import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import iv # Modified Bessel function of the first kind
from scipy.integrate import quad

# --- Configuration ---
# You can change these parameters as needed
INPUT_FILE = 'data/traces.txt'
OUTPUT_JSON_FILE = 'data/phase_data.json'
OUTPUT_PLOT_FILE = 'data/phase__plot.png'
VMF_COMPONENTS = 3# Number of von Mises-Fisher distributions to fit
EM_MAX_ITER = 1000   # Maximum iterations for the EM algorithm
EM_TOLERANCE = 1e-6 # Convergence tolerance for the EM algorithm

# --- Data Loading ---
def load_data(filepath):
    """
    Loads space-separated data from a text file (x y z throughput).

    Args:
        filepath (str): The path to the data file.

    Returns:
        tuple: A tuple containing:
            - np.ndarray: A (n, 3) array of normalized direction vectors.
            - np.ndarray: A (n,) array of path throughputs (weights).
    """
    try:
        data = np.loadtxt(filepath)
        if data.ndim == 1: # Handle case with a single line of data
            data = data.reshape(1, -1)
        
        directions = data[:, 0:3]
        # Normalize direction vectors to be safe
        norms = np.linalg.norm(directions, axis=1, keepdims=True)
        # Avoid division by zero for null vectors
        norms[norms == 0] = 1
        directions /= norms
        
        weights = data[:, 3]
        print(f"Successfully loaded {len(weights)} data points from '{filepath}'.")
        return directions, weights
    except IOError:
        print(f"Error: Could not find or read the file at '{filepath}'.")
        return None, None
    except ValueError:
        print(f"Error: Data in '{filepath}' is not in the correct format (four numbers per line).")
        return None, None


# --- Henyey-Greenstein (HG) Fitting ---

def fit_hg(directions, weights):
    """
    Fits the data to a Henyey-Greenstein phase function.
    The asymmetry parameter 'g' is the weighted average of cos(theta).
    cos(theta) is the dot product with the incoming direction (0,0,1), which is the z-component.

    Args:
        directions (np.ndarray): Array of outgoing direction vectors.
        weights (np.ndarray): Array of corresponding weights.

    Returns:
        float: The fitted asymmetry parameter 'g'.
    """
    # cos(theta) is the z-component of the direction vector
    cos_thetas = directions[:, 2]
    g = np.sum(weights * cos_thetas) / np.sum(weights)
    return g

def hg_pdf_theta(theta, g):
    """
    Henyey-Greenstein PDF as a function of theta.
    """
    cos_theta = np.cos(theta)
    numerator = 1 - g**2
    denominator = 2 * (1 + g**2 - 2 * g * cos_theta)**(1.5)
    # The Jacobian for spherical coordinates integration is sin(theta)
    return np.sin(theta) * numerator / denominator


# --- von Mises-Fisher (vMF) Mixture Model Fitting (EM Algorithm) ---

def _vmf_kappa_approximation(R_bar, D=3):
    """
    Approximation for kappa (concentration parameter) for vMF.
    """
    if R_bar < 1e-5:
        return 0.0
    return R_bar * (D - R_bar**2) / (1 - R_bar**2)

def fit_vmf_mixture(directions, data_weights, n_components):
    """
    Fits a mixture of von Mises-Fisher distributions using Expectation-Maximization.

    Args:
        directions (np.ndarray): Array of 3D direction vectors.
        data_weights (np.ndarray): Weights for each data point.
        n_components (int): The number of vMF components to fit.

    Returns:
        list: A list of dictionaries, each containing the parameters
              'mu', 'kappa', and 'weight' for a vMF component.
    """
    print(f"Starting EM algorithm to fit {n_components} vMF components...")
    n_points = len(directions)
    X = directions # Use direction vectors directly

    # 1. Initialization
    weights = np.ones(n_components) / n_components
    random_indices = np.random.choice(n_points, n_components, replace=False)
    mus = X[random_indices, :]
    kappas = np.full(n_components, 10.0)
    log_likelihood_old = 0
    
    for i in range(EM_MAX_ITER):
        # 2. E-Step: Calculate responsibilities
        responsibilities = np.zeros((n_points, n_components))
        for k in range(n_components):
            kappa = kappas[k]
            mu = mus[k]
            # Numerically stable calculation of log PDF
            if kappa > 0:
                log_C_kappa = np.log(kappa) - np.log(4 * np.pi) - np.log(np.sinh(kappa))
            else:
                log_C_kappa = -np.log(4 * np.pi)
            
            log_pdf = log_C_kappa + kappa * (X @ mu)
            responsibilities[:, k] = weights[k] * np.exp(log_pdf)
        
        # Check for invalid values and normalize
        responsibilities[np.isnan(responsibilities)] = 0
        row_sums = np.sum(responsibilities, axis=1, keepdims=True)
        log_likelihood_new = np.sum(data_weights * np.log(row_sums.flatten()))
        # Avoid division by zero
        responsibilities = np.divide(responsibilities, row_sums, out=np.zeros_like(responsibilities), where=row_sums!=0)

        weighted_responsibilities = responsibilities * data_weights[:, np.newaxis]

        # 3. M-Step: Update parameters
        weights = np.sum(weighted_responsibilities, axis=0) / np.sum(data_weights)

        for k in range(n_components):
            r_vec = np.sum(weighted_responsibilities[:, k, np.newaxis] * X, axis=0)
            r_norm = np.linalg.norm(r_vec)
            mus[k] = r_vec / r_norm if r_norm > 0 else np.array([0, 0, 1])
            N_k = np.sum(weighted_responsibilities[:, k])
            R_bar = r_norm / N_k if N_k > 0 else 0
            kappas[k] = _vmf_kappa_approximation(R_bar)

        # 4. Convergence Check
        if i > 0 and np.abs(log_likelihood_new - log_likelihood_old) < EM_TOLERANCE:
            print(f"EM algorithm converged after {i+1} iterations.")
            break
        log_likelihood_old = log_likelihood_new
        if i == EM_MAX_ITER - 1:
            print("EM algorithm reached max iterations.")

    results = []
    for k in range(n_components):
        results.append({
            "mu": mus[k].tolist(),
            "kappa": float(kappas[k]),
            "weight": float(weights[k])
        })
    return results

def vmf_pdf_theta(theta, kappa, mu):
    """
    Calculates the marginal vMF PDF as a function of theta (angle from z-axis)
    by integrating over the azimuthal angle phi.
    """
    if kappa == 0:
        return np.sin(theta) / 2.0

    mu_x, mu_y, mu_z = mu
    
    # Integrand for the integral over phi
    integrand = lambda phi: np.exp(kappa * (
        np.sin(theta) * (mu_x * np.cos(phi) + mu_y * np.sin(phi)) +
        np.cos(theta) * mu_z
    ))
    
    integral_val, _ = quad(integrand, 0, 2 * np.pi)
    
    C_kappa = kappa / (4 * np.pi * np.sinh(kappa))
    
    # The Jacobian sin(theta) is included to make it a density over theta
    return C_kappa * np.sin(theta) * integral_val

def vmf_mixture_pdf_theta(theta_range, vmf_params):
    """
    Calculates the combined PDF for a mixture of vMF distributions.
    """
    total_pdf = np.zeros_like(theta_range, dtype=float)
    for params in vmf_params:
        # Calculate the PDF for this component across all theta values
        # using a list comprehension, which is clearer and avoids the vectorization error.
        component_pdf = np.array([
            vmf_pdf_theta(theta, params['kappa'], np.array(params['mu'])) for theta in theta_range
        ])
        total_pdf += params['weight'] * component_pdf
    return total_pdf

def save_phase_data(vmf_params,g , filepath):
    """ Saves the fitted vMF parameters to a JSON file. """
    with open(filepath, 'w') as f:
        json.dump({"vmf" : vmf_params, "hg_g": g}, f, indent=4)
    print(f"vMF parameters saved to '{filepath}'.")

# --- Plotting ---

def plot_results(directions, weights, hg_g, vmf_params):
    """ Generates a plot comparing the data and the fitted distributions. """
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(12, 8))

    # Calculate thetas from direction z-component for histogram
    thetas = np.arccos(np.clip(directions[:, 2], -1.0, 1.0))

    ax.hist(thetas, bins=100, weights=weights, density=True,
            alpha=0.6, label='Original Data (Histogram)', color='gray')

    plot_theta = np.linspace(0, np.pi, 200) # Fewer points for faster plotting

    hg_y = hg_pdf_theta(plot_theta, hg_g)
    ax.plot(plot_theta, hg_y, 'r--', linewidth=2,
            label=f'Henyey-Greenstein Fit (g = {hg_g:.4f})')

    vmf_y = vmf_mixture_pdf_theta(plot_theta, vmf_params)
    ax.plot(plot_theta, vmf_y, 'b-', linewidth=2,
            label=f'{len(vmf_params)}-Component vMF Mixture Fit')
            
    for i, params in enumerate(vmf_params):
        mu = np.array(params['mu'])
        theta_mu = np.degrees(np.arccos(mu[2]))
        phi_mu = np.degrees(np.arctan2(mu[1], mu[0]))

        # Calculate component PDF using a list comprehension to avoid vectorization issues.
        comp_y = params['weight'] * np.array([
            vmf_pdf_theta(theta, params['kappa'], mu) for theta in plot_theta
        ])
        ax.plot(plot_theta, comp_y, linestyle=':', linewidth=1.5,
                label=f'vMF Comp {i+1} (κ={params["kappa"]:.1f}, w={params["weight"]:.2f}, μ(θ,φ)={theta_mu:.1f}°,{phi_mu:.1f}°)')

    ax.set_title('Phase Function Fitting Comparison', fontsize=16)
    ax.set_xlabel('Theta (radians) - Angle from (0,0,1)', fontsize=12)
    ax.set_ylabel('Probability Density', fontsize=12)
    ax.legend(fontsize=10)
    ax.set_xlim(0, np.pi)
    ax.set_ylim(bottom=0)
    
    plt.tight_layout()
    plt.savefig(OUTPUT_PLOT_FILE)
    print(f"Plot saved to '{OUTPUT_PLOT_FILE}'.")
    plt.show()

# --- Main Execution ---

if __name__ == "__main__":
    directions, weights = load_data(INPUT_FILE)

    if directions is not None:
        g = fit_hg(directions, weights)
        print("\n--- Henyey-Greenstein Fit ---")
        print(f"Asymmetry Parameter (g): {g:.6f}")

        vmf_params = fit_vmf_mixture(directions, weights, n_components=VMF_COMPONENTS)
        print("\n--- von Mises-Fisher Mixture Fit ---")
        for i, p in enumerate(vmf_params):
            print(f"Component {i+1}:")
            print(f"  Direction (mu): [{p['mu'][0]:.4f}, {p['mu'][1]:.4f}, {p['mu'][2]:.4f}]")
            print(f"  Concentration (kappa): {p['kappa']:.4f}")
            print(f"  Weight: {p['weight']:.4f}")
        
        save_phase_data(vmf_params, g, OUTPUT_JSON_FILE)

        plot_results(directions, weights, g, vmf_params)


