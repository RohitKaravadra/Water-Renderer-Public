#pragma once


struct OpticalProperties
{
	float ior;			// index of refraction
	float ior_imag;		// imaginary part of index of refraction
	float sigma_a;		// absorption coefficient
	float sigma_s;		// scattering coefficient
	float sigma_t;		// extinction coefficient
	float g;			// asymmetry parameter for Henyey-Greenstein phase function
};

// A: Physically realistic (coastal plankton, visible light ~440–550 nm)
constexpr OpticalProperties Plankton_Physical{
	1.04f,        // ior (real part)
	0.0005f,      // ior_imag
	1e-7f,        // sigma_a (0.1 m^-1)
	5e-7f,        // sigma_s (0.5 m^-1)
	6e-7f,        // sigma_t
	0.94f         // g
};

// B: Bloom conditions (stronger absorption & scattering)
constexpr OpticalProperties Plankton_Bloom{
	1.05f,        // ior (real part)
	0.001f,       // ior_imag
	5e-7f,        // sigma_a (0.5 m^-1)
	2e-6f,        // sigma_s (2.0 m^-1)
	2.5e-6f,      // sigma_t
	0.96f         // g
};

// C: Renderer-friendly scaled values (not physically accurate,
//     but closer to your original constants for visuals)
constexpr OpticalProperties Plankton_Scaled{
	1.04f,        // ior
	0.0005f,      // ior_imag
	0.02f,        // sigma_a
	0.005f,       // sigma_s
	0.025f,       // sigma_t
	0.90f         // g
};