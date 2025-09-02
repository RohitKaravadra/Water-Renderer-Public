# save as compare_deltaE2000_static.py
import os
import cv2
import numpy as np
import matplotlib.pyplot as plt
from skimage.color import rgb2lab, deltaE_ciede2000

def load_rgb(path):
    img = cv2.imread(path, cv2.IMREAD_UNCHANGED)
    if img is None:
        raise FileNotFoundError(f"Can't open: {path}")
    if img.dtype == np.uint8:
        img = img.astype(np.float32) / 255.0
    else:
        img = img.astype(np.float32)
        if img.max() > 1.0:
            img = img / img.max()
    return cv2.cvtColor(img, cv2.COLOR_BGR2RGB)

def normalize_map(m, clip_percent=(1,99)):
    lo, hi = np.percentile(m, clip_percent)
    m = np.clip(m, lo, hi)
    return (m - lo) / (hi - lo + 1e-9)

def compare_deltaE2000(img1_path, img2_path, out_dir="out"):
    os.makedirs(out_dir, exist_ok=True)

    img1 = load_rgb(img1_path)
    img2 = load_rgb(img2_path)

    if img1.shape != img2.shape:
        img2 = cv2.resize(img2, (img1.shape[1], img1.shape[0]), interpolation=cv2.INTER_LINEAR)

    lab1 = rgb2lab(img1)
    lab2 = rgb2lab(img2)
    deltaE = deltaE_ciede2000(lab1, lab2)

    mean_val = np.mean(deltaE)
    max_val  = np.max(deltaE)
    print(f"Mean ΔE2000: {mean_val:.4f}")
    print(f"Max  ΔE2000: {max_val:.4f}")

    # normalize for visualization
    deltaE_norm = normalize_map(deltaE)

    # --- save heatmap with colorbar ---
    fig, ax = plt.subplots()
    im = ax.imshow(deltaE_norm, cmap="inferno")
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Normalized ΔE2000")
    ax.set_title("ΔE2000 Heatmap")
    ax.axis("off")
    heatmap_path = os.path.join(out_dir, "deltaE2000_heatmap.png")
    plt.savefig(heatmap_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    # --- save overlay ---
    cmap = plt.get_cmap("inferno")
    heat_rgb = cmap(deltaE_norm)[...,:3]
    overlay = (0.6 * heat_rgb + 0.4 * img1)
    overlay = np.clip(overlay, 0, 1)
    overlay_path = os.path.join(out_dir, "deltaE2000_overlay.png")
    plt.imsave(overlay_path, overlay)

    print(f"Saved heatmap: {heatmap_path}")
    print(f"Saved overlay: {overlay_path}")

    # --- show results at the end ---
    plt.figure(figsize=(15,5))
    plt.subplot(1,3,1); plt.imshow(img1); plt.title("Image 1"); plt.axis("off")
    plt.subplot(1,3,2); plt.imshow(img2); plt.title("Image 2"); plt.axis("off")
    plt.subplot(1,3,3); plt.imshow(deltaE_norm, cmap="inferno"); plt.title("ΔE2000 Heatmap"); plt.axis("off")
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    folder = input("Enter folder path (where both images are located): ").strip()
    file1  = input("Enter filename for image 1: ").strip()
    file2  = input("Enter filename for image 2: ").strip()
    path1 = os.path.join(folder, file1)
    path2 = os.path.join(folder, file2)
    compare_deltaE2000(path1, path2, out_dir=os.path.join(folder, "out"))
