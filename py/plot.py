import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import splrep, splev

# --- 1. Load and Prepare Data ---
try:
    df = pd.read_csv('Draft1.txt', sep=r'\s+')
except FileNotFoundError:
    print("Error: Draft1.txt not found. Please ensure the data file is in the same directory.")
    exit()

df.rename(columns={'V(grid)': 'Vg', 'V(plate)': 'Vp'}, inplace=True)
# CRITICAL FIX: Remove duplicate Vg points to ensure valid differentiation
df = df.groupby('Vg', as_index=False).mean()
df.sort_values(by='Vg', inplace=True)

Vg_data = df['Vg'].values
Vp_data = df['Vp'].values

# --- 2. Define User-Specified Segments and Boundaries ---
boundaries = [-8.0, -1.0, 1.0, 10.0]
poly_deg = 5

# --- 3. Global Smooth Spline Fit ---
# Filter data to the specified modeling range
spline_mask = (Vg_data >= boundaries[0]) & (Vg_data <= boundaries[-1])
# splrep creates a B-spline representation. s is a smoothing factor.
# A small s value ensures the spline fits the data points closely.
tck = splrep(Vg_data[spline_mask], Vp_data[spline_mask], s=0.5)

# --- 4. Generate Polynomials from the Smooth Spline ---
all_coeffs = []
for i in range(len(boundaries) - 1):
    start, end = boundaries[i], boundaries[i+1]
    # Create a dense grid of Vg points for the current segment
    vg_dense = np.linspace(start, end, 200)
    # Evaluate the smooth spline on this dense grid to get "perfect" data
    vp_smooth = splev(vg_dense, tck)
    # Fit a polynomial to this smooth, dense data
    coeffs = np.polyfit(vg_dense, vp_smooth, poly_deg)
    all_coeffs.append(coeffs)

coeffs_s1, coeffs_s2, coeffs_s3 = all_coeffs
p1 = np.poly1d(coeffs_s1)
p2 = np.poly1d(coeffs_s2)
p3 = np.poly1d(coeffs_s3)

# --- 5. Define Clamp Values and Full Model ---
V_clamp_low = p1(boundaries[0])
V_clamp_high = p3(boundaries[-1])

def model_12ax7(vg_in):
    vg = np.asarray(vg_in, dtype=float)
    vp_out = np.zeros_like(vg)
    
    vp_out[vg < boundaries[0]] = V_clamp_low
    mask1 = (vg >= boundaries[0]) & (vg < boundaries[1])
    vp_out[mask1] = p1(vg[mask1])
    mask2 = (vg >= boundaries[1]) & (vg < boundaries[2])
    vp_out[mask2] = p2(vg[mask2])
    mask3 = (vg >= boundaries[2]) & (vg <= boundaries[3])
    vp_out[mask3] = p3(vg[mask3])
    vp_out[vg > boundaries[3]] = V_clamp_high
    
    return vp_out

# --- 6. Plot for Verification ---
plt.style.use('seaborn-v0_8-whitegrid')
fig, ax = plt.subplots(figsize=(14, 8))

vg_full_range = np.linspace(-24, 24, 1000)
ax.plot(Vg_data, Vp_data, 'o', markersize=3, label='SPICE Data (Cleaned)', color='gray', alpha=0.6)
ax.plot(vg_full_range, model_12ax7(vg_full_range), label='Final Piecewise Polynomial Model', color='red', linewidth=2.5)
for b in boundaries:
    ax.axvline(b, color='k', linestyle=':', alpha=0.7)

ax.set_xlabel('Grid Voltage (Vg)')
ax.set_ylabel('Plate Voltage (Vp)')
ax.set_title('12AX7 Model with User-Defined Segments and C2 Continuity')
ax.legend()
ax.grid(True)
ax.set_xlim(-24, 24)
plt.show()

# --- 7. Print Coefficients and Verification Data ---
print("\n--- Final Polynomial Coefficients (a*x^5 + b*x^4 + ...) ---")
names = ["Segment 1: -8V to -1V", "Segment 2: -1V to 1V", "Segment 3: 1V to 10V"]
for name, coeffs in zip(names, all_coeffs):
    print(f"\n{name}")
    for c in coeffs: print(f"{c:e}")

print(f"\nClamp Low (< -8.0V): {V_clamp_low:.4f} V")
print(f"Clamp High (> 10.0V): {V_clamp_high:.4f} V")

# --- 8. Boundary Continuity Check ---
print("\n--- Boundary Continuity Check ---")
boundaries_to_check = [-1.0, 1.0]
polys = [p1, p2, p3]
for i, b in enumerate(boundaries_to_check):
    p_left, p_right = polys[i], polys[i+1]
    print(f"\nAt Vg = {b:.1f} V (Boundary between S{i+1} and S{i+2}):")
    print(f"  Vp (Left):     {p_left(b):.4f} V")
    print(f"  Vp (Right):    {p_right(b):.4f} V")
    print(f"  Gain (Left):   {p_left.deriv(1)(b):.4f}")
    print(f"  Gain (Right):  {p_right.deriv(1)(b):.4f}")
    print(f"  Curve (Left):  {p_left.deriv(2)(b):.4f}")
    print(f"  Curve (Right): {p_right.deriv(2)(b):.4f}")
