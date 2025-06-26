import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import splrep, splev
from scipy.optimize import minimize

# --- 1. Load and Prepare Data (Corrected) ---
try:
    # Step 1: Load data from the file
    df = pd.read_csv('Draft1.txt', sep=r'\s+')
    # Step 2: Rename the columns immediately
    df.rename(columns={'V(grid)': 'Vg', 'V(plate)': 'Vp'}, inplace=True)
    # Step 3: Now, perform the cleaning operations on the correctly named columns
    df = df.groupby('Vg', as_index=False).mean().sort_values(by='Vg')
except FileNotFoundError:
    print("Error: Draft1.txt not found. Please ensure the data file is in the same directory.")
    exit()

Vg_data, Vp_data = df['Vg'].values, df['Vp'].values

# --- 2. Define Final Segments and Boundaries ---
boundaries = [-8.0, -1.0, 1.0, 15.68266]
poly_deg = 5
B_PLUS_LIMIT = 405.0

# --- 3. Global Smooth Spline Fit for Target Data ---
spline_mask = (Vg_data >= boundaries[0]) & (Vg_data <= boundaries[-1])
tck = splrep(Vg_data[spline_mask], Vp_data[spline_mask], s=0.5)

# --- 4. Fit Segments 2 and 3 using Spline-to-Poly method ---
all_coeffs = [None] * 3
for i in [1, 2]: # Fit S2 and S3 first
    start, end = boundaries[i], boundaries[i+1]
    vg_dense = np.linspace(start, end, 200)
    vp_smooth = splev(vg_dense, tck)
    all_coeffs[i] = np.polyfit(vg_dense, vp_smooth, poly_deg)

coeffs_s2, coeffs_s3 = all_coeffs[1], all_coeffs[2]
p2 = np.poly1d(coeffs_s2)

# --- 5. Fit Segment 1 using Constrained Optimization ---
start_s1, end_s1 = boundaries[0], boundaries[1]
vg_s1_dense = np.linspace(start_s1, end_s1, 200)
vp_s1_target = splev(vg_s1_dense, tck)

def poly_model(coeffs, vg):
    return np.polyval(coeffs, vg)

def objective_func(coeffs, vg, vp):
    return np.sum((poly_model(coeffs, vg) - vp)**2)

p2_deriv1 = p2.deriv(1)
p2_deriv2 = p2.deriv(2)
target_val_at_boundary = p2(end_s1)
target_d1_at_boundary = p2_deriv1(end_s1)
target_d2_at_boundary = p2_deriv2(end_s1)

# CORRECTED CODE
constraints = [
    # C0 constraint (already correct)
    {'type': 'eq', 'fun': lambda c: np.polyval(c, end_s1) - target_val_at_boundary},
    # C1 constraint: Evaluate the 1st derivative polynomial at the boundary
    {'type': 'eq', 'fun': lambda c: np.polyval(np.polyder(c, 1), end_s1) - target_d1_at_boundary},
    # C2 constraint: Evaluate the 2nd derivative polynomial at the boundary
    {'type': 'eq', 'fun': lambda c: np.polyval(np.polyder(c, 2), end_s1) - target_d2_at_boundary},
    # B+ limit constraint (already correct)
    {'type': 'ineq', 'fun': lambda c: B_PLUS_LIMIT - np.polyval(c, vg_s1_dense)}
]

initial_guess = np.polyfit(vg_s1_dense, vp_s1_target, poly_deg)

opt_result = minimize(
    objective_func,
    initial_guess,
    args=(vg_s1_dense, vp_s1_target),
    method='SLSQP',
    constraints=constraints,
    options={'maxiter': 500}
)
coeffs_s1 = opt_result.x
all_coeffs[0] = coeffs_s1
p1 = np.poly1d(coeffs_s1)
p3 = np.poly1d(coeffs_s3)

# --- 6. Assemble Final Model and Clamps ---
V_clamp_low = p1(boundaries[0])
V_clamp_high = p3(boundaries[-1])

def model_12ax7(vg_in):
    vg = np.asarray(vg_in, dtype=float)
    vp_out = np.zeros_like(vg)
    vp_out[vg < boundaries[0]] = V_clamp_low
    vp_out[(vg >= boundaries[0]) & (vg < boundaries[1])] = p1(vg[(vg >= boundaries[0]) & (vg < boundaries[1])])
    vp_out[(vg >= boundaries[1]) & (vg < boundaries[2])] = p2(vg[(vg >= boundaries[1]) & (vg < boundaries[2])])
    vp_out[(vg >= boundaries[2]) & (vg <= boundaries[3])] = p3(vg[(vg >= boundaries[2]) & (vg <= boundaries[3])])
    vp_out[vg > boundaries[3]] = V_clamp_high
    return vp_out

# --- 7. Plot for Verification ---
fig, axes = plt.subplots(2, 1, figsize=(15, 12), gridspec_kw={'height_ratios': [3, 1]})
ax1 = axes[0]
vg_full_range = np.linspace(-24, 24, 2000)
ax1.plot(Vg_data, Vp_data, 'o', markersize=3, label='SPICE Data', color='gray', alpha=0.6)
ax1.plot(vg_full_range, model_12ax7(vg_full_range), label='Final Constrained Model', color='red', linewidth=2.5)
for b in boundaries: ax1.axvline(b, color='k', linestyle=':', alpha=0.7)
ax1.axhline(B_PLUS_LIMIT, color='blue', linestyle='--', label=f'B+ Limit ({B_PLUS_LIMIT}V)')
ax1.set_ylabel('Plate Voltage (Vp)'), ax1.set_title('Final 12AX7 Model with Constraints'), ax1.legend(), ax1.grid(True)
ax1.set_xlim(-24, 24)

ax2 = axes[1]
ax2.plot(Vg_data, Vp_data, 'o', markersize=4, color='gray', alpha=0.6)
ax2.plot(vg_full_range, model_12ax7(vg_full_range), 'r-', linewidth=2.5)
ax2.axhline(B_PLUS_LIMIT, color='blue', linestyle='--')
ax2.set_xlim(-12, -4), ax2.set_ylim(395, 406)
ax2.set_title('Zoomed View: "Knee" Region (No Ripple/Overshoot)')
ax2.set_xlabel('Grid Voltage (Vg)'), ax2.set_ylabel('Plate Voltage (Vp)'), ax2.grid(True)
plt.tight_layout()
plt.show()

# --- 8. Print Coefficients and Verification Data ---
print("\n--- Final, Constrained Polynomial Coefficients ---")
names = ["Segment 1: -8V to -1V (Constrained)", "Segment 2: -1V to 1V", "Segment 3: 1V to 15.68V"]
for name, coeffs in zip(names, all_coeffs):
    print(f"\n{name}")
    for c in coeffs: print(f"{c:e}")
print(f"\nClamp Low (< -8.0V): {V_clamp_low:.4f} V"), print(f"Clamp High (> 15.68V): {V_clamp_high:.4f} V")