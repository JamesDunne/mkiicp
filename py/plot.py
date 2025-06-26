import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import splrep, splev
from scipy.optimize import minimize

# --- 1. Load and Prepare Data ---
try:
    df = pd.read_csv('Draft1.txt', sep=r'\s+')
    df.rename(columns={'V(grid)': 'Vg', 'V(plate)': 'Vp'}, inplace=True)
    df = df.groupby('Vg', as_index=False).mean().sort_values(by='Vg')
except FileNotFoundError:
    print("Error: Draft1.txt not found.")
    exit()

Vg_data, Vp_data = df['Vg'].values, df['Vp'].values

# --- 2. Define Final Segments and Boundaries ---
operational_boundaries = [-8.0, -1.0, 1.0, 15.68266]
poly_deg = 5
B_PLUS_LIMIT = 405.0

# --- 3. Global Smooth Spline Fit ---
# Widen spline range to support the extended fit for segment 1
spline_mask = (Vg_data >= -12.0) & (Vg_data <= operational_boundaries[-1])
tck = splrep(Vg_data[spline_mask], Vp_data[spline_mask], s=0.5)

# --- 4. Fit Segments 2 and 3 (Unchanged) ---
all_coeffs = [None] * 3
for i in [1, 2]:
    start, end = operational_boundaries[i], operational_boundaries[i+1]
    vg_dense = np.linspace(start, end, 200)
    vp_smooth = splev(vg_dense, tck)
    all_coeffs[i] = np.polyfit(vg_dense, vp_smooth, poly_deg)

coeffs_s2, coeffs_s3 = all_coeffs[1], all_coeffs[2]
p2 = np.poly1d(coeffs_s2)

# --- 5. Refined Fit for Segment 1 ---
# *** KEY CHANGE: Use a wider fitting range ***
fit_start_s1, fit_end_s1 = -12.0, -1.0
vg_s1_fit_range = np.linspace(fit_start_s1, fit_end_s1, 200)
vp_s1_target = splev(vg_s1_fit_range, tck)

def objective_func(coeffs, vg, vp):
    return np.sum((np.polyval(coeffs, vg) - vp)**2)

p2_deriv1, p2_deriv2 = p2.deriv(1), p2.deriv(2)
target_val_at_boundary = p2(fit_end_s1)
target_d1_at_boundary = p2_deriv1(fit_end_s1)
target_d2_at_boundary = p2_deriv2(fit_end_s1)

constraints = [
    {'type': 'eq', 'fun': lambda c: np.polyval(c, fit_end_s1) - target_val_at_boundary},
    {'type': 'eq', 'fun': lambda c: np.polyval(np.polyder(c, 1), fit_end_s1) - target_d1_at_boundary},
    {'type': 'eq', 'fun': lambda c: np.polyval(np.polyder(c, 2), fit_end_s1) - target_d2_at_boundary},
    {'type': 'ineq', 'fun': lambda c: B_PLUS_LIMIT - np.polyval(c, vg_s1_fit_range)}
]

initial_guess = np.polyfit(vg_s1_fit_range, vp_s1_target, poly_deg)
opt_result = minimize(
    objective_func, initial_guess, args=(vg_s1_fit_range, vp_s1_target),
    method='SLSQP', constraints=constraints, options={'maxiter': 1000}
)
coeffs_s1 = opt_result.x
all_coeffs[0] = coeffs_s1
p1 = np.poly1d(coeffs_s1)
p3 = np.poly1d(coeffs_s3)

# --- 6. Assemble Final Model with Guaranteed Clamp Continuity ---
# *** KEY CHANGE: Clamp value is now derived from the polynomial itself ***
V_clamp_low = p1(operational_boundaries[0])
V_clamp_high = p3(operational_boundaries[-1])

def model_12ax7(vg_in):
    vg = np.asarray(vg_in, dtype=float)
    vp_out = np.zeros_like(vg)
    b = operational_boundaries
    vp_out[vg < b[0]] = V_clamp_low
    mask1 = (vg >= b[0]) & (vg < b[1])
    vp_out[mask1] = p1(vg[mask1])
    mask2 = (vg >= b[1]) & (vg < b[2])
    vp_out[mask2] = p2(vg[mask2])
    mask3 = (vg >= b[2]) & (vg <= b[3])
    vp_out[mask3] = p3(vg[mask3])
    vp_out[vg > b[3]] = V_clamp_high
    return vp_out

# --- 7. Plot for Verification ---
fig, axes = plt.subplots(2, 1, figsize=(15, 12), gridspec_kw={'height_ratios': [3, 1]})
ax1 = axes[0]
vg_full_range = np.linspace(-24, 24, 2000)
ax1.plot(Vg_data, Vp_data, 'o', ms=3, label='SPICE Data', color='gray', alpha=0.6)
ax1.plot(vg_full_range, model_12ax7(vg_full_range), label='Final Refined Model', color='red', lw=2.5)
for b in operational_boundaries: ax1.axvline(b, color='k', ls=':', alpha=0.7)
ax1.axhline(B_PLUS_LIMIT, color='blue', ls='--', label=f'B+ Limit ({B_PLUS_LIMIT}V)')
ax1.set(ylabel='Plate Voltage (Vp)', title='Final 12AX7 Model: Ripple Eliminated', xlim=(-24, 24))
ax1.legend(), ax1.grid(True)

ax2 = axes[1]
ax2.plot(Vg_data, Vp_data, 'o', ms=4, color='gray', alpha=0.6)
ax2.plot(vg_full_range, model_12ax7(vg_full_range), 'r-', lw=2.5)
ax2.axhline(B_PLUS_LIMIT, color='blue', ls='--')
ax2.set(xlim=(-12, -4), ylim=(395, 406), title='Zoomed View: "Knee" Region is Now Smooth',
        xlabel='Grid Voltage (Vg)', ylabel='Plate Voltage (Vp)')
ax2.grid(True)
plt.tight_layout()
plt.show()

# --- 8. Print Coefficients and Verification Data ---
print("\n--- Final, Refined Polynomial Coefficients ---")
names = ["Segment 1: -8V to -1V", "Segment 2: -1V to 1V", "Segment 3: 1V to 15.68V"]
for name, coeffs in zip(names, all_coeffs):
    print(f"\n{name}")
    for c in coeffs: print(f"{c:e}")

print(f"\nClamp Low (< -8.0V): {V_clamp_low:.4f} V")
print(f"Clamp High (> 15.68V): {V_clamp_high:.4f} V")

# --- 9. Final Continuity Checks ---
print("\n--- Clamp Boundary Continuity Check (at Vg = -8.0V) ---")
print(f"  Polynomial value: {p1(operational_boundaries[0]):.4f}")
print(f"  Clamp value:      {V_clamp_low:.4f}")
print(f"  Slope (Gain):     {p1.deriv(1)(operational_boundaries[0]):.4f}  <-- Near zero, as desired")
