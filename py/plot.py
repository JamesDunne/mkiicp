import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import splrep, splev
from scipy.optimize import minimize, minimize_scalar

# --- 1. Load and Prepare Data ---
try:
    df = pd.read_csv('Draft1.txt', sep=r'\s+').rename(columns={'V(grid)': 'Vg', 'V(plate)': 'Vp'})
    df = df.groupby('Vg', as_index=False).mean().sort_values(by='Vg')
except FileNotFoundError:
    print("Error: Draft1.txt not found.")
    exit()

Vg_data, Vp_data = df['Vg'].values, df['Vp'].values

# --- 2. Define Final Segments and Boundaries ---
b = operational_boundaries = [-8.0, -1.0, 1.0, 15.68266]
poly_deg = 5
B_PLUS_LIMIT = 405.0

# --- 3. Master Spline and Key Point Identification ---
spline_mask = (Vg_data >= -12.0) & (Vg_data <= b[-1])
tck = splrep(Vg_data[spline_mask], Vp_data[spline_mask], s=0.5)

res = minimize_scalar(lambda x: splev(x, tck, der=1), bounds=(-8, -2), method='bounded')
knee_vg = res.x
knee_vp = splev(knee_vg, tck)
print(f"Identified Knee Pinning Point: Vg={knee_vg:.4f}V, Vp={knee_vp:.4f}V")

def objective_func(c, vg, vp): return np.sum((np.polyval(c, vg) - vp)**2)
all_coeffs = [None] * 3

# --- 4. Stage 1: Fit Segment 3 (Saturation, C1-Constrained) ---
s3_start, s3_end = b[2], b[3]
vg_s3_fit = np.linspace(s3_start, s3_end, 200)
vp_s3_target = splev(vg_s3_fit, tck)

# --- CHANGE: Relax C2 constraint to eliminate ripple ---
s3_constraints = [
    {'type': 'eq', 'fun': lambda c: np.polyval(c, s3_start) - splev(s3_start, tck)},       # C0
    {'type': 'eq', 'fun': lambda c: np.polyval(np.polyder(c, 1), s3_start) - splev(s3_start, tck, der=1)}, # C1
    {'type': 'eq', 'fun': lambda c: np.polyval(c, s3_end) - splev(s3_end, tck)},
    {'type': 'eq', 'fun': lambda c: np.polyval(np.polyder(c, 1), s3_end)}
]
s3_res = minimize(objective_func, np.polyfit(vg_s3_fit, vp_s3_target, poly_deg), args=(vg_s3_fit, vp_s3_target),
                  method='SLSQP', constraints=s3_constraints, options={'maxiter': 1000})
all_coeffs[2] = s3_res.x
p3 = np.poly1d(all_coeffs[2])

# --- 5. Stage 2: Fit Segment 2 (Central) ---
s2_start, s2_end = b[1], b[2]
vg_s2_fit = np.linspace(s2_start, s2_end, 200)
vp_s2_target = splev(vg_s2_fit, tck)

s2_constraints = [
    {'type': 'eq', 'fun': lambda c: np.polyval(c, s2_start) - splev(s2_start, tck)},
    {'type': 'eq', 'fun': lambda c: np.polyval(np.polyder(c, 1), s2_start) - splev(s2_start, tck, der=1)},
    {'type': 'eq', 'fun': lambda c: np.polyval(np.polyder(c, 2), s2_start) - splev(s2_start, tck, der=2)},
    {'type': 'eq', 'fun': lambda c: np.polyval(c, s2_end) - p3(s2_end)}, # Match new S3
    {'type': 'eq', 'fun': lambda c: np.polyval(np.polyder(c, 1), s2_end) - p3.deriv(1)(s2_end)}
]
s2_res = minimize(objective_func, np.polyfit(vg_s2_fit, vp_s2_target, poly_deg), args=(vg_s2_fit, vp_s2_target),
                  method='SLSQP', constraints=s2_constraints, options={'maxiter': 1000})
all_coeffs[1] = s2_res.x
p2 = np.poly1d(all_coeffs[1])

# --- 6. Stage 3: Fit Segment 1 (Knee) ---
s1_start, s1_end = -12.0, b[1]
# --- FIX: Use a much denser grid for B+ limit check ---
vg_s1_fit = np.linspace(s1_start, s1_end, 1000)
vp_s1_target = splev(vg_s1_fit, tck)

s1_constraints = [
    {'type': 'eq', 'fun': lambda c: np.polyval(c, s1_end) - p2(s1_end)},
    {'type': 'eq', 'fun': lambda c: np.polyval(np.polyder(c, 1), s1_end) - p2.deriv(1)(s1_end)},
    {'type': 'eq', 'fun': lambda c: np.polyval(np.polyder(c, 2), s1_end) - p2.deriv(2)(s1_end)},
    {'type': 'eq', 'fun': lambda c: np.polyval(c, b[0]) - splev(b[0], tck)},
    {'type': 'eq', 'fun': lambda c: np.polyval(c, knee_vg) - knee_vp},
    {'type': 'ineq', 'fun': lambda c: B_PLUS_LIMIT - np.polyval(c, vg_s1_fit)}
]
s1_res = minimize(objective_func, np.polyfit(vg_s1_fit, vp_s1_target, poly_deg), args=(vg_s1_fit, vp_s1_target),
                  method='SLSQP', constraints=s1_constraints, options={'maxiter': 1000})
all_coeffs[0] = s1_res.x
p1 = np.poly1d(all_coeffs[0])

# --- 7. Assemble Final Model and Clamps ---
V_clamp_low = p1(b[0])
V_clamp_high = p3(b[3])

def model_12ax7(vg_in):
    vg = np.asarray(vg_in, dtype=float)
    vp_out = np.zeros_like(vg)
    vp_out[vg < b[0]] = V_clamp_low
    mask1 = (vg >= b[0]) & (vg < b[1]); vp_out[mask1] = p1(vg[mask1])
    mask2 = (vg >= b[1]) & (vg < b[2]); vp_out[mask2] = p2(vg[mask2])
    mask3 = (vg >= b[2]) & (vg <= b[3]); vp_out[mask3] = p3(vg[mask3])
    vp_out[vg > b[3]] = V_clamp_high
    return vp_out

# --- 8. Plot for Verification ---
fig, axes = plt.subplots(2, 1, figsize=(16, 9), gridspec_kw={'height_ratios': [1, 1]})
ax1, ax2 = axes
vg_plot_range = np.linspace(-10, 18, 2000)
ax1.plot(Vg_data, Vp_data, 'o', ms=3, label='SPICE Data', color='gray', alpha=0.6)
ax1.plot(vg_plot_range, model_12ax7(vg_plot_range), label='Definitive Model', color='red', lw=2.5)
ax1.plot([b[0], knee_vg], [p1(b[0]), knee_vp], 'x', color='blue', markersize=10, mew=2, label='Pinned Points')
for bound in b: ax1.axvline(bound, color='k', ls=':', alpha=0.7)
ax1.set(ylabel='Plate Voltage (Vp)', title='Definitive 12AX7 Model - Fully Constrained & Ripple-Free', xlim=(-10, 18))
ax1.legend(), ax1.grid(True)

ax2.plot(Vg_data, Vp_data, 'o', ms=4, color='gray', alpha=0.6)
ax2.plot(vg_plot_range, model_12ax7(vg_plot_range), 'r-', lw=2.5)
ax2.set(xlim=(0, 18), ylim=(10, 140), title='Zoomed View: Saturation Region is Now Smooth',
        xlabel='Grid Voltage (Vg)', ylabel='Plate Voltage (Vp)')
ax2.grid(True)
plt.tight_layout()
plt.show()

# --- 9. Print Coefficients and Verification Data ---
print("\n--- Definitive, Final Polynomial Coefficients ---")
names = ["Segment 1: -8V to -1V", "Segment 2: -1V to 1V", "Segment 3: 1V to 15.68V"]
for name, coeffs in zip(names, all_coeffs):
    print(f"\n{name}\n" + "\n".join([f"{c:e}" for c in coeffs]))
print(f"\nClamp Low (< -8.0V): {V_clamp_low:.4f} V"), print(f"Clamp High (> 15.68V): {V_clamp_high:.4f} V")
