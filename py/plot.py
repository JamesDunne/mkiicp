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

# Find Vg for the "cutoff knee"
res_cutoff_knee = minimize_scalar(lambda x: splev(x, tck, der=1), bounds=(-8, -2), method='bounded')
cutoff_knee_vg, cutoff_knee_vp = res_cutoff_knee.x, splev(res_cutoff_knee.x, tck)

# Find Vg for the "saturation shoulder" (point of max curvature)
res_sat_knee = minimize_scalar(lambda x: -splev(x, tck, der=2), bounds=(1, 10), method='bounded')
sat_knee_vg, sat_knee_vp = res_sat_knee.x, splev(res_sat_knee.x, tck)

print(f"Cutoff Knee Pin: Vg={cutoff_knee_vg:.4f}, Vp={cutoff_knee_vp:.4f}")
print(f"Saturation Shoulder Pin: Vg={sat_knee_vg:.4f}, Vp={sat_knee_vp:.4f}")

def objective_func(c, vg, vp): return np.sum((np.polyval(c, vg) - vp)**2)
all_coeffs = [None] * 3

# --- 4. Stage 1: Fit Segment 3 (Saturation) ---
s3_start, s3_end = b[2], b[3]
vg_s3_fit = np.linspace(s3_start, s3_end, 200)
vp_s3_target = splev(vg_s3_fit, tck)
s3_constraints = [
    {'type': 'eq', 'fun': lambda c: np.polyval(c, s3_start) - splev(s3_start, tck)},
    {'type': 'eq', 'fun': lambda c: np.polyval(np.polyder(c, 1), s3_start) - splev(s3_start, tck, der=1)},
    {'type': 'eq', 'fun': lambda c: np.polyval(c, sat_knee_vg) - sat_knee_vp}, # Pin at saturation shoulder
    {'type': 'eq', 'fun': lambda c: np.polyval(c, s3_end) - splev(s3_end, tck)},
    {'type': 'eq', 'fun': lambda c: np.polyval(np.polyder(c, 1), s3_end)}
]
s3_res = minimize(objective_func, np.polyfit(vg_s3_fit, vp_s3_target, poly_deg), args=(vg_s3_fit, vp_s3_target),
                  method='SLSQP', constraints=s3_constraints)
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
    {'type': 'eq', 'fun': lambda c: np.polyval(c, s2_end) - p3(s2_end)},
    {'type': 'eq', 'fun': lambda c: np.polyval(np.polyder(c, 1), s2_end) - p3.deriv(1)(s2_end)}
]
s2_res = minimize(objective_func, np.polyfit(vg_s2_fit, vp_s2_target, poly_deg), args=(vg_s2_fit, vp_s2_target),
                  method='SLSQP', constraints=s2_constraints)
all_coeffs[1] = s2_res.x
p2 = np.poly1d(all_coeffs[1])

# --- 6. Stage 3: Fit Segment 1 (Knee) ---
s1_start_fit, s1_end_fit = -12.0, b[1]
vg_s1_fit = np.linspace(s1_start_fit, s1_end_fit, 200)
vp_s1_target = splev(vg_s1_fit, tck)

# --- Definitive B+ Constraint Function ---
def robust_b_plus_check(coeffs):
    p_deriv = np.polyder(coeffs, 1)
    extrema_candidates = np.roots(p_deriv)
    # Filter for real roots within the segment
    real_extrema_vg = extrema_candidates[np.isreal(extrema_candidates)].real
    valid_extrema_vg = real_extrema_vg[(real_extrema_vg >= s1_start_fit) & (real_extrema_vg <= s1_end_fit)]
    if len(valid_extrema_vg) > 0:
        max_val = np.max(np.polyval(coeffs, valid_extrema_vg))
        return B_PLUS_LIMIT - max_val
    return B_PLUS_LIMIT # No local maxima in range, constraint is met

s1_constraints = [
    {'type': 'eq', 'fun': lambda c: np.polyval(c, s1_end_fit) - p2(s1_end_fit)},
    {'type': 'eq', 'fun': lambda c: np.polyval(np.polyder(c, 1), s1_end_fit) - p2.deriv(1)(s1_end_fit)},
    {'type': 'eq', 'fun': lambda c: np.polyval(np.polyder(c, 2), s1_end_fit) - p2.deriv(2)(s1_end_fit)},
    {'type': 'eq', 'fun': lambda c: np.polyval(c, b[0]) - splev(b[0], tck)},
    {'type': 'eq', 'fun': lambda c: np.polyval(c, cutoff_knee_vg) - cutoff_knee_vp},
    {'type': 'ineq', 'fun': robust_b_plus_check}
]
s1_res = minimize(objective_func, np.polyfit(vg_s1_fit, vp_s1_target, poly_deg), args=(vg_s1_fit, vp_s1_target),
                  method='SLSQP', constraints=s1_constraints)
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
fig, ax = plt.subplots(figsize=(15, 9))
vg_plot_range = np.linspace(-12, 18, 2000)
ax.plot(Vg_data, Vp_data, 'o', ms=3, label='SPICE Data', color='gray', alpha=0.6)
ax.plot(vg_plot_range, model_12ax7(vg_plot_range), label='Definitive Model', color='red', lw=2.5)
# Add markers for all pinned points
pin_vgs = [b[0], cutoff_knee_vg, sat_knee_vg]
pin_vps = [p1(b[0]), p1(cutoff_knee_vg), p3(sat_knee_vg)]
ax.plot(pin_vgs, pin_vps, 'x', color='blue', markersize=10, mew=2, label='Pinned Points')
for bound in b: ax.axvline(bound, color='k', ls=':', alpha=0.7)
ax.axhline(B_PLUS_LIMIT, color='purple', ls='--', label=f'B+ Limit ({B_PLUS_LIMIT}V)')
ax.set(ylabel='Plate Voltage (Vp)', title='Definitive 12AX7 Model - Fully Constrained, Pinned, and Ripple-Free', xlim=(-12, 18))
ax.legend(), ax.grid(True)
plt.show()

# --- 9. Print Coefficients and Verification Data ---
print("\n--- Definitive, Final Polynomial Coefficients ---")
names = ["Segment 1: -8V to -1V", "Segment 2: -1V to 1V", "Segment 3: 1V to 15.68V"]
for name, coeffs in zip(names, all_coeffs):
    print(f"\n{name}\n" + "\n".join([f"{c:e}" for c in coeffs]))
print(f"\nClamp Low (< -8.0V): {V_clamp_low:.4f} V"), print(f"Clamp High (> 15.68V): {V_clamp_high:.4f} V")
