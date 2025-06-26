import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import splrep, splev
from scipy.optimize import minimize

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

# Find true apexes by finding zero-crossings of the 3rd derivative
spline_d3 = splev(np.linspace(-10, 15, 1000), tck, der=3)
zero_crossings = np.where(np.diff(np.sign(spline_d3)))[0]
vg_grid = np.linspace(-10, 15, 1000)
cutoff_apex_vg = vg_grid[zero_crossings[0]]
sat_apex_vg = vg_grid[zero_crossings[1]]
cutoff_apex_vp = splev(cutoff_apex_vg, tck)
sat_apex_vp = splev(sat_apex_vg, tck)

print(f"Cutoff Knee Apex (3rd deriv zero): Vg={cutoff_apex_vg:.4f}, Vp={cutoff_apex_vp:.4f}")
print(f"Saturation Shoulder Apex (3rd deriv zero): Vg={sat_apex_vg:.4f}, Vp={sat_apex_vp:.4f}")

def objective_func(c, vg, vp): return np.sum((np.polyval(c, vg) - vp)**2)
all_coeffs = [None] * 3

# --- 4. Independent Constrained Fit for Each Segment ---

# Segment 3 (Saturation)
s3_start, s3_end = b[2], b[3]
vg_s3_fit = np.linspace(s3_start, s3_end, 200)
vp_s3_target = splev(vg_s3_fit, tck)
s3_constraints = [
    {'type': 'eq', 'fun': lambda c: np.polyval(c, s3_start) - splev(s3_start, tck)},
    {'type': 'eq', 'fun': lambda c: np.polyval(np.polyder(c, 1), s3_start) - splev(s3_start, tck, der=1)},
    {'type': 'eq', 'fun': lambda c: np.polyval(c, sat_apex_vg) - sat_apex_vp},
    {'type': 'eq', 'fun': lambda c: np.polyval(c, s3_end) - splev(s3_end, tck)},
    {'type': 'eq', 'fun': lambda c: np.polyval(np.polyder(c, 1), s3_end)}
]
all_coeffs[2] = minimize(objective_func, np.polyfit(vg_s3_fit, vp_s3_target, poly_deg), args=(vg_s3_fit, vp_s3_target), method='SLSQP', constraints=s3_constraints).x
p3 = np.poly1d(all_coeffs[2])

# Segment 2 (Central Bridge)
s2_start, s2_end = b[1], b[2]
vg_s2_fit = np.linspace(s2_start, s2_end, 200)
vp_s2_target = splev(vg_s2_fit, tck)
s2_constraints = [
    {'type': 'eq', 'fun': lambda c: np.polyval(c, s2_start) - splev(s2_start, tck)},
    {'type': 'eq', 'fun': lambda c: np.polyval(np.polyder(c, 1), s2_start) - splev(s2_start, tck, der=1)},
    {'type': 'eq', 'fun': lambda c: np.polyval(c, s2_end) - splev(s2_end, tck)},
    {'type': 'eq', 'fun': lambda c: np.polyval(np.polyder(c, 1), s2_end) - splev(s2_end, tck, der=1)}
]
all_coeffs[1] = minimize(objective_func, np.polyfit(vg_s2_fit, vp_s2_target, poly_deg), args=(vg_s2_fit, vp_s2_target), method='SLSQP', constraints=s2_constraints).x
p2 = np.poly1d(all_coeffs[1])

# Segment 1 (Knee)
s1_start_fit, s1_end_fit = -12.0, b[1]
def robust_b_plus_check(coeffs):
    p_deriv_roots = np.roots(np.polyder(coeffs, 1))
    real_extrema_vg = p_deriv_roots[np.isreal(p_deriv_roots)].real
    valid_extrema_vg = real_extrema_vg[(real_extrema_vg >= s1_start_fit) & (real_extrema_vg <= s1_end_fit)]
    if len(valid_extrema_vg) == 0: return B_PLUS_LIMIT
    return B_PLUS_LIMIT - np.max(np.polyval(coeffs, valid_extrema_vg))
vg_s1_fit = np.linspace(s1_start_fit, s1_end_fit, 200)
vp_s1_target = splev(vg_s1_fit, tck)
s1_constraints = [
    {'type': 'eq', 'fun': lambda c: np.polyval(c, s1_end_fit) - splev(s1_end_fit, tck)},
    {'type': 'eq', 'fun': lambda c: np.polyval(np.polyder(c, 1), s1_end_fit) - splev(s1_end_fit, tck, der=1)},
    {'type': 'eq', 'fun': lambda c: np.polyval(c, b[0]) - splev(b[0], tck)},
    {'type': 'eq', 'fun': lambda c: np.polyval(c, cutoff_apex_vg) - cutoff_apex_vp},
    {'type': 'ineq', 'fun': robust_b_plus_check}
]
all_coeffs[0] = minimize(objective_func, np.polyfit(vg_s1_fit, vp_s1_target, poly_deg), args=(vg_s1_fit, vp_s1_target), method='SLSQP', constraints=s1_constraints).x
p1 = np.poly1d(all_coeffs[0])

# --- 5. Assemble Final Model ---
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

# --- 6. Plot for Verification ---
fig, ax = plt.subplots(figsize=(15, 9))
vg_plot_range = np.linspace(-12, 18, 2000)
ax.plot(Vg_data, Vp_data, 'o', ms=3, label='SPICE Data', color='gray', alpha=0.6)
ax.plot(vg_plot_range, model_12ax7(vg_plot_range), label='Definitive Production Model', color='red', lw=2.5)
pin_vgs = [b[0], cutoff_apex_vg, sat_apex_vg, b[3]]
pin_vps = [p1(b[0]), p1(cutoff_apex_vg), p3(sat_apex_vg), p3(b[3])]
ax.plot(pin_vgs, pin_vps, 'x', color='blue', markersize=10, mew=2, label='Pinned Points')
for bound in b: ax.axvline(bound, color='k', ls=':', alpha=0.7)
ax.axhline(B_PLUS_LIMIT, color='purple', ls='--', label=f'B+ Limit ({B_PLUS_LIMIT}V)')
ax.set(ylabel='Plate Voltage (Vp)', title='Definitive 12AX7 Model - Production Quality', xlim=(-12, 18))
ax.legend(), ax.grid(True)
plt.show()

# --- 7. Print Coefficients ---
print("\n--- Definitive, Final Polynomial Coefficients ---")
names = ["Segment 1: -8V to -1V", "Segment 2: -1V to 1V", "Segment 3: 1V to 15.68V"]
for name, coeffs in zip(names, all_coeffs):
    print(f"\n{name}\n" + "\n".join([f"{c:e}" for c in coeffs]))
print(f"\nClamp Low (< -8.0V): {V_clamp_low:.4f} V"), print(f"Clamp High (> 15.68V): {V_clamp_high:.4f} V")
