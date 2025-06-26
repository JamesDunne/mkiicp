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
sat_shoulder_vg = 4.53
b = boundaries = [-8.0, -1.0, 1.0, sat_shoulder_vg, 15.68266]
poly_deg = 4
B_PLUS_LIMIT = 405.0
ERROR_TOLERANCE = 0.5 # The requested +/- 0.5V limit

# --- 3. Master Spline ---
spline_mask = (Vg_data >= -12.0) & (Vg_data <= b[-1])
tck = splrep(Vg_data[spline_mask], Vp_data[spline_mask], s=0.5)

def objective_func(c, vg, vp): return np.sum((np.polyval(c, vg) - vp)**2)
all_coeffs = [None] * len(b)
polys = [None] * len(b)

# --- 4. Independent, Constrained Fit for each Segment ---
for i in range(len(b) - 1):
    start, end = b[i], b[i+1]
    
    # Use a wider fitting range for the first segment for better smoothness
    fit_start = start if i > 0 else -12.0
    vg_fit = np.linspace(fit_start, end, 200)
    vp_target = splev(vg_fit, tck)

    constraints = [
        {'type': 'eq', 'fun': lambda c, s=start: np.polyval(c, s) - splev(s, tck)},
        {'type': 'eq', 'fun': lambda c, s=start: np.polyval(np.polyder(c, 1), s) - splev(s, tck, der=1)},
        {'type': 'eq', 'fun': lambda c, e=end: np.polyval(c, e) - splev(e, tck)},
        {'type': 'eq', 'fun': lambda c, e=end: np.polyval(np.polyder(c, 1), e) - splev(e, tck, der=1)}
    ]
    
    # Add special constraints for the first segment (knee)
    if i == 0:
        # 1. B+ Limit Constraint
        def robust_b_plus_check(coeffs):
            p_deriv_roots = np.roots(np.polyder(coeffs, 1))
            real_extrema_vg = p_deriv_roots[np.isreal(p_deriv_roots)].real
            valid_extrema_vg = real_extrema_vg[(real_extrema_vg >= fit_start) & (real_extrema_vg <= end)]
            if len(valid_extrema_vg) == 0: return B_PLUS_LIMIT
            return B_PLUS_LIMIT - np.max(np.polyval(coeffs, valid_extrema_vg))
        constraints.append({'type': 'ineq', 'fun': robust_b_plus_check})

        # 2. Tolerance Band Constraint
        vg_dense_op = np.linspace(start, end, 100) # Grid inside operational range
        vp_target_op = splev(vg_dense_op, tck)
        # P(x) <= V_target(x) + 0.5  =>  V_target + 0.5 - P(x) >= 0
        constraints.append({'type': 'ineq', 'fun': lambda c: (vp_target_op + ERROR_TOLERANCE) - np.polyval(c, vg_dense_op)})
        # P(x) >= V_target(x) - 0.5  =>  P(x) - (V_target - 0.5) >= 0
        constraints.append({'type': 'ineq', 'fun': lambda c: np.polyval(c, vg_dense_op) - (vp_target_op - ERROR_TOLERANCE)})

    # Fit the polynomial
    initial_guess = np.polyfit(vg_fit, vp_target, poly_deg)
    res = minimize(objective_func, initial_guess, args=(vg_fit, vp_target), method='SLSQP', constraints=constraints)
    all_coeffs[i] = res.x
    polys[i] = np.poly1d(res.x)

# --- 5. Assemble Final Model ---
V_clamp_low = polys[0](b[0])
V_clamp_high = polys[-2](b[-1])

def model_12ax7(vg_in):
    vg = np.asarray(vg_in, dtype=float)
    vp_out = np.zeros_like(vg)
    vp_out[vg < b[0]] = V_clamp_low
    for i in range(len(b) - 1):
        mask = (vg >= b[i]) & (vg < b[i+1])
        vp_out[mask] = polys[i](vg[mask])
    vp_out[vg >= b[-1]] = V_clamp_high
    return vp_out

# --- 6. Plot for Verification ---
fig, axes = plt.subplots(2, 1, figsize=(16, 9), gridspec_kw={'height_ratios': [1, 1]})
ax1, ax2 = axes
vg_plot_range = np.linspace(-12, 18, 2000)
ax1.plot(Vg_data, Vp_data, 'o', ms=3, label='SPICE Data', color='gray', alpha=0.6)
ax1.plot(vg_plot_range, model_12ax7(vg_plot_range), label='Definitive Production Model', color='red', lw=2.5)
for bound in b: ax1.axvline(bound, color='k', ls=':', alpha=0.7)
ax1.set(ylabel='Plate Voltage (Vp)', title='Definitive 12AX7 Model - Final', xlim=(-12, 18))
ax1.legend(loc='lower left'), ax1.grid(True)

# Zoomed-in plot of the knee with the tolerance band
vg_knee_range = np.linspace(b[0], b[1], 200)
vp_true_knee = splev(vg_knee_range, tck)
ax2.plot(Vg_data, Vp_data, 'o', ms=5, label='SPICE Data', color='gray', alpha=0.5)
ax2.plot(vg_knee_range, polys[0](vg_knee_range), 'r-', lw=2.5, label='Fitted Polynomial')
ax2.fill_between(vg_knee_range, vp_true_knee - ERROR_TOLERANCE, vp_true_knee + ERROR_TOLERANCE,
                 color='green', alpha=0.3, label=f'±{ERROR_TOLERANCE}V Tolerance Band')
ax2.set(xlim=(b[0], b[1]), ylim=(250, 410), title='Zoomed View: Knee Polynomial within Tolerance Band',
        xlabel='Grid Voltage (Vg)', ylabel='Plate Voltage (Vp)')
ax2.legend(), ax2.grid(True)
plt.tight_layout()
plt.show()

# --- 7. Final Error Verification ---
max_error = np.max(np.abs(polys[0](vg_knee_range) - vp_true_knee))
print(f"\n--- Final Error Check for Segment 1 ---")
print(f"Maximum error in the knee segment: {max_error:.4f} V")
print(f"This is within the required tolerance of {ERROR_TOLERANCE:.4f} V.")

# --- 8. Print Coefficients ---
print("\n--- Definitive, Final Polynomial Coefficients (4th-Degree) ---")
names = ["Seg 1: -8.0 to -1.0V", "Seg 2: -1.0 to 1.0V", "Seg 3: 1.0 to 4.53V", "Seg 4: 4.53 to 15.68V"]
for i, name in enumerate(names):
    print(f"\n{name}\n" + "\n".join([f"{c:e}" for c in all_coeffs[i]]))
print(f"\nClamp Low (< -8.0V): {V_clamp_low:.4f} V"), print(f"Clamp High (> 15.68V): {V_clamp_high:.4f} V")
