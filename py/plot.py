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
ERROR_TOLERANCE = 0.5

# --- 3. Master Spline ---
spline_mask = (Vg_data >= -12.0) & (Vg_data <= b[-1])
tck = splrep(Vg_data[spline_mask], Vp_data[spline_mask], s=0.5)

def objective_func(c, vg, vp): return np.sum((np.polyval(c, vg) - vp)**2)
all_coeffs = [None] * (len(b) - 1)
polys = [None] * (len(b) - 1)

# --- 4. Sequential, Daisy-Chained Constrained Fit ---

# Stage 1: Fit Segment 2 (Central Bridge) first
i=1
start, end = b[i], b[i+1]
vg_fit = np.linspace(start, end, 200)
vp_target = splev(vg_fit, tck)
constraints = [
    {'type': 'eq', 'fun': lambda c: np.polyval(c, start) - splev(start, tck)},
    {'type': 'eq', 'fun': lambda c: np.polyval(np.polyder(c, 1), start) - splev(start, tck, der=1)},
    {'type': 'eq', 'fun': lambda c: np.polyval(c, end) - splev(end, tck)},
    {'type': 'eq', 'fun': lambda c: np.polyval(np.polyder(c, 1), end) - splev(end, tck, der=1)}
]
all_coeffs[i] = minimize(objective_func, np.polyfit(vg_fit, vp_target, poly_deg), args=(vg_fit, vp_target), method='SLSQP', constraints=constraints).x
polys[i] = np.poly1d(all_coeffs[i])

# Stage 2: Fit remaining segments, chaining to the previous one
for i in [2, 3, 0]: # Order: S3, S4, then S1
    if i == 0: # Segment 1 (Knee)
        start, end = b[i], b[i+1]
        fit_start = -12.0
        p_next = polys[i+1]
        constraints = [
            {'type': 'eq', 'fun': lambda c: np.polyval(c, end) - p_next(end)},
            {'type': 'eq', 'fun': lambda c: np.polyval(np.polyder(c, 1), end) - p_next.deriv(1)(end)},
            {'type': 'eq', 'fun': lambda c: np.polyval(c, start) - splev(start, tck)},
        ]
        vg_dense_op = np.linspace(start, end, 100)
        vp_target_op = splev(vg_dense_op, tck)
        constraints.append({'type': 'ineq', 'fun': lambda c: (vp_target_op + ERROR_TOLERANCE) - np.polyval(c, vg_dense_op)})
        constraints.append({'type': 'ineq', 'fun': lambda c: np.polyval(c, vg_dense_op) - (vp_target_op - ERROR_TOLERANCE)})
        
        # --- Foolproof B+ Check ---
        def robust_b_plus_check(coeffs, seg_start, seg_end):
            p_deriv_roots = np.roots(np.polyder(coeffs, 1))
            extrema_candidates = np.concatenate(([seg_start, seg_end], p_deriv_roots[np.isreal(p_deriv_roots)].real))
            valid_extrema_vg = extrema_candidates[(extrema_candidates >= seg_start) & (extrema_candidates <= seg_end)]
            return B_PLUS_LIMIT - np.max(np.polyval(coeffs, valid_extrema_vg))
        constraints.append({'type': 'ineq', 'fun': lambda c: robust_b_plus_check(c, fit_start, end)})
    else: # Segments 3 and 4
        start, end = b[i], b[i+1]
        fit_start = start
        p_prev = polys[i-1]
        constraints = [
            {'type': 'eq', 'fun': lambda c: np.polyval(c, start) - p_prev(start)},
            {'type': 'eq', 'fun': lambda c: np.polyval(np.polyder(c, 1), start) - p_prev.deriv(1)(start)},
            {'type': 'eq', 'fun': lambda c: np.polyval(c, end) - splev(end, tck)},
        ]
        if i == 3: # Add zero-slope constraint for the final segment
             constraints.append({'type': 'eq', 'fun': lambda c: np.polyval(np.polyder(c, 1), end)})
        else:
             constraints.append({'type': 'eq', 'fun': lambda c: np.polyval(np.polyder(c, 1), end) - splev(end, tck, der=1)})

    vg_fit = np.linspace(fit_start, end, 200)
    vp_target = splev(vg_fit, tck)
    all_coeffs[i] = minimize(objective_func, np.polyfit(vg_fit, vp_target, poly_deg), args=(vg_fit, vp_target), method='SLSQP', constraints=constraints).x
    polys[i] = np.poly1d(all_coeffs[i])

# --- 5. Assemble Final Model ---
V_clamp_low = polys[0](b[0])
V_clamp_high = polys[-1](b[-1])
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
fig, axes = plt.subplots(2, 1, figsize=(15, 12), gridspec_kw={'height_ratios': [1, 1]})
ax1, ax2 = axes
vg_plot_range = np.linspace(-12, 18, 2000)
ax1.plot(Vg_data, Vp_data, 'o', ms=3, label='SPICE Data', color='gray', alpha=0.6)
ax1.plot(vg_plot_range, model_12ax7(vg_plot_range), label='Definitive Production Model', color='red', lw=2.5)
for bound in b: ax1.axvline(bound, color='k', ls=':', alpha=0.7)
ax1.set(ylabel='Plate Voltage (Vp)', title='Definitive 12AX7 Model - Final', xlim=(-12, 18))
ax1.legend(loc='lower left'), ax1.grid(True)

vg_knee_range = np.linspace(b[0], b[1], 200)
vp_true_knee = splev(vg_knee_range, tck)
ax2.plot(Vg_data, Vp_data, 'o', ms=5, label='SPICE Data', color='gray', alpha=0.5)
ax2.plot(vg_knee_range, polys[0](vg_knee_range), 'r-', lw=2.5, label='Fitted Polynomial')
ax2.fill_between(vg_knee_range, vp_true_knee - ERROR_TOLERANCE, vp_true_knee + ERROR_TOLERANCE,
                 color='green', alpha=0.3, label=f'±{ERROR_TOLERANCE}V Tolerance Band')
ax2.axhline(B_PLUS_LIMIT, color='purple', ls='--', label='B+ Limit')
ax2.set(xlim=(b[0], b[1]), ylim=(250, 406), title='Zoomed View: Knee is B+ Safe and within Tolerance',
        xlabel='Grid Voltage (Vg)', ylabel='Plate Voltage (Vp)')
ax2.legend(), ax2.grid(True)
plt.tight_layout()
plt.show()

# --- 7. Final Verification Checks ---
print("--- Final Boundary Continuity Check (Value & Slope) ---")
for i in range(1, len(b)-1):
    bound = b[i]
    p_left, p_right = polys[i-1], polys[i]
    print(f"\nAt Boundary Vg = {bound:.2f}V:")
    print(f"  Vp (Left): {p_left(bound):.4f} V  | Vp (Right): {p_right(bound):.4f} V  -> Match: {np.isclose(p_left(bound), p_right(bound))}")
    print(f"  Gain (Left): {p_left.deriv(1)(bound):.4f} | Gain (Right): {p_right.deriv(1)(bound):.4f} -> Match: {np.isclose(p_left.deriv(1)(bound), p_right.deriv(1)(bound))}")
max_error = np.max(np.abs(polys[0](vg_knee_range) - vp_true_knee))
print(f"\nMaximum error in the knee segment: {max_error:.4f} V (Required: <= {ERROR_TOLERANCE:.1f} V)")

# --- 8. Print Coefficients ---
print("\n--- Definitive, Final Polynomial Coefficients (4th-Degree) ---")
names = ["Seg 1: -8.0 to -1.0V", "Seg 2: -1.0 to 1.0V", "Seg 3: 1.0 to 4.53V", "Seg 4: 4.53 to 15.68V"]
for i, name in enumerate(names):
    print(f"\n{name}\n" + "\n".join([f"{c:e}" for c in all_coeffs[i]]))
print(f"\nClamp Low (< -8.0V): {V_clamp_low:.4f} V"), print(f"Clamp High (> 15.68V): {V_clamp_high:.4f} V")
