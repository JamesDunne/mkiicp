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
poly_deg = 4 # Use more stable 4th-degree polynomials
B_PLUS_LIMIT = 405.0

# --- 3. Master Spline ---
spline_mask = (Vg_data >= -12.0) & (Vg_data <= b[-1])
tck = splrep(Vg_data[spline_mask], Vp_data[spline_mask], s=0.5)

def objective_func(c, vg, vp): return np.sum((np.polyval(c, vg) - vp)**2)
all_coeffs = [None] * len(b)
polys = [None] * len(b)

# --- 4. Independent, C1-Constrained Fit for each Segment ---
for i in range(len(b) - 1):
    start, end = b[i], b[i+1]
    
    # Use a wider fitting range for the first segment
    fit_start = start if i > 0 else -12.0
    vg_fit = np.linspace(fit_start, end, 200)
    vp_target = splev(vg_fit, tck)

    constraints = [
        {'type': 'eq', 'fun': lambda c, s=start: np.polyval(c, s) - splev(s, tck)},
        {'type': 'eq', 'fun': lambda c, s=start: np.polyval(np.polyder(c, 1), s) - splev(s, tck, der=1)},
        {'type': 'eq', 'fun': lambda c, e=end: np.polyval(c, e) - splev(e, tck)},
        {'type': 'eq', 'fun': lambda c, e=end: np.polyval(np.polyder(c, 1), e) - splev(e, tck, der=1)}
    ]
    
    # Add B+ limit constraint only for the first segment
    if i == 0:
        def robust_b_plus_check(coeffs):
            p_deriv_roots = np.roots(np.polyder(coeffs, 1))
            real_extrema_vg = p_deriv_roots[np.isreal(p_deriv_roots)].real
            valid_extrema_vg = real_extrema_vg[(real_extrema_vg >= fit_start) & (real_extrema_vg <= end)]
            if len(valid_extrema_vg) == 0: return B_PLUS_LIMIT
            return B_PLUS_LIMIT - np.max(np.polyval(coeffs, valid_extrema_vg))
        constraints.append({'type': 'ineq', 'fun': robust_b_plus_check})

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
    vp_out[vg >= b[-1]] = V_clamp_high # Handle the very last point
    return vp_out

# --- 6. Plot for Verification ---
fig, ax = plt.subplots(figsize=(15, 9))
vg_plot_range = np.linspace(-12, 18, 2000)
ax.plot(Vg_data, Vp_data, 'o', ms=3, label='SPICE Data', color='gray', alpha=0.6)
ax.plot(vg_plot_range, model_12ax7(vg_plot_range), label='Definitive Production Model', color='red', lw=2.5)
for bound in b: ax.axvline(bound, color='k', ls=':', alpha=0.7)
ax.axhline(B_PLUS_LIMIT, color='purple', ls='--', label=f'B+ Limit ({B_PLUS_LIMIT}V)')
ax.set(ylabel='Plate Voltage (Vp)', title='Definitive 12AX7 Model - Ripple-Free & C1 Continuous', xlim=(-12, 18))
ax.legend(loc='lower left'), ax.grid(True)
plt.show()

# --- 7. Final Verification Checks ---
print("--- Final Boundary Continuity Check (Value & Slope) ---")
for i in range(1, len(b)-1):
    bound = b[i]
    p_left, p_right = polys[i-1], polys[i]
    print(f"\nAt Boundary Vg = {bound:.2f}V:")
    print(f"  Vp (Left): {p_left(bound):.4f} V  | Vp (Right): {p_right(bound):.4f} V")
    print(f"  Gain (Left): {p_left.deriv(1)(bound):.4f} | Gain (Right): {p_right.deriv(1)(bound):.4f}")

# --- 8. Print Coefficients ---
print("\n--- Definitive, Final Polynomial Coefficients (4th-Degree) ---")
names = ["Seg 1: -8.0 to -1.0V", "Seg 2: -1.0 to 1.0V", "Seg 3: 1.0 to 4.53V", "Seg 4: 4.53 to 15.68V"]
for i, name in enumerate(names):
    print(f"\n{name}\n" + "\n".join([f"{c:e}" for c in all_coeffs[i]]))
print(f"\nClamp Low (< -8.0V): {V_clamp_low:.4f} V"), print(f"Clamp High (> 15.68V): {V_clamp_high:.4f} V")
