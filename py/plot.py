import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# --- 1. Data Loading and Preparation ---

# Load the data from the file
try:
    df = pd.read_csv('Draft1.txt', sep='\t')
except FileNotFoundError:
    print("Error: Draft1.txt not found. Please make sure the file is in the same directory.")
    exit()

# Rename columns and convert to numpy arrays
df.rename(columns={'V(grid)': 'Vg', 'V(plate)': 'Vp'}, inplace=True)
vg_orig = df['Vg'].to_numpy()
vp_orig = df['Vp'].to_numpy()

# --- 2. Clamping and Segmentation ---

# Define boundaries for clamping and segmentation
CLAMP_LOW_VG = -8.0
CLAMP_HIGH_VG = 17.62846
CLAMP_HIGH_VP = 405.0
CLAMP_LOW_VP = 15.68266

# Create the target Vp curve with clamping
vp_target = vp_orig.copy()
vp_target[vg_orig <= CLAMP_LOW_VG] = CLAMP_HIGH_VP
vp_target[vg_orig >= CLAMP_HIGH_VG] = CLAMP_LOW_VP

# Define segment boundaries and overlap for smoother fitting
b_seg1_end = -1.0
b_seg2_end = 1.0
overlap = 0.5  # Use 0.5V of overlap on each side of a segment for fitting

# --- 3. Constrained Polynomial Fitting with scipy.optimize.minimize ---

# Helper functions for the optimizer
def poly_func(p, x):
    """Calculates polynomial value. p are coefficients, highest power first."""
    return np.polyval(p, x)

def poly_deriv_func(p, x):
    """Calculates polynomial derivative value."""
    deriv_p = np.polyder(p)
    return np.polyval(deriv_p, x)

def objective_func(p, x, y):
    """Objective function to minimize: Sum of Squared Errors."""
    return np.sum((poly_func(p, x) - y)**2)

print("--- Starting Constrained Optimization ---")

# == Segment 1: Cut-off Knee (Vg <= -1.0) ==
# Degree: 4 (5 coefficients). Constraints: 3. DoF: 2.
poly1_deg = 4
p0_1 = np.zeros(poly1_deg + 1)
seg1_mask = (vg_orig >= CLAMP_LOW_VG) & (vg_orig <= b_seg1_end + overlap)
vg_seg1, vp_seg1 = vg_orig[seg1_mask], vp_target[seg1_mask]

# Constraints:
# 1. Value at Vg=-8.0 must be 405V
# 2. Derivative at Vg=-8.0 must be 0 (tangent)
# 3. Value at the specified constraint point
constraints1 = [
    {'type': 'eq', 'fun': lambda p: poly_func(p, CLAMP_LOW_VG) - CLAMP_HIGH_VP},
    {'type': 'eq', 'fun': lambda p: poly_deriv_func(p, CLAMP_LOW_VG) - 0.0},
    {'type': 'eq', 'fun': lambda p: poly_func(p, -5.340293) - 400.6191}
]

res1 = minimize(objective_func, p0_1, args=(vg_seg1, vp_seg1),
                method='SLSQP', constraints=constraints1)
poly1 = np.poly1d(res1.x)
print(f"Segment 1 Fit Success: {res1.success} ({res1.message})")


# == Segment 2: Critical Region (-1.0 to 1.0) ==
# Must connect smoothly to Segment 1 at Vg = -1.0
val_at_join1 = poly1(b_seg1_end)
deriv_at_join1 = np.polyder(poly1)(b_seg1_end)

# Degree: 5 (6 coefficients). Constraints: 2. DoF: 4.
poly2_deg = 5
p0_2 = np.zeros(poly2_deg + 1)
seg2_mask = (vg_orig >= b_seg1_end - overlap) & (vg_orig <= b_seg2_end + overlap)
vg_seg2, vp_seg2 = vg_orig[seg2_mask], vp_target[seg2_mask]

# Constraints:
# 1. Value at Vg=-1.0 must match end of poly1 (C0 continuity)
# 2. Derivative at Vg=-1.0 must match end of poly1 (C1 continuity)
constraints2 = [
    {'type': 'eq', 'fun': lambda p: poly_func(p, b_seg1_end) - val_at_join1},
    {'type': 'eq', 'fun': lambda p: poly_deriv_func(p, b_seg1_end) - deriv_at_join1}
]

res2 = minimize(objective_func, p0_2, args=(vg_seg2, vp_seg2),
                method='SLSQP', constraints=constraints2)
poly2 = np.poly1d(res2.x)
print(f"Segment 2 Fit Success: {res2.success} ({res2.message})")


# == Segment 3: Saturation Knee (Vg >= 1.0) ==
# Must connect smoothly to Segment 2 at Vg = 1.0
val_at_join2 = poly2(b_seg2_end)
deriv_at_join2 = np.polyder(poly2)(b_seg2_end)

# Degree: 6 (7 coefficients). Constraints: 5. DoF: 2.
poly3_deg = 6
p0_3 = np.zeros(poly3_deg + 1)
seg3_mask = (vg_orig >= b_seg2_end - overlap) & (vg_orig <= CLAMP_HIGH_VG)
vg_seg3, vp_seg3 = vg_orig[seg3_mask], vp_target[seg3_mask]

# Constraints:
# 1. Value at Vg=1.0 must match end of poly2 (C0)
# 2. Derivative at Vg=1.0 must match end of poly2 (C1)
# 3. Value at the specified constraint point
# 4. Value at Vg=17.6... must be clamp value
# 5. Derivative at Vg=17.6... must be 0 (tangent)
constraints3 = [
    {'type': 'eq', 'fun': lambda p: poly_func(p, b_seg2_end) - val_at_join2},
    {'type': 'eq', 'fun': lambda p: poly_deriv_func(p, b_seg2_end) - deriv_at_join2},
    {'type': 'eq', 'fun': lambda p: poly_func(p, 5.253457) - 37.7934},
    {'type': 'eq', 'fun': lambda p: poly_func(p, CLAMP_HIGH_VG) - CLAMP_LOW_VP},
    {'type': 'eq', 'fun': lambda p: poly_deriv_func(p, CLAMP_HIGH_VG) - 0.0}
]

res3 = minimize(objective_func, p0_3, args=(vg_seg3, vp_seg3),
                method='SLSQP', constraints=constraints3, options={'maxiter': 500})
poly3 = np.poly1d(res3.x)
print(f"Segment 3 Fit Success: {res3.success} ({res3.message})")


print("\n--- Polynomial Coefficients (highest power first) ---")
print(f"Segment 1 (Cut-off, Vg in [-8.0, -1.0]):\n{poly1}\n")
print(f"Segment 2 (Critical, Vg in [-1.0, 1.0]):\n{poly2}\n")
print(f"Segment 3 (Saturation, Vg in [1.0, 17.6]):\n{poly3}\n")

# --- 4. Final Model Construction and Verification ---

def model_vp(vg):
    vg = np.asarray(vg)
    vp = np.zeros_like(vg, dtype=float)

    conds = [
        vg <= CLAMP_LOW_VG,
        (vg > CLAMP_LOW_VG) & (vg <= b_seg1_end),
        (vg > b_seg1_end) & (vg <= b_seg2_end),
        (vg > b_seg2_end) & (vg <= CLAMP_HIGH_VG),
        vg > CLAMP_HIGH_VG
    ]
    # Note the change in boundary names
    funcs = [
        lambda x: CLAMP_HIGH_VP,
        lambda x: poly1(x),
        lambda x: poly2(x),
        lambda x: poly3(x),
        lambda x: CLAMP_LOW_VP
    ]
    
    for cond, func in zip(conds, funcs):
        # Clip the polynomial outputs just in case of numerical instability
        if func in [lambda x: poly1(x), lambda x: poly2(x), lambda x: poly3(x)]:
             vp[cond] = np.clip(func(vg[cond]), CLAMP_LOW_VP, CLAMP_HIGH_VP)
        else:
             vp[cond] = func(vg[cond])
             
    return vp

vg_model_fine = np.linspace(-24, 24, 2000)
vp_model_fine = model_vp(vg_model_fine)

fit_range_mask = (vg_orig >= CLAMP_LOW_VG) & (vg_orig <= CLAMP_HIGH_VG)
vg_fit_range = vg_orig[fit_range_mask]
vp_target_fit_range = vp_target[fit_range_mask]
vp_model_fit_range = model_vp(vg_fit_range)
error = vp_model_fit_range - vp_target_fit_range

# --- 5. Visualization ---

plt.style.use('seaborn-v0_8-whitegrid')
fig1, ax1 = plt.subplots(figsize=(16, 9))

# Plot data
ax1.plot(vg_orig, vp_orig, 'o', color='lightgray', markersize=4, label='Original SPICE Data')
ax1.plot(vg_orig, vp_target, '-', color='royalblue', linewidth=3, alpha=0.6, label='Target Curve (with Clamping)')
ax1.plot(vg_model_fine, vp_model_fine, '--', color='red', linewidth=2.5, label='C1-Continuous Polynomial Model')

# Mark boundaries and constraints
ax1.axvline(x=CLAMP_LOW_VG, color='k', linestyle=':', label='Segment/Clamp Boundaries')
ax1.axvline(x=b_seg1_end, color='k', linestyle=':')
ax1.axvline(x=b_seg2_end, color='k', linestyle=':')
ax1.axvline(x=CLAMP_HIGH_VG, color='k', linestyle=':')
ax1.plot(-5.340293, 400.6191, 'X', color='darkorange', markersize=12, mew=2, label='Constraint Point (Cut-off)')
ax1.plot(5.253457, 37.7934, 'P', color='purple', markersize=12, mew=2, label='Constraint Point (Saturation)')

# Add annotations for tangency points
ax1.annotate('Tangent to 405V', xy=(CLAMP_LOW_VG, CLAMP_HIGH_VP), xytext=(-15, 380),
             arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=8),
             fontsize=12, ha='center')
ax1.annotate(f'Tangent to {CLAMP_LOW_VP:.2f}V', xy=(CLAMP_HIGH_VG, CLAMP_LOW_VP), xytext=(12, 50),
             arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=8),
             fontsize=12, ha='center')


ax1.set_title('12AX7 Triode Model: Vp vs. Vg with C1-Continuous Piecewise Polynomial Fit', fontsize=16)
ax1.set_xlabel('Grid Voltage (Vg) [V]', fontsize=12)
ax1.set_ylabel('Plate Voltage (Vp) [V]', fontsize=12)
ax1.legend(fontsize=11)
ax1.set_xlim(-25, 25)
ax1.set_ylim(0, 425)
ax1.grid(True, which='both')

plt.tight_layout()
plt.show()

# Plot 2: Error Analysis
fig2, ax2 = plt.subplots(figsize=(16, 9))
ax2.plot(vg_fit_range, error, 'o', color='crimson', markersize=5, label='Model Error')
ax2.axhline(y=0.25, color='green', linestyle='--', label='Tolerance (+/- 0.25V)')
ax2.axhline(y=-0.25, color='green', linestyle='--')
ax2.axhline(y=0, color='black', linestyle='-', linewidth=0.7)
ax2.set_title('Model Error (Polynomial Model vs. Target Curve)', fontsize=16)
ax2.set_xlabel('Grid Voltage (Vg) [V]', fontsize=12)
ax2.set_ylabel('Error (Vp_model - Vp_target) [V]', fontsize=12)
ax2.set_xlim(CLAMP_LOW_VG, CLAMP_HIGH_VG)
ax2.set_ylim(-0.3, 0.3)
ax2.legend(fontsize=11)
ax2.grid(True, which='both')

print(f"\n--- Error Analysis ---")
print(f"Maximum positive error: {np.max(error):.4f} V")
print(f"Maximum negative error: {np.min(error):.4f} V")
if np.all(np.abs(error) <= 0.25):
    print("Success: All points are within the +/- 0.25V tolerance.")
else:
    print("Warning: Some points exceed the +/- 0.25V tolerance.")

plt.tight_layout()
plt.show()
