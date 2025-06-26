import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# --- 1. Define All Constants and Target Curve ---

# Define boundaries, constraints, and tolerance
CLAMP_LOW_VG = -8.0
CLAMP_HIGH_VG = 17.62846
CLAMP_HIGH_VP = 405.0
CLAMP_LOW_VP = 15.68266
TOLERANCE = 0.5

# Define the HARD boundaries for the final model
b_seg1_end = -1.0
b_seg2_end = 1.0

# Load data
df = pd.read_csv('Draft1.txt', sep='\t', usecols=['V(grid)', 'V(plate)'])
df.rename(columns={'V(grid)': 'Vg', 'V(plate)': 'Vp'}, inplace=True)
vg_orig, vp_orig = df['Vg'].to_numpy(), df['Vp'].to_numpy()

# Hard constraints for the polynomials
P1_CONSTRAINT = (-5.340293, 400.6191)
P3_CONSTRAINT_1 = (5.253457, 37.7934)
# NEW: Interpolate from original data to get Vp for the new constraints
P3_CONSTRAINT_2 = (8.0, np.interp(8.0, vg_orig, vp_orig))
P3_CONSTRAINT_3 = (13.5, np.interp(13.5, vg_orig, vp_orig))
print("--- Added Pass-Through Constraints for Saturation Curve ---")
print(f"Constraint at Vg=+8.0V -> Vp={P3_CONSTRAINT_2[1]:.4f}V")
print(f"Constraint at Vg=+13.5V -> Vp={P3_CONSTRAINT_3[1]:.4f}V")


# Create the target curve
vp_target = vp_orig.copy()
vp_target[vg_orig <= CLAMP_LOW_VG] = CLAMP_HIGH_VP
vp_target[vg_orig >= CLAMP_HIGH_VG] = CLAMP_LOW_VP

# --- 2. Solve and Fit Polynomials with Targeted Degrees ---

# Helper functions
def poly_func(p, x): return np.polyval(p, x)
def poly_deriv_func(p, x): return np.polyval(np.polyder(p), x)
def objective_func(p, x, y): return np.sum((poly_func(p, x) - y)**2)

print("\n--- Starting Final High-Accuracy Modeling (Degrees: 4, 2, 6) ---")

# == Segment 1: Quartic (unchanged) ==
poly1_deg = 4
fit_mask1 = (vg_orig >= CLAMP_LOW_VG - 1.0) & (vg_orig <= b_seg1_end + 1.0)
vg_fit1, vp_fit1 = vg_orig[fit_mask1], vp_target[fit_mask1]
p0_1 = np.polyfit(vg_fit1, vp_fit1, poly1_deg)
constraints1 = [
    {'type': 'eq', 'fun': lambda p: poly_func(p, CLAMP_LOW_VG) - CLAMP_HIGH_VP},
    {'type': 'eq', 'fun': lambda p: poly_deriv_func(p, CLAMP_LOW_VG) - 0.0},
    {'type': 'eq', 'fun': lambda p: poly_func(p, P1_CONSTRAINT[0]) - P1_CONSTRAINT[1]}
]
res1 = minimize(objective_func, p0_1, args=(vg_fit1, vp_fit1), method='SLSQP', constraints=constraints1)
poly1 = np.poly1d(res1.x)
print(f"Segment 1 (Quartic) Fit Success: {res1.success}")

# == Segment 2: Quadratic (unchanged) ==
val_at_join1 = poly1(b_seg1_end)
deriv_at_join1 = np.polyder(poly1)(b_seg1_end)
poly2_deg = 2
fit_mask2 = (vg_orig >= b_seg1_end - 1.0) & (vg_orig <= b_seg2_end + 1.0)
vg_fit2, vp_fit2 = vg_orig[fit_mask2], vp_target[fit_mask2]
p0_2 = np.polyfit(vg_fit2, vp_fit2, poly2_deg)
constraints2 = [
    {'type': 'eq', 'fun': lambda p: poly_func(p, b_seg1_end) - val_at_join1},
    {'type': 'eq', 'fun': lambda p: poly_deriv_func(p, b_seg1_end) - deriv_at_join1}
]
res2 = minimize(objective_func, p0_2, args=(vg_fit2, vp_fit2), method='SLSQP', constraints=constraints2)
poly2 = np.poly1d(res2.x)
print(f"Segment 2 (Quadratic) Fit Success: {res2.success}")

# == Segment 3: 6th-Degree Polynomial via Direct Solver ==
print("Segment 3 (6th-Degree): Solving directly for a unique ripple-free solution...")
val_at_join2 = poly2(b_seg2_end)
deriv_at_join2 = np.polyder(poly2)(b_seg2_end)
# Define the 7 constraints
c = [
    (b_seg2_end, val_at_join2),
    (P3_CONSTRAINT_1[0], P3_CONSTRAINT_1[1]),
    (P3_CONSTRAINT_2[0], P3_CONSTRAINT_2[1]),
    (P3_CONSTRAINT_3[0], P3_CONSTRAINT_3[1]),
    (CLAMP_HIGH_VG, CLAMP_LOW_VP),
]
dc = [
    (b_seg2_end, deriv_at_join2),
    (CLAMP_HIGH_VG, 0.0)
]
# Build the A matrix (7x7) and B vector (7x1) for A*p = B
poly3_deg = 6
A = []
B = []
# Add value constraint rows
for x, y in c:
    A.append([x**(poly3_deg-i) for i in range(poly3_deg+1)])
    B.append(y)
# Add derivative constraint rows
for x, dy in dc:
    A.append([(poly3_deg-i)*x**(poly3_deg-i-1) if (poly3_deg-i)>0 else 0 for i in range(poly3_deg+1)])
    B.append(dy)
# Solve the system
try:
    p3 = np.linalg.solve(np.array(A), np.array(B))
    poly3 = np.poly1d(p3)
    print("Segment 3 Direct Solution Success.")
except np.linalg.LinAlgError:
    print("FATAL: Could not solve for Segment 3 coefficients.")
    exit()

print("\n--- Final Polynomial Coefficients ---")
print(f"Segment 1 (Deg 4):\n{poly1}\n")
print(f"Segment 2 (Deg 2):\n{poly2}\n")
print(f"Segment 3 (Deg 6):\n{poly3}\n")

# --- 3. Final Model Construction and Verification ---
def model_vp(vg):
    vg = np.asarray(vg)
    vp = np.zeros_like(vg, dtype=float)
    conds = [
        vg <= CLAMP_LOW_VG, (vg > CLAMP_LOW_VG) & (vg <= b_seg1_end),
        (vg > b_seg1_end) & (vg <= b_seg2_end), (vg > b_seg2_end) & (vg <= CLAMP_HIGH_VG),
        vg > CLAMP_HIGH_VG
    ]
    funcs = [lambda x: CLAMP_HIGH_VP, poly1, poly2, poly3, lambda x: CLAMP_LOW_VP]
    for cond, func in zip(conds, funcs): vp[cond] = func(vg[cond])
    return vp

# --- 4. Plotting and Error Analysis ---
vg_model_fine = np.linspace(-24, 24, 2000)
vp_model_fine = model_vp(vg_model_fine)
fit_range_mask = (vg_orig >= CLAMP_LOW_VG) & (vg_orig <= CLAMP_HIGH_VG)
vg_fit_range = vg_orig[fit_range_mask]
vp_target_fit_range = vp_target[fit_range_mask]
vp_model_fit_range = model_vp(vg_fit_range)
error = vp_model_fit_range - vp_target_fit_range

# Main Plot
plt.style.use('seaborn-v0_8-whitegrid')
fig1, ax1 = plt.subplots(figsize=(16, 9))
ax1.plot(vg_orig, vp_orig, 'o', color='lightgray', markersize=4, label='Original SPICE Data')
ax1.plot(vg_orig, vp_target, '-', color='royalblue', linewidth=3, alpha=0.6, label='Target Curve')
ax1.plot(vg_model_fine, vp_model_fine, '--', color='red', linewidth=2.5, label='Final Model (Degs: 4, 2, 6)')
b_list = [CLAMP_LOW_VG, b_seg1_end, b_seg2_end, CLAMP_HIGH_VG]
ax1.axvline(x=b_list[0], color='k', linestyle=':', label='Segment Boundaries')
for b in b_list[1:]: ax1.axvline(x=b, color='k', linestyle=':')
p_list = [P1_CONSTRAINT, P3_CONSTRAINT_1, P3_CONSTRAINT_2, P3_CONSTRAINT_3]
ax1.plot(p_list[0][0], p_list[0][1], 'X', color='darkorange', markersize=12, mew=2, label='Constraint Points')
for p in p_list[1:]: ax1.plot(p[0], p[1], 'X', color='darkorange', markersize=12, mew=2)
ax1.legend(fontsize=11)
ax1.set_title('Final High-Accuracy Model with Pinned Saturation Curve', fontsize=16)
ax1.set_xlabel('Grid Voltage (Vg) [V]', fontsize=12)
ax1.set_ylabel('Plate Voltage (Vp) [V]', fontsize=12)
ax1.set_xlim(-25, 25)
ax1.set_ylim(0, 425)
ax1.grid(True, which='both')
plt.tight_layout()
plt.show()

# Error Plot
fig2, ax2 = plt.subplots(figsize=(16, 9))
ax2.plot(vg_fit_range, error, 'o', color='crimson', markersize=5, label='Model Error')
ax2.axhline(y=TOLERANCE, color='green', linestyle='--', label=f'Tolerance (+/- {TOLERANCE}V)')
ax2.axhline(y=-TOLERANCE, color='green', linestyle='--')
ax2.axhline(y=0, color='black', linestyle='-', linewidth=0.7)
ax2.set_title('Model Error (Final Model vs. Target Curve)', fontsize=16)
ax2.set_xlabel('Grid Voltage (Vg) [V]', fontsize=12)
ax2.set_ylabel('Error (Vp_model - Vp_target) [V]', fontsize=12)
ax2.set_xlim(CLAMP_LOW_VG, CLAMP_HIGH_VG)
ax2.set_ylim(-TOLERANCE * 1.2, TOLERANCE * 1.2)
ax2.legend(fontsize=11)
ax2.grid(True, which='both')
print(f"\n--- Final Error Analysis (Tolerance: +/- {TOLERANCE}V) ---")
print(f"Maximum positive error: {np.max(error):.4f} V")
print(f"Maximum negative error: {np.min(error):.4f} V")
if np.all(np.abs(error) <= TOLERANCE):
    print("Success: All points are within the specified tolerance.")
else:
    print("Warning: Some points exceed the specified tolerance.")
plt.tight_layout()
plt.show()
