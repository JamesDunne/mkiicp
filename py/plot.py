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

# Hard constraints for the polynomials
P1_CONSTRAINT_1 = (CLAMP_LOW_VG, CLAMP_HIGH_VP) # Value at clamp
P1_CONSTRAINT_2 = (CLAMP_LOW_VG, 0.0)           # Derivative at clamp
P1_CONSTRAINT_3 = (-5.340293, 400.6191)
P1_CONSTRAINT_4 = (-1.027793, 252.1438) # New constraint, defines end of segment 1

P3_CONSTRAINT_1 = (5.253457, 37.7934)

# New segment boundaries based on constraints
b_seg1_end = P1_CONSTRAINT_4[0] # -1.027793V
b_seg2_end = 1.0

# Load data and create the target curve
df = pd.read_csv('Draft1.txt', sep='\t', usecols=['V(grid)', 'V(plate)'])
df.rename(columns={'V(grid)': 'Vg', 'V(plate)': 'Vp'}, inplace=True)
vg_orig, vp_orig = df['Vg'].to_numpy(), df['Vp'].to_numpy()
vp_target = vp_orig.copy()
vp_target[vg_orig <= CLAMP_LOW_VG] = CLAMP_HIGH_VP
vp_target[vg_orig >= CLAMP_HIGH_VG] = CLAMP_LOW_VP

# --- 2. Solve and Fit Polynomials with New Strategy ---

# Helper functions
def poly_func(p, x): return np.polyval(p, x)
def poly_deriv_func(p, x): return np.polyval(np.polyder(p), x)
def objective_func(p, x, y): return np.sum((poly_func(p, x) - y)**2)

print("--- Starting Final Robust Modeling ---")

# == Segment 1: Direct Solution for a Cubic Polynomial ==
# A cubic (ax^3+bx^2+cx+d) is uniquely defined by 4 constraints. We solve A*p=B.
print("Segment 1: Solving directly for a unique cubic polynomial...")
x1, y1 = P1_CONSTRAINT_1
x2, dy2 = P1_CONSTRAINT_2
x3, y3 = P1_CONSTRAINT_3
x4, y4 = P1_CONSTRAINT_4

A = np.array([
    [x1**3, x1**2, x1, 1],  # P(x1) = y1
    [3*x2**2, 2*x2, 1, 0],  # P'(x2) = dy2
    [x3**3, x3**2, x3, 1],  # P(x3) = y3
    [x4**3, x4**2, x4, 1]   # P(x4) = y4
])
B = np.array([y1, dy2, y3, y4])

try:
    p1 = np.linalg.solve(A, B)
    poly1 = np.poly1d(p1)
    # Verification: Check if it exceeds the physical limit
    vg_fine_seg1 = np.linspace(CLAMP_LOW_VG, b_seg1_end, 100)
    max_vp_seg1 = np.max(poly1(vg_fine_seg1))
    print(f"Segment 1 solved. Max Vp in range: {max_vp_seg1:.4f}V (Limit: {CLAMP_HIGH_VP}V)")
    if max_vp_seg1 > CLAMP_HIGH_VP + 1e-6: # Add tolerance for float error
        print("WARNING: Segment 1 slightly exceeds the 405V physical limit.")
except np.linalg.LinAlgError:
    print("FATAL: Could not solve for Segment 1 coefficients. The constraints may be ill-conditioned.")
    exit()

# == Segment 2: Constrained Fit for a Cubic Polynomial (No Ripple) ==
print("Segment 2: Fitting a stable cubic polynomial...")
val_at_join1 = poly1(b_seg1_end)
deriv_at_join1 = np.polyder(poly1)(b_seg1_end)
poly2_deg = 3
seg2_mask = (vg_orig >= b_seg1_end) & (vg_orig <= b_seg2_end)
vg_seg2, vp_seg2 = vg_orig[seg2_mask], vp_target[seg2_mask]
p0_2 = np.polyfit(vg_seg2, vp_seg2, poly2_deg) # Smart initial guess
constraints2 = [
    {'type': 'eq', 'fun': lambda p: poly_func(p, b_seg1_end) - val_at_join1},
    {'type': 'eq', 'fun': lambda p: poly_deriv_func(p, b_seg1_end) - deriv_at_join1}
]
res2 = minimize(objective_func, p0_2, args=(vg_seg2, vp_seg2), method='SLSQP', constraints=constraints2)
poly2 = np.poly1d(res2.x)

# == Segment 3: Constrained Fit for a Quintic Polynomial ==
print("Segment 3: Fitting a quintic polynomial for the saturation curve...")
val_at_join2 = poly2(b_seg2_end)
deriv_at_join2 = np.polyder(poly2)(b_seg2_end)
poly3_deg = 5
seg3_mask = (vg_orig >= b_seg2_end) & (vg_orig <= CLAMP_HIGH_VG)
vg_seg3, vp_seg3 = vg_orig[seg3_mask], vp_target[seg3_mask]
p0_3 = np.polyfit(vg_seg3, vp_seg3, poly3_deg) # Smart initial guess
constraints3 = [
    {'type': 'eq', 'fun': lambda p: poly_func(p, b_seg2_end) - val_at_join2},
    {'type': 'eq', 'fun': lambda p: poly_deriv_func(p, b_seg2_end) - deriv_at_join2},
    {'type': 'eq', 'fun': lambda p: poly_func(p, P3_CONSTRAINT_1[0]) - P3_CONSTRAINT_1[1]},
    {'type': 'eq', 'fun': lambda p: poly_func(p, CLAMP_HIGH_VG) - CLAMP_LOW_VP},
    {'type': 'eq', 'fun': lambda p: poly_deriv_func(p, CLAMP_HIGH_VG) - 0.0}
]
res3 = minimize(objective_func, p0_3, args=(vg_seg3, vp_seg3), method='SLSQP', constraints=constraints3)
poly3 = np.poly1d(res3.x)

print("\n--- Final Polynomial Coefficients (highest power first) ---")
print(f"Segment 1 (Deg 3, Vg in [-8.0, {b_seg1_end}]):\n{poly1}\n")
print(f"Segment 2 (Deg 3, Vg in [{b_seg1_end}, {b_seg2_end}]):\n{poly2}\n")
print(f"Segment 3 (Deg 5, Vg in [{b_seg2_end}, {CLAMP_HIGH_VG}]):\n{poly3}\n")

# --- 3. Final Model Construction and Verification ---

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
    funcs = [lambda x: CLAMP_HIGH_VP, poly1, poly2, poly3, lambda x: CLAMP_LOW_VP]
    for cond, func in zip(conds, funcs):
        vp[cond] = func(vg[cond])
    return vp

vg_model_fine = np.linspace(-24, 24, 2000)
vp_model_fine = model_vp(vg_model_fine)
fit_range_mask = (vg_orig >= CLAMP_LOW_VG) & (vg_orig <= CLAMP_HIGH_VG)
vg_fit_range = vg_orig[fit_range_mask]
vp_target_fit_range = vp_target[fit_range_mask]
vp_model_fit_range = model_vp(vg_fit_range)
error = vp_model_fit_range - vp_target_fit_range

# --- 4. Visualization ---

plt.style.use('seaborn-v0_8-whitegrid')
fig1, ax1 = plt.subplots(figsize=(16, 9))
ax1.plot(vg_orig, vp_orig, 'o', color='lightgray', markersize=4, label='Original SPICE Data')
ax1.plot(vg_orig, vp_target, '-', color='royalblue', linewidth=3, alpha=0.6, label='Target Curve (with Clamping)')
ax1.plot(vg_model_fine, vp_model_fine, '--', color='red', linewidth=2.5, label='Final Model (Degs: 3, 3, 5)')
# Mark boundaries and constraints
b_list = [CLAMP_LOW_VG, b_seg1_end, b_seg2_end, CLAMP_HIGH_VG]
ax1.axvline(x=b_list[0], color='k', linestyle=':', label='Segment Boundaries')
for b in b_list[1:]: ax1.axvline(x=b, color='k', linestyle=':')
p_list = [P1_CONSTRAINT_3, P1_CONSTRAINT_4, P3_CONSTRAINT_1]
for p in p_list: ax1.plot(p[0], p[1], 'X', color='darkorange', markersize=12, mew=2)
ax1.get_legend_handles_labels()[0][-1] = 'Constraint Points' # Relabel for clarity
ax1.legend()
ax1.set_title('12AX7 Triode Model: Final, Stable, and Accurate Fit', fontsize=16)
ax1.set_xlabel('Grid Voltage (Vg) [V]', fontsize=12)
ax1.set_ylabel('Plate Voltage (Vp) [V]', fontsize=12)
ax1.set_xlim(-25, 25)
ax1.set_ylim(0, 425)
ax1.grid(True, which='both')
plt.tight_layout()
plt.show()

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
