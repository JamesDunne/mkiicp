import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

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

# Define segment boundaries
b1 = CLAMP_LOW_VG  # -8.0V
b2 = -1.0
b3 = 1.0
b4 = CLAMP_HIGH_VG # 17.62846V

# Filter data for each segment from the target curve
seg1_mask = (vg_orig >= b1) & (vg_orig <= b2)
vg_seg1, vp_seg1 = vg_orig[seg1_mask], vp_target[seg1_mask]

seg2_mask = (vg_orig >= b2) & (vg_orig <= b3)
vg_seg2, vp_seg2 = vg_orig[seg2_mask], vp_target[seg2_mask]

seg3_mask = (vg_orig >= b3) & (vg_orig <= b4)
vg_seg3, vp_seg3 = vg_orig[seg3_mask], vp_target[seg3_mask]


# --- 3. Constrained Polynomial Fitting ---

# Constraint points from the prompt
p1_constraint = (-5.340293, 400.6191)
p3_given_constraint = (5.253457, 37.7934)

# == Segment 1: Cut-off Knee (-8V to -1V) ==
# Polynomial P(x) = (x - xc) * (a*x^2 + b*x + c) + yc
# This forces P(xc) = yc. We fit for a, b, c.
def constrained_poly_1(x, a, b, c):
    xc, yc = p1_constraint
    return (x - xc) * (a * x**2 + b * x + c) + yc

popt1, _ = curve_fit(constrained_poly_1, vg_seg1, vp_seg1)
a1, b1_fit, c1 = popt1
poly1 = np.poly1d(np.polyadd(np.polymul([a1, b1_fit, c1], [1, -p1_constraint[0]]), [0, 0, 0, p1_constraint[1]]))

# == Segment 2: Critical Region (-1V to +1V) ==
# Ensure continuity: The start of poly2 must match the end of poly1
p2_constraint = (b2, poly1(b2))

# Using a higher degree (5th) for better fit in the critical region
# P(x) = (x - xc) * (a*x^4 + b*x^3 + c*x^2 + d*x + e) + yc
def constrained_poly_2(x, a, b, c, d, e):
    xc, yc = p2_constraint
    return (x - xc) * (a * x**4 + b * x**3 + c * x**2 + d * x + e) + yc

popt2, _ = curve_fit(constrained_poly_2, vg_seg2, vp_seg2, maxfev=10000)
a2, b2_fit, c2, d2, e2 = popt2
poly2 = np.poly1d(np.polyadd(np.polymul([a2, b2_fit, c2, d2, e2], [1, -p2_constraint[0]]), [0, 0, 0, 0, 0, p2_constraint[1]]))

# == Segment 3: Saturation Knee (+1V to +17.6V) ==
# This segment has TWO constraints:
# 1. Continuity with Segment 2: (b3, poly2(b3))
# 2. The point given in the prompt.
p3_constraint1 = (b3, poly2(b3))
p3_constraint2 = p3_given_constraint

# Define the line L(x) that passes through the two constraint points
x1, y1 = p3_constraint1
x2, y2 = p3_constraint2
line_slope = (y2 - y1) / (x2 - x1)
line_intercept = y1 - line_slope * x1
line_poly = np.poly1d([line_slope, line_intercept])

# P(x) = (x - x1)*(x - x2)*(a*x^2 + b*x + c) + L(x)
# This forces P(x1)=y1 and P(x2)=y2. We fit for a, b, c.
def constrained_poly_3(x, a, b, c):
    # R(x) is a 2nd degree polynomial, so P(x) is 4th degree
    r_poly = a * x**2 + b * x + c
    return (x - x1) * (x - x2) * r_poly + line_poly(x)

popt3, _ = curve_fit(constrained_poly_3, vg_seg3, vp_seg3, maxfev=10000)
a3, b3_fit, c3 = popt3
r_poly3 = np.poly1d([a3, b3_fit, c3])
poly3 = np.polyadd(np.polymul(r_poly3, [1, -(x1 + x2), x1 * x2]), line_poly)

print("--- Polynomial Coefficients (highest power first) ---")
print(f"Segment 1 (Cut-off, Vg in [-8.0, -1.0]):\n{poly1}\n")
print(f"Segment 2 (Critical, Vg in [-1.0, 1.0]):\n{poly2}\n")
print(f"Segment 3 (Saturation, Vg in [1.0, 17.6]):\n{poly3}\n")


# --- 4. Final Model Construction and Verification ---

def model_vp(vg):
    """Applies the piecewise polynomial model to a Vg value or array."""
    vg = np.asarray(vg)
    vp = np.zeros_like(vg, dtype=float)

    # Apply clamping and polynomials based on Vg value
    conds = [
        vg <= b1,
        (vg > b1) & (vg <= b2),
        (vg > b2) & (vg <= b3),
        (vg > b3) & (vg <= b4),
        vg > b4
    ]
    funcs = [
        lambda x: CLAMP_HIGH_VP,
        lambda x: poly1(x),
        lambda x: poly2(x),
        lambda x: poly3(x),
        lambda x: CLAMP_LOW_VP
    ]
    
    for cond, func in zip(conds, funcs):
        vp[cond] = func(vg[cond])

    # Final check on the physical B+ limit
    vp[vp > CLAMP_HIGH_VP] = CLAMP_HIGH_VP

    return vp

# Generate data for plotting the model
vg_model_fine = np.linspace(-24, 24, 2000)
vp_model_fine = model_vp(vg_model_fine)

# Calculate error within the fitted range
fit_range_mask = (vg_orig >= b1) & (vg_orig <= b4)
vg_fit_range = vg_orig[fit_range_mask]
vp_target_fit_range = vp_target[fit_range_mask]
vp_model_fit_range = model_vp(vg_fit_range)

error = vp_model_fit_range - vp_target_fit_range


# --- 5. Visualization ---

plt.style.use('seaborn-v0_8-whitegrid')

# Plot 1: Main Curve Fit
fig1, ax1 = plt.subplots(figsize=(16, 9))

# Plot the raw and target data
ax1.plot(vg_orig, vp_orig, 'o', color='lightgray', markersize=4, label='Original SPICE Data')
ax1.plot(vg_orig, vp_target, '-', color='royalblue', linewidth=2.5, label='Target Curve (with Clamping)')

# Plot the piecewise polynomial model
ax1.plot(vg_model_fine, vp_model_fine, '--', color='red', linewidth=2, label='Piecewise Polynomial Model')

# Mark segment boundaries
ax1.axvline(x=b1, color='k', linestyle=':', label='Segment Boundaries')
ax1.axvline(x=b2, color='k', linestyle=':')
ax1.axvline(x=b3, color='k', linestyle=':')
ax1.axvline(x=b4, color='k', linestyle=':')

# Mark constraint points
ax1.plot(*p1_constraint, 'X', color='darkorange', markersize=10, mew=2, label='Constraint Point (Cut-off)')
ax1.plot(*p3_given_constraint, 'P', color='purple', markersize=10, mew=2, label='Constraint Point (Saturation)')

ax1.set_title('12AX7 Triode Model: Vp vs. Vg with Piecewise Polynomial Fit', fontsize=16)
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

# Plot tolerance bands
ax2.axhline(y=0.25, color='green', linestyle='--', label='Tolerance (+/- 0.25V)')
ax2.axhline(y=-0.25, color='green', linestyle='--')
ax2.axhline(y=0, color='black', linestyle='-', linewidth=0.7)


ax2.set_title('Model Error (Polynomial Model vs. Target Curve)', fontsize=16)
ax2.set_xlabel('Grid Voltage (Vg) [V]', fontsize=12)
ax2.set_ylabel('Error (Vp_model - Vp_target) [V]', fontsize=12)
ax2.set_xlim(b1, b4)
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
