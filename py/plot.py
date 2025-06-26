import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import minimize

# --- 1. Data Loading and Preparation ---

# Load the data from the text file
try:
    with open('Draft1.txt', 'r') as f:
        data_content = f.read()
    
    # Use pandas to read the tab-delimited data
    from io import StringIO
    df = pd.read_csv(StringIO(data_content), sep='\t')
    
    # Rename columns as requested
    df = df.rename(columns={'V(grid)': 'Vg', 'V(plate)': 'Vp'})

    # Convert to NumPy arrays for numerical operations
    vg_data = df['Vg'].values
    vp_data = df['Vp'].values

    print("Data loaded successfully.")
    print(f"Loaded {len(df)} data points.")

except FileNotFoundError:
    print("Error: Draft1.txt not found. Please ensure the file is in the same directory.")
    exit()

# Create an interpolation function to get original Vp values at any Vg
get_original_vp = interp1d(vg_data, vp_data, bounds_error=False, fill_value="extrapolate")

# --- 2. Define Constants and Segments ---

# Polynomial degree (cubic is a good choice for curves)
POLY_DEGREE = 3
COEFF_COUNT = POLY_DEGREE + 1

# Segment boundaries
VG_START = -8.0
VG_CUTOFF_END = -1.0
VG_LINEAR_END = 1.0
VG_SAT_END = 17.62846

# Clamping values
CLAMP_HIGH_V = 405.0
# The saturation clamp value is the Vp at the end of the saturation curve
CLAMP_LOW_V = get_original_vp(VG_SAT_END) 

# Specific constraint points
CUTOFF_CONSTRAINT_VG = -5.340293
CUTOFF_CONSTRAINT_VP = 400.6191
SATURATION_CONSTRAINT_VG = 5.253457
SATURATION_CONSTRAINT_VP = 37.7934

# Overlap for fitting windows to prevent ripple at edges
# This extends the data used for fitting beyond the segment boundaries
FIT_OVERLAP = 2.0 

# --- 3. Setup for Optimization ---

# Helper function to evaluate a polynomial from a coefficient list
def eval_poly(coeffs, x):
    return np.polyval(coeffs, x)

# Helper function to evaluate the derivative of a polynomial
def eval_poly_deriv(coeffs, x):
    deriv_coeffs = np.polyder(coeffs)
    return np.polyval(deriv_coeffs, x)

# Objective function: Minimize the total squared error across all segments
def objective(p, vg_data, vp_data):
    # Unpack the 12 coefficients into 3 sets of 4
    p_cutoff = p[0:COEFF_COUNT]
    p_linear = p[COEFF_COUNT:2*COEFF_COUNT]
    p_saturation = p[2*COEFF_COUNT:3*COEFF_COUNT]

    # --- Define fitting windows (wider than segments) ---
    # Cutoff fit window
    idx_cutoff = (vg_data >= VG_START - FIT_OVERLAP) & (vg_data <= VG_CUTOFF_END + FIT_OVERLAP)
    # Linear fit window
    idx_linear = (vg_data >= VG_CUTOFF_END - FIT_OVERLAP) & (vg_data <= VG_LINEAR_END + FIT_OVERLAP)
    # Saturation fit window
    idx_saturation = (vg_data >= VG_LINEAR_END - FIT_OVERLAP) & (vg_data <= VG_SAT_END + FIT_OVERLAP)
    
    # --- Calculate errors ---
    err_cutoff = np.sum((vp_data[idx_cutoff] - eval_poly(p_cutoff, vg_data[idx_cutoff]))**2)
    err_linear = np.sum((vp_data[idx_linear] - eval_poly(p_linear, vg_data[idx_linear]))**2)
    err_saturation = np.sum((vp_data[idx_saturation] - eval_poly(p_saturation, vg_data[idx_saturation]))**2)
    
    return err_cutoff + err_linear + err_saturation

# --- 4. Define Constraints for the Optimizer ---

# Each constraint function takes the coefficient vector `p` and returns a value.
# For 'eq' type, the value must be 0. For 'ineq' type, it must be >= 0.

constraints = []

# --- A. Equality Constraints ('eq') ---

# C0 Continuity: Polynomials must meet at boundaries
# Constraint 1: Cutoff and Linear polys meet at VG_CUTOFF_END
constraints.append({'type': 'eq', 'fun': lambda p: eval_poly(p[0:COEFF_COUNT], VG_CUTOFF_END) - eval_poly(p[COEFF_COUNT:2*COEFF_COUNT], VG_CUTOFF_END)})
# Constraint 2: Linear and Saturation polys meet at VG_LINEAR_END
constraints.append({'type': 'eq', 'fun': lambda p: eval_poly(p[COEFF_COUNT:2*COEFF_COUNT], VG_LINEAR_END) - eval_poly(p[2*COEFF_COUNT:3*COEFF_COUNT], VG_LINEAR_END)})

# C1 Continuity: Derivatives must match at boundaries for smoothness
# Constraint 3: Cutoff and Linear derivatives match
constraints.append({'type': 'eq', 'fun': lambda p: eval_poly_deriv(p[0:COEFF_COUNT], VG_CUTOFF_END) - eval_poly_deriv(p[COEFF_COUNT:2*COEFF_COUNT], VG_CUTOFF_END)})
# Constraint 4: Linear and Saturation derivatives match
constraints.append({'type': 'eq', 'fun': lambda p: eval_poly_deriv(p[COEFF_COUNT:2*COEFF_COUNT], VG_LINEAR_END) - eval_poly_deriv(p[2*COEFF_COUNT:3*COEFF_COUNT], VG_LINEAR_END)})

# Boundary Value Constraints: Match original curve at segment start/end
# Constraint 5: Cutoff poly must start at the correct Vp value
constraints.append({'type': 'eq', 'fun': lambda p: eval_poly(p[0:COEFF_COUNT], VG_START) - get_original_vp(VG_START)})
# Constraint 6: Saturation poly must end at the correct Vp value
constraints.append({'type': 'eq', 'fun': lambda p: eval_poly(p[2*COEFF_COUNT:3*COEFF_COUNT], VG_SAT_END) - get_original_vp(VG_SAT_END)})

# Specific Point Constraints
# Constraint 7: Cutoff knee must pass through the specified point
constraints.append({'type': 'eq', 'fun': lambda p: eval_poly(p[0:COEFF_COUNT], CUTOFF_CONSTRAINT_VG) - CUTOFF_CONSTRAINT_VP})
# Constraint 8: Saturation knee must pass through its specified point
constraints.append({'type': 'eq', 'fun': lambda p: eval_poly(p[2*COEFF_COUNT:3*COEFF_COUNT], SATURATION_CONSTRAINT_VG) - SATURATION_CONSTRAINT_VP})


# --- B. Inequality Constraints ('ineq') ---

# Constraint: Polynomials must not exceed B+ voltage (405V)
# We check this on a fine grid of points within each segment
def check_b_plus(p):
    # Unpack coefficients
    p_cutoff = p[0:COEFF_COUNT]
    p_linear = p[COEFF_COUNT:2*COEFF_COUNT]
    p_saturation = p[2*COEFF_COUNT:3*COEFF_COUNT]
    
    # Create test grids
    vg_grid_cutoff = np.linspace(VG_START, VG_CUTOFF_END, 20)
    vg_grid_linear = np.linspace(VG_CUTOFF_END, VG_LINEAR_END, 20)
    vg_grid_saturation = np.linspace(VG_LINEAR_END, VG_SAT_END, 20)
    
    # Calculate predicted Vp
    vp_pred_cutoff = eval_poly(p_cutoff, vg_grid_cutoff)
    vp_pred_linear = eval_poly(p_linear, vg_grid_linear)
    vp_pred_saturation = eval_poly(p_saturation, vg_grid_saturation)
    
    # The constraint is CLAMP_HIGH_V - Vp >= 0.
    # We return the minimum of this value. If it's >= 0, all points are valid.
    min_diff_cutoff = np.min(CLAMP_HIGH_V - vp_pred_cutoff)
    min_diff_linear = np.min(CLAMP_HIGH_V - vp_pred_linear)
    min_diff_saturation = np.min(CLAMP_HIGH_V - vp_pred_saturation)
    
    return min(min_diff_cutoff, min_diff_linear, min_diff_saturation)

constraints.append({'type': 'ineq', 'fun': check_b_plus})


# --- 5. Run the Optimization ---

# Generate a good initial guess (p0) by doing a simple polyfit on each segment
p0_cutoff = np.polyfit(vg_data[(vg_data >= VG_START) & (vg_data < VG_CUTOFF_END)], 
                       vp_data[(vg_data >= VG_START) & (vg_data < VG_CUTOFF_END)], POLY_DEGREE)
p0_linear = np.polyfit(vg_data[(vg_data >= VG_CUTOFF_END) & (vg_data < VG_LINEAR_END)], 
                       vp_data[(vg_data >= VG_CUTOFF_END) & (vg_data < VG_LINEAR_END)], POLY_DEGREE)
p0_saturation = np.polyfit(vg_data[(vg_data >= VG_LINEAR_END) & (vg_data <= VG_SAT_END)], 
                           vp_data[(vg_data >= VG_LINEAR_END) & (vg_data <= VG_SAT_END)], POLY_DEGREE)

p0 = np.concatenate([p0_cutoff, p0_linear, p0_saturation])

print("Starting optimization... This may take a moment.")
result = minimize(objective, p0, args=(vg_data, vp_data), method='SLSQP', constraints=constraints, options={'maxiter': 1000})

if result.success:
    print("Optimization successful!")
    p_opt = result.x
    # Unpack final optimized coefficients
    p_opt_cutoff = p_opt[0:COEFF_COUNT]
    p_opt_linear = p_opt[COEFF_COUNT:2*COEFF_COUNT]
    p_opt_saturation = p_opt[2*COEFF_COUNT:3*COEFF_COUNT]
else:
    print("Optimization failed.")
    print(result.message)
    exit()

# --- 6. Create the Final Piecewise Model Function ---

def get_model_vp(vg):
    vg = np.asarray(vg)
    # Create an empty array for the results
    vp_out = np.zeros_like(vg, dtype=float)

    # Define conditions for each segment
    cond_clamp_high = vg < VG_START
    cond_cutoff = (vg >= VG_START) & (vg < VG_CUTOFF_END)
    cond_linear = (vg >= VG_CUTOFF_END) & (vg < VG_LINEAR_END)
    cond_saturation = (vg >= VG_LINEAR_END) & (vg <= VG_SAT_END)
    cond_clamp_low = vg > VG_SAT_END

    # Apply the correct function based on the condition
    vp_out[cond_clamp_high] = CLAMP_HIGH_V
    vp_out[cond_cutoff] = eval_poly(p_opt_cutoff, vg[cond_cutoff])
    vp_out[cond_linear] = eval_poly(p_opt_linear, vg[cond_linear])
    vp_out[cond_saturation] = eval_poly(p_opt_saturation, vg[cond_saturation])
    vp_out[cond_clamp_low] = CLAMP_LOW_V
    
    return vp_out


# --- 7. Visualization ---

print("Generating plot...")

# Create a dense grid of Vg points for smooth plotting
vg_plot = np.linspace(-24, 24, 2000)
vp_model_plot = get_model_vp(vg_plot)

# Create the plot
plt.style.use('seaborn-v0_8-whitegrid')
fig, ax = plt.subplots(figsize=(16, 9))

# Plot original data
ax.plot(vg_data, vp_data, 'o', color='skyblue', label='Original SPICE Data', markersize=4, alpha=0.7, zorder=1)

# Plot the final stitched model
ax.plot(vg_plot, vp_model_plot, '-', color='red', label='Stitched Polynomial Model', linewidth=2.5, zorder=2)

# Highlight the constraint points and boundaries to verify the fit
constraint_vgs = [VG_START, VG_CUTOFF_END, VG_LINEAR_END, VG_SAT_END, CUTOFF_CONSTRAINT_VG, SATURATION_CONSTRAINT_VG]
constraint_vps = get_model_vp(constraint_vgs)
ax.plot(constraint_vgs, constraint_vps, 'x', color='black', markersize=10, mew=2, label='Segment Boundaries & Constraints', zorder=3)

# Add B+ line
ax.axhline(y=CLAMP_HIGH_V, color='gray', linestyle='--', label=f'B+ Voltage Limit ({CLAMP_HIGH_V}V)')

# Formatting
ax.set_title('12AX7 Triode Vp vs. Vg Model (SPICE Data vs. Constrained Polynomial Fit)', fontsize=16)
ax.set_xlabel('Grid Voltage (Vg) [V]', fontsize=12)
ax.set_ylabel('Plate Voltage (Vp) [V]', fontsize=12)
ax.legend(fontsize=11)
ax.grid(True, which='both', linestyle='--', linewidth=0.5)

# Set plot limits to focus on the important area
ax.set_xlim(-25, 25)
ax.set_ylim(-10, CLAMP_HIGH_V + 20)

plt.tight_layout()
plt.show()
