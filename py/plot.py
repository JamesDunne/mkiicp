import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import minimize

# --- 1. Data Loading and Preparation ---
try:
    with open('Draft1.txt', 'r') as f:
        data_content = f.read()
    from io import StringIO
    df = pd.read_csv(StringIO(data_content), sep='\t')
    df = df.rename(columns={'V(grid)': 'Vg', 'V(plate)': 'Vp'})
    vg_data = df['Vg'].values
    vp_data = df['Vp'].values
    print("Data loaded successfully.")
    print(f"Loaded {len(df)} data points.")
except FileNotFoundError:
    print("Error: Draft1.txt not found. Please ensure the file is in the same directory.")
    exit()

get_original_vp = interp1d(vg_data, vp_data, bounds_error=False, fill_value="extrapolate")

# --- 2. Define Constants and Segments ---

# *** FINAL 5-SEGMENT STRATEGY ***
# Two cubics for the sharp cutoff, three quadratics for the rest.
DEGREES = [3, 3, 2, 2, 2] 
COEFF_COUNTS = [d + 1 for d in DEGREES]
# Calculate indices for slicing the master coefficient vector 'p'
p_indices = np.cumsum([0] + COEFF_COUNTS)

# Segment boundaries - added a new boundary at -5.3V
BOUNDARIES = [-8.0, -5.3, -1.0, 1.0, 6.0, 17.62846]

# Clamping values
CLAMP_HIGH_V = 405.0
CLAMP_LOW_V = get_original_vp(BOUNDARIES[-1])

# --- 3. Setup for Optimization ---
def eval_poly(coeffs, x):
    return np.polyval(coeffs, x)

def eval_poly_deriv(coeffs, x):
    return np.polyval(np.polyder(coeffs), x)

def objective(p, vg_data, vp_data):
    total_error = 0
    # Loop through each segment to calculate its contribution to the error
    for i in range(len(DEGREES)):
        p_segment = p[p_indices[i]:p_indices[i+1]]
        vg_start, vg_end = BOUNDARIES[i], BOUNDARIES[i+1]
        
        # Use a wider window for fitting to ensure smooth transitions
        idx_segment_fit = (vg_data >= vg_start - 1.0) & (vg_data <= vg_end + 1.0)

        if np.any(idx_segment_fit):
            error = np.sum((vp_data[idx_segment_fit] - eval_poly(p_segment, vg_data[idx_segment_fit]))**2)
            total_error += error
            
    return total_error

# --- 4. Define Constraints for the Optimizer ---
constraints = []

# C0 and C1 Continuity at all internal junctions
for i in range(len(DEGREES) - 1):
    junction_vg = BOUNDARIES[i+1]
    
    # Define lambdas to capture the correct indices for each iteration
    get_p_before = (lambda i_cap: lambda p: p[p_indices[i_cap]:p_indices[i_cap+1]])(i)
    get_p_after = (lambda i_cap: lambda p: p[p_indices[i_cap+1]:p_indices[i_cap+2]])(i)
    
    # C0 constraint: values must be equal
    c0_constraint = {'type': 'eq', 'fun': lambda p, v=junction_vg, pb=get_p_before, pa=get_p_after: eval_poly(pb(p), v) - eval_poly(pa(p), v)}
    constraints.append(c0_constraint)
    
    # C1 constraint: derivatives must be equal
    c1_constraint = {'type': 'eq', 'fun': lambda p, v=junction_vg, pb=get_p_before, pa=get_p_after: eval_poly_deriv(pb(p), v) - eval_poly_deriv(pa(p), v)}
    constraints.append(c1_constraint)

# Pin the absolute start and end points of the entire curve to the original data
p_first_segment = lambda p: p[p_indices[0]:p_indices[1]]
p_last_segment = lambda p: p[p_indices[-2]:p_indices[-1]]
constraints.append({'type': 'eq', 'fun': lambda p: eval_poly(p_first_segment(p), BOUNDARIES[0]) - get_original_vp(BOUNDARIES[0])})
constraints.append({'type': 'eq', 'fun': lambda p: eval_poly(p_last_segment(p), BOUNDARIES[-1]) - get_original_vp(BOUNDARIES[-1])})

# B+ Inequality Constraint for all segments
def check_b_plus(p):
    min_diff = float('inf')
    for i in range(len(DEGREES)):
        p_segment = p[p_indices[i]:p_indices[i+1]]
        grid = np.linspace(BOUNDARIES[i], BOUNDARIES[i+1], 20)
        min_diff = min(min_diff, np.min(CLAMP_HIGH_V - eval_poly(p_segment, grid)))
    return min_diff

constraints.append({'type': 'ineq', 'fun': check_b_plus})


# --- 5. Run the Optimization ---
# Generate a good initial guess (p0) for all five segments
p0 = []
for i in range(len(DEGREES)):
    # Ensure there is data in the segment for the initial fit
    idx_seg = (vg_data >= BOUNDARIES[i]) & (vg_data <= BOUNDARIES[i+1])
    if np.sum(idx_seg) > DEGREES[i]:
        p0.append(np.polyfit(vg_data[idx_seg], vp_data[idx_seg], DEGREES[i]))
    else: # Fallback for very small segments
        p0.append(np.zeros(DEGREES[i] + 1))
p0 = np.concatenate(p0)


print("Starting final 5-segment optimization...")
result = minimize(objective, p0, args=(vg_data, vp_data), method='SLSQP', constraints=constraints, options={'maxiter': 4000, 'ftol': 1e-10})

if result.success:
    print("Optimization successful!")
    p_opt = result.x
else:
    print("Optimization failed.")
    print(result.message)
    exit()

# --- 6. Create the Final Piecewise Model Function ---
def get_model_vp(vg):
    vg = np.asarray(vg)
    vp_out = np.zeros_like(vg, dtype=float)

    vp_out[vg < BOUNDARIES[0]] = CLAMP_HIGH_V
    vp_out[vg > BOUNDARIES[-1]] = CLAMP_LOW_V
    
    for i in range(len(DEGREES)):
        p_segment = p_opt[p_indices[i]:p_indices[i+1]]
        vg_start, vg_end = BOUNDARIES[i], BOUNDARIES[i+1]
        
        cond = (vg >= vg_start) & (vg <= vg_end)
        # On the last segment, make sure the condition is inclusive of the endpoint
        if i == len(DEGREES) - 1:
            cond = (vg >= vg_start) & (vg <= vg_end)
        else:
            cond = (vg >= vg_start) & (vg < vg_end)

        vp_out[cond] = eval_poly(p_segment, vg[cond])
        
    return vp_out

# --- 7. Visualization ---
print("Generating plot...")
vg_plot = np.linspace(-24, 24, 2000)
vp_model_plot = get_model_vp(vg_plot)

plt.style.use('seaborn-v0_8-whitegrid')
fig, ax = plt.subplots(figsize=(16, 9))

ax.plot(vg_data, vp_data, 'o', color='skyblue', label='Original SPICE Data', markersize=4, alpha=0.6, zorder=1)
ax.plot(vg_plot, vp_model_plot, '-', color='red', label='Final 5-Segment Model', linewidth=2.5, zorder=2)

constraint_vgs = BOUNDARIES
constraint_vps = get_model_vp(constraint_vgs)
ax.plot(constraint_vgs, constraint_vps, 'x', color='black', markersize=10, mew=2, label='Segment Boundaries', zorder=3)

ax.axhline(y=CLAMP_HIGH_V, color='gray', linestyle='--', label=f'B+ Voltage Limit ({CLAMP_HIGH_V}V)')

ax.set_title('12AX7 Triode Vp vs. Vg Model (Final 5-Segment High-Fidelity Fit)', fontsize=16)
ax.set_xlabel('Grid Voltage (Vg) [V]', fontsize=12)
ax.set_ylabel('Plate Voltage (Vp) [V]', fontsize=12)
ax.legend(fontsize=11)
ax.grid(True, which='both', linestyle='--', linewidth=0.5)

ax.set_xlim(-25, 25)
ax.set_ylim(-10, CLAMP_HIGH_V + 20)

plt.tight_layout()
plt.show()
