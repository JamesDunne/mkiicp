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

# *** FINAL MODEL: 5-SEGMENTS WITH OPTIMIZED DEGREES ***
# Updated last segment to quadratic as requested.
DEGREES = [2, 3, 2, 2, 2] 
COEFF_COUNTS = [d + 1 for d in DEGREES]
TOTAL_COEFFS = sum(COEFF_COUNTS)
NUM_INTERNAL_BOUNDARIES = 4

# Define FIXED outer boundaries
VG_START = -8.0
VG_END = 17.62846

# Clamping values
CLAMP_HIGH_V = 405.0
CLAMP_LOW_V = get_original_vp(VG_END)

# --- 3. Setup the Combined Optimization Problem ---
def objective_with_boundaries(x, vg_data, vp_data):
    internal_boundaries = x[:NUM_INTERNAL_BOUNDARIES]
    coeffs = x[NUM_INTERNAL_BOUNDARIES:]
    
    all_boundaries = [VG_START] + list(internal_boundaries) + [VG_END]
    p_indices = np.cumsum([0] + COEFF_COUNTS)

    total_error = 0
    for i in range(len(DEGREES)):
        p_segment = coeffs[p_indices[i]:p_indices[i+1]]
        vg_start, vg_end = all_boundaries[i], all_boundaries[i+1]
        
        idx_segment = (vg_data >= vg_start) & (vg_data <= vg_end)
        if np.any(idx_segment):
            error = np.sum((vp_data[idx_segment] - np.polyval(p_segment, vg_data[idx_segment]))**2)
            total_error += error
            
    return total_error

# --- 4. Define Constraints for the Combined Problem ---
constraints = []
p_indices = np.cumsum([0] + COEFF_COUNTS)

# Constraint 1: Boundary ordering (with a small margin to prevent collapse)
constraints.append({'type': 'ineq', 'fun': lambda x: x[0] - VG_START - 0.2})
constraints.append({'type': 'ineq', 'fun': lambda x: x[1] - x[0] - 0.2})
constraints.append({'type': 'ineq', 'fun': lambda x: x[2] - x[1] - 0.2})
constraints.append({'type': 'ineq', 'fun': lambda x: x[3] - x[2] - 0.2})
constraints.append({'type': 'ineq', 'fun': lambda x: VG_END - x[3] - 0.2})

# Constraint 2: C0 and C1 Continuity
for i in range(NUM_INTERNAL_BOUNDARIES):
    get_p_before = (lambda i_cap: lambda x: x[NUM_INTERNAL_BOUNDARIES + p_indices[i_cap] : NUM_INTERNAL_BOUNDARIES + p_indices[i_cap+1]])(i)
    get_p_after = (lambda i_cap: lambda x: x[NUM_INTERNAL_BOUNDARIES + p_indices[i_cap+1] : NUM_INTERNAL_BOUNDARIES + p_indices[i_cap+2]])(i)
    
    constraints.append({'type': 'eq', 'fun': lambda x, i=i, pb=get_p_before, pa=get_p_after: np.polyval(pb(x), x[i]) - np.polyval(pa(x), x[i])})
    constraints.append({'type': 'eq', 'fun': lambda x, i=i, pb=get_p_before, pa=get_p_after: np.polyval(np.polyder(pb(x)), x[i]) - np.polyval(np.polyder(pa(x)), x[i])})

# Constraint 3: Pin the absolute start and end points
p_first = lambda x: x[NUM_INTERNAL_BOUNDARIES + p_indices[0] : NUM_INTERNAL_BOUNDARIES + p_indices[1]]
p_last = lambda x: x[NUM_INTERNAL_BOUNDARIES + p_indices[-2] : NUM_INTERNAL_BOUNDARIES + p_indices[-1]]
constraints.append({'type': 'eq', 'fun': lambda x: np.polyval(p_first(x), VG_START) - get_original_vp(VG_START)})
constraints.append({'type': 'eq', 'fun': lambda x: np.polyval(p_last(x), VG_END) - get_original_vp(VG_END)})

# Constraint 4: B+ Inequality Constraint
def check_b_plus(x):
    internal_boundaries = x[:NUM_INTERNAL_BOUNDARIES]
    coeffs = x[NUM_INTERNAL_BOUNDARIES:]
    all_boundaries = [VG_START] + list(internal_boundaries) + [VG_END]

    min_diff = float('inf')
    for i in range(len(DEGREES)):
        p_segment = coeffs[p_indices[i]:p_indices[i+1]]
        grid = np.linspace(all_boundaries[i], all_boundaries[i+1], 20)
        min_diff = min(min_diff, np.min(CLAMP_HIGH_V - np.polyval(p_segment, grid)))
    return min_diff

constraints.append({'type': 'ineq', 'fun': check_b_plus})

# --- 5. Run the Final Optimization ---
initial_boundaries = [-5.3, -1.0, 1.0, 6.0]
initial_all_boundaries = [VG_START] + initial_boundaries + [VG_END]
initial_coeffs = []
for i in range(len(DEGREES)):
    idx_seg = (vg_data >= initial_all_boundaries[i]) & (vg_data <= initial_all_boundaries[i+1])
    if np.sum(idx_seg) > DEGREES[i]:
        initial_coeffs.append(np.polyfit(vg_data[idx_seg], vp_data[idx_seg], DEGREES[i]))
    else:
        initial_coeffs.append(np.zeros(DEGREES[i] + 1))
x0 = np.concatenate([initial_boundaries, np.concatenate(initial_coeffs)])

print("Starting final optimization of boundaries and coefficients...")
result = minimize(objective_with_boundaries, x0, args=(vg_data, vp_data), method='SLSQP', constraints=constraints, options={'maxiter': 6000, 'ftol': 1e-11})

if result.success:
    print("Optimization successful!")
    x_opt = result.x
else:
    print("Optimization failed or did not fully converge, but the result is likely an improvement.")
    print(result.message)
    x_opt = result.x

# --- 6. Unpack and Display Final Optimized Model ---
optimized_boundaries = x_opt[:NUM_INTERNAL_BOUNDARIES]
optimized_coeffs = x_opt[NUM_INTERNAL_BOUNDARIES:]
final_all_boundaries = [VG_START] + list(optimized_boundaries) + [VG_END]

print("\n\n" + "="*80)
print("--- Final Optimized Model Parameters for C++ Implementation ---")
print("="*80)
print(f"\nUpper Clamp (for Vg < {final_all_boundaries[0]:.6f}): {CLAMP_HIGH_V:.6f}")
print(f"Lower Clamp (for Vg >= {final_all_boundaries[-1]:.6f}): {CLAMP_LOW_V:.6f}\n")

degree_names = {1: "Linear (ax + b)", 2: "Quadratic (ax^2 + bx + c)", 3: "Cubic (ax^3 + bx^2 + cx + d)"}

for i in range(len(DEGREES)):
    p_segment = optimized_coeffs[p_indices[i]:p_indices[i+1]]
    vg_start, vg_end = final_all_boundaries[i], final_all_boundaries[i+1]
    
    # Determine segment range string
    if i < len(DEGREES) - 1:
        range_str = f"Vg in [{vg_start:.6f}, {vg_end:.6f})"
    else: # Last segment is inclusive of the end boundary
        range_str = f"Vg in [{vg_start:.6f}, {vg_end:.6f}]"

    print(f"--- Segment {i+1}: {range_str} ---")
    print(f"  Type: {degree_names[DEGREES[i]]}")
    # Note on coefficient order
    print("  Coefficients (in descending power order, e.g., a, b, c, ...):")
    # Format each coefficient in scientific notation for precision
    formatted_coeffs = [f"{c:.8e}" for c in p_segment]
    print(f"  {formatted_coeffs}\n")

print("="*80)


# --- 7. Visualization ---
def get_model_vp(vg):
    vg = np.asarray(vg)
    vp_out = np.zeros_like(vg, dtype=float)
    vp_out[vg < final_all_boundaries[0]] = CLAMP_HIGH_V
    vp_out[vg > final_all_boundaries[-1]] = CLAMP_LOW_V
    
    for i in range(len(DEGREES)):
        p_segment = optimized_coeffs[p_indices[i]:p_indices[i+1]]
        vg_start, vg_end = final_all_boundaries[i], final_all_boundaries[i+1]
        cond = (vg >= vg_start) & (vg <= vg_end if i == len(DEGREES) - 1 else vg < vg_end)
        vp_out[cond] = np.polyval(p_segment, vg[cond])
        
    return vp_out

print("\nGenerating final verification plot...")
vg_plot = np.linspace(-24, 24, 2000)
vp_model_plot = get_model_vp(vg_plot)

plt.style.use('seaborn-v0_8-whitegrid')
fig, ax = plt.subplots(figsize=(16, 9))
ax.plot(vg_data, vp_data, 'o', color='skyblue', label='Original SPICE Data', markersize=4, alpha=0.6, zorder=1)
ax.plot(vg_plot, vp_model_plot, '-', color='red', label='Final Optimized Model', linewidth=2.5, zorder=2)
ax.plot(final_all_boundaries, get_model_vp(final_all_boundaries), 'x', color='black', markersize=10, mew=2, label='Optimized Segment Boundaries', zorder=3)
ax.axhline(y=CLAMP_HIGH_V, color='gray', linestyle='--', label=f'B+ Voltage Limit ({CLAMP_HIGH_V}V)')
ax.set_title('12AX7 Triode Vp vs. Vg Model (Final Boundary-Optimized Fit)', fontsize=16)
ax.set_xlabel('Grid Voltage (Vg) [V]', fontsize=12)
ax.set_ylabel('Plate Voltage (Vp) [V]', fontsize=12)
ax.legend(fontsize=11)
ax.grid(True, which='both', linestyle='--', linewidth=0.5)
ax.set_xlim(-25, 25)
ax.set_ylim(-10, CLAMP_HIGH_V + 20)
plt.tight_layout()
plt.show()
