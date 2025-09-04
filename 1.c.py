import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

# Given Data
Q = 10  # Discharge (m³/s)
n = 0.013  # Manning's roughness coefficient
g = 9.81  # Gravity (m/s²)

# Trapezoidal Channel Data
B_trap = 5  # Bottom width (m)
z_trap = 1.5  # Side slope (H:V)
S0_trap = 0.0008  # Bed slope
L_trap = 200  # Length (m)

# Rectangular Channel Data
B_rect = 4  # Width (m)
S0_rect1 = 0.0012  # Slope (upper reach)
L_rect1 = 150  # Length (m)
S0_rect2 = 0.0024  # Slope (lower reach)
L_rect2 = 250  # Length (m)

# Function to Compute Normal Depth using Manning's Equation
def normal_depth(B, z, S0, n, Q, channel_type):
    y = sp.Symbol('y')
    if channel_type == "trapezoidal":
        A = B * y + z * y**2
        P = B + 2 * y * (1 + z**2)**0.5
    else:  # Rectangular
        A = B * y
        P = B + 2 * y
    R = A / P
    eqn = (1/n) * A * R**(2/3) * S0**0.5 - Q
    y_n = sp.nsolve(eqn, y, 1)
    return float(y_n)

# Function to Compute Critical Depth using Energy Equation
def critical_depth(B, z, Q, channel_type):
    y = sp.Symbol('y', positive=True, real=True)
    if channel_type == "trapezoidal":
        A = B * y + z * y**2
        T = B + 2 * z * y
    else:  # Rectangular
        A = B * y
        T = B
    eqn = Q**2 - (A**3 * g / T)
    y_c = sp.nsolve(eqn, y, 1)
    return float(y_c)

# Compute depths
y_n_trap = normal_depth(B_trap, z_trap, S0_trap, n, Q, "trapezoidal")
y_c_trap = critical_depth(B_trap, z_trap, Q, "trapezoidal")

y_n_rect1 = normal_depth(B_rect, 0, S0_rect1, n, Q, "rectangular")
y_c_rect1 = critical_depth(B_rect, 0, Q, "rectangular")

y_n_rect2 = normal_depth(B_rect, 0, S0_rect2, n, Q, "rectangular")
y_c_rect2 = critical_depth(B_rect, 0, Q, "rectangular")

# Print results
print("Normal Depths (m):")
print(f"Trapezoidal Channel: {y_n_trap:.3f}")
print(f"Upper Rectangular Channel: {y_n_rect1:.3f}")
print(f"Lower Rectangular Channel: {y_n_rect2:.3f}")

print("\nCritical Depths (m):")
print(f"Trapezoidal Channel: {y_c_trap:.3f}")
print(f"Upper Rectangular Channel: {y_c_rect1:.3f}")
print(f"Lower Rectangular Channel: {y_c_rect2:.3f}")

channels = [
    {"type": "trapezoidal", "B": 5, "z": 1.5, "S0": 0.0008, "L": 200},
    {"type": "rectangular", "B": 4, "z": 0, "S0": 0.0012, "L": 150},
    {"type": "rectangular", "B": 4, "z": 0, "S0": 0.0024, "L": 250}
]

dx = 10  # Step size (m)

# Compute Normal and Critical Depths
y_n_values = []
y_c_values = []
for ch in channels:
    y_n_values.append(normal_depth(ch["B"], ch["z"], ch["S0"], n, Q, ch["type"]))
    y_c_values.append(critical_depth(ch["B"], ch["z"], Q, ch["type"]))

def flow_properties(y, B, z, channel_type):
    if channel_type == "trapezoidal":
        A = B * y + z * y**2
        P = B + 2 * y * (1 + z**2)**0.5
    else:  # Rectangular
        A = B * y
        P = B + 2 * y
    R = A / P  # Hydraulic radius
    return A, R

# Iterative Standard-Step Method Function
def compute_water_profile(channel, y_start, direction="upstream"):
    B, z, S0, L, channel_type = channel["B"], channel["z"], channel["S0"], channel["L"], channel["type"]
    steps = int(L / dx) + 1
    x = np.linspace(0, L, steps)
    y = np.zeros(steps)
    y[0] = y_start  # Set initial depth

    for i in range(1, steps):
        y_prev = y[i-1]
        A_prev, R_prev = flow_properties(y_prev, B, z, channel_type)
        V_prev = Q / A_prev
        EGL_prev = y_prev + (V_prev**2 / (2 * g))  # Energy Grade Line
            
        # Solve for y_next using the energy equation
        y_next = sp.Symbol('y')
        A_next, R_next = flow_properties(y_next, B, z, channel_type)
        V_next = Q / A_next
        EGL_next = y_next + (V_next**2 / (2 * g))
        eqn = EGL_prev - EGL_next
            
        try:
            y[i] = float(sp.nsolve(eqn, y_next, y_prev))  # Solve for new depth
        except:
            y[i] = y_prev  # Keeping previous depth if no solution found
        
    if direction == "downstream":
        x = x[::-1]  # Reverse if computing downstream

    return x, y

x_total = []
y_total = []
offset = 0
cumulative_bottom = 0
bottom_total = []
y_start = y_c_values[-1]  # Assume flow starts at critical depth (free overfall)
for i, ch in enumerate(reversed(channels)):  # Compute upstream
    x, y = compute_water_profile(ch, y_start, direction="upstream")
    x_total.extend(x + offset)
    y_total.extend(y)
    offset += ch["L"]
    y_start = y[-1]
def energy_grade_line(y, Q, z):
    if z > 0:
        A = 5 * y + (1.5 * y**2)
    else:  # Rectangular
        A = 4 * y   
    V = Q / A
    EGL = y + (V**2 / (2 * g))
    return EGL

# Define Sections
x_values = np.array([0, L_trap, L_trap + L_rect1, L_trap + L_rect1 + L_rect2])
y_n_values = np.array([y_n_trap, y_n_trap, y_n_rect1, y_n_rect2])
y_c_values = np.array([y_c_trap, y_c_trap, y_c_rect1, y_c_rect2])
channel_types = ['trapezoidal' if x <= L_trap else 'rectangular' for x in x_values]
EGL_values = []
current_offset = 0
for ch in reversed(channels):
    x_ch = np.linspace(0, ch["L"], int(ch["L"] / dx) + 1)
    y_ch = y_total[current_offset:current_offset + len(x_ch)]
    EGL_values.extend([energy_grade_line(y, Q, ch["z"]) for y in y_ch])
    current_offset += len(x_ch)

# Channel Bottom Elevation Profile (Corrected)
bottom_elevations = []
current_elevation = S0_trap * L_trap + S0_rect1 * L_rect1 + S0_rect2 * L_rect2 # starting elevation
offset = 0
for ch in reversed(channels):
    bottom_elevations.extend(np.linspace(current_elevation, current_elevation - ch["S0"] * ch["L"], int(ch["L"] / dx) + 1))
    current_elevation -= ch["S0"] * ch["L"]
    offset += ch["L"]

bottom_elevations = np.array(bottom_elevations)
x_total = np.array(x_total)

# Water Surface Elevation (Corrected)
water_surface = np.array(y_total) + bottom_elevations

# Plot Water Surface Profile
plt.figure(figsize=(12, 6))
plt.plot(x_total, bottom_elevations, label="Channel Bottom", color="black", linestyle="-")
plt.plot(x_total, water_surface, label="Water Surface Elevation", color="blue")

# Calculate x_values for normal and critical depth lines
x_values = [0, L_trap, L_trap + L_rect1, L_trap + L_rect1 + L_rect2]
bottom_elevation_at_x = [bottom_elevations[0], bottom_elevations[int(L_trap/dx)], bottom_elevations[int((L_trap + L_rect1)/dx)], bottom_elevations[-1]]

plt.plot(x_values, [bottom_elevation_at_x[i] + y_c_values[i] for i in range(4)], label="Critical Depth Line", color="red")
plt.plot(x_values, [bottom_elevation_at_x[i] + y_n_values[i] for i in range(4)], label="Normal Depth Line", color="green")
plt.plot(x_total, [y + 0.6 for y in EGL_values], label="Energy Grade Line (EGL)", color="purple")

# Labels and Titles
plt.xlabel("Channel Length (m)")
plt.ylabel("Elevation (m)")
plt.title("Water Surface Profile using Standard-Step Method")
plt.legend()
plt.grid()

# Annotate sections
plt.text(50, water_surface[0], "Trapezoidal Channel")
plt.text(275, water_surface[1], "Upper Rectangular Channel")
plt.text(450, water_surface[2], "Lower Rectangular Channel")
plt.show()