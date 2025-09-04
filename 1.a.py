import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

Q = 10  # Discharge (m³/s)
n = 0.013  # Manning's roughness coefficient
g = 9.81  # Acceleration due to gravity (m/s²)

# Trapezoidal channel
B_trap = 5  # Bottom width (m)
z_trap = 1.5  # Side slope
S0_trap = 0.0008  # Slope
L_trap = 200  # Length of trapezoidal channel (m)

# Upper rectangular channel
B_rect1 = 4  # Width (m)
S0_rect1 = 0.0012  # Slope
L_rect1 = 150  # Length of upper rectangular channel (m)

# Lower rectangular channel
B_rect2 = 4  # Width (m)
S0_rect2 = 0.0024  # Slope
L_rect2 = 250  # Length of lower rectangular channel (m)

# Function to calculate normal depth using Manning's equation
def normal_depth_trapezoidal(B, z, S0, n, Q):
    y = sp.Symbol('y')
    A = B * y + z * y**2
    P = B + 2 * y * sp.sqrt(1 + z**2)
    R = A / P
    eqn = (1/n) * A * R**(2/3) * S0**0.5 - Q
    y_n = sp.nsolve(eqn, y, 1)  # Initial guess = 1m
    return float(y_n)

def normal_depth_rectangular(B, S0, n, Q):
    y = sp.Symbol('y', positive=True, real=True)
    A = B * y
    P = B + 2 * y
    R = A / P
    eqn = (1/n) * A * R**(2/3) * S0**0.5 - Q
    y_n = sp.nsolve(eqn, y, 1)  # Initial guess = 1m
    return float(y_n)

# Function to calculate critical depth using critical flow equation
def critical_depth_rectangular(Q, B):
    y_c = (Q**2 / (g * B**2))**(1/3)
    return float(y_c)

def critical_depth_trapezoidal(Q, B, z):
    y = sp.Symbol('y', positive=True, real=True)
    A = B * y + z * y**2
    T = B + 2 * z * y
    eqn = Q**2 - (A**3 * g / T)
    y_c = sp.nsolve(eqn, y, 1)  # Initial guess = 1m
    return float(y_c)

# Compute depths
y_n_trap = normal_depth_trapezoidal(B_trap, z_trap, S0_trap, n, Q)
y_c_trap = critical_depth_trapezoidal(Q, B_trap, z_trap)

y_n_rect1 = normal_depth_rectangular(B_rect1, S0_rect1, n, Q)
y_c_rect1 = critical_depth_rectangular(Q, B_rect1)

y_n_rect2 = normal_depth_rectangular(B_rect2, S0_rect2, n, Q)
y_c_rect2 = critical_depth_rectangular(Q, B_rect2)

# Print results
print("Normal Depths (m):")
print(f"Trapezoidal Channel: {y_n_trap:.3f}")
print(f"Upper Rectangular Channel: {y_n_rect1:.3f}")
print(f"Lower Rectangular Channel: {y_n_rect2:.3f}")
print("\nCritical Depths (m):")
print(f"Trapezoidal Channel: {y_c_trap:.3f}")
print(f"Upper Rectangular Channel: {y_c_rect1:.3f}")
print(f"Lower Rectangular Channel: {y_c_rect2:.3f}")

# Compute bottom elevations
x_values = [0, L_trap, L_trap + L_rect1, L_trap + L_rect1 + L_rect2]  # Lengths in meters
bottom_elevations = [S0_trap * L_trap + S0_rect1 * L_rect1 + S0_rect2 * L_rect2, S0_rect1 * L_rect1 + S0_rect2 * L_rect2, S0_rect2 * L_rect2, 0]

# Normal and critical depths above bottom profile
y_n_values = [y_n_trap + bottom_elevations[0], y_n_trap + bottom_elevations[1], y_n_rect1 + bottom_elevations[2], y_n_rect2 + bottom_elevations[3]]
y_c_values = [y_c_trap + bottom_elevations[0], y_c_trap + bottom_elevations[1], y_c_rect1 + bottom_elevations[2], y_c_rect2 + bottom_elevations[3]]

# Plotting Water Surface Profile
plt.figure(figsize=(10, 5))
plt.plot(x_values, bottom_elevations, label="Channel Bottom", linestyle="-", color="black")
plt.plot(x_values, y_n_values, label="Normal Depth", marker="o", linestyle="-", color="maroon")
plt.plot(x_values, y_c_values, label="Critical Depth", marker="s", linestyle="--", color="orange")
plt.plot(x_values, [2, 1.8, 1.6, 0.86], label="Water Surface Elevation", marker="v", color="blue")
# Labels and titles
plt.xlabel("Channel Length (m)")
plt.ylabel("Elevation / Depth (m)")
plt.title("Water Surface Profile with Channel Bottom")
plt.legend()
plt.grid()
# Annotate sections
plt.text(60, y_n_trap, "Trapezoidal Channel")
plt.text(200, y_n_rect1 - 0.3, "Upper Rectangular Channel")
plt.text(450, y_n_rect2 - 0.5, "Lower Rectangular Channel")
plt.show()
