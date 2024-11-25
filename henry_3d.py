import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import warnings

# Suppress integration warnings

# Define parameters
Dg = 1.0    # Diffusivity inside the bubble (gas phase) [m^2/s]
Dl = 4.0    # Diffusivity outside the bubble (liquid phase) [m^2/s]
R0 = 1.0    # Bubble radius [m]
d0 = 2 * R0  # Bubble diameter [m]
alpha = 1.0  # Dimensionless Henry's law constant
cg0 = 1.0    # Initial concentration inside the bubble [mol/m³]

# Maximum integration limit (since infinity is not practical)
x_max = 1000  # Adjust as needed for convergence

# Time array for plotting (from t=0.01 to t=5 seconds)
t_values = np.linspace(0.01, 5.0, num=200)

# Positions
r_inside = 0.5 * R0  # r/d0 = 0.25 -> r = 0.5 * R0
r_outside = 1.5 * R0 # r/d0 = 0.75 -> r = 1.5 * R0

# Define lambda_g(-x) and lambda_l(-x)
def lambda_g(x):
    return 1j * np.sqrt(x / Dg)

def lambda_l(x):
    return 1j * np.sqrt(x / Dl)

# Define xi(-x) and zeta(-x)
def xi(x):
    lgR0 = lambda_g(x) * R0
    llR0 = lambda_l(x) * R0
    numerator = 2 * Dg * (lgR0 * np.cosh(lgR0) - np.sinh(lgR0))
    denominator = Dl * np.exp(-llR0) * (1 + llR0)
    return numerator / denominator

def zeta(x):
    llR0 = lambda_l(x) * R0
    lgR0 = lambda_g(x) * R0
    term1 = xi(x) * np.exp(-llR0) / (alpha * R0)
    term2 = 2 / R0 * np.sinh(lgR0)
    return term1 + term2

# Define the integrand for cg(t, r)
def integrand_cg(x, r, t):
    lg_r = lambda_g(x) * r
    numerator = np.sinh(lg_r)
    denominator = zeta(x)
    integrand = np.imag((numerator / denominator) * np.exp(-x * t)) / x
    return integrand

# Compute cg(t, r)
def cg(t, r):
    integral, _ = quad(integrand_cg, 0, x_max, args=(r, t), limit=500)
    return -2 * cg0 / (np.pi * r) * integral

# Define the integrand for cl(t, r)
def integrand_cl(x, r, t):
    ll_r = lambda_l(x) * r
    numerator = xi(x) * np.exp(-ll_r)
    denominator = zeta(x)
    integrand = np.imag((numerator / denominator) * np.exp(-x * t)) / x
    return integrand

# Compute cl(t, r)
def cl(t, r):
    integral, _ = quad(integrand_cl, 0, x_max, args=(r, t), limit=500)
    return cg0 / (np.pi * r) * integral

# (a) Compute cg(t) at r/d0 = 0.25 (r = 0.5 R0)
cg_values = []
for t in t_values:
    c_g_value = cg(t, r_inside)
    cg_values.append(c_g_value)

# (b) Compute cl(t) at r/d0 = 0.75 (r = 1.5 R0)
cl_values = []
for t in t_values:
    c_l_value = cl(t, r_outside)
    cl_values.append(c_l_value)

# (c) Radial concentration profile at time t = 1.2 s
t_fixed = 1.2  # Fixed time
# Radial positions inside and outside the bubble
r_values_inside = np.linspace(1e-3, R0, num=100)  # Avoid r=0
r_values_outside = np.linspace(R0, 5 * R0, num=200)

# Compute concentration profiles inside the bubble
cg_profile = []
for r in r_values_inside:
    c_g_value = cg(t_fixed, r)
    cg_profile.append(c_g_value)

# Compute concentration profiles outside the bubble
cl_profile = []
for r in r_values_outside:
    c_l_value = cl(t_fixed, r)
    cl_profile.append(c_l_value)

# Combine radial positions and concentrations for plotting
r_values = np.concatenate((r_values_inside, r_values_outside))
c_values = np.concatenate((cg_profile, cl_profile))

# Plot all three plots in a single figure with subplots
fig, axs = plt.subplots(2, 2, figsize=(12, 10))

# Plot (a) Concentration inside the bubble at r/d0 = 0.25
axs[0, 0].plot(t_values, cg_values, label=f'c_g(t) at r/d₀ = 0.25', color='blue')
axs[0, 0].set_xlabel('Time t (s)')
axs[0, 0].set_ylabel('Concentration c_g (mol/m³)')
axs[0, 0].set_title('Concentration inside the Bubble at r/d₀ = 0.25')
axs[0, 0].legend()
axs[0, 0].grid(True)

# Plot (b) Concentration outside the bubble at r/d0 = 0.75
axs[0, 1].plot(t_values, cl_values, label=f'c_l(t) at r/d₀ = 0.75', color='green')
axs[0, 1].set_xlabel('Time t (s)')
axs[0, 1].set_ylabel('Concentration c_l (mol/m³)')
axs[0, 1].set_title('Concentration outside the Bubble at r/d₀ = 0.75')
axs[0, 1].legend()
axs[0, 1].grid(True)

# Plot (c) Radial concentration profile at t = 1.2 s
axs[1, 0].plot(r_values / d0, c_values, label=f'Concentration profile at t = {t_fixed} s', color='red')
axs[1, 0].axvline(x=R0 / d0, color='k', linestyle='--', label='Bubble boundary (r/d₀ = 0.5)')
axs[1, 0].set_xlabel('Radial position r/d₀')
axs[1, 0].set_ylabel('Concentration c (mol/m³)')
axs[1, 0].set_title('Radial Concentration Profile at t = 1.2 s')
axs[1, 0].legend()
axs[1, 0].grid(True)

# Remove the empty subplot (bottom right)
fig.delaxes(axs[1, 1])

# Adjust layout
plt.tight_layout()

# Display the combined figure
plt.show()
