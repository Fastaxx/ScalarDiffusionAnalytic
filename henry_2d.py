import numpy as np
from scipy.special import j0, j1, y0, y1
from scipy.integrate import quad
import matplotlib.pyplot as plt

# Physical parameters
Dg = 1.0    # Gas diffusivity
Dl = 4.0    # Liquid diffusivity
R0 = 2.0    # Bubble radius
cg0 = 1.0   # Initial concentration inside bubble
He = 2.0    # Henry

D = np.sqrt(Dg / Dl)

# Define Φ(u) and Ψ(u)
def Phi(u):
    term1 = Dg * np.sqrt(Dl) * j1(u * R0) * y0(D * u * R0)
    term2 = He * Dl * np.sqrt(Dg) * j0(u * R0) * y1(D * u * R0)
    return term1 - term2

def Psi(u):
    term1 = Dg * np.sqrt(Dl) * j1(u * R0) * j0(D * u * R0)
    term2 = He * Dl * np.sqrt(Dg) * j0(u * R0) * j1(D * u * R0)
    return term1 - term2

# Integrand for c_g(r, t)
def cg_integrand(u, r, t):
    Phi_u = Phi(u)
    Psi_u = Psi(u)
    denominator = u**2 * (Phi_u**2 + Psi_u**2)
    numerator = np.exp(-Dg * u**2 * t) * j0(u * r) * j1(u * R0)
    # Handle division by zero
    if denominator == 0:
        return 0.0
    else:
        return numerator / denominator

# Integrand for c_l(r, t)
def cl_integrand(u, r, t):
    Phi_u = Phi(u)
    Psi_u = Psi(u)
    denominator = u * (Phi_u**2 + Psi_u**2)
    term1 = j0(D * u * r) * Phi_u
    term2 = y0(D * u * r) * Psi_u
    numerator = np.exp(-Dg * u**2 * t) * j1(u * R0) * (term1 - term2)
    # Handle division by zero
    if denominator == 0:
        return 0.0
    else:
        return numerator / denominator

# Compute c_g(r, t)
def compute_cg(r_values, t_values):
    prefactor = (4 * cg0 * Dg * Dl * Dl * He) / (np.pi**2 * R0)
    cg_results = np.zeros((len(t_values), len(r_values)))
    for i, t in enumerate(t_values):
        Umax = 5.0 / np.sqrt(Dg * t)
        for j, r in enumerate(r_values):
            integrand = lambda u: cg_integrand(u, r, t)
            integral, error = quad(integrand, 0, Umax, limit=500, epsabs=1e-6, epsrel=1e-6)
            cg_results[i, j] = prefactor * integral
    return cg_results

# Compute c_l(r, t)
def compute_cl(r_values, t_values):
    prefactor = (2 * cg0 * Dg * np.sqrt(Dl)*He) / np.pi
    cl_results = np.zeros((len(t_values), len(r_values)))
    for i, t in enumerate(t_values):
        Umax = 5.0 / np.sqrt(Dg * t)
        for j, r in enumerate(r_values):
            integrand = lambda u: cl_integrand(u, r, t)
            integral, error = quad(integrand, 0, Umax, limit=500, epsabs=1e-6, epsrel=1e-6)
            cl_results[i, j] = prefactor * integral
    return cl_results

# Main script
if __name__ == "__main__":
    # Spatial and temporal points
    r_values_inside = np.linspace(1e-6, R0, 100)  # Inside the bubble (avoid r=0)
    r_values_outside = np.linspace(R0, 4 * R0, 100)  # Outside the bubble
    t_values = [0.01,0.25,0.5,0.75,1.0]

    # Compute c_g(r, t)
    print("Computing c_g(r, t)...")
    cg_results = compute_cg(r_values_inside, t_values)

    # Compute c_l(r, t)
    print("Computing c_l(r, t)...")
    cl_results = compute_cl(r_values_outside, t_values)

    # Plotting c_g(r, t)
    plt.figure(figsize=(12, 5))

    for i, t in enumerate(t_values):
        plt.plot(r_values_inside / R0, cg_results[i], label=f'c₉ (t={t})', linestyle='-')
        plt.plot(r_values_outside / R0, cl_results[i], label=f'cₗ (t={t})', linestyle='--')

    plt.xlabel('r / R₀')
    plt.ylabel('Concentration')
    plt.title('Concentration Profiles (Gas and Liquid)')
    plt.legend()
    plt.grid(True)
    plt.show()
