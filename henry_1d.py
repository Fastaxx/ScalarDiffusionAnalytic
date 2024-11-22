import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf, erfc

def analytical_solution(x, t, k1, k2, He, x0, C0):
    """
    Calculate analytical solution for the 1D diffusion problem
    
    Parameters:
    -----------
    x : array-like
        Spatial coordinates
    t : float
        Time
    k1, k2 : float
        Diffusion coefficients
    He : float
        Henry's constant
    x0 : float
        Interface position
    """
    T1 = np.zeros_like(x)
    T2 = np.zeros_like(x)
    
    # Domain 1 (x < x0)
    mask1 = x > x0
    T1[mask1] = C0/(1 + He*np.sqrt(k2/k1)) * (1 + He * np.sqrt(k2/k1) * erf(x[mask1]/(2 * np.sqrt(k1 * t))))
    
    # Domain 2 (x > x0)
    mask2 = x <= x0
    T2[mask2] = (He*C0)/(1+He*np.sqrt(k2/k1)) * erfc(np.abs(x[mask2])/(2 * np.sqrt(k2 * t)))
    
    return T1, T2

def plot_solutions(x, t, k1, k2, He, x0, C0):
    """Plot both analytical and numerical solutions"""
    # Calculate solutions
    T1_analytical, T2_analytical = analytical_solution(x, t, k1, k2, He, x0, C0)

    # Create plot
    plt.figure(figsize=(12, 6))
    
    # Plot analytical solution
    plt.plot(x[x > x0], T1_analytical[x > x0], 'b-', label='T1 Analytical')
    plt.plot(x[x <= x0], T2_analytical[x <= x0], 'r-', label='T2 Analytical')
    
    plt.axvline(x=x0, color='k', linestyle=':', label='Interface')
    plt.grid(True)
    plt.xlabel('Position (x)')
    plt.ylabel('Temperature (T)')
    plt.title(f'Temperature Distribution at t = {t}')
    plt.legend()
    plt.show()

# Example usage
if __name__ == "__main__":
    # Parameters
    k1 = 1.0  # Diffusion coefficient in domain 1
    k2 = 4.0  # Diffusion coefficient in domain 2
    He = 2.0  # Henry's constant
    C0 = 1.0
    x0 = 0.0  # Interface position
    
    # Domain
    x = np.linspace(-5, 5, 2000)
    t = 0.1
    
    # Plot solutions
    plot_solutions(x, t, k1, k2, He, x0, C0)