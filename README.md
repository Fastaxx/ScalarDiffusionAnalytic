# ScalarDiffusionAnalytic

Python scripts to compute the analytical solutions of the 2-phase diffusion problem in 1D, 2D, and 3D.

---

## Table of Contents

- [Introduction](#introduction)
- [Features](#features)
- [Theoretical Background](#theoretical-background)
  - [1D Diffusion](#1d-diffusion)
  - [2D Diffusion](#2d-diffusion)
  - [3D Diffusion](#3d-diffusion)
- [Installation](#installation)
- [Usage](#usage)
  - [Running the Scripts](#running-the-scripts)
  - [Example Outputs](#example-outputs)
- [Dependencies](#dependencies)
- [Contributing](#contributing)
- [License](#license)
- [Acknowledgments](#acknowledgments)

---

## Introduction

**ScalarDiffusionAnalytic** is a collection of Python scripts that compute analytical solutions for two-phase diffusion problems in one-dimensional (1D), two-dimensional (2D), and three-dimensional (3D) geometries. These scripts are designed to visualize the diffusion processes between two different media with varying diffusivities and interface conditions.

---

## Features

- **Analytical Solutions**: Provides exact solutions derived using methods such as Laplace transforms, separation of variables, and integral transforms.
- **Dimensional Flexibility**: Includes scripts for 1D, 2D (cylindrical coordinates), and 3D (spherical coordinates) diffusion problems.
- **Customizable Parameters**: Allows users to input physical parameters such as diffusivities, interface conditions, initial concentrations, and geometrical dimensions.
- **Visualization**: Generates plots of concentration profiles over time and space to aid in the visualization of diffusion processes.
- **Henry's Law Interface Conditions**: Incorporates different interface conditions, including Henry's Law with adjustable constants.

---

## Theoretical Background

### 1D Diffusion

**Problem Statement**:

We consider a one-dimensional diffusion problem where two semi-infinite media meet at an interface located at \( x = x_0 \). The diffusivities in the two media are $$ D_1 $$ and $$ D_2 $$, respectively.

**Governing Equations**:

- For x < x0 :
  $$
  \frac{\partial c_1}{\partial t} = D_1 \frac{\partial^2 c_1}{\partial x^2}
  $$
- For \( x > x_0 \):
  $$
  \frac{\partial c_2}{\partial t} = D_2 \frac{\partial^2 c_2}{\partial x^2}
  $$

**Interface Conditions**:

- **Concentration Discontinuity (Henry's Law)**:
  $$
  c_1(x_0, t) = \alpha c_2(x_0, t)
  $$
- **Flux Continuity**:
  $$
  D_1 \left. \frac{\partial c_1}{\partial x} \right|_{x = x_0^-} = D_2 \left. \frac{\partial c_2}{\partial x} \right|_{x = x_0^+}
  $$

**Solution Approach**:

- Apply Laplace transforms to eliminate time dependence.
- Solve the resulting ordinary differential equations (ODEs).
- Apply boundary and interface conditions to determine constants.
- Perform inverse Laplace transform to obtain the time-dependent solution.

### 2D Diffusion

**Problem Statement**:

We analyze diffusion from a circular bubble (in 2D) of radius $$ R_0 $$ into the surrounding medium.

**Governing Equations**:

- Inside the bubble (gas phase):
  $$
  \frac{\partial c_g}{\partial t} = D_g \left( \frac{\partial^2 c_g}{\partial r^2} + \frac{1}{r} \frac{\partial c_g}{\partial r} \right)
  $$
- Outside the bubble (liquid phase):
  $$
  \frac{\partial c_l}{\partial t} = D_l \left( \frac{\partial^2 c_l}{\partial r^2} + \frac{1}{r} \frac{\partial c_l}{\partial r} \right)
  $$

**Interface Conditions**:

- **Henry's Law**:
  $$
  c_g(R_0, t) = \alpha c_l(R_0, t)
  $$
- **Flux Continuity**:
  $$
  D_g \left. \frac{\partial c_g}{\partial r} \right|_{r = R_0^-} = D_l \left. \frac{\partial c_l}{\partial r} \right|_{r = R_0^+}
  $$

**Solution Approach**:

- Use Laplace transforms and Bessel functions to solve the radial diffusion equations.
- Apply boundary conditions to find constants.
- Invert transforms to obtain solutions in the time domain.

### 3D Diffusion

**Problem Statement**:

We consider diffusion from a spherical bubble (in 3D) of radius $$ R_0 $$ into the surrounding liquid.

**Governing Equations**:

- Inside the bubble:
  $$
  \frac{\partial c_g}{\partial t} = D_g \left( \frac{\partial^2 c_g}{\partial r^2} + \frac{2}{r} \frac{\partial c_g}{\partial r} \right)
  $$
- Outside the bubble:
  $$
  \frac{\partial c_l}{\partial t} = D_l \left( \frac{\partial^2 c_l}{\partial r^2} + \frac{2}{r} \frac{\partial c_l}{\partial r} \right)
  $$

**Interface Conditions**:

- **Henry's Law**:
  $$
  c_g(R_0, t) = \alpha c_l(R_0, t)
  $$
- **Flux Continuity**:
  $$
  D_g \left. \frac{\partial c_g}{\partial r} \right|_{r = R_0^-} = D_l \left. \frac{\partial c_l}{\partial r} \right|_{r = R_0^+}
  $$

**Solution Approach**:

- Similar to the 2D case but involves spherical Bessel functions.
- The solutions are derived using separation of variables and series expansions.

---

## Installation

1. **Clone the Repository**:

   ```bash
   git clone https://github.com/yourusername/ScalarDiffusionAnalytic.git
   ```

2. **Navigate to the Directory**:

   ```bash
   cd ScalarDiffusionAnalytic
   ```

3. **Install Dependencies**:

   Ensure you have Python 3 installed along with the necessary libraries.

   ```bash
   pip install numpy scipy matplotlib
   ```

---

## Usage

### Running the Scripts

Each dimension has its own script:

- **1D Diffusion**: `diffusion_1d.py`
- **2D Diffusion**: `diffusion_2d.py`
- **3D Diffusion**: `diffusion_3d.py`

To run a script, use the following command:

```bash
python diffusion_1d.py
```

### Configuring Parameters

At the beginning of each script, you can set the physical and numerical parameters:

```python
# Physical parameters
D1 = 1.0     # Diffusivity in medium 1
D2 = 0.5     # Diffusivity in medium 2
alpha = 0.8  # Henry's law constant
c0 = 1.0     # Initial concentration

# Simulation parameters
x0 = 0.0     # Interface position
t_values = [0.1, 0.5, 1.0, 2.0, 5.0]  # Time points
x_values = np.linspace(-5, 5, 200)    # Spatial domain
```

Adjust these parameters according to your problem.

### Example Outputs

After running a script, it will display plots of concentration profiles:

- **Concentration vs. Position**: Shows how concentration varies in space at different times.
- **Concentration vs. Time**: Shows how concentration at a fixed position evolves over time.

---

## Dependencies

The scripts require the following Python libraries:

- **NumPy**: For numerical computations.
- **SciPy**: For special functions and integration.
- **Matplotlib**: For plotting results.

Install them using:

```bash
pip install numpy scipy matplotlib
```

---

## Contributing

Contributions are welcome! If you have suggestions or improvements, please submit a pull request or open an issue.

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Acknowledgments

- **Scientific Computing Resources**: The scripts utilize computational methods and special functions available in scientific libraries.
- **Mathematical References**:
  - Crank, J. *The Mathematics of Diffusion*. Oxford University Press.
  - Palas Kumar Farsoiya, Stéphane Popinet, and Luc Deike. Bubble-mediated transfer of dilute gas in turbulence. Journal of Fluid Mechanics, 920(A34), June 2021.
- **Community Support**: Thanks to the open-source community for providing tools and libraries that make this project possible.
