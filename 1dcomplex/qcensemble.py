#
#  Ensemble quantum simulation of dynamics of 1d complex lattice
#
# Ref: " Exponential Quantum Speedup for Simulating Classical Lattice Dynamics" 
# 
#
# Xiantao Li, Penn State U, March 2025
#

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm

# -----------------------------
# Set Matplotlib RCParams 
# -----------------------------
plt.rcParams.update({
    "font.family": "serif",
    "font.size": 14,
    "axes.linewidth": 1.5,
    "lines.linewidth": 1.5,
    "legend.fontsize": 12,
    "xtick.major.size": 6,
    "xtick.minor.size": 3,
    "ytick.major.size": 6,
    "ytick.minor.size": 3,
})

# -----------------------------
# Simulation Parameters
# -----------------------------
N = 127                 # Number of atoms (alternating mass distribution)
num_atoms = N           # Alias for clarity
dt = 0.05               # Time step
num_simulations = 1024  # Number of independent runs
equil_dist = 1.0        # Equilibrium distance between atoms

# Total simulation times to compare.
T_total_values = [2, 48, 96]  # Total simulation time values

# Temperature parameters for the Gaussian bump on the left
T_peak = 2.0           # Maximum temperature at the left boundary
T_base = 0.05           # Baseline (low) temperature for the rest of the chain
sigma_temp = 6         # Width of the Gaussian bump

# -----------------------------
# Masses for Diatomic Chain (alternating 1 and 1.5)
# -----------------------------
masses = np.array([1 if i % 2 == 0 else 1.5 for i in range(num_atoms)], dtype=float)

# Define B matrix from Fejér-Riesz factors; B=Q^T
q0 = 1
q1 = -1
B = np.zeros((N, N + 1))
for i in range(N):
    B[i, i] = q0
    B[i, i + 1] = q1

# Introduce the diagonal mass matrix M and its inverse square root.
M_inv_sqrt = np.diag(1 / np.sqrt(masses))
B = M_inv_sqrt @ B

# Construct Hermitian matrix H using the scaled B.
zero_block = np.zeros((N + 1, N + 1))
H_upper = np.hstack([np.zeros((N, N)), 1j * B])
H_lower = np.hstack([-1j * B.T, zero_block])
H = np.vstack([H_upper, H_lower])

# Precompute the exponential operator for a single time step (dt)
U = expm(-1j * H * dt)

# Container for storing the average temperature profile for each T_total value
avg_temps = {}

for T_total in T_total_values:
    steps = int(T_total / dt)
    temperature_profiles = np.zeros((num_simulations, num_atoms))
    
    for sim in range(num_simulations):
        # Initialize positions: equally spaced with fixed endpoints.
        positions = np.linspace(0, (num_atoms - 1) * equil_dist, num_atoms)
    
        # Define a Gaussian temperature profile on the left.
        indices = np.arange(num_atoms) - 16
        T_profile = T_base + (T_peak - T_base) * np.exp(-indices**2 / (2 * sigma_temp**2))
    
        # Initialize velocities with mass dependence:
        velocities = np.array([np.random.normal(0.0, np.sqrt(T_profile[i] / masses[i]))
                                 for i in range(num_atoms)])
    
        # Initialize displacement for the degrees of freedom (length = num_atoms)
        displacement = np.array([np.random.normal(0.0, np.sqrt(T_profile[i]))
                                  for i in range(num_atoms)])
    
        # Enforce fixed boundaries on velocities (first and last atoms)
        velocities[0] = 0.0
        velocities[-1] = 0.0
    
        # Construct psi vector: first part is velocities, second part is B.T @ displacement.
        psi = np.concatenate([velocities, B.T @ displacement])
    
        # Time evolution using matrix exponential for the given number of steps.
        for _ in range(steps):
            psi = U @ psi
    
        # After time evolution, extract the updated velocities (first num_atoms components)
        velocities_updated = psi[:num_atoms].real
    
        # Compute kinetic energy and local temperature:
        kinetic_energy = 0.5 * masses * velocities_updated**2
        local_temperature = 2 * kinetic_energy  # using 0.5*m*v^2 = 0.5*T (with k_B = 1)
        temperature_profiles[sim, :] = local_temperature
    
    # Store the ensemble-averaged temperature profile.
    avg_temps[T_total] = np.mean(temperature_profiles, axis=0)

# -----------------------------
# Plot the Averaged Temperature Profiles
# -----------------------------
x = np.arange(num_atoms)
plt.figure(figsize=(8, 4))

markers=['o', '+', '^']
for i, T_total in enumerate(T_total_values):
    plt.plot(x, avg_temps[T_total], marker=markers[i], linestyle='-', label=f"T_total = {T_total}")

    

plt.xlabel('Atom Index')
plt.ylabel('Average Kinetic Energy')
plt.legend(frameon=False)
plt.tight_layout()
plt.show()
