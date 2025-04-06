#
#  Quantum simulation of lattice wave on a 1d lattice
#  with nearest & next nearest neighbor interactions 
#
# Ref: " Exponential Quantum Speedup for Simulating Classical Lattice Dynamics" 
# 
#
# Xiantao Li, Penn State U, March 2025
#

import numpy as np
import matplotlib.pyplot as plt
from qiskit import Aer
from qiskit.opflow import MatrixEvolution, StateFn, OperatorBase
from qiskit.opflow.primitive_ops import MatrixOp

# Parameters: numer of atoms, step size, total simulation time, and # of time steps
N = 63
dt = 0.05
T = 60
steps = int(T / dt)

# Initial wave packet parameters: center, width and wave number of the wave packet
x0 = N // 2
sigma = 6
k0 = 1.2

# Dispersion relation for initial velocity
omega_k = np.sqrt(2 * (1 - np.cos(k0)) - (1/3) * (1 - np.cos(2*k0)))

# Initial displacement and velocity
x = np.arange(N)
gaussian_envelope = np.exp(-((x - x0) ** 2) / (2 * sigma ** 2))
nu_initial = gaussian_envelope * np.cos(k0 * x)
v_initial = -omega_k * gaussian_envelope * np.sin(k0 * x)

# Construct the matrix Q for the factorization of D
q0 = (1 + 1/np.sqrt(3)) / 2
q1 = -1
q2 = (1 - 1/np.sqrt(3)) / 2

Q = np.zeros((N, N + 2))
for i in range(N):
    Q[i, i] = q0
    Q[i, i + 1] = q1
    Q[i, i + 2] = q2

# construct the Hamiltonian matrix H     
zero_block = np.zeros((N+2, N+2))
H_upper = np.hstack([np.zeros((N, N)), 1j * Q])
H_lower = np.hstack([-1j * Q.T, zero_block])
H = np.vstack([H_upper, H_lower])

# Initial state |psi(0)>; assume m=1
psi = np.concatenate([v_initial, Q.T @ nu_initial])

# Qiskit operator and state preparation
H_op: OperatorBase = MatrixOp(H)
state = StateFn(psi)
evo = MatrixEvolution().convert((H_op * dt).exp_i())

backend = Aer.get_backend('statevector_simulator')

velocity_history = [psi[:N].real.copy()]

# Quantum time evolution
for _ in range(steps):
    state = evo @ state
    psi = np.array(state.eval().primitive).flatten()
    if _ % 20 == 0:
        velocity_history.append(psi[:N].real.copy())

# Plot velocity propagation
plt.figure(figsize=(10, 6))
for i, v in enumerate(velocity_history):
    plt.clf()
    plt.plot(x, v, label=f'Time: {i * 20 * dt:.2f}')
    plt.ylim(-1, 1)
    plt.xlabel('Atom index')
    plt.ylabel('Velocity')
    plt.title('Wave Packet Velocity (Quantum Simulation with Qiskit)')
    plt.legend()
    plt.pause(0.05)

plt.show()
