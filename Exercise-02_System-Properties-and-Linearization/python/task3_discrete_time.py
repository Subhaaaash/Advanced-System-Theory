"""Task 3 — evolution of the discrete-time system x[k+1] = A x[k] + B u[k].

A = [[0, -2], [1, -3]],  B = [1, 1]^T,  y[k] = x[k]

* (a) x[1], x[2], x[3] by iteration, and the closed form
      x[k] = (-1)^k [2, 1]^T - (-2)^k [1, 1]^T            (x[0] = [1, 0]^T)
* (b) eigen-decomposition, A^k in closed form, convergence check
* (c) forced response x[k] = A^k x0 + sum_i (-2)^(k-1-i) u[i] [1, 1]^T
* (d) impulse response x[k] = (-2)^(k-1) [1, 1]^T  (k >= 1)  -> unbounded

Run:  python task3_discrete_time.py      (writes figures/task3_*.png)
"""
import os

import numpy as np
import sympy as sp
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)

A = np.array([[0, -2], [1, -3]])
B = np.array([1, 1])

# (a) iteration -------------------------------------------------------------
x = np.array([1, 0])
print("(a) iteration with u = 0, x[0] = [1, 0]")
traj = [x]
for k in range(1, 4):
    x = A @ x
    traj.append(x)
    print(f"    x[{k}] = {x}")

closed = lambda k: (-1) ** k * np.array([2, 1]) - (-2) ** k * np.array([1, 1])
for k in range(0, 15):
    assert np.array_equal(np.linalg.matrix_power(A, k) @ np.array([1, 0]), closed(k))
print("    closed form (-1)^k [2,1] - (-2)^k [1,1] matches A^k x0 for k = 0..14")

# (b) eigen-decomposition and A^k --------------------------------------------
As = sp.Matrix([[0, -2], [1, -3]])
k = sp.symbols("k", integer=True, nonnegative=True)
print("\n(b) characteristic polynomial:", sp.factor(As.charpoly().as_expr()))
for val, mult, vecs in As.eigenvects():
    print(f"    lambda = {val}: eigenvector {list(vecs[0])}")
Ak = sp.Matrix([[2 * (-1) ** k - (-2) ** k, -2 * (-1) ** k + 2 * (-2) ** k],
                [(-1) ** k - (-2) ** k, -(-1) ** k + 2 * (-2) ** k]])
for kk in range(0, 12):
    assert Ak.subs(k, kk) == As ** kk
print("    A^k closed form verified for k = 0..11:", Ak.tolist())
print("    spectral radius =", max(abs(np.linalg.eigvals(A))), "> 1  -> not internally stable")
v1 = np.array([2, 1])
print("    x0 = c*[2,1] gives x[k] = (-1)^k c [2,1]: bounded, oscillating, NOT convergent:",
      [(np.linalg.matrix_power(A, j) @ v1).tolist() for j in range(4)])

# (c) forced response --------------------------------------------------------
print("\n(c) A^j B = (-2)^j B  (B is an eigenvector for lambda = -2):",
      all(np.array_equal(np.linalg.matrix_power(A, j) @ B, (-2) ** j * B) for j in range(10)))
rng = np.random.default_rng(1)
u = rng.normal(size=10)
x0 = rng.normal(size=2)
x = x0.copy()
for kk in range(10):
    x = A @ x + B * u[kk]
formula = np.linalg.matrix_power(A, 10) @ x0 + sum((-2.0) ** (10 - 1 - i) * u[i] for i in range(10)) * B
print("    simulation vs closed-form forced response at k = 10 agree:", np.allclose(x, formula))

# (d) impulse response -------------------------------------------------------
print("\n(d) impulse input, x[0] = 0")
x = np.zeros(2)
ys = [x.copy()]
for kk in range(10):
    x = A @ x + B * (1.0 if kk == 0 else 0.0)
    ys.append(x.copy())
for kk in range(6):
    print(f"    y[{kk}] = {ys[kk]}")
assert all(np.array_equal(ys[kk], (-2.0) ** (kk - 1) * B) for kk in range(1, 11))
print("    y[k] = (-2)^(k-1) [1,1] for k >= 1  -> |y[k]| = 2^(k-1) grows without bound")

# figure ---------------------------------------------------------------------
K = np.arange(0, 11)
fig, axs = plt.subplots(1, 2, figsize=(10, 3.8))
xa = np.array([np.linalg.matrix_power(A, j) @ np.array([1, 0]) for j in K])
xb = np.array([np.linalg.matrix_power(A, j) @ v1 for j in K])
axs[0].plot(K, xa[:, 0], "o-", label="$x_1[k]$, $x[0]=[1,0]^T$")
axs[0].plot(K, xb[:, 0], "s--", label="$x_1[k]$, $x[0]=[2,1]^T$ (only $\\lambda=-1$ mode)")
axs[0].set_yscale("symlog", linthresh=3)
axs[0].set_title("Free response (u = 0)", fontsize=10)
axs[0].set_xlabel("k"); axs[0].grid(alpha=0.3); axs[0].legend(fontsize=8)
yi = np.array(ys)
axs[1].stem(K, yi[:, 0])
axs[1].set_yscale("symlog", linthresh=3)
axs[1].set_title("Impulse response $y_1[k] = y_2[k]$ ($x[0]=0$)", fontsize=10)
axs[1].set_xlabel("k"); axs[1].grid(alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "task3_responses.png"), dpi=140)
print("\nwrote figures/task3_responses.png")
