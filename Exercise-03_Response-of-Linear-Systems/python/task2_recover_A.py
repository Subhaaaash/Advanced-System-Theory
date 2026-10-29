"""Task 2 — recover A and Phi(t,0) from two known free trajectories.

Idea: x(t) = Phi(t) x(0) for every initial state, so stacking the two
trajectories as columns gives  X(t) = Phi(t) X0  =>  Phi(t) = X(t) X0^-1,
and  A = dPhi/dt at t = 0.

Run:  python task2_recover_A.py
"""
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from common import check, save

t = sp.symbols("t", real=True)
e1, e2 = sp.exp(-t), sp.exp(-2 * t)
xa = sp.Matrix([4 * e1 - 3 * e2, -2 * e1 + 3 * e2])   # from x(0) = [1, 1]
xb = sp.Matrix([e2, -e2])                             # from x(0) = [1, -1]
X = sp.Matrix.hstack(xa, xb)
X0 = X.subs(t, 0)
check("X(0) = [[1,1],[1,-1]]", X0 == sp.Matrix([[1, 1], [1, -1]]))
check("X0^-1 = 1/2 [[1,1],[1,-1]]", X0.inv() == sp.Matrix([[1, 1], [1, -1]]) / 2)

Phi = sp.simplify(sp.expand(X * X0.inv()))
Phi_ref = sp.Matrix([[2 * e1 - e2, 2 * e1 - 2 * e2], [-e1 + e2, -e1 + 2 * e2]])
check("Phi(t) = [[2e^-t - e^-2t, 2e^-t - 2e^-2t], [-e^-t + e^-2t, -e^-t + 2e^-2t]]", sp.simplify(Phi - Phi_ref) == sp.zeros(2))
check("Phi(0) = I", Phi_ref.subs(t, 0) == sp.eye(2))

A = sp.simplify(Phi_ref.diff(t).subs(t, 0))
check("A = dPhi/dt(0) = [[0,2],[-1,-3]]", A == sp.Matrix([[0, 2], [-1, -3]]))
check("exp(At) = Phi(t)", sp.simplify((A * t).exp() - Phi_ref) == sp.zeros(2))
check("dPhi/dt = A Phi for all t", sp.simplify(Phi_ref.diff(t) - A * Phi_ref) == sp.zeros(2))
check("both given trajectories solve x' = Ax", all(sp.simplify(x.diff(t) - A * x) == sp.zeros(2, 1) for x in (xa, xb)))
ev = A.eigenvects()
print("  eigen-pairs of A:", [(v, list(vec[0])) for v, _, vec in ev])
check("[1,-1] is the eigenvector of -2 (second trajectory is a pure mode)", sp.simplify(A * sp.Matrix([1, -1]) + 2 * sp.Matrix([1, -1])) == sp.zeros(2, 1))
check("[2,-1] is the eigenvector of -1", sp.simplify(A * sp.Matrix([2, -1]) + sp.Matrix([2, -1])) == sp.zeros(2, 1))

# numerical cross-check: simulate x' = Ax from both initial states
An = np.array(A, float)
tt = np.linspace(0, 5, 200)
fa, fb = sp.lambdify(t, xa, "numpy"), sp.lambdify(t, xb, "numpy")
fig, ax = plt.subplots(figsize=(5.2, 4.4))
for x0, f, name in (([1, 1], fa, "$x(0)=[1,1]^T$"), ([1, -1], fb, "$x(0)=[1,-1]^T$")):
    sim = solve_ivp(lambda _, x: An @ x, (0, 5), x0, t_eval=tt, rtol=1e-10, atol=1e-12)
    given = np.array([np.broadcast_to(np.asarray(v, float).ravel(), tt.shape) for v in f(tt)]).reshape(2, -1)
    check(f"simulation with recovered A reproduces trajectory from {x0}", np.allclose(sim.y, given, atol=1e-7))
    ax.plot(given[0], given[1], lw=2, label=name)
    ax.plot(*x0, "ko", ms=4)
ax.plot(0, 0, "r*", ms=10)
ax.set_xlabel("$x_1$"); ax.set_ylabel("$x_2$"); ax.grid(alpha=0.3); ax.legend()
ax.set_title("Task 2: the two given trajectories (both reach the origin)", fontsize=10)
save(fig, "task2_trajectories.png")
