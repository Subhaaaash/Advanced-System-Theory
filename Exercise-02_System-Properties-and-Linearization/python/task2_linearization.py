"""Task 2 — equilibria and Jacobian linearisation (symbolic + numerical checks).

* Task 2a: equilibria (k*pi, (-1)^k), Jacobians and eigenvalues
* Task 2b: equilibria (0,0) and (0,-1) — symbolic check + a dense numerical
           search that confirms there are no others
* Task 2c: I/O relation  y'' + (3 + y'^2) y' + (1 + y + y^2) u = 0,
           state space and linearisation at the origin
* Figure : nonlinear vs linearised trajectories near the equilibria

Run:  python task2_linearization.py      (writes figures/task2_*.png)
"""
import os

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)

x1, x2, u = sp.symbols("x1 x2 u", real=True)
X = sp.Matrix([x1, x2])

# ---------------------------------------------------------------- Task 2a
print("=" * 70, "\nTask 2a:  x1' = -sin x1,  x2' = -x2 + cos x1")
f_a = sp.Matrix([-sp.sin(x1), -x2 + sp.cos(x1)])
J_a = f_a.jacobian(X)
print("Jacobian df/dx =", J_a.tolist())
for k in range(-2, 3):
    eq = {x1: k * sp.pi, x2: sp.Integer(-1) ** k}
    assert sp.simplify(f_a.subs(eq)) == sp.zeros(2, 1)
    A = J_a.subs(eq)
    print(f"  k={k:+d}: x0 = ({k}*pi, {(-1)**(k % 2):+d})  A = {A.tolist()}  eig = {list(A.eigenvals())}")

# ---------------------------------------------------------------- Task 2b
print("=" * 70, "\nTask 2b")
f_b = sp.Matrix([
    x1**2 + x2**2 + x2 * sp.cos(x1),
    (1 + x1) * x1 + (1 + x2) * x2 + x1 * sp.sin(x2),
])
J_b = f_b.jacobian(X)
print("Jacobian df/dx =", J_b.tolist())
for eq in ({x1: 0, x2: 0}, {x1: 0, x2: -1}):
    assert sp.simplify(f_b.subs(eq)) == sp.zeros(2, 1)
    A = J_b.subs(eq)
    ev = [complex(sp.N(e)) for e in A.eigenvals()]
    print(f"  x0 = ({eq[x1]}, {eq[x2]}):  A = {A.tolist()}")
    print(f"      eigenvalues = {[f'{e.real:+.4f}' for e in ev]}  (exact: {list(A.eigenvals())})")

fb_num = sp.lambdify((x1, x2), list(f_b), "numpy")
found = set()
for a in np.linspace(-3, 3, 121):
    for b in np.linspace(-3, 3, 121):
        s, info, ier, _ = fsolve(lambda z: fb_num(*z), [a, b], full_output=True, xtol=1e-14)
        if ier == 1 and np.max(np.abs(fb_num(*s))) < 1e-11:
            found.add((round(float(s[0]), 8) + 0.0, round(float(s[1]), 8) + 0.0))
print("  numerical search on [-3,3]^2 finds the equilibria:", sorted(found))

# ---------------------------------------------------------------- Task 2c
print("=" * 70, "\nTask 2c:  y'' + (3 + y'^2) y' + (1 + y + y^2) u = 0")
f_c = sp.Matrix([x2, -(3 + x2**2) * x2 - (1 + x1 + x1**2) * u])
g_c = sp.Matrix([x1])
op = {x1: 0, x2: 0, u: 0}
assert f_c.subs(op) == sp.zeros(2, 1)
A = f_c.jacobian(X).subs(op)
B = f_c.jacobian(sp.Matrix([u])).subs(op)
C = g_c.jacobian(X).subs(op)
D = g_c.jacobian(sp.Matrix([u])).subs(op)
print("  A =", A.tolist(), " B =", B.tolist(), " C =", C.tolist(), " D =", D.tolist())
print("  eigenvalues of A:", list(A.eigenvals()))
s = sp.symbols("s")
G = sp.simplify((C * (s * sp.eye(2) - A).inv() * B + D)[0])
print("  transfer function of the linear model G(s) =", G)
print("  1 + x1 + x1^2 = 0 has real roots?", sp.solve(1 + x1 + x1**2, x1), "(complex only -> equilibria need u* = 0)")

# ---------------------------------------------------------------- figure
def compare(f_np, A_np, x0_eq, offsets, T, title, fname):
    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    for i, d in enumerate(offsets):
        z0 = np.array(x0_eq) + np.array(d)
        nl = solve_ivp(lambda t, z: f_np(*z), (0, T), z0, max_step=0.02)
        lin = solve_ivp(lambda t, dz: A_np @ dz, (0, T), np.array(d), max_step=0.02)
        ax.plot(nl.y[0], nl.y[1], color="tab:blue", lw=1.6, label="nonlinear" if i == 0 else None)
        ax.plot(lin.y[0] + x0_eq[0], lin.y[1] + x0_eq[1], "--", color="tab:orange", lw=1.4,
                label="linearised" if i == 0 else None)
        ax.plot(*z0, "k.", ms=5)
    ax.plot(*x0_eq, "r*", ms=12, label="equilibrium")
    ax.set_xlabel("$x_1$"); ax.set_ylabel("$x_2$"); ax.set_title(title, fontsize=10)
    ax.grid(alpha=0.3); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, fname), dpi=140); plt.close(fig)
    print("  wrote", os.path.join("figures", fname))


print("=" * 70, "\nFigures")
fa_np = sp.lambdify((x1, x2), list(f_a), "numpy")
ring = [(0.4 * np.cos(a), 0.4 * np.sin(a)) for a in np.linspace(0, 2 * np.pi, 8, endpoint=False)]
compare(fa_np, np.array([[-1.0, 0], [0, -1]]), (0.0, 1.0), ring, 6,
        "Task 2a, k = 0: stable node at (0, 1)", "task2a_k0_stable.png")
compare(fa_np, np.array([[1.0, 0], [0, -1]]), (np.pi, -1.0), [(0.05 * np.cos(a), 0.4 * np.sin(a)) for a in np.linspace(0, 2*np.pi, 8, endpoint=False)], 2.2,
        "Task 2a, k = 1: saddle at ($\\pi$, -1)", "task2a_k1_saddle.png")
A_b2 = np.array(J_b.subs({x1: 0, x2: -1}).evalf(), dtype=float)
compare(fb_num, A_b2, (0.0, -1.0), [(0.15 * np.cos(a), 0.15 * np.sin(a)) for a in np.linspace(0, 2*np.pi, 8, endpoint=False)], 12,
        "Task 2b: stable node at (0, -1)", "task2b_stable_node.png")
