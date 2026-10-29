"""Task 3 — find x(0) such that the natural response is y(t) = t e^{-t}.

A = [[-1,1,0],[0,-1,0],[0,0,2]] (a 2x2 Jordan block for -1 plus eigenvalue 2),
C = [1 1 1].

Run:  python task3_initial_condition.py
"""
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from common import check, save

t = sp.symbols("t", real=True)
a, b, c = sp.symbols("a b c")
A = sp.Matrix([[-1, 1, 0], [0, -1, 0], [0, 0, 2]])
C = sp.Matrix([[1, 1, 1]])

N = A[:2, :2] + sp.eye(2)
check("Jordan block J = -I + N with N^2 = 0", N**2 == sp.zeros(2))
Phi_ref = sp.Matrix([[sp.exp(-t), t * sp.exp(-t), 0], [0, sp.exp(-t), 0], [0, 0, sp.exp(2 * t)]])
check("exp(At) = [[e^-t, t e^-t, 0],[0, e^-t, 0],[0, 0, e^2t]]", sp.simplify((A * t).exp() - Phi_ref) == sp.zeros(3))

y = sp.expand((C * Phi_ref * sp.Matrix([a, b, c]))[0])
check("y(t) = (a+b) e^-t + b t e^-t + c e^2t",
      sp.simplify(y - ((a + b) * sp.exp(-t) + b * t * sp.exp(-t) + c * sp.exp(2 * t))) == 0)
sol = sp.solve([a + b, b - 1, c], [a, b, c], dict=True)
check("unique solution a=-1, b=1, c=0", sol == [{a: -1, b: 1, c: 0}])
x0 = sp.Matrix([-1, 1, 0])
check("y(t) = t e^-t for x(0) = [-1,1,0]", sp.simplify((C * Phi_ref * x0)[0] - t * sp.exp(-t)) == 0)

tt = np.linspace(0, 8, 400)
sim = solve_ivp(lambda _, x: np.array(A, float) @ x, (0, 8), [-1, 1, 0], t_eval=tt, rtol=1e-10, atol=1e-12)
y_sim = np.array(C, float) @ sim.y
check("simulation gives y = t e^-t", np.allclose(y_sim.ravel(), tt * np.exp(-tt), atol=1e-7))

fig, ax = plt.subplots(figsize=(6.4, 3.6))
ax.plot(tt, y_sim.ravel(), lw=2, label="$y(t)$ from $x(0)=[-1,1,0]^T$")
ax.plot(tt, sim.y[0], "--", label="$x_1(t) = (t-1)e^{-t}$")
ax.plot(tt, sim.y[1], ":", label="$x_2(t) = e^{-t}$")
ax.axhline(0, color="grey", lw=0.6); ax.set_xlabel("t"); ax.grid(alpha=0.3); ax.legend(fontsize=8)
ax.set_title("Task 3: natural response $y(t) = te^{-t}$ ($x_3 = 0$, so $e^{2t}$ is never excited)", fontsize=10)
save(fig, "task3_natural_response.png")
