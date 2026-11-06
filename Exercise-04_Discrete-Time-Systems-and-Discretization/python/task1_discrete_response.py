"""Task 1 — response of x[k+1] = A x[k] + B u[k], A = [[0,1],[-1,0]], B = [1,1],
C = [1 0], x[0] = [1,-1], u[k] = sigma[k].

Closed form: x[k] = [1 - sin(k pi/2), -cos(k pi/2)],  y[k] = 1 - sin(k pi/2)  (period 4).

Run:  python task1_discrete_response.py      (writes figures/task1_response.png)
"""
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from common import check, save

A = sp.Matrix([[0, 1], [-1, 0]]); B = sp.Matrix([1, 1]); C = sp.Matrix([[1, 0]]); x0 = sp.Matrix([1, -1])

xs = [x0]
for n in range(12):
    xs.append(A * xs[-1] + B)
check("x[1] = [0,0], x[2] = [1,1], x[3] = [2,0]", xs[1:4] == [sp.Matrix([0, 0]), sp.Matrix([1, 1]), sp.Matrix([2, 0])])
check("x[4] = x[0]: the state sequence has period 4", xs[4] == x0 and xs[8] == x0)
check("A^2 = -I, A^4 = I", A**2 == -sp.eye(2) and A**4 == sp.eye(2))
check("eigenvalues +-j", set(A.eigenvals()) == {sp.I, -sp.I})
xstar = (sp.eye(2) - A).inv() * B
check("fixed point x* = (I-A)^-1 B = [1,0]", xstar == sp.Matrix([1, 0]))
for n in range(13):
    Ak = sp.cos(n * sp.pi / 2) * sp.eye(2) + sp.sin(n * sp.pi / 2) * A
    assert Ak == A**n
    formula = A**n * x0 + sum((A**(n - 1 - j) * B for j in range(n)), sp.zeros(2, 1))
    closed = sp.Matrix([1 - sp.sin(n * sp.pi / 2), -sp.cos(n * sp.pi / 2)])
    assert formula == xs[n] == closed == A**n * (x0 - xstar) + xstar
check("A^k = cos(k pi/2) I + sin(k pi/2) A, and x[k] = A^k(x0-x*) + x* = [1-sin, -cos], k = 0..12", True)
ys = [(C * x)[0] for x in xs]
check("y = 1, 0, 1, 2, 1, 0, 1, 2, ...", ys[:8] == [1, 0, 1, 2, 1, 0, 1, 2])

K = np.arange(13)
fig, ax = plt.subplots(figsize=(6.4, 3.2))
ax.stem(K, [float(v) for v in ys], basefmt=" ")
ax.set_xlabel("k"); ax.set_ylabel("y[k]"); ax.grid(alpha=0.3)
ax.set_title(r"Task 1: $y[k] = 1 - \sin(k\pi/2)$ (period 4, never settles)", fontsize=10)
save(fig, "task1_response.png")
