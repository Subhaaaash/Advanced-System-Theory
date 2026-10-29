"""Task 1 — exp(At) by diagonalization, modes, impulse response.

A = [[0,-2],[1,-1]], B = [1,1]^T, C = [1 0], x(0) = 0.
Eigenvalues lambda_{1,2} = -1/2 +- j sqrt(7)/2 (complex pair, Re < 0).

Run:  python task1_modes_impulse.py      (writes figures/task1_*.png)
"""
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from scipy.signal import StateSpace, impulse
from common import check, save

t = sp.symbols("t", real=True)
j = sp.I
w = sp.sqrt(7) / 2
A = sp.Matrix([[0, -2], [1, -1]]); B = sp.Matrix([1, 1]); C = sp.Matrix([[1, 0]])

print("(a) eigenvalues / eigenvectors")
lam1, lam2 = -sp.Rational(1, 2) + j * w, -sp.Rational(1, 2) - j * w
check("char. polynomial = lambda^2 + lambda + 2", sp.expand(A.charpoly().as_expr()) == sp.Symbol("lambda")**2 + sp.Symbol("lambda") + 2)
p1 = sp.Matrix([1 + j * sp.sqrt(7), 2]); p2 = sp.Matrix([1 - j * sp.sqrt(7), 2])
check("A p1 = lambda1 p1", sp.simplify(A * p1 - lam1 * p1) == sp.zeros(2, 1))
check("A p2 = lambda2 p2", sp.simplify(A * p2 - lam2 * p2) == sp.zeros(2, 1))
P = sp.Matrix.hstack(p1, p2); Pi = sp.simplify(P.inv())
print("  P^-1 =", Pi.tolist())
Phi = sp.simplify(sp.expand_complex(P * sp.diag(sp.exp(lam1 * t), sp.exp(lam2 * t)) * Pi))
Phi_ref = sp.exp(-t / 2) * sp.Matrix([
    [sp.cos(w * t) + sp.sin(w * t) / sp.sqrt(7), -4 * sp.sin(w * t) / sp.sqrt(7)],
    [2 * sp.sin(w * t) / sp.sqrt(7), sp.cos(w * t) - sp.sin(w * t) / sp.sqrt(7)],
])
check("P exp(Lambda t) P^-1 = real closed form", sp.simplify(Phi - Phi_ref) == sp.zeros(2))
check("closed form = sympy expm", sp.simplify(sp.expand_complex((A * t).exp()) - Phi_ref) == sp.zeros(2))
check("Phi(0) = I", sp.simplify(Phi_ref.subs(t, 0)) == sp.eye(2))

print("(b) modes")
M1 = sp.simplify(p1 * Pi[0, :]); M2 = sp.simplify(p2 * Pi[1, :])
print("  M1 = p1 p~1 =", M1.tolist())
N = sp.Matrix([[-1, 4], [-2, 1]])
check("M1 = 1/2 (I + j/sqrt7 * [[-1,4],[-2,1]])", sp.simplify(M1 - (sp.eye(2) + j / sp.sqrt(7) * N) / 2) == sp.zeros(2))
check("real form exp(At) = e^{-t/2}[cos(wt) I + sin(wt)(2A+I)/sqrt7]",
      sp.simplify(Phi_ref - sp.exp(-t / 2) * (sp.cos(w * t) * sp.eye(2) + sp.sin(w * t) * (2 * A + sp.eye(2)) / sp.sqrt(7))) == sp.zeros(2))
check("M2 = conj(M1)", sp.simplify(M2 - M1.conjugate()) == sp.zeros(2))
check("M1 + M2 = I", sp.simplify(M1 + M2) == sp.eye(2))
check("M1 e^{l1 t} + M2 e^{l2 t} = exp(At)",
      sp.simplify(sp.expand_complex(M1 * sp.exp(lam1 * t) + M2 * sp.exp(lam2 * t)) - Phi_ref) == sp.zeros(2))
print("  Re(lambda) = -1/2 for both modes -> |e^{lambda t}| = e^{-t/2} -> 0: no mode diverges")

print("(c) impulse response, x(0) = 0")
y = sp.simplify((C * Phi_ref * B)[0])
y_ref = sp.exp(-t / 2) * (sp.cos(w * t) - 3 * sp.sin(w * t) / sp.sqrt(7))
check("y(t) = C exp(At) B = e^{-t/2}(cos wt - 3/sqrt7 sin wt)", sp.simplify(y - y_ref) == 0)
check("y(0) = CB = 1", y_ref.subs(t, 0) == 1)

tt = np.linspace(0, 12, 1200)
_, y_sim = impulse(StateSpace(np.array(A, float), np.array(B, float), np.array(C, float), [[0]]), T=tt)
y_num = sp.lambdify(t, y_ref, "numpy")(tt)
check("analytic y(t) matches scipy impulse()", np.allclose(y_sim, y_num, atol=1e-6))

fig, ax = plt.subplots(figsize=(6.4, 3.6))
ax.plot(tt, y_num, lw=2, label="$y(t)$ (impulse response)")
ax.plot(tt, np.sqrt(16 / 7) * np.exp(-tt / 2), "k--", lw=1, label=r"envelope $\pm\frac{4}{\sqrt{7}}e^{-t/2}$")
ax.plot(tt, -np.sqrt(16 / 7) * np.exp(-tt / 2), "k--", lw=1)
ax.axhline(0, color="grey", lw=0.6); ax.set_xlabel("t"); ax.grid(alpha=0.3); ax.legend()
ax.set_title("Task 1(c): damped oscillation, $\\omega = \\sqrt{7}/2$, decay $e^{-t/2}$", fontsize=10)
save(fig, "task1_impulse_response.png")
