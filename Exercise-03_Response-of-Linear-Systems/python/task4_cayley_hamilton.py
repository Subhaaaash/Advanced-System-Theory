"""Task 4 — exp(At) with Cayley-Hamilton, response for a step applied at t0 > 0.

A = [[-6,-1],[5,0]], B = [1,0]^T, C = [1 1], x(t0) = [1,2]^T, u(t) = sigma(t - t0).

Run:  python task4_cayley_hamilton.py      (writes figures/task4_response.png)
"""
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from common import check, save

t, tau, t0 = sp.symbols("t tau t0", real=True)
a0, a1 = sp.symbols("alpha0 alpha1")
A = sp.Matrix([[-6, -1], [5, 0]]); B = sp.Matrix([1, 0]); C = sp.Matrix([[1, 1]]); xt0 = sp.Matrix([1, 2])
lam = sp.Symbol("lambda")

print("(a) Cayley-Hamilton")
check("char. polynomial = (lambda+1)(lambda+5)", sp.factor(A.charpoly(lam).as_expr()) == (lam + 1) * (lam + 5))
check("Cayley-Hamilton: A^2 + 6A + 5I = 0", A**2 + 6 * A + 5 * sp.eye(2) == sp.zeros(2))
sol = sp.solve([a0 - a1 - sp.exp(-t), a0 - 5 * a1 - sp.exp(-5 * t)], [a0, a1])
check("alpha0 = (5e^-t - e^-5t)/4", sp.simplify(sol[a0] - (5 * sp.exp(-t) - sp.exp(-5 * t)) / 4) == 0)
check("alpha1 = (e^-t - e^-5t)/4", sp.simplify(sol[a1] - (sp.exp(-t) - sp.exp(-5 * t)) / 4) == 0)
Phi = sp.simplify(sol[a0] * sp.eye(2) + sol[a1] * A)
Phi_ref = sp.Matrix([[-sp.exp(-t) + 5 * sp.exp(-5 * t), -sp.exp(-t) + sp.exp(-5 * t)],
                     [5 * sp.exp(-t) - 5 * sp.exp(-5 * t), 5 * sp.exp(-t) - sp.exp(-5 * t)]]) / 4
check("exp(At) = 1/4 [[-e^-t+5e^-5t, -e^-t+e^-5t],[5e^-t-5e^-5t, 5e^-t-e^-5t]]", sp.simplify(Phi - Phi_ref) == sp.zeros(2))
check("matches sympy expm", sp.simplify((A * t).exp() - Phi_ref) == sp.zeros(2))

print("(b) response")
v5 = sp.Matrix([1, -1])
check("[1,-1] is the eigenvector of -5 and C [1,-1] = 0 (output never sees e^-5t)", A * v5 == -5 * v5 and (C * v5)[0] == 0)
x_nat = sp.simplify(Phi_ref.subs(t, t - t0) * xt0)
y_nat = sp.simplify((C * x_nat)[0])
check("natural part C e^{A(t-t0)} x(t0) = 3 e^{-(t-t0)}", sp.simplify(y_nat - 3 * sp.exp(-(t - t0))) == 0)
x_forced = sp.simplify(sp.integrate(Phi_ref.subs(t, t - tau) * B, (tau, t0, t)))
check("forced state = A^-1 (e^{A(t-t0)} - I) B", sp.simplify(x_forced - A.inv() * (Phi_ref.subs(t, t - t0) - sp.eye(2)) * B) == sp.zeros(2, 1))
y_forced = sp.simplify((C * x_forced)[0])
check("forced part = 1 - e^{-(t-t0)}", sp.simplify(y_forced - (1 - sp.exp(-(t - t0)))) == 0)
y_tot = sp.simplify(y_nat + y_forced)
check("t >= t0: y(t) = 1 + 2 e^{-(t-t0)}", sp.simplify(y_tot - (1 + 2 * sp.exp(-(t - t0)))) == 0)
T = t - t0
x_tot = sp.Matrix([-sp.exp(-T) / 2 + sp.Rational(3, 2) * sp.exp(-5 * T), 1 + sp.Rational(5, 2) * sp.exp(-T) - sp.Rational(3, 2) * sp.exp(-5 * T)])
check("state for t >= t0 = [-e/2 + 3/2 e5, 1 + 5/2 e - 3/2 e5]", sp.simplify(x_nat + x_forced - x_tot) == sp.zeros(2, 1))
check("continuity at t0: both branches give 3", y_nat.subs(t, t0) == 3 and sp.simplify(y_tot.subs(t, t0)) == 3)
check("steady state = -C A^-1 B = 1", (-C * A.inv() * B)[0] == 1)

# numerical check with t0 = 2: integrate backwards to 0 and forwards to 8
T0 = 2.0
An, Bn, Cn = np.array(A, float), np.array(B, float).ravel(), np.array(C, float).ravel()
back = solve_ivp(lambda _, x: An @ x, (T0, 0), [1, 2], dense_output=True, rtol=1e-10, atol=1e-12)
fwd = solve_ivp(lambda _, x: An @ x + Bn * 1.0, (T0, 8), [1, 2], dense_output=True, rtol=1e-10, atol=1e-12)
tb, tf = np.linspace(0, T0, 200), np.linspace(T0, 8, 400)
yb, yf = Cn @ back.sol(tb), Cn @ fwd.sol(tf)
check("simulation, 0 <= t <= t0: y = 3 e^{-(t-t0)}", np.allclose(yb, 3 * np.exp(-(tb - T0)), rtol=1e-7))
check("simulation, t >= t0: y = 1 + 2 e^{-(t-t0)}", np.allclose(yf, 1 + 2 * np.exp(-(tf - T0)), rtol=1e-7))

fig, ax = plt.subplots(figsize=(6.6, 3.8))
ax.plot(tb, yb, lw=2, label=r"$0\leq t\leq t_0$: $3e^{-(t-t_0)}$ (no input)")
ax.plot(tf, yf, lw=2, label=r"$t\geq t_0$: $1 + 2e^{-(t-t_0)}$ (step on)")
ax.axvline(T0, color="k", ls=":", lw=1); ax.text(T0 + 0.05, 20, "$t_0 = 2$ (example)", fontsize=8)
ax.axhline(1, color="grey", ls="--", lw=1, label="steady state $-CA^{-1}B = 1$")
ax.set_yscale("log"); ax.set_xlabel("t"); ax.set_ylabel("y(t)"); ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=8)
ax.set_title("Task 4(b): output response, step applied at $t_0$", fontsize=10)
save(fig, "task4_response.png")
