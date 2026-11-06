"""Task 5 — y'' + y' = u: state space, exact and Euler discretization, step response
for Ts = 0.5 and Ts = 2*pi.

Closed forms (x[0] = 0, u = sigma[k]):
    exact :  y[k] = k Ts - 1 + (e^{-Ts})^k      (= samples of y(t) = t - 1 + e^{-t})
    Euler :  y[k] = k Ts - 1 + (1 - Ts)^k

Run:  python task5_discretization.py      (writes figures/task5_step_Ts*.png)
"""
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from common import check, save

a, T = sp.symbols("alpha T_s", positive=True)
A = sp.Matrix([[0, 1], [0, -1]]); B = sp.Matrix([0, 1]); C = sp.Matrix([[1, 0]])

print("(a)")
E = sp.simplify((A * a).exp())
check("exp(At) = [[1, 1-e^-t],[0, e^-t]]", sp.simplify(E - sp.Matrix([[1, 1 - sp.exp(-a)], [0, sp.exp(-a)]])) == sp.zeros(2))
Ad = E.subs(a, T)
Bd = sp.simplify(sp.integrate(E * B, (a, 0, T)))
check("Bd = [Ts - 1 + e^-Ts, 1 - e^-Ts]", sp.simplify(Bd - sp.Matrix([T - 1 + sp.exp(-T), 1 - sp.exp(-T)])) == sp.zeros(2, 1))
Ae, Be = sp.eye(2) + T * A, T * B
check("Euler: Ad~ = [[1,Ts],[0,1-Ts]], Bd~ = [0,Ts]", Ae == sp.Matrix([[1, T], [0, 1 - T]]))

print("(b)")
for Tn, name in ((0.5, "0p5"), (2 * np.pi, "2pi")):
    Adn, Bdn = np.array(Ad.subs(T, Tn), float), np.array(Bd.subs(T, Tn), float).ravel()
    Aen, Ben = np.array(Ae.subs(T, Tn), float), np.array(Be.subs(T, Tn), float).ravel()
    K = np.arange(0, int(np.ceil(10 / Tn)) + 1 if Tn < 1 else 6)
    def sim(Ad_, Bd_):
        x, ys = np.zeros(2), []
        for _ in K:
            ys.append(x[0]); x = Ad_ @ x + Bd_
        return np.array(ys)
    ye, yE = sim(Adn, Bdn), sim(Aen, Ben)
    check(f"Ts={Tn:.3g}: exact y[k] = k Ts - 1 + e^(-k Ts)", np.allclose(ye, K * Tn - 1 + np.exp(-K * Tn)))
    check(f"Ts={Tn:.3g}: Euler y[k] = k Ts - 1 + (1-Ts)^k", np.allclose(yE, K * Tn - 1 + (1 - Tn)**K))
    print(f"    Euler eigenvalues: 1 and {1 - Tn:.3f} -> {'stable' if abs(1 - Tn) < 1 else 'UNSTABLE'}")
    tt = np.linspace(0, K[-1] * Tn, 600)
    fig, ax = plt.subplots(figsize=(6.6, 3.8))
    ax.plot(tt, tt - 1 + np.exp(-tt), "k", lw=1.5, label="continuous $y(t) = t - 1 + e^{-t}$")
    ax.plot(K * Tn, ye, "o", ms=6, label="exact discretization")
    ax.plot(K * Tn, yE, "s:", ms=6, mfc="none", lw=1, label="forward Euler")
    if Tn > 1:
        ax.set_yscale("symlog", linthresh=10)
    ax.set_xlabel("t = k $T_s$"); ax.set_ylabel("y"); ax.grid(alpha=0.3); ax.legend(fontsize=8)
    ax.set_title(f"Task 5(b): step response, $T_s$ = {'0.5' if Tn < 1 else '2π'}", fontsize=10)
    save(fig, f"task5_step_Ts{name}.png")
