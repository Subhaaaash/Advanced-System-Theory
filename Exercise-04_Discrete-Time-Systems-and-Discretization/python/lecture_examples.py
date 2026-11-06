"""Week 4 lecture examples — symbolic/numerical checks of every slide example.

Run:  python lecture_examples.py      (writes figures/lecture_discretization_*.png)
"""
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from scipy.linalg import expm
from common import check, save

k = sp.symbols("k", integer=True, nonnegative=True)
t, Ts = sp.symbols("t T_s", positive=True)

print("Slide 10: A = [[-1,2],[0,1]], A^2 = I")
A = sp.Matrix([[-1, 2], [0, 1]])
check("A^2 = I", A**2 == sp.eye(2))
Phi = sp.Matrix([[(-1)**k, 1 - (-1)**k], [0, 1]])
check("Phi[k] = [[(-1)^k, 1-(-1)^k],[0,1]] for k = 0..9", all(Phi.subs(k, n) == A**n for n in range(10)))
check("x0=[2,1]: odd k -> [0,1], even k -> [2,1]", A * sp.Matrix([2, 1]) == sp.Matrix([0, 1]))

print("Slide 13: impulse response, A = [[0,1],[0,-1]], B=[0,1], C=[1,0]")
A = sp.Matrix([[0, 1], [0, -1]]); B = sp.Matrix([0, 1]); C = sp.Matrix([[1, 0]])
check("A^k = [[0,(-1)^(k-1)],[0,(-1)^k]] for k >= 1", all(A**n == sp.Matrix([[0, (-1)**(n - 1)], [0, (-1)**n]]) for n in range(1, 10)))
H = [(C * A**(n - 1) * B)[0] for n in range(1, 10)]
check("H[1] = 0, H[k] = (-1)^(k-2) for k >= 2", H == [0] + [(-1)**(n - 2) for n in range(2, 10)])

print("Slides 22-29: discretization of A = [[0,-1],[1/2,-1]], B = [0,1]")
A = sp.Matrix([[0, -1], [sp.Rational(1, 2), -1]]); B = sp.Matrix([0, 1])
e, c, s = sp.exp(-t / 2), sp.cos(t / 2), sp.sin(t / 2)
Phi = sp.Matrix([[e * (c + s), -2 * e * s], [e * s, e * (c - s)]])
check("char. poly = lambda^2 + lambda + 1/2", sp.expand(A.charpoly().as_expr()) == sp.Symbol("lambda")**2 + sp.Symbol("lambda") + sp.Rational(1, 2))
check("exp(At) of slide 23 is correct", sp.simplify(sp.expand_complex((A * t).exp()) - Phi) == sp.zeros(2))
check("A^-1 = [[-2,2],[-1,0]]", A.inv() == sp.Matrix([[-2, 2], [-1, 0]]))
Bd = A.inv() * (Phi - sp.eye(2)) * B
check("Bd = [2e(c+s)-2, 2es] (slide 28)", sp.simplify(Bd - sp.Matrix([2 * e * (c + s) - 2, 2 * e * s])) == sp.zeros(2, 1))
check("Bd = integral_0^Ts exp(A a) B da", sp.simplify(sp.integrate(Phi * B, (t, 0, t)) - Bd) == sp.zeros(2, 1))
xt = Phi * sp.Matrix([1, 0]) + Bd
check("total response x(t) = [3e(c+s)-2, 3es] (slide 27)", sp.simplify(xt - sp.Matrix([3 * e * (c + s) - 2, 3 * e * s])) == sp.zeros(2, 1))
An, Bn, T = np.array(A, float), np.array(B, float).ravel(), 0.5
Ad = expm(An * T); Bdn = np.linalg.solve(An, (Ad - np.eye(2)) @ Bn)
check("Ts=0.5: Ad = [[0.947,-0.385],[0.193,0.562]]", np.allclose(Ad, [[0.947, -0.385], [0.193, 0.562]], atol=6e-4))
check("Ts=0.5: Bd = [-0.106, 0.385]", np.allclose(Bdn, [-0.106, 0.385], atol=6e-4))
Ae, Be = np.eye(2) + T * An, T * Bn
check("Euler: Ad~ = [[1,-0.5],[0.25,0.5]], Bd~ = [0,0.5]", np.allclose(Ae, [[1, -0.5], [0.25, 0.5]]) and np.allclose(Be, [0, 0.5]))

tt = np.linspace(0, 10, 400); K = np.arange(0, 21)
fx = sp.lambdify(t, xt, "numpy")
xc = np.array(fx(tt)).reshape(2, -1)
def sim(Ad_, Bd_):
    x = np.array([1.0, 0.0]); out = [x]
    for _ in K[1:]:
        x = Ad_ @ x + Bd_; out.append(x)
    return np.array(out)
xd, xe = sim(Ad, Bdn), sim(Ae, Be)
check("exact discretization = samples of continuous response", np.allclose(xd.T, np.array(fx(K * T)).reshape(2, -1), atol=1e-9))
fig, axs = plt.subplots(2, 1, figsize=(6.6, 5), sharex=True)
for i, ax in enumerate(axs):
    ax.plot(tt, xc[i], "k", lw=1.5, label="continuous $x(t)$")
    ax.plot(K * T, xd[:, i], "o", ms=5, label="exact $A_d, B_d$")
    ax.plot(K * T, xe[:, i], "s", ms=4, mfc="none", label="forward Euler")
    ax.set_ylabel(f"$x_{i+1}$"); ax.grid(alpha=0.3)
axs[0].legend(fontsize=8); axs[1].set_xlabel("t")
axs[0].set_title("Slides 29-30: step response, $x_0=[1,0]^T$, $T_s = 0.5$", fontsize=10)
save(fig, "lecture_discretization_step.png")
