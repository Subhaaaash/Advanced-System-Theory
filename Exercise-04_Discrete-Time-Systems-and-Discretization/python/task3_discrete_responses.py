"""Task 3 — total responses of two discrete systems to impulse and step, x[0] = [1,-1].

System 1: A = diag(-1,-2), B = [1,-1], C = [1 1]
System 2: A = [[0,1],[-2,-3]], B = [0,1], C = [1 0]

Run:  python task3_discrete_responses.py      (writes figures/task3_responses.png)
"""
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from common import check, save

R = sp.Rational
systems = {
    1: (sp.diag(-1, -2), sp.Matrix([1, -1]), sp.Matrix([[1, 1]])),
    2: (sp.Matrix([[0, 1], [-2, -3]]), sp.Matrix([0, 1]), sp.Matrix([[1, 0]])),
}
x0 = sp.Matrix([1, -1])
N = 12


def simulate(A, B, C, u):
    x, ys = x0, []
    for n in range(N):
        ys.append((C * x)[0]); x = A * x + B * u(n)
    return ys


imp = lambda n: 1 if n == 0 else 0
step = lambda n: 1
expected = {
    (1, "zi"): lambda n: (-1)**n - (-2)**n,
    (2, "zi"): lambda n: (-1)**n,
    (1, "imp"): lambda n: 0 if n == 0 else (-2)**(n - 1),
    (2, "imp"): lambda n: 1 if n == 0 else -(-2)**(n - 1),
    (1, "step"): lambda n: R(1, 6) + R(1, 2) * (-1)**n - R(2, 3) * (-2)**n,
    (2, "step"): lambda n: R(1, 6) + R(1, 2) * (-1)**n + R(1, 3) * (-2)**n,
}
res = {}
for i, (A, B, C) in systems.items():
    check(f"system {i}: eigenvalues -1, -2", set(A.eigenvals()) == {-1, -2})
    H = [(C * A**(n - 1) * B)[0] for n in range(1, N)]
    check(f"system {i}: H[k] = (-1)^(k-1) - (-2)^(k-1), k >= 1", H == [(-1)**(n - 1) - (-2)**(n - 1) for n in range(1, N)])
    check(f"system {i}: zero-input response", simulate(A, B, C, lambda n: 0) == [expected[(i, "zi")](n) for n in range(N)])
    for name, u in (("imp", imp), ("step", step)):
        ys = simulate(A, B, C, u); res[(i, name)] = ys
        check(f"system {i}: {name} total response matches closed form", ys == [expected[(i, name)](n) for n in range(N)])
check("system 2: x0 = [1,-1] is the eigenvector of -1", systems[2][0] * x0 == -x0)

K = np.arange(9)
fig, axs = plt.subplots(1, 2, figsize=(10, 3.6))
for ax, name, title in zip(axs, ("imp", "step"), ("impulse $u=\\delta[k]$", "step $u=\\sigma[k]$")):
    for i, m in ((1, "o"), (2, "s")):
        ax.plot(K, [float(v) for v in res[(i, name)][:9]], m + "-", ms=5, lw=0.8, label=f"system {i}")
    ax.set_yscale("symlog", linthresh=2); ax.set_xlabel("k"); ax.grid(alpha=0.3); ax.legend(fontsize=8)
    ax.set_title(f"Task 3: total response, {title}, $x[0]=[1,-1]^T$", fontsize=9)
save(fig, "task3_responses.png")
