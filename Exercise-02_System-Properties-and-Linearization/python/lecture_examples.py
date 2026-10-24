"""Lecture 2 examples — symbolic checks.

* Pendulum  m l^2 th'' + b th' + m g l sin th = 0: equilibria (k pi, 0) and
  the linearisation at the hanging (k = 0) and upright (k = 1) positions
* Cart-pendulum (slides 24-26): Jacobians at x0 = 0, u0 = 0

Run:  python lecture_examples.py
"""
import sympy as sp

m, l, g, b, M = sp.symbols("m l g b M", positive=True)
x1, x2, x3, x4, u = sp.symbols("x1 x2 x3 x4 u", real=True)

print("Pendulum")
f = sp.Matrix([x2, -b / (m * l**2) * x2 - g / l * sp.sin(x1)])
J = f.jacobian([x1, x2])
for k in (0, 1):
    A = J.subs({x1: k * sp.pi, x2: 0})
    print(f"  k={k}: A = {A.tolist()}")
    print(f"        char. poly: {sp.expand(A.charpoly().as_expr())}")

print("\nCart-pendulum at x0 = 0, u0 = 0")
den = M + m * sp.sin(x2) ** 2
f = sp.Matrix([
    x3,
    x4,
    (m * l * sp.sin(x2) * x4**2 + m * g * sp.sin(x2) * sp.cos(x2) + u) / den,
    -(m * l * sp.sin(x2) * sp.cos(x2) * x4**2 + (M + m) * g * sp.sin(x2) + sp.cos(x2) * u) / (l * den),
])
op = {x1: 0, x2: 0, x3: 0, x4: 0, u: 0}
A = sp.simplify(f.jacobian([x1, x2, x3, x4]).subs(op))
B = sp.simplify(f.jacobian([u]).subs(op))
print("  A =", A.tolist())
print("  B =", B.tolist())
