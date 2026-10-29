"""Week 3 lecture examples — symbolic checks of every exp(At) on the slides.

Note: on slide 28 the coefficients alpha_0, alpha_1, alpha_2 are correct, but
the printed Phi(t) is not; the correct matrix is printed below.

Run:  python lecture_examples.py
"""
import sympy as sp
from common import check

t, tau = sp.symbols("t tau", real=True)


def expm(A):
    return sp.simplify((A * t).exp())


print("Similarity transform (slides 12-13): A = [[2,1],[1,2]]")
A = sp.Matrix([[2, 1], [1, 2]])
P = sp.Matrix([[1, 1], [-1, 1]])
L = sp.diag(1, 3)
check("P^-1 A P = diag(1,3)", sp.simplify(P.inv() * A * P - L) == sp.zeros(2))
Phi = P * sp.diag(sp.exp(t), sp.exp(3 * t)) * P.inv()
check("P exp(Lt) P^-1 = exp(At)", sp.simplify(Phi - expm(A)) == sp.zeros(2))
check("x0 = [1,1] -> x(t) = [e^3t, e^3t]", sp.simplify(Phi * sp.Matrix([1, 1]) - sp.Matrix([sp.exp(3 * t)] * 2)) == sp.zeros(2, 1))

print("Nilpotent (slide 14): A = [[1,-1/2],[2,-1]]")
A = sp.Matrix([[1, sp.Rational(-1, 2)], [2, -1]])
check("A^2 = 0", A**2 == sp.zeros(2))
check("exp(At) = I + At", sp.simplify(expm(A) - (sp.eye(2) + A * t)) == sp.zeros(2))

print("Modes (slides 16-17): A = [[-1,1],[0,1]]")
A = sp.Matrix([[-1, 1], [0, 1]])
P = sp.Matrix([[1, 1], [0, 2]])
Pi = P.inv()
modes = P[:, 0] * Pi[0, :] * sp.exp(-t) + P[:, 1] * Pi[1, :] * sp.exp(t)
check("sum of modes = exp(At)", sp.simplify(modes - expm(A)) == sp.zeros(2))
check("x0 = p1 excites only e^-t", sp.simplify(expm(A) * P[:, 0] - P[:, 0] * sp.exp(-t)) == sp.zeros(2, 1))

print("Cayley-Hamilton, complex eigenvalues (slide 24): A = [[0,1],[-1,0]]")
A = sp.Matrix([[0, 1], [-1, 0]])
check("exp(At) = cos t I + sin t A", sp.simplify(expm(A) - (sp.cos(t) * sp.eye(2) + sp.sin(t) * A)) == sp.zeros(2))

print("Cayley-Hamilton, double eigenvalue 0 (slide 25): A = [[0,1],[0,0]]")
A = sp.Matrix([[0, 1], [0, 0]])
check("exp(At) = [[1,t],[0,1]]", expm(A) == sp.Matrix([[1, t], [0, 1]]))

print("Cayley-Hamilton, repeated eigenvalue (slides 27-28)")
A = sp.Matrix([[0, 0, 2], [0, 1, 0], [-1, 0, 3]])
a0 = -2 * t * sp.exp(t) + sp.exp(2 * t)
a1 = (2 + 3 * t) * sp.exp(t) - 2 * sp.exp(2 * t)
a2 = -(t + 1) * sp.exp(t) + sp.exp(2 * t)
Phi = sp.simplify(a0 * sp.eye(3) + a1 * A + a2 * A**2)
check("alpha_0..2 from the slide are correct", sp.simplify(Phi - expm(A)) == sp.zeros(3))
print("  correct Phi(t) =", sp.simplify(sp.expand(Phi)).tolist())

print("Total response (slides 39-41)")
A = sp.Matrix([[0, 1], [0, 0]]); B = sp.Matrix([0, 1]); C = sp.Matrix([[0, 1]]); x0 = sp.Matrix([1, -1])
yh = (C * expm(A) * x0)[0]
yp = sp.integrate((C * expm(A).subs(t, t - tau) * B)[0], (tau, 0, t))
check("natural response y_h = -1", yh == -1)
check("forced step response y_p = t", sp.simplify(yp - t) == 0)
check("impulse response H(t) = 1", (C * expm(A) * B)[0] == 1)
