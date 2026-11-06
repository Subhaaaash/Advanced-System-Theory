"""Task 2 — reachable states of x[k+1] = A x[k] + B u[k], A = [[0,1],[-0.5,-1]], B = [1,2].

Run:  python task2_reachability.py
"""
import sympy as sp
from common import check

h = sp.Rational(1, 2)
A = sp.Matrix([[0, 1], [-h, -1]]); B = sp.Matrix([1, 2]); lam = sp.Symbol("lambda")

print("(a)")
check("char. poly = lambda^2 + lambda + 1/2", sp.expand(A.charpoly(lam).as_expr()) == lam**2 + lam + h)
check("Cayley-Hamilton: A^2 = -A - 1/2 I", A**2 == -A - h * sp.eye(2))
R = sp.Matrix.hstack(B, A * B)
check("AB = [2, -5/2]", A * B == sp.Matrix([2, -sp.Rational(5, 2)]))
check("det[B, AB] = -13/2 != 0 (full rank)", R.det() == -sp.Rational(13, 2))
for n in range(2, 8):
    c = R.solve(A**n * B)
    assert R * c == A**n * B
check("A^k B lies in span{B, AB} for k = 2..7", True)

print("(b)")
u0, u1 = sp.symbols("u0 u1")
sol = sp.solve(list(A * B * u0 + B * u1 - sp.Matrix([4, 2])), [u0, u1])
check("u[0] = 12/13, u[1] = 28/13", sol == {u0: sp.Rational(12, 13), u1: sp.Rational(28, 13)})
x1 = B * sol[u0]; x2 = A * x1 + B * sol[u1]
check("x[1] = [12/13, 24/13], x[2] = [4,2]", x1 == sp.Matrix([12, 24]) / 13 and x2 == sp.Matrix([4, 2]))
check("one step is not enough: [4,2] is not a multiple of B", sp.Matrix.hstack(B, sp.Matrix([4, 2])).det() != 0)

print("(d) counterexample with x[0] != 0")
A2 = sp.diag(2, 1); B2 = sp.Matrix([0, 1]); x0 = sp.Matrix([1, 0])
check("rank[B, AB] = 1 < n: not reachable", sp.Matrix.hstack(B2, A2 * B2).rank() == 1)
check("x1[k] = 2^k regardless of u, so [8, 0] is reached at k = 3 > n = 2 only",
      all((A2**n * x0)[0] == 2**n for n in range(6)))
