"""Task 4 — exact (zero-order hold) and forward-Euler discretization.

(a) double integrator  y'' = u :  A = [[0,1],[0,0]], B = [0,1], C = [1 0]
(b) A = [[-1,0],[1,0]], B = [1,0], C = [0 1]

Both A are singular, so Bd = A^-1 (Ad - I) B cannot be used; Bd = int_0^Ts e^{A a} B da.

Run:  python task4_discretization.py
"""
import numpy as np
import sympy as sp
from scipy.signal import cont2discrete
from common import check

a, T = sp.symbols("alpha T_s", positive=True)


def exact(A, B):
    Ad = sp.simplify((A * T).exp())
    Bd = sp.simplify(sp.integrate((A * a).exp() * B, (a, 0, T)))
    return Ad, Bd


def against_scipy(A, B, C, Ad, Bd, label):
    for Tn in (0.1, 0.5, 2.0, 5.0):
        Adn, Bdn, *_ = cont2discrete((np.array(A, float), np.array(B, float), np.array(C, float), [[0]]), Tn, method="zoh")
        assert np.allclose(Adn, np.array(Ad.subs(T, Tn), float)) and np.allclose(Bdn, np.array(Bd.subs(T, Tn), float))
    check(f"{label}: matches scipy cont2discrete(zoh) for Ts = 0.1, 0.5, 2, 5", True)


print("(a) double integrator")
A = sp.Matrix([[0, 1], [0, 0]]); B = sp.Matrix([0, 1]); C = sp.Matrix([[1, 0]])
check("A nilpotent (A^2 = 0), A singular", A**2 == sp.zeros(2) and A.det() == 0)
Ad, Bd = exact(A, B)
check("Ad = [[1,Ts],[0,1]]", Ad == sp.Matrix([[1, T], [0, 1]]))
check("Bd = [Ts^2/2, Ts]", Bd == sp.Matrix([T**2 / 2, T]))
check("Euler: Ad~ = I + Ts A = Ad (exact, since A^2 = 0), Bd~ = [0, Ts]", sp.eye(2) + T * A == Ad)
against_scipy(A, B, C, Ad, Bd, "(a)")

print("(b)")
A = sp.Matrix([[-1, 0], [1, 0]]); B = sp.Matrix([1, 0]); C = sp.Matrix([[0, 1]])
check("A singular, eigenvalues -1 and 0", A.det() == 0 and set(A.eigenvals()) == {-1, 0})
Ad, Bd = exact(A, B)
check("Ad = [[e^-Ts, 0],[1-e^-Ts, 1]]", sp.simplify(Ad - sp.Matrix([[sp.exp(-T), 0], [1 - sp.exp(-T), 1]])) == sp.zeros(2))
check("Bd = [1-e^-Ts, Ts-1+e^-Ts]", sp.simplify(Bd - sp.Matrix([1 - sp.exp(-T), T - 1 + sp.exp(-T)])) == sp.zeros(2, 1))
check("Euler: Ad~ = [[1-Ts,0],[Ts,1]], Bd~ = [Ts,0]", sp.eye(2) + T * A == sp.Matrix([[1 - T, 0], [T, 1]]))
check("first-order Taylor of exact Ad, Bd equals Euler",
      sp.simplify(sp.Matrix(Ad.applyfunc(lambda f: sp.series(f, T, 0, 2).removeO())) - (sp.eye(2) + T * A)) == sp.zeros(2)
      and sp.simplify(Bd.applyfunc(lambda f: sp.series(f, T, 0, 2).removeO()) - T * B) == sp.zeros(2, 1))
against_scipy(A, B, C, Ad, Bd, "(b)")
