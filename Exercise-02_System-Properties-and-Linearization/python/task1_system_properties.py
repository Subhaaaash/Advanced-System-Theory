"""Task 1 — numerical checks of causality, time invariance and linearity.

Each system is written as a function that takes an input signal u (a Python
function of t) and returns the output signal y (again a function of t).
The three definitions of the lecture are then tested on random test signals:

* linearity (Def. 1.1)       : u1 -> y1, u2 -> y2  =>  a1*u1 + a2*u2 -> a1*y1 + a2*y2
* time-invariance (Def. 1.2) : y(u(t), t - T) == y(u(t - T), t)   for T > 0
* causality (Def. 1.3)       : u1(t) = u2(t) for all t <= T  =>  y1(t) = y2(t) for all t <= T

A numerical test can only *disprove* a property (one counterexample is
enough). Passing the test is a sanity check, the proof is in solutions.md.

Run:  python task1_system_properties.py
"""
import numpy as np

rng = np.random.default_rng(0)
A, F0, PHI = 2.0, 0.7, 0.3          # parameters of systems c) and d)
t_grid = np.linspace(-3, 3, 601)


def random_signal():
    """Smooth random test signal: a sum of a few sinusoids."""
    amp = rng.normal(size=4)
    freq = rng.uniform(0.2, 2.0, size=4)
    ph = rng.uniform(0, 2 * np.pi, size=4)
    return lambda t: sum(a * np.sin(2 * np.pi * f * t + p) for a, f, p in zip(amp, freq, ph))


systems = {
    "a) y(t) = u(t) - 0.5 u(t+1)":       lambda u: (lambda t: u(t) - 0.5 * u(t + 1)),
    "b) y(t) = u(2t)":                   lambda u: (lambda t: u(2 * t)),
    "c) y(t) = A cos(2 pi f0 t + u(t))": lambda u: (lambda t: A * np.cos(2 * np.pi * F0 * t + u(t))),
    "d) y(t) = u(t) cos(2 pi f0 t + phi)": lambda u: (lambda t: u(t) * np.cos(2 * np.pi * F0 * t + PHI)),
}


def is_linear(S, trials=20):
    for _ in range(trials):
        u1, u2 = random_signal(), random_signal()
        a, b = rng.normal(size=2)
        lhs = S(lambda t: a * u1(t) + b * u2(t))(t_grid)
        rhs = a * S(u1)(t_grid) + b * S(u2)(t_grid)
        if not np.allclose(lhs, rhs, atol=1e-9):
            return False
    return True


def is_time_invariant(S, trials=20):
    for _ in range(trials):
        u = random_signal()
        T = rng.uniform(0.1, 1.5)
        shifted_input = S(lambda t: u(t - T))(t_grid)      # shift first, then system
        shifted_output = S(u)(t_grid - T)                  # system first, then shift
        if not np.allclose(shifted_input, shifted_output, atol=1e-9):
            return False
    return True


def is_causal(S, trials=20):
    """Build u2 = u1 for t <= T0 and different for t > T0, compare outputs on t <= T0."""
    for _ in range(trials):
        u1, bump = random_signal(), random_signal()
        T0 = rng.uniform(-1, 1)
        u2 = lambda t, u1=u1, bump=bump, T0=T0: u1(t) + np.where(t > T0, 1.0 + bump(t) ** 2, 0.0)
        past = t_grid[t_grid <= T0]
        if not np.allclose(S(u1)(past), S(u2)(past), atol=1e-9):
            return False
    return True


if __name__ == "__main__":
    print(f"{'system':40s} {'causal':>7s} {'TI':>5s} {'linear':>7s}")
    for name, S in systems.items():
        yn = lambda b: "Y" if b else "N"
        print(f"{name:40s} {yn(is_causal(S)):>7s} {yn(is_time_invariant(S)):>5s} {yn(is_linear(S)):>7s}")

    # The concrete counterexamples used in the written solution
    print("\nCounterexamples used in the solutions:")
    # c) zero input does not give zero output -> not linear
    yc = systems["c) y(t) = A cos(2 pi f0 t + u(t))"](lambda t: 0 * t)
    print(f"  c) u = 0 gives y(0) = {yc(0.0):.3f}  (should be 0 for a linear system)")
    # b) shift test with a unit step u(t) = 1 for t >= 0
    step = lambda t: np.where(t >= 0, 1.0, 0.0)
    Sb = systems["b) y(t) = u(2t)"]
    T, t = 1.0, 0.75
    print(f"  b) T=1, t=0.75: S[u(.-T)](t) = u(2t-T) = {Sb(lambda s: step(s - T))(t):.0f},"
          f"  y(t-T) = u(2t-2T) = {Sb(step)(t - T):.0f}  -> time-variant")
