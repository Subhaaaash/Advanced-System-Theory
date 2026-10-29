# Exercise 3 — Response of Linear Systems

Lecture: Week 03 — Response of linear systems. Problem sheet: Exercise 3 (keep a local copy in [`../course-material/`](../course-material/); it is not published here).

## Files

| File | What it contains |
|---|---|
| [`concepts.md`](concepts.md) · [PDF](pdf/concepts.pdf) | Intuitive overview (what a system response is and why we need it, the lecture in one picture), scalar warm-up, complex numbers and Dirac impulse, homogeneous/inhomogeneous ODEs, state transition matrix, matrix exponential, the four methods (series/nilpotent, diagonalization incl. complex eigenvalues, Jordan form, Cayley–Hamilton incl. repeated eigenvalues), modes, convolution integral, impulse and step response, inverse problems, cheat sheet |
| [`solutions.md`](solutions.md) · [PDF](pdf/solutions.pdf) | Tasks 1–4, fully worked |
| [`python/lecture_examples.py`](python/lecture_examples.py) | All slide examples checked symbolically (incl. the corrected matrix of slide 28) |
| [`python/task1_modes_impulse.py`](python/task1_modes_impulse.py) | Task 1: complex diagonalization, mode matrices, impulse response vs. `scipy.signal.impulse` |
| [`python/task2_recover_A.py`](python/task2_recover_A.py) | Task 2: $\boldsymbol{\Phi} = \boldsymbol{X}(t)\boldsymbol{X}_0^{-1}$, $\boldsymbol{A} = \dot{\boldsymbol{\Phi}}(0)$, simulation check |
| [`python/task3_initial_condition.py`](python/task3_initial_condition.py) | Task 3: Jordan-block exponential, coefficient matching, simulation check |
| [`python/task4_cayley_hamilton.py`](python/task4_cayley_hamilton.py) | Task 4: Cayley–Hamilton coefficients, response before and after $t_0$, backward/forward simulation |

Run any script with `python python/<name>.py` (needs `numpy scipy matplotlib sympy`).

## Progress

- [x] Week 3 concepts
- [x] Task 1 — Modes and responses of LTI systems
- [x] Task 2 — Recovering the state transition matrix
- [x] Task 3 — Determining initial conditions
- [x] Task 4 — Cayley–Hamilton method and response of LTI systems

## Key results

- **Task 1:** $\lambda_{1,2} = -\tfrac12\pm j\tfrac{\sqrt7}{2}$, $e^{\boldsymbol{A}t} = e^{-t/2}\big[\cos\omega t\,\boldsymbol{I} + \tfrac{\sin\omega t}{\sqrt7}(2\boldsymbol{A} + \boldsymbol{I})\big]$, $\omega = \tfrac{\sqrt7}{2}$; no mode diverges; $y(t) = e^{-t/2}\big(\cos\omega t - \tfrac{3}{\sqrt7}\sin\omega t\big)$.
- **Task 2:** $\boldsymbol{A} = \begin{bmatrix}0 & 2\\ -1 & -3\end{bmatrix}$, $\boldsymbol{\Phi}(t,0) = \begin{bmatrix}2e^{-t} - e^{-2t} & 2e^{-t} - 2e^{-2t}\\ -e^{-t} + e^{-2t} & -e^{-t} + 2e^{-2t}\end{bmatrix}$.
- **Task 3:** $\boldsymbol{x}(0) = [-1,\ 1,\ 0]^\top$.
- **Task 4:** $e^{\boldsymbol{A}t} = \tfrac{5e^{-t} - e^{-5t}}{4}\boldsymbol{I} + \tfrac{e^{-t} - e^{-5t}}{4}\boldsymbol{A}$; $y(t) = 3e^{-(t-t_0)}$ for $0\le t\le t_0$, $y(t) = 1 + 2e^{-(t-t_0)}$ for $t\ge t_0$.

![Task 1: impulse response](python/figures/task1_impulse_response.png)
