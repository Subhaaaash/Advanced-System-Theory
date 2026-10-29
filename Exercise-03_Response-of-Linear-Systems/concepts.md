---
title: "Exercise 3 — Concepts"
subtitle: "Advanced System Theory · Week 03: Response of Linear Systems"
---

# Exercise 3 — Concepts

**Course:** Advanced System Theory (Paderborn University, Automatic Control Group RAT, Prof. Dr. Erdal Kayacan), winter semester 2026/27
**Lecture:** Week 03 — Response of linear systems
**Covers:** everything needed for Exercise 3 — Task 1 (diagonalization, modes, impulse response), Task 2 (recovering $\boldsymbol{A}$ and $\boldsymbol{\Phi}$ from trajectories), Task 3 (choosing initial conditions), Task 4 (Cayley–Hamilton, response to a delayed step)

This sheet starts with the **intuition** — what a system response is and why we need it — then builds every tool from the ground up. Each method is given as a recipe with a worked example. The full solutions are in [`solutions.md`](solutions.md); the scripts in [`python/`](python/) check every number.

**Notation (as in the lecture slides).** Bold symbols are vectors and matrices: state $\boldsymbol{x}(t)\in\mathbb{R}^n$, input $\boldsymbol{u}(t)$, output $\boldsymbol{y}(t)$, system matrices $\boldsymbol{A},\boldsymbol{B},\boldsymbol{C},\boldsymbol{D}$, forcing term $\boldsymbol{r}(t)$. $\boldsymbol{\Phi}(t, t_0)$ is the state transition matrix, $\exp(\boldsymbol{A}t) = e^{\boldsymbol{A}t}$ the matrix exponential, $\boldsymbol{P}$ the eigenvector matrix and $\boldsymbol{\Lambda} = \operatorname{diag}(\lambda_1,\dots,\lambda_n)$, $\boldsymbol{p}_i$ the $i$-th column of $\boldsymbol{P}$ and $\tilde{\boldsymbol{p}}_i$ the $i$-th row of $\boldsymbol{P}^{-1}$, $\alpha_i(t)$ the Cayley–Hamilton coefficients, $\sigma(t)$ the unit step, $\delta(t)$ the Dirac impulse, $\boldsymbol{H}(t)$ the impulse response matrix. Subscripts $h$ and $p$ mark the homogeneous and particular solution. $j = \sqrt{-1}$.

---

## Contents

0. [The big picture: what this lecture is about](#0-the-big-picture-what-this-lecture-is-about)
1. [Warm-up: the scalar case $\dot x = ax + bu$](#1-warm-up-the-scalar-case-dot-x--ax--bu)
2. [Mathematical toolkit](#2-mathematical-toolkit)
3. [Homogeneous and inhomogeneous equations](#3-homogeneous-and-inhomogeneous-equations)
4. [The state transition matrix](#4-the-state-transition-matrix)
5. [The matrix exponential](#5-the-matrix-exponential)
6. [Method 1 — series expansion and nilpotent matrices](#6-method-1--series-expansion-and-nilpotent-matrices)
7. [Method 2 — similarity transform (diagonalization)](#7-method-2--similarity-transform-diagonalization)
8. [Modes](#8-modes)
9. [Method 3 — Jordan form (non-diagonalizable A)](#9-method-3--jordan-form-non-diagonalizable-a)
10. [Method 4 — Cayley–Hamilton theorem](#10-method-4--cayleyhamilton-theorem)
11. [Which method to choose?](#11-which-method-to-choose)
12. [Response with an input: the convolution integral](#12-response-with-an-input-the-convolution-integral)
13. [Standard test signals, impulse and step response](#13-standard-test-signals-impulse-and-step-response)
14. [Worked lecture example: total response](#14-worked-lecture-example-total-response)
15. [Two "inverse" problems from the exercise](#15-two-inverse-problems-from-the-exercise)
16. [Cheat sheet](#16-cheat-sheet)

---

## 0. The big picture: what this lecture is about

### 0.1 What is a "system response"?

In Week 2 we learned to *write down* a system: $\dot{\boldsymbol{x}} = \boldsymbol{A}\boldsymbol{x} + \boldsymbol{B}\boldsymbol{u}$, $\boldsymbol{y} = \boldsymbol{C}\boldsymbol{x} + \boldsymbol{D}\boldsymbol{u}$. That equation only says **how fast the state is changing right now**. It does not directly say **where the state will be in 5 seconds**.

The **system response** answers exactly that question:

> *Given where the system starts ($\boldsymbol{x}_0$ at time $t_0$) and what we do to it (the input $\boldsymbol{u}(t)$), what are the state $\boldsymbol{x}(t)$ and the output $\boldsymbol{y}(t)$ at every later time?*

In other words: **solving the differential equation**. The state equation is the *rule*; the response is the *movie* that the rule produces.

An everyday picture: a car. The rule is "acceleration = (engine force − friction)/mass". The response is the actual speed-versus-time curve when you start at 30 km/h and press the pedal halfway for 10 seconds.

The response always consists of **two parts** (thanks to the superposition principle of Week 2, Definition 1.1):

| Part | Also called | Cause | Car picture |
|---|---|---|---|
| $\boldsymbol{x}_h$, homogeneous solution | **zero-input / natural / free** response | the initial state alone ($\boldsymbol{u} = \boldsymbol{0}$) | you take your foot off the pedal: the car rolls out from 30 km/h |
| $\boldsymbol{x}_p$, particular solution | **zero-state / forced** response | the input alone ($\boldsymbol{x}_0 = \boldsymbol{0}$) | starting from standstill, you press the pedal |

$$
\text{total response} = \text{natural response} + \text{forced response}.
$$

Because the system is linear, the two can be computed **separately** and simply **added**.

### 0.2 Why do we need it?

- **Prediction and simulation.** Knowing $\boldsymbol{x}(t)$ tells you whether a robot arm reaches its target, whether a tank overflows, how long a motor needs to reach speed.
- **Performance.** Rise time, overshoot, oscillation, settling time — all are read off the response to a standard test signal (step, impulse).
- **Stability.** Does the response die out, stay bounded or blow up? The answer sits in the exponentials $e^{\lambda_i t}$ of the response (formalised in Weeks 9–10).
- **Understanding.** The response decomposes into **modes** — the natural "building blocks" of motion the system can perform (Section 8). Knowing them tells you which initial states or inputs excite slow, fast, oscillating or unstable behaviour.
- **Foundation for the rest of the course.** The state transition matrix $e^{\boldsymbol{A}t}$ reappears in discretization (Week 4: $\boldsymbol{A}_d = e^{\boldsymbol{A}T}$), in the Laplace transform (Week 5: $\mathcal{L}\{e^{\boldsymbol{A}t}\} = (s\boldsymbol{I} - \boldsymbol{A})^{-1}$), in stability (Weeks 9–10), and in controllability/observability (Weeks 11–12).

**The key promise of this lecture (slide 43):** for linear time-invariant (LTI) systems a **general closed-form solution always exists**. For nonlinear systems it usually does not — which is why Week 2 linearised them first.

### 0.3 The lecture in one picture

$$
\boxed{\;\boldsymbol{y}(t) = \underbrace{\boldsymbol{C}\,e^{\boldsymbol{A}(t - t_0)}\boldsymbol{x}_0}_{\text{natural response}} \;+\; \underbrace{\boldsymbol{C}\int_{t_0}^{t} e^{\boldsymbol{A}(t-\tau)}\boldsymbol{B}\,\boldsymbol{u}(\tau)\,\mathrm{d}\tau + \boldsymbol{D}\boldsymbol{u}(t)}_{\text{forced response}}\;}
$$

Everything in the lecture serves this one formula:

| Topic | Intuition in one line | Section |
|---|---|---|
| Homogeneous vs. inhomogeneous ODE | "no input" vs. "with input" — superposition lets us treat them separately | 3 |
| State transition matrix $\boldsymbol{\Phi}(t,t_0)$ | a **time machine** for the state: multiply the state at $t_0$ by $\boldsymbol{\Phi}$ and you get the state at $t$ | 4 |
| Matrix exponential $e^{\boldsymbol{A}t}$ | the matrix version of $e^{at}$ from school — for LTI systems it *is* $\boldsymbol{\Phi}$ | 5 |
| Series / nilpotent | add up $\boldsymbol{I} + \boldsymbol{A}t + \boldsymbol{A}^2t^2/2 + \dots$; sometimes it stops after a few terms | 6 |
| Similarity transform | change coordinates so that the states **decouple** — then each one is a scalar $e^{\lambda t}$ | 7 |
| Modes | the response is a **sum of simple motions** $e^{\lambda_i t}$ in fixed directions $\boldsymbol{p}_i$ — like the notes in a chord | 8 |
| Jordan form | if there are not enough eigenvectors, terms like $t\,e^{\lambda t}$ appear | 9 |
| Cayley–Hamilton | any function of an $n\times n$ matrix is a polynomial of degree $n-1$ — only $n$ unknown coefficients | 10 |
| Convolution integral | the forced response **adds up the after-effects of all past inputs**, each weighted by how long ago it happened | 12 |
| Impulse response $\boldsymbol{H}(t)$ | the system's **fingerprint**: hit it once, record the output; every other response follows from it | 13 |

### 0.4 Which section is needed for which task?

| Exercise task | Needs |
|---|---|
| Task 1(a) diagonalization of $e^{\boldsymbol{A}t}$ | §2.1 complex numbers, §5, §7 (incl. §7.4 complex eigenvalues) |
| Task 1(b) modes, divergence | §8 |
| Task 1(c) impulse response | §2.4 Dirac impulse, §12, §13.2 |
| Task 2 recover $\boldsymbol{A}$ and $\boldsymbol{\Phi}$ | §4.3, §5.2, §15.1 |
| Task 3 initial condition for $y(t) = te^{-t}$ | §2.3 linear independence, §9, §15.2 |
| Task 4(a) Cayley–Hamilton | §10 |
| Task 4(b) response, step at $t_0 > 0$ | §4.2 (backwards in time), §12, §13.3 |

---

## 1. Warm-up: the scalar case $\dot x = ax + bu$

Everything in this lecture is the matrix version of what happens with **one** state. Understanding the scalar case first makes the rest feel natural.

### 1.1 No input: $\dot x = ax$

Try $x(t) = e^{a(t - t_0)}x_0$. Its derivative is $a\,e^{a(t-t_0)}x_0 = a\,x(t)$ ✓, and at $t = t_0$ it gives $e^0x_0 = x_0$ ✓. So

$$
x(t) = e^{a(t - t_0)}\,x_0 .
$$

- $a < 0$: the state decays — the system "forgets" its initial condition.
- $a > 0$: the state grows — unstable.
- $\lvert a\rvert$ sets the **speed**; $1/\lvert a\rvert$ is the **time constant**.

The factor $e^{a(t - t_0)}$ is the scalar "state transition": it carries the state from time $t_0$ to time $t$.

### 1.2 With input: integrating factor

$\dot x - ax = bu$. Multiply by $e^{-at}$; the left side becomes a perfect derivative:

$$
e^{-at}\dot x - a e^{-at}x = \frac{\mathrm{d}}{\mathrm{d}t}\big(e^{-at}x\big) = e^{-at}\,b\,u(t).
$$

Integrate from $t_0$ to $t$:

$$
e^{-at}x(t) - e^{-at_0}x_0 = \int_{t_0}^{t}e^{-a\tau}b\,u(\tau)\,\mathrm{d}\tau
\quad\Longrightarrow\quad
x(t) = \underbrace{e^{a(t - t_0)}x_0}_{\text{natural}} + \underbrace{\int_{t_0}^{t}e^{a(t - \tau)}\,b\,u(\tau)\,\mathrm{d}\tau}_{\text{forced}} .
$$

**Reading the integral.** The input at an earlier time $\tau$ gives the state a little "kick" $b\,u(\tau)\,\mathrm{d}\tau$. By time $t$ that kick has evolved for $t - \tau$ seconds, so it has been multiplied by $e^{a(t - \tau)}$. The integral **adds up all these past kicks**. For $a < 0$, old kicks ($t - \tau$ large) have faded — the system has a **fading memory**.

### 1.3 Step input

$u(t) = 1$ for $t\ge t_0$, $x_0 = 0$:

$$
x(t) = b\int_{t_0}^{t}e^{a(t-\tau)}\,\mathrm{d}\tau = \frac{b}{a}\big(e^{a(t - t_0)} - 1\big)\ \xrightarrow{\ a<0,\ t\to\infty\ }\ -\frac{b}{a}.
$$

The matrix version of exactly this formula is used in Task 4 (§13.3).

**The rest of the lecture replaces the number $a$ by the matrix $\boldsymbol{A}$.** The only real difficulty is: *what is $e^{\boldsymbol{A}t}$, and how do we compute it?*

---

## 2. Mathematical toolkit

### 2.1 Complex numbers and Euler's formula (needed for Task 1)

A complex number is $z = \sigma + j\omega$ with real part $\sigma$, imaginary part $\omega$ and $j^2 = -1$. Its **conjugate** is $\bar z = \sigma - j\omega$.

- Multiplication: $(a + jb)(c + jd) = (ac - bd) + j(ad + bc)$.
- Division: multiply numerator and denominator by the conjugate of the denominator, e.g. $\dfrac{1}{a + jb} = \dfrac{a - jb}{a^2 + b^2}$.
- $z + \bar z = 2\operatorname{Re}z$ (real), $z - \bar z = 2j\operatorname{Im}z$.

**Euler's formula:**

$$
e^{j\omega t} = \cos(\omega t) + j\sin(\omega t),\qquad
\cos(\omega t) = \frac{e^{j\omega t} + e^{-j\omega t}}{2},\qquad
\sin(\omega t) = \frac{e^{j\omega t} - e^{-j\omega t}}{2j}.
$$

So for a complex eigenvalue $\lambda = \sigma + j\omega$:

$$
e^{\lambda t} = e^{\sigma t}e^{j\omega t} = e^{\sigma t}\big(\cos\omega t + j\sin\omega t\big),\qquad \lvert e^{\lambda t}\rvert = e^{\sigma t}.
$$

**Interpretation:** the **real part $\sigma$ decides growth or decay** (envelope $e^{\sigma t}$), the **imaginary part $\omega$ decides the oscillation frequency**. A real matrix always has complex eigenvalues in **conjugate pairs** $\sigma\pm j\omega$, and their contributions combine into real $\cos$/$\sin$ terms.

**Roots of $\lambda^2 + a_1\lambda + a_0 = 0$** with negative discriminant: $\lambda = -\dfrac{a_1}{2} \pm j\sqrt{a_0 - \dfrac{a_1^2}{4}}$.

### 2.2 Matrix facts

- **2×2 inverse:** $\begin{bmatrix}a & b\\ c & d\end{bmatrix}^{-1} = \dfrac{1}{ad - bc}\begin{bmatrix}d & -b\\ -c & a\end{bmatrix}$ (works for complex entries too).
- **Eigenvalues:** roots of $\det(\lambda\boldsymbol{I} - \boldsymbol{A}) = 0$; for 2×2, $\lambda^2 - \operatorname{tr}(\boldsymbol{A})\lambda + \det(\boldsymbol{A}) = 0$.
- **Eigenvectors:** solve $(\lambda_i\boldsymbol{I} - \boldsymbol{A})\boldsymbol{p} = \boldsymbol{0}$; any nonzero multiple is also an eigenvector.
- **Triangular matrices:** eigenvalues are the diagonal entries.
- **Block-diagonal matrices** $\operatorname{diag}(\boldsymbol{A}_1, \boldsymbol{A}_2)$: everything (eigenvalues, powers, exponentials) can be done block by block.
- **Matrix product order matters:** $\boldsymbol{A}\boldsymbol{B}\ne\boldsymbol{B}\boldsymbol{A}$ in general.

### 2.3 Linear independence of exponential functions (needed for Tasks 2 and 3)

The functions $e^{\lambda_1 t}, e^{\lambda_2 t},\dots$ (distinct $\lambda_i$), together with $t\,e^{\lambda t}, t^2e^{\lambda t},\dots$, are **linearly independent**: if

$$
c_1e^{\lambda_1t} + c_2\,t\,e^{\lambda_1 t} + c_3e^{\lambda_2t} + \dots = 0\quad\text{for all } t,
$$

then **every coefficient is zero**. Consequence — **coefficient matching**: two such sums are equal for all $t$ only if the coefficients of each function agree. E.g. $(a + b)e^{-t} + b\,te^{-t} + ce^{2t} = te^{-t}$ for all $t$ forces $a + b = 0$, $b = 1$, $c = 0$.

### 2.4 The unit step and the Dirac impulse

- **Unit step:** $\sigma(t) = 1$ for $t\ge0$, $0$ otherwise. A delayed step $\sigma(t - t_0)$ switches on at $t_0$.
- **Rectangle:** $\Pi(t) = 1/\epsilon$ for $0\le t\le\epsilon$, $0$ otherwise — area always $1$.
- **Dirac impulse:** $\delta(t) = \lim_{\epsilon\to0}\Pi(t)$ — infinitely short, infinitely high, **area 1**. It is not an ordinary function; it is defined by what it does inside an integral, the **sifting property**:

$$
\int_{a}^{b} f(\tau)\,\delta(\tau - t_1)\,\mathrm{d}\tau = f(t_1)\qquad\text{if } a < t_1 < b .
$$

The impulse "picks out" the value of $f$ at the instant of the impulse. Physical picture: a hammer blow — a huge force for a very short time, delivering a finite momentum.

**Convention for an impulse at $t = 0$:** the impulse acts at the very beginning, so the integral $\int_0^t$ is understood as starting just before it ($0^-$) and the whole area counts: $\int_{0^-}^{t} f(\tau)\delta(\tau)\,\mathrm{d}\tau = f(0)$ for $t > 0$.

### 2.5 Integrating exponentials

$\displaystyle\int e^{c\tau}\,\mathrm{d}\tau = \frac{1}{c}e^{c\tau}$ ($c\ne0$), and therefore $\displaystyle\int_{t_0}^{t}e^{a(t - \tau)}\,\mathrm{d}\tau = \frac{1}{a}\big(e^{a(t-t_0)} - 1\big)$. For matrices: §13.3.

---

## 3. Homogeneous and inhomogeneous equations

### 3.1 Definition 1.1

We consider the system of $n$ linear **homogeneous** ODEs

$$
\dot{\boldsymbol{x}}(t) = \boldsymbol{A}(t)\boldsymbol{x}(t)
$$

and the system of $n$ linear **inhomogeneous** ODEs

$$
\dot{\boldsymbol{x}}(t) = \boldsymbol{A}(t)\boldsymbol{x}(t) + \boldsymbol{r}(t),
$$

where $\boldsymbol{A}(t)$ is the dynamic matrix, $\boldsymbol{r}(t)\in\mathbb{R}^n$ the system input (forcing term) and $\boldsymbol{x}(t)\in\mathbb{R}^n$ the state. For specified initial conditions $\boldsymbol{x}(t_0) = \boldsymbol{x}_0$ these systems possess **unique** solutions for every $(t_0, \boldsymbol{x}_0)$.

For a state-space model the forcing term is $\boldsymbol{r}(t) = \boldsymbol{B}\boldsymbol{u}(t)$ (slide 33).

### 3.2 Structure of the solution (slide 6)

$$
\boldsymbol{x}(t, t_0, \boldsymbol{x}_0) = \underbrace{\boldsymbol{\Phi}(t, t_0)\boldsymbol{x}_0}_{\boldsymbol{x}_h} + \underbrace{\int_{t_0}^{t}\boldsymbol{\Phi}(t, \tau)\boldsymbol{r}(\tau)\,\mathrm{d}\tau}_{\boldsymbol{x}_p}.
$$

- **Homogeneous solution $\boldsymbol{x}_h$:** nonzero initial conditions, zero forcing ⇒ **zero-input response**.
- **Particular solution $\boldsymbol{x}_p$:** zero initial conditions, nonzero forcing ⇒ **zero-state response**.
- For given $\boldsymbol{x}_0$ and given forcing, $\boldsymbol{x}(t)$ is known for **all** $t$. *This is why we call $\boldsymbol{x}$ the state:* it contains everything about the past that matters for the future.

---

## 4. The state transition matrix

### 4.1 Definition 1.2

The **state transition matrix** $\boldsymbol{\Phi}(t, t_0)$ is the unique (nonsingular) solution of the matrix equation

$$
\dot{\boldsymbol{\Phi}}(t, t_0) = \boldsymbol{A}(t)\boldsymbol{\Phi}(t, t_0)\qquad\text{with}\qquad \boldsymbol{\Phi}(t_0, t_0) = \boldsymbol{I}.
$$

It is uniquely determined by $\boldsymbol{A}(t)$, and the solution of the homogeneous ODE is

$$
\boldsymbol{x}(t, t_0, \boldsymbol{x}_0) = \boldsymbol{\Phi}(t, t_0)\,\boldsymbol{x}_0\qquad\forall t.
$$

**Intuition — a time machine for the state.** Give $\boldsymbol{\Phi}(t, t_0)$ the state at time $t_0$ and it returns the state at time $t$. *Why does $\boldsymbol{\Phi}$ satisfy the same ODE as $\boldsymbol{x}$?* Because each **column** of $\boldsymbol{\Phi}$ is itself a solution: the $i$-th column is the trajectory that starts at the $i$-th unit vector $\boldsymbol{e}_i$ (since $\boldsymbol{\Phi}(t_0,t_0)\boldsymbol{e}_i = \boldsymbol{e}_i$). Any other initial state is a combination $\boldsymbol{x}_0 = \sum x_{0,i}\boldsymbol{e}_i$, so by superposition its trajectory is the same combination of the columns — that is $\boldsymbol{\Phi}\boldsymbol{x}_0$.

### 4.2 Properties (slide 8)

For any $t, \tau, \sigma\in\mathbb{R}$:

1. **Composition:** $\boldsymbol{\Phi}(t, \tau) = \boldsymbol{\Phi}(t, \sigma)\boldsymbol{\Phi}(\sigma, \tau)$ — travelling from $\tau$ to $\sigma$ and then from $\sigma$ to $t$ is the same as travelling from $\tau$ to $t$ directly.
2. **Invertibility:** $\det\boldsymbol{\Phi}(t, t_0)\ne0$ and
   $$\big(\boldsymbol{\Phi}(t, t_0)\big)^{-1} = \boldsymbol{\Phi}(t_0, t).$$
   The inverse travels **backwards**. In particular, $\boldsymbol{\Phi}(t, t_0)$ with $t < t_0$ is perfectly meaningful: if you know the state at $t_0$, you can compute the state at **earlier** times. *(Needed in Task 4, where the state is given at $t_0 > 0$ and the response is asked for $0\le t\le t_0$.)*
3. Closed-form expressions of $\boldsymbol{\Phi}(t, t_0)$ exist only for special classes of $\boldsymbol{A}(t)$ — most importantly for **constant** $\boldsymbol{A}$ (next section).

### 4.3 Building $\boldsymbol{\Phi}$ from known trajectories (needed for Task 2)

Suppose you know $n$ trajectories $\boldsymbol{x}^{(1)}(t),\dots,\boldsymbol{x}^{(n)}(t)$ of $\dot{\boldsymbol{x}} = \boldsymbol{A}\boldsymbol{x}$, starting from **linearly independent** initial states $\boldsymbol{x}^{(i)}(0)$. Each satisfies $\boldsymbol{x}^{(i)}(t) = \boldsymbol{\Phi}(t, 0)\boldsymbol{x}^{(i)}(0)$. Stack them as columns:

$$
\underbrace{\big[\boldsymbol{x}^{(1)}(t)\ \cdots\ \boldsymbol{x}^{(n)}(t)\big]}_{\boldsymbol{X}(t)} = \boldsymbol{\Phi}(t, 0)\underbrace{\big[\boldsymbol{x}^{(1)}(0)\ \cdots\ \boldsymbol{x}^{(n)}(0)\big]}_{\boldsymbol{X}_0}
\quad\Longrightarrow\quad
\boxed{\boldsymbol{\Phi}(t, 0) = \boldsymbol{X}(t)\,\boldsymbol{X}_0^{-1}}
$$

$\boldsymbol{X}_0$ is invertible because its columns are independent. ($\boldsymbol{X}(t)$ is called a **fundamental matrix**.) Then $\boldsymbol{A}$ follows from Definition 1.2 at $t = 0$: $\dot{\boldsymbol{\Phi}}(0,0) = \boldsymbol{A}\boldsymbol{\Phi}(0,0) = \boldsymbol{A}$, i.e.

$$
\boxed{\boldsymbol{A} = \dot{\boldsymbol{\Phi}}(t, 0)\big|_{t=0}.}
$$

---

## 5. The matrix exponential

### 5.1 Theorem 1.1 — state transition matrix for autonomous LTI systems

For time-invariant $\boldsymbol{A}$, the state transition matrix is

$$
\boldsymbol{\Phi}(t, t_0) = \exp\big[\boldsymbol{A}(t - t_0)\big],\qquad\text{in particular}\qquad
\boldsymbol{\Phi}(t, 0) = \exp(\boldsymbol{A}t) = \sum_{k=0}^{\infty}\frac{t^k}{k!}\boldsymbol{A}^k = \boldsymbol{I} + \boldsymbol{A}t + \frac{\boldsymbol{A}^2t^2}{2!} + \frac{\boldsymbol{A}^3t^3}{3!} + \cdots
$$

This is the same power series as the scalar $e^{at} = 1 + at + \frac{a^2t^2}{2!} + \cdots$, with $a$ replaced by $\boldsymbol{A}$ and $1$ by $\boldsymbol{I}$. The series always converges.

**Proof (slide 9)** — differentiate term by term:

$$
\frac{\mathrm{d}}{\mathrm{d}t}\exp(\boldsymbol{A}t) = \sum_{k=1}^{\infty}\frac{k\,t^{k-1}}{k\,(k-1)!}\boldsymbol{A}^k = \sum_{k=0}^{\infty}\frac{t^k}{k!}\boldsymbol{A}^{k+1} = \boldsymbol{A}\exp(\boldsymbol{A}t),
$$

so $\dot{\boldsymbol{\Phi}} = \boldsymbol{A}\boldsymbol{\Phi}$; and at $t = 0$ only the $k = 0$ term survives: $\exp(\boldsymbol{0}) = \boldsymbol{I}$, $\det\boldsymbol{I} = 1\ne0$. Both conditions of Definition 1.2 hold.

### 5.2 Useful properties

| Property | Formula | Comment |
|---|---|---|
| start value | $e^{\boldsymbol{A}\cdot0} = \boldsymbol{I}$ | always check this first |
| derivative | $\frac{\mathrm{d}}{\mathrm{d}t}e^{\boldsymbol{A}t} = \boldsymbol{A}e^{\boldsymbol{A}t} = e^{\boldsymbol{A}t}\boldsymbol{A}$ | at $t = 0$ this gives $\boldsymbol{A}$ (Task 2) |
| addition of times | $e^{\boldsymbol{A}(t+s)} = e^{\boldsymbol{A}t}e^{\boldsymbol{A}s}$ | the composition property for LTI systems |
| inverse | $\big(e^{\boldsymbol{A}t}\big)^{-1} = e^{-\boldsymbol{A}t}$ | backwards in time |
| diagonal $\boldsymbol{A}$ | $e^{\operatorname{diag}(\lambda_i)t} = \operatorname{diag}(e^{\lambda_it})$ | element by element **only** for diagonal matrices |
| block-diagonal | $e^{\operatorname{diag}(\boldsymbol{A}_1,\boldsymbol{A}_2)t} = \operatorname{diag}(e^{\boldsymbol{A}_1t}, e^{\boldsymbol{A}_2t})$ | used in Task 3 |
| warning | $e^{(\boldsymbol{A}+\boldsymbol{B})t}\ne e^{\boldsymbol{A}t}e^{\boldsymbol{B}t}$ in general | equality only if $\boldsymbol{A}\boldsymbol{B} = \boldsymbol{B}\boldsymbol{A}$ |

**Common mistake:** $e^{\boldsymbol{A}t}$ is **not** the matrix of $e^{a_{ij}t}$. For $\boldsymbol{A} = \begin{bmatrix}0 & 1\\ 0 & 0\end{bmatrix}$ the element-wise exponential would be $\begin{bmatrix}1 & e^t\\ 1 & 1\end{bmatrix}$, but the true one is $\begin{bmatrix}1 & t\\ 0 & 1\end{bmatrix}$.

Computing the infinite series by hand is impractical, so the lecture gives four methods (Sections 6, 7, 9, 10).

---

## 6. Method 1 — series expansion and nilpotent matrices

A matrix is **nilpotent** if $\boldsymbol{A}^k = \boldsymbol{0}$ for some $k$. Then all later powers vanish too and the series **terminates**:

$$
\boldsymbol{A}^2 = \boldsymbol{0}\quad\Longrightarrow\quad e^{\boldsymbol{A}t} = \boldsymbol{I} + \boldsymbol{A}t .
$$

**Slide 14:** $\boldsymbol{A} = \begin{bmatrix}1 & -1/2\\ 2 & -1\end{bmatrix}$, $\boldsymbol{A}^2 = \begin{bmatrix}1 - 1 & -\frac12 + \frac12\\ 2 - 2 & -1 + 1\end{bmatrix} = \boldsymbol{0}$, so

$$
\boldsymbol{\Phi}(t, t_0) = \boldsymbol{I} + (t - t_0)\boldsymbol{A} = \begin{bmatrix}1 + (t - t_0) & -\frac12(t - t_0)\\ 2(t - t_0) & 1 - (t - t_0)\end{bmatrix}.
$$

A nilpotent matrix has **all eigenvalues equal to 0** (if $\boldsymbol{A}\boldsymbol{v} = \lambda\boldsymbol{v}$ then $\boldsymbol{0} = \boldsymbol{A}^k\boldsymbol{v} = \lambda^k\boldsymbol{v}$, so $\lambda = 0$). For non-nilpotent matrices, truncating the series gives only an approximation, useful for small $t$.

---

## 7. Method 2 — similarity transform (diagonalization)

### 7.1 Idea: change coordinates so the states decouple

If $\boldsymbol{A}$ were diagonal, every state would evolve on its own: $\dot x_i = \lambda_ix_i$, so $x_i(t) = e^{\lambda_it}x_i(0)$ — $n$ independent scalar problems. Most matrices are not diagonal, but many become diagonal **in the right coordinates**: the coordinates along the eigenvectors.

Two observations (slide 11):

1. If $\boldsymbol{A}' = \boldsymbol{P}^{-1}\boldsymbol{A}\boldsymbol{P}$ ($\boldsymbol{A}'$ and $\boldsymbol{A}$ are **similar**), then $\exp(\boldsymbol{A}'t) = \boldsymbol{P}^{-1}\exp(\boldsymbol{A}t)\boldsymbol{P}$.
   *Reason:* $(\boldsymbol{P}^{-1}\boldsymbol{A}\boldsymbol{P})^k = \boldsymbol{P}^{-1}\boldsymbol{A}^k\boldsymbol{P}$ (all inner $\boldsymbol{P}\boldsymbol{P}^{-1}$ cancel), so every term of the series transforms the same way.
2. If $\boldsymbol{A}$ is **diagonalizable**, i.e. $\boldsymbol{P}^{-1}\boldsymbol{A}\boldsymbol{P} = \boldsymbol{\Lambda} = \operatorname{diag}(\lambda_1,\dots,\lambda_n)$, then

$$
\boxed{\exp(\boldsymbol{A}t) = \boldsymbol{P}\exp(\boldsymbol{\Lambda}t)\boldsymbol{P}^{-1} = \boldsymbol{P}\operatorname{diag}\big(e^{\lambda_1t},\dots,e^{\lambda_nt}\big)\boldsymbol{P}^{-1}.}
$$

- $\boldsymbol{A}$ is diagonalizable if it has $n$ linearly independent eigenvectors $\boldsymbol{u}_i$; then $\boldsymbol{P} = [\boldsymbol{u}_1\ \cdots\ \boldsymbol{u}_n]$. **Distinct eigenvalues always give independent eigenvectors.**
- Matrices with repeated eigenvalues *may* not be diagonalizable (→ Jordan form, §9). Real symmetric matrices always are.

**Reading the formula right to left:** $\boldsymbol{P}^{-1}$ converts the initial state into eigen-coordinates, $e^{\boldsymbol{\Lambda}t}$ lets each eigen-coordinate evolve on its own, $\boldsymbol{P}$ converts back.

### 7.2 Recipe

1. Eigenvalues: solve $\det(\lambda\boldsymbol{I} - \boldsymbol{A}) = 0$.
2. Eigenvector for each $\lambda_i$: solve $(\lambda_i\boldsymbol{I} - \boldsymbol{A})\boldsymbol{u}_i = \boldsymbol{0}$.
3. $\boldsymbol{P} = [\boldsymbol{u}_1\ \cdots\ \boldsymbol{u}_n]$, compute $\boldsymbol{P}^{-1}$.
4. Multiply $\boldsymbol{P}\operatorname{diag}(e^{\lambda_it})\boldsymbol{P}^{-1}$.
5. **Check:** at $t = 0$ the result must be $\boldsymbol{I}$; its derivative at $t = 0$ must be $\boldsymbol{A}$.

### 7.3 Worked lecture example (slides 12–13)

$\boldsymbol{A} = \begin{bmatrix}2 & 1\\ 1 & 2\end{bmatrix}$: $\det(\lambda\boldsymbol{I} - \boldsymbol{A}) = (\lambda - 2)^2 - 1 = (\lambda - 1)(\lambda - 3)$, so $\lambda_1 = 1$, $\lambda_2 = 3$.

- $\lambda_1 = 1$: $\lambda_1\boldsymbol{I} - \boldsymbol{A} = \begin{bmatrix}-1 & -1\\ -1 & -1\end{bmatrix}$ ⇒ $u_1 + u_2 = 0$ ⇒ $\boldsymbol{u}_1 = [1, -1]^\top$.
- $\lambda_2 = 3$: $\lambda_2\boldsymbol{I} - \boldsymbol{A} = \begin{bmatrix}1 & -1\\ -1 & 1\end{bmatrix}$ ⇒ $u_1 = u_2$ ⇒ $\boldsymbol{u}_2 = [1, 1]^\top$.

$$
\boldsymbol{P} = \begin{bmatrix}1 & 1\\ -1 & 1\end{bmatrix},\quad \boldsymbol{\Lambda} = \begin{bmatrix}1 & 0\\ 0 & 3\end{bmatrix},\quad \boldsymbol{P}^{-1} = \frac12\begin{bmatrix}1 & -1\\ 1 & 1\end{bmatrix},
$$

$$
\boldsymbol{\Phi}(t, t_0) = \boldsymbol{P}\begin{bmatrix}e^{t - t_0} & 0\\ 0 & e^{3(t - t_0)}\end{bmatrix}\boldsymbol{P}^{-1} = \frac12\begin{bmatrix}e^{t - t_0} + e^{3(t - t_0)} & -e^{t - t_0} + e^{3(t - t_0)}\\ -e^{t - t_0} + e^{3(t - t_0)} & e^{t - t_0} + e^{3(t - t_0)}\end{bmatrix}.
$$

For $t_0 = 0$, $\boldsymbol{x}_0 = [1, 1]^\top$: $\boldsymbol{x}(t) = [e^{3t}, e^{3t}]^\top$ — the initial state is exactly the eigenvector $\boldsymbol{u}_2$, so only $e^{3t}$ appears (a preview of modes, §8).

### 7.4 Complex eigenvalues (needed for Task 1)

If $\det(\lambda\boldsymbol{I} - \boldsymbol{A}) = 0$ has complex roots $\lambda_{1,2} = \sigma\pm j\omega$, **the recipe works unchanged** — the eigenvectors and $\boldsymbol{P}$ are complex, but the final $e^{\boldsymbol{A}t}$ is real. Practical points:

- **Eigenvectors come in conjugate pairs:** if $\boldsymbol{u}_1$ belongs to $\lambda_1$, then $\boldsymbol{u}_2 = \bar{\boldsymbol{u}}_1$ belongs to $\lambda_2 = \bar\lambda_1$. Compute one, conjugate for the other.
- **Finding $\boldsymbol{u}_1$:** $\lambda_1\boldsymbol{I} - \boldsymbol{A}$ is singular, so one row suffices — pick the simpler row.
- **Simplify with Euler:** the result contains $e^{\lambda_1t}$ and $e^{\lambda_2t}$ in conjugate pairs; group them as $e^{\sigma t}\cos\omega t$ and $e^{\sigma t}\sin\omega t$ (§2.1).

**Shortcut for any real 2×2 matrix with eigenvalues $\sigma\pm j\omega$** ($\omega\ne0$):

$$
\boxed{e^{\boldsymbol{A}t} = e^{\sigma t}\Big[\cos(\omega t)\,\boldsymbol{I} + \frac{\sin(\omega t)}{\omega}\big(\boldsymbol{A} - \sigma\boldsymbol{I}\big)\Big]}
$$

*Check:* at $t = 0$ it gives $\boldsymbol{I}$; its derivative at $t = 0$ is $\sigma\boldsymbol{I} + (\boldsymbol{A} - \sigma\boldsymbol{I}) = \boldsymbol{A}$ ✓. (It follows from Cayley–Hamilton, §10.) Use it to **verify** the diagonalization result. *Example (slide 24):* $\boldsymbol{A} = \begin{bmatrix}0 & 1\\ -1 & 0\end{bmatrix}$, $\sigma = 0$, $\omega = 1$: $e^{\boldsymbol{A}t} = \cos t\,\boldsymbol{I} + \sin t\,\boldsymbol{A} = \begin{bmatrix}\cos t & \sin t\\ -\sin t & \cos t\end{bmatrix}$ — a rotation: the state circles the origin forever.

---

## 8. Modes

### 8.1 The modal decomposition (slide 15)

Write $\boldsymbol{P} = [\boldsymbol{p}_1\ \cdots\ \boldsymbol{p}_n]$ (columns = eigenvectors) and $\boldsymbol{P}^{-1}$ with rows $\tilde{\boldsymbol{p}}_1,\dots,\tilde{\boldsymbol{p}}_n$. Multiplying out $\boldsymbol{P}\exp(\boldsymbol{\Lambda}t)\boldsymbol{P}^{-1}$ column-times-row gives a **sum of $n$ modes**:

$$
\boxed{\exp(\boldsymbol{A}t) = \sum_{i=1}^{n}\boldsymbol{p}_i\tilde{\boldsymbol{p}}_i\,e^{\lambda_it}}
$$

Each $\boldsymbol{p}_i\tilde{\boldsymbol{p}}_i$ is an $n\times n$ matrix (column × row). **Each term is called a mode** of $\dot{\boldsymbol{x}} = \boldsymbol{A}\boldsymbol{x}$. In general the solution contains contributions of all modes.

Useful facts about the mode matrices $\boldsymbol{M}_i = \boldsymbol{p}_i\tilde{\boldsymbol{p}}_i$:

- $\tilde{\boldsymbol{p}}_i\boldsymbol{p}_i = 1$ and $\tilde{\boldsymbol{p}}_j\boldsymbol{p}_i = 0$ for $i\ne j$ (because $\boldsymbol{P}^{-1}\boldsymbol{P} = \boldsymbol{I}$).
- $\sum_i\boldsymbol{M}_i = \boldsymbol{I}$ (set $t = 0$) — a quick check.
- For a complex pair, $\boldsymbol{M}_2 = \bar{\boldsymbol{M}}_1$, so the two modes add up to a real term $2\operatorname{Re}\big(\boldsymbol{M}_1e^{\lambda_1t}\big)$.

### 8.2 Exciting a single mode (slides 16–17)

If the initial state lies along one eigenvector, $\boldsymbol{x}(0) = \alpha\boldsymbol{p}_i$, then only that mode appears:

$$
\boldsymbol{x}(t) = \exp(\boldsymbol{A}t)\boldsymbol{x}(0) = \sum_k\boldsymbol{p}_k\underbrace{\tilde{\boldsymbol{p}}_k\boldsymbol{p}_i}_{1\text{ if }k = i,\ 0\text{ else}}\alpha\,e^{\lambda_kt} = \alpha\boldsymbol{p}_ie^{\lambda_it}.
$$

The state moves **along the straight line through $\boldsymbol{p}_i$**, growing or shrinking like $e^{\lambda_it}$.

*Lecture example:* $\boldsymbol{A} = \begin{bmatrix}-1 & 1\\ 0 & 1\end{bmatrix}$, $\lambda_1 = -1$, $\lambda_2 = 1$, $\boldsymbol{P} = \begin{bmatrix}1 & 1\\ 0 & 2\end{bmatrix}$, $\boldsymbol{P}^{-1} = \begin{bmatrix}1 & -1/2\\ 0 & 1/2\end{bmatrix}$:

$$
\boldsymbol{\Phi}(t) = \begin{bmatrix}1 & -1/2\\ 0 & 0\end{bmatrix}e^{-t} + \begin{bmatrix}0 & 1/2\\ 0 & 1\end{bmatrix}e^{t}.
$$

For $\boldsymbol{x}_0 = \alpha\boldsymbol{p}_1 = [\alpha, 0]^\top$: $\boldsymbol{x}(t) = [\alpha, 0]^\top e^{-t}$ — the unstable mode $e^{t}$ is **not present**, even though the system is unstable. Any tiny error in $\boldsymbol{x}_0$ would excite it, though.

### 8.3 Which mode diverges?

The size of mode $i$ is governed by $\lvert e^{\lambda_it}\rvert = e^{\operatorname{Re}(\lambda_i)\,t}$:

| $\operatorname{Re}\lambda_i$ | $\operatorname{Im}\lambda_i$ | Mode behaviour as $t\to\infty$ |
|---|---|---|
| $< 0$ | $0$ | decays monotonically |
| $< 0$ | $\ne0$ | decaying oscillation |
| $= 0$ | $0$ | constant |
| $= 0$ | $\ne0$ | sustained oscillation (bounded) |
| $> 0$ | any | **diverges** ⇒ instability |

**A mode introduces divergence exactly when the real part of its eigenvalue is positive.** If all real parts are negative, every mode decays and $e^{\boldsymbol{A}t}\to\boldsymbol{0}$. (Repeated eigenvalues with $\operatorname{Re}\lambda = 0$ can produce terms like $t$ that diverge slowly — §9.)

---

## 9. Method 3 — Jordan form (non-diagonalizable A)

### 9.1 When eigenvectors are missing

If an eigenvalue is repeated, there may be fewer independent eigenvectors than its multiplicity; then $\boldsymbol{A}$ is **not diagonalizable**. One can still transform to the **Jordan canonical form** $\boldsymbol{J} = \boldsymbol{P}^{-1}\boldsymbol{A}\boldsymbol{P}$ (slide 18): block-diagonal, with a diagonal part $\boldsymbol{J}_0$ and **Jordan blocks**

$$
\boldsymbol{J}_i = \begin{bmatrix}\lambda & 1 & & 0\\ & \lambda & \ddots & \\ & & \ddots & 1\\ 0 & & & \lambda\end{bmatrix} = \lambda\boldsymbol{I} + \boldsymbol{N}_i,
$$

where $\boldsymbol{N}_i$ has ones on the superdiagonal and is **nilpotent** ($\boldsymbol{N}_i^{n_i} = \boldsymbol{0}$).

### 9.2 Exponential of a Jordan block (slide 19)

$\lambda\boldsymbol{I}$ and $\boldsymbol{N}_i$ commute, so $e^{\boldsymbol{J}_it} = e^{\lambda t}e^{\boldsymbol{N}_it}$, and the series for $e^{\boldsymbol{N}_it}$ terminates:

$$
e^{\boldsymbol{J}_it} = e^{\lambda t}\begin{bmatrix}1 & t & \cdots & \frac{t^{n_i-1}}{(n_i-1)!}\\ 0 & 1 & \cdots & \frac{t^{n_i-2}}{(n_i-2)!}\\ \vdots & & \ddots & \vdots\\ 0 & 0 & \cdots & 1\end{bmatrix}.
\qquad\text{For a 2×2 block:}\quad
\exp\!\left(\begin{bmatrix}\lambda & 1\\ 0 & \lambda\end{bmatrix}t\right) = e^{\lambda t}\begin{bmatrix}1 & t\\ 0 & 1\end{bmatrix}.
$$

**Consequence:** non-diagonalizable matrices produce terms $t\,e^{\lambda t}$, $t^2e^{\lambda t}$, … in the response. For $\lambda < 0$ they still decay (the exponential wins over the polynomial), but they can first **rise** — e.g. $t\,e^{-t}$ grows to $e^{-1}$ at $t = 1$ before decaying.

*Example (slide 25):* $\boldsymbol{A} = \begin{bmatrix}0 & 1\\ 0 & 0\end{bmatrix}$ is a Jordan block with $\lambda = 0$: $e^{\boldsymbol{A}t} = \begin{bmatrix}1 & t\\ 0 & 1\end{bmatrix}$, and $\boldsymbol{x}_0 = [\alpha, \beta]^\top$ gives $\boldsymbol{x}(t) = [\alpha + \beta t,\ \beta]^\top$ — a mass moving at constant speed $\beta$. Here the term $t$ diverges even though $\operatorname{Re}\lambda = 0$.

**If $\boldsymbol{A}$ is already in Jordan form** (as in Task 3), no transformation is needed — just exponentiate block by block.

---

## 10. Method 4 — Cayley–Hamilton theorem

### 10.1 The theorem

Every square matrix satisfies its own characteristic equation: if $p(\lambda) = \det(\lambda\boldsymbol{I} - \boldsymbol{A}) = \lambda^n + c_{n-1}\lambda^{n-1} + \dots + c_0$, then

$$
p(\boldsymbol{A}) = \boldsymbol{A}^n + c_{n-1}\boldsymbol{A}^{n-1} + \dots + c_0\boldsymbol{I} = \boldsymbol{0}.
$$

*Example:* $\boldsymbol{A} = \begin{bmatrix}-6 & -1\\ 5 & 0\end{bmatrix}$, $p(\lambda) = \lambda^2 + 6\lambda + 5$. Then $\boldsymbol{A}^2 = \begin{bmatrix}31 & 6\\ -30 & -5\end{bmatrix}$ and $\boldsymbol{A}^2 + 6\boldsymbol{A} + 5\boldsymbol{I} = \begin{bmatrix}31 - 36 + 5 & 6 - 6\\ -30 + 30 & -5 + 0 + 5\end{bmatrix} = \boldsymbol{0}$ ✓.

**Why this helps:** $\boldsymbol{A}^n$ can be written in terms of $\boldsymbol{I},\boldsymbol{A},\dots,\boldsymbol{A}^{n-1}$; so can $\boldsymbol{A}^{n+1} = \boldsymbol{A}\cdot\boldsymbol{A}^n$, and so on. Hence **every power — and every power series, including $e^{\boldsymbol{A}t}$ — collapses to a polynomial of degree $n - 1$** (slide 21):

$$
\exp(\boldsymbol{A}t) = \alpha_0(t)\boldsymbol{I} + \alpha_1(t)\boldsymbol{A} + \dots + \alpha_{n-1}(t)\boldsymbol{A}^{n-1}.
$$

Only the $n$ scalar functions $\alpha_k(t)$ are unknown. For a 2×2 matrix: $e^{\boldsymbol{A}t} = \alpha_0(t)\boldsymbol{I} + \alpha_1(t)\boldsymbol{A}$.

### 10.2 Finding the coefficients: distinct eigenvalues (slides 22–23)

Multiply the polynomial identity by an eigenvector $\boldsymbol{u}_i$: since $\boldsymbol{A}^k\boldsymbol{u}_i = \lambda_i^k\boldsymbol{u}_i$, both sides become multiples of $\boldsymbol{u}_i$, giving a **scalar** equation for every eigenvalue:

$$
e^{\lambda_it} = \sum_{k=0}^{n-1}\alpha_k(t)\lambda_i^k,\qquad i = 1,\dots,n.
$$

In matrix form, with the **Vandermonde matrix** $\boldsymbol{V}$:

$$
\begin{bmatrix}e^{\lambda_1t}\\ \vdots\\ e^{\lambda_nt}\end{bmatrix} = \underbrace{\begin{bmatrix}1 & \lambda_1 & \cdots & \lambda_1^{n-1}\\ \vdots & & & \vdots\\ 1 & \lambda_n & \cdots & \lambda_n^{n-1}\end{bmatrix}}_{\boldsymbol{V}}\begin{bmatrix}\alpha_0(t)\\ \vdots\\ \alpha_{n-1}(t)\end{bmatrix}.
$$

$\boldsymbol{V}$ is invertible **iff all eigenvalues are distinct**; then $[\alpha_0,\dots,\alpha_{n-1}]^\top = \boldsymbol{V}^{-1}[e^{\lambda_1t},\dots,e^{\lambda_nt}]^\top$. **No eigenvectors are needed** — only eigenvalues.

**2×2 with distinct real eigenvalues $\lambda_1\ne\lambda_2$** — solving the two equations $e^{\lambda_it} = \alpha_0 + \alpha_1\lambda_i$:

$$
\alpha_1(t) = \frac{e^{\lambda_1t} - e^{\lambda_2t}}{\lambda_1 - \lambda_2},\qquad
\alpha_0(t) = \frac{\lambda_1e^{\lambda_2t} - \lambda_2e^{\lambda_1t}}{\lambda_1 - \lambda_2}.
$$

(Subtract the two equations to get $\alpha_1$; substitute back for $\alpha_0$.)

### 10.3 Recipe

1. Characteristic polynomial and eigenvalues.
2. Write $e^{\boldsymbol{A}t} = \alpha_0\boldsymbol{I} + \alpha_1\boldsymbol{A} + \dots$ (degree $n - 1$).
3. One equation $e^{\lambda_it} = \sum\alpha_k\lambda_i^k$ per eigenvalue (derivative equations for repeated ones, §10.5).
4. Solve for the $\alpha_k(t)$.
5. Insert into the matrix polynomial and multiply out.
6. **Check** $t = 0$ ⇒ $\boldsymbol{I}$.

### 10.4 Complex eigenvalues (slide 24)

$\boldsymbol{A} = \begin{bmatrix}0 & 1\\ -1 & 0\end{bmatrix}$, $\lambda^2 + 1 = 0$, $\lambda_{1,2} = \pm j$:

$$
e^{jt} = \cos t + j\sin t = \alpha_0 + j\alpha_1,\qquad e^{-jt} = \cos t - j\sin t = \alpha_0 - j\alpha_1
\ \Longrightarrow\ \alpha_0 = \cos t,\ \alpha_1 = \sin t,
$$

(add / subtract the two equations), so $\boldsymbol{\Phi}(t) = \cos t\,\boldsymbol{I} + \sin t\,\boldsymbol{A} = \begin{bmatrix}\cos t & \sin t\\ -\sin t & \cos t\end{bmatrix}$. For general $\sigma\pm j\omega$ the same steps give the shortcut of §7.4.

### 10.5 Repeated eigenvalues (slides 26–28)

If $\lambda_i$ has algebraic multiplicity $\mu_i\ge2$, the Vandermonde matrix has identical rows and is singular. **Fix:** differentiate the scalar identity $e^{\lambda t} = \sum\alpha_k\lambda^k$ with respect to $\lambda$ (the $\alpha_k$ do not depend on $\lambda$) to get $\mu_i - 1$ additional equations:

$$
t\,e^{\lambda_it} = \alpha_1 + 2\alpha_2\lambda_i + \dots + (n-1)\alpha_{n-1}\lambda_i^{n-2},\qquad
t^2e^{\lambda_it} = 2\alpha_2 + \dots,\quad\dots
$$

These replace the duplicate rows.

*Lecture example:* $\boldsymbol{A} = \begin{bmatrix}0 & 0 & 2\\ 0 & 1 & 0\\ -1 & 0 & 3\end{bmatrix}$, $\lambda_1 = 2$, $\lambda_{2,3} = 1$. Ansatz $e^{\lambda t} = \alpha_0 + \alpha_1\lambda + \alpha_2\lambda^2$:

$$
\begin{aligned}
\lambda = 1:&\quad e^{t} = \alpha_0 + \alpha_1 + \alpha_2\\
\tfrac{\mathrm{d}}{\mathrm{d}\lambda}\text{ at }\lambda = 1:&\quad te^{t} = \alpha_1 + 2\alpha_2\\
\lambda = 2:&\quad e^{2t} = \alpha_0 + 2\alpha_1 + 4\alpha_2
\end{aligned}
\qquad\Longrightarrow\qquad
\begin{aligned}
\alpha_0 &= -2te^{t} + e^{2t}\\
\alpha_1 &= (2 + 3t)e^{t} - 2e^{2t}\\
\alpha_2 &= -(t + 1)e^{t} + e^{2t}.
\end{aligned}
$$

Then $\boldsymbol{\Phi}(t) = \alpha_0\boldsymbol{I} + \alpha_1\boldsymbol{A} + \alpha_2\boldsymbol{A}^2$ with $\boldsymbol{A}^2 = \begin{bmatrix}-2 & 0 & 6\\ 0 & 1 & 0\\ -3 & 0 & 7\end{bmatrix}$:

$$
\boldsymbol{\Phi}(t) = \begin{bmatrix}2e^{t} - e^{2t} & 0 & 2e^{2t} - 2e^{t}\\ 0 & e^{t} & 0\\ e^{t} - e^{2t} & 0 & 2e^{2t} - e^{t}\end{bmatrix}.
$$

**Note:** the $\alpha_k$ on slide 28 are correct, but the matrix printed there contains errors in most entries (e.g. terms like $-12t$, which cannot appear because the $t$-terms cancel). The matrix above is the correct one: it satisfies $\boldsymbol{\Phi}(0) = \boldsymbol{I}$ and $\dot{\boldsymbol{\Phi}} = \boldsymbol{A}\boldsymbol{\Phi}$, and agrees with a direct computation (`python/lecture_examples.py`). That the $t\,e^{t}$ terms cancel means this $\boldsymbol{A}$ is actually diagonalizable (the repeated eigenvalue 1 has two independent eigenvectors).

---

## 11. Which method to choose?

| Method | Main idea | When is it useful? |
|---|---|---|
| Series expansion | evaluate $\sum t^k\boldsymbol{A}^k/k!$ | nilpotent $\boldsymbol{A}$ ($\boldsymbol{A}^k = \boldsymbol{0}$) — the series terminates exactly; otherwise a small-$t$ approximation |
| Diagonalization | transform to $\boldsymbol{\Lambda}$, exponentiate the diagonal | $n$ independent, easy eigenvectors; **also gives the modes** (Task 1) |
| Jordan form | exponentiate each Jordan block via its nilpotent part | $\boldsymbol{A}$ not diagonalizable, or already in Jordan form (Task 3) |
| Cayley–Hamilton | $e^{\boldsymbol{A}t}$ as a polynomial of degree $n - 1$ | small matrices, eigenvalues known but eigenvectors inconvenient (Task 4); repeated eigenvalues via derivatives |

Whatever method you use, always run the two **checks**: $\boldsymbol{\Phi}(0) = \boldsymbol{I}$ and $\dot{\boldsymbol{\Phi}}(0) = \boldsymbol{A}$.

---

## 12. Response with an input: the convolution integral

### 12.1 Variation of parameters (slides 31–32)

**Homogeneous part.** $\boldsymbol{x}_h(t) = \exp[\boldsymbol{A}(t - t_0)]\boldsymbol{x}_0$ solves $\dot{\boldsymbol{x}} = \boldsymbol{A}\boldsymbol{x}$ (§5.1).

**Particular part.** Guess a solution of the form "homogeneous solution with a varying coefficient":

$$
\boldsymbol{x}_p(t) = \exp(\boldsymbol{A}t)\,\boldsymbol{p}(t),\qquad \boldsymbol{p}(t)\ \text{unknown}.
$$

Insert into $\dot{\boldsymbol{x}} = \boldsymbol{A}\boldsymbol{x} + \boldsymbol{r}$ (product rule on the left):

$$
\boldsymbol{A}\exp(\boldsymbol{A}t)\boldsymbol{p}(t) + \exp(\boldsymbol{A}t)\dot{\boldsymbol{p}}(t) = \boldsymbol{A}\exp(\boldsymbol{A}t)\boldsymbol{p}(t) + \boldsymbol{r}(t)
\ \Longrightarrow\ \dot{\boldsymbol{p}}(t) = \exp(-\boldsymbol{A}t)\,\boldsymbol{r}(t).
$$

Integrate with $\boldsymbol{p}(t_0) = \boldsymbol{0}$ (the initial condition is already in $\boldsymbol{x}_h$):

$$
\boldsymbol{p}(t) = \int_{t_0}^{t}\exp(-\boldsymbol{A}\tau)\boldsymbol{r}(\tau)\,\mathrm{d}\tau
\ \Longrightarrow\
\boldsymbol{x}_p(t) = \int_{t_0}^{t}\exp\big(\boldsymbol{A}(t - \tau)\big)\boldsymbol{r}(\tau)\,\mathrm{d}\tau = \int_{t_0}^{t}\boldsymbol{\Phi}(t, \tau)\boldsymbol{r}(\tau)\,\mathrm{d}\tau .
$$

This is the matrix version of the integrating-factor result in §1.2.

### 12.2 Total response (slide 33)

With $\boldsymbol{r}(t) = \boldsymbol{B}\boldsymbol{u}(t)$ and the output equation $\boldsymbol{y} = \boldsymbol{C}\boldsymbol{x} + \boldsymbol{D}\boldsymbol{u}$:

$$
\boxed{\;\boldsymbol{x}(t, t_0, \boldsymbol{x}_0) = \boldsymbol{\Phi}(t, t_0)\boldsymbol{x}_0 + \int_{t_0}^{t}\boldsymbol{\Phi}(t, \tau)\boldsymbol{B}\boldsymbol{u}(\tau)\,\mathrm{d}\tau\;}
$$

$$
\boxed{\;\boldsymbol{y}(t) = \underbrace{\boldsymbol{C}\boldsymbol{\Phi}(t, t_0)\boldsymbol{x}_0}_{\text{zero-input / natural response}} + \underbrace{\boldsymbol{C}\int_{t_0}^{t}\boldsymbol{\Phi}(t, \tau)\boldsymbol{B}\boldsymbol{u}(\tau)\,\mathrm{d}\tau + \boldsymbol{D}\boldsymbol{u}(t)}_{\text{zero-state / forced response}}\;}
$$

with $\boldsymbol{\Phi}(t, \tau) = e^{\boldsymbol{A}(t - \tau)}$ for LTI systems. $\boldsymbol{B}$ represents the **actuator structure** — how the input can act on the states.

### 12.3 Reading the convolution integral

Exactly as in the scalar case: the input at time $\tau$ gives a kick $\boldsymbol{B}\boldsymbol{u}(\tau)\,\mathrm{d}\tau$ to the state; by time $t$ that kick has been propagated by $e^{\boldsymbol{A}(t - \tau)}$; the integral sums all past kicks. Only inputs **up to time $t$** enter — the response is **causal** (Week 2).

### 12.4 Change of variables (slides 34–35)

Substituting $\theta = t - \tau$ ($\tau = t - \theta$, $\mathrm{d}\tau = -\mathrm{d}\theta$, limits $\tau = 0\mapsto\theta = t$, $\tau = t\mapsto\theta = 0$) gives the equivalent form

$$
\boldsymbol{x}_p(t) = \int_0^t e^{\boldsymbol{A}(t-\tau)}\boldsymbol{B}\boldsymbol{u}(\tau)\,\mathrm{d}\tau = \int_0^t e^{\boldsymbol{A}\theta}\boldsymbol{B}\boldsymbol{u}(t - \theta)\,\mathrm{d}\theta .
$$

$\theta$ is "how long ago" — the second form weights the input from $\theta$ seconds ago by $e^{\boldsymbol{A}\theta}$. Use whichever form makes the integral easier.

---

## 13. Standard test signals, impulse and step response

### 13.1 Test signals (slide 36)

To compare systems, engineers apply standard inputs: the **unit step** $\sigma(t)$ (switch something on and keep it on), the **rectangle** $\Pi(t)$ and the **Dirac impulse** $\delta(t) = \lim_{\epsilon\to0}\Pi(t)$ (a short hit). See §2.4.

### 13.2 Impulse response (slides 37–38)

Consider only the zero-state response and write the feedthrough also as an integral, $\boldsymbol{D}\boldsymbol{u}(t) = \int\boldsymbol{D}\boldsymbol{u}(\tau)\delta(t - \tau)\,\mathrm{d}\tau$ (sifting property):

$$
\boldsymbol{y}(t) = \int_{t_0}^{t}\big[\boldsymbol{C}\boldsymbol{\Phi}(t, \tau)\boldsymbol{B} + \boldsymbol{D}\delta(t - \tau)\big]\boldsymbol{u}(\tau)\,\mathrm{d}\tau = \int_{t_0}^{t}\boldsymbol{H}(t, \tau)\boldsymbol{u}(\tau)\,\mathrm{d}\tau .
$$

**Definition 1.3 — impulse response matrix:**

$$
\boldsymbol{H}(t, \tau) = \begin{cases}\boldsymbol{C}\boldsymbol{\Phi}(t, \tau)\boldsymbol{B} + \boldsymbol{D}\delta(t - \tau), & t\ge\tau\\ \boldsymbol{0}, & t < \tau\end{cases}
$$

For LTI systems it depends only on $t - \tau$, so with $\tau = 0$:

$$
\boxed{\boldsymbol{H}(t) = \begin{cases}\boldsymbol{C}e^{\boldsymbol{A}t}\boldsymbol{B} + \boldsymbol{D}\delta(t), & t\ge0\\ \boldsymbol{0}, & t < 0\end{cases}}
\qquad\text{and}\qquad
\boldsymbol{y}(t) = \boldsymbol{C}e^{\boldsymbol{A}(t - t_0)}\boldsymbol{x}_0 + \boldsymbol{H}(t) * \boldsymbol{u}(t).
$$

**Why it is called "impulse response":** apply $u(t) = \delta(t)$ with $\boldsymbol{x}(0) = \boldsymbol{0}$. The sifting property collapses the integral to the value of the integrand at $\tau = 0$:

$$
\boldsymbol{x}(t) = \int_{0^-}^{t}e^{\boldsymbol{A}(t-\tau)}\boldsymbol{B}\,\delta(\tau)\,\mathrm{d}\tau = e^{\boldsymbol{A}t}\boldsymbol{B},\qquad
\boldsymbol{y}(t) = \boldsymbol{C}e^{\boldsymbol{A}t}\boldsymbol{B} + \boldsymbol{D}\delta(t) = \boldsymbol{H}(t)\qquad(t\ge0).
$$

**Intuition:** the impulse instantly sets the state to $\boldsymbol{B}$ ("the hammer blow puts the state there"), after which the system evolves freely: $\boldsymbol{x}(t) = e^{\boldsymbol{A}t}\boldsymbol{B}$. So **an impulse response is a natural response started from $\boldsymbol{x}(0^+) = \boldsymbol{B}$**. The entry $h_{ij}(t)$ is a dynamic input-to-output sensitivity: how much output $i$ depends on what input $j$ was $t$ seconds ago (slide 41). Knowing $\boldsymbol{H}(t)$ gives the response to *any* input via the convolution — it is the system's fingerprint.

### 13.3 Step response — and steps that start at $t_0$ (needed for Task 4)

For a constant input $\boldsymbol{u}(\tau) = \boldsymbol{u}_c$ on $[t_0, t]$ and invertible $\boldsymbol{A}$:

$$
\int_{t_0}^{t}e^{\boldsymbol{A}(t - \tau)}\,\mathrm{d}\tau = \boldsymbol{A}^{-1}\big(e^{\boldsymbol{A}(t - t_0)} - \boldsymbol{I}\big)
\qquad\Longrightarrow\qquad
\boldsymbol{x}_p(t) = \boldsymbol{A}^{-1}\big(e^{\boldsymbol{A}(t - t_0)} - \boldsymbol{I}\big)\boldsymbol{B}\boldsymbol{u}_c .
$$

*Derivation:* with $\theta = t - \tau$ the integral is $\int_0^{t - t_0}e^{\boldsymbol{A}\theta}\,\mathrm{d}\theta$; since $\frac{\mathrm{d}}{\mathrm{d}\theta}\big(\boldsymbol{A}^{-1}e^{\boldsymbol{A}\theta}\big) = e^{\boldsymbol{A}\theta}$, the antiderivative is $\boldsymbol{A}^{-1}e^{\boldsymbol{A}\theta}$, evaluated between $0$ and $t - t_0$. (Matrix version of §1.3; $\boldsymbol{A}^{-1}$ and $e^{\boldsymbol{A}\theta}$ commute.)

**Steady state** (all eigenvalues with $\operatorname{Re}\lambda < 0$, so $e^{\boldsymbol{A}(t-t_0)}\to\boldsymbol{0}$):

$$
\boldsymbol{x}(\infty) = -\boldsymbol{A}^{-1}\boldsymbol{B}\boldsymbol{u}_c,\qquad \boldsymbol{y}(\infty) = \big(-\boldsymbol{C}\boldsymbol{A}^{-1}\boldsymbol{B} + \boldsymbol{D}\big)\boldsymbol{u}_c .
$$

Sanity check: at steady state $\dot{\boldsymbol{x}} = \boldsymbol{0} = \boldsymbol{A}\boldsymbol{x}(\infty) + \boldsymbol{B}\boldsymbol{u}_c$ gives the same $\boldsymbol{x}(\infty)$ ✓.

**Delayed step $u(t) = \sigma(t - t_0)$ with the state given at $t_0$.** Split the time axis:

- **$t\le t_0$:** the input is zero, so only the natural motion acts. Since the state is known at $t_0$, use the state transition matrix **backwards** (§4.2): $\boldsymbol{x}(t) = \boldsymbol{\Phi}(t, t_0)\boldsymbol{x}(t_0) = e^{\boldsymbol{A}(t - t_0)}\boldsymbol{x}(t_0)$ with $t - t_0\le0$.
- **$t\ge t_0$:** natural response from $\boldsymbol{x}(t_0)$ plus forced response of a constant input switched on at $t_0$:
  $\boldsymbol{x}(t) = e^{\boldsymbol{A}(t - t_0)}\boldsymbol{x}(t_0) + \boldsymbol{A}^{-1}\big(e^{\boldsymbol{A}(t - t_0)} - \boldsymbol{I}\big)\boldsymbol{B}$.
- Both expressions give $\boldsymbol{x}(t_0)$ at $t = t_0$ — the state (and, with $\boldsymbol{D} = \boldsymbol{0}$, the output) is **continuous**.

**Time invariance pays off:** the result depends only on $t - t_0$. The response to a step at $t_0$ is the response to a step at $0$, shifted by $t_0$ (Week 2, Definition 1.2).

---

## 14. Worked lecture example: total response

Slides 39–41: $\boldsymbol{A} = \begin{bmatrix}0 & 1\\ 0 & 0\end{bmatrix}$, $\boldsymbol{B} = \begin{bmatrix}0\\ 1\end{bmatrix}$, $\boldsymbol{C} = \begin{bmatrix}0 & 1\end{bmatrix}$, $\boldsymbol{D} = 0$, $\boldsymbol{x}_0 = [1, -1]^\top$, $t_0 = 0$, $u(t) = \sigma(t)$. (A mass: $x_1$ position, $x_2$ velocity, input = force, output = velocity.)

1. **State transition matrix.** $\boldsymbol{A}$ is nilpotent, so $\boldsymbol{\Phi}(t) = \boldsymbol{I} + t\boldsymbol{A} = \begin{bmatrix}1 & t\\ 0 & 1\end{bmatrix}$.
2. **Natural response.** $y_h(t) = \boldsymbol{C}\boldsymbol{\Phi}(t)\boldsymbol{x}_0 = \begin{bmatrix}0 & 1\end{bmatrix}\begin{bmatrix}1 - t\\ -1\end{bmatrix} = -1$ — the velocity stays $-1$ without a force.
3. **Forced response.** $y_p(t) = \boldsymbol{C}\int_0^t\boldsymbol{\Phi}(\tau)\boldsymbol{B}\,\mathrm{d}\tau = \begin{bmatrix}0 & 1\end{bmatrix}\begin{bmatrix}\int_0^t\tau\,\mathrm{d}\tau\\ \int_0^t1\,\mathrm{d}\tau\end{bmatrix} = \begin{bmatrix}0 & 1\end{bmatrix}\begin{bmatrix}t^2/2\\ t\end{bmatrix} = t$ — a unit force accelerates the mass at 1 m/s². (The slide uses the $\theta$-form of §12.4.)
4. **Total response.** $y(t) = y_h + y_p = t - 1$.
5. **Impulse response.** $H(t) = \boldsymbol{C}\boldsymbol{\Phi}(t)\boldsymbol{B} = \begin{bmatrix}0 & 1\end{bmatrix}\begin{bmatrix}t\\ 1\end{bmatrix} = 1$: a hammer blow gives the mass a constant velocity forever — the sensitivity is constant in time (special case).

---

## 15. Two "inverse" problems from the exercise

Usually we are given $\boldsymbol{A}$ and $\boldsymbol{x}_0$ and compute the response. Tasks 2 and 3 go the other way.

### 15.1 Given trajectories → find $\boldsymbol{A}$ and $\boldsymbol{\Phi}$ (Task 2)

1. Read off the initial states $\boldsymbol{x}^{(i)}(0)$ by setting $t = 0$ in each trajectory (check they match the given ones).
2. $\boldsymbol{X}(t) = [\boldsymbol{x}^{(1)}(t)\ \cdots\ \boldsymbol{x}^{(n)}(t)]$, $\boldsymbol{X}_0 = \boldsymbol{X}(0)$.
3. $\boldsymbol{\Phi}(t, 0) = \boldsymbol{X}(t)\boldsymbol{X}_0^{-1}$ (§4.3).
4. $\boldsymbol{A} = \dot{\boldsymbol{\Phi}}(0)$.
5. Checks: $\boldsymbol{\Phi}(0) = \boldsymbol{I}$, $\dot{\boldsymbol{\Phi}} = \boldsymbol{A}\boldsymbol{\Phi}$; each given trajectory satisfies $\dot{\boldsymbol{x}} = \boldsymbol{A}\boldsymbol{x}$.

**Bonus insight:** the exponentials appearing in the trajectories are the modes $e^{\lambda_it}$, so the eigenvalues can be read off directly. A trajectory containing only **one** exponential started on an eigenvector (§8.2).

### 15.2 Given a desired natural response → find $\boldsymbol{x}_0$ (Task 3)

1. Compute $e^{\boldsymbol{A}t}$ (any method).
2. Write $\boldsymbol{x}_0 = [a, b, c,\dots]^\top$ with unknowns and compute $y(t) = \boldsymbol{C}e^{\boldsymbol{A}t}\boldsymbol{x}_0$ — a combination of $e^{\lambda_it}$, $te^{\lambda_it}$, … with coefficients linear in $a, b, c$.
3. **Match coefficients** with the desired $y(t)$ (§2.3): one linear equation per function.
4. Solve; check by substitution.

**What is possible?** Only functions that appear as modes of $\boldsymbol{A}$ (and are seen through $\boldsymbol{C}$) can be produced. $y(t) = te^{-t}$ is possible only because $\boldsymbol{A}$ has a 2×2 Jordan block for $\lambda = -1$; a system with only simple eigenvalues could never produce a $t\,e^{\lambda t}$ term. And to keep the unstable mode out of $y(t)$, the initial state must **not excite it**.

---

## 16. Cheat sheet

**State transition matrix**

$$
\dot{\boldsymbol{\Phi}} = \boldsymbol{A}\boldsymbol{\Phi},\ \ \boldsymbol{\Phi}(t_0, t_0) = \boldsymbol{I},\qquad
\boldsymbol{\Phi}(t, t_0) = e^{\boldsymbol{A}(t - t_0)}\ (\text{LTI}),\qquad
\boldsymbol{\Phi}(t, \tau) = \boldsymbol{\Phi}(t, \sigma)\boldsymbol{\Phi}(\sigma, \tau),\qquad
\boldsymbol{\Phi}(t, t_0)^{-1} = \boldsymbol{\Phi}(t_0, t).
$$

**Computing $e^{\boldsymbol{A}t}$**

| Situation | Formula |
|---|---|
| nilpotent ($\boldsymbol{A}^2 = \boldsymbol{0}$) | $\boldsymbol{I} + \boldsymbol{A}t$ |
| diagonalizable | $\boldsymbol{P}\operatorname{diag}(e^{\lambda_it})\boldsymbol{P}^{-1} = \sum_i\boldsymbol{p}_i\tilde{\boldsymbol{p}}_ie^{\lambda_it}$ |
| 2×2, $\lambda = \sigma\pm j\omega$ | $e^{\sigma t}\big[\cos\omega t\,\boldsymbol{I} + \frac{\sin\omega t}{\omega}(\boldsymbol{A} - \sigma\boldsymbol{I})\big]$ |
| Jordan block $\begin{bmatrix}\lambda & 1\\ 0 & \lambda\end{bmatrix}$ | $e^{\lambda t}\begin{bmatrix}1 & t\\ 0 & 1\end{bmatrix}$ |
| Cayley–Hamilton, 2×2, $\lambda_1\ne\lambda_2$ | $\alpha_0\boldsymbol{I} + \alpha_1\boldsymbol{A}$, $\alpha_1 = \frac{e^{\lambda_1t} - e^{\lambda_2t}}{\lambda_1 - \lambda_2}$, $\alpha_0 = \frac{\lambda_1e^{\lambda_2t} - \lambda_2e^{\lambda_1t}}{\lambda_1 - \lambda_2}$ |
| repeated eigenvalue (CH) | add $\frac{\mathrm{d}}{\mathrm{d}\lambda}$: $te^{\lambda t} = \alpha_1 + 2\alpha_2\lambda + \dots$ |

Always check $e^{\boldsymbol{A}\cdot0} = \boldsymbol{I}$ and $\frac{\mathrm{d}}{\mathrm{d}t}e^{\boldsymbol{A}t}\big|_0 = \boldsymbol{A}$.

**Response**

$$
\boldsymbol{y}(t) = \boldsymbol{C}e^{\boldsymbol{A}(t - t_0)}\boldsymbol{x}_0 + \boldsymbol{C}\int_{t_0}^{t}e^{\boldsymbol{A}(t - \tau)}\boldsymbol{B}\boldsymbol{u}(\tau)\,\mathrm{d}\tau + \boldsymbol{D}\boldsymbol{u}(t)
$$

- Impulse ($\boldsymbol{x}_0 = \boldsymbol{0}$): $\boldsymbol{x}(t) = e^{\boldsymbol{A}t}\boldsymbol{B}$, $\boldsymbol{H}(t) = \boldsymbol{C}e^{\boldsymbol{A}t}\boldsymbol{B} + \boldsymbol{D}\delta(t)$.
- Step from $t_0$: $\boldsymbol{x}_p = \boldsymbol{A}^{-1}(e^{\boldsymbol{A}(t - t_0)} - \boldsymbol{I})\boldsymbol{B}$; steady state $\boldsymbol{y}(\infty) = -\boldsymbol{C}\boldsymbol{A}^{-1}\boldsymbol{B} + \boldsymbol{D}$.
- State known at $t_0$, earlier times: $\boldsymbol{x}(t) = e^{\boldsymbol{A}(t - t_0)}\boldsymbol{x}(t_0)$, $t < t_0$.

**Modes:** mode $i$ diverges ⇔ $\operatorname{Re}\lambda_i > 0$. $\boldsymbol{x}_0 = \alpha\boldsymbol{p}_i$ excites only mode $i$. Trajectories from data: $\boldsymbol{\Phi} = \boldsymbol{X}(t)\boldsymbol{X}_0^{-1}$, $\boldsymbol{A} = \dot{\boldsymbol{\Phi}}(0)$.
