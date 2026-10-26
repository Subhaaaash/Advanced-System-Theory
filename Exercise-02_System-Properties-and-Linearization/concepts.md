---
title: "Exercise 2 — Concepts"
subtitle: "Advanced System Theory · Week 02: System Properties and Linearization"
---

# Exercise 2 — Concepts

**Course:** Advanced System Theory (Paderborn University, Automatic Control Group RAT, Prof. Dr. Erdal Kayacan), winter semester 2026/27
**Lecture:** Week 02 — System properties and linearization of nonlinear systems
**Covers:** everything needed for Exercise 2 — Task 1 (system properties), Task 2 (equilibria and linearization), Task 3 (evolution of a discrete-time system, internal and BIBO stability)

This sheet builds every idea from the ground up and assumes only school calculus and a little linear algebra. Each test is given as a recipe you can apply mechanically, followed by worked examples. The full solutions are in [`solutions.md`](solutions.md); the scripts in [`python/`](python/) check every number.

**Notation (as in the lecture slides).**

- As on slide 5: **bold symbols are non-scalar quantities** (vectors and matrices): state $\boldsymbol{x}\in\mathbb{R}^n$, input $\boldsymbol{u}\in\mathbb{R}^m$, output $\boldsymbol{y}$, vector field $\boldsymbol{f}$, output map $\boldsymbol{g}$, system matrices $\boldsymbol{A},\boldsymbol{B},\boldsymbol{C},\boldsymbol{D}$. Plain symbols are scalars: components $x_1, x_2$, scalar signals $u(t), y(t)$ (as in Task 1), numbers $\alpha_1,\alpha_2,\beta_1,\beta_2\in\mathbb{R}$.
- $y(\boldsymbol{x}_0, \boldsymbol{u}, t)$ is the output at time $t$ for initial state $\boldsymbol{x}(t_0) = \boldsymbol{x}_0$ and input $\boldsymbol{u}$ (Definitions 1.1, 1.2). $u_1(t)\to y_1(t)$ reads "input $u_1$ produces output $y_1$" (slide 12).
- $\frac{\mathrm{d}}{\mathrm{d}t}\boldsymbol{x} = \dot{\boldsymbol{x}}$. Operating/equilibrium point $(\boldsymbol{x}_0,\boldsymbol{u}_0)$ with deviations $\Delta\boldsymbol{x} = \boldsymbol{x} - \boldsymbol{x}_0$, $\Delta\boldsymbol{u} = \boldsymbol{u} - \boldsymbol{u}_0$, $\Delta\boldsymbol{y} = \boldsymbol{y} - \boldsymbol{y}_0$; $\cdot|_0$ means evaluation at $(\boldsymbol{x}_0,\boldsymbol{u}_0)$ (slide 25). Task 2(c) of the sheet writes the operating point as $(\boldsymbol{x}^\ast, u^\ast)$.
- Continuous time $t\in\mathbb{R}$ with round brackets $x(t)$; discrete time $k, n\in\mathbb{Z}$ with square brackets $\boldsymbol{x}[k]$, $y[n]$.

---

## Contents

1. [Signals, systems and how to manipulate signals](#1-signals-systems-and-how-to-manipulate-signals)
2. [Linearity](#2-linearity)
3. [Time-invariance](#3-time-invariance)
4. [Causality](#4-causality)
5. [Invertibility and stability](#5-invertibility-and-stability)
6. [Worked lecture exercise: the property table](#6-worked-lecture-exercise-the-property-table)
7. [Nonlinear state-space models](#7-nonlinear-state-space-models)
8. [Equilibrium points and how to find them](#8-equilibrium-points-and-how-to-find-them)
9. [Linearization with Jacobians](#9-linearization-with-jacobians)
10. [What the linear model says about stability (2×2 toolkit)](#10-what-the-linear-model-says-about-stability-22-toolkit)
11. [Discrete-time LTI systems: solution of the state equation](#11-discrete-time-lti-systems-solution-of-the-state-equation)
12. [Computing $\boldsymbol{A}^k$ with eigenvalues and eigenvectors](#12-computing-ak-with-eigenvalues-and-eigenvectors)
13. [Modes, internal stability and BIBO stability](#13-modes-internal-stability-and-bibo-stability)
14. [Cheat sheet](#14-cheat-sheet)

---

## 1. Signals, systems and how to manipulate signals

### 1.1 What a "system" is here

A system maps an input signal to an output signal. With the arrow notation of slide 12,

$$
u(t)\ \to\ y(t),
$$

and, when the system has a state, the output also depends on the initial state: $y(\boldsymbol{x}_0, \boldsymbol{u}, t)$ (slide 5). The property definitions compare **whole signals** — the values at all times $t$ — not single values.

In Task 1 the systems are given as explicit **I/O relationships** $y(t) = \dots$ with scalar $u$ and $y$. They have no initial state, so the argument $\boldsymbol{x}_0$ is simply dropped. Some are **memoryless (static)**: $y(t)$ depends only on $u$ at the same instant $t$ (e.g. $y(t) = u(t)\cos(2\pi f_0 t + \phi)$). Others have **memory**: $y(t)$ depends on $u$ at other times (e.g. $u(t+1)$ or $u(2t)$).

In Tasks 2 and 3 the systems are given as **state-space models** — differential or difference equations with a state $\boldsymbol{x}$ (§7, §11).

### 1.2 Time shift, time scaling, and the "replace $t$" rule

Three operations appear constantly. Let $u(t)$ be a signal and $T > 0$.

| Operation | New signal | Picture |
|---|---|---|
| Delay by $T$ | $u(t - T)$ | graph moves **right** by $T$ (events happen later) |
| Advance by $T$ | $u(t + T)$ | graph moves **left** by $T$ (events happen earlier) |
| Compression by 2 | $u(2t)$ | graph squeezed towards $t=0$; runs twice as fast |

*Why does $u(t-T)$ move right?* The value that $u$ had at time $0$ is now found where $t - T = 0$, i.e. at $t = T$.

**The "replace $t$" rule (the single most important trick for Task 1).** If $y(t) = u(2t)$ and you want the output at another time, say $y(t - T)$, you replace *every* $t$ in the formula by $t - T$:

$$
y(t-T) = u\big(2(t - T)\big) = u(2t - 2T).
$$

The common mistake is to write $u(2t - T)$. That is a different signal: it is $u$ first delayed by $T$ and *then* compressed. Keep the two orders strictly apart:

- **Delay the input first, then apply the system:** feed the delayed input $u(t - T)$ into the system. In $y(t) = u(2t)$ the system reads its input at time $2t$, so the output is $u(2t - T)$. In the lecture notation this is $y(u(t-T), t)$.
- **Apply the system first, then delay the output:** $y(t - T) = u(2(t - T)) = u(2t - 2T)$. In the lecture notation this is $y(u(t), t - T)$.

The safest way to avoid errors is to give the delayed input its own name ($\tilde u(t) := u(t - T)$) and only then substitute.

### 1.3 How properties are proved and disproved

All properties in this sheet are "for all" statements: *for all* inputs, *for all* shifts $T > 0$, *for all* $\alpha_1,\alpha_2\in\mathbb{R}$.

- To **prove** a property you must show it in general (with symbols, not numbers).
- To **disprove** it, a **single counterexample** is enough: one particular input, one particular $T$, one particular $\alpha$ for which it fails.

This is why the solutions show general derivations for "Y" and one concrete counterexample for "N".

---

## 2. Linearity

### 2.1 Definition 1.1 — the superposition principle

For a linear system, the following conditions for $t > t_0$ and $\alpha_1,\alpha_2,\beta_1,\beta_2\in\mathbb{R}$ are fulfilled for all (permissible) initial states $\boldsymbol{x}(t_0) = \boldsymbol{x}_0$, input signals $\boldsymbol{u}(t)$ and initial time $t_0 > 0$:

$$
\begin{aligned}
\boldsymbol{y}(\boldsymbol{0},\ \alpha_1\boldsymbol{u}_1 + \alpha_2\boldsymbol{u}_2,\ t) &= \alpha_1\boldsymbol{y}(\boldsymbol{0},\boldsymbol{u}_1,t) + \alpha_2\boldsymbol{y}(\boldsymbol{0},\boldsymbol{u}_2,t),\\
\boldsymbol{y}(\beta_1\boldsymbol{x}_{0,1} + \beta_2\boldsymbol{x}_{0,2},\ \boldsymbol{0},\ t) &= \beta_1\boldsymbol{y}(\boldsymbol{x}_{0,1},\boldsymbol{0},t) + \beta_2\boldsymbol{y}(\boldsymbol{x}_{0,2},\boldsymbol{0},t),\\
\boldsymbol{y}(\boldsymbol{x}_0,\boldsymbol{u},t) &= \boldsymbol{y}(\boldsymbol{0},\boldsymbol{u},t) + \boldsymbol{y}(\boldsymbol{x}_0,\boldsymbol{0},t).
\end{aligned}
$$

Line by line:

1. linear in the **input** (with zero initial state);
2. linear in the **initial state** (with zero input);
3. total response $=$ **zero-state (forced) response** $+$ **zero-input (free) response**.

The superposition principle comprises **additivity** (the response to a sum is the sum of the responses) and **homogeneity** (scaling the input scales the output) with respect to the initial state and the input. You will see the split of line 3 exactly in Task 3(c): $\boldsymbol{x}[k] = \underbrace{\boldsymbol{A}^k\boldsymbol{x}[0]}_{\text{free}} + \underbrace{\sum\boldsymbol{A}^{k-1-i}\boldsymbol{B}u[i]}_{\text{forced}}$.

### 2.2 The form used for I/O relationships (slide 12)

Without an initial state only line 1 remains. The lecture writes it with arrows:

$$
\text{CT:}\quad x_1(t)\to y_1(t)\ \text{ and }\ x_2(t)\to y_2(t)\ \Longrightarrow\ a\,x_1(t) + b\,x_2(t)\ \to\ a\,y_1(t) + b\,y_2(t),
$$

$$
\text{DT:}\quad x_1[n]\to y_1[n]\ \text{ and }\ x_2[n]\to y_2[n]\ \Longrightarrow\ a\,x_1[n] + b\,x_2[n]\ \to\ a\,y_1[n] + b\,y_2[n].
$$

(On slides 11–13 the input is called $x$; in Task 1 it is called $u$.) In Task 1 notation, with $\alpha_1,\alpha_2$ as in Definition 1.1:

$$
u_1(t)\to y_1(t),\ \ u_2(t)\to y_2(t)\quad\Longrightarrow\quad \alpha_1u_1(t) + \alpha_2u_2(t)\ \to\ \alpha_1y_1(t) + \alpha_2y_2(t).
$$

Every LTI state-space model $\dot{\boldsymbol{x}}(t) = \boldsymbol{A}\boldsymbol{x}(t) + \boldsymbol{B}\boldsymbol{u}(t)$, $\boldsymbol{y}(t) = \boldsymbol{C}\boldsymbol{x}(t) + \boldsymbol{D}\boldsymbol{u}(t)$ is linear in the sense of Definition 1.1.

### 2.3 Recipe: testing an I/O relationship for linearity

1. **Quick necessary test — zero in, zero out.** Put $u(t) = 0$ for all $t$. A linear system must give $y(t) = 0$ for all $t$, because homogeneity with $\alpha_1 = 0$ gives $0\cdot u_1\to 0\cdot y_1 = 0$. If $y\not\equiv 0$, the system is **not linear** — done.
2. **General test.** Let $u_1(t)\to y_1(t)$, $u_2(t)\to y_2(t)$ and define $u_3(t) = \alpha_1u_1(t) + \alpha_2u_2(t)$. Compute its output $y_3(t)$ by inserting $u_3$ into the formula wherever $u$ appears (with the *same* time argument). Check whether $y_3(t) = \alpha_1y_1(t) + \alpha_2y_2(t)$.
3. **If they differ for some choice**, give one concrete counterexample.

### 2.4 Typical patterns

| Relationship | Linear? | Reason |
|---|---|---|
| $y(t) = c(t)\,u(t - \tau)$, any known function $c(t)$ | Y | $u$ appears to the first power, multiplied only by things that do not depend on $u$ |
| $y(t) = u(\psi(t))$, any time map $\psi$ (e.g. $u(2t)$, $u(-t)$) | Y | changing the *time argument* is linear: $u_3(\psi(t)) = \alpha_1u_1(\psi(t)) + \alpha_2u_2(\psi(t))$ |
| sums, integrals, derivatives of $u$ | Y | these operations are linear |
| $y = u^2$, $3x(t)x(t-1)$, $\sqrt{u}$, $\lvert u\rvert$ | N | not homogeneous: $(2u)^2 = 4u^2 \neq 2u^2$ |
| $y = \cos(u)$, $\sin(\dots + u)$, $e^{u}$ | N | nonlinear function *of the input value* |
| $y = 3u + 1$ (affine) | N | fails zero-in-zero-out |

**Key distinction (needed for Task 1 c, d):** a function of *time* that multiplies the input, like $u(t)\cos(2\pi f_0 t + \phi)$, keeps the system linear. A nonlinear function *applied to the input value*, like $A\cos(2\pi f_0 t + u(t))$, destroys linearity.

---

## 3. Time-invariance

### 3.1 Definition 1.2

A system is called **time-invariant** if, for every admissible initial state $\boldsymbol{x}_0$, every admissible input $\boldsymbol{u}(t)$, and every time shift $T > 0$,

$$
\boldsymbol{y}\big(\boldsymbol{x}_0,\ \boldsymbol{u}(t),\ t - T\big) = \boldsymbol{y}\big(\boldsymbol{x}_0,\ \boldsymbol{u}(t - T),\ t\big).
$$

- **Left side:** the original output, observed $T$ later — i.e. the output delayed by $T$.
- **Right side:** the output at time $t$ when the input is delayed by $T$.

In words: *the behaviour and characteristics of the system are fixed over time* — the same experiment carried out later gives the same result, just later. Otherwise the system is **time-variant**. In the lecture's solution (slide 12): a time shift of the input signal results in the same time shift of the output signal.

A system is **linear and time-invariant (LTI)** if the system matrices $\boldsymbol{A},\boldsymbol{B},\boldsymbol{C},\boldsymbol{D}$ are constant ($\boldsymbol{A}(t) = \boldsymbol{A}\ \forall\, t > 0$ and so on):

$$
\dot{\boldsymbol{x}}(t) = \boldsymbol{A}\boldsymbol{x}(t) + \boldsymbol{B}\boldsymbol{u}(t),\qquad \boldsymbol{y}(t) = \boldsymbol{C}\boldsymbol{x}(t) + \boldsymbol{D}\boldsymbol{u}(t).
$$

### 3.2 Recipe: compare both sides of Definition 1.2

For an I/O relationship $y(t) = \dots$ (no $\boldsymbol{x}_0$):

1. **Right side $y(u(t - T), t)$ — delay the input, then apply the system.** In the formula, replace every $u(\text{something})$ by $u(\text{something} - T)$. Leave every other $t$ untouched.
2. **Left side $y(u(t), t - T)$ — apply the system, then delay the output.** Take the original formula $y(t)$ and replace **every** $t$ by $t - T$ (the "replace $t$" rule of §1.2).
3. **Compare.** If both sides agree for all $u$, $t$ and $T > 0$, the system is time-invariant; otherwise it is time-variant.

### 3.3 Red flags and patterns

- **$t$ appears explicitly outside the input** (as in $\cos(2\pi f_0 t)$, $t\,u(t)$, $0.5e^tu(t)$): almost always time-variant. On the right side the explicit $t$ stays $t$; on the left side it becomes $t - T$.
  *Exception:* if the explicit time function is constant (e.g. $f_0 = 0$), it disappears and the system can be time-invariant.
- **Time scaling or reversal of the input** ($u(2t)$, $x[4 - n]$): time-variant, because a shift $T$ of the input becomes a shift of $T/2$ (or $-T$) at the output.
- **Pure shifts of the input** ($u(t-1)$, $u(t+1)$): time-invariant.
- **State space:** constant $\boldsymbol{A},\boldsymbol{B},\boldsymbol{C},\boldsymbol{D}$ ⇒ time-invariant; $\boldsymbol{A}(t)$, … depending on time ⇒ generally time-variant.

*Example (time-invariant).* $y(t) = u(t-1)$. Right side: $y(u(t-T),t) = u(t - 1 - T)$. Left side: $y(u(t), t-T) = u(t - T - 1)$. Equal ⇒ time-invariant.

*Example (time-variant).* $y(t) = t\,u(t)$. Right side: $y(u(t-T),t) = t\,u(t - T)$. Left side: $y(u(t),t-T) = (t - T)\,u(t - T)$. These differ (e.g. $u(t) = 1$, $T = 1$: $t$ vs $t - 1$) ⇒ time-variant.

---

## 4. Causality

### 4.1 Definition 1.3

A system is called **causal** if and only if for any input signal pair $u_1(t)$ and $u_2(t)$ and any $T\in\lbrace\mathbb{R}\mid 0 < T\rbrace$ such that

$$
u_1(t) = u_2(t),\quad \forall t\le T,
$$

the corresponding outputs satisfy

$$
y_1(t) = y_2(t),\quad \forall t\le T.
$$

In simple words: **a system is causal if the output depends on past and current inputs but not on future inputs.** Every physical system that runs in real time is causal: it cannot react to something that has not happened yet.

### 4.2 Recipe

For the output at time $t$, list every time instant at which the input is read ("read times").

- If **all** read times are $\le t$ **for all** $t$ ⇒ **causal**.
- If **some** read time is $> t$ for **some** $t$ ⇒ **non-causal**. Counterexample: take $u_1, u_2$ that agree for $t\le T$ but differ just after $T$ (Definition 1.3 with that $T$).

Memoryless systems read only time $t$, so **every memoryless system is causal**.

### 4.3 Patterns

| Read time | Causal? | Comment |
|---|---|---|
| $u(t)$, $x(t - 1)$, $\sum_{k=-\infty}^{n}x[k]$ | Y | present and past only |
| $u(t + 1)$ | N | needs the input one second in the future |
| $u(2t)$ | N | for $t > 0$: $2t > t$, i.e. the future (for $t<0$ it reads the past, but causality must hold for *all* $t$) |
| $x[4 - n]$ | N | for $n = 0$ it reads $x[4]$, the future |

Non-causal does not mean useless: in offline signal processing (a recorded file) the "future" samples are available. It only means the system cannot run in real time.

### 4.4 State-space models are causal

$\frac{\mathrm{d}}{\mathrm{d}t}\boldsymbol{x} = \boldsymbol{f}(\boldsymbol{x},\boldsymbol{u})$ is integrated forwards: $\boldsymbol{x}(t) = \boldsymbol{x}(t_0) + \int_{t_0}^t \boldsymbol{f}(\boldsymbol{x}(\tau), \boldsymbol{u}(\tau))\,\mathrm{d}\tau$ uses inputs only up to $t$. As the lecture states, **the LTI state-space representation is causal by definition.**

---

## 5. Invertibility and stability

These two are not asked for in Task 1, but they are part of the lecture and of the lecture exercise in §6.

### 5.1 Definition 1.4 — invertible system

A system is said to be **invertible** if distinct inputs lead to distinct outputs.

- $y(t) = 2u(t)$ has the inverse system $w(t) = 0.5\,y(t)$ (slide 8).
- $y(t) = u(t)^2$ is not invertible ($u$ and $-u$ give the same output).
- $y(t) = u(t)\cos(2\pi f_0t + \phi)$ is not invertible: at the zeros of the cosine the input value is lost.

### 5.2 Definition 1.5 — stability

A system is said to be **stable** if small perturbations to the system lead to responses that do not diverge. In layman's terms (slide 10): a stable system tends to return to a particular state after being temporarily disturbed (ball in a valley); an unstable one moves away (ball on a hilltop).

Two precise versions are used in the course:

- **BIBO stability (external, week 10)** — used in the lecture's solution (slide 12): every **b**ounded **i**nput produces a **b**ounded **o**utput: $\lvert u(t)\rvert \le M_u < \infty\ \forall t \;\Rightarrow\; \lvert y(t)\rvert \le M_y < \infty\ \forall t$.
- **Internal stability (week 9)** — with zero input, the state returns to (or stays near) the equilibrium after a perturbation of the initial state. For linear systems this is decided by the eigenvalues of $\boldsymbol{A}$ (§10, §13).

*Lecture question (slide 9):* $y(t) = 0.5e^tu(t)$. The bounded input $u(t) = 1$ gives $y(t) = 0.5e^t\to\infty$. Not BIBO stable.

**BIBO recipe for I/O relationships.** Assume $\lvert u\rvert\le M$. Try to bound $\lvert y\rvert$ using $\lvert\cos\rvert\le1$, the triangle inequality $\lvert a + b\rvert\le\lvert a\rvert + \lvert b\rvert$, etc. If a bound exists ⇒ BIBO stable. If you find one bounded input with an unbounded output (sums that grow, factors like $e^t$ or $t$) ⇒ not BIBO stable.

---

## 6. Worked lecture exercise: the property table

Slides 11–13 ask for causality, stability, time-invariance and linearity of four systems (input $x$, output $y$). Working them through is the best preparation for Task 1.

**(1) $y[n] = x[n - 4] + x[4 - n]$**

- *Causal?* Read times $n - 4$ and $4 - n$. For $n = 0$: $4 - n = 4 > 0$ (future). **N.**
- *Stable?* $\lvert y[n]\rvert \le \lvert x[n-4]\rvert + \lvert x[4-n]\rvert \le 2M$. **Y.**
- *Time-invariant?* Input delayed: $x[n - 4 - T] + x[4 - n - T]$. Output delayed: $y[n - T] = x[n - T - 4] + x[4 - n + T]$. The second terms differ ($-T$ vs $+T$) because of the time reversal. **N.**
- *Linear?* Only time arguments are changed. **Y.**

**(2) $y(t) = 3x(t) + x(t - 1)$** — present and past ⇒ causal; $\lvert y\rvert\le 4M$ ⇒ stable; pure shifts with constant coefficients ⇒ time-invariant; first powers only ⇒ linear. **Y Y Y Y.**

**(3) $y(t) = 3x(t)x(t - 1)$** — causal, stable ($\lvert y\rvert\le 3M^2$), time-invariant; but the input $2x$ gives $3(2x(t))(2x(t-1)) = 4\cdot3x(t)x(t-1)\ne 2\,y(t)$ ⇒ **not linear**. **Y Y Y N.**

**(4) $y[n] = \sum_{k=-\infty}^{n} x[k]$ (accumulator)** — past and present ⇒ causal. The bounded input $x[k] = 1$ for $k\ge0$ (0 before) gives $y[n] = n + 1\to\infty$ ⇒ **not stable**. Input delayed: $\sum_{k\le n}x[k - T] = \sum_{j\le n-T}x[j]$; output delayed: $y[n - T] = \sum_{k\le n - T}x[k]$; equal ⇒ time-invariant. A sum is linear ⇒ linear. **Y N Y Y.**

| | causal | stable | time-invariant | linear |
|---|:---:|:---:|:---:|:---:|
| $y[n] = x[n-4] + x[4-n]$ | N | Y | N | Y |
| $y(t) = 3x(t) + x(t-1)$ | Y | Y | Y | Y |
| $y(t) = 3x(t)x(t-1)$ | Y | Y | Y | N |
| $y[n] = \sum_{k=-\infty}^{n}x[k]$ | Y | N | Y | Y |

---

## 7. Nonlinear state-space models

### 7.1 The general form (slide 16)

A nonlinear but time-invariant dynamic system can usually be represented by a set of nonlinear differential equations in the form

$$
\frac{\mathrm{d}}{\mathrm{d}t}\boldsymbol{x} = \boldsymbol{f}(\boldsymbol{x},\boldsymbol{u}),
$$

where $\boldsymbol{f}:\mathbb{R}^n\times\mathbb{R}^m\to\mathbb{R}^n$ is a continuously differentiable nonlinear vector field, $\boldsymbol{x}\in\mathbb{R}^n$ is the state vector and $\boldsymbol{u}\in\mathbb{R}^m$ is the input vector. If there is an output, it is given by a time-invariant output equation $\boldsymbol{y} = \boldsymbol{g}(\boldsymbol{x},\boldsymbol{u})$ (slide 22). Written out for $n = 2$:

$$
\dot x_1(t) = f_1\big(x_1(t), x_2(t), u(t)\big),\qquad \dot x_2(t) = f_2\big(x_1(t), x_2(t), u(t)\big).
$$

A system **without input** ($\frac{\mathrm{d}}{\mathrm{d}t}\boldsymbol{x} = \boldsymbol{f}(\boldsymbol{x})$, as in Task 2a and 2b) is called **autonomous**.

### 7.2 From an $n$-th order ODE to state space

Many models arrive as one high-order differential equation, like the pendulum $ml^2\ddot\theta + b\dot\theta + mgl\sin\theta = 0$ or the I/O relation of Task 2(c). The recipe:

1. **Solve for the highest derivative.** For a second-order equation: $\ddot y = h(y, \dot y, u)$.
2. **Choose the states** as the output and its derivatives up to one order below the highest: $x_1 = y$, $x_2 = \dot y$, i.e. $\boldsymbol{x} = [x_1, x_2]^\top = [y, \dot y]^\top$ (for order $n$: $x_i = y^{(i-1)}$, $i = 1,\dots,n$).
3. **Write one first-order equation per state:**
   - $\dot x_1 = \dot y = x_2$ (this is just the definition of $x_2$),
   - $\dot x_2 = \ddot y = h(x_1, x_2, u)$ (the ODE, with $y\to x_1$, $\dot y\to x_2$).
4. **Output equation:** $y = x_1$, so $g(\boldsymbol{x},u) = x_1$.

*Pendulum example (slide 18).* With $x_1 = \theta$ and $x_2 = \dot\theta$, the state-space equation is

$$
\begin{bmatrix}\dot x_1\\ \dot x_2\end{bmatrix} = \begin{bmatrix} x_2\\ -\frac{b}{ml^2}x_2 - \frac{g}{l}\sin x_1\end{bmatrix}.
$$

---

## 8. Equilibrium points and how to find them

### 8.1 Definition 1.6 — equilibrium

For a nonlinear state-space model

$$
\frac{\mathrm{d}}{\mathrm{d}t}\boldsymbol{x} = \boldsymbol{f}(\boldsymbol{x},\boldsymbol{u}),
$$

an equilibrium point $(\boldsymbol{x}_0,\boldsymbol{u}_0)$ satisfies

$$
\boldsymbol{f}(\boldsymbol{x}_0,\boldsymbol{u}_0) = \boldsymbol{0}.
$$

At an equilibrium point the state remains constant ($\dot{\boldsymbol{x}} = \boldsymbol{0}$): if the state starts at $\boldsymbol{x}_0$ and the input is held at $\boldsymbol{u}_0$, it stays at $\boldsymbol{x}_0$ for all future time. For an autonomous system the condition is $\boldsymbol{f}(\boldsymbol{x}_0) = \boldsymbol{0}$.

Finding equilibria means solving $n$ (usually nonlinear) equations $f_1 = 0,\dots,f_n = 0$ for the $n$ unknowns. There may be none, one, several, infinitely many isolated ones, or a whole continuum.

### 8.2 Strategy

1. **Look for an equation with only one unknown** (or that factorises) and solve it first; substitute into the others.
2. **Use trigonometric facts:**
   - $\sin z = 0 \iff z = k\pi$, $k\in\mathbb{Z}$.
   - $\cos(k\pi) = (-1)^k$ ($+1$ for even $k$, $-1$ for odd $k$); $\sin(k\pi) = 0$.
   - $\lvert\sin z\rvert\le 1$, $\lvert\cos z\rvert\le1$; $\sin$ is odd, $\cos$ is even.
   - $1 - \cos z = 2\sin^2(z/2) \ge 0$, and $1 - \cos z \le z^2/2$.
3. **Always check the obvious candidate $\boldsymbol{x} = \boldsymbol{0}$.**
4. **Combine equations** (add, subtract) to cancel the hard terms.
5. **Bound the solutions** when the equations cannot be solved in closed form: e.g. if $a^2 + b^2 = -b\cos a$, then $a^2 + b^2 \le \lvert b\rvert$, which confines $(a,b)$ to a small disk. Inside a small region, signs of terms become known (e.g. $\cos a > 0$ for $\lvert a\rvert < \pi/2$), and this often rules out all but a few points.
6. **Verify** each candidate by substituting it back into **all** equations.
7. **For a system with input**, the equilibrium condition involves $\boldsymbol{u}_0$ too. Often one asks "for which constant $\boldsymbol{u}_0$ do equilibria exist?" (see Task 2c).

### 8.3 Examples from the lecture

- *Pendulum (slide 18):* the equilibrium points are given by $x_2 = 0$ and $\sin(x_1) = 0$. This leads to an infinite number of equilibrium points $(k\pi, 0)$, $k\in\mathbb{Z}$. Even $k$: hanging down; odd $k$: upright.
- *Cart-pendulum (slide 26):* $x_3 = x_4 = 0$ (velocities zero); with $u_0 = 0$ the third and fourth components of $\boldsymbol{f}$ need $\sin x_2 = 0$; the cart position $x_1$ is free. The lecture chooses $\boldsymbol{x}_0 = \boldsymbol{0}$, $u_0 = 0$.

---

## 9. Linearization with Jacobians

### 9.1 Why do we linearize? (slide 15)

- Most real-world dynamical systems are nonlinear, and nonlinear systems are often difficult to analyse and control directly: no general closed-form solution, behaviour depending strongly on the initial condition, multiple equilibria, saturation.
- A linear model gives access to powerful tools: eigenvalues and modes, local stability analysis, controllability and observability, controller design (LQR, pole placement, observers).

**A nonlinear system is not replaced by a linear system globally; it is approximated locally around an operating point.**

### 9.2 One variable: the tangent line (slide 20)

A scalar function $f:\mathbb{R}\to\mathbb{R}$ is approximated around $x_0$ by the first two terms of the Taylor series:

$$
f(x) \approx f(x_0) + f'(x_0)(x - x_0).
$$

The neglected terms are of order $(x - x_0)^2$: the tangent line is a good approximation only for $x\approx x_0$, and the approximation error grows as we move away from the operating point. Example: $\sin x\approx x$ near $0$ (since $\sin 0 = 0$, $\cos 0 = 1$) — the small-angle approximation of slide 23; $\sin x \approx -(x - \pi)$ near $\pi$ (since $\cos\pi = -1$).

### 9.3 Several variables: partial derivatives and the Jacobian matrices (slide 19)

A **partial derivative** $\frac{\partial f_i}{\partial x_j}$ is the ordinary derivative of $f_i$ with respect to $x_j$, treating all other variables as constants.

Rules you need (with $x_2$ treated as a constant when differentiating w.r.t. $x_1$, and vice versa):

| Expression | $\frac{\partial}{\partial x_1}$ | $\frac{\partial}{\partial x_2}$ |
|---|---|---|
| $x_1^2 + x_2^2$ | $2x_1$ | $2x_2$ |
| $x_2\cos x_1$ | $-x_2\sin x_1$ | $\cos x_1$ |
| $x_1\sin x_2$ | $\sin x_2$ | $x_1\cos x_2$ |
| $(1 + x_1)x_1 = x_1 + x_1^2$ | $1 + 2x_1$ | $0$ |
| $(3 + x_2^2)x_2 = 3x_2 + x_2^3$ | $0$ | $3 + 3x_2^2$ |
| $(1 + x_1 + x_1^2)\,u$ | $(1 + 2x_1)\,u$ | $0$ |

(Product rule: $(ab)' = a'b + ab'$; chain rule: $(\sin(h))' = \cos(h)\,h'$.) Tip: multiply out products before differentiating — fewer mistakes.

We can linearize $\boldsymbol{f}(\boldsymbol{x},\boldsymbol{u})$ around the operating point $(\boldsymbol{x}_0,\boldsymbol{u}_0)$ by expanding it into a Taylor series and cutting off all higher-order terms:

$$
\frac{\mathrm{d}}{\mathrm{d}t}\boldsymbol{x} = \boldsymbol{f}(\boldsymbol{x},\boldsymbol{u}) \approx \boldsymbol{f}(\boldsymbol{x}_0,\boldsymbol{u}_0) + \left.\frac{\partial\boldsymbol{f}}{\partial\boldsymbol{x}}\right|_{(\boldsymbol{x}_0,\boldsymbol{u}_0)}(\boldsymbol{x} - \boldsymbol{x}_0) + \left.\frac{\partial\boldsymbol{f}}{\partial\boldsymbol{u}}\right|_{(\boldsymbol{x}_0,\boldsymbol{u}_0)}(\boldsymbol{u} - \boldsymbol{u}_0)
$$

with the Jacobian matrices (row $i$ belongs to $f_i$, column $j$ to the variable)

$$
\frac{\partial\boldsymbol{f}}{\partial\boldsymbol{x}} = \begin{bmatrix}
\frac{\partial f_1}{\partial x_1} & \cdots & \frac{\partial f_1}{\partial x_n}\\
\vdots & \ddots & \vdots\\
\frac{\partial f_n}{\partial x_1} & \cdots & \frac{\partial f_n}{\partial x_n}
\end{bmatrix}
\qquad\text{and}\qquad
\frac{\partial\boldsymbol{f}}{\partial\boldsymbol{u}} = \begin{bmatrix}
\frac{\partial f_1}{\partial u_1} & \cdots & \frac{\partial f_1}{\partial u_m}\\
\vdots & \ddots & \vdots\\
\frac{\partial f_n}{\partial u_1} & \cdots & \frac{\partial f_n}{\partial u_m}
\end{bmatrix}.
$$

### 9.4 The linearized state-space model (slides 21–22)

Consider the deviations from the equilibrium point $\Delta\boldsymbol{x} = \boldsymbol{x} - \boldsymbol{x}_0$ and $\Delta\boldsymbol{u} = \boldsymbol{u} - \boldsymbol{u}_0$. Since $\boldsymbol{x}_0$ is constant, $\frac{\mathrm{d}}{\mathrm{d}t}\boldsymbol{x}_0 = \boldsymbol{0}$, and

$$
\begin{aligned}
\frac{\mathrm{d}}{\mathrm{d}t}(\Delta\boldsymbol{x}) &= \frac{\mathrm{d}}{\mathrm{d}t}\boldsymbol{x} - \frac{\mathrm{d}}{\mathrm{d}t}\boldsymbol{x}_0\\
&\approx \left[\boldsymbol{f}(\boldsymbol{x}_0,\boldsymbol{u}_0) + \left.\frac{\partial\boldsymbol{f}}{\partial\boldsymbol{x}}\right|_{(\boldsymbol{x}_0,\boldsymbol{u}_0)}(\boldsymbol{x} - \boldsymbol{x}_0) + \left.\frac{\partial\boldsymbol{f}}{\partial\boldsymbol{u}}\right|_{(\boldsymbol{x}_0,\boldsymbol{u}_0)}(\boldsymbol{u} - \boldsymbol{u}_0)\right] - \frac{\mathrm{d}}{\mathrm{d}t}\boldsymbol{x}_0\\
&= \boldsymbol{0} + \left.\frac{\partial\boldsymbol{f}}{\partial\boldsymbol{x}}\right|_{(\boldsymbol{x}_0,\boldsymbol{u}_0)}\Delta\boldsymbol{x} + \left.\frac{\partial\boldsymbol{f}}{\partial\boldsymbol{u}}\right|_{(\boldsymbol{x}_0,\boldsymbol{u}_0)}\Delta\boldsymbol{u} - \boldsymbol{0}\\
&= \boldsymbol{A}\Delta\boldsymbol{x} + \boldsymbol{B}\Delta\boldsymbol{u}.
\end{aligned}
$$

The first $\boldsymbol{0}$ is $\boldsymbol{f}(\boldsymbol{x}_0,\boldsymbol{u}_0) = \boldsymbol{0}$ — this is why the equilibrium is the natural choice: **the resulting deviation model has no constant offset term** (slide 17). At a non-equilibrium operating point the constant $\boldsymbol{f}(\boldsymbol{x}_0,\boldsymbol{u}_0)$ would remain.

If the system also has a time-invariant output equation $\boldsymbol{y} = \boldsymbol{g}(\boldsymbol{x},\boldsymbol{u})$, then with $\Delta\boldsymbol{y} = \boldsymbol{y} - \boldsymbol{y}_0 = \boldsymbol{g}(\boldsymbol{x},\boldsymbol{u}) - \boldsymbol{g}(\boldsymbol{x}_0,\boldsymbol{u}_0)$:

$$
\Delta\boldsymbol{y} \approx \left.\frac{\partial\boldsymbol{g}}{\partial\boldsymbol{x}}\right|_{(\boldsymbol{x}_0,\boldsymbol{u}_0)}\Delta\boldsymbol{x} + \left.\frac{\partial\boldsymbol{g}}{\partial\boldsymbol{u}}\right|_{(\boldsymbol{x}_0,\boldsymbol{u}_0)}\Delta\boldsymbol{u} = \boldsymbol{C}\Delta\boldsymbol{x} + \boldsymbol{D}\Delta\boldsymbol{u},
$$

with the Jacobians $\frac{\partial\boldsymbol{g}}{\partial\boldsymbol{x}}\in\mathbb{R}^{p\times n}$ and $\frac{\partial\boldsymbol{g}}{\partial\boldsymbol{u}}\in\mathbb{R}^{p\times m}$ built like those of $\boldsymbol{f}$. Summary:

$$
\boxed{\;\frac{\mathrm{d}}{\mathrm{d}t}(\Delta\boldsymbol{x}) = \boldsymbol{A}\Delta\boldsymbol{x} + \boldsymbol{B}\Delta\boldsymbol{u},\quad \Delta\boldsymbol{y} = \boldsymbol{C}\Delta\boldsymbol{x} + \boldsymbol{D}\Delta\boldsymbol{u},\quad
\boldsymbol{A} = \left.\frac{\partial\boldsymbol{f}}{\partial\boldsymbol{x}}\right|_0,\ \boldsymbol{B} = \left.\frac{\partial\boldsymbol{f}}{\partial\boldsymbol{u}}\right|_0,\ \boldsymbol{C} = \left.\frac{\partial\boldsymbol{g}}{\partial\boldsymbol{x}}\right|_0,\ \boldsymbol{D} = \left.\frac{\partial\boldsymbol{g}}{\partial\boldsymbol{u}}\right|_0\;}
$$

### 9.5 Recipe

1. Write the model as $\frac{\mathrm{d}}{\mathrm{d}t}\boldsymbol{x} = \boldsymbol{f}(\boldsymbol{x},\boldsymbol{u})$, $\boldsymbol{y} = \boldsymbol{g}(\boldsymbol{x},\boldsymbol{u})$ (§7.2).
2. Find the equilibria $(\boldsymbol{x}_0,\boldsymbol{u}_0)$ from $\boldsymbol{f}(\boldsymbol{x}_0,\boldsymbol{u}_0) = \boldsymbol{0}$ (§8).
3. Compute the Jacobians **symbolically** (as functions of $\boldsymbol{x}$, $\boldsymbol{u}$).
4. **Evaluate** at each equilibrium ($\cdot|_0$) to get the numerical $\boldsymbol{A},\boldsymbol{B},\boldsymbol{C},\boldsymbol{D}$. (Differentiate first, evaluate second — never the other way round.)
5. State the model in deviation variables, and say which equilibrium it belongs to.
6. Optionally interpret: eigenvalues of $\boldsymbol{A}$ (§10).

### 9.6 Lecture examples

*Pendulum at $(k\pi, 0)$.* $\frac{\partial\boldsymbol{f}}{\partial\boldsymbol{x}} = \begin{bmatrix}0 & 1\\ -\frac{g}{l}\cos x_1 & -\frac{b}{ml^2}\end{bmatrix}$, and $\cos(k\pi) = (-1)^k$:

$$
\boldsymbol{A}_{\text{down}} = \begin{bmatrix}0 & 1\\ -\frac{g}{l} & -\frac{b}{ml^2}\end{bmatrix}\ (k\ \text{even}),\qquad
\boldsymbol{A}_{\text{up}} = \begin{bmatrix}0 & 1\\ +\frac{g}{l} & -\frac{b}{ml^2}\end{bmatrix}\ (k\ \text{odd}).
$$

For the undamped pendulum of slide 23 ($b = 0$, written with $\theta$ and $\omega = \dot\theta$: $\frac{\mathrm{d}}{\mathrm{d}t}\theta = \omega$, $\frac{\mathrm{d}}{\mathrm{d}t}\omega = -\frac{g}{l}\sin\theta$), the small-angle approximation $\sin\theta\approx\theta$ gives exactly $\boldsymbol{A}_{\text{down}}$:

$$
\frac{\mathrm{d}}{\mathrm{d}t}\begin{bmatrix}\Delta\theta\\ \Delta\omega\end{bmatrix} = \begin{bmatrix}0 & 1\\ -g/l & 0\end{bmatrix}\begin{bmatrix}\Delta\theta\\ \Delta\omega\end{bmatrix}.
$$

*Cart-pendulum (slides 24–26).* $\boldsymbol{x}(t)^\top = \big[\zeta(t), \phi(t), \dot\zeta(t), \dot\phi(t)\big]$, $u(t) = F(t)$, $y(t) = x_2(t)$. At $x_2 = 0$: $\sin x_2 = 0$, $\cos x_2 = 1$, and every term with $x_4^2$ has zero derivative at $x_4 = 0$ (derivative $2x_4 = 0$). With $\boldsymbol{x}_0 = \boldsymbol{0}$, $u_0 = 0$:

$$
\boldsymbol{A} = \begin{bmatrix}0&0&1&0\\0&0&0&1\\0&\frac{mg}{M}&0&0\\0&-\frac{(M+m)g}{lM}&0&0\end{bmatrix},\quad
\boldsymbol{B} = \begin{bmatrix}0\\0\\\frac1M\\-\frac1{lM}\end{bmatrix},\quad \boldsymbol{C} = \begin{bmatrix}0&1&0&0\end{bmatrix},\quad \boldsymbol{D} = 0.
$$

### 9.7 Validity

The linear approximation is only valid **locally** around the operating point and works well when $\Delta\boldsymbol{x}$ and $\Delta\boldsymbol{u}$ are sufficiently small (slide 21). Each equilibrium has its *own* linear model, and they can behave completely differently (pendulum down: oscillation; up: falls over). The figures in [`python/figures/`](python/figures/) show nonlinear and linearized trajectories lying on top of each other near the equilibrium and separating further away.

---

## 10. What the linear model says about stability (2×2 toolkit)

Task 2 asks only for the linearization, but interpreting it is what linearization is for, and you need eigenvalues for Task 3 anyway.

### 10.1 Eigenvalues of a 2×2 matrix

The eigenvalues $\lambda$ of $\boldsymbol{A} = \begin{bmatrix}a_{11} & a_{12}\\ a_{21} & a_{22}\end{bmatrix}$ solve the **characteristic equation**

$$
\det(\lambda\boldsymbol{I} - \boldsymbol{A}) = (\lambda - a_{11})(\lambda - a_{22}) - a_{12}a_{21} = \lambda^2 - \operatorname{tr}(\boldsymbol{A})\,\lambda + \det(\boldsymbol{A}) = 0,
$$

with $\operatorname{tr}\boldsymbol{A} = a_{11} + a_{22}$ and $\det\boldsymbol{A} = a_{11}a_{22} - a_{12}a_{21}$. So

$$
\lambda_{1,2} = \frac{\operatorname{tr}\boldsymbol{A} \pm\sqrt{(\operatorname{tr}\boldsymbol{A})^2 - 4\det\boldsymbol{A}}}{2},\qquad \lambda_1 + \lambda_2 = \operatorname{tr}\boldsymbol{A},\quad \lambda_1\lambda_2 = \det\boldsymbol{A}.
$$

For a **diagonal or triangular** matrix, the eigenvalues are simply the diagonal entries.

### 10.2 Continuous time: sign of the real parts

For $\frac{\mathrm{d}}{\mathrm{d}t}(\Delta\boldsymbol{x}) = \boldsymbol{A}\Delta\boldsymbol{x}$ the solution is built from terms $e^{\lambda t}$ (week 3). Hence:

| Eigenvalues of $\boldsymbol{A}$ | Behaviour of $\Delta\boldsymbol{x}$ | Name (2D) |
|---|---|---|
| all $\operatorname{Re}\lambda < 0$ | decays to $\boldsymbol{0}$ | stable node / focus |
| some $\operatorname{Re}\lambda > 0$ | grows | unstable (saddle if the other is $<0$) |
| some $\operatorname{Re}\lambda = 0$, rest $< 0$ | neither decays nor necessarily grows | marginal case |

Quick 2×2 test: **both eigenvalues have negative real part ⇔ $\operatorname{tr}\boldsymbol{A} < 0$ and $\det\boldsymbol{A} > 0$.** $\det\boldsymbol{A} < 0$ ⇔ real eigenvalues of opposite sign ⇔ saddle.

### 10.3 Linearization and local stability of the nonlinear system (preview of week 9)

**Lyapunov's indirect method:** if all eigenvalues of $\boldsymbol{A} = \frac{\partial\boldsymbol{f}}{\partial\boldsymbol{x}}\big|_0$ have $\operatorname{Re}\lambda < 0$, the equilibrium of the *nonlinear* system is locally asymptotically stable; if at least one has $\operatorname{Re}\lambda > 0$, it is unstable. If some eigenvalue has $\operatorname{Re}\lambda = 0$ (and none positive), the linearization is **inconclusive** — the neglected higher-order terms decide. Task 2(c) is such a case.

---

## 11. Discrete-time LTI systems: solution of the state equation

### 11.1 The model

$$
\boldsymbol{x}[k + 1] = \boldsymbol{A}\boldsymbol{x}[k] + \boldsymbol{B}u[k],\qquad \boldsymbol{y}[k] = \boldsymbol{C}\boldsymbol{x}[k] + \boldsymbol{D}u[k],\qquad k\in\mathbb{Z}.
$$

Unlike a differential equation, a difference equation is a **recipe you can run**: given $\boldsymbol{x}[0]$ and the inputs, compute $\boldsymbol{x}[1]$, then $\boldsymbol{x}[2]$, and so on.

Reminder of matrix–vector multiplication: $\begin{bmatrix}a_{11} & a_{12}\\ a_{21} & a_{22}\end{bmatrix}\begin{bmatrix}p\\ q\end{bmatrix} = \begin{bmatrix}a_{11}p + a_{12}q\\ a_{21}p + a_{22}q\end{bmatrix}$ (row times column).

### 11.2 Free response ($u[k] = 0$ for all $k$)

$\boldsymbol{x}[1] = \boldsymbol{A}\boldsymbol{x}[0]$, $\boldsymbol{x}[2] = \boldsymbol{A}\boldsymbol{x}[1] = \boldsymbol{A}^2\boldsymbol{x}[0]$, …, so

$$
\boldsymbol{x}[k] = \boldsymbol{A}^k \boldsymbol{x}[0],\qquad \boldsymbol{A}^0 = \boldsymbol{I}.
$$

### 11.3 Full response — derivation

Iterate with the input:

$$
\begin{aligned}
\boldsymbol{x}[1] &= \boldsymbol{A}\boldsymbol{x}[0] + \boldsymbol{B}u[0]\\
\boldsymbol{x}[2] &= \boldsymbol{A}\boldsymbol{x}[1] + \boldsymbol{B}u[1] = \boldsymbol{A}^2\boldsymbol{x}[0] + \boldsymbol{A}\boldsymbol{B}u[0] + \boldsymbol{B}u[1]\\
\boldsymbol{x}[3] &= \boldsymbol{A}\boldsymbol{x}[2] + \boldsymbol{B}u[2] = \boldsymbol{A}^3\boldsymbol{x}[0] + \boldsymbol{A}^2\boldsymbol{B}u[0] + \boldsymbol{A}\boldsymbol{B}u[1] + \boldsymbol{B}u[2].
\end{aligned}
$$

The pattern: the input applied at step $i$ has been "processed" by $\boldsymbol{A}$ exactly $k - 1 - i$ times by step $k$. Hence

$$
\boxed{\;\boldsymbol{x}[k] = \underbrace{\boldsymbol{A}^k\boldsymbol{x}[0]}_{\text{free (zero-input) response}} + \underbrace{\sum_{i=0}^{k-1}\boldsymbol{A}^{k-1-i}\boldsymbol{B}\,u[i]}_{\text{forced (zero-state) response}}\;},\qquad k\ge1 .
$$

*Proof by induction.* For $k = 1$ it reads $\boldsymbol{x}[1] = \boldsymbol{A}\boldsymbol{x}[0] + \boldsymbol{B}u[0]$ ✓. If it holds for $k$, then
$\boldsymbol{x}[k+1] = \boldsymbol{A}\boldsymbol{x}[k] + \boldsymbol{B}u[k] = \boldsymbol{A}^{k+1}\boldsymbol{x}[0] + \sum_{i=0}^{k-1}\boldsymbol{A}^{k-i}\boldsymbol{B}u[i] + \boldsymbol{A}^0\boldsymbol{B}u[k] = \boldsymbol{A}^{k+1}\boldsymbol{x}[0] + \sum_{i=0}^{k}\boldsymbol{A}^{k-i}\boldsymbol{B}u[i]$ ✓.

This is exactly the split $\boldsymbol{y}(\boldsymbol{x}_0,\boldsymbol{u},t) = \boldsymbol{y}(\boldsymbol{0},\boldsymbol{u},t) + \boldsymbol{y}(\boldsymbol{x}_0,\boldsymbol{0},t)$ of Definition 1.1 (§2.1). The sum is a **discrete convolution** of the input with $\boldsymbol{A}^{j}\boldsymbol{B}$ (week 3).

### 11.4 Impulse response

With $\boldsymbol{x}[0] = \boldsymbol{0}$ and the unit impulse $u[k] = \delta[k]$ ($1$ for $k = 0$, $0$ for $k\ne0$), only the $i = 0$ term survives:

$$
\boldsymbol{x}[k] = \boldsymbol{A}^{k-1}\boldsymbol{B}\ (k\ge1),\quad \boldsymbol{x}[0] = \boldsymbol{0},\qquad \boldsymbol{y}[k] = \boldsymbol{C}\boldsymbol{A}^{k-1}\boldsymbol{B}\ (k\ge1),\quad \boldsymbol{y}[0] = \boldsymbol{D}.
$$

---

## 12. Computing $\boldsymbol{A}^k$ with eigenvalues and eigenvectors

### 12.1 Eigenvalues and eigenvectors

$\boldsymbol{v}\ne\boldsymbol{0}$ is an **eigenvector** of $\boldsymbol{A}$ with **eigenvalue** $\lambda$ if

$$
\boldsymbol{A}\boldsymbol{v} = \lambda\boldsymbol{v}.
$$

$\boldsymbol{A}$ acts on $\boldsymbol{v}$ just by scaling it. Applying $\boldsymbol{A}$ $k$ times: $\boldsymbol{A}^k\boldsymbol{v} = \lambda^k\boldsymbol{v}$. **This is the whole trick.**

*How to compute them:*

1. Eigenvalues: solve $\det(\lambda\boldsymbol{I} - \boldsymbol{A}) = 0$ (§10.1).
2. For each $\lambda$: solve $(\lambda\boldsymbol{I} - \boldsymbol{A})\boldsymbol{v} = \boldsymbol{0}$ (the same matrix as in step 1; $(\boldsymbol{A} - \lambda\boldsymbol{I})\boldsymbol{v} = \boldsymbol{0}$ is equivalent — every sign flips, the solutions don't change). The matrix $\lambda\boldsymbol{I} - \boldsymbol{A}$ is singular, so its rows are multiples of each other; one row $\begin{bmatrix}p & q\end{bmatrix}$ gives $pv_1 + qv_2 = 0$, and $\boldsymbol{v} = \begin{bmatrix}q\\ -p\end{bmatrix}$ (or any multiple) works.
3. Check: $\boldsymbol{A}\boldsymbol{v} = \lambda\boldsymbol{v}$.

### 12.2 Method 1 — decompose the initial state

If $\boldsymbol{A}\in\mathbb{R}^{2\times2}$ has two distinct eigenvalues $\lambda_1\ne\lambda_2$ with eigenvectors $\boldsymbol{v}_1, \boldsymbol{v}_2$, these form a basis: every $\boldsymbol{x}[0]$ can be written as $\boldsymbol{x}[0] = c_1\boldsymbol{v}_1 + c_2\boldsymbol{v}_2$ (solve a 2×2 linear system for $c_1,c_2$). Then

$$
\boldsymbol{x}[k] = \boldsymbol{A}^k\boldsymbol{x}[0] = c_1\lambda_1^k\,\boldsymbol{v}_1 + c_2\lambda_2^k\,\boldsymbol{v}_2 .
$$

### 12.3 Method 2 — diagonalization

Put the eigenvectors as columns in $\boldsymbol{V} = [\boldsymbol{v}_1\ \boldsymbol{v}_2]$ and $\boldsymbol{\Lambda} = \operatorname{diag}(\lambda_1,\lambda_2)$. Then $\boldsymbol{A}\boldsymbol{V} = \boldsymbol{V}\boldsymbol{\Lambda}$, so $\boldsymbol{A} = \boldsymbol{V}\boldsymbol{\Lambda}\boldsymbol{V}^{-1}$ and

$$
\boldsymbol{A}^k = (\boldsymbol{V}\boldsymbol{\Lambda}\boldsymbol{V}^{-1})(\boldsymbol{V}\boldsymbol{\Lambda}\boldsymbol{V}^{-1})\cdots(\boldsymbol{V}\boldsymbol{\Lambda}\boldsymbol{V}^{-1}) = \boldsymbol{V}\boldsymbol{\Lambda}^k\boldsymbol{V}^{-1},\qquad \boldsymbol{\Lambda}^k = \begin{bmatrix}\lambda_1^k & 0\\ 0 & \lambda_2^k\end{bmatrix}.
$$

(All inner $\boldsymbol{V}^{-1}\boldsymbol{V} = \boldsymbol{I}$ cancel.) **2×2 inverse:** $\begin{bmatrix}a_{11} & a_{12}\\ a_{21} & a_{22}\end{bmatrix}^{-1} = \dfrac{1}{a_{11}a_{22} - a_{12}a_{21}}\begin{bmatrix}a_{22} & -a_{12}\\ -a_{21} & a_{11}\end{bmatrix}$.

Equivalently, with $\boldsymbol{w}_i^\top$ the rows of $\boldsymbol{V}^{-1}$: $\boldsymbol{A}^k = \lambda_1^k\,\boldsymbol{v}_1\boldsymbol{w}_1^\top + \lambda_2^k\,\boldsymbol{v}_2\boldsymbol{w}_2^\top$ (sum of rank-one "mode" matrices).

Always **check** a closed form at $k = 0$ (must give $\boldsymbol{I}$ or $\boldsymbol{x}[0]$) and $k = 1$ (must give $\boldsymbol{A}$ or $\boldsymbol{A}\boldsymbol{x}[0]$).

### 12.4 Powers of negative numbers

$(-1)^k$ alternates $1, -1, 1, -1, \dots$ (bounded, never converges). $(-2)^k$ alternates sign and doubles in size: $1, -2, 4, -8, 16,\dots$ In general $\lvert\lambda^k\rvert = \lvert\lambda\rvert^k$.

### 12.5 A useful shortcut

If $\boldsymbol{B}$ (or $\boldsymbol{x}[0]$) happens to be an eigenvector, i.e. $\boldsymbol{A}\boldsymbol{B} = \lambda\boldsymbol{B}$, then $\boldsymbol{A}^j\boldsymbol{B} = \lambda^j\boldsymbol{B}$ for every $j$ — no matrix powers needed. Look for this before computing (Task 3 c, d).

---

## 13. Modes, internal stability and BIBO stability

### 13.1 Modes in discrete time

From §12.2, the free response is a sum of **modes** $c_i\lambda_i^k \boldsymbol{v}_i$. Each mode behaves according to its eigenvalue:

| $\lambda$ (real) | $\lambda^k$ | mode |
|---|---|---|
| $0 < \lambda < 1$ | $\to 0$ monotonically | decays |
| $-1 < \lambda < 0$ | $\to0$ with alternating sign | decays, oscillating |
| $\lambda = 1$ | $1$ | constant |
| $\lambda = -1$ | $\pm1$ alternating | bounded, does **not** converge |
| $\lambda > 1$ | $\to\infty$ | diverges |
| $\lambda < -1$ | $\pm\infty$ alternating | diverges, oscillating |

(For complex $\lambda$, $\lvert\lambda\rvert$ decides the growth and the angle gives a rotation.) **The boundary in discrete time is the unit circle $\lvert\lambda\rvert = 1$**, whereas in continuous time it is the imaginary axis $\operatorname{Re}\lambda = 0$.

### 13.2 Internal stability

The exercise sheet defines it as the behaviour of the state as $k\to\infty$ in the absence of an external input, i.e. of $\boldsymbol{x}[k] = \boldsymbol{A}^k\boldsymbol{x}[0]$:

- **Asymptotically stable:** $\boldsymbol{x}[k]\to\boldsymbol{0}$ for every $\boldsymbol{x}[0]$ ⇔ all $\lvert\lambda_i\rvert < 1$.
- **Marginally stable:** $\boldsymbol{x}[k]$ bounded for every $\boldsymbol{x}[0]$ — requires all $\lvert\lambda_i\rvert \le 1$ (and eigenvalues on the circle not "repeated with a Jordan block" — week 9).
- **Unstable:** some $\boldsymbol{x}[0]$ makes $\boldsymbol{x}[k]$ unbounded ⇔ some $\lvert\lambda_i\rvert>1$ (or a defective eigenvalue on the circle).

For a **specific** $\boldsymbol{x}[0]$, convergence depends on which modes it excites: if $c_i = 0$ for every bad eigenvalue, those modes never appear. The answer to "does it converge for $\boldsymbol{x}[0] = \dots$?" therefore needs the decomposition $\boldsymbol{x}[0] = \sum c_i\boldsymbol{v}_i$.

### 13.3 BIBO stability

BIBO stability is about the input–output behaviour from $\boldsymbol{x}[0] = \boldsymbol{0}$: every bounded input must give a bounded output. For an LTI system it holds if and only if the impulse response is absolutely summable, $\sum_{k}\lVert \boldsymbol{C}\boldsymbol{A}^{k-1}\boldsymbol{B}\rVert < \infty$ (week 10). In particular: **if the response to the bounded input $\delta[k]$ is unbounded, the system is not BIBO stable.**

**Internal vs BIBO.** Internal (asymptotic) stability implies BIBO stability, but not conversely: an unstable mode that the input cannot excite (or the output cannot see) leaves the I/O behaviour bounded. In Task 3 both fail, because $\boldsymbol{B}$ lies exactly on the unstable eigenvector.

---

## 14. Cheat sheet

**System properties (I/O relationship $y(t) = \dots$)**

| Property | Test (lecture definition) | Quick signs of failure |
|---|---|---|
| Linear (Def. 1.1) | $u_1\to y_1,\ u_2\to y_2\Rightarrow\alpha_1u_1 + \alpha_2u_2\to\alpha_1y_1 + \alpha_2y_2$; first check $u = 0\Rightarrow y = 0$ | $u^2$, $x(t)x(t-1)$, $\cos(u)$, $e^u$, constant offsets |
| Time-invariant (Def. 1.2) | $y(u(t), t - T) = y(u(t - T), t)$; on the left replace **every** $t$ by $t - T$ | explicit $t$ outside $u$ (e.g. $\cos 2\pi f_0t$), time scaling $u(2t)$, reversal $x[4-n]$ |
| Causal (Def. 1.3) | $u_1 = u_2\ \forall t\le T\Rightarrow y_1 = y_2\ \forall t\le T$; all read times $\le t$ | $u(t + 1)$, $u(2t)$, $x[4 - n]$ |
| BIBO stable | $\lvert u\rvert\le M\Rightarrow\lvert y\rvert\le M'$ | accumulation, $e^t$, $t\,u(t)$ |

Memoryless ⇒ causal. Time-only factors ($c(t)u(t)$) keep linearity. Changing only time arguments keeps linearity.

**Linearization**

$$
\boldsymbol{f}(\boldsymbol{x}_0,\boldsymbol{u}_0) = \boldsymbol{0},\qquad \frac{\mathrm{d}}{\mathrm{d}t}(\Delta\boldsymbol{x}) = \boldsymbol{A}\Delta\boldsymbol{x} + \boldsymbol{B}\Delta\boldsymbol{u},\quad \Delta\boldsymbol{y} = \boldsymbol{C}\Delta\boldsymbol{x} + \boldsymbol{D}\Delta\boldsymbol{u},\qquad
\boldsymbol{A} = \left.\frac{\partial\boldsymbol{f}}{\partial\boldsymbol{x}}\right|_0,\ \boldsymbol{B} = \left.\frac{\partial\boldsymbol{f}}{\partial\boldsymbol{u}}\right|_0,\ \boldsymbol{C} = \left.\frac{\partial\boldsymbol{g}}{\partial\boldsymbol{x}}\right|_0,\ \boldsymbol{D} = \left.\frac{\partial\boldsymbol{g}}{\partial\boldsymbol{u}}\right|_0 .
$$

Differentiate first, evaluate second. $\sin(k\pi) = 0$, $\cos(k\pi) = (-1)^k$. 2×2: $\lambda^2 - \operatorname{tr}\boldsymbol{A}\,\lambda + \det\boldsymbol{A} = 0$; stable ⇔ $\operatorname{tr}\boldsymbol{A}<0$, $\det\boldsymbol{A}>0$; saddle ⇔ $\det\boldsymbol{A} < 0$.

**Discrete time**

$$
\boldsymbol{x}[k] = \boldsymbol{A}^k\boldsymbol{x}[0] + \sum_{i=0}^{k-1}\boldsymbol{A}^{k-1-i}\boldsymbol{B}u[i],\qquad \boldsymbol{A}^k = \boldsymbol{V}\boldsymbol{\Lambda}^k\boldsymbol{V}^{-1},\qquad \boldsymbol{A}^k\boldsymbol{v} = \lambda^k\boldsymbol{v}.
$$

Asymptotically stable ⇔ all $\lvert\lambda\rvert < 1$. Impulse response $\boldsymbol{x}[k] = \boldsymbol{A}^{k-1}\boldsymbol{B}$; unbounded ⇒ not BIBO stable.
