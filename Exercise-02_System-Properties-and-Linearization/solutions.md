---
title: "Exercise 2 — Step-by-Step Solutions"
subtitle: "Advanced System Theory · System Properties and Linearization"
---

# Exercise 2 — Step-by-Step Solutions

**Course:** Advanced System Theory (Paderborn University, Automatic Control Group RAT, Prof. Dr. Erdal Kayacan)
**Covers:** Task 1 (a–d), Task 2 (a, b and the I/O-relation part, called "2c" here because the sheet labels it "b" twice), Task 3 (a–d)
**Checking:** no official solution was available when this was written. Every result is verified by the scripts in [`python/`](python/) (symbolic algebra with SymPy, numerical search, simulation).

Every step is written out. References like "Concepts §3.2" point to [`concepts.md`](concepts.md). The notation is that of the lecture slides: bold symbols are vectors and matrices ($\boldsymbol{x}$, $\boldsymbol{f}$, $\boldsymbol{A}$, …), $(\boldsymbol{x}_0,\boldsymbol{u}_0)$ is the equilibrium, $\Delta\boldsymbol{x} = \boldsymbol{x} - \boldsymbol{x}_0$, and $\cdot|_0$ means evaluation at the equilibrium. Symbols given on the exercise sheet ($\phi$, $(\boldsymbol{x}^\ast, u^\ast)$, $\delta[k]$, …) are kept as they are.

---

## Contents

- [Task 1 — System properties](#task-1--system-properties)
- [Task 2a — Linearization](#task-2a--linearization)
- [Task 2b — Linearization of the second nonlinear system](#task-2b--linearization-of-the-second-nonlinear-system)
- [Task 2c — State space and linearization of the I/O relation](#task-2c--state-space-and-linearization-of-the-io-relation)
- [Task 3 — Evolution of a discrete-time system](#task-3--evolution-of-a-discrete-time-system)
- [Summary of results](#summary-of-results)

---

## Task 1 — System properties

For each I/O relationship we test the three definitions of the lecture:

- **Causality (Definition 1.3):** $u_1(t) = u_2(t)\ \forall t\le T \;\Rightarrow\; y_1(t) = y_2(t)\ \forall t\le T$, for every $T > 0$ (Concepts §4.2).
- **Time-invariance (Definition 1.2):** $y\big(u(t), t - T\big) = y\big(u(t - T), t\big)$ for every $T > 0$. The I/O relationships have no initial state, so $\boldsymbol{x}_0$ is dropped. Left side: output delayed by $T$ (replace every $t$ by $t - T$). Right side: output for the delayed input $u(t - T)$ (Concepts §3.2).
- **Linearity (Definition 1.1 / slide 12):** $u_1(t)\to y_1(t)$, $u_2(t)\to y_2(t)$ $\Rightarrow$ $\alpha_1u_1(t) + \alpha_2u_2(t)\to\alpha_1y_1(t) + \alpha_2y_2(t)$ (Concepts §2.3).

A "Y" is proved in general; an "N" is shown with one counterexample (Concepts §1.3). Throughout, $A\ne0$ and $f_0\ne0$ (special cases are noted at the end).

### 1(a) $\;y(t) = u(t) - 0.5u(t + 1)$

**Causality.** The output at time $t$ reads the input at $t$ and at $t + 1$. The read time $t + 1 > t$ lies in the future for every $t$.

*Counterexample.* Take any $T > 0$. Let $u_1(t) = 0$ for all $t$, and $u_2(t) = 0$ for $t\le T$, $u_2(t) = 1$ for $t > T$. Then $u_1(t) = u_2(t)\ \forall t\le T$. At $t = T$:
$y_1(T) = 0$, but $y_2(T) = u_2(T) - 0.5u_2(T+1) = 0 - 0.5 = -0.5$. So $y_1(t)\ne y_2(t)$ for a $t\le T$ ⇒ **not causal**.

**Time-invariance.**

- Right side, input delayed: $y\big(u(t - T), t\big) = u(t - T) - 0.5u(t + 1 - T)$.
- Left side, output delayed (replace every $t$ by $t - T$): $y\big(u(t), t - T\big) = u(t - T) - 0.5u(t - T + 1)$.
- Identical for every $u$, $t$ and $T > 0$ ⇒ **time-invariant**.

**Linearity.** Let $u_1(t)\to y_1(t)$, $u_2(t)\to y_2(t)$ and apply $\alpha_1u_1(t) + \alpha_2u_2(t)$:

$$
\begin{aligned}
&\alpha_1u_1(t) + \alpha_2u_2(t) - 0.5\big(\alpha_1u_1(t+1) + \alpha_2u_2(t+1)\big)\\
&\quad= \alpha_1\underbrace{\big(u_1(t) - 0.5u_1(t+1)\big)}_{y_1(t)} + \alpha_2\underbrace{\big(u_2(t) - 0.5u_2(t+1)\big)}_{y_2(t)} = \alpha_1y_1(t) + \alpha_2y_2(t).
\end{aligned}
$$

⇒ **linear**.

*Interpretation:* a system that "looks one second ahead" — perfectly linear and time-invariant, but it cannot be built as a real-time device.

### 1(b) $\;y(t) = u(2t)$

**Causality.** The read time is $2t$. For $t > 0$, $2t > t$: the output needs future input values. (For $t < 0$ it reads the past, but Definition 1.3 must hold for all $T > 0$.)

*Counterexample.* $T = 1$, $u_1(t) = 0$ for all $t$, $u_2(t) = 0$ for $t\le1$ and $u_2(t) = 1$ for $t > 1$. Then $u_1(t) = u_2(t)\ \forall t\le1$. At $t = 0.75\le1$: $y_1(0.75) = u_1(1.5) = 0$, $y_2(0.75) = u_2(1.5) = 1$ ⇒ **not causal**.

**Time-invariance.** Careful use of the "replace $t$" rule (Concepts §1.2):

- Right side, input delayed: the system reads its input at time $2t$, and the input is now $u(t - T)$, so
  $y\big(u(t - T), t\big) = u(2t - T)$.
- Left side, output delayed: $y\big(u(t), t - T\big) = u\big(2(t - T)\big) = u(2t - 2T)$.
- These differ: the delay at the output is $T/2$ instead of $T$ (the compression also compresses the delay).

*Counterexample.* $u(t)$ = unit step ($u(t) = 1$ for $t\ge0$, else $0$), $T = 1$, $t = 0.75$: the right side gives $u(0.5) = 1$, the left side gives $u(-0.5) = 0$ ⇒ **time-variant**.

**Linearity.** $\alpha_1u_1(t) + \alpha_2u_2(t)$ produces $\alpha_1u_1(2t) + \alpha_2u_2(2t) = \alpha_1y_1(t) + \alpha_2y_2(t)$ ⇒ **linear**. (Changing the time argument never destroys linearity.)

*Interpretation:* playing a recording at double speed. It is linear, but a 1 s delay of the recording shows up as only 0.5 s at the output, and to play at double speed in real time you would need the future of the recording.

### 1(c) $\;y(t) = A\cos(2\pi f_0t + u(t))$

This is **phase modulation**: the input sets the phase of a carrier. ($A$ is the scalar amplitude here, not a matrix.)

**Causality.** Memoryless: $y(t)$ uses $u$ only at the same instant $t$ ⇒ **causal**. (If $u_1(t) = u_2(t)\ \forall t\le T$, then $y_1(t) = A\cos(2\pi f_0t + u_1(t)) = A\cos(2\pi f_0t + u_2(t)) = y_2(t)\ \forall t\le T$.)

**Time-invariance.**

- Right side: $y\big(u(t - T), t\big) = A\cos\big(2\pi f_0t + u(t - T)\big)$ — the carrier term $2\pi f_0t$ is *not* shifted, because it is part of the system, not of the input.
- Left side: $y\big(u(t), t - T\big) = A\cos\big(2\pi f_0(t - T) + u(t - T)\big)$.

*Counterexample.* $u(t) = 0$ for all $t$, $T = \frac{1}{4f_0}$ (a quarter period), $t = 0$: the right side gives $A\cos 0 = A$; the left side gives $A\cos(-\pi/2) = 0$. Different ⇒ **time-variant**. The explicit $t$ in the carrier is the red flag (Concepts §3.3).

**Linearity.** Zero-in-zero-out test (Concepts §2.3): $u(t) = 0$ gives $y(t) = A\cos(2\pi f_0t)\not\equiv0$ ⇒ **not linear**. (Homogeneity also fails in general: $2u(t)$ produces $A\cos(2\pi f_0t + 2u(t))\ne 2A\cos(2\pi f_0t + u(t))$; the input sits inside the cosine.)

### 1(d) $\;y(t) = u(t)\cos(2\pi f_0t + \phi)$

This is **amplitude modulation** (a mixer): the input multiplies a carrier.

**Causality.** Memoryless ⇒ **causal**.

**Time-invariance.**

- Right side: $y\big(u(t - T), t\big) = u(t - T)\cos(2\pi f_0t + \phi)$.
- Left side: $y\big(u(t), t - T\big) = u(t - T)\cos\big(2\pi f_0(t - T) + \phi\big)$.

*Counterexample.* $u(t) = 1$ for all $t$, $T = \frac{1}{4f_0}$, $t = -\frac{\phi}{2\pi f_0}$ (so that $2\pi f_0t + \phi = 0$): the right side gives $\cos0 = 1$, the left side gives $\cos(-\pi/2) = 0$ ⇒ **time-variant**.

**Linearity.** The factor $\cos(2\pi f_0t + \phi)$ is a known function of time that does not depend on $u$:

$$
\big(\alpha_1u_1(t) + \alpha_2u_2(t)\big)\cos(2\pi f_0t + \phi) = \alpha_1u_1(t)\cos(2\pi f_0t + \phi) + \alpha_2u_2(t)\cos(2\pi f_0t + \phi) = \alpha_1y_1(t) + \alpha_2y_2(t)
$$

⇒ **linear** (a linear, time-variant system).

**Compare (c) and (d):** both contain a carrier with an explicit $t$, so both are time-variant. In (d) the input multiplies the cosine (linear); in (c) it sits inside the cosine (nonlinear).

### Summary table and special cases

| | causal | time-invariant | linear |
|---|:---:|:---:|:---:|
| a) $y(t) = u(t) - 0.5u(t + 1)$ | N | Y | Y |
| b) $y(t) = u(2t)$ | N | N | Y |
| c) $y(t) = A\cos(2\pi f_0t + u(t))$ | Y | N | N |
| d) $y(t) = u(t)\cos(2\pi f_0t + \phi)$ | Y | N | Y |

*Special parameter values.* If $f_0 = 0$, the explicit time dependence disappears: (c) becomes $A\cos(u(t))$ and (d) becomes $u(t)\cos\phi$, both time-invariant; (c) stays nonlinear unless $A = 0$ (then $y(t) = 0$, which is trivially linear and time-invariant).

*Extra (not asked):* all four are BIBO stable: $\lvert y\rvert\le1.5M$, $\lvert y\rvert\le M$, $\lvert y\rvert\le\lvert A\rvert$ and $\lvert y\rvert\le M$ respectively.

*Check:* `python/task1_system_properties.py` tests all three properties numerically on random signals and reproduces the table and both counterexamples.

---

## Task 2a — Linearization

$$
\dot x_1(t) = -\sin x_1(t),\qquad \dot x_2(t) = -x_2(t) + \cos x_1(t),
\qquad\text{i.e.}\qquad
\frac{\mathrm{d}}{\mathrm{d}t}\boldsymbol{x} = \boldsymbol{f}(\boldsymbol{x}) = \begin{bmatrix}f_1(\boldsymbol{x})\\ f_2(\boldsymbol{x})\end{bmatrix} = \begin{bmatrix}-\sin x_1\\ -x_2 + \cos x_1\end{bmatrix}.
$$

The system is autonomous (no input), so the linear model is $\frac{\mathrm{d}}{\mathrm{d}t}(\Delta\boldsymbol{x}) = \boldsymbol{A}\Delta\boldsymbol{x}$.

### Step 1 — Equilibrium points (Definition 1.6)

Solve $\boldsymbol{f}(\boldsymbol{x}_0) = \boldsymbol{0}$ (Concepts §8.2). The first component contains only $x_1$:

$$
-\sin x_1 = 0 \iff x_1 = k\pi,\quad k\in\mathbb{Z}.
$$

Substitute into the second: $-x_2 + \cos(k\pi) = 0 \iff x_2 = \cos(k\pi) = (-1)^k$.

$$
\boxed{\;\boldsymbol{x}_{0} = \begin{bmatrix}k\pi\\ (-1)^k\end{bmatrix},\qquad k\in\mathbb{Z}\;}
$$

So, as for the pendulum of slide 18, there is an infinite number of equilibrium points, in two families:

- $k$ even: $(0, 1), (\pm2\pi, 1), (\pm4\pi, 1),\dots$
- $k$ odd: $(\pm\pi, -1), (\pm3\pi, -1),\dots$

*Check $k = 1$:* $f_1 = -\sin\pi = 0$, $f_2 = -(-1) + \cos\pi = 1 - 1 = 0$ ✓.

### Step 2 — Jacobian matrix (symbolic)

$$
\frac{\partial\boldsymbol{f}}{\partial\boldsymbol{x}} = \begin{bmatrix}
\frac{\partial f_1}{\partial x_1} & \frac{\partial f_1}{\partial x_2}\\[2pt]
\frac{\partial f_2}{\partial x_1} & \frac{\partial f_2}{\partial x_2}
\end{bmatrix} = \begin{bmatrix} -\cos x_1 & 0\\ -\sin x_1 & -1\end{bmatrix}.
$$

($\frac{\partial}{\partial x_1}(-\sin x_1) = -\cos x_1$; $f_1$ does not contain $x_2$; $\frac{\partial}{\partial x_1}\cos x_1 = -\sin x_1$; $\frac{\partial}{\partial x_2}(-x_2) = -1$.)

### Step 3 — Evaluate at the equilibrium points

Use $\cos(k\pi) = (-1)^k$, $\sin(k\pi) = 0$:

$$
\boldsymbol{A} = \left.\frac{\partial\boldsymbol{f}}{\partial\boldsymbol{x}}\right|_{\boldsymbol{x}_0} = \begin{bmatrix} -(-1)^k & 0\\ 0 & -1\end{bmatrix}
\quad\Longrightarrow\quad
\boldsymbol{A}_{\text{even}} = \begin{bmatrix}-1 & 0\\ 0 & -1\end{bmatrix},\qquad
\boldsymbol{A}_{\text{odd}} = \begin{bmatrix}1 & 0\\ 0 & -1\end{bmatrix}.
$$

### Step 4 — Linearized models

With $\Delta\boldsymbol{x} = \boldsymbol{x} - \boldsymbol{x}_0$:

$$
\boxed{\;k\ \text{even}:\ \frac{\mathrm{d}}{\mathrm{d}t}\begin{bmatrix}\Delta x_1\\ \Delta x_2\end{bmatrix} = \begin{bmatrix}-1 & 0\\ 0 & -1\end{bmatrix}\begin{bmatrix}\Delta x_1\\ \Delta x_2\end{bmatrix},\qquad
k\ \text{odd}:\ \frac{\mathrm{d}}{\mathrm{d}t}\begin{bmatrix}\Delta x_1\\ \Delta x_2\end{bmatrix} = \begin{bmatrix}1 & 0\\ 0 & -1\end{bmatrix}\begin{bmatrix}\Delta x_1\\ \Delta x_2\end{bmatrix}\;}
$$

where, e.g. for $k = 0$, $\Delta x_1 = x_1$, $\Delta x_2 = x_2 - 1$; for $k = 1$, $\Delta x_1 = x_1 - \pi$, $\Delta x_2 = x_2 + 1$.

### Step 5 — Interpretation (Concepts §10)

The matrices are diagonal, so the eigenvalues are the diagonal entries.

- $k$ even: $\lambda = -1, -1$ ⇒ both negative ⇒ **locally asymptotically stable** (stable node). Small deviations decay like $e^{-t}$.
- $k$ odd: $\lambda = +1, -1$ ⇒ **saddle, unstable**. Deviations in $x_1$ grow like $e^{t}$; deviations in $x_2$ decay.

The structure is the same as for the pendulum: $\dot x_1 = -\sin x_1$ pushes $x_1$ towards the nearest even multiple of $\pi$ and away from odd multiples; $x_2$ simply follows $\cos x_1$ with a first-order lag.

*Check:* `python/task2_linearization.py` (symbolic Jacobians for $k = -2,\dots,2$); figures [`task2a_k0_stable.png`](python/figures/task2a_k0_stable.png) and [`task2a_k1_saddle.png`](python/figures/task2a_k1_saddle.png) overlay nonlinear and linearized trajectories.

![Saddle at (π, −1): nonlinear vs linearized](python/figures/task2a_k1_saddle.png)

---

## Task 2b — Linearization of the second nonlinear system

$$
\begin{aligned}
\dot x_1(t) &= x_1(t)^2 + x_2(t)^2 + x_2(t)\cos x_1(t) &&=: f_1(\boldsymbol{x}),\\
\dot x_2(t) &= (1 + x_1(t))x_1(t) + (1 + x_2(t))x_2(t) + x_1(t)\sin x_2(t) &&=: f_2(\boldsymbol{x}).
\end{aligned}
$$

Multiplied out: $f_2 = x_1 + x_1^2 + x_2 + x_2^2 + x_1\sin x_2$.

### Step 1 — Equilibrium points

We need $\boldsymbol{f}(\boldsymbol{x}_0) = \boldsymbol{0}$:

$$
\text{(I)}\quad x_1^2 + x_2^2 + x_2\cos x_1 = 0,\qquad
\text{(II)}\quad x_1^2 + x_2^2 + x_1 + x_2 + x_1\sin x_2 = 0 .
$$

**Obvious candidate.** $\boldsymbol{x} = \boldsymbol{0}$: (I) $= 0$, (II) $= 0$ ✓.

**Look along $x_1 = 0$.** Then (I) becomes $x_2^2 + x_2 = x_2(x_2 + 1) = 0$ and (II) becomes $x_2^2 + x_2 = 0$ — the same equation. So $x_2 = 0$ or $x_2 = -1$:

$$
(0, 0)\quad\text{and}\quad(0, -1).
$$

*Check $(0,-1)$:* (I) $= 0 + 1 + (-1)\cos0 = 0$ ✓; (II) $= 0 + 1 + 0 - 1 + 0 = 0$ ✓.

**Are there others?** We show there are none (Concepts §8.2, steps 4–5).

*(i) Confine the solutions.* From (I), $x_1^2 + x_2^2 = -x_2\cos x_1 \le \lvert x_2\rvert$ (because $\lvert\cos\rvert\le1$). Hence $x_2^2\le\lvert x_2\rvert$, i.e. $\lvert x_2\rvert\le1$, and $x_1^2\le\lvert x_2\rvert - x_2^2\le\frac14$ (the maximum of $s - s^2$ is $\frac14$), i.e. $\lvert x_1\rvert\le\frac12$. Then $\cos x_1 > 0$, and since $-x_2\cos x_1 = x_1^2 + x_2^2\ge0$ we get $x_2\le0$. So every equilibrium lies in $\lvert x_1\rvert\le\frac12$, $-1\le x_2\le0$.

*(ii) Combine the equations.* Subtract (I) from (II):

$$
x_1 + x_2 + x_1\sin x_2 - x_2\cos x_1 = 0
\quad\Longleftrightarrow\quad
x_1\,(1 + \sin x_2) = -x_2\,(1 - \cos x_1). \tag{III}
$$

In the region from (i): $1 + \sin x_2\ge 1 - \sin 1\approx0.159 > 0$, $-x_2\ge0$, $1 - \cos x_1\ge0$. So the right side is $\ge0$, hence $x_1\ge0$.

*(iii) Rule out $x_1 > 0$.* Suppose $x_1 > 0$. Using $1 - \cos x_1\le x_1^2/2$ and $\lvert x_2\rvert\le1$ in (III):
$x_1(1 + \sin x_2)\le x_1^2/2$, so $1 + \sin x_2\le x_1/2\le\frac14$, which forces $\sin x_2\le-\frac34$, i.e. $x_2\le\arcsin(-0.75)\approx-0.848$. Feeding this back into $x_1^2\le\lvert x_2\rvert(1 - \lvert x_2\rvert)$ (from (i)) gives $x_1^2\le0.848\cdot0.152\approx0.129$, so $x_1\le0.359$; then $1 + \sin x_2\le0.180$, $\sin x_2\le -0.820$, $x_2\le-0.961$, so $x_1^2\le0.961\cdot0.039\approx0.037$, $x_1\le0.194$; then $1 + \sin x_2\le0.097$, i.e. $\sin x_2\le-0.903$. (Each time we use that $s(1 - s)$ decreases for $s > \frac12$, so its largest value on $[\lvert x_2\rvert_{\min}, 1]$ is at the left end.) But on $[-1, 0]$, $\sin x_2\ge\sin(-1)\approx-0.841$. Contradiction. Hence $x_1 = 0$, and we already solved that case.

$$
\boxed{\;\text{equilibrium points: }\ \boldsymbol{x}_{0} = \begin{bmatrix}0\\ 0\end{bmatrix}\ \text{ and }\ \boldsymbol{x}_{0} = \begin{bmatrix}0\\ -1\end{bmatrix}\;}
$$

(The script additionally runs a dense numerical root search on $[-3,3]^2$ and finds exactly these two points.)

### Step 2 — Jacobian matrix (symbolic)

Differentiate term by term (Concepts §9.3 table):

$$
\frac{\partial f_1}{\partial x_1} = 2x_1 - x_2\sin x_1,\quad
\frac{\partial f_1}{\partial x_2} = 2x_2 + \cos x_1,\quad
\frac{\partial f_2}{\partial x_1} = 1 + 2x_1 + \sin x_2,\quad
\frac{\partial f_2}{\partial x_2} = 1 + 2x_2 + x_1\cos x_2 .
$$

$$
\frac{\partial\boldsymbol{f}}{\partial\boldsymbol{x}} = \begin{bmatrix}2x_1 - x_2\sin x_1 & 2x_2 + \cos x_1\\ 1 + 2x_1 + \sin x_2 & 1 + 2x_2 + x_1\cos x_2\end{bmatrix}.
$$

### Step 3 — Evaluate at $\boldsymbol{x}_0 = [0,\ 0]^\top$

$\sin0 = 0$, $\cos0 = 1$:

$$
\boldsymbol{A} = \left.\frac{\partial\boldsymbol{f}}{\partial\boldsymbol{x}}\right|_{(0,0)} = \begin{bmatrix}0 - 0 & 0 + 1\\ 1 + 0 + 0 & 1 + 0 + 0\end{bmatrix} = \begin{bmatrix}0 & 1\\ 1 & 1\end{bmatrix},
\qquad \boxed{\;\frac{\mathrm{d}}{\mathrm{d}t}(\Delta\boldsymbol{x}) = \begin{bmatrix}0 & 1\\ 1 & 1\end{bmatrix}\Delta\boldsymbol{x},\quad \Delta\boldsymbol{x} = \boldsymbol{x}\;}
$$

Eigenvalues: $\operatorname{tr}\boldsymbol{A} = 1$, $\det\boldsymbol{A} = 0\cdot1 - 1\cdot1 = -1$, so $\lambda^2 - \lambda - 1 = 0$, $\lambda = \frac{1\pm\sqrt5}{2}\approx1.618,\ -0.618$. $\det\boldsymbol{A} < 0$ ⇒ **saddle, unstable**.

### Step 4 — Evaluate at $\boldsymbol{x}_0 = [0,\ -1]^\top$

$\sin0 = 0$, $\cos0 = 1$, $\sin(-1) = -\sin1$, $\cos(-1) = \cos1$:

$$
\boldsymbol{A} = \left.\frac{\partial\boldsymbol{f}}{\partial\boldsymbol{x}}\right|_{(0,-1)} = \begin{bmatrix}0 - (-1)\cdot0 & -2 + 1\\ 1 + 0 - \sin1 & 1 - 2 + 0\end{bmatrix}
= \begin{bmatrix}0 & -1\\ 1 - \sin1 & -1\end{bmatrix}
\approx\begin{bmatrix}0 & -1\\ 0.1585 & -1\end{bmatrix}.
$$

$$
\boxed{\;\frac{\mathrm{d}}{\mathrm{d}t}(\Delta\boldsymbol{x}) = \begin{bmatrix}0 & -1\\ 1 - \sin1 & -1\end{bmatrix}\Delta\boldsymbol{x},\qquad \Delta\boldsymbol{x} = \begin{bmatrix}x_1\\ x_2 + 1\end{bmatrix}\;}
$$

Eigenvalues: $\operatorname{tr}\boldsymbol{A} = -1$, $\det\boldsymbol{A} = 0\cdot(-1) - (-1)(1 - \sin1) = 1 - \sin1\approx0.1585$:

$$
\lambda = \frac{-1\pm\sqrt{1 - 4(1 - \sin1)}}{2} = \frac{-1\pm\sqrt{4\sin1 - 3}}{2}\approx -0.198,\ -0.802 .
$$

$\operatorname{tr}\boldsymbol{A} < 0$, $\det\boldsymbol{A} > 0$ ⇒ both negative real ⇒ **stable node, locally asymptotically stable**.

*Check:* `python/task2_linearization.py`; figure [`task2b_stable_node.png`](python/figures/task2b_stable_node.png).

---

## Task 2c — State space and linearization of the I/O relation

$$
\ddot y(t) + (3 + \dot y(t)^2)\,\dot y(t) + (1 + y(t) + y(t)^2)\,u(t) = 0,\qquad \boldsymbol{x} = [x_1,\ x_2]^\top = [y,\ \dot y]^\top.
$$

### Step 1 — Nonlinear state-space representation (Concepts §7.2)

Solve for the highest derivative: $\ddot y = -(3 + \dot y^2)\dot y - (1 + y + y^2)u$. With $\dot x_1 = \dot y = x_2$ and $\dot x_2 = \ddot y$:

$$
\boxed{\;
\frac{\mathrm{d}}{\mathrm{d}t}\boldsymbol{x} = \boldsymbol{f}(\boldsymbol{x},u) = \begin{bmatrix} x_2\\ -(3 + x_2^2)\,x_2 - (1 + x_1 + x_1^2)\,u\end{bmatrix},\qquad
y = g(\boldsymbol{x},u) = x_1\;}
$$

Multiplied out: $f_2 = -3x_2 - x_2^3 - u - x_1u - x_1^2u$.

### Step 2 — Check that the origin is an equilibrium point

$f_1(\boldsymbol{0},0) = 0$, $f_2(\boldsymbol{0},0) = -(3 + 0)\cdot0 - (1 + 0 + 0)\cdot0 = 0$ ✓. So $(\boldsymbol{x}^\ast, u^\ast) = (\boldsymbol{0}, 0)$ is an equilibrium point and $y^\ast = g(\boldsymbol{0},0) = 0$.

*(Bonus: all equilibrium points.)* $f_1 = 0\Rightarrow x_2 = 0$; then $f_2 = -(1 + x_1 + x_1^2)u = 0$. The factor $1 + x_1 + x_1^2 = (x_1 + \tfrac12)^2 + \tfrac34 > 0$ never vanishes, so $u^\ast = 0$. Hence the equilibrium points are $\boldsymbol{x}^\ast = [x_1^\ast,\ 0]^\top$ with **any** $x_1^\ast$ and $u^\ast = 0$ — a whole line. Constant nonzero inputs have no equilibrium.

### Step 3 — Jacobian matrices

$$
\frac{\partial\boldsymbol{f}}{\partial\boldsymbol{x}} = \begin{bmatrix}0 & 1\\ -(1 + 2x_1)\,u & -3 - 3x_2^2\end{bmatrix},\qquad
\frac{\partial\boldsymbol{f}}{\partial u} = \begin{bmatrix}0\\ -(1 + x_1 + x_1^2)\end{bmatrix},\qquad
\frac{\partial g}{\partial\boldsymbol{x}} = \begin{bmatrix}1 & 0\end{bmatrix},\qquad
\frac{\partial g}{\partial u} = 0 .
$$

(For $\frac{\partial f_2}{\partial x_2}$: $\frac{\partial}{\partial x_2}(-3x_2 - x_2^3) = -3 - 3x_2^2$. For $\frac{\partial f_2}{\partial x_1}$: $\frac{\partial}{\partial x_1}\big(-(1 + x_1 + x_1^2)u\big) = -(1 + 2x_1)u$.)

### Step 4 — Evaluate at $(\boldsymbol{x}^\ast,u^\ast) = (\boldsymbol{0},0)$

$$
\boxed{\;
\boldsymbol{A} = \left.\frac{\partial\boldsymbol{f}}{\partial\boldsymbol{x}}\right|_{(\boldsymbol{x}^\ast,u^\ast)} = \begin{bmatrix}0 & 1\\ 0 & -3\end{bmatrix},\quad
\boldsymbol{B} = \left.\frac{\partial\boldsymbol{f}}{\partial u}\right|_{(\boldsymbol{x}^\ast,u^\ast)} = \begin{bmatrix}0\\ -1\end{bmatrix},\quad
\boldsymbol{C} = \begin{bmatrix}1 & 0\end{bmatrix},\quad
\boldsymbol{D} = 0\;}
$$

$$
\frac{\mathrm{d}}{\mathrm{d}t}(\Delta\boldsymbol{x}) = \begin{bmatrix}0 & 1\\ 0 & -3\end{bmatrix}\Delta\boldsymbol{x} + \begin{bmatrix}0\\ -1\end{bmatrix}\Delta u,\qquad \Delta y = \begin{bmatrix}1 & 0\end{bmatrix}\Delta\boldsymbol{x}.
$$

Here $\Delta\boldsymbol{x} = \boldsymbol{x}$, $\Delta u = u$, $\Delta y = y$ because the operating point is the origin. Note that $\frac{\partial f_2}{\partial x_1} = -(1 + 2x_1)u$ vanishes because $u^\ast = 0$: around the origin the dependence of the gain on $y$ is a second-order effect.

**Equivalent linear I/O relation.** $\Delta\dot x_1 = \Delta x_2$, $\Delta\dot x_2 = -3\Delta x_2 - \Delta u$, $\Delta y = \Delta x_1$ give

$$
\Delta\ddot y(t) + 3\,\Delta\dot y(t) + \Delta u(t) = 0 .
$$

This is exactly what you get by dropping all products of small quantities ($\dot y^3$, $yu$, $y^2u$) in the original equation — a good sanity check.

### Step 5 — Interpretation

$\boldsymbol{A}$ is upper triangular, so its eigenvalues are $0$ and $-3$. The eigenvalue $0$ reflects the line of equilibria: shifting $y$ by a constant (with $u = 0$, $\dot y = 0$) gives another equilibrium, so a perturbation of $y$ neither decays nor grows. Because one eigenvalue has zero real part, the linearization alone cannot decide the stability of the nonlinear system (Concepts §10.3).

*Check:* `python/task2_linearization.py` computes $\boldsymbol{A},\boldsymbol{B},\boldsymbol{C},\boldsymbol{D}$ and the eigenvalues symbolically.

---

## Task 3 — Evolution of a discrete-time system

$$
\boldsymbol{x}[k+1] = \underbrace{\begin{bmatrix}0 & -2\\ 1 & -3\end{bmatrix}}_{\boldsymbol{A}}\boldsymbol{x}[k] + \underbrace{\begin{bmatrix}1\\ 1\end{bmatrix}}_{\boldsymbol{B}}u[k],\qquad \boldsymbol{y}[k] = \boldsymbol{x}[k]\quad (\boldsymbol{C} = \boldsymbol{I},\ \boldsymbol{D} = \boldsymbol{0}),\qquad k\in\mathbb{Z}.
$$

### Preparation — eigenvalues and eigenvectors of $\boldsymbol{A}$ (Concepts §12.1)

$\operatorname{tr}\boldsymbol{A} = -3$, $\det\boldsymbol{A} = 0\cdot(-3) - (-2)\cdot1 = 2$:

$$
\det(\lambda\boldsymbol{I} - \boldsymbol{A}) = \lambda^2 + 3\lambda + 2 = (\lambda + 1)(\lambda + 2) = 0\ \Longrightarrow\ \lambda_1 = -1,\quad \lambda_2 = -2 .
$$

Each eigenvector solves $(\lambda_i\boldsymbol{I} - \boldsymbol{A})\boldsymbol{v} = \boldsymbol{0}$ — the same matrix as in the characteristic equation. (Using $\boldsymbol{A} - \lambda_i\boldsymbol{I}$ instead flips every sign and gives the same eigenvectors.)

- $\lambda_1 = -1$: $\lambda_1\boldsymbol{I} - \boldsymbol{A} = \begin{bmatrix}-1 - 0 & 0 + 2\\ 0 - 1 & -1 + 3\end{bmatrix} = \begin{bmatrix}-1 & 2\\ -1 & 2\end{bmatrix}$, the row gives $-v_1 + 2v_2 = 0$ ⇒ $\boldsymbol{v}_1 = \begin{bmatrix}2\\ 1\end{bmatrix}$. Check: $\boldsymbol{A}\boldsymbol{v}_1 = \begin{bmatrix}0 - 2\\ 2 - 3\end{bmatrix} = \begin{bmatrix}-2\\ -1\end{bmatrix} = -\boldsymbol{v}_1$ ✓.
- $\lambda_2 = -2$: $\lambda_2\boldsymbol{I} - \boldsymbol{A} = \begin{bmatrix}-2 - 0 & 0 + 2\\ 0 - 1 & -2 + 3\end{bmatrix} = \begin{bmatrix}-2 & 2\\ -1 & 1\end{bmatrix}$, the row gives $-v_1 + v_2 = 0$ ⇒ $\boldsymbol{v}_2 = \begin{bmatrix}1\\ 1\end{bmatrix}$. Check: $\boldsymbol{A}\boldsymbol{v}_2 = \begin{bmatrix}-2\\ 1 - 3\end{bmatrix} = -2\boldsymbol{v}_2$ ✓.

**Notice:** $\boldsymbol{B} = \boldsymbol{v}_2$. The input enters exactly along the eigenvector of $\lambda_2 = -2$. This will matter in (c) and (d).

### 3(a) — Iteration and closed form ($u[k] = 0$ for all $k$, $\boldsymbol{x}[0] = [1,\ 0]^\top$)

**Iteration** ($\boldsymbol{x}[k+1] = \boldsymbol{A}\boldsymbol{x}[k]$, row times column):

$$
\boldsymbol{x}[1] = \begin{bmatrix}0 & -2\\ 1 & -3\end{bmatrix}\begin{bmatrix}1\\ 0\end{bmatrix} = \begin{bmatrix}0\cdot1 + (-2)\cdot0\\ 1\cdot1 + (-3)\cdot0\end{bmatrix} = \begin{bmatrix}0\\ 1\end{bmatrix},
$$

$$
\boldsymbol{x}[2] = \boldsymbol{A}\begin{bmatrix}0\\ 1\end{bmatrix} = \begin{bmatrix}0 - 2\\ 0 - 3\end{bmatrix} = \begin{bmatrix}-2\\ -3\end{bmatrix},\qquad
\boldsymbol{x}[3] = \boldsymbol{A}\begin{bmatrix}-2\\ -3\end{bmatrix} = \begin{bmatrix}0 + 6\\ -2 + 9\end{bmatrix} = \begin{bmatrix}6\\ 7\end{bmatrix}.
$$

**Closed form, method 1 (decompose $\boldsymbol{x}[0]$, Concepts §12.2).** Write $\boldsymbol{x}[0] = c_1\boldsymbol{v}_1 + c_2\boldsymbol{v}_2$:

$$
\begin{bmatrix}1\\ 0\end{bmatrix} = c_1\begin{bmatrix}2\\ 1\end{bmatrix} + c_2\begin{bmatrix}1\\ 1\end{bmatrix}
\iff \begin{cases}2c_1 + c_2 = 1\\ c_1 + c_2 = 0\end{cases}
\iff c_1 = 1,\ c_2 = -1 .
$$

Since $\boldsymbol{A}^k\boldsymbol{v}_i = \lambda_i^k\boldsymbol{v}_i$:

$$
\boxed{\;\boldsymbol{x}[k] = (-1)^k\begin{bmatrix}2\\ 1\end{bmatrix} - (-2)^k\begin{bmatrix}1\\ 1\end{bmatrix} = \begin{bmatrix}2(-1)^k - (-2)^k\\ (-1)^k - (-2)^k\end{bmatrix}\;}
$$

*Check:* $k = 0$: $[2 - 1,\ 1 - 1]^\top = [1, 0]^\top$ ✓; $k = 1$: $[-2 + 2,\ -1 + 2]^\top = [0, 1]^\top$ ✓; $k = 2$: $[2 - 4,\ 1 - 4]^\top = [-2, -3]^\top$ ✓; $k = 3$: $[-2 + 8,\ -1 + 8]^\top = [6, 7]^\top$ ✓.

**Closed form, method 2 (general $\boldsymbol{A}^k$, Concepts §12.3).** $\boldsymbol{V} = [\boldsymbol{v}_1\ \boldsymbol{v}_2] = \begin{bmatrix}2 & 1\\ 1 & 1\end{bmatrix}$, $\det\boldsymbol{V} = 1$, $\boldsymbol{V}^{-1} = \begin{bmatrix}1 & -1\\ -1 & 2\end{bmatrix}$, $\boldsymbol{\Lambda} = \operatorname{diag}(-1, -2)$:

$$
\boldsymbol{A}^k = \boldsymbol{V}\boldsymbol{\Lambda}^k\boldsymbol{V}^{-1} = (-1)^k\begin{bmatrix}2\\ 1\end{bmatrix}\begin{bmatrix}1 & -1\end{bmatrix} + (-2)^k\begin{bmatrix}1\\ 1\end{bmatrix}\begin{bmatrix}-1 & 2\end{bmatrix}
= (-1)^k\begin{bmatrix}2 & -2\\ 1 & -1\end{bmatrix} + (-2)^k\begin{bmatrix}-1 & 2\\ -1 & 2\end{bmatrix},
$$

$$
\boxed{\;\boldsymbol{A}^k = \begin{bmatrix}2(-1)^k - (-2)^k & -2(-1)^k + 2(-2)^k\\ (-1)^k - (-2)^k & -(-1)^k + 2(-2)^k\end{bmatrix}\;}
$$

*Check:* $k = 0$ gives $\boldsymbol{I}$; $k = 1$ gives $\begin{bmatrix}-2 + 2 & 2 - 4\\ -1 + 2 & 1 - 4\end{bmatrix} = \boldsymbol{A}$ ✓. Multiplying by $\boldsymbol{x}[0] = [1, 0]^\top$ picks the first column — the same result as method 1.

### 3(b) — Convergence and internal stability

**For $\boldsymbol{x}[0] = [1, 0]^\top$:** $\boldsymbol{x}[k]$ contains $-(-2)^k[1,1]^\top$, whose size is $2^k\to\infty$ with alternating sign ($1, -2, 4, -8, \dots$). The other part $(-1)^k[2,1]^\top$ stays bounded but keeps flipping sign. So $\lVert\boldsymbol{x}[k]\rVert\to\infty$: **the state does not converge; it diverges with alternating sign** (as $\boldsymbol{x}[3] = [6, 7]^\top$ already suggests).

**For arbitrary $\boldsymbol{x}[0]\in\mathbb{R}^2$:** decompose $\boldsymbol{x}[0] = c_1\boldsymbol{v}_1 + c_2\boldsymbol{v}_2$ (always possible, $\boldsymbol{v}_1, \boldsymbol{v}_2$ are a basis), so $\boldsymbol{x}[k] = c_1(-1)^k\boldsymbol{v}_1 + c_2(-2)^k\boldsymbol{v}_2$.

| Case | Behaviour of $\boldsymbol{x}[k]$ |
|---|---|
| $c_2\ne0$ (i.e. $\boldsymbol{x}[0]$ not a multiple of $[2,1]^\top$) | diverges ($\lvert\lambda_2\rvert = 2 > 1$) |
| $c_2 = 0$, $c_1\ne0$ ($\boldsymbol{x}[0] = c_1[2,1]^\top$) | bounded but **not convergent**: alternates between $c_1[2,1]^\top$ and $-c_1[2,1]^\top$ ($\lambda_1 = -1$ on the unit circle) |
| $\boldsymbol{x}[0] = \boldsymbol{0}$ | stays $\boldsymbol{0}$ |

Hence **$\boldsymbol{x}[k]$ converges only for $\boldsymbol{x}[0] = \boldsymbol{0}$**. In the language of Concepts §13.2: the largest eigenvalue magnitude is $\max\lvert\lambda_i\rvert = 2 > 1$, so the system is **not internally stable (unstable)**. Even without $\lambda_2$, the eigenvalue $\lambda_1 = -1$ alone would prevent asymptotic stability.

### 3(c) — Closed form with input

From Concepts §11.3:

$$
\boldsymbol{x}[k] = \boldsymbol{A}^k\boldsymbol{x}[0] + \sum_{i=0}^{k-1}\boldsymbol{A}^{k-1-i}\boldsymbol{B}\,u[i],\qquad k\ge1 .
$$

Because $\boldsymbol{B} = \boldsymbol{v}_2$ is an eigenvector with $\lambda_2 = -2$, $\boldsymbol{A}^j\boldsymbol{B} = (-2)^j\boldsymbol{B}$ (Concepts §12.5), so the sum simplifies dramatically:

$$
\boxed{\;\boldsymbol{x}[k] = \boldsymbol{A}^k\boldsymbol{x}[0] + \Big(\sum_{i=0}^{k-1}(-2)^{k-1-i}\,u[i]\Big)\begin{bmatrix}1\\ 1\end{bmatrix},\qquad \boldsymbol{A}^k\ \text{as in 3(a)}\;}
$$

and $\boldsymbol{y}[k] = \boldsymbol{x}[k]$. *Check $k = 1$:* $\boldsymbol{x}[1] = \boldsymbol{A}\boldsymbol{x}[0] + (-2)^0u[0]\boldsymbol{B} = \boldsymbol{A}\boldsymbol{x}[0] + \boldsymbol{B}u[0]$ ✓. *Check $k = 2$:* $\boldsymbol{A}^2\boldsymbol{x}[0] + (-2u[0] + u[1])\boldsymbol{B}$, and directly $\boldsymbol{x}[2] = \boldsymbol{A}\boldsymbol{x}[1] + \boldsymbol{B}u[1] = \boldsymbol{A}^2\boldsymbol{x}[0] + \boldsymbol{A}\boldsymbol{B}u[0] + \boldsymbol{B}u[1] = \boldsymbol{A}^2\boldsymbol{x}[0] - 2\boldsymbol{B}u[0] + \boldsymbol{B}u[1]$ ✓.

For $\boldsymbol{x}[0] = [1,0]^\top$ this reads componentwise:

$$
\boldsymbol{x}[k] = \begin{bmatrix}2(-1)^k - (-2)^k + s_k\\ (-1)^k - (-2)^k + s_k\end{bmatrix},\qquad s_k = \sum_{i=0}^{k-1}(-2)^{k-1-i}u[i].
$$

*Interpretation:* the input only ever excites the $\lambda = -2$ mode; the $\lambda = -1$ mode is untouched by $u$ (the difference $x_1[k] - x_2[k] = (-1)^k(x_1[0] - x_2[0])$ evolves independently of the input).

### 3(d) — Impulse response and BIBO stability

With $\boldsymbol{x}[0] = \boldsymbol{0}$ and $u[k] = \delta[k]$, only the $i = 0$ term of the sum survives:

$$
\boldsymbol{x}[0] = \boldsymbol{0},\qquad \boldsymbol{x}[k] = \boldsymbol{A}^{k-1}\boldsymbol{B} = (-2)^{k-1}\begin{bmatrix}1\\ 1\end{bmatrix}\quad (k\ge1),\qquad \boldsymbol{y}[k] = \boldsymbol{x}[k].
$$

Values: $\boldsymbol{y}[1] = [1, 1]^\top$, $\boldsymbol{y}[2] = [-2, -2]^\top$, $\boldsymbol{y}[3] = [4, 4]^\top$, $\boldsymbol{y}[4] = [-8, -8]^\top$, …, $\lvert y_i[k]\rvert = 2^{k-1}\to\infty$.

The input $\delta[k]$ is bounded ($\lvert u[k]\rvert\le1$), but the output is unbounded. **The output does not remain bounded as $k\to\infty$ ⇒ the system is not BIBO stable.**

**Link to (b).** The unstable eigenvalue $\lambda_2 = -2$ is excited directly by the input (since $\boldsymbol{B} = \boldsymbol{v}_2$) and seen by the output ($\boldsymbol{C} = \boldsymbol{I}$), so internal instability shows up in the I/O behaviour. If $\boldsymbol{B}$ had been along $\boldsymbol{v}_1$ instead, the impulse response would be $(-1)^{k-1}\boldsymbol{v}_1$ — bounded, though not decaying — even though the system is still internally unstable. This is why internal stability and BIBO stability are different concepts (Concepts §13.3).

*Check:* `python/task3_discrete_time.py` verifies all numbers, both closed forms, the forced response against a simulation with random inputs, and plots the responses:

![Free and impulse responses](python/figures/task3_responses.png)

---

## Summary of results

- **Task 1** (causal / time-invariant / linear): a) N Y Y; b) N N Y; c) Y N N; d) Y N Y.
- **Task 2a:** equilibrium points $\boldsymbol{x}_0 = [k\pi,\ (-1)^k]^\top$; $\boldsymbol{A} = \operatorname{diag}(-1, -1)$ for even $k$ (stable), $\boldsymbol{A} = \operatorname{diag}(1, -1)$ for odd $k$ (saddle).
- **Task 2b:** equilibrium points $[0,0]^\top$ and $[0,-1]^\top$ only; $\boldsymbol{A} = \begin{bmatrix}0 & 1\\ 1 & 1\end{bmatrix}$ (saddle, $\lambda = \frac{1\pm\sqrt5}{2}$) and $\boldsymbol{A} = \begin{bmatrix}0 & -1\\ 1 - \sin1 & -1\end{bmatrix}$ (stable node, $\lambda\approx-0.198, -0.802$).
- **Task 2c:** $\dot x_1 = x_2$, $\dot x_2 = -(3 + x_2^2)x_2 - (1 + x_1 + x_1^2)u$, $y = x_1$; at $(\boldsymbol{x}^\ast,u^\ast) = (\boldsymbol{0},0)$: $\boldsymbol{A} = \begin{bmatrix}0 & 1\\ 0 & -3\end{bmatrix}$, $\boldsymbol{B} = \begin{bmatrix}0\\ -1\end{bmatrix}$, $\boldsymbol{C} = [1\ 0]$, $\boldsymbol{D} = 0$, i.e. $\Delta\ddot y + 3\Delta\dot y + \Delta u = 0$.
- **Task 3:** $\boldsymbol{x}[1] = [0,1]^\top$, $\boldsymbol{x}[2] = [-2,-3]^\top$, $\boldsymbol{x}[3] = [6,7]^\top$; $\boldsymbol{x}[k] = (-1)^k[2,1]^\top - (-2)^k[1,1]^\top$; diverges, converges only for $\boldsymbol{x}[0] = \boldsymbol{0}$ (eigenvalues $-1, -2$); $\boldsymbol{x}[k] = \boldsymbol{A}^k\boldsymbol{x}[0] + \sum_{i<k}(-2)^{k-1-i}u[i]\,[1,1]^\top$; impulse response $(-2)^{k-1}[1,1]^\top$ is unbounded ⇒ not BIBO stable.
