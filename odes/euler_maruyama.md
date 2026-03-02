---
authors:
  - name: Andreas Buttenschoen
---

# The Euler–Maruyama Method

:::{tip} Big Idea
Euler–Maruyama extends forward Euler to **stochastic differential equations** (SDEs):
at each step we add a deterministic contribution (the *drift*) and a random
contribution (the *diffusion*). The method is

$$X_{n+1} = X_n + a(t_n, X_n)\,h + b(t_n, X_n)\,\Delta W_n$$

where $\Delta W_n = \sqrt{h}\,Z_n$ with $Z_n \sim \mathcal{N}(0,1)$.
Each simulation produces one possible **sample path**, and
convergence is measured in expectation rather than pointwise.
Setting $b = 0$ recovers forward Euler exactly.
:::

## Why Stochastic?

Many systems we want to model are not purely deterministic:

- **Finance:** Stock prices fluctuate due to unpredictable market forces. The
  Black–Scholes model for asset prices is an SDE, specifically geometric Brownian motion.
- **Biology:** Populations are subject to environmental noise. Deterministic
  logistic growth $dX = rX(1 - X/K)\,dt$ becomes stochastic when we add
  $\sigma X\,dW$ to model random environmental fluctuations.
- **Physics:** A pollen grain suspended in water gets buffeted by water
  molecules, the original Brownian motion observed by Robert Brown in 1827.
  An ODE for the grain's position would predict it stays still (no net force).
  The actual path is a random, jittery curve.

---

## Brownian Motion

### The Drunkard's Walk

Imagine a person standing at the origin who, at regular time intervals
$\delta t$, takes a step of random size. Let $\xi_j$ denote the $j$-th
step, drawn from a probability distribution with density $p(\xi)$
(so that $\mathbb{P}(a \leq \xi_j \leq b) = \int_a^b p(\xi)\,d\xi$).
Suppose the steps are independent with mean zero and variance $\sigma^2$:

$$\mathbb{E}[\xi_j] = \int_{-\infty}^{\infty} \xi\,p(\xi)\,d\xi = 0,
\qquad \operatorname{Var}(\xi_j)
= \int_{-\infty}^{\infty} \xi^2\,p(\xi)\,d\xi = \sigma^2.$$

Here $\mathbb{E}[X]$ denotes the **expected value** (or mean) of a random
variable $X$: the average over all possible outcomes, weighted by their
probabilities.

After $n$ steps the position is $S_n = \sum_{j=1}^n \xi_j$. What are the
mean and variance of $S_n$? Writing out the expected value against the
joint density $p(\xi_1, \ldots, \xi_n)$:

$$\mathbb{E}[S_n]
= \int \cdots \int \left(\sum_{j=1}^n \xi_j\right)
  p(\xi_1, \ldots, \xi_n)\,d\xi_1 \cdots d\xi_n.$$

The integral is linear in the sum. Since the steps are independent, the
joint density factors as
$p(\xi_1, \ldots, \xi_n) = p(\xi_1)\cdots p(\xi_n)$, and the
cross-integrals collapse: $\int p(\xi_k)\,d\xi_k = 1$ for $k \neq j$.
This gives

$$\mathbb{E}[S_n]
= \sum_{j=1}^n \int \xi_j\,p(\xi_j)\,d\xi_j
= \sum_{j=1}^n \mathbb{E}[\xi_j]
= 0.$$

For the variance, independence is essential. Recall that
$\operatorname{Var}(X) = \mathbb{E}[X^2] - (\mathbb{E}[X])^2$. Since
$\mathbb{E}[S_n] = 0$, the variance reduces to
$\operatorname{Var}(S_n) = \mathbb{E}[S_n^2]$, the **mean-square
displacement**: the average of the squared distance from the origin,
taken over many realizations of the walk. Expanding $S_n^2$:

$$\mathbb{E}[S_n^2]
= \mathbb{E}\!\left[\left(\sum_{j=1}^n \xi_j\right)^{\!2}\right]
= \sum_{j=1}^n \mathbb{E}[\xi_j^2]
+ 2\sum_{i < j} \mathbb{E}[\xi_i\,\xi_j].$$

Independence means $\mathbb{E}[\xi_i\,\xi_j] = \mathbb{E}[\xi_i]\,\mathbb{E}[\xi_j] = 0$
for $i \neq j$, so all cross terms vanish and

$$\operatorname{Var}(S_n)
= \sum_{j=1}^n \operatorname{Var}(\xi_j)
= n\sigma^2.$$

After time $t = n\,\delta t$ the mean-square displacement is
$\mathbb{E}[S_n^2] = (t / \delta t)\,\sigma^2$. For this to depend
only on $t$ and not on the discretization $\delta t$, we need
$\sigma^2 \propto \delta t$, i.e. each step has standard deviation
proportional to $\sqrt{\delta t}$.

This is precisely what is observed experimentally: if you track a pollen
grain suspended in water under a microscope, the mean-square displacement
grows linearly with time. Einstein (1905) explained this by the argument
above, connecting the microscopic randomness of molecular collisions to the
macroscopic diffusion rate.

### From the Random Walk to Continuous Noise

Setting $\sigma^2 = \delta t$ and taking the steps to be Gaussian
$\xi_j = \sqrt{\delta t}\,Z_j$ with $Z_j \sim \mathcal{N}(0,1)$, the
random walk updates as

$$W_{j} = W_{j-1} + \sqrt{\delta t}\; Z_j$$

The $\sqrt{\delta t}$ scaling is precisely what Einstein's argument
requires. When $\delta t \to 0$, this
random walk converges to a continuous random function $W(t)$ called
**Brownian motion** (or a **Wiener process**).

::::{prf:definition} Brownian Motion
:label: def-brownian-motion

A stochastic process $W(t)$ is a (standard) **Brownian motion** if it
satisfies:

1. $W(0) = 0$.
2. **Independent increments:** $W(t_3)-W(t_2)$ is independent of
   $W(t_1)-W(t_0)$ for $0 \leq t_0 < t_1 \leq t_2 < t_3$.
3. **Normal increments:** $W(t+h)-W(t) \sim \mathcal{N}(0,h)$.
4. **Continuous paths:** $W(t)$ is a continuous function of $t$, with
   probability 1.

:::{prf:remark} Properties and computation
:class: dropdown

**Simulation.** Property 3 tells us exactly how to simulate $W(t)$ on a
computer. If we want the increment over a time step of size $h$, we
generate $\Delta W = \sqrt{h}\,Z$ where $Z \sim \mathcal{N}(0,1)$. Each
increment is one call to the random number generator, scaled by $\sqrt{h}$.

**Increment statistics.** From property 3, the increment
$\Delta W = W(t+h) - W(t) \sim \mathcal{N}(0,h)$ is a Gaussian random
variable with mean $0$ and variance $h$:

$$\mathbb{E}[\Delta W] = 0 \qquad \text{(the noise has no preferred direction)}$$

$$\mathbb{E}[(\Delta W)^2] = \operatorname{Var}(\Delta W) = h
\qquad \text{(the magnitude scales with the time step)}$$

The second property also follows directly from the drunkard's walk
calculation: the variance of position grows linearly in time.

**Independence.** By property 2, consecutive increments
$\Delta W_0, \Delta W_1, \ldots$ are independent: knowing one tells you
nothing about any other.

:::
::::

<!-- This is a comment that will be hidden in the rendered output.

### Why $\sqrt{h}$ and Not $h$?

This is the single most important scaling in the subject. Consider the
random walk after $n$ independent steps, each of size $\sqrt{\delta t}\,Z_j$.
The position is $\sum_{j=1}^n \sqrt{\delta t}\,Z_j$, and its variance is

$$\operatorname{Var}\!\left(\sum_{j=1}^n \sqrt{\delta t}\,Z_j\right)
= n \cdot \delta t \cdot 1 = n\,\delta t = T$$

where $T = n\,\delta t$ is the total time. The variance grows *linearly
in time*, regardless of the discretization. If we used $\delta t$ instead
of $\sqrt{\delta t}$, the variance would be $n \cdot \delta t^2 =
T \cdot \delta t$, which vanishes as $\delta t \to 0$, so the noise
would disappear in the limit.

:::{prf:remark} Why Stochastic Calculus Differs
:label: rmk-sqrt-scaling

The $\sqrt{\delta t}$ scaling is why stochastic calculus differs from
ordinary calculus. In deterministic numerics, a term of size $\delta t$ is
"first order" and $(\delta t)^2$ is negligible. In stochastic numerics,
$\Delta W \sim \sqrt{\delta t}$, so $(\Delta W)^2 \sim \delta t$, which is
*not* negligible.
:::

-->

### Computing with Stochastic Processes

A stochastic process like $W(t)$ gives a different function on every
realization. The best we can do is describe *statistics* of the process:
expected values, variances, and correlations across realizations. When we
solve an SDE numerically, the solution $X_n$ at each time step is itself a
random variable that depends on all the random increments
$\Delta W_0, \Delta W_1, \ldots, \Delta W_{n-1}$.

In practice, we rarely know the density of $X_n$ explicitly, so we cannot
evaluate $\mathbb{E}[X_n]$ as an integral. Instead, we approximate it by
**Monte Carlo** sampling: run $M$ independent simulations and average the
results,

$$\mathbb{E}[X]
= \int_{-\infty}^{\infty} x\,p(x)\,dx
\approx \frac{1}{M}\sum_{i=1}^{M} X^{(i)},
\qquad X^{(i)} \sim p.$$

This is a quadrature rule for the integral $\int x\,p(x)\,dx$, where the
sample points $X^{(i)}$ are drawn from the distribution $p$ rather than
placed on a deterministic grid. See the companion notebook
[](../notebooks/monte-carlo.ipynb) for Monte Carlo integration examples,
importance sampling, and the Metropolis--Hastings algorithm.

---

## From Euler's Method to Euler–Maruyama

### The Integral Form of an ODE (Duhamel's Principle)

Start from $\frac{dx}{dt} = f(t,x)$. Integrate both sides from $t_n$ to
$t_{n+1}$:

$$x(t_{n+1}) = x(t_n) + \int_{t_n}^{t_{n+1}} f(s, x(s))\,ds$$

This is **Duhamel's principle**. Euler's method approximates the integral by
freezing the integrand at the left endpoint (a rectangle rule):

$$\int_{t_n}^{t_{n+1}} f(s, x(s))\,ds \approx f(t_n, x_n) \cdot h$$

### The SDE and Its Integral Form

:::{prf:definition} Stochastic Differential Equation
:label: def-sde

A stochastic differential equation (SDE) has the form

$$dX = a(t,X)\,dt + b(t,X)\,dW(t), \qquad X(0) = x$$

where $a(t,x)$ is the **drift**, $b(t,x)$ is the **diffusion**, and $W(t)$
is Brownian motion.
:::

This notation is shorthand for an integral equation, Duhamel's principle
with two integrals. Over one time step $[t_n, t_{n+1}]$:

$$X(t_{n+1}) = X(t_n)
+ \underbrace{\int_{t_n}^{t_{n+1}} a(s, X(s))\,ds}_{\text{deterministic integral}}
+ \underbrace{\int_{t_n}^{t_{n+1}} b(s, X(s))\,dW(s)}_{\text{stochastic integral}}$$

We approximate both integrals by the rectangle-rule logic:

- **Deterministic integral:**
  $\int_{t_n}^{t_{n+1}} a(s, X(s))\,ds \approx a(t_n, X_n) \cdot h$.
  This is Euler's approximation.
- **Stochastic integral:**
  $\int_{t_n}^{t_{n+1}} b(s, X(s))\,dW(s) \approx b(t_n, X_n) \cdot \Delta W_n$.
  We freeze $b$ at the left endpoint; what remains is
  $\int_{t_n}^{t_{n+1}} dW(s) = W(t_{n+1}) - W(t_n) = \Delta W_n$.

### The Euler–Maruyama Scheme

:::{prf:algorithm} Euler–Maruyama Method
:label: alg-euler-maruyama

**Input:** Drift $a(t,x)$, diffusion $b(t,x)$, initial value $x$, step
size $h$, number of steps $N$.

**Output:** Approximations $X_0, X_1, \ldots, X_N$.

1. Set $X_0 = x$.
2. For $n = 0, 1, \ldots, N-1$:
   1. Generate $Z_n \sim \mathcal{N}(0,1)$.
   2. Set $\Delta W_n = \sqrt{h}\,Z_n$.
   3. Set $X_{n+1} = X_n + a(t_n, X_n)\,h + b(t_n, X_n)\,\Delta W_n$.
:::

:::{prf:example} Geometric Brownian Motion
:label: ex-gbm
:class: dropdown

The SDE

$$dX = \lambda X\,dt + \mu X\,dW, \qquad X(0) = X_0$$

where $\lambda, \mu$ are constants, is called **geometric Brownian motion**
(GBM). Both drift and diffusion are proportional to $X$. This is the SDE
underlying the Black–Scholes model for asset prices.

**Modelling interpretation.** The relative change in price over a short
interval has two components:
- A deterministic trend $\lambda\,dt$ (expected rate of return).
- A random fluctuation $\mu\,dW$ (volatility).

So $dX/X = \lambda\,dt + \mu\,dW$, i.e. $dX = \lambda X\,dt + \mu X\,dW$.
The multiplicative noise ensures prices stay positive and fluctuations scale
with price level.

**Euler–Maruyama:** With $a(t,X) = \lambda X$ and $b(t,X) = \mu X$:

$$X_{n+1} = X_n + \lambda X_n\,h + \mu X_n\,\Delta W_n
= X_n(1 + \lambda h + \mu\,\Delta W_n)$$

**Exact solution.** This is one of the rare SDEs with a closed-form
solution, making it an ideal test problem: we can compare the
Euler–Maruyama approximation against the true solution on the same
Brownian path. The exact solution is

$$X(t) = X_0 \exp\!\left[(\lambda - \tfrac{1}{2}\mu^2)\,t + \mu\,W(t)\right]$$

The $-\frac{1}{2}\mu^2$ correction is a mathematical consequence of the Itô
correction term (see [](#thm-ito-formula) below), not a modelling choice.
The derivation uses **Itô's formula**, the stochastic chain rule derived in
the next section: set $Y = \ln X$ and apply Itô's formula with
$\varphi(x) = \ln x$ (so $\varphi' = 1/x$, $\varphi'' = -1/x^2$). Here
$f(X) = \lambda X$ and $g(X) = \mu X$, so

$$dY = \frac{1}{X}(\lambda X\,dt + \mu X\,dW)
+ \frac{1}{2}\!\left(-\frac{1}{X^2}\right)(\mu X)^2\,dt
= (\lambda - \tfrac{1}{2}\mu^2)\,dt + \mu\,dW.$$

This has constant coefficients. Integrating from $0$ to $t$:

$$Y(t) - Y(0) = (\lambda - \tfrac{1}{2}\mu^2)\,t + \mu\,W(t).$$

Exponentiating ($X = e^Y$) gives the result.
:::

---

## Convergence of Euler–Maruyama

Both $X_n$ (numerical) and $X(t_n)$ (exact) are random variables. On any
particular run, the error $|X_n - X(t_n)|$ depends on the Brownian path.
We need convergence concepts that account for randomness. Now that we have
the expected value at our disposal, there are two natural ways to measure
the error.

### Strong Convergence: Mean of the Error

Strong convergence asks: **how close is each numerical path to the true
path, on average?** We compute the pathwise error $|X_n - X(t_n)|$ on
each realization, then take the expected value across all realizations.
This matters when individual trajectories are important, for example
when simulating a specific stock price path or a particular particle
trajectory.

:::{prf:definition} Strong Convergence
:label: def-strong-convergence

A method has **strong order of convergence** $\gamma$ if

$$\mathbb{E}\left|X_n - X(t_n)\right| \leq C\,h^{\gamma}$$

for some constant $C$ and sufficiently small $h$.
:::

::::{prf:theorem} Strong Convergence of Euler–Maruyama
:label: thm-em-strong

Under global Lipschitz conditions on the drift $a$ and diffusion $b$,
Euler–Maruyama has **strong order $\gamma = 1/2$**.

:::{prf:proof}
:class: dropdown

We work with the scalar SDE $dX = f(X)\,dt + g(X)\,dW$ assuming $f$ and $g$
are globally Lipschitz with constant $L$ and satisfy a linear growth bound.
Define $e_n = X_n - X(t_n)$.

**Stage 1: Local error.** Suppose we start exactly at $X(t_n)$ and take one
EM step. The local error is

$$\ell_n = \int_{t_n}^{t_{n+1}} [f(X(t_n)) - f(X(s))]\,ds
+ \int_{t_n}^{t_{n+1}} [g(X(t_n)) - g(X(s))]\,dW(s)$$

Using $(a+b)^2 \leq 2a^2 + 2b^2$, Cauchy–Schwarz on the deterministic
integral, and the **Itô isometry** ([](#prop-ito-integral))
$\mathbb{E}|\int \phi\,dW|^2 = \int \mathbb{E}|\phi|^2\,ds$ on the
stochastic integral:

- Deterministic integral contribution: $O(h^3)$ in mean-square.
- Stochastic integral contribution: $O(h^2)$ in mean-square.

The stochastic term dominates:
$\mathbb{E}|\ell_n|^2 \leq C_4\,h^2$.

**Stage 2: Global error accumulation.** Over one step:

$$e_{n+1} = e_n + h[f(X_n) - f(X(t_n))] + [g(X_n) - g(X(t_n))]\Delta W_n + \ell_n$$

Since $\Delta W_n$ is independent of $e_n$ and has zero mean, cross-terms
vanish. Using Lipschitz bounds:

$$\mathbb{E}|e_{n+1}|^2 \leq (1 + C_5 h)\,\mathbb{E}|e_n|^2 + C_6\,h^2$$

By the **discrete Grönwall inequality**, unrolling from $e_0 = 0$:

$$\mathbb{E}|e_n|^2 \leq C\,h$$

By Jensen's inequality $(\mathbb{E}|e_n|)^2 \leq \mathbb{E}|e_n|^2$:

$$\mathbb{E}|X_n - X(t_n)| \leq \sqrt{C}\,h^{1/2} \qquad \square$$
:::
::::

:::{prf:remark} Why Order $1/2$?
:label: rmk-why-half

The bottleneck is the stochastic integral's local error: $O(h^2)$ in
mean-square vs $O(h^3)$ for the deterministic part. The Itô isometry converts
$\int g\,dW$ into $\int |g|^2\,ds$, gaining only one power of $h$ instead
of two. This is a direct consequence of $(dW)^2 = dt$.

Contrast with Euler for ODEs: order 1 (halve the step, halve the error).
For EM: order $1/2$ (halve the step, reduce error by factor $\sqrt{2}$).
:::

### Weak Convergence: Error of the Means

Strong convergence is demanding: it requires every individual path to be
accurate. Often we only care about **statistics** of the solution, for
instance the expected payoff of a financial derivative or the mean
concentration of a chemical species. Weak convergence asks a different
question: **how well does the method reproduce the expected value?**
Rather than averaging the error, we look at the error of the average:
$|\mathbb{E}[X_n] - \mathbb{E}[X(t_n)]|$.

:::{prf:definition} Weak Convergence
:label: def-weak-convergence

A method has **weak order of convergence** $\gamma$ if

$$\left|\mathbb{E}[X_n] - \mathbb{E}[X(t_n)]\right|
\leq C\,h^{\gamma}$$

for some constant $C$ and sufficiently small $h$.
:::

::::{prf:theorem} Weak Convergence of Euler–Maruyama
:label: thm-em-weak

Under sufficient smoothness of $a$ and $b$,
Euler–Maruyama has **weak order $\gamma = 1$**.

:::{prf:proof}
:class: dropdown

**Local error.** Fix $x = X(t_n)$. The exact solution satisfies

$$\mathbb{E}[X(t_{n+1}) \mid X(t_n) = x]
= x + f(x)\,h + O(h^2)$$

(the $O(h^2)$ term comes from the drift's variation over $[t_n, t_{n+1}]$).
The Euler--Maruyama step gives $\hat{X} = x + f(x)\,h + g(x)\,\Delta W$,
so

$$\mathbb{E}[\hat{X}] = x + f(x)\,h$$

since $\mathbb{E}[\Delta W] = 0$. The local weak error is $O(h^2)$.

The pathwise error is dominated by the $g(x)\,\Delta W$ term
($O(\sqrt{h})$ in magnitude), but this has **zero mean** and washes out.

**Global accumulation.** Summing $n = T/h$ local errors of size $O(h^2)$:

$$|\mathbb{E}[X_n] - \mathbb{E}[X(t_n)]|
\leq n \cdot C\,h^2 = CT\,h \qquad \square$$
:::
::::

:::{prf:remark} Why Weak Order Is Higher
:label: rmk-weak-vs-strong

Strong convergence pays the full price of the noise; weak convergence only
pays for the systematic bias. The leading pathwise error has zero mean and
cancels in expectation, giving a full order improvement:
**strong $1/2$, weak $1$**.

More generally, weak convergence holds for
$|\mathbb{E}[\varphi(X_n)] - \mathbb{E}[\varphi(X(t_n))]| \leq C\,h^\gamma$
for any smooth test function $\varphi$. Taking $\varphi(x) = x^k$ gives
convergence of higher moments; the identity $\varphi(x) = x$ is the case
stated above.
:::

---

## Stochastic Calculus and Itô's Formula

To understand where the exact GBM solution comes from, and in particular
the $-\frac{1}{2}\mu^2$ correction, we need a stochastic version of
the chain rule. This requires first making sense of integration with
respect to Brownian motion.

### The Itô Stochastic Integral

We have already seen that the integral form of an SDE involves an expression
$\int_0^T H(s)\,dW(s)$, where $H$ is some process. The **Itô integral**
defines this as the limit of left-endpoint Riemann sums:

:::{prf:definition} Itô Integral
:label: def-ito-integral

For a process $H(t)$ that is adapted (i.e., $H(t)$ depends only on
information up to time $t$), the **Itô integral** is

$$\int_0^T H(t)\,dW(t)
= \lim_{N\to\infty} \sum_{j=0}^{N-1} H(t_j)\,\bigl(W(t_{j+1}) - W(t_j)\bigr)$$

where $0 = t_0 < t_1 < \cdots < t_N = T$ is a partition with mesh going
to zero, and the integrand is evaluated at the **left endpoint** $t_j$ of
each subinterval.
:::

The left-endpoint choice is what makes this the *Itô* integral (as opposed
to a midpoint or right-endpoint convention). It is also the choice that
Euler–Maruyama naturally implements.

Two properties of the Itô integral are essential for everything that follows:

:::{prf:property} Properties of the Itô Integral
:label: prop-ito-integral

1. **Zero mean (martingale property):**
   $\mathbb{E}\!\left[\int_0^T H\,dW\right] = 0$.

2. **Itô isometry:**
   $\mathbb{E}\!\left[\left(\int_0^T H\,dW\right)^{\!2}\right]
   = \int_0^T \mathbb{E}\bigl[H(s)^2\bigr]\,ds$.
:::

### Itô's Formula

Now we can derive the stochastic chain rule. Suppose $X(t)$ satisfies the
SDE $dX = f(X)\,dt + g(X)\,dW$, and let $\varphi$ be a twice continuously
differentiable function. We want the SDE satisfied by $Y(t) = \varphi(X(t))$.

::::{prf:theorem} Itô's Formula
:label: thm-ito-formula

If $dX = f(X)\,dt + g(X)\,dW$ and $\varphi \in C^2$, then

$$d\varphi = \left[\varphi'(X)\,f(X) + \tfrac{1}{2}\,g(X)^2\,\varphi''(X)\right]dt
+ \varphi'(X)\,g(X)\,dW$$

This is the ordinary chain rule $\varphi'(X)\,dX$ plus a correction term
$\frac{1}{2}g^2\,\varphi''\,dt$ that has no deterministic analogue.

:::{prf:proof}
:class: dropdown

Partition $[0,T]$ into subintervals of width $\delta t$. On each subinterval,
write $\Delta X_j = X(t_{j+1}) - X(t_j)$ and Taylor-expand $\varphi$:

$$\varphi(X(t_{j+1})) - \varphi(X(t_j))
= \varphi'(X(t_j))\,\Delta X_j
+ \tfrac{1}{2}\,\varphi''(X(t_j))\,(\Delta X_j)^2
+ O(|\Delta X_j|^3)$$

Sum over all subintervals:

$$\varphi(X(T)) - \varphi(X(0))
= \sum_j \varphi'(X(t_j))\,\Delta X_j
+ \tfrac{1}{2}\sum_j \varphi''(X(t_j))\,(\Delta X_j)^2
+ \text{h.o.t.}$$

**The first sum** is a left-endpoint Riemann sum for $\int_0^T \varphi'(X)\,dX$.
Since $dX = f\,dt + g\,dW$, this converges to

$$\int_0^T \varphi'(X)\,f(X)\,dt + \int_0^T \varphi'(X)\,g(X)\,dW$$

**The second sum** requires care. From the SDE,

$$(\Delta X_j)^2 = \bigl(f(X(t_j))\,\delta t + g(X(t_j))\,\Delta W_j\bigr)^2$$

$$= f^2\,(\delta t)^2 + 2fg\,\delta t\,\Delta W_j + g^2\,(\Delta W_j)^2$$

The first two terms are $O((\delta t)^2)$ and $O((\delta t)^{3/2})$
respectively, both vanishing in the limit. The third term involves
$(\Delta W_j)^2$. Since each $(\Delta W_j)^2$ has mean $\delta t$
(from the Brownian increment properties),
$\sum_j g(X(t_j))^2\,(\Delta W_j)^2$ converges to
$\int_0^T g(X)^2\,dt$ (not zero!). Therefore:

$$\tfrac{1}{2}\sum_j \varphi''(X(t_j))\,(\Delta X_j)^2
\;\longrightarrow\; \tfrac{1}{2}\int_0^T \varphi''(X)\,g(X)^2\,dt$$

**The higher-order terms** involve $|\Delta X_j|^3 \sim (\delta t)^{3/2}$,
so $\sum_j O(|\Delta X_j|^3) = O((\delta t)^{1/2}) \to 0$.

Combining, and writing the result in differential form:

$$d\varphi = \varphi'(X)\,\bigl(f\,dt + g\,dW\bigr)
+ \tfrac{1}{2}\,\varphi''(X)\,g^2\,dt$$

which gives the stated formula. $\square$
:::
::::

:::{prf:remark} The Itô Correction
:label: rmk-ito-correction

In ordinary calculus, the chain rule is $d\varphi = \varphi'(X)\,dX$; the
second-order Taylor term $(dX)^2$ is negligible. For SDEs, $(dX)^2$ contains
a $g^2(dW)^2 = g^2\,dt$ piece that is *first order* and survives. This
produces the correction $\frac{1}{2}g^2\,\varphi''\,dt$.

This is the origin of the $-\frac{1}{2}\mu^2$ term in the GBM solution
([](#ex-gbm)): applying $\varphi = \ln$ gives
$\varphi'' = -1/X^2$, and the correction is
$\frac{1}{2}(\mu X)^2 \cdot (-1/X^2) = -\frac{1}{2}\mu^2$.
:::

---

## Summary

- Euler–Maruyama extends Euler's method to SDEs.
- Numerical solutions are **random variables**; each simulation gives one
  *sample path*.
- The $\sqrt{h}$ scaling of Brownian increments is why $(dW)^2 = dt$, why
  the chain rule gains a correction (Itô's formula), and why strong
  convergence is order $1/2$ instead of $1$.
- **Strong order $1/2$, weak order $1$.**

:::{seealso}
The companion notebook
[](../notebooks/euler-maruyama.ipynb)
contains computational experiments: Brownian motion simulation, Itô vs
Stratonovich integrals, Euler–Maruyama for geometric Brownian motion, and
strong/weak convergence verification. For the underlying algorithms, see
{cite}`Higham2001`.
:::
