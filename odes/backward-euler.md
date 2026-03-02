---
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/backward-euler.pdf
    id: odes-backward-euler-pdf
downloads:
  - id: odes-backward-euler-pdf
    title: Download PDF
---

# Backward Euler Method

:::{tip} Big Idea
Backward Euler evaluates $f$ at the *next* time step instead of the
current one, making it **implicit**: $u_{n+1}$ appears on both sides,
requiring a solve at each step. This extra cost buys **unconditional
stability**: the method is stable for any step size when
$\operatorname{Re}(\lambda) < 0$. For stiff problems, where explicit
methods waste enormous effort satisfying stability constraints, this
is essential.
:::

---

## Derivation

Integrate $u' = f(t,u)$ over one time step:

$$
u(t_{n+1}) = u(t_n) + \int_{t_n}^{t_{n+1}} f(s, u(s))\,ds
$$

Forward Euler ([](#def-forward-euler)) approximates this integral by
freezing the integrand at the **left** endpoint $t_n$ (a left rectangle
rule). Backward Euler freezes it at the **right** endpoint $t_{n+1}$
instead:

$$
\int_{t_n}^{t_{n+1}} f(s, u(s))\,ds \approx h\,f(t_{n+1}, u(t_{n+1}))
$$

This gives the **backward Euler** (or implicit Euler) method:

:::{prf:definition} Backward Euler Method
:label: def-backward-euler

$$
u_{n+1} = u_n + h f(t_{n+1}, u_{n+1})
$$

The value $u_{n+1}$, the next approximation of the solution, appears on
both sides: it is defined *implicitly* as the solution to

$$
u_{n+1} - u_n - h f(t_{n+1}, u_{n+1}) = 0
$$
:::

---

## The Cost of Implicitness

Each step of backward Euler requires solving an equation for $u_{n+1}$.

**Linear problems.** For $f(t, u) = \lambda u + g(t)$, we can solve
algebraically:

$$
u_{n+1} = u_n + h(\lambda u_{n+1} + g(t_{n+1}))
\quad\implies\quad
u_{n+1} = \frac{u_n + hg(t_{n+1})}{1 - h\lambda}
$$

No iteration needed.

**Nonlinear problems.** For general $f$, we must solve the nonlinear
equation

$$
G(u_{n+1}) = u_{n+1} - u_n - hf(t_{n+1}, u_{n+1}) = 0
$$

at every time step, using any of the root-finding methods from the
[nonlinear equations chapter](../nonlinear-equations/index.md) (e.g.,
Newton's method or a fixed-point iteration). The previous time step
$u_n$ provides a natural initial guess, and for small $h$ convergence
is typically rapid (1--2 iterations).

See the [Euler's method notebook](../notebooks/euler-method.ipynb) for an
interactive implementation.

---

## Error Analysis and Convergence

Following the same analysis as for
[forward Euler](forward-euler.md), we obtain the same accuracy order.

::::{prf:theorem} Local Truncation Error of Backward Euler
:label: thm-lte-backward-euler

For backward Euler applied to $u' = f(t,u)$ with $u \in C^2$, the local
truncation error is

$$
\tau_n = -\frac{h}{2}u''(t_{n+1}) + O(h^2) = O(h)
$$

Backward Euler has consistency order $p = 1$.

:::{prf:proof}
:class: dropdown

Taylor expand the exact solution about $t_{n+1}$ (not $t_n$, since
backward Euler evaluates $f$ at the next time step):

$$
u(t_n) = u(t_{n+1}) - hu'(t_{n+1}) + \frac{h^2}{2}u''(t_{n+1}) + O(h^3)
$$

Since $u'(t_{n+1}) = f(t_{n+1}, u(t_{n+1}))$, substituting into the
scheme gives

$$
\tau_n = \frac{u(t_{n+1}) - u(t_n)}{h} - f(t_{n+1}, u(t_{n+1}))
= -\frac{h}{2}u''(t_{n+1}) + O(h^2)
$$

The sign of the leading term is opposite to forward Euler's
($+\frac{h}{2}u''$). This affects stability but not the accuracy order.
:::
::::

::::{prf:theorem} Convergence of Backward Euler
:label: thm-convergence-backward-euler

Under the same Lipschitz conditions as [](#thm-one-step-convergence),
the global error of backward Euler satisfies

$$
|E_n| \leq \frac{C}{L}\left(e^{L(t_n - t_0)} - 1\right) h
$$

Backward Euler converges with order 1.

:::{prf:proof}
:class: dropdown

The argument is identical to the forward Euler convergence proof
([](#thm-one-step-convergence)). Backward Euler is a one-step method
with consistency order $p = 1$, so the same error recurrence, Lipschitz
bound, and discrete Grönwall argument apply.
:::
::::

---

## Stability Analysis

Apply backward Euler to the test equation $u' = \lambda u$:

$$
u_{n+1} = u_n + h\lambda\,u_{n+1}
\quad\implies\quad
u_{n+1} = \frac{1}{1 - h\lambda}\,u_n
$$

The stability function is $R(z) = \frac{1}{1-z}$ where $z = h\lambda$.
Compare forward Euler's $R(z) = 1 + z$
([](#def-stability-region-fwd-euler)).

:::{prf:definition} Stability Region of Backward Euler
:label: def-stability-region-bwd-euler

The **stability region** of backward Euler is

$$
\mathcal{S} = \left\{z \in \mathbb{C} : \left|\frac{1}{1-z}\right| \leq 1\right\}
= \{z \in \mathbb{C} : |1-z| \geq 1\}
$$

This is everything **outside** the disk of radius 1 centered at $(1, 0)$.
In particular, it contains the entire left half-plane
$\{\operatorname{Re}(z) \leq 0\}$.


```{figure} ../img/stability_bwd_euler.png
:width: 50%
:align: center

Stability region of backward Euler: everything **outside** the disk $|1 - z| \geq 1$. The entire left half-plane is included, so backward Euler is stable for all $h > 0$ whenever $\operatorname{Re}(\lambda) < 0$.
```

:::

For $\operatorname{Re}(\lambda) < 0$, backward Euler is stable for
**all** $h > 0$: **unconditionally stable**.

Compare forward Euler: $|1+z| \leq 1$ restricts $h$ to a small disk
around $(-1,0)$. Backward Euler: $|1-z| \geq 1$ includes the entire
left half-plane. No matter how large $|\lambda|$ is, we can choose $h$
based on accuracy alone.

---

## Stiffness

Stiff problems ([](#def-stiff-problem)) arise whenever the ODE has
widely separated time scales. The fast components demand small steps for
explicit stability, even after those components have decayed to
negligible levels. This is where backward Euler's unconditional stability
pays off.

:::{prf:example} The Prothero--Robinson Equation
:label: ex-prothero-robinson
:class: dropdown

Consider $u'(t) = \lambda(u - \cos t) - \sin t$ with $u(0) = 1$.

The exact solution is $u(t) = \cos t$, independent of $\lambda$.

With $\lambda = -1000$, forward Euler requires $h < 0.002$ for stability,
but the solution varies on a scale of $O(1)$. The method is forced to
take thousands of tiny steps tracking a transient that has already decayed.

Backward Euler with $h = 0.1$ gives
$|R(z)| = |1/(1 - h\lambda)| = 1/101 \ll 1$: stable and accurate with
100x fewer steps.
:::

:::{prf:example} Heat Equation Semi-Discretization
:label: ex-heat-stiffness
:class: dropdown

Discretize $u_t = u_{xx}$ on $[0,1]$ with spacing $\Delta x$:

$$
\frac{du_j}{dt} = \frac{u_{j+1} - 2u_j + u_{j-1}}{(\Delta x)^2}
$$

This is a system $\mathbf{u}' = A\mathbf{u}$ where $A$ is tridiagonal
with eigenvalues

$$
\lambda_k = -\frac{4}{(\Delta x)^2}\sin^2\!\left(\frac{k\pi}{2(N+1)}\right)
$$

The smallest eigenvalue is $\lambda_1 \approx -\pi^2$ (slow modes) and
the largest is $\lambda_N \approx -4/(\Delta x)^2$ (fast modes). The
ratio $|\lambda_N|/|\lambda_1| \to \infty$ as $\Delta x \to 0$.

**Forward Euler** requires $h \leq (\Delta x)^2/2$. With
$\Delta x = 0.01$: $h \leq 5 \times 10^{-5}$, requiring 20,000 time
steps to reach $T = 1$.

**Backward Euler** has no stability restriction. Choose $h$ based on
accuracy alone; perhaps 100 steps suffice.
:::
