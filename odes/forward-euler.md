# Forward Euler Method

:::{tip} Big Idea
Forward Euler is the simplest time-stepping method: replace the derivative with
a forward difference. It's **explicit**—we compute $u_{n+1}$ directly from
$u_n$—making it cheap per step. But this simplicity comes at a cost: forward
Euler is only **conditionally stable**, requiring step sizes small enough that
$|1 + h\lambda| \leq 1$. For stiff problems, this constraint becomes
impractical.
:::

---

## The Initial Value Problem

We seek to approximate the solution of

$$
\frac{du}{dt} = f(t, u), \quad u(t_0) = u_0
$$

on a grid of times $t_0 < t_1 < t_2 < \cdots < t_N = T$ with step size $h = t_{n+1} - t_n$.

Let $u_n$ denote our approximation to $u(t_n)$.

---

## Derivation

Replace $du/dt$ with a forward difference:

$$
\frac{u(t_{n+1}) - u(t_n)}{h} \approx f(t_n, u(t_n))
$$

This gives the **forward Euler** (or explicit Euler) method:

:::{prf:definition} Forward Euler Method
:label: def-forward-euler

$$
u_{n+1} = u_n + h f(t_n, u_n)
$$

Starting from $u_0$, we can march forward in time: given $u_n$, we compute $u_{n+1}$ directly.
:::

:::{prf:remark} Derivation via Duhamel's Principle
:label: rmk-duhamel-fwd-euler
:class: dropdown

An equivalent derivation starts from the integral form of the ODE.
Integrate $u' = f(t,u)$ over one step:

$$
u(t_{n+1}) = u(t_n) + \int_{t_n}^{t_{n+1}} f(s, u(s))\,ds
$$

Approximate the integral by a **left rectangle rule**, freezing the
integrand at $s = t_n$:

$$
\int_{t_n}^{t_{n+1}} f(s, u(s))\,ds \approx h\,f(t_n, u(t_n))
$$

This recovers forward Euler. Freezing at the **right** endpoint
$s = t_{n+1}$ instead gives [backward Euler](backward-euler.md).
Higher-order quadrature rules (midpoint, Simpson, etc.) lead to
higher-order ODE methods.
:::

### Geometric Interpretation

Forward Euler follows the tangent line:
1. At $(t_n, u_n)$, the slope is $f(t_n, u_n)$
2. Follow this slope for time $h$
3. Arrive at $u_{n+1} = u_n + h \cdot \text{slope}$

Each step follows a tangent to a (possibly different) solution curve.

```{figure} ../img/euler.png
:width: 95%
:align: center

**Top:** Forward Euler solution of $u' = 2u$. Each step follows a tangent line to a different solution curve (colored), corresponding to the ODE solution with a different initial condition passing through $(t_n, u_n)$. The numerical solution (black) underestimates the exact solution (blue) because the tangent line falls below the convex exponential. **Bottom:** Backward Euler for the same problem.
```

### Implementation

See the [Euler's Method notebook](../notebooks/euler-method.ipynb) for an interactive implementation with visualizations of convergence, stability, and stiffness.

---

## Sensitivity of the Initial Value Problem

Before analyzing the algorithm's accuracy, we should ask: how sensitive is the *problem* itself?

The colored curves in the figure above are exact solutions of the same ODE with different initial conditions. How quickly do nearby solutions diverge or converge?

### Lipschitz Continuity

:::{prf:definition} Lipschitz Continuity
:label: def-ode-lipschitz

A function $f$ is **Lipschitz continuous** on a domain $D$ if there exists a constant $L \geq 0$ such that

$$
|f(u) - f(v)| \leq L|u - v| \quad \text{for all } u, v \in D
$$

The smallest such $L$ is the **Lipschitz constant**.
:::

:::{prf:remark} Connection to Existence and Uniqueness
:label: rmk-lipschitz-picard
:class: dropdown

This is the same condition that the Picard-Lindelöf theorem requires to guarantee existence and uniqueness of ODE solutions. It is not a coincidence: the convergence proof essentially re-derives the Picard iteration in discrete form. The Lipschitz constant $L$ controls both how quickly nearby solution curves diverge *and* how quickly numerical errors can grow.
:::

:::{prf:example} $C^1$ Functions are Lipschitz
:label: ex-c1-lipschitz
:class: dropdown

If $f(u)$ is continuously differentiable on a bounded domain, then $f$ is Lipschitz with constant $L = \max |f'|$. This follows from the mean value theorem:

$$
|f(u) - f(v)| = |f'(\xi)| \, |u - v| \leq L|u - v|
$$

for some $\xi$ between $u$ and $v$.

For example, $f(u) = \lambda u$ has Lipschitz constant $L = |\lambda|$.
:::

### Condition Number of the IVP

Suppose $u(t)$ and $v(t)$ both solve $u' = f(t,u)$ but with different initial conditions $u(0) = u_0$ and $v(0) = v_0$. Their difference satisfies

$$
\frac{d}{dt}(u - v) = f(t, u) - f(t, v)
$$

If $f$ is Lipschitz continuous in $u$ with constant $L$ ([](#def-ode-lipschitz)), then $\frac{d}{dt}|u - v| \leq L\,|u - v|$. By Grönwall's inequality:

$$
|u(T) - v(T)| \leq e^{LT}\,|u_0 - v_0|
$$

The factor $e^{LT}$ is the **condition number of the IVP**. It is a property of the *problem*, not the algorithm:

- If $L > 0$: nearby solution curves **diverge**. Small perturbations grow exponentially, and $e^{LT}$ can be very large. The problem is **ill-conditioned** over long times.
- If $L < 0$: nearby solution curves **converge**. Perturbations are damped, and $e^{LT} < 1$. The problem is **well-conditioned**.
- If $L = 0$ (e.g., $u' = g(t)$): solution curves are parallel. Errors neither grow nor shrink.

No numerical method can do better than the problem allows. We will see this same factor reappear in the convergence bound ([](#thm-one-step-convergence)).

---

## Local Truncation Error

If we take a single step of forward Euler starting from the *exact* solution, how far off are we? The **local truncation error** measures this per-step accuracy: the error committed in one step, assuming we started with the exact answer.

:::{prf:definition} Local Truncation Error
:label: def-local-truncation-error

Substitute the $\underline{\text{exact solution}}$ $u(t)$ into the numerical scheme. The **local truncation error** $\tau_n$ is the residual:

```{math}
:label: eq-lte-definition
\frac{u(t_{n+1}) - u(t_n)}{h} = f(t_n, u(t_n)) + \tau_n
```

Equivalently, rearranging: if we start at the exact value $u(t_n)$ and take one Euler step, the **one-step error** $\mathcal{L}_n = h\tau_n$ is:

$$
\underbrace{u(t_{n+1}) - \underbrace{\bigl[u(t_n) + h f(t_n, u(t_n))\bigr]}_{\text{one Euler step}}}_{\text{one-step error } \mathcal{L}_n} = h\tau_n
$$

**Warning:** Notice that [](#eq-lte-definition) divides by $h$, so $\tau_n$ is a *rate* (error per unit time), not the actual error committed in one step. The one-step error $\mathcal{L}_n = h\tau_n$ is the actual distance you miss by after one step. For forward Euler: $\tau_n = O(h)$ but $\mathcal{L}_n = O(h^2)$.
:::

:::{prf:proposition} Local Truncation Error of Forward Euler
:label: prop-lte-forward-euler

For forward Euler applied to $u' = f(t,u)$ with $u \in C^2$, the local truncation error is

$$
\tau_n = \frac{h}{2}u''(t_n) + O(h^2) = O(h)
$$

:::{prf:proof}
:class: dropdown

Taylor expand the exact solution about $t_n$:

$$
u(t_{n+1}) = u(t_n) + h u'(t_n) + \frac{h^2}{2}u''(t_n) + O(h^3)
$$

Since $u'(t_n) = f(t_n, u(t_n))$, substituting into [](#eq-lte-definition) gives:

$$
\tau_n = \frac{u(t_{n+1}) - u(t_n)}{h} - f(t_n, u(t_n)) = \frac{h}{2}u''(t_n) + O(h^2)
$$
:::

### Consistency Order

:::{prf:definition} Consistency Order
:label: def-consistency-order

A single-step method has **consistency order** $p$ if for sufficiently smooth $f$:

$$
\|\tau_n\| \leq Ch^{p} \quad \text{for all } h \in (0, h_0] \text{ and all } n
$$

Forward Euler has consistency order $p = 1$.
:::

This is the same notion of "order of accuracy" we saw for [finite differences](../approximation-theory/numerical-differentiation.md): the forward difference $\frac{f(x+h) - f(x)}{h}$ approximates $f'(x)$ with error $O(h)$. Forward Euler *is* a forward difference applied to $du/dt$, so it inherits that same first-order accuracy.

---

## Global Error and Convergence

The local truncation error tells us how accurate a *single step* is. But we take many steps, so the real question is: **do these per-step errors accumulate controllably, or do they blow up?**

:::{prf:definition} Global Error
:label: def-global-error

The **global error** at step $n$ is

$$
E_n = u_n - u(t_n)
$$

the difference between the numerical solution and the exact solution.
:::

The global error is the cumulative effect of all one-step errors. Two things determine its size:

1. **Error accumulation:** Recall that the one-step error is $\mathcal{L}_n = h\tau_n = O(h^{p+1})$. Over $N = T/h$ steps, these pile up: $N \times O(h^{p+1}) = O(h^p)$.
2. **Stability:** Local errors don't just add up; each one gets *amplified* by subsequent steps. If the amplification factor is bounded, errors accumulate gently ($O(h^p)$). If it exceeds 1, errors grow exponentially and the method is useless.


::::{prf:example} Error Propagation for the Linear Test Problem
:label: ex-error-propagation
:class: dropdown

Consider the linear test problem $u' = \lambda u + g(t)$. The exact solution satisfies:

$$
u(t_{n+1}) = (1 + h\lambda)u(t_n) + hg(t_n) + h\tau_n
$$

The numerical solution satisfies:

$$
u_{n+1} = (1 + h\lambda)u_n + hg(t_n)
$$

Subtracting the numerical solution from the exact solution gives the error recurrence:

$$
E_{n+1} = (1 + h\lambda)E_n - h\tau_n
$$

Unrolling this recurrence:

$$
E_n = (1 + h\lambda)^n E_0 - h\sum_{j=0}^{n-1}(1 + h\lambda)^{n-1-j}\tau_j
$$

This is the **discrete Duhamel principle**: the global error at step $n$ equals the initial error propagated forward, plus all local truncation errors propagated by the appropriate powers of the amplification factor $R = 1 + h\lambda$.

Both ingredients are visible here:
- **Accumulation:** Each $\tau_j$ contributes to $E_n$. The sum has $n$ terms.
- **Stability:** Each $\tau_j$ is amplified by $R^{n-1-j}$. If $|R| \leq 1$, these powers are bounded and errors stay controlled. If $|R| > 1$, the earliest errors get amplified the most and the sum explodes.

::::

### Convergence Theorem

:::{prf:theorem} Convergence of One-Step Methods
:label: thm-one-step-convergence

Consider forward Euler ([](#def-forward-euler)) applied to $u' = f(t,u)$, where $f$ is Lipschitz continuous in $u$ with constant $L$ ([](#def-ode-lipschitz)). If the method has consistency order $p$, i.e.,

$$
|\tau_n| \leq Ch^{p} \quad \text{for all } h \in (0, h_0] \text{ and all } n
$$

then the global discretization error satisfies:

$$
|E_n| \leq \frac{C}{L}\left(e^{L(t_n - t_0)} - 1\right) h^p
$$

:::{prf:proof}
:class: dropdown

**Step 1: Error recurrence.**
The numerical method gives $u_{n+1} = u_n + hf(t_n, u_n)$.
The exact solution satisfies $u(t_{n+1}) = u(t_n) + hf(t_n, u(t_n)) + h\tau_n$.

Subtracting the numerical solution from the exact solution:
$$
E_{n+1} = E_n + h[f(t_n, u_n) - f(t_n, u(t_n))] - h\tau_n
$$

**Step 2: Lipschitz bound.**
Since $f$ is Lipschitz with constant $L$:
$$
|f(t_n, u_n) - f(t_n, u(t_n))| \leq L|E_n|
$$

Thus:
$$
|E_{n+1}| \leq |E_n| + hL|E_n| + h|\tau_n| = (1 + hL)|E_n| + h|\tau_n|
$$

**Step 3: Unroll the recurrence.**
Let $\tau_{\max} = \max_j |\tau_j|$. Then:
$$
|E_n| \leq (1 + hL)^n|E_0| + h\tau_{\max}\sum_{j=0}^{n-1}(1 + hL)^j
$$

**Step 4: Bound the geometric sum.**
The inequality $(1 + hL)^k \leq e^{khL}$ follows from the elementary
bound $1 + x \leq e^x$ (valid for all $x \in \mathbb{R}$), raised to
the $k$-th power. The sum is a finite geometric series with ratio
$r = 1 + hL$:
$\sum_{j=0}^{n-1} r^j = \frac{r^n - 1}{r - 1} = \frac{(1+hL)^n - 1}{hL}$.
Applying both estimates:
$$
|E_n| \leq e^{nhL}|E_0| + \frac{e^{nhL} - 1}{L}\tau_{\max}
$$

**Step 5: Substitute $nh = t_n - t_0$ and $\tau_{\max} = O(h^p)$.**
With $E_0 = 0$ (exact initial condition):
$$
|E_n| \leq \frac{e^{L(t_n - t_0)} - 1}{L} \cdot Ch^p = O(h^p)
$$
:::

::::{prf:corollary} Consistency Implies Convergence
:label: cor-consistency-convergence

For one-step methods applied to an ODE with Lipschitz continuous $f$:

$$
\text{Consistency of order } p \implies \text{Convergence of order } p
$$

:::{prf:proof}
:class: dropdown

This follows directly from [](#thm-one-step-convergence). If the method has consistency order $p$, then $|\tau_n| \leq Ch^p$. The theorem gives:

$$
|E_n| \leq \frac{C}{L}\left(e^{LT} - 1\right) h^p = O(h^p)
$$

The key point: the constant $\frac{C}{L}(e^{LT} - 1)$ depends on the problem (through $L$ and $T$) but not on $h$. Thus consistency of order $p$ implies convergence of order $p$.
:::
::::

:::{prf:remark} Backward Error, Forward Error, and the Condition Number
:label: rmk-fwd-bwd-error-ode
:class: dropdown

The convergence theorem is an instance of a pattern that runs through all of numerical analysis:

$$
\underbrace{\text{Global error}}_{\text{what we care about}} \;\leq\; \underbrace{e^{LT}}_{\text{problem sensitivity}} \;\times\; \underbrace{O(h^p)}_{\text{per-step accuracy}}
$$

The three ingredients have names that we will formalize when we study
[linear systems](../qr-least-squares/forward-backward-error.md):

- **Backward error** = the local truncation error $\tau_n = O(h^p)$.
  The numerical solution does not satisfy the original ODE, but it
  *does* satisfy a nearby, perturbed ODE $u' = f(t,u) + \tau(t)$. A
  method with small $\tau$ solves a problem *close* to the original.
- **Forward error** = the global error $E_n = u_n - u(t_n)$, the
  quantity we actually care about: how far is our answer from the truth?
- **Condition number** = $e^{LT}$, the amplification factor. The
  Lipschitz constant $L$ measures how sensitive the ODE is to
  perturbations. This is a property of the *problem*, not the algorithm.

Small backward error does not automatically guarantee small forward
error. It depends on whether the problem amplifies or damps per-step
perturbations.
:::

---

## Stability Analysis

The convergence theorem guarantees that errors vanish as $h \to 0$. But in practice we use a *fixed* step size. Can errors grow from step to step?

The error propagation example ([](#ex-error-propagation)) already answered this for the linear problem $u' = \lambda u + g(t)$: the error at each step is multiplied by the **amplification factor** $R = 1 + h\lambda$. If $|R| > 1$, errors grow exponentially and the method is useless.

This is why the scalar **test equation** $u' = \lambda u$ plays a central role. Any smooth ODE is locally linear (by Taylor expansion), so the behavior of a method on $u' = \lambda u$ reveals its stability properties in general. Applying forward Euler gives

$$
u_{n+1} = (1 + h\lambda)\,u_n
$$

so the numerical solution is $u_n = (1 + h\lambda)^n u_0$. For the solution to remain bounded, we need $|1 + h\lambda| \leq 1$. Writing $z = h\lambda$, this defines the **stability region**.

:::{prf:definition} Stability Region
:label: def-stability-region-fwd-euler

The **stability region** of forward Euler is:

$$
\mathcal{S} = \{z \in \mathbb{C} : |1 + z| \leq 1\}
$$

This is a disk of radius 1 centered at $(-1, 0)$ in the complex plane.

```{figure} ../img/stability_fwd_euler.png
:width: 50%
:align: center

Stability region of forward Euler: the disk $|1 + z| \leq 1$. For stability, $z = h\lambda$ must lie inside this region.
```

:::

For real $\lambda < 0$, stability requires $-2 \leq h\lambda \leq 0$, giving the step-size restriction

$$
h \leq \frac{2}{|\lambda|}
$$

Forward Euler is **conditionally stable**: the step size is restricted by the problem's eigenvalues.

### Stiffness

When $|\lambda|$ is large, the stability constraint $h \leq 2/|\lambda|$
forces extremely small step sizes even when the solution itself varies
slowly.

:::{prf:definition} Stiff Problem
:label: def-stiff-problem

A problem is **stiff** when an explicit method must use a step size
far smaller than what accuracy alone would require. Two step sizes are
in play:

- $h_{\text{accuracy}}$: the largest $h$ that resolves the solution to
  the desired tolerance. This is set by the truncation error and how
  rapidly the solution varies.
- $h_{\text{stability}}$: the largest $h$ for which the method does not
  blow up. This is set by the stability region and the eigenvalues of
  the problem.

The problem is stiff when

$$
h_{\text{stability}} \ll h_{\text{accuracy}}
$$

so that stability, not accuracy, dictates the step size.
:::

:::{prf:example} Stiffness in Action
:label: ex-stiffness-fwd-euler

Consider $u' = \lambda(u - \cos t) - \sin t$ with $u(0) = 1$ and
$\lambda = -1000$. The exact solution is $u(t) = \cos t$ regardless of
$\lambda$: a smooth, slowly varying function that a step size of
$h = 0.1$ would track accurately ($h_{\text{accuracy}} \approx 0.1$).
Yet forward Euler's stability constraint demands
$h < 2/|\lambda| = 0.002$ ($h_{\text{stability}} = 0.002$). To
integrate to $T = 1$ we need at least 500 steps, all dictated by
stability rather than accuracy.

:::

```{figure} ../img/stiffness.png
:width: 95%
:align: center

Forward Euler applied to $u' = \lambda(u - \cos t) - \sin t$ with fixed $h = 10^{-3}$. The exact solution is $u(t) = \cos t$ (dashed) for all $\lambda$. **Left/center:** $\lambda = 0$ and $\lambda = -10$ satisfy the stability condition $h < 2/|\lambda|$ and track the solution accurately. **Right:** $\lambda = -2100$ violates the stability condition ($h|\lambda| = 2.1 > 2$) and the numerical solution oscillates with exponentially growing amplitude.
```
---

## Limitations

**Low order.** Forward Euler is only $O(h)$ accurate. Just as central
differences improved on forward differences by averaging, we can build
higher-order ODE methods by evaluating $f$ at cleverly chosen intermediate
points and combining the results. This is the idea behind **Runge-Kutta
methods**; see the
[optional section on Runge-Kutta and adaptive methods](runge-kutta.md).

**Fixed step size.** For a nonlinear ODE, the effective $\lambda$ changes along
the solution. A step size that is stable in one region may be unstable in
another, and a step size that is accurate during rapid transients may be
wasteful during slow evolution. With a fixed $h$, there is no way to adapt. This
motivates **adaptive time-stepping**, where we estimate the local truncation
error at each step and adjust $h$ automatically; see the
[optional section on Runge-Kutta and adaptive methods](runge-kutta.md).

