---
kernelspec:
  name: python3
  display_name: Python 3
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/globalization.pdf
    id: nonlinear-systems-globalization-pdf
downloads:
  - id: nonlinear-systems-globalization-pdf
    title: Download PDF
---

# Globalization Strategies

:::{tip} Big Idea
Newton's method converges quadratically, but only if you start close enough.
**Globalization** modifies Newton to converge from farther away by controlling
step sizes. The two main approaches differ in *what they monitor*: the residual
$\|\mathbf{F}(\mathbf{x})\|$ (backward error) or the Newton correction
$\|\Delta\mathbf{x}\|$ (forward error estimate).
:::

## The Problem with Pure Newton

Newton's method has a fundamental tension:

- **Near the solution:** quadratic convergence, take the full Newton step.
- **Far from the solution:** the linear model may be poor, and the full step
  could make things worse.

:::{prf:example} Newton Diverging
:label: ex-newton-diverging-arctan

Consider $f(x) = \arctan(x)$ with root at $x^* = 0$.

Newton's method: $x_{k+1} = x_k - (1+x_k^2)\arctan(x_k)$

Starting from $x_0 = 2$: $x_1 \approx -3.5$, $x_2 \approx 10.2$,
$x_3 \to \infty$. The full Newton step overshoots dramatically.
:::

```{code-cell} python
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt

f = lambda x: np.arctan(x)
fp = lambda x: 1 / (1 + x**2)

# Newton iterates
x_n = [2.0]
for _ in range(3):
    xk = x_n[-1]
    x_n.append(xk - f(xk) / fp(xk))

fig, ax = plt.subplots(figsize=(9, 4.5))
x = np.linspace(-12, 12, 500)
ax.plot(x, f(x), 'b-', linewidth=2)
ax.axhline(0, color='k', linewidth=0.8)

# Plot tangent lines and iterates
colors = ['#d62728', '#ff7f0e', '#2ca02c']
for k in range(min(3, len(x_n) - 1)):
    xk = x_n[k]
    xk1 = x_n[k + 1]
    # Tangent line
    slope = fp(xk)
    tangent = f(xk) + slope * (x - xk)
    ax.plot(x, tangent, '--', color=colors[k], linewidth=1.2, alpha=0.7)
    # Mark x_k on curve
    ax.plot(xk, f(xk), 'o', color=colors[k], markersize=7, zorder=5)
    # Vertical line from x_k to x-axis
    ax.plot([xk, xk], [0, f(xk)], ':', color=colors[k], linewidth=1, alpha=0.5)
    # Mark x_k on x-axis
    ax.plot(xk, 0, 's', color=colors[k], markersize=6, zorder=5)
    # Label
    ax.annotate(f'$x_{k}$', (xk, -0.25), ha='center', fontsize=10, color=colors[k])

# Mark the root
ax.plot(0, 0, 'ko', markersize=8, zorder=5)
ax.annotate('$x^* = 0$', (0.3, -0.25), fontsize=10)

ax.set_xlabel('$x$')
ax.set_ylabel('$f(x)$')
ax.set_title(r"Newton diverges for $f(x) = \arctan(x)$: tangent lines overshoot")
ax.set_xlim(-12, 12)
ax.set_ylim(-2.5, 2.5)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

print("Newton iterates:")
for k, xk in enumerate(x_n):
    print(f"  x_{k} = {xk:.4f}")
```

The fix: **damped Newton**. Instead of taking the full Newton step
$\Delta\mathbf{x}$, take a fraction:

$$
\mathbf{x}_{k+1} = \mathbf{x}_k + \lambda_k \Delta\mathbf{x}_k
$$

where $\lambda_k \in (0, 1]$ is the **damping factor**. The question is how to
choose $\lambda_k$, and this is where the two approaches diverge.

## Approach 1: Backtracking Line Search

The simplest strategy: define a **merit function** that measures progress, and
require each step to decrease it.

For root-finding, the natural merit function is the squared residual:

$$
\phi(\mathbf{x}) = \frac{1}{2}\|\mathbf{F}(\mathbf{x})\|_2^2
$$

Recall from the [error analysis](../nonlinear-equations/error-analysis.md) that
the residual $\|\mathbf{F}(\hat{\mathbf{x}})\|$ is the **backward error**: it
measures how well $\hat{\mathbf{x}}$ satisfies the equation. So backtracking on
$\phi$ amounts to requiring each step to reduce the backward error.

The gradient of $\phi$ is
$\nabla\phi(\mathbf{x}) = D\mathbf{F}(\mathbf{x})^T\mathbf{F}(\mathbf{x})$, connecting root-finding to
optimization.

### The Armijo Condition

Start with $\lambda = 1$ (full Newton step) and reduce until we see
**sufficient decrease** in $\phi$.

:::{prf:definition} The Armijo Condition
:label: def-armijo-condition

Accept step size $\lambda$ if:

$$
\phi(\mathbf{x}_k + \lambda\Delta\mathbf{x}) \leq \phi(\mathbf{x}_k) + c\lambda \nabla\phi(\mathbf{x}_k)^T\Delta\mathbf{x}
$$

where $c \in (0, 1)$ is a small constant (typically $c = 10^{-4}$).
:::

The right-hand side is a linear prediction of $\phi$ along the search
direction, scaled by $c$. Since the Newton direction is a descent direction
($\nabla\phi^T\Delta\mathbf{x} < 0$ when $D\mathbf{F}$ is nonsingular), this
line has negative slope. The small $c$ makes the condition easy to satisfy: we
just need *some* decrease, not optimal decrease.

:::{prf:algorithm} Backtracking Line Search
:label: alg-backtracking-line-search

**Input:** Current point $\mathbf{x}_k$, Newton direction $\Delta\mathbf{x}$,
parameters $c \in (0, 1)$, $\rho \in (0, 1)$

**Output:** Step size $\lambda$

1. $\lambda \gets 1$
2. **while** $\phi(\mathbf{x}_k + \lambda\Delta\mathbf{x}) > \phi(\mathbf{x}_k) + c\lambda \nabla\phi(\mathbf{x}_k)^T\Delta\mathbf{x}$:
3. $\qquad \lambda \gets \rho \lambda$ *(reduce step size)*
4. **return** $\lambda$
:::

Typical values: $c = 10^{-4}$, $\rho = 0.5$.

### Properties

- **Far from solution:** $\lambda_k$ is small, preventing overshoot.
- **Near solution:** $\lambda_k = 1$ (full step accepted), recovering quadratic
  convergence.
- **Guaranteed progress:** each iteration decreases $\phi$.

:::{prf:theorem} Convergence of Backtracking Newton
:label: thm-global-convergence-damped-newton

Under mild conditions on $\mathbf{F}$, damped Newton with Armijo backtracking
satisfies:

$$
\lim_{k \to \infty} \|\nabla\phi(\mathbf{x}_k)\| = 0
$$

If additionally the Jacobian is nonsingular at limit points, the method
converges to a root.
:::

:::{prf:proof}
:class: dropdown

Each step decreases $\phi$ by at least a fixed fraction (the Armijo condition
guarantees this). Since $\phi \geq 0$, the decreases must eventually become
small. This happens only when $\nabla\phi \to 0$, which for our merit function
means $D\mathbf{F}^T\mathbf{F} \to 0$.
:::

:::{prf:remark} What "Global Convergence" Actually Means
:label: rmk-global-convergence-caveat
:class: dropdown

This result is sometimes called "global convergence," meaning convergence from
any starting point (as opposed to Newton's local convergence which requires
starting near the root). However, the guarantee is weaker than it sounds:
$\nabla\phi \to 0$ means $D\mathbf{F}^T\mathbf{F} \to 0$, which could be
satisfied at a stationary point of $\phi$ where $\|\mathbf{F}\| > 0$ (a local
minimum of the residual that is not a root). True convergence to a root requires
the additional condition that $D\mathbf{F}$ is nonsingular at the limit point.
:::

### Limitations

Backtracking on $\|\mathbf{F}\|$ has a subtle flaw: it is **not affine
invariant**. To understand why this matters, we first need to define what affine
invariance means.

## Affine Invariance

:::{prf:definition} Affine Invariance
:label: def-affine-invariance

An iterative method for solving $\mathbf{F}(\mathbf{x}) = \mathbf{0}$ is
**affine invariant** if it produces the same iterates (up to coordinate change)
when applied to the transformed problem
$\tilde{\mathbf{F}}(\mathbf{x}) = A\mathbf{F}(\mathbf{x})$ for any
nonsingular matrix $A$.
:::

Why should we care? The matrix $A$ represents an arbitrary scaling or rotation
of the equations. If we multiply the first equation by 1000, the mathematical
problem is the same, but a non-invariant method may behave completely
differently. An affine invariant method is insensitive to such choices.

:::{prf:remark} Newton's Method is Affine Invariant
:label: rmk-newton-affine-invariant

Pure Newton's method (with $\lambda = 1$) is affine invariant. Under
$\tilde{\mathbf{F}} = A\mathbf{F}$:

$$
D\tilde{\mathbf{F}}^{-1}\tilde{\mathbf{F}}
= (AD\mathbf{F})^{-1}(A\mathbf{F})
= D\mathbf{F}^{-1}\mathbf{F}
$$

The Newton correction $\Delta\mathbf{x} = -D\mathbf{F}^{-1}\mathbf{F}$ is
unchanged. But the residual norm $\|\tilde{\mathbf{F}}\| = \|A\mathbf{F}\|
\neq \|\mathbf{F}\|$ in general.

So **Newton itself is affine invariant, but Armijo backtracking on
$\|\mathbf{F}\|$ breaks this invariance.**
:::

## Approach 2: Error-Oriented Step Control (NLEQ-ERR)

Deuflhard's insight: instead of monitoring the residual $\|\mathbf{F}\|$
(backward error), monitor the **Newton correction**
$\|\Delta\mathbf{x}\|$ (forward error estimate). Since the Newton correction is
affine invariant, this preserves the invariance of the overall method.

### The Natural Level Function

The natural level function for error-oriented Newton is:

$$
h(\mathbf{x}) = \|D\mathbf{F}(\mathbf{x})^{-1}\mathbf{F}(\mathbf{x})\| = \|\Delta\mathbf{x}\|
$$

To see why this estimates the forward error, note that by the mean value theorem:

$$
\mathbf{F}(\mathbf{x}) = \mathbf{F}(\mathbf{x}) - \mathbf{F}(\mathbf{x}^*) = D\mathbf{F}(\xi)(\mathbf{x} - \mathbf{x}^*)
$$

for some $\xi$ between $\mathbf{x}$ and $\mathbf{x}^*$. Inverting gives:

$$
\mathbf{x} - \mathbf{x}^* = D\mathbf{F}(\xi)^{-1}\mathbf{F}(\mathbf{x}) \approx D\mathbf{F}(\mathbf{x})^{-1}\mathbf{F}(\mathbf{x}) = -\Delta\mathbf{x}
$$

So the Newton correction is literally the linearized estimate of the error
$\mathbf{x} - \mathbf{x}^*$, and this estimate improves as
$\mathbf{x} \to \mathbf{x}^*$. Monitoring $\|\Delta\mathbf{x}\|$ is the right
choice because:

1. It estimates what we actually care about (distance to the root).
2. It is affine invariant ($A$ cancels in $D\mathbf{F}^{-1}\mathbf{F}$).
3. Near the root, $\|\Delta\mathbf{x}\| \to 0$ is exactly the convergence
   criterion.

### The Contraction Factor

After a damped step $\mathbf{x}_{k+1} = \mathbf{x}_k + \lambda_k\Delta\mathbf{x}_k$ with $\lambda_k \in (0, 1]$,
compute the **simplified Newton correction** at the trial point using the *old*
Jacobian:

$$
\overline{\Delta\mathbf{x}}_{k+1} = -D\mathbf{F}(\mathbf{x}_k)^{-1}\mathbf{F}(\mathbf{x}_{k+1})
$$

The contraction factor is:

$$
\theta_k = \frac{\|\overline{\Delta\mathbf{x}}_{k+1}\|}{\|\Delta\mathbf{x}_k\|}
$$

This measures how much the error estimate shrank in one step. From the
[Banach fixed point theorem](../nonlinear-equations/fixed-point.md#thm-banach),
convergence requires $\theta_k < 1$ (the iteration must be a contraction). Deuflhard uses a
**restricted monotonicity test**:

$$
\theta_k < 1 - \frac{\lambda_k}{4}
$$

which is stronger than $\theta_k < 1$ and ensures well-controlled convergence.

### Adaptive Step Size Prediction

A key feature of NLEQ_ERR: **predict the optimal $\lambda$ for the next step**
rather than just halving when the test fails.

From the deviation between the actual and predicted corrections:

$$
\mu_k = \frac{\frac{1}{2}\|\Delta\mathbf{x}_k\|\lambda_k^2}{\|\overline{\Delta\mathbf{x}}_{k+1} - (1 - \lambda_k)\Delta\mathbf{x}_k\|}
$$

The next damping factor is $\lambda_{k+1} = \min(1, \mu_k)$. This adapts
aggressively: when the linearization is accurate ($\mu_k$ is large), $\lambda$
ramps up quickly toward 1. When it is poor, $\lambda$ stays small.

### The Algorithm

:::{prf:algorithm} NLEQ_ERR (Simplified)
:label: alg-nleq-err

**Input:** $\mathbf{F}$, $D\mathbf{F}$, initial guess $\mathbf{x}_0$,
tolerance $\varepsilon$

**Output:** Approximate root $\mathbf{x}$

1. Choose initial $\lambda_0 \in (0, 1]$ (based on nonlinearity estimate: $\lambda_0 = 1$ for mildly nonlinear, $\lambda_0 = 10^{-2}$ for highly nonlinear, $\lambda_0 = 10^{-4}$ for extremely nonlinear)
2. **for** $k = 0, 1, 2, \ldots$:
3. $\qquad$ Solve $D\mathbf{F}(\mathbf{x}_k)\Delta\mathbf{x}_k = -\mathbf{F}(\mathbf{x}_k)$
4. $\qquad$ **if** $\|\Delta\mathbf{x}_k\| < \varepsilon$: **return** $\mathbf{x}_k + \Delta\mathbf{x}_k$
5. $\qquad$ **Damping loop:**
6. $\qquad\qquad$ $\mathbf{x}^{\text{trial}} \gets \mathbf{x}_k + \lambda_k\Delta\mathbf{x}_k$
7. $\qquad\qquad$ Solve $D\mathbf{F}(\mathbf{x}_k)\overline{\Delta\mathbf{x}} = -\mathbf{F}(\mathbf{x}^{\text{trial}})$ *(reuse Jacobian)*
8. $\qquad\qquad$ Compute $\theta_k = \|\overline{\Delta\mathbf{x}}\|/\|\Delta\mathbf{x}_k\|$
9. $\qquad\qquad$ **if** $\theta_k < 1 - \lambda_k/4$: accept step, predict $\lambda_{k+1}$ from $\mu_k$
10. $\qquad\qquad$ **else:** $\lambda_k \gets \min(\mu_k, \lambda_k/2)$, repeat damping loop
11. $\qquad$ $\mathbf{x}_{k+1} \gets \mathbf{x}^{\text{trial}}$
:::

Note that step 7 reuses the Jacobian from step 3. The factored Jacobian is
already available, so computing $\overline{\Delta\mathbf{x}}$ requires only a
forward/backward substitution (cheap compared to a new factorization).

:::{prf:theorem} Local Quadratic Convergence Recovery
:label: thm-eventual-unit-step-size

For both methods: if $\mathbf{x}_k \to \mathbf{x}^*$ where
$\mathbf{F}(\mathbf{x}^*) = \mathbf{0}$ and $D\mathbf{F}(\mathbf{x}^*)$ is
nonsingular, then $\lambda_k = 1$ for all $k$ sufficiently large.
Quadratic convergence is recovered near the root.
:::

:::{prf:proof}
:class: dropdown

Near the solution, the Newton step becomes very accurate. For backtracking, the
full step easily satisfies the Armijo condition. For NLEQ_ERR, the contraction
factor $\theta_k \to 0$ (since Newton converges quadratically), so the
monotonicity test is satisfied with $\lambda = 1$.
:::

## Visualization: Newton Fractals

The difference between these approaches is dramatically visible in the complex
plane. Applying Newton's method to $f(z) = z^3 - 1$ from a grid of starting
points and coloring by which root the iteration converges to produces the
classic **Newton fractal**. Backtracking on $|f(z)|$ does not eliminate the
fractal (because $|f|$ has multiple valleys), while NLEQ-ERR's error-oriented
monitoring significantly smooths the basin boundaries.

:::{seealso}
[Newton's Method in the Complex Plane](../notebooks/deflated-continuation.ipynb):
interactive comparison of standard Newton, residual-based damping (NLEQ-RES),
and error-oriented damping (NLEQ-ERR) applied to $z^3 - 1 = 0$.
:::

**Reference:** P. Deuflhard, *Newton Methods for Nonlinear Problems: Affine
Invariance and Adaptive Algorithms*, Springer, 2011.
