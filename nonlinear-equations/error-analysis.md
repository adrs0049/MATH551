---
kernelspec:
  name: python3
  display_name: Python 3
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/error-analysis.pdf
    id: nonlinear-equations-error-analysis-pdf
downloads:
  - id: nonlinear-equations-error-analysis-pdf
    title: Download PDF
---

# Error Analysis for Root Finding

:::{tip} Big Idea
A computed root is never exact. The question is: how wrong can it be? **Condition numbers** tell us how sensitive the root is to perturbations in the problem, and **forward/backward error** gives us a framework for interpreting what our algorithm actually computed.
:::

## The Absolute Condition Number

Before analyzing specific root-finding methods, we need to understand a more
fundamental question: **how sensitive is the root itself to perturbations in the
problem?** We introduced [condition numbers](#def-condition-number) and
[stability](#def-stable-algorithm) in the error and stability chapter. Here we
apply these ideas specifically to root finding. This is independent of the algorithm used. Even with exact arithmetic,
if the problem data is slightly wrong (due to measurement error, rounding, or
modelling), the root will shift. The condition number quantifies this
sensitivity.

For root finding, we work with the **absolute condition number** rather than the
relative condition number. The reason is fundamental: the input to a root-finding
problem is the function value $f(x^*) = 0$, and relative error is undefined when
the quantity is zero. Since we are perturbing around zero, only absolute
perturbations make sense.

:::{prf:definition} Absolute Condition Number
:label: def-absolute-condition-number

Consider a problem that maps input data $x$ to output $y = F(x)$. The
**absolute condition number** is the worst-case ratio of absolute output change
to absolute input change for small perturbations:

$$
\hat{\kappa} = \lim_{\delta \to 0} \sup_{0 < |\Delta x| \leq \delta} \frac{|F(x + \Delta x) - F(x)|}{|\Delta x|}
$$

When $F$ is differentiable, this is exactly $\hat{\kappa} = |F'(x)|$.

A problem is **well-conditioned** if $\hat{\kappa}$ is small and
**ill-conditioned** if $\hat{\kappa}$ is large.
:::

:::{prf:remark}
:label: rmk-abs-cond-limit
:class: dropdown

Unlike the [relative condition number](#def-condition-number), the limit here
cannot be avoided. In the relative case, the perturbation $|\Delta x|$ cancels
between numerator and denominator, keeping the ratio bounded. For the absolute
condition number, no such cancellation occurs: if $F$ is nonlinear, the ratio
$|F(x + \Delta x) - F(x)|/|\Delta x|$ can grow without bound as $|\Delta x|$
increases. The limit restricts attention to infinitesimal perturbations, where
the ratio reduces to $|F'(x)|$.
:::

## The Condition Number for Roots

Now we apply this to root finding. Suppose $f(x^*) = 0$ and we perturb $f$
slightly, say to $\tilde{f} = f + \epsilon g$ for some bounded function $g$.
The perturbed equation $\tilde{f}(\tilde{x}^*) = 0$ has a shifted root
$\tilde{x}^*$. How large is the shift $|\tilde{x}^* - x^*|$ relative to the
perturbation size $\epsilon$?

:::{prf:theorem} Condition Number for a Simple Root
:label: thm-condition-number-simple-root

If $f(x^*) = 0$ and $f'(x^*) \neq 0$ (a **simple root**), then the absolute condition number of the root with respect to perturbations of $f$ is:

$$
\hat{\kappa} = \frac{1}{|f'(x^*)|}
$$
:::

:::{prf:proof}
:class: dropdown

Perturb $f$ to $\tilde{f} = f + \epsilon g$. The perturbed root $\tilde{x}^*$ satisfies $\tilde{f}(\tilde{x}^*) = 0$, so:

$$
f(\tilde{x}^*) + \epsilon g(\tilde{x}^*) = 0
$$

Expand $f(\tilde{x}^*)$ in a Taylor series around $x^*$:

$$
f(x^*) + f'(x^*)(\tilde{x}^* - x^*) + \mathcal{O}((\tilde{x}^* - x^*)^2) + \epsilon g(\tilde{x}^*) = 0
$$

Since $f(x^*) = 0$:

$$
f'(x^*)(\tilde{x}^* - x^*) \approx -\epsilon g(\tilde{x}^*)
$$

Therefore:

$$
\tilde{x}^* - x^* \approx -\frac{\epsilon g(x^*)}{f'(x^*)}
$$

Taking absolute values:

$$
|\tilde{x}^* - x^*| \approx \frac{|g(x^*)|}{|f'(x^*)|} \cdot \epsilon
$$

The perturbation to $f$ has amplitude $\epsilon |g(x^*)|$, so the condition
number is the ratio of the output change $|\tilde{x}^* - x^*|$ to the input
change $\epsilon |g(x^*)|$:

$$
\hat{\kappa} = \frac{|\tilde{x}^* - x^*|}{\epsilon |g(x^*)|} \approx \frac{1}{|f'(x^*)|}
$$

Note that $g$ cancels entirely: the condition number depends only on $f'(x^*)$,
not on the shape of the perturbation.
:::

:::{prf:remark} Connection to the Lipschitz Constant
:label: rmk-condition-lipschitz
:class: dropdown

The quantity $|f'(x^*)|$ is the local Lipschitz constant of $f$ near the root:
it controls how much $f$ changes per unit change in $x$. The condition number
$\hat{\kappa} = 1/|f'(x^*)|$ is its reciprocal, which is the Lipschitz constant
of the *inverse* map from perturbations in $f$ to shifts in the root. A large
Lipschitz constant for $f$ (steep crossing) means a small condition number. A
small Lipschitz constant for $f$ (shallow crossing) means a large condition
number.
:::

:::{prf:remark} Geometric Interpretation
:label: rmk-condition-number-geometry
:class: dropdown

The condition number has a simple geometric meaning:

- **Well-conditioned** ($\kappa$ small): $f$ crosses zero steeply. A small vertical perturbation of the curve barely moves the root horizontally.
- **Ill-conditioned** ($\kappa$ large): $f$ crosses zero at a shallow angle. The same small vertical perturbation slides the root far along the $x$-axis.

This is why nearly-tangent crossings are dangerous: the root sits on a nearly flat part of $f$, so it's easy to perturb.
:::

```{code-cell} python
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
eps = 0.3

# --- Left: Well-conditioned (steep crossing) ---
ax = axes[0]
x = np.linspace(-1, 3, 300)
f_well = lambda x: 2*(x - 1)
ax.plot(x, f_well(x), 'b-', linewidth=2)
ax.plot(x, f_well(x) + eps, 'b--', linewidth=1.5, alpha=0.6)
ax.plot(x, f_well(x) - eps, 'b--', linewidth=1.5, alpha=0.6)
ax.axhline(0, color='k', linewidth=0.8)

# Original root
root = 1.0
ax.plot(root, 0, 'ko', markersize=8, zorder=5)

# Perturbed roots
root_plus = 1 - eps/2
root_minus = 1 + eps/2
ax.plot(root_plus, 0, 'rs', markersize=7, zorder=5)
ax.plot(root_minus, 0, 'rs', markersize=7, zorder=5)

# Show horizontal shift
ax.annotate('', xy=(root_minus, -0.4), xytext=(root, -0.4),
            arrowprops=dict(arrowstyle='<->', color='red', lw=1.5))
ax.text((root + root_minus)/2, -0.7, r'$\Delta x^*$', ha='center',
        fontsize=11, color='red')

# Show vertical perturbation
ax.annotate('', xy=(-0.3, eps), xytext=(-0.3, 0),
            arrowprops=dict(arrowstyle='<->', color='gray', lw=1.2))
ax.text(-0.55, eps/2, r'$\epsilon$', ha='center', fontsize=11, color='gray')

ax.fill_between(x, f_well(x) - eps, f_well(x) + eps, alpha=0.08, color='blue')
ax.set_xlim(-1, 3)
ax.set_ylim(-3, 3)
ax.set_xlabel('$x$')
ax.set_ylabel('$f(x)$')
ax.set_title(r"Well-conditioned: $|f'(x^*)| = 2$, $\hat{\kappa} = 0.5$"
             "\nsteep crossing, small horizontal shift")
ax.grid(True, alpha=0.3)

# --- Right: Ill-conditioned (shallow crossing) ---
ax = axes[1]
f_ill = lambda x: 0.3*(x - 1)
ax.plot(x, f_ill(x), 'b-', linewidth=2)
ax.plot(x, f_ill(x) + eps, 'b--', linewidth=1.5, alpha=0.6)
ax.plot(x, f_ill(x) - eps, 'b--', linewidth=1.5, alpha=0.6)
ax.axhline(0, color='k', linewidth=0.8)

# Original root
ax.plot(root, 0, 'ko', markersize=8, zorder=5)

# Perturbed roots
root_plus = 1 - eps/0.3
root_minus = 1 + eps/0.3
ax.plot(root_plus, 0, 'rs', markersize=7, zorder=5)
ax.plot(root_minus, 0, 'rs', markersize=7, zorder=5)

# Show horizontal shift
ax.annotate('', xy=(root_minus, -0.4), xytext=(root, -0.4),
            arrowprops=dict(arrowstyle='<->', color='red', lw=1.5))
ax.text((root + root_minus)/2, -0.7, r'$\Delta x^*$', ha='center',
        fontsize=11, color='red')

# Show vertical perturbation
ax.annotate('', xy=(-0.3, eps), xytext=(-0.3, 0),
            arrowprops=dict(arrowstyle='<->', color='gray', lw=1.2))
ax.text(-0.55, eps/2, r'$\epsilon$', ha='center', fontsize=11, color='gray')

ax.fill_between(x, f_ill(x) - eps, f_ill(x) + eps, alpha=0.08, color='blue')
ax.set_xlim(-1, 3)
ax.set_ylim(-3, 3)
ax.set_xlabel('$x$')
ax.set_ylabel('$f(x)$')
ax.set_title(r"Ill-conditioned: $|f'(x^*)| = 0.3$, $\hat{\kappa} \approx 3.3$"
             "\nshallow crossing, large horizontal shift")
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()
```

(fig-condition-number-geometry)=
**Conditioning of a root.** Both panels show the same vertical perturbation $\epsilon$ (gray) applied to $f$ (blue band between dashed curves). *Left:* When $f$ crosses zero steeply ($|f'(x^*)| = 2$), the perturbed roots (red squares) barely move from the true root (black dot). *Right:* When $f$ crosses zero at a shallow angle ($|f'(x^*)| = 0.3$), the same perturbation slides the root much further along the $x$-axis. The horizontal shift $\Delta x^*$ is amplified by the condition number $\hat{\kappa} = 1/|f'(x^*)|$.

### Examples

::::{prf:example} Well-Conditioned Root
:label: ex-well-conditioned-root
:class: dropdown

Consider $f(x) = x - 1$ with root $x^* = 1$.

Since $f'(x) = 1$, we have $f'(x^*) = 1$, so:

$$
\kappa = \frac{1}{|f'(1)|} = 1
$$

This is perfectly conditioned. A perturbation of size $\epsilon$ in $f$ causes a perturbation of size $\epsilon$ in the root.
::::

::::{prf:example} Ill-Conditioned Root
:label: ex-ill-conditioned-root
:class: dropdown

Consider $f(x) = (x-1)^2 - 10^{-10}$ with roots near $x = 1$.

The roots are at $x^* = 1 \pm 10^{-5}$.

At these roots, $f'(x^*) = 2(x^* - 1) = \pm 2 \times 10^{-5}$, so:

$$
\kappa = \frac{1}{2 \times 10^{-5}} = 5 \times 10^{4}
$$

A perturbation of size $10^{-16}$ (machine epsilon) could move the root by $5 \times 10^{-12}$. We can expect to lose about 5 digits of accuracy!
::::

## Forward and Backward Error

The condition number tells us how sensitive the *problem* is, but in practice we
face another problem: our algorithm returns an approximation $\hat{x}$ and
we need to judge how good it is. The true root $x^*$ is unknown, so we cannot
directly compute $|\hat{x} - x^*|$. What *can* we compute? We can evaluate
$f(\hat{x})$, which tells us how well $\hat{x}$ satisfies the equation. But
does a small $f(\hat{x})$ guarantee that $\hat{x}$ is close to $x^*$?

The answer, as the condition number analysis suggests, is: *it depends on the
conditioning*. This leads to two different ways of measuring error. One measures
how far $\hat{x}$ is from the true answer. The other measures how well $\hat{x}$
satisfies the equation it is supposed to solve.

:::{prf:definition} Forward and Backward Error for Root Finding
:label: def-forward-backward-error-root-finding

Let $x^*$ be the true root of $f(x) = 0$ and let $\hat{x}$ be an approximation.

- The **forward error** is $|\hat{x} - x^*|$: how far $\hat{x}$ is from the
  true root.
- The **backward error** (also called the **residual**) is $|f(\hat{x})|$: how
  well $\hat{x}$ satisfies the equation.
:::

The forward error is what we care about but cannot compute (we don't know
$x^*$). The backward error is easy to compute, just evaluate $f(\hat{x})$.
As we will see, they are connected through the condition number.

:::{prf:theorem} Relating Forward and Backward Error
:label: thm-relating-forward-backward-error

Suppose $f \in \mathcal{C}^1$, $f(x^*) = 0$ is a simple root, and $\hat{x}$ is
sufficiently close to $x^*$. Then:

$$
\underbrace{|\hat{x} - x^*|}_{\text{forward error}} \leq \hat{\kappa} \cdot \underbrace{|f(\hat{x})|}_{\text{backward error}}
$$

where $\hat{\kappa} = 1/|f'(x^*)|$ is the condition number of the root.
:::

:::{prf:proof}
:class: dropdown

By the mean value theorem:

$$
f(\hat{x}) = f(\hat{x}) - f(x^*) = f'(\xi)(\hat{x} - x^*)
$$

for some $\xi$ in the interval $I = [x^*, \hat{x}]$ (or $[\hat{x}, x^*]$).
Define $m = \min_{x \in I} |f'(x)|$ and $M = \max_{x \in I} |f'(x)|$. Since
$\xi \in I$, we have $m \leq |f'(\xi)| \leq M$, so:

$$
m \, |\hat{x} - x^*| \leq |f(\hat{x})| \leq M \, |\hat{x} - x^*|
$$

Rearranging gives bounds on the forward error in both directions:

$$
\frac{1}{M} \cdot |f(\hat{x})| \leq |\hat{x} - x^*| \leq \frac{1}{m} \cdot |f(\hat{x})|
$$

As $I$ shrinks around $x^*$, both $m$ and $M$ converge to $|f'(x^*)|$, so the
upper bound becomes $|\hat{x} - x^*| \leq |f(\hat{x})| / |f'(x^*)| = \hat{\kappa} \cdot |f(\hat{x})|$.
:::

This is a fundamental principle that appears throughout numerical analysis:

$$
\text{forward error} \leq \text{condition number} \times \text{backward error}
$$

A small residual guarantees a small forward error *only when the problem is
well-conditioned*. We see it here for root finding with
$\hat{\kappa} = 1/|f'(x^*)|$. We will see it again for
[linear systems](../direct-methods/index.md) in the form
$\|\Delta\mathbf{x}\|/\|\mathbf{x}\| \leq \kappa(\mathbf{A}) \cdot \|\Delta\mathbf{b}\|/\|\mathbf{b}\|$.
The same structure governs eigenvalue problems, least squares, differential
equations, and essentially every problem in scientific computing.

```{code-cell} python
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

residual = 0.5
x_star = 0.0
xlim = (-1, 4)
ylim = (-1, 2)

for ax, slope, title in zip(axes,
    [2.0, 0.2],
    [r"$|f'(x^*)| = 2$, $\hat{\kappa} = 0.5$",
     r"$|f'(x^*)| = 0.2$, $\hat{\kappa} = 5$"]):

    x = np.linspace(xlim[0], xlim[1], 400)
    ax.plot(x, slope * x, 'b-', linewidth=2.5)
    ax.axhline(0, color='k', linewidth=0.8)

    x_hat = residual / slope

    # Backward error (green): vertical line from (x_hat, 0) to (x_hat, residual)
    ax.plot([x_hat, x_hat], [0, residual], 'g-', linewidth=4, solid_capstyle='round', zorder=4)

    # Forward error (red): horizontal line from x* to x_hat on the x-axis
    ax.plot([x_star, x_hat], [0, 0], 'r-', linewidth=4, solid_capstyle='round', zorder=3)

    # Mark points
    ax.plot(x_star, 0, 'ko', markersize=9, zorder=5, label=r'true root $x^*$')
    ax.plot(x_hat, 0, 'rs', markersize=8, zorder=5, label=r'approximation $\hat{x}$')
    ax.legend(fontsize=9, loc='upper left')

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel('$x$', fontsize=12)
    ax.set_ylabel('$f(x)$', fontsize=12)
    ax.set_title(title, fontsize=11)
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()
```

(fig-forward-backward-error)=
**Forward vs backward error.** The green arrow is the backward error (residual $|f(\hat{x})|$), the red arrow is the forward error ($|\hat{x} - x^*|$). Both panels have the same residual of $0.5$. *Left:* The steep crossing has a small condition number ($\hat{\kappa} = 0.5$), so the forward error is $0.25$, smaller than the residual. *Right:* The shallower crossing has a larger condition number ($\hat{\kappa} = 5$), and the same residual now gives a forward error of $2.5$, ten times larger. A larger condition number means more amplification from backward to forward error.

### Stopping Criteria

This relationship has immediate consequences for when to stop iterating. So far,
our algorithms for fixed point iteration and Newton's method included a stopping
test $|x_{n+1} - x_n| < \varepsilon$, but we never justified *why* a small step
size means we are close to the answer. The step size $|x_{n+1} - x_n|$ is not
the forward error $|x_n - x^*|$. So why is it a useful proxy?

:::{prf:lemma} Step Size Bounds the Forward Error
:label: lem-step-size-error-bound

Suppose an iterative method produces a sequence $\{x_n\}$ that is contractive:

$$
|x_{n+1} - x_n| \leq \theta |x_n - x_{n-1}| \quad \text{for some } \theta < 1
$$

Then the sequence converges to some limit $x^*$, and the forward error satisfies:

$$
|x_n - x^*| \leq \frac{|x_{n+1} - x_n|}{1 - \theta}
$$
:::

:::{prf:proof}
:class: dropdown

Since the sequence converges, we can write the error as a telescoping sum:

$$
x_n - x^* = x_n - \lim_{k \to \infty} x_k = \sum_{k=0}^{\infty} (x_{n+k} - x_{n+k+1})
$$

Taking absolute values and using the contraction property
$|x_{n+k} - x_{n+k+1}| \leq \theta^k |x_n - x_{n+1}|$:

$$
|x_n - x^*| \leq \sum_{k=0}^{\infty} \theta^k |x_n - x_{n+1}| = \frac{|x_{n+1} - x_n|}{1 - \theta}
$$
:::

This is a consequence of the
[Banach fixed point theorem](fixed-point.md#thm-banach) and applies to any
contractive iteration (fixed point, Jacobi, Gauss-Seidel, Newton near the
root). The step test $|x_{n+1} - x_n| < \varepsilon$ therefore guarantees
$|x_n - x^*| < \varepsilon / (1 - \theta)$. When $\theta$ is small (fast
convergence), this is close to $\varepsilon$. When $\theta$ is close to 1 (slow
convergence), the factor $1/(1-\theta)$ can be large and the step test becomes
unreliable.

:::{prf:remark} Two Stopping Criteria
:label: rmk-stopping-criteria
:class: dropdown

In practice, we use **two** stopping criteria:

1. **Residual test:** Stop when $|f(x_n)| < \varepsilon_1$ (small backward error)
2. **Step test:** Stop when $|x_{n+1} - x_n| < \varepsilon_2$ (approximate forward error)

For **well-conditioned** problems ($\hat{\kappa} \approx 1$), these are
equivalent: a small residual guarantees a small forward error.

For **ill-conditioned** problems ($\hat{\kappa} \gg 1$), they can disagree
dramatically. A small residual $|f(\hat{x})| = 10^{-15}$ with $\hat{\kappa} =
10^5$ means the forward error could be as large as $10^{-10}$.

**Rule of thumb:** always check both. If the residual is small but the iterates
are still moving, the problem is ill-conditioned.
:::

## Repeated Roots

We now turn to a class of problems where ill-conditioning is intrinsic.

A root $x^*$ is called **repeated** (or a root of higher multiplicity) when $f$
not only vanishes at $x^*$ but also touches zero "flatly," meaning its first
several derivatives vanish there as well.

:::{prf:definition} Root Multiplicity
:label: def-root-multiplicity

A root $x^*$ of $f$ has **multiplicity** $m$ if:

$$
f(x) = (x - x^*)^m h(x), \quad h(x^*) \neq 0
$$

Equivalently, $f(x^*) = f'(x^*) = \cdots = f^{(m-1)}(x^*) = 0$ but $f^{(m)}(x^*) \neq 0$.

A root with $m = 1$ is called **simple**. A root with $m \geq 2$ is called
**repeated**.
:::

```{code-cell} python
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(8, 4))
x = np.linspace(-1.5, 1.5, 400)

for m, color, label in [(1, '#1f77b4', r'$f(x) = x$ (simple, $m=1$)'),
                         (2, '#ff7f0e', r'$f(x) = x^2$ (double, $m=2$)'),
                         (4, '#2ca02c', r'$f(x) = x^4$ (quadruple, $m=4$)')]:
    ax.plot(x, x**m, color=color, linewidth=2.5, label=label)

ax.axhline(0, color='k', linewidth=0.8)
ax.plot(0, 0, 'ko', markersize=8, zorder=5)

ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-0.5, 1.5)
ax.set_xlabel('$x$', fontsize=12)
ax.set_ylabel('$f(x)$', fontsize=12)
ax.set_title('Simple vs repeated roots: higher multiplicity means flatter near the root')
ax.legend(fontsize=9, loc='upper center')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()
```

(fig-repeated-roots-intro)=
**Simple vs repeated roots.** A simple root ($m=1$, blue) crosses zero at a nonzero angle. A double root ($m=2$, orange) merely touches zero and turns around. A quadruple root ($m=4$, green) is flatter still. The higher the multiplicity, the flatter $f$ is near the root, which is precisely what makes repeated roots ill-conditioned.

The distinction matters for two reasons: conditioning and convergence speed.

### Conditioning of repeated roots

For a simple root, the condition number is $\hat{\kappa} = 1/|f'(x^*)|$, which
is finite. For a repeated root, $f'(x^*) = 0$, so the formula gives
$\hat{\kappa} = \infty$. But what does "infinite condition number" actually mean
in practice? It does *not* mean the root is impossible to compute. It means
that the relationship between perturbation size and root shift is no longer
linear.

:::{prf:theorem} Sensitivity of Repeated Roots
:label: thm-ill-conditioning-repeated-roots

If $x^*$ is a root of multiplicity $m$, then a perturbation $f \to f + \epsilon$
shifts the root by $\mathcal{O}(\epsilon^{1/m})$.
:::

:::{prf:proof}
:class: dropdown

Near a root of multiplicity $m$, we have $f(x) \approx c(x - x^*)^m$ for some
constant $c = f^{(m)}(x^*)/m! \neq 0$. The perturbed equation
$f(\tilde{x}^*) + \epsilon = 0$ gives:

$$
c(\tilde{x}^* - x^*)^m \approx -\epsilon \implies \tilde{x}^* - x^* \approx \left(\frac{\epsilon}{|c|}\right)^{1/m}
$$
:::

The exponent $1/m$ is the key. When $\epsilon$ is small (as in floating-point
rounding), $\epsilon^{1/m} \gg \epsilon$. For example, with $\epsilon = 10^{-16}$
(machine epsilon), a simple root shifts by $10^{-16}$, a double root shifts by
$10^{-8}$, and a triple root shifts by $10^{-16/3} \approx 10^{-5.3}$.
The following figure shows *why* this happens geometrically.

```{code-cell} python
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

x_star = 1.0
eps = 0.2
x = np.linspace(0.2, 1.8, 500)

titles = [r'Simple root ($m=1$): shift $= \epsilon$',
          r'Double root ($m=2$): shift $= \sqrt{\epsilon}$',
          r'Triple root ($m=3$): shift $= \epsilon^{1/3}$']

for ax, m, title in zip(axes, [1, 2, 3], titles):
    f_vals = (x - x_star) ** m
    ax.plot(x, f_vals, 'b-', linewidth=2.5, label=r'$f(x)$')
    # For m=1, shift up; for m>=2 (non-negative near root), shift down
    sign = 1 if m == 1 else -1
    ax.plot(x, f_vals + sign * eps, 'b--', linewidth=2, alpha=0.6,
            label=rf'$f(x) {"+" if sign > 0 else "-"} \epsilon$')
    ax.axhline(0, color='k', linewidth=0.8)

    # Perturbed root
    shift = eps ** (1.0 / m)
    x_tilde = x_star - shift if m == 1 else x_star + shift

    # Forward error (red): horizontal line showing the root shift
    ax.plot([x_tilde, x_star], [0, 0], 'r-', linewidth=3.5,
            solid_capstyle='round', zorder=3)

    # Mark points
    ax.plot(x_star, 0, 'ko', markersize=9, zorder=5)
    ax.plot(x_tilde, 0, 'rs', markersize=8, zorder=5)

    ax.set_xlim(0.2, 1.8)
    ax.set_ylim(-0.15, 0.25)
    ax.set_xlabel('$x$', fontsize=11)
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=8, loc='upper left')
    ax.grid(True, alpha=0.3)

axes[0].set_ylabel('$f(x)$', fontsize=11)
plt.tight_layout()
plt.show()
```

(fig-repeated-root-scaling)=
**Root shift under perturbation for $f(x) = (x - x^*)^m$.** The solid blue curve
is the original function, the dashed blue curve is the perturbed function. The
black dot is the root of $f$, the red square is the root of the perturbed
function, and the red line is the root shift (forward error). *Left:* A simple
root crosses zero steeply, so the perturbation barely moves the root. *Center:*
A double root is flat near $x^*$, so the perturbed curve's zero crossing slides
much further: the shift is $\sqrt{\epsilon}$. *Right:* A triple root is flatter
still, and the shift grows to $\epsilon^{1/3}$. Note that for even-multiplicity
roots ($m = 2$), we must shift *down* to create a new root. A perturbation in
the other direction would lift the curve away from zero entirely, destroying the
root. This is another way even-multiplicity roots are fragile: they can simply
disappear under perturbation.

### Convergence of Newton's method at repeated roots

Beyond the conditioning issue, Newton's method also *converges more slowly* at
repeated roots.

:::{prf:proposition} Reduced Convergence at Repeated Roots
:label: prop-repeated-roots-convergence

When $x^*$ is a root of multiplicity $m \geq 2$, Newton's method converges only
**linearly** with rate $\rho = 1 - 1/m$.
:::

:::{prf:proof}
:class: dropdown

Write $f(x) = (x - x^*)^m h(x)$ with $h(x^*) \neq 0$. Newton's iteration
function is $g(x) = x - f(x)/f'(x)$. A calculation shows:

$$
g'(x^*) = 1 - \frac{1}{m}
$$

Since $g'(x^*) \neq 0$, convergence is linear with rate
$\rho = |g'(x^*)| = 1 - 1/m$.

For a double root ($m = 2$), $\rho = 1/2$. For a triple root ($m = 3$),
$\rho = 2/3$. As the multiplicity increases, convergence slows.
:::

So repeated roots cause a double penalty: fewer accurate digits *and* slower
convergence. Both problems can be mitigated if the multiplicity is known.

:::{prf:remark} Modified Newton for Repeated Roots
:label: rmk-modified-newton
:class: dropdown

If the multiplicity $m$ is known, the modified iteration:

$$
x_{n+1} = x_n - m\frac{f(x_n)}{f'(x_n)}
$$

restores quadratic convergence. This works because the modified iteration
function $g(x) = x - mf(x)/f'(x)$ satisfies $g'(x^*) = 0$.

In practice, $m$ is rarely known in advance. An alternative is to apply Newton's
method to $u(x) = f(x)/f'(x)$, which has a simple root at $x^*$ regardless of
the multiplicity of $f$. This avoids the need to know $m$ but requires
computing both $f$ and $f'$ at each step.
:::
