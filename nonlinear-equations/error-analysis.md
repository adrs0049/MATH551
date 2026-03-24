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

:::{figure}
:label: fig-condition-number-geometry

**Conditioning of a root.** Both panels show the same vertical perturbation $\epsilon$ (gray) applied to $f$ (blue band between dashed curves). *Left:* When $f$ crosses zero steeply ($|f'(x^*)| = 2$), the perturbed roots (red squares) barely move from the true root (black dot). *Right:* When $f$ crosses zero at a shallow angle ($|f'(x^*)| = 0.3$), the same perturbation slides the root much further along the $x$-axis. The horizontal shift $\Delta x^*$ is amplified by the condition number $\hat{\kappa} = 1/|f'(x^*)|$.
:::

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

## Repeated Roots

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
rounding), $\epsilon^{1/m} \gg \epsilon$:

| Multiplicity | Root shift from $\epsilon = 10^{-16}$ | Digits lost |
|:---:|:---:|:---:|
| $m = 1$ (simple) | $\sim 10^{-16}$ | 0 |
| $m = 2$ (double) | $\sim 10^{-8}$ | 8 |
| $m = 3$ (triple) | $\sim 10^{-5.3}$ | 11 |

So if $f$ itself carries uncertainty of size $\epsilon$ (e.g., from measured
coefficients or upstream rounding), a double root can be located to at most
about 8 digits in double precision, and a triple root to about 5. This is a
limitation of the *problem*, not the algorithm.

However, if $f$ is given exactly (e.g., a polynomial with known coefficients),
the situation is better. Methods like modified Newton (see below) can converge
to a repeated root to nearly full precision, because the function is not
perturbed; only the floating-point arithmetic introduces errors at each step.

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

## Forward and Backward Error

The condition number tells us how sensitive the *problem* is, but in practice we
face a different question: our algorithm returns an approximation $\hat{x}$ and
we need to judge how good it is. The true root $x^*$ is unknown, so we cannot
directly compute $|\hat{x} - x^*|$. What *can* we compute? We can evaluate
$f(\hat{x})$, which tells us how well $\hat{x}$ satisfies the equation. But
does a small $f(\hat{x})$ guarantee that $\hat{x}$ is close to $x^*$?

The answer, as the condition number analysis suggests, is: *it depends on the
conditioning*. This leads to two complementary ways of measuring error.

:::{prf:definition} Forward and Backward Error for Root Finding
:label: def-forward-backward-error-root-finding

Given an approximate root $\hat{x}$ of $f(x) = 0$:

- The **forward error** is $|\hat{x} - x^*|$: how far $\hat{x}$ is from the
  true root.
- The **backward error** is $|f(\hat{x})|$: how well $\hat{x}$ satisfies the
  equation.
:::

The forward error is what we care about but cannot compute (we don't know
$x^*$). The backward error is easy to compute, just evaluate $f(\hat{x})$.
The condition number connects them.

### Relating Forward and Backward Error

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

This is the fundamental relationship of numerical root finding: the forward
error is bounded by the condition number times the backward error. A small
residual guarantees a small forward error *only when the problem is
well-conditioned* ($\hat{\kappa}$ small). This relationship reappears throughout
scientific computing: in [linear systems](../direct-methods/index.md), the
analogous result is $\|\Delta\mathbf{x}\|/\|\mathbf{x}\| \leq \kappa(\mathbf{A}) \cdot \|\Delta\mathbf{b}\|/\|\mathbf{b}\|$,
and in [globalization strategies](../nonlinear-systems/globalization.md), the
choice of *what to monitor* (residual vs. correction) determines the robustness
of the solver.

### Stopping Criteria

This relationship has immediate consequences for when to stop iterating. But
first, we need to understand why the step size $|x_{n+1} - x_n|$ is a useful
proxy for the forward error.

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

