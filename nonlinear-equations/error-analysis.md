---
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

## Condition Number for Roots

The condition number measures how sensitive the root $x^*$ is to perturbations in the function $f$.

:::{admonition} Setup
:class: note

Suppose $f(x^*) = 0$ and we perturb $f$ to $\tilde{f} = f + \epsilon g$. The perturbed root $\tilde{x}^*$ satisfies $\tilde{f}(\tilde{x}^*) = 0$.
:::

::::::{prf:theorem} Condition Number for a Simple Root
:label: thm-condition-number-simple-root

If $f(x^*) = 0$ and $f'(x^*) \neq 0$ (a **simple root**), then the absolute condition number is:

$$
\kappa = \frac{1}{|f'(x^*)|}
$$

The root is **well-conditioned** when $|f'(x^*)|$ is large (steep crossing) and **ill-conditioned** when $|f'(x^*)|$ is small (tangent crossing).

::::{dropdown} Derivation
Starting from $\tilde{f}(\tilde{x}^*) = 0$:

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

The sensitivity to the perturbation is $1/|f'(x^*)|$.
::::
::::::

:::{prf:remark} Geometric Interpretation
:label: rmk-condition-number-geometry

The condition number has a simple geometric meaning:

- **Well-conditioned** ($\kappa$ small): $f$ crosses zero steeply. A small vertical perturbation of the curve barely moves the root horizontally.
- **Ill-conditioned** ($\kappa$ large): $f$ crosses zero at a shallow angle. The same small vertical perturbation slides the root far along the $x$-axis.

This is why nearly-tangent crossings are dangerous: the root sits on a nearly flat part of $f$, so it's easy to perturb.
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

## Multiple Roots

Multiple roots deserve special attention because they are **both ill-conditioned and slow to converge** — a double penalty.

:::{prf:definition} Root Multiplicity
:label: def-root-multiplicity

A root $x^*$ of $f$ has **multiplicity** $m$ if:

$$
f(x) = (x - x^*)^m h(x), \quad h(x^*) \neq 0
$$

Equivalently, $f(x^*) = f'(x^*) = \cdots = f^{(m-1)}(x^*) = 0$ but $f^{(m)}(x^*) \neq 0$.
:::

::::::{prf:theorem} Ill-Conditioning at Multiple Roots
:label: thm-ill-conditioning-multiple-roots

If $x^*$ is a root of multiplicity $m \geq 2$, the condition number is **infinite**: the root is infinitely sensitive to perturbations.

::::{prf:proof}
:class: dropdown

At a multiple root, $f'(x^*) = 0$, so $\kappa = 1/|f'(x^*)| = \infty$.

More precisely, near a double root $f(x) \approx c(x - x^*)^2$. A perturbation $f \to f + \epsilon$ shifts the root by:

$$
c(\tilde{x}^* - x^*)^2 \approx -\epsilon \implies \tilde{x}^* - x^* \approx \pm\sqrt{\epsilon/c}
$$

The root moves by $\mathcal{O}(\sqrt{\epsilon})$ instead of $\mathcal{O}(\epsilon)$ — much worse. For a root of multiplicity $m$, the sensitivity is $\mathcal{O}(\epsilon^{1/m})$.
::::
::::::

:::{prf:proposition} Reduced Convergence for Multiple Roots
:label: prop-multiple-roots

When $x^*$ is a root of multiplicity $m \geq 2$, Newton's method converges only **linearly** with rate $\rho = 1 - 1/m$.
:::

:::{prf:proof}
:class: dropdown

Write $f(x) = (x - x^*)^m h(x)$ with $h(x^*) \neq 0$. Newton's iteration function is $g(x) = x - f(x)/f'(x)$. A calculation shows:

$$
g'(x^*) = 1 - \frac{1}{m}
$$

Since $g'(x^*) \neq 0$, convergence is linear with rate $\rho = |g'(x^*)| = 1 - 1/m$.

For a double root ($m = 2$), $\rho = 1/2$. For a triple root ($m = 3$), $\rho = 2/3$. As the multiplicity increases, convergence slows dramatically.
:::

:::{prf:remark} Modified Newton for Multiple Roots
:label: rmk-modified-newton
:class: dropdown

If the multiplicity $m$ is known, the modified iteration:
$$
x_{n+1} = x_n - m\frac{f(x_n)}{f'(x_n)}
$$
restores quadratic convergence. This works because the modified iteration function $g(x) = x - mf(x)/f'(x)$ satisfies $g'(x^*) = 0$.

In practice, $m$ is rarely known in advance. An alternative is to apply Newton's method to $u(x) = f(x)/f'(x)$, which has a simple root at $x^*$ regardless of the multiplicity of $f$.
:::

## Forward and Backward Error

Root-finding algorithms return an approximation $\hat{x} \approx x^*$. There are two natural ways to measure the quality of this approximation.

:::{prf:definition} Forward Error
:label: def-forward-error-root-finding

The **forward error** is the distance from the computed root to the true root:

$$
\Delta x = |\hat{x} - x^*|
$$
:::

:::{prf:definition} Backward Error
:label: def-backward-error-root-finding

The **backward error** is the residual — how well $\hat{x}$ satisfies the equation:

$$
|f(\hat{x})|
$$

This answers: *for what perturbation of the problem is $\hat{x}$ the exact answer?*
:::

The forward error is what we care about but can't compute (we don't know $x^*$). The backward error is easy to compute — just evaluate $f(\hat{x})$.

### The Connection

::::::{prf:theorem} Forward Error $\approx$ Condition Number $\times$ Backward Error
:label: thm-relating-forward-backward-error

For root finding near a simple root:

$$
\underbrace{|\hat{x} - x^*|}_{\text{forward error}} \approx \underbrace{\frac{1}{|f'(x^*)|}}_{\kappa} \cdot \underbrace{|f(\hat{x})|}_{\text{backward error}}
$$

::::{prf:proof}
:class: dropdown

By Taylor's theorem:

$$
f(\hat{x}) = f(x^*) + f'(\xi)(\hat{x} - x^*) = f'(\xi)(\hat{x} - x^*)
$$

for some $\xi$ between $\hat{x}$ and $x^*$. If $\hat{x} \approx x^*$, then $f'(\xi) \approx f'(x^*)$, so:

$$
|\hat{x} - x^*| \approx \frac{|f(\hat{x})|}{|f'(x^*)|}
$$
::::
::::::

This is the fundamental relationship of numerical root finding, and it has immediate practical consequences.

### Stopping Criteria

:::{prf:remark} Two Stopping Criteria
:label: rmk-stopping-criteria

In practice, we use **two** stopping criteria:

1. **Residual test:** Stop when $|f(x_n)| < \varepsilon_1$ (small backward error)
2. **Step test:** Stop when $|x_{n+1} - x_n| < \varepsilon_2$ (approximate forward error)

For **well-conditioned** problems ($\kappa \approx 1$), these are equivalent: a small residual guarantees a small forward error.

For **ill-conditioned** problems ($\kappa \gg 1$), they can disagree dramatically. A small residual $|f(\hat{x})| = 10^{-15}$ with $\kappa = 10^5$ means the forward error could be as large as $10^{-10}$.

**Rule of thumb:** always check both. If the residual is small but the iterates are still moving, the problem is ill-conditioned.
:::

