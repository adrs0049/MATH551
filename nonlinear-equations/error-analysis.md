# Error Analysis for Root Finding

:::{tip} Big Idea
Understanding *how* errors behave in root finding connects three ideas: **condition numbers** tell us how sensitive the root is to perturbations, **forward/backward error** gives us a framework for interpreting computed results, and **convergence rates** tell us how quickly our iterates approach the root.
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

### Multiple Roots: A Warning

::::::{prf:remark} Ill-Conditioning at Multiple Roots
:label: rmk-ill-conditioning-multiple-roots

If $f(x^*) = 0$ and $f'(x^*) = 0$ (a **multiple root**), the condition number is infinite. Multiple roots are inherently ill-conditioned.

::::{dropdown} Why?
At a double root, $f(x) \approx c(x - x^*)^2$ near $x^*$. A small perturbation $\epsilon$ to the function can:
- Split the double root into two simple roots
- Eliminate the root entirely
- Move the root by $\mathcal{O}(\sqrt{\epsilon})$ instead of $\mathcal{O}(\epsilon)$

This $\sqrt{\epsilon}$ sensitivity means the condition number is effectively $1/\sqrt{\epsilon} \to \infty$.
::::
::::::

### Examples

::::{admonition} Example: Well-Conditioned Root
:class: note

Consider $f(x) = x - 1$ with root $x^* = 1$.

:::{dropdown} Analysis
Since $f'(x) = 1$, we have $f'(x^*) = 1$, so:

$$
\kappa = \frac{1}{|f'(1)|} = 1
$$

This is perfectly conditioned. A perturbation of size $\epsilon$ in $f$ causes a perturbation of size $\epsilon$ in the root.
:::
::::

::::{admonition} Example: Ill-Conditioned Root
:class: note

Consider $f(x) = (x-1)^2 - 10^{-10}$ with roots near $x = 1$.

:::{dropdown} Analysis
The roots are at $x^* = 1 \pm 10^{-5}$.

At these roots, $f'(x^*) = 2(x^* - 1) = \pm 2 \times 10^{-5}$, so:

$$
\kappa = \frac{1}{2 \times 10^{-5}} = 5 \times 10^{4}
$$

A perturbation of size $10^{-16}$ (machine epsilon) could move the root by $5 \times 10^{-12}$. We can expect to lose about 5 digits of accuracy!
:::
::::

## Forward and Backward Error

The [forward/backward error framework](../qr-least-squares/forward-backward-error.md) from earlier applies naturally to root finding.

### Forward Error

:::{prf:definition} Forward Error for Root Finding
:label: def-forward-error-root-finding

Given an approximate root $\hat{x}$, the **forward error** is:

$$
\text{Forward Error} = |\hat{x} - x^*|
$$

This measures how far our approximation is from the true root.
:::

The forward error is what we ultimately care about, but it's hard to compute directly (we don't know $x^*$!).

### Backward Error

:::{prf:definition} Backward Error for Root Finding
:label: def-backward-error-root-finding

Given an approximate root $\hat{x}$, the **backward error** is:

$$
\text{Backward Error} = |f(\hat{x})|
$$

This measures how well $\hat{x}$ satisfies the equation $f(x) = 0$.
:::

The backward error is easy to compute: just evaluate $f(\hat{x})$.

### The Connection

::::::{prf:theorem} Relating Forward and Backward Error
:label: thm-relating-forward-backward-error

For root finding, the forward and backward errors are related by:

$$
|\hat{x} - x^*| \approx \frac{|f(\hat{x})|}{|f'(x^*)|} = \kappa \cdot |f(\hat{x})|
$$

::::{dropdown} Derivation
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

### Practical Implications

:::{admonition} Stopping Criteria
:class: note

This relationship explains why we use **two stopping criteria** in practice:

1. **Residual test:** Stop when $|f(x_n)| < \epsilon_1$ (small backward error)
2. **Step test:** Stop when $|x_{n+1} - x_n| < \epsilon_2$ (approximate forward error)

For well-conditioned problems ($\kappa \approx 1$), these are equivalent. For ill-conditioned problems, a small residual doesn't guarantee a small forward error!
:::

## Convergence Rates

Different root-finding methods converge at different rates. Understanding these rates helps us choose the right method.

### Linear Convergence

::::::{prf:definition} Linear (First-Order) Convergence
:label: def-linear-convergence

A sequence $\{x_n\}$ converges **linearly** to $x^*$ if:

$$
|x_{n+1} - x^*| \leq L \, |x_n - x^*|
$$

for some constant $L < 1$.

::::{dropdown} What this means
Each iteration reduces the error by a constant factor $L$. If $L = 0.5$, we gain about 0.3 digits per iteration ($-\log_{10}(0.5) \approx 0.3$).

After $n$ iterations:
$$
|x_n - x^*| \leq L^n |x_0 - x^*|
$$

To gain $d$ digits of accuracy, we need approximately $n \approx d / \log_{10}(1/L)$ iterations.
::::
::::::

**Examples of linear convergence:**
- **Bisection:** $L = 0.5$ (exactly halves error each step)
- **Fixed-point iteration:** $L = |g'(x^*)|$ when $|g'(x^*)| < 1$

### Quadratic Convergence

::::::{prf:definition} Quadratic (Second-Order) Convergence
:label: def-quadratic-convergence

A sequence $\{x_n\}$ converges **quadratically** to $x^*$ if:

$$
|x_{n+1} - x^*| \leq C \, |x_n - x^*|^2
$$

for some constant $C$.

::::{dropdown} What this means
The number of correct digits approximately **doubles** each iteration!

If the error is $10^{-3}$, the next error is roughly $10^{-6}$, then $10^{-12}$.

To go from 3 correct digits to 12 correct digits takes only 2 iterations.
::::
::::::

**Examples of quadratic convergence:**
- **Newton's method** (for simple roots with good initial guess)

### Superlinear Convergence

::::::{prf:definition} Superlinear Convergence
:label: def-superlinear-convergence

A sequence $\{x_n\}$ converges **superlinearly** to $x^*$ if:

$$
\lim_{n \to \infty} \frac{|x_{n+1} - x^*|}{|x_n - x^*|} = 0
$$

::::{dropdown} What this means
The ratio of consecutive errors goes to zero—faster than any geometric rate, but not necessarily quadratic.

Superlinear convergence is "between" linear and quadratic: the effective contraction factor improves as we get closer to the root.
::::
::::::

**Examples of superlinear convergence:**
- **Secant method:** Order $\approx 1.618$ (the golden ratio)
- **Quasi-Newton methods:** Achieve superlinear without computing exact derivatives

### General Order of Convergence

::::::{prf:definition} Order of Convergence
:label: def-order-convergence

A sequence $\{x_n\}$ converging to $x^*$ has **order of convergence** $p \geq 1$ if:

$$
\lim_{n \to \infty} \frac{|x_{n+1} - x^*|}{|x_n - x^*|^p} = C
$$

for some constant $C > 0$ (the **asymptotic error constant**).

| Order $p$ | Name | Digits gained per iteration |
|-----------|------|----------------------------|
| $p = 1$ | Linear | $-\log_{10}(C)$ (constant) |
| $1 < p < 2$ | Superlinear | Increasing |
| $p = 2$ | Quadratic | Doubles |
| $p = 3$ | Cubic | Triples |
::::::

::::::{prf:example} Orders of Common Methods
:label: ex-convergence-orders

| Method | Order $p$ | Asymptotic Error Constant |
|--------|-----------|---------------------------|
| Bisection | 1 | $C = 0.5$ |
| Fixed-point ($g'(x^*) \neq 0$) | 1 | $C = \|g'(x^*)\|$ |
| Secant | $\phi \approx 1.618$ | Depends on $f$ |
| Newton (simple root) | 2 | $C = \|f''(x^*)\|/(2\|f'(x^*)\|)$ |
| Newton (double root) | 1 | Degrades to linear! |
| Halley | 3 | Uses second derivatives |

:::{dropdown} Why the secant method has order $\phi$
The secant method uses $x_{n+1} = x_n - f(x_n) \frac{x_n - x_{n-1}}{f(x_n) - f(x_{n-1})}$.

The error satisfies $e_{n+1} \approx C \cdot e_n \cdot e_{n-1}$. Assuming $e_n \sim e_0^{p^n}$ leads to $p^2 = p + 1$, giving $p = \frac{1 + \sqrt{5}}{2} = \phi$.
:::
::::::

### Comparison

| Method | Convergence | Rate | Iterations for 10 digits |
|--------|-------------|------|--------------------------|
| Bisection | Linear | $L = 0.5$ | ~34 |
| Fixed-point (typical) | Linear | $L = 0.9$ | ~220 |
| Newton's method | Quadratic | doubles digits | ~4-5 |

:::{admonition} The Trade-off
:class: tip

Faster convergence comes at a cost:
- **Bisection:** Slow but guaranteed (only needs $f$ continuous, sign change)
- **Fixed-point:** Moderate, needs $|g'| < 1$
- **Newton:** Fast but needs $f'$, good initial guess, and $f'(x^*) \neq 0$

In practice, a common strategy is: use bisection to get close, then switch to Newton for fast final convergence.
:::

## Summary

| Concept | Formula | Interpretation |
|---------|---------|----------------|
| Condition number | $\kappa = 1/\|f'(x^*)\|$ | Sensitivity of root to perturbations |
| Forward error | $\|\hat{x} - x^*\|$ | Distance from true root |
| Backward error | $\|f(\hat{x})\|$ | Residual |
| Error relationship | Forward $\approx \kappa \times$ Backward | Why both stopping criteria matter |
| Linear convergence | $e_{n+1} \leq L \cdot e_n$ | Constant factor reduction |
| Superlinear convergence | $e_{n+1}/e_n \to 0$ | Ratio improves each step |
| Quadratic convergence | $e_{n+1} \leq C \cdot e_n^2$ | Digits double each step |
| Order $p$ convergence | $e_{n+1} \sim C \cdot e_n^p$ | General framework |
