# Second-Order Finite Differences

:::{tip} Big Idea
By using information on both sides of a point, the central difference formula achieves second-order accuracy ($\mathcal{O}(h^2)$) compared to the first-order accuracy ($\mathcal{O}(h)$) of the forward difference.
:::

## Motivation

In the previous section, we developed a forward-difference approximation with error proportional to $h$. Can we do better? That is, can we develop a finite difference approximation with a **smaller error**?

From geometric intuition (draw a picture!), using points on both sides of $x_0$ should give a better approximation to the tangent line.

## The Central Difference Formula

:::{prf:definition} Central Difference Approximation
:label: def-central-difference

$$
f'(x_0) \approx \frac{f(x_0 + h) - f(x_0 - h)}{2h}
$$
:::

## Error Analysis

We replace both $f(x_0 \pm h)$ using Taylor series with remainder:

$$
f(x_0 \pm h) = f(x_0) \pm h f'(x_0) + \frac{h^2}{2}f''(x_0) \pm \frac{h^3}{6} f'''(\xi^{\pm})
$$

Substituting into the central difference formula:

$$
\begin{split}
\frac{f(x_0 + h) - f(x_0 - h)}{2h}
&= \frac{1}{2h}\bigg(
    f(x_0) + h f'(x_0) + \frac{h^2}{2}f''(x_0) + \frac{h^3}{6} f'''(\xi^+) \\
    &\qquad\qquad - f(x_0) + h f'(x_0) - \frac{h^2}{2}f''(x_0) + \frac{h^3}{6} f'''(\xi^-) \bigg) \\
&= \frac{1}{2h} \left( 2hf'(x_0) + \frac{h^3}{6}\left( f'''(\xi^+) + f'''(\xi^-) \right)\right)
\end{split}
$$

Thus:

$$
\frac{f(x_0 + h) - f(x_0 - h)}{2h} = f'(x_0) + \frac{h^2}{12} \left( f^{(3)}(\xi^-) + f^{(3)}(\xi^+) \right)
$$

:::{prf:proposition} Central Difference Error
:label: prop-central-difference-error

The error of the central difference approximation is:

$$
E(h) = \frac{h^2}{12} \left( f^{(3)}(\xi^-) + f^{(3)}(\xi^+) \right) = \mathcal{O}(h^2)
$$

The error scales with the **square** of the discretization size, making it smaller than the forward-difference scheme for the same $h$.
:::

## Comparison: Forward vs Central Difference

| Method | Formula | Error | Order |
|--------|---------|-------|-------|
| Forward Difference | $\frac{f(x_0 + h) - f(x_0)}{h}$ | $\mathcal{O}(h)$ | 1st order |
| Central Difference | $\frac{f(x_0 + h) - f(x_0 - h)}{2h}$ | $\mathcal{O}(h^2)$ | 2nd order |

:::{prf:example} Comparing Forward and Central Difference Errors
:label: ex-fd-comparison

**Comparing errors for $h = 0.01$:**

- Forward difference error: $\sim 0.01$
- Central difference error: $\sim 0.0001$

The central difference is **100 times more accurate** for the same step size!
:::

## Computational Exercise

:::{admonition} Exercise
:class: seealso

1. Implement the central difference scheme on a computer.

2. Plot a log-log plot of the error vs. the step size.

3. Verify that the slope of the line is approximately **2** (indicating second-order convergence).
:::

## Optimal Step Size for Central Differences

Following the same analysis as for forward differences, but now with second-order approximation error:

$$
E(h) = \frac{M_0 \eta}{h} + M_1 h^2
$$

The optimal step size is:

$$
h_{\text{opt}} \approx \eta^{1/3} \approx 6 \times 10^{-6}
$$

for double precision arithmetic.

:::{note}
**Practical takeaway:** For central differences with double precision, use $h \approx 10^{-5}$ to $10^{-6}$. This is larger than the optimal $h$ for forward differences because the approximation error decreases faster.
:::

## Summary

The central difference formula demonstrates a key principle in numerical analysis: **using more information can improve accuracy**. By sampling the function on both sides of the point of interest, we achieve quadratic rather than linear convergence.

This idea generalizes: higher-order finite difference formulas use more points and achieve even higher accuracy, though with diminishing returns due to round-off error and increased computational cost.
