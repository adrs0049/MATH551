# Numerical Differentiation

:::{tip} Big Idea
Derivatives can be approximated using function values at discrete points. Taylor's theorem tells us exactly how accurate these approximations are and reveals the fundamental trade-off between truncation error and round-off error.
:::

## The Forward Difference Formula

Suppose we have a function $f(x)$ and we want to compute its derivative at a point $x_0$ (i.e., the slope of the tangent line to the curve $y = f(x)$ at $x_0$).

Recall the limit definition of the derivative:

$$
f'(x_0) = \lim_{h\to 0} \frac{f(x_0 + h) - f(x_0)}{h}
$$

The left-hand side is the slope of the tangent line, while the right-hand side is the limit of the slope of the secant line connecting the points $(x_0, f(x_0))$ and $(x_0 + h, f(x_0 + h))$.

This suggests we can approximate the derivative using:

$$
f'(x_0) \approx \frac{f(x_0 + h) - f(x_0)}{h}
$$

for $h$ "small". A major question is: **how accurate is this approximation?**

### Error Analysis

:::{prf:proposition} Forward Difference Error
:label: prop-forward-diff

The forward difference approximation satisfies:

$$
\frac{f(x_0 + h) - f(x_0)}{h} = f'(x_0) + E(h)
$$

where the error is:

$$
E(h) = -\frac{h}{2} f''(\xi) = \mathcal{O}(h)
$$

for some $\xi \in (x_0, x_0 + h)$. The error is **first-order** in $h$.
:::

:::{prf:proof}
:class: dropdown

Using Taylor's theorem with $k = 1$:

$$
f(x_0 + h) = f(x_0) + f'(x_0) h + \frac{f''(\xi)}{2} h^2
$$

where $\xi \in (x_0, x_0 + h)$.

Substituting into our secant line approximation:

$$
\begin{split}
\frac{f(x_0 + h) - f(x_0)}{h} &= \frac{1}{h}\left( f(x_0) + h f'(x_0) + \frac{1}{2} h^2 f''(\xi) - f(x_0) \right) \\
&= f'(x_0) + \frac{1}{2} h f''(\xi)
\end{split}
$$
:::

:::{prf:remark} Computational Verification
:label: rmk-fd-verification
:class: dropdown

To verify the error theory: if $E(h) = Ch$, then $\log E = \log C + \log h$.

Plot $\log E$ vs. $\log h$ — you should see a **straight line with slope 1**.
:::

```{figure} /img/FWD_abs_error.jpg
:width: 75%
:align: center

The absolute error of the forward difference approximation of the derivative for $\sin(x_0)$ at $x_0 = 1.2$.
```

## The Round-off Error Trade-off

The plot above shows something surprising: the error **increases** for very small $h$!

This is our first encounter with a fundamental phenomenon in scientific computing: **round-off error**. When $h$ becomes very small, we're subtracting two nearly equal numbers ($f(x_0 + h) \approx f(x_0)$), and then dividing by a tiny $h$. This amplifies the small errors inherent in floating-point arithmetic.

:::{tip}
In the [next chapter on floating-point arithmetic](../error-stability/floating-point.md), we'll learn:
- Why computers can't represent most real numbers exactly
- What **machine epsilon** is and why it's approximately $10^{-16}$ for double precision
- How to find the **optimal step size** that balances truncation and round-off error

For now, the practical takeaway is: **smaller $h$ is not always better**. For first-order finite differences, $h \approx 10^{-8}$ is often optimal.
:::

## The Central Difference Formula

The forward difference only uses information to the *right* of $x_0$. From geometric intuition, using points on **both sides** should give a better approximation to the tangent line—the secant line through $(x_0 - h, f(x_0 - h))$ and $(x_0 + h, f(x_0 + h))$ is more symmetric.

:::{prf:proposition} Central Difference Approximation
:label: prop-central-diff

The central difference approximation:

$$
f'(x_0) \approx \frac{f(x_0 + h) - f(x_0 - h)}{2h}
$$

has error $E(h) = \mathcal{O}(h^2)$ — **second-order** in $h$. Squaring $h$ squares the error!
:::

:::{prf:proof}
:class: dropdown

We expand both $f(x_0 + h)$ and $f(x_0 - h)$ using Taylor series:

$$
f(x_0 \pm h) = f(x_0) \pm h f'(x_0) + \frac{h^2}{2}f''(x_0) \pm \frac{h^3}{6} f'''(\xi^{\pm}) + O(h^4)
$$

Substituting into the central difference formula:

$$
\begin{split}
\frac{f(x_0 + h) - f(x_0 - h)}{2h}
&= \frac{1}{2h}\left( 2hf'(x_0) + \frac{h^3}{6}\left( f'''(\xi^+) + f'''(\xi^-) \right) + O(h^5)\right) \\
&= f'(x_0) + \frac{h^2}{12} \left( f'''(\xi^-) + f'''(\xi^+) \right) + O(h^4)
\end{split}
$$

**Why the improvement?** The even-powered terms ($h^2, h^4, \ldots$) **cancel** due to symmetry.
:::

```{figure} cdf.png
:width: 80%
:align: center

Central finite difference error for approximating $\tan'(0.2)$. The error follows $O(h^2)$ (red line) until roundoff error dominates at $h \approx 10^{-6}$.
```

### Comparison: Forward vs Central

| Method | Formula | Error | Order |
|--------|---------|-------|-------|
| Forward | $\frac{f(x_0 + h) - f(x_0)}{h}$ | $O(h)$ | 1st |
| Central | $\frac{f(x_0 + h) - f(x_0 - h)}{2h}$ | $O(h^2)$ | 2nd |

:::{prf:example} Accuracy Comparison
:label: ex-fd-comparison
:class: dropdown

For $h = 0.01$:
- Forward difference error: $\sim 0.01$
- Central difference error: $\sim 0.0001$

The central difference is **100× more accurate** for the same step size!
:::

:::{prf:remark} Optimal Step Size for Central Differences
:label: rmk-optimal-h-central
:class: dropdown

The total error (truncation + roundoff) is:

$$
E(h) = \frac{M_0 \mu}{h} + M_1 h^2
$$

Minimizing gives $h_{\text{opt}} \approx \mu^{1/3} \approx 6 \times 10^{-6}$ for double precision.

This is *larger* than for forward differences because truncation error decreases faster ($h^2$ vs $h$).
:::

## Second Derivative Approximation

:::{prf:proposition} Second Derivative Formula
:label: prop-second-deriv

The central difference approximation for the second derivative:

$$
f''(x_0) \approx \frac{f(x_0 + h) - 2f(x_0) + f(x_0 - h)}{h^2}
$$

has error $O(h^2)$.
:::

:::{prf:proof}
:class: dropdown

Adding the Taylor expansions (instead of subtracting):

$$
f(x_0 + h) + f(x_0 - h) = 2f(x_0) + h^2 f''(x_0) + O(h^4)
$$

Solving for $f''(x_0)$ gives the formula above.
:::

This formula is fundamental in numerical PDEs—it's the **discrete Laplacian** in one dimension.

## Summary

| Method | Formula | Error Order |
|--------|---------|-------------|
| Forward difference | $\frac{f(x+h) - f(x)}{h}$ | $O(h)$ |
| Central difference | $\frac{f(x+h) - f(x-h)}{2h}$ | $O(h^2)$ |
| Second derivative | $\frac{f(x+h) - 2f(x) + f(x-h)}{h^2}$ | $O(h^2)$ |

**Key principle:** Using symmetric stencils improves accuracy, but there are diminishing returns due to round-off error. The optimal step size balances truncation and round-off errors.
