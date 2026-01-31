# Chebyshev Polynomials

:::{tip} Big Idea
Chebyshev polynomials solve the minimax problem for polynomial interpolation: they tell us where to place nodes to minimize the worst-case interpolation error.
:::

## The Minimax Problem

Recall the interpolation error formula:
$$
|f(x) - p(x)| = \left|\frac{f^{(n)}(\xi)}{n!}\right| \cdot |(x-x_1)(x-x_2) \cdots (x-x_n)|
$$

We can't control the first factor (it depends on $f$), but we can minimize the second by choosing nodes wisely.

**The Minimax Problem:** Choose $x_1, \ldots, x_n \in [-1, 1]$ to minimize
$$
\max_{x \in [-1,1]} |(x-x_1)(x-x_2) \cdots (x-x_n)|
$$

## The Solution: Chebyshev Roots

:::{prf:theorem} Optimality of Chebyshev Nodes
:label: thm-chebyshev-optimal

The points minimizing the minimax problem on $[-1, 1]$ are the **Chebyshev roots**:
$$
x_k = \cos\left(\frac{(2k-1)\pi}{2n}\right), \quad k = 1, \ldots, n
$$

The minimum value is $\frac{1}{2^{n-1}}$.
:::

For $n = 10$, the Chebyshev roots are:
$$
x_k = \cos\left(\frac{(2k-1)\pi}{20}\right), \quad k = 1, \ldots, 10
$$

These points cluster near the endpoints of the interval—exactly where Runge's phenomenon causes problems with equally spaced nodes!

:::{prf:remark} Geometric Interpretation
Chebyshev nodes arise from projecting equally spaced points on the unit circle onto the $x$-axis. Place $n$ points at angles $\theta_k = (2k-1)\pi/(2n)$ on the circle; their $x$-coordinates are the Chebyshev nodes $x_k = \cos\theta_k$. This projection naturally clusters points near $x = \pm 1$.
:::

## Chebyshev Polynomials

```{figure} ../img/cheb.png
:width: 85%
:align: center

The first six Chebyshev polynomials $T_0(x)$ through $T_5(x)$. Note that each $T_k$ oscillates between $-1$ and $+1$ exactly $k+1$ times on $[-1, 1]$.
```

### Definition via Complex Exponentials

The most illuminating definition connects Chebyshev polynomials to Fourier series and complex analysis.

:::{prf:definition} Chebyshev Polynomials
:label: def-chebyshev-poly

Let $z = e^{i\theta}$ be a point on the unit circle with $x = \Re(z) = \cos\theta$. Then:
$$
T_k(x) = \frac{1}{2}(z^k + z^{-k}) = \cos(k\theta)
$$
where $\theta = \arccos(x)$.
:::

This definition directly links **Chebyshev series** to **Fourier series** and **Laurent series**—three fundamental tools connected through the change of variables $x = \cos\theta$.

### Why They're Polynomials

From the definition, we verify the first few:
- $T_0(x) = \frac{1}{2}(z^0 + z^0) = 1$
- $T_1(x) = \frac{1}{2}(z + z^{-1}) = \frac{1}{2}(e^{i\theta} + e^{-i\theta}) = \cos\theta = x$
- $T_2(x) = \cos(2\theta) = 2\cos^2\theta - 1 = 2x^2 - 1$
- $T_3(x) = \cos(3\theta) = 4\cos^3\theta - 3\cos\theta = 4x^3 - 3x$
- $T_4(x) = 8x^4 - 8x^2 + 1$

### Recurrence Relation

The recurrence follows from the complex exponential definition:

$$
\frac{1}{2}(z + z^{-1})(z^k + z^{-k}) = \frac{1}{2}(z^{k+1} + z^{-k-1}) + \frac{1}{2}(z^{k-1} + z^{-k+1})
$$

which gives:

:::{prf:property} Three-Term Recurrence
:label: prop-chebyshev-recurrence

$$
T_{k+1}(x) = 2x T_k(x) - T_{k-1}(x)
$$
:::

By induction, each $T_k(x)$ is a polynomial of degree exactly $k$.

### Key Properties

1. **Roots:** $T_n$ has $n$ roots at the Chebyshev nodes
2. **Extrema:** $T_n$ oscillates between $-1$ and $+1$ exactly $n+1$ times on $[-1, 1]$
3. **Leading coefficient:** The coefficient of $x^n$ in $T_n$ is $2^{n-1}$ (for $n \geq 1$)
4. **Minimax property:** $\frac{1}{2^{n-1}} T_n(x)$ is the monic polynomial of degree $n$ with smallest maximum on $[-1, 1]$

## Proof of Optimality

The proof uses a contradiction argument based on the alternation property.

:::{prf:proof}
:class: dropdown

Let $p(x) = (x - x_1) \cdots (x - x_n)$ be any monic polynomial of degree $n$, and suppose
$$
\max_{x \in [-1,1]} |p(x)| < \frac{1}{2^{n-1}}
$$

Let $q(x) = \frac{1}{2^{n-1}} T_n(x)$, which is monic and alternates between $\pm\frac{1}{2^{n-1}}$ at $n+1$ points $y_1, \ldots, y_{n+1}$.

At each $y_i$:
- If $q(y_i) = +\frac{1}{2^{n-1}}$, then $q(y_i) - p(y_i) > 0$
- If $q(y_i) = -\frac{1}{2^{n-1}}$, then $q(y_i) - p(y_i) < 0$

So $q - p$ alternates sign at least $n+1$ times, meaning it has at least $n$ roots.

But $q - p$ has degree at most $n - 1$ (both are monic of degree $n$), so it can have at most $n - 1$ roots. Contradiction!
:::

## Transformation to Arbitrary Intervals

The Chebyshev nodes are defined for $[-1, 1]$. For an arbitrary interval $[a, b]$, we transform:

$$
x_k = \frac{b-a}{2} \cos\left(\frac{(2k-1)\pi}{2n}\right) + \frac{a+b}{2}
$$

This:
1. Scales the Chebyshev points by $\frac{b-a}{2}$ (the ratio of interval lengths)
2. Shifts the center from 0 to $\frac{a+b}{2}$

The minimum value of the node polynomial on $[a, b]$ becomes $\frac{(b-a)^n}{2^{2n-1}}$.

## Chebyshev vs. Equally Spaced: Runge's Function

Consider $f(x) = \frac{1}{1 + 25x^2}$ on $[-1, 1]$.

```{figure} ../img/cheb_conv.png
:width: 95%
:align: center

**Chebyshev vs. equally spaced nodes:** Convergence comparison for polynomial interpolation. Chebyshev nodes achieve exponential convergence for smooth functions, while equally spaced nodes fail for Runge's function.
```

| Nodes | Error behavior as $n \to \infty$ |
|-------|----------------------------------|
| Equally spaced | Error **grows** without bound |
| Chebyshev | Error **decreases** exponentially |

This is the power of optimal node placement!

## When to Use Chebyshev Interpolation

Chebyshev interpolation is particularly valuable for:

1. **Spectral methods** for differential equations
2. **Gaussian quadrature** for numerical integration
3. **Approximating smooth functions** with high accuracy

However, for most practical data fitting, **piecewise polynomials (splines)** are often preferred because:
- They're simpler to work with
- They don't require special node placement
- They naturally handle non-smooth data

## Chebyshev Series

Beyond interpolation, Chebyshev polynomials provide a powerful **series representation** for functions, analogous to Taylor series but with superior convergence properties.

### Lipschitz Continuity

:::{prf:definition} Lipschitz Continuous
:label: def-lipschitz

A function $f$ is **Lipschitz continuous** on $[-1, 1]$ if there exists a constant $C > 0$ such that:
$$
|f(x) - f(y)| \leq C|x - y| \quad \text{for all } x, y \in [-1, 1]
$$
:::

Any function $f \in C^1[-1, 1]$ is Lipschitz (by the mean value theorem), but not every continuous function is Lipschitz.

### The Chebyshev Series Theorem

:::{prf:theorem} Existence of Chebyshev Series
:label: thm-chebyshev-series

If $f$ is Lipschitz continuous on $[-1, 1]$, then it has a unique representation:
$$
f(x) = \sum_{k=0}^{\infty} c_k T_k(x)
$$
which converges absolutely and uniformly.
:::

The coefficients are given by the integral formula:

:::{prf:property} Chebyshev Coefficients
:label: prop-chebyshev-coefficients

$$
c_k = \frac{2}{\pi} \int_{-1}^{1} \frac{f(x) T_k(x)}{\sqrt{1 - x^2}} \, dx \quad \text{for } k \geq 1
$$
and $c_0$ by the same formula with factor $\frac{1}{\pi}$ instead of $\frac{2}{\pi}$.
:::

:::{dropdown} Connection to Fourier Series
Under the change of variables $x = \cos\theta$, we have $T_k(x) = \cos(k\theta)$ and:
$$
\frac{dx}{\sqrt{1-x^2}} = -d\theta
$$
So the Chebyshev coefficients are precisely Fourier cosine coefficients of $f(\cos\theta)$!
:::

### Why Chebyshev Over Taylor?

| Aspect | Taylor Series | Chebyshev Series |
|--------|--------------|------------------|
| Center | Single point $x_0$ | Entire interval $[-1, 1]$ |
| Convergence | Disk of convergence | Whole interval (uniform) |
| Truncation | Best near $x_0$ | Near-best everywhere |
| Complexity | $O(n^3)$ via Vandermonde | $O(n \log n)$ via DCT |

For functions on an interval, Chebyshev series are almost always preferable.

## From Interpolation to Approximation

So far, we've treated interpolation as a data-fitting problem: given points $(x_i, f_i)$, find a polynomial passing through them. But there's a more powerful perspective.

**Key insight:** If we have a *function* $f$, we can sample it at $n+1$ Chebyshev nodes to get data, then interpolate. The natural question becomes: *how well does the interpolant approximate the original function?*

This shift—from fitting given data to approximating a known function—is the foundation of **spectral methods**. The Chebyshev interpolant $p_n$ becomes a computational stand-in for $f$: we differentiate, integrate, or solve equations using $p_n$ instead of $f$.

## Convergence Rates: A Preview

The smoothness of $f$ determines the convergence rate:

| Smoothness | Coefficient Decay | Approximation Error |
|------------|-------------------|---------------------|
| $f \in C^k$ | $O(n^{-k-1})$ | $O(n^{-k})$ (algebraic) |
| $f$ analytic | $O(\rho^{-n})$ | $O(\rho^{-n})$ (exponential) |

The jump from algebraic to exponential convergence when $f$ becomes analytic is dramatic—this is what makes spectral methods so powerful for smooth problems.

See the [Spectral Accuracy](spectral-accuracy.md) chapter for the precise theorems, including:
- The **Bernstein ellipse** and how complex singularities control convergence
- Coefficient decay rates and their relationship to smoothness
- Examples showing how to diagnose function smoothness from coefficients

