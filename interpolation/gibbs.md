# The Gibbs Phenomenon

:::{tip} Big Idea
Polynomial approximations to discontinuous functions exhibit persistent **overshoot** near the discontinuity. This overshoot (~9% for step functions) does not decrease as the polynomial degree increases—it's a fundamental property of polynomial approximation, not a numerical artifact.
:::

## The Phenomenon

Consider the sign function $f(x) = \text{sign}(x)$ on $[-1, 1]$:
$$
\text{sign}(x) = \begin{cases} -1 & x < 0 \\ 0 & x = 0 \\ +1 & x > 0 \end{cases}
$$

When we interpolate this discontinuous function with a polynomial:

```python
import numpy as np

def chebpts(n):
    """Chebyshev points of the second kind."""
    return np.cos(np.pi * np.arange(n) / (n - 1))

def bary(xs, f, x):
    """Barycentric interpolation."""
    n = len(x)
    w = np.ones(n)
    w[0] = 0.5
    w[-1] = 0.5
    w[1::2] *= -1

    result = np.zeros_like(xs)
    for i, xeval in enumerate(xs):
        if xeval in x:
            result[i] = f[np.where(x == xeval)[0][0]]
        else:
            terms = w / (xeval - x)
            result[i] = np.dot(terms, f) / np.sum(terms)
    return result

xs = np.linspace(-1, 1, 10000)

for n in [10, 50, 100]:
    x = chebpts(n)
    f = np.sign(x)
    p = bary(xs, f, x)

    print(f"n={n}: max overshoot = {np.max(p):.4f}")
    # Output: always around 1.09 regardless of n!
```

## The Overshoot Persists

```{figure} ../img/gibbs.png
:width: 95%
:align: center

The Gibbs phenomenon: polynomial interpolation of a step function shows persistent ~9% overshoot near the discontinuity, regardless of the polynomial degree.
```

| Degree $n$ | Maximum of interpolant |
|------------|----------------------|
| 10 | ~1.09 |
| 50 | ~1.09 |
| 100 | ~1.09 |
| 1000 | ~1.09 |

The overshoot **does not decrease** as $n \to \infty$!

More precisely, for the sign function:
$$
\lim_{n \to \infty} \max_{x \in [-1,1]} p_n(x) = \frac{2}{\pi} \int_0^\pi \frac{\sin t}{t} \, dt \approx 1.0895
$$

This is approximately **9% overshoot**.

## Why Does This Happen?

The Gibbs phenomenon occurs because:

1. **Discontinuities create slowly-decaying coefficients:** For $\text{sign}(x)$, the Chebyshev coefficients decay as $O(k^{-1})$—too slow for pointwise convergence at the jump.

2. **Convergence is non-uniform:** While $\|f - p_n\|_2 \to 0$, the $L^\infty$ error near the discontinuity does not vanish.

3. **Oscillations concentrate:** As $n$ increases, the oscillations don't disappear—they become more numerous but more localized near the jump.

## Decay Away from Discontinuities

The oscillations do decay as we move away from the discontinuity:

- At distance $d$ from the jump: oscillations decay like $O(1/(nd))$
- For $f'$ having the jump: decay like $O(1/(nd)^2)$

So the "ringing" is **localized** but **persistent** at the discontinuity itself.

## Practical Solutions

### 1. Piecewise Interpolation

If you know where the discontinuities are, interpolate each smooth piece separately:

```python
def piecewise_interp(f, discontinuities, n_per_piece):
    """Interpolate between known discontinuities."""
    pieces = []
    bounds = [-1] + discontinuities + [1]

    for i in range(len(bounds) - 1):
        a, b = bounds[i], bounds[i+1]
        # Interpolate on [a, b] separately
        x_local = (b-a)/2 * chebpts(n_per_piece) + (a+b)/2
        f_local = f(x_local)
        pieces.append((a, b, x_local, f_local))

    return pieces
```

### 2. Filtering/Smoothing

Apply a filter to damp high-frequency oscillations (at the cost of some accuracy).

### 3. Edge Detection

Use the coefficient decay pattern to automatically detect discontinuities, then apply piecewise methods.

## Connection to Fourier Series

The Gibbs phenomenon was first observed in **Fourier series** (by Henry Wilbraham in 1848, later rediscovered by J. Willard Gibbs in 1899).

The same 9% overshoot appears when approximating a step function with Fourier series:
$$
\text{sign}(x) \approx \frac{4}{\pi}\left(\sin x + \frac{\sin 3x}{3} + \frac{\sin 5x}{5} + \cdots\right)
$$

Since Chebyshev series are related to Fourier series (via $x = \cos\theta$), the same phenomenon appears.

## Summary

:::{prf:remark} Key Takeaways
:label: rmk-gibbs-takeaways

1. **The overshoot is real:** ~9% for step functions, and it persists as $n \to \infty$
2. **It's not a bug:** This is a fundamental property of polynomial/Fourier approximation
3. **Oscillations localize:** They become more concentrated near the discontinuity
4. **Practical solution:** Use piecewise interpolation when discontinuity locations are known
:::

**The lesson:** Polynomial interpolation is designed for **smooth** functions. When approximating discontinuous functions, either:
- Accept the oscillations, or
- Break the problem into smooth pieces
