# Spectral Accuracy and Coefficient Decay

:::{tip} Big Idea
The rate at which Chebyshev coefficients decay reveals the smoothness of a function. For analytic functions, coefficients decay exponentially—this is **spectral accuracy**. For functions with limited smoothness, coefficients decay algebraically, and the rate tells us exactly how many derivatives exist.
:::

## Total Variation

:::{prf:definition} Total Variation
:label: def-total-variation

The **total variation** of a function $f$ on $[-1, 1]$ is:
$$
V = \|f'\|_1 = \int_{-1}^{1} |f'(x)| \, dx
$$
:::

For discontinuous $f'$, interpret this as the sum of $|f'|$ integrals plus the jump magnitudes.

## The Fundamental Theorem

:::{prf:theorem} Coefficient Decay
:label: thm-coefficient-decay

Let $f$ and its derivatives through $f^{(p-1)}$ be absolutely continuous on $[-1, 1]$, and suppose $f^{(p)}$ has bounded variation $V$. Then for $k \geq p+1$, the Chebyshev coefficients satisfy:
$$
|c_k| \leq \frac{2V}{\pi (k-p)^{p+1}}
$$
:::

**Interpretation:**
- $f$ continuous ($p=0$): $|c_k| = O(k^{-1})$
- $f$ has one derivative ($p=1$): $|c_k| = O(k^{-2})$
- $f$ has $p$ derivatives: $|c_k| = O(k^{-p-1})$
- $f$ analytic: $|c_k| = O(\rho^{-k})$ for some $\rho > 1$

## Analytic Functions: The Bernstein Ellipse

For analytic functions, we can give a precise geometric characterization of the convergence rate.

:::{prf:theorem} Exponential Convergence for Analytic Functions
:label: thm-cheb-exponential-convergence

Let $f$ be analytic on $[-1, 1]$ and analytically continuable to the **Bernstein ellipse** $\mathcal{E}_\rho$ with parameter $\rho > 1$. Then the Chebyshev interpolant $p_n$ satisfies:
$$
\|f - p_n\|_\infty \leq \frac{2M}{\rho^n(\rho - 1)}
$$
where $M = \max_{z \in \mathcal{E}_\rho} |f(z)|$.
:::

:::{prf:definition} Bernstein Ellipse
:label: def-bernstein-ellipse

The **Bernstein ellipse** $\mathcal{E}_\rho$ is the image of the circle $|z| = \rho$ under the Joukowski map $z \mapsto \frac{1}{2}(z + z^{-1})$. It has:
- Foci at $\pm 1$
- Semi-major axis $a = \frac{1}{2}(\rho + \rho^{-1})$
- Semi-minor axis $b = \frac{1}{2}(\rho - \rho^{-1})$
:::

The largest $\rho$ for which $f$ is analytic inside $\mathcal{E}_\rho$ determines the convergence rate. **Singularities in the complex plane—even far from $[-1, 1]$—limit $\rho$ and slow convergence.**

:::{prf:example} Runge's Function Revisited
:label: ex-runge-bernstein

For $f(x) = \frac{1}{1 + 25x^2}$, the singularities are at $x = \pm i/5$. The Bernstein ellipse passing through these points has $\rho \approx 1.2$, giving convergence rate $O(1.2^{-n})$—exponential, but slow.

For $f(x) = \frac{1}{1 + x^2}$, the singularities at $x = \pm i$ are farther away, giving $\rho \approx 2.4$ and much faster convergence.
:::

## Approximation Error

The error in truncating the Chebyshev series at $n$ terms:

$$
\|f - p_n\|_\infty \leq \sum_{k=n+1}^{\infty} |c_k|
$$

Using the coefficient decay:
- $f \in C^p$: $\|f - p_n\|_\infty = O(n^{-p})$ (algebraic)
- $f$ analytic: $\|f - p_n\|_\infty = O(\rho^{-n})$ (exponential)

| Smoothness | Coefficient Decay | Approximation Error |
|------------|-------------------|---------------------|
| Discontinuous | $O(k^{-1})$ | $O(n^{-1/2})$ |
| Continuous | $O(k^{-1})$ | $O(n^{-1})$ |
| $C^1$ | $O(k^{-2})$ | $O(n^{-1})$ |
| $C^p$ | $O(k^{-p-1})$ | $O(n^{-p})$ |
| Analytic | $O(\rho^{-k})$ | $O(\rho^{-n})$ |

## Example: The Sign Function

$f(x) = \text{sign}(x)$ is discontinuous at $x = 0$.

```python
n = 2000
x = chebpts(n)
f = np.sign(x)
c = vals2coeffs(f)

# Coefficients decay as O(k^{-1})
plt.semilogy(np.abs(c))
# Only odd coefficients are nonzero (antisymmetric function)
```

Approximation error: $\|f - p_n\|_\infty = O(n^{-1})$ (not $O(n^{-1/2})$ because Chebyshev avoids Gibbs phenomenon in max norm).

## Example: $|x|$

$f(x) = |x|$ is continuous but not differentiable at $x = 0$.

```python
n = 1000
x = chebpts(n)
f = np.abs(x)
c = vals2coeffs(f)

# Coefficients decay as O(k^{-2})
# Approximation error: O(n^{-1})
```

The error bound:
$$
\|f - p_n\|_\infty \leq \frac{8}{\pi(n-1)}
$$

## Example: $|\sin(5x)|^3$

This function has three continuous derivatives (jumps in $f^{(3)}$ at roots of $\sin(5x)$).

```python
f = lambda x: np.abs(np.sin(5*x))**3

n = 10000
x = chebpts(n)
c = vals2coeffs(f(x))

# Coefficients decay as O(k^{-4})
# Approximation error: O(n^{-3})
```

## Example: Analytic Function

$f(x) = \sin(6x) + \sin(60 e^x)$ is entire (analytic everywhere).

```python
f = lambda x: np.sin(6*x) + np.sin(60*np.exp(x))

n = 250
x = chebpts(n)
c = vals2coeffs(f(x))

# Coefficients decay EXPONENTIALLY
# Around n=150, coefficients hit machine precision
```

The coefficients plateau at $\sim 10^{-15}$ (machine precision), indicating the series has converged.

## Visualizing Coefficient Decay

```python
fig, axes = plt.subplots(3, 2, figsize=(10, 8))

functions = [
    (np.sign, r"$\mathrm{sign}(x)$", "Discontinuous"),
    (np.abs, r"$|x|$", "Continuous, not $C^1$"),
    (lambda x: np.abs(np.sin(5*x))**3, r"$|\sin(5x)|^3$", "$C^3$"),
]

for i, (f, label, desc) in enumerate(functions):
    x = chebpts(1000)
    vals = f(x)
    coeffs = vals2coeffs(vals)

    axes[i, 0].plot(x, vals, 'k')
    axes[i, 0].set_title(f"{label} — {desc}")

    axes[i, 1].semilogy(np.abs(coeffs), 'k.')
    axes[i, 1].set_ylabel("$|c_k|$")
    axes[i, 1].set_ylim([1e-16, 1e1])
```

## The Gibbs Phenomenon

For discontinuous functions, polynomial approximations exhibit **overshoot** near discontinuities:

- The maximum of the interpolant exceeds the function by about 9% (for step functions)
- This overshoot persists as $n \to \infty$
- It's a fundamental property of polynomial approximation, not a numerical artifact

```python
xs = np.linspace(-1, 1, 10000)

for n in [10, 50, 100]:
    x = chebpts(n)
    f = np.sign(x)
    p = bary(xs, f, x)

    print(f"n={n}: max overshoot = {np.max(p):.4f}")
    # Always around 1.09 regardless of n
```

## Adaptive Resolution

The coefficient decay principle enables **adaptive** algorithms:

```python
def adaptive_approx(f, tol=1e-12):
    """Find n such that coefficients decay below tol."""
    for n in [16, 32, 64, 128, 256, 512, 1024]:
        x = chebpts(n)
        c = vals2coeffs(f(x))

        # Check if last few coefficients are small
        if np.max(np.abs(c[-5:])) < tol * np.max(np.abs(c)):
            # Trim to significant coefficients
            k = np.where(np.abs(c) > tol * np.max(np.abs(c)))[0][-1]
            return c[:k+1], n

    raise ValueError("Function requires more than 1024 points")
```

This is the core idea behind Chebfun's automatic degree selection.

## Locating Non-Smoothness

The coefficient decay pattern reveals **where** smoothness breaks down:

- Smooth everywhere: exponential decay
- Singularity at one point: algebraic decay
- Multiple singularities: decay rate limited by worst singularity

For functions with boundary layers or internal layers, the coefficients reveal the layer location and sharpness.

## Lebesgue Constants

The **Lebesgue constant** quantifies how close interpolation is to best approximation.

:::{prf:definition} Lebesgue Constant
:label: def-lebesgue-constant

$$
\Lambda_n = \sup_{f \neq 0} \frac{\|p_n\|_\infty}{\|f\|_\infty}
$$
where $p_n$ is the interpolant of $f$ at $n+1$ nodes.
:::

:::{prf:theorem} Near-Best Approximation
:label: thm-near-best-approximation

Let $p_n$ be the Chebyshev interpolant and $p_n^*$ be the best polynomial approximation of degree $n$. Then:
$$
\|f - p_n\|_\infty \leq (\Lambda_n + 1) \|f - p_n^*\|_\infty
$$
:::

The Lebesgue constant tells us how much worse interpolation can be compared to best approximation.

### Growth Rates by Node Type

| Nodes | $\Lambda_n$ growth |
|-------|-------------------|
| Chebyshev | $\frac{2}{\pi}\log(n+1) + O(1)$ (logarithmic) |
| Equispaced | $\frac{2^{n+1}}{en\log n}$ (exponential!) |
| Lower bound (any nodes) | $\frac{2}{\pi}\log(n+1) + \gamma$ |

Chebyshev nodes are **nearly optimal**—their Lebesgue constant grows only logarithmically.

## The Weierstrass Approximation Theorem

A foundational result guaranteeing that polynomials can approximate any continuous function:

:::{prf:theorem} Weierstrass (1885)
:label: thm-weierstrass

Let $f$ be continuous on $[-1, 1]$. For any $\varepsilon > 0$, there exists a polynomial $p$ such that:
$$
\|f - p\|_\infty < \varepsilon
$$
:::

This theorem applies even to **pathological** continuous functions (like those that are continuous but nowhere differentiable).

### The Faber-Bernstein Result

However, Weierstrass does not guarantee *uniform* convergence of interpolants:

:::{prf:theorem} Faber-Bernstein Theorem
:label: thm-faber-bernstein

There exists **no** set of interpolation nodes such that polynomial interpolation converges for **all** continuous functions.
:::

This means: for any node choice, there exists some continuous function whose interpolants diverge. However, for the **vast majority of functions** encountered in practice (those with some smoothness), Chebyshev interpolation works phenomenally well.

## Practical Guidelines

1. **Check coefficient decay** to verify your approximation is resolved
2. **Exponential decay** to machine precision indicates an analytic function
3. **Algebraic decay** suggests limited smoothness—count the rate to find $p$
4. **Plateauing coefficients** at some level indicates numerical noise or rounding

## Summary

| What You See | What It Means |
|--------------|---------------|
| $\|c_k\| \sim \rho^{-k}$ | Analytic function |
| $\|c_k\| \sim k^{-p-1}$ | Function has $p$ derivatives |
| $\|c_k\| \sim k^{-1}$ | Discontinuous function |
| Plateau at $10^{-15}$ | Machine precision reached |
| No decay | Function not resolved, need more points |

**The key insight:** Chebyshev coefficients encode smoothness. Spectral methods achieve their remarkable accuracy because smooth functions have rapidly decaying coefficients.
