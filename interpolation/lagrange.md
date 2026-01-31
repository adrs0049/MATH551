# Lagrange Interpolation

:::{tip} Big Idea
Given data points $(x_i, y_i)$, there exists a unique polynomial of degree at most $n$ passing through all points. The Lagrange form gives an explicit formula, but **barycentric interpolation** is the numerically stable way to evaluate it.
:::

## The Interpolation Problem

**Problem:** Given $n+1$ distinct points $(x_0, y_0), (x_1, y_1), \ldots, (x_n, y_n)$, find a polynomial $p(x)$ such that:
$$
p(x_i) = y_i \quad \text{for } i = 0, 1, \ldots, n
$$

:::{prf:theorem} Existence and Uniqueness of Polynomial Interpolation
:label: thm-interp-existence

Given $n+1$ points with distinct $x$-coordinates, there exists a **unique** polynomial of degree at most $n$ passing through all points.
:::

:::{prf:proof}
:class: dropdown

**Existence:** The Lagrange form (below) explicitly constructs such a polynomial.

**Uniqueness:** Suppose $p(x)$ and $q(x)$ are both polynomials of degree at most $n$ interpolating the same data. Consider the difference $d(x) = p(x) - q(x)$.

- $d(x)$ is a polynomial of degree at most $n$
- $d(x_i) = p(x_i) - q(x_i) = y_i - y_i = 0$ for all $i = 0, 1, \ldots, n$

So $d(x)$ has $n+1$ roots. But a nonzero polynomial of degree $n$ can have at most $n$ roots. Therefore $d(x) \equiv 0$, which means $p(x) = q(x)$.
:::

## Lagrange Basis Polynomials

For nodes $x_0, x_1, \ldots, x_n$, define the **Lagrange basis polynomial**:

$$
L_j(x) = \prod_{i=0, i \neq j}^{n} \frac{x - x_i}{x_j - x_i}
$$

**Key property:**
$$
L_j(x_i) = \begin{cases} 1 & \text{if } i = j \\ 0 & \text{if } i \neq j \end{cases}
$$

The interpolating polynomial is then:
$$
p(x) = \sum_{j=0}^{n} y_j L_j(x)
$$

:::{prf:example} Lagrange Interpolation Through Three Points
:label: ex-lagrange-three-points

Find the Lagrange interpolation function through the points $(0, 1), (1, 0), (2, 3)$.

:::{dropdown} Solution

**Basis polynomials:**
$$
L_0(x) = \frac{(x-1)(x-2)}{(0-1)(0-2)} = \frac{(x-1)(x-2)}{2}
$$
$$
L_1(x) = \frac{(x-0)(x-2)}{(1-0)(1-2)} = -x(x-2)
$$
$$
L_2(x) = \frac{(x-0)(x-1)}{(2-0)(2-1)} = \frac{x(x-1)}{2}
$$

**Interpolant:**
$$
p(x) = 1 \cdot L_0(x) + 0 \cdot L_1(x) + 3 \cdot L_2(x) = 2x^2 - 3x + 1
$$

**Verify:** $p(0) = 1$, $p(1) = 0$, $p(2) = 3$
:::

:::{prf:remark} Limitations of Direct Lagrange Evaluation
:label: rmk-lagrange-limitations

The Lagrange formula is elegant but has practical limitations:

| Issue | Problem |
|-------|---------|
| **Cost** | $O(n^2)$ operations per evaluation point |
| **Stability** | Products of many terms can overflow/underflow |
| **Flexibility** | Adding a new point requires recomputing all basis polynomials |

See the [Interpolation Methods notebook](./interpolation-methods.ipynb) for a computational comparison.
:::

## Barycentric Interpolation

Can we do better? Yes—by algebraically reorganizing the Lagrange formula, we obtain the **barycentric form**, which resolves all three issues while computing the same polynomial.

### The Node Polynomial

Define the **node polynomial**:
$$
\ell(x) = \prod_{k=0}^{n} (x - x_k) = (x - x_0)(x - x_1) \cdots (x - x_n)
$$

Its derivative at a node $x_j$ is:
$$
\ell'(x_j) = \prod_{k \neq j} (x_j - x_k)
$$

This lets us rewrite the Lagrange basis as:
$$
L_j(x) = \frac{\ell(x)}{\ell'(x_j)(x - x_j)}
$$

### Barycentric Weights

Define the **barycentric weights**:
$$
\lambda_j = \frac{1}{\ell'(x_j)} = \frac{1}{\prod_{k \neq j} (x_j - x_k)}
$$

Then the Lagrange basis becomes:
$$
L_j(x) = \ell(x) \cdot \frac{\lambda_j}{x - x_j}
$$

### First Barycentric Formula

Substituting into the interpolant gives the **first barycentric formula**:

:::{prf:proposition} First Barycentric Formula
:label: prop-bary-first

$$
p(x) = \ell(x) \sum_{j=0}^{n} \frac{\lambda_j}{x - x_j} f_j
$$
:::

Once the weights $\lambda_j$ are known (computed once in $O(n^2)$ time), evaluating $p(x)$ costs only $O(n)$.

### Second Barycentric Formula

An even more elegant formula comes from the identity $\sum_{j=0}^{n} L_j(x) = 1$ (the Lagrange polynomials form a partition of unity—they interpolate the constant function 1).

Dividing the interpolant by this identity:

:::{prf:proposition} Second Barycentric Formula
:label: prop-bary-second

$$
p(x) = \frac{\displaystyle\sum_{j=0}^{n} \frac{\lambda_j}{x - x_j} f_j}{\displaystyle\sum_{j=0}^{n} \frac{\lambda_j}{x - x_j}}
$$

with the convention $p(x_j) = f_j$ when $x = x_j$.
:::

:::{dropdown} Why $\sum L_j(x) = 1$?
The function $g(x) = 1$ (constant) is interpolated by itself. Since $g(x_j) = 1$ for all $j$:
$$
1 = \sum_{j=0}^{n} 1 \cdot L_j(x) = \sum_{j=0}^{n} L_j(x)
$$
This follows from the uniqueness of polynomial interpolation.
:::

### Advantages

1. **$O(n)$ per evaluation** after $O(n^2)$ preprocessing for weights
2. **Numerically stable** even for large $n$
3. **Adding a point** only requires updating weights

### Implementation

```python
def bary_weights(x):
    """Compute barycentric weights for nodes x."""
    n = len(x)
    w = np.ones(n)
    for j in range(n):
        for i in range(n):
            if i != j:
                w[j] /= (x[j] - x[i])
    return w

def bary_interp(xeval, x, y, w):
    """Evaluate interpolant at xeval using barycentric formula."""
    # Handle evaluation at nodes
    for j, xj in enumerate(x):
        if np.isclose(xeval, xj):
            return y[j]

    terms = w / (xeval - x)
    return np.dot(terms, y) / np.sum(terms)
```

## Chebyshev Points and Barycentric Weights

:::{prf:definition} Chebyshev Points
:label: def-chebyshev-points

The **Chebyshev points** (of the first kind) on $[-1, 1]$ are:
$$
x_k = \cos\left(\frac{k\pi}{n}\right), \quad k = 0, 1, \ldots, n
$$

These $n+1$ points are the projections of equally spaced points on the unit circle onto the $x$-axis. They cluster near the endpoints $\pm 1$—this uneven distribution will be important later.
:::

For Chebyshev points, the barycentric weights have a remarkably simple closed form:

:::{prf:theorem} Chebyshev Barycentric Weights
:label: thm-cheb-bary-weights

For the $n+1$ Chebyshev points, the barycentric weights are:
$$
\lambda_k = (-1)^k \delta_k, \quad \text{where } \delta_k = \begin{cases} 1/2 & k = 0 \text{ or } k = n \\ 1 & \text{otherwise} \end{cases}
$$
:::

This makes Chebyshev interpolation especially efficient—no $O(n^2)$ weight computation needed!

```python
def cheb_bary_weights(n):
    """Barycentric weights for Chebyshev points."""
    w = np.ones(n+1)
    w[0] = 0.5
    w[-1] = 0.5
    w[1::2] *= -1
    return w
```

## Numerical Stability

The two barycentric formulas have different stability properties. The following results are due to {cite:t}`Higham2004`.

### Backward Stability of the First Formula

:::{prf:theorem} Backward Stability of Modified Lagrange Formula
:label: thm-bary-backward-stable

The first barycentric formula (modified Lagrange form)
$$
p(x) = \ell(x) \sum_{j=0}^{n} \frac{\lambda_j}{x - x_j} f_j
$$
is **backward stable**. The computed value $\hat{p}(x)$ satisfies
$$
\hat{p}(x) = \sum_{j=0}^{n} f_j(1 + \delta_j) L_j(x)
$$
where $|\delta_j| \lesssim nu$ and $u$ is the unit roundoff ($u \approx 10^{-16}$ in double precision).

:::{dropdown} Interpretation
Backward stability means the computed result is the *exact* interpolant of slightly perturbed data. The perturbations are at the level of machine precision—the best we can hope for.
:::

### Forward Stability of the Second Formula

:::{prf:theorem} Forward Stability of Barycentric Formula
:label: thm-bary-forward-stable

The second barycentric formula
$$
p(x) = \frac{\sum_{j=0}^{n} \frac{\lambda_j}{x - x_j} f_j}{\sum_{j=0}^{n} \frac{\lambda_j}{x - x_j}}
$$
is **forward stable** when the Lebesgue constant $\Lambda_n$ is small. Specifically:
$$
\frac{|\hat{p}(x) - p(x)|}{|p(x)|} \lesssim \Lambda_n \cdot nu
$$

For Chebyshev points, $\Lambda_n = O(\log n)$, so the formula is stable for all practical $n$.

:::{dropdown} Why the Lebesgue constant?
The Lebesgue constant measures how much the interpolation process amplifies errors in the data. For well-chosen nodes (Chebyshev), $\Lambda_n$ grows slowly. For equispaced nodes, $\Lambda_n$ grows exponentially, making even the barycentric formula unreliable.
:::

### Comparison

| Formula | Stability Type | Best For |
|---------|---------------|----------|
| First (modified Lagrange) | Backward stable | Extrapolation, any nodes |
| Second (barycentric) | Forward stable | Interpolation with good nodes |

:::{prf:remark} Practical Recommendation
:label: rmk-which-formula

For $x \in [-1, 1]$ with Chebyshev nodes, both formulas give comparable accuracy. The second formula is preferred because:
1. Scale invariance allows rescaling weights to avoid overflow
2. No need to compute $\ell(x)$

For extrapolation or equispaced nodes, use the first formula.

See {cite:t}`BerrutTrefethen2004` and {cite:t}`WebbTrefethenGonnet2012` for further analysis.
:::

### Why Scale Invariance Matters

The second barycentric formula is **scale invariant**: we can multiply all weights by any nonzero constant without changing the result:
$$
\frac{\sum_j \frac{c\lambda_j}{x - x_j} f_j}{\sum_j \frac{c\lambda_j}{x - x_j}} = \frac{\sum_j \frac{\lambda_j}{x - x_j} f_j}{\sum_j \frac{\lambda_j}{x - x_j}}
$$

This means we can rescale weights to avoid overflow, which is critical for large $n$.

### Warning: Equispaced Points

For equispaced points, the weights grow like:
$$
|\lambda_j| \sim \frac{2^n}{n!} \binom{n}{j}
$$

This grows exponentially, making polynomial interpolation through equispaced points numerically unstable—even with the barycentric formula. This is separate from (but compounds) Runge's phenomenon.

## Runge's Phenomenon: A Warning

Consider $f(x) = \frac{1}{1 + 25x^2}$ on $[-1, 1]$.

```{figure} ../img/runge.png
:width: 95%
:align: center

**Runge's phenomenon:** Polynomial interpolation of $f(x) = 1/(1+25x^2)$ through equally spaced nodes. As the degree increases, the interpolant develops large oscillations near the boundaries, despite the function being smooth.
```

```{figure} ../img/error_interp.png
:width: 95%
:align: center

**Interpolation error comparison:** (Left) Equidistant points cause error to grow exponentially with $n$. (Center) Chebyshev points cluster near endpoints. (Right) Chebyshev interpolation error decreases exponentially—the hallmark of spectral accuracy.
```

| Nodes | Error as $n \to \infty$ |
|-------|------------------------|
| Equally spaced | **Grows** without bound |
| Chebyshev | **Decreases** exponentially |

```python
# Equally spaced: error GROWS
x_eq = np.linspace(-1, 1, n)
f_eq = 1/(1 + 25*x_eq**2)
# Interpolant oscillates wildly near boundaries!

# Chebyshev: error DECREASES
x_cheb = np.cos(np.pi * np.arange(n) / (n-1))
f_cheb = 1/(1 + 25*x_cheb**2)
# Smooth, accurate approximation
```

The lesson: **node placement matters**. Chebyshev nodes cluster near the endpoints, exactly where equally spaced nodes cause trouble.

**Use barycentric interpolation.** It's the standard in modern software (Chebfun, etc.).

## References

```{bibliography}
:filter: docname in docnames
```
