# Numerical Integration

:::{tip} Big Idea
Integrals can be approximated by replacing the integrand with simple functions (polynomials) and integrating those exactly. Unlike differentiation, integration is a *smoothing* operation—errors tend to average out rather than amplify.
:::

## The Trapezoidal Rule

The simplest approach: approximate $f(x)$ by a straight line and integrate that.

### Single Interval

:::{prf:proposition} Trapezoidal Rule (Single Interval)
:label: prop-trapezoidal

$$
\int_a^b f(x)\,dx \approx \frac{b-a}{2}\left(f(a) + f(b)\right)
$$

Connect $(a, f(a))$ to $(b, f(b))$ with a line; compute the trapezoid's area.
:::

### Local Error Analysis

What is the error of this approximation? Let $h = b - a$ and use Taylor's theorem.

:::{prf:theorem} Local Truncation Error for Trapezoidal Rule
:label: thm-trap-local-error

For $f \in C^2[a, b]$, the **local truncation error** (error on a single interval) is:

$$
\int_a^b f(x)\,dx - \frac{h}{2}\left(f(a) + f(b)\right) = -\frac{h^3}{12}f''(\xi)
$$

for some $\xi \in (a, b)$, where $h = b - a$.
:::

:::{prf:proof}
:class: dropdown

Expand $f(x)$ around the midpoint $c = (a+b)/2$ using Taylor's theorem:

$$
f(x) = f(c) + f'(c)(x - c) + \frac{f''(\eta(x))}{2}(x - c)^2
$$

Integrating from $a$ to $b$:

$$
\int_a^b f(x)\,dx = hf(c) + \frac{f''(\eta)}{2}\int_a^b (x-c)^2\,dx = hf(c) + \frac{h^3}{24}f''(\eta)
$$

where we used that $\int_a^b (x-c)\,dx = 0$ by symmetry, and $\int_a^b (x-c)^2\,dx = h^3/12$.

For the trapezoidal approximation, expand $f(a)$ and $f(b)$ around $c$:

$$
f(a) = f(c) - \frac{h}{2}f'(c) + \frac{h^2}{8}f''(\xi_1)
$$

$$
f(b) = f(c) + \frac{h}{2}f'(c) + \frac{h^2}{8}f''(\xi_2)
$$

Adding:

$$
\frac{h}{2}(f(a) + f(b)) = hf(c) + \frac{h^3}{16}\cdot\frac{f''(\xi_1) + f''(\xi_2)}{2}
$$

Taking the difference and using the intermediate value theorem to combine the $f''$ terms:

$$
\int_a^b f(x)\,dx - \frac{h}{2}(f(a) + f(b)) = -\frac{h^3}{12}f''(\xi)
$$

for some $\xi \in (a, b)$.
:::

**Key observation:** The local error is $O(h^3)$—cubic in the interval width.

## Composite Trapezoidal Rule

For better accuracy, divide $[a, b]$ into $n$ subintervals of equal width $h = (b-a)/n$, with nodes $x_k = a + kh$ for $k = 0, 1, \ldots, n$.

:::{prf:definition} Composite Trapezoidal Rule
:label: def-composite-trapezoidal

$$
T_n(f) = h\left(\frac{f(x_0) + f(x_n)}{2} + \sum_{k=1}^{n-1} f(x_k)\right)
$$

This applies the trapezoidal rule to each subinterval and sums the results.
:::

### From Local to Global Error

The **global error** is the total error when approximating the integral over $[a, b]$.

:::{prf:theorem} Global Error for Composite Trapezoidal Rule
:label: thm-trap-global-error

For $f \in C^2[a, b]$:

$$
\int_a^b f(x)\,dx - T_n(f) = -\frac{(b-a)h^2}{12}f''(\xi) = O(h^2)
$$

for some $\xi \in (a, b)$.
:::

:::{prf:proof}
:class: dropdown

On each subinterval $[x_{k-1}, x_k]$, the local error is:

$$
\int_{x_{k-1}}^{x_k} f(x)\,dx - \frac{h}{2}(f(x_{k-1}) + f(x_k)) = -\frac{h^3}{12}f''(\xi_k)
$$

for some $\xi_k \in (x_{k-1}, x_k)$.

Summing over all $n$ subintervals:

$$
\int_a^b f(x)\,dx - T_n(f) = -\frac{h^3}{12}\sum_{k=1}^{n} f''(\xi_k)
$$

Since $f'' \in C[a,b]$, the sum $\frac{1}{n}\sum_{k=1}^{n} f''(\xi_k)$ lies between $\min f''$ and $\max f''$. By the intermediate value theorem, there exists $\xi \in (a, b)$ such that:

$$
\frac{1}{n}\sum_{k=1}^{n} f''(\xi_k) = f''(\xi)
$$

Therefore:

$$
\int_a^b f(x)\,dx - T_n(f) = -\frac{h^3}{12} \cdot n \cdot f''(\xi) = -\frac{h^3 n}{12}f''(\xi)
$$

Since $n = (b-a)/h$:

$$
\int_a^b f(x)\,dx - T_n(f) = -\frac{(b-a)h^2}{12}f''(\xi)
$$
:::

### Understanding Local vs Global Error

| Error Type | Definition | Trapezoidal Rule |
|------------|------------|------------------|
| **Local** | Error on one subinterval of width $h$ | $O(h^3)$ |
| **Global** | Total error over $[a, b]$ | $O(h^2)$ |

**Why does the order drop from 3 to 2?**

The global error accumulates local errors from $n \sim 1/h$ subintervals:

$$
\text{Global error} \sim n \times \text{Local error} \sim \frac{1}{h} \times h^3 = h^2
$$

This is the typical pattern: **global order = local order − 1**.

:::{prf:remark} Python Implementation
:label: rmk-trapezoidal-code
:class: dropdown

```python
def trapezoidal(f, a, b, n):
    """Composite trapezoidal rule with n subintervals."""
    h = (b - a) / n
    x = np.linspace(a, b, n + 1)
    y = f(x)
    return h * (0.5 * y[0] + np.sum(y[1:-1]) + 0.5 * y[-1])
```
:::

## Higher-Order Methods

By using higher-degree polynomial approximations, we can achieve better accuracy.

**Simpson's rule** uses a quadratic polynomial through three points $(a, f(a))$, $(m, f(m))$, $(b, f(b))$ where $m = (a+b)/2$:

$$
\int_a^b f(x)\,dx \approx \frac{h}{6}\left(f(a) + 4f(m) + f(b)\right)
$$

where $h = b - a$. This has local error $O(h^5)$ and global error $O(h^4)$—two orders better than trapezoidal.

Even higher-order **Newton-Cotes formulas** exist (using more equally-spaced points), though they become unstable for high orders. The optimal approach—**Gaussian quadrature**—chooses both the nodes and weights optimally and achieves remarkable efficiency. We will explore this in the chapter on [interpolation](../interpolation/integration.md).

## Why Integration is Easier Than Differentiation

:::{prf:remark} Smoothing vs. Roughening
:label: rmk-integration-easier
:class: dropdown

| Operation | Error behavior |
|-----------|---------------|
| Differentiation | Errors **amplify** (dividing by small $h$) |
| Integration | Errors **average out** (summing many terms) |

This is why numerical integration is generally more stable than numerical differentiation. Integration "smooths," differentiation "roughens."
:::

From the perspective of conditioning:
- **Differentiation** amplifies high-frequency noise
- **Integration** damps high-frequency components

This is why we can often integrate noisy data reliably, but differentiating noisy data is notoriously difficult.

## Summary

| Rule | Local Error | Global Error |
|------|-------------|--------------|
| Trapezoidal | $O(h^3)$ | $O(h^2)$ |
| Simpson's | $O(h^5)$ | $O(h^4)$ |

**Key principle:** Global order = local order − 1, because we sum $O(1/h)$ local errors.
