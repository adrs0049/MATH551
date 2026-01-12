# Condition Numbers

:::{tip} Big Idea
The condition number measures how sensitive a problem is to small perturbations in the input. Large condition numbers signal potential trouble **regardless of the algorithm used**—it's a property of the mathematics, not the implementation.
:::

## Why Condition Numbers Matter

Before we discuss algorithms and their errors, we need to understand: **some problems are inherently harder than others**.

Consider computing $f(x) = \tan(x)$ near $x = \pi/2$. Even with infinite precision arithmetic and a perfect algorithm, tiny uncertainties in $x$ cause huge changes in $\tan(x)$. This isn't a flaw in our algorithm—it's the nature of the tangent function near its pole.

The **condition number** quantifies this sensitivity. It answers: *if my input is slightly wrong, how wrong might my output be?*

## Absolute and Relative Error

Before measuring problem sensitivity, we need precise definitions of error.

:::{prf:definition} Absolute and Relative Error
:label: def-error-types

Let $x$ be the true value and $x^*$ be a numerical approximation.

**Absolute error:**
$$
\text{Abs}(x) = |x - x^*|
$$

**Relative error:**
$$
\text{Rel}(x) = \frac{|x - x^*|}{|x|}, \quad x \neq 0
$$
:::

### Why Relative Error Matters

:::{prf:example} Same Absolute Error, Different Quality
:label: ex-rel-error
:class: dropdown

| | $x = 1$, $x^* = 2$ | $x = 10^6$, $x^* = 10^6 + 1$ |
|--|--|--|
| Absolute error | $1$ | $1$ |
| Relative error | $100\%$ | $0.0001\%$ |
| Quality | Terrible | Excellent |

Both have the same absolute error, but the first is off by 100% while the second is nearly perfect. **Relative error captures what matters.**
:::

## Condition of $f$ at $x$

When we want to evaluate $f(x)$, even if $x^*$ is close to $x$, the value $f(x^*)$ may be far from $f(x)$. The **condition number** quantifies this sensitivity.

:::{prf:definition} Condition Number
:label: def-condition-number

Given a function $f(x)$ and a fixed location $x$, the condition number is:

$$
\kappa := \sup_{|x - x^*| \neq 0} \left| \frac{\text{Rel}(f(x))}{\text{Rel}(x)} \right| = \sup_{|x - x^*| \neq 0} \left| \frac{\frac{f(x) - f(x^*)}{f(x)}}{\frac{x - x^*}{x}} \right|
$$

- If $\kappa \gg 1$ (large), then evaluating $f$ at $x$ is **ill-conditioned**
- If $\kappa \approx \mathcal{O}(1)$, then evaluating $f$ at $x$ is **well-conditioned**
:::

## Simplified Formula via Taylor's Theorem

:::{prf:proposition} Condition Number Formula
:label: prop-condition-formula

For a differentiable function $f$, the condition number at $x$ is:

$$
\kappa = \left| \frac{x f'(x)}{f(x)} \right|
$$
:::

:::{prf:proof}
:class: dropdown

Since $x^*$ is close to $x$, we can write $x^* = x + h$ and apply Taylor's theorem:

$$
f(x^*) = f(x + h) \approx f(x) + f'(x)(x - x^*)
$$

This gives us:

$$
f(x) - f(x^*) \approx f'(x)(x - x^*)
$$

Substituting into the condition number definition:

$$
\kappa = \sup_x \left| \frac{x f'(x)}{f(x)} \right|
$$
:::

## Examples

:::{prf:example} Square Root (Well-Conditioned)
:label: ex-sqrt-condition
:class: dropdown

Consider $f(x) = \sqrt{x}$ with $f'(x) = \frac{1}{2\sqrt{x}}$.

Near $\bar{x} = 1$, we have $\bar{y} = f(\bar{x}) = 1$. Using a Taylor approximation:

$$
y - \bar{y} = \sqrt{x} - 1 \approx \frac{1}{2}(x - 1) = \frac{1}{2}(x - \bar{x})
$$

Variations in $y$ are always smaller than variations in $x$. Computing the condition number:

$$
\kappa = \left| \frac{x f'(x)}{f(x)} \right| = \left| \frac{x \cdot \frac{1}{2\sqrt{x}}}{\sqrt{x}} \right| = \frac{1}{2}
$$

**Result:** $\kappa = 1/2$ — evaluating $\sqrt{x}$ is **well-conditioned**.
:::

:::{prf:example} Tangent Near π/2 (Ill-Conditioned)
:label: ex-tan-condition
:class: dropdown

Consider $f(x) = \tan(x)$ near $x = \frac{\pi}{2}$.

Take two points:
$$
x_1 = \frac{\pi}{2} - 0.001, \quad x_2 = \frac{\pi}{2} - 0.002
$$

Then $|x_1 - x_2| = 0.001$, but $|f(x_1) - f(x_2)| \approx 500$. A tiny input change causes a huge output change!

Computing the condition number:

$$
\kappa = \left| \frac{x f'(x)}{f(x)} \right| = \left| \frac{x}{\cos(x)\sin(x)} \right| = |2x \csc(2x)| \to \infty \text{ as } x \to \pi/2
$$

**Result:** $\kappa \to \infty$ — evaluating $\tan(x)$ near $\frac{\pi}{2}$ is **ill-conditioned**.
:::

:::{prf:example} Logarithm Near 1 (Ill-Conditioned)
:label: ex-log-condition
:class: dropdown

Consider $f(x) = \ln(x)$ near $x = 1$.

$$
\kappa = \left| \frac{x f'(x)}{f(x)} \right| = \left| \frac{x \cdot (1/x)}{\ln(x)} \right| = \frac{1}{|\ln(x)|}
$$

As $x \to 1$, we have $\ln(x) \to 0$, so $\kappa \to \infty$.

**Result:** $\kappa \to \infty$ — evaluating $\ln(x)$ near $x = 1$ is **ill-conditioned**.
:::

## Summary

| Problem | Condition Number | Well/Ill-Conditioned |
|---------|-----------------|---------------------|
| $f(x) = \sqrt{x}$ | $\kappa = 1/2$ | Well-conditioned |
| $f(x) = x^2$ | $\kappa = 2$ | Well-conditioned |
| $f(x) = \tan(x)$ near $\pi/2$ | $\kappa \to \infty$ | Ill-conditioned |
| $f(x) = \ln(x)$ near $1$ | $\kappa \to \infty$ | Ill-conditioned |

The condition number is a property of the **mathematical problem**—it tells us the best possible accuracy any algorithm could achieve. But even well-conditioned problems can give bad results if we use a bad algorithm. That's the subject of the [next section](forward-backward-error.md).
