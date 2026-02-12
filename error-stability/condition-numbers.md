# Condition Numbers

:::{tip} Big Idea
Every computation — whether evaluating a mathematical function or performing a step in an algorithm — amplifies input errors by some factor. The **condition number** $\kappa$ measures that amplification: *if my input has relative error $\varepsilon$, how large is the relative error in the output?* When $\kappa$ is large, the computation is **ill-conditioned** and accuracy is lost.
:::

In the [previous section](floating-point.md), we saw that floating-point inputs carry small [relative errors](floating-point.md#def-error-types) bounded by [machine epsilon](floating-point.md#def-machine-epsilon), and that subtracting nearly equal numbers [amplifies those errors catastrophically](floating-point.md#rmk-cancellation-subtraction). The condition number makes this precise — and applies to *any* computation, not just subtraction.

## Condition of $f$ at $x$

When we evaluate $f(x)$ but only have an approximation $x^*$ to the input, the output $f(x^*)$ may be far from $f(x)$. The **condition number** is the worst-case ratio of [relative error](floating-point.md#def-error-types) in the output to relative error in the input.

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

:::{prf:proposition} Condition Number for Function Evaluation
:label: prop-condition-formula

Consider the problem of **evaluating** a differentiable function $f$ at a point $x$. The condition number of this problem is:

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

But note *why* this happens: the condition number blows up because $f(x) = \ln(x) \to 0$, so the relative error $|f(x) - f(x^*)|/|f(x)|$ has a denominator going to zero. The absolute error is perfectly well-behaved — $|f'(1)| = 1$, so a small input perturbation produces a small absolute output perturbation. The "ill-conditioning" here is really that **relative error is undefined at a zero of $f$**. This happens for any function near one of its roots, and is an artifact of measuring error in relative terms rather than a genuine sensitivity of the computation.
:::

## Condition Number of Subtraction

The condition number formula $\kappa = |xf'(x)/f(x)|$ works for functions of one variable. For a function of **two** variables $g(a, b) = a - b$, we generalize: the condition number measures how the relative error in $g$ relates to the relative errors in $a$ and $b$.

:::{prf:example} Subtraction is Ill-Conditioned When $a \approx b$
:label: ex-subtraction-condition
:class: dropdown

Consider $g(a, b) = a - b$ with small perturbations $\tilde{a} = a(1 + \varepsilon_1)$ and $\tilde{b} = b(1 + \varepsilon_2)$:

$$
\tilde{g} - g = a\varepsilon_1 - b\varepsilon_2
$$

The relative error in the result is:

$$
\frac{|\tilde{g} - g|}{|g|} \leq \frac{|a| + |b|}{|a - b|} \cdot \max(|\varepsilon_1|, |\varepsilon_2|)
$$

The condition number of subtraction is therefore:

$$
\kappa = \frac{|a| + |b|}{|a - b|}
$$

- When $a$ and $b$ are well-separated: $\kappa \approx \mathcal{O}(1)$ — subtraction is well-conditioned.
- When $a \approx b$: $\kappa \to \infty$ — subtraction is **ill-conditioned**.

This is precisely the [catastrophic cancellation](floating-point.md#rmk-cancellation-subtraction) we saw in the floating-point chapter, now explained through the lens of condition numbers.
:::

## Condition Number of the Finite Difference

We can now give a **condition number explanation** of the [finite difference trade-off](floating-point.md#solving-the-mystery-the-finite-difference-trade-off).

:::{prf:example} Finite Difference Condition Number
:label: ex-fd-condition
:class: dropdown

The forward difference $\frac{f(x_0 + h) - f(x_0)}{h}$ requires computing $a - b$ where $a = f(x_0 + h)$ and $b = f(x_0)$.

Applying the subtraction condition number with $a \approx b \approx f(x_0)$:

$$
\kappa \approx \frac{|f(x_0 + h)| + |f(x_0)|}{|f(x_0 + h) - f(x_0)|} \approx \frac{2|f(x_0)|}{h|f'(x_0)|}
$$

As $h \to 0$, this condition number grows like $1/h$. Each input carries relative error $\mu$ (machine epsilon), so the round-off contribution to the finite difference is:

$$
\text{round-off error} \approx \kappa \cdot \mu \cdot |f'(x_0)| \approx \frac{2\mu|f(x_0)|}{h}
$$

This matches exactly the round-off term we derived from the floating-point analysis — but now we see it as a **conditioning problem**: the subtraction step is ill-conditioned for small $h$.
:::

## Summary

The condition number $\kappa$ measures how much a problem amplifies input errors.
It is a property of the **mathematical problem**, not the algorithm — if $\kappa$
is large, no algorithm can avoid accuracy loss. But what if $\kappa$ is small and
we still get poor results? Then the problem isn't to blame — the algorithm is.
In the [next section](stability.md), we make this precise with the concept of
**stable and unstable algorithms**.

