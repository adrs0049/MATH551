---
kernelspec:
  name: python3
  display_name: Python 3
---

# Errors and Floating-Point Arithmetic

:::{tip} Big Idea
Every computation on a computer introduces small errors because real numbers must be stored in finite precision. The key guarantee: the relative error of representing any real number is bounded by **machine epsilon**. This simple fact explains the mysterious finite difference behavior from [Chapter 1](../approximation-theory/numerical-differentiation.md) — and understanding it is the first step toward writing robust numerical code.
:::

## Absolute and Relative Error

Before we can analyze floating-point arithmetic, we need precise definitions of error.

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

## Floating-Point Numbers as Approximations

We need to represent all real numbers — from Avogadro's number ($6.02 \times 10^{23}$) to Planck's constant ($6.63 \times 10^{-34}$) — using finite resources. The solution is **floating-point arithmetic**: a finite set of numbers $\mathbb{F}$ that approximates the real line $\mathbb{R}$, with a guaranteed bound on the relative error of the approximation.

The key fact is this: for any real number $x$, its floating-point representation $\text{fl}(x)$ satisfies

```{math}
:label: eq-fp-guarantee
\text{fl}(x) = x(1 + \varepsilon), \quad |\varepsilon| \leq \mu
```

where $\mu$ is **machine epsilon** (or **unit roundoff**). This says that floating-point representation introduces a small *relative* error — the same relative accuracy whether $x$ is tiny or huge.

:::{prf:definition} Machine Epsilon
:label: def-machine-epsilon

The **machine epsilon** $\mu$ is the smallest number such that

$$
\frac{|\text{fl}(x) - x|}{|x|} \leq \mu \quad \text{for all } x \in \mathbb{R}
$$

For a floating-point system with $t$ mantissa bits:

$$
\mu = \frac{1}{2} \times 2^{-t}
$$

| Precision | Mantissa bits | Machine epsilon | Decimal digits |
|-----------|--------------|-----------------|----------------|
| Single (float32) | 23 | $\approx 5.96 \times 10^{-8}$ | ~7 |
| Double (float64) | 52 | $\approx 1.11 \times 10^{-16}$ | ~16 |
:::

:::{prf:remark} Practical Implication
:label: rmk-precision-digits

Single precision gives about **7 decimal digits** of accuracy; double precision gives about **16 decimal digits**. This means that two real numbers that agree to 16 significant digits will have the *same* double-precision representation.
:::

## The Floating-Point Number Line

Every floating-point number is a rational number of the form $m \cdot 2^e$, so the set of representable numbers $\mathbb{F} \subset \mathbb{Q}$ is finite. But unlike a uniform grid, floating-point numbers are **not evenly spaced**. Within each exponent range $[2^e, 2^{e+1})$, the $2^t$ mantissa values divide the interval into equally spaced points, with gap size:

$$
\Delta = 2^{e - t}
$$

where $t$ is the number of mantissa bits. When the exponent increases by 1, the gap doubles. The result is a logarithmic distribution: dense near zero, sparse for large magnitudes.

The figure below illustrates this for a toy system with 3-bit mantissa and exponents $e \in \{0, 1, 2, 3\}$. Notice how the spacing doubles at each power of 2.

```{code-cell} python
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt

# Toy floating-point system: 3-bit mantissa, exponents 0..3
t = 3  # mantissa bits
floats = []
for e in range(4):  # exponents 0, 1, 2, 3
    for m in range(2**t):  # mantissa values 0, 1, ..., 7
        x = (1 + m / 2**t) * 2**e
        floats.append(x)

floats = sorted(set(floats))

fig, ax = plt.subplots(figsize=(10, 1.5))
ax.scatter(floats, np.zeros_like(floats), marker='|', s=400, linewidths=2.5, color='C0')

# Mark powers of 2 to show exponent boundaries
for e in range(5):
    ax.axvline(2**e, color='gray', linestyle='--', linewidth=0.7, alpha=0.5)
    ax.text(2**e, 0.25, f'$2^{e}$', ha='center', fontsize=10, color='gray')

ax.set_xlim(0.8, 17)
ax.set_ylim(-0.3, 0.45)
ax.set_xlabel('Value')
ax.set_title('Floating-point numbers (3-bit mantissa): spacing doubles at each power of 2')
ax.set_yticks([])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
plt.tight_layout()
plt.show()
```

:::{prf:remark} Relative Spacing is Constant
:label: rmk-relative-spacing

The absolute gap $\Delta = 2^{e-t}$ varies with magnitude, but the **relative** gap is nearly constant:

$$
\frac{\Delta}{2^e} = 2^{-t}
$$

This is why machine epsilon — a *relative* error bound — characterizes floating-point precision. No matter how large or small a number is, the nearest representable neighbor is within a fixed relative distance.
:::

:::{prf:remark} Floating Point as a Finite Approximation to $\mathbb{R}$
:label: rmk-fp-approximation-reals
:class: dropdown

In analysis, the real numbers $\mathbb{R}$ are constructed as the **completion** of the rationals $\mathbb{Q}$ — every gap between rationals is filled by taking limits of Cauchy sequences. The reals are the smallest set where every convergent sequence has a limit.

Floating-point arithmetic attempts something analogous with finite resources: approximate $\mathbb{R}$ using a finite subset $\mathbb{F} \subset \mathbb{Q}$, chosen so that every real number has a nearby representative. The logarithmic spacing ensures this works across many orders of magnitude — but unlike $\mathbb{R}$, the gaps never close. Machine epsilon is the price we pay for finiteness: no matter how many bits we use, there is always a smallest relative gap below which distinct real numbers become indistinguishable.
:::

## Solving the Mystery: The Finite Difference Trade-off

In the [previous chapter](../approximation-theory/numerical-differentiation.md), we observed that finite difference errors **increase** for very small step sizes. Now we have the tools to explain why.

### Round-off Error in Subtraction

Using the floating-point guarantee {eq}`eq-fp-guarantee`, consider computing $a - b$ when both values carry representation error:

$$
\text{fl}(a) - \text{fl}(b) = a(1 + \varepsilon_1) - b(1 + \varepsilon_2) = (a - b) + \underbrace{(a\varepsilon_1 - b\varepsilon_2)}_{\text{error}}
$$

The absolute error of the subtraction is bounded by:

$$
|a\varepsilon_1 - b\varepsilon_2| \leq |a|\mu + |b|\mu = (|a| + |b|)\mu
$$

:::{prf:remark} Catastrophic Cancellation
:label: rmk-cancellation-subtraction

When $a \approx b$, the absolute error $(|a| + |b|)\mu$ stays roughly constant, but $|a - b|$ becomes tiny. The **relative** error of the subtraction explodes:

$$
\frac{(|a| + |b|)\mu}{|a - b|} \to \infty \quad \text{as } a \to b
$$

Subtracting nearly equal numbers destroys relative accuracy.
:::

### Applying This to Finite Differences

In the forward difference $\frac{f(x+h) - f(x)}{h}$, we subtract $a = f(x+h)$ and $b = f(x)$. When $h$ is small, both are close to $f(x)$, so the subtraction error is:

$$
(|f(x+h)| + |f(x)|)\mu \approx 2|f(x)|\mu
$$

Dividing by $h$ gives the round-off contribution to the finite difference:

$$
\text{round-off error} \approx \frac{2\mu|f(x)|}{h}
$$

Combined with the truncation error from [Taylor's theorem](../approximation-theory/taylor.md), the total error is:

$$
E(h) = \underbrace{\frac{2\mu|f(x)|}{h}}_{\text{round-off}} + \underbrace{\frac{h}{2}|f''(\xi)|}_{\text{truncation}}
$$

:::{note}
The two terms compete: truncation error **decreases** as $h \to 0$, while round-off error **increases** as $h \to 0$.
:::

### The Optimal Step Size

Minimizing $E(h)$ by setting $dE/dh = 0$:

$$
-\frac{2\mu|f(x)|}{h^2} + \frac{|f''(\xi)|}{2} = 0 \quad \implies \quad h_{\text{opt}} = 2\sqrt{\frac{\mu|f(x)|}{|f''(\xi)|}} \approx \sqrt{\mu}
$$

where the last approximation assumes $|f(x)| \approx |f''(\xi)|$.

:::{prf:remark} Optimal Step Sizes
:label: rmk-optimal-step-sizes

| Precision | Machine epsilon $\mu$ | Optimal $h$ (forward diff) |
|-----------|----------------------|---------------------------|
| Single | $\sim 10^{-8}$ | $\sim 10^{-4}$ |
| Double | $\sim 10^{-16}$ | $\sim 10^{-8}$ |
:::
