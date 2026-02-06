---
kernelspec:
  name: python3
  display_name: Python 3
execute:
  enabled: true
---

# Floating Point Numbers

:::{tip} Big Idea
Computers represent real numbers using finite precision floating point arithmetic. Understanding this representation explains why certain computations lose accuracy and helps us write more robust numerical code.
:::

## Number Systems

### Integers in Base 10

We typically write integers in base 10. A number $d_N d_{N-1} \dots d_1 d_0$ (where $d_i \in \{0, 1, \dots, 9\}$) represents:

$$
\# = d_0 \cdot 10^0 + d_1 \cdot 10^1 + \cdots + d_N \cdot 10^N
$$

### Binary Representation

Computers use base 2 (binary), where digits $d_i \in \{0, 1\}$. A binary number $d_N d_{N-1} \dots d_1 d_0$ represents:

$$
\# = \sum_{i=0}^{N} d_i \cdot 2^i
$$

:::{prf:example} Binary to Decimal
:label: ex-binary-decimal
:class: dropdown

The binary number $110_2 = 1 \cdot 2^2 + 1 \cdot 2^1 + 0 \cdot 2^0 = 6$.

**Bit shifting:**
- Left shift: $110_2 \to 1100_2 = 12$ (doubles the number)
- Right shift: $110_2 \to 11_2 = 3$ (halves the number)
:::

### Signed Integers: Two's Complement

Negative integers use **two's complement** representation. For an $N$-bit signed integer $a = d_{N-1}d_{N-2}\dots d_1 d_0$:

$$
a = -d_{N-1} \cdot 2^{N-1} + \sum_{i=0}^{N-2} d_i \cdot 2^i
$$

:::{prf:example} Two's Complement Negation
:label: ex-twos-complement
:class: dropdown

To get $-5$ from $5$ in 8-bit two's complement:

1. Start with $5 = 0000\,0101_2$
2. Invert all bits: $1111\,1010_2$
3. Add one: $1111\,1011_2 = -5$
:::

## Fixed Point Notation

To represent fractions, we allow digits after a radix point:

$$
d_N \dots d_1 d_0 . d_{-1} d_{-2} \dots d_{-M}
$$

represents:

$$
\# = \sum_{i=-M}^{N} d_i \cdot b^i
$$

where $b$ is the base.

:::{prf:example} Fixed Point Binary
:label: ex-fixed-point
:class: dropdown

In binary: $1.01_2 = 1 + 0 \cdot 2^{-1} + 1 \cdot 2^{-2} = 1.25$
:::

With finitely many digits, some numbers cannot be represented exactly (e.g., $\frac{1}{3}$ or $\pi$).

## Floating Point Numbers

For scientific computing, we need to represent numbers of vastly different magnitudes—from Avogadro's number ($6.02 \times 10^{23}$) to Planck's constant ($6.63 \times 10^{-34}$).

**Scientific notation** allows the radix point to "float":

$$
245000 = 2.45 \times 10^5
$$

In binary:
$$
11000_2 = 1.1_2 \times 2^4, \quad 0.0101_2 = 1.01_2 \times 2^{-2}
$$

Note: In normalized binary scientific notation, the digit before the radix point is always 1, so we don't need to store it!

## IEEE 754 Standard

A floating point number consists of three parts:
- **Sign bit** $S$: 0 for positive, 1 for negative
- **Exponent** $E$: Shifted to allow negative exponents
- **Mantissa/Fraction** $m$: The significant digits

### Single Precision (32-bit)

| Component | Bits |
|-----------|------|
| Sign | 1 |
| Exponent | 8 |
| Mantissa | 23 |

The value represented is:

$$
\text{fl}(x) = \pm \left(1 + \frac{d_0}{2^1} + \frac{d_1}{2^2} + \cdots + \frac{d_{22}}{2^{23}}\right) \times 2^{E - 127}
$$

where $d_i$ are the mantissa bits and $E$ is the stored exponent.

### Double Precision (64-bit)

| Component | Bits |
|-----------|------|
| Sign | 1 |
| Exponent | 11 |
| Mantissa | 52 |

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

## Rounding Error (Machine Epsilon)

Given a real number $x$, its floating point representation $\text{fl}(x)$ satisfies:

$$
\frac{|\text{fl}(x) - x|}{|x|} \leq \mu
$$

where $\mu$ is the **machine epsilon** (or **unit roundoff**):

$$
\mu = \frac{1}{2} \times 2^{-t}
$$

with $t$ being the number of mantissa bits.

| Precision | Mantissa bits | Machine epsilon |
|-----------|--------------|-----------------|
| Single (float) | 23 | $\approx 5.96 \times 10^{-8}$ |
| Double | 52 | $\approx 1.11 \times 10^{-16}$ |

:::{prf:remark} Practical Implication
:label: rmk-precision-digits

- Single precision gives about **7 decimal digits** of accuracy
- Double precision gives about **16 decimal digits** of accuracy
:::

## Application: The Finite Difference Trade-off

In the [previous chapter](../approximation-theory/numerical-differentiation.md), we observed that finite difference errors **increase** for very small step sizes. Machine epsilon explains why.

### Round-off Error in Subtraction

The machine epsilon bound can be rewritten in multiplicative form: for any real number $a$,

$$
\text{fl}(a) = a(1 + \varepsilon), \quad |\varepsilon| \leq \mu
$$

Now consider computing $a - b$ when both values carry representation error:

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

