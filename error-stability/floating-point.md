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

## Integer Interpretation of Floating Point

The same 32 bits can be interpreted as either a float or an integer. Given the floating point representation $(S_x, E_x, m_x)$, the integer value is:

$$
\text{Int}(x) = 2^{31} S_x + 2^{23} E_x + M_x
$$

where $M_x = 2^{23} m_x$ is the integer value of the mantissa bits.

This dual interpretation is exploited in fast numerical algorithms like the [fast inverse square root](fast-inverse-sqrt.md).

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

In the [approximation theory chapter](../approximation-theory/finite-differences.md), we observed that finite difference errors **increase** for very small step sizes. Now we can explain why.

### The Total Error

When computing the forward difference approximation:

$$
f'(x) \approx \frac{f(x+h) - f(x)}{h}
$$

we make **two types of errors**:

1. **Truncation error** from Taylor's theorem: $\frac{h}{2}|f''(\xi)|$
2. **Round-off error** from floating-point arithmetic

For the round-off error: when $h$ is small, $f(x+h) \approx f(x)$, so we're subtracting two nearly equal numbers. If both values have relative error $\mu$, the subtraction has absolute error roughly $2\mu|f(x)|$. Dividing by $h$ amplifies this to:

$$
\text{Round-off error} \approx \frac{2\mu|f(x)|}{h}
$$

The total error is:

$$
E(h) = \underbrace{\frac{2\mu|f(x)|}{h}}_{\text{round-off}} + \underbrace{\frac{h}{2}|f''(\xi)|}_{\text{truncation}}
$$

### The Optimal Step Size

To minimize $E(h)$, differentiate and set to zero:

$$
\frac{dE}{dh} = -\frac{2\mu|f(x)|}{h^2} + \frac{|f''(\xi)|}{2} = 0
$$

Solving (and assuming $|f(x)| \approx |f''(\xi)|$ for simplicity):

$$
h_{\text{opt}} = 2\sqrt{\mu} \approx \sqrt{\mu}
$$

:::{prf:remark} Optimal Step Sizes
:label: rmk-optimal-step-sizes

| Precision | Machine epsilon $\mu$ | Optimal $h$ for FD |
|-----------|----------------------|-------------------|
| Single | $\sim 10^{-8}$ | $\sim 10^{-4}$ |
| Double | $\sim 10^{-16}$ | $\sim 10^{-8}$ |

For double precision, the optimal forward difference step is $h \approx 10^{-8}$.
:::

### Connection to the Error Framework

This is a perfect illustration of our error analysis framework:

- **Backward error** (truncation): How much did we perturb the mathematical problem? We approximated $f'$ by a secant slope—error $O(h)$.
- **Forward error** (round-off amplification): How much did floating-point errors affect our answer? Error $O(\mu/h)$.

The condition number of the operation "subtract then divide by $h$" grows like $1/h$, which is why round-off errors get amplified for small $h$.

:::{admonition} Key Lesson
:class: tip

**Smaller step sizes are not always better.** There is a sweet spot where truncation and round-off errors balance. This principle appears throughout scientific computing—in finite differences, numerical integration, and iterative solvers.
:::
