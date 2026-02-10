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

Single precision gives about **7 decimal digits** of accuracy; double precision gives about **16 decimal digits**. This means that two real numbers that agree to 16 significant digits will have the *same* double-precision representation. We will see exactly how these bits encode numbers when we study [IEEE 754 representation](fast-inverse-sqrt.md).
:::

## Solving the Mystery: The Finite Difference Trade-off

In the [previous chapter](../approximation-theory/numerical-differentiation.md), we observed that finite difference errors **increase** for very small step sizes. Now we have the tools to explain why.

### The Computation on a Machine

On a computer, function values carry round-off error. By the floating-point guarantee {eq}`eq-fp-guarantee`, the computed values satisfy

$$
\widetilde{f}(x+h) = f(x+h)(1 + \varepsilon_1), \quad \widetilde{f}(x) = f(x)(1 + \varepsilon_2), \qquad |\varepsilon_i| \leq \mu.
$$

So the computed forward difference is

$$
\widetilde{D}_h f(x) \;:=\; \frac{\widetilde{f}(x+h) - \widetilde{f}(x)}{h}.
$$

### Separating the Two Errors

Expanding:

$$
\widetilde{D}_h f(x) = \underbrace{\frac{f(x+h) - f(x)}{h}}_{\text{exact finite difference}} \;+\; \underbrace{\frac{f(x+h)\,\varepsilon_1 - f(x)\,\varepsilon_2}{h}}_{\text{round-off perturbation}}.
$$

Applying [Taylor's theorem](../approximation-theory/taylor.md) to the first term:

$$
\frac{f(x+h) - f(x)}{h} = f'(x) + \frac{h}{2}f''(\xi),
$$

the total error is

$$
\widetilde{D}_h f(x) - f'(x) = \underbrace{\frac{f(x+h)\,\varepsilon_1 - f(x)\,\varepsilon_2}{h}}_{\text{round-off error}} \;+\; \underbrace{\frac{h}{2}f''(\xi)}_{\text{truncation error}}.
$$

Since $|\varepsilon_i| \leq \mu$ and $f(x+h) \approx f(x)$ for small $h$:

```{math}
:label: eq-fd-total-error
\left|\widetilde{D}_h f(x) - f'(x)\right| \;\leq\; \underbrace{\frac{2\mu\,|f(x)|}{h}}_{\text{round-off}} \;+\; \underbrace{\frac{h}{2}\,|f''(\xi)|}_{\text{truncation}}
```

:::{note}
The two terms compete: truncation error **decreases** as $h \to 0$, while round-off error **increases** as $h \to 0$.
:::

:::{prf:remark} Catastrophic Cancellation
:label: rmk-cancellation-subtraction

The round-off term $2\mu|f(x)|/h$ grows as $h \to 0$ because we are subtracting nearly equal numbers. The absolute error of $\widetilde{f}(x+h) - \widetilde{f}(x)$ stays roughly constant at $(|f(x+h)| + |f(x)|)\mu$, but dividing by $h$ amplifies it. More generally, subtracting nearly equal floating-point numbers destroys relative accuracy — this is called **catastrophic cancellation**.
:::

### The Optimal Step Size

Minimizing {eq}`eq-fd-total-error` by setting $dE/dh = 0$:

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
