# Approximation Theory

:::{tip} Big Idea
Given a function $f$, how can we represent it on a computer? We approximate $f$ by a polynomial, and work in two equivalent representations: **values** at special points, or **coefficients** in a polynomial basis. Transforming between these views—via the DCT/FFT—unlocks spectral accuracy.
:::

## Overview

Computers can only perform arithmetic. To work with functions like $\sin(x)$ or solutions to differential equations, we need approximations. Polynomials are the natural choice:

- **Easy to evaluate** (Horner's method, barycentric formulas)
- **Easy to differentiate and integrate**
- **Converge rapidly** for smooth functions

The key insight from Trefethen's *Approximation Theory and Approximation Practice* is that we can work in two dual spaces:

| Representation | Description | Advantages |
|---------------|-------------|------------|
| **Value space** | $f(x_0), f(x_1), \ldots, f(x_n)$ | Direct sampling, differentiation matrices |
| **Coefficient space** | $c_0, c_1, \ldots, c_n$ in $p(x) = \sum c_k T_k(x)$ | Integration, coefficient decay analysis |

The **DCT** (Discrete Cosine Transform) converts between them in $O(n \log n)$ operations.

## Why Chebyshev?

For polynomial interpolation on $[-1, 1]$, node placement matters enormously:

- **Equally spaced nodes**: Runge phenomenon—error grows exponentially!
- **Chebyshev nodes**: Error decreases exponentially for smooth functions

Chebyshev polynomials $T_n(x) = \cos(n \arccos x)$ are optimal in a minimax sense.

## Chapter Contents

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} Lagrange Interpolation
:link: lagrange
:link-type: doc

The interpolation problem and barycentric formulas for stable evaluation.
:::

:::{grid-item-card} Chebyshev Polynomials
:link: chebyshev
:link-type: doc

Optimal node placement, minimax properties, and defeating Runge's phenomenon.
:::

:::{grid-item-card} Values and Coefficients
:link: values-coefficients
:link-type: doc

The dual representations and the DCT/FFT transforms between them.
:::

:::{grid-item-card} Differentiation Matrices
:link: differentiation
:link-type: doc

Computing derivatives in value space—spectral differentiation.
:::

:::{grid-item-card} Integration
:link: integration
:link-type: doc

Quadrature via Chebyshev coefficients vs. trapezoidal rule.
:::

:::{grid-item-card} Spectral Accuracy
:link: spectral-accuracy
:link-type: doc

Why smooth functions achieve exponential convergence—coefficient decay theory.
:::

::::

## Learning Outcomes

After completing this chapter, you should be able to:

- **L7.1:** State existence/uniqueness for polynomial interpolation.
- **L7.2:** Construct Lagrange basis polynomials and use barycentric evaluation.
- **L7.3:** Explain Runge's phenomenon and why Chebyshev nodes avoid it.
- **L7.4:** Define Chebyshev polynomials via $T_n(x) = \cos(n\arccos x)$.
- **L7.5:** Transform between values and coefficients using DCT.
- **L7.6:** Construct and apply Chebyshev differentiation matrices.
- **L7.7:** Integrate using Chebyshev coefficients.
- **L7.8:** Relate coefficient decay to function smoothness.
- **L7.9:** Explain spectral vs algebraic convergence rates.

## Key Formulas

| Concept | Formula |
|---------|---------|
| Chebyshev polynomial | $T_n(x) = \cos(n \arccos x)$ |
| Chebyshev nodes | $x_k = \cos\left(\frac{k\pi}{n}\right)$, $k = 0, \ldots, n$ |
| Barycentric weights | $w_k = (-1)^k \cdot \begin{cases} 1/2 & k = 0, n \\ 1 & \text{otherwise} \end{cases}$ |
| Values → Coefficients | DCT (Type I) |
| Coefficients → Values | Inverse DCT |
| Chebyshev integral | $\int_{-1}^{1} T_k(x)\,dx = \frac{2}{1-k^2}$ (even $k$) |
| Coefficient decay | $\|c_k\| = O(k^{-p-1})$ if $f^{(p)}$ has bounded variation |

## Applications

- **Function approximation**: Represent $\sin$, $\exp$, special functions
- **Spectral methods**: Solve ODEs/PDEs with exponential accuracy
- **Quadrature**: High-accuracy numerical integration
- **Root finding**: Chebyshev companion matrix methods

## Connection to Fourier

Under the substitution $x = \cos\theta$:
- Chebyshev polynomials become cosines: $T_n(\cos\theta) = \cos(n\theta)$
- The DCT is essentially a real FFT
- Chebyshev interpolation $\leftrightarrow$ cosine interpolation

This connection is why Chebyshev methods achieve "spectral accuracy"—the same exponential convergence as Fourier methods, but for non-periodic functions on $[-1, 1]$.
