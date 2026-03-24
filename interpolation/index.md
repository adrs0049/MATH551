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

