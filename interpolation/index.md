# Approximation Theory

:::{tip} Big Idea
This chapter extends the linear algebra you have learned to *function
spaces*. Polynomials of degree at most $n$ form an $(n+1)$-dimensional
vector space. A function is then a vector. Sampling at $n+1$ nodes is a
choice of basis. Differentiation and integration become a linear operator
and a linear functional on this space. Approximation theory studies which
basis to pick, where to put the nodes, and how the smoothness of the
underlying function controls the error.
:::

## What This Chapter Does

We start from polynomial interpolation as a change of basis between values
and coefficients, and identify which bases are numerically usable. We then
ask where to place the nodes, with the answer being Chebyshev points.
Coefficient decay reveals the smoothness of the underlying function, and
that same dictionary controls the accuracy of every downstream operation:
differentiation by the matrix $D_N$, integration by Clenshaw–Curtis
quadrature, and the solution of linear boundary value problems by spectral
collocation. The chapter ends with adaptivity, where the coefficient tail
itself tells us when we have used enough basis functions.

## Learning Outcomes

After completing this chapter you should be able to:

- **L7.1** State existence and uniqueness of polynomial interpolation, and
  identify the Vandermonde map between values and coefficients in any
  basis.
- **L7.2** Use the barycentric formula to evaluate a Chebyshev interpolant
  in $O(n)$.
- **L7.3** Explain Runge's phenomenon, quantify it via the Lebesgue
  constant, and describe why Chebyshev nodes are nearly optimal.
- **L7.4** Read a coefficient-magnitude plot and infer the regularity
  class of the underlying function.
- **L7.5** Build the Chebyshev differentiation matrix $D_N$, use it to
  approximate $f'$ at the nodes, and relate its conditioning to $N$.
- **L7.6** Apply Clenshaw–Curtis quadrature, derive its weights from the
  DCT, and predict its convergence rate from the smoothness of the
  integrand.
- **L7.7** Discretize a linear BVP via Chebyshev collocation, impose
  boundary conditions by row replacement, and explain the trade-off that
  motivates ultraspherical methods.
- **L7.8** Implement an adaptive stopping rule based on coefficient-tail
  decay, and recognize the failure mode for non-smooth $f$.
