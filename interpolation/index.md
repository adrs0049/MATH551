# Approximation Theory

:::{tip} Big Idea
Polynomial interpolation at Chebyshev nodes is the standard tool for
one-dimensional approximation. The textbook warning that "high-order
interpolation is unstable" is a statement about *equispaced nodes*, not
about polynomial interpolation as such. With Chebyshev nodes and the
barycentric evaluation formula, the interpolant converges to any
sufficiently smooth function at a rate determined by the smoothness of
the function itself.

The conceptual tool behind everything that follows is linear algebra on
function spaces. Polynomials of degree at most $n$ form an
$(n+1)$-dimensional vector space; a function is a vector, sampling is a
choice of basis, and differentiation and integration are a linear
operator and a linear functional on that space. Approximation theory
studies which basis to pick, where to put the nodes, and how the
smoothness of the underlying function controls the error. The remainder
of the chapter formalises this picture and develops the differentiation,
integration, root-finding, and boundary-value-problem solvers it
supports.
:::

## Why Polynomials?

The starting point is a 19th-century theorem that justifies everything
that follows.

:::{prf:theorem} Weierstrass Approximation Theorem
:label: thm-weierstrass

Let $f \in C([a,b])$ and $\varepsilon > 0$. There exists a polynomial $p$
such that
$$
\|f - p\|_\infty \;<\; \varepsilon.
$$
:::

So polynomials are *dense* in the space of continuous functions on a
bounded interval. Every continuous $f$ can be approximated arbitrarily
well, in the supremum norm, by some polynomial of some degree. This is
the existence statement that licenses the rest of the chapter. The harder
question, and the one we actually answer, is constructive: given $f$, how
do we compute a near-best polynomial approximation, and at what cost?

A subtlety from the same era warns us not to be naive about it. Faber and
Bernstein showed that there is *no fixed sequence of nodes*
$\{x_j^{(n)}\}$ for which the polynomial interpolant converges to $f$ for
every continuous $f$. Convergence depends on both the nodes and the
regularity of $f$. The two themes of this chapter follow from those two
ingredients: [where to put the nodes](point-choice.md), and how the
[smoothness of $f$](regularity-and-decay.md) controls how fast a truncated
Chebyshev series converges. The good news is that for the functions we
actually meet in applications, analytic, $C^r$ with bounded $r$-th
derivative, or just Lipschitz, the rates are excellent.

## What This Chapter Does

We start from polynomial interpolation as a change of basis between values
and coefficients, and identify which bases are numerically usable. We then
ask where to place the nodes, with the answer being Chebyshev points.
Coefficient decay reveals the smoothness of the underlying function, and
that same dictionary controls the accuracy of every downstream operation:
differentiation by the matrix $D_N$, integration by Clenshaw–Curtis
quadrature, root-finding via the colleague-matrix eigenvalue problem, and
the solution of linear boundary value problems by spectral collocation.
The chapter ends with adaptivity, where the coefficient tail itself tells
us when we have used enough basis functions.

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
- **L7.7** Compute the roots of a smooth function on $[-1, 1]$ by forming
  the Chebyshev colleague matrix and computing its eigenvalues, and
  explain why the colleague formulation is better conditioned than the
  monomial companion.
- **L7.8** Discretize a linear BVP via Chebyshev collocation, impose
  boundary conditions by row replacement, and explain the trade-off that
  motivates ultraspherical methods.
- **L7.9** Implement an adaptive stopping rule based on coefficient-tail
  decay, and recognize the failure mode for non-smooth $f$.
