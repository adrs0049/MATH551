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

The chapter is in two parts.

**Part A: 1D, deterministic.** We start from polynomial interpolation
as a change of basis between values and coefficients, and identify
which bases are numerically usable. We then ask where to place the
nodes, with the answer being Chebyshev points. Coefficient decay
reveals the smoothness of the underlying function, and that same
dictionary controls the accuracy of every downstream operation:
differentiation by the matrix $D_N$, integration by Clenshaw–Curtis
quadrature, root-finding via the colleague-matrix eigenvalue problem,
and the solution of linear boundary value problems by spectral
collocation. Part A ends with adaptivity, where the coefficient tail
itself tells us when we have used enough basis functions.

**Part B: high-dimensional, Monte Carlo.** Tensor-product Chebyshev in
$d$ dimensions inherits the 1D rate but pays an $n^d$ cost for the
grid. The way out is to pick a *parameterised* basis (ridge functions
$\sigma(w \cdot x + b)$) and a different regularity notion (the
Barron norm, an L¹ moment of the Fourier transform). The
universal-approximation theorem is the modern Weierstrass for this
basis, and Barron's theorem is the analogous regularity-rate
dictionary, with rate $O(C_f/\sqrt n)$ in any dimension. Part B closes
with PyTorch examples that read both the density theorem and the
dimension-free rate off training curves.

## Learning Outcomes

After completing this chapter you should be able to:

- **L7.1** State existence and uniqueness of polynomial interpolation, and
  identify the Vandermonde map between values and coefficients in any
  basis.
- **L7.2** Use the barycentric formula to evaluate a Chebyshev interpolant
  in $O(n)$.
- **L7.3** Explain Runge's phenomenon and describe why Chebyshev nodes
  are nearly optimal.
- **L7.4** Read a coefficient-magnitude plot and infer the regularity
  class of the underlying function.
- **L7.5** Differentiate $f$ in the Chebyshev coefficient basis using
  the $T \to U$ recurrence between Chebyshev coefficients of $f$ and
  $f'$, and predict the accuracy from the coefficient decay.
- **L7.6** Integrate $f$ in the Chebyshev coefficient basis as a
  linear functional on $\{c_k\}$ (inner product against fixed weights),
  and predict the quadrature error from the coefficient decay.
- **L7.7** Discretize a linear BVP via Chebyshev collocation, impose
  boundary conditions by row replacement, and explain the trade-off
  that motivates ultraspherical methods.
- **L7.8** Implement an adaptive stopping rule based on coefficient-
  tail decay, and recognize the failure mode for non-smooth $f$.
- **L7.9** Quantify the curse of dimensionality for tensor-product
  Chebyshev, and identify which features of $f$ control how badly $d$
  shows up in the cost.

**Optional:**

- **L7.10** Distinguish the Chebyshev *projection* (truncated series,
  $L^2$-optimal) from the Chebyshev *interpolant* (sampled values,
  near-best in $L^\infty$), quantify the gap via the Lebesgue
  constant, and apply the near-best approximation theorem.
- **L7.11** Build the Chebyshev differentiation matrix $D_N$ in the
  *value* basis, use it to approximate $f'$ at the nodes, and relate
  its conditioning to $N$.
- **L7.12** Apply Clenshaw–Curtis quadrature in the *value* basis,
  derive its weights from the DCT, and predict its convergence rate
  from the smoothness of the integrand.
- **L7.13** Compute the roots of a smooth function on $[-1, 1]$ by
  forming the Chebyshev colleague matrix and computing its
  eigenvalues, and explain why the colleague formulation is better
  conditioned than the monomial companion.
- **L7.14** State the universal-approximation theorem and identify it
  as the Weierstrass theorem for the ridge-function basis.
- **L7.15** State Barron's theorem, identify the Barron norm as a
  spectral-side regularity notion, and explain why a Monte Carlo
  representation produces a dimension-free rate.
