---
kernelspec:
  name: python3
  display_name: Python 3
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/lagrange.pdf
    id: interpolation-lagrange-pdf
downloads:
  - id: interpolation-lagrange-pdf
    title: Download PDF
---

# Polynomial Interpolation: Values and Coefficients

:::{tip} Big Idea
Polynomial interpolation is a *change of basis* between two equivalent
representations of the same polynomial: its values $\{f_j\}$ at $n+1$ nodes
and its coefficients $\{c_j\}$ in some basis of $\mathbb{P}_n$. The
Vandermonde matrix is that change-of-basis map. In the **monomial** basis
this map is catastrophically ill-conditioned. In the **Lagrange** basis it is
the identity, but naive evaluation is unstable. The **barycentric** formula
is the only practical way to evaluate a Lagrange interpolant.
:::

## The Interpolation Problem

Given $n+1$ distinct nodes $x_0 < x_1 < \cdots < x_n$ and values $f_0, \ldots, f_n$,
find a polynomial $p_n \in \mathbb{P}_n$ such that

$$
p_n(x_j) = f_j, \qquad j = 0, 1, \ldots, n.
$$

:::{prf:theorem} Existence and Uniqueness
:label: thm-interp-existence

For any $n+1$ distinct nodes there exists a unique $p_n \in \mathbb{P}_n$
interpolating the data.
:::

:::{prf:proof}
:class: dropdown

**Uniqueness.** If $p, q \in \mathbb{P}_n$ both interpolate the data then
$d = p - q \in \mathbb{P}_n$ has $n+1$ roots, so $d \equiv 0$.

**Existence.** The Lagrange formula below constructs one.
:::

## Three Bases for $\mathbb{P}_n$

Polynomials of degree $\le n$ form an $(n+1)$-dimensional vector space.
Picking a basis $\{\phi_0, \ldots, \phi_n\}$ writes every polynomial as
$p_n(x) = \sum_j c_j\, \phi_j(x)$. The interpolation conditions
$p_n(x_i) = f_i$ then become a linear system $\Phi\, \mathbf{c} =
\mathbf{f}$ where $\Phi_{ij} = \phi_j(x_i)$. The cost of interpolation,
its conditioning, and how easily we can extend it later all depend on
*which* basis we pick. We discuss three.

### Monomial basis

The basis you already know, $\{1, x, x^2, \ldots, x^n\}$. The motivation
is familiarity: every polynomial is written in this basis when you first
meet it. The coefficients $c_j$ are then the usual coefficients of $x^j$,
and the values↔coefficients system $\Phi\, \mathbf{c} = \mathbf{f}$ takes
the explicit form

$$
\underbrace{\begin{pmatrix}
1 & x_0 & x_0^2 & \cdots & x_0^n \\
1 & x_1 & x_1^2 & \cdots & x_1^n \\
\vdots & & & & \vdots \\
1 & x_n & x_n^2 & \cdots & x_n^n
\end{pmatrix}}_{V}
\begin{pmatrix} c_0 \\ c_1 \\ \vdots \\ c_n \end{pmatrix}
= \begin{pmatrix} f_0 \\ f_1 \\ \vdots \\ f_n \end{pmatrix}.
$$

$V$ is the **Vandermonde matrix**. Each row records the powers of one
node; each column records one monomial sampled at all nodes. The system
is solvable whenever the nodes are distinct, but solving it is
*numerically a disaster*: the condition number $\kappa(V)$ grows
exponentially in $n$.

```{code-cell} python
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt

ns = np.arange(4, 41, 2)
kappa_eq = []
kappa_ch = []
for n in ns:
    x_eq = np.linspace(-1, 1, n+1)
    x_ch = np.cos(np.pi * np.arange(n+1) / n)
    kappa_eq.append(np.linalg.cond(np.vander(x_eq, increasing=True)))
    kappa_ch.append(np.linalg.cond(np.vander(x_ch, increasing=True)))

fig, ax = plt.subplots(figsize=(7, 4.2))
ax.semilogy(ns, kappa_eq, 'o-', label='equispaced nodes')
ax.semilogy(ns, kappa_ch, 's-', label='Chebyshev nodes')
ax.axhline(1/np.finfo(float).eps, color='k', ls=':', alpha=0.6,
           label=r'$1/\varepsilon_{\mathrm{mach}}$')
ax.set_xlabel('$n$')
ax.set_ylabel(r'$\kappa(V)$')
ax.set_title('Vandermonde condition number')
ax.legend()
plt.tight_layout()
plt.show()
```

For equispaced nodes $\kappa(V)$ blows past $1/\varepsilon_{\text{mach}}$
already by $n \approx 25$. Solving $V\mathbf{c} = \mathbf{f}$ via a
*general* linear solver in double precision then loses essentially all
significant digits. Chebyshev nodes (anticipated in
[](point-choice.md)) help, but even then $\kappa(V)$ grows.

This blow-up reflects the choice of *algorithm* applied to the
values↔coefficients map, not the conditioning of polynomial interpolation
itself. Specialized solvers that exploit the Vandermonde structure, such
as the Björck–Pereyra algorithm {cite:p}`BjorckPereyra1970`, recover the
coefficients to far higher accuracy than $\kappa(V)$ alone would suggest.
A more robust strategy is to choose a basis in which the
values↔coefficients matrix $\Phi$ has *trivial* structure: triangular,
identity, or orthogonal. The conditioning then ceases to be a problem,
because $\kappa(\Phi) = 1$ for the identity and for any orthogonal matrix,
and triangular systems are solved by direct substitution. The Newton basis
(below) makes $\Phi$ lower triangular. The Lagrange basis after that makes
it the identity. The Chebyshev basis at Chebyshev nodes, introduced in
[](point-choice.md), makes $\Phi$ (a scaling of) the Discrete Cosine
Transform matrix, which is orthogonal up to that scaling. From the linear
algebra you have already developed in
[](../qr-least-squares/qr-factorization.md) we know that orthogonal
factors do not amplify error. That is exactly why these alternative bases
sidestep the Vandermonde blow-up: they replace the ill-conditioned dense
solve by a multiplication against a perfectly-conditioned matrix, which
the FFT moreover computes in $O(n \log n)$ operations.

### Newton basis

The motivation for Newton is *incrementality*: we want a basis where the
values↔coefficients matrix $\Phi$ is lower triangular, so that adding a
new data point appends one equation and one unknown without disturbing any
earlier coefficient. That forces $\phi_k$ to vanish at $x_0, \ldots,
x_{k-1}$, giving the basis

$$
\{1,\; (x-x_0),\; (x-x_0)(x-x_1),\; \ldots,\; (x-x_0)\cdots(x-x_{n-1})\}.
$$

The coefficients in this basis are the **divided differences** of $f$, and
the interpolant of degree $n+1$ is obtained from the interpolant of degree
$n$ by adding *one* new term. Newton is therefore the basis of choice when
the number of nodes is not known in advance, or when nodes are added
adaptively. We will not need this construction explicitly in what follows.
Its modern replacement is the adaptive evaluation discussed in
[](adaptive-qr.md), which gets the same incremental property in a
different way.

### Lagrange basis

The motivation here is hands-on: build the interpolant as an explicit
product of linear factors that automatically passes through the data.
For each node $x_j$, the polynomial

$$
\ell_j(x) = \prod_{i \ne j} \frac{x - x_i}{x_j - x_i}
$$

is, by construction, zero at every other node $x_i$ (one of the factors
in the product vanishes there) and equal to $1$ at $x_j$ (the numerator
matches the denominator). So $\ell_j(x_i) = \delta_{ij}$, and the linear
combination

$$
p_n(x) = \sum_{j=0}^n f_j\, \ell_j(x)
$$

automatically satisfies $p_n(x_i) = f_i$ at every node. No system to
solve, no inversion: just a product of monomials weighted by the data.
As a bonus, the values↔coefficients matrix $\Phi$ is now the identity, so
the coefficients of $p_n$ in this basis are literally the data values.

The cost is paid in *evaluation*: the formula above is $O(n^2)$ per
evaluation point and is prone to overflow as $n$ grows. The barycentric
reformulation below fixes both issues without leaving the basis.

```{code-cell} python
:tags: [hide-input]

n = 5
x_nodes = np.cos(np.pi * np.arange(n+1) / n)
xx = np.linspace(-1, 1, 400)

def lagrange_basis(j, x, nodes):
    out = np.ones_like(x)
    xj = nodes[j]
    for i, xi in enumerate(nodes):
        if i == j:
            continue
        out *= (x - xi) / (xj - xi)
    return out

fig, ax = plt.subplots(figsize=(7, 4))
for j in range(n+1):
    ax.plot(xx, lagrange_basis(j, xx, x_nodes), label=f'$\\ell_{j}$')
ax.plot(x_nodes, np.zeros_like(x_nodes), 'ko', ms=5)
ax.axhline(1, color='gray', lw=0.5)
ax.axhline(0, color='gray', lw=0.5)
ax.set_xlabel('$x$')
ax.set_ylim([-0.7, 1.1])
ax.set_title(f'Lagrange basis polynomials, $n={n}$')
ax.legend(ncol=3, fontsize=9, loc='lower center')
plt.tight_layout()
plt.show()
```

Each $\ell_j$ is the unique degree-$n$ polynomial that is $1$ at $x_j$ and
$0$ at every other node. The interpolant $\sum f_j \ell_j$ inherits the
required values automatically.

## The Barycentric Formula

We can reorganize the Lagrange formula algebraically into a form that
evaluates in $O(n)$ and is numerically stable. Define the **node polynomial**
and **barycentric weights**

$$
\ell(x) = \prod_{k=0}^n (x - x_k), \qquad
\lambda_j = \frac{1}{\prod_{k \ne j}(x_j - x_k)}.
$$

Then $\ell_j(x) = \ell(x)\, \lambda_j / (x - x_j)$. Substituting and using the
fact that the $\ell_j$ are a partition of unity ($\sum_j \ell_j(x) = 1$) gives
the **second barycentric formula**:

:::{prf:proposition} Barycentric Interpolation Formula
:label: prop-bary

$$
p_n(x) = \frac{\displaystyle\sum_{j=0}^n \frac{\lambda_j}{x - x_j}\, f_j}
              {\displaystyle\sum_{j=0}^n \frac{\lambda_j}{x - x_j}},
\qquad p_n(x_j) = f_j.
$$
:::

Once the weights $\lambda_j$ are computed once in $O(n^2)$, each evaluation
costs $O(n)$. For Chebyshev nodes (next section) the weights have a
trivial closed form $\lambda_j = (-1)^j \delta_j$ with $\delta_j = \tfrac12$
at the endpoints and $1$ otherwise. No precomputation at all.

### Naive Lagrange vs. Barycentric

```{code-cell} python
:tags: [hide-input]

def lagrange_naive(xeval, x, f):
    n = len(x)
    out = np.zeros_like(xeval)
    for j in range(n):
        Lj = np.ones_like(xeval)
        for i in range(n):
            if i != j:
                Lj *= (xeval - x[i]) / (x[j] - x[i])
        out += f[j] * Lj
    return out

def bary_weights_cheb(n):
    w = np.ones(n+1)
    w[0] = 0.5; w[-1] = 0.5
    w[1::2] *= -1
    return w

def bary(xeval, x, f, w):
    diff = xeval[:, None] - x[None, :]
    at_node = np.isclose(diff, 0.0)
    diff[at_node] = 1.0
    terms = w / diff
    out = (terms * f).sum(axis=1) / terms.sum(axis=1)
    r, c = np.where(at_node)
    out[r] = f[c]
    return out

ns = np.arange(10, 201, 10)
xe = np.linspace(-0.999, 0.999, 200)
true = np.exp(np.sin(5*xe))
err_naive, err_bary = [], []
for n in ns:
    xn = np.cos(np.pi * np.arange(n+1) / n)
    fn = np.exp(np.sin(5*xn))
    w = bary_weights_cheb(n)
    err_naive.append(np.max(np.abs(lagrange_naive(xe, xn, fn) - true)))
    err_bary.append(np.max(np.abs(bary(xe, xn, fn, w) - true)))

fig, ax = plt.subplots(figsize=(7, 4))
ax.semilogy(ns, err_naive, 'o-', label='naive Lagrange')
ax.semilogy(ns, err_bary, 's-', label='barycentric')
ax.axhline(np.finfo(float).eps, color='k', ls=':', alpha=0.6,
           label=r'$\varepsilon_{\mathrm{mach}}$')
ax.set_xlabel('$n$')
ax.set_ylabel(r'max error  $|p_n - f|$')
ax.set_title(r'Evaluating the Chebyshev interpolant of $e^{\sin 5x}$')
ax.legend()
plt.tight_layout()
plt.show()
```

The mathematical interpolant is the same; only the *evaluation algorithm*
differs. The naive formula loses digits as $n$ grows because the products
$\prod_i (x - x_i)$ overflow and underflow against each other. The
barycentric formula is backward stable; see {cite:t}`Higham2004` and
{cite:t}`BerrutTrefethen2004` for the analysis.

## The Interpolation Error

How well does $p_n$ approximate a function $f$ from which the values
$f_j = f(x_j)$ are sampled?

:::{prf:theorem} Interpolation Error Formula
:label: thm-interp-error

If $f \in C^{n+1}[a,b]$ and $p_n$ interpolates $f$ at $x_0, \ldots, x_n \in
[a,b]$, then for every $x \in [a,b]$ there exists $\xi \in [a,b]$ with

$$
f(x) - p_n(x) = \frac{f^{(n+1)}(\xi)}{(n+1)!} \prod_{j=0}^n (x - x_j).
$$
:::

:::{prf:proof}
:class: dropdown

Fix $x \notin \{x_j\}$ and define
$g(t) = f(t) - p_n(t) - K \prod_j (t - x_j)$ with $K$ chosen so $g(x) = 0$.
Then $g$ has $n+2$ zeros in $[a,b]$ (the $n+1$ nodes plus $x$). By Rolle's
theorem, $g^{(n+1)}$ has a zero $\xi$. Differentiating: $g^{(n+1)}(t) =
f^{(n+1)}(t) - K(n+1)!$, so $K = f^{(n+1)}(\xi)/(n+1)!$.
:::

Two factors control the error:

1. **Smoothness of $f$**, through $f^{(n+1)}(\xi)$.
2. **Node placement**, through the **node polynomial** $\omega(x) = \prod (x - x_j)$.

The node polynomial is the lever we control. Choosing the nodes to make
$\max_x |\omega(x)|$ small is the entire point of the next section: equispaced
nodes leave $|\omega|$ enormous near the endpoints (Runge), while Chebyshev
nodes minimize it.

## References

```{bibliography}
:filter: docname in docnames
```
