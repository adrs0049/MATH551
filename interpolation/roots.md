---
kernelspec:
  name: python3
  display_name: Python 3
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/roots.pdf
    id: interpolation-roots-pdf
downloads:
  - id: interpolation-roots-pdf
    title: Download PDF
---

# Roots as Eigenvalues: The Colleague Matrix

:::{tip} Big Idea
Finding the roots of a smooth function $f$ on $[-1, 1]$ reduces to a
linear-algebra problem in the same value↔coefficient picture that
handled differentiation and integration. Interpolate $f$ at Chebyshev
nodes, read off the coefficients via the DCT, and form the **colleague
matrix**, a tridiagonal-plus-rank-one $n \times n$ matrix whose
eigenvalues are exactly the roots of the truncated Chebyshev series.
Root-finding is then one call to a dense eigensolver.
:::

## The Reduction

Given a degree-$n$ Chebyshev expansion

$$
p_n(x) = \sum_{k=0}^n c_k T_k(x), \qquad c_n \ne 0,
$$

we want its roots in $[-1, 1]$. The standard trick for polynomials in the
monomial basis is to form the **companion matrix**: its characteristic
polynomial is (up to sign) the original polynomial, so its eigenvalues
are the roots. The Chebyshev-basis analogue is the **colleague matrix**.

The construction is a compact encoding of the multiplication-by-$x$ map
on the quotient ring $\mathbb{P}_n / (p_n)$, written in the Chebyshev
basis. Two ingredients:

1. *Three-term recurrence.* From $T_{k+1}(x) = 2 x T_k(x) - T_{k-1}(x)$,
   $$
   x\, T_0(x) = T_1(x), \qquad
   x\, T_k(x) = \tfrac{1}{2}\big(T_{k-1}(x) + T_{k+1}(x)\big)\ \text{for}\ k \ge 1.
   $$
   Multiplication by $x$ is therefore tridiagonal on the basis
   $\{T_k\}$.

2. *Reduction modulo $p_n$.* The product $x T_{n-1}$ produces a $T_n$
   that sits outside the basis $\{T_0, \ldots, T_{n-1}\}$. But on the
   roots of $p_n$, $T_n$ can be written in terms of the lower-index
   polynomials: from $p_n(x) = 0$,
   $$
   T_n(x) \;=\; -\sum_{k=0}^{n-1} \frac{c_k}{c_n} T_k(x).
   $$
   Substituting into $x T_{n-1} = \tfrac{1}{2}(T_{n-2} + T_n)$ gives a
   rank-one correction concentrated in the last column of the
   multiplication matrix.

Collecting both ingredients:

:::{prf:definition} Colleague matrix
:label: def-colleague

For $p_n(x) = \sum_{k=0}^n c_k T_k(x)$ with $c_n \ne 0$, the **colleague
matrix** $A \in \mathbb{R}^{n \times n}$ is

$$
A \;=\;
\begin{pmatrix}
0    & \tfrac12 &          &          &          \\
1    & 0        & \tfrac12 &          &          \\
     & \tfrac12 & 0        & \tfrac12 &          \\
     &          & \ddots   & \ddots   & \ddots   \\
     &          &          & \tfrac12 & 0
\end{pmatrix}
\;-\;
\frac{1}{2 c_n}
\begin{pmatrix}
0      & \cdots & 0 & c_0     \\
0      & \cdots & 0 & c_1     \\
\vdots &        & \vdots & \vdots \\
0      & \cdots & 0 & c_{n-1}
\end{pmatrix}.
$$

The tridiagonal part encodes multiplication by $x$ on $\{T_0, \ldots,
T_{n-1}\}$; the rank-one correction in the last column enforces
reduction modulo $p_n$.
:::

:::{prf:theorem} Roots are eigenvalues
:label: thm-colleague

The eigenvalues of $A$ are exactly the roots of $p_n$, counted with
multiplicity.
:::

:::{prf:proof}
:class: dropdown

The matrix $A$ represents multiplication by $x$ as a linear map on
the quotient ring $\mathbb{P}_n / (p_n)$, using the basis $\{T_0, \ldots,
T_{n-1}\}$. The eigenvalues of multiplication-by-$x$ on $\mathbb{P}_n /
(p_n)$ are the zeros of $p_n$ (this is the Chinese Remainder Theorem
for polynomial rings, written concretely). So the characteristic
polynomial of $A$ equals $p_n$ up to a constant, and the eigenvalues of
$A$ are the roots of $p_n$.
:::

## Why Colleague Rather Than Companion

The same reduction applied in the monomial basis produces the standard
**companion matrix**, whose eigenvalues are again the roots of the
polynomial. The two formulations differ in their conditioning, not in
what they compute.

From [](lagrange.md) we already know that the monomial-basis
values↔coefficients map, the Vandermonde matrix $V$, has condition
number growing exponentially in $n$, so recovering the monomial
coefficients from nodal samples of a smooth $f$ already loses digits.
Passing those coefficients to a companion-matrix eigensolver compounds
the loss: monomial companion eigenvalue problems are themselves poorly
conditioned for polynomials with clustered or large-modulus roots.

The Chebyshev basis sidesteps both failures. The values↔coefficients
map is a DCT with $\kappa_2 = \sqrt{2}$ (see [](point-choice.md)), so
the $c_k$ are obtained to working precision. The colleague matrix is
tridiagonal plus rank one, with entries bounded by $O(1)$ when $|c_k| /
|c_n|$ is bounded, and its eigenvalue problem is well-conditioned for
roots in $[-1, 1]$.

## Algorithm

For a function $f$ given symbolically or as a black box:

:::{prf:algorithm} Chebyshev root-finding on $[-1, 1]$
:label: alg-cheb-roots

**Input.** Function $f$; tolerance $\mathrm{tol}$.

**Output.** Approximate roots of $f$ in $[-1, 1]$.

1. Resolve $f$ by Chebyshev interpolation: call
   [](#alg-adaptive-cheb) to obtain coefficients $(c_0, \ldots, c_n)$
   with $n$ adaptively chosen.
2. Form the colleague matrix $A$ from the $c_k$.
3. Compute the eigenvalues $\{\lambda_j\}$ of $A$ with a dense
   eigensolver (e.g.\ `numpy.linalg.eigvals`).
4. Filter: retain $\lambda_j$ with $|\mathrm{Im}\,\lambda_j| <
   \mathrm{tol}$ and $\mathrm{Re}\,\lambda_j \in [-1, 1]$.
5. (Optional) Polish each surviving root with one Newton step on $f$.
:::

Step 5 is usually unnecessary: the colleague eigenvalues are already
accurate to near machine precision for smooth $f$. It becomes useful
when $f$ has a root on the boundary $\pm 1$ or when very tight accuracy
is needed.

### Demo: $f(x) = \sin(5 x) - x^2$

Four roots on $[-1, 1]$. The colleague method finds all of them at
once; compare against `brentq` applied to each sign-change bracket.

```{code-cell} python
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt
import scipy.fft as fft
from scipy.optimize import brentq

def chebpts(N):
    return np.cos(np.pi * np.arange(N+1) / N)

def cheb_coeffs(v):
    n = len(v) - 1
    c = fft.dct(v, type=1, norm='forward')
    c[1:n] *= 2
    return c

def colleague_matrix(c):
    n = len(c) - 1
    A = np.zeros((n, n))
    for i in range(n - 1):
        A[i, i+1] = 0.5
        A[i+1, i] = 0.5
    A[1, 0] = 1.0
    A[:n, n-1] -= c[:n] / (2 * c[n])
    return A

def cheb_roots(f, tol=1e-10, Nmax=2048):
    N = 16
    while N <= Nmax:
        c = cheb_coeffs(f(chebpts(N)))
        if np.max(np.abs(c[int(0.9*N):])) < tol:
            break
        N *= 2
    else:
        raise RuntimeError("failed to resolve")
    while len(c) > 2 and abs(c[-1]) < tol:
        c = c[:-1]
    A = colleague_matrix(c)
    ev = np.linalg.eigvals(A)
    real_ev = ev[np.abs(ev.imag) < 1e-6].real
    real_ev = real_ev[(real_ev > -1.001) & (real_ev < 1.001)]
    return np.sort(np.clip(real_ev, -1.0, 1.0))

f = lambda x: np.sin(5*x) - x**2
roots_cheb = cheb_roots(f)

# brentq validation
xx_fine = np.linspace(-1, 1, 4000)
vals = f(xx_fine)
sc = np.where(np.diff(np.sign(vals)))[0]
roots_bq = np.array(sorted(brentq(f, xx_fine[i], xx_fine[i+1]) for i in sc))

xx = np.linspace(-1, 1, 1000)
fig, ax = plt.subplots(figsize=(7.5, 4.2))
ax.axhline(0, color='gray', lw=0.6)
ax.plot(xx, f(xx), 'k', lw=1, label=r'$f(x) = \sin(5x) - x^2$')
ax.plot(roots_cheb, f(roots_cheb), 'o', color='C3', ms=8,
        label='colleague eigenvalues')
ax.plot(roots_bq, f(roots_bq), 'x', color='C0', ms=8,
        label='brentq (validation)')
ax.set_xlabel('$x$'); ax.set_ylabel('$f(x)$')
ax.set_title('Roots of $\\sin(5x) - x^2$ on $[-1, 1]$')
ax.legend(fontsize=9); plt.tight_layout(); plt.show()

print(f'colleague roots:  {roots_cheb}')
print(f'brentq roots:     {roots_bq}')
print(f'max discrepancy:  {np.max(np.abs(roots_cheb - roots_bq)):.2e}')
```

The two methods agree to better than $10^{-10}$. The colleague call
returns all four roots from a single eigensolve; `brentq` required
four bracket hunts plus four bisections. For a handful of roots the
cost is comparable. For a function with many roots, or when the
bracketing intervals are not easy to produce, the colleague approach
wins outright.

## Extensions

- **Other intervals.** For $[a, b]$, apply the affine change of
  variable $x = \frac{b - a}{2}\, t + \frac{b + a}{2}$, compute roots
  in $t \in [-1, 1]$, and map back.
- **Many roots, deflation.** If $n$ is large, solving a dense
  $n \times n$ eigenproblem costs $O(n^3)$. The *recursive* algorithm in
  {cite:t}`Boyd2002` subdivides the interval until each subinterval
  holds few roots, turning the total cost into $O(n \log^2 n)$. This is
  the strategy used by the Chebfun `roots` routine
  {cite:p}`BattlesTrefethen2004`.
- **Complex roots.** The unfiltered eigenvalues of $A$ include complex
  roots of $p_n$, which approximate complex zeros of $f$ provided they
  lie inside the Bernstein ellipse where the Chebyshev expansion of $f$
  converges. See {cite:t}`Trefethen2013`, Chapter 18.

```{bibliography}
:filter: docname in docnames
```
