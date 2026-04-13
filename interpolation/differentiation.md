---
kernelspec:
  name: python3
  display_name: Python 3
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/differentiation.pdf
    id: interpolation-differentiation-pdf
downloads:
  - id: interpolation-differentiation-pdf
    title: Download PDF
---

# Differentiation as a Linear Operator

:::{tip} Big Idea
Differentiation is *linear*. On the finite-dimensional space $\mathbb{P}_n$
it is therefore a single matrix, and we can choose to write that matrix
either in the **value basis** (acting on nodal values) or in the
**coefficient basis** (acting on Chebyshev coefficients). Both
representations differentiate the same polynomial. The value-space matrix
$D$ is dense and reads "differentiate the global Lagrange interpolant".
The coefficient-space matrix is banded and reads "shift Chebyshev indices
down by one with weight $2k$". Different operations are easier in
different representations, and we move freely between them through the
DCT.
:::

## Two Pictures of Differentiation

We have already met the local picture: [forward differences](../numerical-algorithms/finite-differences.md)
give $O(h)$ accuracy and centered differences give $O(h^2)$. Each
derivative estimate uses 2–3 nearby values; the rest of the data is
ignored.

The spectral picture is the opposite. The interpolant $p_n$ already
encodes all the data globally. Its derivative $p_n'$ is just another
polynomial in $\mathbb{P}_{n-1} \subset \mathbb{P}_n$, which we can in
turn represent by its values at the same Chebyshev nodes or by its
Chebyshev coefficients. So differentiation is a linear map from
$\mathbb{P}_n$ to itself, and the question is: what does its matrix look
like in each of the two bases?

A naming convention before we start. Throughout this section and the
next we use "value space" and "coefficient space" as shorthand for two
specific bases of $\mathbb{P}_n$:

- **Value space** is the **Lagrange basis** $\{\ell_j\}$. The coordinates
  of a polynomial in this basis are literally its nodal values
  $(f_0, \ldots, f_n)$, since $p_n(x) = \sum_j f_j\, \ell_j(x)$.
- **Coefficient space** is the **Chebyshev basis** $\{T_k\}$. The
  coordinates are the Chebyshev coefficients $(c_0, \ldots, c_n)$.

Both are bases of the same $(n+1)$-dimensional space $\mathbb{P}_n$. The
DCT from [§2](point-choice.md) is precisely the change-of-basis matrix
between them. The "value-space differentiation matrix" $D$ below
represents differentiation in the Lagrange basis; the "coefficient-space"
operator $\mathcal{D}$ represents the same map in the Chebyshev basis.

## Differentiation in the Value Basis

In the Lagrange (value) basis the polynomial is

$$
p_n(x) = \sum_{j=0}^n f_j\, \ell_j(x),
$$

so

$$
p_n'(x) = \sum_{j=0}^n f_j\, \ell_j'(x).
$$

Evaluating $p_n'$ at the same nodes $x_i$ gives

$$
p_n'(x_i) \;=\; \sum_{j=0}^n \ell_j'(x_i)\, f_j.
$$

This is a matrix-vector product. Define the **value-space differentiation
matrix** $D \in \mathbb{R}^{(n+1)\times(n+1)}$ by

:::{prf:definition} Value-space differentiation matrix
:label: def-diff-matrix

$$
D_{ij} \;=\; \ell_j'(x_i),
\qquad
\big(D\,\mathbf{f}\big)_i \;=\; p_n'(x_i).
$$
:::

The entries of $D$ are determined entirely by the *nodes*, not by $f$.
Compute them once and you can differentiate any function sampled on those
nodes by a single matrix-vector multiplication.

### Closed-form entries from the barycentric formula

Recall from [§1](lagrange.md) that the barycentric form of the Lagrange
basis is

$$
\ell_j(x) \;=\; \ell(x)\, \frac{\lambda_j}{x - x_j},
\qquad
\ell(x) = \prod_{k=0}^n (x - x_k),
\qquad
\lambda_j = \frac{1}{\prod_{k \ne j}(x_j - x_k)}.
$$

Differentiating $\ell_j$ by the product rule and evaluating at $x = x_i$
with $i \ne j$, the term containing $\ell(x_i) = 0$ drops out and one
obtains, after some algebra,

$$
\ell_j'(x_i) \;=\; \frac{\lambda_j}{\lambda_i}\,\frac{1}{x_i - x_j},
\qquad i \ne j.
$$

The diagonal entries follow from a separate trick: since the constant
function $1$ is its own interpolant, $\sum_j \ell_j(x) \equiv 1$, so
differentiating gives $\sum_j \ell_j'(x) \equiv 0$. Evaluated at $x =
x_i$, this says the rows of $D$ sum to zero:

$$
D_{ii} \;=\; -\sum_{k \ne i} D_{ik}.
$$

:::{prf:proposition} Entries of $D$
:label: prop-diff-entries

For any distinct nodes $x_0, \ldots, x_n$ with barycentric weights
$\lambda_j$,

$$
D_{ij} \;=\;
\begin{cases}
\dfrac{\lambda_j / \lambda_i}{x_i - x_j}, & i \ne j, \\[1ex]
\displaystyle -\sum_{k \ne i} D_{ik}, & i = j.
\end{cases}
$$
:::

For Chebyshev nodes $x_j = \cos(j\pi/n)$ the barycentric weights are
$\lambda_j = (-1)^j \delta_j$ with $\delta_0 = \delta_n = \tfrac12$ and
$\delta_j = 1$ otherwise, so the formula evaluates without any extra
work. The value-space matrix $D$ is dense: every output value $p_n'(x_i)$
depends on every input value $f_j$. That is the price of using a global
interpolant.

```{code-cell} python
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt

def cheb_diff_matrix(N):
    if N == 0:
        return np.zeros((1, 1)), np.zeros(1)
    x = np.cos(np.pi * np.arange(N+1) / N)
    c = np.ones(N+1); c[0] = 2; c[N] = 2; c[1::2] *= -1
    X = np.outer(x, np.ones(N+1))
    dX = X - X.T + np.eye(N+1)
    D = np.outer(c, 1/c) / dX
    D -= np.diag(D.sum(axis=1))
    return D, x

f  = lambda x: np.exp(np.sin(5*x))
fp = lambda x: 5*np.cos(5*x) * np.exp(np.sin(5*x))

fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))

D, _ = cheb_diff_matrix(16)
im = axes[0].imshow(np.abs(D), cmap='viridis',
                    norm=plt.matplotlib.colors.LogNorm())
axes[0].set_title(r'Value-space matrix $|D|$, $N = 16$')
plt.colorbar(im, ax=axes[0], fraction=0.046)

Ns = np.arange(4, 81, 2)
err_spec, err_fd2 = [], []
for N in Ns:
    D, x = cheb_diff_matrix(N)
    err_spec.append(np.max(np.abs(D @ f(x) - fp(x))))
    xu = np.linspace(-1, 1, N+1); h = xu[1] - xu[0]
    fu = f(xu); df = np.empty_like(fu)
    df[1:-1] = (fu[2:] - fu[:-2]) / (2*h)
    df[0]  = (fu[1] - fu[0]) / h
    df[-1] = (fu[-1] - fu[-2]) / h
    err_fd2.append(np.max(np.abs(df - fp(xu))))

axes[1].semilogy(Ns, err_spec, 'o-', label='Chebyshev $D$')
axes[1].semilogy(Ns, err_fd2, 's-', label='2nd-order centred FD')
axes[1].set_xlabel('$N$'); axes[1].set_ylabel(r"$\|D \mathbf{f} - f'\|_\infty$")
axes[1].set_title(r"Convergence on $e^{\sin 5x}$")
axes[1].legend()
plt.tight_layout(); plt.show()
```

The left panel is the structure of $D$ at $N = 16$: every entry is
nonzero, and the corner entries dominate (the $(0,0)$ and $(N,N)$ entries
equal $\pm (2N^2 + 1)/6$). Those large boundary entries are exactly what
allows boundary conditions to influence the interior in spectral BVP
solvers ([§6](spectral-bvp.md)).

The right panel is convergence on $f(x) = e^{\sin 5x}$, an entire
function. Centered differences drop algebraically as $N^{-2}$. Chebyshev
differentiation drops *geometrically*: it reaches machine precision near
$N \approx 30$ and stays there. By the
[regularity-decay](regularity-and-decay.md) dictionary, the
differentiation error inherits the convergence rate of $p_n \to f$, with
at most a single power of $N$ lost in differentiation. Same dictionary
also predicts that on a non-smooth $f$ both schemes drop to algebraic and
the spectral edge disappears.

## Differentiation in the Coefficient Basis

The same operator looks completely different when written in the
Chebyshev basis $\{T_k\}$. Suppose

$$
p_n(x) \;=\; \sum_{k=0}^n c_k\, T_k(x).
$$

We want the coefficients $c'_k$ of $p_n'$ in the *same* basis:

$$
p_n'(x) \;=\; \sum_{k=0}^{n-1} c'_k\, T_k(x).
$$

To get a recurrence, start from the trigonometric identity
$T_k(\cos\theta) = \cos(k\theta)$. Differentiating $T_{k+1} - T_{k-1} =
2T_k$ (which itself follows from the product-to-sum formulas) and using
$T_k'(x) = k\, U_{k-1}(x)$ for the Chebyshev polynomials of the second
kind, one derives the relation

$$
\frac{T_{k+1}'(x)}{k+1} - \frac{T_{k-1}'(x)}{k-1} \;=\; 2\,T_k(x).
$$

Reading this as an equation between the coefficient vectors of $T_k$ and
$T_k'$ gives a backward two-term recurrence for the derivative
coefficients.

:::{prf:proposition} Coefficient-space differentiation
:label: prop-coeff-diff

Given the Chebyshev coefficients $c_0, c_1, \ldots, c_n$ of $p_n$, the
coefficients $c'_0, c'_1, \ldots, c'_{n-1}$ of $p_n'$ are computed by

$$
c'_n = 0, \qquad c'_{n-1} = 2 n\, c_n,
\qquad
c'_{k-1} \;=\; c'_{k+1} + 2 k\, c_k \quad \text{for } k = n-1, n-2, \ldots, 1,
$$

and finally $c'_0$ is halved if one uses the convention with a doubled
$c_0$ entry. The recurrence is exact for $p_n \in \mathbb{P}_n$.
:::

This is **differentiation in the coefficient basis**. In matrix terms,
the operator is the strictly upper-triangular, *banded* matrix

$$
\mathcal{D} \;=\;
\begin{pmatrix}
0 & 1 & 0 & 3 & 0 & 5 & \cdots \\
0 & 0 & 4 & 0 & 8 & 0 & \cdots \\
0 & 0 & 0 & 6 & 0 & 10 & \cdots \\
  &   &   &   & \ddots &  &
\end{pmatrix},
$$

with non-trivial entries only on every other super-diagonal. Compare to
$D$ in the value basis: there every entry is nonzero. The same
operation, two bases, two completely different sparsity patterns.

```{code-cell} python
:tags: [hide-input]

import scipy.fft as fft

def chebpts(N): return np.cos(np.pi * np.arange(N+1)/N)

def vals2coeffs(v):
    n = len(v) - 1
    c = fft.dct(v[::-1], type=1, norm='forward')
    c[1:n] *= 2
    return c

def coeffs2vals(c):
    n = len(c) - 1
    cs = c.copy(); cs[1:n] /= 2
    v = fft.idct(cs, type=1, norm='forward')
    return v[::-1]

def coeff_diff_matrix(N):
    """Banded operator that maps Chebyshev coefficients of p_n
    to Chebyshev coefficients of p_n' (length N+1, last entry zero)."""
    M = np.zeros((N+1, N+1))
    for k in range(N, 0, -1):
        # c'_{k-1} = c'_{k+1} + 2k c_k
        M[k-1, k] += 2*k
        if k+1 <= N:
            M[k-1, :] += M[k+1, :]
    M[0, :] /= 2
    return M

N = 16
D_value,  _ = cheb_diff_matrix(N)
D_coeff     = coeff_diff_matrix(N)

# Verify both compute the same derivative
x = chebpts(N); v = f(x)
deriv_via_value  = D_value @ v
c = vals2coeffs(v)
cprime = D_coeff @ c
deriv_via_coeff = coeffs2vals(cprime)

fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
axes[0].spy(D_value, markersize=6)
axes[0].set_title(r'Value-space $D$ (dense), $N = 16$')
axes[1].spy(D_coeff, markersize=6)
axes[1].set_title(r'Coefficient-space $\mathcal{D}$ (banded)')
plt.tight_layout(); plt.show()

err = np.max(np.abs(deriv_via_value - deriv_via_coeff))
print(f'max |D f - C^(-1) D_coeff C f|  =  {err:.2e}')
```

Two basis choices, two completely different sparsity patterns, identical
answer to within rounding error. The value-space matrix has every entry
nonzero, so applying it costs $O(N^2)$. The coefficient-space matrix is
strictly upper-triangular with non-trivial entries only on every other
super-diagonal, so applying it costs $O(N)$. Wrapping it with two DCTs
(values $\to$ coefficients, then back) gives a total of $O(N \log N)$ for
the same derivative.

## Choose the Basis That Makes the Operation Easy

The bigger point of this section is one we will keep returning to.
Polynomial calculus offers us a *choice of representation* for every
operation, and we are free to pick whichever is cheapest. Sampling a
function and multiplying two polynomials pointwise are easy in the value
basis. Reading off smoothness from coefficient decay, truncating to
degree $n$, and (as we will see in [§5](integration.md)) integrating
against the Chebyshev weight are easy in the coefficient basis.
Translating between the two costs a single DCT, that is $O(N \log N)$,
so we move freely between them. The choice between value-space
collocation and coefficient-space (ultraspherical) discretisation in
[§6](spectral-bvp.md) is the same idea applied to BVPs, and it is what
determines whether the resulting linear system is dense or banded.
