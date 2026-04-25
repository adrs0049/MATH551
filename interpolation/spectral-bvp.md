---
kernelspec:
  name: python3
  display_name: Python 3
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/spectral-bvp.pdf
    id: interpolation-spectral-bvp-pdf
downloads:
  - id: interpolation-spectral-bvp-pdf
    title: Download PDF
---

# Spectral Methods for Boundary Value Problems

:::{tip} Big Idea
With $D_N$ in hand from the [previous section](differentiation.md), a linear
BVP $\mathcal{L}u = f$ on $[-1, 1]$ becomes a finite linear system
$L_N \mathbf{u} = \mathbf{f}$. Boundary conditions are imposed by replacing
the appropriate rows. For smooth solutions the error decays geometrically
in $N$. The price is that $L_N$ is **dense and ill-conditioned**: the
*ultraspherical* method recasts the same problem in a basis where the
operator matrix is sparse and well-conditioned.
:::

## Collocation in One Picture

Take a linear two-point BVP

$$
a(x)\, u''(x) + b(x)\, u'(x) + c(x)\, u(x) = f(x),
\qquad u(-1) = \alpha,\ u(1) = \beta.
$$

Let $\mathbf{u} = (u(x_0), \ldots, u(x_N))^T$ be the *unknown* values at the
Chebyshev nodes. Differentiation at the nodes is just $D_N \mathbf{u}$ and
$D_N^2 \mathbf{u}$. The differential operator becomes the matrix

$$
L_N = \mathrm{diag}(a)\, D_N^2 + \mathrm{diag}(b)\, D_N + \mathrm{diag}(c).
$$

Demanding $L_N \mathbf{u} = \mathbf{f}$ at every interior node and replacing
the first and last rows with the boundary conditions ($u_0 = \beta$,
$u_N = \alpha$, with the convention $x_0 = 1$, $x_N = -1$) gives a square
linear system to solve.

That's all there is to *physical-space* spectral collocation.

## Worked Example

Solve $-u'' + (1 + x^2)\, u = f$ on $[-1, 1]$ with homogeneous Dirichlet
boundary conditions, manufacturing $f$ from $u_{\mathrm{exact}}(x) =
(1 - x^2)\, e^{\sin 5x}$.

```{code-cell} python
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la

def cheb_diff_matrix(N):
    x = np.cos(np.pi * np.arange(N+1) / N)
    c = np.ones(N+1); c[0] = 2; c[N] = 2; c[1::2] *= -1
    X = np.outer(x, np.ones(N+1)); dX = X - X.T + np.eye(N+1)
    D = np.outer(c, 1/c) / dX
    D -= np.diag(D.sum(axis=1))
    return D, x

def solve_bvp(a, b, c, f, alpha, beta, N):
    D, x = cheb_diff_matrix(N)
    L = np.diag(a(x)) @ D @ D + np.diag(b(x)) @ D + np.diag(c(x))
    rhs = f(x)
    L[0, :] = 0; L[0, 0] = 1; rhs[0] = beta
    L[N, :] = 0; L[N, N] = 1; rhs[N] = alpha
    return x, la.solve(L, rhs)

ue = lambda x: (1 - x**2) * np.exp(np.sin(5*x))
def rhs_f(x):
    s, c5 = np.sin(5*x), np.cos(5*x); E = np.exp(s)
    # u'  = -2x E + (1-x^2)(5 c5 E)
    # u'' = -2 E + (-2x)(5 c5 E) + (-2x)(5 c5 E) + (1-x^2)(-25 s E + 25 c5^2 E)
    upp = -2*E - 20*x*c5*E + (1-x**2)*((25*c5**2 - 25*s)*E)
    return -upp + (1 + x**2)*ue(x)

a = lambda x: -np.ones_like(x)
b = lambda x:  np.zeros_like(x)
c = lambda x:  1 + x**2

xx = np.linspace(-1, 1, 600)
x16, u16 = solve_bvp(a, b, c, rhs_f, 0, 0, N=16)

fig, (axL, axR) = plt.subplots(1, 2, figsize=(11, 4.2))
axL.plot(xx, ue(xx), 'k', lw=2, label='exact')
axL.plot(x16, u16, 'o', ms=5, label='spectral, $N=16$')
axL.set_title('BVP solution'); axL.legend(); axL.set_xlabel('$x$')

Ns = np.arange(8, 81, 4)
err = []
for N in Ns:
    x, u = solve_bvp(a, b, c, rhs_f, 0, 0, N)
    err.append(np.max(np.abs(u - ue(x))))
axR.semilogy(Ns, err, 'o-')
axR.set_xlabel('$N$'); axR.set_ylabel(r'$\|u_N - u\|_\infty$')
axR.set_title('Spectral convergence')
plt.tight_layout(); plt.show()
```

The error halves the digit count between $N = 8$ and $N = 32$. This is
geometric convergence, exactly the
[coefficient-decay](regularity-and-decay.md) result transported through the
linear solve.

## The Conditioning Problem

The same plot for finite differences would show $O(N^{-2})$ instead of
geometric. The reward is dramatic. The cost is two-fold:

1. **Density.** $L_N$ is dense, so each solve costs $O(N^3)$ instead of
   $O(N)$ for tridiagonal FD.
2. **Conditioning.** From [§5](differentiation.md), $\kappa(D_N^2) \sim N^4$.
   Backward stability of `solve` then bounds the achievable accuracy at
   $\kappa(L_N) \cdot \varepsilon_{\mathrm{mach}} \sim N^4 \cdot 10^{-16}$.
   Pushing $N$ much above $\sim 10^3$ in double precision is hopeless.

These are not small problems for stiff or high-order operators (think
$u^{(4)}$ in beam equations). They have driven a substantial body of
research, of which the cleanest answer is *ultraspherical* methods.

## Ultraspherical Spectral Methods (Sketch)

The standard collocation method writes everything in the *value* basis at
Chebyshev nodes. The same basis is used for $u$, $u'$, and $u''$, and that
is exactly why $D_N^2$ is dense and ill-conditioned. The
**ultraspherical** approach of {cite:t}`OlverTownsend2013` (see also
{cite:t}`Trefethen2013`) uses *different bases for different orders of
derivative*. The unknown $u$ is expanded in the Chebyshev basis $T_k$. Its
first derivative $u'$ lives naturally in the basis of Chebyshev
polynomials of the second kind $U_k = C_k^{(1)}$, its second derivative
in $C_k^{(2)}$, and so on. The price for doing this is *bookkeeping*:
operators that convert between adjacent bases. The reward is a
**banded** linear system.

### A worked example: $u'(x) = f(x)$, $u(-1) = \alpha$

To see how this plays out, consider the simplest non-trivial BVP, the
first-order initial-value problem

$$
u'(x) = f(x), \qquad u(-1) = \alpha,
$$

solved on $[-1, 1]$. Write the unknown in the Chebyshev basis,

$$
u(x) = \sum_{k=0}^{N} c_k\, T_k(x),
$$

so the unknowns are the coefficients $\mathbf{c} = (c_0, \ldots, c_N)$.

**Step 1: Differentiation operator.** From the identity $T_k'(x) = k\,
U_{k-1}(x)$ for $k \ge 1$,

$$
u'(x) \;=\; \sum_{k=1}^{N} k\, c_k\, U_{k-1}(x).
$$

So the *coefficients of $u'$ in the $U$ basis* are obtained from
$\mathbf{c}$ by the rectangular sparse matrix

$$
\mathcal{D}_0 \;=\;
\begin{pmatrix}
0 & 1 & 0 & 0 & \cdots & 0 \\
0 & 0 & 2 & 0 & \cdots & 0 \\
0 & 0 & 0 & 3 & \cdots & 0 \\
\vdots & & & & \ddots & \vdots \\
0 & 0 & 0 & 0 & \cdots & N
\end{pmatrix}
\;\in\; \mathbb{R}^{N \times (N+1)},
$$

a single super-diagonal of strictly increasing entries. Applying it costs
$O(N)$.

**Step 2: Conversion operator.** The right-hand side $f$ is given as a
Chebyshev expansion $f = \sum f_k\, T_k$, but the equation $u' = f$ now
asks us to compare two coefficient vectors *in different bases*: the
$U$-basis coefficients of $u'$ on the left, the $T$-basis coefficients
of $f$ on the right. The bridge is the identity

$$
T_0 = U_0, \qquad
T_1 = \tfrac{1}{2} U_1, \qquad
T_k = \tfrac{1}{2}\big(U_k - U_{k-2}\big) \text{ for } k \ge 2,
$$

so the conversion matrix from $T$- to $U$-coefficients is

$$
\mathcal{S}_0 \;=\;
\begin{pmatrix}
1 & 0   & -\tfrac{1}{2} & 0 & 0 & \cdots \\
0 & \tfrac{1}{2} & 0   & -\tfrac{1}{2} & 0 & \cdots \\
0 & 0   & \tfrac{1}{2} & 0   & -\tfrac{1}{2} & \cdots \\
\vdots & & & & \ddots & \\
\end{pmatrix},
$$

a tridiagonal matrix (bandwidth two: a main diagonal and a second
super-diagonal). It costs $O(N)$ to apply.

**Step 3: Boundary condition.** Since $T_k(-1) = (-1)^k$, the condition
$u(-1) = \alpha$ becomes one linear constraint on $\mathbf{c}$:

$$
\mathbf{b}^\top \mathbf{c} = \alpha,
\qquad
\mathbf{b} = (1,\, -1,\, 1,\, -1,\, \ldots,\, (-1)^N).
$$

This is a single *dense row*.

**Step 4: Assemble.** The full discretisation is

$$
\underbrace{\begin{pmatrix} \mathbf{b}^\top \\[0.4ex] \mathcal{D}_0 \end{pmatrix}}_{\mathcal{L} \in \mathbb{R}^{(N+1) \times (N+1)}}
\;\mathbf{c} \;=\;
\begin{pmatrix} \alpha \\[0.4ex] \mathcal{S}_0 \mathbf{f}_T \end{pmatrix},
$$

with $\mathbf{f}_T$ the Chebyshev coefficients of $f$. The operator
$\mathcal{L}$ has *one* dense row on top of a sparse banded block. This is
the **almost-banded** structure that makes ultraspherical solves
efficient: an $O(N)$ banded factorisation handles the bulk, with a
small dense correction for the boundary row.

```{code-cell} python
:tags: [hide-input]

N = 32

# Physical-space collocation D for u'
D_phys, _ = cheb_diff_matrix(N)

# Ultraspherical operator: dense boundary row on top of banded D_0
def ultra_first_order(N):
    L = np.zeros((N+1, N+1))
    L[0, :] = (-1.0)**np.arange(N+1)        # u(-1) = sum (-1)^k c_k
    for k in range(1, N+1):
        L[k, k] = k                          # D_0
    return L

L_ultra = ultra_first_order(N)

fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
axes[0].spy(D_phys, markersize=3)
axes[0].set_title(f'Physical-space $D$, dense, $N = {N}$')
axes[1].spy(L_ultra, markersize=6)
axes[1].set_title(r'Ultraspherical $\mathcal{L}$: 1 dense row + banded $\mathcal{D}_0$')
plt.tight_layout(); plt.show()
```

### Multiplication operator for variable coefficients

For a variable-coefficient problem $a(x)\, u'(x) + \cdots = f(x)$, we also
need a *multiplication operator* that takes the coefficient vector of $u'$
and returns the coefficient vector of $a u'$ in the same basis. Suppose
$a$ has the short Chebyshev expansion

$$
a(x) \;=\; \sum_{k=0}^d a_k\, T_k(x),
$$

with $d$ small. We need the matrix $M_1[a]$ that acts on $U$-basis
coefficients, since $u'$ lives there. The basic building block is the
multiplication-by-$x$ operator, derived from the recurrence
$x\, U_k = \tfrac{1}{2}(U_{k+1} + U_{k-1})$ (with $U_{-1} = 0$):

$$
M_1[x] \;=\;
\frac{1}{2}\!
\begin{pmatrix}
0 & 1 & & \\
1 & 0 & 1 & \\
  & 1 & 0 & 1 \\
  &   & \ddots & \ddots & \ddots
\end{pmatrix},
$$

a tridiagonal matrix. From this and the three-term Chebyshev recurrence,
$M_1[T_k]$ is built by polynomial recursion: $M_1[T_0] = I$,
$M_1[T_1] = M_1[x]$, and

$$
M_1[T_{k+1}] \;=\; 2\, M_1[x]\, M_1[T_k] - M_1[T_{k-1}].
$$

The result is a Toeplitz-plus-Hankel matrix with **bandwidth equal to the
degree $k$**. Therefore $M_1[a] = \sum_k a_k\, M_1[T_k]$ has bandwidth at
most $d$. Coefficients $a$ with short Chebyshev expansion give *banded*
multiplication operators. For example $a(x) = 1 + x$ gives the
tridiagonal $M_1[1+x] = I + M_1[x]$.

The full operator $M_1[a]\,\mathcal{D}_0 + \mathcal{S}_0$ for the equation
$a(x)\, u'(x) + u(x) = f(x)$ inherits the almost-banded pattern: a few
dense rows from boundary conditions, sparse banded everything else.
Higher-order problems work the same way, with $\mathcal{D}_\lambda$ and
$\mathcal{S}_\lambda$ for the next basis up. We use this construction in
[](adaptive-qr.md) to drive the visualisation of the QR sweep on a
genuinely variable-coefficient problem.

The pay-off is dramatic. Solving the almost-banded system by an
adapted-bandwidth QR factorisation costs $O(N)$. The conditioning of
$\mathcal{L}$ stays $O(1)$ in $N$ instead of $O(N^{2k})$ for a $k$-th
order problem. The cost is the basis bookkeeping above and the loss of a
direct nodal interpretation. For the moderate-$N$ collocation problems
of this chapter, the dense version is fine. For the
[adaptive QR](adaptive-qr.md) algorithm in the next section, the
almost-banded structure is essential.

```{bibliography}
:filter: docname in docnames
```
