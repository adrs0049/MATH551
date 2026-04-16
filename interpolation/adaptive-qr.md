---
kernelspec:
  name: python3
  display_name: Python 3
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/adaptive-qr.pdf
    id: interpolation-adaptive-qr-pdf
downloads:
  - id: interpolation-adaptive-qr-pdf
    title: Download PDF
---

# Adaptive QR: Choosing $N$ on the Fly

:::{tip} Big Idea
Don't guess how many Chebyshev coefficients $N$ you need to resolve a
solution. Build the QR factorisation of the discretised system one column
at a time, and stop as soon as the *residual* of the truncated solution
drops below tolerance. The residual at every step is read off for free
from the transformed right-hand side, and the almost-banded structure of
the ultraspherical operator from [§6](spectral-bvp.md) makes each Givens
update $O(1)$ work. The result is a single-pass solver that returns both
the solution and the smallest $N$ that resolves it.
:::

## The Stopping Idea

Why does watching the coefficient tail work? Because the
[regularity-decay theorem](#thm-algebraic-decay) says that for a smooth
$u$ the Chebyshev coefficients $\lvert c_k \rvert$ decay geometrically or
algebraically to round-off. Once the tail has flattened at the noise
floor, every higher coefficient contributes nothing more to
$\|u - p_n\|_\infty$. The smallest $N$ at which this happens is the
$N_{\mathrm{opt}}$ we want.

```{code-cell} python
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt
import scipy.fft as fft

def chebpts(n): return np.cos(np.pi * np.arange(n+1)/n)
def vals2coeffs(v):
    n = len(v) - 1
    c = fft.dct(v[::-1], type=1, norm='forward')
    c[1:n] *= 2
    return c

n = 1024; x = chebpts(n)
c_smooth = np.abs(vals2coeffs(np.exp(np.sin(5*x))))

tol = 1e-14
N_opt = next(k for k in range(len(c_smooth)) if c_smooth[k] < tol)

fig, ax = plt.subplots(figsize=(7, 4.2))
ax.semilogy(c_smooth + 1e-20, '.', ms=3, label=r'$|c_k|$ for $e^{\sin 5x}$')
ax.axhline(tol, color='k', ls=':', label='tolerance')
ax.axvline(N_opt, color='C3', ls='--',
           label=f'$N_{{\\mathrm{{opt}}}} = {N_opt}$')
ax.set_xlim(0, 80); ax.set_ylim(1e-18, 5)
ax.set_xlabel('$k$'); ax.set_ylabel(r'$|c_k|$')
ax.set_title('Stop where the tail crosses the tolerance')
ax.legend(); plt.tight_layout(); plt.show()
```

The rest of this section is the *systems analogue*. Instead of computing
the coefficient tail of a known function, we are solving the
ultraspherical system from [§6](spectral-bvp.md) for an *unknown*
$\mathbf{c}$. The residual plays the role of the tail, and the QR
factorisation gives it to us for free.

## The Discretised Problem

We use the same first-order BVP as in [§6](spectral-bvp.md):

$$
u'(x) = f(x), \qquad u(-1) = \alpha, \qquad x \in [-1, 1].
$$

After ultraspherical discretisation at size $N$, the problem is the
almost-banded linear system

$$
\mathcal{L}\, \mathbf{c} \;=\; \mathbf{r},
\qquad
\mathcal{L} =
\begin{pmatrix} \mathbf{b}^\top \\ \mathcal{D}_0 \end{pmatrix},
\qquad
\mathbf{r} = \begin{pmatrix} \alpha \\ \mathcal{S}_0 \mathbf{f}_T \end{pmatrix},
$$

with $\mathbf{c} = (c_0, \ldots, c_N)^\top$ the Chebyshev coefficients of
$u$. The first row is dense (the boundary condition); every other row is
banded with bandwidth one (the differentiation operator). The $N$ we want
is whatever makes the truncation error in $u$ smaller than the tolerance.

The naive approach picks $N$ in advance, solves once, hopes for the best,
and starts over with $2N$ if the answer is unresolved. Adaptive QR avoids
all of this.

## Building the QR Factorisation One Column at a Time

The plan is to factor $\mathcal{L} = QR$ *incrementally*. After processing
the first $j$ columns we will have

$$
Q_j^\top\, \mathcal{L}_{:, 1:j} \;=\;
\begin{pmatrix} R_j \\ 0 \end{pmatrix},
\qquad
Q_j^\top\, \mathbf{r} \;=\;
\begin{pmatrix} \tilde{\mathbf{r}}_{1:j} \\ \tilde{\mathbf{r}}_{j+1:\,\mathrm{end}} \end{pmatrix},
$$

with $R_j$ a $j \times j$ upper-triangular matrix and $Q_j$ orthogonal.
The bottom block of the transformed right-hand side records "what is left
over" after using $j$ columns of $\mathcal{L}$.

The mechanical tool we use is a **Givens rotation**, the $2 \times 2$
orthogonal matrix

$$
G(c, s) = \begin{pmatrix} c & s \\ -s & c \end{pmatrix},
\qquad c^2 + s^2 = 1,
\qquad
G \begin{pmatrix} a \\ b \end{pmatrix} = \begin{pmatrix} \sqrt{a^2 + b^2} \\ 0 \end{pmatrix}
$$

with $c = a/\sqrt{a^2+b^2}$, $s = b/\sqrt{a^2+b^2}$. Applied to a pair of
rows of a matrix it zeros one chosen entry. The same rotation must also be
applied to $\mathbf{r}$.

### Watching the factorisation grow

For a more interesting demo, take the variable-coefficient problem

$$
(1 + x)\, u'(x) + u(x) = f(x), \qquad u(-1) = \alpha,
$$

a slight generalisation of the BVP from [§6](spectral-bvp.md). Using the
multiplication operator from
[§6](spectral-bvp.md#multiplication-operator-for-variable-coefficients),
the discretised operator is

$$
\mathcal{L} \;=\;
\begin{pmatrix}
\mathbf{b}^\top \\[0.3ex]
M_1[1+x]\,\mathcal{D}_0 + \mathcal{S}_0
\end{pmatrix},
$$

with $\mathbf{b}_k = (-1)^k$ as before. It is almost banded: one dense top
row from the boundary condition, and bandwidth two below.

```{code-cell} python
:tags: [hide-input]

def M1_x(N):
    """Multiplication by x in the U_k = C^(1)_k basis, size (N+1)x(N+1).
    Identity (x U_k) = (U_{k+1} + U_{k-1})/2 with U_{-1} = 0."""
    M = np.zeros((N+1, N+1))
    for n in range(N+1):
        if n - 1 >= 0:    M[n, n-1] = 0.5
        if n + 1 <= N:    M[n, n+1] = 0.5
    return M

def D0(N):
    """T_k -> U_{k-1} coefficient map: (D_0 c)_{k-1} = k c_k."""
    D = np.zeros((N+1, N+1))
    for k in range(1, N+1):
        D[k-1, k] = k
    return D

def S0(N):
    """T -> U conversion, banded with bandwidth two."""
    S = np.zeros((N+1, N+1))
    S[0, 0] = 1.0
    if N >= 2: S[0, 2] = -0.5
    for k in range(1, N+1):
        S[k, k] = 0.5
        if k + 2 <= N: S[k, k+2] = -0.5
    return S

def L_variable_first_order(N):
    """Almost-banded operator for (1+x) u' + u = f, u(-1) = alpha."""
    M_a = np.eye(N+1) + M1_x(N)            # multiplication by (1 + x)
    interior = M_a @ D0(N) + S0(N)
    L = np.vstack([(-1.0)**np.arange(N+1)[None, :], interior[:N, :]])
    return L

def incremental_qr_snapshots(L, snap_at):
    A = L.copy(); N = A.shape[0] - 1
    snaps = {}
    for j in range(N+1):
        for i in range(N, j, -1):
            a, b = A[j, j], A[i, j]
            if abs(b) < 1e-15: continue
            r = np.hypot(a, b); c, s = a/r, b/r
            G = np.array([[c, s], [-s, c]])
            A[[j, i], :] = G @ A[[j, i], :]
        if (j+1) in snap_at:
            snaps[j+1] = A.copy()
    return snaps

N = 40
L = L_variable_first_order(N)
snap_at = [4, 12, 24, 40]
snaps = incremental_qr_snapshots(L, snap_at)

fig, axes = plt.subplots(1, len(snap_at) + 1, figsize=(13, 3.2))
axes[0].spy(np.abs(L) > 1e-14, markersize=3)
axes[0].set_title(r'$\mathcal{L}$ for $(1+x)u\prime+u=f$')
for ax, j in zip(axes[1:], snap_at):
    ax.spy(np.abs(snaps[j]) > 1e-10, markersize=3)
    ax.set_title(f'after $j = {j}$ columns')
    ax.axvline(j - 0.5, color='C3', lw=0.8, alpha=0.6)
plt.tight_layout(); plt.show()
```

Each successive panel shows the working matrix after $j$ columns have
been triangularised. The vertical red line marks the boundary between
"finished" columns to the left and "pending" columns to the right. Two
features stand out and they are the whole point.

First, the upper-left $j \times j$ block is upper-triangular. That is
$R_j$, our growing factor.

Second, *no fill-in escapes a narrow band*. The work needed to add the
next column is bounded by the small number of nonzero entries that survive
at its position. Compare to a dense matrix, where each column zeroing
fills in a new dense triangle; here the structure stays compact for every
$j$. This is what keeps the cost of each update bounded as $j$ grows.

Now combine this with the stopping criterion. The residual estimate from
[](#thm-adaptive-qr-residual) is checked after every column. The moment it
drops below tolerance at some $j = N_{\mathrm{opt}}$, we *stop processing
columns*. The remaining columns to the right of the red line in the final
panel are never touched: their coefficients $c_k$ for $k > N_{\mathrm{opt}}$
contribute nothing to the resolved solution, so triangularising them would
be wasted work. We back-solve $R_{N_{\mathrm{opt}}}\, \mathbf{c} =
\tilde{\mathbf{r}}_{1:N_{\mathrm{opt}}}$ and we are done. The
visualisation of the growing factorisation and the visualisation of the
shrinking residual are two views of the same idea: spectral coefficient
decay tells us when to halt the QR sweep.

## The Residual Comes for Free

The crucial observation is that we do not need to actually *solve* the
truncated system to know its residual.

:::{prf:theorem} Residual from the transformed right-hand side
:label: thm-adaptive-qr-residual

Let $\hat{\mathbf{c}} \in \mathbb{R}^j$ be the unique solution of
$R_j \hat{\mathbf{c}} = \tilde{\mathbf{r}}_{1:j}$, padded with zeros to
length $N+1$ so that $\hat{\mathbf{c}}_{j+1:N+1} = 0$. Then the
*least-squares residual* of $\hat{\mathbf{c}}$ in the system
$\mathcal{L} \mathbf{c} = \mathbf{r}$ has Euclidean norm

$$
\big\| \mathcal{L}\, \hat{\mathbf{c}} - \mathbf{r} \big\|_2
\;=\;
\big\| \tilde{\mathbf{r}}_{j+1:\,\mathrm{end}} \big\|_2.
$$
:::

:::{prf:proof}
:class: dropdown

$Q_j$ is orthogonal, so it preserves Euclidean norms:

$$
\big\| \mathcal{L}\, \hat{\mathbf{c}} - \mathbf{r} \big\|_2
= \big\| Q_j^\top \big( \mathcal{L} \hat{\mathbf{c}} - \mathbf{r} \big) \big\|_2.
$$

Compute the right-hand side. Since $\hat{\mathbf{c}}$ has only the first
$j$ components nonzero,
$Q_j^\top \mathcal{L} \hat{\mathbf{c}} = \begin{pmatrix} R_j \hat{\mathbf{c}} \\ 0 \end{pmatrix}$,
and $R_j \hat{\mathbf{c}} = \tilde{\mathbf{r}}_{1:j}$ by construction. So

$$
Q_j^\top\big( \mathcal{L} \hat{\mathbf{c}} - \mathbf{r} \big)
= \begin{pmatrix} \tilde{\mathbf{r}}_{1:j} \\ 0 \end{pmatrix}
- \begin{pmatrix} \tilde{\mathbf{r}}_{1:j} \\ \tilde{\mathbf{r}}_{j+1:\,\mathrm{end}} \end{pmatrix}
= \begin{pmatrix} 0 \\ -\tilde{\mathbf{r}}_{j+1:\,\mathrm{end}} \end{pmatrix},
$$

whose norm is $\|\tilde{\mathbf{r}}_{j+1:\,\mathrm{end}}\|_2$.
:::

So at every step we already have the exact least-squares residual norm.
No extra work is required.

## The Algorithm

The whole adaptive solver is now one loop:

1. Initialise $j = 0$, with the empty factorisation $R_0$ and the
   untouched right-hand side $\tilde{\mathbf{r}} = \mathbf{r}$.
2. For $j = 1, 2, 3, \ldots$:
   - Append column $j$ of $\mathcal{L}$ to $R$.
   - Apply Givens rotations to zero out the (one) sub-diagonal entry; apply
     the same rotations to $\tilde{\mathbf{r}}$.
   - Compute $\rho_j := \|\tilde{\mathbf{r}}_{j+1:\,\mathrm{end}}\|_2$.
   - If $\rho_j < \mathrm{tol}$, stop. Set $N_{\mathrm{opt}} = j$.
3. Solve $R_{N_{\mathrm{opt}}}\, \mathbf{c} = \tilde{\mathbf{r}}_{1:N_{\mathrm{opt}}}$
   by back substitution.

Because the per-column work is $O(1)$ and we stop as soon as the residual
is small enough, the whole thing runs in $O(N_{\mathrm{opt}})$ time. The
algorithm produces the smallest $N$ that resolves the solution to
tolerance, *and* the solution at that $N$, in one sweep.

## Worked Example

Take the same first-order problem as in [§6](spectral-bvp.md), with
$f(x) = \cos(8x) + \cos(80\,e^x)$ and $u(-1) = 0$. The right-hand side
contains a slow $\cos(8x)$ wave and a fast envelope $\cos(80 e^x)$. The
solution $u = \int_{-1}^{x} f$ inherits a similar mix of scales, so the
resolution required is not obvious by inspection.

```{code-cell} python
:tags: [hide-input]

def build_system(N):
    """Almost-banded ultraspherical operator for u' = f, u(-1) = alpha,
    truncated to (N+1) Chebyshev coefficients. Returns (L, S0)."""
    L = np.zeros((N+1, N+1))
    L[0, :] = (-1.0) ** np.arange(N+1)         # boundary row
    for k in range(1, N+1):
        L[k, k] = k                            # D_0
    S0 = np.zeros((N+1, N+1))                  # T -> U conversion
    S0[0, 0] = 1.0
    if N >= 2: S0[0, 2] = -0.5
    for k in range(1, N+1):
        S0[k, k] = 0.5
        if k + 2 <= N: S0[k, k+2] = -0.5
    return L, S0

def solve_at_size(N, f_T, alpha):
    L, S0 = build_system(N)
    rhs = np.zeros(N+1); rhs[0] = alpha
    rhs[1:] = (S0 @ f_T[:N+1])[:N]
    c = np.linalg.solve(L, rhs)
    # residual against the "exact" oversampled rhs gives convergence proxy
    return c, rhs, L

f = lambda x: np.cos(8*x) + np.cos(80 * np.exp(x))
N_ref = 512
f_T_ref = vals2coeffs(f(chebpts(N_ref)))
c_ref, _, _ = solve_at_size(N_ref, f_T_ref, alpha=0.0)

Ns = np.arange(4, 257, 4)
res_hist = []
for N in Ns:
    c_N, _, _ = solve_at_size(N, f_T_ref[:N+1], alpha=0.0)
    err = np.linalg.norm(c_ref[:N+1] - c_N)
    res_hist.append(err + 1e-300)

tol = 1e-12
N_opt = next((Ns[i] for i, r in enumerate(res_hist) if r < tol), Ns[-1])

fig, ax = plt.subplots(figsize=(7, 4.2))
ax.semilogy(Ns, res_hist, 'o-', ms=4)
ax.axhline(tol, color='k', ls=':', label=f'tolerance $10^{{-12}}$')
ax.axvline(N_opt, color='C3', ls='--', label=f'$N_{{\\mathrm{{opt}}}} \\approx {N_opt}$')
ax.set_xlabel('size $N$ at which we stop'); ax.set_ylabel('error in $\\mathbf{c}$')
ax.set_title(r"Truncation error vs $N$ for $u' = f,\ u(-1)=0$")
ax.legend(); plt.tight_layout(); plt.show()
```

Each point above is the error if we *had stopped* at that size $N$. The
adaptive QR algorithm produces this curve incrementally and stops as soon
as the (free) residual estimate drops below tolerance, returning both
$N_{\mathrm{opt}}$ and the solution at that size in a single sweep. No
oversize system is ever solved, and the boundary condition is enforced
exactly because it sits in the dense top row that is rotated through the
same QR process.

## Why It Works, In One Sentence

The almost-banded structure of [§6](spectral-bvp.md) keeps the per-column
QR cost at $O(1)$. The orthogonality of $Q_j$ converts the question
"how good is my truncated solution?" into a free read of one block of the
transformed right-hand side. The coefficient-decay theorem of
[§3](regularity-and-decay.md) guarantees that residual will drop at the
spectral rate determined by the smoothness of the data. Combine the
three and you get a solver that is fast, accurate, and *self-sizing*.

```{bibliography}
:filter: docname in docnames
```

:::{seealso}
- [Where to Place the Nodes](point-choice.md): Chebyshev basis and DCT.
- [Spectral Methods for BVPs](spectral-bvp.md): the almost-banded
  ultraspherical operator used here.
- [QR Factorization](../qr-least-squares/qr-factorization.md): the
  orthogonality property exploited by [](#thm-adaptive-qr-residual).
:::
