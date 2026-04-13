---
kernelspec:
  name: python3
  display_name: Python 3
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/point-choice.pdf
    id: interpolation-point-choice-pdf
downloads:
  - id: interpolation-point-choice-pdf
    title: Download PDF
---

# Where to Place the Nodes

:::{tip} Big Idea
The interpolation error has two factors: $f^{(n+1)}(\xi)/(n+1)!$, which we
cannot control, and the **node polynomial** $\omega(x) = \prod (x - x_j)$,
which we can. Equispaced nodes leave $|\omega|$ exponentially large near the
endpoints (Runge's phenomenon). **Chebyshev nodes** flatten it. The Lebesgue
constant turns this into a sharp quantitative comparison: equispaced nodes
amplify errors by $\sim 2^n$, Chebyshev nodes by $\sim \log n$.
:::

## Runge's Phenomenon

The notorious cautionary example. Take the smooth function
$f(x) = 1/(1 + 25x^2)$ on $[-1,1]$ and interpolate it through equispaced nodes.

```{code-cell} python
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt

f = lambda x: 1/(1 + 25*x**2)

def bary_eval(xeval, x, fvals, w):
    diff = xeval[:, None] - x[None, :]
    at_node = np.isclose(diff, 0.0)
    diff[at_node] = 1.0
    terms = w / diff
    out = (terms * fvals).sum(axis=1) / terms.sum(axis=1)
    r, c = np.where(at_node)
    out[r] = fvals[c]
    return out

def equispaced_weights(x):
    n = len(x)
    w = np.ones(n)
    for j in range(n):
        for i in range(n):
            if i != j:
                w[j] /= (x[j] - x[i])
    return w

xx = np.linspace(-1, 1, 1000)
fig, ax = plt.subplots(figsize=(7.5, 4.2))
ax.plot(xx, f(xx), 'k', lw=2, label='$f(x) = 1/(1+25x^2)$')
for n, color in zip([8, 16, 24], ['C0', 'C1', 'C3']):
    xe = np.linspace(-1, 1, n+1)
    fe = f(xe)
    we = equispaced_weights(xe)
    pe = bary_eval(xx, xe, fe, we)
    ax.plot(xx, pe, color=color, lw=1.2, label=f'$p_{{{n}}}$ (equispaced)')
ax.set_ylim(-2, 2.5)
ax.set_xlabel('$x$'); ax.legend(fontsize=9, loc='upper right')
ax.set_title("Runge: equispaced interpolation diverges near $\\pm 1$")
plt.tight_layout(); plt.show()
```

The interpolant *passes through every node*, but between nodes near the
endpoints it oscillates with amplitude that grows without bound. The
function is analytic on the whole real line; the divergence is entirely a
property of the **node placement**.

## Diagnosis: The Node Polynomial

Recall the [error formula](#thm-interp-error):

$$
f(x) - p_n(x) = \frac{f^{(n+1)}(\xi)}{(n+1)!} \, \omega(x), \qquad
\omega(x) = \prod_{j=0}^n (x - x_j).
$$

The first factor depends on $f$ and is fixed. Everything we control sits in
$\omega(x)$. The minimax question is then:

> **Choose $\{x_j\}$ to minimize $\max_{x \in [-1,1]} |\omega(x)|$.**

Look at $|\omega|$ for two choices of $n+1$ nodes:

```{code-cell} python
:tags: [hide-input]

def omega(x, nodes):
    out = np.ones_like(x)
    for xj in nodes:
        out *= (x - xj)
    return out

xx = np.linspace(-1, 1, 2000)
fig, ax = plt.subplots(figsize=(7.5, 4.2))
for n, ls in zip([10, 20, 40], ['-', '--', ':']):
    x_eq = np.linspace(-1, 1, n+1)
    x_ch = np.cos(np.pi * np.arange(n+1) / n)
    w_eq = np.abs(omega(xx, x_eq))**(1.0/n)
    w_ch = np.abs(omega(xx, x_ch))**(1.0/n)
    ax.plot(xx, w_eq, color='C0', ls=ls, label=f'equispaced, $n={n}$')
    ax.plot(xx, w_ch, color='C1', ls=ls, label=f'Chebyshev, $n={n}$')
ax.set_xlabel('$x$'); ax.set_ylabel(r'$|\omega(x)|^{1/n}$')
ax.set_title('Node polynomial, normalised by $n$')
ax.legend(ncol=2, fontsize=9); plt.tight_layout(); plt.show()
```

The right normalisation is $|\omega|^{1/n}$, the *per-node geometric
mean*. As $n \to \infty$ the curves converge to a fixed *equilibrium
potential* of the limiting node distribution. For Chebyshev nodes the
limit is a constant ($\tfrac12$), so $|\omega(x)|$ is spatially uniform up
to a single global factor. For equispaced nodes the limit bulges in the
interior and drops at the endpoints, creating a *ratio* between interior
and endpoint values of $|\omega|$ that grows geometrically with $n$.

That spatial non-uniformity is what causes Runge. The error formula pairs
$\omega(x)$ with $f^{(n+1)}(\xi)/(n+1)!$. For functions with complex
singularities close to $[-1,1]$, that derivative factor itself grows
geometrically. With Chebyshev nodes the two effects balance and the
product is small everywhere. With equispaced nodes the derivative growth
combines with the endpoint amplification of $|\omega|$, and the error
diverges near $\pm 1$.

## Chebyshev Nodes

:::{prf:definition} Chebyshev Nodes (of the second kind)
:label: def-chebyshev-points

The **Chebyshev points** on $[-1, 1]$ are
$$
x_j = \cos\!\left(\frac{j\pi}{n}\right), \qquad j = 0, 1, \ldots, n.
$$
:::

### Geometric origin

Place $n+1$ equispaced points on the upper unit semicircle and project them
straight down to the $x$-axis. Equispaced angles $\theta_j = j\pi/n$, so
$x_j = \cos\theta_j$. That is the picture.

```{code-cell} python
:tags: [hide-input]

n = 12
theta = np.pi * np.arange(n+1) / n
x_ch = np.cos(theta)

fig, ax = plt.subplots(figsize=(7, 3.6))
phi = np.linspace(0, np.pi, 200)
ax.plot(np.cos(phi), np.sin(phi), 'k', lw=1)
ax.plot([-1, 1], [0, 0], 'k', lw=1)
ax.plot(x_ch, np.sin(theta), 'o', color='C0', ms=6, label='equispaced on circle')
for xc, yc in zip(x_ch, np.sin(theta)):
    ax.plot([xc, xc], [0, yc], color='C0', lw=0.6, alpha=0.5)
ax.plot(x_ch, np.zeros_like(x_ch), 's', color='C3', ms=6, label='Chebyshev nodes')
ax.set_aspect('equal'); ax.axis('off'); ax.legend(loc='upper right')
ax.set_title('Chebyshev nodes are projections of equispaced points on the circle')
plt.tight_layout(); plt.show()
```

The clustering near $\pm 1$ is exactly what flattens $|\omega|$ at the
endpoints, where equispaced nodes leave it largest.

### Chebyshev polynomials and the recurrence

The definition $T_k(x) = \cos(k \arccos x)$ does not look like a
polynomial. Yet using the trigonometric identity
$\cos((k{+}1)\theta) + \cos((k{-}1)\theta) = 2\cos\theta\,\cos(k\theta)$ and
substituting $x = \cos\theta$ gives the three-term recurrence

$$
T_0(x) = 1, \quad T_1(x) = x, \quad T_{k+1}(x) = 2x\,T_k(x) - T_{k-1}(x),
$$

which produces a sequence of *honest polynomials* in $x$. The first few
are

$$
T_0 = 1, \quad T_1 = x, \quad T_2 = 2x^2 - 1, \quad
T_3 = 4x^3 - 3x, \quad T_4 = 8x^4 - 8x^2 + 1.
$$

So the Chebyshev polynomials are just a re-encoding of the same
$(n+1)$-dimensional space $\mathbb{P}_n$ as the monomials $1, x, \ldots,
x^n$, related by an upper-triangular change of basis. The recurrence is
also what makes the **Clenshaw algorithm** for evaluating
$\sum_k c_k T_k(x)$ work in $O(n)$.

## The Chebyshev Series

So far the Chebyshev polynomials have appeared as a tool for placing nodes.
They are also a *basis* of $\mathbb{P}_n$, and more importantly an
orthogonal basis for an entire function space on $[-1,1]$.

:::{prf:definition} Chebyshev series
:label: def-chebyshev-series

The **Chebyshev series** of a Lipschitz function $f$ on $[-1,1]$ is

$$
f(x) = \sum_{k=0}^\infty c_k\, T_k(x),
\qquad
c_k = \frac{2}{\pi} \int_{-1}^{1} \frac{f(x)\, T_k(x)}{\sqrt{1 - x^2}}\, dx,
$$

with $c_0$ obtained from the same formula but with prefactor $1/\pi$. The
series converges uniformly on $[-1,1]$. The coefficients $c_k$ are the
**Chebyshev coefficients** of $f$.
:::

The integral defining $c_k$ is a Fourier cosine integral in disguise: under
$x = \cos\theta$,

$$
c_k = \frac{2}{\pi} \int_0^\pi f(\cos\theta) \cos(k\theta)\, d\theta.
$$

Truncating the series at index $n$ gives the *projection* of $f$ onto
$\mathbb{P}_n$, which is the best $L^2_w$ approximation with respect to the
Chebyshev weight $w(x) = (1 - x^2)^{-1/2}$. This is *not* the same
polynomial as the Chebyshev interpolant, but the next two lemmas show that
the two agree to spectral accuracy and the interpolant is what we actually
compute.

## Computing the Coefficients

We rarely have a closed form for $c_k$. Instead we sample $f$ at the
Chebyshev nodes and use a discrete formula. Two facts make this efficient:
discrete orthogonality of the $T_k$ at those nodes, and the cosine
structure that lets the FFT do the work.

### The values↔coefficients matrix

Just as the monomial basis produced the Vandermonde system
$V \mathbf{c} = \mathbf{f}$, the Chebyshev basis at the Chebyshev nodes
produces its own values↔coefficients map. With $x_j = \cos(j\pi/n)$ and
$T_k(x_j) = \cos(jk\pi/n)$, the equations $p_n(x_j) = f_j$ become

$$
\underbrace{\begin{pmatrix}
T_0(x_0) & T_1(x_0) & \cdots & T_n(x_0) \\
T_0(x_1) & T_1(x_1) & \cdots & T_n(x_1) \\
\vdots   &          &        & \vdots \\
T_0(x_n) & T_1(x_n) & \cdots & T_n(x_n)
\end{pmatrix}}_{\textstyle T \;=\; \big(\cos(jk\pi/n)\big)_{jk}}
\begin{pmatrix} c_0 \\ c_1 \\ \vdots \\ c_n \end{pmatrix}
\;=\;
\begin{pmatrix} f_0 \\ f_1 \\ \vdots \\ f_n \end{pmatrix}.
$$

Compare to the Vandermonde matrix in the previous section: same layout
(rows = nodes, columns = basis functions), but the entries are *cosines*
rather than monomial powers. This change of basis fixes the ill
conditioning. The matrix $T$ satisfies

$$
T^\top T \;=\; D, \qquad D = \mathrm{diag}(n, n/2, \ldots, n/2, n),
$$

which is the discrete orthogonality relation for cosines. So $D^{-1/2} T$
is genuinely orthogonal, and $\kappa_2(D^{-1/2} T) = 1$. The full matrix
$T$ then has condition number

$$
\kappa_2(T) \;=\; \sqrt{\kappa_2(D)} \;=\; \sqrt{2},
$$

independent of $n$. Contrast this with $\kappa_2(V)$ for the Vandermonde
matrix, which grew exponentially in $n$ and exceeded $1/\varepsilon_{\text{mach}}$
already by $n \approx 25$. Inverting $T$ therefore costs no linear solve at
all: multiplying by $D^{-1} T^\top$ recovers $\mathbf{c}$ in one shot, and
the FFT does that multiplication in $O(n \log n)$. The two lemmas below
make these statements precise.

:::{prf:lemma} Discrete coefficient formula
:label: lem-discrete-cheb

For $f$ sampled at the Chebyshev nodes $x_j = \cos(j\pi/n)$, $j = 0,
\ldots, n$, the **discrete Chebyshev coefficients**

$$
\tilde c_k \;=\; \frac{2}{n} \sum_{j=0}^{n}{}'' f(x_j)\, \cos\!\left(\frac{jk\pi}{n}\right),
\qquad k = 0, 1, \ldots, n,
$$

(with $\sum''$ denoting that the $j = 0$ and $j = n$ terms are halved) are
the unique coefficients of the Chebyshev interpolant
$p_n = \sum_{k=0}^n \tilde c_k T_k$ satisfying $p_n(x_j) = f(x_j)$ for all
$j$. If $f \in \mathbb{P}_n$ then $\tilde c_k = c_k$ exactly.
:::

:::{prf:proof}
:class: dropdown

The $T_k$ satisfy a discrete orthogonality relation at the Chebyshev nodes:

$$
\sum_{j=0}^{n}{}'' T_k(x_j)\, T_l(x_j) \;=\;
\begin{cases}
n,   & k = l \in \{0, n\}, \\
n/2, & k = l, \quad 0 < k < n, \\
0,   & k \ne l, \quad 0 \le k, l \le n.
\end{cases}
$$

This follows from $T_k(x_j) = \cos(jk\pi/n)$ and the standard discrete
orthogonality of cosines. Imposing $p_n(x_j) = f(x_j)$ and applying the
orthogonality relation isolates $\tilde c_k$ as stated. If $f \in
\mathbb{P}_n$ the interpolant *is* $f$, so $\tilde c_k = c_k$.
:::

For functions $f \notin \mathbb{P}_n$ the discrete coefficients $\tilde
c_k$ approximate the true $c_k$ to within an *aliasing error* that decays
at the same rate as $c_k$ itself, so for our purposes we treat the two
interchangeably and write $c_k$ for both.

:::{prf:lemma} DCT computation
:label: lem-dct

The discrete coefficients of [](#lem-discrete-cheb) are the output of the
type-I **Discrete Cosine Transform** of the value vector $(f(x_0),
\ldots, f(x_n))$. Equivalently, they are the result of solving
$T \mathbf{c} = \mathbf{f}$. The Cooley–Tukey **FFT** algorithm
{cite:p}`CooleyTukey1965` carries out this basis change in $O(n \log n)$
operations.
:::

:::{prf:proof}
:class: dropdown

The kernel $\cos(jk\pi/n)$ is precisely the DCT-I kernel; embedding the
cosine sum into a length-$2n$ DFT by even symmetric reflection turns it
into a real FFT. The FFT itself is a structured factorisation of the
orthogonal DFT matrix into $O(\log n)$ sparse factors (the
butterflies), each of which can be applied in $O(n)$ work. The total
cost is $O(n \log n)$. The "solve" here is not Gaussian elimination on
$T$; it is the application of an explicit factorisation of $T^{-1}$
that exploits the recursive structure of the cosine kernel.
:::

The combination of the two lemmas is the punchline: starting from $n+1$
nodal samples of a smooth $f$, one FFT call performs the basis change
from values to Chebyshev coefficients. Inverting the same transform
reconstructs the values. This $O(n^2) \to O(n \log n)$ acceleration is
the reason spectral methods are practical, and we use it freely in the
[differentiation](differentiation.md) and [integration](integration.md)
sections that follow.

## (Optional) Sharper Bounds: Minimax and Lebesgue Constants

The growth-factor analysis of the previous sections is enough for the rest
of this chapter. The two refinements below sharpen *what exactly* makes
Chebyshev nodes good. They are standard results worth knowing about, but
not needed downstream.

### Minimax property

The node polynomial $\omega$ for Chebyshev nodes is, up to a constant, the
Chebyshev polynomial $T_n$. The polynomial $T_n(x) = \cos(n \arccos x)$
oscillates between $\pm 1$ exactly $n+1$ times on $[-1,1]$. This is the
*equioscillation* property, and it turns out to be optimal.

:::{prf:theorem} Chebyshev minimax
:label: thm-cheb-minimax

Among all monic polynomials of degree $n$ on $[-1,1]$, the rescaled
Chebyshev polynomial $2^{1-n} T_n(x)$ has the smallest maximum value, equal
to $2^{1-n}$. Equivalently, the Chebyshev nodes minimize $\max |\omega(x)|$.
:::

:::{prf:proof}
:class: dropdown

Suppose some monic $p_n$ has $\max |p_n| < 2^{1-n}$. Then $q = 2^{1-n} T_n -
p_n$ is a polynomial of degree at most $n-1$ that takes the same sign as
$T_n$ at the $n+1$ extrema of $T_n$. So $q$ alternates sign $n+1$ times,
hence has $n$ roots, which is impossible for degree $n-1$.
:::

### Lebesgue constant

The minimax property tells us about $|\omega|$. To compare interpolation
schemes head-on we want a single number that says *how much worse than the
best polynomial approximation can interpolation be*.

:::{prf:definition} Lebesgue constant
:label: def-lebesgue

Let $\ell_j$ be the Lagrange basis for nodes $x_0, \ldots, x_n$. The
**Lebesgue function** is $\Lambda_n(x) = \sum_j |\ell_j(x)|$ and the
**Lebesgue constant** is $\Lambda_n = \max_x \Lambda_n(x)$.
:::

:::{prf:theorem} Near-best approximation
:label: thm-near-best

If $p_n$ interpolates $f$ at the nodes and $p_n^*$ is the best $L^\infty$
polynomial approximation of degree $n$, then
$$
\|f - p_n\|_\infty \le (1 + \Lambda_n)\,\|f - p_n^*\|_\infty.
$$
:::

So $\Lambda_n$ is the amplification factor between best approximation and
interpolation. Look at the Lebesgue function for the two node families:

```{code-cell} python
:tags: [hide-input]

def lebesgue_func(xeval, nodes):
    n = len(nodes)
    out = np.zeros_like(xeval)
    for j in range(n):
        Lj = np.ones_like(xeval)
        for i in range(n):
            if i != j:
                Lj *= (xeval - nodes[i]) / (nodes[j] - nodes[i])
        out += np.abs(Lj)
    return out

n = 20
x_eq = np.linspace(-1, 1, n+1)
x_ch = np.cos(np.pi * np.arange(n+1) / n)
xx = np.linspace(-1, 1, 4000)

fig, ax = plt.subplots(figsize=(7.5, 4.2))
ax.semilogy(xx, lebesgue_func(xx, x_eq), label='equispaced')
ax.semilogy(xx, lebesgue_func(xx, x_ch), label='Chebyshev')
ax.set_xlabel('$x$'); ax.set_ylabel(r'$\Lambda_n(x) = \sum_j|\ell_j(x)|$')
ax.set_title(f'Lebesgue function, $n = {n}$')
ax.legend(); plt.tight_layout(); plt.show()
```

For equispaced nodes the Lebesgue function spikes near $\pm 1$ with peaks
that grow exponentially in $n$. For Chebyshev nodes it stays $O(\log n)$
across the whole interval. Tracking the maximum vs $n$:

```{code-cell} python
:tags: [hide-input]

ns = np.arange(4, 41, 2)
L_eq, L_ch = [], []
for n in ns:
    xe = np.linspace(-1, 1, n+1)
    xc = np.cos(np.pi * np.arange(n+1) / n)
    xx = np.linspace(-1, 1, 4000)
    L_eq.append(lebesgue_func(xx, xe).max())
    L_ch.append(lebesgue_func(xx, xc).max())

fig, ax = plt.subplots(figsize=(7.5, 4.2))
ax.semilogy(ns, L_eq, 'o-', label='equispaced')
ax.semilogy(ns, L_ch, 's-', label='Chebyshev')
ax.semilogy(ns, 2/np.pi * np.log(ns) + 1, 'k:', label=r'$\frac{2}{\pi}\log n$')
ax.set_xlabel('$n$'); ax.set_ylabel(r'$\Lambda_n$')
ax.set_title('Lebesgue constant grows exponentially vs logarithmically')
ax.legend(); plt.tight_layout(); plt.show()
```

The asymptotics ({cite:t}`Trefethen2013`):

$$
\Lambda_n^{\text{equi}} \sim \frac{2^{n+1}}{e\,n\log n}, \qquad
\Lambda_n^{\text{Cheb}} = \frac{2}{\pi}\log(n+1) + O(1).
$$

The Chebyshev growth is, up to a constant, optimal: no interpolation scheme
can do better than $\frac{2}{\pi}\log n$.

```{bibliography}
:filter: docname in docnames
```
