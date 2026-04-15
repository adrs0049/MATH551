---
kernelspec:
  name: python3
  display_name: Python 3
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/qr-factorization.pdf
    id: qr-least-squares-qr-factorization-pdf
downloads:
  - id: qr-least-squares-qr-factorization-pdf
    title: Download PDF
---

# QR Factorization

:::{tip} Big Idea
The QR factorization decomposes any matrix into an orthogonal matrix times an
upper triangular matrix. Householder reflections provide a numerically stable
way to compute it, the gold standard for numerical linear algebra.
:::

## The QR Factorization

:::{prf:definition} QR Factorization
:label: def-qr-factorization

For $A \in \mathbb{R}^{m \times n}$ with $m \geq n$, the **QR factorization**
is:

$$
A = QR
$$

where $Q \in \mathbb{R}^{m \times m}$ is orthogonal ($Q^T Q = I$) and
$R \in \mathbb{R}^{m \times n}$ is upper triangular.
:::

:::{prf:theorem} Existence and Uniqueness of QR
:label: thm-qr-existence

Every $A \in \mathbb{R}^{m \times n}$ with $m \geq n$ has a QR factorization.

If $A$ has full column rank, the reduced QR is **unique** (up to signs of
columns of $Q$) when we require $R$ to have positive diagonal entries.
:::

:::{prf:proof}
:class: dropdown

The [Gram-Schmidt algorithm](gram-schmidt.md#alg-classical-gs) constructs an
orthonormal basis $q_1, \ldots, q_n$ for the column space of $A$ and expresses
each column $a_j$ as a linear combination of $q_1, \ldots, q_j$:

$$
a_j = \sum_{i=1}^{j} r_{ij} q_i
$$

where $r_{ij} = \langle a_j, q_i \rangle$ for $i < j$ and $r_{jj} = \|v_j\|$.
In matrix form this is $A = QR$ with $Q$ orthogonal and $R$ upper triangular.

Uniqueness: if $A = Q_1 R_1 = Q_2 R_2$ with positive diagonals, then
$Q_2^T Q_1 = R_2 R_1^{-1}$ is both orthogonal and upper triangular, hence a
diagonal matrix with $\pm 1$ entries. The positive diagonal constraint forces
it to be $I$.
:::

## Computing QR: Householder Reflections

Gram-Schmidt constructs QR by orthogonalizing columns via subtraction, which
is numerically unstable (see [loss of orthogonality](gram-schmidt.md#ex-gs-nearly-parallel)).
The stable alternative is the same idea as column elimination in
[Gaussian elimination](../direct-methods/gaussian-elimination.md), with one
substitution: instead of eliminating below the diagonal by subtracting
multiples of rows (recorded as triangular factors $L_k$), we eliminate by
applying orthogonal reflectors $H_k$ that zero entries below the diagonal
of each column. The triangular factors of GE are cheap but can amplify
errors; the orthogonal reflectors of Householder preserve norms and do not.

Concretely, we triangularize $A$ by multiplying with orthogonal matrices:

$$
Q_n \cdots Q_2 Q_1 A = R \quad \Rightarrow \quad A = \underbrace{Q_1^T Q_2^T \cdots Q_n^T}_{Q} R
$$

Each $Q_k$ is a **Householder reflector** that zeros out entries below the
diagonal in column $k$.

### Householder Reflectors

:::{prf:definition} Householder Reflector
:label: def-householder

Given a unit vector $v$, the **Householder reflector** is:

$$
H = I - 2vv^T
$$

This matrix reflects any vector across the hyperplane perpendicular to $v$.
:::

:::{prf:proposition} Properties of Householder Reflectors
:label: prop-householder-properties

A Householder reflector $H = I - 2vv^T$ satisfies:

1. $H$ is symmetric: $H^T = H$
2. $H$ is orthogonal: $H^T H = I$
3. $H$ is an involution: $H^2 = I$
4. $Hx$ reflects $x$ across the hyperplane $\{y : v^T y = 0\}$
:::

:::{prf:proof}
:class: dropdown

(1) $(I - 2vv^T)^T = I - 2vv^T$ since $(vv^T)^T = vv^T$.

(2) $H^TH = H^2 = (I - 2vv^T)(I - 2vv^T) = I - 4vv^T + 4v(v^Tv)v^T = I$
since $v^Tv = 1$.

(3) Follows from (2) and $H^T = H$.

(4) Decompose $x = (v^Tx)v + (x - (v^Tx)v)$. The first term is the component
along $v$, the second is in the hyperplane. Then
$Hx = -(v^Tx)v + (x - (v^Tx)v)$: the component along $v$ is negated.
:::

```{code-cell} python
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Left: reflection of a single vector
ax = axes[0]
v = np.array([1, 1]) / np.sqrt(2)  # reflection normal
x = np.array([2, 0.5])
H = np.eye(2) - 2 * np.outer(v, v)
Hx = H @ x

# Hyperplane (perpendicular to v)
perp = np.array([-v[1], v[0]])
t = np.linspace(-2.5, 2.5, 100)
ax.plot(t * perp[0], t * perp[1], 'k-', lw=1, alpha=0.4, label='hyperplane')

# Vectors
ax.annotate('', xy=x, xytext=[0, 0],
            arrowprops=dict(arrowstyle='->', color='#1f77b4', lw=2))
ax.annotate('', xy=Hx, xytext=[0, 0],
            arrowprops=dict(arrowstyle='->', color='#d62728', lw=2))
ax.annotate('', xy=v * 1.5, xytext=[0, 0],
            arrowprops=dict(arrowstyle='->', color='gray', lw=1.5, ls='--'))

# Dashed line showing reflection
ax.plot([x[0], Hx[0]], [x[1], Hx[1]], 'k:', lw=1, alpha=0.5)

ax.text(x[0] + 0.1, x[1] + 0.1, r'$\mathbf{x}$', fontsize=13, color='#1f77b4')
ax.text(Hx[0] + 0.1, Hx[1] - 0.25, r'$H\mathbf{x}$', fontsize=13, color='#d62728')
ax.text(v[0] * 1.5 + 0.1, v[1] * 1.5 + 0.1, r'$\mathbf{v}$', fontsize=13, color='gray')

ax.set_xlim(-2.5, 2.5)
ax.set_ylim(-2.5, 2.5)
ax.set_aspect('equal')
ax.grid(True, alpha=0.3)
ax.set_title(r'Householder reflection: $H = I - 2\mathbf{v}\mathbf{v}^T$')

# Right: choosing v to map x to ||x||e1
ax = axes[1]
x2 = np.array([1.5, 2.0])
norm_x = np.linalg.norm(x2)
target = norm_x * np.array([1, 0])

# Householder vector
v2 = x2 - target
v2 = v2 / np.linalg.norm(v2)
H2 = np.eye(2) - 2 * np.outer(v2, v2)
Hx2 = H2 @ x2

# Hyperplane
perp2 = np.array([-v2[1], v2[0]])
ax.plot(t * perp2[0], t * perp2[1], 'k-', lw=1, alpha=0.4, label='hyperplane')

ax.annotate('', xy=x2, xytext=[0, 0],
            arrowprops=dict(arrowstyle='->', color='#1f77b4', lw=2))
ax.annotate('', xy=Hx2, xytext=[0, 0],
            arrowprops=dict(arrowstyle='->', color='#d62728', lw=2))

# Show that ||x|| = ||Hx|| (circle)
theta = np.linspace(0, 2 * np.pi, 100)
ax.plot(norm_x * np.cos(theta), norm_x * np.sin(theta), 'k:', lw=0.8, alpha=0.3)

ax.text(x2[0] + 0.1, x2[1] + 0.1, r'$\mathbf{x}$', fontsize=13, color='#1f77b4')
ax.text(Hx2[0] + 0.1, Hx2[1] - 0.25, r'$H\mathbf{x} = \|\mathbf{x}\|\mathbf{e}_1$',
        fontsize=12, color='#d62728')

ax.set_xlim(-1, 3)
ax.set_ylim(-2, 2.5)
ax.set_aspect('equal')
ax.grid(True, alpha=0.3)
ax.set_title(r'Choosing $\mathbf{v}$ to map $\mathbf{x} \to \|\mathbf{x}\|\mathbf{e}_1$')

plt.tight_layout()
plt.show()
```

### Computing the Householder Vector

To zero out everything below the first entry of a vector $a$, we need
$Ha = \|a\|_2 e_1$. Since $H$ is a reflection and preserves norms, the target
$\|a\|_2 e_1$ lies on the same sphere as $a$ (right panel above). The
Householder vector $v$ is the bisector of $a$ and the target:

$$
v = a + \text{sign}(a_1) \|a\|_2 e_1, \qquad v \leftarrow v / \|v\|_2
$$

The sign choice avoids cancellation when $a$ is nearly parallel to $e_1$.

### The Algorithm

:::{prf:algorithm} Householder QR
:label: alg-householder-qr

**Input:** $A \in \mathbb{R}^{m \times n}$

**Output:** Orthogonal $Q$ and upper triangular $R$ with $A = QR$

1. **for** $k = 1, 2, \ldots, n$:
2. $\qquad$ Let $a = A_{k:m, k}$ (the subcolumn below the diagonal)
3. $\qquad$ $v = a + \text{sign}(a_1)\|a\|_2 e_1$, $\quad v \leftarrow v/\|v\|$
4. $\qquad$ $A_{k:m, k:n} \leftarrow A_{k:m, k:n} - 2v(v^T A_{k:m, k:n})$
5. $Q = H_1 H_2 \cdots H_n$
:::

Note that step 4 applies the reflection only to the submatrix, and the
matrix-vector product $v^T A_{k:m,k:n}$ avoids forming $H$ explicitly.

For $A \in \mathbb{R}^{m \times n}$, the factorization costs
$2mn^2 - \tfrac{2}{3}n^3$ flops; forming $Q$ explicitly costs an additional
$2mn^2$, while applying $Q^T$ to a single vector costs only $4mn - 2n^2$. In
the square case $m = n$, the factorization cost reduces to $\tfrac{4}{3}n^3$
flops, twice the cost of LU factorization. For
least squares we usually do not need $Q$ explicitly, just $Q^T b$, which can be
computed by applying the Householder reflectors sequentially.

## Numerical Stability

We developed the framework of [forward and backward stability](condition-stability.md)
in the previous chapter. Recall that [forward stability](condition-stability.md#def-forward-stable)
is unattainable ([](#prop-forward-stability-impossible)), and the right standard
is [backward stability](condition-stability.md#def-backward-stable-linalg): the
computed result is the exact result for a slightly perturbed input. Householder
QR achieves this.

What does this mean for a *factorization*? The algorithm takes $A$ and returns
computed factors $\hat{Q}, \hat{R}$. These will not satisfy $\hat{Q}\hat{R} = A$
exactly. Backward stability asks for the next best thing: there exists a small
perturbation $\delta A$ such that $\hat{Q}\hat{R} = A + \delta A$ *exactly*. In
other words, the computed factors are the exact QR factorization of a nearby
matrix. The size of $\delta A$ is the backward error, and we want
$\|\delta A\|/\|A\| = O(\varepsilon_{\text{mach}})$.

:::{prf:theorem} Backward Stability of Householder QR
:label: thm-householder-stability

The computed factors $\hat{Q}, \hat{R}$ satisfy:

$$
\hat{Q}\hat{R} = A + \delta A \quad \text{where} \quad \frac{\|\delta A\|}{\|A\|} = O(\varepsilon_{\text{mach}})
$$

Moreover, $\hat{Q}$ is orthogonal to machine precision:
$\|\hat{Q}^T\hat{Q} - I\| = O(\varepsilon_{\text{mach}})$.
:::

:::{prf:proof}
:class: dropdown

(Sketch, following Higham, *Accuracy and Stability of Numerical Algorithms*, Ch. 19.)

The full proof tracks rounding errors through each Householder step. We outline
the three key ingredients and how they combine.

**Step 1: One Householder reflection is backward stable.**

At step $k$, we apply a Householder reflector $H_k = I - 2v_kv_k^T$ to the
current matrix. The key operation is $H_k B$ for a submatrix $B$. In floating
point, the computed result satisfies

$$
\text{fl}(H_k B) = (H_k + \delta H_k) B, \quad \|\delta H_k\| = O(\varepsilon_{\text{mach}})
$$

To see why $\delta H_k$ is small, look at how $H_k B = B - 2 v_k (v_k^T B)$
is actually computed. The two intermediate quantities are the row vector
$w^T = v_k^T B$ and the rank-one update $2 v_k w^T$. Apply the standard
floating-point error model:

- $w^T$ is a row of dot products. Each computed entry $\hat{w}_j$ satisfies
  $|\hat{w}_j - w_j| \le c\,\varepsilon\, \|v_k\|_2\, \|B e_j\|_2$ for a small
  constant $c$. Since $\|v_k\|_2 = 1$, this is $O(\varepsilon)\,\|B\|_2$.
- The subtraction $B - 2 v_k \hat{w}^T$ adds another $O(\varepsilon)\,\|B\|_2$
  per entry from the multiply-add.

Collecting both, the computed result equals $(H_k + \delta H_k) B$ with
$\|\delta H_k\|_2 = O(\varepsilon)$. The constant does **not** depend on
$\|B\|$ because $\|v_k\|_2 = 1$ and $H_k$ does not stretch any direction.

Contrast this with Gaussian elimination. The elementary matrix
$L_k = I + \ell_k e_k^T$ contains the multipliers
$\ell_k = a_{ik}/a_{kk}$, which can be arbitrarily large when $|a_{kk}|$ is
small. The analogous bound there carries a $\|L_k\|$ factor, which is
exactly the growth factor that pivoting tries to control. Orthogonal $H_k$
has $\|H_k\|_2 = 1$ unconditionally, so no such amplification appears.

**Step 2: The perturbations compose additively, not multiplicatively.**

After $n$ steps, the computed upper triangular factor $\hat{R}$ satisfies

$$
\hat{R} = (H_n + \delta H_n) \cdots (H_1 + \delta H_1) A.
$$

We claim the perturbed product equals $Q^T + E$ with
$\|E\|_2 = O(n\varepsilon_{\text{mach}})$, where $Q^T = H_n \cdots H_1$. The
key point is the $n$, not $2^n$.

Start with two factors. Writing $\hat{H}_k = H_k + \delta H_k$ with
$\|\delta H_k\|_2 = O(\varepsilon)$,

$$
\hat{H}_2 \hat{H}_1 = H_2 H_1 + H_2\, \delta H_1 + \delta H_2\, H_1 + \delta H_2\, \delta H_1.
$$

The last term is $O(\varepsilon^2)$ and we drop it. The two first-order terms
each contain a single $\delta H_k$ flanked by orthogonal factors. Since
$\|H_k\|_2 = 1$, the flanking does not change the size:
$\|H_2\, \delta H_1\|_2 = \|\delta H_1\|_2 = O(\varepsilon)$, and likewise
for the other term. The total error is bounded by
$\|\delta H_1\|_2 + \|\delta H_2\|_2 = O(\varepsilon)$, the **sum** of the
factor errors rather than a product.

If $H_1, H_2$ were not orthogonal the sandwich $H_2\, \delta H_1$ would
scale by $\|H_2\|_2$, and with $n$ factors the analogous bound carries a
$\prod \|H_j\|_2$ amplification. Orthogonality is exactly what kills this
multiplicative blow-up.

Extending to $n$ factors gives the same conclusion: each of the $n$
first-order terms is a single $\delta H_k$ flanked by orthogonal matrices,
contributing $O(\varepsilon)$, so the total is $O(n\varepsilon)$.

**Step 3: Rewriting as a backward error on $A$.**

The computed $\hat{Q}$ is implicitly defined as the inverse of the product
of computed reflectors:
$\hat{Q}^T := \hat{H}_n \cdots \hat{H}_1 = Q^T + E$. Transposing,
$\hat{Q} = Q + E^T$. Substituting into $\hat{Q}\hat{R}$ and using
$\hat{R} = \hat{Q}^T A = (Q^T + E)A$ from Step 2,

$$
\hat{Q}\hat{R} = (Q + E^T)(Q^T + E) A
= A + (QE + E^T Q^T) A + O(\varepsilon^2),
$$

where we used $QQ^T = I$. Setting
$\delta A = (QE + E^T Q^T) A$, the relative backward error satisfies

$$
\frac{\|\delta A\|}{\|A\|} \leq 2\|Q\|_2\, \|E\|_2 = O(n\varepsilon_{\text{mach}}),
$$

since $\|Q\|_2 = 1$.

**Takeaway.** The proof works because orthogonal matrices have three
properties that prevent error growth: unit norm ($\|Q\|_2 = 1$, no
amplification), norm preservation ($\|Qx\|_2 = \|x\|_2$, no cancellation),
and unit condition number ($\kappa_2(Q) = 1$, perturbations stay small under
composition).
:::

The orthogonality of the computed $\hat{Q}$ is the key difference from
Gram-Schmidt. For Gram-Schmidt, the loss of orthogonality is
$\|\hat{Q}^T\hat{Q} - I\| = O(\kappa(A) \cdot \varepsilon_{\text{mach}})$,
which can be catastrophic for ill-conditioned matrices. Householder maintains
$O(\varepsilon_{\text{mach}})$ regardless of conditioning.

## Why Orthogonal Transformations Are Ideal

The stability of Householder QR is not accidental. It reflects a deeper
principle: **errors under composition of orthogonal matrices add, they do
not multiply.**

For a single orthogonal $Q$ this is obvious. If $x$ carries an error $e$,
then

$$
\|Q(x+e) - Qx\|_2 = \|Qe\|_2 = \|e\|_2.
$$

The error is transported, never magnified. The composition property is the
consequence we cashed in during Step 2 of
{prf:ref}`thm-householder-stability`. The perturbed product expands to

$$
(H_n + \delta H_n) \cdots (H_1 + \delta H_1)
= H_n \cdots H_1
+ \sum_{k=1}^{n} H_n \cdots H_{k+1}\, \delta H_k\, H_{k-1} \cdots H_1
+ O(\varepsilon^2).
$$

For orthogonal $H_k$ each summand has norm $\|\delta H_k\|_2$ because the
flanking factors have unit norm, so the total first-order error is bounded
by

$$
\sum_{k=1}^{n} \|\delta H_k\|_2 = O(n\varepsilon).
$$

For general invertible factors $M_k$ the same expansion holds, but the
sandwich inequality gives

$$
\| M_n \cdots M_{k+1}\, \delta M_k\, M_{k-1} \cdots M_1 \|_2
\leq \|\delta M_k\|_2 \prod_{j \neq k} \|M_j\|_2,
$$

so the total error is bounded by

$$
\sum_{k=1}^{n} \|\delta M_k\|_2 \prod_{j \neq k} \|M_j\|_2
= O\!\left( \varepsilon \prod_{j=1}^{n} \|M_j\|_2 \right),
$$

multiplicative in the factor norms. If each $\|M_j\|_2 = \|M\|_2$, this is
$O(\varepsilon \|M\|_2^n)$, the exponential blow-up that orthogonality
avoids.

This is why orthogonal transformations are the tool of choice throughout
numerical linear algebra, from QR factorization to eigenvalue algorithms:
they are the only family of matrices for which long sequences of operations
do not amplify rounding errors.

## Solving Linear Systems with QR

Given a square system $A\mathbf{x} = \mathbf{b}$ with $A \in \mathbb{R}^{n \times n}$ invertible, the QR factorization reduces it to a triangular solve:

$$
A\mathbf{x} = \mathbf{b} \quad \Rightarrow \quad QR\mathbf{x} = \mathbf{b} \quad \Rightarrow \quad R\mathbf{x} = Q^T\mathbf{b}
$$

The steps are:

1. **Factor:** compute $A = QR$ via Householder ($\frac{4}{3}n^3$ flops for square $A$)
2. **Transform:** compute $\tilde{\mathbf{b}} = Q^T\mathbf{b}$ ($\mathcal{O}(n^2)$ flops)
3. **Solve:** back substitution on $R\mathbf{x} = \tilde{\mathbf{b}}$ ($\mathcal{O}(n^2)$ flops)

:::{prf:example} Solving a Linear System with QR
:label: ex-solve-qr
:class: dropdown

Consider:

$$
A = \begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}, \quad \mathbf{b} = \begin{pmatrix} 3 \\ 1 \end{pmatrix}
$$

**Step 1: QR factorization.** The first column has norm $\|a_1\| = \sqrt{2}$.
The Householder vector maps $a_1$ to $\sqrt{2}\,e_1$. After applying the
reflector to the full matrix:

$$
Q = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}, \quad
R = \begin{pmatrix} \sqrt{2} & 0 \\ 0 & \sqrt{2} \end{pmatrix}
$$

(In this case $A$ is already orthogonal up to scaling.)

**Step 2: Transform the right-hand side.**

$$
\tilde{\mathbf{b}} = Q^T\mathbf{b} = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}\begin{pmatrix} 3 \\ 1 \end{pmatrix} = \begin{pmatrix} 4/\sqrt{2} \\ 2/\sqrt{2} \end{pmatrix} = \begin{pmatrix} 2\sqrt{2} \\ \sqrt{2} \end{pmatrix}
$$

**Step 3: Back substitution.**

$$
R\mathbf{x} = \tilde{\mathbf{b}} \implies \begin{pmatrix} \sqrt{2} & 0 \\ 0 & \sqrt{2} \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 2\sqrt{2} \\ \sqrt{2} \end{pmatrix} \implies \mathbf{x} = \begin{pmatrix} 2 \\ 1 \end{pmatrix}
$$
:::

Compared to LU, the QR approach costs roughly twice as many flops
($\frac{4}{3}n^3$ for QR vs $\frac{2}{3}n^3$ for LU). For square systems where
stability is not a concern, LU is preferred. The real advantage of QR appears
for [least squares problems](least-squares.md), where the system is
overdetermined and LU does not apply.
