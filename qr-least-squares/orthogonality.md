# Orthogonality and Projections

:::{tip} Big Idea
Orthogonality is the key to understanding least squares. The best approximation to a vector within a subspace is characterized by the residual being **orthogonal** to that subspace. This principle—formulated through **inner products**—extends far beyond $\mathbb{R}^n$ to function spaces (Hilbert spaces), where it underlies Fourier series, wavelets, and PDE theory.
:::

## Inner Product Spaces

The concept of "angle" and "orthogonality" generalizes beyond $\mathbb{R}^n$ through the abstraction of an **inner product**.

:::{prf:definition} Inner Product
:label: def-inner-product

An **inner product** on a vector space $V$ (over $\mathbb{R}$) is a function $\langle \cdot, \cdot \rangle: V \times V \to \mathbb{R}$ satisfying:

1. **Symmetry:** $\langle x, y \rangle = \langle y, x \rangle$
2. **Linearity:** $\langle \alpha x + \beta y, z \rangle = \alpha\langle x, z \rangle + \beta\langle y, z \rangle$
3. **Positive definiteness:** $\langle x, x \rangle \geq 0$, with equality iff $x = 0$

A vector space with an inner product is called an **inner product space**.
:::

Every inner product induces a norm: $\|x\| = \sqrt{\langle x, x \rangle}$.

:::{prf:example} Inner Products
:label: ex-inner-products
:class: dropdown

| Space | Inner Product | Induced Norm |
|-------|---------------|--------------|
| $\mathbb{R}^n$ | $\langle x, y \rangle = x^T y = \sum_i x_i y_i$ | Euclidean norm $\|x\|_2$ |
| $L^2[a,b]$ | $\langle f, g \rangle = \int_a^b f(x)g(x)\,dx$ | $L^2$ norm $\|f\|_2$ |
| $\ell^2$ (sequences) | $\langle x, y \rangle = \sum_{i=1}^\infty x_i y_i$ | $\ell^2$ norm |

The same theorems we prove for $\mathbb{R}^n$—Pythagorean theorem, best approximation, Gram-Schmidt—work in **any** inner product space. When the space is complete (Cauchy sequences converge), it's called a **Hilbert space**.

Examples of Hilbert spaces:
- $\mathbb{R}^n$ with the dot product
- $L^2[a,b]$ — the natural setting for Fourier series
- Sobolev spaces $H^k$ — the natural setting for PDEs

The finite-dimensional theory you learn here is the template for infinite-dimensional analysis.
:::

## Orthogonality

:::{prf:definition} Orthogonality
:label: def-orthogonal

Vectors $x, y$ in an inner product space are **orthogonal** (written $x \perp y$) if:

$$
\langle x, y \rangle = 0
$$
:::

In $\mathbb{R}^n$: $\langle x, y \rangle = \|x\|_2 \|y\|_2 \cos(\theta)$, so orthogonal means $\theta = \pm \pi/2$—**at right angles**.

## The Pythagorean Theorem (Generalized)

:::{prf:theorem} Pythagorean Theorem
:label: thm-pythagorean

Let $v, w$ be vectors in an inner product space. If $v \perp w$, then:

$$
\|v\|^2 + \|w\|^2 = \|v - w\|^2
$$
:::

:::{prf:proof}
:class: dropdown

Using only the algebraic properties of the inner product:

$$
\begin{align}
\|v - w\|^2 &= \langle v-w, v-w \rangle \\
&= \langle v, v \rangle - \langle v, w \rangle - \langle w, v \rangle + \langle w, w \rangle \\
&= \|v\|^2 + \|w\|^2 \quad \text{(since } \langle v, w \rangle = 0\text{)}
\end{align}
$$
:::

This 2D picture captures the high-dimensional truth:

```
        v
        |
        |
        +------w
```

## Subspaces

:::{prf:definition} Subspace
:label: def-subspace

A **subspace** $U$ of an inner product space $V$ is a subset that is closed under addition and scalar multiplication:
- If $u, v \in U$, then $u + v \in U$
- If $u \in U$ and $\alpha \in \mathbb{R}$, then $\alpha u \in U$
:::

**Examples in $\mathbb{R}^n$:**
- $\{(x, y, 0) : x, y \in \mathbb{R}\}$ is a subspace of $\mathbb{R}^3$ (the $xy$-plane)
- $\text{span}\{v_1, \ldots, v_k\} = \{\alpha_1 v_1 + \cdots + \alpha_k v_k\}$
- The **range** of a matrix: $\text{R}(A) = \{Ax : x \in \mathbb{R}^n\}$

**Examples in $L^2[a,b]$:**
- Polynomials of degree $\leq n$
- Trigonometric polynomials $\text{span}\{1, \cos x, \sin x, \cos 2x, \sin 2x, \ldots\}$

## Orthogonal Complements

:::{prf:definition} Orthogonal Complement
:label: def-orthogonal-complement

The **orthogonal complement** of a subspace $U$ is:

$$
U^\perp = \{v \in V : \langle u, v \rangle = 0 \text{ for all } u \in U\}
$$
:::

This is the set of all vectors orthogonal to everything in $U$.

**Example in $\mathbb{R}^3$:** If $U = \{(x, y, 0)\}$ (the $xy$-plane), then $U^\perp = \{(0, 0, z)\}$ (the $z$-axis).

## Orthogonal Projection

The **orthogonal projection** of $v$ onto a unit vector $u$ is:

$$
\text{proj}_u v = \langle v, u \rangle \, u
$$

This gives the component of $v$ in the direction of $u$.

For projection onto a subspace $U$ with orthonormal basis $\{u_1, \ldots, u_m\}$:

$$
\text{proj}_U v = \sum_{i=1}^m \langle v, u_i \rangle \, u_i
$$

## The Best Approximation Theorem

:::{prf:theorem} Best Approximation
:label: thm-best-approximation

Let $U$ be a subspace of an inner product space $V$ and $x \in V$. Then $z \in U$ is the **best approximation** to $x$ in $U$ (minimizing $\|x - u\|$ over all $u \in U$) if and only if:

$$
x - z \in U^\perp
$$

That is, the error vector is orthogonal to the subspace.
:::

:::{prf:proof}
:class: dropdown

Suppose $z \in U$ and $x - z \in U^\perp$. For any $u \in U$:

Since $z - u \in U$ (subspaces are closed under subtraction) and $x - z \perp U$:

$$
\langle x - z, z - u \rangle = 0
$$

By Pythagoras:

$$
\|x - z\|^2 + \|z - u\|^2 = \|x - u\|^2
$$

Since $\|z - u\|^2 \geq 0$, we have $\|x - z\| \leq \|x - u\|$ for all $u \in U$. ∎
:::

**Geometric Picture:**

```
        x
       /|
      / |  (x - z) ⟂ U
     /  |
    z---+---- U (subspace)
```

## Orthonormal Bases

:::{prf:definition} Orthonormal Basis
:label: def-orthonormal-basis

A basis $\{u_1, \ldots, u_m\}$ for a subspace $U$ is **orthonormal** if:

1. $\langle u_i, u_j \rangle = 0$ for all $i \neq j$ (orthogonal)
2. $\|u_i\| = 1$ for all $i$ (normalized)

Equivalently: $\langle u_i, u_j \rangle = \delta_{ij}$ where $\delta_{ij}$ is the Kronecker delta.
:::

**Why orthonormal bases are useful:**

1. **Coefficients are easy:** $v = \sum_i c_i u_i$ implies $c_i = \langle v, u_i \rangle$
2. **Projections are simple:** $\text{proj}_U v = \sum_i \langle v, u_i \rangle \, u_i$
3. **Condition number is 1:** No amplification of errors

:::{prf:remark} The Bigger Picture: Fourier Series
:label: rmk-fourier-connection

In $L^2[-\pi, \pi]$, the functions $\{1, \cos x, \sin x, \cos 2x, \sin 2x, \ldots\}$ form an orthogonal basis. The Fourier coefficients $a_n = \langle f, \cos(nx) \rangle$ are just projections onto basis elements—the same formula works in $\mathbb{R}^n$ and in function spaces!
:::

## Orthogonal Matrices

:::{prf:definition} Orthogonal Matrix
:label: def-orthogonal-matrix

A square matrix $Q \in \mathbb{R}^{n \times n}$ is **orthogonal** if its columns form an orthonormal basis for $\mathbb{R}^n$.

Equivalently: $Q^T Q = I$ (and thus $Q^{-1} = Q^T$).
:::

**Key properties of orthogonal matrices:**

| Property | Meaning |
|----------|---------|
| $Q^T Q = I$ | Columns are orthonormal |
| $Q Q^T = I$ | Rows are orthonormal |
| $Q^{-1} = Q^T$ | Inverse is just transpose |
| $\|Qx\|_2 = \|x\|_2$ | Preserves lengths (isometry) |
| $\kappa_2(Q) = 1$ | Perfect conditioning |
