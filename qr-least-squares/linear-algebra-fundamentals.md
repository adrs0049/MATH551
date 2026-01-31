# Linear Algebra Fundamentals

:::{tip} Big Idea
**Vectors** are the objects we compute with; **matrices** are linear functions that act on vectors. To analyze algorithms, we need to measure size: **norms** for vectors, **induced norms** for matrices. The **condition number** $\kappa(A) = \|A\|\|A^{-1}\|$ tells us how sensitive a linear system is to perturbations.
:::

## Vectors and Vector Spaces

The fundamental objects in numerical linear algebra are **vectors**—elements of a vector space.

:::{prf:definition} Vector Space
:label: def-vector-space

A **vector space** $V$ over $\mathbb{R}$ is a set with two operations:
- **Addition:** $\mathbf{x} + \mathbf{y} \in V$ for $\mathbf{x}, \mathbf{y} \in V$
- **Scalar multiplication:** $\alpha\mathbf{x} \in V$ for $\alpha \in \mathbb{R}$, $\mathbf{x} \in V$

satisfying the usual axioms (associativity, commutativity, distributivity, zero element, inverses).
:::

The canonical example is $\mathbb{R}^n$—column vectors with $n$ real components. But vector spaces are far more general:

:::::{ prf:example} Examples of Vector Spaces
:label: ex-vector-spaces

| Space | Elements | Dimension |
|-------|----------|-----------|
| $\mathbb{R}^n$ | Column vectors | $n$ |
| $\mathbb{C}^n$ | Complex vectors | $n$ (over $\mathbb{C}$) |
| $\mathcal{P}_n$ | Polynomials of degree $\leq n$ | $n+1$ |
| $C[a,b]$ | Continuous functions on $[a,b]$ | $\infty$ |
| $L^2[a,b]$ | Square-integrable functions | $\infty$ |

::::{dropdown} Why This Matters
The same linear algebra concepts—basis, dimension, linear maps, norms—apply to **all** these spaces. When you solve a PDE numerically, you're doing linear algebra in a function space. The finite-dimensional theory ($\mathbb{R}^n$) is the template for the infinite-dimensional theory (functional analysis).

This is why we emphasize the abstract structure: vectors are elements of a vector space, matrices are linear maps. The specifics of $\mathbb{R}^n$ are just one instance.
::::
:::::

## Vector Norms

To analyze errors and convergence, we need to measure the "size" of vectors.

:::{prf:definition} Norm
:label: def-norm

A **norm** $\|\cdot\|: V \to \mathbb{R}$ on a vector space $V$ satisfies:

1. $\|\mathbf{x}\| \geq 0$ with equality iff $\mathbf{x} = \mathbf{0}$ (positive definiteness)
2. $\|\alpha\mathbf{x}\| = |\alpha|\|\mathbf{x}\|$ (homogeneity)
3. $\|\mathbf{x} + \mathbf{y}\| \leq \|\mathbf{x}\| + \|\mathbf{y}\|$ (triangle inequality)
:::

A vector space equipped with a norm is called a **normed vector space**. If it's also complete (Cauchy sequences converge), it's a **Banach space**—the natural setting for analysis.

### The $p$-Norms on $\mathbb{R}^n$

$$
\|\mathbf{x}\|_p = \left(\sum_{i=1}^n |x_i|^p\right)^{1/p}
$$

| Name | Formula | Interpretation |
|------|---------|---------------|
| 1-norm | $\|\mathbf{x}\|_1 = \sum_i \lvert x_i \rvert$ | Manhattan distance |
| 2-norm | $\|\mathbf{x}\|_2 = \sqrt{\sum_i x_i^2}$ | Euclidean length |
| $\infty$-norm | $\|\mathbf{x}\|_\infty = \max_i \lvert x_i \rvert$ | Maximum component |

### Norm Equivalence

:::{prf:theorem} Norm Equivalence
:label: thm-norm-equivalence

All norms on $\mathbb{R}^n$ are **equivalent**: for any two norms $\|\cdot\|_a$ and $\|\cdot\|_b$, there exist constants $c, C > 0$ such that:

$$
c\|\mathbf{x}\|_a \leq \|\mathbf{x}\|_b \leq C\|\mathbf{x}\|_a \quad \text{for all } \mathbf{x}
$$
:::

This is a **finite-dimensional phenomenon**. In infinite dimensions (function spaces), different norms can give genuinely different notions of convergence—a key subtlety in PDE theory.

### Function Space Norms

The same idea extends to functions:

| Space | Norm | Formula |
|-------|------|---------|
| $C[a,b]$ | Supremum norm | $\|f\|_\infty = \max_{x \in [a,b]} \lvert f(x) \rvert$ |
| $L^2[a,b]$ | $L^2$ norm | $\|f\|_2 = \sqrt{\int_a^b \lvert f(x) \rvert^2 dx}$ |
| $L^p[a,b]$ | $L^p$ norm | $\|f\|_p = \left(\int_a^b \lvert f(x) \rvert^p dx\right)^{1/p}$ |

These are the continuous analogs of the discrete $p$-norms—sums become integrals.

## Matrices as Linear Maps

**Matrices** are linear functions between vector spaces. A matrix $A \in \mathbb{R}^{m \times n}$ defines:

$$
T_A: \mathbb{R}^n \to \mathbb{R}^m, \qquad T_A(\mathbf{x}) = A\mathbf{x}
$$

**Linearity** means:
- $T_A(\mathbf{x} + \mathbf{y}) = T_A(\mathbf{x}) + T_A(\mathbf{y})$
- $T_A(\alpha\mathbf{x}) = \alpha T_A(\mathbf{x})$

Every linear map $\mathbb{R}^n \to \mathbb{R}^m$ corresponds to a unique $m \times n$ matrix, and vice versa.

:::{admonition} The Bigger Picture
:class: note
In infinite dimensions, linear maps between function spaces are called **operators**. Differential operators ($d/dx$), integral operators ($\int K(x,y) f(y) dy$), and solution operators for PDEs are all linear maps—the infinite-dimensional analogs of matrices.
:::

### The Matrix-Vector Product

Given $A \in \mathbb{R}^{m \times n}$ and $\mathbf{x} \in \mathbb{R}^n$:

$$
(A\mathbf{x})_i = \sum_{j=1}^{n} a_{ij} x_j, \quad i = 1, \ldots, m
$$

**Cost:** $2mn$ floating-point operations.

**Two views:**

| Row View | Column View |
|----------|-------------|
| Each $(A\mathbf{x})_i$ is a dot product: $\mathbf{a}_i^T \cdot \mathbf{x}$ | $A\mathbf{x}$ is a linear combination: $\sum_j x_j \mathbf{a}^{(j)}$ |

The **column view** reveals that $A\mathbf{x}$ lives in the **column space** (range) of $A$.

### Geometric Interpretation

| Matrix Type | Geometric Effect |
|-------------|------------------|
| Diagonal | Scaling along coordinate axes |
| Orthogonal ($Q^TQ = I$) | Rotation and/or reflection |
| Symmetric | Scaling along eigenvector directions |

## Matrix Norms

Since matrices are linear maps, we measure their size by how much they "stretch" vectors.

:::{prf:definition} Induced (Operator) Norm
:label: def-operator-norm

$$
\|A\| = \max_{\mathbf{x} \neq \mathbf{0}} \frac{\|A\mathbf{x}\|}{\|\mathbf{x}\|} = \max_{\|\mathbf{x}\| = 1} \|A\mathbf{x}\|
$$

The maximum stretching factor over all unit vectors.
:::

This definition works for **any** linear map between normed spaces—it's how we measure operators in functional analysis too.

| Name | Formula | Computation |
|------|---------|-------------|
| 1-norm | $\|A\|_1 = \max_j \sum_i \lvert a_{ij} \rvert$ | Maximum column sum |
| $\infty$-norm | $\|A\|_\infty = \max_i \sum_j \lvert a_{ij} \rvert$ | Maximum row sum |
| 2-norm | $\|A\|_2 = \sigma_{\max}(A)$ | Largest singular value |

**Key properties:**
- $\|A\mathbf{x}\| \leq \|A\| \cdot \|\mathbf{x}\|$ — the defining inequality
- $\|AB\| \leq \|A\| \cdot \|B\|$ — submultiplicativity
- $\|I\| = 1$

:::{admonition} Note
:class: note
The 1-norm and $\infty$-norm are cheap (just sums). The 2-norm requires singular values—more expensive but geometrically natural.
:::

