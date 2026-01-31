# Matrix-Vector Products

:::{tip} Big Idea
The matrix-vector product $A\mathbf{x}$ is the fundamental operation in numerical linear algebra. Understanding it geometrically—as a linear transformation—is key to understanding everything that follows.
:::

## Definition

Given $A \in \mathbb{R}^{m \times n}$ and $\mathbf{x} \in \mathbb{R}^n$, the product $\mathbf{y} = A\mathbf{x}$ is:

$$
y_i = \sum_{j=1}^{n} a_{ij} x_j, \quad i = 1, \ldots, m
$$

**Cost:** $2mn$ floating-point operations (flops).

## Two Ways to View $A\mathbf{x}$

### Row View: Dot Products

Each component of $\mathbf{y}$ is a dot product of a row of $A$ with $\mathbf{x}$:

$$
\mathbf{y} = \begin{pmatrix} \mathbf{a}_1^T \cdot \mathbf{x} \\ \mathbf{a}_2^T \cdot \mathbf{x} \\ \vdots \\ \mathbf{a}_m^T \cdot \mathbf{x} \end{pmatrix}
$$

### Column View: Linear Combination

$\mathbf{y}$ is a linear combination of the columns of $A$:

$$
A\mathbf{x} = x_1 \mathbf{a}^{(1)} + x_2 \mathbf{a}^{(2)} + \cdots + x_n \mathbf{a}^{(n)}
$$

where $\mathbf{a}^{(j)}$ is the $j$-th column of $A$.

:::{admonition} Key Insight
:class: important
The column view is often more useful: the columns of $A$ span the **column space** (or **range**) of $A$, and $A\mathbf{x}$ lives in this space.
:::

## Matrix as Linear Transformation

A matrix $A$ defines a linear map $T: \mathbb{R}^n \to \mathbb{R}^m$ via $T(\mathbf{x}) = A\mathbf{x}$.

**Properties:**
- $T(\mathbf{x} + \mathbf{y}) = T(\mathbf{x}) + T(\mathbf{y})$
- $T(\alpha \mathbf{x}) = \alpha T(\mathbf{x})$

**Geometric actions of matrices:**

| Matrix Type | Geometric Effect |
|-------------|------------------|
| Diagonal | Scaling along axes |
| Orthogonal | Rotation/reflection |
| Upper triangular | Shearing |
| Symmetric | Scaling along eigenvector directions |

## Matrix-Matrix Products

For $A \in \mathbb{R}^{m \times p}$ and $B \in \mathbb{R}^{p \times n}$:

$$
(AB)_{ij} = \sum_{k=1}^{p} a_{ik} b_{kj}
$$

**Cost:** $2mpn$ flops.

**Key properties:**
- $(AB)C = A(BC)$ (associative)
- $AB \neq BA$ in general (not commutative!)
- $(AB)^T = B^T A^T$
- $(AB)^{-1} = B^{-1} A^{-1}$ (if both invertible)

## Sparse Matrices

Many applications produce matrices with mostly zero entries:

- **Finite difference discretizations** of PDEs
- **Cell-cell interaction networks** (cells only interact with neighbors)
- **Graph Laplacians** from connectivity data

:::{admonition} Example: Cell Friction Network
In cell migration models, the friction matrix $F$ has entry $F_{ij} \neq 0$ only if cells $i$ and $j$ are in contact. For $N$ cells with $k$ neighbors each, $F$ has only $O(kN)$ nonzeros instead of $O(N^2)$.
:::

For sparse matrices:
- Store only nonzero entries (CSR, CSC formats)
- Matrix-vector product costs $O(\text{nnz})$ instead of $O(mn)$
- Specialized algorithms preserve sparsity

## Computational Considerations

### Memory Access Patterns

Matrix storage order matters for performance:

- **Row-major (C/Python):** Store rows contiguously
- **Column-major (Fortran/MATLAB):** Store columns contiguously

Access data in storage order for cache efficiency!

### Level 2 BLAS

Matrix-vector operations are **Level 2 BLAS** (Basic Linear Algebra Subprograms):

```python
import numpy as np

# Matrix-vector product
y = A @ x          # Preferred syntax
y = np.dot(A, x)   # Also works

# In scipy
from scipy.linalg import blas
y = blas.dgemv(1.0, A, x)  # Direct BLAS call
```

Modern libraries (NumPy, SciPy) use optimized BLAS implementations automatically.
