# Triangular Systems and Matrix Factorizations

:::{tip} Big Idea
Some linear systems are cheap to solve: diagonal systems cost $\mathcal{O}(n)$,
triangular systems cost $\mathcal{O}(n^2)$. Since a general system costs
$\mathcal{O}(n^3)$, the strategy for solving $A\mathbf{x} = \mathbf{b}$ is to
**factorize** $A$ into triangular (or orthogonal) factors, then solve the
resulting easy systems.
:::

## Diagonal Systems

For a diagonal matrix $\mathcal{D}$:

$$
\mathcal{D}\mathbf{x} = \begin{pmatrix} d_1 & & \\ & \ddots & \\ & & d_n \end{pmatrix} \begin{pmatrix} x_1 \\ \vdots \\ x_n \end{pmatrix} = \begin{pmatrix} b_1 \\ \vdots \\ b_n \end{pmatrix}
$$

The solution is immediate: $x_i = b_i / d_i$ for all $i$.

**Cost:** $n$ divisions = $\mathcal{O}(n)$ flops.

## Upper Triangular Systems (Back Substitution)

For an upper triangular matrix $\mathcal{U}$:

$$
\mathcal{U}\mathbf{x} = \begin{pmatrix} u_{11} & u_{12} & \cdots & u_{1n} \\ & u_{22} & \cdots & u_{2n} \\ & & \ddots & \vdots \\ & & & u_{nn} \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ \vdots \\ x_n \end{pmatrix} = \begin{pmatrix} b_1 \\ b_2 \\ \vdots \\ b_n \end{pmatrix}
$$

The last equation gives $x_n$ directly:
$$
u_{nn} x_n = b_n \implies x_n = \frac{b_n}{u_{nn}}
$$

Then work upward:
$$
x_i = \frac{b_i - \sum_{j=i+1}^{n} u_{ij} x_j}{u_{ii}}, \quad i = n-1, n-2, \ldots, 1
$$

:::{prf:algorithm} Back Substitution
:label: alg-back-substitution

**Input:** Upper triangular $\mathcal{U} \in \mathbb{R}^{n \times n}$,
right-hand side $\mathbf{b} \in \mathbb{R}^n$

**Output:** Solution $\mathbf{x}$ of $\mathcal{U}\mathbf{x} = \mathbf{b}$

1. $x_n \gets b_n / u_{nn}$
2. **for** $i = n-1, n-2, \ldots, 1$:
3. $\qquad x_i \gets \left(b_i - \sum_{j=i+1}^{n} u_{ij} x_j\right) / u_{ii}$
:::

**Cost:** $\mathcal{O}(n^2)$ flops.

## Lower Triangular Systems (Forward Substitution)

For a lower triangular matrix $\mathcal{L}$, the first equation gives $x_1$:
$$
l_{11} x_1 = b_1 \implies x_1 = \frac{b_1}{l_{11}}
$$

Then work downward:
$$
x_i = \frac{b_i - \sum_{j=1}^{i-1} l_{ij} x_j}{l_{ii}}, \quad i = 2, 3, \ldots, n
$$

:::{prf:algorithm} Forward Substitution
:label: alg-forward-substitution

**Input:** Lower triangular $\mathcal{L} \in \mathbb{R}^{n \times n}$,
right-hand side $\mathbf{b} \in \mathbb{R}^n$

**Output:** Solution $\mathbf{x}$ of $\mathcal{L}\mathbf{x} = \mathbf{b}$

1. $x_1 \gets b_1 / l_{11}$
2. **for** $i = 2, 3, \ldots, n$:
3. $\qquad x_i \gets \left(b_i - \sum_{j=1}^{i-1} l_{ij} x_j\right) / l_{ii}$
:::

**Cost:** $\mathcal{O}(n^2)$ flops.

:::{prf:example} Back Substitution
:label: ex-back-substitution
:class: dropdown

Solve $\mathcal{U}\mathbf{x} = \mathbf{b}$ where:

$$
\begin{pmatrix} 2 & 1 & 3 \\ 0 & 2 & 1 \\ 0 & 0 & 4 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} 1 \\ -1 \\ 4 \end{pmatrix}
$$

Working from the bottom:
- $x_3 = 4/4 = 1$
- $x_2 = (-1 - 1 \cdot 1)/2 = -1$
- $x_1 = (1 - 1 \cdot (-1) - 3 \cdot 1)/2 = -1/2$

Solution: $\mathbf{x} = (-1/2, -1, 1)^T$.
:::

## Why This Matters: Matrix Factorizations

The key observation: triangular systems cost $\mathcal{O}(n^2)$ to solve, while
a general system costs $\mathcal{O}(n^3)$. This motivates the central strategy
of numerical linear algebra: **factorize** the matrix into simpler pieces.

The naive approach to solving $A\mathbf{x} = \mathbf{b}$ is to compute
$A^{-1}$ and then $\mathbf{x} = A^{-1}\mathbf{b}$. This is a bad idea:
computing $A^{-1}$ costs $\mathcal{O}(n^3)$ operations, is numerically less
stable than alternatives, and does not exploit structure. Instead, we write $A$
as a product of matrices that are easy to solve with.

### The Main Factorizations

**LU factorization** writes $A = PLU$ where $P$ is a permutation matrix (row
swaps), $L$ is lower triangular (with ones on the diagonal), and $U$ is upper
triangular. You have seen this before: it is
[Gaussian elimination](https://www.buttenschoen.ca/MATH545/systems/lu.html) (Math 235)
organized as a matrix product. Solving $A\mathbf{x} = \mathbf{b}$
becomes $PLU \mathbf{x} = \mathbf{b}$ which reduces to:

1. Permute: $P\mathbf{b}$
2. Forward substitution: solve $L\mathbf{y} = P\mathbf{b}$
3. Back substitution: solve $U\mathbf{x} = \mathbf{y}$

Steps 2 and 3 are each $\mathcal{O}(n^2)$. The expensive part is computing the
factorization itself ($\mathcal{O}(n^3)$), but this is done only once. Solving
for additional right-hand sides costs only $\mathcal{O}(n^2)$ each. See the
[appendix](../direct-methods/index.md) for details on LU decomposition.

**QR factorization** writes $A = QR$ where $Q$ is orthogonal and $R$ is upper
triangular. This is the [Gram-Schmidt](https://www.buttenschoen.ca/MATH545/orthogonality/qr.html)
process (Math 235) rewritten as a matrix
factorization. Then $A\mathbf{x} = \mathbf{b}$ becomes
$R\mathbf{x} = Q^T\mathbf{b}$: a matrix-vector multiply followed by back
substitution. QR is the focus of this chapter because it is more numerically
stable than LU and extends naturally to rectangular matrices and least squares
problems.

**SVD** ([singular value decomposition](https://www.buttenschoen.ca/MATH545/eigenvalues/svd.html))
writes $A = U\Sigma V^T$ where $U$ and
$V$ are orthogonal and $\Sigma$ is diagonal. This is connected to eigenvector
diagonalization and the fundamental theorem of linear algebra
([Math 545](https://www.buttenschoen.ca/MATH545)). The
SVD is the most powerful factorization for analysis (it reveals rank,
conditioning, and the geometry of the linear map), but it is also the most
expensive to compute.
