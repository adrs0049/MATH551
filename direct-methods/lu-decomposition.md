# LU Decomposition

:::{tip} Big Idea
LU decomposition factors $\mathcal{A} = \mathcal{L}\mathcal{U}$ where
$\mathcal{L}$ is lower triangular and $\mathcal{U}$ is upper triangular. Once
computed, we can solve $\mathcal{A}\mathbf{x} = \mathbf{b}$ for any
$\mathbf{b}$ cheaply using two triangular solves.
:::

## The Factorization

Given an invertible matrix $\mathcal{A}$, we seek:

$$
\mathcal{A} = \mathcal{L}\mathcal{U}
$$

where $\mathcal{L}$ is lower triangular with ones on the diagonal and
$\mathcal{U}$ is upper triangular.

## Solving Linear Systems with LU

To solve $\mathcal{A}\mathbf{x} = \mathbf{b}$:

1. **Decompose:** compute $\mathcal{A} = \mathcal{L}\mathcal{U}$
   (cost: $\frac{2}{3}n^3$)
2. **Forward solve:** solve $\mathcal{L}\mathbf{y} = \mathbf{b}$ for
   $\mathbf{y}$ (cost: $n^2$)
3. **Back solve:** solve $\mathcal{U}\mathbf{x} = \mathbf{y}$ for
   $\mathbf{x}$ (cost: $n^2$)

The expensive decomposition is done only once. Solving for additional
right-hand sides costs only $\mathcal{O}(n^2)$ each.

## Connection to Gaussian Elimination

Gaussian elimination implicitly computes the LU decomposition:
- $\mathcal{U}$ is the final upper triangular matrix
- $\mathcal{L}$ stores the multipliers $l_{ik} = a_{ik}/a_{kk}$

$$
\mathcal{L} = \begin{pmatrix} 1 & & & \\ l_{21} & 1 & & \\ l_{31} & l_{32} & 1 & \\ \vdots & & \ddots & \ddots \\ l_{n1} & l_{n2} & \cdots & l_{n,n-1} & 1 \end{pmatrix}
$$

:::{prf:example} LU Decomposition
:label: ex-lu-decomposition

For $\mathcal{A} = \begin{pmatrix} 2 & 1 & 3 \\ 4 & 4 & 7 \\ 2 & 5 & 9 \end{pmatrix}$,
the Gaussian elimination from [](#ex-gaussian-elimination) gives multipliers
$l_{21} = 2$, $l_{31} = 1$, $l_{32} = 2$, so:

$$
\mathcal{L} = \begin{pmatrix} 1 & 0 & 0 \\ 2 & 1 & 0 \\ 1 & 2 & 1 \end{pmatrix}, \quad
\mathcal{U} = \begin{pmatrix} 2 & 1 & 3 \\ 0 & 2 & 1 \\ 0 & 0 & 4 \end{pmatrix}
$$

One can verify that $\mathcal{L}\mathcal{U} = \mathcal{A}$.
:::

## When Does LU Exist?

:::{prf:theorem} LU Existence and Uniqueness
:label: thm-lu-existence

If all leading principal minors of $\mathcal{A}$ are nonzero (i.e.,
$\det(\mathcal{A}_{1:k,1:k}) \neq 0$ for $k = 1, \ldots, n$), then the LU
decomposition exists and is unique.
:::

When this condition fails, we need [pivoting](pivoting.md) to proceed.
