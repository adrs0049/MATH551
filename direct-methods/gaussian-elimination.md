# Gaussian Elimination

:::{tip} Big Idea
Gaussian elimination transforms a general matrix into upper triangular form
using row operations. Combined with back substitution, this solves any
invertible linear system.
:::

## Elementary Row Operations

Three operations that preserve the solution of a linear system:

1. **Multiply a row by a nonzero scalar:** $cR_i \to R_i$
2. **Add a multiple of one row to another:** $R_i - c R_j \to R_i$
3. **Interchange two rows:** $R_i \leftrightarrow R_j$

## The Algorithm

Transform the augmented matrix $[\mathcal{A} | \mathbf{b}]$ to upper triangular
form $[\mathcal{U} | \tilde{\mathbf{b}}]$ by eliminating entries below the
diagonal, column by column.

:::{prf:algorithm} Gaussian Elimination
:label: alg-gaussian-elimination

**Input:** $\mathcal{A} \in \mathbb{R}^{n \times n}$,
$\mathbf{b} \in \mathbb{R}^n$

**Output:** Upper triangular $\mathcal{U}$, modified $\tilde{\mathbf{b}}$

1. **for** $k = 1, 2, \ldots, n-1$:
2. $\qquad$ **for** $i = k+1, \ldots, n$:
3. $\qquad\qquad$ $l_{ik} \gets a_{ik} / a_{kk}$ *(multiplier)*
4. $\qquad\qquad$ **for** $j = k+1, \ldots, n$:
5. $\qquad\qquad\qquad$ $a_{ij} \gets a_{ij} - l_{ik} \cdot a_{kj}$
6. $\qquad\qquad$ $b_i \gets b_i - l_{ik} \cdot b_k$
:::

:::{prf:example} Gaussian Elimination by Hand
:label: ex-gaussian-elimination

Solve:

$$
\begin{aligned}
2x_1 + x_2 + 3x_3 &= 1 \\
4x_1 + 4x_2 + 7x_3 &= 1 \\
2x_1 + 5x_2 + 9x_3 &= 3
\end{aligned}
$$

**Step 1:** Form augmented matrix

$$
\left[\begin{array}{ccc|c} 2 & 1 & 3 & 1 \\ 4 & 4 & 7 & 1 \\ 2 & 5 & 9 & 3 \end{array}\right]
$$

**Step 2:** Eliminate below first pivot ($R_2 - 2R_1$, $R_3 - R_1$)

$$
\left[\begin{array}{ccc|c} 2 & 1 & 3 & 1 \\ 0 & 2 & 1 & -1 \\ 0 & 4 & 6 & 2 \end{array}\right]
$$

**Step 3:** Eliminate below second pivot ($R_3 - 2R_2$)

$$
\left[\begin{array}{ccc|c} 2 & 1 & 3 & 1 \\ 0 & 2 & 1 & -1 \\ 0 & 0 & 4 & 4 \end{array}\right]
$$

**Step 4:** Back substitution gives $\mathbf{x} = (-1/2, -1, 1)^T$.
:::

## Computational Cost

At step $k$, the elimination updates a $(n-k) \times (n-k)$ submatrix. The
total cost is:

$$
\sum_{k=1}^{n-1} 2(n-k)^2 = \frac{2}{3}n^3 + \mathcal{O}(n^2)
$$

Gaussian elimination is $\mathcal{O}(n^3)$: cubic in the matrix size. The
subsequent back substitution is only $\mathcal{O}(n^2)$, so the elimination
dominates.

:::{seealso}
[Gaussian Elimination and LU Decomposition](../notebooks/gaussian-elimination-notebook.ipynb): implementation of GE, elementary matrices, the connection to LU, and visualization of structured LU factors.
:::
