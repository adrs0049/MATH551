# Pivoting

:::{tip} Big Idea
Pivoting reorders rows (and sometimes columns) to avoid division by zero and
minimize round-off error. Partial pivoting, swapping rows to use the largest
pivot, is essential for numerical stability.
:::

## The Problem

Gaussian elimination fails if a pivot element $a_{kk}^{(k)} = 0$. Even when
not zero, small pivots cause problems: division by small numbers amplifies
errors and the algorithm becomes unstable.

## Partial Pivoting

At each step $k$, find the largest element in the current column (below the
diagonal) and swap rows before eliminating.

:::{prf:algorithm} Gaussian Elimination with Partial Pivoting
:label: alg-ge-partial-pivoting

**Input:** $\mathcal{A} \in \mathbb{R}^{n \times n}$,
$\mathbf{b} \in \mathbb{R}^n$

**Output:** Upper triangular $\mathcal{U}$, modified $\tilde{\mathbf{b}}$,
permutation record

1. **for** $k = 1, 2, \ldots, n-1$:
2. $\qquad$ Find $p = \arg\max_{i \geq k} |a_{ik}|$
3. $\qquad$ Swap rows $k$ and $p$ in $\mathcal{A}$ and $\mathbf{b}$
4. $\qquad$ **for** $i = k+1, \ldots, n$:
5. $\qquad\qquad$ $l_{ik} \gets a_{ik} / a_{kk}$
6. $\qquad\qquad$ **for** $j = k+1, \ldots, n$:
7. $\qquad\qquad\qquad$ $a_{ij} \gets a_{ij} - l_{ik} \cdot a_{kj}$
8. $\qquad\qquad$ $b_i \gets b_i - l_{ik} \cdot b_k$
:::

## PA = LU Decomposition

With pivoting, the factorization becomes:

$$
\mathcal{P}\mathcal{A} = \mathcal{L}\mathcal{U}
$$

where $\mathcal{P}$ is a **permutation matrix** that records the row swaps.

To solve $\mathcal{A}\mathbf{x} = \mathbf{b}$:

1. Apply permutation: $\mathcal{P}\mathbf{b}$
2. Forward solve: $\mathcal{L}\mathbf{y} = \mathcal{P}\mathbf{b}$
3. Back solve: $\mathcal{U}\mathbf{x} = \mathbf{y}$

:::{prf:example} Why Pivoting Matters
:label: ex-pivoting-matters

Consider in 3-digit arithmetic:

$$
\begin{pmatrix} 0.0001 & 1 \\ 1 & 1 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 1 \\ 2 \end{pmatrix}
$$

**Without pivoting:** the multiplier is $l_{21} = 1/0.0001 = 10000$. The
elimination produces $a_{22}^{(1)} = 1 - 10000 \cdot 1 = -9999$, but in
3-digit arithmetic this is rounded to $-10000$, and the subsequent back
substitution gives $x_2 \approx 1.000$ and $x_1 \approx 0.000$. The true
solution is $x_1 \approx 1.0001$, $x_2 \approx 0.9999$.

**With pivoting:** swap the rows first:

$$
\begin{pmatrix} 1 & 1 \\ 0.0001 & 1 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 2 \\ 1 \end{pmatrix}
$$

Now the multiplier is $l_{21} = 0.0001$, which is safe. The result is accurate
even in limited precision.
:::

## Types of Pivoting

| Strategy | Description | Cost |
|----------|-------------|------|
| **Partial pivoting** | Swap rows to get largest pivot in column | $\mathcal{O}(n^2)$ comparisons |
| **Complete pivoting** | Swap both rows and columns | $\mathcal{O}(n^3)$ comparisons |
| **Scaled pivoting** | Account for row scaling before choosing pivot | $\mathcal{O}(n^2)$ |

Partial pivoting is the standard choice: the overhead is $\mathcal{O}(n^2)$
comparisons, negligible compared to $\mathcal{O}(n^3)$ elimination, and it is
sufficient for nearly all practical problems.

:::{prf:remark} Pivoting in Practice
:label: rmk-pivoting-practice

All standard libraries (NumPy, SciPy, LAPACK) use partial pivoting by default.
When you call `scipy.linalg.solve(A, b)` or `scipy.linalg.lu(A)`, partial
pivoting is applied automatically.
:::

:::{seealso}
[Stability of Gaussian Elimination](../notebooks/stability-ge.ipynb): growth factors for random matrices and the structure of LU factors with partial pivoting.
:::
