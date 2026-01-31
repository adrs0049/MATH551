# Linear Algebra and QR Factorization

:::{tip} Big Idea
Linear systems $A\mathbf{x} = \mathbf{b}$ are fundamental to scientific computing. The QR factorization provides a numerically stable way to solve these systems—and handles least squares problems elegantly where other methods fail.
:::

## Linear Systems

Consider a system of $n$ equations in $n$ unknowns:

$$
\begin{aligned}
a_{11}x_1 + a_{12}x_2 + \cdots + a_{1n}x_n &= b_1 \\
a_{21}x_1 + a_{22}x_2 + \cdots + a_{2n}x_n &= b_2 \\
&\vdots \\
a_{n1}x_1 + a_{n2}x_2 + \cdots + a_{nn}x_n &= b_n
\end{aligned}
$$

In matrix form: $A\mathbf{x} = \mathbf{b}$ where $A \in \mathbb{R}^{n \times n}$, and $\mathbf{x}, \mathbf{b} \in \mathbb{R}^n$.

## The Fundamental Theorem of Linear Algebra

:::{prf:theorem} Invertibility Conditions
:label: thm-invertibility-conditions

For a matrix $A \in \mathbb{R}^{n \times n}$, the following are equivalent:

1. $A$ is **invertible** (i.e., $A^{-1}$ exists)
2. $A\mathbf{x} = \mathbf{b}$ has a **unique solution** for each $\mathbf{b}$
3. $A\mathbf{x} = \mathbf{0}$ has only the **trivial solution** $\mathbf{x} = \mathbf{0}$
4. $\det(A) \neq 0$
5. All **eigenvalues** of $A$ are non-zero
:::

## Matrix Factorizations: The Key to Solving Linear Systems

The naive approach to solving $A\mathbf{x} = \mathbf{b}$ is to compute $A^{-1}$ and then $\mathbf{x} = A^{-1}\mathbf{b}$. This is a **bad idea**:
- Computing $A^{-1}$ is expensive ($O(n^3)$ operations)
- It's numerically less stable than alternatives
- It doesn't exploit structure

Instead, we **factorize** the matrix: write $A$ as a product of simpler matrices that are easy to invert or solve with.

### Why Factorizations?

1. **Solve efficiently:** Triangular systems are cheap ($O(n^2)$)
2. **Reuse work:** Once factored, solve for multiple right-hand sides cheaply
3. **Reveal structure:** Factorizations expose rank, conditioning, eigenvalues

### The Big Three Factorizations

| Factorization | Form | Best For |
|---------------|------|----------|
| **LU** | $A = PLU$ | Square systems, multiple RHS |
| **QR** | $A = QR$ | Least squares, stability |
| **SVD** | $A = U\Sigma V^T$ | Rank-deficient problems, analysis |

**This chapter focuses on QR**, which handles rectangular matrices and least squares problems where LU doesn't apply.

:::{prf:remark} LU Factorization
:label: rmk-lu-factorization

LU factorization is essentially **Gaussian elimination** organized as a matrix product. Given $A$, we find:
- $P$ = permutation matrix (row swaps for pivoting)
- $L$ = lower triangular (multipliers from elimination)
- $U$ = upper triangular (the row-echelon form)

Then $A\mathbf{x} = \mathbf{b}$ becomes $L(U\mathbf{x}) = P^T\mathbf{b}$: solve $L\mathbf{y} = P^T\mathbf{b}$ (forward substitution), then $U\mathbf{x} = \mathbf{y}$ (back substitution).

**LU will be covered in homework** — the key ideas are in your linear algebra background.
:::

## Learning Outcomes

After completing this chapter, you should be able to:

- **L4.1:** State the fundamental theorem of linear algebra.
- **L4.2:** Define and compute vector and matrix norms.
- **L4.3:** Define condition number and explain its significance.
- **L4.4:** Define orthogonality and compute inner products.
- **L4.5:** State the best approximation theorem.
- **L4.6:** Describe Gram-Schmidt and its instability.
- **L4.7:** Explain Householder reflections.
- **L4.8:** Write the QR factorization (full and reduced).
- **L4.9:** Derive the normal equations.
- **L4.10:** Explain why $\kappa(A^T A) = \kappa(A)^2$.
- **L4.11:** Solve least squares using QR.
- **L4.12:** Compare normal equations vs QR stability.

