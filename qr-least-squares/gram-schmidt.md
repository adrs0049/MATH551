# Gram-Schmidt Orthogonalization

:::{tip} Big Idea
The Gram-Schmidt algorithm transforms any basis into an orthonormal basis—and in doing so, computes the QR factorization. Just as Gaussian elimination organized as matrix operations gives LU, Gram-Schmidt organized as matrix operations gives QR.
:::

:::{admonition} Prerequisite: Math 235
:class: note
You learned the Gram-Schmidt process in linear algebra (Math 235). Here we revisit it with two new perspectives: (1) understanding it as computing a matrix factorization, and (2) analyzing its numerical behavior.
:::

## The Problem

Given linearly independent vectors $a_1, \ldots, a_n$, find an orthonormal basis $q_1, \ldots, q_n$ for the same subspace.

## The Algorithm

The idea is simple: process vectors one at a time, subtracting off their projections onto previously computed orthonormal vectors.

:::{prf:algorithm} Classical Gram-Schmidt
:label: alg-classical-gs

**Input:** Linearly independent vectors $a_1, \ldots, a_n$

**Output:** Orthonormal vectors $q_1, \ldots, q_n$

1. **for** $j = 1, 2, \ldots, n$:
2. $\qquad v_j = a_j$
3. $\qquad$ **for** $i = 1, 2, \ldots, j-1$:
4. $\qquad\qquad v_j = v_j - \langle a_j, q_i \rangle \, q_i$ $\quad$ *(subtract projection)*
5. $\qquad q_j = v_j / \|v_j\|$ $\quad$ *(normalize)*
:::

:::{prf:example} Step-by-Step Gram-Schmidt
:label: ex-gs-step-by-step
:class: dropdown

Given $a_1 = \begin{pmatrix} 1 \\ 1 \\ 0 \end{pmatrix}$, $a_2 = \begin{pmatrix} 1 \\ 0 \\ 1 \end{pmatrix}$:

**Step 1:** Normalize $a_1$:
$$
q_1 = \frac{a_1}{\|a_1\|} = \frac{1}{\sqrt{2}} \begin{pmatrix} 1 \\ 1 \\ 0 \end{pmatrix}
$$

**Step 2:** Compute projection coefficient:
$$
\langle a_2, q_1 \rangle = \begin{pmatrix} 1 \\ 0 \\ 1 \end{pmatrix} \cdot \frac{1}{\sqrt{2}} \begin{pmatrix} 1 \\ 1 \\ 0 \end{pmatrix} = \frac{1}{\sqrt{2}}
$$

**Step 3:** Subtract projection:
$$
v_2 = a_2 - \langle a_2, q_1 \rangle \, q_1 = \begin{pmatrix} 1 \\ 0 \\ 1 \end{pmatrix} - \frac{1}{\sqrt{2}} \cdot \frac{1}{\sqrt{2}} \begin{pmatrix} 1 \\ 1 \\ 0 \end{pmatrix} = \begin{pmatrix} 1/2 \\ -1/2 \\ 1 \end{pmatrix}
$$

**Step 4:** Normalize:
$$
\|v_2\| = \sqrt{1/4 + 1/4 + 1} = \sqrt{3/2}
$$

$$
q_2 = \frac{v_2}{\|v_2\|} = \sqrt{\frac{2}{3}} \begin{pmatrix} 1/2 \\ -1/2 \\ 1 \end{pmatrix} = \begin{pmatrix} 1/\sqrt{6} \\ -1/\sqrt{6} \\ 2/\sqrt{6} \end{pmatrix}
$$

**Verify orthogonality:** $\langle q_1, q_2 \rangle = \frac{1}{\sqrt{2}} \cdot \frac{1}{\sqrt{6}} - \frac{1}{\sqrt{2}} \cdot \frac{1}{\sqrt{6}} + 0 = 0$ ✓
:::

## From Gram-Schmidt to QR Factorization

Gram-Schmidt is actually computing a QR factorization! Rearranging the algorithm:

$$
a_j = \sum_{i=1}^{j-1} \langle a_j, q_i \rangle \, q_i + \|v_j\| \, q_j
$$

This means $a_j$ is a linear combination of $q_1, \ldots, q_j$.

In matrix form:

$$
\underbrace{\begin{pmatrix} | & | & & | \\ a_1 & a_2 & \cdots & a_n \\ | & | & & | \end{pmatrix}}_A = \underbrace{\begin{pmatrix} | & | & & | \\ q_1 & q_2 & \cdots & q_n \\ | & | & & | \end{pmatrix}}_Q \underbrace{\begin{pmatrix} r_{11} & r_{12} & \cdots & r_{1n} \\ 0 & r_{22} & \cdots & r_{2n} \\ \vdots & & \ddots & \vdots \\ 0 & 0 & \cdots & r_{nn} \end{pmatrix}}_R
$$

where:
- $r_{ij} = \langle a_j, q_i \rangle$ for $i < j$ (the projection coefficients)
- $r_{jj} = \|v_j\|$ (the normalization factor)

:::{prf:remark} The Parallel: Gaussian Elimination → LU, Gram-Schmidt → QR
:label: rmk-gs-lu-parallel
:class: dropdown

| Algorithm | What it does | Matrix factorization |
|-----------|--------------|---------------------|
| Gaussian elimination | Row operations to get upper triangular | $A = LU$ |
| Gram-Schmidt | Orthogonalization to get orthonormal basis | $A = QR$ |

In both cases:
- The algorithm performs a sequence of operations on columns/rows
- Recording these operations as matrices gives the factorization
- $L$ stores the elimination multipliers; $R$ stores the projection coefficients
- $U$ is the reduced form; $Q$ is the orthonormal basis

**Making This Concrete:**

**Gaussian elimination:** Each step subtracts a multiple of one row from another. The multipliers $\ell_{ij}$ go into $L$; the result is $U$.

**Gram-Schmidt:** Each step subtracts projections onto previous vectors. The projection coefficients $r_{ij} = \langle a_j, q_i \rangle$ go into $R$; the orthonormal vectors form $Q$.

The triangular structure of $R$ reflects the sequential nature: $a_j$ only projects onto $q_1, \ldots, q_{j-1}$, so $r_{ij} = 0$ for $i > j$.
:::

## Numerical Problems with Classical Gram-Schmidt

Despite its elegance, classical Gram-Schmidt has serious numerical issues.

### Catastrophic Cancellation

:::{admonition} Warning
:class: warning

When the input vectors are nearly linearly dependent, Gram-Schmidt suffers from **catastrophic cancellation**.
:::

:::{prf:example} Nearly Parallel Columns
:label: ex-gs-nearly-parallel
:class: dropdown

Consider (from Trefethen & Bau):
$$
A = \begin{pmatrix} 0.70000 & 0.70711 \\ 0.70001 & 0.70711 \end{pmatrix}
$$

The columns are nearly parallel.

**Why This Fails:**

Computing $q_2$:

$$
q_2 = \frac{a_2 - \langle a_2, q_1 \rangle \, q_1}{\|a_2 - \langle a_2, q_1 \rangle \, q_1\|}
$$

The numerator involves subtracting two nearly equal vectors:

$$
a_2 - \langle a_2, q_1 \rangle \, q_1 \approx \begin{pmatrix} 0.70711 \\ 0.70711 \end{pmatrix} - 1.00000 \begin{pmatrix} 0.70710 \\ 0.70711 \end{pmatrix}
$$

**Cancellation!** We're subtracting numbers that agree in many digits, losing precision. The result is dominated by rounding errors.
:::

### Loss of Orthogonality

In exact arithmetic, Gram-Schmidt produces perfectly orthogonal vectors. In floating-point arithmetic, the computed $\hat{q}_i$ may **not** be orthogonal:

$$
\langle \hat{q}_i, \hat{q}_j \rangle \neq 0 \quad \text{for } i \neq j
$$

The loss of orthogonality scales with the condition number: for ill-conditioned matrices, the computed "orthonormal" vectors can be far from orthogonal. This defeats the entire purpose of orthogonalization!

## Modified Gram-Schmidt

A simple reordering of operations improves stability:

:::{prf:algorithm} Modified Gram-Schmidt
:label: alg-modified-gs

**Input:** Linearly independent vectors $a_1, \ldots, a_n$

**Output:** Orthonormal vectors $q_1, \ldots, q_n$ and upper triangular $R$

1. **for** $j = 1, 2, \ldots, n$:
2. $\qquad v_j = a_j$
3. $\qquad$ **for** $i = 1, 2, \ldots, j-1$:
4. $\qquad\qquad r_{ij} = \langle v_j, q_i \rangle$ $\quad$ *(use current $v_j$, not original $a_j$)*
5. $\qquad\qquad v_j = v_j - r_{ij} \, q_i$
6. $\qquad r_{jj} = \|v_j\|$
7. $\qquad q_j = v_j / r_{jj}$
:::

### Why Is Modified Gram-Schmidt Better?

The key difference is **when** we compute the projection coefficients:

| Classical | Modified |
|-----------|----------|
| $r_{ij} = \langle a_j, q_i \rangle$ | $r_{ij} = \langle v_j, q_i \rangle$ |
| Uses original $a_j$ | Uses current $v_j$ (already partially orthogonalized) |

:::{prf:remark} The Improvement in Modified GS
:label: rmk-mgs-improvement
:class: dropdown

In classical Gram-Schmidt, we compute all projections using the original vector $a_j$, then subtract them all at once. Rounding errors in the projections accumulate.

In modified Gram-Schmidt, we subtract each projection immediately, then compute the next projection using the updated vector. Each projection is computed against a vector that is **already more orthogonal** to the previous $q_i$'s.

**Analogy:** Think of it like this: if you're trying to make a vector orthogonal to $q_1$ and $q_2$, classical GS computes both corrections using the original vector and subtracts them together. Modified GS first makes the vector orthogonal to $q_1$, *then* computes how to make *that result* orthogonal to $q_2$.

The second approach is more accurate because the second projection starts from a better approximation.
:::

**However:** Modified Gram-Schmidt still has orthogonality loss proportional to $\kappa(A) \cdot \varepsilon_{\text{mach}}$. For ill-conditioned matrices, this can still be unacceptable. We need a different approach.
