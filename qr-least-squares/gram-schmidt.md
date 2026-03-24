---
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/gram-schmidt.pdf
    id: qr-least-squares-gram-schmidt-pdf
downloads:
  - id: qr-least-squares-gram-schmidt-pdf
    title: Download PDF
---

# Gram-Schmidt Orthogonalization

:::{tip} Big Idea
The Gram-Schmidt algorithm transforms any basis into an orthonormal basis—and in
doing so, computes the QR factorization. Just as Gaussian elimination organized
as matrix operations gives LU, Gram-Schmidt organized as matrix operations gives
QR.
:::

:::{admonition} Prerequisite: Math 235
:class: note
You learned the Gram-Schmidt process in linear algebra (Math 235). Here we
revisit it with two new perspectives: (1) understanding it as computing a matrix
factorization, and (2) analyzing its numerical behavior.
:::

## The Problem

Given linearly independent vectors $a_1, \ldots, a_n$, find an orthonormal basis
$q_1, \ldots, q_n$ for the same subspace.

## The Algorithm

The idea is simple: process vectors one at a time, subtracting off their
projections onto previously computed orthonormal vectors.

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

**Verify orthogonality:** $\langle q_1, q_2 \rangle = \frac{1}{\sqrt{2}} \cdot
\frac{1}{\sqrt{6}} - \frac{1}{\sqrt{2}} \cdot \frac{1}{\sqrt{6}} + 0 = 0$
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

## Numerical Problems with Classical Gram-Schmidt

Despite its elegance, classical Gram-Schmidt has serious numerical issues.

When the columns of $A$ are nearly linearly dependent, the subtraction
$v_j = a_j - \sum \langle a_j, q_i \rangle q_i$ involves subtracting nearly
equal vectors. This **catastrophic cancellation** means the computed $v_j$ is
dominated by rounding errors, and the resulting $\hat{q}_j$ is not orthogonal
to the previous vectors.

The loss of orthogonality scales with the condition number of $A$: the computed
vectors satisfy $|\langle \hat{q}_i, \hat{q}_j \rangle| = \mathcal{O}(\kappa(A) \cdot \varepsilon_{\text{mach}})$.
For ill-conditioned matrices, this can be unacceptably large.

:::{prf:example} Loss of Orthogonality
:label: ex-gs-nearly-parallel
:class: dropdown

Consider (from Trefethen and Bau):

$$
A = \begin{pmatrix} 0.70000 & 0.70711 \\ 0.70001 & 0.70711 \end{pmatrix}
$$

The columns are nearly parallel. Computing $q_2$ requires:

$$
v_2 = a_2 - \langle a_2, q_1 \rangle \, q_1 \approx \begin{pmatrix} 0.70711 \\ 0.70711 \end{pmatrix} - 1.00000 \begin{pmatrix} 0.70710 \\ 0.70711 \end{pmatrix}
$$

The subtraction cancels most significant digits. The result is dominated by
rounding errors, so $\hat{q}_2$ is far from orthogonal to $q_1$.
:::

This motivates the [Householder approach](qr-factorization.md): instead of
orthogonalizing columns by subtraction (which cancels), we triangularize $A$ by
multiplication with orthogonal matrices, which preserves norms and avoids
cancellation entirely.
