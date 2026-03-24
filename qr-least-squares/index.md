# Linear Algebra and QR Factorization

:::{tip} Big Idea
Linear systems $A\mathbf{x} = \mathbf{b}$ are fundamental to scientific
computing. The QR factorization provides a numerically stable way to solve
these systems and handles least squares problems elegantly where other methods
fail.
:::

We begin by reviewing the linear algebra prerequisites (norms, inner products)
and introducing triangular systems, which are cheap to solve. This motivates
matrix factorizations: we decompose $A$ into simpler factors and solve the
resulting triangular systems. The chapter focuses on the **QR factorization**
$A = QR$, built from orthogonal matrices, which is numerically stable and
extends naturally to rectangular (overdetermined) systems and least squares
problems.

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
