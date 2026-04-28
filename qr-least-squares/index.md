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

**Linear algebra prerequisites.**

- **L4.1:** Define and compute vector and matrix norms.
- **L4.2:** Define orthogonality, compute inner products, and use orthogonal projections.

**Errors and conditioning.**

- **L4.3:** Distinguish forward error ($\|\hat{\mathbf{x}} - \mathbf{x}\|$) from backward error (residual or perturbed input), and explain why numerical analysis judges algorithms by the latter.
- **L4.4:** Define the matrix condition number $\kappa(A) = \|A\|\,\|A^{-1}\|$ and apply the **golden rule**: forward error $\lesssim \kappa(A) \cdot$ backward error.
- **L4.5:** Define **backward** and **forward stability**. State the central principle: a backward-stable algorithm applied to a well-conditioned problem produces a small forward error.
- **L4.6:** Distinguish *conditioning* (a property of the problem) from *stability* (a property of the algorithm). Recognise that a small residual alone does not imply a small error; the condition number is the bridge.
- **L4.7:** Produce an *a-posteriori* error bound on any computed solution to $A\mathbf{x} = \mathbf{b}$ from a condition-number estimate and the residual.

**Factorisations and how to use them.**

- **L4.8:** Reduce $A\mathbf{x} = \mathbf{b}$ to triangular form via factorisation and solve by forward/back substitution.
- **L4.9:** Apply the **LU factorisation** (with pivoting) to solve square linear systems $A\mathbf{x} = \mathbf{b}$.
- **L4.10:** Describe Gram-Schmidt and explain its loss of orthogonality at moderate $\kappa$.
- **L4.11:** Describe Householder reflections and explain why they yield a backward-stable QR.
- **L4.12:** Apply the **QR factorisation** to solve square linear systems $A\mathbf{x} = \mathbf{b}$.

**Optional: least squares.**

- **L4.13:** Use QR to solve overdetermined least-squares problems and compare with the normal equations; explain why $\kappa(A^\top A) = \kappa(A)^2$ matters in practice.
