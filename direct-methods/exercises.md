# Exercises

Practice implementing and analyzing direct methods for linear systems.

---

### Q5.1: Triangular Solves

- **(a)** Solve by back substitution:
$$
\begin{pmatrix} 3 & 1 & 2 \\ 0 & 2 & 4 \\ 0 & 0 & 5 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} 9 \\ 8 \\ 10 \end{pmatrix}
$$

- **(b)** Solve by forward substitution:
$$
\begin{pmatrix} 2 & 0 & 0 \\ 1 & 3 & 0 \\ 4 & 2 & 1 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} 4 \\ 5 \\ 11 \end{pmatrix}
$$

---

### Q5.2: Gaussian Elimination

Solve the following system using Gaussian elimination (show all steps):
$$
\begin{aligned}
x_1 + 2x_2 + 3x_3 &= 6 \\
2x_1 + 3x_2 + x_3 &= 6 \\
3x_1 + x_2 + 2x_3 &= 6
\end{aligned}
$$

---

### Q5.3: LU Decomposition

- **(a)** Find the LU decomposition of:
$$
\mathcal{A} = \begin{pmatrix} 2 & 1 & 1 \\ 4 & 3 & 3 \\ 8 & 7 & 9 \end{pmatrix}
$$

- **(b)** Use your decomposition to solve $\mathcal{A}\mathbf{x} = (4, 10, 24)^T$.

- **(c)** Solve $\mathcal{A}\mathbf{x} = (3, 7, 17)^T$ (same LU, different RHS).

---

### Q5.4: Flop Counting

- **(a)** Verify that back substitution requires exactly $n^2$ flops.

- **(b)** Verify that Gaussian elimination requires $\frac{2}{3}n^3 + \mathcal{O}(n^2)$ flops.

- **(c)** For $n = 1000$, how many times faster is solving with a precomputed LU versus fresh Gaussian elimination?

---

### Q5.5: Pivoting

Consider:
$$
\begin{pmatrix} 0 & 1 \\ 1 & 1 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 1 \\ 2 \end{pmatrix}
$$

- **(a)** Why does Gaussian elimination without pivoting fail?

- **(b)** Apply Gaussian elimination with partial pivoting.

- **(c)** Write out the $\mathcal{P}\mathcal{A} = \mathcal{L}\mathcal{U}$ decomposition.

---

### Q5.6: Implementation

Implement the following in Python:

- **(a)** Back substitution for upper triangular systems.

- **(b)** LU decomposition without pivoting.

- **(c)** LU decomposition with partial pivoting.

Test on random $100 \times 100$ matrices and compare with `scipy.linalg.lu`.

---

### Q5.7: Stability Investigation

- **(a)** Create a $10 \times 10$ matrix where pivoting makes a significant difference. Solve the system with and without pivoting and compare errors.

- **(b)** Investigate the Hilbert matrix: compare the computed solution to the exact solution for increasing matrix sizes.

---

## Self-Assessment Questions

Test your understanding with these conceptual questions:

1. **Complexity:** How many flops does Gaussian elimination require for an $n \times n$ matrix?

2. **LU Advantage:** If you need to solve $\mathcal{A}\mathbf{x} = \mathbf{b}$ for 100 different $\mathbf{b}$ vectors, is it better to use Gaussian elimination 100 times or compute LU once?

3. **Pivoting:** Why does dividing by small numbers cause problems in floating point arithmetic?

4. **Back Substitution:** Why is back substitution $\mathcal{O}(n^2)$ rather than $\mathcal{O}(n^3)$?
