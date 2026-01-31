# Exercises

Practice QR factorization and least squares.

---

### Q9.1: Orthogonality

Let $u = \begin{pmatrix} 1 \\ 2 \\ 2 \end{pmatrix}$ and $v = \begin{pmatrix} 2 \\ 1 \\ -2 \end{pmatrix}$.

- **(a)** Compute $u \cdot v$. Are $u$ and $v$ orthogonal?
- **(b)** Normalize $u$ to get a unit vector $\hat{u}$.
- **(c)** Compute the projection of $v$ onto $u$.

---

### Q9.2: Gram-Schmidt by Hand

Apply the Gram-Schmidt algorithm to the vectors:
$$
a_1 = \begin{pmatrix} 1 \\ 1 \\ 0 \end{pmatrix}, \quad a_2 = \begin{pmatrix} 1 \\ 0 \\ 1 \end{pmatrix}, \quad a_3 = \begin{pmatrix} 0 \\ 1 \\ 1 \end{pmatrix}
$$

- **(a)** Compute the orthonormal basis $q_1, q_2, q_3$.
- **(b)** Verify that $q_i \cdot q_j = \delta_{ij}$.
- **(c)** Write the QR factorization $A = QR$ where $A = [a_1 | a_2 | a_3]$.

---

### Q9.3: QR Factorization

For the matrix $A = \begin{pmatrix} 1 & 1 \\ 1 & 0 \\ 0 & 1 \end{pmatrix}$:

- **(a)** Compute the reduced QR factorization $A = \hat{Q}\hat{R}$ by hand.
- **(b)** Verify that $\hat{Q}^T \hat{Q} = I$.
- **(c)** Use Python to check your answer with `numpy.linalg.qr`.

---

### Q9.4: Normal Equations

Consider the least squares problem $\min_x \|Ax - b\|_2$ with:
$$
A = \begin{pmatrix} 1 & 1 \\ 1 & 2 \\ 1 & 3 \end{pmatrix}, \quad b = \begin{pmatrix} 1 \\ 2 \\ 2 \end{pmatrix}
$$

- **(a)** Form the normal equations $A^T A \hat{x} = A^T b$.
- **(b)** Solve for $\hat{x}$.
- **(c)** Compute the residual $r = b - A\hat{x}$ and verify that $A^T r = 0$.

---

### Q9.5: QR for Least Squares

Solve the same least squares problem from Q9.4 using QR factorization:

- **(a)** Compute the QR factorization of $A$.
- **(b)** Solve $\hat{R}\hat{x} = \hat{Q}^T b$ by back substitution.
- **(c)** Verify you get the same answer as Q9.4.

---

### Q9.6: Condition Number Comparison

For the polynomial fitting problem with design matrix:
$$
X_{ij} = t_i^{j-1}, \quad t_i = \frac{i-1}{n-1}, \quad i = 1, \ldots, n, \quad j = 1, \ldots, k+1
$$

- **(a)** Write Python code to construct $X$ for $n = 20$ and $k = 10$.
- **(b)** Compute $\kappa_2(X)$ and $\kappa_2(X^T X)$ using `numpy.linalg.cond`.
- **(c)** Verify that $\kappa_2(X^T X) \approx \kappa_2(X)^2$.

---

### Q9.7: Ill-Conditioned Regression

Implement the "dramatic example" from the notes:

```python
def fake_regression():
    N = 100
    k = 14
    t = np.linspace(0, 1, N)
    X = np.column_stack([t**i for i in range(k+1)])
    y = np.exp(np.sin(4*t))
    y = y / 2006.787453080206
    return X, y
```

- **(a)** Solve for regression coefficients using the normal equations.
- **(b)** Solve using `numpy.linalg.qr` and back substitution.
- **(c)** Compare $\hat{\beta}_{15}$ (Python index 14) to the true value of 1. Which method is more accurate?
- **(d)** Explain the difference in terms of condition numbers.

---

### Q9.8: Orthogonal Matrix Properties

Let $Q$ be an orthogonal matrix. Prove the following:

- **(a)** $\|Qx\|_2 = \|x\|_2$ for all $x$.
- **(b)** If $Q_1$ and $Q_2$ are orthogonal, then $Q_1 Q_2$ is orthogonal.
- **(c)** The eigenvalues of $Q$ have absolute value 1.

---

### Q9.9: Householder Reflection

- **(a)** For $v = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$, write the Householder matrix $H = I - 2vv^T$.
- **(b)** Show that $H$ reflects vectors across the $y$-axis.
- **(c)** Find the Householder reflector that maps $\begin{pmatrix} 3 \\ 4 \end{pmatrix}$ to $\begin{pmatrix} 5 \\ 0 \end{pmatrix}$.

---

### Q9.10: Polynomial Curve Fitting

Given data points $(0, 1), (1, 2), (2, 5), (3, 10), (4, 17)$:

- **(a)** Set up the least squares problem for fitting a quadratic $p(t) = c_0 + c_1 t + c_2 t^2$.
- **(b)** Solve using QR factorization (by hand or using Python).
- **(c)** Plot the data points and your fitted curve.
- **(d)** What is the sum of squared residuals?

---

### Q9.11: Operator Norm of Orthogonal Matrices

- **(a)** Show that $\|Q\|_2 = 1$ for any orthogonal matrix $Q$.
- **(b)** Show that $\|Q^{-1}\|_2 = 1$.
- **(c)** Conclude that $\kappa_2(Q) = 1$.

---

### Q9.12: Comparing Methods (Implementation)

Write Python functions to solve least squares three ways:

```python
def solve_normal_equations(A, b):
    """Solve via A^T A x = A^T b"""
    pass

def solve_qr(A, b):
    """Solve via QR factorization"""
    pass

def solve_lstsq(A, b):
    """Solve using numpy.linalg.lstsq"""
    pass
```

Test all three on the ill-conditioned polynomial problem and compare accuracy.

---

### Q9.13: Householder QR Implementation *(optional)*

Implement QR factorization using Householder reflections.

- **(a)** Write a function to compute the Householder vector $v$ such that $H = I - 2vv^T/\|v\|^2$ maps $x$ to $\|x\|_2 e_1$:
  ```python
  def householder_vector(x):
      """Return Householder vector v such that Hx = ||x||*e1."""
      pass
  ```

- **(b)** Write the full Householder QR factorization:
  ```python
  def householder_qr(A):
      """Compute QR factorization using Householder reflections.

      Returns:
          Q: Orthogonal matrix
          R: Upper triangular matrix
      """
      pass
  ```

- **(c)** Test your implementation against `numpy.linalg.qr` on random matrices.

- **(d)** Compare the numerical stability of your Householder QR to classical Gram-Schmidt on the ill-conditioned polynomial fitting problem from Q9.7.

:::{dropdown} Hint
For numerical stability, choose the sign of the Householder vector to avoid cancellation:
$$v = x + \text{sign}(x_1)\|x\|_2 e_1$$
If $x_1 = 0$, use the positive sign.
:::

---

### Q9.14: Givens Rotations *(optional)*

Givens rotations zero out individual matrix entries. The rotation matrix:

$$
G(i, j, \theta) = \begin{pmatrix}
1 & & & & \\
& \cos\theta & & -\sin\theta & \\
& & 1 & & \\
& \sin\theta & & \cos\theta & \\
& & & & 1
\end{pmatrix}
$$

rotates in the $(i,j)$ plane by angle $\theta$.

- **(a)** Given $a$ and $b$, derive formulas for $c = \cos\theta$ and $s = \sin\theta$ such that:
  $$
  \begin{pmatrix} c & -s \\ s & c \end{pmatrix} \begin{pmatrix} a \\ b \end{pmatrix} = \begin{pmatrix} r \\ 0 \end{pmatrix}
  $$

- **(b)** Implement a function to compute $c$ and $s$:
  ```python
  def givens_rotation(a, b):
      """Compute c, s such that [c -s; s c] @ [a; b] = [r; 0]."""
      pass
  ```

- **(c)** Implement QR factorization using Givens rotations:
  ```python
  def givens_qr(A):
      """Compute QR factorization using Givens rotations."""
      pass
  ```

- **(d)** When are Givens rotations preferred over Householder? (Hint: think about sparse matrices.)

:::{dropdown} Hint
For numerical stability, use:
$$
r = \sqrt{a^2 + b^2}, \quad c = a/r, \quad s = -b/r
$$
Handle the case $r = 0$ specially (set $c = 1$, $s = 0$).
:::

---

### Q9.15: Partial Pivoting *(optional)*

In Gaussian elimination, **partial pivoting** swaps rows to put the largest (in absolute value) element on the diagonal. This improves numerical stability.

- **(a)** Implement LU factorization with partial pivoting:
  ```python
  def lu_pivot(A):
      """LU factorization with partial pivoting.

      Returns:
          P: Permutation matrix (or permutation vector)
          L: Lower triangular with 1s on diagonal
          U: Upper triangular

      Such that PA = LU.
      """
      pass
  ```

- **(b)** Test on the matrix:
  $$
  A = \begin{pmatrix} 10^{-20} & 1 \\ 1 & 1 \end{pmatrix}
  $$
  Compare LU without pivoting vs. with pivoting. What happens to the computed solution of $Ax = b$ for $b = (1, 2)^T$?

- **(c)** Explain why the element $u_{22}$ becomes very large without pivoting, leading to large errors.

- **(d)** The **growth factor** is $\rho = \max_{ij} |u_{ij}| / \max_{ij} |a_{ij}|$. For the matrix above, compute $\rho$ with and without pivoting.

:::{dropdown} Hint
Without pivoting, the first step gives:
$$
L = \begin{pmatrix} 1 & 0 \\ 10^{20} & 1 \end{pmatrix}, \quad
U = \begin{pmatrix} 10^{-20} & 1 \\ 0 & 1 - 10^{20} \end{pmatrix}
$$

With pivoting, we swap rows first, giving much more reasonable values.
:::

---

## Self-Assessment Questions

Test your understanding with these conceptual questions:

1. **Orthogonality:** If $q_1, q_2$ are orthonormal, what is $q_1 \cdot q_2$? What is $\|q_1\|_2$?

2. **Best approximation:** Why is the residual $b - A\hat{x}$ orthogonal to the column space of $A$?

3. **Gram-Schmidt instability:** When does catastrophic cancellation occur in Gram-Schmidt?

4. **QR factorization:** If $A = QR$, what are the shapes of $Q$ and $R$ for $A \in \mathbb{R}^{100 \times 5}$?

5. **Condition number:** If $\kappa(A) = 10^6$, what is $\kappa(A^T A)$? Why is this bad?

6. **QR vs normal equations:** Why does QR give 7 correct digits when normal equations give none?

7. **Solving via QR:** Given $A = \hat{Q}\hat{R}$, what system do you solve for least squares?



### Q4.1: Vector Norms

Compute the 1-norm, 2-norm, and $\infty$-norm of the following vectors:

- **(a)** $\mathbf{x} = (3, -4)^T$
- **(b)** $\mathbf{x} = (1, 1, 1, 1)^T$
- **(c)** $\mathbf{x} = (1, 2, 3, 4, 5)^T$

---

### Q4.2: Matrix Norms

Compute $\|\mathcal{A}\|_1$ and $\|\mathcal{A}\|_\infty$ for:

- **(a)** $\mathcal{A} = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$

- **(b)** $\mathcal{A} = \begin{pmatrix} 2 & -1 & 0 \\ -1 & 2 & -1 \\ 0 & -1 & 2 \end{pmatrix}$

---

### Q4.3: Condition Number Computation

- **(a)** Compute the condition number $\kappa_1(\mathcal{A})$ for $\mathcal{A} = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}$.

- **(b)** Use NumPy to compute the condition number of the $10 \times 10$ Hilbert matrix.

- **(c)** For what size Hilbert matrix does the condition number exceed $1/\varepsilon_{\text{mach}}$ in double precision?

---

### Q4.4: Sensitivity Analysis

Consider the system $\mathcal{A}\mathbf{x} = \mathbf{b}$ where $\mathcal{A} = \begin{pmatrix} 1 & 1 \\ 1 & 1.0001 \end{pmatrix}$ and $\mathbf{b} = \begin{pmatrix} 2 \\ 2 \end{pmatrix}$.

- **(a)** Solve the system.
- **(b)** Compute the condition number.
- **(c)** Perturb $b_2$ by 0.0001 and solve again. How much does the solution change?
- **(d)** Is this consistent with the condition number?

---

### Q4.5: Iterative Refinement (Newton's Method in Disguise)

When solving $A\mathbf{x} = \mathbf{b}$, **iterative refinement** improves an approximate solution $\hat{\mathbf{x}}$ by repeatedly correcting it:

```
1. Compute residual: r = b - A*x̂
2. Solve for correction: A*δ = r
3. Update solution: x̂ ← x̂ + δ
4. Repeat until ||r|| is small enough
```

- **(a)** Let $F(\mathbf{x}) = A\mathbf{x} - \mathbf{b}$. Show that finding a root of $F(\mathbf{x}) = \mathbf{0}$ is equivalent to solving $A\mathbf{x} = \mathbf{b}$.

- **(b)** Write Newton's method for solving $F(\mathbf{x}) = \mathbf{0}$. Recall that Newton's method is:
  $$\mathbf{x}_{k+1} = \mathbf{x}_k - [F'(\mathbf{x}_k)]^{-1} F(\mathbf{x}_k)$$
  What is the Jacobian $F'(\mathbf{x})$ for a linear function?

- **(c)** Show that Newton's method applied to $F(\mathbf{x}) = A\mathbf{x} - \mathbf{b}$ gives exactly the iterative refinement algorithm above.

- **(d)** **Implementation:** Write Python code for iterative refinement. Use the LU factorization of $A$ (computed once) to solve for each correction $\delta$.

  ```python
  def iterative_refinement(A, b, x0, tol=1e-12, max_iter=10):
      """Improve solution x0 using iterative refinement."""
      # Your code here
      pass
  ```

- **(e)** Test on an ill-conditioned system. Create a $10 \times 10$ Hilbert matrix, solve $A\mathbf{x} = \mathbf{b}$ where $\mathbf{b} = A \cdot \mathbf{1}$ (so the true solution is all ones), and apply iterative refinement. How many iterations are needed?

:::{dropdown} Hint
For part (d), compute the LU factorization once with `scipy.linalg.lu_factor`, then use `scipy.linalg.lu_solve` for each correction. The key insight is that the expensive $O(n^3)$ factorization is reused, making each refinement step only $O(n^2)$.

For part (e), computing the residual $\mathbf{r} = \mathbf{b} - A\hat{\mathbf{x}}$ in higher precision (or carefully) is important for ill-conditioned problems. In practice, LAPACK's `xGERFS` routine does this.
:::

---

### Q4.6: Componentwise Condition Number (Skeel/Bauer)

The standard condition number $\kappa(A) = \|A\| \|A^{-1}\|$ can be overly pessimistic for **badly scaled** matrices. Consider:

$$
A = \begin{pmatrix} 10^{-6} & 0 \\ 0 & 1 \end{pmatrix}, \quad \mathbf{b} = \begin{pmatrix} 10^{-6} \\ 1 \end{pmatrix}
$$

- **(a)** Compute $\kappa_\infty(A) = \|A\|_\infty \|A^{-1}\|_\infty$. Is this matrix "ill-conditioned" by this measure?

- **(b)** Solve $A\mathbf{x} = \mathbf{b}$. What is the exact solution? Now perturb $A$ by relative errors of size $\varepsilon$ in each entry (i.e., replace $a_{ij}$ with $a_{ij}(1 + \varepsilon_{ij})$ where $|\varepsilon_{ij}| \leq \varepsilon$). How much does $\mathbf{x}$ change?

- **(c)** The **componentwise relative condition number** is defined as:
  $$\kappa_{\text{CR}}(A) = \| |A^{-1}| \cdot |A| \|$$
  where $|A|$ denotes the matrix of absolute values. Compute $\kappa_{\text{CR}}(A)$ for the diagonal matrix above.

- **(d)** **General principle:** Show that for *any* diagonal matrix $D$, we have $\kappa_{\text{CR}}(D) = 1$. This explains why diagonal systems are always "easy" even when $\kappa(D)$ is large.

- **(e)** More generally, if $D$ is diagonal and $B$ is any invertible matrix, show that:
  $$\kappa_{\text{CR}}(DB) = \kappa_{\text{CR}}(B)$$
  This means row scaling doesn't affect the componentwise condition number.

:::{dropdown} Hint
For part (c), note that for a diagonal matrix $A = \text{diag}(d_1, \ldots, d_n)$:
- $A^{-1} = \text{diag}(1/d_1, \ldots, 1/d_n)$
- $|A^{-1}| \cdot |A| = \text{diag}(|1/d_1| \cdot |d_1|, \ldots) = I$

For part (e), use $(DB)^{-1} = B^{-1}D^{-1}$ and the fact that $|D^{-1}| \cdot |D| = I$ for diagonal $D$.
:::

:::{admonition} Why This Matters
:class: note

Demmel's *Applied Numerical Linear Algebra* (Section 2.2.1) shows that when perturbations satisfy $|\delta A| \leq \varepsilon |A|$ (relative perturbations in each entry), the error bound becomes:

$$\frac{\|\delta \mathbf{x}\|}{\|\mathbf{x}\|} \lesssim \kappa_{\text{CR}}(A) \cdot \varepsilon$$

This is often much smaller than the bound using $\kappa(A)$, especially for badly scaled systems. LAPACK routines like `xGESVX` compute this tighter bound.
:::

---

### Q4.7: Row Equilibration

**Row equilibration** scales each row of $A$ so that all rows have comparable norms. Specifically, if we define the diagonal matrix:

$$
D_r = \text{diag}\left(\frac{1}{\|a_1\|_\infty}, \frac{1}{\|a_2\|_\infty}, \ldots, \frac{1}{\|a_n\|_\infty}\right)
$$

where $a_i^T$ is the $i$-th row of $A$, then $\tilde{A} = D_r A$ has $\|\tilde{a}_i\|_\infty = 1$ for all rows.

- **(a)** Consider the badly scaled matrix:
  $$
  A = \begin{pmatrix} 10^{-8} & 10^{-8} \\ 1 & 2 \end{pmatrix}
  $$
  Compute $\kappa_\infty(A)$. Then compute the row equilibration matrix $D_r$, form $\tilde{A} = D_r A$, and compute $\kappa_\infty(\tilde{A})$.

- **(b)** **Proof:** Show that for any invertible diagonal matrix $D$, the scaled system $DA\mathbf{x} = D\mathbf{b}$ has the same solution as $A\mathbf{x} = \mathbf{b}$. Why does this mean row equilibration doesn't change the mathematical problem?

- **(c)** **Proof:** Show that row equilibration does *not* change the componentwise condition number:
  $$
  \kappa_{\text{CR}}(D_r A) = \kappa_{\text{CR}}(A)
  $$
  (Hint: Use the result from Q4.6(e).)

- **(d)** **Implementation:** Write a Python function for row equilibration:
  ```python
  def row_equilibrate(A):
      """
      Compute row equilibration of A.

      Returns:
          A_eq: Row-equilibrated matrix
          D_r: Diagonal scaling matrix (as 1D array)
      """
      # Your code here
      pass
  ```

- **(e)** Test your implementation on the matrix from part (a). Solve $A\mathbf{x} = \mathbf{b}$ where $\mathbf{b} = A \cdot (1, 1)^T$ both with and without equilibration. Which gives a more accurate solution?

:::{dropdown} Hint
For part (c), recall that $\kappa_{\text{CR}}(A) = \||A^{-1}| \cdot |A|\|$. For a diagonal matrix $D$ with positive diagonal entries, we have $|D| = D$ and $|D^{-1}| = D^{-1}$, so $|D^{-1}||D| = I$.
:::

---

### Q4.8: Row and Column Equilibration

For highly ill-conditioned matrices, we can apply scaling to both rows *and* columns:

$$
\tilde{A} = D_r A D_c
$$

where $D_r$ scales rows and $D_c$ scales columns. The transformed system becomes:

$$
\tilde{A}\mathbf{y} = \tilde{\mathbf{b}}, \quad \text{where } \tilde{\mathbf{b}} = D_r \mathbf{b} \text{ and } \mathbf{x} = D_c \mathbf{y}
$$

- **(a)** **Proof:** Verify that if $\tilde{A}\mathbf{y} = \tilde{\mathbf{b}}$ with $\tilde{A} = D_r A D_c$ and $\tilde{\mathbf{b}} = D_r \mathbf{b}$, then $\mathbf{x} = D_c \mathbf{y}$ satisfies $A\mathbf{x} = \mathbf{b}$.

- **(b)** Consider the matrix:
  $$
  A = \begin{pmatrix} 1 & 10^6 \\ 10^{-6} & 1 \end{pmatrix}
  $$
  Compute $\kappa_\infty(A)$. Now apply row equilibration only (as in Q4.7) and compute the new condition number. Finally, apply column equilibration to the result and compute the final condition number. How much did equilibration help?

- **(c)** **Proof:** Unlike row scaling, column scaling *can* change the componentwise condition number. Show by example that $\kappa_{\text{CR}}(AD_c) \neq \kappa_{\text{CR}}(A)$ in general.

  (Hint: Try $A = \begin{pmatrix} 1 & 1 \\ 1 & 2 \end{pmatrix}$ and $D_c = \text{diag}(1, 100)$.)

- **(d)** **Implementation:** Write a complete equilibration solver:
  ```python
  def solve_equilibrated(A, b):
      """
      Solve Ax = b using row and column equilibration.

      Steps:
      1. Compute row scaling D_r (as in Q4.7)
      2. Apply row scaling: A1 = D_r @ A, b1 = D_r @ b
      3. Compute column scaling D_c for A1
      4. Apply column scaling: A2 = A1 @ D_c
      5. Solve A2 @ y = b1
      6. Recover x = D_c @ y

      Returns:
          x: Solution to original system
          kappa_original: Condition number of A
          kappa_equilibrated: Condition number of equilibrated A
      """
      # Your code here
      pass
  ```

- **(e)** **Experiment:** Generate a random ill-conditioned matrix using:
  ```python
  def make_ill_conditioned(n, target_kappa):
      """Create an n×n matrix with condition number ≈ target_kappa."""
      U, _ = np.linalg.qr(np.random.randn(n, n))
      V, _ = np.linalg.qr(np.random.randn(n, n))
      s = np.logspace(0, -np.log10(target_kappa), n)
      return U @ np.diag(s) @ V.T
  ```

  For $n = 50$ and `target_kappa = 1e12`:
  - Solve $A\mathbf{x} = \mathbf{b}$ directly using `numpy.linalg.solve`
  - Solve using your equilibration routine
  - Compare the errors when the true solution is $\mathbf{x} = (1, 1, \ldots, 1)^T$

:::{dropdown} Hint
For column equilibration, you can use the $\infty$-norm of each column:
```python
D_c = np.diag(1.0 / np.max(np.abs(A), axis=0))
```

For part (c), you need to compute $|A^{-1}||A|$ before and after column scaling. The column scaling changes the product in a non-trivial way.
:::

:::{admonition} Why This Matters
:class: note

LAPACK's expert driver routines (`xGESVX`, `xPOSVX`) automatically perform equilibration before solving. The routine returns:
- The equilibrated condition number
- Error bounds based on the equilibrated problem
- A flag indicating if the matrix was badly scaled

The goal of equilibration is to transform the problem so that round-off errors in the solve are as small as possible, while not changing the mathematical answer.
:::

---

## Self-Assessment Questions

Test your understanding with these conceptual questions:

1. **Invertibility:** If $\det(\mathcal{A}) = 0$, what can you say about the solution of $\mathcal{A}\mathbf{x} = \mathbf{b}$?

2. **Norms:** For $\mathbf{x} = (1, -2, 3)^T$, compute $\|\mathbf{x}\|_1$, $\|\mathbf{x}\|_2$, and $\|\mathbf{x}\|_\infty$.

3. **Condition Number:** If $\kappa(\mathcal{A}) = 10^8$, approximately how many digits of accuracy can you expect in double precision?

4. **Equilibration:** Why does row scaling not change the componentwise condition number, but column scaling can?

5. **Iterative Refinement:** Why do we compute the residual in higher precision (or carefully) for ill-conditioned problems?
