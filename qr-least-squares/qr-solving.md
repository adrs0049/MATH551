# Solving Least Squares via QR

:::{tip} Big Idea
Using QR factorization to solve least squares problems avoids squaring the condition number, making it far more stable than the normal equations. The difference can be dramatic!
:::

## The Problem with Normal Equations

Recall the normal equations: $A^T A \hat{x} = A^T b$

**The issue:** $\kappa(A^T A) = \kappa(A)^2$

If $\kappa(A) = 10^8$, then $\kappa(A^T A) = 10^{16}$—near machine precision!

## Reduced vs Full QR

For least squares, we need to distinguish two forms of the QR factorization:

:::{prf:remark} Full vs Reduced QR
:label: rmk-full-reduced-qr
:class: dropdown

**Full QR:** $A = QR$ where $Q \in \mathbb{R}^{m \times m}$, $R \in \mathbb{R}^{m \times n}$

$$
\underbrace{A}_{m \times n} = \underbrace{Q}_{m \times m} \underbrace{\begin{pmatrix} R_1 \\ 0 \end{pmatrix}}_{m \times n}
$$

**Reduced (thin) QR:** $A = \hat{Q}\hat{R}$ where $\hat{Q} \in \mathbb{R}^{m \times n}$, $\hat{R} \in \mathbb{R}^{n \times n}$

$$
\underbrace{A}_{m \times n} = \underbrace{\hat{Q}}_{m \times n} \underbrace{\hat{R}}_{n \times n}
$$

**Which to Use?**
- **Reduced QR** is more memory-efficient (stores only $mn$ elements for $\hat{Q}$ vs $m^2$ for $Q$)
- **Reduced QR** is sufficient for solving the least squares problem
- **Full QR** is useful when you need to compute the residual norm or analyze the orthogonal complement
:::

## The QR Solution

Given the reduced QR factorization $A = \hat{Q}\hat{R}$:

$$
A\hat{x} = b \quad \Rightarrow \quad \hat{Q}\hat{R}\hat{x} = b
$$

The least squares solution satisfies $A\hat{x} = \text{proj}_{\text{R}(A)} b = \hat{Q}\hat{Q}^T b$.

Therefore:
$$
\hat{Q}\hat{R}\hat{x} = \hat{Q}\hat{Q}^T b
$$

Multiply both sides by $\hat{Q}^T$:
$$
\hat{R}\hat{x} = \hat{Q}^T b
$$

:::{prf:theorem} QR Solution to Least Squares
:label: thm-qr-least-squares

Given $A = \hat{Q}\hat{R}$ (reduced QR), the least squares solution is:

$$
\hat{R}\hat{x} = \hat{Q}^T b
$$

This is an upper triangular system—solve by back substitution!
:::

## Why This Works

**Step 1:** Project $b$ onto the column space of $A$:
$$
\text{proj}_{\text{R}(A)} b = \hat{Q}\hat{Q}^T b
$$

**Step 2:** Since $\text{R}(\hat{Q}) = \text{R}(A)$, there exists $\hat{x}$ with:
$$
A\hat{x} = \hat{Q}\hat{R}\hat{x} = \hat{Q}\hat{Q}^T b
$$

**Step 3:** Multiply by $\hat{Q}^T$ (using $\hat{Q}^T\hat{Q} = I$):
$$
\hat{R}\hat{x} = \hat{Q}^T b
$$

## Algorithm

```python
import numpy as np
from scipy.linalg import solve_triangular

def least_squares_qr(A, b):
    """Solve least squares via QR factorization."""
    Q, R = np.linalg.qr(A, mode='reduced')
    c = Q.T @ b
    x = solve_triangular(R, c, lower=False)
    return x
```

**Steps:**
1. Compute $\hat{Q}, \hat{R} = \text{qr}(A)$
2. Compute $c = \hat{Q}^T b$
3. Solve $\hat{R}\hat{x} = c$ by back substitution

## Condition Number Comparison

| Method | Effective Condition Number |
|--------|---------------------------|
| Normal equations | $\kappa(A)^2$ |
| QR factorization | $\kappa(A)$ |

The QR method never squares the condition number!

## Dramatic Example

This example from Trefethen and Bau shows the difference:

```python
import numpy as np

def fake_regression():
    """Ill-conditioned polynomial regression problem."""
    N = 100
    k = 14
    t = np.linspace(0, 1, N)
    X = np.column_stack([t**i for i in range(k+1)])
    y = np.exp(np.sin(4*t))
    y = y / 2006.787453080206  # normalize
    return X, y

X, y = fake_regression()

# Method 1: Normal equations (BAD!)
XtX = X.T @ X
Xty = X.T @ y
beta_normal = np.linalg.solve(XtX, Xty)

# Method 2: QR factorization (GOOD!)
Q, R = np.linalg.qr(X, mode='reduced')
beta_qr = np.linalg.solve(R, Q.T @ y)

print(f"Normal equations: beta[14] = {beta_normal[14]}")
print(f"QR factorization: beta[14] = {beta_qr[14]}")
print(f"True value:       beta[14] ≈ 1.0")
```

**Typical output:**
```
Normal equations: beta[14] = 0.3847...  (WRONG!)
QR factorization: beta[14] = 0.9999993... (correct to 7 figures)
```

## Why Such a Big Difference?

For this problem:
- $\kappa(X) \approx 10^9$
- $\kappa(X^T X) \approx 10^{18}$

With machine epsilon $\epsilon \approx 10^{-16}$:

| Method | Expected relative error |
|--------|------------------------|
| Normal equations | $\kappa(X^T X) \cdot \epsilon \approx 10^{2}$ (useless!) |
| QR factorization | $\kappa(X) \cdot \epsilon \approx 10^{-7}$ (7 digits) |

## The Residual

After solving, the residual is:

$$
r = b - A\hat{x} = b - \hat{Q}\hat{Q}^T b = (I - \hat{Q}\hat{Q}^T)b
$$

The residual norm squared (sum of squared errors):

$$
\|r\|_2^2 = \|b\|_2^2 - \|\hat{Q}^T b\|_2^2
$$

## Full QR Viewpoint

Using the full QR, $A = QR$ where $Q \in \mathbb{R}^{m \times m}$:

$$
Q^T b = \begin{pmatrix} c_1 \\ c_2 \end{pmatrix}
$$

where $c_1 \in \mathbb{R}^n$ and $c_2 \in \mathbb{R}^{m-n}$.

The least squares solution satisfies $R_1 \hat{x} = c_1$, and:

$$
\|r\|_2 = \|c_2\|_2
$$

The residual norm is determined by the "extra" part $c_2$.

## Comparison Summary

| Aspect | Normal Equations | QR Method |
|--------|-----------------|-----------|
| Condition | $\kappa(A)^2$ | $\kappa(A)$ |
| Form $A^T A$? | Yes | No |
| Matrix-matrix product | Yes ($O(mn^2)$) | No |
| Memory | $n \times n$ matrix | $m \times n$ factors |
| Stability | Poor | Excellent |
| Cost | $mn^2 + \frac{1}{3}n^3$ | $2mn^2 - \frac{2}{3}n^3$ |

## When to Use Each Method

**Use normal equations when:**
- $A$ is well-conditioned ($\kappa(A) < 10^7$ or so)
- $A^T A$ is already available
- Speed is critical and accuracy less so

**Use QR when:**
- $A$ may be ill-conditioned
- High accuracy is needed
- $A$ is sparse (specialized QR methods exist)
- You're not sure about the conditioning

**Rule of thumb:** When in doubt, use QR.

## Python Best Practices

```python
import numpy as np

# BEST: Use lstsq (automatically chooses good method)
x, residuals, rank, s = np.linalg.lstsq(A, b, rcond=None)

# GOOD: Explicit QR
Q, R = np.linalg.qr(A)
x = np.linalg.solve(R, Q.T @ b)

# AVOID: Normal equations
# x = np.linalg.solve(A.T @ A, A.T @ b)  # Don't do this!
```

## Summary

:::{admonition} Key Takeaways
:class: tip

1. **Normal equations** square the condition number—often disastrous
2. **QR method** preserves the original conditioning
3. **Solve $\hat{R}\hat{x} = \hat{Q}^T b$** by back substitution
4. **Use `np.linalg.lstsq`** in practice—it does the right thing
5. **The difference is dramatic** for ill-conditioned problems
:::

## Further Reading

- Trefethen & Bau, *Numerical Linear Algebra*, Lecture 11
- Golub & Van Loan, *Matrix Computations*, Chapter 5
