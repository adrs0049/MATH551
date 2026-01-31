# Least Squares Problems

:::{tip} Big Idea
When a linear system has more equations than unknowns (overdetermined), we can't solve it exactly. Instead, we find the solution that minimizes the residual—this is the least squares solution.
:::

## The Problem

Given $A \in \mathbb{R}^{m \times n}$ with $m > n$ and $b \in \mathbb{R}^m$, the system $Ax = b$ typically has **no solution**.

**Why?** The vector $b$ usually doesn't lie in the column space of $A$.

**Goal:** Find $\hat{x}$ that minimizes the residual:

$$
\hat{x} = \arg\min_{x \in \mathbb{R}^n} \|Ax - b\|_2
$$

## Linear Regression Example

The most common application: fitting a model to data.

**Example:** Given $N$ data points $(t_i, y_i)$, fit a polynomial $p(t) = c_0 + c_1 t + c_2 t^2$:

$$
\begin{pmatrix} 1 & t_1 & t_1^2 \\ 1 & t_2 & t_2^2 \\ \vdots & \vdots & \vdots \\ 1 & t_N & t_N^2 \end{pmatrix} \begin{pmatrix} c_0 \\ c_1 \\ c_2 \end{pmatrix} = \begin{pmatrix} y_1 \\ y_2 \\ \vdots \\ y_N \end{pmatrix}
$$

With $N > 3$ data points, this system is overdetermined.

## Geometric Interpretation

The least squares solution finds the point in $\text{R}(A)$ closest to $b$:

```
        b
       /|
      / |  residual r = b - Ax̂
     /  |  is perpendicular to R(A)
    Ax̂--+---- R(A) (column space)
```

**Key insight:** The residual $r = b - A\hat{x}$ is orthogonal to the column space of $A$.

## The Normal Equations

From the orthogonality condition $r \perp \text{R}(A)$, we can derive:

:::{prf:theorem} Normal Equations
:label: thm-normal-equations

The least squares solution $\hat{x}$ satisfies:

$$
A^T A \hat{x} = A^T b
$$

If $A$ has full column rank, then $A^T A$ is invertible and:

$$
\hat{x} = (A^T A)^{-1} A^T b
$$
:::

**Derivation:** We need $(b - A\hat{x}) \perp Av$ for all $v \in \mathbb{R}^n$.

$$
(Av)^T(b - A\hat{x}) = 0 \quad \forall v
$$

$$
v^T A^T(b - A\hat{x}) = 0 \quad \forall v
$$

This requires $A^T(b - A\hat{x}) = 0$, giving $A^T A \hat{x} = A^T b$.

## Why Are They Called "Normal" Equations?

The name comes from the fact that the residual is **normal** (perpendicular) to the column space—not because they're "standard" equations.

## Properties of $A^T A$

When $A$ has full column rank:

| Property | Statement |
|----------|-----------|
| Symmetric | $(A^T A)^T = A^T A$ |
| Positive definite | $x^T A^T A x = \|Ax\|_2^2 > 0$ for $x \neq 0$ |
| Invertible | Follows from positive definiteness |
| Condition number | $\kappa(A^T A) = \kappa(A)^2$ ⚠️ |

:::{prf:proposition} Condition Number Squaring
:label: prop-condition-squaring

The condition number of $A^T A$ is the **square** of the condition number of $A$:

$$
\kappa_2(A^T A) = \kappa_2(A)^2
$$

This is why solving normal equations can be unstable.
:::

## Minimization Viewpoint

The least squares problem is equivalent to minimizing:

$$
f(x) = \|Ax - b\|_2^2 = x^T A^T A x - 2b^T A x + b^T b
$$

Taking the gradient and setting to zero:

$$
\nabla f(x) = 2A^T A x - 2A^T b = 0
$$

yields the normal equations.

**Observation:** This is a quadratic function with positive definite Hessian $2A^T A$, so there's a unique global minimum.

## Multiple Linear Regression

In statistics notation, the least squares problem for regression:

$$
\hat{\beta} = \arg\min_{\beta} \|y - X\beta\|_2^2
$$

where:
- $X \in \mathbb{R}^{N \times (k+1)}$ is the design matrix (observations of explanatory variables)
- $y \in \mathbb{R}^N$ is the response vector
- $\beta \in \mathbb{R}^{k+1}$ are the regression coefficients

The columns of $X$ typically include a column of ones (for the intercept).

## The Pseudoinverse

The matrix $(A^T A)^{-1} A^T$ is called the **Moore-Penrose pseudoinverse** $A^+$:

$$
A^+ = (A^T A)^{-1} A^T
$$

So $\hat{x} = A^+ b$.

**Properties:**
- $A^+ A = I_n$ (left inverse)
- $A A^+ \neq I_m$ in general (not a right inverse)
- $A A^+$ projects onto $\text{R}(A)$

## Numerical Considerations

:::{admonition} Don't Solve Normal Equations Directly!
:class: danger

Forming $A^T A$ and solving $(A^T A)\hat{x} = A^T b$:

1. **Squares the condition number:** $\kappa(A^T A) = \kappa(A)^2$
2. **Loses precision:** Errors amplified by squared condition number
3. **Requires forming product:** Extra $O(mn^2)$ operations

**Better approach:** Use QR factorization (next section)!
:::

## Example: Polynomial Fitting

```python
import numpy as np

# Data points
t = np.array([0, 1, 2, 3, 4])
y = np.array([1.0, 2.1, 3.9, 8.2, 15.8])

# Design matrix for quadratic fit
X = np.column_stack([np.ones_like(t), t, t**2])

# Normal equations (don't do this!)
# beta_bad = np.linalg.solve(X.T @ X, X.T @ y)

# Better: use lstsq which uses QR internally
beta, residuals, rank, s = np.linalg.lstsq(X, y, rcond=None)
```

## Residual Analysis

After finding $\hat{x}$, the residual is:

$$
r = b - A\hat{x}
$$

**Properties:**
- $\|r\|_2^2$ is the sum of squared errors (SSE)
- $r \perp \text{R}(A)$
- $A^T r = 0$

The residual measures how well the model fits the data.

## Summary

| Concept | Formula |
|---------|---------|
| Least squares problem | $\min_x \|Ax - b\|_2$ |
| Normal equations | $A^T A \hat{x} = A^T b$ |
| Solution (if full rank) | $\hat{x} = (A^T A)^{-1} A^T b$ |
| Residual | $r = b - A\hat{x} \perp \text{R}(A)$ |
| Condition number issue | $\kappa(A^T A) = \kappa(A)^2$ |

:::{admonition} Key Takeaway
:class: tip

The normal equations are mathematically correct but numerically dangerous. Use QR factorization instead—it avoids squaring the condition number.
:::
