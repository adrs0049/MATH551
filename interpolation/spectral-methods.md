# Spectral Methods for Boundary Value Problems

:::{tip} Big Idea
Spectral methods solve differential equations by expanding solutions in global basis functions (Fourier or Chebyshev). For smooth problems, they achieve exponential accuracy—far surpassing finite difference methods. Spectral collocation enforces the differential equation at grid points, yielding a linear system that can be solved directly.
:::

## Overview

**Spectral methods** represent the solution to a differential equation as:
$$
u(x) \approx \sum_{k=0}^{N} a_k \phi_k(x)
$$

where $\{\phi_k\}$ are global basis functions:
- **Fourier:** $e^{ikx}$ for periodic problems
- **Chebyshev:** $T_k(x)$ for non-periodic problems

The key insight: differentiation becomes algebraic operations on coefficients!

## Why Spectral Methods?

For smooth problems, spectral methods offer:

| Method | Error decay | Points for 10-digit accuracy |
|--------|-------------|------------------------------|
| 2nd-order finite difference | $O(h^2)$ | ~$10^5$ |
| 4th-order finite difference | $O(h^4)$ | ~$10^{2.5} \approx 300$ |
| Spectral (smooth function) | $O(e^{-\alpha N})$ | ~20-50 |

**Exponential convergence** means spectral methods can achieve high accuracy with remarkably few points.

## Spectral Collocation

**Idea:** Approximate the solution as a polynomial that:
1. Satisfies the ODE at interior grid points
2. Satisfies boundary conditions at endpoints

Given BVP:
$$
\mathcal{L}u = f, \quad u(a) = \alpha, \quad u(b) = \beta
$$

At Chebyshev points $x_0, x_1, \ldots, x_N$ (with $x_0 = 1$ and $x_N = -1$ for $[-1, 1]$):

- **Interior equations:** $(\mathcal{L}u)_j = f_j$ for $j = 1, \ldots, N-1$
- **Boundary conditions:** $u_0 = \alpha$, $u_N = \beta$

## Second-Order BVPs

Consider:
$$
u''(x) = f(x), \quad x \in [-1, 1], \quad u(-1) = \alpha, \quad u(1) = \beta
$$

Let $D$ be the $(N+1) \times (N+1)$ Chebyshev differentiation matrix.

Then $D^2 = D \cdot D$ approximates the second derivative.

### Setting Up the System

The system $D^2 \mathbf{u} = \mathbf{f}$ with boundary conditions becomes:

$$
\begin{pmatrix}
1 & 0 & 0 & \cdots & 0 \\
D^2_{10} & D^2_{11} & D^2_{12} & \cdots & D^2_{1N} \\
\vdots & & & & \vdots \\
D^2_{N-1,0} & D^2_{N-1,1} & & \cdots & D^2_{N-1,N} \\
0 & 0 & 0 & \cdots & 1
\end{pmatrix}
\begin{pmatrix}
u_0 \\ u_1 \\ \vdots \\ u_{N-1} \\ u_N
\end{pmatrix}
=
\begin{pmatrix}
\beta \\ f_1 \\ \vdots \\ f_{N-1} \\ \alpha
\end{pmatrix}
$$

The first and last rows enforce boundary conditions.

### Implementation

```python
import numpy as np
import scipy.linalg as la

def chebpts(N):
    """Chebyshev points of the second kind."""
    return np.cos(np.pi * np.arange(N+1) / N)

def cheb_diff_matrix(N):
    """Chebyshev differentiation matrix."""
    x = chebpts(N)
    c = np.ones(N+1)
    c[0] = c[N] = 2
    c = c * ((-1.0) ** np.arange(N+1))

    X = np.tile(x, (N+1, 1))
    dX = X - X.T

    D = np.outer(c, 1.0/c) / (dX + np.eye(N+1))
    D = D - np.diag(D.sum(axis=1))
    return D

def solve_bvp(f, alpha, beta, N):
    """
    Solve u'' = f on [-1, 1] with u(-1) = alpha, u(1) = beta.
    """
    x = chebpts(N)
    D = cheb_diff_matrix(N)
    D2 = D @ D

    # Right-hand side
    rhs = f(x)

    # Apply boundary conditions
    D2[0, :] = 0
    D2[0, 0] = 1
    D2[N, :] = 0
    D2[N, N] = 1
    rhs[0] = beta   # u(1) = beta (x_0 = 1)
    rhs[N] = alpha  # u(-1) = alpha (x_N = -1)

    # Solve
    u = la.solve(D2, rhs)
    return x, u
```

## Example: Poisson Equation

Solve $u'' = e^{4x}$ on $[-1, 1]$ with $u(\pm 1) = 0$.

```python
f = lambda x: np.exp(4*x)
x, u = solve_bvp(f, 0, 0, N=16)

# Compare to exact solution
u_exact = (np.exp(4*x) - x*np.sinh(4) - np.cosh(4)) / 16
error = np.max(np.abs(u - u_exact))
print(f"Max error with N=16: {error:.2e}")  # Around 1e-13!
```

With just 16 points, we achieve near-machine-precision accuracy!

## Variable Coefficient Problems

For $a(x)u'' + b(x)u' + c(x)u = f(x)$:

$$
\text{diag}(a) \cdot D^2 + \text{diag}(b) \cdot D + \text{diag}(c)
$$

where $\text{diag}(a)$ is the diagonal matrix with $a(x_j)$ on the diagonal.

```python
def solve_variable_coeff_bvp(a, b, c, f, alpha, beta, N):
    """Solve a(x)u'' + b(x)u' + c(x)u = f(x)."""
    x = chebpts(N)
    D = cheb_diff_matrix(N)
    D2 = D @ D

    # Operator matrix
    A = np.diag(a(x)) @ D2 + np.diag(b(x)) @ D + np.diag(c(x))
    rhs = f(x)

    # Boundary conditions
    A[0, :] = 0; A[0, 0] = 1; rhs[0] = beta
    A[N, :] = 0; A[N, N] = 1; rhs[N] = alpha

    return x, la.solve(A, rhs)
```

## Neumann Boundary Conditions

For $u'(-1) = \alpha$ instead of $u(-1) = \alpha$:

Replace the last row of the system with the last row of $D$:

```python
# Instead of A[N, :] = [0, ..., 0, 1]
A[N, :] = D[N, :]  # Row for u'(x_N) = u'(-1)
rhs[N] = alpha     # u'(-1) = alpha
```

:::{prf:remark} Mixed Boundary Conditions
For problems with one Dirichlet and one Neumann condition, modify only the appropriate row. For two Neumann conditions, the problem may be ill-posed (the solution is only determined up to a constant).
:::

## Eigenvalue Problems

For $-u'' = \lambda u$ with $u(\pm 1) = 0$:

The interior equations give a generalized eigenvalue problem:
$$
-D^2_{int} \mathbf{u}_{int} = \lambda \mathbf{u}_{int}
$$

where the subscript "int" means we remove boundary rows/columns.

```python
def laplacian_eigenvalues(N):
    """Eigenvalues of -d^2/dx^2 on [-1,1] with Dirichlet BCs."""
    x = chebpts(N)
    D = cheb_diff_matrix(N)
    D2 = D @ D

    # Remove boundary points
    D2_int = D2[1:N, 1:N]

    # Eigenvalues
    eigs = la.eigvals(-D2_int)
    return np.sort(np.real(eigs))

# Compare to exact: lambda_k = (k*pi/2)^2
N = 20
computed = laplacian_eigenvalues(N)
exact = [(k * np.pi / 2)**2 for k in range(1, N)]
```

The exact eigenvalues are $\lambda_k = (k\pi/2)^2$ for $k = 1, 2, 3, \ldots$

## Convergence: Spectral vs. Finite Difference

For a smooth problem, compare spectral vs. finite difference:

```python
def solve_bvp_fd(f, alpha, beta, N):
    """Solve u'' = f using 2nd-order finite differences."""
    h = 2.0 / N
    x = np.linspace(-1, 1, N+1)

    # Tridiagonal matrix for u''
    A = np.diag(-2 * np.ones(N-1)) + np.diag(np.ones(N-2), 1) + np.diag(np.ones(N-2), -1)
    A /= h**2

    rhs = f(x[1:N])
    rhs[0] -= alpha / h**2
    rhs[-1] -= beta / h**2

    u_int = la.solve(A, rhs)
    return x, np.concatenate([[alpha], u_int, [beta]])
```

| $N$ | FD Error | Spectral Error |
|-----|----------|----------------|
| 8 | $10^{-2}$ | $10^{-4}$ |
| 16 | $10^{-3}$ | $10^{-9}$ |
| 32 | $10^{-4}$ | $10^{-14}$ |

Spectral methods improve exponentially; finite differences improve algebraically.

