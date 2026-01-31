# Iterative Methods for Linear Systems

:::{tip} Big Idea
Many numerical problems reduce to finding $\mathbf{x}$ satisfying $\mathbf{x} = M\mathbf{x} + \mathbf{c}$. This is a fixed point problem. Three perspectives—**Banach contraction**, **Neumann series**, and **matrix splitting**—give you the same answer and reveal deep connections across numerical analysis.
:::

## The Unifying Framework

Consider solving the linear system $A\mathbf{x} = \mathbf{b}$. Instead of direct factorization, we can reformulate as a fixed point problem:

$$
\mathbf{x} = M\mathbf{x} + \mathbf{c}
$$

and iterate:

$$
\mathbf{x}^{(k+1)} = M\mathbf{x}^{(k)} + \mathbf{c}
$$

**When does this converge?** Three equivalent perspectives:

## Perspective 1: Banach Fixed Point Theorem

:::{prf:theorem} Convergence of Fixed Point Iteration
:label: thm-matrix-fixed-point

If $\rho(M) < 1$ (spectral radius less than 1), the iteration $\mathbf{x}^{(k+1)} = M\mathbf{x}^{(k)} + \mathbf{c}$ converges to a unique fixed point $\mathbf{x}^*$ from any starting point.
:::

:::{prf:remark} Connection to Scalar Case
:label: rmk-matrix-scalar-connection
:class: dropdown

This is the **matrix version** of the Banach Fixed Point Theorem from nonlinear equations!

| Scalar | Matrix |
|--------|--------|
| $x_{n+1} = g(x_n)$ | $\mathbf{x}^{(k+1)} = M\mathbf{x}^{(k)} + \mathbf{c}$ |
| Converges if $|g'(x)| < 1$ | Converges if $\rho(M) < 1$ |
| Rate: $|e_n| \leq q^n |e_0|$ | Rate: $\|\mathbf{e}^{(k)}\| \approx C \rho(M)^k$ |

The unifying principle: **contractions converge geometrically**.
:::

## Perspective 2: Neumann Series

Rearrange $\mathbf{x} = M\mathbf{x} + \mathbf{c}$ to get $(I - M)\mathbf{x} = \mathbf{c}$, so:

$$
\mathbf{x} = (I - M)^{-1}\mathbf{c}
$$

When does this inverse exist? The **Neumann series** gives the answer:

:::{prf:theorem} Neumann Series
:label: thm-neumann-series

If $\rho(M) < 1$, then $(I - M)$ is invertible and:

$$
(I - M)^{-1} = \sum_{k=0}^{\infty} M^k = I + M + M^2 + M^3 + \cdots
$$
:::

:::{prf:proof}
:class: dropdown

Let $S_n = I + M + M^2 + \cdots + M^n$. Then:
$$
(I - M)S_n = I - M^{n+1}
$$
If $\rho(M) < 1$, then $M^{n+1} \to 0$ as $n \to \infty$, so $S_n \to (I-M)^{-1}$.

This is the matrix analog of the geometric series:
$$
\frac{1}{1-r} = 1 + r + r^2 + \cdots \quad \text{for } |r| < 1
$$
:::

## Perspective 3: Iterates Are Neumann Partial Sums

Here's the beautiful connection: **the fixed point iterates compute the Neumann series term by term**.

:::{prf:theorem} Iterates = Partial Sums
:label: thm-iterates-partial-sums

The fixed point iteration starting from $\mathbf{x}^{(0)}$ gives:

$$
\mathbf{x}^{(k)} = M^k \mathbf{x}^{(0)} + \left(\sum_{j=0}^{k-1} M^j\right)\mathbf{c}
$$
:::

:::{prf:proof}
:class: dropdown

Expand the recurrence $\mathbf{x}^{(k+1)} = M\mathbf{x}^{(k)} + \mathbf{c}$:

$$
\mathbf{x}^{(1)} = M\mathbf{x}^{(0)} + \mathbf{c}
$$

$$
\mathbf{x}^{(2)} = M\mathbf{x}^{(1)} + \mathbf{c} = M^2\mathbf{x}^{(0)} + M\mathbf{c} + \mathbf{c}
$$

$$
\mathbf{x}^{(3)} = M\mathbf{x}^{(2)} + \mathbf{c} = M^3\mathbf{x}^{(0)} + M^2\mathbf{c} + M\mathbf{c} + \mathbf{c}
$$

By induction: $\mathbf{x}^{(k)} = M^k \mathbf{x}^{(0)} + \sum_{j=0}^{k-1} M^j \mathbf{c}$.
:::

As $k \to \infty$ with $\rho(M) < 1$:
- $M^k \mathbf{x}^{(0)} \to \mathbf{0}$ — the initial condition **washes out**
- $\sum_{j=0}^{k-1} M^j \to (I - M)^{-1}$

**Special case:** If $\mathbf{x}^{(0)} = \mathbf{0}$, then $\mathbf{x}^{(k)} = \left(\sum_{j=0}^{k-1} M^j\right)\mathbf{c}$ — the iterates *are* the Neumann partial sums applied to $\mathbf{c}$.

## Matrix Splittings for $A\mathbf{x} = \mathbf{b}$

To solve $A\mathbf{x} = \mathbf{b}$, split $A = N - P$ where $N$ is easy to invert:

$$
N\mathbf{x} = P\mathbf{x} + \mathbf{b} \implies \mathbf{x} = N^{-1}P\mathbf{x} + N^{-1}\mathbf{b} = M\mathbf{x} + \mathbf{c}
$$

where $M = N^{-1}P$ is the **iteration matrix**.

**Convergence** ⟺ $\rho(M) < 1$ ⟺ Neumann series for $(I - M)^{-1}$ converges.

### Jacobi Method

Split $A = D - (L + U)$ where $D$ = diagonal, $L$ = strict lower, $U$ = strict upper.

- $N = D$, $P = L + U$
- $M_J = D^{-1}(L + U)$

:::{prf:algorithm} Jacobi Iteration
:label: alg-jacobi

**Input:** $A$, $\mathbf{b}$, initial guess $\mathbf{x}^{(0)}$, tolerance $\varepsilon$

**Output:** Approximate solution $\mathbf{x}$

1. **for** $k = 0, 1, 2, \ldots$:
2. $\qquad$ **for** $i = 1, \ldots, n$:
3. $\qquad\qquad$ $x_i^{(k+1)} = \frac{1}{a_{ii}}\left(b_i - \sum_{j \neq i} a_{ij}x_j^{(k)}\right)$
4. $\qquad$ **if** $\|\mathbf{x}^{(k+1)} - \mathbf{x}^{(k)}\| < \varepsilon$: **return** $\mathbf{x}^{(k+1)}$
:::

**Key property:** All components use values from the **previous** iteration — **parallelizable**.

### Gauss-Seidel Method

Use updated values immediately as they become available.

- $N = D - L$, $P = U$
- $M_{GS} = (D - L)^{-1}U$

:::{prf:algorithm} Gauss-Seidel Iteration
:label: alg-gauss-seidel

**Input:** $A$, $\mathbf{b}$, initial guess $\mathbf{x}^{(0)}$, tolerance $\varepsilon$

**Output:** Approximate solution $\mathbf{x}$

1. **for** $k = 0, 1, 2, \ldots$:
2. $\qquad$ **for** $i = 1, \ldots, n$:
3. $\qquad\qquad$ $x_i^{(k+1)} = \frac{1}{a_{ii}}\left(b_i - \sum_{j < i} a_{ij}x_j^{(k+1)} - \sum_{j > i} a_{ij}x_j^{(k)}\right)$
4. $\qquad$ **if** $\|\mathbf{x}^{(k+1)} - \mathbf{x}^{(k)}\| < \varepsilon$: **return** $\mathbf{x}^{(k+1)}$
:::

**Key property:** Uses **new** values for $j < i$ — typically faster convergence but **not parallelizable**.

### Comparison

| Method | Iteration Matrix $M$ | Parallelizable | Typical Convergence |
|--------|---------------------|----------------|---------------------|
| Jacobi | $D^{-1}(L + U)$ | Yes | Slower |
| Gauss-Seidel | $(D-L)^{-1}U$ | No | Faster |
| SOR | Modified G-S with $\omega$ | No | Tunable |

## Convergence Conditions

### Spectral Radius

:::{prf:theorem} Spectral Radius Criterion
:label: thm-spectral-radius

The iteration $\mathbf{x}^{(k+1)} = M\mathbf{x}^{(k)} + \mathbf{c}$ converges for **all** starting points if and only if $\rho(M) < 1$.

The **convergence rate** is:
$$
\|\mathbf{x}^{(k)} - \mathbf{x}^*\| \approx C \cdot \rho(M)^k
$$

Number of iterations for accuracy $\varepsilon$: $k \approx \frac{\log(1/\varepsilon)}{\log(1/\rho(M))}$
:::

### Sufficient Conditions

Computing $\rho(M)$ is expensive. These conditions are easier to check:

:::{prf:proposition} Strictly Diagonally Dominant
:label: prop-diag-dominant

If $|a_{ii}| > \sum_{j \neq i} |a_{ij}|$ for all rows $i$, then **both** Jacobi and Gauss-Seidel converge.
:::

:::{prf:proof}
:class: dropdown

Diagonal dominance implies $\|M\|_\infty < 1$, and $\rho(M) \leq \|M\|$ for any matrix norm.
:::

:::{prf:proposition} Symmetric Positive Definite Convergence
:label: prop-spd-convergence

If $A$ is symmetric positive definite, then **Gauss-Seidel converges**.

(Jacobi may or may not converge for SPD matrices.)
:::

## Connection to ODE Discretization

The same framework appears in time-stepping for ODEs!

For $\mathbf{y}' = A\mathbf{y}$, **implicit Euler** gives:

$$
\mathbf{y}_{n+1} = \mathbf{y}_n + h A\mathbf{y}_{n+1}
$$

Rearranging:

$$
(I - hA)\mathbf{y}_{n+1} = \mathbf{y}_n \implies \mathbf{y}_{n+1} = (I - hA)^{-1}\mathbf{y}_n
$$

**Stability** depends on $(I - hA)^{-1}$ being bounded — the Neumann series tells you when!

:::{prf:remark} The Stability Connection
:label: rmk-stability-connection
:class: dropdown

| Context | Fixed Point Form | Convergence/Stability |
|---------|-----------------|----------------------|
| Iterative solver | $\mathbf{x} = M\mathbf{x} + \mathbf{c}$ | $\rho(M) < 1$ |
| Implicit Euler | $\mathbf{y}_{n+1} = (I-hA)^{-1}\mathbf{y}_n$ | Eigenvalues in stability region |
| ODE equilibrium | $\mathbf{0} = A\mathbf{x}^*$ | $\text{Re}(\lambda_i) < 0$ |

The same spectral analysis underlies all three!
:::

## The Deeper Structure: Fredholm Theory

For **compact operators** $K$ (think: integral equations, discretized PDEs), the equation $(I - K)\mathbf{x} = \mathbf{c}$ has special structure:

- The spectrum of $K$ is **countable** with 0 as the only accumulation point
- $(I - K)^{-1}$ exists **except** at isolated values

When the inverse fails to exist, the **Fredholm alternative** applies: either a unique solution exists, or the kernel is **finite-dimensional** and solvability depends on $\mathbf{c}$ being orthogonal to the cokernel.

:::{admonition} The Punchline
:class: tip

The identity $I$ is the trivial operator. Everything in numerical analysis asks: **"How far can we perturb $I$ and still invert?"**

- **Neumann series** gives the explicit answer
- **Banach contraction** gives the abstract guarantee
- **Iterative methods** are the computational realization
- **Fredholm theory** describes what happens at the boundary of invertibility
:::

## Beyond Classical Methods

### Successive Over-Relaxation (SOR)

Gauss-Seidel can be accelerated with a **relaxation parameter** $\omega \in (0, 2)$:

$$
x_i^{(k+1)} = (1-\omega)x_i^{(k)} + \omega \cdot (\text{Gauss-Seidel update})
$$

With $\omega = 1$ you get standard Gauss-Seidel; $\omega > 1$ (over-relaxation) can dramatically accelerate convergence if $\omega$ is chosen well; $\omega < 1$ (under-relaxation) increases stability. The optimal $\omega$ depends on the spectrum of $M_{GS}$.

### Modern Krylov Subspace Methods

Classical methods (Jacobi, Gauss-Seidel, SOR) are **rarely used** in modern large-scale computing. **Krylov subspace methods** dominate:

| Method | For | Key Property |
|--------|-----|--------------|
| Conjugate Gradient (CG) | SPD systems | Optimal in energy norm |
| GMRES | General systems | Minimizes residual |
| BiCGSTAB | Non-symmetric | Memory efficient |

These methods search for solutions in the **Krylov subspace** $\mathcal{K}_k(A, \mathbf{b}) = \text{span}\{\mathbf{b}, A\mathbf{b}, A^2\mathbf{b}, \ldots, A^{k-1}\mathbf{b}\}$—finding good approximations using only matrix-vector products, never forming $A^{-1}$.

### When to Use What

| Situation | Method |
|-----------|--------|
| Small/medium dense $A$ | Direct (LU, QR) |
| Large sparse $A$, SPD | Conjugate Gradient |
| Large sparse $A$, general | GMRES, BiCGSTAB |
| Multiple right-hand sides | Direct (reuse factorization) |
| Only $A\mathbf{v}$ available | Krylov methods |

The classical methods matter for understanding the theory—the Neumann series perspective, the role of $\rho(M)$, the connection to fixed points and ODE stability. Modern methods build on these ideas with more sophisticated polynomial approximations in Krylov subspaces.
