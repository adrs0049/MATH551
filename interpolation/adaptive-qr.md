# Adaptive QR: Automatic Accuracy for Spectral Methods

:::{tip} Big Idea
When solving differential equations with spectral methods, how do we know how many basis functions (Chebyshev polynomials) to use? Too few gives a bad answer; too many wastes computation. The **adaptive QR algorithm** elegantly solves both problems at once—it finds the right number of terms *while* computing the solution, and tells us exactly when it's accurate enough.
:::

## The Problem: How Many Terms?

Recall that spectral methods approximate the solution to a boundary value problem as:
$$
u(x) \approx \sum_{k=0}^{n-1} c_k T_k(x)
$$

where $T_k$ are Chebyshev polynomials and $c_k$ are unknown coefficients we need to find.

**The challenge:** We don't know $n$ in advance!

- Choose $n$ too small → the solution is inaccurate
- Choose $n$ too large → we waste computation solving a bigger system

### The Naive Approach: Double and Retry

The obvious strategy:
1. Try $n = 16$, solve the system, check if accurate enough
2. If not, try $n = 32$, solve again from scratch
3. If not, try $n = 64$, solve again from scratch
4. ...

**Problem:** Each time we double $n$, we throw away all the work we did before!

:::{prf:remark} Why This Matters
For smooth solutions, spectral methods converge *exponentially* fast—often $n = 20$ to $50$ is enough for 14 digits of accuracy. But finding the right $n$ by trial-and-error means solving multiple systems, which is wasteful.
:::

## The Key Insight: Build the Solution Incrementally

Instead of solving separate problems for $n = 16$, then $n = 32$, etc., the adaptive QR algorithm:

1. Starts solving with $n = 1$ term
2. Adds one term at a time, updating the solution
3. Checks at each step: "Is this accurate enough?"
4. Stops as soon as the answer is "yes"

**The magic:** Step 2 can be done efficiently using **Givens rotations**, and step 3 gives us the residual *for free* as a byproduct of the computation!

## From ODE to Linear System

When we discretize a boundary value problem using spectral collocation, we get a linear system:
$$
A \mathbf{c} = \mathbf{b}
$$

where:
- $\mathbf{c} = [c_0, c_1, \ldots, c_{n-1}]^T$ are the Chebyshev coefficients we want
- $A$ encodes the differential operator (e.g., second derivative) in coefficient space
- $\mathbf{b}$ encodes the right-hand side and boundary conditions

### A Special Structure: "Almost Banded"

Here's the key observation that makes everything work. The matrix $A$ has a very special structure:

```
┌─────────────────────────────────────────┐
│   K rows from boundary conditions       │  ← Dense rows (usually 2)
├─────────────────────────────────────────┤
│                                         │
│        Banded part from the ODE         │  ← Only a few diagonals
│                                         │
└─────────────────────────────────────────┘
```

- The **top few rows** (one per boundary condition) are dense
- The **rest of the matrix** is banded—only a few diagonals are nonzero

:::{prf:definition} Almost Banded Matrix
:label: def-almost-banded

A matrix is **almost banded** if it has the form:
$$
A = \begin{pmatrix} B \\ C \end{pmatrix}
$$
where $B$ has $K$ dense rows (from boundary conditions) and $C$ is banded (from the differential operator).
:::

**Why is this useful?** Banded matrices can be factored in $O(n)$ operations instead of $O(n^3)$. The "almost banded" structure preserves most of this efficiency!

## QR Factorization: A Quick Review

Recall that any matrix $A$ can be factored as:
$$
A = QR
$$
where $Q$ is orthogonal ($Q^T Q = I$) and $R$ is upper triangular.

To solve $A\mathbf{c} = \mathbf{b}$:
1. Compute $QR = A$
2. Compute $Q^T \mathbf{b}$ (easy, just matrix-vector multiply)
3. Solve $R\mathbf{c} = Q^T \mathbf{b}$ by back-substitution

**The key property of $Q$:** Since $Q$ is orthogonal, it preserves lengths:
$$
\|Q\mathbf{x}\|_2 = \|\mathbf{x}\|_2 \quad \text{for any vector } \mathbf{x}
$$

This will be crucial for understanding the residual!

## Givens Rotations: Building QR One Column at a Time

A **Givens rotation** is a simple $2 \times 2$ orthogonal matrix that zeros out one element:

$$
G = \begin{pmatrix} c & s \\ -s & c \end{pmatrix}
$$

where $c^2 + s^2 = 1$ (so $c = \cos\theta$ and $s = \sin\theta$ for some angle $\theta$).

When we apply $G$ to a 2-vector:

$$
\begin{pmatrix} c & s \\ -s & c \end{pmatrix} \begin{pmatrix} a \\ b \end{pmatrix} = \begin{pmatrix} \sqrt{a^2 + b^2} \\ 0 \end{pmatrix}
$$

By choosing $c = a/\sqrt{a^2 + b^2}$ and $s = b/\sqrt{a^2 + b^2}$, we can zero out $b$ while preserving the length $\sqrt{a^2 + b^2}$.

### Using Givens to Zero a Column

To create zeros below the diagonal in column $j$ of a matrix:

```
Before column j:          After applying Givens:
     column j                  column j
        ↓                         ↓
[  × × × × ×  ]           [  × × × × ×  ]
[  0 × × × ×  ]           [  0 × × × ×  ]
[  0 0 × × ×  ]    →      [  0 0 × × ×  ]   ← R[j,j]
[  0 0 a × ×  ]           [  0 0 0 × ×  ]   ← zeroed!
[  0 0 b × ×  ]           [  0 0 0 × ×  ]   ← zeroed!
```

We apply a sequence of Givens rotations, each zeroing one entry below the diagonal.

**Important:** We apply the same rotations to the right-hand side $\mathbf{b}$!

## The Adaptive QR Algorithm

Now we can explain the algorithm. Here's the beautiful idea:

### Step 1: Process Columns One at a Time

Instead of building the full matrix and then factoring it, we:
1. Start with just column 1
2. Apply Givens rotations to make it upper triangular
3. Add column 2, apply more Givens rotations
4. Add column 3, apply more Givens rotations
5. ...continue...

At each step, we have a partial QR factorization of the columns processed so far.

### Step 2: Check the Residual (For Free!)

Here's the magic. After processing $j$ columns, the transformed system looks like:

$$
Q_j A = \begin{pmatrix} R_j & * \\ 0 & * \end{pmatrix}, \quad
Q_j \mathbf{b} = \begin{pmatrix} \mathbf{r}_{1:j} \\ \mathbf{r}_{j+1:\text{end}} \end{pmatrix}
$$

- The top $j$ rows form an upper triangular system
- The bottom rows have zeros in the first $j$ columns

**The key theorem:**

:::{prf:theorem} Residual from Transformed RHS
:label: thm-residual-qr

If we truncate the solution to $j$ terms (setting $c_k = 0$ for $k \geq j$) and solve the top $j \times j$ system, then:
$$
\|\text{residual}\|_2 = \|\mathbf{r}_{j+1:\text{end}}\|_2
$$

The residual norm equals the norm of the "leftover" components of the transformed right-hand side!
:::

:::{prf:proof}
:class: dropdown

Let $\hat{\mathbf{c}}$ be the solution with only $j$ nonzero coefficients, found by solving $R_j \hat{\mathbf{c}}_{1:j} = \mathbf{r}_{1:j}$.

The residual is:
$$
\text{residual} = \mathbf{b} - A\hat{\mathbf{c}}
$$

Multiply both sides by $Q_j$ (which preserves length since $Q_j$ is orthogonal):
$$
Q_j \cdot \text{residual} = Q_j \mathbf{b} - Q_j A \hat{\mathbf{c}}
$$

The right-hand side is:

$$
\begin{pmatrix} \mathbf{r}_{1:j} \\ \mathbf{r}_{j+1:\text{end}} \end{pmatrix} - \begin{pmatrix} R_j & * \\ 0 & * \end{pmatrix} \begin{pmatrix} \hat{\mathbf{c}}_{1:j} \\ 0 \end{pmatrix}
= \begin{pmatrix} \mathbf{r}_{1:j} - R_j \hat{\mathbf{c}}_{1:j} \\ \mathbf{r}_{j+1:\text{end}} \end{pmatrix}
= \begin{pmatrix} \mathbf{0} \\ \mathbf{r}_{j+1:\text{end}} \end{pmatrix}
$$

The top part is zero because $\hat{\mathbf{c}}_{1:j}$ was chosen to satisfy $R_j \hat{\mathbf{c}}_{1:j} = \mathbf{r}_{1:j}$.

Since $Q_j$ preserves lengths:
$$
\|\text{residual}\|_2 = \|Q_j \cdot \text{residual}\|_2 = \|\mathbf{r}_{j+1:\text{end}}\|_2
$$
:::

### Step 3: Stop When Accurate Enough

The algorithm becomes:

```
for j = 1, 2, 3, ...
    1. Add column j to the system
    2. Apply Givens rotations to zero entries below diagonal
    3. Apply same rotations to right-hand side b
    4. Check: is ||r_{j+1:end}|| < tolerance?
       - If YES: stop, n_opt = j
       - If NO: continue to j+1

Solve R_{n_opt} c = r_{1:n_opt} by back-substitution
```

**We get the optimal $n$ and the solution simultaneously!**

## Why "Almost Banded" Matters

Remember the matrix structure:
```
┌─────────────────────────────────────────┐
│   K dense boundary rows                 │
├─────────────────────────────────────────┤
│        Banded part (bandwidth m)        │
└─────────────────────────────────────────┘
```

When we apply Givens rotations:
- Each column only has about $m$ entries below the diagonal (because of banding)
- So we only need about $m$ Givens rotations per column
- Total work: $O(m^2 n)$ instead of $O(n^3)$

The dense boundary rows cause some "fill-in" (new nonzeros created), but this is limited and can be tracked efficiently.

:::{prf:remark} Memory Efficiency
The fill-in from boundary rows affects at most $K$ columns of each row, where $K$ is the number of boundary conditions (usually 2). The algorithm stores only the banded part plus this limited fill-in, using $O(mn)$ memory instead of $O(n^2)$.
:::

## Example: Second-Order BVP

Consider the problem:
$$
u''(x) + u(x) = f(x), \quad u(-1) = \alpha, \quad u(1) = \beta
$$

The almost-banded structure arises because:
1. **Boundary rows:** $u(-1) = \alpha$ and $u(1) = \beta$ involve all Chebyshev coefficients (dense)
2. **Interior rows:** The relation between Chebyshev coefficients of $u''$ and $u$ is local (banded)

For the ultraspherical spectral method (a sophisticated variant), the bandwidth is small:
- Second derivative: bandwidth 3
- Fourth derivative: bandwidth 5

This means we can solve BVPs with spectral accuracy in just $O(n)$ operations!

## Putting It All Together

The adaptive QR algorithm combines several ideas we've seen:

| Concept | Role in Algorithm |
|---------|------------------|
| Chebyshev polynomials | Basis functions with exponential convergence |
| Spectral collocation | Turns ODE into linear system |
| Almost-banded structure | Enables $O(n)$ factorization |
| QR factorization | Stable way to solve the system |
| Givens rotations | Builds QR incrementally |
| Orthogonal transformations preserve norms | Gives us the residual for free |

**The result:** An algorithm that automatically finds how many Chebyshev coefficients are needed, solves the system, and verifies the accuracy—all in a single pass through the columns!

## Practical Impact

This algorithm (from Olver & Townsend, 2013) is implemented in software like Chebfun and enables:

- Solving ODEs to 15 digits of accuracy automatically
- Handling problems where the required resolution varies (e.g., boundary layers)
- Systems of coupled ODEs with the same framework

:::{prf:example} Automatic Accuracy
:class: dropdown

Consider solving $u'' = e^{4x}$ on $[-1, 1]$ with $u(\pm 1) = 0$.

With traditional spectral methods, you might:
1. Try $n = 16$, check residual: too big
2. Try $n = 32$, check residual: okay!

With adaptive QR:
1. Process columns 1, 2, 3, ...
2. At column 24, residual drops below $10^{-14}$
3. Stop and return the solution

No wasted work, and the residual check was free!
:::

## Looking Ahead

This algorithm illustrates a powerful principle in scientific computing: **structure-exploiting algorithms**. By understanding the mathematical structure of the problem (almost-banded matrices, orthogonal transformations, coefficient decay), we can design algorithms that are both:

- **Fast:** $O(n)$ instead of $O(n^3)$
- **Reliable:** Automatic accuracy verification

This same philosophy—understanding structure to design better algorithms—appears throughout numerical analysis and is what makes scientific computing both an art and a science.

## References

- Olver, S. and Townsend, A. (2013). "A Fast and Well-Conditioned Spectral Method." *SIAM Review*, 55(3), 462-489.

:::{seealso}
- [Chebyshev Polynomials](chebyshev.md) - The basis functions
- [Spectral Methods](spectral-methods.md) - Collocation for BVPs
- [Orthogonality and Projections](../qr-least-squares/orthogonality.md) - Why orthogonal matrices preserve norms
:::
