# Condition Numbers and Stability for Linear Systems

:::{tip} Big Idea
The **condition number** $\kappa(A) = \|A\|\|A^{-1}\|$ tells us how sensitive a linear system is to perturbations. It's an intrinsic property of the problem—no algorithm can do better. When $\kappa(A) \gtrsim 1/\varepsilon_{\text{mach}}$, the matrix is **numerically indistinguishable from singular**.
:::

## Recap: Forward and Backward Error

We introduced [forward and backward error](forward-backward-error.md) earlier. The [golden rule](forward-backward-error.md#golden-rule) connects them:

$$
\text{relative forward error} \lesssim \kappa \times \text{relative backward error}
$$

**The key questions for linear systems:**
1. What is the condition number $\kappa(A)$?
2. When is a system ill-conditioned?
3. Which algorithms achieve small backward error?

## Sensitivity of Linear Systems

How sensitive is the solution $\mathbf{x}$ to perturbations in $A$ and $\mathbf{b}$?

For a function $f(x)$, we perturbed the input $x$. But a linear system $A\mathbf{x} = \mathbf{b}$ has **two inputs**: the matrix $A$ and the vector $\mathbf{b}$. Both are subject to errors:

- **$\mathbf{b}$ comes from measurements** — always has some noise
- **$A$ comes from a model** — coefficients may be uncertain, or stored with roundoff error

So we must understand how errors in *both* $A$ and $\mathbf{b}$ propagate to errors in $\mathbf{x}$.

:::::{ prf:theorem} Sensitivity of Linear Systems
:label: thm-sensitivity-linear-systems

For the linear system $A\mathbf{x} = \mathbf{b}$:

$$
\frac{\|\delta\mathbf{x}\|}{\|\mathbf{x}\|} \lesssim \kappa(A) \left(\frac{\|\delta A\|}{\|A\|} + \frac{\|\delta\mathbf{b}\|}{\|\mathbf{b}\|}\right)
$$

The quantity $\kappa(A) = \|A\| \|A^{-1}\|$ is the **amplification factor** from relative input perturbation to relative output error.

::::{dropdown} Derivation
**Start simple: perturb only $\mathbf{b}$.**

If $A\mathbf{x} = \mathbf{b}$, then $\mathbf{x} = A^{-1}\mathbf{b}$. Perturb $\mathbf{b} \to \mathbf{b} + \delta\mathbf{b}$:

$$
\mathbf{x} + \delta\mathbf{x} = A^{-1}(\mathbf{b} + \delta\mathbf{b}) \quad \Rightarrow \quad \delta\mathbf{x} = A^{-1}\delta\mathbf{b}
$$

Taking norms: $\|\delta\mathbf{x}\| \leq \|A^{-1}\| \|\delta\mathbf{b}\|$

For **relative error**, we want $\|\delta\mathbf{x}\|/\|\mathbf{x}\|$ in terms of $\|\delta\mathbf{b}\|/\|\mathbf{b}\|$.

The trick: use $\|\mathbf{b}\| = \|A\mathbf{x}\| \leq \|A\| \|\mathbf{x}\|$, so $1/\|\mathbf{x}\| \leq \|A\|/\|\mathbf{b}\|$.

$$
\frac{\|\delta\mathbf{x}\|}{\|\mathbf{x}\|} \leq \|A^{-1}\| \|\delta\mathbf{b}\| \cdot \frac{\|A\|}{\|\mathbf{b}\|} = \underbrace{\|A\| \|A^{-1}\|}_{\kappa(A)} \cdot \frac{\|\delta\mathbf{b}\|}{\|\mathbf{b}\|}
$$

**The condition number emerges naturally!**

**Now perturb $A$:** Suppose $(A + \delta A)(\mathbf{x} + \delta\mathbf{x}) = \mathbf{b}$. Expanding:

$$
\cancel{A\mathbf{x}} + A\delta\mathbf{x} + \delta A \cdot \mathbf{x} + \underbrace{\delta A \cdot \delta\mathbf{x}}_{\approx 0} = \cancel{\mathbf{b}}
$$

Solving: $\delta\mathbf{x} = -A^{-1}(\delta A \cdot \mathbf{x})$

Taking norms and forming relative errors:

$$
\frac{\|\delta\mathbf{x}\|}{\|\mathbf{x}\|} \leq \|A^{-1}\| \|A\| \cdot \frac{\|\delta A\|}{\|A\|} = \kappa(A) \cdot \frac{\|\delta A\|}{\|A\|}
$$
::::
:::::

## The Condition Number

The sensitivity theorem motivates the following definition:

:::{prf:definition} Condition Number
:label: def-condition-number-matrix

$$
\kappa(A) = \|A\| \cdot \|A^{-1}\|
$$

**Properties:**
- $\kappa(A) \geq 1$ always (equality for orthogonal matrices)
- $\kappa(I) = 1$
- $\kappa(A) = \infty$ if $A$ is singular (since $\sigma_{\min} = 0$)
- For the 2-norm: $\kappa_2(A) = \sigma_{\max}/\sigma_{\min}$ (ratio of largest to smallest singular value)
:::

**Rule of thumb:** Expect to lose $\log_{10}\kappa(A)$ digits of accuracy.

| Condition Number | Digits Lost |
|-----------------|-------------|
| $\kappa \approx 10^k$ | ~$k$ digits |
| $\kappa \gtrsim 1/\varepsilon_{\text{mach}} \approx 10^{16}$ | All digits |

But what does it *mean* for a matrix to be ill-conditioned? The next section provides the key insight.

## The Deep Insight: Numerically Singular Matrices

:::::{admonition} Extension: When Ill-Conditioned Means Singular (Demmel)
:class: note

A matrix with $\kappa(A) \gtrsim 1/\varepsilon_{\text{mach}}$ is **numerically indistinguishable from a singular matrix**.

::::{dropdown} Why? (Demmel's Perspective)
Consider the "distance to singularity"—how much do we need to perturb $A$ to make it singular?

For the 2-norm, this distance is exactly $\sigma_{\min}(A)$, the smallest singular value. In relative terms:

$$
\frac{\text{distance to nearest singular matrix}}{\|A\|} = \frac{\sigma_{\min}}{\sigma_{\max}} = \frac{1}{\kappa(A)}
$$

Now consider floating-point arithmetic. Every matrix $A$ is stored with relative error $\sim \varepsilon_{\text{mach}}$. The computer doesn't see $A$—it sees $A + E$ where $\|E\|/\|A\| \sim \varepsilon_{\text{mach}}$.

**If $\kappa(A) \gtrsim 1/\varepsilon_{\text{mach}}$:**
- The distance to singularity is $\lesssim \varepsilon_{\text{mach}}$
- The storage error is $\sim \varepsilon_{\text{mach}}$
- The computer cannot distinguish $A$ from a singular matrix!

This is why ill-conditioned systems are fundamentally hard—not because of bad algorithms, but because the *problem itself* is on the edge of being unsolvable.
::::
:::::

## Residuals and Backward Error

For linear systems, backward error has a simple form:

- **Residual:** $\mathbf{r} = \mathbf{b} - A\hat{\mathbf{x}}$

The computed solution $\hat{\mathbf{x}}$ exactly solves $A\hat{\mathbf{x}} = \mathbf{b} - \mathbf{r}$. The relative backward error is $\|\mathbf{r}\|/\|\mathbf{b}\|$.

:::{admonition} Warning: Small Residual ≠ Small Error
:class: warning

$$
\frac{\|\hat{\mathbf{x}} - \mathbf{x}\|}{\|\mathbf{x}\|} \lesssim \kappa(A) \cdot \frac{\|\mathbf{r}\|}{\|\mathbf{b}\|}
$$

For ill-conditioned $A$, a tiny residual can hide a large error. The residual measures backward error; you need to multiply by $\kappa(A)$ to bound forward error.
:::

## Practical Guideline: Always Check the Condition Number

Since the forward error satisfies

$$
\frac{\|\hat{\mathbf{x}} - \mathbf{x}\|}{\|\mathbf{x}\|} \lesssim \kappa(A) \cdot \varepsilon_{\text{mach}}
$$

we should **always estimate $\kappa(A)$ before trusting the solution**. But how?

### The Challenge

Computing $\kappa(A) = \|A\| \|A^{-1}\|$ exactly requires $A^{-1}$, which costs $O(n^3)$ operations—as expensive as solving the system! We need a cheaper approach.

### Hager's Algorithm: A Clever Trick

The key insight (Hager, 1984; refined by Higham) is that we can **estimate** $\|A^{-1}\|$ using only a few solves with the already-factored matrix.

:::{prf:algorithm} Hager's 1-Norm Estimator
:label: alg-hager

**Input:** LU factorization of $A$

**Output:** Estimate of $\|A^{-1}\|_1$

1. $x = \mathbf{1}/n$
2. **repeat**
3. $\qquad$ Solve $A^T y = x$ $\quad$ *(reuses the LU factorization!)*
4. $\qquad$ **if** $\|y\|_1 \leq \|y_{\text{prev}}\|_1$ **then break**
5. $\qquad$ $\xi = \text{sign}(y)$
6. $\qquad$ Solve $Az = \xi$
7. $\qquad$ $j = \arg\max_i |z_i|$
8. $\qquad$ **if** $\|z\|_\infty \leq z^T x$ **then break**
9. $\qquad$ $x = e_j$ $\quad$ *(unit vector with 1 in position $j$)*
10. **return** $\|y\|_1$
:::

**Cost:** Each iteration requires two triangular solves. Typically converges in 2–5 iterations, so total cost is $O(n^2)$—much cheaper than the $O(n^3)$ factorization.

See the [Condition Number Estimation notebook](./CondNumberEst.ipynb) for a Python implementation.

### What LAPACK Does

LAPACK's routines (e.g., `dgecon`) implement this estimation automatically. When you call `np.linalg.cond(A)` or use SciPy's linear solvers with condition estimation, this is what happens behind the scenes.

:::{admonition} Practical Workflow
:class: tip

**Before solving $Ax = b$:**
1. Factor $A = LU$ (or $A = QR$)
2. Estimate $\kappa(A)$ using Hager's algorithm — costs only $O(n^2)$ extra
3. Check: if $\kappa(A) \cdot \varepsilon_{\text{mach}} \gtrsim 1$, the answer may be meaningless
4. Solve the system
5. The forward error bound is $\lesssim \kappa(A) \cdot \varepsilon_{\text{mach}}$

This gives you a **reliable error bound** essentially for free!
:::

## Stability of Algorithms: A Preview

An algorithm is **backward stable** if it produces a solution with backward error $\sim \varepsilon_{\text{mach}}$. Combined with the golden rule:

$$
\text{forward error} \lesssim \kappa(A) \cdot \varepsilon_{\text{mach}}
$$

This is the best we can hope for—any algorithm must contend with the condition number.

:::{admonition} What About Forward Stability?
:class: note

An algorithm is **forward stable** if it achieves forward error $\lesssim \kappa(A) \cdot \varepsilon_{\text{mach}}$ directly.

Backward stability *implies* forward stability (via the golden rule), but is a stronger guarantee: it tells us the computed answer is the *exact* answer to a nearby problem. This is more useful because:
- It separates algorithm quality from problem sensitivity
- It applies even when we can't compute the true answer
:::

**Coming up:** We'll see that:
- **Householder QR** is backward stable (the gold standard)
- **LU with partial pivoting** is backward stable in practice (with caveats)
- **Classical Gram-Schmidt** is *not* stable—orthogonality loss scales with $\kappa(A)$
