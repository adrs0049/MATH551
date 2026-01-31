# Core Concepts and Fundamentals

:::{tip} Big Idea
These are the fundamental concepts that unify the entire course: approximation theorems justify *what* we're doing, while error analysis tells us *how well* we're doing it.
:::

## Foundational Theorems

:::{prf:theorem} Taylor's Theorem
:label: thm-taylor-core

If $f$ is $(n+1)$-times differentiable, then:
$$
f(x) = \sum_{k=0}^{n} \frac{f^{(k)}(a)}{k!}(x-a)^k + \frac{f^{(n+1)}(\xi)}{(n+1)!}(x-a)^{n+1}
$$
for some $\xi$ between $a$ and $x$.

**Why it matters:** This is THE foundation of numerical analysis.
- **Truncation error:** The remainder term tells us the error when we truncate
- **Finite differences:** Derived by Taylor expanding $f(x+h)$
- **Newton's method:** Comes from linearizing $f$ (first-order Taylor)
- **Convergence order:** Higher-order methods use more Taylor terms
:::

:::{prf:theorem} Weierstrass Approximation Theorem
:label: thm-weierstrass-core

Let $f$ be continuous on $[a, b]$. For any $\varepsilon > 0$, there exists a polynomial $p$ such that:
$$
\|f - p\|_\infty < \varepsilon
$$

**Why it matters:** Justifies the entire polynomial approximation enterprise. Continuous functions CAN be approximated by polynomials. The question is how efficiently.

:::{warning}
Weierstrass guarantees existence but says nothing about:
- Which polynomial degree is needed
- Which nodes to use for interpolation
- How to construct the polynomial efficiently

These are the questions numerical analysis answers.
:::

:::{prf:theorem} Best Approximation Theorem
:label: thm-best-approximation-core

Let $V$ be an inner product space and $W \subset V$ a finite-dimensional subspace. For any $v \in V$, there exists a unique $w^* \in W$ minimizing $\|v - w\|$.

The best approximation is characterized by **orthogonality**:
$$
\langle v - w^*, w \rangle = 0 \quad \text{for all } w \in W
$$

**Why it matters:**
- **Least squares:** $A^T(b - Ax) = 0$ is exactly this orthogonality condition
- **Fourier/Chebyshev:** Coefficients come from projecting onto basis functions
- **Why QR works:** Orthogonal projection onto $\text{Col}(A)$
:::

## Forward and Backward Error

:::{prf:definition} Forward and Backward Error
:label: def-forward-backward-error

For a problem $f: X \to Y$ and computed solution $\hat{y}$:

**Forward error:** How far is $\hat{y}$ from the true answer?
$$
\text{forward error} = \|\hat{y} - f(x)\|
$$

**Backward error:** For what perturbed input would $\hat{y}$ be the exact answer?
$$
\text{backward error} = \min \{ \|\delta x\| : f(x + \delta x) = \hat{y} \}
$$
:::

:::{prf:definition} Condition Number
:label: def-condition-number-core

The **condition number** measures sensitivity of the problem to input perturbations:
$$
\kappa = \lim_{\delta \to 0} \sup_{\|\delta x\| \leq \delta} \frac{\|f(x + \delta x) - f(x)\| / \|f(x)\|}{\|\delta x\| / \|x\|}
$$

For an invertible matrix $A$:
$$
\kappa(A) = \|A\| \cdot \|A^{-1}\|
$$

For the 2-norm: $\kappa_2(A) = \sigma_{\max}(A) / \sigma_{\min}(A)$.
:::

:::{prf:theorem} The Fundamental Error Relationship
:label: thm-fundamental-error

$$
\frac{\|\text{forward error}\|}{\|\text{true solution}\|} \leq \kappa \cdot \frac{\|\text{backward error}\|}{\|\text{input}\|}
$$

Or more simply:
$$
\text{relative forward error} \lesssim \kappa \times \text{relative backward error}
$$

**Backward stable algorithm + well-conditioned problem = accurate answer.**
:::

:::{prf:remark} The Error Analysis Philosophy
:label: rmk-error-philosophy

This decomposition separates concerns:
- **Algorithm designers** aim for backward stability (small backward error)
- **Problem formulators** aim for well-conditioned problems (small $\kappa$)
- **Users** should check conditioning before trusting results

A backward stable algorithm gives the exact answer to a *slightly wrong* problem. If the problem is well-conditioned, "slightly wrong problem" means "nearly correct answer."
:::

## The Contraction Mapping Principle

:::{prf:theorem} Banach Fixed Point Theorem
:label: thm-banach-core

Let $(X, d)$ be a **complete metric space** and $T: X \to X$ a **contraction mapping**, i.e., there exists $q \in [0, 1)$ such that
$$
d(T(x), T(y)) \leq q \cdot d(x, y) \quad \text{for all } x, y \in X.
$$

Then:
1. $T$ has a **unique fixed point** $x^* \in X$
2. For any starting point $x_0 \in X$, the iteration $x_{n+1} = T(x_n)$ converges to $x^*$
3. The convergence is **geometric**: $d(x_n, x^*) \leq q^n \cdot d(x_0, x^*)$
4. **A priori bound:** $d(x_n, x^*) \leq \frac{q^n}{1-q} d(x_0, x_1)$
5. **A posteriori bound:** $d(x_n, x^*) \leq \frac{q}{1-q} d(x_{n-1}, x_n)$
:::

:::{prf:proof}
:class: dropdown

**Existence and convergence:** The sequence $\{x_n\}$ is Cauchy. For $m > n$:
$$
d(x_n, x_m) \leq \sum_{k=n}^{m-1} d(x_k, x_{k+1}) \leq \sum_{k=n}^{m-1} q^k d(x_0, x_1) \leq \frac{q^n}{1-q} d(x_0, x_1)
$$
Since $q < 1$, this tends to 0 as $n \to \infty$. By completeness, $x^* = \lim_{n \to \infty} x_n$ exists.

Taking limits in $x_{n+1} = T(x_n)$ and using continuity of $T$ (contractions are Lipschitz, hence continuous):
$$
x^* = \lim_{n \to \infty} x_{n+1} = \lim_{n \to \infty} T(x_n) = T\left(\lim_{n \to \infty} x_n\right) = T(x^*)
$$

**Uniqueness:** If $T(x^*) = x^*$ and $T(y^*) = y^*$, then
$$
d(x^*, y^*) = d(T(x^*), T(y^*)) \leq q \cdot d(x^*, y^*)
$$
Since $q < 1$, this forces $d(x^*, y^*) = 0$, so $x^* = y^*$.

**Convergence rate:** $d(x_n, x^*) = d(T(x_{n-1}), T(x^*)) \leq q \cdot d(x_{n-1}, x^*) \leq q^n d(x_0, x^*)$.

**A posteriori bound:** Let $m \to \infty$ in the Cauchy estimate with $n$ replaced by $n-1$.
:::

:::{prf:remark} Where We Use the Banach Fixed Point Theorem
:label: rmk-banach-applications

This theorem is THE foundation for iterative methods throughout scientific computing:

| Application | Space $X$ | Contraction $T$ |
|-------------|-----------|-----------------|
| **Fixed-point iteration** | $\mathbb{R}$ or $[a,b]$ | $g(x)$ with $\|g'\|_\infty < 1$ |
| **Newton's method** (local) | Neighborhood of root | Newton map near simple root |
| **Picard iteration for ODEs** | $C([0,T])$ | Integral operator |
| **Iterative linear solvers** | $\mathbb{R}^n$ | $x \mapsto Mx + c$ with $\rho(M) < 1$ |
| **Implicit function theorem** | Banach space | Contraction from IFT proof |

The a posteriori bound is particularly useful: it gives a computable error estimate from consecutive iterates.
:::

## The Neumann Series

:::{prf:theorem} Neumann Series
:label: thm-neumann-core

Let $\mathcal{A}$ be a **unital Banach algebra** with identity $1$ and norm $\|\cdot\|$. If $a \in \mathcal{A}$ satisfies $\|a\| < 1$, then $(1 - a)$ is invertible and
$$
(1 - a)^{-1} = \sum_{k=0}^{\infty} a^k = 1 + a + a^2 + a^3 + \cdots
$$

More generally, if $\mathcal{B}(X)$ denotes the bounded linear operators on a Banach space $X$, and $A \in \mathcal{B}(X)$ has **spectral radius** $\rho(A) < 1$, then $(I - A)^{-1} = \sum_{k=0}^\infty A^k$.
:::

:::{prf:proof}
:class: dropdown

**Convergence:** Since $\|a^k\| \leq \|a\|^k$ (submultiplicativity of the norm in a Banach algebra) and $\|a\| < 1$, the series $\sum_{k=0}^\infty a^k$ converges absolutely, hence converges in the complete space $\mathcal{A}$.

**Verification:** The partial sums $S_n = \sum_{k=0}^n a^k$ satisfy
$$
(1 - a) S_n = S_n (1 - a) = 1 - a^{n+1}
$$
As $n \to \infty$, $\|a^{n+1}\| \leq \|a\|^{n+1} \to 0$, so $(1 - a) \sum_{k=0}^\infty a^k = 1$.

**Spectral radius version:** Gelfand's formula gives $\rho(A) = \lim_{k \to \infty} \|A^k\|^{1/k}$. If $\rho(A) < 1$, then for any $r$ with $\rho(A) < r < 1$, we have $\|A^k\| \leq C r^k$ for large $k$, ensuring convergence.
:::

:::{prf:remark} Where We Use the Neumann Series
:label: rmk-neumann-applications

The Neumann series is the operator-theoretic generalization of the geometric series $\frac{1}{1-r} = \sum_{k=0}^\infty r^k$ for $|r| < 1$.

| Application | Banach Algebra | Element $a$ |
|-------------|----------------|-------------|
| **Iterative linear solvers** | $\mathbb{R}^{n \times n}$ | Iteration matrix $M$ |
| **Perturbation of inverses** | $\mathcal{B}(X)$ | $A^{-1}\Delta A$ for perturbed $A + \Delta A$ |
| **Integral equations** | $\mathcal{B}(L^2)$ | Integral operator $K$ |
| **ODE stability** | $\mathbb{C}^{n \times n}$ | $hA$ in implicit methods |
| **Resolvent** | $\mathcal{B}(X)$ | $\lambda^{-1}A$ for $|\lambda| > \rho(A)$ |

**Key insight:** The condition $\rho(M) < 1$ (not $\|M\| < 1$) is necessary and sufficient for convergence. This is why spectral radius, not norm, controls iterative methods.
:::

:::{prf:corollary} Perturbation of Inverses
:label: cor-perturbation-inverse

If $A$ is invertible and $\|A^{-1}\| \cdot \|\Delta A\| < 1$, then $A + \Delta A$ is invertible and
$$
\|(A + \Delta A)^{-1} - A^{-1}\| \leq \frac{\|A^{-1}\|^2 \|\Delta A\|}{1 - \|A^{-1}\| \|\Delta A\|}
$$
:::

:::{prf:proof}
:class: dropdown

Write $A + \Delta A = A(I + A^{-1}\Delta A)$. Since $\|A^{-1}\Delta A\| \leq \|A^{-1}\| \|\Delta A\| < 1$, the Neumann series gives
$$
(A + \Delta A)^{-1} = (I + A^{-1}\Delta A)^{-1} A^{-1} = \sum_{k=0}^\infty (-A^{-1}\Delta A)^k A^{-1}
$$
The bound follows from estimating $\|(I + A^{-1}\Delta A)^{-1} - I\|$.
:::

