# Fixed Point Iteration

:::{tip} Big Idea
Root finding problems can be reformulated as fixed point problems. The key insight: if $|g'(x)| < 1$ near a fixed point, the iteration $x_{n+1} = g(x_n)$ converges. This provides a unified framework for analyzing iterative methods.
:::

## Fixed Points

:::{prf:definition} Fixed Point
:label: def-fixed-point

$x_0 \in \mathbb{R}$ is a **fixed point** of $g(x)$ if $g(x_0) = x_0$.
:::

The basic idea: reformulate a root-finding problem as a fixed-point problem.

$$
\text{Root finding: } f(x) = 0 \quad \longrightarrow \quad \text{Fixed point: } x = g(x)
$$

There are many ways to do this! Given $f(x) = 0$, you could write:
- $g(x) = x - f(x)$
- $g(x) = x - cf(x)$ for any nonzero $c$
- $g(x) = x - f(x)/f'(x)$ (this gives Newton's method!)

**The choice of $g$ matters enormously** for convergence—as we'll see.

## The Algorithm

The fixed-point iteration is beautifully simple.

:::{prf:algorithm} Fixed Point Iteration
:label: alg-fixed-point

**Input:** Function $g$, initial guess $x_0$, tolerance $\varepsilon$, max iterations $N$

**Output:** Approximate fixed point $x$

1. **for** $n = 0, 1, 2, \ldots, N-1$:
2. $\qquad x_{n+1} \gets g(x_n)$
3. $\qquad$ **if** $|x_{n+1} - x_n| < \varepsilon$: **return** $x_{n+1}$
4. **return** $x_N$ *(or indicate failure)*
:::

That's it. But whether this converges—and how fast—depends entirely on properties of $g$.

## Why the Choice of $g$ Matters

Consider finding the root of $f(x) = x^2 - 3 = 0$ (i.e., finding $\sqrt{3}$).

Here are three valid reformulations:

- **G1:** From $x^2 = 3$, we get $g_1(x) = 3/x$
- **G2:** Add $x$ to both sides: $g_2(x) = x - (x^2 - 3)/2$
- **G3:** Divide by $2x$: $g_3(x) = (x^2 + 3)/(2x)$

All three have $\sqrt{3}$ as a fixed point. But their behavior is dramatically different:
- **G1** cycles forever, never converging
- **G2** converges slowly (linearly)
- **G3** converges rapidly (quadratically—this is Newton's method!)

:::{seealso}
[Fixed Point Iteration Demo](fixed-point-iteration.ipynb) — Compare G1, G2, G3 convergence behavior
:::

## Existence and Uniqueness

When does a fixed point exist? When is it unique?

:::{prf:theorem} Fixed Point Existence and Uniqueness
:label: thm-fixed-point-existence

Given $g: [a, b] \to \mathbb{R}$:

**Existence:** If $g$ is continuous and maps $[a,b]$ into itself (i.e., $g([a,b]) \subseteq [a,b]$), then $g$ has at least one fixed point in $[a,b]$.

**Uniqueness:** If additionally $|g'(x)| \leq \rho < 1$ for all $x \in [a,b]$, then the fixed point is unique.
:::

:::{prf:proof} Existence
:class: dropdown

Define $h(x) = x - g(x)$.

If $g(a) = a$ or $g(b) = b$, we're done—we've found a fixed point.

Otherwise, since $g([a,b]) \subseteq [a,b]$:
- $g(a) > a$, so $h(a) = a - g(a) < 0$
- $g(b) < b$, so $h(b) = b - g(b) > 0$

By the Intermediate Value Theorem, there exists $c \in (a,b)$ with $h(c) = 0$, i.e., $g(c) = c$.
:::

:::{prf:proof} Uniqueness
:class: dropdown

Suppose two fixed points $c_1 < c_2$ exist. By the Mean Value Theorem:

$$
|c_1 - c_2| = |g(c_1) - g(c_2)| = |g'(\xi)||c_1 - c_2| \leq \rho|c_1 - c_2|
$$

for some $\xi \in (c_1, c_2)$.

This implies $(1-\rho)|c_1 - c_2| \leq 0$. But $\rho < 1$ and $c_1 \neq c_2$, so $(1-\rho)|c_1 - c_2| > 0$—a contradiction.
:::

## Convergence

:::{prf:theorem} Convergence of Fixed Point Iteration
:label: thm-fixed-point-convergence

If $g([a,b]) \subseteq [a,b]$ and $|g'(x)| \leq \rho < 1$ on $[a,b]$, then for any $x_0 \in [a,b]$, the sequence $x_{n+1} = g(x_n)$ converges to the unique fixed point $c$.

Moreover, the convergence is geometric:
$$
|x_n - c| \leq \rho^n |x_0 - c|
$$
:::

:::{prf:proof}
:class: dropdown

Since $c$ is a fixed point, $g(c) = c$. Using the Mean Value Theorem:

$$
|x_{n+1} - c| = |g(x_n) - g(c)| = |g'(\xi_n)||x_n - c| \leq \rho|x_n - c|
$$

Applying this recursively:

$$
|x_n - c| \leq \rho|x_{n-1} - c| \leq \rho^2|x_{n-2} - c| \leq \cdots \leq \rho^n|x_0 - c|
$$

Since $\rho < 1$, we have $\rho^n \to 0$, so $x_n \to c$.
:::

:::{prf:remark} Local Convergence
:label: rmk-local-convergence

Even if $|g'(x)| < 1$ only *at* the fixed point (not on the whole interval), the iteration still converges—provided we start close enough. This is called **local convergence**.

Specifically: if $g \in \mathcal{C}^1$, $g(c) = c$, and $|g'(c)| < 1$, then there exists $\delta > 0$ such that the iteration converges for any $x_0$ with $|x_0 - c| < \delta$.
:::

## The Derivative at the Fixed Point

The key insight: **$|g'(c)|$ determines everything**.

For our three reformulations of $x^2 - 3 = 0$:

- $g_1(x) = 3/x$: We have $g_1'(x) = -3/x^2$, so $|g_1'(\sqrt{3})| = 1$. Right on the boundary—no convergence guaranteed (and indeed, it fails).

- $g_2(x) = x - (x^2-3)/2$: We have $g_2'(x) = 1 - x$, so $|g_2'(\sqrt{3})| = |1 - \sqrt{3}| \approx 0.73$. Linear convergence with rate $\rho \approx 0.73$.

- $g_3(x) = (x^2+3)/(2x)$: We have $g_3'(x) = 1/2 - 3/(2x^2)$, so $g_3'(\sqrt{3}) = 0$. The derivative vanishes—this signals faster-than-linear convergence.

## Order of Convergence

:::{prf:definition} Order of Convergence
:label: def-convergence-order

A sequence $\{x_n\}$ converging to $\ell$ has **order** $p$ if:

$$
\lim_{n\to\infty} \frac{|x_{n+1} - \ell|}{|x_n - \ell|^p} = C
$$

for some constant $C > 0$.

- $p = 1$: **Linear** convergence (error shrinks by constant factor)
- $p = 2$: **Quadratic** convergence (digits of accuracy double each step)
:::

:::{prf:theorem} Order of Fixed Point Iteration
:label: thm-fixed-point-order

If $g(c) = c$ and $g'(c) = g''(c) = \cdots = g^{(p-1)}(c) = 0$ but $g^{(p)}(c) \neq 0$, then the iteration converges with order $p$.
:::

:::{prf:proof}
:class: dropdown

Taylor expand $g(x_n)$ around $c$:

$$
x_{n+1} = g(x_n) = g(c) + g'(c)(x_n - c) + \cdots + \frac{g^{(p)}(\xi)}{p!}(x_n - c)^p
$$

Since $g(c) = c$ and the first $p-1$ derivatives vanish:

$$
x_{n+1} - c = \frac{g^{(p)}(\xi)}{p!}(x_n - c)^p
$$

Thus $|x_{n+1} - c| \approx C|x_n - c|^p$ with $C = |g^{(p)}(c)|/p!$.
:::

This explains why $g_3$ converges so fast: $g_3'(\sqrt{3}) = 0$ means at least quadratic convergence.

## The Banach Fixed Point Theorem

The convergence results above are special cases of a fundamental principle that appears throughout mathematics.

(theoretical-foundation-the-banach-fixed-point-theorem)=
:::{prf:theorem} Banach Fixed Point Theorem
:label: thm-banach

Let $(X, d)$ be a **complete metric space** and $T: X \to X$ be a **contraction**:

$$
d(T(x), T(y)) \leq q \cdot d(x, y) \quad \text{for all } x, y \in X
$$

for some $q < 1$. Then:
1. $T$ has a **unique fixed point** $x^*$
2. The iteration $x_{n+1} = T(x_n)$ **converges** from any starting point
3. Convergence is **geometric**: $d(x_n, x^*) \leq q^n \cdot d(x_0, x^*)$
:::

:::{prf:remark} Why the Banach FPT Matters
:label: rmk-banach-applications
:class: dropdown

The Banach FPT is not just about scalar equations. The same principle governs:

- **Newton's method for systems**: The iteration is a contraction near the solution
- **Picard iteration for ODEs**: Proves existence and uniqueness for $y' = f(t,y)$
- **Iterative linear solvers**: Jacobi and Gauss-Seidel converge when the iteration matrix is a contraction

Whenever you see an iterative method that "works," there's often a contraction hiding underneath.
:::

## Advantages and Disadvantages

**Advantages:**
- **Simple** — Just iterate $x_{n+1} = g(x_n)$
- **Unifying framework** — Newton's method is a fixed-point iteration in disguise
- **Flexible** — Many choices for $g$; can design for fast convergence
- **Generalizes** — Extends to systems, ODEs, infinite dimensions (Banach FPT)

**Disadvantages:**
- **Not guaranteed** — Can diverge if $|g'(c)| \geq 1$
- **Sensitive** — Different reformulations give wildly different behavior
- **Slow** — Linear convergence when $|g'(c)|$ is close to 1

**Design principle:** Make $|g'(c)|$ small. If $g'(c) = 0$, you get quadratic convergence—this is Newton's method.
