# Newton's Method

:::{tip} Big Idea
Newton's method approximates a function by its tangent line at each iteration. This simple idea yields quadratic convergence—the number of correct digits roughly doubles each step. It's the workhorse of nonlinear equation solving.
:::

## Derivation

Given $f(x) = 0$ and an initial guess $x_0$, approximate $f$ by its tangent line at $x_0$:

$$
y = f(x_0) + f'(x_0)(x - x_0)
$$

Setting $y = 0$ and solving for $x$ gives the next approximation:

$$
x_1 = x_0 - \frac{f(x_0)}{f'(x_0)}
$$

Repeating this process:

$$
x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}
$$

## Newton's Method as Fixed Point Iteration

Newton's method is a fixed point iteration with:

$$
g(x) = x - \frac{f(x)}{f'(x)}
$$

If $x_n \to c$, then:
$$
c = c - \frac{f(c)}{f'(c)} \implies f(c) = 0
$$

So the fixed point is indeed a root of $f$.

## The Algorithm

:::{prf:algorithm} Newton's Method
:label: alg-newton

**Input:** Functions $f$ and $f'$, initial guess $x_0$, tolerance $\varepsilon$, max iterations $N$

**Output:** Approximate root $x$

1. **for** $n = 0, 1, 2, \ldots, N-1$:
2. $\qquad x_{n+1} \gets x_n - f(x_n)/f'(x_n)$
3. $\qquad$ **if** $|x_{n+1} - x_n| < \varepsilon$: **return** $x_{n+1}$
4. **return** $x_N$ *(or indicate failure)*
:::

## Convergence Analysis

:::{prf:theorem} Local Convergence of Newton's Method
:label: thm-newton-convergence

Suppose $f \in \mathcal{C}^2([a,b])$, $f(c) = 0$ for some $c \in (a,b)$, and $c$ is a **simple root** (i.e., $f'(c) \neq 0$). Then for $x_0$ sufficiently close to $c$, Newton's method converges to $c$.
:::

:::{prf:proof}
:class: dropdown

The iteration function is $g(x) = x - \frac{f(x)}{f'(x)}$.

Computing the derivative:
$$
g'(x) = 1 - \frac{(f')^2 - ff''}{(f')^2} = \frac{ff''}{(f')^2}
$$

At the root:
$$
g'(c) = \frac{f(c)f''(c)}{(f'(c))^2} = 0
$$

Since $g'(c) = 0$ and $g$ is continuous, there exists $\delta > 0$ such that $|g'(x)| < 1$ for $x \in (c-\delta, c+\delta)$. By the fixed point convergence theorem, the iteration converges.
:::

:::{prf:theorem} Quadratic Convergence
:label: thm-newton-quadratic

Under the same assumptions as [](#thm-newton-convergence), Newton's method converges with order at least $p = 2$ (quadratic convergence). Specifically:

$$
|x_{n+1} - c| \leq M |x_n - c|^2
$$

for some constant $M > 0$ depending on $f$.
:::

:::{prf:proof}
:class: dropdown

This follows from the fixed point order theorem: since $g'(c) = 0$ but generally $g''(c) \neq 0$, the convergence is quadratic.

More precisely, Taylor expansion of $g$ around $c$ gives:
$$
x_{n+1} - c = g(x_n) - g(c) = \underbrace{g'(c)}_{=0}(x_n - c) + \frac{g''(\xi)}{2}(x_n - c)^2
$$
for some $\xi$ between $x_n$ and $c$.
:::

:::{prf:remark} Newton's Method as a Local Contraction
:label: rmk-newton-contraction
:class: dropdown

The convergence proof reveals why Newton's method is "local": it's a **contraction mapping** only near the root!

Recall from the [Banach Fixed Point Theorem](fixed-point.md#theoretical-foundation-the-banach-fixed-point-theorem) that contractions converge geometrically. For Newton's iteration $g(x) = x - f(x)/f'(x)$:

- At the root: $|g'(c)| = 0 < 1$ ✓ (strongly contractive)
- Far from the root: $|g'(x)|$ may be large (not a contraction!)

This explains Newton's behavior:
- **Good initial guess:** You're in the contraction region → rapid convergence
- **Bad initial guess:** You're outside the contraction region → possible divergence

The condition $|g'(x)| < 1$ from fixed-point theory tells us exactly when Newton iterates move closer to the root. Newton's special feature is that $g'(c) = 0$, giving **quadratic** rather than just linear convergence.
:::

## Example: Computing $\sqrt{3}$

:::{prf:example} Babylonian Method
:label: ex-babylonian
:class: dropdown

For $f(x) = x^2 - 3$, Newton's method gives:

$$
x_{n+1} = x_n - \frac{x_n^2 - 3}{2x_n} = \frac{x_n^2 + 3}{2x_n} = \frac{1}{2}\left(x_n + \frac{3}{x_n}\right)
$$

This is the **Babylonian method** for square roots! The number of correct digits roughly doubles each iteration.

Starting from $x_0 = 2$:

| $n$ | $x_n$ | Error |
|-----|-------|-------|
| 0 | 2.0 | 0.27 |
| 1 | 1.75 | 0.018 |
| 2 | 1.732143 | $9 \times 10^{-5}$ |
| 3 | 1.7320508 | $2 \times 10^{-9}$ |
:::

:::{seealso}
[Fixed Point Iteration Demo](fixed-point-iteration.ipynb) — Compares different iteration functions for $\sqrt{3}$
:::

## Newton's Method in Higher Dimensions

For a vector-valued function $\mathbf{f}: \mathbb{R}^n \to \mathbb{R}^n$, Newton's method becomes:

$$
\mathbf{x}_{n+1} = \mathbf{x}_n - [D\mathbf{f}(\mathbf{x}_n)]^{-1}\mathbf{f}(\mathbf{x}_n)
$$

where $D\mathbf{f}$ is the **Jacobian matrix**.

:::{prf:definition} Jacobian Matrix
:label: def-jacobian

For $\mathbf{f}: \mathbb{R}^n \to \mathbb{R}^m$, the **Jacobian matrix** is:

$$
D\mathbf{f}(\mathbf{x}) = \begin{pmatrix}
\frac{\partial f_1}{\partial x_1} & \cdots & \frac{\partial f_1}{\partial x_n} \\
\vdots & \ddots & \vdots \\
\frac{\partial f_m}{\partial x_1} & \cdots & \frac{\partial f_m}{\partial x_n}
\end{pmatrix}
$$
:::

In 2D, for $\mathbf{f}(x,y) = \begin{pmatrix} f_1(x,y) \\ f_2(x,y) \end{pmatrix}$:

$$
\begin{pmatrix} x_{n+1} \\ y_{n+1} \end{pmatrix} = \begin{pmatrix} x_n \\ y_n \end{pmatrix} - \begin{pmatrix} \frac{\partial f_1}{\partial x} & \frac{\partial f_1}{\partial y} \\ \frac{\partial f_2}{\partial x} & \frac{\partial f_2}{\partial y} \end{pmatrix}^{-1} \begin{pmatrix} f_1(x_n, y_n) \\ f_2(x_n, y_n) \end{pmatrix}
$$

Each Newton step requires solving a linear system—this motivates the linear algebra chapters!

## Multiple Roots

:::{prf:proposition} Reduced Convergence for Multiple Roots
:label: prop-multiple-roots

When $c$ is a root of multiplicity $p \geq 2$:
$$
f(x) = (x-c)^p h(x), \quad h(c) \neq 0
$$

Newton's method only converges **linearly** with rate $C = 1 - 1/p$.
:::

:::{prf:remark} Modified Newton for Multiple Roots
:label: rmk-modified-newton
:class: dropdown

If the multiplicity $p$ is known, the modified iteration:
$$
x_{n+1} = x_n - p\frac{f(x_n)}{f'(x_n)}
$$
restores quadratic convergence.
:::

## Advantages and Disadvantages

**Advantages:**
- **Fast** — Quadratic convergence; digits double each iteration
- **Generalizes** — Extends naturally to systems in $\mathbb{R}^n$
- **Foundation** — Basis for many advanced optimization methods

**Disadvantages:**
- **Requires derivative** — Need $f'(x)$; can approximate with finite differences
- **Local only** — May diverge with bad initial guess
- **Singular points** — Fails where $f'(x) \approx 0$

**Variants:**
- **Secant method:** Approximates $f'(x_n)$ using finite differences. Order $\approx 1.618$ (the golden ratio!)
- **Quasi-Newton:** For systems, approximate the Jacobian to reduce cost
