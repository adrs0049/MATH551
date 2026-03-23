---
kernelspec:
  name: python3
  display_name: Python 3
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/fixed-point.pdf
    id: nonlinear-equations-fixed-point-pdf
downloads:
  - id: nonlinear-equations-fixed-point-pdf
    title: Download PDF
---

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

That's it. But this simplicity raises three fundamental questions:

1. **Does a fixed point even exist?** And if so, is it unique?
2. **Does the iteration converge?** Starting from $x_0$, will the sequence $\{x_n\}$ actually reach the fixed point?
3. **How fast does it converge?** Can we do better than bisection's linear rate?

The answer to all three turns on a single quantity: **$|g'(x)|$ near the fixed point**. When $|g'| < 1$, the map $g$ is a *contraction* — it pulls nearby points closer together, and the iteration spirals inward. When $|g'| \geq 1$, it pushes points apart, and the iteration diverges. The choice of $g$ controls $g'$, which is why different reformulations of the same problem can behave so differently.

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

When does a fixed point exist? When is it unique? The answers come from two conditions, each with a clear geometric meaning.

**Existence** requires that $g$ maps $[a,b]$ into itself: $g([a,b]) \subseteq [a,b]$. Geometrically, this traps the graph of $g$ inside the box $[a,b] \times [a,b]$. Since the graph enters the box with $g(a) \geq a$ and exits with $g(b) \leq b$, it *must* cross the diagonal $y = x$ somewhere — that crossing is a fixed point.

**Uniqueness** requires that $|g'(x)| < 1$ on the interval — the graph of $g$ is *flatter* than the diagonal. A curve flatter than the line it crosses can only cross it once.

```{code-cell} python
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

a, b = 0.5, 2.5
x = np.linspace(a, b, 200)

# --- Left panel: Existence (g maps [a,b] into [a,b]) ---
g1 = lambda x: 0.4*np.sin(3*x) + 1.5  # maps [a,b] into [a,b], multiple fixed points
y1 = g1(x)

# The box
ax1.plot([a, b, b, a, a], [a, a, b, b, a], 'k-', linewidth=1.5, label=r'Box $[a,b] \times [a,b]$')
ax1.fill_between([a, b], a, b, alpha=0.08, color='blue')
ax1.plot(x, x, 'k--', linewidth=1, alpha=0.6, label=r'$y = x$')
ax1.plot(x, y1, 'b-', linewidth=2.5, label=r'$y = g(x)$')

# Mark fixed points (where g(x) = x)
from scipy.optimize import brentq
# Find crossings
crossings = []
h1 = y1 - x
for i in range(len(h1)-1):
    if h1[i] * h1[i+1] < 0:
        root = brentq(lambda t: g1(t) - t, x[i], x[i+1])
        crossings.append(root)

for c in crossings:
    ax1.plot(c, c, 'ro', markersize=8, zorder=5)

ax1.set_xlim(a - 0.15, b + 0.15)
ax1.set_ylim(a - 0.15, b + 0.15)
ax1.set_xlabel('$x$')
ax1.set_ylabel('$y$')
ax1.set_title('Existence: $g([a,b]) \\subseteq [a,b]$\ngraph must cross $y = x$')
ax1.legend(fontsize=8, loc='upper left')
ax1.set_aspect('equal')

# --- Right panel: Uniqueness (|g'| < 1 means flatter than diagonal) ---
g2 = lambda x: 0.4*x + 0.9   # |g'| = 0.4 < 1, unique fixed point
y2 = g2(x)

ax2.plot([a, b, b, a, a], [a, a, b, b, a], 'k-', linewidth=1.5)
ax2.fill_between([a, b], a, b, alpha=0.08, color='blue')
ax2.plot(x, x, 'k--', linewidth=1, alpha=0.6, label=r'$y = x$ (slope 1)')
ax2.plot(x, y2, 'b-', linewidth=2.5, label=r"$y = g(x)$, $|g'| = 0.4 < 1$")

# Mark the unique fixed point
c2 = 0.9 / 0.6  # g(x) = x => 0.4x + 0.9 = x => x = 1.5
ax2.plot(c2, c2, 'ro', markersize=8, zorder=5, label='Unique fixed point')

# Show slope comparison
xm = 1.8
ax2.annotate('', xy=(xm + 0.5, g2(xm) + 0.2), xytext=(xm, g2(xm)),
            arrowprops=dict(arrowstyle='->', color='blue', lw=1.5))
ax2.annotate('', xy=(xm + 0.5, xm + 0.5), xytext=(xm, xm),
            arrowprops=dict(arrowstyle='->', color='gray', lw=1.5))
ax2.text(xm + 0.52, g2(xm) + 0.05, "slope $< 1$", fontsize=9, color='blue')
ax2.text(xm + 0.52, xm + 0.35, "slope $= 1$", fontsize=9, color='gray')

ax2.set_xlim(a - 0.15, b + 0.15)
ax2.set_ylim(a - 0.15, b + 0.15)
ax2.set_xlabel('$x$')
ax2.set_ylabel('$y$')
ax2.set_title('Uniqueness: $|g\'| < 1$\nflatter than diagonal $\\Rightarrow$ one crossing')
ax2.legend(fontsize=8, loc='upper left')
ax2.set_aspect('equal')

plt.tight_layout()
plt.show()
```

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

:::{prf:remark} Understanding Geometric Convergence
:label: rmk-geometric-convergence

The bound $|x_n - c| \leq \rho^n |x_0 - c|$ is called **geometric** (or **linear**) convergence. At each step, the error is multiplied by at most $\rho$:

$$
|e_{n+1}| \leq \rho\, |e_n|
$$

The value of $\rho$ — the **contraction factor** — controls the speed:
- $\rho$ close to $0$: fast convergence (each step eliminates most of the error)
- $\rho$ close to $1$: painfully slow (each step barely improves the approximation)
- $\rho = 1$: no guarantee of convergence at all

Taking logarithms of the error bound gives $\log|e_n| \leq n \log \rho + \log|e_0|$, so **geometric convergence appears as a straight line on a log-scale plot**, with slope $\log \rho$. Smaller $\rho$ means a steeper downward slope.
:::

```{code-cell} python
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(8, 4.5))

n = np.arange(0, 30)
e0 = 1.0

for rho, color, style in [(0.9, '#d62728', '-'),
                           (0.7, '#ff7f0e', '-'),
                           (0.5, '#2ca02c', '-'),
                           (0.2, '#1f77b4', '-')]:
    errors = e0 * rho**n
    ax.semilogy(n, errors, style, color=color, linewidth=2,
                marker='o', markersize=3,
                label=rf'$\rho = {rho}$  (slope $\log_{{10}}\rho = {np.log10(rho):.2f}$)')

ax.set_xlabel('Iteration $n$', fontsize=12)
ax.set_ylabel('Error bound $\\rho^n |e_0|$', fontsize=12)
ax.set_title('Geometric convergence: smaller $\\rho$ = steeper descent on log scale')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3, which='both')
ax.set_xlim(0, 29)
ax.set_ylim(1e-15, 2)

plt.tight_layout()
plt.show()
```

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

:::{prf:remark} What Higher-Order Convergence Really Means
:label: rmk-higher-order-convergence

The error relation $|e_{n+1}| \approx C|e_n|^p$ has dramatic consequences as $p$ increases.

**Linear** ($p = 1$): Each step multiplies the error by $\rho$. To gain one digit of accuracy, you always need a fixed number of iterations. Progress is steady but slow.

**Quadratic** ($p = 2$): Each step *squares* the error. If $e_n \approx 10^{-3}$, then $e_{n+1} \approx 10^{-6}$, then $e_{n+2} \approx 10^{-12}$. The number of correct digits **doubles** at every step. This is why Newton's method feels like it "suddenly" converges — the first few iterations seem slow, then accuracy explodes.

**Cubic** ($p = 3$): Each step *cubes* the error. The number of correct digits **triples** each step — even faster, but rarely worth the extra cost per iteration in practice.
:::

```{code-cell} python
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(8, 5))

# Simulate error sequences with |e_{n+1}| = C * |e_n|^p
e0 = 0.5
n_max = 12

# Linear: e_{n+1} = rho * e_n, rho = 0.5
rho = 0.5
e_linear = [e0]
for _ in range(n_max - 1):
    e_linear.append(rho * e_linear[-1])

# Quadratic: e_{n+1} = C * e_n^2, C = 1 (so e_{n+1} = e_n^2)
e_quad = [e0]
for _ in range(n_max - 1):
    e_next = e_quad[-1]**2
    if e_next < 1e-16:
        e_quad.append(1e-16)
    else:
        e_quad.append(e_next)

# Cubic: e_{n+1} = C * e_n^3, C = 1
e_cubic = [e0]
for _ in range(n_max - 1):
    e_next = e_cubic[-1]**3
    if e_next < 1e-16:
        e_cubic.append(1e-16)
    else:
        e_cubic.append(e_next)

n = np.arange(n_max)

ax.semilogy(n, e_linear, 'o-', color='#d62728', linewidth=2, markersize=6,
            label=r'Linear ($p=1$, $\rho=0.5$): $e_{n+1} = 0.5\, e_n$')
ax.semilogy(n, e_quad, 's-', color='#1f77b4', linewidth=2, markersize=6,
            label=r'Quadratic ($p=2$): $e_{n+1} = e_n^2$')
ax.semilogy(n, e_cubic, 'D-', color='#2ca02c', linewidth=2, markersize=6,
            label=r'Cubic ($p=3$): $e_{n+1} = e_n^3$')

ax.axhline(y=1e-15, color='gray', linestyle=':', linewidth=1, label='Machine epsilon')

# Annotate the quadratic "doubling digits" effect
ax.annotate('1 digit', xy=(1, e_quad[1]), xytext=(1.6, 1e-1),
            fontsize=9, color='#1f77b4',
            arrowprops=dict(arrowstyle='->', color='#1f77b4', lw=1.2))
ax.annotate('2 digits', xy=(2, e_quad[2]), xytext=(2.6, 3e-3),
            fontsize=9, color='#1f77b4',
            arrowprops=dict(arrowstyle='->', color='#1f77b4', lw=1.2))
ax.annotate('4 digits', xy=(3, e_quad[3]), xytext=(3.6, 3e-5),
            fontsize=9, color='#1f77b4',
            arrowprops=dict(arrowstyle='->', color='#1f77b4', lw=1.2))
ax.annotate('8 digits', xy=(4, e_quad[4]), xytext=(4.6, 3e-7),
            fontsize=9, color='#1f77b4',
            arrowprops=dict(arrowstyle='->', color='#1f77b4', lw=1.2))

ax.set_xlabel('Iteration $n$', fontsize=12)
ax.set_ylabel('Error $|e_n|$', fontsize=12)
ax.set_title('Linear vs quadratic vs cubic convergence')
ax.legend(fontsize=9, loc='lower left')
ax.grid(True, alpha=0.3, which='both')
ax.set_xlim(-0.3, n_max - 0.5)
ax.set_ylim(1e-16, 2)

plt.tight_layout()
plt.show()
```

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
:label: rmk-banach-fpt-applications
:class: dropdown

The Banach FPT is not just about scalar equations. The same principle governs:

- **Newton's method for systems**: The iteration is a contraction near the solution
- **Picard iteration for ODEs**: Proves existence and uniqueness for $y' = f(t,y)$
- **Iterative linear solvers**: Jacobi and Gauss-Seidel converge when the iteration matrix is a contraction

Whenever you see an iterative method that "works," there's often a contraction hiding underneath.
:::

## Summary: How to Choose $g(x)$

Everything in this section comes down to one question: **how do you pick $g$?**

Given $f(x) = 0$, there are infinitely many ways to rewrite it as $x = g(x)$. The right choice is governed by $|g'(c)|$ at the fixed point:

| $|g'(c)|$ | Behavior | Convergence |
|-----------|----------|-------------|
| $> 1$ | Iteration diverges | None |
| $= 1$ | Borderline — may or may not converge | Uncertain |
| $0 < \rho < 1$ | Linear convergence with rate $\rho$ | Slow if $\rho \approx 1$, fast if $\rho \approx 0$ |
| $= 0$ | At least quadratic convergence | Fast — digits of accuracy double each step |

**The design principle:** make $|g'(c)|$ as small as possible. The optimal choice is $g'(c) = 0$, which yields quadratic convergence. This is exactly what Newton's method achieves by setting $g(x) = x - f(x)/f'(x)$ — as we'll see in the next section.
