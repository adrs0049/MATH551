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

**The choice of $g$ matters enormously** for convergence, as we'll see.

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

The answer to all three turns on a single quantity: **$|g'(x)|$ near the fixed point**. When $|g'| < 1$, the map $g$ is a *contraction*: it pulls nearby points closer together, and the iteration spirals inward. When $|g'| \geq 1$, it pushes points apart, and the iteration diverges. The choice of $g$ controls $g'$, which is why different reformulations of the same problem can behave so differently.

## Why the Choice of $g$ Matters

Consider finding the root of $f(x) = x^2 - 3 = 0$ (i.e., finding $\sqrt{3}$).

Here are three valid reformulations:

- **G1:** From $x^2 = 3$, we get $g_1(x) = 3/x$
- **G2:** Add $x$ to both sides: $g_2(x) = x - (x^2 - 3)/2$
- **G3:** Divide by $2x$: $g_3(x) = (x^2 + 3)/(2x)$

All three have $\sqrt{3}$ as a fixed point. But their behavior is dramatically different:
- **G1** cycles forever, never converging
- **G2** converges slowly (linearly)
- **G3** converges rapidly (quadratically, this is Newton's method!)

:::{seealso}
[Fixed Point Iteration Demo](fixed-point-iteration.ipynb): Compare G1, G2, G3 convergence behavior
:::

## Existence, Uniqueness and Convergence

There are three separate questions to answer, each requiring its own condition:

**1. Existence:**  *Does a fixed point exist at all?* This only requires that
$g$ maps $[a,b]$ into itself: $g([a,b]) \subseteq [a,b]$. Geometrically, this
traps the graph of $g$ inside the box $[a,b] \times [a,b]$. Since $g(a) \geq a$
and $g(b) \leq b$, the graph *must* cross the diagonal $y = x$. That crossing
is a fixed point (left panel). But there could be many crossings.

**2. Uniqueness:**  *Is there only one fixed point?* This requires $|g'(x)| \leq
\rho < 1$ on the interval. The map is then a *contraction*: it shrinks
distances. If two fixed points $c_1, c_2$ existed, the MVT gives $|c_1 - c_2|
= |g(c_1) - g(c_2)| = |g'(\xi)||c_1 - c_2| \leq \rho|c_1 - c_2|$,
a contradiction since $\rho < 1$. Geometrically, $|g'| < 1$ means the graph is
everywhere less steep than the diagonal, so it can only cross once (center
panel). When $|g'| > 1$, the graph can be steeper than the diagonal, allowing
multiple crossings (right panel).

**3. Convergence:**  *Does the iteration $x_{n+1} = g(x_n)$ actually reach the
fixed point?* This is a question about the *dynamics* of iterating $g$, not just
its graph. The same contraction condition $|g'| \leq \rho < 1$ guarantees
convergence: each step shrinks the distance to the fixed point, so the cobweb
spirals inward (center panel). When $|g'| > 1$, the iteration *amplifies* errors.
The cobweb diverges away from the fixed point even if you start nearby (right
panel).

```{code-cell} python
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
ax1, ax2, ax3 = axes

a, b = 0.5, 2.5
x = np.linspace(a, b, 200)

from scipy.optimize import brentq

def plot_box(ax):
    ax.plot([a, b, b, a, a], [a, a, b, b, a], 'k-', linewidth=1.5)
    ax.fill_between([a, b], a, b, alpha=0.08, color='blue')
    ax.plot(x, x, 'k--', linewidth=1, alpha=0.6)

def find_and_mark_crossings(ax, g):
    """Find and mark fixed points; return list of roots."""
    y = g(x)
    h = y - x
    roots = []
    for i in range(len(h)-1):
        if h[i] * h[i+1] < 0:
            root = brentq(lambda t: g(t) - t, x[i], x[i+1])
            roots.append(root)
            ax.plot(root, root, 'ro', markersize=8, zorder=5)
    return roots

def setup_ax(ax, title):
    ax.set_xlim(a - 0.15, b + 0.15)
    ax.set_ylim(a - 0.15, b + 0.15)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    ax.set_title(title)
    ax.set_aspect('equal')

# --- Panel 1: Existence (g maps [a,b] into [a,b]) ---
g1 = lambda x: 0.4*np.sin(3*x) + 1.5
plot_box(ax1)
ax1.plot(x, g1(x), 'b-', linewidth=2.5)
find_and_mark_crossings(ax1, g1)
setup_ax(ax1, 'Existence: $g([a,b]) \\subseteq [a,b]$\ngraph must cross $y = x$')

# --- Panel 2: Contraction (|g'| < 1) => unique + converges ---
g2 = lambda x: 0.4*x + 0.9   # |g'| = 0.4 < 1
plot_box(ax2)
ax2.plot(x, g2(x), 'b-', linewidth=2.5)
find_and_mark_crossings(ax2, g2)

# Cobweb showing contraction
x_cob = 2.2
ax2.plot(x_cob, x_cob, 'k*', markersize=12, zorder=6)
for _ in range(6):
    x_next = g2(x_cob)
    ax2.plot([x_cob, x_cob], [x_cob, x_next], 'g-', alpha=0.5, linewidth=1)
    ax2.plot([x_cob, x_next], [x_next, x_next], 'g-', alpha=0.5, linewidth=1)
    x_cob = x_next
ax2.plot(x_cob, x_cob, 'rs', markersize=8, zorder=6)

setup_ax(ax2, "Contraction: $|g'| < 1$\ncobweb spirals in, unique fixed point")

# --- Panel 3: Nonlinear with |g'| > 1 in middle => multiple fixed points ---
# Steep S-curve: maps [a,b] into [a,b] but |g'| > 1 in the middle
g3 = lambda x: 1.5 + 0.95*np.tanh(2.5*(x - 1.5))
plot_box(ax3)
ax3.plot(x, g3(x), 'r-', linewidth=2.5)
roots = find_and_mark_crossings(ax3, g3)

# Cobweb starting near the central fixed point, showing divergence away from it
# The middle root is the unstable one
middle_root = sorted(roots)[1]
x_cob = middle_root + 0.05
ax3.plot(x_cob, x_cob, 'k*', markersize=12, zorder=6)
for _ in range(12):
    x_next = g3(x_cob)
    if x_next < a - 0.3 or x_next > b + 0.3:
        break
    ax3.plot([x_cob, x_cob], [x_cob, x_next], 'g-', alpha=0.5, linewidth=1)
    ax3.plot([x_cob, x_next], [x_next, x_next], 'g-', alpha=0.5, linewidth=1)
    x_cob = x_next
ax3.plot(x_cob, x_cob, 'rs', markersize=8, zorder=6)

setup_ax(ax3, "Not a contraction: $|g'| > 1$\ncobweb diverges from central fixed point")

plt.tight_layout()
plt.show()
```

:::{figure} #
:label: fig-fixed-point-geometry

**Geometric intuition for fixed point iteration.** In the cobweb diagrams (center and right), the black star ($\star$) marks the starting point $x_0$ and the red square marks the final iterate. *Left:* If $g$ maps $[a,b]$ into itself, its graph stays inside the blue box and must cross the diagonal $y = x$ at least once, guaranteeing existence of a fixed point. *Center:* When $|g'| < 1$ everywhere on $[a,b]$, the graph is less steep than the diagonal, so there is exactly one crossing. The cobweb diagram shows the iteration spiraling inward from the start ($\star$) to the unique fixed point. *Right:* When $g([a,b]) \subseteq [a,b]$ but $|g'| > 1$ near the center, the graph is steep enough to cross the diagonal multiple times, creating three fixed points. The cobweb starting near the central (unstable) fixed point ($\star$) diverges away from it, eventually settling at one of the stable outer fixed points where $|g'| < 1$.
:::

:::{prf:theorem} Fixed Point Existence, Uniqueness, and Convergence
:label: thm-fixed-point-existence

Given $g: [a, b] \to \mathbb{R}$ continuous:

**Existence:** If $g$ maps $[a,b]$ into itself (i.e., $g([a,b]) \subseteq [a,b]$), then $g$ has at least one fixed point in $[a,b]$.

**Uniqueness:** If additionally $|g'(x)| \leq \rho < 1$ for all $x \in [a,b]$, then the fixed point is unique.

**Convergence:** Under the same conditions, for any $x_0 \in [a,b]$, the sequence $x_{n+1} = g(x_n)$ converges to the unique fixed point $c$, and the convergence is geometric:
$$
|x_n - c| \leq \rho^n |x_0 - c|
$$
:::

:::{prf:proof} Existence
:class: dropdown

Define $h(x) = x - g(x)$.

If $g(a) = a$ or $g(b) = b$, we're done. We've found a fixed point.

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

This implies $(1-\rho)|c_1 - c_2| \leq 0$. But $\rho < 1$ and $c_1 \neq c_2$, so $(1-\rho)|c_1 - c_2| > 0$. This is a contradiction.
:::

:::{prf:proof} Convergence
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
:class: dropdown

The bound $|x_n - c| \leq \rho^n |x_0 - c|$ is called **geometric** (or **linear**) convergence. At each step, the error is multiplied by at most $\rho$:

$$
|e_{n+1}| \leq \rho\, |e_n|
$$

The value of $\rho$, called the **contraction factor**, controls the speed:
- $\rho$ close to $0$: fast convergence (each step eliminates most of the error).
- $\rho$ close to $1$: painfully slow (each step barely improves the approximation).
- $\rho = 1$: no convergence guarantee.

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
:class: dropdown

Even if $|g'(x)| < 1$ only *at* the fixed point (not on the whole interval), the
iteration still converges, provided we start close enough. This is called
**local convergence**.

Specifically: if $g \in \mathcal{C}^1$, $g(c) = c$, and $|g'(c)| < 1$, then
there exists $\delta > 0$ such that the iteration converges for any $x_0$ with
$|x_0 - c| < \delta$.
:::

## The Derivative at the Fixed Point

The key insight: **$|g'(c)|$ determines everything**.

For our three reformulations of $x^2 - 3 = 0$:

- $g_1(x) = 3/x$: We have $g_1'(x) = -3/x^2$, so $|g_1'(\sqrt{3})| = 1$. Right on the boundary, so no convergence is guaranteed (and indeed, it fails).

- $g_2(x) = x - (x^2-3)/2$: We have $g_2'(x) = 1 - x$, so $|g_2'(\sqrt{3})| = |1 - \sqrt{3}| \approx 0.73$. Linear convergence with rate $\rho \approx 0.73$.

- $g_3(x) = (x^2+3)/(2x)$: We have $g_3'(x) = 1/2 - 3/(2x^2)$, so $g_3'(\sqrt{3}) = 0$. The derivative vanishes, which signals faster-than-linear convergence.

## Order of Convergence

We know that fixed point iteration converges when $|g'(c)| < 1$, but *how fast*?
The contraction factor $\rho = |g'(c)|$ controls linear convergence, but what
happens when $g'(c) = 0$, as with $g_3$? We need a finer notion of convergence
speed.

:::{prf:definition} Order of Convergence
:label: def-convergence-order

A sequence $\{x_n\}$ converging to $\ell$ has **order** $p$ if:

$$
\lim_{n\to\infty} \frac{|x_{n+1} - \ell|}{|x_n - \ell|^p} = C
$$

for some constant $C > 0$.

- $p = 1$: **Linear** convergence (error shrinks by a constant factor each step)
- $p = 2$: **Quadratic** convergence (digits of accuracy double each step)
- $p = 3$: **Cubic** convergence (digits of accuracy triple each step)
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

This explains why $g_3$ converges so fast. We already know $g_3'(\sqrt{3}) = 0$.
Computing the second derivative:

$$
g_3'(x) = \frac{1}{2} - \frac{3}{2x^2}, \qquad g_3''(x) = \frac{3}{x^3}
$$

so $g_3''(\sqrt{3}) = 3/(3\sqrt{3}) = 1/\sqrt{3} \neq 0$. Since $g_3'(c) = 0$
and $g_3''(c) \neq 0$, the theorem gives $p = 2$: **quadratic convergence** with
asymptotic error constant $C = |g_3''(\sqrt{3})|/2! = 1/(2\sqrt{3}) \approx 0.289$.

:::{prf:remark} What Higher-Order Convergence Means
:label: rmk-higher-order-convergence
:class: dropdown

The error relation $|e_{n+1}| \approx C|e_n|^p$ has dramatic consequences as $p$ increases.

**Linear** ($p = 1$): Each step multiplies the error by $\rho$. To gain one
digit of accuracy, you always need a fixed number of iterations. Progress is
steady but slow.

**Quadratic** ($p = 2$): Each step *squares* the error. If $e_n \approx
10^{-3}$, then $e_{n+1} \approx 10^{-6}$, then $e_{n+2} \approx 10^{-12}$. The
number of correct digits **doubles** at every step. This is characteristic of
Newton's method: the initial iterations show modest improvement, but once the
error becomes small the doubling effect produces rapid convergence to machine
precision in just a few steps.

**Cubic** ($p = 3$): Each step *cubes* the error. The number of correct digits **triples** each step, even faster.
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

**The design principle:** make $|g'(c)|$ as small as possible. The optimal
choice is $g'(c) = 0$, which yields quadratic convergence. This is exactly what
Newton's method achieves by setting $g(x) = x - f(x)/f'(x)$.

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


