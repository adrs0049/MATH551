---
kernelspec:
  name: python3
  display_name: Python 3
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/practical-root-finding.pdf
    id: nonlinear-equations-practical-root-finding-pdf
downloads:
  - id: nonlinear-equations-practical-root-finding-pdf
    title: Download PDF
---

# Practical Root Finding

:::{tip} Big Idea
Real-world problems often have multiple roots, and we need to find all of them
reliably. This requires combining **bracketing** (to isolate roots) with
**fast solvers** (to compute them accurately). Production solvers like MATLAB's
`fzero` use hybrid methods that get the best of both worlds.
:::

## The Challenge: Multiple Roots

So far we have studied methods that converge to *one* root:

- **Bisection** is safe but slow, and requires a bracket $[a, b]$ with
  $f(a)f(b) < 0$.
- **Newton's method** is fast but local: it converges to whichever root is
  nearest (roughly), and may miss others entirely.

In practice, a function can have many roots, and we need a systematic way to
find them all.

:::{prf:example} Multiple Roots
:label: ex-multiple-roots-practical

Consider $f(x) = x^3 - 6x^2 + 11x - 6 = (x-1)(x-2)(x-3)$. There are three
roots at $x = 1, 2, 3$. Starting Newton's method from $x_0 = 0$ converges to
$x = 1$. Starting from $x_0 = 2.5$ converges to $x = 3$ (not the nearest root
$x = 2$). There is no single starting point that finds all three.
:::

## Step 1: Bracketing All Roots

The first step is to **isolate** each root in its own interval. The simplest
approach: evaluate $f$ on a grid and look for sign changes.

:::{prf:algorithm} Sign-Change Scanning
:label: alg-sign-change-scanning

**Input:** Function $f$, interval $[a, b]$, number of sample points $N$

**Output:** List of brackets $[a_i, b_i]$ with $f(a_i)f(b_i) < 0$

1. Set $h = (b - a)/N$
2. **for** $i = 0, 1, \ldots, N-1$:
3. $\qquad$ $x_i = a + ih$, $x_{i+1} = a + (i+1)h$
4. $\qquad$ **if** $f(x_i) \cdot f(x_{i+1}) < 0$: record bracket $[x_i, x_{i+1}]$
:::

By the Intermediate Value Theorem, each bracket contains at least one root.

:::{prf:remark} Limitations of Sign-Change Scanning
:label: rmk-sign-change-limitations

This approach can miss roots in two situations:

1. **Even-multiplicity roots** (e.g., $f(x) = x^2$): the function touches zero
   but does not change sign, so no bracket is detected.
2. **Closely spaced roots**: if two roots are closer than the grid spacing $h$,
   they may fall in the same interval. The sign change detects only one bracket,
   and bisection/Newton applied to that bracket may find only one of the two
   roots.

A finer grid reduces the chance of missing roots but increases the cost.
:::

```{code-cell} python
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt

f = lambda x: x**3 - 6*x**2 + 11*x - 6

x = np.linspace(-0.5, 4, 500)
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(x, f(x), 'b-', linewidth=2)
ax.axhline(0, color='k', linewidth=0.8)

# Sign-change scanning
N = 20
x_grid = np.linspace(-0.5, 4, N + 1)
f_grid = f(x_grid)

brackets = []
for i in range(N):
    if f_grid[i] * f_grid[i + 1] < 0:
        brackets.append((x_grid[i], x_grid[i + 1]))
        ax.axvspan(x_grid[i], x_grid[i + 1], alpha=0.2, color='green')

# Mark roots
roots = [1, 2, 3]
for r in roots:
    ax.plot(r, 0, 'ro', markersize=8, zorder=5)

ax.set_xlabel('$x$')
ax.set_ylabel('$f(x)$')
ax.set_title(f'Sign-change scanning: {len(brackets)} brackets found (green)')
ax.grid(True, alpha=0.3)
ax.set_ylim(-3, 3)
plt.tight_layout()
plt.show()
```

## Step 2: Solve Each Bracket

Once we have brackets, apply a fast solver to each one independently. This is
the **bracket-then-solve** strategy.

```{code-cell} python
:tags: [hide-input]

from scipy.optimize import brentq

print("Brackets found and roots computed:")
for a_i, b_i in brackets:
    root = brentq(f, a_i, b_i)
    print(f"  [{a_i:.3f}, {b_i:.3f}] -> root = {root:.12f}")
```

## Hybrid Methods: The Best of Both Worlds

The most effective solvers combine bracketing (safety) with Newton-like methods
(speed). The idea: maintain a bracket $[a, b]$ that always contains the root,
but use faster steps when possible.

### Brent's Method

MATLAB's `fzero` and SciPy's `brentq` implement **Brent's method** (1973),
which combines three strategies:

1. **Bisection**: always safe, always shrinks the bracket by half.
2. **Secant method**: uses the two most recent points to extrapolate (no
   derivative needed, superlinear convergence).
3. **Inverse quadratic interpolation**: fits a parabola through the three most
   recent points for even faster convergence.

At each step, Brent's method tries the fastest option first (inverse quadratic
interpolation), falls back to secant if that fails, and uses bisection as the
last resort. The bracket is always maintained, so convergence is guaranteed.

:::{prf:remark} Why Not Newton Inside a Bracket?
:label: rmk-why-not-newton-bracket

One could also maintain a bracket and use Newton steps when they stay inside the
bracket, falling back to bisection otherwise. This is called the
**Newton-bisection hybrid** and works well in practice. Brent's method avoids
Newton because it does not require the derivative $f'$, which makes it more
broadly applicable.
:::

### What MATLAB's `fzero` Does

MATLAB's `fzero(f, x0)` works in two phases:

1. **Phase 1 (bracket finding):** starting from the initial guess `x0`, it
   searches outward (expanding intervals) until it finds a sign change. This
   gives a bracket $[a, b]$.
2. **Phase 2 (bracket solving):** apply Brent's method to the bracket.

If called as `fzero(f, [a, b])` with a bracket directly, it skips Phase 1.

SciPy's `brentq(f, a, b)` requires the bracket up front and goes directly to
Brent's method.

## Finding All Roots Systematically

Combining everything, a robust strategy for finding all roots of $f$ on $[a, b]$:

:::{prf:algorithm} Find All Roots
:label: alg-find-all-roots

**Input:** Function $f$, interval $[a, b]$, grid size $N$, tolerance $\varepsilon$

**Output:** List of all detected roots

1. **Scan:** evaluate $f$ on a grid of $N+1$ points in $[a, b]$.
2. **Bracket:** identify all intervals $[x_i, x_{i+1}]$ where $f$ changes sign.
3. **Solve:** apply Brent's method (or Newton-bisection) to each bracket.
4. **Deduplicate:** remove any roots that are within $\varepsilon$ of each other.
:::

```{code-cell} python
:tags: [hide-input]

def find_all_roots(f, a, b, N=1000, tol=1e-12):
    """Find all sign-change roots of f on [a, b]."""
    x_grid = np.linspace(a, b, N + 1)
    f_grid = f(x_grid)

    roots = []
    for i in range(N):
        if f_grid[i] * f_grid[i + 1] < 0:
            root = brentq(f, x_grid[i], x_grid[i + 1], xtol=tol)
            # Deduplicate
            if not roots or abs(root - roots[-1]) > 10 * tol:
                roots.append(root)
    return roots

# Example: find all roots of a more complex function
g = lambda x: np.sin(5*x) * np.exp(-x/3) - 0.1

roots_found = find_all_roots(g, 0, 10)

fig, ax = plt.subplots(figsize=(10, 4))
x = np.linspace(0, 10, 1000)
ax.plot(x, g(x), 'b-', linewidth=2)
ax.axhline(0, color='k', linewidth=0.8)
for r in roots_found:
    ax.plot(r, 0, 'ro', markersize=8, zorder=5)
ax.set_xlabel('$x$')
ax.set_ylabel('$f(x)$')
ax.set_title(f'Found {len(roots_found)} roots of $\\sin(5x) e^{{-x/3}} - 0.1$')
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

print("Roots found:")
for i, r in enumerate(roots_found):
    print(f"  x_{i} = {r:.10f},  f(x_{i}) = {g(r):.2e}")
```

## Deflation: An Alternative Approach

Instead of bracketing all roots up front, **deflation** finds roots one at a
time: after finding a root $r_1$, divide it out and solve the deflated problem.

For a polynomial $f(x)$ with a known root $r_1$, define:

$$
f_1(x) = \frac{f(x)}{x - r_1}
$$

The roots of $f_1$ are the remaining roots of $f$. Apply Newton's method (or
any solver) to $f_1$ to find the next root, then deflate again.

:::{prf:remark} Deflation Pitfalls
:label: rmk-deflation-pitfalls

Deflation works well for polynomials but has issues for general functions:

1. **Numerical instability**: dividing by $(x - r_1)$ amplifies errors near
   $r_1$. If $r_1$ is only approximate, the deflated function has a pole near
   (but not at) the true root.
2. **Error accumulation**: each deflation step introduces errors that propagate
   to later roots. Finding roots in order of increasing magnitude helps.
3. **Non-polynomial functions**: for transcendental functions like $\sin(x)$,
   deflation is not straightforward since there are infinitely many roots.

For these reasons, the bracket-then-solve approach is usually preferred for
general functions, while deflation is mainly used for polynomials.
:::

