---
kernelspec:
  name: python3
  display_name: Python 3
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/condition-stability.pdf
    id: qr-least-squares-condition-stability-pdf
downloads:
  - id: qr-least-squares-condition-stability-pdf
    title: Download PDF
---

# Stability in Linear Systems

:::{tip} Big Idea
A numerical algorithm cannot avoid the sensitivity inherent in the problem it
solves. The best it can do is avoid making things worse. **Backward stability**
is the formal statement of this goal: the computed answer is the exact solution
of a slightly perturbed problem. Combined with the condition number, this gives
the tightest error bound we can hope for.
:::

## Forward and Backward Error

(forward-backward-error)=

When an algorithm computes an approximate solution $\hat{x}$ to a problem
$f(x) = y$, there are two ways to measure how good $\hat{x}$ is.

:::{prf:definition} Forward and Backward Error
:label: def-forward-backward-error-linalg

For a problem $f: X \to Y$ with true solution $y = f(x)$ and computed solution
$\hat{y}$:

**Forward error:** how far is $\hat{y}$ from the true answer?

$$
\text{forward error} = \|\hat{y} - f(x)\|
$$

**Backward error:** what is the smallest perturbation $\delta x$ such that
$\hat{y}$ is the *exact* answer to the perturbed problem?

$$
\text{backward error} = \min \{ \|\delta x\| : f(x + \delta x) = \hat{y} \}
$$
:::

**For linear systems,** the problem is $f(A, \mathbf{b}) = A^{-1}\mathbf{b}$.
An algorithm produces $\hat{\mathbf{x}}$. Then:

- **Forward error:** $\|\hat{\mathbf{x}} - \mathbf{x}\|$, how far the computed solution is from the true one.
- **Backward error:** the smallest $(\delta A, \delta\mathbf{b})$ such that $(A + \delta A)\hat{\mathbf{x}} = \mathbf{b} + \delta\mathbf{b}$. The computed solution *exactly* solves a nearby system.

Forward error is what we care about. Backward error is what we can control.

(golden-rule)=

:::{prf:theorem} The Golden Rule
:label: thm-golden-rule

$$
\text{relative forward error} \lesssim \kappa \times \text{relative backward error}
$$

where $\kappa$ is the condition number of the problem.
:::

This separates concerns cleanly:

- **The problem** determines $\kappa$ (how sensitive the answer is to perturbations).
- **The algorithm** determines the backward error (how much perturbation it introduces).

Neither can compensate for the other. A stable algorithm applied to an
ill-conditioned problem still gives a poor answer. A brilliant problem
formulation is wasted on an unstable algorithm.

## Sensitivity of Linear Systems

How sensitive is the solution $\mathbf{x}$ to perturbations in $A$ and $\mathbf{b}$?

A linear system $A\mathbf{x} = \mathbf{b}$ has **two inputs**: the matrix $A$
and the vector $\mathbf{b}$. Both are subject to errors:

- **$\mathbf{b}$ comes from measurements** and always has some noise
- **$A$ comes from a model** and its coefficients may be uncertain or stored with roundoff error

So we must understand how errors in *both* $A$ and $\mathbf{b}$ propagate to errors in $\mathbf{x}$.

:::::{prf:theorem} Sensitivity of Linear Systems
:label: thm-sensitivity-linear-systems

For the linear system $A\mathbf{x} = \mathbf{b}$:

$$
\frac{\|\delta\mathbf{x}\|}{\|\mathbf{x}\|} \lesssim \kappa(A) \left(\frac{\|\delta A\|}{\|A\|} + \frac{\|\delta\mathbf{b}\|}{\|\mathbf{b}\|}\right)
$$

The quantity $\kappa(A) = \|A\| \|A^{-1}\|$ is the **amplification factor** from relative input perturbation to relative output error.

::::{prf:proof}
:class: dropdown

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

## Numerically Singular Matrices

Recall that the [condition number](linear-algebra-fundamentals.md#def-condition-number)
$\kappa(A) = \|A\|\|A^{-1}\|$ measures the ratio of maximum to minimum stretch
of a unit vector. The sensitivity theorem tells us what this means in practice.

**Rule of thumb:** Expect to lose $\log_{10}\kappa(A)$ digits of accuracy.

| Condition Number | Digits Lost |
|-----------------|-------------|
| $\kappa \approx 10^k$ | ~$k$ digits |
| $\kappa \gtrsim 1/\varepsilon_{\text{mach}} \approx 10^{16}$ | All digits |

See the [Condition Numbers and Lost Digits notebook](../notebooks/condition-number-digits.ipynb) for concrete demonstrations with diagonal, Hilbert, and Vandermonde matrices.

But what does it *mean* geometrically for a matrix to be ill-conditioned? How
close is it to being singular? The answer involves the geometry of the set of
singular matrices in the operator norm.

### The singular set is closed (optional)

::::{prf:lemma} $\Sigma$ is closed
:label: lem-singular-set-closed

The set of singular matrices
$\Sigma = \{M \in \mathbb{R}^{n\times n} : \det(M) = 0\}$ is closed in the
operator norm.

:::{prf:proof}
:class: dropdown

The determinant $\det : \mathbb{R}^{n \times n} \to \mathbb{R}$ is a
polynomial in the $n^2$ entries, so it is continuous. The singleton $\{0\}$ is
closed in $\mathbb{R}$, so $\Sigma = \det^{-1}(\{0\})$ is closed as the
preimage of a closed set under a continuous map.
:::
::::

Because $\Sigma$ is closed, every invertible matrix $A$ sits at some strictly
positive distance from $\Sigma$.

How far is an invertible matrix from $\Sigma$? The following theorem gives the
exact answer.

:::::::{prf:theorem} Distance to Singularity
:label: thm-distance-singularity

For any invertible $A \in \mathbb{R}^{n \times n}$ with the 2-norm,

$$
\operatorname{dist}_2(A, \Sigma)
\;=\; \min_{\substack{E \in \mathbb{R}^{n\times n} \\ A + E \text{ singular}}} \|E\|_2
\;=\; \sigma_{\min}(A)
\;=\; \frac{\|A\|_2}{\kappa_2(A)}.
$$

::::::{prf:proof}
:class: dropdown

Recall that $\sigma_{\min}(A) = \min_{\|\mathbf{x}\|=1} \|A\mathbf{x}\|$.

**Lower bound.** Let $E$ be any perturbation such that $A + E$ is singular.
Then there exists a unit vector $\mathbf{z}$ with $(A + E)\mathbf{z} = \mathbf{0}$,
so $E\mathbf{z} = -A\mathbf{z}$. Therefore

$$
\|E\|_2 \geq \|E\mathbf{z}\|_2 = \|A\mathbf{z}\|_2 \geq \sigma_{\min}(A).
$$

**Upper bound.** Let $\mathbf{v}$ be a unit vector achieving
$\sigma_{\min}(A) = \|A\mathbf{v}\|_2$. Set $E^* = -A\mathbf{v}\,\mathbf{v}^T$.
Then $\|E^*\|_2 = \|A\mathbf{v}\|_2 = \sigma_{\min}(A)$, and

$$
(A + E^*)\mathbf{v} = A\mathbf{v} - A\mathbf{v}(\mathbf{v}^T\mathbf{v}) = \mathbf{0},
$$

so $A + E^*$ is singular.

Combining: $\operatorname{dist}(A, \Sigma) = \sigma_{\min}(A)$.
The identity $\sigma_{\min} = \|A\|_2 / \kappa_2(A)$ follows from
$\kappa_2(A) = \sigma_{\max}/\sigma_{\min}$.
::::::
:::::::

### Demmel's insight

:::::{prf:example} Numerically singular matrices
:label: ex-demmel-numerically-singular

A matrix with $\kappa(A) \gtrsim 1/\varepsilon_{\text{mach}}$ is
**numerically indistinguishable from a singular matrix**.

By [](#thm-distance-singularity), the relative distance from $A$ to
the nearest singular matrix is

$$
\frac{\operatorname{dist}(A, \Sigma)}{\|A\|}
= \frac{\sigma_{\min}}{\sigma_{\max}}
= \frac{1}{\kappa(A)}.
$$

Now consider floating-point storage. The computer does not see $A$. It sees
$\tilde{A} = A + E$ where $\|E\|/\|A\| \sim \varepsilon_{\text{mach}}$. The
computer's view of $A$ is a matrix somewhere inside a ball of radius
$\varepsilon_{\text{mach}}\|A\|$ centered at $A$.

**If $\kappa(A) \gtrsim 1/\varepsilon_{\text{mach}}$:**

- The distance to the nearest singular matrix is
  $\sigma_{\min} = \|A\|/\kappa(A) \lesssim \varepsilon_{\text{mach}}\|A\|$.
- The floating-point storage error is $\sim \varepsilon_{\text{mach}}\|A\|$.
- The uncertainty ball overlaps $\Sigma$. The computer cannot tell whether $A$
  is singular or not.

This is why ill-conditioned systems are fundamentally hard. It is not a failure
of the algorithm. The problem itself is on the boundary of being unsolvable.
:::::

The figure below illustrates this geometry. Each matrix $A_i$ sits at distance
$\sigma_{\min}(A_i)$ from $\Sigma$. The dashed circles show the floating-point
uncertainty ball of radius $\varepsilon_{\text{mach}}\|A\|$. For $A_3$, the
ball overlaps $\Sigma$.

```{code-cell} python
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch

fig, ax = plt.subplots(1, 1, figsize=(9, 5.5))

# --- Singular set as a smooth curve (codimension-1 surface in the schematic) ---
t = np.linspace(-1.0, 6.0, 500)
# A gently curving "wall" representing Σ = {det(A) = 0}
curve_x = 1.8 + 0.25 * np.sin(0.9 * t)
curve_y = t
ax.fill_betweenx(curve_y, -0.5, curve_x,
                 color='#d62728', alpha=0.08)
ax.plot(curve_x, curve_y, color='#d62728', lw=2.5, zorder=3)
ax.text(0.55, 5.2, r'$\Sigma = \{\det(A) = 0\}$',
        fontsize=12, color='#d62728', ha='center', fontweight='bold')

# Floating-point ball radius (same for all, since ||A|| ~ same scale)
eps_r = 0.38

# --- A1: well-conditioned, far from Σ ---
A1 = np.array([6.0, 4.0])
# Distance to curve at this y-value
x_curve_at_A1 = 1.8 + 0.25 * np.sin(0.9 * A1[1])
dist1 = A1[0] - x_curve_at_A1

ax.plot(*A1, 'o', color='#1f77b4', ms=9, zorder=5)
ax.text(A1[0] + 0.15, A1[1] + 0.35, r'$A_1$',
        fontsize=13, color='#1f77b4', fontweight='bold')
# Uncertainty ball
c1 = plt.Circle(A1, eps_r, fill=True, facecolor='#1f77b4', alpha=0.12,
                edgecolor='#1f77b4', ls='--', lw=1.5, zorder=4)
ax.add_patch(c1)
# Distance bracket
ax.annotate('', xy=(x_curve_at_A1 + 0.05, A1[1]), xytext=(A1[0] - 0.05, A1[1]),
            arrowprops=dict(arrowstyle='<->', color='#1f77b4', lw=1.5,
                            shrinkA=0, shrinkB=0))
ax.text((A1[0] + x_curve_at_A1) / 2, A1[1] - 0.35,
        r'$\sigma_{\min}(A_1)$', fontsize=10, color='#1f77b4', ha='center')

# --- A2: ill-conditioned, close to Σ ---
A2 = np.array([2.85, 2.0])
x_curve_at_A2 = 1.8 + 0.25 * np.sin(0.9 * A2[1])
dist2 = A2[0] - x_curve_at_A2

ax.plot(*A2, 'o', color='#ff7f0e', ms=9, zorder=5)
ax.text(A2[0] + 0.05, A2[1] + 0.4, r'$A_2$',
        fontsize=13, color='#ff7f0e', fontweight='bold')
# Uncertainty ball
c2 = plt.Circle(A2, eps_r, fill=True, facecolor='#ff7f0e', alpha=0.12,
                edgecolor='#ff7f0e', ls='--', lw=1.5, zorder=4)
ax.add_patch(c2)
# Distance bracket
ax.annotate('', xy=(x_curve_at_A2 + 0.05, A2[1]), xytext=(A2[0] - 0.05, A2[1]),
            arrowprops=dict(arrowstyle='<->', color='#ff7f0e', lw=1.5,
                            shrinkA=0, shrinkB=0))
ax.text((A2[0] + x_curve_at_A2) / 2, A2[1] - 0.35,
        r'$\sigma_{\min}(A_2)$', fontsize=10, color='#ff7f0e', ha='center')

# --- A3: numerically singular, ball overlaps Σ ---
y3 = 0.5
x_curve_at_A3 = 1.8 + 0.25 * np.sin(0.9 * y3)
# sigma_min is small but visible; the eps-ball (radius 0.55) clearly overlaps Σ
dist3 = 0.35
eps_r3 = 0.55  # larger uncertainty ball for this matrix
A3 = np.array([x_curve_at_A3 + dist3, y3])

ax.plot(*A3, 'o', color='#2ca02c', ms=9, zorder=5)
ax.text(A3[0] + 0.55, A3[1] - 0.05, r'$A_3$',
        fontsize=13, color='#2ca02c', fontweight='bold')
# Uncertainty ball -- this one overlaps Σ
c3 = plt.Circle(A3, eps_r3, fill=True, facecolor='#2ca02c', alpha=0.15,
                edgecolor='#2ca02c', ls='--', lw=1.5, zorder=4)
ax.add_patch(c3)
# Distance arrow (horizontal, at the point's y-level)
ax.annotate('', xy=(x_curve_at_A3 + 0.03, A3[1]),
            xytext=(A3[0] - 0.03, A3[1]),
            arrowprops=dict(arrowstyle='<->', color='#2ca02c', lw=1.5,
                            shrinkA=0, shrinkB=0))
# Label above the arrow, shifted right to avoid overlapping the boundary
ax.text(A3[0] + 0.1, A3[1] + 0.25,
        r'$\sigma_{\min}$', fontsize=10, color='#2ca02c', ha='left')

# --- Annotations ---
# Legend-style text box
legend_text = (
    r'$A_1$: well-conditioned ($\kappa \approx 1$)' + '\n'
    r'$A_2$: ill-conditioned ($\kappa \gg 1$)' + '\n'
    r'$A_3$: numerically singular ($\kappa \gtrsim 1/\varepsilon_{\mathrm{mach}}$)'
)
ax.text(7.8, 1.3, legend_text, fontsize=9.5, va='center', ha='right',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='#f7f7f7',
                  edgecolor='#cccccc', alpha=0.95),
        linespacing=1.8)

# Ball label
ax.text(7.8, -0.2,
        r'Dashed circles = $\varepsilon_{\mathrm{mach}}\|A\|$'
        '\n(floating-point uncertainty)',
        fontsize=9, color='#666666', ha='right', va='top',
        linespacing=1.6)

# Axis label
ax.text(4.5, -0.8, r'Space of matrices $\mathbb{R}^{n \times n}$  (schematic)',
        fontsize=11, ha='center', color='#444444')

ax.set_xlim(-0.5, 8.2)
ax.set_ylim(-1.0, 5.8)
ax.set_aspect('equal')
ax.axis('off')
ax.set_title('Distance to singularity in the operator norm', fontsize=13, pad=12)

plt.tight_layout()
plt.show()
```


## Residuals and Backward Error

:::{prf:definition} Residual
:label: def-residual

Given a linear system $A\mathbf{x} = \mathbf{b}$ and a computed solution
$\hat{\mathbf{x}}$, the **residual** is

$$
\mathbf{r} = \mathbf{b} - A\hat{\mathbf{x}}
$$
:::

The residual is cheap to compute and directly measures the backward error.

:::{prf:proposition} The Residual Measures Backward Error
:label: prop-residual-backward-error

The relative backward error (perturbing $\mathbf{b}$ only) equals
$\|\mathbf{r}\|/\|\mathbf{b}\|$. That is, $\hat{\mathbf{x}}$ is the exact
solution of the perturbed system $A\hat{\mathbf{x}} = \mathbf{b} + \delta\mathbf{b}$
with $\delta\mathbf{b} = -\mathbf{r}$.
:::

:::{prf:proof}
:class: dropdown

By definition, $\mathbf{r} = \mathbf{b} - A\hat{\mathbf{x}}$, so rearranging:

$$
A\hat{\mathbf{x}} = \mathbf{b} - \mathbf{r}
$$

Setting $\delta\mathbf{b} = -\mathbf{r}$, we have
$A\hat{\mathbf{x}} = \mathbf{b} + \delta\mathbf{b}$ exactly. The computed
solution $\hat{\mathbf{x}}$ is the *exact* solution of this perturbed system,
with perturbation size

$$
\frac{\|\delta\mathbf{b}\|}{\|\mathbf{b}\|} = \frac{\|\mathbf{r}\|}{\|\mathbf{b}\|}
$$

This is the smallest perturbation to $\mathbf{b}$ (with $A$ fixed) that makes
$\hat{\mathbf{x}}$ an exact solution, so $\|\mathbf{r}\|/\|\mathbf{b}\|$ is
the relative backward error.
:::

:::{admonition} Warning: Small Residual $\neq$ Small Error
:class: warning

Combining the residual with the sensitivity theorem
([](#thm-sensitivity-linear-systems)):

$$
\frac{\|\hat{\mathbf{x}} - \mathbf{x}\|}{\|\mathbf{x}\|} \lesssim \kappa(A) \cdot \frac{\|\mathbf{r}\|}{\|\mathbf{b}\|}
$$

For ill-conditioned $A$, a tiny residual can hide a large error. The residual
measures backward error; you need to multiply by $\kappa(A)$ to bound forward
error.
:::

## Forward Stability

Now we turn to algorithms. We have two candidates for what "good algorithm"
should mean. The first is the most natural wish.

:::{prf:definition} Forward Stable Algorithm
:label: def-forward-stable

An algorithm for solving $A\mathbf{x} = \mathbf{b}$ is **forward stable** if
the computed solution $\hat{\mathbf{x}}$ satisfies:

$$
\frac{\|\hat{\mathbf{x}} - \mathbf{x}\|}{\|\mathbf{x}\|} = O(\varepsilon_{\text{mach}})
$$

That is, the relative forward error is of order machine epsilon, regardless
of the conditioning of $A$.
:::

This sounds ideal. Why not demand it?

:::{prf:proposition} Forward Stability Is Unattainable
:label: prop-forward-stability-impossible

No algorithm for solving $A\mathbf{x} = \mathbf{b}$ in floating-point
arithmetic can be forward stable for all nonsingular $A$.

:::

:::{prf:proof}
:class: dropdown

The obstacle is not algorithmic; it is inherent in floating-point
representation.

**The input is already perturbed.** When we store $A$ and $\mathbf{b}$ in
floating point, we do not have the true $A$ and $\mathbf{b}$. We have
$\tilde{A}$ and $\tilde{\mathbf{b}}$ with

$$
\frac{\|\tilde{A} - A\|}{\|A\|} = O(\varepsilon_{\text{mach}}), \qquad
\frac{\|\tilde{\mathbf{b}} - \mathbf{b}\|}{\|\mathbf{b}\|} = O(\varepsilon_{\text{mach}})
$$

Even before any computation begins, the problem has been perturbed. By the
sensitivity theorem ([](#thm-sensitivity-linear-systems)), the exact solution
of the perturbed problem already differs from $\mathbf{x}$ by

$$
\frac{\|\tilde{\mathbf{x}} - \mathbf{x}\|}{\|\mathbf{x}\|} = O(\kappa(A) \cdot \varepsilon_{\text{mach}})
$$

No algorithm can undo this. The input perturbation gets amplified by $\kappa(A)$,
and no amount of clever arithmetic can recover the lost information. An
algorithm that achieved $O(\varepsilon_{\text{mach}})$ forward error regardless
of $\kappa(A)$ would have to "know" the true $A$ despite only having access to
$\tilde{A}$.

**Forward stability requires $\kappa(A) = O(1)$**, which is a condition on the
problem, not the algorithm. For well-conditioned problems, forward and backward
stability coincide. For ill-conditioned problems, forward stability is
unachievable by any algorithm.
:::

## Backward Stability: The Right Standard

Since we cannot control forward error directly, we control what we can: the
backward error. This leads to the Higham standard.

:::{prf:definition} Backward Stable Algorithm
:label: def-backward-stable-linalg

An algorithm for solving $A\mathbf{x} = \mathbf{b}$ is **backward stable** if
the computed solution $\hat{\mathbf{x}}$ satisfies:

$$
(A + \delta A)\hat{\mathbf{x}} = \mathbf{b} + \delta\mathbf{b}, \quad
\frac{\|\delta A\|}{\|A\|} = O(\varepsilon_{\text{mach}}), \quad
\frac{\|\delta\mathbf{b}\|}{\|\mathbf{b}\|} = O(\varepsilon_{\text{mach}})
$$

That is, $\hat{\mathbf{x}}$ is the *exact* solution of a *nearby* problem.
:::

**Why this is the right standard.** A backward stable algorithm introduces
perturbations no larger than those already present from storing $A$ and
$\mathbf{b}$ in floating point. The algorithm does not make the situation any
worse than it already is. You cannot ask for more than that.

Combined with the sensitivity theorem ([](#thm-sensitivity-linear-systems)):

$$
\frac{\|\hat{\mathbf{x}} - \mathbf{x}\|}{\|\mathbf{x}\|} \lesssim \kappa(A) \cdot \varepsilon_{\text{mach}}
$$

This is the **optimal forward error bound**. It is optimal because:

1. The factor $\varepsilon_{\text{mach}}$ comes from the backward error, which is as small as possible.
2. The factor $\kappa(A)$ is intrinsic to the problem. No algorithm can remove it.

A backward stable algorithm achieves the best forward error that any algorithm
can achieve on that problem. If the answer is not accurate enough, the problem
is to blame, not the algorithm.

:::{prf:remark} Which Algorithms Are Backward Stable?
:label: rmk-backward-stable-algorithms
:class: dropdown

- **Householder QR** is backward stable (see
  [](#thm-householder-stability)). This is the gold standard.
- **LU with partial pivoting** is backward stable in practice, though the
  theoretical worst case allows growth factor $2^{n-1}$ (see the
  [stability notebook](../notebooks/stability-ge.ipynb)).
- **Classical Gram-Schmidt** is *not* backward stable: the loss of
  orthogonality scales with $\kappa(A)$, not $\varepsilon_{\text{mach}}$.
:::

## Practical Guideline: Always Check the Condition Number

Since the forward error satisfies

$$
\frac{\|\hat{\mathbf{x}} - \mathbf{x}\|}{\|\mathbf{x}\|} \lesssim \kappa(A) \cdot \varepsilon_{\text{mach}}
$$

we should **always estimate $\kappa(A)$ before trusting the solution**. But how?

### The Challenge

Computing $\kappa(A) = \|A\| \|A^{-1}\|$ exactly requires $A^{-1}$, which costs $O(n^3)$ operations, as expensive as solving the system! We need a cheaper approach.

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

**Cost:** Each iteration requires two triangular solves. Typically converges in 2-5 iterations, so total cost is $O(n^2)$, much cheaper than the $O(n^3)$ factorization.

See the [Condition Number Estimation notebook](./CondNumberEst.ipynb) for a Python implementation.

:::{prf:remark} Why the 1-norm is sufficient
:label: rmk-1-norm-sufficient
:class: dropdown

Hager's algorithm estimates $\kappa_1(A)$, not $\kappa_2(A)$. This is fine for
two reasons.

First, all matrix norms on $\mathbb{R}^{n \times n}$ are equivalent. The sharp
bounds between the 1-norm and 2-norm are:

$$
\frac{1}{\sqrt{n}}\|A\|_1 \leq \|A\|_2 \leq \sqrt{n}\,\|A\|_1
$$

so $\kappa_1(A)$ and $\kappa_2(A)$ can differ by at most a factor of $n$. Since
the forward error bound ([](#prop-residual-backward-error)) in the 1-norm is

$$
\frac{\|\hat{\mathbf{x}} - \mathbf{x}\|_1}{\|\mathbf{x}\|_1} \lesssim \kappa_1(A) \cdot \frac{\|\mathbf{r}\|_1}{\|\mathbf{b}\|_1}
$$

we can convert to a 2-norm bound using $\|\mathbf{v}\|_2 \leq \|\mathbf{v}\|_1
\leq \sqrt{n}\,\|\mathbf{v}\|_2$ for vectors and the matrix norm equivalence
above. This gives

$$
\frac{\|\hat{\mathbf{x}} - \mathbf{x}\|_2}{\|\mathbf{x}\|_2} \lesssim n \cdot \kappa_1(A) \cdot \frac{\|\mathbf{r}\|_2}{\|\mathbf{b}\|_2}
$$

The extra factor of $n$ is pessimistic (it comes from worst-casing the norm
conversion at every step). In practice the 1-norm bound is already a good
estimate of the 2-norm error. For the "digits lost" rule of thumb, a factor
of $n$ changes $\log_{10}\kappa$ by at most $\log_{10} n$, which is negligible.

Second, the 1-norm is cheap: $\|A\|_1 = \max_j \sum_i |a_{ij}|$ costs $O(n^2)$
and requires no singular value computation. The 2-norm $\|A\|_2 = \sigma_{\max}$
is more expensive. Since Hager's algorithm is a practical tool for cheap error
bounds, the 1-norm is the natural choice.
:::

:::{admonition} Practical Workflow
:class: tip

**Solving $A\mathbf{x} = \mathbf{b}$ responsibly:**

1. **Factor:** Compute $A = LU$ (or $A = QR$).
2. **Solve:** Use the factorization to obtain $\hat{\mathbf{x}}$.
3. **Estimate $\kappa(A)$:** Apply Hager's algorithm to the factorization. This costs only $O(n^2)$ extra. For QR, the same idea works: estimate $\|R^{-1}\|$ via triangular solves, then $\kappa(A) \approx \|A\| \cdot \|R^{-1}\|$.
4. **Compute the residual:** $\mathbf{r} = \mathbf{b} - A\hat{\mathbf{x}}$.
5. **Forward error bound:** By [](#prop-residual-backward-error) and [](#thm-sensitivity-linear-systems):

$$
\frac{\|\hat{\mathbf{x}} - \mathbf{x}\|}{\|\mathbf{x}\|} \lesssim \kappa(A) \cdot \frac{\|\mathbf{r}\|}{\|\mathbf{b}\|}
$$

If $\kappa(A) \cdot \varepsilon_{\text{mach}} \gtrsim 1$, the answer may be meaningless.

Steps 3-5 are essentially free compared to the $O(n^3)$ factorization, and they give you a **computable forward error bound**.
:::
