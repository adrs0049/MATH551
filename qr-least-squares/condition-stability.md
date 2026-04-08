---
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

## Consequences of the Condition Number

Recall that the [condition number](linear-algebra-fundamentals.md#def-condition-number)
$\kappa(A) = \|A\|\|A^{-1}\|$ measures the ratio of maximum to minimum stretch
of a unit vector. The sensitivity theorem tells us what this means in practice.

**Rule of thumb:** Expect to lose $\log_{10}\kappa(A)$ digits of accuracy.

| Condition Number | Digits Lost |
|-----------------|-------------|
| $\kappa \approx 10^k$ | ~$k$ digits |
| $\kappa \gtrsim 1/\varepsilon_{\text{mach}} \approx 10^{16}$ | All digits |

See the [Condition Numbers and Lost Digits notebook](../notebooks/condition-number-digits.ipynb) for concrete demonstrations with diagonal, Hilbert, and Vandermonde matrices.

But what does it *mean* for a matrix to be ill-conditioned? The next section provides the key insight.

## The Deep Insight: Numerically Singular Matrices

:::::{admonition} Extension: When Ill-Conditioned Means Singular (Demmel)
:class: note

A matrix with $\kappa(A) \gtrsim 1/\varepsilon_{\text{mach}}$ is **numerically indistinguishable from a singular matrix**.

::::{dropdown} Why? (Demmel's Perspective)
Consider the "distance to singularity": how much do we need to perturb $A$ to make it singular?

For the 2-norm, this distance is exactly $\sigma_{\min}(A)$, the smallest singular value. In relative terms:

$$
\frac{\text{distance to nearest singular matrix}}{\|A\|} = \frac{\sigma_{\min}}{\sigma_{\max}} = \frac{1}{\kappa(A)}
$$

Now consider floating-point arithmetic. Every matrix $A$ is stored with relative error $\sim \varepsilon_{\text{mach}}$. The computer does not see $A$. It sees $A + E$ where $\|E\|/\|A\| \sim \varepsilon_{\text{mach}}$.

**If $\kappa(A) \gtrsim 1/\varepsilon_{\text{mach}}$:**
- The distance to singularity is $\lesssim \varepsilon_{\text{mach}}$
- The storage error is $\sim \varepsilon_{\text{mach}}$
- The computer cannot distinguish $A$ from a singular matrix!

This is why ill-conditioned systems are fundamentally hard. It is not because of
bad algorithms. The *problem itself* is on the edge of being unsolvable.
::::
:::::

The picture below illustrates this geometry. The set of singular matrices is
closed in the operator norm topology. Every matrix $A$ sits at distance
$\sigma_{\min}(A) = \|A\|/\kappa(A)$ from this set. When the condition number
is large, $A$ sits close to the boundary. When $\kappa(A) \gtrsim
1/\varepsilon_{\text{mach}}$, the floating-point "uncertainty ball" around $A$
overlaps the singular set.

```{code-cell} python
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

fig, ax = plt.subplots(1, 1, figsize=(8, 5))

# Draw the singular set as a filled region on the left
singular_x = np.linspace(-0.5, 2.0, 300)
singular_upper = 4.5 + 0.3 * np.sin(3 * singular_x)
singular_lower = -0.5 - 0.3 * np.sin(3 * singular_x + 1)
ax.fill_between(singular_x, singular_lower, singular_upper,
                color='#d62728', alpha=0.15, label='Singular matrices')
ax.plot(singular_x, singular_upper, color='#d62728', lw=2)
ax.plot(singular_x, singular_lower, color='#d62728', lw=2)
# Boundary label
ax.text(0.75, 4.9, 'Singular matrices', fontsize=11, color='#d62728',
        ha='center', style='italic')

# Well-conditioned matrix: far from the boundary
A_well = np.array([5.5, 2.0])
dist_well = 3.2
eps_ball_well = 0.4
ax.plot(*A_well, 'o', color='#1f77b4', ms=8, zorder=5)
ax.text(A_well[0] + 0.15, A_well[1] + 0.3, r'$A_1$: $\kappa \approx 1$',
        fontsize=11, color='#1f77b4')
circle_well = plt.Circle(A_well, eps_ball_well, fill=False,
                          color='#1f77b4', ls='--', lw=1.5, alpha=0.7)
ax.add_patch(circle_well)
# Distance arrow
ax.annotate('', xy=(2.0, 2.0), xytext=A_well,
            arrowprops=dict(arrowstyle='<->', color='#1f77b4', lw=1.5))
ax.text(3.75, 2.5, r'$\sigma_{\min} \gg \varepsilon_{\mathrm{mach}}\|A\|$',
        fontsize=9, color='#1f77b4', ha='center')

# Ill-conditioned matrix: close to the boundary
A_ill = np.array([2.7, 3.2])
eps_ball_ill = 0.4
ax.plot(*A_ill, 'o', color='#ff7f0e', ms=8, zorder=5)
ax.text(A_ill[0] + 0.15, A_ill[1] + 0.3, r'$A_2$: $\kappa \gg 1$',
        fontsize=11, color='#ff7f0e')
circle_ill = plt.Circle(A_ill, eps_ball_ill, fill=False,
                          color='#ff7f0e', ls='--', lw=1.5, alpha=0.7)
ax.add_patch(circle_ill)
# Distance arrow
ax.annotate('', xy=(2.0, 3.2), xytext=A_ill,
            arrowprops=dict(arrowstyle='<->', color='#ff7f0e', lw=1.5))
ax.text(2.35, 3.55, r'$\sigma_{\min}$', fontsize=9, color='#ff7f0e', ha='center')

# Numerically singular matrix: ball overlaps singular set
A_num = np.array([2.25, 1.0])
eps_ball_num = 0.5
ax.plot(*A_num, 'o', color='#2ca02c', ms=8, zorder=5)
ax.text(A_num[0] + 0.2, A_num[1] - 0.45,
        r'$A_3$: $\kappa \gtrsim 1/\varepsilon_{\mathrm{mach}}$',
        fontsize=11, color='#2ca02c')
circle_num = plt.Circle(A_num, eps_ball_num, fill=False,
                          color='#2ca02c', ls='--', lw=1.5, alpha=0.7)
ax.add_patch(circle_num)

# Dashed balls label
ax.text(6.2, 4.3, r'Dashed circles: $\varepsilon_{\mathrm{mach}}\|A\|$',
        fontsize=9, color='gray', ha='center',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='gray', alpha=0.8))
ax.text(6.2, 3.8, '(floating-point uncertainty)',
        fontsize=9, color='gray', ha='center')

ax.set_xlim(-0.7, 8)
ax.set_ylim(-0.8, 5.5)
ax.set_aspect('equal')
ax.set_xlabel(r'$\mathbb{R}^{n \times n}$ (schematic)', fontsize=11)
ax.set_title('Distance to singularity in the operator norm', fontsize=12)
ax.axis('off')

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
