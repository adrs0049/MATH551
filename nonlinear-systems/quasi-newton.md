---
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/quasi-newton.pdf
    id: nonlinear-systems-quasi-newton-pdf
downloads:
  - id: nonlinear-systems-quasi-newton-pdf
    title: Download PDF
---

# Quasi-Newton Methods

:::{tip} Big Idea
Computing the Jacobian costs $\mathcal{O}(n^2)$ evaluations per iteration.
**Quasi-Newton methods** build an approximation
$B_k \approx D\mathbf{F}(\mathbf{x}_k)$ and update it cheaply using
information from the iteration. The trade-off: superlinear convergence instead
of quadratic, but much cheaper iterations.
:::

Each Newton iteration requires evaluating $n^2$ partial derivatives and an
$\mathcal{O}(n^3)$ LU factorization. For large $n$, or when derivatives are
expensive, this is prohibitive. Quasi-Newton methods replace the exact Jacobian
with an approximation $B_k \approx D\mathbf{F}(\mathbf{x}_k)$ that is updated
cheaply ($\mathcal{O}(n^2)$ per step) using information from the iteration
itself.

## From the Secant Method to Systems

In one dimension, the **secant method** avoids computing $f'$ by approximating
it from two previous iterates:

$$
f'(x_k) \approx \frac{f(x_k) - f(x_{k-1})}{x_k - x_{k-1}}
$$

This replaces the tangent line (Newton) with a secant line through the two most
recent points. The secant method converges with order $p \approx 1.618$ (the
golden ratio), which is between linear and quadratic.

For a system $\mathbf{F}: \mathbb{R}^n \to \mathbb{R}^n$, we want to
generalize this idea: build an approximation $B_{k+1}$ to the Jacobian using
the step $\mathbf{s}_k = \mathbf{x}_{k+1} - \mathbf{x}_k$ and the
corresponding change in $\mathbf{F}$,
$\mathbf{y}_k = \mathbf{F}(\mathbf{x}_{k+1}) - \mathbf{F}(\mathbf{x}_k)$.
The natural requirement is that $B_{k+1}$ should reproduce this finite
difference exactly.

:::{prf:definition} The Secant Condition
:label: def-secant-condition

The approximation $B_{k+1}$ must satisfy:

$$
B_{k+1}\mathbf{s}_k = \mathbf{y}_k
$$

where $\mathbf{s}_k = \mathbf{x}_{k+1} - \mathbf{x}_k$ and
$\mathbf{y}_k = \mathbf{F}(\mathbf{x}_{k+1}) - \mathbf{F}(\mathbf{x}_k)$.
:::

This says: the approximation should be exact in the direction we just moved.
But the secant condition gives only $n$ equations for $n^2$ unknowns (the
entries of $B_{k+1}$), so infinitely many matrices satisfy it. We need an
additional criterion to select one.

## Broyden's Method

The simplest and most natural choice: among all matrices satisfying the secant
condition, pick the one that changes $B_k$ as little as possible. This leads to
**Broyden's method** (1965).

:::{prf:definition} Broyden's Update
:label: def-broyden-update

Given the current approximation $B_k$, the updated approximation is:

$$
B_{k+1} = B_k + \frac{(\mathbf{y}_k - B_k\mathbf{s}_k)\mathbf{s}_k^T}{\mathbf{s}_k^T\mathbf{s}_k}
$$
:::

:::{prf:proof} Derivation
:class: dropdown

We solve the constrained minimization problem:

$$
\min_{B_{k+1}} \|B_{k+1} - B_k\|_F \quad \text{subject to} \quad B_{k+1}\mathbf{s}_k = \mathbf{y}_k
$$

The constraint requires
$(B_{k+1} - B_k)\mathbf{s}_k = \mathbf{y}_k - B_k\mathbf{s}_k$. The minimum
Frobenius norm perturbation satisfying a single linear constraint is a rank-1
matrix:

$$
B_{k+1} - B_k = \frac{(\mathbf{y}_k - B_k\mathbf{s}_k)\mathbf{s}_k^T}{\mathbf{s}_k^T\mathbf{s}_k}
$$

To verify: multiply both sides by $\mathbf{s}_k$ on the right to recover the
constraint. The factor $\mathbf{s}_k^T/(\mathbf{s}_k^T\mathbf{s}_k)$ is a
projection that ensures the perturbation acts only in the direction
$\mathbf{s}_k$, leaving $B_k$ unchanged in all orthogonal directions.
:::

The update is a **rank-1 perturbation** of $B_k$, which means:
- $B_{k+1}$ agrees with $B_k$ on all directions orthogonal to $\mathbf{s}_k$.
- $B_{k+1}$ is corrected only in the direction we just moved, where we have
  new information.

### The Algorithm

:::{prf:algorithm} Broyden's Method
:label: alg-broyden-method

**Input:** $\mathbf{F}$, initial guess $\mathbf{x}_0$, initial
$B_0 \approx D\mathbf{F}(\mathbf{x}_0)$, tolerance $\varepsilon$

**Output:** Approximate root $\mathbf{x}$

1. **for** $k = 0, 1, 2, \ldots$:
2. $\qquad$ Solve $B_k \Delta\mathbf{x} = -\mathbf{F}(\mathbf{x}_k)$
3. $\qquad$ $\mathbf{x}_{k+1} \gets \mathbf{x}_k + \Delta\mathbf{x}$
4. $\qquad$ **if** $\|\Delta\mathbf{x}\| < \varepsilon$: **return** $\mathbf{x}_{k+1}$
5. $\qquad$ $\mathbf{s}_k \gets \Delta\mathbf{x}$
6. $\qquad$ $\mathbf{y}_k \gets \mathbf{F}(\mathbf{x}_{k+1}) - \mathbf{F}(\mathbf{x}_k)$
7. $\qquad$ $B_{k+1} \gets B_k + \frac{(\mathbf{y}_k - B_k\mathbf{s}_k)\mathbf{s}_k^T}{\mathbf{s}_k^T\mathbf{s}_k}$
:::

### Efficient Implementation via Sherman-Morrison

Each iteration requires solving $B_k \Delta\mathbf{x} = -\mathbf{F}$, which
naively costs $\mathcal{O}(n^3)$ for LU factorization. Since $B_{k+1}$ differs
from $B_k$ by a rank-1 matrix, we can instead maintain $H_k = B_k^{-1}$
directly and update it using the **Sherman-Morrison formula**:

:::{prf:proposition} Sherman-Morrison Formula
:label: prop-sherman-morrison

If $A$ is invertible and $A + \mathbf{u}\mathbf{v}^T$ is invertible, then:

$$
(A + \mathbf{u}\mathbf{v}^T)^{-1} = A^{-1} - \frac{A^{-1}\mathbf{u}\mathbf{v}^T A^{-1}}{1 + \mathbf{v}^T A^{-1}\mathbf{u}}
$$
:::

Applied to Broyden's update, the inverse approximation $H_k = B_k^{-1}$
satisfies:

$$
H_{k+1} = H_k + \frac{(\mathbf{s}_k - H_k\mathbf{y}_k)\mathbf{s}_k^T H_k}{\mathbf{s}_k^T H_k \mathbf{y}_k}
$$

With this approach, each iteration costs $\mathcal{O}(n^2)$ (matrix-vector
products for the update, and $\Delta\mathbf{x} = -H_k\mathbf{F}$ is a
matrix-vector multiply instead of a linear solve).

## Convergence

:::{prf:theorem} Superlinear Convergence of Broyden's Method
:label: thm-superlinear-convergence-broyden

Suppose $\mathbf{F}(\mathbf{x}^*) = \mathbf{0}$,
$D\mathbf{F}(\mathbf{x}^*)$ is nonsingular, and $D\mathbf{F}$ is Lipschitz
continuous near $\mathbf{x}^*$. If $\mathbf{x}_0$ is sufficiently close to
$\mathbf{x}^*$ and $B_0$ is sufficiently close to
$D\mathbf{F}(\mathbf{x}^*)$, then Broyden's method converges
**superlinearly**:

$$
\lim_{k \to \infty} \frac{\|\mathbf{x}_{k+1} - \mathbf{x}^*\|}{\|\mathbf{x}_k - \mathbf{x}^*\|} = 0
$$
:::

:::{prf:proof}
:class: dropdown

Write $\mathbf{e}_k = \mathbf{x}_k - \mathbf{x}^*$ and
$E_k = B_k - D\mathbf{F}(\mathbf{x}^*)$. The Newton step with exact Jacobian
would give $\mathbf{e}_{k+1}^{\text{Newton}} = \mathcal{O}(\|\mathbf{e}_k\|^2)$
(quadratic). With Broyden's approximation, the error satisfies:

$$
\mathbf{e}_{k+1} = \mathbf{e}_k + \Delta\mathbf{x}_k
= \mathbf{e}_k - B_k^{-1}\mathbf{F}(\mathbf{x}_k)
$$

Adding and subtracting $D\mathbf{F}(\mathbf{x}^*)^{-1}\mathbf{F}(\mathbf{x}_k)$:

$$
\mathbf{e}_{k+1} = \underbrace{\mathbf{e}_k - D\mathbf{F}(\mathbf{x}^*)^{-1}\mathbf{F}(\mathbf{x}_k)}_{\mathcal{O}(\|\mathbf{e}_k\|^2) \text{ (Newton residual)}} + \underbrace{(D\mathbf{F}(\mathbf{x}^*)^{-1} - B_k^{-1})\mathbf{F}(\mathbf{x}_k)}_{\text{Jacobian approximation error}}
$$

The first term is $\mathcal{O}(\|\mathbf{e}_k\|^2)$ from the standard Newton
analysis. The second term is bounded by
$\|D\mathbf{F}(\mathbf{x}^*)^{-1}\| \cdot \|E_k\| \cdot \|\mathbf{e}_k\|$
(since $\|\mathbf{F}(\mathbf{x}_k)\| = \mathcal{O}(\|\mathbf{e}_k\|)$).

The key result (Dennis and Moré, 1974) is that the Broyden update satisfies:

$$
\|E_{k+1}\mathbf{s}_k\| \leq \|E_k\mathbf{s}_k\| + \mathcal{O}(\|\mathbf{e}_k\|^2)
$$

The error in $B_k$ does not grow in the direction $\mathbf{s}_k$, and the
$\mathcal{O}(\|\mathbf{e}_k\|^2)$ term shrinks as we converge. A careful
analysis (the Dennis-Moré characterization theorem) shows that the Jacobian
approximation error satisfies:

$$
\frac{\|E_k \mathbf{s}_k\|}{\|\mathbf{s}_k\|} \to 0
$$

This means $B_k$ becomes increasingly accurate *in the directions that matter*
(the step directions), even though $\|E_k\|$ may not converge to zero overall.
Combined with the $\mathcal{O}(\|\mathbf{e}_k\|^2)$ Newton residual, this gives
$\|\mathbf{e}_{k+1}\|/\|\mathbf{e}_k\| \to 0$, i.e., superlinear convergence.
:::

The convergence is superlinear but not quadratic because $B_k$ is only an
approximation to the Jacobian. Newton's method uses the exact Jacobian at each
step, giving $\mathcal{O}(\|\mathbf{e}_k\|^2)$ error. Broyden's method adds a
Jacobian approximation error that is $o(\|\mathbf{e}_k\|)$ (vanishing relative
to the error) but not $\mathcal{O}(\|\mathbf{e}_k\|^2)$. The approximation
improves in the step directions as the iteration progresses, which is enough for
superlinear but not quadratic convergence.

## Broyden's "Good" and "Bad" Methods

The update above is called **Broyden's first method** (or "good" Broyden). It
updates $B_k$ to match the secant condition $B_{k+1}\mathbf{s}_k = \mathbf{y}_k$.

**Broyden's second method** (or "bad" Broyden) instead updates the inverse
$H_k = B_k^{-1}$ to satisfy the *inverse* secant condition
$H_{k+1}\mathbf{y}_k = \mathbf{s}_k$, using the analogous minimum-change
principle:

$$
H_{k+1} = H_k + \frac{(\mathbf{s}_k - H_k\mathbf{y}_k)\mathbf{y}_k^T}{\mathbf{y}_k^T\mathbf{y}_k}
$$

The names are historical and misleading: neither method is universally better.
Broyden's first method tends to work better when $D\mathbf{F}$ is well
conditioned; the second method can be preferable when $D\mathbf{F}$ is ill
conditioned, since it avoids inverting the approximation.

## Error-Oriented Quasi-Newton (QNERR)

Broyden's method is only locally convergent: it has no mechanism to detect or
recover from a bad approximation. In Deuflhard's framework (see
[Globalization](globalization.md)), the natural combination of quasi-Newton with
globalization is **QNERR**: the error-oriented quasi-Newton method.

In [NLEQ_ERR](globalization.md#alg-nleq-err), once the contraction factor
$\theta_k$ drops below a threshold $\theta_{\max}$ (typically 0.5), the
iteration is converging well and recomputing the Jacobian is wasteful. QNERR
**freezes the Jacobian** and applies rank-1 corrections using the iteration
history, while continuing to monitor $\theta_k$. If the contraction degrades
($\theta_k > \theta_{\max}$), QNERR hands control back to NLEQ_ERR, which
recomputes the Jacobian.

This gives a hybrid strategy:

- **Phase 1 (NLEQ_ERR):** full Newton with damping until contraction is
  established ($\theta < \theta_{\max}$).
- **Phase 2 (QNERR):** freeze the Jacobian, apply rank-1 corrections, monitor
  $\theta$. If contraction degrades, switch back to Phase 1.

Unlike standalone Broyden, QNERR always starts from the true Jacobian (computed
in the Newton phase) and has a built-in safety net through $\theta_k$
monitoring. It also inherits the
[affine invariance](globalization.md#def-affine-invariance) of NLEQ_ERR, since
it monitors $\|\Delta\mathbf{x}\|$ rather than $\|\mathbf{F}\|$.
