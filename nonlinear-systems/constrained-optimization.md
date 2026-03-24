---
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/constrained-optimization.pdf
    id: nonlinear-systems-constrained-optimization-pdf
downloads:
  - id: nonlinear-systems-constrained-optimization-pdf
    title: Download PDF
---

# Constrained Optimization and Lagrange Multipliers

:::{tip} Big Idea
Constrained optimization problems become nonlinear systems via **Lagrange multipliers**. The KKT conditions—necessary conditions for optimality—are just a system $\mathbf{F}(\mathbf{x}, \boldsymbol{\lambda}) = \mathbf{0}$. Newton's method solves them directly.
:::

## The Constrained Optimization Problem

We want to minimize a function subject to constraints:

$$
\min_{\mathbf{x} \in \mathbb{R}^n} f(\mathbf{x}) \quad \text{subject to} \quad \mathbf{g}(\mathbf{x}) = \mathbf{0}
$$

where $\mathbf{g}: \mathbb{R}^n \to \mathbb{R}^m$ represents $m$ equality constraints.

::::{admonition} Example: Closest Point on a Curve
:class: note

Find the point on the circle $x^2 + y^2 = 1$ closest to $(2, 1)$.

**Objective:** $f(x, y) = (x-2)^2 + (y-1)^2$ (squared distance)

**Constraint:** $g(x, y) = x^2 + y^2 - 1 = 0$

We can't just set $\nabla f = 0$—the minimum of $f$ is at $(2, 1)$, which violates the constraint!
::::

## The Geometric Insight

At a constrained minimum, the gradient of $f$ must be perpendicular to the constraint surface. Why? If $\nabla f$ had a component tangent to the constraint, we could move along the constraint and decrease $f$.

::::::{prf:theorem} The Lagrange Condition
:label: thm-lagrange-condition

At a constrained minimum $\mathbf{x}^*$, there exist multipliers $\boldsymbol{\lambda} = (\lambda_1, \ldots, \lambda_m)^T$ such that:

$$
\nabla f(\mathbf{x}^*) = \sum_{i=1}^m \lambda_i \nabla g_i(\mathbf{x}^*) = D\mathbf{g}(\mathbf{x}^*)^T \boldsymbol{\lambda}
$$

The gradient of $f$ is a linear combination of the constraint gradients.

::::{dropdown} Geometric Interpretation
The constraint surface $\mathbf{g}(\mathbf{x}) = \mathbf{0}$ has tangent space orthogonal to $\nabla g_i$ at each point. At the optimum:
- $\nabla f$ must be orthogonal to the tangent space (no improving direction)
- This means $\nabla f$ lies in the span of $\{\nabla g_1, \ldots, \nabla g_m\}$
::::
::::::

## The Lagrangian

Define the **Lagrangian function**:

$$
\mathcal{L}(\mathbf{x}, \boldsymbol{\lambda}) = f(\mathbf{x}) - \boldsymbol{\lambda}^T \mathbf{g}(\mathbf{x}) = f(\mathbf{x}) - \sum_{i=1}^m \lambda_i g_i(\mathbf{x})
$$

The necessary conditions for a constrained minimum are:

$$
\nabla_{\mathbf{x}} \mathcal{L} = \nabla f - D\mathbf{g}^T \boldsymbol{\lambda} = \mathbf{0}
$$
$$
\nabla_{\boldsymbol{\lambda}} \mathcal{L} = -\mathbf{g}(\mathbf{x}) = \mathbf{0}
$$

**This is a system of $n + m$ equations in $n + m$ unknowns!**

## KKT Conditions as a Nonlinear System

The Karush-Kuhn-Tucker (KKT) conditions for equality-constrained optimization:

::::::{prf:definition} KKT System (Equality Constraints)
:label: def-kkt-system-equality

Find $(\mathbf{x}^*, \boldsymbol{\lambda}^*)$ such that:

$$
\mathbf{F}(\mathbf{x}, \boldsymbol{\lambda}) = \begin{pmatrix} \nabla f(\mathbf{x}) - D\mathbf{g}(\mathbf{x})^T \boldsymbol{\lambda} \\ \mathbf{g}(\mathbf{x}) \end{pmatrix} = \mathbf{0}
$$

This is an $(n+m)$-dimensional nonlinear system.
::::::

::::{admonition} Example: Circle Problem as KKT System
:class: note

For the closest-point problem:
- $f(x,y) = (x-2)^2 + (y-1)^2$
- $g(x,y) = x^2 + y^2 - 1$

The KKT system is:

$$
\mathbf{F}(x, y, \lambda) = \begin{pmatrix} 2(x-2) - 2\lambda x \\ 2(y-1) - 2\lambda y \\ x^2 + y^2 - 1 \end{pmatrix} = \mathbf{0}
$$

Three equations in three unknowns $(x, y, \lambda)$.
::::

## Newton's Method for KKT Systems

Apply Newton's method to the KKT system. The Jacobian of $\mathbf{F}(\mathbf{x}, \boldsymbol{\lambda})$ has a special structure:

$$
D\mathbf{F} = \begin{pmatrix} \nabla^2_{xx}\mathcal{L} & -D\mathbf{g}^T \\ D\mathbf{g} & \mathbf{0} \end{pmatrix}
$$

where $\nabla^2_{xx}\mathcal{L} = \nabla^2 f - \sum_i \lambda_i \nabla^2 g_i$ is the Hessian of the Lagrangian with respect to $\mathbf{x}$.

:::{prf:algorithm} Newton for Constrained Optimization
:label: alg-newton-constrained-optimization

**Input:** $f$, $\mathbf{g}$, initial guess $(\mathbf{x}_0, \boldsymbol{\lambda}_0)$, tolerance $\varepsilon$

**Output:** Approximate solution $(\mathbf{x}^*, \boldsymbol{\lambda}^*)$

1. **for** $k = 0, 1, 2, \ldots$:
2. $\qquad$ Evaluate residual $\mathbf{r}_k = \begin{pmatrix} \nabla f - D\mathbf{g}^T\boldsymbol{\lambda}_k \\ \mathbf{g}(\mathbf{x}_k) \end{pmatrix}$
3. $\qquad$ **if** $\|\mathbf{r}_k\| < \varepsilon$: **return** $(\mathbf{x}_k, \boldsymbol{\lambda}_k)$
4. $\qquad$ Form Jacobian $\begin{pmatrix} H_k & -A_k^T \\ A_k & 0 \end{pmatrix}$ where $H_k = \nabla^2_{xx}\mathcal{L}$, $A_k = D\mathbf{g}$
5. $\qquad$ Solve $\begin{pmatrix} H_k & -A_k^T \\ A_k & 0 \end{pmatrix} \begin{pmatrix} \Delta\mathbf{x} \\ \Delta\boldsymbol{\lambda} \end{pmatrix} = -\mathbf{r}_k$
6. $\qquad$ $(\mathbf{x}_{k+1}, \boldsymbol{\lambda}_{k+1}) \gets (\mathbf{x}_k + \Delta\mathbf{x}, \boldsymbol{\lambda}_k + \Delta\boldsymbol{\lambda})$
:::

### The KKT Matrix Structure

The matrix in step 5 is called the **KKT matrix**:

$$
K = \begin{pmatrix} H & -A^T \\ A & 0 \end{pmatrix}
$$

This is a **saddle point system**—it's symmetric but indefinite (has both positive and negative eigenvalues). Special solvers exploit this structure.

:::::{admonition} Properties of the KKT Matrix
:class: note

The KKT matrix is nonsingular if:
1. $H$ is positive definite on the null space of $A$ (sufficient for local minimum)
2. $A$ has full row rank (constraints are linearly independent)

::::{dropdown} Why Indefinite?
Consider the $2 \times 2$ case with $H = 1$, $A = 1$:
$$
K = \begin{pmatrix} 1 & -1 \\ 1 & 0 \end{pmatrix}
$$
Eigenvalues: $\frac{1 \pm \sqrt{5}}{2}$ — one positive, one negative.
::::
:::::

## Inequality Constraints

Real optimization often includes inequality constraints:

$$
\min_{\mathbf{x}} f(\mathbf{x}) \quad \text{subject to} \quad \mathbf{g}(\mathbf{x}) = \mathbf{0}, \quad \mathbf{h}(\mathbf{x}) \leq \mathbf{0}
$$

The KKT conditions become more complex:

::::::{prf:theorem} KKT Conditions (Inequality Constraints)
:label: thm-kkt-conditions-inequality

At a constrained minimum, there exist multipliers $\boldsymbol{\lambda}$, $\boldsymbol{\mu}$ such that:

1. **Stationarity:** $\nabla f = D\mathbf{g}^T\boldsymbol{\lambda} + D\mathbf{h}^T\boldsymbol{\mu}$
2. **Primal feasibility:** $\mathbf{g}(\mathbf{x}) = \mathbf{0}$, $\mathbf{h}(\mathbf{x}) \leq \mathbf{0}$
3. **Dual feasibility:** $\boldsymbol{\mu} \geq \mathbf{0}$
4. **Complementary slackness:** $\mu_i h_i(\mathbf{x}) = 0$ for all $i$

::::{dropdown} Complementary Slackness
The condition $\mu_i h_i = 0$ means:
- If constraint $i$ is inactive ($h_i < 0$), then $\mu_i = 0$ (constraint doesn't affect optimum)
- If $\mu_i > 0$, then $h_i = 0$ (constraint is active/binding)

This is the tricky part—we don't know which constraints are active beforehand!
::::
::::::

### Active Set Methods

One approach: guess which constraints are active, solve the equality-constrained problem, check if the guess was correct.

:::{prf:algorithm} Active Set Method (Sketch)
:label: alg-active-set-method

1. Start with a guess of active constraints $\mathcal{A}$
2. Solve equality-constrained problem with $h_i = 0$ for $i \in \mathcal{A}$
3. Check KKT conditions:
   - If $\mu_i < 0$ for some $i \in \mathcal{A}$: remove $i$ from $\mathcal{A}$
   - If $h_j > 0$ for some $j \notin \mathcal{A}$: add $j$ to $\mathcal{A}$
4. Repeat until convergence
:::

### Interior Point Methods

A more modern approach: replace the inequality constraints with a barrier function:

$$
\min_{\mathbf{x}} f(\mathbf{x}) - \mu \sum_i \log(-h_i(\mathbf{x}))
$$

As $\mu \to 0$, solutions approach the constrained optimum. This leads to smooth systems solvable by Newton's method.

## Alternative Approaches

### Penalty Methods

Replace the constrained problem with a sequence of unconstrained problems:

$$
\min_{\mathbf{x}} f(\mathbf{x}) + \frac{\rho}{2}\|\mathbf{g}(\mathbf{x})\|^2
$$

As $\rho \to \infty$, solutions approach feasibility.

**Problem:** Large $\rho$ makes the problem ill-conditioned.

### Augmented Lagrangian

Combine Lagrange multipliers with a penalty:

$$
\mathcal{L}_\rho(\mathbf{x}, \boldsymbol{\lambda}) = f(\mathbf{x}) - \boldsymbol{\lambda}^T\mathbf{g}(\mathbf{x}) + \frac{\rho}{2}\|\mathbf{g}(\mathbf{x})\|^2
$$

Iterate:
1. Minimize $\mathcal{L}_\rho$ with respect to $\mathbf{x}$ (unconstrained)
2. Update $\boldsymbol{\lambda} \gets \boldsymbol{\lambda} - \rho\mathbf{g}(\mathbf{x})$

**Advantage:** Achieves exact feasibility without $\rho \to \infty$.

## The Big Picture

| Problem | Necessary Conditions | System Structure |
|---------|---------------------|------------------|
| Unconstrained: $\min f$ | $\nabla f = 0$ | $n$ equations |
| Root-finding: $\mathbf{F} = 0$ | Given directly | $n$ equations |
| Constrained: $\min f$ s.t. $\mathbf{g} = 0$ | KKT conditions | $(n+m)$ equations |
| Optimization with $\nabla^2 f$ | Newton on $\nabla f = 0$ | Uses Hessian |
| Constrained optimization | Newton on KKT | Uses Hessian + Jacobian |

Everything reduces to **solving a nonlinear system**. The techniques we've developed—Newton's method, globalization, quasi-Newton—apply directly.

## Example: Constrained Least Squares

:::::{admonition} Example: Regression with Constraints
:class: note

Fit a line $y = ax + b$ to data, but require $a + b = 1$ (perhaps a normalization constraint).

**Objective:** $f(a, b) = \sum_i (y_i - ax_i - b)^2$

**Constraint:** $g(a, b) = a + b - 1 = 0$

The Lagrangian:
$$
\mathcal{L}(a, b, \lambda) = \sum_i (y_i - ax_i - b)^2 - \lambda(a + b - 1)
$$

KKT conditions:
$$
\frac{\partial \mathcal{L}}{\partial a} = -2\sum_i x_i(y_i - ax_i - b) - \lambda = 0
$$
$$
\frac{\partial \mathcal{L}}{\partial b} = -2\sum_i (y_i - ax_i - b) - \lambda = 0
$$
$$
a + b - 1 = 0
$$

This is a linear system (since $f$ is quadratic and $g$ is linear)!
:::::

## The Bigger Picture: From KKT to Optimal Control

The Lagrange multiplier framework extends far beyond finite-dimensional optimization. The same ideas—with the multiplier becoming an **adjoint** or **costate** variable—appear throughout applied mathematics.

### Calculus of Variations

Instead of minimizing $f(\mathbf{x})$, consider minimizing a **functional**:

$$
J[y] = \int_a^b L(x, y(x), y'(x)) \, dx
$$

over functions $y(x)$ with fixed endpoints $y(a) = \alpha$, $y(b) = \beta$.

::::::{prf:theorem} Euler-Lagrange Equation
:label: thm-euler-lagrange

A necessary condition for $y^*$ to minimize $J$ is:

$$
\frac{\partial L}{\partial y} - \frac{d}{dx}\frac{\partial L}{\partial y'} = 0
$$

::::{dropdown} Connection to Lagrange Multipliers
Think of the constraint "$y$ passes through all intermediate points" as infinitely many constraints. The Euler-Lagrange equation is the continuous analog of $\nabla f = \lambda \nabla g$.

More formally: discretize the integral and take the limit. The finite-dimensional KKT conditions become the Euler-Lagrange equation.
::::
::::::

### Optimal Control

Now add dynamics: minimize a cost while a system evolves according to an ODE.

$$
\min_{u(t)} \int_0^T L(\mathbf{x}(t), u(t)) \, dt + \Phi(\mathbf{x}(T))
$$

subject to the **state equation**:

$$
\dot{\mathbf{x}}(t) = \mathbf{f}(\mathbf{x}(t), u(t)), \qquad \mathbf{x}(0) = \mathbf{x}_0
$$

Here $\mathbf{x}(t)$ is the state, $u(t)$ is the control, and the ODE is a constraint!

::::::{prf:theorem} Pontryagin's Maximum Principle
:label: thm-pontryagin-maximum-principle

Define the **Hamiltonian**:

$$
H(\mathbf{x}, u, \boldsymbol{\lambda}) = L(\mathbf{x}, u) + \boldsymbol{\lambda}^T \mathbf{f}(\mathbf{x}, u)
$$

The optimal control satisfies:

1. **State equation:** $\dot{\mathbf{x}} = \nabla_{\boldsymbol{\lambda}} H = \mathbf{f}(\mathbf{x}, u)$ (forward in time)
2. **Adjoint equation:** $\dot{\boldsymbol{\lambda}} = -\nabla_{\mathbf{x}} H$ (backward in time)
3. **Optimality:** $\nabla_u H = 0$ (or $H$ minimized over $u$)

::::{dropdown} The Adjoint = Lagrange Multiplier
The adjoint variable $\boldsymbol{\lambda}(t)$ is exactly a continuous-time Lagrange multiplier enforcing the ODE constraint. Compare:

| Finite-dimensional | Optimal Control |
|-------------------|-----------------|
| $\nabla f = D\mathbf{g}^T\boldsymbol{\lambda}$ | $\dot{\boldsymbol{\lambda}} = -\nabla_{\mathbf{x}} H$ |
| $\mathbf{g}(\mathbf{x}) = \mathbf{0}$ | $\dot{\mathbf{x}} = \mathbf{f}(\mathbf{x}, u)$ |
| Solve simultaneously | Solve forward-backward |
::::
::::::

### Neural ODEs and the Adjoint Method

Neural ODEs model a neural network as a continuous dynamical system:

$$
\frac{d\mathbf{h}}{dt} = f_\theta(\mathbf{h}(t), t)
$$

where $\theta$ are learnable parameters. To train, we need gradients of a loss $L(\mathbf{h}(T))$ with respect to $\theta$.

:::::{admonition} The Adjoint Sensitivity Method
:class: note

Instead of backpropagating through ODE solver steps (memory-intensive), use the adjoint:

$$
\frac{d\mathbf{a}}{dt} = -\mathbf{a}^T \frac{\partial f_\theta}{\partial \mathbf{h}}, \qquad \mathbf{a}(T) = \frac{\partial L}{\partial \mathbf{h}(T)}
$$

where $\mathbf{a}(t) = \partial L / \partial \mathbf{h}(t)$ is the adjoint state.

The parameter gradient is:

$$
\frac{dL}{d\theta} = -\int_T^0 \mathbf{a}(t)^T \frac{\partial f_\theta}{\partial \theta} \, dt
$$

::::{dropdown} Why This Works
This is Pontryagin's principle applied to the "control problem" of choosing $\theta$ to minimize loss! The adjoint equation propagates sensitivity backward through time—exactly like backpropagation, but continuous.

**Key insight:** Backprop through an ODE is just the adjoint method. The "Lagrange multiplier" is the gradient flowing backward.
::::
:::::

### The Unifying Pattern

| Setting | Constraint | Multiplier | Optimality Condition |
|---------|-----------|------------|---------------------|
| Finite-dim optimization | $\mathbf{g}(\mathbf{x}) = 0$ | $\boldsymbol{\lambda}$ | KKT system |
| Calculus of variations | $y$ smooth | (implicit) | Euler-Lagrange |
| Optimal control | $\dot{\mathbf{x}} = \mathbf{f}$ | $\boldsymbol{\lambda}(t)$ (adjoint) | Pontryagin |
| Neural ODEs | $\dot{\mathbf{h}} = f_\theta$ | $\mathbf{a}(t)$ (adjoint) | Adjoint sensitivity |
| PDE-constrained opt | PDE constraint | Adjoint PDE | Lagrangian + adjoints |

**The pattern:** Introduce multipliers to enforce constraints, derive necessary conditions, solve the resulting system.

### Numerical Methods

All of these lead to systems that must be solved numerically:

- **Finite-dim KKT:** Newton's method (as we've seen)
- **Optimal control:** Shooting methods, collocation, direct transcription
- **Neural ODEs:** ODE solvers + adjoint ODE solvers
- **PDE-constrained:** Discretize-then-optimize or optimize-then-discretize

The globalization and quasi-Newton ideas from this chapter apply throughout!

:::::{admonition} Direct Transcription
:class: note

A powerful approach: **discretize** the optimal control problem to get a finite-dimensional NLP, then apply standard optimization.

Discretize $\mathbf{x}(t)$ at times $t_0, t_1, \ldots, t_N$:
- Variables: $\mathbf{x}_0, \mathbf{x}_1, \ldots, \mathbf{x}_N, u_0, u_1, \ldots, u_{N-1}$
- Constraints: $\mathbf{x}_{k+1} = \mathbf{x}_k + h\mathbf{f}(\mathbf{x}_k, u_k)$ (or better integrators)
- Objective: $\sum_k L(\mathbf{x}_k, u_k) + \Phi(\mathbf{x}_N)$

Now it's just KKT! Newton's method on the discretized system.
:::::

## Summary

**Key insight:** Constrained optimization = nonlinear system via Lagrange multipliers.

**Advantages of Newton on KKT:**
- Quadratic convergence (locally)
- Handles equality constraints directly
- Natural framework for inequality constraints (with extensions)

**Challenges:**
- KKT matrix is indefinite (saddle point structure)
- Inequality constraints require active set or interior point methods
- Still need good initial guess

**The unified view:** Whether you're finding roots, minimizing functions, optimizing with constraints, solving optimal control problems, or training neural ODEs—it all comes down to:

1. Formulate necessary conditions (introduce Lagrange multipliers / adjoints)
2. Solve the resulting system with Newton-type methods
3. Use globalization to ensure convergence from reasonable initial guesses

The techniques in this chapter—Newton's method, Jacobian computation, globalization, quasi-Newton updates—are the foundation for computational methods across applied mathematics, from engineering optimization to machine learning.
