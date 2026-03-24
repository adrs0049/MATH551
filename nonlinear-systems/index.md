# Nonlinear Systems of Equations

:::{admonition} Optional Chapter
:class: warning

This chapter covers material beyond the core MATH 551 syllabus. It extends the
ideas from the [nonlinear equations](../nonlinear-equations/index.md) chapter to
systems of equations in $\mathbb{R}^n$ and introduces globalization strategies
for robust convergence. This material is included for interested students and as
preparation for MATH 552.
:::

:::{tip} Big Idea
Newton's method extends naturally to systems: replace division by $f'(x)$ with solving a linear system involving the Jacobian. But in higher dimensions, **globalization** becomes essential—local quadratic convergence means nothing if you can't get close to the root. This connects nonlinear equation solving to optimization.
:::

## Overview

We seek $\mathbf{x}^* \in \mathbb{R}^n$ such that $\mathbf{F}(\mathbf{x}^*) = \mathbf{0}$, where $\mathbf{F}: \mathbb{R}^n \to \mathbb{R}^n$.

This problem appears everywhere:
- **Engineering:** Equilibrium equations, circuit analysis, structural mechanics
- **Optimization:** Finding critical points where $\nabla f(\mathbf{x}) = \mathbf{0}$
- **Differential equations:** Implicit methods require solving nonlinear systems at each step
- **Machine learning:** Training neural networks (finding zeros of the gradient)

### The Optimization Connection

There's a deep connection between root-finding and optimization. If we want to minimize a function $\phi: \mathbb{R}^n \to \mathbb{R}$, necessary conditions require:

$$
\nabla \phi(\mathbf{x}^*) = \mathbf{0}
$$

This is a system of $n$ nonlinear equations! Newton's method for optimization (finding minima) and Newton's method for root-finding are closely related—the Hessian $\nabla^2 \phi$ plays the role of the Jacobian.

### The Challenge in Higher Dimensions

In one dimension, Newton's method is straightforward: good initial guess → quadratic convergence. In $n$ dimensions, things get harder:

1. **Finding a good initial guess** is much harder in high-dimensional spaces
2. **The Jacobian may be singular** or nearly singular
3. **Computing the Jacobian** requires $n^2$ partial derivatives
4. **Each iteration solves** an $n \times n$ linear system: $\mathcal{O}(n^3)$ cost

These challenges motivate **globalization strategies**: modifications that make Newton's method converge from farther away.

## Learning Outcomes

After completing this chapter, you should be able to:

**Core:**
- **L6.1:** Write systems as $\mathbf{F}(\mathbf{x}) = \mathbf{0}$.
- **L6.2:** Compute the Jacobian matrix.
- **L6.3:** State Newton's method for systems.
- **L6.4:** Perform Newton iterations by hand (2D examples).
- **L6.5:** Implement Newton for systems.
- **L6.6:** Explain the cost per iteration ($\mathcal{O}(n^3)$).
- **L6.7:** Explain why globalization is needed.
- **L6.8:** State the Armijo condition for line search.
- **L6.9:** Implement damped Newton with backtracking.

**Optional:**
- **L6.10:** State local quadratic convergence conditions.
- **L6.11:** Explain the idea behind quasi-Newton methods.
- **L6.12:** Formulate constrained optimization as a KKT system.
- **L6.13:** Explain the role of Lagrange multipliers geometrically.

