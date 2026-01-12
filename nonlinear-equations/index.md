# Nonlinear Equations

:::{tip} Big Idea
This chapter is about **iterative methods**: algorithms that generate sequences of approximations converging to a solution. Newton's method—which you know from calculus—is the workhorse, but understanding *when* and *why* it works requires careful analysis.
:::

## Motivation: The Quake Fast Inverse Square Root

Before diving into theory, consider this legendary piece of code from Quake III Arena (1999):

```c
float Q_rsqrt(float x) {
    long i;
    float x2, y;
    x2 = x * 0.5F;
    y  = x;
    i  = *(long*)&y;                    // interpret float bits as integer
    i  = 0x5f3759df - (i >> 1);         // magic!
    y  = *(float*)&i;
    y  = y * (1.5F - (x2 * y * y));     // Newton refinement
    return y;
}
```

This computes $y = 1/\sqrt{x}$—a function that cannot be evaluated directly using just $+, -, \times, \div$. How do we turn this into something computable?

### Translating to a Root-Finding Problem

::::::{prf:remark} Key Insight
:label: rmk-inverse-sqrt-root-finding

Finding $y = 1/\sqrt{x}$ is equivalent to solving $g(y) = 0$ where $g(y) = 1/y^2 - x$.

:::{dropdown} Why?
If $y = 1/\sqrt{x}$, then $y^2 = 1/x$, so $1/y^2 = x$, which means $1/y^2 - x = 0$.
:::
::::::

**Why translate to root finding?** Because you already know an algorithm for finding roots: **Newton's method** from calculus.

### Newton's Method Refresher

Recall: to solve $g(y) = 0$, approximate $g$ by its tangent line and find where it crosses zero:

$$
y_{n+1} = y_n - \frac{g(y_n)}{g'(y_n)}
$$

::::{dropdown} Applying to $g(y) = 1/y^2 - x$
For $g(y) = 1/y^2 - x$, we have $g'(y) = -2/y^3$. Substituting:

$$
y_{n+1} = y_n - \frac{1/y_n^2 - x}{-2/y_n^3} = y_n + \frac{y_n}{2}\left(1 - xy_n^2\right) = y_n\left(\frac{3}{2} - \frac{x}{2}y_n^2\right)
$$

This is exactly the last line of the Quake code!
::::

The [floating-point bit tricks](../error-stability/fast-inverse-sqrt.md) provide a clever initial guess, and one Newton iteration refines it to sufficient accuracy.

## The Root-Finding Problem

Given a nonlinear function $f(x)$, find $x^*$ such that $f(x^*) = 0$.

**Examples:**
- $f(x) = e^x - 2$ has root $x^* = \ln 2$
- $f(x) = x^2 - 2$ has roots $x^* = \pm\sqrt{2}$
- $f(x) = \cos(x) - x$ has a root near $x^* \approx 0.739$

For most nonlinear equations, we cannot find closed-form solutions. We need a different approach.

## The Big Picture: Iterative Methods

The Quake example illustrates a fundamental pattern in scientific computing:

:::{prf:definition} Iterative Methods
:label: def-iterative-methods

Many problems cannot be solved in a finite number of arithmetic operations. Instead, we construct a **sequence** $x_0, x_1, x_2, \ldots$ that converges to the solution:

$$
\lim_{n \to \infty} x_n = x^*
$$

In practice, we stop when $|x_{n+1} - x_n|$ or some other criterion is sufficiently small.
:::

The key questions for any iterative method are:
1. **Does it converge?** Under what conditions?
2. **How fast?** How many iterations do we need?
3. **How sensitive is the answer?** (This is where condition numbers come in.)

### A Unifying Principle

Newton's method, bisection, and fixed-point iteration may seem like three unrelated techniques. But they share a common structure: each generates a sequence that (hopefully) converges to the root.

The key insight is that convergence happens when the iteration *shrinks* errors:

:::{admonition} Contraction → Convergence
:class: tip

If each step of an iteration reduces the error by a factor $L < 1$:

$$
|x_{n+1} - x^*| \leq L \, |x_n - x^*|
$$

then the sequence converges geometrically to $x^*$.
:::

This principle—formalized as the **Banach Fixed Point Theorem** in [Fixed Point Iteration](fixed-point.md)—is one of the most important ideas in numerical analysis. It explains:
- Why Newton's method converges (when it does)
- Why bisection always converges (it halves the interval each step)
- Why iterative linear solvers like Jacobi and Gauss-Seidel work

The same idea recurs throughout this course: in ODE integrators, optimization algorithms, and beyond.

## Learning Outcomes

After completing this chapter, you should be able to:

- **L3.1:** Derive Newton's method from Taylor series.
- **L3.2:** Implement Newton's method in code.
- **L3.3:** Explain quadratic convergence.
- **L3.4:** Identify when Newton's method fails.
- **L3.5:** State the Intermediate Value Theorem.
- **L3.6:** Derive the bisection algorithm.
- **L3.7:** Calculate required bisection iterations.
- **L3.8:** Define fixed points and restate root-finding.
- **L3.9:** State the Banach Fixed Point Theorem.
- **L3.10:** Apply convergence conditions.
- **L3.11:** Compute condition numbers for roots.
