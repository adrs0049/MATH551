# Approximation Theory

:::{tip} Big Idea
Computers can only do arithmetic. Everything else—$e^x$, $\sin(x)$, derivatives, integrals—must be approximated using additions and multiplications. Taylor's theorem tells us how to do this and how much error we incur.
:::

## Overview

Suppose you want to write a program to compute $e^x$, or solve a differential equation, or evaluate an integral. How would you do it?

The challenge is that computers are simple machines. They can only perform basic arithmetic: $+$, $-$, $\times$, $\div$. But scientific computing demands much more:

- How do you compute $e^x$? You can't—not exactly—using only arithmetic.
- How do you compute $\sin(x)$, $\ln(x)$, $\sqrt{x}$? Same problem.
- How do you compute a derivative $f'(x)$? You'd need a limit—not a finite operation.
- How do you compute an integral $\int_a^b f(x)\,dx$? Again, not directly computable.

**The solution:** Approximate these objects using only arithmetic. Taylor polynomials let us write:

$$
e^x \approx 1 + x + \frac{x^2}{2} + \frac{x^3}{6} + \cdots
$$

Now we *can* compute the right-hand side—it's just additions and multiplications. The cost is an approximation error. Taylor's theorem tells us exactly how large that error is:

$$
f(x) = \underbrace{f(a) + f'(a)(x-a) + \cdots + \frac{f^{(n)}(a)}{n!}(x-a)^n}_{\text{computable}} + \underbrace{R_n(x)}_{\text{error}}
$$

This is the starting point for all of numerical analysis.

## Learning Outcomes

After completing this chapter, you should be able to:

- **L1.1:** Write Taylor expansions with Lagrange remainder.
- **L1.2:** Derive finite difference formulas from Taylor series.
- **L1.3:** Analyze truncation error order ($O(h)$, $O(h^2)$, etc.).
- **L1.4:** Explain the trade-off between truncation and roundoff error.
- **L1.5:** Choose appropriate step sizes for numerical differentiation.


