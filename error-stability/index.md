# Floating Point and Error Analysis

:::{tip} Big Idea
Numbers on a computer are approximations. Every floating-point computation introduces small errors, and some mathematical operations amplify those errors dramatically. Understanding floating-point representation and condition numbers lets us predict when computations will lose accuracy — and when clever tricks can exploit the representation itself.
:::

## A Teaser: The Fastest Inverse Square Root

Here are 6 lines of C that compute $1/\sqrt{x}$ faster than any library call in the 1990s:

```c
float Q_rsqrt( float x )
{
    long i;
    float x2, y;
    const float threehalfs = 1.5F;

    x2 = x * 0.5F;
    y  = x;
    i  = * ( long * ) &y;               // evil floating point bit level hacking
    i  = 0x5f3759df - ( i >> 1 );       // what the fuck?
    y  = * ( float * ) &i;
    y  = y * ( threehalfs - ( x2 * y * y ) );   // 1st iteration
    return y;
}
```

This code — from the [Quake III Arena source](https://github.com/id-Software/Quake-III-Arena) — powered real-time 3D graphics on 1990s hardware. It looks like black magic, but every line has a precise mathematical explanation.

**By the end of this chapter, you'll understand every line.**

## Overview

This chapter is driven by a mystery from [Chapter 1](../approximation-theory/numerical-differentiation.md): *why does the finite difference error increase when we make the step size smaller?* By the end, you'll not only understand the answer, but also see how floating-point representation can be exploited for creative algorithmic tricks.

1. **Errors and floating-point arithmetic** — We define absolute and relative error, then introduce floating-point numbers as approximations with a relative error guarantee (machine epsilon). This immediately solves the finite difference mystery: subtracting nearly equal numbers amplifies relative error, creating a trade-off between truncation and round-off.

2. **Condition numbers** — Could we have predicted the finite difference problem before running any code? Yes — the condition number $\kappa$ measures how much a computation amplifies input errors. We'll see that subtraction is ill-conditioned when operands are close, giving a deeper theoretical explanation of the trade-off.

3. **Fast inverse square root** — Now we look *inside* the floating-point representation (IEEE 754 bit layout) and discover that reinterpreting float bits as an integer gives an approximate logarithm. This insight, combined with Newton's method, produces one of the most famous algorithms in computer science.

## Learning Outcomes

After completing this chapter, you should be able to:

- **L2.1:** Define absolute and relative error, and explain why relative error is the appropriate measure for floating-point computations.
- **L2.2:** State the floating-point guarantee ($\text{fl}(x) = x(1+\varepsilon)$) and define machine epsilon.
- **L2.3:** Explain the finite difference error trade-off between truncation and round-off.
- **L2.4:** Compute condition numbers for simple functions and identify ill-conditioned problems.
- **L2.5:** Explain how the fast inverse square root exploits IEEE 754 bit representation.
