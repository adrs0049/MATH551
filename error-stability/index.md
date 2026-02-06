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

This chapter builds up the ideas needed to explain the fast inverse square root:

1. **Floating-point representation** — How computers store numbers using IEEE 754. Every number has a small representation error bounded by machine epsilon. We'll see how subtracting nearly equal numbers destroys accuracy, and derive the optimal step size for finite differences.

2. **Condition numbers** — Some mathematical problems amplify input errors. The condition number $\kappa$ measures this sensitivity. We'll see that subtraction is ill-conditioned when operands are close — giving a deeper explanation of the finite difference trade-off.

3. **Fast inverse square root** — Now we can explain the bit manipulation (it exploits the floating-point representation to approximate a logarithm) and the Newton refinement step (which bridges naturally to root-finding methods in the [next chapter](../nonlinear-equations/index.md)).

## Learning Outcomes

After completing this chapter, you should be able to:

- **L2.1:** Explain IEEE 754 floating-point representation.
- **L2.2:** Define machine epsilon and its significance.
- **L2.3:** Compute condition numbers for simple functions.
- **L2.4:** Identify sources of numerical instability.
- **L2.5:** Reformulate expressions to avoid cancellation.
