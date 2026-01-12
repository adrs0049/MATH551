# Error and Stability

:::{tip} Big Idea
Numbers on a computer are approximations. Forward error measures how wrong our
answer is; backward error measures how much we'd have to perturb the input to
make our answer exact. A stable algorithm has small backward error.
:::

## Overview

Every floating-point computation introduces error:
- **Representation error:** Most real numbers can't be stored exactly
- **Rounding error:** Arithmetic operations round their results
- **Accumulation:** Errors compound through sequences of operations

Understanding these errors—and designing algorithms that control them—is essential for reliable scientific computing.
This chapter builds up the error analysis framework in four steps:

1. **Floating-point representation** — Computers can't store most real numbers exactly. Every number has a small representation error bounded by machine epsilon.

2. **The Quake fast inverse square root** — A famous algorithm showing how deep understanding of floating-point enables creative numerical tricks. It bridges to Newton's method in the next chapter.

3. **Condition numbers** — Some mathematical problems amplify errors. The condition number $\kappa$ measures this sensitivity. This is a property of the **problem**, not the algorithm.

4. **Forward and backward error** — Algorithms introduce additional errors. Backward error measures algorithm quality; forward error is what we actually care about.

## The Golden Rule

$$
\boxed{\text{Forward error} \leq \text{Condition number} \times \text{Backward error}}
$$

This cleanly separates:
- **Problem sensitivity** ($\kappa$) — intrinsic to the mathematics, we can't change it
- **Algorithm quality** (backward error) — what we can control through careful implementation

A **stable algorithm** produces answers with small backward error.
An **ill-conditioned problem** amplifies any error, no matter how good the algorithm.

## Learning Outcomes

After completing this chapter, you should be able to:

- **L2.1:** Explain IEEE 754 floating-point representation.
- **L2.2:** Define machine epsilon and its significance.
- **L2.3:** Distinguish forward and backward error.
- **L2.4:** Compute condition numbers for simple functions.
- **L2.5:** Identify sources of numerical instability.
- **L2.6:** Reformulate expressions to avoid cancellation.

