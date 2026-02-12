# Stable and Unstable Algorithms

:::{tip} Big Idea
The condition number tells us about the **problem**. But the same problem can be solved by different algorithms — and some amplify errors while others don't. A **stable** algorithm solves a nearby problem exactly; an **unstable** algorithm magnifies rounding errors, producing results far worse than the problem's conditioning would suggest.
:::

In the [previous section](condition-numbers.md), we saw that the [condition number](condition-numbers.md#def-condition-number) measures how sensitive a problem is to input perturbations. But sensitivity of the problem is only half the story — the **algorithm** we choose matters just as much.

## Definitions

:::{prf:definition} Stable Algorithm
:label: def-stable-algorithm

A **stable** algorithm for a given problem yields a numerical approximation that is the exact solution of an only slightly perturbed problem. Equivalently, a stable algorithm does not amplify rounding errors beyond what the problem's conditioning requires.

An **unstable** algorithm magnifies rounding errors, producing results far worse than the problem's conditioning would suggest.
:::

:::{prf:remark}
:label: rmk-stability-conditioning-pair

The two key questions in any numerical computation are:

1. **Is the problem well-conditioned?** (condition number — intrinsic to the mathematics)
2. **Is the algorithm stable?** (algorithm design — what we can control)

If both answers are yes, the computed result will be accurate. If the problem is ill-conditioned, no algorithm can help. But if the problem is well-conditioned and the result is poor, the algorithm is to blame.
:::

## Example: Two Algorithms for the Same Function

:::{prf:example} Stable vs. Unstable Evaluation of $\sqrt{x+1} - \sqrt{x}$
:label: ex-stable-unstable-sqrt

Consider $f(x) = \sqrt{x+1} - \sqrt{x}$. The condition number satisfies $\kappa \leq 1/2$ — this is a **well-conditioned problem**.

**Algorithm 1 (unstable):** Evaluate directly using 4-digit decimal arithmetic:

$$
f(1984) = \sqrt{1985} - \sqrt{1984} = 0.4455 \times 10^2 - 0.4454 \times 10^2 = 0.0001 \times 10^2 = 0.1000 \times 10^{-1}
$$

The true value is $0.1122 \times 10^{-1}$. We lost two significant digits — **10% relative error** from a well-conditioned problem! The subtraction of nearly equal numbers is [ill-conditioned](condition-numbers.md#ex-subtraction-condition), with $\kappa = (|a|+|b|)/|a-b| \approx 0.8909 \times 10^2 / 0.0001 \times 10^2 \approx 10^4$, causing catastrophic cancellation.

**Algorithm 2 (stable):** Rationalize the expression:

$$
\sqrt{x+1} - \sqrt{x} = \frac{(\sqrt{x+1} - \sqrt{x})(\sqrt{x+1} + \sqrt{x})}{\sqrt{x+1} + \sqrt{x}} = \frac{1}{\sqrt{x+1} + \sqrt{x}}
$$

Now:
$$
f(1984) = \frac{1}{\sqrt{1985} + \sqrt{1984}} = \frac{1}{0.4455 \times 10^2 + 0.4454 \times 10^2} = \frac{1}{0.8909 \times 10^2} = 0.1122 \times 10^{-1} \quad \checkmark
$$

Same problem, different algorithm, full accuracy. Algorithm 1 is **unstable** because it introduces an ill-conditioned subtraction step; Algorithm 2 avoids the subtraction and is **stable**.
:::

## Example: The Finite Difference Revisited

The same pattern — well-conditioned problem, unstable algorithm — appears in a computation we already know well.

:::{prf:example} Computing $f'(x)$ — Well-Conditioned Problem, Unstable Algorithm
:label: ex-fd-stability

Consider computing $f'(x)$ for $f(x) = \sin(x)$ at $x = 1$. The problem's condition number is:

$$
\kappa_{\text{problem}} = \left|\frac{x\,f''(x)}{f'(x)}\right| = \left|\frac{-\sin(1)}{\cos(1)}\right| = |\tan(1)| \approx 1.56
$$

The *problem* is well-conditioned. But the forward difference algorithm computes $f'(x)$ via the subtraction $f(x+h) - f(x)$, whose [condition number](condition-numbers.md#ex-subtraction-condition) is:

$$
\kappa_{\text{subtraction}} \approx \frac{2|f(x)|}{h|f'(x)|} \sim \frac{1}{h} \to \infty \quad \text{as } h \to 0
$$

The algorithm introduces an ill-conditioned subtraction step that the problem itself does not require. A different method (e.g., automatic differentiation) avoids this subtraction entirely and is stable.
:::

Recognizing this pattern — *is the problem sensitive, or is the algorithm choosing a sensitive path?* — is one of the central skills in numerical analysis.

## Floating-Point Operations to Watch For

Not every source of error is an ill-conditioning issue. The three operations below can each cause trouble, but for **different reasons**:

1. **Subtracting nearly equal numbers** ($x - y$ when $x \approx y$): Subtraction of nearly equal numbers is genuinely **ill-conditioned** — the [condition number of subtraction](condition-numbers.md#ex-subtraction-condition) is $\kappa = (|a|+|b|)/|a-b| \to \infty$ as $a \to b$. The key insight is that an algorithm may *choose* to subtract nearly equal numbers as an intermediate step, even when the overall problem doesn't require it. If the subtraction can be avoided by algebraic reformulation, the algorithm becomes stable.

2. **Dividing by a small number** ($x / y$ when $|y| \ll 1$): This is *not* ill-conditioned — the condition number of $f(y) = x/y$ with respect to $y$ is $\kappa = 1$. The issue is that small absolute errors in $y$ become large absolute errors in the result. If $y$ itself was computed with some absolute error $\delta$, then $x/(y+\delta) - x/y \approx -x\delta/y^2$, which is large when $y$ is small. **Remedy:** rewrite using logarithms or rescale the computation.

3. **Adding numbers of vastly different magnitude** ($x + y$ when $|x| \gg |y|$): This is also *not* an ill-conditioning issue. Subtraction is not involved, so there is no cancellation. The problem is that floating-point numbers have finite precision: if $y$ is smaller than the spacing between consecutive floats near $x$, then $\text{fl}(x + y) = x$ and $y$ is simply lost. **Remedy:** sort and sum smallest first, or use compensated summation (Kahan summation).

## Looking Ahead

We will formalize these ideas using **forward and backward error** when we study [linear systems](../qr-least-squares/forward-backward-error.md). There, we will see that backward stability — the precise version of "solves a nearby problem exactly" — is the gold standard for numerical algorithms. For now, the key takeaway is: **distinguish between the problem (conditioning) and the algorithm (stability)**.
