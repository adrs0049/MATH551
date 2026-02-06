# Forward and Backward Error

:::{tip} Big Idea
Forward error asks: "How wrong is my answer?" Backward error asks: "What problem did I actually solve?" A **stable algorithm** solves a nearby problem exactly—and the condition number determines whether "nearby problem" means "nearby answer."
:::

## A Motivating Example

Consider computing $f(x) = \sqrt{x+1} - \sqrt{x}$ for large $x$.

The condition number is $\kappa \leq 1/2$ — this is a **well-conditioned problem**!

Yet using 4-digit decimal arithmetic:

$$
f(1984) = \sqrt{1985} - \sqrt{1984} = 44.55 - 44.54 = 0.01
$$

The true value is $0.01122$. We lost two significant digits — **10% relative error** from a well-conditioned problem!

**What went wrong?** The algorithm subtracts two nearly-equal numbers, amplifying their individual roundoff errors. This is called **subtractive cancellation**.

**A stable alternative:** Rationalize the expression:

$$
\sqrt{x+1} - \sqrt{x} = \frac{(\sqrt{x+1} - \sqrt{x})(\sqrt{x+1} + \sqrt{x})}{\sqrt{x+1} + \sqrt{x}} = \frac{1}{\sqrt{x+1} + \sqrt{x}} = \frac{1}{89.09} = 0.01122 \quad \checkmark
$$

Same problem, different algorithm, full accuracy.

This example reveals something fundamental: **not all errors are the problem's fault**. We need a framework to separate algorithm quality from problem sensitivity.

## Forward Error

Consider computing $y = f(x)$. Due to roundoff or approximation, we actually compute $\tilde{y}$.

:::{prf:definition} Forward Error
:label: def-forward-error

$$
\text{Forward error} = |\tilde{y} - y| = |\tilde{y} - f(x)|
$$

**Question answered:** How far is our answer from the truth?
:::

Forward error is the most intuitive measure—it directly answers "how wrong is my answer?"

### Absolute vs. Relative Error

For meaningful comparisons, we often use **relative error**:

:::{prf:definition} Relative Forward Error
:label: def-relative-forward-error

$$
\text{Relative forward error} = \frac{|\tilde{y} - y|}{|y|}
$$

This measures error as a fraction of the true value.
:::

**Example:** An error of 0.001 means different things for:
- $y = 1$: relative error is 0.1% (excellent)
- $y = 0.001$: relative error is 100% (useless)

Relative error captures this distinction. When we say "accurate to 6 digits," we mean relative error $\approx 10^{-6}$.

:::{prf:definition} Forward Stable Algorithm
:label: def-forward-stable

An algorithm is **forward stable** if:

$$
\frac{|\tilde{y} - y|}{|y|} = O(\varepsilon_{\text{mach}} \cdot \kappa)
$$

The forward error is bounded by machine precision times the condition number.
:::

### The Catch

Forward error seems like the natural measure—we want to know how wrong we are! But there's a problem:

**Forward error requires knowing the true answer $y = f(x)$.**

If we knew the true answer, we wouldn't need to compute it. This is where backward error comes in.

## Backward Error

Instead of asking "how wrong is my answer?", backward error asks a different question.

:::{prf:definition} Backward Error
:label: def-backward-error

$$
\text{Backward error} = \min\{|\tilde{x} - x| : f(\tilde{x}) = \tilde{y}\}
$$

**Question answered:** What input would give our output exactly?
:::

In other words: our computed $\tilde{y}$ might not equal $f(x)$, but it *does* equal $f(\tilde{x})$ for some nearby $\tilde{x}$. The backward error measures how far $\tilde{x}$ is from $x$.

### Why Backward Error?

At first, backward error seems less useful than forward error. But it turns out to be more fundamental for algorithm analysis. Here's why:

**1. Backward error is computable**

Forward error requires knowing the true answer—but if we knew that, we wouldn't need to compute it!

Backward error, on the other hand, is often **directly computable**:
- **Linear systems:** Given computed $\tilde{x}$, the residual $r = b - A\tilde{x}$ tells us backward error immediately—we solved $A\tilde{x} = b - r$ exactly.
- **Root finding:** Given approximate root $\tilde{r}$, the residual $f(\tilde{r})$ measures backward error.

**2. Backward error separates algorithm from problem**

There are two sources of error:
- **Algorithm quality:** Does the algorithm solve *some* nearby problem accurately?
- **Problem sensitivity:** How much does the answer change when the problem changes slightly?

Backward error isolates the first question. If an algorithm has small backward error, it's doing its job well—any remaining forward error is the problem's fault (ill-conditioning), not the algorithm's.

**3. Backward error accounts for input uncertainty**

Real-world inputs have measurement error. If your input $x$ has uncertainty $\delta$, then:
- You don't actually know you're solving $f(x)$—you might be solving $f(x + \epsilon)$ for some $|\epsilon| \leq \delta$
- An algorithm with backward error $\leq \delta$ is **as good as exact** for your purposes

**Principle:** Don't ask for more accuracy than your inputs justify.

### Residual as Backward Error

For many problems, backward error equals the **residual**—a directly computable quantity.

| Problem | Residual | Interpretation |
|---------|----------|----------------|
| Linear system $Ax = b$ | $r = b - A\tilde{x}$ | $\tilde{x}$ solves $Ax = b - r$ exactly |
| Root finding $f(x) = 0$ | $f(\tilde{r})$ | $\tilde{r}$ is a root of $f(x) - f(\tilde{r})$ |
| Least squares | $\|b - A\tilde{x}\|$ | How well $\tilde{x}$ fits the data |

:::{prf:definition} Backward Stable Algorithm
:label: def-backward-stable

An algorithm is **backward stable** if for all inputs $x$, the computed output $\tilde{y}$ satisfies:

$$
\tilde{y} = f(\tilde{x}) \quad \text{for some } \tilde{x} \text{ with } \frac{|\tilde{x} - x|}{|x|} = O(\varepsilon_{\text{mach}})
$$

The computed answer is the exact answer to a slightly perturbed problem.
:::

**Key insight:** Backward stability implies forward stability (for well-conditioned problems). This is why we focus on backward stability—it's the stronger, more useful guarantee.

(golden-rule)=
## The Relationship: Forward Error ≤ κ × Backward Error

Now we connect backward error to forward error using the **condition number** from the [error-stability chapter](../error-stability/condition-numbers.md).

### Derivation

Suppose we compute $\tilde{y}$ instead of the true $y = f(x)$. If the algorithm is backward stable:

$$
\tilde{y} = f(\tilde{x}) \quad \text{for some } \tilde{x} \text{ near } x
$$

The forward error is the difference $|\tilde{y} - y| = |f(\tilde{x}) - f(x)|$.

By Taylor expansion:

$$
|f(\tilde{x}) - f(x)| \approx |f'(x)| \cdot |\tilde{x} - x|
$$

Converting to relative errors:

$$
\frac{|\tilde{y} - y|}{|y|} \approx \frac{|f'(x)|}{|f(x)|} \cdot |\tilde{x} - x| = \underbrace{\left|\frac{x f'(x)}{f(x)}\right|}_{\kappa} \cdot \frac{|\tilde{x} - x|}{|x|}
$$

::::{prf:theorem} The Golden Rule
:label: thm-golden-rule

$$
\underbrace{\frac{|\tilde{y} - y|}{|y|}}_{\text{relative forward error}} \lesssim \underbrace{\kappa}_{\text{condition number}} \cdot \underbrace{\frac{|\tilde{x} - x|}{|x|}}_{\text{relative backward error}}
$$

**Forward error ≤ Condition number × Backward error**

This is the central equation of numerical analysis. It cleanly separates:
- **Problem sensitivity** (κ) — intrinsic to the mathematics
- **Algorithm quality** (backward error) — what we can control
::::

### What This Means

| Component | Controlled by | Can we improve it? |
|-----------|--------------|-------------------|
| Forward error | Both | Indirectly |
| Condition number $\kappa$ | The problem | No (it's math) |
| Backward error | The algorithm | Yes! |

**Practical implications:**
- **Well-conditioned + backward stable → accurate:** Small $\kappa$ and small backward error guarantee small forward error
- **Ill-conditioned → trouble:** Even perfect algorithms give poor forward error when $\kappa$ is large
- **Focus on backward stability:** Write algorithms with small backward error; that's all you can control

## Unstable Algorithms: Common Pitfalls

Even well-conditioned problems can produce terrible results with unstable algorithms. Three operations are particularly dangerous:

| Operation | Problem | Stable Alternative |
|-----------|---------|-------------------|
| Subtracting close numbers | Relative error explodes | Reformulate algebraically |
| Dividing by small numbers | Amplifies numerator error | Use logarithms or rescale |
| Adding numbers of vastly different magnitude | Small values lost | Sort and sum smallest first |

::::{dropdown} Example: The Quadratic Formula
:color: danger
:icon: alert

Solve $x^2 - 10^6 x + 1 = 0$. The roots are $x_1 \approx 10^6$ and $x_2 \approx 10^{-6}$.

**Standard formula for the small root:**

$$
x_2 = \frac{10^6 - \sqrt{10^{12} - 4}}{2} \approx \frac{10^6 - 10^6}{2} = 0
$$

Catastrophic cancellation! The true value is $10^{-6}$, so we have **100% relative error**.

**Stable alternative:** Use $x_1 x_2 = c/a = 1$, so $x_2 = 1/x_1$. Compute $x_1$ (no cancellation), then $x_2 = 1/x_1$ with full precision.
::::

::::{dropdown} Example: Computing $e^x - 1$
:color: danger
:icon: alert

For small $x$, direct computation fails:

```python
x = 1e-15
result = np.exp(x) - 1  # Returns ~1.1e-15 (only 1 correct digit!)
```

**Why?** $e^{10^{-15}} \approx 1.000000000000001$, stored as $1.0$ due to limited precision. Then $1.0 - 1.0 \approx 0$.

**Stable alternative:** Use the dedicated function:

```python
result = np.expm1(1e-15)  # Returns 1.0e-15 (full precision!)
```
::::

::::{dropdown} Example: Kahan Summation
:color: success
:icon: check-circle

When summing many numbers, roundoff errors accumulate. **Kahan summation** tracks and compensates for these errors:

```python
def kahan_sum(x):
    """Compensated summation — much more accurate than naive sum."""
    s = 0.0  # Running sum
    c = 0.0  # Compensation for lost low-order bits
    for xi in x:
        y = xi - c        # Compensated value to add
        t = s + y         # New sum (may lose precision)
        c = (t - s) - y   # Recover what was lost
        s = t
    return s
```

This achieves $O(\varepsilon_{\text{mach}})$ error instead of $O(n \cdot \varepsilon_{\text{mach}})$ for naive summation.
::::

## Summary

| Concept | Measures | Controlled by | Computable? |
|---------|----------|---------------|-------------|
| Forward error | Distance to true answer | Problem + algorithm | Usually not |
| Backward error | Perturbation to input | Algorithm | Often yes (residual) |
| Condition number | Sensitivity of problem | Problem (math) | Sometimes |

**The workflow:**
1. Check if the problem is well-conditioned (small $\kappa$)
2. Use a backward-stable algorithm (small backward error)
3. Then forward error will be small automatically

**If forward error is large despite backward stability, the problem is ill-conditioned—reformulate if possible.**
