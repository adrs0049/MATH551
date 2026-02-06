# Exercises

Practice analyzing errors, stability, and conditioning.

## Self-Assessment Questions

Test your understanding with these conceptual questions:

1. **Machine Epsilon:** What happens when you add a number smaller than $\varepsilon_{\text{mach}}$ to $1.0$?

2. **Conditioning:** What makes evaluating $\tan(x)$ near $x = \pi/2$ ill-conditioned? Can any algorithm fix this?

3. **Stability:** Why is computing $\sqrt{x+1} - \sqrt{x}$ directly unstable for large $x$? What's a stable alternative?

4. **Machine Epsilon:** What is the approximate machine epsilon for single and double precision? What does it tell us about the accuracy of floating point computations?

5. **Binary Representation:** How do you convert a positive integer to its negative two's complement representation?

6. **Subtractive Cancellation:** When subtracting two numbers $a$ and $b$, under what conditions is the result reliable?

7. **Trade-offs:** In finite differences, why can't we just use an extremely small $h$ to get arbitrary accuracy?


## Floating-Point Arithmetic

---

### Q2.1: Machine Epsilon

- **(a)** Write a program to experimentally determine machine epsilon for single (32-bit) and double (64-bit) precision.
- **(b)** Verify your answers against the theoretical values $2^{-24}$ and $2^{-53}$.
- **(c)** Compute $1 + \varepsilon_{\text{mach}}/2$ in double precision. What do you get? Why?

---

### Q2.2: Representation Errors

- **(a)** Show that $0.1$ cannot be represented exactly in binary floating-point.
- **(b)** Compute `0.1 + 0.1 + 0.1 - 0.3` in Python. Explain the result.
- **(c)** Why might `sum([0.1] * 10)` not equal `1.0`?

---

### Q2.3: Summation Order

Consider summing $s = \sum_{k=1}^{n} \frac{1}{k}$ for $n = 10^7$.

- **(a)** Compute the sum from $k = 1$ to $n$ (forward).
- **(b)** Compute the sum from $k = n$ to $1$ (backward).
- **(c)** Which is more accurate? Why? (Hint: think about the relative sizes of terms being added.)

---

### Q2.4: Kahan Summation

Kahan (1965) proposed a **compensated summation** algorithm that tracks roundoff error:

```python
def kahan_sum(xs):
    s = 0.0  # Running sum
    c = 0.0  # Compensation for lost low-order bits

    for x in xs:
        y = x - c        # Compensate for previous error
        t = s + y        # Add to sum
        c = (t - s) - y  # Recover the roundoff error
        s = t

    return s
```

- **(a)** Implement this algorithm.
- **(b)** Test with the sum $S = 10000 + \pi + e$ using 6-digit decimal arithmetic. Trace through by hand.
- **(c)** What is the role of the variable `c`? Why does it improve accuracy?
- **(d)** Compare Kahan summation to naive summation for the integral:
  $$I = \int_1^{2^{22}} \frac{dx}{\sqrt{x}} = 2(2^{11} - 1)$$
  using the trapezoidal rule with $N = 10^6$ points, in single precision (`np.float32`).

---

## Condition Numbers

---

### Q2.8: Condition Number Computation

Compute the condition number $\kappa = |xf'(x)/f(x)|$ for:

- **(a)** $f(x) = e^x$ at $x = 1$
- **(b)** $f(x) = \ln(x)$ at $x = 1$ and $x = e$
- **(c)** $f(x) = \sin(x)$ at $x = 0$ and $x = \pi$
- **(d)** $f(x) = x^n$ at $x = 1$ for integer $n$

Which evaluations are well-conditioned? Which are ill-conditioned?

:::{dropdown} Solution for (a)
For $f(x) = e^x$, we have $f'(x) = e^x$.

$$
K = \left|\frac{x f'(x)}{f(x)}\right| = \left|\frac{x e^x}{e^x}\right| = |x|
$$

At $x = 1$: $K = 1$. This is well-conditioned.
:::

---

### Q2.9: Subtractive Cancellation

The expression $\sqrt{x+1} - \sqrt{x}$ suffers from cancellation for large $x$.

- **(a)** Compute this directly for $x = 10^{16}$.
- **(b)** Rationalize to get $\frac{1}{\sqrt{x+1} + \sqrt{x}}$ and compute again.
- **(c)** What is the relative error of the direct method?

---

### Q2.10: Trigonometric Reformulation

For small $x$, the expression $1 - \cos(x)$ loses accuracy.

- **(a)** Compute $1 - \cos(10^{-8})$ directly.
- **(b)** Use the identity $1 - \cos(x) = 2\sin^2(x/2)$ and compute again.
- **(c)** Compare to the Taylor approximation $1 - \cos(x) \approx x^2/2$.

---

## Stability Analysis

---

### Q2.11: Stable Reformulation

For each of the following expressions, identify the potential numerical issue and propose a stable reformulation:

- **(a)** $\sqrt{x^2 + 1} - x$ for large positive $x$

- **(b)** $\frac{1 - \cos(x)}{x^2}$ for small $x$

- **(c)** $e^x - 1$ for small $x$

- **(d)** $\ln(x) - \ln(y)$ when $x \approx y$

:::{dropdown} Solution for (a)
**Problem:** Subtractive cancellation when $\sqrt{x^2+1} \approx x$.

**Stable form:** Multiply by $\frac{\sqrt{x^2+1}+x}{\sqrt{x^2+1}+x}$:

$$
\sqrt{x^2+1} - x = \frac{(x^2+1) - x^2}{\sqrt{x^2+1}+x} = \frac{1}{\sqrt{x^2+1}+x}
$$
:::

---

## Computational Exercises

---

### Q2.12: Trapezoidal Rule with Round-off

Consider computing the integral:

$$
I = \int_1^{2^k} \frac{dx}{\sqrt{x}} = 2(2^{k/2} - 1)
$$

using the trapezoidal rule with $k = 22$. The approximation is:

$$
I \approx \Delta x \left(\sum_{i=1}^{N-1} f(x_i) + \frac{f(x_0) + f(x_N)}{2}\right)
$$

where $f(x) = x^{-1/2}$ and $\Delta x = \frac{2^k - 1}{N}$.

- **(a)** Implement this using naive summation with single precision (32-bit) floats. Plot the relative error vs. $N$ for $N = 10^p$, $p \in [2, 8]$. Explain what you observe.

- **(b)** Implement using Kahan summation. Plot the relative error vs. $N$. Explain what you observe.

:::{dropdown} Hint
To force single precision in NumPy:
```python
s = np.asarray(0.0, dtype=np.float32)
```
:::

---

### Q2.13: Fast Inverse Square Root Accuracy

Refer to the fast inverse square root code.

- **(a)** Compute the relative error of the fast inverse square root **without** Newton iteration (just the bit manipulation). Plot the relative error for inputs $x \in [10^0, 10^5]$.

- **(b)** Compute the relative error **with** Newton iterations. How many iterations would you use for a game engine? Support your answer with error plots for 1, 2, and 3 iterations.

---
