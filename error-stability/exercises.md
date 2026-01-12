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

7. **Forward vs. Backward:** If an algorithm has small backward error but large forward error, what does this tell you about the problem?

8. **Trade-offs:** In finite differences, why can't we just use an extremely small $h$ to get arbitrary accuracy?


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

## Forward and Backward Error

---

### Q2.5: Quadratic Formula

For $ax^2 + bx + c = 0$ with $a = 1$, $b = -10^8$, $c = 1$:

- **(a)** Compute both roots using the standard quadratic formula.
- **(b)** Identify which root suffers from cancellation.
- **(c)** Use the identity $x_1 x_2 = c/a$ to compute the small root accurately.
- **(d)** Compute the backward error (residual) for both methods.

---

### Q2.6: Exponential Minus One

- **(a)** For $x = 10^{-k}$, $k = 1, \ldots, 16$, compute $e^x - 1$ using direct subtraction.
- **(b)** Compare to `numpy.expm1(x)`.
- **(c)** Plot the relative error of the direct method vs. $x$.
- **(d)** At what $x$ does the direct method lose all significant digits?

---

### Q2.7: Variance Computation

Consider computing variance: $\sigma^2 = \frac{1}{n}\sum(x_i - \bar{x})^2$.

An alternative "one-pass" formula is: $\sigma^2 = \frac{1}{n}\sum x_i^2 - \bar{x}^2$.

- **(a)** For $x = [10^8, 10^8 + 1, 10^8 + 2]$, compute variance both ways.
- **(b)** Which method is more accurate? Why?
- **(c)** Research Welford's algorithm. Why is it preferred?

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

### Q2.12: Residuals and Forward Error

Solve the system:

$$
\begin{pmatrix} 1 & 1 \\ 1 & 1.0001 \end{pmatrix}
\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} =
\begin{pmatrix} 2 \\ 2 \end{pmatrix}
$$

- **(a)** Solve using Gaussian elimination.
- **(b)** Perturb $b_2$ by $0.0001$ and solve again.
- **(c)** Compute the condition number of $A$.
- **(d)** Verify that forward error $\approx \kappa \times$ relative perturbation.

---

## Computational Exercises

---

### Q2.13: Unstable Recursion

Consider evaluating the integrals

$$
y_n = \int_0^1 \frac{x^n}{10 + x}\, dx
$$

for $n = 1, 2, \dots, 30$.

- **(a)** Show that $y_n + 10y_{n-1} = \frac{1}{n}$.

- **(b)** Show that $y_0 = \log 11 - \log 10$, then use the recursion
  $$
  y_n = \frac{1}{n} - 10y_{n-1}
  $$
  to numerically generate $y_1$ through $y_{30}$. Be sure to use $y_0$ in the form **above** (not $\log(11/10)$). Include a table showing $n$ and $y_n$. Briefly explain what you observe.

- **(c)** Show that for $n \geq 0$, we have $0 \leq y_n \leq 1$. Discuss the results from (b) in light of this bound.

- **(d)** Derive a formula for computing $y_{n-1}$ given $y_n$ (the backward recursion).

- **(e)** Show that $y_n$ equals an infinite series.

:::{dropdown} Hint for (e)
Use the backward recursion $y_n = \frac{1}{10}\left(\frac{1}{n+1} - y_{n+1}\right)$ assuming $y_\infty = 0$. Repeatedly substitute to find the pattern.
:::

- **(f)** Show that for any $\varepsilon > 0$ and positive integer $n_0$, there exists $n_1 \geq n_0$ such that taking $y_{n_1} = 0$ as a starting value produces integral evaluations $y_n$ with absolute error smaller than $\varepsilon$ for all $0 < n \leq n_0$.

- **(g)** Explain why the backward algorithm is stable.

:::{dropdown} Hint for (g)
Include a small round-off error at each step. What happens to this error as you iterate backward?
:::

- **(h)** Write a computer program that computes $y_{20}$ within an absolute error of at most $10^{-5}$. Explain how you chose $n_1$.

:::{dropdown} Solution for (a)
Using integration by parts or direct manipulation:
$$
y_n + 10y_{n-1} = \int_0^1 \frac{x^n + 10x^{n-1}}{10+x}\,dx = \int_0^1 \frac{x^{n-1}(x+10)}{10+x}\,dx = \int_0^1 x^{n-1}\,dx = \frac{1}{n}
$$
:::

---

### Q2.14: Trapezoidal Rule with Round-off

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

### Q2.15: Fast Inverse Square Root Accuracy

Refer to the fast inverse square root code.

- **(a)** Compute the relative error of the fast inverse square root **without** Newton iteration (just the bit manipulation). Plot the relative error for inputs $x \in [10^0, 10^5]$.

- **(b)** Compute the relative error **with** Newton iterations. How many iterations would you use for a game engine? Support your answer with error plots for 1, 2, and 3 iterations.

---
