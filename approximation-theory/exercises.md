# Exercises

Practice applying Taylor's theorem and implementing finite difference methods. These exercises build both theoretical understanding and computational skills.

## Self-Assessment Questions

Test your understanding with these conceptual questions:

1. **Taylor's Theorem:** Can you write down Taylor's theorem with the Lagrange form of the remainder? What conditions on $f$ are required?

2. **Error Analysis:** Given a finite difference formula, can you derive its error using Taylor expansions?

3. **Order of Accuracy:** What does it mean for a method to be "first-order" or "second-order" accurate? How do you verify this computationally?

4. **Implementation:** Can you implement forward and central difference approximations in Python?

5. **Log-Log Plots:** Why do we use log-log plots to verify convergence rates? What slope should you expect for first-order and second-order methods?

6. **Round-off Error:** Why does the error start to increase for very small $h$? What is the optimal step size for forward differences with double precision?

## Computational Exercises

---

### Q1.1: Taylor Series Visualization

Consider the program `taylor.py` or the equivalent Jupyter notebook posted on Moodle.

- **(a)** Setup your computer so that you can execute either program.

- **(b)** Create a plot showing the first four Taylor series expansions of $f(x) = e^x$:
  - Use a different color for each line.
  - Plot NumPy's exponential function in dashed black.
  - Make sure your name appears in the plot title.
  - Put the absolute error of the approximation of $e$ for each graph in the plot's legend. Use two significant figures and scientific notation.

:::{dropdown} Hint
Use `np.exp(x)` for the exact function. For Taylor polynomials, you can compute them iteratively by adding terms.
:::

---

### Q1.2: Taylor Polynomials of Quadratics

Consider the polynomial $f(x) = x^2 - x - 2$.

- **(a)** Find $P_1(x)$, $P_2(x)$ and $P_3(x)$ for $f(x)$ centered at $x_0 = 0$. What is the relation between $P_3(x)$ and $f(x)$? Why?

- **(b)** Find $P_1(x)$, $P_2(x)$ and $P_3(x)$ for $f(x)$ centered at $x_0 = 2$. What is the relation between $P_3(x)$ and $f(x)$? Why?

- **(c)** In general, given a polynomial $f(x)$ with degree $\leq m$, what can you say about $f(x) - P_n(x)$ for $n \geq m$?

:::{dropdown} Solution
**(a)** At $x_0 = 0$:
- $f(0) = -2$, $f'(0) = -1$, $f''(0) = 2$, $f'''(0) = 0$
- $P_1(x) = -2 - x$
- $P_2(x) = -2 - x + x^2$
- $P_3(x) = -2 - x + x^2 = f(x)$

$P_3(x) = f(x)$ because $f$ is a degree 2 polynomial, and $P_n(x) = f(x)$ for all $n \geq 2$.

**(b)** At $x_0 = 2$:
- $f(2) = 0$, $f'(2) = 3$, $f''(2) = 2$
- $P_1(x) = 3(x-2)$
- $P_2(x) = 3(x-2) + (x-2)^2 = x^2 - x - 2 = f(x)$

**(c)** For any polynomial $f$ of degree $\leq m$ and $n \geq m$: $f(x) - P_n(x) = 0$.
:::

---

### Q1.3: Taylor Approximation of tan(x)

Find both $P_2(x)$ and $P_3(x)$ for $f(x) = \tan x$ about $x_0 = 0$, and use them to approximate $\tan(0.2)$.

*Hint: Plots may be helpful to answer the following questions.*

- **(a)** Show that in each case the remainder term provides an upper bound for the true error.

- **(b)** Using the remainder term, find a minimum value of $k$ for $P_k(x)$ to approximate $f(x)$ to within $5 \times 10^{-2}$ on $[0, 0.5]$.

:::{dropdown} Solution
For $f(x) = \tan x$ at $x_0 = 0$:
- $f(0) = 0$
- $f'(x) = \sec^2 x$, so $f'(0) = 1$
- $f''(x) = 2\sec^2 x \tan x$, so $f''(0) = 0$
- $f'''(x) = 2\sec^4 x + 4\sec^2 x \tan^2 x$, so $f'''(0) = 2$

Therefore:
- $P_2(x) = x$
- $P_3(x) = x + \frac{x^3}{3}$

Approximations: $P_2(0.2) = 0.2$, $P_3(0.2) = 0.2027$

The true value is $\tan(0.2) \approx 0.2027$, so $P_3$ is much more accurate.
:::

---

## Finite Difference Exercises

---

### Q1.4: Central Difference Analysis

Consider the central finite difference scheme to approximate $f'(x_0)$:

$$
f'(x_0) \approx \frac{f(x_0 + h) - f(x_0 - h)}{2h}
$$

- **(a)** Show that the approximation error is $E(h) = Ch^q$, where $C$ and $q$ are constants that you will determine.

- **(b)** Consider the function $f(x) = \tan(x)$ at $x_0 = 0.2$. Develop a function that approximates $f'(x_0)$ using the central difference scheme. Using that function, create a log-log plot of the absolute error against the step size $h$ for $h = 10^{-N}$, $N \in [0, 20]$. Also plot the function $E(h) = Ch^q$. Does the slope of this line match the slope of the absolute error? Explain what you observe.

:::{dropdown} Solution
**(a)** From the derivation in the notes: $q = 2$ and $C = \frac{1}{12}(f'''(\xi^+) + f'''(\xi^-))$.

**(b)** The log-log plot should show slope 2 for moderate values of $h$. For very small $h$, round-off error dominates and the error increases.
:::

---

### Q1.5: Second Derivative Approximation

Given a function $f(x)$, use Taylor's theorem to derive a second-order approximation of $f''(x_0)$.

*Hint: Your finite difference approximation for $f''(x_0)$ should include the terms $f(x_0 \pm h)$ and $f(x_0)$.*

- **(a)** What is the precise form of the error term?

- **(b)** Create the following plots, and for each explain what you observe:
  1. Plot the value $f''(1)$ against the discretization size $h$.
  2. Plot the absolute error vs. the discretization size $h$.
  3. Plot the absolute error divided by $h^2$ vs. the discretization size $h$. Verify that the values converge to the value derived from the error term.

:::{dropdown} Solution
The second derivative approximation is:

$$
f''(x_0) \approx \frac{f(x_0 + h) - 2f(x_0) + f(x_0 - h)}{h^2}
$$

The error is $E(h) = -\frac{h^2}{12}f^{(4)}(\xi)$ for some $\xi$ near $x_0$.

This is a second-order method ($\mathcal{O}(h^2)$).
:::

---

### Q1.6: One-Sided Second-Order Approximation

Given a function $f(x)$, use Taylor approximations to derive a **second-order** (i.e., the error should be $E(h) = Ch^2$) approximation to $f'(x_0)$.

*Hint: Your approximation should include the terms $f(x_0)$, $f(x_0 + h)$, and $f(x_0 + 2h)$.*

- **(a)** What is the precise form of the error term?

- **(b)** Using your formula, approximate $f'(0)$ where $f(x) = \sin(x)$ for $h = 2^{-p}$, $p \in [1, 15]$. Create:
  1. A semilog plot of $f'(0)$ against $h$.
  2. A log-log plot of the absolute error vs. $h$.
  3. A plot of the absolute error divided by $h^2$ vs. $h$.

:::{dropdown} Hint
Expand $f(x_0 + h)$ and $f(x_0 + 2h)$ in Taylor series, then find coefficients $a, b, c$ such that:

$$
af(x_0) + bf(x_0 + h) + cf(x_0 + 2h) = f'(x_0) + \mathcal{O}(h^2)
$$
:::

:::{dropdown} Solution
The formula is:

$$
f'(x_0) \approx \frac{-3f(x_0) + 4f(x_0+h) - f(x_0+2h)}{2h}
$$

This is useful when you can only sample to one side of $x_0$ (e.g., at a boundary).
:::

---

## Taylor Series Practice

---

### Q1.7: Computing Taylor Series

Compute Taylor series for the following functions. Report the first three nonzero terms and show your work.

- **(a)** $\sqrt{x}$ at $x = a$ (for $a > 0$)

- **(b)** $x^{3/2}$ at $x = a$ (for $a > 0$)

- **(c)** $e^x$ at $x = 0, 1, 2$

- **(d)** $\cos(x)$ at $x = 0, \pi/4, \pi/2$

- **(e)** $\sin(x)$ at $x = 0, \pi/4, \pi/2$

- **(f)** $\tan(x)$ at $x = 0, \pi/4$. Can you compute a Taylor series near $x = \pi/2$? Why or why not?

:::{dropdown} Solution for (f)
At $x = 0$: $\tan x = x + \frac{x^3}{3} + \frac{2x^5}{15} + \cdots$

At $x = \pi/4$: $\tan x = 1 + 2(x - \pi/4) + 2(x - \pi/4)^2 + \cdots$

At $x = \pi/2$: **No Taylor series exists** because $\tan(x)$ has a vertical asymptote at $x = \pi/2$ (it is not continuous there, let alone differentiable).
:::

---

