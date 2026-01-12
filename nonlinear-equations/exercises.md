# Exercises

Practice implementing and analyzing root-finding algorithms.

## Self-Assessment Questions

Test your understanding with these conceptual questions:

1. **IVT:** If $f(0) = -2$ and $f(1) = 3$, what can you conclude about roots of $f$ in $[0,1]$?

2. **Bisection:** How many bisection iterations are needed to reduce a starting interval of width 1 to width $10^{-6}$?

3. **Fixed Points:** Given $g(x) = \cos(x)$, does the iteration $x_{n+1} = \cos(x_n)$ converge? What is the fixed point?

4. **Convergence Rate:** If a sequence converges linearly with rate $C = 0.5$, approximately how many iterations to reduce error by a factor of $10^{-6}$?

5. **Newton's Method:** Why does Newton's method fail or converge slowly for $f(x) = x^3$ near $x = 0$?

---

## Bisection Method

---

### Q3.1: Basic Bisection

Consider the function $f(x) = \sqrt{x} - 1.1$ on the interval $[0, 2]$.

- **(a)** Use the bisection method to find a zero with $\texttt{atol} = 10^{-8}$. How many iterations are required? Does the iteration count match expectations from the convergence analysis?

- **(b)** What is the resulting absolute error? Could this error be predicted by the convergence analysis?

---

### Q3.2: Finding Roots with Bisection

Using bisection, determine the following roots to within $\texttt{atol} = 10^{-8}$. For each, report the solution (8 significant figures) and number of iterations.

- **(a)** The positive root of $x^2 - 3 = 0$

- **(b)** The real root of $x^3 - x^2 - x - 1 = 0$

- **(c)** The smallest positive root of $\cos(x) - \sin(x) - \frac{1}{2} = 0$

- **(d)** The smallest positive root of $\tan(x) - x = 0$

- **(e)** The root of $\tan(x) - x = 0$ closest to $x = 100$

---

## Fixed Point Iteration

---

### Q3.3: Cosine Fixed Point

Prove that for any $x_0 \in \mathbb{R}$, the iteration $x_{n+1} = \cos(x_n)$ converges to a unique fixed point $c$. Find $c$ to at least 12 decimal places.

:::{dropdown} Hint
It may be helpful to consider $h(x) = (g \circ g)(x) = \cos(\cos(x))$.
:::

---

### Q3.4: Sine Iteration

Consider the sequence $x_{n+1} = \sin(x_n)$ with $x_0 = 1$.

- **(a)** Show that there is a unique fixed point.

- **(b)** Show that the sequence is monotone decreasing using $|\sin(x)| < |x|$. Show it's bounded below by zero. Conclude convergence.

- **(c)** Plot the first 500 values on a semilogx plot. Estimate the slope as iteration number increases. What does this mean for the rate of convergence?

---

### Q3.5: Convergence Analysis

For each iteration, determine if it converges to the indicated $c$ when $x_0$ is sufficiently close. If it converges, find the order of convergence using both graphical analysis and Taylor series expansion.

- **(a)** $x_{n+1} = \frac{15x_n^2 - 24x_n + 13}{4x_n}$, $c = 1$

- **(b)** $x_{n+1} = \frac{3}{4}x_n + \frac{1}{x_n^3}$, $c = \sqrt{2}$

---

## Newton's Method

---

### Q3.6: Newton's Method Applications

Write programs to approximate the following values using Newton's method with $\texttt{atol} = 10^{-12}$. Use Taylor's theorem to find a good initial guess.

- **(a)** $\sqrt{2}$ and $\sqrt{3}$

- **(b)** $\pi$ from a root-finding problem

- **(c)** The positive root of $\cos(x) = x$

- **(d)** The smallest positive root of $\tan(x) = x$

:::{dropdown} Hint for (d)
You need an expansion that handles the singularity at $\pi/2$.
:::

---

### Q3.7: Newton's Method Basins of Attraction

Apply Newton's method to $f(x) = 4x^4 - 6x^2 - \frac{11}{4}$.

- **(a)** Plot the function. Use Newton's method with starting guesses close to the roots to find the two real roots to 8 digits.

- **(b)** For starting values `x = np.linspace(-2, 2, 200)`, determine which converge to positive/negative roots and which fail to converge.

- **(c)** Extend to complex starting values $z = x + yi$ for $-2 \leq x, y \leq 2$. Find all four roots. Create a plot showing regions of convergence for each root (Newton fractal).

---

### Q3.8: Multiple Roots

The function $f(x) = (x-1)^2 e^x$ has a double root at $x = 1$.

- **(a)** Derive Newton's method for this function. Show the iteration is well-defined for $x_n \neq -1$, and the convergence rate is similar to bisection.

- **(b)** Implement Newton's method starting from $x_0 = 2$. Observe its performance.

- **(c)** Could you use bisection for this problem? Explain.

---

### Q3.9: Division Without Division

Suppose your calculator's division button is broken (only +, −, × work).

- **(a)** Given $b \neq 0$, suggest a quadratically convergent iteration to compute $1/b$.

- **(b)** Implement your algorithm with convergence criterion $|x_n - x_{n-1}| < 10^{-10}$.

- **(c)** Apply to $b = \pi$ with initial guesses $x_0 = 1$ and $x_0 = 0.1$. Explain results.

- **(d)** Use bit operations (like in fast inverse square root) to propose an algorithm for a good initial guess.

- **(e)** Demonstrate numerically that your initial guess guarantees convergence.

:::{dropdown} Hint for (a)
To find $1/b$, solve $f(x) = 0$ where $f(x) = 1/x - b = 0$, or equivalently $f(x) = bx - 1 = 0$... but that uses division. Try $f(x) = 1 - bx$.
:::

---

### Q3.10: Newton in 2D

Consider the system:
$$
f(x,y) = 2x^2 - 2xy + 2y^2 - x - y = 0
$$
$$
g(x,y) = 4x - y + 2 = 0
$$

Find an approximate solution by taking **2 steps** of Newton's method starting from $(x_0, y_0) = (2, 0)$. Compute $(x_1, y_1)$ and $(x_2, y_2)$.

---

## Theory Questions

---

### Q3.11: Newton for Multiple Roots

If $c$ is a root of multiplicity $p \geq 2$, then $f(x) = (x-c)^p h(x)$ with $h(c) \neq 0$.

- **(a)** Write out Newton's iteration function $g(x)$ in terms of $h(x)$ and $h'(x)$.

- **(b)** Show that $g'(c) = 1 - 1/p \neq 0$. Explain why this implies only linear convergence.

:::{dropdown} Solution for (b)
Since $g'(c) = 1 - 1/p$ is nonzero (for $p \geq 2$), the fixed point theorem gives only linear convergence with rate $C = |1 - 1/p|$.

For a double root ($p = 2$): $C = 1/2$
For a triple root ($p = 3$): $C = 2/3$

The convergence slows as multiplicity increases.
:::

---

