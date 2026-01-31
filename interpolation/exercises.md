# Exercises

Practice polynomial interpolation and spectral methods.

## Part I: Polynomial Interpolation

---

### Q1: Lagrange Interpolation

Find the Lagrange interpolating polynomial for:

- **(a)** $(0, 1), (1, 3)$
- **(b)** $(-1, 1), (0, 0), (1, 1)$
- **(c)** $(0, 1), (1, 1), (2, 7)$

---

### Q2: Divided Differences

Compute the divided difference table and Newton interpolant for $(1, 1), (2, 4), (3, 9), (4, 16)$.

What polynomial do you get? Is this expected?

---

### Q3: Error Estimation

For $f(x) = e^x$ interpolated at $x_0 = 0, x_1 = 1$:

- **(a)** Find the interpolating polynomial.
- **(b)** Use the error formula to bound $|f(0.5) - p(0.5)|$.
- **(c)** Compare with the actual error.

---

### Q4: Runge's Phenomenon

Interpolate $f(x) = \frac{1}{1 + 25x^2}$ on $[-1, 1]$ using:

- **(a)** $n = 5$ equally spaced nodes
- **(b)** $n = 10$ equally spaced nodes
- **(c)** $n = 20$ equally spaced nodes
- **(d)** $n = 20$ Chebyshev nodes

Plot each interpolant along with $f(x)$. What do you observe?

---

### Q5: Implementation

Implement:

- **(a)** Lagrange interpolation
- **(b)** Newton interpolation with divided differences
- **(c)** Horner's method for evaluating Newton form

Compare timings for evaluating at 1000 points with $n = 50$ nodes.

---

### Q6: Adding Points

You have the Newton interpolant for $(0, 1), (1, 0), (2, 3)$.

- **(a)** Add the point $(3, 10)$ efficiently.
- **(b)** What is the new polynomial?

---

### Q7: Uniqueness Proof

Given a set of points $x_0 < x_1 < \cdots < x_n$, show that there exists a **unique** polynomial of degree at most $n$ that interpolates a function at these points.

*Hint:* Proceed by contradiction. Assume two different interpolating polynomials exist and consider their difference $d(x)$. How many roots does $d(x)$ have?

---

### Q8: Barycentric Formula Derivation

This problem fills in the details of the barycentric interpolation derivation.

**(a)** Given the node polynomial $\ell(x) = \prod_{k=0}^{n}(x - x_k)$, show that:
$$
\ell_j(x) = \frac{\ell(x)}{\ell'(x_j)(x - x_j)}
$$

**(b)** Show that $\sum_{j=0}^{n} \ell_j(x) = 1$.

*Hint:* Use the uniqueness of polynomial interpolation—what function do the Lagrange polynomials interpolate when all data values are 1?

---

### Q9: Chebyshev Polynomials

The Chebyshev polynomials are defined as $T_k(x) = \frac{1}{2}(z^k + z^{-k})$ where $z = e^{i\theta}$ and $x = \cos\theta$.

**(a)** Show that $T_k(x) = \cos(k\theta)$ for an appropriate choice of $\theta$.

**(b)** Compute $T_0$, $T_1$, $T_2$, and $T_3$ explicitly as polynomials in $x$.

**(c)** Write $x^5$ as a sum of Chebyshev polynomials.

**(d)** Prove the recurrence relation: $T_{k+1}(x) = 2xT_k(x) - T_{k-1}(x)$.

---

### Q10: Quadrature Comparison

Consider the integral $I = \int_{-1}^{1} |\sin(5x)|^3 \, dx$. The exact value is:
$$
I = \frac{24 + 9\cos(5) - \cos(15)}{30}
$$

**(a)** Implement the trapezoidal rule with $n+1$ equispaced points.

**(b)** Implement Chebyshev quadrature (Clenshaw-Curtis) with $n+1$ Chebyshev points.

**(c)** For $n = 2^k$ with $k = 1, 2, \ldots, 12$, compute the absolute error for both methods. Plot on a log-log scale. What convergence rates do you observe?

**(d)** Repeat for the Gaussian integral:
$$
I = \frac{1}{\sqrt{\pi}\varepsilon}\int_{-1}^{1} e^{-(x/\varepsilon)^2} \, dx = \text{erf}(1/\varepsilon)
$$
with $\varepsilon = 0.5$. How do the methods compare for this smooth function?

---

## Part II: Spectral Differentiation

---

### Q11: Chebyshev Differentiation Matrix

- **(a)** Implement the Chebyshev differentiation matrix for $N+1$ points.
- **(b)** Verify that $D \cdot \mathbf{1} = \mathbf{0}$ (derivative of constant is zero).
- **(c)** For $u(x) = x^3$, compute $D\mathbf{u}$ and compare to exact $3x^2$.

---

### Q12: Spectral vs. Finite Difference

For $u(x) = e^{\sin(\pi x)}$ on $[-1, 1]$:

- **(a)** Compute the derivative using a 2nd-order centered finite difference.
- **(b)** Compute the derivative using the Chebyshev differentiation matrix.
- **(c)** Plot max error vs. $N$ for both methods on a semilog plot.
- **(d)** How many points does each method need for 10 digits of accuracy?

---

### Q13: Differentiation Matrix Properties

Consider the function $f(x) = \sin(\pi x)$ on $[-1, 1]$.

**(a)** Implement the Chebyshev differentiation matrix construction.

**(b)** Compute the derivative at the Chebyshev points using $\mathbf{f}' = D\mathbf{f}$.

**(c)** Compare with the exact derivative $f'(x) = \pi\cos(\pi x)$ at the same points.

**(d)** Implement a second-order finite difference approximation on an equispaced grid with the same number of points. Compare the errors.

---

## Part III: Spectral Methods for BVPs

---

### Q14: Simple BVP

Solve $u'' = -\pi^2 \sin(\pi x)$ on $[-1, 1]$ with $u(\pm 1) = 0$.

- **(a)** Find the exact solution.
- **(b)** Implement spectral collocation.
- **(c)** Plot error vs. $N$. Is convergence exponential?

---

### Q15: Variable Coefficient BVP

Solve $(1 + x^2)u'' + 2xu' - u = 0$ on $[-1, 1]$ with $u(-1) = 1$ and $u(1) = 2$.

- **(a)** Set up the spectral collocation system.
- **(b)** Solve for $N = 8, 16, 32$.
- **(c)** How can you assess convergence when you don't know the exact solution?

---

### Q16: Neumann Boundary Conditions

Solve $u'' = \cos(x)$ on $[-1, 1]$ with $u'(-1) = 0$ and $u(1) = 0$.

- **(a)** Modify the spectral collocation setup for mixed boundary conditions.
- **(b)** Solve and compare to the exact solution.
- **(c)** What changes if both conditions are Neumann? Is the problem well-posed?

---

### Q17: Eigenvalue Problem

For the Laplacian eigenvalue problem $-u'' = \lambda u$ on $[-1, 1]$ with $u(\pm 1) = 0$:

- **(a)** The exact eigenvalues are $\lambda_k = (k\pi/2)^2$. Compute them spectrally.
- **(b)** Plot the error in the first 10 eigenvalues vs. $N$.
- **(c)** Why are higher eigenvalues less accurate?

---

### Q18: Helmholtz Equation

Solve $u'' + k^2 u = f(x)$ on $[-1, 1]$ with $u(\pm 1) = 0$, where $k = 10$ and $f(x) = 1$.

- **(a)** Solve using spectral collocation.
- **(b)** What happens as $k$ approaches an eigenvalue of $-d^2/dx^2$?
- **(c)** How many points are needed for a given accuracy as $k$ increases?

---

### Q19: Comparison with scipy.solve_bvp

Use `scipy.integrate.solve_bvp` to solve:
$$
u'' = e^u, \quad u(-1) = u(1) = 0
$$

- **(a)** Solve using `solve_bvp` (finite difference based).
- **(b)** Solve using spectral collocation (need Newton iteration for nonlinear!).
- **(c)** Compare accuracy and number of points needed.

---

### Q20: Fourth-Order Problem

Solve the beam equation $u'''' = 1$ on $[-1, 1]$ with clamped boundary conditions: $u(\pm 1) = u'(\pm 1) = 0$.

- **(a)** Form $D^4 = D \cdot D \cdot D \cdot D$.
- **(b)** How do you impose 4 boundary conditions?
- **(c)** Solve and compare to the exact solution.

---

### Q21: Periodic Problem with Fourier

Solve $u'' + u = \cos(2x)$ on $[0, 2\pi]$ with periodic boundary conditions.

- **(a)** Use FFT-based differentiation instead of Chebyshev.
- **(b)** Compare to the exact solution.
- **(c)** Why is Fourier better than Chebyshev for periodic problems?

---

## Self-Assessment Questions

Test your understanding with these conceptual questions:

### Polynomial Interpolation

1. **Uniqueness:** If you have 4 data points, what is the maximum degree of the interpolating polynomial?

2. **Lagrange:** What is $L_2(x_2)$? What is $L_2(x_0)$?

3. **Error:** Why does the interpolation error formula look similar to Taylor's theorem error?

4. **Runge:** Why do equally spaced nodes cause problems for some functions?

5. **Barycentric:** Why is the second barycentric formula "scale invariant"? Why does this matter numerically?

6. **Chebyshev:** Why do Chebyshev points cluster near the endpoints of $[-1, 1]$?

7. **Coefficient Decay:** If a function has a jump discontinuity, what decay rate do you expect for its Chebyshev coefficients?

8. **Gibbs:** Does increasing the polynomial degree eliminate the overshoot near a discontinuity? Why or why not?

### Spectral Methods

9. **Differentiation Matrix:** Why is the Chebyshev differentiation matrix dense while finite difference matrices are sparse?

10. **Accuracy:** Why do spectral methods achieve exponential convergence for smooth problems?

11. **Boundary Conditions:** How do you modify the spectral system to impose Dirichlet vs. Neumann boundary conditions?

12. **Eigenvalues:** Why are the high-frequency eigenvalues of $D^2$ less accurate than low-frequency ones?

13. **Limitations:** When would you choose finite differences over spectral methods?

14. **Grid Points:** Why do Chebyshev points cluster near the boundaries?

15. **Periodic vs. Non-periodic:** When should you use Fourier methods vs. Chebyshev methods?
