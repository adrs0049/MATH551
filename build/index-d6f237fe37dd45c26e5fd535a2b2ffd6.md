# MATH 551 Introduction to Scientific Computing

## Instructor

This is a course developed by Brian van Koten and [Andreas Buttenschoen](https://www.buttenschoen.ca) at the University of Massachusetts Amherst.

## Description

Math 551 is a first course in numerical analysis that introduces students to the
fundamental concepts and techniques of scientific computing. We will study
numerical algorithms, floating-point arithmetic and round-off errors, root-finding
methods for nonlinear equations, and numerical linear algebra including direct
and iterative methods for solving linear systems. The course emphasizes both
theoretical understanding and practical implementation using Python.

:::{note}
**These notes are a work in progress.** If you find bugs, missing references, typos, or other issues, please [file an issue on GitHub](https://github.com/adrs0049/MATH551/issues).
:::

## Learning Goals

* Understand floating-point arithmetic and sources of numerical error
* Apply Taylor's theorem to analyze the accuracy of numerical approximations
* Implement and analyze root-finding algorithms (bisection, Newton's method, fixed-point iteration)
* Solve linear systems using direct methods (LU factorization, Gaussian elimination)
* Apply iterative methods (Jacobi, Gauss-Seidel) for large sparse systems
* Construct and use polynomial interpolants (Lagrange, Newton, Chebyshev)
* Perform matrix computations using Python, NumPy, and SciPy

## Prerequisites

* MATH 233 (Multivariate Calculus)
* MATH 235 (Introduction to Linear Algebra)
* A scientific programming course: COMPSCI 121, E&C-ENG 122, PHYSICS 281, or E&C-ENG 242

## How to Use This Book

**Learning Outcomes** appear at the start of each chapter. Use them to check your understanding—if you can do what each outcome describes, you've mastered that material.

**Exercises** come in three types:
- *Self-assessment questions* — quick conceptual checks
- *Pencil-and-paper problems* — work through by hand to build intuition
- *Computational exercises* — implement algorithms in Python

Problems marked *(optional)* go beyond the core material.

**Notebooks** in the [Interactive Notebooks](notebooks/index.md) section contain runnable code. Click the Colab badge to open in Google Colaboratory, or use the launch button for Binder/JupyterHub.

## Course Topics

1. **Numerical Algorithms** - Taylor polynomials, finite difference approximations, error analysis
2. **Round-off Errors** - Floating point arithmetic, machine epsilon, catastrophic cancellation
3. **Nonlinear Equations** - Bisection method, Newton's method, fixed point iteration
4. **Linear Algebra Background** - Matrix norms, condition numbers, sensitivity analysis
5. **Direct Methods** - LU factorization, Gaussian elimination, Cholesky decomposition
6. **Iterative Methods** - Jacobi method, Gauss-Seidel, convergence analysis
7. **Polynomial Interpolation** - Lagrange interpolation, Newton form, Chebyshev polynomials

## Mathematical Python Resources

* [Mathematical Python](https://patrickwalls.github.io/mathematicalpython/) - Great resource to get started with mathematical Python
* [NumPy Documentation](https://numpy.org/doc/)
* [SciPy Documentation](https://docs.scipy.org/doc/scipy/)

## Course Grading Tools

* [PLOM Grading](https://plomgrading.org/) will be used to grade quizzes and exams.

## License

This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.

<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a>
