# Norms and Condition Numbers

:::{tip} Big Idea
Norms measure the "size" of vectors and matrices. The condition number measures how sensitive a linear system is to perturbations—a large condition number means small input changes can cause large output changes.
:::

## Vector Norms

:::{prf:definition} Vector Norm
:label: def-vector-norm

A **norm** on $\mathbb{R}^n$ is a function $\|\cdot\|: \mathbb{R}^n \to \mathbb{R}$ satisfying:

1. $\|\mathbf{x}\| \geq 0$ with equality iff $\mathbf{x} = \mathbf{0}$ (positive definiteness)
2. $\|\alpha\mathbf{x}\| = |\alpha|\|\mathbf{x}\|$ for all scalars $\alpha$ (homogeneity)
3. $\|\mathbf{x} + \mathbf{y}\| \leq \|\mathbf{x}\| + \|\mathbf{y}\|$ (triangle inequality)
:::

### Common Vector Norms

The **$p$-norms** are defined as:

$$
\|\mathbf{x}\|_p = \left(\sum_{i=1}^n |x_i|^p\right)^{1/p}
$$

| Name | Formula | Interpretation |
|------|---------|---------------|
| 1-norm | $\|\mathbf{x}\|_1 = \sum_i |x_i|$ | Sum of absolute values |
| 2-norm (Euclidean) | $\|\mathbf{x}\|_2 = \sqrt{\sum_i x_i^2}$ | Length of vector |
| $\infty$-norm | $\|\mathbf{x}\|_\infty = \max_i |x_i|$ | Maximum component |

```python
import numpy as np

x = np.array([1, -2, 3])
print(np.linalg.norm(x, 1))    # 6.0
print(np.linalg.norm(x, 2))    # 3.74...
print(np.linalg.norm(x, np.inf))  # 3.0
```

## Matrix Norms

For matrices, we use **induced norms** (also called operator norms):

$$
\|A\| = \max_{\mathbf{x} \neq \mathbf{0}} \frac{\|A\mathbf{x}\|}{\|\mathbf{x}\|} = \max_{\|\mathbf{x}\| = 1} \|A\mathbf{x}\|
$$

This measures the maximum "stretching" the matrix does to any vector.

### Common Matrix Norms

| Name | Formula | Computation |
|------|---------|-------------|
| 1-norm | $\|A\|_1 = \max_j \sum_i |a_{ij}|$ | Max column sum |
| $\infty$-norm | $\|A\|_\infty = \max_i \sum_j |a_{ij}|$ | Max row sum |
| 2-norm (spectral) | $\|A\|_2 = \sqrt{\lambda_{\max}(A^T A)}$ | Largest singular value |
| Frobenius | $\|A\|_F = \sqrt{\sum_{i,j} a_{ij}^2}$ | Like vector 2-norm |

:::{admonition} Note
The 1-norm and $\infty$-norm are easy to compute (just sums). The 2-norm requires computing singular values, making it more expensive.
:::

### Key Properties

For any induced matrix norm:

1. $\|A\mathbf{x}\| \leq \|A\| \cdot \|\mathbf{x}\|$
2. $\|AB\| \leq \|A\| \cdot \|B\|$ (submultiplicativity)
3. $\|I\| = 1$

## Condition Number

:::{prf:definition} Condition Number
:label: def-matrix-condition-number

The **condition number** of an invertible matrix $A$ is:

$$
\kappa(A) = \|A\| \cdot \|A^{-1}\|
$$
:::

**Properties:**
- $\kappa(A) \geq 1$ for any matrix norm
- $\kappa(I) = 1$
- For the 2-norm: $\kappa_2(A) = \sigma_{\max}/\sigma_{\min}$ (ratio of largest to smallest singular values)

## Sensitivity of Linear Systems

Consider the perturbed system $(A + \delta A)(\mathbf{x} + \delta\mathbf{x}) = \mathbf{b} + \delta\mathbf{b}$.

:::{prf:theorem} Perturbation Bound for Linear Systems
:label: thm-perturbation-bound

If $\|\delta A\| < \|A^{-1}\|^{-1}$, then:

$$
\frac{\|\delta\mathbf{x}\|}{\|\mathbf{x}\|} \leq \kappa(A) \left(\frac{\|\delta A\|}{\|A\|} + \frac{\|\delta\mathbf{b}\|}{\|\mathbf{b}\|}\right) + O(\text{higher order})
$$
:::

**Interpretation:** The relative error in the solution can be amplified by a factor of $\kappa(A)$.

This connects directly to the forward/backward error framework:

$$
\text{forward error} \lesssim \kappa \times \text{backward error}
$$

## Well-Conditioned vs Ill-Conditioned

| Condition Number | Classification | Digits Lost |
|-----------------|----------------|-------------|
| $\kappa \approx 1$ | Well-conditioned | ~0 |
| $\kappa \approx 10^k$ | Moderate | ~$k$ digits |
| $\kappa \approx 1/\varepsilon_{\text{mach}}$ | Ill-conditioned | All digits |

:::{prf:example} Hilbert Matrix
:label: ex-hilbert-matrix
:class: dropdown

The Hilbert matrix $H_{ij} = \frac{1}{i+j-1}$ is notoriously ill-conditioned:

| Size | $\kappa_2(H)$ |
|------|--------------|
| $5 \times 5$ | $4.8 \times 10^5$ |
| $10 \times 10$ | $1.6 \times 10^{13}$ |
| $15 \times 15$ | $\approx 10^{17}$ |

For $n \geq 12$, the Hilbert matrix is effectively singular in double precision!
:::

## Computing Condition Numbers

```python
import numpy as np

A = np.array([[1, 2], [3, 4]])

# 2-norm condition number (default)
kappa = np.linalg.cond(A)

# Other norms
kappa_1 = np.linalg.cond(A, 1)
kappa_inf = np.linalg.cond(A, np.inf)
```

:::{admonition} Warning
:class: warning
Computing $\kappa(A)$ requires computing $A^{-1}$ or singular values—expensive for large matrices! Often we estimate condition numbers instead.
:::

## Why Condition Numbers Matter for QR

Recall from forward/backward error analysis:

$$
\frac{\|\hat{x} - x\|}{\|x\|} \lesssim \kappa(A) \cdot \varepsilon_{\text{mach}}
$$

For least squares with the normal equations $A^T A \hat{x} = A^T b$:

$$
\kappa(A^T A) = \kappa(A)^2
$$

This **squaring** of the condition number is disastrous! If $\kappa(A) = 10^6$, then $\kappa(A^T A) = 10^{12}$.

**QR factorization avoids this** by never forming $A^T A$.

## Practical Guidelines

1. **Always check condition numbers** before trusting numerical results
2. **Rescaling** (preconditioning) can improve conditioning
3. **No algorithm can overcome** inherent ill-conditioning—it's a property of the problem
4. **Prefer QR over normal equations** for least squares to avoid squaring $\kappa$

## Application: Cell Friction Networks

In cell migration models, the friction matrix $F$ relates cell velocities to forces:

$$
F \mathbf{v} = \mathbf{f}
$$

The condition number of $F$ depends on:
- **Network topology**: Well-connected networks are better conditioned
- **Friction coefficient ratios**: Large variations increase $\kappa$
- **Near-singular configurations**: Cells about to detach

Monitoring $\kappa(F)$ can detect numerical instabilities before they cause simulation failures.
