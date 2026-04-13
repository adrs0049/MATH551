---
kernelspec:
  name: python3
  display_name: Python 3
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/integration.pdf
    id: interpolation-integration-pdf
downloads:
  - id: interpolation-integration-pdf
    title: Download PDF
---

# Integration as a Linear Functional

:::{tip} Big Idea
$\int_{-1}^{1} f\,dx$ is a *linear functional* on $f$. Restricted to
$\mathbb{P}_n$ it is just a row vector of weights $\mathbf{w}$ acting on the
nodal values: $\int p_n = \mathbf{w}^T \mathbf{f}$. Every interpolatory
quadrature is determined by its node placement and the corresponding weight
vector. Picking Chebyshev nodes and using the closed-form Chebyshev
integrals gives **Clenshaw–Curtis quadrature**. Its weights are computable
from a DCT, and it has spectral accuracy on smooth integrands.
:::

## Quadrature as Weights on Values

The simplest quadrature you know is the [trapezoidal
rule](../approximation-theory/numerical-integration.md), which integrates
the *piecewise linear* interpolant of $f$ exactly. That is the model: pick
a basis, integrate the basis functions analytically, and the integral of
the interpolant is

$$
\int_{-1}^{1} p_n(x)\, dx
= \int_{-1}^{1} \sum_{j=0}^n f_j\, \ell_j(x)\, dx
= \sum_{j=0}^n \underbrace{\Big(\int_{-1}^{1} \ell_j(x)\, dx\Big)}_{w_j} f_j
= \sum_{j=0}^n w_j\, f_j
= \mathbf{w}^\top \mathbf{f}
= \langle \mathbf{w},\, \mathbf{f}\rangle.
$$

The continuous integral $\int p_n$ has become a *finite-dimensional inner
product* on $\mathbb{R}^{n+1}$: the value vector $\mathbf{f}$ tested
against the weight vector $\mathbf{w}$. Once we fix the nodes (and hence
the $\ell_j$), the **quadrature weights** $w_j = \int \ell_j$ are the
integrals of the basis functions. Trapezoidal uses the hat basis on
an equispaced grid, gets weights $w_j = h$ (with $h/2$ at the endpoints),
and integrates affine functions exactly. Equispaced higher-order rules
(Simpson, Newton–Cotes) generalize this. For large $n$ they suffer from
the same Runge instability we saw in [§2](point-choice.md). The weights
become large with mixed signs. Don't go that route.

The fix is the same one as before: keep the framework, change the nodes.

## Clenshaw–Curtis Quadrature

Use Chebyshev nodes and the Chebyshev basis. The Chebyshev expansion
$f \approx \sum_k c_k T_k(x)$ integrates term by term using

$$
\int_{-1}^{1} T_k(x)\, dx =
\begin{cases}
0, & k \text{ odd}, \\
\dfrac{2}{1 - k^2}, & k \text{ even}.
\end{cases}
$$

So

$$
\int_{-1}^{1} p_n(x)\, dx = 2 c_0 + \sum_{\substack{k = 2 \\ k \text{ even}}}^{n} \frac{2 c_k}{1 - k^2}.
$$

Two views, one calculation:

- *Coefficient view.* Compute $\mathbf{c} = $ DCT$(\mathbf{f})$ in $O(N \log N)$,
  then dot with the Chebyshev integral weights above.
- *Value view.* Convolve the integral weights with the DCT and you get a
  weight vector $\mathbf{w}$ such that $\int p_n = \mathbf{w}^T \mathbf{f}$
  directly. These are the **Clenshaw–Curtis weights**.

### What the weights look like

```{code-cell} python
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt
import scipy.fft as fft

def chebpts(N):
    return np.cos(np.pi * np.arange(N+1) / N)

def cc_weights(N):
    """Clenshaw–Curtis weights on Chebyshev points x_j = cos(j*pi/N).
    Trefethen, "Spectral Methods in MATLAB", clencurt.m."""
    theta = np.pi * np.arange(N+1) / N
    w = np.zeros(N+1)
    v = np.ones(N-1)
    if N % 2 == 0:
        w[0] = 1.0 / (N**2 - 1); w[N] = w[0]
        for k in range(1, N//2):
            v -= 2 * np.cos(2*k*theta[1:N]) / (4*k**2 - 1)
        v -= np.cos(N*theta[1:N]) / (N**2 - 1)
    else:
        w[0] = 1.0 / N**2; w[N] = w[0]
        for k in range(1, (N-1)//2 + 1):
            v -= 2 * np.cos(2*k*theta[1:N]) / (4*k**2 - 1)
    w[1:N] = 2 * v / N
    return w

N = 32
w = cc_weights(N); x = chebpts(N)
fig, ax = plt.subplots(figsize=(7, 3.8))
ax.stem(np.arange(N+1), w, basefmt=' ')
ax.set_xlabel('node index $j$'); ax.set_ylabel('$w_j$')
ax.set_title(f'Clenshaw–Curtis weights, $N = {N}$ (all positive)')
plt.tight_layout(); plt.show()
```

All weights are positive (no cancellation in the sum) and smoothly
distributed: the endpoint weights are small to compensate for the
endpoint clustering of Chebyshev nodes. Compare to Newton–Cotes at large
$n$, where weights alternate sign and grow exponentially.

## Convergence

The quadrature error inherits the truncation error of the Chebyshev
expansion, so the [regularity-decay](regularity-and-decay.md) dictionary
applies directly:

| $f$ regularity | $\lvert I - I_N\rvert$ |
|---|---|
| Analytic in a Bernstein ellipse with parameter $\rho$ | $O(\rho^{-N})$ |
| $C^r$ with bounded variation of $f^{(r)}$ | $O(N^{-r})$ |
| Discontinuous | $O(N^{-1})$ at best |

```{code-cell} python
:tags: [hide-input]

from scipy.special import erf

def trapezoid(f, N):
    x = np.linspace(-1, 1, N+1); h = x[1] - x[0]
    fx = f(x)
    return h * (fx[0]/2 + fx[-1]/2 + fx[1:-1].sum())

def cc(f, N):
    x = chebpts(N); return cc_weights(N) @ f(x)

# (a) Analytic integrand
f1 = lambda x: np.exp(-x**2)
I1 = np.sqrt(np.pi) * erf(1.0)

# (b) C^2 integrand: |x|^3
f2 = lambda x: np.abs(x)**3
I2 = 0.5

Ns = 2**np.arange(2, 11)
err_t1 = [abs(trapezoid(f1, N) - I1) for N in Ns]
err_c1 = [abs(cc(f1, N)        - I1) for N in Ns]
err_t2 = [abs(trapezoid(f2, N) - I2) for N in Ns]
err_c2 = [abs(cc(f2, N)        - I2) for N in Ns]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
ax1.loglog(Ns, err_t1, 'o-', label='trapezoid')
ax1.loglog(Ns, err_c1, 's-', label='Clenshaw–Curtis')
ax1.set_title(r'analytic: $\int_{-1}^{1} e^{-x^2} dx$')
ax1.set_xlabel('$N$'); ax1.set_ylabel('error'); ax1.legend()

ax2.loglog(Ns, err_t2, 'o-', label='trapezoid')
ax2.loglog(Ns, err_c2, 's-', label='Clenshaw–Curtis')
ax2.loglog(Ns, 1.0/Ns**4, 'k:', label=r'$N^{-4}$')
ax2.set_title(r'$C^2$: $\int_{-1}^{1} |x|^3 dx$')
ax2.set_xlabel('$N$'); ax2.legend()
plt.tight_layout(); plt.show()
```

The analytic integrand drops to machine precision around $N \approx 30$
under Clenshaw–Curtis; trapezoid plods along at $O(N^{-2})$. For the
$C^2$ integrand both rules eventually drop algebraically, but
Clenshaw–Curtis benefits from each extra derivative the integrand has,
while trapezoid is permanently stuck at second order.

## References

```{bibliography}
:filter: docname in docnames
```
