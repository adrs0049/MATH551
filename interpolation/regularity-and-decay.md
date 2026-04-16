---
kernelspec:
  name: python3
  display_name: Python 3
exports:
  - format: pdf
    template: ./_templates/plain_narrow
    output: exports/regularity-and-decay.pdf
    id: interpolation-regularity-and-decay-pdf
downloads:
  - id: interpolation-regularity-and-decay-pdf
    title: Download PDF
---

# Regularity and Coefficient Decay

:::{tip} Big Idea
The smoothness of $f$ controls the decay rate of its Chebyshev coefficients
$c_k$. Analytic functions give *geometric* decay $|c_k| = O(\rho^{-k})$.
Functions with $r$ continuous derivatives give *algebraic* decay
$|c_k| = O(k^{-r-1})$. Discontinuities give no decay at all and produce
**Gibbs oscillations** that do not vanish as $n \to \infty$. Reading a
coefficient-decay plot is therefore reading the regularity of $f$.
:::

## The Coefficient Decay Theorem

For Lipschitz $f$ on $[-1,1]$ the Chebyshev series

$$
f(x) = \sum_{k=0}^\infty c_k T_k(x), \qquad
c_k = \frac{2}{\pi} \int_{-1}^{1} \frac{f(x) T_k(x)}{\sqrt{1-x^2}}\, dx
\quad (k \ge 1)
$$

converges uniformly. The interesting question is *how fast*.

:::{prf:theorem} Algebraic decay
:label: thm-algebraic-decay

If $f, f', \ldots, f^{(r-1)}$ are absolutely continuous on $[-1,1]$ and
$f^{(r)}$ has bounded variation $V$, then for $k > r$
$$
|c_k| \le \frac{2V}{\pi (k-r)^{r+1}}.
$$
:::

:::{prf:theorem} Geometric decay
:label: thm-geometric-decay

If $f$ is analytic on $[-1,1]$ and extends analytically to the **Bernstein
ellipse** $\mathcal{E}_\rho$ (foci $\pm 1$, semi-axis sum $\rho > 1$), with
$|f| \le M$ there, then
$$
|c_k| \le 2 M \rho^{-k}.
$$
:::

The proofs are integrations-by-parts (algebraic case) and contour shifts
(analytic case); see {cite:t}`Trefethen2013`, Ch. 7–8.

## Reading Decay Plots

Take four functions of decreasing regularity and plot $|c_k|$.

```{code-cell} python
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt
import scipy.fft as fft

def chebpts(n):
    return np.cos(np.pi * np.arange(n+1) / n)

def vals2coeffs(v):
    n = len(v) - 1
    c = fft.dct(v[::-1], type=1, norm='forward')
    c[1:n] *= 2
    return c

n = 1024
x = chebpts(n)
fig, axes = plt.subplots(2, 2, figsize=(10, 6.5))

cases = [
    (np.sin(5*x),                  r'$\sin(5x)$ (entire)',          'C0'),
    (1/(1 + 25*x**2),              r'$1/(1+25x^2)$ (poles at $\pm i/5$)', 'C2'),
    (np.abs(x)**3,                 r'$|x|^3$ ($C^2$)',              'C1'),
    (np.sign(x),                   r'$\mathrm{sign}(x)$ (jump)',    'C3'),
]

for ax, (v, label, color) in zip(axes.flat, cases):
    c = vals2coeffs(v)
    ax.semilogy(np.abs(c) + 1e-20, '.', ms=3, color=color)
    ax.set_title(label)
    ax.set_xlabel('$k$'); ax.set_ylabel(r'$|c_k|$')
    ax.set_ylim(1e-18, 5)
plt.tight_layout(); plt.show()
```

What to read from each panel:

- **$\sin(5x)$** is entire (analytic everywhere). The coefficients drop
  faster than any geometric rate and hit machine precision near $k \sim 30$.
- **$1/(1+25x^2)$** is analytic on $[-1,1]$ but has complex poles at
  $\pm i/5$, *close* to the real interval. So the Bernstein-ellipse
  parameter $\rho$ is small ($\rho \approx 1.22$), and decay is geometric
  but slow. A clean straight line on log-y axes.
- **$|x|^3$** has two continuous derivatives ($r = 2$); the third has a
  jump. Coefficients decay like $k^{-3}$, a straight line on log-log axes,
  not on log-y.
- **$\mathrm{sign}(x)$** is discontinuous; coefficients decay only as
  $k^{-1}$, the slowest possible rate compatible with convergence in any
  averaged sense.

So a quick look at $|c_k|$ tells you which regularity class $f$ belongs to.

## A View from Sobolev Theory

This regularity↔decay correspondence is not specific to Chebyshev series.
It is the *discrete* face of a duality you may have seen in functional
analysis: smoothness in physical space corresponds to localization in
frequency space.

:::{prf:remark} The MATH 725 picture
:label: rmk-sobolev

For periodic $f$ on $[0, 2\pi]$ with Fourier coefficients $\hat f_k$,
membership in the Sobolev space $H^s$ is *equivalent* to
$\sum_k |k|^{2s} |\hat f_k|^2 < \infty$. Roughly, $f$ has $s$ derivatives in
$L^2$ iff $|\hat f_k|$ decays faster than $|k|^{-s-1/2}$. Under
$x = \cos\theta$ the Chebyshev series becomes a Fourier cosine series and
inherits the same dictionary:

| Physical-space regularity | Frequency-space decay |
|---|---|
| Analytic on a strip / Bernstein ellipse | Geometric: $\lvert c_k\rvert \le C\rho^{-k}$ |
| $C^r$ with bounded $r$-th variation | Algebraic: $\lvert c_k\rvert = O(k^{-r-1})$ |
| Continuous, no derivatives | $\lvert c_k\rvert = O(k^{-1})$ |
| Jump discontinuity | No decay (only conditional convergence) |

This is one face of the **uncertainty principle**: a function that is
extremely *localized* in physical space (a delta, a jump) is extremely
*spread out* in frequency, and vice versa. Spectral methods exploit this
duality every time they truncate a coefficient tail.
:::

## The Approximation Error

Truncation of the series gives the standard error bound

$$
\|f - p_n\|_\infty \le \sum_{k=n+1}^\infty |c_k|.
$$

Combined with the decay theorems:

| $f$ regularity | Approximation error |
|---|---|
| Analytic (Bernstein parameter $\rho$) | $O(\rho^{-n})$ |
| $C^r$ with bounded variation of $f^{(r)}$ | $O(n^{-r})$ |
| Continuous only | $O(n^{-1/2})$ in $L^2$, $O(1)$ in $L^\infty$ |

The jump from algebraic to geometric convergence is what makes spectral
methods spectacular for smooth problems and unhelpful for non-smooth ones.

## The Gibbs Phenomenon

When $f$ has a jump, the truncated series exhibits a persistent overshoot
near the discontinuity. The overshoot magnitude does not shrink as
$n \to \infty$. Only its *width* shrinks.

```{code-cell} python
:tags: [hide-input]

def coeffs2vals(c):
    n = len(c) - 1
    cs = c.copy()
    cs[1:n] /= 2
    v = fft.idct(cs, type=1, norm='forward')
    return v[::-1]

def cheb_interp(xeval, n):
    x = chebpts(n)
    f = np.sign(x)
    c = vals2coeffs(f)
    out = np.zeros_like(xeval)
    Tkm1 = np.ones_like(xeval); Tk = xeval.copy()
    out += c[0] * Tkm1
    if n >= 1:
        out += c[1] * Tk
    for k in range(2, n+1):
        Tkp1 = 2*xeval*Tk - Tkm1
        out += c[k] * Tkp1
        Tkm1, Tk = Tk, Tkp1
    return out

xx = np.linspace(-0.5, 0.5, 4000)
fig, ax = plt.subplots(figsize=(7.5, 4.2))
for n, color in zip([16, 64, 256], ['C0', 'C1', 'C3']):
    ax.plot(xx, cheb_interp(xx, n), color=color, lw=1, label=f'$n = {n}$')
ax.plot([-0.5, 0, 0, 0.5], [-1, -1, 1, 1], 'k--', lw=1, label='$\\mathrm{sign}(x)$')
ax.axhline(1.0895, color='gray', ls=':', lw=1, alpha=0.7,
           label=r'asymptotic peak $\approx 1.09$')
ax.set_xlabel('$x$'); ax.set_ylim(-1.4, 1.4)
ax.set_title('Gibbs phenomenon: overshoot persists, only narrows')
ax.legend(fontsize=9, loc='lower right'); plt.tight_layout(); plt.show()
```

The peak overshoot of the truncated series approaches

$$
\frac{2}{\pi} \int_0^{\pi} \frac{\sin t}{t}\, dt \approx 1.0895,
$$

an $\sim 8.9\%$ overshoot above the true jump height. This is the same
constant that appears for Fourier series. Under $x = \cos\theta$ the two
phenomena are identical.

The connection to coefficient decay is direct: $\lvert c_k\rvert = O(k^{-1})$
is too slow for the partial sums to converge uniformly near the jump, so
the limiting overshoot has nowhere to go.

The Sobolev picture from earlier sharpens this. Membership in $H^s$ is
equivalent to $\sum_k k^{2s}\,\lvert c_k\rvert^2 < \infty$. For a function
with $\lvert c_k\rvert \sim 1/k$ the sum is

$$
\sum_k k^{2s} \cdot k^{-2} \;=\; \sum_k k^{2s - 2},
$$

which converges iff $2s - 2 < -1$, that is $s < \tfrac{1}{2}$. So a
function with a jump sits exactly at the Sobolev threshold $s = \tfrac12$
and belongs to $H^s$ only for $s < \tfrac12$. The threshold $s = \tfrac12$
is also precisely the embedding boundary: $H^s \hookrightarrow C^0$ holds
for $s > \tfrac12$ but fails at and below. So a function below the
threshold cannot be continuous, the partial sums cannot converge uniformly
to it, and the persistent Gibbs overshoot is the geometric record of that
failure. The uncertainty principle made it inevitable: a function that
is sharply localised in physical space (the jump) must be spread out in
frequency, slowly enough that no finite truncation can reassemble it
without overshoot.

**If you must approximate discontinuous data, split the domain at the
jumps and do piecewise interpolation on each smooth piece.**

```{bibliography}
:filter: docname in docnames
```
