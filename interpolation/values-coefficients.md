# Values and Coefficients

:::{tip} Big Idea
A polynomial can be represented by its **values** at Chebyshev points or by its **coefficients** in the Chebyshev basis. The DCT transforms between these representations in $O(n \log n)$ time. Each view has advantages: values for pointwise operations, coefficients for integration and smoothness analysis. This duality is the foundation of modern spectral methods.
:::

## The Two Faces of Polynomial Interpolants

Given $n+1$ data points $(x_j, f_j)$ where $x_j$ are Chebyshev points, the polynomial interpolant $p_n(x)$ is unique. But we can represent it in two equivalent ways:

### Value Representation

Using Lagrange polynomials $\ell_j(x)$:
$$
p_n(x) = \sum_{j=0}^{n} f_j \ell_j(x) \quad \longleftrightarrow \quad \mathbf{p}_n \equiv \begin{pmatrix} f_0 \\ f_1 \\ \vdots \\ f_n \end{pmatrix}
$$

The polynomial is uniquely determined by the **sampled function values** $\{f_j\}$.

### Coefficient Representation

Using Chebyshev polynomials $T_k(x)$:
$$
p_n(x) = \sum_{k=0}^{n} c_k T_k(x) \quad \longleftrightarrow \quad \mathbf{p}_n \equiv \begin{pmatrix} c_0 \\ c_1 \\ \vdots \\ c_n \end{pmatrix}
$$

The polynomial is uniquely determined by its **Chebyshev coefficients** $\{c_k\}$.

### The Connection: A Linear Map

Evaluating the Chebyshev series at the data points gives:
$$
f_i = p_n(x_i) = \sum_{k=0}^{n} c_k T_k(x_i)
$$

This is a matrix-vector product:

$$
\begin{pmatrix} f_0 \\ f_1 \\ \vdots \\ f_n \end{pmatrix}
=
\begin{pmatrix}
T_0(x_0) & T_1(x_0) & \cdots & T_n(x_0) \\
T_0(x_1) & T_1(x_1) & \cdots & T_n(x_1) \\
\vdots & \vdots & \ddots & \vdots \\
T_0(x_n) & T_1(x_n) & \cdots & T_n(x_n)
\end{pmatrix}
\begin{pmatrix} c_0 \\ c_1 \\ \vdots \\ c_n \end{pmatrix}
$$

Since $T_k(x_j) = T_k(\cos(j\pi/n)) = \cos(jk\pi/n)$, this is exactly the **Discrete Cosine Transform** matrix!

## Two Representations

Consider a polynomial $p(x)$ of degree at most $n$.

### Value Space

Store the values at the $n+1$ Chebyshev points:
$$
\mathbf{f} = \begin{pmatrix} f_0 \\ f_1 \\ \vdots \\ f_n \end{pmatrix} = \begin{pmatrix} p(x_0) \\ p(x_1) \\ \vdots \\ p(x_n) \end{pmatrix}
$$

where $x_k = \cos(k\pi/n)$ for $k = 0, 1, \ldots, n$.

**Advantages:**
- Direct sampling: just evaluate $f(x_k)$
- Differentiation via matrices: $p'(x_k) = (D\mathbf{f})_k$
- Nonlinear operations easy: $\sin(f)$ is just componentwise

### Coefficient Space

Store the Chebyshev coefficients:
$$
\mathbf{c} = \begin{pmatrix} c_0 \\ c_1 \\ \vdots \\ c_n \end{pmatrix} \quad \text{where} \quad p(x) = \sum_{k=0}^{n} c_k T_k(x)
$$

**Advantages:**
- Integration has a simple formula
- Coefficient decay reveals smoothness
- Truncation gives best polynomial approximation

## The Discrete Cosine Transform

The key insight: under $x = \cos\theta$, Chebyshev polynomials become cosines:
$$
T_k(\cos\theta) = \cos(k\theta)
$$

At Chebyshev points $x_j = \cos(j\pi/n)$, we have $\theta_j = j\pi/n$, so:
$$
T_k(x_j) = \cos\left(\frac{jk\pi}{n}\right)
$$

This is exactly the DCT-I matrix!

### Values to Coefficients

Given values $\mathbf{f}$, the coefficients are:
$$
c_k = \frac{2}{n} \sum_{j=0}^{n}{}'' f_j \cos\left(\frac{jk\pi}{n}\right)
$$

where $\sum''$ means the first and last terms are halved.

```python
import scipy.fft as fft

def vals2coeffs(values):
    """Convert values at Chebyshev points to coefficients."""
    n = len(values) - 1
    if n == 0:
        return values.copy()

    # Use DCT-I
    coeffs = fft.dct(values[::-1], type=1, norm='forward')

    # Scale interior coefficients
    coeffs[1:n] *= 2.0
    return coeffs
```

### Coefficients to Values

Given coefficients $\mathbf{c}$, the values are:
$$
f_j = \sum_{k=0}^{n} c_k \cos\left(\frac{jk\pi}{n}\right)
$$

```python
def coeffs2vals(coeffs):
    """Convert coefficients to values at Chebyshev points."""
    n = len(coeffs) - 1
    if n == 0:
        return coeffs.copy()

    # Undo scaling
    coeffs_scaled = coeffs.copy()
    coeffs_scaled[1:n] /= 2.0

    # Use inverse DCT-I
    values = fft.idct(coeffs_scaled, type=1, norm='forward')
    return values[::-1]
```

## Example: $f(x) = x^3$

```python
n = 4
x = np.cos(np.pi * np.arange(n) / (n-1))  # Chebyshev points
f = x**3                                    # Values

c = vals2coeffs(f)
# c = [0, 0.75, 0, 0.25]
```

The polynomial $x^3$ has the exact Chebyshev expansion:
$$
x^3 = \frac{3}{4}T_1(x) + \frac{1}{4}T_3(x)
$$

:::{dropdown} The Linear System

Let's derive this by setting up and solving the linear system explicitly.

**Step 1: Chebyshev points** ($n = 4$, so $\theta_j = j\pi/3$):

| $j$ | $\theta_j$ | $x_j = \cos\theta_j$ | $f_j = x_j^3$ |
|-----|------------|----------------------|---------------|
| 0 | $0$ | $1$ | $1$ |
| 1 | $\pi/3$ | $1/2$ | $1/8$ |
| 2 | $2\pi/3$ | $-1/2$ | $-1/8$ |
| 3 | $\pi$ | $-1$ | $-1$ |

**Step 2: Build the Chebyshev matrix** $T_{jk} = T_k(x_j)$:

Using $T_0 = 1$, $T_1 = x$, $T_2 = 2x^2 - 1$, $T_3 = 4x^3 - 3x$:

$$
T = \begin{pmatrix}
T_0(1) & T_1(1) & T_2(1) & T_3(1) \\
T_0(\tfrac{1}{2}) & T_1(\tfrac{1}{2}) & T_2(\tfrac{1}{2}) & T_3(\tfrac{1}{2}) \\
T_0(-\tfrac{1}{2}) & T_1(-\tfrac{1}{2}) & T_2(-\tfrac{1}{2}) & T_3(-\tfrac{1}{2}) \\
T_0(-1) & T_1(-1) & T_2(-1) & T_3(-1)
\end{pmatrix}
= \begin{pmatrix}
1 & 1 & 1 & 1 \\
1 & \tfrac{1}{2} & -\tfrac{1}{2} & -1 \\
1 & -\tfrac{1}{2} & -\tfrac{1}{2} & 1 \\
1 & -1 & 1 & -1
\end{pmatrix}
$$

**Step 3: Set up the system** $\mathbf{f} = T\mathbf{c}$:

$$
\begin{pmatrix} 1 \\ 1/8 \\ -1/8 \\ -1 \end{pmatrix}
= \begin{pmatrix}
1 & 1 & 1 & 1 \\
1 & \tfrac{1}{2} & -\tfrac{1}{2} & -1 \\
1 & -\tfrac{1}{2} & -\tfrac{1}{2} & 1 \\
1 & -1 & 1 & -1
\end{pmatrix}
\begin{pmatrix} c_0 \\ c_1 \\ c_2 \\ c_3 \end{pmatrix}
$$

**Step 4: Solve** (or use the DCT, which exploits the structure of $T$):

$$
\mathbf{c} = T^{-1}\mathbf{f} = \begin{pmatrix} 0 \\ 3/4 \\ 0 \\ 1/4 \end{pmatrix}
$$

**Verification:** $T\mathbf{c} = \mathbf{f}$?

- Row 0: $0(1) + \tfrac{3}{4}(1) + 0(1) + \tfrac{1}{4}(1) = 1$ ✓
- Row 1: $0(1) + \tfrac{3}{4}(\tfrac{1}{2}) + 0(-\tfrac{1}{2}) + \tfrac{1}{4}(-1) = \tfrac{3}{8} - \tfrac{1}{4} = \tfrac{1}{8}$ ✓
- Row 2: $0(1) + \tfrac{3}{4}(-\tfrac{1}{2}) + 0(-\tfrac{1}{2}) + \tfrac{1}{4}(1) = -\tfrac{3}{8} + \tfrac{1}{4} = -\tfrac{1}{8}$ ✓
- Row 3: $0(1) + \tfrac{3}{4}(-1) + 0(1) + \tfrac{1}{4}(-1) = -1$ ✓

:::

:::{dropdown} Why the DCT is faster
The matrix $T$ has special structure: $T_{jk} = \cos(jk\pi/(n-1))$. This is exactly the **DCT-I matrix**, and the FFT algorithm exploits this structure to compute $T^{-1}\mathbf{f}$ in $O(n \log n)$ instead of the $O(n^3)$ required for general matrix inversion.
:::

:::{dropdown} Algebraic verification
Using $T_1(x) = x$ and $T_3(x) = 4x^3 - 3x$, we can verify directly:
$$
\frac{3}{4}T_1(x) + \frac{1}{4}T_3(x) = \frac{3}{4}x + \frac{1}{4}(4x^3 - 3x) = \frac{3}{4}x + x^3 - \frac{3}{4}x = x^3 \quad \checkmark
$$
:::

## Example: $f(x) = x^5$

Similarly:
$$
x^5 = \frac{10}{16}T_1(x) + \frac{5}{16}T_3(x) + \frac{1}{16}T_5(x)
$$

The pattern: monomials have sparse Chebyshev representations (only odd or even terms appear).

## The Clenshaw Algorithm

To evaluate $p(x) = \sum_{k=0}^{n} c_k T_k(x)$ without computing each $T_k$, use **Clenshaw's algorithm**:

```python
def clenshaw(x, coeffs):
    """Evaluate Chebyshev series at x using Clenshaw algorithm."""
    n = len(coeffs) - 1
    if n == 0:
        return coeffs[0] * np.ones_like(x)

    b_k1 = np.zeros_like(x)
    b_k2 = np.zeros_like(x)

    for k in range(n, 1, -1):
        b_k2, b_k1 = b_k1, coeffs[k] + 2*x*b_k1 - b_k2

    return coeffs[0] + x*b_k1 - b_k2
```

This is analogous to Horner's method for monomial expansions.

## Why Two Representations?

Different operations are natural in each space:

| Operation | Value Space | Coefficient Space |
|-----------|-------------|-------------------|
| **Sample $f$** | Direct: $f(x_k)$ | Need inverse DCT |
| **Differentiate** | Matrix multiply: $D\mathbf{f}$ | Recurrence relation |
| **Integrate** | Need DCT first | Direct formula |
| **Multiply $f \cdot g$** | Componentwise | Convolution (expensive) |
| **Assess smoothness** | Not easy | Coefficient decay |
| **Truncate** | Not easy | Drop small coefficients |

## The Freedom to Choose

The key insight is that **translation between representations is cheap** ($O(n \log n)$), so we can work in whichever space is most convenient for each operation.

**Example: Computing $\int_{-1}^{1} f(x)^2 \, dx$**

Given the interpolant $p_n(x)$ for $f(x)$:
1. **Square in value space:** $(f_0^2, f_1^2, \ldots, f_n^2)$ — just componentwise!
2. **Transform to coefficient space:** Use DCT
3. **Integrate in coefficient space:** Use the closed-form Chebyshev integral formula

This hybrid approach is often optimal.

## Extension to Infinite Dimensions

For a general Lipschitz continuous function $f(x)$, the Chebyshev series is infinite:
$$
f(x) = \sum_{k=0}^{\infty} c_k T_k(x), \quad f \equiv (c_0, c_1, c_2, \ldots)^T
$$

This "infinite vector" is called a **quasivector** in the literature. While working with infinitely many coefficients requires care (functional analysis!), the practical reality is that:

:::{tip} Key Insight
For functions encountered in applications, the coefficients $c_k$ decay rapidly. Truncating to $n$ terms gives accurate approximations, and we can compute with finite vectors.
:::

## The Chebfun Philosophy

Modern software like **Chebfun** (MATLAB) and **ApproxFun** (Julia) represent functions as their Chebyshev coefficients, automatically:

1. Sampling the function at Chebyshev points
2. Computing coefficients via DCT
3. Adaptively choosing $n$ until coefficients decay to machine precision

This gives "numerical functions" that can be manipulated like symbolic functions but with guaranteed accuracy.

```python
# Conceptual chebfun-style workflow
def chebfun(f, tol=1e-14):
    """Create adaptive Chebyshev approximation."""
    for n in [16, 32, 64, 128, 256, 512, 1024]:
        x = chebpts(n)
        vals = f(x)
        coeffs = vals2coeffs(vals)

        # Check if coefficients have decayed
        if np.max(np.abs(coeffs[-3:])) < tol * np.max(np.abs(coeffs)):
            return coeffs  # Converged!

    raise ValueError("Function too complex or non-smooth")
```

## Complexity Summary

| Operation | Cost |
|-----------|------|
| Values → Coefficients | $O(n \log n)$ via DCT |
| Coefficients → Values | $O(n \log n)$ via inverse DCT |
| Evaluate at one point | $O(n)$ via Clenshaw or barycentric |
| Differentiate (value space) | $O(n^2)$ matrix-vector product |
| Integrate (coefficient space) | $O(n)$ |

The DCT/FFT connection makes transforming between representations cheap, so we can use whichever is more convenient for each operation.

## Aliasing: The Sampling Pitfall

When we sample a function at $N$ discrete points, we **cannot distinguish** high-frequency components from low-frequency ones. This is **aliasing**—the bane of spectral methods.

### The Nyquist Limit

:::{prf:theorem} Nyquist-Shannon Sampling Theorem
:label: thm-nyquist

To uniquely recover a signal with frequencies up to $f_{\max}$, we must sample at rate:
$$
f_s \geq 2 f_{\max}
$$

Equivalently: with $N$ sample points, we can only resolve frequencies $k \leq N/2$.
:::

For polynomial interpolation at $N+1$ points: we can represent polynomials of degree at most $N$. Anything higher gets **aliased** to lower frequencies.

### Aliasing in Action

Consider sampling $\cos(10\theta)$ at $N = 8$ points ($\theta_j = j\pi/8$):

```python
N = 8
theta = np.pi * np.arange(N) / N
f_high = np.cos(10 * theta)  # Frequency 10 > N/2 = 4

# These samples are IDENTICAL to cos(6*theta)!
f_low = np.cos(6 * theta)
np.allclose(f_high, f_low)  # True!
```

The frequency-10 wave "folds back" and appears as frequency 6. We can't tell them apart from samples alone.

:::{prf:remark} Aliasing Formula
:label: rmk-aliasing-formula

On a grid of $N$ points, frequency $k$ aliases to frequency $k \mod N$ (mapped to $[-N/2, N/2]$). For Chebyshev grids, the folding is around $N$: frequency $k$ aliases to $|2N - k|$ when $k > N$.
:::

### Why Aliasing Matters

Aliasing causes **silent errors**:

1. **No warning:** The computed coefficients look reasonable
2. **Wrong answer:** High-frequency content corrupts low-frequency coefficients
3. **Hard to detect:** Everything "works" until you compare with a finer grid

### The 2/3 Rule for Nonlinear Terms

When computing products like $u \cdot v$ in spectral methods, aliasing becomes critical.

**The problem:** If $u$ and $v$ each have frequencies up to $N/2$, then $u \cdot v$ has frequencies up to $N$—exceeding our resolution!

:::{prf:proposition} The 2/3 Dealiasing Rule
:label: prop-two-thirds-rule

To compute $u \cdot v$ without aliasing on an $N$-point grid:
1. **Zero-pad:** Extend both signals to $M = 3N/2$ points
2. **Transform:** Use inverse FFT to get values on the finer grid
3. **Multiply:** Compute pointwise product in physical space
4. **Transform back:** Use FFT and truncate to $N$ coefficients

Alternatively: keep only the lowest $2N/3$ modes of $u$ and $v$ before multiplying.
:::

```python
def dealias_product(u_hat, v_hat):
    """Compute u*v with 2/3 dealiasing."""
    N = len(u_hat)
    M = 3 * N // 2  # Padding factor

    # Zero-pad to finer grid
    u_pad = np.zeros(M, dtype=complex)
    v_pad = np.zeros(M, dtype=complex)
    u_pad[:N//2] = u_hat[:N//2]
    u_pad[-(N//2):] = u_hat[-(N//2):]
    v_pad[:N//2] = v_hat[:N//2]
    v_pad[-(N//2):] = v_hat[-(N//2):]

    # Multiply in physical space
    uv = np.fft.ifft(u_pad) * np.fft.ifft(v_pad)

    # Transform back and truncate
    uv_hat = np.fft.fft(uv)
    result = np.zeros(N, dtype=complex)
    result[:N//2] = uv_hat[:N//2]
    result[-(N//2):] = uv_hat[-(N//2):]
    return result * (M / N)  # Normalization
```

:::{warning}
**Ignoring aliasing in nonlinear problems leads to instability!** The classic example: solving Burgers' equation $u_t + uu_x = 0$ without dealiasing produces spurious oscillations that grow unboundedly.
:::

### Practical Guidelines

| Situation | Recommendation |
|-----------|----------------|
| Linear problems | Aliasing less critical—errors stay bounded |
| Nonlinear products | Always dealias (2/3 rule or zero-padding) |
| Function evaluation | Check convergence by comparing $N$ and $2N$ |
| Unknown smoothness | Monitor coefficient decay |

### Detecting Aliasing

How to know if your solution is aliased:

1. **Coefficient plateau:** If $|c_k|$ doesn't decay to machine precision, you may need more points
2. **Resolution test:** Double $N$ and compare—if answer changes significantly, you were underresolved
3. **Energy in high modes:** If significant energy near $k = N/2$, aliasing is likely

```python
def check_resolution(coeffs, tol=1e-10):
    """Warn if coefficients suggest underresolution."""
    if np.max(np.abs(coeffs[-5:])) > tol * np.max(np.abs(coeffs)):
        print("Warning: coefficients not fully resolved—aliasing possible!")
```

### The Moral

> *"Aliasing is the price of discretization."*

When sampling continuous functions at discrete points, information is lost. Nyquist tells us exactly how much bandwidth we can capture. For spectral methods to work reliably:

- **Linear problems:** Ensure solution is smooth enough for the grid
- **Nonlinear problems:** Dealias religiously
- **Unknown problems:** Always check convergence
