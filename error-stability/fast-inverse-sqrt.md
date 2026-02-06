# Fast Inverse Square Root

:::{tip} Big Idea
The famous fast inverse square root algorithm from Quake III demonstrates how understanding floating point representation enables creative numerical tricks. It combines bit-level manipulation with Newton's method for refinement.
:::

## Motivation

Computing unit normal vectors is essential for:
- Light reflections in games
- Collision detection
- Computer-aided design

Given a vector $\mathbf{v} = (v_1, v_2, v_3)$, its unit vector is:

$$
\hat{\mathbf{v}} = \frac{\mathbf{v}}{\sqrt{v_1^2 + v_2^2 + v_3^2}}
$$

This requires evaluating the inverse square root: $f(x) = \frac{1}{\sqrt{x}}$

In the 1990s when computational power was limited, a fast implementation was crucial for real-time graphics.

## The Famous Code

From the [Quake III Arena source code](https://github.com/id-Software/Quake-III-Arena/blob/dbe4ddb10315479fc00086f08e25d968b4b43c49/code/game/q_math.c#L552):

```c
float Q_rsqrt( float x )
{
    long i;
    float x2, y;
    const float threehalfs = 1.5F;

    x2 = x * 0.5F;
    y  = x;
    i  = * ( long * ) &y;               // evil floating point bit level hacking
    i  = 0x5f3759df - ( i >> 1 );       // what the fuck?
    y  = * ( float * ) &i;
    y  = y * ( threehalfs - ( x2 * y * y ) );   // 1st iteration
    return y;
}
```

Let's understand what this does.

## The Strategy: Turn Multiplication into Addition

We want to compute $y = x^{-1/2}$. Taking $\log_2$ of both sides:

$$
\log_2(y) = -\tfrac{1}{2}\log_2(x)
$$

If we had a cheap way to compute $\log_2$, we could replace the expensive power $x^{-1/2}$ with a simple multiply-by-$(-1/2)$, then exponentiate back. Of course, computing $\log_2$ is itself expensive — unless the hardware is already doing it for us.

## Integer Interpretation of Float Bits

Recall that a 32-bit float stores sign bit $S_x$, exponent $E_x$, and mantissa bits $m_x$. If we read those same 32 bits as an integer, we get:

$$
\text{Int}(x) = 2^{31} S_x + 2^{23} E_x + M_x
$$

where $M_x = 2^{23} m_x$ is the integer value of the mantissa bits. This is what the C line `i = *(long *)&y` does — it reinterprets the float's bits as an integer without changing them.

## Why Integer Bits $\approx$ Logarithm

From the [IEEE 754 representation](floating-point.md), a positive float represents:

$$
x = 2^{E_x - 127}(1 + m_x)
$$

Taking $\log_2$:

$$
\log_2(x) = (E_x - 127) + \log_2(1 + m_x)
$$

Since $m_x \in [0, 1)$, we approximate $\log_2(1 + m_x) \approx m_x + \sigma$ where $\sigma \approx 0.0430$. Substituting $m_x = M_x / 2^{23}$:

$$
\log_2(x) \approx \frac{1}{2^{23}}\underbrace{(M_x + 2^{23} E_x)}_{\text{Int}(x)} + (\sigma - 127)
$$

That is:

$$
\log_2(x) \approx \frac{1}{2^{23}} \cdot \text{Int}(x) + (\sigma - 127)
$$

**The integer interpretation of float bits is a scaled-and-shifted logarithm.** This is the key insight that makes the entire algorithm work.

## The Magic Formula

Now apply this to $\log_2(y) = -\frac{1}{2}\log_2(x)$. Writing both sides in terms of integer bit patterns:

$$
\frac{1}{2^{23}} \text{Int}(y) + (\sigma - 127) = -\frac{1}{2}\left[\frac{1}{2^{23}} \text{Int}(x) + (\sigma - 127)\right]
$$

Solving for $\text{Int}(y)$:

$$
\text{Int}(y) = \underbrace{\tfrac{3}{2} \cdot 2^{23}(127 - \sigma)}_{\texttt{0x5f3759df}} - \frac{1}{2}\,\text{Int}(x)
$$

This is exactly the "what the fuck?" line:

```c
i = 0x5f3759df - (i >> 1);
```

- `i >> 1`: integer right shift divides $\text{Int}(x)$ by 2
- Subtract from the magic constant: completes the formula

The magic number `0x5f3759df` is just $\frac{3}{2} \cdot 2^{23}(127 - \sigma)$ evaluated for the optimal $\sigma$.

## Newton's Method Refinement

The bit manipulation gives a rough approximation. To improve accuracy, we apply Newton's method to:

$$
g(y) = \frac{1}{y^2} - x
$$

The root of $g$ is $y = \frac{1}{\sqrt{x}}$.

Newton's iteration:

$$
y_{n+1} = y_n - \frac{g(y_n)}{g'(y_n)} = y_n\left(\frac{3}{2} - \frac{x}{2}y_n^2\right)
$$

This is exactly line 12 of the code:
```c
y = y * (threehalfs - (x2 * y * y));
```

## Why It Works

1. **Bit manipulation** exploits the logarithmic relationship between a float's value and its bit pattern to get a rough initial guess
2. **Newton's method** rapidly refines this guess (one iteration often suffices for graphics applications)

The algorithm is remarkably accurate while using only:
- One multiplication
- One subtraction
- One bit shift
- One Newton iteration

No expensive division or square root operations are needed.

## Modern Relevance

Modern CPUs have dedicated instructions for inverse square root (e.g., `rsqrtss` on x86), making this trick less necessary.

