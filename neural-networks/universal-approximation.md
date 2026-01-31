# Universal Approximation Theorem

:::{tip} Big Idea
Neural networks with a single hidden layer can approximate any continuous function on a compact set to arbitrary accuracy. This is a qualitative existence result—it doesn't tell us how many neurons we need.
:::

## What is a Neural Network?

A **single hidden layer neural network** computes:

$$
f_n(x) = \sum_{k=1}^{n} a_k \sigma(w_k \cdot x + b_k)
$$

where:
- $x \in \mathbb{R}^d$ is the input
- $w_k \in \mathbb{R}^d$ are weight vectors
- $b_k \in \mathbb{R}$ are biases
- $a_k \in \mathbb{R}$ are output weights
- $\sigma: \mathbb{R} \to \mathbb{R}$ is an **activation function**
- $n$ is the number of neurons (width)

### Common Activation Functions

| Name | Formula | Properties |
|------|---------|------------|
| Sigmoid | $\sigma(t) = \frac{1}{1 + e^{-t}}$ | Smooth, bounded |
| Tanh | $\sigma(t) = \tanh(t)$ | Smooth, bounded, zero-centered |
| ReLU | $\sigma(t) = \max(0, t)$ | Simple, unbounded |
| Softplus | $\sigma(t) = \log(1 + e^t)$ | Smooth approximation to ReLU |

## The Theorem

:::{prf:theorem} Universal Approximation Theorem (Cybenko 1989, Hornik 1991)
:label: thm-universal-approximation

Let $\sigma: \mathbb{R} \to \mathbb{R}$ be a non-constant, bounded, continuous function. Let $K \subset \mathbb{R}^d$ be compact.

For any $f \in C(K)$ and any $\epsilon > 0$, there exists a neural network $f_n$ such that:
$$
\|f - f_n\|_\infty < \epsilon
$$
:::

**In words:** Neural networks are **dense** in $C(K)$—just like polynomials (Weierstrass) and trigonometric functions.

## Connection to the Weierstrass Theorem

The Universal Approximation Theorem is a modern descendant of one of the most beautiful results in classical analysis.

:::{prf:theorem} Weierstrass Approximation Theorem (1885)
:label: thm-weierstrass

Let $f: [a,b] \to \mathbb{R}$ be continuous. For any $\epsilon > 0$, there exists a polynomial $p$ such that:
$$
\|f - p\|_\infty < \epsilon
$$
:::

The structural parallel is striking:

| Weierstrass (1885) | Universal Approximation (1989) |
|--------------------|-------------------------------|
| Polynomials $\sum a_k x^k$ | Neural networks $\sum a_k \sigma(w_k x + b_k)$ |
| Monomials $\{1, x, x^2, \ldots\}$ as basis | Ridge functions $\{\sigma(w \cdot x + b)\}$ as basis |
| Dense in $C([a,b])$ | Dense in $C(K)$ |
| Existence only—no degree bound | Existence only—no width bound |

Both theorems answer the same fundamental question: **Can this family of functions approximate anything?** Both say "yes" without telling us how many terms we need.

### The Crucial Caveat: Existence vs. Convergence of a Fixed Basis

As Trefethen emphasizes, there's an important subtlety here. The Weierstrass theorem says: for each $f$ and $\epsilon$, there *exists* some polynomial that works. But this polynomial depends on $f$—it's not saying that expanding $f$ in a fixed basis (like Chebyshev or Fourier) will converge.

For a **fixed orthogonal basis** to converge, we need more regularity:

| Basis | Convergence requirement |
|-------|------------------------|
| Fourier series | Bounded variation for pointwise convergence; smoothness for rates |
| Chebyshev expansion | Bounded variation; Lipschitz or better for good rates |
| Legendre expansion | Similar regularity requirements |

**The Weierstrass theorem allows us to tailor the polynomial to the function.** This is analogous to how the Universal Approximation Theorem lets us choose the weights $w_k, b_k, a_k$ to match the target function.

:::{note}
The parallel deepens: just as Chebyshev/Fourier require bounded variation or smoothness for guaranteed convergence rates, Barron's theorem (next section) shows that neural networks achieve dimension-independent rates only when the target function has bounded "Barron norm"—a spectral regularity condition.
:::

**Key difference in high dimensions:** In $d$ dimensions, polynomials of degree $n$ have $\binom{n+d}{d}$ terms—exponential in $d$. Neural networks with the right activation can sometimes escape this curse for functions with appropriate spectral structure, which is why the Universal Approximation Theorem opened a new chapter in approximation theory.

## Proof Sketch

The key insight: ridge functions $\sigma(w \cdot x + b)$ form a "rich enough" set.

### For Sigmoid Activation

1. As the weight magnitude $\|w\|$ increases, $\sigma(w \cdot x + b)$ approaches a step function
2. Step functions can approximate indicator functions of half-spaces
3. Indicator functions of half-spaces can approximate any continuous function (by chopping the domain into pieces)

### For General Activations

The proof uses:
- Hahn-Banach theorem (functional analysis)
- The fact that $\sigma$ is non-polynomial (or bounded and non-constant)
- Density arguments in $C(K)$

## What the Theorem Does NOT Say

:::{warning}
The universal approximation theorem is an **existence** result. It does not tell us:

1. **How many neurons** are needed for a given accuracy
2. **How to find** the weights (training is a separate problem)
3. **Which functions** are easy vs. hard to approximate
4. **How depth** (multiple layers) helps
:::

This is analogous to:
- Weierstrass says polynomials are dense, but doesn't say which degree you need
- Stone-Weierstrass is qualitative, not quantitative

## The Curse of Dimensionality

For polynomials of degree $n$ in $d$ dimensions, we have $\binom{n+d}{d}$ terms.

For $d = 10$ and $n = 10$: approximately $10^6$ terms!

The universal approximation theorem doesn't save us from this. Generic continuous functions still require exponentially many neurons as dimension grows.

**The key question:** Are there function classes where neural networks beat this exponential barrier?

**Answer:** Yes! These are the **Barron spaces**.

## Extensions

### ReLU Networks

The original theorem required bounded activation. Extensions show:

:::{prf:theorem} Universal Approximation for ReLU
:label: thm-relu-universal

Networks with ReLU activation $\sigma(t) = \max(0, t)$ are also universal approximators, despite being unbounded and non-smooth.
:::

### Deep Networks

The universal approximation theorem applies to single hidden layer networks. What about deeper networks?

**Depth helps expressivity:** Some functions can be represented exactly with $O(\log n)$ neurons in a deep network, but require exponentially many neurons in a shallow network.

Example: The function $f(x) = x^{2^L}$ can be computed by:
- Shallow network: needs $O(2^L)$ neurons
- Deep network with $L$ layers: needs $O(L)$ neurons

## Comparison with Classical Approximation

| Aspect | Polynomials | Fourier | Neural Networks |
|--------|-------------|---------|-----------------|
| Universal in $C(K)$? | Yes | Yes (periodic) | Yes |
| Constructive? | Yes (Lagrange, etc.) | Yes (FFT) | Sort of (training) |
| Convergence rate? | Known | Known | Depends on function class |
| High-$d$? | Curse | Curse | Can escape for Barron |

## Summary

The universal approximation theorem tells us:
- ✅ Neural networks **can** approximate any continuous function
- ❌ It doesn't say **how** or **how efficiently**
- 🔑 The real power comes from Barron's theorem (next section)

The theorem is like knowing "there exists a polynomial approximation"—useful, but we need quantitative bounds to do practical computation.
