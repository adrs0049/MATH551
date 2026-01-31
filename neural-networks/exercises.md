# Exercises

Practice with neural network approximation theory.

---

### Q13.1: Simple Neural Network

Consider a single-neuron network: $f(x) = a \cdot \sigma(wx + b)$ where $\sigma(t) = \tanh(t)$.

- **(a)** For $a = 1, w = 1, b = 0$, plot $f(x)$ on $[-5, 5]$.
- **(b)** For $a = 1, w = 10, b = 0$, plot $f(x)$. What happens as $w \to \infty$?
- **(c)** Explain how adjusting $w$ and $b$ shifts and scales the transition region.

---

### Q13.2: Approximating a Step Function

Use a sum of sigmoid neurons to approximate the step function $f(x) = \mathbf{1}_{x > 0}$.

- **(a)** Use a single neuron $\sigma(wx)$ with large $w$. How close can you get?
- **(b)** Use the difference $\sigma(w(x - a)) - \sigma(w(x - b))$ to create a "bump." Plot for various $w$.
- **(c)** Add several such bumps to approximate a piecewise constant function.

---

### Q13.3: ReLU Networks

For ReLU activation $\sigma(t) = \max(0, t)$:

- **(a)** Show that $|t| = \sigma(t) + \sigma(-t)$.
- **(b)** How can you build a "tent function" using two ReLU neurons?
- **(c)** Use 4 ReLU neurons to approximate a function that is zero outside $[-1, 1]$ and has a peak at $x = 0$.

---

### Q13.4: Universal Approximation in Action

Using `numpy` and optimization (e.g., `scipy.optimize.minimize`), fit a neural network with $n = 5$ sigmoid neurons to:

- **(a)** $f(x) = \sin(2\pi x)$ on $[0, 1]$
- **(b)** $f(x) = |x|$ on $[-1, 1]$
- **(c)** Compare errors for $n = 5, 10, 20$ neurons.

---

### Q13.5: Dimension Dependence

Consider approximating $f(x_1, \ldots, x_d) = \exp(-\|x\|^2)$ on $[-1, 1]^d$.

- **(a)** For polynomial approximation of degree $n$, how many terms are needed?
- **(b)** A theorem says Gaussians have finite Barron norm. What does this predict about neural network approximation?
- **(c)** Why can't we numerically verify this for $d = 100$?

---

### Q13.6: Barron Norm Computation

For the 1D function $f(x) = e^{-x^2/2}$:

- **(a)** Compute the Fourier transform $\hat{f}(\omega)$.
- **(b)** Verify that $C_f = \int |\omega||\hat{f}(\omega)|d\omega < \infty$.
- **(c)** What does Barron's theorem predict for the approximation rate?

---

### Q13.7: Comparing Approximation Methods

For $f(x) = \frac{1}{1 + 25x^2}$ on $[-1, 1]$:

- **(a)** Approximate with degree-$n$ polynomial using Chebyshev nodes. Plot error vs. $n$.
- **(b)** Approximate with a neural network with $n$ neurons (use optimization). Plot error vs. $n$.
- **(c)** Which achieves 6-digit accuracy with fewer parameters?

---

### Q13.8: Deep vs. Shallow

Consider $f(x) = x^{2^L}$ for $L = 4$ (so $f(x) = x^{16}$).

- **(a)** A single hidden layer network needs $O(16)$ neurons to represent this exactly with ReLU. Why?
- **(b)** A deep network with $L$ layers needs $O(L) = O(4)$ neurons. Explain by building it layer by layer: $x \to x^2 \to x^4 \to x^8 \to x^{16}$.
- **(c)** This is an example where depth helps exponentially. Can you think of other such examples?

---

### Q13.9: Non-Approximable Functions

Consider $f(x) = \sin(1000x)$ on $[0, 1]$.

- **(a)** How many polynomial terms (Chebyshev) are needed to capture this oscillation?
- **(b)** How does the Barron norm scale with frequency?
- **(c)** Is this function "hard" for neural networks? Why or why not?

---

### Q13.10: Neural Networks for PDEs

Consider the ODE $u'' = -\pi^2 u$ with $u(0) = 0, u(1) = 0$.

- **(a)** Write a loss function that measures how well a neural network $u_\theta(x)$ satisfies the ODE and boundary conditions.
- **(b)** Implement this (known as a "physics-informed neural network" or PINN).
- **(c)** Compare accuracy to spectral collocation. Which is more efficient for this smooth problem?

---

## Self-Assessment Questions

Test your understanding with these conceptual questions:

1. **Universal Approximation:** What does the universal approximation theorem say? What does it NOT say?

2. **Activation Functions:** Why do we need nonlinear activation functions? What happens with $\sigma(t) = t$?

3. **Curse of Dimensionality:** For polynomial approximation in $d$ dimensions with degree $n$, how many terms are needed?

4. **Barron Norm:** What is the Barron norm? What property of a function does it measure?

5. **Dimension Independence:** Why is Barron's $O(1/\sqrt{n})$ rate remarkable compared to polynomial rates?

6. **Existence vs. Computation:** Barron's theorem proves good approximations exist. Why doesn't this immediately solve machine learning?

7. **Deep vs. Shallow:** Give an example where depth helps. Does Barron's theorem apply to deep networks?
