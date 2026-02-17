---
authors:
  - name: Thejani Gamage
---

# Numerical Simulation of Stochastic Differential Equations

:::{tip} Big Idea
Ordinary differential equations (ODEs) describe the evolution of systems when
the dynamics are deterministic. In many real-world applications — such as finance,
physics, and biology — systems are influenced by **random noise**. Stochastic
differential equations (SDEs) model the evolution of such systems.

The Euler–Maruyama method extends forward Euler to SDEs: at each step we add
a deterministic contribution (the *drift*) and a random contribution (the
*diffusion*). Each simulation produces one possible **sample path**, and
convergence is measured in expectation rather than pointwise.
:::

## Why Stochastic?

Many systems we want to model are not purely deterministic:

- **Finance:** Stock prices fluctuate due to unpredictable market forces. The
  Black–Scholes model for asset prices is an SDE — geometric Brownian motion.
- **Biology:** Populations are subject to environmental noise. Deterministic
  logistic growth $dX = rX(1 - X/K)\,dt$ becomes stochastic when we add
  $\sigma X\,dW$ to model random environmental fluctuations.
- **Physics:** A particle suspended in fluid is buffeted by random molecular
  collisions — the original Brownian motion observed by Robert Brown in 1827.

In each case, a deterministic ODE captures the average behavior, but misses
the inherent randomness. SDEs provide the mathematical framework to include it.

## Stochastic Differential Equation (SDE)

A simple SDE has the form
$$dX(t) = a(t,X(t))\,dt + b(t,X(t))\,dW(t),$$
$$X(0)=x$$
where:
- $a(t,x)$ is the **drift** (represents deterministic changes),
- $b(t,x)$ is the **diffusion** (represents random changes),
- $W(t)$ is **Brownian motion** (Wiener process).
- $x$ is the initial point (usually known).

---

## Brownian Motion

Brownian motion $W(t)$ is a stochastic process. For each time $t$, $W(t)$ is a random variable and it satisfies:

1. $W(0) = 0$
2. Increments are independent; $W(t_3)-W(t_2)$ is independent from  $W(t_1)-W(t_0)$ for $0 \leq t_0 < t_1 \leq t_2 < t_3$
3. $W(t+h)-W(t) \sim \mathcal{N}(0,h)$
4. $W(t)$ is a continuous function of $t$, with probability $1$.

---

## Why Euler’s Method Needs Modification

Formally integrating the SDE over $[t,t+h]$ we obtain:

$$
X(t+h) = X(t)
+ \int_t^{t+h} a(s,X(s))\,ds
+ \int_t^{t+h} b(s,X(s))\,dW(s).
$$

We approximate:
- the deterministic integral as in Euler’s method,
- the stochastic integral using properties of Brownian motion.

---

## Euler–Maruyama Method

Let $t_n = nh$. Define
$$
\Delta W_n = W(t_{n+1}) - W(t_n).
$$

### Scheme

$$
X_{n+1}
=
X_n
+ a(t_n,X_n)\,h
+ b(t_n,X_n)\,\Delta W_n
$$

where
$$
\Delta W_n \sim \mathcal{N}(0,h).
$$
and $$X_0=x.$$

This is the **Euler–Maruyama method**.

## Convergence Analysis

Brownian paths are:
- nowhere differentiable,
- extremely irregular.

Therefore:
- pointwise convergence is too strong,
- convergence in terms of expected values are the natural quantities to study.
- Convergence guarantees depend on the properties of the drift and diffusion function


---

## Example: Geometric Brownian Motion

$$
dX_t = a X_t\,dt + b X_t\,dW_t.
$$
where a,b are constants.
Euler–Maruyama gives

$$
X_{n+1}
=
X_n
+ aX_n h
+ b X_n \Delta W_n.$$

where
$$
\Delta W_n \sim \mathcal{N}(0,h).
$$
and
$$
X_0=x.
$$


This is the SDE underlying the Black–Scholes model for asset prices in financial mathematics.

The exact solution is known:

$$X(t) = x \exp\!\left[\left(a - \tfrac{b^2}{2}\right)t + b\, W(t)\right].$$

This is one of the rare SDEs with a closed-form solution, making it an ideal
test problem: we can compare the Euler–Maruyama approximation against the true
solution evaluated at the same Brownian path. Most SDEs encountered in practice
cannot be solved explicitly.

## Summary

- Euler–Maruyama extends Euler’s method to SDEs.
- Random increments satisfy
  $$
  \mathbb{E}[\Delta W]=0,\quad \mathbb{E}[(\Delta W)^2]=h.
  $$
- Numerical solutions are **random variables**.
- Accuracy is measured in **expectation**, not pointwise.

:::{seealso}
For a complete algorithmic treatment with code examples and convergence
experiments, see {cite}`Higham2001`.
:::
