# Numerical Methods for ODEs

:::{tip} Big Idea
Solving ODEs numerically is **time-stepping**: we repeatedly approximate the future from the present. The fundamental trade-off is between **explicit methods** (cheap but stability-restricted) and **implicit methods** (expensive but unconditionally stable). For **stiff problems**—where stability forces impractically small step sizes—implicit methods are essential.
:::

## Discretization

We consider the initial value problem:

$$
\frac{du}{dt} = f(t, u), \quad u(t_0) = u_0
$$

where $u$ may be a vector. We assume $f$ is Lipschitz continuous, guaranteeing a unique solution.

The solution $u(t)$ is a continuous function, but computers work with discrete data. We discretize the time interval $[t_0, t_e]$ into a **lattice**:

$$
t_0 < t_1 < t_2 < \cdots < t_N = t_e
$$

with step sizes $h_n = t_{n+1} - t_n$. For simplicity, we often use uniform spacing $h_n = h$.

A **discretization method** associates to each lattice a **lattice function**:

$$
u_n \approx u(t_n), \quad n = 0, 1, \ldots, N
$$

We store only the values $\{u_0, u_1, \ldots, u_N\}$—a finite amount of data representing the continuous solution.

The simplest methods replace derivatives with finite differences:

$$
\frac{du}{dt} \approx \frac{u_{n+1} - u_n}{h}
$$

This immediately gives **forward Euler**: $u_{n+1} = u_n + h f(t_n, u_n)$.

## The Error Framework for ODEs

The error analysis for ODE solvers mirrors the structure we have seen throughout the course:

| Concept | ODE Numerics |
|---------|--------------|
| **Per-step accuracy** | Local truncation error $\tau_n = O(h^p)$ |
| **Global error** | $E_n = u_n - u(t_n)$ |
| **Problem sensitivity** | Condition number $e^{LT}$ (from Lipschitz constant $L$) |
| **Stability** | Amplification factor $\lvert 1 + h\lambda \rvert \leq 1$ |

Every numerical method introduces **local truncation error** at each step. The central question: does this error accumulate controllably, or does it explode?

$$
\text{Global error} \;\leq\; \underbrace{e^{LT}}_{\text{problem sensitivity}} \;\times\; \underbrace{O(h^p)}_{\text{per-step accuracy}}
$$

The factor $e^{LT}$ depends on the **problem** (its Lipschitz constant $L$ and integration time $T$). Whether the method amplifies or damps these errors depends on the **amplification factor**, which must stay bounded for the method to be useful. This is the **stability** requirement.

For **stiff problems**, where $|\lambda|$ is large, explicit methods require impractically small step sizes for stability even when accuracy would permit larger steps. Implicit methods handle such problems gracefully.

## Learning Outcomes

After completing this chapter, you should be able to:

- **L5.1:** Derive forward Euler from the forward difference approximation and interpret it geometrically
- **L5.2:** Compute the local truncation error $\tau_n$ and distinguish it from the one-step error $h\tau_n$
- **L5.3:** State the Lipschitz condition and explain its role in both existence/uniqueness and convergence
- **L5.4:** Derive the condition number $e^{LT}$ of an initial value problem from Gronwall's inequality
- **L5.5:** State and interpret the convergence theorem: consistency of order $p$ implies convergence of order $p$
- **L5.6:** Determine the stability region of forward Euler and compute the step-size restriction $h \leq 2/|\lambda|$
- **L5.7:** Recognize stiff problems and explain why they require implicit methods or adaptive step sizes
