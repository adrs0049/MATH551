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

The concepts from our error analysis carry over directly:

| Concept | General Setting | ODE Numerics |
|---------|-----------------|--------------|
| **Backward error** | Residual $\|f(\tilde{x}) - b\|$ | Local truncation error $\tau_n$ |
| **Forward error** | $\|\tilde{x} - x\|$ | Global error $\|u_n - u(t_n)\|$ |
| **Condition number** | Sensitivity to perturbations | Lipschitz constant $L$ of $f$ |
| **Stability** | Backward stable algorithm | Absolutely stable method |

Every numerical method introduces **local truncation error** at each step—this is the backward error. The central question: does this error accumulate controllably, or does it explode?

$$
\text{Global error} \lesssim \text{(Amplification factor)}^n \times \text{Local truncation error}
$$

The amplification factor depends on both the **method** (its stability properties) and the **problem** (the Lipschitz constant). A method is useful when errors remain bounded—this is the **stability** requirement, the ODE analog of backward stability.

For **stiff problems**, where eigenvalues span many orders of magnitude, explicit methods require impractically small step sizes for stability even when accuracy would permit larger steps. Implicit methods—the "backward stable" algorithms of ODE numerics—handle such problems gracefully.

The **stiffness ratio** $R = \max|\lambda_i| / \min|\lambda_i|$ can be understood as a condition number for method choice: it measures how ill-conditioned the "use an explicit method" approach is. When $R$ is large, the ratio of work required for stability versus work required for accuracy becomes enormous—this is the computational signature of stiffness.

## Learning Outcomes

After completing this chapter, you should be able to (THIS NEEDS UPDATES):

- **L5.1:** Derive and implement forward and backward Euler methods
- **L5.2:** Compute local truncation error and determine consistency order
- **L5.3:** Prove convergence of one-step methods using Lipschitz bounds
- **L5.4:** Determine stability regions and identify A-stable methods
- **L5.5:** Recognize stiff problems and explain why they require implicit methods
- **L5.6:** Derive RK2 methods by matching Taylor series coefficients
- **L5.7:** Implement adaptive time-stepping using embedded pairs
