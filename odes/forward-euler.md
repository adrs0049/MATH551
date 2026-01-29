# Forward Euler Method

:::{tip} Big Idea
Forward Euler is the simplest time-stepping method: replace the derivative with a forward difference. It's **explicit**—we compute $u_{n+1}$ directly from $u_n$—making it cheap per step. But this simplicity comes at a cost: forward Euler is only **conditionally stable**, requiring step sizes small enough that $|1 + h\lambda| \leq 1$. For stiff problems, this constraint becomes impractical.
:::

## The Initial Value Problem

We seek to approximate the solution of

$$
\frac{du}{dt} = f(t, u), \quad u(t_0) = u_0
$$

on a grid of times $t_0 < t_1 < t_2 < \cdots < t_N = T$ with step size $h = t_{n+1} - t_n$.

Let $u_n$ denote our approximation to $u(t_n)$.

## Derivation

Replace $du/dt$ with a forward difference:

$$
\frac{u(t_{n+1}) - u(t_n)}{h} \approx f(t_n, u(t_n))
$$

This gives the **forward Euler** (or explicit Euler) method:

:::{prf:definition} Forward Euler Method
:label: def-forward-euler

$$
u_{n+1} = u_n + h f(t_n, u_n)
$$

Starting from $u_0$, we can march forward in time: given $u_n$, we compute $u_{n+1}$ directly.
:::

### Geometric Interpretation

Forward Euler follows the tangent line:
1. At $(t_n, u_n)$, the slope is $f(t_n, u_n)$
2. Follow this slope for time $h$
3. Arrive at $u_{n+1} = u_n + h \cdot \text{slope}$

Each step follows a tangent to a (possibly different) solution curve.

### Implementation

```python
def forward_euler(f, u0, t0, T, h):
    """
    Solve u' = f(t, u) from t0 to T with step size h.
    """
    t = t0
    u = u0
    times = [t]
    values = [u]

    while t < T:
        u = u + h * f(t, u)
        t = t + h
        times.append(t)
        values.append(u)

    return np.array(times), np.array(values)
```

## Local Truncation Error

The **local truncation error** measures how well the exact solution satisfies the numerical scheme.

:::{prf:definition} Local Truncation Error
:label: def-local-truncation-error

Substitute the exact solution $u(t)$ into the finite difference scheme. The local truncation error $\tau_n$ is defined by:

$$
\frac{u(t_{n+1}) - u(t_n)}{h} = f(t_n, u(t_n)) + \tau_n
$$
:::

For forward Euler, Taylor expand the exact solution:

$$
u(t_{n+1}) = u(t_n) + h u'(t_n) + \frac{h^2}{2}u''(t_n) + O(h^3)
$$

Since $u'(t_n) = f(t_n, u(t_n))$:

$$
\tau_n = \frac{h}{2}u''(t_n) + O(h^2) = O(h)
$$

:::{prf:definition} Consistency Order
:label: def-consistency-order

A single-step method has **consistency order** $p$ if for sufficiently smooth $f$:

$$
\|\tau(t + h)\| \leq Ch^{p} \quad \text{for all } h \in (0, h_0]
$$

Forward Euler has consistency order $p = 1$.
:::

:::{prf:remark} Local Truncation Error as Backward Error
:label: rmk-lte-backward-error

The local truncation error is the **backward error** of one time step. When we compute $u_{n+1}$ from $u_n$, we're solving a *perturbed* problem:

$$
\frac{u(t_{n+1}) - u(t_n)}{h} = f(t_n, u(t_n)) + \tau_n
$$

The residual $\tau_n$ measures how much we've perturbed the differential equation. A method with small $\tau_n$ solves a problem *close* to the original—this is backward stability in the ODE setting.
:::

## Global Error and Convergence

The **global error** $E_n = u_n - u(t_n)$ accumulates local errors over all steps.

### Error Propagation: Building Intuition

To understand how local errors accumulate, consider the linear test problem $u' = \lambda u + g(t)$. The exact solution satisfies:

$$
u(t_{n+1}) = (1 + h\lambda)u(t_n) + hg(t_n) + h\tau_n
$$

The numerical solution satisfies:

$$
u_{n+1} = (1 + h\lambda)u_n + hg(t_n)
$$

Subtracting gives the error recurrence:

$$
E_{n+1} = (1 + h\lambda)E_n - h\tau_n
$$

Unrolling this recurrence:

$$
E_n = (1 + h\lambda)^n E_0 - h\sum_{j=0}^{n-1}(1 + h\lambda)^{n-1-j}\tau_j
$$

This is the **discrete Duhamel principle**: the global error at step $n$ equals the initial error propagated forward, plus all local truncation errors propagated by the appropriate powers of the amplification factor.

:::{prf:remark} Connection to Neumann Series
:label: rmk-duhamel-neumann
:class: dropdown

This structure is identical to the **Neumann series** from iterative methods! Compare:

| Context | Recurrence | Solution |
|---------|------------|----------|
| Iterative methods | $x_{k+1} = Ax_k + b$ | $x_n = A^n x_0 + \sum_{j=0}^{n-1} A^{n-1-j}b$ |
| ODE error | $E_{n+1} = RE_n - h\tau_n$ | $E_n = R^n E_0 - h\sum_{j=0}^{n-1} R^{n-1-j}\tau_j$ |

where $R = 1 + h\lambda$ is the amplification factor.

Both are **discrete convolutions with a geometric kernel**. The convergence/stability condition is the same: $|R| < 1$ (ODE stability) parallels $\|A\| < 1$ (Neumann series convergence).

This also discretizes the **continuous variation of constants formula**:
$$
u(t) = e^{\lambda t}u_0 + \int_0^t e^{\lambda(t-s)}g(s)\,ds
$$

The exponential $e^{\lambda t}$ becomes the power $R^n = (1+h\lambda)^n$, and the integral becomes a sum.
:::

### Convergence Theorem

:::{prf:theorem} Convergence of One-Step Methods
:label: thm-one-step-convergence

Let $u(t)$ be the solution to the initial value problem with $f$ sufficiently smooth. Let $\Psi$ be the step function of a one-step method that is Lipschitz continuous in $u$ with constant $L_\Psi$. If the method has consistency order $p$, i.e.,

$$
\|\tau(t + h)\| \leq Ch^{p} \quad \text{for } h \in (0, h_0], \; t \in [t_0, t_e - h]
$$

then the global discretization error satisfies:

$$
\|E_n\| \leq \frac{C}{L_\Psi}\left(e^{L_\Psi(t_n - t_0)} - 1\right) h_{\max}^p
$$
:::

:::{prf:proof}
:class: dropdown

We prove this for forward Euler on a nonlinear problem; the general case is similar.

**Step 1: Error recurrence.**
The numerical method gives $u_{n+1} = u_n + hf(t_n, u_n)$.
The exact solution satisfies $u(t_{n+1}) = u(t_n) + hf(t_n, u(t_n)) + h\tau_n$.

Subtracting:
$$
E_{n+1} = E_n + h[f(t_n, u_n) - f(t_n, u(t_n))] - h\tau_n
$$

**Step 2: Lipschitz bound.**
Since $f$ is Lipschitz with constant $L$:
$$
|f(t_n, u_n) - f(t_n, u(t_n))| \leq L|E_n|
$$

Thus:
$$
|E_{n+1}| \leq |E_n| + hL|E_n| + h|\tau_n| = (1 + hL)|E_n| + h|\tau_n|
$$

**Step 3: Unroll the recurrence.**
Let $\tau_{\max} = \max_j |\tau_j|$. Then:
$$
|E_n| \leq (1 + hL)^n|E_0| + h\tau_{\max}\sum_{j=0}^{n-1}(1 + hL)^j
$$

**Step 4: Bound the geometric sum.**
Using $(1 + hL)^k \leq e^{khL}$ and $\sum_{j=0}^{n-1}(1+hL)^j \leq \frac{(1+hL)^n - 1}{hL}$:
$$
|E_n| \leq e^{nhL}|E_0| + \frac{e^{nhL} - 1}{L}\tau_{\max}
$$

**Step 5: Substitute $nh = t_n - t_0$ and $\tau_{\max} = O(h^p)$.**
With $E_0 = 0$ (exact initial condition):
$$
|E_n| \leq \frac{e^{L(t_n - t_0)} - 1}{L} \cdot Ch^p = O(h^p)
$$
:::

:::{prf:corollary} Consistency Implies Convergence
:label: cor-consistency-convergence

For one-step methods with Lipschitz continuous step function:

$$
\text{Consistency of order } p \implies \text{Convergence of order } p
$$
:::

:::{prf:proof}
:class: dropdown

This follows directly from [](#thm-one-step-convergence). If the method has consistency order $p$, then $\|\tau\| \leq Ch^p$. The theorem gives:

$$
\|E_n\| \leq \frac{C}{L_\Psi}\left(e^{L_\Psi T} - 1\right) h^p = O(h^p)
$$

The key point: the constant $\frac{C}{L_\Psi}(e^{L_\Psi T} - 1)$ depends on the problem (through $L_\Psi$ and $T$) but not on $h$. Thus consistency of order $p$ implies convergence of order $p$.
:::

:::{prf:remark} The Forward/Backward Error Framework for ODEs
:label: rmk-fwd-bwd-error-ode
:class: dropdown

The convergence theorem is the ODE version of our fundamental error relationship:

$$
\text{Forward error} \leq \text{Condition number} \times \text{Backward error}
$$

In the ODE setting:
- **Backward error** = Local truncation error $\tau_n = O(h^p)$
- **Forward error** = Global error $\|E_n\| = O(h^p)$
- **Condition number** = $e^{L(t_n - t_0)}$ (Lipschitz amplification over time)

The Lipschitz constant $L$ measures how sensitive the ODE is to perturbations. The factor $e^{LT}$ is the **condition number for time evolution**—a property of the *problem*, not the algorithm.

A consistent method (small backward error per step) achieves convergence (controlled forward error) as long as the problem isn't too ill-conditioned (moderate $L$).
:::

## Stability Analysis

Convergence (as $h \to 0$) doesn't guarantee useful results for fixed $h > 0$. We need **absolute stability**.

The convergence theorem tells us that errors are bounded as $h \to 0$. But in practice, we use a *fixed* step size $h > 0$. The question becomes: does the error recurrence

$$
E_{n+1} = (1 + h\lambda)E_n - h\tau_n
$$

produce bounded errors for our chosen $h$? This depends on the **amplification factor** $R = 1 + h\lambda$. If $|R| > 1$, errors grow exponentially with each step—the method is unstable. Stability is the condition that ensures our backward errors (local truncation errors) don't get amplified catastrophically as we march forward in time.

### The Test Equation

Apply forward Euler to $u' = \lambda u$ with $\text{Re}(\lambda) < 0$ (decaying solutions):

$$
u_{n+1} = (1 + h\lambda)u_n = R(z)u_n, \quad z = h\lambda
$$

where $R(z) = 1 + z$ is the **stability function**.

### Stability Region

For the numerical solution to remain bounded, we need $|R(z)| = |1 + z| \leq 1$.

:::{prf:definition} Stability Region
:label: def-stability-region-fwd-euler

The **stability region** of forward Euler is:

$$
\mathcal{S} = \{z \in \mathbb{C} : |1 + z| \leq 1\}
$$

This is a disk of radius 1 centered at $(-1, 0)$.
:::

For real $\lambda < 0$: stability requires $-2 \leq h\lambda \leq 0$, so:

$$
h \leq \frac{2}{|\lambda|}
$$

Forward Euler is **conditionally stable**—the step size is restricted by the problem's eigenvalues.

:::{prf:remark} The Exponential Factor
:label: rmk-exponential-factor
:class: dropdown

The factor $e^{L(t_n - t_0)}$ in the convergence bound is unavoidable—it reflects the **problem's sensitivity** (condition number), not the algorithm. The Lipschitz constant $L$ measures how quickly nearby solution curves diverge or converge:

- If $L = 0$ (e.g., $u' = g(t)$): solution curves are parallel, errors don't grow
- If $L > 0$: solution curves diverge, errors can grow exponentially
- If $L < 0$: solution curves converge, errors are damped
:::

## Limitations

Forward Euler's conditional stability becomes problematic for **stiff problems**—where large eigenvalues force tiny step sizes even when the solution varies slowly.

For the test equation with $\lambda = -1000$:
- Stability requires $h < 0.002$
- But the solution $e^{-1000t}$ decays to negligible levels by $t = 0.01$

We're forced to take thousands of tiny steps tracking a transient we don't care about. This motivates **implicit methods** like backward Euler, which we develop next.
