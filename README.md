# A covariance-matched Langevin integrator preserving harmonic equilibrium


## Files

| File | Purpose |
|---|---|
| `integrator_COVLD_setup.m` | Checks the inputs and time-step admissibility, constructs the covariance-matched noise factor, and factorizes the mechanical operators once. |
| `integrator_COVLD.m` | Advances displacement and stored endpoint velocity by one predictor-corrector step. |
| `example_harmonic_oscillator.m` | Runs a harmonic oscillator ensemble and compares trajectory-based, discrete, continuous, and canonical statistics. |

## Quick start

Place the three `.m` files in the same folder, make that folder the MATLAB current folder, and run

```matlab
example_harmonic_oscillator
```

For direct use,

```matlab
op = integrator_COVLD_setup(M, K, Z, kBT, dt);
[U, V] = integrator_COVLD(op, U, V);
```

## Model and inputs

The physical displacement and velocity are collected in `U` and `V`.
The constant physical mass, stiffness, and friction matrices are `M`, `K`, and `Z`.
The model is

```math
\mathbf M\ddot{\mathbf U}
+
\mathbf Z\dot{\mathbf U}
+
\mathbf K\mathbf U
=
\mathbf R(t)
+
\boldsymbol{\beta}(t),
\qquad
\left\langle
\boldsymbol{\beta}(t)\boldsymbol{\beta}(t')^{\mathrm T}
\right\rangle
=
2k_BT\,\mathbf Z\,\delta(t-t').
```

Use consistent physical units.
`kBT` is the thermal energy.
Do not mass-normalize `K` or `Z` before supplying them with the physical `M`.

`K` defines the number of degrees of freedom `N`.
It may be a scalar, a vector of diagonal entries, or an `N`-by-`N` matrix.
`M` and `Z` may additionally be scalars representing multiples of the identity.
The current setup requires symmetric positive-definite `M`, `K`, and `Z`, positive `kBT`, and a numerically resolved strictly admissible `dt`.

For `L` independent trajectories, `U` and `V` are `N`-by-`L` arrays with one trajectory per column.

```matlab
M = [2.0, 0.2;
     0.2, 1.0];
K = [3.0, -1.0;
     -1.0, 2.0];
Z = [1.0, 0.1;
     0.1, 0.8];

op = integrator_COVLD_setup(M, K, Z, 1.0, 0.05);
U = zeros(2,1);
V = zeros(2,1);

[U, V] = integrator_COVLD(op, U, V);
```

Optional fourth and fifth arguments supply the physical external forces at the current and next time nodes.
An optional sixth argument supplies the independent standard normal variables.

```matlab
R_n = zeros(op.N,1);
R_next = zeros(op.N,1);
eta = randn(2*op.N,size(U,2));

[U, V] = integrator_COVLD(op,U,V,R_n,R_next,eta);
```

When `eta` is omitted, the function generates it with MATLAB `randn`.
The displacement and velocity noise blocks use the same `eta` so that their required cross covariance is retained.

## Harmonic oscillator example

The default example uses `M = K = Z = kBT = 1`, `dt = 0.1`, and `t_end = 10` in reduced units.
It advances 20,000 independent trajectories from `U = V = 0` and stores ten complete sample paths.
The initial ensemble is therefore deliberately out of equilibrium.

### Figs. 1 and 2

`fig_harmonic_displacement` and `fig_harmonic_velocity` show ten independent stochastic realizations of the same COVLD integrator.
The ten curves are different trajectories, not different numerical methods.

For trajectory `r`,

```math
\mathbf x_{n+1}^{(r)}
=
\mathbf A\mathbf x_n^{(r)}
+
\mathbf G\boldsymbol{\eta}_{n+1}^{(r)},
\qquad
\mathbf x_n^{(r)}
=
\begin{bmatrix}
u_n^{(r)}\\
v_n^{(r)}
\end{bmatrix}.
```

The Gaussian vectors are independent between trajectories and time steps.
Figs. 1 and 2 therefore show individual thermal trajectories in displacement and velocity.

### Figs. 3 and 4

`fig_harmonic_displacement_variance` and `fig_harmonic_velocity_variance` compare four covariance quantities.

**Canonical value — black solid line.**  
This is the equilibrium target.

```math
\boldsymbol{\Sigma}_{\mathrm{eq}}
=
\begin{bmatrix}
k_BT/K & 0\\
0 & k_BT/M
\end{bmatrix}.
```

The displacement and velocity variances are therefore

```math
\operatorname{Var}_{\mathrm{eq}}(u)=\frac{k_BT}{K},
\qquad
\operatorname{Var}_{\mathrm{eq}}(v)=\frac{k_BT}{M}.
```

**Ensemble — blue solid line.**  
This is the finite-ensemble covariance measured directly from the 20,000 simulated trajectories.
With `L` trajectories,

```math
\overline{\mathbf x}_n
=
\frac{1}{L}
\sum_{r=1}^{L}\mathbf x_n^{(r)},
```

```math
\widehat{\boldsymbol{\Sigma}}_n
=
\frac{1}{L-1}
\sum_{r=1}^{L}
\left(
\mathbf x_n^{(r)}-\overline{\mathbf x}_n
\right)
\left(
\mathbf x_n^{(r)}-\overline{\mathbf x}_n
\right)^{\mathrm T}.
```

Thus, the ensemble curve is obtained from the stochastic trajectories themselves.

**Discrete covariance recursion — red dashed line.**  
This is the covariance of the discrete COVLD update computed directly, without trajectory sampling.

```math
\boldsymbol{\Sigma}_{n+1}^{\mathrm d}
=
\mathbf A
\boldsymbol{\Sigma}_n^{\mathrm d}
\mathbf A^{\mathrm T}
+
\mathbf Q,
```

where the covariance-matched increment satisfies

```math
\mathbf Q
=
\boldsymbol{\Sigma}_{\mathrm{eq}}
-
\mathbf A
\boldsymbol{\Sigma}_{\mathrm{eq}}
\mathbf A^{\mathrm T}.
```

For this example, `U = V = 0` initially, so

```math
\boldsymbol{\Sigma}_0^{\mathrm d}=\mathbf 0.
```

**Continuous reference — green dotted line.**  
This is the exact covariance of the continuous linear Langevin system.
For the scalar oscillator,

```math
\mathbf A_{\mathrm c}
=
\begin{bmatrix}
0 & 1\\
-K/M & -Z/M
\end{bmatrix},
\qquad
\mathbf E
=
\exp\!\left(
\Delta t\,\mathbf A_{\mathrm c}
\right).
```

The exact one-step covariance increment is

```math
\mathbf Q_{\mathrm c}
=
\boldsymbol{\Sigma}_{\mathrm{eq}}
-
\mathbf E
\boldsymbol{\Sigma}_{\mathrm{eq}}
\mathbf E^{\mathrm T},
```

and the continuous reference is propagated as

```math
\boldsymbol{\Sigma}_{n+1}^{\mathrm c}
=
\mathbf E
\boldsymbol{\Sigma}_n^{\mathrm c}
\mathbf E^{\mathrm T}
+
\mathbf Q_{\mathrm c}.
```

The variance plots are normalized by the corresponding canonical values,

```math
\frac{\mathrm{Var}(u)}{k_BT/K},
\qquad
\frac{\operatorname{Var}(v)}{k_BT/M},
```

so the canonical line is equal to one.

The difference between **Ensemble** and **Discrete covariance recursion** is finite-ensemble sampling error.
The difference between **Discrete covariance recursion** and **Continuous reference** is time-discretization error.
At long times, the discrete covariance approaches the same canonical target by construction at an admissible time step.

Each run saves MAT/CSV data, source snapshots, a log, and FIG/EPS/PNG figures under

```text
example_harmonic_oscillator/YYYYMMDD_HHMMSS/
```

## Requirements

MATLAB R2016b or later.
No additional toolbox is required.

## Citation

> J. Y. Lee, *A covariance-matched Langevin integrator preserving harmonic equilibrium*, submitted, 2026.
