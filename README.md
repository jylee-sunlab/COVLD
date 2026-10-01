# A covariance-matched Langevin integrator preserving harmonic equilibrium

## Files

| File | Purpose |
|---|---|
<<<<<<< HEAD
| `integrator_PGLGV_setup.m` | Checks the inputs and time-step admissibility, constructs the noise factor, and factorizes the mechanical operators once. |
| `integrator_PGLGV.m` | Advances displacement and stored endpoint velocity by one predictor–corrector step. |
=======
| `integrator_COVLD_setup.m` | Checks the inputs and time-step admissibility, constructs the noise factor, and factorizes the mechanical operators once. |
| `integrator_COVLD.m` | Advances displacement and stored endpoint velocity by one predictor–corrector step. |
>>>>>>> 731e7f0 (update)
| `example_harmonic_oscillator.m` | Runs a harmonic oscillator ensemble, compares its moments with discrete and continuous reference covariances, and saves the results. |

## Quick start

<<<<<<< HEAD
Place the four `.m` files in the same folder, make that folder the MATLAB current folder, and run
=======
Place the three `.m` files in the same folder, make that folder the MATLAB current folder, and run
>>>>>>> 731e7f0 (update)

```matlab
example_harmonic_oscillator
```

<<<<<<< HEAD

=======
>>>>>>> 731e7f0 (update)
## Model and inputs

The physical displacement and velocity vectors are denoted by $\boldsymbol U$ and $\boldsymbol V$.
The constant mass, stiffness and friction matrices are $\boldsymbol M$, $\boldsymbol K$ and $\boldsymbol Z$.
The prescribed deterministic external force is $\boldsymbol R(t)$, and the zero-mean Gaussian white thermal force is $\boldsymbol\beta(t)$.
The model is

$$
\boldsymbol M\ddot{\boldsymbol U}+\boldsymbol Z\dot{\boldsymbol U}+\boldsymbol K\boldsymbol U=\boldsymbol R(t)+\boldsymbol\beta(t),
\qquad
\langle\boldsymbol\beta(t)\boldsymbol\beta(t')^{\mathrm T}\rangle=2k_BT\boldsymbol Z\delta(t-t').
$$

Use consistent physical units throughout.
`kBT` is the thermal energy, not the temperature in kelvin.
For a scalar oscillator, the mass-normalized friction is `Z/M` and the undamped frequency is `sqrt(K/M)`.
Do not mass-normalize `K` or `Z` before supplying them with the physical `M`.

`K` defines the number of degrees of freedom `N`.
It can be a scalar, a vector of diagonal entries, or an `N`-by-`N` matrix.
`M` and `Z` can additionally be scalars representing multiples of the identity.
All inputs must be finite, real double-precision values.
<<<<<<< HEAD
=======
`M`, `K` and `Z` must be symmetric positive definite, with positive `kBT` and `dt`.
Only numerically resolved strictly admissible steps are accepted, without covariance clipping or automatic time-step adjustment.
>>>>>>> 731e7f0 (update)

For `L` independent trajectories, use `N`-by-`L` state arrays.
A multi-degree-of-freedom state is a column, not a row.
For example, the same functions can be used without modal decoupling as follows.

```matlab
<<<<<<< HEAD
M = [2.0, 0.2; 0.2, 1.0];
K = [3.0, -1.0; -1.0, 2.0];
Z = [1.0, 0.1; 0.1, 0.8];
op = integrator_PGLGV_setup(M, K, Z, 1.0, 0.05);
U = zeros(2, 1);
V = zeros(2, 1);
[U, V] = integrator_PGLGV(op, U, V);
=======
M = [2.0, 0.2;
     0.2, 1.0];
K = [3.0, -1.0;
     -1.0, 2.0];
Z = [1.0, 0.1;
     0.1, 0.8];
op = integrator_COVLD_setup(M, K, Z, 1.0, 0.05);
U = zeros(2, 1);
V = zeros(2, 1);
[U, V] = integrator_COVLD(op, U, V);
>>>>>>> 731e7f0 (update)
```

### Prescribed forces and random variables

Optional fourth and fifth arguments supply the physical external forces at the current and next time nodes.
Each force can be `N`-by-1, shared by all trajectories, or `N`-by-`L`.

```matlab
R_n = zeros(op.N, 1);
R_next = zeros(op.N, 1);
<<<<<<< HEAD
[U, V] = integrator_PGLGV(op, U, V, R_n, R_next);
=======
[U, V] = integrator_COVLD(op, U, V, R_n, R_next);
>>>>>>> 731e7f0 (update)
```

An omitted or empty `R_n` is zero.
An omitted or empty `R_next` equals `R_n`.
The internal force is evaluated as `op.K*U` inside the step function.

```matlab
eta = randn(2*op.N, size(U, 2));
<<<<<<< HEAD
[U, V] = integrator_PGLGV(op, U, V, R_n, R_next, eta);
=======
[U, V] = integrator_COVLD(op, U, V, R_n, R_next, eta);
>>>>>>> 731e7f0 (update)
```

When `eta` is omitted or empty, the function uses MATLAB `randn` and the current random-number-generator state.
A supplied `eta` does not change the random-number-generator state.
For the stated statistical properties, the draws must be independent across steps and trajectories and independent of the initial state.
The two displacement and velocity noise blocks use the **same** `eta` to retain the required cross covariance.

<<<<<<< HEAD

=======
>>>>>>> 731e7f0 (update)
## Harmonic oscillator example

The default example uses `M = K = Z = kBT = 1`, `dt = 0.1`, and `t_end = 20` in reduced units.
It advances 20,000 independent trajectories from zero displacement and velocity and stores four complete sample paths.
The initial ensemble is deliberately out of equilibrium.

<<<<<<< HEAD
The example calculates the sample mean and centered covariance, the covariance of the discrete update, and an exact continuous reference obtained with `expm`.
The discrete covariance recursion does not use Monte Carlo sampling.
A small covariance-balance residual checks equilibrium preservation, while the difference between discrete and continuous transient covariances measures time-discretization error.
Finite ensemble estimates retain sampling error even when the discrete equilibrium covariance is preserved.


## Requirements and validation

The scripts are written for MATLAB R2016b or later and use base MATLAB functions only.


## Citation

> J. Y. Lee, *A covariance-matched Langevin integrator preserving harmonic equilibrium*, submitted, 2026.

=======
### Figures 1 and 2

`fig_harmonic_displacement` and `fig_harmonic_velocity` show four independent stochastic realizations of the same integrator.
The four curves are not four different methods.
For realization $r$, the discrete harmonic state is

$$
\boldsymbol x_{n+1}^{(r)}
=
\boldsymbol A\boldsymbol x_n^{(r)}+\boldsymbol G\boldsymbol\eta_{n+1}^{(r)},
\qquad
\boldsymbol x_n^{(r)}=
\begin{bmatrix}u_n^{(r)}\\v_n^{(r)}\end{bmatrix},
$$

where the Gaussian vectors $\boldsymbol\eta_{n+1}^{(r)}$ are independent between trajectories and time steps.
The figures therefore illustrate thermal fluctuations of individual paths.

### Figures 3 and 4

`fig_harmonic_displacement_variance` and `fig_harmonic_velocity_variance` compare four different covariance quantities.
Let $L=20000$ denote the number of trajectories.

1. **Ensemble** is the finite-ensemble estimate obtained directly from the simulated trajectories,

$$
\widehat{\boldsymbol\Sigma}_n
=
\frac{1}{L-1}
\sum_{r=1}^{L}
\left(\boldsymbol x_n^{(r)}-\overline{\boldsymbol x}_n\right)
\left(\boldsymbol x_n^{(r)}-\overline{\boldsymbol x}_n\right)^{\mathrm T}.
$$

2. **Discrete covariance recursion** is the exact covariance of the discrete integrator and does not require trajectory sampling,

$$
\boldsymbol\Sigma_{n+1}^{\mathrm d}
=
\boldsymbol A\boldsymbol\Sigma_n^{\mathrm d}\boldsymbol A^{\mathrm T}+\boldsymbol Q,
\qquad
\boldsymbol Q
=
\boldsymbol\Sigma_{\mathrm{eq}}
-
\boldsymbol A\boldsymbol\Sigma_{\mathrm{eq}}\boldsymbol A^{\mathrm T}.
$$

3. **Continuous reference** is the exact covariance of the continuous linear Langevin system.
With

$$
\boldsymbol A_{\mathrm c}
=
\begin{bmatrix}
0 & 1\\
-K/M & -Z/M
\end{bmatrix},
\qquad
\boldsymbol E=e^{\Delta t\boldsymbol A_{\mathrm c}},
$$

its one-step covariance recursion is

$$
\boldsymbol\Sigma_{n+1}^{\mathrm c}
=
\boldsymbol E\boldsymbol\Sigma_n^{\mathrm c}\boldsymbol E^{\mathrm T}
+
\boldsymbol Q_{\mathrm c},
\qquad
\boldsymbol Q_{\mathrm c}
=
\boldsymbol\Sigma_{\mathrm{eq}}
-
\boldsymbol E\boldsymbol\Sigma_{\mathrm{eq}}\boldsymbol E^{\mathrm T}.
$$

4. **Canonical value** is the equilibrium target,

$$
\boldsymbol\Sigma_{\mathrm{eq}}
=
\begin{bmatrix}
k_BT/K & 0\\
0 & k_BT/M
\end{bmatrix}.
$$

The variance plots are normalized by the corresponding canonical values, so the canonical line is equal to one.
The difference between **Ensemble** and **Discrete covariance recursion** is finite-sample error, while the difference between **Discrete covariance recursion** and **Continuous reference** is time-discretization error.
Both discrete and continuous covariances relax toward the same canonical target.

Each run saves MAT/CSV data, source snapshots, a log, and FIG/EPS/PNG figures under `example_harmonic_oscillator/YYYYMMDD_HHMMSS/`.

## Requirements

The scripts are written for MATLAB R2016b or later and use base MATLAB functions only.

## Citation

> J. Y. Lee, *A covariance-matched Langevin integrator preserving harmonic equilibrium*, submitted, 2026.
>>>>>>> 731e7f0 (update)
