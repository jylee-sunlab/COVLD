# A covariance-matched Langevin integrator preserving harmonic equilibrium

## Files

| File | Purpose |
|---|---|
| `integrator_PGLGV_setup.m` | Checks the inputs and time-step admissibility, constructs the noise factor, and factorizes the mechanical operators once. |
| `integrator_PGLGV.m` | Advances displacement and stored endpoint velocity by one predictor–corrector step. |
| `example_harmonic_oscillator.m` | Runs a harmonic oscillator ensemble, compares its moments with discrete and continuous reference covariances, and saves the results. |

## Quick start

Place the four `.m` files in the same folder, make that folder the MATLAB current folder, and run

```matlab
example_harmonic_oscillator
```


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

For `L` independent trajectories, use `N`-by-`L` state arrays.
A multi-degree-of-freedom state is a column, not a row.
For example, the same functions can be used without modal decoupling as follows.

```matlab
M = [2.0, 0.2; 0.2, 1.0];
K = [3.0, -1.0; -1.0, 2.0];
Z = [1.0, 0.1; 0.1, 0.8];
op = integrator_PGLGV_setup(M, K, Z, 1.0, 0.05);
U = zeros(2, 1);
V = zeros(2, 1);
[U, V] = integrator_PGLGV(op, U, V);
```

### Prescribed forces and random variables

Optional fourth and fifth arguments supply the physical external forces at the current and next time nodes.
Each force can be `N`-by-1, shared by all trajectories, or `N`-by-`L`.

```matlab
R_n = zeros(op.N, 1);
R_next = zeros(op.N, 1);
[U, V] = integrator_PGLGV(op, U, V, R_n, R_next);
```

An omitted or empty `R_n` is zero.
An omitted or empty `R_next` equals `R_n`.
The internal force is evaluated as `op.K*U` inside the step function.

```matlab
eta = randn(2*op.N, size(U, 2));
[U, V] = integrator_PGLGV(op, U, V, R_n, R_next, eta);
```

When `eta` is omitted or empty, the function uses MATLAB `randn` and the current random-number-generator state.
A supplied `eta` does not change the random-number-generator state.
For the stated statistical properties, the draws must be independent across steps and trajectories and independent of the initial state.
The two displacement and velocity noise blocks use the **same** `eta` to retain the required cross covariance.


## Harmonic oscillator example

The default example uses `M = K = Z = kBT = 1`, `dt = 0.1`, and `t_end = 20` in reduced units.
It advances 20,000 independent trajectories from zero displacement and velocity and stores four complete sample paths.
The initial ensemble is deliberately out of equilibrium.

The example calculates the sample mean and centered covariance, the covariance of the discrete update, and an exact continuous reference obtained with `expm`.
The discrete covariance recursion does not use Monte Carlo sampling.
A small covariance-balance residual checks equilibrium preservation, while the difference between discrete and continuous transient covariances measures time-discretization error.
Finite ensemble estimates retain sampling error even when the discrete equilibrium covariance is preserved.

Every run creates a distinct directory of the form

```text
example_harmonic_oscillator/YYYYMMDD_HHMMSS/
```

The directory contains `inputs.mat`, `results.mat`, `harmonic_moments.csv`, `run.log`, and copies of the executed script and integrator functions.
The saved variables include the settings, operator, selected paths, final ensemble, moment histories, covariance errors and random-number-generator states.
Four figures show displacement paths, stored velocity paths, displacement variance and velocity variance.
Each figure is saved as `.fig`, `.eps` and `.png` with the same base filename.
The variables remain available in the MATLAB workspace.
The example is a small demonstration, not a script reproducing every figure in the manuscript.

The regression script similarly writes its results under `test_integrator_PGLGV/YYYYMMDD_HHMMSS/`.
Generated run directories are output data and need not be uploaded as repository source files.

## Requirements and validation

The scripts are written for MATLAB R2016b or later and use base MATLAB functions only.


## Citation

> J. Y. Lee, *A covariance-matched Langevin integrator preserving harmonic equilibrium*, submitted, 2026.

