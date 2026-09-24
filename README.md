# Cahn–Hilliard Chemomechanical Particle Model

A MATLAB finite-volume model of phase-separating transport coupled to mechanical deformation in a spherical electrode particle, with optional elastic or elastoplastic surface-shell mechanics.

The model combines a radially symmetric Cahn–Hilliard transport equation with quasistatic mechanical equilibrium. Concentration generates chemical strain, while the hydrostatic stress contribution enters the chemical potential. An optional reduced-order shell modifies the radial traction at the particle surface.

> **Project status:** The numerical solver, post-processing pipeline, and core visualisation routines are implemented. Analytical verification, convergence studies, automated tests, and publication-quality example figures are still being developed.

## Features

- Spherically symmetric Cahn–Hilliard transport
- Concentration-dependent chemical strain
- Stress-dependent chemical potential
- Quasistatic particle mechanics
- Conservative spherical finite-volume discretisation
- Implicit differential-algebraic equation solution using `ode15i`
- Consistent initial conditions using `decic`
- Traction-free, elastic-shell, and elastoplastic-shell modes
- Event-driven detection of initial shell yielding
- Separate numerical solution, post-processing, and visualisation stages

## Shell modes

The surface model is selected using `params.shellModel`.

| Mode | Setting | Surface behaviour |
|---|---|---|
| Traction-free | `"traction-free"` | Zero radial surface traction |
| Elastic shell | `"elastic"` | Shell traction remains elastic throughout the simulation |
| Elastoplastic shell | `"elastoplastic"` | The shell begins elastically and switches to a yield-limited plastic response |

Even in traction-free mode, the particle retains bulk chemomechanical coupling: concentration affects deformation and stress, and stress affects the chemical potential.

## Model formulation

The discrete state vector is

```text
w = [u_1, ..., u_N, c_1, ..., c_N, ep]
```

where:

- `u` is radial displacement;
- `c` is normalised concentration;
- `ep` is the reduced shell plastic-strain variable.

The chemical potential is implemented in the form

```text
mu = mu_chem(c) - kappa*laplacian(c) - gamma*trace(sigma)/3
```

The homogeneous chemical contribution is currently

```text
mu_chem(c) = log(c/(1-c)) + 3*(1-2*c)
```

The transport equation is discretised in conservative spherical finite-volume form. Mechanical equilibrium and the shell constraint are algebraic equations, producing a coupled DAE system.

## Numerical workflow

The code is split into three stages:

```text
parameters
    ↓
numerical_solver
    ↓
postProcessing
    ↓
visualisation
```

1. `numerical_solver` integrates the coupled DAE and returns the raw state history.
2. `postProcessing` reconstructs concentration, displacement, stress, chemical potential, flux-related quantities, and average concentration.
3. `visualisation` plots simulation summaries, radial profiles, and shell-specific quantities.

This separation allows saved simulations to be analysed and replotted without rerunning the numerical solver.

## Repository structure

```text
.
├── master.m              # Main example and model configuration
├── parameters.m          # Geometry, material, shell, and loading parameters
├── numerical_solver.m    # Loading segments and raw state assembly
├── solve_segment.m       # DAE integration and shell-regime switching
├── RHS_particle.m        # Coupled transport-mechanics residual
├── elastEvent.m          # Elastic-shell yield event
├── mu_chem.m             # Homogeneous chemical potential
├── postProcessing.m      # Derived fields and scalar quantities
├── visualisation.m       # Summary, profile, and shell plots
├── validation/           # Analytical comparisons (planned)
├── tests/                # Regression and conservation tests (planned)
└── figures/              # Selected exported figures (planned)
```

## Requirements

- MATLAB
- A MATLAB release containing `ode15i`, `decic`, and `tiledlayout`

No third-party MATLAB packages are currently required.

## Running a simulation

Clone the repository, open its root directory in MATLAB, and select a shell model in `master.m`:

```matlab
params = parameters();

params.shellModel = "traction-free";
% params.shellModel = "elastic";
% params.shellModel = "elastoplastic";
```

Run the full workflow:

```matlab
raw = numerical_solver(params);
sol = postProcessing(raw, params);
figs = visualisation(sol, params);
```

Alternatively, run:

```matlab
master
```

The current `master.m` also creates a reduced `sol_save` structure for storing selected outputs at lower radial resolution.

## Main parameters

The principal controls are collected in `parameters.m`.

| Group | Parameters |
|---|---|
| Mesh | `R`, `nbin`, `dr`, `r_bin`, `r_edge` |
| Particle mechanics | `E_si`, `nu_si`, `alph`, `nu` |
| Phase-field transport | `D0`, `kappa`, `gamm`, `c_ref` |
| Shell model | `shellModel`, `A`, `delta`, `sig_y` |
| Physical constants | `Far`, `R_gT`, `V_c`, `c_m` |
| Loading protocol | `durations`, `currents` |

The current default mesh uses 400 radial control volumes. Applied loading is supplied as matching `durations` and `currents` arrays, allowing multiple piecewise-constant current segments.

## Raw solver output

`numerical_solver` returns a structure named `raw` containing:

| Field | Description |
|---|---|
| `raw.sol_total` | Full DAE state history |
| `raw.t_total` | Simulation time vector |
| `raw.I` | Applied current/flux history |
| `raw.t_plast` | Detected yield-transition times |
| `raw.t_elast` | Recorded return-to-elastic times |

## Post-processed output

`postProcessing` returns a structure named `sol` containing:

| Field | Description |
|---|---|
| `sol.t_total` | Simulation time vector |
| `sol.SOC` | Radially averaged concentration |
| `sol.c` | Concentration field |
| `sol.u` | Radial displacement field |
| `sol.ep` | Shell plastic-strain history |
| `sol.sigma_rr` | Radial stress field |
| `sol.sigma_tt` | Hoop stress field |
| `sol.sigma_kk` | Stress trace |
| `sol.mu` | Chemical-potential field |
| `sol.I` | Applied current/flux history |

Additional strain, elastic-work, and flux outputs are calculated internally and can be enabled in `postProcessing.m` if required.

## Visualisation

`visualisation.m` automatically adapts to the selected shell mode. It produces:

- average concentration against time;
- surface concentration and displacement histories;
- applied current/flux history;
- surface chemical-potential history;
- radial concentration and displacement profiles;
- radial chemical-potential and stress profiles;
- shell traction, trial stress, and plastic strain when a shell is active.

The function returns figure handles in the `figs` structure:

```matlab
figs = visualisation(sol, params);
```

Plot limits, labels, and formatting can be adjusted in the relevant plotting block inside `visualisation.m`.

## Planned validation

The following verification and comparison cases are planned:

- lithium conservation against the imposed surface flux;
- mechanical-equilibrium residual checks;
- traction-free surface-stress verification;
- spatially homogeneous solutions;
- analytical displacement and stress solutions;
- analytical homogeneous and phase-separated chemical potentials;
- mesh-convergence studies;
- solver-tolerance studies;
- regression tests for each shell mode;
- direct comparison of traction-free, elastic, and elastoplastic responses.

These files will be placed in `validation/` and `tests/` rather than mixed with the core solver.

## Current limitations

- The model assumes spherical symmetry.
- Particle mechanics is quasistatic.
- The shell is represented by a reduced surface constitutive law rather than a spatially resolved shell domain.
- The present implementation uses nondimensional model equations, while some dimensional material parameters are retained for scaling and future development.
- Voltage and overpotential are not included in the current post-processed output.
- The current elastoplastic implementation detects initial yielding, but the unloading and reverse-yield treatment is still under development.
- Analytical validation and automated regression tests have not yet been added to the public workflow.
- The code is research software and should not yet be used as a validated predictive engineering tool.

## Reproducibility

The current initial concentration is spatially uniform. If a random perturbation is reintroduced to initiate phase separation, set and record the MATLAB random-number seed:

```matlab
rng(1);
```

For reproducible results, also record the MATLAB release, complete parameter structure, mesh resolution, loading protocol, and solver tolerances.

## Development roadmap

- Add analytical-comparison functions
- Add mass-conservation and mechanical-residual diagnostics
- Add automated tests for all three shell modes
- Implement and verify unloading/reverse-yield events
- Add voltage and overpotential reconstruction
- Add convergence studies
- Export a small set of representative figures
- Add a software citation file and archived release

## Citation

A `CITATION.cff` file will be added when the associated publication or thesis citation is finalised.

## Licence

A licence should be added before public release. The choice must be consistent with any university, collaborator, publication, or funder requirements.

## Author

**Joseph Petrassem de Sousa**

Research interests include phase-transforming battery electrodes, continuum thermodynamics, chemomechanical coupling, and numerical modelling of lithium-ion battery materials.
