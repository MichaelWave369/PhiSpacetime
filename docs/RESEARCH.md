# Research charter and acceptance gates

## Hypothesis H0
A nine-way spatial pocket hierarchy may improve accuracy per resource over uniform discretization or conventional adaptive refinement for some geometric field problems. No physical energy reduction or numerical superiority is assumed.

## R0 scope
Implement a regularized annular coordinate map, principal stretches, transformed relative constitutive tensor in a local orthonormal basis (vacuum starting medium), and an analytic scalar quadrature benchmark. The one-dimensional triadic subdivision is a precursor to a 3x3 spatial pocket tree, not an integration of existing code. The coordinate map is prescribed, not solved from Maxwell or Einstein equations. Exterior ray agreement follows by construction and is not experimental evidence of invisibility.

## R1: scalar wave baseline
Implement a verified 2D scalar Helmholtz or time-domain wave solver. Use manufactured solutions, documented boundaries, convergence studies, vacuum and uncloaked controls. State polarization, dispersion, material loss, discretization and absorbing boundary assumptions. Record complex exterior field residual and energy balance; a scalar result cannot establish full vector cloaking.

## R2: comparative refinement
Compare fixed mesh, established binary/quadtree adaptive mesh and 3x3 pocket refinement using identical equations, tolerances, initial/boundary conditions and hardware. Report error, degrees of freedom, field evaluations, runtime and peak memory. Pre-register tolerances, use multiple wavelengths and withheld geometries, publish failures and resource limits. Only then integrate audited NBG Pocket Mesh and PhiMycelium adapters.

## R3–R8 roadmap
| Track | Benchmark and guardrail |
|---|---|
| Metric | Minkowski and Schwarzschild curvature; coordinate invariance |
| Wormhole | Morris–Thorne throat in proper radial coordinate; null energy condition |
| Warp | Alcubierre metric and independent stress-energy calculation; no propulsion claim |
| Vacuum | parallel-plate Casimir idealization; no macroscopic negative energy claim |
| NullMesh | documented interface matching, flux conservation and causal structure |
| Altermath | boundary observation operator, reconstruction error and non-identifiability |
| Validator | independent convergence, dimensional consistency, reproducibility and provenance |

## Integration policy
Nested Bubble Gear, NBG Pocket Mesh, Infinite Holographic Altermath, PhiMycelium and PhiNullMesh are candidate computational abstractions; their APIs and physical relevance must be audited before use. Record source repository, exact commit, license, adapter tests and a conventional baseline. No special significance is assigned to 3/6/9 or golden-ratio parameters absent controlled evidence.

## Evidence ladder
Analytic identity -> verified numerical reproduction -> benchmark comparison -> independent replication -> physical experiment. Never promote a result beyond the level supported. Each experiment records assumptions, units, parameters, code revision, environment and uncertainty. A report's government provenance is not scientific validation.
