# R1c: interfaces and matched observation

## Why revisit R1b?

At n=48, t=.27, the published R1b branch report found vacuum-relative errors of 0.000384411 (bare), 0.107570 (transformed), 0.097334 (density-only), and 1.8184e-12 (virtual). These are unfavorable to the transformed medium. The initial pulse support overlapped the outer annulus by .01 and observation was early; each material used different stable time steps. These factors confound interpretation. Source: https://github.com/MichaelWave369/PhiSpacetime/actions/runs/45917730266. No historical result is removed.

## Aligned transmission verification

A piecewise linear periodic coordinate X(x) has slopes 0.5 and 1.5 across the two half domains. The coefficients Axx=1/X', Ayy=rho=X' yield the exact standing mode sin(2 pi X(x)) cos(2 pi t). Both displacement and normal flux are continuous at the interior and periodic interfaces. The numerical field converges against this exact solution on aligned meshes (16,32,64); gate: observed order >1.5 and finest amplitude error <.01. This is an aligned-interface solution check, not an independently reconstructed numerical flux jump test. Curved geometry remains stair-stepped.

## Independent temporal verification

At fixed n=16, use the known discrete Fourier eigenmode of the vacuum mass-lumped operator. Verify the eigenvector residual and compare leapfrog against the exact semidiscrete cosine evolution, halving dt twice; order in [1.9,2.1]. This isolates time integration from spatial discretization. It does not validate temporal error for every transformed material.

## Matched exterior observations

Pulse support lies strictly outside the annulus, identical physical/virtual initial conditions there. All five configurations use the same time steps per run, selected from the most restrictive stability limit; same observation nodes outside R=.25. Report vacuum and virtual relative L2 errors, energy drift, mask counts, and time steps. The finite regularization makes the virtual obstacle the relevant reference rather than perfect vacuum. The main sweep uses n=32,48,64 and t=.38, with additional half-dt, later t=.46 and off-axis illumination cases. No ranking is a pass criterion; results are exploratory.

Periodic reentry, coarse curved masks, quadrature across discontinuities, and finite-time observation norm can bias comparisons. A rigorous fitted-interface or cut-cell discretization and absorbing/open domain are still required for a physical scattering interpretation. The next gate is a curved-interface convergence study with geometry-aware assembly and boundary observation histories. Adaptive pocket integration should not proceed on the assumption that the current cloak candidate succeeds.

Run python -m phispacetime.interface_benchmark for the complete JSON report. Record all cases, including failures.
