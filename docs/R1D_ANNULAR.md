# R1d: fitted circular interfaces

R1c showed mixed results, with transformed vacuum-relative error falling from .214768 at n=32 to .100950 at n=64, versus .137940 for the bare disk at n=64. Layered and time convergence orders were about two. See https://github.com/MichaelWave369/PhiSpacetime/actions/runs/45919392073. This is not proof of cloaking.

## Fitted annular formulation

Use polar Q1 elements with exact circular coordinate boundaries, periodic theta, homogeneous Dirichlet wave boundaries at the inner and outer radii, and matched radial break at R=.22. Physical inner radius .12, virtual .04, outer .45. The virtual radial nodes are mapped into the physical domain by the R0 map; outside the transformation they coincide. Each mesh has n radial nodes and n angular nodes. The weak form uses radial stiffness A_r R, angular stiffness A_theta/R and mass rho R under dR dtheta. Two-point Gauss integration per coordinate, lumped mass, existing conservative leapfrog. No stair-stepped masks within the mesh.

Paired virtual/transformed matrices must agree up to floating point because their coordinate maps are exactly related. This is an implementation covariance check, not an independent numerical or physical discovery; it cannot establish reduced scattering for nonmapped meshes. The bare physical annulus uses the same physical nodes with unit coefficients.

## Independent static reference

Solve the Laplace Dirichlet problem with zero inner value and cos(theta) at the outer boundary. The virtual exact solution is (r-epsilon^2/r)/(L-epsilon^2/L)*cos(theta). The bare solution uses its own inner radius. The transformed exact solution is the virtual one evaluated at the inverse map. Conjugate gradient residual and resolution convergence are checked. This static test is independent of the time integrator but does not verify all dynamic interface fluxes.

## Dynamic comparison

Compact radial pulse times cos(theta), initially outside the transformation, common time step from the most restrictive stability limit, homogeneous fixed inner and outer boundaries. Observe common exterior nodes outside R=.25. Report full matrix/mass discrepancy, static and dynamic exterior errors, static residual and wave energy drift. The finite outer boundary reflects waves; this is a cavity, not open scattering. The virtual reference is not vacuum.

Run python -m phispacetime.annular_benchmark. Sweep n=17,33,65. Preserve R1c Cartesian results as an independent approximate comparison.

## Next gate

Compare independently generated physical and virtual meshes using common observation interpolation, test nonradial geometry and interface flux, add open-boundary verification, then compare uniform/quadtree/triadic adaptive discretizations. Paired coordinate invariance alone is not evidence that our pocket architecture improves numerical performance. No upstream integration or physical experiment is claimed.
