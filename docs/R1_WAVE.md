# R1a: scalar propagation control

## Equation and discretization

On the nondimensional periodic unit square, solve u_tt = u_xx + u_yy with wave speed 1. Spatial derivatives use the five-point second-order stencil; time uses centered leapfrog with Taylor startup. The Courant number dt/h is limited to 1/sqrt(2). A fixed-zero disk is imposed by zeroing masked degrees of freedom, with initial displacement and velocity compatible with the constraint. Grid-aligned disk approximation is not a convergence-validated boundary method.

The conserved discrete energy at half steps is the squared time difference plus the cross-time spatial gradient inner product (with the usual factor 1/2). This is not instantaneous continuum energy. Conservation is checked for unforced vacuum and obstacle controls.

## Preregistered checks

- Standing mode sin(2 pi x) cos(2 pi y) cos(2 pi sqrt(2) t), zero initial velocity; final t=0.25.
- Resolutions 24, 48, 96; Courant number 0.5; relative discrete L2 error and observed order log2(E_h/E_h/2).
- Acceptance: observed orders 1.9–2.1 and finest relative L2 below 0.002.
- Relative discrete energy drift below 1e-12 at n=32 and in both obstacle/vacuum controls.
- Obstacle comparison at n=64, t=0.35, disk radius 0.12, center (0.65,0.5); compact pulse radius 0.12, center (0.25,0.5). Identical initial conditions, exterior observation radius 0.18. Fixed nodes exactly zero; exterior discrepancy greater than 0.001.
- Constant solution, periodic translation symmetry, invalid CFL and input rejection.

Run python -m phispacetime.wave_benchmark for machine-readable results. These are deterministic verification cases, not general performance estimates. No absorbing boundary, driven scattering, material anisotropy, physical cloak, or electromagnetic polarization is modeled. Exterior discrepancy is a finite-time field difference, not a scattering cross section. The periodic domain is explicit; results must not be interpreted as open-space radiation.

## Next: R1b

Implement a separately verified divergence-form scalar wave solver with symmetric positive-definite anisotropic tensor and mass density. Derive transformation coefficients for the chosen scalar polarization; test manufactured variable-coefficient solutions, interface flux and mesh convergence. Compare vacuum, bare obstacle, regularized transformation and deliberately incorrect material controls, with an explicit inner-boundary model. Then evaluate an open-domain scattering configuration with documented absorbing boundaries. Do not silently substitute ray mapping for wave propagation.

R2 begins only once the same PDE and comparable error/resource metrics can be applied to uniform and adaptive spatial discretizations.
