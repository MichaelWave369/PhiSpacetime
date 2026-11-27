# R1b: transformed scalar media

Solve rho u_tt = div(A grad u) on a periodic unit square. Q1 finite elements, tensor-product two-point Gauss quadrature, row-sum lumped mass, centered explicit leapfrog with Taylor startup. The stiffness matrix is symmetric; the constant mode has zero stiffness. Dirichlet obstacles constrain nodes, not exact curved boundaries. The discrete half-step energy uses the cross-time stiffness bilinear form. Gershgorin bounds the maximum eigenvalue of M^-1 K; time steps use 0.65 of the resulting stability limit or smaller.

For x=F(X), with J=dx/dX and virtual coefficients equal to unity, the weak-form pushforward is A=JJ^T/det J and rho=1/det J in two dimensions. For the R0 annular mapping the radial and tangential stretches are a=dR/dr and b=R/r. Thus principal A entries are a/b,b/a, rotated into Cartesian axes, and rho=1/(ab). This is a scalar transformed-wave equation, not a vector Maxwell constitutive model. The material is regularized by a finite virtual inner radius; vacuum extension inside the core plus masked nodes is an approximation, not exact enforcement of the curved interface.

## Verification gates

- Matrix symmetry, zero row sum, positive mass and known integral of density.
- Reject non-SPD tensors and nonpositive or nonfinite density.
- Constant field and stable time step.
- Two analytic periodic solutions: constant off-diagonal anisotropy and transversely varying mass with Axx=rho(y). Resolutions 16,32,64. Observed initial-amplitude-normalized L2 error order in [1.7,2.3], finest error below .01 and discrete energy drift below 1e-11. The anisotropic standing mode at t=.13 is near a zero crossing; normalizing error by the initial mode norm avoids an ill-conditioned instantaneous denominator. Both error measures are reported; the original instantaneous relative criterion failed and is not a robust amplitude accuracy metric at this time.
- Radial coefficient rotation, determinant and exterior identity.
- Finite-time exterior comparisons against vacuum for bare disk, transformed disk, density-only negative control and smaller virtual disk; energy and fixed-node checks. No required ranking among the controls.

Run python -m phispacetime.medium_benchmark. The default exploratory comparison uses n=48, end=.27, a=.12, b=.22, epsilon=.04; a compact pulse centered at (.15,.5), radius .09. Observe nodes outside radius .25. Different coefficient fields impose different stable time steps; all runs end at the same physical time but temporal discretization errors need separate convergence tests. No superiority assertion follows from one resolution or pulse. A density-only control is not an optimized conventional cloak.

## Next gate

R1c: resolve curved interfaces, verify interface flux using manufactured piecewise coefficients and independent spatial/time refinement; compare the transformed medium against its virtual obstacle with matching outer observations. Repeat across mesh sizes, pulse profiles, and times, and publish unfavorable outcomes. Add a separately verified open-boundary formulation before interpreting results as free-space scattering. R2 adaptive nine-way pockets must use the same PDE, error criteria and resource accounting. No existing Phi repository has been imported.
