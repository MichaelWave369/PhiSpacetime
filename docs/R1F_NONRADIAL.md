# R1f: compact nonradial coordinate map

R1e independent fitted polar grids showed decreasing transformed exterior wave error (.1490, .04349, .00927 at n=17,33,65). Later observation at n=33 gave .05697 transformed versus .48280 bare; not all settings favor transformation. Source: https://github.com/MichaelWave369/PhiSpacetime/actions/runs/45921209815.

R1f removes radial symmetry, not the scalar-model restriction. On a periodic unit square prescribe F(x,y)=(x,y+a b(x)b(y)), b=sin^4(2 pi (x-.25)) in (.25,.75), zero elsewhere. F is identity outside the support square, has a positive determinant for the allowed amplitude bound, and is inverted by monotone bisection. J=[[1,0],[s,d]], A=JJ^T/d, rho=1/d. The physical and virtual Cartesian grids use the same fixed coordinate locations, not mapped node pairs. The coefficients are smooth enough for the verification problem; no singular map, obstacle or cloak is present.

Analytic reference: sin(2 pi x)cos(2 pi y)cos(2 pi sqrt(2)t) evaluated at F inverse; zero initial velocity. Compare initial-amplitude-normalized L2 error at 16,32,64, requiring order >1.4 and finest error <.02; energy drift <1e-10. Verify Jacobian by finite differences, inverse round trip, determinant, SPD and identity support. These checks are numerical verification, not new physics.

Exterior experiments: identical compact pulse outside the deformation, vacuum, full transformed medium and density-only negative control, shared time steps, exterior nodes with buffer .03. Report six settings including opposite deformation sign and changed pulse location, half time step and later observation. No ranking gate. The zero-amplitude control must yield exact agreement. The finite periodic domain can exhibit reentry; no free-space scattering interpretation. Results are scalar, not vector Maxwell. No NBG integration yet.

Run python -m phispacetime.nonradial_benchmark. Next: absorbing/open-domain verification, then a controlled adaptive mesh study with equal PDE and error metrics. Preserve earlier unfavorable controls. This is a coordinate-deformation test, not a noncircular obstacle cloak.
