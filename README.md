# ΦSpacetime

A recursive computational laboratory for geometric physics, from PHI369 Labs.

**R0: coordinate-transformation optics and numerical refinement.** This project is not evidence of physical spacetime engineering, invisibility, faster-than-light travel, or a novel gravitational theory.

## Run locally

Requires Python 3.10+; no third-party dependencies.

```sh
python -m unittest discover -s tests -v
python -m phispacetime
```

The CLI emits JSON comparing uniform Simpson quadrature and binary/triadic adaptive Simpson quadrature against an analytic integral of a regularized radial transformation's tangential stretch. This is a scalar numerical benchmark, not an electromagnetic scattering simulation. Ray mapping is geometric optics only; the regularized inner boundary is not a perfect cloak.

See [research charter](docs/RESEARCH.md) and [source register](docs/SOURCES.md).

## Nine research tracks

Metric, Wormhole, Warp, Cloak, Vacuum, Mycelium Field, NullMesh, Altermath Physics, Validator. Only the R0 coordinate-map and quadrature primitives are implemented. Future integration must pin and verify upstream versions, not infer equivalence from names.

Reproduce. Simulate. Measure. Falsify. Refine. Repeat.
