"""Emit R0 benchmark results; no Maxwell or gravity solve."""
import json
from .optics import Annulus
from .quadrature import adaptive, uniform


def benchmark():
    ann = Annulus(1.0, 3.0, 0.1)
    exact = ann.stretch_integral()
    f = lambda r: ann.stretches(r)[1]
    rows = []
    for tolerance in (1e-4, 1e-6, 1e-8):
        for arity in (2, 3):
            estimate, n = adaptive(f, ann.virtual_inner, ann.outer, tolerance, arity)
            rows.append(dict(method=f"adaptive-{arity}", tolerance=tolerance,
                             value=estimate, evaluations=n, abs_error=abs(estimate-exact)))
        # Uniform control: first doubling that meets the same measured error target.
        panels = 2
        total = 0
        while panels <= 131072:
            estimate, n = uniform(f, ann.virtual_inner, ann.outer, panels)
            total += n
            if abs(estimate-exact) <= tolerance:
                break
            panels *= 2
        else:
            raise RuntimeError("uniform benchmark budget exhausted")
        rows.append(dict(method="uniform", tolerance=tolerance, value=estimate,
                         evaluations=n, cumulative_search_evaluations=total,
                         abs_error=abs(estimate-exact)))
    return dict(scope="scalar stretch integral; not electromagnetic scattering",
                radii=dict(inner=ann.inner, outer=ann.outer, virtual_inner=ann.virtual_inner),
                exact=exact, rows=rows)


if __name__ == "__main__":
    print(json.dumps(benchmark(), indent=2, allow_nan=False))
