"""Deterministic uniform and adaptive Simpson rules, with explicit budget failure."""
import math


def _validate(a, b):
    if not (math.isfinite(a) and math.isfinite(b) and a < b):
        raise ValueError("finite increasing bounds required")


def uniform(f, a, b, panels):
    _validate(a, b)
    if not isinstance(panels, int) or panels < 2 or panels % 2:
        raise ValueError("positive even panel count required")
    h = (b - a) / panels
    values = [f(a + i * h) for i in range(panels + 1)]
    if not all(map(math.isfinite, values)):
        raise ValueError("nonfinite integrand")
    return h / 3 * (values[0] + values[-1] + sum(
        (4 if i % 2 else 2) * values[i] for i in range(1, panels))), len(values)


def adaptive(f, a, b, tolerance=1e-8, arity=2, max_depth=20, max_evals=100000):
    _validate(a, b)
    if (arity not in (2, 3) or not math.isfinite(tolerance) or tolerance <= 0
            or not isinstance(max_depth, int) or max_depth < 0 or max_depth > 64
            or not isinstance(max_evals, int) or max_evals < 3):
        raise ValueError("invalid refinement parameters")
    cache = {}

    def sample(x):
        if x not in cache:
            if len(cache) >= max_evals:
                raise RuntimeError("evaluation budget exhausted")
            value = f(x)
            if not math.isfinite(value):
                raise ValueError("nonfinite integrand")
            cache[x] = value
        return cache[x]

    def simpson(l, r):
        m = (l + r) / 2
        return (r - l) / 6 * (sample(l) + 4 * sample(m) + sample(r))

    def refine(l, r, coarse, tol, depth):
        bounds = [l] + [l + (r - l) * j / arity for j in range(1, arity)] + [r]
        segments = list(zip(bounds, bounds[1:]))
        fine = [simpson(x, y) for x, y in segments]
        total = math.fsum(fine)
        # Simpson error scales as h^4 for sufficiently smooth integrands.
        correction = (total - coarse) / (arity ** 4 - 1)
        if abs(correction) <= tol:
            return total + correction
        if depth >= max_depth:
            raise RuntimeError("refinement depth exhausted")
        return math.fsum(refine(x, y, s, tol / arity, depth + 1)
                         for (x, y), s in zip(segments, fine))

    result = refine(a, b, simpson(a, b), tolerance, 0)
    return result, len(cache)
