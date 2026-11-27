"""Regularized annular radial map. Lengths share arbitrary consistent units."""
from dataclasses import dataclass
import math


@dataclass(frozen=True)
class Annulus:
    inner: float
    outer: float
    virtual_inner: float

    def __post_init__(self):
        if not all(math.isfinite(v) for v in (self.inner, self.outer, self.virtual_inner)):
            raise ValueError("finite radii required")
        if not 0 < self.virtual_inner < self.inner < self.outer:
            raise ValueError("require 0 < virtual_inner < inner < outer")

    @property
    def slope(self):
        return (self.outer - self.inner) / (self.outer - self.virtual_inner)

    def radius(self, r):
        if not math.isfinite(r) or r < self.virtual_inner:
            raise ValueError("outside mapped domain")
        if r >= self.outer:
            return r
        return self.inner + self.slope * (r - self.virtual_inner)

    def inverse(self, R):
        if not math.isfinite(R) or R < self.inner:
            raise ValueError("outside physical mapped domain")
        if R >= self.outer:
            return R
        return self.virtual_inner + (R - self.inner) / self.slope

    def stretches(self, r):
        R = self.radius(r)
        # Interior one-sided value at outer interface; exterior is identity.
        if r > self.outer:
            return (1.0, 1.0, 1.0)
        return (self.slope, R / r, 1.0)

    def relative_tensor(self, r):
        """Diagonal epsilon_r and mu_r in physical orthonormal radial basis.
        J J^T / det(J), virtual vacuum, no material realization implied.
        """
        a, b, c = self.stretches(r)
        determinant = a * b * c
        return tuple(s * s / determinant for s in (a, b, c))

    def map_point(self, x, y):
        r = math.hypot(x, y)
        R = self.radius(r)
        return (x * R / r, y * R / r)

    def mapped_ray(self, impact, extent, samples):
        """Map a horizontal virtual ray; reject rays entering excluded core."""
        if (not math.isfinite(impact) or not math.isfinite(extent)
                or abs(impact) < self.virtual_inner or extent <= self.outer
                or not isinstance(samples, int) or samples < 2):
            raise ValueError("ray intersects excluded core or invalid sampling")
        return [self.map_point(-extent + 2 * extent * i / (samples - 1), impact)
                for i in range(samples)]

    def stretch_integral(self):
        """Integral of tangential stretch over virtual annulus."""
        k = self.inner - self.slope * self.virtual_inner
        return self.slope * (self.outer - self.virtual_inner) + k * math.log(
            self.outer / self.virtual_inner)
