import math
import unittest
from phispacetime.optics import Annulus
from phispacetime.quadrature import adaptive, uniform
from phispacetime.__main__ import benchmark


class OpticsTests(unittest.TestCase):
    def setUp(self):
        self.a = Annulus(1, 3, 0.1)

    def test_invalid(self):
        for args in [(1, 3, 0), (1, 1, .1), (1, 3, 2), (math.nan, 3, .1)]:
            with self.assertRaises(ValueError):
                Annulus(*args)

    def test_boundaries(self):
        self.assertAlmostEqual(self.a.radius(.1), 1)
        self.assertEqual(self.a.radius(3), 3)
        self.assertEqual(self.a.radius(4), 4)
        self.assertAlmostEqual(self.a.radius(3-1e-8), 3, places=7)

    def test_inverse(self):
        for r in (.1, .2, 1, 2, 3, 4):
            self.assertAlmostEqual(self.a.inverse(self.a.radius(r)), r)

    def test_jacobian(self):
        r, h = .8, 1e-6
        derivative = (self.a.radius(r+h)-self.a.radius(r-h))/(2*h)
        self.assertAlmostEqual(derivative, self.a.stretches(r)[0], places=8)
        self.assertAlmostEqual(self.a.stretches(r)[1], self.a.radius(r)/r)

    def test_tensor(self):
        for r in (.1, .5, 2, 3, 4):
            er, et, ez = self.a.relative_tensor(r)
            self.assertTrue(all(math.isfinite(t) and t > 0 for t in (er, et, ez)))
            self.assertAlmostEqual(er*et, 1)
        self.assertEqual(self.a.relative_tensor(4), (1, 1, 1))

    def test_ray(self):
        points = self.a.mapped_ray(.25, 5, 101)
        self.assertEqual(points[0], (-5, .25))
        self.assertEqual(points[-1], (5, .25))
        self.assertTrue(all(math.hypot(x, y) >= 1 for x, y in points))
        with self.assertRaises(ValueError):
            self.a.mapped_ray(.01, 5, 101)


class QuadratureTests(unittest.TestCase):
    def test_polynomial(self):
        self.assertAlmostEqual(uniform(lambda x: x**3, 0, 1, 10)[0], .25)
        for arity in (2, 3):
            self.assertAlmostEqual(adaptive(lambda x: x**3, 0, 1, arity=arity)[0], .25)

    def test_budget_and_invalid(self):
        with self.assertRaises(RuntimeError):
            adaptive(math.exp, 0, 1, max_evals=3)
        with self.assertRaises(ValueError):
            uniform(math.exp, 0, 1, 3)
        with self.assertRaises(ValueError):
            adaptive(lambda x: math.nan, 0, 1)

    def test_benchmark(self):
        result = benchmark()
        self.assertEqual(len(result["rows"]), 9)
        for row in result["rows"]:
            self.assertLessEqual(row["abs_error"], row["tolerance"])
            self.assertGreater(row["evaluations"], 0)


if __name__ == "__main__":
    unittest.main()
