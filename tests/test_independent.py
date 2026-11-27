import math
import unittest
from phispacetime.independent import independent_nodes, interpolate, compare
from phispacetime.annular import mesh, assemble_annulus, harmonic


class IndependentTests(unittest.TestCase):
    def test_independent_geometry(self):
        v,p,a=mesh(17)
        r,phase=independent_nodes(17,"transformed")
        self.assertNotEqual(r,p)
        self.assertNotEqual(phase,0)
        self.assertAlmostEqual(r[8],a.outer)
        self.assertNotEqual(r[10],v[10])
        with self.assertRaises(ValueError):
            assemble_annulus(17,"transformed",list(reversed(r)))

    def test_interpolation(self):
        r,phase=independent_nodes(17,"transformed")
        n=len(r)
        f=[2*x+1 for x in r for j in range(n)]
        self.assertAlmostEqual(interpolate(f,r,phase,.31,-.01),1.62)
        self.assertAlmostEqual(interpolate(f,r,phase,.31,2*math.pi-.01),1.62)
        for i in (0,5,16):
            self.assertAlmostEqual(interpolate(f,r,phase,r[i],phase),f[i*n])

    def test_static_mode(self):
        for case in ("virtual","transformed","bare"):
            r,phase=independent_nodes(33,case)
            _,d=harmonic(33,case,r,phase,2)
            self.assertLess(d["relative_l2"],.01)
            self.assertLess(d["relative_residual"],1e-9)

    def test_common_probes(self):
        result=compare(17)
        self.assertEqual(result["probe_count"],259)
        self.assertEqual(result["cases"]["virtual"]["wave_error"],0.)
        self.assertEqual(result["cases"]["virtual"]["static_error"],0.)
        for c in result["cases"].values():
            self.assertEqual(c["wave"]["steps"],result["steps"])
            self.assertLess(c["wave"]["relative_energy_drift"],1e-10)
            self.assertTrue(math.isfinite(c["flux"]["jump"]))
        self.assertGreater(result["cases"]["bare"]["static_error"],.001)

    def test_static_refinement(self):
        rows=[compare(n) for n in (17,33)]
        for name in ("virtual","transformed"):
            a,b=[r["cases"][name]["static"]["relative_l2"] for r in rows]
            self.assertLess(b,a)
        self.assertLess(rows[-1]["cases"]["transformed"]["static_error"],.01)


if __name__=="__main__":
    unittest.main()
