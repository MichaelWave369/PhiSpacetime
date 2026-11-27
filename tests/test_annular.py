import math
import unittest
from phispacetime.annular import mesh, assemble_annulus, harmonic, comparison


class AnnularTests(unittest.TestCase):
    def test_mesh(self):
        v,p,a=mesh(17)
        self.assertAlmostEqual(v[0],.04)
        self.assertAlmostEqual(p[0],.12)
        self.assertAlmostEqual(v[8],.22)
        self.assertEqual(v[8:],p[8:])
        self.assertAlmostEqual(p[-1],.45)
        with self.assertRaises(ValueError):
            mesh(16)

    def test_mass_and_symmetry(self):
        for case in ("virtual","transformed","bare"):
            op,r=assemble_annulus(17,case)
            self.assertTrue(all(m>0 for m in op.mass))
            inner=.12 if case=="bare" else .04
            self.assertAlmostEqual(sum(op.mass),math.pi*(.45**2-inner**2),places=12)
            for i,row in enumerate(op.rows):
                self.assertAlmostEqual(sum(v for _,v in row),0.,places=11)
                for j,v in row:
                    self.assertAlmostEqual(v,dict(op.rows[j])[i],places=12)

    def test_harmonic_convergence(self):
        errors=[harmonic(n,"virtual")[1]["relative_l2"] for n in (17,33,65)]
        for a,b in zip(errors,errors[1:]):
            self.assertGreater(math.log2(a/b),1.5)
        self.assertLess(errors[-1],.002)

    def test_paired_covariance(self):
        result=comparison(17)
        self.assertLess(result["max_matrix_difference"],1e-11)
        self.assertLess(result["max_mass_difference"],1e-13)
        t=result["cases"]["transformed"]
        self.assertLess(t["wave_exterior_error"],1e-10)
        self.assertLess(t["static_exterior_error"],1e-9)
        self.assertGreater(result["cases"]["bare"]["static_exterior_error"],.001)
        for case in result["cases"].values():
            self.assertLess(case["wave"]["relative_energy_drift"],1e-10)
            self.assertLess(case["static"]["relative_residual"],1e-9)


if __name__=="__main__":
    unittest.main()
