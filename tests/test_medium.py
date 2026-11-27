import math
import unittest
from phispacetime.medium import assemble, evolve, manufactured, radial_medium, experiment


class MediumTests(unittest.TestCase):
    def test_matrix(self):
        op=assemble(8,lambda x,y:(2.,.3,1.,1.+.2*y))
        self.assertAlmostEqual(sum(op.mass),1.1)
        for i,row in enumerate(op.rows):
            self.assertAlmostEqual(sum(v for _,v in row),0.,places=12)
            for j,v in row:
                self.assertAlmostEqual(v,dict(op.rows[j])[i],places=12)
        self.assertTrue(all(m>0 for m in op.mass))

    def test_invalid_material(self):
        for value in ((1.,2.,1.,1.),(1.,0.,1.,0.),(math.nan,0.,1.,1.)):
            with self.assertRaises(ValueError):
                assemble(4,lambda x,y:value)

    def test_constant(self):
        op=assemble(8)
        u,diag=evolve(op,[2.]*64,[0.]*64,.1)
        self.assertLess(max(abs(v-2.) for v in u),1e-12)
        self.assertLess(diag["dt"],diag["dt_limit"])

    def test_manufactured(self):
        for variable in (False,True):
            a,b,c=[manufactured(n,variable) for n in (16,32,64)]
            for p,q in ((a,b),(b,c)):
                order=math.log2(p["amplitude_normalized_l2"]/q["amplitude_normalized_l2"])
                self.assertGreater(order,1.7)
                self.assertLess(order,2.3)
            self.assertLess(c["amplitude_normalized_l2"],.01)
            self.assertLess(c["relative_energy_drift"],1e-11)

    def test_transform(self):
        mat=radial_medium(.12,.22,.04)
        self.assertEqual(mat(.9,.5),(1.,0.,1.,1.))
        a,b,d,rho=mat(.66,.5)
        self.assertAlmostEqual(b,0.)
        self.assertAlmostEqual(a*d-b*b,1.)
        self.assertGreater(rho,0.)
        aa,bb,dd,rr=mat(.5,.66)
        self.assertAlmostEqual(aa,d)
        self.assertAlmostEqual(dd,a)
        self.assertAlmostEqual(rr,rho)

    def test_controls(self):
        result=experiment(32)
        cases=result["cases"]
        self.assertEqual(cases["vacuum"]["exterior_error"],0.)
        for name in ("bare","transformed","density_only"):
            self.assertEqual(cases[name]["fixed_max"],0.)
            self.assertGreater(cases[name]["exterior_error"],0.)
        for row in cases.values():
            self.assertTrue(math.isfinite(row["exterior_error"]))
            self.assertLess(row["relative_energy_drift"],1e-11)


if __name__=="__main__":
    unittest.main()
