import math
import unittest
from phispacetime.medium import assemble
from phispacetime.absorbing import sponge, coordinates, wall_mask, evolve_damped, reflection, transformed_comparison


class AbsorbingTests(unittest.TestCase):
    def test_nonperiodic_assembly(self):
        op=assemble(9,periodic=False,length=2.,origin=(-.5,-.5))
        self.assertAlmostEqual(sum(op.mass),4.)
        self.assertNotIn(8,dict(op.rows[0]))
        self.assertNotIn(72,dict(op.rows[0]))
        for i,row in enumerate(op.rows):
            self.assertAlmostEqual(sum(v for _,v in row),0.,places=12)
            for j,v in row:
                self.assertAlmostEqual(v,dict(op.rows[j])[i],places=12)

    def test_sponge(self):
        self.assertEqual(sponge(.5,.5),0.)
        self.assertEqual(sponge(0.,.5),35.)
        self.assertGreater(sponge(.05,.5),sponge(.1,.5))
        with self.assertRaises(ValueError):
            sponge(.5,.5,strength=-1.)

    def test_energy_balance(self):
        n=17
        op=assemble(n,periodic=False)
        xy=coordinates(n)
        mask=wall_mask(n)
        u=[0. if m else math.sin(math.pi*x)*math.sin(math.pi*y)
           for (x,y),m in zip(xy,mask)]
        gamma=[sponge(x,y) for x,y in xy]
        out,d=evolve_damped(op,u,[0.]*(n*n),.5,32,gamma,mask)
        self.assertGreater(d["dissipated"],0.)
        self.assertLess(d["final_energy"],d["initial_energy"])
        self.assertLess(d["relative_balance_error"],1e-10)
        self.assertLess(d["max_relative_energy_increase"],1e-10)
        self.assertTrue(all(v==0 for v,m in zip(out,mask) if m))
        with self.assertRaises(ValueError):
            evolve_damped(op,u,[0.]*(n*n),.5,2,gamma,mask)

    def test_reference_comparison(self):
        row=reflection(25)
        self.assertGreater(row["probes"],0)
        self.assertEqual(row["reference"]["dissipated"],0.)
        for case in row["cases"].values():
            self.assertEqual(case["steps"],row["reference"]["steps"])
            self.assertLess(case["relative_balance_error"],1e-10)
            self.assertTrue(math.isfinite(case["return_error"]))
        # A functioning absorber must improve this fixed reference case.
        self.assertLess(row["cases"]["sponge"]["return_error"],
                        row["cases"]["reflecting"]["return_error"])

    def test_transformed_controls(self):
        row=transformed_comparison(25)
        self.assertGreater(row["probe_count"],0)
        self.assertEqual(row["cases"]["vacuum"]["exterior_error"],0.)
        for case in row["cases"].values():
            self.assertLess(case["relative_balance_error"],1e-10)
            self.assertTrue(math.isfinite(case["exterior_error"]))


if __name__=="__main__":
    unittest.main()
