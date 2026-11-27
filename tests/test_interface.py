import math
import unittest
from phispacetime.interface import layered_mode, layered_material, mapped_x, temporal_mode, matched_experiment, evolve_steps
from phispacetime.medium import assemble


class InterfaceTests(unittest.TestCase):
    def test_exact_flux(self):
        # Analytic one-sided physical derivatives: slope times virtual derivative.
        for x,left,right in ((.5,.5,1.5),(0.,1.5,.5)):
            derivative=2*math.pi*math.cos(2*math.pi*mapped_x(x))
            self.assertAlmostEqual((1/left)*left*derivative,(1/right)*right*derivative)
        self.assertAlmostEqual(mapped_x(.5),.25)
        self.assertAlmostEqual(mapped_x(1.),1.)
        self.assertNotEqual(layered_material(.25,0),layered_material(.75,0))

    def test_layer_convergence(self):
        rows=[layered_mode(n) for n in (16,32,64)]
        for a,b in zip(rows,rows[1:]):
            self.assertGreater(math.log2(a["amplitude_error"]/b["amplitude_error"]),1.5)
        self.assertLess(rows[-1]["amplitude_error"],.01)
        for row in rows:
            self.assertLess(row["relative_energy_drift"],1e-11)

    def test_temporal_convergence(self):
        result=temporal_mode()
        self.assertLess(result["eigen_residual"],1e-12)
        for a,b in zip(result["rows"],result["rows"][1:]):
            order=math.log2(a["error"]/b["error"])
            self.assertGreater(order,1.9)
            self.assertLess(order,2.1)

    def test_matched_controls(self):
        row=matched_experiment(32)
        self.assertGreater(row["initial_support_gap"],0.)
        self.assertGreater(row["virtual_nodes"],0)
        self.assertGreater(row["physical_nodes"],row["virtual_nodes"])
        self.assertEqual(row["cases"]["vacuum"]["vacuum_error"],0.)
        self.assertEqual(row["cases"]["virtual"]["virtual_error"],0.)
        self.assertGreater(row["cases"]["bare"]["vacuum_error"],0.)
        for case in row["cases"].values():
            self.assertEqual(case["steps"],row["steps"])
            self.assertEqual(case["dt"],row["end"]/row["steps"])
            self.assertLess(case["relative_energy_drift"],1e-11)
            self.assertTrue(math.isfinite(case["virtual_error"]))

    def test_reject_unsafe_steps(self):
        with self.assertRaises(ValueError):
            evolve_steps(assemble(8),[0.]*64,1.,2)


if __name__=="__main__":
    unittest.main()
