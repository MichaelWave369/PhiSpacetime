import math
import unittest
from phispacetime.wave import Grid, solve, standing_mode, obstacle_control, relative_l2


class WaveTests(unittest.TestCase):
    def test_cfl(self):
        for dt in (0.1, 0, math.nan):
            with self.assertRaises(ValueError):
                Grid(16,dt)

    def test_constant_solution(self):
        g = Grid(8,0.05)
        field,energy = solve(g,[2.0]*64,[0.0]*64,10)
        self.assertEqual(field,[2.0]*64)
        self.assertEqual(energy["max_absolute_energy_drift"],0)
        self.assertIsNone(energy["max_relative_energy_drift"])

    def test_periodic_translation(self):
        g=Grid(8,0.05)
        f=g.sample(lambda x,y: math.sin(2*math.pi*x))
        shifted=[f[j*8+(i-1)%8] for j in range(8) for i in range(8)]
        a,_=solve(g,f,[0.0]*64,8)
        b,_=solve(g,shifted,[0.0]*64,8)
        for j in range(8):
            for i in range(8):
                self.assertAlmostEqual(b[j*8+i],a[j*8+(i-1)%8],places=12)

    def test_second_order_convergence(self):
        a,b,c = [standing_mode(n) for n in (24,48,96)]
        for p,q in ((a,b),(b,c)):
            order=math.log2(p["relative_l2"]/q["relative_l2"])
            self.assertGreater(order,1.9)
            self.assertLess(order,2.1)
        self.assertLess(c["relative_l2"],0.002)

    def test_energy(self):
        result=standing_mode(32)
        self.assertGreater(result["initial_energy"],0)
        self.assertLess(result["max_relative_energy_drift"],1e-12)

    def test_obstacle(self):
        result=obstacle_control()
        self.assertGreater(result["obstacle_nodes"],0)
        self.assertEqual(result["obstacle_max_abs"],0)
        self.assertGreater(result["exterior_relative_l2"],1e-3)
        for key in ("vacuum_energy","obstacle_energy"):
            self.assertLess(result[key]["max_relative_energy_drift"],1e-12)

    def test_invalid_fields(self):
        g=Grid(8,0.05)
        for args in (([0.0]*63,[0.0]*64,1),([0.0]*64,[0.0]*64,0),
                     ([math.nan]*64,[0.0]*64,1)):
            with self.assertRaises(ValueError):
                solve(g,*args)
        with self.assertRaises(ValueError):
            solve(g,[1.0]*64,[0.0]*64,1,[True]*64)
        with self.assertRaises(ValueError):
            relative_l2([1.0],[0.0])


if __name__ == "__main__":
    unittest.main()
