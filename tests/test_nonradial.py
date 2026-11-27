import math
import unittest
from phispacetime.nonradial import Deformation, manufactured, comparison


class NonradialTests(unittest.TestCase):
    def test_inverse_and_identity(self):
        for a in (-.035,0.,.035):
            d=Deformation(a)
            for x in (.1,.25,.37,.5,.73,.9):
                for y in (.1,.25,.39,.5,.7,.9):
                    X,Y=d.forward(x,y)
                    u,v=d.inverse(X,Y)
                    self.assertAlmostEqual(u,x,places=14)
                    self.assertAlmostEqual(v,y,places=14)
                    if x<=.25 or x>=.75 or y<=.25 or y>=.75:
                        self.assertEqual((X,Y),(x,y))
                        self.assertEqual(d.material(X,Y),(1.,0.,1.,1.))
        with self.assertRaises(ValueError):
            Deformation(.1)

    def test_jacobian_and_spd(self):
        d=Deformation()
        h=1e-6
        for x,y in ((.38,.41),(.59,.65),(.5,.37)):
            s,j=d.jacobian(x,y)
            self.assertGreater(j,0.)
            self.assertAlmostEqual(s,(d.forward(x+h,y)[1]-d.forward(x-h,y)[1])/(2*h),places=8)
            self.assertAlmostEqual(j,(d.forward(x,y+h)[1]-d.forward(x,y-h)[1])/(2*h),places=8)
            a,b,c,rho=d.material(*d.forward(x,y))
            self.assertGreater(a,0.)
            self.assertGreater(rho,0.)
            self.assertAlmostEqual(a*c-b*b,1.,places=12)

    def test_manufactured(self):
        rows=[manufactured(n) for n in (16,32,64)]
        for a,b in zip(rows,rows[1:]):
            self.assertGreater(math.log2(a["amplitude_error"]/b["amplitude_error"]),1.4)
        self.assertLess(rows[-1]["amplitude_error"],.02)
        for row in rows:
            self.assertLess(row["relative_energy_drift"],1e-10)

    def test_exterior_controls(self):
        row=comparison(24)
        self.assertGreater(row["exterior_nodes"],0)
        self.assertEqual(row["cases"]["vacuum"]["exterior_error"],0.)
        for d in row["cases"].values():
            self.assertEqual(d["steps"],row["steps"])
            self.assertLess(d["relative_energy_drift"],1e-10)
            self.assertTrue(math.isfinite(d["exterior_error"]))

    def test_identity_control(self):
        row=comparison(24,amplitude=0.)
        for d in row["cases"].values():
            self.assertEqual(d["exterior_error"],0.)


if __name__=="__main__":
    unittest.main()
