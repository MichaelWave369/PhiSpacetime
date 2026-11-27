import math
import unittest
from collections import Counter
from phispacetime.pockets import leaves,triangulate,assemble,geometry,run,exact,source


class PocketTests(unittest.TestCase):
    def test_partition_and_budget(self):
        for arity in (2,3):
            scale,cells=leaves(arity,100)
            self.assertLessEqual(len(cells),100)
            self.assertEqual(sum(c.h*c.h for c in cells),scale*scale)
            self.assertEqual((len(cells)-1)%(arity**2-1),0)
            self.assertEqual((scale,cells),leaves(arity,100))

    def test_conforming_edges(self):
        for arity in (2,3):
            scale,cells=leaves(arity,100)
            xy,tri,fixed=triangulate(scale,cells)
            edges=Counter(tuple(sorted((a,b))) for t in tri for a,b in
                          zip(t,t[1:]+t[:1]))
            self.assertAlmostEqual(sum(geometry(xy,t)[0] for t in tri),1.)
            for (i,j),count in edges.items():
                a,b=xy[i],xy[j]
                boundary=any(a[k]==b[k] and a[k] in (0.,1.) for k in (0,1))
                self.assertEqual(count,1 if boundary else 2)

    def test_operator(self):
        scale,cells=leaves(3,64)
        xy,tri,fixed=triangulate(scale,cells)
        rows,mass,load=assemble(xy,tri)
        self.assertAlmostEqual(sum(mass),1.)
        self.assertTrue(all(v>0 for v in mass))
        for i,row in enumerate(rows):
            self.assertAlmostEqual(sum(v for _,v in row),0.,places=10)
            for j,v in row:
                self.assertAlmostEqual(v,dict(rows[j])[i],places=10)
            if not fixed[i]:
                self.assertAlmostEqual(sum(v*xy[j][0] for j,v in row),0.,places=10)

    def test_manufactured_source(self):
        h=1e-4
        for x,y in ((.36,.61),(.45,.57),(.1,.2)):
            u=exact(x,y)
            lap=(exact(x+h,y)+exact(x-h,y)+exact(x,y+h)+exact(x,y-h)-4*u)/h**2
            self.assertAlmostEqual(-lap,source(x,y),delta=0.001)

    def test_uniform_convergence(self):
        a,b,c=[run(2,k,False,False) for k in (16,64,256)]
        self.assertLess(b["relative_l2"],a["relative_l2"])
        self.assertLess(c["relative_l2"],b["relative_l2"])
        self.assertGreater(math.log2(b["relative_l2"]/c["relative_l2"]),1.5)
        self.assertLess(c["relative_residual"],1e-9)

    def test_adaptive_controls(self):
        for arity in (2,3):
            a,b=[run(arity,k) for k in (64,256)]
            self.assertLess(b["relative_l2"],a["relative_l2"])
            self.assertLess(b["relative_residual"],1e-9)
            self.assertAlmostEqual(b["area"],1.)
            self.assertGreater(b["free_dofs"],0)


if __name__=="__main__":
    unittest.main()
