"""Periodic Q1 finite-element scalar wave solver with lumped mass.

rho*u_tt = div(A grad u), symmetric positive definite A, positive rho.
Two-point Gauss quadrature per axis. All quantities nondimensional.
"""
import math
from dataclasses import dataclass
from .optics import Annulus
from .wave import relative_l2


def vacuum(x, y):
    return 1., 0., 1., 1.


def radial_medium(inner, outer, virtual_inner, center=(.5,.5)):
    ann = Annulus(inner, outer, virtual_inner)
    cx, cy = center

    def material(x, y):
        dx, dy = x-cx, y-cy
        R = math.hypot(dx,dy)
        # Core material is a numerical extension, not a physical cloak.
        # The fixed-node mask approximates a Dirichlet inner boundary.
        if R < inner or R >= outer:
            return vacuum(x,y)
        r = ann.inverse(R)
        a, b, _ = ann.stretches(r)
        det = a*b
        ar, at = a*a/det, b*b/det
        nx, ny = dx/R, dy/R
        return (ar*nx*nx+at*ny*ny, (ar-at)*nx*ny,
                ar*ny*ny+at*nx*nx, 1/det)

    return material


@dataclass
class Operator:
    n: int
    mass: list
    rows: list

    def apply(self, field):
        return [math.fsum(a*field[j] for j,a in row) for row in self.rows]

    @property
    def dt_limit(self):
        # Gershgorin upper bound on lambda_max(M^-1 K).
        bound = max(math.fsum(abs(a) for _,a in row)/m
                    for row,m in zip(self.rows,self.mass))
        return 2/math.sqrt(bound)


def assemble(n, material=vacuum):
    if not isinstance(n,int) or n < 4:
        raise ValueError("n must be an integer >= 4")
    h=1/n
    rows=[{} for _ in range(n*n)]
    mass=[0.]*(n*n)
    g=1/math.sqrt(3)
    signs=((-1,-1),(1,-1),(1,1),(-1,1))
    for iy in range(n):
        for ix in range(n):
            ids=[iy*n+ix,iy*n+(ix+1)%n,
                 ((iy+1)%n)*n+(ix+1)%n,((iy+1)%n)*n+ix]
            for xi in (-g,g):
                for eta in (-g,g):
                    x=(ix+(xi+1)/2)*h
                    y=(iy+(eta+1)/2)*h
                    ax,b,ay,rho=material(x,y)
                    if (not all(math.isfinite(v) for v in (ax,b,ay,rho))
                            or ax<=0 or ay<=0 or ax*ay-b*b<=0 or rho<=0):
                        raise ValueError("material must be finite, SPD, positive density")
                    shape=[(1+s*xi)*(1+t*eta)/4 for s,t in signs]
                    gx=[s*(1+t*eta)/(2*h) for s,t in signs]
                    gy=[t*(1+s*xi)/(2*h) for s,t in signs]
                    w=h*h/4
                    for p,i in enumerate(ids):
                        mass[i]+=w*rho*shape[p]
                        for q,j in enumerate(ids):
                            val=w*(ax*gx[p]*gx[q]+b*(gx[p]*gy[q]+gy[p]*gx[q])
                                   +ay*gy[p]*gy[q])
                            rows[i][j]=rows[i].get(j,0.)+val
    return Operator(n,mass,[sorted(row.items()) for row in rows])


def evolve(op, initial, velocity, end, fixed=None, safety=.65):
    count=op.n**2
    if (not math.isfinite(end) or end<=0 or not 0<safety<1
            or len(initial)!=count or len(velocity)!=count
            or not all(math.isfinite(v) for v in (*initial,*velocity))):
        raise ValueError("invalid evolution parameters")
    mask=[False]*count if fixed is None else list(fixed)
    if len(mask)!=count or any(type(v) is not bool for v in mask):
        raise ValueError("invalid fixed mask")
    if any(m and (u!=0 or v!=0) for m,u,v in zip(mask,initial,velocity)):
        raise ValueError("incompatible fixed-node data")
    steps=max(1,math.ceil(end/(safety*op.dt_limit)))
    dt=end/steps
    prev=list(initial)
    k=op.apply(prev)
    curr=[0. if mask[i] else u+dt*v-.5*dt**2*k[i]/op.mass[i]
          for i,(u,v) in enumerate(zip(initial,velocity))]

    def energy(a,b):
        kb=op.apply(b)
        return .5*math.fsum(m*((v-u)/dt)**2+u*z
                            for m,u,v,z in zip(op.mass,a,b,kb))

    e0=energy(prev,curr)
    drift=0.
    for _ in range(1,steps):
        k=op.apply(curr)
        nxt=[0. if mask[i] else 2*curr[i]-prev[i]-dt**2*k[i]/op.mass[i]
             for i in range(count)]
        prev,curr=curr,nxt
        e=energy(prev,curr)
        if not math.isfinite(e):
            raise RuntimeError("nonfinite discrete energy")
        drift=max(drift,abs(e-e0))
    return curr,dict(steps=steps,dt=dt,dt_limit=op.dt_limit,energy=e0,
                     relative_energy_drift=drift/abs(e0) if e0 else None)


def sample(n,f):
    return [f(i/n,j/n) for j in range(n) for i in range(n)]


def manufactured(n, variable=False, end=.13):
    if variable:
        # rho=Axx varies transversely; u depends only on x.
        def material(x,y):
            r=1+.3*math.cos(2*math.pi*y)
            return r,0.,2.,r
        init=sample(n,lambda x,y: math.sin(2*math.pi*x))
        omega=2*math.pi
    else:
        def material(x,y):
            return 2.,.35,1.,1.
        init=sample(n,lambda x,y: math.sin(2*math.pi*(x+y)))
        omega=2*math.pi*math.sqrt(2+2*.35+1)
    op=assemble(n,material)
    out,info=evolve(op,init,[0.]*(n*n),end)
    exact=[u*math.cos(omega*end) for u in init]
    return dict(n=n,amplitude_normalized_l2=math.sqrt(
        math.fsum((a-b)**2 for a,b in zip(out,exact))/math.fsum(u*u for u in init)),
        instantaneous_relative_l2=relative_l2(out,exact),**info)


def experiment(n=48,end=.27):
    """Exploratory finite-time scalar exterior comparison, no cloak claim."""
    inner,outer,virtual=.12,.22,.04
    radius=sample(n,lambda x,y: math.hypot(x-.5,y-.5))
    fixed=[v<=inner for v in radius]
    exterior=[v>=outer+.03 for v in radius]
    def pulse(x,y):
        q=((x-.15)**2+(y-.5)**2)/.09**2
        return max(0.,1-q)**4
    initial=sample(n,pulse)
    zeros=[0.]*(n*n)
    base=assemble(n)
    transformed=radial_medium(inner,outer,virtual)
    def wrong(x,y):
        ax,b,ay,rho=transformed(x,y)
        return 1.,0.,1.,rho
    ref,info=evolve(base,initial,zeros,end)
    cases={"vacuum":dict(exterior_error=0.,**info)}
    for name,operator in (("bare",base),("transformed",assemble(n,transformed)),
                          ("density_only",assemble(n,wrong))):
        out,diag=evolve(operator,initial,zeros,end,fixed)
        cases[name]=dict(exterior_error=relative_l2(out,ref,exterior),
                         fixed_max=max(abs(u) for u,m in zip(out,fixed) if m),**diag)
    # Small virtual obstacle control, approximating the regularized map's
    # reference problem; geometry and discretization errors remain.
    small=[v<=virtual for v in radius]
    out,diag=evolve(base,initial,zeros,end,small)
    cases["virtual_obstacle"]=dict(exterior_error=relative_l2(out,ref,exterior),**diag)
    return dict(n=n,end=end,radii=[inner,outer,virtual],cases=cases)
