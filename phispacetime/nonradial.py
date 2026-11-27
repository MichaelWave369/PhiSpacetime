"""Compactly supported nonradial coordinate deformation of a periodic square.

This is a smooth invertible scalar-wave coordinate map, without an obstacle
or singular cloak. The physical and virtual Cartesian grids are not mapped
to one another. All quantities nondimensional.
"""
from dataclasses import dataclass
import math
from .medium import assemble, sample, vacuum
from .interface import evolve_steps
from .wave import relative_l2


def bump(x):
    if not .25 < x < .75:
        return 0.,0.
    t=2*math.pi*(x-.25)
    # sin^4 has three continuous derivatives at the support boundary.
    return math.sin(t)**4,8*math.pi*math.sin(t)**3*math.cos(t)


@dataclass(frozen=True)
class Deformation:
    amplitude: float=.035

    def __post_init__(self):
        # Conservative global bound |bump'| <= 8*pi.
        if not math.isfinite(self.amplitude) or abs(self.amplitude)*8*math.pi>=.95:
            raise ValueError("amplitude violates invertibility bound")

    def forward(self,x,y):
        bx,_=bump(x)
        by,_=bump(y)
        return x,y+self.amplitude*bx*by

    def jacobian(self,x,y):
        bx,dx=bump(x)
        by,dy=bump(y)
        return self.amplitude*dx*by,1+self.amplitude*bx*dy

    def inverse(self,X,Y):
        # Monotone vertical map; identity outside the support.
        if not (.25<X<.75 and .25<Y<.75):
            return X,Y
        lo,hi=.25,.75
        for _ in range(53):
            mid=(lo+hi)/2
            if self.forward(X,mid)[1]<Y:
                lo=mid
            else:
                hi=mid
        return X,(lo+hi)/2

    def material(self,X,Y):
        x,y=self.inverse(X,Y)
        s,d=self.jacobian(x,y)
        return 1/d,s/d,(s*s+d*d)/d,1/d


def matched_steps(ops,end,factor=1):
    if not isinstance(factor,int) or factor<1:
        raise ValueError("positive integer time refinement factor required")
    return max(2,math.ceil(end/(.65*min(op.dt_limit for op in ops))))*factor


def manufactured(n,amplitude=.035,end=.19):
    deformation=Deformation(amplitude)
    op=assemble(n,deformation.material)
    def mode(X,Y):
        x,y=deformation.inverse(X,Y)
        return math.sin(2*math.pi*x)*math.cos(2*math.pi*y)
    initial=sample(n,mode)
    exact=[u*math.cos(2*math.pi*math.sqrt(2)*end) for u in initial]
    steps=matched_steps([op],end)
    out,info=evolve_steps(op,initial,end,steps)
    error=math.sqrt(math.fsum((a-b)**2 for a,b in zip(out,exact))/
                    math.fsum(u*u for u in initial))
    return dict(n=n,amplitude=amplitude,amplitude_error=error,**info)


def comparison(n=48,amplitude=.035,end=.35,profile=0,factor=1):
    if profile not in (0,1):
        raise ValueError("unknown profile")
    deformation=Deformation(amplitude)
    def wrong(X,Y):
        return 1.,0.,1.,deformation.material(X,Y)[3]
    ops={name:assemble(n,material) for name,material in
         (("vacuum",vacuum),("transformed",deformation.material),("density_only",wrong))}
    cx,cy=(.12,.5) if profile==0 else (.12,.37)
    radius=.075
    def pulse(x,y):
        q=((x-cx)**2+(y-cy)**2)/radius**2
        return max(0.,1-q)**4
    initial=sample(n,pulse)
    # Exclude a buffer around the transformation's support.
    exterior=sample(n,lambda x,y:x<.22 or x>.78 or y<.22 or y>.78)
    steps=matched_steps(list(ops.values()),end,factor)
    fields={};diags={}
    for name,op in ops.items():
        fields[name],diags[name]=evolve_steps(op,initial,end,steps)
    return dict(n=n,amplitude=amplitude,end=end,profile=profile,factor=factor,steps=steps,
                exterior_nodes=sum(exterior),
                cases={name:dict(exterior_error=relative_l2(field,fields["vacuum"],exterior),
                                 **diags[name]) for name,field in fields.items()})
