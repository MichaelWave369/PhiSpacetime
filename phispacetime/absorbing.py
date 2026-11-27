"""Nonperiodic scalar FEM with a dissipative sponge and fixed outer wall.

A finite sponge is an approximate absorbing treatment, not a PML or
an exact radiation condition. Energy accounting is discrete.
"""
import math
from .medium import assemble, vacuum
from .nonradial import Deformation
from .wave import relative_l2


def coordinates(n,length=1.,origin=(0.,0.)):
    return [(origin[0]+length*i/(n-1),origin[1]+length*j/(n-1))
            for j in range(n) for i in range(n)]


def wall_mask(n):
    return [i in (0,n-1) or j in (0,n-1)
            for j in range(n) for i in range(n)]


def sponge(x,y,width=.18,strength=35.):
    if not 0<width<.5 or not math.isfinite(strength) or strength<0:
        raise ValueError("invalid sponge parameters")
    d=min(x,y,1-x,1-y)
    return strength*max(0.,1-d/width)**3


def evolve_damped(op,initial,velocity,end,steps,damping,fixed):
    count=op.n**2
    if (not isinstance(steps,int) or steps<2 or not math.isfinite(end) or end<=0
            or any(len(v)!=count for v in (initial,velocity,damping,fixed))
            or not all(math.isfinite(v) for v in (*initial,*velocity,*damping))
            or any(v<0 for v in damping) or any(type(m) is not bool for m in fixed)):
        raise ValueError("invalid damped evolution")
    dt=end/steps
    if dt>=op.dt_limit:
        raise ValueError("unsafe step")
    if any(m and (u!=0 or v!=0) for m,u,v in zip(fixed,initial,velocity)):
        raise ValueError("incompatible wall data")
    prev=list(initial)
    k=op.apply(prev)
    curr=[0. if fixed[i] else u+dt*v-.5*dt*dt*(k[i]/op.mass[i]+damping[i]*v)
          for i,(u,v) in enumerate(zip(initial,velocity))]

    def energy(a,b):
        return .5*math.fsum(m*((v-u)/dt)**2+u*z
                            for m,u,v,z in zip(op.mass,a,b,op.apply(b)))

    e0=energy(prev,curr)
    last=e0
    dissipated=0.
    max_balance=0.
    max_increase=0.
    for _ in range(1,steps):
        k=op.apply(curr)
        nxt=[0. if fixed[i] else
             (2*curr[i]-(1-.5*dt*damping[i])*prev[i]-dt*dt*k[i]/op.mass[i])/
             (1+.5*dt*damping[i]) for i in range(count)]
        dissipated+=dt*math.fsum(m*s*((z-a)/(2*dt))**2
                                 for m,s,z,a in zip(op.mass,damping,nxt,prev))
        prev,curr=curr,nxt
        e=energy(prev,curr)
        if not math.isfinite(e):
            raise RuntimeError("nonfinite energy")
        max_increase=max(max_increase,e-last)
        max_balance=max(max_balance,abs(e+dissipated-e0))
        last=e
    scale=abs(e0) or 1.
    return curr,dict(steps=steps,dt=dt,initial_energy=e0,final_energy=last,
                     dissipated=dissipated,relative_balance_error=max_balance/scale,
                     max_relative_energy_increase=max_increase/scale)


def pulse(x,y):
    q=((x-.5)**2+(y-.5)**2)/.11**2
    return max(0.,1-q)**4


def matched_end(ops,end,factor=1):
    return max(2,math.ceil(end/(.6*min(op.dt_limit for op in ops))))*factor


def reflection(n=25,end=.95,width=.18,strength=35.,factor=1):
    """Small box versus larger undamped box, same h and dt.

    Larger box spans [-.5,1.5]^2. Its walls cannot affect the central
    observation square by the requested default time in the continuum.
    Numerical dispersion and finite-domain errors remain.
    """
    if n<9 or n%2!=1 or end<=0 or end>1.1:
        raise ValueError("odd n >=9 and 0<end<=1.1 required")
    big=2*(n-1)+1
    xy=coordinates(n)
    refxy=coordinates(big,2.,(-.5,-.5))
    op=assemble(n,periodic=False)
    refop=assemble(big,periodic=False,length=2.,origin=(-.5,-.5))
    steps=matched_end([op,refop],end,factor)
    init=[pulse(x,y) for x,y in xy]
    refinit=[pulse(x,y) for x,y in refxy]
    ref,rd=evolve_damped(refop,refinit,[0.]*(big*big),end,steps,
                          [0.]*(big*big),wall_mask(big))
    ids=[k for k,(x,y) in enumerate(xy) if .3<=x<=.7 and .3<=y<=.7]
    shift=(n-1)//2
    refprobe=[ref[(k//n+shift)*big+(k%n+shift)] for k in ids]
    # Normalize by initial pulse norm on the same probes, not tiny late-time reference.
    norm=math.sqrt(math.fsum(init[k]**2 for k in ids))
    cases={}
    for name,gamma in (("reflecting",[0.]*(n*n)),
                       ("sponge",[sponge(x,y,width,strength) for x,y in xy])):
        out,diag=evolve_damped(op,init,[0.]*(n*n),end,steps,gamma,wall_mask(n))
        error=math.sqrt(math.fsum((out[k]-v)**2 for k,v in zip(ids,refprobe)))/norm
        cases[name]=dict(return_error=error,**diag)
    return dict(n=n,reference_n=big,end=end,width=width,strength=strength,
                probes=len(ids),reference=rd,cases=cases)


def transformed_comparison(n=41,end=.65,factor=1):
    """Matched open-boundary approximation for compact nonradial deformation."""
    d=Deformation()
    xy=coordinates(n)
    def wrong(x,y):
        return 1.,0.,1.,d.material(x,y)[3]
    ops={name:assemble(n,mat,periodic=False) for name,mat in
         (("vacuum",vacuum),("transformed",d.material),("density_only",wrong))}
    def initial(x,y):
        q=((x-.205)**2+(y-.5)**2)/.04**2
        return max(0.,1-q)**4
    # Support lies outside the deformation but partly in the sponge;
    # the same damping applies in all configurations and is outside map support.
    u=[initial(x,y) for x,y in xy]
    gamma=[sponge(x,y) for x,y in xy]
    steps=matched_end(list(ops.values()),end,factor)
    probes=[.19<=x<=.81 and .19<=y<=.81 and
            (x<.24 or x>.76 or y<.24 or y>.76) for x,y in xy]
    fields={};diags={}
    for name,op in ops.items():
        fields[name],diags[name]=evolve_damped(op,u,[0.]*(n*n),end,steps,gamma,wall_mask(n))
    return dict(n=n,end=end,steps=steps,probe_count=sum(probes),
                cases={name:dict(exterior_error=relative_l2(v,fields["vacuum"],probes),**diags[name])
                       for name,v in fields.items()})
