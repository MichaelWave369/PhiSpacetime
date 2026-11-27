"""R1c: interface and matched-observation verification for scalar waves."""
import math
from .medium import assemble, evolve, radial_medium, sample, vacuum
from .wave import relative_l2


def mapped_x(x):
    """Continuous periodic-circle coordinate map with two linear layers."""
    return .5*x if x < .5 else .25+1.5*(x-.5)


def layered_material(x,y):
    s=.5 if x < .5 else 1.5
    return 1/s,0.,s,s


def layered_mode(n,end=.19):
    """Exact solution sin(2*pi*X(x))*cos(2*pi*t).
    Axx*u_x is continuous at both interfaces, including periodic seam.
    """
    if n%2:
        raise ValueError("interface must align with mesh")
    op=assemble(n,layered_material)
    initial=sample(n,lambda x,y:math.sin(2*math.pi*mapped_x(x)))
    exact=[u*math.cos(2*math.pi*end) for u in initial]
    out,info=evolve(op,initial,[0.]*(n*n),end)
    return dict(n=n,amplitude_error=math.sqrt(
        math.fsum((a-b)**2 for a,b in zip(out,exact))/math.fsum(u*u for u in initial)),**info)


def evolve_steps(op,initial,end,steps,fixed=None):
    """Choose the existing integrator's safety factor to enforce exact step count."""
    if not isinstance(steps,int) or steps<2:
        raise ValueError("at least two steps required")
    safety=end/((steps-.25)*op.dt_limit)
    if not 0<safety<1:
        raise ValueError("requested time step violates conservative bound")
    out,info=evolve(op,initial,[0.]*(op.n**2),end,fixed,safety=safety)
    if info["steps"]!=steps:
        raise RuntimeError("time-step matching failed")
    return out,info


def temporal_mode(n=16,end=.19):
    """Isolate time error against the exact semidiscrete periodic eigenmode."""
    op=assemble(n)
    initial=sample(n,lambda x,y:math.sin(2*math.pi*x))
    k=op.apply(initial)
    lam=math.fsum(u*v for u,v in zip(initial,k))/math.fsum(
        m*u*u for m,u in zip(op.mass,initial))
    residual=max(abs(v-lam*m*u) for m,u,v in zip(op.mass,initial,k))
    exact=[u*math.cos(math.sqrt(lam)*end) for u in initial]
    minimum=max(2,math.ceil(end/(.65*op.dt_limit)))
    rows=[]
    for steps in (minimum,2*minimum,4*minimum):
        out,info=evolve_steps(op,initial,end,steps)
        rows.append(dict(steps=steps,error=math.sqrt(
            math.fsum((a-b)**2 for a,b in zip(out,exact))/math.fsum(u*u for u in initial)),
            relative_energy_drift=info["relative_energy_drift"]))
    return dict(n=n,eigenvalue=lam,eigen_residual=residual,rows=rows)


def matched_experiment(n=48,end=.38,profile=0,time_factor=1):
    """Compare all cases on exactly matched time steps and exterior nodes.

    The compact pulse is entirely outside the transformed material. A finite
    virtual disk is the correct reference target, not perfect vacuum.
    """
    if profile not in (0,1) or time_factor not in (1,2,4):
        raise ValueError("invalid experiment variant")
    inner,outer,virtual=.12,.22,.04
    cx,cy=(.14,.5) if profile==0 else (.15,.61)
    radius=.07
    if math.hypot(cx-.5,cy-.5)-radius <= outer:
        raise ValueError("initial pulse overlaps transformed material")
    def pulse(x,y):
        q=((x-cx)**2+(y-cy)**2)/radius**2
        return max(0.,1-q)**4
    initial=sample(n,pulse)
    radii=sample(n,lambda x,y:math.hypot(x-.5,y-.5))
    fixed=[v<=inner for v in radii]
    small=[v<=virtual for v in radii]
    exterior=[v>=outer+.03 for v in radii]
    base=assemble(n)
    transform=radial_medium(inner,outer,virtual)
    def wrong(x,y):
        return 1.,0.,1.,transform(x,y)[3]
    operators={"vacuum":(base,None),"bare":(base,fixed),
               "virtual":(base,small),"transformed":(assemble(n,transform),fixed),
               "density_only":(assemble(n,wrong),fixed)}
    steps=max(2,math.ceil(end/(.65*min(op.dt_limit for op,_ in operators.values()))))*time_factor
    fields={}
    diagnostics={}
    for name,(op,mask) in operators.items():
        fields[name],diagnostics[name]=evolve_steps(op,initial,end,steps,mask)
    cases={}
    for name,field in fields.items():
        cases[name]=dict(vacuum_error=relative_l2(field,fields["vacuum"],exterior),
                         virtual_error=relative_l2(field,fields["virtual"],exterior),
                         **diagnostics[name])
    return dict(n=n,end=end,profile=profile,time_factor=time_factor,steps=steps,
                initial_support_gap=math.hypot(cx-.5,cy-.5)-radius-outer,
                exterior_nodes=sum(exterior),physical_nodes=sum(fixed),
                virtual_nodes=sum(small),cases=cases)
