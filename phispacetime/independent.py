"""Independent radial distributions and angular phases; common physical probes."""
import bisect
import math
from .annular import mesh, assemble_annulus, harmonic, boundary_mask
from .interface import evolve_steps
from .wave import relative_l2


def independent_nodes(n,case):
    v,p,ann=mesh(n)
    if case=="virtual":
        return v,0.
    if case not in ("transformed","bare"):
        raise ValueError("unknown case")
    half=(n-1)//2
    # Independent physical placement: nonlinear inner spacing, different
    # exterior spacing and a rotated angular lattice. Not mapped virtual nodes.
    nodes=[p[0]+(ann.outer-p[0])*(i/half)**1.2 for i in range(half)]
    nodes += [ann.outer+(p[-1]-ann.outer)*(i/half)**1.1 for i in range(half+1)]
    return nodes,math.pi/(3*n)


def interpolate(field,radii,phase,R,theta):
    n=len(radii)
    if len(field)!=n*n or not radii[0]<=R<=radii[-1]:
        raise ValueError("probe outside grid")
    i=min(n-2,max(0,bisect.bisect_right(radii,R)-1))
    a=(R-radii[i])/(radii[i+1]-radii[i])
    t=((theta-phase)/(2*math.pi)*n)%n
    j=int(t)%n
    b=t-int(t)
    k=(j+1)%n
    return ((1-a)*((1-b)*field[i*n+j]+b*field[i*n+k])
            +a*((1-b)*field[(i+1)*n+j]+b*field[(i+1)*n+k]))


def probes(field,radii,phase):
    return [interpolate(field,radii,phase,R,2*math.pi*(j+.23)/37)
            for R in (.255,.285,.315,.345,.375,.405,.435) for j in range(37)]


def interface_flux(field,radii,phase,case,mode):
    """One-sided secant flux estimates versus exact static flux at R=.22.
    First-order reconstruction, not the weak discrete FEM flux.
    """
    _,_,ann=mesh(len(radii))
    n=len(radii)
    i=min(range(n),key=lambda k:abs(radii[k]-ann.outer))
    R=ann.outer
    e=radii[0] if case=="bare" else ann.virtual_inner
    L=radii[-1]
    derivative=mode*((R/L)**mode+(e*e/(R*L))**mode)/(R*(1-(e/L)**(2*mode)))
    ar=ann.slope if case=="transformed" else 1.
    left=[];right=[];exact=[]
    for j in range(n):
        left.append(ar*(field[i*n+j]-field[(i-1)*n+j])/(R-radii[i-1]))
        right.append((field[(i+1)*n+j]-field[i*n+j])/(radii[i+1]-R))
        exact.append(derivative*math.cos(mode*(2*math.pi*j/n+phase)))
    return dict(left_error=relative_l2(left,exact),right_error=relative_l2(right,exact),
                jump=relative_l2(left,right))


def compare(n=33,mode=1,end=.24,time_factor=1):
    if mode not in (1,2) or time_factor not in (1,2):
        raise ValueError("invalid comparison parameters")
    grids={case:independent_nodes(n,case) for case in ("virtual","transformed","bare")}
    ops={case:assemble_annulus(n,case,radii)[0] for case,(radii,_) in grids.items()}
    steps=max(2,math.ceil(end/(.65*min(op.dt_limit for op in ops.values()))))*time_factor
    waves={};statics={};diagnostics={}
    for case,(radii,phase) in grids.items():
        initial=[max(0.,1-((r-.34)/.065)**2)**4*math.cos(mode*(2*math.pi*j/n+phase))
                 for r in radii for j in range(n)]
        out,info=evolve_steps(ops[case],initial,end,steps,boundary_mask(n))
        steady,static=harmonic(n,case,radii,phase,mode)
        waves[case]=probes(out,radii,phase)
        statics[case]=probes(steady,radii,phase)
        diagnostics[case]=dict(wave=info,static=static,
                               flux=interface_flux(steady,radii,phase,case,mode))
    return dict(n=n,mode=mode,end=end,steps=steps,time_factor=time_factor,probe_count=len(waves["virtual"]),
                cases={case:dict(wave_error=relative_l2(waves[case],waves["virtual"]),
                                 static_error=relative_l2(statics[case],statics["virtual"]),
                                 **diagnostics[case]) for case in grids})
