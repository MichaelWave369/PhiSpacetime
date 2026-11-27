"""Fitted polar Q1 scalar FEM on a bounded annulus.

Circular interfaces are exact in polar coordinates. Virtual and transformed
meshes are deliberately paired by the prescribed map; their agreement is
a discrete covariance check, not independent evidence of cloaking.
"""
import math
from .optics import Annulus
from .medium import Operator, evolve
from .interface import evolve_steps
from .wave import relative_l2


def mesh(n, inner=.12, interface=.22, outer=.45, epsilon=.04):
    if not isinstance(n,int) or n<9 or n%2!=1:
        raise ValueError("odd n >= 9 required")
    if not 0<epsilon<inner<interface<outer:
        raise ValueError("invalid annulus radii")
    half=(n-1)//2
    virtual=[epsilon+(interface-epsilon)*i/half for i in range(half)]
    virtual += [interface+(outer-interface)*i/half for i in range(half+1)]
    ann=Annulus(inner,interface,epsilon)
    physical=[ann.radius(r) for r in virtual]
    return virtual,physical,ann


def assemble_annulus(n,case,radial_nodes=None):
    virtual,physical,ann=mesh(n)
    if case not in ("virtual","transformed","bare"):
        raise ValueError("unknown case")
    radii=virtual if case=="virtual" else physical
    if radial_nodes is not None:
        radii=list(radial_nodes)
        expected=virtual if case=="virtual" else physical
        if (len(radii)!=n or not all(math.isfinite(v) for v in radii)
                or any(a>=b for a,b in zip(radii,radii[1:]))
                or abs(radii[0]-expected[0])>1e-14 or abs(radii[-1]-expected[-1])>1e-14
                or not any(abs(v-ann.outer)<1e-14 for v in radii)):
            raise ValueError("invalid fitted radial nodes")
    mass=[0.]*(n*n)
    rows=[{} for _ in mass]
    dtheta=2*math.pi/n
    g=1/math.sqrt(3)
    signs=((-1,-1),(1,-1),(1,1),(-1,1))
    for i in range(n-1):
        lo,hi=radii[i:i+2]
        h=hi-lo
        for j in range(n):
            jp=(j+1)%n
            ids=[i*n+j,(i+1)*n+j,(i+1)*n+jp,i*n+jp]
            for xi in (-g,g):
                R=lo+h*(xi+1)/2
                ar,at,rho=1.,1.,1.
                if case=="transformed" and R<ann.outer:
                    rv=ann.inverse(R)
                    a,b,_=ann.stretches(rv)
                    ar,at,rho=a/b,b/a,1/(a*b)
                for eta in (-g,g):
                    N=[(1+s*xi)*(1+t*eta)/4 for s,t in signs]
                    dR=[s*(1+t*eta)/(2*h) for s,t in signs]
                    dT=[t*(1+s*xi)/(2*dtheta) for s,t in signs]
                    w=h*dtheta/4
                    for p,k in enumerate(ids):
                        mass[k]+=w*rho*R*N[p]
                        for q,l in enumerate(ids):
                            value=w*(ar*R*dR[p]*dR[q]+at/R*dT[p]*dT[q])
                            rows[k][l]=rows[k].get(l,0.)+value
    return Operator(n,mass,[sorted(row.items()) for row in rows]),radii


def boundary_mask(n):
    return [i==0 or i==n-1 for i in range(n) for _ in range(n)]


def harmonic(n,case,radial_nodes=None,phase=0.,mode=1):
    """Solve the static Dirichlet problem by CG; exact virtual m=1 mode."""
    op,radii=assemble_annulus(n,case,radial_nodes)
    if not isinstance(mode,int) or mode<1 or mode>=n/2 or not math.isfinite(phase):
        raise ValueError("invalid harmonic mode or phase")
    mask=boundary_mask(n)
    boundary=[math.cos(mode*(2*math.pi*j/n+phase)) if i==n-1 else 0.
              for i in range(n) for j in range(n)]
    rhs=[-v if not m else 0. for v,m in zip(op.apply(boundary),mask)]
    def apply(x):
        return [0. if m else v for v,m in zip(op.apply(x),mask)]
    def dot(a,b):
        return math.fsum(x*y for x,y in zip(a,b))
    x=[0.]*(n*n)
    residual=rhs[:]
    p=residual[:]
    rr=dot(residual,residual)
    target=rr*1e-22
    count=0
    while rr>target and count<20*n*n:
        ap=apply(p)
        denom=dot(p,ap)
        if denom<=0:
            raise RuntimeError("nonpositive CG curvature")
        alpha=rr/denom
        x=[a+alpha*b for a,b in zip(x,p)]
        residual=[a-alpha*b for a,b in zip(residual,ap)]
        new=dot(residual,residual)
        p=[a+(new/rr)*b for a,b in zip(residual,p)]
        rr=new
        count+=1
    if rr>target:
        raise RuntimeError("CG failed to converge")
    field=[a+b for a,b in zip(x,boundary)]
    _,_,ann=mesh(n)
    reference_radii=[ann.inverse(r) for r in radii] if case=="transformed" else radii
    e=reference_radii[0]
    L=reference_radii[-1]
    exact=[((r/L)**mode-(e*e/(r*L))**mode)/(1-(e/L)**(2*mode))*math.cos(mode*(2*math.pi*j/n+phase))
           for r in reference_radii for j in range(n)]
    residual_actual=math.sqrt(dot([a-b for a,b in zip(rhs,apply(x))],
                                  [a-b for a,b in zip(rhs,apply(x))]))
    return field,dict(relative_l2=relative_l2(field,exact),iterations=count,
                      relative_residual=residual_actual/math.sqrt(dot(rhs,rhs)))


def comparison(n=33,end=.3):
    vr,pr,ann=mesh(n)
    cases={c:assemble_annulus(n,c) for c in ("virtual","transformed","bare")}
    # Shared exterior initial data; radial compact support excludes all
    # transformed material and both fixed boundaries.
    initial=[max(0.,1-((r-.34)/.065)**2)**4*math.cos(2*math.pi*j/n)
             for r in vr for j in range(n)]
    mask=boundary_mask(n)
    exterior=[r>=ann.outer+.03 for r in vr for _ in range(n)]
    steps=max(2,math.ceil(end/(.65*min(op.dt_limit for op,_ in cases.values()))))
    wave={}
    static={}
    for name,(op,_) in cases.items():
        wave[name]=evolve_steps(op,initial,end,steps,mask)
        static[name]=harmonic(n,name)
    ref=cases["virtual"][0]
    trans=cases["transformed"][0]
    stiffness_error=max(abs(a-dict(trans.rows[i])[j]) for i,row in enumerate(ref.rows) for j,a in row)
    mass_error=max(abs(a-b) for a,b in zip(ref.mass,trans.mass))
    return dict(n=n,end=end,steps=steps,exterior_nodes=sum(exterior),
                max_matrix_difference=stiffness_error,max_mass_difference=mass_error,
                cases={name:dict(wave_exterior_error=relative_l2(wave[name][0],wave["virtual"][0],exterior),
                                 static_exterior_error=relative_l2(static[name][0],static["virtual"][0],exterior),
                                 wave=wave[name][1],static=static[name][1])
                       for name in cases})
