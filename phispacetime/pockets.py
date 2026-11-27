"""R2a: conforming P1 finite elements on hierarchical square leaves.

Integer lattice topology, shared edge breakpoints and center-fan triangles.
Static forced scalar PDE -laplacian(u)=f, homogeneous Dirichlet boundary.
No upstream pocket API or time-dependent adaptive solver is claimed.
"""
from dataclasses import dataclass
import math
import heapq
from collections import defaultdict


@dataclass(frozen=True)
class Cell:
    x: int
    y: int
    h: int
    depth: int


def bump(x,c=.36,w=.19):
    z=(x-c)/w
    if abs(z)>=1:
        return 0.,0.
    t=1-z*z
    return t**4,(-8*t**3+48*z*z*t*t)/(w*w)


def exact(x,y,localized=True):
    return math.sin(math.pi*x)*math.sin(math.pi*y)+(
        .5*bump(x)[0]*bump(y,.61)[0] if localized else 0.)


def source(x,y,localized=True):
    f=2*math.pi**2*math.sin(math.pi*x)*math.sin(math.pi*y)
    if localized:
        a,aa=bump(x)
        b,bb=bump(y,.61)
        f-=.5*(aa*b+a*bb)
    return f


def leaves(arity,budget,adaptive=True,localized=True,max_depth=5):
    if arity not in (2,3) or not isinstance(budget,int) or budget<1 or not 1<=max_depth<=7:
        raise ValueError("invalid hierarchy configuration")
    scale=arity**max_depth
    def priority(c):
        if not adaptive:
            return -c.depth
        h=c.h/scale
        return h*h*max(abs(source((c.x+c.h*i/4)/scale,(c.y+c.h*j/4)/scale,localized))
                       for i in (1,2,3) for j in (1,2,3))
    heap=[]
    serial=0
    def push(c):
        nonlocal serial
        heapq.heappush(heap,(-priority(c),serial,c))
        serial+=1
    push(Cell(0,0,scale,0))
    frozen=[]
    while heap and len(heap)+len(frozen)+arity**2-1<=budget:
        _,_,c=heapq.heappop(heap)
        if c.depth==max_depth:
            frozen.append(c)
            continue
        h=c.h//arity
        for j in range(arity):
            for i in range(arity):
                push(Cell(c.x+i*h,c.y+j*h,h,c.depth+1))
    return scale,sorted(frozen+[item[2] for item in heap],key=lambda c:(c.y,c.x,c.depth))


def triangulate(scale,cells):
    """Conforming mesh even at hanging vertices; no interpolated constraints."""
    horizontal=defaultdict(set)
    vertical=defaultdict(set)
    for c in cells:
        x,y,h=2*c.x,2*c.y,2*c.h
        for a,b in ((x,y),(x+h,y),(x+h,y+h),(x,y+h)):
            horizontal[b].add(a)
            vertical[a].add(b)
    points=[];ids={};triangles=[]
    def node(p):
        if p not in ids:
            ids[p]=len(points)
            points.append(p)
        return ids[p]
    for c in cells:
        x,y,h=2*c.x,2*c.y,2*c.h
        bottom=[(a,y) for a in sorted(horizontal[y]) if x<=a<x+h]
        right=[(x+h,b) for b in sorted(vertical[x+h]) if y<=b<y+h]
        top=[(a,y+h) for a in sorted(horizontal[y+h],reverse=True) if x<a<=x+h]
        left=[(x,b) for b in sorted(vertical[x],reverse=True) if y<b<=y+h]
        ring=bottom+right+top+left
        center=node((x+h//2,y+h//2))
        for a,b in zip(ring,ring[1:]+ring[:1]):
            triangles.append((center,node(a),node(b)))
    xy=[(x/(2*scale),y/(2*scale)) for x,y in points]
    fixed=[x in (0,2*scale) or y in (0,2*scale) for x,y in points]
    return xy,triangles,fixed


QUAD=((2/3,1/6,1/6),(1/6,2/3,1/6),(1/6,1/6,2/3))


def geometry(xy,tri):
    (x,y),(a,b),(c,d)=[xy[i] for i in tri]
    det=(a-x)*(d-y)-(c-x)*(b-y)
    if det<=0:
        raise ValueError("nonpositive triangle")
    return det/2,((b-d,c-a),(d-y,x-c),(y-b,a-x)),det


def assemble(xy,triangles,localized=True):
    rows=[{} for _ in xy]
    mass=[0.]*len(xy);load=[0.]*len(xy)
    for tri in triangles:
        area,grad,det=geometry(xy,tri)
        for p,i in enumerate(tri):
            mass[i]+=area/3
            for q,j in enumerate(tri):
                value=area*(grad[p][0]*grad[q][0]+grad[p][1]*grad[q][1])/det**2
                rows[i][j]=rows[i].get(j,0.)+value
        for w in QUAD:
            x=sum(w[k]*xy[i][0] for k,i in enumerate(tri))
            y=sum(w[k]*xy[i][1] for k,i in enumerate(tri))
            f=source(x,y,localized)
            for k,i in enumerate(tri):
                load[i]+=area/3*w[k]*f
    return [sorted(row.items()) for row in rows],mass,load


def solve(rows,load,fixed):
    """Jacobi-preconditioned conjugate gradient on unconstrained nodes."""
    def dot(a,b):
        return math.fsum(x*y for x,y in zip(a,b))
    def apply(u):
        return [0. if m else math.fsum(v*u[j] for j,v in row)
                for row,m in zip(rows,fixed)]
    rhs=[0. if m else v for v,m in zip(load,fixed)]
    diag=[dict(row)[i] for i,row in enumerate(rows)]
    x=[0.]*len(rows);r=rhs[:]
    norm=math.sqrt(dot(r,r))
    if not norm:
        return x,0,0.
    z=[v/d for v,d in zip(r,diag)]
    p=z[:];rz=dot(r,z)
    iterations=0
    while math.sqrt(dot(r,r))>1e-10*norm and iterations<10*len(rows):
        ap=apply(p)
        denom=dot(p,ap)
        if denom<=0:
            raise RuntimeError("nonpositive stiffness curvature")
        a=rz/denom
        x=[u+a*v for u,v in zip(x,p)]
        r=[u-a*v for u,v in zip(r,ap)]
        z=[v/d for v,d in zip(r,diag)]
        new=dot(r,z)
        p=[u+new/rz*v for u,v in zip(z,p)]
        rz=new;iterations+=1
    actual=math.sqrt(dot([a-b for a,b in zip(rhs,apply(x))],
                         [a-b for a,b in zip(rhs,apply(x))]))/norm
    if actual>1e-9:
        raise RuntimeError("linear solve failed")
    return x,iterations,actual


def errors(xy,triangles,field,localized=True):
    # Independent four-point-per-axis Duffy quadrature for the error norm.
    nodes=(.0694318442029737,.330009478207572,.669990521792428,.930568155797026)
    weights=(.173927422568727,.326072577431273,.326072577431273,.173927422568727)
    numerator=denominator=0.
    for tri in triangles:
        area,_,_=geometry(xy,tri)
        for a,wa in zip(nodes,weights):
            for b,wb in zip(nodes,weights):
                w=(1-a,a*(1-b),a*b)
                x=sum(w[k]*xy[i][0] for k,i in enumerate(tri))
                y=sum(w[k]*xy[i][1] for k,i in enumerate(tri))
                truth=exact(x,y,localized)
                u=sum(w[k]*field[i] for k,i in enumerate(tri))
                q=2*area*a*wa*wb
                numerator+=q*(u-truth)**2
                denominator+=q*truth**2
    return math.sqrt(numerator/denominator)


def run(arity,budget,adaptive=True,localized=True):
    scale,cells=leaves(arity,budget,adaptive,localized)
    xy,tri,fixed=triangulate(scale,cells)
    rows,mass,load=assemble(xy,tri,localized)
    field,it,residual=solve(rows,load,fixed)
    return dict(arity=arity,adaptive=adaptive,localized=localized,budget=budget,
                leaves=len(cells),nodes=len(xy),free_dofs=len(xy)-sum(fixed),
                triangles=len(tri),matrix_nonzeros=sum(len(row) for row in rows),
                iterations=it,relative_residual=residual,relative_l2=errors(xy,tri,field,localized),
                area=math.fsum(mass),max_depth=max(c.depth for c in cells))
