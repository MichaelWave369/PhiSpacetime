"""Conservative 2-D scalar wave equation on a periodic unit square.

u_tt = Laplacian(u), with nondimensional wave speed 1. Optional fixed
zero-valued grid nodes represent a scalar Dirichlet obstacle, not a
general electromagnetic conductor or a transformed material.
"""
from dataclasses import dataclass
import math


@dataclass(frozen=True)
class Grid:
    n: int
    dt: float

    def __post_init__(self):
        if (not isinstance(self.n, int) or self.n < 4
                or not math.isfinite(self.dt) or self.dt <= 0
                or self.dt * self.n > 1 / math.sqrt(2)):
            raise ValueError("invalid grid or 2-D CFL limit exceeded")

    @property
    def h(self):
        return 1 / self.n

    def sample(self, f):
        return [f(i*self.h, j*self.h)
                for j in range(self.n) for i in range(self.n)]

    def neighbors(self):
        n = self.n
        return [((i-1)%n+j*n, (i+1)%n+j*n,
                 i+((j-1)%n)*n, i+((j+1)%n)*n)
                for j in range(n) for i in range(n)]


def solve(grid, initial, velocity, steps, fixed=None):
    """Return final field and discrete energy diagnostics.

    Leapfrog uses centered differences; Taylor initialization is
    second-order. Obstacles remain fixed at zero at all time levels.
    """
    n, dt = grid.n, grid.dt
    count = n*n
    if not isinstance(steps, int) or steps < 1:
        raise ValueError("positive integer step count required")
    if len(initial) != count or len(velocity) != count:
        raise ValueError("field size mismatch")
    if not all(math.isfinite(v) for v in (*initial, *velocity)):
        raise ValueError("nonfinite initial data")
    mask = [False]*count if fixed is None else list(fixed)
    if len(mask) != count or any(type(v) is not bool for v in mask):
        raise ValueError("invalid fixed-node mask")
    if any(mask[k] and (initial[k] != 0 or velocity[k] != 0)
           for k in range(count)):
        raise ValueError("initial data incompatible with fixed nodes")
    nbr = grid.neighbors()
    lam = (dt*n)**2
    prev = list(initial)

    def lap(field, k):
        a, b, c, d = nbr[k]
        return field[a]+field[b]+field[c]+field[d]-4*field[k]

    current = [0.0 if mask[k] else prev[k]+dt*velocity[k]+0.5*lam*lap(prev,k)
               for k in range(count)]

    def energy(a, b):
        # Conserved leapfrog energy at a half step. Cross-time potential
        # term, not the instantaneous continuum energy.
        kinetic = math.fsum(((b[k]-a[k])/dt)**2 for k in range(count))*grid.h**2
        potential = math.fsum(
            (b[nbr[k][1]]-b[k])*(a[nbr[k][1]]-a[k])+
            (b[nbr[k][3]]-b[k])*(a[nbr[k][3]]-a[k])
            for k in range(count))
        return 0.5*(kinetic+potential)

    e0 = energy(prev,current)
    max_drift = 0.0
    for _ in range(1,steps):
        nxt = [0.0 if mask[k] else 2*current[k]-prev[k]+lam*lap(current,k)
               for k in range(count)]
        prev,current = current,nxt
        e = energy(prev,current)
        if not math.isfinite(e):
            raise RuntimeError("nonfinite energy")
        max_drift = max(max_drift, abs(e-e0))
    return current, dict(initial_energy=e0, final_energy=energy(prev,current),
                         max_absolute_energy_drift=max_drift,
                         max_relative_energy_drift=(max_drift/abs(e0) if e0 else None))


def relative_l2(actual, reference, include=None):
    if len(actual) != len(reference) or (include is not None and len(include) != len(actual)):
        raise ValueError("field size mismatch")
    indices = [k for k in range(len(actual)) if include is None or include[k]]
    denom = math.fsum(reference[k]**2 for k in indices)
    if denom <= 0:
        raise ValueError("zero reference norm")
    return math.sqrt(math.fsum((actual[k]-reference[k])**2 for k in indices)/denom)


def standing_mode(n, end=0.25):
    """Exact periodic standing mode with zero initial velocity."""
    steps = math.ceil(end*n/0.5)
    grid = Grid(n,end/steps)
    initial = grid.sample(lambda x,y: math.sin(2*math.pi*x)*math.cos(2*math.pi*y))
    field, energy = solve(grid,initial,[0.0]*(n*n),steps)
    omega = 2*math.pi*math.sqrt(2)
    exact = [u*math.cos(omega*end) for u in initial]
    return dict(n=n,steps=steps,dt=grid.dt,relative_l2=relative_l2(field,exact),**energy)


def obstacle_control(n=64,end=0.35):
    """Compare vacuum with fixed-zero disk, using identical initial data."""
    steps = math.ceil(end*n/0.5)
    grid = Grid(n,end/steps)
    # Compact smooth-enough pulse initially disjoint from the disk.
    def pulse(x,y):
        q = ((x-0.25)**2+(y-0.5)**2)/0.12**2
        return max(0.0,1-q)**4
    radii = grid.sample(lambda x,y: (x-0.65)**2+(y-0.5)**2)
    fixed = [r <= 0.12**2 for r in radii]
    exterior = [r >= 0.18**2 for r in radii]
    initial = grid.sample(pulse)
    velocity = [0.0]*(n*n)
    vacuum,ev = solve(grid,initial,velocity,steps)
    obstacle,eo = solve(grid,initial,velocity,steps,fixed)
    return dict(n=n,steps=steps,dt=grid.dt,obstacle_nodes=sum(fixed),
                exterior_nodes=sum(exterior),exterior_relative_l2=relative_l2(obstacle,vacuum,exterior),
                obstacle_max_abs=max(abs(v) for v,m in zip(obstacle,fixed) if m),
                vacuum_energy=ev,obstacle_energy=eo)
