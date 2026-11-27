"""R1a reproducible benchmark; prints JSON."""
import json
import math
from .wave import standing_mode, obstacle_control


def benchmark():
    modes = [standing_mode(n) for n in (24,48,96)]
    orders = [math.log2(a["relative_l2"]/b["relative_l2"])
              for a,b in zip(modes,modes[1:])]
    return dict(scope="periodic 2-D scalar wave equation; no cloak or vector Maxwell solve",
                standing_modes=modes, observed_orders=orders,
                obstacle=obstacle_control())


if __name__ == "__main__":
    print(json.dumps(benchmark(),indent=2,allow_nan=False))
