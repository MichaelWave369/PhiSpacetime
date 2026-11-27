import json
import math
from .annular import comparison


def benchmark():
    rows=[comparison(n) for n in (17,33,65)]
    errors=[r["cases"]["virtual"]["static"]["relative_l2"] for r in rows]
    return dict(scope="bounded fitted polar scalar FEM; paired-map covariance, not independent cloak verification",
                rows=rows,static_orders=[math.log2(a/b) for a,b in zip(errors,errors[1:])])


if __name__=="__main__":
    print(json.dumps(benchmark(),indent=2,allow_nan=False))
