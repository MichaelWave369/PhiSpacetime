import json
import math
from .nonradial import manufactured, comparison


def benchmark():
    rows=[manufactured(n) for n in (16,32,64)]
    return dict(scope="periodic scalar wave; compact nonradial deformation; no obstacle or cloak",
                manufactured=rows,orders=[math.log2(a["amplitude_error"]/b["amplitude_error"])
                                           for a,b in zip(rows,rows[1:])],
                comparisons=[comparison(n,a,t,p,f) for n,a,t,p,f in
                             ((24,.035,.35,0,1),(48,.035,.35,0,1),(72,.035,.35,0,1),
                              (48,-.025,.35,1,1),(48,.035,.35,0,2),(48,.035,.45,0,1))])


if __name__=="__main__":
    print(json.dumps(benchmark(),indent=2,allow_nan=False))
