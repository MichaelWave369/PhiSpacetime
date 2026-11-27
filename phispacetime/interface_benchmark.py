"""R1c deterministic report; preserves favorable and unfavorable outcomes."""
import json
import math
from .interface import layered_mode, temporal_mode, matched_experiment


def orders(values):
    return [math.log2(a/b) for a,b in zip(values,values[1:])]


def benchmark():
    layers=[layered_mode(n) for n in (16,32,64)]
    time=temporal_mode()
    time["orders"]=orders([r["error"] for r in time["rows"]])
    # Paired spatial and temporal refinement, plus a second illumination.
    comparisons=[matched_experiment(n,end,profile,factor) for n,end,profile,factor in
                 ((32,.38,0,1),(48,.38,0,1),(64,.38,0,1),
                  (48,.38,0,2),(48,.46,0,1),(48,.38,1,1))]
    return dict(scope="periodic scalar FEM; exploratory curved interfaces, not physical cloaking",
                layered=dict(rows=layers,orders=orders([r["amplitude_error"] for r in layers])),
                temporal=time,comparisons=comparisons)


if __name__=="__main__":
    print(json.dumps(benchmark(),indent=2,allow_nan=False))
