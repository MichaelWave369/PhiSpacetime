"""R1b numerical verification and exploratory transformed scalar medium."""
import json
import math
from .medium import manufactured, experiment


def benchmark():
    results={}
    for label,variable in (("constant_anisotropy",False),("variable_density",True)):
        rows=[manufactured(n,variable) for n in (16,32,64)]
        results[label]=dict(rows=rows,orders=[
            math.log2(a["relative_l2"]/b["relative_l2"]) for a,b in zip(rows,rows[1:])])
    return dict(scope="scalar periodic finite elements; exploratory, not Maxwell or physical cloaking",
                verification=results,experiment=experiment())


if __name__=="__main__":
    print(json.dumps(benchmark(),indent=2,allow_nan=False))
