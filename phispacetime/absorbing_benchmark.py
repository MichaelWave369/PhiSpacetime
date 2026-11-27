import json
from .absorbing import reflection, transformed_comparison


def benchmark():
    return dict(scope="finite scalar sponge and larger-box reference; not exact open space or Maxwell",
                reflection=[reflection(n,width=w,strength=s,factor=f) for n,w,s,f in
                            ((17,.18,35.,1),(25,.18,35.,1),(33,.18,35.,1),
                             (25,.12,35.,1),(25,.18,15.,1),(25,.18,35.,2))],
                transformed=[transformed_comparison(n,end,f) for n,end,f in
                             ((25,.55,1),(41,.55,1),(41,.65,1),(41,.55,2))])


if __name__=="__main__":
    print(json.dumps(benchmark(),indent=2,allow_nan=False))
