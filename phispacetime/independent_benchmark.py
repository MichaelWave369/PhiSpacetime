import json
from .independent import compare


def benchmark():
    return dict(scope="independent fitted polar scalar meshes; bounded domain, not free-space scattering",
                rows=[compare(n,mode,end,factor) for n,mode,end,factor in
                      ((17,1,.24,1),(33,1,.24,1),(65,1,.24,1),
                       (33,2,.24,1),(33,1,.24,2),(33,1,.33,1))])


if __name__=="__main__":
    print(json.dumps(benchmark(),indent=2,allow_nan=False))
