import json
from .pockets import run


def benchmark():
    return dict(scope="static forced scalar PDE; source-informed refinement; no dynamic AMR",
                rows=[dict(method=name,**run(arity,budget,adaptive,localized))
                      for localized in (False,True) for budget in (16,64,256,1024)
                      for name,arity,adaptive in (("uniform",2,False),("quadtree",2,True),("nine_way",3,True))])


if __name__=="__main__":
    print(json.dumps(benchmark(),indent=2,allow_nan=False))
