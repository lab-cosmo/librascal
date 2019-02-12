from ..lib import NeighbourList

_neighbourlist_list = ["centers","neighbourlist","strict","maxorder","halflist","fulllist"]

_neighbourlists = {}
for k,v in NeighbourList.__dict__.items():
    if "make_adapted_manager" in k or "make_structure_manager" in k:
        name = k.lower().replace('make_adapted_manager_','').replace('make_structure_manager_','')
        _neighbourlists[name] = v

def NeighbourListFactory(name,*args):
    if name not in _neighbourlists:
        raise NameError('The neighbourlist factory {} has not been registered. The available combinations are: {}'.format(name,list(_neighbourlists.keys())))
    return _neighbourlists[name](*args)