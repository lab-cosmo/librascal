# from ..lib.Kernels import CosineKernel
from ..lib._rascal.Models.Kernels import Kernel as Kernelcpp
import json

class Kernel(object):
    def __init__(self, representation, name='Cosine', target_type='structure', **kwargs):

        hypers = dict(name=name,target_type=target_type)
        hypers.update(**kwargs)
        hypers_str = json.dumps(hypers)
        self.representation = representation._representation

        self._kernel = Kernelcpp(hypers_str)


    def __call__(self, X, Y=None):
        if Y is None:
            Y = X
        return self._kernel.compute(self.representation, X, Y)
