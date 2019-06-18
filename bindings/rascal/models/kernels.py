# from ..lib.Kernels import CosineKernel
from ..lib._rascal.Models.Kernels import CosineKernel
import json

class Kernel(object):
    def __init__(self, representation, name='cosine', target_type='structure', **kwargs):

        hypers = dict(target_type=target_type)
        hypers.update(**kwargs)
        hypers_str = json.dumps(hypers)
        self.representation = representation._representation
        if name == 'cosine':
            self._kernel = CosineKernel(hypers_str)
        else:
            raise NotImplementedError("The kernel: {}, has not been implemented".format(name))

    def __call__(self, X, Y=None):
        if Y is None:
            Y = X
        return self._kernel.compute(self.representation, X, Y)
