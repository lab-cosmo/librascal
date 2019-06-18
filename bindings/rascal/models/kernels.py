from ..lib.Kernels import CosineKernel

class Kernel(object):
    def __init__(self, representation, name='cosine', target_type='structure', **kwargs):

        hypers = dict(target_type=target_type)
        hypers.update(**kwargs)
        self.representation = representation
        if name == 'cosine':
            self._kernel = CosineKernel(hypers)
        else:
            raise NotImplementedError("The kernel: {}, has not been implemented".format(name))

    def __call__(self, X, Y=None):
        if Y is None:
            Y = X
        return self._kernel(self.representation, X, Y)
