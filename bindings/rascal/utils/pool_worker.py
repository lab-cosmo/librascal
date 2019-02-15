from multiprocessing.dummy import Pool as ThreadPool

def FactoryPool(method='thread',n_workers=1,disable_pbar=False):
    if n_workers < 2:
        pool = SerialPool()
    elif method == 'thread':
        pool = ThreadPool(n_workers)
    else:
        raise NameError('Worker pool {} is not implemented'.format(method))
    return pool


class SerialPool(object):
    def __init__(self,*args,**kwargs):
        pass

    @staticmethod
    def starmap(func, itarables):
        results = []
        for itarable in itarables:
            results.append(func(*itarable))
        return results

    @staticmethod
    def map(func, itarable):
        results = []
        for item in itarable:
            results.append(func(item))
        return results

    def close(self):
        pass

    def join(self):
        pass