from .fps import fps
from ..lib._rascal.utils import ostream_redirect
from ..lib import utils
from .pool_worker import FactoryPool

ostream_redirect = utils.__dict__['ostream_redirect']


def is_notebook():
    from IPython import get_ipython
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return True   # Jupyter notebook or qtconsole
        elif shell == 'TerminalInteractiveShell':
            return False  # Terminal running IPython
        else:
            return False  # Other type (?)
    except NameError:
        return False      # Probably standard Python interpreter


# if is_notebook():
#     from tqdm import tqdm_notebook as tqdm_cs
#     ascii = False
# else:
#     from tqdm import tqdm as tqdm_cs
#     ascii = True
