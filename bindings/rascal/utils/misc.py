"""Miscellaneous useful utilities"""

import logging

LOGGER = logging.getLogger(__name__)


def is_notebook():
    """Is this being run inside an IPython Notebook?"""
    from IPython import get_ipython

    try:
        shell = get_ipython().__class__.__name__
        if shell == "ZMQInteractiveShell":
            return True  # Jupyter notebook or qtconsole
        elif shell == "TerminalInteractiveShell":
            return False  # Terminal running IPython
        else:
            return False  # Other type (?)
    except NameError:
        return False  # Probably standard Python interpreter


class tqdm_nop:
    """A simple no-op class to replace tqdm if it is not available"""

    def __init__(self, iterable, **kwargs):
        LOGGER.warn("tqdm not available")
        LOGGER.warn("(tried to call tqdm with args: {:s})".format(str(kwargs)))
        self.iterable = iterable

    def __iter__(self):
        return iter(self.iterable)

    def update(self):
        pass

    def close(self):
        pass
