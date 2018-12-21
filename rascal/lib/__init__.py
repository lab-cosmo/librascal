import sys
import os

sys.path.append(os.path.join(
    os.path.dirname(__file__), "../../bindings/"))

"""This is not an elegant solution and not pep8 conform (comment of Till)
bindings/ and /rascal/ should be the same.
"""

import _rascal
from _rascal import (StructureManager, FeatureManager, Adaptor,
                     RepresentationManager, utils, math)

from _rascal.utils import sparsification
