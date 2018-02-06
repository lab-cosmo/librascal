import sys
from glob import glob
sys.path.insert(0,glob('./build/lib.*')[0]+'/')

import proteus as p

#assert p.__version__ == '0.0.1'
assert p.add(1, 2) == 3
assert p.subtract(1, 2) == -1
