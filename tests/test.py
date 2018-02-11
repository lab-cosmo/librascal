import sys
from glob import glob
sys.path.insert(0,glob('./build/lib.*')[0]+'/')
import numpy as np
from scipy.spatial.distance import cdist
import proteus as p

aa = np.array([[0.66196274, 0.91671878, 0.85604587],
       [0.06124809, 0.30113883, 0.5895188 ],
       [0.71365461, 0.69070415, 0.02928229],
       [0.41993124, 0.88230171, 0.82275001],
       [0.14505813, 0.32507217, 0.02294067]])
bb = np.array([[0.95522114, 0.66015242, 0.68513688],
       [0.53203634, 0.32235253, 0.38465644]])


assert np.allclose(cdist(aa,bb,metric='euclidean'),p.cdist(p.Matrix(aa),p.Matrix(bb)))


