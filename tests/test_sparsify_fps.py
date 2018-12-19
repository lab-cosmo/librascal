# please change this to your local librascal folder
import sys, os
path = '/home/michele/lavoro/code/librascal/'
sys.path.insert(0, os.path.join(path, 'build/'))

from rascal.utils import fps #, fps_voronoi
import numpy as np

#np.random.seed(1234)
x = np.random.normal(size=(10000,10))

# returns indices of the selected structures, FPS distance list, Hausdorff distances of all the inputs to the FPS selections
ifps, dfps, lfps = fps(x, 2000, 0, method="simple")
print("ifps-direct", ifps[:3], ifps[-3:], ifps[998:1002])

# now we do this in two stages, by restarting
# I first select 1000
ifps, dfps, lfps = fps(x, 1000, 0, method="simple")
print("ifps-half", ifps[:3], ifps[-3:])

# I build a tuple containing the return arrays
res_tuple = (ifps, dfps, lfps)
# I ask for 1000 but provide the restart tuple. Selection will continue from point 1001
ifps, dfps, lfps = fps(x, 2000, 0, method="simple", restart=res_tuple)
print("ifps-restart", ifps[:3], ifps[-3:], ifps[998:1002])

# You can also provide only the selected indices
ifps, dfps, lfps = fps(x, 1000, 0, method="simple")
res_tuple = (ifps, np.zeros(0,float), np.zeros(0,float))
ifps, dfps, lfps = fps(x, 2000, 0, method="simple", restart=res_tuple)
print("ifps-restart", ifps[:3], ifps[-3:], ifps[998:1002])

# You can also add new points (in this case we reuse the same feature matrix)
# and continue the selection
ifps, dfps, lfps = fps(x, 1000, 0, method="simple")
xx = np.concatenate((x[ifps], x))
res_tuple = (np.asarray(range(1000)), np.zeros(0,float), np.zeros(0,float))
ifps, dfps, lfps = fps(xx, 2000, 0, method="simple", restart=res_tuple)
print("ifps-restart", ifps[:3], ifps[-3:], ifps[998:1002])
