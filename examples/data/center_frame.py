from ase.io import read, write
import numpy as np
import sys
filename=sys.argv[1]
frame = read(filename)
frame.set_positions(np.subtract(frame.positions, np.mean(frame.positions, axis=0)))
write(filename, frame)
