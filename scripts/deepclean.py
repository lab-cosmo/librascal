from os import listdir, remove
from os.path import isfile, join
from shutil import rmtree
for f in listdir('./'):
    if (isfile(f)):
        remove(f)
    elif (f != 'external') and (f != 'third-party'):
        rmtree(f)
