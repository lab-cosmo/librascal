import sys
from os import listdir, remove, getcwd
from os.path import isfile, join
from shutil import rmtree
def deepclean():
    for f in listdir('./'):
        if (isfile(f)):
            remove(f)
        elif (f != 'external') and (f != 'third-party'):
            rmtree(f)

if __name__ == '__main__':
    if len(sys.argv) > 1:
        if (sys.argv[1] == "deepclean"):
            deepclean()
        else:
            print("Command "+sys.argv[1]+" not known.")
    else:
        print("Please provide an argument which functionality you want to use.")

