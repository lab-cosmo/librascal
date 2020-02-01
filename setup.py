import sys
import os
from os.path import join
from distutils.sysconfig import get_python_lib

try:
    from skbuild import setup
    from skbuild.constants import CMAKE_INSTALL_DIR

except ImportError:
    print('scikit-build is required to build from source.', file=sys.stderr)
    print('Please run:', file=sys.stderr)
    print('', file=sys.stderr)
    print('  python -m pip install scikit-build')
    sys.exit(1)


setup(
    name="rascal",
    version="0.3.2",
    description="""A versatile and scalable computation of representations of
atomic structures for machine learning.""",
    author='Félix Musil',
    license="LGPL-3",
    # packages=['rascal'],
    # packages_dir={'rascal':'./bindings/rascal'},
    cmake_args=[
      '-DINSTALL_PATH:STRING='+join(os.getcwd(),CMAKE_INSTALL_DIR()),
      '-DBUILD_EXAMPLES:BOOL=OFF'
    ]
)