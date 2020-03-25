import sys
import os
from os.path import join
from distutils.sysconfig import get_python_lib
from setuptools import find_packages

try:
    from skbuild import setup
    from skbuild.constants import CMAKE_INSTALL_DIR

except ImportError:
    print('scikit-build is required to build from source.', file=sys.stderr)
    print('Please run:', file=sys.stderr)
    print('', file=sys.stderr)
    print('  python3 -m pip install scikit-build')
    sys.exit(1)

with open('./requirements_pip.txt', 'r') as fp:
    requirements = list(filter(lambda x: '#' not in x, (line.strip() for line in fp)))

setup(
    name="rascal",
    version="0.3.9",
    description="""A versatile and scalable computation of representations of
atomic structures for machine learning.""",
    author='Félix Musil',
    license="LGPL-3",
    cmake_args=[
      '-DINSTALL_PATH:STRING='+join(os.getcwd(),CMAKE_INSTALL_DIR()),
      '-DBUILD_EXAMPLES:BOOL=OFF'
    ],
    package_dir={"": "bindings"},
    packages=find_packages(where='bindings'),
    install_requires=requirements
)