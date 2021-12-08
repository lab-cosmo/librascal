import sys
import os
from os.path import join
from distutils.sysconfig import get_python_lib

from setuptools import find_packages
from skbuild import setup
from skbuild.constants import CMAKE_INSTALL_DIR
import re
from bindings.librascal import __version__

with open("AUTHORS.txt") as fp:
    author = ", ".join(fp.readlines()).replace("\n", "").replace("\n", "")

from pkg_resources import parse_requirements
with open("./requirements/minimal.txt", "r") as fp:
    install_requires  = [str(requirement) for requirement in parse_requirements(fp.read())]

from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.rst").read_text()

with open('README.rst') as f:
    long_description = f.read()

setup(
    name="librascal",
    version=__version__,
    long_description_content_type='text/x-rst',
    long_description=long_description,
    description="""A versatile and scalable computation of representations of \
    atomic structures for machine learning.""",
    author=author,
    author_email="michele.ceriotti@epfl.ch",
    license="LGPL-3.0-or-later",
    cmake_args=[
        "-DINSTALL_PATH:STRING=" + join(os.getcwd(), CMAKE_INSTALL_DIR()),
        "-DBUILD_EXAMPLES:BOOL=OFF",
    ],
    url="https://github.com/cosmo-epfl/librascal/",
    package_dir={"": "bindings"},
    packages=find_packages(where="bindings"),
    install_requires=install_requires,
    # include_package_data=True,
    package_data={"": ["lib/librascal.*"]},
    zip_safe=False,
)
