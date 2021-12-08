import sys
import os
import re
from os.path import join
from distutils.sysconfig import get_python_lib
from setuptools import find_packages
from pkg_resources import parse_requirements
from skbuild import setup
from skbuild.constants import CMAKE_INSTALL_DIR

version = re.search(
    r'__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
    open("bindings/librascal/__init__.py").read(),
).group(1)

with open("README.rst") as f:
    long_description = f.read()

with open("AUTHORS.txt") as fp:
    author = ", ".join(fp.readlines()).replace("\n", "").replace("\n", "")

with open("./requirements/minimal.txt", "r") as fp:
    install_requires = [
        str(requirement) for requirement in parse_requirements(fp.read())
    ]

setup(
    name="librascal",
    version=version,
    long_description_content_type="text/x-rst",
    long_description=long_description,
    description="""A versatile and scalable computation of representations of \
    atomic structures for machine learning.""",
    author=author,
    author_email="michele.ceriotti@epfl.ch",
    cmake_args=[
        "-DINSTALL_PATH:STRING=" + join(os.getcwd(), CMAKE_INSTALL_DIR()),
        "-DBUILD_EXAMPLES:BOOL=OFF",
    ],
    url="https://github.com/cosmo-epfl/librascal/",
    package_dir={"": "bindings"},
    packages=find_packages(where="bindings"),
    install_requires=install_requires,
    package_data={"": ["lib/librascal.*"]},
    zip_safe=False,
    license="LGPLv3+",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: C++",
        "Topic :: Scientific/Engineering",
        "Operating System :: Unix",
    ],
)
