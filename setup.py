import sys
from os.path import join
from distutils.sysconfig import get_python_lib
from setuptools import find_packages
from pkg_resources import parse_requirements
from skbuild import setup
from skbuild.constants import CMAKE_INSTALL_DIR

from bindings.librascal import __version__

with open("AUTHORS.txt") as fp:
    author = ", ".join(fp.readlines()).replace("\n", "").replace("\n", "")


with open("./requirements/minimal.txt", "r") as fp:
    install_requires = [
        str(requirement) for requirement in parse_requirements(fp.read())
    ]

with open("README.rst") as f:
    long_description = f.read()

setup(
    name="librascal",
    version=__version__,
    long_description_content_type="text/x-rst",
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
