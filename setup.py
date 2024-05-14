import sys
import platform
from distutils.ccompiler import new_compiler
from distutils.sysconfig import customize_compiler
from setuptools import setup, find_packages
from setuptools.extension import Extension
from subprocess import getoutput


# Check python version
if sys.version_info[:2] < (3, 8):
    raise RuntimeError("Python version >= 3.8 required.")

# The recommendation
# (https://cython.readthedocs.io/en/latest/src/userguide/source_files_and_compilation.html#distributing-cython-modules)
# is to distribute Cython compiled source files with the package. This is why
# the Cython compilation step is disabled by default here.
USE_CYTHON = False
ext = '.pyx' if USE_CYTHON else '.cpp'
language = "c++"
includedirs = []
cpp_extra_link_args = []

# These default compile flags mimic the flags used in the Zeo++/Voro++ Makefile
cpp_extra_compile_args = [
    "-Wall",
    "-ansi",
    "-pedantic",
    "-O3",
]

def using_clang():
    """Will we be using a clang compiler?
    """
    compiler = new_compiler()
    customize_compiler(compiler)
    compiler_ver = getoutput("{0} -v".format(compiler.compiler[0]))
    return "clang" in compiler_ver

# Needed to specify C++ runtime library on OSX. This solution is replicated
# from the setup.py of mdanalysis
if platform.system() == "Darwin" and using_clang():
    cpp_extra_compile_args.append("-stdlib=libc++")
    cpp_extra_compile_args.append("-mmacosx-version-min=10.9")
    cpp_extra_link_args.append("-stdlib=libc++")
    cpp_extra_link_args.append("-mmacosx-version-min=10.7")

extensions = [
    Extension(
        "pyzeo.extension", 
        sources=[
            'src/pyzeo/extension'+ext,
            'src/area_and_volume.cc',
            'src/channel.cc',
            'src/cluster.cc',
            'src/cycle.cc',
            'src/grid.cc',
            'src/geometry.cc',
            'src/graphstorage.cc',
            'src/voro++/src/voro++.cc',
            'src/net.cc',
            'src/networkaccessibility.cc',
            'src/networkanalysis.cc',
            'src/networkstorage.cc',
            'src/networkinfo.cc',
            'src/network.cc',
            'src/networkio.cc', 
            'src/material.cc',
            'src/mindist.cc',
            'src/OMS.cc',
            'src/psd.cc',
            'src/sphere_approx.cc',
            'src/string_additions.cc',
            'src/symbcalc.cc',
            'src/symmetry.cc',
            'src/ray.cc',
            'src/rmsd.cc',
            'src/voronoicell.cc',
            'src/v_network.cc',
        ],
        include_dirs=includedirs,
        extra_compile_args=cpp_extra_compile_args,
        extra_link_args=cpp_extra_link_args,
        language=language
    )
]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

setup(
    name='pyzeo',
    version='0.1.5',
    description="Python interface to Zeo++",
    long_description="Python interface to Zeo++",
    url="https://github.com/nomad-coe/pyzeo",
    author="Lauri Himanen",
    license="",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Programming Language :: C++",
        "Programming Language :: Cython",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3 :: Only",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering"
    ],
    # In the Cython build guide it is said that zip_safe should be disabled
    # when building with setuptools
    packages=find_packages('src'),
    package_dir={'': 'src'},
    zip_safe=False,
    ext_modules=extensions,
    keywords="zeo++ porous materials science",
    python_requires=">=3.8",
)
