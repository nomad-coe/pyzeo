import sys
import platform
from distutils.ccompiler import new_compiler
from distutils.sysconfig import customize_compiler
from setuptools import setup, find_packages
from setuptools.extension import Extension
from subprocess import getoutput


# Check python version
if sys.version_info[:2] < (3, 6):
    raise RuntimeError("Python version >= 3.6 required.")

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
# with the addition of the -fPIC flag which is needed to build a shared library.
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

common_srcfiles = [
    'src/voro++/src/voro++.cc',
    'src/networkanalysis.cc',
    'src/networkstorage.cc',
    'src/networkinfo.cc',
    'src/network.cc',
    'src/net.cc',
    'src/mindist.cc',
    'src/geometry.cc',
    'src/OMS.cc',
    'src/voronoicell.cc',
    'src/v_network.cc',
    'src/graphstorage.cc',
    'src/channel.cc',
	'src/symmetry.cc',
    'src/ray.cc',
    'src/rmsd.cc',
    'src/material.cc',
    'src/psd.cc',
    'src/area_and_volume.cc',
    'src/networkaccessibility.cc',
    'src/string_additions.cc',
]
netstorage_srcfiles = [
    'src/zeoplusplus/netstorage'+ext,
    'src/networkio.cc',
    'src/grid.cc',
    'src/symbcalc.cc',
] + common_srcfiles
voronoicell_srcfiles = [
    'src/zeoplusplus/voronoicell'+ext,
] + common_srcfiles
highaccuracy_srcfiles = [
    'src/zeoplusplus/high_accuracy'+ext,
	'src/networkio.cc',
	'src/grid.cc',
	'src/symbcalc.cc',
    'src/sphere_approx.cc',
] + common_srcfiles
areavol_srcfiles = [
    'src/zeoplusplus/area_volume'+ext,
    'src/networkio.cc',
    'src/grid.cc',
    'src/symbcalc.cc',
] + common_srcfiles
cluster_srcfiles = [
    'src/zeoplusplus/cluster'+ext,
    'src/cluster.cc',
    'src/sphere_approx.cc',
] + common_srcfiles
cycle_srcfiles = [
    'src/zeoplusplus/cycle'+ext,
    'src/cycle.cc',
    'src/sphere_approx.cc',
] + common_srcfiles
graphstorage_srcfiles = [
    'src/zeoplusplus/graphstorage'+ext,
] + common_srcfiles
netio_srcfiles = [
    'src/zeoplusplus/netio'+ext,
    'src/networkio.cc', 
    'src/networkinfo.cc',
    'src/string_additions.cc',
    'src/grid.cc',
    'src/mindist.cc',
    'src/symbcalc.cc',
    'src/symmetry.cc',
    'src/networkstorage.cc',
    'src/geometry.cc',
    'src/net.cc',
    'src/rmsd.cc',
]
channel_srcfiles = [
    'src/zeoplusplus/channel'+ext,
    'src/channel.cc',
]
geometry_srcfiles = [
    'src/zeoplusplus/geometry'+ext,
    'src/geometry.cc'
]
netinfo_srcfiles = [
    'src/zeoplusplus/netinfo'+ext,
    'src/networkinfo.cc',
]
psd_srcfiles = [
    'src/zeoplusplus/psd'+ext,
    'src/psd.cc',
]
extensions = [
    Extension(
        "zeoplusplus.netstorage",
        sources=netstorage_srcfiles, 
        include_dirs=includedirs,
        extra_compile_args=cpp_extra_compile_args,
        extra_link_args=cpp_extra_link_args,
        language=language
    ),
    Extension(
        "zeoplusplus.voronoicell",
        sources=voronoicell_srcfiles,
        include_dirs=includedirs,
        extra_compile_args=cpp_extra_compile_args,
        extra_link_args=cpp_extra_link_args,
        language=language
    ),
    Extension(
        "zeoplusplus.high_accuracy", 
        sources=highaccuracy_srcfiles,
        include_dirs=includedirs,
        extra_compile_args=cpp_extra_compile_args,
        extra_link_args=cpp_extra_link_args,
        language=language
    ),
    Extension(
        "zeoplusplus.area_volume", 
        sources=areavol_srcfiles,
        include_dirs=includedirs,
        extra_compile_args=cpp_extra_compile_args,
        extra_link_args=cpp_extra_link_args,
        language=language
    ),
    Extension(
        "zeoplusplus.cluster", 
        sources=cluster_srcfiles,
        include_dirs=includedirs,
        extra_compile_args=cpp_extra_compile_args,
        extra_link_args=cpp_extra_link_args,
        language=language
    ),
    Extension(
        "zeoplusplus.cycle", 
        sources=cycle_srcfiles,
        include_dirs=includedirs,
        extra_compile_args=cpp_extra_compile_args,
        extra_link_args=cpp_extra_link_args,
        language=language
    ),
    Extension(
        "zeoplusplus.graphstorage",
        sources=graphstorage_srcfiles,
        include_dirs=includedirs,
        extra_compile_args=cpp_extra_compile_args,
        extra_link_args=cpp_extra_link_args,
        language=language
    ),
    Extension(
        "zeoplusplus.netio",
        sources=netio_srcfiles,
        include_dirs=includedirs,
        extra_compile_args=cpp_extra_compile_args,
        extra_link_args=cpp_extra_link_args,
        language=language
    ),
    Extension(
        "zeoplusplus.channel", 
        sources=channel_srcfiles,
        include_dirs=includedirs,
        extra_compile_args=cpp_extra_compile_args,
        extra_link_args=cpp_extra_link_args,
        language=language
    ),
    Extension(
        "zeoplusplus.geometry", 
        sources=geometry_srcfiles,
        extra_compile_args=cpp_extra_compile_args,
        extra_link_args=cpp_extra_link_args,
        language=language
    ),
    Extension(
        "zeoplusplus.netinfo", 
        sources=netinfo_srcfiles,
        extra_compile_args=cpp_extra_compile_args,
        extra_link_args=cpp_extra_link_args,
        language=language
    ),
    Extension(
        "zeoplusplus.psd",
        sources=psd_srcfiles,
        include_dirs=includedirs,
        extra_compile_args=cpp_extra_compile_args,
        extra_link_args=cpp_extra_link_args,
        language=language
    ),
]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

setup(
    name='zeoplusplus',
    version='0.1.0',
    description="Python interface to Zeo++",
    long_description="Python interface to Zeo++",
    url="https://github.com/lauri-codes/zeoplusplus",
    author="Lauri Himanen",
    license="",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Programming Language :: C++",
        "Programming Language :: Cython",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
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
    python_requires=">=3.6",
)
