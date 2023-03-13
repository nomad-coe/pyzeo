# pyzeo
Python wrapper for the [Zeo++ library](http://zeoplusplus.org). Based on the latest released version 0.3.

## Installation

### pip
```sh
pip install pyzeo
```

### From source

```sh
git clone https://github.com/nomad-coe/pyzeo
cd pyzeo
pip install .
```

### Cython
By default the package comes with pre-built Cython binding code. The cython
wrapper definitions live inside src/pyzeo. These bindings can be recreated
by first setting `USE_CYTHON=True` in setup.py, and then recreating the bindings
with:

```sh
python setup.py build_ext --inplace --force
```

Remember to disable cython afterwards by setting `USE_CYTHON=False` in setup.py.

### Common issues

- **fatal error: Python.h: No such file or directory**: The package depends on
   C/C++ extensions that are compiled during the setup. For the compilation to
   work you will need to install the *pythonX.X-dev*-package, where X.X is the
   python version you use. E.g. for python 3.9 on Ubuntu this package could be
   installed with:

   .. code-block:: sh

       sudo apt install python3.9-dev

 - **Installation errors on MacOS**: The package depends on C++ extensions that
   are compiled during the setup. If experiencing problems with setup on MacOS,
   you may need to install the Xcode Command Line tools package. This can be
   done with:

   .. code-block:: sh

       xcode-select --install
