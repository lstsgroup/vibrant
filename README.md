# vibrant

## Installation 

### Quick install

The CMake build framework is used in order to compile the `vibrant` code. To install it type the following in the root directory of the repository and make sure that you have a working GNU fortran compiler (e.g. `gfortran`):
```bash 
mkdir build && cd build 
cmake ..
make -j 5
```
Optional: Test if the installation is working as expected by typing 
```bash
ctest 
```
in the build directory to run the regression tests.

### Install options

Options can be passed to cmake via an `initial_cache.cmake` file:
```bash 
cmake -C ../initial_cache.cmake ..
```
e.g to specify a dependency path (for more options see `initial_cache.example.cmake`):
```cmake 
set(BLAS_LIBRARY_DIR "/path/to/BLAS" CACHE PATH "Path to BLAS library directory")
```
Within this file the dependency paths, the compiler, the compiler flags etc. can be specified. For more options see `initial_cache.example.cmake`.

### Dependencies

| Dependency | Version | Comment |
|---|---|---|
LAPACK (e.g. OpenBLAS) | (>0.3) | Required, if no path provided then cmake tries to find the library |
[FFTW](https://www.fftw.org/) | >3.3 | Required, if no path provided then installed together with vibrant | 
GreenX | v2.0 | Required, if no path provided then installed together with vibrant |
Python | >3.10 | Required, only if `ENABLE_REGRESSION_TESTS` is `ON` |


## Regression test 
The regression tests can be run by running the following in the build directory after `vibrant` was build:
```
ctest 
```

## code style

fprettify config:
```
indent = 4
strict-indent = True
case = [2, 2, 2, 2]
whitespace = 2
whitespace-relational = False
line-length = 300
```

## Website

#### install sphinx 
For building the website, the following python modules are required:
```bash
pip install -U sphinx
pip install sphinx-wagtail-theme
pip install --upgrade myst-parser
```

#### build the website 
enable website generation in cmake:
```bash 
mkdir build && cd build
cmake -DENABLE_DOCS=ON ..
```
build the website with:
```bash 
make docs
```
Display the website by opening the file `vibrant/build/html/index.html` in a browser.


## Contributing

Contributions to vibrant are highly appreciated! If you consider contributing [check out CONTRIBUTING.md]()

## License note

Vibrant is distributed under the Apache 2.0 license. Please note that, while the source code is licensed permissively under Apache 2.0, any compiled binary is subject to the GPL 2.0 license due to its dependency on FFTW.
