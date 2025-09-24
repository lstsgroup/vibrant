<img src="docs/_static/all_img/logo_vibrant.svg" height="120">

---

[![Regression Tests](https://github.com/lstsgroup/vibrant/actions/workflows/action.yml/badge.svg?branch=main)](https://github.com/lstsgroup/vibrant/actions/workflows/action.yml)


## Quick install

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

For more a more fine-grained installation please refer to the [documentation of vibrant](https://lstsgroup.github.io/vibrant/Installation.html).

## Documentation

For more Information about the vibrant code please consult the extensive [documentation website](https://lstsgroup.github.io/vibrant/).

## Contributing

Contributions to vibrant are highly appreciated! If you consider contributing [check out CONTRIBUTING.md](https://github.com/lstsgroup/vibrant/blob/main/CONTRIBUTING.md)

## License note

Vibrant is distributed under the Apache 2.0 license. Please note that, while the source code is licensed permissively under Apache 2.0, any compiled binary is subject to the GPL 2.0 license due to its dependency on FFTW.
