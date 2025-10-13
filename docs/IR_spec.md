# IR Spectra

The infrared (IR) spectrum calculation in Vibrant is based on the post-processing of the dipole moments or Wannier centers. The IR intensities are computed using the following equation:

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

## Install options

Options can be passed to cmake via an `initial_cache.cmake` file:
```bash 
cmake -C ../initial_cache.cmake ..
```
e.g to specify a dependency path (for more options see `initial_cache.example.cmake`):
```cmake 
set(BLAS_LIBRARY_DIR "/path/to/BLAS" CACHE PATH "Path to BLAS library directory")
```
Within this file the dependency paths, the compiler, the compiler flags etc. can be specified. For more options see `initial_cache.example.cmake`.

## Dependencies

| Dependency | Version | Comment |
|---|---|---|
LAPACK (e.g. OpenBLAS) | (>0.3) | Required, if no path provided then cmake tries to find the library |
[FFTW](https://www.fftw.org/) | >3.3 | Required, if no path provided then installed together with vibrant | 
[GreenX](https://github.com/nomad-coe/greenX) | v2.0 | Required, if no path provided then installed together with vibrant |
Python | >3.10 | Required, only if `ENABLE_REGRESSION_TESTS` is `ON` |


## Regression test 
The regression tests can be run by running the following in the build directory after `vibrant` was build:
```
ctest 
```
