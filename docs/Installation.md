# Installation 

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

## Usage 

The compiled binary is now available in the build directory as `build/vibrant` and can be invoked e.g. from the root directory by the command
```bash
build/vibrant input.txt
```
with `input.txt` as a input file with the `vibrant` settings.




## Installing on macOS
> This was tested on macOS Tahoe 26.2 (M4 Chip), Homebrew 5.1.11, GCC 15.2.0_1, OpenBLAS 0.3.33, FFTW 3.3.11, CMake 4.3.2

Use paket manager `homebrew` to install GCC (including the `gfortran` compiler), OpenBLAS, FFTW and cmake:
```bash
brew install gcc openblas fftw cmake
```
Uncomment the following lines in the file `initial_cache.example.cmake` to let the build system know where the libraries are:
```cmake
# specifying the fortran compiler
set(CMAKE_Fortran_COMPILER "gfortran" CACHE FILEPATH "Fortran compiler")

# specifying library paths installed with Homebrew:
set(FFTW_LIBRARY_DIR "/opt/homebrew/opt/fftw/lib/" CACHE PATH "Path to FFTW library directory")
set(FFTW_INCLUDE_DIR "/opt/homebrew/opt/fftw/include/" CACHE PATH "Path to FFTW include directory")
set(BLA_VENDOR "OpenBLAS" CACHE STRING "BLAS vendor")
set(CMAKE_PREFIX_PATH "/opt/homebrew/opt/openblas" CACHE PATH "Path to OpenBLAS installation")
```
Now you can install vibrant with the following commands:
```bash
mkdir build && cd build
cmake -C ../initial_cache.example.cmake ..
make -j 5
```
Optionally you can run the application tests to verify a successfull installation:
```
ctest
``` 


