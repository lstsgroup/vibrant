# vibrant

## Installation 

The CMake build framework is used in order to compile the `vibrant` code. Currently it has two dependencies (`fftw-3.3.10`, `greenX`) that are fetched and compiled automatically when the `vibrant` software is installed. To install it type the following in the root directory of the repository:
```bash 
mkdir build && cd build 
cmake ..
make -j 10
```
To change to a specific compiler, run the above CMake command like this:
```bash 
FC=gfortran cmake .. 
```

## Regression test 
The regression tests can be run by running the following in the build directory after `vibrant` was build:
```
ctest 
```