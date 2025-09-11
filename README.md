# vibrant

## Installation 

The CMake build framework is used in order to compile the `vibrant` code. Currently it has two dependencies (`fftw-3.3.10`, `greenX`) that are fetched and compiled automatically when the `vibrant` software is installed. To install it type the following in the root directory of the repository:
```bash 
mkdir build && cd build 
cmake ..
make -j 10
```
The default compiler is `gfortran`. To change to a specific compiler, run the above CMake command like this:
```bash 
cmake -DCMAKE_Fortran_COMPILER=f90 .. 
```

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