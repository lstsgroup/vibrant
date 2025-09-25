.. vibrant documentation master file, created by
   sphinx-quickstart on Mon Sep 23 12:49:43 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. raw:: html

   <div style="text-align: left; margin-bottom: 1em;">

     <picture>
       <source media="(prefers-color-scheme: dark)" srcset="_static/all_img/logo_vibrant_dark.svg">
       <source media="(prefers-color-scheme: light)" srcset="_static/all_img/logo_vibrant.svg">
       <img alt="Vibrant logo" src="_static/all_img/logo_vibrant.svg" style="height:140px; display: block;">
     </picture>

   </div>

Documentation of Ekin's very cool code.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Installation
   Usage
   Cite
   Contribute


Quick Start 
===========

The CMake build framework is used in order to compile the ``vibrant`` code. To install it type the following in the root directory of the repository and make sure that you have a working GNU Fortran compiler (e.g. ``gfortran``):

.. code-block:: bash

   mkdir build && cd build
   cmake ..
   make -j 5
   
Optional: Test if the installation is working as expected by typing 

.. code-block:: bash
   
   ctest 

in the build directory to run the regression tests.

For more a more fine-grained installation please refer to the `documentation of vibrant <https://lstsgroup.github.io/vibrant/Installation.html>`.
