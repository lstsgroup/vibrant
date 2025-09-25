.. vibrant documentation master file, created by
   sphinx-quickstart on Mon Sep 23 12:49:43 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. raw:: html

   <style>
     @media (prefers-color-scheme: dark) {
       .vibrant-logo-light { display: none !important; }
       .vibrant-logo-dark { display: inline !important; }
     }
     @media (prefers-color-scheme: light), (prefers-color-scheme: no-preference) {
       .vibrant-logo-dark { display: none !important; }
       .vibrant-logo-light { display: inline !important; }
     }
   </style>
   <img src="_static/all_img/logo_vibrant.svg" alt="Vibrant Logo" class="vibrant-logo-light" style="height:140px; float:left;">
   <img src="_static/all_img/logo_vibrant_dark.svg" alt="Vibrant Logo (Dark)" class="vibrant-logo-dark" style="height:140px; float:left;">

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
