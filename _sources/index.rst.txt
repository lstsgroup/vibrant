.. vibrant documentation master file, created by
   sphinx-quickstart on Mon Sep 23 12:49:43 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.



   

Overview 
===========
The ``vibrant`` code is a post-processing tool for calculating vibrational spectra from atomistic simulation data.

.. raw:: html

   <div style="text-align: center; margin-bottom: 1em;">
      <picture>
         <img alt="Vibrant logo" src="_static/all_img/logo_vibrant.svg" style="height:140px; display: inline-block;">
      </picture>
   </div>
Features
===========
   * **Vibrational frequencies, IR and Raman spectra** from static calculations or molecular dynamics (MD) simulations.  
   * **Resonance Raman spectra** via efficient Padé approximation or explicit frequency sampling.  
   * Calculation of **Absorption spectra** alongside vibrational data for comprehensive analysis. 
   * **Dissection of MD-based spectra** into user-defined molecular fragments.
   * **Fast and efficient** post-processing of large systems and long trajectories.
   * Applicable to **single molecules, liquids and materials**.

Compatibility and Applications
===========

Provided that the required quantities (e.g., forces, dipole moments, polarizabilities) are available in the expected format (see :doc:`File Formats <file_formats>`) , ``vibrant`` is compatible with outputs from multiple electronic structure packages.

One example of a popular electronic structure package is `CP2K <https://www.cp2k.org/about>`_. For users performing calculations with CP2K, tutorials and documentation for setting up the required electronic structure calculations can be found in the official CP2K `manual <https://manual.cp2k.org/trunk/>`_ and `exercises <https://www.cp2k.org/exercises>`_.

So far, ``vibrant`` has been used in several publications to generate vibrational or electronic spectra, where the underlying electronic structure calculations were carried out using either the CP2K or FHI-aims program packages. Examples include:

1. Bas, E. E.; Garcia Alvarez, K. M.; Schneemann, A.; Heine, T.; Golze, D. Robust computation and analysis of vibrational spectra of layered framework materials including host–guest interactions. *J. Chem. Theory Comput.* 2024, *20* (21), 9547–9561. https://doi.org/10.1021/acs.jctc.4c01021
2. Feuerstein, L.; Bas, E. E.; Golze, D.; Heine, T.; Oschatz, M.; Weidinger, I. M. Nitrile groups as build-in molecular sensors for interfacial effects at electrocatalytically active carbon–nitrogen materials. *ACS Appl. Mater. Interfaces* 2025, *17* (16), 23996–24004. https://doi.org/10.1021/acsami.5c02366
3. Leucke, M.; Panadés-Barrueta, R. L.; Bas, E. E.; Golze, D. Analytic continuation component of the GreenX library: Robust padé approximants with symmetry constraints. *J. Open Source Softw.* 2025, *10* (109), 7859. https://doi.org/10.21105/joss.07859

For publication 1, the input and output files of the electronic structure calculations carried out with CP2K/FHI-aims are available in the corresponding `Zenodo repository <https://zenodo.org/records/11065943>`_.

Quick Start
===========
To get started with ``vibrant``, see the :doc:`Quick Start <Quick_Start>` section of the documentation.


Documentation
===========
.. toctree::
   :maxdepth: 1

   Quick_Start
   Installation
   Usage
   Keyword_Glossary
   Examples
   Cite
   Contribute
