.. vibrant documentation master file, created by
   sphinx-quickstart on Mon Sep 23 12:49:43 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. raw:: html

   <div style="text-align: center; margin-bottom: 1em;">
   <picture>
      <img alt="Vibrant logo" src="_static/all_img/logo_vibrant.svg" style="height:140px; display: inline-block;">
   </picture>
   </div>

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   Installation
   Usage
   Keyword_Glossary
   Cite
   Contribute

Quick Start 
===========
Welcome to the ``vibrant`` documentation webpage.

The ``vibrant`` code is a post-processing tool for calculating vibrational spectra from ab initio data.

Installation 
-------------
See the :doc:`Installation <Installation>` section for fine-grained installation instructions.

After installation, a compiled binary should be available in the build directory as ``build/vibrant`` can be invoked e.g. from the root directory by the command

.. code-block:: bash
   
   build/vibrant input.txt

For more details on the ``input.txt`` format, see the :doc:`Usage <Usage>` section and the :doc:`Keyword Glossary <Keyword_Glossary>`.
Information about the required data and file formats can be found on the :doc:`File Formats <file_formats>` (in Usage).

Example static IR Spectrum of o-Nitrophenol
-------------

For the calculation of the IR spectrum, we need to define the ``input.txt`` file.
For o-NP, it can look as follows.
.. code-block:: input.txt

    &global
      spectra IR
      temperature 300
      fwhm 5
   &end global
   &system
      filename o-NP.xyz
   &end system
   &static
      diag_hessian y
      &hessian
         force_file o-NP-force.data3
      &end hessian
      displacement 0.001
   &end static
   &dipoles
      type_dipole berry
      dip_file dipole_o-NP_free_static.xyz
   &end dipoles
In addtion we need the data files. For calculating the static IR spectrum here we need the XYZ coordinates of the structrue, 
the force data for computing the Hessian, and the dipole data for calculating the intensities. Check :doc:`IR spectrum documentation <IR_spec>` for more theoretical background.
Download datafiles
.. raw:: html
   <div style="display: flex; align-items: flex-start; justify-content: space-between; margin-bottom: 1em;">
   <div style="flex: 0 0 250px; text-align: right;">
      <img src="_static/index/o-NP.png" alt="o-NP" style="max-width: 100%; height: auto; border-radius: 6px;">
   </div>
   </div>
After running the vibrant exectuable we get mutiple outputs the spectra information is written to result_static_ir an can be plotten with a simple python scirpt e.g

.. code-block:: plotting
   import numpy as np
   import matplotlib.pyplot as plt

   filename = "result_static_ir.txt"
   data = np.loadtxt(filename, skiprows=1)
   x, y = data[:, 0], data[:, 1]

   # === Plot ===
   plt.figure(figsize=(6, 4))
   plt.plot(x, y, color="red", linewidth=1)
   plt.xlabel("Wavenumber (cm$^{-1}$)")
   plt.ylabel("Intensity")
   plt.title("IR Spectrum")
   ax = plt.gca()
   ax.set_yticks([])      
   ax.set_yticklabels([])
   plt.tight_layout()
   plt.show()
   
Resuting spectra

.. raw:: html
   <div style="display: flex; align-items: flex-start; justify-content: space-between; margin-bottom: 1em;">
   <div style="flex: 0 0 250px; text-align: right;">
      <img src="_static/index/ir_spectrum_o-NP.png" alt="o-NP" style="max-width: 100%; height: auto; border-radius: 6px;">
   </div>
   </div>
--> more examples
