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

Example Static IR Spectrum of o-Nitrophenol
-------------

For the calculation of the IR spectrum, we need to define the ``input.txt`` file.
For o-Nitrophenol (o-NP), it can look as follows.

.. code-block:: input
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

In addition, we need the data to feed to vibrant. For calculating the static IR spectrum we need the coordinates of o-NP, 
the forces for computing the Hessian, and the dipole data for calculating the intensities. 
Check :doc:`IR spectrum documentation <IR_spec>` for more theoretical background.

.. raw:: html

   <div style="text-align: center; margin: 1em 0;">
     <img src="_static/index/o-NP.png" alt="o-NP" style="max-width: 300px; height: auto; border-radius: 6px;">
     <p style="font-size: 0.9em; color: #555;">Figure: o-Nitrophenol structure used for IR spectrum calculation</p>
   </div>

Data files for o-NP can be downloaded here:

* :download:`Geometry file <../test/IR_Static_NMA/o-NP.xyz>`
* :download:`Force data <../test/IR_Static_NMA/o-NP-force.data3>`
* :download:`Dipole data <../test/IR_Static_NMA/dipole_o-NP_free_static.xyz>`

After running ``vibrant``, multiple output files are produced.
The spectral information is written to ``result_static_ir.txt`` and can be plotted
with a simple Python script, for example:

.. code-block:: python

   import numpy as np
   import matplotlib.pyplot as plt

   filename = "result_static_ir.txt"
   data = np.loadtxt(filename, skiprows=1)
   x, y = data[:, 0], data[:, 1]
   # === Plot ===
   plt.figure(figsize=(6, 4))
   plt.plot(x, y, color="red", linewidth=1)
   plt.xlabel("Wavenumber (cm$^{-1}$)")
   plt.title("")
   ax = plt.gca()
   ax.yaxis.set_visible(False) 
   ax.spines['top'].set_visible(False)    
   ax.spines['right'].set_visible(False) 
   ax.spines['left'].set_visible(False)   

   plt.tight_layout()
   plt.show()

.. raw:: html

   <div style="text-align: center; margin-bottom: 1em;">
   <picture>
      <img alt="IR spectrum o-NP" src="_static/index/ir_spectrum_o-NP.png" style="max-width: auto; height: auto; display: inline-block;">
   </picture>
   </div>

--> more examples
