---
title: 'Vibrant: A Post-Processing Tool for Static and MD-based Vibrational Spectroscopy Simulations'
tags:
  - Fortran
  - spectroscopy
  - Raman
  - IR
authors:
  - name:
      given-names: Ekin
      surname: Winogradow
    affiliation: '1'

  - name: 
      given-names: Johannes
      surname: Scheffler
    affiliation: '1'

  - name: 
      given-names: Moritz
      surname: Leucke
    affiliation: '1'

  - name: 
      given-names: Thomas
      surname: Heine
    affiliation: '1, 3, 4'

  - name: 
      given-names: Dorothea
      surname: Golze
    affiliation: '1, 2'
affiliations:
  - index: 1
    name: Technische Universität Dresden, Germany
    ror: 042aqky30
  - index: 2
    name: Julius-Maximilians-Universität Würzburg, Germany
    ror: 00fbnyb24
  - index: 3
    name: Helmholtz-Zentrum Dresden-Rossendorf, Centrum for Advanced Systems Understanding, CASUS, Germany
  - index: 4
    name: Yonsei University and IBS center for nanomedicine, Republic of Korea
date: 15 September 2025
bibliography: paper.bib
---

# Summary 

In this work, we present Vibrant, a vibrational analysis program written in Fortran 2008 with interfaces to the FFTW[@frigo1998fftw] and GreenX (GX-AC)[@leucke2025analytic] libraries and parallelized using OpenMP for maximum efficiency. Vibrant enables the generation of vibrational spectra for a wide range of materials through post-processing atomic positions, velocities, forces, polarizability tensors and static or time-dependent dipole moments. Its functionalities include the computation of both static and molecular dynamics (MD)-based infrared (IR), Raman, resonance Raman (RR) and absorption spectra, as well as frequency analyses such as normal mode analysis and power spectrum calculations. Vibrant is distributed under the Apache 2.0 license and is available as open-source software on GitHub.

# Statement of need

Vibrational spectroscopy, particularly infrared (IR) and Raman, is a powerful and widely used technique for characterizing a wide range of materials from gaseous system or liquids to solids such as two-dimensional (2D)-materials. An important extension of Raman spectroscopy is resonance Raman (RR) spectroscopy, where the incident laser frequency matches the energy of an electronic transition in the system[@long2002raman]. As a result, intensities of the vibrational modes that are coupled to these electronic transitions get enhanced, enabling a more feasible characterization of certain vibrational modes in materials[@meinke2014uv; @goetz2019influence].

Computational simulation of vibrational spectroscopy is often far from straightforward. Most quantum chemistry packages that perform electronic structure calculations provide only the dipole moments or polarizability tensors, rather than the final spectrum itself. The challenge becomes greater when one goes beyond the static regime to compute vibrational spectra, such as when performing molecular dynamics simulations to capture anharmonic modes or solvent effects, or real-time time-dependent density functional theory (RT-TDDFT) simulations to access excited state dynamics, which is essential to calculate a RR spectrum.

Post-processing programs that are capable of handling these properties are often limited in scope and scattered across different software packages. In this work, we introduce Vibrant, a computational tool which bridges this gap by performing a series of post-processing procedures on dipole moment and polarizability data to generate both MD-based and static vibrational spectra. Figure&#160;1 provides an overview of the major capabilities of Vibrant.

The first panel (top left) emphasizes Vibrant's key feature, its ability to perform both static calculations based on the harmonic approximation and the post-processing of MD-based trajectories. Static spectra are computed via evaluating finite displacement data and performing normal mode analysis to obtain vibrational frequencies and normal mode coordinates, the latter is being used to produce the static vibrational spectra[@thomas2016theoretical]. In the MD-based approach, the normal mode derivatives are replaced by the Fourier transforms of the time-autocorrelation functions[@thomas2013computing; @gordon1965molecular].

The second panel (top right) features the main strengths of Vibrant, highlighting its capabilities for both static and MD-based vibrational spectroscopy, particularly IR, Raman and RR spectra. The latter is calculated by processing the time-dependent dipole trajectories, for example those generated from RT-TDDFT simulations. Although not depicted on the panel, Vibrant can also compute absorption spectra from time-dependent dipole data.

The third panel (bottom left) demonstrates Vibrant's capability to dissect MD-based spectra and compute the subspectra for user-defined molecular fragments. This feature facilitates the analysis and characterization of vibrational peaks, and is especially useful for identifying spectral features in guest-host systems or evaluating solvent contributions, where vibrational spectroscopy is of great importance[@ekin2024; guo2010cyclodextrin; sun2020covalent].

The last panel (bottom right) demonstrates the use of Padé approximants in absorption or RR spectra calculations, implemented in Vibrant through the integration of the GreenX library’s analytic continuation component[@leucke2025analytic]. In RT-TDDFT simulations, frequency resolution strongly depends on the trajectory length, and shorter trajectories often struggle with poor spectral resolution.[@bruner2016accelerated] Applying Padé interpolation is particularly useful for achieving spectral convergence from short RT-TDDFT trajectories[@bruner2016accelerated; @mattiat2018efficient]. As shown in the figure, Padé interpolation enables convergence and improved frequency resolution even for trajectories as short as 20 fs, yielding up to a fivefold reduction in overall computational cost.

 ![Overview of Vibrant's functionalities. The top left panel illustrates its ability to perform both static and MD-based vibrational spectroscopy calculations, showing examples of the position autocorrelation function and the second derivative of the energy with respect to finite displacements, which are used to compute vibrational frequencies frequencies from two different methods. The top right panel presents the major functionalities available in Vibrant. The bottom left panel highlights Vibrant's feature for calculating subspectra of user-defined molecular fragments. The bottom right panel demonstrates its integration with Padé approximants to achieve finer frequency resolution from time-dependent dipole data. More information about the theoretical background, calculation procedures and input file formats are available on the [Vibrant website](https://lstsgroup.github.io/vibrant/index.html).](flowchart.png){label="overview"}

# State of the field

Currently, there are only a limited number of computational programs available for calculating molecular dynamics (MD)-based and static vibrational spectra. A prominent example is TRAVIS[@brehm2020travis; brehm2011travis], which can process a wide range of properties obtained from MD simulations to compute various types of spectra and correlation functions, including power, IR, Raman[@thomas2013computing; @thomas2015voronoi], resonance Raman[@brehm2019computing], vibrational circular dichroism[@thomas2016classical], and vibrational optical activity[@brehm2017computing], in addition to distribution functions such as the radial distribution function. However, TRAVIS is mainly designed for liquid systems, and in order to compute aforementioned spectra, it processes either Wannier centers, or uses the electronic density to perform Voronoi integration[@thomas2015voronoi]. Other than TRAVIS, there are Python-based scripts, for example, those that process CP2K[@kuhne2020cp2k] polarizability tensors[@mattiat2018efficient] to generate static Raman spectra[@beat_hubmann_2020_4026342], or the scripts distributed with FHI-aims[@blum2009ab; abbott2025roadmap], which can read FHI-aims output files to perform normal mode analysis or compute static IR and Raman spectra. Another widely used tool is Phonopy[@phonopy-phono3py-JPCM; @phonopy-phono3py-JPSJ], which is highly effective for calculating phonon dispersion relations and vibrational densities of state; however, it is not designed for post-processing time-dependent quantities such as dipole moments or polarizabilities that come from MD or RT-TDDFT simulations.

To the best of our knowledge, there is currently no computational tool that offers the flexibility to post-process both Wannier centers and external static, time-dependent or field-induced Berry phase dipole moments[@vanderbilt2018berry] and density functional perturbation theory (DFPT) polarizability tensors[@gonze1989density; @gonze1995adiabatic], derived from either normal mode analysis, MD simulations or RT-TDDFT calculations.

Vibrant provides a flexible and efficient environment for post-processing different types of properties and computing static and MD-based vibrational spectra for gas-phase, liquid, and solid-state systems. Furthermore, it enables the computation of subspectra for user-defined molecular fragments, a feature which is not available in other codes. Another significant aspect of Vibrant is that it features the implementation of Padé approximants and is integrated with OpenMP parallelization. This feature promotes the computation of a well-resolved resonance Raman or absorption spectra even from relatively short RT-TDDFT trajectories.

# Acknowledgements

The Emmy Noether Programme of the German Research Foundation (Project No. 453275048) is gratefully acknowledged for funding. The authors also acknowledge the German Research Foundation for providing support within CRC 1415 (Chemistry of Synthetic Two-Dimensional Materials, project number 417590517) and RTG 2861 (Research Training Group, project number 491865171). The ZIH of the TU Dresden, the Jülich Supercomputing Centre, the Paderborn Center for Parallel Computing (PC2) as well as the Finnish CSC - IT Center for Science are acknowledged with appreciation for providing computational resources.

# References
