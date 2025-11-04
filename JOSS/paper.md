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

In this work, we present Vibrant, a vibrational analysis program written in Fortran 2008. It interfaces to the FFTW [@frigo1998fftw] and GreenX (GX-AC) [@leucke2025analytic] libraries and is parallelized using OpenMP for computational efficiency. Vibrant enables the generation of vibrational spectra for gases, liquids and materials through post-processing atomic positions, velocities, forces, dipole moments and polarizability tensors. Its functionalities include the computation of vibrational frequencies, as well as infrared (IR) and Raman intensities, using either static or molecular dynamics (MD)-based approaches. In addition, Vibrant enables the computation of resonance Raman (RR) and absorption spectra based on real-time propagation trajectories. Vibrant is distributed under the Apache 2.0 license and is available as open-source software on GitHub.

# Statement of need

Vibrational spectroscopy, particularly infrared (IR) and Raman techniques, provides valuable insight into the structure of gaseous, liquid, and solid materials. While in conventional Raman spectroscopy the excitation is non-resonant with electronic transitions, resonance Raman (RR) spectroscopy occurs when the excitation wavelength matches such a transition [@long2002raman]. The resulting enhancement of specific vibrational modes allows, for example, in-situ monitoring of reactants or products through characteristic marker bands during synthesis [@reichmayr2025].

Computational simulation of vibrational spectroscopy is often not straightforward. Most quantum chemistry packages that perform electronic structure calculations provide only the dipole moments or polarizability tensors, rather than the final spectrum itself. The challenge becomes greater when one goes beyond the static regime to compute vibrational spectra, such as when performing MD simulations to capture anharmonic modes or solvent effects, or real-time time-dependent density functional theory (RT-TDDFT) simulations to access excited state dynamics, relevant for RR spectroscopy.

Post-processing programs that can handle these properties are often limited in scope and scattered across different software packages. In this work, we introduce Vibrant, a computational tool which bridges this gap by performing a series of post-processing procedures on dipole moment and polarizability data to generate MD-based or static vibrational spectra. Our code can provide these spectra for gas-phase, liquids and a wide range of materials. Figure&#160;1 provides an overview of Vibrant's major functionalities.

The first panel (top-left) summarizes Vibrant’s frequency analysis capabilities, which include static calculations within the harmonic approximation and post-processing of MD trajectories. In the static approach, conventional normal mode analysis is used to obtain vibrational frequencies and normal mode coordinates, while in the MD-based approach, normal mode derivatives are replaced by Fourier transforms of time-autocorrelation functions, yielding so-called power spectra [@thomas2013computing; @gordon1965molecular].

The second panel (top-right) highlights Vibrant’s functionality for computing IR, non-resonant Raman and RR spectral intensities, which can be obtained from either static or MD-based calculations. For IR and non-resonant Raman spectroscopy, Vibrant can handle different types of dipole moments and polarizabilities. Dipole moments may be provided as Berry-phase dipoles or generated from Wannier centers, while polarizabilities can be taken from density functional perturbation theory (DFPT) or computed from induced dipole moments via a finite-difference approach. RR intensities are obtained by post-processing time-dependent (Berry-phase) dipole trajectories RT-TDDFT simulations. Although not shown in the panel, Vibrant also supports the computation of absorption spectra from RT-TDDFT dipole data. 

The third panel (bottom-left) displays Vibrant's ability to dissect MD-based spectra via computing the subspectra for user-defined molecular fragments. This feature facilitates the analysis of the spectra, and is especially useful for characterizing guest-host systems, evaluating solvent contributions and molecular materials [@bas2024robust; @guo2010cyclodextrin; @sun2020covalent].

The last panel (bottom-right) demonstrates the use of Padé approximants in absorption or RR spectra calculations, implemented in Vibrant through integrating the GreenX library’s analytic continuation component [@leucke2025analytic]. In RT-TDDFT simulations, frequency resolution strongly depends on the trajectory length, and shorter trajectories often struggle with poor spectral resolution [@bruner2016accelerated]. Applying Padé interpolation is particularly useful for achieving spectral convergence from short RT-TDDFT trajectories [@bruner2016accelerated; @mattiat2018efficient]. As demonstrated in the bottom-right panel, Padé interpolation enables convergence and improved frequency resolution even for trajectories as short as 20 fs, resulting in a fivefold reduction of overall computational cost.

 ![Overview of Vibrant's functionalities. The top-left panel illustrates its ability to perform both static and MD-based vibrational spectroscopy calculations, where the vibrational frequencies are calculated from the position autocorrelation functions or the energy derivatives along displaced coordinates. The top-right panel summarizes the available corresponding spectral intensities. The bottom-left panel highlights Vibrant's feature of calculating subspectra for user-defined molecular fragments. The bottom-right panel demonstrates its integration with Padé approximants to achieve finer frequency resolution from time-dependent dipole data. More information about the theoretical background and calculation procedures are available on the [Vibrant website](https://lstsgroup.github.io/vibrant/index.html).](flowchart.png){label="overview"}

# State of the field

Currently, there are only a limited number of computational programs available for calculating MD-based and static vibrational spectra. A prominent example is TRAVIS [@brehm2020travis; @brehm2011travis], which processes various properties obtained from MD simulations to compute different types of spectra, including power, IR, Raman [@thomas2013computing; @thomas2015voronoi], RR [@brehm2019computing], vibrational circular dichroism [@thomas2016classical], and vibrational optical activity [@brehm2017computing]. However, TRAVIS is mainly designed for liquids, computing molecular dipole moments by post-processes either Wannier centers, or using the electronic density to perform Voronoi integration [@thomas2015voronoi]. Other than TRAVIS, there are Python-based scripts, for example, those that process CP2K [@kuhne2020cp2k] polarizabilities [@mattiat2018efficient] to generate static  Raman spectra [@beat_hubmann_2020_4026342], or the scripts distributed with FHI-aims [@blum2009ab; @abbott2025roadmap], which read FHI-aims output files to perform normal mode analysis or compute static IR and Raman spectra. A popular tool for materials is Phonopy [@phonopy-phono3py-JPCM; @phonopy-phono3py-JPSJ], which can calculate phonon dispersion relations and vibrational densities of state; however, it is not designed for MD-based spectral analysis.

To the best of our knowledge, no existing computational tool provides the same level of flexibility as Vibrant, which supports a wide range of vibrational spectra, both static and MD-based methodologies, diverse dipole and polarizability types, and applicability across gaseous, liquid, and solid-state systems.

# Acknowledgements

The Emmy Noether Programme of the German Research Foundation (Project No. 453275048) is gratefully acknowledged for funding. The authors acknowledge the German Research Foundation for support within CRC 1415 (Chemistry of Synthetic Two-Dimensional Materials, Project No. 417590517) and RTG 2861 (Research Training Group, Project No. 491865171). The ZIH of the TU Dresden, the Jülich Supercomputing Centre, the Paderborn Center for Parallel Computing (PC2) as well as the Finnish CSC - IT Center for Science are acknowledged for providing computational resources.

# References
