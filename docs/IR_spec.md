# IR spectra

The infrared (IR) spectrum calculation in Vibrant is based on the post-processing of either Berry-phase dipole moments or Wannier centers. In principle, Vibrant can compute both static or molecular-dynamics (MD)-based IR spectra. More details are provided in the following sections. 

## a) Static IR intensities

The static IR calculation in Vibrant relies on the normal mode analysis within the harmonic approximation. The static IR intensities are computed from the derivatives of the dipole moments $\boldsymbol{\mu}$ with respect to the mass-weighted normal coordinates $Q_{p}$. The following equation is used to generate the final static IR spectrum [static IR spectrum](https://onlinelibrary.wiley.com/doi/full/10.1002/jcc.10089):

$$ 
I_{\textnormal{IR}}(\tilde{\nu}_p) = \frac{N_{A} }{12\varepsilon _{0}c^{2}}\left ( \frac{\partial \boldsymbol{\mu} }{\partial Q_{p}} \right )^{2}
$$

where $\varepsilon_{0}$ is the vacuum permittivity, $N_{A}$ is the Avogadro constant and $c$ is the speed of light. $\tilde{\nu_{p}}$ is the wavenumber defined as $\tilde{\nu_{p}}=\nu_{p}/c$.

For each displaced structure (see Section [Frequencies](frequency.md) for details) the user must provide the dipole moments appended in a single file (see File Formats for more details.) For the derivatives along the normal modes, the user can either provide the normal mode coordinates or alternatively the forces for each displaced structure, together with the cartesian coordinates of the optimized geometry. An example input section for the static section may look like:

```bash
&global
 spectra IR
 fwhm {$full_width_at_half_maximum}
&end global
&system
 filename {$optimized_geometry_file_name}
&end system
&static
 diag_hessian y
 &hessian
  force_file {$FORCE_FILE_NAME}
 &end hessian
 displacement {$amount_of_displacement}
&end static
``` 

or if the user want to skip the Hessian diagonalization and give the normal mode frequencies and eigenvectors externally, the input section becomes:

```bash
&global
 spectra IR
 fwhm {$full_width_at_half_maximum}
&end global
&system
 filename {$optimized_geometry_file_name}
&end system
&static
 diag_hessian n
 normal_freq_file {$normal_mode_frequency_file}
 normal_displ_file {$normal_mode_displacement_file}
 displacement {$amount_of_displacement}
&end static
``` 

Vibrant applies Gaussian broadening to the final discrete set of frequencies and intensities. The `fwhm` keyword controls the full width at half maximum (FWHM) value in cm ${^{-1}}$ and if not specified, it is set to 5 cm ${^{-1}}$. The final static IR intensities are reported in km/mol. 

## b) MD-based IR intensities

In the MD-based IR calculation, the dipole moment derivatives along the normal mode coordinates $Q_{p}$ are replaced with the [time derivatives of the dipole moment autocorrelation functions](https://pubs.rsc.org/en/content/articlehtml/2013/cp/c3cp44302g). We applied fast Fourier transformations to the time-domain dipole autocorrelations using the [FFTW](https://www.fftw.org/) library to convert them into the frequency-domain. The following equation is used to generate the final MD-based IR spectrum:

$$
I_{\textnormal{IR}} (\tilde{\nu}) = \frac{2N_{A}}{12\varepsilon _{0}ck_{B}T}\int_{-\infty }^{\infty }\left \langle \dot{\boldsymbol{\mu} }(\tau )\cdot \dot{\boldsymbol{\mu} }(\tau+t ) \right \rangle_{\tau}\exp(-2\pi c\tilde{\nu}t)dt 
$$

where $\varepsilon_{0}$ is the vacuum permittivity, $N_{A}$ is the Avogadro constant, $c$ is the speed of light, $k_{B}$ is the Boltzmann constant, $T$ is the temperature and $\tilde{\nu}$ is the wavenumber. The final MD-based IR intensities are reported in K.cm.km/mol.

The MD-based IR intensities can be computed either from the Berry phase dipole moments or Wannier centers. The latter can also be used to extract spectra of user-specified molecular blocks or fragments.

> ⚠️ **Warning:** The MD-based spectrum should contain the `cell` subsection alongside the `md` section. Because the dipole moment in periodic systems is [ill-defined](https://www.cambridge.org/de/universitypress/subjects/physics/condensed-matter-physics-nanoscience-and-mesoscopic-physics/berry-phases-electronic-structure-theory-electric-polarization-orbital-magnetization-and-topological-insulators?format=HB&isbn=9781107157651) (i.e. only defined modulo a polarization quantum), it may exhibit discontinuous jumps. Vibrant automatically detects and corrects these jumps by subtracting the integer multiples of polarization quantum, ensuring continuous dipole trajectories. `cell` subsection provides the cell parameters which are necessary to calculate the polarization quantum.

An example input section may look like:

```bash
&global
 spectra MD-IR
&end global
&system
 &cell
  box_x {$x_vector_of_cell}
  box_y {$y_vector_of_cell}
  box_z {$z_vector_of_cell}
  angle_alpha {$angle_between_x_and_z}
  angle_beta {$angle_between_y_and_z}
  angle_gamma {$angle_between_x_and_y}
 &end cell
&end system
&md
 time_step {$MD_time_step}
 correlation_depth {$correlation_depth}
&end md
``` 

### Processing of the final spectra

Vibrant applies several procedures to improve the spectral resolution. It applies [finite difference correction](https://www.sciencedirect.com/science/article/pii/S0377042700003484) via dividing the intensities by the sinc function $\left(\frac{\sin(\tilde{\nu} \Delta t)}{\tilde{\nu} t} \right )^{2}$. Additionally, in order to increase the spectral resolution, Vibrant applies ["data mirroring"](https://link.springer.com/book/10.1007/978-3-319-49628-3?utm_medium=referral&utm_source=google_books&utm_campaign=3_pier05_buy_print&utm_content=en_08082017) to the autocorrelation functions. This done by appending the data set in reverse order to the original data set. Vibrant also employs [Hann Window](https://link.springer.com/book/10.1007/978-3-319-49628-3?utm_medium=referral&utm_source=google_books&utm_campaign=3_pier05_buy_print&utm_content=en_08082017) function to further improve the resolution. This is achieved by multiplying the autocorrelation data with $\displaystyle\cos^{2}[\pi/(2\tau_{\textnormal{max}}-2)]$ before performing the Fourier transformations. $\tau_{\text{max}}$ corresponds to the correlation depth.

## c) Dipole moment types

As stated above, the MD-based IR intensities can be computed either from the Berry phase dipole moments (for the whole supercell) or Wannier centers. The [Berry phase dipole moment](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.47.1651) describes the overall polarization of a periodic system by tracking how the electronic wavefunctions change across the supercell, giving a consistent way to define dipoles in materials with periodic boundaries.

Wannier functions are localized orbitals obtained from a unitary transformation of Bloch states; [maximally localized Wannier functions](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.56.12847) (MLWFs) apply an optimization procedure to achieve the most localized and unique representation. Wannier centers are the cartesian coordinates of the MLWFs in real-space. The dipole moments can be computed through Wannier centers employing:

$$
\boldsymbol{\mu}^{\text{MLWF}}=-2e\sum_{i}^{}\mathbf{\tilde{r}}_{i}^{\mathrm{WC}} + e\sum_{i}Z_{i}\mathbf{\tilde{R}}_{i}
$$

where the distance vectors $\mathbf{\tilde{r}}_{i}^{\mathrm{WC}} = \mathbf{r}_i^{\mathrm{WC}}-\mathbf{r}_0$ and $\mathbf{\tilde{R}}_{i} =  \mathbf{R}_i-\mathbf{r}_0$ are implemented with periodic boundary conditions (PBCs) and $\mathbf{r}_0$ is a reference point, which is the center of mass in Vibrant. $\mathbf{R}_{i}$ and $Z_{i}$ denote the atomic positions and charges, respectively.

 > 💡 **Note:** The static IR intensities currently can only be computed from the Berry phase dipole moments provided for the whole system. 
 
 The selection of the dipole moment type is controlled with the keyword `type_dipole`. An example input section for the specification of Berry phase dipole moments may look like:

```bash 
&dipoles
 type_dipole berry
 dip_file {$DIPOLE_FILE_NAME}
&end dipoles
```
And for the Wannier centers, it would be:

```bash 
&dipoles
 type_dipole wannier
 dip_file {$DIPOLE_FILE_NAME}
&end dipoles
```

The Wannier centers can be used to compute the dipole moment of the whole supercell, however they can also be used to extract spectra of user-specified molecular blocks or fragments, which is discussed in Section [Subspectra for MD-based calculations](fragments.md).

More information on the all available keywords can be found on Section .. and all complete example input files are available on ....

