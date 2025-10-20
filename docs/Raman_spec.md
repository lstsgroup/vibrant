# Raman spectra

The Raman spectrum calculation in Vibrant is based either on the post-processing of polarizability tensors, for example calculated by the density functional theory (DFPT), of through the induced dipole moments. The latter can process either Berry phase dipole moments or Wannier centers under applied electric field. In principle, Vibrant can compute both static or molecular-dynamics (MD)-based Raman spectra.

Vibrant can calculate the Raman intensities for a range of incident laser wavelengths $\tilde{\nu}_{in}$ specified in the input file via the `laser_in` keyword in the `raman` section. Every Raman spectrum calculated for a different laser wavenumber are appended as additional columns.

More details on the differences of static and MD-based Raman spectra are provided in the following sections.

## a) Static Raman intensities

The static Raman calculation in Vibrant relies on the normal mode analysis within the harmonic approximation. The static Raman intensities are computed from the derivatives of the polarizability tensors $\boldsymbol{\alpha}$ with respect to the mass-weighted normal coordinates $Q_{p}$. Following [Placzek's polarization theory](https://onlinelibrary.wiley.com/doi/book/10.1002/0470845767), the [Raman intensity](https://link.springer.com/book/10.1007/978-3-319-49628-3?utm_medium=referral&utm_source=google_books&utm_campaign=3_pier05_buy_print&utm_content=en_08082017) for any normal mode $p$ is given by:

$$
I_{\textnormal{Raman}}(\tilde{\nu}_p)=\frac{h}{8\epsilon_{0}^{2}c} \frac{(\tilde{\nu}_{in}-\tilde{\nu}_{p})^{4}}{\tilde{\nu}_{p}} \frac{1}{1-\exp\left(-\frac{hc\tilde{\nu}_{p}}{k_{B}T}\right)}\frac{X\beta^{2}_{p}+Y\gamma^{2}_{p} }{45} 
$$

where $\tilde{\nu}_{in}$ denotes the wavenumber of the incident light, $\varepsilon_{0}$ is the vacuum permittivity, $c$ is the speed of light,  
$h$ is Planck’s constant, $k_{B}$ is the Boltzmann constant and $T$ is the temperature which is specified by the user with the keyword `temperature` in `global` section.  $\tilde{\nu}_p$ is the wavenumber defined as $\tilde{\nu}_p=\nu_p/c$. $X$ and $Y$ are functions of the scattering geometry and polarization and will be discussed below. The isotropic contribution $\beta^{2}_{p}$ is obtained from the derivatives of the diagonal components of the polarizability tensor $\boldsymbol{\alpha}$ with respect to $Q_{p}$:

$$
\beta^{2}_{p}=\frac{1}{9}\left ( \frac{\partial \alpha _{xx}}{\partial Q_{p}}+\frac{\partial\alpha _{yy} }{\partial Q_{p}} +\frac{\partial \alpha _{zz}}{\partial Q_{p}}\right )^{2}
$$

and the anisotropic contribution $\gamma^{2}_{p}$ contains also the derivatives of the off-diagonal elements of $\boldsymbol{\alpha}$:

$$
\gamma_{p}^{2}=
\frac{1}{2}\left (\frac{\partial \alpha_{xx}}{\partial Q_{p}}-\frac{\partial \alpha_{yy}}{\partial Q_{p}}\right )^{2}
+\frac{1}{2}\left (\frac{\partial \alpha_{yy}}{\partial Q_{p}}-\frac{\partial \alpha_{zz}}{\partial Q_{p}}\right )^{2}
+\frac{1}{2}\left (\frac{\partial \alpha_{zz}}{\partial Q_{p}}-\frac{\partial \alpha_{xx}}{\partial Q_{p}}\right )^{2}
+3\left (\frac{\partial \alpha_{xy}}{\partial Q_{p}}\right )^{2}
+3\left (\frac{\partial \alpha_{yz}}{\partial Q_{p}}\right )^{2}
+3\left (\frac{\partial \alpha_{zx}}{\partial Q_{p}}\right )^{2}
$$

[$X$ and $Y$](https://onlinelibrary.wiley.com/doi/book/10.1002/0470845767) can take the values listed below:

| Intensities | X | Y |
|---|---|---|
Orthogonal ($I_{\perp }(\tilde{\nu})$) | 0 | 3 |
Parallel ($I_{\parallel }(\tilde{\nu})$) | 45 | 4 | 
Unpolarized ($I(\tilde{\nu}) = I_{\perp }(\tilde{\nu}) + I_{\parallel }(\tilde{\nu})$) | 45 | 7 |

By default, Vibrant prints only the unpolarized (total) Raman intensities, but the other types of spectra, together with the depolarization ratio ($\rho (\tilde{\nu})=\frac{I_{\perp }(\tilde{\nu})}{I_{\parallel  }(\tilde{\nu})}$), can be printed out by setting the `spectra_verbosity` keyword to `high`.

For each displaced structure (see Section [Frequencies](frequency.md) for details) the user must provide the polarizability tensors appended in a single file (see [File Formats](file_formats.md) for more details.) For the derivatives along the normal modes, the user can either provide the normal mode coordinates or alternatively the forces for each displaced structure, together with the cartesian coordinates of the optimized geometry. An example input section for the static Raman section may look like:

```bash
&global
 spectra R
 temperature {$temperature_in_K}
 fwhm {$full_width_at_half_maximum}
&end global
...
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
 spectra R
 temperature {$temperature_in_K}
 fwhm {$full_width_at_half_maximum}
&end global
...
&static
 diag_hessian n
 normal_freq_file {$normal_mode_frequency_file}
 normal_displ_file {$normal_mode_displacement_file}
 displacement {$amount_of_displacement}
 write_mol_file
&end static
``` 

Vibrant applies Gaussian broadening to the final discrete set of frequencies and intensities. The `fwhm` keyword controls the full width at half maximum (FWHM) value in cm ${^{-1}}$ and if not specified, it is set to 5 cm ${^{-1}}$. 

 ```{note}
  The keyword `write_mol_file` is optional and it executes the printing of a `{$filename}.mol` file, which includes the optimized geometry, normal mode frequencies, normal mode coordinates and the non-broadened Raman intensities. The `{$filename}.mol` file can be opened with [MOLDEN](https://www.theochem.ru.nl/molden/) to visualize the normal modes alongside the Raman spectrum. If multiple incident laser frequencies are requested, Vibrant generates a separate `{$filename}.mol` file for each frequency.
```

The final static Raman intensities are reported in 10 $^{-30}$ cm $^2$ /sr.

```{warning}
 Currently, Vibrant does not support the use of induced dipoles for calculating static Raman spectra, as the finite-difference error can be significant when using the limited number of structures available in a static calculation compared to an MD trajectory. Therefore, static Raman spectra can only be computed from the provided polarizabilities.
 ```

## b) MD-based Raman intensities

In the MD-based Raman calculation, the polarizability derivatives along the normal mode coordinates $Q_{p}$ are replaced with the [time derivatives of the polarizability autocorrelation functions](https://pubs.rsc.org/en/content/articlehtml/2013/cp/c3cp44302g). We applied fast Fourier transformations to the time-domain polarizability autocorrelations using the [FFTW](https://www.fftw.org/) library to convert them into the frequency-domain. The [MD-based Raman intensities](https://link.springer.com/book/10.1007/978-3-319-49628-3?utm_medium=referral&utm_source=google_books&utm_campaign=3_pier05_buy_print&utm_content=en_08082017) are calculated as: 

$$
I_{\textnormal{Raman}} (\tilde{\nu}) = \frac{2h}{8\epsilon_{0}^{2}k_{B}T} \frac{(\tilde{\nu}_{in}-\tilde{\nu})^{4}}{\tilde{\nu}} \frac{1}{1-\exp\left(-\frac{hc\tilde{\nu}}{k_{B}T}\right)}\frac{X\delta^{2} (\tilde{\nu})+Y\epsilon^{2} (\tilde{\nu}) }{45} 
$$

where $h$ is the Planck’s constant, $\tilde{\nu}_{in}$ is the wavenumber of the incident laser, $\varepsilon_{0}$ is the vacuum permittivity, $k_{B}$ is the Boltzmann constant, $c$ is the speed of light and $T$ is the temperature which is specified by the user with the keyword `temperature` in `global` section. $X$ and $Y$ are functions of the scattering geometry and polarization as discussed in the section above. $\delta^{2} (\tilde{\nu})$ and $\epsilon^{2} (\tilde{\nu})$ refer to the [isotropic and anisotropic contributions](https://books.google.de/books/about/Quantum_Chemistry.html?id=zzxLTIljQB4C&redir_esc=y):

$$
\delta^{2} (\tilde{\nu}) = 
\int_{-\infty}^{\infty}
\left\langle 
\frac{\dot\alpha_{xx}(\tau) + \dot\alpha_{yy}(\tau) + \dot\alpha_{zz}(\tau)}{3}
\times
\frac{\dot\alpha_{xx}(\tau+t) + \dot\alpha_{yy}(\tau+t) + \dot\alpha_{zz}(\tau+t)}{3}
\right\rangle_{\tau}
\, \exp(-2\pi i c \tilde{\nu} t)\, dt
$$

where $\tau$ and $t$ refer to the times for two different snapshots, $\dot{\boldsymbol{\alpha}}$ denotes the time derivative of the polarizability tensor $\boldsymbol{\alpha}$. The anisotropic contribution $\epsilon^{2} (\tilde{\nu})$ is given by:

$$
\begin{aligned}
\epsilon^{2} (\tilde{\nu}) = {} &
\int_{-\infty}^{\infty} \Biggl[
\frac{1}{2}\left\langle
\left\{\dot\alpha_{xx}(\tau)-\dot\alpha_{yy}(\tau)\right\}
\left\{\dot\alpha_{xx}(\tau+t)-\dot\alpha_{yy}(\tau+t)\right\}
\right\rangle_{\tau} \\
& + \frac{1}{2}\left\langle
\left\{\dot\alpha_{yy}(\tau)-\dot\alpha_{zz}(\tau)\right\}
\left\{\dot\alpha_{yy}(\tau+t)-\dot\alpha_{zz}(\tau+t)\right\}
\right\rangle_{\tau} \\
& + \frac{1}{2}\left\langle
\left\{\dot\alpha_{zz}(\tau)-\dot\alpha_{xx}(\tau)\right\}
\left\{\dot\alpha_{zz}(\tau+t)-\dot\alpha_{xx}(\tau+t)\right\}
\right\rangle_{\tau} \\
& + \frac{3}{4}\left\langle
\left\{\dot\alpha_{xy}(\tau)+\dot\alpha_{yx}(\tau)\right\}
\left\{\dot\alpha_{xy}(\tau+t)+\dot\alpha_{yx}(\tau+t)\right\}
\right\rangle_{\tau} \\
& + \frac{3}{4}\left\langle
\left\{\dot\alpha_{yz}(\tau)+\dot\alpha_{zy}(\tau)\right\}
\left\{\dot\alpha_{yz}(\tau+t)+\dot\alpha_{zy}(\tau+t)\right\}
\right\rangle_{\tau} \\
& + \frac{3}{4}\left\langle
\left\{\dot\alpha_{zx}(\tau)+\dot\alpha_{xz}(\tau)\right\}
\left\{\dot\alpha_{zx}(\tau+t)+\dot\alpha_{xz}(\tau+t)\right\}
\right\rangle_{\tau}
\Biggr] \\
& \times \exp(-2\pi i c \tilde{\nu} t)\, dt
\end{aligned}
$$

The `md` section should be included for the MD-based Raman calculation:

```bash
&global
 spectra MD-R
 temperature {$temperature_in_K}
&end global
...
&md
 time_step {$MD_time_step}
 correlation_depth {$correlation_depth}
&end md
```

 ```{note}
  Vibrant applies the same processing to the final MD-based Raman intensities as given in the Section [IR Spectra](IR_spec.md), including the application of data mirroring and Hann Window function to the autocorrelation data and the application of the sinc function to the final intensities. 
```

 The final MD-based Raman intensities are reported in m $^2$ K cm 10 $^{-30}$ . The polarizability/dipole moment types that can be processed by Vibrant for computing Raman spectra are discussed in the next section.

## c) Different polarizability tensors

The MD-based Raman intensities can be computed either directly from the DFPT polarizabilities or from an induced dipole approach where Berry phase dipole moments or Wannier centers are processed. [CP2K](https://www.cp2k.org/) program package enables the [calculation of the polarizability tensors via DFPT](https://pubs.aip.org/aip/jcp/article/141/9/094503/193754), which can be processed by Vibrant. [DFPT](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.39.13120) analytically evaluates the response of a system to external electric or magnetic fields. The DFPT polarizabilities are obtained for the whole supercell, and are obtained directly from a single linear-response calculation, without any field adjustments. The DFPT polarizabilities are uniquely defined unlike the dipole moments, meaning no phase indeterminacy is involved. The DFPT polarizabilities can be specified in the Vibrant input file using the `type_dipole` keyword:

```bash
...
&dipoles
 type_dipole dfpt
 dip_x_file {$dipole_file_name_x_field}
 dip_y_file {$dipole_file_name_y_field}
 dip_z_file {$dipole_file_name_z_field}
&end dipoles 
```

Polarizabilities can be also calculated numerically through [induced dipole moments](https://pubs.rsc.org/en/content/articlehtml/2013/cp/c3cp44302g) $\boldsymbol{\mu}_{\text{ind}}$. The polarizability tensor $\boldsymbol{\alpha}$ is related to the change in the dipole moment under the applied external electric field $\mathbf{E}=(E_x, E_y, E_z)$ by:

$$
\boldsymbol{\mu}_{\textnormal{ind}} =\boldsymbol{\alpha} \mathbf{E} 
$$

For each selected MD snapshot, three single-point calculations must be performed with a finite periodic field applied along the $x$, $y$ or $z$ direction in addition to a field-free (unperturbed) calculation. The induced dipole moment is then calculated via forward finite differences as:

$$
\boldsymbol{\mu}_{\textnormal{ind}}^k = \boldsymbol{\mu}^{k}-\boldsymbol{\mu}^{0} \quad \textnormal{with} \quad k=x,y,z
$$

where $\boldsymbol{\mu}^{0}$ is the dipole moment vector obtained from the field-free calculations. The components $\alpha_{\alpha\beta}$ of the polarizability tensor are then given by $ \alpha_{\alpha,\beta=k} = \mu^k_{\textnormal{ind},\alpha}/E$ .

Dipoles can be computed using the Berry-phase or MLWF formalism, enabling consistent treatment of IR and Raman spectra and local Raman analyses. The field strength must be optimized to balance finite-difference accuracy and noise. More details on the Berry phase method and maximally localized Wannier functions (MLWFs) are provided in Section [IR Spectra](IR_spec.md). The induced dipole section must also include the file that contains the field-free dipole moments and the `field_strength` keyword:

```bash
...
&dipoles
 type_dipole berry #or `wannier`
 field_strength {$field_strength_in_au}
 dip_file {$dipole_file_name_field_free}
 dip_x_file {$dipole_file_name_x_field}
 dip_y_file {$dipole_file_name_y_field}
 dip_z_file {$dipole_file_name_z_field}
&end dipoles
```

More details on the format of the dipole moment files can be found at ...

```{warning}
 Similar to the MD-based IR spectra, the MD-based spectrum should contain the `cell` subsection alongside the `md` section if the induced dipole method is selected. For more details, please Section [IR Spectra](IR_spec.md).
```

```{note}
 The Wannier centers can be used to compute the dipole moment and consequently the polarizability of the whole supercell, however they can also be used to extract the spectra of user-specified molecular blocks or fragments, which is discussed in Section [Subspectra for MD-based calculations](fragments.md).
```

More information on the all available keywords can be found on Section .. and all complete example input files are available on ....
