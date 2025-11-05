# Frequencies

`vibrant` can calculate the frequencies either statically from normal mode analysis or via processing MD trajectories. It can also produce a Molden output upon request for the visualization of the vibrational modes.

## a) Normal mode analysis

`vibrant` can carry out a normal mode analysis using the atomic forces. In a polyatomic molecule, the potential energy surface can be approximated by a Taylor series expansion around minimum:

$$
V = V_{0} + \sum_{i}^{3N} \left(\frac{\partial V}{\partial u_{i}}\right)_ {0} u_{i} + \frac{1}{2}\sum_{i,j}^{3N} \left(\frac{\partial^{2} V}{\partial u_{i} u_{j}}\right)_ {0} u_{i} u_{j} + \frac{1}{6} \sum_{i,j,k}^{3N} \left(\frac{\partial^{3} V}{\partial u_{i} u_{j} u_{k}}\right)_ {0} u_{i} u_{j} u_{k} + \cdots
$$

where ${u_{i} u_{j}}$ are the mass-weighted Cartesian displacements from the equilibrium position and $N$ is the number of atoms. The first two terms, corresponding to the energy at the current geometry ($V_{0}$) and the first derivative of the energy along the coordinate $u_{i}$ ($\frac{\partial V}{\partial u_{i}}$) vanish. This leaves us with the third term, that is also called Hessian ($ H_{jk} = \left(\frac{\partial^{2} V}{\partial u_{i} u_{j}}\right)_ {0}$).

The term $-\frac{\partial V}{\partial u_{i}}$ corresponds to the atomic forces. The Hessian matrix can therefore be constructed from the derivatives of atomic forces along the coordinates for each atom $N$ which are slightly displaced in $\pm x$, $\pm y$ and $\pm z$ directions, yielding $6N$ displaced structures. Diagonalizing the $6N\times6N$ mass-weighted Hessian matrix ($\frac{H_{ij}}{\sqrt{m_{i} m_{j}}}$) yields $3N$ eigenvalues $\lambda$ and $3N$ eigenvectors. Frequencies are then obtained from:

$$
\omega_{k} =\sqrt{\lambda_{k}}  \qquad \nu_{k} = \frac{\omega_{k}}{2\pi}
$$

where $\omega$ is the angular frequency  (measured in radians per second), and $\nu$ is the ordinary frequency measured in hertz. `vibrant` prints all the normal mode frequencies, meaning that the first 6 (or 5 if the molecule is linear) normal modes belong to the rotational and translational modes, and the rest are the vibrational modes. The unit of the frequencies are cm $^{-1}$.

The eigenvectors $\mathbf{\tilde{u}}_ {p}$ provide the amplitude of the normal mode vibrations, and eventually provide the normal mode coordinates $Q_{p} = \mathbf{\tilde{u}}_ {p}^{T} \mathbf {u}_ {p}$. The normal mode coordinates are also printed by `vibrant`. Their mass-weighted versions can later be used in computing static IR, Raman or resonance Raman intensities.

For each displaced structure, the user must provide the atomic forces appended in a single file (see [File Formats](file_formats.md#21-atomic-forces) for more details). Hessian diagonalization is performed using the [LAPACK](https://www.netlib.org/lapack/) library.

An example input section for performing the normal mode analysis is shown below:

```bash
&global
 spectra NMA
&end global
&system
 filename <optimized_geometry_file_name>
&end system
&static
  &hessian
   force_file <force_file_name>
  &end hessian
  displacement <amount_of_displacement>
  write_mol_file
&end static
```

where the `hessian` section invokes the Hessian diagonalization, and the `displacement` refers to the displacement between the shifted geometries (given in Å).

 ```{note}
 The `write_mol_file` keyword is optional and it executes the printing of a `<filename>.mol` file, which includes the optimized geometry, normal mode frequencies and normal mode coordinates. The `<filename>.mol` file can be opened with [MOLDEN](https://www.theochem.ru.nl/molden/) to visualize the normal modes.  
 ```

## b) Power spectrum

The MD-based equivalent of obtaining all frequencies irrespective of IR or Raman selection rules is to compute the power spectrum. The power spectrum can be computed from the autocorrelation functions of the time derivatives of particle velocities, which are usually available in an MD simulation. The general equation for an [autocorrelation function](https://academic.oup.com/book/27866) is as follows:

$$
\left \langle A(\tau) A(\tau + t) \right \rangle_{\tau}
= \int A(\tau) A(\tau + t)\, d\tau
\approx \frac{1}{\tau_{\textnormal{max}}} \sum_{\tau=1}^{\tau_{\textnormal{max}}} A(\tau) A(\tau + t)
$$

where $\tau$ and $t$ refer to the times for two different snapshots and $\tau_{\textnormal{max}}$ is the correlation depth, which is an input parameter (specified with the `correlation_depth` keyword in the `vibrant` input file).

Alternatively, particle velocities can also be obtained from the derivatives of the position vectors with respect to the time step. This means that simply providing an MD position trajectory would be enough to generate a Power spectrum. Applying fast Fourier transformation using the [FFTW](https://www.fftw.org/) library to the sum of velocity autocorrelation functions in a system gives the power spectrum for that system, as shown below with the [equation](https://link.springer.com/book/10.1007/978-3-319-49628-3):

$$
P(\omega) = \frac{2c}{3k_{B}} \int_{-\infty}^{\infty} \left< \dot{r}(\tau)\dot{r}(t+\tau)\right>_{\tau} e^{-i\omega t} dt
$$

where $\omega$ is the angular frequency, $c$ is the speed of light, $\dot{r}$ is the time derivative of the position and $\tau$ and $t$ stand for the times for two different snapshots. The time-derivative of any property $x$ is evaluated through central finite differences:

$$
\dot{x} = \frac{\delta x(t)}{\delta t} \approx \frac{x(t+\Delta t)-x(t-\Delta t)}{2\Delta t}
$$

where $\Delta t$ refers to the time step between the selected equidistant MD snapshots.

Power spectrum calculation in `vibrant` can be invoked by adding the section:

```bash
...
&global
 spectra P
 ...
&end global
&system
 type_traj <type_traj>
 filename <trajectory_file_name>
 mass_weighting <mass_weighting_flag>
&end system
&md
...
&end md
```

where `type_traj` can be specified by user as `pos` or `vel` standing for "positions" or "velocities". Upon request, `vibrant` can compute the mass-weighted power spectrum as well, which is controlled by the keyword `mass_weighting`, which be specified as `y` or `n` standing for "yes" or "no".

The final power intensities are given in K cm .

### Spectral refinement

`vibrant` applies several procedures to improve the frequency resolution of MD-based spectra. It applies [finite difference correction](https://www.sciencedirect.com/science/article/pii/S0377042700003484) via dividing the intensities by the sinc function $\left(\frac{\sin(\tilde{\nu} \Delta t)}{\tilde{\nu} t} \right )^{2}$. Additionally, in order to increase the spectral resolution, `vibrant` applies ["data mirroring"](https://link.springer.com/book/10.1007/978-3-319-49628-3?utm_medium=referral&utm_source=google_books&utm_campaign=3_pier05_buy_print&utm_content=en_08082017) to the autocorrelation functions. This is done by appending the data set in reverse order to the original data set. `vibrant` also employs [Hann Window](https://link.springer.com/book/10.1007/978-3-319-49628-3?utm_medium=referral&utm_source=google_books&utm_campaign=3_pier05_buy_print&utm_content=en_08082017) function to further improve the resolution. This is achieved by multiplying the autocorrelation data with $\displaystyle\cos^{2}[\pi/(2\tau_{\textnormal{max}}-2)]$ before performing the Fourier transformations. $\tau_{\text{max}}$ corresponds to the correlation depth.


More information on the all available keywords can be found on [Keyword Glossary](Keyword_Glossary.rst) and all complete example input files are available on [Examples](Examples.md).
