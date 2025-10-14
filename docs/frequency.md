# Frequencies

Vibrant can calculate the frequencies either statically from normal mode analysis or via processing MD trajectories. It can also produce a Molden output upon request for the visualization of the vibrational modes.

## a) Normal mode analysis

Vibrant can carry out a normal mode analysis using the atomic forces. In a polyatomic molecule, the potential energy surface can be approximated by a Taylor series expansion around minimum:

$$ 
V = V_0 + \sum_i^{3N} \left(\frac{\partial V}{\partial u_i}\right)_0 u_i + \frac{1}{2}\sum_{i,j}^{3N} \left(\frac{\partial^2 V}{\partial u_i u_j}\right)_0 u_i u_j + \frac{1}{6} \sum_{i,j,k}^{3N} \left(\frac{\partial^3 V}{\partial u_i u_j u_k}\right)_0 u_i u_j u_k + \cdots
$$

where ${u_i u_j}$ are the mass-weighted Cartesian displacements from the equilibrium position and $N$ is the number of atoms. The first two terms, corresponding to the energy at the current geometry ($V_0$) and the first derivative of the energy along the coordinate $u_i$ ($\frac{\partial V}{\partial u_i}$) vanish. This leaves us with the third term, that is also called Hessian ($ H_{jk} = \left(\frac{\partial^2 V}{\partial u_i u_j}\right)_0$).

The term $\frac{\partial V}{\partial u_i}$ corresponds to the atomic forces. Hessian matrix can therefore be constructed from the derivatives of atomic forces along the coordinates for each atom $N$ which are slightly displaced in $\pm x$, $\pm y$ and $\pm z$ directions, yielding $6N$ displaced structures. Diagonalizing the $6N\times6N$ mass-weighted Hessian matrix ($\frac{H_{ij}}{\sqrt{m_{i} m_{j}}}$) yields $3N$ eigenvalues $\lambda$ and $3N$ eigenvectors. Frequencies are then obtained from:

$$
\omega_k =\sqrt{\lambda_k}  \qquad \nu_k = \frac{\omega_k}{2\pi}
$$

where $\omega$ is the angular frequency  (measured in radians per second), and $\nu$ is the ordinary frequency measured in hertz. Vibrant per default prints only the vibrational frequencies, eliminating the first 6 (or 5 for linear molecules) frequencies which belong to the rotational and translational modes. The unit of the frequencies are cm $^{-1}$.

The eigenvectors $\mathbf{\tilde{u}}_{p}$ provide the amplitude of the normal mode vibrations, and eventually provide the normal mode coordinates $Q_p=\mathbf{\tilde{u}}_p^T\mathbf{u}_p$. The normal mode coordinates are also printed by Vibrant. Their mass-weighted versions can later be used in computing static IR (see Section [IR Spectra](IR_spec.md)) or Raman intensities. 

For each displaced structure, the user must provide the atomic forces appended in a single file (see File Formats for more details). Hessian diagonalization in Vibrant is performed using the [LAPACK](https://www.netlib.org/lapack/) library.

An example input section for performing the normal mode analysis is shown below:

```bash
&global
 spectra NMA
&end global
&system
 filename {$optimized_geometry_file_name}
&end system
&static
  &hessian
   force_file {$force_file_name}
  &end hessian
  displacement 0.001
&end static
```

where the `hessian` section revokes the hessian diagonalization, and the `displacement` refers to the displacement between the shifted geometries (given in Angstrom).

## b) Power spectrum

The MD-based equivalent of obtaining all frequencies irrespective of IR or Raman selection rules is to compute the power spectrum. The power spectrum can be computed from the autocorrelation functions of the time derivatives of particle velocities, which are usually available in an MD simulation. Alternatively, particle velocities can also be obtained from the derivatives of the position vectors with respect to the time step. This means that simply providing an MD position trajectory would be enough to generate a Power spectrum. Applying fast Fourier transformation using the [FFTW](https://www.fftw.org/) library to the sum of velocity autocorrelation functions in a system gives the power spectrum for that system, as shown below with the [equation](https://link.springer.com/book/10.1007/978-3-319-49628-3):

$$
P(\omega ) =  \frac{2c}{3k_{B}}\int_{-\infty }^{\infty } \left< \dot{r}(\tau )\dot{r}(t+\tau )\right>_{\tau }e^{-i\omega t}dt
$$

where $\omega$ is the angular frequency, $c$ is the speed of light, $\dot{r}$ is the time derivative of the position and $\tau$ and $t$ stand for the times for two different snapshots. The final power intensities are given in K.cm.

Power spectrum calculation in Vibrant can be revoked by adding the section:

```bash
...
&global
 spectra P
 ...
&end global
&system
 type_traj {$type_traj}
 filename {$trajectory_file_name}
 mass_weighting {$mass_weighting_flag}
&end system
&md
...
&end md
``` 

where `$type_traj` can be specified by user as `pos` or `vel` standing for "positions" or "velocities". Upon request, Vibrant can compute the mass-weighted power spectrum as well, which is controlled by the keyword `mass_weighting`. `$mass_weighting_flag` can be specified as `y` or `n` standing for "yes" or "no".
