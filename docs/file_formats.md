# File Formats and Example Scripts

`vibrant` processes multiple data files, such as the dipole moments, polarizabilities, and forces, in order to compute the desired vibrational spectrum. In principle, these data can be obtained from any quantum chemistry package, provided that the appended data files follow the input formats expected by `vibrant`, which are described in this section.

More information on all available keywords can be found on [Keyword Glossary](Keyword_Glossary.rst) and all complete example input files are available on [Examples](Examples.md).

## Example scripts

Although the _ab initio_ data can, in principle, be obtained using various quantum chemistry packages, most of our DFT calculations are performed with CP2K or FHI-aims. Therefore, we provide a set of utility Bash scripts that can read selected CP2K and FHI-aims output files and append the required data in the appropriate format for `vibrant`. These scripts are intended both to assist users of these packages and to serve as example starting points for users working with other electronic structure codes. The folder containing the example scripts and data files can be downloaded [here](_static/example_scripts_vibrant.tar.gz). The folder contains three directories:

1. `static_calculations` → Contains four Bash scripts for static spectrum calculations: one for appending CP2K forces, one for FHI-aims forces, one for CP2K dipole moments, and one for CP2K DFPT polarizabilities.  Additionally, directories containing CP2K and FHI-aims outputs for finite-displaced structures of the water molecule are provided to demonstrate how the data can be correctly appended by looping over these structures. In practice, these directories would also contain the corresponding displaced geometries required to perform the DFT calculations. Such directory structures can either be created manually by the user or generated automatically, for example, using the `get_vibrations.py` utility of [FHI-aims](https://fhi-aims.org/uploads/manual/Ch4/S6.html), prior to the DFT calculations. 
2. `md_append` → Contains two subdirectories that provide Bash scripts for reading MD-based Berry-phase dipole moments and DFPT polarizability data from CP2K outputs, along with example CP2K output files.
3. `rt-tddft_dipoles` → Contains a Bash script for reading and appending dipole moments from a CP2K RT-TDDFT output file. The provided script loops over MD frames; however, it can easily be adapted to loop over displaced structures for static calculations, following a similar structure to that described in Step 1 in this list. An example CP2K RT-TDDFT output file is also provided.

## Formatting of the input data files

### 1) Coordinates of atomic positions (.xyz files)

Static spectral calculations that are based on the normal mode analysis require the cartesian coordinates (in Å) of each atom in the system to be provided in a standard `.xyz` [format](https://en.wikipedia.org/wiki/XYZ_file_format).

Similarly, the power spectrum calculations requires that the user provides the atomic position (in Å) or velocity vectors (in bohr/au_time, where au_time is the time given in a.u.) for each MD snapshot in the form of an `.xyz` file format:

```bash
<total_number_of_atoms> #first MD frame
                       #leave 1 blank line
<element_symbol> <x_coord> <y_coord> <z_coord>
...
...
<total_number_of_atoms> #second MD frame
                       #leave 1 blank line
...

```

An example file can be downloaded from [here](_static/force.data).

### 2) Files for normal mode analysis

#### 2.1 Atomic forces

To carry out a normal mode analysis, the user must provide the atomic force vectors in a.u. for each $6N$ displaced structure (for more information, see Section [Frequencies](frequency.md#a-normal-mode-analysis). The file should have the format (ignoring the comment lines):

```bash
#start with forces coming from Atom 1 displaced in +x direction
<force_x <force_y> <force_z> #first atom
<force_x> <force_y> <force_z> #second atom
...
<force_x> <force_y> <force_z> #last atom
#forces coming from Atom 1 displaced in +y direction
<force_x> <force_y> <force_z> #first atom
...
#forces coming from Atom 1 displaced in +z direction
<force_x> <force_y> <force_z> #first atom
...
#forces coming from Atom 2 displaced in +x direction
<force_x> <force_y> <force_z> #first atom
...
... 
#finish looping over all atoms displaces in +x, +y, +z directions.
#continue with forces coming from Atom 1 displaced in -x direction
<force_x> <force_y> <force_z> #first atom
...
...
#finish with the last atom displaced in -z direction
<force_x> <force_y> <force_z> #first atom
...
<force_x> <force_y> <force_z> #last atom
```

As an example, [this bash script](_static/force.sh) enters the directories of each displaced structure and for example, appends the [CP2K force files](_static/force.xyz) accordingly.

#### 2.2 Normal mode frequencies and coordinates

The user can provide the normal mode frequencies and mass-weighted eigenvectors instead of forces. In that case, the normal mode frequencies should be given as:

```bash
<frequency> #for the first normal mode
<frequency> #for the second normal mode
...
<frequency> #for the last normal mode
```

An example file can be downloaded from [here](_static/normal_mode_freq.dat).

And the mass-weighted eigenvectors should be given as the atomic coordinates for each vibrational mode, such as:

```bash
<x_coord> <y_coord> <z_coord> # coordinates of Atom 1 for normal mode 1
<x_coord> <y_coord> <z_coord> # coordinates of Atom 2 for normal mode 1
...
<x_coord> <y_coord> <z_coord> # coordinates of the last atom for normal mode 1
<x_coord> <y_coord> <z_coord> # coordinates of Atom 1 for normal mode 2
...
...
<x_coord> <y_coord> <z_coord> # coordinates of the last atom for the last normal mode
```

An example file can be downloaded from [here](_static/normal_mode_displ.dat).

### 3) Berry-phase dipole moments

#### 3.1 MD-based spectra

The Berry phase dipole moment vectors for the MD-based spectra should be given in Debye in the `.xyz` format, combined for all MD snapshots in a single file such as:

```bash
1 # because there is only 1 dipole vector for the whole system
                       #leave 1 blank line
<dummy_value> <dipole_x> <dipole_y> <dipole_z> #first MD frame
1 
                       #leave 1 blank line
<dummy_value> <dipole_x> <dipole_y> <dipole_z> #second MD frame
...

                       #leave 1 blank line
<dummy_value> <dipole_x> <dipole_y> <dipole_z> #last MD frame
```

An example file can be downloaded from [here](_static/berry_dipole.xyz).

#### 3.2 Static spectra based on normal modes

The Berry phase dipole moment vectors for a static calculation should be given in Debye for each $6N$ displaced structure (for more information, see Section [Frequencies](frequency.md#a-normal-mode-analysis) combined in a single file, such as (ignoring the comment lines):

```bash
#dipoles coming from Atom 1 displaced in +x direction
1 # because there is only 1 dipole vector for the whole system
                       #leave 1 blank line
<dummy_value> <dipole_x> <dipole_y> <dipole_z> 
#dipoles coming from Atom 1 displaced in +y direction
1 
                       
<dummy_value> <dipole_x> <dipole_y> <dipole_z> 
#dipoles coming from Atom 1 displaced in +z direction
1 
                       
<dummy_value> <dipole_x> <dipole_y> <dipole_z> 
#dipoles coming from Atom 2 displaced in +x direction
1 
                       
<dummy_value> <dipole_x> <dipole_y> <dipole_z> 
...
...
#dipoles coming from Atom 2 displaced in +x direction
1 
                       
<dummy_value> <dipole_x> <dipole_y> <dipole_z> 
...
#dipoles coming from Atom 1 displaced in -x direction
1 
                       
<dummy_value> <dipole_x> <dipole_y> <dipole_z> 
...
#dipoles coming from the last atom displaced in -z direction
1 
                       
<dummy_value> <dipole_x> <dipole_y> <dipole_z> 
```

An example file can be downloaded from [here](_static/dipoles_water.xyz).

### 4) Wannier centers

If any MD-based spectrum is being calculated from Wannier centers, then the `wannier.xyz` for the whole system (atoms + Wannier centers labeled as "X") should be given for each MD frame in the `.xyz` file format. 

An example could be seen [here](_static/COF-1_solv_wannier_free.xyz).

### 5) DFPT polarizabilities

#### 5.1 MD-based spectra

The DFPT polarizabilities for the MD-based spectra should be given in Å $^3$ in the `.xyz` format. There should be three separate files in total, such that the first file contains $\alpha_{xx}$, $\alpha_{xy}$, $\alpha_{xz}$; the second file contains $\alpha_{yx}$, $\alpha_{yy}$, $\alpha_{yz}$; and the third file contains $\alpha_{zx}$, $\alpha_{zy}$, $\alpha_{zz}$. Each file should have the format:


```bash
1 # because there is only 1 polarizability tensor for the whole system
                       #leave 1 blank line
<dummy_value> <alpha_ix> <alpha_iy> <alpha_iz> #first MD frame
1 
                       #leave 1 blank line
<dummy_value> <alpha_ix> <alpha_iy> <alpha_iz> #second MD frame
...
1 
                       #leave 1 blank line
<dummy_value> <alpha_ix> <alpha_iy> <alpha_iz> #last MD frame
```

where $i$ is $x$, $y$, or $z$.

An example file can be downloaded from [here](_static/alpha.xyz).

#### 5.2 Static spectra based on normal modes

The DFPT polarizabilities for a static calculation should be given in Å $^3$ for each $6N$ displaced structure (for more information, see Section [Frequencies](frequency.md#a-normal-mode-analysis) combined in a single file, such as (ignoring the comment lines):

```bash
#leave 6 blank lines
#Atom 1 displaced in +x direction
<dummy_value> <alpha_xx> <alpha_yy> <alpha_zz> 
<dummy_value> <alpha_xy> <alpha_xz> <alpha_yz> 
<dummy_value> <alpha_yx> <alpha_zx> <alpha_zy> 
#Atom 1 displaced in +y direction
#leave 6 blank lines
<dummy_value> <alpha_xx> <alpha_yy> <alpha_zz> 
<dummy_value> <alpha_xy> <alpha_xz> <alpha_yz> 
<dummy_value> <alpha_yx> <alpha_zx> <alpha_zy> 
#Atom 1 displaced in +z direction
#leave 6 blank lines
<dummy_value> <alpha_xx> <alpha_yy> <alpha_zz> 
<dummy_value> <alpha_xy> <alpha_xz> <alpha_yz> 
<dummy_value> <alpha_yx> <alpha_zx> <alpha_zy> 
...
#Atom 2 displaced in +x direction
#leave 6 blank lines
<dummy_value> <alpha_xx> <alpha_yy> <alpha_zz> 
<dummy_value> <alpha_xy> <alpha_xz> <alpha_yz> 
<dummy_value> <alpha_yx> <alpha_zx> <alpha_zy> 
...
#Atom 1 displaced in -x direction
#leave 6 blank lines
<dummy_value> <alpha_xx> <alpha_yy> <alpha_zz> 
<dummy_value> <alpha_xy> <alpha_xz> <alpha_yz> 
<dummy_value> <alpha_yx> <alpha_zx> <alpha_zy> 
...
#Last atom displaced in -z direction
#leave 6 blank lines
<dummy_value> <alpha_xx> <alpha_yy> <alpha_zz> 
<dummy_value> <alpha_xy> <alpha_xz> <alpha_yz> 
<dummy_value> <alpha_yx> <alpha_zx> <alpha_zy> 
```

The reason for leaving six blank spaces is because, for example, the user can directly append the [CP2K DFPT polarizabilities](https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/PROPERTIES/LINRES/POLAR.html#CP2K_INPUT.FORCE_EVAL.PROPERTIES.LINRES.POLAR.DO_RAMAN). An example is shown below:

```bash

POLARIZABILITY TENSOR (atomic units):
xx,yy,zz        8.25528        6.66403        5.31055
xy,xz,yz       -0.00001        0.00194       -0.00116
yx,zx,zy       -0.00000        0.00174       -0.00112
POLARIZABILITY TENSOR (Angstrom^3):
xx,yy,zz        1.22331        0.98751        0.78694
xy,xz,yz       -0.00000        0.00029       -0.00017
yx,zx,zy       -0.00000        0.00026       -0.00017

```

`vibrant` only reads the polarizabilities in Å $^3$, ignoring the a.u. block.

An example file can be downloaded from [here](_static/polarizabilities.dat).

### 6) Time-dependent Berry phase dipole moments coming from RT-TDDFT

#### 6.1 MD-based spectra

The time-dependent Berry phase dipole moment vectors for the MD-based resonance Raman (RR) or absorption spectra should be given in Debye in the `.xyz` format, combined for all MD snapshots in a single file such as:

```bash
<dummy_value> #first MD frame
                       #leave 1 blank line
<dummy_value> <dipole_x> <dipole_y> <dipole_z> #unperturbed dipoles
<dummy_value> <dipole_x> <dipole_y> <dipole_z> #RT-TDDFT step 1
<dummy_value> <dipole_x> <dipole_y> <dipole_z> #RT-TDDFT step 2
...
<dummy_value> <dipole_x> <dipole_y> <dipole_z> #last RT-TDDFT step
<dummy_value> #second MD frame
                       #leave 1 blank line
<dummy_value> <dipole_x> <dipole_y> <dipole_z> #unperturbed dipoles
<dummy_value> <dipole_x> <dipole_y> <dipole_z> #RT-TDDFT step 1
...
...
<dummy_value> #last MD frame
                       #leave 1 blank line
<dummy_value> <dipole_x> <dipole_y> <dipole_z> #unperturbed dipoles
<dummy_value> <dipole_x> <dipole_y> <dipole_z> #RT-TDDFT step 1
...
<dummy_value> <dipole_x> <dipole_y> <dipole_z> #last RT-TDDFT step
```

Here the unperturbed dipole moments refer to $\mu^{0}_{i}$ described in Section [Resonance Raman (RR) and absorption spectra](Resonance_Raman_spec.md). Three files should be given in total, each containing the unperturbed dipole moments + time-dependent dipole moments perturbed under electric field applied in $x$, $y$ or $z$ direction.

An example file can be downloaded from [here](_static/RTP_dipoles_MD.xyz).

#### 6.2 Static spectra based on normal modes

The time-dependent Berry phase dipole moment vectors for a static RR or absorption spectrum calculation should be given in Debye for each $6N$ displaced structure (for more information, see Section [Frequencies](frequency.md#a-normal-mode-analysis) combined in a single file, such as (ignoring the comment lines):

```bash
#dipoles coming from Atom 1 displaced in +x direction
<dummy_value>
                       #leave 1 blank line
<dummy_value> <dipole_x> <dipole_y> <dipole_z> #unperturbed dipoles
<dummy_value> <dipole_x> <dipole_y> <dipole_z> #RT-TDDFT step 1
<dummy_value> <dipole_x> <dipole_y> <dipole_z> #RT-TDDFT step 2
...
<dummy_value> <dipole_x> <dipole_y> <dipole_z> #last RT-TDDFT step
#dipoles coming from Atom 1 displaced in +y direction
<dummy_value> 
                       #leave 1 blank line
<dummy_value> <dipole_x> <dipole_y> <dipole_z> #unperturbed dipoles
<dummy_value> <dipole_x> <dipole_y> <dipole_z> #RT-TDDFT step 1
...
#dipoles coming from Atom 1 displaced in +z direction
<dummy_value> 
                       #leave 1 blank line
<dummy_value> <dipole_x> <dipole_y> <dipole_z> #unperturbed dipoles
<dummy_value> <dipole_x> <dipole_y> <dipole_z> #RT-TDDFT step 1
...
#dipoles coming from Atom 2 displaced in +x direction
<dummy_value>
                       #leave 1 blank line
<dummy_value> <dipole_x> <dipole_y> <dipole_z> #unperturbed dipoles
<dummy_value> <dipole_x> <dipole_y> <dipole_z> #RT-TDDFT step 1
...
...
#dipoles coming from Atom 1 displaced in -x direction
<dummy_value>
                       #leave 1 blank line
<dummy_value> <dipole_x> <dipole_y> <dipole_z> #unperturbed dipoles
<dummy_value> <dipole_x> <dipole_y> <dipole_z> #RT-TDDFT step 1
...
...
#dipoles coming from the last atom displaced in -z direction
<dummy_value>
                       #leave 1 blank line
<dummy_value> <dipole_x> <dipole_y> <dipole_z> #unperturbed dipoles
<dummy_value> <dipole_x> <dipole_y> <dipole_z> #RT-TDDFT step 1
...
<dummy_value> <dipole_x> <dipole_y> <dipole_z> #last RT-TDDFT step
```

Similarly to the MD-based spectra, here the unperturbed dipole moments refer to $\mu^{0}_{i}$ described in Section [Resonance Raman (RR) and absorption spectra](Resonance_Raman_spec.md). Three files should be given in total, each containing the unperturbed dipole moments + time-dependent dipole moments perturbed under electric field applied in $x$, $y$ or $z$ direction.

 ```{note}
If the user only wants to calculate an absorption spectrum, then it is sufficient to provide the time-dependent dipole moments of only one geometry.

An example file for only one geometry can be downloaded from [here](_static/RTP_dipoles_static.xyz).
 ```
