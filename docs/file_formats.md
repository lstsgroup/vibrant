# File Formats

Vibrant processes multiple data files, such as dipole moments, polarizabilities, and forces, to compute the desired vibrational spectrum. These files must follow specific formats which we describe in this section.

More information on the all available keywords can be found on Section .. and all complete example input files are available on ....

## 1) Coordinates of atomic positions (.xyz files)

Static spectral calculations that are based on the normal mode analysis require the cartesian coordinates (in Å) of each atom in the system to be provided in a standard `.xyz` [format](https://en.wikipedia.org/wiki/XYZ_file_format).

Similarly, the power spectrum calculations requires that the user provides the atomic position or velocity vectors for each MD snapshot in the form of an `.xyz` file format, such as:

```bash
$total_number_of_atoms #first MD frame
                       #leave 1 blank line
$element_symbol $x_coord $y_coord $z_coor
...
...
$total_number_of_atoms #second MD frame
                       #leave 1 blank line
...

```

## 2) Files for normal mode analysis

### 2.1 Atomic forces

To carry out a normal mode analysis, the user must provide the atomic force vectors in a.u. for each $6N$ displaced structure (for more information, see Section [Frequencies](frequency.md). The file should have the format (ignoring the comment lines):

```bash
#start with forces coming from Atom 1 displaced in +x direction
$force_x $force_y $force_z #first atom
$force_x $force_y $force_z #second atom
...
$force_x $force_y $force_z #last atom
#forces coming from Atom 1 displaced in +y direction
$force_x $force_y $force_z #first atom
...
#forces coming from Atom 1 displaced in +z direction
$force_x $force_y $force_z #first atom
...
#forces coming from Atom 2 displaced in +x direction
$force_x $force_y $force_z #first atom
...
... 
#finish looping over all atoms displaces in +x, +y, +z directions.
#continue with forces coming from Atom 1 displaced in -x direction
$force_x $force_y $force_z #first atom
...
...
#finish with the last atom displaced in -z direction
$force_x $force_y $force_z #first atom
...
$force_x $force_y $force_z #last atom
```

As an example, [this bash script](_static/force.sh) enters directories of each displaced structure and appends the [CP2K force files](_static/force.xyz) accordingly.

### 2.2 Normal mode frequencies and coordinates

The user can provide the normal mode frequencies and coordinates instead of forces. In that case, the normal mode frequencies should be given as:

```bash
$frequency #for the first normal mode
$frequency #for the second normal mode
...
$frequency #for the last normal mode
```

And the normal mode coordinates should be given as the atomic cartesian coordinates (in. a.u.) for each vibrational mode, such as:

```bash
$x_coord $y_coord $z_coord # coordinates of Atom 1 for normal mode 1
$x_coord $y_coord $z_coord # coordinates of Atom 2 for normal mode 1
...
$x_coord $y_coord $z_coord # coordinates of the last atom for normal mode 1
$x_coord $y_coord $z_coord # coordinates of Atom 1 for normal mode 2
...
...
$x_coord $y_coord $z_coord # coordinates of the last atom for the last normal mode
```

## 3) Berry-phase dipole moments

### 3.1 MD-based spectra

The Berry phase dipole moment vectors for the MD-based spectra should be given in Debye in the `.xyz` format, combined for all MD snapshots in a single file such as:

```bash
1 # because there is only 1 dipole vector for the whole system
                       #leave 1 blank line
'dummy_value' $dipole_x $dipole_y $dipole_z #first MD frame
1 
                       #leave 1 blank line
'dummy_value' $dipole_x $dipole_y $dipole_z #second MD frame
...
1 
                       #leave 1 blank line
'dummy_value' $dipole_x $dipole_y $dipole_z #last MD frame
```

### 3.2 Static spectra based on normal modes

The Berry phase dipole moment vectors for a static calculation should be given in Debye for each $6N$ displaced structure (for more information, see Section [Frequencies](frequency.md) combined in a single file, such as (ignoring the comment lines):

```bash
#dipoles coming from Atom 1 displaced in +x direction
1 # because there is only 1 dipole vector for the whole system
                       #leave 1 blank line
'dummy_value' $dipole_x $dipole_y $dipole_z 
#dipoles coming from Atom 1 displaced in +y direction
1 
                       
'dummy_value' $dipole_x $dipole_y $dipole_z 
#dipoles coming from Atom 1 displaced in +z direction
1 
                       
'dummy_value' $dipole_x $dipole_y $dipole_z 
#dipoles coming from Atom 2 displaced in +x direction
1 
                       
'dummy_value' $dipole_x $dipole_y $dipole_z 
...
...
#dipoles coming from Atom 2 displaced in +x direction
1 
                       
'dummy_value' $dipole_x $dipole_y $dipole_z 
...
#dipoles coming from Atom 1 displaced in -x direction
1 
                       
'dummy_value' $dipole_x $dipole_y $dipole_z 
...
#dipoles coming from the last atom displaced in -z direction
1 
                       
'dummy_value' $dipole_x $dipole_y $dipole_z 
```

## 4) Wannier centers

If any MD-based spectrum is being calculated from Wannier centers, then the `wannier.xyz` for whole system (atoms + wannier centers labeled as "X") should be given for each MD frame in the `.xyz` file format. An example could be seen [here](_static/COF-1_solv_wannier_free.xyz).

## 5) DFPT polarizabilities

### 5.1 MD-based spectra

The DFPT polarizabilities for the MD-based spectra should be given in Å $^3$ in the `.xyz` format. There should be three separate files in total, such that the first file contains $\alpha_{xx}$, $\alpha_{xy}$, $\alpha_{xz}$; the second file contains $\alpha_{yx}$, $\alpha_{yy}$, $\alpha_{yz}$; and the third file contains $\alpha_{zx}$, $\alpha_{zy}$, $\alpha_{zz}$. Each file should have the format:


```bash
1 # because there is only 1 polarizability tensor for the whole system
                       #leave 1 blank line
'dummy_value' $alpha_ix $alpha_iy $alpha_iz #first MD frame
1 
                       #leave 1 blank line
'dummy_value' $alpha_ix $alpha_iy $alpha_iz #second MD frame
...
1 
                       #leave 1 blank line
'dummy_value' $alpha_ix $alpha_iy $alpha_iz #last MD frame
```

where $i$ is $x$, $y$, or $z$.

### 5.2 Static spectra based on normal modes

The DFPT polarizabilities for a static calculation should be given in Å $^3$ for each $6N$ displaced structure (for more information, see Section [Frequencies](frequency.md) combined in a single file, such as (ignoring the comment lines):

```bash
#leave 6 blank lines
#Atom 1 displaced in +x direction
'dummy_value' $alpha_xx $alpha_yy $alpha_zz 
'dummy_value' $alpha_xy $alpha_xz $alpha_yz 
'dummy_value' $alpha_yx $alpha_zx $alpha_zy 
#Atom 1 displaced in +y direction
#leave 6 blank lines
'dummy_value' $alpha_xx $alpha_yy $alpha_zz 
'dummy_value' $alpha_xy $alpha_xz $alpha_yz 
'dummy_value' $alpha_yx $alpha_zx $alpha_zy 
#Atom 1 displaced in +z direction
#leave 6 blank lines
'dummy_value' $alpha_xx $alpha_yy $alpha_zz 
'dummy_value' $alpha_xy $alpha_xz $alpha_yz 
'dummy_value' $alpha_yx $alpha_zx $alpha_zy 
...
#Atom 2 displaced in +x direction
#leave 6 blank lines
'dummy_value' $alpha_xx $alpha_yy $alpha_zz 
'dummy_value' $alpha_xy $alpha_xz $alpha_yz 
'dummy_value' $alpha_yx $alpha_zx $alpha_zy 
...
#Atom 1 displaced in -x direction
#leave 6 blank lines
'dummy_value' $alpha_xx $alpha_yy $alpha_zz 
'dummy_value' $alpha_xy $alpha_xz $alpha_yz 
'dummy_value' $alpha_yx $alpha_zx $alpha_zy 
...
#Last atom displaced in -z direction
#leave 6 blank lines
'dummy_value' $alpha_xx $alpha_yy $alpha_zz 
'dummy_value' $alpha_xy $alpha_xz $alpha_yz 
'dummy_value' $alpha_yx $alpha_zx $alpha_zy 
```

The reason for leaving six blank spaces is because the user can directly append the CP2K DFPT polarizabilities, an example is shown below:

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

Vibrant only reads the polarizabilities in Å $^3$, ignoring the a.u. block.

## 6) Time-dependent Berry phase dipole moments coming from RT-TDDFT

### 6.1 MD-based spectra

The time-dependent Berry phase dipole moment vectors for the MD-based resonance Raman (RR) or absorption spectra should be given in Debye in the `.xyz` format, combined for all MD snapshots in a single file such as:

```bash
'dummy_value' #first MD frame
                       #leave 1 blank line
'dummy_value' $dipole_x $dipole_y $dipole_z #unperturbed dipoles
'dummy_value' $dipole_x $dipole_y $dipole_z #RT-TDDFT step 1
'dummy_value' $dipole_x $dipole_y $dipole_z #RT-TDDFT step 2
...
'dummy_value' $dipole_x $dipole_y $dipole_z #last RT-TDDFT step
'dummy_value' #second MD frame
                       #leave 1 blank line
'dummy_value' $dipole_x $dipole_y $dipole_z #unperturbed dipoles
'dummy_value' $dipole_x $dipole_y $dipole_z #RT-TDDFT step 1
...
...
'dummy_value' #last MD frame
                       #leave 1 blank line
'dummy_value' $dipole_x $dipole_y $dipole_z #unperturbed dipoles
'dummy_value' $dipole_x $dipole_y $dipole_z #RT-TDDFT step 1
...
'dummy_value' $dipole_x $dipole_y $dipole_z #last RT-TDDFT step
```

Here the unperturbed dipole moments refer to $\mu^{0}_{i}$ described in Section [Resonance Raman (RR) and absorption spectra](Resonance_Raman_spec.md). Three files should be given in total, each containing the unperturbed dipole moments + time-dependent dipole moments perturbed under electric field applied in $x$, $y$ or $z$ direction.

### 6.2 Static spectra based on normal modes

The time-dependent Berry phase dipole moment vectors for a static RR or absorption spectrum calculation should be given in Debye for each $6N$ displaced structure (for more information, see Section [Frequencies](frequency.md) combined in a single file, such as (ignoring the comment lines):

```bash
#dipoles coming from Atom 1 displaced in +x direction
'dummy_value'
                       #leave 1 blank line
'dummy_value' $dipole_x $dipole_y $dipole_z #unperturbed dipoles
'dummy_value' $dipole_x $dipole_y $dipole_z #RT-TDDFT step 1
'dummy_value' $dipole_x $dipole_y $dipole_z #RT-TDDFT step 2
...
'dummy_value' $dipole_x $dipole_y $dipole_z #last RT-TDDFT step
#dipoles coming from Atom 1 displaced in +y direction
'dummy_value' 
                       #leave 1 blank line
'dummy_value' $dipole_x $dipole_y $dipole_z #unperturbed dipoles
'dummy_value' $dipole_x $dipole_y $dipole_z #RT-TDDFT step 1
...
#dipoles coming from Atom 1 displaced in +z direction
'dummy_value' 
                       #leave 1 blank line
'dummy_value' $dipole_x $dipole_y $dipole_z #unperturbed dipoles
'dummy_value' $dipole_x $dipole_y $dipole_z #RT-TDDFT step 1
...
#dipoles coming from Atom 2 displaced in +x direction
'dummy_value'
                       #leave 1 blank line
'dummy_value' $dipole_x $dipole_y $dipole_z #unperturbed dipoles
'dummy_value' $dipole_x $dipole_y $dipole_z #RT-TDDFT step 1
...
...
#dipoles coming from Atom 1 displaced in -x direction
'dummy_value'
                       #leave 1 blank line
'dummy_value' $dipole_x $dipole_y $dipole_z #unperturbed dipoles
'dummy_value' $dipole_x $dipole_y $dipole_z #RT-TDDFT step 1
...
...
#dipoles coming from the last atom displaced in -z direction
'dummy_value'
                       #leave 1 blank line
'dummy_value' $dipole_x $dipole_y $dipole_z #unperturbed dipoles
'dummy_value' $dipole_x $dipole_y $dipole_z #RT-TDDFT step 1
...
'dummy_value' $dipole_x $dipole_y $dipole_z #last RT-TDDFT step
```

Similarly to the MD-based spectra, here the unperturbed dipole moments refer to $\mu^{0}_{i}$ described in Section [Resonance Raman (RR) and absorption spectra](Resonance_Raman_spec.md). Three files should be given in total, each containing the unperturbed dipole moments + time-dependent dipole moments perturbed under electric field applied in $x$, $y$ or $z$ direction.

