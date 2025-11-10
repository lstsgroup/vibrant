
Keyword Glossary
================

Block: global
-------------

.. code-block:: bash

    &global
        spectra 
        temperature
        fwhm
        spectra_verbosity
    &end global

.. keyword:: spectra
    :section: global
    :type: string
    :default: None

    Type of spectra to calculate. Possible values are:

    * **NMA**, **IR**, **R** — Normal Mode Analysis, Static IR, Static Raman
    * **P**, **MD-IR**, **MD-R** — Power Spectrum, Dynamic IR, Dynamic Raman 
    * **ABS**, **RR**, **MD-RR** — Absorption Spectrum, Static Resonance Raman, Dynamic Resonance Raman


.. keyword:: temperature
    :section: global
    :type: float
    :default: 300.0
    :unit: K

    This is the Temperature.

.. keyword:: fwhm
    :section: global
    :type: float
    :default: 10
    :unit: :math:`\mathrm{cm}^{-1}`

    Full width at half-maximum (fwhm) value for Gaussian broadening of static spectra. 

.. keyword:: spectra_verbosity
    :section: global
    :type: string
    :default: normal

    This is a verbosity setting for writing raman spectra files.


Block: system
-------------

.. code-block:: bash

    &system
        filename 
        type_traj 
        mass_weighting 
    &end system
    
.. keyword:: filename
    :section: system
    :unit:  Ang
    :type: string

    This is the path of the file containing the system coordinates. When velocities are provided, the unit Ang/fs is assumed.

.. keyword:: type_traj
    :section: system
    :type: string

    Define type of provided trajectory. Possible values are: 

    * **pos** — for a trajectory of positions
    * **vel** — for a trajectory of velocities

.. keyword:: mass_weighting
    :section: system
    :type: string

    Define if mass weighting should be applied for spectra computation. Possible values are: 

    * **y** — yes, apply mass weighting 
    * **n** — no, dont apply mass weighting 

Subblock: cell
~~~~~~~~~~~~~~~~
.. code-block:: bash

    &cell
        # Option 1
        box_xyz
        angle_alpha_beta_gamma 

        # Option 2
        lattice_x
        lattice_y
        lattice_z
    &end cell

.. keyword:: box_xyz
    :section: system/cell
    :unit: Ang
    :type: list[float]

    Length of the box lenght :math:`\mathbf{x}` :math:`\mathbf{y}` :math:`\mathbf{z}`. The box lengths can be specified individually using the keywords box_x, box_y, and box_z.

.. keyword:: angle_alpha_beta_gamma
    :section: system/cell
    :unit: deg
    :type: list[float]

    Cell angle  :math:`\alpha` :math:`\beta` :math:`\gamma`. The cell angle can be specified individually using the keywords angle_alpha, angle_beta, and angle_gamma.

.. keyword:: lattice_x
    :section: system/cell
    :unit: Ang
    :type: list[float]

    Lattice vector in x direction of system.

.. keyword:: lattice_y
    :section: system/cell
    :unit: Ang
    :type: list[float]


    Lattice vector in y direction of system.

.. keyword:: lattice_z
    :section: system/cell
    :unit: Ang
    :type: list[float]


    Lattice vector in z direction of system.

Subblock: fragment
~~~~~~~~~~~~~~~~
.. code-block:: bash

    &fragment
        atom_list 
    &end fragment

.. keyword:: atom_list
    :section: fragment
    :type: list[int]

    List of atomic indices that belong to a fragment (one entry per fragment).


Block: md
-------------

.. code-block:: bash

    &md
        time_step 
        correlation_depth
        write_acf_file 
    &end md

.. keyword:: time_step
    :section: md
    :type: float
    :unit: fs

    Time step of molecular dynamics trajectory.

.. keyword:: correlation_depth
    :section: md
    :type: int

    Correlation depth used in autocorrelation function.

.. keyword:: write_acf_file
    :section: md
    :type: string

    Flag for writing the time-autocorrelation data.

Block: static
-------------

.. code-block:: bash
    
    &static
        displacement
        diag_hessian 
        normal_freq_file
        normal_displ_file
        write_mol_file
    &end static


.. keyword:: displacement
   :section: static
   :type: float
   :unit: Ang

   This is the finite difference displacement.

.. keyword:: diag_hessian
    :section: static
    :type: float

    Define if the Hessian should be diagonalized. Possible values are:

    * **y** — yes, diagonalize the Hessian --> provide `force_file <#keyword-force_file>`__
    * **n** — no, load existing data --> provide `normal_freq_file <#keyword-normal_freq_file>`__ and `normal_displ_file <#keyword-normal_displ_file>`__

.. keyword:: normal_freq_file
    :section: static
    :unit: :math:`\mathrm{cm}^{-1}`
    :type: string

    This is the path of the file containing the systems normal mode frequencies.

.. keyword:: normal_displ_file
    :section: static
    :unit: -
    :type: string

    This is the path of the file containing the systems normal mode displacements.

.. keyword:: write_mol_file
    :section: static
    :type: string

    Flag for writing a Molden file for analyzing static spectra.

Block: hessian
~~~~~~~~~~~~~~~~
.. code-block:: bash
    
    &hessian
        force_file
    &end hessian

.. keyword:: force_file
    :section: static
    :unit: Hartree/Bohr (a.u.)
    :type: string

    This is the path of the file containing the system forces to build hessian.

Block: dipoles
-------------

.. code-block:: bash
    
    &dipoles
        type_dipole
        dip_file 
        dip_x_file
        dip_y_fil
        dip_z_file
    &end dipoles

.. keyword:: type_dipole
    :section: dipoles
    :type: string
    :default: berry

    Defines type of the provided dipole moments. Possible values are: 

    * **berry** — for Berry phase dipole moments
    * **wannier** — for Wannier centers
    
.. keyword:: dip_file
    :section: dipoles
    :unit:  Debye
    :type: string
    
    This is the path of the file containing the systems dipole moments (if only a single file is used).
    
.. keyword:: dip_x_file
    :section: dipoles
    :unit:  Debye
    :type: string
    
    This is the path of the file containing the systems dipole moments obtained under an electric field in the x-direction.
    
.. keyword:: dip_y_file
    :section: dipoles
    :unit:  Debye
    :type: string
    
    This is the path of the file containing the systems dipole moments obtained under an electric field in the y-direction.
    
.. keyword:: dip_z_file
    :section: dipoles
    :unit:  Debye
    :type: string
    
    This is the path of the file containing the systems dipole moments obtained under an electric field in the z-direction.
    


Block: polarizabilities
-------------

.. code-block:: bash
    
    &dipoles
        type_pol
        field_strength
        static_pol_file
    &end dipoles

.. keyword:: type_pol
    :section: polarizabilities
    :type: string
    :default: -

    Defines type of the provided polarizabilities. Possible values are: 

    * **induced** — induced dipole moments 
    * **analytical** — direct polarizabilities

.. keyword:: field_strength
    :section: dipoles
    :type: float
    :unit: Hartree/[:math:`e` Bohr] (a.u.)
    
    Strength of the applied electric field in atomic units.

.. keyword:: static_pol_file
    :section: dipoles
    :unit:  :math:`\mathrm{Ang}³`
    :type: string
    
    This is the path of the file containing the systems static polarizabilities.

Block: rtp
-------------

.. code-block:: bash

    &rtp
        rtp_time_step 
        rtp_framecount 
        check_pade
        pade_framecount 
        damping_constant 
    &end rtp

.. keyword:: rtp_time_step
    :section: rtp
    :type: float
    :unit: fs

    Time step used for the RTP (real-time propagation) analysis.

.. keyword:: rtp_framecount
    :section: rtp
    :type: int

    Number of frames used for RTP analysis.

.. keyword:: check_pade
    :section: rtp
    :type: string
    :default: 'n'

    Define if Padé interpolation should be used. Possible values are:

    * **y** — yes, apply Padé interpolation
    * **n** — no, dont apply Padé interpolation

.. keyword:: pade_framecount
    :section: rtp
    :type: int
    :default: 80000

    Number of frames used after Padé interpolation.

.. keyword:: damping_constant
    :section: rtp
    :type: float
    :unit: eV
    :default: 0.1

    Damping constant used in RTP spectrum calculations.

Block: raman
-------------

.. code-block:: bash

    &raman
        laser_in 
    &end raman

.. keyword:: laser_in
    :section: raman
    :type: float or list[float]
    :unit: eV
    :default: 0.5

    Incoming laser energy in eV. Multiple values (max 4) may be specified.
