
Keyword Glossary
================

link to a keyword: `See temperature keyword <#keyword-temperature>`__


Block: global
-------------

.. code-block:: text

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

.. code-block:: text

    &system
        filename 
        type_traj 
        mass_weighting 
    &end system
    
.. keyword:: filename
   :section: system
   :type: string
    
    This is the path of the file containing the system coordinates.

.. keyword:: type_traj
   :section: system
   :type: string
    
    Define type of provided trajectory. Possible values are: 

    * **pos** for a trajectory of positions
    * **vel** or a trajectory of velocities

.. keyword:: mass_weighting
   :section: system
   :type: string
    
   Define if mass weighting should be applied for spectra computation. Possible values are: 

    * **y** yes, apply mass weighting 
    * **n** no, dont apply mass weighting 

Subblock: cell
~~~~~~~~~~~~~~~~
.. code-block:: text

    &cell
        box_x 
        box_y 
        box_z 
        angle_alpha 
        angle_beta 
        angle_gamma
    &end cell

.. keyword:: cell_type
    :section: system/cell
    :type: string

    Define type of the simulation cell. Possible values are: 
    * **orthorhombic** yes, apply mass weighting 
    * **hexagonal** no, dont apply mass weighting 
    * **triclinic** no, dont apply mass weighting 
    CHECK AGAIN

.. keyword:: box_x
    :section: system/cell
    :type: float

    Length of the box lenght :math:`\mathbf{a}`. Uses the same unit as coordinates.

.. keyword:: box_y
    :section: system/cell
    :type: float

    Length of the box lenght :math:`\mathbf{b}`. Uses the same unit as coordinates.

.. keyword:: box_z
    :section: system/cell
    :type: float

    Length of the box lenght :math:`\mathbf{c}`. Uses the same unit as coordinates.

.. keyword:: angle_alpha
    :section: system/cell
    :type: float
    :unit: deg

    Cell angle  :math:`\alpha`.

.. keyword:: angle_beta
    :section: system/cell
    :type: float
    :unit: deg

    Cell angle :math:`\beta`.

.. keyword:: angle_gamma
    :section: system/cell
    :type: float
    :unit: deg

    Cell angle :math:`\gamma`.

Subblock: fragment
~~~~~~~~~~~~~~~~
.. code-block:: text

    &fragment
        atom_list 
    &end fragment

.. keyword:: atom_list
    :section: fragment
    :type: list[int]

    List of atomic indices that belong to a fragment (one entry per fragment).


Block: md
-------------

.. code-block:: text

    &md
        time_step 
        correlation_depth 
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

Block: static
-------------

.. code-block:: text
    
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
   :unit: Angstrom

   This is the finite difference displacement.

.. keyword:: diag_hessian
   :section: static
   :type: float
   :unit: Angstrom

   Define if the Hessian should be diagonalized. Possible values are:

   * **y** — yes, diagonalize the Hessian --> provide `force_file <#keyword-force_file>`__
   * **n** — no, load existing data --> provide `normal_freq_file <#keyword-normal_freq_file>`__ and `normal_displ_file <#keyword-normal_displ_file>`__

   .. keyword:: normal_freq_file
    :section: static
    :type: string

    This is the path of the file containing the systems normal mode frequencies.

.. keyword:: normal_displ_file
    :section: static
    :type: string

    This is the path of the file containing the systems normal mode displacements.

.. keyword:: write_mol_file
    :section: static
    :type: string

    Flag for writing a Molden file for analyzing static spectra.

Block: hessian
~~~~~~~~~~~~~~~~
.. code-block:: text
    
    &hessian
        force_file
    &end hessian
.. keyword:: force_file
   :section: static
   :type: string

   This is the path of the file containing the system forces to build hessian.

Block: dipoles
-------------

.. code-block:: text
    
    &dipoles
        type_dipole
        field_strength
        dip_file 
        dip_x_file
        dip_y_fil
        dip_z_file
        static_pol_file
    &end dipoles

.. keyword:: type_dipole
    :section: dipoles
    :type: string
    :default: berry
    
    Define type of provided dipole moments. Possible values are: 

    * **berry** for berry phase dipolemoments
    * **wannier** for wannier dipolemoments
    * **dfpt** for density functional perturbation theory
    check again

.. keyword:: field_strength
    :section: dipoles
    :type: float
    :unit: a.u.
    
    Strength of the applied electric field in atomic units.
    
.. keyword:: dip_file
    :section: dipoles
    :type: string
    
    This is the path of the file containing the systems dipole moments (if only a single file is used).
    
.. keyword:: dip_x_file
    :section: dipoles
    :type: string
    
    This is the path of the file containing the systems dipole moments obtained under an electric field in the x-direction.
    
.. keyword:: dip_y_file
    :section: dipoles
    :type: string
    
    This is the path of the file containing the systems dipole moments obtained under an electric field in the y-direction.
    
.. keyword:: dip_z_file
    :section: dipoles
    :type: string
    
    This is the path of the file containing the systems dipole moments obtained under an electric field in the z-direction.
    
.. keyword:: static_pol_file
    :section: dipoles
    :type: string
    
    This is the path of the file containing the systems static polarizabilities.

Block: rtp
-------------

.. code-block:: text

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

    Flag to enable Padé interpolation (“y” or “n”).

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

.. code-block:: text

    &raman
        laser_in 
    &end raman

.. keyword:: laser_in
    :section: raman
    :type: float or list[float]
    :unit: eV or cm-1 ?? not clear from readinput. check
    :default: 0.5

    Incoming laser energy in eV. Multiple values (max 10) may be specified.
