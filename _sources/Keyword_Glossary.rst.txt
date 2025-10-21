
Keyword Glossary
================

link to a keyword: `See temperature keyword <#keyword-temperature>`__

All Keywords and sections are listed here:

.. code-block:: text

    &global
        spectra MD-RR
        temperature 300
    &end global
    &dipoles
        ype_dipole berry
        ield_strength 0.005
        ip_x_file o-NP_RTP_dipoles_X_256.xyz
        ip_y_file o-NP_RTP_dipoles_Y_256.xyz
        ip_z_file o-NP_RTP_dipoles_Z_256.xyz
    &end dipoles
    &md
        time_step 2.5
        correlation_depth 64
    &end md
    &rtp
        rtp_time_step 0.0625
        rtp_framecount 256
        check_pade y
        pade_framecount 400
        damping_constant 0.1
    &end rtp
    &raman
        laser_in 1.17
    &end raman




Block: global
-------------

.. keyword:: spectra
   :section: global
   :type: string
   :default: None

   This is the spectra type. Possible values are: MD-RR, MD-RT, RR, RT, FT, ABS

.. keyword:: temperature
   :section: global
   :type: float
   :default: 300.0
   :unit: K

   This is the Temperature. with formula :math:`\omega = 2\pi f`.

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

   This is a verbosity setting for creating raman spectra files.


Block: system
-------------
.. keyword:: filename
   :section: system
   :type: string
    
    This is the name of coordinates file.

.. keyword:: type_traj
   :section: system
   :type: string
    
   Define type of provided trajectory. Values: 'pos' for positions, 'vel' for velocities.

.. keyword:: mass_weighting
   :section: system
   :type: string
    
   Define if mass weighting should be applied for spectra computation. Values: 'y' or 'n'

Subblock: cell
-------------
.. keyword:: cell_type
    :section: cell
    :type: string

    Crystal type of the simulation cell (e.g., orthorhombic, hexagonal, triclinic).

.. keyword:: box_x
    :section: cell
    :type: float

    Length of the cell vector :math:\mathbf{a} (x component). Uses the same unit as coordinates.

.. keyword:: box_y
    :section: cell
    :type: float

    Length of the cell vector :math:\mathbf{b} (y component). Uses the same unit as coordinates.

.. keyword:: box_z
    :section: cell
    :type: float

    Length of the cell vector :math:\mathbf{c} (z component). Uses the same unit as coordinates.

.. keyword:: angle_alpha
    :section: cell
    :type: float
    :unit: deg

    Cell angle :math:\alpha=\angle(\mathbf{b},\mathbf{c}).

.. keyword:: angle_beta
    :section: cell
    :type: float
    :unit: deg

    Cell angle :math:\beta=\angle(\mathbf{a},\mathbf{c}).

.. keyword:: angle_gamma
    :section: cell
    :type: float
    :unit: deg

    Cell angle :math:\gamma=\angle(\mathbf{a},\mathbf{b}).

Block: fragment
-------------
.. keyword:: atom_list
    :section: fragment
    :type: list[int]

List of atomic indices that belong to a fragment (one entry per fragment).


Block: md
-------------
.. keyword:: time_step
    :section: md
    :type: float
    :unit: fs

Time step of molecular dynamics trajectory.

.. keyword:: correlation_depth
    :section: md
    :type: int

    Correlation depth used in autocorrelationfunction

Block: static
-------------
.. keyword:: force_file
   :section: static
   :type: string

   This is the name of force file.

.. keyword:: displacement
   :section: static
   :type: float
   :unit: Angstrom

   This is the finite difference displacement.

.. keyword:: diag_hessian
   :section: static
   :type: float
   :unit: Angstrom

   This is the finite difference displacement.

   .. keyword:: normal_freq_file
    :section: static
    :type: string

File containing normal mode frequencies for static or mode-based analyses.

.. keyword:: normal_displ_file
    :section: static
    :type: string

File containing normal mode displacements or eigenvectors.

.. keyword:: write_mol_file
    :section: static
    :type: string

Flag to write a molecular structure file (values: “y” or “n”).

Block: hessian
-------------

Block: dipoles
-------------

.. keyword:: type_dipole
    :section: dipoles
    :type: string
    :default: berry
    
    Specifies the method used to obtain dipole moments (e.g., “berry”, “wannier”, “dfpt”).
    This value can be dfpt 1 berry or wannier? check again
.. keyword:: field_strength
    :section: dipoles
    :type: float
    :unit: a.u.
    
    Strength of the applied electric field in atomic units.
    
.. keyword:: dip_file
    :section: dipoles
    :type: string
    
    File containing dipole moments (if only a single file is used).
    
.. keyword:: dip_x_file
    :section: dipoles
    :type: string
    
    File containing dipole moments obtained under an electric field in the x-direction.
    
.. keyword:: dip_y_file
    :section: dipoles
    :type: string
    
    File containing dipole moments obtained under an electric field in the y-direction.
    
.. keyword:: dip_z_file
    :section: dipoles
    :type: string
    
    File containing dipole moments obtained under an electric field in the z-direction.
    
.. keyword:: static_pol_file
    :section: dipoles
    :type: string
    
    File containing static polarizabilities.

Block: rtp
-------------
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

.. keyword:: laser_in
    :section: raman
    :type: float or list[float]
    :unit: eV or cm-1 ?? not clear from readinput. check
    :default: 0.5

Incoming laser energy in eV. Multiple values (max 10) may be specified.
