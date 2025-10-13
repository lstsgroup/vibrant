
Keyword Glossary
================

.. code::

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



Block: system
-------------

.. keyword:: temperature
   :section: system
   :type: float
   :default: 300.0
   :unit: Kelvin

   Temperature in Kelvin used for sampling.

.. keyword:: pressure
   :section: system
   :type: float
   :default: 23.2
   :unit: kPa

   This is the Pressure. with formula $pV=nRT$.