! **************************************************************************************************
!> \brief Definition of physical constants:
!>
!>      speed_light : Speed of light in vacuum [cm/s]
!>      const_planck : Planck constant [m^2*kg/s] or [J.s]
!>      const_permit : F*m^−1 ?permittivity of vacuum [F/m]?
!>      pi : pi
!>      const_charge :  Charge ?Elementary charge [C]?
!>      const_boltz : Boltzmann constant [m^2*kg*s^-2*K-1] or [J/K]
!>      damping_constant : Damping Constant [ev]
!>      joule_unit : J
!>      debye : debye
!>      ev_unit : ev
!>      action_unit : J.s
!>      bohr2ang : Conversion factors [Bohr] -> [Angstrom]
!>      hartreebohr2evang : Conversion factors [Hartree/Bohr] -> [eV/Angstrom]
!>      at_u : Conversion factors atomic mass unit-kilogram relationship
!>      ang : Conversion factors [Angstrom] ->  (wave numbers)
!>      fs2s : Conversion factors [fs] -> [s]
!>      reccm2ev : Conversion factors [1/cm] -> [eV]
!>      t_cor : correlation depth?
!>      temp : Temperature in K
!>      hessian_factor : ???
!> \note
!> \par History
!>      - Style adapted for from in CP2K common/physcon.F
!> \author Ekin Bas and Johannes Scheffler
! **************************************************************************************************

MODULE constants

    USE kinds, ONLY: dp

    IMPLICIT NONE
    
    PRIVATE 
    
    PUBLIC ::   speed_light, const_planck, const_permit, pi, const_charge, const_boltz, damping_constant, joule_unit, debye, &
    ev_unit, action_unit, bohr2ang, hartreebohr2evang, at_u, ang, fs2s, reccm2ev, t_cor, temp, hessian_factor
    
    ! Constants
    
    ! Speed of light in vacuum [cm/s]
    REAL(kind=dp), PARAMETER                            :: speed_light = 2.9979246e+10_dp  
    
    ! Planck constant [m^2*kg/s] or [J.s]
    REAL(kind=dp), PARAMETER                            :: const_planck = 6.62607015e-34_dp 

    !Permittivity of vacuum [F/m]?
    REAL(kind=dp), PARAMETER                            :: const_permit = 8.8541878128e-12

    ! Pi
    REAL(kind=dp), PARAMETER                            :: pi = 3.14159_dp

    ! Elementary charge [C] and [eV] -> [J]
    REAL(kind=dp), PARAMETER                            :: const_charge = 1.602176565e-19_dp

    ! Boltzmann constant [m^2*kg*s^-2*K-1] or [J/K]
    REAL(kind=dp), PARAMETER                            :: const_boltz = 1.380649e-23_dp 


    ! Conversion factors
    
    ! [Hartree] -> [J]
    REAL(kind=dp), PARAMETER                            :: joule_unit = 4.359744722e-18

    ! [a.u. of action] -> [J*s]
    REAL(kind=dp), PARAMETER                            :: action_unit = 1.054571817e-34_dp 
    
    ! [Debye] -> [a.u.]
    REAL(kind=dp), PARAMETER                            :: debye = 0.393456_dp
    
    ! [a.u.] -> [eV]
    REAL(kind=dp), PARAMETER                            :: ev_unit = 27.211386_dp

    ! [Bohr] -> [Angstrom]
    REAL(kind=dp), PARAMETER                            :: bohr2ang = 0.5291772109_dp 

    ! [Hartree/Bohr] -> [eV/Angstrom]
    REAL(kind=dp), PARAMETER                            :: hartreebohr2evang = 51.42208619083232_dp 
    
    ! [a.m.u.] -> [kg]
    REAL(kind=dp), PARAMETER                            :: at_u = 1.6605390666e-27_dp 

    ! [Angstrom] -> [m]
    REAL(kind=dp), PARAMETER                            :: ang = 1.0e-10_dp

    ! [fs] -> [s]
    REAL(kind=dp), PARAMETER                            :: fs2s = 1.0e-15_dp

    ! [1/cm] -> [eV]
    REAL(kind=dp), PARAMETER                            :: reccm2ev = 0.000124_dp 
    
    ! [eV/(a.m.u*Ang^2)] -> [J/(kg*m^2)]
    REAL(kind=dp), PARAMETER                            :: hessian_factor = REAL(const_charge/(at_u*ang*ang), kind=dp)


    ! Input parameters 
    
    ! Damping Constant [eV]
    REAL(Kind=dp), PARAMETER                            :: damping_constant = 0.10_dp 

    ! Correlation depth for MD
    INTEGER, PARAMETER                                  :: t_cor = 1024

    ! Temperature [K]
    REAL(kind=dp), PARAMETER                            :: temp = 300.0_dp       
    

END MODULE constants
