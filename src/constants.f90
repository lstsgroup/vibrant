!
!   Copyright 2025 Ekin E. Winogradow, Johannes Scheffler, Moritz Leucke, Dorothea Golze
!
!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!
!       http://www.apache.org/licenses/LICENSE-2.0
!
!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.
!

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

   PUBLIC ::   speed_light, const_planck, const_permit, pi, const_charge, const_boltz, joule_unit, debye, &
             debye2cm, ev_unit, action_unit, bohr2ang, hartreebohr2evang, am_u, at_u, ang, fs2s, reccm2ev, &
             hessian_factor, avo_num, au2vm, cm2m, a3_to_debye_per_e, speed_light_au, reccm2au

   ! Constants

   ! Speed of light in vacuum [cm/s]
   REAL(kind=dp), PARAMETER                            :: speed_light = 2.9979246e+10_dp
   
   ! Speed of light in vacuum [a.u.]
   REAL(kind=dp), PARAMETER                            :: speed_light_au = 137_dp

   ! Planck constant [m^2*kg/s] or [J.s]
   REAL(kind=dp), PARAMETER                            :: const_planck = 6.62607015e-34_dp

   !Permittivity of vacuum [F/m]?
   REAL(kind=dp), PARAMETER                            :: const_permit = 8.8541878128e-12_dp

   ! Pi
   REAL(kind=dp), PARAMETER                            :: pi = 3.14159_dp

   ! Elementary charge [C] and [eV] -> [J]
   REAL(kind=dp), PARAMETER                            :: const_charge = 1.602176565e-19_dp

   ! Boltzmann constant [m^2*kg*s^-2*K-1] or [J/K]
   REAL(kind=dp), PARAMETER                            :: const_boltz = 1.380649e-23_dp
   
   ! Avogadro's number [mol^-1]
   REAL(kind=dp), PARAMETER                            :: avo_num = 6.02214e+23_dp

   ! Conversion factors

   ! [Hartree] -> [J]
   REAL(kind=dp), PARAMETER                            :: joule_unit = 4.359744722e-18

   ! [a.u. of action] -> [J*s]
   REAL(kind=dp), PARAMETER                            :: action_unit = 1.054571817e-34_dp

   ! [Debye] -> [a.u.]
   REAL(kind=dp), PARAMETER                            :: debye = 0.393456_dp
   
   ! [Debye] -> [C*m]
   REAL(kind=dp), PARAMETER                            :: debye2cm = 3.33564e-30_dp
   
   ! [Debye] -> [C*m]
   REAL(kind=dp), PARAMETER                            :: reccm2au = 4.556335e-6_dp 
   
   ! [cm] -> [m]
   REAL(kind=dp), PARAMETER                            :: cm2m = 0.01_dp
   
   ! [a.u.] -> [V/m]
   REAL(kind=dp), PARAMETER                            :: au2vm = 5.14220675112e+11_dp
   
   ! [A^3] -> [Debye/E (a.u.)]
   REAL(kind=dp), PARAMETER                            :: a3_to_debye_per_e = 1.713005_dp

   ! [a.u.] -> [eV]
   REAL(kind=dp), PARAMETER                            :: ev_unit = 27.211386_dp

   ! [Bohr] -> [Angstrom]
   REAL(kind=dp), PARAMETER                            :: bohr2ang = 0.5291772109_dp

   ! [Hartree/Bohr] -> [eV/Angstrom]
   REAL(kind=dp), PARAMETER                            :: hartreebohr2evang = 51.42208619083232_dp

   ! [a.m.u.] -> [kg]
   REAL(kind=dp), PARAMETER                            :: am_u = 1.6605390666e-27_dp
   
   ! [a.t.u.] -> [s]
   REAL(kind=dp), PARAMETER                            :: at_u = 2.4188843265864e-17_dp

   ! [Angstrom] -> [m]
   REAL(kind=dp), PARAMETER                            :: ang = 1.0e-10_dp

   ! [fs] -> [s]
   REAL(kind=dp), PARAMETER                            :: fs2s = 1.0e-15_dp

   ! [1/cm] -> [eV]
   REAL(kind=dp), PARAMETER                            :: reccm2ev = 0.000124_dp

   ! [eV/(a.m.u*Ang^2)] -> [J/(kg*m^2)]
   REAL(kind=dp), PARAMETER                            :: hessian_factor = REAL(const_charge/(am_u*ang*ang), kind=dp)

    !! MAGIC NUMBERS
    ! spec_ir
    REAL(kind=dp), PARAMETER                            :: ir_factor = 42.256_dp

    ! raman
    REAL(kind=dp), PARAMETER                            :: r_factor = -1.438777_dp

    ! refquencey
    REAL(kind=dp), PARAMETER                            :: frq_factor = 1.883652d-4

    ! refquencey
    REAL(kind=dp), PARAMETER                            :: power_factor = 7.211349d-9

    ! refquencey
    REAL(kind=dp), PARAMETER                            :: rtp_factor = 1.23984198e-4

    ! refquencey
    REAL(kind=dp), PARAMETER                            :: int_rman_factor = 1d-29*0.421_dp

    ! refquencey
    REAL(kind=dp), PARAMETER                            :: int_ir_factor = 3047.2310_dp

END MODULE constants
