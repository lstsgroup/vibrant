MODULE types
 
 USE kinds,    ONLY: dp, default_string_length
 
 IMPLICIT NONE

 PRIVATE

 PUBLIC :: section_type

 TYPE :: cell_section
    REAL(KIND=dp) :: abc(3)
    REAL(KIND=dp) :: alpha_beta_gamma(3)
    LOGICAL :: present = .false.
 END TYPE
 
 TYPE :: coordinates_section
    CHARACTER(LEN=default_string_length) :: xyz_filename
    LOGICAL :: present = .false.
 END TYPE 
!**********************************************************
 TYPE :: fragments_section
    INTEGER :: num_lists
    INTEGER, ALLOCATABLE:: atom_list(:,:)
    CHARACTER(LEN=default_string_length) :: pdb_filename
    LOGICAL :: present = .false.
 END TYPE
!**********************************************************

 TYPE :: system_section
    TYPE(cell_section) :: cell
    TYPE(coordinates_section) :: coordinates
    TYPE(fragments_section) :: fragments
 END TYPE
!**********************************************************

 TYPE :: md_section
    CHARACTER(LEN=default_string_length) :: trajectory_file
    CHARACTER(LEN=default_string_length) :: velocity_file
    REAL(KIND=dp) :: snapshot_time_step
    REAL(KIND=dp) :: correlation_depth
 END TYPE

!**********************************************************

 TYPE section_type
  TYPE(system_section) :: system
  TYPE(md_section)     :: md
 END TYPE section_type

CONTAINS

 SUBROUTINE  deallocate_sections()

  !... potentially some stuff
 END SUBROUTINE deallocate_sections

END MODULE types
