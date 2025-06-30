program vibrant_input_trail

  USE kinds,      ONLY: dp, default_string_length
  USE types,      ONLY: section_type
  USE read_input, ONLY: parse_input,&
                        parse_command_line

  IMPLICIT NONE

  TYPE(section_type)                            :: input_sections
  CHARACTER(LEN=default_string_length)          :: input_file_name
  
  CALL parse_command_line(input_sections, input_file_name)
  CALL parse_input(input_sections, input_file_name)

  !*** TEST
  write(*,*) "abc", input_sections%system%cell%abc(1:3)
  write(*,*) "alpha_beta_gamma", input_sections%system%cell%alpha_beta_gamma(1:3)
  write(*,*) "coordinate_file: ", input_sections%system%coordinates%xyz_filename
  write(*,*) "md files: ", input_sections%md%trajectory_file
  write(*,*) "velocity files:  ", input_sections%md%velocity_file
  write(*,*) "snapshot_time_step  ", input_sections%md%snapshot_time_step
  write(*,*) "correlation_depth  ", input_sections%md%correlation_depth
  
end program
