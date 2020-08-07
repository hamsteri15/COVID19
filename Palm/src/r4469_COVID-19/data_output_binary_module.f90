!> @file data_output_binary_module.f90
!--------------------------------------------------------------------------------------------------!
! This file is part of the PALM model system.
!
! PALM is free software: you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software
! Foundation, either version 3 of the License, or (at your option) any later
! version.
!
! PALM is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
! A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with
! PALM. If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 2019-2019 Leibniz Universitaet Hannover
!--------------------------------------------------------------------------------------------------!
!
! Current revisions:
! ------------------
! 
! 
! Former revisions:
! -----------------
! $Id: data_output_binary_module.f90 4408 2020-02-14 10:04:39Z gronemeier $
! Enable character-array output
!
! 4232 2019-09-20 09:34:22Z knoop
! Bugfix: INCLUDE "mpif.h" must be placed after IMPLICIT NONE statement
!
! 4147 2019-08-07 09:42:31Z gronemeier
! corrected indentation according to coding standard
!
! 4141 2019-08-05 12:24:51Z gronemeier
! Initial revision
!
!
! Authors:
! --------
!> @author: Tobias Gronemeier
!
! Description:
! ------------
!> Binary output module to write output data into binary files.
!>
!> @todo Get iostat value of write statements.
!--------------------------------------------------------------------------------------------------!
 MODULE data_output_binary_module

    USE kinds

#if defined( __parallel ) && !defined( __mpifh )
    USE MPI
#endif

    IMPLICIT NONE

#if defined( __parallel ) && defined( __mpifh )
    INCLUDE "mpif.h"
#endif

    INTEGER, PARAMETER ::  charlen = 100  !< maximum length of character variables

    CHARACTER(LEN=*), PARAMETER ::  config_file_name = 'BINARY_TO_NETCDF_CONFIG'  !< name of config file
    CHARACTER(LEN=*), PARAMETER ::  mode_binary = 'binary'                        !< string to select operation mode of module
    CHARACTER(LEN=*), PARAMETER ::  file_prefix = 'BIN_'                          !< file prefix for binary files

    CHARACTER(LEN=charlen)      ::  file_suffix = ''             !< file suffix added to each file name
    CHARACTER(LEN=800)          ::  internal_error_message = ''  !< string containing the last error message
    CHARACTER(LEN=800)          ::  temp_string                  !< dummy string

    INTEGER ::  binary_file_lowest_unit = 1000  !< lowest unit number of all binary files created by this module
    INTEGER ::  config_file_unit                !< unit number of config file
    INTEGER ::  debug_output_unit               !< Fortran Unit Number of the debug-output file
    INTEGER ::  global_id_in_file = -1          !< value of global ID within a file
    INTEGER ::  master_rank                     !< master rank for tasks to be executed by single PE only
    INTEGER ::  next_available_unit             !< next unit number available for new file
    INTEGER ::  output_group_comm               !< MPI communicator addressing all MPI ranks which participate in output

    INTEGER, DIMENSION(:), ALLOCATABLE ::  files_highest_variable_id  !< highest assigned ID of variable or dimension in a file

    LOGICAL ::  binary_open_file_first_call = .TRUE.  !< true if binary_open_file routine was not called yet
    LOGICAL ::  config_file_open = .FALSE.            !< true if config file is opened and not closed
    LOGICAL ::  print_debug_output = .FALSE.          !< if true, debug output is printed

    SAVE

    PRIVATE

    INTERFACE binary_init_module
       MODULE PROCEDURE binary_init_module
    END INTERFACE binary_init_module

    INTERFACE binary_open_file
       MODULE PROCEDURE binary_open_file
    END INTERFACE binary_open_file

    INTERFACE binary_init_dimension
       MODULE PROCEDURE binary_init_dimension
    END INTERFACE binary_init_dimension

    INTERFACE binary_init_variable
       MODULE PROCEDURE binary_init_variable
    END INTERFACE binary_init_variable

    INTERFACE binary_write_attribute
       MODULE PROCEDURE binary_write_attribute
    END INTERFACE binary_write_attribute

    INTERFACE binary_stop_file_header_definition
       MODULE PROCEDURE binary_stop_file_header_definition
    END INTERFACE binary_stop_file_header_definition

    INTERFACE binary_write_variable
       MODULE PROCEDURE binary_write_variable
    END INTERFACE binary_write_variable

    INTERFACE binary_finalize
       MODULE PROCEDURE binary_finalize
    END INTERFACE binary_finalize

    INTERFACE binary_get_error_message
       MODULE PROCEDURE binary_get_error_message
    END INTERFACE binary_get_error_message

    PUBLIC &
       binary_finalize, &
       binary_get_error_message, &
       binary_init_dimension, &
       binary_stop_file_header_definition, &
       binary_init_module, &
       binary_init_variable, &
       binary_open_file, &
       binary_write_attribute, &
       binary_write_variable


 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize data-output module.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE binary_init_module( file_suffix_of_output_group, mpi_comm_of_output_group, &
                                master_output_rank,                                    &
                                program_debug_output_unit, debug_output, dom_global_id )

    CHARACTER(LEN=*), INTENT(IN) ::  file_suffix_of_output_group  !> file-name suffix added to each file;
                                                                  !> must be unique for each output group

    INTEGER, INTENT(IN) ::  dom_global_id              !< global id within a file defined by DOM
    INTEGER, INTENT(IN) ::  master_output_rank         !< MPI rank executing tasks which must be executed by a single PE
    INTEGER, INTENT(IN) ::  mpi_comm_of_output_group   !< MPI communicator specifying the rank group participating in output
    INTEGER, INTENT(IN) ::  program_debug_output_unit  !< file unit number for debug output

    LOGICAL, INTENT(IN) ::  debug_output  !< if true, debug output is printed


    file_suffix = file_suffix_of_output_group
    output_group_comm = mpi_comm_of_output_group
    master_rank = master_output_rank

    debug_output_unit = program_debug_output_unit
    print_debug_output = debug_output

    global_id_in_file = dom_global_id

 END SUBROUTINE binary_init_module

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Open binary file.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE binary_open_file( mode, file_name, file_id, return_value )

    CHARACTER(LEN=charlen)             ::  bin_filename = ''  !< actual name of binary file
    CHARACTER(LEN=charlen), INTENT(IN) ::  file_name          !< name of file
    CHARACTER(LEN=7)                   ::  my_rank_char       !< string containing value of my_rank with leading zeros
    CHARACTER(LEN=*),       INTENT(IN) ::  mode               !< operation mode

    CHARACTER(LEN=*), PARAMETER ::  routine_name = 'binary_open_file'  !< name of this routine

    INTEGER, INTENT(OUT) ::  file_id       !< file ID
    INTEGER              ::  my_rank       !< MPI rank of local processor
    INTEGER              ::  nranks        !< number of MPI ranks participating in output
    INTEGER, INTENT(OUT) ::  return_value  !< return value

    INTEGER, DIMENSION(:), ALLOCATABLE ::  files_highest_variable_id_tmp  !< temporary list of given variable IDs in file

    LOGICAL ::  file_exists  !< true if file to be opened already exists


    return_value = 0

#if defined( __parallel )
    CALL MPI_COMM_SIZE( output_group_comm, nranks, return_value )
    IF ( return_value == 0 )  CALL MPI_COMM_RANK( output_group_comm, my_rank, return_value )
    IF ( return_value == 0 )  THEN
       WRITE( my_rank_char, '("_",I6.6)' )  my_rank
    ELSE
       CALL internal_message( 'error', routine_name // ': MPI error' )
    ENDIF
#else
    nranks = 1
    my_rank = master_rank
    WRITE( my_rank_char, '("_",I6.6)' )  my_rank
#endif
!
!-- Check mode (not required, added for compatibility reasons)
    IF ( TRIM( mode ) == mode_binary )  CONTINUE
!
!-- Open binary config file for combining script
    IF ( return_value == 0  .AND.  binary_open_file_first_call )  THEN

       binary_open_file_first_call = .FALSE.
       config_file_unit = binary_file_lowest_unit

       IF ( my_rank == master_rank )  THEN
!
!--       Remove any pre-existing file
          INQUIRE( FILE=TRIM( config_file_name ) // TRIM( file_suffix ), &
                   EXIST=file_exists )

          IF ( file_exists )  THEN
             CALL internal_message( 'debug', routine_name //     &
                                    ': Remove existing file ' // &
                                    TRIM( config_file_name ) // TRIM( file_suffix ) )
             !> @note Fortran2008 feature 'EXECUTE_COMMAND_LINE' not yet supported by
             !>       PGI 18.10 compiler. Hence, non-standard 'SYSTEM' call must be used
             ! CALL EXECUTE_COMMAND_LINE(                                                &
             !         COMMAND='rm ' // TRIM( config_file_name ) // TRIM( file_suffix ), &
             !         WAIT=.TRUE., EXITSTAT=return_value )
             CALL SYSTEM( 'rm ' // TRIM( config_file_name ) // TRIM( file_suffix ) )
          ENDIF

          OPEN( config_file_unit, FILE=TRIM( config_file_name ) // TRIM( file_suffix ), &
                FORM='UNFORMATTED', STATUS='NEW', IOSTAT=return_value )

          IF ( return_value == 0 )  THEN

             config_file_open = .TRUE.
!
!--          Write some general information to config file
             WRITE( config_file_unit )  nranks
             WRITE( config_file_unit )  master_rank
             WRITE( config_file_unit )  LEN( file_prefix )
             WRITE( config_file_unit )  file_prefix
             WRITE( config_file_unit )  charlen
             WRITE( config_file_unit )  global_id_in_file

          ELSE

             return_value = 1
             CALL internal_message( 'error', routine_name // ': could not create config' )

          ENDIF

       ENDIF

       next_available_unit = binary_file_lowest_unit + 1

    ENDIF
!
!-- Initialize output file: open, write header, initialize variable/dimension IDs
    IF ( return_value == 0 )  THEN

       bin_filename = file_prefix // TRIM( file_name ) // TRIM( file_suffix ) // my_rank_char
!
!--    Remove any pre-existing file
       INQUIRE( FILE=TRIM( bin_filename ), EXIST=file_exists )

       IF ( file_exists )  THEN
          CALL internal_message( 'debug', routine_name // &
                                 ': remove existing file ' // TRIM( bin_filename ) )
          !> @note Fortran2008 feature 'EXECUTE_COMMAND_LINE' not yet supported by
          !>       PGI 18.10 compiler. Hence, non-standard 'SYSTEM' call must be used
          ! CALL EXECUTE_COMMAND_LINE( COMMAND='rm ' // TRIM( bin_filename ), &
          !                            WAIT=.TRUE., EXITSTAT=return_value )
          CALL SYSTEM( 'rm ' // TRIM( bin_filename ) )
       ENDIF
!
!--    Open binary file
       CALL internal_message( 'debug', routine_name // ': open file ' // TRIM( bin_filename ) )
       OPEN ( next_available_unit, FILE=TRIM( bin_filename ), &
              FORM='UNFORMATTED', STATUS='NEW', IOSTAT=return_value )

       IF ( return_value == 0 )  THEN
!
!--       Add file_name to config file
          IF ( my_rank == master_rank )  THEN
             WRITE( config_file_unit )  file_name
          ENDIF
!
!--       Save file ID and increase next file unit number
          file_id = next_available_unit
          next_available_unit = next_available_unit + 1
!
!--       Write some meta data to file
          WRITE ( file_id )  charlen
          WRITE ( file_id )  file_id
          WRITE ( file_id )  file_name
!
!--       Extend file-variable/dimension-ID list by 1 and set it to 0 for new file.
          IF ( ALLOCATED( files_highest_variable_id ) )  THEN
             ALLOCATE( files_highest_variable_id_tmp(SIZE( files_highest_variable_id )) )
             files_highest_variable_id_tmp = files_highest_variable_id
             DEALLOCATE( files_highest_variable_id )
             ALLOCATE( files_highest_variable_id(binary_file_lowest_unit+1:file_id) )
             files_highest_variable_id(:file_id-1) = files_highest_variable_id_tmp
             DEALLOCATE( files_highest_variable_id_tmp )
          ELSE
             ALLOCATE( files_highest_variable_id(binary_file_lowest_unit+1:file_id) )
          ENDIF
          files_highest_variable_id(file_id) = 0

       ELSE
          return_value = 1
          CALL internal_message( 'error', routine_name // &
                                 ': could not open file "' // TRIM( file_name ) // '"')
       ENDIF

    ENDIF

 END SUBROUTINE binary_open_file

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write attribute to file.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE binary_write_attribute( file_id, variable_id, attribute_name, &
               value_char, value_int8, value_int16, value_int32,          &
               value_real32, value_real64, return_value )

    CHARACTER(LEN=*), PARAMETER ::  routine_name = 'binary_write_attribute'  !< name of this routine

    CHARACTER(LEN=charlen), INTENT(IN)           ::  attribute_name        !< name of attribute
    CHARACTER(LEN=charlen), INTENT(IN), OPTIONAL ::  value_char  !< value of attribute
    CHARACTER(LEN=charlen)                       ::  attribute_type        !< data type of attribute
    CHARACTER(LEN=charlen)                       ::  output_string         !< output string

    INTEGER(KIND=1), INTENT(IN), OPTIONAL ::  value_int8   !< value of attribute
    INTEGER(KIND=2), INTENT(IN), OPTIONAL ::  value_int16  !< value of attribute
    INTEGER(KIND=4), INTENT(IN), OPTIONAL ::  value_int32  !< value of attribute

    INTEGER, INTENT(IN)  ::  file_id       !< file ID
    INTEGER, INTENT(IN)  ::  variable_id   !< variable ID
    INTEGER, INTENT(OUT) ::  return_value  !< return value

    REAL(KIND=4), INTENT(IN), OPTIONAL ::  value_real32  !< value of attribute
    REAL(KIND=8), INTENT(IN), OPTIONAL ::  value_real64  !< value of attribute


    return_value = 0

    CALL internal_message( 'debug', TRIM( routine_name ) // &
                           ': write attribute ' // TRIM( attribute_name ) )
!
!-- Write attribute to file
    output_string = 'attribute'
    WRITE( file_id )  output_string

    WRITE( file_id )  variable_id
    WRITE( file_id )  attribute_name

    IF ( PRESENT( value_char ) )  THEN
       attribute_type = 'char'
       WRITE( file_id )  attribute_type
       WRITE( file_id )  value_char
    ELSEIF ( PRESENT( value_int8 ) )  THEN
       attribute_type = 'int8'
       WRITE( file_id )  attribute_type
       WRITE( file_id )  value_int8
    ELSEIF ( PRESENT( value_int16 ) )  THEN
       attribute_type = 'int16'
       WRITE( file_id )  attribute_type
       WRITE( file_id )  value_int16
    ELSEIF ( PRESENT( value_int32 ) )  THEN
       attribute_type = 'int32'
       WRITE( file_id )  attribute_type
       WRITE( file_id )  value_int32
    ELSEIF ( PRESENT( value_real32 ) )  THEN
       attribute_type = 'real32'
       WRITE( file_id )  attribute_type
       WRITE( file_id )  value_real32
    ELSEIF ( PRESENT( value_real64 ) )  THEN
       attribute_type = 'real64'
       WRITE( file_id )  attribute_type
       WRITE( file_id )  value_real64
    ELSE
       return_value = 1
       CALL internal_message( 'error', TRIM( routine_name ) // &
                              ': no value given for attribute "' // TRIM( attribute_name ) // '"' )
    ENDIF

 END SUBROUTINE binary_write_attribute

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize dimension. Write information in file header
!> and save dimension values to be later written to file.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE binary_init_dimension( mode, file_id, dimension_id, variable_id, &
               dimension_name, dimension_type, dimension_length, return_value )

    CHARACTER(LEN=charlen), INTENT(IN) ::  dimension_name  !< name of dimension
    CHARACTER(LEN=charlen), INTENT(IN) ::  dimension_type  !< data type of dimension
    CHARACTER(LEN=charlen)             ::  output_string   !< output string
    CHARACTER(LEN=*),       INTENT(IN) ::  mode            !< operation mode

    CHARACTER(LEN=*), PARAMETER ::  routine_name = 'binary_init_dimension'  !< name of this routine

    INTEGER, INTENT(OUT) ::  dimension_id      !< dimension ID
    INTEGER, INTENT(IN)  ::  dimension_length  !< length of dimension
    INTEGER, INTENT(IN)  ::  file_id           !< file ID
    INTEGER, INTENT(OUT) ::  return_value      !< return value
    INTEGER, INTENT(OUT) ::  variable_id       !< variable ID


    return_value = 0

    CALL internal_message( 'debug', routine_name // ': init dimension ' // TRIM( dimension_name ) )
!
!-- Check mode (not required, added for compatibility reasons only)
    IF ( TRIM( mode ) == mode_binary )  CONTINUE
!
!-- Assign dimension ID
    dimension_id = files_highest_variable_id( file_id ) + 1
    files_highest_variable_id( file_id ) = dimension_id
!
!-- Define dimension in file
    output_string = 'dimension'
    WRITE( file_id )  output_string
    WRITE( file_id )  dimension_name
    WRITE( file_id )  dimension_id
    WRITE( file_id )  dimension_type
    WRITE( file_id )  dimension_length
!
!-- Define variable associated with dimension
    CALL binary_init_variable( mode, file_id, variable_id, dimension_name, dimension_type, &
                               (/ dimension_id /), is_global=.TRUE., return_value=return_value )
    IF ( return_value /= 0 )  THEN
       CALL internal_message( 'error', routine_name // &
                              ': init dimension "' // TRIM( dimension_name ) // '"' )
    ENDIF

 END SUBROUTINE binary_init_dimension

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize variable. Write information of variable into file header.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE binary_init_variable( mode, file_id, variable_id, variable_name, variable_type, &
                                  dimension_ids, is_global, return_value )

    CHARACTER(LEN=charlen)             ::  output_string   !< output string
    CHARACTER(LEN=charlen), INTENT(IN) ::  variable_name   !< name of variable
    CHARACTER(LEN=charlen), INTENT(IN) ::  variable_type   !< data type of variable
    CHARACTER(LEN=*),       INTENT(IN) ::  mode            !< operation mode

    CHARACTER(LEN=*), PARAMETER ::  routine_name = 'binary_init_variable'  !< name of this routine

    INTEGER, INTENT(IN)  ::  file_id       !< file ID
    INTEGER, INTENT(OUT) ::  variable_id   !< variable ID
    INTEGER, INTENT(OUT) ::  return_value  !< return value

    INTEGER, DIMENSION(:), INTENT(IN) ::  dimension_ids  !< list of dimension IDs used by variable

    LOGICAL, INTENT(IN)  ::  is_global  !< true if variable is global (same on all PE)


    return_value = 0

    CALL internal_message( 'debug', routine_name // ': init variable ' // TRIM( variable_name ) )
!
!-- Check mode (not required, added for compatibility reasons only)
    IF ( TRIM( mode ) == mode_binary )  CONTINUE
!
!-- Check if variable is global (not required, added for compatibility reasons only)
    IF ( is_global )  CONTINUE
!
!-- Assign variable ID
    variable_id = files_highest_variable_id( file_id ) + 1
    files_highest_variable_id( file_id ) = variable_id
!
!-- Write variable information in file
    output_string = 'variable'
    WRITE( file_id )  output_string
    WRITE( file_id )  variable_name
    WRITE( file_id )  variable_id
    WRITE( file_id )  variable_type
    WRITE( file_id )  SIZE( dimension_ids )
    WRITE( file_id )  dimension_ids

 END SUBROUTINE binary_init_variable

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Leave file definition state.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE binary_stop_file_header_definition( file_id, return_value )

    CHARACTER(LEN=charlen) ::  output_string  !< output string

    CHARACTER(LEN=*), PARAMETER ::  routine_name = 'binary_stop_file_header_definition'  !< name of this routine

    INTEGER, INTENT(IN)  ::  file_id       !< file ID
    INTEGER, INTENT(OUT) ::  return_value  !< return value


    return_value = 0

    WRITE( temp_string, * ) file_id
    CALL internal_message( 'debug', routine_name // &
                           ': finalize file definition (file_id=' // TRIM( temp_string ) // ')' )

    output_string = '*** end file header ***'
    WRITE( file_id )  output_string

 END SUBROUTINE binary_stop_file_header_definition

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write variable to file.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE binary_write_variable(                                                     &
               file_id, variable_id, bounds_start, value_counts, bounds_origin,        &
               is_global,                                                              &
               values_char_0d,   values_char_1d,   values_char_2d,   values_char_3d,   &
               values_int8_0d,   values_int8_1d,   values_int8_2d,   values_int8_3d,   &
               values_int16_0d,  values_int16_1d,  values_int16_2d,  values_int16_3d,  &
               values_int32_0d,  values_int32_1d,  values_int32_2d,  values_int32_3d,  &
               values_intwp_0d,  values_intwp_1d,  values_intwp_2d,  values_intwp_3d,  &
               values_real32_0d, values_real32_1d, values_real32_2d, values_real32_3d, &
               values_real64_0d, values_real64_1d, values_real64_2d, values_real64_3d, &
               values_realwp_0d, values_realwp_1d, values_realwp_2d, values_realwp_3d, &
               return_value )

    CHARACTER(LEN=charlen) ::  output_string  !< output string

    CHARACTER(LEN=*), PARAMETER ::  routine_name = 'binary_write_variable'  !< name of this routine

    CHARACTER(LEN=1), POINTER,             INTENT(IN), OPTIONAL                   ::  values_char_0d  !< output variable
    CHARACTER(LEN=1), POINTER, CONTIGUOUS, INTENT(IN), OPTIONAL, DIMENSION(:)     ::  values_char_1d  !< output variable
    CHARACTER(LEN=1), POINTER, CONTIGUOUS, INTENT(IN), OPTIONAL, DIMENSION(:,:)   ::  values_char_2d  !< output variable
    CHARACTER(LEN=1), POINTER, CONTIGUOUS, INTENT(IN), OPTIONAL, DIMENSION(:,:,:) ::  values_char_3d  !< output variable

    INTEGER, INTENT(IN)  ::  file_id       !< file ID
    INTEGER, INTENT(OUT) ::  return_value  !< return value
    INTEGER, INTENT(IN)  ::  variable_id   !< variable ID

    INTEGER, DIMENSION(:), INTENT(IN) ::  bounds_origin  !< starting index of each dimension
    INTEGER, DIMENSION(:), INTENT(IN) ::  bounds_start   !< starting index of variable
    INTEGER, DIMENSION(:), INTENT(IN) ::  value_counts   !< count of values along each dimension to be written

    INTEGER(KIND=1), POINTER,             INTENT(IN), OPTIONAL                   ::  values_int8_0d   !< output variable
    INTEGER(KIND=2), POINTER,             INTENT(IN), OPTIONAL                   ::  values_int16_0d  !< output variable
    INTEGER(KIND=4), POINTER,             INTENT(IN), OPTIONAL                   ::  values_int32_0d  !< output variable
    INTEGER(iwp),    POINTER,             INTENT(IN), OPTIONAL                   ::  values_intwp_0d  !< output variable
    INTEGER(KIND=1), POINTER, CONTIGUOUS, INTENT(IN), OPTIONAL, DIMENSION(:)     ::  values_int8_1d   !< output variable
    INTEGER(KIND=2), POINTER, CONTIGUOUS, INTENT(IN), OPTIONAL, DIMENSION(:)     ::  values_int16_1d  !< output variable
    INTEGER(KIND=4), POINTER, CONTIGUOUS, INTENT(IN), OPTIONAL, DIMENSION(:)     ::  values_int32_1d  !< output variable
    INTEGER(iwp),    POINTER, CONTIGUOUS, INTENT(IN), OPTIONAL, DIMENSION(:)     ::  values_intwp_1d  !< output variable
    INTEGER(KIND=1), POINTER, CONTIGUOUS, INTENT(IN), OPTIONAL, DIMENSION(:,:)   ::  values_int8_2d   !< output variable
    INTEGER(KIND=2), POINTER, CONTIGUOUS, INTENT(IN), OPTIONAL, DIMENSION(:,:)   ::  values_int16_2d  !< output variable
    INTEGER(KIND=4), POINTER, CONTIGUOUS, INTENT(IN), OPTIONAL, DIMENSION(:,:)   ::  values_int32_2d  !< output variable
    INTEGER(iwp),    POINTER, CONTIGUOUS, INTENT(IN), OPTIONAL, DIMENSION(:,:)   ::  values_intwp_2d  !< output variable
    INTEGER(KIND=1), POINTER, CONTIGUOUS, INTENT(IN), OPTIONAL, DIMENSION(:,:,:) ::  values_int8_3d   !< output variable
    INTEGER(KIND=2), POINTER, CONTIGUOUS, INTENT(IN), OPTIONAL, DIMENSION(:,:,:) ::  values_int16_3d  !< output variable
    INTEGER(KIND=4), POINTER, CONTIGUOUS, INTENT(IN), OPTIONAL, DIMENSION(:,:,:) ::  values_int32_3d  !< output variable
    INTEGER(iwp),    POINTER, CONTIGUOUS, INTENT(IN), OPTIONAL, DIMENSION(:,:,:) ::  values_intwp_3d  !< output variable

    LOGICAL, INTENT(IN) ::  is_global  !< true if variable is global (same on all PE)

    REAL(KIND=4), POINTER,             INTENT(IN), OPTIONAL                   ::  values_real32_0d  !< output variable
    REAL(KIND=8), POINTER,             INTENT(IN), OPTIONAL                   ::  values_real64_0d  !< output variable
    REAL(wp),     POINTER,             INTENT(IN), OPTIONAL                   ::  values_realwp_0d  !< output variable
    REAL(KIND=4), POINTER, CONTIGUOUS, INTENT(IN), OPTIONAL, DIMENSION(:)     ::  values_real32_1d  !< output variable
    REAL(KIND=8), POINTER, CONTIGUOUS, INTENT(IN), OPTIONAL, DIMENSION(:)     ::  values_real64_1d  !< output variable
    REAL(wp),     POINTER, CONTIGUOUS, INTENT(IN), OPTIONAL, DIMENSION(:)     ::  values_realwp_1d  !< output variable
    REAL(KIND=4), POINTER, CONTIGUOUS, INTENT(IN), OPTIONAL, DIMENSION(:,:)   ::  values_real32_2d  !< output variable
    REAL(KIND=8), POINTER, CONTIGUOUS, INTENT(IN), OPTIONAL, DIMENSION(:,:)   ::  values_real64_2d  !< output variable
    REAL(wp),     POINTER, CONTIGUOUS, INTENT(IN), OPTIONAL, DIMENSION(:,:)   ::  values_realwp_2d  !< output variable
    REAL(KIND=4), POINTER, CONTIGUOUS, INTENT(IN), OPTIONAL, DIMENSION(:,:,:) ::  values_real32_3d  !< output variable
    REAL(KIND=8), POINTER, CONTIGUOUS, INTENT(IN), OPTIONAL, DIMENSION(:,:,:) ::  values_real64_3d  !< output variable
    REAL(wp),     POINTER, CONTIGUOUS, INTENT(IN), OPTIONAL, DIMENSION(:,:,:) ::  values_realwp_3d  !< output variable


    return_value = 0

    WRITE( temp_string, '(": write variable ",I6," into file ",I6)' ) variable_id, file_id
    CALL internal_message( 'debug', routine_name // TRIM( temp_string ) )

    IF ( is_global )  CONTINUE  ! reqired to prevent compiler warning

    IF ( .NOT. ANY( value_counts == 0 ) )  THEN
       WRITE( file_id )  variable_id
       WRITE( file_id )  bounds_start
       WRITE( file_id )  value_counts
       WRITE( file_id )  bounds_origin
!
!--    character output
       IF ( PRESENT( values_char_0d ) )  THEN
          output_string = 'char'
          WRITE( file_id )  output_string
          WRITE( file_id )  values_char_0d
       ELSEIF ( PRESENT( values_char_1d ) )  THEN
          output_string = 'char'
          WRITE( file_id )  output_string
          WRITE( file_id )  values_char_1d
       ELSEIF ( PRESENT( values_char_2d ) )  THEN
          output_string = 'char'
          WRITE( file_id )  output_string
          WRITE( file_id )  values_char_2d
       ELSEIF ( PRESENT( values_char_3d ) )  THEN
          output_string = 'char'
          WRITE( file_id )  output_string
          WRITE( file_id )  values_char_3d
!
!--    8bit integer output
       ELSEIF ( PRESENT( values_int8_0d ) )  THEN
          output_string = 'int8'
          WRITE( file_id )  output_string
          WRITE( file_id )  values_int8_0d
       ELSEIF ( PRESENT( values_int8_1d ) )  THEN
          output_string = 'int8'
          WRITE( file_id )  output_string
          WRITE( file_id )  values_int8_1d
       ELSEIF ( PRESENT( values_int8_2d ) )  THEN
          output_string = 'int8'
          WRITE( file_id )  output_string
          WRITE( file_id )  values_int8_2d
       ELSEIF ( PRESENT( values_int8_3d ) )  THEN
          output_string = 'int8'
          WRITE( file_id )  output_string
          WRITE( file_id )  values_int8_3d
!
!--    16bit integer output
       ELSEIF ( PRESENT( values_int16_0d ) )  THEN
          output_string = 'int16'
          WRITE( file_id )  output_string
          WRITE( file_id )  values_int16_0d
       ELSEIF ( PRESENT( values_int16_1d ) )  THEN
          output_string = 'int16'
          WRITE( file_id )  output_string
          WRITE( file_id )  values_int16_1d
       ELSEIF ( PRESENT( values_int16_2d ) )  THEN
          output_string = 'int16'
          WRITE( file_id )  output_string
          WRITE( file_id )  values_int16_2d
       ELSEIF ( PRESENT( values_int16_3d ) )  THEN
          output_string = 'int16'
          WRITE( file_id )  output_string
          WRITE( file_id )  values_int16_3d
!
!--    32bit integer output
       ELSEIF ( PRESENT( values_int32_0d ) )  THEN
          output_string = 'int32'
          WRITE( file_id )  output_string
          WRITE( file_id )  values_int32_0d
       ELSEIF ( PRESENT( values_int32_1d ) )  THEN
          output_string = 'int32'
          WRITE( file_id )  output_string
          WRITE( file_id )  values_int32_1d
       ELSEIF ( PRESENT( values_int32_2d ) )  THEN
          output_string = 'int32'
          WRITE( file_id )  output_string
          WRITE( file_id )  values_int32_2d
       ELSEIF ( PRESENT( values_int32_3d ) )  THEN
          output_string = 'int32'
          WRITE( file_id )  output_string
          WRITE( file_id )  values_int32_3d
!
!--    working-precision integer output
       ELSEIF ( PRESENT( values_intwp_0d ) )  THEN
          output_string = 'intwp'
          WRITE( file_id )  output_string
          WRITE( file_id )  values_intwp_0d
       ELSEIF ( PRESENT( values_intwp_1d ) )  THEN
          output_string = 'intwp'
          WRITE( file_id )  output_string
          WRITE( file_id )  values_intwp_1d
       ELSEIF ( PRESENT( values_intwp_2d ) )  THEN
          output_string = 'intwp'
          WRITE( file_id )  output_string
          WRITE( file_id )  values_intwp_2d
       ELSEIF ( PRESENT( values_intwp_3d ) )  THEN
          output_string = 'intwp'
          WRITE( file_id )  output_string
          WRITE( file_id )  values_intwp_3d
!
!--    32bit real output
       ELSEIF ( PRESENT( values_real32_0d ) )  THEN
          output_string = 'real32'
          WRITE( file_id )  output_string
          WRITE( file_id )  values_real32_0d
       ELSEIF ( PRESENT( values_real32_1d ) )  THEN
          output_string = 'real32'
          WRITE( file_id )  output_string
          WRITE( file_id )  values_real32_1d
       ELSEIF ( PRESENT( values_real32_2d ) )  THEN
          output_string = 'real32'
          WRITE( file_id )  output_string
          WRITE( file_id )  values_real32_2d
       ELSEIF ( PRESENT( values_real32_3d ) )  THEN
          output_string = 'real32'
          WRITE( file_id )  output_string
          WRITE( file_id )  values_real32_3d
!
!--    64bit real output
       ELSEIF ( PRESENT( values_real64_0d ) )  THEN
          output_string = 'real64'
          WRITE( file_id )  output_string
          WRITE( file_id )  values_real64_0d
       ELSEIF ( PRESENT( values_real64_1d ) )  THEN
          output_string = 'real64'
          WRITE( file_id )  output_string
          WRITE( file_id )  values_real64_1d
       ELSEIF ( PRESENT( values_real64_2d ) )  THEN
          output_string = 'real64'
          WRITE( file_id )  output_string
          WRITE( file_id )  values_real64_2d
       ELSEIF ( PRESENT( values_real64_3d ) )  THEN
          output_string = 'real64'
          WRITE( file_id )  output_string
          WRITE( file_id )  values_real64_3d
!
!--    working-precision real output
       ELSEIF ( PRESENT( values_realwp_0d ) )  THEN
          output_string = 'realwp'
          WRITE( file_id )  output_string
          WRITE( file_id )  values_realwp_0d
       ELSEIF ( PRESENT( values_realwp_1d ) )  THEN
          output_string = 'realwp'
          WRITE( file_id )  output_string
          WRITE( file_id )  values_realwp_1d
       ELSEIF ( PRESENT( values_realwp_2d ) )  THEN
          output_string = 'realwp'
          WRITE( file_id )  output_string
          WRITE( file_id )  values_realwp_2d
       ELSEIF ( PRESENT( values_realwp_3d ) )  THEN
          output_string = 'realwp'
          WRITE( file_id )  output_string
          WRITE( file_id )  values_realwp_3d
       ELSE
          return_value = 1
          CALL internal_message( 'error', routine_name // ': no values given' )
       ENDIF

    ENDIF

 END SUBROUTINE binary_write_variable

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Close opened files.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE binary_finalize( file_id, return_value )

    CHARACTER(LEN=charlen) ::  output_string  !< output string

    CHARACTER(LEN=*), PARAMETER ::  routine_name = 'binary_finalize'  !< name of this routine

    INTEGER, INTENT(IN)  ::  file_id       !< file ID
    INTEGER, INTENT(OUT) ::  return_value  !< return value


    IF ( config_file_open )  THEN

       output_string = '*** end config file ***'
       WRITE( config_file_unit )  output_string

       CLOSE( config_file_unit, IOSTAT=return_value )

       IF ( return_value /= 0 )  THEN
          CALL internal_message( 'error', routine_name // ': cannot close configuration file' )
       ELSE
          config_file_open = .FALSE.
       ENDIF

    ELSE

       return_value = 0

    ENDIF

    IF ( return_value == 0 )  THEN

       WRITE( temp_string, * ) file_id
       CALL internal_message( 'debug', routine_name // &
                              ': close file (file_id=' // TRIM( temp_string ) // ')' )

       CLOSE( file_id, IOSTAT=return_value )
       IF ( return_value /= 0 )  THEN
          WRITE( temp_string, * ) file_id
          CALL internal_message( 'error', routine_name // &
                                 ': cannot close file (file_id=' // TRIM( temp_string ) // ')' )
       ENDIF

    ENDIF

 END SUBROUTINE binary_finalize

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Message routine writing debug information into the debug file
!> or creating the error message string.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE internal_message( level, string )

    CHARACTER(LEN=*), INTENT(IN) ::  level   !< message importance level
    CHARACTER(LEN=*), INTENT(IN) ::  string  !< message string


    IF ( TRIM( level ) == 'error' )  THEN

       WRITE( internal_error_message, '(A,A)' ) ': ', string

    ELSEIF ( TRIM( level ) == 'debug'  .AND.  print_debug_output )  THEN

       WRITE( debug_output_unit, '(A,A)' ) 'DOM DEBUG: ', string
       FLUSH( debug_output_unit )

    ENDIF

 END SUBROUTINE internal_message

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Return the last created error message.
!--------------------------------------------------------------------------------------------------!
 FUNCTION binary_get_error_message() RESULT( error_message )

    CHARACTER(LEN=800) ::  error_message  !< return error message to main program


    error_message = TRIM( internal_error_message )

    internal_error_message = ''

 END FUNCTION binary_get_error_message

 END MODULE data_output_binary_module
