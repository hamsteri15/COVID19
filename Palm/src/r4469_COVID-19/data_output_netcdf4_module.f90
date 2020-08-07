!> @file data_output_netcdf4_module.f90
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
! $Id: data_output_netcdf4_module.f90 4429 2020-02-27 15:24:30Z raasch $
! bugfix: cpp-directive moved to avoid compile error due to unused dummy argument
! 
! 4408 2020-02-14 10:04:39Z gronemeier
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
!> NetCDF output module to write data to NetCDF files.
!> This is either done in parallel mode via parallel NetCDF4 I/O or in serial mode only by PE0.
!--------------------------------------------------------------------------------------------------!
 MODULE data_output_netcdf4_module

    USE kinds

#if defined( __parallel ) && !defined( __mpifh )
    USE MPI
#endif

#if defined( __netcdf4 )
    USE NETCDF
#endif

    IMPLICIT NONE

#if defined( __parallel ) && defined( __mpifh )
    INCLUDE "mpif.h"
#endif

    CHARACTER(LEN=800) ::  internal_error_message = ''  !< string containing the last error message
    CHARACTER(LEN=100) ::  file_suffix = ''             !< file suffix added to each file name
    CHARACTER(LEN=800) ::  temp_string                  !< dummy string

    CHARACTER(LEN=*), PARAMETER ::  mode_parallel = 'parallel'  !< string selecting netcdf4 parallel mode
    CHARACTER(LEN=*), PARAMETER ::  mode_serial   = 'serial'    !< string selecting netcdf4 serial mode

    INTEGER ::  debug_output_unit       !< Fortran Unit Number of the debug-output file
    INTEGER ::  global_id_in_file = -1  !< value of global ID within a file
    INTEGER ::  master_rank             !< master rank for tasks to be executed by single PE only
    INTEGER ::  output_group_comm       !< MPI communicator addressing all MPI ranks which participate in output

    LOGICAL ::  print_debug_output = .FALSE.  !< if true, debug output is printed

    SAVE

    PRIVATE

    INTERFACE netcdf4_init_module
       MODULE PROCEDURE netcdf4_init_module
    END INTERFACE netcdf4_init_module

    INTERFACE netcdf4_open_file
       MODULE PROCEDURE netcdf4_open_file
    END INTERFACE netcdf4_open_file

    INTERFACE netcdf4_init_dimension
       MODULE PROCEDURE netcdf4_init_dimension
    END INTERFACE netcdf4_init_dimension

    INTERFACE netcdf4_init_variable
       MODULE PROCEDURE netcdf4_init_variable
    END INTERFACE netcdf4_init_variable

    INTERFACE netcdf4_write_attribute
       MODULE PROCEDURE netcdf4_write_attribute
    END INTERFACE netcdf4_write_attribute

    INTERFACE netcdf4_stop_file_header_definition
       MODULE PROCEDURE netcdf4_stop_file_header_definition
    END INTERFACE netcdf4_stop_file_header_definition

    INTERFACE netcdf4_write_variable
       MODULE PROCEDURE netcdf4_write_variable
    END INTERFACE netcdf4_write_variable

    INTERFACE netcdf4_finalize
       MODULE PROCEDURE netcdf4_finalize
    END INTERFACE netcdf4_finalize

    INTERFACE netcdf4_get_error_message
       MODULE PROCEDURE netcdf4_get_error_message
    END INTERFACE netcdf4_get_error_message

    PUBLIC &
       netcdf4_finalize, &
       netcdf4_get_error_message, &
       netcdf4_init_dimension, &
       netcdf4_stop_file_header_definition, &
       netcdf4_init_module, &
       netcdf4_init_variable, &
       netcdf4_open_file, &
       netcdf4_write_attribute, &
       netcdf4_write_variable


 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize data-output module.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE netcdf4_init_module( file_suffix_of_output_group, mpi_comm_of_output_group, &
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

 END SUBROUTINE netcdf4_init_module

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Open netcdf file.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE netcdf4_open_file( mode, file_name, file_id, return_value )

    CHARACTER(LEN=*), INTENT(IN) ::  file_name  !< name of file
    CHARACTER(LEN=*), INTENT(IN) ::  mode       !< operation mode (either parallel or serial)

    CHARACTER(LEN=*), PARAMETER ::  routine_name = 'netcdf4_open_file'  !< name of this routine

    INTEGER, INTENT(OUT) ::  file_id       !< file ID
    INTEGER              ::  my_rank       !< MPI rank of processor
    INTEGER              ::  nc_stat       !< netcdf return value
    INTEGER, INTENT(OUT) ::  return_value  !< return value


    return_value = 0
    file_id = -1
!
!-- Open new file
    CALL internal_message( 'debug', routine_name // ': create file "' // TRIM( file_name ) // '"' )

    IF ( TRIM( mode ) == mode_serial )  THEN

#if defined( __netcdf4 )
#if defined( __parallel )
       CALL MPI_COMM_RANK( output_group_comm, my_rank, return_value )
       IF ( return_value /= 0 )  THEN
          CALL internal_message( 'error', routine_name // ': MPI_COMM_RANK error' )
       ENDIF
       IF ( my_rank /= master_rank )  THEN
          return_value = 1
          CALL internal_message( 'error', routine_name // &
                                 ': trying to define a NetCDF file in serial mode by an MPI ' // &
                                 'rank other than the master output rank. Serial NetCDF ' // &
                                 'files can only be defined by the master output rank!' )
       ENDIF
#else
       my_rank = master_rank
       return_value = 0
#endif

       IF ( return_value == 0 )  &
          nc_stat = NF90_CREATE( TRIM( file_name ) // TRIM( file_suffix ), &
                                 IOR( NF90_NOCLOBBER, NF90_NETCDF4 ),      &
                                 file_id )
#else
       nc_stat = 0
       return_value = 1
       CALL internal_message( 'error', routine_name //                               &
                              ': pre-processor directive "__netcdf4" not given. ' // &
                              'Using NetCDF4 output not possible' )
#endif

    ELSEIF ( TRIM( mode ) == mode_parallel )  THEN

#if defined( __parallel ) && defined( __netcdf4 ) && defined( __netcdf4_parallel )
       nc_stat = NF90_CREATE( TRIM( file_name ) // TRIM( file_suffix ),               &
                              IOR( NF90_NOCLOBBER, IOR( NF90_NETCDF4, NF90_MPIIO ) ), &
                              file_id, COMM = output_group_comm, INFO = MPI_INFO_NULL )
#else
       nc_stat = 0
       return_value = 1
       CALL internal_message( 'error', routine_name //                                 &
                              ': pre-processor directives "__parallel" and/or ' //     &
                              '"__netcdf4" and/or "__netcdf4_parallel" not given. ' // &
                              'Using parallel NetCDF4 output not possible' )
#endif

    ELSE
       nc_stat = 0
       return_value = 1
       CALL internal_message( 'error', routine_name // ': selected mode "' //  &
                                       TRIM( mode ) // '" must be either "' // &
                                       mode_serial // '" or "' // mode_parallel // '"' )
    ENDIF

#if defined( __netcdf4 )
    IF ( nc_stat /= NF90_NOERR  .AND.  return_value == 0 )  THEN
       return_value = 1
       CALL internal_message( 'error', routine_name //                 &
                              ': NetCDF error while opening file "' // &
                              TRIM( file_name ) // '": ' // NF90_STRERROR( nc_stat ) )
    ENDIF
#endif

 END SUBROUTINE netcdf4_open_file

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write attribute to netcdf file.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE netcdf4_write_attribute( file_id, variable_id, attribute_name, &
                  value_char, value_int8, value_int16, value_int32,        &
                  value_real32, value_real64, return_value )

    CHARACTER(LEN=*), INTENT(IN)           ::  attribute_name  !< name of attribute
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL ::  value_char      !< value of attribute

    CHARACTER(LEN=*), PARAMETER ::  routine_name = 'netcdf4_write_attribute'  !< name of this routine

    INTEGER ::  nc_stat    !< netcdf return value
    INTEGER ::  target_id  !< ID of target which gets attribute (either global or variable_id)

    INTEGER, INTENT(IN)  ::  file_id       !< file ID
    INTEGER, INTENT(OUT) ::  return_value  !< return value
    INTEGER, INTENT(IN)  ::  variable_id   !< variable ID

    INTEGER(KIND=1), INTENT(IN), OPTIONAL ::  value_int8   !< value of attribute
    INTEGER(KIND=2), INTENT(IN), OPTIONAL ::  value_int16  !< value of attribute
    INTEGER(KIND=4), INTENT(IN), OPTIONAL ::  value_int32  !< value of attribute

    REAL(KIND=4), INTENT(IN), OPTIONAL ::  value_real32  !< value of attribute
    REAL(KIND=8), INTENT(IN), OPTIONAL ::  value_real64  !< value of attribute


#if defined( __netcdf4 )
    return_value = 0

    IF ( variable_id == global_id_in_file )  THEN
       target_id = NF90_GLOBAL
    ELSE
       target_id = variable_id
    ENDIF

    CALL internal_message( 'debug', routine_name // &
                           ': write attribute "' // TRIM( attribute_name ) // '"' )

    IF ( PRESENT( value_char ) )  THEN
       nc_stat = NF90_PUT_ATT( file_id, target_id, TRIM( attribute_name ), TRIM( value_char ) )
    ELSEIF ( PRESENT( value_int8 ) )  THEN
       nc_stat = NF90_PUT_ATT( file_id, target_id, TRIM( attribute_name ), value_int8 )
    ELSEIF ( PRESENT( value_int16 ) )  THEN
       nc_stat = NF90_PUT_ATT( file_id, target_id, TRIM( attribute_name ), value_int16 )
    ELSEIF ( PRESENT( value_int32 ) )  THEN
       nc_stat = NF90_PUT_ATT( file_id, target_id, TRIM( attribute_name ), value_int32 )
    ELSEIF ( PRESENT( value_real32 ) )  THEN
       nc_stat = NF90_PUT_ATT( file_id, target_id, TRIM( attribute_name ), value_real32 )
    ELSEIF ( PRESENT( value_real64 ) )  THEN
       nc_stat = NF90_PUT_ATT( file_id, target_id, TRIM( attribute_name ), value_real64 )
    ELSE
       return_value = 1
       CALL internal_message( 'error', routine_name // &
                              ': no value given for attribute "' // TRIM( attribute_name ) // '"' )
    ENDIF

    IF ( return_value == 0 )  THEN
       IF ( nc_stat /= NF90_NOERR )  THEN
          return_value = 1
          CALL internal_message( 'error', routine_name //                      &
                                 ': NetCDF error while writing attribute "' // &
                                 TRIM( attribute_name ) // '": ' // NF90_STRERROR( nc_stat ) )
       ENDIF
    ENDIF
#else
    return_value = 1
#endif

 END SUBROUTINE netcdf4_write_attribute

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize dimension.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE netcdf4_init_dimension( mode, file_id, dimension_id, variable_id, &
               dimension_name, dimension_type, dimension_length, return_value )

    CHARACTER(LEN=*), INTENT(IN) ::  dimension_name  !< name of dimension
    CHARACTER(LEN=*), INTENT(IN) ::  dimension_type  !< data type of dimension
    CHARACTER(LEN=*), INTENT(IN) ::  mode            !< operation mode (either parallel or serial)

    CHARACTER(LEN=*), PARAMETER ::  routine_name = 'netcdf4_init_dimension'  !< name of this routine

    INTEGER, INTENT(OUT) ::  dimension_id         !< dimension ID
    INTEGER, INTENT(IN)  ::  dimension_length     !< length of dimension
    INTEGER, INTENT(IN)  ::  file_id              !< file ID
    INTEGER              ::  nc_dimension_length  !< length of dimension
    INTEGER              ::  nc_stat              !< netcdf return value
    INTEGER, INTENT(OUT) ::  return_value         !< return value
    INTEGER, INTENT(OUT) ::  variable_id          !< variable ID


#if defined( __netcdf4 )
    return_value = 0
    variable_id = -1

    CALL internal_message( 'debug', routine_name // &
                           ': init dimension "' // TRIM( dimension_name ) // '"' )
!
!-- Check if dimension is unlimited
    IF ( dimension_length < 0 )  THEN
       nc_dimension_length = NF90_UNLIMITED
    ELSE
       nc_dimension_length = dimension_length
    ENDIF
!
!-- Define dimension in file
    nc_stat = NF90_DEF_DIM( file_id, dimension_name, nc_dimension_length, dimension_id )

    IF ( nc_stat == NF90_NOERR )  THEN
!
!--    Define variable holding dimension values in file
       CALL netcdf4_init_variable( mode, file_id, variable_id, dimension_name, dimension_type, &
                                   (/ dimension_id /), is_global=.TRUE., return_value=return_value )

    ELSE
       return_value = 1
       CALL internal_message( 'error', routine_name //                           &
                              ': NetCDF error while initializing dimension "' // &
                              TRIM( dimension_name ) // '": ' // NF90_STRERROR( nc_stat ) )
    ENDIF
#else
    return_value = 1
    variable_id = -1
    dimension_id = -1
#endif

 END SUBROUTINE netcdf4_init_dimension

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize variable.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE netcdf4_init_variable( mode, file_id, variable_id, variable_name, variable_type, &
                                   dimension_ids, is_global, return_value )

    CHARACTER(LEN=*), INTENT(IN) ::  mode           !< operation mode (either parallel or serial)
    CHARACTER(LEN=*), INTENT(IN) ::  variable_name  !< name of variable
    CHARACTER(LEN=*), INTENT(IN) ::  variable_type  !< data type of variable

    CHARACTER(LEN=*), PARAMETER ::  routine_name = 'netcdf4_init_variable'  !< name of this routine

    INTEGER, INTENT(IN)  ::  file_id           !< file ID
    INTEGER              ::  nc_stat           !< netcdf return value
    INTEGER              ::  nc_variable_type  !< netcdf data type
    INTEGER, INTENT(OUT) ::  return_value      !< return value
    INTEGER, INTENT(OUT) ::  variable_id       !< variable ID

    INTEGER, DIMENSION(:), INTENT(IN) ::  dimension_ids  !< list of dimension IDs used by variable

    LOGICAL, INTENT(IN)  ::  is_global  !< true if variable is global (same on all PE)


#if defined( __netcdf4 )
    return_value = 0

    WRITE( temp_string, * ) is_global
    CALL internal_message( 'debug', routine_name //                        &
                           ': init variable "' // TRIM( variable_name ) // &
                           '" ( is_global = ' // TRIM( temp_string ) // ')' )

    nc_variable_type = get_netcdf_data_type( variable_type )

    IF ( nc_variable_type /= -1 )  THEN
!
!--    Define variable in file
       nc_stat = NF90_DEF_VAR( file_id, variable_name, nc_variable_type, dimension_ids, variable_id )

!
!--    Define how variable can be accessed by PEs in parallel netcdf file
       IF ( nc_stat == NF90_NOERR  .AND.  TRIM( mode ) == mode_parallel )  THEN
#if defined( __netcdf4_parallel )
          IF ( is_global )  THEN
             nc_stat = NF90_VAR_PAR_ACCESS( file_id, variable_id, NF90_INDEPENDENT )
          ELSE
             nc_stat = NF90_VAR_PAR_ACCESS( file_id, variable_id, NF90_COLLECTIVE )
          ENDIF
#else
          CONTINUE
#endif
       ENDIF

       IF ( nc_stat /= NF90_NOERR )  THEN
          return_value = 1
          CALL internal_message( 'error', routine_name //                          &
                                 ': NetCDF error while initializing variable "' // &
                                 TRIM( variable_name ) // '": ' // NF90_STRERROR( nc_stat ) )
       ENDIF

    ELSE
       return_value = 1
    ENDIF

#else
    return_value = 1
    variable_id = -1
#endif

 END SUBROUTINE netcdf4_init_variable

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Leave file definition state.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE netcdf4_stop_file_header_definition( file_id, return_value )

    CHARACTER(LEN=*), PARAMETER ::  routine_name = 'netcdf4_stop_file_header_definition'  !< name of this routine

    INTEGER, INTENT(IN)  ::  file_id        !< file ID
    INTEGER              ::  nc_stat        !< netcdf return value
    INTEGER              ::  old_fill_mode  !< previous netcdf fill mode
    INTEGER, INTENT(OUT) ::  return_value   !< return value


#if defined( __netcdf4 )
    return_value = 0

    WRITE( temp_string, * ) file_id
    CALL internal_message( 'debug', routine_name // &
                           ': finalize file definition (file_id=' // TRIM( temp_string ) // ')' )
!
!-- Set general no fill, otherwise the performance drops significantly
    nc_stat = NF90_SET_FILL( file_id, NF90_NOFILL, old_fill_mode )

    IF ( nc_stat == NF90_NOERR )  THEN
       nc_stat = NF90_ENDDEF( file_id )
    ENDIF

    IF ( nc_stat /= NF90_NOERR )  THEN
       return_value = 1
       CALL internal_message( 'error', routine_name // &
                              ': NetCDF error: ' // NF90_STRERROR( nc_stat ) )
    ENDIF
#else
    return_value = 1
#endif

 END SUBROUTINE netcdf4_stop_file_header_definition

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write variable of different kind into netcdf file.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE netcdf4_write_variable(                                                    &
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

    CHARACTER(LEN=*), PARAMETER ::  routine_name = 'netcdf4_write_variable'  !< name of this routine

    CHARACTER(LEN=1), POINTER,             INTENT(IN), OPTIONAL                   ::  values_char_0d  !< output variable
    CHARACTER(LEN=1), POINTER, CONTIGUOUS, INTENT(IN), OPTIONAL, DIMENSION(:)     ::  values_char_1d  !< output variable
    CHARACTER(LEN=1), POINTER, CONTIGUOUS, INTENT(IN), OPTIONAL, DIMENSION(:,:)   ::  values_char_2d  !< output variable
    CHARACTER(LEN=1), POINTER, CONTIGUOUS, INTENT(IN), OPTIONAL, DIMENSION(:,:,:) ::  values_char_3d  !< output variable

    INTEGER              ::  d             !< loop index
    INTEGER, INTENT(IN)  ::  file_id       !< file ID
    INTEGER              ::  my_rank       !< MPI rank of processor
    INTEGER              ::  nc_stat       !< netcdf return value
    INTEGER              ::  ndims         !< number of dimensions of variable in file
    INTEGER, INTENT(OUT) ::  return_value  !< return value
    INTEGER, INTENT(IN)  ::  variable_id   !< variable ID

    INTEGER, DIMENSION(:),              INTENT(IN)  ::  bounds_origin      !< starting index of each dimension
    INTEGER, DIMENSION(:),              INTENT(IN)  ::  bounds_start       !< starting index of variable
    INTEGER, DIMENSION(:), ALLOCATABLE              ::  dimension_ids      !< IDs of dimensions of variable in file
    INTEGER, DIMENSION(:), ALLOCATABLE              ::  dimension_lengths  !< length of dimensions of variable in file
    INTEGER, DIMENSION(:),              INTENT(IN)  ::  value_counts       !< count of values along each dimension to be written

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


#if defined( __netcdf4 )

#if defined( __parallel )
    CALL MPI_COMM_RANK( output_group_comm, my_rank, return_value )
    IF ( return_value /= 0 )  THEN
       CALL internal_message( 'error', routine_name // ': MPI_COMM_RANK error' )
    ENDIF
#else
    my_rank = master_rank
    return_value = 0
#endif

    IF ( return_value == 0  .AND.  ( .NOT. is_global  .OR.  my_rank == master_rank ) )  THEN

       WRITE( temp_string, * ) variable_id
       CALL internal_message( 'debug', routine_name // ': write variable ' // TRIM( temp_string ) )

       ndims = SIZE( bounds_start )

!
!--    character output
       IF ( PRESENT( values_char_0d ) )  THEN
          nc_stat = NF90_PUT_VAR( file_id, variable_id, (/ values_char_0d /), &
                                  start = bounds_start - bounds_origin + 1,   &
                                  count = value_counts )
       ELSEIF ( PRESENT( values_char_1d ) )  THEN
          nc_stat = NF90_PUT_VAR( file_id, variable_id, values_char_1d,     &
                                  start = bounds_start - bounds_origin + 1, &
                                  count = value_counts )
       ELSEIF ( PRESENT( values_char_2d ) )  THEN
          nc_stat = NF90_PUT_VAR( file_id, variable_id, values_char_2d,     &
                                  start = bounds_start - bounds_origin + 1, &
                                  count = value_counts )
       ELSEIF ( PRESENT( values_char_3d ) )  THEN
          nc_stat = NF90_PUT_VAR( file_id, variable_id, values_char_3d,     &
                                  start = bounds_start - bounds_origin + 1, &
                                  count = value_counts )
!
!--    8bit integer output
       ELSEIF ( PRESENT( values_int8_0d ) )  THEN
          nc_stat = NF90_PUT_VAR( file_id, variable_id, (/ values_int8_0d /), &
                                  start = bounds_start - bounds_origin + 1,   &
                                  count = value_counts )
       ELSEIF ( PRESENT( values_int8_1d ) )  THEN
          nc_stat = NF90_PUT_VAR( file_id, variable_id, values_int8_1d,     &
                                  start = bounds_start - bounds_origin + 1, &
                                  count = value_counts )
       ELSEIF ( PRESENT( values_int8_2d ) )  THEN
          nc_stat = NF90_PUT_VAR( file_id, variable_id, values_int8_2d,     &
                                  start = bounds_start - bounds_origin + 1, &
                                  count = value_counts )
       ELSEIF ( PRESENT( values_int8_3d ) )  THEN
          nc_stat = NF90_PUT_VAR( file_id, variable_id, values_int8_3d,     &
                                  start = bounds_start - bounds_origin + 1, &
                                  count = value_counts )
!
!--    16bit integer output
       ELSEIF ( PRESENT( values_int16_0d ) )  THEN
          nc_stat = NF90_PUT_VAR( file_id, variable_id, (/ values_int16_0d /), &
                                  start = bounds_start - bounds_origin + 1,    &
                                  count = value_counts )
       ELSEIF ( PRESENT( values_int16_1d ) )  THEN
          nc_stat = NF90_PUT_VAR( file_id, variable_id, values_int16_1d,    &
                                  start = bounds_start - bounds_origin + 1, &
                                  count = value_counts )
       ELSEIF ( PRESENT( values_int16_2d ) )  THEN
          nc_stat = NF90_PUT_VAR( file_id, variable_id, values_int16_2d,    &
                                  start = bounds_start - bounds_origin + 1, &
                                  count = value_counts )
       ELSEIF ( PRESENT( values_int16_3d ) )  THEN
          nc_stat = NF90_PUT_VAR( file_id, variable_id, values_int16_3d,    &
                                  start = bounds_start - bounds_origin + 1, &
                                  count = value_counts )
!
!--    32bit integer output
       ELSEIF ( PRESENT( values_int32_0d ) )  THEN
          nc_stat = NF90_PUT_VAR( file_id, variable_id, (/ values_int32_0d /),  &
                                  start = bounds_start - bounds_origin + 1,     &
                                  count = value_counts )
       ELSEIF ( PRESENT( values_int32_1d ) )  THEN
          nc_stat = NF90_PUT_VAR( file_id, variable_id, values_int32_1d,    &
                                  start = bounds_start - bounds_origin + 1, &
                                  count = value_counts )
       ELSEIF ( PRESENT( values_int32_2d ) )  THEN
          nc_stat = NF90_PUT_VAR( file_id, variable_id, values_int32_2d,    &
                                  start = bounds_start - bounds_origin + 1, &
                                  count = value_counts )
       ELSEIF ( PRESENT( values_int32_3d ) )  THEN
          nc_stat = NF90_PUT_VAR( file_id, variable_id, values_int32_3d,    &
                                  start = bounds_start - bounds_origin + 1, &
                                  count = value_counts )
!
!--    working-precision integer output
       ELSEIF ( PRESENT( values_intwp_0d ) )  THEN
          nc_stat = NF90_PUT_VAR( file_id, variable_id, (/ values_intwp_0d /),  &
                                  start = bounds_start - bounds_origin + 1,     &
                                  count = value_counts )
       ELSEIF ( PRESENT( values_intwp_1d ) )  THEN
          nc_stat = NF90_PUT_VAR( file_id, variable_id, values_intwp_1d,    &
                                  start = bounds_start - bounds_origin + 1, &
                                  count = value_counts )
       ELSEIF ( PRESENT( values_intwp_2d ) )  THEN
          nc_stat = NF90_PUT_VAR( file_id, variable_id, values_intwp_2d,    &
                                  start = bounds_start - bounds_origin + 1, &
                                  count = value_counts )
       ELSEIF ( PRESENT( values_intwp_3d ) )  THEN
          nc_stat = NF90_PUT_VAR( file_id, variable_id, values_intwp_3d,    &
                                  start = bounds_start - bounds_origin + 1, &
                                  count = value_counts )
!
!--    32bit real output
       ELSEIF ( PRESENT( values_real32_0d ) )  THEN
          nc_stat = NF90_PUT_VAR( file_id, variable_id, (/ values_real32_0d /), &
                                  start = bounds_start - bounds_origin + 1,     &
                                  count = value_counts )
       ELSEIF ( PRESENT( values_real32_1d ) )  THEN
          nc_stat = NF90_PUT_VAR( file_id, variable_id, values_real32_1d,   &
                                  start = bounds_start - bounds_origin + 1, &
                                  count = value_counts )
       ELSEIF ( PRESENT( values_real32_2d ) )  THEN
          nc_stat = NF90_PUT_VAR( file_id, variable_id, values_real32_2d,   &
                                  start = bounds_start - bounds_origin + 1, &
                                  count = value_counts )
       ELSEIF ( PRESENT( values_real32_3d ) )  THEN
          nc_stat = NF90_PUT_VAR( file_id, variable_id, values_real32_3d,   &
                                  start = bounds_start - bounds_origin + 1, &
                                  count = value_counts )
!
!--    64bit real output
       ELSEIF ( PRESENT( values_real64_0d ) )  THEN
          nc_stat = NF90_PUT_VAR( file_id, variable_id, (/ values_real64_0d /), &
                                  start = bounds_start - bounds_origin + 1,     &
                                  count = value_counts )
       ELSEIF ( PRESENT( values_real64_1d ) )  THEN
          nc_stat = NF90_PUT_VAR( file_id, variable_id, values_real64_1d,   &
                                  start = bounds_start - bounds_origin + 1, &
                                  count = value_counts )
       ELSEIF ( PRESENT( values_real64_2d ) )  THEN
          nc_stat = NF90_PUT_VAR( file_id, variable_id, values_real64_2d,   &
                                  start = bounds_start - bounds_origin + 1, &
                                  count = value_counts )
       ELSEIF ( PRESENT( values_real64_3d ) )  THEN
          nc_stat = NF90_PUT_VAR( file_id, variable_id, values_real64_3d,   &
                                  start = bounds_start - bounds_origin + 1, &
                                  count = value_counts )
!
!--    working-precision real output
       ELSEIF ( PRESENT( values_realwp_0d ) )  THEN
          nc_stat = NF90_PUT_VAR( file_id, variable_id, (/ values_realwp_0d /), &
                                  start = bounds_start - bounds_origin + 1,     &
                                  count = value_counts )
       ELSEIF ( PRESENT( values_realwp_1d ) )  THEN
          nc_stat = NF90_PUT_VAR( file_id, variable_id, values_realwp_1d,   &
                                  start = bounds_start - bounds_origin + 1, &
                                  count = value_counts )
       ELSEIF ( PRESENT( values_realwp_2d ) )  THEN
          nc_stat = NF90_PUT_VAR( file_id, variable_id, values_realwp_2d,   &
                                  start = bounds_start - bounds_origin + 1, &
                                  count = value_counts )
       ELSEIF ( PRESENT( values_realwp_3d ) )  THEN
          nc_stat = NF90_PUT_VAR( file_id, variable_id, values_realwp_3d,   &
                                  start = bounds_start - bounds_origin + 1, &
                                  count = value_counts )
       ELSE
          return_value = 1
          nc_stat = NF90_NOERR
          WRITE( temp_string, '(": variable_id=",I6,", file_id=",I6)' ) variable_id, file_id
          CALL internal_message( 'error', routine_name // &
                                 ': no output values given ' // TRIM( temp_string ) )
       ENDIF
!
!--    Check for errors
       IF ( nc_stat /= NF90_NOERR )  THEN
          return_value = 1

          IF ( nc_stat == NF90_EEDGE .OR. nc_stat == NF90_EINVALCOORDS )  THEN
!
!--          If given bounds exceed dimension bounds, get information of bounds in file
             WRITE( temp_string, * )  NF90_STRERROR( nc_stat )

             ALLOCATE( dimension_ids(ndims) )
             ALLOCATE( dimension_lengths(ndims) )

             nc_stat = NF90_INQUIRE_VARIABLE( file_id, variable_id, dimids=dimension_ids )

             d = 1
             DO WHILE ( d <= ndims .AND. nc_stat == NF90_NOERR )
                nc_stat = NF90_INQUIRE_DIMENSION( file_id, dimension_ids(d), &
                                                  LEN=dimension_lengths(d) )
                d = d + 1
             ENDDO

             IF ( nc_stat == NF90_NOERR )  THEN
                WRITE( temp_string, * )  TRIM( temp_string ) // '; given variable bounds: ' //  &
                   'start=', bounds_start, ', count=', value_counts, ', origin=', bounds_origin
                CALL internal_message( 'error', routine_name //     &
                                       ': error while writing: ' // TRIM( temp_string ) )
             ELSE
!
!--             Error occured during NF90_INQUIRE_VARIABLE or NF90_INQUIRE_DIMENSION
                CALL internal_message( 'error', routine_name //            &
                                       ': error while accessing file: ' // &
                                        NF90_STRERROR( nc_stat ) )
             ENDIF

          ELSE
!
!--          Other NetCDF error
             CALL internal_message( 'error', routine_name //     &
                                    ': error while writing: ' // NF90_STRERROR( nc_stat ) )
          ENDIF
       ENDIF

    ENDIF
#else
    return_value = 1
#endif

 END SUBROUTINE netcdf4_write_variable

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Close netcdf file.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE netcdf4_finalize( file_id, return_value )

    CHARACTER(LEN=*), PARAMETER ::  routine_name = 'netcdf4_finalize'  !< name of routine

    INTEGER, INTENT(IN)  ::  file_id       !< file ID
    INTEGER              ::  nc_stat       !< netcdf return value
    INTEGER, INTENT(OUT) ::  return_value  !< return value


#if defined( __netcdf4 )
    WRITE( temp_string, * ) file_id
    CALL internal_message( 'debug', routine_name // &
                           ': close file (file_id=' // TRIM( temp_string ) // ')' )

    nc_stat = NF90_CLOSE( file_id )
    IF ( nc_stat == NF90_NOERR )  THEN
       return_value = 0
    ELSE
       return_value = 1
       CALL internal_message( 'error', routine_name // &
                              ': NetCDF error: ' // NF90_STRERROR( nc_stat ) )
    ENDIF
#else
    return_value = 1
#endif

 END SUBROUTINE netcdf4_finalize

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Convert data_type string into netcdf data type value.
!--------------------------------------------------------------------------------------------------!
 FUNCTION get_netcdf_data_type( data_type ) RESULT( return_value )

    CHARACTER(LEN=*), INTENT(IN) ::  data_type  !< requested data type

    CHARACTER(LEN=*), PARAMETER ::  routine_name = 'get_netcdf_data_type'  !< name of this routine

    INTEGER ::  return_value  !< netcdf data type


    SELECT CASE ( TRIM( data_type ) )

#if defined( __netcdf4 )
       CASE ( 'char' )
          return_value = NF90_CHAR

       CASE ( 'int8' )
          return_value = NF90_BYTE

       CASE ( 'int16' )
          return_value = NF90_SHORT

       CASE ( 'int32' )
          return_value = NF90_INT

       CASE ( 'real32' )
          return_value = NF90_FLOAT

       CASE ( 'real64' )
          return_value = NF90_DOUBLE
#endif

       CASE DEFAULT
          CALL internal_message( 'error', routine_name // &
                                 ': data type unknown (' // TRIM( data_type ) // ')' )
          return_value = -1

    END SELECT

 END FUNCTION get_netcdf_data_type

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
 FUNCTION netcdf4_get_error_message() RESULT( error_message )

    CHARACTER(LEN=800) ::  error_message  !< return error message to main program


    error_message = TRIM( internal_error_message )

    internal_error_message = ''

 END FUNCTION netcdf4_get_error_message


 END MODULE data_output_netcdf4_module
