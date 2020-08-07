!> @file data_output_module.f90
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
! $Id: data_output_module.f90 4408 2020-02-14 10:04:39Z gronemeier $
! Enable character-array output
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
!> @author Tobias Gronemeier
!> @author Helge Knoop
!
!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Data-output module to handle output of variables into output files.
!>
!> The module first creates an interal database containing all meta data of all output quantities.
!> After defining all meta data, the output files are initialized and prepared for writing. When
!> writing is finished, files can be finalized and closed.
!> The order of calls are as follows:
!>   1. Initialize the module via
!>      'dom_init'
!>   2. Define output files via (multiple calls of)
!>      'dom_def_file', 'dom_def_att', 'dom_def_dim', 'dom_def_var'
!>   3. Leave definition stage via
!>      'dom_def_end'
!>   4. Write output data into file via
!>      'dom_write_var'
!>   5. Finalize the output via
!>      'dom_finalize_output'
!> If any routine exits with a non-zero return value, the error message of the last encountered
!> error can be fetched via 'dom_get_error_message'.
!> For debugging purposes, the content of the database can be written to the debug output via
!> 'dom_database_debug_output'.
!>
!> @todo Convert variable if type of given values do not fit specified type.
!--------------------------------------------------------------------------------------------------!
 MODULE data_output_module

    USE kinds

    USE data_output_netcdf4_module, &
       ONLY: netcdf4_init_dimension, &
             netcdf4_get_error_message, &
             netcdf4_stop_file_header_definition, &
             netcdf4_init_module, &
             netcdf4_init_variable, &
             netcdf4_finalize, &
             netcdf4_open_file, &
             netcdf4_write_attribute, &
             netcdf4_write_variable

    USE data_output_binary_module, &
       ONLY: binary_finalize, &
             binary_get_error_message, &
             binary_init_dimension, &
             binary_stop_file_header_definition, &
             binary_init_module, &
             binary_init_variable, &
             binary_open_file, &
             binary_write_attribute, &
             binary_write_variable

    IMPLICIT NONE

    INTEGER, PARAMETER ::  charlen = 100  !< maximum length of character variables
    INTEGER, PARAMETER ::  no_id = -1     !< default ID if no ID was assigned

    TYPE attribute_type
       CHARACTER(LEN=charlen) ::  data_type = ''  !< data type
       CHARACTER(LEN=charlen) ::  name            !< attribute name
       CHARACTER(LEN=charlen) ::  value_char      !< attribute value if character
       INTEGER(KIND=1)        ::  value_int8      !< attribute value if 8bit integer
       INTEGER(KIND=2)        ::  value_int16     !< attribute value if 16bit integer
       INTEGER(KIND=4)        ::  value_int32     !< attribute value if 32bit integer
       REAL(KIND=4)           ::  value_real32    !< attribute value if 32bit real
       REAL(KIND=8)           ::  value_real64    !< attribute value if 64bit real
    END TYPE attribute_type

    TYPE variable_type
       CHARACTER(LEN=charlen)                            ::  data_type = ''       !< data type
       CHARACTER(LEN=charlen)                            ::  name                 !< variable name
       INTEGER                                           ::  id = no_id           !< id within file
       LOGICAL                                           ::  is_global = .FALSE.  !< true if global variable
       CHARACTER(LEN=charlen), DIMENSION(:), ALLOCATABLE ::  dimension_names      !< list of dimension names used by variable
       INTEGER,                DIMENSION(:), ALLOCATABLE ::  dimension_ids        !< list of dimension ids used by variable
       TYPE(attribute_type),   DIMENSION(:), ALLOCATABLE ::  attributes           !< list of attributes
    END TYPE variable_type

    TYPE dimension_type
       CHARACTER(LEN=charlen)                     ::  data_type = ''        !< data type
       CHARACTER(LEN=charlen)                     ::  name                  !< dimension name
       INTEGER                                    ::  id = no_id            !< dimension id within file
       INTEGER                                    ::  length                !< length of dimension
       INTEGER                                    ::  length_mask           !< length of masked dimension
       INTEGER                                    ::  variable_id = no_id   !< associated variable id within file
       LOGICAL                                    ::  is_masked = .FALSE.   !< true if masked
       INTEGER,         DIMENSION(2)              ::  bounds                !< lower and upper bound of dimension
       INTEGER,         DIMENSION(:), ALLOCATABLE ::  masked_indices        !< list of masked indices of dimension
       INTEGER(KIND=1), DIMENSION(:), ALLOCATABLE ::  masked_values_int8    !< masked dimension values if 16bit integer
       INTEGER(KIND=2), DIMENSION(:), ALLOCATABLE ::  masked_values_int16   !< masked dimension values if 16bit integer
       INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE ::  masked_values_int32   !< masked dimension values if 32bit integer
       INTEGER(iwp),    DIMENSION(:), ALLOCATABLE ::  masked_values_intwp   !< masked dimension values if working-precision int
       INTEGER(KIND=1), DIMENSION(:), ALLOCATABLE ::  values_int8           !< dimension values if 16bit integer
       INTEGER(KIND=2), DIMENSION(:), ALLOCATABLE ::  values_int16          !< dimension values if 16bit integer
       INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE ::  values_int32          !< dimension values if 32bit integer
       INTEGER(iwp),    DIMENSION(:), ALLOCATABLE ::  values_intwp          !< dimension values if working-precision integer
       LOGICAL,         DIMENSION(:), ALLOCATABLE ::  mask                  !< mask
       REAL(KIND=4),    DIMENSION(:), ALLOCATABLE ::  masked_values_real32  !< masked dimension values if 32bit real
       REAL(KIND=8),    DIMENSION(:), ALLOCATABLE ::  masked_values_real64  !< masked dimension values if 64bit real
       REAL(wp),        DIMENSION(:), ALLOCATABLE ::  masked_values_realwp  !< masked dimension values if working-precision real
       REAL(KIND=4),    DIMENSION(:), ALLOCATABLE ::  values_real32         !< dimension values if 32bit real
       REAL(KIND=8),    DIMENSION(:), ALLOCATABLE ::  values_real64         !< dimension values if 64bit real
       REAL(wp),        DIMENSION(:), ALLOCATABLE ::  values_realwp         !< dimension values if working-precision real
       TYPE(attribute_type), DIMENSION(:), ALLOCATABLE ::  attributes       !< list of attributes
    END TYPE dimension_type

    TYPE file_type
       CHARACTER(LEN=charlen)                          ::  format = ''        !< file format
       CHARACTER(LEN=charlen)                          ::  name = ''          !< file name
       INTEGER                                         ::  id = no_id         !< id of file
       LOGICAL                                         ::  is_init = .FALSE.  !< true if initialized
       TYPE(attribute_type), DIMENSION(:), ALLOCATABLE ::  attributes         !< list of attributes
       TYPE(dimension_type), DIMENSION(:), ALLOCATABLE ::  dimensions         !< list of dimensions
       TYPE(variable_type),  DIMENSION(:), ALLOCATABLE ::  variables          !< list of variables
    END TYPE file_type


    CHARACTER(LEN=charlen) ::  output_file_suffix = ''      !< file suffix added to each file name
    CHARACTER(LEN=800)     ::  internal_error_message = ''  !< string containing the last error message
    CHARACTER(LEN=800)     ::  temp_string                  !< dummy string

    INTEGER ::  debug_output_unit  !< Fortran Unit Number of the debug-output file
    INTEGER ::  nfiles = 0         !< number of files
    INTEGER ::  master_rank = 0    !< master rank for tasks to be executed by single PE only
    INTEGER ::  output_group_comm  !< MPI communicator addressing all MPI ranks which participate in output

    LOGICAL ::  print_debug_output = .FALSE.  !< if true, debug output is printed

    TYPE(file_type), DIMENSION(:), ALLOCATABLE ::  files  !< file list

    SAVE

    PRIVATE

    !> Initialize the data-output module
    INTERFACE dom_init
       MODULE PROCEDURE dom_init
    END INTERFACE dom_init

    !> Add files to database
    INTERFACE dom_def_file
       MODULE PROCEDURE dom_def_file
    END INTERFACE dom_def_file

    !> Add dimensions to database
    INTERFACE dom_def_dim
       MODULE PROCEDURE dom_def_dim
    END INTERFACE dom_def_dim

    !> Add variables to database
    INTERFACE dom_def_var
       MODULE PROCEDURE dom_def_var
    END INTERFACE dom_def_var

    !> Add attributes to database
    INTERFACE dom_def_att
       MODULE PROCEDURE dom_def_att_char
       MODULE PROCEDURE dom_def_att_int8
       MODULE PROCEDURE dom_def_att_int16
       MODULE PROCEDURE dom_def_att_int32
       MODULE PROCEDURE dom_def_att_real32
       MODULE PROCEDURE dom_def_att_real64
    END INTERFACE dom_def_att

    !> Prepare for output: evaluate database and create files
    INTERFACE dom_def_end
       MODULE PROCEDURE dom_def_end
    END INTERFACE dom_def_end

    !> Write variables to file
    INTERFACE dom_write_var
       MODULE PROCEDURE dom_write_var
    END INTERFACE dom_write_var

    !> Last actions required for output befor termination
    INTERFACE dom_finalize_output
       MODULE PROCEDURE dom_finalize_output
    END INTERFACE dom_finalize_output

    !> Return error message
    INTERFACE dom_get_error_message
       MODULE PROCEDURE dom_get_error_message
    END INTERFACE dom_get_error_message

    !> Write database to debug output
    INTERFACE dom_database_debug_output
       MODULE PROCEDURE dom_database_debug_output
    END INTERFACE dom_database_debug_output

    PUBLIC &
       dom_init, &
       dom_def_file, &
       dom_def_dim, &
       dom_def_var, &
       dom_def_att, &
       dom_def_end, &
       dom_write_var, &
       dom_finalize_output, &
       dom_get_error_message, &
       dom_database_debug_output

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize data-output module.
!> Provide some general information of the main program.
!> The optional argument 'file_suffix_of_output_group' defines a file suffix which is added to all
!> output files. If multiple output groups (groups of MPI ranks, defined by
!> 'mpi_comm_of_output_group') exist, a unique file suffix must be given for each group. This
!> prevents that multiple groups try to open and write to the same output file.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dom_init( file_suffix_of_output_group, mpi_comm_of_output_group, master_output_rank, &
                      program_debug_output_unit, debug_output )

    CHARACTER(LEN=*), INTENT(IN), OPTIONAL ::  file_suffix_of_output_group  !< file-name suffix added to each file;
                                                                            !> must be unique for each output group

    INTEGER, INTENT(IN), OPTIONAL ::  master_output_rank         !< MPI rank executing tasks which must
                                                                 !> be executed by a single PE only
    INTEGER, INTENT(IN)           ::  mpi_comm_of_output_group   !< MPI communicator specifying the MPI group
                                                                 !> which participate in the output
    INTEGER, INTENT(IN)           ::  program_debug_output_unit  !< file unit number for debug output

    LOGICAL, INTENT(IN)           ::  debug_output               !< if true, debug output is printed


    IF ( PRESENT( file_suffix_of_output_group ) )  output_file_suffix = file_suffix_of_output_group
    IF ( PRESENT( master_output_rank ) )  master_rank = master_output_rank

    output_group_comm = mpi_comm_of_output_group

    debug_output_unit = program_debug_output_unit
    print_debug_output = debug_output

    CALL binary_init_module( output_file_suffix, output_group_comm, master_rank, &
                             debug_output_unit, debug_output, no_id )

    CALL netcdf4_init_module( output_file_suffix, output_group_comm, master_rank, &
                             debug_output_unit, debug_output, no_id )

 END SUBROUTINE dom_init

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Define output file.
!> Example call:
!>   status = dom_def_file( 'my_output_file_name', 'binary' )
!--------------------------------------------------------------------------------------------------!
 FUNCTION dom_def_file( file_name, file_format ) RESULT( return_value )

    CHARACTER(LEN=*), INTENT(IN) ::  file_name    !< name of file to be created
    CHARACTER(LEN=*), INTENT(IN) ::  file_format  !< format of file to be created

    CHARACTER(LEN=*), PARAMETER ::  routine_name = 'dom_def_file'  !< name of this routine

    INTEGER ::  f             !< loop index
    INTEGER ::  return_value  !< return value

    TYPE(file_type), DIMENSION(:), ALLOCATABLE ::  files_tmp  !< temporary file list


    return_value = 0

    CALL internal_message( 'debug', routine_name // ': define file "' // TRIM( file_name ) // '"' )
!
!-- Allocate file list or extend it by 1
    IF ( .NOT. ALLOCATED( files ) ) THEN

       nfiles = 1
       ALLOCATE( files(nfiles) )

    ELSE

       nfiles = SIZE( files )
!
!--    Check if file already exists
       DO  f = 1, nfiles
          IF ( files(f)%name == TRIM( file_name ) )  THEN
             return_value = 1
             CALL internal_message( 'error', routine_name // &
                     ': file "' // TRIM( file_name ) // '" already exists' )
             EXIT
          ENDIF
       ENDDO
!
!--    Extend file list
       IF ( return_value == 0 )  THEN
          ALLOCATE( files_tmp(nfiles) )
          files_tmp = files
          DEALLOCATE( files )
          nfiles = nfiles + 1
          ALLOCATE( files(nfiles) )
          files(:nfiles-1) = files_tmp
          DEALLOCATE( files_tmp )
       ENDIF

    ENDIF
!
!-- Add new file to database
    IF ( return_value == 0 )  THEN
       files(nfiles)%name = TRIM( file_name )
       files(nfiles)%format = TRIM( file_format )
    ENDIF

 END FUNCTION dom_def_file

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Define dimension.
!> Dimensions can either be limited (a lower and upper bound is given) or unlimited (only a lower
!> bound is given). Also, instead of providing all values of the dimension, a single value can be
!> given which is then used to fill the entire dimension.
!> An optional mask can be given to mask limited dimensions.
!> Example call:
!>   - fixed dimension with 100 entries (values known):
!>       status = dom_def_dim( file_name='my_output_file_name', dimension_name='my_dimension', &
!>                             output_type='real32', bounds=(/1,100/), &
!>                             values_real32=my_dim(1:100), mask=my_dim_mask(1:100) )
!>   - fixed dimension with 50 entries (values not yet known):
!>       status = dom_def_dim( file_name='my_output_file_name', dimension_name='my_dimension', &
!>                             output_type='int32', bounds=(/0,49/), &
!>                             values_int32=(/fill_value/) )
!>   - masked dimension with 75 entries:
!>       status = dom_def_dim( file_name='my_output_file_name', dimension_name='my_dimension', &
!>                             output_type='real64', bounds=(/101,175/), &
!>                             values_real64=my_dim(1:75), mask=my_dim_mask(1:75) )
!>   - unlimited dimension:
!>       status = dom_def_dim( file_name='my_output_file_name', dimension_name='my_dimension', &
!>                             output_type='real32', bounds=(/1/), &
!>                             values_real32=(/fill_value/) )
!>
!> @todo Convert given values into selected output_type.
!--------------------------------------------------------------------------------------------------!
 FUNCTION dom_def_dim( file_name, dimension_name, output_type, bounds,        &
                       values_int8, values_int16, values_int32, values_intwp, &
                       values_real32, values_real64, values_realwp,           &
                       mask ) RESULT( return_value )

    CHARACTER(LEN=*), INTENT(IN) ::  file_name       !< name of file
    CHARACTER(LEN=*), INTENT(IN) ::  dimension_name  !< name of dimension
    CHARACTER(LEN=*), INTENT(IN) ::  output_type     !< data type of dimension variable in output file

    CHARACTER(LEN=*), PARAMETER ::  routine_name = 'dom_def_dim'  !< name of this routine

    INTEGER ::  d             !< loop index
    INTEGER ::  f             !< loop index
    INTEGER ::  i             !< loop index
    INTEGER ::  j             !< loop index
    INTEGER ::  ndims         !< number of dimensions in file
    INTEGER ::  return_value  !< return value

    INTEGER,         DIMENSION(:), INTENT(IN)           ::  bounds         !< lower and upper bound of dimension variable
    INTEGER(KIND=1), DIMENSION(:), INTENT(IN), OPTIONAL ::  values_int8    !< values of dimension
    INTEGER(KIND=2), DIMENSION(:), INTENT(IN), OPTIONAL ::  values_int16   !< values of dimension
    INTEGER(KIND=4), DIMENSION(:), INTENT(IN), OPTIONAL ::  values_int32   !< values of dimension
    INTEGER(iwp),    DIMENSION(:), INTENT(IN), OPTIONAL ::  values_intwp   !< values of dimension

    LOGICAL,         DIMENSION(:), INTENT(IN), OPTIONAL ::  mask           !< mask of dimesion

    REAL(KIND=4),    DIMENSION(:), INTENT(IN), OPTIONAL ::  values_real32  !< values of dimension
    REAL(KIND=8),    DIMENSION(:), INTENT(IN), OPTIONAL ::  values_real64  !< values of dimension
    REAL(wp),        DIMENSION(:), INTENT(IN), OPTIONAL ::  values_realwp  !< values of dimension

    TYPE(dimension_type)                            ::  dimension       !< new dimension
    TYPE(dimension_type), DIMENSION(:), ALLOCATABLE ::  dimensions_tmp  !< temporary dimension list


    return_value = 0

    CALL internal_message( 'debug', routine_name //                    &
                           ': define dimension ' //                    &
                           '(dimension "' // TRIM( dimension_name ) // &
                           '", file "' // TRIM( file_name ) // '")' )

    dimension%name      = TRIM( dimension_name )
    dimension%data_type = TRIM( output_type )
!
!-- Check dimension bounds and allocate dimension according to bounds
    IF ( SIZE( bounds ) == 1 )  THEN
!
!--    Dimension has only lower bound, which means it changes its size
!--    during simulation.
!--    Set length to -1 as indicator.
       dimension%bounds(:) = bounds(1)
       dimension%length    = -1

       IF ( PRESENT( mask ) )  THEN
          return_value = 1
          CALL internal_message( 'error', routine_name //                      &
                                 ': unlimited dimensions cannot be masked ' // &
                                 '(dimension "' // TRIM( dimension_name ) //   &
                                 '", file "' // TRIM( file_name ) // '")!' )
       ENDIF

    ELSEIF ( SIZE( bounds ) == 2 )  THEN

       dimension%bounds = bounds
       dimension%length = bounds(2) - bounds(1) + 1
!
!--    Save dimension values
       IF ( PRESENT( values_int8 ) )  THEN
          ALLOCATE( dimension%values_int8(dimension%bounds(1):dimension%bounds(2)) )
          IF ( SIZE( values_int8 ) == dimension%length )  THEN
             dimension%values_int8 = values_int8
          ELSEIF ( SIZE( values_int8 ) == 1 )  THEN
             dimension%values_int8(:) = values_int8(1)
          ELSE
             return_value = 2
          ENDIF
       ELSEIF( PRESENT( values_int16 ) )  THEN
          ALLOCATE( dimension%values_int16(dimension%bounds(1):dimension%bounds(2)) )
          IF ( SIZE( values_int16 ) == dimension%length )  THEN
             dimension%values_int16 = values_int16
          ELSEIF ( SIZE( values_int16 ) == 1 )  THEN
             dimension%values_int16(:) = values_int16(1)
          ELSE
             return_value = 2
          ENDIF
       ELSEIF( PRESENT( values_int32 ) )  THEN
          ALLOCATE( dimension%values_int32(dimension%bounds(1):dimension%bounds(2)) )
          IF ( SIZE( values_int32 ) == dimension%length )  THEN
             dimension%values_int32 = values_int32
          ELSEIF ( SIZE( values_int32 ) == 1 )  THEN
             dimension%values_int32(:) = values_int32(1)
          ELSE
             return_value = 2
          ENDIF
       ELSEIF( PRESENT( values_intwp ) )  THEN
          ALLOCATE( dimension%values_intwp(dimension%bounds(1):dimension%bounds(2)) )
          IF ( SIZE( values_intwp ) == dimension%length )  THEN
             dimension%values_intwp = values_intwp
          ELSEIF ( SIZE( values_intwp ) == 1 )  THEN
             dimension%values_intwp(:) = values_intwp(1)
          ELSE
             return_value = 2
          ENDIF
       ELSEIF( PRESENT( values_real32 ) )  THEN
          ALLOCATE( dimension%values_real32(dimension%bounds(1):dimension%bounds(2)) )
          IF ( SIZE( values_real32 ) == dimension%length )  THEN
             dimension%values_real32 = values_real32
          ELSEIF ( SIZE( values_real32 ) == 1 )  THEN
             dimension%values_real32(:) = values_real32(1)
          ELSE
             return_value = 2
          ENDIF
       ELSEIF( PRESENT( values_real64 ) )  THEN
          ALLOCATE( dimension%values_real64(dimension%bounds(1):dimension%bounds(2)) )
          IF ( SIZE( values_real64 ) == dimension%length )  THEN
             dimension%values_real64 = values_real64
          ELSEIF ( SIZE( values_real64 ) == 1 )  THEN
             dimension%values_real64(:) = values_real64(1)
          ELSE
             return_value = 2
          ENDIF
       ELSEIF( PRESENT( values_realwp ) )  THEN
          ALLOCATE( dimension%values_realwp(dimension%bounds(1):dimension%bounds(2)) )
          IF ( SIZE( values_realwp ) == dimension%length )  THEN
             dimension%values_realwp = values_realwp
          ELSEIF ( SIZE( values_realwp ) == 1 )  THEN
             dimension%values_realwp(:) = values_realwp(1)
          ELSE
             return_value = 2
          ENDIF
       ELSE
          return_value = 1
          CALL internal_message( 'error', routine_name //                    &
                                 ': no values given ' //                     &
                                 '(dimension "' // TRIM( dimension_name ) // &
                                 '", file "' // TRIM( file_name ) // '")!' )
       ENDIF

       IF ( return_value == 2 )  THEN
          return_value = 1
          CALL internal_message( 'error', routine_name //                               &
                                 ': number of values and given bounds do not match ' // &
                                 '(dimension "' // TRIM( dimension_name ) //            &
                                 '", file "' // TRIM( file_name ) // '")!' )
       ENDIF
!
!--    Initialize mask
       IF ( PRESENT( mask )  .AND.  return_value == 0 )  THEN

          IF ( dimension%length == SIZE( mask ) )  THEN

             IF ( ALL( mask ) )  THEN

                CALL internal_message( 'debug', routine_name //                              &
                                       ': mask contains only TRUE values. Ignoring mask ' // &
                                       '(dimension "' // TRIM( dimension_name ) //           &
                                       '", file "' // TRIM( file_name ) // '")!' )

             ELSE

                dimension%is_masked = .TRUE.
                dimension%length_mask = COUNT( mask )

                ALLOCATE( dimension%mask(dimension%bounds(1):dimension%bounds(2)) )
                ALLOCATE( dimension%masked_indices(0:dimension%length_mask-1) )

                dimension%mask = mask
!
!--             Save masked positions and masked values
                IF ( ALLOCATED( dimension%values_int8 ) )  THEN

                   ALLOCATE( dimension%masked_values_int8(0:dimension%length_mask-1) )
                   j = 0
                   DO  i = dimension%bounds(1), dimension%bounds(2)
                      IF ( dimension%mask(i) )  THEN
                         dimension%masked_values_int8(j) = dimension%values_int8(i)
                         dimension%masked_indices(j) = i
                         j = j + 1
                      ENDIF
                   ENDDO

                ELSEIF ( ALLOCATED( dimension%values_int16 ) )  THEN

                   ALLOCATE( dimension%masked_values_int16(0:dimension%length_mask-1) )
                   j = 0
                   DO  i = dimension%bounds(1), dimension%bounds(2)
                      IF ( dimension%mask(i) )  THEN
                         dimension%masked_values_int16(j) = dimension%values_int16(i)
                         dimension%masked_indices(j) = i
                         j = j + 1
                      ENDIF
                   ENDDO

                ELSEIF ( ALLOCATED( dimension%values_int32 ) )  THEN

                   ALLOCATE( dimension%masked_values_int32(0:dimension%length_mask-1) )
                   j = 0
                   DO  i =dimension%bounds(1), dimension%bounds(2)
                      IF ( dimension%mask(i) )  THEN
                         dimension%masked_values_int32(j) = dimension%values_int32(i)
                         dimension%masked_indices(j) = i
                         j = j + 1
                      ENDIF
                   ENDDO

                ELSEIF ( ALLOCATED( dimension%values_intwp ) )  THEN

                   ALLOCATE( dimension%masked_values_intwp(0:dimension%length_mask-1) )
                   j = 0
                   DO  i = dimension%bounds(1), dimension%bounds(2)
                      IF ( dimension%mask(i) )  THEN
                         dimension%masked_values_intwp(j) = dimension%values_intwp(i)
                         dimension%masked_indices(j) = i
                         j = j + 1
                      ENDIF
                   ENDDO

                ELSEIF ( ALLOCATED( dimension%values_real32 ) )  THEN

                   ALLOCATE( dimension%masked_values_real32(0:dimension%length_mask-1) )
                   j = 0
                   DO  i = dimension%bounds(1), dimension%bounds(2)
                      IF ( dimension%mask(i) )  THEN
                         dimension%masked_values_real32(j) = dimension%values_real32(i)
                         dimension%masked_indices(j) = i
                         j = j + 1
                      ENDIF
                   ENDDO

                ELSEIF ( ALLOCATED(dimension%values_real64) )  THEN

                   ALLOCATE( dimension%masked_values_real64(0:dimension%length_mask-1) )
                   j = 0
                   DO  i = dimension%bounds(1), dimension%bounds(2)
                      IF ( dimension%mask(i) )  THEN
                         dimension%masked_values_real64(j) = dimension%values_real64(i)
                         dimension%masked_indices(j) = i
                         j = j + 1
                      ENDIF
                   ENDDO

                ELSEIF ( ALLOCATED(dimension%values_realwp) )  THEN

                   ALLOCATE( dimension%masked_values_realwp(0:dimension%length_mask-1) )
                   j = 0
                   DO  i = dimension%bounds(1), dimension%bounds(2)
                      IF ( dimension%mask(i) )  THEN
                         dimension%masked_values_realwp(j) = dimension%values_realwp(i)
                         dimension%masked_indices(j) = i
                         j = j + 1
                      ENDIF
                   ENDDO

                ENDIF

             ENDIF  ! if not all mask = true

          ELSE
             return_value = 1
             CALL internal_message( 'error', routine_name //                           &
                                    ': size of mask and given bounds do not match ' // &
                                    '(dimension "' // TRIM( dimension_name ) //        &
                                    '", file "' // TRIM( file_name ) // '")!' )
          ENDIF

       ENDIF

    ELSE

       return_value = 1
       CALL internal_message( 'error', routine_name //                                       &
                              ': at least one but no more than two bounds must be given ' // &
                              '(dimension "' // TRIM( dimension_name ) //                    &
                              '", file "' // TRIM( file_name ) // '")!' )

    ENDIF
!
!-- Add dimension to database
    IF ( return_value == 0 )  THEN

       DO  f = 1, nfiles

          IF ( TRIM( file_name ) == files(f)%name )  THEN

             IF ( files(f)%is_init )  THEN

                return_value = 1
                CALL internal_message( 'error', routine_name //                      &
                                       ': file already initialized. ' //             &
                                       'No further dimension definition allowed ' // &
                                       '(dimension "' // TRIM( dimension_name ) //   &
                                       '", file "' // TRIM( file_name ) // '")!' )
                EXIT

             ELSEIF ( .NOT. ALLOCATED( files(f)%dimensions ) )  THEN

                ndims = 1
                ALLOCATE( files(f)%dimensions(ndims) )

             ELSE
!
!--             Check if any variable of the same name as the new dimension is already defined
                IF ( ALLOCATED( files(f)%variables ) )  THEN
                   DO  i = 1, SIZE( files(f)%variables )
                      IF ( files(f)%variables(i)%name == dimension%name )  THEN
                         return_value = 1
                         CALL internal_message( 'error', routine_name //                    &
                                 ': file already has a variable of this name defined. ' //  &
                                 'Defining a dimension of the same name is not allowed ' // &
                                 '(dimension "' // TRIM( dimension_name ) //                &
                                 '", file "' // TRIM( file_name ) // '")!' )
                         EXIT
                      ENDIF
                   ENDDO
                ENDIF

                IF ( return_value == 0 )  THEN
!
!--                Check if dimension already exists in file
                   ndims = SIZE( files(f)%dimensions )

                   DO  d = 1, ndims
                      IF ( files(f)%dimensions(d)%name == dimension%name )  THEN
                         return_value = 1
                         CALL internal_message( 'error', routine_name //     &
                                 ': dimension already exists in file ' //    &
                                 '(dimension "' // TRIM( dimension_name ) // &
                                 '", file "' // TRIM( file_name ) // '")!' )
                         EXIT
                      ENDIF
                   ENDDO
!
!--                Extend dimension list
                   IF ( return_value == 0 )  THEN
                      ALLOCATE( dimensions_tmp(ndims) )
                      dimensions_tmp = files(f)%dimensions
                      DEALLOCATE( files(f)%dimensions )
                      ndims = ndims + 1
                      ALLOCATE( files(f)%dimensions(ndims) )
                      files(f)%dimensions(:ndims-1) = dimensions_tmp
                      DEALLOCATE( dimensions_tmp )
                   ENDIF
                ENDIF

             ENDIF
!
!--          Add new dimension to database
             IF ( return_value == 0 )  files(f)%dimensions(ndims) = dimension

             EXIT

          ENDIF
       ENDDO

       IF ( f > nfiles )  THEN
          return_value = 1
          CALL internal_message( 'error', routine_name //                                     &
                                 ': file not found (dimension "' // TRIM( dimension_name ) // &
                                 '", file "' // TRIM( file_name ) // '")!' )
       ENDIF

    ENDIF

 END FUNCTION dom_def_dim

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Add variable to database.
!> If a variable is identical for each MPI rank, the optional argument 'is_global' should be set to
!> TRUE. This flags the variable to be a global variable and is later only written once by the
!> master output rank.
!> Example call:
!>   dom_def_var( file_name =  'my_output_file_name', &
!>                variable_name = 'u', &
!>                dimension_names = (/'x   ', 'y   ', 'z   ', 'time'/), &
!>                output_type = 'real32' )
!> @note The order of dimensions must match in reversed order to the dimensions of the
!>       corresponding variable array. The last given dimension can also be non-existent within the
!>       variable array if at any given call of 'dom_write_var' for this variable, the last
!>       dimension has only a single index.
!>       Hence, the array 'u' must be allocated with dimension 'x' as its last dimension, preceded
!>       by 'y', then 'z', and 'time' being the first dimension. If at any given write statement,
!>       only a single index of dimension 'time' is to be written, the dimension can be non-present
!>       in the variable array leaving dimension 'z' as the first dimension.
!>       So, the variable array needs to be allocated like either:
!>          ALLOCATE( u(<time>,<z>,<y>,<x>) )
!>       or
!>          ALLOCATE( u(<z>,<y>,<x>) )
!--------------------------------------------------------------------------------------------------!
 FUNCTION dom_def_var( file_name, variable_name, dimension_names, output_type, is_global ) &
             RESULT( return_value )

    CHARACTER(LEN=*), INTENT(IN) ::  file_name      !< name of file
    CHARACTER(LEN=*), INTENT(IN) ::  variable_name  !< name of variable
    CHARACTER(LEN=*), INTENT(IN) ::  output_type    !< data type of variable

    CHARACTER(LEN=*), PARAMETER ::  routine_name = 'dom_def_var'  !< name of this routine

    CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) ::  dimension_names  !< list of dimension names

    INTEGER ::  d             !< loop index
    INTEGER ::  f             !< loop index
    INTEGER ::  i             !< loop index
    INTEGER ::  nvars         !< number of variables in file
    INTEGER ::  return_value  !< return value

    LOGICAL                       ::  found      !< true if requested dimension is defined in file
    LOGICAL, INTENT(IN), OPTIONAL ::  is_global  !< true if variable is global (same on all PE)

    TYPE(variable_type)                            ::  variable       !< new variable
    TYPE(variable_type), DIMENSION(:), ALLOCATABLE ::  variables_tmp  !< temporary variable list


    return_value = 0
    found = .FALSE.

    CALL internal_message( 'debug', routine_name //                                     &
                           ': define variable (variable "' // TRIM( variable_name ) //  &
                           '", file "' // TRIM( file_name ) // '")' )

    variable%name = TRIM( variable_name )

    ALLOCATE( variable%dimension_names(SIZE( dimension_names )) )
    ALLOCATE( variable%dimension_ids(SIZE( dimension_names )) )

    variable%dimension_names = dimension_names
    variable%dimension_ids = -1
    variable%data_type = TRIM( output_type )

    IF ( PRESENT( is_global ) )  THEN
       variable%is_global = is_global
    ELSE
       variable%is_global = .FALSE.
    ENDIF
!
!-- Add variable to database
    DO  f = 1, nfiles

       IF ( TRIM( file_name ) == files(f)%name )  THEN

          IF ( files(f)%is_init )  THEN

             return_value = 1
             CALL internal_message( 'error', routine_name //                                  &
                     ': file already initialized. No further variable definition allowed ' // &
                     '(variable "' // TRIM( variable_name ) //                                &
                     '", file "' // TRIM( file_name ) // '")!' )
             EXIT

          ELSEIF ( ALLOCATED( files(f)%dimensions ) )  THEN
!
!--          Check if any dimension of the same name as the new variable is already defined
             DO  d = 1, SIZE( files(f)%dimensions )
                IF ( files(f)%dimensions(d)%name == variable%name )  THEN
                   return_value = 1
                   CALL internal_message( 'error', routine_name //                    &
                           ': file already has a dimension of this name defined. ' // &
                           'Defining a variable of the same name is not allowed ' //  &
                           '(variable "' // TRIM( variable_name ) //                  &
                           '", file "' // TRIM( file_name ) // '")!' )
                   EXIT
                ENDIF
             ENDDO
!
!--          Check if dimensions assigned to variable are defined within file
             IF ( return_value == 0 )  THEN
                DO  i = 1, SIZE( variable%dimension_names )
                   found = .FALSE.
                   DO  d = 1, SIZE( files(f)%dimensions )
                      IF ( files(f)%dimensions(d)%name == variable%dimension_names(i) )  THEN
                         found = .TRUE.
                         EXIT
                      ENDIF
                   ENDDO
                   IF ( .NOT. found )  THEN
                      return_value = 1
                      CALL internal_message( 'error', routine_name //                            &
                              ': required dimension "'//  TRIM( variable%dimension_names(i) ) // &
                              '" for variable is not defined ' //                                &
                              '(variable "' // TRIM( variable_name ) //                          &
                              '", file "' // TRIM( file_name ) // '")!' )
                      EXIT
                   ENDIF
                ENDDO
             ENDIF

          ELSE

             return_value = 1
             CALL internal_message( 'error', routine_name //                      &
                     ': no dimensions defined in file. Cannot define variable '// &
                     '(variable "' // TRIM( variable_name ) //                    &
                     '", file "' // TRIM( file_name ) // '")!' )

          ENDIF

          IF ( return_value == 0 )  THEN
!
!--          Check if variable already exists
             IF ( .NOT. ALLOCATED( files(f)%variables ) )  THEN

                nvars = 1
                ALLOCATE( files(f)%variables(nvars) )

             ELSE

                nvars = SIZE( files(f)%variables )
                DO  i = 1, nvars
                   IF ( files(f)%variables(i)%name == variable%name )  THEN
                      return_value = 1
                      CALL internal_message( 'error', routine_name //   &
                              ': variable already exists '//            &
                              '(variable "' // TRIM( variable_name ) // &
                              '", file "' // TRIM( file_name ) // '")!' )
                      EXIT
                   ENDIF
                ENDDO

                IF ( return_value == 0 )  THEN
!
!--                Extend variable list
                   ALLOCATE( variables_tmp(nvars) )
                   variables_tmp = files(f)%variables
                   DEALLOCATE( files(f)%variables )
                   nvars = nvars + 1
                   ALLOCATE( files(f)%variables(nvars) )
                   files(f)%variables(:nvars-1) = variables_tmp
                   DEALLOCATE( variables_tmp )
                ENDIF

             ENDIF
!
!--          Add new variable to database
             IF ( return_value == 0 )  files(f)%variables(nvars) = variable

          ENDIF

          EXIT

       ENDIF

    ENDDO

    IF ( f > nfiles )  THEN
       return_value = 1
       CALL internal_message( 'error', routine_name //                                   &
                              ': file not found (variable "' // TRIM( variable_name ) // &
                              '", file "' // TRIM( file_name ) // '")!' )
    ENDIF

 END FUNCTION dom_def_var

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Create attribute with value of type character.
!> If the optional argument 'variable_name' is given, the attribute is added to the respective
!> variable or dimension of that name. Otherwise, the attribute is added as a global attribute to
!> the file itself.
!> If an attribute of similar name already exists, it is updated (overwritten) with the new value.
!> If the optional argument 'append' is set TRUE, the value of an already existing attribute of
!> similar name is appended by the new value instead of overwritten.
!> Example call:
!>   - define a global file attribute:
!>      dom_def_att( file_name='my_output_file_name', &
!>                   attribute_name='my_attribute', &
!>                   value='This is the attribute value' )
!>   - define a variable attribute:
!>      dom_def_att( file_name='my_output_file_name', &
!>                   variable_name='my_variable', &
!>                   attribute_name='my_attribute', &
!>                   value='This is the attribute value' )
!>   - append an attribute:
!>      dom_def_att( file_name='my_output_file_name', &
!>                   attribute_name='my_attribute', &
!>                   value=' and this part was appended', append=.TRUE. )
!--------------------------------------------------------------------------------------------------!
 FUNCTION dom_def_att_char( file_name, variable_name, attribute_name, value, append ) &
             RESULT( return_value )

    CHARACTER(LEN=*),      INTENT(IN)           ::  file_name               !< name of file
    CHARACTER(LEN=*),      INTENT(IN)           ::  attribute_name          !< name of attribute
    CHARACTER(LEN=*),      INTENT(IN)           ::  value                   !< attribute value
    CHARACTER(LEN=*),      INTENT(IN), OPTIONAL ::  variable_name           !< name of variable
    CHARACTER(LEN=charlen)                      ::  variable_name_internal  !< internal copy of variable_name

    ! CHARACTER(LEN=*), PARAMETER ::  routine_name = 'dom_def_att_char'  !< name of routine

    INTEGER ::  return_value  !< return value

    LOGICAL                       ::  append_internal  !< same as 'append'
    LOGICAL, INTENT(IN), OPTIONAL ::  append           !< if true, append value to existing value

    TYPE(attribute_type) ::  attribute  !< new attribute


    return_value = 0

    IF ( PRESENT( append ) )  THEN
       append_internal = append
    ELSE
       append_internal = .FALSE.
    ENDIF

    attribute%name       = TRIM( attribute_name )
    attribute%data_type  = 'char'
    attribute%value_char = TRIM( value )

    IF ( PRESENT( variable_name ) )  THEN
       variable_name_internal = TRIM( variable_name )
    ELSE
       variable_name_internal = ''
    ENDIF

    return_value = save_attribute_in_database( file_name=TRIM( file_name ), &
                      variable_name=TRIM( variable_name_internal ),         &
                      attribute=attribute, append=append_internal )

 END FUNCTION dom_def_att_char

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Create attribute with value of type int8.
!> If the optional argument 'variable_name' is given, the attribute is added to the respective
!> variable or dimension of that name. Otherwise, the attribute is added as a global attribute to
!> the file itself.
!> Numerical attributes cannot be appended, only updated (append=.TRUE. will cause an error).
!> Example call:
!>   - define a global file attribute:
!>      dom_def_att( file_name='my_output_file_name', &
!>                   attribute_name='my_attribute', &
!>                   value=0_1 )
!>   - define a variable attribute:
!>      dom_def_att( file_name='my_output_file_name', &
!>                   variable_name='my_variable', &
!>                   attribute_name='my_attribute', &
!>                   value=1_1 )
!--------------------------------------------------------------------------------------------------!
 FUNCTION dom_def_att_int8( file_name, variable_name, attribute_name, value, append ) &
             RESULT( return_value )

    CHARACTER(LEN=*),      INTENT(IN)           ::  file_name               !< name of file
    CHARACTER(LEN=*),      INTENT(IN)           ::  attribute_name          !< name of attribute
    CHARACTER(LEN=*),      INTENT(IN), OPTIONAL ::  variable_name           !< name of variable
    CHARACTER(LEN=charlen)                      ::  variable_name_internal  !< internal copy of variable_name

    CHARACTER(LEN=*), PARAMETER ::  routine_name = 'dom_def_att_int8'  !< name of routine

    INTEGER(KIND=1), INTENT(IN) ::  value  !< attribute value

    INTEGER ::  return_value  !< return value

    LOGICAL                       ::  append_internal  !< same as 'append'
    LOGICAL, INTENT(IN), OPTIONAL ::  append           !< if true, append value to existing value

    TYPE(attribute_type) ::  attribute  !< new attribute


    return_value = 0

    IF ( PRESENT( variable_name ) )  THEN
       variable_name_internal = TRIM( variable_name )
    ELSE
       variable_name_internal = ''
    ENDIF

    IF ( PRESENT( append ) )  THEN
       IF ( append )  THEN
          return_value = 1
          CALL internal_message( 'error', routine_name //                             &
                                 ': numeric attribute cannot be appended ' //         &
                                 '(attribute "' // TRIM( attribute_name ) //          &
                                 '", variable "' // TRIM( variable_name_internal ) // &
                                 '", file "' // TRIM( file_name ) // '")!' )
       ENDIF
    ENDIF

    IF ( return_value == 0 )  THEN
       append_internal = .FALSE.

       attribute%name       = TRIM( attribute_name )
       attribute%data_type  = 'int8'
       attribute%value_int8 = value

       return_value = save_attribute_in_database( file_name=TRIM( file_name ), &
                         variable_name=TRIM( variable_name_internal ),         &
                         attribute=attribute, append=append_internal )
    ENDIF

 END FUNCTION dom_def_att_int8

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Create attribute with value of type int16.
!> If the optional argument 'variable_name' is given, the attribute is added to the respective
!> variable or dimension of that name. Otherwise, the attribute is added as a global attribute to
!> the file itself.
!> Numerical attributes cannot be appended, only updated (append=.TRUE. will cause an error).
!> Example call:
!>   - define a global file attribute:
!>      dom_def_att( file_name='my_output_file_name', &
!>                   attribute_name='my_attribute', &
!>                   value=0_2 )
!>   - define a variable attribute:
!>      dom_def_att( file_name='my_output_file_name', &
!>                   variable_name='my_variable', &
!>                   attribute_name='my_attribute', &
!>                   value=1_2 )
!--------------------------------------------------------------------------------------------------!
 FUNCTION dom_def_att_int16( file_name, variable_name, attribute_name, value, append ) &
             RESULT( return_value )

    CHARACTER(LEN=*),      INTENT(IN)           ::  file_name               !< name of file
    CHARACTER(LEN=*),      INTENT(IN)           ::  attribute_name          !< name of attribute
    CHARACTER(LEN=*),      INTENT(IN), OPTIONAL ::  variable_name           !< name of variable
    CHARACTER(LEN=charlen)                      ::  variable_name_internal  !< internal copy of variable_name

    CHARACTER(LEN=*), PARAMETER ::  routine_name = 'dom_def_att_int16'  !< name of routine

    INTEGER(KIND=2), INTENT(IN) ::  value  !< attribute value

    INTEGER ::  return_value  !< return value

    LOGICAL                       ::  append_internal  !< same as 'append'
    LOGICAL, INTENT(IN), OPTIONAL ::  append           !< if true, append value to existing value

    TYPE(attribute_type) ::  attribute  !< new attribute


    return_value = 0

    IF ( PRESENT( variable_name ) )  THEN
       variable_name_internal = TRIM( variable_name )
    ELSE
       variable_name_internal = ''
    ENDIF

    IF ( PRESENT( append ) )  THEN
       IF ( append )  THEN
          return_value = 1
          CALL internal_message( 'error', routine_name //                             &
                                 ': numeric attribute cannot be appended ' //         &
                                 '(attribute "' // TRIM( attribute_name ) //          &
                                 '", variable "' // TRIM( variable_name_internal ) // &
                                 '", file "' // TRIM( file_name ) // '")!' )
       ENDIF
    ENDIF

    IF ( return_value == 0 )  THEN
       append_internal = .FALSE.

       attribute%name        = TRIM( attribute_name )
       attribute%data_type   = 'int16'
       attribute%value_int16 = value

       return_value = save_attribute_in_database( file_name=TRIM( file_name ), &
                         variable_name=TRIM( variable_name_internal ),         &
                         attribute=attribute, append=append_internal )
    ENDIF

 END FUNCTION dom_def_att_int16

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Create attribute with value of type int32.
!> If the optional argument 'variable_name' is given, the attribute is added to the respective
!> variable or dimension of that name. Otherwise, the attribute is added as a global attribute to
!> the file itself.
!> Numerical attributes cannot be appended, only updated (append=.TRUE. will cause an error).
!> Example call:
!>   - define a global file attribute:
!>      dom_def_att( file_name='my_output_file_name', &
!>                   attribute_name='my_attribute', &
!>                   value=0_4 )
!>   - define a variable attribute:
!>      dom_def_att( file_name='my_output_file_name', &
!>                   variable_name='my_variable', &
!>                   attribute_name='my_attribute', &
!>                   value=1_4 )
!--------------------------------------------------------------------------------------------------!
 FUNCTION dom_def_att_int32( file_name, variable_name, attribute_name, value, append ) &
             RESULT( return_value )

    CHARACTER(LEN=*),      INTENT(IN)           ::  file_name               !< name of file
    CHARACTER(LEN=*),      INTENT(IN)           ::  attribute_name          !< name of attribute
    CHARACTER(LEN=*),      INTENT(IN), OPTIONAL ::  variable_name           !< name of variable
    CHARACTER(LEN=charlen)                      ::  variable_name_internal  !< internal copy of variable_name

    CHARACTER(LEN=*), PARAMETER ::  routine_name = 'dom_def_att_int32'  !< name of routine

    INTEGER(KIND=4), INTENT(IN) ::  value  !< attribute value

    INTEGER ::  return_value  !< return value

    LOGICAL                       ::  append_internal  !< same as 'append'
    LOGICAL, INTENT(IN), OPTIONAL ::  append           !< if true, append value to existing value

    TYPE(attribute_type) ::  attribute  !< new attribute


    return_value = 0

    IF ( PRESENT( variable_name ) )  THEN
       variable_name_internal = TRIM( variable_name )
    ELSE
       variable_name_internal = ''
    ENDIF

    IF ( PRESENT( append ) )  THEN
       IF ( append )  THEN
          return_value = 1
          CALL internal_message( 'error', routine_name //                             &
                                 ': numeric attribute cannot be appended ' //         &
                                 '(attribute "' // TRIM( attribute_name ) //          &
                                 '", variable "' // TRIM( variable_name_internal ) // &
                                 '", file "' // TRIM( file_name ) // '")!' )
       ENDIF
    ENDIF

    IF ( return_value == 0 )  THEN
       append_internal = .FALSE.

       attribute%name        = TRIM( attribute_name )
       attribute%data_type   = 'int32'
       attribute%value_int32 = value

       return_value = save_attribute_in_database( file_name=TRIM( file_name ), &
                         variable_name=TRIM( variable_name_internal ),         &
                         attribute=attribute, append=append_internal )
    ENDIF

 END FUNCTION dom_def_att_int32

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Create attribute with value of type real32.
!> If the optional argument 'variable_name' is given, the attribute is added to the respective
!> variable or dimension of that name. Otherwise, the attribute is added as a global attribute to
!> the file itself.
!> Numerical attributes cannot be appended, only updated (append=.TRUE. will cause an error).
!> Example call:
!>   - define a global file attribute:
!>      dom_def_att( file_name='my_output_file_name', &
!>                   attribute_name='my_attribute', &
!>                   value=1.0_4 )
!>   - define a variable attribute:
!>      dom_def_att( file_name='my_output_file_name', &
!>                   variable_name='my_variable', &
!>                   attribute_name='my_attribute', &
!>                   value=1.0_4 )
!--------------------------------------------------------------------------------------------------!
 FUNCTION dom_def_att_real32( file_name, variable_name, attribute_name, value, append ) &
             RESULT( return_value )

    CHARACTER(LEN=*),      INTENT(IN)           ::  file_name               !< name of file
    CHARACTER(LEN=*),      INTENT(IN)           ::  attribute_name          !< name of attribute
    CHARACTER(LEN=*),      INTENT(IN), OPTIONAL ::  variable_name           !< name of variable
    CHARACTER(LEN=charlen)                      ::  variable_name_internal  !< internal copy of variable_name

    CHARACTER(LEN=*), PARAMETER ::  routine_name = 'dom_def_att_real32'  !< name of routine

    INTEGER ::  return_value  !< return value

    LOGICAL                       ::  append_internal  !< same as 'append'
    LOGICAL, INTENT(IN), OPTIONAL ::  append           !< if true, append value to existing value

    REAL(KIND=4), INTENT(IN) ::  value  !< attribute value

    TYPE(attribute_type) ::  attribute  !< new attribute


    return_value = 0

    IF ( PRESENT( variable_name ) )  THEN
       variable_name_internal = TRIM( variable_name )
    ELSE
       variable_name_internal = ''
    ENDIF

    IF ( PRESENT( append ) )  THEN
       IF ( append )  THEN
          return_value = 1
          CALL internal_message( 'error', routine_name //                             &
                                 ': numeric attribute cannot be appended ' //         &
                                 '(attribute "' // TRIM( attribute_name ) //          &
                                 '", variable "' // TRIM( variable_name_internal ) // &
                                 '", file "' // TRIM( file_name ) // '")!' )
       ENDIF
    ENDIF

    IF ( return_value == 0 )  THEN
       append_internal = .FALSE.

       attribute%name         = TRIM( attribute_name )
       attribute%data_type    = 'real32'
       attribute%value_real32 = value

       return_value = save_attribute_in_database( file_name=TRIM( file_name ), &
                         variable_name=TRIM( variable_name_internal ),         &
                         attribute=attribute, append=append_internal )
    ENDIF

 END FUNCTION dom_def_att_real32

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Create attribute with value of type real64.
!> If the optional argument 'variable_name' is given, the attribute is added to the respective
!> variable or dimension of that name. Otherwise, the attribute is added as a global attribute to
!> the file itself.
!> Numerical attributes cannot be appended, only updated (append=.TRUE. will cause an error).
!> Example call:
!>   - define a global file attribute:
!>      dom_def_att( file_name='my_output_file_name', &
!>                   attribute_name='my_attribute', &
!>                   value=0.0_8 )
!>   - define a variable attribute:
!>      dom_def_att( file_name='my_output_file_name', &
!>                   variable_name='my_variable', &
!>                   attribute_name='my_attribute', &
!>                   value=1.0_8 )
!--------------------------------------------------------------------------------------------------!
 FUNCTION dom_def_att_real64( file_name, variable_name, attribute_name, value, append ) &
             RESULT( return_value )

    CHARACTER(LEN=*),      INTENT(IN)           ::  file_name               !< name of file
    CHARACTER(LEN=*),      INTENT(IN)           ::  attribute_name          !< name of attribute
    CHARACTER(LEN=*),      INTENT(IN), OPTIONAL ::  variable_name           !< name of variable
    CHARACTER(LEN=charlen)                      ::  variable_name_internal  !< internal copy of variable_name

    CHARACTER(LEN=*), PARAMETER ::  routine_name = 'dom_def_att_real64'  !< name of routine

    INTEGER ::  return_value  !< return value

    LOGICAL                       ::  append_internal  !< same as 'append'
    LOGICAL, INTENT(IN), OPTIONAL ::  append           !< if true, append value to existing value

    REAL(KIND=8), INTENT(IN) ::  value  !< attribute value

    TYPE(attribute_type) ::  attribute  !< new attribute


    return_value = 0

    IF ( PRESENT( variable_name ) )  THEN
       variable_name_internal = TRIM( variable_name )
    ELSE
       variable_name_internal = ''
    ENDIF

    IF ( PRESENT( append ) )  THEN
       IF ( append )  THEN
          return_value = 1
          CALL internal_message( 'error', routine_name //                             &
                                 ': numeric attribute cannot be appended ' //         &
                                 '(attribute "' // TRIM( attribute_name ) //          &
                                 '", variable "' // TRIM( variable_name_internal ) // &
                                 '", file "' // TRIM( file_name ) // '")!' )
       ENDIF
    ENDIF

    IF ( return_value == 0 )  THEN
       append_internal = .FALSE.

       attribute%name         = TRIM( attribute_name )
       attribute%data_type    = 'real64'
       attribute%value_real64 = value

       return_value = save_attribute_in_database( file_name=TRIM( file_name ), &
                         variable_name=TRIM( variable_name_internal ),         &
                         attribute=attribute, append=append_internal )
    ENDIF

 END FUNCTION dom_def_att_real64

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> End output definition.
!> The database is cleared from unused files and dimensions. Then, the output files are initialized
!> and prepared for writing output values to them. The saved values of the dimensions are written
!> to the files.
!--------------------------------------------------------------------------------------------------!
 FUNCTION dom_def_end() RESULT( return_value )

    CHARACTER(LEN=*), PARAMETER ::  routine_name = 'dom_def_end'  !< name of routine

    INTEGER ::  d             !< loop index
    INTEGER ::  f             !< loop index
    INTEGER ::  return_value  !< return value

    INTEGER(KIND=1), DIMENSION(:), ALLOCATABLE, TARGET ::  values_int8           !< target array for dimension values
    INTEGER(KIND=2), DIMENSION(:), ALLOCATABLE, TARGET ::  values_int16          !< target array for dimension values
    INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE, TARGET ::  values_int32          !< target array for dimension values
    INTEGER(iwp),    DIMENSION(:), ALLOCATABLE, TARGET ::  values_intwp          !< target array for dimension values

    INTEGER(KIND=1), DIMENSION(:), POINTER, CONTIGUOUS ::  values_int8_pointer   !< pointer to target array
    INTEGER(KIND=2), DIMENSION(:), POINTER, CONTIGUOUS ::  values_int16_pointer  !< pointer to target array
    INTEGER(KIND=4), DIMENSION(:), POINTER, CONTIGUOUS ::  values_int32_pointer  !< pointer to target array
    INTEGER(iwp),    DIMENSION(:), POINTER, CONTIGUOUS ::  values_intwp_pointer  !< pointer to target array

    REAL(KIND=4), DIMENSION(:), ALLOCATABLE, TARGET ::  values_real32            !< target array for dimension values
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE, TARGET ::  values_real64            !< target array for dimension values
    REAL(wp),     DIMENSION(:), ALLOCATABLE, TARGET ::  values_realwp            !< target array for dimension values

    REAL(KIND=4), DIMENSION(:), POINTER, CONTIGUOUS ::  values_real32_pointer    !< pointer to target array
    REAL(KIND=8), DIMENSION(:), POINTER, CONTIGUOUS ::  values_real64_pointer    !< pointer to target array
    REAL(wp),     DIMENSION(:), POINTER, CONTIGUOUS ::  values_realwp_pointer    !< pointer to target array


    return_value = 0
    CALL internal_message( 'debug', routine_name // ': start' )
!
!-- Clear database from empty files and unused dimensions
    IF ( nfiles > 0 )  return_value = cleanup_database()

    IF ( return_value == 0 )  THEN
       DO  f = 1, nfiles
!
!--       Skip initialization if file is already initialized
          IF ( files(f)%is_init )  CYCLE

          CALL internal_message( 'debug', routine_name // ': initialize file "' // &
                                 TRIM( files(f)%name ) // '"' )
!
!--       Open file
          CALL open_output_file( files(f)%format, files(f)%name, files(f)%id, &
                                 return_value=return_value )
!
!--       Initialize file header:
!--       define dimensions and variables and write attributes
          IF ( return_value == 0 )  &
             CALL init_file_header( files(f), return_value=return_value )
!
!--       End file definition
          IF ( return_value == 0 )  &
             CALL stop_file_header_definition( files(f)%format, files(f)%id, &
                                               files(f)%name, return_value )

          IF ( return_value == 0 )  THEN
!
!--          Flag file as initialized
             files(f)%is_init = .TRUE.
!
!--          Write dimension values into file
             DO  d = 1, SIZE( files(f)%dimensions )
                IF ( ALLOCATED( files(f)%dimensions(d)%values_int8 ) )  THEN
                   ALLOCATE( values_int8(files(f)%dimensions(d)%bounds(1): &
                                         files(f)%dimensions(d)%bounds(2)) )
                   values_int8 = files(f)%dimensions(d)%values_int8
                   values_int8_pointer => values_int8
                   return_value = dom_write_var( files(f)%name, files(f)%dimensions(d)%name, &
                                     bounds_start=(/ files(f)%dimensions(d)%bounds(1) /),    &
                                     bounds_end  =(/ files(f)%dimensions(d)%bounds(2) /),    &
                                     values_int8_1d=values_int8_pointer )
                   DEALLOCATE( values_int8 )
                ELSEIF ( ALLOCATED( files(f)%dimensions(d)%values_int16 ) )  THEN
                   ALLOCATE( values_int16(files(f)%dimensions(d)%bounds(1): &
                                          files(f)%dimensions(d)%bounds(2)) )
                   values_int16 = files(f)%dimensions(d)%values_int16
                   values_int16_pointer => values_int16
                   return_value = dom_write_var( files(f)%name, files(f)%dimensions(d)%name, &
                                     bounds_start=(/ files(f)%dimensions(d)%bounds(1) /),    &
                                     bounds_end  =(/ files(f)%dimensions(d)%bounds(2) /),    &
                                     values_int16_1d=values_int16_pointer )
                   DEALLOCATE( values_int16 )
                ELSEIF ( ALLOCATED( files(f)%dimensions(d)%values_int32 ) )  THEN
                   ALLOCATE( values_int32(files(f)%dimensions(d)%bounds(1): &
                                          files(f)%dimensions(d)%bounds(2)) )
                   values_int32 = files(f)%dimensions(d)%values_int32
                   values_int32_pointer => values_int32
                   return_value = dom_write_var( files(f)%name, files(f)%dimensions(d)%name, &
                                     bounds_start=(/ files(f)%dimensions(d)%bounds(1) /),    &
                                     bounds_end  =(/ files(f)%dimensions(d)%bounds(2) /),    &
                                     values_int32_1d=values_int32_pointer )
                   DEALLOCATE( values_int32 )
                ELSEIF ( ALLOCATED( files(f)%dimensions(d)%values_intwp ) )  THEN
                   ALLOCATE( values_intwp(files(f)%dimensions(d)%bounds(1): &
                                          files(f)%dimensions(d)%bounds(2)) )
                   values_intwp = files(f)%dimensions(d)%values_intwp
                   values_intwp_pointer => values_intwp
                   return_value = dom_write_var( files(f)%name, files(f)%dimensions(d)%name, &
                                     bounds_start=(/ files(f)%dimensions(d)%bounds(1) /),    &
                                     bounds_end  =(/ files(f)%dimensions(d)%bounds(2) /),    &
                                     values_intwp_1d=values_intwp_pointer )
                   DEALLOCATE( values_intwp )
                ELSEIF ( ALLOCATED( files(f)%dimensions(d)%values_real32 ) )  THEN
                   ALLOCATE( values_real32(files(f)%dimensions(d)%bounds(1): &
                                           files(f)%dimensions(d)%bounds(2)) )
                   values_real32 = files(f)%dimensions(d)%values_real32
                   values_real32_pointer => values_real32
                   return_value = dom_write_var( files(f)%name, files(f)%dimensions(d)%name, &
                                     bounds_start=(/ files(f)%dimensions(d)%bounds(1) /),    &
                                     bounds_end  =(/ files(f)%dimensions(d)%bounds(2) /),    &
                                     values_real32_1d=values_real32_pointer )
                   DEALLOCATE( values_real32 )
                ELSEIF ( ALLOCATED( files(f)%dimensions(d)%values_real64 ) )  THEN
                   ALLOCATE( values_real64(files(f)%dimensions(d)%bounds(1): &
                                           files(f)%dimensions(d)%bounds(2)) )
                   values_real64 = files(f)%dimensions(d)%values_real64
                   values_real64_pointer => values_real64
                   return_value = dom_write_var( files(f)%name, files(f)%dimensions(d)%name, &
                                     bounds_start=(/ files(f)%dimensions(d)%bounds(1) /),    &
                                     bounds_end  =(/ files(f)%dimensions(d)%bounds(2) /),    &
                                     values_real64_1d=values_real64_pointer )
                   DEALLOCATE( values_real64 )
                ELSEIF ( ALLOCATED( files(f)%dimensions(d)%values_realwp ) )  THEN
                   ALLOCATE( values_realwp(files(f)%dimensions(d)%bounds(1): &
                                           files(f)%dimensions(d)%bounds(2)) )
                   values_realwp = files(f)%dimensions(d)%values_realwp
                   values_realwp_pointer => values_realwp
                   return_value = dom_write_var( files(f)%name, files(f)%dimensions(d)%name, &
                                     bounds_start=(/ files(f)%dimensions(d)%bounds(1) /),    &
                                     bounds_end  =(/ files(f)%dimensions(d)%bounds(2) /),    &
                                     values_realwp_1d=values_realwp_pointer )
                   DEALLOCATE( values_realwp )
                ENDIF
                IF ( return_value /= 0 )  EXIT
             ENDDO

          ENDIF

          IF ( return_value /= 0 )  EXIT

       ENDDO
    ENDIF

    CALL internal_message( 'debug', routine_name // ': finished' )

 END FUNCTION dom_def_end

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write variable to file.
!> Example call:
!>   dom_write_var( file_name = 'my_output_file_name', &
!>                  name = 'u', &
!>                  bounds_start = (/nxl, nys, nzb, time_step/), &
!>                  bounds_end = (/nxr, nyn, nzt, time_step/), &
!>                  values_real64_3d = u )
!> @note The order of dimension bounds must match to the order of dimensions given in call
!>       'dom_def_var'. I.e., the corresponding variable definition should be like:
!>          dom_def_var( file_name =  'my_output_file_name', &
!>                       name = 'u', &
!>                       dimension_names = (/'x   ', 'y   ', 'z   ', 'time'/), &
!>                       output_type = <desired-output-type> )
!> @note The values given do not need to be of the same data type as was defined in the
!>       corresponding 'dom_def_var' call. If the output format 'netcdf' was chosen, the values are
!>       automatically converted to the data type given during the definition. If 'binary' was
!>       chosen, the values are written to file as given in the 'dom_write_var' call.
!--------------------------------------------------------------------------------------------------!
 FUNCTION dom_write_var( file_name, variable_name, bounds_start, bounds_end,         &
             values_char_0d,   values_char_1d,   values_char_2d,   values_char_3d,   &
             values_int8_0d,   values_int8_1d,   values_int8_2d,   values_int8_3d,   &
             values_int16_0d,  values_int16_1d,  values_int16_2d,  values_int16_3d,  &
             values_int32_0d,  values_int32_1d,  values_int32_2d,  values_int32_3d,  &
             values_intwp_0d,  values_intwp_1d,  values_intwp_2d,  values_intwp_3d,  &
             values_real32_0d, values_real32_1d, values_real32_2d, values_real32_3d, &
             values_real64_0d, values_real64_1d, values_real64_2d, values_real64_3d, &
             values_realwp_0d, values_realwp_1d, values_realwp_2d, values_realwp_3d  &
             ) RESULT( return_value )

    CHARACTER(LEN=charlen)            ::  file_format    !< file format chosen for file
    CHARACTER(LEN=*),      INTENT(IN) ::  file_name      !< name of file
    CHARACTER(LEN=*),      INTENT(IN) ::  variable_name  !< name of variable

    CHARACTER(LEN=*), PARAMETER ::  routine_name = 'dom_write_var'  !< name of routine

    CHARACTER(LEN=1), POINTER, INTENT(IN), OPTIONAL                   ::  values_char_0d           !< output variable
    CHARACTER(LEN=1), POINTER, INTENT(IN), OPTIONAL, DIMENSION(:)     ::  values_char_1d           !< output variable
    CHARACTER(LEN=1), POINTER, INTENT(IN), OPTIONAL, DIMENSION(:,:)   ::  values_char_2d           !< output variable
    CHARACTER(LEN=1), POINTER, INTENT(IN), OPTIONAL, DIMENSION(:,:,:) ::  values_char_3d           !< output variable

    CHARACTER(LEN=1), TARGET, ALLOCATABLE, DIMENSION(:)               ::  values_char_1d_resorted  !< resorted output variable
    CHARACTER(LEN=1), TARGET, ALLOCATABLE, DIMENSION(:,:)             ::  values_char_2d_resorted  !< resorted output variable
    CHARACTER(LEN=1), TARGET, ALLOCATABLE, DIMENSION(:,:,:)           ::  values_char_3d_resorted  !< resorted output variable

    CHARACTER(LEN=1), POINTER                                         ::  values_char_0d_pointer   !< pointer to resortet array
    CHARACTER(LEN=1), POINTER, CONTIGUOUS, DIMENSION(:)               ::  values_char_1d_pointer   !< pointer to resortet array
    CHARACTER(LEN=1), POINTER, CONTIGUOUS, DIMENSION(:,:)             ::  values_char_2d_pointer   !< pointer to resortet array
    CHARACTER(LEN=1), POINTER, CONTIGUOUS, DIMENSION(:,:,:)           ::  values_char_3d_pointer   !< pointer to resortet array

    INTEGER ::  file_id              !< file ID
    INTEGER ::  i                    !< loop index
    INTEGER ::  j                    !< loop index
    INTEGER ::  k                    !< loop index
    INTEGER ::  output_return_value  !< return value of a called output routine
    INTEGER ::  return_value         !< return value
    INTEGER ::  variable_id          !< variable ID

    INTEGER, DIMENSION(:),   INTENT(IN)  ::  bounds_end             !< end index per dimension of variable
    INTEGER, DIMENSION(:),   INTENT(IN)  ::  bounds_start           !< start index per dimension of variable
    INTEGER, DIMENSION(:),   ALLOCATABLE ::  bounds_origin          !< first index of each dimension
    INTEGER, DIMENSION(:),   ALLOCATABLE ::  bounds_start_internal  !< start index per dim. for output after masking
    INTEGER, DIMENSION(:),   ALLOCATABLE ::  value_counts           !< count of indices to be written per dimension
    INTEGER, DIMENSION(:,:), ALLOCATABLE ::  masked_indices         !< list containing all output indices along a dimension

    LOGICAL ::  do_output  !< true if any data lies within given range of masked dimension
    LOGICAL ::  is_global  !< true if variable is global

    INTEGER(KIND=1), POINTER, INTENT(IN), OPTIONAL                   ::  values_int8_0d             !< output variable
    INTEGER(KIND=2), POINTER, INTENT(IN), OPTIONAL                   ::  values_int16_0d            !< output variable
    INTEGER(KIND=4), POINTER, INTENT(IN), OPTIONAL                   ::  values_int32_0d            !< output variable
    INTEGER(iwp),    POINTER, INTENT(IN), OPTIONAL                   ::  values_intwp_0d            !< output variable
    INTEGER(KIND=1), POINTER, INTENT(IN), OPTIONAL, DIMENSION(:)     ::  values_int8_1d             !< output variable
    INTEGER(KIND=2), POINTER, INTENT(IN), OPTIONAL, DIMENSION(:)     ::  values_int16_1d            !< output variable
    INTEGER(KIND=4), POINTER, INTENT(IN), OPTIONAL, DIMENSION(:)     ::  values_int32_1d            !< output variable
    INTEGER(iwp),    POINTER, INTENT(IN), OPTIONAL, DIMENSION(:)     ::  values_intwp_1d            !< output variable
    INTEGER(KIND=1), POINTER, INTENT(IN), OPTIONAL, DIMENSION(:,:)   ::  values_int8_2d             !< output variable
    INTEGER(KIND=2), POINTER, INTENT(IN), OPTIONAL, DIMENSION(:,:)   ::  values_int16_2d            !< output variable
    INTEGER(KIND=4), POINTER, INTENT(IN), OPTIONAL, DIMENSION(:,:)   ::  values_int32_2d            !< output variable
    INTEGER(iwp),    POINTER, INTENT(IN), OPTIONAL, DIMENSION(:,:)   ::  values_intwp_2d            !< output variable
    INTEGER(KIND=1), POINTER, INTENT(IN), OPTIONAL, DIMENSION(:,:,:) ::  values_int8_3d             !< output variable
    INTEGER(KIND=2), POINTER, INTENT(IN), OPTIONAL, DIMENSION(:,:,:) ::  values_int16_3d            !< output variable
    INTEGER(KIND=4), POINTER, INTENT(IN), OPTIONAL, DIMENSION(:,:,:) ::  values_int32_3d            !< output variable
    INTEGER(iwp),    POINTER, INTENT(IN), OPTIONAL, DIMENSION(:,:,:) ::  values_intwp_3d            !< output variable

    INTEGER(KIND=1), TARGET, ALLOCATABLE, DIMENSION(:)               ::  values_int8_1d_resorted    !< resorted output variable
    INTEGER(KIND=2), TARGET, ALLOCATABLE, DIMENSION(:)               ::  values_int16_1d_resorted   !< resorted output variable
    INTEGER(KIND=4), TARGET, ALLOCATABLE, DIMENSION(:)               ::  values_int32_1d_resorted   !< resorted output variable
    INTEGER(iwp),    TARGET, ALLOCATABLE, DIMENSION(:)               ::  values_intwp_1d_resorted   !< resorted output variable
    INTEGER(KIND=1), TARGET, ALLOCATABLE, DIMENSION(:,:)             ::  values_int8_2d_resorted    !< resorted output variable
    INTEGER(KIND=2), TARGET, ALLOCATABLE, DIMENSION(:,:)             ::  values_int16_2d_resorted   !< resorted output variable
    INTEGER(KIND=4), TARGET, ALLOCATABLE, DIMENSION(:,:)             ::  values_int32_2d_resorted   !< resorted output variable
    INTEGER(iwp),    TARGET, ALLOCATABLE, DIMENSION(:,:)             ::  values_intwp_2d_resorted   !< resorted output variable
    INTEGER(KIND=1), TARGET, ALLOCATABLE, DIMENSION(:,:,:)           ::  values_int8_3d_resorted    !< resorted output variable
    INTEGER(KIND=2), TARGET, ALLOCATABLE, DIMENSION(:,:,:)           ::  values_int16_3d_resorted   !< resorted output variable
    INTEGER(KIND=4), TARGET, ALLOCATABLE, DIMENSION(:,:,:)           ::  values_int32_3d_resorted   !< resorted output variable
    INTEGER(iwp),    TARGET, ALLOCATABLE, DIMENSION(:,:,:)           ::  values_intwp_3d_resorted   !< resorted output variable

    INTEGER(KIND=1), POINTER                                         ::  values_int8_0d_pointer     !< pointer to resortet array
    INTEGER(KIND=2), POINTER                                         ::  values_int16_0d_pointer    !< pointer to resortet array
    INTEGER(KIND=4), POINTER                                         ::  values_int32_0d_pointer    !< pointer to resortet array
    INTEGER(iwp),    POINTER                                         ::  values_intwp_0d_pointer    !< pointer to resortet array
    INTEGER(KIND=1), POINTER, CONTIGUOUS, DIMENSION(:)               ::  values_int8_1d_pointer     !< pointer to resortet array
    INTEGER(KIND=2), POINTER, CONTIGUOUS, DIMENSION(:)               ::  values_int16_1d_pointer    !< pointer to resortet array
    INTEGER(KIND=4), POINTER, CONTIGUOUS, DIMENSION(:)               ::  values_int32_1d_pointer    !< pointer to resortet array
    INTEGER(iwp),    POINTER, CONTIGUOUS, DIMENSION(:)               ::  values_intwp_1d_pointer    !< pointer to resortet array
    INTEGER(KIND=1), POINTER, CONTIGUOUS, DIMENSION(:,:)             ::  values_int8_2d_pointer     !< pointer to resortet array
    INTEGER(KIND=2), POINTER, CONTIGUOUS, DIMENSION(:,:)             ::  values_int16_2d_pointer    !< pointer to resortet array
    INTEGER(KIND=4), POINTER, CONTIGUOUS, DIMENSION(:,:)             ::  values_int32_2d_pointer    !< pointer to resortet array
    INTEGER(iwp),    POINTER, CONTIGUOUS, DIMENSION(:,:)             ::  values_intwp_2d_pointer    !< pointer to resortet array
    INTEGER(KIND=1), POINTER, CONTIGUOUS, DIMENSION(:,:,:)           ::  values_int8_3d_pointer     !< pointer to resortet array
    INTEGER(KIND=2), POINTER, CONTIGUOUS, DIMENSION(:,:,:)           ::  values_int16_3d_pointer    !< pointer to resortet array
    INTEGER(KIND=4), POINTER, CONTIGUOUS, DIMENSION(:,:,:)           ::  values_int32_3d_pointer    !< pointer to resortet array
    INTEGER(iwp),    POINTER, CONTIGUOUS, DIMENSION(:,:,:)           ::  values_intwp_3d_pointer    !< pointer to resortet array

    REAL(KIND=4), POINTER, INTENT(IN), OPTIONAL                      ::  values_real32_0d           !< output variable
    REAL(KIND=8), POINTER, INTENT(IN), OPTIONAL                      ::  values_real64_0d           !< output variable
    REAL(wp),     POINTER, INTENT(IN), OPTIONAL                      ::  values_realwp_0d           !< output variable
    REAL(KIND=4), POINTER, INTENT(IN), OPTIONAL, DIMENSION(:)        ::  values_real32_1d           !< output variable
    REAL(KIND=8), POINTER, INTENT(IN), OPTIONAL, DIMENSION(:)        ::  values_real64_1d           !< output variable
    REAL(wp),     POINTER, INTENT(IN), OPTIONAL, DIMENSION(:)        ::  values_realwp_1d           !< output variable
    REAL(KIND=4), POINTER, INTENT(IN), OPTIONAL, DIMENSION(:,:)      ::  values_real32_2d           !< output variable
    REAL(KIND=8), POINTER, INTENT(IN), OPTIONAL, DIMENSION(:,:)      ::  values_real64_2d           !< output variable
    REAL(wp),     POINTER, INTENT(IN), OPTIONAL, DIMENSION(:,:)      ::  values_realwp_2d           !< output variable
    REAL(KIND=4), POINTER, INTENT(IN), OPTIONAL, DIMENSION(:,:,:)    ::  values_real32_3d           !< output variable
    REAL(KIND=8), POINTER, INTENT(IN), OPTIONAL, DIMENSION(:,:,:)    ::  values_real64_3d           !< output variable
    REAL(wp),     POINTER, INTENT(IN), OPTIONAL, DIMENSION(:,:,:)    ::  values_realwp_3d           !< output variable

    REAL(KIND=4), TARGET, ALLOCATABLE, DIMENSION(:)                  ::  values_real32_1d_resorted  !< resorted output variable
    REAL(KIND=8), TARGET, ALLOCATABLE, DIMENSION(:)                  ::  values_real64_1d_resorted  !< resorted output variable
    REAL(wp),     TARGET, ALLOCATABLE, DIMENSION(:)                  ::  values_realwp_1d_resorted  !< resorted output variable
    REAL(KIND=4), TARGET, ALLOCATABLE, DIMENSION(:,:)                ::  values_real32_2d_resorted  !< resorted output variable
    REAL(KIND=8), TARGET, ALLOCATABLE, DIMENSION(:,:)                ::  values_real64_2d_resorted  !< resorted output variable
    REAL(wp),     TARGET, ALLOCATABLE, DIMENSION(:,:)                ::  values_realwp_2d_resorted  !< resorted output variable
    REAL(KIND=4), TARGET, ALLOCATABLE, DIMENSION(:,:,:)              ::  values_real32_3d_resorted  !< resorted output variable
    REAL(KIND=8), TARGET, ALLOCATABLE, DIMENSION(:,:,:)              ::  values_real64_3d_resorted  !< resorted output variable
    REAL(wp),     TARGET, ALLOCATABLE, DIMENSION(:,:,:)              ::  values_realwp_3d_resorted  !< resorted output variable

    REAL(KIND=4), POINTER                                            ::  values_real32_0d_pointer   !< pointer to resortet array
    REAL(KIND=8), POINTER                                            ::  values_real64_0d_pointer   !< pointer to resortet array
    REAL(wp),     POINTER                                            ::  values_realwp_0d_pointer   !< pointer to resortet array
    REAL(KIND=4), POINTER, CONTIGUOUS, DIMENSION(:)                  ::  values_real32_1d_pointer   !< pointer to resortet array
    REAL(KIND=8), POINTER, CONTIGUOUS, DIMENSION(:)                  ::  values_real64_1d_pointer   !< pointer to resortet array
    REAL(wp),     POINTER, CONTIGUOUS, DIMENSION(:)                  ::  values_realwp_1d_pointer   !< pointer to resortet array
    REAL(KIND=4), POINTER, CONTIGUOUS, DIMENSION(:,:)                ::  values_real32_2d_pointer   !< pointer to resortet array
    REAL(KIND=8), POINTER, CONTIGUOUS, DIMENSION(:,:)                ::  values_real64_2d_pointer   !< pointer to resortet array
    REAL(wp),     POINTER, CONTIGUOUS, DIMENSION(:,:)                ::  values_realwp_2d_pointer   !< pointer to resortet array
    REAL(KIND=4), POINTER, CONTIGUOUS, DIMENSION(:,:,:)              ::  values_real32_3d_pointer   !< pointer to resortet array
    REAL(KIND=8), POINTER, CONTIGUOUS, DIMENSION(:,:,:)              ::  values_real64_3d_pointer   !< pointer to resortet array
    REAL(wp),     POINTER, CONTIGUOUS, DIMENSION(:,:,:)              ::  values_realwp_3d_pointer   !< pointer to resortet array

    TYPE(dimension_type), DIMENSION(:), ALLOCATABLE ::  dimension_list  !< list of used dimensions of variable


    return_value = 0
    output_return_value = 0

    CALL internal_message( 'debug', routine_name // ': write ' // TRIM( variable_name ) // &
                           ' into file ' // TRIM( file_name ) )
!
!-- Search for variable within file
    CALL find_var_in_file( file_name, variable_name, file_format, file_id, variable_id, &
                           is_global, dimension_list, return_value=return_value  )

    IF ( return_value == 0 )  THEN
!
!--    Check if the correct amount of variable bounds were given
       IF ( SIZE( bounds_start ) /= SIZE( dimension_list )  .OR.  &
            SIZE( bounds_end ) /= SIZE( dimension_list ) )  THEN
          return_value = 1
          CALL internal_message( 'error', routine_name //                  &
                                 ': number bounds do not match with ' //   &
                                 'number of dimensions of variable ' //    &
                                 '(variable "' // TRIM( variable_name ) // &
                                 '", file "' // TRIM( file_name ) // '")!' )
       ENDIF

    ENDIF

    IF ( return_value == 0 )  THEN
!
!--    Save starting index (lower bounds) of each dimension
       ALLOCATE( bounds_origin(SIZE( dimension_list )) )
       ALLOCATE( bounds_start_internal(SIZE( dimension_list )) )
       ALLOCATE( value_counts(SIZE( dimension_list )) )

       WRITE( temp_string, * ) bounds_start
       CALL internal_message( 'debug', routine_name //                    &
                              ': file "' // TRIM( file_name ) //          &
                              '", variable "' // TRIM( variable_name ) // &
                              '", bounds_start =' // TRIM( temp_string ) )
       WRITE( temp_string, * ) bounds_end
       CALL internal_message( 'debug', routine_name //                    &
                              ': file "' // TRIM( file_name ) //          &
                              '", variable "' // TRIM( variable_name ) // &
                              '", bounds_end =' // TRIM( temp_string ) )
!
!--    Get bounds for masking
       CALL get_masked_indices_and_masked_dimension_bounds( dimension_list,                  &
               bounds_start, bounds_end, bounds_start_internal, value_counts, bounds_origin, &
               masked_indices )

       do_output = .NOT. ANY( value_counts == 0 )

       WRITE( temp_string, * ) bounds_start_internal
       CALL internal_message( 'debug', routine_name //                    &
                              ': file "' // TRIM( file_name ) //          &
                              '", variable "' // TRIM( variable_name ) // &
                              '", bounds_start_internal =' // TRIM( temp_string ) )
       WRITE( temp_string, * ) value_counts
       CALL internal_message( 'debug', routine_name //                    &
                              ': file "' // TRIM( file_name ) //          &
                              '", variable "' // TRIM( variable_name ) // &
                              '", value_counts =' // TRIM( temp_string ) )
!
!--    Mask and resort variable
!--    character output
       IF ( PRESENT( values_char_0d ) )  THEN
          values_char_0d_pointer => values_char_0d
       ELSEIF ( PRESENT( values_char_1d ) )  THEN
          IF ( do_output ) THEN
             ALLOCATE( values_char_1d_resorted(0:value_counts(1)-1) )
             !$OMP PARALLEL PRIVATE (i)
             !$OMP DO
             DO  i = 0, value_counts(1) - 1
                values_char_1d_resorted(i) = values_char_1d(masked_indices(1,i))
             ENDDO
             !$OMP END PARALLEL
          ELSE
             ALLOCATE( values_char_1d_resorted(1) )
             values_char_1d_resorted = ' '
          ENDIF
          values_char_1d_pointer => values_char_1d_resorted
       ELSEIF ( PRESENT( values_char_2d ) )  THEN
          IF ( do_output ) THEN
             ALLOCATE( values_char_2d_resorted(0:value_counts(1)-1, &
                                               0:value_counts(2)-1) )
             !$OMP PARALLEL PRIVATE (i,j)
             !$OMP DO
             DO  i = 0, value_counts(1) - 1
                DO  j = 0, value_counts(2) - 1
                   values_char_2d_resorted(i,j) = values_char_2d(masked_indices(2,j), &
                                                                 masked_indices(1,i)  )
                ENDDO
             ENDDO
             !$OMP END PARALLEL
          ELSE
             ALLOCATE( values_char_2d_resorted(1,1) )
             values_char_2d_resorted = ' '
          ENDIF
          values_char_2d_pointer => values_char_2d_resorted
       ELSEIF ( PRESENT( values_char_3d ) )  THEN
          IF ( do_output ) THEN
             ALLOCATE( values_char_3d_resorted(0:value_counts(1)-1, &
                                               0:value_counts(2)-1, &
                                               0:value_counts(3)-1) )
             !$OMP PARALLEL PRIVATE (i,j,k)
             !$OMP DO
             DO  i = 0, value_counts(1) - 1
                DO  j = 0, value_counts(2) - 1
                   DO  k = 0, value_counts(3) - 1
                      values_char_3d_resorted(i,j,k) = values_char_3d(masked_indices(3,k), &
                                                                      masked_indices(2,j), &
                                                                      masked_indices(1,i)  )
                   ENDDO
                ENDDO
             ENDDO
             !$OMP END PARALLEL
          ELSE
             ALLOCATE( values_char_3d_resorted(1,1,1) )
             values_char_3d_resorted = ' '
          ENDIF
          values_char_3d_pointer => values_char_3d_resorted
!
!--    8bit integer output
       ELSEIF ( PRESENT( values_int8_0d ) )  THEN
          values_int8_0d_pointer => values_int8_0d
       ELSEIF ( PRESENT( values_int8_1d ) )  THEN
          IF ( do_output ) THEN
             ALLOCATE( values_int8_1d_resorted(0:value_counts(1)-1) )
             !$OMP PARALLEL PRIVATE (i)
             !$OMP DO
             DO  i = 0, value_counts(1) - 1
                values_int8_1d_resorted(i) = values_int8_1d(masked_indices(1,i))
             ENDDO
             !$OMP END PARALLEL
          ELSE
             ALLOCATE( values_int8_1d_resorted(1) )
             values_int8_1d_resorted = 0_1
          ENDIF
          values_int8_1d_pointer => values_int8_1d_resorted
       ELSEIF ( PRESENT( values_int8_2d ) )  THEN
          IF ( do_output ) THEN
             ALLOCATE( values_int8_2d_resorted(0:value_counts(1)-1, &
                                               0:value_counts(2)-1) )
             !$OMP PARALLEL PRIVATE (i,j)
             !$OMP DO
             DO  i = 0, value_counts(1) - 1
                DO  j = 0, value_counts(2) - 1
                   values_int8_2d_resorted(i,j) = values_int8_2d(masked_indices(2,j), &
                                                                 masked_indices(1,i)  )
                ENDDO
             ENDDO
             !$OMP END PARALLEL
          ELSE
             ALLOCATE( values_int8_2d_resorted(1,1) )
             values_int8_2d_resorted = 0_1
          ENDIF
          values_int8_2d_pointer => values_int8_2d_resorted
       ELSEIF ( PRESENT( values_int8_3d ) )  THEN
          IF ( do_output ) THEN
             ALLOCATE( values_int8_3d_resorted(0:value_counts(1)-1, &
                                               0:value_counts(2)-1, &
                                               0:value_counts(3)-1) )
             !$OMP PARALLEL PRIVATE (i,j,k)
             !$OMP DO
             DO  i = 0, value_counts(1) - 1
                DO  j = 0, value_counts(2) - 1
                   DO  k = 0, value_counts(3) - 1
                      values_int8_3d_resorted(i,j,k) = values_int8_3d(masked_indices(3,k), &
                                                                      masked_indices(2,j), &
                                                                      masked_indices(1,i)  )
                   ENDDO
                ENDDO
             ENDDO
             !$OMP END PARALLEL
          ELSE
             ALLOCATE( values_int8_3d_resorted(1,1,1) )
             values_int8_3d_resorted = 0_1
          ENDIF
          values_int8_3d_pointer => values_int8_3d_resorted
!
!--    16bit integer output
       ELSEIF ( PRESENT( values_int16_0d ) )  THEN
          values_int16_0d_pointer => values_int16_0d
       ELSEIF ( PRESENT( values_int16_1d ) )  THEN
          IF ( do_output ) THEN
             ALLOCATE( values_int16_1d_resorted(0:value_counts(1)-1) )
             !$OMP PARALLEL PRIVATE (i)
             !$OMP DO
             DO  i = 0, value_counts(1) - 1
                values_int16_1d_resorted(i) = values_int16_1d(masked_indices(1,i))
             ENDDO
             !$OMP END PARALLEL
          ELSE
             ALLOCATE( values_int16_1d_resorted(1) )
             values_int16_1d_resorted = 0_1
          ENDIF
          values_int16_1d_pointer => values_int16_1d_resorted
       ELSEIF ( PRESENT( values_int16_2d ) )  THEN
          IF ( do_output ) THEN
             ALLOCATE( values_int16_2d_resorted(0:value_counts(1)-1, &
                                                0:value_counts(2)-1) )
             !$OMP PARALLEL PRIVATE (i,j)
             !$OMP DO
             DO  i = 0, value_counts(1) - 1
                DO  j = 0, value_counts(2) - 1
                   values_int16_2d_resorted(i,j) = values_int16_2d(masked_indices(2,j), &
                                                                   masked_indices(1,i))
                ENDDO
             ENDDO
             !$OMP END PARALLEL
          ELSE
             ALLOCATE( values_int16_2d_resorted(1,1) )
             values_int16_2d_resorted = 0_1
          ENDIF
          values_int16_2d_pointer => values_int16_2d_resorted
       ELSEIF ( PRESENT( values_int16_3d ) )  THEN
          IF ( do_output ) THEN
             ALLOCATE( values_int16_3d_resorted(0:value_counts(1)-1, &
                                                0:value_counts(2)-1, &
                                                0:value_counts(3)-1) )
             !$OMP PARALLEL PRIVATE (i,j,k)
             !$OMP DO
             DO  i = 0, value_counts(1) - 1
                DO  j = 0, value_counts(2) - 1
                   DO  k = 0, value_counts(3) - 1
                      values_int16_3d_resorted(i,j,k) = values_int16_3d(masked_indices(3,k), &
                                                                        masked_indices(2,j), &
                                                                        masked_indices(1,i)  )
                   ENDDO
                ENDDO
             ENDDO
             !$OMP END PARALLEL
          ELSE
             ALLOCATE( values_int16_3d_resorted(1,1,1) )
             values_int16_3d_resorted = 0_1
          ENDIF
          values_int16_3d_pointer => values_int16_3d_resorted
!
!--    32bit integer output
       ELSEIF ( PRESENT( values_int32_0d ) )  THEN
          values_int32_0d_pointer => values_int32_0d
       ELSEIF ( PRESENT( values_int32_1d ) )  THEN
          IF ( do_output ) THEN
             ALLOCATE( values_int32_1d_resorted(0:value_counts(1)-1) )
             !$OMP PARALLEL PRIVATE (i)
             !$OMP DO
             DO  i = 0, value_counts(1) - 1
                values_int32_1d_resorted(i) = values_int32_1d(masked_indices(1,i))
             ENDDO
             !$OMP END PARALLEL
          ELSE
             ALLOCATE( values_int32_1d_resorted(1) )
             values_int32_1d_resorted = 0_1
          ENDIF
          values_int32_1d_pointer => values_int32_1d_resorted
       ELSEIF ( PRESENT( values_int32_2d ) )  THEN
          IF ( do_output ) THEN
             ALLOCATE( values_int32_2d_resorted(0:value_counts(1)-1, &
                                                0:value_counts(2)-1) )
             !$OMP PARALLEL PRIVATE (i,j)
             !$OMP DO
             DO  i = 0, value_counts(1) - 1
                DO  j = 0, value_counts(2) - 1
                   values_int32_2d_resorted(i,j) = values_int32_2d(masked_indices(2,j), &
                                                                   masked_indices(1,i)  )
                ENDDO
             ENDDO
             !$OMP END PARALLEL
          ELSE
             ALLOCATE( values_int32_2d_resorted(1,1) )
             values_int32_2d_resorted = 0_1
          ENDIF
          values_int32_2d_pointer => values_int32_2d_resorted
       ELSEIF ( PRESENT( values_int32_3d ) )  THEN
          IF ( do_output ) THEN
             ALLOCATE( values_int32_3d_resorted(0:value_counts(1)-1, &
                                                0:value_counts(2)-1, &
                                                0:value_counts(3)-1) )
             !$OMP PARALLEL PRIVATE (i,j,k)
             !$OMP DO
             DO  i = 0, value_counts(1) - 1
                DO  j = 0, value_counts(2) - 1
                   DO  k = 0, value_counts(3) - 1
                      values_int32_3d_resorted(i,j,k) = values_int32_3d(masked_indices(3,k), &
                                                                        masked_indices(2,j), &
                                                                        masked_indices(1,i)  )
                   ENDDO
                ENDDO
             ENDDO
             !$OMP END PARALLEL
          ELSE
             ALLOCATE( values_int32_3d_resorted(1,1,1) )
             values_int32_3d_resorted = 0_1
          ENDIF
          values_int32_3d_pointer => values_int32_3d_resorted
!
!--    working-precision integer output
       ELSEIF ( PRESENT( values_intwp_0d ) )  THEN
          values_intwp_0d_pointer => values_intwp_0d
       ELSEIF ( PRESENT( values_intwp_1d ) )  THEN
          IF ( do_output ) THEN
             ALLOCATE( values_intwp_1d_resorted(0:value_counts(1)-1) )
             !$OMP PARALLEL PRIVATE (i)
             !$OMP DO
             DO  i = 0, value_counts(1) - 1
                values_intwp_1d_resorted(i) = values_intwp_1d(masked_indices(1,i))
             ENDDO
             !$OMP END PARALLEL
          ELSE
             ALLOCATE( values_intwp_1d_resorted(1) )
             values_intwp_1d_resorted = 0_1
          ENDIF
          values_intwp_1d_pointer => values_intwp_1d_resorted
       ELSEIF ( PRESENT( values_intwp_2d ) )  THEN
          IF ( do_output ) THEN
             ALLOCATE( values_intwp_2d_resorted(0:value_counts(1)-1, &
                                                0:value_counts(2)-1) )
             !$OMP PARALLEL PRIVATE (i,j)
             !$OMP DO
             DO  i = 0, value_counts(1) - 1
                DO  j = 0, value_counts(2) - 1
                   values_intwp_2d_resorted(i,j) = values_intwp_2d(masked_indices(2,j), &
                                                                   masked_indices(1,i)  )
                ENDDO
             ENDDO
             !$OMP END PARALLEL
          ELSE
             ALLOCATE( values_intwp_2d_resorted(1,1) )
             values_intwp_2d_resorted = 0_1
          ENDIF
          values_intwp_2d_pointer => values_intwp_2d_resorted
       ELSEIF ( PRESENT( values_intwp_3d ) )  THEN
          IF ( do_output ) THEN
             ALLOCATE( values_intwp_3d_resorted(0:value_counts(1)-1, &
                                                0:value_counts(2)-1, &
                                                0:value_counts(3)-1) )
             !$OMP PARALLEL PRIVATE (i,j,k)
             !$OMP DO
             DO  i = 0, value_counts(1) - 1
                DO  j = 0, value_counts(2) - 1
                   DO  k = 0, value_counts(3) - 1
                      values_intwp_3d_resorted(i,j,k) = values_intwp_3d(masked_indices(3,k), &
                                                                        masked_indices(2,j), &
                                                                        masked_indices(1,i)  )
                   ENDDO
                ENDDO
             ENDDO
             !$OMP END PARALLEL
          ELSE
             ALLOCATE( values_intwp_3d_resorted(1,1,1) )
             values_intwp_3d_resorted = 0_1
          ENDIF
          values_intwp_3d_pointer => values_intwp_3d_resorted
!
!--    32bit real output
       ELSEIF ( PRESENT( values_real32_0d ) )  THEN
          values_real32_0d_pointer => values_real32_0d
       ELSEIF ( PRESENT( values_real32_1d ) )  THEN
          IF ( do_output ) THEN
             ALLOCATE( values_real32_1d_resorted(0:value_counts(1)-1) )
             !$OMP PARALLEL PRIVATE (i)
             !$OMP DO
             DO  i = 0, value_counts(1) - 1
                values_real32_1d_resorted(i) = values_real32_1d(masked_indices(1,i))
             ENDDO
             !$OMP END PARALLEL
          ELSE
             ALLOCATE( values_real32_1d_resorted(1) )
             values_real32_1d_resorted = 0_1
          ENDIF
          values_real32_1d_pointer => values_real32_1d_resorted
       ELSEIF ( PRESENT( values_real32_2d ) )  THEN
          IF ( do_output ) THEN
             ALLOCATE( values_real32_2d_resorted(0:value_counts(1)-1, &
                                                 0:value_counts(2)-1) )
             !$OMP PARALLEL PRIVATE (i,j)
             !$OMP DO
             DO  i = 0, value_counts(1) - 1
                DO  j = 0, value_counts(2) - 1
                   values_real32_2d_resorted(i,j) = values_real32_2d(masked_indices(2,j), &
                                                                     masked_indices(1,i)  )
                ENDDO
             ENDDO
             !$OMP END PARALLEL
          ELSE
             ALLOCATE( values_real32_2d_resorted(1,1) )
             values_real32_2d_resorted = 0_1
          ENDIF
          values_real32_2d_pointer => values_real32_2d_resorted
       ELSEIF ( PRESENT( values_real32_3d ) )  THEN
          IF ( do_output ) THEN
             ALLOCATE( values_real32_3d_resorted(0:value_counts(1)-1, &
                                                 0:value_counts(2)-1, &
                                                 0:value_counts(3)-1) )
             !$OMP PARALLEL PRIVATE (i,j,k)
             !$OMP DO
             DO  i = 0, value_counts(1) - 1
                DO  j = 0, value_counts(2) - 1
                   DO  k = 0, value_counts(3) - 1
                      values_real32_3d_resorted(i,j,k) = values_real32_3d(masked_indices(3,k), &
                                                                          masked_indices(2,j), &
                                                                          masked_indices(1,i)  )
                   ENDDO
                ENDDO
             ENDDO
             !$OMP END PARALLEL
          ELSE
             ALLOCATE( values_real32_3d_resorted(1,1,1) )
             values_real32_3d_resorted = 0_1
          ENDIF
          values_real32_3d_pointer => values_real32_3d_resorted
!
!--    64bit real output
       ELSEIF ( PRESENT( values_real64_0d ) )  THEN
          values_real64_0d_pointer => values_real64_0d
       ELSEIF ( PRESENT( values_real64_1d ) )  THEN
          IF ( do_output ) THEN
             ALLOCATE( values_real64_1d_resorted(0:value_counts(1)-1) )
             !$OMP PARALLEL PRIVATE (i)
             !$OMP DO
             DO  i = 0, value_counts(1) - 1
                values_real64_1d_resorted(i) = values_real64_1d(masked_indices(1,i))
             ENDDO
             !$OMP END PARALLEL
          ELSE
             ALLOCATE( values_real64_1d_resorted(1) )
             values_real64_1d_resorted = 0_1
          ENDIF
          values_real64_1d_pointer => values_real64_1d_resorted
       ELSEIF ( PRESENT( values_real64_2d ) )  THEN
          IF ( do_output ) THEN
             ALLOCATE( values_real64_2d_resorted(0:value_counts(1)-1, &
                                                 0:value_counts(2)-1) )
             !$OMP PARALLEL PRIVATE (i,j)
             !$OMP DO
             DO  i = 0, value_counts(1) - 1
                DO  j = 0, value_counts(2) - 1
                   values_real64_2d_resorted(i,j) = values_real64_2d(masked_indices(2,j), &
                                                                     masked_indices(1,i)  )
                ENDDO
             ENDDO
             !$OMP END PARALLEL
          ELSE
             ALLOCATE( values_real64_2d_resorted(1,1) )
             values_real64_2d_resorted = 0_1
          ENDIF
          values_real64_2d_pointer => values_real64_2d_resorted
       ELSEIF ( PRESENT( values_real64_3d ) )  THEN
          IF ( do_output ) THEN
             ALLOCATE( values_real64_3d_resorted(0:value_counts(1)-1, &
                                                 0:value_counts(2)-1, &
                                                 0:value_counts(3)-1) )
             !$OMP PARALLEL PRIVATE (i,j,k)
             !$OMP DO
             DO  i = 0, value_counts(1) - 1
                DO  j = 0, value_counts(2) - 1
                   DO  k = 0, value_counts(3) - 1
                      values_real64_3d_resorted(i,j,k) = values_real64_3d(masked_indices(3,k), &
                                                                          masked_indices(2,j), &
                                                                          masked_indices(1,i)  )
                   ENDDO
                ENDDO
             ENDDO
             !$OMP END PARALLEL
          ELSE
             ALLOCATE( values_real64_3d_resorted(1,1,1) )
             values_real64_3d_resorted = 0_1
          ENDIF
          values_real64_3d_pointer => values_real64_3d_resorted
!
!--    working-precision real output
       ELSEIF ( PRESENT( values_realwp_0d ) )  THEN
          values_realwp_0d_pointer => values_realwp_0d
       ELSEIF ( PRESENT( values_realwp_1d ) )  THEN
          IF ( do_output ) THEN
             ALLOCATE( values_realwp_1d_resorted(0:value_counts(1)-1) )
             !$OMP PARALLEL PRIVATE (i)
             !$OMP DO
             DO  i = 0, value_counts(1) - 1
                values_realwp_1d_resorted(i) = values_realwp_1d(masked_indices(1,i))
             ENDDO
             !$OMP END PARALLEL
          ELSE
             ALLOCATE( values_realwp_1d_resorted(1) )
             values_realwp_1d_resorted = 0_1
          ENDIF
          values_realwp_1d_pointer => values_realwp_1d_resorted
       ELSEIF ( PRESENT( values_realwp_2d ) )  THEN
          IF ( do_output ) THEN
             ALLOCATE( values_realwp_2d_resorted(0:value_counts(1)-1, &
                                                 0:value_counts(2)-1) )
             !$OMP PARALLEL PRIVATE (i,j)
             !$OMP DO
             DO  i = 0, value_counts(1) - 1
                DO  j = 0, value_counts(2) - 1
                   values_realwp_2d_resorted(i,j) = values_realwp_2d(masked_indices(2,j), &
                                                                     masked_indices(1,i)  )
                ENDDO
             ENDDO
             !$OMP END PARALLEL
          ELSE
             ALLOCATE( values_realwp_2d_resorted(1,1) )
             values_realwp_2d_resorted = 0_1
          ENDIF
          values_realwp_2d_pointer => values_realwp_2d_resorted
       ELSEIF ( PRESENT( values_realwp_3d ) )  THEN
          IF ( do_output ) THEN
             ALLOCATE( values_realwp_3d_resorted(0:value_counts(1)-1, &
                                                 0:value_counts(2)-1, &
                                                 0:value_counts(3)-1) )
             !$OMP PARALLEL PRIVATE (i,j,k)
             !$OMP DO
             DO  i = 0, value_counts(1) - 1
                DO  j = 0, value_counts(2) - 1
                   DO  k = 0, value_counts(3) - 1
                      values_realwp_3d_resorted(i,j,k) = values_realwp_3d(masked_indices(3,k), &
                                                                          masked_indices(2,j), &
                                                                          masked_indices(1,i)  )
                   ENDDO
                ENDDO
             ENDDO
             !$OMP END PARALLEL
          ELSE
             ALLOCATE( values_realwp_3d_resorted(1,1,1) )
             values_realwp_3d_resorted = 0_1
          ENDIF
          values_realwp_3d_pointer => values_realwp_3d_resorted

       ELSE
          return_value = 1
          CALL internal_message( 'error', routine_name //                  &
                                 ': no output values given ' //            &
                                 '(variable "' // TRIM( variable_name ) // &
                                 '", file "' // TRIM( file_name ) // '")!'  )
       ENDIF

       DEALLOCATE( masked_indices )

    ENDIF  ! Check for error

    IF ( return_value == 0 )  THEN
!
!--    Write variable into file
       SELECT CASE ( TRIM( file_format ) )

          CASE ( 'binary' )
!
!--          character output
             IF ( PRESENT( values_char_0d ) )  THEN
                 CALL binary_write_variable( file_id, variable_id,                      &
                         bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_char_0d=values_char_0d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_char_1d ) )  THEN
                CALL binary_write_variable( file_id, variable_id,                      &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_char_1d=values_char_1d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_char_2d ) )  THEN
                CALL binary_write_variable( file_id, variable_id,                      &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_char_2d=values_char_2d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_char_3d ) )  THEN
                CALL binary_write_variable( file_id, variable_id,                      &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_char_3d=values_char_3d_pointer, return_value=output_return_value )
!
!--          8bit integer output
             ELSEIF ( PRESENT( values_int8_0d ) )  THEN
                CALL binary_write_variable( file_id, variable_id,                      &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_int8_0d=values_int8_0d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_int8_1d ) )  THEN
                CALL binary_write_variable( file_id, variable_id,                      &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_int8_1d=values_int8_1d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_int8_2d ) )  THEN
                CALL binary_write_variable( file_id, variable_id,                      &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_int8_2d=values_int8_2d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_int8_3d ) )  THEN
                CALL binary_write_variable( file_id, variable_id,                      &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_int8_3d=values_int8_3d_pointer, return_value=output_return_value )
!
!--          16bit integer output
             ELSEIF ( PRESENT( values_int16_0d ) )  THEN
                CALL binary_write_variable( file_id, variable_id,                      &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_int16_0d=values_int16_0d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_int16_1d ) )  THEN
                CALL binary_write_variable( file_id, variable_id,                      &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_int16_1d=values_int16_1d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_int16_2d ) )  THEN
                CALL binary_write_variable( file_id, variable_id,                      &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_int16_2d=values_int16_2d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_int16_3d ) )  THEN
                CALL binary_write_variable( file_id, variable_id,                      &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_int16_3d=values_int16_3d_pointer, return_value=output_return_value )
!
!--          32bit integer output
             ELSEIF ( PRESENT( values_int32_0d ) )  THEN
                CALL binary_write_variable( file_id, variable_id,                      &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_int32_0d=values_int32_0d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_int32_1d ) )  THEN
                CALL binary_write_variable( file_id, variable_id,                      &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_int32_1d=values_int32_1d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_int32_2d ) )  THEN
                CALL binary_write_variable( file_id, variable_id,                      &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_int32_2d=values_int32_2d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_int32_3d ) )  THEN
                CALL binary_write_variable( file_id, variable_id,                      &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_int32_3d=values_int32_3d_pointer, return_value=output_return_value )
!
!--          working-precision integer output
             ELSEIF ( PRESENT( values_intwp_0d ) )  THEN
                CALL binary_write_variable( file_id, variable_id,                      &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_intwp_0d=values_intwp_0d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_intwp_1d ) )  THEN
                CALL binary_write_variable( file_id, variable_id,                      &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_intwp_1d=values_intwp_1d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_intwp_2d ) )  THEN
                CALL binary_write_variable( file_id, variable_id,                      &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_intwp_2d=values_intwp_2d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_intwp_3d ) )  THEN
                CALL binary_write_variable( file_id, variable_id,                      &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_intwp_3d=values_intwp_3d_pointer, return_value=output_return_value )
!
!--          32bit real output
             ELSEIF ( PRESENT( values_real32_0d ) )  THEN
                CALL binary_write_variable( file_id, variable_id,                      &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_real32_0d=values_real32_0d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_real32_1d ) )  THEN
                CALL binary_write_variable( file_id, variable_id,                      &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_real32_1d=values_real32_1d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_real32_2d ) )  THEN
                CALL binary_write_variable( file_id, variable_id,                      &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_real32_2d=values_real32_2d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_real32_3d ) )  THEN
                CALL binary_write_variable( file_id, variable_id,                      &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_real32_3d=values_real32_3d_pointer, return_value=output_return_value )
!
!--          64bit real output
             ELSEIF ( PRESENT( values_real64_0d ) )  THEN
                CALL binary_write_variable( file_id, variable_id,                      &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_real64_0d=values_real64_0d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_real64_1d ) )  THEN
                CALL binary_write_variable( file_id, variable_id,                      &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_real64_1d=values_real64_1d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_real64_2d ) )  THEN
                CALL binary_write_variable( file_id, variable_id,                      &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_real64_2d=values_real64_2d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_real64_3d ) )  THEN
                CALL binary_write_variable( file_id, variable_id,                      &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_real64_3d=values_real64_3d_pointer, return_value=output_return_value )
!
!--          working-precision real output
             ELSEIF ( PRESENT( values_realwp_0d ) )  THEN
                CALL binary_write_variable( file_id, variable_id,                      &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_realwp_0d=values_realwp_0d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_realwp_1d ) )  THEN
                CALL binary_write_variable( file_id, variable_id,                      &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_realwp_1d=values_realwp_1d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_realwp_2d ) )  THEN
                CALL binary_write_variable( file_id, variable_id,                      &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_realwp_2d=values_realwp_2d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_realwp_3d ) )  THEN
                CALL binary_write_variable( file_id, variable_id,                      &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_realwp_3d=values_realwp_3d_pointer, return_value=output_return_value )
             ELSE
                return_value = 1
                CALL internal_message( 'error', routine_name //                          &
                                       ': output_type not supported by file format "' // &
                                       TRIM( file_format ) // '" ' //                    &
                                       '(variable "' // TRIM( variable_name ) //         &
                                       '", file "' // TRIM( file_name ) // '")!' )
             ENDIF

          CASE ( 'netcdf4-parallel', 'netcdf4-serial' )
!
!--          character output
             IF ( PRESENT( values_char_0d ) )  THEN
                 CALL netcdf4_write_variable( file_id, variable_id,                     &
                         bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_char_0d=values_char_0d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_char_1d ) )  THEN
                CALL netcdf4_write_variable( file_id, variable_id,                     &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_char_1d=values_char_1d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_char_2d ) )  THEN
                CALL netcdf4_write_variable( file_id, variable_id,                     &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_char_2d=values_char_2d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_char_3d ) )  THEN
                CALL netcdf4_write_variable( file_id, variable_id,                     &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_char_3d=values_char_3d_pointer, return_value=output_return_value )
!
!--          8bit integer output
             ELSEIF ( PRESENT( values_int8_0d ) )  THEN
                CALL netcdf4_write_variable( file_id, variable_id,                     &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_int8_0d=values_int8_0d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_int8_1d ) )  THEN
                CALL netcdf4_write_variable( file_id, variable_id,                     &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_int8_1d=values_int8_1d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_int8_2d ) )  THEN
                CALL netcdf4_write_variable( file_id, variable_id,                     &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_int8_2d=values_int8_2d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_int8_3d ) )  THEN
                CALL netcdf4_write_variable( file_id, variable_id,                     &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_int8_3d=values_int8_3d_pointer, return_value=output_return_value )
!
!--          16bit integer output
             ELSEIF ( PRESENT( values_int16_0d ) )  THEN
                CALL netcdf4_write_variable( file_id, variable_id,                     &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_int16_0d=values_int16_0d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_int16_1d ) )  THEN
                CALL netcdf4_write_variable( file_id, variable_id,                     &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_int16_1d=values_int16_1d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_int16_2d ) )  THEN
                CALL netcdf4_write_variable( file_id, variable_id,                     &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_int16_2d=values_int16_2d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_int16_3d ) )  THEN
                CALL netcdf4_write_variable( file_id, variable_id,                     &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_int16_3d=values_int16_3d_pointer, return_value=output_return_value )
!
!--          32bit integer output
             ELSEIF ( PRESENT( values_int32_0d ) )  THEN
                CALL netcdf4_write_variable( file_id, variable_id,                     &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_int32_0d=values_int32_0d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_int32_1d ) )  THEN
                CALL netcdf4_write_variable( file_id, variable_id,                     &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_int32_1d=values_int32_1d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_int32_2d ) )  THEN
                CALL netcdf4_write_variable( file_id, variable_id,                     &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_int32_2d=values_int32_2d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_int32_3d ) )  THEN
                CALL netcdf4_write_variable( file_id, variable_id,                     &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_int32_3d=values_int32_3d_pointer, return_value=output_return_value )
!
!--          working-precision integer output
             ELSEIF ( PRESENT( values_intwp_0d ) )  THEN
                CALL netcdf4_write_variable( file_id, variable_id,                     &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_intwp_0d=values_intwp_0d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_intwp_1d ) )  THEN
                CALL netcdf4_write_variable( file_id, variable_id,                     &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_intwp_1d=values_intwp_1d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_intwp_2d ) )  THEN
                CALL netcdf4_write_variable( file_id, variable_id,                     &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_intwp_2d=values_intwp_2d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_intwp_3d ) )  THEN
                CALL netcdf4_write_variable( file_id, variable_id,                     &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_intwp_3d=values_intwp_3d_pointer, return_value=output_return_value )
!
!--          32bit real output
             ELSEIF ( PRESENT( values_real32_0d ) )  THEN
                CALL netcdf4_write_variable( file_id, variable_id,                     &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_real32_0d=values_real32_0d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_real32_1d ) )  THEN
                CALL netcdf4_write_variable( file_id, variable_id,                     &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_real32_1d=values_real32_1d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_real32_2d ) )  THEN
                CALL netcdf4_write_variable( file_id, variable_id,                     &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_real32_2d=values_real32_2d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_real32_3d ) )  THEN
                CALL netcdf4_write_variable( file_id, variable_id,                     &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_real32_3d=values_real32_3d_pointer, return_value=output_return_value )
!
!--          64bit real output
             ELSEIF ( PRESENT( values_real64_0d ) )  THEN
                CALL netcdf4_write_variable( file_id, variable_id,                     &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_real64_0d=values_real64_0d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_real64_1d ) )  THEN
                CALL netcdf4_write_variable( file_id, variable_id,                     &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_real64_1d=values_real64_1d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_real64_2d ) )  THEN
                CALL netcdf4_write_variable( file_id, variable_id,                     &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_real64_2d=values_real64_2d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_real64_3d ) )  THEN
                CALL netcdf4_write_variable( file_id, variable_id,                     &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_real64_3d=values_real64_3d_pointer, return_value=output_return_value )
!
!--          working-precision real output
             ELSEIF ( PRESENT( values_realwp_0d ) )  THEN
                CALL netcdf4_write_variable( file_id, variable_id,                     &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_realwp_0d=values_realwp_0d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_realwp_1d ) )  THEN
                CALL netcdf4_write_variable( file_id, variable_id,                     &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_realwp_1d=values_realwp_1d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_realwp_2d ) )  THEN
                CALL netcdf4_write_variable( file_id, variable_id,                     &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_realwp_2d=values_realwp_2d_pointer, return_value=output_return_value )
             ELSEIF ( PRESENT( values_realwp_3d ) )  THEN
                CALL netcdf4_write_variable( file_id, variable_id,                     &
                        bounds_start_internal, value_counts, bounds_origin, is_global, &
                        values_realwp_3d=values_realwp_3d_pointer, return_value=output_return_value )
             ELSE
                return_value = 1
                CALL internal_message( 'error', routine_name //                          &
                                       ': output_type not supported by file format "' // &
                                       TRIM( file_format ) // '" ' //                    &
                                       '(variable "' // TRIM( variable_name ) //         &
                                       '", file "' // TRIM( file_name ) // '")!' )
             ENDIF

          CASE DEFAULT
             return_value = 1
             CALL internal_message( 'error', routine_name //                    &
                                    ': file format "' // TRIM( file_format ) // &
                                    '" not supported ' //                       &
                                    '(variable "' // TRIM( variable_name ) //   &
                                    '", file "' // TRIM( file_name ) // '")!' )

       END SELECT

       IF ( return_value == 0  .AND.  output_return_value /= 0 )  THEN
          return_value = 1
          CALL internal_message( 'error', routine_name //                  &
                                 ': error while writing variable ' //      &
                                 '(variable "' // TRIM( variable_name ) // &
                                 '", file "' // TRIM( file_name ) // '")!' )
       ENDIF

    ENDIF

 END FUNCTION dom_write_var

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Finalize output.
!> All necessary steps are carried out to close all output files. If a file could not be closed,
!> this is noted in the error message.
!>
!> @bug if multiple files failed to be closed, only the last failure is given in the error message.
!--------------------------------------------------------------------------------------------------!
 FUNCTION dom_finalize_output() RESULT( return_value )

    CHARACTER(LEN=*), PARAMETER ::  routine_name = 'dom_finalize_output'  !< name of routine

    INTEGER ::  f                      !< loop index
    INTEGER ::  output_return_value    !< return value from called routines
    INTEGER ::  return_value           !< return value
    INTEGER ::  return_value_internal  !< error code after closing a single file


    return_value = 0

    DO  f = 1, nfiles

       IF ( files(f)%is_init )  THEN

          output_return_value = 0
          return_value_internal = 0

          SELECT CASE ( TRIM( files(f)%format ) )

             CASE ( 'binary' )
                CALL binary_finalize( files(f)%id, output_return_value )

             CASE ( 'netcdf4-parallel', 'netcdf4-serial' )
                CALL netcdf4_finalize( files(f)%id, output_return_value )

             CASE DEFAULT
                return_value_internal = 1

          END SELECT

          IF ( output_return_value /= 0 )  THEN
             return_value = output_return_value
             CALL internal_message( 'error', routine_name //             &
                                    ': error while finalizing file "' // &
                                    TRIM( files(f)%name ) // '"' )
          ELSEIF ( return_value_internal /= 0 )  THEN
             return_value = return_value_internal
             CALL internal_message( 'error', routine_name //                     &
                                    ': unsupported file format "' //             &
                                    TRIM( files(f)%format ) // '" for file "' // &
                                    TRIM( files(f)%name ) // '"' )
          ENDIF

       ENDIF

    ENDDO

 END FUNCTION dom_finalize_output

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Return the last created error message.
!--------------------------------------------------------------------------------------------------!
 FUNCTION dom_get_error_message() RESULT( error_message )

    CHARACTER(LEN=800) ::  error_message  !< return error message to main program


    error_message = TRIM( internal_error_message )

    error_message = TRIM( error_message ) // TRIM( binary_get_error_message() )

    error_message = TRIM( error_message ) // TRIM( netcdf4_get_error_message() )

    internal_error_message = ''

 END FUNCTION dom_get_error_message

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Add attribute to database.
!>
!> @todo Try to combine similar code parts and shorten routine.
!--------------------------------------------------------------------------------------------------!
 FUNCTION save_attribute_in_database( file_name, variable_name, attribute, append ) &
             RESULT( return_value )

    CHARACTER(LEN=*), INTENT(IN) ::  file_name      !< name of file
    CHARACTER(LEN=*), INTENT(IN) ::  variable_name  !< name of variable

    CHARACTER(LEN=*), PARAMETER ::  routine_name = 'save_attribute_in_database'  !< name of routine

    INTEGER ::  a             !< loop index
    INTEGER ::  d             !< loop index
    INTEGER ::  f             !< loop index
    INTEGER ::  natts         !< number of attributes
    INTEGER ::  return_value  !< return value

    LOGICAL             ::  found   !< true if variable or dimension of name 'variable_name' found
    LOGICAL, INTENT(IN) ::  append  !< if true, append value to existing value

    TYPE(attribute_type), INTENT(IN) ::  attribute  !< new attribute

    TYPE(attribute_type), DIMENSION(:), ALLOCATABLE ::  atts_tmp  !< temporary attribute list


    return_value = 0
    found = .FALSE.

    CALL internal_message( 'debug', routine_name //                            &
                           ': define attribute "' // TRIM( attribute%name ) // &
                           '" of variable "' // TRIM( variable_name ) //       &
                           '" in file "' // TRIM( file_name ) // '"' )

    DO  f = 1, nfiles

       IF ( TRIM( file_name ) == files(f)%name )  THEN

          IF ( files(f)%is_init )  THEN
             return_value = 1
             CALL internal_message( 'error', routine_name // ': file "' // TRIM( file_name ) // &
                     '" is already initialized. No further attribute definition allowed!' )
             EXIT
          ENDIF
!
!--       Add attribute to file
          IF ( TRIM( variable_name ) == '' )  THEN
!
!--          Initialize first file attribute
             IF ( .NOT. ALLOCATED( files(f)%attributes ) )  THEN
                natts = 1
                ALLOCATE( files(f)%attributes(natts) )
             ELSE
                natts = SIZE( files(f)%attributes )
!
!--             Check if attribute already exists
                DO  a = 1, natts
                   IF ( files(f)%attributes(a)%name == attribute%name )  THEN
                      IF ( append )  THEN
!
!--                      Append existing string attribute
                         files(f)%attributes(a)%value_char =             &
                            TRIM( files(f)%attributes(a)%value_char ) // &
                            TRIM( attribute%value_char )
                      ELSE
                         files(f)%attributes(a) = attribute
                      ENDIF
                      found = .TRUE.
                      EXIT
                   ENDIF
                ENDDO
!
!--             Extend attribute list by 1
                IF ( .NOT. found )  THEN
                   ALLOCATE( atts_tmp(natts) )
                   atts_tmp = files(f)%attributes
                   DEALLOCATE( files(f)%attributes )
                   natts = natts + 1
                   ALLOCATE( files(f)%attributes(natts) )
                   files(f)%attributes(:natts-1) = atts_tmp
                   DEALLOCATE( atts_tmp )
                ENDIF
             ENDIF
!
!--          Save new attribute to the end of the attribute list
             IF ( .NOT. found )  THEN
                files(f)%attributes(natts) = attribute
                found = .TRUE.
             ENDIF

             EXIT

          ELSE
!
!--          Add attribute to dimension
             IF ( ALLOCATED( files(f)%dimensions ) )  THEN

                DO  d = 1, SIZE( files(f)%dimensions )

                   IF ( files(f)%dimensions(d)%name == TRIM( variable_name ) )  THEN

                      IF ( .NOT. ALLOCATED( files(f)%dimensions(d)%attributes ) )  THEN
!
!--                      Initialize first attribute
                         natts = 1
                         ALLOCATE( files(f)%dimensions(d)%attributes(natts) )
                      ELSE
                         natts = SIZE( files(f)%dimensions(d)%attributes )
!
!--                      Check if attribute already exists
                         DO  a = 1, natts
                            IF ( files(f)%dimensions(d)%attributes(a)%name == attribute%name ) &
                            THEN
                               IF ( append )  THEN
!
!--                               Append existing character attribute
                                  files(f)%dimensions(d)%attributes(a)%value_char =             &
                                     TRIM( files(f)%dimensions(d)%attributes(a)%value_char ) // &
                                     TRIM( attribute%value_char )
                               ELSE
!
!--                               Update existing attribute
                                  files(f)%dimensions(d)%attributes(a) = attribute
                               ENDIF
                               found = .TRUE.
                               EXIT
                            ENDIF
                         ENDDO
!
!--                      Extend attribute list
                         IF ( .NOT. found )  THEN
                            ALLOCATE( atts_tmp(natts) )
                            atts_tmp = files(f)%dimensions(d)%attributes
                            DEALLOCATE( files(f)%dimensions(d)%attributes )
                            natts = natts + 1
                            ALLOCATE( files(f)%dimensions(d)%attributes(natts) )
                            files(f)%dimensions(d)%attributes(:natts-1) = atts_tmp
                            DEALLOCATE( atts_tmp )
                         ENDIF
                      ENDIF
!
!--                   Add new attribute to database
                      IF ( .NOT. found )  THEN
                         files(f)%dimensions(d)%attributes(natts) = attribute
                         found = .TRUE.
                      ENDIF

                      EXIT

                   ENDIF  ! dimension found

                ENDDO  ! loop over dimensions

             ENDIF  ! dimensions exist in file
!
!--          Add attribute to variable
             IF ( .NOT. found  .AND.  ALLOCATED( files(f)%variables) )  THEN

                DO  d = 1, SIZE( files(f)%variables )

                   IF ( files(f)%variables(d)%name == TRIM( variable_name ) )  THEN

                      IF ( .NOT. ALLOCATED( files(f)%variables(d)%attributes ) )  THEN
!
!--                      Initialize first attribute
                         natts = 1
                         ALLOCATE( files(f)%variables(d)%attributes(natts) )
                      ELSE
                         natts = SIZE( files(f)%variables(d)%attributes )
!
!--                      Check if attribute already exists
                         DO  a = 1, natts
                            IF ( files(f)%variables(d)%attributes(a)%name == attribute%name )  &
                            THEN
                               IF ( append )  THEN
!
!--                               Append existing character attribute
                                  files(f)%variables(d)%attributes(a)%value_char =             &
                                     TRIM( files(f)%variables(d)%attributes(a)%value_char ) // &
                                     TRIM( attribute%value_char )
                               ELSE
!
!--                               Update existing attribute
                                  files(f)%variables(d)%attributes(a) = attribute
                               ENDIF
                               found = .TRUE.
                               EXIT
                            ENDIF
                         ENDDO
!
!--                      Extend attribute list
                         IF ( .NOT. found )  THEN
                            ALLOCATE( atts_tmp(natts) )
                            atts_tmp = files(f)%variables(d)%attributes
                            DEALLOCATE( files(f)%variables(d)%attributes )
                            natts = natts + 1
                            ALLOCATE( files(f)%variables(d)%attributes(natts) )
                            files(f)%variables(d)%attributes(:natts-1) = atts_tmp
                            DEALLOCATE( atts_tmp )
                         ENDIF

                      ENDIF
!
!--                   Add new attribute to database
                      IF ( .NOT. found )  THEN
                         files(f)%variables(d)%attributes(natts) = attribute
                         found = .TRUE.
                      ENDIF

                      EXIT

                   ENDIF  ! variable found

                ENDDO  ! loop over variables

             ENDIF  ! variables exist in file

             IF ( .NOT. found )  THEN
                return_value = 1
                CALL internal_message( 'error',                                        &
                        routine_name //                                                &
                        ': requested dimension/variable "' // TRIM( variable_name ) // &
                        '" for attribute "' // TRIM( attribute%name ) //               &
                        '" does not exist in file "' // TRIM( file_name ) // '"' )
             ENDIF

             EXIT

          ENDIF  ! variable_name not empty

       ENDIF  ! check file_name

    ENDDO  ! loop over files

    IF ( .NOT. found  .AND.  return_value == 0 )  THEN
       return_value = 1
       CALL internal_message( 'error',                                         &
                              routine_name //                                  &
                              ': requested file "' // TRIM( file_name ) //     &
                              '" for attribute "' // TRIM( attribute%name ) // &
                              '" does not exist' )
    ENDIF

 END FUNCTION save_attribute_in_database

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check database and delete any unused dimensions and empty files (i.e. files
!> without variables).
!--------------------------------------------------------------------------------------------------!
 FUNCTION cleanup_database() RESULT( return_value )

    ! CHARACTER(LEN=*), PARAMETER ::  routine_name = 'cleanup_database'  !< name of routine

    INTEGER ::  d             !< loop index
    INTEGER ::  f             !< loop index
    INTEGER ::  i             !< loop index
    INTEGER ::  ndims         !< number of dimensions in a file
    INTEGER ::  ndims_used    !< number of used dimensions in a file
    INTEGER ::  nfiles_used   !< number of used files
    INTEGER ::  nvars         !< number of variables in a file
    INTEGER ::  return_value  !< return value

    LOGICAL, DIMENSION(1:nfiles)             ::  file_is_used       !< true if file contains variables
    LOGICAL, DIMENSION(:),       ALLOCATABLE ::  dimension_is_used  !< true if dimension is used by any variable

    TYPE(dimension_type), DIMENSION(:), ALLOCATABLE ::  used_dimensions  !< list of used dimensions

    TYPE(file_type), DIMENSION(:), ALLOCATABLE ::  used_files  !< list of used files


    return_value = 0
!
!-- Flag files which contain output variables as used
    file_is_used(:) = .FALSE.
    DO  f = 1, nfiles
       IF ( ALLOCATED( files(f)%variables ) )  THEN
          file_is_used(f) = .TRUE.
       ENDIF
    ENDDO
!
!-- Copy flagged files into temporary list
    nfiles_used = COUNT( file_is_used )
    ALLOCATE( used_files(nfiles_used) )
    i = 0
    DO  f = 1, nfiles
       IF ( file_is_used(f) )  THEN
          i = i + 1
          used_files(i) = files(f)
       ENDIF
    ENDDO
!
!-- Replace file list with list of used files
    DEALLOCATE( files )
    nfiles = nfiles_used
    ALLOCATE( files(nfiles) )
    files = used_files
    DEALLOCATE( used_files )
!
!-- Check every file for unused dimensions
    DO  f = 1, nfiles
!
!--    If a file is already initialized, it was already checked previously
       IF ( files(f)%is_init )  CYCLE
!
!--    Get number of defined dimensions
       ndims = SIZE( files(f)%dimensions )
       ALLOCATE( dimension_is_used(ndims) )
!
!--    Go through all variables and flag all used dimensions
       nvars = SIZE( files(f)%variables )
       DO  d = 1, ndims
          DO  i = 1, nvars
             dimension_is_used(d) = &
                ANY( files(f)%dimensions(d)%name == files(f)%variables(i)%dimension_names )
             IF ( dimension_is_used(d) )  EXIT
          ENDDO
       ENDDO
!
!--    Copy used dimensions to temporary list
       ndims_used = COUNT( dimension_is_used )
       ALLOCATE( used_dimensions(ndims_used) )
       i = 0
       DO  d = 1, ndims
          IF ( dimension_is_used(d) )  THEN
             i = i + 1
             used_dimensions(i) = files(f)%dimensions(d)
          ENDIF
       ENDDO
!
!--    Replace dimension list with list of used dimensions
       DEALLOCATE( files(f)%dimensions )
       ndims = ndims_used
       ALLOCATE( files(f)%dimensions(ndims) )
       files(f)%dimensions = used_dimensions
       DEALLOCATE( used_dimensions )
       DEALLOCATE( dimension_is_used )

    ENDDO

 END FUNCTION cleanup_database

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Open requested output file.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE open_output_file( file_format, file_name, file_id, return_value )

    CHARACTER(LEN=*), INTENT(IN) ::  file_format  !< file format chosen for file
    CHARACTER(LEN=*), INTENT(IN) ::  file_name    !< name of file to be checked

    CHARACTER(LEN=*), PARAMETER ::  routine_name = 'open_output_file'  !< name of routine

    INTEGER, INTENT(OUT) ::  file_id              !< file ID
    INTEGER              ::  output_return_value  !< return value of a called output routine
    INTEGER, INTENT(OUT) ::  return_value         !< return value


    return_value = 0
    output_return_value = 0

    SELECT CASE ( TRIM( file_format ) )

       CASE ( 'binary' )
          CALL binary_open_file( 'binary', file_name, file_id, output_return_value )

       CASE ( 'netcdf4-serial' )
          CALL netcdf4_open_file( 'serial', file_name, file_id, output_return_value )

       CASE ( 'netcdf4-parallel' )
          CALL netcdf4_open_file( 'parallel', file_name, file_id, output_return_value )

       CASE DEFAULT
          return_value = 1

    END SELECT

    IF ( output_return_value /= 0 )  THEN
       return_value = output_return_value
       CALL internal_message( 'error', routine_name // &
                              ': error while opening file "' // TRIM( file_name ) // '"' )
    ELSEIF ( return_value /= 0 )  THEN
       CALL internal_message( 'error', routine_name //                     &
                              ': file "' // TRIM( file_name ) //           &
                              '": file format "' // TRIM( file_format ) // &
                              '" not supported' )
    ENDIF

 END SUBROUTINE open_output_file

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize attributes, dimensions and variables in a file.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE init_file_header( file, return_value )

    ! CHARACTER(LEN=*), PARAMETER ::  routine_name = 'init_file_header'  !< name of routine

    INTEGER              ::  a             !< loop index
    INTEGER              ::  d             !< loop index
    INTEGER, INTENT(OUT) ::  return_value  !< return value

    TYPE(file_type), INTENT(INOUT) ::  file  !< initialize header of this file


    return_value  = 0
!
!-- Write file attributes
    IF ( ALLOCATED( file%attributes ) )  THEN
       DO  a = 1, SIZE( file%attributes )
          return_value = write_attribute( file%format, file%id, file%name,     &
                                          variable_id=no_id, variable_name='', &
                                          attribute=file%attributes(a) )
          IF ( return_value /= 0 )  EXIT
       ENDDO
    ENDIF

    IF ( return_value == 0 )  THEN
!
!--    Initialize file dimensions
       DO  d = 1, SIZE( file%dimensions )

          IF ( .NOT. file%dimensions(d)%is_masked )  THEN
!
!--          Initialize non-masked dimension
             CALL init_file_dimension( file%format, file%id, file%name,       &
                     file%dimensions(d)%id, file%dimensions(d)%name,          &
                     file%dimensions(d)%data_type, file%dimensions(d)%length, &
                     file%dimensions(d)%variable_id, return_value )

          ELSE
!
!--          Initialize masked dimension
             CALL init_file_dimension( file%format, file%id, file%name,            &
                     file%dimensions(d)%id, file%dimensions(d)%name,               &
                     file%dimensions(d)%data_type, file%dimensions(d)%length_mask, &
                     file%dimensions(d)%variable_id, return_value )

          ENDIF

          IF ( return_value == 0  .AND.  ALLOCATED( file%dimensions(d)%attributes ) )  THEN
!
!--          Write dimension attributes
             DO  a = 1, SIZE( file%dimensions(d)%attributes )
                return_value = write_attribute( file%format, file%id, file%name, &
                                  variable_id=file%dimensions(d)%variable_id,    &
                                  variable_name=file%dimensions(d)%name,         &
                                  attribute=file%dimensions(d)%attributes(a) )
                IF ( return_value /= 0 )  EXIT
             ENDDO
          ENDIF

          IF ( return_value /= 0 )  EXIT

       ENDDO
!
!--    Save dimension IDs for variables wihtin database
       IF ( return_value == 0 )  &
          CALL collect_dimesion_ids_for_variables( file%name, file%variables, file%dimensions, &
                                                   return_value )
!
!--    Initialize file variables
       IF ( return_value == 0 )  THEN
          DO  d = 1, SIZE( file%variables )

             CALL init_file_variable( file%format, file%id, file%name,                          &
                     file%variables(d)%id, file%variables(d)%name, file%variables(d)%data_type, &
                     file%variables(d)%dimension_ids,                                           &
                     file%variables(d)%is_global, return_value )

             IF ( return_value == 0  .AND.  ALLOCATED( file%variables(d)%attributes ) )  THEN
!
!--             Write variable attributes
                DO  a = 1, SIZE( file%variables(d)%attributes )
                   return_value = write_attribute( file%format, file%id, file%name, &
                                     variable_id=file%variables(d)%id,              &
                                     variable_name=file%variables(d)%name,          &
                                     attribute=file%variables(d)%attributes(a) )
                   IF ( return_value /= 0 )  EXIT
                ENDDO
             ENDIF

             IF ( return_value /= 0 )  EXIT

          ENDDO
       ENDIF

    ENDIF

 END SUBROUTINE init_file_header

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize dimension in file.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE init_file_dimension( file_format, file_id, file_name,              &
               dimension_id, dimension_name, dimension_type, dimension_length, &
               variable_id, return_value )

    CHARACTER(LEN=*), INTENT(IN) ::  dimension_name  !< name of dimension
    CHARACTER(LEN=*), INTENT(IN) ::  dimension_type  !< data type of dimension
    CHARACTER(LEN=*), INTENT(IN) ::  file_format     !< file format chosen for file
    CHARACTER(LEN=*), INTENT(IN) ::  file_name       !< name of file

    CHARACTER(LEN=*), PARAMETER ::  routine_name = 'init_file_dimension'  !< file format chosen for file

    INTEGER, INTENT(OUT) ::  dimension_id         !< dimension ID
    INTEGER, INTENT(IN)  ::  dimension_length     !< length of dimension
    INTEGER, INTENT(IN)  ::  file_id              !< file ID
    INTEGER              ::  output_return_value  !< return value of a called output routine
    INTEGER, INTENT(OUT) ::  return_value         !< return value
    INTEGER, INTENT(OUT) ::  variable_id          !< associated variable ID


    return_value = 0
    output_return_value = 0

    temp_string = '(file "' // TRIM( file_name ) // &
                  '", dimension "' // TRIM( dimension_name ) // '")'

    SELECT CASE ( TRIM( file_format ) )

       CASE ( 'binary' )
          CALL binary_init_dimension( 'binary', file_id, dimension_id, variable_id, &
                  dimension_name, dimension_type, dimension_length,                 &
                  return_value=output_return_value )

       CASE ( 'netcdf4-serial' )
          CALL netcdf4_init_dimension( 'serial', file_id, dimension_id, variable_id, &
                  dimension_name, dimension_type, dimension_length,                  &
                  return_value=output_return_value )

       CASE ( 'netcdf4-parallel' )
          CALL netcdf4_init_dimension( 'parallel', file_id, dimension_id, variable_id, &
                  dimension_name, dimension_type, dimension_length,                    &
                  return_value=output_return_value )

       CASE DEFAULT
          return_value = 1
          CALL internal_message( 'error', routine_name //                    &
                                 ': file format "' // TRIM( file_format ) // &
                                 '" not supported ' // TRIM( temp_string ) )

    END SELECT

    IF ( output_return_value /= 0 )  THEN
       return_value = output_return_value
       CALL internal_message( 'error', routine_name // &
                              ': error while defining dimension ' // TRIM( temp_string ) )
    ENDIF

 END SUBROUTINE init_file_dimension

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize variable.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE init_file_variable( file_format, file_id, file_name,        &
                                variable_id, variable_name, variable_type, dimension_ids, &
                                is_global, return_value )

    CHARACTER(LEN=*), INTENT(IN) ::  file_format    !< file format chosen for file
    CHARACTER(LEN=*), INTENT(IN) ::  file_name      !< file name
    CHARACTER(LEN=*), INTENT(IN) ::  variable_name  !< name of variable
    CHARACTER(LEN=*), INTENT(IN) ::  variable_type  !< data type of variable

    CHARACTER(LEN=*), PARAMETER ::  routine_name = 'init_file_variable'  !< file format chosen for file

    INTEGER, INTENT(IN)  ::  file_id              !< file ID
    INTEGER              ::  output_return_value  !< return value of a called output routine
    INTEGER, INTENT(OUT) ::  return_value         !< return value
    INTEGER, INTENT(OUT) ::  variable_id          !< variable ID

    INTEGER, DIMENSION(:), INTENT(IN) ::  dimension_ids  !< list of dimension IDs used by variable

    LOGICAL, INTENT(IN)  ::  is_global  !< true if variable is global


    return_value = 0
    output_return_value = 0

    temp_string = '(file "' // TRIM( file_name ) // &
                  '", variable "' // TRIM( variable_name ) // '")'

    SELECT CASE ( TRIM( file_format ) )

       CASE ( 'binary' )
          CALL binary_init_variable( 'binary', file_id, variable_id, variable_name, &
                  variable_type, dimension_ids, is_global, return_value=output_return_value )

       CASE ( 'netcdf4-serial' )
          CALL netcdf4_init_variable( 'serial', file_id, variable_id, variable_name, &
                  variable_type, dimension_ids, is_global, return_value=output_return_value )

       CASE ( 'netcdf4-parallel' )
          CALL netcdf4_init_variable( 'parallel', file_id, variable_id, variable_name, &
                  variable_type, dimension_ids, is_global, return_value=output_return_value )

       CASE DEFAULT
          return_value = 1
          CALL internal_message( 'error', routine_name //                    &
                                 ': file format "' // TRIM( file_format ) // &
                                 '" not supported ' // TRIM( temp_string ) )

    END SELECT

    IF ( output_return_value /= 0 )  THEN
       return_value = output_return_value
       CALL internal_message( 'error', routine_name // &
                              ': error while defining variable ' // TRIM( temp_string ) )
    ENDIF

 END SUBROUTINE init_file_variable

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write attribute to file.
!--------------------------------------------------------------------------------------------------!
 FUNCTION write_attribute( file_format, file_id, file_name,        &
                           variable_id, variable_name, attribute ) RESULT( return_value )

    CHARACTER(LEN=*), INTENT(IN) ::  file_format    !< file format chosen for file
    CHARACTER(LEN=*), INTENT(IN) ::  file_name      !< file name
    CHARACTER(LEN=*), INTENT(IN) ::  variable_name  !< variable name

    CHARACTER(LEN=*), PARAMETER ::  routine_name = 'write_attribute'  !< file format chosen for file

    INTEGER, INTENT(IN) ::  file_id              !< file ID
    INTEGER             ::  return_value         !< return value
    INTEGER             ::  output_return_value  !< return value of a called output routine
    INTEGER, INTENT(IN) ::  variable_id          !< variable ID

    TYPE(attribute_type), INTENT(IN) ::  attribute  !< attribute to be written


    return_value = 0
    output_return_value = 0
!
!-- Prepare for possible error message
    temp_string = '(file "' // TRIM( file_name ) //           &
                  '", variable "' // TRIM( variable_name ) // &
                  '", attribute "' // TRIM( attribute%name ) // '")'
!
!-- Write attribute to file
    SELECT CASE ( TRIM( file_format ) )

       CASE ( 'binary' )

          SELECT CASE ( TRIM( attribute%data_type ) )

             CASE( 'char' )
                CALL binary_write_attribute( file_id=file_id, variable_id=variable_id,  &
                        attribute_name=attribute%name, value_char=attribute%value_char, &
                        return_value=output_return_value )

             CASE( 'int8' )
                CALL binary_write_attribute( file_id=file_id, variable_id=variable_id,  &
                        attribute_name=attribute%name, value_int8=attribute%value_int8, &
                        return_value=output_return_value )

             CASE( 'int16' )
                CALL binary_write_attribute( file_id=file_id, variable_id=variable_id,    &
                        attribute_name=attribute%name, value_int16=attribute%value_int16, &
                        return_value=output_return_value )

             CASE( 'int32' )
                CALL binary_write_attribute( file_id=file_id, variable_id=variable_id,    &
                        attribute_name=attribute%name, value_int32=attribute%value_int32, &
                        return_value=output_return_value )

             CASE( 'real32' )
                CALL binary_write_attribute( file_id=file_id, variable_id=variable_id,      &
                        attribute_name=attribute%name, value_real32=attribute%value_real32, &
                        return_value=output_return_value )

             CASE( 'real64' )
                CALL binary_write_attribute( file_id=file_id, variable_id=variable_id,      &
                        attribute_name=attribute%name, value_real64=attribute%value_real64, &
                        return_value=output_return_value )

             CASE DEFAULT
                return_value = 1
                CALL internal_message( 'error', routine_name //                     &
                                       ': file format "' // TRIM( file_format ) //  &
                                       '" does not support attribute data type "'// &
                                       TRIM( attribute%data_type ) //               &
                                       '" ' // TRIM( temp_string ) )

          END SELECT

       CASE ( 'netcdf4-parallel', 'netcdf4-serial' )

          SELECT CASE ( TRIM( attribute%data_type ) )

             CASE( 'char' )
                CALL netcdf4_write_attribute( file_id=file_id, variable_id=variable_id, &
                        attribute_name=attribute%name, value_char=attribute%value_char, &
                        return_value=output_return_value )

             CASE( 'int8' )
                CALL netcdf4_write_attribute( file_id=file_id, variable_id=variable_id, &
                        attribute_name=attribute%name, value_int8=attribute%value_int8, &
                        return_value=output_return_value )

             CASE( 'int16' )
                CALL netcdf4_write_attribute( file_id=file_id, variable_id=variable_id,   &
                        attribute_name=attribute%name, value_int16=attribute%value_int16, &
                        return_value=output_return_value )

             CASE( 'int32' )
                CALL netcdf4_write_attribute( file_id=file_id, variable_id=variable_id,   &
                        attribute_name=attribute%name, value_int32=attribute%value_int32, &
                        return_value=output_return_value )

             CASE( 'real32' )
                CALL netcdf4_write_attribute( file_id=file_id, variable_id=variable_id,     &
                        attribute_name=attribute%name, value_real32=attribute%value_real32, &
                        return_value=output_return_value )

             CASE( 'real64' )
                CALL netcdf4_write_attribute( file_id=file_id, variable_id=variable_id,     &
                        attribute_name=attribute%name, value_real64=attribute%value_real64, &
                        return_value=output_return_value )

             CASE DEFAULT
                return_value = 1
                CALL internal_message( 'error', routine_name //                     &
                                       ': file format "' // TRIM( file_format ) //  &
                                       '" does not support attribute data type "'// &
                                       TRIM( attribute%data_type ) //               &
                                       '" ' // TRIM( temp_string ) )

          END SELECT

       CASE DEFAULT
          return_value = 1
          CALL internal_message( 'error', routine_name //                                &
                                 ': unsupported file format "' // TRIM( file_format ) // &
                                 '" ' // TRIM( temp_string ) )

    END SELECT

    IF ( output_return_value /= 0 )  THEN
       return_value = output_return_value
       CALL internal_message( 'error', routine_name // &
                              ': error while writing attribute ' // TRIM( temp_string ) )
    ENDIF

 END FUNCTION write_attribute

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Get dimension IDs and save them to variables.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE collect_dimesion_ids_for_variables( file_name, variables, dimensions, return_value )

    CHARACTER(LEN=*), INTENT(IN) ::  file_name !< name of file

    CHARACTER(LEN=*), PARAMETER ::  routine_name = 'collect_dimesion_ids_for_variables'  !< file format chosen for file

    INTEGER              ::  d             !< loop index
    INTEGER              ::  i             !< loop index
    INTEGER              ::  j             !< loop index
    INTEGER              ::  ndims         !< number of dimensions
    INTEGER              ::  nvars         !< number of variables
    INTEGER, INTENT(OUT) ::  return_value  !< return value

    LOGICAL ::  found  !< true if dimension required by variable was found in dimension list

    TYPE(dimension_type), DIMENSION(:), INTENT(IN) ::  dimensions  !< list of dimensions in file

    TYPE(variable_type), DIMENSION(:), INTENT(INOUT) ::  variables  !< list of variables in file


    return_value  = 0
    ndims = SIZE( dimensions )
    nvars = SIZE( variables )

    DO  i = 1, nvars
       DO  j = 1, SIZE( variables(i)%dimension_names )
          found = .FALSE.
          DO  d = 1, ndims
             IF ( variables(i)%dimension_names(j) == dimensions(d)%name )  THEN
                variables(i)%dimension_ids(j) = dimensions(d)%id
                found = .TRUE.
                EXIT
             ENDIF
          ENDDO
          IF ( .NOT. found )  THEN
             return_value = 1
             CALL internal_message( 'error', routine_name //                                &
                     ': required dimension "' // TRIM( variables(i)%dimension_names(j) ) // &
                     '" is undefined (variable "' // TRIM( variables(i)%name ) //           &
                     '", file "' // TRIM( file_name ) // '")!' )
             EXIT
          ENDIF
       ENDDO
       IF ( .NOT. found )  EXIT
    ENDDO

 END SUBROUTINE collect_dimesion_ids_for_variables

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Leave file definition/initialization.
!>
!> @todo Do we need an MPI barrier at the end?
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE stop_file_header_definition( file_format, file_id, file_name, return_value )

    CHARACTER(LEN=*), INTENT(IN) ::  file_format  !< file format
    CHARACTER(LEN=*), INTENT(IN) ::  file_name    !< file name

    CHARACTER(LEN=*), PARAMETER ::  routine_name = 'stop_file_header_definition'  !< name of routine

    INTEGER, INTENT(IN)  ::  file_id              !< file id
    INTEGER              ::  output_return_value  !< return value of a called output routine
    INTEGER, INTENT(OUT) ::  return_value         !< return value


    return_value = 0
    output_return_value = 0

    temp_string = '(file "' // TRIM( file_name ) // '")'

    SELECT CASE ( TRIM( file_format ) )

       CASE ( 'binary' )
          CALL binary_stop_file_header_definition( file_id, output_return_value )

       CASE ( 'netcdf4-parallel', 'netcdf4-serial' )
          CALL netcdf4_stop_file_header_definition( file_id, output_return_value )

       CASE DEFAULT
          return_value = 1
          CALL internal_message( 'error', routine_name //                    &
                                 ': file format "' // TRIM( file_format ) // &
                                 '" not supported ' // TRIM( temp_string ) )

    END SELECT

    IF ( output_return_value /= 0 )  THEN
       return_value = output_return_value
       CALL internal_message( 'error', routine_name //                          &
                              ': error while leaving file-definition state ' // &
                              TRIM( temp_string ) )
    ENDIF

    ! CALL MPI_Barrier( MPI_COMM_WORLD, return_value )

 END SUBROUTINE stop_file_header_definition

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Find a requested variable 'variable_name' and its used dimensions in requested file 'file_name'.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE find_var_in_file( file_name, variable_name, file_format, file_id, variable_id, &
                              is_global, dimensions, return_value )

    CHARACTER(LEN=charlen), INTENT(OUT) ::  file_format    !< file format chosen for file
    CHARACTER(LEN=*),       INTENT(IN)  ::  file_name      !< name of file
    CHARACTER(LEN=*),       INTENT(IN)  ::  variable_name  !< name of variable

    CHARACTER(LEN=*), PARAMETER ::  routine_name = 'find_var_in_file'  !< name of routine

    INTEGER              ::  d             !< loop index
    INTEGER              ::  dd            !< loop index
    INTEGER              ::  f             !< loop index
    INTEGER, INTENT(OUT) ::  file_id       !< file ID
    INTEGER, INTENT(OUT) ::  return_value  !< return value
    INTEGER, INTENT(OUT) ::  variable_id   !< variable ID

    INTEGER, DIMENSION(:), ALLOCATABLE ::  dimension_ids  !< list of dimension IDs used by variable

    LOGICAL              ::  found      !< true if requested variable found in requested file
    LOGICAL, INTENT(OUT) ::  is_global  !< true if variable is global

    TYPE(dimension_type), DIMENSION(:), ALLOCATABLE, INTENT(OUT) ::  dimensions  !< list of dimensions used by variable


    return_value = 0
    found = .FALSE.

    DO  f = 1, nfiles
       IF ( TRIM( file_name ) == TRIM( files(f)%name ) )  THEN

          IF ( .NOT. files(f)%is_init )  THEN
             return_value = 1
             CALL internal_message( 'error', routine_name //                     &
                                    ': file not initialized. ' //                &
                                    'Writing variable to file is impossible ' // &
                                    '(variable "' // TRIM( variable_name ) //    &
                                    '", file "' // TRIM( file_name ) // '")!' )
             EXIT
          ENDIF

          file_id     = files(f)%id
          file_format = files(f)%format
!
!--       Search for variable in file
          DO  d = 1, SIZE( files(f)%variables )
             IF ( TRIM( variable_name ) == TRIM( files(f)%variables(d)%name ) )  THEN

                variable_id    = files(f)%variables(d)%id
                is_global = files(f)%variables(d)%is_global

                ALLOCATE( dimension_ids(SIZE( files(f)%variables(d)%dimension_ids )) )
                ALLOCATE( dimensions(SIZE( files(f)%variables(d)%dimension_ids )) )

                dimension_ids = files(f)%variables(d)%dimension_ids

                found = .TRUE.
                EXIT

             ENDIF
          ENDDO

          IF ( found )  THEN
!
!--          Get list of dimensions used by variable
             DO  d = 1, SIZE( files(f)%dimensions )
                DO  dd = 1, SIZE( dimension_ids )
                   IF ( dimension_ids(dd) == files(f)%dimensions(d)%id )  THEN
                      dimensions(dd) = files(f)%dimensions(d)
                      EXIT
                   ENDIF
                ENDDO
             ENDDO

          ELSE
!
!--          If variable was not found, search for a dimension instead
             DO  d = 1, SIZE( files(f)%dimensions )
                IF ( TRIM( variable_name ) == TRIM( files(f)%dimensions(d)%name ) )  THEN

                   variable_id    = files(f)%dimensions(d)%variable_id
                   is_global = .TRUE.

                   ALLOCATE( dimensions(1) )

                   dimensions(1) = files(f)%dimensions(d)

                   found = .TRUE.
                   EXIT

                ENDIF
             ENDDO

          ENDIF
!
!--       If variable was not found in requested file, return an error
          IF ( .NOT. found )  THEN
             return_value = 1
             CALL internal_message( 'error', routine_name //                  &
                                    ': variable not found in file ' //        &
                                    '(variable "' // TRIM( variable_name ) // &
                                    '", file "' // TRIM( file_name ) // '")!' )
          ENDIF

          EXIT

       ENDIF  ! file found
    ENDDO  ! loop over files

    IF ( .NOT. found  .AND.  return_value == 0 )  THEN
       return_value = 1
       CALL internal_message( 'error', routine_name //                  &
                              ': file not found ' //                    &
                              '(variable "' // TRIM( variable_name ) // &
                              '", file "' // TRIM( file_name ) // '")!' )
    ENDIF

 END SUBROUTINE find_var_in_file

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Search for masked indices of dimensions within the given bounds ('bounds_start' and
!> 'bounds_end'). Return the masked indices ('masked_indices') of the dimensions, the first index
!> of the masked dimensions containing these indices ('bounds_masked_start'), the count of masked
!> indices within given bounds ('value_counts') and the origin index of each dimension
!> ('bounds_origin'). If, for any dimension, no masked index lies within the given bounds, counts,
!> starts and origins are set to zero for all dimensions.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE get_masked_indices_and_masked_dimension_bounds(                             &
               dimensions, bounds_start, bounds_end, bounds_masked_start, value_counts, &
               bounds_origin, masked_indices )

    ! CHARACTER(LEN=*), PARAMETER ::  routine_name = 'get_masked_indices_and_masked_dimension_bounds'  !< name of routine

    INTEGER ::  d  !< loop index
    INTEGER ::  i  !< loop index

    INTEGER, DIMENSION(:), INTENT(IN)  ::  bounds_end           !< upper bonuds to be searched in
    INTEGER, DIMENSION(:), INTENT(OUT) ::  bounds_masked_start  !< lower bounds of masked dimensions within given bounds
    INTEGER, DIMENSION(:), INTENT(OUT) ::  bounds_origin        !< first index of each dimension, 0 if dimension is masked
    INTEGER, DIMENSION(:), INTENT(IN)  ::  bounds_start         !< lower bounds to be searched in
    INTEGER, DIMENSION(:), INTENT(OUT) ::  value_counts         !< count of indices per dimension to be output

    INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) ::  masked_indices  !< masked indices within given bounds

    TYPE(dimension_type), DIMENSION(:), INTENT(IN) ::  dimensions  !< dimensions to be searched for masked indices


    ALLOCATE( masked_indices(SIZE( dimensions ),0:MAXVAL( bounds_end - bounds_start + 1 )) )
    masked_indices = -HUGE( 0 )
!
!-- Check for masking and update lower and upper bounds if masked
    DO  d = 1, SIZE( dimensions )

       IF ( dimensions(d)%is_masked )  THEN

          bounds_origin(d) = 0

          bounds_masked_start(d) = -HUGE( 0 )
!
!--       Find number of masked values within given variable bounds
          value_counts(d) = 0
          DO  i = LBOUND( dimensions(d)%masked_indices, DIM=1 ), &
                  UBOUND( dimensions(d)%masked_indices, DIM=1 )
!
!--          Is masked index within given bounds?
             IF ( dimensions(d)%masked_indices(i) >= bounds_start(d)  .AND.  &
                  dimensions(d)%masked_indices(i) <= bounds_end(d)           )  THEN
!
!--             Save masked index
                masked_indices(d,value_counts(d)) = dimensions(d)%masked_indices(i)
                value_counts(d) = value_counts(d) + 1
!
!--             Save bounds of mask within given bounds
                IF ( bounds_masked_start(d) == -HUGE( 0 ) )  bounds_masked_start(d) = i

             ENDIF

          ENDDO
!
!--       Set masked bounds to zero if no masked index lies within bounds
          IF ( value_counts(d) == 0 )  THEN
             bounds_origin(:) = 0
             bounds_masked_start(:) = 0
             value_counts(:) = 0
             EXIT
          ENDIF

       ELSE
!
!--       If dimension is not masked, save all indices within bounds for output
          bounds_origin(d) = dimensions(d)%bounds(1)
          bounds_masked_start(d) = bounds_start(d)
          value_counts(d) = bounds_end(d) - bounds_start(d) + 1

          DO  i = 0, value_counts(d) - 1
             masked_indices(d,i) = bounds_start(d) + i
          ENDDO

       ENDIF

    ENDDO

 END SUBROUTINE get_masked_indices_and_masked_dimension_bounds

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

       WRITE( internal_error_message, '(A,A)' ) 'DOM ERROR: ', string

    ELSEIF ( TRIM( level ) == 'debug'  .AND.  print_debug_output )  THEN

       WRITE( debug_output_unit, '(A,A)' ) 'DOM DEBUG: ', string
       FLUSH( debug_output_unit )

    ENDIF

 END SUBROUTINE internal_message

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Print contents of the created database to debug_output_unit. This routine can be called at any
!> stage after the call to 'dom_init'. Multiple calls are possible.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dom_database_debug_output

    CHARACTER(LEN=*), PARAMETER ::  separation_string = '---'                   !< string separating blocks in output
    CHARACTER(LEN=50)           ::  write_format1                               !< format for write statements
    CHARACTER(LEN=*), PARAMETER ::  routine_name = 'dom_database_debug_output'  !< name of this routine

    INTEGER            ::  f                       !< loop index
    INTEGER, PARAMETER ::  indent_depth = 3        !< space per indentation
    INTEGER            ::  indent_level            !< indentation level
    INTEGER, PARAMETER ::  max_keyname_length = 6  !< length of longest key name
    INTEGER            ::  natts                   !< number of attributes
    INTEGER            ::  ndims                   !< number of dimensions
    INTEGER            ::  nvars                   !< number of variables


    CALL internal_message( 'debug', routine_name // ': write database to debug output' )

    WRITE( debug_output_unit, '(A)' ) 'DOM database:'
    WRITE( debug_output_unit, '(A)' ) REPEAT( separation_string // ' ', 4 )

    IF ( .NOT. ALLOCATED( files ) .OR. nfiles == 0 )  THEN

       WRITE( debug_output_unit, '(A)' ) 'database is empty'

    ELSE

       indent_level = 1
       WRITE( write_format1, '(A,I3,A,I3,A)' ) '(T', indent_level * indent_depth + 1, ',A,T',  &
                                         indent_level * indent_depth + 1 + max_keyname_length, &
                                         ',(": ")'

       DO  f = 1, nfiles

          natts = 0
          ndims = 0
          nvars = 0
          IF ( ALLOCATED( files(f)%attributes ) ) natts = SIZE( files(f)%attributes )
          IF ( ALLOCATED( files(f)%dimensions ) ) ndims = SIZE( files(f)%dimensions )
          IF ( ALLOCATED( files(f)%variables  ) ) nvars = SIZE( files(f)%variables  )

          WRITE( debug_output_unit, '(A)' ) 'file:'
          WRITE( debug_output_unit, TRIM( write_format1 ) // ',A)' ) 'name', TRIM( files(f)%name )
          WRITE( debug_output_unit, TRIM( write_format1 ) // ',A)' ) 'format', TRIM(files(f)%format)
          WRITE( debug_output_unit, TRIM( write_format1 ) // ',I7)' ) 'id', files(f)%id
          WRITE( debug_output_unit, TRIM( write_format1 ) // ',L1)' ) 'is init', files(f)%is_init
          WRITE( debug_output_unit, TRIM( write_format1 ) // ',I7)' ) '#atts', natts
          WRITE( debug_output_unit, TRIM( write_format1 ) // ',I7)' ) '#dims', ndims
          WRITE( debug_output_unit, TRIM( write_format1 ) // ',I7)' ) '#vars', nvars

          IF ( natts /= 0 )  CALL print_attributes( indent_level, files(f)%attributes )
          IF ( ndims /= 0 )  CALL print_dimensions( indent_level, files(f)%dimensions )
          IF ( nvars /= 0 )  CALL print_variables( indent_level, files(f)%variables )

          WRITE( debug_output_unit, '(/A/)' ) REPEAT( separation_string // ' ', 4 )

       ENDDO

    ENDIF

    CONTAINS

!--------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Print list of attributes.
!--------------------------------------------------------------------------------------------!
    SUBROUTINE print_attributes( indent_level, attributes )

       CHARACTER(LEN=50) ::  write_format1  !< format for write statements
       CHARACTER(LEN=50) ::  write_format2  !< format for write statements

       INTEGER             ::  i                       !< loop index
       INTEGER, INTENT(IN) ::  indent_level            !< indentation level
       INTEGER, PARAMETER  ::  max_keyname_length = 6  !< length of longest key name
       INTEGER             ::  nelement                !< number of elements to print

       TYPE(attribute_type), DIMENSION(:), INTENT(IN) ::  attributes  !< list of attributes


       WRITE( write_format1, '(A,I3,A)' ) '(T', indent_level * indent_depth + 1, ',A'
       WRITE( write_format2, '(A,I3,A,I3,A)' ) &
          '(T', ( indent_level + 1 ) * indent_depth + 1, ',A,T', &
          ( indent_level + 1 ) * indent_depth + 1 + max_keyname_length, ',(": ")'

       WRITE( debug_output_unit, TRIM( write_format1 ) // ',A)' ) &
          REPEAT( separation_string // ' ', 4 )
       WRITE( debug_output_unit, TRIM( write_format1 ) // ')' ) 'attributes:'

       nelement = SIZE( attributes )
       DO  i = 1, nelement
          WRITE( debug_output_unit, TRIM( write_format2 ) // ',A)' ) &
             'name', TRIM( attributes(i)%name )
          WRITE( debug_output_unit, TRIM( write_format2 ) // ',A)' ) &
             'type', TRIM( attributes(i)%data_type )

          IF ( TRIM( attributes(i)%data_type ) == 'char' )  THEN
             WRITE( debug_output_unit, TRIM( write_format2 ) // ',A)' ) &
                'value', TRIM( attributes(i)%value_char )
          ELSEIF ( TRIM( attributes(i)%data_type ) == 'int8' )  THEN
             WRITE( debug_output_unit, TRIM( write_format2 ) // ',I4)' ) &
                'value', attributes(i)%value_int8
          ELSEIF ( TRIM( attributes(i)%data_type ) == 'int16' )  THEN
             WRITE( debug_output_unit, TRIM( write_format2 ) // ',I6)' ) &
                'value', attributes(i)%value_int16
          ELSEIF ( TRIM( attributes(i)%data_type ) == 'int32' )  THEN
             WRITE( debug_output_unit, TRIM( write_format2 ) // ',I11)' ) &
                'value', attributes(i)%value_int32
          ELSEIF ( TRIM( attributes(i)%data_type ) == 'real32' )  THEN
             WRITE( debug_output_unit, TRIM( write_format2 ) // ',E14.7)' ) &
                'value', attributes(i)%value_real32
          ELSEIF (  TRIM(attributes(i)%data_type) == 'real64' )  THEN
             WRITE( debug_output_unit, TRIM( write_format2 ) // ',E22.15)' ) &
                'value', attributes(i)%value_real64
          ENDIF
          IF ( i < nelement )  &
             WRITE( debug_output_unit, TRIM( write_format1 ) // ')' ) separation_string
       ENDDO

    END SUBROUTINE print_attributes

!--------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Print list of dimensions.
!--------------------------------------------------------------------------------------------!
    SUBROUTINE print_dimensions( indent_level, dimensions )

       CHARACTER(LEN=50) ::  write_format1  !< format for write statements
       CHARACTER(LEN=50) ::  write_format2  !< format for write statements

       INTEGER             ::  i                        !< loop index
       INTEGER, INTENT(IN) ::  indent_level             !< indentation level
       INTEGER             ::  j                        !< loop index
       INTEGER, PARAMETER  ::  max_keyname_length = 15  !< length of longest key name
       INTEGER             ::  nelement                 !< number of elements to print

       LOGICAL ::  is_masked  !< true if dimension is masked

       TYPE(dimension_type), DIMENSION(:), INTENT(IN) ::  dimensions  !< list of dimensions


       WRITE( write_format1, '(A,I3,A)' ) '(T', indent_level * indent_depth + 1, ',A'
       WRITE( write_format2, '(A,I3,A,I3,A)' ) &
          '(T', ( indent_level + 1 ) * indent_depth + 1, ',A,T', &
          ( indent_level + 1 ) * indent_depth + 1 + max_keyname_length, ',(": ")'

       WRITE( debug_output_unit, TRIM( write_format1 ) // ',A)' ) &
          REPEAT( separation_string // ' ', 4 )
       WRITE( debug_output_unit, TRIM( write_format1 ) // ')' ) 'dimensions:'

       nelement = SIZE( dimensions )
       DO  i = 1, nelement
          is_masked = dimensions(i)%is_masked
!
!--       Print general information
          WRITE( debug_output_unit, TRIM( write_format2 ) // ',A)' ) &
             'name', TRIM( dimensions(i)%name )
          WRITE( debug_output_unit, TRIM( write_format2 ) // ',A)' ) &
             'type', TRIM( dimensions(i)%data_type )
          WRITE( debug_output_unit, TRIM( write_format2 ) // ',I7)' ) &
             'id', dimensions(i)%id
          WRITE( debug_output_unit, TRIM( write_format2 ) // ',I7)' ) &
             'length', dimensions(i)%length
          WRITE( debug_output_unit, TRIM( write_format2 ) // ',I7,A,I7)' ) &
             'bounds', dimensions(i)%bounds(1), ',', dimensions(i)%bounds(2)
          WRITE( debug_output_unit, TRIM( write_format2 ) // ',L1)' ) &
             'is masked', dimensions(i)%is_masked
!
!--       Print information about mask
          IF ( is_masked )  THEN
             WRITE( debug_output_unit, TRIM( write_format2 ) // ',I7)' ) &
                'masked length', dimensions(i)%length_mask

             WRITE( debug_output_unit, TRIM( write_format2 ) // ',L1)', ADVANCE='no' ) &
                'mask', dimensions(i)%mask(dimensions(i)%bounds(1))
             DO  j = dimensions(i)%bounds(1)+1, dimensions(i)%bounds(2)
                WRITE( debug_output_unit, '(A,L1)', ADVANCE='no' ) ',', dimensions(i)%mask(j)
             ENDDO
             WRITE( debug_output_unit, '(A)' )  ''  ! write line-end

             WRITE( debug_output_unit, TRIM( write_format2 ) // ',I6)', ADVANCE='no' ) &
                'masked indices', dimensions(i)%masked_indices(0)
             DO  j = 1, dimensions(i)%length_mask-1
                WRITE( debug_output_unit, '(A,I6)', ADVANCE='no' ) &
                   ',', dimensions(i)%masked_indices(j)
             ENDDO
             WRITE( debug_output_unit, '(A)' )  ''  ! write line-end
          ENDIF
!
!--       Print saved values
          IF ( ALLOCATED( dimensions(i)%values_int8 ) )  THEN

             WRITE( debug_output_unit, TRIM( write_format2 ) // ',I4)', ADVANCE='no' ) &
                'values', dimensions(i)%values_int8(dimensions(i)%bounds(1))
             DO  j = dimensions(i)%bounds(1)+1, dimensions(i)%bounds(2)
                WRITE( debug_output_unit, '(A,I4)', ADVANCE='no' ) &
                   ',', dimensions(i)%values_int8(j)
             ENDDO
             WRITE( debug_output_unit, '(A)' )  ''  ! write line-end
             IF ( is_masked )  THEN
                WRITE( debug_output_unit, TRIM( write_format2 ) // ',I4)', ADVANCE='no' ) &
                   'masked values', dimensions(i)%masked_values_int8(0)
                DO  j = 1, dimensions(i)%length_mask-1
                   WRITE( debug_output_unit, '(A,I4)', ADVANCE='no' ) &
                      ',', dimensions(i)%masked_values_int8(j)
                ENDDO
                WRITE( debug_output_unit, '(A)' )  ''  ! write line-end
             ENDIF

          ELSEIF ( ALLOCATED( dimensions(i)%values_int16 ) )  THEN

             WRITE( debug_output_unit, TRIM( write_format2 ) // ',I6)', ADVANCE='no' ) &
                'values', dimensions(i)%values_int16(dimensions(i)%bounds(1))
             DO  j = dimensions(i)%bounds(1)+1, dimensions(i)%bounds(2)
                WRITE( debug_output_unit, '(A,I6)', ADVANCE='no' ) &
                   ',', dimensions(i)%values_int16(j)
             ENDDO
             WRITE( debug_output_unit, '(A)' )  ''  ! write line-end
             IF ( is_masked )  THEN
                WRITE( debug_output_unit, TRIM( write_format2 ) // ',I6)', ADVANCE='no' ) &
                   'masked values', dimensions(i)%masked_values_int16(0)
                DO  j = 1, dimensions(i)%length_mask-1
                   WRITE( debug_output_unit, '(A,I6)', ADVANCE='no' ) &
                      ',', dimensions(i)%masked_values_int16(j)
                ENDDO
                WRITE( debug_output_unit, '(A)' )  ''  ! write line-end
             ENDIF

          ELSEIF ( ALLOCATED( dimensions(i)%values_int32 ) )  THEN

             WRITE( debug_output_unit, TRIM( write_format2 ) // ',I11)', ADVANCE='no' ) &
                'values', dimensions(i)%values_int32(dimensions(i)%bounds(1))
             DO  j = dimensions(i)%bounds(1)+1, dimensions(i)%bounds(2)
                WRITE( debug_output_unit, '(A,I11)', ADVANCE='no' ) &
                   ',', dimensions(i)%values_int32(j)
             ENDDO
             WRITE( debug_output_unit, '(A)' )  ''  ! write line-end
             IF ( is_masked )  THEN
                WRITE( debug_output_unit, TRIM( write_format2 ) // ',I11)', ADVANCE='no' ) &
                   'masked values', dimensions(i)%masked_values_int32(0)
                DO  j = 1, dimensions(i)%length_mask-1
                   WRITE( debug_output_unit, '(A,I11)', ADVANCE='no' ) &
                      ',', dimensions(i)%masked_values_int32(j)
                ENDDO
                WRITE( debug_output_unit, '(A)' )  ''  ! write line-end
             ENDIF

          ELSEIF ( ALLOCATED( dimensions(i)%values_intwp ) )  THEN

             WRITE( debug_output_unit, TRIM( write_format2 ) // ',I11)', ADVANCE='no' ) &
                'values', dimensions(i)%values_intwp(dimensions(i)%bounds(1))
             DO  j = dimensions(i)%bounds(1)+1, dimensions(i)%bounds(2)
                WRITE( debug_output_unit, '(A,I11)', ADVANCE='no' ) &
                   ',', dimensions(i)%values_intwp(j)
             ENDDO
             WRITE( debug_output_unit, '(A)' )  ''  ! write line-end
             IF ( is_masked )  THEN
                WRITE( debug_output_unit, TRIM( write_format2 ) // ',I11)', ADVANCE='no' ) &
                   'masked values', dimensions(i)%masked_values_intwp(0)
                DO  j = 1, dimensions(i)%length_mask-1
                   WRITE( debug_output_unit, '(A,I11)', ADVANCE='no' ) &
                      ',', dimensions(i)%masked_values_intwp(j)
                ENDDO
                WRITE( debug_output_unit, '(A)' )  ''  ! write line-end
             ENDIF

          ELSEIF ( ALLOCATED( dimensions(i)%values_real32 ) )  THEN

             WRITE( debug_output_unit, TRIM( write_format2 ) // ',E14.7)', ADVANCE='no' ) &
                'values', dimensions(i)%values_real32(dimensions(i)%bounds(1))
             DO  j = dimensions(i)%bounds(1)+1, dimensions(i)%bounds(2)
                WRITE( debug_output_unit, '(A,E14.7)', ADVANCE='no' ) &
                   ',', dimensions(i)%values_real32(j)
             ENDDO
             WRITE( debug_output_unit, '(A)' )  ''  ! write line-end
             IF ( is_masked )  THEN
                WRITE( debug_output_unit, TRIM( write_format2 ) // ',E14.7)', ADVANCE='no' ) &
                   'masked values', dimensions(i)%masked_values_real32(0)
                DO  j = 1, dimensions(i)%length_mask-1
                   WRITE( debug_output_unit, '(A,E14.7)', ADVANCE='no' ) &
                      ',', dimensions(i)%masked_values_real32(j)
                ENDDO
                WRITE( debug_output_unit, '(A)' )  ''  ! write line-end
             ENDIF

          ELSEIF ( ALLOCATED( dimensions(i)%values_real64 ) )  THEN

             WRITE( debug_output_unit, TRIM( write_format2 ) // ',E22.15)', ADVANCE='no' ) &
                'values', dimensions(i)%values_real64(dimensions(i)%bounds(1))
             DO  j = dimensions(i)%bounds(1)+1, dimensions(i)%bounds(2)
                WRITE( debug_output_unit, '(A,E22.15)', ADVANCE='no' ) &
                   ',', dimensions(i)%values_real64(j)
             ENDDO
             WRITE( debug_output_unit, '(A)' )  ''  ! write line-end
             IF ( is_masked )  THEN
                WRITE( debug_output_unit, TRIM( write_format2 ) // ',E22.15)', ADVANCE='no' ) &
                   'masked values', dimensions(i)%masked_values_real64(0)
                DO  j = 1, dimensions(i)%length_mask-1
                   WRITE( debug_output_unit, '(A,E22.15)', ADVANCE='no' ) &
                      ',', dimensions(i)%masked_values_real64(j)
                ENDDO
                WRITE( debug_output_unit, '(A)' )  ''  ! write line-end
             ENDIF

          ELSEIF ( ALLOCATED( dimensions(i)%values_realwp ) )  THEN

             WRITE( debug_output_unit, TRIM( write_format2 ) // ',E22.15)', ADVANCE='no' ) &
                'values', dimensions(i)%values_realwp(dimensions(i)%bounds(1))
             DO  j = dimensions(i)%bounds(1)+1, dimensions(i)%bounds(2)
                WRITE( debug_output_unit, '(A,E22.15)', ADVANCE='no' ) &
                   ',', dimensions(i)%values_realwp(j)
             ENDDO
             WRITE( debug_output_unit, '(A)' )  ''  ! write line-end
             IF ( is_masked )  THEN
                WRITE( debug_output_unit, TRIM( write_format2 ) // ',E22.15)', ADVANCE='no' ) &
                   'masked values', dimensions(i)%masked_values_realwp(0)
                DO  j = 1, dimensions(i)%length_mask-1
                   WRITE( debug_output_unit, '(A,E22.15)', ADVANCE='no' ) &
                      ',', dimensions(i)%masked_values_realwp(j)
                ENDDO
                WRITE( debug_output_unit, '(A)' )  ''  ! write line-end
             ENDIF

          ENDIF

          IF ( ALLOCATED( dimensions(i)%attributes ) )  &
             CALL print_attributes( indent_level+1, dimensions(i)%attributes )

          IF ( i < nelement )  &
             WRITE( debug_output_unit, TRIM( write_format1 ) // ')' ) separation_string
       ENDDO

    END SUBROUTINE print_dimensions

!--------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Print list of variables.
!--------------------------------------------------------------------------------------------!
    SUBROUTINE print_variables( indent_level, variables )

       CHARACTER(LEN=50) ::  write_format1  !< format for write statements
       CHARACTER(LEN=50) ::  write_format2  !< format for write statements

       INTEGER             ::  i                        !< loop index
       INTEGER, INTENT(IN) ::  indent_level             !< indentation level
       INTEGER             ::  j                        !< loop index
       INTEGER, PARAMETER  ::  max_keyname_length = 16  !< length of longest key name
       INTEGER             ::  nelement                 !< number of elements to print

       TYPE(variable_type), DIMENSION(:), INTENT(IN) ::  variables  !< list of variables


       WRITE( write_format1, '(A,I3,A)' ) '(T', indent_level * indent_depth + 1, ',A'
       WRITE( write_format2, '(A,I3,A,I3,A)' ) &
          '(T', ( indent_level + 1 ) * indent_depth + 1, ',A,T', &
          ( indent_level + 1 ) * indent_depth + 1 + max_keyname_length, ',(": ")'

       WRITE( debug_output_unit, TRIM( write_format1 ) // ',A)' ) &
          REPEAT( separation_string // ' ', 4 )
       WRITE( debug_output_unit, TRIM( write_format1 ) // ')' ) 'variables:'

       nelement = SIZE( variables )
       DO  i = 1, nelement
          WRITE( debug_output_unit, TRIM( write_format2 ) // ',A)' ) &
             'name', TRIM( variables(i)%name )
          WRITE( debug_output_unit, TRIM( write_format2 ) // ',A)' ) &
             'type', TRIM( variables(i)%data_type )
          WRITE( debug_output_unit, TRIM( write_format2 ) // ',I7)' ) &
             'id', variables(i)%id
          WRITE( debug_output_unit, TRIM( write_format2 ) // ',L1)' ) &
             'is global', variables(i)%is_global

          WRITE( debug_output_unit, TRIM( write_format2 ) // ',A)', ADVANCE='no' ) &
             'dimension names', TRIM( variables(i)%dimension_names(1) )
          DO  j = 2, SIZE( variables(i)%dimension_names )
             WRITE( debug_output_unit, '(A,A)', ADVANCE='no' ) &
                ',', TRIM( variables(i)%dimension_names(j) )
          ENDDO
          WRITE( debug_output_unit, '(A)' )  ''  ! write line-end

          WRITE( debug_output_unit, TRIM( write_format2 ) // ',I7)', ADVANCE='no' ) &
             'dimension ids', variables(i)%dimension_ids(1)
          DO  j = 2, SIZE( variables(i)%dimension_names )
             WRITE( debug_output_unit, '(A,I7)', ADVANCE='no' ) &
                ',', variables(i)%dimension_ids(j)
          ENDDO
          WRITE( debug_output_unit, '(A)' )  ''  ! write line-end

          IF ( ALLOCATED( variables(i)%attributes ) )  &
             CALL print_attributes( indent_level+1, variables(i)%attributes )
          IF ( i < nelement )  &
             WRITE( debug_output_unit, TRIM( write_format1 ) // ')' ) separation_string
       ENDDO

    END SUBROUTINE print_variables

 END SUBROUTINE dom_database_debug_output

 END MODULE data_output_module
