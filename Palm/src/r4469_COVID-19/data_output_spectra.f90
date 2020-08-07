!> @file data_output_spectra.f90
!------------------------------------------------------------------------------!
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
! Copyright 1997-2020 Leibniz Universitaet Hannover
!------------------------------------------------------------------------------!
!
! Current revisions:
! ------------------
! 
! 
! Former revisions:
! -----------------
! $Id: data_output_spectra.f90 4360 2020-01-07 11:25:50Z suehring $
! Corrected "Former revisions" section
! 
! 3655 2019-01-07 16:51:22Z knoop
! variables documented
!
! Revision 1.1  2001/01/05 15:14:20  raasch
! Initial revision
!
!
! Description:
! ------------
!> Writing spectra data on file, using a special format which allows
!> plotting of these data with PROFIL-graphic-software
!------------------------------------------------------------------------------!
 SUBROUTINE data_output_spectra
 
#if defined( __netcdf )
    USE control_parameters,                                                    &
        ONLY:  message_string, time_since_reference_point

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point

    USE kinds

    USE NETCDF

    USE netcdf_interface,                                                      &
        ONLY:  id_set_sp, id_var_time_sp, nc_stat, netcdf_handle_error

    USE pegrid

    USE spectra_mod,                                                           &
        ONLY:  average_count_sp, averaging_interval_sp, comp_spectra_level,    &
               data_output_sp, dosp_time_count, spectra_direction, spectrum_x, &
               spectrum_y


    IMPLICIT NONE

    INTEGER(iwp) ::  m       !< running index over spectra output
    INTEGER(iwp) ::  pr      !< index used to assign default quantities to data output
    
    CALL cpu_log( log_point(31), 'data_output_spectra', 'start' )

!
!-- Check if user gave any levels for spectra to be calculated
    IF ( comp_spectra_level(1) == 999999 )  RETURN

!
!-- Output is only performed on PE0
    IF ( myid == 0 )  THEN

!
!--    Open file for spectra output in NetCDF format
       CALL check_open( 107 )

!
!--    Increment the counter for number of output times
       dosp_time_count = dosp_time_count + 1

!
!--    Update the spectra time axis
       nc_stat = NF90_PUT_VAR( id_set_sp, id_var_time_sp,        &
                               (/ time_since_reference_point /), &
                               start = (/ dosp_time_count /), count = (/ 1 /) )
       CALL netcdf_handle_error( 'data_output_spectra', 47 )

!
!--    If necessary, calculate time average and reset average counter
       IF ( average_count_sp == 0 )  THEN
           message_string = 'no spectra data available'
           CALL message( 'data_output_spectra', 'PA0186', 0, 0, 0, 6, 0 )
       ENDIF
       IF ( average_count_sp /= 1 )  THEN
          spectrum_x = spectrum_x / REAL( average_count_sp, KIND=wp )
          spectrum_y = spectrum_y / REAL( average_count_sp, KIND=wp )
          average_count_sp = 0
       ENDIF

!
!--    Loop over all spectra defined by the user
       m = 1
       DO WHILE ( data_output_sp(m) /= ' '  .AND.  m <= 10 )

          SELECT CASE ( TRIM( data_output_sp(m) ) )

             CASE ( 'u' )
                pr = 1

             CASE ( 'v' )
                pr = 2

             CASE ( 'w' )
                pr = 3

             CASE ( 'theta' )
                pr = 4

             CASE ( 'q' )
                pr = 5

             CASE ( 's' )
                pr = 6

             CASE DEFAULT
!
!--             The DEFAULT case is reached either if the parameter 
!--             data_output_sp(m) contains a wrong character string or if the 
!--             user has coded a special case in the user interface. There, the 
!--             subroutine user_spectra checks which of these two conditions 
!--             applies.
                CALL user_spectra( 'data_output', m, pr )

          END SELECT

!
!--       Output of spectra in NetCDF format
!--       Output of x-spectra
          IF ( INDEX( spectra_direction(m), 'x' ) /= 0 ) THEN
             CALL output_spectra_netcdf( m, 'x' )
          ENDIF
!
!--       Output of y-spectra
          IF ( INDEX( spectra_direction(m), 'y' ) /= 0 ) THEN
             CALL output_spectra_netcdf( m, 'y' )
          ENDIF

!
!--       Increase counter for next spectrum
          m = m + 1

       ENDDO

!
!--    Reset spectra values
       spectrum_x = 0.0_wp; spectrum_y = 0.0_wp

    ENDIF

    CALL cpu_log( log_point(31), 'data_output_spectra', 'stop' )

#if defined( __parallel )
!    CALL MPI_BARRIER( comm2d, ierr )  ! really necessary
#endif

#endif
 END SUBROUTINE data_output_spectra


!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
 SUBROUTINE output_spectra_netcdf( nsp, direction )
#if defined( __netcdf )

    USE basic_constants_and_equations_mod,                                     &
        ONLY:  pi

    USE grid_variables,                                                        &
        ONLY:  dx, dy

    USE indices,                                                               &
        ONLY:  nx, ny

    USE kinds

    USE NETCDF

    USE netcdf_interface,                                                      &
        ONLY:  id_set_sp, id_var_dospx, id_var_dospy, nc_stat,                 &
               netcdf_handle_error

    USE spectra_mod,                                                           &
        ONLY:  dosp_time_count, n_sp_x, n_sp_y, spectrum_x, spectrum_y


    IMPLICIT NONE

    CHARACTER (LEN=1), INTENT(IN) ::  direction     !< directio of spectra evaluation

    INTEGER(iwp), INTENT(IN)      ::  nsp           !< number of spectrum

    INTEGER(iwp)                  ::  i             !< running index in frequency space
    INTEGER(iwp)                  ::  k             !< running index over number of spectrum

    REAL(wp)                      ::  frequency     !< wavenumber

    REAL(wp), DIMENSION(nx/2)     ::  netcdf_data_x !< normalized wavenumber along x written into NetCDF file
    REAL(wp), DIMENSION(ny/2)     ::  netcdf_data_y !< normalized wavenumber along y written into NetCDF file


    IF ( direction == 'x' )  THEN

       DO  k = 1, n_sp_x

          DO  i = 1, nx/2
             frequency = 2.0_wp * pi * i / ( dx * ( nx + 1 ) )
             netcdf_data_x(i) = frequency * spectrum_x(i,k,nsp)
          ENDDO

          nc_stat = NF90_PUT_VAR( id_set_sp, id_var_dospx(nsp), netcdf_data_x, &
                                  start = (/ 1, k, dosp_time_count /), &
                                  count = (/ nx/2, 1, 1 /) )
          CALL netcdf_handle_error( 'data_output_spectra', 348 )

       ENDDO

    ENDIF

    IF ( direction == 'y' )  THEN

       DO  k = 1, n_sp_y

          DO  i = 1, ny/2
             frequency = 2.0_wp * pi * i / ( dy * ( ny + 1 ) )
             netcdf_data_y(i) = frequency * spectrum_y(i,k,nsp)
          ENDDO

          nc_stat = NF90_PUT_VAR( id_set_sp, id_var_dospy(nsp), netcdf_data_y, &
                                  start = (/ 1, k, dosp_time_count /), &
                                  count = (/ ny/2, 1, 1 /) )
          CALL netcdf_handle_error( 'data_output_spectra', 349 )

       ENDDO

    ENDIF

#endif
 END SUBROUTINE output_spectra_netcdf
