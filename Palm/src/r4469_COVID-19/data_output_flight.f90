!> @file data_output_flight.f90
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
! $Id: data_output_flight.f90 4360 2020-01-07 11:25:50Z suehring $
! Corrected "Former revisions" section
! 
! 3655 2019-01-07 16:51:22Z knoop
! unused variables removed
! 
! 1957 2016-07-07 10:43:48Z suehring
! Initial revision
!
!
! Description:
! ------------
!> Writing data from flight measurements on file.
!------------------------------------------------------------------------------!
 SUBROUTINE data_output_flight
 
#if defined( __netcdf )
    USE control_parameters,                                                    &
        ONLY:  num_leg, num_var_fl, time_since_reference_point, virtual_flight

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point

    USE flight_mod,                                                            &
        ONLY:  sensor, x_pos, y_pos, z_pos

    USE kinds

    USE NETCDF

    USE netcdf_interface,                                                      &
        ONLY:  dofl_time_count, id_set_fl, id_var_dofl, id_var_time_fl,        &
               id_var_x_fl, id_var_y_fl, id_var_z_fl, netcdf_handle_error,     &
               nc_stat 
              
    USE pegrid

    IMPLICIT NONE

    INTEGER(iwp) ::  i !< loop variable for output quantities
    INTEGER(iwp) ::  l !< loop variable for flight legs
    INTEGER(iwp) ::  k !< internal count for variable labels and units

    CALL cpu_log( log_point(64), 'data_output_flight', 'start' )

!
!-- Check if virtual flights are carried out
    IF ( .NOT. virtual_flight )  RETURN

!
!-- Output is only performed on PE0
    IF ( myid == 0 )  THEN

!
!--    Open file for flight output in NetCDF format
       CALL check_open( 199 )

!
!--    Increment the counter for number of output times
       dofl_time_count = dofl_time_count + 1

!
!--    Update the flight-output time and spatial coordinates
       nc_stat = NF90_PUT_VAR( id_set_fl, id_var_time_fl,                      &
                               (/ time_since_reference_point /),               &
                               start = (/ dofl_time_count /),                  &
                               count = (/ 1 /) )
       CALL netcdf_handle_error( 'data_output_flight', 554 )

       DO  l = 1, num_leg
          nc_stat = NF90_PUT_VAR( id_set_fl, id_var_x_fl(l),                   &
                                  (/ x_pos(l) /),                              &
                                  start = (/ dofl_time_count /),               &
                                  count = (/ 1 /) )
          nc_stat = NF90_PUT_VAR( id_set_fl, id_var_y_fl(l),                   &
                                  (/ y_pos(l) /),                              &
                                  start = (/ dofl_time_count /),               &
                                  count = (/ 1 /) )
          nc_stat = NF90_PUT_VAR( id_set_fl, id_var_z_fl(l),                   &
                                  (/ z_pos(l) /),                              &
                                  start = (/ dofl_time_count /),               &
                                  count = (/ 1 /) )
          CALL netcdf_handle_error( 'data_output_flight', 555 )
       ENDDO
!
!--    Output measured quantities
       k = 1
       DO  l = 1, num_leg
          DO i = 1, num_var_fl
             nc_stat = NF90_PUT_VAR( id_set_fl, id_var_dofl(k),                &
                                     (/ sensor(i,l) /),                        &
                                     start = (/ dofl_time_count /),            &
                                     count = (/ 1 /) )

             CALL netcdf_handle_error( 'data_output_flight', 556 )
          
             k = k + 1
          ENDDO
       ENDDO
    ENDIF
    
    CALL cpu_log( log_point(64), 'data_output_flight', 'stop' )
    
#endif
 END SUBROUTINE data_output_flight
