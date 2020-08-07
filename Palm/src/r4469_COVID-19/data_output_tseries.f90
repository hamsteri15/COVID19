!> @file data_output_tseries.f90
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
! -----------------
! 
! 
! Former revisions:
! -----------------
! $Id: data_output_tseries.f90 4360 2020-01-07 11:25:50Z suehring $
! Corrected "Former revisions" section
! 
! 3655 2019-01-07 16:51:22Z knoop
! unused format removed
!
! Revision 1.1  1998/03/03 08:00:13  raasch
! Initial revision
!
!
! Description:
! ------------
!> Time series output for PROFIL. Always all time series are stored. A selection
!> can be applied via the PROFIL-parameters in close_file.
!------------------------------------------------------------------------------!
 SUBROUTINE data_output_tseries
 

    USE control_parameters,                                                    &
        ONLY:  dots_time_count, time_since_reference_point

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point  

    USE kinds

#if defined( __netcdf )
    USE NETCDF
#endif
    USE netcdf_interface,                                                      &
        ONLY:  dots_num, id_set_ts, id_var_dots, id_var_time_ts, nc_stat,      &
               netcdf_handle_error

    USE pegrid

    USE profil_parameter
    
    USE statistics,                                                            &
        ONLY:  flow_statistics_called, statistic_regions, ts_value

    IMPLICIT NONE


    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  sr !<


!
!-- If required, compute statistics.
    IF ( .NOT. flow_statistics_called )  CALL flow_statistics

!
!-- Flow_statistics has its own cpu-time measuring.
    CALL cpu_log( log_point(21), 'data_output_tseries', 'start' )

    IF ( myid == 0 )  THEN

!
!--    Open file for time series output in NetCDF format
       CALL check_open( 105 )
       
!--    Increment the counter for number of output times
!      CAUTION: The following line has to be after the call of the subroutine
!               check_open, since check_open resets the counter dots_time_count
!               to 0, if a new file is opened
       dots_time_count = dots_time_count + 1
       
#if defined( __netcdf )
!
!--    Update the time series time axis
       nc_stat = NF90_PUT_VAR( id_set_ts, id_var_time_ts,        &
                               (/ time_since_reference_point /), &
                               start = (/ dots_time_count /),    &
                               count = (/ 1 /) )
       CALL netcdf_handle_error( 'data_output_tseries', 350 )
#endif

!
!--    Time series output for the total domain (and each subregion, if
!--    applicable)
       DO  sr = 0, statistic_regions

#if defined( __netcdf )
          DO  i = 1, dots_num
             nc_stat = NF90_PUT_VAR( id_set_ts, id_var_dots(i,sr),  &
                                     (/ ts_value(i,sr) /),          &
                                     start = (/ dots_time_count /), &
                                     count = (/ 1 /) )
             CALL netcdf_handle_error( 'data_output_tseries', 351 )
          ENDDO
#endif

       ENDDO

    ENDIF


    CALL cpu_log( log_point(21), 'data_output_tseries', 'stop' )

 END SUBROUTINE data_output_tseries
