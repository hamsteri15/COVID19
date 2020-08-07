!> @file data_output_profiles.f90
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
! $Id: data_output_profiles.f90 4360 2020-01-07 11:25:50Z suehring $
! Corrected "Former revisions" section
! 
! 3655 2019-01-07 16:51:22Z knoop
! add variable description
!
! Revision 1.1  1997/09/12 06:28:48  raasch
! Initial revision
!
!
! Description:
! ------------
!> Plot output of 1D-profiles for PROFIL
!------------------------------------------------------------------------------!
 SUBROUTINE data_output_profiles
 

    USE control_parameters,                                                    &
        ONLY:  average_count_pr, averaging_interval_pr, coupling_start_time,   &
               dopr_n, dopr_time_count, normalizing_region,                    &
               time_since_reference_point

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point

    USE indices,                                                               &
        ONLY:  nzb, nzt

    USE kinds

#if defined( __netcdf )
    USE NETCDF
#endif

    USE netcdf_interface,                                                      &
        ONLY:  id_set_pr, id_var_dopr, id_var_norm_dopr, id_var_time_pr,       &
               nc_stat, netcdf_handle_error, output_for_t0

    USE pegrid

    USE profil_parameter

    USE statistics,                                                            &
        ONLY:  flow_statistics_called, hom, hom_sum, pr_palm, statistic_regions

    IMPLICIT NONE


    INTEGER(iwp) ::  i  !< loop index
    INTEGER(iwp) ::  sr !< statistic region index

!
!-- If required, compute statistics
    IF ( .NOT. flow_statistics_called )  CALL flow_statistics

!
!-- Flow_statistics has its own CPU time measurement
    CALL cpu_log( log_point(15), 'data_output_profiles', 'start' )

!
!-- If required, compute temporal average
    IF ( averaging_interval_pr == 0.0_wp )  THEN
       hom_sum(:,:,:) = hom(:,1,:,:)
    ELSE
       IF ( average_count_pr > 0 )  THEN
          hom_sum = hom_sum / REAL( average_count_pr, KIND=wp )
       ELSE
!
!--       This case may happen if dt_dopr is changed in the 
!--       runtime_parameters-list of a restart run
          RETURN
       ENDIF
    ENDIF

    
    IF ( myid == 0 )  THEN

!
!--    Plot-output for each (sub-)region

!
!--    Open file for profile output in NetCDF format
       CALL check_open( 104 )

!
!--    Increment the counter for number of output times
       dopr_time_count = dopr_time_count + 1

!
!--    Output of initial profiles
       IF ( dopr_time_count == 1 )  THEN
        
          IF ( .NOT. output_for_t0 ) THEN 

#if defined( __netcdf )          
!
!--          Store initial time to time axis, but only if an output
!--          is required for at least one of the profiles. The initial time
!--          is either 0, or, in case of a prerun for coupled atmosphere-ocean
!--          runs, has a negative value
             DO  i = 1, dopr_n
                IF ( dopr_initial_index(i) /= 0 )  THEN
                   nc_stat = NF90_PUT_VAR( id_set_pr, id_var_time_pr,          &
                                        (/ -coupling_start_time /),            &
                                        start = (/ 1 /), count = (/ 1 /) )
                   CALL netcdf_handle_error( 'data_output_profiles', 329 )
                   output_for_t0 = .TRUE.
                   EXIT
                ENDIF
             ENDDO

!
!--          Store normalization factors
             nc_stat = NF90_PUT_VAR( id_set_pr, id_var_norm_dopr(1), & ! wpt0
                                  (/ hom_sum(nzb,18,normalizing_region) /), &
                                     start = (/ 1 /), count = (/ 1 /) )
             CALL netcdf_handle_error( 'data_output_profiles', 330 )

             nc_stat = NF90_PUT_VAR( id_set_pr, id_var_norm_dopr(2), & ! ws2
                        (/ hom_sum(nzb+8,pr_palm,normalizing_region)**2 /), &
                                     start = (/ 1 /), count = (/ 1 /) )
             CALL netcdf_handle_error( 'data_output_profiles', 331 )
             nc_stat = NF90_PUT_VAR( id_set_pr, id_var_norm_dopr(3), & ! tsw2
                        (/ hom_sum(nzb+3,pr_palm,normalizing_region)**2 /), &
                                  start = (/ 1 /), count = (/ 1 /) )
             CALL netcdf_handle_error( 'data_output_profiles', 332 )
             nc_stat = NF90_PUT_VAR( id_set_pr, id_var_norm_dopr(4), & ! ws3
                        (/ hom_sum(nzb+8,pr_palm,normalizing_region)**3 /), &
                                     start = (/ 1 /), count = (/ 1 /) )
             CALL netcdf_handle_error( 'data_output_profiles', 333 )

             nc_stat = NF90_PUT_VAR( id_set_pr, id_var_norm_dopr(5), &!ws2tsw
                        (/ hom_sum(nzb+8,pr_palm,normalizing_region)**3 *   &
                           hom_sum(nzb+3,pr_palm,normalizing_region)    /), &
                                     start = (/ 1 /), count = (/ 1 /) )
             CALL netcdf_handle_error( 'data_output_profiles', 334 )

             nc_stat = NF90_PUT_VAR( id_set_pr, id_var_norm_dopr(6), &!wstsw2
                        (/ hom_sum(nzb+8,pr_palm,normalizing_region) *      &
                           hom_sum(nzb+3,pr_palm,normalizing_region)**2 /), &
                                     start = (/ 1 /), count = (/ 1 /) )
             CALL netcdf_handle_error( 'data_output_profiles', 335 )

             nc_stat = NF90_PUT_VAR( id_set_pr, id_var_norm_dopr(7), & ! z_i
                           (/ hom_sum(nzb+6,pr_palm,normalizing_region) /), &
                                     start = (/ 1 /), count = (/ 1 /) )
             CALL netcdf_handle_error( 'data_output_profiles', 336 )
             
#endif
!
!--          Loop over all 1D variables
             DO  i = 1, dopr_n

                IF ( dopr_initial_index(i) /= 0 )  THEN

!
!--                Output for the individual (sub-)regions
                   DO  sr = 0, statistic_regions

#if defined( __netcdf )
!
!--                   Write data to netcdf file
                      nc_stat = NF90_PUT_VAR( id_set_pr, id_var_dopr(i,sr),    &
                                    hom(nzb:nzt+1,1,dopr_initial_index(i),sr), &
                                              start = (/ 1, 1 /),              &
                                              count = (/ nzt-nzb+2, 1 /) )
                      CALL netcdf_handle_error( 'data_output_profiles', 337 )
#endif

                   ENDDO

                ENDIF   ! Initial profile available

             ENDDO   ! Loop over dopr_n for initial profiles

             IF ( output_for_t0 )  THEN
                dopr_time_count = dopr_time_count + 1
             ENDIF

          END IF
       ENDIF   ! Initial profiles

#if defined( __netcdf )

!
!--    Store time to time axis
       nc_stat = NF90_PUT_VAR( id_set_pr, id_var_time_pr,        &
                               (/ time_since_reference_point /), &
                               start = (/ dopr_time_count /),    &
                               count = (/ 1 /) )
       CALL netcdf_handle_error( 'data_output_profiles', 338 )

!
!--    Store normalization factors
       nc_stat = NF90_PUT_VAR( id_set_pr, id_var_norm_dopr(1), &  ! wpt0
                               (/ hom_sum(nzb,18,normalizing_region) /), &
                               start = (/ dopr_time_count /),               &
                               count = (/ 1 /) )
       CALL netcdf_handle_error( 'data_output_profiles', 339 )

       nc_stat = NF90_PUT_VAR( id_set_pr, id_var_norm_dopr(2), &  ! ws2
                     (/ hom_sum(nzb+8,pr_palm,normalizing_region)**2 /), &
                               start = (/ dopr_time_count /),               &
                               count = (/ 1 /) )
       CALL netcdf_handle_error( 'data_output_profiles', 340 )

       nc_stat = NF90_PUT_VAR( id_set_pr, id_var_norm_dopr(3), &  ! tsw2
                     (/ hom_sum(nzb+3,pr_palm,normalizing_region)**2 /), &
                               start = (/ dopr_time_count /),               &
                               count = (/ 1 /) )
       CALL netcdf_handle_error( 'data_output_profiles', 341 )

       nc_stat = NF90_PUT_VAR( id_set_pr, id_var_norm_dopr(4), &  ! ws3
                     (/ hom_sum(nzb+8,pr_palm,normalizing_region)**3 /), &
                               start = (/ dopr_time_count /),               &
                               count = (/ 1 /) )
       CALL netcdf_handle_error( 'data_output_profiles', 342 )

       nc_stat = NF90_PUT_VAR( id_set_pr, id_var_norm_dopr(5), &  ! ws2tsw
                     (/ hom_sum(nzb+8,pr_palm,normalizing_region)**3 *   &
                        hom_sum(nzb+3,pr_palm,normalizing_region)    /), &
                               start = (/ dopr_time_count /),               &
                               count = (/ 1 /) )
       CALL netcdf_handle_error( 'data_output_profiles', 343 )

       nc_stat = NF90_PUT_VAR( id_set_pr, id_var_norm_dopr(6), &  ! wstsw2
                     (/ hom_sum(nzb+8,pr_palm,normalizing_region) *      &
                        hom_sum(nzb+3,pr_palm,normalizing_region)**2 /), &
                               start = (/ dopr_time_count /),               &
                               count = (/ 1 /) )
       CALL netcdf_handle_error( 'data_output_profiles', 344 )

       nc_stat = NF90_PUT_VAR( id_set_pr, id_var_norm_dopr(7), &  ! z_i
                        (/ hom_sum(nzb+6,pr_palm,normalizing_region) /), &
                               start = (/ dopr_time_count /),               &
                               count = (/ 1 /) )
       CALL netcdf_handle_error( 'data_output_profiles', 345 )
#endif

!
!--    Output of the individual (non-initial) profiles
       DO  i = 1, dopr_n

!
!--       Output for the individual (sub-)domains
          DO  sr = 0, statistic_regions

#if defined( __netcdf )
!
!--          Write data to netcdf file
             nc_stat = NF90_PUT_VAR( id_set_pr, id_var_dopr(i,sr),          &
                                     hom_sum(nzb:nzt+1,dopr_index(i),sr),&
                                     start = (/ 1, dopr_time_count /),      &
                                     count = (/ nzt-nzb+2, 1 /) )
             CALL netcdf_handle_error( 'data_output_profiles', 346 )
#endif

          ENDDO

       ENDDO

    ENDIF  ! Output on PE0

!
!-- If averaging has been done above, the summation counter must be re-set.
    IF ( averaging_interval_pr /= 0.0_wp )  THEN
       average_count_pr = 0
    ENDIF

    CALL cpu_log( log_point(15), 'data_output_profiles','stop' )

 END SUBROUTINE data_output_profiles
