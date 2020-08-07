!> @file diagnostic_output_quantities_mod.f90
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
! $Id: diagnostic_output_quantities_mod.f90 4457 2020-03-11 14:20:43Z raasch $
! use statement for exchange horiz added,
! bugfix for call of exchange horiz 2d
! 
! 4431 2020-02-27 23:23:01Z gronemeier
! added wspeed and wdir output; bugfix: set fill_value in case of masked output
!
! 4360 2020-01-07 11:25:50Z suehring
! added output of wu, wv, wtheta and wq to enable covariance calculation
! according to temporal EC method
!
! 4346 2019-12-18 11:55:56Z motisi
! Introduction of wall_flags_total_0, which currently sets bits based on static
! topography information used in wall_flags_static_0
!
! 4331 2019-12-10 18:25:02Z suehring
! - Modularize 2-m potential temperature output
! - New output for 10-m wind speed
!
! 4329 2019-12-10 15:46:36Z motisi
! Renamed wall_flags_0 to wall_flags_static_0
!
! 4182 2019-08-22 15:20:23Z scharf
! Corrected "Former revisions" section
!
! 4167 2019-08-16 11:01:48Z suehring
! Changed behaviour of masked output over surface to follow terrain and ignore
! buildings (J.Resler, T.Gronemeier)
!
! 4157 2019-08-14 09:19:12Z suehring
! Initialization restructured, in order to work also when data output during
! spin-up is enabled.
!
! 4132 2019-08-02 12:34:17Z suehring
! Bugfix in masked data output
!
! 4069 2019-07-01 14:05:51Z Giersch
! Masked output running index mid has been introduced as a local variable to
! avoid runtime error (Loop variable has been modified) in time_integration
!
! 4039 2019-06-18 10:32:41Z suehring
! - Add output of uu, vv, ww to enable variance calculation according temporal
!   EC method
! - Allocate arrays only when they are required
! - Formatting adjustment
! - Rename subroutines
! - Further modularization
!
! 3998 2019-05-23 13:38:11Z suehring
! Bugfix in gathering all output strings
!
! 3995 2019-05-22 18:59:54Z suehring
! Avoid compiler warnings about unused variable and fix string operation which
! is not allowed with PGI compiler
!
! 3994 2019-05-22 18:08:09Z suehring
! Initial revision
!
! Authors:
! --------
! @author Farah Kanani-Suehring
!
!
! Description:
! ------------
!> ...
!------------------------------------------------------------------------------!
 MODULE diagnostic_output_quantities_mod

    USE arrays_3d,                                                             &
        ONLY:  ddzu,                                                           &
               pt,                                                             &
               q,                                                              &
               u,                                                              &
               v,                                                              &
               w,                                                              &
               zu,                                                             &
               zw

    USE basic_constants_and_equations_mod,                                     &
        ONLY:  kappa, pi

    USE control_parameters,                                                    &
        ONLY:  current_timestep_number,                                        &
               data_output,                                                    &
               message_string,                                                 &
               varnamelength
!
!     USE cpulog,                                                                &
!         ONLY:  cpu_log, log_point

    USE exchange_horiz_mod,                                                    &
        ONLY:  exchange_horiz_2d

   USE grid_variables,                                                         &
        ONLY:  ddx, ddy

    USE indices,                                                               &
        ONLY:  nbgp,                                                           &
               nxl,                                                            &
               nxlg,                                                           &
               nxr,                                                            &
               nxrg,                                                           &
               nyn,                                                            &
               nyng,                                                           &
               nys,                                                            &
               nysg,                                                           &
               nzb,                                                            &
               nzt,                                                            &
               wall_flags_total_0

    USE kinds

    USE surface_mod,                                                           &
        ONLY:  surf_def_h,                                                     &
               surf_lsm_h,                                                     &
               surf_type,                                                      &
               surf_usm_h


    IMPLICIT NONE

    CHARACTER(LEN=varnamelength), DIMENSION(500) ::  do_all = ' '

    INTEGER(iwp) ::  timestep_number_at_prev_calc = 0  !< ...at previous diagnostic output calculation

    LOGICAL ::  initialized_diagnostic_output_quantities = .FALSE. !< flag indicating whether output is initialized
    LOGICAL ::  prepared_diagnostic_output_quantities = .FALSE.    !< flag indicating whether output is p

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  pt_2m     !< 2-m air potential temperature
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  pt_2m_av  !< averaged 2-m air potential temperature
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  uv_10m    !< horizontal wind speed at 10m
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  uv_10m_av !< averaged horizontal wind speed at 10m

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  ti     !< rotation(u,v,w) aka turbulence intensity
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  ti_av  !< avg. rotation(u,v,w) aka turbulence intensity
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  u_center     !< u at center of grid box
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  u_center_av  !< mean of u_center
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  uu           !< uu
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  uu_av        !< mean of uu
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  wspeed       !< horizontal wind speed
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  wspeed_av    !< mean of horizotal wind speed
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  v_center     !< v at center of grid box
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  v_center_av  !< mean of v_center
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  vv           !< vv
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  vv_av        !< mean of vv
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  wdir         !< wind direction
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  wdir_av      !< mean wind direction
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  ww           !< ww
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  ww_av        !< mean of ww
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  wu           !< wu
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  wu_av        !< mean of wu
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  wv           !< wv
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  wv_av        !< mean of wv
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  wtheta       !< wtheta
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  wtheta_av    !< mean of wtheta
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  wq           !< wq
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  wq_av        !< mean of wq


    SAVE

    PRIVATE

!
!-- Public variables
    PUBLIC do_all,                                                             &
           initialized_diagnostic_output_quantities,                           &
           prepared_diagnostic_output_quantities,                              &
           timestep_number_at_prev_calc,                                       &
           pt_2m_av,                                                           &
           ti_av,                                                              &
           u_center_av,                                                        &
           uu_av,                                                              &
           uv_10m_av,                                                          &
           v_center_av,                                                        &
           vv_av,                                                              &
           wdir_av,                                                            &
           wspeed_av,                                                          &
           ww_av
!
!-- Public routines
    PUBLIC doq_3d_data_averaging,                                              &
           doq_calculate,                                                      &
           doq_check_data_output,                                              &
           doq_define_netcdf_grid,                                             &
           doq_output_2d,                                                      &
           doq_output_3d,                                                      &
           doq_output_mask,                                                    &
           doq_init,                                                           &
           doq_wrd_local
!          doq_rrd_local,                                                      &


    INTERFACE doq_3d_data_averaging
       MODULE PROCEDURE doq_3d_data_averaging
    END INTERFACE doq_3d_data_averaging

    INTERFACE doq_calculate
       MODULE PROCEDURE doq_calculate
    END INTERFACE doq_calculate

    INTERFACE doq_check_data_output
       MODULE PROCEDURE doq_check_data_output
    END INTERFACE doq_check_data_output

    INTERFACE doq_define_netcdf_grid
       MODULE PROCEDURE doq_define_netcdf_grid
    END INTERFACE doq_define_netcdf_grid

    INTERFACE doq_output_2d
       MODULE PROCEDURE doq_output_2d
    END INTERFACE doq_output_2d

    INTERFACE doq_output_3d
       MODULE PROCEDURE doq_output_3d
    END INTERFACE doq_output_3d

    INTERFACE doq_output_mask
       MODULE PROCEDURE doq_output_mask
    END INTERFACE doq_output_mask

    INTERFACE doq_init
       MODULE PROCEDURE doq_init
    END INTERFACE doq_init

    INTERFACE doq_prepare
       MODULE PROCEDURE doq_prepare
    END INTERFACE doq_prepare

!     INTERFACE doq_rrd_local
!        MODULE PROCEDURE doq_rrd_local
!     END INTERFACE doq_rrd_local

    INTERFACE doq_wrd_local
       MODULE PROCEDURE doq_wrd_local
    END INTERFACE doq_wrd_local


 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Sum up and time-average diagnostic output quantities as well as allocate
!> the array necessary for storing the average.
!------------------------------------------------------------------------------!
 SUBROUTINE doq_3d_data_averaging( mode, variable )

    USE control_parameters,                                                    &
        ONLY:  average_count_3d

    CHARACTER (LEN=*) ::  mode     !<
    CHARACTER (LEN=*) ::  variable !<

    INTEGER(iwp) ::  i !<
    INTEGER(iwp) ::  j !<
    INTEGER(iwp) ::  k !<

    IF ( mode == 'allocate' )  THEN

       SELECT CASE ( TRIM( variable ) )

          CASE ( 'ti' )
             IF ( .NOT. ALLOCATED( ti_av ) )  THEN
                ALLOCATE( ti_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
             ENDIF
             ti_av = 0.0_wp

          CASE ( 'uu' )
             IF ( .NOT. ALLOCATED( uu_av ) )  THEN
                ALLOCATE( uu_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
             ENDIF
             uu_av = 0.0_wp

          CASE ( 'vv' )
             IF ( .NOT. ALLOCATED( vv_av ) )  THEN
                ALLOCATE( vv_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
             ENDIF
             vv_av = 0.0_wp

          CASE ( 'ww' )
             IF ( .NOT. ALLOCATED( ww_av ) )  THEN
                ALLOCATE( ww_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
             ENDIF
             ww_av = 0.0_wp

           CASE ( 'wu' )
             IF ( .NOT. ALLOCATED( wu_av ) )  THEN
                ALLOCATE( wu_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
             ENDIF
             wu_av = 0.0_wp

           CASE ( 'wv' )
             IF ( .NOT. ALLOCATED( wv_av ) )  THEN
                ALLOCATE( wv_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
             ENDIF
             wv_av = 0.0_wp

           CASE ( 'wtheta' )
             IF ( .NOT. ALLOCATED( wtheta_av ) )  THEN
                ALLOCATE( wtheta_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
             ENDIF
             wtheta_av = 0.0_wp

           CASE ( 'wq' )
             IF ( .NOT. ALLOCATED( wq_av ) )  THEN
                ALLOCATE( wq_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
             ENDIF
             wq_av = 0.0_wp

          CASE ( 'theta_2m*' )
             IF ( .NOT. ALLOCATED( pt_2m_av ) )  THEN
                ALLOCATE( pt_2m_av(nysg:nyng,nxlg:nxrg) )
             ENDIF
             pt_2m_av = 0.0_wp

          CASE ( 'wspeed_10m*' )
             IF ( .NOT. ALLOCATED( uv_10m_av ) )  THEN
                ALLOCATE( uv_10m_av(nysg:nyng,nxlg:nxrg) )
             ENDIF
             uv_10m_av = 0.0_wp

          CASE ( 'wspeed' )
             IF ( .NOT. ALLOCATED( wspeed_av ) )  THEN
                ALLOCATE( wspeed_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
             ENDIF
             wspeed_av = 0.0_wp

          CASE ( 'wdir' )
             IF ( .NOT. ALLOCATED( u_center_av ) )  THEN
                ALLOCATE( u_center_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
             ENDIF
             IF ( .NOT. ALLOCATED( v_center_av ) )  THEN
                ALLOCATE( v_center_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
             ENDIF
             u_center_av = 0.0_wp
             v_center_av = 0.0_wp

          CASE DEFAULT
             CONTINUE

       END SELECT

    ELSEIF ( mode == 'sum' )  THEN

       SELECT CASE ( TRIM( variable ) )

          CASE ( 'ti' )
             IF ( ALLOCATED( ti_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         ti_av(k,j,i) = ti_av(k,j,i) + ti(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'uu' )
             IF ( ALLOCATED( uu_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         uu_av(k,j,i) = uu_av(k,j,i) + uu(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'vv' )
             IF ( ALLOCATED( vv_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         vv_av(k,j,i) = vv_av(k,j,i) + vv(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'ww' )
             IF ( ALLOCATED( ww_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         ww_av(k,j,i) = ww_av(k,j,i) + ww(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'wu' )
             IF ( ALLOCATED( wu_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         wu_av(k,j,i) = wu_av(k,j,i) + wu(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'wv' )
             IF ( ALLOCATED( wv_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         wv_av(k,j,i) = wv_av(k,j,i) + wv(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'wtheta' )
             IF ( ALLOCATED( wtheta_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         wtheta_av(k,j,i) = wtheta_av(k,j,i) + wtheta(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'wq' )
             IF ( ALLOCATED( wq_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         wq_av(k,j,i) = wq_av(k,j,i) + wq(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'theta_2m*' )
             IF ( ALLOCATED( pt_2m_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      pt_2m_av(j,i) = pt_2m_av(j,i) + pt_2m(j,i)
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'wspeed_10m*' )
             IF ( ALLOCATED( uv_10m_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      uv_10m_av(j,i) = uv_10m_av(j,i) + uv_10m(j,i)
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'wspeed' )
            IF ( ALLOCATED( wspeed_av ) ) THEN
               DO  i = nxl, nxr
                  DO  j = nys, nyn
                     DO  k = nzb, nzt+1
                         wspeed_av(k,j,i) = wspeed_av(k,j,i) + wspeed(k,j,i)
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF

          CASE ( 'wdir' )
             IF ( ALLOCATED( u_center_av )  .AND.  ALLOCATED( v_center_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                        u_center_av(k,j,i) = u_center_av(k,j,i) + u_center(k,j,i)
                        v_center_av(k,j,i) = v_center_av(k,j,i) + v_center(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE DEFAULT
             CONTINUE

       END SELECT

    ELSEIF ( mode == 'average' )  THEN

       SELECT CASE ( TRIM( variable ) )

          CASE ( 'ti' )
             IF ( ALLOCATED( ti_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         ti_av(k,j,i) = ti_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'uu' )
             IF ( ALLOCATED( uu_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         uu_av(k,j,i) = uu_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'vv' )
             IF ( ALLOCATED( vv_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         vv_av(k,j,i) = vv_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'ww' )
             IF ( ALLOCATED( ww_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         ww_av(k,j,i) = ww_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'wu' )
             IF ( ALLOCATED( wu_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         wu_av(k,j,i) = wu_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'wv' )
             IF ( ALLOCATED( wv_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         wv_av(k,j,i) = wv_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'wtheta' )
             IF ( ALLOCATED( wtheta_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         wtheta_av(k,j,i) = wtheta_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'wq' )
             IF ( ALLOCATED( wq_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         wq_av(k,j,i) = wq_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

         CASE ( 'theta_2m*' )
            IF ( ALLOCATED( pt_2m_av ) ) THEN
               DO  i = nxlg, nxrg
                  DO  j = nysg, nyng
                     pt_2m_av(j,i) = pt_2m_av(j,i) / REAL( average_count_3d, KIND=wp )
                  ENDDO
               ENDDO
               CALL exchange_horiz_2d( pt_2m_av )
            ENDIF

         CASE ( 'wspeed_10m*' )
            IF ( ALLOCATED( uv_10m_av ) ) THEN
               DO  i = nxlg, nxrg
                  DO  j = nysg, nyng
                     uv_10m_av(j,i) = uv_10m_av(j,i) / REAL( average_count_3d, KIND=wp )
                  ENDDO
               ENDDO
               CALL exchange_horiz_2d( uv_10m_av )
            ENDIF

         CASE ( 'wspeed' )
             IF ( ALLOCATED( wspeed_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         wspeed_av(k,j,i) = wspeed_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'wdir' )
             IF ( ALLOCATED( u_center_av )  .AND.  ALLOCATED( v_center_av ) ) THEN

                IF ( .NOT. ALLOCATED( wdir_av ) )  THEN
                   ALLOCATE( wdir_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
                ENDIF
                wdir_av = 0.0_wp

                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         u_center_av(k,j,i) = u_center_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                         v_center_av(k,j,i) = v_center_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                         wdir_av(k,j,i) = ATAN2( u_center_av(k,j,i), v_center_av(k,j,i) ) &
                                        / pi * 180.0_wp + 180.0_wp
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

       END SELECT

    ENDIF


 END SUBROUTINE doq_3d_data_averaging

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check data output for diagnostic output
!------------------------------------------------------------------------------!
 SUBROUTINE doq_check_data_output( var, unit, i, ilen, k )

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  unit  !<
    CHARACTER (LEN=*) ::  var   !<

    INTEGER(iwp), OPTIONAL, INTENT(IN) ::  i     !< Current element of data_output
    INTEGER(iwp), OPTIONAL, INTENT(IN) ::  ilen  !< Length of current entry in data_output
    INTEGER(iwp), OPTIONAL, INTENT(IN) ::  k     !< Output is xy mode? 0 = no, 1 = yes

    SELECT CASE ( TRIM( var ) )

       CASE ( 'ti' )
          unit = '1/s'

       CASE ( 'uu' )
          unit = 'm2/s2'

       CASE ( 'vv' )
          unit = 'm2/s2'

       CASE ( 'ww' )
          unit = 'm2/s2'

       CASE ( 'wu' )
          unit = 'm2/s2'

       CASE ( 'wv' )
          unit = 'm2/s2'

       CASE ( 'wtheta' )
          unit = 'Km/s'

       CASE ( 'wq' )
          unit = 'm/s'

       CASE ( 'wspeed' )
          unit = 'm/s'

       CASE ( 'wdir' )
          unit = 'degree'
!
!--    Treat horizotal cross-section output quanatities
       CASE ( 'theta_2m*', 'wspeed_10m*' )
!
!--       Check if output quantity is _xy only.
          IF ( k == 0  .OR.  data_output(i)(ilen-2:ilen) /= '_xy' )  THEN
             message_string = 'illegal value for data_output: "' //            &
                              TRIM( var ) // '" & only 2d-horizontal ' //      &
                              'cross sections are allowed for this value'
             CALL message( 'diagnostic_output', 'PA0111', 1, 2, 0, 6, 0 )
          ENDIF

          IF ( TRIM( var ) == 'theta_2m*'   )  unit = 'K'
          IF ( TRIM( var ) == 'wspeed_10m*' )  unit = 'm/s'

       CASE DEFAULT
          unit = 'illegal'

    END SELECT


 END SUBROUTINE doq_check_data_output

!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining appropriate grid for netcdf variables.
!------------------------------------------------------------------------------!
 SUBROUTINE doq_define_netcdf_grid( variable, found, grid_x, grid_y, grid_z )

    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT(IN)  ::  variable    !<
    LOGICAL, INTENT(OUT)           ::  found       !<
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_x      !<
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_y      !<
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_z      !<

    found  = .TRUE.

    SELECT CASE ( TRIM( variable ) )
!
!--    s grid
       CASE ( 'ti', 'ti_xy', 'ti_xz', 'ti_yz',                                 &
              'wspeed', 'wspeed_xy', 'wspeed_xz', 'wspeed_yz',                 &
              'wdir', 'wdir_xy', 'wdir_xz', 'wdir_yz',                         &
              'wu', 'wu_xy', 'wu_xz', 'wu_yz',                                 &
              'wv', 'wv_xy', 'wv_xz', 'wv_yz',                                 &
              'wtheta', 'wtheta_xy', 'wtheta_xz', 'wtheta_yz',                 &
              'wq', 'wq_xy', 'wq_xz', 'wq_yz')

          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zu'
!
!--    s grid surface variables
       CASE ( 'theta_2m*_xy', 'wspeed_10m*' )

          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zu'
!
!--    u grid
       CASE ( 'uu', 'uu_xy', 'uu_xz', 'uu_yz' )

          grid_x = 'xu'
          grid_y = 'y'
          grid_z = 'zu'
!
!--    v grid
       CASE ( 'vv', 'vv_xy', 'vv_xz', 'vv_yz'  )

          grid_x = 'x'
          grid_y = 'yv'
          grid_z = 'zu'
!
!--    w grid
       CASE ( 'ww', 'ww_xy', 'ww_xz', 'ww_yz'  )

          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zw'

       CASE DEFAULT
          found  = .FALSE.
          grid_x = 'none'
          grid_y = 'none'
          grid_z = 'none'

    END SELECT


 END SUBROUTINE doq_define_netcdf_grid

!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining 2D output variables
!------------------------------------------------------------------------------!
 SUBROUTINE doq_output_2d( av, variable, found, grid,                          &
                           mode, local_pf, two_d, nzb_do, nzt_do, fill_value )


    IMPLICIT NONE

    CHARACTER (LEN=*) ::  grid     !<
    CHARACTER (LEN=*) ::  mode     !<
    CHARACTER (LEN=*) ::  variable !<

    INTEGER(iwp) ::  av       !< value indicating averaged or non-averaged output
    INTEGER(iwp) ::  flag_nr  !< number of the topography flag (0: scalar, 1: u, 2: v, 3: w)
    INTEGER(iwp) ::  i        !< grid index x-direction
    INTEGER(iwp) ::  j        !< grid index y-direction
    INTEGER(iwp) ::  k        !< grid index z-direction
    INTEGER(iwp) ::  nzb_do   !<
    INTEGER(iwp) ::  nzt_do   !<

    LOGICAL ::  found             !< true if variable is in list
    LOGICAL ::  resorted          !< true if array is resorted
    LOGICAL ::  two_d             !< flag parameter that indicates 2D variables (horizontal cross sections)

    REAL(wp) ::  fill_value       !< value for the _FillValue attribute

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf !<
    REAL(wp), DIMENSION(:,:,:), POINTER ::                 to_be_resorted  !< points to array which needs to be resorted for output

    flag_nr  = 0
    found    = .TRUE.
    resorted = .FALSE.
    two_d    = .FALSE.

    SELECT CASE ( TRIM( variable ) )

       CASE ( 'ti_xy', 'ti_xz', 'ti_yz' )
           IF ( av == 0 )  THEN
              to_be_resorted => ti
           ELSE
              IF ( .NOT. ALLOCATED( ti_av ) ) THEN
                 ALLOCATE( ti_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
                 ti_av = REAL( fill_value, KIND = wp )
              ENDIF
              to_be_resorted => ti_av
           ENDIF
           flag_nr = 0

           IF ( mode == 'xy' )  grid = 'zu'

       CASE ( 'uu_xy', 'uu_xz', 'uu_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => uu
          ELSE
             IF ( .NOT. ALLOCATED( uu_av ) ) THEN
                ALLOCATE( uu_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
                uu_av = REAL( fill_value, KIND = wp )
             ENDIF
             to_be_resorted => uu_av
          ENDIF
          flag_nr = 1

          IF ( mode == 'xy' )  grid = 'zu'

       CASE ( 'vv_xy', 'vv_xz', 'vv_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => vv
          ELSE
             IF ( .NOT. ALLOCATED( vv_av ) ) THEN
                ALLOCATE( vv_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
                vv_av = REAL( fill_value, KIND = wp )
             ENDIF
             to_be_resorted => vv_av
          ENDIF
          flag_nr = 2

          IF ( mode == 'xy' )  grid = 'zu'

       CASE ( 'ww_xy', 'ww_xz', 'ww_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => ww
          ELSE
             IF ( .NOT. ALLOCATED( ww_av ) ) THEN
                ALLOCATE( ww_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
                ww_av = REAL( fill_value, KIND = wp )
             ENDIF
             to_be_resorted => ww_av
          ENDIF
          flag_nr = 3

          IF ( mode == 'xy' )  grid = 'zw'

       CASE ( 'wu_xy', 'wu_xz', 'wu_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => wu
          ELSE
             IF ( .NOT. ALLOCATED( wu_av ) ) THEN
                ALLOCATE( wu_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
                wu_av = REAL( fill_value, KIND = wp )
             ENDIF
             to_be_resorted => wu_av
          ENDIF
          flag_nr = 0

          IF ( mode == 'xy' )  grid = 'zw'

       CASE ( 'wv_xy', 'wv_xz', 'wv_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => wv
          ELSE
             IF ( .NOT. ALLOCATED( wv_av ) ) THEN
                ALLOCATE( wv_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
                wv_av = REAL( fill_value, KIND = wp )
             ENDIF
             to_be_resorted => wv_av
          ENDIF
          flag_nr = 0

          IF ( mode == 'xy' )  grid = 'zw'

       CASE ( 'wtheta_xy', 'wtheta_xz', 'wtheta_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => wtheta
          ELSE
             IF ( .NOT. ALLOCATED( wtheta_av ) ) THEN
                ALLOCATE( wtheta_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
                wtheta_av = REAL( fill_value, KIND = wp )
             ENDIF
             to_be_resorted => wtheta_av
          ENDIF
          flag_nr = 0

          IF ( mode == 'xy' )  grid = 'zw'

       CASE ( 'wq_xy', 'wq_xz', 'wq_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => wq
          ELSE
             IF ( .NOT. ALLOCATED( wq_av ) ) THEN
                ALLOCATE( wq_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
                wq_av = REAL( fill_value, KIND = wp )
             ENDIF
             to_be_resorted => wq_av
          ENDIF
          flag_nr = 0

          IF ( mode == 'xy' )  grid = 'zw'

       CASE ( 'theta_2m*_xy' )        ! 2d-array
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = pt_2m(j,i)
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( pt_2m_av ) ) THEN
                ALLOCATE( pt_2m_av(nysg:nyng,nxlg:nxrg) )
                pt_2m_av = REAL( fill_value, KIND = wp )
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = pt_2m_av(j,i)
                ENDDO
             ENDDO
          ENDIF
          resorted = .TRUE.
          two_d    = .TRUE.
          grid     = 'zu1'

       CASE ( 'wspeed_10m*_xy' )        ! 2d-array
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = uv_10m(j,i)
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( uv_10m_av ) ) THEN
                ALLOCATE( uv_10m_av(nysg:nyng,nxlg:nxrg) )
                uv_10m_av = REAL( fill_value, KIND = wp )
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = uv_10m_av(j,i)
                ENDDO
             ENDDO
          ENDIF
          resorted = .TRUE.
          two_d    = .TRUE.
          grid     = 'zu1'

       CASE ( 'wspeed_xy', 'wspeed_xz', 'wspeed_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => wspeed
          ELSE
             IF ( .NOT. ALLOCATED( wspeed_av ) ) THEN
                ALLOCATE( wspeed_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
                wspeed_av = REAL( fill_value, KIND = wp )
             ENDIF
             to_be_resorted => wspeed_av
          ENDIF
          flag_nr = 0

          IF ( mode == 'xy' )  grid = 'zu'

       CASE ( 'wdir_xy', 'wdir_xz', 'wdir_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => wdir
          ELSE
             IF ( .NOT. ALLOCATED( wdir_av ) ) THEN
                ALLOCATE( wdir_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
                wdir_av = REAL( fill_value, KIND = wp )
             ENDIF
             to_be_resorted => wdir_av
          ENDIF
          flag_nr = 0

          IF ( mode == 'xy' )  grid = 'zu'

       CASE DEFAULT
          found = .FALSE.
          grid  = 'none'

    END SELECT

    IF ( found  .AND.  .NOT. resorted )  THEN
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb_do, nzt_do
                local_pf(i,j,k) = MERGE( to_be_resorted(k,j,i),                &
                                     REAL( fill_value, KIND = wp ),            &
                                     BTEST( wall_flags_total_0(k,j,i), flag_nr ) )
             ENDDO
          ENDDO
       ENDDO
    ENDIF

 END SUBROUTINE doq_output_2d


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining 3D output variables
!------------------------------------------------------------------------------!
 SUBROUTINE doq_output_3d( av, variable, found, local_pf, fill_value, nzb_do,  &
                           nzt_do )

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  variable !<

    INTEGER(iwp) ::  av       !< index indicating averaged or instantaneous output
    INTEGER(iwp) ::  flag_nr  !< number of the topography flag (0: scalar, 1: u, 2: v, 3: w)
    INTEGER(iwp) ::  i        !< index variable along x-direction
    INTEGER(iwp) ::  j        !< index variable along y-direction
    INTEGER(iwp) ::  k        !< index variable along z-direction
    INTEGER(iwp) ::  nzb_do   !< lower limit of the data output (usually 0)
    INTEGER(iwp) ::  nzt_do   !< vertical upper limit of the data output (usually nz_do3d)

    LOGICAL ::  found             !< true if variable is in list
    LOGICAL ::  resorted          !< true if array is resorted

    REAL(wp) ::  fill_value       !< value for the _FillValue attribute

    REAL(sp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf        !<
    REAL(wp), DIMENSION(:,:,:), POINTER ::                 to_be_resorted  !< points to array which needs to be resorted for output

    flag_nr  = 0
    found    = .TRUE.
    resorted = .FALSE.

    SELECT CASE ( TRIM( variable ) )

       CASE ( 'ti' )
          IF ( av == 0 )  THEN
             to_be_resorted => ti
          ELSE
             IF ( .NOT. ALLOCATED( ti_av ) ) THEN
                ALLOCATE( ti_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
                ti_av = REAL( fill_value, KIND = wp )
             ENDIF
             to_be_resorted => ti_av
          ENDIF
          flag_nr = 0

       CASE ( 'uu' )
          IF ( av == 0 )  THEN
             to_be_resorted => uu
          ELSE
             IF ( .NOT. ALLOCATED( uu_av ) ) THEN
                ALLOCATE( uu_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
                uu_av = REAL( fill_value, KIND = wp )
             ENDIF
             to_be_resorted => uu_av
          ENDIF
          flag_nr = 1

       CASE ( 'vv' )
          IF ( av == 0 )  THEN
             to_be_resorted => vv
          ELSE
             IF ( .NOT. ALLOCATED( vv_av ) ) THEN
                ALLOCATE( vv_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
                vv_av = REAL( fill_value, KIND = wp )
             ENDIF
             to_be_resorted => vv_av
          ENDIF
          flag_nr = 2

       CASE ( 'ww' )
          IF ( av == 0 )  THEN
             to_be_resorted => ww
          ELSE
             IF ( .NOT. ALLOCATED( ww_av ) ) THEN
                ALLOCATE( ww_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
                ww_av = REAL( fill_value, KIND = wp )
             ENDIF
             to_be_resorted => ww_av
          ENDIF
          flag_nr = 3

       CASE ( 'wu' )
          IF ( av == 0 )  THEN
             to_be_resorted => wu
          ELSE
             IF ( .NOT. ALLOCATED( wu_av ) ) THEN
                ALLOCATE( wu_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
                wu_av = REAL( fill_value, KIND = wp )
             ENDIF
             to_be_resorted => wu_av
          ENDIF
          flag_nr = 0

       CASE ( 'wv' )
          IF ( av == 0 )  THEN
             to_be_resorted => wv
          ELSE
             IF ( .NOT. ALLOCATED( wv_av ) ) THEN
                ALLOCATE( wv_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
                wv_av = REAL( fill_value, KIND = wp )
             ENDIF
             to_be_resorted => wv_av
          ENDIF
          flag_nr = 0

       CASE ( 'wtheta' )
          IF ( av == 0 )  THEN
             to_be_resorted => wtheta
          ELSE
             IF ( .NOT. ALLOCATED( wtheta_av ) ) THEN
                ALLOCATE( wtheta_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
                wtheta_av = REAL( fill_value, KIND = wp )
             ENDIF
             to_be_resorted => wtheta_av
          ENDIF
          flag_nr = 0

       CASE ( 'wq' )
          IF ( av == 0 )  THEN
             to_be_resorted => wq
          ELSE
             IF ( .NOT. ALLOCATED( wq_av ) ) THEN
                ALLOCATE( wq_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
                wq_av = REAL( fill_value, KIND = wp )
             ENDIF
             to_be_resorted => wq_av
          ENDIF
          flag_nr = 0

       CASE ( 'wspeed' )
          IF ( av == 0 )  THEN
             to_be_resorted => wspeed
          ELSE
             IF ( .NOT. ALLOCATED( wspeed_av ) ) THEN
                ALLOCATE( wspeed_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
                wspeed_av = REAL( fill_value, KIND = wp )
             ENDIF
             to_be_resorted => wspeed_av
          ENDIF
          flag_nr = 0

       CASE ( 'wdir' )
          IF ( av == 0 )  THEN
             to_be_resorted => wdir
          ELSE
             IF ( .NOT. ALLOCATED( wdir_av ) ) THEN
                ALLOCATE( wdir_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
                wdir_av = REAL( fill_value, KIND = wp )
             ENDIF
             to_be_resorted => wdir_av
          ENDIF
          flag_nr = 0

       CASE DEFAULT
          found = .FALSE.

    END SELECT

    IF ( found  .AND.  .NOT. resorted )  THEN
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb_do, nzt_do
                local_pf(i,j,k) = MERGE( to_be_resorted(k,j,i),                &
                                     REAL( fill_value, KIND = wp ),            &
                                     BTEST( wall_flags_total_0(k,j,i), flag_nr ) )
             ENDDO
          ENDDO
       ENDDO
    ENDIF

 END SUBROUTINE doq_output_3d

! Description:
! ------------
!> Resorts the user-defined output quantity with indices (k,j,i) to a
!> temporary array with indices (i,j,k) for masked data output.
!------------------------------------------------------------------------------!
 SUBROUTINE doq_output_mask( av, variable, found, local_pf, mid )

    USE control_parameters

    USE indices

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  variable   !<
    CHARACTER (LEN=5) ::  grid       !< flag to distinquish between staggered grids

    INTEGER(iwp) ::  av              !< index indicating averaged or instantaneous output
    INTEGER(iwp) ::  flag_nr         !< number of the topography flag (0: scalar, 1: u, 2: v, 3: w)
    INTEGER(iwp) ::  i               !< index variable along x-direction
    INTEGER(iwp) ::  j               !< index variable along y-direction
    INTEGER(iwp) ::  k               !< index variable along z-direction
    INTEGER(iwp) ::  im              !< loop index for masked variables
    INTEGER(iwp) ::  jm              !< loop index for masked variables
    INTEGER(iwp) ::  kk              !< masked output index variable along z-direction
    INTEGER(iwp) ::  mid             !< masked output running index
    INTEGER(iwp) ::  ktt             !< k index of highest horizontal surface

    LOGICAL      ::  found           !< true if variable is in list
    LOGICAL      ::  resorted        !< true if array is resorted

    REAL(wp),                                                                  &
       DIMENSION(mask_size_l(mid,1),mask_size_l(mid,2),mask_size_l(mid,3)) ::  &
          local_pf   !<
    REAL(wp), DIMENSION(:,:,:), POINTER ::  to_be_resorted  !< points to array which needs to be resorted for output

    REAL(wp), PARAMETER   ::  fill_value = -9999.0_wp       !< value for the _FillValue attribute

    flag_nr  = 0
    found    = .TRUE.
    resorted = .FALSE.
    grid     = 's'

    SELECT CASE ( TRIM( variable ) )

       CASE ( 'ti' )
          IF ( av == 0 )  THEN
             to_be_resorted => ti
          ELSE
             to_be_resorted => ti_av
          ENDIF
          grid = 's'
          flag_nr = 0

       CASE ( 'uu' )
          IF ( av == 0 )  THEN
             to_be_resorted => uu
          ELSE
             to_be_resorted => uu_av
          ENDIF
          grid = 'u'
          flag_nr = 1

       CASE ( 'vv' )
          IF ( av == 0 )  THEN
             to_be_resorted => vv
          ELSE
             to_be_resorted => vv_av
          ENDIF
          grid = 'v'
          flag_nr = 2

       CASE ( 'ww' )
          IF ( av == 0 )  THEN
             to_be_resorted => ww
          ELSE
             to_be_resorted => ww_av
          ENDIF
          grid = 'w'
          flag_nr = 3

       CASE ( 'wu' )
          IF ( av == 0 )  THEN
             to_be_resorted => wu
          ELSE
             to_be_resorted => wu_av
          ENDIF
          grid = 's'
          flag_nr = 0

       CASE ( 'wv' )
          IF ( av == 0 )  THEN
             to_be_resorted => wv
          ELSE
             to_be_resorted => wv_av
          ENDIF
          grid = 's'
          flag_nr = 0

       CASE ( 'wtheta' )
          IF ( av == 0 )  THEN
             to_be_resorted => wtheta
          ELSE
             to_be_resorted => wtheta_av
          ENDIF
          grid = 's'
          flag_nr = 0

       CASE ( 'wq' )
          IF ( av == 0 )  THEN
             to_be_resorted => wq
          ELSE
             to_be_resorted => wq_av
          ENDIF
          grid = 's'
          flag_nr = 0

       CASE ( 'wspeed' )
          IF ( av == 0 )  THEN
             to_be_resorted => wspeed
          ELSE
             to_be_resorted => wspeed_av
          ENDIF
          grid = 's'
          flag_nr = 0

       CASE ( 'wdir' )
          IF ( av == 0 )  THEN
             to_be_resorted => wdir
          ELSE
             to_be_resorted => wdir_av
          ENDIF
          grid = 's'
          flag_nr = 0

       CASE DEFAULT
          found = .FALSE.

    END SELECT

    IF ( found  .AND.  .NOT. resorted )  THEN
       IF ( .NOT. mask_surface(mid) )  THEN
!
!--       Default masked output
          DO  i = 1, mask_size_l(mid,1)
             DO  j = 1, mask_size_l(mid,2)
                DO  k = 1, mask_size_l(mid,3)
                   local_pf(i,j,k) = MERGE( to_be_resorted(mask_k(mid,k),  &
                                                           mask_j(mid,j),  &
                                                           mask_i(mid,i)), &
                                            REAL( fill_value, KIND = wp ), &
                                            BTEST( wall_flags_total_0(     &
                                                           mask_k(mid,k),  &
                                                           mask_j(mid,j),  &
                                                           mask_i(mid,i)), &
                                                   flag_nr ) )
                ENDDO
             ENDDO
          ENDDO

       ELSE
!
!--       Terrain-following masked output
          DO  i = 1, mask_size_l(mid,1)
             DO  j = 1, mask_size_l(mid,2)
!
!--             Get k index of the highest terraing surface
                im = mask_i(mid,i)
                jm = mask_j(mid,j)
                ktt = MINLOC( MERGE( 1, 0, BTEST( wall_flags_total_0(:,jm,im), 5 )), &
                              DIM = 1 ) - 1
                DO  k = 1, mask_size_l(mid,3)
                   kk = MIN( ktt+mask_k(mid,k), nzt+1 )
!
!--                Set value if not in building
                   IF ( BTEST( wall_flags_total_0(kk,jm,im), 6 ) )  THEN
                      local_pf(i,j,k) = fill_value
                   ELSE
                      local_pf(i,j,k) = to_be_resorted(kk,jm,im)
                   ENDIF
                ENDDO
             ENDDO
          ENDDO

       ENDIF
    ENDIF

 END SUBROUTINE doq_output_mask

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate required arrays
!------------------------------------------------------------------------------!
 SUBROUTINE doq_init

    IMPLICIT NONE

    INTEGER(iwp) ::  ivar   !< loop index over all 2d/3d/mask output quantities

!
!-- Next line is to avoid compiler warnings about unused variables
    IF ( timestep_number_at_prev_calc == 0 )  CONTINUE
!
!-- Preparatory steps and initialization of output arrays
    IF ( .NOT.  prepared_diagnostic_output_quantities )  CALL doq_prepare

    initialized_diagnostic_output_quantities = .FALSE.

    ivar = 1

    DO  WHILE ( ivar <= SIZE( do_all ) )

       SELECT CASE ( TRIM( do_all(ivar) ) )
!
!--       Allocate array for 'turbulence intensity'
          CASE ( 'ti' )
             IF ( .NOT. ALLOCATED( ti ) )  THEN
                ALLOCATE( ti(nzb:nzt+1,nys:nyn,nxl:nxr) )
                ti = 0.0_wp
             ENDIF
!
!--       Allocate array for uu
          CASE ( 'uu' )
             IF ( .NOT. ALLOCATED( uu ) )  THEN
                ALLOCATE( uu(nzb:nzt+1,nys:nyn,nxl:nxr) )
                uu = 0.0_wp
             ENDIF
!
!--       Allocate array for vv
          CASE ( 'vv' )
             IF ( .NOT. ALLOCATED( vv ) )  THEN
                ALLOCATE( vv(nzb:nzt+1,nys:nyn,nxl:nxr) )
                vv = 0.0_wp
             ENDIF
!
!--       Allocate array for ww
          CASE ( 'ww' )
             IF ( .NOT. ALLOCATED( ww ) )  THEN
                ALLOCATE( ww(nzb:nzt+1,nys:nyn,nxl:nxr) )
                ww = 0.0_wp
             ENDIF
!
!--       Allocate array for wu
          CASE ( 'wu' )
             IF ( .NOT. ALLOCATED( wu ) )  THEN
                ALLOCATE( wu(nzb:nzt+1,nys:nyn,nxl:nxr) )
                wu = 0.0_wp
             ENDIF
!
!--       Allocate array for wv
          CASE ( 'wv' )
             IF ( .NOT. ALLOCATED( wv ) )  THEN
                ALLOCATE( wv(nzb:nzt+1,nys:nyn,nxl:nxr) )
                wv = 0.0_wp
             ENDIF
!
!--       Allocate array for wtheta
          CASE ( 'wtheta' )
             IF ( .NOT. ALLOCATED( wtheta ) )  THEN
                ALLOCATE( wtheta(nzb:nzt+1,nys:nyn,nxl:nxr) )
                wtheta = 0.0_wp
             ENDIF
!
!--       Allocate array for wq
          CASE ( 'wq' )
             IF ( .NOT. ALLOCATED( wq ) )  THEN
                ALLOCATE( wq(nzb:nzt+1,nys:nyn,nxl:nxr) )
                wq = 0.0_wp
             ENDIF
!
!--       Allocate array for 2-m potential temperature
          CASE ( 'theta_2m*' )
             IF ( .NOT. ALLOCATED( pt_2m ) )  THEN
                ALLOCATE( pt_2m(nys:nyn,nxl:nxr) )
                pt_2m = 0.0_wp
             ENDIF
!
!--       Allocate array for 10-m wind speed
          CASE ( 'wspeed_10m*' )
             IF ( .NOT. ALLOCATED( uv_10m ) )  THEN
                ALLOCATE( uv_10m(nys:nyn,nxl:nxr) )
                uv_10m = 0.0_wp
             ENDIF
!
!--       Allocate array for wspeed
          CASE ( 'wspeed' )
             IF ( .NOT. ALLOCATED( wspeed ) )  THEN
                ALLOCATE( wspeed(nzb:nzt+1,nys:nyn,nxl:nxr) )
                wspeed = 0.0_wp
             ENDIF

!
!--       Allocate array for wdir
          CASE ( 'wdir' )
             IF ( .NOT. ALLOCATED( u_center ) )  THEN
                ALLOCATE( u_center(nzb:nzt+1,nys:nyn,nxl:nxr) )
                u_center = 0.0_wp
             ENDIF
             IF ( .NOT. ALLOCATED( v_center ) )  THEN
                ALLOCATE( v_center(nzb:nzt+1,nys:nyn,nxl:nxr) )
                v_center = 0.0_wp
             ENDIF
             IF ( .NOT. ALLOCATED( wdir ) )  THEN
                ALLOCATE( wdir(nzb:nzt+1,nys:nyn,nxl:nxr) )
                wdir = 0.0_wp
             ENDIF

       END SELECT

       ivar = ivar + 1
    ENDDO

    initialized_diagnostic_output_quantities = .TRUE.

 END SUBROUTINE doq_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate diagnostic quantities
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE doq_calculate

    IMPLICIT NONE

    INTEGER(iwp) ::  i          !< grid index x-dimension
    INTEGER(iwp) ::  j          !< grid index y-dimension
    INTEGER(iwp) ::  k          !< grid index z-dimension
    INTEGER(iwp) ::  ivar       !< loop index over all 2d/3d/mask output quantities

    TYPE(surf_type), POINTER ::  surf     !< surf-type array, used to generalize subroutines


!     CALL cpu_log( log_point(41), 'calculate_quantities', 'start' )

!
!-- Save timestep number to check in time_integration if doq_calculate
!-- has been called already, since the CALL occurs at two locations, but the calculations need to be
!-- done only once per timestep.
    timestep_number_at_prev_calc = current_timestep_number

    ivar = 1

    DO  WHILE ( ivar <= SIZE( do_all ) )

       SELECT CASE ( TRIM( do_all(ivar) ) )
!
!--       Calculate 'turbulence intensity' from rot[(u,v,w)] at scalar grid point
          CASE ( 'ti' )
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt
                      ti(k,j,i) = 0.25_wp * SQRT(                              &
                        (   (   w(k,j+1,i) + w(k-1,j+1,i)                      &
                              - w(k,j-1,i) - w(k-1,j-1,i) ) * ddy              &
                          - (   v(k+1,j,i) + v(k+1,j+1,i)                      &
                              - v(k-1,j,i) - v(k-1,j+1,i) ) * ddzu(k) )**2     &
                      + (   (   u(k+1,j,i) + u(k+1,j,i+1)                      &
                              - u(k-1,j,i) - u(k-1,j,i+1) ) * ddzu(k)          &
                          - (   w(k,j,i+1) + w(k-1,j,i+1)                      &
                              - w(k,j,i-1) - w(k-1,j,i-1) ) * ddx     )**2     &
                      + (   (   v(k,j,i+1) + v(k,j+1,i+1)                      &
                              - v(k,j,i-1) - v(k,j+1,i-1) ) * ddx              &
                          - (   u(k,j+1,i) + u(k,j+1,i+1)                      &
                              - u(k,j-1,i) - u(k,j-1,i+1) ) * ddy     )**2  )  &
                       * MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0) )
                   ENDDO
                ENDDO
             ENDDO
!
!--       uu
          CASE ( 'uu' )
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt
                      uu(k,j,i) = u(k,j,i) * u(k,j,i)                          &
                       * MERGE( 1.0_wp, 0.0_wp,                                &
                       BTEST( wall_flags_total_0(k,j,i), 1) )
                   ENDDO
                ENDDO
             ENDDO
!
!--       vv
          CASE ( 'vv' )
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt
                      vv(k,j,i) = v(k,j,i) * v(k,j,i)                          &
                       * MERGE( 1.0_wp, 0.0_wp,                                &
                       BTEST( wall_flags_total_0(k,j,i), 2) )
                   ENDDO
                ENDDO
             ENDDO
!
!--       ww
          CASE ( 'ww' )
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt-1
                      ww(k,j,i) = w(k,j,i) * w(k,j,i)                          &
                       * MERGE( 1.0_wp, 0.0_wp,                                &
                       BTEST( wall_flags_total_0(k,j,i), 3) )
                   ENDDO
                ENDDO
             ENDDO
!
!--       wu
          CASE ( 'wu' )
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt-1
                      wu(k,j,i) = w(k,j,i) * u(k,j,i)                          &
                       * MERGE( 1.0_wp, 0.0_wp,                                &
                       BTEST( wall_flags_total_0(k,j,i), 0) )
                   ENDDO
                ENDDO
             ENDDO
!
!--       wv
          CASE ( 'wv' )
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt-1
                      wv(k,j,i) = w(k,j,i) * v(k,j,i)                          &
                       * MERGE( 1.0_wp, 0.0_wp,                                &
                       BTEST( wall_flags_total_0(k,j,i), 0) )
                   ENDDO
                ENDDO
             ENDDO
!
!--       wtheta
          CASE ( 'wtheta' )
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt-1
                      wtheta(k,j,i) = w(k,j,i) * pt(k,j,i)                     &
                       * MERGE( 1.0_wp, 0.0_wp,                                &
                       BTEST( wall_flags_total_0(k,j,i), 0) )
                   ENDDO
                ENDDO
             ENDDO
!
!--       wq
          CASE ( 'wq' )
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt-1
                      wq(k,j,i) = w(k,j,i) * q(k,j,i)                          &
                       * MERGE( 1.0_wp, 0.0_wp,                                &
                       BTEST( wall_flags_total_0(k,j,i), 0) )
                   ENDDO
                ENDDO
             ENDDO
!
!--       2-m potential temperature
          CASE ( 'theta_2m*' )
!
!--          2-m potential temperature is caluclated from surface arrays. In
!--          case the 2m level is below the first grid point, MOST is applied,
!--          else, linear interpolation between two vertical grid levels is
!--          applied. To access all surfaces, iterate over all horizontally-
!--          upward facing surface types.
             surf => surf_def_h(0)
             CALL calc_pt_2m
             surf => surf_lsm_h
             CALL calc_pt_2m
             surf => surf_usm_h
             CALL calc_pt_2m
!
!--       10-m wind speed
          CASE ( 'wspeed_10m*' )
!
!--          10-m wind speed is caluclated from surface arrays. In
!--          case the 10m level is below the first grid point, MOST is applied,
!--          else, linear interpolation between two vertical grid levels is
!--          applied. To access all surfaces, iterate over all horizontally-
!--          upward facing surface types.
             surf => surf_def_h(0)
             CALL calc_wind_10m
             surf => surf_lsm_h
             CALL calc_wind_10m
             surf => surf_usm_h
             CALL calc_wind_10m
!
!--       horizontal wind speed
          CASE ( 'wspeed' )
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb, nzt+1
                      wspeed(k,j,i) = SQRT( ( 0.5_wp * ( u(k,j,i) + u(k,j,i+1) ) )**2             &
                                          + ( 0.5_wp * ( v(k,j,i) + v(k,j+1,i) ) )**2 )           &
                                    * MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0) )
                   ENDDO
                ENDDO
             ENDDO

!
!--       horizontal wind direction
          CASE ( 'wdir' )
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb, nzt+1
                      u_center(k,j,i) = 0.5_wp * ( u(k,j,i) + u(k,j,i+1) )
                      v_center(k,j,i) = 0.5_wp * ( v(k,j,i) + v(k,j+1,i) )

                      wdir(k,j,i) = ATAN2( u_center(k,j,i), v_center(k,j,i) ) &
                                  / pi * 180.0_wp + 180.0_wp
                   ENDDO
                ENDDO
             ENDDO

       END SELECT

       ivar = ivar + 1
    ENDDO

!     CALL cpu_log( log_point(41), 'calculate_quantities', 'stop' )

!
!-- The following block contains subroutines to calculate diagnostic
!-- surface quantities.
    CONTAINS
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculation of 2-m potential temperature.
!------------------------------------------------------------------------------!
       SUBROUTINE calc_pt_2m

          USE surface_layer_fluxes_mod,                                        &
              ONLY:  psi_h

          IMPLICIT NONE

          INTEGER(iwp) ::  kk     !< running index along the z-dimension
          INTEGER(iwp) ::  m      !< running index for surface elements

          DO  m = 1, surf%ns

             i = surf%i(m)
             j = surf%j(m)
             k = surf%k(m)
!
!--          If 2-m level is below the first grid level, MOST is
!--          used for calculation of 2-m temperature.
             IF ( surf%z_mo(m) > 2.0_wp )  THEN
                pt_2m(j,i) = surf%pt_surface(m) + surf%ts(m) / kappa           &
                                * ( LOG( 2.0_wp /  surf%z0h(m) )               &
                                  - psi_h( 2.0_wp / surf%ol(m) )               &
                                  + psi_h( surf%z0h(m) / surf%ol(m) ) )
!
!--          If 2-m level is above the first grid level, 2-m temperature
!--          is linearly interpolated between the two nearest vertical grid
!--          levels. Note, since 2-m temperature is only computed for
!--          horizontal upward-facing surfaces, only a vertical
!--          interpolation is necessary.
             ELSE
!
!--             zw(k-1) defines the height of the surface.
                kk = k
                DO WHILE ( zu(kk) - zw(k-1) < 2.0_wp  .AND.  kk <= nzt )
                   kk = kk + 1
                ENDDO
!
!--             kk defines the index of the first grid level >= 2m.
                pt_2m(j,i) = pt(kk-1,j,i) +                                    &
                              ( zw(k-1) + 2.0_wp - zu(kk-1)     ) *            &
                              ( pt(kk,j,i)       - pt(kk-1,j,i) ) /            &
                              ( zu(kk)           - zu(kk-1)     )
             ENDIF

          ENDDO

       END SUBROUTINE calc_pt_2m

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculation of 10-m wind speed.
!------------------------------------------------------------------------------!
       SUBROUTINE calc_wind_10m

          USE surface_layer_fluxes_mod,                                        &
              ONLY:  psi_m

          IMPLICIT NONE

          INTEGER(iwp) ::  kk     !< running index along the z-dimension
          INTEGER(iwp) ::  m      !< running index for surface elements

          REAL(wp) ::  uv_l !< wind speed at lower grid point
          REAL(wp) ::  uv_u !< wind speed at upper grid point

          DO  m = 1, surf%ns

             i = surf%i(m)
             j = surf%j(m)
             k = surf%k(m)
!
!--          If 10-m level is below the first grid level, MOST is
!--          used for calculation of 10-m temperature.
             IF ( surf%z_mo(m) > 10.0_wp )  THEN
                uv_10m(j,i) = surf%us(m) / kappa                               &
                          * ( LOG( 10.0_wp /  surf%z0(m) )                     &
                              - psi_m( 10.0_wp    / surf%ol(m) )               &
                              + psi_m( surf%z0(m) / surf%ol(m) ) )
!
!--          If 10-m level is above the first grid level, 10-m wind speed
!--          is linearly interpolated between the two nearest vertical grid
!--          levels. Note, since 10-m temperature is only computed for
!--          horizontal upward-facing surfaces, only a vertical
!--          interpolation is necessary.
             ELSE
!
!--             zw(k-1) defines the height of the surface.
                kk = k
                DO WHILE ( zu(kk) - zw(k-1) < 10.0_wp  .AND.  kk <= nzt )
                   kk = kk + 1
                ENDDO
!
!--             kk defines the index of the first grid level >= 10m.
                uv_l = SQRT( ( 0.5_wp * ( u(kk-1,j,i) + u(kk-1,j,i+1) ) )**2   &
                           + ( 0.5_wp * ( v(kk-1,j,i) + v(kk-1,j+1,i) ) )**2 )

                uv_u = SQRT( ( 0.5_wp * ( u(kk,j,i)   + u(kk,j,i+1)   ) )**2   &
                           + ( 0.5_wp * ( v(kk,j,i)   + v(kk,j+1,i)   ) )**2 )

                uv_10m(j,i) = uv_l + ( zw(k-1) + 10.0_wp - zu(kk-1) ) *        &
                                     ( uv_u              - uv_l     ) /        &
                                     ( zu(kk)            - zu(kk-1) )

             ENDIF

          ENDDO

       END SUBROUTINE calc_wind_10m

 END SUBROUTINE doq_calculate


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Preparation of the diagnostic output, counting of the module-specific
!> output quantities and gathering of the output names.
!------------------------------------------------------------------------------!
 SUBROUTINE doq_prepare


    USE control_parameters,                                                    &
        ONLY:  do2d, do3d, domask, masks

    IMPLICIT NONE

    CHARACTER (LEN=varnamelength), DIMENSION(0:1,500) ::  do2d_var = ' '  !<
                                                          !< label array for 2d output quantities

    INTEGER(iwp) ::  av         !< index defining type of output, av=0 instantaneous, av=1 averaged
    INTEGER(iwp) ::  ivar       !< loop index
    INTEGER(iwp) ::  ivar_all   !< loop index
    INTEGER(iwp) ::  l          !< index for cutting string
    INTEGER(iwp) ::  mid          !< masked output running index

    prepared_diagnostic_output_quantities = .FALSE.

    ivar     = 1
    ivar_all = 1

    DO  av = 0, 1
!
!--    Remove _xy, _xz, or _yz from string
       l = MAX( 3, LEN_TRIM( do2d(av,ivar) ) )
       do2d_var(av,ivar)(1:l-3) = do2d(av,ivar)(1:l-3)
!
!--    Gather 2d output quantity names.
!--    Check for double occurrence of output quantity, e.g. by _xy,
!--    _yz, _xz.
       DO  WHILE ( do2d_var(av,ivar)(1:1) /= ' ' )
          IF ( .NOT.  ANY( do_all == do2d_var(av,ivar) ) )  THEN
             do_all(ivar_all) = do2d_var(av,ivar)
          ENDIF
          ivar = ivar + 1
          ivar_all = ivar_all + 1
          l = MAX( 3, LEN_TRIM( do2d(av,ivar) ) )
          do2d_var(av,ivar)(1:l-3) = do2d(av,ivar)(1:l-3)
       ENDDO

       ivar = 1
!
!--    Gather 3d output quantity names
       DO  WHILE ( do3d(av,ivar)(1:1) /= ' ' )
          do_all(ivar_all) = do3d(av,ivar)
          ivar = ivar + 1
          ivar_all = ivar_all + 1
       ENDDO

       ivar = 1
!
!--    Gather masked output quantity names. Also check for double output
!--    e.g. by different masks.
       DO  mid = 1, masks
          DO  WHILE ( domask(mid,av,ivar)(1:1) /= ' ' )
             IF ( .NOT.  ANY( do_all == domask(mid,av,ivar) ) )  THEN
                do_all(ivar_all) = domask(mid,av,ivar)
             ENDIF

             ivar = ivar + 1
             ivar_all = ivar_all + 1
          ENDDO
          ivar = 1
       ENDDO

    ENDDO

    prepared_diagnostic_output_quantities = .TRUE.

 END SUBROUTINE doq_prepare

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine reads local (subdomain) restart data
!> Note: With the current structure reading of non-standard array is not
!> possible
!------------------------------------------------------------------------------!
!  SUBROUTINE doq_rrd_local( k, nxlf, nxlc, nxl_on_file, nxrf, nxrc,             &
!                            nxr_on_file, nynf, nync, nyn_on_file, nysf,         &
!                            nysc, nys_on_file, tmp_3d_non_standard, found )
!
!
!     USE control_parameters
!
!     USE indices
!
!     USE kinds
!
!     USE pegrid
!
!
!     IMPLICIT NONE
!
!     INTEGER(iwp) ::  k               !<
!     INTEGER(iwp) ::  nxlc            !<
!     INTEGER(iwp) ::  nxlf            !<
!     INTEGER(iwp) ::  nxl_on_file     !<
!     INTEGER(iwp) ::  nxrc            !<
!     INTEGER(iwp) ::  nxrf            !<
!     INTEGER(iwp) ::  nxr_on_file     !<
!     INTEGER(iwp) ::  nync            !<
!     INTEGER(iwp) ::  nynf            !<
!     INTEGER(iwp) ::  nyn_on_file     !<
!     INTEGER(iwp) ::  nysc            !<
!     INTEGER(iwp) ::  nysf            !<
!     INTEGER(iwp) ::  nys_on_file     !<
!
!     LOGICAL, INTENT(OUT)  :: found
!
!     REAL(wp), DIMENSION(:,:,:), ALLOCATABLE  ::  tmp_3d_non_standard !< temporary array for storing 3D data with non standard dimensions
! !
! !-- If temporary non-standard array for reading is already allocated,
! !-- deallocate it.
!     IF ( ALLOCATED( tmp_3d_non_standard ) )  DEALLOCATE( tmp_3d_non_standard )
!
!     found = .TRUE.
!
!     SELECT CASE ( restart_string(1:length) )
!
!        CASE ( 'ti_av' )
!           IF ( .NOT. ALLOCATED( ti_av ) )  THEN
!              ALLOCATE( ti_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
!           ENDIF
!           IF ( k == 1 )  THEN
!              ALLOCATE( tmp_3d_non_standard(nzb:nzt+1,nys_on_file:nyn_on_file,  &
!                                            nxl_on_file:nxr_on_file) )
!              READ ( 13 )  tmp_3d_non_standard
!           ENDIF
!           ti_av(:,nysc:nync,nxlc:nxrc) = tmp_3d_non_standard(:,nysf:nynf,nxlf:nxrf)
!
!        CASE ( 'uu_av' )
!           IF ( .NOT. ALLOCATED( uu_av ) )  THEN
!              ALLOCATE( uu_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
!           ENDIF
!           IF ( k == 1 )  THEN
!              ALLOCATE( tmp_3d_non_standard(nzb:nzt+1,nys_on_file:nyn_on_file,  &
!                                            nxl_on_file:nxr_on_file) )
!              READ ( 13 )  tmp_3d_non_standard
!           ENDIF
!           uu_av(:,nysc:nync,nxlc:nxrc) = tmp_3d_non_standard(:,nysf:nynf,nxlf:nxrf)
!
!        CASE ( 'vv_av' )
!           IF ( .NOT. ALLOCATED( vv_av ) )  THEN
!              ALLOCATE( vv_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
!           ENDIF
!           IF ( k == 1 )  THEN
!              ALLOCATE( tmp_3d_non_standard(nzb:nzt+1,nys_on_file:nyn_on_file,  &
!                                            nxl_on_file:nxr_on_file) )
!              READ ( 13 )  tmp_3d_non_standard
!           ENDIF
!           vv_av(:,nysc:nync,nxlc:nxrc) = tmp_3d_non_standard(:,nysf:nynf,nxlf:nxrf)
!
!        CASE ( 'ww_av' )
!           IF ( .NOT. ALLOCATED( ww_av ) )  THEN
!              ALLOCATE( ww_av(nzb:nzt+1,nys:nyn,nxl:nxr) )
!           ENDIF
!           IF ( k == 1 )  THEN
!              ALLOCATE( tmp_3d_non_standard(nzb:nzt+1,nys_on_file:nyn_on_file,  &
!                                            nxl_on_file:nxr_on_file) )
!              READ ( 13 )  tmp_3d_non_standard
!           ENDIF
!           ww_av(:,nysc:nync,nxlc:nxrc) = tmp_3d_non_standard(:,nysf:nynf,nxlf:nxrf)
!
!
!        CASE DEFAULT
!
!           found = .FALSE.
!
!     END SELECT
!
!  END SUBROUTINE doq_rrd_local

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine writes local (subdomain) restart data
!------------------------------------------------------------------------------!
 SUBROUTINE doq_wrd_local


    IMPLICIT NONE

    IF ( ALLOCATED( pt_2m_av ) )  THEN
       CALL wrd_write_string( 'pt_2m_av' )
       WRITE ( 14 )  pt_2m_av
    ENDIF

    IF ( ALLOCATED( ti_av ) )  THEN
       CALL wrd_write_string( 'ti_av' )
       WRITE ( 14 )  ti_av
    ENDIF

    IF ( ALLOCATED( uu_av ) )  THEN
       CALL wrd_write_string( 'uu_av' )
       WRITE ( 14 )  uu_av
    ENDIF

    IF ( ALLOCATED( uv_10m_av ) )  THEN
       CALL wrd_write_string( 'uv_10m_av' )
       WRITE ( 14 )  uv_10m_av
    ENDIF

    IF ( ALLOCATED( vv_av ) )  THEN
       CALL wrd_write_string( 'vv_av' )
       WRITE ( 14 )  vv_av
    ENDIF

    IF ( ALLOCATED( ww_av ) )  THEN
       CALL wrd_write_string( 'ww_av' )
       WRITE ( 14 )  ww_av
    ENDIF

    IF ( ALLOCATED( wu_av ) )  THEN
       CALL wrd_write_string( 'wu_av' )
       WRITE ( 14 )  wu_av
    ENDIF

    IF ( ALLOCATED( wv_av ) )  THEN
       CALL wrd_write_string( 'wv_av' )
       WRITE ( 14 )  wv_av
    ENDIF

    IF ( ALLOCATED( wtheta_av ) )  THEN
       CALL wrd_write_string( 'wtheta_av' )
       WRITE ( 14 )  wtheta_av
    ENDIF

    IF ( ALLOCATED( wq_av ) )  THEN
       CALL wrd_write_string( 'wq_av' )
       WRITE ( 14 )  wq_av
    ENDIF

    IF ( ALLOCATED( wspeed_av ) )  THEN
       CALL wrd_write_string( 'wspeed_av' )
       WRITE ( 14 )  wspeed_av
    ENDIF

    IF ( ALLOCATED( u_center_av ) )  THEN
       CALL wrd_write_string( 'u_center_av' )
       WRITE ( 14 )  u_center_av
    ENDIF

    IF ( ALLOCATED( v_center_av ) )  THEN
       CALL wrd_write_string( 'v_center_av' )
       WRITE ( 14 )  v_center_av
    ENDIF

 END SUBROUTINE doq_wrd_local



 END MODULE diagnostic_output_quantities_mod
