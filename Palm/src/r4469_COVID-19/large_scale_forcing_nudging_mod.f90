!> @file large_scale_forcing_nudging_mod.f90
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
! $Id: large_scale_forcing_nudging_mod.f90 4360 2020-01-07 11:25:50Z suehring $
! Introduction of wall_flags_total_0, which currently sets bits based on static
! topography information used in wall_flags_static_0
! 
! 4329 2019-12-10 15:46:36Z motisi
! Renamed wall_flags_0 to wall_flags_static_0
! 
! 4182 2019-08-22 15:20:23Z scharf
! Corrected "Former revisions" section
! 
! 3719 2019-02-06 13:10:18Z kanani
! Removed USE cpulog (unused)
! 
! 3655 2019-01-07 16:51:22Z knoop
! unused variables removed
! 2320 2017-07-21 12:47:43Z suehring
! initial revision
!
! Description:
! ------------
!> Calculates large scale forcings (geostrophic wind and subsidence velocity) as 
!> well as surfaces fluxes dependent on time given in an external file (LSF_DATA).
!> Moreover, module contains nudging routines, where u, v, pt and q are nudged 
!> to given profiles on a relaxation timescale tnudge. 
!> Profiles are read in from NUDGING_DATA. 
!> Code is based on Neggers et al. (2012) and also in parts on DALES and UCLA-LES.
!> @todo: Revise reading of ASCII-files
!> @todo: Remove unused variables and control flags
!> @todo: Revise large-scale facing of surface variables
!> @todo: Revise control flags lsf_exception, lsf_surf, lsf_vert, etc. 
!--------------------------------------------------------------------------------!
 MODULE lsf_nudging_mod

    USE arrays_3d,                                                             &
        ONLY:  dzw, e, diss, heatflux_input_conversion, pt, pt_init, q,        &
               q_init, s, tend, u, u_init, ug, v, v_init, vg, w, w_subs,       &
               waterflux_input_conversion, zu, zw                  

    USE control_parameters,                                                    &
        ONLY:  bc_lr, bc_ns, bc_pt_b, bc_q_b, constant_diffusion,              &
               constant_heatflux, constant_waterflux,                          &
               data_output_pr, dt_3d, end_time,                                &
               humidity, initializing_actions, intermediate_timestep_count,    &
               ibc_pt_b, ibc_q_b,                                              &
               large_scale_forcing, large_scale_subsidence, lsf_surf, lsf_vert,&
               lsf_exception, message_string, neutral,                         &
               nudging, passive_scalar, pt_surface, ocean_mode, q_surface,     &
               surface_heatflux, surface_pressure, surface_waterflux,          &
               topography, use_subsidence_tendencies
               
    USE grid_variables

    USE indices,                                                               &
        ONLY:  nbgp, ngp_sums_ls, nx, nxl, nxlg, nxlu, nxr, nxrg, ny, nys,     &
               nysv, nysg, nyn, nyng, nzb, nz, nzt, wall_flags_total_0

    USE kinds

    USE pegrid

    USE surface_mod,                                                           &
        ONLY:  surf_def_h, surf_lsm_h, surf_usm_h

    USE statistics,                                                            &
        ONLY:  hom, statistic_regions, sums_ls_l, weight_substep

    INTEGER(iwp) ::  nlsf = 1000                       !< maximum number of profiles in LSF_DATA (large scale forcing)
    INTEGER(iwp) ::  ntnudge = 1000                    !< maximum number of profiles in NUDGING_DATA (nudging)

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ptnudge     !< vertical profile of pot. temperature interpolated to vertical grid (nudging)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  qnudge      !< vertical profile of water vapor mixing ratio interpolated to vertical grid (nudging) 
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  tnudge      !< vertical profile of nudging time scale interpolated to vertical grid (nudging)  
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  td_lsa_lpt  !< temperature tendency due to large scale advection (large scale forcing)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  td_lsa_q    !< water vapor mixing ratio tendency due to large scale advection (large scale forcing)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  td_sub_lpt  !< temperature tendency due to subsidence/ascent (large scale forcing)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  td_sub_q    !< water vapor mixing ratio tendency due to subsidence/ascent (large scale forcing)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ug_vert     !< vertical profile of geostrophic wind component in x-direction interpolated to vertical grid (large scale forcing)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  unudge      !< vertical profile of wind component in x-direction interpolated to vertical grid (nudging) 
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  vnudge      !< vertical profile of wind component in y-direction interpolated to vertical grid (nudging) 
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  vg_vert     !< vertical profile of geostrophic wind component in y-direction interpolated to vertical grid (large scale forcing)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  wnudge      !< vertical profile of subsidence/ascent velocity interpolated to vertical grid (nudging) ???
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  wsubs_vert  !< vertical profile of wind component in z-direction interpolated to vertical grid (nudging) ???

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  shf_surf      !< time-dependent surface sensible heat flux (large scale forcing)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  timenudge     !< times at which vertical profiles are defined in NUDGING_DATA (nudging)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  time_surf     !< times at which surface values/fluxes are defined in LSF_DATA (large scale forcing)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  time_vert     !< times at which vertical profiles are defined in LSF_DATA (large scale forcing)

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  tmp_tnudge    !< current nudging time scale

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  p_surf        !< time-dependent surface pressure (large scale forcing)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  pt_surf       !< time-dependent surface temperature (large scale forcing)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  qsws_surf     !< time-dependent surface latent heat flux (large scale forcing)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  q_surf        !< time-dependent surface water vapor mixing ratio (large scale forcing)

    SAVE
    PRIVATE
!
!-- Public subroutines
    PUBLIC calc_tnudge, ls_forcing_surf, ls_forcing_vert, ls_advec, lsf_init,  &
           lsf_nudging_check_parameters, nudge_init,                           &
           lsf_nudging_check_data_output_pr, lsf_nudging_header,               &
           nudge, nudge_ref
           
!
!-- Public variables
    PUBLIC qsws_surf, shf_surf, td_lsa_lpt, td_lsa_q, td_sub_lpt,              &
           td_sub_q, time_vert


    INTERFACE ls_advec
       MODULE PROCEDURE ls_advec
       MODULE PROCEDURE ls_advec_ij
    END INTERFACE ls_advec

    INTERFACE nudge
       MODULE PROCEDURE nudge
       MODULE PROCEDURE nudge_ij
    END INTERFACE nudge

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE lsf_nudging_check_parameters

       IMPLICIT NONE
!
!--    Check nudging and large scale forcing from external file
       IF ( nudging  .AND.  (  .NOT.  large_scale_forcing ) )  THEN
          message_string = 'Nudging requires large_scale_forcing = .T.. &'//   &
                        'Surface fluxes and geostrophic wind should be &'//    &
                        'prescribed in file LSF_DATA'
          CALL message( 'check_parameters', 'PA0374', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( large_scale_forcing  .AND.  ( bc_lr /= 'cyclic'  .OR.              &
                                          bc_ns /= 'cyclic' ) )  THEN
          message_string = 'Non-cyclic lateral boundaries do not allow for &'//&
                        'the usage of large scale forcing from external file.'
          CALL message( 'check_parameters', 'PA0375', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( large_scale_forcing  .AND.  (  .NOT.  humidity ) )  THEN
          message_string = 'The usage of large scale forcing from external &'//&
                        'file LSF_DATA requires humidity = .T..'
          CALL message( 'check_parameters', 'PA0376', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( large_scale_forcing  .AND.  passive_scalar )  THEN
          message_string = 'The usage of large scale forcing from external &'// &
                        'file LSF_DATA is not implemented for passive scalars'
          CALL message( 'check_parameters', 'PA0440', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( large_scale_forcing  .AND.  topography /= 'flat'                   &
                              .AND.  .NOT.  lsf_exception )  THEN
          message_string = 'The usage of large scale forcing from external &'//&
                        'file LSF_DATA is not implemented for non-flat topography'
          CALL message( 'check_parameters', 'PA0377', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( large_scale_forcing  .AND.  ocean_mode )  THEN
          message_string = 'The usage of large scale forcing from external &'//&
                        'file LSF_DATA is not implemented for ocean mode'
          CALL message( 'check_parameters', 'PA0378', 1, 2, 0, 6, 0 )
       ENDIF

    END SUBROUTINE lsf_nudging_check_parameters

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check data output of profiles for land surface model
!------------------------------------------------------------------------------!
    SUBROUTINE lsf_nudging_check_data_output_pr( variable, var_count, unit,    &
                                                 dopr_unit )
 
       USE profil_parameter

       IMPLICIT NONE
   
       CHARACTER (LEN=*) ::  unit      !< 
       CHARACTER (LEN=*) ::  variable  !< 
       CHARACTER (LEN=*) ::  dopr_unit !< local value of dopr_unit
 
       INTEGER(iwp) ::  var_count     !< 

       SELECT CASE ( TRIM( variable ) )
       

          CASE ( 'td_lsa_thetal' )
             IF (  .NOT.  large_scale_forcing )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) //          &
                                 ' is not implemented for ' //                 &
                                 'large_scale_forcing = .FALSE.'
                CALL message( 'lsf_nudging_check_data_output_pr', 'PA0393',    &
                               1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 81
                dopr_unit             = 'K/s'
                unit                  = 'K/s'
                hom(:,2,81,:) = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'td_lsa_q' )
             IF (  .NOT.  large_scale_forcing )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) //          &
                                 ' is not implemented for ' //                 &
                                 'large_scale_forcing = .FALSE.'
                CALL message( 'lsf_nudging_check_data_output_pr', 'PA0393',    &
                               1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 82
                dopr_unit             = 'kg/kgs'
                unit                  = 'kg/kgs'
                hom(:,2,82,:) = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF
          CASE ( 'td_sub_thetal' )
             IF (  .NOT.  large_scale_forcing )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) //          &
                                 ' is not implemented for ' //                 &
                                 'large_scale_forcing = .FALSE.'
                CALL message( 'lsf_nudging_check_data_output_pr', 'PA0393',    &
                               1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 83
                dopr_unit             = 'K/s'
                unit                  = 'K/s'
                hom(:,2,83,:) = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'td_sub_q' )
             IF (  .NOT.  large_scale_forcing )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) //          &
                                 ' is not implemented for ' //                 &
                                 'large_scale_forcing = .FALSE.'
                CALL message( 'lsf_nudging_check_data_output_pr', 'PA0393',    &
                               1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 84
                dopr_unit             = 'kg/kgs'
                unit                  = 'kg/kgs'
                hom(:,2,84,:) = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'td_nud_thetal' )
             IF (  .NOT.  nudging )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) //          &
                                 ' is not implemented for ' //                 &
                                 'nudging = .FALSE.'
                CALL message( 'lsf_nudging_check_data_output_pr', 'PA0394',    &
                               1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 85
                dopr_unit             = 'K/s'
                unit                  = 'K/s'
                hom(:,2,85,:) = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'td_nud_q' )
             IF (  .NOT.  nudging )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) //          &
                                 ' is not implemented for ' //                 &
                                 'nudging = .FALSE.'
                CALL message( 'lsf_nudging_check_data_output_pr', 'PA0394',    &
                               1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 86
                dopr_unit             = 'kg/kgs'
                unit                  = 'kg/kgs'
                hom(:,2,86,:) = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'td_nud_u' )
             IF (  .NOT.  nudging )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) //          &
                                 ' is not implemented for ' //                 &
                                 'nudging = .FALSE.'
                CALL message( 'lsf_nudging_check_data_output_pr', 'PA0394',    &
                               1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 87
                dopr_unit             = 'm/s2'
                unit                  = 'm/s2'
                hom(:,2,87,:) = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'td_nud_v' )
             IF (  .NOT.  nudging )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) //          &
                                 ' is not implemented for ' //                 &
                                 'nudging = .FALSE.'
                CALL message( 'lsf_nudging_check_data_output_pr', 'PA0394',    &
                               1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 88
                dopr_unit             = 'm/s2'
                unit                  = 'm/s2'
                hom(:,2,88,:) = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF


          CASE DEFAULT
             unit = 'illegal'
   
       END SELECT

    END SUBROUTINE lsf_nudging_check_data_output_pr

!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE lsf_nudging_header ( io )

       IMPLICIT NONE

       INTEGER(iwp), INTENT(IN) ::  io !< Unit of the output file

       WRITE ( io, 1 )
       IF ( large_scale_forcing )  THEN
          WRITE ( io, 3 )
          WRITE ( io, 4 )

          IF ( large_scale_subsidence )  THEN
             IF ( .NOT. use_subsidence_tendencies )  THEN
                WRITE ( io, 5 )
             ELSE
                WRITE ( io, 6 )
             ENDIF 
          ENDIF

          IF ( bc_pt_b == 'dirichlet' )  THEN
             WRITE ( io, 12 )
          ELSEIF ( bc_pt_b == 'neumann' )  THEN
             WRITE ( io, 13 )
          ENDIF

          IF ( bc_q_b == 'dirichlet' )  THEN
             WRITE ( io, 14 )
          ELSEIF ( bc_q_b == 'neumann' )  THEN
             WRITE ( io, 15 )
          ENDIF

          WRITE ( io, 7 )
          IF ( nudging )  THEN
             WRITE ( io, 10 )
          ENDIF
       ELSE
          WRITE ( io, 2 )
          WRITE ( io, 11 )
       ENDIF
       IF ( large_scale_subsidence )  THEN
          WRITE ( io, 8 )
          WRITE ( io, 9 )
       ENDIF


1 FORMAT (//' Large scale forcing and nudging:'/ &
              ' -------------------------------'/)
2 FORMAT (' --> No large scale forcing from external is used (default) ')
3 FORMAT (' --> Large scale forcing from external file LSF_DATA is used: ')
4 FORMAT ('     - large scale advection tendencies ')
5 FORMAT ('     - large scale subsidence velocity w_subs ')
6 FORMAT ('     - large scale subsidence tendencies ')
7 FORMAT ('     - and geostrophic wind components ug and vg')
8 FORMAT (' --> Large-scale vertical motion is used in the ', &
                  'prognostic equation(s) for')
9 FORMAT ('     the scalar(s) only')
10 FORMAT (' --> Nudging is used')
11 FORMAT (' --> No nudging is used (default) ')
12 FORMAT ('     - prescribed surface values for temperature')
13 FORMAT ('     - prescribed surface fluxes for temperature')
14 FORMAT ('     - prescribed surface values for humidity')
15 FORMAT ('     - prescribed surface fluxes for humidity')

    END SUBROUTINE lsf_nudging_header 

!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE lsf_init

       IMPLICIT NONE

       CHARACTER(100) ::  chmess      !<
       CHARACTER(1)   ::  hash        !<

       INTEGER(iwp) ::  ierrn         !<
       INTEGER(iwp) ::  finput = 90   !<
       INTEGER(iwp) ::  k             !<
       INTEGER(iwp) ::  nt            !<

       REAL(wp) ::  fac               !<
       REAL(wp) ::  highheight        !<
       REAL(wp) ::  highug_vert       !<
       REAL(wp) ::  highvg_vert       !<
       REAL(wp) ::  highwsubs_vert    !<
       REAL(wp) ::  lowheight         !<
       REAL(wp) ::  lowug_vert        !<
       REAL(wp) ::  lowvg_vert        !<
       REAL(wp) ::  lowwsubs_vert     !<
       REAL(wp) ::  high_td_lsa_lpt   !<
       REAL(wp) ::  low_td_lsa_lpt    !<
       REAL(wp) ::  high_td_lsa_q     !<
       REAL(wp) ::  low_td_lsa_q      !<
       REAL(wp) ::  high_td_sub_lpt   !<
       REAL(wp) ::  low_td_sub_lpt    !<
       REAL(wp) ::  high_td_sub_q     !<
       REAL(wp) ::  low_td_sub_q      !<
       REAL(wp) ::  r_dummy           !<

       ALLOCATE( p_surf(0:nlsf), pt_surf(0:nlsf), q_surf(0:nlsf),              &
                 qsws_surf(0:nlsf), shf_surf(0:nlsf),                          &
                 td_lsa_lpt(nzb:nzt+1,0:nlsf), td_lsa_q(nzb:nzt+1,0:nlsf),     &
                 td_sub_lpt(nzb:nzt+1,0:nlsf), td_sub_q(nzb:nzt+1,0:nlsf),     &
                 time_vert(0:nlsf), time_surf(0:nlsf),                         &
                 ug_vert(nzb:nzt+1,0:nlsf), vg_vert(nzb:nzt+1,0:nlsf),         &
                 wsubs_vert(nzb:nzt+1,0:nlsf) )

       p_surf = 0.0_wp; pt_surf = 0.0_wp; q_surf = 0.0_wp; qsws_surf = 0.0_wp
       shf_surf = 0.0_wp; time_vert = 0.0_wp; td_lsa_lpt = 0.0_wp
       td_lsa_q = 0.0_wp; td_sub_lpt = 0.0_wp; td_sub_q = 0.0_wp
       time_surf = 0.0_wp; ug_vert = 0.0_wp; vg_vert = 0.0_wp
       wsubs_vert = 0.0_wp

!
!--    Array for storing large scale forcing and nudging tendencies at each 
!--    timestep for data output
       ALLOCATE( sums_ls_l(nzb:nzt+1,0:7) )
       sums_ls_l = 0.0_wp

       ngp_sums_ls = (nz+2)*6

       OPEN ( finput, FILE='LSF_DATA', STATUS='OLD', &
              FORM='FORMATTED', IOSTAT=ierrn )

       IF ( ierrn /= 0 )  THEN
          message_string = 'file LSF_DATA does not exist'
          CALL message( 'ls_forcing', 'PA0368', 1, 2, 0, 6, 0 )
       ENDIF

       ierrn = 0
!
!--    First three lines of LSF_DATA contain header
       READ ( finput, FMT='(a100)', IOSTAT=ierrn ) chmess
       READ ( finput, FMT='(a100)', IOSTAT=ierrn ) chmess
       READ ( finput, FMT='(a100)', IOSTAT=ierrn ) chmess

       IF ( ierrn /= 0 )  THEN
          message_string = 'errors in file LSF_DATA'
          CALL message( 'ls_forcing', 'PA0369', 1, 2, 0, 6, 0 )
       ENDIF

!
!--    Surface values are read in
       nt     = 0
       ierrn = 0

       DO WHILE ( time_surf(nt) < end_time )
          nt = nt + 1
          READ ( finput, *, IOSTAT = ierrn ) time_surf(nt), shf_surf(nt),      &
                                             qsws_surf(nt), pt_surf(nt),       &
                                             q_surf(nt), p_surf(nt)

          IF ( ierrn /= 0 )  THEN
            WRITE ( message_string, * ) 'No time dependent surface ' //        &
                              'variables in & LSF_DATA for end of run found'

             CALL message( 'ls_forcing', 'PA0363', 1, 2, 0, 6, 0 )
          ENDIF
       ENDDO

       IF ( time_surf(1) > end_time )  THEN
          WRITE ( message_string, * ) 'Time dependent surface variables in ' //&
                                      '&LSF_DATA set in after end of ' ,       &
                                      'simulation - lsf_surf is set to FALSE'
          CALL message( 'ls_forcing', 'PA0371', 0, 0, 0, 6, 0 )
          lsf_surf = .FALSE.
       ENDIF

!
!--    Go to the end of the list with surface variables
       DO WHILE ( ierrn == 0 )
          READ ( finput, *, IOSTAT = ierrn ) r_dummy
       ENDDO

!
!--    Profiles of ug, vg and w_subs are read in (large scale forcing)

       nt = 0
       DO WHILE ( time_vert(nt) < end_time )
          nt = nt + 1
          hash = "#"
          ierrn = 1 ! not zero
!
!--       Search for the next line consisting of "# time", 
!--       from there onwards the profiles will be read
          DO WHILE ( .NOT. ( hash == "#" .AND. ierrn == 0 ) ) 
             READ ( finput, *, IOSTAT=ierrn ) hash, time_vert(nt)
             IF ( ierrn < 0 )  THEN 
                WRITE( message_string, * ) 'No time dependent vertical profiles',&
                                 ' in & LSF_DATA for end of run found'
                CALL message( 'ls_forcing', 'PA0372', 1, 2, 0, 6, 0 )
             ENDIF
          ENDDO

          IF ( nt == 1 .AND. time_vert(nt) > end_time ) EXIT

          READ ( finput, *, IOSTAT=ierrn ) lowheight, lowug_vert, lowvg_vert,  &
                                           lowwsubs_vert, low_td_lsa_lpt,      &
                                           low_td_lsa_q, low_td_sub_lpt,       &
                                           low_td_sub_q
          IF ( ierrn /= 0 )  THEN
             message_string = 'errors in file LSF_DATA'
             CALL message( 'ls_forcing', 'PA0369', 1, 2, 0, 6, 0 )
          ENDIF

          READ ( finput, *, IOSTAT=ierrn ) highheight, highug_vert,            &
                                           highvg_vert, highwsubs_vert,        &
                                           high_td_lsa_lpt, high_td_lsa_q,     &
                                           high_td_sub_lpt, high_td_sub_q
       
          IF ( ierrn /= 0 )  THEN
             message_string = 'errors in file LSF_DATA'
             CALL message( 'ls_forcing', 'PA0369', 1, 2, 0, 6, 0 )
          ENDIF


          DO  k = nzb, nzt+1
             IF ( highheight < zu(k) )  THEN
                lowheight      = highheight
                lowug_vert     = highug_vert
                lowvg_vert     = highvg_vert
                lowwsubs_vert  = highwsubs_vert
                low_td_lsa_lpt = high_td_lsa_lpt
                low_td_lsa_q   = high_td_lsa_q
                low_td_sub_lpt = high_td_sub_lpt
                low_td_sub_q   = high_td_sub_q

                ierrn = 0
                READ ( finput, *, IOSTAT=ierrn ) highheight, highug_vert,      &
                                                 highvg_vert, highwsubs_vert,  &
                                                 high_td_lsa_lpt,              &
                                                 high_td_lsa_q,                &
                                                 high_td_sub_lpt, high_td_sub_q

                IF ( ierrn /= 0 )  THEN
                   WRITE( message_string, * ) 'zu(',k,') = ', zu(k), 'm ',     &
                        'is higher than the maximum height in LSF_DATA ',      &
                        'which is ', lowheight, 'm. Interpolation on PALM ',   &
                        'grid is not possible.'
                   CALL message( 'ls_forcing', 'PA0395', 1, 2, 0, 6, 0 )
                ENDIF

             ENDIF

!
!--          Interpolation of prescribed profiles in space 
             fac = (highheight-zu(k))/(highheight - lowheight)

             ug_vert(k,nt)    = fac * lowug_vert                               &
                                + ( 1.0_wp - fac ) * highug_vert
             vg_vert(k,nt)    = fac * lowvg_vert                               &
                                + ( 1.0_wp - fac ) * highvg_vert
             wsubs_vert(k,nt) = fac * lowwsubs_vert                            &
                                + ( 1.0_wp - fac ) * highwsubs_vert

             td_lsa_lpt(k,nt) = fac * low_td_lsa_lpt                           &
                                + ( 1.0_wp - fac ) * high_td_lsa_lpt
             td_lsa_q(k,nt)   = fac * low_td_lsa_q                             &
                                + ( 1.0_wp - fac ) * high_td_lsa_q
             td_sub_lpt(k,nt) = fac * low_td_sub_lpt                           &
                                + ( 1.0_wp - fac ) * high_td_sub_lpt
             td_sub_q(k,nt)   = fac * low_td_sub_q                             &
                                + ( 1.0_wp - fac ) * high_td_sub_q

          ENDDO

       ENDDO 

!
!--    Large scale vertical velocity has to be zero at the surface
       wsubs_vert(nzb,:) = 0.0_wp
    
       IF ( time_vert(1) > end_time )  THEN
          WRITE ( message_string, * ) 'Time dependent large scale profile ',   &
                             'forcing from&LSF_DATA sets in after end of ' ,   &
                             'simulation - lsf_vert is set to FALSE'
          CALL message( 'ls_forcing', 'PA0373', 0, 0, 0, 6, 0 )
          lsf_vert = .FALSE.
       ENDIF

       CLOSE( finput )

    END SUBROUTINE lsf_init

!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE ls_forcing_surf ( time )

       IMPLICIT NONE

       INTEGER(iwp) ::  nt                     !<

       REAL(wp)             :: dum_surf_flux  !<
       REAL(wp)             :: fac            !<
       REAL(wp), INTENT(in) :: time           !<

!
!--    Interpolation in time of LSF_DATA at the surface
       nt = 1
       DO WHILE ( time > time_surf(nt) )
          nt = nt + 1
       ENDDO
       IF ( time /= time_surf(nt) )  THEN
         nt = nt - 1
       ENDIF

       fac = ( time -time_surf(nt) ) / ( time_surf(nt+1) - time_surf(nt) )

       IF ( ibc_pt_b == 0 )  THEN
!
!--       In case of Dirichlet boundary condition shf must not 
!--       be set - it is calculated via MOST in prandtl_fluxes
          pt_surface = pt_surf(nt) + fac * ( pt_surf(nt+1) - pt_surf(nt) )

       ELSEIF ( ibc_pt_b == 1 )  THEN
!
!--       In case of Neumann boundary condition pt_surface is needed for 
!--       calculation of reference density
          dum_surf_flux = ( shf_surf(nt) + fac *                               &
                            ( shf_surf(nt+1) - shf_surf(nt) )                  &
                          ) * heatflux_input_conversion(nzb)
!
!--       Save surface sensible heat flux on default, natural and urban surface
!--       type, if required 
          IF ( surf_def_h(0)%ns >= 1 )  surf_def_h(0)%shf(:) = dum_surf_flux
          IF ( surf_lsm_h%ns    >= 1 )  surf_lsm_h%shf(:)    = dum_surf_flux
          IF ( surf_usm_h%ns    >= 1 )  surf_usm_h%shf(:)    = dum_surf_flux

          pt_surface    = pt_surf(nt) + fac * ( pt_surf(nt+1) - pt_surf(nt) )

       ENDIF

       IF ( ibc_q_b == 0 )  THEN
!
!--       In case of Dirichlet boundary condition qsws must not 
!--       be set - it is calculated via MOST in prandtl_fluxes
          q_surface = q_surf(nt) + fac * ( q_surf(nt+1) - q_surf(nt) )

       ELSEIF ( ibc_q_b == 1 )  THEN
          dum_surf_flux = ( qsws_surf(nt) + fac *                              &
                             ( qsws_surf(nt+1) - qsws_surf(nt) )               &
                             ) * waterflux_input_conversion(nzb)
!
!--       Save surface sensible heat flux on default, natural and urban surface
!--       type, if required 
          IF ( surf_def_h(0)%ns >= 1 )  surf_def_h(0)%qsws(:) = dum_surf_flux
          IF ( surf_lsm_h%ns    >= 1 )  surf_lsm_h%qsws(:)    = dum_surf_flux
          IF ( surf_usm_h%ns    >= 1 )  surf_usm_h%qsws(:)    = dum_surf_flux

       ENDIF
!
!--    Surface heat- and waterflux will be written later onto surface elements 
       IF ( .NOT.  neutral  .AND.  constant_heatflux  .AND.                    &
            TRIM( initializing_actions ) /= 'read_restart_data' )  THEN
             surface_heatflux = shf_surf(1)
       ENDIF
       IF ( humidity  .AND.  constant_waterflux  .AND.                         &
            TRIM( initializing_actions ) /= 'read_restart_data' )  THEN
             surface_waterflux = qsws_surf(1)
       ENDIF

       surface_pressure = p_surf(nt) + fac * ( p_surf(nt+1) - p_surf(nt) )

    END SUBROUTINE ls_forcing_surf 




!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE ls_forcing_vert ( time )


       IMPLICIT NONE

       INTEGER(iwp) ::  nt                     !<

       REAL(wp)             ::  fac           !<
       REAL(wp), INTENT(in) ::  time          !<

!
!--    Interpolation in time of LSF_DATA for ug, vg and w_subs
       nt = 1
       DO WHILE ( time > time_vert(nt) )
          nt = nt + 1
       ENDDO
       IF ( time /= time_vert(nt) )  THEN
         nt = nt - 1
       ENDIF

       fac = ( time-time_vert(nt) ) / ( time_vert(nt+1)-time_vert(nt) )

       ug     = ug_vert(:,nt) + fac * ( ug_vert(:,nt+1) - ug_vert(:,nt) )
       vg     = vg_vert(:,nt) + fac * ( vg_vert(:,nt+1) - vg_vert(:,nt) )

       IF ( large_scale_subsidence )  THEN
          w_subs = wsubs_vert(:,nt)                                            &
                   + fac * ( wsubs_vert(:,nt+1) - wsubs_vert(:,nt) )
       ENDIF

    END SUBROUTINE ls_forcing_vert


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE ls_advec ( time, prog_var )
      

       IMPLICIT NONE

       CHARACTER (LEN=*) ::  prog_var   !< 

       REAL(wp), INTENT(in)  :: time    !< 
       REAL(wp) :: fac                  !<  

       INTEGER(iwp) ::  i               !< 
       INTEGER(iwp) ::  j               !< 
       INTEGER(iwp) ::  k               !< 
       INTEGER(iwp) ::  nt               !< 

!
!--    Interpolation in time of LSF_DATA 
       nt = 1
       DO WHILE ( time > time_vert(nt) )
          nt = nt + 1
       ENDDO
       IF ( time /= time_vert(nt) )  THEN
         nt = nt - 1
       ENDIF

       fac = ( time-time_vert(nt) ) / ( time_vert(nt+1)-time_vert(nt) )

!
!--    Add horizontal large scale advection tendencies of pt and q 
       SELECT CASE ( prog_var )

          CASE ( 'pt' )

             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt
                      tend(k,j,i) = tend(k,j,i) + td_lsa_lpt(k,nt) + fac *     &
                                    ( td_lsa_lpt(k,nt+1) - td_lsa_lpt(k,nt) ) *&
                                        MERGE( 1.0_wp, 0.0_wp,                 &
                                        BTEST( wall_flags_total_0(k,j,i), 0 ) )
                   ENDDO
                ENDDO
             ENDDO

          CASE ( 'q' )

             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt
                      tend(k,j,i) = tend(k,j,i) + td_lsa_q(k,nt) + fac *       &
                                    ( td_lsa_q(k,nt+1) - td_lsa_q(k,nt) ) *    &
                                        MERGE( 1.0_wp, 0.0_wp,                 &
                                        BTEST( wall_flags_total_0(k,j,i), 0 ) )
                   ENDDO
                ENDDO
             ENDDO

       END SELECT

!
!--    Subsidence of pt and q with prescribed subsidence tendencies
       IF ( large_scale_subsidence .AND. use_subsidence_tendencies )  THEN

          SELECT CASE ( prog_var )

             CASE ( 'pt' )

                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         tend(k,j,i) = tend(k,j,i) + td_sub_lpt(k,nt) + fac *  &
                                     ( td_sub_lpt(k,nt+1) - td_sub_lpt(k,nt) )*&
                                        MERGE( 1.0_wp, 0.0_wp,                 &
                                        BTEST( wall_flags_total_0(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
  
             CASE ( 'q' )

                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         tend(k,j,i) = tend(k,j,i) + td_sub_q(k,nt) + fac *    &
                                       ( td_sub_q(k,nt+1) - td_sub_q(k,nt) ) * &
                                        MERGE( 1.0_wp, 0.0_wp,                 &
                                        BTEST( wall_flags_total_0(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO

          END SELECT

       ENDIF

    END SUBROUTINE ls_advec


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE ls_advec_ij ( i, j, time, prog_var )

       IMPLICIT NONE

       CHARACTER (LEN=*) ::  prog_var   !< 

       REAL(wp), INTENT(in)  :: time    !< 
       REAL(wp) :: fac                  !< 

       INTEGER(iwp) ::  i               !< 
       INTEGER(iwp) ::  j               !< 
       INTEGER(iwp) ::  k               !< 
       INTEGER(iwp) ::  nt               !< 

!
!--    Interpolation in time of LSF_DATA 
       nt = 1
       DO WHILE ( time > time_vert(nt) )
          nt = nt + 1
       ENDDO
       IF ( time /= time_vert(nt) )  THEN
         nt = nt - 1
       ENDIF

       fac = ( time-time_vert(nt) ) / ( time_vert(nt+1)-time_vert(nt) )

!
!--    Add horizontal large scale advection tendencies of pt and q 
       SELECT CASE ( prog_var )

          CASE ( 'pt' )

             DO  k = nzb+1, nzt
                tend(k,j,i) = tend(k,j,i) + td_lsa_lpt(k,nt)                   &
                             + fac * ( td_lsa_lpt(k,nt+1) - td_lsa_lpt(k,nt) )*&
                                        MERGE( 1.0_wp, 0.0_wp,                 &
                                        BTEST( wall_flags_total_0(k,j,i), 0 ) )
             ENDDO

          CASE ( 'q' )

             DO  k = nzb+1, nzt
                tend(k,j,i) = tend(k,j,i) + td_lsa_q(k,nt)                     &
                              + fac * ( td_lsa_q(k,nt+1) - td_lsa_q(k,nt) ) *  &
                                        MERGE( 1.0_wp, 0.0_wp,                 &
                                        BTEST( wall_flags_total_0(k,j,i), 0 ) )
             ENDDO

       END SELECT

!
!--    Subsidence of pt and q with prescribed profiles
       IF ( large_scale_subsidence .AND. use_subsidence_tendencies )  THEN

          SELECT CASE ( prog_var )

             CASE ( 'pt' )

                DO  k = nzb+1, nzt
                   tend(k,j,i) = tend(k,j,i) + td_sub_lpt(k,nt) + fac *        &
                                 ( td_sub_lpt(k,nt+1) - td_sub_lpt(k,nt) ) *   &
                                        MERGE( 1.0_wp, 0.0_wp,                 &
                                        BTEST( wall_flags_total_0(k,j,i), 0 ) )
                ENDDO
  
             CASE ( 'q' )

                DO  k = nzb+1, nzt
                   tend(k,j,i) = tend(k,j,i) + td_sub_q(k,nt) + fac *          &
                                 ( td_sub_q(k,nt+1) - td_sub_q(k,nt) ) *       &
                                        MERGE( 1.0_wp, 0.0_wp,                 &
                                        BTEST( wall_flags_total_0(k,j,i), 0 ) )
                ENDDO

          END SELECT

       ENDIF

    END SUBROUTINE ls_advec_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE nudge_init

       IMPLICIT NONE


       INTEGER(iwp) ::  finput = 90  !<
       INTEGER(iwp) ::  ierrn        !<
       INTEGER(iwp) ::  k            !<
       INTEGER(iwp) ::  nt            !<

       CHARACTER(1) ::  hash     !<

       REAL(wp) ::  highheight   !<
       REAL(wp) ::  highqnudge   !<
       REAL(wp) ::  highptnudge  !<
       REAL(wp) ::  highunudge   !<
       REAL(wp) ::  highvnudge   !<
       REAL(wp) ::  highwnudge   !<
       REAL(wp) ::  hightnudge   !<

       REAL(wp) ::  lowheight    !<
       REAL(wp) ::  lowqnudge    !<
       REAL(wp) ::  lowptnudge   !<
       REAL(wp) ::  lowunudge    !<
       REAL(wp) ::  lowvnudge    !<
       REAL(wp) ::  lowwnudge    !<
       REAL(wp) ::  lowtnudge    !<

       REAL(wp) ::  fac          !<

       ALLOCATE( ptnudge(nzb:nzt+1,1:ntnudge), qnudge(nzb:nzt+1,1:ntnudge), &
                 tnudge(nzb:nzt+1,1:ntnudge), unudge(nzb:nzt+1,1:ntnudge),  &
                 vnudge(nzb:nzt+1,1:ntnudge), wnudge(nzb:nzt+1,1:ntnudge)  )

       ALLOCATE( tmp_tnudge(nzb:nzt) )

       ALLOCATE( timenudge(0:ntnudge) )

       ptnudge = 0.0_wp; qnudge = 0.0_wp; tnudge = 0.0_wp; unudge = 0.0_wp
       vnudge = 0.0_wp; wnudge = 0.0_wp; timenudge = 0.0_wp
!
!--    Initialize array tmp_nudge with a current nudging time scale of 6 hours
       tmp_tnudge = 21600.0_wp

       nt = 0
       OPEN ( finput, FILE='NUDGING_DATA', STATUS='OLD', &
              FORM='FORMATTED', IOSTAT=ierrn )

       IF ( ierrn /= 0 )  THEN
          message_string = 'file NUDGING_DATA does not exist'
          CALL message( 'nudging', 'PA0365', 1, 2, 0, 6, 0 )
       ENDIF

       ierrn = 0

 rloop:DO
          nt = nt + 1
          hash = "#"
          ierrn = 1 ! not zero
!
!--       Search for the next line consisting of "# time", 
!--       from there onwards the profiles will be read
          DO WHILE ( .NOT. ( hash == "#" .AND. ierrn == 0 ) ) 
          
            READ ( finput, *, IOSTAT=ierrn ) hash, timenudge(nt)
            IF ( ierrn < 0 )  EXIT rloop

          ENDDO

          ierrn = 0
          READ ( finput, *, IOSTAT=ierrn ) lowheight, lowtnudge, lowunudge,   &
                                           lowvnudge, lowwnudge , lowptnudge, &
                                           lowqnudge

          IF ( ierrn /= 0 )  THEN
             message_string = 'errors in file NUDGING_DATA'
             CALL message( 'nudging', 'PA0366', 1, 2, 0, 6, 0 )
          ENDIF

          ierrn = 0
          READ ( finput, *, IOSTAT=ierrn ) highheight, hightnudge, highunudge,   &
                                           highvnudge, highwnudge , highptnudge, &
                                           highqnudge

          IF ( ierrn /= 0 )  THEN
             message_string = 'errors in file NUDGING_DATA'
             CALL message( 'nudging', 'PA0366', 1, 2, 0, 6, 0 )
          ENDIF

          DO  k = nzb, nzt+1
             DO WHILE ( highheight < zu(k) )
                lowheight  = highheight
                lowtnudge  = hightnudge
                lowunudge  = highunudge
                lowvnudge  = highvnudge
                lowwnudge  = highwnudge
                lowptnudge = highptnudge
                lowqnudge  = highqnudge
 
                ierrn = 0
                READ ( finput, *, IOSTAT=ierrn )  highheight , hightnudge ,    &
                                                  highunudge , highvnudge ,    &
                                                  highwnudge , highptnudge,    &
                                                  highqnudge
                IF (ierrn /= 0 )  THEN
                   WRITE( message_string, * ) 'zu(',k,') = ', zu(k), 'm is ',  &
                        'higher than the maximum height in NUDING_DATA which ',&
                        'is ', lowheight, 'm. Interpolation on PALM ',         &
                        'grid is not possible.'
                   CALL message( 'nudging', 'PA0364', 1, 2, 0, 6, 0 )
                ENDIF
             ENDDO

!
!--          Interpolation of prescribed profiles in space 

             fac = ( highheight - zu(k) ) / ( highheight - lowheight )

             tnudge(k,nt)  = fac * lowtnudge  + ( 1.0_wp - fac ) * hightnudge
             unudge(k,nt)  = fac * lowunudge  + ( 1.0_wp - fac ) * highunudge
             vnudge(k,nt)  = fac * lowvnudge  + ( 1.0_wp - fac ) * highvnudge
             wnudge(k,nt)  = fac * lowwnudge  + ( 1.0_wp - fac ) * highwnudge
             ptnudge(k,nt) = fac * lowptnudge + ( 1.0_wp - fac ) * highptnudge
             qnudge(k,nt)  = fac * lowqnudge  + ( 1.0_wp - fac ) * highqnudge
          ENDDO

       ENDDO rloop

       CLOSE ( finput )

!
!--    Overwrite initial profiles in case of nudging
       IF ( nudging )  THEN
          pt_init = ptnudge(:,1)
          u_init  = unudge(:,1)
          v_init  = vnudge(:,1)
          IF ( humidity  )  THEN ! is passive_scalar correct???
             q_init = qnudge(:,1)
          ENDIF

          WRITE( message_string, * ) 'Initial profiles of u, v, pt and q ',    &
                                     'from NUDGING_DATA are used.'
          CALL message( 'large_scale_forcing_nudging', 'PA0370', 0, 0, 0, 6, 0 )
       ENDIF


    END SUBROUTINE nudge_init

!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE calc_tnudge ( time )

       IMPLICIT NONE


       REAL(wp) ::  dtm         !<
       REAL(wp) ::  dtp         !<
       REAL(wp) ::  time        !<

       INTEGER(iwp) ::  k   !<
       INTEGER(iwp) ::  nt  !<

       nt = 1
       DO WHILE ( time > timenudge(nt) )
         nt = nt+1
       ENDDO
       IF ( time /= timenudge(1) ) THEN
         nt = nt-1
       ENDIF

       dtm = ( time - timenudge(nt) ) / ( timenudge(nt+1) - timenudge(nt) )
       dtp = ( timenudge(nt+1) - time ) / ( timenudge(nt+1) - timenudge(nt) )

       DO  k = nzb, nzt
          tmp_tnudge(k) = MAX( dt_3d, tnudge(k,nt) * dtp + tnudge(k,nt+1) * dtm )
       ENDDO

    END SUBROUTINE calc_tnudge

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE nudge ( time, prog_var )

       IMPLICIT NONE

       CHARACTER (LEN=*) ::  prog_var  !<

       REAL(wp) ::  tmp_tend    !<
       REAL(wp) ::  dtm         !<
       REAL(wp) ::  dtp         !<
       REAL(wp) ::  time        !<

       INTEGER(iwp) ::  i  !<
       INTEGER(iwp) ::  j  !<
       INTEGER(iwp) ::  k  !<
       INTEGER(iwp) ::  nt  !<


       nt = 1
       DO WHILE ( time > timenudge(nt) )
         nt = nt+1
       ENDDO
       IF ( time /= timenudge(1) ) THEN
         nt = nt-1
       ENDIF

       dtm = ( time - timenudge(nt) ) / ( timenudge(nt+1) - timenudge(nt) )
       dtp = ( timenudge(nt+1) - time ) / ( timenudge(nt+1) - timenudge(nt) )

       SELECT CASE ( prog_var )

          CASE ( 'u' )

             DO  i = nxl, nxr
                DO  j = nys, nyn

                   DO  k = nzb+1, nzt

                      tmp_tend = - ( hom(k,1,1,0) - ( unudge(k,nt) * dtp +     &
                                     unudge(k,nt+1) * dtm ) ) / tmp_tnudge(k)

                      tend(k,j,i) = tend(k,j,i) + tmp_tend *                   &
                                        MERGE( 1.0_wp, 0.0_wp,                 &
                                        BTEST( wall_flags_total_0(k,j,i), 1 ) )

                      sums_ls_l(k,6) = sums_ls_l(k,6) + tmp_tend *             &
                                     weight_substep(intermediate_timestep_count)
                   ENDDO
                  
                   sums_ls_l(nzt+1,6) = sums_ls_l(nzt,6)
 
                ENDDO
            ENDDO

          CASE ( 'v' )

             DO  i = nxl, nxr
                DO  j = nys, nyn

                   DO  k = nzb+1, nzt

                      tmp_tend = - ( hom(k,1,2,0) - ( vnudge(k,nt) * dtp +     &
                                     vnudge(k,nt+1) * dtm ) ) / tmp_tnudge(k)

                      tend(k,j,i) = tend(k,j,i) + tmp_tend *                   &
                                        MERGE( 1.0_wp, 0.0_wp,                 &
                                        BTEST( wall_flags_total_0(k,j,i), 2 ) )

                      sums_ls_l(k,7) = sums_ls_l(k,7) + tmp_tend *             &
                                     weight_substep(intermediate_timestep_count)
                   ENDDO
                  
                   sums_ls_l(nzt+1,7) = sums_ls_l(nzt,7)

                ENDDO
            ENDDO

          CASE ( 'pt' )

             DO  i = nxl, nxr
                DO  j = nys, nyn

                   DO  k = nzb+1, nzt

                      tmp_tend = - ( hom(k,1,4,0) - ( ptnudge(k,nt) * dtp +    &
                                     ptnudge(k,nt+1) * dtm ) ) / tmp_tnudge(k)

                      tend(k,j,i) = tend(k,j,i) + tmp_tend *                   &
                                        MERGE( 1.0_wp, 0.0_wp,                 &
                                        BTEST( wall_flags_total_0(k,j,i), 0 ) )

                      sums_ls_l(k,4) = sums_ls_l(k,4) + tmp_tend *             &
                                     weight_substep(intermediate_timestep_count)
                   ENDDO

                   sums_ls_l(nzt+1,4) = sums_ls_l(nzt,4)

                ENDDO
            ENDDO

          CASE ( 'q' )

             DO  i = nxl, nxr
                DO  j = nys, nyn

                   DO  k = nzb+1, nzt

                      tmp_tend = - ( hom(k,1,41,0) - ( qnudge(k,nt) * dtp +    &
                                     qnudge(k,nt+1) * dtm ) ) / tmp_tnudge(k)

                      tend(k,j,i) = tend(k,j,i) + tmp_tend *                   &
                                        MERGE( 1.0_wp, 0.0_wp,                 &
                                        BTEST( wall_flags_total_0(k,j,i), 0 ) )

                      sums_ls_l(k,5) = sums_ls_l(k,5) + tmp_tend *             &
                                     weight_substep(intermediate_timestep_count)
                   ENDDO
                  
                   sums_ls_l(nzt+1,5) = sums_ls_l(nzt,5)

                ENDDO
            ENDDO

          CASE DEFAULT
             message_string = 'unknown prognostic variable "' // prog_var // '"'
             CALL message( 'nudge', 'PA0367', 1, 2, 0, 6, 0 )

       END SELECT

    END SUBROUTINE nudge


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!------------------------------------------------------------------------------!

    SUBROUTINE nudge_ij( i, j, time, prog_var )

       IMPLICIT NONE


       CHARACTER (LEN=*) ::  prog_var  !<

       REAL(wp) ::  tmp_tend    !<
       REAL(wp) ::  dtm         !<
       REAL(wp) ::  dtp         !<
       REAL(wp) ::  time        !<

       INTEGER(iwp) ::  i  !<
       INTEGER(iwp) ::  j  !<
       INTEGER(iwp) ::  k  !<
       INTEGER(iwp) ::  nt  !<


       nt = 1
       DO WHILE ( time > timenudge(nt) )
         nt = nt+1
       ENDDO
       IF ( time /= timenudge(1) )  THEN
         nt = nt-1
       ENDIF

       dtm = ( time - timenudge(nt) ) / ( timenudge(nt+1) - timenudge(nt) )
       dtp = ( timenudge(nt+1) - time ) / ( timenudge(nt+1) - timenudge(nt) )

       SELECT CASE ( prog_var )

          CASE ( 'u' )

             DO  k = nzb+1, nzt

                tmp_tend = - ( hom(k,1,1,0) - ( unudge(k,nt) * dtp +           &
                               unudge(k,nt+1) * dtm ) ) / tmp_tnudge(k)

                tend(k,j,i) = tend(k,j,i) + tmp_tend *                         &
                                        MERGE( 1.0_wp, 0.0_wp,                 &
                                        BTEST( wall_flags_total_0(k,j,i), 1 ) )

                sums_ls_l(k,6) = sums_ls_l(k,6) + tmp_tend                     &
                                 * weight_substep(intermediate_timestep_count)
             ENDDO

             sums_ls_l(nzt+1,6) = sums_ls_l(nzt,6)

          CASE ( 'v' )

             DO  k = nzb+1, nzt

                tmp_tend = - ( hom(k,1,2,0) - ( vnudge(k,nt) * dtp +           &
                               vnudge(k,nt+1) * dtm ) ) / tmp_tnudge(k)

                tend(k,j,i) = tend(k,j,i) + tmp_tend *                         &
                                        MERGE( 1.0_wp, 0.0_wp,                 &
                                        BTEST( wall_flags_total_0(k,j,i), 2 ) )

                sums_ls_l(k,7) = sums_ls_l(k,7) + tmp_tend                     &
                                 * weight_substep(intermediate_timestep_count)
             ENDDO

             sums_ls_l(nzt+1,7) = sums_ls_l(nzt,7)

          CASE ( 'pt' )

             DO  k = nzb+1, nzt

                tmp_tend = - ( hom(k,1,4,0) - ( ptnudge(k,nt) * dtp +          &
                               ptnudge(k,nt+1) * dtm ) ) / tmp_tnudge(k)

                tend(k,j,i) = tend(k,j,i) + tmp_tend *                         &
                                        MERGE( 1.0_wp, 0.0_wp,                 &
                                        BTEST( wall_flags_total_0(k,j,i), 0 ) )

                sums_ls_l(k,4) = sums_ls_l(k,4) + tmp_tend                     &
                                 * weight_substep(intermediate_timestep_count)
             ENDDO

             sums_ls_l(nzt+1,4) = sums_ls_l(nzt,4)


          CASE ( 'q' )

             DO  k = nzb+1, nzt

                tmp_tend = - ( hom(k,1,41,0) - ( qnudge(k,nt) * dtp +          &
                               qnudge(k,nt+1) * dtm ) ) / tmp_tnudge(k)

                tend(k,j,i) = tend(k,j,i) + tmp_tend *                         &
                                        MERGE( 1.0_wp, 0.0_wp,                 &
                                        BTEST( wall_flags_total_0(k,j,i), 0 ) )

                sums_ls_l(k,5) = sums_ls_l(k,5) + tmp_tend                     &
                                 * weight_substep(intermediate_timestep_count)
             ENDDO

             sums_ls_l(nzt+1,5) = sums_ls_l(nzt,5)

          CASE DEFAULT
             message_string = 'unknown prognostic variable "' // prog_var // '"'
             CALL message( 'nudge', 'PA0367', 1, 2, 0, 6, 0 )

       END SELECT


    END SUBROUTINE nudge_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE nudge_ref ( time )

       IMPLICIT NONE

       INTEGER(iwp) ::  nt                    !<

       REAL(wp)             ::  fac           !<
       REAL(wp), INTENT(in) ::  time          !<

!
!--    Interpolation in time of NUDGING_DATA for pt_init and q_init. This is 
!--    needed for correct upper boundary conditions for pt and q and in case that 
!      large scale subsidence as well as scalar Rayleigh-damping are used
       nt = 1
       DO WHILE ( time > time_vert(nt) )
          nt = nt + 1
       ENDDO
       IF ( time /= time_vert(nt) )  THEN
        nt = nt - 1
       ENDIF

       fac = ( time-time_vert(nt) ) / ( time_vert(nt+1)-time_vert(nt) )

       pt_init = ptnudge(:,nt) + fac * ( ptnudge(:,nt+1) - ptnudge(:,nt) )
       q_init  = qnudge(:,nt) + fac * ( qnudge(:,nt+1) - qnudge(:,nt) )
       u_init  = unudge(:,nt) + fac * ( unudge(:,nt+1) - unudge(:,nt) )
       v_init  = vnudge(:,nt) + fac * ( vnudge(:,nt+1) - vnudge(:,nt) )

    END SUBROUTINE nudge_ref


 END MODULE lsf_nudging_mod
