!> @file turbulence_closure_mod.f90
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
! Copyright 2017-2018 Leibniz Universitaet Hannover
!--------------------------------------------------------------------------------!
!
! Current revisions:
! -----------------
! 
! 
! Former revisions:
! -----------------
! $Id: turbulence_closure_mod.f90 4457 2020-03-11 14:20:43Z raasch $
! use statement for exchange horiz added
! 
! 4433 2020-02-28 22:14:43Z gronemeier
! remove warning for newly implemented RANS mode
!
! 4370 2020-01-10 14:00:44Z raasch
! vector directives added to force vectorization on Intel19 compiler
!
! 4346 2019-12-18 11:55:56Z motisi
! Introduction of wall_flags_total_0, which currently sets bits based on static
! topography information used in wall_flags_static_0
!
! 4329 2019-12-10 15:46:36Z motisi
! Renamed wall_flags_0 to wall_flags_static_0
!
! 4182 2019-08-22 15:20:23Z scharf
! Corrected "Former revisions" section
!
! 4177 2019-08-20 14:32:34Z gronemeier
! add comment
!
! 4170 2019-08-19 17:12:31Z gronemeier
! - add performance optimizations according to K. Ketelsen
!   to diffusion_e and tcm_diffusivities_default
! - bugfix in calculating l_wall for vertical walls
! - bugfix in using l_wall in initialization (consider wall_adjustment_factor)
! - always initialize diss and save the dissipation to that array
!
! 4168 2019-08-16 13:50:17Z suehring
! Replace function get_topography_top_index by topo_top_ind
!
! 4110 2019-07-22 17:05:21Z suehring
! pass integer flag array as well as boundary flags to WS scalar advection
! routine
!
! 4109 2019-07-22 17:00:34Z suehring
! - Modularize setting of boundary conditions for TKE and dissipation
! - Neumann boundary condition for TKE at model top is set also in child domain
! - Revise setting of Neumann boundary conditions at non-cyclic lateral
!   boundaries
! - Bugfix, set Neumann boundary condition for TKE at vertical wall instead of
!   an implicit Dirichlet boundary condition which implied a sink of TKE
!   at vertical walls
!
! 4048 2019-06-21 21:00:21Z knoop
! write out preprocessor directives; remove tailing whitespaces
!
! 3775 2019-03-04 12:40:20Z gronemeier
! removed unused variables
!
! 3724 2019-02-06 16:28:23Z kanani
! Correct double-used log_point_s units
!
! 3719 2019-02-06 13:10:18Z kanani
! Changed log_point to log_point_s, otherwise this overlaps with
! 'all progn.equations' cpu measurement.
!
! 3684 2019-01-20 20:20:58Z knoop
! Remove unused variable simulated_time
!
! 2696 2017-12-14 17:12:51Z kanani
! Initial revision
!
!
! Authors:
! --------
! @author Tobias Gronemeier
! @author Hauke Wurps
!
! Description:
! ------------
!> This module contains the available turbulence closures for PALM.
!>
!>
!> @todo test initialization for all possibilities
!> @todo add OpenMP directives whereever possible
!> @todo Check for random disturbances
!> @note <Enter notes on the module>
!-----------------------------------------------------------------------------!
 MODULE turbulence_closure_mod


    USE arrays_3d,                                                            &
        ONLY:  diss, diss_1, diss_2, diss_3, diss_p, dzu, e, e_1, e_2, e_3,   &
               e_p, kh, km, mean_inflow_profiles, prho, pt, tdiss_m,          &
               te_m, tend, u, v, vpt, w

    USE basic_constants_and_equations_mod,                                    &
        ONLY:  g, kappa, lv_d_cp, lv_d_rd, rd_d_rv

    USE control_parameters,                                                   &
        ONLY:  bc_dirichlet_l,                                                &
               bc_dirichlet_n,                                                &
               bc_dirichlet_r,                                                &
               bc_dirichlet_s,                                                &
               bc_radiation_l,                                                &
               bc_radiation_n,                                                &
               bc_radiation_r,                                                &
               bc_radiation_s,                                                &
               child_domain,                                                  &
               constant_diffusion, dt_3d, e_init, humidity,                   &
               initializing_actions, intermediate_timestep_count,             &
               intermediate_timestep_count_max, km_constant,                  &
               les_dynamic, les_mw, ocean_mode, plant_canopy, prandtl_number, &
               pt_reference, rans_mode, rans_tke_e, rans_tke_l,               &
               timestep_scheme, turbulence_closure,                           &
               turbulent_inflow, use_upstream_for_tke, vpt_reference,         &
               ws_scheme_sca, current_timestep_number

    USE advec_ws,                                                             &
        ONLY:  advec_s_ws

    USE advec_s_bc_mod,                                                       &
        ONLY:  advec_s_bc

    USE advec_s_pw_mod,                                                       &
        ONLY:  advec_s_pw

    USE advec_s_up_mod,                                                       &
        ONLY:  advec_s_up

    USE cpulog,                                                               &
        ONLY:  cpu_log, log_point_s

    USE indices,                                                              &
        ONLY:  advc_flags_s,                                                  &
               nbgp, nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg, nzb, nzt,    &
               topo_top_ind,                                                  &
               wall_flags_total_0

    USE kinds

    USE ocean_mod,                                                            &
        ONLY:  prho_reference

    USE pegrid

    USE plant_canopy_model_mod,                                               &
        ONLY:  pcm_tendency

    USE statistics,                                                           &
        ONLY:  hom, hom_sum, statistic_regions

    USE surface_mod,                                                          &
        ONLY:  bc_h,                                                          &
               bc_v,                                                          &
               surf_def_h,                                                    &
               surf_def_v,                                                    &
               surf_lsm_h,                                                    &
               surf_lsm_v,                                                    &
               surf_usm_h,                                                    &
               surf_usm_v

    IMPLICIT NONE


    REAL(wp) ::  c_0                !< constant used for diffusion coefficient and dissipation (dependent on mode RANS/LES)
    REAL(wp) ::  c_1                !< model constant for RANS mode
    REAL(wp) ::  c_2                !< model constant for RANS mode
    REAL(wp) ::  c_3                !< model constant for RANS mode
    REAL(wp) ::  c_4                !< model constant for RANS mode
    REAL(wp) ::  l_max              !< maximum length scale for Blackadar mixing length
    REAL(wp) ::  dsig_e = 1.0_wp    !< factor to calculate Ke from Km (1/sigma_e)
    REAL(wp) ::  dsig_diss = 1.0_wp !< factor to calculate K_diss from Km (1/sigma_diss)

    REAL(wp), DIMENSION(0:4) :: rans_const_c = &       !< model constants for RANS mode (namelist param)
       (/ 0.55_wp, 1.44_wp, 1.92_wp, 1.44_wp, 0.0_wp /) !> default values fit for standard-tke-e closure

    REAL(wp), DIMENSION(2) :: rans_const_sigma = &     !< model constants for RANS mode, sigma values (sigma_e, sigma_diss) (namelist param)
       (/ 1.0_wp, 1.30_wp /)

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  l_black    !< mixing length according to Blackadar

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  l_grid     !< geometric mean of grid sizes dx, dy, dz

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  l_wall !< near-wall mixing length

!
!-- Public variables
    PUBLIC c_0, rans_const_c, rans_const_sigma

    SAVE

    PRIVATE
!
!-- Public subroutines
    PUBLIC                                                                     &
       tcm_boundary_conds,                                                     &
       tcm_check_parameters,                                                   &
       tcm_check_data_output,                                                  &
       tcm_define_netcdf_grid,                                                 &
       tcm_init_arrays,                                                        &
       tcm_init,                                                               &
       tcm_actions,                                                            &
       tcm_prognostic_equations,                                               &
       tcm_swap_timelevel,                                                     &
       tcm_3d_data_averaging,                                                  &
       tcm_data_output_2d,                                                     &
       tcm_data_output_3d,                                                     &
       tcm_diffusivities

!
!-- PALM interfaces:
!-- Boundary conditions for subgrid TKE and dissipation
    INTERFACE tcm_boundary_conds
       MODULE PROCEDURE tcm_boundary_conds
    END INTERFACE tcm_boundary_conds
!
!-- Input parameter checks to be done in check_parameters
    INTERFACE tcm_check_parameters
       MODULE PROCEDURE tcm_check_parameters
    END INTERFACE tcm_check_parameters

!
!-- Data output checks for 2D/3D data to be done in check_parameters
    INTERFACE tcm_check_data_output
       MODULE PROCEDURE tcm_check_data_output
    END INTERFACE tcm_check_data_output

!
!-- Definition of data output quantities
    INTERFACE tcm_define_netcdf_grid
       MODULE PROCEDURE tcm_define_netcdf_grid
    END INTERFACE tcm_define_netcdf_grid

!
!-- Initialization of arrays
    INTERFACE tcm_init_arrays
       MODULE PROCEDURE tcm_init_arrays
    END INTERFACE tcm_init_arrays

!
!-- Initialization actions
    INTERFACE tcm_init
       MODULE PROCEDURE tcm_init
    END INTERFACE tcm_init

!
!-- Location specific actions
    INTERFACE tcm_actions
       MODULE PROCEDURE tcm_actions
       MODULE PROCEDURE tcm_actions_ij
    END INTERFACE tcm_actions

!
!-- Prognostic equations for TKE and TKE dissipation rate
    INTERFACE tcm_prognostic_equations
       MODULE PROCEDURE tcm_prognostic_equations
       MODULE PROCEDURE tcm_prognostic_equations_ij
    END INTERFACE tcm_prognostic_equations

!
!-- Swapping of time levels (required for prognostic variables)
    INTERFACE tcm_swap_timelevel
       MODULE PROCEDURE tcm_swap_timelevel
    END INTERFACE tcm_swap_timelevel

!
!-- Averaging of 3D data for output
    INTERFACE tcm_3d_data_averaging
       MODULE PROCEDURE tcm_3d_data_averaging
    END INTERFACE tcm_3d_data_averaging

!
!-- Data output of 2D quantities
    INTERFACE tcm_data_output_2d
       MODULE PROCEDURE tcm_data_output_2d
    END INTERFACE tcm_data_output_2d

!
!-- Data output of 3D data
    INTERFACE tcm_data_output_3d
       MODULE PROCEDURE tcm_data_output_3d
    END INTERFACE tcm_data_output_3d

!
!-- Call tcm_diffusivities_default and tcm_diffusivities_dynamic
    INTERFACE tcm_diffusivities
       MODULE PROCEDURE tcm_diffusivities
    END INTERFACE tcm_diffusivities


 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check parameters routine for turbulence closure module.
!------------------------------------------------------------------------------!
 SUBROUTINE tcm_boundary_conds

    USE pmc_interface,                                                         &
        ONLY : rans_mode_parent

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !< grid index x direction
    INTEGER(iwp) ::  j  !< grid index y direction
    INTEGER(iwp) ::  k  !< grid index z direction
    INTEGER(iwp) ::  l  !< running index boundary type, for up- and downward-facing walls
    INTEGER(iwp) ::  m  !< running index surface elements
!
!-- Boundary conditions for TKE.
    IF ( .NOT. constant_diffusion )  THEN
!
!--    In LES mode, Neumann conditions with de/x_i=0 are assumed at solid walls.
!--    Note, only TKE is prognostic in this case and dissipation is only
!--    a diagnostic quantity.
       IF ( .NOT. rans_mode )  THEN
!
!--       Horizontal walls, upward- and downward-facing
          DO  l = 0, 1
             !$OMP PARALLEL DO PRIVATE( i, j, k )
             !$ACC PARALLEL LOOP PRIVATE(i, j, k) &
             !$ACC PRESENT(bc_h, e_p)
             DO  m = 1, bc_h(l)%ns
                i = bc_h(l)%i(m)
                j = bc_h(l)%j(m)
                k = bc_h(l)%k(m)
                e_p(k+bc_h(l)%koff,j,i) = e_p(k,j,i)
             ENDDO
          ENDDO
!
!--       Vertical walls
          DO  l = 0, 3
!
!--          Note concerning missing ACC directive for this loop: Even though
!--          the data structure bc_v is present, it may not contain any
!--          allocated arrays in the flat but also in a topography case,
!--          leading to a runtime error. Therefore, omit ACC directives
!--          for this loop, in contrast to the bc_h loop.
             !$OMP PARALLEL DO PRIVATE( i, j, k )
             DO  m = 1, bc_v(l)%ns
                i = bc_v(l)%i(m)
                j = bc_v(l)%j(m)
                k = bc_v(l)%k(m)
                e_p(k,j+bc_v(l)%joff,i+bc_v(l)%ioff) = e_p(k,j,i)
             ENDDO
          ENDDO
!
!--    In RANS mode, wall function is used as boundary condition for TKE
       ELSE
!
!--       Use wall function within constant-flux layer
!--       Note, grid points listed in bc_h are not included in any calculations in RANS mode and
!--       are therefore not set here.
!
!--       Upward-facing surfaces
!--       Default surfaces
          DO  m = 1, surf_def_h(0)%ns
             i = surf_def_h(0)%i(m)
             j = surf_def_h(0)%j(m)
             k = surf_def_h(0)%k(m)
             e_p(k,j,i) = ( surf_def_h(0)%us(m) / c_0 )**2
          ENDDO
!
!--       Natural surfaces
          DO  m = 1, surf_lsm_h%ns
             i = surf_lsm_h%i(m)
             j = surf_lsm_h%j(m)
             k = surf_lsm_h%k(m)
             e_p(k,j,i) = ( surf_lsm_h%us(m) / c_0 )**2
          ENDDO
!
!--       Urban surfaces
          DO  m = 1, surf_usm_h%ns
             i = surf_usm_h%i(m)
             j = surf_usm_h%j(m)
             k = surf_usm_h%k(m)
             e_p(k,j,i) = ( surf_usm_h%us(m) / c_0 )**2
          ENDDO
!
!--       Vertical surfaces
          DO  l = 0, 3
!
!--          Default surfaces
             DO  m = 1, surf_def_v(l)%ns
                i = surf_def_v(l)%i(m)
                j = surf_def_v(l)%j(m)
                k = surf_def_v(l)%k(m)
                e_p(k,j,i) = ( surf_def_v(l)%us(m) / c_0 )**2
             ENDDO
!
!--          Natural surfaces
             DO  m = 1, surf_lsm_v(l)%ns
                i = surf_lsm_v(l)%i(m)
                j = surf_lsm_v(l)%j(m)
                k = surf_lsm_v(l)%k(m)
                e_p(k,j,i) = ( surf_lsm_v(l)%us(m) / c_0 )**2
             ENDDO
!
!--          Urban surfaces
             DO  m = 1, surf_usm_v(l)%ns
                i = surf_usm_v(l)%i(m)
                j = surf_usm_v(l)%j(m)
                k = surf_usm_v(l)%k(m)
                e_p(k,j,i) = ( surf_usm_v(l)%us(m) / c_0 )**2
             ENDDO
          ENDDO
       ENDIF
!
!--    Set Neumann boundary condition for TKE at model top. Do this also
!--    in case of a nested run.
       !$ACC KERNELS PRESENT(e_p)
       e_p(nzt+1,:,:) = e_p(nzt,:,:)
       !$ACC END KERNELS
!
!--    Nesting case: if parent operates in RANS mode and child in LES mode,
!--    no TKE is transfered. This case, set Neumann conditions at lateral and
!--    top child boundaries.
!--    If not ( both either in RANS or in LES mode ), TKE boundary condition
!--    is treated in the nesting.
       If ( child_domain )  THEN
          IF ( rans_mode_parent  .AND.  .NOT. rans_mode )  THEN

             e_p(nzt+1,:,:) = e_p(nzt,:,:)
             IF ( bc_dirichlet_l )  e_p(:,:,nxl-1) = e_p(:,:,nxl)
             IF ( bc_dirichlet_r )  e_p(:,:,nxr+1) = e_p(:,:,nxr)
             IF ( bc_dirichlet_s )  e_p(:,nys-1,:) = e_p(:,nys,:)
             IF ( bc_dirichlet_n )  e_p(:,nyn+1,:) = e_p(:,nyn,:)

          ENDIF
       ENDIF
!
!--    At in- and outflow boundaries also set Neumann boundary conditions
!--    for the SGS-TKE. An exception is made for the child domain if
!--    both parent and child operate in RANS mode. This case no
!--    lateral Neumann boundary conditions will be set but Dirichlet
!--    conditions will be set in the nesting.
       IF ( .NOT. child_domain  .AND.  .NOT. rans_mode_parent  .AND.           &
            .NOT. rans_mode )  THEN
          IF ( bc_dirichlet_s  .OR.  bc_radiation_s )  THEN
             e_p(:,nys-1,:) = e_p(:,nys,:)
             IF ( rans_tke_e )  diss_p(:,nys-1,:) = diss_p(:,nys,:)
          ENDIF
          IF ( bc_dirichlet_n  .OR.  bc_radiation_n )  THEN
             e_p(:,nyn+1,:) = e_p(:,nyn,:)
             IF ( rans_tke_e )  diss_p(:,nyn+1,:) = diss_p(:,nyn,:)
          ENDIF
          IF ( bc_dirichlet_l  .OR.  bc_radiation_l )  THEN
             e_p(:,:,nxl-1) = e_p(:,:,nxl)
             IF ( rans_tke_e )  diss_p(:,nyn+1,:) = diss_p(:,nyn,:)
          ENDIF
          IF ( bc_dirichlet_r  .OR.  bc_radiation_r )  THEN
             e_p(:,:,nxr+1) = e_p(:,:,nxr)
             IF ( rans_tke_e )  diss_p(:,nyn+1,:) = diss_p(:,nyn,:)
          ENDIF
       ENDIF
    ENDIF

!
!-- Boundary conditions for TKE dissipation rate in RANS mode.
    IF ( rans_tke_e )  THEN
!
!--    Use wall function within constant-flux layer
!--    Upward-facing surfaces
!--    Default surfaces
       DO  m = 1, surf_def_h(0)%ns
          i = surf_def_h(0)%i(m)
          j = surf_def_h(0)%j(m)
          k = surf_def_h(0)%k(m)
          diss_p(k,j,i) = surf_def_h(0)%us(m)**3          &
                        / ( kappa * surf_def_h(0)%z_mo(m) )
       ENDDO
!
!--    Natural surfaces
       DO  m = 1, surf_lsm_h%ns
          i = surf_lsm_h%i(m)
          j = surf_lsm_h%j(m)
          k = surf_lsm_h%k(m)
          diss_p(k,j,i) = surf_lsm_h%us(m)**3          &
                        / ( kappa * surf_lsm_h%z_mo(m) )
       ENDDO
!
!--    Urban surfaces
       DO  m = 1, surf_usm_h%ns
          i = surf_usm_h%i(m)
          j = surf_usm_h%j(m)
          k = surf_usm_h%k(m)
          diss_p(k,j,i) = surf_usm_h%us(m)**3          &
                        / ( kappa * surf_usm_h%z_mo(m) )
       ENDDO
!
!--    Vertical surfaces
       DO  l = 0, 3
!
!--       Default surfaces
          DO  m = 1, surf_def_v(l)%ns
             i = surf_def_v(l)%i(m)
             j = surf_def_v(l)%j(m)
             k = surf_def_v(l)%k(m)
             diss_p(k,j,i) = surf_def_v(l)%us(m)**3          &
                           / ( kappa * surf_def_v(l)%z_mo(m) )
          ENDDO
!
!--       Natural surfaces
          DO  m = 1, surf_lsm_v(l)%ns
             i = surf_lsm_v(l)%i(m)
             j = surf_lsm_v(l)%j(m)
             k = surf_lsm_v(l)%k(m)
             diss_p(k,j,i) = surf_lsm_v(l)%us(m)**3          &
                           / ( kappa * surf_lsm_v(l)%z_mo(m) )
          ENDDO
!
!--       Urban surfaces
          DO  m = 1, surf_usm_v(l)%ns
             i = surf_usm_v(l)%i(m)
             j = surf_usm_v(l)%j(m)
             k = surf_usm_v(l)%k(m)
             diss_p(k,j,i) = surf_usm_v(l)%us(m)**3          &
                           / ( kappa * surf_usm_v(l)%z_mo(m) )
          ENDDO
       ENDDO
!
!--    Limit change of diss to be between -90% and +100%. Also, set an absolute
!--    minimum value
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb, nzt+1
                diss_p(k,j,i) = MAX( MIN( diss_p(k,j,i),          &
                                          2.0_wp * diss(k,j,i) ), &
                                     0.1_wp * diss(k,j,i),        &
                                     0.0001_wp )
             ENDDO
          ENDDO
       ENDDO

       diss_p(nzt+1,:,:) = diss_p(nzt,:,:)

    ENDIF

 END SUBROUTINE tcm_boundary_conds

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check parameters routine for turbulence closure module.
!------------------------------------------------------------------------------!
 SUBROUTINE tcm_check_parameters

    USE control_parameters,                                                    &
        ONLY:  message_string, turbulent_inflow, turbulent_outflow

    IMPLICIT NONE

!
!-- Define which turbulence closure is going to be used
    SELECT CASE ( TRIM( turbulence_closure ) )

       CASE ( 'dynamic' )
          les_dynamic = .TRUE.

       CASE ( 'Moeng_Wyngaard' )
          les_mw = .TRUE.

       CASE ( 'TKE-l' )
          rans_tke_l = .TRUE.
          rans_mode = .TRUE.

       CASE ( 'TKE-e' )
          rans_tke_e = .TRUE.
          rans_mode = .TRUE.

       CASE DEFAULT
          message_string = 'Unknown turbulence closure: ' //                &
                           TRIM( turbulence_closure )
          CALL message( 'tcm_check_parameters', 'PA0500', 1, 2, 0, 6, 0 )

    END SELECT
!
!-- Set variables for RANS mode or LES mode
    IF ( rans_mode )  THEN
!
!--    Assign values to constants for RANS mode
       dsig_e    = 1.0_wp / rans_const_sigma(1)
       dsig_diss = 1.0_wp / rans_const_sigma(2)

       c_0 = rans_const_c(0)
       c_1 = rans_const_c(1)
       c_2 = rans_const_c(2)
       c_3 = rans_const_c(3)   !> @todo clarify how to switch between different models
       c_4 = rans_const_c(4)

       IF ( turbulent_inflow .OR. turbulent_outflow )  THEN
          message_string = 'turbulent inflow/outflow is not yet '//            &
                           'implemented for RANS mode'
          CALL message( 'tcm_check_parameters', 'PA0501', 1, 2, 0, 6, 0 )
       ENDIF

    ELSE
!
!--    LES mode
       c_0 = 0.1_wp    !according to Lilly (1967) and Deardorff (1980)

       dsig_e = 1.0_wp !assure to use K_m to calculate TKE instead
                       !of K_e which is used in RANS mode

    ENDIF

 END SUBROUTINE tcm_check_parameters

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check data output.
!------------------------------------------------------------------------------!
 SUBROUTINE tcm_check_data_output( var, unit )

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  unit     !< unit of output variable
    CHARACTER (LEN=*) ::  var      !< name of output variable


    SELECT CASE ( TRIM( var ) )

       CASE ( 'diss' )
          unit = 'm2/s3'

       CASE ( 'kh', 'km' )
          unit = 'm2/s'

       CASE DEFAULT
          unit = 'illegal'

    END SELECT

 END SUBROUTINE tcm_check_data_output


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Define appropriate grid for netcdf variables.
!> It is called out from subroutine netcdf.
!------------------------------------------------------------------------------!
 SUBROUTINE tcm_define_netcdf_grid( var, found, grid_x, grid_y, grid_z )

    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT(OUT) ::  grid_x   !< x grid of output variable
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_y   !< y grid of output variable
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_z   !< z grid of output variable
    CHARACTER (LEN=*), INTENT(IN)  ::  var      !< name of output variable

    LOGICAL, INTENT(OUT) ::  found   !< flag if output variable is found

    found  = .TRUE.

!
!-- Check for the grid
    SELECT CASE ( TRIM( var ) )

       CASE ( 'diss', 'diss_xy', 'diss_xz', 'diss_yz' )
          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zu'

       CASE ( 'kh', 'kh_xy', 'kh_xz', 'kh_yz' )
          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zu'

       CASE ( 'km', 'km_xy', 'km_xz', 'km_yz' )
          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zu'

       CASE DEFAULT
          found  = .FALSE.
          grid_x = 'none'
          grid_y = 'none'
          grid_z = 'none'

    END SELECT

 END SUBROUTINE tcm_define_netcdf_grid


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Average 3D data.
!------------------------------------------------------------------------------!
 SUBROUTINE tcm_3d_data_averaging( mode, variable )


    USE averaging,                                                             &
        ONLY:  diss_av, kh_av, km_av

    USE control_parameters,                                                    &
        ONLY:  average_count_3d

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  mode       !< flag defining mode 'allocate', 'sum' or 'average'
    CHARACTER (LEN=*) ::  variable   !< name of variable

    INTEGER(iwp) ::  i   !< loop index
    INTEGER(iwp) ::  j   !< loop index
    INTEGER(iwp) ::  k   !< loop index

    IF ( mode == 'allocate' )  THEN

       SELECT CASE ( TRIM( variable ) )

          CASE ( 'diss' )
             IF ( .NOT. ALLOCATED( diss_av ) )  THEN
                ALLOCATE( diss_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             diss_av = 0.0_wp

          CASE ( 'kh' )
             IF ( .NOT. ALLOCATED( kh_av ) )  THEN
                ALLOCATE( kh_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             kh_av = 0.0_wp

          CASE ( 'km' )
             IF ( .NOT. ALLOCATED( km_av ) )  THEN
                ALLOCATE( km_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             km_av = 0.0_wp

          CASE DEFAULT
             CONTINUE

       END SELECT

    ELSEIF ( mode == 'sum' )  THEN

       SELECT CASE ( TRIM( variable ) )

          CASE ( 'diss' )
             IF ( ALLOCATED( diss_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         diss_av(k,j,i) = diss_av(k,j,i) + diss(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'kh' )
             IF ( ALLOCATED( kh_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         kh_av(k,j,i) = kh_av(k,j,i) + kh(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'km' )
             IF ( ALLOCATED( km_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         km_av(k,j,i) = km_av(k,j,i) + km(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE DEFAULT
             CONTINUE

       END SELECT

    ELSEIF ( mode == 'average' )  THEN

       SELECT CASE ( TRIM( variable ) )

          CASE ( 'diss' )
             IF ( ALLOCATED( diss_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         diss_av(k,j,i) = diss_av(k,j,i)                       &
                                        / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'kh' )
             IF ( ALLOCATED( kh_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         kh_av(k,j,i) = kh_av(k,j,i)                           &
                                        / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'km' )
             IF ( ALLOCATED( km_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         km_av(k,j,i) = km_av(k,j,i)                           &
                                        / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

       END SELECT

    ENDIF

 END SUBROUTINE tcm_3d_data_averaging


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Define 2D output variables.
!------------------------------------------------------------------------------!
 SUBROUTINE tcm_data_output_2d( av, variable, found, grid, mode, local_pf,     &
                                nzb_do, nzt_do )

    USE averaging,                                                             &
        ONLY:  diss_av, kh_av, km_av

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  grid       !< name of vertical grid
    CHARACTER (LEN=*) ::  mode       !< either 'xy', 'xz' or 'yz'
    CHARACTER (LEN=*) ::  variable   !< name of variable

    INTEGER(iwp) ::  av        !< flag for (non-)average output
    INTEGER(iwp) ::  flag_nr   !< number of masking flag
    INTEGER(iwp) ::  i         !< loop index
    INTEGER(iwp) ::  j         !< loop index
    INTEGER(iwp) ::  k         !< loop index
    INTEGER(iwp) ::  nzb_do    !< vertical output index (bottom)
    INTEGER(iwp) ::  nzt_do    !< vertical output index (top)

    LOGICAL ::  found     !< flag if output variable is found
    LOGICAL ::  resorted  !< flag if output is already resorted

    REAL(wp) ::  fill_value = -9999.0_wp  !< value for the _FillValue attribute

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf !< local
       !< array to which output data is resorted to

    REAL(wp), DIMENSION(:,:,:), POINTER ::  to_be_resorted  !< points to selected output variable

    found = .TRUE.
    resorted = .FALSE.
!
!-- Set masking flag for topography for not resorted arrays
    flag_nr = 0

    SELECT CASE ( TRIM( variable ) )

       CASE ( 'diss_xy', 'diss_xz', 'diss_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => diss
          ELSE
             IF ( .NOT. ALLOCATED( diss_av ) ) THEN
                ALLOCATE( diss_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                diss_av = REAL( fill_value, KIND = wp )
             ENDIF
             to_be_resorted => diss_av
          ENDIF
          IF ( mode == 'xy' ) grid = 'zu'

       CASE ( 'kh_xy', 'kh_xz', 'kh_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => kh
          ELSE
             IF ( .NOT. ALLOCATED( kh_av ) ) THEN
                ALLOCATE( kh_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                kh_av = REAL( fill_value, KIND = wp )
             ENDIF
             to_be_resorted => kh_av
          ENDIF
          IF ( mode == 'xy' ) grid = 'zu'

       CASE ( 'km_xy', 'km_xz', 'km_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => km
          ELSE
             IF ( .NOT. ALLOCATED( km_av ) ) THEN
                ALLOCATE( km_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                km_av = REAL( fill_value, KIND = wp )
             ENDIF
             to_be_resorted => km_av
          ENDIF
          IF ( mode == 'xy' ) grid = 'zu'

       CASE DEFAULT
          found = .FALSE.
          grid  = 'none'

    END SELECT

    IF ( found .AND. .NOT. resorted )  THEN
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb_do, nzt_do
                local_pf(i,j,k) = MERGE( to_be_resorted(k,j,i),                &
                                  REAL( fill_value, KIND = wp ),               &
                                  BTEST( wall_flags_total_0(k,j,i), flag_nr ) )
             ENDDO
          ENDDO
       ENDDO
    ENDIF

 END SUBROUTINE tcm_data_output_2d


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Define 3D output variables.
!------------------------------------------------------------------------------!
 SUBROUTINE tcm_data_output_3d( av, variable, found, local_pf, nzb_do, nzt_do )


    USE averaging,                                                             &
        ONLY:  diss_av, kh_av, km_av

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  variable   !< name of variable

    INTEGER(iwp) ::  av        !< flag for (non-)average output
    INTEGER(iwp) ::  flag_nr   !< number of masking flag
    INTEGER(iwp) ::  i         !< loop index
    INTEGER(iwp) ::  j         !< loop index
    INTEGER(iwp) ::  k         !< loop index
    INTEGER(iwp) ::  nzb_do    !< lower limit of the data output (usually 0)
    INTEGER(iwp) ::  nzt_do    !< vertical upper limit of the data output (usually nz_do3d)

    LOGICAL ::  found     !< flag if output variable is found
    LOGICAL ::  resorted  !< flag if output is already resorted

    REAL(wp) ::  fill_value = -9999.0_wp  !< value for the _FillValue attribute

    REAL(sp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf   !< local
       !< array to which output data is resorted to

    REAL(wp), DIMENSION(:,:,:), POINTER ::  to_be_resorted  !< points to selected output variable

    found = .TRUE.
    resorted = .FALSE.
!
!-- Set masking flag for topography for not resorted arrays
    flag_nr = 0

    SELECT CASE ( TRIM( variable ) )

       CASE ( 'diss' )
          IF ( av == 0 )  THEN
             to_be_resorted => diss
          ELSE
             IF ( .NOT. ALLOCATED( diss_av ) ) THEN
                ALLOCATE( diss_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                diss_av = REAL( fill_value, KIND = wp )
             ENDIF
             to_be_resorted => diss_av
          ENDIF

       CASE ( 'kh' )
          IF ( av == 0 )  THEN
             to_be_resorted => kh
          ELSE
             IF ( .NOT. ALLOCATED( kh_av ) ) THEN
                ALLOCATE( kh_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                kh_av = REAL( fill_value, KIND = wp )
             ENDIF
             to_be_resorted => kh_av
          ENDIF

       CASE ( 'km' )
          IF ( av == 0 )  THEN
             to_be_resorted => km
          ELSE
             IF ( .NOT. ALLOCATED( km_av ) ) THEN
                ALLOCATE( km_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                km_av = REAL( fill_value, KIND = wp )
             ENDIF
             to_be_resorted => km_av
          ENDIF

       CASE DEFAULT
          found = .FALSE.

    END SELECT


    IF ( found .AND. .NOT. resorted )  THEN
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb_do, nzt_do
                local_pf(i,j,k) = MERGE(                                 &
                                   to_be_resorted(k,j,i),                &
                                   REAL( fill_value, KIND = wp ),        &
                                   BTEST( wall_flags_total_0(k,j,i), flag_nr ) )
             ENDDO
          ENDDO
       ENDDO
       resorted = .TRUE.
    ENDIF

 END SUBROUTINE tcm_data_output_3d


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate arrays and assign pointers.
!------------------------------------------------------------------------------!
 SUBROUTINE tcm_init_arrays

    USE bulk_cloud_model_mod,                                                  &
        ONLY:  collision_turbulence

    USE pmc_interface,                                                         &
        ONLY:  nested_run

    IMPLICIT NONE

    ALLOCATE( kh(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( km(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

    ALLOCATE( e_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( e_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( e_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

!
!-- Allocate arrays required for dissipation.
!-- Please note, if it is a nested run, arrays need to be allocated even if
!-- they do not necessarily need to be transferred, which is attributed to
!-- the design of the model coupler which allocates memory for each variable.
    ALLOCATE( diss_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

    IF ( rans_tke_e  .OR.  nested_run )  THEN
       ALLOCATE( diss_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       ALLOCATE( diss_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ENDIF

!
!-- Initial assignment of pointers
    e  => e_1;   e_p  => e_2;   te_m  => e_3

    diss => diss_1
    IF ( rans_tke_e  .OR.  nested_run )  THEN
       diss_p => diss_2; tdiss_m => diss_3
    ENDIF

 END SUBROUTINE tcm_init_arrays


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of turbulence closure module.
!------------------------------------------------------------------------------!
 SUBROUTINE tcm_init

    USE control_parameters,                                                    &
        ONLY:  bc_dirichlet_l, complex_terrain, topography

    USE model_1d_mod,                                                          &
        ONLY:  e1d, kh1d, km1d

    IMPLICIT NONE

    INTEGER(iwp) :: i            !< loop index
    INTEGER(iwp) :: j            !< loop index
    INTEGER(iwp) :: k            !< loop index
    INTEGER(iwp) :: nz_s_shift   !< lower shift index for scalars
    INTEGER(iwp) :: nz_s_shift_l !< local lower shift index in case of turbulent inflow

!
!-- Initialize mixing length
    CALL tcm_init_mixing_length

!
!-- Actions for initial runs
    IF ( TRIM( initializing_actions ) /= 'read_restart_data'  .AND.            &
         TRIM( initializing_actions ) /= 'cyclic_fill' )  THEN

       IF ( INDEX( initializing_actions, 'set_1d-model_profiles' ) /= 0 )  THEN

          IF ( .NOT. rans_tke_e ) THEN
!
!--          Transfer initial profiles to the arrays of the 3D model
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   e(:,j,i)  = e1d
                   kh(:,j,i) = kh1d
                   km(:,j,i) = km1d
                ENDDO
             ENDDO

             IF ( constant_diffusion )  THEN
                e = 0.0_wp
             ENDIF

             diss = 0.0_wp

          ELSE
!
!--          In case of TKE-e closure in RANS mode, do not use e, diss, and km
!--          profiles from 1D model. Instead, initialize with constant profiles
             IF ( constant_diffusion )  THEN
                km = km_constant
                kh = km / prandtl_number
                e  = 0.0_wp
             ELSEIF ( e_init > 0.0_wp )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb+1, nzt
                         km(k,j,i) = c_0 * l_wall(k,j,i) * SQRT( e_init )
                      ENDDO
                   ENDDO
                ENDDO
                km(nzb,:,:)   = km(nzb+1,:,:)
                km(nzt+1,:,:) = km(nzt,:,:)
                kh = km / prandtl_number
                e  = e_init
             ELSE
                IF ( .NOT. ocean_mode )  THEN
                   kh   = 0.01_wp   ! there must exist an initial diffusion, because
                   km   = 0.01_wp   ! otherwise no TKE would be produced by the
                                    ! production terms, as long as not yet
                                    ! e = (u*/cm)**2 at k=nzb+1
                ELSE
                   kh   = 0.00001_wp
                   km   = 0.00001_wp
                ENDIF
                e    = 0.0_wp
             ENDIF

             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb+1, nzt
                      diss(k,j,i) = c_0**4 * e(k,j,i)**2 / km(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
             diss(nzb,:,:) = diss(nzb+1,:,:)
             diss(nzt+1,:,:) = diss(nzt,:,:)

          ENDIF

       ELSEIF ( INDEX( initializing_actions, 'set_constant_profiles' ) /= 0 .OR. &
                INDEX( initializing_actions, 'inifor' ) /= 0 )  THEN

          IF ( constant_diffusion )  THEN
             km = km_constant
             kh = km / prandtl_number
             e  = 0.0_wp
          ELSEIF ( e_init > 0.0_wp )  THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb+1, nzt
                      km(k,j,i) = c_0 * l_wall(k,j,i) * SQRT( e_init )
                   ENDDO
                ENDDO
             ENDDO
             km(nzb,:,:)   = km(nzb+1,:,:)
             km(nzt+1,:,:) = km(nzt,:,:)
             kh = km / prandtl_number
             e  = e_init
          ELSE
             IF ( .NOT. ocean_mode )  THEN
                kh   = 0.01_wp   ! there must exist an initial diffusion, because
                km   = 0.01_wp   ! otherwise no TKE would be produced by the
                                 ! production terms, as long as not yet
                                 ! e = (u*/cm)**2 at k=nzb+1
             ELSE
                kh   = 0.00001_wp
                km   = 0.00001_wp
             ENDIF
             e    = 0.0_wp
          ENDIF

          IF ( rans_tke_e )  THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb+1, nzt
                      diss(k,j,i) = c_0**4 * e(k,j,i)**2 / km(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
             diss(nzb,:,:) = diss(nzb+1,:,:)
             diss(nzt+1,:,:) = diss(nzt,:,:)
          ELSE
             diss = 0.0_wp
          ENDIF

       ENDIF
!
!--    Store initial profiles for output purposes etc.
       hom(:,1,23,:) = SPREAD( km(:,nys,nxl), 2, statistic_regions+1 )
       hom(:,1,24,:) = SPREAD( kh(:,nys,nxl), 2, statistic_regions+1 )
!
!--    Initialize old and new time levels.
       te_m = 0.0_wp
       e_p = e
       IF ( rans_tke_e )  THEN
          tdiss_m = 0.0_wp
          diss_p = diss
       ENDIF

    ELSEIF ( TRIM( initializing_actions ) == 'read_restart_data'  .OR.         &
             TRIM( initializing_actions ) == 'cyclic_fill' )                   &
    THEN

!
!--    In case of complex terrain and cyclic fill method as initialization,
!--    shift initial data in the vertical direction for each point in the
!--    x-y-plane depending on local surface height
       IF ( complex_terrain  .AND.                                             &
            TRIM( initializing_actions ) == 'cyclic_fill' )  THEN
          DO  i = nxlg, nxrg
             DO  j = nysg, nyng
                nz_s_shift = topo_top_ind(j,i,0)

                e(nz_s_shift:nzt+1,j,i)  =  e(0:nzt+1-nz_s_shift,j,i)
                km(nz_s_shift:nzt+1,j,i) = km(0:nzt+1-nz_s_shift,j,i)
                kh(nz_s_shift:nzt+1,j,i) = kh(0:nzt+1-nz_s_shift,j,i)
             ENDDO
          ENDDO
          IF ( rans_tke_e )  THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   nz_s_shift = topo_top_ind(j,i,0)

                   diss(nz_s_shift:nzt+1,j,i) = diss(0:nzt+1-nz_s_shift,j,i)
                ENDDO
             ENDDO
          ELSE
             diss = 0.0_wp
          ENDIF
       ENDIF

!
!--    Initialization of the turbulence recycling method
       IF ( TRIM( initializing_actions ) == 'cyclic_fill'  .AND.               &
            turbulent_inflow )  THEN
          mean_inflow_profiles(:,5) = hom_sum(:,8,0)   ! e
!
!--       In case of complex terrain, determine vertical displacement at inflow
!--       boundary and adjust mean inflow profiles
          IF ( complex_terrain )  THEN
             IF ( nxlg <= 0 .AND. nxrg >= 0 .AND.  &
                  nysg <= 0 .AND. nyng >= 0        )  THEN
                nz_s_shift_l = topo_top_ind(0,0,0)
             ELSE
                nz_s_shift_l = 0
             ENDIF
#if defined( __parallel )
             CALL MPI_ALLREDUCE(nz_s_shift_l, nz_s_shift, 1, MPI_INTEGER,      &
                                MPI_MAX, comm2d, ierr)
#else
             nz_s_shift = nz_s_shift_l
#endif
             mean_inflow_profiles(nz_s_shift:nzt+1,5) =  &
                hom_sum(0:nzt+1-nz_s_shift,8,0)  ! e
          ENDIF
!
!--       Use these mean profiles at the inflow (provided that Dirichlet
!--       conditions are used)
          IF ( bc_dirichlet_l )  THEN
             DO  j = nysg, nyng
                DO  k = nzb, nzt+1
                   e(k,j,nxlg:-1)  = mean_inflow_profiles(k,5)
                ENDDO
             ENDDO
          ENDIF
       ENDIF
!
!--    Inside buildings set TKE back to zero
       IF ( TRIM( initializing_actions ) == 'cyclic_fill' .AND.                &
            topography /= 'flat' )  THEN
!
!--       Inside buildings set TKE back to zero.
!--       Other scalars (km, kh,...) are ignored at present,
!--       maybe revise later.
          DO  i = nxlg, nxrg
             DO  j = nysg, nyng
                DO  k = nzb, nzt
                   e(k,j,i)     = MERGE( e(k,j,i), 0.0_wp,                     &
                                         BTEST( wall_flags_total_0(k,j,i), 0 ) )
                ENDDO
             ENDDO
          ENDDO

          IF ( rans_tke_e )  THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt
                      diss(k,j,i)    = MERGE( diss(k,j,i), 0.0_wp,             &
                                         BTEST( wall_flags_total_0(k,j,i), 0 ) )
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDIF
!
!--    Initialize new time levels (only done in order to set boundary values
!--    including ghost points)
       e_p = e
!
!--    Allthough tendency arrays are set in prognostic_equations, they have
!--    to be predefined here because there they are used (but multiplied with 0)
!--    before they are set.
       te_m = 0.0_wp

       IF ( rans_tke_e )  THEN
          diss_p = diss
          tdiss_m = 0.0_wp
       ENDIF

    ENDIF

 END SUBROUTINE tcm_init


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Pre-computation of grid-dependent and near-wall mixing length.
!> @todo consider walls in horizontal direction at a distance further than a
!>       single grid point (RANS mode)
!------------------------------------------------------------------------------!
 SUBROUTINE tcm_init_mixing_length

    USE arrays_3d,                                                             &
        ONLY:  dzw, ug, vg, zu, zw

    USE control_parameters,                                                    &
        ONLY:  f, message_string, wall_adjustment, wall_adjustment_factor

    USE exchange_horiz_mod,                                                    &
        ONLY:  exchange_horiz, exchange_horiz_int

    USE grid_variables,                                                        &
        ONLY:  dx, dy

    USE indices,                                                               &
        ONLY:  nbgp, nx, nxl, nxlg, nxr, nxrg, ny, nyn, nyng, nys, nysg, nzb,  &
               nzt, wall_flags_total_0

    USE kinds


    IMPLICIT NONE

    INTEGER(iwp) :: dist_dx     !< found distance devided by dx
    INTEGER(iwp) :: i           !< index variable along x
    INTEGER(iwp) :: ii          !< index variable along x
    INTEGER(iwp) :: j           !< index variable along y
    INTEGER(iwp) :: k           !< index variable along z
    INTEGER(iwp) :: k_max_topo  !< index of maximum topography height
    INTEGER(iwp) :: kk          !< index variable along z
    INTEGER(iwp) :: rad_i       !< search radius in grid points along x
    INTEGER(iwp) :: rad_i_l     !< possible search radius to the left
    INTEGER(iwp) :: rad_i_r     !< possible search radius to the right
    INTEGER(iwp) :: rad_j       !< search radius in grid points along y
    INTEGER(iwp) :: rad_j_n     !< possible search radius to north
    INTEGER(iwp) :: rad_j_s     !< possible search radius to south
    INTEGER(iwp) :: rad_k       !< search radius in grid points along z
    INTEGER(iwp) :: rad_k_b     !< search radius in grid points along negative z
    INTEGER(iwp) :: rad_k_t     !< search radius in grid points along positive z

    INTEGER(KIND=1), DIMENSION(:,:), ALLOCATABLE :: vic_yz !< contains a quarter of a single yz-slice of vicinity

    INTEGER(KIND=1), DIMENSION(:,:,:), ALLOCATABLE :: vicinity !< contains topography information of the vicinity of (i/j/k)

    REAL(wp) :: radius          !< search radius in meter

    ALLOCATE( l_grid(1:nzt) )
    ALLOCATE( l_wall(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
!
!-- Initialize the mixing length in case of an LES-simulation
    IF ( .NOT. rans_mode )  THEN
!
!--    Compute the grid-dependent mixing length.
       DO  k = 1, nzt
          l_grid(k)  = ( dx * dy * dzw(k) )**0.33333333333333_wp
       ENDDO
!
!--    Initialize near-wall mixing length l_wall only in the vertical direction
!--    for the moment, multiplication with wall_adjustment_factor further below
       l_wall(nzb,:,:)   = l_grid(1)
       DO  k = nzb+1, nzt
          l_wall(k,:,:)  = l_grid(k)
       ENDDO
       l_wall(nzt+1,:,:) = l_grid(nzt)

       IF ( wall_adjustment )  THEN

          DO  k = 1, nzt
             IF ( l_grid(k) > 1.5_wp * dx * wall_adjustment_factor .OR.            &
                  l_grid(k) > 1.5_wp * dy * wall_adjustment_factor )  THEN
                WRITE( message_string, * ) 'grid anisotropy exceeds ',             &
                                           'threshold given by only local',        &
                                           ' &horizontal reduction of near_wall ', &
                                           'mixing length l_wall',                 &
                                           ' &starting from height level k = ', k, &
                                           '.'
                CALL message( 'init_grid', 'PA0202', 0, 1, 0, 6, 0 )
                EXIT
             ENDIF
          ENDDO
!
!--       In case of topography: limit near-wall mixing length l_wall further:
!--       Go through all points of the subdomain one by one and look for the closest
!--       surface.
!--       Is this correct in the ocean case?
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
!
!--                Check if current gridpoint belongs to the atmosphere
                   IF ( BTEST( wall_flags_total_0(k,j,i), 0 ) )  THEN
!
!--                   Check for neighbouring grid-points.
!--                   Vertical distance, down
                      IF ( .NOT. BTEST( wall_flags_total_0(k-1,j,i), 0 ) )             &
                         l_wall(k,j,i) = MIN( l_grid(k), zu(k) - zw(k-1) )
!
!--                   Vertical distance, up
                      IF ( .NOT. BTEST( wall_flags_total_0(k+1,j,i), 0 ) )             &
                         l_wall(k,j,i) = MIN( l_wall(k,j,i), l_grid(k), zw(k) - zu(k) )
!
!--                   y-distance
                      IF ( .NOT. BTEST( wall_flags_total_0(k,j-1,i), 0 )  .OR.         &
                           .NOT. BTEST( wall_flags_total_0(k,j+1,i), 0 ) )             &
                         l_wall(k,j,i) = MIN( l_wall(k,j,i), l_grid(k), 0.5_wp * dy )
!
!--                   x-distance
                      IF ( .NOT. BTEST( wall_flags_total_0(k,j,i-1), 0 )  .OR.         &
                           .NOT. BTEST( wall_flags_total_0(k,j,i+1), 0 ) )             &
                         l_wall(k,j,i) = MIN( l_wall(k,j,i), l_grid(k), 0.5_wp * dx )
!
!--                   yz-distance (vertical edges, down)
                      IF ( .NOT. BTEST( wall_flags_total_0(k-1,j-1,i), 0 )  .OR.       &
                           .NOT. BTEST( wall_flags_total_0(k-1,j+1,i), 0 )  )          &
                         l_wall(k,j,i) = MIN( l_wall(k,j,i), l_grid(k),                &
                                              SQRT( 0.25_wp * dy**2 +                  &
                                             ( zu(k) - zw(k-1) )**2 ) )
!
!--                   yz-distance (vertical edges, up)
                      IF ( .NOT. BTEST( wall_flags_total_0(k+1,j-1,i), 0 )  .OR.       &
                           .NOT. BTEST( wall_flags_total_0(k+1,j+1,i), 0 )  )          &
                         l_wall(k,j,i) = MIN( l_wall(k,j,i), l_grid(k),                &
                                              SQRT( 0.25_wp * dy**2 +                  &
                                             ( zw(k) - zu(k) )**2 ) )
!
!--                   xz-distance (vertical edges, down)
                      IF ( .NOT. BTEST( wall_flags_total_0(k-1,j,i-1), 0 )  .OR.       &
                           .NOT. BTEST( wall_flags_total_0(k-1,j,i+1), 0 )  )          &
                         l_wall(k,j,i) = MIN( l_wall(k,j,i), l_grid(k),                &
                                              SQRT( 0.25_wp * dx**2 +                  &
                                             ( zu(k) - zw(k-1) )**2 ) )
!
!--                   xz-distance (vertical edges, up)
                      IF ( .NOT. BTEST( wall_flags_total_0(k+1,j,i-1), 0 )  .OR.       &
                           .NOT. BTEST( wall_flags_total_0(k+1,j,i+1), 0 )  )          &
                         l_wall(k,j,i) = MIN( l_wall(k,j,i), l_grid(k),                &
                                              SQRT( 0.25_wp * dx**2 +                  &
                                             ( zw(k) - zu(k) )**2 ) )
!
!--                   xy-distance (horizontal edges)
                      IF ( .NOT. BTEST( wall_flags_total_0(k,j-1,i-1), 0 )  .OR.        &
                           .NOT. BTEST( wall_flags_total_0(k,j+1,i-1), 0 )  .OR.        &
                           .NOT. BTEST( wall_flags_total_0(k,j-1,i+1), 0 )  .OR.        &
                           .NOT. BTEST( wall_flags_total_0(k,j+1,i+1), 0 ) )            &
                         l_wall(k,j,i) = MIN( l_wall(k,j,i), l_grid(k),                 &
                                              SQRT( 0.25_wp * ( dx**2 + dy**2 ) ) )
!
!--                   xyz distance (vertical and horizontal edges, down)
                      IF ( .NOT. BTEST( wall_flags_total_0(k-1,j-1,i-1), 0 )  .OR.      &
                           .NOT. BTEST( wall_flags_total_0(k-1,j+1,i-1), 0 )  .OR.      &
                           .NOT. BTEST( wall_flags_total_0(k-1,j-1,i+1), 0 )  .OR.      &
                           .NOT. BTEST( wall_flags_total_0(k-1,j+1,i+1), 0 ) )          &
                         l_wall(k,j,i) = MIN( l_wall(k,j,i), l_grid(k),                 &
                                              SQRT( 0.25_wp * ( dx**2 + dy**2 )         &
                                                    +  ( zu(k) - zw(k-1) )**2  ) )
!
!--                   xyz distance (vertical and horizontal edges, up)
                      IF ( .NOT. BTEST( wall_flags_total_0(k+1,j-1,i-1), 0 )  .OR.      &
                           .NOT. BTEST( wall_flags_total_0(k+1,j+1,i-1), 0 )  .OR.      &
                           .NOT. BTEST( wall_flags_total_0(k+1,j-1,i+1), 0 )  .OR.      &
                           .NOT. BTEST( wall_flags_total_0(k+1,j+1,i+1), 0 ) )          &
                         l_wall(k,j,i) = MIN( l_wall(k,j,i), l_grid(k),                 &
                                              SQRT( 0.25_wp * ( dx**2 + dy**2 )         &
                                                    +  ( zw(k) - zu(k) )**2  ) )

                   ENDIF
!
!--                Adjust mixing length by wall-adjustment factor and limit it by l_grid
                   l_wall(k,j,i) = MIN( l_wall(k,j,i) * wall_adjustment_factor, l_grid(k) )

                ENDDO  !k loop
             ENDDO  !j loop
          ENDDO  !i loop

       ENDIF  !if wall_adjustment

    ELSE
!
!--    Initialize the mixing length in case of a RANS simulation
       ALLOCATE( l_black(nzb:nzt+1) )

!
!--    Calculate mixing length according to Blackadar (1962)
       IF ( f /= 0.0_wp )  THEN
          l_max = 2.7E-4_wp * SQRT( ug(nzt+1)**2 + vg(nzt+1)**2 ) /            &
                  ABS( f ) + 1.0E-10_wp
       ELSE
          l_max = 30.0_wp
       ENDIF

       DO  k = nzb, nzt
          l_black(k) = kappa * zu(k) / ( 1.0_wp + kappa * zu(k) / l_max )
       ENDDO

       l_black(nzt+1) = l_black(nzt)

!
!--    Get height level of highest topography within local subdomain
       k_max_topo = 0
       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
             DO  k = nzb+1, nzt-1
                IF ( .NOT. BTEST( wall_flags_total_0(k,j,i), 0 ) .AND.  &
                     k > k_max_topo )  &
                   k_max_topo = k
             ENDDO
          ENDDO
       ENDDO

       l_wall(nzb,:,:) = l_black(nzb)
       l_wall(nzt+1,:,:) = l_black(nzt+1)
!
!--    Limit mixing length to either nearest wall or Blackadar mixing length.
!--    For that, analyze each grid point (i/j/k) ("analysed grid point") and
!--    search within its vicinity for the shortest distance to a wall by cal-
!--    culating the distance between the analysed grid point and the "viewed
!--    grid point" if it contains a wall (belongs to topography).
       DO  k = nzb+1, nzt

          radius = l_black(k)  ! radius within walls are searched
!
!--       Set l_wall to its default maximum value (l_back)
          l_wall(k,:,:) = radius

!
!--       Compute search radius as number of grid points in all directions
          rad_i = CEILING( radius / dx )
          rad_j = CEILING( radius / dy )

          DO  kk = 0, nzt-k
             rad_k_t = kk
!
!--          Limit upward search radius to height of maximum topography
             IF ( zu(k+kk)-zu(k) >= radius .OR. k+kk >= k_max_topo )  EXIT
          ENDDO

          DO  kk = 0, k
             rad_k_b = kk
             IF ( zu(k)-zu(k-kk) >= radius )  EXIT
          ENDDO

!
!--       Get maximum vertical radius; necessary for defining arrays
          rad_k = MAX( rad_k_b, rad_k_t )
!
!--       When analysed grid point lies above maximum topography, set search
!--       radius to 0 if the distance between the analysed grid point and max
!--       topography height is larger than the maximum search radius
          IF ( zu(k-rad_k_b) > zu(k_max_topo) )  rad_k_b = 0
!
!--       Search within vicinity only if the vertical search radius is >0
          IF ( rad_k_b /= 0 .OR. rad_k_t /= 0 )  THEN

             !> @note shape of vicinity is larger in z direction
             !>   Shape of vicinity is two grid points larger than actual search
             !>   radius in vertical direction. The first and last grid point is
             !>   always set to 1 to asure correct detection of topography. See
             !>   function "shortest_distance" for details.
             !>   2018-03-16, gronemeier
             ALLOCATE( vicinity(-rad_k-1:rad_k+1,-rad_j:rad_j,-rad_i:rad_i) )
             ALLOCATE( vic_yz(0:rad_k+1,0:rad_j) )

             vicinity = 1

             DO  i = nxl, nxr
                DO  j = nys, nyn
!
!--                Start search only if (i/j/k) belongs to atmosphere
                   IF ( BTEST( wall_flags_total_0(k,j,i), 0 )  )  THEN
!
!--                   Reset topography within vicinity
                      vicinity(-rad_k:rad_k,:,:) = 0
!
!--                   Copy area surrounding analysed grid point into vicinity.
!--                   First, limit size of data copied to vicinity by the domain
!--                   border
                      !> @note limit copied area to 1 grid point in hor. dir.
                      !>   Ignore walls in horizontal direction which are
                      !>   further away than a single grid point. This allows to
                      !>   only search within local subdomain without the need
                      !>   of global topography information.
                      !>   The error made by this assumption are acceptable at
                      !>   the moment.
                      !>   2018-10-01, gronemeier
                      rad_i_l = MIN( 1, rad_i, i )
                      rad_i_r = MIN( 1, rad_i, nx-i )

                      rad_j_s = MIN( 1, rad_j, j )
                      rad_j_n = MIN( 1, rad_j, ny-j )

                      CALL copy_into_vicinity( k, j, i,           &
                                               -rad_k_b, rad_k_t, &
                                               -rad_j_s, rad_j_n, &
                                               -rad_i_l, rad_i_r  )
                      !> @note in case of cyclic boundaries, those parts of the
                      !>   topography which lies beyond the domain borders but
                      !>   still within the search radius still needs to be
                      !>   copied into 'vicinity'. As the effective search
                      !>   radius is limited to 1 at the moment, no further
                      !>   copying is needed. Old implementation (prior to
                      !>   2018-10-01) had this covered but used a global array.
                      !>   2018-10-01, gronemeier

!
!--                   Search for walls only if there is any within vicinity
                      IF ( MAXVAL( vicinity(-rad_k:rad_k,:,:) ) /= 0 )  THEN
!
!--                      Search within first half (positive x)
                         dist_dx = rad_i
                         DO  ii = 0, dist_dx
!
!--                         Search along vertical direction only if below
!--                         maximum topography
                            IF ( rad_k_t > 0 ) THEN
!
!--                            Search for walls within octant (+++)
                               vic_yz = vicinity(0:rad_k+1,0:rad_j,ii)
                               l_wall(k,j,i) = MIN( l_wall(k,j,i),             &
                                       shortest_distance( vic_yz, .TRUE., ii ) )
!
!--                            Search for walls within octant (+-+)
!--                            Switch order of array so that the analysed grid
!--                            point is always located at (0/0) (required by
!--                            shortest_distance").
                               vic_yz = vicinity(0:rad_k+1,0:-rad_j:-1,ii)
                               l_wall(k,j,i) = MIN( l_wall(k,j,i),             &
                                       shortest_distance( vic_yz, .TRUE., ii ) )

                            ENDIF
!
!--                         Search for walls within octant (+--)
                            vic_yz = vicinity(0:-rad_k-1:-1,0:-rad_j:-1,ii)
                            l_wall(k,j,i) = MIN( l_wall(k,j,i),                &
                                      shortest_distance( vic_yz, .FALSE., ii ) )
!
!--                         Search for walls within octant (++-)
                            vic_yz = vicinity(0:-rad_k-1:-1,0:rad_j,ii)
                            l_wall(k,j,i) = MIN( l_wall(k,j,i),                &
                                      shortest_distance( vic_yz, .FALSE., ii ) )
!
!--                         Reduce search along x by already found distance
                            dist_dx = CEILING( l_wall(k,j,i) / dx )

                         ENDDO
!
!-                       Search within second half (negative x)
                         DO  ii = 0, -dist_dx, -1
!
!--                         Search along vertical direction only if below
!--                         maximum topography
                            IF ( rad_k_t > 0 ) THEN
!
!--                            Search for walls within octant (-++)
                               vic_yz = vicinity(0:rad_k+1,0:rad_j,ii)
                               l_wall(k,j,i) = MIN( l_wall(k,j,i),             &
                                      shortest_distance( vic_yz, .TRUE., -ii ) )
!
!--                            Search for walls within octant (--+)
!--                            Switch order of array so that the analysed grid
!--                            point is always located at (0/0) (required by
!--                            shortest_distance").
                               vic_yz = vicinity(0:rad_k+1,0:-rad_j:-1,ii)
                               l_wall(k,j,i) = MIN( l_wall(k,j,i),             &
                                      shortest_distance( vic_yz, .TRUE., -ii ) )

                            ENDIF
!
!--                         Search for walls within octant (---)
                            vic_yz = vicinity(0:-rad_k-1:-1,0:-rad_j:-1,ii)
                            l_wall(k,j,i) = MIN( l_wall(k,j,i),                &
                                     shortest_distance( vic_yz, .FALSE., -ii ) )
!
!--                         Search for walls within octant (-+-)
                            vic_yz = vicinity(0:-rad_k-1:-1,0:rad_j,ii)
                            l_wall(k,j,i) = MIN( l_wall(k,j,i),                &
                                     shortest_distance( vic_yz, .FALSE., -ii ) )
!
!--                         Reduce search along x by already found distance
                            dist_dx = CEILING( l_wall(k,j,i) / dx )

                         ENDDO

                      ENDIF  !Check for any walls within vicinity

                   ELSE  !Check if (i,j,k) belongs to atmosphere

                      l_wall(k,j,i) = l_black(k)

                   ENDIF

                ENDDO  !j loop
             ENDDO  !i loop

             DEALLOCATE( vicinity )
             DEALLOCATE( vic_yz )

          ENDIF  !check vertical size of vicinity

       ENDDO  !k loop

       !$ACC ENTER DATA COPYIN(l_black(nzb:nzt+1))

    ENDIF  !LES or RANS mode

!
!-- Set lateral boundary conditions for l_wall
    CALL exchange_horiz( l_wall, nbgp )

    !$ACC ENTER DATA COPYIN(l_grid(nzb:nzt+1)) &
    !$ACC COPYIN(l_wall(nzb:nzt+1,nysg:nyng,nxlg:nxrg))

    CONTAINS
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate the shortest distance between position (i/j/k)=(0/0/0) and
!> (pos_i/jj/kk), where (jj/kk) is the position of the maximum of 'array'
!> closest to the origin (0/0) of 'array'.
!------------------------------------------------------------------------------!
    REAL(wp) FUNCTION shortest_distance( array, orientation, pos_i )

       IMPLICIT NONE

       LOGICAL, INTENT(IN) :: orientation    !< flag if array represents an array oriented upwards (true) or downwards (false)

       INTEGER(iwp), INTENT(IN) :: pos_i     !< x position of the yz-plane 'array'

       INTEGER(iwp) :: a                     !< loop index
       INTEGER(iwp) :: b                     !< loop index
       INTEGER(iwp) :: jj                    !< loop index

       INTEGER(KIND=1) :: maximum            !< maximum of array along z dimension

       INTEGER(iwp), DIMENSION(0:rad_j) :: loc_k !< location of closest wall along vertical dimension

       INTEGER(KIND=1), DIMENSION(0:rad_k+1,0:rad_j), INTENT(IN) :: array !< array containing a yz-plane at position pos_i

!
!--    Get coordinate of first maximum along vertical dimension
!--    at each y position of array (similar to function maxloc but more stable).
       DO  a = 0, rad_j
          loc_k(a) = rad_k+1
          maximum = MAXVAL( array(:,a) )
          DO  b = 0, rad_k+1
             IF ( array(b,a) == maximum )  THEN
                loc_k(a) = b
                EXIT
             ENDIF
          ENDDO
       ENDDO
!
!--    Set distance to the default maximum value (=search radius)
       shortest_distance = radius
!
!--    Calculate distance between position (0/0/0) and
!--    position (pos_i/jj/loc(jj)) and only save the shortest distance.
       IF ( orientation ) THEN  !if array is oriented upwards
          DO  jj = 0, rad_j
             shortest_distance =                                               &
                MIN( shortest_distance,                                        &
                     SQRT( MAX(REAL(pos_i, KIND=wp)*dx-0.5_wp*dx, 0.0_wp)**2   &
                         + MAX(REAL(jj, KIND=wp)*dy-0.5_wp*dy, 0.0_wp)**2      &
                         + MAX(zw(loc_k(jj)+k-1)-zu(k), 0.0_wp)**2             &
                         )                                                     &
                   )
          ENDDO
       ELSE  !if array is oriented downwards
          !> @note MAX within zw required to circumvent error at domain border
          !>   At the domain border, if non-cyclic boundary is present, the
          !>   index for zw could be -1, which will be errorneous (zw(-1) does
          !>   not exist). The MAX function limits the index to be at least 0.
          DO  jj = 0, rad_j
             shortest_distance =                                               &
                MIN( shortest_distance,                                        &
                     SQRT( MAX(REAL(pos_i, KIND=wp)*dx-0.5_wp*dx, 0.0_wp)**2   &
                         + MAX(REAL(jj, KIND=wp)*dy-0.5_wp*dy, 0.0_wp)**2      &
                         + MAX(zu(k)-zw(MAX(k-loc_k(jj),0_iwp)), 0.0_wp)**2    &
                         )                                                     &
                   )
          ENDDO
       ENDIF

    END FUNCTION

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Copy a subarray of size (kb:kt,js:jn,il:ir) centered around grid point
!> (kp,jp,ip) containing the first bit of wall_flags_total_0 into the array
!> 'vicinity'. Only copy first bit as this indicates the presence of topography.
!------------------------------------------------------------------------------!
    SUBROUTINE copy_into_vicinity( kp, jp, ip, kb, kt, js, jn, il, ir )

       IMPLICIT NONE

       INTEGER(iwp), INTENT(IN) :: il !< left loop boundary
       INTEGER(iwp), INTENT(IN) :: ip !< center position in x-direction
       INTEGER(iwp), INTENT(IN) :: ir !< right loop boundary
       INTEGER(iwp), INTENT(IN) :: jn !< northern loop boundary
       INTEGER(iwp), INTENT(IN) :: jp !< center position in y-direction
       INTEGER(iwp), INTENT(IN) :: js !< southern loop boundary
       INTEGER(iwp), INTENT(IN) :: kb !< bottom loop boundary
       INTEGER(iwp), INTENT(IN) :: kp !< center position in z-direction
       INTEGER(iwp), INTENT(IN) :: kt !< top loop boundary

       INTEGER(iwp) :: i   !< loop index
       INTEGER(iwp) :: j   !< loop index
       INTEGER(iwp) :: k   !< loop index

       DO  i = il, ir
          DO  j = js, jn
             DO  k = kb, kt
                vicinity(k,j,i) = MERGE( 0, 1,               &
                       BTEST( wall_flags_total_0(kp+k,jp+j,ip+i), 0 ) )
             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE copy_into_vicinity

 END SUBROUTINE tcm_init_mixing_length


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize virtual velocities used later in production_e.
!------------------------------------------------------------------------------!
 SUBROUTINE production_e_init

    USE arrays_3d,                                                             &
        ONLY:  drho_air_zw, zu

    USE control_parameters,                                                    &
        ONLY:  constant_flux_layer

    USE surface_layer_fluxes_mod,                                              &
        ONLY:  phi_m

    IMPLICIT NONE

    INTEGER(iwp) ::  i      !< grid index x-direction
    INTEGER(iwp) ::  j      !< grid index y-direction
    INTEGER(iwp) ::  k      !< grid index z-direction
    INTEGER(iwp) ::  m      !< running index surface elements

    REAL(wp) ::  km_sfc     !< diffusion coefficient, used to compute virtual velocities

    IF ( constant_flux_layer )  THEN
!
!--    Calculate a virtual velocity at the surface in a way that the
!--    vertical velocity gradient at k = 1 (u(k+1)-u_0) matches the
!--    Prandtl law (-w'u'/km). This gradient is used in the TKE shear
!--    production term at k=1 (see production_e_ij).
!--    The velocity gradient has to be limited in case of too small km
!--    (otherwise the timestep may be significantly reduced by large
!--    surface winds).
!--    not available in case of non-cyclic boundary conditions.
!--    Default surfaces, upward-facing
       !$OMP PARALLEL DO PRIVATE(i,j,k,m)
       !$ACC PARALLEL LOOP PRIVATE(i, j, k, m, km_sfc) &
       !$ACC PRESENT(surf_def_h(0), u, v, drho_air_zw, zu)
       DO  m = 1, surf_def_h(0)%ns

          i = surf_def_h(0)%i(m)
          j = surf_def_h(0)%j(m)
          k = surf_def_h(0)%k(m)
!
!--       Note, calculation of u_0 and v_0 is not fully accurate, as u/v
!--       and km are not on the same grid. Actually, a further
!--       interpolation of km onto the u/v-grid is necessary. However, the
!--       effect of this error is negligible.
          km_sfc = kappa * surf_def_h(0)%us(m) * surf_def_h(0)%z_mo(m) /       &
                   phi_m( surf_def_h(0)%z_mo(m) / surf_def_h(0)%ol(m) )

          surf_def_h(0)%u_0(m) = u(k+1,j,i) + surf_def_h(0)%usws(m) *          &
                                     drho_air_zw(k-1)               *          &
                                     ( zu(k+1) - zu(k-1)    )       /          &
                                     ( km_sfc  + 1.0E-20_wp )
          surf_def_h(0)%v_0(m) = v(k+1,j,i) + surf_def_h(0)%vsws(m) *          &
                                     drho_air_zw(k-1)               *          &
                                     ( zu(k+1) - zu(k-1)    )       /          &
                                     ( km_sfc  + 1.0E-20_wp )

          IF ( ABS( u(k+1,j,i) - surf_def_h(0)%u_0(m) )  >                     &
               ABS( u(k+1,j,i) - u(k-1,j,i)           )                        &
             )  surf_def_h(0)%u_0(m) = u(k-1,j,i)

          IF ( ABS( v(k+1,j,i) - surf_def_h(0)%v_0(m) )  >                     &
               ABS( v(k+1,j,i) - v(k-1,j,i)           )                        &
             )  surf_def_h(0)%v_0(m) = v(k-1,j,i)

       ENDDO
!
!--    Default surfaces, downward-facing surfaces
       !$OMP PARALLEL DO PRIVATE(i,j,k,m)
       !$ACC PARALLEL LOOP PRIVATE(i, j, k, m, km_sfc) &
       !$ACC PRESENT(surf_def_h(1), u, v, drho_air_zw, zu, km)
       DO  m = 1, surf_def_h(1)%ns

          i = surf_def_h(1)%i(m)
          j = surf_def_h(1)%j(m)
          k = surf_def_h(1)%k(m)
!
!--       Note, calculation of u_0 and v_0 is not fully accurate, as u/v
!--       and km are not on the same grid. Actually, a further
!--       interpolation of km onto the u/v-grid is necessary. However, the
!--       effect of this error is negligible.
          surf_def_h(1)%u_0(m) = u(k-1,j,i) - surf_def_h(1)%usws(m) *          &
                                     drho_air_zw(k-1) *                        &
                                     ( zu(k+1)    - zu(k-1)    )  /            &
                                     ( km(k,j,i)  + 1.0E-20_wp )
          surf_def_h(1)%v_0(m) = v(k-1,j,i) - surf_def_h(1)%vsws(m) *          &
                                     drho_air_zw(k-1) *                        &
                                     ( zu(k+1)    - zu(k-1)    )  /            &
                                     ( km(k,j,i)  + 1.0E-20_wp )

          IF ( ABS( surf_def_h(1)%u_0(m) - u(k-1,j,i) )  >                     &
               ABS( u(k+1,j,i)           - u(k-1,j,i) )                        &
             )  surf_def_h(1)%u_0(m) = u(k+1,j,i)

          IF ( ABS( surf_def_h(1)%v_0(m) - v(k-1,j,i) )  >                     &
               ABS( v(k+1,j,i)           - v(k-1,j,i) )                        &
             )  surf_def_h(1)%v_0(m) = v(k+1,j,i)

       ENDDO
!
!--    Natural surfaces, upward-facing
       !$OMP PARALLEL DO PRIVATE(i,j,k,m)
       !$ACC PARALLEL LOOP PRIVATE(i, j, k, m, km_sfc) &
       !$ACC PRESENT(surf_lsm_h, u, v, drho_air_zw, zu)
       DO  m = 1, surf_lsm_h%ns

          i = surf_lsm_h%i(m)
          j = surf_lsm_h%j(m)
          k = surf_lsm_h%k(m)
!
!--       Note, calculation of u_0 and v_0 is not fully accurate, as u/v
!--       and km are not on the same grid. Actually, a further
!--       interpolation of km onto the u/v-grid is necessary. However, the
!--       effect of this error is negligible.
          km_sfc = kappa * surf_lsm_h%us(m) * surf_lsm_h%z_mo(m) /             &
                   phi_m( surf_lsm_h%z_mo(m) / surf_lsm_h%ol(m) )

          surf_lsm_h%u_0(m) = u(k+1,j,i) + surf_lsm_h%usws(m)    *             &
                                        drho_air_zw(k-1)         *             &
                                        ( zu(k+1) - zu(k-1)    ) /             &
                                        ( km_sfc  + 1.0E-20_wp )
          surf_lsm_h%v_0(m) = v(k+1,j,i) + surf_lsm_h%vsws(m)    *             &
                                        drho_air_zw(k-1)         *             &
                                        ( zu(k+1) - zu(k-1)    ) /             &
                                        ( km_sfc  + 1.0E-20_wp )

          IF ( ABS( u(k+1,j,i) - surf_lsm_h%u_0(m) )  >                        &
               ABS( u(k+1,j,i) - u(k-1,j,i)   )                                &
             )  surf_lsm_h%u_0(m) = u(k-1,j,i)

          IF ( ABS( v(k+1,j,i) - surf_lsm_h%v_0(m) )  >                        &
               ABS( v(k+1,j,i) - v(k-1,j,i)   )                                &
             )  surf_lsm_h%v_0(m) = v(k-1,j,i)

       ENDDO
!
!--    Urban surfaces, upward-facing
       !$OMP PARALLEL DO PRIVATE(i,j,k,m)
       !$ACC PARALLEL LOOP PRIVATE(i, j, k, m, km_sfc) &
       !$ACC PRESENT(surf_usm_h, u, v, drho_air_zw, zu)
       DO  m = 1, surf_usm_h%ns

          i = surf_usm_h%i(m)
          j = surf_usm_h%j(m)
          k = surf_usm_h%k(m)
!
!--       Note, calculation of u_0 and v_0 is not fully accurate, as u/v
!--       and km are not on the same grid. Actually, a further
!--       interpolation of km onto the u/v-grid is necessary. However, the
!--       effect of this error is negligible.
          km_sfc = kappa * surf_usm_h%us(m) * surf_usm_h%z_mo(m) /             &
                   phi_m( surf_usm_h%z_mo(m) / surf_usm_h%ol(m) )

          surf_usm_h%u_0(m) = u(k+1,j,i) + surf_usm_h%usws(m)    *             &
                                        drho_air_zw(k-1)         *             &
                                        ( zu(k+1) - zu(k-1)    ) /             &
                                        ( km_sfc  + 1.0E-20_wp )
          surf_usm_h%v_0(m) = v(k+1,j,i) + surf_usm_h%vsws(m)    *             &
                                        drho_air_zw(k-1)         *             &
                                        ( zu(k+1) - zu(k-1)    ) /             &
                                        ( km_sfc  + 1.0E-20_wp )

          IF ( ABS( u(k+1,j,i) - surf_usm_h%u_0(m) )  >                        &
               ABS( u(k+1,j,i) - u(k-1,j,i)   )                                &
             )  surf_usm_h%u_0(m) = u(k-1,j,i)

          IF ( ABS( v(k+1,j,i) - surf_usm_h%v_0(m) )  >                        &
               ABS( v(k+1,j,i) - v(k-1,j,i)   )                                &
             )  surf_usm_h%v_0(m) = v(k-1,j,i)

       ENDDO

    ENDIF

 END SUBROUTINE production_e_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Execute module-specific actions for all grid points
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE tcm_actions( location )


    CHARACTER (LEN=*) ::  location !<

!    INTEGER(iwp) ::  i !<
!    INTEGER(iwp) ::  j !<
!    INTEGER(iwp) ::  k !<

!
!-- Here the module-specific actions follow
!-- No calls for single grid points are allowed at locations before and
!-- after the timestep, since these calls are not within an i,j-loop
    SELECT CASE ( location )

       CASE ( 'before_timestep' )


       CASE ( 'before_prognostic_equations' )

          IF ( .NOT. constant_diffusion )  CALL production_e_init


       CASE ( 'after_integration' )


       CASE ( 'after_timestep' )


       CASE ( 'u-tendency' )


       CASE ( 'v-tendency' )


       CASE ( 'w-tendency' )


       CASE ( 'pt-tendency' )


       CASE ( 'sa-tendency' )


       CASE ( 'e-tendency' )


       CASE ( 'q-tendency' )


       CASE ( 's-tendency' )


       CASE DEFAULT
          CONTINUE

    END SELECT

 END SUBROUTINE tcm_actions


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Execute module-specific actions for grid point i,j
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE tcm_actions_ij( i, j, location )


    CHARACTER (LEN=*) ::  location

    INTEGER(iwp) ::  i
    INTEGER(iwp) ::  j

!
!-- Here the module-specific actions follow
    SELECT CASE ( location )

       CASE ( 'u-tendency' )

!--       Next line is to avoid compiler warning about unused variables. Please remove.
          IF ( i +  j < 0 )  CONTINUE

       CASE ( 'v-tendency' )


       CASE ( 'w-tendency' )


       CASE ( 'pt-tendency' )


       CASE ( 'sa-tendency' )


       CASE ( 'e-tendency' )


       CASE ( 'q-tendency' )


       CASE ( 's-tendency' )


       CASE DEFAULT
          CONTINUE

    END SELECT

 END SUBROUTINE tcm_actions_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Prognostic equation for subgrid-scale TKE and TKE dissipation rate.
!> Vector-optimized version
!------------------------------------------------------------------------------!
 SUBROUTINE tcm_prognostic_equations

    USE control_parameters,                                                    &
        ONLY:  scalar_advec, tsc

    IMPLICIT NONE

    INTEGER(iwp) ::  i       !< loop index
    INTEGER(iwp) ::  j       !< loop index
    INTEGER(iwp) ::  k       !< loop index

    REAL(wp)     ::  sbt     !< wheighting factor for sub-time step

!
!-- If required, compute prognostic equation for turbulent kinetic
!-- energy (TKE)
    IF ( .NOT. constant_diffusion )  THEN

       CALL cpu_log( log_point_s(67), 'tke-equation', 'start' )

       sbt = tsc(2)
       IF ( .NOT. use_upstream_for_tke )  THEN
          IF ( scalar_advec == 'bc-scheme' )  THEN

             IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--             Bott-Chlond scheme always uses Euler time step. Thus:
                sbt = 1.0_wp
             ENDIF
             tend = 0.0_wp
             CALL advec_s_bc( e, 'e' )

          ENDIF
       ENDIF

!
!--    TKE-tendency terms with no communication
       IF ( scalar_advec /= 'bc-scheme'  .OR.  use_upstream_for_tke )  THEN
          IF ( use_upstream_for_tke )  THEN
             tend = 0.0_wp
             CALL advec_s_up( e )
          ELSE
             !$ACC KERNELS PRESENT(tend)
             tend = 0.0_wp
             !$ACC END KERNELS
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( ws_scheme_sca )  THEN
                   CALL advec_s_ws( advc_flags_s, e, 'e',                      &
                                    bc_dirichlet_l  .OR.  bc_radiation_l,      &
                                    bc_dirichlet_n  .OR.  bc_radiation_n,      &
                                    bc_dirichlet_r  .OR.  bc_radiation_r,      &
                                    bc_dirichlet_s  .OR.  bc_radiation_s )
                ELSE
                   CALL advec_s_pw( e )
                ENDIF
             ELSE
                CALL advec_s_up( e )
             ENDIF
          ENDIF
       ENDIF

       CALL production_e( .FALSE. )

       IF ( .NOT. humidity )  THEN
          IF ( ocean_mode )  THEN
             CALL diffusion_e( prho, prho_reference )
          ELSE
             CALL diffusion_e( pt, pt_reference )
          ENDIF
       ELSE
          CALL diffusion_e( vpt, pt_reference )
       ENDIF

!
!--    Additional sink term for flows through plant canopies
       IF ( plant_canopy )  CALL pcm_tendency( 6 )

!       CALL user_actions( 'e-tendency' ) ToDo: find general solution for circular dependency between modules

!
!--    Prognostic equation for TKE.
!--    Eliminate negative TKE values, which can occur due to numerical
!--    reasons in the course of the integration. In such cases the old TKE
!--    value is reduced by 90%.
       !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i, j, k) &
       !$ACC PRESENT(e, tend, te_m, wall_flags_total_0) &
       !$ACC PRESENT(tsc(3:3)) &
       !$ACC PRESENT(e_p)
       DO  i = nxl, nxr
          DO  j = nys, nyn
             !following directive is required to vectorize on Intel19
             !DIR$ IVDEP
             DO  k = nzb+1, nzt
                e_p(k,j,i) = e(k,j,i) + ( dt_3d * ( sbt * tend(k,j,i) +        &
                                                 tsc(3) * te_m(k,j,i) )        &
                                        )                                      &
                                   * MERGE( 1.0_wp, 0.0_wp,                    &
                                       BTEST( wall_flags_total_0(k,j,i), 0 )   &
                                          )
                IF ( e_p(k,j,i) < 0.0_wp )  e_p(k,j,i) = 0.1_wp * e(k,j,i)
             ENDDO
          ENDDO
       ENDDO

!
!--    Calculate tendencies for the next Runge-Kutta step
       IF ( timestep_scheme(1:5) == 'runge' )  THEN
          IF ( intermediate_timestep_count == 1 )  THEN
             !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i, j, k) &
             !$ACC PRESENT(tend, te_m)
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt
                      te_m(k,j,i) = tend(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSEIF ( intermediate_timestep_count < &
                   intermediate_timestep_count_max )  THEN
             !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i, j, k) &
             !$ACC PRESENT(tend, te_m)
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt
                      te_m(k,j,i) =   -9.5625_wp * tend(k,j,i)                 &
                                     + 5.3125_wp * te_m(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDIF

       CALL cpu_log( log_point_s(67), 'tke-equation', 'stop' )

    ENDIF   ! TKE equation

!
!-- If required, compute prognostic equation for TKE dissipation rate
    IF ( rans_tke_e )  THEN

       CALL cpu_log( log_point_s(64), 'diss-equation', 'start' )

       sbt = tsc(2)
       IF ( .NOT. use_upstream_for_tke )  THEN
          IF ( scalar_advec == 'bc-scheme' )  THEN

             IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--             Bott-Chlond scheme always uses Euler time step. Thus:
                sbt = 1.0_wp
             ENDIF
             tend = 0.0_wp
             CALL advec_s_bc( diss, 'diss' )

          ENDIF
       ENDIF

!
!--    dissipation-tendency terms with no communication
       IF ( scalar_advec /= 'bc-scheme'  .OR.  use_upstream_for_tke )  THEN
          IF ( use_upstream_for_tke )  THEN
             tend = 0.0_wp
             CALL advec_s_up( diss )
          ELSE
             tend = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( ws_scheme_sca )  THEN
                   CALL advec_s_ws( advc_flags_s, diss, 'diss',                &
                                    bc_dirichlet_l  .OR.  bc_radiation_l,      &
                                    bc_dirichlet_n  .OR.  bc_radiation_n,      &
                                    bc_dirichlet_r  .OR.  bc_radiation_r,      &
                                    bc_dirichlet_s  .OR.  bc_radiation_s )
                ELSE
                   CALL advec_s_pw( diss )
                ENDIF
             ELSE
                CALL advec_s_up( diss )
             ENDIF
          ENDIF
       ENDIF
!
!--    Production of TKE dissipation rate
       CALL production_e( .TRUE. )
!
!--    Diffusion term of TKE dissipation rate
       CALL diffusion_diss
!
!--    Additional sink term for flows through plant canopies
!        IF ( plant_canopy )  CALL pcm_tendency( ? )         !> @todo not yet implemented

!       CALL user_actions( 'e-tendency' ) ToDo: find general solution for circular dependency between modules

!
!--    Prognostic equation for TKE dissipation.
!--    Eliminate negative dissipation values, which can occur due to numerical
!--    reasons in the course of the integration. In such cases the old
!--    dissipation value is reduced by 90%.
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                diss_p(k,j,i) = diss(k,j,i) + ( dt_3d * ( sbt * tend(k,j,i) +  &
                                                 tsc(3) * tdiss_m(k,j,i) )     &
                                        )                                      &
                                   * MERGE( 1.0_wp, 0.0_wp,                    &
                                       BTEST( wall_flags_total_0(k,j,i), 0 )   &
                                          )
                IF ( diss_p(k,j,i) < 0.0_wp )                                  &
                   diss_p(k,j,i) = 0.1_wp * diss(k,j,i)
             ENDDO
          ENDDO
       ENDDO

!
!--    Calculate tendencies for the next Runge-Kutta step
       IF ( timestep_scheme(1:5) == 'runge' )  THEN
          IF ( intermediate_timestep_count == 1 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt
                      tdiss_m(k,j,i) = tend(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSEIF ( intermediate_timestep_count < &
                   intermediate_timestep_count_max )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt
                      tdiss_m(k,j,i) =   -9.5625_wp * tend(k,j,i)              &
                                        + 5.3125_wp * tdiss_m(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDIF

       CALL cpu_log( log_point_s(64), 'diss-equation', 'stop' )

    ENDIF

 END SUBROUTINE tcm_prognostic_equations


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Prognostic equation for subgrid-scale TKE and TKE dissipation rate.
!> Cache-optimized version
!------------------------------------------------------------------------------!
 SUBROUTINE tcm_prognostic_equations_ij( i, j, i_omp, tn )

    USE arrays_3d,                                                             &
        ONLY:  diss_l_diss, diss_l_e, diss_s_diss, diss_s_e, flux_l_diss,      &
               flux_l_e, flux_s_diss, flux_s_e

    USE control_parameters,                                                    &
        ONLY:  tsc

    IMPLICIT NONE

    INTEGER(iwp) ::  i       !< loop index x direction
    INTEGER(iwp) ::  i_omp   !< first loop index of i-loop in prognostic_equations
    INTEGER(iwp) ::  j       !< loop index y direction
    INTEGER(iwp) ::  k       !< loop index z direction
    INTEGER(iwp) ::  tn      !< task number of openmp task

!
!-- If required, compute prognostic equation for turbulent kinetic
!-- energy (TKE)
    IF ( .NOT. constant_diffusion )  THEN

!
!--    Tendency-terms for TKE
       tend(:,j,i) = 0.0_wp
       IF ( timestep_scheme(1:5) == 'runge'  &
           .AND.  .NOT. use_upstream_for_tke )  THEN
           IF ( ws_scheme_sca )  THEN
               CALL advec_s_ws( advc_flags_s,                                  &
                                i, j, e, 'e', flux_s_e, diss_s_e,              &
                                flux_l_e, diss_l_e , i_omp, tn,                &
                                bc_dirichlet_l  .OR.  bc_radiation_l,          &
                                bc_dirichlet_n  .OR.  bc_radiation_n,          &
                                bc_dirichlet_r  .OR.  bc_radiation_r,          &
                                bc_dirichlet_s  .OR.  bc_radiation_s )
           ELSE
               CALL advec_s_pw( i, j, e )
           ENDIF
       ELSE
          CALL advec_s_up( i, j, e )
       ENDIF

       CALL production_e_ij( i, j, .FALSE. )

       IF ( .NOT. humidity )  THEN
          IF ( ocean_mode )  THEN
             CALL diffusion_e_ij( i, j, prho, prho_reference )
          ELSE
             CALL diffusion_e_ij( i, j, pt, pt_reference )
          ENDIF
       ELSE
          CALL diffusion_e_ij( i, j, vpt, pt_reference )
       ENDIF

!
!--    Additional sink term for flows through plant canopies
       IF ( plant_canopy )  CALL pcm_tendency( i, j, 6 )

!       CALL user_actions( i, j, 'e-tendency' ) ToDo: find general solution for circular dependency between modules

!
!--    Prognostic equation for TKE.
!--    Eliminate negative TKE values, which can occur due to numerical
!--    reasons in the course of the integration. In such cases the old
!--    TKE value is reduced by 90%.
       DO  k = nzb+1, nzt
          e_p(k,j,i) = e(k,j,i) + ( dt_3d * ( tsc(2) * tend(k,j,i) +           &
                                              tsc(3) * te_m(k,j,i) )           &
                                  )                                            &
                                 * MERGE( 1.0_wp, 0.0_wp,                      &
                                    BTEST( wall_flags_total_0(k,j,i), 0 )      &
                                        )
          IF ( e_p(k,j,i) <= 0.0_wp )  e_p(k,j,i) = 0.1_wp * e(k,j,i)
       ENDDO

!
!--    Calculate tendencies for the next Runge-Kutta step
       IF ( timestep_scheme(1:5) == 'runge' )  THEN
          IF ( intermediate_timestep_count == 1 )  THEN
             DO  k = nzb+1, nzt
                te_m(k,j,i) = tend(k,j,i)
             ENDDO
          ELSEIF ( intermediate_timestep_count < &
                   intermediate_timestep_count_max )  THEN
             DO  k = nzb+1, nzt
                te_m(k,j,i) =   -9.5625_wp * tend(k,j,i) +                     &
                                 5.3125_wp * te_m(k,j,i)
             ENDDO
          ENDIF
       ENDIF

    ENDIF   ! TKE equation

!
!-- If required, compute prognostic equation for TKE dissipation rate
    IF ( rans_tke_e )  THEN
!
!--    Tendency-terms for dissipation
       tend(:,j,i) = 0.0_wp
       IF ( timestep_scheme(1:5) == 'runge'  &
           .AND.  .NOT. use_upstream_for_tke )  THEN
           IF ( ws_scheme_sca )  THEN
               CALL advec_s_ws( advc_flags_s,                                  &
                                i, j, diss, 'diss', flux_s_diss, diss_s_diss,  &
                                flux_l_diss, diss_l_diss, i_omp, tn,           &
                                bc_dirichlet_l  .OR.  bc_radiation_l,          &
                                bc_dirichlet_n  .OR.  bc_radiation_n,          &
                                bc_dirichlet_r  .OR.  bc_radiation_r,          &
                                bc_dirichlet_s  .OR.  bc_radiation_s )
           ELSE
               CALL advec_s_pw( i, j, diss )
           ENDIF
       ELSE
          CALL advec_s_up( i, j, diss )
       ENDIF
!
!--    Production of TKE dissipation rate
       CALL production_e_ij( i, j, .TRUE. )
!
!--    Diffusion term of TKE dissipation rate
       CALL diffusion_diss_ij( i, j )
!
!--    Additional sink term for flows through plant canopies
!        IF ( plant_canopy )  CALL pcm_tendency( i, j, ? )     !> @todo not yet implemented

!       CALL user_actions( i, j, 'diss-tendency' ) ToDo: find general solution for circular dependency between modules

!
!--    Prognostic equation for TKE dissipation
!--    Eliminate negative dissipation values, which can occur due to
!--    numerical reasons in the course of the integration. In such cases
!--    the old dissipation value is reduced by 90%.
       DO  k = nzb+1, nzt
          diss_p(k,j,i) = diss(k,j,i) + ( dt_3d * ( tsc(2) * tend(k,j,i) +     &
                                                    tsc(3) * tdiss_m(k,j,i) )  &
                                        )                                      &
                                        * MERGE( 1.0_wp, 0.0_wp,               &
                                          BTEST( wall_flags_total_0(k,j,i), 0 )&
                                               )
       ENDDO

!
!--    Calculate tendencies for the next Runge-Kutta step
       IF ( timestep_scheme(1:5) == 'runge' )  THEN
          IF ( intermediate_timestep_count == 1 )  THEN
             DO  k = nzb+1, nzt
                tdiss_m(k,j,i) = tend(k,j,i)
             ENDDO
          ELSEIF ( intermediate_timestep_count < &
                   intermediate_timestep_count_max )  THEN
             DO  k = nzb+1, nzt
                tdiss_m(k,j,i) =   -9.5625_wp * tend(k,j,i) +                  &
                                    5.3125_wp * tdiss_m(k,j,i)
             ENDDO
          ENDIF
       ENDIF

    ENDIF   ! dissipation equation

 END SUBROUTINE tcm_prognostic_equations_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Production terms (shear + buoyancy) of the TKE.
!> Vector-optimized version
!> @warning The case with constant_flux_layer = F and use_surface_fluxes = T is
!>          not considered well!
!------------------------------------------------------------------------------!
 SUBROUTINE production_e( diss_production )

    USE arrays_3d,                                                             &
        ONLY:  ddzw, dd2zu, drho_air_zw, q, ql, d_exner, exner

    USE control_parameters,                                                    &
        ONLY:  cloud_droplets, constant_flux_layer, neutral,                   &
               rho_reference, use_single_reference_value, use_surface_fluxes,  &
               use_top_fluxes

    USE grid_variables,                                                        &
        ONLY:  ddx, dx, ddy, dy

    USE bulk_cloud_model_mod,                                                  &
        ONLY:  bulk_cloud_model

    IMPLICIT NONE

    LOGICAL :: diss_production

    INTEGER(iwp) ::  i       !< running index x-direction
    INTEGER(iwp) ::  j       !< running index y-direction
    INTEGER(iwp) ::  k       !< running index z-direction
    INTEGER(iwp) ::  l       !< running index for different surface type orientation
    INTEGER(iwp) ::  m       !< running index surface elements
    INTEGER(iwp) ::  surf_e  !< end index of surface elements at given i-j position
    INTEGER(iwp) ::  surf_s  !< start index of surface elements at given i-j position
    INTEGER(iwp) ::  flag_nr !< number of masking flag

    REAL(wp)     ::  def         !< ( du_i/dx_j + du_j/dx_i ) * du_i/dx_j
    REAL(wp)     ::  flag        !< flag to mask topography
    REAL(wp)     ::  k1          !< temporary factor
    REAL(wp)     ::  k2          !< temporary factor
    REAL(wp)     ::  km_neutral  !< diffusion coefficient assuming neutral conditions - used to compute shear production at surfaces
    REAL(wp)     ::  theta       !< virtual potential temperature
    REAL(wp)     ::  temp        !< theta * Exner-function
    REAL(wp)     ::  sign_dir    !< sign of wall-tke flux, depending on wall orientation
    REAL(wp)     ::  usvs        !< momentum flux u"v"
    REAL(wp)     ::  vsus        !< momentum flux v"u"
    REAL(wp)     ::  wsus        !< momentum flux w"u"
    REAL(wp)     ::  wsvs        !< momentum flux w"v"

    REAL(wp), DIMENSION(nzb+1:nzt) ::  dudx  !< Gradient of u-component in x-direction
    REAL(wp), DIMENSION(nzb+1:nzt) ::  dudy  !< Gradient of u-component in y-direction
    REAL(wp), DIMENSION(nzb+1:nzt) ::  dudz  !< Gradient of u-component in z-direction
    REAL(wp), DIMENSION(nzb+1:nzt) ::  dvdx  !< Gradient of v-component in x-direction
    REAL(wp), DIMENSION(nzb+1:nzt) ::  dvdy  !< Gradient of v-component in y-direction
    REAL(wp), DIMENSION(nzb+1:nzt) ::  dvdz  !< Gradient of v-component in z-direction
    REAL(wp), DIMENSION(nzb+1:nzt) ::  dwdx  !< Gradient of w-component in x-direction
    REAL(wp), DIMENSION(nzb+1:nzt) ::  dwdy  !< Gradient of w-component in y-direction
    REAL(wp), DIMENSION(nzb+1:nzt) ::  dwdz  !< Gradient of w-component in z-direction
    REAL(wp), DIMENSION(nzb+1:nzt) ::  tmp_flux  !< temporary flux-array in z-direction



!
!-- Calculate TKE production by shear. Calculate gradients at all grid
!-- points first, gradients at surface-bounded grid points will be
!-- overwritten further below.
    !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(i, j, l) &
    !$ACC PRIVATE(surf_s, surf_e) &
    !$ACC PRIVATE(dudx(:), dudy(:), dudz(:), dvdx(:), dvdy(:), dvdz(:), dwdx(:), dwdy(:), dwdz(:)) &
    !$ACC PRESENT(e, u, v, w, diss, dd2zu, ddzw, km, wall_flags_total_0) &
    !$ACC PRESENT(tend) &
    !$ACC PRESENT(surf_def_h(0:1), surf_def_v(0:3)) &
    !$ACC PRESENT(surf_lsm_h, surf_lsm_v(0:3)) &
    !$ACC PRESENT(surf_usm_h, surf_usm_v(0:3))
    DO  i = nxl, nxr
       DO  j = nys, nyn
          !$ACC LOOP PRIVATE(k)
          DO  k = nzb+1, nzt

             dudx(k) =           ( u(k,j,i+1) - u(k,j,i)     ) * ddx
             dudy(k) = 0.25_wp * ( u(k,j+1,i) + u(k,j+1,i+1) -                 &
                                   u(k,j-1,i) - u(k,j-1,i+1) ) * ddy
             dudz(k) = 0.5_wp  * ( u(k+1,j,i) + u(k+1,j,i+1) -                 &
                                   u(k-1,j,i) - u(k-1,j,i+1) ) * dd2zu(k)

             dvdx(k) = 0.25_wp * ( v(k,j,i+1) + v(k,j+1,i+1) -                 &
                                   v(k,j,i-1) - v(k,j+1,i-1) ) * ddx
             dvdy(k) =           ( v(k,j+1,i) - v(k,j,i)     ) * ddy
             dvdz(k) = 0.5_wp  * ( v(k+1,j,i) + v(k+1,j+1,i) -                 &
                                     v(k-1,j,i) - v(k-1,j+1,i) ) * dd2zu(k)

             dwdx(k) = 0.25_wp * ( w(k,j,i+1) + w(k-1,j,i+1) -                 &
                                   w(k,j,i-1) - w(k-1,j,i-1) ) * ddx
             dwdy(k) = 0.25_wp * ( w(k,j+1,i) + w(k-1,j+1,i) -                 &
                                   w(k,j-1,i) - w(k-1,j-1,i) ) * ddy
             dwdz(k) =           ( w(k,j,i)   - w(k-1,j,i)   ) * ddzw(k)

          ENDDO


          flag_nr = 29


          IF ( constant_flux_layer )  THEN
!

             flag_nr = 0

!--          Position beneath wall
!--          (2) - Will allways be executed.
!--          'bottom and wall: use u_0,v_0 and wall functions'
!
!--          Compute gradients at north- and south-facing surfaces.
!--          First, for default surfaces, then for urban surfaces.
!--          Note, so far no natural vertical surfaces implemented
             DO  l = 0, 1
                surf_s = surf_def_v(l)%start_index(j,i)
                surf_e = surf_def_v(l)%end_index(j,i)
                !$ACC LOOP PRIVATE(m, k, usvs, wsvs, km_neutral, sign_dir)
                DO  m = surf_s, surf_e
                   k           = surf_def_v(l)%k(m)
                   usvs        = surf_def_v(l)%mom_flux_tke(0,m)
                   wsvs        = surf_def_v(l)%mom_flux_tke(1,m)

                   km_neutral = kappa * ( usvs**2 + wsvs**2 )**0.25_wp         &
                                   * 0.5_wp * dy
!
!--                -1.0 for right-facing wall, 1.0 for left-facing wall
                   sign_dir = MERGE( 1.0_wp, -1.0_wp,                          &
                              BTEST( wall_flags_total_0(k,j-1,i), flag_nr ) )
                   dudy(k) = sign_dir * usvs / ( km_neutral + 1E-10_wp )
                   dwdy(k) = sign_dir * wsvs / ( km_neutral + 1E-10_wp )
                ENDDO
!
!--             Natural surfaces
                surf_s = surf_lsm_v(l)%start_index(j,i)
                surf_e = surf_lsm_v(l)%end_index(j,i)
                !$ACC LOOP PRIVATE(m, k, usvs, wsvs, km_neutral, sign_dir)
                DO  m = surf_s, surf_e
                   k           = surf_lsm_v(l)%k(m)
                   usvs        = surf_lsm_v(l)%mom_flux_tke(0,m)
                   wsvs        = surf_lsm_v(l)%mom_flux_tke(1,m)

                   km_neutral = kappa * ( usvs**2 + wsvs**2 )**0.25_wp         &
                                   * 0.5_wp * dy
!
!--                -1.0 for right-facing wall, 1.0 for left-facing wall
                   sign_dir = MERGE( 1.0_wp, -1.0_wp,                          &
                              BTEST( wall_flags_total_0(k,j-1,i), flag_nr ) )
                   dudy(k) = sign_dir * usvs / ( km_neutral + 1E-10_wp )
                   dwdy(k) = sign_dir * wsvs / ( km_neutral + 1E-10_wp )
                ENDDO
!
!--             Urban surfaces
                surf_s = surf_usm_v(l)%start_index(j,i)
                surf_e = surf_usm_v(l)%end_index(j,i)
                !$ACC LOOP PRIVATE(m, k, usvs, wsvs, km_neutral, sign_dir)
                DO  m = surf_s, surf_e
                   k           = surf_usm_v(l)%k(m)
                   usvs        = surf_usm_v(l)%mom_flux_tke(0,m)
                   wsvs        = surf_usm_v(l)%mom_flux_tke(1,m)

                   km_neutral = kappa * ( usvs**2 + wsvs**2 )**0.25_wp         &
                                   * 0.5_wp * dy
!
!--                -1.0 for right-facing wall, 1.0 for left-facing wall
                   sign_dir = MERGE( 1.0_wp, -1.0_wp,                          &
                              BTEST( wall_flags_total_0(k,j-1,i), flag_nr ) )
                   dudy(k) = sign_dir * usvs / ( km_neutral + 1E-10_wp )
                   dwdy(k) = sign_dir * wsvs / ( km_neutral + 1E-10_wp )
                ENDDO
             ENDDO
!
!--          Compute gradients at east- and west-facing walls
             DO  l = 2, 3
                surf_s = surf_def_v(l)%start_index(j,i)
                surf_e = surf_def_v(l)%end_index(j,i)
                !$ACC LOOP PRIVATE(m, k, vsus, wsus, km_neutral, sign_dir)
                DO  m = surf_s, surf_e
                   k     = surf_def_v(l)%k(m)
                   vsus  = surf_def_v(l)%mom_flux_tke(0,m)
                   wsus  = surf_def_v(l)%mom_flux_tke(1,m)

                   km_neutral = kappa * ( vsus**2 + wsus**2 )**0.25_wp         &
                                      * 0.5_wp * dx
!
!--                -1.0 for right-facing wall, 1.0 for left-facing wall
                   sign_dir = MERGE( 1.0_wp, -1.0_wp,                          &
                              BTEST( wall_flags_total_0(k,j,i-1), flag_nr ) )
                   dvdx(k) = sign_dir * vsus / ( km_neutral + 1E-10_wp )
                   dwdx(k) = sign_dir * wsus / ( km_neutral + 1E-10_wp )
                ENDDO
!
!--             Natural surfaces
                surf_s = surf_lsm_v(l)%start_index(j,i)
                surf_e = surf_lsm_v(l)%end_index(j,i)
                !$ACC LOOP PRIVATE(m, k, vsus, wsus, km_neutral, sign_dir)
                DO  m = surf_s, surf_e
                   k     = surf_lsm_v(l)%k(m)
                   vsus  = surf_lsm_v(l)%mom_flux_tke(0,m)
                   wsus  = surf_lsm_v(l)%mom_flux_tke(1,m)

                   km_neutral = kappa * ( vsus**2 + wsus**2 )**0.25_wp         &
                                      * 0.5_wp * dx
!
!--                -1.0 for right-facing wall, 1.0 for left-facing wall
                   sign_dir = MERGE( 1.0_wp, -1.0_wp,                          &
                              BTEST( wall_flags_total_0(k,j,i-1), flag_nr ) )
                   dvdx(k) = sign_dir * vsus / ( km_neutral + 1E-10_wp )
                   dwdx(k) = sign_dir * wsus / ( km_neutral + 1E-10_wp )
                ENDDO
!
!--             Urban surfaces
                surf_s = surf_usm_v(l)%start_index(j,i)
                surf_e = surf_usm_v(l)%end_index(j,i)
                !$ACC LOOP PRIVATE(m, k, vsus, wsus, km_neutral, sign_dir)
                DO  m = surf_s, surf_e
                   k     = surf_usm_v(l)%k(m)
                   vsus  = surf_usm_v(l)%mom_flux_tke(0,m)
                   wsus  = surf_usm_v(l)%mom_flux_tke(1,m)

                   km_neutral = kappa * ( vsus**2 + wsus**2 )**0.25_wp         &
                                      * 0.5_wp * dx
!
!--                -1.0 for right-facing wall, 1.0 for left-facing wall
                   sign_dir = MERGE( 1.0_wp, -1.0_wp,                          &
                              BTEST( wall_flags_total_0(k,j,i-1), flag_nr ) )
                   dvdx(k) = sign_dir * vsus / ( km_neutral + 1E-10_wp )
                   dwdx(k) = sign_dir * wsus / ( km_neutral + 1E-10_wp )
                ENDDO
             ENDDO
!
!--          Compute gradients at upward-facing surfaces
             surf_s = surf_def_h(0)%start_index(j,i)
             surf_e = surf_def_h(0)%end_index(j,i)
             !$ACC LOOP PRIVATE(m, k)
             DO  m = surf_s, surf_e
                k = surf_def_h(0)%k(m)
!
!--             Please note, actually, an interpolation of u_0 and v_0
!--             onto the grid center would be required. However, this
!--             would require several data transfers between 2D-grid and
!--             wall type. The effect of this missing interpolation is
!--             negligible. (See also production_e_init).
                dudz(k) = ( u(k+1,j,i) - surf_def_h(0)%u_0(m) ) * dd2zu(k)
                dvdz(k) = ( v(k+1,j,i) - surf_def_h(0)%v_0(m) ) * dd2zu(k)

             ENDDO
!
!--          Natural surfaces
             surf_s = surf_lsm_h%start_index(j,i)
             surf_e = surf_lsm_h%end_index(j,i)
             !$ACC LOOP PRIVATE(m, k)
             DO  m = surf_s, surf_e
                k = surf_lsm_h%k(m)

                dudz(k) = ( u(k+1,j,i) - surf_lsm_h%u_0(m) ) * dd2zu(k)
                dvdz(k) = ( v(k+1,j,i) - surf_lsm_h%v_0(m) ) * dd2zu(k)

             ENDDO
!
!--          Urban surfaces
             surf_s = surf_usm_h%start_index(j,i)
             surf_e = surf_usm_h%end_index(j,i)
             !$ACC LOOP PRIVATE(m, k)
             DO  m = surf_s, surf_e
                k = surf_usm_h%k(m)

                dudz(k) = ( u(k+1,j,i) - surf_usm_h%u_0(m) ) * dd2zu(k)
                dvdz(k) = ( v(k+1,j,i) - surf_usm_h%v_0(m) ) * dd2zu(k)

             ENDDO
!
!--          Compute gradients at downward-facing walls, only for
!--          non-natural default surfaces
             surf_s = surf_def_h(1)%start_index(j,i)
             surf_e = surf_def_h(1)%end_index(j,i)
             !$ACC LOOP PRIVATE(m, k)
             DO  m = surf_s, surf_e
                k = surf_def_h(1)%k(m)

                dudz(k) = ( surf_def_h(1)%u_0(m) - u(k-1,j,i) ) * dd2zu(k)
                dvdz(k) = ( surf_def_h(1)%v_0(m) - v(k-1,j,i) ) * dd2zu(k)

             ENDDO


          ENDIF


          !$ACC LOOP PRIVATE(k, def, flag)
          DO  k = nzb+1, nzt

             def = 2.0_wp * ( dudx(k)**2 + dvdy(k)**2 + dwdz(k)**2 ) +         &
                              dudy(k)**2 + dvdx(k)**2 + dwdx(k)**2 +           &
                              dwdy(k)**2 + dudz(k)**2 + dvdz(k)**2 +           &
                   2.0_wp * ( dvdx(k)*dudy(k) + dwdx(k)*dudz(k) +              &
                              dwdy(k)*dvdz(k) )

             IF ( def < 0.0_wp )  def = 0.0_wp

             flag = MERGE( 1.0_wp, 0.0_wp, BTEST(wall_flags_total_0(k,j,i),flag_nr) )

             IF ( .NOT. diss_production )  THEN

!--             Compute tendency for TKE-production from shear
                tend(k,j,i) = tend(k,j,i) + km(k,j,i) * def * flag

             ELSE

!--             RANS mode: Compute tendency for dissipation-rate-production from shear
                tend(k,j,i) = tend(k,j,i) + km(k,j,i) * def * flag *           &
                              diss(k,j,i)/( e(k,j,i) + 1.0E-20_wp ) * c_1

             ENDIF

          ENDDO


       ENDDO
    ENDDO

!
!-- If required, calculate TKE production by buoyancy
    IF ( .NOT. neutral )  THEN

       IF ( .NOT. humidity )  THEN

          IF ( ocean_mode )  THEN
!
!--          So far in the ocean no special treatment of density flux
!--          in the bottom and top surface layer
             DO  i = nxl, nxr
                DO  j = nys, nyn

                   DO  k = nzb+1, nzt
                      tmp_flux(k) = kh(k,j,i) * ( prho(k+1,j,i) - prho(k-1,j,i) ) * dd2zu(k)
                   ENDDO
!
!--                Treatment of near-surface grid points, at up- and down-
!--                ward facing surfaces
                   IF ( use_surface_fluxes )  THEN
                      DO  l = 0, 1
                         surf_s = surf_def_h(l)%start_index(j,i)
                         surf_e = surf_def_h(l)%end_index(j,i)
                         DO  m = surf_s, surf_e
                            k = surf_def_h(l)%k(m)
                            tmp_flux(k) = drho_air_zw(k-1) * surf_def_h(l)%shf(m)
                         ENDDO
                      ENDDO
                   ENDIF

                   IF ( use_top_fluxes )  THEN
                      surf_s = surf_def_h(2)%start_index(j,i)
                      surf_e = surf_def_h(2)%end_index(j,i)
                      DO  m = surf_s, surf_e
                         k = surf_def_h(2)%k(m)
                         tmp_flux(k) = drho_air_zw(k) * surf_def_h(2)%shf(m)
                      ENDDO
                   ENDIF

                   IF ( .NOT. diss_production )  THEN

!--                   Compute tendency for TKE-production from shear
                      DO  k = nzb+1, nzt
                         flag = MERGE( 1.0_wp, 0.0_wp, BTEST(wall_flags_total_0(k,j,i),0) )
                         tend(k,j,i) = tend(k,j,i) + flag * tmp_flux(k) * ( g / &
                                       MERGE( rho_reference, prho(k,j,i),       &
                                              use_single_reference_value ) )
                      ENDDO

                   ELSE

!--                   RANS mode: Compute tendency for dissipation-rate-production from shear
                      DO  k = nzb+1, nzt
                         flag = MERGE( 1.0_wp, 0.0_wp, BTEST(wall_flags_total_0(k,j,i),0) )
                         tend(k,j,i) = tend(k,j,i) + flag * tmp_flux(k) * ( g / &
                                       MERGE( rho_reference, prho(k,j,i),       &
                                              use_single_reference_value ) ) *  &
                                       diss(k,j,i)/( e(k,j,i) + 1.0E-20_wp ) *  &
                                       c_3
                      ENDDO

                   ENDIF

                ENDDO
             ENDDO

          ELSE ! or IF ( .NOT. ocean_mode )  THEN

             !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(i, j) &
             !$ACC PRIVATE(surf_s, surf_e) &
             !$ACC PRIVATE(tmp_flux(nzb+1:nzt)) &
             !$ACC PRESENT(e, diss, kh, pt, dd2zu, drho_air_zw, wall_flags_total_0) &
             !$ACC PRESENT(tend) &
             !$ACC PRESENT(surf_def_h(0:2)) &
             !$ACC PRESENT(surf_lsm_h) &
             !$ACC PRESENT(surf_usm_h)
             DO  i = nxl, nxr
                DO  j = nys, nyn

                   !$ACC LOOP PRIVATE(k)
                   DO  k = nzb+1, nzt
                      tmp_flux(k) = -1.0_wp * kh(k,j,i) * ( pt(k+1,j,i) - pt(k-1,j,i) ) * dd2zu(k)
                   ENDDO

                   IF ( use_surface_fluxes )  THEN
!
!--                   Default surfaces, up- and downward-facing
                      DO  l = 0, 1
                         surf_s = surf_def_h(l)%start_index(j,i)
                         surf_e = surf_def_h(l)%end_index(j,i)
                         !$ACC LOOP PRIVATE(m, k)
                         DO  m = surf_s, surf_e
                            k = surf_def_h(l)%k(m)
                            tmp_flux(k) = drho_air_zw(k-1) * surf_def_h(l)%shf(m)
                         ENDDO
                      ENDDO
!
!--                   Natural surfaces
                      surf_s = surf_lsm_h%start_index(j,i)
                      surf_e = surf_lsm_h%end_index(j,i)
                      !$ACC LOOP PRIVATE(m, k)
                      DO  m = surf_s, surf_e
                         k = surf_lsm_h%k(m)
                         tmp_flux(k) = drho_air_zw(k-1) * surf_lsm_h%shf(m)
                      ENDDO
!
!--                   Urban surfaces
                      surf_s = surf_usm_h%start_index(j,i)
                      surf_e = surf_usm_h%end_index(j,i)
                      !$ACC LOOP PRIVATE(m, k)
                      DO  m = surf_s, surf_e
                         k = surf_usm_h%k(m)
                         tmp_flux(k) = drho_air_zw(k-1) * surf_usm_h%shf(m)
                      ENDDO
                   ENDIF

                   IF ( use_top_fluxes )  THEN
                      surf_s = surf_def_h(2)%start_index(j,i)
                      surf_e = surf_def_h(2)%end_index(j,i)
                      !$ACC LOOP PRIVATE(m, k)
                      DO  m = surf_s, surf_e
                         k = surf_def_h(2)%k(m)
                         tmp_flux(k) = drho_air_zw(k) * surf_def_h(2)%shf(m)
                      ENDDO
                   ENDIF

                   IF ( .NOT. diss_production )  THEN

!--                   Compute tendency for TKE-production from shear
                     !$ACC LOOP PRIVATE(k, flag)
                      DO  k = nzb+1, nzt
                         flag = MERGE( 1.0_wp, 0.0_wp, BTEST(wall_flags_total_0(k,j,i),0) )
                         tend(k,j,i) = tend(k,j,i) + flag * tmp_flux(k) * ( g / &
                                       MERGE( pt_reference, pt(k,j,i),          &
                                              use_single_reference_value ) )
                      ENDDO

                   ELSE

!--                   RANS mode: Compute tendency for dissipation-rate-production from shear
                      DO  k = nzb+1, nzt
                         flag = MERGE( 1.0_wp, 0.0_wp, BTEST(wall_flags_total_0(k,j,i),0) )
                         tend(k,j,i) = tend(k,j,i) + flag * tmp_flux(k) * ( g / &
                                       MERGE( pt_reference, pt(k,j,i),          &
                                              use_single_reference_value ) ) *  &
                                       diss(k,j,i)/( e(k,j,i) + 1.0E-20_wp ) *  &
                                       c_3
                      ENDDO

                   ENDIF

                ENDDO
             ENDDO

          ENDIF ! from IF ( .NOT. ocean_mode )

       ELSE ! or IF ( humidity )  THEN

          DO  i = nxl, nxr
             DO  j = nys, nyn

                DO  k = nzb+1, nzt

                   IF ( .NOT. bulk_cloud_model .AND. .NOT. cloud_droplets ) THEN
                      k1 = 1.0_wp + 0.61_wp * q(k,j,i)
                      k2 = 0.61_wp * pt(k,j,i)
                      tmp_flux(k) = -1.0_wp * kh(k,j,i) *                      &
                                      ( k1 * ( pt(k+1,j,i) - pt(k-1,j,i) ) +   &
                                        k2 * ( q(k+1,j,i)  - q(k-1,j,i) )      &
                                      ) * dd2zu(k)
                   ELSE IF ( bulk_cloud_model )  THEN
                      IF ( ql(k,j,i) == 0.0_wp )  THEN
                         k1 = 1.0_wp + 0.61_wp * q(k,j,i)
                         k2 = 0.61_wp * pt(k,j,i)
                      ELSE
                         theta = pt(k,j,i) + d_exner(k) * lv_d_cp * ql(k,j,i)
                         temp  = theta * exner(k)
                         k1 = ( 1.0_wp - q(k,j,i) + 1.61_wp *                  &
                                       ( q(k,j,i) - ql(k,j,i) ) *              &
                              ( 1.0_wp + rd_d_rv * lv_d_rd / temp ) ) /        &
                              ( 1.0_wp + rd_d_rv * lv_d_rd * lv_d_cp *         &
                              ( q(k,j,i) - ql(k,j,i) ) / ( temp * temp ) )
                         k2 = theta * ( lv_d_cp / temp * k1 - 1.0_wp )
                      ENDIF
                      tmp_flux(k) = -1.0_wp * kh(k,j,i) *                      &
                                      ( k1 * ( pt(k+1,j,i) - pt(k-1,j,i) ) +   &
                                        k2 * ( q(k+1,j,i)  - q(k-1,j,i) )      &
                                      ) * dd2zu(k)
                   ELSE IF ( cloud_droplets )  THEN
                      k1 = 1.0_wp + 0.61_wp * q(k,j,i) - ql(k,j,i)
                      k2 = 0.61_wp * pt(k,j,i)
                      tmp_flux(k) = -1.0_wp * kh(k,j,i) * &
                                      ( k1 * ( pt(k+1,j,i) - pt(k-1,j,i) ) +   &
                                        k2 * ( q(k+1,j,i)  - q(k-1,j,i) ) -    &
                                        pt(k,j,i) * ( ql(k+1,j,i) -            &
                                        ql(k-1,j,i) ) ) * dd2zu(k)
                   ENDIF

                ENDDO

                IF ( use_surface_fluxes )  THEN
!
!--                Treat horizontal default surfaces
                   DO  l = 0, 1
                      surf_s = surf_def_h(l)%start_index(j,i)
                      surf_e = surf_def_h(l)%end_index(j,i)
                      DO  m = surf_s, surf_e
                         k = surf_def_h(l)%k(m)

                         IF ( .NOT. bulk_cloud_model .AND. .NOT. cloud_droplets ) THEN
                            k1 = 1.0_wp + 0.61_wp * q(k,j,i)
                            k2 = 0.61_wp * pt(k,j,i)
                         ELSE IF ( bulk_cloud_model )  THEN
                            IF ( ql(k,j,i) == 0.0_wp )  THEN
                               k1 = 1.0_wp + 0.61_wp * q(k,j,i)
                               k2 = 0.61_wp * pt(k,j,i)
                            ELSE
                               theta = pt(k,j,i) + d_exner(k) * lv_d_cp * ql(k,j,i)
                               temp  = theta * exner(k)
                               k1 = ( 1.0_wp - q(k,j,i) + 1.61_wp *            &
                                          ( q(k,j,i) - ql(k,j,i) ) *           &
                                 ( 1.0_wp + rd_d_rv * lv_d_rd / temp ) ) /     &
                                 ( 1.0_wp + rd_d_rv * lv_d_rd * lv_d_cp *      &
                                 ( q(k,j,i) - ql(k,j,i) ) / ( temp * temp ) )
                               k2 = theta * ( lv_d_cp / temp * k1 - 1.0_wp )
                            ENDIF
                         ELSE IF ( cloud_droplets )  THEN
                            k1 = 1.0_wp + 0.61_wp * q(k,j,i) - ql(k,j,i)
                            k2 = 0.61_wp * pt(k,j,i)
                         ENDIF

                         tmp_flux(k) = ( k1 * surf_def_h(l)%shf(m) +           &
                                         k2 * surf_def_h(l)%qsws(m)            &
                                       ) * drho_air_zw(k-1)
                      ENDDO
                   ENDDO
!
!--                Treat horizontal natural surfaces
                   surf_s = surf_lsm_h%start_index(j,i)
                   surf_e = surf_lsm_h%end_index(j,i)
                   DO  m = surf_s, surf_e
                      k = surf_lsm_h%k(m)

                      IF ( .NOT. bulk_cloud_model .AND. .NOT. cloud_droplets ) THEN
                         k1 = 1.0_wp + 0.61_wp * q(k,j,i)
                         k2 = 0.61_wp * pt(k,j,i)
                      ELSE IF ( bulk_cloud_model )  THEN
                         IF ( ql(k,j,i) == 0.0_wp )  THEN
                            k1 = 1.0_wp + 0.61_wp * q(k,j,i)
                            k2 = 0.61_wp * pt(k,j,i)
                         ELSE
                            theta = pt(k,j,i) + d_exner(k) * lv_d_cp * ql(k,j,i)
                            temp  = theta * exner(k)
                            k1 = ( 1.0_wp - q(k,j,i) + 1.61_wp *               &
                                          ( q(k,j,i) - ql(k,j,i) ) *           &
                                 ( 1.0_wp + rd_d_rv * lv_d_rd / temp ) ) /     &
                                 ( 1.0_wp + rd_d_rv * lv_d_rd * lv_d_cp *      &
                                 ( q(k,j,i) - ql(k,j,i) ) / ( temp * temp ) )
                            k2 = theta * ( lv_d_cp / temp * k1 - 1.0_wp )
                         ENDIF
                      ELSE IF ( cloud_droplets )  THEN
                         k1 = 1.0_wp + 0.61_wp * q(k,j,i) - ql(k,j,i)
                         k2 = 0.61_wp * pt(k,j,i)
                      ENDIF

                      tmp_flux(k) = ( k1 * surf_lsm_h%shf(m) +                 &
                                      k2 * surf_lsm_h%qsws(m)                  &
                                    ) * drho_air_zw(k-1)
                   ENDDO
!
!--                Treat horizontal urban surfaces
                   surf_s = surf_usm_h%start_index(j,i)
                   surf_e = surf_usm_h%end_index(j,i)
                   DO  m = surf_s, surf_e
                      k = surf_usm_h%k(m)

                      IF ( .NOT. bulk_cloud_model .AND. .NOT. cloud_droplets ) THEN
                         k1 = 1.0_wp + 0.61_wp * q(k,j,i)
                         k2 = 0.61_wp * pt(k,j,i)
                      ELSE IF ( bulk_cloud_model )  THEN
                         IF ( ql(k,j,i) == 0.0_wp )  THEN
                            k1 = 1.0_wp + 0.61_wp * q(k,j,i)
                            k2 = 0.61_wp * pt(k,j,i)
                         ELSE
                            theta = pt(k,j,i) + d_exner(k) * lv_d_cp * ql(k,j,i)
                            temp  = theta * exner(k)
                            k1 = ( 1.0_wp - q(k,j,i) + 1.61_wp *               &
                                          ( q(k,j,i) - ql(k,j,i) ) *           &
                                 ( 1.0_wp + rd_d_rv * lv_d_rd / temp ) ) /     &
                                 ( 1.0_wp + rd_d_rv * lv_d_rd * lv_d_cp *      &
                                 ( q(k,j,i) - ql(k,j,i) ) / ( temp * temp ) )
                            k2 = theta * ( lv_d_cp / temp * k1 - 1.0_wp )
                         ENDIF
                      ELSE IF ( cloud_droplets )  THEN
                         k1 = 1.0_wp + 0.61_wp * q(k,j,i) - ql(k,j,i)
                         k2 = 0.61_wp * pt(k,j,i)
                      ENDIF

                      tmp_flux(k) = ( k1 * surf_usm_h%shf(m) +                 &
                                      k2 * surf_usm_h%qsws(m)                  &
                                    ) * drho_air_zw(k-1)
                   ENDDO

                ENDIF ! from IF ( use_surface_fluxes )  THEN

                IF ( use_top_fluxes )  THEN

                   surf_s = surf_def_h(2)%start_index(j,i)
                   surf_e = surf_def_h(2)%end_index(j,i)
                   DO  m = surf_s, surf_e
                      k = surf_def_h(2)%k(m)

                      IF ( .NOT. bulk_cloud_model .AND. .NOT. cloud_droplets ) THEN
                         k1 = 1.0_wp + 0.61_wp * q(k,j,i)
                         k2 = 0.61_wp * pt(k,j,i)
                      ELSE IF ( bulk_cloud_model )  THEN
                         IF ( ql(k,j,i) == 0.0_wp )  THEN
                            k1 = 1.0_wp + 0.61_wp * q(k,j,i)
                            k2 = 0.61_wp * pt(k,j,i)
                         ELSE
                            theta = pt(k,j,i) + d_exner(k) * lv_d_cp * ql(k,j,i)
                            temp  = theta * exner(k)
                            k1 = ( 1.0_wp - q(k,j,i) + 1.61_wp *               &
                                       ( q(k,j,i) - ql(k,j,i) ) *              &
                              ( 1.0_wp + rd_d_rv * lv_d_rd / temp ) ) /        &
                              ( 1.0_wp + rd_d_rv * lv_d_rd * lv_d_cp *         &
                              ( q(k,j,i) - ql(k,j,i) ) / ( temp * temp ) )
                            k2 = theta * ( lv_d_cp / temp * k1 - 1.0_wp )
                         ENDIF
                      ELSE IF ( cloud_droplets )  THEN
                         k1 = 1.0_wp + 0.61_wp * q(k,j,i) - ql(k,j,i)
                         k2 = 0.61_wp * pt(k,j,i)
                      ENDIF

                      tmp_flux(k) = ( k1 * surf_def_h(2)%shf(m) +              &
                                      k2 * surf_def_h(2)%qsws(m)               &
                                    ) * drho_air_zw(k)

                   ENDDO

                ENDIF ! from IF ( use_top_fluxes )  THEN

                IF ( .NOT. diss_production )  THEN

!--                Compute tendency for TKE-production from shear
                   DO  k = nzb+1, nzt
                      flag = MERGE( 1.0_wp, 0.0_wp, BTEST(wall_flags_total_0(k,j,i),0) )
                      tend(k,j,i) = tend(k,j,i) + flag * tmp_flux(k) * ( g / &
                                    MERGE( vpt_reference, vpt(k,j,i),          &
                                           use_single_reference_value ) )
                   ENDDO

                ELSE

!--                RANS mode: Compute tendency for dissipation-rate-production from shear
                   DO  k = nzb+1, nzt
                      flag = MERGE( 1.0_wp, 0.0_wp, BTEST(wall_flags_total_0(k,j,i),0) )
                      tend(k,j,i) = tend(k,j,i) + flag * tmp_flux(k) * ( g /   &
                                    MERGE( vpt_reference, vpt(k,j,i),          &
                                           use_single_reference_value ) ) *    &
                                    diss(k,j,i)/( e(k,j,i) + 1.0E-20_wp ) *    &
                                    c_3
                   ENDDO

                ENDIF

             ENDDO
          ENDDO

       ENDIF

    ENDIF

 END SUBROUTINE production_e


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Production terms (shear + buoyancy) of the TKE.
!> Cache-optimized version
!> @warning The case with constant_flux_layer = F and use_surface_fluxes = T is
!>          not considered well!
!------------------------------------------------------------------------------!
 SUBROUTINE production_e_ij( i, j, diss_production )

    USE arrays_3d,                                                             &
        ONLY:  ddzw, dd2zu, drho_air_zw, q, ql, d_exner, exner

    USE control_parameters,                                                    &
        ONLY:  cloud_droplets, constant_flux_layer, neutral,                   &
               rho_reference, use_single_reference_value, use_surface_fluxes,  &
               use_top_fluxes

    USE grid_variables,                                                        &
        ONLY:  ddx, dx, ddy, dy

    USE bulk_cloud_model_mod,                                                  &
        ONLY:  bulk_cloud_model

    IMPLICIT NONE

    LOGICAL :: diss_production

    INTEGER(iwp) ::  i       !< running index x-direction
    INTEGER(iwp) ::  j       !< running index y-direction
    INTEGER(iwp) ::  k       !< running index z-direction
    INTEGER(iwp) ::  l       !< running index for different surface type orientation
    INTEGER(iwp) ::  m       !< running index surface elements
    INTEGER(iwp) ::  surf_e  !< end index of surface elements at given i-j position
    INTEGER(iwp) ::  surf_s  !< start index of surface elements at given i-j position
    INTEGER(iwp) ::  flag_nr !< number of masking flag

    REAL(wp)     ::  def         !< ( du_i/dx_j + du_j/dx_i ) * du_i/dx_j
    REAL(wp)     ::  flag        !< flag to mask topography
    REAL(wp)     ::  k1          !< temporary factor
    REAL(wp)     ::  k2          !< temporary factor
    REAL(wp)     ::  km_neutral  !< diffusion coefficient assuming neutral conditions - used to compute shear production at surfaces
    REAL(wp)     ::  theta       !< virtual potential temperature
    REAL(wp)     ::  temp        !< theta * Exner-function
    REAL(wp)     ::  sign_dir    !< sign of wall-tke flux, depending on wall orientation
    REAL(wp)     ::  usvs        !< momentum flux u"v"
    REAL(wp)     ::  vsus        !< momentum flux v"u"
    REAL(wp)     ::  wsus        !< momentum flux w"u"
    REAL(wp)     ::  wsvs        !< momentum flux w"v"

    REAL(wp), DIMENSION(nzb+1:nzt) ::  dudx  !< Gradient of u-component in x-direction
    REAL(wp), DIMENSION(nzb+1:nzt) ::  dudy  !< Gradient of u-component in y-direction
    REAL(wp), DIMENSION(nzb+1:nzt) ::  dudz  !< Gradient of u-component in z-direction
    REAL(wp), DIMENSION(nzb+1:nzt) ::  dvdx  !< Gradient of v-component in x-direction
    REAL(wp), DIMENSION(nzb+1:nzt) ::  dvdy  !< Gradient of v-component in y-direction
    REAL(wp), DIMENSION(nzb+1:nzt) ::  dvdz  !< Gradient of v-component in z-direction
    REAL(wp), DIMENSION(nzb+1:nzt) ::  dwdx  !< Gradient of w-component in x-direction
    REAL(wp), DIMENSION(nzb+1:nzt) ::  dwdy  !< Gradient of w-component in y-direction
    REAL(wp), DIMENSION(nzb+1:nzt) ::  dwdz  !< Gradient of w-component in z-direction
    REAL(wp), DIMENSION(nzb+1:nzt) ::  tmp_flux  !< temporary flux-array in z-direction



!
!-- Calculate TKE production by shear. Calculate gradients at all grid
!-- points first, gradients at surface-bounded grid points will be
!-- overwritten further below.
    DO  k = nzb+1, nzt

       dudx(k) =           ( u(k,j,i+1) - u(k,j,i)     ) * ddx
       dudy(k) = 0.25_wp * ( u(k,j+1,i) + u(k,j+1,i+1) -                 &
                             u(k,j-1,i) - u(k,j-1,i+1) ) * ddy
       dudz(k) = 0.5_wp  * ( u(k+1,j,i) + u(k+1,j,i+1) -                 &
                             u(k-1,j,i) - u(k-1,j,i+1) ) * dd2zu(k)

       dvdx(k) = 0.25_wp * ( v(k,j,i+1) + v(k,j+1,i+1) -                 &
                             v(k,j,i-1) - v(k,j+1,i-1) ) * ddx
       dvdy(k) =           ( v(k,j+1,i) - v(k,j,i)     ) * ddy
       dvdz(k) = 0.5_wp  * ( v(k+1,j,i) + v(k+1,j+1,i) -                 &
                             v(k-1,j,i) - v(k-1,j+1,i) ) * dd2zu(k)

       dwdx(k) = 0.25_wp * ( w(k,j,i+1) + w(k-1,j,i+1) -                 &
                             w(k,j,i-1) - w(k-1,j,i-1) ) * ddx
       dwdy(k) = 0.25_wp * ( w(k,j+1,i) + w(k-1,j+1,i) -                 &
                             w(k,j-1,i) - w(k-1,j-1,i) ) * ddy
       dwdz(k) =           ( w(k,j,i)   - w(k-1,j,i)   ) * ddzw(k)

    ENDDO

    flag_nr = 29

    IF ( constant_flux_layer )  THEN

       flag_nr = 0

!--    Position beneath wall
!--    (2) - Will allways be executed.
!--    'bottom and wall: use u_0,v_0 and wall functions'
!
!--    Compute gradients at north- and south-facing surfaces.
!--    First, for default surfaces, then for urban surfaces.
!--    Note, so far no natural vertical surfaces implemented
       DO  l = 0, 1
          surf_s = surf_def_v(l)%start_index(j,i)
          surf_e = surf_def_v(l)%end_index(j,i)
          DO  m = surf_s, surf_e
             k           = surf_def_v(l)%k(m)
             usvs        = surf_def_v(l)%mom_flux_tke(0,m)
             wsvs        = surf_def_v(l)%mom_flux_tke(1,m)

             km_neutral = kappa * ( usvs**2 + wsvs**2 )**0.25_wp         &
                             * 0.5_wp * dy
!
!--          -1.0 for right-facing wall, 1.0 for left-facing wall
             sign_dir = MERGE( 1.0_wp, -1.0_wp,                          &
                        BTEST( wall_flags_total_0(k,j-1,i), flag_nr ) )
             dudy(k) = sign_dir * usvs / ( km_neutral + 1E-10_wp )
             dwdy(k) = sign_dir * wsvs / ( km_neutral + 1E-10_wp )
          ENDDO
!
!--       Natural surfaces
          surf_s = surf_lsm_v(l)%start_index(j,i)
          surf_e = surf_lsm_v(l)%end_index(j,i)
          DO  m = surf_s, surf_e
             k           = surf_lsm_v(l)%k(m)
             usvs        = surf_lsm_v(l)%mom_flux_tke(0,m)
             wsvs        = surf_lsm_v(l)%mom_flux_tke(1,m)

             km_neutral = kappa * ( usvs**2 + wsvs**2 )**0.25_wp         &
                             * 0.5_wp * dy
!
!--          -1.0 for right-facing wall, 1.0 for left-facing wall
             sign_dir = MERGE( 1.0_wp, -1.0_wp,                          &
                        BTEST( wall_flags_total_0(k,j-1,i), flag_nr ) )
             dudy(k) = sign_dir * usvs / ( km_neutral + 1E-10_wp )
             dwdy(k) = sign_dir * wsvs / ( km_neutral + 1E-10_wp )
          ENDDO
!
!--       Urban surfaces
          surf_s = surf_usm_v(l)%start_index(j,i)
          surf_e = surf_usm_v(l)%end_index(j,i)
          DO  m = surf_s, surf_e
             k           = surf_usm_v(l)%k(m)
             usvs        = surf_usm_v(l)%mom_flux_tke(0,m)
             wsvs        = surf_usm_v(l)%mom_flux_tke(1,m)

             km_neutral = kappa * ( usvs**2 + wsvs**2 )**0.25_wp         &
                             * 0.5_wp * dy
!
!--          -1.0 for right-facing wall, 1.0 for left-facing wall
             sign_dir = MERGE( 1.0_wp, -1.0_wp,                          &
                        BTEST( wall_flags_total_0(k,j-1,i), flag_nr ) )
             dudy(k) = sign_dir * usvs / ( km_neutral + 1E-10_wp )
             dwdy(k) = sign_dir * wsvs / ( km_neutral + 1E-10_wp )
          ENDDO
       ENDDO
!
!--    Compute gradients at east- and west-facing walls
       DO  l = 2, 3
          surf_s = surf_def_v(l)%start_index(j,i)
          surf_e = surf_def_v(l)%end_index(j,i)
          DO  m = surf_s, surf_e
             k     = surf_def_v(l)%k(m)
             vsus  = surf_def_v(l)%mom_flux_tke(0,m)
             wsus  = surf_def_v(l)%mom_flux_tke(1,m)

             km_neutral = kappa * ( vsus**2 + wsus**2 )**0.25_wp         &
                                * 0.5_wp * dx
!
!--          -1.0 for right-facing wall, 1.0 for left-facing wall
             sign_dir = MERGE( 1.0_wp, -1.0_wp,                          &
                        BTEST( wall_flags_total_0(k,j,i-1), flag_nr ) )
             dvdx(k) = sign_dir * vsus / ( km_neutral + 1E-10_wp )
             dwdx(k) = sign_dir * wsus / ( km_neutral + 1E-10_wp )
          ENDDO
!
!--       Natural surfaces
          surf_s = surf_lsm_v(l)%start_index(j,i)
          surf_e = surf_lsm_v(l)%end_index(j,i)
          DO  m = surf_s, surf_e
             k     = surf_lsm_v(l)%k(m)
             vsus  = surf_lsm_v(l)%mom_flux_tke(0,m)
             wsus  = surf_lsm_v(l)%mom_flux_tke(1,m)

             km_neutral = kappa * ( vsus**2 + wsus**2 )**0.25_wp         &
                                * 0.5_wp * dx
!
!--          -1.0 for right-facing wall, 1.0 for left-facing wall
             sign_dir = MERGE( 1.0_wp, -1.0_wp,                          &
                        BTEST( wall_flags_total_0(k,j,i-1), flag_nr ) )
             dvdx(k) = sign_dir * vsus / ( km_neutral + 1E-10_wp )
             dwdx(k) = sign_dir * wsus / ( km_neutral + 1E-10_wp )
          ENDDO
!
!--       Urban surfaces
          surf_s = surf_usm_v(l)%start_index(j,i)
          surf_e = surf_usm_v(l)%end_index(j,i)
          DO  m = surf_s, surf_e
             k     = surf_usm_v(l)%k(m)
             vsus  = surf_usm_v(l)%mom_flux_tke(0,m)
             wsus  = surf_usm_v(l)%mom_flux_tke(1,m)

             km_neutral = kappa * ( vsus**2 + wsus**2 )**0.25_wp         &
                                * 0.5_wp * dx
!
!--          -1.0 for right-facing wall, 1.0 for left-facing wall
             sign_dir = MERGE( 1.0_wp, -1.0_wp,                          &
                        BTEST( wall_flags_total_0(k,j,i-1), flag_nr ) )
             dvdx(k) = sign_dir * vsus / ( km_neutral + 1E-10_wp )
             dwdx(k) = sign_dir * wsus / ( km_neutral + 1E-10_wp )
          ENDDO
       ENDDO
!
!--    Compute gradients at upward-facing surfaces
       surf_s = surf_def_h(0)%start_index(j,i)
       surf_e = surf_def_h(0)%end_index(j,i)
       DO  m = surf_s, surf_e
          k = surf_def_h(0)%k(m)
!
!--       Please note, actually, an interpolation of u_0 and v_0
!--       onto the grid center would be required. However, this
!--       would require several data transfers between 2D-grid and
!--       wall type. The effect of this missing interpolation is
!--       negligible. (See also production_e_init).
          dudz(k) = ( u(k+1,j,i) - surf_def_h(0)%u_0(m) ) * dd2zu(k)
          dvdz(k) = ( v(k+1,j,i) - surf_def_h(0)%v_0(m) ) * dd2zu(k)

       ENDDO
!
!--    Natural surfaces
       surf_s = surf_lsm_h%start_index(j,i)
       surf_e = surf_lsm_h%end_index(j,i)
       DO  m = surf_s, surf_e
          k = surf_lsm_h%k(m)

          dudz(k) = ( u(k+1,j,i) - surf_lsm_h%u_0(m) ) * dd2zu(k)
          dvdz(k) = ( v(k+1,j,i) - surf_lsm_h%v_0(m) ) * dd2zu(k)

       ENDDO
!
!--    Urban surfaces
       surf_s = surf_usm_h%start_index(j,i)
       surf_e = surf_usm_h%end_index(j,i)
       DO  m = surf_s, surf_e
          k = surf_usm_h%k(m)

          dudz(k) = ( u(k+1,j,i) - surf_usm_h%u_0(m) ) * dd2zu(k)
          dvdz(k) = ( v(k+1,j,i) - surf_usm_h%v_0(m) ) * dd2zu(k)

       ENDDO
!
!--    Compute gradients at downward-facing walls, only for
!--    non-natural default surfaces
       surf_s = surf_def_h(1)%start_index(j,i)
       surf_e = surf_def_h(1)%end_index(j,i)
       DO  m = surf_s, surf_e
          k = surf_def_h(1)%k(m)

          dudz(k) = ( surf_def_h(1)%u_0(m) - u(k-1,j,i) ) * dd2zu(k)
          dvdz(k) = ( surf_def_h(1)%v_0(m) - v(k-1,j,i) ) * dd2zu(k)

       ENDDO

    ENDIF

    DO  k = nzb+1, nzt

       def = 2.0_wp * ( dudx(k)**2 + dvdy(k)**2 + dwdz(k)**2 ) +         &
                        dudy(k)**2 + dvdx(k)**2 + dwdx(k)**2 +           &
                        dwdy(k)**2 + dudz(k)**2 + dvdz(k)**2 +           &
             2.0_wp * ( dvdx(k)*dudy(k) + dwdx(k)*dudz(k) +              &
                        dwdy(k)*dvdz(k) )

       IF ( def < 0.0_wp )  def = 0.0_wp

       flag = MERGE( 1.0_wp, 0.0_wp, BTEST(wall_flags_total_0(k,j,i),flag_nr) )

       IF ( .NOT. diss_production )  THEN

!--       Compute tendency for TKE-production from shear
          tend(k,j,i) = tend(k,j,i) + km(k,j,i) * def * flag

       ELSE

!--       RANS mode: Compute tendency for dissipation-rate-production from shear
          tend(k,j,i) = tend(k,j,i) + km(k,j,i) * def * flag *           &
                        diss(k,j,i)/( e(k,j,i) + 1.0E-20_wp ) * c_1

       ENDIF

    ENDDO

!
!-- If required, calculate TKE production by buoyancy
    IF ( .NOT. neutral )  THEN

       IF ( .NOT. humidity )  THEN

          IF ( ocean_mode )  THEN
!
!--          So far in the ocean no special treatment of density flux
!--          in the bottom and top surface layer
             DO  k = nzb+1, nzt
                tmp_flux(k) = kh(k,j,i) * ( prho(k+1,j,i) - prho(k-1,j,i) ) * dd2zu(k)
             ENDDO
!
!--          Treatment of near-surface grid points, at up- and down-
!--          ward facing surfaces
             IF ( use_surface_fluxes )  THEN
                DO  l = 0, 1
                   surf_s = surf_def_h(l)%start_index(j,i)
                   surf_e = surf_def_h(l)%end_index(j,i)
                   DO  m = surf_s, surf_e
                      k = surf_def_h(l)%k(m)
                      tmp_flux(k) = drho_air_zw(k-1) * surf_def_h(l)%shf(m)
                   ENDDO
                ENDDO
             ENDIF

             IF ( use_top_fluxes )  THEN
                surf_s = surf_def_h(2)%start_index(j,i)
                surf_e = surf_def_h(2)%end_index(j,i)
                DO  m = surf_s, surf_e
                   k = surf_def_h(2)%k(m)
                   tmp_flux(k) = drho_air_zw(k) * surf_def_h(2)%shf(m)
                ENDDO
             ENDIF

             IF ( .NOT. diss_production )  THEN

!--             Compute tendency for TKE-production from shear
                DO  k = nzb+1, nzt
                   flag = MERGE( 1.0_wp, 0.0_wp, BTEST(wall_flags_total_0(k,j,i),0) )
                   tend(k,j,i) = tend(k,j,i) + flag * tmp_flux(k) * ( g / &
                                 MERGE( rho_reference, prho(k,j,i),       &
                                        use_single_reference_value ) )
                ENDDO

             ELSE

!--             RANS mode: Compute tendency for dissipation-rate-production from shear
                DO  k = nzb+1, nzt
                   flag = MERGE( 1.0_wp, 0.0_wp, BTEST(wall_flags_total_0(k,j,i),0) )
                   tend(k,j,i) = tend(k,j,i) + flag * tmp_flux(k) * ( g / &
                                 MERGE( rho_reference, prho(k,j,i),       &
                                        use_single_reference_value ) ) *  &
                                 diss(k,j,i)/( e(k,j,i) + 1.0E-20_wp ) *  &
                                 c_3
                ENDDO

             ENDIF


          ELSE ! or IF ( .NOT. ocean_mode )  THEN

             DO  k = nzb+1, nzt
                tmp_flux(k) = -1.0_wp * kh(k,j,i) * ( pt(k+1,j,i) - pt(k-1,j,i) ) * dd2zu(k)
             ENDDO

             IF ( use_surface_fluxes )  THEN
!
!--             Default surfaces, up- and downward-facing
                DO  l = 0, 1
                   surf_s = surf_def_h(l)%start_index(j,i)
                   surf_e = surf_def_h(l)%end_index(j,i)
                   DO  m = surf_s, surf_e
                      k = surf_def_h(l)%k(m)
                      tmp_flux(k) = drho_air_zw(k-1) * surf_def_h(l)%shf(m)
                   ENDDO
                ENDDO
!
!--             Natural surfaces
                surf_s = surf_lsm_h%start_index(j,i)
                surf_e = surf_lsm_h%end_index(j,i)
                DO  m = surf_s, surf_e
                   k = surf_lsm_h%k(m)
                   tmp_flux(k) = drho_air_zw(k-1) * surf_lsm_h%shf(m)
                ENDDO
!
!--             Urban surfaces
                surf_s = surf_usm_h%start_index(j,i)
                surf_e = surf_usm_h%end_index(j,i)
                DO  m = surf_s, surf_e
                   k = surf_usm_h%k(m)
                   tmp_flux(k) = drho_air_zw(k-1) * surf_usm_h%shf(m)
                ENDDO
             ENDIF

             IF ( use_top_fluxes )  THEN
                surf_s = surf_def_h(2)%start_index(j,i)
                surf_e = surf_def_h(2)%end_index(j,i)
                DO  m = surf_s, surf_e
                   k = surf_def_h(2)%k(m)
                   tmp_flux(k) = drho_air_zw(k) * surf_def_h(2)%shf(m)
                ENDDO
             ENDIF

             IF ( .NOT. diss_production )  THEN

!--             Compute tendency for TKE-production from shear
                DO  k = nzb+1, nzt
                   flag = MERGE( 1.0_wp, 0.0_wp, BTEST(wall_flags_total_0(k,j,i),0) )
                   tend(k,j,i) = tend(k,j,i) + flag * tmp_flux(k) * ( g / &
                                 MERGE( pt_reference, pt(k,j,i),          &
                                        use_single_reference_value ) )
                ENDDO

             ELSE

!--             RANS mode: Compute tendency for dissipation-rate-production from shear
                DO  k = nzb+1, nzt
                   flag = MERGE( 1.0_wp, 0.0_wp, BTEST(wall_flags_total_0(k,j,i),0) )
                   tend(k,j,i) = tend(k,j,i) + flag * tmp_flux(k) * ( g / &
                                 MERGE( pt_reference, pt(k,j,i),          &
                                        use_single_reference_value ) ) *  &
                                 diss(k,j,i)/( e(k,j,i) + 1.0E-20_wp ) *  &
                                 c_3
                ENDDO

             ENDIF

          ENDIF ! from IF ( .NOT. ocean_mode )

       ELSE ! or IF ( humidity )  THEN

          DO  k = nzb+1, nzt

             IF ( .NOT. bulk_cloud_model .AND. .NOT. cloud_droplets ) THEN
                k1 = 1.0_wp + 0.61_wp * q(k,j,i)
                k2 = 0.61_wp * pt(k,j,i)
                tmp_flux(k) = -1.0_wp * kh(k,j,i) *                      &
                                ( k1 * ( pt(k+1,j,i) - pt(k-1,j,i) ) +   &
                                  k2 * ( q(k+1,j,i)  - q(k-1,j,i) )      &
                                ) * dd2zu(k)
             ELSE IF ( bulk_cloud_model )  THEN
                IF ( ql(k,j,i) == 0.0_wp )  THEN
                   k1 = 1.0_wp + 0.61_wp * q(k,j,i)
                   k2 = 0.61_wp * pt(k,j,i)
                ELSE
                   theta = pt(k,j,i) + d_exner(k) * lv_d_cp * ql(k,j,i)
                   temp  = theta * exner(k)
                   k1 = ( 1.0_wp - q(k,j,i) + 1.61_wp *                  &
                                 ( q(k,j,i) - ql(k,j,i) ) *              &
                        ( 1.0_wp + rd_d_rv * lv_d_rd / temp ) ) /        &
                        ( 1.0_wp + rd_d_rv * lv_d_rd * lv_d_cp *         &
                        ( q(k,j,i) - ql(k,j,i) ) / ( temp * temp ) )
                   k2 = theta * ( lv_d_cp / temp * k1 - 1.0_wp )
                ENDIF
                tmp_flux(k) = -1.0_wp * kh(k,j,i) *                      &
                                ( k1 * ( pt(k+1,j,i) - pt(k-1,j,i) ) +   &
                                  k2 * ( q(k+1,j,i)  - q(k-1,j,i) )      &
                                ) * dd2zu(k)
             ELSE IF ( cloud_droplets )  THEN
                k1 = 1.0_wp + 0.61_wp * q(k,j,i) - ql(k,j,i)
                k2 = 0.61_wp * pt(k,j,i)
                tmp_flux(k) = -1.0_wp * kh(k,j,i) * &
                                ( k1 * ( pt(k+1,j,i) - pt(k-1,j,i) ) +   &
                                  k2 * ( q(k+1,j,i)  - q(k-1,j,i) ) -    &
                                  pt(k,j,i) * ( ql(k+1,j,i) -            &
                                  ql(k-1,j,i) ) ) * dd2zu(k)
             ENDIF

          ENDDO

          IF ( use_surface_fluxes )  THEN
!
!--          Treat horizontal default surfaces
             DO  l = 0, 1
                surf_s = surf_def_h(l)%start_index(j,i)
                surf_e = surf_def_h(l)%end_index(j,i)
                DO  m = surf_s, surf_e
                   k = surf_def_h(l)%k(m)

                   IF ( .NOT. bulk_cloud_model .AND. .NOT. cloud_droplets ) THEN
                      k1 = 1.0_wp + 0.61_wp * q(k,j,i)
                      k2 = 0.61_wp * pt(k,j,i)
                   ELSE IF ( bulk_cloud_model )  THEN
                      IF ( ql(k,j,i) == 0.0_wp )  THEN
                         k1 = 1.0_wp + 0.61_wp * q(k,j,i)
                         k2 = 0.61_wp * pt(k,j,i)
                      ELSE
                         theta = pt(k,j,i) + d_exner(k) * lv_d_cp * ql(k,j,i)
                         temp  = theta * exner(k)
                         k1 = ( 1.0_wp - q(k,j,i) + 1.61_wp *            &
                                    ( q(k,j,i) - ql(k,j,i) ) *           &
                           ( 1.0_wp + rd_d_rv * lv_d_rd / temp ) ) /     &
                           ( 1.0_wp + rd_d_rv * lv_d_rd * lv_d_cp *      &
                           ( q(k,j,i) - ql(k,j,i) ) / ( temp * temp ) )
                         k2 = theta * ( lv_d_cp / temp * k1 - 1.0_wp )
                      ENDIF
                   ELSE IF ( cloud_droplets )  THEN
                      k1 = 1.0_wp + 0.61_wp * q(k,j,i) - ql(k,j,i)
                      k2 = 0.61_wp * pt(k,j,i)
                   ENDIF

                   tmp_flux(k) = ( k1 * surf_def_h(l)%shf(m) +           &
                                   k2 * surf_def_h(l)%qsws(m)            &
                                 ) * drho_air_zw(k-1)
                ENDDO
             ENDDO
!
!--          Treat horizontal natural surfaces
             surf_s = surf_lsm_h%start_index(j,i)
             surf_e = surf_lsm_h%end_index(j,i)
             DO  m = surf_s, surf_e
                k = surf_lsm_h%k(m)

                IF ( .NOT. bulk_cloud_model .AND. .NOT. cloud_droplets ) THEN
                   k1 = 1.0_wp + 0.61_wp * q(k,j,i)
                   k2 = 0.61_wp * pt(k,j,i)
                ELSE IF ( bulk_cloud_model )  THEN
                   IF ( ql(k,j,i) == 0.0_wp )  THEN
                      k1 = 1.0_wp + 0.61_wp * q(k,j,i)
                      k2 = 0.61_wp * pt(k,j,i)
                   ELSE
                      theta = pt(k,j,i) + d_exner(k) * lv_d_cp * ql(k,j,i)
                      temp  = theta * exner(k)
                      k1 = ( 1.0_wp - q(k,j,i) + 1.61_wp *               &
                                    ( q(k,j,i) - ql(k,j,i) ) *           &
                           ( 1.0_wp + rd_d_rv * lv_d_rd / temp ) ) /     &
                           ( 1.0_wp + rd_d_rv * lv_d_rd * lv_d_cp *      &
                           ( q(k,j,i) - ql(k,j,i) ) / ( temp * temp ) )
                      k2 = theta * ( lv_d_cp / temp * k1 - 1.0_wp )
                   ENDIF
                ELSE IF ( cloud_droplets )  THEN
                   k1 = 1.0_wp + 0.61_wp * q(k,j,i) - ql(k,j,i)
                   k2 = 0.61_wp * pt(k,j,i)
                ENDIF

                tmp_flux(k) = ( k1 * surf_lsm_h%shf(m) +                 &
                                k2 * surf_lsm_h%qsws(m)                  &
                              ) * drho_air_zw(k-1)
             ENDDO
!
!--          Treat horizontal urban surfaces
             surf_s = surf_usm_h%start_index(j,i)
             surf_e = surf_usm_h%end_index(j,i)
             DO  m = surf_s, surf_e
                k = surf_usm_h%k(m)

                IF ( .NOT. bulk_cloud_model .AND. .NOT. cloud_droplets ) THEN
                   k1 = 1.0_wp + 0.61_wp * q(k,j,i)
                   k2 = 0.61_wp * pt(k,j,i)
                ELSE IF ( bulk_cloud_model )  THEN
                   IF ( ql(k,j,i) == 0.0_wp )  THEN
                      k1 = 1.0_wp + 0.61_wp * q(k,j,i)
                      k2 = 0.61_wp * pt(k,j,i)
                   ELSE
                      theta = pt(k,j,i) + d_exner(k) * lv_d_cp * ql(k,j,i)
                      temp  = theta * exner(k)
                      k1 = ( 1.0_wp - q(k,j,i) + 1.61_wp *               &
                                    ( q(k,j,i) - ql(k,j,i) ) *           &
                           ( 1.0_wp + rd_d_rv * lv_d_rd / temp ) ) /     &
                           ( 1.0_wp + rd_d_rv * lv_d_rd * lv_d_cp *      &
                           ( q(k,j,i) - ql(k,j,i) ) / ( temp * temp ) )
                      k2 = theta * ( lv_d_cp / temp * k1 - 1.0_wp )
                   ENDIF
                ELSE IF ( cloud_droplets )  THEN
                   k1 = 1.0_wp + 0.61_wp * q(k,j,i) - ql(k,j,i)
                   k2 = 0.61_wp * pt(k,j,i)
                ENDIF

                tmp_flux(k) = ( k1 * surf_usm_h%shf(m) +                 &
                                k2 * surf_usm_h%qsws(m)                  &
                              ) * drho_air_zw(k-1)
             ENDDO

          ENDIF ! from IF ( use_surface_fluxes )  THEN

          IF ( use_top_fluxes )  THEN

             surf_s = surf_def_h(2)%start_index(j,i)
             surf_e = surf_def_h(2)%end_index(j,i)
             DO  m = surf_s, surf_e
                k = surf_def_h(2)%k(m)

                IF ( .NOT. bulk_cloud_model .AND. .NOT. cloud_droplets ) THEN
                   k1 = 1.0_wp + 0.61_wp * q(k,j,i)
                   k2 = 0.61_wp * pt(k,j,i)
                ELSE IF ( bulk_cloud_model )  THEN
                   IF ( ql(k,j,i) == 0.0_wp )  THEN
                      k1 = 1.0_wp + 0.61_wp * q(k,j,i)
                      k2 = 0.61_wp * pt(k,j,i)
                   ELSE
                      theta = pt(k,j,i) + d_exner(k) * lv_d_cp * ql(k,j,i)
                      temp  = theta * exner(k)
                      k1 = ( 1.0_wp - q(k,j,i) + 1.61_wp *               &
                                 ( q(k,j,i) - ql(k,j,i) ) *              &
                        ( 1.0_wp + rd_d_rv * lv_d_rd / temp ) ) /        &
                        ( 1.0_wp + rd_d_rv * lv_d_rd * lv_d_cp *         &
                        ( q(k,j,i) - ql(k,j,i) ) / ( temp * temp ) )
                      k2 = theta * ( lv_d_cp / temp * k1 - 1.0_wp )
                   ENDIF
                ELSE IF ( cloud_droplets )  THEN
                   k1 = 1.0_wp + 0.61_wp * q(k,j,i) - ql(k,j,i)
                   k2 = 0.61_wp * pt(k,j,i)
                ENDIF

                tmp_flux(k) = ( k1 * surf_def_h(2)%shf(m) +              &
                                k2 * surf_def_h(2)%qsws(m)               &
                              ) * drho_air_zw(k)

             ENDDO

          ENDIF ! from IF ( use_top_fluxes )  THEN

          IF ( .NOT. diss_production )  THEN

!--          Compute tendency for TKE-production from shear
             DO  k = nzb+1, nzt
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST(wall_flags_total_0(k,j,i),0) )
                tend(k,j,i) = tend(k,j,i) + flag * tmp_flux(k) * ( g / &
                              MERGE( vpt_reference, vpt(k,j,i),          &
                                     use_single_reference_value ) )
             ENDDO

          ELSE

!--          RANS mode: Compute tendency for dissipation-rate-production from shear
             DO  k = nzb+1, nzt
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST(wall_flags_total_0(k,j,i),0) )
                tend(k,j,i) = tend(k,j,i) + flag * tmp_flux(k) * ( g /   &
                              MERGE( vpt_reference, vpt(k,j,i),          &
                                     use_single_reference_value ) ) *    &
                              diss(k,j,i)/( e(k,j,i) + 1.0E-20_wp ) *    &
                              c_3
             ENDDO

          ENDIF

       ENDIF

    ENDIF

 END SUBROUTINE production_e_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Diffusion and dissipation terms for the TKE.
!> Vector-optimized version
!> @todo Try to avoid the usage of the 3d-array 'diss' where possible (case les
!>       and rans_tke_l if not wang_kernel, use_sgs_for_particles, or
!>       collision_turbulence).
!------------------------------------------------------------------------------!
 SUBROUTINE diffusion_e( var, var_reference )

    USE arrays_3d,                                                             &
        ONLY:  dd2zu, ddzu, ddzw, drho_air, rho_air_zw

    USE control_parameters,                                                    &
        ONLY:  atmos_ocean_sign, use_single_reference_value

    USE grid_variables,                                                        &
        ONLY:  ddx2, ddy2

    USE bulk_cloud_model_mod,                                                  &
        ONLY:  collision_turbulence

    USE particle_attributes,                                                   &
        ONLY:  use_sgs_for_particles, wang_kernel

    IMPLICIT NONE

    INTEGER(iwp) ::  i              !< running index x direction
    INTEGER(iwp) ::  j              !< running index y direction
    INTEGER(iwp) ::  k              !< running index z direction
    INTEGER(iwp) ::  m              !< running index surface elements

    REAL(wp)     ::  duv2_dz2       !< squared vertical gradient of wind vector
    REAL(wp)     ::  dvar_dz        !< vertical gradient of var
    REAL(wp)     ::  l              !< mixing length
    REAL(wp)     ::  var_reference  !< reference temperature

    REAL(wp), DIMENSION(nzb+1:nzt) ::  l_stable  !< mixing length according to stratification
    REAL(wp), DIMENSION(nzb+1:nzt) ::  rif       !< Richardson flux number

    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg), INTENT(IN) ::  var  !< temperature

!
!-- Calculate the dissipation
    IF ( les_dynamic .OR. les_mw )  THEN

       !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(i, j) &
       !$ACC PRIVATE(l, l_stable, dvar_dz) &
       !$ACC PRESENT(diss, e, var, wall_flags_total_0) &
       !$ACC PRESENT(dd2zu, l_grid, l_wall)
       DO  i = nxl, nxr
          DO  j = nys, nyn
             !$ACC LOOP PRIVATE(k)
             DO  k = nzb+1, nzt

                dvar_dz = atmos_ocean_sign * ( var(k+1,j,i) - var(k-1,j,i) ) * dd2zu(k)
                IF ( dvar_dz > 0.0_wp ) THEN
                   IF ( use_single_reference_value )  THEN
                     l_stable(k) = 0.76_wp * SQRT( e(k,j,i) )                  &
                                 / SQRT( g / var_reference * dvar_dz ) + 1E-5_wp
                   ELSE
                     l_stable(k) = 0.76_wp * SQRT( e(k,j,i) )               &
                                 / SQRT( g / var(k,j,i) * dvar_dz ) + 1E-5_wp
                   ENDIF
                ELSE
                   l_stable(k) = l_grid(k)
                ENDIF

             ENDDO

             !$ACC LOOP PRIVATE(k)
             !DIR$ IVDEP
             DO  k = nzb+1, nzt

                l  = MIN( l_wall(k,j,i), l_stable(k) )

                diss(k,j,i) = ( 0.19_wp + 0.74_wp * l / l_wall(k,j,i) )                &
                              * e(k,j,i) * SQRT( e(k,j,i) ) / l                        &
                              * MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )

             ENDDO

          ENDDO
       ENDDO

    ELSEIF ( rans_tke_l )  THEN

      !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(i, j) &
      !$ACC PRIVATE(l_stable, duv2_dz2, rif, dvar_dz) &
      !$ACC PRESENT(diss, e, u, v, var, wall_flags_total_0) &
      !$ACC PRESENT(dd2zu, l_black, l_wall)
       DO  i = nxl, nxr
          DO  j = nys, nyn
!
!--          Calculate Richardson-flux number
             IF ( use_single_reference_value )  THEN
                !$ACC LOOP PRIVATE(k)
                DO  k = nzb+1, nzt

                   dvar_dz = atmos_ocean_sign * ( var(k+1,j,i) - var(k-1,j,i) ) * dd2zu(k)

                   duv2_dz2 =   ( ( u(k+1,j,i) - u(k-1,j,i) ) * dd2zu(k) )**2  &
                              + ( ( v(k+1,j,i) - v(k-1,j,i) ) * dd2zu(k) )**2  &
                              + 1E-30_wp

                   rif(k) = MIN( MAX( g / var_reference * dvar_dz / duv2_dz2, -5.0_wp ),  1.0_wp )
                ENDDO
             ELSE
                !$ACC LOOP PRIVATE(k)
                DO  k = nzb+1, nzt

                   dvar_dz = atmos_ocean_sign * ( var(k+1,j,i) - var(k-1,j,i) ) * dd2zu(k)

                   duv2_dz2 =   ( ( u(k+1,j,i) - u(k-1,j,i) ) * dd2zu(k) )**2  &
                              + ( ( v(k+1,j,i) - v(k-1,j,i) ) * dd2zu(k) )**2  &
                              + 1E-30_wp

                   rif(k) = MIN( MAX( g / var(k,j,i) * dvar_dz / duv2_dz2, -5.0_wp ),  1.0_wp )
                ENDDO
             ENDIF
!
!--          Calculate diabatic mixing length using Dyer-profile functions
             !$ACC LOOP PRIVATE(k)
             DO  k = nzb+1, nzt
                IF ( rif(k) >= 0.0_wp )  THEN
                   l_stable(k) = MIN( l_black(k) / ( 1.0_wp + 5.0_wp * rif(k) ), l_wall(k,j,i) )
                ELSE
                   l_stable(k) = l_wall(k,j,i) * SQRT( 1.0_wp - 16.0_wp * rif(k) )
                ENDIF
             ENDDO

             !$ACC LOOP PRIVATE(k)
             !DIR$ IVDEP
             DO  k = nzb+1, nzt

                diss(k,j,i) = c_0**3 * e(k,j,i) * SQRT( e(k,j,i) ) / l_stable(k)     &
                            * MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )

             ENDDO

          ENDDO
       ENDDO

!-- Note, in case of rans_tke_e, the dissipation is already calculated
!-- in prognostic_equations
    ENDIF

    !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i, j, k) &
    !$ACC PRESENT(diss, e, km, tend, wall_flags_total_0) &
    !$ACC PRESENT(ddzu, ddzw, rho_air_zw, drho_air)
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt

             tend(k,j,i) = tend(k,j,i) + (                                            &
                                           (                                          &
                       ( km(k,j,i)+km(k,j,i+1) ) * ( e(k,j,i+1)-e(k,j,i) )            &
                     - ( km(k,j,i)+km(k,j,i-1) ) * ( e(k,j,i)-e(k,j,i-1) )            &
                                           ) * ddx2                                   &
                                         + (                                          &
                       ( km(k,j,i)+km(k,j+1,i) ) * ( e(k,j+1,i)-e(k,j,i) )            &
                     - ( km(k,j,i)+km(k,j-1,i) ) * ( e(k,j,i)-e(k,j-1,i) )            &
                                           ) * ddy2                                   &
                                         + (                                          &
            ( km(k,j,i)+km(k+1,j,i) ) * ( e(k+1,j,i)-e(k,j,i) ) * ddzu(k+1)           &
                                                          * rho_air_zw(k)             &
          - ( km(k,j,i)+km(k-1,j,i) ) * ( e(k,j,i)-e(k-1,j,i) ) * ddzu(k)             &
                                                          * rho_air_zw(k-1)           &
                                           ) * ddzw(k) * drho_air(k)                  &
                                         ) * dsig_e                                   &
                                           * MERGE( 1.0_wp, 0.0_wp,                   &
                                             BTEST( wall_flags_total_0(k,j,i), 0 ) )  &
                          - diss(k,j,i)

          ENDDO
       ENDDO
    ENDDO

!
!-- Neumann boundary condition for dissipation diss(nzb,:,:) = diss(nzb+1,:,:).
!-- Note, bc cannot be set in tcm_boundary conditions as the dissipation
!-- in LES mode is only a diagnostic quantity.
    IF ( .NOT. rans_tke_e .AND. ( use_sgs_for_particles  .OR.                  &
         wang_kernel  .OR.  collision_turbulence  ) )  THEN
!
!--    Upward facing surfaces
       DO  m = 1, bc_h(0)%ns
          i = bc_h(0)%i(m)
          j = bc_h(0)%j(m)
          k = bc_h(0)%k(m)
          diss(k-1,j,i) = diss(k,j,i)
       ENDDO
!
!--    Downward facing surfaces
       DO  m = 1, bc_h(1)%ns
          i = bc_h(1)%i(m)
          j = bc_h(1)%j(m)
          k = bc_h(1)%k(m)
          diss(k+1,j,i) = diss(k,j,i)
       ENDDO

    ENDIF

 END SUBROUTINE diffusion_e


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Diffusion and dissipation terms for the TKE.
!> Cache-optimized version
!> @todo Try to avoid the usage of the 3d-array 'diss' where possible (case les
!>       and rans_tke_l if not wang_kernel, use_sgs_for_particles, or
!>       collision_turbulence).
!------------------------------------------------------------------------------!
 SUBROUTINE diffusion_e_ij( i, j, var, var_reference )

    USE arrays_3d,                                                             &
        ONLY:  dd2zu, ddzu, ddzw, drho_air, rho_air_zw

    USE control_parameters,                                                    &
        ONLY:  atmos_ocean_sign, use_single_reference_value

    USE grid_variables,                                                        &
        ONLY:  ddx2, ddy2

    USE bulk_cloud_model_mod,                                                  &
        ONLY:  collision_turbulence

    USE particle_attributes,                                                   &
        ONLY:  use_sgs_for_particles, wang_kernel

    IMPLICIT NONE

    INTEGER(iwp) ::  i              !< running index x direction
    INTEGER(iwp) ::  j              !< running index y direction
    INTEGER(iwp) ::  k              !< running index z direction
    INTEGER(iwp) ::  m              !< running index surface elements
    INTEGER(iwp) ::  surf_e         !< End index of surface elements at (j,i)-gridpoint
    INTEGER(iwp) ::  surf_s         !< Start index of surface elements at (j,i)-gridpoint

    REAL(wp)     ::  duv2_dz2       !< squared vertical gradient of wind vector
    REAL(wp)     ::  dvar_dz        !< vertical gradient of var
    REAL(wp)     ::  l              !< mixing length
    REAL(wp)     ::  var_reference  !< reference temperature

    REAL(wp), DIMENSION(nzb+1:nzt) ::  l_stable  !< mixing length according to stratification
    REAL(wp), DIMENSION(nzb+1:nzt) ::  rif       !< Richardson flux number

    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg), INTENT(IN) ::  var  !< temperature

!
!-- Calculate the mixing length and dissipation...
!-- ...in case of LES
    IF ( les_dynamic .OR. les_mw )  THEN

       DO  k = nzb+1, nzt
          dvar_dz = atmos_ocean_sign * ( var(k+1,j,i) - var(k-1,j,i) ) * dd2zu(k)
          IF ( dvar_dz > 0.0_wp ) THEN
             IF ( use_single_reference_value )  THEN
                l_stable(k) = 0.76_wp * SQRT( e(k,j,i) )                  &
                            / SQRT( g / var_reference * dvar_dz ) + 1E-5_wp
             ELSE
                l_stable(k) = 0.76_wp * SQRT( e(k,j,i) )               &
                            / SQRT( g / var(k,j,i) * dvar_dz ) + 1E-5_wp
             ENDIF
          ELSE
             l_stable(k) = l_grid(k)
          ENDIF
       ENDDO

       !DIR$ IVDEP
       DO  k = nzb+1, nzt
          l  = MIN( l_wall(k,j,i), l_stable(k) )

          diss(k,j,i) = ( 0.19_wp + 0.74_wp * l / l_wall(k,j,i) )              &
                        * e(k,j,i) * SQRT( e(k,j,i) ) / l                      &
                        * MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )
       ENDDO

!
!-- ...in case of RANS
    ELSEIF ( rans_tke_l )  THEN

!
!--    Calculate Richardson-flux number
       IF ( use_single_reference_value )  THEN
          DO  k = nzb+1, nzt
             dvar_dz = atmos_ocean_sign * ( var(k+1,j,i) - var(k-1,j,i) ) * dd2zu(k)

             duv2_dz2 =   ( ( u(k+1,j,i) - u(k-1,j,i) ) * dd2zu(k) )**2           &
                        + ( ( v(k+1,j,i) - v(k-1,j,i) ) * dd2zu(k) )**2           &
                        + 1E-30_wp

             rif(k) = MIN( MAX( g / var_reference * dvar_dz / duv2_dz2, -5.0_wp ),  1.0_wp )
          ENDDO
       ELSE
          DO  k = nzb+1, nzt
             dvar_dz = atmos_ocean_sign * ( var(k+1,j,i) - var(k-1,j,i) ) * dd2zu(k)

             duv2_dz2 =   ( ( u(k+1,j,i) - u(k-1,j,i) ) * dd2zu(k) )**2           &
                        + ( ( v(k+1,j,i) - v(k-1,j,i) ) * dd2zu(k) )**2           &
                        + 1E-30_wp

             rif(k) = MIN( MAX( g / var(k,j,i) * dvar_dz / duv2_dz2, -5.0_wp ), 1.0_wp )
          ENDDO
       ENDIF
!
!--    Calculate diabatic mixing length using Dyer-profile functions
       DO  k = nzb+1, nzt
          IF ( rif(k) >= 0.0_wp )  THEN
             l_stable(k) = MIN( l_black(k) / ( 1.0_wp + 5.0_wp * rif(k) ), l_wall(k,j,i) )
          ELSE
             l_stable(k) = l_wall(k,j,i) * SQRT( 1.0_wp - 16.0_wp * rif(k) )
          ENDIF

       ENDDO

       !DIR$ IVDEP
       DO  k = nzb+1, nzt
          diss(k,j,i) = c_0**3 * e(k,j,i) * SQRT( e(k,j,i) ) / l_stable(k)     &
                      * MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )
       ENDDO

!-- Note, in case of rans_tke_e, the dissipation is already calculated
!-- in prognostic_equations
    ENDIF

!
!-- Calculate the tendency term
    !DIR$ IVDEP
    DO  k = nzb+1, nzt

       tend(k,j,i) = tend(k,j,i) + (                                           &
                                      (                                        &
                      ( km(k,j,i)+km(k,j,i+1) ) * ( e(k,j,i+1)-e(k,j,i) )      &
                    - ( km(k,j,i)+km(k,j,i-1) ) * ( e(k,j,i)-e(k,j,i-1) )      &
                                      ) * ddx2                                 &
                                    + (                                        &
                      ( km(k,j,i)+km(k,j+1,i) ) * ( e(k,j+1,i)-e(k,j,i) )      &
                    - ( km(k,j,i)+km(k,j-1,i) ) * ( e(k,j,i)-e(k,j-1,i) )      &
                                      ) * ddy2                                 &
                                    + (                                        &
           ( km(k,j,i)+km(k+1,j,i) ) * ( e(k+1,j,i)-e(k,j,i) ) * ddzu(k+1)     &
                                                         * rho_air_zw(k)       &
         - ( km(k,j,i)+km(k-1,j,i) ) * ( e(k,j,i)-e(k-1,j,i) ) * ddzu(k)       &
                                                         * rho_air_zw(k-1)     &
                                      ) * ddzw(k) * drho_air(k)                &
                                   ) * dsig_e                                  &
                                     * MERGE( 1.0_wp, 0.0_wp,                  &
                                        BTEST( wall_flags_total_0(k,j,i), 0 ) )&
                                 - diss(k,j,i)

    ENDDO

!
!-- Set boundary conditions of dissipation if needed for calculating the sgs
!-- particle velocities.
!-- Neumann boundary condition for dissipation diss(nzb,:,:) = diss(nzb+1,:,:)
!-- For each surface type determine start and end index (in case of elevated
!-- topography several up/downward facing surfaces may exist.
!-- Note, bc cannot be set in tcm_boundary conditions as the dissipation
!-- in LES mode is only a diagnostic quantity.
    IF ( .NOT. rans_tke_e .AND.  ( use_sgs_for_particles  .OR.  wang_kernel    &
          .OR.  collision_turbulence ) )  THEN
       surf_s = bc_h(0)%start_index(j,i)
       surf_e = bc_h(0)%end_index(j,i)
       DO  m = surf_s, surf_e
          k             = bc_h(0)%k(m)
          diss(k-1,j,i) = diss(k,j,i)
       ENDDO
!
!--    Downward facing surfaces
       surf_s = bc_h(1)%start_index(j,i)
       surf_e = bc_h(1)%end_index(j,i)
       DO  m = surf_s, surf_e
          k             = bc_h(1)%k(m)
          diss(k+1,j,i) = diss(k,j,i)
       ENDDO
    ENDIF

 END SUBROUTINE diffusion_e_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Diffusion term for the TKE dissipation rate
!> Vector-optimized version
!------------------------------------------------------------------------------!
 SUBROUTINE diffusion_diss
    USE arrays_3d,                                                             &
        ONLY:  ddzu, ddzw, drho_air, rho_air_zw

    USE grid_variables,                                                        &
        ONLY:  ddx2, ddy2

    IMPLICIT NONE

    INTEGER(iwp) ::  i              !< running index x direction
    INTEGER(iwp) ::  j              !< running index y direction
    INTEGER(iwp) ::  k              !< running index z direction

    REAL(wp)     ::  flag           !< flag to mask topography

!
!-- Calculate the tendency terms
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt

!
!--          Predetermine flag to mask topography
             flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )

             tend(k,j,i) = tend(k,j,i) +                                       &
                         (       (                                             &
                 ( km(k,j,i)+km(k,j,i+1) ) * ( diss(k,j,i+1)-diss(k,j,i) )     &
               - ( km(k,j,i)+km(k,j,i-1) ) * ( diss(k,j,i)-diss(k,j,i-1) )     &
                                 ) * ddx2                                      &
                               + (                                             &
                 ( km(k,j,i)+km(k,j+1,i) ) * ( diss(k,j+1,i)-diss(k,j,i) )     &
               - ( km(k,j,i)+km(k,j-1,i) ) * ( diss(k,j,i)-diss(k,j-1,i) )     &
                                 ) * ddy2                                      &
                               + (                                             &
      ( km(k,j,i)+km(k+1,j,i) ) * ( diss(k+1,j,i)-diss(k,j,i) ) * ddzu(k+1)    &
                                                    * rho_air_zw(k)            &
    - ( km(k,j,i)+km(k-1,j,i) ) * ( diss(k,j,i)-diss(k-1,j,i) ) * ddzu(k)      &
                                                    * rho_air_zw(k-1)          &
                                 ) * ddzw(k) * drho_air(k)                     &
                         ) * flag * dsig_diss                                  &
                         - c_2 * diss(k,j,i)**2                                &
                               / ( e(k,j,i) + 1.0E-20_wp ) * flag

          ENDDO
       ENDDO
    ENDDO

 END SUBROUTINE diffusion_diss


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Diffusion term for the TKE dissipation rate
!> Cache-optimized version
!------------------------------------------------------------------------------!
 SUBROUTINE diffusion_diss_ij( i, j )

    USE arrays_3d,                                                             &
        ONLY:  ddzu, ddzw, drho_air, rho_air_zw

    USE grid_variables,                                                        &
        ONLY:  ddx2, ddy2

    IMPLICIT NONE

    INTEGER(iwp) ::  i         !< running index x direction
    INTEGER(iwp) ::  j         !< running index y direction
    INTEGER(iwp) ::  k         !< running index z direction

    REAL(wp)     ::  flag      !< flag to mask topography

!
!-- Calculate the mixing length (for dissipation)
    DO  k = nzb+1, nzt

!
!--    Predetermine flag to mask topography
       flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )

!
!--    Calculate the tendency term
       tend(k,j,i) =  tend(k,j,i) +                                            &
                   (            (                                              &
                ( km(k,j,i)+km(k,j,i+1) ) * ( diss(k,j,i+1)-diss(k,j,i) )      &
              - ( km(k,j,i)+km(k,j,i-1) ) * ( diss(k,j,i)-diss(k,j,i-1) )      &
                                ) * ddx2                                       &
                              + (                                              &
                ( km(k,j,i)+km(k,j+1,i) ) * ( diss(k,j+1,i)-diss(k,j,i) )      &
              - ( km(k,j,i)+km(k,j-1,i) ) * ( diss(k,j,i)-diss(k,j-1,i) )      &
                                ) * ddy2                                       &
                              + (                                              &
     ( km(k,j,i)+km(k+1,j,i) ) * ( diss(k+1,j,i)-diss(k,j,i) ) * ddzu(k+1)     &
                                                   * rho_air_zw(k)             &
   - ( km(k,j,i)+km(k-1,j,i) ) * ( diss(k,j,i)-diss(k-1,j,i) ) * ddzu(k)       &
                                                   * rho_air_zw(k-1)           &
                                ) * ddzw(k) * drho_air(k)                      &
                   ) * flag * dsig_diss                                        &
                   - c_2 * diss(k,j,i)**2 / ( e(k,j,i) + 1.0E-20_wp ) * flag

    ENDDO

 END SUBROUTINE diffusion_diss_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Computation of the turbulent diffusion coefficients for momentum and heat.
!> @bug unstable stratification is not properly considered for kh in rans mode.
!------------------------------------------------------------------------------!
 SUBROUTINE tcm_diffusivities( var, var_reference )

    USE control_parameters,                                                    &
        ONLY:  bc_radiation_l, bc_radiation_n, bc_radiation_r, bc_radiation_s, &
               e_min

    USE exchange_horiz_mod,                                                    &
        ONLY:  exchange_horiz

    USE surface_layer_fluxes_mod,                                              &
        ONLY:  phi_m

    INTEGER(iwp) ::  i          !< loop index
    INTEGER(iwp) ::  j          !< loop index
    INTEGER(iwp) ::  k          !< loop index
    INTEGER(iwp) ::  m          !< loop index
    INTEGER(iwp) ::  n          !< loop index

    REAL(wp) ::  var_reference  !< reference temperature

    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg), INTENT(IN) ::  var  !< temperature


!
!-- Introduce an optional minimum tke
    IF ( e_min > 0.0_wp )  THEN
       !$OMP PARALLEL DO PRIVATE(i,j,k)
       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
             DO  k = nzb+1, nzt
                e(k,j,i) = MAX( e(k,j,i), e_min ) *                            &
                        MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )
             ENDDO
          ENDDO
       ENDDO
    ENDIF

!
!-- Call default diffusivities routine. This is always used to calculate kh.
    CALL tcm_diffusivities_default( var, var_reference )
!
!-- Call dynamic subgrid model to calculate km.
    IF ( les_dynamic )  THEN
       CALL tcm_diffusivities_dynamic
    ENDIF

!
!-- In RANS mode, use MOST to calculate km and kh within the surface layer.
    IF ( rans_tke_e )  THEN
!
!--    Upward facing surfaces
!--    Default surfaces
       n = 0
       !$OMP PARALLEL DO PRIVATE(i,j,k,m)
       DO  m = 1, surf_def_h(0)%ns
          i = surf_def_h(0)%i(m)
          j = surf_def_h(0)%j(m)
          k = surf_def_h(0)%k(m)
          km(k,j,i) = kappa * surf_def_h(0)%us(m) * surf_def_h(0)%z_mo(m) /    &
                      phi_m( surf_def_h(0)%z_mo(m) / surf_def_h(0)%ol(m) )
          kh(k,j,i) = 1.35_wp * km(k,j,i)
       ENDDO
!
!--    Natural surfaces
       !$OMP PARALLEL DO PRIVATE(i,j,k,m)
       DO  m = 1, surf_lsm_h%ns
          i = surf_lsm_h%i(m)
          j = surf_lsm_h%j(m)
          k = surf_lsm_h%k(m)
          km(k,j,i) = kappa * surf_lsm_h%us(m) * surf_lsm_h%z_mo(m) /          &
                      phi_m( surf_lsm_h%z_mo(m) / surf_lsm_h%ol(m) )
          kh(k,j,i) = 1.35_wp * km(k,j,i)
       ENDDO
!
!--    Urban surfaces
       !$OMP PARALLEL DO PRIVATE(i,j,k,m)
       DO  m = 1, surf_usm_h%ns
          i = surf_usm_h%i(m)
          j = surf_usm_h%j(m)
          k = surf_usm_h%k(m)
          km(k,j,i) = kappa * surf_usm_h%us(m) * surf_usm_h%z_mo(m) /          &
                      phi_m( surf_usm_h%z_mo(m) / surf_usm_h%ol(m) )
          kh(k,j,i) = 1.35_wp * km(k,j,i)
       ENDDO

!
!--    North-, south-, west and eastward facing surfaces
!--    Do not consider stratification at these surfaces.
       DO  n = 0, 3
!
!--       Default surfaces
          !$OMP PARALLEL DO PRIVATE(i,j,k,m)
          DO  m = 1, surf_def_v(n)%ns
             i = surf_def_v(n)%i(m)
             j = surf_def_v(n)%j(m)
             k = surf_def_v(n)%k(m)
             km(k,j,i) = kappa * surf_def_v(n)%us(m) * surf_def_v(n)%z_mo(m)
             kh(k,j,i) = 1.35_wp * km(k,j,i)
          ENDDO
!
!--       Natural surfaces
          !$OMP PARALLEL DO PRIVATE(i,j,k,m)
          DO  m = 1, surf_lsm_v(n)%ns
             i = surf_lsm_v(n)%i(m)
             j = surf_lsm_v(n)%j(m)
             k = surf_lsm_v(n)%k(m)
             km(k,j,i) = kappa * surf_lsm_v(n)%us(m) * surf_lsm_v(n)%z_mo(m)
             kh(k,j,i) = 1.35_wp * km(k,j,i)
          ENDDO
!
!--       Urban surfaces
          !$OMP PARALLEL DO PRIVATE(i,j,k,m)
          DO  m = 1, surf_usm_v(n)%ns
             i = surf_usm_v(n)%i(m)
             j = surf_usm_v(n)%j(m)
             k = surf_usm_v(n)%k(m)
             km(k,j,i) = kappa * surf_usm_v(n)%us(m) * surf_usm_v(n)%z_mo(m)
             kh(k,j,i) = 1.35_wp * km(k,j,i)
          ENDDO
       ENDDO

       CALL exchange_horiz( km, nbgp )
       CALL exchange_horiz( kh, nbgp )

    ENDIF
!
!-- Set boundary values (Neumann conditions)
!-- Downward facing surfaces
    !$OMP PARALLEL DO PRIVATE(i,j,k)
    !$ACC PARALLEL LOOP PRIVATE(i,j,k) &
    !$ACC PRESENT(bc_h(1), kh, km)
    DO  m = 1, bc_h(1)%ns
       i = bc_h(1)%i(m)
       j = bc_h(1)%j(m)
       k = bc_h(1)%k(m)
       km(k+1,j,i) = km(k,j,i)
       kh(k+1,j,i) = kh(k,j,i)
    ENDDO
!
!-- Downward facing surfaces
    !$OMP PARALLEL DO PRIVATE(i,j,k)
    !$ACC PARALLEL LOOP PRIVATE(i,j,k) &
    !$ACC PRESENT(bc_h(0), kh, km)
    DO  m = 1, bc_h(0)%ns
       i = bc_h(0)%i(m)
       j = bc_h(0)%j(m)
       k = bc_h(0)%k(m)
       km(k-1,j,i) = km(k,j,i)
       kh(k-1,j,i) = kh(k,j,i)
    ENDDO
!
!-- Model top
    !$OMP PARALLEL DO
    !$ACC PARALLEL LOOP COLLAPSE(2) &
    !$ACC PRESENT(kh, km)
    DO  i = nxlg, nxrg
       DO  j = nysg, nyng
          km(nzt+1,j,i) = km(nzt,j,i)
          kh(nzt+1,j,i) = kh(nzt,j,i)
       ENDDO
    ENDDO

!
!-- Set Neumann boundary conditions at the outflow boundaries in case of
!-- non-cyclic lateral boundaries
    IF ( bc_radiation_l )  THEN
       km(:,:,nxl-1) = km(:,:,nxl)
       kh(:,:,nxl-1) = kh(:,:,nxl)
    ENDIF
    IF ( bc_radiation_r )  THEN
       km(:,:,nxr+1) = km(:,:,nxr)
       kh(:,:,nxr+1) = kh(:,:,nxr)
    ENDIF
    IF ( bc_radiation_s )  THEN
       km(:,nys-1,:) = km(:,nys,:)
       kh(:,nys-1,:) = kh(:,nys,:)
    ENDIF
    IF ( bc_radiation_n )  THEN
       km(:,nyn+1,:) = km(:,nyn,:)
       kh(:,nyn+1,:) = kh(:,nyn,:)
    ENDIF

 END SUBROUTINE tcm_diffusivities


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Computation of the turbulent diffusion coefficients for momentum and heat
!> according to Prandtl-Kolmogorov.
!------------------------------------------------------------------------------!
 SUBROUTINE tcm_diffusivities_default( var, var_reference )

    USE arrays_3d,                                                             &
        ONLY:  dd2zu

    USE control_parameters,                                                    &
        ONLY:  atmos_ocean_sign, use_single_reference_value

    USE statistics,                                                            &
        ONLY :  rmask, sums_l_l

    IMPLICIT NONE

    INTEGER(iwp) ::  i                   !< loop index
    INTEGER(iwp) ::  j                   !< loop index
    INTEGER(iwp) ::  k                   !< loop index
!$  INTEGER(iwp) ::  omp_get_thread_num  !< opemmp function to get thread number
    INTEGER(iwp) ::  sr                  !< statistic region
    INTEGER(iwp) ::  tn                  !< thread number

    REAL(wp)     ::  duv2_dz2            !< squared vertical gradient of wind vector
    REAL(wp)     ::  dvar_dz             !< vertical gradient of var
    REAL(wp)     ::  l                   !< mixing length (single height)
    REAL(wp)     ::  var_reference       !< reference temperature

    !DIR$ ATTRIBUTES ALIGN:64:: l_v, l_stable, rif
    REAL(wp), DIMENSION(nzb+1:nzt) ::  l_v       !< mixing length (all heights)
    REAL(wp), DIMENSION(nzb+1:nzt) ::  l_stable  !< mixing length according to stratification
    REAL(wp), DIMENSION(nzb+1:nzt) ::  rif       !< Richardson flux number

    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg), INTENT(IN) ::  var  !< temperature

!
!-- Default thread number in case of one thread
    tn = 0

!
!-- Initialization for calculation of the mixing length profile
    !$ACC KERNELS PRESENT(sums_l_l)
    sums_l_l = 0.0_wp
    !$ACC END KERNELS

!
!-- Compute the turbulent diffusion coefficient for momentum
    !$OMP PARALLEL PRIVATE (i,j,k,l,sr,tn)
!$  tn = omp_get_thread_num()

    IF ( les_dynamic .OR. les_mw )  THEN
       !$OMP DO
       !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(i, j) &
       !$ACC PRIVATE(dvar_dz, l, l_stable, l_v) &
       !$ACC PRESENT(wall_flags_total_0, var, dd2zu, e, l_wall, l_grid, rmask) &
       !$ACC PRESENT(kh, km, sums_l_l)
       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
             !$ACC LOOP PRIVATE(k)
             DO  k = nzb+1, nzt
!
!--             Determine the mixing length
!--             @note The following code cannot be transferred to a subroutine
!--             due to errors when using OpenACC directives. The execution
!--             crashes reliably if a subroutine is called at this point (the
!--             reasong for this behaviour is unknown, however).
                dvar_dz = atmos_ocean_sign * ( var(k+1,j,i) - var(k-1,j,i) ) * dd2zu(k)
                IF ( dvar_dz > 0.0_wp ) THEN
                   IF ( use_single_reference_value )  THEN
                      l_stable(k) = 0.76_wp * SQRT( e(k,j,i) )                  &
                                  / SQRT( g / var_reference * dvar_dz ) + 1E-5_wp
                   ELSE
                      l_stable(k) = 0.76_wp * SQRT( e(k,j,i) )               &
                                  / SQRT( g / var(k,j,i) * dvar_dz ) + 1E-5_wp
                   ENDIF
                ELSE
                   l_stable(k) = l_grid(k)
                ENDIF

             ENDDO

             !$ACC LOOP PRIVATE(k)
             !DIR$ IVDEP
             DO  k = nzb+1, nzt

                l_v(k) = MIN( l_wall(k,j,i), l_stable(k) )                      &
                       * MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )
                l = l_v(k)
!
!--             Compute diffusion coefficients for momentum and heat
!
!--             COVID-19 specific code: Viscosity of air added by Mikko.
                km(k,j,i) = c_0 * l * SQRT( e(k,j,i) ) + 1.4086E-5_wp
!
!--             COVID-19 specific code ends
                kh(k,j,i) = ( 1.0_wp + 2.0_wp * l / l_wall(k,j,i) ) * km(k,j,i)

             ENDDO
!
!--          Summation for averaged profile (cf. flow_statistics)
             !$ACC LOOP PRIVATE(sr, k)
             DO  sr = 0, statistic_regions
                DO  k = nzb+1, nzt
                   sums_l_l(k,sr,tn) = sums_l_l(k,sr,tn) + l_v(k) * rmask(j,i,sr)
                ENDDO
             ENDDO

          ENDDO
       ENDDO

    ELSEIF ( rans_tke_l )  THEN

       !$OMP DO
       !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(i, j) &
       !$ACC PRIVATE(dvar_dz, duv2_dz2, l_stable, l_v, rif) &
       !$ACC PRESENT(wall_flags_total_0, var, dd2zu, e, u, v, l_wall, l_black, rmask) &
       !$ACC PRESENT(kh, km, sums_l_l)
       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
!
!--          Calculate Richardson-flux number
             IF ( use_single_reference_value )  THEN
                !$ACC LOOP PRIVATE(k)
                !DIR$ IVDEP
                DO  k = nzb+1, nzt
                   dvar_dz = atmos_ocean_sign * ( var(k+1,j,i) - var(k-1,j,i) ) * dd2zu(k)

                   duv2_dz2 =   ( ( u(k+1,j,i) - u(k-1,j,i) ) * dd2zu(k) )**2  &
                              + ( ( v(k+1,j,i) - v(k-1,j,i) ) * dd2zu(k) )**2  &
                              + 1E-30_wp

                   rif(k) = MIN( MAX( g / var_reference * dvar_dz / duv2_dz2, -5.0_wp ),  1.0_wp )
                ENDDO
             ELSE
                !$ACC LOOP PRIVATE(k)
                !DIR$ IVDEP
                DO  k = nzb+1, nzt
                   dvar_dz = atmos_ocean_sign * ( var(k+1,j,i) - var(k-1,j,i) ) * dd2zu(k)

                   duv2_dz2 =   ( ( u(k+1,j,i) - u(k-1,j,i) ) * dd2zu(k) )**2  &
                              + ( ( v(k+1,j,i) - v(k-1,j,i) ) * dd2zu(k) )**2  &
                              + 1E-30_wp

                   rif(k) = MIN( MAX( g / var(k,j,i) * dvar_dz / duv2_dz2, -5.0_wp ),  1.0_wp )
                ENDDO
             ENDIF
!
!--          Calculate diabatic mixing length using Dyer-profile functions
!--          In case of unstable stratification, use mixing length of neutral case
             !$ACC LOOP PRIVATE(k)
             DO  k = nzb+1, nzt
                IF ( rif(k) >= 0.0_wp )  THEN
                   l_stable(k)  = MIN( l_black(k) / ( 1.0_wp + 5.0_wp * rif(k) ), l_wall(k,j,i) )
                ELSE
                   l_stable(k)  = l_wall(k,j,i)
                ENDIF

             ENDDO
!
!--          Compute diffusion coefficients for momentum and heat
             !$ACC LOOP PRIVATE(k)
             !DIR$ IVDEP
             DO  k = nzb+1, nzt
                l_v(k)    = l_stable(k) * MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )
                km(k,j,i) = c_0 * l_v(k) * SQRT( e(k,j,i) )
                kh(k,j,i) = km(k,j,i) / prandtl_number
             ENDDO
!
!--          Summation for averaged profile (cf. flow_statistics)
             !$ACC LOOP PRIVATE(sr, k)
             DO  sr = 0, statistic_regions
                DO  k = nzb+1, nzt
                   sums_l_l(k,sr,tn) = sums_l_l(k,sr,tn) + l_v(k) * rmask(j,i,sr)
                ENDDO
             ENDDO

          ENDDO
       ENDDO

    ELSEIF ( rans_tke_e )  THEN

       !$OMP DO
       !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(i, j) &
       !$ACC PRIVATE(l_v) &
       !$ACC PRESENT(wall_flags_total_0, e, diss, rmask) &
       !$ACC PRESENT(kh, km, sums_l_l)
       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
!
!--          Compute diffusion coefficients for momentum and heat
             !$ACC LOOP PRIVATE(k)
             !DIR$ IVDEP
             DO  k = nzb+1, nzt

                l_v(k) = c_0**3 * e(k,j,i) * SQRT(e(k,j,i)) / ( diss(k,j,i) + 1.0E-30_wp ) &
                       * MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )

                km(k,j,i) = c_0 * SQRT( e(k,j,i) ) * l_v(k)
                kh(k,j,i) = km(k,j,i) / prandtl_number

             ENDDO
!
!--          Summation for averaged profile of mixing length (cf. flow_statistics)
             !$ACC LOOP PRIVATE(sr, k)
             DO  sr = 0, statistic_regions
                DO  k = nzb+1, nzt
                   sums_l_l(k,sr,tn) = sums_l_l(k,sr,tn) + l_v(k) * rmask(j,i,sr)
                ENDDO
             ENDDO

          ENDDO
       ENDDO

    ENDIF

    !$ACC KERNELS PRESENT(sums_l_l)
    sums_l_l(nzt+1,:,tn) = sums_l_l(nzt,:,tn)   ! quasi boundary-condition for
                                                ! data output
    !$ACC END KERNELS
!$OMP END PARALLEL

 END SUBROUTINE tcm_diffusivities_default


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculates the eddy viscosity dynamically using the linear dynamic model
!> according to
!> Heinz, Stefan. "Realizability of dynamic subgrid-scale stress models via
!> stochastic analysis."
!> Monte Carlo Methods and Applications 14.4 (2008): 311-329.
!>
!> Furthermore dynamic bounds are used to limit the absolute value of c* as
!> described in
!> Mokhtarpoor, Reza, and Stefan Heinz. "Dynamic large eddy simulation:
!> Stability via realizability." Physics of Fluids 29.10 (2017): 105104.
!>
!> @author Hauke Wurps
!> @author Björn Maronga
!------------------------------------------------------------------------------!
 SUBROUTINE tcm_diffusivities_dynamic

    USE arrays_3d,                                                             &
        ONLY:  ddzw, dzw, dd2zu, w, ug, vg

    USE grid_variables,                                                        &
        ONLY : ddx, ddy, dx, dy

    IMPLICIT NONE

    INTEGER(iwp) ::  i           !< running index x-direction
    INTEGER(iwp) ::  j           !< running index y-direction
    INTEGER(iwp) ::  k           !< running index z-direction
    INTEGER(iwp) ::  l           !< running index
    INTEGER(iwp) ::  m           !< running index

    REAL(wp)     ::  dudx        !< Gradient of u-component in x-direction
    REAL(wp)     ::  dudy        !< Gradient of u-component in y-direction
    REAL(wp)     ::  dudz        !< Gradient of u-component in z-direction
    REAL(wp)     ::  dvdx        !< Gradient of v-component in x-direction
    REAL(wp)     ::  dvdy        !< Gradient of v-component in y-direction
    REAL(wp)     ::  dvdz        !< Gradient of v-component in z-direction
    REAL(wp)     ::  dwdx        !< Gradient of w-component in x-direction
    REAL(wp)     ::  dwdy        !< Gradient of w-component in y-direction
    REAL(wp)     ::  dwdz        !< Gradient of w-component in z-direction

    REAL(wp)     ::  flag        !< topography flag

    REAL(wp)     ::  uc(-1:1,-1:1)  !< u on grid center
    REAL(wp)     ::  vc(-1:1,-1:1)  !< v on grid center
    REAL(wp)     ::  wc(-1:1,-1:1)  !< w on grid center

    REAL(wp)     ::  ut(nzb:nzt+1,nysg:nyng,nxlg:nxrg)  !< test filtered u
    REAL(wp)     ::  vt(nzb:nzt+1,nysg:nyng,nxlg:nxrg)  !< test filtered v
    REAL(wp)     ::  wt(nzb:nzt+1,nysg:nyng,nxlg:nxrg)  !< test filtered w

    REAL(wp)     ::  uct         !< test filtered u on grid center
    REAL(wp)     ::  vct         !< test filtered v on grid center
    REAL(wp)     ::  wct         !< test filtered w on grid center
    REAL(wp)     ::  u2t         !< test filtered u**2 on grid center
    REAL(wp)     ::  v2t         !< test filtered v**2 on grid center
    REAL(wp)     ::  w2t         !< test filtered w**2 on grid center
    REAL(wp)     ::  uvt         !< test filtered u*v on grid center
    REAL(wp)     ::  uwt         !< test filtered u*w on grid center
    REAL(wp)     ::  vwt         !< test filtered v*w on grid center

    REAL(wp)     ::  sd11        !< deviatoric shear tensor
    REAL(wp)     ::  sd22        !< deviatoric shear tensor
    REAL(wp)     ::  sd33        !<f deviatoric shear tensor
    REAL(wp)     ::  sd12        !< deviatoric shear tensor
    REAL(wp)     ::  sd13        !< deviatoric shear tensor
    REAL(wp)     ::  sd23        !< deviatoric shear tensor

    REAL(wp)     ::  sd2         !< sum: sd_ij*sd_ij

    REAL(wp)     ::  sdt11       !< filtered deviatoric shear tensor
    REAL(wp)     ::  sdt22       !< filtered deviatoric shear tensor
    REAL(wp)     ::  sdt33       !< filtered deviatoric shear tensor
    REAL(wp)     ::  sdt12       !< filtered deviatoric shear tensor
    REAL(wp)     ::  sdt13       !< filtered deviatoric shear tensor
    REAL(wp)     ::  sdt23       !< filtered deviatoric shear tensor

    REAL(wp)     ::  sdt2        !< sum: sdt_ij*sdt_ij

    REAL(wp)     ::  ld11        !< deviatoric stress tensor
    REAL(wp)     ::  ld22        !< deviatoric stress tensor
    REAL(wp)     ::  ld33        !< deviatoric stress tensor
    REAL(wp)     ::  ld12        !< deviatoric stress tensor
    REAL(wp)     ::  ld13        !< deviatoric stress tensor
    REAL(wp)     ::  ld23        !< deviatoric stress tensor

    REAL(wp)     ::  lnn         !< sum ld_nn
    REAL(wp)     ::  ldsd        !< sum: ld_ij*sd_ij

    REAL(wp)     ::  delta       !< grid size
    REAL(wp)     ::  cst         !< c*
    REAL(wp)     ::  cstnust_t   !< product c*nu*
    REAL(wp)     ::  cst_max     !< bounds of c*

    REAL(wp), PARAMETER :: fac_cmax = 23.0_wp/(24.0_wp*sqrt(3.0_wp))  !< constant

!
!-- velocities on grid centers:
    CALL tcm_box_filter_2d_array( u, ut )
    CALL tcm_box_filter_2d_array( v, vt )
    CALL tcm_box_filter_2d_array( w, wt )

    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt

             flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )

!
!--          Compute the deviatoric shear tensor s_ij on grid centers:
!--          s_ij =  0.5 * ( du_i/dx_j + du_j/dx_i )
             dudx  =           ( u(k,j,i+1) - u(k,j,i)     ) * ddx
             dudy  = 0.25_wp * ( u(k,j+1,i) + u(k,j+1,i+1) - &
                                 u(k,j-1,i) - u(k,j-1,i+1) ) * ddy
             dudz  = 0.5_wp  * ( u(k+1,j,i) + u(k+1,j,i+1) - &
                                 u(k-1,j,i) - u(k-1,j,i+1) ) * dd2zu(k)

             dvdx  = 0.25_wp * ( v(k,j,i+1) + v(k,j+1,i+1) - &
                                 v(k,j,i-1) - v(k,j+1,i-1) ) * ddx
             dvdy  =           ( v(k,j+1,i) - v(k,j,i)     ) * ddy
             dvdz  = 0.5_wp  * ( v(k+1,j,i) + v(k+1,j+1,i) - &
                                 v(k-1,j,i) - v(k-1,j+1,i) ) * dd2zu(k)

             dwdx  = 0.25_wp * ( w(k,j,i+1) + w(k-1,j,i+1) - &
                                 w(k,j,i-1) - w(k-1,j,i-1) ) * ddx
             dwdy  = 0.25_wp * ( w(k,j+1,i) + w(k-1,j+1,i) - &
                                 w(k,j-1,i) - w(k-1,j-1,i) ) * ddy
             dwdz  =           ( w(k,j,i)   - w(k-1,j,i)   ) * ddzw(k)

             sd11 = dudx
             sd22 = dvdy
             sd33 = dwdz
             sd12 = 0.5_wp * ( dudy + dvdx )
             sd13 = 0.5_wp * ( dudz + dwdx )
             sd23 = 0.5_wp * ( dvdz + dwdy )
!
!--          sum: sd_ij*sd_ij
             sd2 = sd11**2 + sd22**2 + sd33**2                     &
                   + 2.0_wp * ( sd12**2 + sd13**2 + sd23**2 )
!
!--          The filtered velocities are needed to calculate the filtered shear
!--          tensor: sdt_ij = 0.5 * ( dut_i/dx_j + dut_j/dx_i )
             dudx  =           ( ut(k,j,i+1) - ut(k,j,i)     ) * ddx
             dudy  = 0.25_wp * ( ut(k,j+1,i) + ut(k,j+1,i+1) - &
                                 ut(k,j-1,i) - ut(k,j-1,i+1) ) * ddy
             dudz  = 0.5_wp  * ( ut(k+1,j,i) + ut(k+1,j,i+1) - &
                                 ut(k-1,j,i) - ut(k-1,j,i+1) ) * dd2zu(k)

             dvdx  = 0.25_wp * ( vt(k,j,i+1) + vt(k,j+1,i+1) - &
                                 vt(k,j,i-1) - vt(k,j+1,i-1) ) * ddx
             dvdy  =           ( vt(k,j+1,i) - vt(k,j,i)     ) * ddy
             dvdz  = 0.5_wp  * ( vt(k+1,j,i) + vt(k+1,j+1,i) - &
                                 vt(k-1,j,i) - vt(k-1,j+1,i) ) * dd2zu(k)

             dwdx  = 0.25_wp * ( wt(k,j,i+1) + wt(k-1,j,i+1) - &
                                 wt(k,j,i-1) - wt(k-1,j,i-1) ) * ddx
             dwdy  = 0.25_wp * ( wt(k,j+1,i) + wt(k-1,j+1,i) - &
                                 wt(k,j-1,i) - wt(k-1,j-1,i) ) * ddy
             dwdz  =           ( wt(k,j,i)   - wt(k-1,j,i)   ) * ddzw(k)

             sdt11 = dudx
             sdt22 = dvdy
             sdt33 = dwdz
             sdt12 = 0.5_wp * ( dudy + dvdx )
             sdt13 = 0.5_wp * ( dudz + dwdx )
             sdt23 = 0.5_wp * ( dvdz + dwdy )
!
!--          sum: sd_ij*sd_ij
             sdt2 = sdt11**2 + sdt22**2 + sdt33**2         &
                   + 2.0_wp * ( sdt12**2 + sdt13**2 + sdt23**2 )
!
!--          Need filtered velocities and filtered squared velocities on grid
!--          centers. Substraction of geostrophic velocity helps to avoid
!--          numerical errors in the expression <u**2> - <u>*<u>, which can be
!--          very small (<...> means filtered).
             DO  l = -1, 1
                DO  m = -1, 1
                   uc(l,m) = 0.5_wp * ( u(k,j+l,i+m) + u(k,j+l,i+m+1) ) - ug(k)
                   vc(l,m) = 0.5_wp * ( v(k,j+l,i+m) + v(k,j+l+1,i+m) ) - vg(k)
                   wc(l,m) = 0.5_wp * ( w(k-1,j+l,i+m) + w(k,j+l,i+m) )
                ENDDO
             ENDDO

             CALL tcm_box_filter_2d_single( uc, uct )
             CALL tcm_box_filter_2d_single( vc, vct )
             CALL tcm_box_filter_2d_single( wc, wct )
             CALL tcm_box_filter_2d_single( uc**2, u2t )
             CALL tcm_box_filter_2d_single( vc**2, v2t )
             CALL tcm_box_filter_2d_single( wc**2, w2t )
             CALL tcm_box_filter_2d_single( uc*vc, uvt )
             CALL tcm_box_filter_2d_single( uc*wc, uwt )
             CALL tcm_box_filter_2d_single( vc*wc, vwt )

             ld11 = u2t - uct*uct
             ld22 = v2t - vct*vct
             ld33 = w2t - wct*wct
             ld12 = uvt - uct*vct
             ld13 = uwt - uct*wct
             ld23 = vwt - vct*wct

             lnn = ld11 + ld22 + ld33
!
!--          Substract trace to get deviatoric resolved stress
             ld11 = ld11 - lnn / 3.0_wp
             ld22 = ld22 - lnn / 3.0_wp
             ld33 = ld33 - lnn / 3.0_wp

             ldsd = ld11 * sdt11 + ld22 * sdt22 + ld33 * sdt33 + &
                    2.0_wp * ( ld12 * sdt12 + ld13 * sdt13 + ld23 * sdt23 )
!
!--          c* nu*^T is SGS viscosity on test filter level:
             cstnust_t = -ldsd / ( sdt2 + 1.0E-20_wp )
!
!--          The model was only tested for an isotropic grid. The following
!--          expression was a recommendation of Stefan Heinz.
             delta = MAX( dx, dy, dzw(k) )

             IF ( lnn <= 0.0_wp ) THEN
                cst = 0.0_wp
             ELSE
                cst = cstnust_t /                                              &
                   ( 4.0_wp * delta * SQRT( lnn / 2.0_wp ) + 1.0E-20_wp )
             ENDIF

!
!--          Calculate border according to Mokhtarpoor and Heinz (2017)
             cst_max = fac_cmax * SQRT( e(k,j,i) ) /                           &
                       ( delta * SQRT( 2.0_wp * sd2 ) + 1.0E-20_wp )

             IF ( ABS( cst ) > cst_max )  THEN
                cst = cst_max * cst / ABS( cst )
             ENDIF

             km(k,j,i) = cst * delta * SQRT( e(k,j,i) ) * flag

          ENDDO
       ENDDO
    ENDDO

 END SUBROUTINE tcm_diffusivities_dynamic


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This subroutine acts as a box filter with filter width 2 * dx.
!> Output is only one point.
!------------------------------------------------------------------------------!
 SUBROUTINE tcm_box_filter_2d_single( var, var_fil )

    IMPLICIT NONE

    REAL(wp)     ::  var(-1:1,-1:1)      !< variable to be filtered
    REAL(wp)     ::  var_fil             !< filtered variable
!
!-- It is assumed that a box with a side length of 2 * dx and centered at the
!-- variable*s position contains one half of the four closest neigbours and one
!-- forth of the diagonally closest neighbours.
    var_fil = 0.25_wp * ( var(0,0) +                                           &
                      0.5_wp  * ( var(0,1)  + var(1,0)   +                     &
                                  var(0,-1) + var(-1,0)  ) +                   &
                      0.25_wp * ( var(1,1)  + var(1,-1)  +                     &
                                  var(-1,1) + var(-1,-1) ) )

 END SUBROUTINE tcm_box_filter_2d_single

!------------------------------------------------------------------------------!
! Description:
! ------------
!> This subroutine acts as a box filter with filter width 2 * dx.
!> The filtered variable var_fil is on the same grid as var.
!------------------------------------------------------------------------------!
 SUBROUTINE tcm_box_filter_2d_array( var, var_fil )

    IMPLICIT NONE

    INTEGER(iwp) ::  i    !< running index x-direction
    INTEGER(iwp) ::  j    !< running index y-direction
    INTEGER(iwp) ::  k    !< running index z-direction

    REAL(wp)     ::  var(nzb:nzt+1,nysg:nyng,nxlg:nxrg)      !< variable to be filtered
    REAL(wp)     ::  var_fil(nzb:nzt+1,nysg:nyng,nxlg:nxrg)  !< filtered variable
!
!-- It is assumed that a box with a side length of 2 * dx and centered at the
!-- variable's position contains one half of the four closest neigbours and one
!-- forth of the diagonally closest neighbours.
    DO  i = nxlg+1, nxrg-1
       DO  j = nysg+1, nyng-1
          DO  k = nzb, nzt+1
             var_fil(k,j,i) = 0.25_wp * ( var(k,j,i) +                         &
                              0.5_wp  * ( var(k,j,i+1)   + var(k,j+1,i)   +    &
                                          var(k,j,i-1)   + var(k,j-1,i)     ) +&
                              0.25_wp * ( var(k,j+1,i+1) + var(k,j+1,i-1) +    &
                                          var(k,j-1,i+1) + var(k,j-1,i-1)   )  )
          END DO
       END DO
    END DO

 END SUBROUTINE tcm_box_filter_2d_array


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Swapping of timelevels.
!------------------------------------------------------------------------------!
 SUBROUTINE tcm_swap_timelevel ( mod_count )

    IMPLICIT NONE


    INTEGER, INTENT(IN) ::  mod_count  !< flag defining where pointers point to


    SELECT CASE ( mod_count )

       CASE ( 0 )

          IF ( .NOT. constant_diffusion )  THEN
             e => e_1;    e_p => e_2
          ENDIF

          IF ( rans_tke_e )  THEN
             diss => diss_1;    diss_p => diss_2
          ENDIF

       CASE ( 1 )

          IF ( .NOT. constant_diffusion )  THEN
             e => e_2;    e_p => e_1
          ENDIF

          IF ( rans_tke_e )  THEN
             diss => diss_2;    diss_p => diss_1
          ENDIF

    END SELECT

 END SUBROUTINE tcm_swap_timelevel


 END MODULE turbulence_closure_mod
