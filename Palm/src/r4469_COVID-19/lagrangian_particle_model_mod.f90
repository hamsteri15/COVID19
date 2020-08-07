!> @file lagrangian_particle_model_mod.f90
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
! $Id: lagrangian_particle_model_mod.f90 4457 2020-03-11 14:20:43Z raasch $
! use statement for exchange horiz added
! 
! 4444 2020-03-05 15:59:50Z raasch
! bugfix: cpp-directives for serial mode added
! 
! 4430 2020-02-27 18:02:20Z suehring
! - Bugfix in logarithmic interpolation of near-ground particle speed (density
!   was not considered).
! - Revise CFL-check when SGS particle speeds are considered.
! - In nested case with SGS particle speeds in the child domain, do not give 
!   warning that particles are on domain boundaries. At the end of the particle
!   time integration these will be transferred to the parent domain anyhow.
! 
! 4360 2020-01-07 11:25:50Z suehring
! Introduction of wall_flags_total_0, which currently sets bits based on static
! topography information used in wall_flags_static_0
! 
! 4336 2019-12-13 10:12:05Z raasch
! bugfix: wrong header output of particle group features (density ratio) in case
! of restarts corrected
! 
! 4329 2019-12-10 15:46:36Z motisi
! Renamed wall_flags_0 to wall_flags_static_0
! 
! 4282 2019-10-29 16:18:46Z schwenkel
! Bugfix of particle timeseries in case of more than one particle group
! 
! 4277 2019-10-28 16:53:23Z schwenkel
! Bugfix: Added first_call_lpm in use statement
! 
! 4276 2019-10-28 16:03:29Z schwenkel
! Modularize lpm: Move conditions in time intergration to module
! 
! 4275 2019-10-28 15:34:55Z schwenkel
! Change call of simple predictor corrector method, i.e. two divergence free
! velocitiy fields are now used.
!
! 4232 2019-09-20 09:34:22Z knoop
! Removed INCLUDE "mpif.h", as it is not needed because of USE pegrid
! 
! 4195 2019-08-28 13:44:27Z schwenkel
! Bugfix for simple_corrector interpolation method in case of ocean runs and
! output particle advection interpolation method into header
! 
! 4182 2019-08-22 15:20:23Z scharf
! Corrected "Former revisions" section
! 
! 4168 2019-08-16 13:50:17Z suehring
! Replace function get_topography_top_index by topo_top_ind
! 
! 4145 2019-08-06 09:55:22Z schwenkel
! Some reformatting
! 
! 4144 2019-08-06 09:11:47Z raasch
! relational operators .EQ., .NE., etc. replaced by ==, /=, etc.
! 
! 4143 2019-08-05 15:14:53Z schwenkel
! Rename variable and change select case to if statement
! 
! 4122 2019-07-26 13:11:56Z schwenkel
! Implement reset method as bottom boundary condition
! 
! 4121 2019-07-26 10:01:22Z schwenkel
! Implementation of an simple method for interpolating the velocities to
! particle position
! 
! 4114 2019-07-23 14:09:27Z schwenkel
! Bugfix: Added working precision for if statement
! 
! 4054 2019-06-27 07:42:18Z raasch
! bugfix for calculating the minimum particle time step
! 
! 4044 2019-06-19 12:28:27Z schwenkel
! Bugfix in case of grid strecting: corrected calculation of k-Index
!
! 4043 2019-06-18 16:59:00Z schwenkel
! Remove min_nr_particle, Add lpm_droplet_interactions_ptq into module
! 
! 4028 2019-06-13 12:21:37Z schwenkel
! Further modularization of particle code components
! 
! 4020 2019-06-06 14:57:48Z schwenkel
! Removing submodules 
! 
! 4018 2019-06-06 13:41:50Z eckhard
! Bugfix for former revision
! 
! 4017 2019-06-06 12:16:46Z schwenkel
! Modularization of all lagrangian particle model code components
! 
! 3655 2019-01-07 16:51:22Z knoop 
! bugfix to guarantee correct particle releases in case that the release
! interval is smaller than the model timestep
!
! Revision 1.1  1999/11/25 16:16:06  raasch
! Initial revision
!
!
! Description:
! ------------
!> The embedded LPM allows for studying transport and dispersion processes within
!> turbulent flows. This model including passive particles that do not show any
!> feedback on the turbulent flow. Further also particles with inertia and
!> cloud droplets ca be simulated explicitly.
!>
!> @todo test lcm
!>       implement simple interpolation method for subgrid scale velocites
!> @note <Enter notes on the module>
!> @bug  <Enter bug on the module>
!------------------------------------------------------------------------------!
 MODULE lagrangian_particle_model_mod

    USE, INTRINSIC ::  ISO_C_BINDING

    USE arrays_3d,                                                             &
        ONLY:  de_dx, de_dy, de_dz,                                            &
               d_exner,                                                        &
               drho_air_zw,                                                    &
               dzw, zu, zw,  ql_c, ql_v, ql_vp, hyp,                           &
               pt, q, exner, ql, diss, e, u, v, w, km, ql_1, ql_2, pt_p, q_p
 
    USE averaging,                                                             &
        ONLY:  ql_c_av, pr_av, pc_av, ql_vp_av, ql_v_av

    USE basic_constants_and_equations_mod,                                     &
        ONLY: molecular_weight_of_solute, molecular_weight_of_water, magnus,   &
              pi, rd_d_rv, rho_l, r_v, rho_s, vanthoff, l_v, kappa, g, lv_d_cp

    USE control_parameters,                                                    &
        ONLY:  bc_dirichlet_l, bc_dirichlet_n, bc_dirichlet_r, bc_dirichlet_s, &
               child_domain,                                                   &
               cloud_droplets, constant_flux_layer, current_timestep_number,   &
               dt_3d, dt_3d_reached, first_call_lpm, humidity,                 &
               dt_3d_reached_l, dt_dopts, dz, initializing_actions,            &
               intermediate_timestep_count, intermediate_timestep_count_max,   &
               message_string, molecular_viscosity, ocean_mode,                &
               particle_maximum_age, iran,                                     & 
               simulated_time, topography, dopts_time_count,                   &
               time_since_reference_point, rho_surface, u_gtrans, v_gtrans,    &
               dz_stretch_level, dz_stretch_level_start

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point, log_point_s

    USE indices,                                                               &
        ONLY:  nx, nxl, nxlg, nxrg, nxr, ny, nyn, nys, nyng, nysg, nz, nzb,    &
               nzb_max, nzt,nbgp, ngp_2dh_outer,                               &
               topo_top_ind,                                                   &
               wall_flags_total_0

    USE kinds

    USE pegrid

    USE particle_attributes

#if defined( __parallel )
    USE pmc_particle_interface,                                                &
        ONLY: pmcp_c_get_particle_from_parent, pmcp_p_fill_particle_win,       &
              pmcp_c_send_particle_to_parent, pmcp_p_empty_particle_win,       &
              pmcp_p_delete_particles_in_fine_grid_area, pmcp_g_init,          &
              pmcp_g_print_number_of_particles
#endif

    USE pmc_interface,                                                         &
        ONLY: nested_run

    USE grid_variables,                                                        &
        ONLY:  ddx, dx, ddy, dy

    USE netcdf_interface,                                                      &
        ONLY:  netcdf_data_format, netcdf_deflate, dopts_num, id_set_pts,      &
               id_var_dopts, id_var_time_pts, nc_stat,                         &
               netcdf_handle_error

    USE random_function_mod,                                                   &
        ONLY:  random_function

    USE statistics,                                                            &
        ONLY:  hom

    USE surface_mod,                                                           &
        ONLY:  bc_h,                                                           &
               surf_def_h,                                                     &
               surf_lsm_h,                                                     &
               surf_usm_h

#if defined( __parallel )  &&  !defined( __mpifh )
    USE MPI
#endif

#if defined( __netcdf )
    USE NETCDF
#endif

    IMPLICIT NONE

    CHARACTER(LEN=15) ::  aero_species = 'nacl'                   !< aerosol species
    CHARACTER(LEN=15) ::  aero_type    = 'maritime'               !< aerosol type
    CHARACTER(LEN=15) ::  bc_par_lr    = 'cyclic'                 !< left/right boundary condition
    CHARACTER(LEN=15) ::  bc_par_ns    = 'cyclic'                 !< north/south boundary condition
    CHARACTER(LEN=15) ::  bc_par_b     = 'reflect'                !< bottom boundary condition
    CHARACTER(LEN=15) ::  bc_par_t     = 'absorb'                 !< top boundary condition
    CHARACTER(LEN=15) ::  collision_kernel   = 'none'             !< collision kernel

    CHARACTER(LEN=5)  ::  splitting_function = 'gamma'            !< function for calculation critical weighting factor
    CHARACTER(LEN=5)  ::  splitting_mode     = 'const'            !< splitting mode

    CHARACTER(LEN=25) ::  particle_advection_interpolation = 'trilinear' !< interpolation method for calculatin the particle

    INTEGER(iwp) ::  deleted_particles = 0                        !< number of deleted particles per time step    
    INTEGER(iwp) ::  i_splitting_mode                             !< dummy for splitting mode
    INTEGER(iwp) ::  iran_part = -1234567                         !< number for random generator    
    INTEGER(iwp) ::  max_number_particles_per_gridbox = 100       !< namelist parameter (see documentation)
    INTEGER(iwp) ::  isf                                          !< dummy for splitting function
    INTEGER(iwp) ::  number_particles_per_gridbox = -1            !< namelist parameter (see documentation)
    INTEGER(iwp) ::  number_of_sublayers = 20                     !< number of sublayers for particle velocities betwenn surface and first grid level
    INTEGER(iwp) ::  offset_ocean_nzt = 0                         !< in case of oceans runs, the vertical index calculations need an offset
    INTEGER(iwp) ::  offset_ocean_nzt_m1 = 0                      !< in case of oceans runs, the vertical index calculations need an offset
    INTEGER(iwp) ::  particles_per_point = 1                      !< namelist parameter (see documentation)
    INTEGER(iwp) ::  radius_classes = 20                          !< namelist parameter (see documentation)
    INTEGER(iwp) ::  splitting_factor = 2                         !< namelist parameter (see documentation)
    INTEGER(iwp) ::  splitting_factor_max = 5                     !< namelist parameter (see documentation)
    INTEGER(iwp) ::  step_dealloc = 100                           !< namelist parameter (see documentation)
    INTEGER(iwp) ::  total_number_of_particles                    !< total number of particles in the whole model domain
    INTEGER(iwp) ::  trlp_count_sum                               !< parameter for particle exchange of PEs
    INTEGER(iwp) ::  trlp_count_recv_sum                          !< parameter for particle exchange of PEs
    INTEGER(iwp) ::  trrp_count_sum                               !< parameter for particle exchange of PEs
    INTEGER(iwp) ::  trrp_count_recv_sum                          !< parameter for particle exchange of PEs
    INTEGER(iwp) ::  trsp_count_sum                               !< parameter for particle exchange of PEs
    INTEGER(iwp) ::  trsp_count_recv_sum                          !< parameter for particle exchange of PEs
    INTEGER(iwp) ::  trnp_count_sum                               !< parameter for particle exchange of PEs
    INTEGER(iwp) ::  trnp_count_recv_sum                          !< parameter for particle exchange of PEs

    LOGICAL ::  lagrangian_particle_model = .FALSE.       !< namelist parameter (see documentation) 
    LOGICAL ::  curvature_solution_effects = .FALSE.      !< namelist parameter (see documentation)
    LOGICAL ::  deallocate_memory = .TRUE.                !< namelist parameter (see documentation)
    LOGICAL ::  hall_kernel = .FALSE.                     !< flag for collision kernel
    LOGICAL ::  merging = .FALSE.                         !< namelist parameter (see documentation)
    LOGICAL ::  random_start_position = .FALSE.           !< namelist parameter (see documentation)
    LOGICAL ::  read_particles_from_restartfile = .TRUE.  !< namelist parameter (see documentation)
    LOGICAL ::  seed_follows_topography = .FALSE.         !< namelist parameter (see documentation)
    LOGICAL ::  splitting = .FALSE.                       !< namelist parameter (see documentation)
    LOGICAL ::  use_kernel_tables = .FALSE.               !< parameter, which turns on the use of precalculated collision kernels
    LOGICAL ::  write_particle_statistics = .FALSE.       !< namelist parameter (see documentation)
    LOGICAL ::  interpolation_simple_predictor = .FALSE.  !< flag for simple particle advection interpolation with predictor step
    LOGICAL ::  interpolation_simple_corrector = .FALSE.  !< flag for simple particle advection interpolation with corrector step
    LOGICAL ::  interpolation_trilinear = .FALSE.         !< flag for trilinear particle advection interpolation

    LOGICAL, DIMENSION(max_number_of_particle_groups) ::   vertical_particle_advection = .TRUE. !< Switch for vertical particle transport

    REAL(wp) ::  aero_weight = 1.0_wp                      !< namelist parameter (see documentation)
    REAL(wp) ::  dt_min_part = 0.0002_wp                   !< minimum particle time step when SGS velocities are used (s)
    REAL(wp) ::  dt_prel = 9999999.9_wp                    !< namelist parameter (see documentation)
    REAL(wp) ::  dt_write_particle_data = 9999999.9_wp     !< namelist parameter (see documentation)
    REAL(wp) ::  end_time_prel = 9999999.9_wp              !< namelist parameter (see documentation)
    REAL(wp) ::  initial_weighting_factor = 1.0_wp         !< namelist parameter (see documentation)
    REAL(wp) ::  last_particle_release_time = 0.0_wp       !< last time of particle release
    REAL(wp) ::  log_sigma(3) = 1.0_wp                     !< namelist parameter (see documentation)
    REAL(wp) ::  na(3) = 0.0_wp                            !< namelist parameter (see documentation)
    REAL(wp) ::  number_concentration = -1.0_wp            !< namelist parameter (see documentation)
    REAL(wp) ::  radius_merge = 1.0E-7_wp                  !< namelist parameter (see documentation)
    REAL(wp) ::  radius_split = 40.0E-6_wp                 !< namelist parameter (see documentation)
    REAL(wp) ::  rm(3) = 1.0E-6_wp                         !< namelist parameter (see documentation)
    REAL(wp) ::  sgs_wf_part                               !< parameter for sgs
    REAL(wp) ::  time_write_particle_data = 0.0_wp         !< write particle data at current time on file
    REAL(wp) ::  weight_factor_merge = -1.0_wp             !< namelist parameter (see documentation)
    REAL(wp) ::  weight_factor_split = -1.0_wp             !< namelist parameter (see documentation)
    REAL(wp) ::  z0_av_global                              !< horizontal mean value of z0

    REAL(wp) ::  rclass_lbound !<
    REAL(wp) ::  rclass_ubound !<

    REAL(wp), PARAMETER ::  c_0 = 3.0_wp         !< parameter for lagrangian timescale

    REAL(wp), DIMENSION(max_number_of_particle_groups) ::  density_ratio = 9999999.9_wp  !< namelist parameter (see documentation)
    REAL(wp), DIMENSION(max_number_of_particle_groups) ::  pdx = 9999999.9_wp            !< namelist parameter (see documentation)
    REAL(wp), DIMENSION(max_number_of_particle_groups) ::  pdy = 9999999.9_wp            !< namelist parameter (see documentation)
    REAL(wp), DIMENSION(max_number_of_particle_groups) ::  pdz = 9999999.9_wp            !< namelist parameter (see documentation)
    REAL(wp), DIMENSION(max_number_of_particle_groups) ::  psb = 9999999.9_wp            !< namelist parameter (see documentation)
    REAL(wp), DIMENSION(max_number_of_particle_groups) ::  psl = 9999999.9_wp            !< namelist parameter (see documentation)
    REAL(wp), DIMENSION(max_number_of_particle_groups) ::  psn = 9999999.9_wp            !< namelist parameter (see documentation)
    REAL(wp), DIMENSION(max_number_of_particle_groups) ::  psr = 9999999.9_wp            !< namelist parameter (see documentation)
    REAL(wp), DIMENSION(max_number_of_particle_groups) ::  pss = 9999999.9_wp            !< namelist parameter (see documentation)
    REAL(wp), DIMENSION(max_number_of_particle_groups) ::  pst = 9999999.9_wp            !< namelist parameter (see documentation).
    REAL(wp), DIMENSION(max_number_of_particle_groups) ::  radius = 9999999.9_wp         !< namelist parameter (see documentation)

    REAL(wp), DIMENSION(:), ALLOCATABLE     ::  log_z_z0   !< Precalculate LOG(z/z0)  

    INTEGER(iwp), PARAMETER ::  NR_2_direction_move = 10000 !<

#if defined( __parallel )
    INTEGER(iwp)            ::  nr_move_north               !<
    INTEGER(iwp)            ::  nr_move_south               !<

    TYPE(particle_type), DIMENSION(:), ALLOCATABLE ::  move_also_north
    TYPE(particle_type), DIMENSION(:), ALLOCATABLE ::  move_also_south
#endif

    REAL(wp) ::  epsilon_collision !<
    REAL(wp) ::  urms              !<

    REAL(wp), DIMENSION(:),   ALLOCATABLE ::  epsclass  !< dissipation rate class
    REAL(wp), DIMENSION(:),   ALLOCATABLE ::  radclass  !< radius class
    REAL(wp), DIMENSION(:),   ALLOCATABLE ::  winf      !<

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ec        !<
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ecf       !<
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  gck       !<
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  hkernel   !<
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  hwratio   !<

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  ckernel !<
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  u_t   !< u value of old timelevel t
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  v_t   !< v value of old timelevel t
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  w_t   !< w value of old timelevel t


    INTEGER(iwp), PARAMETER         ::  PHASE_INIT    = 1  !<
    INTEGER(iwp), PARAMETER, PUBLIC ::  PHASE_RELEASE = 2  !<

    SAVE

    PRIVATE

    PUBLIC lpm_parin,     &
           lpm_header,    &
           lpm_init_arrays,&
           lpm_init,      &
           lpm_actions,   &
           lpm_data_output_ptseries, &
           lpm_interaction_droplets_ptq, &
           lpm_rrd_local_particles, &
           lpm_wrd_local, &
           lpm_rrd_global, &
           lpm_wrd_global, &
           lpm_rrd_local, &
           lpm_check_parameters

    PUBLIC lagrangian_particle_model

    INTERFACE lpm_check_parameters
       MODULE PROCEDURE lpm_check_parameters
    END INTERFACE lpm_check_parameters

    INTERFACE lpm_parin
       MODULE PROCEDURE lpm_parin
    END INTERFACE lpm_parin

    INTERFACE lpm_header
       MODULE PROCEDURE lpm_header
    END INTERFACE lpm_header

    INTERFACE lpm_init_arrays
       MODULE PROCEDURE lpm_init_arrays
    END INTERFACE lpm_init_arrays
 
    INTERFACE lpm_init
       MODULE PROCEDURE lpm_init
    END INTERFACE lpm_init

    INTERFACE lpm_actions
       MODULE PROCEDURE lpm_actions
    END INTERFACE lpm_actions

    INTERFACE lpm_data_output_ptseries
       MODULE PROCEDURE lpm_data_output_ptseries
    END INTERFACE

    INTERFACE lpm_rrd_local_particles
       MODULE PROCEDURE lpm_rrd_local_particles
    END INTERFACE lpm_rrd_local_particles

    INTERFACE lpm_rrd_global
       MODULE PROCEDURE lpm_rrd_global
    END INTERFACE lpm_rrd_global

    INTERFACE lpm_rrd_local
       MODULE PROCEDURE lpm_rrd_local
    END INTERFACE lpm_rrd_local

    INTERFACE lpm_wrd_local
       MODULE PROCEDURE lpm_wrd_local
    END INTERFACE lpm_wrd_local

    INTERFACE lpm_wrd_global
       MODULE PROCEDURE lpm_wrd_global
    END INTERFACE lpm_wrd_global

    INTERFACE lpm_advec
       MODULE PROCEDURE lpm_advec
    END INTERFACE lpm_advec

    INTERFACE lpm_calc_liquid_water_content
       MODULE PROCEDURE lpm_calc_liquid_water_content
    END INTERFACE

    INTERFACE lpm_interaction_droplets_ptq
       MODULE PROCEDURE lpm_interaction_droplets_ptq
       MODULE PROCEDURE lpm_interaction_droplets_ptq_ij
    END INTERFACE lpm_interaction_droplets_ptq

    INTERFACE lpm_boundary_conds
       MODULE PROCEDURE lpm_boundary_conds
    END INTERFACE lpm_boundary_conds

    INTERFACE lpm_droplet_condensation
       MODULE PROCEDURE lpm_droplet_condensation
    END INTERFACE

    INTERFACE lpm_droplet_collision
       MODULE PROCEDURE lpm_droplet_collision
    END INTERFACE lpm_droplet_collision

    INTERFACE lpm_init_kernels
       MODULE PROCEDURE lpm_init_kernels
    END INTERFACE lpm_init_kernels

    INTERFACE lpm_splitting
       MODULE PROCEDURE lpm_splitting
    END INTERFACE lpm_splitting

    INTERFACE lpm_merging
       MODULE PROCEDURE lpm_merging
    END INTERFACE lpm_merging

    INTERFACE lpm_exchange_horiz
       MODULE PROCEDURE lpm_exchange_horiz
    END INTERFACE lpm_exchange_horiz

    INTERFACE lpm_move_particle
       MODULE PROCEDURE lpm_move_particle
    END INTERFACE lpm_move_particle

    INTERFACE realloc_particles_array
       MODULE PROCEDURE realloc_particles_array
    END INTERFACE realloc_particles_array

    INTERFACE dealloc_particles_array
       MODULE PROCEDURE dealloc_particles_array
    END INTERFACE dealloc_particles_array

    INTERFACE lpm_sort_and_delete
       MODULE PROCEDURE lpm_sort_and_delete
    END INTERFACE lpm_sort_and_delete

    INTERFACE lpm_sort_timeloop_done
       MODULE PROCEDURE lpm_sort_timeloop_done
    END INTERFACE lpm_sort_timeloop_done

    INTERFACE lpm_pack
       MODULE PROCEDURE lpm_pack
    END INTERFACE lpm_pack

 CONTAINS
 

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Parin for &particle_parameters for the Lagrangian particle model
!------------------------------------------------------------------------------!
 SUBROUTINE lpm_parin
 
    CHARACTER (LEN=80) ::  line  !<

    NAMELIST /particles_par/ &
       aero_species, &
       aero_type, &
       aero_weight, &
       alloc_factor, &
       bc_par_b, &
       bc_par_lr, &
       bc_par_ns, &
       bc_par_t, &
       collision_kernel, &
       curvature_solution_effects, &
       deallocate_memory, &
       density_ratio, &
       dissipation_classes, &
       dt_dopts, &
       dt_min_part, &
       dt_prel, &
       dt_write_particle_data, &
       end_time_prel, &
       initial_weighting_factor, &
       log_sigma, &
       max_number_particles_per_gridbox, &
       merging, &
       na, &
       number_concentration, &
       number_of_particle_groups, &
       number_particles_per_gridbox, &
       particles_per_point, &
       particle_advection_start, &
       particle_advection_interpolation, &
       particle_maximum_age, &
       pdx, &
       pdy, &
       pdz, &
       psb, &
       psl, &
       psn, &
       psr, &
       pss, &
       pst, &
       radius, &
       radius_classes, &
       radius_merge, &
       radius_split, &
       random_start_position, &
       read_particles_from_restartfile, &
       rm, &
       seed_follows_topography, &
       splitting, &
       splitting_factor, &
       splitting_factor_max, &
       splitting_function, &
       splitting_mode, &
       step_dealloc, &
       use_sgs_for_particles, &
       vertical_particle_advection, &
       weight_factor_merge, &
       weight_factor_split, &
       write_particle_statistics

       NAMELIST /particle_parameters/ &
       aero_species, &
       aero_type, &
       aero_weight, &
       alloc_factor, &
       bc_par_b, &
       bc_par_lr, &
       bc_par_ns, &
       bc_par_t, &
       collision_kernel, &
       curvature_solution_effects, &
       deallocate_memory, &
       density_ratio, &
       dissipation_classes, &
       dt_dopts, &
       dt_min_part, &
       dt_prel, &
       dt_write_particle_data, &
       end_time_prel, &
       initial_weighting_factor, &
       log_sigma, &
       max_number_particles_per_gridbox, &
       merging, &
       na, &
       number_concentration, &
       number_of_particle_groups, &
       number_particles_per_gridbox, &
       particles_per_point, &
       particle_advection_start, &
       particle_advection_interpolation, &
       particle_maximum_age, &
       pdx, &
       pdy, &
       pdz, &
       psb, &
       psl, &
       psn, &
       psr, &
       pss, &
       pst, &
       radius, &
       radius_classes, &
       radius_merge, &
       radius_split, &
       random_start_position, &
       read_particles_from_restartfile, &
       rm, &
       seed_follows_topography, &
       splitting, &
       splitting_factor, &
       splitting_factor_max, &
       splitting_function, &
       splitting_mode, &
       step_dealloc, &
       use_sgs_for_particles, &
       vertical_particle_advection, &
       weight_factor_merge, &
       weight_factor_split, &
       write_particle_statistics

!
!-- Position the namelist-file at the beginning (it was already opened in
!-- parin), search for the namelist-group of the package and position the
!-- file at this line. Do the same for each optionally used package.
    line = ' '
    
!
!-- Try to find particles package
    REWIND ( 11 )
    line = ' '
    DO   WHILE ( INDEX( line, '&particle_parameters' ) == 0 )
       READ ( 11, '(A)', END=12 )  line
    ENDDO
    BACKSPACE ( 11 )
!
!-- Read user-defined namelist
    READ ( 11, particle_parameters, ERR = 10 )
!
!-- Set flag that indicates that particles are switched on
    particle_advection = .TRUE.
    
    GOTO 14

10  BACKSPACE( 11 )
    READ( 11 , '(A)') line
    CALL parin_fail_message( 'particle_parameters', line )
!
!-- Try to find particles package (old namelist)
12  REWIND ( 11 )
    line = ' '
    DO WHILE ( INDEX( line, '&particles_par' ) == 0 )
       READ ( 11, '(A)', END=14 )  line
    ENDDO
    BACKSPACE ( 11 )
!
!-- Read user-defined namelist
    READ ( 11, particles_par, ERR = 13, END = 14 )

    message_string = 'namelist particles_par is deprecated and will be ' //    &
                     'removed in near future. Please use namelist ' //         &
                     'particle_parameters instead'
    CALL message( 'package_parin', 'PA0487', 0, 1, 0, 6, 0 )

!
!-- Set flag that indicates that particles are switched on
    particle_advection = .TRUE.

    GOTO 14

13    BACKSPACE( 11 )
       READ( 11 , '(A)') line
       CALL parin_fail_message( 'particles_par', line )

14 CONTINUE

 END SUBROUTINE lpm_parin
 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Writes used particle attributes in header file.
!------------------------------------------------------------------------------! 
 SUBROUTINE lpm_header ( io )

    CHARACTER (LEN=40) ::  output_format       !< netcdf format

    INTEGER(iwp) ::  i               !<
    INTEGER(iwp), INTENT(IN) ::  io  !< Unit of the output file


     IF ( humidity  .AND.  cloud_droplets )  THEN
       WRITE ( io, 433 )
       IF ( curvature_solution_effects )  WRITE ( io, 434 )
       IF ( collision_kernel /= 'none' )  THEN
          WRITE ( io, 435 )  TRIM( collision_kernel )
          IF ( collision_kernel(6:9) == 'fast' )  THEN
             WRITE ( io, 436 )  radius_classes, dissipation_classes
          ENDIF
       ELSE
          WRITE ( io, 437 )
       ENDIF
    ENDIF
 
    IF ( particle_advection )  THEN
!
!--    Particle attributes
       WRITE ( io, 480 )  particle_advection_start, TRIM(particle_advection_interpolation), &
                          dt_prel, bc_par_lr, &
                          bc_par_ns, bc_par_b, bc_par_t, particle_maximum_age, &
                          end_time_prel
       IF ( use_sgs_for_particles )  WRITE ( io, 488 )  dt_min_part
       IF ( random_start_position )  WRITE ( io, 481 )
       IF ( seed_follows_topography )  WRITE ( io, 496 )
       IF ( particles_per_point > 1 )  WRITE ( io, 489 )  particles_per_point
       WRITE ( io, 495 )  total_number_of_particles
       IF ( dt_write_particle_data /= 9999999.9_wp )  THEN
          WRITE ( io, 485 )  dt_write_particle_data
          IF ( netcdf_data_format > 1 )  THEN
             output_format = 'netcdf (64 bit offset) and binary'
          ELSE
             output_format = 'netcdf and binary'
          ENDIF
          IF ( netcdf_deflate == 0 )  THEN
             WRITE ( io, 344 )  output_format
          ELSE
             WRITE ( io, 354 )  TRIM( output_format ), netcdf_deflate
          ENDIF
       ENDIF
       IF ( dt_dopts /= 9999999.9_wp )  WRITE ( io, 494 )  dt_dopts
       IF ( write_particle_statistics )  WRITE ( io, 486 )

       WRITE ( io, 487 )  number_of_particle_groups

       DO  i = 1, number_of_particle_groups
          WRITE ( io, 490 )  i, radius(i)
          IF ( density_ratio(i) /= 0.0_wp )  THEN
             WRITE ( io, 491 )  density_ratio(i)
          ELSE
             WRITE ( io, 492 )
          ENDIF
          WRITE ( io, 493 )  psl(i), psr(i), pss(i), psn(i), psb(i), pst(i), &
                             pdx(i), pdy(i), pdz(i)
          IF ( .NOT. vertical_particle_advection(i) )  WRITE ( io, 482 )
       ENDDO

    ENDIF
    
344 FORMAT ('       Output format: ',A/)
354 FORMAT ('       Output format: ',A, '   compressed with level: ',I1/)

433 FORMAT ('    Cloud droplets treated explicitly using the Lagrangian part', &
                 'icle model')
434 FORMAT ('    Curvature and solution effecs are considered for growth of', &
                 ' droplets < 1.0E-6 m')
435 FORMAT ('    Droplet collision is handled by ',A,'-kernel')
436 FORMAT ('       Fast kernel with fixed radius- and dissipation classes ', &
                    'are used'/ &
            '          number of radius classes:       ',I3,'    interval ', &
                       '[1.0E-6,2.0E-4] m'/ &
            '          number of dissipation classes:   ',I2,'    interval ', &
                       '[0,1000] cm**2/s**3')
437 FORMAT ('    Droplet collision is switched off')

480 FORMAT ('    Particles:'/ &
            '    ---------'// &
            '       Particle advection is active (switched on at t = ', F7.1, &
                    ' s)'/ &
            '       Interpolation of particle velocities is done by using ', A, &
                    ' method'/ &
            '       Start of new particle generations every  ',F6.1,' s'/ &
            '       Boundary conditions: left/right: ', A, ' north/south: ', A/&
            '                            bottom:     ', A, ' top:         ', A/&
            '       Maximum particle age:                 ',F9.1,' s'/ &
            '       Advection stopped at t = ',F9.1,' s'/)
481 FORMAT ('       Particles have random start positions'/)
482 FORMAT ('          Particles are advected only horizontally'/)
485 FORMAT ('       Particle data are written on file every ', F9.1, ' s')
486 FORMAT ('       Particle statistics are written on file'/)
487 FORMAT ('       Number of particle groups: ',I2/)
488 FORMAT ('       SGS velocity components are used for particle advection'/ &
            '          minimum timestep for advection:', F8.5/)
489 FORMAT ('       Number of particles simultaneously released at each ', &
                    'point: ', I5/)
490 FORMAT ('       Particle group ',I2,':'/ &
            '          Particle radius: ',E10.3, 'm')
491 FORMAT ('          Particle inertia is activated'/ &
            '             density_ratio (rho_fluid/rho_particle) =',F6.3/)
492 FORMAT ('          Particles are advected only passively (no inertia)'/)
493 FORMAT ('          Boundaries of particle source: x:',F8.1,' - ',F8.1,' m'/&
            '                                         y:',F8.1,' - ',F8.1,' m'/&
            '                                         z:',F8.1,' - ',F8.1,' m'/&
            '          Particle distances:  dx = ',F8.1,' m  dy = ',F8.1, &
                       ' m  dz = ',F8.1,' m'/)
494 FORMAT ('       Output of particle time series in NetCDF format every ', &
                    F8.2,' s'/)
495 FORMAT ('       Number of particles in total domain: ',I10/)
496 FORMAT ('       Initial vertical particle positions are interpreted ', &
                    'as relative to the given topography')
    
 END SUBROUTINE lpm_header
 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Writes used particle attributes in header file.
!------------------------------------------------------------------------------!  
 SUBROUTINE lpm_check_parameters
 
!
!-- Collision kernels:
    SELECT CASE ( TRIM( collision_kernel ) )

       CASE ( 'hall', 'hall_fast' )
          hall_kernel = .TRUE.

       CASE ( 'wang', 'wang_fast' )
          wang_kernel = .TRUE.

       CASE ( 'none' )


       CASE DEFAULT
          message_string = 'unknown collision kernel: collision_kernel = "' // &
                           TRIM( collision_kernel ) // '"'
          CALL message( 'lpm_check_parameters', 'PA0350', 1, 2, 0, 6, 0 )

    END SELECT
    IF ( collision_kernel(6:9) == 'fast' )  use_kernel_tables = .TRUE.

!
!-- Subgrid scale velocites with the simple interpolation method for resolved
!-- velocites is not implemented for passive particles. However, for cloud
!-- it can be combined as the sgs-velocites for active particles are
!-- calculated differently, i.e. no subboxes are needed.
    IF ( .NOT. TRIM( particle_advection_interpolation ) == 'trilinear'  .AND.  &
       use_sgs_for_particles  .AND.  .NOT. cloud_droplets )  THEN
          message_string = 'subrgrid scale velocities in combination with ' // &
                           'simple interpolation method is not '            // &
                           'implemented'
          CALL message( 'lpm_check_parameters', 'PA0659', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( nested_run  .AND.  cloud_droplets )  THEN
       message_string = 'nested runs in combination with cloud droplets ' // &
                        'is not implemented'
          CALL message( 'lpm_check_parameters', 'PA0687', 1, 2, 0, 6, 0 )
    ENDIF


 END SUBROUTINE lpm_check_parameters
 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize arrays for lpm
!------------------------------------------------------------------------------!   
 SUBROUTINE lpm_init_arrays
 
    IF ( cloud_droplets )  THEN
!
!--    Liquid water content, change in liquid water content
       ALLOCATE ( ql_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                         &
                  ql_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
!--    Real volume of particles (with weighting), volume of particles
       ALLOCATE ( ql_v(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                         &
                  ql_vp(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ENDIF


    ALLOCATE( u_t(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                              &
              v_t(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                              &
              w_t(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
!
!-- Initialize values with current time step
    u_t = u
    v_t = v
    w_t = w
!
!--    Initial assignment of the pointers
    IF ( cloud_droplets )  THEN
       ql   => ql_1
       ql_c => ql_2
    ENDIF

 END SUBROUTINE lpm_init_arrays
 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize Lagrangian particle model
!------------------------------------------------------------------------------! 
 SUBROUTINE lpm_init

    INTEGER(iwp) ::  i                           !<
    INTEGER(iwp) ::  j                           !<
    INTEGER(iwp) ::  k                           !<

    REAL(wp) ::  div                             !<
    REAL(wp) ::  height_int                      !<
    REAL(wp) ::  height_p                        !<
    REAL(wp) ::  z_p                             !<
    REAL(wp) ::  z0_av_local                     !<

!
!-- In case of oceans runs, the vertical index calculations need an offset,
!-- because otherwise the k indices will become negative
    IF ( ocean_mode )  THEN
       offset_ocean_nzt    = nzt
       offset_ocean_nzt_m1 = nzt - 1
    ENDIF

!
!-- Define block offsets for dividing a gridcell in 8 sub cells
!-- See documentation for List of subgrid boxes
!-- See pack_and_sort in lpm_pack_arrays.f90 for assignment of the subgrid boxes
    block_offset(0) = block_offset_def ( 0, 0, 0)
    block_offset(1) = block_offset_def ( 0, 0,-1)
    block_offset(2) = block_offset_def ( 0,-1, 0)
    block_offset(3) = block_offset_def ( 0,-1,-1)
    block_offset(4) = block_offset_def (-1, 0, 0)
    block_offset(5) = block_offset_def (-1, 0,-1)
    block_offset(6) = block_offset_def (-1,-1, 0)
    block_offset(7) = block_offset_def (-1,-1,-1)
!
!-- Check the number of particle groups.
    IF ( number_of_particle_groups > max_number_of_particle_groups )  THEN
       WRITE( message_string, * ) 'max_number_of_particle_groups =',           &
                                  max_number_of_particle_groups ,              &
                                  '&number_of_particle_groups reset to ',      &
                                  max_number_of_particle_groups
       CALL message( 'lpm_init', 'PA0213', 0, 1, 0, 6, 0 )
       number_of_particle_groups = max_number_of_particle_groups
    ENDIF
!
!-- Check if downward-facing walls exist. This case, reflection boundary
!-- conditions (as well as subgrid-scale velocities) may do not work
!-- propably (not realized so far).
    IF ( surf_def_h(1)%ns >= 1 )  THEN
       WRITE( message_string, * ) 'Overhanging topography do not work '//      &
                                  'with particles'
       CALL message( 'lpm_init', 'PA0212', 0, 1, 0, 6, 0 )

    ENDIF

!
!-- Set default start positions, if necessary
    IF ( psl(1) == 9999999.9_wp )  psl(1) = 0.0_wp
    IF ( psr(1) == 9999999.9_wp )  psr(1) = ( nx +1 ) * dx
    IF ( pss(1) == 9999999.9_wp )  pss(1) = 0.0_wp
    IF ( psn(1) == 9999999.9_wp )  psn(1) = ( ny +1 ) * dy
    IF ( psb(1) == 9999999.9_wp )  psb(1) = zu(nz/2)
    IF ( pst(1) == 9999999.9_wp )  pst(1) = psb(1)

    IF ( pdx(1) == 9999999.9_wp  .OR.  pdx(1) == 0.0_wp )  pdx(1) = dx
    IF ( pdy(1) == 9999999.9_wp  .OR.  pdy(1) == 0.0_wp )  pdy(1) = dy
    IF ( pdz(1) == 9999999.9_wp  .OR.  pdz(1) == 0.0_wp )  pdz(1) = zu(2) - zu(1)

!
!-- If number_particles_per_gridbox is set, the parametres pdx, pdy and pdz are
!-- calculated diagnostically. Therfore an isotropic distribution is prescribed.
    IF ( number_particles_per_gridbox /= -1 .AND.   &
         number_particles_per_gridbox >= 1 )    THEN
       pdx(1) = (( dx * dy * ( zu(2) - zu(1) ) ) /  &
             REAL(number_particles_per_gridbox))**0.3333333_wp
!
!--    Ensure a smooth value (two significant digits) of distance between
!--    particles (pdx, pdy, pdz).
       div = 1000.0_wp
       DO  WHILE ( pdx(1) < div )
          div = div / 10.0_wp
       ENDDO
       pdx(1) = NINT( pdx(1) * 100.0_wp / div ) * div / 100.0_wp
       pdy(1) = pdx(1)
       pdz(1) = pdx(1)

    ENDIF

    DO  j = 2, number_of_particle_groups
       IF ( psl(j) == 9999999.9_wp )  psl(j) = psl(j-1)
       IF ( psr(j) == 9999999.9_wp )  psr(j) = psr(j-1)
       IF ( pss(j) == 9999999.9_wp )  pss(j) = pss(j-1)
       IF ( psn(j) == 9999999.9_wp )  psn(j) = psn(j-1)
       IF ( psb(j) == 9999999.9_wp )  psb(j) = psb(j-1)
       IF ( pst(j) == 9999999.9_wp )  pst(j) = pst(j-1)
       IF ( pdx(j) == 9999999.9_wp  .OR.  pdx(j) == 0.0_wp )  pdx(j) = pdx(j-1)
       IF ( pdy(j) == 9999999.9_wp  .OR.  pdy(j) == 0.0_wp )  pdy(j) = pdy(j-1)
       IF ( pdz(j) == 9999999.9_wp  .OR.  pdz(j) == 0.0_wp )  pdz(j) = pdz(j-1)
    ENDDO

!
!-- Allocate arrays required for calculating particle SGS velocities.
!-- Initialize prefactor required for stoachastic Weil equation.
    IF ( use_sgs_for_particles  .AND.  .NOT. cloud_droplets )  THEN
       ALLOCATE( de_dx(nzb:nzt+1,nysg:nyng,nxlg:nxrg), &
                 de_dy(nzb:nzt+1,nysg:nyng,nxlg:nxrg), &
                 de_dz(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

       de_dx = 0.0_wp
       de_dy = 0.0_wp
       de_dz = 0.0_wp

       sgs_wf_part = 1.0_wp / 3.0_wp
    ENDIF

!
!-- Allocate array required for logarithmic vertical interpolation of
!-- horizontal particle velocities between the surface and the first vertical
!-- grid level. In order to avoid repeated CPU cost-intensive CALLS of
!-- intrinsic FORTRAN procedure LOG(z/z0), LOG(z/z0) is precalculated for
!-- several heights. Splitting into 20 sublayers turned out to be sufficient.
!-- To obtain exact height levels of particles, linear interpolation is applied
!-- (see lpm_advec.f90).
    IF ( constant_flux_layer )  THEN

       ALLOCATE ( log_z_z0(0:number_of_sublayers) )
       z_p = zu(nzb+1) - zw(nzb)

!
!--    Calculate horizontal mean value of z0 used for logartihmic
!--    interpolation. Note: this is not exact for heterogeneous z0.
!--    However, sensitivity studies showed that the effect is
!--    negligible.
       z0_av_local  = SUM( surf_def_h(0)%z0 ) + SUM( surf_lsm_h%z0 ) +         &
                      SUM( surf_usm_h%z0 )
       z0_av_global = 0.0_wp

#if defined( __parallel )
       CALL MPI_ALLREDUCE(z0_av_local, z0_av_global, 1, MPI_REAL, MPI_SUM, &
                          comm2d, ierr )
#else
       z0_av_global = z0_av_local
#endif

       z0_av_global = z0_av_global  / ( ( ny + 1 ) * ( nx + 1 ) )
!
!--    Horizontal wind speed is zero below and at z0
       log_z_z0(0) = 0.0_wp
!
!--    Calculate vertical depth of the sublayers
       height_int  = ( z_p - z0_av_global ) / REAL( number_of_sublayers, KIND=wp )
!
!--    Precalculate LOG(z/z0)
       height_p    = z0_av_global
       DO  k = 1, number_of_sublayers

          height_p    = height_p + height_int
          log_z_z0(k) = LOG( height_p / z0_av_global )

       ENDDO

    ENDIF

!
!-- Check which particle interpolation method should be used
    IF ( TRIM( particle_advection_interpolation )  ==  'trilinear' )  THEN
       interpolation_simple_corrector = .FALSE.
       interpolation_simple_predictor = .FALSE.
       interpolation_trilinear        = .TRUE.
    ELSEIF ( TRIM( particle_advection_interpolation )  ==  'simple_corrector' )  THEN
       interpolation_simple_corrector = .TRUE.
       interpolation_simple_predictor = .FALSE.
       interpolation_trilinear        = .FALSE.
    ELSEIF ( TRIM( particle_advection_interpolation )  ==  'simple_predictor' )  THEN
       interpolation_simple_corrector = .FALSE.
       interpolation_simple_predictor = .TRUE.
       interpolation_trilinear        = .FALSE.
    ENDIF

!
!-- Check boundary condition and set internal variables
    SELECT CASE ( bc_par_b )

       CASE ( 'absorb' )
          ibc_par_b = 1

       CASE ( 'reflect' )
          ibc_par_b = 2

       CASE ( 'reset' )
          ibc_par_b = 3

       CASE DEFAULT
          WRITE( message_string, * )  'unknown boundary condition ',           &
                                       'bc_par_b = "', TRIM( bc_par_b ), '"'
          CALL message( 'lpm_init', 'PA0217', 1, 2, 0, 6, 0 )

    END SELECT
    SELECT CASE ( bc_par_t )

       CASE ( 'absorb' )
          ibc_par_t = 1

       CASE ( 'reflect' )
          ibc_par_t = 2

       CASE ( 'nested' )
          ibc_par_t = 3

       CASE DEFAULT
          WRITE( message_string, * ) 'unknown boundary condition ',            &
                                     'bc_par_t = "', TRIM( bc_par_t ), '"'
          CALL message( 'lpm_init', 'PA0218', 1, 2, 0, 6, 0 )

    END SELECT
    SELECT CASE ( bc_par_lr )

       CASE ( 'cyclic' )
          ibc_par_lr = 0

       CASE ( 'absorb' )
          ibc_par_lr = 1

       CASE ( 'reflect' )
          ibc_par_lr = 2

       CASE ( 'nested' )
          ibc_par_lr = 3

       CASE DEFAULT
          WRITE( message_string, * ) 'unknown boundary condition ',   &
                                     'bc_par_lr = "', TRIM( bc_par_lr ), '"'
          CALL message( 'lpm_init', 'PA0219', 1, 2, 0, 6, 0 )

    END SELECT
    SELECT CASE ( bc_par_ns )

       CASE ( 'cyclic' )
          ibc_par_ns = 0

       CASE ( 'absorb' )
          ibc_par_ns = 1

       CASE ( 'reflect' )
          ibc_par_ns = 2

       CASE ( 'nested' )
          ibc_par_ns = 3

       CASE DEFAULT
          WRITE( message_string, * ) 'unknown boundary condition ',   &
                                     'bc_par_ns = "', TRIM( bc_par_ns ), '"'
          CALL message( 'lpm_init', 'PA0220', 1, 2, 0, 6, 0 )

    END SELECT
    SELECT CASE ( splitting_mode )

       CASE ( 'const' )
          i_splitting_mode = 1

       CASE ( 'cl_av' )
          i_splitting_mode = 2

       CASE ( 'gb_av' )
          i_splitting_mode = 3

       CASE DEFAULT
          WRITE( message_string, * )  'unknown splitting_mode = "',            &
                                      TRIM( splitting_mode ), '"'
          CALL message( 'lpm_init', 'PA0146', 1, 2, 0, 6, 0 )

    END SELECT
    SELECT CASE ( splitting_function )

       CASE ( 'gamma' )
          isf = 1

       CASE ( 'log' )
          isf = 2

       CASE ( 'exp' )
          isf = 3

       CASE DEFAULT
          WRITE( message_string, * )  'unknown splitting function = "',        &
                                       TRIM( splitting_function ), '"'
          CALL message( 'lpm_init', 'PA0147', 1, 2, 0, 6, 0 )

    END SELECT
!
!-- Initialize collision kernels
    IF ( collision_kernel /= 'none' )  CALL lpm_init_kernels
!
!-- For the first model run of a possible job chain initialize the
!-- particles, otherwise read the particle data from restart file.
    IF ( TRIM( initializing_actions ) == 'read_restart_data'  &
         .AND.  read_particles_from_restartfile )  THEN
       CALL lpm_rrd_local_particles
    ELSE
!
!--    Allocate particle arrays and set attributes of the initial set of
!--    particles, which can be also periodically released at later times.
       ALLOCATE( prt_count(nzb:nzt+1,nysg:nyng,nxlg:nxrg), &
                 grid_particles(nzb+1:nzt,nys:nyn,nxl:nxr) )

       number_of_particles = 0
       prt_count           = 0
!
!--    initialize counter for particle IDs
       grid_particles%id_counter = 1
!
!--    Initialize all particles with dummy values (otherwise errors may
!--    occur within restart runs). The reason for this is still not clear
!--    and may be presumably caused by errors in the respective user-interface.
       zero_particle = particle_type( 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,  &
                                      0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,  &
                                      0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,  &
                                      0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,  &
                                      0, 0, 0_idp, .FALSE., -1 )

       particle_groups = particle_groups_type( 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp )
!
!--    Set values for the density ratio and radius for all particle
!--    groups, if necessary
       IF ( density_ratio(1) == 9999999.9_wp )  density_ratio(1) = 0.0_wp
       IF ( radius(1)        == 9999999.9_wp )  radius(1) = 0.0_wp
       DO  i = 2, number_of_particle_groups
          IF ( density_ratio(i) == 9999999.9_wp )  THEN
             density_ratio(i) = density_ratio(i-1)
          ENDIF
          IF ( radius(i) == 9999999.9_wp )  radius(i) = radius(i-1)
       ENDDO

       DO  i = 1, number_of_particle_groups
          IF ( density_ratio(i) /= 0.0_wp  .AND.  radius(i) == 0 )  THEN
             WRITE( message_string, * ) 'particle group #', i, ' has a',       &
                                        'density ratio /= 0 but radius = 0'
             CALL message( 'lpm_init', 'PA0215', 1, 2, 0, 6, 0 )
          ENDIF
          particle_groups(i)%density_ratio = density_ratio(i)
          particle_groups(i)%radius        = radius(i)
       ENDDO
!
!--    Set a seed value for the random number generator to be exclusively
!--    used for the particle code. The generated random numbers should be
!--    different on the different PEs.
       iran_part = iran_part + myid
!
!--    Create the particle set, and set the initial particles
       CALL lpm_create_particle( phase_init )
       last_particle_release_time = particle_advection_start
!
!--    User modification of initial particles
       CALL user_lpm_init
!
!--    Open file for statistical informations about particle conditions
       IF ( write_particle_statistics )  THEN
          CALL check_open( 80 )
          WRITE ( 80, 8000 )  current_timestep_number, simulated_time,         &
                              number_of_particles
          CALL close_file( 80 )
       ENDIF

    ENDIF

#if defined( __parallel )
    IF ( nested_run )  CALL pmcp_g_init
#endif

!
!-- To avoid programm abort, assign particles array to the local version of
!-- first grid cell
    number_of_particles = prt_count(nzb+1,nys,nxl)
    particles => grid_particles(nzb+1,nys,nxl)%particles(1:number_of_particles)
!
!-- Formats
8000 FORMAT (I6,1X,F7.2,4X,I10,71X,I10)

 END SUBROUTINE lpm_init
 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Create Lagrangian particles
!------------------------------------------------------------------------------!  
 SUBROUTINE lpm_create_particle (phase)

    INTEGER(iwp)               ::  alloc_size  !< relative increase of allocated memory for particles
    INTEGER(iwp)               ::  i           !< loop variable ( particle groups )
    INTEGER(iwp)               ::  ip          !< index variable along x
    INTEGER(iwp)               ::  j           !< loop variable ( particles per point )
    INTEGER(iwp)               ::  jp          !< index variable along y
    INTEGER(iwp)               ::  k           !< index variable along z
    INTEGER(iwp)               ::  k_surf      !< index of surface grid point
    INTEGER(iwp)               ::  kp          !< index variable along z
    INTEGER(iwp)               ::  loop_stride !< loop variable for initialization
    INTEGER(iwp)               ::  n           !< loop variable ( number of particles )
    INTEGER(iwp)               ::  new_size    !< new size of allocated memory for particles

    INTEGER(iwp), INTENT(IN)   ::  phase       !< mode of inititialization

    INTEGER(iwp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  local_count !< start address of new particle
    INTEGER(iwp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  local_start !< start address of new particle

    LOGICAL                    ::  first_stride !< flag for initialization

    REAL(wp)                   ::  pos_x      !< increment for particle position in x
    REAL(wp)                   ::  pos_y      !< increment for particle position in y
    REAL(wp)                   ::  pos_z      !< increment for particle position in z
    REAL(wp)                   ::  rand_contr !< dummy argument for random position

    TYPE(particle_type),TARGET ::  tmp_particle !< temporary particle used for initialization


!
!-- Calculate particle positions and store particle attributes, if
!-- particle is situated on this PE
    DO  loop_stride = 1, 2
       first_stride = (loop_stride == 1)
       IF ( first_stride )   THEN
          local_count = 0           ! count number of particles
       ELSE
          local_count = prt_count   ! Start address of new particles
       ENDIF

!
!--    Calculate initial_weighting_factor diagnostically
       IF ( number_concentration /= -1.0_wp  .AND.  number_concentration > 0.0_wp )  THEN
          initial_weighting_factor =  number_concentration  *                           &
                                      pdx(1) * pdy(1) * pdz(1)
       END IF

       n = 0
       DO  i = 1, number_of_particle_groups
          pos_z = psb(i)
          DO WHILE ( pos_z <= pst(i) )
             IF ( pos_z >= zw(0) .AND.  pos_z < zw(nzt) )  THEN
                pos_y = pss(i)
                DO WHILE ( pos_y <= psn(i) )
                   IF ( pos_y >= nys * dy  .AND.                  &
                        pos_y <  ( nyn + 1 ) * dy  )  THEN
                      pos_x = psl(i)
               xloop: DO WHILE ( pos_x <= psr(i) )
                         IF ( pos_x >= nxl * dx  .AND.            &
                              pos_x <  ( nxr + 1) * dx )  THEN
                            DO  j = 1, particles_per_point
                               n = n + 1
                               tmp_particle%x             = pos_x
                               tmp_particle%y             = pos_y
                               tmp_particle%z             = pos_z
                               tmp_particle%age           = 0.0_wp
                               tmp_particle%age_m         = 0.0_wp
                               tmp_particle%dt_sum        = 0.0_wp
                               tmp_particle%e_m           = 0.0_wp
                               tmp_particle%rvar1         = 0.0_wp
                               tmp_particle%rvar2         = 0.0_wp
                               tmp_particle%rvar3         = 0.0_wp
                               tmp_particle%speed_x       = 0.0_wp
                               tmp_particle%speed_y       = 0.0_wp
                               tmp_particle%speed_z       = 0.0_wp
                               tmp_particle%origin_x      = pos_x
                               tmp_particle%origin_y      = pos_y
                               tmp_particle%origin_z      = pos_z
                               IF ( curvature_solution_effects )  THEN
                                  tmp_particle%aux1      = 0.0_wp    ! dry aerosol radius
                                  tmp_particle%aux2      = dt_3d     ! last Rosenbrock timestep
                               ELSE
                                  tmp_particle%aux1      = 0.0_wp    ! free to use
                                  tmp_particle%aux2      = 0.0_wp    ! free to use
                               ENDIF
                               tmp_particle%radius        = particle_groups(i)%radius
                               tmp_particle%weight_factor = initial_weighting_factor
                               tmp_particle%class         = 1
                               tmp_particle%group         = i
                               tmp_particle%id            = 0_idp
                               tmp_particle%particle_mask = .TRUE.
                               tmp_particle%block_nr      = -1
!
!--                            Determine the grid indices of the particle position
                               ip = INT( tmp_particle%x * ddx )
                               jp = INT( tmp_particle%y * ddy )
!
!--                            In case of stretching the actual k index is found iteratively
                               IF ( dz_stretch_level /= -9999999.9_wp  .OR.           &
                                    dz_stretch_level_start(1) /= -9999999.9_wp )  THEN
                                  kp = MINLOC( ABS( tmp_particle%z - zu ), DIM = 1 ) - 1
                               ELSE
                                  kp = INT( tmp_particle%z / dz(1) + 1 + offset_ocean_nzt )
                               ENDIF
!
!--                            Determine surface level. Therefore, check for
!--                            upward-facing wall on w-grid. 
                               k_surf = topo_top_ind(jp,ip,3)
                               IF ( seed_follows_topography )  THEN
!
!--                               Particle height is given relative to topography
                                  kp = kp + k_surf
                                  tmp_particle%z = tmp_particle%z + zw(k_surf)
!--                               Skip particle release if particle position is
!--                               above model top, or within topography in case
!--                               of overhanging structures.
                                  IF ( kp > nzt  .OR.                          &
                                 .NOT. BTEST( wall_flags_total_0(kp,jp,ip), 0 ) )  THEN
                                     pos_x = pos_x + pdx(i)
                                     CYCLE xloop
                                  ENDIF
!
!--                            Skip particle release if particle position is
!--                            below surface, or within topography in case
!--                            of overhanging structures.
                               ELSEIF ( .NOT. seed_follows_topography .AND.    &
                                         tmp_particle%z <= zw(k_surf)  .OR.    &
                                        .NOT. BTEST( wall_flags_total_0(kp,jp,ip), 0 ) )&
                               THEN
                                  pos_x = pos_x + pdx(i)
                                  CYCLE xloop
                               ENDIF

                               local_count(kp,jp,ip) = local_count(kp,jp,ip) + 1

                               IF ( .NOT. first_stride )  THEN
                                  IF ( ip < nxl  .OR.  jp < nys  .OR.  kp < nzb+1 )  THEN
                                     write(6,*) 'xl ',ip,jp,kp,nxl,nys,nzb+1
                                  ENDIF
                                  IF ( ip > nxr  .OR.  jp > nyn  .OR.  kp > nzt )  THEN
                                     write(6,*) 'xu ',ip,jp,kp,nxr,nyn,nzt
                                  ENDIF
                                  grid_particles(kp,jp,ip)%particles(local_count(kp,jp,ip)) = tmp_particle
                               ENDIF
                            ENDDO
                         ENDIF
                         pos_x = pos_x + pdx(i)
                      ENDDO xloop
                   ENDIF
                   pos_y = pos_y + pdy(i)
                ENDDO
             ENDIF

             pos_z = pos_z + pdz(i)
          ENDDO
       ENDDO

       IF ( first_stride )  THEN
          DO  ip = nxl, nxr
             DO  jp = nys, nyn
                DO  kp = nzb+1, nzt
                   IF ( phase == PHASE_INIT )  THEN
                      IF ( local_count(kp,jp,ip) > 0 )  THEN
                         alloc_size = MAX( INT( local_count(kp,jp,ip) *        &
                            ( 1.0_wp + alloc_factor / 100.0_wp ) ),            &
                            1 )
                      ELSE
                         alloc_size = 1
                      ENDIF
                      ALLOCATE(grid_particles(kp,jp,ip)%particles(1:alloc_size))
                      DO  n = 1, alloc_size
                         grid_particles(kp,jp,ip)%particles(n) = zero_particle
                      ENDDO
                   ELSEIF ( phase == PHASE_RELEASE )  THEN
                      IF ( local_count(kp,jp,ip) > 0 )  THEN
                         new_size   = local_count(kp,jp,ip) + prt_count(kp,jp,ip)
                         alloc_size = MAX( INT( new_size * ( 1.0_wp +          &
                            alloc_factor / 100.0_wp ) ), 1 )
                         IF( alloc_size > SIZE( grid_particles(kp,jp,ip)%particles) )  THEN
                            CALL realloc_particles_array( ip, jp, kp, alloc_size )
                         ENDIF
                      ENDIF
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDIF

    ENDDO

    local_start = prt_count+1
    prt_count   = local_count
!
!-- Calculate particle IDs
    DO  ip = nxl, nxr
       DO  jp = nys, nyn
          DO  kp = nzb+1, nzt
             number_of_particles = prt_count(kp,jp,ip)
             IF ( number_of_particles <= 0 )  CYCLE
             particles => grid_particles(kp,jp,ip)%particles(1:number_of_particles)

             DO  n = local_start(kp,jp,ip), number_of_particles  !only new particles

                particles(n)%id = 10000_idp**3 * grid_particles(kp,jp,ip)%id_counter + &
                                  10000_idp**2 * kp + 10000_idp * jp + ip
!
!--             Count the number of particles that have been released before
                grid_particles(kp,jp,ip)%id_counter =                          &
                                         grid_particles(kp,jp,ip)%id_counter + 1

             ENDDO

          ENDDO
       ENDDO
    ENDDO
!
!-- Initialize aerosol background spectrum
    IF ( curvature_solution_effects )  THEN
       CALL lpm_init_aerosols( local_start )
    ENDIF
!
!-- Add random fluctuation to particle positions.
    IF ( random_start_position )  THEN
       DO  ip = nxl, nxr
          DO  jp = nys, nyn
             DO  kp = nzb+1, nzt
                number_of_particles = prt_count(kp,jp,ip)
                IF ( number_of_particles <= 0 )  CYCLE
                particles => grid_particles(kp,jp,ip)%particles(1:number_of_particles)
!
!--             Move only new particles. Moreover, limit random fluctuation
!--             in order to prevent that particles move more than one grid box,
!--             which would lead to problems concerning particle exchange
!--             between processors in case pdx/pdy are larger than dx/dy,
!--             respectively.
                DO  n = local_start(kp,jp,ip), number_of_particles
                   IF ( psl(particles(n)%group) /= psr(particles(n)%group) )  THEN
                      rand_contr = ( random_function( iran_part ) - 0.5_wp ) * &
                                     pdx(particles(n)%group)
                      particles(n)%x = particles(n)%x +                        &
                              MERGE( rand_contr, SIGN( dx, rand_contr ),       &
                                     ABS( rand_contr ) < dx                    &
                                   )
                   ENDIF
                   IF ( pss(particles(n)%group) /= psn(particles(n)%group) )  THEN
                      rand_contr = ( random_function( iran_part ) - 0.5_wp ) * &
                                     pdy(particles(n)%group)
                      particles(n)%y = particles(n)%y +                        &
                              MERGE( rand_contr, SIGN( dy, rand_contr ),       &
                                     ABS( rand_contr ) < dy                    &
                                   )
                   ENDIF
                   IF ( psb(particles(n)%group) /= pst(particles(n)%group) )  THEN
                      rand_contr = ( random_function( iran_part ) - 0.5_wp ) * &
                                     pdz(particles(n)%group)
                      particles(n)%z = particles(n)%z +                        &
                              MERGE( rand_contr, SIGN( dzw(kp), rand_contr ),  &
                                     ABS( rand_contr ) < dzw(kp)               &
                                   )
                   ENDIF
                ENDDO
!
!--             Identify particles located outside the model domain and reflect
!--             or absorb them if necessary.
                CALL lpm_boundary_conds( 'bottom/top', i, j, k )
!
!--             Furthermore, remove particles located in topography. Note, as
!--             the particle speed is still zero at this point, wall
!--             reflection boundary conditions will not work in this case.
                particles =>                                                   &
                       grid_particles(kp,jp,ip)%particles(1:number_of_particles)
                DO  n = local_start(kp,jp,ip), number_of_particles
                   i = particles(n)%x * ddx
                   j = particles(n)%y * ddy
                   k = particles(n)%z / dz(1) + 1 + offset_ocean_nzt
                   DO WHILE( zw(k) < particles(n)%z )
                      k = k + 1
                   ENDDO
                   DO WHILE( zw(k-1) > particles(n)%z )
                      k = k - 1
                   ENDDO
!
!--                Check if particle is within topography
                   IF ( .NOT. BTEST( wall_flags_total_0(k,j,i), 0 ) )  THEN
                      particles(n)%particle_mask = .FALSE.
                      deleted_particles = deleted_particles + 1
                   ENDIF

                ENDDO
             ENDDO
          ENDDO
       ENDDO
!
!--    Exchange particles between grid cells and processors
       CALL lpm_move_particle
       CALL lpm_exchange_horiz

    ENDIF
!
!-- In case of random_start_position, delete particles identified by
!-- lpm_exchange_horiz and lpm_boundary_conds. Then sort particles into blocks,
!-- which is needed for a fast interpolation of the LES fields on the particle
!-- position.
    CALL lpm_sort_and_delete
!
!-- Determine the current number of particles
    DO  ip = nxl, nxr
       DO  jp = nys, nyn
          DO  kp = nzb+1, nzt
             number_of_particles         = number_of_particles                 &
                                           + prt_count(kp,jp,ip)
          ENDDO
       ENDDO
    ENDDO
!
!-- Calculate the number of particles of the total domain
#if defined( __parallel )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( number_of_particles, total_number_of_particles, 1, &
    MPI_INTEGER, MPI_SUM, comm2d, ierr )
#else
    total_number_of_particles = number_of_particles
#endif

    RETURN

 END SUBROUTINE lpm_create_particle
 
 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine initialize the particles as aerosols with physio-chemical 
!> properties. 
!------------------------------------------------------------------------------!   
 SUBROUTINE lpm_init_aerosols(local_start)

    REAL(wp) ::  afactor            !< curvature effects
    REAL(wp) ::  bfactor            !< solute effects
    REAL(wp) ::  dlogr              !< logarithmic width of radius bin
    REAL(wp) ::  e_a                !< vapor pressure
    REAL(wp) ::  e_s                !< saturation vapor pressure
    REAL(wp) ::  rmin = 0.005e-6_wp !< minimum aerosol radius
    REAL(wp) ::  rmax = 10.0e-6_wp  !< maximum aerosol radius
    REAL(wp) ::  r_mid              !< mean radius of bin
    REAL(wp) ::  r_l                !< left radius of bin
    REAL(wp) ::  r_r                !< right radius of bin
    REAL(wp) ::  sigma              !< surface tension
    REAL(wp) ::  t_int              !< temperature

    INTEGER(iwp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg), INTENT(IN) ::  local_start !<

    INTEGER(iwp) ::  n              !<
    INTEGER(iwp) ::  ip             !<
    INTEGER(iwp) ::  jp             !<
    INTEGER(iwp) ::  kp             !<

!
!-- Set constants for different aerosol species
    IF ( TRIM( aero_species ) == 'nacl' )  THEN
       molecular_weight_of_solute = 0.05844_wp 
       rho_s                      = 2165.0_wp
       vanthoff                   = 2.0_wp
    ELSEIF ( TRIM( aero_species ) == 'c3h4o4' )  THEN
       molecular_weight_of_solute = 0.10406_wp 
       rho_s                      = 1600.0_wp
       vanthoff                   = 1.37_wp
    ELSEIF ( TRIM( aero_species ) == 'nh4o3' )  THEN
       molecular_weight_of_solute = 0.08004_wp 
       rho_s                      = 1720.0_wp
       vanthoff                   = 2.31_wp
    ELSE
       WRITE( message_string, * ) 'unknown aerosol species ',   &
                                'aero_species = "', TRIM( aero_species ), '"'
       CALL message( 'lpm_init', 'PA0470', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- The following typical aerosol spectra are taken from Jaenicke (1993):
!-- Tropospheric aerosols. Published in Aerosol-Cloud-Climate Interactions.
    IF ( TRIM( aero_type ) == 'polar' )  THEN
       na        = (/ 2.17e1, 1.86e-1, 3.04e-4 /) * 1.0E6_wp
       rm        = (/ 0.0689, 0.375, 4.29 /) * 1.0E-6_wp
       log_sigma = (/ 0.245, 0.300, 0.291 /)
    ELSEIF ( TRIM( aero_type ) == 'background' )  THEN
       na        = (/ 1.29e2, 5.97e1, 6.35e1 /) * 1.0E6_wp
       rm        = (/ 0.0036, 0.127, 0.259 /) * 1.0E-6_wp
       log_sigma = (/ 0.645, 0.253, 0.425 /)
    ELSEIF ( TRIM( aero_type ) == 'maritime' )  THEN
       na        = (/ 1.33e2, 6.66e1, 3.06e0 /) * 1.0E6_wp
       rm        = (/ 0.0039, 0.133, 0.29 /) * 1.0E-6_wp
       log_sigma = (/ 0.657, 0.210, 0.396 /)
    ELSEIF ( TRIM( aero_type ) == 'continental' )  THEN
       na        = (/ 3.20e3, 2.90e3, 3.00e-1 /) * 1.0E6_wp
       rm        = (/ 0.01, 0.058, 0.9 /) * 1.0E-6_wp
       log_sigma = (/ 0.161, 0.217, 0.380 /)
    ELSEIF ( TRIM( aero_type ) == 'desert' )  THEN
       na        = (/ 7.26e2, 1.14e3, 1.78e-1 /) * 1.0E6_wp
       rm        = (/ 0.001, 0.0188, 10.8 /) * 1.0E-6_wp
       log_sigma = (/ 0.247, 0.770, 0.438 /)
    ELSEIF ( TRIM( aero_type ) == 'rural' )  THEN
       na        = (/ 6.65e3, 1.47e2, 1.99e3 /) * 1.0E6_wp
       rm        = (/ 0.00739, 0.0269, 0.0419 /) * 1.0E-6_wp
       log_sigma = (/ 0.225, 0.557, 0.266 /)
    ELSEIF ( TRIM( aero_type ) == 'urban' )  THEN
       na        = (/ 9.93e4, 1.11e3, 3.64e4 /) * 1.0E6_wp
       rm        = (/ 0.00651, 0.00714, 0.0248 /) * 1.0E-6_wp
       log_sigma = (/ 0.245, 0.666, 0.337 /)
    ELSEIF ( TRIM( aero_type ) == 'user' )  THEN
       CONTINUE
    ELSE
       WRITE( message_string, * ) 'unknown aerosol type ',   &
                                'aero_type = "', TRIM( aero_type ), '"'
       CALL message( 'lpm_init', 'PA0459', 1, 2, 0, 6, 0 )
    ENDIF

    DO  ip = nxl, nxr
       DO  jp = nys, nyn
          DO  kp = nzb+1, nzt

             number_of_particles = prt_count(kp,jp,ip)
             IF ( number_of_particles <= 0 )  CYCLE
             particles => grid_particles(kp,jp,ip)%particles(1:number_of_particles)

             dlogr   = ( LOG10(rmax) - LOG10(rmin) ) / ( number_of_particles - local_start(kp,jp,ip) + 1 )
!
!--          Initialize the aerosols with a predefined spectral distribution
!--          of the dry radius (logarithmically increasing bins) and a varying
!--          weighting factor
             DO  n = local_start(kp,jp,ip), number_of_particles  !only new particles

                r_l   = 10.0**( LOG10( rmin ) + (n-1) * dlogr )
                r_r   = 10.0**( LOG10( rmin ) + n * dlogr )
                r_mid = SQRT( r_l * r_r )

                particles(n)%aux1          = r_mid
                particles(n)%weight_factor =                                           &
                   ( na(1) / ( SQRT( 2.0_wp * pi ) * log_sigma(1) ) *                     &
                     EXP( - LOG10( r_mid / rm(1) )**2 / ( 2.0_wp * log_sigma(1)**2 ) ) +  &
                     na(2) / ( SQRT( 2.0_wp * pi ) * log_sigma(2) ) *                     &
                     EXP( - LOG10( r_mid / rm(2) )**2 / ( 2.0_wp * log_sigma(2)**2 ) ) +  &
                     na(3) / ( SQRT( 2.0_wp * pi ) * log_sigma(3) ) *                     &
                     EXP( - LOG10( r_mid / rm(3) )**2 / ( 2.0_wp * log_sigma(3)**2 ) )    &
                   ) * ( LOG10(r_r) - LOG10(r_l) ) * ( dx * dy * dzw(kp) )

!
!--             Multiply weight_factor with the namelist parameter aero_weight
!--             to increase or decrease the number of simulated aerosols
                particles(n)%weight_factor = particles(n)%weight_factor * aero_weight

                IF ( particles(n)%weight_factor - FLOOR(particles(n)%weight_factor,KIND=wp) &
                     > random_function( iran_part ) )  THEN
                   particles(n)%weight_factor = FLOOR(particles(n)%weight_factor,KIND=wp) + 1.0_wp
                ELSE
                   particles(n)%weight_factor = FLOOR(particles(n)%weight_factor,KIND=wp)
                ENDIF
!
!--             Unnecessary particles will be deleted
                IF ( particles(n)%weight_factor <= 0.0_wp )  particles(n)%particle_mask = .FALSE.

             ENDDO
!
!--          Set particle radius to equilibrium radius based on the environmental
!--          supersaturation (Khvorostyanov and Curry, 2007, JGR). This avoids
!--          the sometimes lengthy growth toward their equilibrium radius within
!--          the simulation.
             t_int  = pt(kp,jp,ip) * exner(kp)

             e_s = magnus( t_int )
             e_a = q(kp,jp,ip) * hyp(kp) / ( q(kp,jp,ip) + rd_d_rv )

             sigma   = 0.0761_wp - 0.000155_wp * ( t_int - 273.15_wp )
             afactor = 2.0_wp * sigma / ( rho_l * r_v * t_int )

             bfactor = vanthoff * molecular_weight_of_water *    &
                       rho_s / ( molecular_weight_of_solute * rho_l )
!
!--          The formula is only valid for subsaturated environments. For
!--          supersaturations higher than -5 %, the supersaturation is set to -5%.
             IF ( e_a / e_s >= 0.95_wp )  e_a = 0.95_wp * e_s

             DO  n = local_start(kp,jp,ip), number_of_particles  !only new particles
!
!--             For details on this equation, see Eq. (14) of Khvorostyanov and
!--             Curry (2007, JGR)
                particles(n)%radius = bfactor**0.3333333_wp *                  &
                   particles(n)%aux1 / ( 1.0_wp - e_a / e_s )**0.3333333_wp / &
                   ( 1.0_wp + ( afactor / ( 3.0_wp * bfactor**0.3333333_wp *   &
                     particles(n)%aux1 ) ) /                                  &
                     ( 1.0_wp - e_a / e_s )**0.6666666_wp                      &
                   )

             ENDDO

          ENDDO
       ENDDO
    ENDDO

 END SUBROUTINE lpm_init_aerosols


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculates quantities required for considering the SGS velocity fluctuations
!> in the particle transport by a stochastic approach. The respective
!> quantities are: SGS-TKE gradients and horizontally averaged profiles of the
!> SGS TKE and the resolved-scale velocity variances. 
!------------------------------------------------------------------------------!
 SUBROUTINE lpm_init_sgs_tke

    USE exchange_horiz_mod,                                                    &
        ONLY:  exchange_horiz

    USE statistics,                                                            &
        ONLY:  flow_statistics_called, hom, sums, sums_l

    INTEGER(iwp) ::  i      !< index variable along x
    INTEGER(iwp) ::  j      !< index variable along y
    INTEGER(iwp) ::  k      !< index variable along z
    INTEGER(iwp) ::  m      !< running index for the surface elements 

    REAL(wp) ::  flag1      !< flag to mask topography

!
!-- TKE gradient along x and y
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb, nzt+1

             IF ( .NOT. BTEST( wall_flags_total_0(k,j,i-1), 0 )  .AND.         &
                        BTEST( wall_flags_total_0(k,j,i), 0   )  .AND.         &
                        BTEST( wall_flags_total_0(k,j,i+1), 0 ) )              &
             THEN
                de_dx(k,j,i) = 2.0_wp * sgs_wf_part *                          &
                               ( e(k,j,i+1) - e(k,j,i) ) * ddx
             ELSEIF ( BTEST( wall_flags_total_0(k,j,i-1), 0 )  .AND.           &
                      BTEST( wall_flags_total_0(k,j,i), 0   )  .AND.           &
                .NOT. BTEST( wall_flags_total_0(k,j,i+1), 0 ) )                &
             THEN
                de_dx(k,j,i) = 2.0_wp * sgs_wf_part *                          &
                               ( e(k,j,i) - e(k,j,i-1) ) * ddx
             ELSEIF ( .NOT. BTEST( wall_flags_total_0(k,j,i), 22   )  .AND.    &
                      .NOT. BTEST( wall_flags_total_0(k,j,i+1), 22 ) )         &   
             THEN
                de_dx(k,j,i) = 0.0_wp
             ELSEIF ( .NOT. BTEST( wall_flags_total_0(k,j,i-1), 22 )  .AND.    &
                      .NOT. BTEST( wall_flags_total_0(k,j,i), 22   ) )         &
             THEN
                de_dx(k,j,i) = 0.0_wp
             ELSE
                de_dx(k,j,i) = sgs_wf_part * ( e(k,j,i+1) - e(k,j,i-1) ) * ddx
             ENDIF

             IF ( .NOT. BTEST( wall_flags_total_0(k,j-1,i), 0 )  .AND.         &
                        BTEST( wall_flags_total_0(k,j,i), 0   )  .AND.         &
                        BTEST( wall_flags_total_0(k,j+1,i), 0 ) )              &
             THEN
                de_dy(k,j,i) = 2.0_wp * sgs_wf_part *                          &
                               ( e(k,j+1,i) - e(k,j,i) ) * ddy
             ELSEIF ( BTEST( wall_flags_total_0(k,j-1,i), 0 )  .AND.           &
                      BTEST( wall_flags_total_0(k,j,i), 0   )  .AND.           &
                .NOT. BTEST( wall_flags_total_0(k,j+1,i), 0 ) )                &
             THEN
                de_dy(k,j,i) = 2.0_wp * sgs_wf_part *                          &
                               ( e(k,j,i) - e(k,j-1,i) ) * ddy
             ELSEIF ( .NOT. BTEST( wall_flags_total_0(k,j,i), 22   )  .AND.    &
                      .NOT. BTEST( wall_flags_total_0(k,j+1,i), 22 ) )         &   
             THEN
                de_dy(k,j,i) = 0.0_wp
             ELSEIF ( .NOT. BTEST( wall_flags_total_0(k,j-1,i), 22 )  .AND.    &
                      .NOT. BTEST( wall_flags_total_0(k,j,i), 22   ) )         &
             THEN
                de_dy(k,j,i) = 0.0_wp
             ELSE
                de_dy(k,j,i) = sgs_wf_part * ( e(k,j+1,i) - e(k,j-1,i) ) * ddy
             ENDIF

          ENDDO
       ENDDO
    ENDDO

!
!-- TKE gradient along z at topograhy and  including bottom and top boundary conditions
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt-1
!
!--          Flag to mask topography
             flag1 = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0  ) )

             de_dz(k,j,i) = 2.0_wp * sgs_wf_part *                             &
                           ( e(k+1,j,i) - e(k-1,j,i) ) / ( zu(k+1) - zu(k-1) ) &
                                                 * flag1 
          ENDDO
!
!--       upward-facing surfaces
          DO  m = bc_h(0)%start_index(j,i), bc_h(0)%end_index(j,i)
             k            = bc_h(0)%k(m)
             de_dz(k,j,i) = 2.0_wp * sgs_wf_part *                             &
                           ( e(k+1,j,i) - e(k,j,i)   ) / ( zu(k+1) - zu(k) )
          ENDDO
!
!--       downward-facing surfaces
          DO  m = bc_h(1)%start_index(j,i), bc_h(1)%end_index(j,i)
             k            = bc_h(1)%k(m)
             de_dz(k,j,i) = 2.0_wp * sgs_wf_part *                             &
                           ( e(k,j,i) - e(k-1,j,i)   ) / ( zu(k) - zu(k-1) )
          ENDDO

          de_dz(nzb,j,i)   = 0.0_wp
          de_dz(nzt,j,i)   = 0.0_wp
          de_dz(nzt+1,j,i) = 0.0_wp
       ENDDO
    ENDDO
!
!-- Ghost point exchange
    CALL exchange_horiz( de_dx, nbgp )
    CALL exchange_horiz( de_dy, nbgp )
    CALL exchange_horiz( de_dz, nbgp )
    CALL exchange_horiz( diss, nbgp  )
!
!-- Set boundary conditions at non-periodic boundaries. Note, at non-period
!-- boundaries zero-gradient boundary conditions are set for the subgrid TKE.
!-- Thus, TKE gradients normal to the respective lateral boundaries are zero, 
!-- while tangetial TKE gradients then must be the same as within the prognostic
!-- domain.  
    IF ( bc_dirichlet_l )  THEN
       de_dx(:,:,-1) = 0.0_wp
       de_dy(:,:,-1) = de_dy(:,:,0) 
       de_dz(:,:,-1) = de_dz(:,:,0)
    ENDIF
    IF ( bc_dirichlet_r )  THEN
       de_dx(:,:,nxr+1) = 0.0_wp
       de_dy(:,:,nxr+1) = de_dy(:,:,nxr) 
       de_dz(:,:,nxr+1) = de_dz(:,:,nxr)
    ENDIF
    IF ( bc_dirichlet_n )  THEN
       de_dx(:,nyn+1,:) = de_dx(:,nyn,:)
       de_dy(:,nyn+1,:) = 0.0_wp 
       de_dz(:,nyn+1,:) = de_dz(:,nyn,:)
    ENDIF
    IF ( bc_dirichlet_s )  THEN
       de_dx(:,nys-1,:) = de_dx(:,nys,:)
       de_dy(:,nys-1,:) = 0.0_wp 
       de_dz(:,nys-1,:) = de_dz(:,nys,:)
    ENDIF  
!
!-- Calculate the horizontally averaged profiles of SGS TKE and resolved
!-- velocity variances (they may have been already calculated in routine
!-- flow_statistics).
    IF ( .NOT. flow_statistics_called )  THEN

!
!--    First calculate horizontally averaged profiles of the horizontal
!--    velocities.
       sums_l(:,1,0) = 0.0_wp
       sums_l(:,2,0) = 0.0_wp

       DO  i = nxl, nxr
          DO  j =  nys, nyn
             DO  k = nzb, nzt+1
!
!--             Flag indicating vicinity of wall
                flag1 = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 24 ) )

                sums_l(k,1,0)  = sums_l(k,1,0)  + u(k,j,i) * flag1
                sums_l(k,2,0)  = sums_l(k,2,0)  + v(k,j,i) * flag1
             ENDDO
          ENDDO
       ENDDO

#if defined( __parallel )
!
!--    Compute total sum from local sums
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( sums_l(nzb,1,0), sums(nzb,1), nzt+2-nzb, &
                           MPI_REAL, MPI_SUM, comm2d, ierr )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( sums_l(nzb,2,0), sums(nzb,2), nzt+2-nzb, &
                              MPI_REAL, MPI_SUM, comm2d, ierr )
#else
       sums(:,1) = sums_l(:,1,0)
       sums(:,2) = sums_l(:,2,0)
#endif

!
!--    Final values are obtained by division by the total number of grid
!--    points used for the summation.
       hom(:,1,1,0) = sums(:,1) / ngp_2dh_outer(:,0)   ! u
       hom(:,1,2,0) = sums(:,2) / ngp_2dh_outer(:,0)   ! v

!
!--    Now calculate the profiles of SGS TKE and the resolved-scale
!--    velocity variances
       sums_l(:,8,0)  = 0.0_wp
       sums_l(:,30,0) = 0.0_wp
       sums_l(:,31,0) = 0.0_wp
       sums_l(:,32,0) = 0.0_wp
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb, nzt+1
!
!--             Flag indicating vicinity of wall
                flag1 = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 24 ) )

                sums_l(k,8,0)  = sums_l(k,8,0)  + e(k,j,i)                       * flag1
                sums_l(k,30,0) = sums_l(k,30,0) + ( u(k,j,i) - hom(k,1,1,0) )**2 * flag1
                sums_l(k,31,0) = sums_l(k,31,0) + ( v(k,j,i) - hom(k,1,2,0) )**2 * flag1
                sums_l(k,32,0) = sums_l(k,32,0) + w(k,j,i)**2                    * flag1
             ENDDO
          ENDDO
       ENDDO

#if defined( __parallel )
!
!--    Compute total sum from local sums
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( sums_l(nzb,8,0), sums(nzb,8), nzt+2-nzb, &
                           MPI_REAL, MPI_SUM, comm2d, ierr )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( sums_l(nzb,30,0), sums(nzb,30), nzt+2-nzb, &
                           MPI_REAL, MPI_SUM, comm2d, ierr )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( sums_l(nzb,31,0), sums(nzb,31), nzt+2-nzb, &
                           MPI_REAL, MPI_SUM, comm2d, ierr )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( sums_l(nzb,32,0), sums(nzb,32), nzt+2-nzb, &
                           MPI_REAL, MPI_SUM, comm2d, ierr )

#else
       sums(:,8)  = sums_l(:,8,0)
       sums(:,30) = sums_l(:,30,0)
       sums(:,31) = sums_l(:,31,0)
       sums(:,32) = sums_l(:,32,0)
#endif

!
!--    Final values are obtained by division by the total number of grid
!--    points used for the summation.
       hom(:,1,8,0)  = sums(:,8)  / ngp_2dh_outer(:,0)   ! e
       hom(:,1,30,0) = sums(:,30) / ngp_2dh_outer(:,0)   ! u*2
       hom(:,1,31,0) = sums(:,31) / ngp_2dh_outer(:,0)   ! v*2 
       hom(:,1,32,0) = sums(:,32) / ngp_2dh_outer(:,0)   ! w*2

    ENDIF

 END SUBROUTINE lpm_init_sgs_tke
 
 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Sobroutine control lpm actions, i.e. all actions during one time step. 
!------------------------------------------------------------------------------!  
 SUBROUTINE lpm_actions( location )

    USE exchange_horiz_mod,                                                    &
        ONLY:  exchange_horiz

    CHARACTER (LEN=*), INTENT(IN) ::  location !< call location string

    INTEGER(iwp)       ::  i                  !<
    INTEGER(iwp)       ::  ie                 !<
    INTEGER(iwp)       ::  is                 !<
    INTEGER(iwp)       ::  j                  !<
    INTEGER(iwp)       ::  je                 !<
    INTEGER(iwp)       ::  js                 !<
    INTEGER(iwp), SAVE ::  lpm_count = 0      !<
    INTEGER(iwp)       ::  k                  !<
    INTEGER(iwp)       ::  ke                 !<
    INTEGER(iwp)       ::  ks                 !<
    INTEGER(iwp)       ::  m                  !<
    INTEGER(iwp), SAVE ::  steps = 0          !<

    LOGICAL            ::  first_loop_stride  !<


    SELECT CASE ( location )

       CASE ( 'after_pressure_solver' )
!
!--       The particle model is executed if particle advection start is reached and only at the end
!--       of the intermediate time step loop.
          IF ( time_since_reference_point >= particle_advection_start   &
               .AND.  intermediate_timestep_count == intermediate_timestep_count_max )             &
          THEN
             CALL cpu_log( log_point(25), 'lpm', 'start' )
!
!--          Write particle data at current time on file.
!--          This has to be done here, before particles are further processed,
!--          because they may be deleted within this timestep (in case that
!--          dt_write_particle_data = dt_prel = particle_maximum_age).
             time_write_particle_data = time_write_particle_data + dt_3d
             IF ( time_write_particle_data >= dt_write_particle_data )  THEN

                CALL lpm_data_output_particles
!
!--          The MOD function allows for changes in the output interval with restart
!--          runs.
                time_write_particle_data = MOD( time_write_particle_data, &
                                           MAX( dt_write_particle_data, dt_3d ) )
             ENDIF

!
!--          Initialize arrays for marking those particles to be deleted after the
!--          (sub-) timestep
             deleted_particles = 0

!
!--          Initialize variables used for accumulating the number of particles
!--          xchanged between the subdomains during all sub-timesteps (if sgs
!--          velocities are included). These data are output further below on the
!--          particle statistics file.
             trlp_count_sum      = 0
             trlp_count_recv_sum = 0
             trrp_count_sum      = 0
             trrp_count_recv_sum = 0
             trsp_count_sum      = 0
             trsp_count_recv_sum = 0
             trnp_count_sum      = 0
             trnp_count_recv_sum = 0
!
!--          Calculate exponential term used in case of particle inertia for each
!--          of the particle groups
             DO  m = 1, number_of_particle_groups
                IF ( particle_groups(m)%density_ratio /= 0.0_wp )  THEN
                   particle_groups(m)%exp_arg  =                                        &
                             4.5_wp * particle_groups(m)%density_ratio *                &
                             molecular_viscosity / ( particle_groups(m)%radius )**2

                   particle_groups(m)%exp_term = EXP( -particle_groups(m)%exp_arg *     &
                             dt_3d )
                ENDIF
             ENDDO
!
!--          If necessary, release new set of particles
             IF ( ( simulated_time - last_particle_release_time ) >= dt_prel  .AND.     &
                    end_time_prel > simulated_time )  THEN
                DO WHILE ( ( simulated_time - last_particle_release_time ) >= dt_prel )
                   CALL lpm_create_particle( PHASE_RELEASE )
                   last_particle_release_time = last_particle_release_time + dt_prel
                ENDDO
             ENDIF
!
!--          Reset summation arrays
             IF ( cloud_droplets )  THEN
                ql_c  = 0.0_wp
                ql_v  = 0.0_wp
                ql_vp = 0.0_wp
             ENDIF

             first_loop_stride = .TRUE.
             grid_particles(:,:,:)%time_loop_done = .TRUE.
!
!--          Timestep loop for particle advection.
!--          This loop has to be repeated until the advection time of every particle
!--          (within the total domain!) has reached the LES timestep (dt_3d).
!--          In case of including the SGS velocities, the particle timestep may be
!--          smaller than the LES timestep (because of the Lagrangian timescale
!--          restriction) and particles may require to undergo several particle
!--          timesteps, before the LES timestep is reached. Because the number of these
!--          particle timesteps to be carried out is unknown at first, these steps are
!--          carried out in the following infinite loop with exit condition.
             DO
                CALL cpu_log( log_point_s(44), 'lpm_advec', 'start' )
                CALL cpu_log( log_point_s(44), 'lpm_advec', 'pause' )

!
!--             If particle advection includes SGS velocity components, calculate the
!--             required SGS quantities (i.e. gradients of the TKE, as well as
!--             horizontally averaged profiles of the SGS TKE and the resolved-scale
!--             velocity variances)
                IF ( use_sgs_for_particles  .AND.  .NOT. cloud_droplets )  THEN
                   CALL lpm_init_sgs_tke
                ENDIF
!
!--             In case SGS-particle speed is considered, particles may carry out
!--             several particle timesteps. In order to prevent unnecessary
!--             treatment of particles that already reached the final time level,
!--             particles are sorted into contiguous blocks of finished and
!--             not-finished particles, in addition to their already sorting
!--             according to their sub-boxes.
                IF ( .NOT. first_loop_stride  .AND.  use_sgs_for_particles )            &
                   CALL lpm_sort_timeloop_done
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt

                         number_of_particles = prt_count(k,j,i)
!
!--                      If grid cell gets empty, flag must be true
                         IF ( number_of_particles <= 0 )  THEN
                            grid_particles(k,j,i)%time_loop_done = .TRUE.
                            CYCLE
                         ENDIF

                         IF ( .NOT. first_loop_stride  .AND.  &
                              grid_particles(k,j,i)%time_loop_done )  CYCLE

                         particles => grid_particles(k,j,i)%particles(1:number_of_particles)

                         particles(1:number_of_particles)%particle_mask = .TRUE.
!
!--                      Initialize the variable storing the total time that a particle
!--                      has advanced within the timestep procedure
                         IF ( first_loop_stride )  THEN
                            particles(1:number_of_particles)%dt_sum = 0.0_wp
                         ENDIF
!
!--                      Particle (droplet) growth by condensation/evaporation and
!--                      collision
                         IF ( cloud_droplets  .AND.  first_loop_stride)  THEN
!
!--                         Droplet growth by condensation / evaporation
                            CALL lpm_droplet_condensation(i,j,k)
!
!--                         Particle growth by collision
                            IF ( collision_kernel /= 'none' )  THEN
                               CALL lpm_droplet_collision(i,j,k)
                            ENDIF

                         ENDIF
!
!--                      Initialize the switch used for the loop exit condition checked
!--                      at the end of this loop. If at least one particle has failed to
!--                      reach the LES timestep, this switch will be set false in
!--                      lpm_advec.
                         dt_3d_reached_l = .TRUE.

!
!--                      Particle advection
                         CALL lpm_advec( i, j, k )
!
!--                      Particle reflection from walls. Only applied if the particles
!--                      are in the vertical range of the topography. (Here, some
!--                      optimization is still possible.)
                         IF ( topography /= 'flat'  .AND.  k < nzb_max + 2 )  THEN
                            CALL  lpm_boundary_conds( 'walls', i, j, k )
                         ENDIF
!
!--                      User-defined actions after the calculation of the new particle
!--                      position
                         CALL user_lpm_advec( i, j, k )
!
!--                      Apply boundary conditions to those particles that have crossed
!--                      the top or bottom boundary and delete those particles, which are
!--                      older than allowed
                         CALL lpm_boundary_conds( 'bottom/top', i, j, k )
!
!---                     If not all particles of the actual grid cell have reached the
!--                      LES timestep, this cell has to do another loop iteration. Due to
!--                      the fact that particles can move into neighboring grid cells,
!--                      these neighbor cells also have to perform another loop iteration.
!--                      Please note, this realization does not work properly if
!--                      particles move into another subdomain.
                         IF ( .NOT. dt_3d_reached_l )  THEN
                            ks = MAX(nzb+1,k-1)
                            ke = MIN(nzt,k+1)
                            js = MAX(nys,j-1)
                            je = MIN(nyn,j+1)
                            is = MAX(nxl,i-1)
                            ie = MIN(nxr,i+1)
                            grid_particles(ks:ke,js:je,is:ie)%time_loop_done = .FALSE.
                         ELSE
                            grid_particles(k,j,i)%time_loop_done = .TRUE.
                         ENDIF

                      ENDDO
                   ENDDO
                ENDDO
                steps = steps + 1
                dt_3d_reached_l = ALL(grid_particles(:,:,:)%time_loop_done)
!
!--             Find out, if all particles on every PE have completed the LES timestep
!--             and set the switch corespondingly
#if defined( __parallel )
                IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
                CALL MPI_ALLREDUCE( dt_3d_reached_l, dt_3d_reached, 1, MPI_LOGICAL, &
                                    MPI_LAND, comm2d, ierr )
#else
                dt_3d_reached = dt_3d_reached_l
#endif
                CALL cpu_log( log_point_s(44), 'lpm_advec', 'stop' )

!
!--             Apply splitting and merging algorithm
                IF ( cloud_droplets )  THEN
                   IF ( splitting )  THEN
                      CALL lpm_splitting
                   ENDIF
                   IF ( merging )  THEN
                      CALL lpm_merging
                   ENDIF
                ENDIF
!
!--             Move Particles local to PE to a different grid cell
                CALL lpm_move_particle
!
!--             Horizontal boundary conditions including exchange between subdmains
                CALL lpm_exchange_horiz

!
!--             IF .FALSE., lpm_sort_and_delete is done inside pcmp
                IF ( .NOT. dt_3d_reached  .OR.  .NOT. nested_run )   THEN
!
!--                Pack particles (eliminate those marked for deletion),
!--                determine new number of particles
                   CALL lpm_sort_and_delete

!--                Initialize variables for the next (sub-) timestep, i.e., for marking
!--                those particles to be deleted after the timestep
                   deleted_particles = 0
                ENDIF

                IF ( dt_3d_reached )  EXIT

                first_loop_stride = .FALSE.
             ENDDO   ! timestep loop

#if defined( __parallel )
!
!--          in case of nested runs do the transfer of particles after every full model time step
             IF ( nested_run )   THEN
                CALL particles_from_parent_to_child
                CALL particles_from_child_to_parent
                CALL pmcp_p_delete_particles_in_fine_grid_area

                CALL lpm_sort_and_delete

                deleted_particles = 0
             ENDIF
#endif

!
!--          Calculate the new liquid water content for each grid box
             IF ( cloud_droplets )  CALL lpm_calc_liquid_water_content

!
!--          At the end all arrays are exchanged
             IF ( cloud_droplets )  THEN
                CALL exchange_horiz( ql, nbgp )
                CALL exchange_horiz( ql_c, nbgp )
                CALL exchange_horiz( ql_v, nbgp )
                CALL exchange_horiz( ql_vp, nbgp )
             ENDIF

!
!--          Deallocate unused memory
             IF ( deallocate_memory  .AND.  lpm_count == step_dealloc )  THEN
                CALL dealloc_particles_array
                lpm_count = 0
             ELSEIF ( deallocate_memory )  THEN
                lpm_count = lpm_count + 1
             ENDIF

!
!--          Write particle statistics (in particular the number of particles
!--          exchanged between the subdomains) on file
             IF ( write_particle_statistics )  CALL lpm_write_exchange_statistics
!
!--          Execute Interactions of condnesation and evaporation to humidity and
!--          temperature field
             IF ( cloud_droplets )  THEN
                CALL lpm_interaction_droplets_ptq
                CALL exchange_horiz( pt, nbgp )
                CALL exchange_horiz( q, nbgp )
             ENDIF

             CALL cpu_log( log_point(25), 'lpm', 'stop' )

! !
! !--       Output of particle time series
!           IF ( particle_advection )  THEN
!              IF ( time_dopts >= dt_dopts  .OR.                                                        &
!                   ( time_since_reference_point >= particle_advection_start  .AND.                     &
!                    first_call_lpm ) )  THEN
!                 CALL lpm_data_output_ptseries
!                 time_dopts = MOD( time_dopts, MAX( dt_dopts, dt_3d ) )
!              ENDIF
!           ENDIF

!
!--           Set this switch to .false. @todo: maybe find better solution.
              first_call_lpm = .FALSE.
           ENDIF! ENDIF statement of lpm_actions('after_pressure_solver')

       CASE ( 'after_integration' )
!
!--       Call at the end of timestep routine to save particle velocities fields
!--       for the next timestep
          CALL lpm_swap_timelevel_for_particle_advection

       CASE DEFAULT
          CONTINUE

    END SELECT

 END SUBROUTINE lpm_actions
 

#if defined( __parallel )
!------------------------------------------------------------------------------!
! Description:
! ------------
!
!------------------------------------------------------------------------------!
 SUBROUTINE particles_from_parent_to_child

    CALL pmcp_c_get_particle_from_parent                         ! Child actions
    CALL pmcp_p_fill_particle_win                                ! Parent actions

    RETURN

 END SUBROUTINE particles_from_parent_to_child

 
!------------------------------------------------------------------------------!
! Description:
! ------------
!
!------------------------------------------------------------------------------!
 SUBROUTINE particles_from_child_to_parent

    CALL pmcp_c_send_particle_to_parent                         ! Child actions
    CALL pmcp_p_empty_particle_win                              ! Parent actions

    RETURN

 END SUBROUTINE particles_from_child_to_parent
#endif
 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine write exchange statistics of the lpm in a ascii file. 
!------------------------------------------------------------------------------!
 SUBROUTINE lpm_write_exchange_statistics

    INTEGER(iwp) ::  ip         !<
    INTEGER(iwp) ::  jp         !<
    INTEGER(iwp) ::  kp         !<
    INTEGER(iwp) ::  tot_number_of_particles !<

!
!-- Determine the current number of particles
    number_of_particles         = 0
    DO  ip = nxl, nxr
       DO  jp = nys, nyn
          DO  kp = nzb+1, nzt
             number_of_particles = number_of_particles                         &
                                     + prt_count(kp,jp,ip)
          ENDDO
       ENDDO
    ENDDO

    CALL check_open( 80 )
#if defined( __parallel )
    WRITE ( 80, 8000 )  current_timestep_number+1, simulated_time+dt_3d, &
                        number_of_particles, pleft, trlp_count_sum,      &
                        trlp_count_recv_sum, pright, trrp_count_sum,     &
                        trrp_count_recv_sum, psouth, trsp_count_sum,     &
                        trsp_count_recv_sum, pnorth, trnp_count_sum,     &
                        trnp_count_recv_sum
#else
    WRITE ( 80, 8000 )  current_timestep_number+1, simulated_time+dt_3d, &
                        number_of_particles
#endif
    CALL close_file( 80 )

    IF ( number_of_particles > 0 )  THEN
        WRITE(9,*) 'number_of_particles ', number_of_particles,                &
                    current_timestep_number + 1, simulated_time + dt_3d
    ENDIF

#if defined( __parallel )
    CALL MPI_ALLREDUCE( number_of_particles, tot_number_of_particles, 1,       &
                        MPI_INTEGER, MPI_SUM, comm2d, ierr )
#else
    tot_number_of_particles = number_of_particles
#endif

#if defined( __parallel )
    IF ( nested_run )  THEN
       CALL pmcp_g_print_number_of_particles( simulated_time+dt_3d,            &
                                              tot_number_of_particles)
    ENDIF
#endif

!
!-- Formats
8000 FORMAT (I6,1X,F7.2,4X,I10,5X,4(I3,1X,I4,'/',I4,2X),6X,I10)


 END SUBROUTINE lpm_write_exchange_statistics
 

!------------------------------------------------------------------------------! 
! Description:
! ------------
!> Write particle data in FORTRAN binary and/or netCDF format 
!------------------------------------------------------------------------------!
 SUBROUTINE lpm_data_output_particles
 
    INTEGER(iwp) ::  ip !<
    INTEGER(iwp) ::  jp !<
    INTEGER(iwp) ::  kp !<

    CALL cpu_log( log_point_s(40), 'lpm_data_output', 'start' )

!
!-- Attention: change version number for unit 85 (in routine check_open)
!--            whenever the output format for this unit is changed!
    CALL check_open( 85 )

    WRITE ( 85 )  simulated_time
    WRITE ( 85 )  prt_count

    DO  ip = nxl, nxr
       DO  jp = nys, nyn
          DO  kp = nzb+1, nzt
             number_of_particles = prt_count(kp,jp,ip)
             particles => grid_particles(kp,jp,ip)%particles(1:number_of_particles)
             IF ( number_of_particles <= 0 )  CYCLE
             WRITE ( 85 )  particles
          ENDDO
       ENDDO
    ENDDO

    CALL close_file( 85 )


#if defined( __netcdf )
! !
! !-- Output in netCDF format
!     CALL check_open( 108 )
! 
! !
! !-- Update the NetCDF time axis
!     prt_time_count = prt_time_count + 1
! 
!     nc_stat = NF90_PUT_VAR( id_set_prt, id_var_time_prt, &
!                             (/ simulated_time /),        &
!                             start = (/ prt_time_count /), count = (/ 1 /) )
!     CALL netcdf_handle_error( 'lpm_data_output_particles', 1 )
! 
! !
! !-- Output the real number of particles used
!     nc_stat = NF90_PUT_VAR( id_set_prt, id_var_rnop_prt, &
!                             (/ number_of_particles /),   &
!                             start = (/ prt_time_count /), count = (/ 1 /) )
!     CALL netcdf_handle_error( 'lpm_data_output_particles', 2 )
! 
! !
! !-- Output all particle attributes
!     nc_stat = NF90_PUT_VAR( id_set_prt, id_var_prt(1), particles%age,      &
!                             start = (/ 1, prt_time_count /),               &
!                             count = (/ maximum_number_of_particles /) )
!     CALL netcdf_handle_error( 'lpm_data_output_particles', 3 )
! 
!     nc_stat = NF90_PUT_VAR( id_set_prt, id_var_prt(2), particles%user,     &
!                             start = (/ 1, prt_time_count /),               &
!                             count = (/ maximum_number_of_particles /) )
!     CALL netcdf_handle_error( 'lpm_data_output_particles', 4 )
! 
!     nc_stat = NF90_PUT_VAR( id_set_prt, id_var_prt(3), particles%origin_x, &
!                             start = (/ 1, prt_time_count /),               &
!                             count = (/ maximum_number_of_particles /) )
!     CALL netcdf_handle_error( 'lpm_data_output_particles', 5 )
! 
!     nc_stat = NF90_PUT_VAR( id_set_prt, id_var_prt(4), particles%origin_y, &
!                             start = (/ 1, prt_time_count /),               &
!                             count = (/ maximum_number_of_particles /) )
!     CALL netcdf_handle_error( 'lpm_data_output_particles', 6 )
! 
!     nc_stat = NF90_PUT_VAR( id_set_prt, id_var_prt(5), particles%origin_z, &
!                             start = (/ 1, prt_time_count /),               &
!                             count = (/ maximum_number_of_particles /) )
!     CALL netcdf_handle_error( 'lpm_data_output_particles', 7 )
! 
!     nc_stat = NF90_PUT_VAR( id_set_prt, id_var_prt(6), particles%radius,   &
!                             start = (/ 1, prt_time_count /),               &
!                             count = (/ maximum_number_of_particles /) )
!     CALL netcdf_handle_error( 'lpm_data_output_particles', 8 )
! 
!     nc_stat = NF90_PUT_VAR( id_set_prt, id_var_prt(7), particles%speed_x,  &
!                             start = (/ 1, prt_time_count /),               &
!                             count = (/ maximum_number_of_particles /) )
!     CALL netcdf_handle_error( 'lpm_data_output_particles', 9 )
! 
!     nc_stat = NF90_PUT_VAR( id_set_prt, id_var_prt(8), particles%speed_y,  &
!                             start = (/ 1, prt_time_count /),               &
!                             count = (/ maximum_number_of_particles /) )
!     CALL netcdf_handle_error( 'lpm_data_output_particles', 10 )
! 
!     nc_stat = NF90_PUT_VAR( id_set_prt, id_var_prt(9), particles%speed_z,  &
!                             start = (/ 1, prt_time_count /),               &
!                             count = (/ maximum_number_of_particles /) )
!     CALL netcdf_handle_error( 'lpm_data_output_particles', 11 )
! 
!     nc_stat = NF90_PUT_VAR( id_set_prt,id_var_prt(10),                     &
!                             particles%weight_factor,                       &
!                             start = (/ 1, prt_time_count /),               &
!                             count = (/ maximum_number_of_particles /) )
!     CALL netcdf_handle_error( 'lpm_data_output_particles', 12 )
! 
!     nc_stat = NF90_PUT_VAR( id_set_prt, id_var_prt(11), particles%x,       &
!                             start = (/ 1, prt_time_count /),               &
!                             count = (/ maximum_number_of_particles /) )
!     CALL netcdf_handle_error( 'lpm_data_output_particles', 13 )
! 
!     nc_stat = NF90_PUT_VAR( id_set_prt, id_var_prt(12), particles%y,       &
!                             start = (/ 1, prt_time_count /),               &
!                             count = (/ maximum_number_of_particles /) )
!     CALL netcdf_handle_error( 'lpm_data_output_particles', 14 )
! 
!     nc_stat = NF90_PUT_VAR( id_set_prt, id_var_prt(13), particles%z,       &
!                             start = (/ 1, prt_time_count /),               &
!                             count = (/ maximum_number_of_particles /) )
!     CALL netcdf_handle_error( 'lpm_data_output_particles', 15 )
! 
!     nc_stat = NF90_PUT_VAR( id_set_prt, id_var_prt(14), particles%class,   &
!                             start = (/ 1, prt_time_count /),               &
!                             count = (/ maximum_number_of_particles /) )
!     CALL netcdf_handle_error( 'lpm_data_output_particles', 16 )
! 
!     nc_stat = NF90_PUT_VAR( id_set_prt, id_var_prt(15), particles%group,   &
!                             start = (/ 1, prt_time_count /),               &
!                             count = (/ maximum_number_of_particles /) )
!     CALL netcdf_handle_error( 'lpm_data_output_particles', 17 )
! 
!     nc_stat = NF90_PUT_VAR( id_set_prt, id_var_prt(16),                    &
!                             particles%id2,                                 &
!                             start = (/ 1, prt_time_count /),               &
!                             count = (/ maximum_number_of_particles /) )
!     CALL netcdf_handle_error( 'lpm_data_output_particles', 18 )
! 
!     nc_stat = NF90_PUT_VAR( id_set_prt, id_var_prt(17), particles%id1,     &
!                             start = (/ 1, prt_time_count /),               &
!                             count = (/ maximum_number_of_particles /) )
!     CALL netcdf_handle_error( 'lpm_data_output_particles', 19 )
! 
#endif

    CALL cpu_log( log_point_s(40), 'lpm_data_output', 'stop' )

 END SUBROUTINE lpm_data_output_particles
 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine calculates and provide particle timeseries output.
!------------------------------------------------------------------------------!
 SUBROUTINE lpm_data_output_ptseries
 
    INTEGER(iwp) ::  i    !<
    INTEGER(iwp) ::  inum !<
    INTEGER(iwp) ::  j    !<
    INTEGER(iwp) ::  jg   !<
    INTEGER(iwp) ::  k    !<
    INTEGER(iwp) ::  n    !<

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  pts_value   !<
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  pts_value_l !<


    CALL cpu_log( log_point(36), 'data_output_ptseries', 'start' )

    IF ( myid == 0 )  THEN
!
!--    Open file for time series output in NetCDF format
       dopts_time_count = dopts_time_count + 1
       CALL check_open( 109 )
#if defined( __netcdf )
!
!--    Update the particle time series time axis
       nc_stat = NF90_PUT_VAR( id_set_pts, id_var_time_pts,      &
                               (/ time_since_reference_point /), &
                               start = (/ dopts_time_count /), count = (/ 1 /) )
       CALL netcdf_handle_error( 'data_output_ptseries', 391 )
#endif

    ENDIF

    ALLOCATE( pts_value(0:number_of_particle_groups,dopts_num), &
              pts_value_l(0:number_of_particle_groups,dopts_num) )

    pts_value_l = 0.0_wp
    pts_value_l(:,16) = 9999999.9_wp    ! for calculation of minimum radius

!
!-- Calculate or collect the particle time series quantities for all particles
!-- and seperately for each particle group (if there is more than one group)
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb, nzt
             number_of_particles = prt_count(k,j,i)
             IF (number_of_particles <= 0)  CYCLE
             particles => grid_particles(k,j,i)%particles(1:number_of_particles)
             DO  n = 1, number_of_particles

                IF ( particles(n)%particle_mask )  THEN  ! Restrict analysis to active particles

                   pts_value_l(0,1)  = pts_value_l(0,1) + 1.0_wp  ! total # of particles
                   pts_value_l(0,2)  = pts_value_l(0,2) +                      &
                          ( particles(n)%x - particles(n)%origin_x )  ! mean x
                   pts_value_l(0,3)  = pts_value_l(0,3) +                      &
                          ( particles(n)%y - particles(n)%origin_y )  ! mean y
                   pts_value_l(0,4)  = pts_value_l(0,4) +                      &
                          ( particles(n)%z - particles(n)%origin_z )  ! mean z
                   pts_value_l(0,5)  = pts_value_l(0,5) + particles(n)%z        ! mean z (absolute)
                   pts_value_l(0,6)  = pts_value_l(0,6) + particles(n)%speed_x  ! mean u
                   pts_value_l(0,7)  = pts_value_l(0,7) + particles(n)%speed_y  ! mean v
                   pts_value_l(0,8)  = pts_value_l(0,8) + particles(n)%speed_z  ! mean w
                   pts_value_l(0,9)  = pts_value_l(0,9)  + particles(n)%rvar1 ! mean sgsu
                   pts_value_l(0,10) = pts_value_l(0,10) + particles(n)%rvar2 ! mean sgsv
                   pts_value_l(0,11) = pts_value_l(0,11) + particles(n)%rvar3 ! mean sgsw
                   IF ( particles(n)%speed_z > 0.0_wp )  THEN
                      pts_value_l(0,12) = pts_value_l(0,12) + 1.0_wp  ! # of upward moving prts
                      pts_value_l(0,13) = pts_value_l(0,13) +                  &
                                              particles(n)%speed_z ! mean w upw.
                   ELSE
                      pts_value_l(0,14) = pts_value_l(0,14) +                  &
                                              particles(n)%speed_z ! mean w down
                   ENDIF
                   pts_value_l(0,15) = pts_value_l(0,15) + particles(n)%radius ! mean rad
                   pts_value_l(0,16) = MIN( pts_value_l(0,16), particles(n)%radius ) ! minrad
                   pts_value_l(0,17) = MAX( pts_value_l(0,17), particles(n)%radius ) ! maxrad
                   pts_value_l(0,18) = pts_value_l(0,18) + 1.0_wp
                   pts_value_l(0,19) = pts_value_l(0,18) + 1.0_wp
!
!--                Repeat the same for the respective particle group
                   IF ( number_of_particle_groups > 1 )  THEN
                      jg = particles(n)%group

                      pts_value_l(jg,1)  = pts_value_l(jg,1) + 1.0_wp
                      pts_value_l(jg,2)  = pts_value_l(jg,2) +                   &
                           ( particles(n)%x - particles(n)%origin_x )
                      pts_value_l(jg,3)  = pts_value_l(jg,3) +                   &
                           ( particles(n)%y - particles(n)%origin_y )
                      pts_value_l(jg,4)  = pts_value_l(jg,4) +                   &
                           ( particles(n)%z - particles(n)%origin_z )
                      pts_value_l(jg,5)  = pts_value_l(jg,5) + particles(n)%z
                      pts_value_l(jg,6)  = pts_value_l(jg,6) + particles(n)%speed_x
                      pts_value_l(jg,7)  = pts_value_l(jg,7) + particles(n)%speed_y
                      pts_value_l(jg,8)  = pts_value_l(jg,8) + particles(n)%speed_z
                      pts_value_l(jg,9)  = pts_value_l(jg,9)  + particles(n)%rvar1
                      pts_value_l(jg,10) = pts_value_l(jg,10) + particles(n)%rvar2
                      pts_value_l(jg,11) = pts_value_l(jg,11) + particles(n)%rvar3
                      IF ( particles(n)%speed_z > 0.0_wp )  THEN
                         pts_value_l(jg,12) = pts_value_l(jg,12) + 1.0_wp
                         pts_value_l(jg,13) = pts_value_l(jg,13) + particles(n)%speed_z
                      ELSE
                         pts_value_l(jg,14) = pts_value_l(jg,14) + particles(n)%speed_z
                      ENDIF
                      pts_value_l(jg,15) = pts_value_l(jg,15) + particles(n)%radius
                      pts_value_l(jg,16) = MIN( pts_value_l(jg,16), particles(n)%radius )
                      pts_value_l(jg,17) = MAX( pts_value_l(jg,17), particles(n)%radius )
                      pts_value_l(jg,18) = pts_value_l(jg,18) + 1.0_wp
                      pts_value_l(jg,19) = pts_value_l(jg,19) + 1.0_wp
                   ENDIF

                ENDIF

             ENDDO

          ENDDO
       ENDDO
    ENDDO


#if defined( __parallel )
!
!-- Sum values of the subdomains
    inum = number_of_particle_groups + 1

    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( pts_value_l(0,1), pts_value(0,1), 15*inum, MPI_REAL, &
                        MPI_SUM, comm2d, ierr )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( pts_value_l(0,16), pts_value(0,16), inum, MPI_REAL, &
                        MPI_MIN, comm2d, ierr )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( pts_value_l(0,17), pts_value(0,17), inum, MPI_REAL, &
                        MPI_MAX, comm2d, ierr )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( pts_value_l(0,18), pts_value(0,18), inum, MPI_REAL, &
                        MPI_MAX, comm2d, ierr )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( pts_value_l(0,19), pts_value(0,19), inum, MPI_REAL, &
                        MPI_MIN, comm2d, ierr )
#else
    pts_value(:,1:19) = pts_value_l(:,1:19)
#endif

!
!-- Normalize the above calculated quantities (except min/max values) with the
!-- total number of particles
    IF ( number_of_particle_groups > 1 )  THEN
       inum = number_of_particle_groups
    ELSE
       inum = 0
    ENDIF

    DO  j = 0, inum

       IF ( pts_value(j,1) > 0.0_wp )  THEN

          pts_value(j,2:15) = pts_value(j,2:15) / pts_value(j,1)
          IF ( pts_value(j,12) > 0.0_wp  .AND.  pts_value(j,12) < 1.0_wp )  THEN
             pts_value(j,13) = pts_value(j,13) / pts_value(j,12)
             pts_value(j,14) = pts_value(j,14) / ( 1.0_wp - pts_value(j,12) )
          ELSEIF ( pts_value(j,12) == 0.0_wp )  THEN
             pts_value(j,13) = -1.0_wp
          ELSE
             pts_value(j,14) = -1.0_wp
          ENDIF

       ENDIF

    ENDDO

!
!-- Calculate higher order moments of particle time series quantities,
!-- seperately for each particle group (if there is more than one group)
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb, nzt
             number_of_particles = prt_count(k,j,i)
             IF (number_of_particles <= 0)  CYCLE
             particles => grid_particles(k,j,i)%particles(1:number_of_particles)
             DO  n = 1, number_of_particles

                pts_value_l(0,20) = pts_value_l(0,20) + ( particles(n)%x - &
                                    particles(n)%origin_x - pts_value(0,2) )**2 ! x*2
                pts_value_l(0,21) = pts_value_l(0,21) + ( particles(n)%y - &
                                    particles(n)%origin_y - pts_value(0,3) )**2 ! y*2
                pts_value_l(0,22) = pts_value_l(0,22) + ( particles(n)%z - &
                                    particles(n)%origin_z - pts_value(0,4) )**2 ! z*2
                pts_value_l(0,23) = pts_value_l(0,23) + ( particles(n)%speed_x - &
                                                         pts_value(0,6) )**2   ! u*2
                pts_value_l(0,24) = pts_value_l(0,24) + ( particles(n)%speed_y - &
                                                          pts_value(0,7) )**2   ! v*2
                pts_value_l(0,25) = pts_value_l(0,25) + ( particles(n)%speed_z - &
                                                          pts_value(0,8) )**2   ! w*2
                pts_value_l(0,26) = pts_value_l(0,26) + ( particles(n)%rvar1 - &
                                                          pts_value(0,9) )**2   ! u"2
                pts_value_l(0,27) = pts_value_l(0,27) + ( particles(n)%rvar2 - &
                                                          pts_value(0,10) )**2  ! v"2
                pts_value_l(0,28) = pts_value_l(0,28) + ( particles(n)%rvar3 - &
                                                          pts_value(0,11) )**2  ! w"2
!
!--             Repeat the same for the respective particle group
                IF ( number_of_particle_groups > 1 )  THEN
                   jg = particles(n)%group

                   pts_value_l(jg,20) = pts_value_l(jg,20) + ( particles(n)%x - &
                                       particles(n)%origin_x - pts_value(jg,2) )**2
                   pts_value_l(jg,21) = pts_value_l(jg,21) + ( particles(n)%y - &
                                       particles(n)%origin_y - pts_value(jg,3) )**2
                   pts_value_l(jg,22) = pts_value_l(jg,22) + ( particles(n)%z - &
                                       particles(n)%origin_z - pts_value(jg,4) )**2
                   pts_value_l(jg,23) = pts_value_l(jg,23) + ( particles(n)%speed_x - &
                                                             pts_value(jg,6) )**2
                   pts_value_l(jg,24) = pts_value_l(jg,24) + ( particles(n)%speed_y - &
                                                             pts_value(jg,7) )**2
                   pts_value_l(jg,25) = pts_value_l(jg,25) + ( particles(n)%speed_z - &
                                                             pts_value(jg,8) )**2
                   pts_value_l(jg,26) = pts_value_l(jg,26) + ( particles(n)%rvar1 - &
                                                             pts_value(jg,9) )**2
                   pts_value_l(jg,27) = pts_value_l(jg,27) + ( particles(n)%rvar2 - &
                                                             pts_value(jg,10) )**2
                   pts_value_l(jg,28) = pts_value_l(jg,28) + ( particles(n)%rvar3 - &
                                                             pts_value(jg,11) )**2
                ENDIF

             ENDDO
          ENDDO
       ENDDO
    ENDDO

    pts_value_l(0,29) = ( number_of_particles - pts_value(0,1) / numprocs )**2
                                                 ! variance of particle numbers
    IF ( number_of_particle_groups > 1 )  THEN
       DO  j = 1, number_of_particle_groups
          pts_value_l(j,29) = ( pts_value_l(j,1) - &
                                pts_value(j,1) / numprocs )**2
       ENDDO
    ENDIF

#if defined( __parallel )
!
!-- Sum values of the subdomains
    inum = number_of_particle_groups + 1

    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( pts_value_l(0,20), pts_value(0,20), inum*10, MPI_REAL, &
                        MPI_SUM, comm2d, ierr )
#else
    pts_value(:,20:29) = pts_value_l(:,20:29)
#endif

!
!-- Normalize the above calculated quantities with the total number of
!-- particles
    IF ( number_of_particle_groups > 1 )  THEN
       inum = number_of_particle_groups
    ELSE
       inum = 0
    ENDIF

    DO  j = 0, inum

       IF ( pts_value(j,1) > 0.0_wp )  THEN
          pts_value(j,20:28) = pts_value(j,20:28) / pts_value(j,1)
       ENDIF
       pts_value(j,29) = pts_value(j,29) / numprocs

    ENDDO

#if defined( __netcdf )
!
!-- Output particle time series quantities in NetCDF format
    IF ( myid == 0 )  THEN
       DO  j = 0, inum
          DO  i = 1, dopts_num
             nc_stat = NF90_PUT_VAR( id_set_pts, id_var_dopts(i,j),  &
                                     (/ pts_value(j,i) /),           &
                                     start = (/ dopts_time_count /), &
                                     count = (/ 1 /) )
             CALL netcdf_handle_error( 'data_output_ptseries', 392 )
          ENDDO
       ENDDO
    ENDIF
#endif

    DEALLOCATE( pts_value, pts_value_l )

    CALL cpu_log( log_point(36), 'data_output_ptseries', 'stop' )

END SUBROUTINE lpm_data_output_ptseries

 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine reads the respective restart data for the lpm.
!------------------------------------------------------------------------------!
 SUBROUTINE lpm_rrd_local_particles

    CHARACTER (LEN=10) ::  particle_binary_version    !<
    CHARACTER (LEN=10) ::  version_on_file            !<

    INTEGER(iwp) ::  alloc_size !<
    INTEGER(iwp) ::  ip         !<
    INTEGER(iwp) ::  jp         !<
    INTEGER(iwp) ::  kp         !<

    TYPE(particle_type), DIMENSION(:), ALLOCATABLE ::  tmp_particles !<

!
!-- Read particle data from previous model run.
!-- First open the input unit.
    IF ( myid_char == '' )  THEN
       OPEN ( 90, FILE='PARTICLE_RESTART_DATA_IN'//myid_char,                  &
                  FORM='UNFORMATTED' )
    ELSE
       OPEN ( 90, FILE='PARTICLE_RESTART_DATA_IN/'//myid_char,                 &
                  FORM='UNFORMATTED' )
    ENDIF

!
!-- First compare the version numbers
    READ ( 90 )  version_on_file
    particle_binary_version = '4.0'
    IF ( TRIM( version_on_file ) /= TRIM( particle_binary_version ) )  THEN
       message_string = 'version mismatch concerning data from prior ' //      &
                        'run &version on file = "' //                          &
                                      TRIM( version_on_file ) //               &
                        '&version in program = "' //                           &
                                      TRIM( particle_binary_version ) // '"'
       CALL message( 'lpm_read_restart_file', 'PA0214', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- If less particles are stored on the restart file than prescribed by
!-- 1, the remainder is initialized by zero_particle to avoid
!-- errors.
    zero_particle = particle_type( 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,     &
                                   0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,     &
                                   0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,     &
                                   0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,     &
                                   0, 0, 0_idp, .FALSE., -1 )
!
!-- Read some particle parameters and the size of the particle arrays,
!-- allocate them and read their contents.
    READ ( 90 )  bc_par_b, bc_par_lr, bc_par_ns, bc_par_t,                     &
                 last_particle_release_time, number_of_particle_groups,        &
                 particle_groups, time_write_particle_data

    ALLOCATE( prt_count(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                        &
              grid_particles(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

    READ ( 90 )  prt_count

    DO  ip = nxl, nxr
       DO  jp = nys, nyn
          DO  kp = nzb+1, nzt

             number_of_particles = prt_count(kp,jp,ip)
             IF ( number_of_particles > 0 )  THEN
                alloc_size = MAX( INT( number_of_particles *                   &
                             ( 1.0_wp + alloc_factor / 100.0_wp ) ),           &
                             1 )
             ELSE
                alloc_size = 1
             ENDIF

             ALLOCATE( grid_particles(kp,jp,ip)%particles(1:alloc_size) )

             IF ( number_of_particles > 0 )  THEN
                ALLOCATE( tmp_particles(1:number_of_particles) )
                READ ( 90 )  tmp_particles
                grid_particles(kp,jp,ip)%particles(1:number_of_particles) = tmp_particles
                DEALLOCATE( tmp_particles )
                IF ( number_of_particles < alloc_size )  THEN
                   grid_particles(kp,jp,ip)%particles(number_of_particles+1:alloc_size) &
                      = zero_particle
                ENDIF
             ELSE
                grid_particles(kp,jp,ip)%particles(1:alloc_size) = zero_particle
             ENDIF

          ENDDO
       ENDDO
    ENDDO

    CLOSE ( 90 )
!
!-- Must be called to sort particles into blocks, which is needed for a fast
!-- interpolation of the LES fields on the particle position.
    CALL lpm_sort_and_delete


 END SUBROUTINE lpm_rrd_local_particles
 
 
 SUBROUTINE lpm_rrd_local( k, nxlf, nxlc, nxl_on_file, nxrf, nxrc,          &
                              nxr_on_file, nynf, nync, nyn_on_file, nysf,  &
                              nysc, nys_on_file, tmp_3d, found )


   USE control_parameters,                                                 &
       ONLY: length, restart_string

    INTEGER(iwp) ::  k               !<
    INTEGER(iwp) ::  nxlc            !<
    INTEGER(iwp) ::  nxlf            !<
    INTEGER(iwp) ::  nxl_on_file     !<
    INTEGER(iwp) ::  nxrc            !<
    INTEGER(iwp) ::  nxrf            !<
    INTEGER(iwp) ::  nxr_on_file     !<
    INTEGER(iwp) ::  nync            !<
    INTEGER(iwp) ::  nynf            !<
    INTEGER(iwp) ::  nyn_on_file     !<
    INTEGER(iwp) ::  nysc            !<
    INTEGER(iwp) ::  nysf            !<
    INTEGER(iwp) ::  nys_on_file     !<

    LOGICAL, INTENT(OUT)  ::  found

    REAL(wp), DIMENSION(nzb:nzt+1,nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) ::  tmp_3d   !<


    found = .TRUE.

    SELECT CASE ( restart_string(1:length) )

       CASE ( 'iran' ) ! matching random numbers is still unresolved issue
          IF ( k == 1 )  READ ( 13 )  iran, iran_part

        CASE ( 'pc_av' )
           IF ( .NOT. ALLOCATED( pc_av ) )  THEN
              ALLOCATE( pc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
           ENDIF
           IF ( k == 1 )  READ ( 13 )  tmp_3d
           pc_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =          &
              tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

        CASE ( 'pr_av' )
           IF ( .NOT. ALLOCATED( pr_av ) )  THEN
              ALLOCATE( pr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
           ENDIF
           IF ( k == 1 )  READ ( 13 )  tmp_3d
           pr_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =          &
              tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
 
         CASE ( 'ql_c_av' )
            IF ( .NOT. ALLOCATED( ql_c_av ) )  THEN
               ALLOCATE( ql_c_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
            ENDIF
            IF ( k == 1 )  READ ( 13 )  tmp_3d
            ql_c_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =        &
               tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

         CASE ( 'ql_v_av' )
            IF ( .NOT. ALLOCATED( ql_v_av ) )  THEN
               ALLOCATE( ql_v_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
            ENDIF
            IF ( k == 1 )  READ ( 13 )  tmp_3d
            ql_v_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =        &
               tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

         CASE ( 'ql_vp_av' )
            IF ( .NOT. ALLOCATED( ql_vp_av ) )  THEN
               ALLOCATE( ql_vp_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
            ENDIF
            IF ( k == 1 )  READ ( 13 )  tmp_3d
            ql_vp_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =       &
               tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE DEFAULT

             found = .FALSE.

       END SELECT


 END SUBROUTINE lpm_rrd_local
 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine writes the respective restart data for the lpm.
!------------------------------------------------------------------------------!
 SUBROUTINE lpm_wrd_local
 
    CHARACTER (LEN=10) ::  particle_binary_version   !<

    INTEGER(iwp) ::  ip                              !<
    INTEGER(iwp) ::  jp                              !<
    INTEGER(iwp) ::  kp                              !< 
!
!-- First open the output unit.
    IF ( myid_char == '' )  THEN
       OPEN ( 90, FILE='PARTICLE_RESTART_DATA_OUT'//myid_char, &
                  FORM='UNFORMATTED')
    ELSE
       IF ( myid == 0 )  CALL local_system( 'mkdir PARTICLE_RESTART_DATA_OUT' )
#if defined( __parallel )
!
!--    Set a barrier in order to allow that thereafter all other processors
!--    in the directory created by PE0 can open their file
       CALL MPI_BARRIER( comm2d, ierr )
#endif
       OPEN ( 90, FILE='PARTICLE_RESTART_DATA_OUT/'//myid_char, &
                  FORM='UNFORMATTED' )
    ENDIF

!
!-- Write the version number of the binary format.
!-- Attention: After changes to the following output commands the version
!-- ---------  number of the variable particle_binary_version must be
!--            changed! Also, the version number and the list of arrays
!--            to be read in lpm_read_restart_file must be adjusted
!--            accordingly.
    particle_binary_version = '4.0'
    WRITE ( 90 )  particle_binary_version

!
!-- Write some particle parameters, the size of the particle arrays
    WRITE ( 90 )  bc_par_b, bc_par_lr, bc_par_ns, bc_par_t,                    &
                  last_particle_release_time, number_of_particle_groups,       &
                  particle_groups, time_write_particle_data

    WRITE ( 90 )  prt_count
          
    DO  ip = nxl, nxr
       DO  jp = nys, nyn
          DO  kp = nzb+1, nzt
             number_of_particles = prt_count(kp,jp,ip)
             particles => grid_particles(kp,jp,ip)%particles(1:number_of_particles)
             IF ( number_of_particles <= 0 )  CYCLE
             WRITE ( 90 )  particles
          ENDDO
       ENDDO
    ENDDO

    CLOSE ( 90 )

#if defined( __parallel )
       CALL MPI_BARRIER( comm2d, ierr )
#endif

    CALL wrd_write_string( 'iran' ) 
    WRITE ( 14 )  iran, iran_part 


 END SUBROUTINE lpm_wrd_local


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine writes the respective restart data for the lpm.
!------------------------------------------------------------------------------! 
 SUBROUTINE lpm_wrd_global
 
    CALL wrd_write_string( 'curvature_solution_effects' ) 
    WRITE ( 14 )  curvature_solution_effects

    CALL wrd_write_string( 'interpolation_simple_corrector' )
    WRITE ( 14 )  interpolation_simple_corrector

    CALL wrd_write_string( 'interpolation_simple_predictor' )
    WRITE ( 14 )  interpolation_simple_predictor

    CALL wrd_write_string( 'interpolation_trilinear' )
    WRITE ( 14 )  interpolation_trilinear

 END SUBROUTINE lpm_wrd_global
 

!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine writes the respective restart data for the lpm.
!------------------------------------------------------------------------------! 
 SUBROUTINE lpm_rrd_global( found )
 
    USE control_parameters,                            &
        ONLY: length, restart_string

    LOGICAL, INTENT(OUT)  ::  found

    found = .TRUE.

    SELECT CASE ( restart_string(1:length) )

       CASE ( 'curvature_solution_effects' )
          READ ( 13 )  curvature_solution_effects

       CASE ( 'interpolation_simple_corrector' )
          READ ( 13 )  interpolation_simple_corrector

       CASE ( 'interpolation_simple_predictor' )
          READ ( 13 )  interpolation_simple_predictor

       CASE ( 'interpolation_trilinear' )
          READ ( 13 )  interpolation_trilinear

!          CASE ( 'global_paramter' )
!             READ ( 13 )  global_parameter
!          CASE ( 'global_array' )
!             IF ( .NOT. ALLOCATED( global_array ) )  ALLOCATE( global_array(1:10) )
!             READ ( 13 )  global_array

       CASE DEFAULT

          found = .FALSE.

    END SELECT 
    
 END SUBROUTINE lpm_rrd_global


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This is a submodule of the lagrangian particle model. It contains all 
!> dynamic processes of the lpm. This includes the advection (resolved and sub-
!> grid scale) as well as the boundary conditions of particles. As a next step
!> this submodule should be excluded as an own file. 
!------------------------------------------------------------------------------!
 SUBROUTINE lpm_advec (ip,jp,kp)

    LOGICAL ::  subbox_at_wall !< flag to see if the current subgridbox is adjacent to a wall

    INTEGER(iwp) ::  i                           !< index variable along x
    INTEGER(iwp) ::  i_next                      !< index variable along x
    INTEGER(iwp) ::  ip                          !< index variable along x
    INTEGER(iwp) ::  iteration_steps = 1         !< amount of iterations steps for corrector step
    INTEGER(iwp) ::  j                           !< index variable along y
    INTEGER(iwp) ::  j_next                      !< index variable along y
    INTEGER(iwp) ::  jp                          !< index variable along y
    INTEGER(iwp) ::  k                           !< index variable along z
    INTEGER(iwp) ::  k_wall                      !< vertical index of topography top 
    INTEGER(iwp) ::  kp                          !< index variable along z
    INTEGER(iwp) ::  k_next                      !< index variable along z
    INTEGER(iwp) ::  kw                          !< index variable along z
    INTEGER(iwp) ::  kkw                         !< index variable along z
    INTEGER(iwp) ::  n                           !< loop variable over all particles in a grid box
    INTEGER(iwp) ::  nb                          !< block number particles are sorted in
    INTEGER(iwp) ::  particle_end                !< end index for partilce loop
    INTEGER(iwp) ::  particle_start              !< start index for particle loop
    INTEGER(iwp) ::  surf_start                  !< Index on surface data-type for current grid box
    INTEGER(iwp) ::  subbox_end                  !< end index for loop over subboxes in particle advection
    INTEGER(iwp) ::  subbox_start                !< start index for loop over subboxes in particle advection
    INTEGER(iwp) ::  nn                          !< loop variable over iterations steps

    INTEGER(iwp), DIMENSION(0:7) ::  start_index !< start particle index for current block
    INTEGER(iwp), DIMENSION(0:7) ::  end_index   !< start particle index for current block

    REAL(wp) ::  aa                 !< dummy argument for horizontal particle interpolation
    REAL(wp) ::  alpha              !< interpolation facor for x-direction

    REAL(wp) ::  bb                 !< dummy argument for horizontal particle interpolation
    REAL(wp) ::  beta               !< interpolation facor for y-direction
    REAL(wp) ::  cc                 !< dummy argument for horizontal particle interpolation
    REAL(wp) ::  d_z_p_z0           !< inverse of interpolation length for logarithmic interpolation 
    REAL(wp) ::  dd                 !< dummy argument for horizontal particle interpolation 
    REAL(wp) ::  de_dx_int_l        !< x/y-interpolated TKE gradient (x) at particle position at lower vertical level
    REAL(wp) ::  de_dx_int_u        !< x/y-interpolated TKE gradient (x) at particle position at upper vertical level
    REAL(wp) ::  de_dy_int_l        !< x/y-interpolated TKE gradient (y) at particle position at lower vertical level
    REAL(wp) ::  de_dy_int_u        !< x/y-interpolated TKE gradient (y) at particle position at upper vertical level
    REAL(wp) ::  de_dt              !< temporal derivative of TKE experienced by the particle
    REAL(wp) ::  de_dt_min          !< lower level for temporal TKE derivative
    REAL(wp) ::  de_dz_int_l        !< x/y-interpolated TKE gradient (z) at particle position at lower vertical level
    REAL(wp) ::  de_dz_int_u        !< x/y-interpolated TKE gradient (z) at particle position at upper vertical level
    REAL(wp) ::  diameter           !< diamter of droplet
    REAL(wp) ::  diss_int_l         !< x/y-interpolated dissipation at particle position at lower vertical level
    REAL(wp) ::  diss_int_u         !< x/y-interpolated dissipation at particle position at upper vertical level
    REAL(wp) ::  dt_particle_m      !< previous particle time step
    REAL(wp) ::  dz_temp            !< dummy for the vertical grid spacing
    REAL(wp) ::  e_int_l            !< x/y-interpolated TKE at particle position at lower vertical level
    REAL(wp) ::  e_int_u            !< x/y-interpolated TKE at particle position at upper vertical level
    REAL(wp) ::  e_mean_int         !< horizontal mean TKE at particle height
    REAL(wp) ::  exp_arg            !< argument in the exponent - particle radius
    REAL(wp) ::  exp_term           !< exponent term
    REAL(wp) ::  gamma              !< interpolation facor for z-direction
    REAL(wp) ::  gg                 !< dummy argument for horizontal particle interpolation 
    REAL(wp) ::  height_p           !< dummy argument for logarithmic interpolation
    REAL(wp) ::  log_z_z0_int       !< logarithmus used for surface_layer interpolation
    REAL(wp) ::  random_gauss       !< Gaussian-distributed random number used for SGS particle advection
    REAL(wp) ::  RL                 !< Lagrangian autocorrelation coefficient
    REAL(wp) ::  rg1                !< Gaussian distributed random number
    REAL(wp) ::  rg2                !< Gaussian distributed random number
    REAL(wp) ::  rg3                !< Gaussian distributed random number
    REAL(wp) ::  sigma              !< velocity standard deviation
    REAL(wp) ::  u_int_l            !< x/y-interpolated u-component at particle position at lower vertical level
    REAL(wp) ::  u_int_u            !< x/y-interpolated u-component at particle position at upper vertical level
    REAL(wp) ::  unext              !< calculated particle u-velocity of corrector step
    REAL(wp) ::  us_int             !< friction velocity at particle grid box
    REAL(wp) ::  usws_int           !< surface momentum flux (u component) at particle grid box
    REAL(wp) ::  v_int_l            !< x/y-interpolated v-component at particle position at lower vertical level
    REAL(wp) ::  v_int_u            !< x/y-interpolated v-component at particle position at upper vertical level
    REAL(wp) ::  vsws_int           !< surface momentum flux (u component) at particle grid box
    REAL(wp) ::  vnext              !< calculated particle v-velocity of corrector step
    REAL(wp) ::  vv_int             !< dummy to compute interpolated mean SGS TKE, used to scale SGS advection 
    REAL(wp) ::  w_int_l            !< x/y-interpolated w-component at particle position at lower vertical level
    REAL(wp) ::  w_int_u            !< x/y-interpolated w-component at particle position at upper vertical level
    REAL(wp) ::  wnext              !< calculated particle w-velocity of corrector step
    REAL(wp) ::  w_s                !< terminal velocity of droplets
    REAL(wp) ::  x                  !< dummy argument for horizontal particle interpolation 
    REAL(wp) ::  xp                 !< calculated particle position in x of predictor step
    REAL(wp) ::  y                  !< dummy argument for horizontal particle interpolation
    REAL(wp) ::  yp                 !< calculated particle position in y of predictor step
    REAL(wp) ::  z_p                !< surface layer height (0.5 dz)
    REAL(wp) ::  zp                 !< calculated particle position in z of predictor step

    REAL(wp), PARAMETER ::  a_rog = 9.65_wp      !< parameter for fall velocity
    REAL(wp), PARAMETER ::  b_rog = 10.43_wp     !< parameter for fall velocity
    REAL(wp), PARAMETER ::  c_rog = 0.6_wp       !< parameter for fall velocity
    REAL(wp), PARAMETER ::  k_cap_rog = 4.0_wp   !< parameter for fall velocity
    REAL(wp), PARAMETER ::  k_low_rog = 12.0_wp  !< parameter for fall velocity
    REAL(wp), PARAMETER ::  d0_rog = 0.745_wp    !< separation diameter

    REAL(wp), DIMENSION(number_of_particles) ::  term_1_2       !< flag to communicate whether a particle is near topography or not
    REAL(wp), DIMENSION(number_of_particles) ::  dens_ratio     !< ratio between the density of the fluid and the density of the particles 
    REAL(wp), DIMENSION(number_of_particles) ::  de_dx_int      !< horizontal TKE gradient along x at particle position
    REAL(wp), DIMENSION(number_of_particles) ::  de_dy_int      !< horizontal TKE gradient along y at particle position
    REAL(wp), DIMENSION(number_of_particles) ::  de_dz_int      !< horizontal TKE gradient along z at particle position
    REAL(wp), DIMENSION(number_of_particles) ::  diss_int       !< dissipation at particle position
    REAL(wp), DIMENSION(number_of_particles) ::  dt_gap         !< remaining time until particle time integration reaches LES time
    REAL(wp), DIMENSION(number_of_particles) ::  dt_particle    !< particle time step
    REAL(wp), DIMENSION(number_of_particles) ::  e_int          !< TKE at particle position
    REAL(wp), DIMENSION(number_of_particles) ::  fs_int         !< weighting factor for subgrid-scale particle speed
    REAL(wp), DIMENSION(number_of_particles) ::  lagr_timescale !< Lagrangian timescale
    REAL(wp), DIMENSION(number_of_particles) ::  rvar1_temp     !< SGS particle velocity - u-component
    REAL(wp), DIMENSION(number_of_particles) ::  rvar2_temp     !< SGS particle velocity - v-component
    REAL(wp), DIMENSION(number_of_particles) ::  rvar3_temp     !< SGS particle velocity - w-component
    REAL(wp), DIMENSION(number_of_particles) ::  u_int          !< u-component of particle speed
    REAL(wp), DIMENSION(number_of_particles) ::  v_int          !< v-component of particle speed 
    REAL(wp), DIMENSION(number_of_particles) ::  w_int          !< w-component of particle speed
    REAL(wp), DIMENSION(number_of_particles) ::  xv             !< x-position
    REAL(wp), DIMENSION(number_of_particles) ::  yv             !< y-position
    REAL(wp), DIMENSION(number_of_particles) ::  zv             !< z-position

    REAL(wp), DIMENSION(number_of_particles, 3) ::  rg          !< vector of Gaussian distributed random numbers

    CALL cpu_log( log_point_s(44), 'lpm_advec', 'continue' )
!
!-- Determine height of Prandtl layer and distance between Prandtl-layer
!-- height and horizontal mean roughness height, which are required for
!-- vertical logarithmic interpolation of horizontal particle speeds
!-- (for particles below first vertical grid level).
    z_p      = zu(nzb+1) - zw(nzb)
    d_z_p_z0 = 1.0_wp / ( z_p - z0_av_global )

    xv = particles(1:number_of_particles)%x
    yv = particles(1:number_of_particles)%y
    zv = particles(1:number_of_particles)%z
    dt_particle = dt_3d

!
!-- This case uses a simple interpolation method for the particle velocites,
!-- and applying a predictor-corrector method. @note the current time divergence
!-- free time step is denoted with u_t etc.; the velocities of the time level of
!-- t+1 wit u,v, and w, as the model is called after swap timelevel
!-- @attention: for the corrector step the velocities of t(n+1) are required.
!-- Therefore the particle code is executed at the end of the time intermediate
!-- timestep routine. This interpolation method is described in more detail
!-- in Grabowski et al., 2018 (GMD).
    IF ( interpolation_simple_corrector )  THEN
!
!--    Predictor step
       kkw = kp - 1
       DO  n = 1, number_of_particles

          alpha = MAX( MIN( ( particles(n)%x - ip * dx ) * ddx, 1.0_wp ), 0.0_wp )
          u_int(n) = u_t(kp,jp,ip) * ( 1.0_wp - alpha ) + u_t(kp,jp,ip+1) * alpha

          beta  = MAX( MIN( ( particles(n)%y - jp * dy ) * ddy, 1.0_wp ), 0.0_wp )
          v_int(n) = v_t(kp,jp,ip) * ( 1.0_wp - beta ) + v_t(kp,jp+1,ip) * beta

          gamma = MAX( MIN( ( particles(n)%z - zw(kkw) ) /                   &
                            ( zw(kkw+1) - zw(kkw) ), 1.0_wp ), 0.0_wp )
          w_int(n) = w_t(kkw,jp,ip) * ( 1.0_wp - gamma ) + w_t(kkw+1,jp,ip) * gamma

       ENDDO
!
!--    Corrector step
       DO  n = 1, number_of_particles

          IF ( .NOT. particles(n)%particle_mask )  CYCLE

          DO  nn = 1, iteration_steps

!
!--          Guess new position
             xp = particles(n)%x + u_int(n) * dt_particle(n)
             yp = particles(n)%y + v_int(n) * dt_particle(n)
             zp = particles(n)%z + w_int(n) * dt_particle(n)
!
!--          x direction
             i_next = FLOOR( xp * ddx , KIND=iwp)
             alpha  = MAX( MIN( ( xp - i_next * dx ) * ddx, 1.0_wp ), 0.0_wp )
!
!--          y direction
             j_next = FLOOR( yp * ddy )
             beta   = MAX( MIN( ( yp - j_next * dy ) * ddy, 1.0_wp ), 0.0_wp )
!
!--          z_direction
             k_next = MAX( MIN( FLOOR( zp / (zw(kkw+1)-zw(kkw)) + offset_ocean_nzt ), nzt ), 0)
             gamma = MAX( MIN( ( zp - zw(k_next) ) /                      &
                               ( zw(k_next+1) - zw(k_next) ), 1.0_wp ), 0.0_wp )
!
!--          Calculate part of the corrector step
             unext = u(k_next+1, j_next, i_next) * ( 1.0_wp - alpha ) +    &
                     u(k_next+1, j_next,   i_next+1) * alpha

             vnext = v(k_next+1, j_next, i_next) * ( 1.0_wp - beta  ) +    &
                     v(k_next+1, j_next+1, i_next  ) * beta

             wnext = w(k_next,   j_next, i_next) * ( 1.0_wp - gamma ) +    &
                     w(k_next+1, j_next, i_next  ) * gamma

!
!--          Calculate interpolated particle velocity with predictor
!--          corrector step. u_int, v_int and w_int describes the part of
!--          the predictor step. unext, vnext and wnext is the part of the
!--          corrector step. The resulting new position is set below. The
!--          implementation is based on Grabowski et al., 2018 (GMD).
             u_int(n) = 0.5_wp * ( u_int(n) + unext )
             v_int(n) = 0.5_wp * ( v_int(n) + vnext )
             w_int(n) = 0.5_wp * ( w_int(n) + wnext )

          ENDDO
       ENDDO
!
!-- This case uses a simple interpolation method for the particle velocites,
!-- and applying a predictor.
    ELSEIF ( interpolation_simple_predictor )  THEN
!
!--    The particle position for the w velociy is based on the value of kp and kp-1
       kkw = kp - 1
       DO  n = 1, number_of_particles
          IF ( .NOT. particles(n)%particle_mask )  CYCLE

          alpha    = MAX( MIN( ( particles(n)%x - ip * dx ) * ddx, 1.0_wp ), 0.0_wp )
          u_int(n) = u(kp,jp,ip) * ( 1.0_wp - alpha ) + u(kp,jp,ip+1) * alpha

          beta     = MAX( MIN( ( particles(n)%y - jp * dy ) * ddy, 1.0_wp ), 0.0_wp )
          v_int(n) = v(kp,jp,ip) * ( 1.0_wp - beta ) + v(kp,jp+1,ip) * beta

          gamma    = MAX( MIN( ( particles(n)%z - zw(kkw) ) /                   &
                               ( zw(kkw+1) - zw(kkw) ), 1.0_wp ), 0.0_wp )
          w_int(n) = w(kkw,jp,ip) * ( 1.0_wp - gamma ) + w(kkw+1,jp,ip) * gamma
       ENDDO
!
!-- The trilinear interpolation.
    ELSEIF ( interpolation_trilinear )  THEN

       start_index = grid_particles(kp,jp,ip)%start_index
       end_index   = grid_particles(kp,jp,ip)%end_index

       DO  nb = 0, 7
!
!--       Interpolate u velocity-component
          i = ip
          j = jp + block_offset(nb)%j_off
          k = kp + block_offset(nb)%k_off

          DO  n = start_index(nb), end_index(nb)
!
!--          Interpolation of the u velocity component onto particle position.
!--          Particles are interpolation bi-linearly in the horizontal and a
!--          linearly in the vertical. An exception is made for particles below
!--          the first vertical grid level in case of a prandtl layer. In this
!--          case the horizontal particle velocity components are determined using
!--          Monin-Obukhov relations (if branch).
!--          First, check if particle is located below first vertical grid level
!--          above topography (Prandtl-layer height)
!--          Determine vertical index of topography top
             k_wall = topo_top_ind(jp,ip,0)

             IF ( constant_flux_layer  .AND.  zv(n) - zw(k_wall) < z_p )  THEN
!
!--             Resolved-scale horizontal particle velocity is zero below z0.
                IF ( zv(n) - zw(k_wall) < z0_av_global )  THEN
                   u_int(n) = 0.0_wp
                ELSE
!
!--                Determine the sublayer. Further used as index.
                   height_p = ( zv(n) - zw(k_wall) - z0_av_global )            &
                                        * REAL( number_of_sublayers, KIND=wp ) &
                                        * d_z_p_z0
!
!--                Calculate LOG(z/z0) for exact particle height. Therefore,
!--                interpolate linearly between precalculated logarithm.
                   log_z_z0_int = log_z_z0(INT(height_p))                      &
                                    + ( height_p - INT(height_p) )             &
                                    * ( log_z_z0(INT(height_p)+1)              &
                                         - log_z_z0(INT(height_p))             &
                                      )
!
!--                Get friction velocity and momentum flux from new surface data
!--                types.
                   IF ( surf_def_h(0)%start_index(jp,ip) <=                    &
                        surf_def_h(0)%end_index(jp,ip) )  THEN
                      surf_start = surf_def_h(0)%start_index(jp,ip)
!--                   Limit friction velocity. In narrow canyons or holes the
!--                   friction velocity can become very small, resulting in a too
!--                   large particle speed.
                      us_int    = MAX( surf_def_h(0)%us(surf_start), 0.01_wp )
                      usws_int  = surf_def_h(0)%usws(surf_start)               &
                                * drho_air_zw(k_wall)
                   ELSEIF ( surf_lsm_h%start_index(jp,ip) <=                   &
                            surf_lsm_h%end_index(jp,ip) )  THEN
                      surf_start = surf_lsm_h%start_index(jp,ip)
                      us_int    = MAX( surf_lsm_h%us(surf_start), 0.01_wp )
                      usws_int  = surf_lsm_h%usws(surf_start)                  &
                                * drho_air_zw(k_wall)
                   ELSEIF ( surf_usm_h%start_index(jp,ip) <=                   &
                            surf_usm_h%end_index(jp,ip) )  THEN
                      surf_start = surf_usm_h%start_index(jp,ip)
                      us_int    = MAX( surf_usm_h%us(surf_start), 0.01_wp )
                      usws_int  = surf_usm_h%usws(surf_start)                  &
                                * drho_air_zw(k_wall)
                   ENDIF
!
!--                Neutral solution is applied for all situations, e.g. also for
!--                unstable and stable situations. Even though this is not exact
!--                this saves a lot of CPU time since several calls of intrinsic
!--                FORTRAN procedures (LOG, ATAN) are avoided, This is justified
!--                as sensitivity studies revealed no significant effect of
!--                using the neutral solution also for un/stable situations.
                   u_int(n) = -usws_int / ( us_int * kappa + 1E-10_wp )        &
                               * log_z_z0_int - u_gtrans
                ENDIF
!
!--          Particle above the first grid level. Bi-linear interpolation in the
!--          horizontal and linear interpolation in the vertical direction.
             ELSE
                x  = xv(n) - i * dx
                y  = yv(n) + ( 0.5_wp - j ) * dy
                aa = x**2          + y**2
                bb = ( dx - x )**2 + y**2
                cc = x**2          + ( dy - y )**2
                dd = ( dx - x )**2 + ( dy - y )**2
                gg = aa + bb + cc + dd

                u_int_l = ( ( gg - aa ) * u(k,j,i)   + ( gg - bb ) * u(k,j,i+1)   &
                            + ( gg - cc ) * u(k,j+1,i) + ( gg - dd ) *            &
                            u(k,j+1,i+1) ) / ( 3.0_wp * gg ) - u_gtrans

                IF ( k == nzt )  THEN
                   u_int(n) = u_int_l
                ELSE
                   u_int_u = ( ( gg-aa ) * u(k+1,j,i) + ( gg-bb ) * u(k+1,j,i+1)  &
                               + ( gg-cc ) * u(k+1,j+1,i) + ( gg-dd ) *           &
                               u(k+1,j+1,i+1) ) / ( 3.0_wp * gg ) - u_gtrans
                   u_int(n) = u_int_l + ( zv(n) - zu(k) ) / dzw(k+1) *            &
                              ( u_int_u - u_int_l )
                ENDIF
             ENDIF
          ENDDO
!
!--       Same procedure for interpolation of the v velocity-component
          i = ip + block_offset(nb)%i_off
          j = jp
          k = kp + block_offset(nb)%k_off

          DO  n = start_index(nb), end_index(nb)
!
!--          Determine vertical index of topography top
             k_wall = topo_top_ind(jp,ip,0)

             IF ( constant_flux_layer  .AND.  zv(n) - zw(k_wall) < z_p )  THEN
                IF ( zv(n) - zw(k_wall) < z0_av_global )  THEN
!
!--                Resolved-scale horizontal particle velocity is zero below z0.
                   v_int(n) = 0.0_wp
                ELSE
!
!--                Determine the sublayer. Further used as index. Please note,
!--                logarithmus can not be reused from above, as in in case of
!--                topography particle on u-grid can be above surface-layer height,
!--                whereas it can be below on v-grid.
                   height_p = ( zv(n) - zw(k_wall) - z0_av_global )            &
                                     * REAL( number_of_sublayers, KIND=wp )    &
                                     * d_z_p_z0
!
!--                Calculate LOG(z/z0) for exact particle height. Therefore,
!--                interpolate linearly between precalculated logarithm.
                   log_z_z0_int = log_z_z0(INT(height_p))                      &
                                    + ( height_p - INT(height_p) )             &
                                    * ( log_z_z0(INT(height_p)+1)              &
                                         - log_z_z0(INT(height_p))             &
                                      )
!
!--                Get friction velocity and momentum flux from new surface data
!--                types.
                   IF ( surf_def_h(0)%start_index(jp,ip) <=                    &
                        surf_def_h(0)%end_index(jp,ip) )  THEN
                      surf_start = surf_def_h(0)%start_index(jp,ip)
!--                   Limit friction velocity. In narrow canyons or holes the
!--                   friction velocity can become very small, resulting in a too
!--                   large particle speed.
                      us_int    = MAX( surf_def_h(0)%us(surf_start), 0.01_wp )
                      vsws_int  = surf_def_h(0)%vsws(surf_start)               &
                                * drho_air_zw(k_wall)
                   ELSEIF ( surf_lsm_h%start_index(jp,ip) <=                   &
                            surf_lsm_h%end_index(jp,ip) )  THEN
                      surf_start = surf_lsm_h%start_index(jp,ip)
                      us_int    = MAX( surf_lsm_h%us(surf_start), 0.01_wp )
                      vsws_int  = surf_lsm_h%vsws(surf_start)                  &
                                * drho_air_zw(k_wall)
                   ELSEIF ( surf_usm_h%start_index(jp,ip) <=                   &
                            surf_usm_h%end_index(jp,ip) )  THEN
                      surf_start = surf_usm_h%start_index(jp,ip)
                      us_int    = MAX( surf_usm_h%us(surf_start), 0.01_wp )
                      vsws_int  = surf_usm_h%vsws(surf_start)                  &
                                * drho_air_zw(k_wall)
                   ENDIF
!
!--                Neutral solution is applied for all situations, e.g. also for
!--                unstable and stable situations. Even though this is not exact
!--                this saves a lot of CPU time since several calls of intrinsic
!--                FORTRAN procedures (LOG, ATAN) are avoided, This is justified
!--                as sensitivity studies revealed no significant effect of
!--                using the neutral solution also for un/stable situations.
                   v_int(n) = -vsws_int / ( us_int * kappa + 1E-10_wp )        &
                            * log_z_z0_int - v_gtrans

                ENDIF
             ELSE
                x  = xv(n) + ( 0.5_wp - i ) * dx
                y  = yv(n) - j * dy
                aa = x**2          + y**2
                bb = ( dx - x )**2 + y**2
                cc = x**2          + ( dy - y )**2
                dd = ( dx - x )**2 + ( dy - y )**2
                gg = aa + bb + cc + dd

                v_int_l = ( ( gg - aa ) * v(k,j,i)   + ( gg - bb ) * v(k,j,i+1)   &
                          + ( gg - cc ) * v(k,j+1,i) + ( gg - dd ) * v(k,j+1,i+1) &
                          ) / ( 3.0_wp * gg ) - v_gtrans

                IF ( k == nzt )  THEN
                   v_int(n) = v_int_l
                ELSE
                   v_int_u = ( ( gg-aa ) * v(k+1,j,i)   + ( gg-bb ) * v(k+1,j,i+1)   &
                             + ( gg-cc ) * v(k+1,j+1,i) + ( gg-dd ) * v(k+1,j+1,i+1) &
                             ) / ( 3.0_wp * gg ) - v_gtrans
                   v_int(n) = v_int_l + ( zv(n) - zu(k) ) / dzw(k+1) *               &
                                     ( v_int_u - v_int_l )
                ENDIF
             ENDIF
          ENDDO
!
!--       Same procedure for interpolation of the w velocity-component
          i = ip + block_offset(nb)%i_off
          j = jp + block_offset(nb)%j_off
          k = kp - 1

          DO  n = start_index(nb), end_index(nb)
             IF ( vertical_particle_advection(particles(n)%group) )  THEN
                x  = xv(n) + ( 0.5_wp - i ) * dx
                y  = yv(n) + ( 0.5_wp - j ) * dy
                aa = x**2          + y**2
                bb = ( dx - x )**2 + y**2
                cc = x**2          + ( dy - y )**2
                dd = ( dx - x )**2 + ( dy - y )**2
                gg = aa + bb + cc + dd

                w_int_l = ( ( gg - aa ) * w(k,j,i)   + ( gg - bb ) * w(k,j,i+1)   &
                          + ( gg - cc ) * w(k,j+1,i) + ( gg - dd ) * w(k,j+1,i+1) &
                          ) / ( 3.0_wp * gg )

                IF ( k == nzt )  THEN
                   w_int(n) = w_int_l
                ELSE
                   w_int_u = ( ( gg-aa ) * w(k+1,j,i)   + &
                               ( gg-bb ) * w(k+1,j,i+1) + &
                               ( gg-cc ) * w(k+1,j+1,i) + &
                               ( gg-dd ) * w(k+1,j+1,i+1) &
                             ) / ( 3.0_wp * gg )
                   w_int(n) = w_int_l + ( zv(n) - zw(k) ) / dzw(k+1) *            &
                              ( w_int_u - w_int_l )
                ENDIF
             ELSE
                w_int(n) = 0.0_wp
             ENDIF
          ENDDO
       ENDDO
    ENDIF

!-- Interpolate and calculate quantities needed for calculating the SGS
!-- velocities
    IF ( use_sgs_for_particles  .AND.  .NOT. cloud_droplets )  THEN

       DO  nb = 0,7

          subbox_at_wall = .FALSE.
!
!--       In case of topography check if subbox is adjacent to a wall
          IF ( .NOT. topography == 'flat' )  THEN
             i = ip + MERGE( -1_iwp , 1_iwp, BTEST( nb, 2 ) )
             j = jp + MERGE( -1_iwp , 1_iwp, BTEST( nb, 1 ) )
             k = kp + MERGE( -1_iwp , 1_iwp, BTEST( nb, 0 ) )
             IF ( .NOT. BTEST(wall_flags_total_0(k,  jp, ip), 0) .OR.             &
                  .NOT. BTEST(wall_flags_total_0(kp, j,  ip), 0) .OR.             &
                  .NOT. BTEST(wall_flags_total_0(kp, jp, i ), 0) )                &
             THEN
                subbox_at_wall = .TRUE.
             ENDIF
          ENDIF
          IF ( subbox_at_wall )  THEN
             e_int(start_index(nb):end_index(nb))     = e(kp,jp,ip) 
             diss_int(start_index(nb):end_index(nb))  = diss(kp,jp,ip)
             de_dx_int(start_index(nb):end_index(nb)) = de_dx(kp,jp,ip)
             de_dy_int(start_index(nb):end_index(nb)) = de_dy(kp,jp,ip)
             de_dz_int(start_index(nb):end_index(nb)) = de_dz(kp,jp,ip)
!
!--          Set flag for stochastic equation.
             term_1_2(start_index(nb):end_index(nb)) = 0.0_wp
          ELSE
             i = ip + block_offset(nb)%i_off
             j = jp + block_offset(nb)%j_off
             k = kp + block_offset(nb)%k_off

             DO  n = start_index(nb), end_index(nb)
!
!--             Interpolate TKE
                x  = xv(n) + ( 0.5_wp - i ) * dx
                y  = yv(n) + ( 0.5_wp - j ) * dy
                aa = x**2          + y**2
                bb = ( dx - x )**2 + y**2
                cc = x**2          + ( dy - y )**2
                dd = ( dx - x )**2 + ( dy - y )**2
                gg = aa + bb + cc + dd

                e_int_l = ( ( gg-aa ) * e(k,j,i)   + ( gg-bb ) * e(k,j,i+1)   &
                          + ( gg-cc ) * e(k,j+1,i) + ( gg-dd ) * e(k,j+1,i+1) &
                          ) / ( 3.0_wp * gg )

                IF ( k+1 == nzt+1 )  THEN
                   e_int(n) = e_int_l
                ELSE
                   e_int_u = ( ( gg - aa ) * e(k+1,j,i)   + &
                               ( gg - bb ) * e(k+1,j,i+1) + &
                               ( gg - cc ) * e(k+1,j+1,i) + &
                               ( gg - dd ) * e(k+1,j+1,i+1) &
                            ) / ( 3.0_wp * gg )
                   e_int(n) = e_int_l + ( zv(n) - zu(k) ) / dzw(k+1) *            &
                                     ( e_int_u - e_int_l )
                ENDIF
!
!--             Needed to avoid NaN particle velocities (this might not be
!--             required any more)
                IF ( e_int(n) <= 0.0_wp )  THEN
                   e_int(n) = 1.0E-20_wp
                ENDIF
!
!--             Interpolate the TKE gradient along x (adopt incides i,j,k and
!--             all position variables from above (TKE))
                de_dx_int_l = ( ( gg - aa ) * de_dx(k,j,i)   + &
                                ( gg - bb ) * de_dx(k,j,i+1) + &
                                ( gg - cc ) * de_dx(k,j+1,i) + &
                                ( gg - dd ) * de_dx(k,j+1,i+1) &
                               ) / ( 3.0_wp * gg )

                IF ( ( k+1 == nzt+1 )  .OR.  ( k == nzb ) )  THEN
                   de_dx_int(n) = de_dx_int_l
                ELSE
                   de_dx_int_u = ( ( gg - aa ) * de_dx(k+1,j,i)   + &
                                   ( gg - bb ) * de_dx(k+1,j,i+1) + &
                                   ( gg - cc ) * de_dx(k+1,j+1,i) + &
                                   ( gg - dd ) * de_dx(k+1,j+1,i+1) &
                                  ) / ( 3.0_wp * gg )
                   de_dx_int(n) = de_dx_int_l + ( zv(n) - zu(k) ) / dzw(k+1) *    &
                                              ( de_dx_int_u - de_dx_int_l )
                ENDIF
!
!--             Interpolate the TKE gradient along y
                de_dy_int_l = ( ( gg - aa ) * de_dy(k,j,i)   + &
                                ( gg - bb ) * de_dy(k,j,i+1) + &
                                ( gg - cc ) * de_dy(k,j+1,i) + &
                                ( gg - dd ) * de_dy(k,j+1,i+1) &
                               ) / ( 3.0_wp * gg )
                IF ( ( k+1 == nzt+1 )  .OR.  ( k == nzb ) )  THEN
                   de_dy_int(n) = de_dy_int_l
                ELSE
                   de_dy_int_u = ( ( gg - aa ) * de_dy(k+1,j,i)   + &
                                   ( gg - bb ) * de_dy(k+1,j,i+1) + &
                                   ( gg - cc ) * de_dy(k+1,j+1,i) + &
                                   ( gg - dd ) * de_dy(k+1,j+1,i+1) &
                                  ) / ( 3.0_wp * gg )
                      de_dy_int(n) = de_dy_int_l + ( zv(n) - zu(k) ) / dzw(k+1) * &
                                                 ( de_dy_int_u - de_dy_int_l )
                ENDIF

!
!--             Interpolate the TKE gradient along z
                IF ( zv(n) < 0.5_wp * dz(1) )  THEN
                   de_dz_int(n) = 0.0_wp
                ELSE
                   de_dz_int_l = ( ( gg - aa ) * de_dz(k,j,i)   + &
                                   ( gg - bb ) * de_dz(k,j,i+1) + &
                                   ( gg - cc ) * de_dz(k,j+1,i) + &
                                   ( gg - dd ) * de_dz(k,j+1,i+1) &
                                  ) / ( 3.0_wp * gg )

                   IF ( ( k+1 == nzt+1 )  .OR.  ( k == nzb ) )  THEN
                      de_dz_int(n) = de_dz_int_l
                   ELSE
                      de_dz_int_u = ( ( gg - aa ) * de_dz(k+1,j,i)   + &
                                      ( gg - bb ) * de_dz(k+1,j,i+1) + &
                                      ( gg - cc ) * de_dz(k+1,j+1,i) + &
                                      ( gg - dd ) * de_dz(k+1,j+1,i+1) &
                                     ) / ( 3.0_wp * gg )
                      de_dz_int(n) = de_dz_int_l + ( zv(n) - zu(k) ) / dzw(k+1) * &
                                                 ( de_dz_int_u - de_dz_int_l )
                   ENDIF
                ENDIF

!
!--             Interpolate the dissipation of TKE
                diss_int_l = ( ( gg - aa ) * diss(k,j,i)   + &
                               ( gg - bb ) * diss(k,j,i+1) + &
                               ( gg - cc ) * diss(k,j+1,i) + &
                               ( gg - dd ) * diss(k,j+1,i+1) &
                               ) / ( 3.0_wp * gg )

                IF ( k == nzt )  THEN
                   diss_int(n) = diss_int_l
                ELSE
                   diss_int_u = ( ( gg - aa ) * diss(k+1,j,i)   + &
                                  ( gg - bb ) * diss(k+1,j,i+1) + &
                                  ( gg - cc ) * diss(k+1,j+1,i) + &
                                  ( gg - dd ) * diss(k+1,j+1,i+1) &
                                 ) / ( 3.0_wp * gg )
                   diss_int(n) = diss_int_l + ( zv(n) - zu(k) ) / dzw(k+1) *      &
                                            ( diss_int_u - diss_int_l )
                ENDIF

!
!--             Set flag for stochastic equation.
                term_1_2(n) = 1.0_wp
             ENDDO
          ENDIF
       ENDDO

       DO  nb = 0,7
          i = ip + block_offset(nb)%i_off
          j = jp + block_offset(nb)%j_off
          k = kp + block_offset(nb)%k_off

          DO  n = start_index(nb), end_index(nb)
!
!--          Vertical interpolation of the horizontally averaged SGS TKE and
!--          resolved-scale velocity variances and use the interpolated values
!--          to calculate the coefficient fs, which is a measure of the ratio
!--          of the subgrid-scale turbulent kinetic energy to the total amount
!--          of turbulent kinetic energy.
             IF ( k == 0 )  THEN
                e_mean_int = hom(0,1,8,0)
             ELSE
                e_mean_int = hom(k,1,8,0) +                                    &
                                           ( hom(k+1,1,8,0) - hom(k,1,8,0) ) / &
                                           ( zu(k+1) - zu(k) ) *               &
                                           ( zv(n) - zu(k) )
             ENDIF

             kw = kp - 1

             IF ( k == 0 )  THEN
                aa  = hom(k+1,1,30,0)  * ( zv(n) / &
                                         ( 0.5_wp * ( zu(k+1) - zu(k) ) ) )
                bb  = hom(k+1,1,31,0)  * ( zv(n) / &
                                         ( 0.5_wp * ( zu(k+1) - zu(k) ) ) )
                cc  = hom(kw+1,1,32,0) * ( zv(n) / &
                                         ( 1.0_wp * ( zw(kw+1) - zw(kw) ) ) )
             ELSE
                aa  = hom(k,1,30,0) + ( hom(k+1,1,30,0) - hom(k,1,30,0) ) *    &
                           ( ( zv(n) - zu(k) ) / ( zu(k+1) - zu(k) ) )
                bb  = hom(k,1,31,0) + ( hom(k+1,1,31,0) - hom(k,1,31,0) ) *    &
                           ( ( zv(n) - zu(k) ) / ( zu(k+1) - zu(k) ) )
                cc  = hom(kw,1,32,0) + ( hom(kw+1,1,32,0)-hom(kw,1,32,0) ) *   &
                           ( ( zv(n) - zw(kw) ) / ( zw(kw+1)-zw(kw) ) )
             ENDIF

             vv_int = ( 1.0_wp / 3.0_wp ) * ( aa + bb + cc )
!
!--          Needed to avoid NaN particle velocities. The value of 1.0 is just
!--          an educated guess for the given case.
             IF ( vv_int + ( 2.0_wp / 3.0_wp ) * e_mean_int == 0.0_wp )  THEN
                fs_int(n) = 1.0_wp
             ELSE
                fs_int(n) = ( 2.0_wp / 3.0_wp ) * e_mean_int /                 &
                            ( vv_int + ( 2.0_wp / 3.0_wp ) * e_mean_int )
             ENDIF

          ENDDO
       ENDDO

       DO  nb = 0, 7
          DO  n = start_index(nb), end_index(nb)
             rg(n,1) = random_gauss( iran_part, 5.0_wp )
             rg(n,2) = random_gauss( iran_part, 5.0_wp )
             rg(n,3) = random_gauss( iran_part, 5.0_wp )
          ENDDO
       ENDDO

       DO  nb = 0, 7
          DO  n = start_index(nb), end_index(nb)

!
!--          Calculate the Lagrangian timescale according to Weil et al. (2004).
             lagr_timescale(n) = ( 4.0_wp * e_int(n) + 1E-20_wp ) / &
                              ( 3.0_wp * fs_int(n) * c_0 * diss_int(n) + 1E-20_wp )

!
!--          Calculate the next particle timestep. dt_gap is the time needed to
!--          complete the current LES timestep.
             dt_gap(n) = dt_3d - particles(n)%dt_sum
             dt_particle(n) = MIN( dt_3d, 0.025_wp * lagr_timescale(n), dt_gap(n) )
             particles(n)%aux1 = lagr_timescale(n)
             particles(n)%aux2 = dt_gap(n)
!
!--          The particle timestep should not be too small in order to prevent
!--          the number of particle timesteps of getting too large
             IF ( dt_particle(n) < dt_min_part )  THEN
                IF ( dt_min_part < dt_gap(n) )  THEN
                   dt_particle(n) = dt_min_part
                ELSE
                   dt_particle(n) = dt_gap(n)
                ENDIF
             ENDIF

             rvar1_temp(n) = particles(n)%rvar1
             rvar2_temp(n) = particles(n)%rvar2
             rvar3_temp(n) = particles(n)%rvar3
!
!--          Calculate the SGS velocity components
             IF ( particles(n)%age == 0.0_wp )  THEN
!
!--             For new particles the SGS components are derived from the SGS
!--             TKE. Limit the Gaussian random number to the interval
!--             [-5.0*sigma, 5.0*sigma] in order to prevent the SGS velocities
!--             from becoming unrealistically large.
                rvar1_temp(n) = SQRT( 2.0_wp * sgs_wf_part * e_int(n)          &
                                          + 1E-20_wp ) * ( rg(n,1) - 1.0_wp )
                rvar2_temp(n) = SQRT( 2.0_wp * sgs_wf_part * e_int(n)          &
                                          + 1E-20_wp ) * ( rg(n,2) - 1.0_wp )
                rvar3_temp(n) = SQRT( 2.0_wp * sgs_wf_part * e_int(n)          &
                                          + 1E-20_wp ) * ( rg(n,3) - 1.0_wp )
             ELSE
!
!--             Restriction of the size of the new timestep: compared to the 
!--             previous timestep the increase must not exceed 200%. First,
!--             check if age > age_m, in order to prevent that particles get zero
!--             timestep.
                dt_particle_m = MERGE( dt_particle(n),                         &
                                       particles(n)%age - particles(n)%age_m,  &
                                       particles(n)%age - particles(n)%age_m < &
                                       1E-8_wp )
                IF ( dt_particle(n) > 2.0_wp * dt_particle_m )  THEN
                   dt_particle(n) = 2.0_wp * dt_particle_m
                ENDIF

!--             For old particles the SGS components are correlated with the
!--             values from the previous timestep. Random numbers have also to
!--             be limited (see above).
!--             As negative values for the subgrid TKE are not allowed, the
!--             change of the subgrid TKE with time cannot be smaller than
!--             -e_int(n)/dt_particle. This value is used as a lower boundary
!--             value for the change of TKE 
                de_dt_min = - e_int(n) / dt_particle(n)

                de_dt = ( e_int(n) - particles(n)%e_m ) / dt_particle_m

                IF ( de_dt < de_dt_min )  THEN
                   de_dt = de_dt_min
                ENDIF

                CALL weil_stochastic_eq( rvar1_temp(n), fs_int(n), e_int(n),    &
                                        de_dx_int(n), de_dt, diss_int(n),       &
                                        dt_particle(n), rg(n,1), term_1_2(n) )

                CALL weil_stochastic_eq( rvar2_temp(n), fs_int(n), e_int(n),    &
                                        de_dy_int(n), de_dt, diss_int(n),       &
                                        dt_particle(n), rg(n,2), term_1_2(n) )

                CALL weil_stochastic_eq( rvar3_temp(n), fs_int(n), e_int(n),    &
                                        de_dz_int(n), de_dt, diss_int(n),       &
                                        dt_particle(n), rg(n,3), term_1_2(n) )

             ENDIF

          ENDDO
       ENDDO
!
!--    Check if the added SGS velocities result in a violation of the CFL-
!--    criterion. If yes, limt the SGS particle speed to match the
!--    CFL criterion. Note, a re-calculation of the SGS particle speed with
!--    smaller timestep does not necessarily fulfill the CFL criterion as the 
!--    new SGS speed can be even larger (due to the random term with scales with
!--    the square-root of dt_particle, for small dt the random contribution increases).
!--    Thus, we would need to re-calculate the SGS speeds as long as they would
!--    fulfill the requirements, which could become computationally expensive, 
!--    Hence, we just limit them.
       dz_temp = zw(kp)-zw(kp-1)

       DO  nb = 0, 7
          DO  n = start_index(nb), end_index(nb)
             IF ( ABS( u_int(n) + rvar1_temp(n) ) > ( dx      / dt_particle(n) )  .OR.   &
                  ABS( v_int(n) + rvar2_temp(n) ) > ( dy      / dt_particle(n) )  .OR.   &
                  ABS( w_int(n) + rvar3_temp(n) ) > ( dz_temp / dt_particle(n) ) )  THEN
!
!--             If total speed exceeds the allowed speed according to CFL 
!--             criterion, limit the SGS speed to
!--             dx_i / dt_particle - u_resolved_i, considering a safty factor.
                rvar1_temp(n) = MERGE( rvar1_temp(n),                          &
                                       0.9_wp *                                &
                                       SIGN( dx / dt_particle(n)               &
                                           - ABS( u_int(n) ), rvar1_temp(n) ), &
                                       ABS( u_int(n) + rvar1_temp(n) ) <       &
                                       ( dx / dt_particle(n) ) )
                rvar2_temp(n) = MERGE( rvar2_temp(n),                          &
                                       0.9_wp *                                &
                                       SIGN( dy / dt_particle(n)               &
                                           - ABS( v_int(n) ), rvar2_temp(n) ), &
                                       ABS( v_int(n) + rvar2_temp(n) ) <       &
                                       ( dy / dt_particle(n) ) )
                rvar3_temp(n) = MERGE( rvar3_temp(n),                          &
                                       0.9_wp *                                &
                                       SIGN( zw(kp)-zw(kp-1) / dt_particle(n)  &
                                           - ABS( w_int(n) ), rvar3_temp(n) ), &
                                       ABS( w_int(n) + rvar3_temp(n) ) <       &
                                       ( zw(kp)-zw(kp-1) / dt_particle(n) ) )
             ENDIF
!
!--          Update particle velocites 
             particles(n)%rvar1 = rvar1_temp(n)
             particles(n)%rvar2 = rvar2_temp(n)
             particles(n)%rvar3 = rvar3_temp(n)
             u_int(n) = u_int(n) + particles(n)%rvar1
             v_int(n) = v_int(n) + particles(n)%rvar2
             w_int(n) = w_int(n) + particles(n)%rvar3
!
!--          Store the SGS TKE of the current timelevel which is needed for
!--          for calculating the SGS particle velocities at the next timestep
             particles(n)%e_m = e_int(n)
          ENDDO
       ENDDO

    ELSE
!
!--    If no SGS velocities are used, only the particle timestep has to
!--    be set
       dt_particle = dt_3d

    ENDIF

    dens_ratio = particle_groups(particles(1:number_of_particles)%group)%density_ratio
    IF ( ANY( dens_ratio == 0.0_wp ) )  THEN
!
!--    Decide whether the particle loop runs over the subboxes or only over 1,
!--    number_of_particles. This depends on the selected interpolation method.
!--    If particle interpolation method is not trilinear, then the sorting within
!--    subboxes is not required. However, therefore the index start_index(nb) and
!--    end_index(nb) are not defined and the loops are still over
!--    number_of_particles. @todo find a more generic way to write this loop or
!--    delete trilinear interpolation
       IF ( interpolation_trilinear )  THEN
          subbox_start = 0
          subbox_end   = 7
       ELSE
          subbox_start = 1
          subbox_end   = 1
       ENDIF
!
!--    loop over subboxes. In case of simple interpolation scheme no subboxes
!--    are introduced, as they are not required. Accordingly, this loops goes
!--    from 1 to 1.
       DO  nb = subbox_start, subbox_end
          IF ( interpolation_trilinear )  THEN
             particle_start = start_index(nb)
             particle_end   = end_index(nb)
          ELSE
             particle_start = 1
             particle_end   = number_of_particles
          ENDIF
!
!--         Loop from particle start to particle end
            DO  n = particle_start, particle_end

!
!--          Particle advection
             IF ( dens_ratio(n) == 0.0_wp )  THEN
!
!--             Pure passive transport (without particle inertia)
                particles(n)%x = xv(n) + u_int(n) * dt_particle(n)
                particles(n)%y = yv(n) + v_int(n) * dt_particle(n)
                particles(n)%z = zv(n) + w_int(n) * dt_particle(n)

                particles(n)%speed_x = u_int(n)
                particles(n)%speed_y = v_int(n)
                particles(n)%speed_z = w_int(n)

             ELSE
!
!--             Transport of particles with inertia
                particles(n)%x = particles(n)%x + particles(n)%speed_x * &
                                                  dt_particle(n)
                particles(n)%y = particles(n)%y + particles(n)%speed_y * &
                                                  dt_particle(n)
                particles(n)%z = particles(n)%z + particles(n)%speed_z * &
                                                  dt_particle(n)

!
!--             Update of the particle velocity
                IF ( cloud_droplets )  THEN
!
!--                Terminal velocity is computed for vertical direction (Rogers et
!--                al., 1993, J. Appl. Meteorol.)
                   diameter = particles(n)%radius * 2000.0_wp !diameter in mm
                   IF ( diameter <= d0_rog )  THEN
                      w_s = k_cap_rog * diameter * ( 1.0_wp - EXP( -k_low_rog * diameter ) )
                   ELSE
                      w_s = a_rog - b_rog * EXP( -c_rog * diameter )
                   ENDIF

!
!--                If selected, add random velocities following Soelch and Kaercher
!--                (2010, Q. J. R. Meteorol. Soc.)
                   IF ( use_sgs_for_particles )  THEN
                      lagr_timescale(n) = km(kp,jp,ip) / MAX( e(kp,jp,ip), 1.0E-20_wp )
                      RL             = EXP( -1.0_wp * dt_3d / MAX( lagr_timescale(n), &
                                             1.0E-20_wp ) )
                      sigma          = SQRT( e(kp,jp,ip) )

                      rg1 = random_gauss( iran_part, 5.0_wp ) - 1.0_wp
                      rg2 = random_gauss( iran_part, 5.0_wp ) - 1.0_wp
                      rg3 = random_gauss( iran_part, 5.0_wp ) - 1.0_wp

                      particles(n)%rvar1 = RL * particles(n)%rvar1 +              &
                                           SQRT( 1.0_wp - RL**2 ) * sigma * rg1
                      particles(n)%rvar2 = RL * particles(n)%rvar2 +              &
                                           SQRT( 1.0_wp - RL**2 ) * sigma * rg2
                      particles(n)%rvar3 = RL * particles(n)%rvar3 +              &
                                           SQRT( 1.0_wp - RL**2 ) * sigma * rg3

                      particles(n)%speed_x = u_int(n) + particles(n)%rvar1
                      particles(n)%speed_y = v_int(n) + particles(n)%rvar2
                      particles(n)%speed_z = w_int(n) + particles(n)%rvar3 - w_s
                   ELSE
                      particles(n)%speed_x = u_int(n)
                      particles(n)%speed_y = v_int(n)
                      particles(n)%speed_z = w_int(n) - w_s
                   ENDIF

                ELSE

                   IF ( use_sgs_for_particles )  THEN
                      exp_arg  = particle_groups(particles(n)%group)%exp_arg
                      exp_term = EXP( -exp_arg * dt_particle(n) )
                   ELSE
                      exp_arg  = particle_groups(particles(n)%group)%exp_arg
                      exp_term = particle_groups(particles(n)%group)%exp_term
                   ENDIF
                   particles(n)%speed_x = particles(n)%speed_x * exp_term +         &
                                          u_int(n) * ( 1.0_wp - exp_term )
                   particles(n)%speed_y = particles(n)%speed_y * exp_term +         &
                                          v_int(n) * ( 1.0_wp - exp_term )
                   particles(n)%speed_z = particles(n)%speed_z * exp_term +         &
                                          ( w_int(n) - ( 1.0_wp - dens_ratio(n) ) * &
                                          g / exp_arg ) * ( 1.0_wp - exp_term )
                ENDIF

             ENDIF
          ENDDO
       ENDDO

    ELSE
!
!--    Decide whether the particle loop runs over the subboxes or only over 1,
!--    number_of_particles. This depends on the selected interpolation method.
       IF ( interpolation_trilinear )  THEN
          subbox_start = 0
          subbox_end   = 7
       ELSE
          subbox_start = 1
          subbox_end   = 1
       ENDIF
!--    loop over subboxes. In case of simple interpolation scheme no subboxes
!--    are introduced, as they are not required. Accordingly, this loops goes
!--    from 1 to 1.
       DO  nb = subbox_start, subbox_end
          IF ( interpolation_trilinear )  THEN
             particle_start = start_index(nb)
             particle_end   = end_index(nb)
          ELSE
             particle_start = 1
             particle_end   = number_of_particles
          ENDIF
!
!--         Loop from particle start to particle end
            DO  n = particle_start, particle_end

!
!--          Transport of particles with inertia
             particles(n)%x = xv(n) + particles(n)%speed_x * dt_particle(n)
             particles(n)%y = yv(n) + particles(n)%speed_y * dt_particle(n)
             particles(n)%z = zv(n) + particles(n)%speed_z * dt_particle(n)
!
!--          Update of the particle velocity
             IF ( cloud_droplets )  THEN
!
!--             Terminal velocity is computed for vertical direction (Rogers et al.,
!--             1993, J. Appl. Meteorol.)
                diameter = particles(n)%radius * 2000.0_wp !diameter in mm
                IF ( diameter <= d0_rog )  THEN
                   w_s = k_cap_rog * diameter * ( 1.0_wp - EXP( -k_low_rog * diameter ) )
                ELSE
                   w_s = a_rog - b_rog * EXP( -c_rog * diameter )
                ENDIF

!
!--             If selected, add random velocities following Soelch and Kaercher
!--             (2010, Q. J. R. Meteorol. Soc.)
                IF ( use_sgs_for_particles )  THEN
                    lagr_timescale(n) = km(kp,jp,ip) / MAX( e(kp,jp,ip), 1.0E-20_wp )
                     RL             = EXP( -1.0_wp * dt_3d / MAX( lagr_timescale(n), &
                                             1.0E-20_wp ) )
                    sigma          = SQRT( e(kp,jp,ip) )

                    rg1 = random_gauss( iran_part, 5.0_wp ) - 1.0_wp
                    rg2 = random_gauss( iran_part, 5.0_wp ) - 1.0_wp
                    rg3 = random_gauss( iran_part, 5.0_wp ) - 1.0_wp

                    particles(n)%rvar1 = RL * particles(n)%rvar1 +                &
                                         SQRT( 1.0_wp - RL**2 ) * sigma * rg1
                    particles(n)%rvar2 = RL * particles(n)%rvar2 +                &
                                         SQRT( 1.0_wp - RL**2 ) * sigma * rg2
                    particles(n)%rvar3 = RL * particles(n)%rvar3 +                &
                                         SQRT( 1.0_wp - RL**2 ) * sigma * rg3

                    particles(n)%speed_x = u_int(n) + particles(n)%rvar1
                    particles(n)%speed_y = v_int(n) + particles(n)%rvar2
                    particles(n)%speed_z = w_int(n) + particles(n)%rvar3 - w_s
                ELSE
                    particles(n)%speed_x = u_int(n)
                    particles(n)%speed_y = v_int(n)
                    particles(n)%speed_z = w_int(n) - w_s
                ENDIF

             ELSE

                IF ( use_sgs_for_particles )  THEN
                   exp_arg  = particle_groups(particles(n)%group)%exp_arg
                   exp_term = EXP( -exp_arg * dt_particle(n) )
                ELSE
                   exp_arg  = particle_groups(particles(n)%group)%exp_arg
                   exp_term = particle_groups(particles(n)%group)%exp_term
                ENDIF
                particles(n)%speed_x = particles(n)%speed_x * exp_term +             &
                                       u_int(n) * ( 1.0_wp - exp_term )
                particles(n)%speed_y = particles(n)%speed_y * exp_term +             &
                                       v_int(n) * ( 1.0_wp - exp_term )
                particles(n)%speed_z = particles(n)%speed_z * exp_term +             &
                                       ( w_int(n) - ( 1.0_wp - dens_ratio(n) ) * g / &
                                       exp_arg ) * ( 1.0_wp - exp_term )
             ENDIF
          ENDDO
       ENDDO

    ENDIF

!
!-- Store the old age of the particle ( needed to prevent that a
!-- particle crosses several PEs during one timestep, and for the
!-- evaluation of the subgrid particle velocity fluctuations )
    particles(1:number_of_particles)%age_m = particles(1:number_of_particles)%age

!
!--    loop over subboxes. In case of simple interpolation scheme no subboxes
!--    are introduced, as they are not required. Accordingly, this loops goes
!--    from 1 to 1.
!
!-- Decide whether the particle loop runs over the subboxes or only over 1,
!-- number_of_particles. This depends on the selected interpolation method.
    IF ( interpolation_trilinear )  THEN
       subbox_start = 0
       subbox_end   = 7
    ELSE
       subbox_start = 1
       subbox_end   = 1
    ENDIF
    DO  nb = subbox_start, subbox_end
       IF ( interpolation_trilinear )  THEN
          particle_start = start_index(nb)
          particle_end   = end_index(nb)
       ELSE
          particle_start = 1
          particle_end   = number_of_particles
       ENDIF
!
!--    Loop from particle start to particle end and increment the particle 
!--    age and the total time that the particle has advanced within the 
!--    particle timestep procedure.
       DO  n = particle_start, particle_end
          particles(n)%age    = particles(n)%age    + dt_particle(n)
          particles(n)%dt_sum = particles(n)%dt_sum + dt_particle(n)
       ENDDO
!
!--    Particles that leave the child domain during the SGS-timestep loop
!--    must not continue timestepping until they are transferred to the 
!--    parent. Hence, set their dt_sum to dt.
       IF ( child_domain  .AND.  use_sgs_for_particles )  THEN
          DO  n = particle_start, particle_end
             IF ( particles(n)%x < 0.0_wp         .OR.                         &
                  particles(n)%y < 0.0_wp         .OR.                         &
                  particles(n)%x > ( nx+1 ) * dx  .OR.                         &
                  particles(n)%y < ( ny+1 ) * dy )  THEN
                particles(n)%dt_sum = dt_3d
             ENDIF
          ENDDO
       ENDIF
!
!--    Check whether there is still a particle that has not yet completed
!--    the total LES timestep
       DO  n = particle_start, particle_end
          IF ( ( dt_3d - particles(n)%dt_sum ) > 1E-8_wp )                     &
             dt_3d_reached_l = .FALSE.
       ENDDO
    ENDDO

    CALL cpu_log( log_point_s(44), 'lpm_advec', 'pause' )


 END SUBROUTINE lpm_advec

 
!------------------------------------------------------------------------------!  
! Description:
! ------------
!> Calculation of subgrid-scale particle speed using the stochastic model 
!> of Weil et al. (2004, JAS, 61, 2877-2887).
!------------------------------------------------------------------------------! 
 SUBROUTINE weil_stochastic_eq( v_sgs, fs_n, e_n, dedxi_n, dedt_n, diss_n,     &
                                dt_n, rg_n, fac )

    REAL(wp) ::  a1      !< dummy argument
    REAL(wp) ::  dedt_n  !< time derivative of TKE at particle position 
    REAL(wp) ::  dedxi_n !< horizontal derivative of TKE at particle position
    REAL(wp) ::  diss_n  !< dissipation at particle position 
    REAL(wp) ::  dt_n    !< particle timestep
    REAL(wp) ::  e_n     !< TKE at particle position
    REAL(wp) ::  fac     !< flag to identify adjacent topography
    REAL(wp) ::  fs_n    !< weighting factor to prevent that subgrid-scale particle speed becomes too large
    REAL(wp) ::  rg_n    !< random number
    REAL(wp) ::  term1   !< memory term
    REAL(wp) ::  term2   !< drift correction term
    REAL(wp) ::  term3   !< random term
    REAL(wp) ::  v_sgs   !< subgrid-scale velocity component 

!-- At first, limit TKE to a small non-zero number, in order to prevent 
!-- the occurrence of extremely large SGS-velocities in case TKE is zero, 
!-- (could occur at the simulation begin).
    e_n = MAX( e_n, 1E-20_wp )
!
!-- Please note, terms 1 and 2 (drift and memory term, respectively) are 
!-- multiplied by a flag to switch of both terms near topography. 
!-- This is necessary, as both terms may cause a subgrid-scale velocity build up 
!-- if particles are trapped in regions with very small TKE, e.g. in narrow street 
!-- canyons resolved by only a few grid points. Hence, term 1 and term 2 are 
!-- disabled if one of the adjacent grid points belongs to topography. 
!-- Moreover, in this case, the  previous subgrid-scale component is also set
!-- to zero.

    a1 = fs_n * c_0 * diss_n
!
!-- Memory term
    term1 = - a1 * v_sgs * dt_n / ( 4.0_wp * sgs_wf_part * e_n + 1E-20_wp )    &
                 * fac
!
!-- Drift correction term
    term2 = ( ( dedt_n * v_sgs / e_n ) + dedxi_n ) * 0.5_wp * dt_n              &
                 * fac
!
!-- Random term
    term3 = SQRT( MAX( a1, 1E-20_wp ) ) * ( rg_n - 1.0_wp ) * SQRT( dt_n )
!
!-- In cese one of the adjacent grid-boxes belongs to topograhy, the previous
!-- subgrid-scale velocity component is set to zero, in order to prevent a 
!-- velocity build-up. 
!-- This case, set also previous subgrid-scale component to zero. 
    v_sgs = v_sgs * fac + term1 + term2 + term3

 END SUBROUTINE weil_stochastic_eq


!------------------------------------------------------------------------------!
! Description:
! ------------
!> swap timelevel in case of particle advection interpolation 'simple-corrector'
!> This routine is called at the end of one timestep, the velocities are then
!> used for the next timestep
!------------------------------------------------------------------------------!
 SUBROUTINE lpm_swap_timelevel_for_particle_advection

!
!-- save the divergence free velocites of t+1 to use them at the end of the
!-- next time step
    u_t = u
    v_t = v
    w_t = w

 END SUBROUTINE lpm_swap_timelevel_for_particle_advection


!------------------------------------------------------------------------------!  
! Description:
! ------------
!> Boundary conditions for the Lagrangian particles.
!> The routine consists of two different parts. One handles the bottom (flat)
!> and top boundary. In this part, also particles which exceeded their lifetime
!> are deleted.
!> The other part handles the reflection of particles from vertical walls.
!> This part was developed by Jin Zhang during 2006-2007.
!>
!> To do: Code structure for finding the t_index values and for checking the
!> -----  reflection conditions is basically the same for all four cases, so it
!>        should be possible to further simplify/shorten it.
!>
!> THE WALLS PART OF THIS ROUTINE HAS NOT BEEN TESTED FOR OCEAN RUNS SO FAR!!!!
!> (see offset_ocean_*)
!------------------------------------------------------------------------------! 
 SUBROUTINE lpm_boundary_conds( location_bc , i, j, k )

    CHARACTER (LEN=*), INTENT(IN) ::  location_bc !< general mode: boundary conditions at bottom/top of the model domain
                                   !< or at vertical surfaces (buildings, terrain steps)    
    INTEGER(iwp), INTENT(IN) ::  i !< grid index of particle box along x
    INTEGER(iwp), INTENT(IN) ::  j !< grid index of particle box along y
    INTEGER(iwp), INTENT(IN) ::  k !< grid index of particle box along z

    INTEGER(iwp) ::  inc            !< dummy for sorting algorithmus
    INTEGER(iwp) ::  ir             !< dummy for sorting algorithmus
    INTEGER(iwp) ::  i1             !< grid index (x) of old particle position
    INTEGER(iwp) ::  i2             !< grid index (x) of current particle position
    INTEGER(iwp) ::  i3             !< grid index (x) of intermediate particle position
    INTEGER(iwp) ::  index_reset    !< index reset height
    INTEGER(iwp) ::  jr             !< dummy for sorting algorithmus
    INTEGER(iwp) ::  j1             !< grid index (y) of old particle position
    INTEGER(iwp) ::  j2             !< grid index (y) of current particle position
    INTEGER(iwp) ::  j3             !< grid index (y) of intermediate particle position
    INTEGER(iwp) ::  k1             !< grid index (z) of old particle position
    INTEGER(iwp) ::  k2             !< grid index (z) of current particle position
    INTEGER(iwp) ::  k3             !< grid index (z) of intermediate particle position
    INTEGER(iwp) ::  n              !< particle number
    INTEGER(iwp) ::  particles_top  !< maximum reset height
    INTEGER(iwp) ::  t_index        !< running index for intermediate particle timesteps in reflection algorithmus
    INTEGER(iwp) ::  t_index_number !< number of intermediate particle timesteps in reflection algorithmus
    INTEGER(iwp) ::  tmp_x          !< dummy for sorting algorithm
    INTEGER(iwp) ::  tmp_y          !< dummy for sorting algorithm
    INTEGER(iwp) ::  tmp_z          !< dummy for sorting algorithm

    INTEGER(iwp), DIMENSION(0:10) ::  x_ind(0:10) = 0 !< index array (x) of intermediate particle positions
    INTEGER(iwp), DIMENSION(0:10) ::  y_ind(0:10) = 0 !< index array (y) of intermediate particle positions
    INTEGER(iwp), DIMENSION(0:10) ::  z_ind(0:10) = 0 !< index array (z) of intermediate particle positions

    LOGICAL  ::  cross_wall_x    !< flag to check if particle reflection along x is necessary
    LOGICAL  ::  cross_wall_y    !< flag to check if particle reflection along y is necessary
    LOGICAL  ::  cross_wall_z    !< flag to check if particle reflection along z is necessary
    LOGICAL  ::  reflect_x       !< flag to check if particle is already reflected along x
    LOGICAL  ::  reflect_y       !< flag to check if particle is already reflected along y
    LOGICAL  ::  reflect_z       !< flag to check if particle is already reflected along z
    LOGICAL  ::  tmp_reach_x     !< dummy for sorting algorithmus
    LOGICAL  ::  tmp_reach_y     !< dummy for sorting algorithmus
    LOGICAL  ::  tmp_reach_z     !< dummy for sorting algorithmus
    LOGICAL  ::  x_wall_reached  !< flag to check if particle has already reached wall
    LOGICAL  ::  y_wall_reached  !< flag to check if particle has already reached wall
    LOGICAL  ::  z_wall_reached  !< flag to check if particle has already reached wall

    LOGICAL, DIMENSION(0:10) ::  reach_x  !< flag to check if particle is at a yz-wall
    LOGICAL, DIMENSION(0:10) ::  reach_y  !< flag to check if particle is at a xz-wall
    LOGICAL, DIMENSION(0:10) ::  reach_z  !< flag to check if particle is at a xy-wall

    REAL(wp) ::  dt_particle    !< particle timestep
    REAL(wp) ::  eps = 1E-10_wp !< security number to check if particle has reached a wall
    REAL(wp) ::  pos_x          !< intermediate particle position (x)
    REAL(wp) ::  pos_x_old      !< particle position (x) at previous particle timestep
    REAL(wp) ::  pos_y          !< intermediate particle position (y)
    REAL(wp) ::  pos_y_old      !< particle position (y) at previous particle timestep
    REAL(wp) ::  pos_z          !< intermediate particle position (z)
    REAL(wp) ::  pos_z_old      !< particle position (z) at previous particle timestep
    REAL(wp) ::  prt_x          !< current particle position (x)
    REAL(wp) ::  prt_y          !< current particle position (y)
    REAL(wp) ::  prt_z          !< current particle position (z)
    REAL(wp) ::  ran_val        !< location of wall in z
    REAL(wp) ::  reset_top      !< location of wall in z
    REAL(wp) ::  t_old          !< previous reflection time
    REAL(wp) ::  tmp_t          !< dummy for sorting algorithmus
    REAL(wp) ::  xwall          !< location of wall in x
    REAL(wp) ::  ywall          !< location of wall in y
    REAL(wp) ::  zwall          !< location of wall in z

    REAL(wp), DIMENSION(0:10) ::  t  !< reflection time

    SELECT CASE ( location_bc )

       CASE ( 'bottom/top' )

!
!--    Apply boundary conditions to those particles that have crossed the top or
!--    bottom boundary and delete those particles, which are older than allowed
       DO  n = 1, number_of_particles

!
!--       Stop if particles have moved further than the length of one 
!--       PE subdomain (newly released particles have age = age_m!)
          IF ( particles(n)%age /= particles(n)%age_m )  THEN
             IF ( ABS(particles(n)%speed_x) >                                  &
                  ((nxr-nxl+2)*dx)/(particles(n)%age-particles(n)%age_m)  .OR. &
                  ABS(particles(n)%speed_y) >                                  &
                  ((nyn-nys+2)*dy)/(particles(n)%age-particles(n)%age_m) )  THEN

                  WRITE( message_string, * )  'particle too fast.  n = ',  n 
                  CALL message( 'lpm_boundary_conds', 'PA0148', 2, 2, -1, 6, 1 )
             ENDIF
          ENDIF

          IF ( particles(n)%age > particle_maximum_age  .AND.  &
               particles(n)%particle_mask )                              &
          THEN
             particles(n)%particle_mask  = .FALSE.
             deleted_particles = deleted_particles + 1
          ENDIF

          IF ( particles(n)%z >= zw(nz)  .AND.  particles(n)%particle_mask )  THEN
             IF ( ibc_par_t == 1 )  THEN
!
!--             Particle absorption
                particles(n)%particle_mask  = .FALSE.
                deleted_particles = deleted_particles + 1
             ELSEIF ( ibc_par_t == 2 )  THEN
!
!--             Particle reflection
                particles(n)%z       = 2.0_wp * zw(nz) - particles(n)%z
                particles(n)%speed_z = -particles(n)%speed_z
                IF ( use_sgs_for_particles  .AND. &
                     particles(n)%rvar3 > 0.0_wp )  THEN
                   particles(n)%rvar3 = -particles(n)%rvar3
                ENDIF
             ENDIF
          ENDIF

          IF ( particles(n)%z < zw(0)  .AND.  particles(n)%particle_mask )  THEN
             IF ( ibc_par_b == 1 )  THEN
!
!--             Particle absorption
                particles(n)%particle_mask  = .FALSE.
                deleted_particles = deleted_particles + 1
             ELSEIF ( ibc_par_b == 2 )  THEN
!
!--             Particle reflection
                particles(n)%z       = 2.0_wp * zw(0) - particles(n)%z
                particles(n)%speed_z = -particles(n)%speed_z
                IF ( use_sgs_for_particles  .AND. &
                     particles(n)%rvar3 < 0.0_wp )  THEN
                   particles(n)%rvar3 = -particles(n)%rvar3
                ENDIF
             ELSEIF ( ibc_par_b == 3 )  THEN
!
!--             Find reset height. @note this works only in non-strechted cases
                particles_top = INT( pst(1) / dz(1) )
                index_reset = MINLOC( prt_count(nzb+1:particles_top,j,i), DIM = 1 )
                reset_top = zu(index_reset)
                iran_part = iran_part + myid
                ran_val = random_function( iran_part )
                particles(n)%z       = reset_top *  ( 1.0  + ( ran_val / 10.0_wp) )
                particles(n)%speed_z = 0.0_wp
                IF ( curvature_solution_effects )  THEN
                   particles(n)%radius = particles(n)%aux1
                ELSE
                   particles(n)%radius = 1.0E-8
                ENDIF
             ENDIF
          ENDIF
       ENDDO

      CASE ( 'walls' )

       CALL cpu_log( log_point_s(48), 'lpm_wall_reflect', 'start' )

       DO  n = 1, number_of_particles
!
!--       Recalculate particle timestep
          dt_particle = particles(n)%age - particles(n)%age_m
!
!--       Obtain x/y indices for current particle position
          i2 = particles(n)%x * ddx
          j2 = particles(n)%y * ddy
          IF ( zw(k)   < particles(n)%z ) k2 = k + 1
          IF ( zw(k)   > particles(n)%z  .AND.  zw(k-1) < particles(n)%z ) k2 = k
          IF ( zw(k-1) > particles(n)%z ) k2 = k - 1
!
!--       Save current particle positions
          prt_x = particles(n)%x
          prt_y = particles(n)%y
          prt_z = particles(n)%z
!
!--       Recalculate old particle positions
          pos_x_old = particles(n)%x - particles(n)%speed_x * dt_particle
          pos_y_old = particles(n)%y - particles(n)%speed_y * dt_particle
          pos_z_old = particles(n)%z - particles(n)%speed_z * dt_particle
!
!--       Obtain x/y indices for old particle positions
          i1 = i
          j1 = j
          k1 = k
!
!--       Determine horizontal as well as vertical walls at which particle can 
!--       be potentially reflected. 
!--       Start with walls aligned in yz layer.
!--       Wall to the right 
          IF ( prt_x > pos_x_old )  THEN
             xwall = ( i1 + 1 ) * dx
!
!--       Wall to the left
          ELSE
             xwall = i1 * dx
          ENDIF
!
!--       Walls aligned in xz layer
!--       Wall to the north
          IF ( prt_y > pos_y_old )  THEN
             ywall = ( j1 + 1 ) * dy
!--       Wall to the south
          ELSE
             ywall = j1 * dy
          ENDIF

          IF ( prt_z > pos_z_old )  THEN
             zwall = zw(k)
          ELSE
             zwall = zw(k-1)
          ENDIF
!
!--       Initialize flags to check if particle reflection is necessary
          cross_wall_x = .FALSE.
          cross_wall_y = .FALSE.
          cross_wall_z = .FALSE.
!
!--       Initialize flags to check if a wall is reached
          reach_x      = .FALSE.
          reach_y      = .FALSE.
          reach_z      = .FALSE.
!
!--       Initialize flags to check if a particle was already reflected
          reflect_x    = .FALSE.
          reflect_y    = .FALSE.
          reflect_z    = .FALSE.
!
!--       Initialize flags to check if a wall is already crossed.
!--       ( Required to obtain correct indices. )
          x_wall_reached = .FALSE.
          y_wall_reached = .FALSE.
          z_wall_reached = .FALSE.
!
!--       Initialize time array 
          t     = 0.0_wp
!
!--       Check if particle can reach any wall. This case, calculate the 
!--       fractional time needed to reach this wall. Store this fractional
!--       timestep in array t. Moreover, store indices for these grid
!--       boxes where the respective wall belongs to.  
!--       Start with x-direction.
          t_index    = 1
          t(t_index) = ( xwall - pos_x_old )                                   &
                     / MERGE( MAX( prt_x - pos_x_old,  1E-30_wp ),             &
                              MIN( prt_x - pos_x_old, -1E-30_wp ),             &
                              prt_x > pos_x_old )
          x_ind(t_index)   = i2
          y_ind(t_index)   = j1
          z_ind(t_index)   = k1
          reach_x(t_index) = .TRUE.
          reach_y(t_index) = .FALSE.
          reach_z(t_index) = .FALSE.
!
!--       Store these values only if particle really reaches any wall. t must 
!--       be in a interval between [0:1].
          IF ( t(t_index) <= 1.0_wp  .AND.  t(t_index) >= 0.0_wp )  THEN
             t_index      = t_index + 1
             cross_wall_x = .TRUE.
          ENDIF
!
!--       y-direction
          t(t_index) = ( ywall - pos_y_old )                                   &
                     / MERGE( MAX( prt_y - pos_y_old,  1E-30_wp ),             &
                              MIN( prt_y - pos_y_old, -1E-30_wp ),             &
                              prt_y > pos_y_old )
          x_ind(t_index)   = i1
          y_ind(t_index)   = j2
          z_ind(t_index)   = k1
          reach_x(t_index) = .FALSE.
          reach_y(t_index) = .TRUE.
          reach_z(t_index) = .FALSE.
          IF ( t(t_index) <= 1.0_wp  .AND.  t(t_index) >= 0.0_wp )  THEN
             t_index      = t_index + 1
             cross_wall_y = .TRUE.
          ENDIF
!
!--       z-direction
          t(t_index) = (zwall - pos_z_old )                                    &
                     / MERGE( MAX( prt_z - pos_z_old,  1E-30_wp ),             &
                              MIN( prt_z - pos_z_old, -1E-30_wp ),             &
                              prt_z > pos_z_old )

          x_ind(t_index)   = i1
          y_ind(t_index)   = j1
          z_ind(t_index)   = k2
          reach_x(t_index) = .FALSE.
          reach_y(t_index) = .FALSE.
          reach_z(t_index) = .TRUE.
          IF( t(t_index) <= 1.0_wp  .AND.  t(t_index) >= 0.0_wp)  THEN
             t_index      = t_index + 1
             cross_wall_z = .TRUE.
          ENDIF

          t_index_number = t_index - 1
!
!--       Carry out reflection only if particle reaches any wall
          IF ( cross_wall_x  .OR.  cross_wall_y  .OR.  cross_wall_z )  THEN
!
!--          Sort fractional timesteps in ascending order. Also sort the 
!--          corresponding indices and flag according to the time interval a  
!--          particle reaches the respective wall.
             inc = 1
             jr  = 1
             DO WHILE ( inc <= t_index_number )
                inc = 3 * inc + 1
             ENDDO

             DO WHILE ( inc > 1 )
                inc = inc / 3
                DO  ir = inc+1, t_index_number
                   tmp_t       = t(ir)
                   tmp_x       = x_ind(ir)
                   tmp_y       = y_ind(ir)
                   tmp_z       = z_ind(ir)
                   tmp_reach_x = reach_x(ir)
                   tmp_reach_y = reach_y(ir)
                   tmp_reach_z = reach_z(ir)
                   jr    = ir
                   DO WHILE ( t(jr-inc) > tmp_t )
                      t(jr)       = t(jr-inc)
                      x_ind(jr)   = x_ind(jr-inc)
                      y_ind(jr)   = y_ind(jr-inc)
                      z_ind(jr)   = z_ind(jr-inc)
                      reach_x(jr) = reach_x(jr-inc)
                      reach_y(jr) = reach_y(jr-inc)
                      reach_z(jr) = reach_z(jr-inc)
                      jr    = jr - inc
                      IF ( jr <= inc )  EXIT
                   ENDDO
                   t(jr)       = tmp_t
                   x_ind(jr)   = tmp_x
                   y_ind(jr)   = tmp_y
                   z_ind(jr)   = tmp_z
                   reach_x(jr) = tmp_reach_x
                   reach_y(jr) = tmp_reach_y
                   reach_z(jr) = tmp_reach_z
                ENDDO
             ENDDO
!
!--          Initialize temporary particle positions
             pos_x = pos_x_old
             pos_y = pos_y_old
             pos_z = pos_z_old
!
!--          Loop over all times a particle possibly moves into a new grid box
             t_old = 0.0_wp
             DO t_index = 1, t_index_number 
!
!--             Calculate intermediate particle position according to the 
!--             timesteps a particle reaches any wall.
                pos_x = pos_x + ( t(t_index) - t_old ) * dt_particle           &
                                                       * particles(n)%speed_x
                pos_y = pos_y + ( t(t_index) - t_old ) * dt_particle           &
                                                       * particles(n)%speed_y
                pos_z = pos_z + ( t(t_index) - t_old ) * dt_particle           &
                                                       * particles(n)%speed_z
!
!--             Obtain x/y grid indices for intermediate particle position from
!--             sorted index array
                i3 = x_ind(t_index)
                j3 = y_ind(t_index)
                k3 = z_ind(t_index)
!
!--             Check which wall is already reached
                IF ( .NOT. x_wall_reached )  x_wall_reached = reach_x(t_index) 
                IF ( .NOT. y_wall_reached )  y_wall_reached = reach_y(t_index)
                IF ( .NOT. z_wall_reached )  z_wall_reached = reach_z(t_index)
!
!--             Check if a particle needs to be reflected at any yz-wall. If 
!--             necessary, carry out reflection. Please note, a security 
!--             constant is required, as the particle position does not 
!--             necessarily exactly match the wall location due to rounding 
!--             errors.
                IF ( reach_x(t_index)                      .AND.               & 
                     ABS( pos_x - xwall ) < eps            .AND.               &
                     .NOT. BTEST(wall_flags_total_0(k3,j3,i3),0) .AND.         &
                     .NOT. reflect_x )  THEN
! 
! 
!--                Reflection in x-direction. 
!--                Ensure correct reflection by MIN/MAX functions, depending on
!--                direction of particle transport. 
!--                Due to rounding errors pos_x does not exactly match the wall
!--                location, leading to erroneous reflection.             
                   pos_x = MERGE( MIN( 2.0_wp * xwall - pos_x, xwall ),        &
                                  MAX( 2.0_wp * xwall - pos_x, xwall ),        &
                                  particles(n)%x > xwall )
!
!--                Change sign of particle speed                     
                   particles(n)%speed_x = - particles(n)%speed_x
!
!--                Also change sign of subgrid-scale particle speed
                   particles(n)%rvar1 = - particles(n)%rvar1
!
!--                Set flag that reflection along x is already done
                   reflect_x          = .TRUE.
!
!--                As the particle does not cross any further yz-wall during 
!--                this timestep, set further x-indices to the current one.
                   x_ind(t_index:t_index_number) = i1
!
!--             If particle already reached the wall but was not reflected, 
!--             set further x-indices to the new one.
                ELSEIF ( x_wall_reached .AND. .NOT. reflect_x )  THEN
                    x_ind(t_index:t_index_number) = i2
                ENDIF !particle reflection in x direction done

!
!--             Check if a particle needs to be reflected at any xz-wall. If 
!--             necessary, carry out reflection. Please note, a security 
!--             constant is required, as the particle position does not 
!--             necessarily exactly match the wall location due to rounding 
!--             errors. 
                IF ( reach_y(t_index)                      .AND.               & 
                     ABS( pos_y - ywall ) < eps            .AND.               &
                     .NOT. BTEST(wall_flags_total_0(k3,j3,i3),0) .AND.         &
                     .NOT. reflect_y )  THEN
! 
! 
!--                Reflection in y-direction. 
!--                Ensure correct reflection by MIN/MAX functions, depending on
!--                direction of particle transport. 
!--                Due to rounding errors pos_y does not exactly match the wall
!--                location, leading to erroneous reflection.             
                   pos_y = MERGE( MIN( 2.0_wp * ywall - pos_y, ywall ),        &
                                  MAX( 2.0_wp * ywall - pos_y, ywall ),        &
                                  particles(n)%y > ywall )
!
!--                Change sign of particle speed                     
                   particles(n)%speed_y = - particles(n)%speed_y
!
!--                Also change sign of subgrid-scale particle speed
                   particles(n)%rvar2 = - particles(n)%rvar2
!
!--                Set flag that reflection along y is already done
                   reflect_y          = .TRUE.
!
!--                As the particle does not cross any further xz-wall during 
!--                this timestep, set further y-indices to the current one.
                   y_ind(t_index:t_index_number) = j1
!
!--             If particle already reached the wall but was not reflected, 
!--             set further y-indices to the new one.
                ELSEIF ( y_wall_reached .AND. .NOT. reflect_y )  THEN
                    y_ind(t_index:t_index_number) = j2
                ENDIF !particle reflection in y direction done

!
!--             Check if a particle needs to be reflected at any xy-wall. If 
!--             necessary, carry out reflection. Please note, a security 
!--             constant is required, as the particle position does not 
!--             necessarily exactly match the wall location due to rounding 
!--             errors. 
                IF ( reach_z(t_index)                      .AND.               & 
                     ABS( pos_z - zwall ) < eps            .AND.               &
                     .NOT. BTEST(wall_flags_total_0(k3,j3,i3),0) .AND.         &
                     .NOT. reflect_z )  THEN
! 
! 
!--                Reflection in z-direction. 
!--                Ensure correct reflection by MIN/MAX functions, depending on
!--                direction of particle transport. 
!--                Due to rounding errors pos_z does not exactly match the wall
!--                location, leading to erroneous reflection.             
                   pos_z = MERGE( MIN( 2.0_wp * zwall - pos_z, zwall ),        &
                                  MAX( 2.0_wp * zwall - pos_z, zwall ),        &
                                  particles(n)%z > zwall )
!
!--                Change sign of particle speed                     
                   particles(n)%speed_z = - particles(n)%speed_z
!
!--                Also change sign of subgrid-scale particle speed
                   particles(n)%rvar3 = - particles(n)%rvar3
!
!--                Set flag that reflection along z is already done
                   reflect_z          = .TRUE.
!
!--                As the particle does not cross any further xy-wall during 
!--                this timestep, set further z-indices to the current one.
                   z_ind(t_index:t_index_number) = k1
!
!--             If particle already reached the wall but was not reflected, 
!--             set further z-indices to the new one.
                ELSEIF ( z_wall_reached .AND. .NOT. reflect_z )  THEN
                    z_ind(t_index:t_index_number) = k2
                ENDIF !particle reflection in z direction done                

!
!--             Swap time
                t_old = t(t_index)

             ENDDO
!
!--          If a particle was reflected, calculate final position from last
!--          intermediate position.
             IF ( reflect_x  .OR.  reflect_y  .OR.  reflect_z )  THEN

                particles(n)%x = pos_x + ( 1.0_wp - t_old ) * dt_particle      &
                                                         * particles(n)%speed_x
                particles(n)%y = pos_y + ( 1.0_wp - t_old ) * dt_particle      &
                                                         * particles(n)%speed_y
                particles(n)%z = pos_z + ( 1.0_wp - t_old ) * dt_particle      &
                                                         * particles(n)%speed_z

             ENDIF

          ENDIF

       ENDDO

       CALL cpu_log( log_point_s(48), 'lpm_wall_reflect', 'stop' )

       CASE DEFAULT
          CONTINUE

    END SELECT

 END SUBROUTINE lpm_boundary_conds 


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculates change in droplet radius by condensation/evaporation, using
!> either an analytic formula or by numerically integrating the radius growth
!> equation including curvature and solution effects using Rosenbrocks method
!> (see Numerical recipes in FORTRAN, 2nd edition, p. 731).
!> The analytical formula and growth equation follow those given in
!> Rogers and Yau (A short course in cloud physics, 3rd edition, p. 102/103).
!------------------------------------------------------------------------------!
 SUBROUTINE lpm_droplet_condensation (i,j,k)

    INTEGER(iwp), INTENT(IN) ::  i              !<
    INTEGER(iwp), INTENT(IN) ::  j              !<
    INTEGER(iwp), INTENT(IN) ::  k              !<
    INTEGER(iwp) ::  n                          !<

    REAL(wp) ::  afactor                       !< curvature effects
    REAL(wp) ::  arg                           !<
    REAL(wp) ::  bfactor                       !< solute effects
    REAL(wp) ::  ddenom                        !<
    REAL(wp) ::  delta_r                       !<
    REAL(wp) ::  diameter                      !< diameter of cloud droplets
    REAL(wp) ::  diff_coeff                    !< diffusivity for water vapor
    REAL(wp) ::  drdt                          !<
    REAL(wp) ::  dt_ros                        !<
    REAL(wp) ::  dt_ros_sum                    !<
    REAL(wp) ::  d2rdtdr                       !<
    REAL(wp) ::  e_a                           !< current vapor pressure
    REAL(wp) ::  e_s                           !< current saturation vapor pressure
    REAL(wp) ::  error                         !< local truncation error in Rosenbrock
    REAL(wp) ::  k1                            !<
    REAL(wp) ::  k2                            !<
    REAL(wp) ::  r_err                         !< First order estimate of Rosenbrock radius
    REAL(wp) ::  r_ros                         !< Rosenbrock radius
    REAL(wp) ::  r_ros_ini                     !< initial Rosenbrock radius
    REAL(wp) ::  r0                            !< gas-kinetic lengthscale
    REAL(wp) ::  sigma                         !< surface tension of water
    REAL(wp) ::  thermal_conductivity          !< thermal conductivity for water
    REAL(wp) ::  t_int                         !< temperature
    REAL(wp) ::  w_s                           !< terminal velocity of droplets
    REAL(wp) ::  re_p                          !< particle Reynolds number
!
!-- Parameters for Rosenbrock method (see Verwer et al., 1999)
    REAL(wp), PARAMETER ::  prec = 1.0E-3_wp     !< precision of Rosenbrock solution
    REAL(wp), PARAMETER ::  q_increase = 1.5_wp  !< increase factor in timestep
    REAL(wp), PARAMETER ::  q_decrease = 0.9_wp  !< decrease factor in timestep
    REAL(wp), PARAMETER ::  gamma = 0.292893218814_wp !< = 1.0 - 1.0 / SQRT(2.0)
!
!-- Parameters for terminal velocity
    REAL(wp), PARAMETER ::  a_rog = 9.65_wp      !< parameter for fall velocity
    REAL(wp), PARAMETER ::  b_rog = 10.43_wp     !< parameter for fall velocity
    REAL(wp), PARAMETER ::  c_rog = 0.6_wp       !< parameter for fall velocity
    REAL(wp), PARAMETER ::  k_cap_rog = 4.0_wp   !< parameter for fall velocity
    REAL(wp), PARAMETER ::  k_low_rog = 12.0_wp  !< parameter for fall velocity
    REAL(wp), PARAMETER ::  d0_rog = 0.745_wp    !< separation diameter

    REAL(wp), DIMENSION(number_of_particles) ::  ventilation_effect     !<
    REAL(wp), DIMENSION(number_of_particles) ::  new_r                  !<

    CALL cpu_log( log_point_s(42), 'lpm_droplet_condens', 'start' )

!
!-- Absolute temperature
    t_int = pt(k,j,i) * exner(k)
!
!-- Saturation vapor pressure (Eq. 10 in Bolton, 1980)
    e_s = magnus( t_int )
!
!-- Current vapor pressure
    e_a = q(k,j,i) * hyp(k) / ( q(k,j,i) + rd_d_rv )
!
!-- Thermal conductivity for water (from Rogers and Yau, Table 7.1)
    thermal_conductivity = 7.94048E-05_wp * t_int + 0.00227011_wp
!
!-- Moldecular diffusivity of water vapor in air (Hall und Pruppacher, 1976)
    diff_coeff           = 0.211E-4_wp * ( t_int / 273.15_wp )**1.94_wp * &
                           ( 101325.0_wp / hyp(k) )
!
!-- Lengthscale for gas-kinetic effects (from Mordy, 1959, p. 23):
    r0 = diff_coeff / 0.036_wp * SQRT( 2.0_wp * pi / ( r_v * t_int ) )
!
!-- Calculate effects of heat conductivity and diffusion of water vapor on the
!-- diffusional growth process (usually known as 1.0 / (F_k + F_d) )
    ddenom  = 1.0_wp / ( rho_l * r_v * t_int / ( e_s * diff_coeff ) +          &
                         ( l_v / ( r_v * t_int ) - 1.0_wp ) * rho_l *          &
                         l_v / ( thermal_conductivity * t_int )                &
                       )
    new_r = 0.0_wp
!
!-- Determine ventilation effect on evaporation of large drops
    DO  n = 1, number_of_particles

       IF ( particles(n)%radius >= 4.0E-5_wp  .AND.  e_a / e_s < 1.0_wp )  THEN
!
!--       Terminal velocity is computed for vertical direction (Rogers et al.,
!--       1993, J. Appl. Meteorol.)
          diameter = particles(n)%radius * 2000.0_wp !diameter in mm
          IF ( diameter <= d0_rog )  THEN
             w_s = k_cap_rog * diameter * ( 1.0_wp - EXP( -k_low_rog * diameter ) )
          ELSE
             w_s = a_rog - b_rog * EXP( -c_rog * diameter )
          ENDIF
!
!--       Calculate droplet's Reynolds number
          re_p = 2.0_wp * particles(n)%radius * w_s / molecular_viscosity
!
!--       Ventilation coefficient (Rogers and Yau, 1989):
          IF ( re_p > 2.5_wp )  THEN
             ventilation_effect(n) = 0.78_wp + 0.28_wp * SQRT( re_p )
          ELSE
             ventilation_effect(n) = 1.0_wp + 0.09_wp * re_p
          ENDIF
       ELSE
!
!--       For small droplets or in supersaturated environments, the ventilation
!--       effect does not play a role
          ventilation_effect(n) = 1.0_wp
       ENDIF
    ENDDO

    IF( .NOT. curvature_solution_effects )  THEN
!
!--    Use analytic model for diffusional growth including gas-kinetic
!--    effects (Mordy, 1959) but without the impact of aerosols.
       DO  n = 1, number_of_particles
          arg      = ( particles(n)%radius + r0 )**2 + 2.0_wp * dt_3d * ddenom * &
                                                       ventilation_effect(n) *   &
                                                       ( e_a / e_s - 1.0_wp )
          arg      = MAX( arg, ( 0.01E-6 + r0 )**2 )
          new_r(n) = SQRT( arg ) - r0
       ENDDO

    ELSE
!
!--    Integrate the diffusional growth including gas-kinetic (Mordy, 1959),
!--    as well as curvature and solute effects (e.g., Köhler, 1936).
!
!--    Curvature effect (afactor) with surface tension (sigma) by Straka (2009)
       sigma = 0.0761_wp - 0.000155_wp * ( t_int - 273.15_wp )
!
!--    Solute effect (afactor)
       afactor = 2.0_wp * sigma / ( rho_l * r_v * t_int )

       DO  n = 1, number_of_particles
!
!--       Solute effect (bfactor)
          bfactor = vanthoff * rho_s * particles(n)%aux1**3 *                    &
                    molecular_weight_of_water / ( rho_l * molecular_weight_of_solute )

          dt_ros     = particles(n)%aux2  ! use previously stored Rosenbrock timestep
          dt_ros_sum = 0.0_wp

          r_ros     = particles(n)%radius  ! initialize Rosenbrock particle radius
          r_ros_ini = r_ros
!
!--       Integrate growth equation using a 2nd-order Rosenbrock method
!--       (see Verwer et al., 1999, Eq. (3.2)). The Rosenbrock method adjusts
!--       its with internal timestep to minimize the local truncation error.
          DO WHILE ( dt_ros_sum < dt_3d )

             dt_ros = MIN( dt_ros, dt_3d - dt_ros_sum )

             DO

                drdt = ddenom * ventilation_effect(n) * ( e_a / e_s - 1.0_wp - &
                                                          afactor / r_ros +    &
                                                          bfactor / r_ros**3   &
                                                        ) / ( r_ros + r0 )

                d2rdtdr = -ddenom * ventilation_effect(n) * (                  &
                                            (e_a / e_s - 1.0_wp ) * r_ros**4 - &
                                            afactor * r0 * r_ros**2 -          &
                                            2.0_wp * afactor * r_ros**3 +      &
                                            3.0_wp * bfactor * r0 +            &
                                            4.0_wp * bfactor * r_ros           &
                                                            )                  &
                          / ( r_ros**4 * ( r_ros + r0 )**2 )

                k1    = drdt / ( 1.0_wp - gamma * dt_ros * d2rdtdr )

                r_ros = MAX(r_ros_ini + k1 * dt_ros, particles(n)%aux1)
                r_err = r_ros

                drdt  = ddenom * ventilation_effect(n) * ( e_a / e_s - 1.0_wp - &
                                                           afactor / r_ros +    &
                                                           bfactor / r_ros**3   &
                                                         ) / ( r_ros + r0 )

                k2 = ( drdt - dt_ros * 2.0 * gamma * d2rdtdr * k1 ) / &
                     ( 1.0_wp - dt_ros * gamma * d2rdtdr )

                r_ros = MAX(r_ros_ini + dt_ros * ( 1.5_wp * k1 + 0.5_wp * k2), particles(n)%aux1)
   !
   !--          Check error of the solution, and reduce dt_ros if necessary.
                error = ABS(r_err - r_ros) / r_ros
                IF ( error > prec )  THEN
                   dt_ros = SQRT( q_decrease * prec / error ) * dt_ros
                   r_ros  = r_ros_ini
                ELSE
                   dt_ros_sum = dt_ros_sum + dt_ros
                   dt_ros     = q_increase * dt_ros
                   r_ros_ini  = r_ros
                   EXIT
                ENDIF

             END DO

          END DO !Rosenbrock loop
!
!--       Store new particle radius
          new_r(n) = r_ros
!
!--       Store internal time step value for next PALM step
          particles(n)%aux2 = dt_ros

       ENDDO !Particle loop

    ENDIF

    DO  n = 1, number_of_particles
!
!--    Sum up the change in liquid water for the respective grid
!--    box for the computation of the release/depletion of water vapor
!--    and heat.
       ql_c(k,j,i) = ql_c(k,j,i) + particles(n)%weight_factor *          &
                                   rho_l * 1.33333333_wp * pi *                &
                                   ( new_r(n)**3 - particles(n)%radius**3 ) /  &
                                   ( rho_surface * dx * dy * dzw(k) )
!
!--    Check if the increase in liqid water is not too big. If this is the case,
!--    the model timestep might be too long.
       IF ( ql_c(k,j,i) > 100.0_wp )  THEN
          WRITE( message_string, * ) 'k=',k,' j=',j,' i=',i,                &
                       ' ql_c=',ql_c(k,j,i), '&part(',n,')%wf=',            &
                       particles(n)%weight_factor,' delta_r=',delta_r
          CALL message( 'lpm_droplet_condensation', 'PA0143', 2, 2, -1, 6, 1 )
       ENDIF
!
!--    Check if the change in the droplet radius is not too big. If this is the
!--    case, the model timestep might be too long.
       delta_r = new_r(n) - particles(n)%radius
       IF ( delta_r < 0.0_wp  .AND.  new_r(n) < 0.0_wp )  THEN
          WRITE( message_string, * ) '#1 k=',k,' j=',j,' i=',i,                &
                       ' e_s=',e_s, ' e_a=',e_a,' t_int=',t_int,               &
                       '&delta_r=',delta_r,                                    &
                       ' particle_radius=',particles(n)%radius
          CALL message( 'lpm_droplet_condensation', 'PA0144', 2, 2, -1, 6, 1 )
       ENDIF
!
!--    Sum up the total volume of liquid water (needed below for
!--    re-calculating the weighting factors)
       ql_v(k,j,i) = ql_v(k,j,i) + particles(n)%weight_factor * new_r(n)**3
!
!--    Determine radius class of the particle needed for collision
       IF ( use_kernel_tables )  THEN
          particles(n)%class = ( LOG( new_r(n) ) - rclass_lbound ) /           &
                               ( rclass_ubound - rclass_lbound ) *             &
                               radius_classes
          particles(n)%class = MIN( particles(n)%class, radius_classes )
          particles(n)%class = MAX( particles(n)%class, 1 )
       ENDIF
!
!--    Store new radius to particle features
       particles(n)%radius = new_r(n)

    ENDDO

    CALL cpu_log( log_point_s(42), 'lpm_droplet_condens', 'stop' )


 END SUBROUTINE lpm_droplet_condensation


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Release of latent heat and change of mixing ratio due to condensation /
!> evaporation of droplets.
!------------------------------------------------------------------------------!
 SUBROUTINE lpm_interaction_droplets_ptq

    INTEGER(iwp) ::  i    !< running index x direction
    INTEGER(iwp) ::  j    !< running index y direction
    INTEGER(iwp) ::  k    !< running index z direction

    REAL(wp) ::  flag     !< flag to mask topography grid points

    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt
!
!--          Predetermine flag to mask topography
             flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )

             q(k,j,i)  = q_p(k,j,i)  - ql_c(k,j,i) * flag
             pt(k,j,i) = pt(k,j,i) + lv_d_cp * ql_c(k,j,i) * d_exner(k) &
                                                           * flag
          ENDDO
       ENDDO
    ENDDO

 END SUBROUTINE lpm_interaction_droplets_ptq


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Release of latent heat and change of mixing ratio due to condensation /
!> evaporation of droplets. Call for grid point i,j
!------------------------------------------------------------------------------!
 SUBROUTINE lpm_interaction_droplets_ptq_ij( i, j )

    INTEGER(iwp) ::  i    !< running index x direction
    INTEGER(iwp) ::  j    !< running index y direction
    INTEGER(iwp) ::  k    !< running index z direction

    REAL(wp) ::  flag     !< flag to mask topography grid points


    DO  k = nzb+1, nzt
!
!--    Predetermine flag to mask topography
       flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )

       q(k,j,i)  = q(k,j,i)  - ql_c(k,j,i) * flag
       pt(k,j,i) = pt(k,j,i) + lv_d_cp * ql_c(k,j,i) * d_exner(k) * flag
    ENDDO

 END SUBROUTINE lpm_interaction_droplets_ptq_ij


!------------------------------------------------------------------------------! 
! Description:
! ------------
!> Calculate the liquid water content for each grid box.
!------------------------------------------------------------------------------!
 SUBROUTINE lpm_calc_liquid_water_content


    INTEGER(iwp) ::  i   !<
    INTEGER(iwp) ::  j   !<
    INTEGER(iwp) ::  k   !<
    INTEGER(iwp) ::  n   !<

    CALL cpu_log( log_point_s(45), 'lpm_calc_ql', 'start' )

!
!-- Set water content initially to zero
    ql = 0.0_wp;  ql_v = 0.0_wp;  ql_vp = 0.0_wp

!
!-- Calculate for each grid box
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt
             number_of_particles = prt_count(k,j,i)
             IF ( number_of_particles <= 0 )  CYCLE
             particles => grid_particles(k,j,i)%particles(1:number_of_particles)
!
!--          Calculate the total volume in the boxes (ql_v, weighting factor
!--          has to beincluded)
             DO  n = 1, prt_count(k,j,i)
                ql_v(k,j,i)  = ql_v(k,j,i)  + particles(n)%weight_factor *     &
                                              particles(n)%radius**3
             ENDDO
!
!--          Calculate the liquid water content
             IF ( ql_v(k,j,i) /= 0.0_wp )  THEN
                ql(k,j,i) = ql(k,j,i) + rho_l * 1.33333333_wp * pi *           &
                                        ql_v(k,j,i) /                          &
                                        ( rho_surface * dx * dy * dzw(k) )
                IF ( ql(k,j,i) < 0.0_wp )  THEN
                   WRITE( message_string, * )  'LWC out of range: ' , &
                                               ql(k,j,i),i,j,k
                   CALL message( 'lpm_calc_liquid_water_content', '', 2, 2,    &
                                 -1, 6, 1 )
                ENDIF
             ELSE
                ql(k,j,i) = 0.0_wp
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    CALL cpu_log( log_point_s(45), 'lpm_calc_ql', 'stop' )

 END SUBROUTINE lpm_calc_liquid_water_content


!------------------------------------------------------------------------------! 
! Description:
! ------------
!> Calculates change in droplet radius by collision. Droplet collision is
!> calculated for each grid box seperately. Collision is parameterized by
!> using collision kernels. Two different kernels are available:
!> Hall kernel: Kernel from Hall (1980, J. Atmos. Sci., 2486-2507), which
!>              considers collision due to pure gravitational effects.
!> Wang kernel: Beside gravitational effects (treated with the Hall-kernel) also
!>              the effects of turbulence on the collision are considered using
!>              parameterizations of Ayala et al. (2008, New J. Phys., 10,
!>              075015) and Wang and Grabowski (2009, Atmos. Sci. Lett., 10,
!>              1-8). This kernel includes three possible effects of turbulence:
!>              the modification of the relative velocity between the droplets,
!>              the effect of preferential concentration, and the enhancement of
!>              collision efficiencies.
!------------------------------------------------------------------------------! 
 SUBROUTINE lpm_droplet_collision (i,j,k)

    INTEGER(iwp), INTENT(IN) ::  i        !<
    INTEGER(iwp), INTENT(IN) ::  j        !<
    INTEGER(iwp), INTENT(IN) ::  k        !<

    INTEGER(iwp) ::  eclass   !<
    INTEGER(iwp) ::  n        !<
    INTEGER(iwp) ::  m        !<
    INTEGER(iwp) ::  rclass_l !<
    INTEGER(iwp) ::  rclass_s !<

    REAL(wp) ::  collection_probability  !< probability for collection
    REAL(wp) ::  ddV                     !< inverse grid box volume
    REAL(wp) ::  epsilon_collision       !< dissipation rate
    REAL(wp) ::  factor_volume_to_mass   !< 4.0 / 3.0 * pi * rho_l
    REAL(wp) ::  xm                      !< droplet mass of super-droplet m
    REAL(wp) ::  xn                      !< droplet mass of super-droplet n
    REAL(wp) ::  xsm                     !< aerosol mass of super-droplet m
    REAL(wp) ::  xsn                     !< aerosol mass of super-droplet n

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  weight    !< weighting factor
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  mass      !< total mass of super droplet
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  aero_mass !< total aerosol mass of super droplet

    CALL cpu_log( log_point_s(43), 'lpm_droplet_coll', 'start' )

    number_of_particles   = prt_count(k,j,i)
    factor_volume_to_mass = 4.0_wp / 3.0_wp * pi * rho_l
    ddV                   = 1.0_wp / ( dx * dy * dzw(k) )
!
!-- Collision requires at least one super droplet inside the box
    IF ( number_of_particles > 0 )  THEN

       IF ( use_kernel_tables )  THEN
!
!--       Fast method with pre-calculated collection kernels for
!--       discrete radius- and dissipation-classes.
          IF ( wang_kernel )  THEN
             eclass = INT( diss(k,j,i) * 1.0E4_wp / 600.0_wp * &
                           dissipation_classes ) + 1
             epsilon_collision = diss(k,j,i)
          ELSE
             epsilon_collision = 0.0_wp
          ENDIF

          IF ( hall_kernel  .OR.  epsilon_collision * 1.0E4_wp < 0.001_wp )  THEN
             eclass = 0   ! Hall kernel is used
          ELSE
             eclass = MIN( dissipation_classes, eclass )
          ENDIF

       ELSE
!
!--       Collection kernels are re-calculated for every new
!--       grid box. First, allocate memory for kernel table.
!--       Third dimension is 1, because table is re-calculated for
!--       every new dissipation value.
          ALLOCATE( ckernel(1:number_of_particles,1:number_of_particles,1:1) )
!
!--       Now calculate collection kernel for this box. Note that
!--       the kernel is based on the previous time step
          CALL recalculate_kernel( i, j, k )

       ENDIF
!
!--    Temporary fields for total mass of super-droplet, aerosol mass, and
!--    weighting factor are allocated.
       ALLOCATE(mass(1:number_of_particles), weight(1:number_of_particles))
       IF ( curvature_solution_effects )  ALLOCATE(aero_mass(1:number_of_particles))

       mass(1:number_of_particles)   = particles(1:number_of_particles)%weight_factor * &
                                       particles(1:number_of_particles)%radius**3     * &
                                       factor_volume_to_mass

       weight(1:number_of_particles) = particles(1:number_of_particles)%weight_factor

       IF ( curvature_solution_effects )  THEN
          aero_mass(1:number_of_particles) = particles(1:number_of_particles)%weight_factor * &
                                             particles(1:number_of_particles)%aux1**3       * &
                                             4.0_wp / 3.0_wp * pi * rho_s
       ENDIF
!
!--    Calculate collision/coalescence
       DO  n = 1, number_of_particles

          DO  m = n, number_of_particles
!
!--          For collisions, the weighting factor of at least one super-droplet
!--          needs to be larger or equal to one.
             IF ( MIN( weight(n), weight(m) ) < 1.0_wp )  CYCLE
!
!--          Get mass of individual droplets (aerosols)
             xn = mass(n) / weight(n)
             xm = mass(m) / weight(m)
             IF ( curvature_solution_effects )  THEN
                xsn = aero_mass(n) / weight(n)
                xsm = aero_mass(m) / weight(m)
             ENDIF
!
!--          Probability that the necessary collisions take place
             IF ( use_kernel_tables )  THEN
                rclass_l = particles(n)%class
                rclass_s = particles(m)%class

                collection_probability  = MAX( weight(n), weight(m) ) *     &
                                          ckernel(rclass_l,rclass_s,eclass) * ddV * dt_3d
             ELSE
                collection_probability  = MAX( weight(n), weight(m) ) *     &
                                          ckernel(n,m,1) * ddV * dt_3d
             ENDIF
!
!--          Calculate the number of collections and consider multiple collections.
!--          (Accordingly, p_crit will be 0.0, 1.0, 2.0, ...)
             IF ( collection_probability - FLOOR(collection_probability)    &
                  > random_function( iran_part ) )  THEN
                collection_probability = FLOOR(collection_probability) + 1.0_wp
             ELSE
                collection_probability = FLOOR(collection_probability)
             ENDIF

             IF ( collection_probability > 0.0_wp )  THEN
!
!--             Super-droplet n collects droplets of super-droplet m
                IF ( weight(n) < weight(m) )  THEN

                   mass(n)   = mass(n)   + weight(n) * xm * collection_probability
                   weight(m) = weight(m) - weight(n)      * collection_probability
                   mass(m)   = mass(m)   - weight(n) * xm * collection_probability
                   IF ( curvature_solution_effects )  THEN
                      aero_mass(n) = aero_mass(n) + weight(n) * xsm * collection_probability
                      aero_mass(m) = aero_mass(m) - weight(n) * xsm * collection_probability
                   ENDIF

                ELSEIF ( weight(m) < weight(n) )  THEN

                   mass(m)   = mass(m)   + weight(m) * xn * collection_probability
                   weight(n) = weight(n) - weight(m)      * collection_probability
                   mass(n)   = mass(n)   - weight(m) * xn * collection_probability
                   IF ( curvature_solution_effects )  THEN
                      aero_mass(m) = aero_mass(m) + weight(m) * xsn * collection_probability
                      aero_mass(n) = aero_mass(n) - weight(m) * xsn * collection_probability
                   ENDIF

                ELSE
!
!--                Collisions of particles of the same weighting factor.
!--                Particle n collects 1/2 weight(n) droplets of particle m,
!--                particle m collects 1/2 weight(m) droplets of particle n.
!--                The total mass mass changes accordingly.
!--                If n = m, the first half of the droplets coalesces with the
!--                second half of the droplets; mass is unchanged because
!--                xm = xn for n = m.
!--
!--                Note: For m = n this equation is an approximation only
!--                valid for weight >> 1 (which is usually the case). The
!--                approximation is weight(n)-1 = weight(n).
                   mass(n)   = mass(n)   + 0.5_wp * weight(n) * ( xm - xn )
                   mass(m)   = mass(m)   + 0.5_wp * weight(m) * ( xn - xm )
                   IF ( curvature_solution_effects )  THEN
                      aero_mass(n) = aero_mass(n) + 0.5_wp * weight(n) * ( xsm - xsn )
                      aero_mass(m) = aero_mass(m) + 0.5_wp * weight(m) * ( xsn - xsm )
                   ENDIF
                   weight(n) = weight(n) - 0.5_wp * weight(m)
                   weight(m) = weight(n)

                ENDIF

             ENDIF

          ENDDO

          ql_vp(k,j,i) = ql_vp(k,j,i) + mass(n) / factor_volume_to_mass

       ENDDO

       IF ( ANY(weight < 0.0_wp) )  THEN
             WRITE( message_string, * ) 'negative weighting factor'
             CALL message( 'lpm_droplet_collision', 'PA0028',      &
                            2, 2, -1, 6, 1 )
       ENDIF

       particles(1:number_of_particles)%radius = ( mass(1:number_of_particles) /   &
                                                   ( weight(1:number_of_particles) &
                                                     * factor_volume_to_mass       &
                                                   )                               &
                                                 )**0.33333333333333_wp

       IF ( curvature_solution_effects )  THEN
          particles(1:number_of_particles)%aux1 = ( aero_mass(1:number_of_particles) / &
                                                    ( weight(1:number_of_particles)    &
                                                      * 4.0_wp / 3.0_wp * pi * rho_s   &
                                                    )                                  &
                                                  )**0.33333333333333_wp
       ENDIF

       particles(1:number_of_particles)%weight_factor = weight(1:number_of_particles)

       DEALLOCATE( weight, mass )
       IF ( curvature_solution_effects )  DEALLOCATE( aero_mass )
       IF ( .NOT. use_kernel_tables )  DEALLOCATE( ckernel )

!
!--    Check if LWC is conserved during collision process
       IF ( ql_v(k,j,i) /= 0.0_wp )  THEN
          IF ( ql_vp(k,j,i) / ql_v(k,j,i) >= 1.0001_wp  .OR.                      &
               ql_vp(k,j,i) / ql_v(k,j,i) <= 0.9999_wp )  THEN
             WRITE( message_string, * ) ' LWC is not conserved during',           &
                                        ' collision! ',                           &
                                        ' LWC after condensation: ', ql_v(k,j,i), &
                                        ' LWC after collision: ', ql_vp(k,j,i)
             CALL message( 'lpm_droplet_collision', 'PA0040', 2, 2, -1, 6, 1 )
          ENDIF
       ENDIF

    ENDIF

    CALL cpu_log( log_point_s(43), 'lpm_droplet_coll', 'stop' )

 END SUBROUTINE lpm_droplet_collision
 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of the collision efficiency matrix with fixed radius and
!> dissipation classes, calculated at simulation start only.
!------------------------------------------------------------------------------! 
 SUBROUTINE lpm_init_kernels

    INTEGER(iwp) ::  i !<
    INTEGER(iwp) ::  j !<
    INTEGER(iwp) ::  k !<
    
!
!-- Calculate collision efficiencies for fixed radius- and dissipation
!-- classes
    IF ( collision_kernel(6:9) == 'fast' )  THEN

       ALLOCATE( ckernel(1:radius_classes,1:radius_classes,                 &
                 0:dissipation_classes), epsclass(1:dissipation_classes),   &
                 radclass(1:radius_classes) )

!
!--    Calculate the radius class bounds with logarithmic distances
!--    in the interval [1.0E-6, 1000.0E-6] m
       rclass_lbound = LOG( 1.0E-6_wp )
       rclass_ubound = LOG( 1000.0E-6_wp )
       radclass(1)   = EXP( rclass_lbound )
       DO  i = 2, radius_classes
          radclass(i) = EXP( rclass_lbound +                                &
                             ( rclass_ubound - rclass_lbound ) *            &
                             ( i - 1.0_wp ) / ( radius_classes - 1.0_wp ) )
       ENDDO

!
!--    Set the class bounds for dissipation in interval [0.0, 600.0] cm**2/s**3
       DO  i = 1, dissipation_classes
          epsclass(i) = 0.06_wp * REAL( i, KIND=wp ) / dissipation_classes
       ENDDO
!
!--    Calculate collision efficiencies of the Wang/ayala kernel
       ALLOCATE( ec(1:radius_classes,1:radius_classes),  &
                 ecf(1:radius_classes,1:radius_classes), &
                 gck(1:radius_classes,1:radius_classes), &
                 winf(1:radius_classes) )

       DO  k = 1, dissipation_classes

          epsilon_collision = epsclass(k)
          urms    = 2.02_wp * ( epsilon_collision / 0.04_wp )**( 1.0_wp / 3.0_wp )

          CALL turbsd
          CALL turb_enhance_eff
          CALL effic

          DO  j = 1, radius_classes
             DO  i = 1, radius_classes
                ckernel(i,j,k) = ec(i,j) * gck(i,j) * ecf(i,j)
             ENDDO
          ENDDO

       ENDDO

!
!--    Calculate collision efficiencies of the Hall kernel
       ALLOCATE( hkernel(1:radius_classes,1:radius_classes), &
                 hwratio(1:radius_classes,1:radius_classes) )

       CALL fallg
       CALL effic

       DO  j = 1, radius_classes
          DO  i =  1, radius_classes
             hkernel(i,j) = pi * ( radclass(j) + radclass(i) )**2 &
                               * ec(i,j) * ABS( winf(j) - winf(i) )
             ckernel(i,j,0) = hkernel(i,j)  ! hall kernel stored on index 0
           ENDDO
       ENDDO

!
!--    Test output of efficiencies
       IF ( j == -1 )  THEN
          PRINT*, '*** Hall kernel'
          WRITE ( *,'(5X,20(F4.0,1X))' ) ( radclass(i)*1.0E6_wp, &
                                           i = 1,radius_classes )
          DO  j = 1, radius_classes
             WRITE ( *,'(F4.0,1X,20(F8.4,1X))' ) radclass(j),  &
                                       ( hkernel(i,j), i = 1,radius_classes )
          ENDDO

          DO  k = 1, dissipation_classes
             DO  i = 1, radius_classes
                DO  j = 1, radius_classes
                   IF ( hkernel(i,j) == 0.0_wp )  THEN
                      hwratio(i,j) = 9999999.9_wp
                   ELSE
                      hwratio(i,j) = ckernel(i,j,k) / hkernel(i,j)
                   ENDIF
                ENDDO
             ENDDO

             PRINT*, '*** epsilon = ', epsclass(k)
             WRITE ( *,'(5X,20(F4.0,1X))' ) ( radclass(i) * 1.0E6_wp, &
                                              i = 1,radius_classes )
             DO  j = 1, radius_classes
                WRITE ( *,'(F4.0,1X,20(F8.4,1X))' ) radclass(j) * 1.0E6_wp, &
                                       ( hwratio(i,j), i = 1,radius_classes )
             ENDDO
          ENDDO
       ENDIF

       DEALLOCATE( ec, ecf, epsclass, gck, hkernel, winf )

    ENDIF

 END SUBROUTINE lpm_init_kernels
 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculation of collision kernels during each timestep and for each grid box
!------------------------------------------------------------------------------!
 SUBROUTINE recalculate_kernel( i1, j1, k1 )


    INTEGER(iwp) ::  i      !<
    INTEGER(iwp) ::  i1     !<
    INTEGER(iwp) ::  j      !<
    INTEGER(iwp) ::  j1     !<
    INTEGER(iwp) ::  k1     !<


    number_of_particles = prt_count(k1,j1,i1)
    radius_classes      = number_of_particles   ! necessary to use the same
                                                ! subroutines as for 
                                                ! precalculated kernels

    ALLOCATE( ec(1:number_of_particles,1:number_of_particles), &
              radclass(1:number_of_particles), winf(1:number_of_particles) )

!
!-- Store particle radii on the radclass array
    radclass(1:number_of_particles) = particles(1:number_of_particles)%radius

    IF ( wang_kernel )  THEN
       epsilon_collision = diss(k1,j1,i1)   ! dissipation rate in m**2/s**3
    ELSE
       epsilon_collision = 0.0_wp
    ENDIF
    urms    = 2.02_wp * ( epsilon_collision / 0.04_wp )**( 0.33333333333_wp )

    IF ( wang_kernel  .AND.  epsilon_collision > 1.0E-7_wp )  THEN
!
!--    Call routines to calculate efficiencies for the Wang kernel
       ALLOCATE( gck(1:number_of_particles,1:number_of_particles), &
                 ecf(1:number_of_particles,1:number_of_particles) )

       CALL turbsd
       CALL turb_enhance_eff
       CALL effic

       DO  j = 1, number_of_particles
          DO  i =  1, number_of_particles
             ckernel(1+i-1,1+j-1,1) = ec(i,j) * gck(i,j) * ecf(i,j)
          ENDDO
       ENDDO

       DEALLOCATE( gck, ecf )
    ELSE
!
!--    Call routines to calculate efficiencies for the Hall kernel
       CALL fallg
       CALL effic

       DO  j = 1, number_of_particles
          DO  i =  1, number_of_particles
             ckernel(i,j,1) = pi * ( radclass(j) + radclass(i) )**2         &
                                 * ec(i,j) * ABS( winf(j) - winf(i) )
          ENDDO
       ENDDO
    ENDIF

    DEALLOCATE( ec, radclass, winf )

 END SUBROUTINE recalculate_kernel

!------------------------------------------------------------------------------! 
! Description:
! ------------
!> Calculation of effects of turbulence on the geometric collision kernel 
!> (by including the droplets' average radial relative velocities and their 
!> radial distribution function) following the analytic model by Aayala et al. 
!> (2008, New J. Phys.). For details check the second part 2 of the publication,
!> page 37ff.
!>
!> Input parameters, which need to be replaced by PALM parameters:
!>    water density, air density
!------------------------------------------------------------------------------!
 SUBROUTINE turbsd

    INTEGER(iwp) ::  i     !<
    INTEGER(iwp) ::  j     !<

    REAL(wp) ::  ao        !<
    REAL(wp) ::  ao_gr     !<
    REAL(wp) ::  bbb       !<
    REAL(wp) ::  be        !<
    REAL(wp) ::  b1        !<
    REAL(wp) ::  b2        !<
    REAL(wp) ::  ccc       !<
    REAL(wp) ::  c1        !<
    REAL(wp) ::  c1_gr     !<
    REAL(wp) ::  c2        !<
    REAL(wp) ::  d1        !<
    REAL(wp) ::  d2        !<
    REAL(wp) ::  eta       !<
    REAL(wp) ::  e1        !<
    REAL(wp) ::  e2        !<
    REAL(wp) ::  fao_gr    !<
    REAL(wp) ::  fr        !<
    REAL(wp) ::  grfin     !<
    REAL(wp) ::  lambda    !<
    REAL(wp) ::  lambda_re !<
    REAL(wp) ::  lf        !<
    REAL(wp) ::  rc        !<
    REAL(wp) ::  rrp       !<
    REAL(wp) ::  sst       !<
    REAL(wp) ::  tauk      !<
    REAL(wp) ::  tl        !<
    REAL(wp) ::  t2        !<
    REAL(wp) ::  tt        !<
    REAL(wp) ::  t1        !<
    REAL(wp) ::  vk        !<
    REAL(wp) ::  vrms1xy   !<
    REAL(wp) ::  vrms2xy   !<
    REAL(wp) ::  v1        !<
    REAL(wp) ::  v1v2xy    !<
    REAL(wp) ::  v1xysq    !<
    REAL(wp) ::  v2        !<
    REAL(wp) ::  v2xysq    !<
    REAL(wp) ::  wrfin     !<
    REAL(wp) ::  wrgrav2   !<
    REAL(wp) ::  wrtur2xy  !<
    REAL(wp) ::  xx        !<
    REAL(wp) ::  yy        !<
    REAL(wp) ::  z         !<

    REAL(wp), DIMENSION(1:radius_classes) ::  st  !< Stokes number
    REAL(wp), DIMENSION(1:radius_classes) ::  tau !< inertial time scale

    lambda    = urms * SQRT( 15.0_wp * molecular_viscosity / epsilon_collision )
    lambda_re = urms**2 * SQRT( 15.0_wp / epsilon_collision / molecular_viscosity )
    tl        = urms**2 / epsilon_collision
    lf        = 0.5_wp * urms**3 / epsilon_collision
    tauk      = SQRT( molecular_viscosity / epsilon_collision )
    eta       = ( molecular_viscosity**3 / epsilon_collision )**0.25_wp
    vk        = eta / tauk

    ao = ( 11.0_wp + 7.0_wp * lambda_re ) / ( 205.0_wp + lambda_re )
    tt = SQRT( 2.0_wp * lambda_re / ( SQRT( 15.0_wp ) * ao ) ) * tauk

!
!-- Get terminal velocity of droplets
    CALL fallg

    DO  i = 1, radius_classes
       tau(i) = winf(i) / g    ! inertial time scale
       st(i)  = tau(i) / tauk  ! Stokes number
    ENDDO

!
!-- Calculate average radial relative velocity at contact (wrfin)
    z   = tt / tl
    be  = SQRT( 2.0_wp ) * lambda / lf
    bbb = SQRT( 1.0_wp - 2.0_wp * be**2 )
    d1  = ( 1.0_wp + bbb ) / ( 2.0_wp * bbb )
    e1  = lf * ( 1.0_wp + bbb ) * 0.5_wp
    d2  = ( 1.0_wp - bbb ) * 0.5_wp / bbb
    e2  = lf * ( 1.0_wp - bbb ) * 0.5_wp
    ccc = SQRT( 1.0_wp - 2.0_wp * z**2 )
    b1  = ( 1.0_wp + ccc ) * 0.5_wp / ccc
    c1  = tl * ( 1.0_wp + ccc ) * 0.5_wp
    b2  = ( 1.0_wp - ccc ) * 0.5_wp / ccc
    c2  = tl * ( 1.0_wp - ccc ) * 0.5_wp

    DO  i = 1, radius_classes

       v1 = winf(i)
       t1 = tau(i)

       DO  j = 1, i
          rrp = radclass(i) + radclass(j)
          v2  = winf(j)
          t2  = tau(j)

          v1xysq  = b1 * d1 * phi_w(c1,e1,v1,t1) - b1 * d2 * phi_w(c1,e2,v1,t1) &
                  - b2 * d1 * phi_w(c2,e1,v1,t1) + b2 * d2 * phi_w(c2,e2,v1,t1)
          v1xysq  = v1xysq * urms**2 / t1
          vrms1xy = SQRT( v1xysq )

          v2xysq  = b1 * d1 * phi_w(c1,e1,v2,t2) - b1 * d2 * phi_w(c1,e2,v2,t2) &
                  - b2 * d1 * phi_w(c2,e1,v2,t2) + b2 * d2 * phi_w(c2,e2,v2,t2)
          v2xysq  = v2xysq * urms**2 / t2
          vrms2xy = SQRT( v2xysq )

          IF ( winf(i) >= winf(j) )  THEN
             v1 = winf(i)
             t1 = tau(i)
             v2 = winf(j)
             t2 = tau(j)
          ELSE
             v1 = winf(j)
             t1 = tau(j)
             v2 = winf(i)
             t2 = tau(i)
          ENDIF

          v1v2xy   =  b1 * d1 * zhi(c1,e1,v1,t1,v2,t2) - &
                      b1 * d2 * zhi(c1,e2,v1,t1,v2,t2) - &
                      b2 * d1 * zhi(c2,e1,v1,t1,v2,t2) + &
                      b2 * d2* zhi(c2,e2,v1,t1,v2,t2)
          fr       = d1 * EXP( -rrp / e1 ) - d2 * EXP( -rrp / e2 )
          v1v2xy   = v1v2xy * fr * urms**2 / tau(i) / tau(j)
          wrtur2xy = vrms1xy**2 + vrms2xy**2 - 2.0_wp * v1v2xy
          IF ( wrtur2xy < 0.0_wp )  wrtur2xy = 0.0_wp
          wrgrav2  = pi / 8.0_wp * ( winf(j) - winf(i) )**2
          wrfin    = SQRT( ( 2.0_wp / pi ) * ( wrtur2xy + wrgrav2) )

!
!--       Calculate radial distribution function (grfin)
          IF ( st(j) > st(i) )  THEN
             sst = st(j)
          ELSE
             sst = st(i)
          ENDIF

          xx = -0.1988_wp * sst**4 + 1.5275_wp * sst**3 - 4.2942_wp *       &
                sst**2 + 5.3406_wp * sst
          IF ( xx < 0.0_wp )  xx = 0.0_wp
          yy = 0.1886_wp * EXP( 20.306_wp / lambda_re )

          c1_gr  =  xx / ( g / vk * tauk )**yy

          ao_gr  = ao + ( pi / 8.0_wp) * ( g / vk * tauk )**2
          fao_gr = 20.115_wp * SQRT( ao_gr / lambda_re )
          rc     = SQRT( fao_gr * ABS( st(j) - st(i) ) ) * eta

          grfin  = ( ( eta**2 + rc**2 ) / ( rrp**2 + rc**2) )**( c1_gr*0.5_wp )
          IF ( grfin < 1.0_wp )  grfin = 1.0_wp

!
!--       Calculate general collection kernel (without the consideration of
!--       collection efficiencies)
          gck(i,j) = 2.0_wp * pi * rrp**2 * wrfin * grfin
          gck(j,i) = gck(i,j)

       ENDDO
    ENDDO

 END SUBROUTINE turbsd

 REAL(wp) FUNCTION phi_w( a, b, vsett, tau0 )
!
!-- Function used in the Ayala et al. (2008) analytical model for turbulent
!-- effects on the collision kernel
    

    REAL(wp) ::  a     !<
    REAL(wp) ::  aa1   !<
    REAL(wp) ::  b     !<
    REAL(wp) ::  tau0  !<
    REAL(wp) ::  vsett !<

    aa1 = 1.0_wp / tau0 + 1.0_wp / a + vsett / b
    phi_w = 1.0_wp / aa1  - 0.5_wp * vsett / b / aa1**2

 END FUNCTION phi_w

 REAL(wp) FUNCTION zhi( a, b, vsett1, tau1, vsett2, tau2 )
!
!-- Function used in the Ayala et al. (2008) analytical model for turbulent
!-- effects on the collision kernel

    REAL(wp) ::  a      !<
    REAL(wp) ::  aa1    !<
    REAL(wp) ::  aa2    !<
    REAL(wp) ::  aa3    !<
    REAL(wp) ::  aa4    !<
    REAL(wp) ::  aa5    !<
    REAL(wp) ::  aa6    !<
    REAL(wp) ::  b      !<
    REAL(wp) ::  tau1   !<
    REAL(wp) ::  tau2   !<
    REAL(wp) ::  vsett1 !<
    REAL(wp) ::  vsett2 !<

    aa1 = vsett2 / b - 1.0_wp / tau2 - 1.0_wp / a
    aa2 = vsett1 / b + 1.0_wp / tau1 + 1.0_wp / a
    aa3 = ( vsett1 - vsett2 ) / b + 1.0_wp / tau1 + 1.0_wp / tau2
    aa4 = ( vsett2 / b )**2 - ( 1.0_wp / tau2 + 1.0_wp / a )**2
    aa5 = vsett2 / b + 1.0_wp / tau2 + 1.0_wp / a
    aa6 = 1.0_wp / tau1 - 1.0_wp / a + ( 1.0_wp / tau2 + 1.0_wp / a) *      &
          vsett1 / vsett2
    zhi = (1.0_wp / aa1 - 1.0_wp / aa2 ) * ( vsett1 - vsett2 ) * 0.5_wp /   &
          b / aa3**2 + ( 4.0_wp / aa4 - 1.0_wp / aa5**2 - 1.0_wp / aa1**2 ) &
          * vsett2 * 0.5_wp / b /aa6 + ( 2.0_wp * ( b / aa2 - b / aa1 ) -   &
          vsett1 / aa2**2 + vsett2 / aa1**2 ) * 0.5_wp / b / aa3

 END FUNCTION zhi


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Parameterization of terminal velocity following Rogers et al. (1993, J. Appl.
!> Meteorol.)
!------------------------------------------------------------------------------!
 SUBROUTINE fallg

    INTEGER(iwp) ::  j                            !<

    REAL(wp), PARAMETER ::  k_cap_rog = 4.0_wp    !< parameter
    REAL(wp), PARAMETER ::  k_low_rog = 12.0_wp   !< parameter
    REAL(wp), PARAMETER ::  a_rog     = 9.65_wp   !< parameter
    REAL(wp), PARAMETER ::  b_rog     = 10.43_wp  !< parameter
    REAL(wp), PARAMETER ::  c_rog     = 0.6_wp    !< parameter
    REAL(wp), PARAMETER ::  d0_rog    = 0.745_wp  !< seperation diameter

    REAL(wp)            ::  diameter              !< droplet diameter in mm


    DO  j = 1, radius_classes

       diameter = radclass(j) * 2000.0_wp

       IF ( diameter <= d0_rog )  THEN
          winf(j) = k_cap_rog * diameter * ( 1.0_wp -                       &
                                             EXP( -k_low_rog * diameter ) )
       ELSE
          winf(j) = a_rog - b_rog * EXP( -c_rog * diameter )
       ENDIF

    ENDDO

 END SUBROUTINE fallg


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Interpolation of collision efficiencies (Hall, 1980, J. Atmos. Sci.)
!------------------------------------------------------------------------------!
 SUBROUTINE effic
 
    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  iq !<
    INTEGER(iwp) ::  ir !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  k  !<

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  ira !<

    LOGICAL, SAVE ::  first = .TRUE. !<

    REAL(wp) ::  ek              !<
    REAL(wp) ::  particle_radius !<
    REAL(wp) ::  pp              !<
    REAL(wp) ::  qq              !<
    REAL(wp) ::  rq              !<

    REAL(wp), DIMENSION(1:21), SAVE ::  rat        !<
    
    REAL(wp), DIMENSION(1:15), SAVE ::  r0         !<
    
    REAL(wp), DIMENSION(1:15,1:21), SAVE ::  ecoll !<

!
!-- Initial assignment of constants
    IF ( first )  THEN

      first = .FALSE.
      r0  = (/   6.0_wp,   8.0_wp,  10.0_wp, 15.0_wp,  20.0_wp,  25.0_wp,   &
                30.0_wp,  40.0_wp,  50.0_wp, 60.0_wp,  70.0_wp, 100.0_wp,   &
               150.0_wp, 200.0_wp, 300.0_wp /)

      rat = (/ 0.00_wp, 0.05_wp, 0.10_wp, 0.15_wp, 0.20_wp, 0.25_wp,        &
               0.30_wp, 0.35_wp, 0.40_wp, 0.45_wp, 0.50_wp, 0.55_wp,        &
               0.60_wp, 0.65_wp, 0.70_wp, 0.75_wp, 0.80_wp, 0.85_wp,        &
               0.90_wp, 0.95_wp, 1.00_wp /)

      ecoll(:,1)  = (/ 0.001_wp, 0.001_wp, 0.001_wp, 0.001_wp, 0.001_wp,    &
                       0.001_wp, 0.001_wp, 0.001_wp, 0.001_wp, 0.001_wp,    &
                       0.001_wp, 0.001_wp, 0.001_wp, 0.001_wp, 0.001_wp /)
      ecoll(:,2)  = (/ 0.003_wp, 0.003_wp, 0.003_wp, 0.004_wp, 0.005_wp,    &
                       0.005_wp, 0.005_wp, 0.010_wp, 0.100_wp, 0.050_wp,    &
                       0.200_wp, 0.500_wp, 0.770_wp, 0.870_wp, 0.970_wp /)
      ecoll(:,3)  = (/ 0.007_wp, 0.007_wp, 0.007_wp, 0.008_wp, 0.009_wp,    &
                       0.010_wp, 0.010_wp, 0.070_wp, 0.400_wp, 0.430_wp,    &
                       0.580_wp, 0.790_wp, 0.930_wp, 0.960_wp, 1.000_wp /)
      ecoll(:,4)  = (/ 0.009_wp, 0.009_wp, 0.009_wp, 0.012_wp, 0.015_wp,    &
                       0.010_wp, 0.020_wp, 0.280_wp, 0.600_wp, 0.640_wp,    &
                       0.750_wp, 0.910_wp, 0.970_wp, 0.980_wp, 1.000_wp /)
      ecoll(:,5)  = (/ 0.014_wp, 0.014_wp, 0.014_wp, 0.015_wp, 0.016_wp,    &
                       0.030_wp, 0.060_wp, 0.500_wp, 0.700_wp, 0.770_wp,    &
                       0.840_wp, 0.950_wp, 0.970_wp, 1.000_wp, 1.000_wp /)
      ecoll(:,6)  = (/ 0.017_wp, 0.017_wp, 0.017_wp, 0.020_wp, 0.022_wp,    &
                       0.060_wp, 0.100_wp, 0.620_wp, 0.780_wp, 0.840_wp,    &
                       0.880_wp, 0.950_wp, 1.000_wp, 1.000_wp, 1.000_wp /)
      ecoll(:,7)  = (/ 0.030_wp, 0.030_wp, 0.024_wp, 0.022_wp, 0.032_wp,    &
                       0.062_wp, 0.200_wp, 0.680_wp, 0.830_wp, 0.870_wp,    &
                       0.900_wp, 0.950_wp, 1.000_wp, 1.000_wp, 1.000_wp /)
      ecoll(:,8)  = (/ 0.025_wp, 0.025_wp, 0.025_wp, 0.036_wp, 0.043_wp,    &
                       0.130_wp, 0.270_wp, 0.740_wp, 0.860_wp, 0.890_wp,    &
                       0.920_wp, 1.000_wp, 1.000_wp, 1.000_wp, 1.000_wp /)
      ecoll(:,9)  = (/ 0.027_wp, 0.027_wp, 0.027_wp, 0.040_wp, 0.052_wp,    &
                       0.200_wp, 0.400_wp, 0.780_wp, 0.880_wp, 0.900_wp,    &
                       0.940_wp, 1.000_wp, 1.000_wp, 1.000_wp, 1.000_wp /)
      ecoll(:,10) = (/ 0.030_wp, 0.030_wp, 0.030_wp, 0.047_wp, 0.064_wp,    &
                       0.250_wp, 0.500_wp, 0.800_wp, 0.900_wp, 0.910_wp,    &
                       0.950_wp, 1.000_wp, 1.000_wp, 1.000_wp, 1.000_wp /)
      ecoll(:,11) = (/ 0.040_wp, 0.040_wp, 0.033_wp, 0.037_wp, 0.068_wp,    &
                       0.240_wp, 0.550_wp, 0.800_wp, 0.900_wp, 0.910_wp,    &
                       0.950_wp, 1.000_wp, 1.000_wp, 1.000_wp, 1.000_wp /)
      ecoll(:,12) = (/ 0.035_wp, 0.035_wp, 0.035_wp, 0.055_wp, 0.079_wp,    &
                       0.290_wp, 0.580_wp, 0.800_wp, 0.900_wp, 0.910_wp,    &
                       0.950_wp, 1.000_wp, 1.000_wp, 1.000_wp, 1.000_wp /)
      ecoll(:,13) = (/ 0.037_wp, 0.037_wp, 0.037_wp, 0.062_wp, 0.082_wp,    &
                       0.290_wp, 0.590_wp, 0.780_wp, 0.900_wp, 0.910_wp,    &
                       0.950_wp, 1.000_wp, 1.000_wp, 1.000_wp, 1.000_wp /)
      ecoll(:,14) = (/ 0.037_wp, 0.037_wp, 0.037_wp, 0.060_wp, 0.080_wp,    &
                       0.290_wp, 0.580_wp, 0.770_wp, 0.890_wp, 0.910_wp,    &
                       0.950_wp, 1.000_wp, 1.000_wp, 1.000_wp, 1.000_wp /)
      ecoll(:,15) = (/ 0.037_wp, 0.037_wp, 0.037_wp, 0.041_wp, 0.075_wp,    &
                       0.250_wp, 0.540_wp, 0.760_wp, 0.880_wp, 0.920_wp,    &
                       0.950_wp, 1.000_wp, 1.000_wp, 1.000_wp, 1.000_wp /)
      ecoll(:,16) = (/ 0.037_wp, 0.037_wp, 0.037_wp, 0.052_wp, 0.067_wp,    &
                       0.250_wp, 0.510_wp, 0.770_wp, 0.880_wp, 0.930_wp,    &
                       0.970_wp, 1.000_wp, 1.000_wp, 1.000_wp, 1.000_wp /)
      ecoll(:,17) = (/ 0.037_wp, 0.037_wp, 0.037_wp, 0.047_wp, 0.057_wp,    &
                       0.250_wp, 0.490_wp, 0.770_wp, 0.890_wp, 0.950_wp,    &
                       1.000_wp, 1.000_wp, 1.000_wp, 1.000_wp, 1.000_wp /)
      ecoll(:,18) = (/ 0.036_wp, 0.036_wp, 0.036_wp, 0.042_wp, 0.048_wp,    &
                       0.230_wp, 0.470_wp, 0.780_wp, 0.920_wp, 1.000_wp,    &
                       1.020_wp, 1.020_wp, 1.020_wp, 1.020_wp, 1.020_wp /)
      ecoll(:,19) = (/ 0.040_wp, 0.040_wp, 0.035_wp, 0.033_wp, 0.040_wp,    &
                       0.112_wp, 0.450_wp, 0.790_wp, 1.010_wp, 1.030_wp,    &
                       1.040_wp, 1.040_wp, 1.040_wp, 1.040_wp, 1.040_wp /)
      ecoll(:,20) = (/ 0.033_wp, 0.033_wp, 0.033_wp, 0.033_wp, 0.033_wp,    &
                       0.119_wp, 0.470_wp, 0.950_wp, 1.300_wp, 1.700_wp,    &
                       2.300_wp, 2.300_wp, 2.300_wp, 2.300_wp, 2.300_wp /)
      ecoll(:,21) = (/ 0.027_wp, 0.027_wp, 0.027_wp, 0.027_wp, 0.027_wp,    &
                       0.125_wp, 0.520_wp, 1.400_wp, 2.300_wp, 3.000_wp,    &
                       4.000_wp, 4.000_wp, 4.000_wp, 4.000_wp, 4.000_wp /)
    ENDIF

!
!-- Calculate the radius class index of particles with respect to array r
!-- Radius has to be in microns
    ALLOCATE( ira(1:radius_classes) )
    DO  j = 1, radius_classes
       particle_radius = radclass(j) * 1.0E6_wp
       DO  k = 1, 15
          IF ( particle_radius < r0(k) )  THEN
             ira(j) = k
             EXIT
          ENDIF
       ENDDO
       IF ( particle_radius >= r0(15) )  ira(j) = 16
    ENDDO

!
!-- Two-dimensional linear interpolation of the collision efficiency.
!-- Radius has to be in microns
    DO  j = 1, radius_classes
       DO  i = 1, j

          ir = MAX( ira(i), ira(j) )
          rq = MIN( radclass(i) / radclass(j), radclass(j) / radclass(i) )
          iq = INT( rq * 20 ) + 1
          iq = MAX( iq , 2)

          IF ( ir < 16 )  THEN
             IF ( ir >= 2 )  THEN
                pp = ( ( MAX( radclass(j), radclass(i) ) * 1.0E6_wp ) -     &
                       r0(ir-1) ) / ( r0(ir) - r0(ir-1) )
                qq = ( rq - rat(iq-1) ) / ( rat(iq) - rat(iq-1) )
                ec(j,i) = ( 1.0_wp - pp ) * ( 1.0_wp - qq )                 &
                          * ecoll(ir-1,iq-1)                                &
                          + pp * ( 1.0_wp - qq ) * ecoll(ir,iq-1)           &
                          + qq * ( 1.0_wp - pp ) * ecoll(ir-1,iq)           &
                          + pp * qq * ecoll(ir,iq)
             ELSE
                qq = ( rq - rat(iq-1) ) / ( rat(iq) - rat(iq-1) )
                ec(j,i) = ( 1.0_wp - qq ) * ecoll(1,iq-1) + qq * ecoll(1,iq)
             ENDIF
          ELSE
             qq = ( rq - rat(iq-1) ) / ( rat(iq) - rat(iq-1) )
             ek = ( 1.0_wp - qq ) * ecoll(15,iq-1) + qq * ecoll(15,iq)
             ec(j,i) = MIN( ek, 1.0_wp )
          ENDIF

          IF ( ec(j,i) < 1.0E-20_wp )  ec(j,i) = 0.0_wp

          ec(i,j) = ec(j,i)

       ENDDO
    ENDDO

    DEALLOCATE( ira )

 END SUBROUTINE effic


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Interpolation of turbulent enhancement factor for collision efficencies
!> following Wang and Grabowski (2009, Atmos. Sci. Let.)
!------------------------------------------------------------------------------!
 SUBROUTINE turb_enhance_eff

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  iq !<
    INTEGER(iwp) ::  ir !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  k  !<
    INTEGER(iwp) ::  kk !<

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  ira !<

    LOGICAL, SAVE ::  first = .TRUE. !<

    REAL(wp) ::  particle_radius !<
    REAL(wp) ::  pp              !<
    REAL(wp) ::  qq              !<
    REAL(wp) ::  rq              !<
    REAL(wp) ::  y1              !<
    REAL(wp) ::  y2              !<
    REAL(wp) ::  y3              !<

    REAL(wp), DIMENSION(1:11), SAVE ::  rat           !<
    REAL(wp), DIMENSION(1:7), SAVE  ::  r0            !<

    REAL(wp), DIMENSION(1:7,1:11), SAVE ::  ecoll_100 !<
    REAL(wp), DIMENSION(1:7,1:11), SAVE ::  ecoll_400 !<

!
!-- Initial assignment of constants
    IF ( first )  THEN

       first = .FALSE.

       r0  = (/  10.0_wp, 20.0_wp, 30.0_wp, 40.0_wp, 50.0_wp, 60.0_wp,  &
                100.0_wp /)

       rat = (/ 0.0_wp, 0.1_wp, 0.2_wp, 0.3_wp, 0.4_wp, 0.5_wp, 0.6_wp, &
                0.7_wp, 0.8_wp, 0.9_wp, 1.0_wp /)
!
!--    Tabulated turbulent enhancement factor at 100 cm**2/s**3
       ecoll_100(:,1)  = (/  1.74_wp,   1.74_wp,   1.773_wp, 1.49_wp,  &
                             1.207_wp,  1.207_wp,  1.0_wp /)
       ecoll_100(:,2)  = (/  1.46_wp,   1.46_wp,   1.421_wp, 1.245_wp, &
                             1.069_wp,  1.069_wp,  1.0_wp /)
       ecoll_100(:,3)  = (/  1.32_wp,   1.32_wp,   1.245_wp, 1.123_wp, &
                             1.000_wp,  1.000_wp,  1.0_wp /)
       ecoll_100(:,4)  = (/  1.250_wp,  1.250_wp,  1.148_wp, 1.087_wp, &
                             1.025_wp,  1.025_wp,  1.0_wp /)
       ecoll_100(:,5)  = (/  1.186_wp,  1.186_wp,  1.066_wp, 1.060_wp, &
                             1.056_wp,  1.056_wp,  1.0_wp /)
       ecoll_100(:,6)  = (/  1.045_wp,  1.045_wp,  1.000_wp, 1.014_wp, &
                             1.028_wp,  1.028_wp,  1.0_wp /)
       ecoll_100(:,7)  = (/  1.070_wp,  1.070_wp,  1.030_wp, 1.038_wp, &
                             1.046_wp,  1.046_wp,  1.0_wp /)
       ecoll_100(:,8)  = (/  1.000_wp,  1.000_wp,  1.054_wp, 1.042_wp, &
                             1.029_wp,  1.029_wp,  1.0_wp /)
       ecoll_100(:,9)  = (/  1.223_wp,  1.223_wp,  1.117_wp, 1.069_wp, &
                             1.021_wp,  1.021_wp,  1.0_wp /)
       ecoll_100(:,10) = (/  1.570_wp,  1.570_wp,  1.244_wp, 1.166_wp, &
                             1.088_wp,  1.088_wp,  1.0_wp /)
       ecoll_100(:,11) = (/ 20.3_wp,   20.3_wp,   14.6_wp,   8.61_wp,  &
                             2.60_wp,   2.60_wp,   1.0_wp /)
!
!--    Tabulated turbulent enhancement factor at 400 cm**2/s**3
       ecoll_400(:,1)  = (/  4.976_wp,  4.976_wp,  3.593_wp,  2.519_wp, &
                             1.445_wp,  1.445_wp,  1.0_wp /)
       ecoll_400(:,2)  = (/  2.984_wp,  2.984_wp,  2.181_wp,  1.691_wp, &
                             1.201_wp,  1.201_wp,  1.0_wp /)
       ecoll_400(:,3)  = (/  1.988_wp,  1.988_wp,  1.475_wp,  1.313_wp, &
                             1.150_wp,  1.150_wp,  1.0_wp /)
       ecoll_400(:,4)  = (/  1.490_wp,  1.490_wp,  1.187_wp,  1.156_wp, &
                             1.126_wp,  1.126_wp,  1.0_wp /)
       ecoll_400(:,5)  = (/  1.249_wp,  1.249_wp,  1.088_wp,  1.090_wp, &
                             1.092_wp,  1.092_wp,  1.0_wp /)
       ecoll_400(:,6)  = (/  1.139_wp,  1.139_wp,  1.130_wp,  1.091_wp, &
                             1.051_wp,  1.051_wp,  1.0_wp /)
       ecoll_400(:,7)  = (/  1.220_wp,  1.220_wp,  1.190_wp,  1.138_wp, &
                             1.086_wp,  1.086_wp,  1.0_wp /)
       ecoll_400(:,8)  = (/  1.325_wp,  1.325_wp,  1.267_wp,  1.165_wp, &
                             1.063_wp,  1.063_wp,  1.0_wp /)
       ecoll_400(:,9)  = (/  1.716_wp,  1.716_wp,  1.345_wp,  1.223_wp, &
                             1.100_wp,  1.100_wp,  1.0_wp /)
       ecoll_400(:,10) = (/  3.788_wp,  3.788_wp,  1.501_wp,  1.311_wp, &
                             1.120_wp,  1.120_wp,  1.0_wp /)
       ecoll_400(:,11) = (/ 36.52_wp,  36.52_wp,  19.16_wp,  22.80_wp,  &
                            26.0_wp,   26.0_wp,    1.0_wp /)

    ENDIF

!
!-- Calculate the radius class index of particles with respect to array r0
!-- The droplet radius has to be given in microns.
    ALLOCATE( ira(1:radius_classes) )

    DO  j = 1, radius_classes
       particle_radius = radclass(j) * 1.0E6_wp
       DO  k = 1, 7
          IF ( particle_radius < r0(k) )  THEN
             ira(j) = k
             EXIT
          ENDIF
       ENDDO
       IF ( particle_radius >= r0(7) )  ira(j) = 8
    ENDDO

!
!-- Two-dimensional linear interpolation of the turbulent enhancement factor.
!-- The droplet radius has to be given in microns.
    DO  j =  1, radius_classes
       DO  i = 1, j

          ir = MAX( ira(i), ira(j) )
          rq = MIN( radclass(i) / radclass(j), radclass(j) / radclass(i) )

          DO  kk = 2, 11
             IF ( rq <= rat(kk) )  THEN
                iq = kk
                EXIT
             ENDIF
          ENDDO

          y1 = 1.0_wp  ! turbulent enhancement factor at 0 m**2/s**3

          IF ( ir < 8 )  THEN
             IF ( ir >= 2 )  THEN
                pp = ( MAX( radclass(j), radclass(i) ) * 1.0E6_wp -  &
                       r0(ir-1) ) / ( r0(ir) - r0(ir-1) )
                qq = ( rq - rat(iq-1) ) / ( rat(iq) - rat(iq-1) )
                y2 = ( 1.0_wp - pp ) * ( 1.0_wp - qq ) * ecoll_100(ir-1,iq-1) + &
                             pp * ( 1.0_wp - qq ) * ecoll_100(ir,iq-1)        + &
                             qq * ( 1.0_wp - pp ) * ecoll_100(ir-1,iq)        + &
                             pp * qq              * ecoll_100(ir,iq)
                y3 = ( 1.0-pp ) * ( 1.0_wp - qq ) * ecoll_400(ir-1,iq-1)      + &
                             pp * ( 1.0_wp - qq ) * ecoll_400(ir,iq-1)        + &
                             qq * ( 1.0_wp - pp ) * ecoll_400(ir-1,iq)        + &
                             pp * qq              * ecoll_400(ir,iq)
             ELSE
                qq = ( rq - rat(iq-1) ) / ( rat(iq) - rat(iq-1) )
                y2 = ( 1.0_wp - qq ) * ecoll_100(1,iq-1) + qq * ecoll_100(1,iq)
                y3 = ( 1.0_wp - qq ) * ecoll_400(1,iq-1) + qq * ecoll_400(1,iq)
             ENDIF
          ELSE
             qq = ( rq - rat(iq-1) ) / ( rat(iq) - rat(iq-1) )
             y2 = ( 1.0_wp - qq ) * ecoll_100(7,iq-1) + qq * ecoll_100(7,iq)
             y3 = ( 1.0_wp - qq ) * ecoll_400(7,iq-1) + qq * ecoll_400(7,iq)
          ENDIF
!
!--       Linear interpolation of turbulent enhancement factor
          IF ( epsilon_collision <= 0.01_wp )  THEN
             ecf(j,i) = ( epsilon_collision - 0.01_wp ) / ( 0.0_wp  - 0.01_wp ) * y1 &
                      + ( epsilon_collision - 0.0_wp  ) / ( 0.01_wp - 0.0_wp  ) * y2
          ELSEIF ( epsilon_collision <= 0.06_wp )  THEN
             ecf(j,i) = ( epsilon_collision - 0.04_wp ) / ( 0.01_wp - 0.04_wp ) * y2 &
                      + ( epsilon_collision - 0.01_wp ) / ( 0.04_wp - 0.01_wp ) * y3
          ELSE
             ecf(j,i) = ( 0.06_wp - 0.04_wp ) / ( 0.01_wp - 0.04_wp ) * y2 &
                      + ( 0.06_wp - 0.01_wp ) / ( 0.04_wp - 0.01_wp ) * y3
          ENDIF

          IF ( ecf(j,i) < 1.0_wp )  ecf(j,i) = 1.0_wp

          ecf(i,j) = ecf(j,i)

       ENDDO
    ENDDO

 END SUBROUTINE turb_enhance_eff
 
 
 !------------------------------------------------------------------------------!
! Description:
! ------------
! This routine is a part of the Lagrangian particle model. Super droplets which 
! fulfill certain criterion's (e.g. a big weighting factor and a large radius) 
! can be split into several super droplets with a reduced number of 
! represented particles of every super droplet. This mechanism ensures an
! improved representation of the right tail of the drop size distribution with 
! a feasible amount of computational costs. The limits of particle creation 
! should be chosen carefully! The idea of this algorithm is based on 
! Unterstrasser and Soelch, 2014. 
!------------------------------------------------------------------------------!
 SUBROUTINE lpm_splitting

    INTEGER(iwp) ::  i                !< 
    INTEGER(iwp) ::  j                !<
    INTEGER(iwp) ::  jpp              !<
    INTEGER(iwp) ::  k                !<
    INTEGER(iwp) ::  n                !<
    INTEGER(iwp) ::  new_particles_gb !< counter of created particles within one grid box
    INTEGER(iwp) ::  new_size         !< new particle array size
    INTEGER(iwp) ::  np               !< 
    INTEGER(iwp) ::  old_size         !< old particle array size

    INTEGER(iwp), PARAMETER ::  n_max = 100 !< number of radii bin for splitting functions    

    LOGICAL ::  first_loop_stride_sp = .TRUE. !< flag to calculate constants only once

    REAL(wp) ::  diameter                 !< diameter of droplet
    REAL(wp) ::  dlog                     !< factor for DSD calculation
    REAL(wp) ::  factor_volume_to_mass    !< pre calculate factor volume to mass
    REAL(wp) ::  lambda                   !< slope parameter of gamma-distribution
    REAL(wp) ::  lwc                      !< liquid water content of grid box
    REAL(wp) ::  lwc_total                !< average liquid water content of cloud
    REAL(wp) ::  m1                       !< first moment of DSD
    REAL(wp) ::  m1_total                 !< average over all PEs of first moment of DSD
    REAL(wp) ::  m2                       !< second moment of DSD
    REAL(wp) ::  m2_total                 !< average average over all PEs second moment of DSD
    REAL(wp) ::  m3                       !< third moment of DSD
    REAL(wp) ::  m3_total                 !< average average over all PEs third moment of DSD
    REAL(wp) ::  mu                       !< spectral shape parameter of gamma distribution
    REAL(wp) ::  nrclgb                   !< number of cloudy grid boxes (ql >= 1.0E-5 kg/kg) 
    REAL(wp) ::  nrclgb_total             !< average over all PEs of number of cloudy grid boxes
    REAL(wp) ::  nr                       !< number concentration of cloud droplets
    REAL(wp) ::  nr_total                 !< average over all PEs of number of cloudy grid boxes
    REAL(wp) ::  nr0                      !< intercept parameter of gamma distribution
    REAL(wp) ::  pirho_l                  !< pi * rho_l / 6.0
    REAL(wp) ::  ql_crit = 1.0E-5_wp      !< threshold lwc for cloudy grid cells 
                                          !< (Siebesma et al 2003, JAS, 60)
    REAL(wp) ::  rm                       !< volume averaged mean radius
    REAL(wp) ::  rm_total                 !< average over all PEs of volume averaged mean radius
    REAL(wp) ::  r_min = 1.0E-6_wp        !< minimum radius of approximated spectra 
    REAL(wp) ::  r_max = 1.0E-3_wp        !< maximum radius of approximated spectra
    REAL(wp) ::  sigma_log = 1.5_wp       !< standard deviation of the LOG-distribution
    REAL(wp) ::  zeta                     !< Parameter for DSD calculation of Seifert

    REAL(wp), DIMENSION(0:n_max-1) ::  an_spl     !< size dependent critical weight factor
    REAL(wp), DIMENSION(0:n_max-1) ::  r_bin_mid  !< mass weighted mean radius of a bin
    REAL(wp), DIMENSION(0:n_max)   ::  r_bin      !< boundaries of a radius bin

    TYPE(particle_type) ::  tmp_particle   !< temporary particle TYPE

    CALL cpu_log( log_point_s(80), 'lpm_splitting', 'start' )

    IF ( first_loop_stride_sp )  THEN
       IF ( i_splitting_mode == 2  .OR.  i_splitting_mode == 3 )  THEN
          dlog   = ( LOG10(r_max) - LOG10(r_min) ) / ( n_max - 1 )
          DO  i = 0, n_max-1
             r_bin(i) = 10.0_wp**( LOG10(r_min) + i * dlog - 0.5_wp * dlog )
             r_bin_mid(i) = 10.0_wp**( LOG10(r_min) + i * dlog )
          ENDDO
          r_bin(n_max) = 10.0_wp**( LOG10(r_min) + n_max * dlog - 0.5_wp * dlog )
       ENDIF   
       factor_volume_to_mass =  4.0_wp / 3.0_wp * pi * rho_l
       pirho_l  = pi * rho_l / 6.0_wp
       IF ( weight_factor_split == -1.0_wp )  THEN
          weight_factor_split = 0.1_wp * initial_weighting_factor 
       ENDIF
    ENDIF


    IF ( i_splitting_mode == 1 )  THEN

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt

                new_particles_gb = 0
                number_of_particles = prt_count(k,j,i)
                IF ( number_of_particles <= 0  .OR.                            &  
                     ql(k,j,i) < ql_crit )  CYCLE
                particles => grid_particles(k,j,i)%particles(1:number_of_particles)
!
!--             Start splitting operations. Each particle is checked if it
!--             fulfilled the splitting criterion's. In splitting mode 'const'   
!--             a critical radius  (radius_split) a critical weighting factor
!--             (weight_factor_split) and a splitting factor (splitting_factor)
!--             must  be prescribed (see particle_parameters). Super droplets 
!--             which have a larger radius and larger weighting factor are split 
!--             into 'splitting_factor' super droplets. Therefore, the weighting 
!--             factor of  the super droplet and all created clones is reduced 
!--             by the factor of 'splitting_factor'.
                DO  n = 1, number_of_particles
                   IF ( particles(n)%particle_mask  .AND.                      &
                        particles(n)%radius >= radius_split  .AND.             & 
                        particles(n)%weight_factor >= weight_factor_split )    &
                   THEN
!
!--                   Calculate the new number of particles.
                      new_size = prt_count(k,j,i) + splitting_factor - 1
!
!--                   Cycle if maximum number of particles per grid box
!--                   is greater than the allowed maximum number.
                      IF ( new_size >= max_number_particles_per_gridbox )  CYCLE
!
!--                   Reallocate particle array if necessary. 
                      IF ( new_size > SIZE(particles) )  THEN 
                         CALL realloc_particles_array( i, j, k, new_size )
                      ENDIF
                      old_size = prt_count(k,j,i)
!
!--                   Calculate new weighting factor.
                      particles(n)%weight_factor =  & 
                         particles(n)%weight_factor / splitting_factor
                      tmp_particle = particles(n)
!
!--                   Create splitting_factor-1 new particles.
                      DO  jpp = 1, splitting_factor-1
                         grid_particles(k,j,i)%particles(jpp+old_size) =       & 
                            tmp_particle
                      ENDDO  
                      new_particles_gb = new_particles_gb + splitting_factor - 1
!   
!--                   Save the new number of super droplets for every grid box.
                      prt_count(k,j,i) = prt_count(k,j,i) +                    &
                                         splitting_factor - 1
                   ENDIF
                ENDDO

             ENDDO
          ENDDO
       ENDDO

    ELSEIF ( i_splitting_mode == 2 )  THEN 
!
!--    Initialize summing variables.
       lwc          = 0.0_wp
       lwc_total    = 0.0_wp 
       m1           = 0.0_wp
       m1_total     = 0.0_wp
       m2           = 0.0_wp
       m2_total     = 0.0_wp
       m3           = 0.0_wp
       m3_total     = 0.0_wp
       nr           = 0.0_wp
       nrclgb       = 0.0_wp
       nrclgb_total = 0.0_wp
       nr_total     = 0.0_wp
       rm           = 0.0_wp
       rm_total     = 0.0_wp

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                number_of_particles = prt_count(k,j,i)
                IF ( number_of_particles <= 0  .OR.                            & 
                     ql(k,j,i) < ql_crit )  CYCLE
                particles => grid_particles(k,j,i)%particles(1:number_of_particles)
                nrclgb = nrclgb + 1.0_wp
!
!--             Calculate moments of DSD.
                DO  n = 1, number_of_particles
                   IF ( particles(n)%particle_mask  .AND.                      &
                        particles(n)%radius >= r_min )                         &
                   THEN
                      nr  = nr  + particles(n)%weight_factor
                      rm  = rm  + factor_volume_to_mass  *                     &
                                 particles(n)%radius**3  *                     &
                                 particles(n)%weight_factor
                      IF ( isf == 1 )  THEN           
                         diameter   = particles(n)%radius * 2.0_wp
                         lwc = lwc + factor_volume_to_mass *                   &
                                     particles(n)%radius**3 *                  & 
                                     particles(n)%weight_factor 
                         m1  = m1  + particles(n)%weight_factor * diameter
                         m2  = m2  + particles(n)%weight_factor * diameter**2
                         m3  = m3  + particles(n)%weight_factor * diameter**3
                      ENDIF
                   ENDIF
                ENDDO 
             ENDDO
          ENDDO
       ENDDO

#if defined( __parallel )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( nr, nr_total, 1 , &
       MPI_REAL, MPI_SUM, comm2d, ierr )
       CALL MPI_ALLREDUCE( rm, rm_total, 1 , &
       MPI_REAL, MPI_SUM, comm2d, ierr )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( nrclgb, nrclgb_total, 1 , &
       MPI_REAL, MPI_SUM, comm2d, ierr )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( lwc, lwc_total, 1 , &
       MPI_REAL, MPI_SUM, comm2d, ierr )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( m1, m1_total, 1 , &
       MPI_REAL, MPI_SUM, comm2d, ierr )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( m2, m2_total, 1 , &
       MPI_REAL, MPI_SUM, comm2d, ierr )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( m3, m3_total, 1 , &
       MPI_REAL, MPI_SUM, comm2d, ierr )
#endif 

!
!--    Calculate number concentration and mean volume averaged radius.
       nr_total = MERGE( nr_total / nrclgb_total,                              &
                         0.0_wp, nrclgb_total > 0.0_wp                         &
                       )
       rm_total = MERGE( ( rm_total /                                          &
                            ( nr_total * factor_volume_to_mass )               &
                          )**0.3333333_wp, 0.0_wp, nrclgb_total > 0.0_wp       &
                       )
!
!--    Check which function should be used to approximate the DSD.
       IF ( isf == 1 )  THEN
          lwc_total = MERGE( lwc_total / nrclgb_total,                         &
                             0.0_wp, nrclgb_total > 0.0_wp                     &
                           )
          m1_total  = MERGE( m1_total / nrclgb_total,                          &
                             0.0_wp, nrclgb_total > 0.0_wp                     &
                           )
          m2_total  = MERGE( m2_total / nrclgb_total,                          &
                             0.0_wp, nrclgb_total > 0.0_wp                     &
                           )
          m3_total  = MERGE( m3_total / nrclgb_total,                          &
                             0.0_wp, nrclgb_total > 0.0_wp                     &
                           )
          zeta = m1_total * m3_total / m2_total**2
          mu   = MAX( ( ( 1.0_wp - zeta ) * 2.0_wp + 1.0_wp ) /                &
                        ( zeta - 1.0_wp ), 0.0_wp                              &
                    )

          lambda = ( pirho_l * nr_total / lwc_total *                          &
                     ( mu + 3.0_wp ) * ( mu + 2.0_wp ) * ( mu + 1.0_wp )       &
                   )**0.3333333_wp
          nr0 = nr_total / gamma( mu + 1.0_wp ) * lambda**( mu + 1.0_wp ) 

          DO  n = 0, n_max-1
             diameter  = r_bin_mid(n) * 2.0_wp
             an_spl(n) = nr0 * diameter**mu * EXP( -lambda * diameter ) *      & 
                         ( r_bin(n+1) - r_bin(n) ) * 2.0_wp 
          ENDDO
       ELSEIF ( isf == 2 )  THEN
          DO  n = 0, n_max-1
             an_spl(n) = nr_total / ( SQRT( 2.0_wp * pi ) *                    &
                                     LOG(sigma_log) * r_bin_mid(n)             &
                                     ) *                                       &
                         EXP( -( LOG( r_bin_mid(n) / rm_total )**2 ) /         &
                               ( 2.0_wp * LOG(sigma_log)**2 )                  & 
                             ) *                                               & 
                         ( r_bin(n+1) - r_bin(n) )
          ENDDO
       ELSEIF( isf == 3 )  THEN
          DO  n = 0, n_max-1 
             an_spl(n) = 3.0_wp * nr_total * r_bin_mid(n)**2 / rm_total**3  *  &
                         EXP( - ( r_bin_mid(n)**3 / rm_total**3 ) )         *  &
                         ( r_bin(n+1) - r_bin(n) )
          ENDDO
       ENDIF
!
!--    Criterion to avoid super droplets with a weighting factor < 1.0.
       an_spl = MAX(an_spl, 1.0_wp)

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                number_of_particles = prt_count(k,j,i)
                IF ( number_of_particles <= 0  .OR.                            &
                     ql(k,j,i) < ql_crit )  CYCLE
                particles => grid_particles(k,j,i)%particles(1:number_of_particles)
                new_particles_gb = 0
!
!--             Start splitting operations. Each particle is checked if it
!--             fulfilled the splitting criterion's. In splitting mode 'cl_av'
!--             a critical radius (radius_split) and a splitting function must
!--             be prescribed (see particles_par). The critical weighting factor
!--             is calculated while approximating a 'gamma', 'log' or 'exp'-
!--             drop size distribution. In this mode the DSD is calculated as
!--             an average over all cloudy grid boxes. Super droplets which
!--             have a larger radius and larger weighting factor are split into
!--             'splitting_factor' super droplets. In this case the splitting
!--             factor is calculated of weighting factor of the super droplet
!--             and the approximated number concentration for droplet of such
!--             a size. Due to the splitting, the weighting factor of the
!--             super droplet and all created clones is reduced by the factor
!--             of 'splitting_facor'.
                DO  n = 1, number_of_particles
                   DO  np = 0, n_max-1
                      IF ( r_bin(np) >= radius_split  .AND.                    &
                           particles(n)%particle_mask  .AND.                   &
                           particles(n)%radius >= r_bin(np)  .AND.             &
                           particles(n)%radius < r_bin(np+1)  .AND.            &
                           particles(n)%weight_factor >= an_spl(np)  )         &
                      THEN
!
!--                      Calculate splitting factor
                         splitting_factor =                                    & 
                             MIN( INT( particles(n)%weight_factor /            &
                                        an_spl(np)                             &
                                     ), splitting_factor_max                   &
                                )
                         IF ( splitting_factor < 2 )  CYCLE
!
!--                      Calculate the new number of particles.
                         new_size = prt_count(k,j,i) + splitting_factor - 1
!
!--                      Cycle if maximum number of particles per grid box
!--                      is greater than the allowed maximum number.
                         IF ( new_size >= max_number_particles_per_gridbox )   & 
                         CYCLE
!
!--                      Reallocate particle array if necessary. 
                         IF ( new_size > SIZE(particles) )  THEN 
                            CALL realloc_particles_array( i, j, k, new_size )
                         ENDIF
                         old_size  = prt_count(k,j,i)
                         new_particles_gb = new_particles_gb +                 &
                                            splitting_factor - 1
!
!--                      Calculate new weighting factor.
                         particles(n)%weight_factor =                          & 
                            particles(n)%weight_factor / splitting_factor
                         tmp_particle = particles(n)
!
!--                      Create splitting_factor-1 new particles.
                         DO  jpp = 1, splitting_factor-1
                            grid_particles(k,j,i)%particles(jpp+old_size) =    &
                                                                    tmp_particle
                         ENDDO
!
!--                      Save the new number of super droplets. 
                         prt_count(k,j,i) = prt_count(k,j,i) +                 &
                                            splitting_factor - 1
                      ENDIF
                   ENDDO
                ENDDO 

             ENDDO
          ENDDO
       ENDDO

    ELSEIF ( i_splitting_mode == 3 )  THEN 

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt

!
!--             Initialize summing variables.
                lwc = 0.0_wp
                m1  = 0.0_wp
                m2  = 0.0_wp
                m3  = 0.0_wp
                nr  = 0.0_wp
                rm  = 0.0_wp  

                new_particles_gb = 0
                number_of_particles = prt_count(k,j,i)
                IF ( number_of_particles <= 0  .OR.                            & 
                     ql(k,j,i) < ql_crit )  CYCLE
                particles => grid_particles(k,j,i)%particles
!
!--             Calculate moments of DSD.
                DO  n = 1, number_of_particles
                   IF ( particles(n)%particle_mask  .AND.                      &
                        particles(n)%radius >= r_min )                         &
                   THEN
                      nr  = nr + particles(n)%weight_factor
                      rm  = rm + factor_volume_to_mass  *                      &
                                 particles(n)%radius**3  *                     &
                                 particles(n)%weight_factor
                      IF ( isf == 1 )  THEN
                         diameter   = particles(n)%radius * 2.0_wp
                         lwc = lwc + factor_volume_to_mass *                   &
                                     particles(n)%radius**3 *                  &
                                     particles(n)%weight_factor 
                         m1  = m1 + particles(n)%weight_factor * diameter
                         m2  = m2 + particles(n)%weight_factor * diameter**2
                         m3  = m3 + particles(n)%weight_factor * diameter**3
                      ENDIF
                   ENDIF
                ENDDO

                IF ( nr <= 0.0_wp  .OR.  rm <= 0.0_wp )  CYCLE
!
!--             Calculate mean volume averaged radius.
                rm = ( rm / ( nr * factor_volume_to_mass ) )**0.3333333_wp
!
!--             Check which function should be used to approximate the DSD.
                IF ( isf == 1 )  THEN
!
!--                Gamma size distribution to calculate  
!--                critical weight_factor (e.g. Marshall + Palmer, 1948).
                   zeta = m1 * m3 / m2**2
                   mu   = MAX( ( ( 1.0_wp - zeta ) * 2.0_wp + 1.0_wp ) /       &
                                ( zeta - 1.0_wp ), 0.0_wp                      &
                             )   
                   lambda = ( pirho_l * nr / lwc *                             &
                              ( mu + 3.0_wp ) * ( mu + 2.0_wp ) *              &
                              ( mu + 1.0_wp )                                  &
                            )**0.3333333_wp
                   nr0 =  ( nr / (gamma( mu + 1.0_wp ) ) ) *                   &
                          lambda**( mu + 1.0_wp ) 

                   DO  n = 0, n_max-1
                      diameter         = r_bin_mid(n) * 2.0_wp
                      an_spl(n) = nr0 * diameter**mu *                         &
                                  EXP( -lambda * diameter ) *                  & 
                                  ( r_bin(n+1) - r_bin(n) ) * 2.0_wp 
                   ENDDO
                ELSEIF ( isf == 2 )  THEN
!
!--                Lognormal size distribution to calculate critical 
!--                weight_factor (e.g. Levin, 1971, Bradley + Stow, 1974).
                   DO  n = 0, n_max-1
                      an_spl(n) = nr / ( SQRT( 2.0_wp * pi ) *                 &
                                              LOG(sigma_log) * r_bin_mid(n)    &
                                        ) *                                    &
                                  EXP( -( LOG( r_bin_mid(n) / rm )**2 ) /      &
                                        ( 2.0_wp * LOG(sigma_log)**2 )         &
                                      ) *                                      &
                                  ( r_bin(n+1) - r_bin(n) )
                   ENDDO
                ELSEIF ( isf == 3 )  THEN
!
!--                Exponential size distribution to calculate critical 
!--                weight_factor (e.g. Berry + Reinhardt, 1974).  
                   DO  n = 0, n_max-1
                      an_spl(n) = 3.0_wp * nr * r_bin_mid(n)**2 / rm**3 *     &
                                  EXP( - ( r_bin_mid(n)**3 / rm**3 ) ) *      &
                                  ( r_bin(n+1) - r_bin(n) )
                   ENDDO
                ENDIF

!
!--             Criterion to avoid super droplets with a weighting factor < 1.0.
                an_spl = MAX(an_spl, 1.0_wp)
!
!--             Start splitting operations. Each particle is checked if it
!--             fulfilled the splitting criterion's. In splitting mode 'gb_av'
!--             a critical radius (radius_split) and a splitting function must 
!--             be prescribed (see particles_par). The critical weighting factor
!--             is calculated while appoximating a 'gamma', 'log' or 'exp'-
!--             drop size distribution. In this mode a DSD is calculated for 
!--             every cloudy grid box. Super droplets which have a larger 
!--             radius and larger weighting factor are split into
!--             'splitting_factor' super droplets. In this case the splitting  
!--             factor is calculated of weighting factor of the super droplet  
!--             and theapproximated number concentration for droplet of such 
!--             a size. Due to the splitting, the weighting factor of the  
!--             super droplet and all created clones is reduced by the factor  
!--             of 'splitting_facor'.
                DO  n = 1, number_of_particles
                   DO  np = 0, n_max-1
                      IF ( r_bin(np) >= radius_split  .AND.                    &
                           particles(n)%particle_mask  .AND.                   &
                           particles(n)%radius >= r_bin(np)    .AND.           &
                           particles(n)%radius < r_bin(np+1)   .AND.           &
                           particles(n)%weight_factor >= an_spl(np) )          &
                      THEN
!
!--                      Calculate splitting factor.
                         splitting_factor =                                    & 
                             MIN( INT( particles(n)%weight_factor /            &
                                        an_spl(np)                             &
                                     ), splitting_factor_max                   &
                                )
                         IF ( splitting_factor < 2 )  CYCLE

!
!--                      Calculate the new number of particles.
                         new_size = prt_count(k,j,i) + splitting_factor - 1
!
!--                      Cycle if maximum number of particles per grid box
!--                      is greater than the allowed maximum number.
                         IF ( new_size >= max_number_particles_per_gridbox )   &
                         CYCLE
!
!--                      Reallocate particle array if necessary.
                         IF ( new_size > SIZE(particles) )  THEN 
                            CALL realloc_particles_array( i, j, k, new_size )
                         ENDIF
!
!--                      Calculate new weighting factor.
                         particles(n)%weight_factor = & 
                            particles(n)%weight_factor / splitting_factor
                         tmp_particle               = particles(n)
                         old_size                   = prt_count(k,j,i)
!
!--                      Create splitting_factor-1 new particles.
                         DO  jpp = 1, splitting_factor-1
                            grid_particles(k,j,i)%particles( jpp + old_size ) = &
                               tmp_particle
                         ENDDO
!
!--                      Save the new number of droplets for every grid box.
                         prt_count(k,j,i)    = prt_count(k,j,i) +              &
                                               splitting_factor - 1
                         new_particles_gb    = new_particles_gb +              &
                                               splitting_factor - 1
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    CALL cpu_log( log_point_s(80), 'lpm_splitting', 'stop' )

 END SUBROUTINE lpm_splitting
 

!------------------------------------------------------------------------------!
! Description:
! ------------
! This routine is a part of the Lagrangian particle model. Two Super droplets 
! which fulfill certain criterion's (e.g. a big weighting factor and a small
! radius) can be merged into one super droplet with a increased number of 
! represented particles of the super droplet. This mechanism ensures an
! improved a feasible amount of computational costs. The limits of particle 
! creation should be chosen carefully! The idea of this algorithm is based on 
! Unterstrasser and Soelch, 2014. 
!------------------------------------------------------------------------------!
 SUBROUTINE lpm_merging

    INTEGER(iwp) ::  i         !<
    INTEGER(iwp) ::  j         !<
    INTEGER(iwp) ::  k         !<
    INTEGER(iwp) ::  n         !<
    INTEGER(iwp) ::  merge_drp = 0     !< number of merged droplets


    REAL(wp) ::  ql_crit = 1.0E-5_wp  !< threshold lwc for cloudy grid cells 
                                      !< (e.g. Siebesma et al 2003, JAS, 60)

    CALL cpu_log( log_point_s(81), 'lpm_merging', 'start' )

    merge_drp  = 0

    IF ( weight_factor_merge == -1.0_wp )  THEN
       weight_factor_merge = 0.5_wp * initial_weighting_factor 
    ENDIF

    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt

             number_of_particles = prt_count(k,j,i)
             IF ( number_of_particles <= 0  .OR.                               &
                   ql(k,j,i) >= ql_crit )  CYCLE
             particles => grid_particles(k,j,i)%particles(1:number_of_particles)
!
!--          Start merging operations: This routine delete super droplets with
!--          a small radius (radius <= radius_merge) and a low weighting 
!--          factor (weight_factor  <= weight_factor_merge). The number of 
!--          represented particles will be added to the next particle of the 
!--          particle array. Tests showed that this simplified method can be 
!--          used because it will only take place outside of cloudy grid 
!--          boxes where ql <= 1.0E-5 kg/kg. Therefore, especially former cloned
!--          and subsequent evaporated super droplets will be merged.
             DO  n = 1, number_of_particles-1
                IF ( particles(n)%particle_mask                    .AND.       &
                     particles(n+1)%particle_mask                  .AND.       &
                     particles(n)%radius        <= radius_merge    .AND.       &
                     particles(n)%weight_factor <= weight_factor_merge )       &
                THEN
                   particles(n+1)%weight_factor  =                             &
                                       particles(n+1)%weight_factor +          &
                                       ( particles(n)%radius**3     /          &
                                         particles(n+1)%radius**3   *          &
                                         particles(n)%weight_factor            &
                                       )
                   particles(n)%particle_mask = .FALSE.
                   deleted_particles          = deleted_particles + 1 
                   merge_drp                  = merge_drp + 1

                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO


    CALL cpu_log( log_point_s(81), 'lpm_merging', 'stop' )

 END SUBROUTINE lpm_merging

 

 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Exchange between subdomains.
!> As soon as one particle has moved beyond the boundary of the domain, it
!> is included in the relevant transfer arrays and marked for subsequent
!> deletion on this PE.
!> First sweep for crossings in x direction. Find out first the number of
!> particles to be transferred and allocate temporary arrays needed to store
!> them.
!> For a one-dimensional decomposition along y, no transfer is necessary,
!> because the particle remains on the PE, but the particle coordinate has to
!> be adjusted.
!------------------------------------------------------------------------------!
 SUBROUTINE lpm_exchange_horiz

    INTEGER(iwp) ::  ip                !< index variable along x
    INTEGER(iwp) ::  jp                !< index variable along y
    INTEGER(iwp) ::  kp                !< index variable along z
    INTEGER(iwp) ::  n                 !< particle index variable 

#if defined( __parallel )
    INTEGER(iwp) ::  i                 !< grid index (x) of particle positition
    INTEGER(iwp) ::  j                 !< grid index (y) of particle positition
    INTEGER(iwp) ::  par_size          !< Particle size in bytes
    INTEGER(iwp) ::  trlp_count        !< number of particles send to left PE
    INTEGER(iwp) ::  trlp_count_recv   !< number of particles receive from right PE
    INTEGER(iwp) ::  trnp_count        !< number of particles send to north PE
    INTEGER(iwp) ::  trnp_count_recv   !< number of particles receive from south PE
    INTEGER(iwp) ::  trrp_count        !< number of particles send to right PE
    INTEGER(iwp) ::  trrp_count_recv   !< number of particles receive from left PE
    INTEGER(iwp) ::  trsp_count        !< number of particles send to south PE
    INTEGER(iwp) ::  trsp_count_recv   !< number of particles receive from north PE

    TYPE(particle_type), DIMENSION(:), ALLOCATABLE ::  rvlp  !< particles received from right PE
    TYPE(particle_type), DIMENSION(:), ALLOCATABLE ::  rvnp  !< particles received from south PE
    TYPE(particle_type), DIMENSION(:), ALLOCATABLE ::  rvrp  !< particles received from left PE
    TYPE(particle_type), DIMENSION(:), ALLOCATABLE ::  rvsp  !< particles received from north PE
    TYPE(particle_type), DIMENSION(:), ALLOCATABLE ::  trlp  !< particles send to left PE
    TYPE(particle_type), DIMENSION(:), ALLOCATABLE ::  trnp  !< particles send to north PE
    TYPE(particle_type), DIMENSION(:), ALLOCATABLE ::  trrp  !< particles send to right PE
    TYPE(particle_type), DIMENSION(:), ALLOCATABLE ::  trsp  !< particles send to south PE
#endif

    CALL cpu_log( log_point_s(23), 'lpm_exchange_horiz', 'start' )

#if defined( __parallel )

!
!-- Exchange between subdomains.
!-- As soon as one particle has moved beyond the boundary of the domain, it
!-- is included in the relevant transfer arrays and marked for subsequent
!-- deletion on this PE.
!-- First sweep for crossings in x direction. Find out first the number of
!-- particles to be transferred and allocate temporary arrays needed to store
!-- them.
!-- For a one-dimensional decomposition along y, no transfer is necessary,
!-- because the particle remains on the PE, but the particle coordinate has to
!-- be adjusted.
    trlp_count  = 0
    trrp_count  = 0

    trlp_count_recv   = 0
    trrp_count_recv   = 0

    IF ( pdims(1) /= 1 )  THEN
!
!--    First calculate the storage necessary for sending and receiving the data.
!--    Compute only first (nxl) and last (nxr) loop iterration.
       DO  ip = nxl, nxr, nxr - nxl
          DO  jp = nys, nyn
             DO  kp = nzb+1, nzt

                number_of_particles = prt_count(kp,jp,ip)
                IF ( number_of_particles <= 0 )  CYCLE
                particles => grid_particles(kp,jp,ip)%particles(1:number_of_particles)
                DO  n = 1, number_of_particles
                   IF ( particles(n)%particle_mask )  THEN
                      i = particles(n)%x * ddx
!
!--                   Above calculation does not work for indices less than zero
                      IF ( particles(n)%x < 0.0_wp)  i = -1

                      IF ( i < nxl )  THEN
                         trlp_count = trlp_count + 1
                      ELSEIF ( i > nxr )  THEN
                         trrp_count = trrp_count + 1
                      ENDIF
                   ENDIF
                ENDDO

             ENDDO
          ENDDO
       ENDDO

       IF ( trlp_count  == 0 )  trlp_count  = 1
       IF ( trrp_count  == 0 )  trrp_count  = 1

       ALLOCATE( trlp(trlp_count), trrp(trrp_count) )

       trlp = zero_particle
       trrp = zero_particle

       trlp_count  = 0
       trrp_count  = 0

    ENDIF
!
!-- Compute only first (nxl) and last (nxr) loop iterration
    DO  ip = nxl, nxr, nxr-nxl
       DO  jp = nys, nyn
          DO  kp = nzb+1, nzt
             number_of_particles = prt_count(kp,jp,ip)
             IF ( number_of_particles <= 0 )  CYCLE
             particles => grid_particles(kp,jp,ip)%particles(1:number_of_particles)
             DO  n = 1, number_of_particles
!
!--             Only those particles that have not been marked as 'deleted' may
!--             be moved.
                IF ( particles(n)%particle_mask )  THEN

                   i = particles(n)%x * ddx
!
!--                Above calculation does not work for indices less than zero
                   IF ( particles(n)%x < 0.0_wp )  i = -1

                   IF ( i <  nxl )  THEN
                      IF ( i < 0 )  THEN
!
!--                      Apply boundary condition along x
                         IF ( ibc_par_lr == 0 )  THEN
!
!--                         Cyclic condition
                            IF ( pdims(1) == 1 )  THEN
                               particles(n)%x        = ( nx + 1 ) * dx + particles(n)%x
                               particles(n)%origin_x = ( nx + 1 ) * dx + &
                               particles(n)%origin_x
                            ELSE
                               trlp_count = trlp_count + 1
                               trlp(trlp_count)   = particles(n)
                               trlp(trlp_count)%x = ( nx + 1 ) * dx + trlp(trlp_count)%x
                               trlp(trlp_count)%origin_x = trlp(trlp_count)%origin_x + &
                               ( nx + 1 ) * dx
                               particles(n)%particle_mask  = .FALSE.
                               deleted_particles = deleted_particles + 1

                               IF ( trlp(trlp_count)%x >= (nx + 1)* dx - 1.0E-12_wp )  THEN
                                  trlp(trlp_count)%x = trlp(trlp_count)%x - 1.0E-10_wp
                                  !++ why is 1 subtracted in next statement???
                                  trlp(trlp_count)%origin_x = trlp(trlp_count)%origin_x - 1
                               ENDIF

                            ENDIF

                         ELSEIF ( ibc_par_lr == 1 )  THEN
!
!--                         Particle absorption
                            particles(n)%particle_mask = .FALSE.
                            deleted_particles = deleted_particles + 1

                         ELSEIF ( ibc_par_lr == 2 )  THEN
!
!--                         Particle reflection
                            particles(n)%x       = -particles(n)%x
                            particles(n)%speed_x = -particles(n)%speed_x

                         ENDIF
                      ELSE
!
!--                      Store particle data in the transfer array, which will be 
!--                      send to the neighbouring PE
                         trlp_count = trlp_count + 1
                         trlp(trlp_count) = particles(n)
                         particles(n)%particle_mask = .FALSE.
                         deleted_particles = deleted_particles + 1

                      ENDIF

                   ELSEIF ( i > nxr )  THEN
                      IF ( i > nx )  THEN
!
!--                      Apply boundary condition along x
                         IF ( ibc_par_lr == 0 )  THEN
!
!--                         Cyclic condition
                            IF ( pdims(1) == 1 )  THEN
                               particles(n)%x = particles(n)%x - ( nx + 1 ) * dx
                               particles(n)%origin_x = particles(n)%origin_x - &
                               ( nx + 1 ) * dx
                            ELSE
                               trrp_count = trrp_count + 1
                               trrp(trrp_count) = particles(n)
                               trrp(trrp_count)%x = trrp(trrp_count)%x - ( nx + 1 ) * dx
                               trrp(trrp_count)%origin_x = trrp(trrp_count)%origin_x - &
                               ( nx + 1 ) * dx
                               particles(n)%particle_mask = .FALSE.
                               deleted_particles = deleted_particles + 1

                            ENDIF

                         ELSEIF ( ibc_par_lr == 1 )  THEN
!
!--                         Particle absorption
                            particles(n)%particle_mask = .FALSE.
                            deleted_particles = deleted_particles + 1

                         ELSEIF ( ibc_par_lr == 2 )  THEN
!
!--                         Particle reflection
                            particles(n)%x       = 2 * ( nx * dx ) - particles(n)%x
                            particles(n)%speed_x = -particles(n)%speed_x

                         ENDIF
                      ELSE
!
!--                      Store particle data in the transfer array, which will be send
!--                      to the neighbouring PE
                         trrp_count = trrp_count + 1
                         trrp(trrp_count) = particles(n)
                         particles(n)%particle_mask = .FALSE.
                         deleted_particles = deleted_particles + 1

                      ENDIF

                   ENDIF
                ENDIF

             ENDDO
          ENDDO
       ENDDO
    ENDDO

!
!-- STORAGE_SIZE returns the storage size of argument A in bits. However , it 
!-- is needed in bytes. The function C_SIZEOF which produces this value directly
!-- causes problems with gfortran. For this reason the use of C_SIZEOF is avoided
    par_size = STORAGE_SIZE(trlp(1))/8


!
!-- Allocate arrays required for north-south exchange, as these
!-- are used directly after particles are exchange along x-direction.
    ALLOCATE( move_also_north(1:NR_2_direction_move) )
    ALLOCATE( move_also_south(1:NR_2_direction_move) )

    nr_move_north = 0
    nr_move_south = 0
!
!-- Send left boundary, receive right boundary (but first exchange how many
!-- and check, if particle storage must be extended)
    IF ( pdims(1) /= 1 )  THEN

       CALL MPI_SENDRECV( trlp_count,      1, MPI_INTEGER, pleft,  0, &
                          trrp_count_recv, 1, MPI_INTEGER, pright, 0, &
                          comm2d, status, ierr )

       ALLOCATE(rvrp(MAX(1,trrp_count_recv)))

       CALL MPI_SENDRECV( trlp, max(1,trlp_count)*par_size, MPI_BYTE,&
                          pleft, 1, rvrp,                            &
                          max(1,trrp_count_recv)*par_size, MPI_BYTE, pright, 1,&
                          comm2d, status, ierr )

       IF ( trrp_count_recv > 0 )  CALL lpm_add_particles_to_gridcell(rvrp(1:trrp_count_recv))

       DEALLOCATE(rvrp)

!
!--    Send right boundary, receive left boundary
       CALL MPI_SENDRECV( trrp_count,      1, MPI_INTEGER, pright, 0, &
                          trlp_count_recv, 1, MPI_INTEGER, pleft,  0, &
                          comm2d, status, ierr )

       ALLOCATE(rvlp(MAX(1,trlp_count_recv)))
!
!--    This MPI_SENDRECV should work even with odd mixture on 32 and 64 Bit 
!--    variables in structure particle_type (due to the calculation of par_size)
       CALL MPI_SENDRECV( trrp, max(1,trrp_count)*par_size, MPI_BYTE,&
                          pright, 1, rvlp,                           &
                          max(1,trlp_count_recv)*par_size, MPI_BYTE, pleft, 1, &
                          comm2d, status, ierr )

       IF ( trlp_count_recv > 0 )  CALL lpm_add_particles_to_gridcell(rvlp(1:trlp_count_recv))

       DEALLOCATE( rvlp )
       DEALLOCATE( trlp, trrp )

    ENDIF

!
!-- Check whether particles have crossed the boundaries in y direction. Note 
!-- that this case can also apply to particles that have just been received
!-- from the adjacent right or left PE.
!-- Find out first the number of particles to be transferred and allocate
!-- temporary arrays needed to store them.
!-- For a one-dimensional decomposition along y, no transfer is necessary,
!-- because the particle remains on the PE.
    trsp_count  = nr_move_south
    trnp_count  = nr_move_north

    trsp_count_recv   = 0
    trnp_count_recv   = 0

    IF ( pdims(2) /= 1 )  THEN
!
!--    First calculate the storage necessary for sending and receiving the
!--    data
       DO  ip = nxl, nxr
          DO  jp = nys, nyn, nyn-nys    !compute only first (nys) and last (nyn) loop iterration
             DO  kp = nzb+1, nzt
                number_of_particles = prt_count(kp,jp,ip)
                IF ( number_of_particles <= 0 )  CYCLE
                particles => grid_particles(kp,jp,ip)%particles(1:number_of_particles)
                DO  n = 1, number_of_particles
                   IF ( particles(n)%particle_mask )  THEN
                      j = particles(n)%y * ddy
!
!--                   Above calculation does not work for indices less than zero
                      IF ( particles(n)%y < 0.0_wp)  j = -1

                      IF ( j < nys )  THEN
                         trsp_count = trsp_count + 1
                      ELSEIF ( j > nyn )  THEN
                         trnp_count = trnp_count + 1
                      ENDIF
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDDO

       IF ( trsp_count  == 0 )  trsp_count  = 1
       IF ( trnp_count  == 0 )  trnp_count  = 1

       ALLOCATE( trsp(trsp_count), trnp(trnp_count) )

       trsp = zero_particle
       trnp = zero_particle

       trsp_count  = nr_move_south
       trnp_count  = nr_move_north

       trsp(1:nr_move_south) = move_also_south(1:nr_move_south)
       trnp(1:nr_move_north) = move_also_north(1:nr_move_north)

    ENDIF

    DO  ip = nxl, nxr
       DO  jp = nys, nyn, nyn-nys ! compute only first (nys) and last (nyn) loop iterration
          DO  kp = nzb+1, nzt
             number_of_particles = prt_count(kp,jp,ip)
             IF ( number_of_particles <= 0 )  CYCLE
             particles => grid_particles(kp,jp,ip)%particles(1:number_of_particles)
             DO  n = 1, number_of_particles
!
!--             Only those particles that have not been marked as 'deleted' may
!--             be moved.
                IF ( particles(n)%particle_mask )  THEN

                   j = particles(n)%y * ddy
!
!--                Above calculation does not work for indices less than zero
                   IF ( particles(n)%y < 0.0_wp )  j = -1

                   IF ( j < nys )  THEN
                      IF ( j < 0 )  THEN
!
!--                      Apply boundary condition along y
                         IF ( ibc_par_ns == 0 )  THEN
!
!--                         Cyclic condition
                            IF ( pdims(2) == 1 )  THEN
                               particles(n)%y = ( ny + 1 ) * dy + particles(n)%y
                               particles(n)%origin_y = ( ny + 1 ) * dy + &
                                                     particles(n)%origin_y
                            ELSE
                               trsp_count         = trsp_count + 1
                               trsp(trsp_count)   = particles(n)
                               trsp(trsp_count)%y = ( ny + 1 ) * dy + &
                                                 trsp(trsp_count)%y
                               trsp(trsp_count)%origin_y = trsp(trsp_count)%origin_y &
                                                + ( ny + 1 ) * dy
                               particles(n)%particle_mask = .FALSE.
                               deleted_particles = deleted_particles + 1

                               IF ( trsp(trsp_count)%y >= (ny+1)* dy - 1.0E-12_wp )  THEN
                                  trsp(trsp_count)%y = trsp(trsp_count)%y - 1.0E-10_wp
                                  !++ why is 1 subtracted in next statement???
                                  trsp(trsp_count)%origin_y =                        &
                                                  trsp(trsp_count)%origin_y - 1
                               ENDIF

                            ENDIF

                         ELSEIF ( ibc_par_ns == 1 )  THEN
!
!--                         Particle absorption
                            particles(n)%particle_mask = .FALSE.
                            deleted_particles          = deleted_particles + 1

                         ELSEIF ( ibc_par_ns == 2 )  THEN
!
!--                         Particle reflection
                            particles(n)%y       = -particles(n)%y
                            particles(n)%speed_y = -particles(n)%speed_y

                         ENDIF
                      ELSE
!
!--                      Store particle data in the transfer array, which will 
!--                      be send to the neighbouring PE
                         trsp_count = trsp_count + 1
                         trsp(trsp_count) = particles(n)
                         particles(n)%particle_mask = .FALSE.
                         deleted_particles = deleted_particles + 1

                      ENDIF

                   ELSEIF ( j > nyn )  THEN
                      IF ( j > ny )  THEN
!
!--                       Apply boundary condition along y
                         IF ( ibc_par_ns == 0 )  THEN
!
!--                         Cyclic condition
                            IF ( pdims(2) == 1 )  THEN
                               particles(n)%y        = particles(n)%y - ( ny + 1 ) * dy
                               particles(n)%origin_y =                         &
                                          particles(n)%origin_y - ( ny + 1 ) * dy
                            ELSE
                               trnp_count         = trnp_count + 1
                               trnp(trnp_count)   = particles(n)
                               trnp(trnp_count)%y =                            &
                                          trnp(trnp_count)%y - ( ny + 1 ) * dy
                               trnp(trnp_count)%origin_y =                     &
                                         trnp(trnp_count)%origin_y - ( ny + 1 ) * dy
                               particles(n)%particle_mask = .FALSE.
                               deleted_particles          = deleted_particles + 1
                            ENDIF

                         ELSEIF ( ibc_par_ns == 1 )  THEN
!
!--                         Particle absorption
                            particles(n)%particle_mask = .FALSE.
                            deleted_particles = deleted_particles + 1

                         ELSEIF ( ibc_par_ns == 2 )  THEN
!
!--                         Particle reflection
                            particles(n)%y       = 2 * ( ny * dy ) - particles(n)%y
                            particles(n)%speed_y = -particles(n)%speed_y

                         ENDIF
                      ELSE
!
!--                      Store particle data in the transfer array, which will 
!--                      be send to the neighbouring PE
                         trnp_count = trnp_count + 1
                         trnp(trnp_count) = particles(n)
                         particles(n)%particle_mask = .FALSE.
                         deleted_particles = deleted_particles + 1

                      ENDIF

                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO

!
!-- Send front boundary, receive back boundary (but first exchange how many
!-- and check, if particle storage must be extended)
    IF ( pdims(2) /= 1 )  THEN

       CALL MPI_SENDRECV( trsp_count,      1, MPI_INTEGER, psouth, 0, &
                          trnp_count_recv, 1, MPI_INTEGER, pnorth, 0, &
                          comm2d, status, ierr )

       ALLOCATE(rvnp(MAX(1,trnp_count_recv)))
!
!--    This MPI_SENDRECV should work even with odd mixture on 32 and 64 Bit 
!--    variables in structure particle_type (due to the calculation of par_size)
       CALL MPI_SENDRECV( trsp, trsp_count*par_size, MPI_BYTE,      &
                          psouth, 1, rvnp,                             &
                          trnp_count_recv*par_size, MPI_BYTE, pnorth, 1,   &
                          comm2d, status, ierr )

       IF ( trnp_count_recv  > 0 )  CALL lpm_add_particles_to_gridcell(rvnp(1:trnp_count_recv))

       DEALLOCATE(rvnp)

!
!--    Send back boundary, receive front boundary
       CALL MPI_SENDRECV( trnp_count,      1, MPI_INTEGER, pnorth, 0, &
                          trsp_count_recv, 1, MPI_INTEGER, psouth, 0, &
                          comm2d, status, ierr )

       ALLOCATE(rvsp(MAX(1,trsp_count_recv)))
!
!--    This MPI_SENDRECV should work even with odd mixture on 32 and 64 Bit 
!--    variables in structure particle_type (due to the calculation of par_size)
       CALL MPI_SENDRECV( trnp, trnp_count*par_size, MPI_BYTE,      &
                          pnorth, 1, rvsp,                          &
                          trsp_count_recv*par_size, MPI_BYTE, psouth, 1,   &
                          comm2d, status, ierr )

       IF ( trsp_count_recv > 0 )  CALL lpm_add_particles_to_gridcell(rvsp(1:trsp_count_recv))


       DEALLOCATE(rvsp)

       number_of_particles = number_of_particles + trsp_count_recv

       DEALLOCATE( trsp, trnp )

    ENDIF

    DEALLOCATE( move_also_north )
    DEALLOCATE( move_also_south )

#else

    DO  ip = nxl, nxr, nxr-nxl
       DO  jp = nys, nyn
          DO  kp = nzb+1, nzt
             number_of_particles = prt_count(kp,jp,ip)
             IF ( number_of_particles <= 0 )  CYCLE
             particles => grid_particles(kp,jp,ip)%particles(1:number_of_particles)
             DO  n = 1, number_of_particles
!
!--             Apply boundary conditions

                IF ( particles(n)%x < 0.0_wp )  THEN

                   IF ( ibc_par_lr == 0 )  THEN
!
!--                   Cyclic boundary. Relevant coordinate has to be changed.
                      particles(n)%x = ( nx + 1 ) * dx + particles(n)%x
                      particles(n)%origin_x = ( nx + 1 ) * dx + &
                               particles(n)%origin_x
                   ELSEIF ( ibc_par_lr == 1 )  THEN
!
!--                   Particle absorption
                      particles(n)%particle_mask = .FALSE.
                      deleted_particles = deleted_particles + 1

                   ELSEIF ( ibc_par_lr == 2 )  THEN
!
!--                   Particle reflection
                      particles(n)%x       = -dx - particles(n)%x
                      particles(n)%speed_x = -particles(n)%speed_x
                   ENDIF

                ELSEIF ( particles(n)%x >= ( nx + 1) * dx )  THEN

                   IF ( ibc_par_lr == 0 )  THEN
!
!--                   Cyclic boundary. Relevant coordinate has to be changed.
                      particles(n)%x = particles(n)%x - ( nx + 1 ) * dx
                      particles(n)%origin_x = particles(n)%origin_x - &
                               ( nx + 1 ) * dx

                   ELSEIF ( ibc_par_lr == 1 )  THEN
!
!--                   Particle absorption
                      particles(n)%particle_mask = .FALSE.
                      deleted_particles = deleted_particles + 1

                   ELSEIF ( ibc_par_lr == 2 )  THEN
!
!--                   Particle reflection
                      particles(n)%x       = ( nx + 1 ) * dx - particles(n)%x
                      particles(n)%speed_x = -particles(n)%speed_x
                   ENDIF

                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    DO  ip = nxl, nxr
       DO  jp = nys, nyn, nyn-nys
          DO  kp = nzb+1, nzt
             number_of_particles = prt_count(kp,jp,ip)
             IF ( number_of_particles <= 0 )  CYCLE
             particles => grid_particles(kp,jp,ip)%particles(1:number_of_particles)
             DO  n = 1, number_of_particles

                IF ( particles(n)%y < 0.0_wp)  THEN

                   IF ( ibc_par_ns == 0 )  THEN
!
!--                   Cyclic boundary. Relevant coordinate has to be changed.
                      particles(n)%y = ( ny + 1 ) * dy + particles(n)%y
                      particles(n)%origin_y = ( ny + 1 ) * dy + &
                           particles(n)%origin_y

                   ELSEIF ( ibc_par_ns == 1 )  THEN
!
!--                   Particle absorption
                      particles(n)%particle_mask = .FALSE.
                      deleted_particles = deleted_particles + 1

                   ELSEIF ( ibc_par_ns == 2 )  THEN
!
!--                   Particle reflection
                      particles(n)%y       = -dy - particles(n)%y
                      particles(n)%speed_y = -particles(n)%speed_y
                   ENDIF

                ELSEIF ( particles(n)%y >= ( ny + 1) * dy )  THEN

                   IF ( ibc_par_ns == 0 )  THEN
!
!--                   Cyclic boundary. Relevant coordinate has to be changed.
                      particles(n)%y = particles(n)%y - ( ny + 1 ) * dy
                      particles(n)%origin_y = particles(n)%origin_y - &
                                ( ny + 1 ) * dy

                   ELSEIF ( ibc_par_ns == 1 )  THEN
!
!--                   Particle absorption
                      particles(n)%particle_mask = .FALSE.
                      deleted_particles = deleted_particles + 1

                   ELSEIF ( ibc_par_ns == 2 )  THEN
!
!--                   Particle reflection
                      particles(n)%y       = ( ny + 1 ) * dy - particles(n)%y
                      particles(n)%speed_y = -particles(n)%speed_y
                   ENDIF

                ENDIF

             ENDDO
          ENDDO
       ENDDO
    ENDDO
#endif

!
!-- Accumulate the number of particles transferred between the subdomains
#if defined( __parallel )
    trlp_count_sum      = trlp_count_sum      + trlp_count
    trlp_count_recv_sum = trlp_count_recv_sum + trlp_count_recv
    trrp_count_sum      = trrp_count_sum      + trrp_count
    trrp_count_recv_sum = trrp_count_recv_sum + trrp_count_recv
    trsp_count_sum      = trsp_count_sum      + trsp_count
    trsp_count_recv_sum = trsp_count_recv_sum + trsp_count_recv
    trnp_count_sum      = trnp_count_sum      + trnp_count
    trnp_count_recv_sum = trnp_count_recv_sum + trnp_count_recv
#endif

    CALL cpu_log( log_point_s(23), 'lpm_exchange_horiz', 'stop' )

 END SUBROUTINE lpm_exchange_horiz

#if defined( __parallel )
!------------------------------------------------------------------------------!
! Description:
! ------------
!> If a particle moves from one processor to another, this subroutine moves 
!> the corresponding elements from the particle arrays of the old grid cells 
!> to the particle arrays of the new grid cells.
!------------------------------------------------------------------------------!
 SUBROUTINE lpm_add_particles_to_gridcell (particle_array)

    IMPLICIT NONE

    INTEGER(iwp)        ::  ip        !< grid index (x) of particle
    INTEGER(iwp)        ::  jp        !< grid index (x) of particle
    INTEGER(iwp)        ::  kp        !< grid index (x) of particle
    INTEGER(iwp)        ::  n         !< index variable of particle
    INTEGER(iwp)        ::  pindex    !< dummy argument for new number of particles per grid box

    LOGICAL             ::  pack_done !<

    TYPE(particle_type), DIMENSION(:), INTENT(IN)  ::  particle_array !< new particles in a grid box
    TYPE(particle_type), DIMENSION(:), ALLOCATABLE ::  temp_ns        !< temporary particle array for reallocation

    pack_done     = .FALSE.

    DO  n = 1, SIZE(particle_array)

       IF ( .NOT. particle_array(n)%particle_mask )  CYCLE

       ip = particle_array(n)%x * ddx
       jp = particle_array(n)%y * ddy
!
!--    In case of stretching the actual k index must be found
       IF ( dz_stretch_level /= -9999999.9_wp  .OR.         &
            dz_stretch_level_start(1) /= -9999999.9_wp )  THEN
          kp = MINLOC( ABS( particle_array(n)%z - zu ), DIM = 1 ) - 1
       ELSE
          kp = INT( particle_array(n)%z / dz(1) + 1 + offset_ocean_nzt )
       ENDIF

       IF ( ip >= nxl  .AND.  ip <= nxr  .AND.  jp >= nys  .AND.  jp <= nyn    &
            .AND.  kp >= nzb+1  .AND.  kp <= nzt)  THEN ! particle stays on processor
          number_of_particles = prt_count(kp,jp,ip)
          particles => grid_particles(kp,jp,ip)%particles(1:number_of_particles)

          pindex = prt_count(kp,jp,ip)+1
          IF( pindex > SIZE(grid_particles(kp,jp,ip)%particles) )  THEN
             IF ( pack_done )  THEN
                CALL realloc_particles_array ( ip, jp, kp )
             ELSE
                CALL lpm_pack
                prt_count(kp,jp,ip) = number_of_particles
                pindex = prt_count(kp,jp,ip)+1
                IF ( pindex > SIZE(grid_particles(kp,jp,ip)%particles) )  THEN
                   CALL realloc_particles_array ( ip, jp, kp )
                ENDIF
                pack_done = .TRUE.
             ENDIF
          ENDIF
          grid_particles(kp,jp,ip)%particles(pindex) = particle_array(n)
          prt_count(kp,jp,ip) = pindex
       ELSE
          IF ( jp <= nys - 1 )  THEN
             nr_move_south = nr_move_south+1
!
!--          Before particle information is swapped to exchange-array, check 
!--          if enough memory is allocated. If required, reallocate exchange
!--          array.
             IF ( nr_move_south > SIZE(move_also_south) )  THEN
!
!--             At first, allocate further temporary array to swap particle 
!--             information.
                ALLOCATE( temp_ns(SIZE(move_also_south)+NR_2_direction_move) )
                temp_ns(1:nr_move_south-1) = move_also_south(1:nr_move_south-1)
                DEALLOCATE( move_also_south )
                ALLOCATE( move_also_south(SIZE(temp_ns)) )
                move_also_south(1:nr_move_south-1) = temp_ns(1:nr_move_south-1)
                DEALLOCATE( temp_ns )

             ENDIF

             move_also_south(nr_move_south) = particle_array(n)

             IF ( jp == -1 )  THEN
!
!--             Apply boundary condition along y
                IF ( ibc_par_ns == 0 )  THEN
                   move_also_south(nr_move_south)%y =                          &
                      move_also_south(nr_move_south)%y + ( ny + 1 ) * dy
                   move_also_south(nr_move_south)%origin_y =                   &
                      move_also_south(nr_move_south)%origin_y + ( ny + 1 ) * dy
                ELSEIF ( ibc_par_ns == 1 )  THEN
!
!--                Particle absorption
                   move_also_south(nr_move_south)%particle_mask = .FALSE.
                   deleted_particles = deleted_particles + 1

                ELSEIF ( ibc_par_ns == 2 )  THEN
!
!--                Particle reflection
                   move_also_south(nr_move_south)%y       =                    &
                      -move_also_south(nr_move_south)%y
                   move_also_south(nr_move_south)%speed_y =                    &
                      -move_also_south(nr_move_south)%speed_y

                ENDIF
             ENDIF
          ELSEIF ( jp >= nyn+1 )  THEN
             nr_move_north = nr_move_north+1
!
!--          Before particle information is swapped to exchange-array, check 
!--          if enough memory is allocated. If required, reallocate exchange
!--          array.
             IF ( nr_move_north > SIZE(move_also_north) )  THEN
!
!--             At first, allocate further temporary array to swap particle 
!--             information.
                ALLOCATE( temp_ns(SIZE(move_also_north)+NR_2_direction_move) )
                temp_ns(1:nr_move_north-1) = move_also_south(1:nr_move_north-1)
                DEALLOCATE( move_also_north )
                ALLOCATE( move_also_north(SIZE(temp_ns)) )
                move_also_north(1:nr_move_north-1) = temp_ns(1:nr_move_north-1)
                DEALLOCATE( temp_ns )

             ENDIF

             move_also_north(nr_move_north) = particle_array(n)
             IF ( jp == ny+1 )  THEN
!
!--             Apply boundary condition along y
                IF ( ibc_par_ns == 0 )  THEN

                   move_also_north(nr_move_north)%y =                          &
                      move_also_north(nr_move_north)%y - ( ny + 1 ) * dy
                   move_also_north(nr_move_north)%origin_y =                   &
                      move_also_north(nr_move_north)%origin_y - ( ny + 1 ) * dy
                ELSEIF ( ibc_par_ns == 1 )  THEN
!
!--                Particle absorption
                   move_also_north(nr_move_north)%particle_mask = .FALSE.
                   deleted_particles = deleted_particles + 1

                ELSEIF ( ibc_par_ns == 2 )  THEN
!
!--                Particle reflection
                   move_also_north(nr_move_north)%y       =                    &
                      -move_also_north(nr_move_north)%y
                   move_also_north(nr_move_north)%speed_y =                    &
                      -move_also_north(nr_move_north)%speed_y

                ENDIF
             ENDIF
          ELSE
             IF ( .NOT. child_domain )  THEN
                WRITE(0,'(a,8i7)') 'particle out of range ',myid,ip,jp,kp,nxl,nxr,nys,nyn
             ENDIF
          ENDIF
       ENDIF
    ENDDO

 END SUBROUTINE lpm_add_particles_to_gridcell
#endif
 
 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> If a particle moves from one grid cell to another (on the current 
!> processor!), this subroutine moves the corresponding element from the
!> particle array of the old grid cell to the particle array of the new grid 
!> cell.
!------------------------------------------------------------------------------!
 SUBROUTINE lpm_move_particle
 
    INTEGER(iwp)        ::  i           !< grid index (x) of particle position
    INTEGER(iwp)        ::  ip          !< index variable along x
    INTEGER(iwp)        ::  j           !< grid index (y) of particle position
    INTEGER(iwp)        ::  jp          !< index variable along y
    INTEGER(iwp)        ::  k           !< grid index (z) of particle position
    INTEGER(iwp)        ::  kp          !< index variable along z
    INTEGER(iwp)        ::  n           !< index variable for particle array
    INTEGER(iwp)        ::  np_before_move !< number of particles per grid box before moving
    INTEGER(iwp)        ::  pindex      !< dummy argument for number of new particle per grid box

    TYPE(particle_type), DIMENSION(:), POINTER  ::  particles_before_move !< particles before moving

    CALL cpu_log( log_point_s(41), 'lpm_move_particle', 'start' )
    CALL lpm_check_cfl
    DO  ip = nxl, nxr
       DO  jp = nys, nyn
          DO  kp = nzb+1, nzt

             np_before_move = prt_count(kp,jp,ip)
             IF ( np_before_move <= 0 )  CYCLE
             particles_before_move => grid_particles(kp,jp,ip)%particles(1:np_before_move)

             DO  n = 1, np_before_move
                i = particles_before_move(n)%x * ddx
                j = particles_before_move(n)%y * ddy
                k = kp
!
!--             Find correct vertical particle grid box (necessary in case of grid stretching)
!--             Due to the CFL limitations only the neighbouring grid boxes are considered. 
                IF( zw(k)   < particles_before_move(n)%z ) k = k + 1
                IF( zw(k-1) > particles_before_move(n)%z ) k = k - 1 

!--             For lpm_exchange_horiz to work properly particles need to be moved to the outermost gridboxes
!--             of the respective processor. If the particle index is inside the processor the following lines
!--             will not change the index
                i = MIN ( i , nxr )
                i = MAX ( i , nxl )
                j = MIN ( j , nyn )
                j = MAX ( j , nys )

                k = MIN ( k , nzt )
                k = MAX ( k , nzb+1 )

!
!--             Check, if particle has moved to another grid cell.
                IF ( i /= ip  .OR.  j /= jp  .OR.  k /= kp )  THEN
!!
!--                If the particle stays on the same processor, the particle
!--                will be added to the particle array of the new processor.
                   number_of_particles = prt_count(k,j,i)
                   particles => grid_particles(k,j,i)%particles(1:number_of_particles)

                   pindex = prt_count(k,j,i)+1
                   IF (  pindex > SIZE(grid_particles(k,j,i)%particles)  )     &
                   THEN
                      CALL realloc_particles_array( i, j, k )
                   ENDIF

                   grid_particles(k,j,i)%particles(pindex) = particles_before_move(n)
                   prt_count(k,j,i) = pindex

                   particles_before_move(n)%particle_mask = .FALSE.
                ENDIF
             ENDDO

          ENDDO
       ENDDO
    ENDDO

    CALL cpu_log( log_point_s(41), 'lpm_move_particle', 'stop' )

    RETURN

 END SUBROUTINE lpm_move_particle
 

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check CFL-criterion for each particle. If one particle violated the 
!> criterion the particle will be deleted and a warning message is given.
!------------------------------------------------------------------------------!
 SUBROUTINE lpm_check_cfl  

    IMPLICIT NONE

    INTEGER(iwp)  ::  i !< running index, x-direction
    INTEGER(iwp)  ::  j !< running index, y-direction
    INTEGER(iwp)  ::  k !< running index, z-direction
    INTEGER(iwp)  ::  n !< running index, number of particles

    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt
             number_of_particles = prt_count(k,j,i)
             IF ( number_of_particles <= 0 )  CYCLE
             particles => grid_particles(k,j,i)%particles(1:number_of_particles)         
             DO  n = 1, number_of_particles
!
!--             Note, check for CFL does not work at first particle timestep 
!--             when both, age and age_m are zero. 
                IF ( particles(n)%age - particles(n)%age_m > 0.0_wp )  THEN
                   IF( ABS( particles(n)%speed_x ) >                           &
                      ( dx / ( particles(n)%age - particles(n)%age_m) )  .OR.  &
                       ABS( particles(n)%speed_y ) >                           & 
                      ( dy / ( particles(n)%age - particles(n)%age_m) )  .OR.  &
                       ABS( particles(n)%speed_z ) >                           &
                      ( ( zw(k)-zw(k-1) )                                      &
                      / ( particles(n)%age - particles(n)%age_m) ) )  THEN
                      WRITE( message_string, * )                               &
                      'Particle violated CFL-criterion: &particle with id ',   &
                      particles(n)%id, ' will be deleted!'   
                      CALL message( 'lpm_check_cfl', 'PA0475', 0, 1, -1, 6, 0 )

                      particles(n)%particle_mask= .FALSE.
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO   

 END SUBROUTINE lpm_check_cfl
 
 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> If the allocated memory for the particle array do not suffice to add arriving
!> particles from neighbour grid cells, this subrouting reallocates the 
!> particle array to assure enough memory is available. 
!------------------------------------------------------------------------------!
 SUBROUTINE realloc_particles_array ( i, j, k, size_in )

    INTEGER(iwp), INTENT(IN)                       ::  i              !<
    INTEGER(iwp), INTENT(IN)                       ::  j              !<
    INTEGER(iwp), INTENT(IN)                       ::  k              !<
    INTEGER(iwp), INTENT(IN), OPTIONAL             ::  size_in        !<

    INTEGER(iwp)                                   ::  old_size        !<
    INTEGER(iwp)                                   ::  new_size        !<
    TYPE(particle_type), DIMENSION(:), ALLOCATABLE ::  tmp_particles_d !<
    TYPE(particle_type), DIMENSION(500)            ::  tmp_particles_s !<

    old_size = SIZE(grid_particles(k,j,i)%particles)

    IF ( PRESENT(size_in) )   THEN
       new_size = size_in
    ELSE
       new_size = old_size * ( 1.0_wp + alloc_factor / 100.0_wp )
    ENDIF

    new_size = MAX( new_size, 1, old_size + 1 )

    IF ( old_size <= 500 )  THEN

       tmp_particles_s(1:old_size) = grid_particles(k,j,i)%particles(1:old_size)

       DEALLOCATE(grid_particles(k,j,i)%particles)
       ALLOCATE(grid_particles(k,j,i)%particles(new_size))

       grid_particles(k,j,i)%particles(1:old_size)          = tmp_particles_s(1:old_size)
       grid_particles(k,j,i)%particles(old_size+1:new_size) = zero_particle

    ELSE

       ALLOCATE(tmp_particles_d(new_size))
       tmp_particles_d(1:old_size) = grid_particles(k,j,i)%particles

       DEALLOCATE(grid_particles(k,j,i)%particles)
       ALLOCATE(grid_particles(k,j,i)%particles(new_size))

       grid_particles(k,j,i)%particles(1:old_size)          = tmp_particles_d(1:old_size)
       grid_particles(k,j,i)%particles(old_size+1:new_size) = zero_particle

       DEALLOCATE(tmp_particles_d)

    ENDIF
    particles => grid_particles(k,j,i)%particles(1:new_size)

    RETURN
    
 END SUBROUTINE realloc_particles_array
 
 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Not needed but allocated space for particles is dealloced. 
!------------------------------------------------------------------------------!
 SUBROUTINE dealloc_particles_array

 
    INTEGER(iwp) ::  i               !<
    INTEGER(iwp) ::  j               !<
    INTEGER(iwp) ::  k               !<
    INTEGER(iwp) ::  old_size        !<
    INTEGER(iwp) ::  new_size        !<

    LOGICAL ::  dealloc

    TYPE(particle_type), DIMENSION(:), ALLOCATABLE ::  tmp_particles_d !<
    TYPE(particle_type), DIMENSION(500)            ::  tmp_particles_s !<

    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt
!
!--          Determine number of active particles
             number_of_particles = prt_count(k,j,i)
!
!--          Determine allocated memory size
             old_size = SIZE( grid_particles(k,j,i)%particles )
!
!--          Check for large unused memory
             dealloc = ( ( number_of_particles < 1 .AND.         &
                           old_size            > 1 )  .OR.       &
                         ( number_of_particles > 1 .AND.         &
                           old_size - number_of_particles *                    &
                              ( 1.0_wp + 0.01_wp * alloc_factor ) > 0.0_wp ) )

             IF ( dealloc )  THEN
                IF ( number_of_particles < 1 )  THEN
                   new_size = 1
                ELSE 
                   new_size = INT( number_of_particles * ( 1.0_wp + 0.01_wp * alloc_factor ) )
                ENDIF

                IF ( number_of_particles <= 500 )  THEN

                   tmp_particles_s(1:number_of_particles) = grid_particles(k,j,i)%particles(1:number_of_particles)

                   DEALLOCATE(grid_particles(k,j,i)%particles)
                   ALLOCATE(grid_particles(k,j,i)%particles(new_size))

                   grid_particles(k,j,i)%particles(1:number_of_particles)          = tmp_particles_s(1:number_of_particles)
                   grid_particles(k,j,i)%particles(number_of_particles+1:new_size) = zero_particle

                ELSE

                   ALLOCATE(tmp_particles_d(number_of_particles))
                   tmp_particles_d(1:number_of_particles) = grid_particles(k,j,i)%particles(1:number_of_particles)

                   DEALLOCATE(grid_particles(k,j,i)%particles)
                   ALLOCATE(grid_particles(k,j,i)%particles(new_size))

                   grid_particles(k,j,i)%particles(1:number_of_particles)          = tmp_particles_d(1:number_of_particles)
                   grid_particles(k,j,i)%particles(number_of_particles+1:new_size) = zero_particle

                   DEALLOCATE(tmp_particles_d)

                ENDIF

             ENDIF
          ENDDO
       ENDDO
    ENDDO

 END SUBROUTINE dealloc_particles_array 
 
 
!------------------------------------------------------------------------------!
! Description:
! -----------
!> Routine for the whole processor
!> Sort all particles into the 8 respective subgrid boxes (in case of trilinear
!> interpolation method) and free space of particles which has been marked for
!> deletion.
!------------------------------------------------------------------------------!
   SUBROUTINE lpm_sort_and_delete

       INTEGER(iwp) ::  i  !<
       INTEGER(iwp) ::  ip !<
       INTEGER(iwp) ::  is !<
       INTEGER(iwp) ::  j  !<
       INTEGER(iwp) ::  jp !<
       INTEGER(iwp) ::  kp !<
       INTEGER(iwp) ::  m  !<
       INTEGER(iwp) ::  n  !<
       INTEGER(iwp) ::  nn !<
       INTEGER(iwp) ::  sort_index  !<

       INTEGER(iwp), DIMENSION(0:7) ::  sort_count  !<

       TYPE(particle_type), DIMENSION(:,:), ALLOCATABLE ::  sort_particles    !<

       CALL cpu_log( log_point_s(51), 'lpm_sort_and_delete', 'start' )
       IF ( interpolation_trilinear )  THEN
          DO  ip = nxl, nxr
             DO  jp = nys, nyn
                DO  kp = nzb+1, nzt
                   number_of_particles = prt_count(kp,jp,ip)
                   IF ( number_of_particles <= 0 )  CYCLE
                   particles => grid_particles(kp,jp,ip)%particles(1:number_of_particles)
                   nn = 0
                   sort_count = 0
                   ALLOCATE( sort_particles(number_of_particles, 0:7) )

                   DO  n = 1, number_of_particles
                      sort_index = 0

                      IF ( particles(n)%particle_mask )  THEN
                         nn = nn + 1
!
!--                      Sorting particles with a binary scheme
!--                      sort_index=111_2=7_10 -> particle at the left,south,bottom subgridbox
!--                      sort_index=000_2=0_10 -> particle at the right,north,top subgridbox
!--                      For this the center of the gridbox is calculated
                         i = (particles(n)%x + 0.5_wp * dx) * ddx
                         j = (particles(n)%y + 0.5_wp * dy) * ddy

                         IF ( i == ip )  sort_index = sort_index + 4
                         IF ( j == jp )  sort_index = sort_index + 2
                         IF ( zu(kp) > particles(n)%z ) sort_index = sort_index + 1

                         sort_count(sort_index) = sort_count(sort_index) + 1
                         m = sort_count(sort_index)
                         sort_particles(m,sort_index) = particles(n)
                         sort_particles(m,sort_index)%block_nr = sort_index
                      ENDIF
                   ENDDO
!
!--                Delete and resort particles by overwritting and set
!--                the number_of_particles to the actual value.
                   nn = 0
                   DO  is = 0,7
                      grid_particles(kp,jp,ip)%start_index(is) = nn + 1
                      DO  n = 1,sort_count(is)
                         nn = nn + 1
                         particles(nn) = sort_particles(n,is)
                      ENDDO
                      grid_particles(kp,jp,ip)%end_index(is) = nn
                   ENDDO

                   number_of_particles = nn
                   prt_count(kp,jp,ip) = number_of_particles
                   DEALLOCATE(sort_particles)
                ENDDO
             ENDDO
          ENDDO

!--    In case of the simple interpolation method the particles must not
!--    be sorted in subboxes. Particles marked for deletion however, must be
!--    deleted and number of particles must be recalculated as it is also
!--    done for the trilinear particle advection interpolation method.
       ELSE

          DO  ip = nxl, nxr
             DO  jp = nys, nyn
                DO  kp = nzb+1, nzt

                   number_of_particles = prt_count(kp,jp,ip)
                   IF ( number_of_particles <= 0 )  CYCLE
                   particles => grid_particles(kp,jp,ip)%particles(1:number_of_particles)
!
!--                Repack particles array, i.e. delete particles and recalculate
!--                number of particles
                   CALL lpm_pack
                   prt_count(kp,jp,ip) = number_of_particles
                ENDDO
             ENDDO
          ENDDO
       ENDIF
       CALL cpu_log( log_point_s(51), 'lpm_sort_and_delete', 'stop' )

    END SUBROUTINE lpm_sort_and_delete

 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Move all particles not marked for deletion to lowest indices (packing)
!------------------------------------------------------------------------------!
    SUBROUTINE lpm_pack

       INTEGER(iwp) ::  n       !<
       INTEGER(iwp) ::  nn      !<
!
!--    Find out elements marked for deletion and move data from highest index
!--    values to these free indices
       nn = number_of_particles

       DO WHILE ( .NOT. particles(nn)%particle_mask )
          nn = nn-1
          IF ( nn == 0 )  EXIT
       ENDDO

       IF ( nn > 0 )  THEN
          DO  n = 1, number_of_particles
             IF ( .NOT. particles(n)%particle_mask )  THEN
                particles(n) = particles(nn)
                nn = nn - 1
                DO WHILE ( .NOT. particles(nn)%particle_mask )
                   nn = nn-1
                   IF ( n == nn )  EXIT
                ENDDO
             ENDIF
             IF ( n == nn )  EXIT
          ENDDO
       ENDIF

!
!--    The number of deleted particles has been determined in routines
!--    lpm_boundary_conds, lpm_droplet_collision, and lpm_exchange_horiz
       number_of_particles = nn

    END SUBROUTINE lpm_pack 


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Sort particles in each sub-grid box into two groups: particles that already
!> completed the LES timestep, and particles that need further timestepping to
!> complete the LES timestep. 
!------------------------------------------------------------------------------!
    SUBROUTINE lpm_sort_timeloop_done

       INTEGER(iwp) ::  end_index     !< particle end index for each sub-box
       INTEGER(iwp) ::  i             !< index of particle grid box in x-direction
       INTEGER(iwp) ::  j             !< index of particle grid box in y-direction
       INTEGER(iwp) ::  k             !< index of particle grid box in z-direction
       INTEGER(iwp) ::  n             !< running index for number of particles
       INTEGER(iwp) ::  nb            !< index of subgrid boux
       INTEGER(iwp) ::  nf            !< indices for particles in each sub-box that already finalized their substeps
       INTEGER(iwp) ::  nnf           !< indices for particles in each sub-box that need further treatment
       INTEGER(iwp) ::  num_finalized !< number of particles in each sub-box that already finalized their substeps
       INTEGER(iwp) ::  start_index   !< particle start index for each sub-box

       TYPE(particle_type), DIMENSION(:), ALLOCATABLE ::  sort_particles  !< temporary particle array

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt

                number_of_particles = prt_count(k,j,i)
                IF ( number_of_particles <= 0 )  CYCLE
                particles => grid_particles(k,j,i)%particles(1:number_of_particles)

                DO  nb = 0, 7
!
!--                Obtain start and end index for each subgrid box
                   start_index = grid_particles(k,j,i)%start_index(nb)
                   end_index   = grid_particles(k,j,i)%end_index(nb)
!
!--                Allocate temporary array used for sorting. 
                   ALLOCATE( sort_particles(start_index:end_index) )
!
!--                Determine number of particles already completed the LES 
!--                timestep, and write them into a temporary array. 
                   nf = start_index
                   num_finalized = 0
                   DO  n = start_index, end_index
                      IF ( dt_3d - particles(n)%dt_sum < 1E-8_wp )  THEN
                         sort_particles(nf) = particles(n)
                         nf                 = nf + 1
                         num_finalized      = num_finalized + 1
                      ENDIF
                   ENDDO
!
!--                Determine number of particles that not completed the LES 
!--                timestep, and write them into a temporary array. 
                   nnf = nf
                   DO  n = start_index, end_index
                      IF ( dt_3d - particles(n)%dt_sum > 1E-8_wp )  THEN
                         sort_particles(nnf) = particles(n)
                         nnf                 = nnf + 1
                      ENDIF
                   ENDDO
!
!--                Write back sorted particles
                   particles(start_index:end_index) =                          &
                                           sort_particles(start_index:end_index)
!
!--                Determine updated start_index, used to masked already 
!--                completed particles. 
                   grid_particles(k,j,i)%start_index(nb) =                     &
                                      grid_particles(k,j,i)%start_index(nb)    &
                                    + num_finalized
!
!--                Deallocate dummy array
                   DEALLOCATE ( sort_particles )
!
!--                Finally, if number of non-completed particles is non zero 
!--                in any of the sub-boxes, set control flag appropriately. 
                   IF ( nnf > nf )                                             &
                      grid_particles(k,j,i)%time_loop_done = .FALSE.

                ENDDO
             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE lpm_sort_timeloop_done 

END MODULE lagrangian_particle_model_mod
