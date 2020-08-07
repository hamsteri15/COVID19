!> @file parin.f90
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
! $Id: parin.f90 4360 2020-01-07 11:25:50Z suehring $
! removed recycling_yshift
! 
! 4227 2019-09-10 18:04:34Z gronemeier
! implement new palm_date_time_mod
! 
! 4146 2019-08-07 07:47:36Z gronemeier
! added rotation_angle to initialization_parameters
! 
! 4191 2019-08-27 15:45:07Z gronemeier
! bugfix: add recycling_method_for_thermodynamic_quantities to inipar namelist
! 
! 4183 2019-08-23 07:33:16Z oliver.maas
! replaced recycle_absolute_quantities by recycling_method_for_thermodynamic_quantities
! 
! 4182 2019-08-22 15:20:23Z scharf
! Corrected "Former revisions" section
! 
! 4176 2019-08-20 14:10:41Z oliver.maas
! added recycle_absolute_quantities to initialization_parameters namelist
! 
! 4173 2019-08-20 12:04:06Z gronemeier
! add vdi_internal_controls
! 
! 4131 2019-08-02 11:06:18Z monakurppa
! Allocate hom and hom_sum to allow profile output for salsa variables.
! 
! 4079 2019-07-09 18:04:41Z suehring
! +monotonic_limiter_z
! 
! 4022 2019-06-12 11:52:39Z suehring
! Change default top boundary condition for pressure to Neumann in offline
! nesting case
! 
! 4017 2019-06-06 12:16:46Z schwenkel
! Introduce alternative switch for debug output during timestepping
! 
! 3885 2019-04-11 11:29:34Z kanani
! Changes related to global restructuring of location messages and introduction 
! of additional debug messages
! 
! 3806 2019-03-21 12:45:50Z raasch
! additional check for lateral boundary conditions added
! 
! 3747 2019-02-16 15:15:23Z gronemeier
! removed setting of parameter region
! 
! 3746 2019-02-16 12:41:27Z gronemeier
! Removed most_method
! 
! 3649 2019-01-02 16:52:21Z suehring
! Delete debug-print statements
!
! Revision 1.1  1997/07/24 11:22:50  raasch
! Initial revision
!
!
! Description:
! ------------
!> This subroutine reads variables controling the run from the NAMELIST files
!>
!> @todo: Revise max_pr_cs (profiles for chemistry)
!------------------------------------------------------------------------------!
 SUBROUTINE parin
 

    USE arrays_3d,                                                             &
        ONLY:  pt_init, q_init, ref_state, s_init, sa_init,                    &
               ug, u_init, v_init, vg

    USE chem_modules

    USE control_parameters

    USE cpulog,                                                                &
        ONLY:  cpu_log_barrierwait

    USE grid_variables,                                                        &
        ONLY:  dx, dy

    USE indices,                                                               &
        ONLY:  nx, ny, nz

    USE kinds

    USE model_1d_mod,                                                          &
        ONLY:  damp_level_1d, dt_pr_1d, dt_run_control_1d, end_time_1d

    USE module_interface,                                                      &
        ONLY:  module_interface_parin

    USE netcdf_interface,                                                      &
        ONLY:  netcdf_data_format, netcdf_deflate, netcdf_precision

    USE pegrid

    USE pmc_interface,                                                         &
        ONLY:  nested_run, nesting_mode

    USE profil_parameter,                                                      &
        ONLY:  cross_profiles, profile_columns, profile_rows

    USE progress_bar,                                                          &
        ONLY :  progress_bar_disabled

    USE read_restart_data_mod,                                                 &
        ONLY:  rrd_global

    USE statistics,                                                            &
        ONLY:  hom, hom_sum, pr_palm, statistic_regions

    USE turbulence_closure_mod,                                                &
        ONLY:  rans_const_c, rans_const_sigma

    USE vertical_nesting_mod,                                                  &
        ONLY:  vnest_start_time


    IMPLICIT NONE

    CHARACTER (LEN=80) ::  line  !< dummy string that contains the current line of the parameter file 

    INTEGER(iwp) ::  global_id      !< process id with respect to MPI_COMM_WORLD
    INTEGER(iwp) ::  global_procs   !< # of procs with respect to MPI_COMM_WORLD
    INTEGER(iwp) ::  i              !<
    INTEGER(iwp) ::  ioerr          !< error flag for open/read/write

    NAMELIST /inipar/  alpha_surface, approximation, bc_e_b,     &
                       bc_lr, bc_ns, bc_p_b, bc_p_t, bc_pt_b, bc_pt_t, bc_q_b, &
             bc_q_t,bc_s_b, bc_s_t, bc_uv_b, bc_uv_t,                 &
             building_height, building_length_x,          &
             building_length_y, building_wall_left, building_wall_south,       &
             calc_soil_moisture_during_spinup,                                 &
             call_psolver_at_all_substeps,  &
             canyon_height,                                                    &
             canyon_width_x, canyon_width_y, canyon_wall_left,                 &
             canyon_wall_south, cfl_factor, cloud_droplets,   &
             collective_wait, complex_terrain,           &
             conserve_volume_flow,                                             &
             conserve_volume_flow_mode, constant_flux_layer,                   &
             coupling_start_time,             &
             cycle_mg, damp_level_1d,                                          &
             data_output_during_spinup,                                        &
             origin_date_time,                                                 &
             dissipation_1d,                                                   &
             dp_external, dp_level_b, dp_smooth, dpdxy,    &
             dt, dt_pr_1d, dt_run_control_1d, dt_spinup, dx, dy, dz, dz_max,   &
             dz_stretch_factor, dz_stretch_level, dz_stretch_level_start,      &
             dz_stretch_level_end, end_time_1d, ensemble_member_nr, e_init,    &
             e_min, fft_method, flux_input_mode, flux_output_mode,             &
             galilei_transformation, humidity,                                 &
             inflow_damping_height, inflow_damping_width,                      &
             inflow_disturbance_begin, inflow_disturbance_end,                 &
             initializing_actions, km_constant,                                &
             large_scale_forcing, large_scale_subsidence, latitude,            &
             longitude,                                 &
             loop_optimization, lsf_exception, masking_method, mg_cycles,      &
             mg_switch_to_pe0_level, mixing_length_1d, momentum_advec,         &
             monotonic_limiter_z,                                              &
             netcdf_precision, neutral, ngsrb,                                 &
             nsor, nsor_ini, nudging, nx, ny, nz, ocean_mode, omega,           &
             omega_sor, outflow_source_plane, passive_scalar,                  &
             prandtl_number, psolver, pt_damping_factor,        &
             pt_damping_width, pt_reference, pt_surface,                       &
             pt_surface_initial_change, pt_vertical_gradient,                  &
             pt_vertical_gradient_level, q_surface, q_surface_initial_change,  &
             q_vertical_gradient, q_vertical_gradient_level,                   &
             random_generator, random_heatflux, rans_const_c, rans_const_sigma,&
             rayleigh_damping_factor, rayleigh_damping_height,                 &
             recycling_method_for_thermodynamic_quantities, recycling_width,   &
             reference_state, residual_limit,                                  &
             rotation_angle,                                                   &
             roughness_length,                                                 &
             scalar_advec,   &
             scalar_rayleigh_damping,                              &
             spinup_time, spinup_pt_amplitude, spinup_pt_mean,                 &
             statistic_regions, subs_vertical_gradient,                        &
             subs_vertical_gradient_level, surface_heatflux, surface_pressure, &
             surface_scalarflux, surface_waterflux,                            &
             s_surface, s_surface_initial_change, s_vertical_gradient,         &
             s_vertical_gradient_level, timestep_scheme,                       &
             topography, topography_grid_convention, top_heatflux,             &
             top_momentumflux_u, top_momentumflux_v,                           &
             top_scalarflux, transpose_compute_overlap,                        &
             tunnel_height, tunnel_length, tunnel_width_x, tunnel_width_y,     &
             tunnel_wall_depth, turbulence_closure,                            &
             turbulent_inflow, turbulent_outflow,                              &
             use_subsidence_tendencies, ug_surface, ug_vertical_gradient,      &
             use_free_convection_scaling,                                      &
             ug_vertical_gradient_level, use_surface_fluxes, use_cmax,         &
             use_top_fluxes, use_ug_for_galilei_tr, use_upstream_for_tke,      &
             uv_heights, u_bulk, u_profile, vdi_checks, vg_surface,            &
             vg_vertical_gradient,  &
             vg_vertical_gradient_level, v_bulk, v_profile,&
             wall_adjustment, wall_heatflux, wall_humidityflux,                &
             wall_scalarflux, y_shift, zeta_max, zeta_min,  &
             z0h_factor

    NAMELIST /initialization_parameters/  alpha_surface,                       &
             approximation, bc_e_b,                                            &
             bc_lr, bc_ns, bc_p_b, bc_p_t, bc_pt_b, bc_pt_t, bc_q_b,           &
             bc_q_t,bc_s_b, bc_s_t, bc_uv_b, bc_uv_t,                          &
             building_height, building_length_x,                               &
             building_length_y, building_wall_left, building_wall_south,       &
             calc_soil_moisture_during_spinup,                                 &
             call_psolver_at_all_substeps,                                     &
             canyon_height,                                                    &
             canyon_width_x, canyon_width_y, canyon_wall_left,                 &
             canyon_wall_south, cfl_factor, cloud_droplets,                    &
             collective_wait, complex_terrain,                                 &
             conserve_volume_flow,                                             &
             conserve_volume_flow_mode, constant_flux_layer,                   &
             coupling_start_time,                                              &
             cycle_mg, damp_level_1d,                                          &
             data_output_during_spinup,                                        &
             origin_date_time,                                                 &
             dissipation_1d,                                                   &
             dp_external, dp_level_b, dp_smooth, dpdxy,                        &
             dt, dt_pr_1d, dt_run_control_1d, dt_spinup, dx, dy, dz, dz_max,   &
             dz_stretch_factor, dz_stretch_level, dz_stretch_level_start,      &
             dz_stretch_level_end, end_time_1d, ensemble_member_nr, e_init,    &
             e_min, fft_method, flux_input_mode, flux_output_mode,             &
             galilei_transformation, humidity,                                 &
             inflow_damping_height, inflow_damping_width,                      &
             inflow_disturbance_begin, inflow_disturbance_end,                 &
             initializing_actions, km_constant,                                &
             large_scale_forcing, large_scale_subsidence, latitude,            &
             longitude,                                                        &
             loop_optimization, lsf_exception, masking_method, mg_cycles,      &
             mg_switch_to_pe0_level, mixing_length_1d, momentum_advec,         &
             monotonic_limiter_z,                                              &
             netcdf_precision, neutral, ngsrb,                                 &
             nsor, nsor_ini, nudging, nx, ny, nz, ocean_mode, omega,           &
             omega_sor, outflow_source_plane, passive_scalar,                  &
             prandtl_number, psolver, pt_damping_factor,                       &
             pt_damping_width, pt_reference, pt_surface,                       &
             pt_surface_initial_change, pt_vertical_gradient,                  &
             pt_vertical_gradient_level, q_surface, q_surface_initial_change,  &
             q_vertical_gradient, q_vertical_gradient_level,                   &
             random_generator, random_heatflux, rans_const_c, rans_const_sigma,&
             rayleigh_damping_factor, rayleigh_damping_height,                 &
             recycling_method_for_thermodynamic_quantities, recycling_width,   &
             reference_state, residual_limit,                                  &
             rotation_angle,                                                   &
             roughness_length, scalar_advec,                                   &
             scalar_rayleigh_damping,                                          &
             spinup_time, spinup_pt_amplitude, spinup_pt_mean,                 &
             statistic_regions, subs_vertical_gradient,                        &
             subs_vertical_gradient_level, surface_heatflux, surface_pressure, &
             surface_scalarflux, surface_waterflux,                            &
             s_surface, s_surface_initial_change, s_vertical_gradient,         &
             s_vertical_gradient_level, timestep_scheme,                       &
             topography, topography_grid_convention, top_heatflux,             &
             top_momentumflux_u, top_momentumflux_v,                           &
             top_scalarflux, transpose_compute_overlap,                        &
             tunnel_height, tunnel_length, tunnel_width_x, tunnel_width_y,     &
             tunnel_wall_depth, turbulence_closure,                            &
             turbulent_inflow, turbulent_outflow,                              &
             use_subsidence_tendencies, ug_surface, ug_vertical_gradient,      &
             ug_vertical_gradient_level, use_surface_fluxes, use_cmax,         &
             use_top_fluxes, use_ug_for_galilei_tr, use_upstream_for_tke,      &
             use_free_convection_scaling,                                      &
             uv_heights, u_bulk, u_profile, vdi_checks,                        &
             vg_surface, vg_vertical_gradient,  &
             vg_vertical_gradient_level, v_bulk, v_profile,                    &
             wall_adjustment, wall_heatflux, wall_humidityflux,                &
             wall_scalarflux, y_shift, zeta_max, zeta_min, z0h_factor
             
    NAMELIST /d3par/  averaging_interval, averaging_interval_pr,               &
             cpu_log_barrierwait, create_disturbances,                         &
             cross_profiles, data_output, data_output_masks,                   &
             data_output_pr, data_output_2d_on_each_pe,                        &
             debug_output,                                                     &
             debug_output_timestep,                                            &
             disturbance_amplitude,                                            &
             disturbance_energy_limit, disturbance_level_b,                    &
             disturbance_level_t, do2d_at_begin, do3d_at_begin,                &
             dt, dt_averaging_input, dt_averaging_input_pr,                    &
             dt_coupling, dt_data_output, dt_data_output_av, dt_disturb,       &
             dt_domask, dt_dopr, dt_dopr_listing, dt_dots, dt_do2d_xy,         &
             dt_do2d_xz, dt_do2d_yz, dt_do3d, dt_max, dt_restart,              &
             dt_run_control,end_time, force_print_header, mask_k_over_surface, &
             mask_scale_x,                                                     &
             mask_scale_y, mask_scale_z, mask_x, mask_y, mask_z, mask_x_loop,  &
             mask_y_loop, mask_z_loop, netcdf_data_format, netcdf_deflate,     &
             normalizing_region, npex, npey, nz_do3d,                          &
             profile_columns, profile_rows,     &
             restart_time, section_xy, section_xz, section_yz,                 &
             skip_time_data_output, skip_time_data_output_av, skip_time_dopr,  &
             skip_time_do2d_xy, skip_time_do2d_xz, skip_time_do2d_yz,          &
             skip_time_do3d, skip_time_domask, synchronous_exchange,           &
             termination_time_needed, vnest_start_time

    NAMELIST /runtime_parameters/  averaging_interval, averaging_interval_pr,  &
             cpu_log_barrierwait, create_disturbances,                         &
             cross_profiles, data_output, data_output_masks,                   &
             data_output_pr, data_output_2d_on_each_pe,                        &
             debug_output,                                                     &
             debug_output_timestep,                                            &
             disturbance_amplitude,                                            &
             disturbance_energy_limit, disturbance_level_b,                    &
             disturbance_level_t, do2d_at_begin, do3d_at_begin,                &
             dt, dt_averaging_input, dt_averaging_input_pr,                    &
             dt_coupling, dt_data_output, dt_data_output_av, dt_disturb,       &
             dt_domask, dt_dopr, dt_dopr_listing, dt_dots, dt_do2d_xy,         &
             dt_do2d_xz, dt_do2d_yz, dt_do3d, dt_max, dt_restart,              &
             dt_run_control,end_time, force_print_header, mask_k_over_surface, &
             mask_scale_x,                                                     &
             mask_scale_y, mask_scale_z, mask_x, mask_y, mask_z, mask_x_loop,  &
             mask_y_loop, mask_z_loop, netcdf_data_format, netcdf_deflate,     &
             normalizing_region, npex, npey, nz_do3d,                          &
             profile_columns, profile_rows,     &
             restart_time, section_xy, section_xz, section_yz,                 &
             skip_time_data_output, skip_time_data_output_av, skip_time_dopr,  &
             skip_time_do2d_xy, skip_time_do2d_xz, skip_time_do2d_yz,          &
             skip_time_do3d, skip_time_domask, synchronous_exchange,           &
             termination_time_needed, vnest_start_time

    NAMELIST /envpar/  progress_bar_disabled, host,                            &
                       maximum_cpu_time_allowed, maximum_parallel_io_streams,  &
                       read_svf, revision, run_identifier, tasks_per_node,     &
                       write_binary, write_svf

!
!-- First read values of environment variables (this NAMELIST file is
!-- generated by palmrun)
    CALL location_message( 'reading environment parameters from ENVPAR', 'start' )

    OPEN ( 90, FILE='ENVPAR', STATUS='OLD', FORM='FORMATTED', IOSTAT=ioerr )

    IF ( ioerr /= 0 )  THEN
       message_string = 'local file ENVPAR not found' //                       &
                        '&some variables for steering may not be properly set'
       CALL message( 'parin', 'PA0276', 0, 1, 0, 6, 0 )
    ELSE
       READ ( 90, envpar, IOSTAT=ioerr )
       IF ( ioerr < 0 )  THEN
          message_string = 'no envpar-NAMELIST found in local file '  //       &
                           'ENVPAR& or some variables for steering may '  //   &
                           'not be properly set'
          CALL message( 'parin', 'PA0278', 0, 1, 0, 6, 0 )
       ELSEIF ( ioerr > 0 )  THEN
          message_string = 'errors in local file ENVPAR' //                    &
                           '&some variables for steering may not be properly set'
          CALL message( 'parin', 'PA0277', 0, 1, 0, 6, 0 )
       ENDIF
       CLOSE ( 90 )
    ENDIF

    CALL location_message( 'reading environment parameters from ENVPAR', 'finished' )
!
!-- Calculate the number of groups into which parallel I/O is split.
!-- The default for files which are opened by all PEs (or where each
!-- PE opens his own independent file) is, that all PEs are doing input/output
!-- in parallel at the same time. This might cause performance or even more
!-- severe problems depending on the configuration of the underlying file
!-- system.
!-- Calculation of the number of blocks and the I/O group must be based on all
!-- PEs involved in this run. Since myid and numprocs are related to the
!-- comm2d communicator, which gives only a subset of all PEs in case of
!-- nested runs, that information must be inquired again from the global
!-- communicator.
!-- First, set the default:
#if defined( __parallel )
    CALL MPI_COMM_RANK( MPI_COMM_WORLD, global_id, ierr )
    CALL MPI_COMM_SIZE( MPI_COMM_WORLD, global_procs, ierr )
#else
    global_id    = 0
    global_procs = 1
#endif
    IF ( maximum_parallel_io_streams == -1  .OR.                               &
         maximum_parallel_io_streams > global_procs )  THEN
       maximum_parallel_io_streams = global_procs
    ENDIF
!
!-- Now calculate the number of io_blocks and the io_group to which the
!-- respective PE belongs. I/O of the groups is done in serial, but in parallel
!-- for all PEs belonging to the same group.
    io_blocks = global_procs / maximum_parallel_io_streams
    io_group  = MOD( global_id+1, io_blocks )
    
    CALL location_message( 'reading NAMELIST parameters from PARIN', 'start' )
!
!-- Data is read in parallel by groups of PEs
    DO  i = 0, io_blocks-1
       IF ( i == io_group )  THEN

!
!--       Open the NAMELIST-file which is send with this job
          CALL check_open( 11 )

!
!--       Read the control parameters for initialization.
!--       The namelist "inipar" must be provided in the NAMELIST-file.
          READ ( 11, initialization_parameters, ERR=10, END=11 )
          GOTO 14
          
 10       BACKSPACE( 11 )
          READ( 11 , '(A)') line
          CALL parin_fail_message( 'initialization_parameters', line )

 11       REWIND ( 11 )
          READ ( 11, inipar, ERR=12, END=13 )
 
          message_string = 'namelist inipar is deprecated and will be ' //    &
                          'removed in near future. & Please use namelist ' // &
                          'initialization_parameters instead'
          CALL message( 'parin', 'PA0017', 0, 1, 0, 6, 0 )
 
          GOTO 14
 
 12       BACKSPACE( 11 )
          READ( 11 , '(A)') line
          CALL parin_fail_message( 'inipar', line )

 13       message_string = 'no initialization_parameters-namelist found'
          CALL message( 'parin', 'PA0272', 1, 2, 0, 6, 0 )

!
!--       Try to read runtime parameters given by the user for this run
!--       (namelist "runtime_parameters"). The namelist "runtime_parmeters"    
!--       can be omitted. In that case default values are used for the         
!--       parameters.
 14       line = ' '

          REWIND ( 11 )
          line = ' '
          DO WHILE ( INDEX( line, '&runtime_parameters' ) == 0 )
             READ ( 11, '(A)', END=16 )  line
          ENDDO
          BACKSPACE ( 11 )

!
!--       Read namelist
          READ ( 11, runtime_parameters, ERR = 15 )
          GOTO 18

 15       BACKSPACE( 11 )
          READ( 11 , '(A)') line
          CALL parin_fail_message( 'runtime_parameters', line )

 16       REWIND ( 11 )
          line = ' '
          DO WHILE ( INDEX( line, '&d3par' ) == 0 )
             READ ( 11, '(A)', END=18 )  line
          ENDDO
          BACKSPACE ( 11 )

!
!--       Read namelist
          READ ( 11, d3par, ERR = 17, END = 18 )

          message_string = 'namelist d3par is deprecated and will be ' //      &
                          'removed in near future. &Please use namelist ' //   &
                          'runtime_parameters instead'
          CALL message( 'parin', 'PA0487', 0, 1, 0, 6, 0 )

          GOTO 18

 17       BACKSPACE( 11 )
          READ( 11 , '(A)') line
          CALL parin_fail_message( 'd3par', line )

 18       CONTINUE

!
!--       Check for module namelists and read them
          CALL module_interface_parin

!
!--       If required, read control parameters from restart file (produced by
!--       a prior run). All PEs are reading from file created by PE0 (see
!--       check_open)
          IF ( TRIM( initializing_actions ) == 'read_restart_data' )  THEN

             CALL rrd_global
!
!--          Increment the run count
             runnr = runnr + 1
!
!--          In case of a restart run, the number of user-defined profiles on
!--          the restart file (already stored in max_pr_user) has to match the
!--          one given for the current run. max_pr_user_tmp is calculated in 
!--          user_parin and max_pr_user is read in via rrd_global.
             IF ( max_pr_user /= max_pr_user_tmp )  THEN
                WRITE( message_string, * ) 'the number of user-defined ',      &
                      'profiles given in data_output_pr (', max_pr_user_tmp,   &
                      ') does not match the one ',                             &
                      'found in the restart file (', max_pr_user, ')'
                CALL message( 'user_parin', 'UI0009', 1, 2, 0, 6, 0 )
             ENDIF
          ELSE
             max_pr_user = max_pr_user_tmp
          ENDIF
!
!--       Activate spinup
          IF ( land_surface .OR. urban_surface )  THEN
             IF ( spinup_time > 0.0_wp )  THEN
                coupling_start_time = spinup_time
                time_since_reference_point = simulated_time - coupling_start_time
                IF ( spinup_pt_mean == 9999999.9_wp )  THEN
                   spinup_pt_mean = pt_surface
                ENDIF
                end_time = end_time + spinup_time
                IF ( TRIM( initializing_actions ) /= 'read_restart_data' )  THEN
                   spinup = .TRUE.
                ELSEIF (        TRIM( initializing_actions ) == 'read_restart_data'      &
                         .AND.  time_since_reference_point > 0.0_wp )  THEN
                   data_output_during_spinup = .FALSE.  !< required for correct ntdim calculation
                                                        !< in check_parameters for restart run
                ENDIF
             ENDIF
          ENDIF

!
!--       In case of nested runs, explicitly set nesting boundary conditions. 
!--       This will overwrite the user settings and basic defaults.
!--       bc_lr and bc_ns always need to be cyclic for vertical nesting.
          IF ( nested_run )  THEN
             IF ( nesting_mode == 'vertical' )  THEN
                IF (bc_lr /= 'cyclic' .OR. bc_ns /= 'cyclic' )  THEN
                   WRITE ( message_string, *) 'bc_lr and bc_ns were set to ,', &
                        'cyclic for vertical nesting'
                   CALL message( 'parin', 'PA0428', 0, 0, 0, 6, 0 )
                   bc_lr   = 'cyclic'
                   bc_ns   = 'cyclic'
                ENDIF
                IF ( child_domain )  THEN
                   bc_uv_t  = 'nested'
                   bc_pt_t  = 'nested'
                   bc_q_t   = 'nested'
                   bc_s_t   = 'nested'
                   bc_cs_t  = 'nested'
                   bc_p_t   = 'neumann'  
                ENDIF
!
!--          For other nesting modes only set boundary conditions for 
!--          nested domains.
             ELSE 
                IF ( child_domain )  THEN
                   bc_lr    = 'nested'
                   bc_ns    = 'nested'
                   bc_uv_t  = 'nested'
                   bc_pt_t  = 'nested'
                   bc_q_t   = 'nested'
                   bc_s_t   = 'nested'
                   bc_cs_t  = 'nested'
                   bc_p_t   = 'neumann'
                ENDIF
             ENDIF
          ENDIF
!
!--       Set boundary conditions also in case the model is offline-nested in
!--       larger-scale models.
          IF ( nesting_offline )  THEN
             bc_lr    = 'nesting_offline'
             bc_ns    = 'nesting_offline'
             bc_uv_t  = 'nesting_offline'
             bc_pt_t  = 'nesting_offline'
             bc_q_t   = 'nesting_offline'
           !  bc_s_t   = 'nesting_offline'  ! scalar boundary condition is not clear yet
           !  bc_cs_t  = 'nesting_offline'  ! same for chemical species
             bc_p_t   = 'neumann'
          ENDIF

!         
!--       In case of nested runs, make sure that initializing_actions =
!--       'set_constant_profiles' even though the constant-profiles 
!--       initializations for the prognostic variables will be overwritten 
!--       by pmci_child_initialize and pmci_parent_initialize. This is, 
!--       however, important e.g. to make sure that diagnostic variables 
!--       are set properly. An exception is made in case of restart runs and
!--       if user decides to do everything by its own.
          IF ( child_domain  .AND.  .NOT. (                                    &
               TRIM( initializing_actions ) == 'read_restart_data'      .OR.   &
               TRIM( initializing_actions ) == 'set_constant_profiles'  .OR.   &
               TRIM( initializing_actions ) == 'by_user' ) )  THEN
             message_string = 'initializing_actions = ' //                     &
                              TRIM( initializing_actions ) // ' has been ' //  &
                              'changed to set_constant_profiles in child ' //  &
                              'domain.' 
             CALL message( 'parin', 'PA0492', 0, 0, 0, 6, 0 )

             initializing_actions = 'set_constant_profiles'
          ENDIF            
!
!--       Check validity of lateral boundary conditions. This has to be done
!--       here because they are already used in init_pegrid and init_grid and
!--       therefore cannot be check in check_parameters
          IF ( bc_lr /= 'cyclic'  .AND.  bc_lr /= 'dirichlet/radiation'  .AND. &
               bc_lr /= 'radiation/dirichlet'  .AND.  bc_lr /= 'nested'  .AND. &
               bc_lr /= 'nesting_offline' )  THEN
             message_string = 'unknown boundary condition: bc_lr = "' // &
                              TRIM( bc_lr ) // '"'
             CALL message( 'parin', 'PA0049', 1, 2, 0, 6, 0 )
          ENDIF
          IF ( bc_ns /= 'cyclic'  .AND.  bc_ns /= 'dirichlet/radiation'  .AND. &
               bc_ns /= 'radiation/dirichlet'  .AND.  bc_ns /= 'nested'  .AND. &
               bc_ns /= 'nesting_offline' )  THEN
             message_string = 'unknown boundary condition: bc_ns = "' // &
                              TRIM( bc_ns ) // '"'
             CALL message( 'parin', 'PA0050', 1, 2, 0, 6, 0 )
          ENDIF
!
!--       Set internal variables used for speed optimization in if clauses
          IF ( bc_lr /= 'cyclic' )               bc_lr_cyc    = .FALSE.
          IF ( bc_lr == 'dirichlet/radiation' )  bc_lr_dirrad = .TRUE.
          IF ( bc_lr == 'radiation/dirichlet' )  bc_lr_raddir = .TRUE.
          IF ( bc_ns /= 'cyclic' )               bc_ns_cyc    = .FALSE.
          IF ( bc_ns == 'dirichlet/radiation' )  bc_ns_dirrad = .TRUE.
          IF ( bc_ns == 'radiation/dirichlet' )  bc_ns_raddir = .TRUE.
!
!--       Radiation-Dirichlet conditions are allowed along one of the horizontal directions only.
!--       In general, such conditions along x and y may work, but require a) some code changes
!--       (e.g. concerning allocation of c_u, c_v ... arrays), and b) a careful model setup by the
!--       user, in order to guarantee that there is a clearly defined outflow at two sides.
!--       Otherwise, the radiation condition may produce wrong results.
          IF ( ( bc_lr_dirrad .OR. bc_lr_raddir )  .AND.  ( bc_ns_dirrad .OR. bc_ns_raddir ) )  THEN
             message_string = 'bc_lr = "' // TRIM( bc_lr ) // '" and bc_ns = "' //                 &
                              TRIM( bc_ns ) // '" are not allowed to be set at the same time'
             CALL message( 'parin', 'PA0589', 1, 2, 0, 6, 0 )
          ENDIF
!
!--       Check in case of initial run, if the grid point numbers are well
!--       defined and allocate some arrays which are already needed in
!--       init_pegrid or check_parameters. During restart jobs, these arrays
!--       will be allocated in rrd_global. All other arrays are allocated
!--       in init_3d_model.
          IF ( TRIM( initializing_actions ) /= 'read_restart_data' )  THEN

             IF ( nx <= 0 )  THEN
                WRITE( message_string, * ) 'no value or wrong value given',    &
                                           ' for nx: nx=', nx
                CALL message( 'parin', 'PA0273', 1, 2, 0, 6, 0 )
             ENDIF
             IF ( ny <= 0 )  THEN
                WRITE( message_string, * ) 'no value or wrong value given',    &
                                           ' for ny: ny=', ny
                CALL message( 'parin', 'PA0274', 1, 2, 0, 6, 0 )
             ENDIF
             IF ( nz <= 0 )  THEN
                WRITE( message_string, * ) 'no value or wrong value given',    &
                                           ' for nz: nz=', nz
                CALL message( 'parin', 'PA0275', 1, 2, 0, 6, 0 )
             ENDIF

!
!--          As a side condition, routine flow_statistics require at least 14
!--          vertical grid levels (they are used to store time-series data)
!>           @todo   Remove this restriction
             IF ( nz < 14 )  THEN
                WRITE( message_string, * ) 'nz >= 14 is required'
                CALL message( 'parin', 'PA0362', 1, 2, 0, 6, 0 )
             ENDIF

!
!--          ATTENTION: in case of changes to the following statement please
!--                  also check the allocate statement in routine rrd_global
             ALLOCATE( pt_init(0:nz+1), q_init(0:nz+1), s_init(0:nz+1),        &
                       ref_state(0:nz+1), sa_init(0:nz+1), ug(0:nz+1),         &
                       u_init(0:nz+1), v_init(0:nz+1), vg(0:nz+1),             &
                    hom(0:nz+1,2,pr_palm+max_pr_user+max_pr_cs+max_pr_salsa,0:statistic_regions),  &
                    hom_sum(0:nz+1,pr_palm+max_pr_user+max_pr_cs+max_pr_salsa,0:statistic_regions) )

             hom = 0.0_wp

          ENDIF

!
!--       NAMELIST-file is not needed anymore
          CALL close_file( 11 )

       ENDIF
#if defined( __parallel )
       CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
#endif
    ENDDO

    CALL location_message( 'reading NAMELIST parameters from PARIN', 'finished' )

 END SUBROUTINE parin
