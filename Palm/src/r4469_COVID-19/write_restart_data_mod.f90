!> @file write_restart_data_mod.f90
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
! $Id: write_restart_data_mod.f90 4360 2020-01-07 11:25:50Z suehring $
! - Move 2-m potential temperature to diagnostic_output_quantities.
! 
! 4301 2019-11-22 12:09:09Z oliver.maas
! removed recycling_yshift
! 
! 4227 2019-09-10 18:04:34Z gronemeier
! implement new palm_date_time_mod and increased binary version
! 
! 4146 2019-08-07 07:47:36Z gronemeier
! Corrected "Former revisions" section
! 
! 4101 2019-07-17 15:14:26Z gronemeier
! remove old_dt
! 
! 4039 2019-06-18 10:32:41Z suehring
! Modularize diagnostic output
! 
! 4017 2019-06-06 12:16:46Z schwenkel
! output of turbulence intensity added
! 
! 3988 2019-05-22 11:32:37Z kanani
! + time_virtual_measurement (to enable steering of output interval)
! 
! 3936 2019-04-26 15:38:02Z kanani
! Enable time-averaged output of theta_2m* with restarts
! 
! 3668 2019-01-14 12:49:24Z maronga
! Implementation of the PALM module interface
! 
! 2894 2018-03-15 09:17:58Z Giersch
! Initial revision 
!
!
! Description:
! ------------
!> Writes restart data into binary file(s) for restart runs.
!------------------------------------------------------------------------------!
 MODULE write_restart_data_mod


    USE arrays_3d,                                                          &
        ONLY:  inflow_damping_factor, mean_inflow_profiles, pt_init,        &
               q_init, ref_state, s_init, u_init, ug, v_init, vg,           &
               e, kh, km, p, pt, q, ql, s, u, u_m_l, u_m_n, u_m_r,          &
               u_m_s, v, v_m_l, v_m_n, v_m_r, v_m_s, vpt, w, w_m_l, w_m_n,  &
               w_m_r, w_m_s

    USE averaging

    USE control_parameters

    USE grid_variables,                                                        &
        ONLY:  dx, dy

    USE indices,                                                               &
        ONLY:  nx, nxl, nxr, ny, nys, nyn, nz, nzb, nzt

    USE kinds

    USE model_1d_mod,                                                          &
        ONLY:  damp_level_1d, dt_pr_1d, dt_run_control_1d, end_time_1d

    USE module_interface,                                                      &
        ONLY:  module_interface_wrd_global,                                    &
               module_interface_wrd_local

    USE netcdf_interface,                                                      &
        ONLY:  netcdf_precision, output_for_t0

    USE pegrid,                                                                &
        ONLY:  collective_wait, hor_index_bounds, myid, numprocs

    USE radiation_model_mod,                                                   &
        ONLY:  time_radiation

    USE random_function_mod,                                                   &
        ONLY:  random_iv, random_iy

    USE random_generator_parallel,                                             &
        ONLY:  id_random_array, seq_random_array

    USE spectra_mod,                                                           &
        ONLY:  average_count_sp, spectrum_x, spectrum_y

    USE statistics,                                                            &
        ONLY:  statistic_regions, hom, hom_sum, u_max, u_max_ijk, v_max,       &
               v_max_ijk, w_max, w_max_ijk, z_i

    USE surface_mod,                                                           &
        ONLY :  surface_wrd_local

    USE vertical_nesting_mod,                                                  &
        ONLY:  vnest_init

    USE virtual_measurement_mod,                                               &
        ONLY:  time_virtual_measurement


    IMPLICIT NONE


    INTERFACE wrd_global
       MODULE PROCEDURE wrd_global
    END INTERFACE wrd_global

    INTERFACE wrd_local
       MODULE PROCEDURE wrd_local
    END INTERFACE wrd_local


    PUBLIC wrd_local, wrd_global


 CONTAINS


! Description:
! ------------
!> Global data of control variables and arrays is written out for 
!> restarts (binary format).
!> This information is only written to the file opened by PE0.
!------------------------------------------------------------------------------!
    SUBROUTINE wrd_global


       CHARACTER (LEN=10)  ::  binary_version_global   !<


       binary_version_global = '4.9'

       CALL wrd_write_string( 'binary_version_global' )
       WRITE ( 14 )  binary_version_global

       CALL wrd_write_string( 'numprocs' )
       WRITE ( 14 )  numprocs

       CALL wrd_write_string( 'hor_index_bounds' ) 
       WRITE ( 14 )  hor_index_bounds

       CALL wrd_write_string( 'nz' ) 
       WRITE ( 14 )  nz

       CALL wrd_write_string( 'max_pr_user' ) 
       WRITE ( 14 )  max_pr_user

       CALL wrd_write_string( 'statistic_regions' ) 
       WRITE ( 14 )  statistic_regions

!
!-- Caution: After changes in the following parameter-list, the
!-- -------  version number stored in the variable binary_version_global has to
!--          be increased. The same changes must also be done in the parameter-
!--          list in rrd_global.

       CALL wrd_write_string( 'advected_distance_x' ) 
       WRITE ( 14 )  advected_distance_x

       CALL wrd_write_string( 'advected_distance_y' ) 
       WRITE ( 14 )  advected_distance_y

       CALL wrd_write_string( 'alpha_surface' ) 
       WRITE ( 14 )  alpha_surface 

       CALL wrd_write_string( 'average_count_pr' ) 
       WRITE ( 14 )  average_count_pr

       CALL wrd_write_string( 'average_count_sp' ) 
       WRITE ( 14 )  average_count_sp

       CALL wrd_write_string( 'average_count_3d' ) 
       WRITE ( 14 )  average_count_3d

       CALL wrd_write_string( 'bc_e_b' ) 
       WRITE ( 14 )  bc_e_b

       CALL wrd_write_string( 'bc_lr' ) 
       WRITE ( 14 )  bc_lr

       CALL wrd_write_string( 'bc_ns' ) 
       WRITE ( 14 )  bc_ns

       CALL wrd_write_string( 'bc_p_b' ) 
       WRITE ( 14 )  bc_p_b

       CALL wrd_write_string( 'bc_p_t' ) 
       WRITE ( 14 )  bc_p_t

       CALL wrd_write_string( 'bc_pt_b' ) 
       WRITE ( 14 )  bc_pt_b

       CALL wrd_write_string( 'bc_pt_t' ) 
       WRITE ( 14 )  bc_pt_t

       CALL wrd_write_string( 'bc_pt_t_val' ) 
       WRITE ( 14 )  bc_pt_t_val

       CALL wrd_write_string( 'bc_q_b' ) 
       WRITE ( 14 )  bc_q_b

       CALL wrd_write_string( 'bc_q_t' ) 
       WRITE ( 14 )  bc_q_t

       CALL wrd_write_string( 'bc_q_t_val' ) 
       WRITE ( 14 )  bc_q_t_val

       CALL wrd_write_string( 'bc_s_b' ) 
       WRITE ( 14 )  bc_s_b

       CALL wrd_write_string( 'bc_s_t' ) 
       WRITE ( 14 )  bc_s_t

       CALL wrd_write_string( 'bc_uv_b' ) 
       WRITE ( 14 )  bc_uv_b

       CALL wrd_write_string( 'bc_uv_t' ) 
       WRITE ( 14 )  bc_uv_t

       CALL wrd_write_string( 'building_height' ) 
       WRITE ( 14 )  building_height

       CALL wrd_write_string( 'building_length_x' ) 
       WRITE ( 14 )  building_length_x

       CALL wrd_write_string( 'building_length_y' ) 
       WRITE ( 14 )  building_length_y

       CALL wrd_write_string( 'building_wall_left' ) 
       WRITE ( 14 )  building_wall_left

       CALL wrd_write_string( 'building_wall_south' ) 
       WRITE ( 14 )  building_wall_south

       CALL wrd_write_string( 'call_psolver_at_all_substeps' ) 
       WRITE ( 14 )  call_psolver_at_all_substeps

       CALL wrd_write_string( 'canyon_height' ) 
       WRITE ( 14 )  canyon_height 

       CALL wrd_write_string( 'canyon_wall_left' ) 
       WRITE ( 14 )  canyon_wall_left

       CALL wrd_write_string( 'canyon_wall_south' ) 
       WRITE ( 14 )  canyon_wall_south

       CALL wrd_write_string( 'canyon_width_x' ) 
       WRITE ( 14 )  canyon_width_x

       CALL wrd_write_string( 'canyon_width_y' ) 
       WRITE ( 14 )  canyon_width_y

       CALL wrd_write_string( 'cfl_factor' ) 
       WRITE ( 14 )  cfl_factor

       CALL wrd_write_string( 'cloud_droplets' ) 
       WRITE ( 14 )  cloud_droplets

       CALL wrd_write_string( 'collective_wait' ) 
       WRITE ( 14 )  collective_wait

       CALL wrd_write_string( 'conserve_volume_flow' ) 
       WRITE ( 14 )  conserve_volume_flow

       CALL wrd_write_string( 'conserve_volume_flow_mode' ) 
       WRITE ( 14 )  conserve_volume_flow_mode

       CALL wrd_write_string( 'constant_flux_layer' ) 
       WRITE ( 14 )  constant_flux_layer

       CALL wrd_write_string( 'coupling_start_time' ) 
       WRITE ( 14 )  coupling_start_time

       CALL wrd_write_string( 'current_timestep_number' ) 
       WRITE ( 14 )  current_timestep_number

       CALL wrd_write_string( 'cycle_mg' ) 
       WRITE ( 14 )  cycle_mg

       CALL wrd_write_string( 'damp_level_1d' ) 
       WRITE ( 14 )  damp_level_1d

       CALL wrd_write_string( 'origin_date_time' ) 
       WRITE ( 14 )  origin_date_time

       CALL wrd_write_string( 'dissipation_1d' ) 
       WRITE ( 14 )  dissipation_1d

       CALL wrd_write_string( 'do2d_xy_time_count' ) 
       WRITE ( 14 )  do2d_xy_time_count

       CALL wrd_write_string( 'do2d_xz_time_count' ) 
       WRITE ( 14 )  do2d_xz_time_count

       CALL wrd_write_string( 'do2d_yz_time_count' ) 
       WRITE ( 14 )  do2d_yz_time_count

       CALL wrd_write_string( 'do3d_time_count' ) 
       WRITE ( 14 )  do3d_time_count

       CALL wrd_write_string( 'dp_external' ) 
       WRITE ( 14 )  dp_external

       CALL wrd_write_string( 'dp_level_b' ) 
       WRITE ( 14 )  dp_level_b

       CALL wrd_write_string( 'dp_smooth' ) 
       WRITE ( 14 )  dp_smooth

       CALL wrd_write_string( 'dpdxy' ) 
       WRITE ( 14 )  dpdxy

       CALL wrd_write_string( 'dt_3d' ) 
       WRITE ( 14 )  dt_3d

       CALL wrd_write_string( 'dt_pr_1d' ) 
       WRITE ( 14 )  dt_pr_1d

       CALL wrd_write_string( 'dt_run_control_1d' ) 
       WRITE ( 14 )  dt_run_control_1d

       CALL wrd_write_string( 'dx' ) 
       WRITE ( 14 )  dx

       CALL wrd_write_string( 'dy' ) 
       WRITE ( 14 )  dy

       CALL wrd_write_string( 'dz' ) 
       WRITE ( 14 )  dz
       
       CALL wrd_write_string( 'dz_max' ) 
       WRITE ( 14 )  dz_max

       CALL wrd_write_string( 'dz_stretch_factor' ) 
       WRITE ( 14 )  dz_stretch_factor
       
       CALL wrd_write_string( 'dz_stretch_factor_array' ) 
       WRITE ( 14 )  dz_stretch_factor_array
       
       CALL wrd_write_string( 'dz_stretch_level' ) 
       WRITE ( 14 )  dz_stretch_level

       CALL wrd_write_string( 'dz_stretch_level_end' ) 
       WRITE ( 14 )  dz_stretch_level_end
       
       CALL wrd_write_string( 'dz_stretch_level_start' ) 
       WRITE ( 14 )  dz_stretch_level_start
       
       CALL wrd_write_string( 'e_min' ) 
       WRITE ( 14 )  e_min

       CALL wrd_write_string( 'end_time_1d' ) 
       WRITE ( 14 )  end_time_1d

       CALL wrd_write_string( 'fft_method' ) 
       WRITE ( 14 )  fft_method

       CALL wrd_write_string( 'first_call_lpm' ) 
       WRITE ( 14 )  first_call_lpm

       CALL wrd_write_string( 'galilei_transformation' ) 
       WRITE ( 14 )  galilei_transformation

       CALL wrd_write_string( 'hom' ) 
       WRITE ( 14 )  hom

       CALL wrd_write_string( 'hom_sum' ) 
       WRITE ( 14 )  hom_sum

       CALL wrd_write_string( 'humidity' ) 
       WRITE ( 14 )  humidity

       IF ( ALLOCATED( inflow_damping_factor ) )  THEN
          CALL wrd_write_string( 'inflow_damping_factor' ) 
          WRITE ( 14 )  inflow_damping_factor
       ENDIF

       CALL wrd_write_string( 'inflow_damping_height' ) 
       WRITE ( 14 )  inflow_damping_height

       CALL wrd_write_string( 'inflow_damping_width' ) 
       WRITE ( 14 )  inflow_damping_width

       CALL wrd_write_string( 'inflow_disturbance_begin' ) 
       WRITE ( 14 )  inflow_disturbance_begin

       CALL wrd_write_string( 'inflow_disturbance_end' ) 
       WRITE ( 14 )  inflow_disturbance_end

       CALL wrd_write_string( 'km_constant' ) 
       WRITE ( 14 )  km_constant

       CALL wrd_write_string( 'large_scale_forcing' ) 
       WRITE ( 14 )  large_scale_forcing

       CALL wrd_write_string( 'large_scale_subsidence' ) 
       WRITE ( 14 )  large_scale_subsidence

       CALL wrd_write_string( 'latitude' ) 
       WRITE ( 14 )  latitude

       CALL wrd_write_string( 'longitude' )
       WRITE ( 14 )  longitude 

       CALL wrd_write_string( 'loop_optimization' ) 
       WRITE ( 14 )  loop_optimization

       CALL wrd_write_string( 'masking_method' ) 
       WRITE ( 14 )  masking_method

       IF ( ALLOCATED( mean_inflow_profiles ) )  THEN
          CALL wrd_write_string( 'mean_inflow_profiles' ) 
          WRITE ( 14 )  mean_inflow_profiles
       ENDIF

       CALL wrd_write_string( 'mg_cycles' ) 
       WRITE ( 14 )  mg_cycles

       CALL wrd_write_string( 'mg_switch_to_pe0_level' ) 
       WRITE ( 14 )  mg_switch_to_pe0_level

       CALL wrd_write_string( 'mixing_length_1d' ) 
       WRITE ( 14 )  mixing_length_1d

       CALL wrd_write_string( 'momentum_advec' ) 
       WRITE ( 14 )  momentum_advec

       CALL wrd_write_string( 'netcdf_precision' ) 
       WRITE ( 14 )  netcdf_precision

       CALL wrd_write_string( 'neutral' ) 
       WRITE ( 14 )  neutral

       CALL wrd_write_string( 'ngsrb' ) 
       WRITE ( 14 )  ngsrb

       CALL wrd_write_string( 'nsor' ) 
       WRITE ( 14 )  nsor

       CALL wrd_write_string( 'nsor_ini' ) 
       WRITE ( 14 )  nsor_ini

       CALL wrd_write_string( 'nudging' ) 
       WRITE ( 14 )  nudging

       CALL wrd_write_string( 'num_leg' ) 
       WRITE ( 14 )  num_leg

       CALL wrd_write_string( 'nx' ) 
       WRITE ( 14 )  nx

       CALL wrd_write_string( 'ny' ) 
       WRITE ( 14 )  ny

       CALL wrd_write_string( 'ocean_mode' ) 
       WRITE ( 14 )  ocean_mode

       CALL wrd_write_string( 'omega' ) 
       WRITE ( 14 )  omega

       CALL wrd_write_string( 'omega_sor' ) 
       WRITE ( 14 )  omega_sor

       CALL wrd_write_string( 'output_for_t0' ) 
       WRITE ( 14 )  output_for_t0

       CALL wrd_write_string( 'passive_scalar' ) 
       WRITE ( 14 )  passive_scalar

       CALL wrd_write_string( 'prandtl_number' ) 
       WRITE ( 14 )  prandtl_number

       CALL wrd_write_string( 'psolver' )
       WRITE ( 14 )  psolver 

       CALL wrd_write_string( 'pt_damping_factor' ) 
       WRITE ( 14 )  pt_damping_factor

       CALL wrd_write_string( 'pt_damping_width' ) 
       WRITE ( 14 )  pt_damping_width

       CALL wrd_write_string( 'pt_init' ) 
       WRITE ( 14 )  pt_init

       CALL wrd_write_string( 'pt_reference' ) 
       WRITE ( 14 )  pt_reference

       CALL wrd_write_string( 'pt_surface' ) 
       WRITE ( 14 )  pt_surface

       CALL wrd_write_string( 'pt_surface_initial_change' ) 
       WRITE ( 14 )  pt_surface_initial_change

       CALL wrd_write_string( 'pt_vertical_gradient' ) 
       WRITE ( 14 )  pt_vertical_gradient

       CALL wrd_write_string( 'pt_vertical_gradient_level' ) 
       WRITE ( 14 )  pt_vertical_gradient_level

       CALL wrd_write_string( 'pt_vertical_gradient_level_ind' ) 
       WRITE ( 14 )  pt_vertical_gradient_level_ind

       CALL wrd_write_string( 'q_init' ) 
       WRITE ( 14 )  q_init

       CALL wrd_write_string( 'q_surface' ) 
       WRITE ( 14 )  q_surface

       CALL wrd_write_string( 'q_surface_initial_change' ) 
       WRITE ( 14 )  q_surface_initial_change

       CALL wrd_write_string( 'q_vertical_gradient' ) 
       WRITE ( 14 )  q_vertical_gradient

       CALL wrd_write_string( 'q_vertical_gradient_level' ) 
       WRITE ( 14 )  q_vertical_gradient_level

       CALL wrd_write_string( 'q_vertical_gradient_level_ind' ) 
       WRITE ( 14 )  q_vertical_gradient_level_ind

       CALL wrd_write_string( 'random_generator' ) 
       WRITE ( 14 )  random_generator

       CALL wrd_write_string( 'random_heatflux' ) 
       WRITE ( 14 )  random_heatflux

       CALL wrd_write_string( 'rans_mode' ) 
       WRITE ( 14 )  rans_mode

       CALL wrd_write_string( 'rayleigh_damping_factor' ) 
       WRITE ( 14 )  rayleigh_damping_factor

       CALL wrd_write_string( 'rayleigh_damping_height' ) 
       WRITE ( 14 )  rayleigh_damping_height

       CALL wrd_write_string( 'recycling_width' ) 
       WRITE ( 14 )  recycling_width

       CALL wrd_write_string( 'ref_state' ) 
       WRITE ( 14 )  ref_state

       CALL wrd_write_string( 'reference_state' ) 
       WRITE ( 14 )  reference_state

       CALL wrd_write_string( 'residual_limit' ) 
       WRITE ( 14 )  residual_limit

       CALL wrd_write_string( 'roughness_length' ) 
       WRITE ( 14 )  roughness_length

       CALL wrd_write_string( 'run_coupled' ) 
       WRITE ( 14 )  run_coupled

       CALL wrd_write_string( 'runnr' ) 
       WRITE ( 14 )  runnr

       CALL wrd_write_string( 's_init' ) 
       WRITE ( 14 )  s_init

       CALL wrd_write_string( 's_surface' ) 
       WRITE ( 14 )  s_surface

       CALL wrd_write_string( 's_surface_initial_change' ) 
       WRITE ( 14 )  s_surface_initial_change

       CALL wrd_write_string( 's_vertical_gradient' ) 
       WRITE ( 14 )  s_vertical_gradient

       CALL wrd_write_string( 's_vertical_gradient_level' ) 
       WRITE ( 14 )  s_vertical_gradient_level

       CALL wrd_write_string( 's_vertical_gradient_level_ind' ) 
       WRITE ( 14 )  s_vertical_gradient_level_ind

       CALL wrd_write_string( 'scalar_advec' ) 
       WRITE ( 14 )  scalar_advec

       CALL wrd_write_string( 'simulated_time' ) 
       WRITE ( 14 )  simulated_time

       IF ( ALLOCATED( spectrum_x ) )  THEN
          CALL wrd_write_string( 'spectrum_x' ) 
          WRITE ( 14 )  spectrum_x
          CALL wrd_write_string( 'spectrum_y' ) 
          WRITE ( 14 )  spectrum_y
       ENDIF

       CALL wrd_write_string( 'spinup_time ' ) 
       WRITE ( 14 )  spinup_time 

       CALL wrd_write_string( 'surface_heatflux' ) 
       WRITE ( 14 )  surface_heatflux

       CALL wrd_write_string( 'surface_pressure' ) 
       WRITE ( 14 )  surface_pressure

       CALL wrd_write_string( 'surface_scalarflux' ) 
       WRITE ( 14 )  surface_scalarflux

       CALL wrd_write_string( 'surface_waterflux' ) 
       WRITE ( 14 )  surface_waterflux

       CALL wrd_write_string( 'time_coupling' ) 
       WRITE ( 14 )  time_coupling

       CALL wrd_write_string( 'time_disturb' ) 
       WRITE ( 14 )  time_disturb

       CALL wrd_write_string( 'time_do2d_xy' ) 
       WRITE ( 14 )  time_do2d_xy

       CALL wrd_write_string( 'time_do2d_xz' ) 
       WRITE ( 14 )  time_do2d_xz

       CALL wrd_write_string( 'time_do2d_yz' ) 
       WRITE ( 14 )  time_do2d_yz

       CALL wrd_write_string( 'time_do3d' ) 
       WRITE ( 14 )  time_do3d

       CALL wrd_write_string( 'time_do_av' ) 
       WRITE ( 14 )  time_do_av

       CALL wrd_write_string( 'time_do_sla' ) 
       WRITE ( 14 )  time_do_sla

       CALL wrd_write_string( 'time_domask' ) 
       WRITE ( 14 )  time_domask

       CALL wrd_write_string( 'time_dopr' ) 
       WRITE ( 14 )  time_dopr

       CALL wrd_write_string( 'time_dopr_av' ) 
       WRITE ( 14 )  time_dopr_av

       CALL wrd_write_string( 'time_dopr_listing' ) 
       WRITE ( 14 )  time_dopr_listing

       CALL wrd_write_string( 'time_dopts' ) 
       WRITE ( 14 )  time_dopts

       CALL wrd_write_string( 'time_dosp' ) 
       WRITE ( 14 )  time_dosp

       CALL wrd_write_string( 'time_dots' ) 
       WRITE ( 14 )  time_dots

       CALL wrd_write_string( 'time_radiation' ) 
       WRITE ( 14 )  time_radiation

       CALL wrd_write_string( 'time_restart' ) 
       WRITE ( 14 )  time_restart

       CALL wrd_write_string( 'time_run_control' ) 
       WRITE ( 14 )  time_run_control

       CALL wrd_write_string( 'time_since_reference_point' ) 
       WRITE ( 14 )  time_since_reference_point

       CALL wrd_write_string( 'time_virtual_measurement' ) 
       WRITE ( 14 )  time_virtual_measurement

       CALL wrd_write_string( 'timestep_scheme' ) 
       WRITE ( 14 )  timestep_scheme

       CALL wrd_write_string( 'top_heatflux' ) 
       WRITE ( 14 )  top_heatflux

       CALL wrd_write_string( 'top_momentumflux_u' ) 
       WRITE ( 14 )  top_momentumflux_u

       CALL wrd_write_string( 'top_momentumflux_v' ) 
       WRITE ( 14 )  top_momentumflux_v

       CALL wrd_write_string( 'top_scalarflux' ) 
       WRITE ( 14 )  top_scalarflux

       CALL wrd_write_string( 'topography' ) 
       WRITE ( 14 )  topography

       CALL wrd_write_string( 'topography_grid_convention' ) 
       WRITE ( 14 )  topography_grid_convention

       CALL wrd_write_string( 'tsc' ) 
       WRITE ( 14 )  tsc

       CALL wrd_write_string( 'tunnel_height' ) 
       WRITE ( 14 )  tunnel_height

       CALL wrd_write_string( 'tunnel_length' ) 
       WRITE ( 14 )  tunnel_length

       CALL wrd_write_string( 'tunnel_wall_depth' ) 
       WRITE ( 14 )  tunnel_wall_depth

       CALL wrd_write_string( 'tunnel_width_x' ) 
       WRITE ( 14 )  tunnel_width_x

       CALL wrd_write_string( 'tunnel_width_y' ) 
       WRITE ( 14 )  tunnel_width_y

       CALL wrd_write_string( 'turbulence_closure' ) 
       WRITE ( 14 )  turbulence_closure

       CALL wrd_write_string( 'turbulent_inflow' ) 
       WRITE ( 14 )  turbulent_inflow

       CALL wrd_write_string( 'u_bulk' ) 
       WRITE ( 14 )  u_bulk

       CALL wrd_write_string( 'u_init' ) 
       WRITE ( 14 )  u_init

       CALL wrd_write_string( 'u_max' ) 
       WRITE ( 14 )  u_max

       CALL wrd_write_string( 'u_max_ijk' ) 
       WRITE ( 14 )  u_max_ijk

       CALL wrd_write_string( 'ug' ) 
       WRITE ( 14 )  ug

       CALL wrd_write_string( 'ug_surface' ) 
       WRITE ( 14 )  ug_surface

       CALL wrd_write_string( 'ug_vertical_gradient' ) 
       WRITE ( 14 )  ug_vertical_gradient

       CALL wrd_write_string( 'ug_vertical_gradient_level' ) 
       WRITE ( 14 )  ug_vertical_gradient_level

       CALL wrd_write_string( 'ug_vertical_gradient_level_ind' ) 
       WRITE ( 14 )  ug_vertical_gradient_level_ind

       CALL wrd_write_string( 'use_surface_fluxes' ) 
       WRITE ( 14 )  use_surface_fluxes

       CALL wrd_write_string( 'use_top_fluxes' ) 
       WRITE ( 14 )  use_top_fluxes

       CALL wrd_write_string( 'use_ug_for_galilei_tr' ) 
       WRITE ( 14 )  use_ug_for_galilei_tr

       CALL wrd_write_string( 'use_upstream_for_tke' ) 
       WRITE ( 14 )  use_upstream_for_tke

       CALL wrd_write_string( 'v_bulk' ) 
       WRITE ( 14 )  v_bulk

       CALL wrd_write_string( 'v_init' ) 
       WRITE ( 14 )  v_init

       CALL wrd_write_string( 'v_max' ) 
       WRITE ( 14 )  v_max

       CALL wrd_write_string( 'v_max_ijk' ) 
       WRITE ( 14 )  v_max_ijk

       CALL wrd_write_string( 'vg' ) 
       WRITE ( 14 )  vg

       CALL wrd_write_string( 'vg_surface' ) 
       WRITE ( 14 )  vg_surface

       CALL wrd_write_string( 'vg_vertical_gradient' ) 
       WRITE ( 14 ) vg_vertical_gradient

       CALL wrd_write_string( 'vg_vertical_gradient_level' ) 
       WRITE ( 14 )  vg_vertical_gradient_level

       CALL wrd_write_string( 'vg_vertical_gradient_level_ind' ) 
       WRITE ( 14 )  vg_vertical_gradient_level_ind

       CALL wrd_write_string( 'virtual_flight' ) 
       WRITE ( 14 )  virtual_flight

       CALL wrd_write_string( 'vnest_init' ) 
       WRITE ( 14 )  vnest_init

       CALL wrd_write_string( 'volume_flow_area' ) 
       WRITE ( 14 )  volume_flow_area

       CALL wrd_write_string( 'volume_flow_initial' ) 
       WRITE ( 14 )  volume_flow_initial

       CALL wrd_write_string( 'subs_vertical_gradient' ) 
       WRITE ( 14 )  subs_vertical_gradient

       CALL wrd_write_string( 'subs_vertical_gradient_level' ) 
       WRITE ( 14 )  subs_vertical_gradient_level

       CALL wrd_write_string( 'subs_vertical_gradient_level_i' ) 
       WRITE ( 14 )  subs_vertical_gradient_level_i

       CALL wrd_write_string( 'w_max' ) 
       WRITE ( 14 )  w_max

       CALL wrd_write_string( 'w_max_ijk' ) 
       WRITE ( 14 )  w_max_ijk

       CALL wrd_write_string( 'wall_adjustment' ) 
       WRITE ( 14 )  wall_adjustment

       CALL wrd_write_string( 'wall_heatflux' ) 
       WRITE ( 14 )  wall_heatflux

       CALL wrd_write_string( 'wall_humidityflux' ) 
       WRITE ( 14 )  wall_humidityflux

       CALL wrd_write_string( 'wall_scalarflux' ) 
       WRITE ( 14 )  wall_scalarflux

       CALL wrd_write_string( 'y_shift' ) 
       WRITE ( 14 )  y_shift

       CALL wrd_write_string( 'z0h_factor' ) 
       WRITE ( 14 )  z0h_factor

       CALL wrd_write_string( 'zeta_max' ) 
       WRITE ( 14 )  zeta_max

       CALL wrd_write_string( 'zeta_min' ) 
       WRITE ( 14 )  zeta_min

       CALL wrd_write_string( 'z_i' ) 
       WRITE ( 14 )  z_i

!
!--    Write restart data of the other modules
       CALL module_interface_wrd_global


    END SUBROUTINE wrd_global


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Processor specific data of variables and arrays is written out for 
!> restarts (binary format).
!> This information is written to the file opened by each PE.
!------------------------------------------------------------------------------!
    SUBROUTINE wrd_local


       CHARACTER (LEN=10) ::  binary_version_local   !< 

!
!-- Write arrays.
       binary_version_local = '4.7'

       CALL wrd_write_string( 'binary_version_local' )
       WRITE ( 14 )  binary_version_local

       WRITE ( 14 )  numprocs, myid, nxl, nxr, nys, nyn, nzb, nzt

!
!-- Attention: After changes to the following output commands the version number
!-- ---------  of the variable binary_version_local must be changed!
!--            Also, the list of arrays to be read in rrd_local must be
!--            adjusted accordingly.
       CALL wrd_write_string( 'e' ) 
       WRITE ( 14 )  e

       IF ( ALLOCATED( e_av ) )  THEN
          CALL wrd_write_string( 'e_av' ) 
          WRITE ( 14 )  e_av
       ENDIF

       CALL wrd_write_string( 'kh' ) 
       WRITE ( 14 )  kh
  

       IF ( ALLOCATED( kh_av ) )  THEN
          CALL wrd_write_string( 'kh_av' ) 
          WRITE ( 14 )  kh_av
       ENDIF

       CALL wrd_write_string( 'km' ) 
       WRITE ( 14 )  km

       IF ( ALLOCATED( km_av ) )  THEN
          CALL wrd_write_string( 'km_av' ) 
          WRITE ( 14 )  km_av
       ENDIF

       IF ( ALLOCATED( lpt_av ) )  THEN
          CALL wrd_write_string( 'lpt_av' ) 
          WRITE ( 14 )  lpt_av
       ENDIF

       IF ( ALLOCATED( lwp_av ) )  THEN
          CALL wrd_write_string( 'lwp_av' ) 
          WRITE ( 14 )  lwp_av
       ENDIF

       CALL wrd_write_string( 'p' ) 
       WRITE ( 14 )  p

       IF ( ALLOCATED( p_av ) )  THEN
          CALL wrd_write_string( 'p_av' ) 
          WRITE ( 14 )  p_av
       ENDIF

       IF ( ALLOCATED( pc_av ) )  THEN
          CALL wrd_write_string( 'pc_av' ) 
          WRITE ( 14 )  pc_av
       ENDIF

       IF ( ALLOCATED( pr_av ) )  THEN
          CALL wrd_write_string( 'pr_av' ) 
          WRITE ( 14 )  pr_av
       ENDIF

       CALL wrd_write_string( 'pt' ) 
       WRITE ( 14 )  pt

       IF ( ALLOCATED( pt_av ) )  THEN
          CALL wrd_write_string( 'pt_av' ) 
          WRITE ( 14 )  pt_av
       ENDIF

       IF ( humidity )  THEN

          CALL wrd_write_string( 'q' ) 
          WRITE ( 14 )  q

          IF ( ALLOCATED( q_av ) )  THEN
             CALL wrd_write_string( 'q_av' ) 
             WRITE ( 14 )  q_av
          ENDIF

          IF ( cloud_droplets )  THEN

             CALL wrd_write_string( 'ql' )
             WRITE ( 14 )  ql

             IF ( ALLOCATED( ql_av ) )  THEN
                CALL wrd_write_string( 'ql_av' )
                WRITE ( 14 )  ql_av
             ENDIF

          ENDIF

          IF ( ALLOCATED( qsws_av ) )  THEN
             CALL wrd_write_string( 'qsws_av' ) 
             WRITE ( 14 )  qsws_av
          ENDIF

       ENDIF

       IF ( passive_scalar )  THEN

          CALL wrd_write_string( 's' ) 
          WRITE ( 14 )  s

          IF ( ALLOCATED( s_av ) )  THEN
             CALL wrd_write_string( 's_av' ) 
             WRITE ( 14 )  s_av
          ENDIF

          IF ( ALLOCATED( ssws_av ) )  THEN
             CALL wrd_write_string( 'ssws_av' ) 
             WRITE ( 14 )  ssws_av
          ENDIF

       ENDIF
       
       IF ( ALLOCATED( ql_c_av ) )  THEN
          CALL wrd_write_string( 'ql_c_av' ) 
          WRITE ( 14 )  ql_c_av
       ENDIF

       IF ( ALLOCATED( ql_v_av ) )  THEN
          CALL wrd_write_string( 'ql_v_av' ) 
          WRITE ( 14 )  ql_v_av
       ENDIF

       IF ( ALLOCATED( ql_vp_av ) )  THEN
          CALL wrd_write_string( 'ql_vp_av' ) 
          WRITE ( 14 )  ql_vp_av
       ENDIF

       IF ( ALLOCATED( qv_av ) )  THEN
          CALL wrd_write_string( 'qv_av' ) 
          WRITE ( 14 )  qv_av
       ENDIF

       CALL wrd_write_string( 'random_iv' ) 
       WRITE ( 14 )  random_iv
       WRITE ( 14 )  random_iy

       IF ( ALLOCATED( seq_random_array ) )  THEN
          CALL wrd_write_string( 'seq_random_array' ) 
          WRITE ( 14 )  id_random_array
          WRITE ( 14 )  seq_random_array 
       ENDIF

       IF ( ALLOCATED( s_av ) )  THEN
          CALL wrd_write_string( 's_av' ) 
          WRITE ( 14 )  s_av
       ENDIF

       IF ( ALLOCATED( shf_av ) )  THEN
          CALL wrd_write_string( 'shf_av' ) 
          WRITE ( 14 )  shf_av
       ENDIF

       IF ( ALLOCATED( ts_av ) )  THEN
          CALL wrd_write_string( 'ts_av' ) 
          WRITE ( 14 ) ts_av
       ENDIF

       CALL wrd_write_string( 'u' ) 
       WRITE ( 14 )  u

       IF ( ALLOCATED( u_av ) )  THEN
          CALL wrd_write_string( 'u_av' ) 
          WRITE ( 14 )  u_av
       ENDIF

       IF ( ALLOCATED( u_m_l ) )  THEN
          CALL wrd_write_string( 'u_m_l' ) 
          WRITE ( 14 )  u_m_l
       ENDIF

       IF ( ALLOCATED( u_m_n ) )  THEN
          CALL wrd_write_string( 'u_m_n' ) 
          WRITE ( 14 )  u_m_n
       ENDIF

       IF ( ALLOCATED( u_m_r ) )  THEN
          CALL wrd_write_string( 'u_m_r' ) 
          WRITE ( 14 )  u_m_r
       ENDIF

       IF ( ALLOCATED( u_m_s ) )  THEN
          CALL wrd_write_string( 'u_m_s' ) 
          WRITE ( 14 )  u_m_s
       ENDIF

       IF ( ALLOCATED( us_av ) )  THEN
          CALL wrd_write_string( 'us_av' ) 
          WRITE ( 14 )  us_av
       ENDIF

       CALL wrd_write_string( 'v' ) 
       WRITE ( 14 )  v

       IF ( ALLOCATED( v_av ) )  THEN
          CALL wrd_write_string( 'v_av' ) 
          WRITE ( 14 )  v_av
       ENDIF

       IF ( ALLOCATED( v_m_l  ) )  THEN
          CALL wrd_write_string( 'v_m_l' ) 
          WRITE ( 14 )  v_m_l
       ENDIF

       IF ( ALLOCATED( v_m_n ) )  THEN
          CALL wrd_write_string( 'v_m_n' ) 
          WRITE ( 14 )  v_m_n
       ENDIF

       IF ( ALLOCATED( v_m_r ) )  THEN
          CALL wrd_write_string( 'v_m_r' ) 
          WRITE ( 14 )  v_m_r
       ENDIF

       IF ( ALLOCATED( v_m_s ) )  THEN
          CALL wrd_write_string( 'v_m_s' ) 
          WRITE ( 14 )  v_m_s
       ENDIF

       IF ( humidity )  THEN

          CALL wrd_write_string( 'vpt' ) 
          WRITE ( 14 )  vpt

          IF ( ALLOCATED( vpt_av ) )  THEN
             CALL wrd_write_string( 'vpt_av' ) 
             WRITE ( 14 ) vpt_av
          ENDIF

       ENDIF

       CALL wrd_write_string( 'w' ) 
       WRITE ( 14 )  w

       IF ( ALLOCATED( w_av ) )  THEN
          CALL wrd_write_string( 'w_av' ) 
          WRITE ( 14 )  w_av
       ENDIF

       IF ( ALLOCATED( w_m_l ) )  THEN
          CALL wrd_write_string( 'w_m_l' ) 
          WRITE ( 14 )  w_m_l
       ENDIF

       IF ( ALLOCATED( w_m_n ) )  THEN
          CALL wrd_write_string( 'w_m_n' ) 
          WRITE ( 14 )  w_m_n
       ENDIF

       IF ( ALLOCATED( w_m_r ) )  THEN
          CALL wrd_write_string( 'w_m_r' ) 
          WRITE ( 14 )  w_m_r
       ENDIF

       IF ( ALLOCATED( w_m_s ) )  THEN
          CALL wrd_write_string( 'w_m_s' ) 
          WRITE ( 14 )  w_m_s
       ENDIF

       IF ( ALLOCATED( z0_av ) )  THEN
          CALL wrd_write_string( 'z0_av' ) 
          WRITE ( 14 )  z0_av
       ENDIF

       IF ( ALLOCATED( z0h_av ) )  THEN
          CALL wrd_write_string( 'z0h_av' ) 
          WRITE ( 14 )  z0h_av
       ENDIF

       IF ( ALLOCATED( z0q_av ) )  THEN
          CALL wrd_write_string( 'z0q_av' ) 
          WRITE ( 14 )  z0q_av
       ENDIF

       IF ( land_surface  .OR.  urban_surface )  THEN

          IF ( ALLOCATED( ghf_av ) )  THEN
             CALL wrd_write_string( 'ghf_av' ) 
             WRITE ( 14 )  ghf_av
          ENDIF

          IF ( ALLOCATED( r_a_av ) )  THEN
             CALL wrd_write_string( 'r_a_av' ) 
             WRITE ( 14 )  r_a_av
          ENDIF

       ENDIF

       IF ( ALLOCATED( tsurf_av ) )  THEN
          CALL wrd_write_string( 'tsurf_av' ) 
          WRITE ( 14 )  tsurf_av
       ENDIF


       CALL surface_wrd_local

!
!--    Write restart data of other modules
       CALL module_interface_wrd_local

!
!--    Write end label
       CALL wrd_write_string( '*** end ***' )

    END SUBROUTINE wrd_local


 END MODULE write_restart_data_mod
