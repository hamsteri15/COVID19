!> @file time_integration.f90
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
! $Id: time_integration.f90 4466 2020-03-20 16:14:41Z suehring $
! Add advection fluxes to ACC copyin
!
! 4457 2020-03-11 14:20:43Z raasch
! use statement for exchange horiz added
!
! 4444 2020-03-05 15:59:50Z raasch
! bugfix: cpp-directives for serial mode added
!
! 4420 2020-02-24 14:13:56Z maronga
! Added output control for wind turbine model
!
! 4403 2020-02-12 13:08:46Z banzhafs
! Allowing both existing and on-demand emission read modes
!
! 4360 2020-01-07 11:25:50Z suehring
! Bugfix, hour_call_emis uninitialized at first call of time_integration
!
! 4346 2019-12-18 11:55:56Z motisi
! Introduction of wall_flags_total_0, which currently sets bits based on static
! topography information used in wall_flags_static_0
!
! 4329 2019-12-10 15:46:36Z motisi
! Renamed wall_flags_0 to wall_flags_static_0
!
! 4281 2019-10-29 15:15:39Z schwenkel
! Moved boundary conditions to module interface
!
! 4276 2019-10-28 16:03:29Z schwenkel
! Further modularization of lpm code components
!
! 4275 2019-10-28 15:34:55Z schwenkel
! Move call oft lpm to the end of intermediate timestep loop
!
! 4268 2019-10-17 11:29:38Z schwenkel
! Removing module specific boundary conditions an put them into their modules
!
! 4227 2019-09-10 18:04:34Z gronemeier
! implement new palm_date_time_mod
!
! 4226 2019-09-10 17:03:24Z suehring
! Changes in interface for the offline nesting
!
! 4182 2019-08-22 15:20:23Z scharf
! Corrected "Former revisions" section
!
! 4170 2019-08-19 17:12:31Z gronemeier
! copy diss, diss_p, tdiss_m to GPU
!
! 4144 2019-08-06 09:11:47Z raasch
! relational operators .EQ., .NE., etc. replaced by ==, /=, etc.
!
! 4126 2019-07-30 11:09:11Z gronemeier
! renamed routine to calculate uv exposure
!
! 4111 2019-07-22 18:16:57Z suehring
! advc_flags_1 / advc_flags_2 renamed to advc_flags_m / advc_flags_s
!
! 4069 2019-07-01 14:05:51Z Giersch
! Masked output running index mid has been introduced as a local variable to
! avoid runtime error (Loop variable has been modified) in time_integration
!
! 4064 2019-07-01 05:33:33Z gronemeier
! Moved call to radiation module out of intermediate time loop
!
! 4048 2019-06-21 21:00:21Z knoop
! Moved production_e_init call into turbulence_closure_mod
!
! 4047 2019-06-21 18:58:09Z knoop
! Added remainings of swap_timelevel upon its dissolution
!
! 4043 2019-06-18 16:59:00Z schwenkel
! Further LPM modularization
!
! 4039 2019-06-18 10:32:41Z suehring
! Rename subroutines in module for diagnostic quantities
!
! 4029 2019-06-14 14:04:35Z raasch
! exchange of ghost points and boundary conditions separated for chemical species and SALSA module,
! bugfix: decycling of chemistry species after nesting data transfer
!
! 4022 2019-06-12 11:52:39Z suehring
! Call synthetic turbulence generator at last RK3 substep right after boundary
! conditions are updated in offline nesting in order to assure that
! perturbations are always imposed
!
! 4017 2019-06-06 12:16:46Z schwenkel
! Mass (volume) flux correction included to ensure global mass conservation for child domains.
!
! 3994 2019-05-22 18:08:09Z suehring
! output of turbulence intensity added
!
! 3988 2019-05-22 11:32:37Z kanani
! Implement steerable output interval for virtual measurements
!
! 3968 2019-05-13 11:04:01Z suehring
! replace nspec_out with n_matched_vars
!
! 3929 2019-04-24 12:52:08Z banzhafs
! Reverse changes back from revision 3878: use chem_boundary_conds instead of
! chem_boundary_conds_decycle
!
!
! 3885 2019-04-11 11:29:34Z kanani
! Changes related to global restructuring of location messages and introduction
! of additional debug messages
!
! 3879 2019-04-08 20:25:23Z knoop
! Moved wtm_forces to module_interface_actions
!
! 3872 2019-04-08 15:03:06Z knoop
! Modifications made for salsa:
! - Call salsa_emission_update at each time step but do the checks within
!   salsa_emission_update (i.e. skip_time_do_salsa >= time_since_reference_point
!   and next_aero_emission_update <= time_since_reference_point ).
! - Renamed nbins --> nbins_aerosol, ncc_tot --> ncomponents_mass and
!   ngast --> ngases_salsa and loop indices b, c and sg to ib, ic and ig
! - Apply nesting for salsa variables
! - Removed cpu_log calls speciffic for salsa.
!
! 3833 2019-03-28 15:04:04Z forkel
! added USE chem_gasphase_mod, replaced nspec by nspec since fixed compounds are not integrated
!
! 3820 2019-03-27 11:53:41Z forkel
! renamed do_emiss to emissions_anthropogenic (ecc)
!
!
! 3774 2019-03-04 10:52:49Z moh.hefny
! rephrase if statement to avoid unallocated array in case of
! nesting_offline is false (crashing during debug mode)
!
! 3761 2019-02-25 15:31:42Z raasch $
! module section re-formatted and openacc required variables moved to separate section,
! re-formatting to 100 char line width
!
! 3745 2019-02-15 18:57:56Z suehring
! Call indoor model after first timestep
!
! 3744 2019-02-15 18:38:58Z suehring
! - Moved call of bio_calculate_thermal_index_maps from biometeorology module to
! time_integration to make sure averaged input is updated before calculating.
!
! 3739 2019-02-13 08:05:17Z dom_dwd_user
! Removed everything related to "time_bio_results" as this is never used.
!
! 3724 2019-02-06 16:28:23Z kanani
! Correct double-used log_point_s unit
!
! 3719 2019-02-06 13:10:18Z kanani
! - removed wind_turbine cpu measurement, since same time is measured inside
!   wtm_forces subroutine as special measures
! - moved the numerous vnest cpulog to special measures
! - extended radiation cpulog over entire radiation part,
!   moved radiation_interactions cpulog to special measures
! - moved some cpu_log calls to this routine for better overview
!
! 3705 2019-01-29 19:56:39Z suehring
! Data output for virtual measurements added
!
! 3704 2019-01-29 19:51:41Z suehring
! Rename subroutines for surface-data output
!
! 3647 2019-01-02 14:10:44Z kanani
! Bugfix: add time_since_reference_point to IF clause for data_output calls
! (otherwise skip_time_* values don't come into affect with dt_do* = 0.0).
! Clean up indoor_model and biometeorology model call.
!
! Revision 1.1  1997/08/11 06:19:04  raasch
! Initial revision
!
!
! Description:
! ------------
!> Integration in time of the model equations, statistical analysis and graphic
!> output
!------------------------------------------------------------------------------!
 SUBROUTINE time_integration


    USE advec_ws,                                                                                  &
        ONLY:  ws_statistics

    USE arrays_3d,                                                                                 &
        ONLY:  diss, diss_p, dzu, e_p, nc_p, nr_p, prho, pt, pt_p, pt_init, q, qc_p, qr_p, q_init, &
               q_p, ref_state, rho_ocean, sa_p, s_p, tend, u, u_p, v, vpt, v_p, w_p

#if defined( __parallel )  &&  ! defined( _OPENACC )
    USE arrays_3d,                                                                                 &
        ONLY:  e, nc, nr, qc, qr, s, w
#endif

    USE biometeorology_mod,                                                                        &
        ONLY:  bio_calculate_thermal_index_maps, thermal_comfort, bio_calculate_uv_exposure,       &
               uv_exposure

    USE bulk_cloud_model_mod,                                                                      &
        ONLY: bulk_cloud_model, calc_liquid_water_content, collision_turbulence,                   &
              microphysics_morrison, microphysics_seifert

    USE calc_mean_profile_mod,                                                                     &
        ONLY:  calc_mean_profile

    USE chem_emissions_mod,                                                                        &
        ONLY:  chem_emissions_setup, chem_emissions_update_on_demand

    USE chem_gasphase_mod,                                                                         &
        ONLY:  nvar

    USE chem_modules,                                                                              &
        ONLY:  bc_cs_t_val, chem_species, emissions_anthropogenic, emiss_read_legacy_mode,         &
               n_matched_vars

#if defined( __parallel )
    USE chem_modules,                                                                              &
        ONLY:  cs_name
#endif

    USE chemistry_model_mod,                                                                       &
        ONLY:  chem_boundary_conds

    USE control_parameters,                                                                        &
        ONLY:  advected_distance_x, advected_distance_y, air_chemistry, average_count_3d,          &
               averaging_interval, averaging_interval_pr, bc_lr_cyc, bc_ns_cyc, bc_pt_t_val,       &
               bc_q_t_val, biometeorology, call_psolver_at_all_substeps,  child_domain,            &
               constant_flux_layer, constant_heatflux, create_disturbances,        &
               dopr_n, constant_diffusion, coupling_mode, coupling_start_time,                     &
               current_timestep_number, disturbance_created, disturbance_energy_limit, dist_range, &
               do_sum, dt_3d, dt_averaging_input, dt_averaging_input_pr, dt_coupling,              &
               dt_data_output_av, dt_disturb, dt_do2d_xy, dt_do2d_xz, dt_do2d_yz, dt_do3d,         &
               dt_domask,dt_dopts, dt_dopr, dt_dopr_listing, dt_dots, dt_run_control,              &
               end_time, first_call_lpm, first_call_mas, galilei_transformation, humidity,         &
               indoor_model, intermediate_timestep_count, intermediate_timestep_count_max,         &
               land_surface, large_scale_forcing, loop_optimization, lsf_surf, lsf_vert, masks,    &
               multi_agent_system_end, multi_agent_system_start, nesting_offline, neutral,         &
               nr_timesteps_this_run, nudging, ocean_mode, passive_scalar, pt_reference,           &
               pt_slope_offset, random_heatflux, rans_tke_e, run_coupled, salsa,                   &
               simulated_time, simulated_time_chr, skip_time_do2d_xy, skip_time_do2d_xz,           &
               skip_time_do2d_yz, skip_time_do3d, skip_time_domask, skip_time_dopr,                &
               skip_time_data_output_av, sloping_surface, stop_dt, surface_output,                 &
               terminate_coupled, terminate_run, timestep_scheme, time_coupling, time_do2d_xy,     &
               time_do2d_xz, time_do2d_yz, time_do3d, time_domask, time_dopr, time_dopr_av,        &
               time_dopr_listing, time_dopts, time_dosp, time_dosp_av, time_dots, time_do_av,      &
               time_do_sla, time_disturb, time_run_control, time_since_reference_point,            &
               turbulent_inflow, turbulent_outflow, urban_surface,                                 &
               use_initial_profile_as_reference, use_single_reference_value, u_gtrans, v_gtrans,   &
               virtual_flight, virtual_measurement, ws_scheme_mom, ws_scheme_sca, timestep_count

#if defined( __parallel )
    USE control_parameters,                                                                        &
        ONLY:  rans_mode
#endif

    USE cpulog,                                                                                    &
        ONLY:  cpu_log, log_point, log_point_s

    USE diagnostic_output_quantities_mod,                                                          &
        ONLY:  doq_calculate,                                                                      &
               timestep_number_at_prev_calc

    USE exchange_horiz_mod,                                                                        &
        ONLY:  exchange_horiz

    USE flight_mod,                                                                                &
        ONLY:  flight_measurement

    USE indices,                                                                                   &
        ONLY:  nbgp, nx, nxl, nxlg, nxr, nxrg, nzb, nzt
!
!-- COVID-19 specific code
    USE indices,                                                                                   &
        ONLY:  ny, nyn, nyng, nys, nysg
!
!-- COVID-19 specific code ends
    USE indoor_model_mod,                                                                          &
        ONLY:  dt_indoor, im_main_heatcool, time_indoor

    USE interfaces

    USE kinds

    USE land_surface_model_mod,                                                                    &
        ONLY:  lsm_boundary_condition, lsm_energy_balance, lsm_soil_model, skip_time_do_lsm

    USE lagrangian_particle_model_mod,                                                             &
        ONLY:  lpm_data_output_ptseries

    USE lsf_nudging_mod,                                                                           &
        ONLY:  calc_tnudge, ls_forcing_surf, ls_forcing_vert, nudge_ref

    USE module_interface,                                                                          &
        ONLY:  module_interface_actions, module_interface_swap_timelevel,                          &
               module_interface_boundary_conditions

    USE multi_agent_system_mod,                                                                    &
        ONLY:  agents_active, multi_agent_system

    USE nesting_offl_mod,                                                                          &
        ONLY:  nesting_offl_bc,                                                                    &
               nesting_offl_geostrophic_wind,                                                      &
               nesting_offl_input,                                                                 &
               nesting_offl_interpolation_factor,                                                  &
               nesting_offl_mass_conservation

    USE netcdf_data_input_mod,                                                                     &
        ONLY:  chem_emis, chem_emis_att

    USE ocean_mod,                                                                                 &
        ONLY:  prho_reference

    USE palm_date_time_mod,                                                                        &
        ONLY:  get_date_time

    USE particle_attributes,                                                                       &
        ONLY:  particle_advection, particle_advection_start, use_sgs_for_particles, wang_kernel

    USE pegrid

#if defined( __parallel )
    USE pmc_interface,                                                                             &
        ONLY:  nested_run, nesting_mode, pmci_boundary_conds, pmci_datatrans, pmci_synchronize,    &
        pmci_ensure_nest_mass_conservation, pmci_ensure_nest_mass_conservation_vertical,           &
        pmci_set_swaplevel
#endif

    USE progress_bar,                                                                              &
        ONLY:  finish_progress_bar, output_progress_bar

    USE prognostic_equations_mod,                                                                  &
        ONLY:  prognostic_equations_cache, prognostic_equations_vector

    USE radiation_model_mod,                                                                       &
        ONLY: dt_radiation, force_radiation_call, radiation, radiation_control,                    &
              radiation_interaction, radiation_interactions, skip_time_do_radiation, time_radiation

    USE salsa_mod,                                                                                 &
        ONLY: aerosol_number, aerosol_mass, bc_am_t_val, bc_an_t_val, bc_gt_t_val,                 &
              nbins_aerosol, ncomponents_mass, ngases_salsa, salsa_boundary_conds,                 &
              salsa_emission_update, salsa_gas, salsa_gases_from_chem, skip_time_do_salsa

    USE spectra_mod,                                                                               &
        ONLY: average_count_sp, averaging_interval_sp, calc_spectra, dt_dosp, skip_time_dosp

    USE statistics,                                                                                &
        ONLY:  flow_statistics_called, hom, pr_palm, sums_ls_l


    USE surface_layer_fluxes_mod,                                                                  &
        ONLY:  surface_layer_fluxes

    USE surface_data_output_mod,                                                                   &
        ONLY:  average_count_surf, averaging_interval_surf, dt_dosurf, dt_dosurf_av,               &
               surface_data_output, surface_data_output_averaging, skip_time_dosurf,               &
               skip_time_dosurf_av, time_dosurf, time_dosurf_av

    USE surface_mod,                                                                               &
        ONLY:  surf_def_h, surf_lsm_h, surf_usm_h

    USE synthetic_turbulence_generator_mod,                                                        &
        ONLY:  dt_stg_call, dt_stg_adjust, parametrize_inflow_turbulence, stg_adjust, stg_main,    &
               time_stg_adjust, time_stg_call, use_syn_turb_gen

    USE turbulence_closure_mod,                                                                    &
        ONLY:  tcm_diffusivities

    USE urban_surface_mod,                                                                         &
        ONLY:  usm_boundary_condition, usm_material_heat_model, usm_material_model,                &
               usm_surface_energy_balance, usm_green_heat_model

    USE vertical_nesting_mod,                                                                      &
        ONLY:  vnested, vnest_init

#if defined( __parallel )
    USE vertical_nesting_mod,                                                                      &
        ONLY:  vnest_anterpolate, vnest_anterpolate_e, vnest_boundary_conds,                       &
               vnest_boundary_conds_khkm, vnest_deallocate, vnest_init_fine, vnest_start_time
#endif

    USE virtual_measurement_mod,                                                                   &
        ONLY:  dt_virtual_measurement,                                                             &
               time_virtual_measurement,                                                           &
               vm_data_output,                                                                     &
               vm_sampling,                                                                        &
               vm_time_start


    USE wind_turbine_model_mod,                                                                    &
        ONLY:  dt_data_output_wtm, time_wtm, wind_turbine, wtm_data_output

#if defined( _OPENACC )
    USE arrays_3d,                                                                                 &
        ONLY:  d, dd2zu, ddzu, ddzw,                                                               &
               diss_l_u,                                                                           &
               diss_l_v,                                                                           &
               diss_l_w,                                                                           &
               diss_s_u,                                                                           &
               diss_s_v,                                                                           &
               diss_s_w,                                                                           &
               drho_air, drho_air_zw, dzw, e,                                                      &
               flux_l_u,                                                                           &
               flux_l_v,                                                                           &
               flux_l_w,                                                                           &
               flux_s_u,                                                                           &
               flux_s_v,                                                                           &
               flux_s_w,                                                                           &
               heatflux_output_conversion,                                                         &
               kh, km, momentumflux_output_conversion, nc, nr, p, ptdf_x, ptdf_y, qc, qr, rdf,     &
               rdf_sc, rho_air, rho_air_zw, s, tdiss_m, te_m, tpt_m, tu_m, tv_m, tw_m, ug, u_init, &
               u_stokes_zu, vg, v_init, v_stokes_zu, w, zu

    USE control_parameters,                                                                        &
        ONLY:  tsc

    USE indices,                                                                                   &
        ONLY:  advc_flags_m, advc_flags_s, nyn, nyng, nys, nysg, nz, nzb_max, wall_flags_total_0

    USE statistics,                                                                                &
        ONLY:  rmask, statistic_regions, sums_l, sums_l_l, sums_us2_ws_l,                          &
               sums_wsus_ws_l, sums_vs2_ws_l, sums_wsvs_ws_l, sums_ws2_ws_l, sums_wspts_ws_l,      &
               sums_wsqs_ws_l, sums_wssas_ws_l, sums_wsqcs_ws_l, sums_wsqrs_ws_l, sums_wsncs_ws_l, &
               sums_wsnrs_ws_l, sums_wsss_ws_l, weight_substep, sums_salsa_ws_l

    USE surface_mod,                                                                               &
        ONLY:  bc_h, enter_surface_arrays, exit_surface_arrays
#endif


    IMPLICIT NONE

    CHARACTER (LEN=9) ::  time_to_string   !<

    INTEGER(iwp) ::  hour                !< hour of current time
    INTEGER(iwp) ::  hour_call_emis = -1 !< last hour where emission was called
    INTEGER(iwp) ::  ib                  !< index for aerosol size bins
    INTEGER(iwp) ::  ic                  !< index for aerosol mass bins
    INTEGER(iwp) ::  icc                 !< additional index for aerosol mass bins
    INTEGER(iwp) ::  ig                  !< index for salsa gases
    INTEGER(iwp) ::  lsp                 !<
    INTEGER(iwp) ::  mid                 !< masked output running index
#if defined( __parallel )
    INTEGER(iwp) ::  lsp_usr             !<
    INTEGER(iwp) ::  n                   !< loop counter for chemistry species
#endif

    REAL(wp) ::  dt_3d_old  !< temporary storage of timestep to be used for
                            !< steering of run control output interval
    REAL(wp) ::  time_since_reference_point_save  !< original value of
                                                  !< time_since_reference_point


!
!-- Copy data from arrays_3d
!$ACC DATA &
!$ACC COPY(d(nzb+1:nzt,nys:nyn,nxl:nxr)) &
!$ACC COPY(diss(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(e(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(u(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(v(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(w(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(kh(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(km(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(p(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(pt(nzb:nzt+1,nysg:nyng,nxlg:nxrg))

!$ACC DATA &
!$ACC COPYIN(diss_l_u(0:nz+1,nys:nyn,0), flux_l_u(0:nz+1,nys:nyn,0)) &
!$ACC COPYIN(diss_l_v(0:nz+1,nys:nyn,0), flux_l_v(0:nz+1,nys:nyn,0)) &
!$ACC COPYIN(diss_l_w(0:nz+1,nys:nyn,0), flux_l_w(0:nz+1,nys:nyn,0)) &
!$ACC COPYIN(diss_s_u(0:nz+1,0), flux_s_u(0:nz+1,0)) &
!$ACC COPYIN(diss_s_v(0:nz+1,0), flux_s_v(0:nz+1,0)) &
!$ACC COPYIN(diss_s_w(0:nz+1,0), flux_s_w(0:nz+1,0)) &
!$ACC COPY(diss_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(e_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(u_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(v_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(w_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(pt_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(tend(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(tdiss_m(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(te_m(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(tu_m(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(tv_m(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(tw_m(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(tpt_m(nzb:nzt+1,nysg:nyng,nxlg:nxrg))

!$ACC DATA &
!$ACC COPYIN(rho_air(nzb:nzt+1), drho_air(nzb:nzt+1)) &
!$ACC COPYIN(rho_air_zw(nzb:nzt+1), drho_air_zw(nzb:nzt+1)) &
!$ACC COPYIN(zu(nzb:nzt+1)) &
!$ACC COPYIN(dzu(1:nzt+1), dzw(1:nzt+1)) &
!$ACC COPYIN(ddzu(1:nzt+1), dd2zu(1:nzt)) &
!$ACC COPYIN(ddzw(1:nzt+1)) &
!$ACC COPYIN(heatflux_output_conversion(nzb:nzt+1)) &
!$ACC COPYIN(momentumflux_output_conversion(nzb:nzt+1)) &
!$ACC COPYIN(rdf(nzb+1:nzt), rdf_sc(nzb+1:nzt)) &
!$ACC COPYIN(ptdf_x(nxlg:nxrg), ptdf_y(nysg:nyng)) &
!$ACC COPYIN(ref_state(0:nz+1)) &
!$ACC COPYIN(u_init(0:nz+1), v_init(0:nz+1)) &
!$ACC COPYIN(u_stokes_zu(nzb:nzt+1), v_stokes_zu(nzb:nzt+1)) &
!$ACC COPYIN(pt_init(0:nz+1)) &
!$ACC COPYIN(ug(0:nz+1), vg(0:nz+1))

!
!-- Copy data from control_parameters
!$ACC DATA &
!$ACC COPYIN(tsc(1:5))

!
!-- Copy data from indices
!$ACC DATA &
!$ACC COPYIN(advc_flags_m(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPYIN(advc_flags_s(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPYIN(wall_flags_total_0(nzb:nzt+1,nysg:nyng,nxlg:nxrg))

!
!-- Copy data from surface_mod
!$ACC DATA &
!$ACC COPYIN(bc_h(0:1)) &
!$ACC COPYIN(bc_h(0)%i(1:bc_h(0)%ns)) &
!$ACC COPYIN(bc_h(0)%j(1:bc_h(0)%ns)) &
!$ACC COPYIN(bc_h(0)%k(1:bc_h(0)%ns)) &
!$ACC COPYIN(bc_h(1)%i(1:bc_h(1)%ns)) &
!$ACC COPYIN(bc_h(1)%j(1:bc_h(1)%ns)) &
!$ACC COPYIN(bc_h(1)%k(1:bc_h(1)%ns))

!
!-- Copy data from statistics
!$ACC DATA &
!$ACC COPYIN(hom(0:nz+1,1:2,1:4,0)) &
!$ACC COPYIN(rmask(nysg:nyng,nxlg:nxrg,0:statistic_regions)) &
!$ACC COPYIN(weight_substep(1:intermediate_timestep_count_max)) &
!$ACC COPY(sums_l(nzb:nzt+1,1:pr_palm,0)) &
!$ACC COPY(sums_l_l(nzb:nzt+1,0:statistic_regions,0)) &
!$ACC COPY(sums_us2_ws_l(nzb:nzt+1,0)) &
!$ACC COPY(sums_wsus_ws_l(nzb:nzt+1,0)) &
!$ACC COPY(sums_vs2_ws_l(nzb:nzt+1,0)) &
!$ACC COPY(sums_wsvs_ws_l(nzb:nzt+1,0)) &
!$ACC COPY(sums_ws2_ws_l(nzb:nzt+1,0)) &
!$ACC COPY(sums_wspts_ws_l(nzb:nzt+1,0)) &
!$ACC COPY(sums_wssas_ws_l(nzb:nzt+1,0)) &
!$ACC COPY(sums_wsqs_ws_l(nzb:nzt+1,0)) &
!$ACC COPY(sums_wsqcs_ws_l(nzb:nzt+1,0)) &
!$ACC COPY(sums_wsqrs_ws_l(nzb:nzt+1,0)) &
!$ACC COPY(sums_wsncs_ws_l(nzb:nzt+1,0)) &
!$ACC COPY(sums_wsnrs_ws_l(nzb:nzt+1,0)) &
!$ACC COPY(sums_wsss_ws_l(nzb:nzt+1,0)) &
!$ACC COPY(sums_salsa_ws_l(nzb:nzt+1,0))

#if defined( _OPENACC )
    CALL enter_surface_arrays
#endif

!
!-- At beginning determine the first time step
    CALL timestep

#if defined( __parallel )
!
!-- Synchronize the timestep in case of nested run.
    IF ( nested_run )  THEN
!
!--    Synchronization by unifying the time step.
!--    Global minimum of all time-steps is used for all.
       CALL pmci_synchronize
    ENDIF
#endif

!
!-- Determine and print out the run control quantities before the first time
!-- step of this run. For the initial run, some statistics (e.g. divergence)
!-- need to be determined first --> CALL flow_statistics at the beginning of
!-- run_control
    CALL run_control
!
!-- Data exchange between coupled models in case that a call has been omitted
!-- at the end of the previous run of a job chain.
    IF ( coupling_mode /= 'uncoupled'  .AND.  run_coupled  .AND. .NOT. vnested )  THEN
!
!--    In case of model termination initiated by the local model the coupler
!--    must not be called because this would again cause an MPI hang.
       DO WHILE ( time_coupling >= dt_coupling  .AND.  terminate_coupled == 0 )
          CALL surface_coupler
          time_coupling = time_coupling - dt_coupling
       ENDDO
       IF (time_coupling == 0.0_wp  .AND.  time_since_reference_point < dt_coupling )  THEN
          time_coupling = time_since_reference_point
       ENDIF
    ENDIF

    CALL location_message( 'atmosphere (and/or ocean) time-stepping', 'start' )

!
!-- Start of the time loop
    DO  WHILE ( simulated_time < end_time  .AND.  .NOT. stop_dt  .AND. .NOT. terminate_run )

       CALL cpu_log( log_point_s(10), 'timesteps', 'start' )

#if defined( __parallel )
!
!--    Vertical nesting: initialize fine grid
       IF ( vnested ) THEN
          IF ( .NOT. vnest_init  .AND.  simulated_time >= vnest_start_time )  THEN
             CALL cpu_log( log_point_s(22), 'vnest_init', 'start' )
             CALL vnest_init_fine
             vnest_init = .TRUE.
             CALL cpu_log( log_point_s(22), 'vnest_init', 'stop' )
          ENDIF
       ENDIF
#endif

!
!--    Determine ug, vg and w_subs in dependence on data from external file
!--    LSF_DATA
       IF ( large_scale_forcing .AND. lsf_vert )  THEN
           CALL ls_forcing_vert ( simulated_time )
           sums_ls_l = 0.0_wp
       ENDIF

!
!--    Set pt_init and q_init to the current profiles taken from
!--    NUDGING_DATA
       IF ( nudging )  THEN
           CALL nudge_ref ( simulated_time )
!
!--        Store temperature gradient at the top boundary for possible Neumann
!--        boundary condition
           bc_pt_t_val = ( pt_init(nzt+1) - pt_init(nzt) ) / dzu(nzt+1)
           bc_q_t_val  = ( q_init(nzt+1) - q_init(nzt) ) / dzu(nzt+1)
           IF ( air_chemistry )  THEN
              DO  lsp = 1, nvar
                 bc_cs_t_val = (  chem_species(lsp)%conc_pr_init(nzt+1)                            &
                                - chem_species(lsp)%conc_pr_init(nzt) )                            &
                               / dzu(nzt+1)
              ENDDO
           ENDIF
           IF ( salsa  .AND.  time_since_reference_point >= skip_time_do_salsa )  THEN
              DO  ib = 1, nbins_aerosol
                 bc_an_t_val = ( aerosol_number(ib)%init(nzt+1) - aerosol_number(ib)%init(nzt) ) / &
                               dzu(nzt+1)
                 DO  ic = 1, ncomponents_mass
                    icc = ( ic - 1 ) * nbins_aerosol + ib
                    bc_am_t_val = ( aerosol_mass(icc)%init(nzt+1) - aerosol_mass(icc)%init(nzt) ) /&
                                  dzu(nzt+1)
                 ENDDO
              ENDDO
              IF ( .NOT. salsa_gases_from_chem )  THEN
                 DO  ig = 1, ngases_salsa
                    bc_gt_t_val = ( salsa_gas(ig)%init(nzt+1) - salsa_gas(ig)%init(nzt) ) /        &
                                  dzu(nzt+1)
                 ENDDO
              ENDIF
           ENDIF
       ENDIF
!
!--    Input of boundary data.
       IF ( nesting_offline )  CALL nesting_offl_input
!
!--    Execute all other module actions routines
       CALL module_interface_actions( 'before_timestep' )

!--    Start of intermediate step loop
       intermediate_timestep_count = 0
       DO  WHILE ( intermediate_timestep_count < intermediate_timestep_count_max )

          intermediate_timestep_count = intermediate_timestep_count + 1

!
!--       Set the steering factors for the prognostic equations which depend
!--       on the timestep scheme
          CALL timestep_scheme_steering

!
!--       Calculate those variables needed in the tendency terms which need
!--       global communication
          IF ( .NOT. use_single_reference_value  .AND.  .NOT. use_initial_profile_as_reference )   &
          THEN
!
!--          Horizontally averaged profiles to be used as reference state in
!--          buoyancy terms (WARNING: only the respective last call of
!--          calc_mean_profile defines the reference state!)
             IF ( .NOT. neutral )  THEN
                CALL calc_mean_profile( pt, 4 )
                ref_state(:)  = hom(:,1,4,0) ! this is used in the buoyancy term
             ENDIF
             IF ( ocean_mode )  THEN
                CALL calc_mean_profile( rho_ocean, 64 )
                ref_state(:)  = hom(:,1,64,0)
             ENDIF
             IF ( humidity )  THEN
                CALL calc_mean_profile( vpt, 44 )
                ref_state(:)  = hom(:,1,44,0)
             ENDIF
!
!--          Assure that ref_state does not become zero at any level
!--          ( might be the case if a vertical level is completely occupied
!--            with topography ).
             ref_state = MERGE( MAXVAL(ref_state), ref_state, ref_state == 0.0_wp )
          ENDIF

          IF ( ( ws_scheme_mom .OR. ws_scheme_sca )  .AND.  intermediate_timestep_count == 1 )     &
          THEN
             CALL ws_statistics
          ENDIF
!
!--       In case of nudging calculate current nudging time scale and horizontal
!--       means of u, v, pt and q
          IF ( nudging )  THEN
             CALL calc_tnudge( simulated_time )
             CALL calc_mean_profile( u, 1 )
             CALL calc_mean_profile( v, 2 )
             CALL calc_mean_profile( pt, 4 )
             CALL calc_mean_profile( q, 41 )
          ENDIF
!
!--       Execute all other module actions routunes
          CALL module_interface_actions( 'before_prognostic_equations' )
!
!--       Solve the prognostic equations. A fast cache optimized version with
!--       only one single loop is used in case of Piascek-Williams advection
!--       scheme. NEC vector machines use a different version, because
!--       in the other versions a good vectorization is prohibited due to
!--       inlining problems.
          IF ( loop_optimization == 'cache' )  THEN
             CALL prognostic_equations_cache
          ELSEIF ( loop_optimization == 'vector' )  THEN
             CALL prognostic_equations_vector
          ENDIF
!
!--       Movement of agents in multi agent system
          IF ( agents_active  .AND.  time_since_reference_point >= multi_agent_system_start  .AND. &
               time_since_reference_point <= multi_agent_system_end  .AND.                         &
               intermediate_timestep_count == 1 )                                                  &
          THEN
             CALL multi_agent_system
             first_call_mas = .FALSE.
          ENDIF

!
!--       Exchange of ghost points (lateral boundary conditions)
          CALL cpu_log( log_point(26), 'exchange-horiz-progn', 'start' )

          CALL exchange_horiz( u_p, nbgp )
          CALL exchange_horiz( v_p, nbgp )
          CALL exchange_horiz( w_p, nbgp )
          CALL exchange_horiz( pt_p, nbgp )
          IF ( .NOT. constant_diffusion )  CALL exchange_horiz( e_p, nbgp )
          IF ( rans_tke_e  .OR.  wang_kernel  .OR.  collision_turbulence                           &
               .OR.  use_sgs_for_particles )  THEN
             IF ( rans_tke_e )  THEN
                CALL exchange_horiz( diss_p, nbgp )
             ELSE
                CALL exchange_horiz( diss, nbgp )
             ENDIF
          ENDIF
          IF ( ocean_mode )  THEN
             CALL exchange_horiz( sa_p, nbgp )
             CALL exchange_horiz( rho_ocean, nbgp )
             CALL exchange_horiz( prho, nbgp )
          ENDIF
          IF ( humidity )  THEN
             CALL exchange_horiz( q_p, nbgp )
             IF ( bulk_cloud_model .AND. microphysics_morrison )  THEN
                CALL exchange_horiz( qc_p, nbgp )
                CALL exchange_horiz( nc_p, nbgp )
             ENDIF
             IF ( bulk_cloud_model .AND. microphysics_seifert )  THEN
                CALL exchange_horiz( qr_p, nbgp )
                CALL exchange_horiz( nr_p, nbgp )
             ENDIF
          ENDIF
          IF ( passive_scalar )  CALL exchange_horiz( s_p, nbgp )
!
!--       COVID-19 specific code
!--       Set Dirichlet-0 conditions for s_p on left, right, south and north boundaries.
          IF ( passive_scalar )  THEN
             IF ( nxl == 0 )  THEN
                s_p(:,:,nxlg:-1) = 0.0_wp
             ENDIF
             IF ( nxr == nx )  THEN
                s_p(:,:,nx+1:nxrg) = 0.0_wp
             ENDIF
             IF ( nys == 0 )  THEN
                s_p(:,nysg:-1,:) = 0.0_wp
             ENDIF
             IF ( nyn == ny )  THEN
                s_p(:,ny+1:nyng,:) = 0.0_wp
             ENDIF
          ENDIF
!
!--       COVID-19 specific code ends
          IF ( air_chemistry )  THEN
             DO  lsp = 1, nvar
                CALL exchange_horiz( chem_species(lsp)%conc_p, nbgp )
             ENDDO
          ENDIF

          IF ( salsa  .AND.  time_since_reference_point >= skip_time_do_salsa )  THEN
             DO  ib = 1, nbins_aerosol
                CALL exchange_horiz( aerosol_number(ib)%conc_p, nbgp )
                DO  ic = 1, ncomponents_mass
                   icc = ( ic - 1 ) * nbins_aerosol + ib
                   CALL exchange_horiz( aerosol_mass(icc)%conc_p, nbgp )
                ENDDO
             ENDDO
             IF ( .NOT. salsa_gases_from_chem )  THEN
                DO  ig = 1, ngases_salsa
                   CALL exchange_horiz( salsa_gas(ig)%conc_p, nbgp )
                ENDDO
             ENDIF
          ENDIF
          CALL cpu_log( log_point(26), 'exchange-horiz-progn', 'stop' )

!
!--       Boundary conditions for the prognostic quantities (except of the
!--       velocities at the outflow in case of a non-cyclic lateral wall) and
!--       boundary conditions for module-specific variables
          CALL module_interface_boundary_conditions
!
!--       Incrementing timestep counter
          timestep_count = timestep_count + 1

          CALL cpu_log( log_point(28), 'swap_timelevel', 'start' )
!
!--       Set the swap level for all modules
          CALL module_interface_swap_timelevel( MOD( timestep_count, 2) )

#if defined( __parallel )
!
!--       Set the swap level for steering the pmc data transfer
          IF ( nested_run )  CALL pmci_set_swaplevel( MOD( timestep_count, 2) + 1 )  !> @todo: why the +1 ?
#endif

          CALL cpu_log( log_point(28), 'swap_timelevel', 'stop' )

#if defined( __parallel )
!
!--       Vertical nesting: Interpolate fine grid data to the coarse grid
          IF ( vnest_init ) THEN
             CALL cpu_log( log_point_s(37), 'vnest_anterpolate', 'start' )
             CALL vnest_anterpolate
             CALL cpu_log( log_point_s(37), 'vnest_anterpolate', 'stop' )
          ENDIF

          IF ( nested_run )  THEN

             CALL cpu_log( log_point(60), 'nesting', 'start' )
!
!--          Domain nesting. The data transfer subroutines pmci_parent_datatrans
!--          and pmci_child_datatrans are called inside the wrapper
!--          subroutine pmci_datatrans according to the control parameters
!--          nesting_mode and nesting_datatransfer_mode.
!--          TO_DO: why is nesting_mode given as a parameter here?
             CALL pmci_datatrans( nesting_mode )

             IF ( TRIM( nesting_mode ) == 'two-way' .OR.  nesting_mode == 'vertical' )  THEN

                CALL cpu_log( log_point_s(92), 'exchange-horiz-nest', 'start' )
!
!--             Exchange_horiz is needed for all parent-domains after the
!--             anterpolation
                CALL exchange_horiz( u, nbgp )
                CALL exchange_horiz( v, nbgp )
                CALL exchange_horiz( w, nbgp )
                IF ( .NOT. neutral )  CALL exchange_horiz( pt, nbgp )

                IF ( humidity )  THEN

                   CALL exchange_horiz( q, nbgp )

                   IF ( bulk_cloud_model  .AND.  microphysics_morrison )  THEN
                       CALL exchange_horiz( qc, nbgp )
                       CALL exchange_horiz( nc, nbgp )
                   ENDIF
                   IF ( bulk_cloud_model  .AND.  microphysics_seifert )  THEN
                       CALL exchange_horiz( qr, nbgp )
                       CALL exchange_horiz( nr, nbgp )
                   ENDIF

                ENDIF

                IF ( passive_scalar )  CALL exchange_horiz( s, nbgp )
!
!--             COVID-19 specific code
!--             Set Dirichlet-0 conditions for s_p on left, right, south and north boundaries.
!--             Nesting is actually not planned to be used for for the COVID-19 campaign
!--             but this is inserted also here just to be consistent.
                IF ( passive_scalar )  THEN
                   IF ( nxl == 0 )  THEN
                      s(:,:,nxlg:-1) = 0.0_wp
                   ENDIF
                   IF ( nxr == nx )  THEN
                      s(:,:,nx+1:nxrg) = 0.0_wp
                   ENDIF
                   IF ( nys == 0 )  THEN
                      s(:,nysg:-1,:) = 0.0_wp
                   ENDIF
                   IF ( nyn == ny )  THEN
                      s(:,ny+1:nyng,:) = 0.0_wp
                   ENDIF
                ENDIF
!
!--             COVID-19 specific code ends
                IF ( .NOT. constant_diffusion )  CALL exchange_horiz( e, nbgp )

                IF ( .NOT. constant_diffusion  .AND.  rans_mode  .AND.  rans_tke_e )  THEN
                   CALL exchange_horiz( diss, nbgp )
                ENDIF

                IF ( air_chemistry )  THEN
                   DO  n = 1, nvar
                      CALL exchange_horiz( chem_species(n)%conc, nbgp )
                   ENDDO
                ENDIF

                IF ( salsa  .AND. time_since_reference_point >= skip_time_do_salsa )  THEN
                   DO  ib = 1, nbins_aerosol
                      CALL exchange_horiz( aerosol_number(ib)%conc, nbgp )
                      DO  ic = 1, ncomponents_mass
                         icc = ( ic - 1 ) * nbins_aerosol + ib
                         CALL exchange_horiz( aerosol_mass(icc)%conc, nbgp )
                      ENDDO
                   ENDDO
                   IF ( .NOT. salsa_gases_from_chem )  THEN
                      DO  ig = 1, ngases_salsa
                         CALL exchange_horiz( salsa_gas(ig)%conc, nbgp )
                      ENDDO
                   ENDIF
                ENDIF
                CALL cpu_log( log_point_s(92), 'exchange-horiz-nest', 'stop' )

             ENDIF

!
!--          Set boundary conditions again after interpolation and anterpolation.
             CALL pmci_boundary_conds

!
!--          Set chemistry boundary conditions (decycling)
             IF ( air_chemistry )  THEN
                DO  lsp = 1, nvar
                   lsp_usr = 1
                   DO WHILE ( TRIM( cs_name( lsp_usr ) ) /= 'novalue' )
                      IF ( TRIM( chem_species(lsp)%name ) == TRIM( cs_name(lsp_usr) ) )  THEN
                         CALL chem_boundary_conds( chem_species(lsp)%conc,                         &
                                                   chem_species(lsp)%conc_pr_init )
                      ENDIF
                      lsp_usr = lsp_usr + 1
                   ENDDO
                ENDDO
             ENDIF

!
!--          Set SALSA boundary conditions (decycling)
             IF ( salsa  .AND. time_since_reference_point >= skip_time_do_salsa )  THEN
                DO  ib = 1, nbins_aerosol
                   CALL salsa_boundary_conds( aerosol_number(ib)%conc, aerosol_number(ib)%init )
                   DO  ic = 1, ncomponents_mass
                      icc = ( ic - 1 ) * nbins_aerosol + ib
                      CALL salsa_boundary_conds( aerosol_mass(icc)%conc, aerosol_mass(icc)%init )
                   ENDDO
                ENDDO
                IF ( .NOT. salsa_gases_from_chem )  THEN
                   DO  ig = 1, ngases_salsa
                      CALL salsa_boundary_conds( salsa_gas(ig)%conc, salsa_gas(ig)%init )
                   ENDDO
                ENDIF
             ENDIF

             CALL cpu_log( log_point(60), 'nesting', 'stop' )

          ENDIF
#endif

!
!--       Temperature offset must be imposed at cyclic boundaries in x-direction
!--       when a sloping surface is used
          IF ( sloping_surface )  THEN
             IF ( nxl ==  0 )  pt(:,:,nxlg:nxl-1) = pt(:,:,nxlg:nxl-1) - pt_slope_offset
             IF ( nxr == nx )  pt(:,:,nxr+1:nxrg) = pt(:,:,nxr+1:nxrg) + pt_slope_offset
          ENDIF

!
!--       Impose a turbulent inflow using the recycling method
          IF ( turbulent_inflow )  CALL inflow_turbulence

!
!--       Set values at outflow boundary using the special outflow condition
          IF ( turbulent_outflow )  CALL outflow_turbulence

!
!--       Impose a random perturbation on the horizontal velocity field
          IF ( create_disturbances  .AND.  ( call_psolver_at_all_substeps  .AND.                   &
               intermediate_timestep_count == intermediate_timestep_count_max )                    &
               .OR. ( .NOT. call_psolver_at_all_substeps  .AND.  intermediate_timestep_count == 1 ) ) &
          THEN
             time_disturb = time_disturb + dt_3d
             IF ( time_disturb >= dt_disturb )  THEN
                IF ( disturbance_energy_limit /= 0.0_wp  .AND.                                     &
                     hom(nzb+5,1,pr_palm,0) < disturbance_energy_limit )  THEN
                   CALL disturb_field( 'u', tend, u )
                   CALL disturb_field( 'v', tend, v )
                ELSEIF ( ( .NOT. bc_lr_cyc  .OR.  .NOT. bc_ns_cyc )                                &
                         .AND. .NOT. child_domain  .AND.  .NOT.  nesting_offline )                 &
                THEN
!
!--                Runs with a non-cyclic lateral wall need perturbations
!--                near the inflow throughout the whole simulation
                   dist_range = 1
                   CALL disturb_field( 'u', tend, u )
                   CALL disturb_field( 'v', tend, v )
                   dist_range = 0
                ENDIF
                time_disturb = time_disturb - dt_disturb
             ENDIF
          ENDIF

!
!--       Map forcing data derived from larger scale model onto domain
!--       boundaries. Further, update geostrophic wind components.
          IF ( nesting_offline  .AND.  intermediate_timestep_count ==                              &
                                       intermediate_timestep_count_max  )  THEN
!--          Determine interpolation factor before boundary conditions and geostrophic wind
!--          is updated.
             CALL nesting_offl_interpolation_factor
             CALL nesting_offl_bc
             CALL nesting_offl_geostrophic_wind
          ENDIF
!
!--       Impose a turbulent inflow using synthetic generated turbulence.
          IF ( use_syn_turb_gen  .AND.                                                             &
               intermediate_timestep_count == intermediate_timestep_count_max )  THEN
             CALL cpu_log( log_point(57), 'synthetic_turbulence_gen', 'start' )
             CALL stg_main
             CALL cpu_log( log_point(57), 'synthetic_turbulence_gen', 'stop' )
          ENDIF
!
!--       Ensure mass conservation. This need to be done after imposing
!--       synthetic turbulence and top boundary condition for pressure is set to
!--       Neumann conditions.
!--       Is this also required in case of Dirichlet?
          IF ( nesting_offline )  CALL nesting_offl_mass_conservation
!
!--       Reduce the velocity divergence via the equation for perturbation
!--       pressure.
          IF ( intermediate_timestep_count == 1  .OR. &
                call_psolver_at_all_substeps )  THEN

             IF (  vnest_init ) THEN
#if defined( __parallel )
!
!--             Compute pressure in the CG, interpolate top boundary conditions
!--             to the FG and then compute pressure in the FG
                IF ( coupling_mode == 'vnested_crse' )  CALL pres

                CALL cpu_log( log_point_s(30), 'vnest_bc', 'start' )
                CALL vnest_boundary_conds
                CALL cpu_log( log_point_s(30), 'vnest_bc', 'stop' )

                IF ( coupling_mode == 'vnested_fine' )  CALL pres

!--             Anterpolate TKE, satisfy Germano Identity
                CALL cpu_log( log_point_s(28), 'vnest_anter_e', 'start' )
                CALL vnest_anterpolate_e
                CALL cpu_log( log_point_s(28), 'vnest_anter_e', 'stop' )
#else
                CONTINUE
#endif

             ELSE
#if defined( __parallel )
!
!--             Mass (volume) flux correction to ensure global mass conservation for child domains.
                IF ( child_domain )  THEN
                   IF ( nesting_mode == 'vertical' )  THEN
                      CALL pmci_ensure_nest_mass_conservation_vertical
                   ELSE
                      CALL pmci_ensure_nest_mass_conservation
                   ENDIF
                ENDIF
#endif
                CALL pres

             ENDIF

          ENDIF
!
!--       Particle transport/physics with the Lagrangian particle model
!--       (only once during intermediate steps, because it uses an Euler-step)
!--       ### particle model should be moved before prognostic_equations, in order
!--       to regard droplet interactions directly

          CALL module_interface_actions( 'after_pressure_solver' )
!
!--       Interaction of droplets with temperature and mixing ratio.
!--       Droplet condensation and evaporation is calculated within
!--       advec_particles.
!
!--       If required, compute liquid water content
          IF ( bulk_cloud_model )  THEN
             CALL calc_liquid_water_content
          ENDIF
!
!--       If required, compute virtual potential temperature
          IF ( humidity )  THEN
             CALL compute_vpt
          ENDIF

!
!--       Compute the diffusion quantities
          IF ( .NOT. constant_diffusion )  THEN

!
!--          Determine surface fluxes shf and qsws and surface values
!--          pt_surface and q_surface in dependence on data from external
!--          file LSF_DATA respectively
             IF ( ( large_scale_forcing .AND. lsf_surf ) .AND.                                     &
                 intermediate_timestep_count == intermediate_timestep_count_max )                  &
             THEN
                CALL ls_forcing_surf( simulated_time )
             ENDIF

!
!--          First the vertical (and horizontal) fluxes in the surface
!--          (constant flux) layer are computed
             IF ( constant_flux_layer )  THEN
                CALL cpu_log( log_point(19), 'surface_layer_fluxes', 'start' )
                CALL surface_layer_fluxes
                CALL cpu_log( log_point(19), 'surface_layer_fluxes', 'stop' )
             ENDIF
!
!--          If required, solve the energy balance for the surface and run soil
!--          model. Call for horizontal as well as vertical surfaces
             IF ( land_surface .AND. time_since_reference_point >= skip_time_do_lsm)  THEN

                CALL cpu_log( log_point(54), 'land_surface', 'start' )
!
!--             Call for horizontal upward-facing surfaces
                CALL lsm_energy_balance( .TRUE., -1 )
                CALL lsm_soil_model( .TRUE., -1, .TRUE. )
!
!--             Call for northward-facing surfaces
                CALL lsm_energy_balance( .FALSE., 0 )
                CALL lsm_soil_model( .FALSE., 0, .TRUE. )
!
!--             Call for southward-facing surfaces
                CALL lsm_energy_balance( .FALSE., 1 )
                CALL lsm_soil_model( .FALSE., 1, .TRUE. )
!
!--             Call for eastward-facing surfaces
                CALL lsm_energy_balance( .FALSE., 2 )
                CALL lsm_soil_model( .FALSE., 2, .TRUE. )
!
!--             Call for westward-facing surfaces
                CALL lsm_energy_balance( .FALSE., 3 )
                CALL lsm_soil_model( .FALSE., 3, .TRUE. )

!
!--             At the end, set boundary conditons for potential temperature
!--             and humidity after running the land-surface model. This
!--             might be important for the nesting, where arrays are transfered.
                CALL lsm_boundary_condition


                CALL cpu_log( log_point(54), 'land_surface', 'stop' )
             ENDIF
!
!--          If required, solve the energy balance for urban surfaces and run
!--          the material heat model
             IF (urban_surface) THEN
                CALL cpu_log( log_point(74), 'urban_surface', 'start' )

                CALL usm_surface_energy_balance( .FALSE. )
                IF ( usm_material_model )  THEN
                   CALL usm_green_heat_model
                   CALL usm_material_heat_model ( .FALSE. )
                ENDIF

!
!--             At the end, set boundary conditons for potential temperature
!--             and humidity after running the urban-surface model. This
!--             might be important for the nesting, where arrays are transfered.
                CALL usm_boundary_condition

                CALL cpu_log( log_point(74), 'urban_surface', 'stop' )
             ENDIF
!
!--          Compute the diffusion coefficients
             CALL cpu_log( log_point(17), 'diffusivities', 'start' )
             IF ( .NOT. humidity ) THEN
                IF ( ocean_mode )  THEN
                   CALL tcm_diffusivities( prho, prho_reference )
                ELSE
                   CALL tcm_diffusivities( pt, pt_reference )
                ENDIF
             ELSE
                CALL tcm_diffusivities( vpt, pt_reference )
             ENDIF
             CALL cpu_log( log_point(17), 'diffusivities', 'stop' )

#if defined( __parallel )
!
!--          Vertical nesting: set fine grid eddy viscosity top boundary condition
             IF ( vnest_init )  CALL vnest_boundary_conds_khkm
#endif

          ENDIF

       ENDDO   ! Intermediate step loop

!
!--    Will be used at some point by flow_statistics.
       !$ACC UPDATE &
       !$ACC HOST(sums_l_l(nzb:nzt+1,0:statistic_regions,0)) &
       !$ACC HOST(sums_us2_ws_l(nzb:nzt+1,0)) &
       !$ACC HOST(sums_wsus_ws_l(nzb:nzt+1,0)) &
       !$ACC HOST(sums_vs2_ws_l(nzb:nzt+1,0)) &
       !$ACC HOST(sums_wsvs_ws_l(nzb:nzt+1,0)) &
       !$ACC HOST(sums_ws2_ws_l(nzb:nzt+1,0)) &
       !$ACC HOST(sums_wspts_ws_l(nzb:nzt+1,0)) &
       !$ACC HOST(sums_wssas_ws_l(nzb:nzt+1,0)) &
       !$ACC HOST(sums_wsqs_ws_l(nzb:nzt+1,0)) &
       !$ACC HOST(sums_wsqcs_ws_l(nzb:nzt+1,0)) &
       !$ACC HOST(sums_wsqrs_ws_l(nzb:nzt+1,0)) &
       !$ACC HOST(sums_wsncs_ws_l(nzb:nzt+1,0)) &
       !$ACC HOST(sums_wsnrs_ws_l(nzb:nzt+1,0)) &
       !$ACC HOST(sums_wsss_ws_l(nzb:nzt+1,0)) &
       !$ACC HOST(sums_salsa_ws_l(nzb:nzt+1,0))

!
!--    If required, calculate radiative fluxes and heating rates
       IF ( radiation  .AND.  time_since_reference_point > skip_time_do_radiation )  THEN

          time_radiation = time_radiation + dt_3d

          IF ( time_radiation >= dt_radiation  .OR.  force_radiation_call )  THEN

             CALL cpu_log( log_point(50), 'radiation', 'start' )

             IF ( .NOT. force_radiation_call )  THEN
                time_radiation = time_radiation - dt_radiation
             ENDIF

!
!--          Adjust the current time to the time step of the radiation model.
!--          Needed since radiation is pre-calculated and stored only on apparent
!--          solar positions
             time_since_reference_point_save = time_since_reference_point
             time_since_reference_point = REAL( FLOOR( time_since_reference_point /             &
                                                       dt_radiation ), wp ) * dt_radiation

             CALL radiation_control

             IF ( ( urban_surface  .OR.  land_surface )  .AND.  radiation_interactions )  THEN
                CALL cpu_log( log_point_s(46), 'radiation_interaction', 'start' )
                CALL radiation_interaction
                CALL cpu_log( log_point_s(46), 'radiation_interaction', 'stop' )
             ENDIF

!
!--          Return the current time to its original value
             time_since_reference_point = time_since_reference_point_save

             CALL cpu_log( log_point(50), 'radiation', 'stop' )

          ENDIF
       ENDIF


!
!-- 20200203 (ECC)
!-- allows for emission update mode in legacy mode as well as on-demand mode
!-- note that under on-demand mode emission update is no longer restricted to
!-- an hourly frequency, but whenever the simulation time corresponds to an
!-- inrement in emission timestamp value

!
!-- If required, consider chemical emissions

       IF  ( air_chemistry .AND. emissions_anthropogenic )  THEN

          IF  ( emiss_read_legacy_mode )  THEN
!
!-- get hourly index and updates emission data when the hour is passed

             CALL get_date_time( time_since_reference_point, hour=hour )

             IF  ( hour_call_emis /= hour )   THEN

                CALL chem_emissions_setup( chem_emis_att, chem_emis, n_matched_vars )
                hour_call_emis = hour

             ENDIF

          ELSE

             CALL chem_emissions_update_on_demand

          ENDIF

       ENDIF


!
!--    If required, consider aerosol emissions for the salsa model
       IF ( salsa )  THEN
!
!--       Call emission routine to update emissions if needed
          CALL salsa_emission_update

       ENDIF
!
!--    If required, calculate indoor temperature, waste heat, heat flux
!--    through wall, etc.
!--    dt_indoor steers the frequency of the indoor model calculations.
!--    Note, at first timestep indoor model is called, in order to provide
!--    a waste heat flux.
       IF ( indoor_model )  THEN

          time_indoor = time_indoor + dt_3d

          IF ( time_indoor >= dt_indoor  .OR.  current_timestep_number == 0 )  THEN

             time_indoor = time_indoor - dt_indoor

             CALL cpu_log( log_point(76), 'indoor_model', 'start' )
             CALL im_main_heatcool
             CALL cpu_log( log_point(76), 'indoor_model', 'stop' )

          ENDIF
       ENDIF
!
!--    Increase simulation time and output times
       nr_timesteps_this_run      = nr_timesteps_this_run + 1
       current_timestep_number    = current_timestep_number + 1
       simulated_time             = simulated_time   + dt_3d
       time_since_reference_point = simulated_time - coupling_start_time
       simulated_time_chr         = time_to_string( time_since_reference_point )


       IF ( time_since_reference_point >= skip_time_data_output_av )  THEN
          time_do_av         = time_do_av       + dt_3d
       ENDIF
       IF ( time_since_reference_point >= skip_time_do2d_xy )  THEN
          time_do2d_xy       = time_do2d_xy     + dt_3d
       ENDIF
       IF ( time_since_reference_point >= skip_time_do2d_xz )  THEN
          time_do2d_xz       = time_do2d_xz     + dt_3d
       ENDIF
       IF ( time_since_reference_point >= skip_time_do2d_yz )  THEN
          time_do2d_yz       = time_do2d_yz     + dt_3d
       ENDIF
       IF ( time_since_reference_point >= skip_time_do3d    )  THEN
          time_do3d          = time_do3d        + dt_3d
       ENDIF
       DO  mid = 1, masks
          IF ( time_since_reference_point >= skip_time_domask(mid) )  THEN
             time_domask(mid)= time_domask(mid) + dt_3d
          ENDIF
       ENDDO
       IF ( time_since_reference_point >= skip_time_dosp )  THEN
          time_dosp       = time_dosp        + dt_3d
       ENDIF
       time_dots          = time_dots        + dt_3d
       IF ( .NOT. first_call_lpm )  THEN
          time_dopts      = time_dopts       + dt_3d
       ENDIF
       IF ( time_since_reference_point >= skip_time_dopr )  THEN
          time_dopr       = time_dopr        + dt_3d
       ENDIF
       time_dopr_listing  = time_dopr_listing + dt_3d
       time_run_control   = time_run_control + dt_3d
!
!--    Increment time-counter for surface output
       IF ( surface_output )  THEN
          IF ( time_since_reference_point >= skip_time_dosurf )  THEN
             time_dosurf    = time_dosurf + dt_3d
          ENDIF
          IF ( time_since_reference_point >= skip_time_dosurf_av )  THEN
             time_dosurf_av = time_dosurf_av + dt_3d
          ENDIF
       ENDIF
!
!--    Increment time-counter for virtual measurements
       IF ( virtual_measurement  .AND.  vm_time_start <= time_since_reference_point )  THEN
          time_virtual_measurement = time_virtual_measurement + dt_3d
       ENDIF

!
!--    Increment time-counter for wind turbine data output
       IF ( wind_turbine )  THEN
          time_wtm = time_wtm + dt_3d
       ENDIF

!
!--    In case of synthetic turbulence generation and parametrized turbulence
!--    information, update the time counter and if required, adjust the
!--    STG to new atmospheric conditions.
       IF ( use_syn_turb_gen  )  THEN
          IF ( parametrize_inflow_turbulence )  THEN
             time_stg_adjust = time_stg_adjust + dt_3d
             IF ( time_stg_adjust >= dt_stg_adjust )  THEN
                CALL cpu_log( log_point(57), 'synthetic_turbulence_gen', 'start' )
                CALL stg_adjust
                CALL cpu_log( log_point(57), 'synthetic_turbulence_gen', 'stop' )
             ENDIF
          ENDIF
          time_stg_call = time_stg_call + dt_3d
       ENDIF

!
!--    Data exchange between coupled models
       IF ( coupling_mode /= 'uncoupled'  .AND.  run_coupled  .AND.  .NOT. vnested )  THEN
          time_coupling = time_coupling + dt_3d

!
!--       In case of model termination initiated by the local model
!--       (terminate_coupled > 0), the coupler must be skipped because it would
!--       cause an MPI intercomminucation hang.
!--       If necessary, the coupler will be called at the beginning of the
!--       next restart run.
          DO WHILE ( time_coupling >= dt_coupling  .AND.  terminate_coupled == 0 )
             CALL surface_coupler
             time_coupling = time_coupling - dt_coupling
          ENDDO
       ENDIF

!
!--    Biometeorology calculation of stationary thermal indices
!--    Todo (kanani): biometeorology needs own time_... treatment.
!--                   It might be that time_do2d_xy differs from time_do3d,
!--                   and then we might get trouble with the biomet output,
!--                   because we can have 2d and/or 3d biomet output!!
       IF ( biometeorology                                                                         &
            .AND. ( ( time_do3d >= dt_do3d  .AND.  time_since_reference_point >= skip_time_do3d )  &
                  .OR.                                                                             &
            ( time_do2d_xy >= dt_do2d_xy  .AND.  time_since_reference_point >= skip_time_do2d_xy ) &
                    ) )  THEN
!
!--       If required, do thermal comfort calculations
          IF ( thermal_comfort )  THEN
             CALL bio_calculate_thermal_index_maps ( .FALSE. )
          ENDIF
!
!--       If required, do UV exposure calculations
          IF ( uv_exposure )  THEN
             CALL bio_calculate_uv_exposure
          ENDIF
       ENDIF

!
!--    Execute alle other module actions routunes
       CALL module_interface_actions( 'after_integration' )

!
!--    If Galilei transformation is used, determine the distance that the
!--    model has moved so far
       IF ( galilei_transformation )  THEN
          advected_distance_x = advected_distance_x + u_gtrans * dt_3d
          advected_distance_y = advected_distance_y + v_gtrans * dt_3d
       ENDIF

!
!--    Check, if restart is necessary (because cpu-time is expiring or
!--    because it is forced by user) and set stop flag
!--    This call is skipped if the remote model has already initiated a restart.
       IF ( .NOT. terminate_run )  CALL check_for_restart

!
!--    Carry out statistical analysis and output at the requested output times.
!--    The MOD function is used for calculating the output time counters (like
!--    time_dopr) in order to regard a possible decrease of the output time
!--    interval in case of restart runs

!
!--    Set a flag indicating that so far no statistics have been created
!--    for this time step
       flow_statistics_called = .FALSE.

!
!--    If required, call flow_statistics for averaging in time
       IF ( averaging_interval_pr /= 0.0_wp  .AND.                                                 &
            ( dt_dopr - time_dopr ) <= averaging_interval_pr  .AND.                                &
            time_since_reference_point >= skip_time_dopr )  THEN
          time_dopr_av = time_dopr_av + dt_3d
          IF ( time_dopr_av >= dt_averaging_input_pr )  THEN
             do_sum = .TRUE.
             time_dopr_av = MOD( time_dopr_av, MAX( dt_averaging_input_pr, dt_3d ) )
          ENDIF
       ENDIF
       IF ( do_sum )  CALL flow_statistics

!
!--    Sum-up 3d-arrays for later output of time-averaged 2d/3d/masked data
       IF ( averaging_interval /= 0.0_wp  .AND.                                                    &
            ( dt_data_output_av - time_do_av ) <= averaging_interval  .AND.                        &
            time_since_reference_point >= skip_time_data_output_av )                               &
       THEN
          time_do_sla = time_do_sla + dt_3d
          IF ( time_do_sla >= dt_averaging_input )  THEN
             IF ( current_timestep_number > timestep_number_at_prev_calc )                         &
                CALL doq_calculate

             CALL sum_up_3d_data
             average_count_3d = average_count_3d + 1
             time_do_sla = MOD( time_do_sla, MAX( dt_averaging_input, dt_3d ) )
          ENDIF
       ENDIF
!
!--    Average surface data
       IF ( surface_output )  THEN
          IF ( averaging_interval_surf /= 0.0_wp                                                   &
                .AND.  ( dt_dosurf_av - time_dosurf_av ) <= averaging_interval_surf                &
                .AND.  time_since_reference_point >= skip_time_dosurf_av )  THEN
             IF ( time_dosurf_av >= dt_averaging_input )  THEN
                CALL surface_data_output_averaging
                average_count_surf = average_count_surf + 1
             ENDIF
          ENDIF
       ENDIF

!
!--    Calculate spectra for time averaging
       IF ( averaging_interval_sp /= 0.0_wp  .AND. ( dt_dosp - time_dosp ) <= averaging_interval_sp&
            .AND.  time_since_reference_point >= skip_time_dosp )  THEN
          time_dosp_av = time_dosp_av + dt_3d
          IF ( time_dosp_av >= dt_averaging_input_pr )  THEN
             CALL calc_spectra
             time_dosp_av = MOD( time_dosp_av, MAX( dt_averaging_input_pr, dt_3d ) )
          ENDIF
       ENDIF

!
!--    Call flight module and output data
       IF ( virtual_flight )  THEN
          CALL flight_measurement
          CALL data_output_flight
       ENDIF
!
!--    Take virtual measurements
       IF ( virtual_measurement  .AND.  time_virtual_measurement >= dt_virtual_measurement         &
                                 .AND.  vm_time_start <= time_since_reference_point )  THEN
          CALL vm_sampling
          CALL vm_data_output
          time_virtual_measurement = MOD(      time_virtual_measurement,                           &
                                          MAX( dt_virtual_measurement, dt_3d ) )
       ENDIF

!
!--    Output wind turbine data
       IF ( wind_turbine  .AND.  time_wtm >= dt_data_output_wtm )  THEN
          CALL wtm_data_output
          time_wtm = MOD( time_wtm, MAX( dt_data_output_wtm, dt_3d ) )
       ENDIF

!
!--    Profile output (ASCII) on file
       IF ( time_dopr_listing >= dt_dopr_listing )  THEN
          CALL print_1d
          time_dopr_listing = MOD( time_dopr_listing, MAX( dt_dopr_listing, dt_3d ) )
       ENDIF

!
!--    Graphic output for PROFIL
       IF ( time_dopr >= dt_dopr  .AND.  time_since_reference_point >= skip_time_dopr )  THEN
          IF ( dopr_n /= 0 )  CALL data_output_profiles
          time_dopr = MOD( time_dopr, MAX( dt_dopr, dt_3d ) )
          time_dopr_av = 0.0_wp    ! due to averaging (see above)
       ENDIF

!
!--    Graphic output for time series
       IF ( time_dots >= dt_dots )  THEN
          CALL data_output_tseries
          time_dots = MOD( time_dots, MAX( dt_dots, dt_3d ) )
       ENDIF

!
!--    Output of spectra (formatted for use with PROFIL), in case of no
!--    time averaging, spectra has to be calculated before
       IF ( time_dosp >= dt_dosp  .AND.  time_since_reference_point >= skip_time_dosp )  THEN
          IF ( average_count_sp == 0 )  CALL calc_spectra
          CALL data_output_spectra
          time_dosp = MOD( time_dosp, MAX( dt_dosp, dt_3d ) )
       ENDIF

!
!--    2d-data output (cross-sections)
       IF ( time_do2d_xy >= dt_do2d_xy  .AND.  time_since_reference_point >= skip_time_do2d_xy )  THEN
          IF ( current_timestep_number > timestep_number_at_prev_calc )                            &
             CALL doq_calculate

          CALL data_output_2d( 'xy', 0 )
          time_do2d_xy = MOD( time_do2d_xy, MAX( dt_do2d_xy, dt_3d ) )
       ENDIF
       IF ( time_do2d_xz >= dt_do2d_xz  .AND.  time_since_reference_point >= skip_time_do2d_xz )  THEN
          IF ( current_timestep_number > timestep_number_at_prev_calc )                            &

             CALL doq_calculate
          CALL data_output_2d( 'xz', 0 )
          time_do2d_xz = MOD( time_do2d_xz, MAX( dt_do2d_xz, dt_3d ) )
       ENDIF
       IF ( time_do2d_yz >= dt_do2d_yz  .AND.  time_since_reference_point >= skip_time_do2d_yz )  THEN
          IF ( current_timestep_number > timestep_number_at_prev_calc )                            &
             CALL doq_calculate

          CALL data_output_2d( 'yz', 0 )
          time_do2d_yz = MOD( time_do2d_yz, MAX( dt_do2d_yz, dt_3d ) )
       ENDIF

!
!--    3d-data output (volume data)
       IF ( time_do3d >= dt_do3d  .AND.  time_since_reference_point >= skip_time_do3d )  THEN
          IF ( current_timestep_number > timestep_number_at_prev_calc )                            &
             CALL doq_calculate

          CALL data_output_3d( 0 )
          time_do3d = MOD( time_do3d, MAX( dt_do3d, dt_3d ) )
       ENDIF

!
!--    Masked data output
       DO  mid = 1, masks
          IF ( time_domask(mid) >= dt_domask(mid)                                                  &
               .AND.  time_since_reference_point >= skip_time_domask(mid) )  THEN
             IF ( current_timestep_number > timestep_number_at_prev_calc )                         &
                CALL doq_calculate

             CALL data_output_mask( 0, mid )
             time_domask(mid) = MOD( time_domask(mid), MAX( dt_domask(mid), dt_3d ) )
          ENDIF
       ENDDO

!
!--    Output of time-averaged 2d/3d/masked data
       IF ( time_do_av >= dt_data_output_av                                                        &
            .AND.  time_since_reference_point >= skip_time_data_output_av )  THEN
          CALL average_3d_data
!
!--       Udate thermal comfort indices based on updated averaged input
          IF ( biometeorology  .AND.  thermal_comfort )  THEN
             CALL bio_calculate_thermal_index_maps ( .TRUE. )
          ENDIF
          CALL data_output_2d( 'xy', 1 )
          CALL data_output_2d( 'xz', 1 )
          CALL data_output_2d( 'yz', 1 )
          CALL data_output_3d( 1 )
          DO  mid = 1, masks
             CALL data_output_mask( 1, mid )
          ENDDO
          time_do_av = MOD( time_do_av, MAX( dt_data_output_av, dt_3d ) )
       ENDIF
!
!--    Output of surface data, instantaneous and averaged data
       IF ( surface_output )  THEN
          IF ( time_dosurf >= dt_dosurf  .AND.  time_since_reference_point >= skip_time_dosurf )  THEN
             CALL surface_data_output( 0 )
             time_dosurf = MOD( time_dosurf, MAX( dt_dosurf, dt_3d ) )
          ENDIF
          IF ( time_dosurf_av >= dt_dosurf_av  .AND.  time_since_reference_point >= skip_time_dosurf_av )  THEN
             CALL surface_data_output( 1 )
             time_dosurf_av = MOD( time_dosurf_av, MAX( dt_dosurf_av, dt_3d ) )
          ENDIF
       ENDIF

!
!--    Output of particle time series
       IF ( particle_advection )  THEN
          IF ( time_dopts >= dt_dopts  .OR.                                                        &
               ( time_since_reference_point >= particle_advection_start  .AND.                     &
                 first_call_lpm ) )  THEN
             CALL lpm_data_output_ptseries
             time_dopts = MOD( time_dopts, MAX( dt_dopts, dt_3d ) )
          ENDIF
       ENDIF

!
!--    If required, set the heat flux for the next time step to a random value
       IF ( constant_heatflux  .AND.  random_heatflux )  THEN
          IF ( surf_def_h(0)%ns >= 1 )  THEN
             CALL cpu_log( log_point(23), 'disturb_heatflux', 'start' )
             CALL disturb_heatflux( surf_def_h(0) )
             CALL cpu_log( log_point(23), 'disturb_heatflux', 'stop' )
          ENDIF
          IF ( surf_lsm_h%ns    >= 1 )  THEN
             CALL cpu_log( log_point(23), 'disturb_heatflux', 'start' )
             CALL disturb_heatflux( surf_lsm_h    )
             CALL cpu_log( log_point(23), 'disturb_heatflux', 'stop' )
          ENDIF
          IF ( surf_usm_h%ns    >= 1 )  THEN
             CALL cpu_log( log_point(23), 'disturb_heatflux', 'start' )
             CALL disturb_heatflux( surf_usm_h    )
             CALL cpu_log( log_point(23), 'disturb_heatflux', 'stop' )
          ENDIF
       ENDIF

!
!--    Execute alle other module actions routunes
       CALL module_interface_actions( 'after_timestep' )

!
!--    Determine size of next time step. Save timestep dt_3d because it is
!--    newly calculated in routine timestep, but required further below for
!--    steering the run control output interval
       dt_3d_old = dt_3d
       CALL timestep

#if defined( __parallel )
!
!--    Synchronize the timestep in case of nested run.
       IF ( nested_run )  THEN
!
!--       Synchronize by unifying the time step.
!--       Global minimum of all time-steps is used for all.
          CALL pmci_synchronize
       ENDIF
#endif

!
!--    Computation and output of run control parameters.
!--    This is also done whenever perturbations have been imposed
       IF ( time_run_control >= dt_run_control  .OR.                                               &
            timestep_scheme(1:5) /= 'runge'  .OR.  disturbance_created )                           &
       THEN
          CALL run_control
          IF ( time_run_control >= dt_run_control )  THEN
             time_run_control = MOD( time_run_control, MAX( dt_run_control, dt_3d_old ) )
          ENDIF
       ENDIF

!
!--    Output elapsed simulated time in form of a progress bar on stdout
       IF ( myid == 0 )  CALL output_progress_bar

       CALL cpu_log( log_point_s(10), 'timesteps', 'stop' )


    ENDDO   ! time loop

#if defined( _OPENACC )
    CALL exit_surface_arrays
#endif
!$ACC END DATA
!$ACC END DATA
!$ACC END DATA
!$ACC END DATA
!$ACC END DATA
!$ACC END DATA
!$ACC END DATA

#if defined( __parallel )
!
!-- Vertical nesting: Deallocate variables initialized for vertical nesting
    IF ( vnest_init )  CALL vnest_deallocate
#endif

    IF ( myid == 0 )  CALL finish_progress_bar

    CALL location_message( 'atmosphere (and/or ocean) time-stepping', 'finished' )

 END SUBROUTINE time_integration
