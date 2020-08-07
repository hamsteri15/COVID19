!> @file bulk_cloud_model_mod.f90
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
! $Id: bulk_cloud_model_mod.f90 4457 2020-03-11 14:20:43Z raasch $
! use statement for exchange horiz added
! 
! 4418 2020-02-21 09:41:13Z raasch
! bugfix for raindrop number adjustment
! 
! 4370 2020-01-10 14:00:44Z raasch
! vector directives added to force vectorization on Intel19 compiler
! 
! 4360 2020-01-07 11:25:50Z suehring
! Introduction of wall_flags_total_0, which currently sets bits based on static
! topography information used in wall_flags_static_0
! 
! 4329 2019-12-10 15:46:36Z motisi
! Renamed wall_flags_0 to wall_flags_static_0
! 
! 4289 2019-11-05 14:33:41Z knoop
! Removed parameters precipitation and precipitation_amount_interval from namelist
! 
! 4268 2019-10-17 11:29:38Z schwenkel
! Introducing bcm_boundary_conditions
! 
! 4182 2019-08-22 15:20:23Z scharf
! Corrected "Former revisions" section
! 
! 4168 2019-08-16 13:50:17Z suehring
! Replace function get_topography_top_index by topo_top_ind
! 
! 4110 2019-07-22 17:05:21Z suehring
! Pass integer flag array as well as boundary flags to WS scalar advection 
! routine
! 
! 4109 2019-07-22 17:00:34Z suehring
! Added microphyics scheme 'morrision_no_rain'
! 
! 3931 2019-04-24 16:34:28Z schwenkel
! Added bcm_exchange_horiz which is called after non_transport_physics
! 
! 3885 2019-04-11 11:29:34Z kanani
! Changes related to global restructuring of location messages and introduction 
! of additional debug messages
! 
! 3874 2019-04-08 16:53:48Z knoop
! Implemented non_transport_physics module interfaces
! 
! 3870 2019-04-08 13:44:34Z knoop
! Moving prognostic equations of bcm into bulk_cloud_model_mod
! 
! 3869 2019-04-08 11:54:20Z knoop
! moving the furniture around ;-)
! 
! 3786 2019-03-06 16:58:03Z raasch
! unsed variables removed
! 
! 3767 2019-02-27 08:18:02Z raasch
! unused variable for file index removed from rrd-subroutines parameter list
! 
! 3724 2019-02-06 16:28:23Z kanani
! Correct double-used log_point_s unit
! 
! 3700 2019-01-26 17:03:42Z knoop
! nopointer option removed
! 1053 2012-11-13 17:11:03Z hoffmann
! initial revision
!
! Description:
! ------------
!> Calculate bulk cloud microphysics.
!------------------------------------------------------------------------------!
 MODULE bulk_cloud_model_mod


    USE advec_s_bc_mod,                                                        &
        ONLY:  advec_s_bc

    USE advec_s_pw_mod,                                                        &
        ONLY:  advec_s_pw

    USE advec_s_up_mod,                                                        &
        ONLY:  advec_s_up

    USE advec_ws,                                                              &
        ONLY:  advec_s_ws

    USE arrays_3d,                                                             &
        ONLY:  ddzu, diss, dzu, dzw, hyp, hyrho,                               &
               nc, nc_1, nc_2, nc_3, nc_p, nr, nr_1, nr_2, nr_3, nr_p,         &
               precipitation_amount, prr, pt, d_exner, pt_init, q, ql, ql_1,   &
               qc, qc_1, qc_2, qc_3, qc_p, qr, qr_1, qr_2, qr_3, qr_p,         &
               exner, zu, tnc_m, tnr_m, tqc_m, tqr_m, tend, rdf_sc, &
               flux_l_qc, flux_l_qr, flux_l_nc, flux_l_nr, &
               flux_s_qc, flux_s_qr, flux_s_nc, flux_s_nr, &
               diss_l_qc, diss_l_qr, diss_l_nc, diss_l_nr, &
               diss_s_qc, diss_s_qr, diss_s_nc, diss_s_nr

    USE averaging,                                                             &
        ONLY:  nc_av, nr_av, prr_av, qc_av, ql_av, qr_av

    USE basic_constants_and_equations_mod,                                     &
        ONLY:  c_p, g, lv_d_cp, lv_d_rd, l_v, magnus, molecular_weight_of_solute,&
               molecular_weight_of_water, pi, rho_l, rho_s, r_d, r_v, vanthoff,&
               exner_function, exner_function_invers, ideal_gas_law_rho,       &
               ideal_gas_law_rho_pt, barometric_formula, rd_d_rv

    USE control_parameters,                                                    &
        ONLY:  bc_dirichlet_l,                                                 &
               bc_dirichlet_n,                                                 &
               bc_dirichlet_r,                                                 &
               bc_dirichlet_s,                                                 &
               bc_radiation_l,                                                 &
               bc_radiation_n,                                                 &
               bc_radiation_r,                                                 &
               bc_radiation_s,                                                 &
               debug_output,                                                   &
               dt_3d, dt_do2d_xy, intermediate_timestep_count,                 &
               intermediate_timestep_count_max, large_scale_forcing,           &
               lsf_surf, pt_surface, rho_surface, surface_pressure,            &
               time_do2d_xy, message_string, initializing_actions,             &
               ws_scheme_sca, scalar_advec, timestep_scheme, tsc, loop_optimization

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point, log_point_s

    USE diffusion_s_mod,                                                       &
        ONLY:  diffusion_s

    USE grid_variables,                                                        &
        ONLY:  dx, dy

    USE indices,                                                               &
        ONLY:  advc_flags_s,                                                   &
               nbgp, nxl, nxlg, nxr, nxrg, nys, nysg, nyn, nyng, nzb, nzt,     &
               topo_top_ind,                                                   &
               wall_flags_total_0

    USE kinds

    USE pegrid,                                                                &
        ONLY:  threads_per_task

    USE statistics,                                                            &
        ONLY:  weight_pres, weight_substep, sums_wsncs_ws_l, sums_wsnrs_ws_l,  &
               sums_wsqcs_ws_l, sums_wsqrs_ws_l

    USE surface_mod,                                                           &
        ONLY :  bc_h,                                                          &
                surf_bulk_cloud_model,                                         &
                surf_microphysics_morrison, surf_microphysics_seifert, &
                surf_def_h, surf_def_v, surf_lsm_h, surf_lsm_v, surf_usm_h, surf_usm_v

    IMPLICIT NONE

    CHARACTER (LEN=20)   ::  aerosol_bulk = 'nacl'                        !< namelist parameter
    CHARACTER (LEN=20)   ::  cloud_scheme = 'saturation_adjust'           !< namelist parameter

    LOGICAL ::  aerosol_nacl =.TRUE.                             !< nacl aerosol for bulk scheme
    LOGICAL ::  aerosol_c3h4o4 =.FALSE.                          !< malonic acid aerosol for bulk scheme
    LOGICAL ::  aerosol_nh4no3 =.FALSE.                          !< malonic acid aerosol for bulk scheme

    LOGICAL ::  bulk_cloud_model = .FALSE.                       !< namelist parameter

    LOGICAL ::  cloud_water_sedimentation = .FALSE.       !< cloud water sedimentation
    LOGICAL ::  curvature_solution_effects_bulk = .FALSE. !< flag for considering koehler theory
    LOGICAL ::  limiter_sedimentation = .TRUE.            !< sedimentation limiter
    LOGICAL ::  collision_turbulence = .FALSE.            !< turbulence effects
    LOGICAL ::  ventilation_effect = .TRUE.               !< ventilation effect

    LOGICAL ::  call_microphysics_at_all_substeps = .FALSE.      !< namelist parameter
    LOGICAL ::  microphysics_sat_adjust = .FALSE.                !< use saturation adjust bulk scheme?
    LOGICAL ::  microphysics_kessler = .FALSE.                   !< use kessler bulk scheme?
    LOGICAL ::  microphysics_morrison = .FALSE.                  !< use 2-moment Morrison (add. prog. eq. for nc and qc)
    LOGICAL ::  microphysics_seifert = .FALSE.                   !< use 2-moment Seifert and Beheng scheme
    LOGICAL ::  microphysics_morrison_no_rain = .FALSE.          !< use 2-moment Morrison     
    LOGICAL ::  precipitation = .FALSE.                          !< namelist parameter

    REAL(wp) ::  precipitation_amount_interval = 9999999.9_wp    !< namelist parameter

    REAL(wp) ::  a_1 = 8.69E-4_wp          !< coef. in turb. parametrization (cm-2 s3)
    REAL(wp) ::  a_2 = -7.38E-5_wp         !< coef. in turb. parametrization (cm-2 s3)
    REAL(wp) ::  a_3 = -1.40E-2_wp         !< coef. in turb. parametrization
    REAL(wp) ::  a_term = 9.65_wp          !< coef. for terminal velocity (m s-1)
    REAL(wp) ::  a_vent = 0.78_wp          !< coef. for ventilation effect
    REAL(wp) ::  b_1 = 11.45E-6_wp         !< coef. in turb. parametrization (m)
    REAL(wp) ::  b_2 = 9.68E-6_wp          !< coef. in turb. parametrization (m)
    REAL(wp) ::  b_3 = 0.62_wp             !< coef. in turb. parametrization
    REAL(wp) ::  b_term = 9.8_wp           !< coef. for terminal velocity (m s-1)
    REAL(wp) ::  b_vent = 0.308_wp         !< coef. for ventilation effect
    REAL(wp) ::  beta_cc = 3.09E-4_wp      !< coef. in turb. parametrization (cm-2 s3)
    REAL(wp) ::  c_1 = 4.82E-6_wp          !< coef. in turb. parametrization (m)
    REAL(wp) ::  c_2 = 4.8E-6_wp           !< coef. in turb. parametrization (m)
    REAL(wp) ::  c_3 = 0.76_wp             !< coef. in turb. parametrization
    REAL(wp) ::  c_const = 0.93_wp         !< const. in Taylor-microscale Reynolds number
    REAL(wp) ::  c_evap = 0.7_wp           !< constant in evaporation
    REAL(wp) ::  c_term = 600.0_wp         !< coef. for terminal velocity (m-1)
    REAL(wp) ::  diff_coeff_l = 0.23E-4_wp  !< diffusivity of water vapor (m2 s-1)
    REAL(wp) ::  eps_sb = 1.0E-10_wp       !< threshold in two-moments scheme
    REAL(wp) ::  eps_mr = 0.0_wp           !< threshold for morrison scheme
    REAL(wp) ::  k_cc = 9.44E09_wp         !< const. cloud-cloud kernel (m3 kg-2 s-1)
    REAL(wp) ::  k_cr0 = 4.33_wp           !< const. cloud-rain kernel (m3 kg-1 s-1)
    REAL(wp) ::  k_rr = 7.12_wp            !< const. rain-rain kernel (m3 kg-1 s-1)
    REAL(wp) ::  k_br = 1000.0_wp          !< const. in breakup parametrization (m-1)
    REAL(wp) ::  k_st = 1.2E8_wp           !< const. in drizzle parametrization (m-1 s-1)
    REAL(wp) ::  kin_vis_air = 1.4086E-5_wp  !< kin. viscosity of air (m2 s-1)
    REAL(wp) ::  prec_time_const = 0.001_wp  !< coef. in Kessler scheme (s-1)
    REAL(wp) ::  ql_crit = 0.0005_wp       !< coef. in Kessler scheme (kg kg-1)
    REAL(wp) ::  schmidt_p_1d3=0.8921121_wp  !< Schmidt number**0.33333, 0.71**0.33333
    REAL(wp) ::  sigma_gc = 1.3_wp         !< geometric standard deviation cloud droplets
    REAL(wp) ::  thermal_conductivity_l = 2.43E-2_wp  !< therm. cond. air (J m-1 s-1 K-1)
    REAL(wp) ::  w_precipitation = 9.65_wp  !< maximum terminal velocity (m s-1)
    REAL(wp) ::  x0 = 2.6E-10_wp           !< separating drop mass (kg)
!    REAL(wp) ::  xamin = 5.24E-19_wp       !< average aerosol mass (kg) (~ 0.05µm)
    REAL(wp) ::  xcmin = 4.18E-15_wp       !< minimum cloud drop size (kg) (~ 1µm)
    REAL(wp) ::  xrmin = 2.6E-10_wp        !< minimum rain drop size (kg)
    REAL(wp) ::  xrmax = 5.0E-6_wp         !< maximum rain drop site (kg)

    REAL(wp) ::  c_sedimentation = 2.0_wp        !< Courant number of sedimentation process
    REAL(wp) ::  dpirho_l                        !< 6.0 / ( pi * rho_l )
    REAL(wp) ::  dry_aerosol_radius = 0.05E-6_wp !< dry aerosol radius
    REAL(wp) ::  dt_micro                        !< microphysics time step
    REAL(wp) ::  sigma_bulk = 2.0_wp             !< width of aerosol spectrum
    REAL(wp) ::  na_init = 100.0E6_wp            !< Total particle/aerosol concentration (cm-3)
    REAL(wp) ::  nc_const = 70.0E6_wp            !< cloud droplet concentration
    REAL(wp) ::  dt_precipitation = 100.0_wp     !< timestep precipitation (s)
    REAL(wp) ::  sed_qc_const                    !< const. for sedimentation of cloud water
    REAL(wp) ::  pirho_l                         !< pi * rho_l / 6.0;

    REAL(wp) ::  e_s     !< saturation water vapor pressure
    REAL(wp) ::  q_s     !< saturation mixing ratio
    REAL(wp) ::  sat     !< supersaturation
    REAL(wp) ::  t_l     !< actual temperature

    SAVE

    PRIVATE

    PUBLIC bcm_parin, &
           bcm_check_parameters, &
           bcm_check_data_output, &
           bcm_check_data_output_pr, &
           bcm_init_arrays, &
           bcm_init, &
           bcm_header, &
           bcm_actions, &
           bcm_non_advective_processes, &
           bcm_exchange_horiz, &
           bcm_prognostic_equations, &
           bcm_boundary_conditions, &
           bcm_3d_data_averaging, &
           bcm_data_output_2d, &
           bcm_data_output_3d, &
           bcm_swap_timelevel, &
           bcm_rrd_global, &
           bcm_rrd_local, &
           bcm_wrd_global, &
           bcm_wrd_local, &
           calc_liquid_water_content

    PUBLIC call_microphysics_at_all_substeps, &
           cloud_water_sedimentation, &
           bulk_cloud_model, &
           cloud_scheme, &
           collision_turbulence, &
           dt_precipitation, &
           microphysics_morrison, &
           microphysics_morrison_no_rain, &           
           microphysics_sat_adjust, &
           microphysics_seifert, &
           na_init, &
           nc_const, &
           precipitation, &
           sigma_gc


    INTERFACE bcm_parin
       MODULE PROCEDURE bcm_parin
    END INTERFACE bcm_parin

    INTERFACE bcm_check_parameters
       MODULE PROCEDURE bcm_check_parameters
    END INTERFACE bcm_check_parameters

    INTERFACE bcm_check_data_output
       MODULE PROCEDURE bcm_check_data_output
    END INTERFACE bcm_check_data_output

    INTERFACE bcm_check_data_output_pr
       MODULE PROCEDURE bcm_check_data_output_pr
    END INTERFACE bcm_check_data_output_pr

    INTERFACE bcm_init_arrays
       MODULE PROCEDURE bcm_init_arrays
    END INTERFACE bcm_init_arrays

    INTERFACE bcm_init
       MODULE PROCEDURE bcm_init
    END INTERFACE bcm_init

    INTERFACE bcm_header
       MODULE PROCEDURE bcm_header
    END INTERFACE bcm_header

    INTERFACE bcm_actions
       MODULE PROCEDURE bcm_actions
       MODULE PROCEDURE bcm_actions_ij
    END INTERFACE bcm_actions

    INTERFACE bcm_non_advective_processes
       MODULE PROCEDURE bcm_non_advective_processes
       MODULE PROCEDURE bcm_non_advective_processes_ij
    END INTERFACE bcm_non_advective_processes

    INTERFACE bcm_exchange_horiz
       MODULE PROCEDURE bcm_exchange_horiz
    END INTERFACE bcm_exchange_horiz

    INTERFACE bcm_prognostic_equations
       MODULE PROCEDURE bcm_prognostic_equations
       MODULE PROCEDURE bcm_prognostic_equations_ij
    END INTERFACE bcm_prognostic_equations

    INTERFACE bcm_boundary_conditions
       MODULE PROCEDURE bcm_boundary_conditions
    END INTERFACE bcm_boundary_conditions

    INTERFACE bcm_swap_timelevel
       MODULE PROCEDURE bcm_swap_timelevel
    END INTERFACE bcm_swap_timelevel

    INTERFACE bcm_3d_data_averaging
       MODULE PROCEDURE bcm_3d_data_averaging
    END INTERFACE bcm_3d_data_averaging

    INTERFACE bcm_data_output_2d
       MODULE PROCEDURE bcm_data_output_2d
    END INTERFACE bcm_data_output_2d

    INTERFACE bcm_data_output_3d
       MODULE PROCEDURE bcm_data_output_3d
    END INTERFACE bcm_data_output_3d

    INTERFACE bcm_rrd_global
       MODULE PROCEDURE bcm_rrd_global
    END INTERFACE bcm_rrd_global

    INTERFACE bcm_rrd_local
       MODULE PROCEDURE bcm_rrd_local
    END INTERFACE bcm_rrd_local

    INTERFACE bcm_wrd_global
       MODULE PROCEDURE bcm_wrd_global
    END INTERFACE bcm_wrd_global

    INTERFACE bcm_wrd_local
       MODULE PROCEDURE bcm_wrd_local
    END INTERFACE bcm_wrd_local

    INTERFACE calc_liquid_water_content
       MODULE PROCEDURE calc_liquid_water_content
    END INTERFACE calc_liquid_water_content

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Parin for &bulk_cloud_parameters for the bulk cloud module
!------------------------------------------------------------------------------!
    SUBROUTINE bcm_parin


       IMPLICIT NONE

       CHARACTER (LEN=80)  ::  line  !< dummy string that contains the current line of the parameter file

       NAMELIST /bulk_cloud_parameters/  &
          aerosol_bulk, &
          c_sedimentation, &
          call_microphysics_at_all_substeps, &
          bulk_cloud_model, &
          cloud_scheme, &
          cloud_water_sedimentation, &
          collision_turbulence, &
          curvature_solution_effects_bulk, &
          dry_aerosol_radius, &
          limiter_sedimentation, &
          na_init, &
          nc_const, &
          sigma_bulk, &
          ventilation_effect

       line = ' '
!
!--    Try to find bulk cloud module namelist
       REWIND ( 11 )
       line = ' '
       DO   WHILE ( INDEX( line, '&bulk_cloud_parameters' ) == 0 )
          READ ( 11, '(A)', END=10 )  line
       ENDDO
       BACKSPACE ( 11 )
!
!--    Read user-defined namelist
       READ ( 11, bulk_cloud_parameters )
!
!--    Set flag that indicates that the bulk cloud module is switched on
       !bulk_cloud_model = .TRUE.

10     CONTINUE


    END SUBROUTINE bcm_parin


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check parameters routine for bulk cloud module
!------------------------------------------------------------------------------!
    SUBROUTINE bcm_check_parameters


       IMPLICIT NONE
!
!--    Check cloud scheme
!--    This scheme considers only saturation adjustment,
!--    i.e. water vapor surplus is converted into liquid
!--    water. No other microphysical processes are considered
       IF ( cloud_scheme == 'saturation_adjust' )  THEN
          microphysics_sat_adjust = .TRUE.
          microphysics_seifert    = .FALSE.
          microphysics_kessler    = .FALSE.
          precipitation           = .FALSE.
          microphysics_morrison_no_rain = .FALSE.
!
!--    This scheme includes all process of the seifert 
!--    beheng scheme (2001,2006). Especially rain processes are 
!--    considered with prognostic quantities of qr and nr. 
!--    Cloud droplet concentration is assumed to be constant and 
!--    qc is diagnostic.
!--    Technical remark: The switch 'microphysics_seifert' allocates
!--    fields of qr and nr and enables all rain processes.
       ELSEIF ( cloud_scheme == 'seifert_beheng' )  THEN
          microphysics_sat_adjust = .FALSE.
          microphysics_seifert    = .TRUE.
          microphysics_kessler    = .FALSE.
          microphysics_morrison  = .FALSE.
          precipitation           = .TRUE.
          microphysics_morrison_no_rain = .FALSE.
!
!--    The kessler scheme is a simplified scheme without any
!--    prognostic quantities for microphyical variables but 
!--    considering autoconversion.
       ELSEIF ( cloud_scheme == 'kessler' )  THEN
          microphysics_sat_adjust = .FALSE.
          microphysics_seifert    = .FALSE.
          microphysics_kessler    = .TRUE.
          microphysics_morrison   = .FALSE.
          precipitation           = .TRUE.
          microphysics_morrison_no_rain = .FALSE.
!
!--    The morrison scheme is an extension of the seifer beheng scheme
!--    including also relevant processes for cloud droplet size particles
!--    such as activation and an diagnostic mehtod for diffusional growth.
!--    I.e. here all processes of Seifert and Beheng as well as of the 
!--    morrision scheme are used. Therefore, ztis includes prognostic 
!--    quantities for qc and nc.
!--    Technical remark: The switch 'microphysics_morrison' allocates
!--    fields of qc and nc and enables diagnostic diffusional growth and 
!--    activation.
       ELSEIF ( cloud_scheme == 'morrison' )  THEN
          microphysics_sat_adjust = .FALSE.
          microphysics_seifert    = .TRUE.
          microphysics_kessler    = .FALSE.
          microphysics_morrison   = .TRUE.
          precipitation           = .TRUE.
          microphysics_morrison_no_rain = .FALSE.    
!          
!--    The 'morrision_no_rain' scheme includes only processes of morrision scheme
!--    without the rain processes of seifert beheng. Therfore, the prog. quantities
!--    of qr and nr remain unallocated. This might be appropiate for cloud in which 
!--    the size distribution is narrow, e.g. fog. 
       ELSEIF ( cloud_scheme == 'morrison_no_rain' )  THEN
          microphysics_sat_adjust = .FALSE.
          microphysics_seifert    = .FALSE.
          microphysics_kessler    = .FALSE.
          microphysics_morrison   = .TRUE.
          microphysics_morrison_no_rain = .TRUE.          
          precipitation           = .FALSE.     
       ELSE
          message_string = 'unknown cloud microphysics scheme cloud_scheme ="' // &
                           TRIM( cloud_scheme ) // '"'
          CALL message( 'check_parameters', 'PA0357', 1, 2, 0, 6, 0 )
       ENDIF



!
!--    Set the default value for the integration interval of precipitation amount
       IF ( microphysics_seifert  .OR.  microphysics_kessler )  THEN
          IF ( precipitation_amount_interval == 9999999.9_wp )  THEN
             precipitation_amount_interval = dt_do2d_xy
          ELSE
             IF ( precipitation_amount_interval > dt_do2d_xy )  THEN
                WRITE( message_string, * )  'precipitation_amount_interval = ',   &
                    precipitation_amount_interval, ' must not be larger than ',   &
                    'dt_do2d_xy = ', dt_do2d_xy
                CALL message( 'check_parameters', 'PA0090', 1, 2, 0, 6, 0 )
             ENDIF
          ENDIF
       ENDIF

       ! TODO: find better sollution for circular dependency problem
       surf_bulk_cloud_model = bulk_cloud_model
       surf_microphysics_morrison = microphysics_morrison
       surf_microphysics_seifert = microphysics_seifert

!
!--    Check aerosol
       IF ( aerosol_bulk == 'nacl' )  THEN
          aerosol_nacl   = .TRUE.
          aerosol_c3h4o4 = .FALSE.
          aerosol_nh4no3 = .FALSE.
       ELSEIF ( aerosol_bulk == 'c3h4o4' )  THEN
          aerosol_nacl   = .FALSE.
          aerosol_c3h4o4 = .TRUE.
          aerosol_nh4no3 = .FALSE.
       ELSEIF ( aerosol_bulk == 'nh4no3' )  THEN
          aerosol_nacl   = .FALSE.
          aerosol_c3h4o4 = .FALSE.
          aerosol_nh4no3 = .TRUE.
       ELSE
          message_string = 'unknown aerosol = "' // TRIM( aerosol_bulk ) // '"'
          CALL message( 'check_parameters', 'PA0469', 1, 2, 0, 6, 0 )
       ENDIF


    END SUBROUTINE bcm_check_parameters

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check data output for bulk cloud module
!------------------------------------------------------------------------------!
    SUBROUTINE bcm_check_data_output( var, unit )

       IMPLICIT NONE

       CHARACTER (LEN=*) ::  unit  !<
       CHARACTER (LEN=*) ::  var   !<

       SELECT CASE ( TRIM( var ) )

          CASE ( 'nc' )
             IF ( .NOT.  microphysics_morrison )  THEN
                message_string = 'output of "' // TRIM( var ) // '" ' //       &
                                 'requires ' //                                &
                                 'cloud_scheme = "morrison"'
                CALL message( 'check_parameters', 'PA0359', 1, 2, 0, 6, 0 )
             ENDIF
             unit = '1/m3'

          CASE ( 'nr' )
             IF ( .NOT.  microphysics_seifert )  THEN
                message_string = 'output of "' // TRIM( var ) // '" ' //       &
                                 'requires ' //                                &
                                 'cloud_scheme = "seifert_beheng"'
                CALL message( 'check_parameters', 'PA0359', 1, 2, 0, 6, 0 )
             ENDIF
             unit = '1/m3'

          CASE ( 'prr' )
             IF ( microphysics_sat_adjust )  THEN
                message_string = 'output of "' // TRIM( var ) // '" ' //       &
                                 'is not available for ' //                    &
                                 'cloud_scheme = "saturation_adjust"'
                CALL message( 'check_parameters', 'PA0423', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'kg/kg m/s'

          CASE ( 'qc' )
             unit = 'kg/kg'

          CASE ( 'qr' )
             IF ( .NOT.  microphysics_seifert ) THEN
                message_string = 'output of "' // TRIM( var ) // '" ' //       &
                                 'requires ' //                                &
                                 'cloud_scheme = "seifert_beheng"'
                CALL message( 'check_parameters', 'PA0359', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'kg/kg'

          CASE ( 'pra*' )
             IF ( .NOT. microphysics_kessler  .AND.                            &
                  .NOT. microphysics_seifert )  THEN
                message_string = 'output of "' // TRIM( var ) // '" ' //       &
                                 'requires ' //                                &
                                 'cloud_scheme = "kessler" or "seifert_beheng"'
                CALL message( 'check_parameters', 'PA0112', 1, 2, 0, 6, 0 )
             ENDIF
! TODO: find sollution (maybe connected to flow_statistics redesign?)
!              IF ( j == 1 )  THEN
!                 message_string = 'temporal averaging of precipitation ' //     &
!                           'amount "' // TRIM( var ) // '" is not possible'
!                 CALL message( 'check_parameters', 'PA0113', 1, 2, 0, 6, 0 )
!              ENDIF
             unit = 'mm'

          CASE ( 'prr*' )
             IF ( .NOT. microphysics_kessler  .AND.                            &
                  .NOT. microphysics_seifert )  THEN
                message_string = 'output of "' // TRIM( var ) // '"' //        &
                         ' requires' //                                        &
                         ' cloud_scheme = "kessler" or "seifert_beheng"'
                CALL message( 'check_parameters', 'PA0112', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'mm/s'

          CASE DEFAULT
             unit = 'illegal'

       END SELECT


    END SUBROUTINE bcm_check_data_output


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check data output of profiles for bulk cloud module
!------------------------------------------------------------------------------!
    SUBROUTINE bcm_check_data_output_pr( variable, var_count, unit, dopr_unit )

       USE arrays_3d,                                                          &
           ONLY: zu

       USE control_parameters,                                                 &
           ONLY: data_output_pr

       USE profil_parameter,                                                   &
           ONLY: dopr_index

       USE statistics,                                                         &
           ONLY: hom, statistic_regions

       IMPLICIT NONE

       CHARACTER (LEN=*) ::  unit      !<
       CHARACTER (LEN=*) ::  variable  !<
       CHARACTER (LEN=*) ::  dopr_unit !< local value of dopr_unit

       INTEGER(iwp) ::  var_count      !<
       INTEGER(iwp) ::  pr_index       !<

       SELECT CASE ( TRIM( variable ) )

! TODO: make index generic: pr_index = pr_palm+1

          CASE ( 'nc' )
             IF ( .NOT.  microphysics_morrison )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) //          &
                                 ' is not implemented for' //                  &
                                 ' cloud_scheme /= morrison'
                CALL message( 'check_parameters', 'PA0358', 1, 2, 0, 6, 0 )
             ENDIF
             pr_index = 123
             dopr_index(var_count) = pr_index
             dopr_unit     = '1/m3'
             unit = dopr_unit
             hom(:,2,pr_index,:)  = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'nr' )
             IF ( .NOT.  microphysics_seifert )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) //          &
                                 ' is not implemented for' //                  &
                                 ' cloud_scheme /= seifert_beheng'
                CALL message( 'check_parameters', 'PA0358', 1, 2, 0, 6, 0 )
             ENDIF
             pr_index = 73
             dopr_index(var_count) = pr_index
             dopr_unit     = '1/m3'
             unit = dopr_unit
             hom(:,2,pr_index,:)  = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'prr' )
             IF ( microphysics_sat_adjust )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) //          &
                                 ' is not available for' //                    &
                                 ' cloud_scheme = saturation_adjust'
                CALL message( 'check_parameters', 'PA0422', 1, 2, 0, 6, 0 )
             ENDIF 
             pr_index = 76
             dopr_index(var_count) = pr_index
             dopr_unit     = 'kg/kg m/s'
             unit = dopr_unit
             hom(:,2,pr_index,:)  = SPREAD( zu, 2, statistic_regions+1 )
          CASE ( 'qc' )
             pr_index = 75
             dopr_index(var_count) = pr_index
             dopr_unit     = 'kg/kg'
             unit = dopr_unit
             hom(:,2,pr_index,:)  = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'qr' )
             IF ( .NOT.  microphysics_seifert )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) //          &
                                 ' is not implemented for' //                  &
                                 ' cloud_scheme /= seifert_beheng'
                CALL message( 'check_parameters', 'PA0358', 1, 2, 0, 6, 0 )
             ENDIF
             pr_index = 74
             dopr_index(var_count) = pr_index
             dopr_unit     = 'kg/kg'
             unit = dopr_unit
             hom(:,2,pr_index,:)  = SPREAD( zu, 2, statistic_regions+1 )

          CASE DEFAULT
             unit = 'illegal'

       END SELECT

    END SUBROUTINE bcm_check_data_output_pr


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate bulk cloud module arrays and define pointers
!------------------------------------------------------------------------------!
    SUBROUTINE bcm_init_arrays

       USE indices,                                                            &
           ONLY:  nxlg, nxrg, nysg, nyng, nzb, nzt


       IMPLICIT NONE

!
!--    Liquid water content
       ALLOCATE ( ql_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

!
!--    3D-cloud water content
       IF ( .NOT. microphysics_morrison )  THEN
          ALLOCATE( qc_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       ENDIF
!
!--    Precipitation amount and rate (only needed if output is switched)
       ALLOCATE( precipitation_amount(nysg:nyng,nxlg:nxrg) )

!
!--    3d-precipitation rate
       ALLOCATE( prr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

       IF ( microphysics_morrison )  THEN
!
!--       3D-cloud drop water content, cloud drop concentration arrays
          ALLOCATE( nc_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                    &
                    nc_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                    &
                    nc_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                    &
                    qc_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                    &
                    qc_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                    &
                    qc_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       ENDIF

       IF ( microphysics_seifert )  THEN
!
!--       3D-rain water content, rain drop concentration arrays
          ALLOCATE( nr_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                    &
                    nr_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                    &
                    nr_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                    &
                    qr_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                    &
                    qr_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                    &
                    qr_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       ENDIF

       IF ( ws_scheme_sca )  THEN

          IF ( microphysics_morrison )  THEN
             ALLOCATE( sums_wsqcs_ws_l(nzb:nzt+1,0:threads_per_task-1) )
             ALLOCATE( sums_wsncs_ws_l(nzb:nzt+1,0:threads_per_task-1) )
             sums_wsqcs_ws_l = 0.0_wp
             sums_wsncs_ws_l = 0.0_wp
          ENDIF

          IF ( microphysics_seifert )  THEN
             ALLOCATE( sums_wsqrs_ws_l(nzb:nzt+1,0:threads_per_task-1) )
             ALLOCATE( sums_wsnrs_ws_l(nzb:nzt+1,0:threads_per_task-1) )
             sums_wsqrs_ws_l = 0.0_wp
             sums_wsnrs_ws_l = 0.0_wp
          ENDIF

       ENDIF

!
!--    Arrays needed for reasons of speed optimization for cache version.
!--    For the vector version the buffer arrays are not necessary,
!--    because the the fluxes can swapped directly inside the loops of the
!--    advection routines.
       IF ( loop_optimization /= 'vector' )  THEN

          IF ( ws_scheme_sca )  THEN

             IF ( microphysics_morrison )  THEN
                ALLOCATE( flux_s_qc(nzb+1:nzt,0:threads_per_task-1),           &
                          diss_s_qc(nzb+1:nzt,0:threads_per_task-1),           &
                          flux_s_nc(nzb+1:nzt,0:threads_per_task-1),           &
                          diss_s_nc(nzb+1:nzt,0:threads_per_task-1) )
                ALLOCATE( flux_l_qc(nzb+1:nzt,nys:nyn,0:threads_per_task-1),   &
                          diss_l_qc(nzb+1:nzt,nys:nyn,0:threads_per_task-1),   &
                          flux_l_nc(nzb+1:nzt,nys:nyn,0:threads_per_task-1),   &
                          diss_l_nc(nzb+1:nzt,nys:nyn,0:threads_per_task-1) )
             ENDIF

             IF ( microphysics_seifert )  THEN
                ALLOCATE( flux_s_qr(nzb+1:nzt,0:threads_per_task-1),           &
                          diss_s_qr(nzb+1:nzt,0:threads_per_task-1),           &
                          flux_s_nr(nzb+1:nzt,0:threads_per_task-1),           &
                          diss_s_nr(nzb+1:nzt,0:threads_per_task-1) )
                ALLOCATE( flux_l_qr(nzb+1:nzt,nys:nyn,0:threads_per_task-1),   &
                          diss_l_qr(nzb+1:nzt,nys:nyn,0:threads_per_task-1),   &
                          flux_l_nr(nzb+1:nzt,nys:nyn,0:threads_per_task-1),   &
                          diss_l_nr(nzb+1:nzt,nys:nyn,0:threads_per_task-1) )
             ENDIF

          ENDIF

       ENDIF

!
!--    Initial assignment of the pointers
       ql => ql_1
       IF ( .NOT. microphysics_morrison )  THEN
          qc => qc_1
       ENDIF
       IF ( microphysics_morrison )  THEN
          qc => qc_1;  qc_p  => qc_2;  tqc_m  => qc_3
          nc => nc_1;  nc_p  => nc_2;  tnc_m  => nc_3
       ENDIF
       IF ( microphysics_seifert )  THEN
          qr => qr_1;  qr_p  => qr_2;  tqr_m  => qr_3
          nr => nr_1;  nr_p  => nr_2;  tnr_m  => nr_3
       ENDIF


    END SUBROUTINE bcm_init_arrays


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of the bulk cloud module
!------------------------------------------------------------------------------!
    SUBROUTINE bcm_init

       IMPLICIT NONE

       INTEGER(iwp) ::  i !<
       INTEGER(iwp) ::  j !<

       IF ( debug_output )  CALL debug_message( 'bcm_init', 'start' )

       IF ( bulk_cloud_model )  THEN
          IF ( TRIM( initializing_actions ) /= 'read_restart_data' )  THEN
!
!--          Initialize the remaining quantities
             IF ( microphysics_morrison )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      qc(:,j,i) = 0.0_wp
                      nc(:,j,i) = 0.0_wp
                   ENDDO
                ENDDO
             ENDIF

             IF ( microphysics_seifert )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      qr(:,j,i) = 0.0_wp
                      nr(:,j,i) = 0.0_wp
                   ENDDO
                ENDDO
             ENDIF
!
!--          Liquid water content and precipitation amount
!--          are zero at beginning of the simulation
             ql = 0.0_wp
             qc = 0.0_wp
             precipitation_amount = 0.0_wp
             prr = 0.0_wp
!
!--          Initialize old and new time levels.
             IF ( microphysics_morrison )  THEN
                tqc_m = 0.0_wp
                tnc_m = 0.0_wp
                qc_p  = qc
                nc_p  = nc
             ENDIF
             IF ( microphysics_seifert )  THEN
                tqr_m = 0.0_wp
                tnr_m = 0.0_wp
                qr_p  = qr
                nr_p  = nr
             ENDIF
          ENDIF ! Only if not read_restart_data
!
!--       constant for the sedimentation of cloud water (2-moment cloud physics)
          sed_qc_const = k_st * ( 3.0_wp / ( 4.0_wp * pi * rho_l )                &
                             )**( 2.0_wp / 3.0_wp ) *                             &
                         EXP( 5.0_wp * LOG( sigma_gc )**2 )

!
!--       Calculate timestep according to precipitation
          IF ( microphysics_seifert )  THEN
             dt_precipitation = c_sedimentation * MINVAL( dzu(nzb+2:nzt) ) /      &
                                w_precipitation
          ENDIF

!
!--       Set constants for certain aerosol type
          IF ( microphysics_morrison )  THEN
             IF ( aerosol_nacl ) THEN
                molecular_weight_of_solute = 0.05844_wp
                rho_s                      = 2165.0_wp
                vanthoff                   = 2.0_wp
             ELSEIF ( aerosol_c3h4o4 ) THEN
                molecular_weight_of_solute = 0.10406_wp
                rho_s                      = 1600.0_wp
                vanthoff                   = 1.37_wp
             ELSEIF ( aerosol_nh4no3 ) THEN
                molecular_weight_of_solute = 0.08004_wp
                rho_s                      = 1720.0_wp
                vanthoff                   = 2.31_wp
             ENDIF
          ENDIF

!
!--       Pre-calculate frequently calculated fractions of pi and rho_l
          pirho_l  = pi * rho_l / 6.0_wp
          dpirho_l = 1.0_wp / pirho_l

          IF ( debug_output )  CALL debug_message( 'bcm_init', 'end' )

       ELSE

          IF ( debug_output )  CALL debug_message( 'bcm_init skipped', 'end' )

       ENDIF

    END SUBROUTINE bcm_init


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Header output for bulk cloud module
!------------------------------------------------------------------------------!
    SUBROUTINE bcm_header ( io )


       IMPLICIT NONE

       INTEGER(iwp), INTENT(IN) ::  io  !< Unit of the output file

!
!--    Write bulk cloud module header
       WRITE ( io, 1 )

       WRITE ( io, 2 )
       WRITE ( io, 3 )

       IF ( microphysics_kessler )  THEN
          WRITE ( io, 4 ) 'Kessler-Scheme'
       ENDIF

       IF ( microphysics_seifert )  THEN
          WRITE ( io, 4 ) 'Seifert-Beheng-Scheme'
          IF ( cloud_water_sedimentation )  WRITE ( io, 5 )
          IF ( collision_turbulence )  WRITE ( io, 6 )
          IF ( ventilation_effect )  WRITE ( io, 7 )
          IF ( limiter_sedimentation )  WRITE ( io, 8 )
       ENDIF

       WRITE ( io, 20 )
       WRITE ( io, 21 ) surface_pressure
       WRITE ( io, 22 ) r_d
       WRITE ( io, 23 ) rho_surface
       WRITE ( io, 24 ) c_p
       WRITE ( io, 25 ) l_v

       IF ( microphysics_seifert )  THEN
          WRITE ( io, 26 ) 1.0E-6_wp * nc_const
          WRITE ( io, 27 ) c_sedimentation
       ENDIF


1   FORMAT ( //' Bulk cloud module information:'/ &
               ' ------------------------------------------'/ )
2   FORMAT ( '--> Bulk scheme with liquid water potential temperature and'/ &
             '    total water content is used.' )
3   FORMAT ( '--> Condensation is parameterized via 0% - or 100% scheme.' )
4   FORMAT ( '--> Precipitation parameterization via ', A )

5   FORMAT ( '--> Cloud water sedimentation parameterization via Stokes law' )
6   FORMAT ( '--> Turbulence effects on precipitation process' )
7   FORMAT ( '--> Ventilation effects on evaporation of rain drops' )
8   FORMAT ( '--> Slope limiter used for sedimentation process' )

20  FORMAT ( '--> Essential parameters:' )
21  FORMAT ( '       Surface pressure             :   p_0   = ', F7.2, ' hPa')
22  FORMAT ( '       Gas constant                 :   R     = ', F5.1, ' J/(kg K)')
23  FORMAT ( '       Density of air               :   rho_0 = ', F6.3, ' kg/m**3')
24  FORMAT ( '       Specific heat cap.           :   c_p   = ', F6.1, ' J/(kg K)')
25  FORMAT ( '       Vapourization heat           :   L_v   = ', E9.2, ' J/kg')
26  FORMAT ( '       Droplet density              :   N_c   = ', F6.1, ' 1/cm**3' )
27  FORMAT ( '       Sedimentation Courant number :   C_s   = ', F4.1 )


    END SUBROUTINE bcm_header


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Control of microphysics for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE bcm_actions( location )


    CHARACTER (LEN=*), INTENT(IN) ::  location !< call location string

    SELECT CASE ( location )

       CASE ( 'before_timestep' )

          IF ( ws_scheme_sca )  THEN

             IF ( microphysics_morrison )  THEN
                sums_wsqcs_ws_l = 0.0_wp
                sums_wsncs_ws_l = 0.0_wp
             ENDIF
             IF ( microphysics_seifert )  THEN
                sums_wsqrs_ws_l = 0.0_wp
                sums_wsnrs_ws_l = 0.0_wp
             ENDIF

          ENDIF

       CASE DEFAULT
          CONTINUE

    END SELECT

    END SUBROUTINE bcm_actions


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Control of microphysics for grid points i,j
!------------------------------------------------------------------------------!

    SUBROUTINE bcm_actions_ij( i, j, location )


    INTEGER(iwp),      INTENT(IN) ::  i         !< grid index in x-direction
    INTEGER(iwp),      INTENT(IN) ::  j         !< grid index in y-direction
    CHARACTER (LEN=*), INTENT(IN) ::  location  !< call location string
    INTEGER(iwp)  ::  dummy  !< call location string

    IF ( bulk_cloud_model    )   dummy = i + j

    SELECT CASE ( location )

       CASE ( 'before_timestep' )

          IF ( ws_scheme_sca )  THEN

             IF ( microphysics_morrison )  THEN
                sums_wsqcs_ws_l = 0.0_wp
                sums_wsncs_ws_l = 0.0_wp
             ENDIF
             IF ( microphysics_seifert )  THEN
                sums_wsqrs_ws_l = 0.0_wp
                sums_wsnrs_ws_l = 0.0_wp
             ENDIF

          ENDIF

       CASE DEFAULT
          CONTINUE

    END SELECT


    END SUBROUTINE bcm_actions_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Control of microphysics for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE bcm_non_advective_processes


       CALL cpu_log( log_point(51), 'microphysics', 'start' )

       IF ( .NOT. microphysics_sat_adjust  .AND.         &
            ( intermediate_timestep_count == 1  .OR.                              &
              call_microphysics_at_all_substeps ) )                               &
       THEN

          IF ( large_scale_forcing  .AND.  lsf_surf ) THEN
!
!--          Calculate vertical profile of the hydrostatic pressure (hyp)
             hyp    = barometric_formula(zu, pt_surface * exner_function(surface_pressure * 100.0_wp), surface_pressure * 100.0_wp)
             d_exner = exner_function_invers(hyp)
             exner = 1.0_wp / exner_function_invers(hyp)
             hyrho  = ideal_gas_law_rho_pt(hyp, pt_init)
!
!--          Compute reference density
             rho_surface = ideal_gas_law_rho(surface_pressure * 100.0_wp, pt_surface * exner_function(surface_pressure * 100.0_wp))
          ENDIF

!
!--       Compute length of time step
          IF ( call_microphysics_at_all_substeps )  THEN
             dt_micro = dt_3d * weight_pres(intermediate_timestep_count)
          ELSE
             dt_micro = dt_3d
          ENDIF

!
!--       Reset precipitation rate
          IF ( intermediate_timestep_count == 1 )  prr = 0.0_wp

!
!--       Compute cloud physics
!--       Here the the simple kessler scheme is used.
          IF ( microphysics_kessler )  THEN
             CALL autoconversion_kessler
             IF ( cloud_water_sedimentation )  CALL sedimentation_cloud
             
!
!--       Here the seifert beheng scheme is used. Cloud concentration is assumed to
!--       a constant value an qc a diagnostic value.
          ELSEIF ( microphysics_seifert  .AND.  .NOT. microphysics_morrison )  THEN
             CALL adjust_cloud
             CALL autoconversion
             CALL accretion
             CALL selfcollection_breakup
             CALL evaporation_rain
             CALL sedimentation_rain
             IF ( cloud_water_sedimentation )  CALL sedimentation_cloud
             
!
!--       Here the morrison scheme is used. No rain processes are considered and qr and nr 
!--       are not allocated
          ELSEIF ( microphysics_morrison_no_rain  .AND.  .NOT. microphysics_seifert )  THEN
             CALL activation
             CALL condensation
             IF ( cloud_water_sedimentation )  CALL sedimentation_cloud   
             
!
!--       Here the full morrison scheme is used and all processes of Seifert and Beheng are 
!--       included
          ELSEIF ( microphysics_morrison  .AND.  microphysics_seifert )  THEN
             CALL adjust_cloud
             CALL activation
             CALL condensation
             CALL autoconversion
             CALL accretion
             CALL selfcollection_breakup
             CALL evaporation_rain
             CALL sedimentation_rain
             IF ( cloud_water_sedimentation )  CALL sedimentation_cloud            

          ENDIF

          CALL calc_precipitation_amount

       ENDIF

       CALL cpu_log( log_point(51), 'microphysics', 'stop' )

    END SUBROUTINE bcm_non_advective_processes


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Control of microphysics for grid points i,j
!------------------------------------------------------------------------------!

    SUBROUTINE bcm_non_advective_processes_ij( i, j )


       INTEGER(iwp) ::  i                 !<
       INTEGER(iwp) ::  j                 !<

       IF ( .NOT. microphysics_sat_adjust  .AND.         &
            ( intermediate_timestep_count == 1  .OR.                              &
              call_microphysics_at_all_substeps ) )                               &
       THEN

          IF ( large_scale_forcing  .AND.  lsf_surf ) THEN
!
!--          Calculate vertical profile of the hydrostatic pressure (hyp)
             hyp    = barometric_formula(zu, pt_surface * exner_function(surface_pressure * 100.0_wp), surface_pressure * 100.0_wp)
             d_exner = exner_function_invers(hyp)
             exner = 1.0_wp / exner_function_invers(hyp)
             hyrho  = ideal_gas_law_rho_pt(hyp, pt_init)
!
!--          Compute reference density
             rho_surface = ideal_gas_law_rho(surface_pressure * 100.0_wp, pt_surface * exner_function(surface_pressure * 100.0_wp))
          ENDIF

!
!--       Compute length of time step
          IF ( call_microphysics_at_all_substeps )  THEN
             dt_micro = dt_3d * weight_pres(intermediate_timestep_count)
          ELSE
             dt_micro = dt_3d
          ENDIF
!
!--       Reset precipitation rate
          IF ( intermediate_timestep_count == 1 )  prr(:,j,i) = 0.0_wp

!
!--       Compute cloud physics
!--       Here the the simple kessler scheme is used.
          IF( microphysics_kessler )  THEN
             CALL autoconversion_kessler_ij( i,j )
             IF ( cloud_water_sedimentation )  CALL sedimentation_cloud_ij( i,j )

!
!--       Here the seifert beheng scheme is used. Cloud concentration is assumed to
!--       a constant value an qc a diagnostic value.             
          ELSEIF ( microphysics_seifert  .AND.  .NOT. microphysics_morrison )  THEN
             CALL adjust_cloud_ij( i,j )
             CALL autoconversion_ij( i,j )
             CALL accretion_ij( i,j )
             CALL selfcollection_breakup_ij( i,j )
             CALL evaporation_rain_ij( i,j )
             CALL sedimentation_rain_ij( i,j )
             IF ( cloud_water_sedimentation )  CALL sedimentation_cloud_ij( i,j )
             
!
!--       Here the morrison scheme is used. No rain processes are considered and qr and nr 
!--       are not allocated
          ELSEIF ( microphysics_morrison_no_rain  .AND.  .NOT. microphysics_seifert )  THEN
             CALL activation_ij( i,j )
             CALL condensation_ij( i,j )
             IF ( cloud_water_sedimentation )  CALL sedimentation_cloud_ij( i,j ) 
             
!
!--       Here the full morrison scheme is used and all processes of Seifert and Beheng are 
!--       included             
          ELSEIF ( microphysics_morrison  .AND.  microphysics_seifert )  THEN
             CALL adjust_cloud_ij( i,j )
             CALL activation_ij( i,j )
             CALL condensation_ij( i,j )
             CALL autoconversion_ij( i,j )
             CALL accretion_ij( i,j )
             CALL selfcollection_breakup_ij( i,j )
             CALL evaporation_rain_ij( i,j )
             CALL sedimentation_rain_ij( i,j )
             IF ( cloud_water_sedimentation )  CALL sedimentation_cloud_ij( i,j )           

          ENDIF

          CALL calc_precipitation_amount_ij( i,j )

       ENDIF

    END SUBROUTINE bcm_non_advective_processes_ij
    
    
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Control of microphysics for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE bcm_exchange_horiz

       USE exchange_horiz_mod,                                                  &
           ONLY:  exchange_horiz


       IF ( .NOT. microphysics_sat_adjust  .AND.                                &
            ( intermediate_timestep_count == 1  .OR.                            &
              call_microphysics_at_all_substeps ) )                             &
       THEN
          IF ( microphysics_morrison )  THEN
             CALL exchange_horiz( nc, nbgp )
             CALL exchange_horiz( qc, nbgp )          
          ENDIF
          IF ( microphysics_seifert ) THEN
             CALL exchange_horiz( qr, nbgp )
             CALL exchange_horiz( nr, nbgp )
          ENDIF
          CALL exchange_horiz( q, nbgp )
          CALL exchange_horiz( pt, nbgp )          
       ENDIF


    END SUBROUTINE bcm_exchange_horiz
    


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Control of microphysics for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE bcm_prognostic_equations


       INTEGER(iwp) ::  i         !< grid index in x-direction
       INTEGER(iwp) ::  j         !< grid index in y-direction
       INTEGER(iwp) ::  k         !< grid index in z-direction

       REAL(wp)     ::  sbt  !<

!
!--    If required, calculate prognostic equations for cloud water content
!--    and cloud drop concentration
       IF ( microphysics_morrison )  THEN

          CALL cpu_log( log_point(67), 'qc-equation', 'start' )

!
!--       Calculate prognostic equation for cloud water content
          sbt = tsc(2)
          IF ( scalar_advec == 'bc-scheme' )  THEN

             IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--             Bott-Chlond scheme always uses Euler time step. Thus:
                sbt = 1.0_wp
             ENDIF
             tend = 0.0_wp
             CALL advec_s_bc( qc, 'qc' )

          ENDIF

!
!--       qc-tendency terms with no communication
          IF ( scalar_advec /= 'bc-scheme' )  THEN
             tend = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( ws_scheme_sca )  THEN
                   CALL advec_s_ws( advc_flags_s, qc, 'qc',                    &
                                    bc_dirichlet_l  .OR.  bc_radiation_l,      &
                                    bc_dirichlet_n  .OR.  bc_radiation_n,      &
                                    bc_dirichlet_r  .OR.  bc_radiation_r,      &
                                    bc_dirichlet_s  .OR.  bc_radiation_s )
                ELSE
                   CALL advec_s_pw( qc )
                ENDIF
             ELSE
                CALL advec_s_up( qc )
             ENDIF
          ENDIF

          CALL diffusion_s( qc,                                                &
                            surf_def_h(0)%qcsws, surf_def_h(1)%qcsws,          &
                            surf_def_h(2)%qcsws,                               &
                            surf_lsm_h%qcsws,    surf_usm_h%qcsws,             &
                            surf_def_v(0)%qcsws, surf_def_v(1)%qcsws,          &
                            surf_def_v(2)%qcsws, surf_def_v(3)%qcsws,          &
                            surf_lsm_v(0)%qcsws, surf_lsm_v(1)%qcsws,          &
                            surf_lsm_v(2)%qcsws, surf_lsm_v(3)%qcsws,          &
                            surf_usm_v(0)%qcsws, surf_usm_v(1)%qcsws,          &
                            surf_usm_v(2)%qcsws, surf_usm_v(3)%qcsws )

!
!--       Prognostic equation for cloud water content
          DO  i = nxl, nxr
             DO  j = nys, nyn
                !following directive is required to vectorize on Intel19
                !DIR$ IVDEP
                DO  k = nzb+1, nzt
                   qc_p(k,j,i) = qc(k,j,i) + ( dt_3d * ( sbt * tend(k,j,i) +   &
                                                      tsc(3) * tqc_m(k,j,i) )  &
                                                    - tsc(5) * rdf_sc(k) *     &
                                                               qc(k,j,i)       &
                                             )                                 &
                                    * MERGE( 1.0_wp, 0.0_wp,                   &
                                       BTEST( wall_flags_total_0(k,j,i), 0 )   &
                                          )
                   IF ( qc_p(k,j,i) < 0.0_wp )  qc_p(k,j,i) = 0.0_wp
                ENDDO
             ENDDO
          ENDDO

!
!--       Calculate tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         tqc_m(k,j,i) = tend(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ELSEIF ( intermediate_timestep_count < &
                      intermediate_timestep_count_max )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         tqc_m(k,j,i) =   -9.5625_wp * tend(k,j,i)             &
                                         + 5.3125_wp * tqc_m(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ENDIF

          CALL cpu_log( log_point(67), 'qc-equation', 'stop' )

          CALL cpu_log( log_point(68), 'nc-equation', 'start' )
!
!--       Calculate prognostic equation for cloud drop concentration
          sbt = tsc(2)
          IF ( scalar_advec == 'bc-scheme' )  THEN

             IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--             Bott-Chlond scheme always uses Euler time step. Thus:
                sbt = 1.0_wp
             ENDIF
             tend = 0.0_wp
             CALL advec_s_bc( nc, 'nc' )

          ENDIF

!
!--       nc-tendency terms with no communication
          IF ( scalar_advec /= 'bc-scheme' )  THEN
             tend = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( ws_scheme_sca )  THEN
                   CALL advec_s_ws( advc_flags_s, nc, 'nc',                    &
                                    bc_dirichlet_l  .OR.  bc_radiation_l,      &
                                    bc_dirichlet_n  .OR.  bc_radiation_n,      &
                                    bc_dirichlet_r  .OR.  bc_radiation_r,      &
                                    bc_dirichlet_s  .OR.  bc_radiation_s )
                ELSE
                   CALL advec_s_pw( nc )
                ENDIF
             ELSE
                CALL advec_s_up( nc )
             ENDIF
          ENDIF

          CALL diffusion_s( nc,                                                &
                            surf_def_h(0)%ncsws, surf_def_h(1)%ncsws,          &
                            surf_def_h(2)%ncsws,                               &
                            surf_lsm_h%ncsws,    surf_usm_h%ncsws,             &
                            surf_def_v(0)%ncsws, surf_def_v(1)%ncsws,          &
                            surf_def_v(2)%ncsws, surf_def_v(3)%ncsws,          &
                            surf_lsm_v(0)%ncsws, surf_lsm_v(1)%ncsws,          &
                            surf_lsm_v(2)%ncsws, surf_lsm_v(3)%ncsws,          &
                            surf_usm_v(0)%ncsws, surf_usm_v(1)%ncsws,          &
                            surf_usm_v(2)%ncsws, surf_usm_v(3)%ncsws )

!
!--       Prognostic equation for cloud drop concentration
          DO  i = nxl, nxr
             DO  j = nys, nyn
                !following directive is required to vectorize on Intel19
                !DIR$ IVDEP
                DO  k = nzb+1, nzt
                   nc_p(k,j,i) = nc(k,j,i) + ( dt_3d * ( sbt * tend(k,j,i) +   &
                                                      tsc(3) * tnc_m(k,j,i) )  &
                                                    - tsc(5) * rdf_sc(k) *     &
                                                               nc(k,j,i)       &
                                             )                                 &
                                   * MERGE( 1.0_wp, 0.0_wp,                    &
                                       BTEST( wall_flags_total_0(k,j,i), 0 )   &
                                          )
                   IF ( nc_p(k,j,i) < 0.0_wp )  nc_p(k,j,i) = 0.0_wp
                ENDDO
             ENDDO
          ENDDO

!
!--       Calculate tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         tnc_m(k,j,i) = tend(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ELSEIF ( intermediate_timestep_count < &
                      intermediate_timestep_count_max )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         tnc_m(k,j,i) =  -9.5625_wp * tend(k,j,i)             &
                                         + 5.3125_wp * tnc_m(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ENDIF

          CALL cpu_log( log_point(68), 'nc-equation', 'stop' )

       ENDIF
!
!--    If required, calculate prognostic equations for rain water content
!--    and rain drop concentration
       IF ( microphysics_seifert )  THEN

          CALL cpu_log( log_point(52), 'qr-equation', 'start' )

!
!--       Calculate prognostic equation for rain water content
          sbt = tsc(2)
          IF ( scalar_advec == 'bc-scheme' )  THEN

             IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--             Bott-Chlond scheme always uses Euler time step. Thus:
                sbt = 1.0_wp
             ENDIF
             tend = 0.0_wp
             CALL advec_s_bc( qr, 'qr' )

          ENDIF

!
!--       qr-tendency terms with no communication
          IF ( scalar_advec /= 'bc-scheme' )  THEN
             tend = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( ws_scheme_sca )  THEN
                   CALL advec_s_ws( advc_flags_s, qr, 'qr',                    &
                                    bc_dirichlet_l  .OR.  bc_radiation_l,      &
                                    bc_dirichlet_n  .OR.  bc_radiation_n,      &
                                    bc_dirichlet_r  .OR.  bc_radiation_r,      &
                                    bc_dirichlet_s  .OR.  bc_radiation_s )
                ELSE
                   CALL advec_s_pw( qr )
                ENDIF
             ELSE
                CALL advec_s_up( qr )
             ENDIF
          ENDIF

          CALL diffusion_s( qr,                                                &
                            surf_def_h(0)%qrsws, surf_def_h(1)%qrsws,          &
                            surf_def_h(2)%qrsws,                               &
                            surf_lsm_h%qrsws,    surf_usm_h%qrsws,             &
                            surf_def_v(0)%qrsws, surf_def_v(1)%qrsws,          &
                            surf_def_v(2)%qrsws, surf_def_v(3)%qrsws,          &
                            surf_lsm_v(0)%qrsws, surf_lsm_v(1)%qrsws,          &
                            surf_lsm_v(2)%qrsws, surf_lsm_v(3)%qrsws,          &
                            surf_usm_v(0)%qrsws, surf_usm_v(1)%qrsws,          &
                            surf_usm_v(2)%qrsws, surf_usm_v(3)%qrsws )

!
!--       Prognostic equation for rain water content
          DO  i = nxl, nxr
             DO  j = nys, nyn
                !following directive is required to vectorize on Intel19
                !DIR$ IVDEP
                DO  k = nzb+1, nzt
                   qr_p(k,j,i) = qr(k,j,i) + ( dt_3d * ( sbt * tend(k,j,i) +   &
                                                      tsc(3) * tqr_m(k,j,i) )  &
                                                    - tsc(5) * rdf_sc(k) *     &
                                                               qr(k,j,i)       &
                                             )                                 &
                                    * MERGE( 1.0_wp, 0.0_wp,                   &
                                       BTEST( wall_flags_total_0(k,j,i), 0 )   &
                                          )
                   IF ( qr_p(k,j,i) < 0.0_wp )  qr_p(k,j,i) = 0.0_wp
                ENDDO
             ENDDO
          ENDDO

!
!--       Calculate tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         tqr_m(k,j,i) = tend(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ELSEIF ( intermediate_timestep_count < &
                      intermediate_timestep_count_max )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         tqr_m(k,j,i) =   -9.5625_wp * tend(k,j,i)             &
                                         + 5.3125_wp * tqr_m(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ENDIF

          CALL cpu_log( log_point(52), 'qr-equation', 'stop' )
          CALL cpu_log( log_point(53), 'nr-equation', 'start' )

!
!--       Calculate prognostic equation for rain drop concentration
          sbt = tsc(2)
          IF ( scalar_advec == 'bc-scheme' )  THEN

             IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--             Bott-Chlond scheme always uses Euler time step. Thus:
                sbt = 1.0_wp
             ENDIF
             tend = 0.0_wp
             CALL advec_s_bc( nr, 'nr' )

          ENDIF

!
!--       nr-tendency terms with no communication
          IF ( scalar_advec /= 'bc-scheme' )  THEN
             tend = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( ws_scheme_sca )  THEN
                   CALL advec_s_ws( advc_flags_s, nr, 'nr',                    &
                                    bc_dirichlet_l  .OR.  bc_radiation_l,      &
                                    bc_dirichlet_n  .OR.  bc_radiation_n,      &
                                    bc_dirichlet_r  .OR.  bc_radiation_r,      &
                                    bc_dirichlet_s  .OR.  bc_radiation_s )
                ELSE
                   CALL advec_s_pw( nr )
                ENDIF
             ELSE
                CALL advec_s_up( nr )
             ENDIF
          ENDIF

          CALL diffusion_s( nr,                                                &
                            surf_def_h(0)%nrsws, surf_def_h(1)%nrsws,          &
                            surf_def_h(2)%nrsws,                               &
                            surf_lsm_h%nrsws,    surf_usm_h%nrsws,             &
                            surf_def_v(0)%nrsws, surf_def_v(1)%nrsws,          &
                            surf_def_v(2)%nrsws, surf_def_v(3)%nrsws,          &
                            surf_lsm_v(0)%nrsws, surf_lsm_v(1)%nrsws,          &
                            surf_lsm_v(2)%nrsws, surf_lsm_v(3)%nrsws,          &
                            surf_usm_v(0)%nrsws, surf_usm_v(1)%nrsws,          &
                            surf_usm_v(2)%nrsws, surf_usm_v(3)%nrsws )

!
!--       Prognostic equation for rain drop concentration
          DO  i = nxl, nxr
             DO  j = nys, nyn
                !following directive is required to vectorize on Intel19
                !DIR$ IVDEP
                DO  k = nzb+1, nzt
                   nr_p(k,j,i) = nr(k,j,i) + ( dt_3d * ( sbt * tend(k,j,i) +   &
                                                      tsc(3) * tnr_m(k,j,i) )  &
                                                    - tsc(5) * rdf_sc(k) *     &
                                                               nr(k,j,i)       &
                                             )                                 &
                                   * MERGE( 1.0_wp, 0.0_wp,                    &
                                       BTEST( wall_flags_total_0(k,j,i), 0 )   &
                                          )
                   IF ( nr_p(k,j,i) < 0.0_wp )  nr_p(k,j,i) = 0.0_wp
                ENDDO
             ENDDO
          ENDDO

!
!--       Calculate tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         tnr_m(k,j,i) = tend(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ELSEIF ( intermediate_timestep_count < &
                      intermediate_timestep_count_max )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         tnr_m(k,j,i) =  -9.5625_wp * tend(k,j,i)             &
                                         + 5.3125_wp * tnr_m(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ENDIF

          CALL cpu_log( log_point(53), 'nr-equation', 'stop' )

       ENDIF

    END SUBROUTINE bcm_prognostic_equations


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Control of microphysics for grid points i,j
!------------------------------------------------------------------------------!

    SUBROUTINE bcm_prognostic_equations_ij( i, j, i_omp_start, tn )


       INTEGER(iwp), INTENT(IN) ::  i            !< grid index in x-direction
       INTEGER(iwp), INTENT(IN) ::  j            !< grid index in y-direction
       INTEGER(iwp)             ::  k            !< grid index in z-direction
       INTEGER(iwp), INTENT(IN) ::  i_omp_start  !< first loop index of i-loop in prognostic_equations
       INTEGER(iwp), INTENT(IN) ::  tn           !< task number of openmp task

!
!--    If required, calculate prognostic equations for cloud water content
!--    and cloud drop concentration
       IF ( microphysics_morrison )  THEN
!
!--       Calculate prognostic equation for cloud water content
          tend(:,j,i) = 0.0_wp
          IF ( timestep_scheme(1:5) == 'runge' ) &
          THEN
             IF ( ws_scheme_sca )  THEN
                CALL advec_s_ws( advc_flags_s, i, j, qc, 'qc', flux_s_qc,      &
                                 diss_s_qc, flux_l_qc, diss_l_qc,              &
                                 i_omp_start, tn,                              &
                                 bc_dirichlet_l  .OR.  bc_radiation_l,         &
                                 bc_dirichlet_n  .OR.  bc_radiation_n,         &
                                 bc_dirichlet_r  .OR.  bc_radiation_r,         &
                                 bc_dirichlet_s  .OR.  bc_radiation_s )
             ELSE
                CALL advec_s_pw( i, j, qc )
             ENDIF
          ELSE
             CALL advec_s_up( i, j, qc )
          ENDIF
          CALL diffusion_s( i, j, qc,                                   &
                            surf_def_h(0)%qcsws, surf_def_h(1)%qcsws,   &
                            surf_def_h(2)%qcsws,                        &
                            surf_lsm_h%qcsws,    surf_usm_h%qcsws,      &
                            surf_def_v(0)%qcsws, surf_def_v(1)%qcsws,   &
                            surf_def_v(2)%qcsws, surf_def_v(3)%qcsws,   &
                            surf_lsm_v(0)%qcsws, surf_lsm_v(1)%qcsws,   &
                            surf_lsm_v(2)%qcsws, surf_lsm_v(3)%qcsws,   &
                            surf_usm_v(0)%qcsws, surf_usm_v(1)%qcsws,   &
                            surf_usm_v(2)%qcsws, surf_usm_v(3)%qcsws )

!
!--       Prognostic equation for cloud water content
          DO  k = nzb+1, nzt
             qc_p(k,j,i) = qc(k,j,i) + ( dt_3d *                         &
                                                ( tsc(2) * tend(k,j,i) + &
                                                  tsc(3) * tqc_m(k,j,i) )&
                                                - tsc(5) * rdf_sc(k)     &
                                                         * qc(k,j,i)     &
                                       )                                 &
                                 * MERGE( 1.0_wp, 0.0_wp,                &
                                    BTEST( wall_flags_total_0(k,j,i), 0 )&
                                        )
             IF ( qc_p(k,j,i) < 0.0_wp )  qc_p(k,j,i) = 0.0_wp
          ENDDO
!
!--       Calculate tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                DO  k = nzb+1, nzt
                   tqc_m(k,j,i) = tend(k,j,i)
                ENDDO
             ELSEIF ( intermediate_timestep_count < &
                      intermediate_timestep_count_max )  THEN
                DO  k = nzb+1, nzt
                   tqc_m(k,j,i) =   -9.5625_wp * tend(k,j,i) +           &
                                     5.3125_wp * tqc_m(k,j,i)
                ENDDO
             ENDIF
          ENDIF

!
!--       Calculate prognostic equation for cloud drop concentration.
          tend(:,j,i) = 0.0_wp
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( ws_scheme_sca )  THEN
                CALL advec_s_ws( advc_flags_s, i, j, nc, 'nc', flux_s_nc,      &
                                 diss_s_nc, flux_l_nc, diss_l_nc,              &
                                 i_omp_start, tn,                              &
                                 bc_dirichlet_l  .OR.  bc_radiation_l,         &
                                 bc_dirichlet_n  .OR.  bc_radiation_n,         &
                                 bc_dirichlet_r  .OR.  bc_radiation_r,         &
                                 bc_dirichlet_s  .OR.  bc_radiation_s )
             ELSE
                CALL advec_s_pw( i, j, nc )
             ENDIF
          ELSE
             CALL advec_s_up( i, j, nc )
          ENDIF
          CALL diffusion_s( i, j, nc,                                    &
                            surf_def_h(0)%ncsws, surf_def_h(1)%ncsws,    &
                            surf_def_h(2)%ncsws,                         &
                            surf_lsm_h%ncsws,    surf_usm_h%ncsws,       &
                            surf_def_v(0)%ncsws, surf_def_v(1)%ncsws,    &
                            surf_def_v(2)%ncsws, surf_def_v(3)%ncsws,    &
                            surf_lsm_v(0)%ncsws, surf_lsm_v(1)%ncsws,    &
                            surf_lsm_v(2)%ncsws, surf_lsm_v(3)%ncsws,    &
                            surf_usm_v(0)%ncsws, surf_usm_v(1)%ncsws,    &
                            surf_usm_v(2)%ncsws, surf_usm_v(3)%ncsws )

!
!--       Prognostic equation for cloud drop concentration
          DO  k = nzb+1, nzt
             nc_p(k,j,i) = nc(k,j,i) + ( dt_3d *                         &
                                                ( tsc(2) * tend(k,j,i) + &
                                                  tsc(3) * tnc_m(k,j,i) )&
                                                - tsc(5) * rdf_sc(k)     &
                                                         * nc(k,j,i)     &
                                       )                                 &
                                 * MERGE( 1.0_wp, 0.0_wp,                &
                                    BTEST( wall_flags_total_0(k,j,i), 0 )&
                                        )
             IF ( nc_p(k,j,i) < 0.0_wp )  nc_p(k,j,i) = 0.0_wp
          ENDDO
!
!--       Calculate tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                DO  k = nzb+1, nzt
                   tnc_m(k,j,i) = tend(k,j,i)
                ENDDO
             ELSEIF ( intermediate_timestep_count < &
                      intermediate_timestep_count_max )  THEN
                DO  k = nzb+1, nzt
                   tnc_m(k,j,i) =   -9.5625_wp * tend(k,j,i) +           &
                                     5.3125_wp * tnc_m(k,j,i)
                ENDDO
             ENDIF
          ENDIF

       ENDIF
!
!--    If required, calculate prognostic equations for rain water content
!--    and rain drop concentration
       IF ( microphysics_seifert )  THEN
!
!--             Calculate prognostic equation for rain water content
          tend(:,j,i) = 0.0_wp
          IF ( timestep_scheme(1:5) == 'runge' ) &
          THEN
             IF ( ws_scheme_sca )  THEN
                CALL advec_s_ws( advc_flags_s, i, j, qr, 'qr', flux_s_qr,      &
                                 diss_s_qr, flux_l_qr, diss_l_qr,              &
                                 i_omp_start, tn,                              &
                                 bc_dirichlet_l  .OR.  bc_radiation_l,         &
                                 bc_dirichlet_n  .OR.  bc_radiation_n,         &
                                 bc_dirichlet_r  .OR.  bc_radiation_r,         &
                                 bc_dirichlet_s  .OR.  bc_radiation_s )
             ELSE
                CALL advec_s_pw( i, j, qr )
             ENDIF
          ELSE
             CALL advec_s_up( i, j, qr )
          ENDIF
          CALL diffusion_s( i, j, qr,                                   &
                            surf_def_h(0)%qrsws, surf_def_h(1)%qrsws,   &
                            surf_def_h(2)%qrsws,                        &
                            surf_lsm_h%qrsws,    surf_usm_h%qrsws,      &
                            surf_def_v(0)%qrsws, surf_def_v(1)%qrsws,   &
                            surf_def_v(2)%qrsws, surf_def_v(3)%qrsws,   &
                            surf_lsm_v(0)%qrsws, surf_lsm_v(1)%qrsws,   &
                            surf_lsm_v(2)%qrsws, surf_lsm_v(3)%qrsws,   &
                            surf_usm_v(0)%qrsws, surf_usm_v(1)%qrsws,   &
                            surf_usm_v(2)%qrsws, surf_usm_v(3)%qrsws )

!
!--       Prognostic equation for rain water content
          DO  k = nzb+1, nzt
             qr_p(k,j,i) = qr(k,j,i) + ( dt_3d *                         &
                                                ( tsc(2) * tend(k,j,i) + &
                                                  tsc(3) * tqr_m(k,j,i) )&
                                                - tsc(5) * rdf_sc(k)     &
                                                         * qr(k,j,i)     &
                                       )                                 &
                                 * MERGE( 1.0_wp, 0.0_wp,                &
                                    BTEST( wall_flags_total_0(k,j,i), 0 )&
                                        )
             IF ( qr_p(k,j,i) < 0.0_wp )  qr_p(k,j,i) = 0.0_wp
          ENDDO
!
!--       Calculate tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                DO  k = nzb+1, nzt
                   tqr_m(k,j,i) = tend(k,j,i)
                ENDDO
             ELSEIF ( intermediate_timestep_count < &
                      intermediate_timestep_count_max )  THEN
                DO  k = nzb+1, nzt
                   tqr_m(k,j,i) =   -9.5625_wp * tend(k,j,i) +           &
                                     5.3125_wp * tqr_m(k,j,i)
                ENDDO
             ENDIF
          ENDIF

!
!--       Calculate prognostic equation for rain drop concentration.
          tend(:,j,i) = 0.0_wp
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( ws_scheme_sca )  THEN
                CALL advec_s_ws( advc_flags_s, i, j, nr, 'nr', flux_s_nr,      &
                                 diss_s_nr, flux_l_nr, diss_l_nr,              &
                                 i_omp_start, tn,                              &
                                 bc_dirichlet_l  .OR.  bc_radiation_l,         &
                                 bc_dirichlet_n  .OR.  bc_radiation_n,         &
                                 bc_dirichlet_r  .OR.  bc_radiation_r,         &
                                 bc_dirichlet_s  .OR.  bc_radiation_s )
             ELSE
                CALL advec_s_pw( i, j, nr )
             ENDIF
          ELSE
             CALL advec_s_up( i, j, nr )
          ENDIF
          CALL diffusion_s( i, j, nr,                                    &
                            surf_def_h(0)%nrsws, surf_def_h(1)%nrsws,    &
                            surf_def_h(2)%nrsws,                         &
                            surf_lsm_h%nrsws,    surf_usm_h%nrsws,       &
                            surf_def_v(0)%nrsws, surf_def_v(1)%nrsws,    &
                            surf_def_v(2)%nrsws, surf_def_v(3)%nrsws,    &
                            surf_lsm_v(0)%nrsws, surf_lsm_v(1)%nrsws,    &
                            surf_lsm_v(2)%nrsws, surf_lsm_v(3)%nrsws,    &
                            surf_usm_v(0)%nrsws, surf_usm_v(1)%nrsws,    &
                            surf_usm_v(2)%nrsws, surf_usm_v(3)%nrsws )

!
!--       Prognostic equation for rain drop concentration
          DO  k = nzb+1, nzt
             nr_p(k,j,i) = nr(k,j,i) + ( dt_3d *                         &
                                                ( tsc(2) * tend(k,j,i) + &
                                                  tsc(3) * tnr_m(k,j,i) )&
                                                - tsc(5) * rdf_sc(k)     &
                                                         * nr(k,j,i)     &
                                       )                                 &
                                 * MERGE( 1.0_wp, 0.0_wp,                &
                                    BTEST( wall_flags_total_0(k,j,i), 0 )&
                                        )
             IF ( nr_p(k,j,i) < 0.0_wp )  nr_p(k,j,i) = 0.0_wp
          ENDDO
!
!--       Calculate tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                DO  k = nzb+1, nzt
                   tnr_m(k,j,i) = tend(k,j,i)
                ENDDO
             ELSEIF ( intermediate_timestep_count < &
                      intermediate_timestep_count_max )  THEN
                DO  k = nzb+1, nzt
                   tnr_m(k,j,i) =   -9.5625_wp * tend(k,j,i) +           &
                                     5.3125_wp * tnr_m(k,j,i)
                ENDDO
             ENDIF
          ENDIF

       ENDIF

    END SUBROUTINE bcm_prognostic_equations_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Swapping of timelevels
!------------------------------------------------------------------------------!
    SUBROUTINE bcm_swap_timelevel ( mod_count )

       IMPLICIT NONE

       INTEGER, INTENT(IN) :: mod_count

       IF ( bulk_cloud_model )  THEN

          SELECT CASE ( mod_count )

             CASE ( 0 )

                IF ( microphysics_morrison )  THEN
                   qc => qc_1;    qc_p => qc_2
                   nc => nc_1;    nc_p => nc_2
                ENDIF
                IF ( microphysics_seifert )  THEN
                   qr => qr_1;    qr_p => qr_2
                   nr => nr_1;    nr_p => nr_2
                ENDIF

             CASE ( 1 )

                IF ( microphysics_morrison )  THEN
                   qc => qc_2;    qc_p => qc_1
                   nc => nc_2;    nc_p => nc_1
                ENDIF
                IF ( microphysics_seifert )  THEN
                   qr => qr_2;    qr_p => qr_1
                   nr => nr_2;    nr_p => nr_1
                ENDIF

          END SELECT

       ENDIF

    END SUBROUTINE bcm_swap_timelevel


!------------------------------------------------------------------------------!
! Description: Boundary conditions of the bulk cloud module variables
!------------------------------------------------------------------------------!
    SUBROUTINE bcm_boundary_conditions

       IMPLICIT NONE

       INTEGER(iwp) ::  i !<
       INTEGER(iwp) ::  j !<
       INTEGER(iwp) ::  k !<
       INTEGER(iwp) ::  m !<
       INTEGER(iwp) ::  l !<

       IF ( microphysics_morrison )  THEN
!
!--       Surface conditions cloud water (Dirichlet)
!--       Run loop over all non-natural and natural walls. Note, in wall-datatype
!--       the k coordinate belongs to the atmospheric grid point, therefore, set
!--       qr_p and nr_p at upward (k-1) and downward-facing (k+1) walls
          DO  l = 0, 1
          !$OMP PARALLEL DO PRIVATE( i, j, k )
             DO  m = 1, bc_h(l)%ns
                i = bc_h(l)%i(m)
                j = bc_h(l)%j(m)
                k = bc_h(l)%k(m)
                qc_p(k+bc_h(l)%koff,j,i) = 0.0_wp
                nc_p(k+bc_h(l)%koff,j,i) = 0.0_wp
             ENDDO
          ENDDO
!
!--       Top boundary condition for cloud water (Dirichlet)
          qc_p(nzt+1,:,:) = 0.0_wp
          nc_p(nzt+1,:,:) = 0.0_wp

       ENDIF

       IF ( microphysics_seifert )  THEN
!
!--       Surface conditions rain water (Dirichlet)
!--       Run loop over all non-natural and natural walls. Note, in wall-datatype
!--       the k coordinate belongs to the atmospheric grid point, therefore, set
!--       qr_p and nr_p at upward (k-1) and downward-facing (k+1) walls
          DO  l = 0, 1
          !$OMP PARALLEL DO PRIVATE( i, j, k )
             DO  m = 1, bc_h(l)%ns
                i = bc_h(l)%i(m)
                j = bc_h(l)%j(m)
                k = bc_h(l)%k(m)
                qr_p(k+bc_h(l)%koff,j,i) = 0.0_wp
                nr_p(k+bc_h(l)%koff,j,i) = 0.0_wp
             ENDDO
          ENDDO
!
!--       Top boundary condition for rain water (Dirichlet)
          qr_p(nzt+1,:,:) = 0.0_wp
          nr_p(nzt+1,:,:) = 0.0_wp

       ENDIF

!
!--    Lateral boundary conditions for scalar quantities at the outflow.
!--    Lateral oundary conditions for TKE and dissipation are set
!--    in tcm_boundary_conds.
       IF ( bc_radiation_s )  THEN
          IF ( microphysics_morrison )  THEN
             qc_p(:,nys-1,:) = qc_p(:,nys,:)
             nc_p(:,nys-1,:) = nc_p(:,nys,:)
          ENDIF
          IF ( microphysics_seifert )  THEN
             qr_p(:,nys-1,:) = qr_p(:,nys,:)
             nr_p(:,nys-1,:) = nr_p(:,nys,:)
          ENDIF
       ELSEIF ( bc_radiation_n )  THEN
          IF ( microphysics_morrison )  THEN
             qc_p(:,nyn+1,:) = qc_p(:,nyn,:)
             nc_p(:,nyn+1,:) = nc_p(:,nyn,:)
          ENDIF
          IF ( microphysics_seifert )  THEN
             qr_p(:,nyn+1,:) = qr_p(:,nyn,:)
             nr_p(:,nyn+1,:) = nr_p(:,nyn,:)
          ENDIF
       ELSEIF ( bc_radiation_l )  THEN
          IF ( microphysics_morrison )  THEN
             qc_p(:,:,nxl-1) = qc_p(:,:,nxl)
             nc_p(:,:,nxl-1) = nc_p(:,:,nxl)
          ENDIF
          IF ( microphysics_seifert )  THEN
             qr_p(:,:,nxl-1) = qr_p(:,:,nxl)
             nr_p(:,:,nxl-1) = nr_p(:,:,nxl)
          ENDIF
       ELSEIF ( bc_radiation_r )  THEN
          IF ( microphysics_morrison )  THEN
             qc_p(:,:,nxr+1) = qc_p(:,:,nxr)
             nc_p(:,:,nxr+1) = nc_p(:,:,nxr)
          ENDIF
          IF ( microphysics_seifert )  THEN
             qr_p(:,:,nxr+1) = qr_p(:,:,nxr)
             nr_p(:,:,nxr+1) = nr_p(:,:,nxr)
          ENDIF
       ENDIF

    END SUBROUTINE bcm_boundary_conditions

!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine for averaging 3D data
!------------------------------------------------------------------------------!
    SUBROUTINE bcm_3d_data_averaging( mode, variable )

       USE control_parameters,                                                 &
           ONLY:  average_count_3d

       IMPLICIT NONE

       CHARACTER (LEN=*) ::  mode     !<
       CHARACTER (LEN=*) ::  variable !<

       INTEGER(iwp) ::  i       !< local index
       INTEGER(iwp) ::  j       !< local index
       INTEGER(iwp) ::  k       !< local index

       IF ( mode == 'allocate' )  THEN

          SELECT CASE ( TRIM( variable ) )

             CASE ( 'nc' )
                IF ( .NOT. ALLOCATED( nc_av ) )  THEN
                   ALLOCATE( nc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                nc_av = 0.0_wp

             CASE ( 'nr' )
                IF ( .NOT. ALLOCATED( nr_av ) )  THEN
                   ALLOCATE( nr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                nr_av = 0.0_wp

             CASE ( 'prr' )
                IF ( .NOT. ALLOCATED( prr_av ) )  THEN
                   ALLOCATE( prr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                prr_av = 0.0_wp

             CASE ( 'qc' )
                IF ( .NOT. ALLOCATED( qc_av ) )  THEN
                   ALLOCATE( qc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                qc_av = 0.0_wp

             CASE ( 'ql' )
                IF ( .NOT. ALLOCATED( ql_av ) )  THEN
                   ALLOCATE( ql_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                ql_av = 0.0_wp

             CASE ( 'qr' )
                IF ( .NOT. ALLOCATED( qr_av ) )  THEN
                   ALLOCATE( qr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                qr_av = 0.0_wp

             CASE DEFAULT
                CONTINUE

          END SELECT

       ELSEIF ( mode == 'sum' )  THEN

          SELECT CASE ( TRIM( variable ) )

             CASE ( 'nc' )
                IF ( ALLOCATED( nc_av ) ) THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            nc_av(k,j,i) = nc_av(k,j,i) + nc(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'nr' )
                IF ( ALLOCATED( nr_av ) ) THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            nr_av(k,j,i) = nr_av(k,j,i) + nr(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'prr' )
                IF ( ALLOCATED( prr_av ) ) THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            prr_av(k,j,i) = prr_av(k,j,i) + prr(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'qc' )
                IF ( ALLOCATED( qc_av ) ) THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            qc_av(k,j,i) = qc_av(k,j,i) + qc(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'ql' )
                IF ( ALLOCATED( ql_av ) ) THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            ql_av(k,j,i) = ql_av(k,j,i) + ql(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'qr' )
                IF ( ALLOCATED( qr_av ) ) THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            qr_av(k,j,i) = qr_av(k,j,i) + qr(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE DEFAULT
                CONTINUE

          END SELECT

       ELSEIF ( mode == 'average' )  THEN

          SELECT CASE ( TRIM( variable ) )

             CASE ( 'nc' )
                IF ( ALLOCATED( nc_av ) ) THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            nc_av(k,j,i) = nc_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'nr' )
                IF ( ALLOCATED( nr_av ) ) THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            nr_av(k,j,i) = nr_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'prr' )
                IF ( ALLOCATED( prr_av ) ) THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            prr_av(k,j,i) = prr_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'qc' )
                IF ( ALLOCATED( qc_av ) ) THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            qc_av(k,j,i) = qc_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'ql' )
                IF ( ALLOCATED( ql_av ) ) THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            ql_av(k,j,i) = ql_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'qr' )
                IF ( ALLOCATED( qr_av ) ) THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            qr_av(k,j,i) = qr_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE DEFAULT
                CONTINUE

          END SELECT

       ENDIF

    END SUBROUTINE bcm_3d_data_averaging


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Define 2D output variables.
!------------------------------------------------------------------------------!
 SUBROUTINE bcm_data_output_2d( av, variable, found, grid, mode, local_pf,     &
                                two_d, nzb_do, nzt_do )


    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT(INOUT) ::  grid       !< name of vertical grid
    CHARACTER (LEN=*), INTENT(IN) ::  mode       !< either 'xy', 'xz' or 'yz'
    CHARACTER (LEN=*), INTENT(IN) ::  variable   !< name of variable

    INTEGER(iwp), INTENT(IN) ::  av        !< flag for (non-)average output
    INTEGER(iwp), INTENT(IN) ::  nzb_do    !< vertical output index (bottom)
    INTEGER(iwp), INTENT(IN) ::  nzt_do    !< vertical output index (top)

    INTEGER(iwp) ::  flag_nr   !< number of masking flag

    INTEGER(iwp) ::  i         !< loop index along x-direction
    INTEGER(iwp) ::  j         !< loop index along y-direction
    INTEGER(iwp) ::  k         !< loop index along z-direction

    LOGICAL, INTENT(INOUT) ::  found   !< flag if output variable is found
    LOGICAL, INTENT(INOUT) ::  two_d !< flag parameter that indicates 2D variables (horizontal cross sections)
    LOGICAL ::  resorted  !< flag if output is already resorted

    REAL(wp), PARAMETER ::  fill_value = -999.0_wp  !< value for the _FillValue attribute

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do), INTENT(INOUT) ::  local_pf !< local
       !< array to which output data is resorted to

    REAL(wp), DIMENSION(:,:,:), POINTER ::  to_be_resorted  !< points to selected output variable

    found = .TRUE.
    resorted = .FALSE.
!
!-- Set masking flag for topography for not resorted arrays
    flag_nr = 0    ! 0 = scalar, 1 = u, 2 = v, 3 = w

    SELECT CASE ( TRIM( variable ) )

       CASE ( 'nc_xy', 'nc_xz', 'nc_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => nc
          ELSE
             IF ( .NOT. ALLOCATED( nc_av ) ) THEN
                ALLOCATE( nc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                nc_av = REAL( fill_value, KIND = wp )
             ENDIF
             to_be_resorted => nc_av
          ENDIF
          IF ( mode == 'xy' ) grid = 'zu'

       CASE ( 'nr_xy', 'nr_xz', 'nr_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => nr
          ELSE
             IF ( .NOT. ALLOCATED( nr_av ) ) THEN
                ALLOCATE( nr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                nr_av = REAL( fill_value, KIND = wp )
             ENDIF
             to_be_resorted => nr_av
          ENDIF
          IF ( mode == 'xy' ) grid = 'zu'

       CASE ( 'pra*_xy' )        ! 2d-array / integral quantity => no av
!                CALL exchange_horiz_2d( precipitation_amount )
          DO  i = nxl, nxr
             DO  j = nys, nyn
                local_pf(i,j,nzb+1) =  precipitation_amount(j,i)
             ENDDO
          ENDDO
          precipitation_amount = 0.0_wp   ! reset for next integ. interval
          resorted = .TRUE.
          two_d = .TRUE.
          IF ( mode == 'xy' ) grid = 'zu1'

       CASE ( 'prr_xy', 'prr_xz', 'prr_yz' )
          IF ( av == 0 )  THEN
!                   CALL exchange_horiz( prr, nbgp )
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb, nzt+1
                      local_pf(i,j,k) = prr(k,j,i) * hyrho(nzb+1)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( prr_av ) ) THEN
                ALLOCATE( prr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                prr_av = REAL( fill_value, KIND = wp )
             ENDIF
!                   CALL exchange_horiz( prr_av, nbgp )
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb, nzt+1
                      local_pf(i,j,k) = prr_av(k,j,i) * hyrho(nzb+1)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          resorted = .TRUE.
          IF ( mode == 'xy' ) grid = 'zu'

       CASE ( 'qc_xy', 'qc_xz', 'qc_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => qc
          ELSE
             IF ( .NOT. ALLOCATED( qc_av ) ) THEN
                ALLOCATE( qc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                qc_av = REAL( fill_value, KIND = wp )
             ENDIF
             to_be_resorted => qc_av
          ENDIF
          IF ( mode == 'xy' ) grid = 'zu'

       CASE ( 'ql_xy', 'ql_xz', 'ql_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => ql
          ELSE
             IF ( .NOT. ALLOCATED( ql_av ) ) THEN
                ALLOCATE( ql_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ql_av = REAL( fill_value, KIND = wp )
             ENDIF
             to_be_resorted => ql_av
          ENDIF
          IF ( mode == 'xy' ) grid = 'zu'

       CASE ( 'qr_xy', 'qr_xz', 'qr_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => qr
          ELSE
             IF ( .NOT. ALLOCATED( qr_av ) ) THEN
                ALLOCATE( qr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                qr_av = REAL( fill_value, KIND = wp )
             ENDIF
             to_be_resorted => qr_av
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
                local_pf(i,j,k) = MERGE( &
                                     to_be_resorted(k,j,i),                      &
                                     REAL( fill_value, KIND = wp ),              &
                                     BTEST( wall_flags_total_0(k,j,i), flag_nr ) &
                                  )
             ENDDO
          ENDDO
       ENDDO
    ENDIF

 END SUBROUTINE bcm_data_output_2d


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Define 3D output variables.
!------------------------------------------------------------------------------!
 SUBROUTINE bcm_data_output_3d( av, variable, found, local_pf, nzb_do, nzt_do )


    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT(IN) ::  variable   !< name of variable

    INTEGER(iwp), INTENT(IN) ::  av        !< flag for (non-)average output
    INTEGER(iwp), INTENT(IN) ::  nzb_do    !< lower limit of the data output (usually 0)
    INTEGER(iwp), INTENT(IN) ::  nzt_do    !< vertical upper limit of the data output (usually nz_do3d)

    INTEGER(iwp) ::  flag_nr   !< number of masking flag

    INTEGER(iwp) ::  i         !< loop index along x-direction
    INTEGER(iwp) ::  j         !< loop index along y-direction
    INTEGER(iwp) ::  k         !< loop index along z-direction

    LOGICAL, INTENT(INOUT) ::  found     !< flag if output variable is found
    LOGICAL ::  resorted  !< flag if output is already resorted

    REAL(wp) ::  fill_value = -999.0_wp  !< value for the _FillValue attribute

    REAL(sp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do), INTENT(INOUT) ::  local_pf   !< local
       !< array to which output data is resorted to

    REAL(wp), DIMENSION(:,:,:), POINTER ::  to_be_resorted  !< points to selected output variable

    found = .TRUE.
    resorted = .FALSE.
!
!-- Set masking flag for topography for not resorted arrays
    flag_nr = 0    ! 0 = scalar, 1 = u, 2 = v, 3 = w

    SELECT CASE ( TRIM( variable ) )

       CASE ( 'nc' )
          IF ( av == 0 )  THEN
             to_be_resorted => nc
          ELSE
             IF ( .NOT. ALLOCATED( nc_av ) ) THEN
                ALLOCATE( nc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                nc_av = REAL( fill_value, KIND = wp )
             ENDIF
             to_be_resorted => nc_av
          ENDIF

       CASE ( 'nr' )
          IF ( av == 0 )  THEN
             to_be_resorted => nr
          ELSE
             IF ( .NOT. ALLOCATED( nr_av ) ) THEN
                ALLOCATE( nr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                nr_av = REAL( fill_value, KIND = wp )
             ENDIF
             to_be_resorted => nr_av
          ENDIF

       CASE ( 'prr' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = prr(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( prr_av ) ) THEN
                ALLOCATE( prr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                prr_av = REAL( fill_value, KIND = wp )
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = prr_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          resorted = .TRUE.

       CASE ( 'qc' )
          IF ( av == 0 )  THEN
             to_be_resorted => qc
          ELSE
             IF ( .NOT. ALLOCATED( qc_av ) ) THEN
                ALLOCATE( qc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                qc_av = REAL( fill_value, KIND = wp )
             ENDIF
             to_be_resorted => qc_av
          ENDIF

       CASE ( 'ql' )
          IF ( av == 0 )  THEN
             to_be_resorted => ql
          ELSE
             IF ( .NOT. ALLOCATED( ql_av ) ) THEN
                ALLOCATE( ql_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ql_av = REAL( fill_value, KIND = wp )
             ENDIF
             to_be_resorted => ql_av
          ENDIF

       CASE ( 'qr' )
          IF ( av == 0 )  THEN
             to_be_resorted => qr
          ELSE
             IF ( .NOT. ALLOCATED( qr_av ) ) THEN
                ALLOCATE( qr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                qr_av = REAL( fill_value, KIND = wp )
             ENDIF
             to_be_resorted => qr_av
          ENDIF

       CASE DEFAULT
          found = .FALSE.

    END SELECT


    IF ( found .AND. .NOT. resorted )  THEN
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb_do, nzt_do
                local_pf(i,j,k) = MERGE(                                         &
                                     to_be_resorted(k,j,i),                      &
                                     REAL( fill_value, KIND = wp ),              &
                                     BTEST( wall_flags_total_0(k,j,i), flag_nr ) &
                                  )
             ENDDO
          ENDDO
       ENDDO
    ENDIF

 END SUBROUTINE bcm_data_output_3d


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine reads the respective restart data for the bulk cloud module.
!------------------------------------------------------------------------------!
    SUBROUTINE bcm_rrd_global( found )


       USE control_parameters,                                                 &
           ONLY: length, restart_string


       IMPLICIT NONE

       LOGICAL, INTENT(OUT)  ::  found


       found = .TRUE.

       SELECT CASE ( restart_string(1:length) )

          CASE ( 'c_sedimentation' )
             READ ( 13 )  c_sedimentation

          CASE ( 'bulk_cloud_model' )
             READ ( 13 )  bulk_cloud_model

          CASE ( 'cloud_scheme' )
             READ ( 13 )  cloud_scheme

          CASE ( 'cloud_water_sedimentation' )
             READ ( 13 )  cloud_water_sedimentation

          CASE ( 'collision_turbulence' )
             READ ( 13 )  collision_turbulence

          CASE ( 'limiter_sedimentation' )
             READ ( 13 )  limiter_sedimentation

          CASE ( 'nc_const' )
             READ ( 13 )  nc_const

          CASE ( 'precipitation' )
             READ ( 13 ) precipitation

          CASE ( 'ventilation_effect' )
             READ ( 13 )  ventilation_effect

          CASE ( 'na_init' )
             READ ( 13 )  na_init

          CASE ( 'dry_aerosol_radius' )
             READ ( 13 )  dry_aerosol_radius

          CASE ( 'sigma_bulk' )
             READ ( 13 )  sigma_bulk

          CASE ( 'aerosol_bulk' )
             READ ( 13 )  aerosol_bulk

          CASE ( 'curvature_solution_effects_bulk' )
             READ ( 13 )  curvature_solution_effects_bulk


!          CASE ( 'global_paramter' )
!             READ ( 13 )  global_parameter
!          CASE ( 'global_array' )
!             IF ( .NOT. ALLOCATED( global_array ) )  ALLOCATE( global_array(1:10) )
!             READ ( 13 )  global_array

          CASE DEFAULT

             found = .FALSE.

       END SELECT


    END SUBROUTINE bcm_rrd_global


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine reads the respective restart data for the bulk cloud module.
!------------------------------------------------------------------------------!
    SUBROUTINE bcm_rrd_local( k, nxlf, nxlc, nxl_on_file, nxrf, nxrc,          &
                              nxr_on_file, nynf, nync, nyn_on_file, nysf,      &
                              nysc, nys_on_file, tmp_2d, tmp_3d, found )


       USE control_parameters

       USE indices

       USE pegrid


       IMPLICIT NONE

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

       REAL(wp), DIMENSION(nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) :: tmp_2d   !<
       REAL(wp), DIMENSION(nzb:nzt+1,nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) :: tmp_3d   !<

!
!-- Here the reading of user-defined restart data follows:
!-- Sample for user-defined output


       found = .TRUE.

       SELECT CASE ( restart_string(1:length) )

          CASE ( 'nc' )
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             nc(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =             &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'nc_av' )
             IF ( .NOT. ALLOCATED( nc_av ) )  THEN
                ALLOCATE( nc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             nc_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =          &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'nr' )
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             nr(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =             &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'nr_av' )
             IF ( .NOT. ALLOCATED( nr_av ) )  THEN
                ALLOCATE( nr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             nr_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =          &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'precipitation_amount' )
             IF ( k == 1 )  READ ( 13 )  tmp_2d
             precipitation_amount(nysc-nbgp:nync+nbgp,                   &
                                  nxlc-nbgp:nxrc+nbgp)  =                &
                   tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'prr' )
             IF ( .NOT. ALLOCATED( prr ) )  THEN
                ALLOCATE( prr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             prr(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =            &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'prr_av' )
             IF ( .NOT. ALLOCATED( prr_av ) )  THEN
                ALLOCATE( prr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             prr_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =         &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'qc' )
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             qc(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =             &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'qc_av' )
             IF ( .NOT. ALLOCATED( qc_av ) )  THEN
                ALLOCATE( qc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             qc_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =          &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'ql' )
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             ql(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =             &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'ql_av' )
             IF ( .NOT. ALLOCATED( ql_av ) )  THEN
                ALLOCATE( ql_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             ql_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =          &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'qr' )
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             qr(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =             &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'qr_av' )
             IF ( .NOT. ALLOCATED( qr_av ) )  THEN
                ALLOCATE( qr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             qr_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =          &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
!
          CASE DEFAULT

             found = .FALSE.

          END SELECT


    END SUBROUTINE bcm_rrd_local


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine writes the respective restart data for the bulk cloud module.
!------------------------------------------------------------------------------!
    SUBROUTINE bcm_wrd_global


       IMPLICIT NONE

       CALL wrd_write_string( 'c_sedimentation' )
       WRITE ( 14 )  c_sedimentation

       CALL wrd_write_string( 'bulk_cloud_model' )
       WRITE ( 14 )  bulk_cloud_model

       CALL wrd_write_string( 'cloud_scheme' )
       WRITE ( 14 )  cloud_scheme

       CALL wrd_write_string( 'cloud_water_sedimentation' )
       WRITE ( 14 )  cloud_water_sedimentation

       CALL wrd_write_string( 'collision_turbulence' )
       WRITE ( 14 )  collision_turbulence

       CALL wrd_write_string( 'limiter_sedimentation' )
       WRITE ( 14 )  limiter_sedimentation

       CALL wrd_write_string( 'nc_const' )
       WRITE ( 14 )  nc_const

       CALL wrd_write_string( 'precipitation' )
       WRITE ( 14 )  precipitation

       CALL wrd_write_string( 'ventilation_effect' )
       WRITE ( 14 )  ventilation_effect

       CALL wrd_write_string( 'na_init' )
       WRITE ( 14 )  na_init

       CALL wrd_write_string( 'dry_aerosol_radius' )
       WRITE ( 14 )  dry_aerosol_radius

       CALL wrd_write_string( 'sigma_bulk' )
       WRITE ( 14 )  sigma_bulk

       CALL wrd_write_string( 'aerosol_bulk' )
       WRITE ( 14 )  aerosol_bulk

       CALL wrd_write_string( 'curvature_solution_effects_bulk' )
       WRITE ( 14 )  curvature_solution_effects_bulk


! needs preceeding allocation if array
!       CALL wrd_write_string( 'global_parameter' )
!       WRITE ( 14 )  global_parameter

!       IF ( ALLOCATED( inflow_damping_factor ) )  THEN
!          CALL wrd_write_string( 'inflow_damping_factor' )
!          WRITE ( 14 )  inflow_damping_factor
!       ENDIF


    END SUBROUTINE bcm_wrd_global


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine writes the respective restart data for the bulk cloud module.
!------------------------------------------------------------------------------!
    SUBROUTINE bcm_wrd_local


       IMPLICIT NONE

       IF ( ALLOCATED( prr ) )  THEN
          CALL wrd_write_string( 'prr' )
          WRITE ( 14 )  prr
       ENDIF

       IF ( ALLOCATED( prr_av ) )  THEN
          CALL wrd_write_string( 'prr_av' )
          WRITE ( 14 )  prr_av
       ENDIF

       IF ( ALLOCATED( precipitation_amount ) )  THEN
          CALL wrd_write_string( 'precipitation_amount' )
          WRITE ( 14 )  precipitation_amount
       ENDIF

       CALL wrd_write_string( 'ql' )
       WRITE ( 14 )  ql

       IF ( ALLOCATED( ql_av ) )  THEN
          CALL wrd_write_string( 'ql_av' )
          WRITE ( 14 )  ql_av
       ENDIF

       CALL wrd_write_string( 'qc' )
       WRITE ( 14 )  qc

       IF ( ALLOCATED( qc_av ) )  THEN
          CALL wrd_write_string( 'qc_av' )
          WRITE ( 14 )  qc_av
       ENDIF

       IF ( microphysics_morrison )  THEN

          CALL wrd_write_string( 'nc' )
          WRITE ( 14 )  nc

          IF ( ALLOCATED( nc_av ) )  THEN
             CALL wrd_write_string( 'nc_av' )
             WRITE ( 14 )  nc_av
          ENDIF

       ENDIF

       IF ( microphysics_seifert )  THEN

          CALL wrd_write_string( 'nr' )
          WRITE ( 14 )  nr

          IF ( ALLOCATED( nr_av ) )  THEN
             CALL wrd_write_string( 'nr_av' )
             WRITE ( 14 )  nr_av
          ENDIF

          CALL wrd_write_string( 'qr' )
          WRITE ( 14 )  qr

          IF ( ALLOCATED( qr_av ) )  THEN
             CALL wrd_write_string( 'qr_av' )
             WRITE ( 14 )  qr_av
          ENDIF

       ENDIF


    END SUBROUTINE bcm_wrd_local

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Adjust number of raindrops to avoid nonlinear effects in sedimentation and
!> evaporation of rain drops due to too small or too big weights
!> of rain drops (Stevens and Seifert, 2008).
!------------------------------------------------------------------------------!
    SUBROUTINE adjust_cloud

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !<
       INTEGER(iwp) ::  j                 !<
       INTEGER(iwp) ::  k                 !<

       REAL(wp) ::  flag                  !< flag to indicate first grid level above

       CALL cpu_log( log_point_s(50), 'adjust_cloud', 'start' )

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )

                IF ( qr(k,j,i) <= eps_sb )  THEN
                   qr(k,j,i) = 0.0_wp
                   nr(k,j,i) = 0.0_wp
                ELSE
                   IF ( nr(k,j,i) * xrmin > qr(k,j,i) * hyrho(k) )  THEN
                      nr(k,j,i) = qr(k,j,i) * hyrho(k) / xrmin * flag
                   ELSEIF ( nr(k,j,i) * xrmax < qr(k,j,i) * hyrho(k) )  THEN
                      nr(k,j,i) = qr(k,j,i) * hyrho(k) / xrmax * flag
                   ENDIF
                ENDIF

                IF ( microphysics_morrison ) THEN
                   IF ( qc(k,j,i) <= eps_sb )  THEN
                      qc(k,j,i) = 0.0_wp
                      nc(k,j,i) = 0.0_wp
                   ELSE
                      IF ( nc(k,j,i) * xcmin > qc(k,j,i) * hyrho(k) )  THEN
                         nc(k,j,i) = qc(k,j,i) * hyrho(k) / xcmin * flag
                      ENDIF
                   ENDIF
                ENDIF

             ENDDO
          ENDDO
       ENDDO

       CALL cpu_log( log_point_s(50), 'adjust_cloud', 'stop' )

    END SUBROUTINE adjust_cloud

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Adjust number of raindrops to avoid nonlinear effects in
!> sedimentation and evaporation of rain drops due to too small or
!> too big weights of rain drops (Stevens and Seifert, 2008).
!> The same procedure is applied to cloud droplets if they are determined
!> prognostically. Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE adjust_cloud_ij( i, j )

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !<
       INTEGER(iwp) ::  j                 !<
       INTEGER(iwp) ::  k                 !<

       REAL(wp) ::  flag                  !< flag to indicate first grid level above surface

       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )

          IF ( qr(k,j,i) <= eps_sb )  THEN
             qr(k,j,i) = 0.0_wp
             nr(k,j,i) = 0.0_wp
          ELSE
!
!--          Adjust number of raindrops to avoid nonlinear effects in
!--          sedimentation and evaporation of rain drops due to too small or
!--          too big weights of rain drops (Stevens and Seifert, 2008).
             IF ( nr(k,j,i) * xrmin > qr(k,j,i) * hyrho(k) )  THEN
                nr(k,j,i) = qr(k,j,i) * hyrho(k) / xrmin * flag
             ELSEIF ( nr(k,j,i) * xrmax < qr(k,j,i) * hyrho(k) )  THEN
                nr(k,j,i) = qr(k,j,i) * hyrho(k) / xrmax * flag
             ENDIF

          ENDIF

          IF ( microphysics_morrison ) THEN
             IF ( qc(k,j,i) <= eps_sb )  THEN
                qc(k,j,i) = 0.0_wp
                nc(k,j,i) = 0.0_wp
             ELSE
                IF ( nc(k,j,i) * xcmin > qc(k,j,i) * hyrho(k) )  THEN
                   nc(k,j,i) = qc(k,j,i) * hyrho(k) / xcmin * flag
                ENDIF
             ENDIF
          ENDIF

       ENDDO

    END SUBROUTINE adjust_cloud_ij

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate number of activated condensation nucleii after simple activation
!> scheme of Twomey, 1959.
!------------------------------------------------------------------------------!
    SUBROUTINE activation

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !<
       INTEGER(iwp) ::  j                 !<
       INTEGER(iwp) ::  k                 !<

       REAL(wp)     ::  activ             !<
       REAL(wp)     ::  afactor           !<
       REAL(wp)     ::  beta_act          !<
       REAL(wp)     ::  bfactor           !<
       REAL(wp)     ::  k_act             !<
       REAL(wp)     ::  n_act             !<
       REAL(wp)     ::  n_ccn             !<
       REAL(wp)     ::  s_0               !<
       REAL(wp)     ::  sat_max           !<
       REAL(wp)     ::  sigma             !<
       REAL(wp)     ::  sigma_act         !<

       REAL(wp) ::  flag                  !< flag to indicate first grid level above

       CALL cpu_log( log_point_s(65), 'activation', 'start' )

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )

!
!--             Call calculation of supersaturation located in subroutine
                CALL supersaturation ( i, j, k )
!
!--             Prescribe parameters for activation
!--             (see: Bott + Trautmann, 2002, Atm. Res., 64)
                k_act  = 0.7_wp
                activ  = 0.0_wp


                IF ( sat > 0.0 .AND. .NOT. curvature_solution_effects_bulk ) THEN
!
!--                Compute the number of activated Aerosols
!--                (see: Twomey, 1959, Pure and applied Geophysics, 43)
                   n_act     = na_init * sat**k_act
!
!--                Compute the number of cloud droplets
!--                (see: Morrison + Grabowski, 2007, JAS, 64)
!                  activ = MAX( n_act - nc(k,j,i), 0.0_wp) / dt_micro

!
!--                Compute activation rate after Khairoutdinov and Kogan
!--                (see: Khairoutdinov + Kogan, 2000, Mon. Wea. Rev., 128)
                   sat_max = 1.0_wp / 100.0_wp
                   activ   = MAX( 0.0_wp, ( (na_init + nc(k,j,i) ) * MIN       &
                      ( 1.0_wp, ( sat / sat_max )**k_act) - nc(k,j,i) ) ) /    &
                       dt_micro
                ELSEIF ( sat > 0.0 .AND. curvature_solution_effects_bulk ) THEN
!
!--                Curvature effect (afactor) with surface tension
!--                parameterization by Straka (2009)
                   sigma = 0.0761_wp - 0.000155_wp * ( t_l - 273.15_wp )
                   afactor = 2.0_wp * sigma / ( rho_l * r_v * t_l )
!
!--                Solute effect (bfactor)
                   bfactor = vanthoff * molecular_weight_of_water *            &
                       rho_s / ( molecular_weight_of_solute * rho_l )

!
!--                Prescribe power index that describes the soluble fraction
!--                of an aerosol particle (beta) 
!--                (see: Morrison + Grabowski, 2007, JAS, 64)
                   beta_act  = 0.5_wp
                   sigma_act = sigma_bulk**( 1.0_wp + beta_act )
!
!--                Calculate mean geometric supersaturation (s_0) with
!--                parameterization by Khvorostyanov and Curry (2006)
                   s_0   = dry_aerosol_radius **(- ( 1.0_wp + beta_act ) ) *    &
                       ( 4.0_wp * afactor**3 / ( 27.0_wp * bfactor ) )**0.5_wp

!
!--                Calculate number of activated CCN as a function of
!--                supersaturation and taking Koehler theory into account
!--                (see: Khvorostyanov + Curry, 2006, J. Geo. Res., 111)
                   n_ccn = ( na_init / 2.0_wp ) * ( 1.0_wp - ERF(              &
                      LOG( s_0 / sat ) / ( SQRT(2.0_wp) * LOG(sigma_act) ) ) )
                   activ = MAX( ( n_ccn - nc(k,j,i) ) / dt_micro, 0.0_wp )
                ENDIF

                nc(k,j,i) = MIN( (nc(k,j,i) + activ * dt_micro * flag), na_init)

             ENDDO
          ENDDO
       ENDDO

       CALL cpu_log( log_point_s(65), 'activation', 'stop' )

    END SUBROUTINE activation

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate number of activated condensation nucleii after simple activation
!> scheme of Twomey, 1959.
!------------------------------------------------------------------------------!
    SUBROUTINE activation_ij( i, j )

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !<
       INTEGER(iwp) ::  j                 !<
       INTEGER(iwp) ::  k                 !<

       REAL(wp)     ::  activ             !<
       REAL(wp)     ::  afactor           !<
       REAL(wp)     ::  beta_act          !<
       REAL(wp)     ::  bfactor           !<
       REAL(wp)     ::  flag              !< flag to indicate first grid level above surface
       REAL(wp)     ::  k_act             !<
       REAL(wp)     ::  n_act             !<
       REAL(wp)     ::  n_ccn             !<
       REAL(wp)     ::  s_0               !<
       REAL(wp)     ::  sat_max           !<
       REAL(wp)     ::  sigma             !<
       REAL(wp)     ::  sigma_act         !<

       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )
!
!--       Call calculation of supersaturation
          CALL supersaturation ( i, j, k )
!
!--       Prescribe parameters for activation
!--       (see: Bott + Trautmann, 2002, Atm. Res., 64)
          k_act  = 0.7_wp
          activ  = 0.0_wp

          IF ( sat > 0.0 .AND. .NOT. curvature_solution_effects_bulk )  THEN
!
!--          Compute the number of activated Aerosols
!--          (see: Twomey, 1959, Pure and applied Geophysics, 43)
             n_act     = na_init * sat**k_act
!
!--          Compute the number of cloud droplets
!--          (see: Morrison + Grabowski, 2007, JAS, 64)
!            activ = MAX( n_act - nc_d1(k), 0.0_wp) / dt_micro

!
!--          Compute activation rate after Khairoutdinov and Kogan
!--          (see: Khairoutdinov + Kogan, 2000, Mon. Wea. Rev., 128)
             sat_max = 0.8_wp / 100.0_wp
             activ   = MAX( 0.0_wp, ( (na_init + nc(k,j,i) ) * MIN             &
                 ( 1.0_wp, ( sat / sat_max )**k_act) - nc(k,j,i) ) ) /         &
                  dt_micro

             nc(k,j,i) = MIN( (nc(k,j,i) + activ * dt_micro), na_init)
          ELSEIF ( sat > 0.0 .AND. curvature_solution_effects_bulk )  THEN
!
!--          Curvature effect (afactor) with surface tension
!--          parameterization by Straka (2009)
             sigma = 0.0761_wp - 0.000155_wp * ( t_l - 273.15_wp )
             afactor = 2.0_wp * sigma / ( rho_l * r_v * t_l )
!
!--          Solute effect (bfactor)
             bfactor = vanthoff * molecular_weight_of_water *                  &
                  rho_s / ( molecular_weight_of_solute * rho_l )

!
!--          Prescribe power index that describes the soluble fraction
!--          of an aerosol particle (beta).
!--          (see: Morrison + Grabowski, 2007, JAS, 64)
             beta_act  = 0.5_wp
             sigma_act = sigma_bulk**( 1.0_wp + beta_act )
!
!--          Calculate mean geometric supersaturation (s_0) with
!--          parameterization by Khvorostyanov and Curry (2006)
             s_0   = dry_aerosol_radius **(- ( 1.0_wp + beta_act ) ) *         &
               ( 4.0_wp * afactor**3 / ( 27.0_wp * bfactor ) )**0.5_wp

!
!--          Calculate number of activated CCN as a function of
!--          supersaturation and taking Koehler theory into account
!--          (see: Khvorostyanov + Curry, 2006, J. Geo. Res., 111)
             n_ccn = ( na_init / 2.0_wp ) * ( 1.0_wp - ERF(                    &
                LOG( s_0 / sat ) / ( SQRT(2.0_wp) * LOG(sigma_act) ) ) )
             activ = MAX( ( n_ccn ) / dt_micro, 0.0_wp )

             nc(k,j,i) = MIN( (nc(k,j,i) + activ * dt_micro * flag), na_init)
          ENDIF

       ENDDO

    END SUBROUTINE activation_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate condensation rate for cloud water content (after Khairoutdinov and
!> Kogan, 2000).
!------------------------------------------------------------------------------!
    SUBROUTINE condensation

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !<
       INTEGER(iwp) ::  j                 !<
       INTEGER(iwp) ::  k                 !<

       REAL(wp)     ::  cond              !<
       REAL(wp)     ::  cond_max          !<
       REAL(wp)     ::  dc                !<
       REAL(wp)     ::  evap              !<
       REAL(wp)     ::  g_fac             !<
       REAL(wp)     ::  nc_0              !<
       REAL(wp)     ::  temp              !<
       REAL(wp)     ::  xc                !<

       REAL(wp) ::  flag                  !< flag to indicate first grid level above

       CALL cpu_log( log_point_s(66), 'condensation', 'start' )

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )
!
!--             Call calculation of supersaturation 
                CALL supersaturation ( i, j, k )
!
!--             Actual temperature:
                IF ( microphysics_seifert ) THEN
                   temp = t_l + lv_d_cp * ( qc(k,j,i) + qr(k,j,i) )
                ELSEIF ( microphysics_morrison_no_rain ) THEN
                   temp = t_l + lv_d_cp * qc(k,j,i) 
                ENDIF  

                g_fac  = 1.0_wp / ( ( l_v / ( r_v * temp ) - 1.0_wp ) *        &
                                    l_v / ( thermal_conductivity_l * temp )    &
                                    + r_v * temp / ( diff_coeff_l * e_s )      &
                                    )
!
!--             Mean weight of cloud drops
                IF ( nc(k,j,i) <= 0.0_wp) CYCLE
                xc = MAX( (hyrho(k) * qc(k,j,i) / nc(k,j,i)), xcmin)
!
!--             Weight averaged diameter of cloud drops:
                dc   = ( xc * dpirho_l )**( 1.0_wp / 3.0_wp )
!
!--             Integral diameter of cloud drops
                nc_0 = nc(k,j,i) * dc
!
!--             Condensation needs only to be calculated in supersaturated regions
                IF ( sat > 0.0_wp )  THEN
!
!--                Condensation rate of cloud water content
!--                after KK scheme.
!--                (see: Khairoutdinov + Kogan, 2000, Mon. Wea. Rev.,128)
                   cond      = 2.0_wp * pi * nc_0 * g_fac * sat / hyrho(k)
                   IF ( microphysics_seifert ) THEN
                      cond_max  = q(k,j,i) - q_s - qc(k,j,i) - qr(k,j,i)
                   ELSEIF ( microphysics_morrison_no_rain ) THEN
                      cond_max  = q(k,j,i) - q_s - qc(k,j,i) 
                   ENDIF
                   cond      = MIN( cond, cond_max / dt_micro )

                   qc(k,j,i) = qc(k,j,i) + cond * dt_micro * flag
                ELSEIF ( sat < 0.0_wp ) THEN
                   evap      = 2.0_wp * pi * nc_0 * g_fac * sat / hyrho(k)
                   evap      = MAX( evap, -qc(k,j,i) / dt_micro )

                   qc(k,j,i) = qc(k,j,i) + evap * dt_micro * flag
                ENDIF
                IF ( nc(k,j,i) * xcmin > qc(k,j,i) * hyrho(k) )  THEN
                   nc(k,j,i) = qc(k,j,i) * hyrho(k) / xcmin
                ENDIF
             ENDDO
          ENDDO
       ENDDO

       CALL cpu_log( log_point_s(66), 'condensation', 'stop' )

    END SUBROUTINE condensation

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate condensation rate for cloud water content (after Khairoutdinov and
!> Kogan, 2000).
!------------------------------------------------------------------------------!
    SUBROUTINE condensation_ij( i, j )

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !<
       INTEGER(iwp) ::  j                 !<
       INTEGER(iwp) ::  k                 !<

       REAL(wp)     ::  cond              !<
       REAL(wp)     ::  cond_max          !<
       REAL(wp)     ::  dc                !<
       REAL(wp)     ::  evap              !<
       REAL(wp)     ::  flag              !< flag to indicate first grid level above surface
       REAL(wp)     ::  g_fac             !<
       REAL(wp)     ::  nc_0              !<
       REAL(wp)     ::  temp              !<
       REAL(wp)     ::  xc                !<


       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )
!
!--       Call calculation of supersaturation
          CALL supersaturation ( i, j, k )
!
!--       Actual temperature:
          IF ( microphysics_seifert ) THEN
             temp = t_l + lv_d_cp * ( qc(k,j,i) + qr(k,j,i) )
          ELSEIF ( microphysics_morrison_no_rain ) THEN
             temp = t_l + lv_d_cp * qc(k,j,i) 
          ENDIF  

          g_fac  = 1.0_wp / ( ( l_v / ( r_v * temp ) - 1.0_wp ) *        &
                              l_v / ( thermal_conductivity_l * temp )    &
                              + r_v * temp / ( diff_coeff_l * e_s )      &
                            )
!
!--       Mean weight of cloud drops
          IF ( nc(k,j,i) <= 0.0_wp) CYCLE
          xc = MAX( (hyrho(k) * qc(k,j,i) / nc(k,j,i)), xcmin)
!
!--       Weight averaged diameter of cloud drops:
          dc   = ( xc * dpirho_l )**( 1.0_wp / 3.0_wp )
!
!--       Integral diameter of cloud drops
          nc_0 = nc(k,j,i) * dc
!
!--       Condensation needs only to be calculated in supersaturated regions
          IF ( sat > 0.0_wp )  THEN
!
!--          Condensation rate of cloud water content
!--          after KK scheme.
!--          (see: Khairoutdinov + Kogan, 2000, Mon. Wea. Rev.,128)
             cond      = 2.0_wp * pi * nc_0 * g_fac * sat / hyrho(k)
             IF ( microphysics_seifert ) THEN
                cond_max  = q(k,j,i) - q_s - qc(k,j,i) - qr(k,j,i)
             ELSEIF ( microphysics_morrison_no_rain ) THEN
                cond_max  = q(k,j,i) - q_s - qc(k,j,i) 
             ENDIF
             cond      = MIN( cond, cond_max / dt_micro )

             qc(k,j,i) = qc(k,j,i) + cond * dt_micro * flag
          ELSEIF ( sat < 0.0_wp ) THEN
             evap      = 2.0_wp * pi * nc_0 * g_fac * sat / hyrho(k)
             evap      = MAX( evap, -qc(k,j,i) / dt_micro )

             qc(k,j,i) = qc(k,j,i) + evap * dt_micro * flag
          ENDIF
       ENDDO

    END SUBROUTINE condensation_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Autoconversion rate (Seifert and Beheng, 2006).
!------------------------------------------------------------------------------!
    SUBROUTINE autoconversion

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !<
       INTEGER(iwp) ::  j                 !<
       INTEGER(iwp) ::  k                 !<

       REAL(wp)     ::  alpha_cc          !<
       REAL(wp)     ::  autocon           !<
       REAL(wp)     ::  dissipation       !<
       REAL(wp)     ::  flag              !< flag to mask topography grid points
       REAL(wp)     ::  k_au              !<
       REAL(wp)     ::  l_mix             !<
       REAL(wp)     ::  nc_auto           !<
       REAL(wp)     ::  nu_c              !<
       REAL(wp)     ::  phi_au            !<
       REAL(wp)     ::  r_cc              !<
       REAL(wp)     ::  rc                !<
       REAL(wp)     ::  re_lambda         !<
       REAL(wp)     ::  sigma_cc          !<
       REAL(wp)     ::  tau_cloud         !<
       REAL(wp)     ::  xc                !<

       CALL cpu_log( log_point_s(47), 'autoconversion', 'start' )

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )

                IF ( microphysics_morrison ) THEN
                   nc_auto = nc(k,j,i)
                ELSE
                   nc_auto = nc_const
                ENDIF

                IF ( qc(k,j,i) > eps_sb  .AND.  nc_auto > eps_mr )  THEN

                   k_au = k_cc / ( 20.0_wp * x0 )
!
!--                Intern time scale of coagulation (Seifert and Beheng, 2006):
!--                (1.0_wp - qc(k,j,i) / ( qc(k,j,i) + qr(k,j,i) ))
                   tau_cloud = MAX( 1.0_wp - qc(k,j,i) / ( qr(k,j,i) +         &
                                    qc(k,j,i) ), 0.0_wp )
!
!--                Universal function for autoconversion process
!--                (Seifert and Beheng, 2006):
                   phi_au = 600.0_wp * tau_cloud**0.68_wp *                    &
                            ( 1.0_wp - tau_cloud**0.68_wp )**3
!
!--                Shape parameter of gamma distribution (Geoffroy et al., 2010):
!--                (Use constant nu_c = 1.0_wp instead?)
                   nu_c = 1.0_wp !MAX( 0.0_wp, 1580.0_wp * hyrho(k) * qc(k,j,i) - 0.28_wp )
!
!--                Mean weight of cloud droplets:
                   xc = MAX( hyrho(k) * qc(k,j,i) / nc_auto, xcmin) 
!
!--                Parameterized turbulence effects on autoconversion (Seifert,
!--                Nuijens and Stevens, 2010)
                   IF ( collision_turbulence )  THEN
!
!--                   Weight averaged radius of cloud droplets:
                      rc = 0.5_wp * ( xc * dpirho_l )**( 1.0_wp / 3.0_wp )

                      alpha_cc = ( a_1 + a_2 * nu_c ) / ( 1.0_wp + a_3 * nu_c )
                      r_cc     = ( b_1 + b_2 * nu_c ) / ( 1.0_wp + b_3 * nu_c )
                      sigma_cc = ( c_1 + c_2 * nu_c ) / ( 1.0_wp + c_3 * nu_c )
!
!--                   Mixing length (neglecting distance to ground and
!--                   stratification)
                      l_mix = ( dx * dy * dzu(k) )**( 1.0_wp / 3.0_wp )
!
!--                   Limit dissipation rate according to Seifert, Nuijens and
!--                   Stevens (2010)
                      dissipation = MIN( 0.06_wp, diss(k,j,i) )
!
!--                   Compute Taylor-microscale Reynolds number:
                      re_lambda = 6.0_wp / 11.0_wp *                           &
                                  ( l_mix / c_const )**( 2.0_wp / 3.0_wp ) *   &
                                  SQRT( 15.0_wp / kin_vis_air ) *              &
                                  dissipation**( 1.0_wp / 6.0_wp )
!
!--                   The factor of 1.0E4 is needed to convert the dissipation
!--                   rate from m2 s-3 to cm2 s-3.
                      k_au = k_au * ( 1.0_wp +                                 &
                                      dissipation * 1.0E4_wp *                 &
                                      ( re_lambda * 1.0E-3_wp )**0.25_wp *     &
                                      ( alpha_cc * EXP( -1.0_wp * ( ( rc -     &
                                                                      r_cc ) / &
                                                        sigma_cc )**2          &
                                                      ) + beta_cc              &
                                      )                                        &
                                    )
                   ENDIF
!
!--                Autoconversion rate (Seifert and Beheng, 2006):
                   autocon = k_au * ( nu_c + 2.0_wp ) * ( nu_c + 4.0_wp ) /    &
                             ( nu_c + 1.0_wp )**2 * qc(k,j,i)**2 * xc**2 *     &
                             ( 1.0_wp + phi_au / ( 1.0_wp - tau_cloud )**2 ) * &
                             rho_surface
                   autocon = MIN( autocon, qc(k,j,i) / dt_micro )

                   qr(k,j,i) = qr(k,j,i) + autocon * dt_micro * flag
                   qc(k,j,i) = qc(k,j,i) - autocon * dt_micro * flag
                   nr(k,j,i) = nr(k,j,i) + autocon / x0 * hyrho(k) * dt_micro  &
                                                              * flag
                   IF ( microphysics_morrison ) THEN
                      nc(k,j,i) = nc(k,j,i) - MIN( nc(k,j,i), 2.0_wp *         &
                                    autocon / x0 * hyrho(k) * dt_micro * flag )
                   ENDIF

                ENDIF

             ENDDO
          ENDDO
       ENDDO

       CALL cpu_log( log_point_s(47), 'autoconversion', 'stop' )

    END SUBROUTINE autoconversion


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Autoconversion rate (Seifert and Beheng, 2006). Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE autoconversion_ij( i, j )

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !<
       INTEGER(iwp) ::  j                 !<
       INTEGER(iwp) ::  k                 !<

       REAL(wp)     ::  alpha_cc          !<
       REAL(wp)     ::  autocon           !<
       REAL(wp)     ::  dissipation       !<
       REAL(wp)     ::  flag              !< flag to indicate first grid level above surface
       REAL(wp)     ::  k_au              !<
       REAL(wp)     ::  l_mix             !<
       REAL(wp)     ::  nc_auto           !<
       REAL(wp)     ::  nu_c              !<
       REAL(wp)     ::  phi_au            !<
       REAL(wp)     ::  r_cc              !<
       REAL(wp)     ::  rc                !<
       REAL(wp)     ::  re_lambda         !<
       REAL(wp)     ::  sigma_cc          !<
       REAL(wp)     ::  tau_cloud         !<
       REAL(wp)     ::  xc                !<

       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )
          IF ( microphysics_morrison ) THEN
             nc_auto = nc(k,j,i)
          ELSE
             nc_auto = nc_const
          ENDIF

          IF ( qc(k,j,i) > eps_sb  .AND.  nc_auto > eps_mr )  THEN

             k_au = k_cc / ( 20.0_wp * x0 )
!
!--          Intern time scale of coagulation (Seifert and Beheng, 2006):
!--          (1.0_wp - qc(k,j,i) / ( qc(k,j,i) + qr(k,j,i) ))
             tau_cloud = MAX( 1.0_wp - qc(k,j,i) / ( qr(k,j,i) + qc(k,j,i) ),     &
                              0.0_wp )
!
!--          Universal function for autoconversion process
!--          (Seifert and Beheng, 2006):
             phi_au = 600.0_wp * tau_cloud**0.68_wp * ( 1.0_wp - tau_cloud**0.68_wp )**3
!
!--          Shape parameter of gamma distribution (Geoffroy et al., 2010):
!--          (Use constant nu_c = 1.0_wp instead?)
             nu_c = 1.0_wp !MAX( 0.0_wp, 1580.0_wp * hyrho(k) * qc(k,j,i) - 0.28_wp )
!
!--          Mean weight of cloud droplets:
             xc = hyrho(k) * qc(k,j,i) / nc_auto
!
!--          Parameterized turbulence effects on autoconversion (Seifert,
!--          Nuijens and Stevens, 2010)
             IF ( collision_turbulence )  THEN
!
!--             Weight averaged radius of cloud droplets:
                rc = 0.5_wp * ( xc * dpirho_l )**( 1.0_wp / 3.0_wp )

                alpha_cc = ( a_1 + a_2 * nu_c ) / ( 1.0_wp + a_3 * nu_c )
                r_cc     = ( b_1 + b_2 * nu_c ) / ( 1.0_wp + b_3 * nu_c )
                sigma_cc = ( c_1 + c_2 * nu_c ) / ( 1.0_wp + c_3 * nu_c )
!
!--             Mixing length (neglecting distance to ground and stratification)
                l_mix = ( dx * dy * dzu(k) )**( 1.0_wp / 3.0_wp )
!
!--             Limit dissipation rate according to Seifert, Nuijens and
!--             Stevens (2010)
                dissipation = MIN( 0.06_wp, diss(k,j,i) )
!
!--             Compute Taylor-microscale Reynolds number:
                re_lambda = 6.0_wp / 11.0_wp *                                 &
                            ( l_mix / c_const )**( 2.0_wp / 3.0_wp ) *         &
                            SQRT( 15.0_wp / kin_vis_air ) *                    &
                            dissipation**( 1.0_wp / 6.0_wp )
!
!--             The factor of 1.0E4 is needed to convert the dissipation rate
!--             from m2 s-3 to cm2 s-3.
                k_au = k_au * ( 1.0_wp +                                       &
                                dissipation * 1.0E4_wp *                       &
                                ( re_lambda * 1.0E-3_wp )**0.25_wp *           &
                                ( alpha_cc * EXP( -1.0_wp * ( ( rc - r_cc ) /  &
                                                  sigma_cc )**2                &
                                                ) + beta_cc                    &
                                )                                              &
                              )
             ENDIF
!
!--          Autoconversion rate (Seifert and Beheng, 2006):
             autocon = k_au * ( nu_c + 2.0_wp ) * ( nu_c + 4.0_wp ) /          &
                       ( nu_c + 1.0_wp )**2 * qc(k,j,i)**2 * xc**2 *            &
                       ( 1.0_wp + phi_au / ( 1.0_wp - tau_cloud )**2 ) *       &
                       rho_surface
             autocon = MIN( autocon, qc(k,j,i) / dt_micro )

             qr(k,j,i) = qr(k,j,i) + autocon * dt_micro                 * flag
             qc(k,j,i) = qc(k,j,i) - autocon * dt_micro                 * flag
             nr(k,j,i) = nr(k,j,i) + autocon / x0 * hyrho(k) * dt_micro * flag
             IF ( microphysics_morrison ) THEN
                nc(k,j,i) = nc(k,j,i) - MIN( nc(k,j,i), 2.0_wp *                &
                                    autocon / x0 * hyrho(k) * dt_micro * flag )
             ENDIF

          ENDIF

       ENDDO

    END SUBROUTINE autoconversion_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Autoconversion process (Kessler, 1969).
!------------------------------------------------------------------------------!
    SUBROUTINE autoconversion_kessler


       IMPLICIT NONE

       INTEGER(iwp) ::  i      !<
       INTEGER(iwp) ::  j      !<
       INTEGER(iwp) ::  k      !<
       INTEGER(iwp) ::  k_wall !< topgraphy top index

       REAL(wp)    ::  dqdt_precip !<
       REAL(wp)    ::  flag        !< flag to mask topography grid points

       DO  i = nxl, nxr
          DO  j = nys, nyn
!
!--          Determine vertical index of topography top
             k_wall = topo_top_ind(j,i,0)
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )

                IF ( qc(k,j,i) > ql_crit )  THEN
                   dqdt_precip = prec_time_const * ( qc(k,j,i) - ql_crit )
                ELSE
                   dqdt_precip = 0.0_wp
                ENDIF

                qc(k,j,i) = qc(k,j,i) - dqdt_precip * dt_micro * flag
                q(k,j,i)  = q(k,j,i)  - dqdt_precip * dt_micro * flag
                pt(k,j,i) = pt(k,j,i) + dqdt_precip * dt_micro * lv_d_cp *     &
                                        d_exner(k)              * flag

!
!--             Compute the rain rate (stored on surface grid point)
                prr(k_wall,j,i) = prr(k_wall,j,i) + dqdt_precip * dzw(k) * flag

             ENDDO
          ENDDO
       ENDDO

   END SUBROUTINE autoconversion_kessler

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Autoconversion process (Kessler, 1969).
!------------------------------------------------------------------------------!
    SUBROUTINE autoconversion_kessler_ij( i, j )


       IMPLICIT NONE

       INTEGER(iwp) ::  i      !<
       INTEGER(iwp) ::  j      !<
       INTEGER(iwp) ::  k      !<
       INTEGER(iwp) ::  k_wall !< topography top index

       REAL(wp)    ::  dqdt_precip !<
       REAL(wp)    ::  flag              !< flag to indicate first grid level above surface

!
!--    Determine vertical index of topography top
       k_wall = topo_top_ind(j,i,0)
       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )

          IF ( qc(k,j,i) > ql_crit )  THEN
             dqdt_precip = prec_time_const * ( qc(k,j,i) - ql_crit )
          ELSE
             dqdt_precip = 0.0_wp
          ENDIF

          qc(k,j,i) = qc(k,j,i) - dqdt_precip * dt_micro * flag
          q(k,j,i)  = q(k,j,i)  - dqdt_precip * dt_micro * flag
          pt(k,j,i) = pt(k,j,i) + dqdt_precip * dt_micro * lv_d_cp * d_exner(k)  &
                                                                 * flag

!
!--       Compute the rain rate (stored on surface grid point)
          prr(k_wall,j,i) = prr(k_wall,j,i) + dqdt_precip * dzw(k) * flag

       ENDDO

    END SUBROUTINE autoconversion_kessler_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Accretion rate (Seifert and Beheng, 2006).
!------------------------------------------------------------------------------!
    SUBROUTINE accretion

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !<
       INTEGER(iwp) ::  j                 !<
       INTEGER(iwp) ::  k                 !<

       REAL(wp)     ::  accr              !<
       REAL(wp)     ::  flag              !< flag to mask topography grid points
       REAL(wp)     ::  k_cr              !<
       REAL(wp)     ::  nc_accr           !<
       REAL(wp)     ::  phi_ac            !<
       REAL(wp)     ::  tau_cloud         !<
       REAL(wp)     ::  xc                !<


       CALL cpu_log( log_point_s(56), 'accretion', 'start' )

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )

                IF ( microphysics_morrison ) THEN
                   nc_accr = nc(k,j,i)
                ELSE
                   nc_accr = nc_const
                ENDIF

                IF ( ( qc(k,j,i) > eps_sb )  .AND.  ( qr(k,j,i) > eps_sb )    &
                                             .AND.  ( nc_accr > eps_mr ) ) THEN
!
!--                Intern time scale of coagulation (Seifert and Beheng, 2006):
                   tau_cloud = 1.0_wp - qc(k,j,i) / ( qc(k,j,i) + qr(k,j,i) )
!
!--                Universal function for accretion process (Seifert and
!--                Beheng, 2001):
                   phi_ac = ( tau_cloud / ( tau_cloud + 5.0E-5_wp ) )**4

!
!--                Mean weight of cloud drops
                   xc = MAX( (hyrho(k) * qc(k,j,i) / nc_accr), xcmin)
!
!--                Parameterized turbulence effects on autoconversion (Seifert,
!--                Nuijens and Stevens, 2010). The factor of 1.0E4 is needed to
!--                convert the dissipation rate (diss) from m2 s-3 to cm2 s-3.
                   IF ( collision_turbulence )  THEN
                      k_cr = k_cr0 * ( 1.0_wp + 0.05_wp *                      &
                                       MIN( 600.0_wp,                          &
                                            diss(k,j,i) * 1.0E4_wp )**0.25_wp  &
                                     )
                   ELSE
                      k_cr = k_cr0
                   ENDIF
!
!--                Accretion rate (Seifert and Beheng, 2006):
                   accr = k_cr * qc(k,j,i) * qr(k,j,i) * phi_ac *              &
                          SQRT( rho_surface * hyrho(k) )
                   accr = MIN( accr, qc(k,j,i) / dt_micro )

                   qr(k,j,i) = qr(k,j,i) + accr * dt_micro * flag
                   qc(k,j,i) = qc(k,j,i) - accr * dt_micro * flag
                   IF ( microphysics_morrison )  THEN
                      nc(k,j,i) = nc(k,j,i) - MIN( nc(k,j,i),                  &
                                         accr / xc * hyrho(k) * dt_micro * flag)
                   ENDIF

                ENDIF

             ENDDO
          ENDDO
       ENDDO

       CALL cpu_log( log_point_s(56), 'accretion', 'stop' )

    END SUBROUTINE accretion

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Accretion rate (Seifert and Beheng, 2006). Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE accretion_ij( i, j )

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !<
       INTEGER(iwp) ::  j                 !<
       INTEGER(iwp) ::  k                 !<

       REAL(wp)     ::  accr              !<
       REAL(wp)     ::  flag              !< flag to indicate first grid level above surface
       REAL(wp)     ::  k_cr              !<
       REAL(wp)     ::  nc_accr           !<
       REAL(wp)     ::  phi_ac            !<
       REAL(wp)     ::  tau_cloud         !<
       REAL(wp)     ::  xc                !<


       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )
          IF ( microphysics_morrison ) THEN
             nc_accr = nc(k,j,i)
          ELSE
             nc_accr = nc_const
          ENDIF

          IF ( ( qc(k,j,i) > eps_sb )  .AND.  ( qr(k,j,i) > eps_sb )  .AND.  &
               ( nc_accr > eps_mr ) )  THEN
!
!--          Intern time scale of coagulation (Seifert and Beheng, 2006):
             tau_cloud = 1.0_wp - qc(k,j,i) / ( qc(k,j,i) + qr(k,j,i) )
!
!--          Universal function for accretion process
!--          (Seifert and Beheng, 2001):
             phi_ac = ( tau_cloud / ( tau_cloud + 5.0E-5_wp ) )**4

!
!--          Mean weight of cloud drops
             xc = MAX( (hyrho(k) * qc(k,j,i) / nc_accr), xcmin)
!
!--          Parameterized turbulence effects on autoconversion (Seifert,
!--          Nuijens and Stevens, 2010). The factor of 1.0E4 is needed to
!--          convert the dissipation rate (diss) from m2 s-3 to cm2 s-3.
             IF ( collision_turbulence )  THEN
                k_cr = k_cr0 * ( 1.0_wp + 0.05_wp *                            &
                                 MIN( 600.0_wp,                                &
                                      diss(k,j,i) * 1.0E4_wp )**0.25_wp        &
                               )
             ELSE
                k_cr = k_cr0
             ENDIF
!
!--          Accretion rate (Seifert and Beheng, 2006):
             accr = k_cr * qc(k,j,i) * qr(k,j,i) * phi_ac *                      &
                    SQRT( rho_surface * hyrho(k) )
             accr = MIN( accr, qc(k,j,i) / dt_micro )

             qr(k,j,i) = qr(k,j,i) + accr * dt_micro * flag
             qc(k,j,i) = qc(k,j,i) - accr * dt_micro * flag
             IF ( microphysics_morrison )  THEN
                nc(k,j,i) = nc(k,j,i) - MIN( nc(k,j,i), accr / xc *               &
                                             hyrho(k) * dt_micro * flag        &
                                         )
             ENDIF


          ENDIF

       ENDDO

    END SUBROUTINE accretion_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Collisional breakup rate (Seifert, 2008).
!------------------------------------------------------------------------------!
    SUBROUTINE selfcollection_breakup

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !<
       INTEGER(iwp) ::  j                 !<
       INTEGER(iwp) ::  k                 !<

       REAL(wp)     ::  breakup           !<
       REAL(wp)     ::  dr                !<
       REAL(wp)     ::  flag              !< flag to mask topography grid points
       REAL(wp)     ::  phi_br            !<
       REAL(wp)     ::  selfcoll          !<

       CALL cpu_log( log_point_s(57), 'selfcollection', 'start' )

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )

                IF ( qr(k,j,i) > eps_sb )  THEN
!
!--                Selfcollection rate (Seifert and Beheng, 2001):
                   selfcoll = k_rr * nr(k,j,i) * qr(k,j,i) *                   &
                              SQRT( hyrho(k) * rho_surface )
!
!--                Weight averaged diameter of rain drops:
                   dr = ( hyrho(k) * qr(k,j,i) /                               &
                          nr(k,j,i) * dpirho_l )**( 1.0_wp / 3.0_wp )
!
!--                Collisional breakup rate (Seifert, 2008):
                   IF ( dr >= 0.3E-3_wp )  THEN
                      phi_br  = k_br * ( dr - 1.1E-3_wp )
                      breakup = selfcoll * ( phi_br + 1.0_wp )
                   ELSE
                      breakup = 0.0_wp
                   ENDIF

                   selfcoll = MAX( breakup - selfcoll, -nr(k,j,i) / dt_micro )
                   nr(k,j,i) = nr(k,j,i) + selfcoll * dt_micro * flag

                ENDIF
             ENDDO
          ENDDO
       ENDDO

       CALL cpu_log( log_point_s(57), 'selfcollection', 'stop' )

    END SUBROUTINE selfcollection_breakup


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Collisional breakup rate (Seifert, 2008). Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE selfcollection_breakup_ij( i, j )

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !<
       INTEGER(iwp) ::  j                 !<
       INTEGER(iwp) ::  k                 !<

       REAL(wp)     ::  breakup           !<
       REAL(wp)     ::  dr                !<
       REAL(wp)     ::  flag              !< flag to indicate first grid level above surface
       REAL(wp)     ::  phi_br            !<
       REAL(wp)     ::  selfcoll          !<

       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )

          IF ( qr(k,j,i) > eps_sb )  THEN
!
!--          Selfcollection rate (Seifert and Beheng, 2001):
             selfcoll = k_rr * nr(k,j,i) * qr(k,j,i) * SQRT( hyrho(k) * rho_surface )
!
!--          Weight averaged diameter of rain drops:
             dr = ( hyrho(k) * qr(k,j,i) / nr(k,j,i) * dpirho_l )**( 1.0_wp / 3.0_wp )
!
!--          Collisional breakup rate (Seifert, 2008):
             IF ( dr >= 0.3E-3_wp )  THEN
                phi_br  = k_br * ( dr - 1.1E-3_wp )
                breakup = selfcoll * ( phi_br + 1.0_wp )
             ELSE
                breakup = 0.0_wp
             ENDIF

             selfcoll = MAX( breakup - selfcoll, -nr(k,j,i) / dt_micro )
             nr(k,j,i) = nr(k,j,i) + selfcoll * dt_micro * flag

          ENDIF
       ENDDO

    END SUBROUTINE selfcollection_breakup_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Evaporation of precipitable water. Condensation is neglected for
!> precipitable water.
!------------------------------------------------------------------------------!
    SUBROUTINE evaporation_rain

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !<
       INTEGER(iwp) ::  j                 !<
       INTEGER(iwp) ::  k                 !<

       REAL(wp)     ::  dr                !<
       REAL(wp)     ::  evap              !<
       REAL(wp)     ::  evap_nr           !<
       REAL(wp)     ::  f_vent            !<
       REAL(wp)     ::  flag              !< flag to mask topography grid points
       REAL(wp)     ::  g_evap            !<
       REAL(wp)     ::  lambda_r          !<
       REAL(wp)     ::  mu_r              !<
       REAL(wp)     ::  mu_r_2            !<
       REAL(wp)     ::  mu_r_5d2          !<
       REAL(wp)     ::  nr_0              !<
       REAL(wp)     ::  temp              !<
       REAL(wp)     ::  xr                !<

       CALL cpu_log( log_point_s(58), 'evaporation', 'start' )

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )

                IF ( qr(k,j,i) > eps_sb )  THEN

!
!--                Call calculation of supersaturation  
                   CALL supersaturation ( i, j, k )
!
!--                Evaporation needs only to be calculated in subsaturated regions
                   IF ( sat < 0.0_wp )  THEN
!
!--                   Actual temperature:
                      temp = t_l + lv_d_cp * ( qc(k,j,i) + qr(k,j,i) )

                      g_evap = 1.0_wp / ( ( l_v / ( r_v * temp ) - 1.0_wp ) *     &
                                          l_v / ( thermal_conductivity_l * temp ) &
                                          + r_v * temp / ( diff_coeff_l * e_s )   &
                                        )
!
!--                   Mean weight of rain drops
                      xr = hyrho(k) * qr(k,j,i) / nr(k,j,i)
!
!--                   Weight averaged diameter of rain drops:
                      dr = ( xr * dpirho_l )**( 1.0_wp / 3.0_wp )
!
!--                   Compute ventilation factor and intercept parameter
!--                   (Seifert and Beheng, 2006; Seifert, 2008):
                      IF ( ventilation_effect )  THEN
!
!--                      Shape parameter of gamma distribution (Milbrandt and Yau,
!--                      2005; Stevens and Seifert, 2008):
                         mu_r = 10.0_wp * ( 1.0_wp + TANH( 1.2E3_wp *          &
                                                          ( dr - 1.4E-3_wp ) ) )
!
!--                      Slope parameter of gamma distribution (Seifert, 2008):
                         lambda_r = ( ( mu_r + 3.0_wp ) * ( mu_r + 2.0_wp ) *  &
                                      ( mu_r + 1.0_wp )                        &
                                    )**( 1.0_wp / 3.0_wp ) / dr

                         mu_r_2   = mu_r + 2.0_wp
                         mu_r_5d2 = mu_r + 2.5_wp

                         f_vent = a_vent * gamm( mu_r_2 ) *                     &
                                  lambda_r**( -mu_r_2 ) + b_vent *              &
                                  schmidt_p_1d3 * SQRT( a_term / kin_vis_air ) *&
                                  gamm( mu_r_5d2 ) * lambda_r**( -mu_r_5d2 ) *  &
                                  ( 1.0_wp -                                    &
                                    0.5_wp * ( b_term / a_term ) *              &
                                    ( lambda_r / ( c_term + lambda_r )          &
                                    )**mu_r_5d2 -                               &
                                    0.125_wp * ( b_term / a_term )**2 *         &
                                    ( lambda_r / ( 2.0_wp * c_term + lambda_r ) &
                                    )**mu_r_5d2 -                               &
                                    0.0625_wp * ( b_term / a_term )**3 *        &
                                    ( lambda_r / ( 3.0_wp * c_term + lambda_r ) &
                                    )**mu_r_5d2 -                               &
                                    0.0390625_wp * ( b_term / a_term )**4 *     &
                                    ( lambda_r / ( 4.0_wp * c_term + lambda_r ) &
                                    )**mu_r_5d2                                 &
                                  )

                         nr_0   = nr(k,j,i) * lambda_r**( mu_r + 1.0_wp ) /    &
                                  gamm( mu_r + 1.0_wp )
                      ELSE
                         f_vent = 1.0_wp
                         nr_0   = nr(k,j,i) * dr
                      ENDIF
   !
   !--                Evaporation rate of rain water content (Seifert and
   !--                Beheng, 2006):
                      evap    = 2.0_wp * pi * nr_0 * g_evap * f_vent * sat /   &
                                hyrho(k)
                      evap    = MAX( evap, -qr(k,j,i) / dt_micro )
                      evap_nr = MAX( c_evap * evap / xr * hyrho(k),            &
                                     -nr(k,j,i) / dt_micro )

                      qr(k,j,i) = qr(k,j,i) + evap    * dt_micro * flag
                      nr(k,j,i) = nr(k,j,i) + evap_nr * dt_micro * flag

                   ENDIF
                ENDIF

             ENDDO
          ENDDO
       ENDDO

       CALL cpu_log( log_point_s(58), 'evaporation', 'stop' )

    END SUBROUTINE evaporation_rain


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Evaporation of precipitable water. Condensation is neglected for
!> precipitable water. Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE evaporation_rain_ij( i, j )

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !<
       INTEGER(iwp) ::  j                 !<
       INTEGER(iwp) ::  k                 !<

       REAL(wp)     ::  dr                !<
       REAL(wp)     ::  evap              !<
       REAL(wp)     ::  evap_nr           !<
       REAL(wp)     ::  f_vent            !<
       REAL(wp)     ::  flag              !< flag to indicate first grid level above surface
       REAL(wp)     ::  g_evap            !<
       REAL(wp)     ::  lambda_r          !<
       REAL(wp)     ::  mu_r              !<
       REAL(wp)     ::  mu_r_2            !<
       REAL(wp)     ::  mu_r_5d2          !<
       REAL(wp)     ::  nr_0              !<
       REAL(wp)     ::  temp              !<
       REAL(wp)     ::  xr                !<

       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )

          IF ( qr(k,j,i) > eps_sb )  THEN
!
!--          Call calculation of supersaturation
             CALL supersaturation ( i, j, k )
!
!--          Evaporation needs only to be calculated in subsaturated regions
             IF ( sat < 0.0_wp )  THEN
!
!--             Actual temperature:
                temp = t_l + lv_d_cp * ( qc(k,j,i) + qr(k,j,i) )

                g_evap = 1.0_wp / ( ( l_v / ( r_v * temp ) - 1.0_wp ) * l_v /  &
                                    ( thermal_conductivity_l * temp ) +        &
                                    r_v * temp / ( diff_coeff_l * e_s )        &
                                  )
!
!--             Mean weight of rain drops
                xr = hyrho(k) * qr(k,j,i) / nr(k,j,i)
!
!--             Weight averaged diameter of rain drops:
                dr = ( xr * dpirho_l )**( 1.0_wp / 3.0_wp )
!
!--             Compute ventilation factor and intercept parameter
!--             (Seifert and Beheng, 2006; Seifert, 2008):
                IF ( ventilation_effect )  THEN
!
!--                Shape parameter of gamma distribution (Milbrandt and Yau, 2005;
!--                Stevens and Seifert, 2008):
                   mu_r = 10.0_wp * ( 1.0_wp + TANH( 1.2E3_wp * ( dr - 1.4E-3_wp ) ) )
!
!--                Slope parameter of gamma distribution (Seifert, 2008):
                   lambda_r = ( ( mu_r + 3.0_wp ) * ( mu_r + 2.0_wp ) *        &
                                ( mu_r + 1.0_wp )                              &
                              )**( 1.0_wp / 3.0_wp ) / dr

                   mu_r_2   = mu_r + 2.0_wp
                   mu_r_5d2 = mu_r + 2.5_wp

                   f_vent = a_vent * gamm( mu_r_2 ) * lambda_r**( -mu_r_2 ) +  &
                            b_vent * schmidt_p_1d3 *                           &
                            SQRT( a_term / kin_vis_air ) * gamm( mu_r_5d2 ) *  &
                            lambda_r**( -mu_r_5d2 ) *                          &
                            ( 1.0_wp -                                         &
                              0.5_wp * ( b_term / a_term ) *                   &
                              ( lambda_r / ( c_term + lambda_r )               &
                              )**mu_r_5d2 -                                    &
                              0.125_wp * ( b_term / a_term )**2 *              &
                              ( lambda_r / ( 2.0_wp * c_term + lambda_r )      &
                              )**mu_r_5d2 -                                    &
                              0.0625_wp * ( b_term / a_term )**3 *             &
                              ( lambda_r / ( 3.0_wp * c_term + lambda_r )      &
                              )**mu_r_5d2 -                                    &
                              0.0390625_wp * ( b_term / a_term )**4 *          &
                              ( lambda_r / ( 4.0_wp * c_term + lambda_r )      &
                              )**mu_r_5d2                                      &
                            )

                   nr_0   = nr(k,j,i) * lambda_r**( mu_r + 1.0_wp ) /           &
                            gamm( mu_r + 1.0_wp )
                ELSE
                   f_vent = 1.0_wp
                   nr_0   = nr(k,j,i) * dr
                ENDIF
!
!--             Evaporation rate of rain water content (Seifert and Beheng, 2006):
                evap    = 2.0_wp * pi * nr_0 * g_evap * f_vent * sat / hyrho(k)
                evap    = MAX( evap, -qr(k,j,i) / dt_micro )
                evap_nr = MAX( c_evap * evap / xr * hyrho(k),                  &
                               -nr(k,j,i) / dt_micro )

                qr(k,j,i) = qr(k,j,i) + evap    * dt_micro * flag
                nr(k,j,i) = nr(k,j,i) + evap_nr * dt_micro * flag

             ENDIF
          ENDIF

       ENDDO

    END SUBROUTINE evaporation_rain_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Sedimentation of cloud droplets (Ackermann et al., 2009, MWR).
!------------------------------------------------------------------------------!
    SUBROUTINE sedimentation_cloud


       IMPLICIT NONE

       INTEGER(iwp) ::  i             !<
       INTEGER(iwp) ::  j             !<
       INTEGER(iwp) ::  k             !<

       REAL(wp) ::  flag              !< flag to mask topography grid points
       REAL(wp) ::  nc_sedi           !<

       REAL(wp), DIMENSION(nzb:nzt+1) ::  sed_qc !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  sed_nc !<


       CALL cpu_log( log_point_s(59), 'sed_cloud', 'start' )

       sed_qc(nzt+1) = 0.0_wp
       sed_nc(nzt+1) = 0.0_wp

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzt, nzb+1, -1
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )

                IF ( microphysics_morrison ) THEN
                   nc_sedi = nc(k,j,i)
                ELSE
                   nc_sedi = nc_const
                ENDIF

!
!--             Sedimentation fluxes for number concentration are only calculated
!--             for cloud_scheme = 'morrison'
                IF ( microphysics_morrison ) THEN
                   IF ( qc(k,j,i) > eps_sb  .AND.  nc(k,j,i) > eps_mr )  THEN
                      sed_nc(k) = sed_qc_const *                               &
                              ( qc(k,j,i) * hyrho(k) )**( 2.0_wp / 3.0_wp ) *  &
                              ( nc(k,j,i) )**( 1.0_wp / 3.0_wp )
                   ELSE
                      sed_nc(k) = 0.0_wp
                   ENDIF

                   sed_nc(k) = MIN( sed_nc(k), hyrho(k) * dzu(k+1) *           &
                                    nc(k,j,i) / dt_micro + sed_nc(k+1)         &
                                  ) * flag

                   nc(k,j,i) = nc(k,j,i) + ( sed_nc(k+1) - sed_nc(k) ) *       &
                                         ddzu(k+1) / hyrho(k) * dt_micro * flag
                ENDIF

                IF ( qc(k,j,i) > eps_sb .AND.  nc_sedi > eps_mr )  THEN
                   sed_qc(k) = sed_qc_const * nc_sedi**( -2.0_wp / 3.0_wp ) *  &
                               ( qc(k,j,i) * hyrho(k) )**( 5.0_wp / 3.0_wp ) * &
                                                                           flag
                ELSE
                   sed_qc(k) = 0.0_wp
                ENDIF

                sed_qc(k) = MIN( sed_qc(k), hyrho(k) * dzu(k+1) * q(k,j,i) /   &
                                            dt_micro + sed_qc(k+1)             &
                               ) * flag

                q(k,j,i)  = q(k,j,i)  + ( sed_qc(k+1) - sed_qc(k) ) *          &
                                        ddzu(k+1) / hyrho(k) * dt_micro * flag
                qc(k,j,i) = qc(k,j,i) + ( sed_qc(k+1) - sed_qc(k) ) *          &
                                        ddzu(k+1) / hyrho(k) * dt_micro * flag
                pt(k,j,i) = pt(k,j,i) - ( sed_qc(k+1) - sed_qc(k) ) *          &
                                        ddzu(k+1) / hyrho(k) * lv_d_cp *       &
                                        d_exner(k) * dt_micro            * flag

!
!--             Compute the precipitation rate due to cloud (fog) droplets
                IF ( call_microphysics_at_all_substeps )  THEN
                   prr(k,j,i) = prr(k,j,i) +  sed_qc(k) / hyrho(k)             &
                                * weight_substep(intermediate_timestep_count)  &
                                * flag
                ELSE
                   prr(k,j,i) = prr(k,j,i) + sed_qc(k) / hyrho(k) * flag
                ENDIF

             ENDDO
          ENDDO
       ENDDO

       CALL cpu_log( log_point_s(59), 'sed_cloud', 'stop' )

    END SUBROUTINE sedimentation_cloud


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Sedimentation of cloud droplets (Ackermann et al., 2009, MWR).
!> Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE sedimentation_cloud_ij( i, j )

       IMPLICIT NONE

       INTEGER(iwp) ::  i             !<
       INTEGER(iwp) ::  j             !<
       INTEGER(iwp) ::  k             !<

       REAL(wp)     ::  flag    !< flag to indicate first grid level above surface
       REAL(wp)     ::  nc_sedi !<

       REAL(wp), DIMENSION(nzb:nzt+1) ::  sed_nc  !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  sed_qc  !<

       sed_qc(nzt+1) = 0.0_wp
       sed_nc(nzt+1) = 0.0_wp


       DO  k = nzt, nzb+1, -1
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )
          IF ( microphysics_morrison ) THEN
             nc_sedi = nc(k,j,i)
          ELSE
             nc_sedi = nc_const
          ENDIF
!
!--       Sedimentation fluxes for number concentration are only calculated
!--       for cloud_scheme = 'morrison'
          IF ( microphysics_morrison ) THEN
             IF ( qc(k,j,i) > eps_sb  .AND.  nc(k,j,i) > eps_mr )  THEN
                sed_nc(k) = sed_qc_const *                                     &
                            ( qc(k,j,i) * hyrho(k) )**( 2.0_wp / 3.0_wp ) *     &
                            ( nc(k,j,i) )**( 1.0_wp / 3.0_wp )
             ELSE
                sed_nc(k) = 0.0_wp
             ENDIF

             sed_nc(k) = MIN( sed_nc(k), hyrho(k) * dzu(k+1) *                 &
                              nc(k,j,i) / dt_micro + sed_nc(k+1)                &
                            ) * flag

             nc(k,j,i) = nc(k,j,i) + ( sed_nc(k+1) - sed_nc(k) ) *               &
                                         ddzu(k+1) / hyrho(k) * dt_micro * flag
          ENDIF

          IF ( qc(k,j,i) > eps_sb  .AND.  nc_sedi > eps_mr )  THEN
             sed_qc(k) = sed_qc_const * nc_sedi**( -2.0_wp / 3.0_wp ) *        &
                         ( qc(k,j,i) * hyrho(k) )**( 5.0_wp / 3.0_wp )  * flag
          ELSE
             sed_qc(k) = 0.0_wp
          ENDIF

          sed_qc(k) = MIN( sed_qc(k), hyrho(k) * dzu(k+1) * q(k,j,i) /          &
                                      dt_micro + sed_qc(k+1)                   &
                         ) * flag

          q(k,j,i)  = q(k,j,i)  + ( sed_qc(k+1) - sed_qc(k) ) * ddzu(k+1) /      &
                                hyrho(k) * dt_micro * flag
          qc(k,j,i) = qc(k,j,i) + ( sed_qc(k+1) - sed_qc(k) ) * ddzu(k+1) /      &
                                hyrho(k) * dt_micro * flag
          pt(k,j,i) = pt(k,j,i) - ( sed_qc(k+1) - sed_qc(k) ) * ddzu(k+1) /      &
                                hyrho(k) * lv_d_cp * d_exner(k) * dt_micro     &
                                                                * flag

!
!--       Compute the precipitation rate of cloud (fog) droplets
          IF ( call_microphysics_at_all_substeps )  THEN
             prr(k,j,i) = prr(k,j,i) + sed_qc(k) / hyrho(k) *                  &
                              weight_substep(intermediate_timestep_count) * flag
          ELSE
             prr(k,j,i) = prr(k,j,i) + sed_qc(k) / hyrho(k) * flag
          ENDIF

       ENDDO

    END SUBROUTINE sedimentation_cloud_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Computation of sedimentation flux. Implementation according to Stevens
!> and Seifert (2008). Code is based on UCLA-LES.
!------------------------------------------------------------------------------!
    SUBROUTINE sedimentation_rain

       IMPLICIT NONE

       INTEGER(iwp) ::  i             !< running index x direction
       INTEGER(iwp) ::  j             !< running index y direction
       INTEGER(iwp) ::  k             !< running index z direction
       INTEGER(iwp) ::  k_run         !<
       INTEGER(iwp) ::  m             !< running index surface elements
       INTEGER(iwp) ::  surf_e        !< End index of surface elements at (j,i)-gridpoint
       INTEGER(iwp) ::  surf_s        !< Start index of surface elements at (j,i)-gridpoint

       REAL(wp)     ::  c_run                      !<
       REAL(wp)     ::  d_max                      !<
       REAL(wp)     ::  d_mean                     !<
       REAL(wp)     ::  d_min                      !<
       REAL(wp)     ::  dr                         !<
       REAL(wp)     ::  flux                       !<
       REAL(wp)     ::  flag                       !< flag to mask topography grid points
       REAL(wp)     ::  lambda_r                   !<
       REAL(wp)     ::  mu_r                       !<
       REAL(wp)     ::  z_run                      !<

       REAL(wp), DIMENSION(nzb:nzt+1) ::  c_nr     !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  c_qr     !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  nr_slope !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  qr_slope !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  sed_nr   !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  sed_qr   !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  w_nr     !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  w_qr     !<

       CALL cpu_log( log_point_s(60), 'sed_rain', 'start' )

!
!--    Compute velocities
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )

                IF ( qr(k,j,i) > eps_sb )  THEN
!
!--                Weight averaged diameter of rain drops:
                   dr = ( hyrho(k) * qr(k,j,i) /                               &
                          nr(k,j,i) * dpirho_l )**( 1.0_wp / 3.0_wp )
!
!--                Shape parameter of gamma distribution (Milbrandt and Yau, 2005;
!--                Stevens and Seifert, 2008):
                   mu_r = 10.0_wp * ( 1.0_wp + TANH( 1.2E3_wp *                &
                                                     ( dr - 1.4E-3_wp ) ) )
!
!--                Slope parameter of gamma distribution (Seifert, 2008):
                   lambda_r = ( ( mu_r + 3.0_wp ) * ( mu_r + 2.0_wp ) *        &
                                ( mu_r + 1.0_wp ) )**( 1.0_wp / 3.0_wp ) / dr

                   w_nr(k) = MAX( 0.1_wp, MIN( 20.0_wp,                        &
                                               a_term - b_term * ( 1.0_wp +    &
                                                  c_term /                     &
                                                  lambda_r )**( -1.0_wp *      &
                                                  ( mu_r + 1.0_wp ) )          &
                                              )                                &
                                ) * flag

                   w_qr(k) = MAX( 0.1_wp, MIN( 20.0_wp,                        &
                                               a_term - b_term * ( 1.0_wp +    &
                                                  c_term /                     &
                                                  lambda_r )**( -1.0_wp *      &
                                                  ( mu_r + 4.0_wp ) )          &
                                             )                                 &
                                ) * flag
                ELSE
                   w_nr(k) = 0.0_wp
                   w_qr(k) = 0.0_wp
                ENDIF
             ENDDO
!
!--          Adjust boundary values using surface data type.
!--          Upward-facing
             surf_s = bc_h(0)%start_index(j,i)
             surf_e = bc_h(0)%end_index(j,i)
             DO  m = surf_s, surf_e
                k         = bc_h(0)%k(m)
                w_nr(k-1) = w_nr(k)
                w_qr(k-1) = w_qr(k)
             ENDDO
!
!--          Downward-facing
             surf_s = bc_h(1)%start_index(j,i)
             surf_e = bc_h(1)%end_index(j,i)
             DO  m = surf_s, surf_e
                k         = bc_h(1)%k(m)
                w_nr(k+1) = w_nr(k)
                w_qr(k+1) = w_qr(k)
             ENDDO
!
!--          Model top boundary value
             w_nr(nzt+1) = 0.0_wp
             w_qr(nzt+1) = 0.0_wp
!
!--          Compute Courant number
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )

                c_nr(k) = 0.25_wp * ( w_nr(k-1) +                              &
                                      2.0_wp * w_nr(k) + w_nr(k+1) ) *         &
                          dt_micro * ddzu(k) * flag
                c_qr(k) = 0.25_wp * ( w_qr(k-1) +                              &
                                      2.0_wp * w_qr(k) + w_qr(k+1) ) *         &
                          dt_micro * ddzu(k) * flag
             ENDDO
!
!--          Limit slopes with monotonized centered (MC) limiter (van Leer, 1977):
             IF ( limiter_sedimentation )  THEN

                DO k = nzb+1, nzt
!
!--                Predetermine flag to mask topography
                   flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )

                   d_mean = 0.5_wp * ( qr(k+1,j,i) - qr(k-1,j,i) )
                   d_min  = qr(k,j,i) - MIN( qr(k+1,j,i), qr(k,j,i), qr(k-1,j,i) )
                   d_max  = MAX( qr(k+1,j,i), qr(k,j,i), qr(k-1,j,i) ) - qr(k,j,i)

                   qr_slope(k) = SIGN(1.0_wp, d_mean) * MIN ( 2.0_wp * d_min,  &
                                                              2.0_wp * d_max,  &
                                                              ABS( d_mean ) )  &
                                                      * flag

                   d_mean = 0.5_wp * ( nr(k+1,j,i) - nr(k-1,j,i) )
                   d_min  = nr(k,j,i) - MIN( nr(k+1,j,i), nr(k,j,i), nr(k-1,j,i) )
                   d_max  = MAX( nr(k+1,j,i), nr(k,j,i), nr(k-1,j,i) ) - nr(k,j,i)

                   nr_slope(k) = SIGN(1.0_wp, d_mean) * MIN ( 2.0_wp * d_min,  &
                                                              2.0_wp * d_max,  &
                                                              ABS( d_mean ) )
                ENDDO

             ELSE

                nr_slope = 0.0_wp
                qr_slope = 0.0_wp

             ENDIF

             sed_nr(nzt+1) = 0.0_wp
             sed_qr(nzt+1) = 0.0_wp
!
!--          Compute sedimentation flux
             DO  k = nzt, nzb+1, -1
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )
!
!--             Sum up all rain drop number densities which contribute to the flux
!--             through k-1/2
                flux  = 0.0_wp
                z_run = 0.0_wp ! height above z(k)
                k_run = k
                c_run = MIN( 1.0_wp, c_nr(k) )
                DO WHILE ( c_run > 0.0_wp  .AND.  k_run <= nzt )
                   flux  = flux + hyrho(k_run) *                               &
                           ( nr(k_run,j,i) + nr_slope(k_run) *                 &
                           ( 1.0_wp - c_run ) * 0.5_wp ) * c_run * dzu(k_run)  &
                                              * flag
                   z_run = z_run + dzu(k_run) * flag
                   k_run = k_run + 1          * flag
                   c_run = MIN( 1.0_wp, c_nr(k_run) - z_run * ddzu(k_run) )    &
                                              * flag
                ENDDO
!
!--             It is not allowed to sediment more rain drop number density than
!--             available
                flux = MIN( flux,                                              &
                            hyrho(k) * dzu(k+1) * nr(k,j,i) + sed_nr(k+1) *    &
                            dt_micro                                           &
                          )

                sed_nr(k) = flux / dt_micro * flag
                nr(k,j,i) = nr(k,j,i) + ( sed_nr(k+1) - sed_nr(k) ) *          &
                                        ddzu(k+1) / hyrho(k) * dt_micro * flag
!
!--             Sum up all rain water content which contributes to the flux
!--             through k-1/2
                flux  = 0.0_wp
                z_run = 0.0_wp ! height above z(k)
                k_run = k
                c_run = MIN( 1.0_wp, c_qr(k) )

                DO WHILE ( c_run > 0.0_wp  .AND.  k_run <= nzt )

                   flux  = flux + hyrho(k_run) * ( qr(k_run,j,i) +             &
                                  qr_slope(k_run) * ( 1.0_wp - c_run ) *       &
                                  0.5_wp ) * c_run * dzu(k_run) * flag
                   z_run = z_run + dzu(k_run)                   * flag
                   k_run = k_run + 1                            * flag
                   c_run = MIN( 1.0_wp, c_qr(k_run) - z_run * ddzu(k_run) )    &
                                                                * flag

                ENDDO
!
!--             It is not allowed to sediment more rain water content than
!--             available
                flux = MIN( flux,                                              &
                            hyrho(k) * dzu(k) * qr(k,j,i) + sed_qr(k+1) *      &
                            dt_micro                                           &
                          )

                sed_qr(k) = flux / dt_micro * flag

                qr(k,j,i) = qr(k,j,i) + ( sed_qr(k+1) - sed_qr(k) ) *          &
                                        ddzu(k+1) / hyrho(k) * dt_micro * flag
                q(k,j,i)  = q(k,j,i)  + ( sed_qr(k+1) - sed_qr(k) ) *          &
                                        ddzu(k+1) / hyrho(k) * dt_micro * flag
                pt(k,j,i) = pt(k,j,i) - ( sed_qr(k+1) - sed_qr(k) ) *          &
                                        ddzu(k+1) / hyrho(k) * lv_d_cp *       &
                                        d_exner(k) * dt_micro            * flag
!
!--             Compute the rain rate
                IF ( call_microphysics_at_all_substeps )  THEN
                   prr(k,j,i) = prr(k,j,i) + sed_qr(k) / hyrho(k)             &
                                * weight_substep(intermediate_timestep_count) &
                                * flag
                ELSE
                   prr(k,j,i) = prr(k,j,i) + sed_qr(k) / hyrho(k) * flag
                ENDIF

             ENDDO
          ENDDO
       ENDDO

       CALL cpu_log( log_point_s(60), 'sed_rain', 'stop' )

    END SUBROUTINE sedimentation_rain


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Computation of sedimentation flux. Implementation according to Stevens
!> and Seifert (2008). Code is based on UCLA-LES. Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE sedimentation_rain_ij( i, j )

       IMPLICIT NONE

       INTEGER(iwp) ::  i             !< running index x direction
       INTEGER(iwp) ::  j             !< running index y direction
       INTEGER(iwp) ::  k             !< running index z direction
       INTEGER(iwp) ::  k_run         !<
       INTEGER(iwp) ::  m             !< running index surface elements
       INTEGER(iwp) ::  surf_e        !< End index of surface elements at (j,i)-gridpoint
       INTEGER(iwp) ::  surf_s        !< Start index of surface elements at (j,i)-gridpoint

       REAL(wp)     ::  c_run                      !<
       REAL(wp)     ::  d_max                      !<
       REAL(wp)     ::  d_mean                     !<
       REAL(wp)     ::  d_min                      !<
       REAL(wp)     ::  dr                         !<
       REAL(wp)     ::  flux                       !<
       REAL(wp)     ::  flag                       !< flag to indicate first grid level above surface
       REAL(wp)     ::  lambda_r                   !<
       REAL(wp)     ::  mu_r                       !<
       REAL(wp)     ::  z_run                      !<

       REAL(wp), DIMENSION(nzb:nzt+1) ::  c_nr     !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  c_qr     !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  nr_slope !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  qr_slope !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  sed_nr   !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  sed_qr   !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  w_nr     !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  w_qr     !<

!
!--    Compute velocities
       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )

          IF ( qr(k,j,i) > eps_sb )  THEN
!
!--          Weight averaged diameter of rain drops:
             dr = ( hyrho(k) * qr(k,j,i) / nr(k,j,i) * dpirho_l )**( 1.0_wp / 3.0_wp )
!
!--          Shape parameter of gamma distribution (Milbrandt and Yau, 2005;
!--          Stevens and Seifert, 2008):
             mu_r = 10.0_wp * ( 1.0_wp + TANH( 1.2E3_wp * ( dr - 1.4E-3_wp ) ) )
!
!--          Slope parameter of gamma distribution (Seifert, 2008):
             lambda_r = ( ( mu_r + 3.0_wp ) * ( mu_r + 2.0_wp ) *              &
                          ( mu_r + 1.0_wp ) )**( 1.0_wp / 3.0_wp ) / dr

             w_nr(k) = MAX( 0.1_wp, MIN( 20.0_wp,                              &
                                         a_term - b_term * ( 1.0_wp +          &
                                            c_term / lambda_r )**( -1.0_wp *   &
                                            ( mu_r + 1.0_wp ) )                &
                                        )                                      &
                          ) * flag
             w_qr(k) = MAX( 0.1_wp, MIN( 20.0_wp,                              &
                                         a_term - b_term * ( 1.0_wp +          &
                                            c_term / lambda_r )**( -1.0_wp *   &
                                            ( mu_r + 4.0_wp ) )                &
                                       )                                       &
                          ) * flag
          ELSE
             w_nr(k) = 0.0_wp
             w_qr(k) = 0.0_wp
          ENDIF
       ENDDO
!
!--    Adjust boundary values using surface data type.
!--    Upward facing non-natural
       surf_s = bc_h(0)%start_index(j,i)
       surf_e = bc_h(0)%end_index(j,i)
       DO  m = surf_s, surf_e
          k         = bc_h(0)%k(m)
          w_nr(k-1) = w_nr(k)
          w_qr(k-1) = w_qr(k)
       ENDDO
!
!--    Downward facing non-natural
       surf_s = bc_h(1)%start_index(j,i)
       surf_e = bc_h(1)%end_index(j,i)
       DO  m = surf_s, surf_e
          k         = bc_h(1)%k(m)
          w_nr(k+1) = w_nr(k)
          w_qr(k+1) = w_qr(k)
       ENDDO
!
!--    Neumann boundary condition at model top
       w_nr(nzt+1) = 0.0_wp
       w_qr(nzt+1) = 0.0_wp
!
!--    Compute Courant number
       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )

          c_nr(k) = 0.25_wp * ( w_nr(k-1) + 2.0_wp * w_nr(k) + w_nr(k+1) ) *   &
                    dt_micro * ddzu(k) * flag
          c_qr(k) = 0.25_wp * ( w_qr(k-1) + 2.0_wp * w_qr(k) + w_qr(k+1) ) *   &
                    dt_micro * ddzu(k) * flag
       ENDDO
!
!--    Limit slopes with monotonized centered (MC) limiter (van Leer, 1977):
       IF ( limiter_sedimentation )  THEN

          DO k = nzb+1, nzt
!
!--          Predetermine flag to mask topography
             flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )

             d_mean = 0.5_wp * ( qr(k+1,j,i) - qr(k-1,j,i) )
             d_min  = qr(k,j,i) - MIN( qr(k+1,j,i), qr(k,j,i), qr(k-1,j,i) )
             d_max  = MAX( qr(k+1,j,i), qr(k,j,i), qr(k-1,j,i) ) - qr(k,j,i)

             qr_slope(k) = SIGN(1.0_wp, d_mean) * MIN ( 2.0_wp * d_min,        &
                                                        2.0_wp * d_max,        &
                                                        ABS( d_mean ) ) * flag

             d_mean = 0.5_wp * ( nr(k+1,j,i) - nr(k-1,j,i) )
             d_min  = nr(k,j,i) - MIN( nr(k+1,j,i), nr(k,j,i), nr(k-1,j,i) )
             d_max  = MAX( nr(k+1,j,i), nr(k,j,i), nr(k-1,j,i) ) - nr(k,j,i)

             nr_slope(k) = SIGN(1.0_wp, d_mean) * MIN ( 2.0_wp * d_min,        &
                                                        2.0_wp * d_max,        &
                                                        ABS( d_mean ) ) * flag
          ENDDO

       ELSE

          nr_slope = 0.0_wp
          qr_slope = 0.0_wp

       ENDIF

       sed_nr(nzt+1) = 0.0_wp
       sed_qr(nzt+1) = 0.0_wp
!
!--    Compute sedimentation flux
       DO  k = nzt, nzb+1, -1
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )
!
!--       Sum up all rain drop number densities which contribute to the flux
!--       through k-1/2
          flux  = 0.0_wp
          z_run = 0.0_wp ! height above z(k)
          k_run = k
          c_run = MIN( 1.0_wp, c_nr(k) )
          DO WHILE ( c_run > 0.0_wp  .AND.  k_run <= nzt )
             flux  = flux + hyrho(k_run) *                                     &
                     ( nr(k_run,j,i) + nr_slope(k_run) * ( 1.0_wp - c_run ) *  &
                     0.5_wp ) * c_run * dzu(k_run) * flag
             z_run = z_run + dzu(k_run)            * flag
             k_run = k_run + 1                     * flag
             c_run = MIN( 1.0_wp, c_nr(k_run) - z_run * ddzu(k_run) ) * flag
          ENDDO
!
!--       It is not allowed to sediment more rain drop number density than
!--       available
          flux = MIN( flux,                                                    &
                      hyrho(k) * dzu(k+1) * nr(k,j,i) + sed_nr(k+1) * dt_micro )

          sed_nr(k) = flux / dt_micro * flag
          nr(k,j,i)  = nr(k,j,i) + ( sed_nr(k+1) - sed_nr(k) ) * ddzu(k+1) /     &
                                    hyrho(k) * dt_micro * flag
!
!--       Sum up all rain water content which contributes to the flux
!--       through k-1/2
          flux  = 0.0_wp
          z_run = 0.0_wp ! height above z(k)
          k_run = k
          c_run = MIN( 1.0_wp, c_qr(k) )

          DO WHILE ( c_run > 0.0_wp  .AND.  k_run <= nzt )

             flux  = flux + hyrho(k_run) *                                     &
                     ( qr(k_run,j,i) + qr_slope(k_run) * ( 1.0_wp - c_run ) *  &
                     0.5_wp ) * c_run * dzu(k_run) * flag
             z_run = z_run + dzu(k_run)            * flag
             k_run = k_run + 1                     * flag
             c_run = MIN( 1.0_wp, c_qr(k_run) - z_run * ddzu(k_run) ) * flag

          ENDDO
!
!--       It is not allowed to sediment more rain water content than available
          flux = MIN( flux,                                                    &
                      hyrho(k) * dzu(k) * qr(k,j,i) + sed_qr(k+1) * dt_micro )

          sed_qr(k) = flux / dt_micro * flag

          qr(k,j,i) = qr(k,j,i) + ( sed_qr(k+1) - sed_qr(k) ) * ddzu(k+1) /      &
                                hyrho(k) * dt_micro * flag
          q(k,j,i)  = q(k,j,i)  + ( sed_qr(k+1) - sed_qr(k) ) * ddzu(k+1) /      &
                                hyrho(k) * dt_micro * flag
          pt(k,j,i) = pt(k,j,i) - ( sed_qr(k+1) - sed_qr(k) ) * ddzu(k+1) /      &
                                hyrho(k) * lv_d_cp * d_exner(k) * dt_micro     &
                                                                * flag
!
!--       Compute the rain rate
          IF ( call_microphysics_at_all_substeps )  THEN
             prr(k,j,i) = prr(k,j,i) + sed_qr(k) / hyrho(k)                    &
                          * weight_substep(intermediate_timestep_count) * flag
          ELSE
             prr(k,j,i) = prr(k,j,i) + sed_qr(k) / hyrho(k) * flag
          ENDIF

       ENDDO

    END SUBROUTINE sedimentation_rain_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Computation of the precipitation amount due to gravitational settling of
!> rain and cloud (fog) droplets
!------------------------------------------------------------------------------!
    SUBROUTINE calc_precipitation_amount

       IMPLICIT NONE

       INTEGER(iwp) ::  i             !< running index x direction
       INTEGER(iwp) ::  j             !< running index y direction
       INTEGER(iwp) ::  k             !< running index y direction
       INTEGER(iwp) ::  m             !< running index surface elements

       IF ( ( dt_do2d_xy - time_do2d_xy ) < precipitation_amount_interval .AND.&
            ( .NOT. call_microphysics_at_all_substeps .OR.                     &
            intermediate_timestep_count == intermediate_timestep_count_max ) ) &
       THEN
!
!--       Run over all upward-facing surface elements, i.e. non-natural,
!--       natural and urban
          DO  m = 1, bc_h(0)%ns
             i = bc_h(0)%i(m)
             j = bc_h(0)%j(m)
             k = bc_h(0)%k(m)
             precipitation_amount(j,i) = precipitation_amount(j,i) +           &
                                               prr(k,j,i) * hyrho(k) * dt_3d
          ENDDO

       ENDIF

    END SUBROUTINE calc_precipitation_amount


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This subroutine computes the precipitation amount due to gravitational
!> settling of rain and cloud (fog) droplets
!------------------------------------------------------------------------------!
    SUBROUTINE calc_precipitation_amount_ij( i, j )

       IMPLICIT NONE

       INTEGER(iwp) ::  i             !< running index x direction
       INTEGER(iwp) ::  j             !< running index y direction
       INTEGER(iwp) ::  k             !< running index z direction
       INTEGER(iwp) ::  m             !< running index surface elements
       INTEGER(iwp) ::  surf_e        !< End index of surface elements at (j,i)-gridpoint
       INTEGER(iwp) ::  surf_s        !< Start index of surface elements at (j,i)-gridpoint

       IF ( ( dt_do2d_xy - time_do2d_xy ) < precipitation_amount_interval .AND.&
            ( .NOT. call_microphysics_at_all_substeps .OR.                     &
            intermediate_timestep_count == intermediate_timestep_count_max ) ) &
       THEN

          surf_s = bc_h(0)%start_index(j,i)
          surf_e = bc_h(0)%end_index(j,i)
          DO  m = surf_s, surf_e
             k                         = bc_h(0)%k(m)
             precipitation_amount(j,i) = precipitation_amount(j,i) +           &
                                               prr(k,j,i) * hyrho(k) * dt_3d
          ENDDO

       ENDIF

    END SUBROUTINE calc_precipitation_amount_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Computation of the diagnostic supersaturation sat, actual liquid water
!< temperature t_l and saturation water vapor mixing ratio q_s
!------------------------------------------------------------------------------!
    SUBROUTINE supersaturation ( i,j,k )

       IMPLICIT NONE

       INTEGER(iwp) ::  i   !< running index
       INTEGER(iwp) ::  j   !< running index
       INTEGER(iwp) ::  k   !< running index

       REAL(wp) ::  alpha   !< correction factor
!
!--    Actual liquid water temperature:
       t_l = exner(k) * pt(k,j,i)
!
!--    Calculate water vapor saturation pressure
       e_s = magnus( t_l )
!
!--    Computation of saturation mixing ratio:
       q_s   = rd_d_rv * e_s / ( hyp(k) - e_s )
!
!--    Correction factor
       alpha = rd_d_rv * lv_d_rd * lv_d_cp / ( t_l * t_l )
!
!--    Correction of the approximated value
!--    (see: Cuijpers + Duynkerke, 1993, JAS, 23)
       q_s   = q_s * ( 1.0_wp + alpha * q(k,j,i) ) / ( 1.0_wp + alpha * q_s )
!
!--    Supersaturation:
!--    Not in case of microphysics_kessler or microphysics_sat_adjust
!--    since qr is unallocated
       IF ( microphysics_seifert ) THEN
          sat = ( q(k,j,i) - qr(k,j,i) - qc(k,j,i) ) / q_s - 1.0_wp
       ELSEIF ( microphysics_morrison_no_rain ) THEN
          sat = ( q(k,j,i) - qc(k,j,i) ) / q_s - 1.0_wp
       ENDIF

    END SUBROUTINE supersaturation


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculation of the liquid water content (0%-or-100%-scheme). This scheme is
!> used by the one and the two moment cloud physics scheme. Using the two moment
!> scheme, this calculation results in the cloud water content.
!------------------------------------------------------------------------------!
    SUBROUTINE calc_liquid_water_content



       IMPLICIT NONE

       INTEGER(iwp) ::  i !<
       INTEGER(iwp) ::  j !<
       INTEGER(iwp) ::  k !<

       REAL(wp)     ::  flag !< flag to indicate first grid level above surface


       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )

!
!--             Call calculation of supersaturation located
                CALL supersaturation( i, j, k )

!
!--             Compute the liquid water content
                IF ( microphysics_seifert  .AND. .NOT. microphysics_morrison ) &
                THEN
                   IF ( ( q(k,j,i) - q_s - qr(k,j,i) ) > 0.0_wp )  THEN
                      qc(k,j,i) = ( q(k,j,i) - q_s - qr(k,j,i) ) * flag
                      ql(k,j,i) = ( qc(k,j,i) + qr(k,j,i) ) * flag
                   ELSE
                      IF ( q(k,j,i) < qr(k,j,i) )  q(k,j,i) = qr(k,j,i)
                      qc(k,j,i) = 0.0_wp
                      ql(k,j,i) = qr(k,j,i) * flag
                   ENDIF
                ELSEIF ( microphysics_morrison  .AND.  microphysics_seifert )  THEN
                   ql(k,j,i) = qc(k,j,i) + qr(k,j,i) * flag
                ELSEIF ( microphysics_morrison  .AND.  .NOT. microphysics_seifert )  THEN
                   ql(k,j,i) = qc(k,j,i)                    
                ELSE
                   IF ( ( q(k,j,i) - q_s ) > 0.0_wp )  THEN
                      qc(k,j,i) = ( q(k,j,i) - q_s ) * flag
                      ql(k,j,i) = qc(k,j,i) * flag
                   ELSE
                      qc(k,j,i) = 0.0_wp
                      ql(k,j,i) = 0.0_wp
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE calc_liquid_water_content

!------------------------------------------------------------------------------!
! Description:
! ------------
!> This function computes the gamma function (Press et al., 1992).
!> The gamma function is needed for the calculation of the evaporation
!> of rain drops.
!------------------------------------------------------------------------------!
    FUNCTION gamm( xx )

       IMPLICIT NONE

       INTEGER(iwp) ::  j            !<

       REAL(wp)     ::  gamm         !<
       REAL(wp)     ::  ser          !<
       REAL(wp)     ::  tmp          !<
       REAL(wp)     ::  x_gamm       !<
       REAL(wp)     ::  xx           !<
       REAL(wp)     ::  y_gamm       !<


       REAL(wp), PARAMETER  ::  stp = 2.5066282746310005_wp               !<
       REAL(wp), PARAMETER  ::  cof(6) = (/ 76.18009172947146_wp,      &
                                           -86.50532032941677_wp,      &
                                            24.01409824083091_wp,      &
                                            -1.231739572450155_wp,     &
                                             0.1208650973866179E-2_wp, &
                                            -0.5395239384953E-5_wp /)     !<

       x_gamm = xx
       y_gamm = x_gamm
       tmp = x_gamm + 5.5_wp
       tmp = ( x_gamm + 0.5_wp ) * LOG( tmp ) - tmp
       ser = 1.000000000190015_wp

       DO  j = 1, 6
          y_gamm = y_gamm + 1.0_wp
          ser    = ser + cof( j ) / y_gamm
       ENDDO

!
!--    Until this point the algorithm computes the logarithm of the gamma
!--    function. Hence, the exponential function is used.
!       gamm = EXP( tmp + LOG( stp * ser / x_gamm ) )
       gamm = EXP( tmp ) * stp * ser / x_gamm

       RETURN

    END FUNCTION gamm

 END MODULE bulk_cloud_model_mod
