!> @file urban_surface_mod.f90
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
! Copyright 2015-2019 Czech Technical University in Prague
! Copyright 2015-2020 Institute of Computer Science of the
!                     Czech Academy of Sciences, Prague
! Copyright 1997-2020 Leibniz Universitaet Hannover
!------------------------------------------------------------------------------!
!
! Current revisions:
! ------------------
! 
! 
! Former revisions:
! -----------------
! $Id: urban_surface_mod.f90 4457 2020-03-11 14:20:43Z raasch $
! use statement for exchange horiz added
! 
! 4442 2020-03-04 19:21:13Z suehring
! Change order of dimension in surface arrays %frac, %emissivity and %albedo 
! to allow for better vectorization in the radiation interactions.
! 
! 4441 2020-03-04 19:20:35Z suehring
! Removed wall_flags_static_0 from USE statements as it's not used within
! the module
! 
! 4329 2019-12-10 15:46:36Z motisi
! Renamed wall_flags_0 to wall_flags_static_0
! 
! 4309 2019-11-26 18:49:59Z suehring
! - Bugfix, include m_liq into restarts
! - Remove unused arrays for liquid water and saturation moisture at vertical
!   walls
!
! 4305 2019-11-25 11:15:40Z suehring
! Revision of some indoor-model parameters
! 
! 4259 2019-10-09 10:05:22Z suehring
! Instead of terminate the job in case the relative wall fractions do not 
! sum-up to one, give only an informative message and normalize the fractions.
! 
! 4258 2019-10-07 13:29:08Z suehring
! - Add checks to ensure that relative fractions of walls, windowns and green
!   surfaces sum-up to one.
! - Revise message calls dealing with local checks.
! 
! 4245 2019-09-30 08:40:37Z pavelkrc
! Initialize explicit per-surface parameters from building_surface_pars
! 
! 4238 2019-09-25 16:06:01Z suehring
! Indoor-model parameters for some building types adjusted in order to avoid
! unrealistically high indoor temperatures (S. Rissmann) 
! 
! 4230 2019-09-11 13:58:14Z suehring
! Bugfix, initialize canopy resistance. Even if no green fraction is set, 
! r_canopy must be initialized for output purposes. 
! 
! 4227 2019-09-10 18:04:34Z gronemeier
! implement new palm_date_time_mod
! 
! 4214 2019-09-02 15:57:02Z suehring
! Bugfix, missing initialization and clearing of soil-moisture tendency 
! (J.Resler)
! 
! 4182 2019-08-22 15:20:23Z scharf
! Corrected "Former revisions" section
! 
! 4168 2019-08-16 13:50:17Z suehring
! Replace function get_topography_top_index by topo_top_ind
! 
! 4148 2019-08-08 11:26:00Z suehring
! - Add anthropogenic heat output factors for heating and cooling to building
!   data base
! - Move definition of building_pars to usm_init_arrays since it is already
!   required in the indoor model
! 
! 4127 2019-07-30 14:47:10Z suehring
! Do not add anthopogenic energy during wall/soil spin-up 
! (merge from branch resler)
! 
! 4077 2019-07-09 13:27:11Z gronemeier
! Set roughness length z0 and z0h/q at ground-floor level to same value as
! those above ground-floor level
! 
! 4051 2019-06-24 13:58:30Z suehring
! Remove work-around for green surface fraction on buildings 
! (do not set it zero) 
! 
! 4050 2019-06-24 13:57:27Z suehring
! In order to avoid confusion with global control parameter, rename the 
! USM-internal flag spinup into during_spinup.
! 
! 3987 2019-05-22 09:52:13Z kanani
! Introduce alternative switch for debug output during timestepping
! 
! 3943 2019-05-02 09:50:41Z maronga
! Removed qsws_eb. Bugfix in calculation of qsws.
! 
! 3933 2019-04-25 12:33:20Z kanani
! Remove allocation of pt_2m, this is done in surface_mod now (surfaces%pt_2m)
! 
! 3921 2019-04-18 14:21:10Z suehring
! Undo accidentally commented initialization  
! 
! 3918 2019-04-18 13:33:11Z suehring
! Set green fraction to zero also at vertical surfaces
! 
! 3914 2019-04-17 16:02:02Z suehring
! In order to obtain correct surface temperature during spinup set window 
! fraction to zero (only during spinup) instead of just disabling
! time-integration of window-surface temperature. 
! 
! 3901 2019-04-16 16:17:02Z suehring
! Workaround - set green fraction to zero ( green-heat model crashes ).
! 
! 3896 2019-04-15 10:10:17Z suehring
! 
! 
! 3896 2019-04-15 10:10:17Z suehring
! Bugfix, wrong index used for accessing building_pars from PIDS 
! 
! 3885 2019-04-11 11:29:34Z kanani
! Changes related to global restructuring of location messages and introduction 
! of additional debug messages
! 
! 3882 2019-04-10 11:08:06Z suehring
! Avoid different type kinds
! Move definition of building-surface properties from declaration block
! to an extra routine
! 
! 3881 2019-04-10 09:31:22Z suehring
! Revise determination of local ground-floor level height. 
! Make level 3 initalization conform with Palm-input-data standard
! Move output of albedo and emissivity to radiation module
! 
! 3832 2019-03-28 13:16:58Z raasch
! instrumented with openmp directives
! 
! 3824 2019-03-27 15:56:16Z pavelkrc
! Remove unused imports
!
! 
! 3814 2019-03-26 08:40:31Z pavelkrc
! unused subroutine commented out
! 
! 3769 2019-02-28 10:16:49Z moh.hefny
! removed unused variables
!
! 3767 2019-02-27 08:18:02Z raasch
! unused variables removed from rrd-subroutines parameter list
! 
! 3748 2019-02-18 10:38:31Z suehring
! Revise conversion of waste-heat flux (do not divide by air density, will
! be done in diffusion_s)
! 
! 3745 2019-02-15 18:57:56Z suehring
! - Remove internal flag indoor_model (is a global control parameter)
! - add waste heat from buildings to the kinmatic heat flux
! - consider waste heat in restart data
! - remove unused USE statements
! 
! 3744 2019-02-15 18:38:58Z suehring
! fixed surface heat capacity in the building parameters
! convert the file back to unix format
!
! 3730 2019-02-11 11:26:47Z moh.hefny
! Formatting and clean-up (rvtils)
! 
! 3710 2019-01-30 18:11:19Z suehring
! Check if building type is set within a valid range.
! 
! 3705 2019-01-29 19:56:39Z suehring
! make nzb_wall public, required for virtual-measurements
! 
! 3704 2019-01-29 19:51:41Z suehring
! Some interface calls moved to module_interface + cleanup
! 
! 3655 2019-01-07 16:51:22Z knoop
! Implementation of the PALM module interface
! 
! 2007 2016-08-24 15:47:17Z kanani
! Initial revision
!
!
! Description:
! ------------
! 2016/6/9 - Initial version of the USM (Urban Surface Model)
!            authors: Jaroslav Resler, Pavel Krc
!                     (Czech Technical University in Prague and Institute of
!                      Computer Science of the Czech Academy of Sciences, Prague)
!            with contributions: Michal Belda, Nina Benesova, Ondrej Vlcek
!            partly inspired by PALM LSM (B. Maronga)
!            parameterizations of Ra checked with TUF3D (E. S. Krayenhoff)
!> Module for Urban Surface Model (USM)
!> The module includes:
!>    1. radiation model with direct/diffuse radiation, shading, reflections
!>       and integration with plant canopy
!>    2. wall and wall surface model
!>    3. surface layer energy balance
!>    4. anthropogenic heat (only from transportation so far)
!>    5. necessary auxiliary subroutines (reading inputs, writing outputs,
!>       restart simulations, ...)
!> It also make use of standard radiation and integrates it into
!> urban surface model.
!>
!> Further work:
!> -------------
!> 1. Remove global arrays surfouts, surfoutl and only keep track of radiosity
!>    from surfaces that are visible from local surfaces (i.e. there is a SVF
!>    where target is local). To do that, radiosity will be exchanged after each
!>    reflection step using MPI_Alltoall instead of current MPI_Allgather.
!>
!> 2. Temporarily large values of surface heat flux can be observed, up to
!>    1.2 Km/s, which seem to be not realistic.
!>
!> @todo Output of _av variables in case of restarts
!> @todo Revise flux conversion in energy-balance solver
!> @todo Check optimizations for RMA operations
!> @todo Alternatives for MPI_WIN_ALLOCATE? (causes problems with openmpi)
!> @todo Check for load imbalances in CPU measures, e.g. for exchange_horiz_prog
!>       factor 3 between min and max time
!> @todo Check divisions in wtend (etc.) calculations for possible division
!>       by zero, e.g. in case fraq(0,m) + fraq(1,m) = 0?!
!> @todo Use unit 90 for OPEN/CLOSE of input files (FK)
!> @todo Move plant canopy stuff into plant canopy code
!------------------------------------------------------------------------------!
 MODULE urban_surface_mod

    USE arrays_3d,                                                             &
        ONLY:  hyp, zu, pt, p, u, v, w, tend, exner, hyrho, prr, q, ql, vpt

    USE calc_mean_profile_mod,                                                 &
        ONLY:  calc_mean_profile

    USE basic_constants_and_equations_mod,                                     &
        ONLY:  c_p, g, kappa, pi, r_d, rho_l, l_v, sigma_sb

    USE control_parameters,                                                    &
        ONLY:  coupling_start_time, topography,                                &
               debug_output, debug_output_timestep, debug_string,              &
               dt_3d, humidity, indoor_model,                                  &
               intermediate_timestep_count, initializing_actions,              &
               intermediate_timestep_count_max, simulated_time, end_time,      &
               timestep_scheme, tsc, coupling_char, io_blocks, io_group,       &
               message_string, time_since_reference_point, surface_pressure,   &
               pt_surface, large_scale_forcing, lsf_surf,                      &
               spinup_pt_mean, spinup_time, time_do3d, dt_do3d,                &
               average_count_3d, varnamelength, urban_surface, dz

    USE bulk_cloud_model_mod,                                                  &
        ONLY: bulk_cloud_model, precipitation
               
    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point, log_point_s

    USE grid_variables,                                                        &
        ONLY:  dx, dy, ddx, ddy, ddx2, ddy2

    USE indices,                                                               &
        ONLY:  nx, ny, nnx, nny, nnz, nxl, nxlg, nxr, nxrg, nyn, nyng, nys,    &
               nysg, nzb, nzt, nbgp, topo_top_ind

    USE, INTRINSIC :: iso_c_binding 

    USE kinds
              
    USE palm_date_time_mod,                                                    &
        ONLY:  get_date_time, seconds_per_hour

    USE pegrid
        
    USE radiation_model_mod,                                                   &
        ONLY:  albedo_type, radiation_interaction,                             &
               radiation, rad_sw_in, rad_lw_in, rad_sw_out, rad_lw_out,        &
               force_radiation_call, iup_u, inorth_u, isouth_u, ieast_u,       &
               iwest_u, iup_l, inorth_l, isouth_l, ieast_l, iwest_l, id,       &
               nz_urban_b, nz_urban_t, unscheduled_radiation_calls

    USE statistics,                                                            &
        ONLY:  hom, statistic_regions

    USE surface_mod,                                                           &
        ONLY:  ind_pav_green, ind_veg_wall, ind_wat_win, surf_usm_h,           &
               surf_usm_v, surface_restore_elements


    IMPLICIT NONE

!
!-- USM model constants

    REAL(wp), PARAMETER ::                     &
              b_ch               = 6.04_wp,    &  !< Clapp & Hornberger exponent
              lambda_h_green_dry = 0.19_wp,    &  !< heat conductivity for dry soil   
              lambda_h_green_sm  = 3.44_wp,    &  !< heat conductivity of the soil matrix
              lambda_h_water     = 0.57_wp,    &  !< heat conductivity of water
              psi_sat            = -0.388_wp,  &  !< soil matrix potential at saturation
              rho_c_soil         = 2.19E6_wp,  &  !< volumetric heat capacity of soil
              rho_c_water        = 4.20E6_wp      !< volumetric heat capacity of water
!               m_max_depth        = 0.0002_wp     ! Maximum capacity of the water reservoir (m)

!
!-- Soil parameters I           alpha_vg,      l_vg_green,    n_vg, gamma_w_green_sat
    REAL(wp), DIMENSION(0:3,1:7), PARAMETER :: soil_pars = RESHAPE( (/     &
                                 3.83_wp,  1.250_wp, 1.38_wp,  6.94E-6_wp, &  !< soil 1
                                 3.14_wp, -2.342_wp, 1.28_wp,  1.16E-6_wp, &  !< soil 2
                                 0.83_wp, -0.588_wp, 1.25_wp,  0.26E-6_wp, &  !< soil 3
                                 3.67_wp, -1.977_wp, 1.10_wp,  2.87E-6_wp, &  !< soil 4
                                 2.65_wp,  2.500_wp, 1.10_wp,  1.74E-6_wp, &  !< soil 5
                                 1.30_wp,  0.400_wp, 1.20_wp,  0.93E-6_wp, &  !< soil 6
                                 0.00_wp,  0.00_wp,  0.00_wp,  0.57E-6_wp  &  !< soil 7
                                 /), (/ 4, 7 /) )

!
!-- Soil parameters II              swc_sat,     fc,   wilt,    swc_res  
    REAL(wp), DIMENSION(0:3,1:7), PARAMETER :: m_soil_pars = RESHAPE( (/ &
                                 0.403_wp, 0.244_wp, 0.059_wp, 0.025_wp, &  !< soil 1
                                 0.439_wp, 0.347_wp, 0.151_wp, 0.010_wp, &  !< soil 2
                                 0.430_wp, 0.383_wp, 0.133_wp, 0.010_wp, &  !< soil 3
                                 0.520_wp, 0.448_wp, 0.279_wp, 0.010_wp, &  !< soil 4
                                 0.614_wp, 0.541_wp, 0.335_wp, 0.010_wp, &  !< soil 5
                                 0.766_wp, 0.663_wp, 0.267_wp, 0.010_wp, &  !< soil 6
                                 0.472_wp, 0.323_wp, 0.171_wp, 0.000_wp  &  !< soil 7
                                 /), (/ 4, 7 /) )
!
!-- value 9999999.9_wp -> generic available or user-defined value must be set
!-- otherwise -> no generic variable and user setting is optional
    REAL(wp) :: alpha_vangenuchten = 9999999.9_wp,      &  !< NAMELIST alpha_vg
                field_capacity = 9999999.9_wp,          &  !< NAMELIST fc
                hydraulic_conductivity = 9999999.9_wp,  &  !< NAMELIST gamma_w_green_sat
                l_vangenuchten = 9999999.9_wp,          &  !< NAMELIST l_vg
                n_vangenuchten = 9999999.9_wp,          &  !< NAMELIST n_vg
                residual_moisture = 9999999.9_wp,       &  !< NAMELIST m_res
                saturation_moisture = 9999999.9_wp,     &  !< NAMELIST m_sat
                wilting_point = 9999999.9_wp               !< NAMELIST m_wilt
    
!
!-- configuration parameters (they can be setup in PALM config)
    LOGICAL ::  usm_material_model = .TRUE.        !< flag parameter indicating wheather the  model of heat in materials is used
    LOGICAL ::  usm_anthropogenic_heat = .FALSE.   !< flag parameter indicating wheather the anthropogenic heat sources
                                                   !< (e.g.transportation) are used
    LOGICAL ::  force_radiation_call_l = .FALSE.   !< flag parameter for unscheduled radiation model calls
    LOGICAL ::  read_wall_temp_3d = .FALSE.
    LOGICAL ::  usm_wall_mod = .FALSE.             !< reduces conductivity of the first 2 wall layers by factor 0.1


    INTEGER(iwp) ::  building_type = 1               !< default building type (preleminary setting)
    INTEGER(iwp) ::  land_category = 2               !< default category for land surface
    INTEGER(iwp) ::  wall_category = 2               !< default category for wall surface over pedestrian zone
    INTEGER(iwp) ::  pedestrian_category = 2         !< default category for wall surface in pedestrian zone
    INTEGER(iwp) ::  roof_category = 2               !< default category for root surface
    REAL(wp)     ::  roughness_concrete = 0.001_wp   !< roughness length of average concrete surface
!
!-- Indices of input attributes in building_pars for (above) ground floor level 
    INTEGER(iwp) ::  ind_alb_wall_agfl     = 38   !< index in input list for albedo_type of wall above ground floor level
    INTEGER(iwp) ::  ind_alb_wall_gfl      = 66   !< index in input list for albedo_type of wall ground floor level
    INTEGER(iwp) ::  ind_alb_wall_r        = 101  !< index in input list for albedo_type of wall roof
    INTEGER(iwp) ::  ind_alb_green_agfl    = 39   !< index in input list for albedo_type of green above ground floor level
    INTEGER(iwp) ::  ind_alb_green_gfl     = 78   !< index in input list for albedo_type of green ground floor level
    INTEGER(iwp) ::  ind_alb_green_r       = 117  !< index in input list for albedo_type of green roof
    INTEGER(iwp) ::  ind_alb_win_agfl      = 40   !< index in input list for albedo_type of window fraction above ground floor level
    INTEGER(iwp) ::  ind_alb_win_gfl       = 77   !< index in input list for albedo_type of window fraction ground floor level
    INTEGER(iwp) ::  ind_alb_win_r         = 115  !< index in input list for albedo_type of window fraction roof
    INTEGER(iwp) ::  ind_c_surface         = 45   !< index in input list for heat capacity wall surface
    INTEGER(iwp) ::  ind_c_surface_green   = 48   !< index in input list for heat capacity green surface
    INTEGER(iwp) ::  ind_c_surface_win     = 47   !< index in input list for heat capacity window surface
    INTEGER(iwp) ::  ind_emis_wall_agfl    = 14   !< index in input list for wall emissivity, above ground floor level
    INTEGER(iwp) ::  ind_emis_wall_gfl     = 32   !< index in input list for wall emissivity, ground floor level
    INTEGER(iwp) ::  ind_emis_wall_r       = 100  !< index in input list for wall emissivity, roof
    INTEGER(iwp) ::  ind_emis_green_agfl   = 15   !< index in input list for green emissivity, above ground floor level
    INTEGER(iwp) ::  ind_emis_green_gfl    = 34   !< index in input list for green emissivity, ground floor level
    INTEGER(iwp) ::  ind_emis_green_r      = 116  !< index in input list for green emissivity, roof
    INTEGER(iwp) ::  ind_emis_win_agfl     = 16   !< index in input list for window emissivity, above ground floor level
    INTEGER(iwp) ::  ind_emis_win_gfl      = 33   !< index in input list for window emissivity, ground floor level
    INTEGER(iwp) ::  ind_emis_win_r        = 113  !< index in input list for window emissivity, roof
    INTEGER(iwp) ::  ind_gflh              = 20   !< index in input list for ground floor level height
    INTEGER(iwp) ::  ind_green_frac_w_agfl = 2    !< index in input list for green fraction on wall, above ground floor level
    INTEGER(iwp) ::  ind_green_frac_w_gfl  = 23   !< index in input list for green fraction on wall, ground floor level
    INTEGER(iwp) ::  ind_green_frac_r_agfl = 3    !< index in input list for green fraction on roof, above ground floor level
    INTEGER(iwp) ::  ind_green_frac_r_gfl  = 24   !< index in input list for green fraction on roof, ground floor level
    INTEGER(iwp) ::  ind_hc1_agfl          = 6    !< index in input list for heat capacity at first wall layer,
                                                  !< above ground floor level
    INTEGER(iwp) ::  ind_hc1_gfl           = 26   !< index in input list for heat capacity at first wall layer, ground floor level
    INTEGER(iwp) ::  ind_hc1_wall_r        = 94   !< index in input list for heat capacity at first wall layer, roof
    INTEGER(iwp) ::  ind_hc1_win_agfl      = 83   !< index in input list for heat capacity at first window layer,
                                                  !< above ground floor level
    INTEGER(iwp) ::  ind_hc1_win_gfl       = 71   !< index in input list for heat capacity at first window layer,
                                                  !< ground floor level
    INTEGER(iwp) ::  ind_hc1_win_r         = 107  !< index in input list for heat capacity at first window layer, roof
    INTEGER(iwp) ::  ind_hc2_agfl          = 7    !< index in input list for heat capacity at second wall layer,
                                                  !< above ground floor level
    INTEGER(iwp) ::  ind_hc2_gfl           = 27   !< index in input list for heat capacity at second wall layer, ground floor level
    INTEGER(iwp) ::  ind_hc2_wall_r        = 95   !< index in input list for heat capacity at second wall layer, roof
    INTEGER(iwp) ::  ind_hc2_win_agfl      = 84   !< index in input list for heat capacity at second window layer,
                                                  !< above ground floor level
    INTEGER(iwp) ::  ind_hc2_win_gfl       = 72   !< index in input list for heat capacity at second window layer,
                                                  !< ground floor level
    INTEGER(iwp) ::  ind_hc2_win_r         = 108  !< index in input list for heat capacity at second window layer, roof
    INTEGER(iwp) ::  ind_hc3_agfl          = 8    !< index in input list for heat capacity at third wall layer,
                                                  !< above ground floor level
    INTEGER(iwp) ::  ind_hc3_gfl           = 28   !< index in input list for heat capacity at third wall layer, ground floor level
    INTEGER(iwp) ::  ind_hc3_wall_r        = 96   !< index in input list for heat capacity at third wall layer, roof
    INTEGER(iwp) ::  ind_hc3_win_agfl      = 85   !< index in input list for heat capacity at third window layer,
                                                  !< above ground floor level
    INTEGER(iwp) ::  ind_hc3_win_gfl       = 73   !< index in input list for heat capacity at third window layer,
                                                  !< ground floor level
    INTEGER(iwp) ::  ind_hc3_win_r         = 109  !< index in input list for heat capacity at third window layer, roof
    INTEGER(iwp) ::  ind_indoor_target_temp_summer = 12
    INTEGER(iwp) ::  ind_indoor_target_temp_winter = 13
    INTEGER(iwp) ::  ind_lai_r_agfl        = 4    !< index in input list for LAI on roof, above ground floor level
    INTEGER(iwp) ::  ind_lai_r_gfl         = 4  !< index in input list for LAI on roof, ground floor level
    INTEGER(iwp) ::  ind_lai_w_agfl        = 5    !< index in input list for LAI on wall, above ground floor level
    INTEGER(iwp) ::  ind_lai_w_gfl         = 25   !< index in input list for LAI on wall, ground floor level
    INTEGER(iwp) ::  ind_lambda_surf       = 46   !< index in input list for thermal conductivity of wall surface
    INTEGER(iwp) ::  ind_lambda_surf_green = 50   !< index in input list for thermal conductivity of green surface
    INTEGER(iwp) ::  ind_lambda_surf_win   = 49   !< index in input list for thermal conductivity of window surface
    INTEGER(iwp) ::  ind_tc1_agfl          = 9    !< index in input list for thermal conductivity at first wall layer,
                                                  !< above ground floor level
    INTEGER(iwp) ::  ind_tc1_gfl           = 29   !< index in input list for thermal conductivity at first wall layer,
                                                  !< ground floor level
    INTEGER(iwp) ::  ind_tc1_wall_r        = 97   !< index in input list for thermal conductivity at first wall layer, roof
    INTEGER(iwp) ::  ind_tc1_win_agfl      = 86   !< index in input list for thermal conductivity at first window layer,
                                                  !< above ground floor level
    INTEGER(iwp) ::  ind_tc1_win_gfl       = 74   !< index in input list for thermal conductivity at first window layer,
                                                  !< ground floor level
    INTEGER(iwp) ::  ind_tc1_win_r         = 110  !< index in input list for thermal conductivity at first window layer, roof
    INTEGER(iwp) ::  ind_tc2_agfl          = 10   !< index in input list for thermal conductivity at second wall layer,
                                                  !< above ground floor level
    INTEGER(iwp) ::  ind_tc2_gfl           = 30   !< index in input list for thermal conductivity at second wall layer,
                                                  !< ground floor level
    INTEGER(iwp) ::  ind_tc2_wall_r        = 98   !< index in input list for thermal conductivity at second wall layer, roof
    INTEGER(iwp) ::  ind_tc2_win_agfl      = 87   !< index in input list for thermal conductivity at second window layer,
                                                  !< above ground floor level
    INTEGER(iwp) ::  ind_tc2_win_gfl       = 75   !< index in input list for thermal conductivity at second window layer,
                                                  !< ground floor level
    INTEGER(iwp) ::  ind_tc2_win_r         = 111  !< index in input list for thermal conductivity at second window layer,
                                                  !< ground floor level
    INTEGER(iwp) ::  ind_tc3_agfl          = 11   !< index in input list for thermal conductivity at third wall layer,
                                                  !< above ground floor level
    INTEGER(iwp) ::  ind_tc3_gfl           = 31   !< index in input list for thermal conductivity at third wall layer,
                                                  !< ground floor level
    INTEGER(iwp) ::  ind_tc3_wall_r        = 99   !< index in input list for thermal conductivity at third wall layer, roof
    INTEGER(iwp) ::  ind_tc3_win_agfl      = 88   !< index in input list for thermal conductivity at third window layer,
                                                  !< above ground floor level
    INTEGER(iwp) ::  ind_tc3_win_gfl       = 76   !< index in input list for thermal conductivity at third window layer,
                                                  !< ground floor level
    INTEGER(iwp) ::  ind_tc3_win_r         = 112  !< index in input list for thermal conductivity at third window layer, roof
    INTEGER(iwp) ::  ind_thick_1_agfl      = 41   !< index for wall layer thickness - 1st layer above ground floor level
    INTEGER(iwp) ::  ind_thick_1_gfl       = 62   !< index for wall layer thickness - 1st layer ground floor level
    INTEGER(iwp) ::  ind_thick_1_wall_r    = 90   !< index for wall layer thickness - 1st layer roof
    INTEGER(iwp) ::  ind_thick_1_win_agfl  = 79   !< index for window layer thickness - 1st layer above ground floor level
    INTEGER(iwp) ::  ind_thick_1_win_gfl   = 67   !< index for window layer thickness - 1st layer ground floor level
    INTEGER(iwp) ::  ind_thick_1_win_r     = 103  !< index for window layer thickness - 1st layer roof
    INTEGER(iwp) ::  ind_thick_2_agfl      = 42   !< index for wall layer thickness - 2nd layer above ground floor level
    INTEGER(iwp) ::  ind_thick_2_gfl       = 63   !< index for wall layer thickness - 2nd layer ground floor level
    INTEGER(iwp) ::  ind_thick_2_wall_r    = 91   !< index for wall layer thickness - 2nd layer roof
    INTEGER(iwp) ::  ind_thick_2_win_agfl  = 80   !< index for window layer thickness - 2nd layer above ground floor level
    INTEGER(iwp) ::  ind_thick_2_win_gfl   = 68   !< index for window layer thickness - 2nd layer ground floor level
    INTEGER(iwp) ::  ind_thick_2_win_r     = 104  !< index for window layer thickness - 2nd layer roof
    INTEGER(iwp) ::  ind_thick_3_agfl      = 43   !< index for wall layer thickness - 3rd layer above ground floor level
    INTEGER(iwp) ::  ind_thick_3_gfl       = 64   !< index for wall layer thickness - 3rd layer ground floor level
    INTEGER(iwp) ::  ind_thick_3_wall_r    = 92   !< index for wall layer thickness - 3rd layer roof
    INTEGER(iwp) ::  ind_thick_3_win_agfl  = 81   !< index for window layer thickness - 3rd layer above ground floor level
    INTEGER(iwp) ::  ind_thick_3_win_gfl   = 69   !< index for window layer thickness - 3rd layer ground floor level  
    INTEGER(iwp) ::  ind_thick_3_win_r     = 105  !< index for window layer thickness - 3rd layer roof
    INTEGER(iwp) ::  ind_thick_4_agfl      = 44   !< index for wall layer thickness - 4th layer above ground floor level
    INTEGER(iwp) ::  ind_thick_4_gfl       = 65   !< index for wall layer thickness - 4th layer ground floor level
    INTEGER(iwp) ::  ind_thick_4_wall_r    = 93   !< index for wall layer thickness - 4st layer roof
    INTEGER(iwp) ::  ind_thick_4_win_agfl  = 82   !< index for window layer thickness - 4th layer above ground floor level
    INTEGER(iwp) ::  ind_thick_4_win_gfl   = 70   !< index for window layer thickness - 4th layer ground floor level
    INTEGER(iwp) ::  ind_thick_4_win_r     = 106  !< index for window layer thickness - 4th layer roof
    INTEGER(iwp) ::  ind_trans_agfl        = 17   !< index in input list for window transmissivity, above ground floor level
    INTEGER(iwp) ::  ind_trans_gfl         = 35   !< index in input list for window transmissivity, ground floor level
    INTEGER(iwp) ::  ind_trans_r           = 114  !< index in input list for window transmissivity, roof
    INTEGER(iwp) ::  ind_wall_frac_agfl    = 0    !< index in input list for wall fraction, above ground floor level
    INTEGER(iwp) ::  ind_wall_frac_gfl     = 21   !< index in input list for wall fraction, ground floor level
    INTEGER(iwp) ::  ind_wall_frac_r       = 89   !< index in input list for wall fraction, roof
    INTEGER(iwp) ::  ind_win_frac_agfl     = 1    !< index in input list for window fraction, above ground floor level
    INTEGER(iwp) ::  ind_win_frac_gfl      = 22   !< index in input list for window fraction, ground floor level
    INTEGER(iwp) ::  ind_win_frac_r        = 102  !< index in input list for window fraction, roof
    INTEGER(iwp) ::  ind_z0_agfl           = 18   !< index in input list for z0, above ground floor level
    INTEGER(iwp) ::  ind_z0_gfl            = 36   !< index in input list for z0, ground floor level
    INTEGER(iwp) ::  ind_z0qh_agfl         = 19   !< index in input list for z0h / z0q, above ground floor level
    INTEGER(iwp) ::  ind_z0qh_gfl          = 37   !< index in input list for z0h / z0q, ground floor level
    INTEGER(iwp) ::  ind_green_type_roof   = 118  !< index in input list for type of green roof
!
!-- Indices of input attributes in building_surface_pars (except for
!-- radiation-related, which are in radiation_model_mod)
    INTEGER(iwp) ::  ind_s_wall_frac                 = 0  !< index for wall fraction (0-1)
    INTEGER(iwp) ::  ind_s_win_frac                  = 1  !< index for window fraction (0-1)
    INTEGER(iwp) ::  ind_s_green_frac_w              = 2  !< index for green fraction on wall (0-1)
    INTEGER(iwp) ::  ind_s_green_frac_r              = 3  !< index for green fraction on roof (0-1)
    INTEGER(iwp) ::  ind_s_lai_r                     = 4  !< index for leaf area index of green fraction
    INTEGER(iwp) ::  ind_s_hc1                       = 5  !< index for heat capacity of wall layer 1
    INTEGER(iwp) ::  ind_s_hc2                       = 6  !< index for heat capacity of wall layer 2
    INTEGER(iwp) ::  ind_s_hc3                       = 7  !< index for heat capacity of wall layer 3
    INTEGER(iwp) ::  ind_s_tc1                       = 8  !< index for thermal conducivity of wall layer 1
    INTEGER(iwp) ::  ind_s_tc2                       = 9  !< index for thermal conducivity of wall layer 2
    INTEGER(iwp) ::  ind_s_tc3                       = 10 !< index for thermal conducivity of wall layer 3
    INTEGER(iwp) ::  ind_s_indoor_target_temp_summer = 11 !< index for indoor target summer temperature
    INTEGER(iwp) ::  ind_s_indoor_target_temp_winter = 12 !< index for indoor target winter temperature
    INTEGER(iwp) ::  ind_s_emis_wall                 = 13 !< index for emissivity of wall fraction (0-1)
    INTEGER(iwp) ::  ind_s_emis_green                = 14 !< index for emissivity of green fraction (0-1)
    INTEGER(iwp) ::  ind_s_emis_win                  = 15 !< index for emissivity o f window fraction (0-1)
    INTEGER(iwp) ::  ind_s_trans                     = 16 !< index for transmissivity of window fraction (0-1)
    INTEGER(iwp) ::  ind_s_z0                        = 17 !< index for roughness length for momentum (m)
    INTEGER(iwp) ::  ind_s_z0qh                      = 18 !< index for roughness length for heat (m)

    REAL(wp)  ::  roof_height_limit = 4.0_wp         !< height for distinguish between land surfaces and roofs
    REAL(wp)  ::  ground_floor_level = 4.0_wp        !< default ground floor level


    CHARACTER(37), DIMENSION(0:7), PARAMETER :: building_type_name = (/     &
                                   'user-defined                         ', &  !< type 0 
                                   'residential - 1950                   ', &  !< type  1 
                                   'residential 1951 - 2000              ', &  !< type  2
                                   'residential 2001 -                   ', &  !< type  3
                                   'office - 1950                        ', &  !< type  4 
                                   'office 1951 - 2000                   ', &  !< type  5
                                   'office 2001 -                        ', &  !< type  6
                                   'bridges                              '  &  !< type  7
                                                                     /)


!
!-- Building facade/wall/green/window properties (partly according to PIDS). 
!-- Initialization of building_pars is outsourced to usm_init_pars. This is
!-- needed because of the huge number of attributes given in building_pars 
!-- (>700), while intel and gfortran compiler have hard limit of continuation 
!-- lines of 511.
    REAL(wp), DIMENSION(0:135,1:7) ::  building_pars
!
!-- Type for surface temperatures at vertical walls. Is not necessary for horizontal walls. 
    TYPE t_surf_vertical
       REAL(wp), DIMENSION(:), ALLOCATABLE         :: t
    END TYPE t_surf_vertical
!
!-- Type for wall temperatures at vertical walls. Is not necessary for horizontal walls. 
    TYPE t_wall_vertical
       REAL(wp), DIMENSION(:,:), ALLOCATABLE       :: t
    END TYPE t_wall_vertical

    TYPE surf_type_usm
       REAL(wp), DIMENSION(:),   ALLOCATABLE ::  var_usm_1d  !< 1D prognostic variable
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  var_usm_2d  !< 2D prognostic variable
    END TYPE surf_type_usm
    
    TYPE(surf_type_usm), POINTER  ::  m_liq_usm_h,        &  !< liquid water reservoir (m), horizontal surface elements 
                                      m_liq_usm_h_p          !< progn. liquid water reservoir (m), horizontal surface elements 

    TYPE(surf_type_usm), TARGET   ::  m_liq_usm_h_1,      &  !<
                                      m_liq_usm_h_2          !<

    TYPE(surf_type_usm), TARGET ::  tm_liq_usm_h_m      !< liquid water reservoir tendency (m), horizontal surface elements 
!
!-- anthropogenic heat sources
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE        ::  aheat             !< daily average of anthropogenic heat (W/m2)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          ::  aheatprof         !< diurnal profiles of anthropogenic heat
                                                                         !< for particular layers
    INTEGER(iwp)                                   ::  naheatlayers = 1  !< number of layers of anthropogenic heat

!
!-- wall surface model
!-- wall surface model constants
    INTEGER(iwp), PARAMETER                        :: nzb_wall = 0       !< inner side of the wall model (to be switched)
    INTEGER(iwp), PARAMETER                        :: nzt_wall = 3       !< outer side of the wall model (to be switched)
    INTEGER(iwp), PARAMETER                        :: nzw = 4            !< number of wall layers (fixed for now)

    REAL(wp), DIMENSION(nzb_wall:nzt_wall)         :: zwn_default        = (/0.0242_wp, 0.0969_wp, 0.346_wp, 1.0_wp /)
    REAL(wp), DIMENSION(nzb_wall:nzt_wall)         :: zwn_default_window = (/0.25_wp,   0.5_wp,    0.75_wp,  1.0_wp /)
    REAL(wp), DIMENSION(nzb_wall:nzt_wall)         :: zwn_default_green  = (/0.25_wp,   0.5_wp,    0.75_wp,  1.0_wp /)
                                                                         !< normalized soil, wall and roof, window and
                                                                         !<green layer depths (m/m)

    REAL(wp)                                       :: wall_inner_temperature   = 295.0_wp    !< temperature of the inner wall
                                                                                             !< surface (~22 degrees C) (K)
    REAL(wp)                                       :: roof_inner_temperature   = 295.0_wp    !< temperature of the inner roof
                                                                                             !< surface (~22 degrees C) (K)
    REAL(wp)                                       :: soil_inner_temperature   = 288.0_wp    !< temperature of the deep soil
                                                                                             !< (~15 degrees C) (K)
    REAL(wp)                                       :: window_inner_temperature = 295.0_wp    !< temperature of the inner window
                                                                                             !< surface (~22 degrees C) (K)

    REAL(wp)                                       :: m_total = 0.0_wp  !< weighted total water content of the soil (m3/m3)
    INTEGER(iwp)                                   :: soil_type

!
!-- surface and material model variables for walls, ground, roofs
    REAL(wp), DIMENSION(:), ALLOCATABLE            :: zwn                !< normalized wall layer depths (m)
    REAL(wp), DIMENSION(:), ALLOCATABLE            :: zwn_window         !< normalized window layer depths (m)
    REAL(wp), DIMENSION(:), ALLOCATABLE            :: zwn_green          !< normalized green layer depths (m)

    REAL(wp), DIMENSION(:), POINTER                :: t_surf_wall_h
    REAL(wp), DIMENSION(:), POINTER                :: t_surf_wall_h_p 
    REAL(wp), DIMENSION(:), POINTER                :: t_surf_window_h
    REAL(wp), DIMENSION(:), POINTER                :: t_surf_window_h_p 
    REAL(wp), DIMENSION(:), POINTER                :: t_surf_green_h
    REAL(wp), DIMENSION(:), POINTER                :: t_surf_green_h_p 

    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET    :: t_surf_wall_h_1
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET    :: t_surf_wall_h_2
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET    :: t_surf_window_h_1
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET    :: t_surf_window_h_2
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET    :: t_surf_green_h_1
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET    :: t_surf_green_h_2

    TYPE(t_surf_vertical), DIMENSION(:), POINTER   ::  t_surf_wall_v
    TYPE(t_surf_vertical), DIMENSION(:), POINTER   ::  t_surf_wall_v_p
    TYPE(t_surf_vertical), DIMENSION(:), POINTER   ::  t_surf_window_v
    TYPE(t_surf_vertical), DIMENSION(:), POINTER   ::  t_surf_window_v_p
    TYPE(t_surf_vertical), DIMENSION(:), POINTER   ::  t_surf_green_v
    TYPE(t_surf_vertical), DIMENSION(:), POINTER   ::  t_surf_green_v_p

    TYPE(t_surf_vertical), DIMENSION(0:3), TARGET  :: t_surf_wall_v_1
    TYPE(t_surf_vertical), DIMENSION(0:3), TARGET  :: t_surf_wall_v_2
    TYPE(t_surf_vertical), DIMENSION(0:3), TARGET  :: t_surf_window_v_1
    TYPE(t_surf_vertical), DIMENSION(0:3), TARGET  :: t_surf_window_v_2
    TYPE(t_surf_vertical), DIMENSION(0:3), TARGET  :: t_surf_green_v_1
    TYPE(t_surf_vertical), DIMENSION(0:3), TARGET  :: t_surf_green_v_2

!
!-- Energy balance variables
!-- parameters of the land, roof and wall surfaces

    REAL(wp), DIMENSION(:,:), POINTER                :: t_wall_h, t_wall_h_p
    REAL(wp), DIMENSION(:,:), ALLOCATABLE, TARGET    :: t_wall_h_1, t_wall_h_2
    REAL(wp), DIMENSION(:,:), POINTER                :: t_window_h, t_window_h_p
    REAL(wp), DIMENSION(:,:), ALLOCATABLE, TARGET    :: t_window_h_1, t_window_h_2
    REAL(wp), DIMENSION(:,:), POINTER                :: t_green_h, t_green_h_p
    REAL(wp), DIMENSION(:,:), ALLOCATABLE, TARGET    :: t_green_h_1, t_green_h_2
    REAL(wp), DIMENSION(:,:), POINTER                :: swc_h, rootfr_h, wilt_h, fc_h, swc_sat_h, swc_h_p, swc_res_h
    REAL(wp), DIMENSION(:,:), ALLOCATABLE, TARGET    :: swc_h_1, rootfr_h_1, &
                                                        wilt_h_1, fc_h_1, swc_sat_h_1, swc_h_2, swc_res_h_1
    

    TYPE(t_wall_vertical), DIMENSION(:), POINTER   :: t_wall_v, t_wall_v_p
    TYPE(t_wall_vertical), DIMENSION(0:3), TARGET  :: t_wall_v_1, t_wall_v_2
    TYPE(t_wall_vertical), DIMENSION(:), POINTER   :: t_window_v, t_window_v_p
    TYPE(t_wall_vertical), DIMENSION(0:3), TARGET  :: t_window_v_1, t_window_v_2
    TYPE(t_wall_vertical), DIMENSION(:), POINTER   :: t_green_v, t_green_v_p
    TYPE(t_wall_vertical), DIMENSION(0:3), TARGET  :: t_green_v_1, t_green_v_2
!
!-- Surface and material parameters classes (surface_type)
!-- albedo, emissivity, lambda_surf, roughness, thickness, volumetric heat capacity, thermal conductivity
    INTEGER(iwp)                                   :: n_surface_types       !< number of the wall type categories
    INTEGER(iwp), PARAMETER                        :: n_surface_params = 9  !< number of parameters for each type of the wall
    INTEGER(iwp), PARAMETER                        :: ialbedo  = 1          !< albedo of the surface
    INTEGER(iwp), PARAMETER                        :: iemiss   = 2          !< emissivity of the surface
    INTEGER(iwp), PARAMETER                        :: ilambdas = 3          !< heat conductivity lambda S between surface
                                                                            !< and material ( W m-2 K-1 )
    INTEGER(iwp), PARAMETER                        :: irough   = 4          !< roughness length z0 for movements
    INTEGER(iwp), PARAMETER                        :: iroughh  = 5          !< roughness length z0h for scalars
                                                                            !< (heat, humidity,...)
    INTEGER(iwp), PARAMETER                        :: icsurf   = 6          !< Surface skin layer heat capacity (J m-2 K-1 )
    INTEGER(iwp), PARAMETER                        :: ithick   = 7          !< thickness of the surface (wall, roof, land)  ( m )
    INTEGER(iwp), PARAMETER                        :: irhoC    = 8          !< volumetric heat capacity rho*C of
                                                                            !< the material ( J m-3 K-1 )
    INTEGER(iwp), PARAMETER                        :: ilambdah = 9          !< thermal conductivity lambda H
                                                                            !< of the wall (W m-1 K-1 )
    CHARACTER(12), DIMENSION(:), ALLOCATABLE       :: surface_type_names    !< names of wall types (used only for reports)
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE        :: surface_type_codes    !< codes of wall types
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          :: surface_params        !< parameters of wall types

!
!-- interfaces of subroutines accessed from outside of this module
    INTERFACE usm_3d_data_averaging
       MODULE PROCEDURE usm_3d_data_averaging
    END INTERFACE usm_3d_data_averaging

    INTERFACE usm_boundary_condition
       MODULE PROCEDURE usm_boundary_condition
    END INTERFACE usm_boundary_condition

    INTERFACE usm_check_data_output
       MODULE PROCEDURE usm_check_data_output
    END INTERFACE usm_check_data_output
    
    INTERFACE usm_check_parameters
       MODULE PROCEDURE usm_check_parameters
    END INTERFACE usm_check_parameters
    
    INTERFACE usm_data_output_3d
       MODULE PROCEDURE usm_data_output_3d
    END INTERFACE usm_data_output_3d
    
    INTERFACE usm_define_netcdf_grid
       MODULE PROCEDURE usm_define_netcdf_grid
    END INTERFACE usm_define_netcdf_grid

    INTERFACE usm_init
       MODULE PROCEDURE usm_init
    END INTERFACE usm_init

    INTERFACE usm_init_arrays
       MODULE PROCEDURE usm_init_arrays
    END INTERFACE usm_init_arrays

    INTERFACE usm_material_heat_model
       MODULE PROCEDURE usm_material_heat_model
    END INTERFACE usm_material_heat_model
    
    INTERFACE usm_green_heat_model
       MODULE PROCEDURE usm_green_heat_model
    END INTERFACE usm_green_heat_model
    
    INTERFACE usm_parin
       MODULE PROCEDURE usm_parin
    END INTERFACE usm_parin

    INTERFACE usm_rrd_local 
       MODULE PROCEDURE usm_rrd_local
    END INTERFACE usm_rrd_local

    INTERFACE usm_surface_energy_balance
       MODULE PROCEDURE usm_surface_energy_balance
    END INTERFACE usm_surface_energy_balance
    
    INTERFACE usm_swap_timelevel
       MODULE PROCEDURE usm_swap_timelevel
    END INTERFACE usm_swap_timelevel
        
    INTERFACE usm_wrd_local
       MODULE PROCEDURE usm_wrd_local
    END INTERFACE usm_wrd_local

    
    SAVE

    PRIVATE 

!
!-- Public functions
    PUBLIC usm_boundary_condition, usm_check_parameters, usm_init,               &
           usm_rrd_local,                                                        & 
           usm_surface_energy_balance, usm_material_heat_model,                  &
           usm_swap_timelevel, usm_check_data_output, usm_3d_data_averaging,     &
           usm_data_output_3d, usm_define_netcdf_grid, usm_parin,                &
           usm_wrd_local, usm_init_arrays

!
!-- Public parameters, constants and initial values
    PUBLIC usm_anthropogenic_heat, usm_material_model, usm_wall_mod, &
           usm_green_heat_model, building_pars,                      &
           nzb_wall, nzt_wall, t_wall_h, t_wall_v,                   &
           t_window_h, t_window_v, building_type



 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> This subroutine creates the necessary indices of the urban surfaces
!> and plant canopy and it allocates the needed arrays for USM
!------------------------------------------------------------------------------!
    SUBROUTINE usm_init_arrays
    
        IMPLICIT NONE
       
        INTEGER(iwp) ::  l

        IF ( debug_output )  CALL debug_message( 'usm_init_arrays', 'start' )

!
!--     Allocate radiation arrays which are part of the new data type. 
!--     For horizontal surfaces.
        ALLOCATE ( surf_usm_h%surfhf(1:surf_usm_h%ns)    )
        ALLOCATE ( surf_usm_h%rad_net_l(1:surf_usm_h%ns) )
!
!--     For vertical surfaces
        DO  l = 0, 3
           ALLOCATE ( surf_usm_v(l)%surfhf(1:surf_usm_v(l)%ns)    )
           ALLOCATE ( surf_usm_v(l)%rad_net_l(1:surf_usm_v(l)%ns) )
        ENDDO

!
!--     Wall surface model
!--     allocate arrays for wall surface model and define pointers
!--     allocate array of wall types and wall parameters
        ALLOCATE ( surf_usm_h%surface_types(1:surf_usm_h%ns)      )
        ALLOCATE ( surf_usm_h%building_type(1:surf_usm_h%ns)      )
        ALLOCATE ( surf_usm_h%building_type_name(1:surf_usm_h%ns) )
        surf_usm_h%building_type      = 0
        surf_usm_h%building_type_name = 'none'
        DO  l = 0, 3
           ALLOCATE ( surf_usm_v(l)%surface_types(1:surf_usm_v(l)%ns)      )
           ALLOCATE ( surf_usm_v(l)%building_type(1:surf_usm_v(l)%ns)      )
           ALLOCATE ( surf_usm_v(l)%building_type_name(1:surf_usm_v(l)%ns) )
           surf_usm_v(l)%building_type      = 0
           surf_usm_v(l)%building_type_name = 'none'
        ENDDO
!
!--     Allocate albedo_type and albedo. Each surface element
!--     has 3 values, 0: wall fraction, 1: green fraction, 2: window fraction.
        ALLOCATE ( surf_usm_h%albedo_type(1:surf_usm_h%ns,0:2) )
        ALLOCATE ( surf_usm_h%albedo(1:surf_usm_h%ns,0:2)      )
        surf_usm_h%albedo_type = albedo_type
        DO  l = 0, 3
           ALLOCATE ( surf_usm_v(l)%albedo_type(1:surf_usm_v(l)%ns,0:2) )
           ALLOCATE ( surf_usm_v(l)%albedo(1:surf_usm_v(l)%ns,0:2)      )
           surf_usm_v(l)%albedo_type = albedo_type
        ENDDO       

!
!--     Allocate indoor target temperature for summer and winter
        ALLOCATE ( surf_usm_h%target_temp_summer(1:surf_usm_h%ns) )
        ALLOCATE ( surf_usm_h%target_temp_winter(1:surf_usm_h%ns) )
        DO  l = 0, 3
           ALLOCATE ( surf_usm_v(l)%target_temp_summer(1:surf_usm_v(l)%ns) )
           ALLOCATE ( surf_usm_v(l)%target_temp_winter(1:surf_usm_v(l)%ns) )
        ENDDO
!
!--     In case the indoor model is applied, allocate memory for waste heat 
!--     and indoor temperature.
        IF ( indoor_model )  THEN
           ALLOCATE ( surf_usm_h%waste_heat(1:surf_usm_h%ns) )
           surf_usm_h%waste_heat = 0.0_wp
           DO  l = 0, 3
              ALLOCATE ( surf_usm_v(l)%waste_heat(1:surf_usm_v(l)%ns) )
              surf_usm_v(l)%waste_heat = 0.0_wp
           ENDDO
        ENDIF
!
!--     Allocate flag indicating ground floor level surface elements
        ALLOCATE ( surf_usm_h%ground_level(1:surf_usm_h%ns) ) 
        DO  l = 0, 3
           ALLOCATE ( surf_usm_v(l)%ground_level(1:surf_usm_v(l)%ns) )
        ENDDO   
!
!--      Allocate arrays for relative surface fraction. 
!--      0 - wall fraction, 1 - green fraction, 2 - window fraction
         ALLOCATE ( surf_usm_h%frac(1:surf_usm_h%ns,0:2) )
         surf_usm_h%frac = 0.0_wp
         DO  l = 0, 3
            ALLOCATE ( surf_usm_v(l)%frac(1:surf_usm_v(l)%ns,0:2) )
            surf_usm_v(l)%frac = 0.0_wp
         ENDDO

!
!--     wall and roof surface parameters. First for horizontal surfaces
        ALLOCATE ( surf_usm_h%isroof_surf(1:surf_usm_h%ns)        )
        ALLOCATE ( surf_usm_h%lambda_surf(1:surf_usm_h%ns)        )
        ALLOCATE ( surf_usm_h%lambda_surf_window(1:surf_usm_h%ns) )
        ALLOCATE ( surf_usm_h%lambda_surf_green(1:surf_usm_h%ns)  )
        ALLOCATE ( surf_usm_h%c_surface(1:surf_usm_h%ns)          )
        ALLOCATE ( surf_usm_h%c_surface_window(1:surf_usm_h%ns)   )
        ALLOCATE ( surf_usm_h%c_surface_green(1:surf_usm_h%ns)    )
        ALLOCATE ( surf_usm_h%transmissivity(1:surf_usm_h%ns)     )
        ALLOCATE ( surf_usm_h%lai(1:surf_usm_h%ns)                )
        ALLOCATE ( surf_usm_h%emissivity(1:surf_usm_h%ns,0:2)     )
        ALLOCATE ( surf_usm_h%r_a(1:surf_usm_h%ns)                )
        ALLOCATE ( surf_usm_h%r_a_green(1:surf_usm_h%ns)          )
        ALLOCATE ( surf_usm_h%r_a_window(1:surf_usm_h%ns)         )
        ALLOCATE ( surf_usm_h%green_type_roof(1:surf_usm_h%ns)    )
        ALLOCATE ( surf_usm_h%r_s(1:surf_usm_h%ns)                )
        
!
!--     For vertical surfaces.
        DO  l = 0, 3
           ALLOCATE ( surf_usm_v(l)%lambda_surf(1:surf_usm_v(l)%ns)        )
           ALLOCATE ( surf_usm_v(l)%c_surface(1:surf_usm_v(l)%ns)          )
           ALLOCATE ( surf_usm_v(l)%lambda_surf_window(1:surf_usm_v(l)%ns) )
           ALLOCATE ( surf_usm_v(l)%c_surface_window(1:surf_usm_v(l)%ns)   )
           ALLOCATE ( surf_usm_v(l)%lambda_surf_green(1:surf_usm_v(l)%ns)  )
           ALLOCATE ( surf_usm_v(l)%c_surface_green(1:surf_usm_v(l)%ns)    )
           ALLOCATE ( surf_usm_v(l)%transmissivity(1:surf_usm_v(l)%ns)     )
           ALLOCATE ( surf_usm_v(l)%lai(1:surf_usm_v(l)%ns)                )
           ALLOCATE ( surf_usm_v(l)%emissivity(1:surf_usm_v(l)%ns,0:2)     )
           ALLOCATE ( surf_usm_v(l)%r_a(1:surf_usm_v(l)%ns)                )
           ALLOCATE ( surf_usm_v(l)%r_a_green(1:surf_usm_v(l)%ns)          )
           ALLOCATE ( surf_usm_v(l)%r_a_window(1:surf_usm_v(l)%ns)         )           
           ALLOCATE ( surf_usm_v(l)%r_s(1:surf_usm_v(l)%ns)                )
        ENDDO

!       
!--     allocate wall and roof material parameters. First for horizontal surfaces
        ALLOCATE ( surf_usm_h%thickness_wall(1:surf_usm_h%ns)                    )
        ALLOCATE ( surf_usm_h%thickness_window(1:surf_usm_h%ns)                  )
        ALLOCATE ( surf_usm_h%thickness_green(1:surf_usm_h%ns)                   )
        ALLOCATE ( surf_usm_h%lambda_h(nzb_wall:nzt_wall,1:surf_usm_h%ns)        )
        ALLOCATE ( surf_usm_h%rho_c_wall(nzb_wall:nzt_wall,1:surf_usm_h%ns)      )
        ALLOCATE ( surf_usm_h%lambda_h_window(nzb_wall:nzt_wall,1:surf_usm_h%ns) )
        ALLOCATE ( surf_usm_h%rho_c_window(nzb_wall:nzt_wall,1:surf_usm_h%ns)    )
        ALLOCATE ( surf_usm_h%lambda_h_green(nzb_wall:nzt_wall,1:surf_usm_h%ns)  )
        ALLOCATE ( surf_usm_h%rho_c_green(nzb_wall:nzt_wall,1:surf_usm_h%ns)     )

        ALLOCATE ( surf_usm_h%rho_c_total_green(nzb_wall:nzt_wall,1:surf_usm_h%ns)    )
        ALLOCATE ( surf_usm_h%n_vg_green(1:surf_usm_h%ns)                             )
        ALLOCATE ( surf_usm_h%alpha_vg_green(1:surf_usm_h%ns)                         )
        ALLOCATE ( surf_usm_h%l_vg_green(1:surf_usm_h%ns)                             )
        ALLOCATE ( surf_usm_h%gamma_w_green_sat(nzb_wall:nzt_wall+1,1:surf_usm_h%ns)  )
        ALLOCATE ( surf_usm_h%lambda_w_green(nzb_wall:nzt_wall,1:surf_usm_h%ns)       )
        ALLOCATE ( surf_usm_h%gamma_w_green(nzb_wall:nzt_wall,1:surf_usm_h%ns)        )
        ALLOCATE ( surf_usm_h%tswc_h_m(nzb_wall:nzt_wall,1:surf_usm_h%ns)             )

!
!--     For vertical surfaces.
        DO  l = 0, 3
           ALLOCATE ( surf_usm_v(l)%thickness_wall(1:surf_usm_v(l)%ns)                    )
           ALLOCATE ( surf_usm_v(l)%thickness_window(1:surf_usm_v(l)%ns)                  )
           ALLOCATE ( surf_usm_v(l)%thickness_green(1:surf_usm_v(l)%ns)                   )
           ALLOCATE ( surf_usm_v(l)%lambda_h(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns)        )
           ALLOCATE ( surf_usm_v(l)%rho_c_wall(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns)      )
           ALLOCATE ( surf_usm_v(l)%lambda_h_window(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns) )
           ALLOCATE ( surf_usm_v(l)%rho_c_window(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns)    )
           ALLOCATE ( surf_usm_v(l)%lambda_h_green(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns)  )
           ALLOCATE ( surf_usm_v(l)%rho_c_green(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns)     )
        ENDDO

!
!--     allocate green wall and roof vegetation and soil parameters. First horizontal surfaces
        ALLOCATE ( surf_usm_h%g_d(1:surf_usm_h%ns)              )
        ALLOCATE ( surf_usm_h%c_liq(1:surf_usm_h%ns)            )
        ALLOCATE ( surf_usm_h%qsws_liq(1:surf_usm_h%ns)         )
        ALLOCATE ( surf_usm_h%qsws_veg(1:surf_usm_h%ns)         )
        ALLOCATE ( surf_usm_h%r_canopy(1:surf_usm_h%ns)         )
        ALLOCATE ( surf_usm_h%r_canopy_min(1:surf_usm_h%ns)     )
        ALLOCATE ( surf_usm_h%pt_10cm(1:surf_usm_h%ns)          ) 

!
!--     For vertical surfaces.
        DO  l = 0, 3
          ALLOCATE ( surf_usm_v(l)%g_d(1:surf_usm_v(l)%ns)              )
          ALLOCATE ( surf_usm_v(l)%c_liq(1:surf_usm_v(l)%ns)            )
          ALLOCATE ( surf_usm_v(l)%qsws_liq(1:surf_usm_v(l)%ns)         )
          ALLOCATE ( surf_usm_v(l)%qsws_veg(1:surf_usm_v(l)%ns)         )
          ALLOCATE ( surf_usm_v(l)%r_canopy(1:surf_usm_v(l)%ns)         )
          ALLOCATE ( surf_usm_v(l)%r_canopy_min(1:surf_usm_v(l)%ns)     )
          ALLOCATE ( surf_usm_v(l)%pt_10cm(1:surf_usm_v(l)%ns)          )
        ENDDO

!
!--     allocate wall and roof layers sizes. For horizontal surfaces.
        ALLOCATE ( zwn(nzb_wall:nzt_wall)                                        )
        ALLOCATE ( surf_usm_h%dz_wall(nzb_wall:nzt_wall+1,1:surf_usm_h%ns)       )
        ALLOCATE ( zwn_window(nzb_wall:nzt_wall)                                 )
        ALLOCATE ( surf_usm_h%dz_window(nzb_wall:nzt_wall+1,1:surf_usm_h%ns)     )
        ALLOCATE ( zwn_green(nzb_wall:nzt_wall)                                  )
        ALLOCATE ( surf_usm_h%dz_green(nzb_wall:nzt_wall+1,1:surf_usm_h%ns)      )
        ALLOCATE ( surf_usm_h%ddz_wall(nzb_wall:nzt_wall+1,1:surf_usm_h%ns)      )
        ALLOCATE ( surf_usm_h%dz_wall_stag(nzb_wall:nzt_wall,1:surf_usm_h%ns)    )
        ALLOCATE ( surf_usm_h%ddz_wall_stag(nzb_wall:nzt_wall,1:surf_usm_h%ns)   )
        ALLOCATE ( surf_usm_h%zw(nzb_wall:nzt_wall,1:surf_usm_h%ns)              )
        ALLOCATE ( surf_usm_h%ddz_window(nzb_wall:nzt_wall+1,1:surf_usm_h%ns)    )
        ALLOCATE ( surf_usm_h%dz_window_stag(nzb_wall:nzt_wall,1:surf_usm_h%ns)  )
        ALLOCATE ( surf_usm_h%ddz_window_stag(nzb_wall:nzt_wall,1:surf_usm_h%ns) )
        ALLOCATE ( surf_usm_h%zw_window(nzb_wall:nzt_wall,1:surf_usm_h%ns)       )
        ALLOCATE ( surf_usm_h%ddz_green(nzb_wall:nzt_wall+1,1:surf_usm_h%ns)     )
        ALLOCATE ( surf_usm_h%dz_green_stag(nzb_wall:nzt_wall,1:surf_usm_h%ns)   )
        ALLOCATE ( surf_usm_h%ddz_green_stag(nzb_wall:nzt_wall,1:surf_usm_h%ns)  )
        ALLOCATE ( surf_usm_h%zw_green(nzb_wall:nzt_wall,1:surf_usm_h%ns)        )

!
!--     For vertical surfaces.
        DO  l = 0, 3
           ALLOCATE ( surf_usm_v(l)%dz_wall(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns)       )
           ALLOCATE ( surf_usm_v(l)%dz_window(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns)     )
           ALLOCATE ( surf_usm_v(l)%dz_green(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns)      )
           ALLOCATE ( surf_usm_v(l)%ddz_wall(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns)      )
           ALLOCATE ( surf_usm_v(l)%dz_wall_stag(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns)    )
           ALLOCATE ( surf_usm_v(l)%ddz_wall_stag(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns)   )
           ALLOCATE ( surf_usm_v(l)%zw(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns)              )
           ALLOCATE ( surf_usm_v(l)%ddz_window(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns)    )
           ALLOCATE ( surf_usm_v(l)%dz_window_stag(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns)  )
           ALLOCATE ( surf_usm_v(l)%ddz_window_stag(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns) )
           ALLOCATE ( surf_usm_v(l)%zw_window(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns)       )
           ALLOCATE ( surf_usm_v(l)%ddz_green(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns)     )
           ALLOCATE ( surf_usm_v(l)%dz_green_stag(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns)   )
           ALLOCATE ( surf_usm_v(l)%ddz_green_stag(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns)  )
           ALLOCATE ( surf_usm_v(l)%zw_green(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns)        )
        ENDDO

!
!--     allocate wall and roof temperature arrays, for horizontal walls
!
!--     Allocate if required. Note, in case of restarts, some of these arrays 
!--     might be already allocated.
        IF ( .NOT. ALLOCATED( t_surf_wall_h_1 ) )                              &
           ALLOCATE ( t_surf_wall_h_1(1:surf_usm_h%ns) )
        IF ( .NOT. ALLOCATED( t_surf_wall_h_2 ) )                              &
           ALLOCATE ( t_surf_wall_h_2(1:surf_usm_h%ns) )
        IF ( .NOT. ALLOCATED( t_wall_h_1 ) )                                   &           
           ALLOCATE ( t_wall_h_1(nzb_wall:nzt_wall+1,1:surf_usm_h%ns) ) 
        IF ( .NOT. ALLOCATED( t_wall_h_2 ) )                                   &           
           ALLOCATE ( t_wall_h_2(nzb_wall:nzt_wall+1,1:surf_usm_h%ns) )          
        IF ( .NOT. ALLOCATED( t_surf_window_h_1 ) )                            &
           ALLOCATE ( t_surf_window_h_1(1:surf_usm_h%ns) )
        IF ( .NOT. ALLOCATED( t_surf_window_h_2 ) )                            &
           ALLOCATE ( t_surf_window_h_2(1:surf_usm_h%ns) )
        IF ( .NOT. ALLOCATED( t_window_h_1 ) )                                 &           
           ALLOCATE ( t_window_h_1(nzb_wall:nzt_wall+1,1:surf_usm_h%ns) ) 
        IF ( .NOT. ALLOCATED( t_window_h_2 ) )                                 &           
           ALLOCATE ( t_window_h_2(nzb_wall:nzt_wall+1,1:surf_usm_h%ns) )          
        IF ( .NOT. ALLOCATED( t_surf_green_h_1 ) )                             &
           ALLOCATE ( t_surf_green_h_1(1:surf_usm_h%ns) )
        IF ( .NOT. ALLOCATED( t_surf_green_h_2 ) )                             &
           ALLOCATE ( t_surf_green_h_2(1:surf_usm_h%ns) )
        IF ( .NOT. ALLOCATED( t_green_h_1 ) )                                  &           
           ALLOCATE ( t_green_h_1(nzb_wall:nzt_wall+1,1:surf_usm_h%ns) ) 
        IF ( .NOT. ALLOCATED( t_green_h_2 ) )                                  &           
           ALLOCATE ( t_green_h_2(nzb_wall:nzt_wall+1,1:surf_usm_h%ns) )          
        IF ( .NOT. ALLOCATED( swc_h_1 ) )                                      &           
           ALLOCATE ( swc_h_1(nzb_wall:nzt_wall+1,1:surf_usm_h%ns) ) 
        IF ( .NOT. ALLOCATED( swc_sat_h_1 ) )                                  &           
           ALLOCATE ( swc_sat_h_1(nzb_wall:nzt_wall+1,1:surf_usm_h%ns) ) 
        IF ( .NOT. ALLOCATED( swc_res_h_1 ) )                                  &           
           ALLOCATE ( swc_res_h_1(nzb_wall:nzt_wall+1,1:surf_usm_h%ns) ) 
        IF ( .NOT. ALLOCATED( swc_h_2 ) )                                      &           
           ALLOCATE ( swc_h_2(nzb_wall:nzt_wall+1,1:surf_usm_h%ns) )
        IF ( .NOT. ALLOCATED( rootfr_h_1 ) )                                   &           
           ALLOCATE ( rootfr_h_1(nzb_wall:nzt_wall+1,1:surf_usm_h%ns) ) 
        IF ( .NOT. ALLOCATED( wilt_h_1 ) )                                     &           
           ALLOCATE ( wilt_h_1(nzb_wall:nzt_wall+1,1:surf_usm_h%ns) ) 
        IF ( .NOT. ALLOCATED( fc_h_1 ) )                                       &           
           ALLOCATE ( fc_h_1(nzb_wall:nzt_wall+1,1:surf_usm_h%ns) ) 

        IF ( .NOT. ALLOCATED( m_liq_usm_h_1%var_usm_1d ) )                     &
           ALLOCATE ( m_liq_usm_h_1%var_usm_1d(1:surf_usm_h%ns) )
        IF ( .NOT. ALLOCATED( m_liq_usm_h_2%var_usm_1d ) )                     &
           ALLOCATE ( m_liq_usm_h_2%var_usm_1d(1:surf_usm_h%ns) )
           
!           
!--     initial assignment of the pointers
        t_wall_h    => t_wall_h_1;   t_wall_h_p   => t_wall_h_2
        t_window_h  => t_window_h_1; t_window_h_p => t_window_h_2
        t_green_h   => t_green_h_1;  t_green_h_p  => t_green_h_2
        t_surf_wall_h   => t_surf_wall_h_1;   t_surf_wall_h_p   => t_surf_wall_h_2           
        t_surf_window_h => t_surf_window_h_1; t_surf_window_h_p => t_surf_window_h_2  
        t_surf_green_h  => t_surf_green_h_1;  t_surf_green_h_p  => t_surf_green_h_2           
        m_liq_usm_h     => m_liq_usm_h_1;     m_liq_usm_h_p     => m_liq_usm_h_2
        swc_h     => swc_h_1; swc_h_p => swc_h_2
        swc_sat_h => swc_sat_h_1
        swc_res_h => swc_res_h_1
        rootfr_h  => rootfr_h_1
        wilt_h    => wilt_h_1
        fc_h      => fc_h_1

!
!--     allocate wall and roof temperature arrays, for vertical walls if required
!
!--     Allocate if required. Note, in case of restarts, some of these arrays 
!--     might be already allocated.
        DO  l = 0, 3
           IF ( .NOT. ALLOCATED( t_surf_wall_v_1(l)%t ) )                      &
              ALLOCATE ( t_surf_wall_v_1(l)%t(1:surf_usm_v(l)%ns) )
           IF ( .NOT. ALLOCATED( t_surf_wall_v_2(l)%t ) )                      &
              ALLOCATE ( t_surf_wall_v_2(l)%t(1:surf_usm_v(l)%ns) )
           IF ( .NOT. ALLOCATED( t_wall_v_1(l)%t ) )                           &           
              ALLOCATE ( t_wall_v_1(l)%t(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns) ) 
           IF ( .NOT. ALLOCATED( t_wall_v_2(l)%t ) )                           &           
              ALLOCATE ( t_wall_v_2(l)%t(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns) )  
           IF ( .NOT. ALLOCATED( t_surf_window_v_1(l)%t ) )                    &
              ALLOCATE ( t_surf_window_v_1(l)%t(1:surf_usm_v(l)%ns) )
           IF ( .NOT. ALLOCATED( t_surf_window_v_2(l)%t ) )                    &
              ALLOCATE ( t_surf_window_v_2(l)%t(1:surf_usm_v(l)%ns) )
           IF ( .NOT. ALLOCATED( t_window_v_1(l)%t ) )                         &           
              ALLOCATE ( t_window_v_1(l)%t(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns) ) 
           IF ( .NOT. ALLOCATED( t_window_v_2(l)%t ) )                         &           
              ALLOCATE ( t_window_v_2(l)%t(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns) )  
           IF ( .NOT. ALLOCATED( t_surf_green_v_1(l)%t ) )                     &
              ALLOCATE ( t_surf_green_v_1(l)%t(1:surf_usm_v(l)%ns) )
           IF ( .NOT. ALLOCATED( t_surf_green_v_2(l)%t ) )                     &
              ALLOCATE ( t_surf_green_v_2(l)%t(1:surf_usm_v(l)%ns) )
           IF ( .NOT. ALLOCATED( t_green_v_1(l)%t ) )                          &           
              ALLOCATE ( t_green_v_1(l)%t(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns) ) 
           IF ( .NOT. ALLOCATED( t_green_v_2(l)%t ) )                          &           
              ALLOCATE ( t_green_v_2(l)%t(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns) )  
        ENDDO
!
!--     initial assignment of the pointers
        t_wall_v        => t_wall_v_1;        t_wall_v_p        => t_wall_v_2
        t_surf_wall_v   => t_surf_wall_v_1;   t_surf_wall_v_p   => t_surf_wall_v_2
        t_window_v      => t_window_v_1;      t_window_v_p      => t_window_v_2
        t_green_v       => t_green_v_1;       t_green_v_p       => t_green_v_2
        t_surf_window_v => t_surf_window_v_1; t_surf_window_v_p => t_surf_window_v_2
        t_surf_green_v  => t_surf_green_v_1;  t_surf_green_v_p  => t_surf_green_v_2

!
!--     Allocate intermediate timestep arrays. For horizontal surfaces.
        ALLOCATE ( surf_usm_h%tt_surface_wall_m(1:surf_usm_h%ns)               )
        ALLOCATE ( surf_usm_h%tt_wall_m(nzb_wall:nzt_wall+1,1:surf_usm_h%ns)   )
        ALLOCATE ( surf_usm_h%tt_surface_window_m(1:surf_usm_h%ns)             )
        ALLOCATE ( surf_usm_h%tt_window_m(nzb_wall:nzt_wall+1,1:surf_usm_h%ns) )
        ALLOCATE ( surf_usm_h%tt_green_m(nzb_wall:nzt_wall+1,1:surf_usm_h%ns)  )
        ALLOCATE ( surf_usm_h%tt_surface_green_m(1:surf_usm_h%ns)              )

!
!--    Allocate intermediate timestep arrays
!--    Horizontal surfaces
       ALLOCATE ( tm_liq_usm_h_m%var_usm_1d(1:surf_usm_h%ns)                   )
       tm_liq_usm_h_m%var_usm_1d = 0.0_wp
!
!--     Set inital values for prognostic quantities
        IF ( ALLOCATED( surf_usm_h%tt_surface_wall_m )   )  surf_usm_h%tt_surface_wall_m   = 0.0_wp
        IF ( ALLOCATED( surf_usm_h%tt_wall_m )           )  surf_usm_h%tt_wall_m           = 0.0_wp
        IF ( ALLOCATED( surf_usm_h%tt_surface_window_m ) )  surf_usm_h%tt_surface_window_m = 0.0_wp
        IF ( ALLOCATED( surf_usm_h%tt_window_m    )      )  surf_usm_h%tt_window_m         = 0.0_wp
        IF ( ALLOCATED( surf_usm_h%tt_green_m    )       )  surf_usm_h%tt_green_m          = 0.0_wp
        IF ( ALLOCATED( surf_usm_h%tt_surface_green_m )  )  surf_usm_h%tt_surface_green_m  = 0.0_wp
!
!--     Now, for vertical surfaces
        DO  l = 0, 3
           ALLOCATE ( surf_usm_v(l)%tt_surface_wall_m(1:surf_usm_v(l)%ns)               )
           ALLOCATE ( surf_usm_v(l)%tt_wall_m(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns)   )
           IF ( ALLOCATED( surf_usm_v(l)%tt_surface_wall_m ) )  surf_usm_v(l)%tt_surface_wall_m = 0.0_wp
           IF ( ALLOCATED( surf_usm_v(l)%tt_wall_m    ) )  surf_usm_v(l)%tt_wall_m    = 0.0_wp
           ALLOCATE ( surf_usm_v(l)%tt_surface_window_m(1:surf_usm_v(l)%ns)             )
           ALLOCATE ( surf_usm_v(l)%tt_window_m(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns) )
           IF ( ALLOCATED( surf_usm_v(l)%tt_surface_window_m ) )  surf_usm_v(l)%tt_surface_window_m = 0.0_wp
           IF ( ALLOCATED( surf_usm_v(l)%tt_window_m  ) )  surf_usm_v(l)%tt_window_m    = 0.0_wp
           ALLOCATE ( surf_usm_v(l)%tt_surface_green_m(1:surf_usm_v(l)%ns)              )
           IF ( ALLOCATED( surf_usm_v(l)%tt_surface_green_m ) )  surf_usm_v(l)%tt_surface_green_m = 0.0_wp
           ALLOCATE ( surf_usm_v(l)%tt_green_m(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns)  )
           IF ( ALLOCATED( surf_usm_v(l)%tt_green_m   ) )  surf_usm_v(l)%tt_green_m    = 0.0_wp
        ENDDO
!
!--     allocate wall heat flux output array and set initial values. For horizontal surfaces
!        ALLOCATE ( surf_usm_h%wshf(1:surf_usm_h%ns)    )  !can be removed 
        ALLOCATE ( surf_usm_h%wshf_eb(1:surf_usm_h%ns) )
        ALLOCATE ( surf_usm_h%wghf_eb(1:surf_usm_h%ns) )
        ALLOCATE ( surf_usm_h%wghf_eb_window(1:surf_usm_h%ns) )
        ALLOCATE ( surf_usm_h%wghf_eb_green(1:surf_usm_h%ns) )
        ALLOCATE ( surf_usm_h%iwghf_eb(1:surf_usm_h%ns) )
        ALLOCATE ( surf_usm_h%iwghf_eb_window(1:surf_usm_h%ns) )
        IF ( ALLOCATED( surf_usm_h%wshf    ) )  surf_usm_h%wshf    = 0.0_wp
        IF ( ALLOCATED( surf_usm_h%wshf_eb ) )  surf_usm_h%wshf_eb = 0.0_wp
        IF ( ALLOCATED( surf_usm_h%wghf_eb ) )  surf_usm_h%wghf_eb = 0.0_wp
        IF ( ALLOCATED( surf_usm_h%wghf_eb_window ) )  surf_usm_h%wghf_eb_window = 0.0_wp
        IF ( ALLOCATED( surf_usm_h%wghf_eb_green ) )  surf_usm_h%wghf_eb_green = 0.0_wp
        IF ( ALLOCATED( surf_usm_h%iwghf_eb ) )  surf_usm_h%iwghf_eb = 0.0_wp
        IF ( ALLOCATED( surf_usm_h%iwghf_eb_window ) )  surf_usm_h%iwghf_eb_window = 0.0_wp
!
!--     Now, for vertical surfaces
        DO  l = 0, 3
!           ALLOCATE ( surf_usm_v(l)%wshf(1:surf_usm_v(l)%ns)    )    ! can be removed
           ALLOCATE ( surf_usm_v(l)%wshf_eb(1:surf_usm_v(l)%ns) )
           ALLOCATE ( surf_usm_v(l)%wghf_eb(1:surf_usm_v(l)%ns) )
           ALLOCATE ( surf_usm_v(l)%wghf_eb_window(1:surf_usm_v(l)%ns) )
           ALLOCATE ( surf_usm_v(l)%wghf_eb_green(1:surf_usm_v(l)%ns) )
           ALLOCATE ( surf_usm_v(l)%iwghf_eb(1:surf_usm_v(l)%ns) )
           ALLOCATE ( surf_usm_v(l)%iwghf_eb_window(1:surf_usm_v(l)%ns) )
           IF ( ALLOCATED( surf_usm_v(l)%wshf    ) )  surf_usm_v(l)%wshf    = 0.0_wp
           IF ( ALLOCATED( surf_usm_v(l)%wshf_eb ) )  surf_usm_v(l)%wshf_eb = 0.0_wp
           IF ( ALLOCATED( surf_usm_v(l)%wghf_eb ) )  surf_usm_v(l)%wghf_eb = 0.0_wp
           IF ( ALLOCATED( surf_usm_v(l)%wghf_eb_window ) )  surf_usm_v(l)%wghf_eb_window = 0.0_wp
           IF ( ALLOCATED( surf_usm_v(l)%wghf_eb_green ) )  surf_usm_v(l)%wghf_eb_green = 0.0_wp
           IF ( ALLOCATED( surf_usm_v(l)%iwghf_eb ) )  surf_usm_v(l)%iwghf_eb = 0.0_wp
           IF ( ALLOCATED( surf_usm_v(l)%iwghf_eb_window ) )  surf_usm_v(l)%iwghf_eb_window = 0.0_wp
        ENDDO
!
!--     Initialize building-surface properties, which are also required by other modules, 
!--     e.g. the indoor model.
        CALL usm_define_pars
        
        IF ( debug_output )  CALL debug_message( 'usm_init_arrays', 'end' )
        
    END SUBROUTINE usm_init_arrays


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Sum up and time-average urban surface output quantities as well as allocate
!> the array necessary for storing the average.
!------------------------------------------------------------------------------!
    SUBROUTINE usm_3d_data_averaging( mode, variable )

        IMPLICIT NONE

        CHARACTER(LEN=*), INTENT(IN) ::  mode
        CHARACTER(LEN=*), INTENT(IN) :: variable
  
        INTEGER(iwp)                                       :: i, j, k, l, m, ids, idsint, iwl, istat  !< runnin indices
        CHARACTER(LEN=varnamelength)                       :: var                                     !< trimmed variable
        INTEGER(iwp), PARAMETER                            :: nd = 5                                  !< number of directions
        CHARACTER(LEN=6), DIMENSION(0:nd-1), PARAMETER     :: dirname = (/ '_roof ', '_south', '_north', '_west ', '_east ' /)
        INTEGER(iwp), DIMENSION(0:nd-1), PARAMETER         :: dirint = (/ iup_u, isouth_u, inorth_u, iwest_u, ieast_u /)

        IF ( variable(1:4) == 'usm_' )  THEN  ! is such a check really rquired?

!
!--     find the real name of the variable
        ids = -1
        l = -1
        var = TRIM(variable)
        DO i = 0, nd-1
            k = len(TRIM(var))
            j = len(TRIM(dirname(i)))
            IF ( TRIM(var(k-j+1:k)) == TRIM(dirname(i)) )  THEN
                ids = i
                idsint = dirint(ids)
                var = var(:k-j)
                EXIT
            ENDIF
        ENDDO
        l = idsint - 2  ! horisontal direction index - terible hack !
        IF ( l < 0 .OR. l > 3 ) THEN
           l = -1
        END IF
        IF ( ids == -1 )  THEN
            var = TRIM(variable)
        ENDIF
        IF ( var(1:11) == 'usm_t_wall_'  .AND.  len(TRIM(var)) >= 12 )  THEN
!
!--          wall layers
            READ(var(12:12), '(I1)', iostat=istat ) iwl
            IF ( istat == 0  .AND.  iwl >= nzb_wall  .AND.  iwl <= nzt_wall )  THEN
                var = var(1:10)
            ELSE
!
!--             wrong wall layer index
                RETURN
            ENDIF
        ENDIF
        IF ( var(1:13) == 'usm_t_window_'  .AND.  len(TRIM(var)) >= 14 )  THEN
!
!--          wall layers
            READ(var(14:14), '(I1)', iostat=istat ) iwl
            IF ( istat == 0  .AND.  iwl >= nzb_wall  .AND.  iwl <= nzt_wall )  THEN
                var = var(1:12)
            ELSE
!
!--             wrong window layer index
                RETURN
            ENDIF
        ENDIF
        IF ( var(1:12) == 'usm_t_green_'  .AND.  len(TRIM(var)) >= 13 )  THEN
!
!--          wall layers
            READ(var(13:13), '(I1)', iostat=istat ) iwl
            IF ( istat == 0  .AND.  iwl >= nzb_wall  .AND.  iwl <= nzt_wall )  THEN
                var = var(1:11)
            ELSE
!
!--             wrong green layer index
                RETURN
            ENDIF
        ENDIF
        IF ( var(1:8) == 'usm_swc_'  .AND.  len(TRIM(var)) >= 9 )  THEN
!
!--          swc layers
            READ(var(9:9), '(I1)', iostat=istat ) iwl
            IF ( istat == 0  .AND.  iwl >= nzb_wall  .AND.  iwl <= nzt_wall )  THEN
                var = var(1:7)
            ELSE
!
!--             wrong swc layer index
                RETURN
            ENDIF
        ENDIF

        IF ( mode == 'allocate' )  THEN
           
           SELECT CASE ( TRIM( var ) )

                CASE ( 'usm_wshf' )
!
!--                 array of sensible heat flux from surfaces
!--                 land surfaces
                    IF ( l == -1 ) THEN
                       IF ( .NOT.  ALLOCATED(surf_usm_h%wshf_eb_av) )  THEN
                          ALLOCATE ( surf_usm_h%wshf_eb_av(1:surf_usm_h%ns) )
                          surf_usm_h%wshf_eb_av = 0.0_wp
                       ENDIF
                    ELSE
                       IF ( .NOT.  ALLOCATED(surf_usm_v(l)%wshf_eb_av) )  THEN
                           ALLOCATE ( surf_usm_v(l)%wshf_eb_av(1:surf_usm_v(l)%ns) )
                           surf_usm_v(l)%wshf_eb_av = 0.0_wp
                       ENDIF
                    ENDIF
                    
                CASE ( 'usm_qsws' )
!
!--                 array of latent heat flux from surfaces
!--                 land surfaces
                    IF ( l == -1 .AND. .NOT.  ALLOCATED(surf_usm_h%qsws_av) )  THEN
                        ALLOCATE ( surf_usm_h%qsws_av(1:surf_usm_h%ns) )
                        surf_usm_h%qsws_av = 0.0_wp
                    ELSE
                       IF ( .NOT.  ALLOCATED(surf_usm_v(l)%qsws_av) )  THEN
                           ALLOCATE ( surf_usm_v(l)%qsws_av(1:surf_usm_v(l)%ns) )
                           surf_usm_v(l)%qsws_av = 0.0_wp
                       ENDIF
                    ENDIF
                    
                CASE ( 'usm_qsws_veg' )
!
!--                 array of latent heat flux from vegetation surfaces
!--                 land surfaces
                    IF ( l == -1 .AND. .NOT.  ALLOCATED(surf_usm_h%qsws_veg_av) )  THEN
                        ALLOCATE ( surf_usm_h%qsws_veg_av(1:surf_usm_h%ns) )
                        surf_usm_h%qsws_veg_av = 0.0_wp
                    ELSE
                       IF ( .NOT.  ALLOCATED(surf_usm_v(l)%qsws_veg_av) )  THEN
                           ALLOCATE ( surf_usm_v(l)%qsws_veg_av(1:surf_usm_v(l)%ns) )
                           surf_usm_v(l)%qsws_veg_av = 0.0_wp
                       ENDIF
                    ENDIF
                    
                CASE ( 'usm_qsws_liq' )
!
!--                 array of latent heat flux from surfaces with liquid
!--                 land surfaces
                    IF ( l == -1 .AND. .NOT.  ALLOCATED(surf_usm_h%qsws_liq_av) )  THEN
                        ALLOCATE ( surf_usm_h%qsws_liq_av(1:surf_usm_h%ns) )
                        surf_usm_h%qsws_liq_av = 0.0_wp
                    ELSE
                       IF ( .NOT.  ALLOCATED(surf_usm_v(l)%qsws_liq_av) )  THEN
                           ALLOCATE ( surf_usm_v(l)%qsws_liq_av(1:surf_usm_v(l)%ns) )
                           surf_usm_v(l)%qsws_liq_av = 0.0_wp
                       ENDIF
                    ENDIF
!
!--             Please note, the following output quantities belongs to the 
!--             individual tile fractions - ground heat flux at wall-, window-, 
!--             and green fraction. Aggregated ground-heat flux is treated
!--             accordingly in average_3d_data, sum_up_3d_data, etc..
                CASE ( 'usm_wghf' )
!
!--                 array of heat flux from ground (wall, roof, land)
                    IF ( l == -1 ) THEN
                       IF ( .NOT.  ALLOCATED(surf_usm_h%wghf_eb_av) )  THEN
                           ALLOCATE ( surf_usm_h%wghf_eb_av(1:surf_usm_h%ns) )
                           surf_usm_h%wghf_eb_av = 0.0_wp
                       ENDIF
                    ELSE
                       IF ( .NOT.  ALLOCATED(surf_usm_v(l)%wghf_eb_av) )  THEN
                           ALLOCATE ( surf_usm_v(l)%wghf_eb_av(1:surf_usm_v(l)%ns) )
                           surf_usm_v(l)%wghf_eb_av = 0.0_wp
                       ENDIF
                    ENDIF

                CASE ( 'usm_wghf_window' )
!
!--                 array of heat flux from window ground (wall, roof, land)
                    IF ( l == -1 ) THEN
                       IF ( .NOT.  ALLOCATED(surf_usm_h%wghf_eb_window_av) )  THEN
                           ALLOCATE ( surf_usm_h%wghf_eb_window_av(1:surf_usm_h%ns) )
                           surf_usm_h%wghf_eb_window_av = 0.0_wp
                       ENDIF
                    ELSE
                       IF ( .NOT.  ALLOCATED(surf_usm_v(l)%wghf_eb_window_av) )  THEN
                           ALLOCATE ( surf_usm_v(l)%wghf_eb_window_av(1:surf_usm_v(l)%ns) )
                           surf_usm_v(l)%wghf_eb_window_av = 0.0_wp
                       ENDIF
                    ENDIF

                CASE ( 'usm_wghf_green' )
!
!--                 array of heat flux from green ground (wall, roof, land)
                    IF ( l == -1 ) THEN
                       IF ( .NOT.  ALLOCATED(surf_usm_h%wghf_eb_green_av) )  THEN
                           ALLOCATE ( surf_usm_h%wghf_eb_green_av(1:surf_usm_h%ns) )
                           surf_usm_h%wghf_eb_green_av = 0.0_wp
                       ENDIF
                    ELSE
                       IF ( .NOT.  ALLOCATED(surf_usm_v(l)%wghf_eb_green_av) )  THEN
                           ALLOCATE ( surf_usm_v(l)%wghf_eb_green_av(1:surf_usm_v(l)%ns) )
                           surf_usm_v(l)%wghf_eb_green_av = 0.0_wp
                       ENDIF
                    ENDIF

                CASE ( 'usm_iwghf' )
!
!--                 array of heat flux from indoor ground (wall, roof, land)
                    IF ( l == -1 ) THEN
                       IF ( .NOT.  ALLOCATED(surf_usm_h%iwghf_eb_av) )  THEN
                           ALLOCATE ( surf_usm_h%iwghf_eb_av(1:surf_usm_h%ns) )
                           surf_usm_h%iwghf_eb_av = 0.0_wp
                       ENDIF
                    ELSE
                       IF ( .NOT.  ALLOCATED(surf_usm_v(l)%iwghf_eb_av) )  THEN
                           ALLOCATE ( surf_usm_v(l)%iwghf_eb_av(1:surf_usm_v(l)%ns) )
                           surf_usm_v(l)%iwghf_eb_av = 0.0_wp
                       ENDIF
                    ENDIF

                CASE ( 'usm_iwghf_window' )
!
!--                 array of heat flux from indoor window ground (wall, roof, land)
                    IF ( l == -1 ) THEN
                       IF ( .NOT.  ALLOCATED(surf_usm_h%iwghf_eb_window_av) )  THEN
                           ALLOCATE ( surf_usm_h%iwghf_eb_window_av(1:surf_usm_h%ns) )
                           surf_usm_h%iwghf_eb_window_av = 0.0_wp
                       ENDIF
                    ELSE
                       IF ( .NOT.  ALLOCATED(surf_usm_v(l)%iwghf_eb_window_av) )  THEN
                           ALLOCATE ( surf_usm_v(l)%iwghf_eb_window_av(1:surf_usm_v(l)%ns) )
                           surf_usm_v(l)%iwghf_eb_window_av = 0.0_wp
                       ENDIF
                    ENDIF

                CASE ( 'usm_t_surf_wall' )
!
!--                 surface temperature for surfaces
                    IF ( l == -1 ) THEN
                       IF ( .NOT.  ALLOCATED(surf_usm_h%t_surf_wall_av) )  THEN
                           ALLOCATE ( surf_usm_h%t_surf_wall_av(1:surf_usm_h%ns) )
                           surf_usm_h%t_surf_wall_av = 0.0_wp
                       ENDIF
                    ELSE
                       IF ( .NOT.  ALLOCATED(surf_usm_v(l)%t_surf_wall_av) )  THEN
                           ALLOCATE ( surf_usm_v(l)%t_surf_wall_av(1:surf_usm_v(l)%ns) )
                           surf_usm_v(l)%t_surf_wall_av = 0.0_wp
                       ENDIF
                    ENDIF

                CASE ( 'usm_t_surf_window' )
!
!--                 surface temperature for window surfaces
                    IF ( l == -1 ) THEN
                       IF ( .NOT.  ALLOCATED(surf_usm_h%t_surf_window_av) )  THEN
                           ALLOCATE ( surf_usm_h%t_surf_window_av(1:surf_usm_h%ns) )
                           surf_usm_h%t_surf_window_av = 0.0_wp
                       ENDIF
                    ELSE
                       IF ( .NOT.  ALLOCATED(surf_usm_v(l)%t_surf_window_av) )  THEN
                           ALLOCATE ( surf_usm_v(l)%t_surf_window_av(1:surf_usm_v(l)%ns) )
                           surf_usm_v(l)%t_surf_window_av = 0.0_wp
                       ENDIF
                    ENDIF
                    
                CASE ( 'usm_t_surf_green' )
!
!--                 surface temperature for green surfaces
                    IF ( l == -1 ) THEN
                       IF ( .NOT.  ALLOCATED(surf_usm_h%t_surf_green_av) )  THEN
                           ALLOCATE ( surf_usm_h%t_surf_green_av(1:surf_usm_h%ns) )
                           surf_usm_h%t_surf_green_av = 0.0_wp
                       ENDIF
                    ELSE
                       IF ( .NOT.  ALLOCATED(surf_usm_v(l)%t_surf_green_av) )  THEN
                           ALLOCATE ( surf_usm_v(l)%t_surf_green_av(1:surf_usm_v(l)%ns) )
                           surf_usm_v(l)%t_surf_green_av = 0.0_wp
                       ENDIF
                    ENDIF
                
                CASE ( 'usm_theta_10cm' )
!
!--                 near surface (10cm) temperature for whole surfaces
                    IF ( l == -1 ) THEN
                       IF ( .NOT.  ALLOCATED(surf_usm_h%pt_10cm_av) )  THEN
                           ALLOCATE ( surf_usm_h%pt_10cm_av(1:surf_usm_h%ns) )
                           surf_usm_h%pt_10cm_av = 0.0_wp
                       ENDIF
                    ELSE
                       IF ( .NOT.  ALLOCATED(surf_usm_v(l)%pt_10cm_av) )  THEN
                           ALLOCATE ( surf_usm_v(l)%pt_10cm_av(1:surf_usm_v(l)%ns) )
                           surf_usm_v(l)%pt_10cm_av = 0.0_wp
                       ENDIF
                    ENDIF
                 
                CASE ( 'usm_t_wall' )
!
!--                 wall temperature for iwl layer of walls and land
                    IF ( l == -1 ) THEN
                       IF ( .NOT.  ALLOCATED(surf_usm_h%t_wall_av) )  THEN
                           ALLOCATE ( surf_usm_h%t_wall_av(nzb_wall:nzt_wall,1:surf_usm_h%ns) )
                           surf_usm_h%t_wall_av = 0.0_wp
                       ENDIF
                    ELSE
                       IF ( .NOT.  ALLOCATED(surf_usm_v(l)%t_wall_av) )  THEN
                           ALLOCATE ( surf_usm_v(l)%t_wall_av(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns) )
                           surf_usm_v(l)%t_wall_av = 0.0_wp
                       ENDIF
                    ENDIF

                CASE ( 'usm_t_window' )
!
!--                 window temperature for iwl layer of walls and land
                    IF ( l == -1 ) THEN
                       IF ( .NOT.  ALLOCATED(surf_usm_h%t_window_av) )  THEN
                           ALLOCATE ( surf_usm_h%t_window_av(nzb_wall:nzt_wall,1:surf_usm_h%ns) )
                           surf_usm_h%t_window_av = 0.0_wp
                       ENDIF
                    ELSE
                       IF ( .NOT.  ALLOCATED(surf_usm_v(l)%t_window_av) )  THEN
                           ALLOCATE ( surf_usm_v(l)%t_window_av(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns) )
                           surf_usm_v(l)%t_window_av = 0.0_wp
                       ENDIF
                    ENDIF

                CASE ( 'usm_t_green' )
!
!--                 green temperature for iwl layer of walls and land
                    IF ( l == -1 ) THEN
                       IF ( .NOT.  ALLOCATED(surf_usm_h%t_green_av) )  THEN
                           ALLOCATE ( surf_usm_h%t_green_av(nzb_wall:nzt_wall,1:surf_usm_h%ns) )
                           surf_usm_h%t_green_av = 0.0_wp
                       ENDIF
                    ELSE
                       IF ( .NOT.  ALLOCATED(surf_usm_v(l)%t_green_av) )  THEN
                           ALLOCATE ( surf_usm_v(l)%t_green_av(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns) )
                           surf_usm_v(l)%t_green_av = 0.0_wp
                       ENDIF
                    ENDIF
                CASE ( 'usm_swc' )
!
!--                 soil water content for iwl layer of walls and land
                    IF ( l == -1 .AND. .NOT.  ALLOCATED(surf_usm_h%swc_av) )  THEN
                        ALLOCATE ( surf_usm_h%swc_av(nzb_wall:nzt_wall,1:surf_usm_h%ns) )
                        surf_usm_h%swc_av = 0.0_wp
                    ELSE
                       IF ( .NOT.  ALLOCATED(surf_usm_v(l)%swc_av) )  THEN
                           ALLOCATE ( surf_usm_v(l)%swc_av(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns) )
                           surf_usm_v(l)%swc_av = 0.0_wp
                       ENDIF
                    ENDIF

               CASE DEFAULT
                   CONTINUE

           END SELECT

        ELSEIF ( mode == 'sum' )  THEN
           
           SELECT CASE ( TRIM( var ) )

                CASE ( 'usm_wshf' )
!
!--                 array of sensible heat flux from surfaces (land, roof, wall)
                    IF ( l == -1 ) THEN
                       DO  m = 1, surf_usm_h%ns
                          surf_usm_h%wshf_eb_av(m) =                              &
                                             surf_usm_h%wshf_eb_av(m) +           &
                                             surf_usm_h%wshf_eb(m)
                       ENDDO
                    ELSE
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%wshf_eb_av(m) =                        &
                                          surf_usm_v(l)%wshf_eb_av(m) +        &
                                          surf_usm_v(l)%wshf_eb(m)
                       ENDDO
                    ENDIF
                    
                CASE ( 'usm_qsws' )
!
!--                 array of latent heat flux from surfaces (land, roof, wall)
                    IF ( l == -1 ) THEN
                    DO  m = 1, surf_usm_h%ns
                       surf_usm_h%qsws_av(m) =                              &
                                          surf_usm_h%qsws_av(m) +           &
                                          surf_usm_h%qsws(m) * l_v
                    ENDDO
                    ELSE
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%qsws_av(m) =                        &
                                          surf_usm_v(l)%qsws_av(m) +        &
                                          surf_usm_v(l)%qsws(m) * l_v
                       ENDDO
                    ENDIF
                    
                CASE ( 'usm_qsws_veg' )
!
!--                 array of latent heat flux from vegetation surfaces (land, roof, wall)
                    IF ( l == -1 ) THEN
                    DO  m = 1, surf_usm_h%ns
                       surf_usm_h%qsws_veg_av(m) =                              &
                                          surf_usm_h%qsws_veg_av(m) +           &
                                          surf_usm_h%qsws_veg(m)
                    ENDDO
                    ELSE
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%qsws_veg_av(m) =                        &
                                          surf_usm_v(l)%qsws_veg_av(m) +        &
                                          surf_usm_v(l)%qsws_veg(m)
                       ENDDO
                    ENDIF
                    
                CASE ( 'usm_qsws_liq' )
!
!--                 array of latent heat flux from surfaces with liquid (land, roof, wall)
                    IF ( l == -1 ) THEN
                    DO  m = 1, surf_usm_h%ns
                       surf_usm_h%qsws_liq_av(m) =                              &
                                          surf_usm_h%qsws_liq_av(m) +           &
                                          surf_usm_h%qsws_liq(m)
                    ENDDO
                    ELSE
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%qsws_liq_av(m) =                        &
                                          surf_usm_v(l)%qsws_liq_av(m) +        &
                                          surf_usm_v(l)%qsws_liq(m)
                       ENDDO
                    ENDIF
                    
                CASE ( 'usm_wghf' )
!
!--                 array of heat flux from ground (wall, roof, land)
                    IF ( l == -1 ) THEN
                       DO  m = 1, surf_usm_h%ns
                          surf_usm_h%wghf_eb_av(m) =                              &
                                             surf_usm_h%wghf_eb_av(m) +           &
                                             surf_usm_h%wghf_eb(m)
                       ENDDO
                    ELSE
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%wghf_eb_av(m) =                        &
                                          surf_usm_v(l)%wghf_eb_av(m) +        &
                                          surf_usm_v(l)%wghf_eb(m)
                       ENDDO
                    ENDIF
                    
                CASE ( 'usm_wghf_window' )
!
!--                 array of heat flux from window ground (wall, roof, land)
                    IF ( l == -1 ) THEN
                       DO  m = 1, surf_usm_h%ns
                          surf_usm_h%wghf_eb_window_av(m) =                              &
                                             surf_usm_h%wghf_eb_window_av(m) +           &
                                             surf_usm_h%wghf_eb_window(m)
                       ENDDO
                    ELSE
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%wghf_eb_window_av(m) =                        &
                                          surf_usm_v(l)%wghf_eb_window_av(m) +        &
                                          surf_usm_v(l)%wghf_eb_window(m)
                       ENDDO
                    ENDIF

                CASE ( 'usm_wghf_green' )
!
!--                 array of heat flux from green ground (wall, roof, land)
                    IF ( l == -1 ) THEN
                       DO  m = 1, surf_usm_h%ns
                          surf_usm_h%wghf_eb_green_av(m) =                              &
                                             surf_usm_h%wghf_eb_green_av(m) +           &
                                             surf_usm_h%wghf_eb_green(m)
                       ENDDO
                    ELSE
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%wghf_eb_green_av(m) =                        &
                                          surf_usm_v(l)%wghf_eb_green_av(m) +        &
                                          surf_usm_v(l)%wghf_eb_green(m)
                       ENDDO
                    ENDIF
                    
                CASE ( 'usm_iwghf' )
!
!--                 array of heat flux from indoor ground (wall, roof, land)
                    IF ( l == -1 ) THEN
                       DO  m = 1, surf_usm_h%ns
                          surf_usm_h%iwghf_eb_av(m) =                              &
                                             surf_usm_h%iwghf_eb_av(m) +           &
                                             surf_usm_h%iwghf_eb(m)
                       ENDDO
                    ELSE
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%iwghf_eb_av(m) =                        &
                                          surf_usm_v(l)%iwghf_eb_av(m) +        &
                                          surf_usm_v(l)%iwghf_eb(m)
                       ENDDO
                    ENDIF
                    
                CASE ( 'usm_iwghf_window' )
!
!--                 array of heat flux from indoor window ground (wall, roof, land)
                    IF ( l == -1 ) THEN
                       DO  m = 1, surf_usm_h%ns
                          surf_usm_h%iwghf_eb_window_av(m) =                              &
                                             surf_usm_h%iwghf_eb_window_av(m) +           &
                                             surf_usm_h%iwghf_eb_window(m)
                       ENDDO
                    ELSE
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%iwghf_eb_window_av(m) =                        &
                                          surf_usm_v(l)%iwghf_eb_window_av(m) +        &
                                          surf_usm_v(l)%iwghf_eb_window(m)
                       ENDDO
                    ENDIF
                    
                CASE ( 'usm_t_surf_wall' )
!
!--                 surface temperature for surfaces
                    IF ( l == -1 ) THEN
                       DO  m = 1, surf_usm_h%ns
                       surf_usm_h%t_surf_wall_av(m) =                               & 
                                          surf_usm_h%t_surf_wall_av(m) +            &
                                          t_surf_wall_h(m)
                       ENDDO
                    ELSE
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%t_surf_wall_av(m) =                         &
                                          surf_usm_v(l)%t_surf_wall_av(m) +         &
                                          t_surf_wall_v(l)%t(m)
                       ENDDO
                    ENDIF
                    
                CASE ( 'usm_t_surf_window' )
!
!--                 surface temperature for window surfaces
                    IF ( l == -1 ) THEN
                       DO  m = 1, surf_usm_h%ns
                          surf_usm_h%t_surf_window_av(m) =                               &
                                             surf_usm_h%t_surf_window_av(m) +            &
                                             t_surf_window_h(m)
                       ENDDO
                    ELSE
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%t_surf_window_av(m) =                         &
                                          surf_usm_v(l)%t_surf_window_av(m) +         &
                                          t_surf_window_v(l)%t(m)
                       ENDDO
                    ENDIF
                    
                CASE ( 'usm_t_surf_green' )
!
!--                 surface temperature for green surfaces
                    IF ( l == -1 ) THEN
                       DO  m = 1, surf_usm_h%ns
                          surf_usm_h%t_surf_green_av(m) =                               &
                                             surf_usm_h%t_surf_green_av(m) +            &
                                             t_surf_green_h(m)
                       ENDDO
                    ELSE
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%t_surf_green_av(m) =                         &
                                          surf_usm_v(l)%t_surf_green_av(m) +         &
                                          t_surf_green_v(l)%t(m)
                       ENDDO
                    ENDIF
                
                CASE ( 'usm_theta_10cm' )
!
!--                 near surface temperature for whole surfaces
                    IF ( l == -1 ) THEN
                       DO  m = 1, surf_usm_h%ns
                          surf_usm_h%pt_10cm_av(m) =                               &
                                             surf_usm_h%pt_10cm_av(m) +            &
                                             surf_usm_h%pt_10cm(m)
                       ENDDO
                    ELSE
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%pt_10cm_av(m) =                         &
                                          surf_usm_v(l)%pt_10cm_av(m) +         &
                                          surf_usm_v(l)%pt_10cm(m)
                       ENDDO
                    ENDIF
                    
                CASE ( 'usm_t_wall' )
!
!--                 wall temperature for  iwl layer of walls and land
                    IF ( l == -1 ) THEN
                       DO  m = 1, surf_usm_h%ns
                          surf_usm_h%t_wall_av(iwl,m) =                           &
                                             surf_usm_h%t_wall_av(iwl,m) +        &
                                             t_wall_h(iwl,m)
                       ENDDO
                    ELSE
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%t_wall_av(iwl,m) =                     &
                                          surf_usm_v(l)%t_wall_av(iwl,m) +     &
                                          t_wall_v(l)%t(iwl,m)
                       ENDDO
                    ENDIF
                    
                CASE ( 'usm_t_window' )
!
!--                 window temperature for  iwl layer of walls and land
                    IF ( l == -1 ) THEN
                       DO  m = 1, surf_usm_h%ns
                          surf_usm_h%t_window_av(iwl,m) =                           &
                                             surf_usm_h%t_window_av(iwl,m) +        &
                                             t_window_h(iwl,m)
                       ENDDO
                    ELSE
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%t_window_av(iwl,m) =                     &
                                          surf_usm_v(l)%t_window_av(iwl,m) +     &
                                          t_window_v(l)%t(iwl,m)
                       ENDDO
                    ENDIF

                CASE ( 'usm_t_green' )
!
!--                 green temperature for  iwl layer of walls and land
                    IF ( l == -1 ) THEN
                       DO  m = 1, surf_usm_h%ns
                          surf_usm_h%t_green_av(iwl,m) =                           &
                                             surf_usm_h%t_green_av(iwl,m) +        &
                                             t_green_h(iwl,m)
                       ENDDO
                    ELSE
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%t_green_av(iwl,m) =                     &
                                          surf_usm_v(l)%t_green_av(iwl,m) +     &
                                          t_green_v(l)%t(iwl,m)
                       ENDDO
                    ENDIF

                CASE ( 'usm_swc' )
!
!--                 soil water content for  iwl layer of walls and land
                    IF ( l == -1 ) THEN
                       DO  m = 1, surf_usm_h%ns
                          surf_usm_h%swc_av(iwl,m) =                           &
                                          surf_usm_h%swc_av(iwl,m) +           &
                                             swc_h(iwl,m)
                       ENDDO
                    ELSE
                    ENDIF

                CASE DEFAULT
                    CONTINUE

           END SELECT

        ELSEIF ( mode == 'average' )  THEN
           
           SELECT CASE ( TRIM( var ) )

                CASE ( 'usm_wshf' )
!
!--                 array of sensible heat flux from surfaces (land, roof, wall)
                    IF ( l == -1 ) THEN
                       DO  m = 1, surf_usm_h%ns
                          surf_usm_h%wshf_eb_av(m) =                              &
                                             surf_usm_h%wshf_eb_av(m) /           &
                                             REAL( average_count_3d, kind=wp )
                       ENDDO
                    ELSE
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%wshf_eb_av(m) =                        &
                                          surf_usm_v(l)%wshf_eb_av(m) /        &
                                          REAL( average_count_3d, kind=wp )
                       ENDDO
                    ENDIF
                    
                CASE ( 'usm_qsws' )
!
!--                 array of latent heat flux from surfaces (land, roof, wall)
                    IF ( l == -1 ) THEN
                    DO  m = 1, surf_usm_h%ns
                       surf_usm_h%qsws_av(m) =                              &
                                          surf_usm_h%qsws_av(m) /           &
                                          REAL( average_count_3d, kind=wp )
                    ENDDO
                    ELSE
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%qsws_av(m) =                        &
                                          surf_usm_v(l)%qsws_av(m) /        &
                                          REAL( average_count_3d, kind=wp )
                       ENDDO
                    ENDIF

                CASE ( 'usm_qsws_veg' )
!
!--                 array of latent heat flux from vegetation surfaces (land, roof, wall)
                    IF ( l == -1 ) THEN
                    DO  m = 1, surf_usm_h%ns
                       surf_usm_h%qsws_veg_av(m) =                              &
                                          surf_usm_h%qsws_veg_av(m) /           &
                                          REAL( average_count_3d, kind=wp )
                    ENDDO
                    ELSE
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%qsws_veg_av(m) =                        &
                                          surf_usm_v(l)%qsws_veg_av(m) /        &
                                          REAL( average_count_3d, kind=wp )
                       ENDDO
                    ENDIF
                    
                CASE ( 'usm_qsws_liq' )
!
!--                 array of latent heat flux from surfaces with liquid (land, roof, wall)
                    IF ( l == -1 ) THEN
                    DO  m = 1, surf_usm_h%ns
                       surf_usm_h%qsws_liq_av(m) =                              &
                                          surf_usm_h%qsws_liq_av(m) /           &
                                          REAL( average_count_3d, kind=wp )
                    ENDDO
                    ELSE
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%qsws_liq_av(m) =                        &
                                          surf_usm_v(l)%qsws_liq_av(m) /        &
                                          REAL( average_count_3d, kind=wp )
                       ENDDO
                    ENDIF
                    
                CASE ( 'usm_wghf' )
!
!--                 array of heat flux from ground (wall, roof, land)
                    IF ( l == -1 ) THEN
                       DO  m = 1, surf_usm_h%ns
                          surf_usm_h%wghf_eb_av(m) =                              &
                                             surf_usm_h%wghf_eb_av(m) /           &
                                             REAL( average_count_3d, kind=wp )
                       ENDDO
                    ELSE
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%wghf_eb_av(m) =                        &
                                          surf_usm_v(l)%wghf_eb_av(m) /        &
                                          REAL( average_count_3d, kind=wp )
                       ENDDO
                    ENDIF
                    
                CASE ( 'usm_wghf_window' )
!
!--                 array of heat flux from window ground (wall, roof, land)
                    IF ( l == -1 ) THEN
                       DO  m = 1, surf_usm_h%ns
                          surf_usm_h%wghf_eb_window_av(m) =                              &
                                             surf_usm_h%wghf_eb_window_av(m) /           &
                                             REAL( average_count_3d, kind=wp )
                       ENDDO
                    ELSE
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%wghf_eb_window_av(m) =                        &
                                          surf_usm_v(l)%wghf_eb_window_av(m) /        &
                                          REAL( average_count_3d, kind=wp )
                       ENDDO
                    ENDIF

                CASE ( 'usm_wghf_green' )
!
!--                 array of heat flux from green ground (wall, roof, land)
                    IF ( l == -1 ) THEN
                       DO  m = 1, surf_usm_h%ns
                          surf_usm_h%wghf_eb_green_av(m) =                              &
                                             surf_usm_h%wghf_eb_green_av(m) /           &
                                             REAL( average_count_3d, kind=wp )
                       ENDDO
                    ELSE
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%wghf_eb_green_av(m) =                        &
                                          surf_usm_v(l)%wghf_eb_green_av(m) /        &
                                          REAL( average_count_3d, kind=wp )
                       ENDDO
                    ENDIF

                CASE ( 'usm_iwghf' )
!
!--                 array of heat flux from indoor ground (wall, roof, land)
                    IF ( l == -1 ) THEN
                       DO  m = 1, surf_usm_h%ns
                          surf_usm_h%iwghf_eb_av(m) =                              &
                                             surf_usm_h%iwghf_eb_av(m) /           &
                                             REAL( average_count_3d, kind=wp )
                       ENDDO
                    ELSE
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%iwghf_eb_av(m) =                        &
                                          surf_usm_v(l)%iwghf_eb_av(m) /        &
                                          REAL( average_count_3d, kind=wp )
                       ENDDO
                    ENDIF
                    
                CASE ( 'usm_iwghf_window' )
!
!--                 array of heat flux from indoor window ground (wall, roof, land)
                    IF ( l == -1 ) THEN
                       DO  m = 1, surf_usm_h%ns
                          surf_usm_h%iwghf_eb_window_av(m) =                              &
                                             surf_usm_h%iwghf_eb_window_av(m) /           &
                                             REAL( average_count_3d, kind=wp )
                       ENDDO
                    ELSE
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%iwghf_eb_window_av(m) =                        &
                                          surf_usm_v(l)%iwghf_eb_window_av(m) /        &
                                          REAL( average_count_3d, kind=wp )
                       ENDDO
                    ENDIF
                    
                CASE ( 'usm_t_surf_wall' )
!
!--                 surface temperature for surfaces
                    IF ( l == -1 ) THEN
                       DO  m = 1, surf_usm_h%ns
                       surf_usm_h%t_surf_wall_av(m) =                               & 
                                          surf_usm_h%t_surf_wall_av(m) /            &
                                             REAL( average_count_3d, kind=wp )
                       ENDDO
                    ELSE
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%t_surf_wall_av(m) =                         &
                                          surf_usm_v(l)%t_surf_wall_av(m) /         &
                                          REAL( average_count_3d, kind=wp )
                       ENDDO
                    ENDIF
                    
                CASE ( 'usm_t_surf_window' )
!
!--                 surface temperature for window surfaces
                    IF ( l == -1 ) THEN
                       DO  m = 1, surf_usm_h%ns
                          surf_usm_h%t_surf_window_av(m) =                               &
                                             surf_usm_h%t_surf_window_av(m) /            &
                                             REAL( average_count_3d, kind=wp )
                       ENDDO
                    ELSE
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%t_surf_window_av(m) =                         &
                                          surf_usm_v(l)%t_surf_window_av(m) /         &
                                          REAL( average_count_3d, kind=wp )
                       ENDDO
                    ENDIF
                    
                CASE ( 'usm_t_surf_green' )
!
!--                 surface temperature for green surfaces
                    IF ( l == -1 ) THEN
                       DO  m = 1, surf_usm_h%ns
                          surf_usm_h%t_surf_green_av(m) =                               &
                                             surf_usm_h%t_surf_green_av(m) /            &
                                             REAL( average_count_3d, kind=wp )
                       ENDDO
                    ELSE
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%t_surf_green_av(m) =                         &
                                          surf_usm_v(l)%t_surf_green_av(m) /         &
                                          REAL( average_count_3d, kind=wp )
                       ENDDO
                    ENDIF
                    
                CASE ( 'usm_theta_10cm' )
!
!--                 near surface temperature for whole surfaces
                    IF ( l == -1 ) THEN
                       DO  m = 1, surf_usm_h%ns
                          surf_usm_h%pt_10cm_av(m) =                               &
                                             surf_usm_h%pt_10cm_av(m) /            &
                                             REAL( average_count_3d, kind=wp )
                       ENDDO
                    ELSE
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%pt_10cm_av(m) =                         &
                                          surf_usm_v(l)%pt_10cm_av(m) /         &
                                          REAL( average_count_3d, kind=wp )
                       ENDDO
                    ENDIF

                    
                CASE ( 'usm_t_wall' )
!
!--                 wall temperature for  iwl layer of walls and land
                    IF ( l == -1 ) THEN
                       DO  m = 1, surf_usm_h%ns
                          surf_usm_h%t_wall_av(iwl,m) =                           &
                                             surf_usm_h%t_wall_av(iwl,m) /        &
                                             REAL( average_count_3d, kind=wp )
                       ENDDO
                    ELSE
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%t_wall_av(iwl,m) =                     &
                                          surf_usm_v(l)%t_wall_av(iwl,m) /     &
                                          REAL( average_count_3d, kind=wp )
                       ENDDO
                    ENDIF

                CASE ( 'usm_t_window' )
!
!--                 window temperature for  iwl layer of walls and land
                    IF ( l == -1 ) THEN
                       DO  m = 1, surf_usm_h%ns
                          surf_usm_h%t_window_av(iwl,m) =                           &
                                             surf_usm_h%t_window_av(iwl,m) /        &
                                             REAL( average_count_3d, kind=wp )
                       ENDDO
                    ELSE
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%t_window_av(iwl,m) =                     &
                                          surf_usm_v(l)%t_window_av(iwl,m) /     &
                                          REAL( average_count_3d, kind=wp )
                       ENDDO
                    ENDIF

                CASE ( 'usm_t_green' )
!
!--                 green temperature for  iwl layer of walls and land
                    IF ( l == -1 ) THEN
                       DO  m = 1, surf_usm_h%ns
                          surf_usm_h%t_green_av(iwl,m) =                           &
                                             surf_usm_h%t_green_av(iwl,m) /        &
                                             REAL( average_count_3d, kind=wp )
                       ENDDO
                    ELSE
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%t_green_av(iwl,m) =                     &
                                          surf_usm_v(l)%t_green_av(iwl,m) /     &
                                          REAL( average_count_3d, kind=wp )
                       ENDDO
                    ENDIF
                    
                CASE ( 'usm_swc' )
!
!--                 soil water content for  iwl layer of walls and land
                    IF ( l == -1 ) THEN
                    DO  m = 1, surf_usm_h%ns
                       surf_usm_h%swc_av(iwl,m) =                           &
                                          surf_usm_h%swc_av(iwl,m) /        &
                                          REAL( average_count_3d, kind=wp )
                    ENDDO
                    ELSE
                       DO  m = 1, surf_usm_v(l)%ns
                          surf_usm_v(l)%swc_av(iwl,m) =                     &
                                          surf_usm_v(l)%swc_av(iwl,m) /     &
                                          REAL( average_count_3d, kind=wp )
                       ENDDO
                    ENDIF


           END SELECT

        ENDIF

        ENDIF

    END SUBROUTINE usm_3d_data_averaging



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Set internal Neumann boundary condition at outer soil grid points 
!> for temperature and humidity. 
!------------------------------------------------------------------------------!
 SUBROUTINE usm_boundary_condition
 
    IMPLICIT NONE

    INTEGER(iwp) :: i      !< grid index x-direction
    INTEGER(iwp) :: ioff   !< offset index x-direction indicating location of soil grid point
    INTEGER(iwp) :: j      !< grid index y-direction
    INTEGER(iwp) :: joff   !< offset index x-direction indicating location of soil grid point
    INTEGER(iwp) :: k      !< grid index z-direction
    INTEGER(iwp) :: koff   !< offset index x-direction indicating location of soil grid point
    INTEGER(iwp) :: l      !< running index surface-orientation
    INTEGER(iwp) :: m      !< running index surface elements

    koff = surf_usm_h%koff
    DO  m = 1, surf_usm_h%ns
       i = surf_usm_h%i(m)
       j = surf_usm_h%j(m)
       k = surf_usm_h%k(m)
       pt(k+koff,j,i) = pt(k,j,i)
    ENDDO

    DO  l = 0, 3
       ioff = surf_usm_v(l)%ioff
       joff = surf_usm_v(l)%joff
       DO  m = 1, surf_usm_v(l)%ns
          i = surf_usm_v(l)%i(m)
          j = surf_usm_v(l)%j(m)
          k = surf_usm_v(l)%k(m)
          pt(k,j+joff,i+ioff) = pt(k,j,i)
       ENDDO
    ENDDO

 END SUBROUTINE usm_boundary_condition


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine checks variables and assigns units.
!> It is called out from subroutine check_parameters.
!------------------------------------------------------------------------------!
    SUBROUTINE usm_check_data_output( variable, unit )

        IMPLICIT NONE

        CHARACTER(LEN=*),INTENT(IN)    ::  variable   !< 
        CHARACTER(LEN=*),INTENT(OUT)   ::  unit       !<

        INTEGER(iwp)                                  :: i,j,l         !< index
        CHARACTER(LEN=2)                              :: ls
        CHARACTER(LEN=varnamelength)                  :: var           !< TRIM(variable)
        INTEGER(iwp), PARAMETER                       :: nl1 = 15      !< number of directional usm variables
        CHARACTER(LEN=varnamelength), DIMENSION(nl1)  :: varlist1 = &  !< list of directional usm variables
                  (/'usm_wshf                      ', &
                    'usm_wghf                      ', &
                    'usm_wghf_window               ', &
                    'usm_wghf_green                ', &
                    'usm_iwghf                     ', &
                    'usm_iwghf_window              ', &
                    'usm_surfz                     ', &
                    'usm_surfwintrans              ', &
                    'usm_surfcat                   ', &
                    'usm_t_surf_wall               ', &
                    'usm_t_surf_window             ', &
                    'usm_t_surf_green              ', &
                    'usm_t_green                   ', &
                    'usm_qsws                      ', &
                    'usm_theta_10cm                '/)

        INTEGER(iwp), PARAMETER                       :: nl2 = 3       !< number of directional layer usm variables
        CHARACTER(LEN=varnamelength), DIMENSION(nl2)  :: varlist2 = &  !< list of directional layer usm variables
                  (/'usm_t_wall                    ', &
                    'usm_t_window                  ', &
                    'usm_t_green                   '/)

        INTEGER(iwp), PARAMETER                       :: nd = 5     !< number of directions
        CHARACTER(LEN=6), DIMENSION(nd), PARAMETER  :: dirname = &  !< direction names
                  (/'_roof ','_south','_north','_west ','_east '/)
        LOGICAL                                       :: lfound     !< flag if the variable is found


        lfound = .FALSE.

        var = TRIM(variable)

!
!--     check if variable exists
!--     directional variables
        DO i = 1, nl1
           DO j = 1, nd
              IF ( TRIM(var) == TRIM(varlist1(i))//TRIM(dirname(j)) ) THEN
                 lfound = .TRUE.
                 EXIT
              ENDIF
              IF ( lfound ) EXIT
           ENDDO
        ENDDO
        IF ( lfound ) GOTO 10
!
!--     directional layer variables
        DO i = 1, nl2
           DO j = 1, nd
              DO l = nzb_wall, nzt_wall
                 WRITE(ls,'(A1,I1)') '_',l
                 IF ( TRIM(var) == TRIM(varlist2(i))//TRIM(ls)//TRIM(dirname(j)) ) THEN
                    lfound = .TRUE.
                    EXIT
                 ENDIF
              ENDDO
              IF ( lfound ) EXIT
           ENDDO
        ENDDO
        IF ( .NOT.  lfound ) THEN
           unit = 'illegal'
           RETURN
        ENDIF
10      CONTINUE

        IF ( var(1:9)  == 'usm_wshf_'  .OR.  var(1:9) == 'usm_wghf_' .OR.                 &
             var(1:16) == 'usm_wghf_window_' .OR. var(1:15) == 'usm_wghf_green_' .OR.     &
             var(1:10) == 'usm_iwghf_' .OR. var(1:17) == 'usm_iwghf_window_'    .OR.      &
             var(1:17) == 'usm_surfwintrans_' .OR.                                        &
             var(1:9)  == 'usm_qsws_'  .OR.  var(1:13)  == 'usm_qsws_veg_'  .OR.          &
             var(1:13) == 'usm_qsws_liq_' ) THEN
            unit = 'W/m2'
        ELSE IF ( var(1:15) == 'usm_t_surf_wall'   .OR.  var(1:10) == 'usm_t_wall' .OR.   &
                  var(1:12) == 'usm_t_window' .OR. var(1:17) == 'usm_t_surf_window' .OR.  &
                  var(1:16) == 'usm_t_surf_green'  .OR.                                   &
                  var(1:11) == 'usm_t_green' .OR.  var(1:7) == 'usm_swc' .OR.             &
                  var(1:14) == 'usm_theta_10cm' )  THEN
            unit = 'K'
        ELSE IF ( var(1:9) == 'usm_surfz'  .OR.  var(1:11) == 'usm_surfcat' )  THEN
            unit = '1'
        ELSE
            unit = 'illegal'
        ENDIF

    END SUBROUTINE usm_check_data_output


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check parameters routine for urban surface model
!------------------------------------------------------------------------------!
    SUBROUTINE usm_check_parameters

       USE control_parameters,                                                 &
           ONLY:  bc_pt_b, bc_q_b, constant_flux_layer, large_scale_forcing,   &
                  lsf_surf, topography

       USE netcdf_data_input_mod,                                             &
            ONLY:  building_type_f

       IMPLICIT NONE

       INTEGER(iwp) ::  i        !< running index, x-dimension
       INTEGER(iwp) ::  j        !< running index, y-dimension

!
!--    Dirichlet boundary conditions are required as the surface fluxes are
!--    calculated from the temperature/humidity gradients in the urban surface
!--    model
       IF ( bc_pt_b == 'neumann'   .OR.   bc_q_b == 'neumann' )  THEN
          message_string = 'urban surface model requires setting of '//        &
                           'bc_pt_b = "dirichlet" and '//                      &
                           'bc_q_b  = "dirichlet"'
          CALL message( 'usm_check_parameters', 'PA0590', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( .NOT.  constant_flux_layer )  THEN
          message_string = 'urban surface model requires '//                   &
                           'constant_flux_layer = .T.'
          CALL message( 'usm_check_parameters', 'PA0084', 1, 2, 0, 6, 0 )
       ENDIF

       IF (  .NOT.  radiation )  THEN
          message_string = 'urban surface model requires '//                   &
                           'the radiation model to be switched on'
          CALL message( 'usm_check_parameters', 'PA0084', 1, 2, 0, 6, 0 )
       ENDIF
!        
!--    Surface forcing has to be disabled for LSF in case of enabled 
!--    urban surface module
       IF ( large_scale_forcing )  THEN
          lsf_surf = .FALSE.
       ENDIF
!
!--    Topography
       IF ( topography == 'flat' )  THEN
          message_string = 'topography /= "flat" is required '//               &
                           'when using the urban surface model'
          CALL message( 'usm_check_parameters', 'PA0592', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    naheatlayers
       IF ( naheatlayers > nzt )  THEN
          message_string = 'number of anthropogenic heat layers '//            &
                           '"naheatlayers" can not be larger than'//           &
                           ' number of domain layers "nzt"'
          CALL message( 'usm_check_parameters', 'PA0593', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Check if building types are set within a valid range.
       IF ( building_type < LBOUND( building_pars, 2 )  .AND.                  &
            building_type > UBOUND( building_pars, 2 ) )  THEN
          WRITE( message_string, * ) 'building_type = ', building_type,        &
                                     ' is out of the valid range'
          CALL message( 'usm_check_parameters', 'PA0529', 2, 2, 0, 6, 0 )
       ENDIF
       IF ( building_type_f%from_file )  THEN
          DO  i = nxl, nxr
             DO  j = nys, nyn
                IF ( building_type_f%var(j,i) /= building_type_f%fill  .AND.   &
              ( building_type_f%var(j,i) < LBOUND( building_pars, 2 )  .OR.    &
                building_type_f%var(j,i) > UBOUND( building_pars, 2 ) ) )      &
                THEN
                   WRITE( message_string, * ) 'building_type = is out of ' //  &
                                              'the valid range at (j,i) = ', j, i
                   CALL message( 'usm_check_parameters', 'PA0529', 2, 2, myid, 6, 0 )
                ENDIF
             ENDDO
          ENDDO
       ENDIF
    END SUBROUTINE usm_check_parameters


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Output of the 3D-arrays in netCDF and/or AVS format 
!> for variables of urban_surface model.
!> It resorts the urban surface module output quantities from surf style
!> indexing into temporary 3D array with indices (i,j,k).
!> It is called from subroutine data_output_3d.
!------------------------------------------------------------------------------!
    SUBROUTINE usm_data_output_3d( av, variable, found, local_pf, nzb_do, nzt_do )
        
        IMPLICIT NONE

        INTEGER(iwp), INTENT(IN)       ::  av        !< flag if averaged
        CHARACTER (len=*), INTENT(IN)  ::  variable  !< variable name
        INTEGER(iwp), INTENT(IN)       ::  nzb_do    !< lower limit of the data output (usually 0)
        INTEGER(iwp), INTENT(IN)       ::  nzt_do    !< vertical upper limit of the data output (usually nz_do3d)
        LOGICAL, INTENT(OUT)           ::  found     !< 
        REAL(sp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf   !< sp - it has to correspond to module data_output_3d
        REAL(sp), DIMENSION(nzb:nzt+1,nys:nyn,nxl:nxr)     ::  temp_pf    !< temp array for urban surface output procedure
        
        CHARACTER (len=varnamelength)                      :: var     !< trimmed variable name
        INTEGER(iwp), PARAMETER                            :: nd = 5  !< number of directions
        CHARACTER(len=6), DIMENSION(0:nd-1), PARAMETER     :: dirname = (/ '_roof ', '_south', '_north', '_west ', '_east ' /)
        INTEGER(iwp), DIMENSION(0:nd-1), PARAMETER         :: dirint =  (/    iup_u, isouth_u, inorth_u,  iwest_u,  ieast_u /)
        INTEGER(iwp), DIMENSION(0:nd-1), PARAMETER         :: diridx =  (/       -1,        1,        0,        3,        2 /)
                                                                      !< index for surf_*_v: 0:3 = (North, South, East, West)
        INTEGER(iwp)                   :: ids,idsint,idsidx
        INTEGER(iwp)                   :: i,j,k,iwl,istat, l, m  !< running indices

        found = .TRUE.
        temp_pf = -1._wp
        
        ids = -1
        var = TRIM(variable)
        DO i = 0, nd-1
            k = len(TRIM(var))
            j = len(TRIM(dirname(i)))
            IF ( TRIM(var(k-j+1:k)) == TRIM(dirname(i)) )  THEN
                ids = i
                idsint = dirint(ids)
                idsidx = diridx(ids)
                var = var(:k-j)
                EXIT
            ENDIF
        ENDDO
        IF ( ids == -1 )  THEN
            var = TRIM(variable)
        ENDIF
        IF ( var(1:11) == 'usm_t_wall_'  .AND.  len(TRIM(var)) >= 12 )  THEN
!
!--         wall layers
            READ(var(12:12), '(I1)', iostat=istat ) iwl
            IF ( istat == 0  .AND.  iwl >= nzb_wall  .AND.  iwl <= nzt_wall )  THEN
                var = var(1:10)
            ENDIF
        ENDIF
        IF ( var(1:13) == 'usm_t_window_'  .AND.  len(TRIM(var)) >= 14 )  THEN
!
!--         window layers
            READ(var(14:14), '(I1)', iostat=istat ) iwl
            IF ( istat == 0  .AND.  iwl >= nzb_wall  .AND.  iwl <= nzt_wall )  THEN
                var = var(1:12)
            ENDIF
        ENDIF
        IF ( var(1:12) == 'usm_t_green_'  .AND.  len(TRIM(var)) >= 13 )  THEN
!
!--         green layers
            READ(var(13:13), '(I1)', iostat=istat ) iwl
            IF ( istat == 0  .AND.  iwl >= nzb_wall  .AND.  iwl <= nzt_wall )  THEN
                var = var(1:11)
            ENDIF
        ENDIF
        IF ( var(1:8) == 'usm_swc_'  .AND.  len(TRIM(var)) >= 9 )  THEN
!
!--         green layers soil water content
            READ(var(9:9), '(I1)', iostat=istat ) iwl
            IF ( istat == 0  .AND.  iwl >= nzb_wall  .AND.  iwl <= nzt_wall )  THEN
                var = var(1:7)
            ENDIF
        ENDIF
        
        SELECT CASE ( TRIM(var) )

          CASE ( 'usm_surfz' )
!
!--           array of surface height (z)
              IF ( idsint == iup_u )  THEN
                 DO  m = 1, surf_usm_h%ns
                    i = surf_usm_h%i(m)
                    j = surf_usm_h%j(m)
                    k = surf_usm_h%k(m)
                    temp_pf(0,j,i) = MAX( temp_pf(0,j,i), REAL( k, KIND = sp) )
                 ENDDO
              ELSE
                 l = idsidx
                 DO  m = 1, surf_usm_v(l)%ns
                    i = surf_usm_v(l)%i(m)
                    j = surf_usm_v(l)%j(m)
                    k = surf_usm_v(l)%k(m)
                    temp_pf(0,j,i) = MAX( temp_pf(0,j,i), REAL( k, KIND = sp) + 1.0_sp )
                 ENDDO
              ENDIF

          CASE ( 'usm_surfcat' )
!
!--           surface category
              IF ( idsint == iup_u )  THEN
                 DO  m = 1, surf_usm_h%ns
                    i = surf_usm_h%i(m)
                    j = surf_usm_h%j(m)
                    k = surf_usm_h%k(m)
                    temp_pf(k,j,i) = surf_usm_h%surface_types(m)
                 ENDDO
              ELSE
                 l = idsidx
                 DO  m = 1, surf_usm_v(l)%ns
                    i = surf_usm_v(l)%i(m)
                    j = surf_usm_v(l)%j(m)
                    k = surf_usm_v(l)%k(m)
                    temp_pf(k,j,i) = surf_usm_v(l)%surface_types(m)
                 ENDDO
              ENDIF
              
          CASE ( 'usm_surfwintrans' )
!
!--           transmissivity window tiles
              IF ( idsint == iup_u )  THEN
                 DO  m = 1, surf_usm_h%ns
                    i = surf_usm_h%i(m)
                    j = surf_usm_h%j(m)
                    k = surf_usm_h%k(m)
                    temp_pf(k,j,i) = surf_usm_h%transmissivity(m)
                 ENDDO
              ELSE
                 l = idsidx
                 DO  m = 1, surf_usm_v(l)%ns
                    i = surf_usm_v(l)%i(m)
                    j = surf_usm_v(l)%j(m)
                    k = surf_usm_v(l)%k(m)
                    temp_pf(k,j,i) = surf_usm_v(l)%transmissivity(m)
                 ENDDO
              ENDIF

          CASE ( 'usm_wshf' )
!
!--           array of sensible heat flux from surfaces
              IF ( av == 0 )  THEN
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%wshf_eb(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%wshf_eb(m)
                    ENDDO
                 ENDIF
              ELSE
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%wshf_eb_av(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%wshf_eb_av(m)
                    ENDDO
                 ENDIF
              ENDIF
              
              
          CASE ( 'usm_qsws' )
!
!--           array of latent heat flux from surfaces
              IF ( av == 0 )  THEN
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%qsws(m) * l_v
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%qsws(m) * l_v
                    ENDDO
                 ENDIF
              ELSE
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%qsws_av(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%qsws_av(m)
                    ENDDO
                 ENDIF
              ENDIF
              
          CASE ( 'usm_qsws_veg' )
!
!--           array of latent heat flux from vegetation surfaces
              IF ( av == 0 )  THEN
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%qsws_veg(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%qsws_veg(m)
                    ENDDO
                 ENDIF
              ELSE
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%qsws_veg_av(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%qsws_veg_av(m)
                    ENDDO
                 ENDIF
              ENDIF
              
          CASE ( 'usm_qsws_liq' )
!
!--           array of latent heat flux from surfaces with liquid
              IF ( av == 0 )  THEN
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%qsws_liq(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%qsws_liq(m)
                    ENDDO
                 ENDIF
              ELSE
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%qsws_liq_av(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%qsws_liq_av(m)
                    ENDDO
                 ENDIF
              ENDIF

          CASE ( 'usm_wghf' )
!
!--           array of heat flux from ground (land, wall, roof)
              IF ( av == 0 )  THEN
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%wghf_eb(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%wghf_eb(m)
                    ENDDO
                 ENDIF
              ELSE
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%wghf_eb_av(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%wghf_eb_av(m)
                    ENDDO
                 ENDIF
              ENDIF

          CASE ( 'usm_wghf_window' )
!
!--           array of heat flux from window ground (land, wall, roof)
              IF ( av == 0 )  THEN
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%wghf_eb_window(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%wghf_eb_window(m)
                    ENDDO
                 ENDIF
              ELSE
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%wghf_eb_window_av(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%wghf_eb_window_av(m)
                    ENDDO
                 ENDIF
              ENDIF

          CASE ( 'usm_wghf_green' )
!
!--           array of heat flux from green ground (land, wall, roof)
              IF ( av == 0 )  THEN
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%wghf_eb_green(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%wghf_eb_green(m)
                    ENDDO
                 ENDIF
              ELSE
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%wghf_eb_green_av(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%wghf_eb_green_av(m)
                    ENDDO
                 ENDIF
              ENDIF

          CASE ( 'usm_iwghf' )
!
!--           array of heat flux from indoor ground (land, wall, roof)
              IF ( av == 0 )  THEN
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%iwghf_eb(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%iwghf_eb(m)
                    ENDDO
                 ENDIF
              ELSE
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%iwghf_eb_av(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%iwghf_eb_av(m)
                    ENDDO
                 ENDIF
              ENDIF

          CASE ( 'usm_iwghf_window' )
!
!--           array of heat flux from indoor window ground (land, wall, roof)
              IF ( av == 0 )  THEN
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%iwghf_eb_window(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%iwghf_eb_window(m)
                    ENDDO
                 ENDIF
              ELSE
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%iwghf_eb_window_av(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%iwghf_eb_window_av(m)
                    ENDDO
                 ENDIF
              ENDIF
              
          CASE ( 'usm_t_surf_wall' )
!
!--           surface temperature for surfaces
              IF ( av == 0 )  THEN
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = t_surf_wall_h(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = t_surf_wall_v(l)%t(m)
                    ENDDO
                 ENDIF
              ELSE
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%t_surf_wall_av(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%t_surf_wall_av(m)
                    ENDDO
                 ENDIF
              ENDIF
              
          CASE ( 'usm_t_surf_window' )
!
!--           surface temperature for window surfaces
              IF ( av == 0 )  THEN
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = t_surf_window_h(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = t_surf_window_v(l)%t(m)
                    ENDDO
                 ENDIF

              ELSE
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%t_surf_window_av(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%t_surf_window_av(m)
                    ENDDO

                 ENDIF

              ENDIF

          CASE ( 'usm_t_surf_green' )
!
!--           surface temperature for green surfaces
              IF ( av == 0 )  THEN
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = t_surf_green_h(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = t_surf_green_v(l)%t(m)
                    ENDDO
                 ENDIF

              ELSE
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%t_surf_green_av(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%t_surf_green_av(m)
                    ENDDO

                 ENDIF

              ENDIF

          CASE ( 'usm_theta_10cm' )
!
!--           near surface temperature for whole surfaces
              IF ( av == 0 )  THEN
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%pt_10cm(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%pt_10cm(m)
                    ENDDO
                 ENDIF
              
              
              ELSE
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%pt_10cm_av(m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%pt_10cm_av(m)
                    ENDDO

                  ENDIF
              ENDIF
             
          CASE ( 'usm_t_wall' )
!
!--           wall temperature for  iwl layer of walls and land
              IF ( av == 0 )  THEN
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = t_wall_h(iwl,m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = t_wall_v(l)%t(iwl,m)
                    ENDDO
                 ENDIF
              ELSE
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%t_wall_av(iwl,m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%t_wall_av(iwl,m)
                    ENDDO
                 ENDIF
              ENDIF
             
          CASE ( 'usm_t_window' )
!
!--           window temperature for iwl layer of walls and land
              IF ( av == 0 )  THEN
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = t_window_h(iwl,m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = t_window_v(l)%t(iwl,m)
                    ENDDO
                 ENDIF
              ELSE
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%t_window_av(iwl,m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%t_window_av(iwl,m)
                    ENDDO
                 ENDIF
              ENDIF

          CASE ( 'usm_t_green' )
!
!--           green temperature for  iwl layer of walls and land
              IF ( av == 0 )  THEN
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = t_green_h(iwl,m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = t_green_v(l)%t(iwl,m)
                    ENDDO
                 ENDIF
              ELSE
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%t_green_av(iwl,m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%t_green_av(iwl,m)
                    ENDDO
                 ENDIF
              ENDIF
              
              CASE ( 'usm_swc' )
!
!--           soil water content for  iwl layer of walls and land
              IF ( av == 0 )  THEN
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = swc_h(iwl,m)
                    ENDDO
                 ELSE

                 ENDIF
              ELSE
                 IF ( idsint == iup_u )  THEN
                    DO  m = 1, surf_usm_h%ns
                       i = surf_usm_h%i(m)
                       j = surf_usm_h%j(m)
                       k = surf_usm_h%k(m)
                       temp_pf(k,j,i) = surf_usm_h%swc_av(iwl,m)
                    ENDDO
                 ELSE
                    l = idsidx
                    DO  m = 1, surf_usm_v(l)%ns
                       i = surf_usm_v(l)%i(m)
                       j = surf_usm_v(l)%j(m)
                       k = surf_usm_v(l)%k(m)
                       temp_pf(k,j,i) = surf_usm_v(l)%swc_av(iwl,m)
                    ENDDO
                 ENDIF
              ENDIF

             
          CASE DEFAULT
              found = .FALSE.
              RETURN
        END SELECT

!
!--     Rearrange dimensions for NetCDF output
!--     FIXME: this may generate FPE overflow upon conversion from DP to SP
        DO  j = nys, nyn
            DO  i = nxl, nxr
                DO  k = nzb_do, nzt_do
                    local_pf(i,j,k) = temp_pf(k,j,i)
                ENDDO
            ENDDO
        ENDDO
        
    END SUBROUTINE usm_data_output_3d
    

!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Soubroutine defines appropriate grid for netcdf variables.
!> It is called out from subroutine netcdf.
!------------------------------------------------------------------------------!
    SUBROUTINE usm_define_netcdf_grid( variable, found, grid_x, grid_y, grid_z )
    
        IMPLICIT NONE

        CHARACTER (len=*), INTENT(IN)  ::  variable    !< 
        LOGICAL, INTENT(OUT)           ::  found       !< 
        CHARACTER (len=*), INTENT(OUT) ::  grid_x      !< 
        CHARACTER (len=*), INTENT(OUT) ::  grid_y      !< 
        CHARACTER (len=*), INTENT(OUT) ::  grid_z      !< 

        CHARACTER (len=varnamelength)  :: var

        var = TRIM(variable)
        IF ( var(1:9) == 'usm_wshf_'  .OR.  var(1:9) == 'usm_wghf_'  .OR.                   &
             var(1:16) == 'usm_wghf_window_'  .OR. var(1:15) == 'usm_wghf_green_' .OR.      &
             var(1:10) == 'usm_iwghf_'  .OR. var(1:17) == 'usm_iwghf_window_' .OR.          &
             var(1:9) == 'usm_qsws_'  .OR.  var(1:13) == 'usm_qsws_veg_'  .OR.              &
             var(1:13) == 'usm_qsws_liq_' .OR.                                              &
             var(1:15) == 'usm_t_surf_wall'  .OR.  var(1:10) == 'usm_t_wall'  .OR.          &
             var(1:17) == 'usm_t_surf_window'  .OR.  var(1:12) == 'usm_t_window'  .OR.      &
             var(1:16) == 'usm_t_surf_green'  .OR. var(1:11) == 'usm_t_green' .OR.          &
             var(1:15) == 'usm_theta_10cm' .OR.                                             &
             var(1:9) == 'usm_surfz'  .OR.  var(1:11) == 'usm_surfcat'  .OR.                &
             var(1:16) == 'usm_surfwintrans'  .OR. var(1:7) == 'usm_swc' ) THEN

            found = .TRUE.
            grid_x = 'x'
            grid_y = 'y'
            grid_z = 'zu'
        ELSE
            found  = .FALSE.
            grid_x = 'none'
            grid_y = 'none'
            grid_z = 'none'
        ENDIF

    END SUBROUTINE usm_define_netcdf_grid
    

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of the wall surface model
!------------------------------------------------------------------------------!
    SUBROUTINE usm_init_material_model

        IMPLICIT NONE

        INTEGER(iwp) ::  k, l, m            !< running indices
        
        IF ( debug_output )  CALL debug_message( 'usm_init_material_model', 'start' )

!
!--     Calculate wall grid spacings. 
!--     Temperature is defined at the center of the wall layers,
!--     whereas gradients/fluxes are defined at the edges (_stag)      
!--     apply for all particular surface grids. First for horizontal surfaces
        DO  m = 1, surf_usm_h%ns

           surf_usm_h%dz_wall(nzb_wall,m) = surf_usm_h%zw(nzb_wall,m)
           DO k = nzb_wall+1, nzt_wall
               surf_usm_h%dz_wall(k,m) = surf_usm_h%zw(k,m) -                  &
                                         surf_usm_h%zw(k-1,m)
           ENDDO
           surf_usm_h%dz_window(nzb_wall,m) = surf_usm_h%zw_window(nzb_wall,m)
           DO k = nzb_wall+1, nzt_wall
               surf_usm_h%dz_window(k,m) = surf_usm_h%zw_window(k,m) -         &
                                         surf_usm_h%zw_window(k-1,m)
           ENDDO
           
           surf_usm_h%dz_wall(nzt_wall+1,m) = surf_usm_h%dz_wall(nzt_wall,m)

           DO k = nzb_wall, nzt_wall-1
               surf_usm_h%dz_wall_stag(k,m) = 0.5 * (                          &
                           surf_usm_h%dz_wall(k+1,m) + surf_usm_h%dz_wall(k,m) )
           ENDDO
           surf_usm_h%dz_wall_stag(nzt_wall,m) = surf_usm_h%dz_wall(nzt_wall,m)
           
           surf_usm_h%dz_window(nzt_wall+1,m) = surf_usm_h%dz_window(nzt_wall,m)

           DO k = nzb_wall, nzt_wall-1
               surf_usm_h%dz_window_stag(k,m) = 0.5 * (                        &
                           surf_usm_h%dz_window(k+1,m) + surf_usm_h%dz_window(k,m) )
           ENDDO
           surf_usm_h%dz_window_stag(nzt_wall,m) = surf_usm_h%dz_window(nzt_wall,m)

           IF (surf_usm_h%green_type_roof(m) == 2.0_wp ) THEN
!
!-- extensive green roof
!-- set ratio of substrate layer thickness, soil-type and LAI
              soil_type = 3
              surf_usm_h%lai(m) = 2.0_wp
             
              surf_usm_h%zw_green(nzb_wall,m)   = 0.05_wp
              surf_usm_h%zw_green(nzb_wall+1,m) = 0.10_wp
              surf_usm_h%zw_green(nzb_wall+2,m) = 0.15_wp
              surf_usm_h%zw_green(nzb_wall+3,m) = 0.20_wp
           ELSE
!
!-- intensiv green roof
!-- set ratio of substrate layer thickness, soil-type and LAI
              soil_type = 6
              surf_usm_h%lai(m) = 4.0_wp
             
              surf_usm_h%zw_green(nzb_wall,m)   = 0.05_wp
              surf_usm_h%zw_green(nzb_wall+1,m) = 0.10_wp
              surf_usm_h%zw_green(nzb_wall+2,m) = 0.40_wp
              surf_usm_h%zw_green(nzb_wall+3,m) = 0.80_wp
           ENDIF
           
           surf_usm_h%dz_green(nzb_wall,m) = surf_usm_h%zw_green(nzb_wall,m)
           DO k = nzb_wall+1, nzt_wall
               surf_usm_h%dz_green(k,m) = surf_usm_h%zw_green(k,m) -           &
                                         surf_usm_h%zw_green(k-1,m)
           ENDDO
           surf_usm_h%dz_green(nzt_wall+1,m) = surf_usm_h%dz_green(nzt_wall,m)

           DO k = nzb_wall, nzt_wall-1
               surf_usm_h%dz_green_stag(k,m) = 0.5 * (                         &
                           surf_usm_h%dz_green(k+1,m) + surf_usm_h%dz_green(k,m) )
           ENDDO
           surf_usm_h%dz_green_stag(nzt_wall,m) = surf_usm_h%dz_green(nzt_wall,m)
           
          IF ( alpha_vangenuchten == 9999999.9_wp )  THEN
             alpha_vangenuchten = soil_pars(0,soil_type)
          ENDIF

          IF ( l_vangenuchten == 9999999.9_wp )  THEN
             l_vangenuchten = soil_pars(1,soil_type)
          ENDIF

          IF ( n_vangenuchten == 9999999.9_wp )  THEN
             n_vangenuchten = soil_pars(2,soil_type)            
          ENDIF

          IF ( hydraulic_conductivity == 9999999.9_wp )  THEN
             hydraulic_conductivity = soil_pars(3,soil_type)            
          ENDIF

          IF ( saturation_moisture == 9999999.9_wp )  THEN
             saturation_moisture = m_soil_pars(0,soil_type)           
          ENDIF

          IF ( field_capacity == 9999999.9_wp )  THEN
             field_capacity = m_soil_pars(1,soil_type)           
          ENDIF

          IF ( wilting_point == 9999999.9_wp )  THEN
             wilting_point = m_soil_pars(2,soil_type)            
          ENDIF

          IF ( residual_moisture == 9999999.9_wp )  THEN
             residual_moisture = m_soil_pars(3,soil_type)       
          ENDIF
          
          DO k = nzb_wall, nzt_wall+1
             swc_h(k,m) = field_capacity
             rootfr_h(k,m) = 0.5_wp
             surf_usm_h%alpha_vg_green(m)      = alpha_vangenuchten
             surf_usm_h%l_vg_green(m)          = l_vangenuchten
             surf_usm_h%n_vg_green(m)          = n_vangenuchten 
             surf_usm_h%gamma_w_green_sat(k,m) = hydraulic_conductivity
             swc_sat_h(k,m)                    = saturation_moisture
             fc_h(k,m)                         = field_capacity
             wilt_h(k,m)                       = wilting_point
             swc_res_h(k,m)                    = residual_moisture
          ENDDO

        ENDDO

        surf_usm_h%ddz_wall        = 1.0_wp / surf_usm_h%dz_wall
        surf_usm_h%ddz_wall_stag   = 1.0_wp / surf_usm_h%dz_wall_stag
        surf_usm_h%ddz_window      = 1.0_wp / surf_usm_h%dz_window
        surf_usm_h%ddz_window_stag = 1.0_wp / surf_usm_h%dz_window_stag
        surf_usm_h%ddz_green       = 1.0_wp / surf_usm_h%dz_green
        surf_usm_h%ddz_green_stag  = 1.0_wp / surf_usm_h%dz_green_stag
!        
!--     For vertical surfaces
        DO  l = 0, 3
           DO  m = 1, surf_usm_v(l)%ns
              surf_usm_v(l)%dz_wall(nzb_wall,m) = surf_usm_v(l)%zw(nzb_wall,m)
              DO k = nzb_wall+1, nzt_wall
                  surf_usm_v(l)%dz_wall(k,m) = surf_usm_v(l)%zw(k,m) -         &
                                               surf_usm_v(l)%zw(k-1,m)
              ENDDO
              surf_usm_v(l)%dz_window(nzb_wall,m) = surf_usm_v(l)%zw_window(nzb_wall,m)
              DO k = nzb_wall+1, nzt_wall
                  surf_usm_v(l)%dz_window(k,m) = surf_usm_v(l)%zw_window(k,m) - &
                                               surf_usm_v(l)%zw_window(k-1,m)
              ENDDO
              surf_usm_v(l)%dz_green(nzb_wall,m) = surf_usm_v(l)%zw_green(nzb_wall,m)
              DO k = nzb_wall+1, nzt_wall
                  surf_usm_v(l)%dz_green(k,m) = surf_usm_v(l)%zw_green(k,m) - &
                                               surf_usm_v(l)%zw_green(k-1,m)
              ENDDO
           
              surf_usm_v(l)%dz_wall(nzt_wall+1,m) =                            &
                                              surf_usm_v(l)%dz_wall(nzt_wall,m)

              DO k = nzb_wall, nzt_wall-1
                  surf_usm_v(l)%dz_wall_stag(k,m) = 0.5 * (                    &
                                                surf_usm_v(l)%dz_wall(k+1,m) + &
                                                surf_usm_v(l)%dz_wall(k,m) )
              ENDDO
              surf_usm_v(l)%dz_wall_stag(nzt_wall,m) =                         &
                                              surf_usm_v(l)%dz_wall(nzt_wall,m)
              surf_usm_v(l)%dz_window(nzt_wall+1,m) =                          &
                                              surf_usm_v(l)%dz_window(nzt_wall,m)

              DO k = nzb_wall, nzt_wall-1
                  surf_usm_v(l)%dz_window_stag(k,m) = 0.5 * (                    &
                                                surf_usm_v(l)%dz_window(k+1,m) + &
                                                surf_usm_v(l)%dz_window(k,m) )
              ENDDO
              surf_usm_v(l)%dz_window_stag(nzt_wall,m) =                         &
                                              surf_usm_v(l)%dz_window(nzt_wall,m)
              surf_usm_v(l)%dz_green(nzt_wall+1,m) =                             &
                                              surf_usm_v(l)%dz_green(nzt_wall,m)

              DO k = nzb_wall, nzt_wall-1
                  surf_usm_v(l)%dz_green_stag(k,m) = 0.5 * (                    &
                                                surf_usm_v(l)%dz_green(k+1,m) + &
                                                surf_usm_v(l)%dz_green(k,m) )
              ENDDO
              surf_usm_v(l)%dz_green_stag(nzt_wall,m) =                         &
                                              surf_usm_v(l)%dz_green(nzt_wall,m)
           ENDDO
           surf_usm_v(l)%ddz_wall        = 1.0_wp / surf_usm_v(l)%dz_wall
           surf_usm_v(l)%ddz_wall_stag   = 1.0_wp / surf_usm_v(l)%dz_wall_stag
           surf_usm_v(l)%ddz_window      = 1.0_wp / surf_usm_v(l)%dz_window
           surf_usm_v(l)%ddz_window_stag = 1.0_wp / surf_usm_v(l)%dz_window_stag
           surf_usm_v(l)%ddz_green       = 1.0_wp / surf_usm_v(l)%dz_green
           surf_usm_v(l)%ddz_green_stag  = 1.0_wp / surf_usm_v(l)%dz_green_stag
        ENDDO      

        
        IF ( debug_output )  CALL debug_message( 'usm_init_material_model', 'end' )

    END SUBROUTINE usm_init_material_model

 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of the urban surface model
!------------------------------------------------------------------------------!
    SUBROUTINE usm_init

        USE arrays_3d,                                                         &
            ONLY:  zw

        USE netcdf_data_input_mod,                                             &
            ONLY:  building_pars_f, building_surface_pars_f, building_type_f,  &
                   terrain_height_f
    
        IMPLICIT NONE

        INTEGER(iwp) ::  i                   !< loop index x-dirction
        INTEGER(iwp) ::  ind_alb_green       !< index in input list for green albedo
        INTEGER(iwp) ::  ind_alb_wall        !< index in input list for wall albedo
        INTEGER(iwp) ::  ind_alb_win         !< index in input list for window albedo
        INTEGER(iwp) ::  ind_emis_wall       !< index in input list for wall emissivity
        INTEGER(iwp) ::  ind_emis_green      !< index in input list for green emissivity
        INTEGER(iwp) ::  ind_emis_win        !< index in input list for window emissivity
        INTEGER(iwp) ::  ind_green_frac_w    !< index in input list for green fraction on wall
        INTEGER(iwp) ::  ind_green_frac_r    !< index in input list for green fraction on roof
        INTEGER(iwp) ::  ind_hc1             !< index in input list for heat capacity at first wall layer
        INTEGER(iwp) ::  ind_hc1_win         !< index in input list for heat capacity at first window layer
        INTEGER(iwp) ::  ind_hc2             !< index in input list for heat capacity at second wall layer
        INTEGER(iwp) ::  ind_hc2_win         !< index in input list for heat capacity at second window layer
        INTEGER(iwp) ::  ind_hc3             !< index in input list for heat capacity at third wall layer
        INTEGER(iwp) ::  ind_hc3_win         !< index in input list for heat capacity at third window layer
        INTEGER(iwp) ::  ind_lai_r           !< index in input list for LAI on roof
        INTEGER(iwp) ::  ind_lai_w           !< index in input list for LAI on wall
        INTEGER(iwp) ::  ind_tc1             !< index in input list for thermal conductivity at first wall layer
        INTEGER(iwp) ::  ind_tc1_win         !< index in input list for thermal conductivity at first window layer
        INTEGER(iwp) ::  ind_tc2             !< index in input list for thermal conductivity at second wall layer
        INTEGER(iwp) ::  ind_tc2_win         !< index in input list for thermal conductivity at second window layer
        INTEGER(iwp) ::  ind_tc3             !< index in input list for thermal conductivity at third wall layer
        INTEGER(iwp) ::  ind_tc3_win         !< index in input list for thermal conductivity at third window layer
        INTEGER(iwp) ::  ind_thick_1         !< index in input list for thickness of first wall layer 
        INTEGER(iwp) ::  ind_thick_1_win     !< index in input list for thickness of first window layer 
        INTEGER(iwp) ::  ind_thick_2         !< index in input list for thickness of second wall layer 
        INTEGER(iwp) ::  ind_thick_2_win     !< index in input list for thickness of second window layer
        INTEGER(iwp) ::  ind_thick_3         !< index in input list for thickness of third wall layer
        INTEGER(iwp) ::  ind_thick_3_win     !< index in input list for thickness of third window layer
        INTEGER(iwp) ::  ind_thick_4         !< index in input list for thickness of fourth wall layer
        INTEGER(iwp) ::  ind_thick_4_win     !< index in input list for thickness of fourth window layer
        INTEGER(iwp) ::  ind_trans           !< index in input list for window transmissivity
        INTEGER(iwp) ::  ind_wall_frac       !< index in input list for wall fraction
        INTEGER(iwp) ::  ind_win_frac        !< index in input list for window fraction
        INTEGER(iwp) ::  ind_z0              !< index in input list for z0
        INTEGER(iwp) ::  ind_z0qh            !< index in input list for z0h / z0q
        INTEGER(iwp) ::  is                  !< loop index input surface element
        INTEGER(iwp) ::  j                   !< loop index y-dirction
        INTEGER(iwp) ::  k                   !< loop index z-dirction
        INTEGER(iwp) ::  l                   !< loop index surface orientation
        INTEGER(iwp) ::  m                   !< loop index surface element
        INTEGER(iwp) ::  st                  !< dummy

        LOGICAL      ::  relative_fractions_corrected !< flag indicating if relative surface fractions require normalization

        REAL(wp)     ::  c, tin, twin
        REAL(wp)     ::  ground_floor_level_l         !< local height of ground floor level
        REAL(wp)     ::  sum_frac                     !< sum of the relative material fractions at a surface element
        REAL(wp)     ::  z_agl                        !< height of the surface element above terrain

        IF ( debug_output )  CALL debug_message( 'usm_init', 'start' )

        CALL cpu_log( log_point_s(78), 'usm_init', 'start' )
!
!--     surface forcing have to be disabled for LSF 
!--     in case of enabled urban surface module
        IF ( large_scale_forcing )  THEN
            lsf_surf = .FALSE.
        ENDIF
!
!--     Flag surface elements belonging to the ground floor level. Therefore, 
!--     use terrain height array from file, if available. This flag is later used
!--     to control initialization of surface attributes.
!--     Todo: for the moment disable initialization of building roofs with
!--     ground-floor-level properties. 
        surf_usm_h%ground_level = .FALSE. 

        DO  l = 0, 3
           surf_usm_v(l)%ground_level = .FALSE.
           DO  m = 1, surf_usm_v(l)%ns
              i = surf_usm_v(l)%i(m) + surf_usm_v(l)%ioff
              j = surf_usm_v(l)%j(m) + surf_usm_v(l)%joff
              k = surf_usm_v(l)%k(m)
!
!--           Determine local ground level. Level 1 - default value,
!--           level 2 - initialization according to building type, 
!--           level 3 - initialization from value read from file.
              ground_floor_level_l = ground_floor_level
              
              IF ( building_type_f%from_file )  THEN
                  ground_floor_level_l =                                       &
                              building_pars(ind_gflh,building_type_f%var(j,i))
              ENDIF
              
              IF ( building_pars_f%from_file )  THEN
                 IF ( building_pars_f%pars_xy(ind_gflh,j,i) /=                 &
                      building_pars_f%fill )                                   &
                    ground_floor_level_l = building_pars_f%pars_xy(ind_gflh,j,i)
              ENDIF
!
!--           Determine height of surface element above ground level. Please 
!--           note, height of surface element is determined with respect to
!--           its height above ground of the reference grid point in atmosphere,
!--           Therefore, substract the offset values when assessing the terrain
!--           height.
              IF ( terrain_height_f%from_file )  THEN
                 z_agl = zw(k) - terrain_height_f%var(j-surf_usm_v(l)%joff,    &
                                                      i-surf_usm_v(l)%ioff)
              ELSE
                 z_agl = zw(k)
              ENDIF
!
!--           Set flag for ground level
              IF ( z_agl <= ground_floor_level_l )                             &
                 surf_usm_v(l)%ground_level(m) = .TRUE.

           ENDDO
        ENDDO
!
!--     Initialization of resistances. 
        DO  m = 1, surf_usm_h%ns
           surf_usm_h%r_a(m)        = 50.0_wp
           surf_usm_h%r_a_green(m)  = 50.0_wp
           surf_usm_h%r_a_window(m) = 50.0_wp
        ENDDO
        DO  l = 0, 3
           DO  m = 1, surf_usm_v(l)%ns
              surf_usm_v(l)%r_a(m)        = 50.0_wp
              surf_usm_v(l)%r_a_green(m)  = 50.0_wp
              surf_usm_v(l)%r_a_window(m) = 50.0_wp
           ENDDO
        ENDDO
        
!
!--    Map values onto horizontal elemements
       DO  m = 1, surf_usm_h%ns
          surf_usm_h%r_canopy(m)     = 200.0_wp !< canopy_resistance
          surf_usm_h%r_canopy_min(m) = 200.0_wp !< min_canopy_resistance
          surf_usm_h%g_d(m)          = 0.0_wp   !< canopy_resistance_coefficient
       ENDDO
!
!--    Map values onto vertical elements, even though this does not make
!--    much sense.
       DO  l = 0, 3
          DO  m = 1, surf_usm_v(l)%ns
             surf_usm_v(l)%r_canopy(m)     = 200.0_wp !< canopy_resistance
             surf_usm_v(l)%r_canopy_min(m) = 200.0_wp !< min_canopy_resistance
             surf_usm_v(l)%g_d(m)          = 0.0_wp   !< canopy_resistance_coefficient
          ENDDO
       ENDDO

!
!--     Initialize urban-type surface attribute. According to initialization in 
!--     land-surface model, follow a 3-level approach. 
!--     Level 1 - initialization via default attributes 
        DO  m = 1, surf_usm_h%ns
!
!--        Now, all horizontal surfaces are roof surfaces (?)
           surf_usm_h%isroof_surf(m)   = .TRUE.
           surf_usm_h%surface_types(m) = roof_category         !< default category for root surface
!
!--        In order to distinguish between ground floor level and 
!--        above-ground-floor level surfaces, set input indices.

           ind_green_frac_r = MERGE( ind_green_frac_r_gfl, ind_green_frac_r_agfl, &
                                     surf_usm_h%ground_level(m) )
           ind_lai_r        = MERGE( ind_lai_r_gfl,        ind_lai_r_agfl,        &
                                     surf_usm_h%ground_level(m) )
           ind_z0           = MERGE( ind_z0_gfl,           ind_z0_agfl,           &
                                     surf_usm_h%ground_level(m) )
           ind_z0qh         = MERGE( ind_z0qh_gfl,         ind_z0qh_agfl,         &
                                     surf_usm_h%ground_level(m) )
!
!--        Store building type and its name on each surface element
           surf_usm_h%building_type(m)      = building_type
           surf_usm_h%building_type_name(m) = building_type_name(building_type)
!
!--        Initialize relatvie wall- (0), green- (1) and window (2) fractions
           surf_usm_h%frac(m,ind_veg_wall)  = building_pars(ind_wall_frac_r,building_type)   
           surf_usm_h%frac(m,ind_pav_green) = building_pars(ind_green_frac_r,building_type)  
           surf_usm_h%frac(m,ind_wat_win)   = building_pars(ind_win_frac_r,building_type)  
           surf_usm_h%lai(m)                = building_pars(ind_lai_r,building_type)  

           surf_usm_h%rho_c_wall(nzb_wall,m)   = building_pars(ind_hc1_wall_r,building_type)  
           surf_usm_h%rho_c_wall(nzb_wall+1,m) = building_pars(ind_hc1_wall_r,building_type)
           surf_usm_h%rho_c_wall(nzb_wall+2,m) = building_pars(ind_hc2_wall_r,building_type)
           surf_usm_h%rho_c_wall(nzb_wall+3,m) = building_pars(ind_hc3_wall_r,building_type)    
           surf_usm_h%lambda_h(nzb_wall,m)   = building_pars(ind_tc1_wall_r,building_type)  
           surf_usm_h%lambda_h(nzb_wall+1,m) = building_pars(ind_tc1_wall_r,building_type) 
           surf_usm_h%lambda_h(nzb_wall+2,m) = building_pars(ind_tc2_wall_r,building_type)
           surf_usm_h%lambda_h(nzb_wall+3,m) = building_pars(ind_tc3_wall_r,building_type)    
           surf_usm_h%rho_c_green(nzb_wall,m)   = rho_c_soil !building_pars(ind_hc1_wall_r,building_type)  
           surf_usm_h%rho_c_green(nzb_wall+1,m) = rho_c_soil !building_pars(ind_hc1_wall_r,building_type)
           surf_usm_h%rho_c_green(nzb_wall+2,m) = rho_c_soil !building_pars(ind_hc2_wall_r,building_type)
           surf_usm_h%rho_c_green(nzb_wall+3,m) = rho_c_soil !building_pars(ind_hc3_wall_r,building_type)    
           surf_usm_h%lambda_h_green(nzb_wall,m)   = lambda_h_green_sm !building_pars(ind_tc1_wall_r,building_type)  
           surf_usm_h%lambda_h_green(nzb_wall+1,m) = lambda_h_green_sm !building_pars(ind_tc1_wall_r,building_type) 
           surf_usm_h%lambda_h_green(nzb_wall+2,m) = lambda_h_green_sm !building_pars(ind_tc2_wall_r,building_type)
           surf_usm_h%lambda_h_green(nzb_wall+3,m) = lambda_h_green_sm !building_pars(ind_tc3_wall_r,building_type)
           surf_usm_h%rho_c_window(nzb_wall,m)   = building_pars(ind_hc1_win_r,building_type)  
           surf_usm_h%rho_c_window(nzb_wall+1,m) = building_pars(ind_hc1_win_r,building_type)
           surf_usm_h%rho_c_window(nzb_wall+2,m) = building_pars(ind_hc2_win_r,building_type)
           surf_usm_h%rho_c_window(nzb_wall+3,m) = building_pars(ind_hc3_win_r,building_type)    
           surf_usm_h%lambda_h_window(nzb_wall,m)   = building_pars(ind_tc1_win_r,building_type)  
           surf_usm_h%lambda_h_window(nzb_wall+1,m) = building_pars(ind_tc1_win_r,building_type) 
           surf_usm_h%lambda_h_window(nzb_wall+2,m) = building_pars(ind_tc2_win_r,building_type)
           surf_usm_h%lambda_h_window(nzb_wall+3,m) = building_pars(ind_tc3_win_r,building_type)    

           surf_usm_h%target_temp_summer(m)  = building_pars(ind_indoor_target_temp_summer,building_type)    
           surf_usm_h%target_temp_winter(m)  = building_pars(ind_indoor_target_temp_winter,building_type)    
!
!--        emissivity of wall-, green- and window fraction 
           surf_usm_h%emissivity(m,ind_veg_wall)  = building_pars(ind_emis_wall_r,building_type)
           surf_usm_h%emissivity(m,ind_pav_green) = building_pars(ind_emis_green_r,building_type)
           surf_usm_h%emissivity(m,ind_wat_win)   = building_pars(ind_emis_win_r,building_type)

           surf_usm_h%transmissivity(m)      = building_pars(ind_trans_r,building_type)

           surf_usm_h%z0(m)                  = building_pars(ind_z0,building_type)
           surf_usm_h%z0h(m)                 = building_pars(ind_z0qh,building_type)
           surf_usm_h%z0q(m)                 = building_pars(ind_z0qh,building_type)
!
!--        albedo type for wall fraction, green fraction, window fraction
           surf_usm_h%albedo_type(m,ind_veg_wall)  = INT( building_pars(ind_alb_wall_r,building_type)  )
           surf_usm_h%albedo_type(m,ind_pav_green) = INT( building_pars(ind_alb_green_r,building_type) )
           surf_usm_h%albedo_type(m,ind_wat_win)   = INT( building_pars(ind_alb_win_r,building_type)   )

           surf_usm_h%zw(nzb_wall,m)         = building_pars(ind_thick_1_wall_r,building_type)
           surf_usm_h%zw(nzb_wall+1,m)       = building_pars(ind_thick_2_wall_r,building_type)
           surf_usm_h%zw(nzb_wall+2,m)       = building_pars(ind_thick_3_wall_r,building_type)
           surf_usm_h%zw(nzb_wall+3,m)       = building_pars(ind_thick_4_wall_r,building_type)
           
           surf_usm_h%zw_green(nzb_wall,m)         = building_pars(ind_thick_1_wall_r,building_type)
           surf_usm_h%zw_green(nzb_wall+1,m)       = building_pars(ind_thick_2_wall_r,building_type)
           surf_usm_h%zw_green(nzb_wall+2,m)       = building_pars(ind_thick_3_wall_r,building_type)
           surf_usm_h%zw_green(nzb_wall+3,m)       = building_pars(ind_thick_4_wall_r,building_type)
           
           surf_usm_h%zw_window(nzb_wall,m)         = building_pars(ind_thick_1_win_r,building_type)
           surf_usm_h%zw_window(nzb_wall+1,m)       = building_pars(ind_thick_2_win_r,building_type)
           surf_usm_h%zw_window(nzb_wall+2,m)       = building_pars(ind_thick_3_win_r,building_type)
           surf_usm_h%zw_window(nzb_wall+3,m)       = building_pars(ind_thick_4_win_r,building_type)

           surf_usm_h%c_surface(m)           = building_pars(ind_c_surface,building_type)  
           surf_usm_h%lambda_surf(m)         = building_pars(ind_lambda_surf,building_type)  
           surf_usm_h%c_surface_green(m)     = building_pars(ind_c_surface_green,building_type)  
           surf_usm_h%lambda_surf_green(m)   = building_pars(ind_lambda_surf_green,building_type)  
           surf_usm_h%c_surface_window(m)    = building_pars(ind_c_surface_win,building_type)  
           surf_usm_h%lambda_surf_window(m)  = building_pars(ind_lambda_surf_win,building_type)  
           
           surf_usm_h%green_type_roof(m)     = building_pars(ind_green_type_roof,building_type)

        ENDDO

        DO  l = 0, 3
           DO  m = 1, surf_usm_v(l)%ns

              surf_usm_v(l)%surface_types(m) = wall_category         !< default category for root surface
!
!--           In order to distinguish between ground floor level and 
!--           above-ground-floor level surfaces, set input indices.
              ind_alb_green    = MERGE( ind_alb_green_gfl,    ind_alb_green_agfl,    &
                                        surf_usm_v(l)%ground_level(m) )
              ind_alb_wall     = MERGE( ind_alb_wall_gfl,     ind_alb_wall_agfl,     &
                                        surf_usm_v(l)%ground_level(m) )
              ind_alb_win      = MERGE( ind_alb_win_gfl,      ind_alb_win_agfl,      &
                                        surf_usm_v(l)%ground_level(m) )
              ind_wall_frac    = MERGE( ind_wall_frac_gfl,    ind_wall_frac_agfl,    &
                                        surf_usm_v(l)%ground_level(m) )
              ind_win_frac     = MERGE( ind_win_frac_gfl,     ind_win_frac_agfl,     &
                                        surf_usm_v(l)%ground_level(m) )
              ind_green_frac_w = MERGE( ind_green_frac_w_gfl, ind_green_frac_w_agfl, &
                                        surf_usm_v(l)%ground_level(m) )
              ind_green_frac_r = MERGE( ind_green_frac_r_gfl, ind_green_frac_r_agfl, &
                                        surf_usm_v(l)%ground_level(m) )
              ind_lai_r        = MERGE( ind_lai_r_gfl,        ind_lai_r_agfl,        &
                                        surf_usm_v(l)%ground_level(m) )
              ind_lai_w        = MERGE( ind_lai_w_gfl,        ind_lai_w_agfl,        &
                                        surf_usm_v(l)%ground_level(m) )
              ind_hc1          = MERGE( ind_hc1_gfl,          ind_hc1_agfl,          &
                                        surf_usm_v(l)%ground_level(m) )
              ind_hc1_win      = MERGE( ind_hc1_win_gfl,      ind_hc1_win_agfl,      &
                                        surf_usm_v(l)%ground_level(m) )
              ind_hc2          = MERGE( ind_hc2_gfl,          ind_hc2_agfl,          &
                                        surf_usm_v(l)%ground_level(m) )
              ind_hc2_win      = MERGE( ind_hc2_win_gfl,      ind_hc2_win_agfl,      &
                                        surf_usm_v(l)%ground_level(m) )
              ind_hc3          = MERGE( ind_hc3_gfl,          ind_hc3_agfl,          &
                                        surf_usm_v(l)%ground_level(m) )
              ind_hc3_win      = MERGE( ind_hc3_win_gfl,      ind_hc3_win_agfl,      &
                                        surf_usm_v(l)%ground_level(m) )
              ind_tc1          = MERGE( ind_tc1_gfl,          ind_tc1_agfl,          &
                                        surf_usm_v(l)%ground_level(m) )
              ind_tc1_win      = MERGE( ind_tc1_win_gfl,      ind_tc1_win_agfl,      &
                                        surf_usm_v(l)%ground_level(m) )
              ind_tc2          = MERGE( ind_tc2_gfl,          ind_tc2_agfl,          &
                                        surf_usm_v(l)%ground_level(m) )
              ind_tc2_win      = MERGE( ind_tc2_win_gfl,      ind_tc2_win_agfl,      &
                                        surf_usm_v(l)%ground_level(m) )
              ind_tc3          = MERGE( ind_tc3_gfl,          ind_tc3_agfl,          &
                                        surf_usm_v(l)%ground_level(m) )
              ind_tc3_win      = MERGE( ind_tc3_win_gfl,      ind_tc3_win_agfl,      &
                                        surf_usm_v(l)%ground_level(m) )
              ind_thick_1      = MERGE( ind_thick_1_gfl,      ind_thick_1_agfl,      &
                                        surf_usm_v(l)%ground_level(m) )
              ind_thick_1_win  = MERGE( ind_thick_1_win_gfl,  ind_thick_1_win_agfl,  &
                                        surf_usm_v(l)%ground_level(m) )
              ind_thick_2      = MERGE( ind_thick_2_gfl,      ind_thick_2_agfl,      &
                                        surf_usm_v(l)%ground_level(m) )
              ind_thick_2_win  = MERGE( ind_thick_2_win_gfl,  ind_thick_2_win_agfl,  &
                                        surf_usm_v(l)%ground_level(m) )
              ind_thick_3      = MERGE( ind_thick_3_gfl,      ind_thick_3_agfl,      &
                                        surf_usm_v(l)%ground_level(m) )
              ind_thick_3_win  = MERGE( ind_thick_3_win_gfl,  ind_thick_3_win_agfl,  &
                                        surf_usm_v(l)%ground_level(m) )
              ind_thick_4      = MERGE( ind_thick_4_gfl,      ind_thick_4_agfl,      &
                                        surf_usm_v(l)%ground_level(m) )
              ind_thick_4_win  = MERGE( ind_thick_4_win_gfl,  ind_thick_4_win_agfl,  &
                                        surf_usm_v(l)%ground_level(m) )
              ind_emis_wall    = MERGE( ind_emis_wall_gfl,    ind_emis_wall_agfl,    &
                                        surf_usm_v(l)%ground_level(m) )
              ind_emis_green   = MERGE( ind_emis_green_gfl,   ind_emis_green_agfl,   &
                                        surf_usm_v(l)%ground_level(m) )
              ind_emis_win     = MERGE( ind_emis_win_gfl,     ind_emis_win_agfl,     &
                                        surf_usm_v(l)%ground_level(m) )
              ind_trans        = MERGE( ind_trans_gfl,       ind_trans_agfl,         &
                                        surf_usm_v(l)%ground_level(m) )
              ind_z0           = MERGE( ind_z0_gfl,           ind_z0_agfl,           &
                                        surf_usm_v(l)%ground_level(m) )
              ind_z0qh         = MERGE( ind_z0qh_gfl,         ind_z0qh_agfl,         &
                                        surf_usm_v(l)%ground_level(m) )
!
!--           Store building type and its name on each surface element
              surf_usm_v(l)%building_type(m)      = building_type
              surf_usm_v(l)%building_type_name(m) = building_type_name(building_type)
!
!--           Initialize relatvie wall- (0), green- (1) and window (2) fractions
              surf_usm_v(l)%frac(m,ind_veg_wall)   = building_pars(ind_wall_frac,building_type)   
              surf_usm_v(l)%frac(m,ind_pav_green)  = building_pars(ind_green_frac_w,building_type) 
              surf_usm_v(l)%frac(m,ind_wat_win)    = building_pars(ind_win_frac,building_type)  
              surf_usm_v(l)%lai(m)                 = building_pars(ind_lai_w,building_type)  

              surf_usm_v(l)%rho_c_wall(nzb_wall,m)   = building_pars(ind_hc1,building_type)  
              surf_usm_v(l)%rho_c_wall(nzb_wall+1,m) = building_pars(ind_hc1,building_type)
              surf_usm_v(l)%rho_c_wall(nzb_wall+2,m) = building_pars(ind_hc2,building_type)
              surf_usm_v(l)%rho_c_wall(nzb_wall+3,m) = building_pars(ind_hc3,building_type)    
              
              surf_usm_v(l)%rho_c_green(nzb_wall,m)   = rho_c_soil !building_pars(ind_hc1,building_type)  
              surf_usm_v(l)%rho_c_green(nzb_wall+1,m) = rho_c_soil !building_pars(ind_hc1,building_type)
              surf_usm_v(l)%rho_c_green(nzb_wall+2,m) = rho_c_soil !building_pars(ind_hc2,building_type)
              surf_usm_v(l)%rho_c_green(nzb_wall+3,m) = rho_c_soil !building_pars(ind_hc3,building_type)    
              
              surf_usm_v(l)%rho_c_window(nzb_wall,m)   = building_pars(ind_hc1_win,building_type)  
              surf_usm_v(l)%rho_c_window(nzb_wall+1,m) = building_pars(ind_hc1_win,building_type)
              surf_usm_v(l)%rho_c_window(nzb_wall+2,m) = building_pars(ind_hc2_win,building_type)
              surf_usm_v(l)%rho_c_window(nzb_wall+3,m) = building_pars(ind_hc3_win,building_type)    

              surf_usm_v(l)%lambda_h(nzb_wall,m)   = building_pars(ind_tc1,building_type)  
              surf_usm_v(l)%lambda_h(nzb_wall+1,m) = building_pars(ind_tc1,building_type) 
              surf_usm_v(l)%lambda_h(nzb_wall+2,m) = building_pars(ind_tc2,building_type)
              surf_usm_v(l)%lambda_h(nzb_wall+3,m) = building_pars(ind_tc3,building_type)    
              
              surf_usm_v(l)%lambda_h_green(nzb_wall,m)   = lambda_h_green_sm !building_pars(ind_tc1,building_type)  
              surf_usm_v(l)%lambda_h_green(nzb_wall+1,m) = lambda_h_green_sm !building_pars(ind_tc1,building_type) 
              surf_usm_v(l)%lambda_h_green(nzb_wall+2,m) = lambda_h_green_sm !building_pars(ind_tc2,building_type)
              surf_usm_v(l)%lambda_h_green(nzb_wall+3,m) = lambda_h_green_sm !building_pars(ind_tc3,building_type)    

              surf_usm_v(l)%lambda_h_window(nzb_wall,m)   = building_pars(ind_tc1_win,building_type)  
              surf_usm_v(l)%lambda_h_window(nzb_wall+1,m) = building_pars(ind_tc1_win,building_type) 
              surf_usm_v(l)%lambda_h_window(nzb_wall+2,m) = building_pars(ind_tc2_win,building_type)
              surf_usm_v(l)%lambda_h_window(nzb_wall+3,m) = building_pars(ind_tc3_win,building_type)    

              surf_usm_v(l)%target_temp_summer(m)  = building_pars(ind_indoor_target_temp_summer,building_type)    
              surf_usm_v(l)%target_temp_winter(m)  = building_pars(ind_indoor_target_temp_winter,building_type)    
!
!--           emissivity of wall-, green- and window fraction 
              surf_usm_v(l)%emissivity(m,ind_veg_wall)  = building_pars(ind_emis_wall,building_type)
              surf_usm_v(l)%emissivity(m,ind_pav_green) = building_pars(ind_emis_green,building_type)
              surf_usm_v(l)%emissivity(m,ind_wat_win)   = building_pars(ind_emis_win,building_type)

              surf_usm_v(l)%transmissivity(m)      = building_pars(ind_trans,building_type)

              surf_usm_v(l)%z0(m)                  = building_pars(ind_z0,building_type)
              surf_usm_v(l)%z0h(m)                 = building_pars(ind_z0qh,building_type)
              surf_usm_v(l)%z0q(m)                 = building_pars(ind_z0qh,building_type)

              surf_usm_v(l)%albedo_type(m,ind_veg_wall)  = INT( building_pars(ind_alb_wall,building_type) )
              surf_usm_v(l)%albedo_type(m,ind_pav_green) = INT( building_pars(ind_alb_green,building_type) )
              surf_usm_v(l)%albedo_type(m,ind_wat_win)   = INT( building_pars(ind_alb_win,building_type) )

              surf_usm_v(l)%zw(nzb_wall,m)         = building_pars(ind_thick_1,building_type)
              surf_usm_v(l)%zw(nzb_wall+1,m)       = building_pars(ind_thick_2,building_type)
              surf_usm_v(l)%zw(nzb_wall+2,m)       = building_pars(ind_thick_3,building_type)
              surf_usm_v(l)%zw(nzb_wall+3,m)       = building_pars(ind_thick_4,building_type)
              
              surf_usm_v(l)%zw_green(nzb_wall,m)         = building_pars(ind_thick_1,building_type)
              surf_usm_v(l)%zw_green(nzb_wall+1,m)       = building_pars(ind_thick_2,building_type)
              surf_usm_v(l)%zw_green(nzb_wall+2,m)       = building_pars(ind_thick_3,building_type)
              surf_usm_v(l)%zw_green(nzb_wall+3,m)       = building_pars(ind_thick_4,building_type)

              surf_usm_v(l)%zw_window(nzb_wall,m)         = building_pars(ind_thick_1_win,building_type)
              surf_usm_v(l)%zw_window(nzb_wall+1,m)       = building_pars(ind_thick_2_win,building_type)
              surf_usm_v(l)%zw_window(nzb_wall+2,m)       = building_pars(ind_thick_3_win,building_type)
              surf_usm_v(l)%zw_window(nzb_wall+3,m)       = building_pars(ind_thick_4_win,building_type)

              surf_usm_v(l)%c_surface(m)           = building_pars(ind_c_surface,building_type)  
              surf_usm_v(l)%lambda_surf(m)         = building_pars(ind_lambda_surf,building_type)
              surf_usm_v(l)%c_surface_green(m)     = building_pars(ind_c_surface_green,building_type)  
              surf_usm_v(l)%lambda_surf_green(m)   = building_pars(ind_lambda_surf_green,building_type)
              surf_usm_v(l)%c_surface_window(m)    = building_pars(ind_c_surface_win,building_type)  
              surf_usm_v(l)%lambda_surf_window(m)  = building_pars(ind_lambda_surf_win,building_type)

           ENDDO
        ENDDO
!
!--     Level 2 - initialization via building type read from file
        IF ( building_type_f%from_file )  THEN
           DO  m = 1, surf_usm_h%ns
              i = surf_usm_h%i(m)
              j = surf_usm_h%j(m)
!
!--           For the moment, limit building type to 6 (to overcome errors in input file).
              st = building_type_f%var(j,i)
              IF ( st /= building_type_f%fill )  THEN

!
!--              In order to distinguish between ground floor level and 
!--              above-ground-floor level surfaces, set input indices.

                 ind_green_frac_r = MERGE( ind_green_frac_r_gfl, ind_green_frac_r_agfl, &
                                           surf_usm_h%ground_level(m) )
                 ind_lai_r        = MERGE( ind_lai_r_gfl,        ind_lai_r_agfl,        &
                                           surf_usm_h%ground_level(m) )
                 ind_z0           = MERGE( ind_z0_gfl,           ind_z0_agfl,           &
                                           surf_usm_h%ground_level(m) )
                 ind_z0qh         = MERGE( ind_z0qh_gfl,         ind_z0qh_agfl,         &
                                           surf_usm_h%ground_level(m) )
!
!--              Store building type and its name on each surface element
                 surf_usm_h%building_type(m)      = st
                 surf_usm_h%building_type_name(m) = building_type_name(st)
!
!--              Initialize relatvie wall- (0), green- (1) and window (2) fractions
                 surf_usm_h%frac(m,ind_veg_wall)  = building_pars(ind_wall_frac_r,st)   
                 surf_usm_h%frac(m,ind_pav_green) = building_pars(ind_green_frac_r,st)  
                 surf_usm_h%frac(m,ind_wat_win)   = building_pars(ind_win_frac_r,st)  
                 surf_usm_h%lai(m)                = building_pars(ind_lai_r,st)  

                 surf_usm_h%rho_c_wall(nzb_wall,m)   = building_pars(ind_hc1_wall_r,st)  
                 surf_usm_h%rho_c_wall(nzb_wall+1,m) = building_pars(ind_hc1_wall_r,st)
                 surf_usm_h%rho_c_wall(nzb_wall+2,m) = building_pars(ind_hc2_wall_r,st)
                 surf_usm_h%rho_c_wall(nzb_wall+3,m) = building_pars(ind_hc3_wall_r,st)    
                 surf_usm_h%lambda_h(nzb_wall,m)   = building_pars(ind_tc1_wall_r,st)  
                 surf_usm_h%lambda_h(nzb_wall+1,m) = building_pars(ind_tc1_wall_r,st) 
                 surf_usm_h%lambda_h(nzb_wall+2,m) = building_pars(ind_tc2_wall_r,st)
                 surf_usm_h%lambda_h(nzb_wall+3,m) = building_pars(ind_tc3_wall_r,st)    
                 
                 surf_usm_h%rho_c_green(nzb_wall,m)   = rho_c_soil !building_pars(ind_hc1_wall_r,st)  
                 surf_usm_h%rho_c_green(nzb_wall+1,m) = rho_c_soil !building_pars(ind_hc1_wall_r,st)
                 surf_usm_h%rho_c_green(nzb_wall+2,m) = rho_c_soil !building_pars(ind_hc2_wall_r,st)
                 surf_usm_h%rho_c_green(nzb_wall+3,m) = rho_c_soil !building_pars(ind_hc3_wall_r,st)    
                 surf_usm_h%lambda_h_green(nzb_wall,m)   = lambda_h_green_sm !building_pars(ind_tc1_wall_r,st)  
                 surf_usm_h%lambda_h_green(nzb_wall+1,m) = lambda_h_green_sm !building_pars(ind_tc1_wall_r,st) 
                 surf_usm_h%lambda_h_green(nzb_wall+2,m) = lambda_h_green_sm !building_pars(ind_tc2_wall_r,st)
                 surf_usm_h%lambda_h_green(nzb_wall+3,m) = lambda_h_green_sm !building_pars(ind_tc3_wall_r,st)    
                
                 surf_usm_h%rho_c_window(nzb_wall,m)   = building_pars(ind_hc1_win_r,st)  
                 surf_usm_h%rho_c_window(nzb_wall+1,m) = building_pars(ind_hc1_win_r,st)
                 surf_usm_h%rho_c_window(nzb_wall+2,m) = building_pars(ind_hc2_win_r,st)
                 surf_usm_h%rho_c_window(nzb_wall+3,m) = building_pars(ind_hc3_win_r,st)    
                 surf_usm_h%lambda_h_window(nzb_wall,m)   = building_pars(ind_tc1_win_r,st)  
                 surf_usm_h%lambda_h_window(nzb_wall+1,m) = building_pars(ind_tc1_win_r,st) 
                 surf_usm_h%lambda_h_window(nzb_wall+2,m) = building_pars(ind_tc2_win_r,st)
                 surf_usm_h%lambda_h_window(nzb_wall+3,m) = building_pars(ind_tc3_win_r,st)    

                 surf_usm_h%target_temp_summer(m)  = building_pars(ind_indoor_target_temp_summer,st)    
                 surf_usm_h%target_temp_winter(m)  = building_pars(ind_indoor_target_temp_winter,st)    
!
!--              emissivity of wall-, green- and window fraction 
                 surf_usm_h%emissivity(m,ind_veg_wall)  = building_pars(ind_emis_wall_r,st)
                 surf_usm_h%emissivity(m,ind_pav_green) = building_pars(ind_emis_green_r,st)
                 surf_usm_h%emissivity(m,ind_wat_win)   = building_pars(ind_emis_win_r,st)

                 surf_usm_h%transmissivity(m)      = building_pars(ind_trans_r,st)

                 surf_usm_h%z0(m)                  = building_pars(ind_z0,st)
                 surf_usm_h%z0h(m)                 = building_pars(ind_z0qh,st)
                 surf_usm_h%z0q(m)                 = building_pars(ind_z0qh,st)
!
!--              albedo type for wall fraction, green fraction, window fraction
                 surf_usm_h%albedo_type(m,ind_veg_wall)  = INT( building_pars(ind_alb_wall_r,st) )
                 surf_usm_h%albedo_type(m,ind_pav_green) = INT( building_pars(ind_alb_green_r,st) )
                 surf_usm_h%albedo_type(m,ind_wat_win)   = INT( building_pars(ind_alb_win_r,st) )

                 surf_usm_h%zw(nzb_wall,m)         = building_pars(ind_thick_1_wall_r,st)
                 surf_usm_h%zw(nzb_wall+1,m)       = building_pars(ind_thick_2_wall_r,st)
                 surf_usm_h%zw(nzb_wall+2,m)       = building_pars(ind_thick_3_wall_r,st)
                 surf_usm_h%zw(nzb_wall+3,m)       = building_pars(ind_thick_4_wall_r,st)
                 
                 surf_usm_h%zw_green(nzb_wall,m)         = building_pars(ind_thick_1_wall_r,st)
                 surf_usm_h%zw_green(nzb_wall+1,m)       = building_pars(ind_thick_2_wall_r,st)
                 surf_usm_h%zw_green(nzb_wall+2,m)       = building_pars(ind_thick_3_wall_r,st)
                 surf_usm_h%zw_green(nzb_wall+3,m)       = building_pars(ind_thick_4_wall_r,st)

                 surf_usm_h%zw_window(nzb_wall,m)         = building_pars(ind_thick_1_win_r,st)
                 surf_usm_h%zw_window(nzb_wall+1,m)       = building_pars(ind_thick_2_win_r,st)
                 surf_usm_h%zw_window(nzb_wall+2,m)       = building_pars(ind_thick_3_win_r,st)
                 surf_usm_h%zw_window(nzb_wall+3,m)       = building_pars(ind_thick_4_win_r,st)

                 surf_usm_h%c_surface(m)           = building_pars(ind_c_surface,st)  
                 surf_usm_h%lambda_surf(m)         = building_pars(ind_lambda_surf,st)
                 surf_usm_h%c_surface_green(m)     = building_pars(ind_c_surface_green,st)  
                 surf_usm_h%lambda_surf_green(m)   = building_pars(ind_lambda_surf_green,st)
                 surf_usm_h%c_surface_window(m)    = building_pars(ind_c_surface_win,st)  
                 surf_usm_h%lambda_surf_window(m)  = building_pars(ind_lambda_surf_win,st)
                 
                 surf_usm_h%green_type_roof(m)     = building_pars(ind_green_type_roof,st)

              ENDIF
           ENDDO

           DO  l = 0, 3
              DO  m = 1, surf_usm_v(l)%ns
                 i = surf_usm_v(l)%i(m) + surf_usm_v(l)%ioff
                 j = surf_usm_v(l)%j(m) + surf_usm_v(l)%joff
!
!--              For the moment, limit building type to 6 (to overcome errors in input file).

                 st = building_type_f%var(j,i)
                 IF ( st /= building_type_f%fill )  THEN

!
!--                 In order to distinguish between ground floor level and 
!--                 above-ground-floor level surfaces, set input indices.
                    ind_alb_green    = MERGE( ind_alb_green_gfl,    ind_alb_green_agfl,    &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_alb_wall     = MERGE( ind_alb_wall_gfl,     ind_alb_wall_agfl,     &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_alb_win      = MERGE( ind_alb_win_gfl,      ind_alb_win_agfl,      &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_wall_frac    = MERGE( ind_wall_frac_gfl,    ind_wall_frac_agfl,    &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_win_frac     = MERGE( ind_win_frac_gfl,     ind_win_frac_agfl,     &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_green_frac_w = MERGE( ind_green_frac_w_gfl, ind_green_frac_w_agfl, &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_green_frac_r = MERGE( ind_green_frac_r_gfl, ind_green_frac_r_agfl, &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_lai_r        = MERGE( ind_lai_r_gfl,        ind_lai_r_agfl,        &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_lai_w        = MERGE( ind_lai_w_gfl,        ind_lai_w_agfl,        &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_hc1          = MERGE( ind_hc1_gfl,          ind_hc1_agfl,          &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_hc1_win      = MERGE( ind_hc1_win_gfl,      ind_hc1_win_agfl,      &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_hc2          = MERGE( ind_hc2_gfl,          ind_hc2_agfl,          &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_hc2_win      = MERGE( ind_hc2_win_gfl,      ind_hc2_win_agfl,      &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_hc3          = MERGE( ind_hc3_gfl,          ind_hc3_agfl,          &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_hc3_win      = MERGE( ind_hc3_win_gfl,      ind_hc3_win_agfl,      &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_tc1          = MERGE( ind_tc1_gfl,          ind_tc1_agfl,          &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_tc1_win      = MERGE( ind_tc1_win_gfl,      ind_tc1_win_agfl,      &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_tc2          = MERGE( ind_tc2_gfl,          ind_tc2_agfl,          &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_tc2_win      = MERGE( ind_tc2_win_gfl,      ind_tc2_win_agfl,      &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_tc3          = MERGE( ind_tc3_gfl,          ind_tc3_agfl,          &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_tc3_win      = MERGE( ind_tc3_win_gfl,      ind_tc3_win_agfl,      &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_thick_1      = MERGE( ind_thick_1_gfl,      ind_thick_1_agfl,      &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_thick_1_win  = MERGE( ind_thick_1_win_gfl,  ind_thick_1_win_agfl,  &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_thick_2      = MERGE( ind_thick_2_gfl,      ind_thick_2_agfl,      &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_thick_2_win  = MERGE( ind_thick_2_win_gfl,  ind_thick_2_win_agfl,  &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_thick_3      = MERGE( ind_thick_3_gfl,      ind_thick_3_agfl,      &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_thick_3_win  = MERGE( ind_thick_3_win_gfl,  ind_thick_3_win_agfl,  &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_thick_4      = MERGE( ind_thick_4_gfl,      ind_thick_4_agfl,      &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_thick_4_win  = MERGE( ind_thick_4_win_gfl,  ind_thick_4_win_agfl,  &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_emis_wall    = MERGE( ind_emis_wall_gfl,    ind_emis_wall_agfl,    &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_emis_green   = MERGE( ind_emis_green_gfl,   ind_emis_green_agfl,   &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_emis_win     = MERGE( ind_emis_win_gfl,     ind_emis_win_agfl,     &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_trans        = MERGE( ind_trans_gfl,       ind_trans_agfl,         &
                                            surf_usm_v(l)%ground_level(m) )
                    ind_z0           = MERGE( ind_z0_gfl,           ind_z0_agfl,           &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_z0qh         = MERGE( ind_z0qh_gfl,         ind_z0qh_agfl,         &
                                              surf_usm_v(l)%ground_level(m) )
!
!--                 Store building type and its name on each surface element
                    surf_usm_v(l)%building_type(m)      = st
                    surf_usm_v(l)%building_type_name(m) = building_type_name(st)
!
!--                 Initialize relatvie wall- (0), green- (1) and window (2) fractions
                    surf_usm_v(l)%frac(m,ind_veg_wall)  = building_pars(ind_wall_frac,st)   
                    surf_usm_v(l)%frac(m,ind_pav_green) = building_pars(ind_green_frac_w,st)  
                    surf_usm_v(l)%frac(m,ind_wat_win)   = building_pars(ind_win_frac,st)   
                    surf_usm_v(l)%lai(m)                = building_pars(ind_lai_w,st)  

                    surf_usm_v(l)%rho_c_wall(nzb_wall,m)   = building_pars(ind_hc1,st)  
                    surf_usm_v(l)%rho_c_wall(nzb_wall+1,m) = building_pars(ind_hc1,st)
                    surf_usm_v(l)%rho_c_wall(nzb_wall+2,m) = building_pars(ind_hc2,st)
                    surf_usm_v(l)%rho_c_wall(nzb_wall+3,m) = building_pars(ind_hc3,st)
                    
                    surf_usm_v(l)%rho_c_green(nzb_wall,m)   = rho_c_soil !building_pars(ind_hc1,st)  
                    surf_usm_v(l)%rho_c_green(nzb_wall+1,m) = rho_c_soil !building_pars(ind_hc1,st)
                    surf_usm_v(l)%rho_c_green(nzb_wall+2,m) = rho_c_soil !building_pars(ind_hc2,st)
                    surf_usm_v(l)%rho_c_green(nzb_wall+3,m) = rho_c_soil !building_pars(ind_hc3,st)
                    
                    surf_usm_v(l)%rho_c_window(nzb_wall,m)   = building_pars(ind_hc1_win,st)  
                    surf_usm_v(l)%rho_c_window(nzb_wall+1,m) = building_pars(ind_hc1_win,st)
                    surf_usm_v(l)%rho_c_window(nzb_wall+2,m) = building_pars(ind_hc2_win,st)
                    surf_usm_v(l)%rho_c_window(nzb_wall+3,m) = building_pars(ind_hc3_win,st)

                    surf_usm_v(l)%lambda_h(nzb_wall,m)   = building_pars(ind_tc1,st)  
                    surf_usm_v(l)%lambda_h(nzb_wall+1,m) = building_pars(ind_tc1,st) 
                    surf_usm_v(l)%lambda_h(nzb_wall+2,m) = building_pars(ind_tc2,st)
                    surf_usm_v(l)%lambda_h(nzb_wall+3,m) = building_pars(ind_tc3,st) 
                    
                    surf_usm_v(l)%lambda_h_green(nzb_wall,m)   = lambda_h_green_sm !building_pars(ind_tc1,st)  
                    surf_usm_v(l)%lambda_h_green(nzb_wall+1,m) = lambda_h_green_sm !building_pars(ind_tc1,st) 
                    surf_usm_v(l)%lambda_h_green(nzb_wall+2,m) = lambda_h_green_sm !building_pars(ind_tc2,st)
                    surf_usm_v(l)%lambda_h_green(nzb_wall+3,m) = lambda_h_green_sm !building_pars(ind_tc3,st) 
                    
                    surf_usm_v(l)%lambda_h_window(nzb_wall,m)   = building_pars(ind_tc1_win,st)  
                    surf_usm_v(l)%lambda_h_window(nzb_wall+1,m) = building_pars(ind_tc1_win,st) 
                    surf_usm_v(l)%lambda_h_window(nzb_wall+2,m) = building_pars(ind_tc2_win,st)
                    surf_usm_v(l)%lambda_h_window(nzb_wall+3,m) = building_pars(ind_tc3_win,st) 

                    surf_usm_v(l)%target_temp_summer(m)  = building_pars(ind_indoor_target_temp_summer,st)    
                    surf_usm_v(l)%target_temp_winter(m)  = building_pars(ind_indoor_target_temp_winter,st)    
!
!--                 emissivity of wall-, green- and window fraction 
                    surf_usm_v(l)%emissivity(m,ind_veg_wall)  = building_pars(ind_emis_wall,st)
                    surf_usm_v(l)%emissivity(m,ind_pav_green) = building_pars(ind_emis_green,st)
                    surf_usm_v(l)%emissivity(m,ind_wat_win)   = building_pars(ind_emis_win,st)

                    surf_usm_v(l)%transmissivity(m)      = building_pars(ind_trans,st)

                    surf_usm_v(l)%z0(m)                  = building_pars(ind_z0,st)
                    surf_usm_v(l)%z0h(m)                 = building_pars(ind_z0qh,st)
                    surf_usm_v(l)%z0q(m)                 = building_pars(ind_z0qh,st)

                    surf_usm_v(l)%albedo_type(m,ind_veg_wall)  = INT( building_pars(ind_alb_wall,st) )
                    surf_usm_v(l)%albedo_type(m,ind_pav_green) = INT( building_pars(ind_alb_green,st) )
                    surf_usm_v(l)%albedo_type(m,ind_wat_win)   = INT( building_pars(ind_alb_win,st) )

                    surf_usm_v(l)%zw(nzb_wall,m)         = building_pars(ind_thick_1,st)
                    surf_usm_v(l)%zw(nzb_wall+1,m)       = building_pars(ind_thick_2,st)
                    surf_usm_v(l)%zw(nzb_wall+2,m)       = building_pars(ind_thick_3,st)
                    surf_usm_v(l)%zw(nzb_wall+3,m)       = building_pars(ind_thick_4,st)
                    
                    surf_usm_v(l)%zw_green(nzb_wall,m)         = building_pars(ind_thick_1,st)
                    surf_usm_v(l)%zw_green(nzb_wall+1,m)       = building_pars(ind_thick_2,st)
                    surf_usm_v(l)%zw_green(nzb_wall+2,m)       = building_pars(ind_thick_3,st)
                    surf_usm_v(l)%zw_green(nzb_wall+3,m)       = building_pars(ind_thick_4,st)
                    
                    surf_usm_v(l)%zw_window(nzb_wall,m)         = building_pars(ind_thick_1_win,st)
                    surf_usm_v(l)%zw_window(nzb_wall+1,m)       = building_pars(ind_thick_2_win,st)
                    surf_usm_v(l)%zw_window(nzb_wall+2,m)       = building_pars(ind_thick_3_win,st)
                    surf_usm_v(l)%zw_window(nzb_wall+3,m)       = building_pars(ind_thick_4_win,st)

                    surf_usm_v(l)%c_surface(m)           = building_pars(ind_c_surface,st)  
                    surf_usm_v(l)%lambda_surf(m)         = building_pars(ind_lambda_surf,st) 
                    surf_usm_v(l)%c_surface_green(m)     = building_pars(ind_c_surface_green,st)  
                    surf_usm_v(l)%lambda_surf_green(m)   = building_pars(ind_lambda_surf_green,st) 
                    surf_usm_v(l)%c_surface_window(m)    = building_pars(ind_c_surface_win,st)  
                    surf_usm_v(l)%lambda_surf_window(m)  = building_pars(ind_lambda_surf_win,st) 


                 ENDIF
              ENDDO
           ENDDO
        ENDIF 
        
!
!--     Level 3 - initialization via building_pars read from file. Note, only
!--     variables that are also defined in the input-standard can be initialized
!--     via file. Other variables will be initialized on level 1 or 2. 
        IF ( building_pars_f%from_file )  THEN
           DO  m = 1, surf_usm_h%ns
              i = surf_usm_h%i(m)
              j = surf_usm_h%j(m)

!
!--           In order to distinguish between ground floor level and 
!--           above-ground-floor level surfaces, set input indices.
              ind_wall_frac    = MERGE( ind_wall_frac_gfl,                     &
                                        ind_wall_frac_agfl,                    &
                                        surf_usm_h%ground_level(m) )
              ind_green_frac_r = MERGE( ind_green_frac_r_gfl,                  &
                                        ind_green_frac_r_agfl,                 &
                                        surf_usm_h%ground_level(m) )
              ind_win_frac     = MERGE( ind_win_frac_gfl,                      &
                                        ind_win_frac_agfl,                     &
                                        surf_usm_h%ground_level(m) )
              ind_lai_r        = MERGE( ind_lai_r_gfl,                         &
                                        ind_lai_r_agfl,                        &
                                        surf_usm_h%ground_level(m) )
              ind_z0           = MERGE( ind_z0_gfl,                            &
                                        ind_z0_agfl,                           &
                                        surf_usm_h%ground_level(m) )
              ind_z0qh         = MERGE( ind_z0qh_gfl,                          &
                                        ind_z0qh_agfl,                         &
                                        surf_usm_h%ground_level(m) )
              ind_hc1          = MERGE( ind_hc1_gfl,                           &
                                        ind_hc1_agfl,                          &
                                        surf_usm_h%ground_level(m) )
              ind_hc2          = MERGE( ind_hc2_gfl,                           &
                                        ind_hc2_agfl,                          &
                                        surf_usm_h%ground_level(m) )
              ind_hc3          = MERGE( ind_hc3_gfl,                           &
                                        ind_hc3_agfl,                          &
                                        surf_usm_h%ground_level(m) )
              ind_tc1          = MERGE( ind_tc1_gfl,                           &
                                        ind_tc1_agfl,                          &
                                        surf_usm_h%ground_level(m) )
              ind_tc2          = MERGE( ind_tc2_gfl,                           &
                                        ind_tc2_agfl,                          &
                                        surf_usm_h%ground_level(m) )
              ind_tc3          = MERGE( ind_tc3_gfl,                           &
                                        ind_tc3_agfl,                          &
                                        surf_usm_h%ground_level(m) )
              ind_emis_wall    = MERGE( ind_emis_wall_gfl,                     &
                                        ind_emis_wall_agfl,                    &
                                        surf_usm_h%ground_level(m) )
              ind_emis_green   = MERGE( ind_emis_green_gfl,                    &
                                        ind_emis_green_agfl,                   &
                                        surf_usm_h%ground_level(m) )
              ind_emis_win     = MERGE( ind_emis_win_gfl,                      &
                                        ind_emis_win_agfl,                     &
                                        surf_usm_h%ground_level(m) )
              ind_trans        = MERGE( ind_trans_gfl,                         &
                                        ind_trans_agfl,                        &
                                        surf_usm_h%ground_level(m) )

!
!--           Initialize relatvie wall- (0), green- (1) and window (2) fractions
              IF ( building_pars_f%pars_xy(ind_wall_frac,j,i) /=               &
                   building_pars_f%fill )                                      &
                 surf_usm_h%frac(m,ind_veg_wall)  =                            &
                                    building_pars_f%pars_xy(ind_wall_frac,j,i)   
                 
              IF ( building_pars_f%pars_xy(ind_green_frac_r,j,i) /=            &          
                   building_pars_f%fill )                                      & 
                 surf_usm_h%frac(m,ind_pav_green) =                            &
                                    building_pars_f%pars_xy(ind_green_frac_r,j,i) 
                 
              IF ( building_pars_f%pars_xy(ind_win_frac,j,i) /=                &
                   building_pars_f%fill )                                      & 
                 surf_usm_h%frac(m,ind_wat_win)   =                            &
                                    building_pars_f%pars_xy(ind_win_frac,j,i)
 
              IF ( building_pars_f%pars_xy(ind_lai_r,j,i) /=                   &
                   building_pars_f%fill )                                      &
                 surf_usm_h%lai(m)  = building_pars_f%pars_xy(ind_lai_r,j,i)

              IF ( building_pars_f%pars_xy(ind_hc1,j,i) /=                     &
                   building_pars_f%fill )  THEN
                 surf_usm_h%rho_c_wall(nzb_wall,m)   =                         &
                                    building_pars_f%pars_xy(ind_hc1,j,i)  
                 surf_usm_h%rho_c_wall(nzb_wall+1,m) =                         &
                                    building_pars_f%pars_xy(ind_hc1,j,i)
              ENDIF
              
              
              IF ( building_pars_f%pars_xy(ind_hc2,j,i) /=                     &
                   building_pars_f%fill )                                      &
                 surf_usm_h%rho_c_wall(nzb_wall+2,m) =                         &
                                    building_pars_f%pars_xy(ind_hc2,j,i)
                 
              IF ( building_pars_f%pars_xy(ind_hc3,j,i) /=                     &
                   building_pars_f%fill )                                      &
                 surf_usm_h%rho_c_wall(nzb_wall+3,m) =                         &
                                    building_pars_f%pars_xy(ind_hc3,j,i)
                 
              IF ( building_pars_f%pars_xy(ind_hc1,j,i) /=                     &
                   building_pars_f%fill )  THEN 
                 surf_usm_h%rho_c_green(nzb_wall,m)   =                        &
                                    building_pars_f%pars_xy(ind_hc1,j,i)  
                 surf_usm_h%rho_c_green(nzb_wall+1,m) =                        &
                                    building_pars_f%pars_xy(ind_hc1,j,i)
              ENDIF
              IF ( building_pars_f%pars_xy(ind_hc2,j,i) /=                     &
                   building_pars_f%fill )                                      &
                 surf_usm_h%rho_c_green(nzb_wall+2,m) =                        &
                                    building_pars_f%pars_xy(ind_hc2,j,i)
                 
              IF ( building_pars_f%pars_xy(ind_hc3,j,i) /=                     &
                   building_pars_f%fill )                                      &
                 surf_usm_h%rho_c_green(nzb_wall+3,m) =                        &
                                    building_pars_f%pars_xy(ind_hc3,j,i)
                 
              IF ( building_pars_f%pars_xy(ind_hc1,j,i) /=                     &
                   building_pars_f%fill )  THEN 
                 surf_usm_h%rho_c_window(nzb_wall,m)   =                       &
                                    building_pars_f%pars_xy(ind_hc1,j,i)  
                 surf_usm_h%rho_c_window(nzb_wall+1,m) =                       &
                                    building_pars_f%pars_xy(ind_hc1,j,i)
              ENDIF
              IF ( building_pars_f%pars_xy(ind_hc2,j,i) /=                     &
                   building_pars_f%fill )                                      &
                 surf_usm_h%rho_c_window(nzb_wall+2,m) =                       &
                                    building_pars_f%pars_xy(ind_hc2,j,i)
                 
              IF ( building_pars_f%pars_xy(ind_hc3,j,i) /=                     &
                   building_pars_f%fill )                                      &
                 surf_usm_h%rho_c_window(nzb_wall+3,m) =                       &
                                    building_pars_f%pars_xy(ind_hc3,j,i)

              IF ( building_pars_f%pars_xy(ind_tc1,j,i) /=                     &
                   building_pars_f%fill )  THEN
                 surf_usm_h%lambda_h(nzb_wall,m)   =                           &
                                    building_pars_f%pars_xy(ind_tc1,j,i)         
                 surf_usm_h%lambda_h(nzb_wall+1,m) =                           &
                                    building_pars_f%pars_xy(ind_tc1,j,i)        
              ENDIF
              IF ( building_pars_f%pars_xy(ind_tc2,j,i) /=                     &  
                   building_pars_f%fill )                                      &
                 surf_usm_h%lambda_h(nzb_wall+2,m) =                           &
                                    building_pars_f%pars_xy(ind_tc2,j,i)
                 
              IF ( building_pars_f%pars_xy(ind_tc3,j,i) /=                     &
                   building_pars_f%fill )                                      & 
                 surf_usm_h%lambda_h(nzb_wall+3,m) =                           &
                                    building_pars_f%pars_xy(ind_tc3,j,i)    
                 
              IF ( building_pars_f%pars_xy(ind_tc1,j,i) /=                     &
                   building_pars_f%fill )  THEN
                 surf_usm_h%lambda_h_green(nzb_wall,m)   =                     &
                                     building_pars_f%pars_xy(ind_tc1,j,i)         
                 surf_usm_h%lambda_h_green(nzb_wall+1,m) =                     &
                                     building_pars_f%pars_xy(ind_tc1,j,i)        
              ENDIF
              IF ( building_pars_f%pars_xy(ind_tc2,j,i) /=                     & 
                   building_pars_f%fill )                                      &
                 surf_usm_h%lambda_h_green(nzb_wall+2,m) =                     &
                                    building_pars_f%pars_xy(ind_tc2,j,i)
                 
              IF ( building_pars_f%pars_xy(ind_tc3,j,i) /=                     &       
                   building_pars_f%fill )                                      &
                 surf_usm_h%lambda_h_green(nzb_wall+3,m) =                     &
                                    building_pars_f%pars_xy(ind_tc3,j,i)    
                 
              IF ( building_pars_f%pars_xy(ind_tc1,j,i) /=                     &
                   building_pars_f%fill )  THEN
                 surf_usm_h%lambda_h_window(nzb_wall,m)   =                    &
                                     building_pars_f%pars_xy(ind_tc1,j,i)         
                 surf_usm_h%lambda_h_window(nzb_wall+1,m) =                    &
                                     building_pars_f%pars_xy(ind_tc1,j,i)        
              ENDIF
              IF ( building_pars_f%pars_xy(ind_tc2,j,i) /=                     &     
                   building_pars_f%fill )                                      &
                 surf_usm_h%lambda_h_window(nzb_wall+2,m) =                    &
                                     building_pars_f%pars_xy(ind_tc2,j,i)
                 
              IF ( building_pars_f%pars_xy(ind_tc3,j,i) /=                     &   
                   building_pars_f%fill )                                      &
                 surf_usm_h%lambda_h_window(nzb_wall+3,m) =                    &
                                    building_pars_f%pars_xy(ind_tc3,j,i)    

              IF ( building_pars_f%pars_xy(ind_indoor_target_temp_summer,j,i) /=&            
                   building_pars_f%fill )                                      & 
                 surf_usm_h%target_temp_summer(m)  =                           &
                      building_pars_f%pars_xy(ind_indoor_target_temp_summer,j,i)   
              IF ( building_pars_f%pars_xy(ind_indoor_target_temp_winter,j,i) /=&            
                   building_pars_f%fill )                                      & 
                 surf_usm_h%target_temp_winter(m)  =                           &
                      building_pars_f%pars_xy(ind_indoor_target_temp_winter,j,i)   

              IF ( building_pars_f%pars_xy(ind_emis_wall,j,i) /=               &   
                   building_pars_f%fill )                                      &
                 surf_usm_h%emissivity(m,ind_veg_wall)  =                      &
                                    building_pars_f%pars_xy(ind_emis_wall,j,i)
                 
              IF ( building_pars_f%pars_xy(ind_emis_green,j,i) /=              &            
                   building_pars_f%fill )                                      &
                 surf_usm_h%emissivity(m,ind_pav_green) =                      &
                                     building_pars_f%pars_xy(ind_emis_green,j,i)
                 
              IF ( building_pars_f%pars_xy(ind_emis_win,j,i) /=                & 
                   building_pars_f%fill )                                      &
                 surf_usm_h%emissivity(m,ind_wat_win)   =                      &
                                     building_pars_f%pars_xy(ind_emis_win,j,i)
                 
              IF ( building_pars_f%pars_xy(ind_trans,j,i) /=                   &    
                   building_pars_f%fill )                                      &
                 surf_usm_h%transmissivity(m) =                                &
                                    building_pars_f%pars_xy(ind_trans,j,i)

              IF ( building_pars_f%pars_xy(ind_z0,j,i) /=                      &          
                   building_pars_f%fill )                                      &
                 surf_usm_h%z0(m) = building_pars_f%pars_xy(ind_z0,j,i)
                 
              IF ( building_pars_f%pars_xy(ind_z0qh,j,i) /=                    &            
                   building_pars_f%fill )                                      &
                 surf_usm_h%z0h(m) = building_pars_f%pars_xy(ind_z0qh,j,i)
              IF ( building_pars_f%pars_xy(ind_z0qh,j,i) /=                    &           
                   building_pars_f%fill )                                      &
                 surf_usm_h%z0q(m) = building_pars_f%pars_xy(ind_z0qh,j,i)

              IF ( building_pars_f%pars_xy(ind_alb_wall_agfl,j,i) /=           &          
                   building_pars_f%fill )                                      & 
                 surf_usm_h%albedo_type(m,ind_veg_wall)  =                     &
                                 building_pars_f%pars_xy(ind_alb_wall_agfl,j,i)
                 
              IF ( building_pars_f%pars_xy(ind_alb_green_agfl,j,i) /=          &           
                   building_pars_f%fill )                                      &
                 surf_usm_h%albedo_type(m,ind_pav_green) =                     &
                                building_pars_f%pars_xy(ind_alb_green_agfl,j,i)
              IF ( building_pars_f%pars_xy(ind_alb_win_agfl,j,i) /=            &          
                   building_pars_f%fill )                                      &
                 surf_usm_h%albedo_type(m,ind_wat_win)   =                     &
                                   building_pars_f%pars_xy(ind_alb_win_agfl,j,i)

              IF ( building_pars_f%pars_xy(ind_thick_1_agfl,j,i) /=            &          
                   building_pars_f%fill )                                      & 
                 surf_usm_h%zw(nzb_wall,m) =                                   &
                                  building_pars_f%pars_xy(ind_thick_1_agfl,j,i)
                 
              IF ( building_pars_f%pars_xy(ind_thick_2_agfl,j,i) /=            &          
                   building_pars_f%fill )                                      &
                 surf_usm_h%zw(nzb_wall+1,m) =                                 &
                                  building_pars_f%pars_xy(ind_thick_2_agfl,j,i)
                 
              IF ( building_pars_f%pars_xy(ind_thick_3_agfl,j,i) /=            &         
                   building_pars_f%fill )                                      &
                 surf_usm_h%zw(nzb_wall+2,m) =                                 &
                                  building_pars_f%pars_xy(ind_thick_3_agfl,j,i)
                 
                 
              IF ( building_pars_f%pars_xy(ind_thick_4_agfl,j,i) /=            &         
                   building_pars_f%fill )                                      & 
                 surf_usm_h%zw(nzb_wall+3,m) =                                 &
                                  building_pars_f%pars_xy(ind_thick_4_agfl,j,i)
                 
              IF ( building_pars_f%pars_xy(ind_thick_1_agfl,j,i) /=            &           
                   building_pars_f%fill )                                      &
                 surf_usm_h%zw_green(nzb_wall,m) =                             &
                                  building_pars_f%pars_xy(ind_thick_1_agfl,j,i)
                 
              IF ( building_pars_f%pars_xy(ind_thick_2_agfl,j,i) /=            &          
                   building_pars_f%fill )                                      &
                 surf_usm_h%zw_green(nzb_wall+1,m) =                           &
                                   building_pars_f%pars_xy(ind_thick_2_agfl,j,i)
                 
              IF ( building_pars_f%pars_xy(ind_thick_3_agfl,j,i) /=            &          
                   building_pars_f%fill )                                      & 
                 surf_usm_h%zw_green(nzb_wall+2,m) =                           &
                                   building_pars_f%pars_xy(ind_thick_3_agfl,j,i)
                 
              IF ( building_pars_f%pars_xy(ind_thick_4_agfl,j,i) /=            &          
                   building_pars_f%fill )                                      &
                 surf_usm_h%zw_green(nzb_wall+3,m) =                           &
                                   building_pars_f%pars_xy(ind_thick_4_agfl,j,i)

              IF ( building_pars_f%pars_xy(ind_c_surface,j,i) /=               &       
                   building_pars_f%fill )                                      & 
                 surf_usm_h%c_surface(m) =                                     &
                                    building_pars_f%pars_xy(ind_c_surface,j,i)
                 
              IF ( building_pars_f%pars_xy(ind_lambda_surf,j,i) /=             &        
                   building_pars_f%fill )                                      &
                 surf_usm_h%lambda_surf(m) =                                   &
                                    building_pars_f%pars_xy(ind_lambda_surf,j,i)
           ENDDO

           DO  l = 0, 3
              DO  m = 1, surf_usm_v(l)%ns
                 i = surf_usm_v(l)%i(m) + surf_usm_v(l)%ioff
                 j = surf_usm_v(l)%j(m) + surf_usm_v(l)%joff
                
!
!--                 In order to distinguish between ground floor level and 
!--                 above-ground-floor level surfaces, set input indices.
                    ind_wall_frac    = MERGE( ind_wall_frac_gfl,               &
                                              ind_wall_frac_agfl,              &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_green_frac_w = MERGE( ind_green_frac_w_gfl,            &
                                              ind_green_frac_w_agfl,           &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_win_frac     = MERGE( ind_win_frac_gfl,                &
                                              ind_win_frac_agfl,               &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_lai_w        = MERGE( ind_lai_w_gfl,                   &
                                              ind_lai_w_agfl,                  &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_z0           = MERGE( ind_z0_gfl,                      &
                                              ind_z0_agfl,                     &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_z0qh         = MERGE( ind_z0qh_gfl,                    &
                                              ind_z0qh_agfl,                   &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_hc1          = MERGE( ind_hc1_gfl,                     &
                                              ind_hc1_agfl,                    &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_hc2          = MERGE( ind_hc2_gfl,                     &
                                              ind_hc2_agfl,                    &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_hc3          = MERGE( ind_hc3_gfl,                     &
                                              ind_hc3_agfl,                    &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_tc1          = MERGE( ind_tc1_gfl,                     &
                                              ind_tc1_agfl,                    &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_tc2          = MERGE( ind_tc2_gfl,                     &
                                              ind_tc2_agfl,                    &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_tc3          = MERGE( ind_tc3_gfl,                     &
                                              ind_tc3_agfl,                    &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_emis_wall    = MERGE( ind_emis_wall_gfl,               &
                                              ind_emis_wall_agfl,              &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_emis_green   = MERGE( ind_emis_green_gfl,              &
                                              ind_emis_green_agfl,             &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_emis_win     = MERGE( ind_emis_win_gfl,                &
                                              ind_emis_win_agfl,               &
                                              surf_usm_v(l)%ground_level(m) )
                    ind_trans        = MERGE( ind_trans_gfl,                   &
                                              ind_trans_agfl,                  &
                                              surf_usm_v(l)%ground_level(m) )
                    
!                   
!--                 Initialize relatvie wall- (0), green- (1) and window (2) fractions
                    IF ( building_pars_f%pars_xy(ind_wall_frac,j,i) /=         &
                         building_pars_f%fill )                                &
                       surf_usm_v(l)%frac(m,ind_veg_wall)  =                   &
                                          building_pars_f%pars_xy(ind_wall_frac,j,i)   
                       
                    IF ( building_pars_f%pars_xy(ind_green_frac_w,j,i) /=      &          
                         building_pars_f%fill )                                & 
                       surf_usm_v(l)%frac(m,ind_pav_green) =                   &
                                  building_pars_f%pars_xy(ind_green_frac_w,j,i) 
                       
                    IF ( building_pars_f%pars_xy(ind_win_frac,j,i) /=          &
                         building_pars_f%fill )                                & 
                       surf_usm_v(l)%frac(m,ind_wat_win)   =                   &
                                       building_pars_f%pars_xy(ind_win_frac,j,i)
                    
                    IF ( building_pars_f%pars_xy(ind_lai_w,j,i) /=             &
                         building_pars_f%fill )                                &
                       surf_usm_v(l)%lai(m)  =                                 &
                                       building_pars_f%pars_xy(ind_lai_w,j,i)
                    
                    IF ( building_pars_f%pars_xy(ind_hc1,j,i) /=               &
                         building_pars_f%fill )  THEN
                       surf_usm_v(l)%rho_c_wall(nzb_wall,m)   =                &
                                          building_pars_f%pars_xy(ind_hc1,j,i) 
                       surf_usm_v(l)%rho_c_wall(nzb_wall+1,m) =                &
                                          building_pars_f%pars_xy(ind_hc1,j,i)
                    ENDIF
                    
                    
                    IF ( building_pars_f%pars_xy(ind_hc2,j,i) /=               &
                         building_pars_f%fill )                                &
                       surf_usm_v(l)%rho_c_wall(nzb_wall+2,m) =                &
                                          building_pars_f%pars_xy(ind_hc2,j,i)
                       
                    IF ( building_pars_f%pars_xy(ind_hc3,j,i) /=               &         
                         building_pars_f%fill )                                &
                       surf_usm_v(l)%rho_c_wall(nzb_wall+3,m) =                &
                                          building_pars_f%pars_xy(ind_hc3,j,i)
                       
                    IF ( building_pars_f%pars_xy(ind_hc1,j,i) /=               &
                         building_pars_f%fill )  THEN 
                       surf_usm_v(l)%rho_c_green(nzb_wall,m)   =               &
                                          building_pars_f%pars_xy(ind_hc1,j,i) 
                       surf_usm_v(l)%rho_c_green(nzb_wall+1,m) =               &
                                          building_pars_f%pars_xy(ind_hc1,j,i)
                    ENDIF
                    IF ( building_pars_f%pars_xy(ind_hc2,j,i) /=               &
                         building_pars_f%fill )                                &
                       surf_usm_v(l)%rho_c_green(nzb_wall+2,m) =               &
                                          building_pars_f%pars_xy(ind_hc2,j,i)
                       
                    IF ( building_pars_f%pars_xy(ind_hc3,j,i) /=               &
                         building_pars_f%fill )                                &
                       surf_usm_v(l)%rho_c_green(nzb_wall+3,m) =               &
                                          building_pars_f%pars_xy(ind_hc3,j,i)
                       
                    IF ( building_pars_f%pars_xy(ind_hc1,j,i) /=               &
                         building_pars_f%fill )  THEN 
                       surf_usm_v(l)%rho_c_window(nzb_wall,m)   =              &
                                          building_pars_f%pars_xy(ind_hc1,j,i) 
                       surf_usm_v(l)%rho_c_window(nzb_wall+1,m) =              &
                                          building_pars_f%pars_xy(ind_hc1,j,i)
                    ENDIF
                    IF ( building_pars_f%pars_xy(ind_hc2,j,i) /=               &
                         building_pars_f%fill )                                &
                       surf_usm_v(l)%rho_c_window(nzb_wall+2,m) =              &
                                          building_pars_f%pars_xy(ind_hc2,j,i)
                       
                    IF ( building_pars_f%pars_xy(ind_hc3,j,i) /=               &
                         building_pars_f%fill )                                &
                       surf_usm_v(l)%rho_c_window(nzb_wall+3,m) =              &
                                          building_pars_f%pars_xy(ind_hc3,j,i)
                    
                    IF ( building_pars_f%pars_xy(ind_tc1,j,i) /=               &
                         building_pars_f%fill )  THEN
                       surf_usm_v(l)%lambda_h(nzb_wall,m)   =                  &
                                          building_pars_f%pars_xy(ind_tc1,j,i)   
                       surf_usm_v(l)%lambda_h(nzb_wall+1,m) =                  &
                                          building_pars_f%pars_xy(ind_tc1,j,i)  
                    ENDIF
                    IF ( building_pars_f%pars_xy(ind_tc2,j,i) /=               &  
                         building_pars_f%fill )                                &
                       surf_usm_v(l)%lambda_h(nzb_wall+2,m) =                  &
                                          building_pars_f%pars_xy(ind_tc2,j,i)
                       
                    IF ( building_pars_f%pars_xy(ind_tc3,j,i) /=               &
                         building_pars_f%fill )                                & 
                       surf_usm_v(l)%lambda_h(nzb_wall+3,m) =                  &
                                          building_pars_f%pars_xy(ind_tc3,j,i) 
                       
                    IF ( building_pars_f%pars_xy(ind_tc1,j,i) /=               &
                         building_pars_f%fill )  THEN
                       surf_usm_v(l)%lambda_h_green(nzb_wall,m)   =            &
                                           building_pars_f%pars_xy(ind_tc1,j,i)   
                       surf_usm_v(l)%lambda_h_green(nzb_wall+1,m) =            &
                                           building_pars_f%pars_xy(ind_tc1,j,i)  
                    ENDIF
                    IF ( building_pars_f%pars_xy(ind_tc2,j,i) /=               & 
                         building_pars_f%fill )                                &
                       surf_usm_v(l)%lambda_h_green(nzb_wall+2,m) =            &
                                          building_pars_f%pars_xy(ind_tc2,j,i)
                       
                    IF ( building_pars_f%pars_xy(ind_tc3,j,i) /=               &       
                         building_pars_f%fill )                                &
                       surf_usm_v(l)%lambda_h_green(nzb_wall+3,m) =            &
                                          building_pars_f%pars_xy(ind_tc3,j,i) 
                       
                    IF ( building_pars_f%pars_xy(ind_tc1,j,i) /=         &
                         building_pars_f%fill )  THEN
                       surf_usm_v(l)%lambda_h_window(nzb_wall,m)   =           &
                                     building_pars_f%pars_xy(ind_tc1,j,i)         
                       surf_usm_v(l)%lambda_h_window(nzb_wall+1,m) =           &
                                     building_pars_f%pars_xy(ind_tc1,j,i)        
                    ENDIF
                    IF ( building_pars_f%pars_xy(ind_tc2,j,i) /=               &     
                         building_pars_f%fill )                                &
                       surf_usm_v(l)%lambda_h_window(nzb_wall+2,m) =           &
                                           building_pars_f%pars_xy(ind_tc2,j,i)
                       
                    IF ( building_pars_f%pars_xy(ind_tc3,j,i) /=               &   
                         building_pars_f%fill )                                &
                       surf_usm_v(l)%lambda_h_window(nzb_wall+3,m) =           &
                                          building_pars_f%pars_xy(ind_tc3,j,i)    
                    
                    IF ( building_pars_f%pars_xy(ind_indoor_target_temp_summer,j,i) /=&            
                         building_pars_f%fill )                                & 
                       surf_usm_v(l)%target_temp_summer(m)  =                  &
                            building_pars_f%pars_xy(ind_indoor_target_temp_summer,j,i)   
                    IF ( building_pars_f%pars_xy(ind_indoor_target_temp_winter,j,i) /=&            
                         building_pars_f%fill )                                & 
                       surf_usm_v(l)%target_temp_winter(m)  =                  &
                            building_pars_f%pars_xy(ind_indoor_target_temp_winter,j,i)   
                    
                    IF ( building_pars_f%pars_xy(ind_emis_wall,j,i) /=         &   
                         building_pars_f%fill )                                &
                       surf_usm_v(l)%emissivity(m,ind_veg_wall)  =             &
                                      building_pars_f%pars_xy(ind_emis_wall,j,i)
                       
                    IF ( building_pars_f%pars_xy(ind_emis_green,j,i) /=        &            
                         building_pars_f%fill )                                &
                       surf_usm_v(l)%emissivity(m,ind_pav_green) =             &
                                      building_pars_f%pars_xy(ind_emis_green,j,i)
                       
                    IF ( building_pars_f%pars_xy(ind_emis_win,j,i) /=          & 
                         building_pars_f%fill )                                &
                       surf_usm_v(l)%emissivity(m,ind_wat_win)   =             &
                                      building_pars_f%pars_xy(ind_emis_win,j,i)
                       
                    IF ( building_pars_f%pars_xy(ind_trans,j,i) /=             &    
                         building_pars_f%fill )                                &
                       surf_usm_v(l)%transmissivity(m) =                       &
                                          building_pars_f%pars_xy(ind_trans,j,i)
                    
                    IF ( building_pars_f%pars_xy(ind_z0,j,i) /=                &          
                         building_pars_f%fill )                                &
                       surf_usm_v(l)%z0(m) = building_pars_f%pars_xy(ind_z0,j,i)
                       
                    IF ( building_pars_f%pars_xy(ind_z0qh,j,i) /=              &            
                         building_pars_f%fill )                                &
                       surf_usm_v(l)%z0h(m) =                                  &
                                       building_pars_f%pars_xy(ind_z0qh,j,i)
                    IF ( building_pars_f%pars_xy(ind_z0qh,j,i) /=              &           
                         building_pars_f%fill )                                &
                       surf_usm_v(l)%z0q(m) =                                  &
                                       building_pars_f%pars_xy(ind_z0qh,j,i)
                    
                    IF ( building_pars_f%pars_xy(ind_alb_wall_agfl,j,i) /=     &          
                         building_pars_f%fill )                                & 
                       surf_usm_v(l)%albedo_type(m,ind_veg_wall)  =            &
                                 building_pars_f%pars_xy(ind_alb_wall_agfl,j,i)
                       
                    IF ( building_pars_f%pars_xy(ind_alb_green_agfl,j,i) /=    &           
                         building_pars_f%fill )                                &
                       surf_usm_v(l)%albedo_type(m,ind_pav_green) =            &
                                 building_pars_f%pars_xy(ind_alb_green_agfl,j,i)
                    IF ( building_pars_f%pars_xy(ind_alb_win_agfl,j,i) /=      &          
                         building_pars_f%fill )                                &
                       surf_usm_v(l)%albedo_type(m,ind_wat_win)   =            &
                                   building_pars_f%pars_xy(ind_alb_win_agfl,j,i)
                    
                    IF ( building_pars_f%pars_xy(ind_thick_1_agfl,j,i) /=      &          
                         building_pars_f%fill )                                & 
                       surf_usm_v(l)%zw(nzb_wall,m) =                          &
                                   building_pars_f%pars_xy(ind_thick_1_agfl,j,i)
                       
                    IF ( building_pars_f%pars_xy(ind_thick_2_agfl,j,i) /=      &          
                         building_pars_f%fill )                                &
                       surf_usm_v(l)%zw(nzb_wall+1,m) =                        &
                                   building_pars_f%pars_xy(ind_thick_2_agfl,j,i)
                       
                    IF ( building_pars_f%pars_xy(ind_thick_3_agfl,j,i) /=      &         
                         building_pars_f%fill )                                &
                       surf_usm_v(l)%zw(nzb_wall+2,m) =                        &
                                   building_pars_f%pars_xy(ind_thick_3_agfl,j,i)
                       
                       
                    IF ( building_pars_f%pars_xy(ind_thick_4_agfl,j,i) /=      &         
                         building_pars_f%fill )                                & 
                       surf_usm_v(l)%zw(nzb_wall+3,m) =                        &
                                   building_pars_f%pars_xy(ind_thick_4_agfl,j,i)
                       
                    IF ( building_pars_f%pars_xy(ind_thick_1_agfl,j,i) /=      &           
                         building_pars_f%fill )                                &
                       surf_usm_v(l)%zw_green(nzb_wall,m) =                    &
                                   building_pars_f%pars_xy(ind_thick_1_agfl,j,i)
                       
                    IF ( building_pars_f%pars_xy(ind_thick_2_agfl,j,i) /=      &          
                         building_pars_f%fill )                                &
                       surf_usm_v(l)%zw_green(nzb_wall+1,m) =                  &
                                   building_pars_f%pars_xy(ind_thick_2_agfl,j,i)
                       
                    IF ( building_pars_f%pars_xy(ind_thick_3_agfl,j,i) /=      &          
                         building_pars_f%fill )                                & 
                       surf_usm_v(l)%zw_green(nzb_wall+2,m) =                  &
                                   building_pars_f%pars_xy(ind_thick_3_agfl,j,i)
                       
                    IF ( building_pars_f%pars_xy(ind_thick_4_agfl,j,i) /=      &          
                         building_pars_f%fill )                                &
                       surf_usm_v(l)%zw_green(nzb_wall+3,m) =                  &
                                   building_pars_f%pars_xy(ind_thick_4_agfl,j,i)
                    
                    IF ( building_pars_f%pars_xy(ind_c_surface,j,i) /=         &       
                         building_pars_f%fill )                                & 
                       surf_usm_v(l)%c_surface(m) =                            &
                                     building_pars_f%pars_xy(ind_c_surface,j,i)
                       
                    IF ( building_pars_f%pars_xy(ind_lambda_surf,j,i) /=       &        
                         building_pars_f%fill )                                &
                       surf_usm_v(l)%lambda_surf(m) =                          &
                                    building_pars_f%pars_xy(ind_lambda_surf,j,i)
                    
              ENDDO
           ENDDO
        ENDIF 
!
!--     Read building surface pars. If present, they override LOD1-LOD3 building
!--     pars where applicable
        IF ( building_surface_pars_f%from_file )  THEN
           DO  m = 1, surf_usm_h%ns
              i = surf_usm_h%i(m)
              j = surf_usm_h%j(m)
              k = surf_usm_h%k(m)
!
!--           Iterate over surfaces in column, check height and orientation
              DO  is = building_surface_pars_f%index_ji(1,j,i), &
                       building_surface_pars_f%index_ji(2,j,i)
                 IF ( building_surface_pars_f%coords(4,is) == -surf_usm_h%koff .AND.            &
                      building_surface_pars_f%coords(1,is) == k )  THEN

                    IF ( building_surface_pars_f%pars(ind_s_wall_frac,is) /=                     &
                         building_surface_pars_f%fill )                                          &
                       surf_usm_h%frac(m,ind_veg_wall) =                                         &
                                building_surface_pars_f%pars(ind_s_wall_frac,is)

                    IF ( building_surface_pars_f%pars(ind_s_green_frac_w,is) /=                  &
                         building_surface_pars_f%fill )                                          &
                       surf_usm_h%frac(m,ind_pav_green) =                                        &
                                building_surface_pars_f%pars(ind_s_green_frac_w,is)

                    IF ( building_surface_pars_f%pars(ind_s_green_frac_r,is) /=                  &
                         building_surface_pars_f%fill )                                          &
                       surf_usm_h%frac(m,ind_pav_green) =                                        &
                                building_surface_pars_f%pars(ind_s_green_frac_r,is)
                                !TODO clarify: why should _w and _r be on the same surface?

                    IF ( building_surface_pars_f%pars(ind_s_win_frac,is) /=                      &
                         building_surface_pars_f%fill )                                          &
                       surf_usm_h%frac(m,ind_wat_win) =                                          &
                                building_surface_pars_f%pars(ind_s_win_frac,is)

                    IF ( building_surface_pars_f%pars(ind_s_lai_r,is) /=                         &
                         building_surface_pars_f%fill )                                          &
                       surf_usm_h%lai(m) =                                                       &
                                building_surface_pars_f%pars(ind_s_lai_r,is)

                    IF ( building_surface_pars_f%pars(ind_s_hc1,is) /=                           &
                         building_surface_pars_f%fill )  THEN
                       surf_usm_h%rho_c_wall(nzb_wall:nzb_wall+1,m) =                            &
                                building_surface_pars_f%pars(ind_s_hc1,is)
                       surf_usm_h%rho_c_green(nzb_wall:nzb_wall+1,m) =                           &
                                building_surface_pars_f%pars(ind_s_hc1,is)
                       surf_usm_h%rho_c_window(nzb_wall:nzb_wall+1,m) =                          &
                                building_surface_pars_f%pars(ind_s_hc1,is)
                    ENDIF

                    IF ( building_surface_pars_f%pars(ind_s_hc2,is) /=                           &
                         building_surface_pars_f%fill )  THEN
                       surf_usm_h%rho_c_wall(nzb_wall+2,m) =                                     &
                                building_surface_pars_f%pars(ind_s_hc2,is)
                       surf_usm_h%rho_c_green(nzb_wall+2,m) =                                    &
                                building_surface_pars_f%pars(ind_s_hc2,is)
                       surf_usm_h%rho_c_window(nzb_wall+2,m) =                                   &
                                building_surface_pars_f%pars(ind_s_hc2,is)
                    ENDIF

                    IF ( building_surface_pars_f%pars(ind_s_hc3,is) /=                           &
                         building_surface_pars_f%fill )  THEN
                       surf_usm_h%rho_c_wall(nzb_wall+3,m) =                                     &
                                building_surface_pars_f%pars(ind_s_hc3,is)
                       surf_usm_h%rho_c_green(nzb_wall+3,m) =                                    &
                                building_surface_pars_f%pars(ind_s_hc3,is)
                       surf_usm_h%rho_c_window(nzb_wall+3,m) =                                   &
                                building_surface_pars_f%pars(ind_s_hc3,is)
                    ENDIF

                    IF ( building_surface_pars_f%pars(ind_s_tc1,is) /=                           &
                         building_surface_pars_f%fill )  THEN
                       surf_usm_h%lambda_h(nzb_wall:nzb_wall+1,m) =                              &
                                building_surface_pars_f%pars(ind_s_tc1,is)
                       surf_usm_h%lambda_h_green(nzb_wall:nzb_wall+1,m) =                        &
                                building_surface_pars_f%pars(ind_s_tc1,is)
                       surf_usm_h%lambda_h_window(nzb_wall:nzb_wall+1,m) =                       &
                                building_surface_pars_f%pars(ind_s_tc1,is)
                    ENDIF

                    IF ( building_surface_pars_f%pars(ind_s_tc2,is) /=                           &
                         building_surface_pars_f%fill )  THEN
                       surf_usm_h%lambda_h(nzb_wall+2,m) =                                       &
                                building_surface_pars_f%pars(ind_s_tc2,is)
                       surf_usm_h%lambda_h_green(nzb_wall+2,m) =                                 &
                                building_surface_pars_f%pars(ind_s_tc2,is)
                       surf_usm_h%lambda_h_window(nzb_wall+2,m) =                                &
                                building_surface_pars_f%pars(ind_s_tc2,is)
                    ENDIF

                    IF ( building_surface_pars_f%pars(ind_s_tc3,is) /=                           &
                         building_surface_pars_f%fill )  THEN
                       surf_usm_h%lambda_h(nzb_wall+3,m) =                                       &
                                building_surface_pars_f%pars(ind_s_tc3,is)
                       surf_usm_h%lambda_h_green(nzb_wall+3,m) =                                 &
                                building_surface_pars_f%pars(ind_s_tc3,is)
                       surf_usm_h%lambda_h_window(nzb_wall+3,m) =                                &
                                building_surface_pars_f%pars(ind_s_tc3,is)
                    ENDIF

                    IF ( building_surface_pars_f%pars(ind_s_indoor_target_temp_summer,is) /=     &
                         building_surface_pars_f%fill )                                          &
                       surf_usm_h%target_temp_summer(m) =                                        &
                                building_surface_pars_f%pars(ind_s_indoor_target_temp_summer,is)

                    IF ( building_surface_pars_f%pars(ind_s_indoor_target_temp_winter,is) /=     &
                         building_surface_pars_f%fill )                                          &
                       surf_usm_h%target_temp_winter(m) =                                        &
                                building_surface_pars_f%pars(ind_s_indoor_target_temp_winter,is)

                    IF ( building_surface_pars_f%pars(ind_s_emis_wall,is) /=                     &
                         building_surface_pars_f%fill )                                          &
                       surf_usm_h%emissivity(m,ind_veg_wall) =                                   &
                                building_surface_pars_f%pars(ind_s_emis_wall,is)

                    IF ( building_surface_pars_f%pars(ind_s_emis_green,is) /=                    &
                         building_surface_pars_f%fill )                                          &
                       surf_usm_h%emissivity(m,ind_pav_green) =                                  &
                                building_surface_pars_f%pars(ind_s_emis_green,is)

                    IF ( building_surface_pars_f%pars(ind_s_emis_win,is) /=                      &
                         building_surface_pars_f%fill )                                          &
                       surf_usm_h%emissivity(m,ind_wat_win) =                                    &
                                building_surface_pars_f%pars(ind_s_emis_win,is)

                    IF ( building_surface_pars_f%pars(ind_s_trans,is) /=                         &
                         building_surface_pars_f%fill )                                          &
                       surf_usm_h%transmissivity(m) =                                            &
                                building_surface_pars_f%pars(ind_s_trans,is)

                    IF ( building_surface_pars_f%pars(ind_s_z0,is) /=                            &
                         building_surface_pars_f%fill )                                          &
                       surf_usm_h%z0(m) =                                                        &
                                building_surface_pars_f%pars(ind_s_z0,is)

                    IF ( building_surface_pars_f%pars(ind_s_z0qh,is) /=                          &
                         building_surface_pars_f%fill )  THEN
                       surf_usm_h%z0q(m) =                                                       &
                                building_surface_pars_f%pars(ind_s_z0qh,is)
                       surf_usm_h%z0h(m) =                                                       &
                                building_surface_pars_f%pars(ind_s_z0qh,is)
                    ENDIF

                    EXIT ! surface was found and processed
                 ENDIF
              ENDDO
           ENDDO

           DO  l = 0, 3
              DO  m = 1, surf_usm_v(l)%ns
                 i = surf_usm_v(l)%i(m)
                 j = surf_usm_v(l)%j(m)
                 k = surf_usm_v(l)%k(m)
!
!--              Iterate over surfaces in column, check height and orientation
                 DO  is = building_surface_pars_f%index_ji(1,j,i), &
                          building_surface_pars_f%index_ji(2,j,i)
                    IF ( building_surface_pars_f%coords(5,is) == -surf_usm_v(l)%joff .AND.    &
                         building_surface_pars_f%coords(6,is) == -surf_usm_v(l)%ioff .AND.    &
                         building_surface_pars_f%coords(1,is) == k )  THEN

                       IF ( building_surface_pars_f%pars(ind_s_wall_frac,is) /=                &
                            building_surface_pars_f%fill )                                     &
                          surf_usm_v(l)%frac(m,ind_veg_wall) =                                 &
                                   building_surface_pars_f%pars(ind_s_wall_frac,is)

                       IF ( building_surface_pars_f%pars(ind_s_green_frac_w,is) /=             &
                            building_surface_pars_f%fill )                                     &
                          surf_usm_v(l)%frac(m,ind_pav_green) =                                &
                                   building_surface_pars_f%pars(ind_s_green_frac_w,is)

                       IF ( building_surface_pars_f%pars(ind_s_green_frac_r,is) /=             &
                            building_surface_pars_f%fill )                                     &
                          surf_usm_v(l)%frac(m,ind_pav_green) =                                &
                                   building_surface_pars_f%pars(ind_s_green_frac_r,is)
                                   !TODO clarify: why should _w and _r be on the same surface?

                       IF ( building_surface_pars_f%pars(ind_s_win_frac,is) /=                 &
                            building_surface_pars_f%fill )                                     &
                          surf_usm_v(l)%frac(m,ind_wat_win) =                                  &
                                   building_surface_pars_f%pars(ind_s_win_frac,is)

                       IF ( building_surface_pars_f%pars(ind_s_lai_r,is) /=                    &
                            building_surface_pars_f%fill )                                     &
                          surf_usm_v(l)%lai(m) =                                               &
                                   building_surface_pars_f%pars(ind_s_lai_r,is)

                       IF ( building_surface_pars_f%pars(ind_s_hc1,is) /=                      &
                            building_surface_pars_f%fill )  THEN
                          surf_usm_v(l)%rho_c_wall(nzb_wall:nzb_wall+1,m) =                    &
                                   building_surface_pars_f%pars(ind_s_hc1,is)
                          surf_usm_v(l)%rho_c_green(nzb_wall:nzb_wall+1,m) =                   &
                                   building_surface_pars_f%pars(ind_s_hc1,is)
                          surf_usm_v(l)%rho_c_window(nzb_wall:nzb_wall+1,m) =                  &
                                   building_surface_pars_f%pars(ind_s_hc1,is)
                       ENDIF

                       IF ( building_surface_pars_f%pars(ind_s_hc2,is) /=                      &
                            building_surface_pars_f%fill )  THEN
                          surf_usm_v(l)%rho_c_wall(nzb_wall+2,m) =                             &
                                   building_surface_pars_f%pars(ind_s_hc2,is)
                          surf_usm_v(l)%rho_c_green(nzb_wall+2,m) =                            &
                                   building_surface_pars_f%pars(ind_s_hc2,is)
                          surf_usm_v(l)%rho_c_window(nzb_wall+2,m) =                           &
                                   building_surface_pars_f%pars(ind_s_hc2,is)
                       ENDIF

                       IF ( building_surface_pars_f%pars(ind_s_hc3,is) /=                      &
                            building_surface_pars_f%fill )  THEN
                          surf_usm_v(l)%rho_c_wall(nzb_wall+3,m) =                             &
                                   building_surface_pars_f%pars(ind_s_hc3,is)
                          surf_usm_v(l)%rho_c_green(nzb_wall+3,m) =                            &
                                   building_surface_pars_f%pars(ind_s_hc3,is)
                          surf_usm_v(l)%rho_c_window(nzb_wall+3,m) =                           &
                                   building_surface_pars_f%pars(ind_s_hc3,is)
                       ENDIF

                       IF ( building_surface_pars_f%pars(ind_s_tc1,is) /=                      &
                            building_surface_pars_f%fill )  THEN
                          surf_usm_v(l)%lambda_h(nzb_wall:nzb_wall+1,m) =                      &
                                   building_surface_pars_f%pars(ind_s_tc1,is)
                          surf_usm_v(l)%lambda_h_green(nzb_wall:nzb_wall+1,m) =                &
                                   building_surface_pars_f%pars(ind_s_tc1,is)
                          surf_usm_v(l)%lambda_h_window(nzb_wall:nzb_wall+1,m) =               &
                                   building_surface_pars_f%pars(ind_s_tc1,is)
                       ENDIF

                       IF ( building_surface_pars_f%pars(ind_s_tc2,is) /=                      &
                            building_surface_pars_f%fill )  THEN
                          surf_usm_v(l)%lambda_h(nzb_wall+2,m) =                               &
                                   building_surface_pars_f%pars(ind_s_tc2,is)
                          surf_usm_v(l)%lambda_h_green(nzb_wall+2,m) =                         &
                                   building_surface_pars_f%pars(ind_s_tc2,is)
                          surf_usm_v(l)%lambda_h_window(nzb_wall+2,m) =                        &
                                   building_surface_pars_f%pars(ind_s_tc2,is)
                       ENDIF

                       IF ( building_surface_pars_f%pars(ind_s_tc3,is) /=                      &
                            building_surface_pars_f%fill )  THEN
                          surf_usm_v(l)%lambda_h(nzb_wall+3,m) =                               &
                                   building_surface_pars_f%pars(ind_s_tc3,is)
                          surf_usm_v(l)%lambda_h_green(nzb_wall+3,m) =                         &
                                   building_surface_pars_f%pars(ind_s_tc3,is)
                          surf_usm_v(l)%lambda_h_window(nzb_wall+3,m) =                        &
                                   building_surface_pars_f%pars(ind_s_tc3,is)
                       ENDIF

                       IF ( building_surface_pars_f%pars(ind_s_indoor_target_temp_summer,is) /=    &
                            building_surface_pars_f%fill )                                         &
                          surf_usm_v(l)%target_temp_summer(m) =                                    &
                                   building_surface_pars_f%pars(ind_s_indoor_target_temp_summer,is)

                       IF ( building_surface_pars_f%pars(ind_s_indoor_target_temp_winter,is) /=    &
                            building_surface_pars_f%fill )                                         &
                          surf_usm_v(l)%target_temp_winter(m) =                                    &
                                   building_surface_pars_f%pars(ind_s_indoor_target_temp_winter,is)

                       IF ( building_surface_pars_f%pars(ind_s_emis_wall,is) /=      &
                            building_surface_pars_f%fill )                           &
                          surf_usm_v(l)%emissivity(m,ind_veg_wall) =                 &
                                   building_surface_pars_f%pars(ind_s_emis_wall,is)

                       IF ( building_surface_pars_f%pars(ind_s_emis_green,is) /=     &
                            building_surface_pars_f%fill )                           &
                          surf_usm_v(l)%emissivity(m,ind_pav_green) =                &
                                   building_surface_pars_f%pars(ind_s_emis_green,is)

                       IF ( building_surface_pars_f%pars(ind_s_emis_win,is) /=       &
                            building_surface_pars_f%fill )                           &
                          surf_usm_v(l)%emissivity(m,ind_wat_win) =                  &
                                   building_surface_pars_f%pars(ind_s_emis_win,is)

                       IF ( building_surface_pars_f%pars(ind_s_trans,is) /=          &
                            building_surface_pars_f%fill )                           &
                          surf_usm_v(l)%transmissivity(m) =                          &
                                   building_surface_pars_f%pars(ind_s_trans,is)

                       IF ( building_surface_pars_f%pars(ind_s_z0,is) /=             &
                            building_surface_pars_f%fill )                           &
                          surf_usm_v(l)%z0(m) =                                      &
                                   building_surface_pars_f%pars(ind_s_z0,is)

                       IF ( building_surface_pars_f%pars(ind_s_z0qh,is) /=           &
                            building_surface_pars_f%fill )  THEN
                          surf_usm_v(l)%z0q(m) =                                     &
                                   building_surface_pars_f%pars(ind_s_z0qh,is)
                          surf_usm_v(l)%z0h(m) =                                     &
                                   building_surface_pars_f%pars(ind_s_z0qh,is)
                       ENDIF

                       EXIT ! surface was found and processed
                    ENDIF
                 ENDDO
              ENDDO
           ENDDO
        ENDIF
!
!--     Run further checks to ensure that the respecitve material fractions are 
!--     prescribed properly. Start with horizontal surfaces (roofs).
        relative_fractions_corrected = .FALSE.
        DO  m = 1, surf_usm_h%ns
           sum_frac = SUM( surf_usm_h%frac(m,:) )
           IF ( sum_frac /= 1.0_wp )  THEN
              relative_fractions_corrected = .TRUE.
!
!--           Normalize relative fractions to 1. Deviations from 1 can
!--           arise, e.g. by rounding errors but also by inconsistent 
!--           driver creation. 
              IF ( sum_frac /= 0.0_wp )  THEN
                 surf_usm_h%frac(m,:) = surf_usm_h%frac(m,:) / sum_frac
!
!--           In case all relative fractions are erroneously set to zero,
!--           set wall fraction to 1.
              ELSE
                 surf_usm_h%frac(m,ind_veg_wall)  = 1.0_wp
                 surf_usm_h%frac(m,ind_wat_win)   = 0.0_wp
                 surf_usm_h%frac(m,ind_pav_green) = 0.0_wp
              ENDIF
           ENDIF
        ENDDO
!
!--     If fractions were normalized, give an informative message.
#if defined( __parallel )
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, relative_fractions_corrected, 1,     &
                            MPI_LOGICAL, MPI_LOR, comm2d, ierr )
#endif
        IF ( relative_fractions_corrected )  THEN
           message_string = 'At some horizotal surfaces the relative '   //    &
                            'material fractions do not sum-up to one . ' //    &
                            'Hence, the respective fractions were normalized.'
           CALL message( 'urban_surface_model_mod', 'PA0686', 0, 0, 0, 6, 0 )
        ENDIF
!
!--     Check relative fractions at vertical surfaces.
        relative_fractions_corrected = .FALSE.
        DO  l = 0, 3
           DO  m = 1, surf_usm_v(l)%ns
              sum_frac = SUM( surf_usm_v(l)%frac(m,:) )
              IF ( sum_frac /= 1.0_wp )  THEN
                 relative_fractions_corrected = .TRUE.
!
!--              Normalize relative fractions to 1.
                 IF ( sum_frac /= 0.0_wp )  THEN
                    surf_usm_v(l)%frac(m,:) = surf_usm_v(l)%frac(m,:) / sum_frac
!
!--              In case all relative fractions are erroneously set to zero,
!--              set wall fraction to 1.
                 ELSE
                    surf_usm_v(l)%frac(m,ind_veg_wall)  = 1.0_wp
                    surf_usm_v(l)%frac(m,ind_wat_win)   = 0.0_wp
                    surf_usm_v(l)%frac(m,ind_pav_green) = 0.0_wp
                 ENDIF
              ENDIF
           ENDDO
        ENDDO
!
!--     Also here, ff fractions were normalized, give an informative message.
#if defined( __parallel )
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, relative_fractions_corrected, 1,     &
                            MPI_LOGICAL, MPI_LOR, comm2d, ierr )
#endif
        IF ( relative_fractions_corrected )  THEN
           message_string = 'At some vertical surfaces the relative '    //    &
                            'material fractions do not sum-up to one . ' //    &
                            'Hence, the respective fractions were normalized.'
           CALL message( 'urban_surface_model_mod', 'PA0686', 0, 0, 0, 6, 0 )
        ENDIF
!
!--     Read the surface_types array.
!--     Please note, here also initialization of surface attributes is done as
!--     long as _urbsurf and _surfpar files are available. Values from above
!--     will be overwritten. This might be removed later, but is still in the 
!--     code to enable compatibility with older model version.
        CALL usm_read_urban_surface_types()
        
        CALL usm_init_material_model()
!        
!--     init anthropogenic sources of heat
        IF ( usm_anthropogenic_heat )  THEN
!
!--         init anthropogenic sources of heat (from transportation for now)
            CALL usm_read_anthropogenic_heat()
        ENDIF

!
!--    Check for consistent initialization.
!--    Check if roughness length for momentum, or heat, exceed surface-layer
!--    height and decrease local roughness length where necessary. 
       DO  m = 1, surf_usm_h%ns
          IF ( surf_usm_h%z0(m) >= surf_usm_h%z_mo(m) )  THEN
          
             surf_usm_h%z0(m) = 0.9_wp * surf_usm_h%z_mo(m)
             
             WRITE( message_string, * ) 'z0 exceeds surface-layer height ' //  &
                            'at horizontal urban surface and is ' //           &
                            'decreased appropriately at grid point (i,j) = ',  &
                            surf_usm_h%i(m), surf_usm_h%j(m)
             CALL message( 'urban_surface_model_mod', 'PA0503',                &
                            0, 0, myid, 6, 0 )
          ENDIF
          IF ( surf_usm_h%z0h(m) >= surf_usm_h%z_mo(m) )  THEN
          
             surf_usm_h%z0h(m) = 0.9_wp * surf_usm_h%z_mo(m)
             surf_usm_h%z0q(m) = 0.9_wp * surf_usm_h%z_mo(m)
             
             WRITE( message_string, * ) 'z0h exceeds surface-layer height ' // &
                            'at horizontal urban surface and is ' //           &
                            'decreased appropriately at grid point (i,j) = ',  &
                            surf_usm_h%i(m), surf_usm_h%j(m)
             CALL message( 'urban_surface_model_mod', 'PA0507',                &
                            0, 0, myid, 6, 0 )
          ENDIF          
       ENDDO
       
       DO  l = 0, 3
          DO  m = 1, surf_usm_v(l)%ns
             IF ( surf_usm_v(l)%z0(m) >= surf_usm_v(l)%z_mo(m) )  THEN
          
                surf_usm_v(l)%z0(m) = 0.9_wp * surf_usm_v(l)%z_mo(m)
             
                WRITE( message_string, * ) 'z0 exceeds surface-layer height '// &
                            'at vertical urban surface and is ' //              &
                            'decreased appropriately at grid point (i,j) = ',   &
                            surf_usm_v(l)%i(m)+surf_usm_v(l)%ioff,              &
                            surf_usm_v(l)%j(m)+surf_usm_v(l)%joff
                CALL message( 'urban_surface_model_mod', 'PA0503',              &
                            0, 0, myid, 6, 0 )
             ENDIF
             IF ( surf_usm_v(l)%z0h(m) >= surf_usm_v(l)%z_mo(m) )  THEN
          
                surf_usm_v(l)%z0h(m) = 0.9_wp * surf_usm_v(l)%z_mo(m)
                surf_usm_v(l)%z0q(m) = 0.9_wp * surf_usm_v(l)%z_mo(m)
             
                WRITE( message_string, * ) 'z0h exceeds surface-layer height '// &
                            'at vertical urban surface and is ' //               &
                            'decreased appropriately at grid point (i,j) = ',    &
                            surf_usm_v(l)%i(m)+surf_usm_v(l)%ioff,               &
                            surf_usm_v(l)%j(m)+surf_usm_v(l)%joff
                CALL message( 'urban_surface_model_mod', 'PA0507',               &
                            0, 0, myid, 6, 0 )
             ENDIF
          ENDDO
       ENDDO
!
!--     Intitialization of the surface and wall/ground/roof temperature
!
!--     Initialization for restart runs
        IF ( TRIM( initializing_actions ) /= 'read_restart_data'  .AND.        &
             TRIM( initializing_actions ) /= 'cyclic_fill' )  THEN

!
!--         At horizontal surfaces. Please note, t_surf_wall_h is defined on a 
!--         different data type, but with the same dimension.
            DO  m = 1, surf_usm_h%ns
               i = surf_usm_h%i(m)            
               j = surf_usm_h%j(m)
               k = surf_usm_h%k(m)

               t_surf_wall_h(m) = pt(k,j,i) * exner(k)
               t_surf_window_h(m) = pt(k,j,i) * exner(k)
               t_surf_green_h(m) = pt(k,j,i) * exner(k)
               surf_usm_h%pt_surface(m) = pt(k,j,i) * exner(k)
            ENDDO
!
!--         At vertical surfaces.
            DO  l = 0, 3
               DO  m = 1, surf_usm_v(l)%ns
                  i = surf_usm_v(l)%i(m)            
                  j = surf_usm_v(l)%j(m)
                  k = surf_usm_v(l)%k(m)

                  t_surf_wall_v(l)%t(m) = pt(k,j,i) * exner(k)
                  t_surf_window_v(l)%t(m) = pt(k,j,i) * exner(k)
                  t_surf_green_v(l)%t(m) = pt(k,j,i) * exner(k)
                  surf_usm_v(l)%pt_surface(m) = pt(k,j,i) * exner(k)
               ENDDO
            ENDDO

!
!--         For the sake of correct initialization, set also q_surface. 
!--         Note, at urban surfaces q_surface is initialized with 0.
            IF ( humidity )  THEN
               DO  m = 1, surf_usm_h%ns
                  surf_usm_h%q_surface(m) = 0.0_wp
               ENDDO
               DO  l = 0, 3
                  DO  m = 1, surf_usm_v(l)%ns
                     surf_usm_v(l)%q_surface(m) = 0.0_wp
                  ENDDO
               ENDDO
            ENDIF
!
!--         initial values for t_wall
!--         outer value is set to surface temperature
!--         inner value is set to wall_inner_temperature
!--         and profile is logaritmic (linear in nz).
!--         Horizontal surfaces
            DO  m = 1, surf_usm_h%ns
!
!--            Roof
               IF ( surf_usm_h%isroof_surf(m) )  THEN
                   tin = roof_inner_temperature
                   twin = window_inner_temperature
!
!--            Normal land surface
               ELSE 
                   tin = soil_inner_temperature
                   twin = window_inner_temperature
               ENDIF

               DO k = nzb_wall, nzt_wall+1
                   c = REAL( k - nzb_wall, wp ) /                              &
                       REAL( nzt_wall + 1 - nzb_wall , wp )

                   t_wall_h(k,m) = ( 1.0_wp - c ) * t_surf_wall_h(m) + c * tin
                   t_window_h(k,m) = ( 1.0_wp - c ) * t_surf_window_h(m) + c * twin
                   t_green_h(k,m) = t_surf_wall_h(m)
                   swc_h(k,m) = 0.5_wp
                   swc_sat_h(k,m) = 0.95_wp
                   swc_res_h(k,m) = 0.05_wp
                   rootfr_h(k,m) = 0.1_wp
                   wilt_h(k,m) = 0.1_wp
                   fc_h(k,m) = 0.9_wp
               ENDDO
            ENDDO
!
!--         Vertical surfaces
            DO  l = 0, 3
               DO  m = 1, surf_usm_v(l)%ns
!
!--               Inner wall
                  tin = wall_inner_temperature
                  twin = window_inner_temperature

                  DO k = nzb_wall, nzt_wall+1
                     c = REAL( k - nzb_wall, wp ) /                            &
                         REAL( nzt_wall + 1 - nzb_wall , wp )
                     t_wall_v(l)%t(k,m) = ( 1.0_wp - c ) * t_surf_wall_v(l)%t(m) + c * tin
                     t_window_v(l)%t(k,m) = ( 1.0_wp - c ) * t_surf_window_v(l)%t(m) + c * twin
                     t_green_v(l)%t(k,m) = t_surf_wall_v(l)%t(m)
                  ENDDO
               ENDDO
            ENDDO
        ENDIF

!
!--     If specified, replace constant wall temperatures with fully 3D values from file
        IF ( read_wall_temp_3d )  CALL usm_read_wall_temperature()

!--
!--     Possibly DO user-defined actions (e.g. define heterogeneous wall surface)
        CALL user_init_urban_surface

!
!--     initialize prognostic values for the first timestep
        t_surf_wall_h_p = t_surf_wall_h
        t_surf_wall_v_p = t_surf_wall_v
        t_surf_window_h_p = t_surf_window_h
        t_surf_window_v_p = t_surf_window_v
        t_surf_green_h_p = t_surf_green_h
        t_surf_green_v_p = t_surf_green_v

        t_wall_h_p = t_wall_h
        t_wall_v_p = t_wall_v
        t_window_h_p = t_window_h
        t_window_v_p = t_window_v
        t_green_h_p = t_green_h
        t_green_v_p = t_green_v

!
!--    Set initial values for prognostic soil quantities
       IF ( TRIM( initializing_actions ) /= 'read_restart_data' )  THEN
          m_liq_usm_h%var_usm_1d  = 0.0_wp
       ENDIF
       m_liq_usm_h_p     = m_liq_usm_h
!
!--    Set initial values for prognostic quantities
!--    Horizontal surfaces
       surf_usm_h%c_liq     = 0.0_wp
       surf_usm_h%qsws_liq  = 0.0_wp
       surf_usm_h%qsws_veg  = 0.0_wp

!
!--    Do the same for vertical surfaces
       DO  l = 0, 3
          surf_usm_v(l)%c_liq     = 0.0_wp
          surf_usm_v(l)%qsws_liq  = 0.0_wp
          surf_usm_v(l)%qsws_veg  = 0.0_wp
       ENDDO



        CALL cpu_log( log_point_s(78), 'usm_init', 'stop' )

        IF ( debug_output )  CALL debug_message( 'usm_init', 'end' )

    END SUBROUTINE usm_init


!------------------------------------------------------------------------------!
! Description:
! ------------
!
!> Wall model as part of the urban surface model. The model predicts vertical
!> and horizontal wall / roof temperatures and window layer temperatures.
!> No window layer temperature calculactions during spinup to increase
!> possible timestep.
!------------------------------------------------------------------------------!
    SUBROUTINE usm_material_heat_model( during_spinup )


        IMPLICIT NONE

        INTEGER(iwp) ::  i,j,k,l,kw, m                      !< running indices

        REAL(wp), DIMENSION(nzb_wall:nzt_wall) :: wtend, wintend  !< tendency
        REAL(wp)     :: win_absorp  !< absorption coefficient from transmissivity
        REAL(wp), DIMENSION(nzb_wall:nzt_wall) :: wall_mod

        LOGICAL      :: during_spinup  !< if true, no calculation of window temperatures


        IF ( debug_output_timestep )  THEN
           WRITE( debug_string, * ) 'usm_material_heat_model | during_spinup: ',&
                                     during_spinup
           CALL debug_message( debug_string, 'start' )
        ENDIF

        !$OMP PARALLEL PRIVATE (m, i, j, k, kw, wtend, wintend, win_absorp, wall_mod)
        wall_mod=1.0_wp
        IF ( usm_wall_mod  .AND.  during_spinup )  THEN
           DO  kw=nzb_wall,nzb_wall+1
               wall_mod(kw)=0.1_wp
           ENDDO
        ENDIF

!
!--     For horizontal surfaces                                   
        !$OMP DO SCHEDULE (STATIC)
        DO  m = 1, surf_usm_h%ns
!
!--        Obtain indices
           i = surf_usm_h%i(m)            
           j = surf_usm_h%j(m)
           k = surf_usm_h%k(m)
!
!--        prognostic equation for ground/roof temperature t_wall_h
           wtend(:) = 0.0_wp
           wtend(nzb_wall) = (1.0_wp / surf_usm_h%rho_c_wall(nzb_wall,m)) *        &
                                       ( surf_usm_h%lambda_h(nzb_wall,m) * wall_mod(nzb_wall) *        &
                                         ( t_wall_h(nzb_wall+1,m)                  &
                                         - t_wall_h(nzb_wall,m) ) *                &
                                         surf_usm_h%ddz_wall(nzb_wall+1,m)         &
                                       + surf_usm_h%frac(m,ind_veg_wall)           &
                                         / (surf_usm_h%frac(m,ind_veg_wall)        &
                                           + surf_usm_h%frac(m,ind_pav_green) )    &
                                         * surf_usm_h%wghf_eb(m)                   &
                                       - surf_usm_h%frac(m,ind_pav_green)          &
                                          / (surf_usm_h%frac(m,ind_veg_wall)       &
                                            + surf_usm_h%frac(m,ind_pav_green) )   &
                                         * ( surf_usm_h%lambda_h_green(nzt_wall,m)* wall_mod(nzt_wall) &
                                           * surf_usm_h%ddz_green(nzt_wall,m)      &
                                           + surf_usm_h%lambda_h(nzb_wall,m) * wall_mod(nzb_wall)      &
                                           * surf_usm_h%ddz_wall(nzb_wall,m) )     &
                                         / ( surf_usm_h%ddz_green(nzt_wall,m)      &
                                           + surf_usm_h%ddz_wall(nzb_wall,m) )     &
                                         * ( t_wall_h(nzb_wall,m)                  &
                                           - t_green_h(nzt_wall,m) ) ) *           &
                                       surf_usm_h%ddz_wall_stag(nzb_wall,m)
!
!-- if indoor model is used inner wall layer is calculated by using iwghf (indoor wall ground heat flux)
           IF ( indoor_model ) THEN 
              DO  kw = nzb_wall+1, nzt_wall-1
                  wtend(kw) = (1.0_wp / surf_usm_h%rho_c_wall(kw,m))              &
                                 * (   surf_usm_h%lambda_h(kw,m) * wall_mod(kw)   &
                                    * ( t_wall_h(kw+1,m) - t_wall_h(kw,m) )       &
                                    * surf_usm_h%ddz_wall(kw+1,m)                 &
                                 - surf_usm_h%lambda_h(kw-1,m) * wall_mod(kw-1)   &
                                    * ( t_wall_h(kw,m) - t_wall_h(kw-1,m) )       &
                                    * surf_usm_h%ddz_wall(kw,m)                   &
                                   ) * surf_usm_h%ddz_wall_stag(kw,m)
              ENDDO
              wtend(nzt_wall) = (1.0_wp / surf_usm_h%rho_c_wall(nzt_wall,m)) *    &
                                         ( -surf_usm_h%lambda_h(nzt_wall-1,m) * wall_mod(nzt_wall-1) * &
                                           ( t_wall_h(nzt_wall,m)                 &
                                           - t_wall_h(nzt_wall-1,m) ) *           &
                                           surf_usm_h%ddz_wall(nzt_wall,m)        &
                                         + surf_usm_h%iwghf_eb(m) ) *             &
                                           surf_usm_h%ddz_wall_stag(nzt_wall,m)
           ELSE
              DO  kw = nzb_wall+1, nzt_wall
                  wtend(kw) = (1.0_wp / surf_usm_h%rho_c_wall(kw,m))              &
                                 * (   surf_usm_h%lambda_h(kw,m)  * wall_mod(kw)  &
                                    * ( t_wall_h(kw+1,m) - t_wall_h(kw,m) )       &
                                    * surf_usm_h%ddz_wall(kw+1,m)                 &
                                 - surf_usm_h%lambda_h(kw-1,m) * wall_mod(kw-1)   &
                                    * ( t_wall_h(kw,m) - t_wall_h(kw-1,m) )       &
                                    * surf_usm_h%ddz_wall(kw,m)                   &
                                   ) * surf_usm_h%ddz_wall_stag(kw,m)
              ENDDO
           ENDIF

           t_wall_h_p(nzb_wall:nzt_wall,m) = t_wall_h(nzb_wall:nzt_wall,m)     &
                                 + dt_3d * ( tsc(2)                            &
                                 * wtend(nzb_wall:nzt_wall) + tsc(3)           &
                                 * surf_usm_h%tt_wall_m(nzb_wall:nzt_wall,m) )   

!
!-- during spinup the tempeature inside window layers is not calculated to make larger timesteps possible
           IF ( .NOT. during_spinup ) THEN
              win_absorp = -log(surf_usm_h%transmissivity(m)) / surf_usm_h%zw_window(nzt_wall,m)
!
!--           prognostic equation for ground/roof window temperature t_window_h
!--           takes absorption of shortwave radiation into account
              wintend(:) = 0.0_wp
              wintend(nzb_wall) = (1.0_wp / surf_usm_h%rho_c_window(nzb_wall,m)) *   &
                                         ( surf_usm_h%lambda_h_window(nzb_wall,m) *  &
                                           ( t_window_h(nzb_wall+1,m)                &
                                           - t_window_h(nzb_wall,m) ) *              &
                                           surf_usm_h%ddz_window(nzb_wall+1,m)       &
                                         + surf_usm_h%wghf_eb_window(m)              &
                                         + surf_usm_h%rad_sw_in(m)                   &
                                           * (1.0_wp - exp(-win_absorp               &
                                           * surf_usm_h%zw_window(nzb_wall,m) ) )    &
                                         ) * surf_usm_h%ddz_window_stag(nzb_wall,m)
   
              IF ( indoor_model ) THEN 
                 DO  kw = nzb_wall+1, nzt_wall-1
                     wintend(kw) = (1.0_wp / surf_usm_h%rho_c_window(kw,m))          &
                                    * (   surf_usm_h%lambda_h_window(kw,m)           &
                                       * ( t_window_h(kw+1,m) - t_window_h(kw,m) )   &
                                       * surf_usm_h%ddz_window(kw+1,m)               &
                                    - surf_usm_h%lambda_h_window(kw-1,m)             &
                                       * ( t_window_h(kw,m) - t_window_h(kw-1,m) )   &
                                       * surf_usm_h%ddz_window(kw,m)                 &
                                    + surf_usm_h%rad_sw_in(m)                        &
                                       * (exp(-win_absorp                            &
                                           * surf_usm_h%zw_window(kw-1,m) )          &
                                           - exp(-win_absorp                         &
                                           * surf_usm_h%zw_window(kw,m) ) )          &
                                      ) * surf_usm_h%ddz_window_stag(kw,m)
   
                 ENDDO
                 wintend(nzt_wall) = (1.0_wp / surf_usm_h%rho_c_window(nzt_wall,m)) *       &
                                            ( -surf_usm_h%lambda_h_window(nzt_wall-1,m) *   &
                                              ( t_window_h(nzt_wall,m)                      &
                                              - t_window_h(nzt_wall-1,m) ) *                &
                                              surf_usm_h%ddz_window(nzt_wall,m)             &
                                            + surf_usm_h%iwghf_eb_window(m)                 &
                                            + surf_usm_h%rad_sw_in(m)                       &
                                              * (exp(-win_absorp                            &
                                              * surf_usm_h%zw_window(nzt_wall-1,m) )        &
                                              - exp(-win_absorp                             &
                                              * surf_usm_h%zw_window(nzt_wall,m) ) )        &
                                            ) * surf_usm_h%ddz_window_stag(nzt_wall,m)
              ELSE
                 DO  kw = nzb_wall+1, nzt_wall
                     wintend(kw) = (1.0_wp / surf_usm_h%rho_c_window(kw,m))          &
                                    * (   surf_usm_h%lambda_h_window(kw,m)           &
                                       * ( t_window_h(kw+1,m) - t_window_h(kw,m) )   &
                                       * surf_usm_h%ddz_window(kw+1,m)               &
                                    - surf_usm_h%lambda_h_window(kw-1,m)             &
                                       * ( t_window_h(kw,m) - t_window_h(kw-1,m) )   &
                                       * surf_usm_h%ddz_window(kw,m)                 &
                                    + surf_usm_h%rad_sw_in(m)                        &
                                       * (exp(-win_absorp                            &
                                           * surf_usm_h%zw_window(kw-1,m) )          &
                                           - exp(-win_absorp                         &
                                           * surf_usm_h%zw_window(kw,m) ) )          &
                                      ) * surf_usm_h%ddz_window_stag(kw,m)
   
                 ENDDO
              ENDIF

              t_window_h_p(nzb_wall:nzt_wall,m) = t_window_h(nzb_wall:nzt_wall,m) &
                                 + dt_3d * ( tsc(2)                               &
                                 * wintend(nzb_wall:nzt_wall) + tsc(3)            &
                                 * surf_usm_h%tt_window_m(nzb_wall:nzt_wall,m) )   

           ENDIF

!
!--        calculate t_wall tendencies for the next Runge-Kutta step
           IF ( timestep_scheme(1:5) == 'runge' )  THEN
               IF ( intermediate_timestep_count == 1 )  THEN
                  DO  kw = nzb_wall, nzt_wall
                     surf_usm_h%tt_wall_m(kw,m) = wtend(kw)
                  ENDDO
               ELSEIF ( intermediate_timestep_count <                          &
                        intermediate_timestep_count_max )  THEN
                   DO  kw = nzb_wall, nzt_wall
                      surf_usm_h%tt_wall_m(kw,m) = -9.5625_wp * wtend(kw) +    &
                                         5.3125_wp * surf_usm_h%tt_wall_m(kw,m)
                   ENDDO
               ENDIF
           ENDIF

           IF ( .NOT. during_spinup )  THEN
!
!--           calculate t_window tendencies for the next Runge-Kutta step
              IF ( timestep_scheme(1:5) == 'runge' )  THEN
                  IF ( intermediate_timestep_count == 1 )  THEN
                     DO  kw = nzb_wall, nzt_wall
                        surf_usm_h%tt_window_m(kw,m) = wintend(kw)
                     ENDDO
                  ELSEIF ( intermediate_timestep_count <                            &
                           intermediate_timestep_count_max )  THEN
                      DO  kw = nzb_wall, nzt_wall
                         surf_usm_h%tt_window_m(kw,m) = -9.5625_wp * wintend(kw) +  &
                                            5.3125_wp * surf_usm_h%tt_window_m(kw,m)
                      ENDDO
                  ENDIF
              ENDIF
           ENDIF

        ENDDO

!
!--     For vertical surfaces     
        !$OMP DO SCHEDULE (STATIC)
        DO  l = 0, 3                              
           DO  m = 1, surf_usm_v(l)%ns
!
!--           Obtain indices
              i = surf_usm_v(l)%i(m)            
              j = surf_usm_v(l)%j(m)
              k = surf_usm_v(l)%k(m)
!
!--           prognostic equation for wall temperature t_wall_v
              wtend(:) = 0.0_wp

              wtend(nzb_wall) = (1.0_wp / surf_usm_v(l)%rho_c_wall(nzb_wall,m)) *    &
                                      ( surf_usm_v(l)%lambda_h(nzb_wall,m) * wall_mod(nzb_wall)  *      &
                                        ( t_wall_v(l)%t(nzb_wall+1,m)                &
                                        - t_wall_v(l)%t(nzb_wall,m) ) *              &
                                        surf_usm_v(l)%ddz_wall(nzb_wall+1,m)         &
                                      + surf_usm_v(l)%frac(m,ind_veg_wall)           &
                                        / (surf_usm_v(l)%frac(m,ind_veg_wall)        &
                                          + surf_usm_v(l)%frac(m,ind_pav_green) )    &
                                        * surf_usm_v(l)%wghf_eb(m)                   &
                                      - surf_usm_v(l)%frac(m,ind_pav_green)          &
                                        / (surf_usm_v(l)%frac(m,ind_veg_wall)        &
                                          + surf_usm_v(l)%frac(m,ind_pav_green) )    &
                                        * ( surf_usm_v(l)%lambda_h_green(nzt_wall,m)* wall_mod(nzt_wall) &
                                          * surf_usm_v(l)%ddz_green(nzt_wall,m)      &
                                          + surf_usm_v(l)%lambda_h(nzb_wall,m)* wall_mod(nzb_wall)       &
                                          * surf_usm_v(l)%ddz_wall(nzb_wall,m) )     &
                                        / ( surf_usm_v(l)%ddz_green(nzt_wall,m)      &
                                          + surf_usm_v(l)%ddz_wall(nzb_wall,m) )     &
                                        * ( t_wall_v(l)%t(nzb_wall,m)                &
                                          - t_green_v(l)%t(nzt_wall,m) ) ) *         &
                                        surf_usm_v(l)%ddz_wall_stag(nzb_wall,m)

              IF ( indoor_model ) THEN 
                 DO  kw = nzb_wall+1, nzt_wall-1
                     wtend(kw) = (1.0_wp / surf_usm_v(l)%rho_c_wall(kw,m))        &
                              * (   surf_usm_v(l)%lambda_h(kw,m)  * wall_mod(kw)  &
                                 * ( t_wall_v(l)%t(kw+1,m) - t_wall_v(l)%t(kw,m) )&
                                 * surf_usm_v(l)%ddz_wall(kw+1,m)                 &
                              - surf_usm_v(l)%lambda_h(kw-1,m)  * wall_mod(kw-1)  &
                                 * ( t_wall_v(l)%t(kw,m) - t_wall_v(l)%t(kw-1,m) )&
                                 * surf_usm_v(l)%ddz_wall(kw,m)                   &
                                 ) * surf_usm_v(l)%ddz_wall_stag(kw,m)
                 ENDDO
                 wtend(nzt_wall) = (1.0_wp / surf_usm_v(l)%rho_c_wall(nzt_wall,m)) * &
                                         ( -surf_usm_v(l)%lambda_h(nzt_wall-1,m) * wall_mod(nzt_wall-1)*    &
                                           ( t_wall_v(l)%t(nzt_wall,m)               &
                                           - t_wall_v(l)%t(nzt_wall-1,m) ) *         &
                                           surf_usm_v(l)%ddz_wall(nzt_wall,m)        &
                                         + surf_usm_v(l)%iwghf_eb(m) ) *             &
                                           surf_usm_v(l)%ddz_wall_stag(nzt_wall,m)
              ELSE
                 DO  kw = nzb_wall+1, nzt_wall
                     wtend(kw) = (1.0_wp / surf_usm_v(l)%rho_c_wall(kw,m))        &
                              * (   surf_usm_v(l)%lambda_h(kw,m) * wall_mod(kw)   &
                                 * ( t_wall_v(l)%t(kw+1,m) - t_wall_v(l)%t(kw,m) )&
                                 * surf_usm_v(l)%ddz_wall(kw+1,m)                 &
                              - surf_usm_v(l)%lambda_h(kw-1,m)  * wall_mod(kw-1)  &
                                 * ( t_wall_v(l)%t(kw,m) - t_wall_v(l)%t(kw-1,m) )&
                                 * surf_usm_v(l)%ddz_wall(kw,m)                   &
                                 ) * surf_usm_v(l)%ddz_wall_stag(kw,m)
                 ENDDO
              ENDIF

              t_wall_v_p(l)%t(nzb_wall:nzt_wall,m) =                           &
                                   t_wall_v(l)%t(nzb_wall:nzt_wall,m)          &
                                 + dt_3d * ( tsc(2)                            &
                                 * wtend(nzb_wall:nzt_wall) + tsc(3)           &
                                 * surf_usm_v(l)%tt_wall_m(nzb_wall:nzt_wall,m) )   

              IF ( .NOT. during_spinup )  THEN
                 win_absorp = -log(surf_usm_v(l)%transmissivity(m)) / surf_usm_v(l)%zw_window(nzt_wall,m)
!
!--              prognostic equation for window temperature t_window_v
                 wintend(:) = 0.0_wp
                 wintend(nzb_wall) = (1.0_wp / surf_usm_v(l)%rho_c_window(nzb_wall,m)) * &
                                         ( surf_usm_v(l)%lambda_h_window(nzb_wall,m) *   &
                                           ( t_window_v(l)%t(nzb_wall+1,m)               &
                                           - t_window_v(l)%t(nzb_wall,m) ) *             &
                                           surf_usm_v(l)%ddz_window(nzb_wall+1,m)        &
                                         + surf_usm_v(l)%wghf_eb_window(m)               &
                                         + surf_usm_v(l)%rad_sw_in(m)                    &
                                           * (1.0_wp - exp(-win_absorp                   &
                                           * surf_usm_v(l)%zw_window(nzb_wall,m) ) )     &
                                         ) * surf_usm_v(l)%ddz_window_stag(nzb_wall,m)
   
                 IF ( indoor_model ) THEN 
                    DO  kw = nzb_wall+1, nzt_wall -1
                        wintend(kw) = (1.0_wp / surf_usm_v(l)%rho_c_window(kw,m))         &
                                 * (   surf_usm_v(l)%lambda_h_window(kw,m)                &
                                    * ( t_window_v(l)%t(kw+1,m) - t_window_v(l)%t(kw,m) ) &
                                    * surf_usm_v(l)%ddz_window(kw+1,m)                    &
                                 - surf_usm_v(l)%lambda_h_window(kw-1,m)                  &
                                    * ( t_window_v(l)%t(kw,m) - t_window_v(l)%t(kw-1,m) ) &
                                    * surf_usm_v(l)%ddz_window(kw,m)                      &
                                 + surf_usm_v(l)%rad_sw_in(m)                             &
                                    * (exp(-win_absorp                                    &
                                       * surf_usm_v(l)%zw_window(kw-1,m)       )          &
                                           - exp(-win_absorp                              &
                                           * surf_usm_v(l)%zw_window(kw,m) ) )            &
                                    ) * surf_usm_v(l)%ddz_window_stag(kw,m)
                     ENDDO
                     wintend(nzt_wall) = (1.0_wp / surf_usm_v(l)%rho_c_window(nzt_wall,m)) *  &
                                             ( -surf_usm_v(l)%lambda_h_window(nzt_wall-1,m) * &
                                               ( t_window_v(l)%t(nzt_wall,m)                  &
                                               - t_window_v(l)%t(nzt_wall-1,m) ) *            &
                                               surf_usm_v(l)%ddz_window(nzt_wall,m)           &
                                             + surf_usm_v(l)%iwghf_eb_window(m)               &
                                             + surf_usm_v(l)%rad_sw_in(m)                     &
                                               * (exp(-win_absorp                             &
                                             * surf_usm_v(l)%zw_window(nzt_wall-1,m) )        &
                                           - exp(-win_absorp                                  &
                                               * surf_usm_v(l)%zw_window(nzt_wall,m) ) )      &
                                             ) * surf_usm_v(l)%ddz_window_stag(nzt_wall,m)
                 ELSE
                    DO  kw = nzb_wall+1, nzt_wall
                        wintend(kw) = (1.0_wp / surf_usm_v(l)%rho_c_window(kw,m))         &
                                 * (   surf_usm_v(l)%lambda_h_window(kw,m)                &
                                    * ( t_window_v(l)%t(kw+1,m) - t_window_v(l)%t(kw,m) ) &
                                    * surf_usm_v(l)%ddz_window(kw+1,m)                    &
                                 - surf_usm_v(l)%lambda_h_window(kw-1,m)                  &
                                    * ( t_window_v(l)%t(kw,m) - t_window_v(l)%t(kw-1,m) ) &
                                    * surf_usm_v(l)%ddz_window(kw,m)                      &
                                 + surf_usm_v(l)%rad_sw_in(m)                             &
                                    * (exp(-win_absorp                                    &
                                       * surf_usm_v(l)%zw_window(kw-1,m)       )          &
                                           - exp(-win_absorp                              &
                                           * surf_usm_v(l)%zw_window(kw,m) ) )            &
                                    ) * surf_usm_v(l)%ddz_window_stag(kw,m)
                    ENDDO
                 ENDIF
   
                 t_window_v_p(l)%t(nzb_wall:nzt_wall,m) =                           &
                                      t_window_v(l)%t(nzb_wall:nzt_wall,m)          &
                                    + dt_3d * ( tsc(2)                              &
                                    * wintend(nzb_wall:nzt_wall) + tsc(3)           &
                                    * surf_usm_v(l)%tt_window_m(nzb_wall:nzt_wall,m) )   
              ENDIF

!
!--           calculate t_wall tendencies for the next Runge-Kutta step
              IF ( timestep_scheme(1:5) == 'runge' )  THEN
                  IF ( intermediate_timestep_count == 1 )  THEN
                     DO  kw = nzb_wall, nzt_wall
                        surf_usm_v(l)%tt_wall_m(kw,m) = wtend(kw)
                     ENDDO
                  ELSEIF ( intermediate_timestep_count <                       &
                           intermediate_timestep_count_max )  THEN
                      DO  kw = nzb_wall, nzt_wall
                         surf_usm_v(l)%tt_wall_m(kw,m) =                       &
                                     - 9.5625_wp * wtend(kw) +                 &
                                       5.3125_wp * surf_usm_v(l)%tt_wall_m(kw,m)
                      ENDDO
                  ENDIF
              ENDIF


              IF ( .NOT. during_spinup )  THEN
!
!--              calculate t_window tendencies for the next Runge-Kutta step
                 IF ( timestep_scheme(1:5) == 'runge' )  THEN
                     IF ( intermediate_timestep_count == 1 )  THEN
                        DO  kw = nzb_wall, nzt_wall
                           surf_usm_v(l)%tt_window_m(kw,m) = wintend(kw)
                        ENDDO
                     ELSEIF ( intermediate_timestep_count <                       &
                              intermediate_timestep_count_max )  THEN
                         DO  kw = nzb_wall, nzt_wall
                            surf_usm_v(l)%tt_window_m(kw,m) =                     &
                                        - 9.5625_wp * wintend(kw) +               &
                                          5.3125_wp * surf_usm_v(l)%tt_window_m(kw,m)
                         ENDDO
                     ENDIF
                 ENDIF
              ENDIF

           ENDDO
        ENDDO
        !$OMP END PARALLEL

        IF ( debug_output_timestep )  THEN
           WRITE( debug_string, * ) 'usm_material_heat_model | during_spinup: ',&
                                    during_spinup
           CALL debug_message( debug_string, 'end' )
        ENDIF

    END SUBROUTINE usm_material_heat_model

!------------------------------------------------------------------------------!
! Description:
! ------------
!
!> Green and substrate model as part of the urban surface model. The model predicts ground
!> temperatures.
!>
!> Important: gree-heat model crashes due to unknown reason. Green fraction
!> is thus set to zero (in favor of wall fraction). 
!------------------------------------------------------------------------------!
    SUBROUTINE usm_green_heat_model


        IMPLICIT NONE

        INTEGER(iwp) ::  i,j,k,l,kw, m              !< running indices

        REAL(wp)     :: ke, lambda_h_green_sat      !< heat conductivity for saturated soil
        REAL(wp)     :: h_vg                        !< Van Genuchten coef. h
        REAL(wp)     :: drho_l_lv                   !< frequently used parameter

        REAL(wp), DIMENSION(nzb_wall:nzt_wall) :: gtend,tend  !< tendency

        REAL(wp), DIMENSION(nzb_wall:nzt_wall) :: root_extr_green

        REAL(wp), DIMENSION(nzb_wall:nzt_wall+1) :: lambda_green_temp  !< temp. lambda
        REAL(wp), DIMENSION(nzb_wall:nzt_wall+1) :: gamma_green_temp   !< temp. gamma

        LOGICAL :: conserve_water_content = .true.


        IF ( debug_output_timestep )  CALL debug_message( 'usm_green_heat_model', 'start' )

        drho_l_lv = 1.0_wp / (rho_l * l_v)

!
!--     For horizontal surfaces.
!--     Set tendency array for soil moisture to zero
        IF ( surf_usm_h%ns > 0 )  THEN
           IF ( intermediate_timestep_count == 1 )  surf_usm_h%tswc_h_m = 0.0_wp
        ENDIF
        
        !$OMP PARALLEL PRIVATE (m, i, j, k, kw, lambda_h_green_sat, ke, lambda_green_temp, gtend,  &
        !$OMP&                  tend, h_vg, gamma_green_temp, m_total, root_extr_green)
        !$OMP DO SCHEDULE (STATIC)
        DO  m = 1, surf_usm_h%ns
           IF (surf_usm_h%frac(m,ind_pav_green) > 0.0_wp) THEN
!
!--           Obtain indices
              i = surf_usm_h%i(m)            
              j = surf_usm_h%j(m)
              k = surf_usm_h%k(m)
   
              DO  kw = nzb_wall, nzt_wall
!
!--              Calculate volumetric heat capacity of the soil, taking 
!--              into account water content 
                 surf_usm_h%rho_c_total_green(kw,m) = (surf_usm_h%rho_c_green(kw,m) * (1.0_wp - swc_sat_h(kw,m)) &
                                      + rho_c_water * swc_h(kw,m))
      
!
!--              Calculate soil heat conductivity at the center of the soil
!--              layers
                 lambda_h_green_sat = lambda_h_green_sm ** (1.0_wp - swc_sat_h(kw,m)) *    &
                                lambda_h_water ** swc_h(kw,m)
      
                 ke = 1.0_wp + LOG10(MAX(0.1_wp,swc_h(kw,m)             &
                      / swc_sat_h(kw,m)))
      
                 lambda_green_temp(kw) = ke * (lambda_h_green_sat - lambda_h_green_dry) +    &
                                  lambda_h_green_dry
   
              ENDDO
              lambda_green_temp(nzt_wall+1) = lambda_green_temp(nzt_wall)
   
   
!
!--           Calculate soil heat conductivity (lambda_h) at the _stag level 
!--           using linear interpolation. For pavement surface, the
!--           true pavement depth is considered
              DO  kw = nzb_wall, nzt_wall
                surf_usm_h%lambda_h_green(kw,m) = ( lambda_green_temp(kw+1) + lambda_green_temp(kw) )  &
                                      * 0.5_wp
              ENDDO

              t_green_h(nzt_wall+1,m) = t_wall_h(nzb_wall,m)
!
!--        prognostic equation for ground/roof temperature t_green_h
              gtend(:) = 0.0_wp
              gtend(nzb_wall) = (1.0_wp / surf_usm_h%rho_c_total_green(nzb_wall,m)) *    &
                                         ( surf_usm_h%lambda_h_green(nzb_wall,m) * &
                                           ( t_green_h(nzb_wall+1,m)               &
                                           - t_green_h(nzb_wall,m) ) *             &
                                           surf_usm_h%ddz_green(nzb_wall+1,m)      &
                                         + surf_usm_h%wghf_eb_green(m) ) *         &
                                           surf_usm_h%ddz_green_stag(nzb_wall,m)
              
               DO  kw = nzb_wall+1, nzt_wall
                   gtend(kw) = (1.0_wp / surf_usm_h%rho_c_total_green(kw,m))       &
                                  * (   surf_usm_h%lambda_h_green(kw,m)            &
                                     * ( t_green_h(kw+1,m) - t_green_h(kw,m) )     &
                                     * surf_usm_h%ddz_green(kw+1,m)                &
                                  - surf_usm_h%lambda_h_green(kw-1,m)              &
                                     * ( t_green_h(kw,m) - t_green_h(kw-1,m) )     &
                                     * surf_usm_h%ddz_green(kw,m)                  &
                                    ) * surf_usm_h%ddz_green_stag(kw,m)
               ENDDO
   
              t_green_h_p(nzb_wall:nzt_wall,m) = t_green_h(nzb_wall:nzt_wall,m)    &
                                    + dt_3d * ( tsc(2)                             &
                                    * gtend(nzb_wall:nzt_wall) + tsc(3)            &
                                    * surf_usm_h%tt_green_m(nzb_wall:nzt_wall,m) )   
   
             
!
!--        calculate t_green tendencies for the next Runge-Kutta step
              IF ( timestep_scheme(1:5) == 'runge' )  THEN
                  IF ( intermediate_timestep_count == 1 )  THEN
                     DO  kw = nzb_wall, nzt_wall
                        surf_usm_h%tt_green_m(kw,m) = gtend(kw)
                     ENDDO
                  ELSEIF ( intermediate_timestep_count <                           &
                           intermediate_timestep_count_max )  THEN
                      DO  kw = nzb_wall, nzt_wall
                         surf_usm_h%tt_green_m(kw,m) = -9.5625_wp * gtend(kw) +    &
                                            5.3125_wp * surf_usm_h%tt_green_m(kw,m)
                      ENDDO
                  ENDIF
              ENDIF

              DO  kw = nzb_wall, nzt_wall

!
!--              Calculate soil diffusivity at the center of the soil layers
                 lambda_green_temp(kw) = (- b_ch * surf_usm_h%gamma_w_green_sat(kw,m) * psi_sat       &
                                   / swc_sat_h(kw,m) ) * ( MAX( swc_h(kw,m),    &
                                   wilt_h(kw,m) ) / swc_sat_h(kw,m) )**(        &
                                   b_ch + 2.0_wp )

!
!--              Parametrization of Van Genuchten
                 IF ( soil_type /= 7 )  THEN
!
!--                 Calculate the hydraulic conductivity after Van Genuchten 
!--                 (1980)
                    h_vg = ( ( (swc_res_h(kw,m) - swc_sat_h(kw,m)) / ( swc_res_h(kw,m) -    &
                               MAX( swc_h(kw,m), wilt_h(kw,m) ) ) )**(      &
                               surf_usm_h%n_vg_green(m) / (surf_usm_h%n_vg_green(m) - 1.0_wp ) ) - 1.0_wp  &
                           )**( 1.0_wp / surf_usm_h%n_vg_green(m) ) / surf_usm_h%alpha_vg_green(m)


                    gamma_green_temp(kw) = surf_usm_h%gamma_w_green_sat(kw,m) * ( ( (1.0_wp +         &
                                    ( surf_usm_h%alpha_vg_green(m) * h_vg )**surf_usm_h%n_vg_green(m))**(  &
                                    1.0_wp - 1.0_wp / surf_usm_h%n_vg_green(m) ) - (        &
                                    surf_usm_h%alpha_vg_green(m) * h_vg )**( surf_usm_h%n_vg_green(m)      &
                                    - 1.0_wp) )**2 )                         &
                                    / ( ( 1.0_wp + ( surf_usm_h%alpha_vg_green(m) * h_vg    &
                                    )**surf_usm_h%n_vg_green(m) )**( ( 1.0_wp  - 1.0_wp     &
                                    / surf_usm_h%n_vg_green(m) ) *( surf_usm_h%l_vg_green(m) + 2.0_wp) ) )

!
!--              Parametrization of Clapp & Hornberger
                 ELSE
                    gamma_green_temp(kw) = surf_usm_h%gamma_w_green_sat(kw,m) * ( swc_h(kw,m)       &
                                    / swc_sat_h(kw,m) )**(2.0_wp * b_ch + 3.0_wp)
                 ENDIF

              ENDDO

!
!--           Prognostic equation for soil moisture content. Only performed,
!--           when humidity is enabled in the atmosphere
              IF ( humidity )  THEN
!
!--              Calculate soil diffusivity (lambda_w) at the _stag level 
!--              using linear interpolation. To do: replace this with
!--              ECMWF-IFS Eq. 8.81
                 DO  kw = nzb_wall, nzt_wall-1
                   
                    surf_usm_h%lambda_w_green(kw,m) = ( lambda_green_temp(kw+1) + lambda_green_temp(kw) )  &
                                      * 0.5_wp
                    surf_usm_h%gamma_w_green(kw,m)  = ( gamma_green_temp(kw+1) + gamma_green_temp(kw) )    &
                                      * 0.5_wp

                 ENDDO

!
!--              In case of a closed bottom (= water content is conserved), 
!--              set hydraulic conductivity to zero to that no water will be 
!--              lost in the bottom layer.
                 IF ( conserve_water_content )  THEN
                    surf_usm_h%gamma_w_green(kw,m) = 0.0_wp
                 ELSE
                    surf_usm_h%gamma_w_green(kw,m) = gamma_green_temp(nzt_wall)
                 ENDIF     

!--              The root extraction (= root_extr * qsws_veg / (rho_l     
!--              * l_v)) ensures the mass conservation for water. The         
!--              transpiration of plants equals the cumulative withdrawals by 
!--              the roots in the soil. The scheme takes into account the 
!--              availability of water in the soil layers as well as the root 
!--              fraction in the respective layer. Layer with moisture below 
!--              wilting point will not contribute, which reflects the 
!--              preference of plants to take water from moister layers.

!
!--              Calculate the root extraction (ECMWF 7.69, the sum of 
!--              root_extr = 1). The energy balance solver guarantees a 
!--              positive transpiration, so that there is no need for an 
!--              additional check.
                 m_total = 0.0_wp
                 DO  kw = nzb_wall, nzt_wall
                     IF ( swc_h(kw,m) > wilt_h(kw,m) )  THEN
                        m_total = m_total + rootfr_h(kw,m) * swc_h(kw,m)
                     ENDIF
                 ENDDO  

                 IF ( m_total > 0.0_wp )  THEN
                    DO  kw = nzb_wall, nzt_wall
                       IF ( swc_h(kw,m) > wilt_h(kw,m) )  THEN
                          root_extr_green(kw) = rootfr_h(kw,m) * swc_h(kw,m)      &
                                                          / m_total
                       ELSE
                          root_extr_green(kw) = 0.0_wp
                       ENDIF
                    ENDDO
                 ENDIF

!
!--              Prognostic equation for soil water content m_soil.
                 tend(:) = 0.0_wp

                 tend(nzb_wall) = ( surf_usm_h%lambda_w_green(nzb_wall,m) * (            &
                          swc_h(nzb_wall+1,m) - swc_h(nzb_wall,m) )    &
                          * surf_usm_h%ddz_green(nzb_wall+1,m) - surf_usm_h%gamma_w_green(nzb_wall,m) - ( &
                             root_extr_green(nzb_wall) * surf_usm_h%qsws_veg(m)          &
!                                + surf_usm_h%qsws_soil_green(m)
                                ) * drho_l_lv )             &
                               * surf_usm_h%ddz_green_stag(nzb_wall,m)

                 DO  kw = nzb_wall+1, nzt_wall-1
                    tend(kw) = ( surf_usm_h%lambda_w_green(kw,m) * ( swc_h(kw+1,m)        &
                              - swc_h(kw,m) ) * surf_usm_h%ddz_green(kw+1,m)              &
                              - surf_usm_h%gamma_w_green(kw,m)                            &
                              - surf_usm_h%lambda_w_green(kw-1,m) * (swc_h(kw,m) -        &
                              swc_h(kw-1,m)) * surf_usm_h%ddz_green(kw,m)                 &
                              + surf_usm_h%gamma_w_green(kw-1,m) - (root_extr_green(kw)   &
                              * surf_usm_h%qsws_veg(m) * drho_l_lv)                       &
                              ) * surf_usm_h%ddz_green_stag(kw,m)

                 ENDDO
                 tend(nzt_wall) = ( - surf_usm_h%gamma_w_green(nzt_wall,m)                  &
                                         - surf_usm_h%lambda_w_green(nzt_wall-1,m)          &
                                         * (swc_h(nzt_wall,m)             &
                                         - swc_h(nzt_wall-1,m))           &
                                         * surf_usm_h%ddz_green(nzt_wall,m)                 &
                                         + surf_usm_h%gamma_w_green(nzt_wall-1,m) - (       &
                                           root_extr_green(nzt_wall)               &
                                         * surf_usm_h%qsws_veg(m) * drho_l_lv  )   &
                                   ) * surf_usm_h%ddz_green_stag(nzt_wall,m)             

                 swc_h_p(nzb_wall:nzt_wall,m) = swc_h(nzb_wall:nzt_wall,m)&
                                                 + dt_3d * ( tsc(2) * tend(:)   &
                                                 + tsc(3) * surf_usm_h%tswc_h_m(:,m) )   
 
!
!--              Account for dry soils (find a better solution here!)
                 DO  kw = nzb_wall, nzt_wall
                    IF ( swc_h_p(kw,m) < 0.0_wp )  swc_h_p(kw,m) = 0.0_wp
                 ENDDO

!
!--              Calculate m_soil tendencies for the next Runge-Kutta step
                 IF ( timestep_scheme(1:5) == 'runge' )  THEN
                    IF ( intermediate_timestep_count == 1 )  THEN
                       DO  kw = nzb_wall, nzt_wall
                          surf_usm_h%tswc_h_m(kw,m) = tend(kw)
                       ENDDO
                    ELSEIF ( intermediate_timestep_count <                   &
                             intermediate_timestep_count_max )  THEN
                       DO  kw = nzb_wall, nzt_wall
                          surf_usm_h%tswc_h_m(kw,m) = -9.5625_wp * tend(kw) + 5.3125_wp&
                                   * surf_usm_h%tswc_h_m(kw,m)
                       ENDDO
                    ENDIF
                 ENDIF
              ENDIF

           ENDIF
           
        ENDDO
        !$OMP END PARALLEL

!
!--     For vertical surfaces     
        DO  l = 0, 3                              
           DO  m = 1, surf_usm_v(l)%ns

              IF (surf_usm_v(l)%frac(m,ind_pav_green) > 0.0_wp) THEN
!
!-- no substrate layer for green walls / only groundbase green walls (ivy i.e.) -> green layers get same
!-- temperature as first wall layer
!-- there fore no temperature calculations for vertical green substrate layers now

! 
! !
! !--              Obtain indices
!                  i = surf_usm_v(l)%i(m)            
!                  j = surf_usm_v(l)%j(m)
!                  k = surf_usm_v(l)%k(m)
!    
!                  t_green_v(l)%t(nzt_wall+1,m) = t_wall_v(l)%t(nzb_wall,m)
! !
! !--              prognostic equation for green temperature t_green_v
!                  gtend(:) = 0.0_wp
!                  gtend(nzb_wall) = (1.0_wp / surf_usm_v(l)%rho_c_green(nzb_wall,m)) * &
!                                          ( surf_usm_v(l)%lambda_h_green(nzb_wall,m) * &
!                                            ( t_green_v(l)%t(nzb_wall+1,m)             &
!                                            - t_green_v(l)%t(nzb_wall,m) ) *           &
!                                            surf_usm_v(l)%ddz_green(nzb_wall+1,m)      &
!                                          + surf_usm_v(l)%wghf_eb(m) ) *               &
!                                            surf_usm_v(l)%ddz_green_stag(nzb_wall,m)
!                
!                  DO  kw = nzb_wall+1, nzt_wall
!                     gtend(kw) = (1.0_wp / surf_usm_v(l)%rho_c_green(kw,m))          &
!                               * (   surf_usm_v(l)%lambda_h_green(kw,m)              &
!                                 * ( t_green_v(l)%t(kw+1,m) - t_green_v(l)%t(kw,m) ) &
!                                 * surf_usm_v(l)%ddz_green(kw+1,m)                   &
!                               - surf_usm_v(l)%lambda_h(kw-1,m)                      &
!                                 * ( t_green_v(l)%t(kw,m) - t_green_v(l)%t(kw-1,m) ) &
!                                 * surf_usm_v(l)%ddz_green(kw,m) )                   &
!                               * surf_usm_v(l)%ddz_green_stag(kw,m)
!                  ENDDO
!    
!                  t_green_v_p(l)%t(nzb_wall:nzt_wall,m) =                              &
!                                       t_green_v(l)%t(nzb_wall:nzt_wall,m)             &
!                                     + dt_3d * ( tsc(2)                                &
!                                     * gtend(nzb_wall:nzt_wall) + tsc(3)               &
!                                     * surf_usm_v(l)%tt_green_m(nzb_wall:nzt_wall,m) )   
!    
! !
! !--              calculate t_green tendencies for the next Runge-Kutta step
!                  IF ( timestep_scheme(1:5) == 'runge' )  THEN
!                      IF ( intermediate_timestep_count == 1 )  THEN
!                         DO  kw = nzb_wall, nzt_wall
!                            surf_usm_v(l)%tt_green_m(kw,m) = gtend(kw)
!                         ENDDO
!                      ELSEIF ( intermediate_timestep_count <                           &
!                               intermediate_timestep_count_max )  THEN
!                          DO  kw = nzb_wall, nzt_wall
!                             surf_usm_v(l)%tt_green_m(kw,m) =                          &
!                                         - 9.5625_wp * gtend(kw) +                     &
!                                           5.3125_wp * surf_usm_v(l)%tt_green_m(kw,m)
!                          ENDDO
!                      ENDIF
!                  ENDIF

                 DO  kw = nzb_wall, nzt_wall+1
                     t_green_v(l)%t(kw,m) = t_wall_v(l)%t(nzb_wall,m)
                 ENDDO
              
              ENDIF

           ENDDO
        ENDDO

        IF ( debug_output_timestep )  CALL debug_message( 'usm_green_heat_model', 'end' )

    END SUBROUTINE usm_green_heat_model

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Parin for &usm_par for urban surface model
!------------------------------------------------------------------------------!
    SUBROUTINE usm_parin

       IMPLICIT NONE

       CHARACTER (LEN=80) ::  line  !< string containing current line of file PARIN

       NAMELIST /urban_surface_par/                                            &
                           building_type,                                      &
                           land_category,                                      &
                           naheatlayers,                                       &
                           pedestrian_category,                                &
                           roughness_concrete,                                 &
                           read_wall_temp_3d,                                  &
                           roof_category,                                      &
                           urban_surface,                                      &
                           usm_anthropogenic_heat,                             &
                           usm_material_model,                                 &
                           wall_category,                                      &
                           wall_inner_temperature,                             &
                           roof_inner_temperature,                             &
                           soil_inner_temperature,                             &
                           window_inner_temperature,                           &
                           usm_wall_mod

       NAMELIST /urban_surface_parameters/                                     &
                           building_type,                                      &
                           land_category,                                      &
                           naheatlayers,                                       &
                           pedestrian_category,                                &
                           roughness_concrete,                                 &
                           read_wall_temp_3d,                                  &
                           roof_category,                                      &
                           urban_surface,                                      &
                           usm_anthropogenic_heat,                             &
                           usm_material_model,                                 &
                           wall_category,                                      &
                           wall_inner_temperature,                             &
                           roof_inner_temperature,                             &
                           soil_inner_temperature,                             &
                           window_inner_temperature,                           &
                           usm_wall_mod
                           
  
!
!--    Try to find urban surface model package
       REWIND ( 11 )
       line = ' '
       DO WHILE ( INDEX( line, '&urban_surface_parameters' ) == 0 )
          READ ( 11, '(A)', END=12 )  line
       ENDDO
       BACKSPACE ( 11 )

!
!--    Read user-defined namelist
       READ ( 11, urban_surface_parameters, ERR = 10 )

!
!--    Set flag that indicates that the urban surface model is switched on
       urban_surface = .TRUE.

       GOTO 14

 10    BACKSPACE( 11 )
       READ( 11 , '(A)') line
       CALL parin_fail_message( 'urban_surface_parameters', line )
!
!--    Try to find old namelist
 12    REWIND ( 11 )
       line = ' '
       DO WHILE ( INDEX( line, '&urban_surface_par' ) == 0 )
          READ ( 11, '(A)', END=14 )  line
       ENDDO
       BACKSPACE ( 11 )

!
!--    Read user-defined namelist
       READ ( 11, urban_surface_par, ERR = 13, END = 14 )

       message_string = 'namelist urban_surface_par is deprecated and will be ' // &
                     'removed in near future. Please use namelist ' //   &
                     'urban_surface_parameters instead'
       CALL message( 'usm_parin', 'PA0487', 0, 1, 0, 6, 0 )

!
!--    Set flag that indicates that the urban surface model is switched on
       urban_surface = .TRUE.

       GOTO 14

 13    BACKSPACE( 11 )
       READ( 11 , '(A)') line
       CALL parin_fail_message( 'urban_surface_par', line )


 14    CONTINUE


    END SUBROUTINE usm_parin

 
!------------------------------------------------------------------------------!
! Description:
! ------------
!
!> This subroutine is part of the urban surface model.
!> It reads daily heat produced by anthropogenic sources
!> and the diurnal cycle of the heat.
!------------------------------------------------------------------------------!
    SUBROUTINE usm_read_anthropogenic_heat
    
        INTEGER(iwp)                  :: i,j,k,ii  !< running indices
        REAL(wp)                      :: heat      !< anthropogenic heat

!
!--     allocation of array of sources of anthropogenic heat and their diural profile
        ALLOCATE( aheat(naheatlayers,nys:nyn,nxl:nxr) )
        ALLOCATE( aheatprof(naheatlayers,0:24) )

!
!--     read daily amount of heat and its daily cycle
        aheat = 0.0_wp
        DO  ii = 0, io_blocks-1
            IF ( ii == io_group )  THEN

!--             open anthropogenic heat file 
                OPEN( 151, file='ANTHROPOGENIC_HEAT'//TRIM(coupling_char), action='read', &
                           status='old', form='formatted', err=11 )
                i = 0
                j = 0
                DO
                    READ( 151, *, err=12, end=13 )  i, j, k, heat
                    IF ( i >= nxl  .AND.  i <= nxr  .AND.  j >= nys  .AND.  j <= nyn )  THEN
                        IF ( k <= naheatlayers  .AND.  k > topo_top_ind(j,i,0) )  THEN
!--                         write heat into the array
                            aheat(k,j,i) = heat
                        ENDIF
                    ENDIF
                    CYCLE
 12                 WRITE(message_string,'(a,2i4)') 'error in file ANTHROPOGENIC_HEAT'//TRIM(coupling_char)//' after line ',i,j
                    CALL message( 'usm_read_anthropogenic_heat', 'PA0515', 0, 1, 0, 6, 0 )
                ENDDO
 13             CLOSE(151)
                CYCLE
 11             message_string = 'file ANTHROPOGENIC_HEAT'//TRIM(coupling_char)//' does not exist'
                CALL message( 'usm_read_anthropogenic_heat', 'PA0516', 1, 2, 0, 6, 0 )
            ENDIF
            
#if defined( __parallel )
            CALL MPI_BARRIER( comm2d, ierr )
#endif
        ENDDO
        
!
!--     read diurnal profiles of heat sources
        aheatprof = 0.0_wp
        DO  ii = 0, io_blocks-1
            IF ( ii == io_group )  THEN
!
!--             open anthropogenic heat profile file 
                OPEN( 151, file='ANTHROPOGENIC_HEAT_PROFILE'//TRIM(coupling_char), action='read', &
                           status='old', form='formatted', err=21 )
                i = 0
                DO
                    READ( 151, *, err=22, end=23 )  i, k, heat
                    IF ( i >= 0  .AND.  i <= 24  .AND.  k <= naheatlayers )  THEN
!--                     write heat into the array
                        aheatprof(k,i) = heat
                    ENDIF
                    CYCLE
 22                 WRITE(message_string,'(a,i4)') 'error in file ANTHROPOGENIC_HEAT_PROFILE'// &
                                                     TRIM(coupling_char)//' after line ',i
                    CALL message( 'usm_read_anthropogenic_heat', 'PA0517', 0, 1, 0, 6, 0 )
                ENDDO
                aheatprof(:,24) = aheatprof(:,0)
 23             CLOSE(151)
                CYCLE
 21             message_string = 'file ANTHROPOGENIC_HEAT_PROFILE'//TRIM(coupling_char)//' does not exist'
                CALL message( 'usm_read_anthropogenic_heat', 'PA0518', 1, 2, 0, 6, 0 )
            ENDIF
            
#if defined( __parallel )
            CALL MPI_BARRIER( comm2d, ierr )
#endif
        ENDDO
        
    END SUBROUTINE usm_read_anthropogenic_heat
   

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Soubroutine reads t_surf and t_wall data from restart files
!------------------------------------------------------------------------------!
    SUBROUTINE usm_rrd_local( k, nxlf, nxlc, nxl_on_file, nxrf, nxr_on_file, nynf, nyn_on_file,    &
                              nysf, nysc, nys_on_file, found )


       USE control_parameters,                                                 &
           ONLY: length, restart_string
           
       IMPLICIT NONE

       INTEGER(iwp)       ::  k                 !< running index over previous input files covering current local domain
       INTEGER(iwp)       ::  l                 !< index variable for surface type
       INTEGER(iwp)       ::  ns_h_on_file_usm  !< number of horizontal surface elements (urban type) on file
       INTEGER(iwp)       ::  nxlc              !< index of left boundary on current subdomain
       INTEGER(iwp)       ::  nxlf              !< index of left boundary on former subdomain 
       INTEGER(iwp)       ::  nxl_on_file       !< index of left boundary on former local domain 
       INTEGER(iwp)       ::  nxrf              !< index of right boundary on former subdomain
       INTEGER(iwp)       ::  nxr_on_file       !< index of right boundary on former local domain 
       INTEGER(iwp)       ::  nynf              !< index of north boundary on former subdomain
       INTEGER(iwp)       ::  nyn_on_file       !< index of north boundary on former local domain 
       INTEGER(iwp)       ::  nysc              !< index of south boundary on current subdomain 
       INTEGER(iwp)       ::  nysf              !< index of south boundary on former subdomain
       INTEGER(iwp)       ::  nys_on_file       !< index of south boundary on former local domain
       
       INTEGER(iwp)       ::  ns_v_on_file_usm(0:3)  !< number of vertical surface elements (urban type) on file
!
!--    Note, the save attribute in the following array declaration is necessary, 
!--    in order to keep the number of urban surface elements on file during
!--    rrd_local calls. 
       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE, SAVE ::  start_index_on_file 
       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE, SAVE ::  end_index_on_file

       LOGICAL, INTENT(OUT)  ::  found 
! MS: Why are there individual temporary arrays that all have the same size?
       REAL(wp), DIMENSION(:), ALLOCATABLE, SAVE ::  tmp_surf_wall_h
       REAL(wp), DIMENSION(:), ALLOCATABLE, SAVE ::  tmp_surf_window_h
       REAL(wp), DIMENSION(:), ALLOCATABLE, SAVE ::  tmp_surf_green_h
       REAL(wp), DIMENSION(:), ALLOCATABLE, SAVE ::  tmp_surf_mliq_h
       REAL(wp), DIMENSION(:), ALLOCATABLE, SAVE ::  tmp_surf_waste_h
       
       REAL(wp), DIMENSION(:,:), ALLOCATABLE, SAVE ::  tmp_wall_h
       REAL(wp), DIMENSION(:,:), ALLOCATABLE, SAVE ::  tmp_window_h
       REAL(wp), DIMENSION(:,:), ALLOCATABLE, SAVE ::  tmp_green_h
       
       TYPE( t_surf_vertical ), DIMENSION(0:3), SAVE ::  tmp_surf_wall_v
       TYPE( t_surf_vertical ), DIMENSION(0:3), SAVE ::  tmp_surf_window_v
       TYPE( t_surf_vertical ), DIMENSION(0:3), SAVE ::  tmp_surf_green_v
       TYPE( t_surf_vertical ), DIMENSION(0:3), SAVE ::  tmp_surf_waste_v
       
       TYPE( t_wall_vertical ), DIMENSION(0:3), SAVE ::  tmp_wall_v
       TYPE( t_wall_vertical ), DIMENSION(0:3), SAVE ::  tmp_window_v
       TYPE( t_wall_vertical ), DIMENSION(0:3), SAVE ::  tmp_green_v


       found = .TRUE.


       SELECT CASE ( restart_string(1:length) ) 

          CASE ( 'ns_h_on_file_usm') 
             IF ( k == 1 )  THEN
                READ ( 13 ) ns_h_on_file_usm
             
                IF ( ALLOCATED( tmp_surf_wall_h ) ) DEALLOCATE( tmp_surf_wall_h )
                IF ( ALLOCATED( tmp_wall_h ) ) DEALLOCATE( tmp_wall_h ) 
                IF ( ALLOCATED( tmp_surf_window_h ) )                       &
                   DEALLOCATE( tmp_surf_window_h ) 
                IF ( ALLOCATED( tmp_window_h) ) DEALLOCATE( tmp_window_h ) 
                IF ( ALLOCATED( tmp_surf_green_h) )                         &
                   DEALLOCATE( tmp_surf_green_h ) 
                IF ( ALLOCATED( tmp_green_h) ) DEALLOCATE( tmp_green_h )
                IF ( ALLOCATED( tmp_surf_mliq_h) )                         &
                   DEALLOCATE( tmp_surf_mliq_h )
                IF ( ALLOCATED( tmp_surf_waste_h) )                         &
                   DEALLOCATE( tmp_surf_waste_h )
  
!
!--             Allocate temporary arrays for reading data on file. Note,
!--             the size of allocated surface elements do not necessarily
!--             need  to match the size of present surface elements on 
!--             current processor, as the number of processors between 
!--             restarts can change. 
                ALLOCATE( tmp_surf_wall_h(1:ns_h_on_file_usm) )
                ALLOCATE( tmp_wall_h(nzb_wall:nzt_wall+1,                   &
                                     1:ns_h_on_file_usm) )
                ALLOCATE( tmp_surf_window_h(1:ns_h_on_file_usm) )
                ALLOCATE( tmp_window_h(nzb_wall:nzt_wall+1,                 &
                                       1:ns_h_on_file_usm) )
                ALLOCATE( tmp_surf_green_h(1:ns_h_on_file_usm) )
                ALLOCATE( tmp_green_h(nzb_wall:nzt_wall+1,                  &
                                      1:ns_h_on_file_usm) )
                ALLOCATE( tmp_surf_mliq_h(1:ns_h_on_file_usm) )
                ALLOCATE( tmp_surf_waste_h(1:ns_h_on_file_usm) )

             ENDIF

          CASE ( 'ns_v_on_file_usm')
             IF ( k == 1 )  THEN
                READ ( 13 ) ns_v_on_file_usm 

                DO  l = 0, 3
                   IF ( ALLOCATED( tmp_surf_wall_v(l)%t ) )                 &
                      DEALLOCATE( tmp_surf_wall_v(l)%t )
                   IF ( ALLOCATED( tmp_wall_v(l)%t ) )                      &
                      DEALLOCATE( tmp_wall_v(l)%t )
                   IF ( ALLOCATED( tmp_surf_window_v(l)%t ) )               & 
                      DEALLOCATE( tmp_surf_window_v(l)%t )
                   IF ( ALLOCATED( tmp_window_v(l)%t ) )                    &
                      DEALLOCATE( tmp_window_v(l)%t )
                   IF ( ALLOCATED( tmp_surf_green_v(l)%t ) )                &
                      DEALLOCATE( tmp_surf_green_v(l)%t )
                   IF ( ALLOCATED( tmp_green_v(l)%t ) )                     &
                      DEALLOCATE( tmp_green_v(l)%t )
                   IF ( ALLOCATED( tmp_surf_waste_v(l)%t ) )                &
                      DEALLOCATE( tmp_surf_waste_v(l)%t )
                ENDDO 

!
!--             Allocate temporary arrays for reading data on file. Note,
!--             the size of allocated surface elements do not necessarily
!--             need to match the size of present surface elements on 
!--             current processor, as the number of processors between 
!--             restarts can change. 
                DO  l = 0, 3
                   ALLOCATE( tmp_surf_wall_v(l)%t(1:ns_v_on_file_usm(l)) )
                   ALLOCATE( tmp_wall_v(l)%t(nzb_wall:nzt_wall+1,           &
                                             1:ns_v_on_file_usm(l) ) )
                   ALLOCATE( tmp_surf_window_v(l)%t(1:ns_v_on_file_usm(l)) )
                   ALLOCATE( tmp_window_v(l)%t(nzb_wall:nzt_wall+1,         & 
                                               1:ns_v_on_file_usm(l) ) )
                   ALLOCATE( tmp_surf_green_v(l)%t(1:ns_v_on_file_usm(l)) )
                   ALLOCATE( tmp_green_v(l)%t(nzb_wall:nzt_wall+1,          &
                                              1:ns_v_on_file_usm(l) ) )
                   ALLOCATE( tmp_surf_waste_v(l)%t(1:ns_v_on_file_usm(l)) )
                ENDDO

             ENDIF    
       
          CASE ( 'usm_start_index_h', 'usm_start_index_v'  )   
             IF ( k == 1 )  THEN

                IF ( ALLOCATED( start_index_on_file ) )                     &
                   DEALLOCATE( start_index_on_file )

                ALLOCATE ( start_index_on_file(nys_on_file:nyn_on_file,     &
                                               nxl_on_file:nxr_on_file) )

                READ ( 13 )  start_index_on_file

             ENDIF
             
          CASE ( 'usm_end_index_h', 'usm_end_index_v' )   
             IF ( k == 1 )  THEN

                IF ( ALLOCATED( end_index_on_file ) )                       &
                   DEALLOCATE( end_index_on_file )

                ALLOCATE ( end_index_on_file(nys_on_file:nyn_on_file,       &
                                             nxl_on_file:nxr_on_file) )

                READ ( 13 )  end_index_on_file

             ENDIF
       
          CASE ( 't_surf_wall_h' )
             IF ( k == 1 )  THEN
                IF ( .NOT.  ALLOCATED( t_surf_wall_h_1 ) )                  &
                   ALLOCATE( t_surf_wall_h_1(1:surf_usm_h%ns) )
                READ ( 13 )  tmp_surf_wall_h
             ENDIF              
             CALL surface_restore_elements(                                 &
                                     t_surf_wall_h_1, tmp_surf_wall_h,      &
                                     surf_usm_h%start_index,                &
                                     start_index_on_file,                   &
                                     end_index_on_file,                     &
                                     nxlc, nysc,                            &
                                     nxlf, nxrf, nysf, nynf,                &
                                     nys_on_file, nyn_on_file,              &
                                     nxl_on_file,nxr_on_file )

          CASE ( 't_surf_wall_v(0)' )
             IF ( k == 1 )  THEN
                IF ( .NOT.  ALLOCATED( t_surf_wall_v_1(0)%t ) )             &
                   ALLOCATE( t_surf_wall_v_1(0)%t(1:surf_usm_v(0)%ns) )
                READ ( 13 )  tmp_surf_wall_v(0)%t
             ENDIF
             CALL surface_restore_elements(                                 &
                                     t_surf_wall_v_1(0)%t, tmp_surf_wall_v(0)%t,      &
                                     surf_usm_v(0)%start_index,             & 
                                     start_index_on_file,                   &
                                     end_index_on_file,                     &
                                     nxlc, nysc,                            &
                                     nxlf, nxrf, nysf, nynf,                &
                                     nys_on_file, nyn_on_file,              &
                                     nxl_on_file,nxr_on_file )
                   
          CASE ( 't_surf_wall_v(1)' )
             IF ( k == 1 )  THEN
                IF ( .NOT.  ALLOCATED( t_surf_wall_v_1(1)%t ) )             &
                   ALLOCATE( t_surf_wall_v_1(1)%t(1:surf_usm_v(1)%ns) )
                READ ( 13 )  tmp_surf_wall_v(1)%t
             ENDIF
             CALL surface_restore_elements(                                 &
                                     t_surf_wall_v_1(1)%t, tmp_surf_wall_v(1)%t,      &
                                     surf_usm_v(1)%start_index,             & 
                                     start_index_on_file,                   &
                                     end_index_on_file,                     &
                                     nxlc, nysc,                            &
                                     nxlf, nxrf, nysf, nynf,                &
                                     nys_on_file, nyn_on_file,              &
                                     nxl_on_file,nxr_on_file )

          CASE ( 't_surf_wall_v(2)' )
             IF ( k == 1 )  THEN
                IF ( .NOT.  ALLOCATED( t_surf_wall_v_1(2)%t ) )             &
                   ALLOCATE( t_surf_wall_v_1(2)%t(1:surf_usm_v(2)%ns) )
                READ ( 13 )  tmp_surf_wall_v(2)%t
             ENDIF
             CALL surface_restore_elements(                                 &
                                     t_surf_wall_v_1(2)%t, tmp_surf_wall_v(2)%t,      &
                                     surf_usm_v(2)%start_index,             &  
                                     start_index_on_file,                   &
                                     end_index_on_file,                     &
                                     nxlc, nysc,                            &
                                     nxlf, nxrf, nysf, nynf,                &
                                     nys_on_file, nyn_on_file,              &
                                     nxl_on_file,nxr_on_file )
                   
          CASE ( 't_surf_wall_v(3)' )
             IF ( k == 1 )  THEN
                IF ( .NOT.  ALLOCATED( t_surf_wall_v_1(3)%t ) )             &
                   ALLOCATE( t_surf_wall_v_1(3)%t(1:surf_usm_v(3)%ns) )
                READ ( 13 )  tmp_surf_wall_v(3)%t
             ENDIF
             CALL surface_restore_elements(                                 &
                                     t_surf_wall_v_1(3)%t, tmp_surf_wall_v(3)%t,      &
                                     surf_usm_v(3)%start_index,             & 
                                     start_index_on_file,                   &
                                     end_index_on_file,                     &
                                     nxlc, nysc,                            &
                                     nxlf, nxrf, nysf, nynf,                &
                                     nys_on_file, nyn_on_file,              &
                                     nxl_on_file,nxr_on_file )

          CASE ( 't_surf_green_h' )
             IF ( k == 1 )  THEN
                IF ( .NOT.  ALLOCATED( t_surf_green_h_1 ) )                 &
                   ALLOCATE( t_surf_green_h_1(1:surf_usm_h%ns) )
                READ ( 13 )  tmp_surf_green_h
             ENDIF
             CALL surface_restore_elements(                                 &
                                     t_surf_green_h_1, tmp_surf_green_h,    &
                                     surf_usm_h%start_index,                &  
                                     start_index_on_file,                   &
                                     end_index_on_file,                     &
                                     nxlc, nysc,                            &
                                     nxlf, nxrf, nysf, nynf,                &
                                     nys_on_file, nyn_on_file,              &
                                     nxl_on_file,nxr_on_file )

          CASE ( 't_surf_green_v(0)' )
             IF ( k == 1 )  THEN
                IF ( .NOT.  ALLOCATED( t_surf_green_v_1(0)%t ) )            &
                   ALLOCATE( t_surf_green_v_1(0)%t(1:surf_usm_v(0)%ns) )
                READ ( 13 )  tmp_surf_green_v(0)%t
             ENDIF
             CALL surface_restore_elements(                                 &
                                     t_surf_green_v_1(0)%t,                 &
                                     tmp_surf_green_v(0)%t,                 &
                                     surf_usm_v(0)%start_index,             & 
                                     start_index_on_file,                   &
                                     end_index_on_file,                     &
                                     nxlc, nysc,                            &
                                     nxlf, nxrf, nysf, nynf,                &
                                     nys_on_file, nyn_on_file,              &
                                     nxl_on_file,nxr_on_file )
                
          CASE ( 't_surf_green_v(1)' )
             IF ( k == 1 )  THEN
                IF ( .NOT.  ALLOCATED( t_surf_green_v_1(1)%t ) )            &
                   ALLOCATE( t_surf_green_v_1(1)%t(1:surf_usm_v(1)%ns) )
                READ ( 13 )  tmp_surf_green_v(1)%t
             ENDIF
             CALL surface_restore_elements(                                 &
                                     t_surf_green_v_1(1)%t,                 &
                                     tmp_surf_green_v(1)%t,                 &
                                     surf_usm_v(1)%start_index,             &  
                                     start_index_on_file,                   &
                                     end_index_on_file,                     &
                                     nxlc, nysc,                            &
                                     nxlf, nxrf, nysf, nynf,                &
                                     nys_on_file, nyn_on_file,              &
                                     nxl_on_file,nxr_on_file )

          CASE ( 't_surf_green_v(2)' )
             IF ( k == 1 )  THEN
                IF ( .NOT.  ALLOCATED( t_surf_green_v_1(2)%t ) )            &
                   ALLOCATE( t_surf_green_v_1(2)%t(1:surf_usm_v(2)%ns) )
                READ ( 13 )  tmp_surf_green_v(2)%t
             ENDIF
             CALL surface_restore_elements(                                 &
                                     t_surf_green_v_1(2)%t,                 &
                                     tmp_surf_green_v(2)%t,                 &
                                     surf_usm_v(2)%start_index,             &  
                                     start_index_on_file,                   &
                                     end_index_on_file,                     &
                                     nxlc, nysc,                            &
                                     nxlf, nxrf, nysf, nynf,                &
                                     nys_on_file, nyn_on_file,              &
                                     nxl_on_file,nxr_on_file )
                
          CASE ( 't_surf_green_v(3)' )
             IF ( k == 1 )  THEN
                IF ( .NOT.  ALLOCATED( t_surf_green_v_1(3)%t ) )            &
                   ALLOCATE( t_surf_green_v_1(3)%t(1:surf_usm_v(3)%ns) )
                READ ( 13 )  tmp_surf_green_v(3)%t
             ENDIF
             CALL surface_restore_elements(                                 &
                                     t_surf_green_v_1(3)%t,                 & 
                                     tmp_surf_green_v(3)%t,                 &
                                     surf_usm_v(3)%start_index,             & 
                                     start_index_on_file,                   &
                                     end_index_on_file,                     &
                                     nxlc, nysc,                            &
                                     nxlf, nxrf, nysf, nynf,                &
                                     nys_on_file, nyn_on_file,              &
                                     nxl_on_file,nxr_on_file )

          CASE ( 't_surf_window_h' )
             IF ( k == 1 )  THEN
                IF ( .NOT.  ALLOCATED( t_surf_window_h_1 ) )                &
                   ALLOCATE( t_surf_window_h_1(1:surf_usm_h%ns) )
                READ ( 13 )  tmp_surf_window_h
             ENDIF
             CALL surface_restore_elements(                                 &
                                     t_surf_window_h_1,                     &
                                     tmp_surf_window_h,                     &
                                     surf_usm_h%start_index,                & 
                                     start_index_on_file,                   &
                                     end_index_on_file,                     &
                                     nxlc, nysc,                            &
                                     nxlf, nxrf, nysf, nynf,                &
                                     nys_on_file, nyn_on_file,              &
                                     nxl_on_file,nxr_on_file )

          CASE ( 't_surf_window_v(0)' )
             IF ( k == 1 )  THEN
                IF ( .NOT.  ALLOCATED( t_surf_window_v_1(0)%t ) )           &
                   ALLOCATE( t_surf_window_v_1(0)%t(1:surf_usm_v(0)%ns) )
                READ ( 13 )  tmp_surf_window_v(0)%t
             ENDIF
             CALL surface_restore_elements(                                 &
                                     t_surf_window_v_1(0)%t,                &
                                     tmp_surf_window_v(0)%t,                &
                                     surf_usm_v(0)%start_index,             & 
                                     start_index_on_file,                   &
                                     end_index_on_file,                     &
                                     nxlc, nysc,                            &
                                     nxlf, nxrf, nysf, nynf,                &
                                     nys_on_file, nyn_on_file,              &
                                     nxl_on_file,nxr_on_file )
                
          CASE ( 't_surf_window_v(1)' )
             IF ( k == 1 )  THEN
                IF ( .NOT.  ALLOCATED( t_surf_window_v_1(1)%t ) )           &
                   ALLOCATE( t_surf_window_v_1(1)%t(1:surf_usm_v(1)%ns) )
                READ ( 13 )  tmp_surf_window_v(1)%t
             ENDIF
             CALL surface_restore_elements(                                 &
                                     t_surf_window_v_1(1)%t,                &
                                     tmp_surf_window_v(1)%t,                &
                                     surf_usm_v(1)%start_index,             & 
                                     start_index_on_file,                   &
                                     end_index_on_file,                     &
                                     nxlc, nysc,                            &
                                     nxlf, nxrf, nysf, nynf,                &
                                     nys_on_file, nyn_on_file,              &
                                     nxl_on_file,nxr_on_file )

          CASE ( 't_surf_window_v(2)' )
             IF ( k == 1 )  THEN
                IF ( .NOT.  ALLOCATED( t_surf_window_v_1(2)%t ) )           &
                   ALLOCATE( t_surf_window_v_1(2)%t(1:surf_usm_v(2)%ns) )
                READ ( 13 )  tmp_surf_window_v(2)%t
             ENDIF
             CALL surface_restore_elements(                                 &
                                     t_surf_window_v_1(2)%t,                & 
                                     tmp_surf_window_v(2)%t,                &
                                     surf_usm_v(2)%start_index,             & 
                                     start_index_on_file,                   &
                                     end_index_on_file,                     &
                                     nxlc, nysc,                            &
                                     nxlf, nxrf, nysf, nynf,                &
                                     nys_on_file, nyn_on_file,              &
                                     nxl_on_file,nxr_on_file )
                
          CASE ( 't_surf_window_v(3)' )
             IF ( k == 1 )  THEN
                IF ( .NOT.  ALLOCATED( t_surf_window_v_1(3)%t ) )           &
                   ALLOCATE( t_surf_window_v_1(3)%t(1:surf_usm_v(3)%ns) )
                READ ( 13 )  tmp_surf_window_v(3)%t
             ENDIF
             CALL surface_restore_elements(                                 &
                                     t_surf_window_v_1(3)%t,                &  
                                     tmp_surf_window_v(3)%t,                &
                                     surf_usm_v(3)%start_index,             & 
                                     start_index_on_file,                   &
                                     end_index_on_file,                     &
                                     nxlc, nysc,                            &
                                     nxlf, nxrf, nysf, nynf,                &
                                     nys_on_file, nyn_on_file,              &
                                     nxl_on_file,nxr_on_file )

          CASE ( 'm_liq_usm_h' )
             IF ( k == 1 )  THEN
                IF ( .NOT.  ALLOCATED( m_liq_usm_h%var_usm_1d ) )           &
                   ALLOCATE( m_liq_usm_h%var_usm_1d(1:surf_usm_h%ns) )
                READ ( 13 )  tmp_surf_mliq_h
             ENDIF              
             CALL surface_restore_elements(                                 &
                                     m_liq_usm_h%var_usm_1d,                &
                                     tmp_surf_mliq_h,                       &
                                     surf_usm_h%start_index,                &
                                     start_index_on_file,                   &
                                     end_index_on_file,                     &
                                     nxlc, nysc,                            &
                                     nxlf, nxrf, nysf, nynf,                &
                                     nys_on_file, nyn_on_file,              &
                                     nxl_on_file,nxr_on_file )

          CASE ( 'waste_heat_h' )
             IF ( k == 1 )  THEN
                IF ( .NOT.  ALLOCATED( surf_usm_h%waste_heat ) )            &
                   ALLOCATE( surf_usm_h%waste_heat(1:surf_usm_h%ns) )
                READ ( 13 )  tmp_surf_waste_h
             ENDIF              
             CALL surface_restore_elements(                                 &
                                     surf_usm_h%waste_heat,                 &
                                     tmp_surf_waste_h,                      &
                                     surf_usm_h%start_index,                &
                                     start_index_on_file,                   &
                                     end_index_on_file,                     &
                                     nxlc, nysc,                            &
                                     nxlf, nxrf, nysf, nynf,                &
                                     nys_on_file, nyn_on_file,              &
                                     nxl_on_file,nxr_on_file )

          CASE ( 'waste_heat_v(0)' )
             IF ( k == 1 )  THEN
                IF ( .NOT.  ALLOCATED( surf_usm_v(0)%waste_heat ) )         &
                   ALLOCATE( surf_usm_v(0)%waste_heat(1:surf_usm_v(0)%ns) )
                READ ( 13 )  tmp_surf_waste_v(0)%t
             ENDIF
             CALL surface_restore_elements(                                 &
                                     surf_usm_v(0)%waste_heat,              &
                                     tmp_surf_waste_v(0)%t,                 &
                                     surf_usm_v(0)%start_index,             & 
                                     start_index_on_file,                   &
                                     end_index_on_file,                     &
                                     nxlc, nysc,                            &
                                     nxlf, nxrf, nysf, nynf,                &
                                     nys_on_file, nyn_on_file,              &
                                     nxl_on_file,nxr_on_file )
                   
          CASE ( 'waste_heat_v(1)' )
             IF ( k == 1 )  THEN
                IF ( .NOT.  ALLOCATED( surf_usm_v(1)%waste_heat ) )         &
                   ALLOCATE( surf_usm_v(1)%waste_heat(1:surf_usm_v(1)%ns) )
                READ ( 13 )  tmp_surf_waste_v(1)%t
             ENDIF
             CALL surface_restore_elements(                                 &
                                     surf_usm_v(1)%waste_heat,              &
                                     tmp_surf_waste_v(1)%t,                 &
                                     surf_usm_v(1)%start_index,             & 
                                     start_index_on_file,                   &
                                     end_index_on_file,                     &
                                     nxlc, nysc,                            &
                                     nxlf, nxrf, nysf, nynf,                &
                                     nys_on_file, nyn_on_file,              &
                                     nxl_on_file,nxr_on_file )

          CASE ( 'waste_heat_v(2)' )
             IF ( k == 1 )  THEN
                IF ( .NOT.  ALLOCATED( surf_usm_v(2)%waste_heat ) )         &
                   ALLOCATE( surf_usm_v(2)%waste_heat(1:surf_usm_v(2)%ns) )
                READ ( 13 )  tmp_surf_waste_v(2)%t
             ENDIF
             CALL surface_restore_elements(                                 &
                                     surf_usm_v(2)%waste_heat,              &
                                     tmp_surf_waste_v(2)%t,                 &
                                     surf_usm_v(2)%start_index,             & 
                                     start_index_on_file,                   &
                                     end_index_on_file,                     &
                                     nxlc, nysc,                            &
                                     nxlf, nxrf, nysf, nynf,                &
                                     nys_on_file, nyn_on_file,              &
                                     nxl_on_file,nxr_on_file )
                   
          CASE ( 'waste_heat_v(3)' )
             IF ( k == 1 )  THEN
                IF ( .NOT.  ALLOCATED( surf_usm_v(3)%waste_heat ) )         &
                   ALLOCATE( surf_usm_v(3)%waste_heat(1:surf_usm_v(3)%ns) )
                READ ( 13 )  tmp_surf_waste_v(3)%t
             ENDIF
             CALL surface_restore_elements(                                 &
                                     surf_usm_v(3)%waste_heat,              &
                                     tmp_surf_waste_v(3)%t,                 &
                                     surf_usm_v(3)%start_index,             & 
                                     start_index_on_file,                   &
                                     end_index_on_file,                     &
                                     nxlc, nysc,                            &
                                     nxlf, nxrf, nysf, nynf,                &
                                     nys_on_file, nyn_on_file,              &
                                     nxl_on_file,nxr_on_file )

          CASE ( 't_wall_h' )
             IF ( k == 1 )  THEN
                IF ( .NOT.  ALLOCATED( t_wall_h_1 ) )                       &
                   ALLOCATE( t_wall_h_1(nzb_wall:nzt_wall+1,                &
                                        1:surf_usm_h%ns) )
                READ ( 13 )  tmp_wall_h
             ENDIF
             CALL surface_restore_elements(                                 &
                                     t_wall_h_1, tmp_wall_h,                &
                                     surf_usm_h%start_index,                & 
                                     start_index_on_file,                   &
                                     end_index_on_file,                     &
                                     nxlc, nysc,                            &
                                     nxlf, nxrf, nysf, nynf,                &
                                     nys_on_file, nyn_on_file,              &
                                     nxl_on_file,nxr_on_file )

          CASE ( 't_wall_v(0)' )
             IF ( k == 1 )  THEN
                IF ( .NOT.  ALLOCATED( t_wall_v_1(0)%t ) )                  &
                   ALLOCATE( t_wall_v_1(0)%t(nzb_wall:nzt_wall+1,           &
                                             1:surf_usm_v(0)%ns) )
                READ ( 13 )  tmp_wall_v(0)%t
             ENDIF
             CALL surface_restore_elements(                                 &
                                     t_wall_v_1(0)%t, tmp_wall_v(0)%t,      &
                                     surf_usm_v(0)%start_index,             &  
                                     start_index_on_file,                   &
                                     end_index_on_file,                     &
                                     nxlc, nysc,                            &
                                     nxlf, nxrf, nysf, nynf,                &
                                     nys_on_file, nyn_on_file,              &
                                     nxl_on_file,nxr_on_file )

          CASE ( 't_wall_v(1)' )
             IF ( k == 1 )  THEN
                IF ( .NOT.  ALLOCATED( t_wall_v_1(1)%t ) )                  &
                   ALLOCATE( t_wall_v_1(1)%t(nzb_wall:nzt_wall+1,           &
                                             1:surf_usm_v(1)%ns) )
                READ ( 13 )  tmp_wall_v(1)%t
             ENDIF
             CALL surface_restore_elements(                                 &
                                     t_wall_v_1(1)%t, tmp_wall_v(1)%t,      &
                                     surf_usm_v(1)%start_index,             &  
                                     start_index_on_file,                   &
                                     end_index_on_file,                     &
                                     nxlc, nysc,                            &
                                     nxlf, nxrf, nysf, nynf,                &
                                     nys_on_file, nyn_on_file,              &
                                     nxl_on_file,nxr_on_file )

          CASE ( 't_wall_v(2)' )
             IF ( k == 1 )  THEN
                IF ( .NOT.  ALLOCATED( t_wall_v_1(2)%t ) )                  &
                   ALLOCATE( t_wall_v_1(2)%t(nzb_wall:nzt_wall+1,           &
                                             1:surf_usm_v(2)%ns) )
                READ ( 13 )  tmp_wall_v(2)%t
             ENDIF
             CALL surface_restore_elements(                                 &
                                     t_wall_v_1(2)%t, tmp_wall_v(2)%t,      &
                                     surf_usm_v(2)%start_index,             & 
                                     start_index_on_file,                   &
                                     end_index_on_file ,                    &
                                     nxlc, nysc,                            &
                                     nxlf, nxrf, nysf, nynf,                &
                                     nys_on_file, nyn_on_file,              &
                                     nxl_on_file,nxr_on_file )

          CASE ( 't_wall_v(3)' )
             IF ( k == 1 )  THEN
                IF ( .NOT.  ALLOCATED( t_wall_v_1(3)%t ) )                  &
                   ALLOCATE( t_wall_v_1(3)%t(nzb_wall:nzt_wall+1,           &
                                             1:surf_usm_v(3)%ns) )
                READ ( 13 )  tmp_wall_v(3)%t
             ENDIF
             CALL surface_restore_elements(                                 &
                                     t_wall_v_1(3)%t, tmp_wall_v(3)%t,      &
                                     surf_usm_v(3)%start_index,             &   
                                     start_index_on_file,                   &
                                     end_index_on_file,                     &
                                     nxlc, nysc,                            &
                                     nxlf, nxrf, nysf, nynf,                &
                                     nys_on_file, nyn_on_file,              &
                                     nxl_on_file,nxr_on_file )

          CASE ( 't_green_h' )
             IF ( k == 1 )  THEN
                IF ( .NOT.  ALLOCATED( t_green_h_1 ) )                      &
                   ALLOCATE( t_green_h_1(nzb_wall:nzt_wall+1,               &
                                         1:surf_usm_h%ns) )
                READ ( 13 )  tmp_green_h
             ENDIF
             CALL surface_restore_elements(                                 &
                                     t_green_h_1, tmp_green_h,              &
                                     surf_usm_h%start_index,                & 
                                     start_index_on_file,                   &
                                     end_index_on_file,                     &
                                     nxlc, nysc,                            &
                                     nxlf, nxrf, nysf, nynf,                &
                                     nys_on_file, nyn_on_file,              &
                                     nxl_on_file,nxr_on_file )

          CASE ( 't_green_v(0)' )
             IF ( k == 1 )  THEN
                IF ( .NOT.  ALLOCATED( t_green_v_1(0)%t ) )                 &
                   ALLOCATE( t_green_v_1(0)%t(nzb_wall:nzt_wall+1,          &
                                              1:surf_usm_v(0)%ns) )
                READ ( 13 )  tmp_green_v(0)%t
             ENDIF
             CALL surface_restore_elements(                                 &
                                     t_green_v_1(0)%t, tmp_green_v(0)%t,    &
                                     surf_usm_v(0)%start_index,             & 
                                     start_index_on_file,                   &
                                     end_index_on_file,                     &
                                     nxlc, nysc,                            &
                                     nxlf, nxrf, nysf, nynf,                &
                                     nys_on_file, nyn_on_file,              &
                                     nxl_on_file,nxr_on_file )

          CASE ( 't_green_v(1)' )
             IF ( k == 1 )  THEN
                IF ( .NOT.  ALLOCATED( t_green_v_1(1)%t ) )                 &
                   ALLOCATE( t_green_v_1(1)%t(nzb_wall:nzt_wall+1,          &
                                              1:surf_usm_v(1)%ns) )
                READ ( 13 )  tmp_green_v(1)%t
             ENDIF
             CALL surface_restore_elements(                                 &
                                     t_green_v_1(1)%t, tmp_green_v(1)%t,    &
                                     surf_usm_v(1)%start_index,             & 
                                     start_index_on_file,                   &
                                     end_index_on_file,                     &
                                     nxlc, nysc,                            &
                                     nxlf, nxrf, nysf, nynf,                &
                                     nys_on_file, nyn_on_file,              &
                                     nxl_on_file,nxr_on_file )

          CASE ( 't_green_v(2)' )
             IF ( k == 1 )  THEN
                IF ( .NOT.  ALLOCATED( t_green_v_1(2)%t ) )                 &
                   ALLOCATE( t_green_v_1(2)%t(nzb_wall:nzt_wall+1,          &
                                              1:surf_usm_v(2)%ns) )
                READ ( 13 )  tmp_green_v(2)%t
             ENDIF
             CALL surface_restore_elements(                                 &
                                     t_green_v_1(2)%t, tmp_green_v(2)%t,    &
                                     surf_usm_v(2)%start_index,             & 
                                     start_index_on_file,                   &
                                     end_index_on_file ,                    &
                                     nxlc, nysc,                            &
                                     nxlf, nxrf, nysf, nynf,                &
                                     nys_on_file, nyn_on_file,              &
                                     nxl_on_file,nxr_on_file )

          CASE ( 't_green_v(3)' )
             IF ( k == 1 )  THEN
                IF ( .NOT.  ALLOCATED( t_green_v_1(3)%t ) )                 &
                   ALLOCATE( t_green_v_1(3)%t(nzb_wall:nzt_wall+1,          &
                                              1:surf_usm_v(3)%ns) )
                READ ( 13 )  tmp_green_v(3)%t
             ENDIF
             CALL surface_restore_elements(                                 &
                                     t_green_v_1(3)%t, tmp_green_v(3)%t,    &
                                     surf_usm_v(3)%start_index,             &  
                                     start_index_on_file,                   &
                                     end_index_on_file,                     &
                                     nxlc, nysc,                            &
                                     nxlf, nxrf, nysf, nynf,                &
                                     nys_on_file, nyn_on_file,              &
                                     nxl_on_file,nxr_on_file )

          CASE ( 't_window_h' )
             IF ( k == 1 )  THEN
                IF ( .NOT.  ALLOCATED( t_window_h_1 ) )                     &
                   ALLOCATE( t_window_h_1(nzb_wall:nzt_wall+1,              &
                                          1:surf_usm_h%ns) )
                READ ( 13 )  tmp_window_h
             ENDIF
             CALL surface_restore_elements(                                 &
                                     t_window_h_1, tmp_window_h,            &
                                     surf_usm_h%start_index,                & 
                                     start_index_on_file,                   &
                                     end_index_on_file,                     &
                                     nxlc, nysc,                            &
                                     nxlf, nxrf, nysf, nynf,                &
                                     nys_on_file, nyn_on_file,              &
                                     nxl_on_file, nxr_on_file )

          CASE ( 't_window_v(0)' )
             IF ( k == 1 )  THEN
                IF ( .NOT.  ALLOCATED( t_window_v_1(0)%t ) )                &
                   ALLOCATE( t_window_v_1(0)%t(nzb_wall:nzt_wall+1,         &
                                               1:surf_usm_v(0)%ns) )
                READ ( 13 )  tmp_window_v(0)%t
             ENDIF
             CALL surface_restore_elements(                                 &
                                     t_window_v_1(0)%t,                     & 
                                     tmp_window_v(0)%t,                     &
                                     surf_usm_v(0)%start_index,             &
                                     start_index_on_file,                   &
                                     end_index_on_file,                     &
                                     nxlc, nysc,                            &
                                     nxlf, nxrf, nysf, nynf,                &
                                     nys_on_file, nyn_on_file,              &
                                     nxl_on_file,nxr_on_file )

          CASE ( 't_window_v(1)' )
             IF ( k == 1 )  THEN
                IF ( .NOT.  ALLOCATED( t_window_v_1(1)%t ) )                &
                   ALLOCATE( t_window_v_1(1)%t(nzb_wall:nzt_wall+1,         &
                                               1:surf_usm_v(1)%ns) )
                READ ( 13 )  tmp_window_v(1)%t
             ENDIF
             CALL surface_restore_elements(                                 &
                                     t_window_v_1(1)%t,                     & 
                                     tmp_window_v(1)%t,                     &
                                     surf_usm_v(1)%start_index,             & 
                                     start_index_on_file,                   &
                                     end_index_on_file,                     &
                                     nxlc, nysc,                            &
                                     nxlf, nxrf, nysf, nynf,                &
                                     nys_on_file, nyn_on_file,              &
                                     nxl_on_file,nxr_on_file )

          CASE ( 't_window_v(2)' )
             IF ( k == 1 )  THEN
                IF ( .NOT.  ALLOCATED( t_window_v_1(2)%t ) )                &
                   ALLOCATE( t_window_v_1(2)%t(nzb_wall:nzt_wall+1,         &
                                               1:surf_usm_v(2)%ns) )
                READ ( 13 )  tmp_window_v(2)%t
             ENDIF
             CALL surface_restore_elements(                                 &
                                     t_window_v_1(2)%t,                     & 
                                     tmp_window_v(2)%t,                     &
                                     surf_usm_v(2)%start_index,             &  
                                     start_index_on_file,                   &
                                     end_index_on_file ,                    &
                                     nxlc, nysc,                            &
                                     nxlf, nxrf, nysf, nynf,                &
                                     nys_on_file, nyn_on_file,              &
                                     nxl_on_file,nxr_on_file )

          CASE ( 't_window_v(3)' )
             IF ( k == 1 )  THEN
                IF ( .NOT.  ALLOCATED( t_window_v_1(3)%t ) )                &
                   ALLOCATE( t_window_v_1(3)%t(nzb_wall:nzt_wall+1,1:surf_usm_v(3)%ns) )
                READ ( 13 )  tmp_window_v(3)%t
             ENDIF
             CALL surface_restore_elements(                                 &
                                     t_window_v_1(3)%t,                     & 
                                     tmp_window_v(3)%t,                     &
                                     surf_usm_v(3)%start_index,             & 
                                     start_index_on_file,                   &
                                     end_index_on_file,                     &
                                     nxlc, nysc,                            &
                                     nxlf, nxrf, nysf, nynf,                &
                                     nys_on_file, nyn_on_file,              &
                                     nxl_on_file,nxr_on_file )

          CASE DEFAULT

             found = .FALSE.

       END SELECT

    END SUBROUTINE usm_rrd_local

    
!------------------------------------------------------------------------------!
! Description:
! ------------
!
!> This subroutine reads walls, roofs and land categories and it parameters
!> from input files.
!------------------------------------------------------------------------------!
    SUBROUTINE usm_read_urban_surface_types
    
        USE netcdf_data_input_mod,                                             &
            ONLY:  building_pars_f, building_type_f

        IMPLICIT NONE

        CHARACTER(12)                                         :: wtn
        INTEGER(iwp)                                          :: wtc
        REAL(wp), DIMENSION(n_surface_params)                 :: wtp
        LOGICAL                                               :: ascii_file = .FALSE.
        INTEGER(iwp), DIMENSION(0:17, nysg:nyng, nxlg:nxrg)   :: usm_par
        REAL(wp), DIMENSION(1:14, nysg:nyng, nxlg:nxrg)       :: usm_val
        INTEGER(iwp)                                          :: k, l, iw, jw, kw, it, ip, ii, ij, m
        INTEGER(iwp)                                          :: i, j
        INTEGER(iwp)                                          :: nz, roof, dirwe, dirsn
        INTEGER(iwp)                                          :: category
        INTEGER(iwp)                                          :: weheight1, wecat1, snheight1, sncat1
        INTEGER(iwp)                                          :: weheight2, wecat2, snheight2, sncat2
        INTEGER(iwp)                                          :: weheight3, wecat3, snheight3, sncat3
        REAL(wp)                                              :: height, albedo, thick
        REAL(wp)                                              :: wealbedo1, wethick1, snalbedo1, snthick1
        REAL(wp)                                              :: wealbedo2, wethick2, snalbedo2, snthick2
        REAL(wp)                                              :: wealbedo3, wethick3, snalbedo3, snthick3


        IF ( debug_output )  CALL debug_message( 'usm_read_urban_surface_types', 'start' )
!
!--     If building_pars or building_type are already read from static input 
!--     file, skip reading ASCII file. 
        IF ( building_type_f%from_file  .OR.  building_pars_f%from_file )      &
           RETURN
!
!--     Check if ASCII input file exists. If not, return and initialize USM
!--     with default settings.
        INQUIRE( FILE = 'SURFACE_PARAMETERS' // coupling_char,                 &
                 EXIST = ascii_file )
                 
        IF ( .NOT. ascii_file )  RETURN

!
!--     read categories of walls and their parameters
        DO  ii = 0, io_blocks-1
            IF ( ii == io_group )  THEN
!
!--             open urban surface file 
                OPEN( 151, file='SURFACE_PARAMETERS'//coupling_char, action='read', &
                           status='old', form='formatted', err=15 )
!
!--             first test and get n_surface_types
                k = 0
                l = 0
                DO
                    l = l+1
                    READ( 151, *, err=11, end=12 )  wtc, wtp, wtn
                    k = k+1
                    CYCLE
 11                 CONTINUE
                ENDDO
 12             n_surface_types = k
                ALLOCATE( surface_type_names(n_surface_types) )
                ALLOCATE( surface_type_codes(n_surface_types) )
                ALLOCATE( surface_params(n_surface_params, n_surface_types) )
!
!--             real reading
                rewind( 151 )
                k = 0
                DO
                    READ( 151, *, err=13, end=14 )  wtc, wtp, wtn
                    k = k+1
                    surface_type_codes(k) = wtc
                    surface_params(:,k) = wtp
                    surface_type_names(k) = wtn
                    CYCLE
13                  WRITE(6,'(i3,a,2i5)') myid, 'readparams2 error k=', k
                    FLUSH(6)
                    CONTINUE
                ENDDO
 14             CLOSE(151)
                CYCLE
 15             message_string = 'file SURFACE_PARAMETERS'//TRIM(coupling_char)//' does not exist'
                CALL message( 'usm_read_urban_surface_types', 'PA0513', 1, 2, 0, 6, 0 )
            ENDIF
        ENDDO
    
!
!--     read types of surfaces
        usm_par = 0
        DO  ii = 0, io_blocks-1
            IF ( ii == io_group )  THEN

!
!--             open csv urban surface file 
                OPEN( 151, file='URBAN_SURFACE'//TRIM(coupling_char), action='read', &
                      status='old', form='formatted', err=23 )
                
                l = 0
                DO
                    l = l+1
!
!--                 i, j, height, nz, roof, dirwe, dirsn, category, soilcat,
!--                 weheight1, wecat1, snheight1, sncat1, weheight2, wecat2, snheight2, sncat2,
!--                 weheight3, wecat3, snheight3, sncat3
                    READ( 151, *, err=21, end=25 )  i, j, height, nz, roof, dirwe, dirsn,            &
                                            category, albedo, thick,                                 &
                                            weheight1, wecat1, wealbedo1, wethick1,                  &
                                            weheight2, wecat2, wealbedo2, wethick2,                  &
                                            weheight3, wecat3, wealbedo3, wethick3,                  &
                                            snheight1, sncat1, snalbedo1, snthick1,                  &
                                            snheight2, sncat2, snalbedo2, snthick2,                  &
                                            snheight3, sncat3, snalbedo3, snthick3

                    IF ( i >= nxlg  .AND.  i <= nxrg  .AND.  j >= nysg  .AND.  j <= nyng )  THEN
!
!--                     write integer variables into array
                        usm_par(:,j,i) = (/1, nz, roof, dirwe, dirsn, category,                      &
                                          weheight1, wecat1, weheight2, wecat2, weheight3, wecat3,   &
                                          snheight1, sncat1, snheight2, sncat2, snheight3, sncat3 /)
!
!--                     write real values into array
                        usm_val(:,j,i) = (/ albedo, thick,                                           &
                                           wealbedo1, wethick1, wealbedo2, wethick2,                 &
                                           wealbedo3, wethick3, snalbedo1, snthick1,                 &
                                           snalbedo2, snthick2, snalbedo3, snthick3 /)
                    ENDIF
                    CYCLE
 21                 WRITE (message_string, "(A,I5)") 'errors in file URBAN_SURFACE'//TRIM(coupling_char)//' on line ', l
                    CALL message( 'usm_read_urban_surface_types', 'PA0512', 0, 1, 0, 6, 0 )
                ENDDO
          
 23             message_string = 'file URBAN_SURFACE'//TRIM(coupling_char)//' does not exist'
                CALL message( 'usm_read_urban_surface_types', 'PA0514', 1, 2, 0, 6, 0 )

 25             CLOSE( 151 )

            ENDIF
#if defined( __parallel )
            CALL MPI_BARRIER( comm2d, ierr )
#endif
        ENDDO
        
!
!--     check completeness and formal correctness of the data
        DO i = nxlg, nxrg
            DO j = nysg, nyng
                IF ( usm_par(0,j,i) /= 0  .AND.  (        &  !< incomplete data,supply default values later
                     usm_par(1,j,i) < nzb  .OR.           &
                     usm_par(1,j,i) > nzt  .OR.           &  !< incorrect height (nz < nzb  .OR.  nz > nzt)
                     usm_par(2,j,i) < 0  .OR.             &
                     usm_par(2,j,i) > 1  .OR.             &  !< incorrect roof sign
                     usm_par(3,j,i) < nzb-nzt  .OR.       & 
                     usm_par(3,j,i) > nzt-nzb  .OR.       &  !< incorrect west-east wall direction sign
                     usm_par(4,j,i) < nzb-nzt  .OR.       &
                     usm_par(4,j,i) > nzt-nzb  .OR.       &  !< incorrect south-north wall direction sign
                     usm_par(6,j,i) < nzb  .OR.           & 
                     usm_par(6,j,i) > nzt  .OR.           &  !< incorrect pedestrian level height for west-east wall
                     usm_par(8,j,i) > nzt  .OR.           &
                     usm_par(10,j,i) > nzt  .OR.          &  !< incorrect wall or roof level height for west-east wall
                     usm_par(12,j,i) < nzb  .OR.          & 
                     usm_par(12,j,i) > nzt  .OR.          &  !< incorrect pedestrian level height for south-north wall
                     usm_par(14,j,i) > nzt  .OR.          &
                     usm_par(16,j,i) > nzt                &  !< incorrect wall or roof level height for south-north wall
                    ) )  THEN
!
!--                 incorrect input data
                    WRITE (message_string, "(A,2I5)") 'missing or incorrect data in file URBAN_SURFACE'// &
                                                       TRIM(coupling_char)//' for i,j=', i,j
                    CALL message( 'usm_read_urban_surface', 'PA0504', 1, 2, 0, 6, 0 )
                ENDIF
                
            ENDDO
        ENDDO
!        
!--     Assign the surface types to the respective data type. 
!--     First, for horizontal upward-facing surfaces. 
!--     Further, set flag indicating that albedo is initialized via ASCII
!--     format, else it would be overwritten in the radiation model.
        surf_usm_h%albedo_from_ascii = .TRUE.
        DO  m = 1, surf_usm_h%ns
           iw = surf_usm_h%i(m)
           jw = surf_usm_h%j(m)
           kw = surf_usm_h%k(m)

           IF ( usm_par(5,jw,iw) == 0 )  THEN

              IF ( zu(kw) >= roof_height_limit )  THEN
                 surf_usm_h%isroof_surf(m)   = .TRUE.
                 surf_usm_h%surface_types(m) = roof_category         !< default category for root surface
              ELSE
                 surf_usm_h%isroof_surf(m)   = .FALSE.
                 surf_usm_h%surface_types(m) = land_category         !< default category for land surface
              ENDIF

              surf_usm_h%albedo(m,:)    = -1.0_wp
              surf_usm_h%thickness_wall(m) = -1.0_wp
              surf_usm_h%thickness_green(m) = -1.0_wp
              surf_usm_h%thickness_window(m) = -1.0_wp
           ELSE
              IF ( usm_par(2,jw,iw)==0 )  THEN
                 surf_usm_h%isroof_surf(m)    = .FALSE.
                 surf_usm_h%thickness_wall(m) = -1.0_wp
                 surf_usm_h%thickness_window(m) = -1.0_wp
                 surf_usm_h%thickness_green(m)  = -1.0_wp
              ELSE
                 surf_usm_h%isroof_surf(m)    = .TRUE.
                 surf_usm_h%thickness_wall(m) = usm_val(2,jw,iw)
                 surf_usm_h%thickness_window(m) = usm_val(2,jw,iw)
                 surf_usm_h%thickness_green(m)  = usm_val(2,jw,iw)
              ENDIF
              surf_usm_h%surface_types(m) = usm_par(5,jw,iw)
              surf_usm_h%albedo(m,:)   = usm_val(1,jw,iw)
              surf_usm_h%transmissivity(m)    = 0.0_wp
           ENDIF
!
!--        Find the type position
           it = surf_usm_h%surface_types(m)
           ip = -99999
           DO k = 1, n_surface_types
              IF ( surface_type_codes(k) == it )  THEN
                 ip = k
                 EXIT
              ENDIF
           ENDDO
           IF ( ip == -99999 )  THEN
!
!--           land/roof category not found
              WRITE (9,"(A,I5,A,3I5)") 'land/roof category ', it,     &
                                       ' not found  for i,j,k=', iw,jw,kw
              FLUSH(9)
              IF ( surf_usm_h%isroof_surf(m) ) THEN
                 category = roof_category
              ELSE
                 category = land_category
              ENDIF
              DO k = 1, n_surface_types
                 IF ( surface_type_codes(k) == roof_category ) THEN
                    ip = k
                    EXIT
                 ENDIF
              ENDDO
              IF ( ip == -99999 )  THEN
!
!--              default land/roof category not found
                 WRITE (9,"(A,I5,A,3I5)") 'Default land/roof category', category, ' not found!'
                 FLUSH(9)
                 ip = 1
              ENDIF
           ENDIF
!
!--        Albedo
           IF ( surf_usm_h%albedo(m,ind_veg_wall) < 0.0_wp )  THEN
              surf_usm_h%albedo(m,:) = surface_params(ialbedo,ip)
           ENDIF
!
!--        Albedo type is 0 (custom), others are replaced later
           surf_usm_h%albedo_type(m,:) = 0
!
!--        Transmissivity
           IF ( surf_usm_h%transmissivity(m) < 0.0_wp )  THEN
              surf_usm_h%transmissivity(m) = 0.0_wp
           ENDIF
!
!--        emissivity of the wall
           surf_usm_h%emissivity(m,:) = surface_params(iemiss,ip)
!            
!--        heat conductivity λS between air and wall ( W m−2 K−1 )
           surf_usm_h%lambda_surf(m) = surface_params(ilambdas,ip)
           surf_usm_h%lambda_surf_window(m) = surface_params(ilambdas,ip)
           surf_usm_h%lambda_surf_green(m)  = surface_params(ilambdas,ip)
!            
!--        roughness length for momentum, heat and humidity
           surf_usm_h%z0(m) = surface_params(irough,ip)
           surf_usm_h%z0h(m) = surface_params(iroughh,ip)
           surf_usm_h%z0q(m) = surface_params(iroughh,ip)
!
!--        Surface skin layer heat capacity (J m−2 K−1 )
           surf_usm_h%c_surface(m) = surface_params(icsurf,ip)
           surf_usm_h%c_surface_window(m) = surface_params(icsurf,ip)
           surf_usm_h%c_surface_green(m)  = surface_params(icsurf,ip)
!            
!--        wall material parameters:
!--        thickness of the wall (m)
!--        missing values are replaced by default value for category
           IF ( surf_usm_h%thickness_wall(m) <= 0.001_wp )  THEN
                surf_usm_h%thickness_wall(m) = surface_params(ithick,ip)
           ENDIF
           IF ( surf_usm_h%thickness_window(m) <= 0.001_wp )  THEN
                surf_usm_h%thickness_window(m) = surface_params(ithick,ip)
           ENDIF
           IF ( surf_usm_h%thickness_green(m) <= 0.001_wp )  THEN
                surf_usm_h%thickness_green(m) = surface_params(ithick,ip)
           ENDIF
!            
!--        volumetric heat capacity rho*C of the wall ( J m−3 K−1 )
           surf_usm_h%rho_c_wall(:,m) = surface_params(irhoC,ip)
           surf_usm_h%rho_c_window(:,m) = surface_params(irhoC,ip)
           surf_usm_h%rho_c_green(:,m)  = surface_params(irhoC,ip)
!            
!--        thermal conductivity λH of the wall (W m−1 K−1 )
           surf_usm_h%lambda_h(:,m) = surface_params(ilambdah,ip)
           surf_usm_h%lambda_h_window(:,m) = surface_params(ilambdah,ip)
           surf_usm_h%lambda_h_green(:,m)  = surface_params(ilambdah,ip)

        ENDDO
!
!--     For vertical surface elements ( 0 -- northward-facing, 1 -- southward-facing,
!--     2 -- eastward-facing, 3 -- westward-facing )
        DO  l = 0, 3
!
!--        Set flag indicating that albedo is initialized via ASCII format. 
!--        Else it would be overwritten in the radiation model.
           surf_usm_v(l)%albedo_from_ascii = .TRUE.
           DO  m = 1, surf_usm_v(l)%ns
              i  = surf_usm_v(l)%i(m)
              j  = surf_usm_v(l)%j(m)
              kw = surf_usm_v(l)%k(m)
              
              IF ( l == 3 )  THEN ! westward facing
                 iw = i
                 jw = j
                 ii = 6
                 ij = 3
              ELSEIF ( l == 2 )  THEN
                 iw = i-1
                 jw = j
                 ii = 6
                 ij = 3
              ELSEIF ( l == 1 )  THEN
                 iw = i
                 jw = j
                 ii = 12
                 ij = 9
              ELSEIF ( l == 0 )  THEN
                 iw = i
                 jw = j-1
                 ii = 12
                 ij = 9
              ENDIF

              IF ( iw < 0 .OR. jw < 0 ) THEN
!
!--              wall on west or south border of the domain - assign default category
                 IF ( kw <= roof_height_limit ) THEN
                     surf_usm_v(l)%surface_types(m) = wall_category   !< default category for wall surface in wall zone
                 ELSE
                     surf_usm_v(l)%surface_types(m) = roof_category   !< default category for wall surface in roof zone
                 END IF
                 surf_usm_v(l)%albedo(m,:)         = -1.0_wp
                 surf_usm_v(l)%thickness_wall(m)   = -1.0_wp
                 surf_usm_v(l)%thickness_window(m) = -1.0_wp
                 surf_usm_v(l)%thickness_green(m)  = -1.0_wp
                 surf_usm_v(l)%transmissivity(m)   = -1.0_wp
              ELSE IF ( kw <= usm_par(ii,jw,iw) )  THEN
!
!--                 pedestrian zone
                 IF ( usm_par(ii+1,jw,iw) == 0 )  THEN
                     surf_usm_v(l)%surface_types(m)  = pedestrian_category   !< default category for wall surface in
                                                                             !<pedestrian zone
                     surf_usm_v(l)%albedo(m,:)         = -1.0_wp
                     surf_usm_v(l)%thickness_wall(m)   = -1.0_wp
                     surf_usm_v(l)%thickness_window(m) = -1.0_wp
                     surf_usm_v(l)%thickness_green(m)  = -1.0_wp
                     surf_usm_v(l)%transmissivity(m)   = -1.0_wp
                 ELSE
                     surf_usm_v(l)%surface_types(m)    = usm_par(ii+1,jw,iw)
                     surf_usm_v(l)%albedo(m,:)         = usm_val(ij,jw,iw)
                     surf_usm_v(l)%thickness_wall(m)   = usm_val(ij+1,jw,iw)
                     surf_usm_v(l)%thickness_window(m) = usm_val(ij+1,jw,iw)
                     surf_usm_v(l)%thickness_green(m)  = usm_val(ij+1,jw,iw)
                     surf_usm_v(l)%transmissivity(m)   = 0.0_wp
                 ENDIF
              ELSE IF ( kw <= usm_par(ii+2,jw,iw) )  THEN
!
!--              wall zone
                 IF ( usm_par(ii+3,jw,iw) == 0 )  THEN
                     surf_usm_v(l)%surface_types(m)    = wall_category         !< default category for wall surface
                     surf_usm_v(l)%albedo(m,:)         = -1.0_wp
                     surf_usm_v(l)%thickness_wall(m)   = -1.0_wp
                     surf_usm_v(l)%thickness_window(m) = -1.0_wp
                     surf_usm_v(l)%thickness_green(m)  = -1.0_wp
                     surf_usm_v(l)%transmissivity(m)   = -1.0_wp
                 ELSE
                     surf_usm_v(l)%surface_types(m)    = usm_par(ii+3,jw,iw)
                     surf_usm_v(l)%albedo(m,:)         = usm_val(ij+2,jw,iw)
                     surf_usm_v(l)%thickness_wall(m)   = usm_val(ij+3,jw,iw)
                     surf_usm_v(l)%thickness_window(m) = usm_val(ij+3,jw,iw)
                     surf_usm_v(l)%thickness_green(m)  = usm_val(ij+3,jw,iw)
                     surf_usm_v(l)%transmissivity(m)   = 0.0_wp
                 ENDIF
              ELSE IF ( kw <= usm_par(ii+4,jw,iw) )  THEN
!
!--              roof zone
                 IF ( usm_par(ii+5,jw,iw) == 0 )  THEN
                     surf_usm_v(l)%surface_types(m)    = roof_category         !< default category for roof surface
                     surf_usm_v(l)%albedo(m,:)         = -1.0_wp
                     surf_usm_v(l)%thickness_wall(m)   = -1.0_wp
                     surf_usm_v(l)%thickness_window(m) = -1.0_wp
                     surf_usm_v(l)%thickness_green(m)  = -1.0_wp
                     surf_usm_v(l)%transmissivity(m)   = -1.0_wp
                 ELSE
                     surf_usm_v(l)%surface_types(m)    = usm_par(ii+5,jw,iw)
                     surf_usm_v(l)%albedo(m,:)         = usm_val(ij+4,jw,iw)
                     surf_usm_v(l)%thickness_wall(m)   = usm_val(ij+5,jw,iw)
                     surf_usm_v(l)%thickness_window(m) = usm_val(ij+5,jw,iw)
                     surf_usm_v(l)%thickness_green(m)  = usm_val(ij+5,jw,iw)
                     surf_usm_v(l)%transmissivity(m)   = 0.0_wp
                 ENDIF
              ELSE
                 WRITE(9,*) 'Problem reading USM data:'
                 WRITE(9,*) l,i,j,kw,topo_top_ind(j,i,0)
                 WRITE(9,*) ii,iw,jw,kw,topo_top_ind(jw,iw,0)
                 WRITE(9,*) usm_par(ii,jw,iw),usm_par(ii+1,jw,iw)
                 WRITE(9,*) usm_par(ii+2,jw,iw),usm_par(ii+3,jw,iw)
                 WRITE(9,*) usm_par(ii+4,jw,iw),usm_par(ii+5,jw,iw)
                 WRITE(9,*) kw,roof_height_limit,wall_category,roof_category
                 FLUSH(9)
!
!--              supply the default category
                 IF ( kw <= roof_height_limit ) THEN
                     surf_usm_v(l)%surface_types(m) = wall_category   !< default category for wall surface in wall zone
                 ELSE
                     surf_usm_v(l)%surface_types(m) = roof_category   !< default category for wall surface in roof zone
                 END IF
                 surf_usm_v(l)%albedo(m,:)         = -1.0_wp
                 surf_usm_v(l)%thickness_wall(m)   = -1.0_wp
                 surf_usm_v(l)%thickness_window(m) = -1.0_wp
                 surf_usm_v(l)%thickness_green(m)  = -1.0_wp
                 surf_usm_v(l)%transmissivity(m)   = -1.0_wp
              ENDIF
!
!--           Find the type position
              it = surf_usm_v(l)%surface_types(m)
              ip = -99999
              DO k = 1, n_surface_types
                 IF ( surface_type_codes(k) == it )  THEN
                    ip = k
                    EXIT
                 ENDIF
              ENDDO
              IF ( ip == -99999 )  THEN
!
!--              wall category not found
                 WRITE (9, "(A,I7,A,3I5)") 'wall category ', it,  &
                                           ' not found  for i,j,k=', iw,jw,kw
                 FLUSH(9)
                 category = wall_category 
                 DO k = 1, n_surface_types
                    IF ( surface_type_codes(k) == category ) THEN
                       ip = k
                       EXIT
                    ENDIF
                 ENDDO
                 IF ( ip == -99999 )  THEN
!
!--                 default wall category not found
                    WRITE (9, "(A,I5,A,3I5)") 'Default wall category', category, ' not found!'
                    FLUSH(9)
                    ip = 1
                 ENDIF
              ENDIF

!
!--           Albedo
              IF ( surf_usm_v(l)%albedo(m,ind_veg_wall) < 0.0_wp )  THEN
                 surf_usm_v(l)%albedo(m,:) = surface_params(ialbedo,ip)
              ENDIF
!--           Albedo type is 0 (custom), others are replaced later
              surf_usm_v(l)%albedo_type(m,:) = 0
!--           Transmissivity of the windows
              IF ( surf_usm_v(l)%transmissivity(m) < 0.0_wp )  THEN
                 surf_usm_v(l)%transmissivity(m) = 0.0_wp
              ENDIF
!
!--           emissivity of the wall
              surf_usm_v(l)%emissivity(:,m) = surface_params(iemiss,ip)
!            
!--           heat conductivity lambda S between air and wall ( W m-2 K-1 )
              surf_usm_v(l)%lambda_surf(m) = surface_params(ilambdas,ip)
              surf_usm_v(l)%lambda_surf_window(m) = surface_params(ilambdas,ip)
              surf_usm_v(l)%lambda_surf_green(m) = surface_params(ilambdas,ip)
!            
!--           roughness length
              surf_usm_v(l)%z0(m) = surface_params(irough,ip)
              surf_usm_v(l)%z0h(m) = surface_params(iroughh,ip)
              surf_usm_v(l)%z0q(m) = surface_params(iroughh,ip)
!            
!--           Surface skin layer heat capacity (J m-2 K-1 )
              surf_usm_v(l)%c_surface(m) = surface_params(icsurf,ip)
              surf_usm_v(l)%c_surface_window(m) = surface_params(icsurf,ip)
              surf_usm_v(l)%c_surface_green(m) = surface_params(icsurf,ip)
!            
!--           wall material parameters:
!--           thickness of the wall (m)
!--           missing values are replaced by default value for category
              IF ( surf_usm_v(l)%thickness_wall(m) <= 0.001_wp )  THEN
                   surf_usm_v(l)%thickness_wall(m) = surface_params(ithick,ip)
              ENDIF
              IF ( surf_usm_v(l)%thickness_window(m) <= 0.001_wp )  THEN
                   surf_usm_v(l)%thickness_window(m) = surface_params(ithick,ip)
              ENDIF
              IF ( surf_usm_v(l)%thickness_green(m) <= 0.001_wp )  THEN
                   surf_usm_v(l)%thickness_green(m) = surface_params(ithick,ip)
              ENDIF
!
!--           volumetric heat capacity rho*C of the wall ( J m-3 K-1 )
              surf_usm_v(l)%rho_c_wall(:,m) = surface_params(irhoC,ip)
              surf_usm_v(l)%rho_c_window(:,m) = surface_params(irhoC,ip)
              surf_usm_v(l)%rho_c_green(:,m) = surface_params(irhoC,ip)
!            
!--           thermal conductivity lambda H of the wall (W m-1 K-1 )
              surf_usm_v(l)%lambda_h(:,m) = surface_params(ilambdah,ip)
              surf_usm_v(l)%lambda_h_window(:,m) = surface_params(ilambdah,ip)
              surf_usm_v(l)%lambda_h_green(:,m) = surface_params(ilambdah,ip)

           ENDDO
        ENDDO 

!
!--     Initialize wall layer thicknesses. Please note, this will be removed
!--     after migration to Palm input data standard.  
        DO k = nzb_wall, nzt_wall
           zwn(k) = zwn_default(k)
           zwn_green(k) = zwn_default_green(k)
           zwn_window(k) = zwn_default_window(k)
        ENDDO
!
!--     apply for all particular surface grids. First for horizontal surfaces
        DO  m = 1, surf_usm_h%ns
           surf_usm_h%zw(:,m) = zwn(:) * surf_usm_h%thickness_wall(m)
           surf_usm_h%zw_green(:,m) = zwn_green(:) * surf_usm_h%thickness_green(m)
           surf_usm_h%zw_window(:,m) = zwn_window(:) * surf_usm_h%thickness_window(m)
        ENDDO
        DO  l = 0, 3
           DO  m = 1, surf_usm_v(l)%ns
              surf_usm_v(l)%zw(:,m) = zwn(:) * surf_usm_v(l)%thickness_wall(m)
              surf_usm_v(l)%zw_green(:,m) = zwn_green(:) * surf_usm_v(l)%thickness_green(m)
              surf_usm_v(l)%zw_window(:,m) = zwn_window(:) * surf_usm_v(l)%thickness_window(m)
           ENDDO
        ENDDO

        IF ( debug_output )  CALL debug_message( 'usm_read_urban_surface_types', 'end' )
   
    END SUBROUTINE usm_read_urban_surface_types


!------------------------------------------------------------------------------!
! Description:
! ------------
!
!> This function advances through the list of local surfaces to find given
!> x, y, d, z coordinates
!------------------------------------------------------------------------------!
    PURE FUNCTION find_surface( x, y, z, d ) result(isurfl)

        INTEGER(iwp), INTENT(in)                :: x, y, z, d
        INTEGER(iwp)                            :: isurfl
        INTEGER(iwp)                            :: isx, isy, isz

        IF ( d == 0 ) THEN
           DO  isurfl = 1, surf_usm_h%ns
              isx = surf_usm_h%i(isurfl)
              isy = surf_usm_h%j(isurfl)
              isz = surf_usm_h%k(isurfl)
              IF ( isx==x .and. isy==y .and. isz==z )  RETURN
           ENDDO
        ELSE
           DO  isurfl = 1, surf_usm_v(d-1)%ns
              isx = surf_usm_v(d-1)%i(isurfl)
              isy = surf_usm_v(d-1)%j(isurfl)
              isz = surf_usm_v(d-1)%k(isurfl)
              IF ( isx==x .and. isy==y .and. isz==z )  RETURN
           ENDDO
        ENDIF
!
!--     coordinate not found
        isurfl = -1

    END FUNCTION


!------------------------------------------------------------------------------!
! Description:
! ------------
!
!> This subroutine reads temperatures of respective material layers in walls,
!> roofs and ground from input files. Data in the input file must be in
!> standard order, i.e. horizontal surfaces first ordered by x, y and then
!> vertical surfaces ordered by x, y, direction, z
!------------------------------------------------------------------------------!
    SUBROUTINE usm_read_wall_temperature

        INTEGER(iwp)                                          :: i, j, k, d, ii, iline  !> running indices
        INTEGER(iwp)                                          :: isurfl
        REAL(wp)                                              :: rtsurf
        REAL(wp), DIMENSION(nzb_wall:nzt_wall+1)              :: rtwall


        IF ( debug_output )  CALL debug_message( 'usm_read_wall_temperature', 'start' )

        DO  ii = 0, io_blocks-1
            IF ( ii == io_group )  THEN
!
!--             open wall temperature file
                OPEN( 152, file='WALL_TEMPERATURE'//coupling_char, action='read', &
                           status='old', form='formatted', err=15 )

                isurfl = 0
                iline = 1
                DO
                    rtwall = -9999.0_wp  !< for incomplete lines
                    READ( 152, *, err=13, end=14 )  i, j, k, d, rtsurf, rtwall

                    IF ( nxl <= i .and. i <= nxr .and. &
                        nys <= j .and. j <= nyn)  THEN  !< local processor
!--                     identify surface id
                        isurfl = find_surface( i, j, k, d )
                        IF ( isurfl == -1 )  THEN
                            WRITE(message_string, '(a,4i5,a,i5,a)') 'Coordinates (xyzd) ', i, j, k, d, &
                                ' on line ', iline, &
                                ' in file WALL_TEMPERATURE are either not present or out of standard order of surfaces.'
                            CALL message( 'usm_read_wall_temperature', 'PA0521', 1, 2, 0, 6, 0 )
                        ENDIF
!
!--                     assign temperatures
                        IF ( d == 0 ) THEN
                           t_surf_wall_h(isurfl) = rtsurf
                           t_wall_h(:,isurfl) = rtwall(:)
                           t_window_h(:,isurfl) = rtwall(:)
                           t_green_h(:,isurfl) = rtwall(:)
                        ELSE
                           t_surf_wall_v(d-1)%t(isurfl) = rtsurf
                           t_wall_v(d-1)%t(:,isurfl) = rtwall(:)
                           t_window_v(d-1)%t(:,isurfl) = rtwall(:)
                           t_green_v(d-1)%t(:,isurfl) = rtwall(:)
                        ENDIF
                    ENDIF

                    iline = iline + 1
                    CYCLE
 13                 WRITE(message_string, '(a,i5,a)') 'Error reading line ', iline, &
                        ' in file WALL_TEMPERATURE.'
                    CALL message( 'usm_read_wall_temperature', 'PA0522', 1, 2, 0, 6, 0 )
                ENDDO
 14             CLOSE(152)
                CYCLE
 15             message_string = 'file WALL_TEMPERATURE'//TRIM(coupling_char)//' does not exist'
                CALL message( 'usm_read_wall_temperature', 'PA0523', 1, 2, 0, 6, 0 )
            ENDIF
#if defined( __parallel )
            CALL MPI_BARRIER( comm2d, ierr )
#endif
        ENDDO

        IF ( debug_output )  CALL debug_message( 'usm_read_wall_temperature', 'end' )

    END SUBROUTINE usm_read_wall_temperature



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Solver for the energy balance at the ground/roof/wall surface.
!> It follows basic ideas and structure of lsm_energy_balance
!> with many simplifications and adjustments.
!> TODO better description
!> No calculation of window surface temperatures during spinup to increase
!> maximum possible timstep
!------------------------------------------------------------------------------!
    SUBROUTINE usm_surface_energy_balance( during_spinup )

        USE exchange_horiz_mod,                                                    &
            ONLY:  exchange_horiz


        IMPLICIT NONE

        INTEGER(iwp)                          :: i, j, k, l, m   !< running indices
        
        INTEGER(iwp) ::  i_off     !< offset to determine index of surface element, seen from atmospheric grid point, for x
        INTEGER(iwp) ::  j_off     !< offset to determine index of surface element, seen from atmospheric grid point, for y
        INTEGER(iwp) ::  k_off     !< offset to determine index of surface element, seen from atmospheric grid point, for z

        LOGICAL                               :: during_spinup      !< flag indicating soil/wall spinup phase
        
        REAL(wp)                              :: frac_win           !< window fraction, used to restore original values during spinup
        REAL(wp)                              :: frac_green         !< green fraction, used to restore original values during spinup
        REAL(wp)                              :: frac_wall          !< wall fraction, used to restore original values during spinup
        REAL(wp)                              :: stend_wall         !< surface tendency
        
        REAL(wp)                              :: stend_window       !< surface tendency
        REAL(wp)                              :: stend_green        !< surface tendency
        REAL(wp)                              :: coef_1             !< first coeficient for prognostic equation
        REAL(wp)                              :: coef_window_1      !< first coeficient for prognostic window equation
        REAL(wp)                              :: coef_green_1       !< first coeficient for prognostic green wall equation
        REAL(wp)                              :: coef_2             !< second  coeficient for prognostic equation
        REAL(wp)                              :: coef_window_2      !< second  coeficient for prognostic window equation
        REAL(wp)                              :: coef_green_2       !< second  coeficient for prognostic green wall equation
        REAL(wp)                              :: rho_cp             !< rho_wall_surface * c_p
        REAL(wp)                              :: f_shf              !< factor for shf_eb
        REAL(wp)                              :: f_shf_window       !< factor for shf_eb window
        REAL(wp)                              :: f_shf_green        !< factor for shf_eb green wall
        REAL(wp)                              :: lambda_surface     !< current value of lambda_surface (heat conductivity
                                                                    !<between air and wall)
        REAL(wp)                              :: lambda_surface_window  !< current value of lambda_surface (heat conductivity
                                                                        !< between air and window)
        REAL(wp)                              :: lambda_surface_green   !< current value of lambda_surface (heat conductivity
                                                                        !< between air and greeb wall)
        
        REAL(wp)                              :: dtime              !< simulated time of day (in UTC)
        INTEGER(iwp)                          :: dhour              !< simulated hour of day (in UTC)
        REAL(wp)                              :: acoef              !< actual coefficient of diurnal profile of anthropogenic heat
        REAL(wp) ::  f1,          &  !< resistance correction term 1
                     f2,          &  !< resistance correction term 2
                     f3,          &  !< resistance correction term 3
                     e,           &  !< water vapour pressure
                     e_s,         &  !< water vapour saturation pressure
                     e_s_dt,      &  !< derivate of e_s with respect to T
                     tend,        &  !< tendency
                     dq_s_dt,     &  !< derivate of q_s with respect to T
                     f_qsws,      &  !< factor for qsws
                     f_qsws_veg,  &  !< factor for qsws_veg
                     f_qsws_liq,  &  !< factor for qsws_liq
                     m_liq_max,   &  !< maxmimum value of the liq. water reservoir
                     qv1,         &  !< specific humidity at first grid level
                     m_max_depth = 0.0002_wp, &  !< Maximum capacity of the water reservoir (m)
                     rho_lv,      &  !< frequently used parameter for green layers
                     drho_l_lv,   &  !< frequently used parameter for green layers
                     q_s             !< saturation specific humidity


        IF ( debug_output_timestep )  THEN
           WRITE( debug_string, * ) 'usm_surface_energy_balance | during_spinup: ',&
                                    during_spinup
           CALL debug_message( debug_string, 'start' )
        ENDIF
!
!--     Index offset of surface element point with respect to adjoining 
!--     atmospheric grid point 
        k_off = surf_usm_h%koff
        j_off = surf_usm_h%joff
        i_off = surf_usm_h%ioff
        
!        
!--     First, treat horizontal surface elements
        !$OMP PARALLEL PRIVATE (m, i, j, k, lambda_surface, lambda_surface_window,                 &
        !$OMP&                  lambda_surface_green, qv1, rho_cp, rho_lv, drho_l_lv, f_shf,       &
        !$OMP&                  f_shf_window, f_shf_green, m_total, f1, f2, e_s, e, f3, f_qsws_veg,&
        !$OMP&                  q_s, f_qsws_liq, f_qsws, e_s_dt, dq_s_dt, coef_1, coef_window_1,   &
        !$OMP&                  coef_green_1, coef_2, coef_window_2, coef_green_2, stend_wall,     &
        !$OMP&                  stend_window, stend_green, tend, m_liq_max)
        !$OMP DO SCHEDULE (STATIC)
        DO  m = 1, surf_usm_h%ns
!
!--       During spinup set green and window fraction to zero and restore
!--       at the end of the loop.
!--       Note, this is a temporary fix and need to be removed later.  
           IF ( during_spinup )  THEN
              frac_win   = surf_usm_h%frac(m,ind_wat_win)
              frac_wall  = surf_usm_h%frac(m,ind_veg_wall)
              frac_green = surf_usm_h%frac(m,ind_pav_green)
              surf_usm_h%frac(m,ind_wat_win)   = 0.0_wp
              surf_usm_h%frac(m,ind_veg_wall)  = 1.0_wp
              surf_usm_h%frac(m,ind_pav_green) = 0.0_wp
           ENDIF
!
!--        Get indices of respective grid point
           i = surf_usm_h%i(m)
           j = surf_usm_h%j(m)
           k = surf_usm_h%k(m)
!
!--        TODO - how to calculate lambda_surface for horizontal surfaces
!--        (lambda_surface is set according to stratification in land surface model)
!--        MS: ???
           IF ( surf_usm_h%ol(m) >= 0.0_wp )  THEN
              lambda_surface = surf_usm_h%lambda_surf(m)
              lambda_surface_window = surf_usm_h%lambda_surf_window(m)
              lambda_surface_green = surf_usm_h%lambda_surf_green(m)
           ELSE
              lambda_surface = surf_usm_h%lambda_surf(m)
              lambda_surface_window = surf_usm_h%lambda_surf_window(m)
              lambda_surface_green = surf_usm_h%lambda_surf_green(m)
           ENDIF

!            pt1  = pt(k,j,i)
           IF ( humidity )  THEN
              qv1 = q(k,j,i)
           ELSE
              qv1 = 0.0_wp
           ENDIF
!
!--        calculate rho * c_p coefficient at surface layer
           rho_cp  = c_p * hyp(k) / ( r_d * surf_usm_h%pt1(m) * exner(k) )

           IF ( surf_usm_h%frac(m,ind_pav_green) > 0.0_wp )  THEN
!
!--           Calculate frequently used parameters
              rho_lv    = rho_cp / c_p * l_v
              drho_l_lv = 1.0_wp / (rho_l * l_v)
           ENDIF

!
!--        Calculate aerodyamic resistance. 
!--        Calculation for horizontal surfaces follows LSM formulation
!--        pt, us, ts are not available for the prognostic time step,
!--        data from the last time step is used here.
!
!--        Workaround: use single r_a as stability is only treated for the
!--        average temperature
           surf_usm_h%r_a(m) = ( surf_usm_h%pt1(m) - surf_usm_h%pt_surface(m) ) /&
                               ( surf_usm_h%ts(m) * surf_usm_h%us(m) + 1.0E-20_wp )   
           surf_usm_h%r_a_window(m) = surf_usm_h%r_a(m)
           surf_usm_h%r_a_green(m)  = surf_usm_h%r_a(m)

!            r_a = ( surf_usm_h%pt1(m) - t_surf_h(m) / exner(k) ) /                              &
!                  ( surf_usm_h%ts(m) * surf_usm_h%us(m) + 1.0E-20_wp )
!            r_a_window = ( surf_usm_h%pt1(m) - t_surf_window_h(m) / exner(k) ) /                &
!                  ( surf_usm_h%ts(m) * surf_usm_h%us(m) + 1.0E-20_wp )
!            r_a_green = ( surf_usm_h%pt1(m) - t_surf_green_h(m) / exner(k) ) /                  &
!                  ( surf_usm_h%ts(m) * surf_usm_h%us(m) + 1.0E-20_wp )
                
!--        Make sure that the resistance does not drop to zero
           IF ( surf_usm_h%r_a(m)        < 1.0_wp )                            &
               surf_usm_h%r_a(m)        = 1.0_wp
           IF ( surf_usm_h%r_a_green(m)  < 1.0_wp )                            &
               surf_usm_h%r_a_green(m)  = 1.0_wp
           IF ( surf_usm_h%r_a_window(m) < 1.0_wp )                            &
               surf_usm_h%r_a_window(m) = 1.0_wp
              
!
!--        Make sure that the resistacne does not exceed a maxmium value in case
!--        of zero velocities
           IF ( surf_usm_h%r_a(m)        > 300.0_wp )                          &
               surf_usm_h%r_a(m)        = 300.0_wp
           IF ( surf_usm_h%r_a_green(m)  > 300.0_wp )                          &
               surf_usm_h%r_a_green(m)  = 300.0_wp
           IF ( surf_usm_h%r_a_window(m) > 300.0_wp )                          &
               surf_usm_h%r_a_window(m) = 300.0_wp                
                
!
!--        factor for shf_eb
           f_shf  = rho_cp / surf_usm_h%r_a(m)
           f_shf_window  = rho_cp / surf_usm_h%r_a_window(m)
           f_shf_green  = rho_cp / surf_usm_h%r_a_green(m)
        

           IF ( surf_usm_h%frac(m,ind_pav_green) > 0.0_wp ) THEN
!--           Adapted from LSM:
!--           Second step: calculate canopy resistance r_canopy
!--           f1-f3 here are defined as 1/f1-f3 as in ECMWF documentation
 
!--           f1: correction for incoming shortwave radiation (stomata close at 
!--           night)
              f1 = MIN( 1.0_wp, ( 0.004_wp * surf_usm_h%rad_sw_in(m) + 0.05_wp ) / &
                               (0.81_wp * (0.004_wp * surf_usm_h%rad_sw_in(m)      &
                                + 1.0_wp)) )
!
!--           f2: correction for soil moisture availability to plants (the 
!--           integrated soil moisture must thus be considered here)
!--           f2 = 0 for very dry soils
              m_total = 0.0_wp
              DO  k = nzb_wall, nzt_wall+1
                  m_total = m_total + rootfr_h(nzb_wall,m)                              &
                            * MAX(swc_h(nzb_wall,m),wilt_h(nzb_wall,m))
              ENDDO 
    
              IF ( m_total > wilt_h(nzb_wall,m)  .AND.  m_total < fc_h(nzb_wall,m) )  THEN
                 f2 = ( m_total - wilt_h(nzb_wall,m) ) / (fc_h(nzb_wall,m) - wilt_h(nzb_wall,m) )
              ELSEIF ( m_total >= fc_h(nzb_wall,m) )  THEN
                 f2 = 1.0_wp
              ELSE
                 f2 = 1.0E-20_wp
              ENDIF
       
!
!--          Calculate water vapour pressure at saturation 
              e_s = 0.01_wp * 610.78_wp * EXP( 17.269_wp * ( t_surf_green_h(m) &
                            - 273.16_wp ) / ( t_surf_green_h(m) - 35.86_wp ) )
!
!--           f3: correction for vapour pressure deficit
              IF ( surf_usm_h%g_d(m) /= 0.0_wp )  THEN
!
!--           Calculate vapour pressure
                 e  = qv1 * surface_pressure / ( qv1 + 0.622_wp )
                 f3 = EXP ( - surf_usm_h%g_d(m) * (e_s - e) )
              ELSE
                 f3 = 1.0_wp
              ENDIF

!
!--           Calculate canopy resistance. In case that c_veg is 0 (bare soils),
!--           this calculation is obsolete, as r_canopy is not used below.
!--           To do: check for very dry soil -> r_canopy goes to infinity
              surf_usm_h%r_canopy(m) = surf_usm_h%r_canopy_min(m) /                   &
                              ( surf_usm_h%lai(m) * f1 * f2 * f3 + 1.0E-20_wp )

!
!--           Calculate the maximum possible liquid water amount on plants and
!--           bare surface. For vegetated surfaces, a maximum depth of 0.2 mm is
!--           assumed, while paved surfaces might hold up 1 mm of water. The 
!--           liquid water fraction for paved surfaces is calculated after 
!--           Noilhan & Planton (1989), while the ECMWF formulation is used for
!--           vegetated surfaces and bare soils.
              m_liq_max = m_max_depth * ( surf_usm_h%lai(m) )

              surf_usm_h%c_liq(m) = MIN( 1.0_wp, ( m_liq_usm_h%var_usm_1d(m) / m_liq_max )**0.67 )
!
!--           Calculate saturation specific humidity
              q_s = 0.622_wp * e_s / ( surface_pressure - e_s )
!
!--           In case of dewfall, set evapotranspiration to zero
!--           All super-saturated water is then removed from the air
              IF ( humidity  .AND.  q_s <= qv1 )  THEN
                 surf_usm_h%r_canopy(m) = 0.0_wp
              ENDIF

!
!--           Calculate coefficients for the total evapotranspiration 
!--           In case of water surface, set vegetation and soil fluxes to zero.
!--           For pavements, only evaporation of liquid water is possible.
              f_qsws_veg  = rho_lv *                                           &
                                ( 1.0_wp        - surf_usm_h%c_liq(m)    ) /   &
                                ( surf_usm_h%r_a_green(m) + surf_usm_h%r_canopy(m) )
              f_qsws_liq  = rho_lv * surf_usm_h%c_liq(m)   /                   &
                                  surf_usm_h%r_a_green(m)
       
              f_qsws = f_qsws_veg + f_qsws_liq
!
!--           Calculate derivative of q_s for Taylor series expansion
              e_s_dt = e_s * ( 17.269_wp / ( t_surf_green_h(m) - 35.86_wp) -   &
                               17.269_wp*( t_surf_green_h(m) - 273.16_wp)      &
                              / ( t_surf_green_h(m) - 35.86_wp)**2 )
       
              dq_s_dt = 0.622_wp * e_s_dt / ( surface_pressure - e_s_dt )
           ENDIF
!
!--        add LW up so that it can be removed in prognostic equation
           surf_usm_h%rad_net_l(m) = surf_usm_h%rad_sw_in(m)  -                &
                                     surf_usm_h%rad_sw_out(m) +                &
                                     surf_usm_h%rad_lw_in(m)  -                &
                                     surf_usm_h%rad_lw_out(m)
!
!--     numerator of the prognostic equation
!--     Todo: Adjust to tile approach. So far, emissivity for wall (element 0)
!--     is used
           coef_1 = surf_usm_h%rad_net_l(m) +                                  & 
                 ( 3.0_wp + 1.0_wp ) * surf_usm_h%emissivity(m,ind_veg_wall) * &
                                       sigma_sb * t_surf_wall_h(m) ** 4 +      &  
                                       f_shf * surf_usm_h%pt1(m) +             &
                                       lambda_surface * t_wall_h(nzb_wall,m)
           IF ( ( .NOT. during_spinup ) .AND. (surf_usm_h%frac(m,ind_wat_win) > 0.0_wp ) ) THEN
              coef_window_1 = surf_usm_h%rad_net_l(m) +                           & 
                      ( 3.0_wp + 1.0_wp ) * surf_usm_h%emissivity(m,ind_wat_win)  &
                                          * sigma_sb * t_surf_window_h(m) ** 4 +  &  
                                          f_shf_window * surf_usm_h%pt1(m) +      &
                                          lambda_surface_window * t_window_h(nzb_wall,m)
           ENDIF                 
           IF ( ( humidity ) .AND. ( surf_usm_h%frac(m,ind_pav_green) > 0.0_wp ) )  THEN
                    coef_green_1 = surf_usm_h%rad_net_l(m) +                                 & 
                   ( 3.0_wp + 1.0_wp ) * surf_usm_h%emissivity(m,ind_pav_green) * sigma_sb * &
                                       t_surf_green_h(m) ** 4 +                  &  
                                          f_shf_green * surf_usm_h%pt1(m) + f_qsws * ( qv1 - q_s    &
                                          + dq_s_dt * t_surf_green_h(m) )        &
                                          +lambda_surface_green * t_green_h(nzb_wall,m)
           ELSE
           coef_green_1 = surf_usm_h%rad_net_l(m) +                            & 
                 ( 3.0_wp + 1.0_wp ) * surf_usm_h%emissivity(m,ind_pav_green) *&
                                       sigma_sb * t_surf_green_h(m) ** 4 +     &  
                                       f_shf_green * surf_usm_h%pt1(m) +       &
                                       lambda_surface_green * t_green_h(nzb_wall,m)
          ENDIF
!
!--        denominator of the prognostic equation
           coef_2 = 4.0_wp * surf_usm_h%emissivity(m,ind_veg_wall) *           &
                             sigma_sb * t_surf_wall_h(m) ** 3                  &
                           + lambda_surface + f_shf / exner(k)
           IF ( ( .NOT. during_spinup ) .AND. ( surf_usm_h%frac(m,ind_wat_win) > 0.0_wp ) ) THEN
              coef_window_2 = 4.0_wp * surf_usm_h%emissivity(m,ind_wat_win) *     &
                                sigma_sb * t_surf_window_h(m) ** 3                &
                              + lambda_surface_window + f_shf_window / exner(k)
           ENDIF
           IF ( ( humidity ) .AND. ( surf_usm_h%frac(m,ind_pav_green) > 0.0_wp ) )  THEN
              coef_green_2 = 4.0_wp * surf_usm_h%emissivity(m,ind_pav_green) * sigma_sb *    &
                                t_surf_green_h(m) ** 3 + f_qsws * dq_s_dt                    &
                              + lambda_surface_green + f_shf_green / exner(k)
           ELSE
           coef_green_2 = 4.0_wp * surf_usm_h%emissivity(m,ind_pav_green) * sigma_sb *    &
                             t_surf_green_h(m) ** 3                                       &
                           + lambda_surface_green + f_shf_green / exner(k)
           ENDIF
!
!--        implicit solution when the surface layer has no heat capacity,
!--        otherwise use RK3 scheme.
           t_surf_wall_h_p(m) = ( coef_1 * dt_3d * tsc(2) +                        &
                             surf_usm_h%c_surface(m) * t_surf_wall_h(m) ) /        & 
                           ( surf_usm_h%c_surface(m) + coef_2 * dt_3d * tsc(2) ) 
           IF (( .NOT. during_spinup ) .AND. (surf_usm_h%frac(m,ind_wat_win) > 0.0_wp)) THEN
              t_surf_window_h_p(m) = ( coef_window_1 * dt_3d * tsc(2) +                        &
                                surf_usm_h%c_surface_window(m) * t_surf_window_h(m) ) /        & 
                              ( surf_usm_h%c_surface_window(m) + coef_window_2 * dt_3d * tsc(2) )
           ENDIF
           t_surf_green_h_p(m) = ( coef_green_1 * dt_3d * tsc(2) +                        &
                             surf_usm_h%c_surface_green(m) * t_surf_green_h(m) ) /        & 
                           ( surf_usm_h%c_surface_green(m) + coef_green_2 * dt_3d * tsc(2) ) 
!
!--        add RK3 term
           t_surf_wall_h_p(m) = t_surf_wall_h_p(m) + dt_3d * tsc(3) *         &
                           surf_usm_h%tt_surface_wall_m(m)

           t_surf_window_h_p(m) = t_surf_window_h_p(m) + dt_3d * tsc(3) *     &
                           surf_usm_h%tt_surface_window_m(m)

           t_surf_green_h_p(m) = t_surf_green_h_p(m) + dt_3d * tsc(3) *       &
                           surf_usm_h%tt_surface_green_m(m)
!
!--        Store surface temperature on pt_surface. Further, in case humidity is used
!--        store also vpt_surface, which is, due to the lack of moisture on roofs simply
!--        assumed to be the surface temperature.
           surf_usm_h%pt_surface(m) = ( surf_usm_h%frac(m,ind_veg_wall) * t_surf_wall_h_p(m)   &
                               + surf_usm_h%frac(m,ind_wat_win) * t_surf_window_h_p(m)         &
                               + surf_usm_h%frac(m,ind_pav_green) * t_surf_green_h_p(m) )      &
                               / exner(k)
                               
           IF ( humidity )  surf_usm_h%vpt_surface(m) =                        &
                                                   surf_usm_h%pt_surface(m)
!
!--        calculate true tendency
           stend_wall = ( t_surf_wall_h_p(m) - t_surf_wall_h(m) - dt_3d * tsc(3) *              &
                     surf_usm_h%tt_surface_wall_m(m)) / ( dt_3d  * tsc(2) )
           stend_window = ( t_surf_window_h_p(m) - t_surf_window_h(m) - dt_3d * tsc(3) *        &
                     surf_usm_h%tt_surface_window_m(m)) / ( dt_3d  * tsc(2) )
           stend_green = ( t_surf_green_h_p(m) - t_surf_green_h(m) - dt_3d * tsc(3) *           &
                     surf_usm_h%tt_surface_green_m(m)) / ( dt_3d  * tsc(2) )
!
!--        calculate t_surf tendencies for the next Runge-Kutta step
           IF ( timestep_scheme(1:5) == 'runge' )  THEN
              IF ( intermediate_timestep_count == 1 )  THEN
                 surf_usm_h%tt_surface_wall_m(m) = stend_wall
                 surf_usm_h%tt_surface_window_m(m) = stend_window
                 surf_usm_h%tt_surface_green_m(m) = stend_green
              ELSEIF ( intermediate_timestep_count <                          &
                        intermediate_timestep_count_max )  THEN
                 surf_usm_h%tt_surface_wall_m(m) = -9.5625_wp * stend_wall +       &
                                     5.3125_wp * surf_usm_h%tt_surface_wall_m(m)
                 surf_usm_h%tt_surface_window_m(m) = -9.5625_wp * stend_window +   &
                                     5.3125_wp * surf_usm_h%tt_surface_window_m(m)
                 surf_usm_h%tt_surface_green_m(m) = -9.5625_wp * stend_green +     &
                                     5.3125_wp * surf_usm_h%tt_surface_green_m(m)
              ENDIF
           ENDIF
!
!--        in case of fast changes in the skin temperature, it is required to
!--        update the radiative fluxes in order to keep the solution stable
           IF ( ( ( ABS( t_surf_wall_h_p(m)   - t_surf_wall_h(m) )   > 1.0_wp )   .OR. &
                (   ABS( t_surf_green_h_p(m)  - t_surf_green_h(m) )  > 1.0_wp )   .OR. &
                (   ABS( t_surf_window_h_p(m) - t_surf_window_h(m) ) > 1.0_wp ) )      &
                   .AND.  unscheduled_radiation_calls  )  THEN
              force_radiation_call_l = .TRUE.
           ENDIF
!
!--        calculate fluxes
!--        rad_net_l is never used!
           surf_usm_h%rad_net_l(m) = surf_usm_h%rad_net_l(m) +                           &
                                     surf_usm_h%frac(m,ind_veg_wall) *                   &
                                     sigma_sb * surf_usm_h%emissivity(m,ind_veg_wall) *  &
                                     ( t_surf_wall_h_p(m)**4 - t_surf_wall_h(m)**4 )     &
                                    + surf_usm_h%frac(m,ind_wat_win) *                   &
                                     sigma_sb * surf_usm_h%emissivity(m,ind_wat_win) *   &
                                     ( t_surf_window_h_p(m)**4 - t_surf_window_h(m)**4 ) &
                                    + surf_usm_h%frac(m,ind_pav_green) *                 &
                                     sigma_sb * surf_usm_h%emissivity(m,ind_pav_green) * &
                                     ( t_surf_green_h_p(m)**4 - t_surf_green_h(m)**4 )

           surf_usm_h%wghf_eb(m)   = lambda_surface *                                    &
                                      ( t_surf_wall_h_p(m) - t_wall_h(nzb_wall,m) )
           surf_usm_h%wghf_eb_green(m)  = lambda_surface_green *                         &
                                          ( t_surf_green_h_p(m) - t_green_h(nzb_wall,m) )
           surf_usm_h%wghf_eb_window(m) = lambda_surface_window *                        &
                                           ( t_surf_window_h_p(m) - t_window_h(nzb_wall,m) )

!
!--        ground/wall/roof surface heat flux
           surf_usm_h%wshf_eb(m)   = - f_shf  * ( surf_usm_h%pt1(m) - t_surf_wall_h_p(m) / exner(k) ) *          &
                                       surf_usm_h%frac(m,ind_veg_wall)         &
                                     - f_shf_window  * ( surf_usm_h%pt1(m) - t_surf_window_h_p(m) / exner(k) ) * &
                                       surf_usm_h%frac(m,ind_wat_win)          &
                                     - f_shf_green  * ( surf_usm_h%pt1(m) - t_surf_green_h_p(m) / exner(k) ) *   &
                                       surf_usm_h%frac(m,ind_pav_green)
!           
!--        store kinematic surface heat fluxes for utilization in other processes
!--        diffusion_s, surface_layer_fluxes,...
           surf_usm_h%shf(m) = surf_usm_h%wshf_eb(m) / c_p
!
!--        If the indoor model is applied, further add waste heat from buildings to the 
!--        kinematic flux.
           IF ( indoor_model )  THEN
              surf_usm_h%shf(m) = surf_usm_h%shf(m) + surf_usm_h%waste_heat(m) / c_p
           ENDIF
      

           IF (surf_usm_h%frac(m,ind_pav_green) > 0.0_wp) THEN
              
           
              IF ( humidity )  THEN
                 surf_usm_h%qsws(m)  = - f_qsws * ( qv1 - q_s + dq_s_dt                     &
                                 * t_surf_green_h(m) - dq_s_dt *               &
                                   t_surf_green_h_p(m) )
       
                 surf_usm_h%qsws_veg(m)  = - f_qsws_veg  * ( qv1 - q_s                      &
                                     + dq_s_dt * t_surf_green_h(m) - dq_s_dt   &
                                     * t_surf_green_h_p(m) )
       
                 surf_usm_h%qsws_liq(m)  = - f_qsws_liq  * ( qv1 - q_s                      &
                                     + dq_s_dt * t_surf_green_h(m) - dq_s_dt   &
                                     * t_surf_green_h_p(m) )
                                     
              ENDIF
 
!
!--           Calculate the true surface resistance
              IF ( .NOT.  humidity )  THEN
                 surf_usm_h%r_s(m) = 1.0E10_wp
              ELSE
                 surf_usm_h%r_s(m) = - rho_lv * ( qv1 - q_s + dq_s_dt                       &
                                 *  t_surf_green_h(m) - dq_s_dt *              &
                                   t_surf_green_h_p(m) ) /                     &
                                   (surf_usm_h%qsws(m) + 1.0E-20)  - surf_usm_h%r_a_green(m)
              ENDIF
 
!
!--           Calculate change in liquid water reservoir due to dew fall or 
!--           evaporation of liquid water
              IF ( humidity )  THEN
!
!--              If precipitation is activated, add rain water to qsws_liq
!--              and qsws_soil according the the vegetation coverage.
!--              precipitation_rate is given in mm.
                 IF ( precipitation )  THEN

!
!--                 Add precipitation to liquid water reservoir, if possible.
!--                 Otherwise, add the water to soil. In case of
!--                 pavements, the exceeding water amount is implicitely removed 
!--                 as runoff as qsws_soil is then not used in the soil model
                    IF ( m_liq_usm_h%var_usm_1d(m) /= m_liq_max )  THEN
                       surf_usm_h%qsws_liq(m) = surf_usm_h%qsws_liq(m)                &
                                        + surf_usm_h%frac(m,ind_pav_green) * prr(k+k_off,j+j_off,i+i_off)&
                                        * hyrho(k+k_off)                              &
                                        * 0.001_wp * rho_l * l_v
                   ENDIF

                 ENDIF

!
!--              If the air is saturated, check the reservoir water level
                 IF ( surf_usm_h%qsws(m) < 0.0_wp )  THEN
!
!--                 Check if reservoir is full (avoid values > m_liq_max)
!--                 In that case, qsws_liq goes to qsws_soil. In this 
!--                 case qsws_veg is zero anyway (because c_liq = 1),       
!--                 so that tend is zero and no further check is needed
                    IF ( m_liq_usm_h%var_usm_1d(m) == m_liq_max )  THEN
!                      surf_usm_h%qsws_soil(m) = surf_usm_h%qsws_soil(m) + surf_usm_h%qsws_liq(m)
                       surf_usm_h%qsws_liq(m)  = 0.0_wp
                    ENDIF

!
!--                 In case qsws_veg becomes negative (unphysical behavior), 
!--                 let the water enter the liquid water reservoir as dew on the
!--                 plant
                    IF ( surf_usm_h%qsws_veg(m) < 0.0_wp )  THEN
                       surf_usm_h%qsws_liq(m) = surf_usm_h%qsws_liq(m) + surf_usm_h%qsws_veg(m)
                       surf_usm_h%qsws_veg(m) = 0.0_wp
                    ENDIF
                 ENDIF                    
  
                 surf_usm_h%qsws(m) = surf_usm_h%qsws(m) / l_v
        
                 tend = - surf_usm_h%qsws_liq(m) * drho_l_lv
                 m_liq_usm_h_p%var_usm_1d(m) = m_liq_usm_h%var_usm_1d(m) + dt_3d *    &
                                               ( tsc(2) * tend +                      &
                                                 tsc(3) * tm_liq_usm_h_m%var_usm_1d(m) )
!
!--             Check if reservoir is overfull -> reduce to maximum
!--             (conservation of water is violated here)
                 m_liq_usm_h_p%var_usm_1d(m) = MIN( m_liq_usm_h_p%var_usm_1d(m),m_liq_max )
 
!
!--             Check if reservoir is empty (avoid values < 0.0)
!--             (conservation of water is violated here)
                 m_liq_usm_h_p%var_usm_1d(m) = MAX( m_liq_usm_h_p%var_usm_1d(m), 0.0_wp )
!
!--             Calculate m_liq tendencies for the next Runge-Kutta step
                 IF ( timestep_scheme(1:5) == 'runge' )  THEN
                    IF ( intermediate_timestep_count == 1 )  THEN
                       tm_liq_usm_h_m%var_usm_1d(m) = tend
                    ELSEIF ( intermediate_timestep_count <                            &
                             intermediate_timestep_count_max )  THEN
                       tm_liq_usm_h_m%var_usm_1d(m) = -9.5625_wp * tend +             &
                                                     5.3125_wp * tm_liq_usm_h_m%var_usm_1d(m)
                    ENDIF
                 ENDIF
 
              ENDIF
           ELSE
              surf_usm_h%r_s(m) = 1.0E10_wp
           ENDIF
!
!--        During spinup green and window fraction are set to zero. Here, the original
!--        values are restored.
           IF ( during_spinup )  THEN
              surf_usm_h%frac(m,ind_wat_win)   = frac_win
              surf_usm_h%frac(m,ind_veg_wall)  = frac_wall
              surf_usm_h%frac(m,ind_pav_green) = frac_green
           ENDIF
 
       ENDDO
!
!--    Now, treat vertical surface elements
       !$OMP DO SCHEDULE (STATIC)
       DO  l = 0, 3
           DO  m = 1, surf_usm_v(l)%ns
!
!--           During spinup set green and window fraction to zero and restore
!--           at the end of the loop.
!--           Note, this is a temporary fix and need to be removed later. 
              IF ( during_spinup )  THEN
                 frac_win   = surf_usm_v(l)%frac(m,ind_wat_win)
                 frac_wall  = surf_usm_v(l)%frac(m,ind_veg_wall)
                 frac_green = surf_usm_v(l)%frac(m,ind_pav_green)
                 surf_usm_v(l)%frac(m,ind_wat_win)   = 0.0_wp
                 surf_usm_v(l)%frac(m,ind_veg_wall)  = 1.0_wp
                 surf_usm_v(l)%frac(m,ind_pav_green) = 0.0_wp
              ENDIF
!
!--          Get indices of respective grid point
              i = surf_usm_v(l)%i(m)
              j = surf_usm_v(l)%j(m)
              k = surf_usm_v(l)%k(m)
 
!
!--          TODO - how to calculate lambda_surface for horizontal (??? do you mean verical ???) surfaces
!--          (lambda_surface is set according to stratification in land surface model).
!--          Please note, for vertical surfaces no ol is defined, since 
!--          stratification is not considered in this case.
              lambda_surface = surf_usm_v(l)%lambda_surf(m)
              lambda_surface_window = surf_usm_v(l)%lambda_surf_window(m)
              lambda_surface_green = surf_usm_v(l)%lambda_surf_green(m)
 
!            pt1  = pt(k,j,i)
              IF ( humidity )  THEN
                 qv1 = q(k,j,i)
              ELSE
                 qv1 = 0.0_wp
              ENDIF
!
!--           calculate rho * c_p coefficient at wall layer
              rho_cp  = c_p * hyp(k) / ( r_d * surf_usm_v(l)%pt1(m) * exner(k) )
              
              IF (surf_usm_v(l)%frac(m,ind_pav_green) > 0.0_wp )  THEN
!
!--              Calculate frequently used parameters
                 rho_lv    = rho_cp / c_p * l_v
                 drho_l_lv = 1.0_wp / (rho_l * l_v)
              ENDIF
 
!--          Calculation of r_a for vertical surfaces
!--
!--          heat transfer coefficient for forced convection along vertical walls
!--          follows formulation in TUF3d model (Krayenhoff & Voogt, 2006)
!--            
!--          H = httc (Tsfc - Tair)
!--          httc = rw * (11.8 + 4.2 * Ueff) - 4.0
!--            
!--                rw: wall patch roughness relative to 1.0 for concrete
!--                Ueff: effective wind speed
!--                - 4.0 is a reduction of Rowley et al (1930) formulation based on
!--                Cole and Sturrock (1977)
!--           
!--                Ucan: Canyon wind speed
!--                wstar: convective velocity
!--                Qs: surface heat flux
!--                zH: height of the convective layer
!--                wstar = (g/Tcan*Qs*zH)**(1./3.)
!--          Effective velocity components must always 
!--          be defined at scalar grid point. The wall normal component is 
!--          obtained by simple linear interpolation. ( An alternative would
!--          be an logarithmic interpolation. )
!--          Parameter roughness_concrete (default value = 0.001) is used
!--          to calculation of roughness relative to concrete
              surf_usm_v(l)%r_a(m) = rho_cp / ( surf_usm_v(l)%z0(m) /           &
                         roughness_concrete * ( 11.8_wp + 4.2_wp *              &
                         SQRT( MAX( ( ( u(k,j,i) + u(k,j,i+1) ) * 0.5_wp )**2 + &
                                    ( ( v(k,j,i) + v(k,j+1,i) ) * 0.5_wp )**2 + &
                                    ( ( w(k,j,i) + w(k-1,j,i) ) * 0.5_wp )**2,  &
                               0.01_wp ) )                                      &
                            )  - 4.0_wp  ) 
!
!--          Limit aerodynamic resistance
              IF ( surf_usm_v(l)%r_a(m) < 1.0_wp )  surf_usm_v(l)%r_a(m) = 1.0_wp   
              
                            
              f_shf         = rho_cp / surf_usm_v(l)%r_a(m)
              f_shf_window  = rho_cp / surf_usm_v(l)%r_a(m)
              f_shf_green   = rho_cp / surf_usm_v(l)%r_a(m)
 

              IF ( surf_usm_v(l)%frac(m,ind_pav_green) > 0.0_wp ) THEN
!
!--             Adapted from LSM:
!--             Second step: calculate canopy resistance r_canopy
!--             f1-f3 here are defined as 1/f1-f3 as in ECMWF documentation
!--             f1: correction for incoming shortwave radiation (stomata close at 
!--             night)
                 f1 = MIN( 1.0_wp, ( 0.004_wp * surf_usm_v(l)%rad_sw_in(m) + 0.05_wp ) / &
                                  (0.81_wp * (0.004_wp * surf_usm_v(l)%rad_sw_in(m)      &
                                   + 1.0_wp)) )
!
!--             f2: correction for soil moisture availability to plants (the 
!--             integrated soil moisture must thus be considered here)
!--             f2 = 0 for very dry soils
 
                 f2=1.0_wp
 
!
!--              Calculate water vapour pressure at saturation 
                 e_s = 0.01_wp * 610.78_wp * EXP( 17.269_wp * (  t_surf_green_v_p(l)%t(m) &
                               - 273.16_wp ) / (  t_surf_green_v_p(l)%t(m) - 35.86_wp ) )
!
!--              f3: correction for vapour pressure deficit
                 IF ( surf_usm_v(l)%g_d(m) /= 0.0_wp )  THEN
!
!--                 Calculate vapour pressure
                    e  = qv1 * surface_pressure / ( qv1 + 0.622_wp )
                    f3 = EXP ( - surf_usm_v(l)%g_d(m) * (e_s - e) )
                 ELSE
                    f3 = 1.0_wp
                 ENDIF
!
!--              Calculate canopy resistance. In case that c_veg is 0 (bare soils),
!--              this calculation is obsolete, as r_canopy is not used below.
!--              To do: check for very dry soil -> r_canopy goes to infinity
                 surf_usm_v(l)%r_canopy(m) = surf_usm_v(l)%r_canopy_min(m) /                  &
                                        ( surf_usm_v(l)%lai(m) * f1 * f2 * f3 + 1.0E-20_wp )
                               
!
!--              Calculate saturation specific humidity
                 q_s = 0.622_wp * e_s / ( surface_pressure - e_s )
!
!--              In case of dewfall, set evapotranspiration to zero
!--              All super-saturated water is then removed from the air
                 IF ( humidity  .AND.  q_s <= qv1 )  THEN
                    surf_usm_v(l)%r_canopy(m) = 0.0_wp
                 ENDIF
 
!
!--              Calculate coefficients for the total evapotranspiration 
!--              In case of water surface, set vegetation and soil fluxes to zero.
!--              For pavements, only evaporation of liquid water is possible.
                 f_qsws_veg  = rho_lv *                                &
                                   ( 1.0_wp        - 0.0_wp ) / & !surf_usm_h%c_liq(m)    ) /   &
                                   ( surf_usm_v(l)%r_a(m) + surf_usm_v(l)%r_canopy(m) )
!                f_qsws_liq  = rho_lv * surf_usm_h%c_liq(m)   /             &
!                              surf_usm_h%r_a_green(m)
          
                 f_qsws = f_qsws_veg! + f_qsws_liq
!
!--              Calculate derivative of q_s for Taylor series expansion
                 e_s_dt = e_s * ( 17.269_wp / ( t_surf_green_v_p(l)%t(m) - 35.86_wp) -   &
                                  17.269_wp*( t_surf_green_v_p(l)%t(m) - 273.16_wp)      &
                                 / ( t_surf_green_v_p(l)%t(m) - 35.86_wp)**2 )
          
                 dq_s_dt = 0.622_wp * e_s_dt / ( surface_pressure - e_s_dt )
              ENDIF

!
!--           add LW up so that it can be removed in prognostic equation
              surf_usm_v(l)%rad_net_l(m) = surf_usm_v(l)%rad_sw_in(m)  -        &
                                           surf_usm_v(l)%rad_sw_out(m) +        &
                                           surf_usm_v(l)%rad_lw_in(m)  -        &
                                           surf_usm_v(l)%rad_lw_out(m)
!
!--           numerator of the prognostic equation
              coef_1 = surf_usm_v(l)%rad_net_l(m) +                             & ! coef +1 corresponds to -lwout
                                                                                  ! included in calculation of radnet_l
              ( 3.0_wp + 1.0_wp ) * surf_usm_v(l)%emissivity(m,ind_veg_wall) *  &
                                      sigma_sb *  t_surf_wall_v(l)%t(m) ** 4 +  &  
                                      f_shf * surf_usm_v(l)%pt1(m) +            &
                                      lambda_surface * t_wall_v(l)%t(nzb_wall,m)
              IF ( ( .NOT. during_spinup ) .AND. ( surf_usm_v(l)%frac(m,ind_wat_win) > 0.0_wp ) ) THEN
                 coef_window_1 = surf_usm_v(l)%rad_net_l(m) +                   & ! coef +1 corresponds to -lwout
                                                                                  ! included in calculation of radnet_l
                ( 3.0_wp + 1.0_wp ) * surf_usm_v(l)%emissivity(m,ind_wat_win) * &
                                      sigma_sb * t_surf_window_v(l)%t(m) ** 4 + &  
                                      f_shf * surf_usm_v(l)%pt1(m) +            &
                                      lambda_surface_window * t_window_v(l)%t(nzb_wall,m)
              ENDIF
              IF ( ( humidity ) .AND. ( surf_usm_v(l)%frac(m,ind_pav_green) > 0.0_wp ) )  THEN
                 coef_green_1 = surf_usm_v(l)%rad_net_l(m) +                      & ! coef +1 corresponds to -lwout
                                                                                    ! included in calculation of radnet_l
                 ( 3.0_wp + 1.0_wp ) * surf_usm_v(l)%emissivity(m,ind_pav_green) * sigma_sb *  &
                                      t_surf_green_v(l)%t(m) ** 4 +               &  
                                      f_shf * surf_usm_v(l)%pt1(m) +     f_qsws * ( qv1 - q_s  &
                                           + dq_s_dt * t_surf_green_v(l)%t(m) ) +              &
                                      lambda_surface_green * t_wall_v(l)%t(nzb_wall,m)
              ELSE
                coef_green_1 = surf_usm_v(l)%rad_net_l(m) +                       & ! coef +1 corresponds to -lwout included
                                                                                    ! in calculation of radnet_l
                ( 3.0_wp + 1.0_wp ) * surf_usm_v(l)%emissivity(m,ind_pav_green) * sigma_sb *  &
                                      t_surf_green_v(l)%t(m) ** 4 +               &  
                                      f_shf * surf_usm_v(l)%pt1(m) +              &
                                      lambda_surface_green * t_wall_v(l)%t(nzb_wall,m)
              ENDIF
                                      
!
!--           denominator of the prognostic equation
              coef_2 = 4.0_wp * surf_usm_v(l)%emissivity(m,ind_veg_wall) * sigma_sb *   &
                                 t_surf_wall_v(l)%t(m) ** 3                             &
                               + lambda_surface + f_shf / exner(k)  
              IF ( ( .NOT. during_spinup ) .AND. ( surf_usm_v(l)%frac(m,ind_wat_win) > 0.0_wp ) ) THEN             
                 coef_window_2 = 4.0_wp * surf_usm_v(l)%emissivity(m,ind_wat_win) * sigma_sb *       &
                                   t_surf_window_v(l)%t(m) ** 3                         &
                                 + lambda_surface_window + f_shf / exner(k)
              ENDIF
              IF ( ( humidity ) .AND. ( surf_usm_v(l)%frac(m,ind_pav_green) > 0.0_wp ) )  THEN
                  coef_green_2 = 4.0_wp * surf_usm_v(l)%emissivity(m,ind_pav_green) * sigma_sb *     &
                                   t_surf_green_v(l)%t(m) ** 3  + f_qsws * dq_s_dt      &
                                 + lambda_surface_green + f_shf / exner(k)
              ELSE
                 coef_green_2 = 4.0_wp * surf_usm_v(l)%emissivity(m,ind_pav_green) * sigma_sb *      &
                                   t_surf_green_v(l)%t(m) ** 3                          &
                                 + lambda_surface_green + f_shf / exner(k)
              ENDIF
!
!--           implicit solution when the surface layer has no heat capacity,
!--           otherwise use RK3 scheme.
              t_surf_wall_v_p(l)%t(m) = ( coef_1 * dt_3d * tsc(2) +                 &
                             surf_usm_v(l)%c_surface(m) * t_surf_wall_v(l)%t(m) ) / & 
                             ( surf_usm_v(l)%c_surface(m) + coef_2 * dt_3d * tsc(2) ) 
              IF ( ( .NOT. during_spinup ) .AND. ( surf_usm_v(l)%frac(m,ind_wat_win) > 0.0_wp ) ) THEN
                 t_surf_window_v_p(l)%t(m) = ( coef_window_1 * dt_3d * tsc(2) +                 &
                                surf_usm_v(l)%c_surface_window(m) * t_surf_window_v(l)%t(m) ) / & 
                              ( surf_usm_v(l)%c_surface_window(m) + coef_window_2 * dt_3d * tsc(2) ) 
              ENDIF
              t_surf_green_v_p(l)%t(m) = ( coef_green_1 * dt_3d * tsc(2) +                 &
                             surf_usm_v(l)%c_surface_green(m) * t_surf_green_v(l)%t(m) ) / & 
                           ( surf_usm_v(l)%c_surface_green(m) + coef_green_2 * dt_3d * tsc(2) ) 
!
!--           add RK3 term
              t_surf_wall_v_p(l)%t(m) = t_surf_wall_v_p(l)%t(m) + dt_3d * tsc(3) *         &
                                surf_usm_v(l)%tt_surface_wall_m(m)
              t_surf_window_v_p(l)%t(m) = t_surf_window_v_p(l)%t(m) + dt_3d * tsc(3) *     &
                                surf_usm_v(l)%tt_surface_window_m(m)
              t_surf_green_v_p(l)%t(m) = t_surf_green_v_p(l)%t(m) + dt_3d * tsc(3) *       &
                                 surf_usm_v(l)%tt_surface_green_m(m)
!
!--           Store surface temperature. Further, in case humidity is used
!--           store also vpt_surface, which is, due to the lack of moisture on roofs simply
!--           assumed to be the surface temperature.     
              surf_usm_v(l)%pt_surface(m) =  ( surf_usm_v(l)%frac(m,ind_veg_wall) * t_surf_wall_v_p(l)%t(m)  &
                                      + surf_usm_v(l)%frac(m,ind_wat_win) * t_surf_window_v_p(l)%t(m)        &
                                      + surf_usm_v(l)%frac(m,ind_pav_green) * t_surf_green_v_p(l)%t(m) )     &
                                      / exner(k)
                                       
              IF ( humidity )  surf_usm_v(l)%vpt_surface(m) =                  &
                                                     surf_usm_v(l)%pt_surface(m)
!
!--           calculate true tendency
              stend_wall = ( t_surf_wall_v_p(l)%t(m) - t_surf_wall_v(l)%t(m) - dt_3d * tsc(3) *      &
                        surf_usm_v(l)%tt_surface_wall_m(m) ) / ( dt_3d  * tsc(2) )
              stend_window = ( t_surf_window_v_p(l)%t(m) - t_surf_window_v(l)%t(m) - dt_3d * tsc(3) *&
                        surf_usm_v(l)%tt_surface_window_m(m) ) / ( dt_3d  * tsc(2) )
              stend_green = ( t_surf_green_v_p(l)%t(m) - t_surf_green_v(l)%t(m) - dt_3d * tsc(3) *   &
                        surf_usm_v(l)%tt_surface_green_m(m) ) / ( dt_3d  * tsc(2) )

!
!--           calculate t_surf_* tendencies for the next Runge-Kutta step
              IF ( timestep_scheme(1:5) == 'runge' )  THEN
                 IF ( intermediate_timestep_count == 1 )  THEN
                    surf_usm_v(l)%tt_surface_wall_m(m) = stend_wall
                    surf_usm_v(l)%tt_surface_window_m(m) = stend_window
                    surf_usm_v(l)%tt_surface_green_m(m) = stend_green
                 ELSEIF ( intermediate_timestep_count <                                 &
                          intermediate_timestep_count_max )  THEN
                    surf_usm_v(l)%tt_surface_wall_m(m) = -9.5625_wp * stend_wall +      &
                                     5.3125_wp * surf_usm_v(l)%tt_surface_wall_m(m)
                    surf_usm_v(l)%tt_surface_green_m(m) = -9.5625_wp * stend_green +    &
                                     5.3125_wp * surf_usm_v(l)%tt_surface_green_m(m)
                    surf_usm_v(l)%tt_surface_window_m(m) = -9.5625_wp * stend_window +  &
                                     5.3125_wp * surf_usm_v(l)%tt_surface_window_m(m)
                 ENDIF
              ENDIF

!
!--           in case of fast changes in the skin temperature, it is required to
!--           update the radiative fluxes in order to keep the solution stable
 
              IF ( ( ( ABS( t_surf_wall_v_p(l)%t(m)   - t_surf_wall_v(l)%t(m) )   > 1.0_wp ) .OR. &
                   (   ABS( t_surf_green_v_p(l)%t(m)  - t_surf_green_v(l)%t(m) )  > 1.0_wp ) .OR. &
                   (   ABS( t_surf_window_v_p(l)%t(m) - t_surf_window_v(l)%t(m) ) > 1.0_wp ) )    &
                      .AND.  unscheduled_radiation_calls )  THEN
                 force_radiation_call_l = .TRUE.
              ENDIF

!
!--           calculate fluxes
!--           prognostic rad_net_l is used just for output!           
              surf_usm_v(l)%rad_net_l(m) = surf_usm_v(l)%frac(m,ind_veg_wall) *                      &
                                           ( surf_usm_v(l)%rad_net_l(m) +                            &
                                           3.0_wp * sigma_sb *                                       &
                                           t_surf_wall_v(l)%t(m)**4 - 4.0_wp * sigma_sb *            &
                                           t_surf_wall_v(l)%t(m)**3 * t_surf_wall_v_p(l)%t(m) )      &
                                         + surf_usm_v(l)%frac(m,ind_wat_win) *                       &
                                           ( surf_usm_v(l)%rad_net_l(m) +                            &
                                           3.0_wp * sigma_sb *                                       &
                                           t_surf_window_v(l)%t(m)**4 - 4.0_wp * sigma_sb *          &
                                           t_surf_window_v(l)%t(m)**3 * t_surf_window_v_p(l)%t(m) )  &
                                         + surf_usm_v(l)%frac(m,ind_pav_green) *                     &
                                           ( surf_usm_v(l)%rad_net_l(m) +                            &
                                           3.0_wp * sigma_sb *                                       &
                                           t_surf_green_v(l)%t(m)**4 - 4.0_wp * sigma_sb *           &
                                           t_surf_green_v(l)%t(m)**3 * t_surf_green_v_p(l)%t(m) )

              surf_usm_v(l)%wghf_eb_window(m) = lambda_surface_window * &
                                                ( t_surf_window_v_p(l)%t(m) - t_window_v(l)%t(nzb_wall,m) )
              surf_usm_v(l)%wghf_eb(m)   = lambda_surface *             &
                                                ( t_surf_wall_v_p(l)%t(m) - t_wall_v(l)%t(nzb_wall,m) )
              surf_usm_v(l)%wghf_eb_green(m)  = lambda_surface_green *  &
                                                ( t_surf_green_v_p(l)%t(m) - t_green_v(l)%t(nzb_wall,m) )

!
!--           ground/wall/roof surface heat flux
              surf_usm_v(l)%wshf_eb(m)   =                                     &
                 - f_shf  * ( surf_usm_v(l)%pt1(m) -                           &
                 t_surf_wall_v_p(l)%t(m) / exner(k) ) * surf_usm_v(l)%frac(m,ind_veg_wall)       &
                 - f_shf_window  * ( surf_usm_v(l)%pt1(m) -                    &
                 t_surf_window_v_p(l)%t(m) / exner(k) ) * surf_usm_v(l)%frac(m,ind_wat_win)&
                 - f_shf_green  * ( surf_usm_v(l)%pt1(m) -                     &
                 t_surf_green_v_p(l)%t(m) / exner(k) ) * surf_usm_v(l)%frac(m,ind_pav_green)

!           
!--           store kinematic surface heat fluxes for utilization in other processes
!--           diffusion_s, surface_layer_fluxes,...
              surf_usm_v(l)%shf(m) = surf_usm_v(l)%wshf_eb(m) / c_p
!
!--           If the indoor model is applied, further add waste heat from buildings to the 
!--           kinematic flux.
              IF ( indoor_model )  THEN
                 surf_usm_v(l)%shf(m) = surf_usm_v(l)%shf(m) +                       &
                                        surf_usm_v(l)%waste_heat(m) / c_p
              ENDIF              

              IF ( surf_usm_v(l)%frac(m,ind_pav_green) > 0.0_wp ) THEN
 

                 IF ( humidity )  THEN
                    surf_usm_v(l)%qsws(m)  = - f_qsws * ( qv1 - q_s + dq_s_dt          &
                                    * t_surf_green_v(l)%t(m) - dq_s_dt *               &
                                      t_surf_green_v_p(l)%t(m) )
          
                    surf_usm_v(l)%qsws(m) = surf_usm_v(l)%qsws(m) / l_v
          
                    surf_usm_v(l)%qsws_veg(m)  = - f_qsws_veg  * ( qv1 - q_s           &
                                        + dq_s_dt * t_surf_green_v(l)%t(m) - dq_s_dt   &
                                        * t_surf_green_v_p(l)%t(m) )
          
!                    surf_usm_h%qsws_liq(m)  = - f_qsws_liq  * ( qv1 - q_s         &
!                                        + dq_s_dt * t_surf_green_h(m) - dq_s_dt   &
!                                        * t_surf_green_h_p(m) )
                 ENDIF
 
!
!--              Calculate the true surface resistance
                 IF ( .NOT.  humidity )  THEN
                    surf_usm_v(l)%r_s(m) = 1.0E10_wp
                 ELSE
                    surf_usm_v(l)%r_s(m) = - rho_lv * ( qv1 - q_s + dq_s_dt             &
                                    *  t_surf_green_v(l)%t(m) - dq_s_dt *               &
                                      t_surf_green_v_p(l)%t(m) ) /                      &
                                      (surf_usm_v(l)%qsws(m) + 1.0E-20)  - surf_usm_v(l)%r_a(m)
                 ENDIF
         
!
!--              Calculate change in liquid water reservoir due to dew fall or 
!--              evaporation of liquid water
                 IF ( humidity )  THEN
!
!--                 If the air is saturated, check the reservoir water level
                    IF ( surf_usm_v(l)%qsws(m) < 0.0_wp )  THEN
       
!
!--                    In case qsws_veg becomes negative (unphysical behavior), 
!--                    let the water enter the liquid water reservoir as dew on the
!--                    plant
                       IF ( surf_usm_v(l)%qsws_veg(m) < 0.0_wp )  THEN
          !                 surf_usm_h%qsws_liq(m) = surf_usm_h%qsws_liq(m) + surf_usm_h%qsws_veg(m)
                          surf_usm_v(l)%qsws_veg(m) = 0.0_wp
                       ENDIF
                    ENDIF
                 
                 ENDIF
              ELSE
                 surf_usm_v(l)%r_s(m) = 1.0E10_wp
              ENDIF
!
!--           During spinup green and window fraction are set to zero. Here, the original
!--           values are restored.
              IF ( during_spinup )  THEN
                 surf_usm_v(l)%frac(m,ind_wat_win)   = frac_win
                 surf_usm_v(l)%frac(m,ind_veg_wall)  = frac_wall
                 surf_usm_v(l)%frac(m,ind_pav_green) = frac_green
              ENDIF

           ENDDO
 
       ENDDO
       !$OMP END PARALLEL

!
!--     Add-up anthropogenic heat, for now only at upward-facing surfaces
         IF ( usm_anthropogenic_heat  .AND.  .NOT. during_spinup  .AND. &
              intermediate_timestep_count == intermediate_timestep_count_max )  THEN
!
!--        application of the additional anthropogenic heat sources
!--        we considere the traffic for now so all heat is absorbed
!--        to the first layer, generalization would be worth. 
!--        calculation of actual profile coefficient
!--        ??? check time_since_reference_point ???
            CALL get_date_time( time_since_reference_point, hour=dhour, second_of_day=dtime )

!--         TO_DO: activate, if testcase is available
!--         !$OMP PARALLEL DO PRIVATE (i, j, k, acoef, rho_cp)
!--         it may also improve performance to move topo_top_ind before the k-loop
            DO i = nxl, nxr
               DO j = nys, nyn
                  DO k = nz_urban_b, min(nz_urban_t,naheatlayers)
                     IF ( k > topo_top_ind(j,i,0) ) THEN
!
!--                    increase of pt in box i,j,k in time dt_3d 
!--                    given to anthropogenic heat aheat*acoef (W*m-2)
!--                    linear interpolation of coeficient
                        acoef = (REAL(dhour+1,wp)-dtime/seconds_per_hour)*aheatprof(k,dhour) + &
                                (dtime/seconds_per_hour-REAL(dhour,wp))*aheatprof(k,dhour+1)
                        IF ( aheat(k,j,i) > 0.0_wp )  THEN
!
!--                       calculate rho * c_p coefficient at layer k
                           rho_cp  = c_p * hyp(k) / ( r_d * pt(k+1,j,i) * exner(k) )
                           pt(k,j,i) = pt(k,j,i) + aheat(k,j,i)*acoef*dt_3d/(exner(k)*rho_cp*dz(1))
                        ENDIF
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
 
         ENDIF
!
!--     pt and shf are defined on nxlg:nxrg,nysg:nyng
!--     get the borders from neighbours
         CALL exchange_horiz( pt, nbgp )
! 
!--     calculation of force_radiation_call:
!--     Make logical OR for all processes.
!--     Force radiation call if at least one processor forces it.
         IF ( intermediate_timestep_count == intermediate_timestep_count_max-1 )&
         THEN
#if defined( __parallel )
           IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
           CALL MPI_ALLREDUCE( force_radiation_call_l, force_radiation_call,    &
                               1, MPI_LOGICAL, MPI_LOR, comm2d, ierr )
#else
           force_radiation_call = force_radiation_call_l
#endif
           force_radiation_call_l = .FALSE.
         ENDIF
 
! !
! !-- Calculate surface specific humidity
!     IF ( humidity )  THEN
!        CALL calc_q_surface_usm
!     ENDIF
 
 
!     CONTAINS
! !------------------------------------------------------------------------------!
! ! Description:
! ! ------------
! !> Calculation of specific humidity of the skin layer (surface). It is assumend
! !> that the skin is always saturated.
! !------------------------------------------------------------------------------!
!        SUBROUTINE calc_q_surface_usm
! 
!           IMPLICIT NONE
! 
!           REAL(wp) :: resistance    !< aerodynamic and soil resistance term
! 
!           DO  m = 1, surf_usm_h%ns
! 
!              i   = surf_usm_h%i(m)            
!              j   = surf_usm_h%j(m)
!              k   = surf_usm_h%k(m)
! 
!!
!!--          Calculate water vapour pressure at saturation
!              e_s = 0.01_wp * 610.78_wp * EXP( 17.269_wp *                  &
!                                     ( t_surf_green_h_p(m) - 273.16_wp ) /  &
!                                     ( t_surf_green_h_p(m) - 35.86_wp  )    &
!                                          )
! 
!!
!!--          Calculate specific humidity at saturation
!              q_s = 0.622_wp * e_s / ( surface_pressure - e_s )
! 
!!              surf_usm_h%r_a_green(m) = ( surf_usm_h%pt1(m) - t_surf_green_h(m) / exner(k) ) /  &
!!                    ( surf_usm_h%ts(m) * surf_usm_h%us(m) + 1.0E-10_wp )
!!                 
!! !--          make sure that the resistance does not drop to zero
!!              IF ( ABS(surf_usm_h%r_a_green(m)) < 1.0E-10_wp )  surf_usm_h%r_a_green(m) = 1.0E-10_wp
! 
!              resistance = surf_usm_h%r_a_green(m) / ( surf_usm_h%r_a_green(m) + surf_usm_h%r_s(m) + 1E-5_wp )
! 
!!
!!--          Calculate specific humidity at surface
!              IF ( bulk_cloud_model )  THEN
!                 q(k,j,i) = resistance * q_s +                   &
!                                            ( 1.0_wp - resistance ) *              &
!                                            ( q(k,j,i) - ql(k,j,i) )
!              ELSE
!                 q(k,j,i) = resistance * q_s +                   &
!                                            ( 1.0_wp - resistance ) *              &
!                                              q(k,j,i)
!              ENDIF
! 
!!
!!--          Update virtual potential temperature
!              vpt(k,j,i) = pt(k,j,i) *         &
!                         ( 1.0_wp + 0.61_wp * q(k,j,i) )
! 
!           ENDDO
!
!!
!!--       Now, treat vertical surface elements
!           DO  l = 0, 3
!              DO  m = 1, surf_usm_v(l)%ns
!!
!!--             Get indices of respective grid point
!                 i = surf_usm_v(l)%i(m)
!                 j = surf_usm_v(l)%j(m)
!                 k = surf_usm_v(l)%k(m)
! 
!!
!!--             Calculate water vapour pressure at saturation
!                 e_s = 0.01_wp * 610.78_wp * EXP( 17.269_wp *                       &
!                                        ( t_surf_green_v_p(l)%t(m) - 273.16_wp ) /  &
!                                        ( t_surf_green_v_p(l)%t(m) - 35.86_wp  )    &
!                                             )
! 
!!
!!--             Calculate specific humidity at saturation
!                 q_s = 0.622_wp * e_s / ( surface_pressure -e_s )
! 
!!
!!--             Calculate specific humidity at surface
!                 IF ( bulk_cloud_model )  THEN
!                    q(k,j,i) = ( q(k,j,i) - ql(k,j,i) )
!                 ELSE
!                    q(k,j,i) = q(k,j,i)
!                 ENDIF
!!
!!--             Update virtual potential temperature
!                 vpt(k,j,i) = pt(k,j,i) *         &
!                            ( 1.0_wp + 0.61_wp * q(k,j,i) )
! 
!              ENDDO
! 
!           ENDDO
! 
!        END SUBROUTINE calc_q_surface_usm

        IF ( debug_output_timestep )  THEN
           WRITE( debug_string, * ) 'usm_surface_energy_balance | during_spinup: ',&
                                    during_spinup
           CALL debug_message( debug_string, 'end' )
        ENDIF

     END SUBROUTINE usm_surface_energy_balance
 
 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Swapping of timelevels for t_surf and t_wall
!> called out from subroutine swap_timelevel
!------------------------------------------------------------------------------!
     SUBROUTINE usm_swap_timelevel( mod_count )
 
        IMPLICIT NONE
 
        INTEGER(iwp), INTENT(IN) ::  mod_count
 
       
        SELECT CASE ( mod_count )
 
           CASE ( 0 )
!
!--          Horizontal surfaces
              t_surf_wall_h    => t_surf_wall_h_1;   t_surf_wall_h_p    => t_surf_wall_h_2
              t_wall_h         => t_wall_h_1;        t_wall_h_p         => t_wall_h_2
              t_surf_window_h  => t_surf_window_h_1; t_surf_window_h_p  => t_surf_window_h_2
              t_window_h       => t_window_h_1;      t_window_h_p       => t_window_h_2
              t_surf_green_h   => t_surf_green_h_1;  t_surf_green_h_p   => t_surf_green_h_2
              t_green_h        => t_green_h_1;       t_green_h_p        => t_green_h_2
!
!--          Vertical surfaces
              t_surf_wall_v    => t_surf_wall_v_1;   t_surf_wall_v_p    => t_surf_wall_v_2
              t_wall_v         => t_wall_v_1;        t_wall_v_p         => t_wall_v_2
              t_surf_window_v  => t_surf_window_v_1; t_surf_window_v_p  => t_surf_window_v_2
              t_window_v       => t_window_v_1;      t_window_v_p       => t_window_v_2
              t_surf_green_v   => t_surf_green_v_1;  t_surf_green_v_p   => t_surf_green_v_2
              t_green_v        => t_green_v_1;       t_green_v_p        => t_green_v_2
           CASE ( 1 )
!
!--          Horizontal surfaces
              t_surf_wall_h    => t_surf_wall_h_2;   t_surf_wall_h_p    => t_surf_wall_h_1
              t_wall_h         => t_wall_h_2;        t_wall_h_p         => t_wall_h_1
              t_surf_window_h  => t_surf_window_h_2; t_surf_window_h_p  => t_surf_window_h_1
              t_window_h       => t_window_h_2;      t_window_h_p       => t_window_h_1
              t_surf_green_h   => t_surf_green_h_2;  t_surf_green_h_p   => t_surf_green_h_1
              t_green_h        => t_green_h_2;       t_green_h_p        => t_green_h_1
!
!--          Vertical surfaces
              t_surf_wall_v    => t_surf_wall_v_2;   t_surf_wall_v_p    => t_surf_wall_v_1
              t_wall_v         => t_wall_v_2;        t_wall_v_p         => t_wall_v_1
              t_surf_window_v  => t_surf_window_v_2; t_surf_window_v_p  => t_surf_window_v_1
              t_window_v       => t_window_v_2;      t_window_v_p       => t_window_v_1
              t_surf_green_v   => t_surf_green_v_2;  t_surf_green_v_p   => t_surf_green_v_1
              t_green_v        => t_green_v_2;       t_green_v_p        => t_green_v_1
        END SELECT
         
     END SUBROUTINE usm_swap_timelevel
 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine writes t_surf and t_wall data into restart files
!------------------------------------------------------------------------------!
     SUBROUTINE usm_wrd_local
 
     
        IMPLICIT NONE
        
        CHARACTER(LEN=1) ::  dum     !< dummy string to create output-variable name  
        INTEGER(iwp)     ::  l       !< index surface type orientation
 
        CALL wrd_write_string( 'ns_h_on_file_usm' )
        WRITE ( 14 )  surf_usm_h%ns
 
        CALL wrd_write_string( 'ns_v_on_file_usm' )
        WRITE ( 14 )  surf_usm_v(0:3)%ns
 
        CALL wrd_write_string( 'usm_start_index_h' )
        WRITE ( 14 )  surf_usm_h%start_index
 
        CALL wrd_write_string( 'usm_end_index_h' )
        WRITE ( 14 )  surf_usm_h%end_index
 
        CALL wrd_write_string( 't_surf_wall_h' )
        WRITE ( 14 )  t_surf_wall_h
 
        CALL wrd_write_string( 't_surf_window_h' )
        WRITE ( 14 )  t_surf_window_h
 
        CALL wrd_write_string( 't_surf_green_h' )
        WRITE ( 14 )  t_surf_green_h

        CALL wrd_write_string( 'm_liq_usm_h' )
        WRITE ( 14 )  m_liq_usm_h%var_usm_1d
!
!--     Write restart data which is especially needed for the urban-surface 
!--     model. In order to do not fill up the restart routines in 
!--     surface_mod. 
!--     Output of waste heat from indoor model. Restart data is required in 
!--     this special case, because the indoor model where waste heat is 
!--     computed is call each hour (current default), so that waste heat would
!--     have zero value until next call of indoor model. 
        IF ( indoor_model )  THEN
           CALL wrd_write_string( 'waste_heat_h' )
           WRITE ( 14 )  surf_usm_h%waste_heat
        ENDIF   
           
        DO  l = 0, 3
 
           CALL wrd_write_string( 'usm_start_index_v' )
           WRITE ( 14 )  surf_usm_v(l)%start_index
 
           CALL wrd_write_string( 'usm_end_index_v' )
           WRITE ( 14 )  surf_usm_v(l)%end_index
 
           WRITE( dum, '(I1)')  l          
 
           CALL wrd_write_string( 't_surf_wall_v(' // dum // ')' )
           WRITE ( 14 )  t_surf_wall_v(l)%t
 
           CALL wrd_write_string( 't_surf_window_v(' // dum // ')' )
           WRITE ( 14 ) t_surf_window_v(l)%t     
 
           CALL wrd_write_string( 't_surf_green_v(' // dum // ')' )
           WRITE ( 14 ) t_surf_green_v(l)%t  
           
           IF ( indoor_model )  THEN
              CALL wrd_write_string( 'waste_heat_v(' // dum // ')' )
              WRITE ( 14 )  surf_usm_v(l)%waste_heat
           ENDIF
           
        ENDDO
 
        CALL wrd_write_string( 'usm_start_index_h' )
        WRITE ( 14 )  surf_usm_h%start_index
 
        CALL wrd_write_string( 'usm_end_index_h' )
        WRITE ( 14 )  surf_usm_h%end_index
 
        CALL wrd_write_string( 't_wall_h' )
        WRITE ( 14 )  t_wall_h
 
        CALL wrd_write_string( 't_window_h' )
        WRITE ( 14 )  t_window_h
 
        CALL wrd_write_string( 't_green_h' )
        WRITE ( 14 )  t_green_h
 
        DO  l = 0, 3
 
           CALL wrd_write_string( 'usm_start_index_v' )
           WRITE ( 14 )  surf_usm_v(l)%start_index
 
           CALL wrd_write_string( 'usm_end_index_v' )
           WRITE ( 14 )  surf_usm_v(l)%end_index
 
           WRITE( dum, '(I1)')  l     
 
           CALL wrd_write_string( 't_wall_v(' // dum // ')' )
           WRITE ( 14 )  t_wall_v(l)%t
 
           CALL wrd_write_string( 't_window_v(' // dum // ')' )
           WRITE ( 14 )  t_window_v(l)%t
 
           CALL wrd_write_string( 't_green_v(' // dum // ')' )
           WRITE ( 14 )  t_green_v(l)%t
        
        ENDDO
        
     END SUBROUTINE usm_wrd_local
     
     
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Define building properties
!------------------------------------------------------------------------------!
     SUBROUTINE usm_define_pars     
!
!--     Define the building_pars
        building_pars(:,1) = (/   &
           0.7_wp,         &  !< parameter 0   - wall fraction above ground floor level
           0.3_wp,         &  !< parameter 1   - window fraction above ground floor level
           0.0_wp,         &  !< parameter 2   - green fraction above ground floor level
           0.0_wp,         &  !< parameter 3   - green fraction roof above ground floor level
           1.5_wp,         &  !< parameter 4   - LAI roof
           1.5_wp,         &  !< parameter 5   - LAI on wall above ground floor level
           2200000.0_wp,   &  !< parameter 6   - heat capacity 1st/2nd wall layer above ground floor level
           1400000.0_wp,   &  !< parameter 7   - heat capacity 3rd wall layer above ground floor level
           1300000.0_wp,   &  !< parameter 8   - heat capacity 4th wall layer above ground floor level
           0.35_wp,        &  !< parameter 9   - thermal conductivity 1st/2nd wall layer above ground floor level
           0.8_wp,         &  !< parameter 10  - thermal conductivity 3rd wall layer above ground floor level
           2.1_wp,         &  !< parameter 11  - thermal conductivity 4th wall layer above ground floor level
           299.15_wp,      &  !< parameter 12  - indoor target summer temperature
           293.15_wp,      &  !< parameter 13  - indoor target winter temperature
           0.93_wp,        &  !< parameter 14  - wall emissivity above ground floor level
           0.86_wp,        &  !< parameter 15  - green emissivity above ground floor level
           0.91_wp,        &  !< parameter 16  - window emissivity above ground floor level
           0.75_wp,        &  !< parameter 17  - window transmissivity above ground floor level
           0.001_wp,       &  !< parameter 18  - z0 roughness above ground floor level
           0.0001_wp,      &  !< parameter 19  - z0h/z0g roughness heat/humidity above ground floor level
           4.0_wp,         &  !< parameter 20  - ground floor level height
           0.75_wp,        &  !< parameter 21  - wall fraction ground floor level
           0.25_wp,        &  !< parameter 22  - window fraction ground floor level
           0.0_wp,         &  !< parameter 23  - green fraction ground floor level
           0.0_wp,         &  !< parameter 24  - green fraction roof ground floor level
           1.5_wp,         &  !< parameter 25  - LAI on wall ground floor level
           2200000.0_wp,   &  !< parameter 26  - heat capacity 1st/2nd wall layer ground floor level
           1400000.0_wp,   &  !< parameter 27  - heat capacity 3rd wall layer ground floor level
           1300000.0_wp,   &  !< parameter 28  - heat capacity 4th wall layer ground floor level
           0.35_wp,        &  !< parameter 29  - thermal conductivity 1st/2nd wall layer ground floor level
           0.8_wp,         &  !< parameter 30  - thermal conductivity 3rd wall layer ground floor level
           2.1_wp,         &  !< parameter 31  - thermal conductivity 4th wall layer ground floor level
           0.93_wp,        &  !< parameter 32  - wall emissivity ground floor level
           0.91_wp,        &  !< parameter 33  - window emissivity ground floor level
           0.86_wp,        &  !< parameter 34  - green emissivity ground floor level
           0.75_wp,        &  !< parameter 35  - window transmissivity ground floor level
           0.001_wp,       &  !< parameter 36  - z0 roughness ground floor level
           0.0001_wp,      &  !< parameter 37  - z0h/z0q roughness heat/humidity
           27.0_wp,        &  !< parameter 38  - wall albedo above ground floor level
           5.0_wp,         &  !< parameter 39  - green albedo above ground floor level
           27.0_wp,        &  !< parameter 40  - window albedo above ground floor level
           0.005_wp,       &  !< parameter 41  - 1st wall layer thickness above ground floor level
           0.01_wp,        &  !< parameter 42  - 2nd wall layer thickness above ground floor level
           0.39_wp,        &  !< parameter 43  - 3rd wall layer thickness above ground floor level
           0.63_wp,        &  !< parameter 44  - 4th wall layer thickness above ground floor level
           20000.0_wp,     &  !< parameter 45  - heat capacity wall surface
           23.0_wp,        &  !< parameter 46  - thermal conductivity of wall surface
           20000.0_wp,     &  !< parameter 47  - heat capacity of window surface
           20000.0_wp,     &  !< parameter 48  - heat capacity of green surface
           23.0_wp,        &  !< parameter 49  - thermal conductivity of window surface
           10.0_wp,        &  !< parameter 50  - thermal conductivty of green surface
           1.0_wp,         &  !< parameter 51  - wall fraction ground plate
           0.005_wp,       &  !< parameter 52  - 1st wall layer thickness ground plate
           0.01_wp,        &  !< parameter 53  - 2nd wall layer thickness ground plate
           0.39_wp,        &  !< parameter 54  - 3rd wall layer thickness ground plate
           0.63_wp,        &  !< parameter 55  - 4th wall layer thickness ground plate
           2200000.0_wp,   &  !< parameter 56  - heat capacity 1st/2nd wall layer ground plate
           1400000.0_wp,   &  !< parameter 57  - heat capacity 3rd wall layer ground plate
           1300000.0_wp,   &  !< parameter 58  - heat capacity 4th wall layer ground plate
           0.35_wp,        &  !< parameter 59  - thermal conductivity 1st/2nd wall layer ground plate
           0.8_wp,         &  !< parameter 60  - thermal conductivity 3rd wall layer ground plate
           2.1_wp,         &  !< parameter 61  - thermal conductivity 4th wall layer ground plate
           0.005_wp,       &  !< parameter 62  - 1st wall layer thickness ground floor level
           0.01_wp,        &  !< parameter 63  - 2nd wall layer thickness ground floor level
           0.39_wp,        &  !< parameter 64  - 3rd wall layer thickness ground floor level
           0.63_wp,        &  !< parameter 65  - 4th wall layer thickness ground floor level
           27.0_wp,        &  !< parameter 66  - wall albedo ground floor level
           0.003_wp,       &  !< parameter 67  - 1st window layer thickness ground floor level
           0.006_wp,       &  !< parameter 68  - 2nd window layer thickness ground floor level
           0.012_wp,       &  !< parameter 69  - 3rd window layer thickness ground floor level
           0.018_wp,       &  !< parameter 70  - 4th window layer thickness ground floor level
           1736000.0_wp,   &  !< parameter 71  - heat capacity 1st/2nd window layer ground floor level
           1736000.0_wp,   &  !< parameter 72  - heat capacity 3rd window layer ground floor level
           1736000.0_wp,   &  !< parameter 73  - heat capacity 4th window layer ground floor level
           0.57_wp,        &  !< parameter 74  - thermal conductivity 1st/2nd window layer ground floor level
           0.57_wp,        &  !< parameter 75  - thermal conductivity 3rd window layer ground floor level
           0.57_wp,        &  !< parameter 76  - thermal conductivity 4th window layer ground floor level
           27.0_wp,        &  !< parameter 77  - window albedo ground floor level
           5.0_wp,         &  !< parameter 78  - green albedo ground floor level
           0.003_wp,       &  !< parameter 79  - 1st window layer thickness above ground floor level
           0.006_wp,       &  !< parameter 80  - 2nd thickness window layer above ground floor level
           0.012_wp,       &  !< parameter 81  - 3rd window layer thickness above ground floor level
           0.018_wp,       &  !< parameter 82  - 4th window layer thickness above ground floor level
           1736000.0_wp,   &  !< parameter 83  - heat capacity 1st/2nd window layer above ground floor level
           1736000.0_wp,   &  !< parameter 84  - heat capacity 3rd window layer above ground floor level
           1736000.0_wp,   &  !< parameter 85  - heat capacity 4th window layer above ground floor level
           0.57_wp,        &  !< parameter 86  - thermal conductivity 1st/2nd window layer above ground floor level
           0.57_wp,        &  !< parameter 87  - thermal conductivity 3rd window layer above ground floor level
           0.57_wp,        &  !< parameter 88  - thermal conductivity 4th window layer above ground floor level
           1.0_wp,         &  !< parameter 89  - wall fraction roof
           0.005_wp,       &  !< parameter 90  - 1st wall layer thickness roof
           0.01_wp,        &  !< parameter 91  - 2nd wall layer thickness roof
           0.31_wp,        &  !< parameter 92  - 3rd wall layer thickness roof
           0.63_wp,        &  !< parameter 93  - 4th wall layer thickness roof
           2200000.0_wp,   &  !< parameter 94  - heat capacity 1st/2nd wall layer roof
           1400000.0_wp,   &  !< parameter 95  - heat capacity 3rd wall layer roof
           1300000.0_wp,   &  !< parameter 96  - heat capacity 4th wall layer roof
           0.35_wp,        &  !< parameter 97  - thermal conductivity 1st/2nd wall layer roof
           0.8_wp,         &  !< parameter 98  - thermal conductivity 3rd wall layer roof
           2.1_wp,         &  !< parameter 99  - thermal conductivity 4th wall layer roof
           0.93_wp,        &  !< parameter 100 - wall emissivity roof
           27.0_wp,        &  !< parameter 101 - wall albedo roof
           0.0_wp,         &  !< parameter 102 - window fraction roof
           0.003_wp,       &  !< parameter 103 - window 1st layer thickness roof
           0.006_wp,       &  !< parameter 104 - window 2nd layer thickness roof
           0.012_wp,       &  !< parameter 105 - window 3rd layer thickness roof
           0.018_wp,       &  !< parameter 106 - window 4th layer thickness roof
           1736000.0_wp,   &  !< parameter 107 - heat capacity 1st/2nd window layer roof
           1736000.0_wp,   &  !< parameter 108 - heat capacity 3rd window layer roof
           1736000.0_wp,   &  !< parameter 109 - heat capacity 4th window layer roof
           0.57_wp,        &  !< parameter 110 - thermal conductivity 1st/2nd window layer roof
           0.57_wp,        &  !< parameter 111 - thermal conductivity 3rd window layer roof
           0.57_wp,        &  !< parameter 112 - thermal conductivity 4th window layer roof
           0.91_wp,        &  !< parameter 113 - window emissivity roof
           0.75_wp,        &  !< parameter 114 - window transmissivity roof
           27.0_wp,        &  !< parameter 115 - window albedo roof
           0.86_wp,        &  !< parameter 116 - green emissivity roof
           5.0_wp,         &  !< parameter 117 - green albedo roof
           0.0_wp,         &  !< parameter 118 - green type roof
           0.8_wp,         &  !< parameter 119 - shading factor
           0.76_wp,        &  !< parameter 120 - g-value windows
           5.0_wp,         &  !< parameter 121 - u-value windows
           0.5_wp,         &  !< parameter 122 - basical airflow without occupancy of the room for - summer 0.5_wp, winter 0.1
           2.0_wp,         &  !< parameter 123 - additional airflow depend of occupancy of the room for - summer 2.0_wp, winter 0.5
           0.0_wp,         &  !< parameter 124 - heat recovery efficiency
           3.5_wp,         &  !< parameter 125 - dynamic parameter specific effective surface
           370000.0_wp,    &  !< parameter 126 - dynamic parameter innner heatstorage
           4.5_wp,         &  !< parameter 127 - ratio internal surface/floor area
           100.0_wp,       &  !< parameter 128 - maximal heating capacity
           0.0_wp,         &  !< parameter 129 - maximal cooling capacity
           2.0_wp,         &  !< parameter 130 - additional internal heat gains dependent on occupancy of the room
           6.0_wp,         &  !< parameter 131 - basic internal heat gains without occupancy of the room
           3.0_wp,         &  !< parameter 132 - storey height
           0.2_wp,         &  !< parameter 133 - ceiling construction height
           0.1_wp,         &  !< parameter 134 - anthropogenic heat output for heating
           1.333_wp        &  !< parameter 135 - anthropogenic heat output for cooling
                            /)
                            
     building_pars(:,2) = (/   &
           0.73_wp,        &  !< parameter 0   - wall fraction above ground floor level
           0.27_wp,        &  !< parameter 1   - window fraction above ground floor level
           0.0_wp,         &  !< parameter 2   - green fraction above ground floor level
           0.0_wp,         &  !< parameter 3   - green fraction roof above ground floor level
           1.5_wp,         &  !< parameter 4   - LAI roof
           1.5_wp,         &  !< parameter 5   - LAI on wall above ground floor level
           2000000.0_wp,   &  !< parameter 6   - heat capacity 1st/2nd wall layer above ground floor level
           103000.0_wp,    &  !< parameter 7   - heat capacity 3rd wall layer above ground floor level
           900000.0_wp,    &  !< parameter 8   - heat capacity 4th wall layer above ground floor level
           0.35_wp,        &  !< parameter 9   - thermal conductivity 1st/2nd wall layer above ground floor level
           0.38_wp,        &  !< parameter 10  - thermal conductivity 3rd wall layer above ground floor level
           0.04_wp,        &  !< parameter 11  - thermal conductivity 4th wall layer above ground floor level
           299.15_wp,      &  !< parameter 12  - indoor target summer temperature
           293.15_wp,      &  !< parameter 13  - indoor target winter temperature
           0.92_wp,        &  !< parameter 14  - wall emissivity above ground floor level
           0.86_wp,        &  !< parameter 15  - green emissivity above ground floor level
           0.87_wp,        &  !< parameter 16  - window emissivity above ground floor level
           0.7_wp,         &  !< parameter 17  - window transmissivity above ground floor level
           0.001_wp,       &  !< parameter 18  - z0 roughness above ground floor level
           0.0001_wp,      &  !< parameter 19  - z0h/z0g roughness heat/humidity above ground floor level
           4.0_wp,         &  !< parameter 20  - ground floor level height
           0.78_wp,        &  !< parameter 21  - wall fraction ground floor level
           0.22_wp,        &  !< parameter 22  - window fraction ground floor level
           0.0_wp,         &  !< parameter 23  - green fraction ground floor level
           0.0_wp,         &  !< parameter 24  - green fraction roof ground floor level
           1.5_wp,         &  !< parameter 25  - LAI on wall ground floor level
           2000000.0_wp,   &  !< parameter 26  - heat capacity 1st/2nd wall layer ground floor level
           103000.0_wp,    &  !< parameter 27  - heat capacity 3rd wall layer ground floor level
           900000.0_wp,    &  !< parameter 28  - heat capacity 4th wall layer ground floor level
           0.35_wp,        &  !< parameter 29  - thermal conductivity 1st/2nd wall layer ground floor level
           0.38_wp,        &  !< parameter 30  - thermal conductivity 3rd wall layer ground floor level
           0.04_wp,        &  !< parameter 31  - thermal conductivity 4th wall layer ground floor level
           0.92_wp,        &  !< parameter 32  - wall emissivity ground floor level
           0.11_wp,        &  !< parameter 33  - window emissivity ground floor level
           0.86_wp,        &  !< parameter 34  - green emissivity ground floor level
           0.7_wp,         &  !< parameter 35  - window transmissivity ground floor level
           0.001_wp,       &  !< parameter 36  - z0 roughness ground floor level
           0.0001_wp,      &  !< parameter 37  - z0h/z0q roughness heat/humidity
           27.0_wp,        &  !< parameter 38  - wall albedo above ground floor level
           5.0_wp,         &  !< parameter 39  - green albedo above ground floor level
           27.0_wp,        &  !< parameter 40  - window albedo above ground floor level
           0.005_wp,       &  !< parameter 41  - 1st wall layer thickness above ground floor level
           0.01_wp,        &  !< parameter 42  - 2nd wall layer thickness above ground floor level
           0.31_wp,        &  !< parameter 43  - 3rd wall layer thickness above ground floor level
           0.43_wp,        &  !< parameter 44  - 4th wall layer thickness above ground floor level
           20000.0_wp,     &  !< parameter 45  - heat capacity wall surface
           23.0_wp,        &  !< parameter 46  - thermal conductivity of wall surface
           20000.0_wp,     &  !< parameter 47  - heat capacity of window surface
           20000.0_wp,     &  !< parameter 48  - heat capacity of green surface
           23.0_wp,        &  !< parameter 49  - thermal conductivity of window surface
           10.0_wp,        &  !< parameter 50  - thermal conductivty of green surface
           1.0_wp,         &  !< parameter 51  - wall fraction ground plate
           0.005_wp,       &  !< parameter 52  - 1st wall layer thickness ground plate
           0.01_wp,        &  !< parameter 53  - 2nd wall layer thickness ground plate
           0.31_wp,        &  !< parameter 54  - 3rd wall layer thickness ground plate
           0.42_wp,        &  !< parameter 55  - 4th wall layer thickness ground plate
           2000000.0_wp,   &  !< parameter 56  - heat capacity 1st/2nd wall layer ground plate
           103000.0_wp,    &  !< parameter 57  - heat capacity 3rd wall layer ground plate
           900000.0_wp,    &  !< parameter 58  - heat capacity 4th wall layer ground plate
           0.35_wp,        &  !< parameter 59  - thermal conductivity 1st/2nd wall layer ground plate
           0.38_wp,        &  !< parameter 60  - thermal conductivity 3rd wall layer ground plate
           0.04_wp,        &  !< parameter 61  - thermal conductivity 4th wall layer ground plate
           0.005_wp,       &  !< parameter 62  - 1st wall layer thickness ground floor level
           0.01_wp,        &  !< parameter 63  - 2nd wall layer thickness ground floor level
           0.31_wp,        &  !< parameter 64  - 3rd wall layer thickness ground floor level
           0.43_wp,        &  !< parameter 65  - 4th wall layer thickness ground floor level
           27.0_wp,        &  !< parameter 66  - wall albedo ground floor level
           0.003_wp,       &  !< parameter 67  - 1st window layer thickness ground floor level
           0.006_wp,       &  !< parameter 68  - 2nd window layer thickness ground floor level
           0.012_wp,       &  !< parameter 69  - 3rd window layer thickness ground floor level
           0.018_wp,       &  !< parameter 70  - 4th window layer thickness ground floor level
           1736000.0_wp,   &  !< parameter 71  - heat capacity 1st/2nd window layer ground floor level
           1736000.0_wp,   &  !< parameter 72  - heat capacity 3rd window layer ground floor level
           1736000.0_wp,   &  !< parameter 73  - heat capacity 4th window layer ground floor level
           0.11_wp,        &  !< parameter 74  - thermal conductivity 1st/2nd window layer ground floor level
           0.11_wp,        &  !< parameter 75  - thermal conductivity 3rd window layer ground floor level
           0.11_wp,        &  !< parameter 76  - thermal conductivity 4th window layer ground floor level
           27.0_wp,        &  !< parameter 77  - window albedo ground floor level
           5.0_wp,         &  !< parameter 78  - green albedo ground floor level
           0.003_wp,       &  !< parameter 79  - 1st window layer thickness above ground floor level
           0.006_wp,       &  !< parameter 80  - 2nd thickness window layer above ground floor level
           0.012_wp,       &  !< parameter 81  - 3rd window layer thickness above ground floor level
           0.018_wp,       &  !< parameter 82  - 4th window layer thickness above ground floor level
           1736000.0_wp,   &  !< parameter 83  - heat capacity 1st/2nd window layer above ground floor level
           1736000.0_wp,   &  !< parameter 84  - heat capacity 3rd window layer above ground floor level
           1736000.0_wp,   &  !< parameter 85  - heat capacity 4th window layer above ground floor level
           0.11_wp,        &  !< parameter 86  - thermal conductivity 1st/2nd window layer above ground floor level
           0.11_wp,        &  !< parameter 87  - thermal conductivity 3rd window layer above ground floor level
           0.11_wp,        &  !< parameter 88  - thermal conductivity 4th window layer above ground floor level
           1.0_wp,         &  !< parameter 89  - wall fraction roof
           0.005_wp,       &  !< parameter 90  - 1st wall layer thickness roof
           0.01_wp,        &  !< parameter 91  - 2nd wall layer thickness roof
           0.5_wp,         &  !< parameter 92  - 3rd wall layer thickness roof
           0.79_wp,        &  !< parameter 93  - 4th wall layer thickness roof
           2000000.0_wp,   &  !< parameter 94  - heat capacity 1st/2nd wall layer roof
           103000.0_wp,    &  !< parameter 95  - heat capacity 3rd wall layer roof
           900000.0_wp,    &  !< parameter 96  - heat capacity 4th wall layer roof
           0.35_wp,        &  !< parameter 97  - thermal conductivity 1st/2nd wall layer roof
           0.38_wp,        &  !< parameter 98  - thermal conductivity 3rd wall layer roof
           0.04_wp,        &  !< parameter 99  - thermal conductivity 4th wall layer roof
           0.93_wp,        &  !< parameter 100 - wall emissivity roof
           27.0_wp,        &  !< parameter 101 - wall albedo roof
           0.0_wp,         &  !< parameter 102 - window fraction roof
           0.003_wp,       &  !< parameter 103 - window 1st layer thickness roof
           0.006_wp,       &  !< parameter 104 - window 2nd layer thickness roof
           0.012_wp,       &  !< parameter 105 - window 3rd layer thickness roof
           0.018_wp,       &  !< parameter 106 - window 4th layer thickness roof
           1736000.0_wp,   &  !< parameter 107 - heat capacity 1st/2nd window layer roof
           1736000.0_wp,   &  !< parameter 108 - heat capacity 3rd window layer roof
           1736000.0_wp,   &  !< parameter 109 - heat capacity 4th window layer roof
           0.11_wp,        &  !< parameter 110 - thermal conductivity 1st/2nd window layer roof
           0.11_wp,        &  !< parameter 111 - thermal conductivity 3rd window layer roof
           0.11_wp,        &  !< parameter 112 - thermal conductivity 4th window layer roof
           0.87_wp,        &  !< parameter 113 - window emissivity roof
           0.7_wp,         &  !< parameter 114 - window transmissivity roof
           27.0_wp,        &  !< parameter 115 - window albedo roof
           0.86_wp,        &  !< parameter 116 - green emissivity roof
           5.0_wp,         &  !< parameter 117 - green albedo roof
           0.0_wp,         &  !< parameter 118 - green type roof
           0.8_wp,         &  !< parameter 119 - shading factor
           0.6_wp,         &  !< parameter 120 - g-value windows
           3.0_wp,         &  !< parameter 121 - u-value windows
           0.5_wp,         &  !< parameter 122 - basical airflow without occupancy of the room for - summer 0.5_wp for winter 0.1
           2.0_wp,         &  !< parameter 123 - additional airflow depend of occupancy of the room for - summer 2.0_wp for winter 0.5
           0.0_wp,         &  !< parameter 124 - heat recovery efficiency
           2.5_wp,         &  !< parameter 125 - dynamic parameter specific effective surface
           165000.0_wp,    &  !< parameter 126 - dynamic parameter innner heatstorage
           4.5_wp,         &  !< parameter 127 - ratio internal surface/floor area
           100.0_wp,       &  !< parameter 128 - maximal heating capacity
           0.0_wp,         &  !< parameter 129 - maximal cooling capacity
           2.0_wp,         &  !< parameter 130 - additional internal heat gains dependent on occupancy of the room
           6.0_wp,         &  !< parameter 131 - basic internal heat gains without occupancy of the room
           3.0_wp,         &  !< parameter 132 - storey height
           0.2_wp,         &  !< parameter 133 - ceiling construction height
           0.1_wp,         &  !< parameter 134 - anthropogenic heat output for heating
           1.333_wp        &  !< parameter 135 - anthropogenic heat output for cooling
                            /)

  building_pars(:,3) = (/   &
           0.7_wp,         &  !< parameter 0   - wall fraction above ground floor level
           0.3_wp,         &  !< parameter 1   - window fraction above ground floor level
           0.0_wp,         &  !< parameter 2   - green fraction above ground floor level
           0.0_wp,         &  !< parameter 3   - green fraction roof above ground floor level
           1.5_wp,         &  !< parameter 4   - LAI roof
           1.5_wp,         &  !< parameter 5   - LAI on wall above ground floor level
           2000000.0_wp,   &  !< parameter 6   - heat capacity 1st/2nd wall layer above ground floor level
           103000.0_wp,    &  !< parameter 7   - heat capacity 3rd wall layer above ground floor level
           900000.0_wp,    &  !< parameter 8   - heat capacity 4th wall layer above ground floor level
           0.35_wp,        &  !< parameter 9   - thermal conductivity 1st/2nd wall layer above ground floor level
           0.14_wp,        &  !< parameter 10  - thermal conductivity 3rd wall layer above ground floor level
           0.035_wp,       &  !< parameter 11  - thermal conductivity 4th wall layer above ground floor level
           299.15_wp,      &  !< parameter 12  - indoor target summer temperature
           293.15_wp,      &  !< parameter 13  - indoor target winter temperature
           0.92_wp,        &  !< parameter 14  - wall emissivity above ground floor level
           0.86_wp,        &  !< parameter 15  - green emissivity above ground floor level
           0.8_wp,         &  !< parameter 16  - window emissivity above ground floor level
           0.6_wp,         &  !< parameter 17  - window transmissivity above ground floor level
           0.001_wp,       &  !< parameter 18  - z0 roughness above ground floor level
           0.0001_wp,      &  !< parameter 19  - z0h/z0g roughness heat/humidity above ground floor level
           3.0_wp,         &  !< parameter 20  - ground floor level height
           0.75_wp,        &  !< parameter 21  - wall fraction ground floor level
           0.25_wp,        &  !< parameter 22  - window fraction ground floor level
           0.0_wp,         &  !< parameter 23  - green fraction ground floor level
           0.0_wp,         &  !< parameter 24  - green fraction roof ground floor level
           1.5_wp,         &  !< parameter 25  - LAI on wall ground floor level
           2000000.0_wp,   &  !< parameter 26  - heat capacity 1st/2nd wall layer ground floor level
           103000.0_wp,    &  !< parameter 27  - heat capacity 3rd wall layer ground floor level
           900000.0_wp,    &  !< parameter 28  - heat capacity 4th wall layer ground floor level
           0.35_wp,        &  !< parameter 29  - thermal conductivity 1st/2nd wall layer ground floor level
           0.14_wp,        &  !< parameter 30  - thermal conductivity 3rd wall layer ground floor level
           0.035_wp,       &  !< parameter 31  - thermal conductivity 4th wall layer ground floor level
           0.92_wp,        &  !< parameter 32  - wall emissivity ground floor level
           0.8_wp,         &  !< parameter 33  - window emissivity ground floor level
           0.86_wp,        &  !< parameter 34  - green emissivity ground floor level
           0.6_wp,         &  !< parameter 35  - window transmissivity ground floor level
           0.001_wp,       &  !< parameter 36  - z0 roughness ground floor level
           0.0001_wp,      &  !< parameter 37  - z0h/z0q roughness heat/humidity
           27.0_wp,        &  !< parameter 38  - wall albedo above ground floor level
           5.0_wp,         &  !< parameter 39  - green albedo above ground floor level
           27.0_wp,        &  !< parameter 40  - window albedo above ground floor level
           0.005_wp,       &  !< parameter 41  - 1st wall layer thickness above ground floor level
           0.01_wp,        &  !< parameter 42  - 2nd wall layer thickness above ground floor level
           0.41_wp,        &  !< parameter 43  - 3rd wall layer thickness above ground floor level
           0.7_wp,         &  !< parameter 44  - 4th wall layer thickness above ground floor level
           20000.0_wp,     &  !< parameter 45  - heat capacity wall surface
           23.0_wp,        &  !< parameter 46  - thermal conductivity of wall surface
           20000.0_wp,     &  !< parameter 47  - heat capacity of window surface
           20000.0_wp,     &  !< parameter 48  - heat capacity of green surface
           23.0_wp,        &  !< parameter 49  - thermal conductivity of window surface
           10.0_wp,        &  !< parameter 50  - thermal conductivty of green surface
           1.0_wp,         &  !< parameter 51  - wall fraction ground plate
           0.005_wp,       &  !< parameter 52  - 1st wall layer thickness ground plate
           0.01_wp,        &  !< parameter 53  - 2nd wall layer thickness ground plate
           0.41_wp,        &  !< parameter 54  - 3rd wall layer thickness ground plate
           0.7_wp,         &  !< parameter 55  - 4th wall layer thickness ground plate
           2000000.0_wp,   &  !< parameter 56  - heat capacity 1st/2nd wall layer ground plate
           103000.0_wp,    &  !< parameter 57  - heat capacity 3rd wall layer ground plate
           900000.0_wp,    &  !< parameter 58  - heat capacity 4th wall layer ground plate
           0.35_wp,        &  !< parameter 59  - thermal conductivity 1st/2nd wall layer ground plate
           0.14_wp,        &  !< parameter 60  - thermal conductivity 3rd wall layer ground plate
           0.035_wp,       &  !< parameter 61  - thermal conductivity 4th wall layer ground plate
           0.005_wp,       &  !< parameter 62  - 1st wall layer thickness ground floor level
           0.01_wp,        &  !< parameter 63  - 2nd wall layer thickness ground floor level
           0.41_wp,        &  !< parameter 64  - 3rd wall layer thickness ground floor level
           0.7_wp,         &  !< parameter 65  - 4th wall layer thickness ground floor level
           27.0_wp,        &  !< parameter 66  - wall albedo ground floor level
           0.003_wp,       &  !< parameter 67  - 1st window layer thickness ground floor level
           0.006_wp,       &  !< parameter 68  - 2nd window layer thickness ground floor level
           0.012_wp,       &  !< parameter 69  - 3rd window layer thickness ground floor level
           0.018_wp,       &  !< parameter 70  - 4th window layer thickness ground floor level
           1736000.0_wp,   &  !< parameter 71  - heat capacity 1st/2nd window layer ground floor level
           1736000.0_wp,   &  !< parameter 72  - heat capacity 3rd window layer ground floor level
           1736000.0_wp,   &  !< parameter 73  - heat capacity 4th window layer ground floor level
           0.037_wp,       &  !< parameter 74  - thermal conductivity 1st/2nd window layer ground floor level
           0.037_wp,       &  !< parameter 75  - thermal conductivity 3rd window layer ground floor level
           0.037_wp,       &  !< parameter 76  - thermal conductivity 4th window layer ground floor level
           27.0_wp,        &  !< parameter 77  - window albedo ground floor level
           5.0_wp,         &  !< parameter 78  - green albedo ground floor level
           0.003_wp,       &  !< parameter 79  - 1st window layer thickness above ground floor level
           0.006_wp,       &  !< parameter 80  - 2nd thickness window layer above ground floor level
           0.012_wp,       &  !< parameter 81  - 3rd window layer thickness above ground floor level
           0.018_wp,       &  !< parameter 82  - 4th window layer thickness above ground floor level
           1736000.0_wp,   &  !< parameter 83  - heat capacity 1st/2nd window layer above ground floor level
           1736000.0_wp,   &  !< parameter 84  - heat capacity 3rd window layer above ground floor level
           1736000.0_wp,   &  !< parameter 85  - heat capacity 4th window layer above ground floor level
           0.037_wp,       &  !< parameter 86  - thermal conductivity 1st/2nd window layer above ground floor level
           0.037_wp,       &  !< parameter 87  - thermal conductivity 3rd window layer above ground floor level
           0.037_wp,       &  !< parameter 88  - thermal conductivity 4th window layer above ground floor level
           1.0_wp,         &  !< parameter 89  - wall fraction roof
           0.005_wp,       &  !< parameter 90  - 1st wall layer thickness roof
           0.01_wp,        &  !< parameter 91  - 2nd wall layer thickness roof
           0.41_wp,        &  !< parameter 92  - 3rd wall layer thickness roof
           0.7_wp,         &  !< parameter 93  - 4th wall layer thickness roof
           2000000.0_wp,   &  !< parameter 94  - heat capacity 1st/2nd wall layer roof
           103000.0_wp,    &  !< parameter 95  - heat capacity 3rd wall layer roof
           900000.0_wp,    &  !< parameter 96  - heat capacity 4th wall layer roof
           0.35_wp,        &  !< parameter 97  - thermal conductivity 1st/2nd wall layer roof
           0.14_wp,        &  !< parameter 98  - thermal conductivity 3rd wall layer roof
           0.035_wp,       &  !< parameter 99  - thermal conductivity 4th wall layer roof
           0.93_wp,        &  !< parameter 100 - wall emissivity roof
           27.0_wp,        &  !< parameter 101 - wall albedo roof
           0.0_wp,         &  !< parameter 102 - window fraction roof
           0.003_wp,       &  !< parameter 103 - window 1st layer thickness roof
           0.006_wp,       &  !< parameter 104 - window 2nd layer thickness roof
           0.012_wp,       &  !< parameter 105 - window 3rd layer thickness roof
           0.018_wp,       &  !< parameter 106 - window 4th layer thickness roof
           1736000.0_wp,   &  !< parameter 107 - heat capacity 1st/2nd window layer roof
           1736000.0_wp,   &  !< parameter 108 - heat capacity 3rd window layer roof
           1736000.0_wp,   &  !< parameter 109 - heat capacity 4th window layer roof
           0.037_wp,       &  !< parameter 110 - thermal conductivity 1st/2nd window layer roof
           0.037_wp,       &  !< parameter 111 - thermal conductivity 3rd window layer roof
           0.037_wp,       &  !< parameter 112 - thermal conductivity 4th window layer roof
           0.8_wp,         &  !< parameter 113 - window emissivity roof
           0.6_wp,         &  !< parameter 114 - window transmissivity roof
           27.0_wp,        &  !< parameter 115 - window albedo roof
           0.86_wp,        &  !< parameter 116 - green emissivity roof
           5.0_wp,         &  !< parameter 117 - green albedo roof
           0.0_wp,         &  !< parameter 118 - green type roof
           0.3_wp,         &  !< parameter 119 - shading factor
           0.5_wp,         &  !< parameter 120 - g-value windows
           1.0_wp,         &  !< parameter 121 - u-value windows
           0.8_wp,         &  !< parameter 122 - basical airflow without occupancy of the room for - summer 0.8_wp, winter 0.1
           2.0_wp,         &  !< parameter 123 - additional airflow depend of occupancy of the room for - summer 2.0_wp, winter 0.5
           0.8_wp,         &  !< parameter 124 - heat recovery efficiency
           2.5_wp,         &  !< parameter 125 - dynamic parameter specific effective surface
           80000.0_wp,     &  !< parameter 126 - dynamic parameter innner heatstorage
           4.5_wp,         &  !< parameter 127 - ratio internal surface/floor area
           100.0_wp,       &  !< parameter 128 - maximal heating capacity
           0.0_wp,         &  !< parameter 129 - maximal cooling capacity
           2.0_wp,         &  !< parameter 130 - additional internal heat gains dependent on occupancy of the room
           6.0_wp,         &  !< parameter 131 - basic internal heat gains without occupancy of the room
           3.0_wp,         &  !< parameter 132 - storey height
           0.2_wp,         &  !< parameter 133 - ceiling construction height
           -2.0_wp,        &  !< parameter 134 - anthropogenic heat output for heating
           1.25_wp         &  !< parameter 135 - anthropogenic heat output for cooling
                            /)    
                            
        building_pars(:,4) = (/   &
           0.5_wp,         &  !< parameter 0   - wall fraction above ground floor level
           0.5_wp,         &  !< parameter 1   - window fraction above ground floor level
           0.0_wp,         &  !< parameter 2   - green fraction above ground floor level
           0.0_wp,         &  !< parameter 3   - green fraction roof above ground floor level
           1.5_wp,         &  !< parameter 4   - LAI roof
           1.5_wp,         &  !< parameter 5   - LAI on wall above ground floor level
           2200000.0_wp,   &  !< parameter 6   - heat capacity 1st/2nd wall layer above ground floor level
           1400000.0_wp,   &  !< parameter 7   - heat capacity 3rd wall layer above ground floor level
           1300000.0_wp,   &  !< parameter 8   - heat capacity 4th wall layer above ground floor level
           0.35_wp,        &  !< parameter 9   - thermal conductivity 1st/2nd wall layer above ground floor level
           0.8_wp,         &  !< parameter 10  - thermal conductivity 3rd wall layer above ground floor level
           2.1_wp,         &  !< parameter 11  - thermal conductivity 4th wall layer above ground floor level
           299.15_wp,      &  !< parameter 12  - indoor target summer temperature
           293.15_wp,      &  !< parameter 13  - indoor target winter temperature
           0.93_wp,        &  !< parameter 14  - wall emissivity above ground floor level
           0.86_wp,        &  !< parameter 15  - green emissivity above ground floor level
           0.91_wp,        &  !< parameter 16  - window emissivity above ground floor level
           0.75_wp,        &  !< parameter 17  - window transmissivity above ground floor level
           0.001_wp,       &  !< parameter 18  - z0 roughness above ground floor level
           0.0001_wp,      &  !< parameter 19  - z0h/z0g roughness heat/humidity above ground floor level
           4.0_wp,         &  !< parameter 20  - ground floor level height
           0.55_wp,        &  !< parameter 21  - wall fraction ground floor level
           0.45_wp,        &  !< parameter 22  - window fraction ground floor level
           0.0_wp,         &  !< parameter 23  - green fraction ground floor level
           0.0_wp,         &  !< parameter 24  - green fraction roof ground floor level
           1.5_wp,         &  !< parameter 25  - LAI on wall ground floor level
           2200000.0_wp,   &  !< parameter 26  - heat capacity 1st/2nd wall layer ground floor level
           1400000.0_wp,   &  !< parameter 27  - heat capacity 3rd wall layer ground floor level
           1300000.0_wp,   &  !< parameter 28  - heat capacity 4th wall layer ground floor level
           0.35_wp,        &  !< parameter 29  - thermal conductivity 1st/2nd wall layer ground floor level
           0.8_wp,         &  !< parameter 30  - thermal conductivity 3rd wall layer ground floor level
           2.1_wp,         &  !< parameter 31  - thermal conductivity 4th wall layer ground floor level
           0.93_wp,        &  !< parameter 32  - wall emissivity ground floor level
           0.91_wp,        &  !< parameter 33  - window emissivity ground floor level
           0.86_wp,        &  !< parameter 34  - green emissivity ground floor level
           0.75_wp,        &  !< parameter 35  - window transmissivity ground floor level
           0.001_wp,       &  !< parameter 36  - z0 roughness ground floor level
           0.0001_wp,      &  !< parameter 37  - z0h/z0q roughness heat/humidity
           27.0_wp,        &  !< parameter 38  - wall albedo above ground floor level
           5.0_wp,         &  !< parameter 39  - green albedo above ground floor level
           27.0_wp,        &  !< parameter 40  - window albedo above ground floor level
           0.005_wp,       &  !< parameter 41  - 1st wall layer thickness above ground floor level
           0.01_wp,        &  !< parameter 42  - 2nd wall layer thickness above ground floor level
           0.39_wp,        &  !< parameter 43  - 3rd wall layer thickness above ground floor level
           0.63_wp,        &  !< parameter 44  - 4th wall layer thickness above ground floor level
           20000.0_wp,     &  !< parameter 45  - heat capacity wall surface
           23.0_wp,        &  !< parameter 46  - thermal conductivity of wall surface
           20000.0_wp,     &  !< parameter 47  - heat capacity of window surface
           20000.0_wp,     &  !< parameter 48  - heat capacity of green surface
           23.0_wp,        &  !< parameter 49  - thermal conductivity of window surface
           10.0_wp,        &  !< parameter 50  - thermal conductivty of green surface
           1.0_wp,         &  !< parameter 51  - wall fraction ground plate
           0.005_wp,       &  !< parameter 52  - 1st wall layer thickness ground plate
           0.01_wp,        &  !< parameter 53  - 2nd wall layer thickness ground plate
           0.39_wp,        &  !< parameter 54  - 3rd wall layer thickness ground plate
           0.63_wp,        &  !< parameter 55  - 4th wall layer thickness ground plate
           2200000.0_wp,   &  !< parameter 56  - heat capacity 1st/2nd wall layer ground plate
           1400000.0_wp,   &  !< parameter 57  - heat capacity 3rd wall layer ground plate
           1300000.0_wp,   &  !< parameter 58  - heat capacity 4th wall layer ground plate
           0.35_wp,        &  !< parameter 59  - thermal conductivity 1st/2nd wall layer ground plate
           0.8_wp,         &  !< parameter 60  - thermal conductivity 3rd wall layer ground plate
           2.1_wp,         &  !< parameter 61  - thermal conductivity 4th wall layer ground plate
           0.005_wp,       &  !< parameter 62  - 1st wall layer thickness ground floor level
           0.01_wp,        &  !< parameter 63  - 2nd wall layer thickness ground floor level
           0.39_wp,        &  !< parameter 64  - 3rd wall layer thickness ground floor level
           0.63_wp,        &  !< parameter 65  - 4th wall layer thickness ground floor level
           27.0_wp,        &  !< parameter 66  - wall albedo ground floor level
           0.003_wp,       &  !< parameter 67  - 1st window layer thickness ground floor level
           0.006_wp,       &  !< parameter 68  - 2nd window layer thickness ground floor level
           0.012_wp,       &  !< parameter 69  - 3rd window layer thickness ground floor level
           0.018_wp,       &  !< parameter 70  - 4th window layer thickness ground floor level
           1736000.0_wp,   &  !< parameter 71  - heat capacity 1st/2nd window layer ground floor level
           1736000.0_wp,   &  !< parameter 72  - heat capacity 3rd window layer ground floor level
           1736000.0_wp,   &  !< parameter 73  - heat capacity 4th window layer ground floor level
           0.57_wp,        &  !< parameter 74  - thermal conductivity 1st/2nd window layer ground floor level
           0.57_wp,        &  !< parameter 75  - thermal conductivity 3rd window layer ground floor level
           0.57_wp,        &  !< parameter 76  - thermal conductivity 4th window layer ground floor level
           27.0_wp,        &  !< parameter 77  - window albedo ground floor level
           5.0_wp,         &  !< parameter 78  - green albedo ground floor level
           0.003_wp,       &  !< parameter 79  - 1st window layer thickness above ground floor level
           0.006_wp,       &  !< parameter 80  - 2nd thickness window layer above ground floor level
           0.012_wp,       &  !< parameter 81  - 3rd window layer thickness above ground floor level
           0.018_wp,       &  !< parameter 82  - 4th window layer thickness above ground floor level
           1736000.0_wp,   &  !< parameter 83  - heat capacity 1st/2nd window layer above ground floor level
           1736000.0_wp,   &  !< parameter 84  - heat capacity 3rd window layer above ground floor level
           1736000.0_wp,   &  !< parameter 85  - heat capacity 4th window layer above ground floor level
           0.57_wp,        &  !< parameter 86  - thermal conductivity 1st/2nd window layer above ground floor level
           0.57_wp,        &  !< parameter 87  - thermal conductivity 3rd window layer above ground floor level
           0.57_wp,        &  !< parameter 88  - thermal conductivity 4th window layer above ground floor level
           1.0_wp,         &  !< parameter 89  - wall fraction roof
           0.005_wp,       &  !< parameter 90  - 1st wall layer thickness roof
           0.01_wp,        &  !< parameter 91  - 2nd wall layer thickness roof
           0.39_wp,        &  !< parameter 92  - 3rd wall layer thickness roof
           0.63_wp,        &  !< parameter 93  - 4th wall layer thickness roof
           2200000.0_wp,   &  !< parameter 94  - heat capacity 1st/2nd wall layer roof
           1400000.0_wp,   &  !< parameter 95  - heat capacity 3rd wall layer roof
           1300000.0_wp,   &  !< parameter 96  - heat capacity 4th wall layer roof
           0.35_wp,        &  !< parameter 97  - thermal conductivity 1st/2nd wall layer roof
           0.8_wp,         &  !< parameter 98  - thermal conductivity 3rd wall layer roof
           2.1_wp,         &  !< parameter 99  - thermal conductivity 4th wall layer roof
           0.93_wp,        &  !< parameter 100 - wall emissivity roof
           27.0_wp,        &  !< parameter 101 - wall albedo roof
           0.0_wp,         &  !< parameter 102 - window fraction roof
           0.003_wp,       &  !< parameter 103 - window 1st layer thickness roof
           0.006_wp,       &  !< parameter 104 - window 2nd layer thickness roof
           0.012_wp,       &  !< parameter 105 - window 3rd layer thickness roof
           0.018_wp,       &  !< parameter 106 - window 4th layer thickness roof
           1736000.0_wp,   &  !< parameter 107 - heat capacity 1st/2nd window layer roof
           1736000.0_wp,   &  !< parameter 108 - heat capacity 3rd window layer roof
           1736000.0_wp,   &  !< parameter 109 - heat capacity 4th window layer roof
           0.57_wp,        &  !< parameter 110 - thermal conductivity 1st/2nd window layer roof
           0.57_wp,        &  !< parameter 111 - thermal conductivity 3rd window layer roof
           0.57_wp,        &  !< parameter 112 - thermal conductivity 4th window layer roof
           0.91_wp,        &  !< parameter 113 - window emissivity roof
           0.75_wp,        &  !< parameter 114 - window transmissivity roof
           27.0_wp,        &  !< parameter 115 - window albedo roof
           0.86_wp,        &  !< parameter 116 - green emissivity roof
           5.0_wp,         &  !< parameter 117 - green albedo roof
           0.0_wp,         &  !< parameter 118 - green type roof
           0.25_wp,        &  !< parameter 119 - shading factor
           0.76_wp,        &  !< parameter 120 - g-value windows
           5.0_wp,         &  !< parameter 121 - u-value windows
           0.1_wp,         &  !< parameter 122 - basical airflow without occupancy of the room for - summer 0.1_wp, winter 0.1
           1.5_wp,         &  !< parameter 123 - additional airflow depend of occupancy of the room for - summer 1.5_wp, winter 1.5
           0.0_wp,         &  !< parameter 124 - heat recovery efficiency
           3.5_wp,         &  !< parameter 125 - dynamic parameter specific effective surface
           370000.0_wp,    &  !< parameter 126 - dynamic parameter innner heatstorage
           4.5_wp,         &  !< parameter 127 - ratio internal surface/floor area
           100.0_wp,       &  !< parameter 128 - maximal heating capacity
           -200.0_wp,      &  !< parameter 129 - maximal cooling capacity
           3.0_wp,         &  !< parameter 130 - additional internal heat gains dependent on occupancy of the room
           10.0_wp,        &  !< parameter 131 - basic internal heat gains without occupancy of the room
           3.0_wp,         &  !< parameter 132 - storey height
           0.2_wp,         &  !< parameter 133 - ceiling construction height
           0.1_wp,         &  !< parameter 134 - anthropogenic heat output for heating
           1.333_wp        &  !< parameter 135 - anthropogenic heat output for cooling
                            /)   
                            
        building_pars(:,5) = (/   &
           0.5_wp,         &  !< parameter 0   - wall fraction above ground floor level
           0.5_wp,         &  !< parameter 1   - window fraction above ground floor level
           0.0_wp,         &  !< parameter 2   - green fraction above ground floor level
           0.0_wp,         &  !< parameter 3   - green fraction roof above ground floor level
           1.5_wp,         &  !< parameter 4   - LAI roof
           1.5_wp,         &  !< parameter 5   - LAI on wall above ground floor level
           2000000.0_wp,   &  !< parameter 6   - heat capacity 1st/2nd wall layer above ground floor level
           103000.0_wp,    &  !< parameter 7   - heat capacity 3rd wall layer above ground floor level
           900000.0_wp,    &  !< parameter 8   - heat capacity 4th wall layer above ground floor level
           0.35_wp,        &  !< parameter 9   - thermal conductivity 1st/2nd wall layer above ground floor level
           0.38_wp,        &  !< parameter 10  - thermal conductivity 3rd wall layer above ground floor level
           0.04_wp,        &  !< parameter 11  - thermal conductivity 4th wall layer above ground floor level
           299.15_wp,      &  !< parameter 12  - indoor target summer temperature
           293.15_wp,      &  !< parameter 13  - indoor target winter temperature
           0.92_wp,        &  !< parameter 14  - wall emissivity above ground floor level
           0.86_wp,        &  !< parameter 15  - green emissivity above ground floor level
           0.87_wp,        &  !< parameter 16  - window emissivity above ground floor level
           0.7_wp,         &  !< parameter 17  - window transmissivity above ground floor level
           0.001_wp,       &  !< parameter 18  - z0 roughness above ground floor level
           0.0001_wp,      &  !< parameter 19  - z0h/z0g roughness heat/humidity above ground floor level
           4.0_wp,         &  !< parameter 20  - ground floor level height
           0.55_wp,        &  !< parameter 21  - wall fraction ground floor level
           0.45_wp,        &  !< parameter 22  - window fraction ground floor level
           0.0_wp,         &  !< parameter 23  - green fraction ground floor level
           0.0_wp,         &  !< parameter 24  - green fraction roof ground floor level
           1.5_wp,         &  !< parameter 25  - LAI on wall ground floor level
           2000000.0_wp,   &  !< parameter 26  - heat capacity 1st/2nd wall layer ground floor level
           103000.0_wp,    &  !< parameter 27  - heat capacity 3rd wall layer ground floor level
           900000.0_wp,    &  !< parameter 28  - heat capacity 4th wall layer ground floor level
           0.35_wp,        &  !< parameter 29  - thermal conductivity 1st/2nd wall layer ground floor level
           0.38_wp,        &  !< parameter 30  - thermal conductivity 3rd wall layer ground floor level
           0.04_wp,        &  !< parameter 31  - thermal conductivity 4th wall layer ground floor level
           0.92_wp,        &  !< parameter 32  - wall emissivity ground floor level
           0.87_wp,        &  !< parameter 33  - window emissivity ground floor level
           0.86_wp,        &  !< parameter 34  - green emissivity ground floor level
           0.7_wp,         &  !< parameter 35  - window transmissivity ground floor level
           0.001_wp,       &  !< parameter 36  - z0 roughness ground floor level
           0.0001_wp,      &  !< parameter 37  - z0h/z0q roughness heat/humidity
           27.0_wp,        &  !< parameter 38  - wall albedo above ground floor level
           5.0_wp,         &  !< parameter 39  - green albedo above ground floor level
           27.0_wp,        &  !< parameter 40  - window albedo above ground floor level
           0.005_wp,       &  !< parameter 41  - 1st wall layer thickness above ground floor level
           0.01_wp,        &  !< parameter 42  - 2nd wall layer thickness above ground floor level
           0.31_wp,        &  !< parameter 43  - 3rd wall layer thickness above ground floor level
           0.43_wp,        &  !< parameter 44  - 4th wall layer thickness above ground floor level
           20000.0_wp,     &  !< parameter 45  - heat capacity wall surface
           23.0_wp,        &  !< parameter 46  - thermal conductivity of wall surface
           20000.0_wp,     &  !< parameter 47  - heat capacity of window surface
           20000.0_wp,     &  !< parameter 48  - heat capacity of green surface
           23.0_wp,        &  !< parameter 49  - thermal conductivity of window surface
           10.0_wp,        &  !< parameter 50  - thermal conductivty of green surface
           1.0_wp,         &  !< parameter 51  - wall fraction ground plate
           0.005_wp,       &  !< parameter 52  - 1st wall layer thickness ground plate
           0.01_wp,        &  !< parameter 53  - 2nd wall layer thickness ground plate
           0.31_wp,        &  !< parameter 54  - 3rd wall layer thickness ground plate
           0.43_wp,        &  !< parameter 55  - 4th wall layer thickness ground plate
           2000000.0_wp,   &  !< parameter 56  - heat capacity 1st/2nd wall layer ground plate
           103000.0_wp,    &  !< parameter 57  - heat capacity 3rd wall layer ground plate
           900000.0_wp,    &  !< parameter 58  - heat capacity 4th wall layer ground plate
           0.35_wp,        &  !< parameter 59  - thermal conductivity 1st/2nd wall layer ground plate
           0.38_wp,        &  !< parameter 60  - thermal conductivity 3rd wall layer ground plate
           0.04_wp,        &  !< parameter 61  - thermal conductivity 4th wall layer ground plate
           0.005_wp,       &  !< parameter 62  - 1st wall layer thickness ground floor level
           0.01_wp,        &  !< parameter 63  - 2nd wall layer thickness ground floor level
           0.31_wp,        &  !< parameter 64  - 3rd wall layer thickness ground floor level
           0.43_wp,        &  !< parameter 65  - 4th wall layer thickness ground floor level
           27.0_wp,        &  !< parameter 66  - wall albedo ground floor level
           0.003_wp,       &  !< parameter 67  - 1st window layer thickness ground floor level
           0.006_wp,       &  !< parameter 68  - 2nd window layer thickness ground floor level
           0.012_wp,       &  !< parameter 69  - 3rd window layer thickness ground floor level
           0.018_wp,       &  !< parameter 70  - 4th window layer thickness ground floor level
           1736000.0_wp,   &  !< parameter 71  - heat capacity 1st/2nd window layer ground floor level
           1736000.0_wp,   &  !< parameter 72  - heat capacity 3rd window layer ground floor level
           1736000.0_wp,   &  !< parameter 73  - heat capacity 4th window layer ground floor level
           0.11_wp,        &  !< parameter 74  - thermal conductivity 1st/2nd window layer ground floor level
           0.11_wp,        &  !< parameter 75  - thermal conductivity 3rd window layer ground floor level
           0.11_wp,        &  !< parameter 76  - thermal conductivity 4th window layer ground floor level
           27.0_wp,        &  !< parameter 77  - window albedo ground floor level
           5.0_wp,         &  !< parameter 78  - green albedo ground floor level
           0.003_wp,       &  !< parameter 79  - 1st window layer thickness above ground floor level
           0.006_wp,       &  !< parameter 80  - 2nd thickness window layer above ground floor level
           0.012_wp,       &  !< parameter 81  - 3rd window layer thickness above ground floor level
           0.018_wp,       &  !< parameter 82  - 4th window layer thickness above ground floor level
           1736000.0_wp,   &  !< parameter 83  - heat capacity 1st/2nd window layer above ground floor level
           1736000.0_wp,   &  !< parameter 84  - heat capacity 3rd window layer above ground floor level
           1736000.0_wp,   &  !< parameter 85  - heat capacity 4th window layer above ground floor level
           0.11_wp,        &  !< parameter 86  - thermal conductivity 1st/2nd window layer above ground floor level
           0.11_wp,        &  !< parameter 87  - thermal conductivity 3rd window layer above ground floor level
           0.11_wp,        &  !< parameter 88  - thermal conductivity 4th window layer above ground floor level
           1.0_wp,         &  !< parameter 89  - wall fraction roof
           0.005_wp,       &  !< parameter 90  - 1st wall layer thickness roof
           0.01_wp,        &  !< parameter 91  - 2nd wall layer thickness roof
           0.31_wp,        &  !< parameter 92  - 3rd wall layer thickness roof
           0.43_wp,        &  !< parameter 93  - 4th wall layer thickness roof
           2000000.0_wp,   &  !< parameter 94  - heat capacity 1st/2nd wall layer roof
           103000.0_wp,    &  !< parameter 95  - heat capacity 3rd wall layer roof
           900000.0_wp,    &  !< parameter 96  - heat capacity 4th wall layer roof
           0.35_wp,        &  !< parameter 97  - thermal conductivity 1st/2nd wall layer roof
           0.38_wp,        &  !< parameter 98  - thermal conductivity 3rd wall layer roof
           0.04_wp,        &  !< parameter 99  - thermal conductivity 4th wall layer roof
           0.91_wp,        &  !< parameter 100 - wall emissivity roof
           27.0_wp,        &  !< parameter 101 - wall albedo roof
           0.0_wp,         &  !< parameter 102 - window fraction roof
           0.003_wp,       &  !< parameter 103 - window 1st layer thickness roof
           0.006_wp,       &  !< parameter 104 - window 2nd layer thickness roof
           0.012_wp,       &  !< parameter 105 - window 3rd layer thickness roof
           0.018_wp,       &  !< parameter 106 - window 4th layer thickness roof
           1736000.0_wp,   &  !< parameter 107 - heat capacity 1st/2nd window layer roof
           1736000.0_wp,   &  !< parameter 108 - heat capacity 3rd window layer roof
           1736000.0_wp,   &  !< parameter 109 - heat capacity 4th window layer roof
           0.11_wp,        &  !< parameter 110 - thermal conductivity 1st/2nd window layer roof
           0.11_wp,        &  !< parameter 111 - thermal conductivity 3rd window layer roof
           0.11_wp,        &  !< parameter 112 - thermal conductivity 4th window layer roof
           0.87_wp,        &  !< parameter 113 - window emissivity roof
           0.7_wp,         &  !< parameter 114 - window transmissivity roof
           27.0_wp,        &  !< parameter 115 - window albedo roof
           0.86_wp,        &  !< parameter 116 - green emissivity roof
           5.0_wp,         &  !< parameter 117 - green albedo roof
           0.0_wp,         &  !< parameter 118 - green type roof
           0.25_wp,        &  !< parameter 119 - shading factor
           0.6_wp,         &  !< parameter 120 - g-value windows
           3.0_wp,         &  !< parameter 121 - u-value windows
           0.1_wp,         &  !< parameter 122 - basical airflow without occupancy of the room for - summer 0.1_wp, winter 0.1
           1.5_wp,         &  !< parameter 123 - additional airflow depend of occupancy of the room for - summer 1.5_wp, winter 1.5
           0.65_wp,        &  !< parameter 124 - heat recovery efficiency
           2.5_wp,         &  !< parameter 125 - dynamic parameter specific effective surface
           165000.0_wp,    &  !< parameter 126 - dynamic parameter innner heatstorage
           4.5_wp,         &  !< parameter 127 - ratio internal surface/floor area
           100.0_wp,       &  !< parameter 128 - maximal heating capacity
           -200.0_wp,      &  !< parameter 129 - maximal cooling capacity
           7.0_wp,         &  !< parameter 130 - additional internal heat gains dependent on occupancy of the room
           20.0_wp,        &  !< parameter 131 - basic internal heat gains without occupancy of the room
           3.0_wp,         &  !< parameter 132 - storey height
           0.2_wp,         &  !< parameter 133 - ceiling construction height
           0.0_wp,         &  !< parameter 134 - anthropogenic heat output for heating
           2.54_wp         &  !< parameter 135 - anthropogenic heat output for cooling
                            /)
                            
        building_pars(:,6) = (/   &
           0.425_wp,       &  !< parameter 0   - wall fraction above ground floor level
           0.575_wp,       &  !< parameter 1   - window fraction above ground floor level
           0.0_wp,         &  !< parameter 2   - green fraction above ground floor level
           0.0_wp,         &  !< parameter 3   - green fraction roof above ground floor level
           1.5_wp,         &  !< parameter 4   - LAI roof
           1.5_wp,         &  !< parameter 5   - LAI on wall above ground floor level
           2000000.0_wp,   &  !< parameter 6   - heat capacity 1st/2nd wall layer above ground floor level
           103000.0_wp,    &  !< parameter 7   - heat capacity 3rd wall layer above ground floor level
           900000.0_wp,    &  !< parameter 8   - heat capacity 4th wall layer above ground floor level
           0.35_wp,        &  !< parameter 9   - thermal conductivity 1st/2nd wall layer above ground floor level
           0.14_wp,        &  !< parameter 10  - thermal conductivity 3rd wall layer above ground floor level
           0.035_wp,       &  !< parameter 11  - thermal conductivity 4th wall layer above ground floor level
           299.15_wp,      &  !< parameter 12  - indoor target summer temperature
           293.15_wp,      &  !< parameter 13  - indoor target winter temperature
           0.92_wp,        &  !< parameter 14  - wall emissivity above ground floor level
           0.86_wp,        &  !< parameter 15  - green emissivity above ground floor level
           0.8_wp,         &  !< parameter 16  - window emissivity above ground floor level
           0.6_wp,         &  !< parameter 17  - window transmissivity above ground floor level
           0.001_wp,       &  !< parameter 18  - z0 roughness above ground floor level
           0.0001_wp,      &  !< parameter 19  - z0h/z0g roughness heat/humidity above ground floor level
           4.0_wp,         &  !< parameter 20  - ground floor level height
           0.475_wp,       &  !< parameter 21  - wall fraction ground floor level
           0.525_wp,       &  !< parameter 22  - window fraction ground floor level
           0.0_wp,         &  !< parameter 23  - green fraction ground floor level
           0.0_wp,         &  !< parameter 24  - green fraction roof ground floor level
           1.5_wp,         &  !< parameter 25  - LAI on wall ground floor level
           2000000.0_wp,   &  !< parameter 26  - heat capacity 1st/2nd wall layer ground floor level
           103000.0_wp,    &  !< parameter 27  - heat capacity 3rd wall layer ground floor level
           900000.0_wp,    &  !< parameter 28  - heat capacity 4th wall layer ground floor level
           0.35_wp,        &  !< parameter 29  - thermal conductivity 1st/2nd wall layer ground floor level
           0.14_wp,        &  !< parameter 30  - thermal conductivity 3rd wall layer ground floor level
           0.035_wp,       &  !< parameter 31  - thermal conductivity 4th wall layer ground floor level
           0.92_wp,        &  !< parameter 32  - wall emissivity ground floor level
           0.8_wp,         &  !< parameter 33  - window emissivity ground floor level
           0.86_wp,        &  !< parameter 34  - green emissivity ground floor level
           0.6_wp,         &  !< parameter 35  - window transmissivity ground floor level
           0.001_wp,       &  !< parameter 36  - z0 roughness ground floor level
           0.0001_wp,      &  !< parameter 37  - z0h/z0q roughness heat/humidity
           27.0_wp,        &  !< parameter 38  - wall albedo above ground floor level
           5.0_wp,         &  !< parameter 39  - green albedo above ground floor level
           27.0_wp,        &  !< parameter 40  - window albedo above ground floor level
           0.005_wp,       &  !< parameter 41  - 1st wall layer thickness above ground floor level
           0.01_wp,        &  !< parameter 42  - 2nd wall layer thickness above ground floor level
           0.41_wp,        &  !< parameter 43  - 3rd wall layer thickness above ground floor level
           0.7_wp,         &  !< parameter 44  - 4th wall layer thickness above ground floor level
           20000.0_wp,     &  !< parameter 45  - heat capacity wall surface
           23.0_wp,        &  !< parameter 46  - thermal conductivity of wall surface
           20000.0_wp,     &  !< parameter 47  - heat capacity of window surface
           20000.0_wp,     &  !< parameter 48  - heat capacity of green surface
           23.0_wp,        &  !< parameter 49  - thermal conductivity of window surface
           10.0_wp,        &  !< parameter 50  - thermal conductivty of green surface
           1.0_wp,         &  !< parameter 51  - wall fraction ground plate
           0.005_wp,       &  !< parameter 52  - 1st wall layer thickness ground plate
           0.01_wp,        &  !< parameter 53  - 2nd wall layer thickness ground plate
           0.41_wp,        &  !< parameter 54  - 3rd wall layer thickness ground plate
           0.7_wp,         &  !< parameter 55  - 4th wall layer thickness ground plate
           2000000.0_wp,   &  !< parameter 56  - heat capacity 1st/2nd wall layer ground plate
           103000.0_wp,    &  !< parameter 57  - heat capacity 3rd wall layer ground plate
           900000.0_wp,    &  !< parameter 58  - heat capacity 4th wall layer ground plate
           0.35_wp,        &  !< parameter 59  - thermal conductivity 1st/2nd wall layer ground plate
           0.14_wp,        &  !< parameter 60  - thermal conductivity 3rd wall layer ground plate
           0.035_wp,       &  !< parameter 61  - thermal conductivity 4th wall layer ground plate
           0.005_wp,       &  !< parameter 62  - 1st wall layer thickness ground floor level
           0.01_wp,        &  !< parameter 63  - 2nd wall layer thickness ground floor level
           0.41_wp,        &  !< parameter 64  - 3rd wall layer thickness ground floor level
           0.7_wp,         &  !< parameter 65  - 4th wall layer thickness ground floor level
           27.0_wp,        &  !< parameter 66  - wall albedo ground floor level
           0.003_wp,       &  !< parameter 67  - 1st window layer thickness ground floor level
           0.006_wp,       &  !< parameter 68  - 2nd window layer thickness ground floor level
           0.012_wp,       &  !< parameter 69  - 3rd window layer thickness ground floor level
           0.018_wp,       &  !< parameter 70  - 4th window layer thickness ground floor level
           1736000.0_wp,   &  !< parameter 71  - heat capacity 1st/2nd window layer ground floor level
           1736000.0_wp,   &  !< parameter 72  - heat capacity 3rd window layer ground floor level
           1736000.0_wp,   &  !< parameter 73  - heat capacity 4th window layer ground floor level
           0.037_wp,       &  !< parameter 74  - thermal conductivity 1st/2nd window layer ground floor level
           0.037_wp,       &  !< parameter 75  - thermal conductivity 3rd window layer ground floor level
           0.037_wp,       &  !< parameter 76  - thermal conductivity 4th window layer ground floor level
           27.0_wp,        &  !< parameter 77  - window albedo ground floor level
           5.0_wp,         &  !< parameter 78  - green albedo ground floor level
           0.003_wp,       &  !< parameter 79  - 1st window layer thickness above ground floor level
           0.006_wp,       &  !< parameter 80  - 2nd thickness window layer above ground floor level
           0.012_wp,       &  !< parameter 81  - 3rd window layer thickness above ground floor level
           0.018_wp,       &  !< parameter 82  - 4th window layer thickness above ground floor level
           1736000.0_wp,   &  !< parameter 83  - heat capacity 1st/2nd window layer above ground floor level
           1736000.0_wp,   &  !< parameter 84  - heat capacity 3rd window layer above ground floor level
           1736000.0_wp,   &  !< parameter 85  - heat capacity 4th window layer above ground floor level
           0.037_wp,       &  !< parameter 86  - thermal conductivity 1st/2nd window layer above ground floor level
           0.037_wp,       &  !< parameter 87  - thermal conductivity 3rd window layer above ground floor level
           0.037_wp,       &  !< parameter 88  - thermal conductivity 4th window layer above ground floor level
           1.0_wp,         &  !< parameter 89  - wall fraction roof
           0.005_wp,       &  !< parameter 90  - 1st wall layer thickness roof
           0.01_wp,        &  !< parameter 91  - 2nd wall layer thickness roof
           0.41_wp,        &  !< parameter 92  - 3rd wall layer thickness roof
           0.7_wp,         &  !< parameter 93  - 4th wall layer thickness roof
           2000000.0_wp,   &  !< parameter 94  - heat capacity 1st/2nd wall layer roof
           103000.0_wp,    &  !< parameter 95  - heat capacity 3rd wall layer roof
           900000.0_wp,    &  !< parameter 96  - heat capacity 4th wall layer roof
           0.35_wp,        &  !< parameter 97  - thermal conductivity 1st/2nd wall layer roof
           0.14_wp,        &  !< parameter 98  - thermal conductivity 3rd wall layer roof
           0.035_wp,       &  !< parameter 99  - thermal conductivity 4th wall layer roof
           0.91_wp,        &  !< parameter 100 - wall emissivity roof
           27.0_wp,        &  !< parameter 101 - wall albedo roof
           0.0_wp,         &  !< parameter 102 - window fraction roof
           0.003_wp,       &  !< parameter 103 - window 1st layer thickness roof
           0.006_wp,       &  !< parameter 104 - window 2nd layer thickness roof
           0.012_wp,       &  !< parameter 105 - window 3rd layer thickness roof
           0.018_wp,       &  !< parameter 106 - window 4th layer thickness roof
           1736000.0_wp,   &  !< parameter 107 - heat capacity 1st/2nd window layer roof
           1736000.0_wp,   &  !< parameter 108 - heat capacity 3rd window layer roof
           1736000.0_wp,   &  !< parameter 109 - heat capacity 4th window layer roof
           0.037_wp,       &  !< parameter 110 - thermal conductivity 1st/2nd window layer roof
           0.037_wp,       &  !< parameter 111 - thermal conductivity 3rd window layer roof
           0.037_wp,       &  !< parameter 112 - thermal conductivity 4th window layer roof
           0.8_wp,         &  !< parameter 113 - window emissivity roof
           0.6_wp,         &  !< parameter 114 - window transmissivity roof
           27.0_wp,        &  !< parameter 115 - window albedo roof
           0.86_wp,        &  !< parameter 116 - green emissivity roof
           5.0_wp,         &  !< parameter 117 - green albedo roof
           0.0_wp,         &  !< parameter 118 - green type roof
           0.25_wp,        &  !< parameter 119 - shading factor
           0.5_wp,         &  !< parameter 120 - g-value windows
           2.5_wp,         &  !< parameter 121 - u-value windows
           0.1_wp,         &  !< parameter 122 - basical airflow without occupancy of the room for - summer 0.1_wp, winter 0.1
           1.5_wp,         &  !< parameter 123 - additional airflow depend of occupancy of the room for - summer 1.5_wp, winter 1.5
           0.9_wp,         &  !< parameter 124 - heat recovery efficiency
           2.5_wp,         &  !< parameter 125 - dynamic parameter specific effective surface
           80000.0_wp,     &  !< parameter 126 - dynamic parameter innner heatstorage
           4.5_wp,         &  !< parameter 127 - ratio internal surface/floor area
           100.0_wp,       &  !< parameter 128 - maximal heating capacity
           -80.0_wp,       &  !< parameter 129 - maximal cooling capacity
           5.0_wp,         &  !< parameter 130 - additional internal heat gains dependent on occupancy of the room
           15.0_wp,        &  !< parameter 131 - basic internal heat gains without occupancy of the room
           3.0_wp,         &  !< parameter 132 - storey height
           0.2_wp,         &  !< parameter 133 - ceiling construction height
           -2.0_wp,        &  !< parameter 134 - anthropogenic heat output for heating
           1.25_wp         &  !< parameter 135 - anthropogenic heat output for cooling
                            /)
                            
        building_pars(:,7) = (/   &
           1.0_wp,         &  !< parameter 0   - wall fraction above ground floor level
           0.0_wp,         &  !< parameter 1   - window fraction above ground floor level
           0.0_wp,         &  !< parameter 2   - green fraction above ground floor level
           0.0_wp,         &  !< parameter 3   - green fraction roof above ground floor level
           1.5_wp,         &  !< parameter 4   - LAI roof
           1.5_wp,         &  !< parameter 5   - LAI on wall above ground floor level
           1950400.0_wp,   &  !< parameter 6   - heat capacity 1st/2nd wall layer above ground floor level
           1848000.0_wp,   &  !< parameter 7   - heat capacity 3rd wall layer above ground floor level
           1848000.0_wp,   &  !< parameter 8   - heat capacity 4th wall layer above ground floor level
           0.7_wp,         &  !< parameter 9   - thermal conductivity 1st/2nd wall layer above ground floor level
           1.0_wp,         &  !< parameter 10  - thermal conductivity 3rd wall layer above ground floor level
           1.0_wp,         &  !< parameter 11  - thermal conductivity 4th wall layer above ground floor level
           299.15_wp,      &  !< parameter 12  - indoor target summer temperature
           293.15_wp,      &  !< parameter 13  - indoor target winter temperature
           0.9_wp,         &  !< parameter 14  - wall emissivity above ground floor level
           0.86_wp,        &  !< parameter 15  - green emissivity above ground floor level
           0.8_wp,         &  !< parameter 16  - window emissivity above ground floor level
           0.6_wp,         &  !< parameter 17  - window transmissivity above ground floor level
           0.001_wp,       &  !< parameter 18  - z0 roughness above ground floor level
           0.0001_wp,      &  !< parameter 19  - z0h/z0g roughness heat/humidity above ground floor level
           4.0_wp,         &  !< parameter 20  - ground floor level height
           1.0_wp,         &  !< parameter 21  - wall fraction ground floor level
           0.0_wp,         &  !< parameter 22  - window fraction ground floor level
           0.0_wp,         &  !< parameter 23  - green fraction ground floor level
           0.0_wp,         &  !< parameter 24  - green fraction roof ground floor level
           1.5_wp,         &  !< parameter 25  - LAI on wall ground floor level
           1950400.0_wp,   &  !< parameter 26  - heat capacity 1st/2nd wall layer ground floor level
           1848000.0_wp,   &  !< parameter 27  - heat capacity 3rd wall layer ground floor level
           1848000.0_wp,   &  !< parameter 28  - heat capacity 4th wall layer ground floor level
           0.7_wp,         &  !< parameter 29  - thermal conductivity 1st/2nd wall layer ground floor level
           1.0_wp,         &  !< parameter 30  - thermal conductivity 3rd wall layer ground floor level
           1.0_wp,         &  !< parameter 31  - thermal conductivity 4th wall layer ground floor level
           0.9_wp,         &  !< parameter 32  - wall emissivity ground floor level
           0.8_wp,         &  !< parameter 33  - window emissivity ground floor level
           0.86_wp,        &  !< parameter 34  - green emissivity ground floor level
           0.6_wp,         &  !< parameter 35  - window transmissivity ground floor level
           0.001_wp,       &  !< parameter 36  - z0 roughness ground floor level
           0.0001_wp,      &  !< parameter 37  - z0h/z0q roughness heat/humidity
           27.0_wp,        &  !< parameter 38  - wall albedo above ground floor level
           5.0_wp,         &  !< parameter 39  - green albedo above ground floor level
           27.0_wp,        &  !< parameter 40  - window albedo above ground floor level
           0.29_wp,        &  !< parameter 41  - 1st wall layer thickness above ground floor level
           0.295_wp,       &  !< parameter 42  - 2nd wall layer thickness above ground floor level
           0.695_wp,       &  !< parameter 43  - 3rd wall layer thickness above ground floor level
           0.985_wp,       &  !< parameter 44  - 4th wall layer thickness above ground floor level
           20000.0_wp,     &  !< parameter 45  - heat capacity wall surface
           23.0_wp,        &  !< parameter 46  - thermal conductivity of wall surface
           20000.0_wp,     &  !< parameter 47  - heat capacity of window surface
           20000.0_wp,     &  !< parameter 48  - heat capacity of green surface
           23.0_wp,        &  !< parameter 49  - thermal conductivity of window surface
           10.0_wp,        &  !< parameter 50  - thermal conductivty of green surface
           1.0_wp,         &  !< parameter 51  - wall fraction ground plate
           0.29_wp,        &  !< parameter 52  - 1st wall layer thickness ground plate
           0.295_wp,       &  !< parameter 53  - 2nd wall layer thickness ground plate
           0.695_wp,       &  !< parameter 54  - 3rd wall layer thickness ground plate
           0.985_wp,       &  !< parameter 55  - 4th wall layer thickness ground plate
           1950400.0_wp,   &  !< parameter 56  - heat capacity 1st/2nd wall layer ground plate
           1848000.0_wp,   &  !< parameter 57  - heat capacity 3rd wall layer ground plate
           1848000.0_wp,   &  !< parameter 58  - heat capacity 4th wall layer ground plate
           0.7_wp,         &  !< parameter 59  - thermal conductivity 1st/2nd wall layer ground plate
           1.0_wp,         &  !< parameter 60  - thermal conductivity 3rd wall layer ground plate
           1.0_wp,         &  !< parameter 61  - thermal conductivity 4th wall layer ground plate
           0.29_wp,        &  !< parameter 62  - 1st wall layer thickness ground floor level
           0.295_wp,       &  !< parameter 63  - 2nd wall layer thickness ground floor level
           0.695_wp,       &  !< parameter 64  - 3rd wall layer thickness ground floor level
           0.985_wp,       &  !< parameter 65  - 4th wall layer thickness ground floor level
           27.0_wp,        &  !< parameter 66  - wall albedo ground floor level
           0.003_wp,       &  !< parameter 67  - 1st window layer thickness ground floor level
           0.006_wp,       &  !< parameter 68  - 2nd window layer thickness ground floor level
           0.012_wp,       &  !< parameter 69  - 3rd window layer thickness ground floor level
           0.018_wp,       &  !< parameter 70  - 4th window layer thickness ground floor level
           1736000.0_wp,   &  !< parameter 71  - heat capacity 1st/2nd window layer ground floor level
           1736000.0_wp,   &  !< parameter 72  - heat capacity 3rd window layer ground floor level
           1736000.0_wp,   &  !< parameter 73  - heat capacity 4th window layer ground floor level
           0.57_wp,        &  !< parameter 74  - thermal conductivity 1st/2nd window layer ground floor level
           0.57_wp,        &  !< parameter 75  - thermal conductivity 3rd window layer ground floor level
           0.57_wp,        &  !< parameter 76  - thermal conductivity 4th window layer ground floor level
           27.0_wp,        &  !< parameter 77  - window albedo ground floor level
           5.0_wp,         &  !< parameter 78  - green albedo ground floor level
           0.003_wp,       &  !< parameter 79  - 1st window layer thickness above ground floor level
           0.006_wp,       &  !< parameter 80  - 2nd thickness window layer above ground floor level
           0.012_wp,       &  !< parameter 81  - 3rd window layer thickness above ground floor level
           0.018_wp,       &  !< parameter 82  - 4th window layer thickness above ground floor level
           1736000.0_wp,   &  !< parameter 83  - heat capacity 1st/2nd window layer above ground floor level
           1736000.0_wp,   &  !< parameter 84  - heat capacity 3rd window layer above ground floor level
           1736000.0_wp,   &  !< parameter 85  - heat capacity 4th window layer above ground floor level
           0.57_wp,        &  !< parameter 86  - thermal conductivity 1st/2nd window layer above ground floor level
           0.57_wp,        &  !< parameter 87  - thermal conductivity 3rd window layer above ground floor level
           0.57_wp,        &  !< parameter 88  - thermal conductivity 4th window layer above ground floor level
           1.0_wp,         &  !< parameter 89  - wall fraction roof
           0.29_wp,        &  !< parameter 90  - 1st wall layer thickness roof
           0.295_wp,       &  !< parameter 91  - 2nd wall layer thickness roof
           0.695_wp,       &  !< parameter 92  - 3rd wall layer thickness roof
           0.985_wp,       &  !< parameter 93  - 4th wall layer thickness roof
           1950400.0_wp,   &  !< parameter 94  - heat capacity 1st/2nd wall layer roof
           1848000.0_wp,   &  !< parameter 95  - heat capacity 3rd wall layer roof
           1848000.0_wp,   &  !< parameter 96  - heat capacity 4th wall layer roof
           0.7_wp,         &  !< parameter 97  - thermal conductivity 1st/2nd wall layer roof
           1.0_wp,         &  !< parameter 98  - thermal conductivity 3rd wall layer roof
           1.0_wp,         &  !< parameter 99  - thermal conductivity 4th wall layer roof
           0.9_wp,         &  !< parameter 100 - wall emissivity roof
           27.0_wp,        &  !< parameter 101 - wall albedo roof
           0.0_wp,         &  !< parameter 102 - window fraction roof
           0.003_wp,       &  !< parameter 103 - window 1st layer thickness roof
           0.006_wp,       &  !< parameter 104 - window 2nd layer thickness roof
           0.012_wp,       &  !< parameter 105 - window 3rd layer thickness roof
           0.018_wp,       &  !< parameter 106 - window 4th layer thickness roof
           1736000.0_wp,   &  !< parameter 107 - heat capacity 1st/2nd window layer roof
           1736000.0_wp,   &  !< parameter 108 - heat capacity 3rd window layer roof
           1736000.0_wp,   &  !< parameter 109 - heat capacity 4th window layer roof
           0.57_wp,        &  !< parameter 110 - thermal conductivity 1st/2nd window layer roof
           0.57_wp,        &  !< parameter 111 - thermal conductivity 3rd window layer roof
           0.57_wp,        &  !< parameter 112 - thermal conductivity 4th window layer roof
           0.8_wp,         &  !< parameter 113 - window emissivity roof
           0.6_wp,         &  !< parameter 114 - window transmissivity roof
           27.0_wp,        &  !< parameter 115 - window albedo roof
           0.86_wp,        &  !< parameter 116 - green emissivity roof
           5.0_wp,         &  !< parameter 117 - green albedo roof
           0.0_wp,         &  !< parameter 118 - green type roof
           0.8_wp,         &  !< parameter 119 - shading factor
           100.0_wp,       &  !< parameter 120 - g-value windows
           100.0_wp,       &  !< parameter 121 - u-value windows
           20.0_wp,        &  !< parameter 122 - basical airflow without occupancy of the room
           20.0_wp,        &  !< parameter 123 - additional airflow depend of occupancy of the room
           0.0_wp,         &  !< parameter 124 - heat recovery efficiency
           1.0_wp,         &  !< parameter 125 - dynamic parameter specific effective surface
           1.0_wp,         &  !< parameter 126 - dynamic parameter innner heatstorage
           4.5_wp,         &  !< parameter 127 - ratio internal surface/floor area
           100000.0_wp,    &  !< parameter 128 - maximal heating capacity
           0.0_wp,         &  !< parameter 129 - maximal cooling capacity
           0.0_wp,         &  !< parameter 130 - additional internal heat gains dependent on occupancy of the room
           0.0_wp,         &  !< parameter 131 - basic internal heat gains without occupancy of the room
           3.0_wp,         &  !< parameter 132 - storey height
           0.2_wp,         &  !< parameter 133 - ceiling construction height
           0.0_wp,         &  !< parameter 134 - anthropogenic heat output for heating
           0.0_wp          &  !< parameter 135 - anthropogenic heat output for cooling
                        /)
                        
     END SUBROUTINE usm_define_pars
 
   
  END MODULE urban_surface_mod
