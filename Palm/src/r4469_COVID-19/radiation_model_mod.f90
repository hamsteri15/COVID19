!> @file radiation_model_mod.f90
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
! Copyright 2015-2020 Institute of Computer Science of the
!                     Czech Academy of Sciences, Prague
! Copyright 2015-2019 Czech Technical University in Prague
! Copyright 1997-2020 Leibniz Universitaet Hannover
!------------------------------------------------------------------------------!
!
! Current revisions:
! ------------------
! 
! 
! Former revisions:
! -----------------
! $Id: radiation_model_mod.f90 4457 2020-03-11 14:20:43Z raasch $
! use statement for exchange horiz added
! 
! 4452 2020-03-10 20:15:32Z suehring
! Bugfix in calc_albedo
!
! 4442 2020-03-04 19:21:13Z suehring
! - Change order of dimension in surface arrays %frac, %emissivity and %albedo
!   to allow for better vectorization in the radiation interactions.
! - Minor formatting issues
!
! 4441 2020-03-04 19:20:35Z suehring
! bugfixes: cpp-directives for serial mode moved, small changes to get serial mode compiled
!
! 4400 2020-02-10 20:32:41Z suehring
! Initialize radiation arrays with zero
!
! 4392 2020-01-31 16:14:57Z pavelkrc
! - Add debug tracing of large radiative fluxes (option trace_fluxes_above)
! - Print exact counts of SVF and CSF if debut_output is enabled
! - Update descriptions of RTM 3.0 and related comments
!
! 4360 2020-01-07 11:25:50Z suehring
! Renamed pc_heating_rate, pc_transpiration_rate, pc_transpiration_rate to
! pcm_heating_rate, pcm_latent_rate, pcm_transpiration_rate
!
! 4340 2019-12-16 08:17:03Z Giersch
! Albedo indices for building_surface_pars are now declared as parameters to
! prevent an error if the gfortran compiler with -Werror=unused-value is used
!
! 4291 2019-11-11 12:36:54Z moh.hefny
! Enabled RTM in case of biometeorology even if there is no vertical
! surfaces or 3D vegetation in the domain
!
! 4286 2019-10-30 16:01:14Z resler
! - Fix wrong treating of time_rad during interpolation in radiation_model_mod
! - Fix wrong checks of time_rad from dynamic driver in radiation_model_mod
! - Add new directional model of human body for MRT: ellipsoid
!
! 4271 2019-10-23 10:46:41Z maronga
! Bugfix: missing parentheses in calculation of snow albedo
!
! 4245 2019-09-30 08:40:37Z pavelkrc
! Initialize explicit per-surface albedos from building_surface_pars
!
! 4238 2019-09-25 16:06:01Z suehring
! Modify check in order to avoid equality comparisons of floating points
!
! 4227 2019-09-10 18:04:34Z gronemeier
! implement new palm_date_time_mod
!
! 4226 2019-09-10 17:03:24Z suehring
! - Netcdf input routine for dimension length renamed
! - Define time variable for external radiation input relative to time_utc_init
!
! 4210 2019-09-02 13:07:09Z suehring
! - Revise steering of splitting diffuse and direct radiation
! - Bugfixes in checks
! - Optimize mapping of radiation components onto 2D arrays, avoid unnecessary
!   operations
!
! 4208 2019-09-02 09:01:07Z suehring
! Bugfix in accessing albedo_pars in the clear-sky branch
! (merge from branch resler)
!
! 4198 2019-08-29 15:17:48Z gronemeier
! Prohibit execution of radiation model if rotation_angle is not zero
!
! 4197 2019-08-29 14:33:32Z suehring
! Revise steering of surface albedo initialization when albedo_pars is provided
!
! 4190 2019-08-27 15:42:37Z suehring
! Implement external radiation forcing also for level-of-detail = 2
! (horizontally 2D radiation)
!
! 4188 2019-08-26 14:15:47Z suehring
! Minor adjustment in error message
!
! 4187 2019-08-26 12:43:15Z suehring
! - Take external radiation from root domain dynamic input if not provided for
!   each nested domain
! - Combine MPI_ALLREDUCE calls to reduce mpi overhead
!
! 4182 2019-08-22 15:20:23Z scharf
! Corrected "Former revisions" section
!
! 4179 2019-08-21 11:16:12Z suehring
! Remove debug prints
!
! 4178 2019-08-21 11:13:06Z suehring
! External radiation forcing implemented.
!
! 4168 2019-08-16 13:50:17Z suehring
! Replace function get_topography_top_index by topo_top_ind
!
! 4157 2019-08-14 09:19:12Z suehring
! Give informative message on raytracing distance only by core zero
!
! 4148 2019-08-08 11:26:00Z suehring
! Comments added
!
! 4134 2019-08-02 18:39:57Z suehring
! Bugfix in formatted write statement
!
! 4127 2019-07-30 14:47:10Z suehring
! Remove unused pch_index (merge from branch resler)
!
! 4089 2019-07-11 14:30:27Z suehring
! - Correct level 2 initialization of spectral albedos in rrtmg branch, long- and
!   shortwave albedos were mixed-up.
! - Change order of albedo_pars so that it is now consistent with the defined
!   order of albedo_pars in PIDS
!
! 4069 2019-07-01 14:05:51Z Giersch
! Masked output running index mid has been introduced as a local variable to
! avoid runtime error (Loop variable has been modified) in time_integration
!
! 4067 2019-07-01 13:29:25Z suehring
! Bugfix, pass dummy string to MPI_INFO_SET (J. Resler)
!
! 4039 2019-06-18 10:32:41Z suehring
! Bugfix for masked data output
!
! 4008 2019-05-30 09:50:11Z moh.hefny
! Bugfix in check variable when a variable's string is less than 3
! characters is processed. All variables now are checked if they
! belong to radiation
!
! 3992 2019-05-22 16:49:38Z suehring
! Bugfix in rrtmg radiation branch in a nested run when the lowest prognistic
! grid points in a child domain are all inside topography
!
! 3987 2019-05-22 09:52:13Z kanani
! Introduce alternative switch for debug output during timestepping
!
! 3943 2019-05-02 09:50:41Z maronga
! Missing blank characteer added.
!
! 3900 2019-04-16 15:17:43Z suehring
! Fixed initialization problem
!
! 3885 2019-04-11 11:29:34Z kanani
! Changes related to global restructuring of location messages and introduction
! of additional debug messages
!
! 3881 2019-04-10 09:31:22Z suehring
! Output of albedo and emissivity moved from USM, bugfixes in initialization
! of albedo
!
! 3861 2019-04-04 06:27:41Z maronga
! Bugfix: factor of 4.0 required instead of 3.0 in calculation of rad_lw_out_change_0
!
! 3859 2019-04-03 20:30:31Z maronga
! Added some descriptions
!
! 3847 2019-04-01 14:51:44Z suehring
! Implement check for dt_radiation (must be > 0)
!
! 3846 2019-04-01 13:55:30Z suehring
! unused variable removed
!
! 3814 2019-03-26 08:40:31Z pavelkrc
! Change zenith(0:0) and others to scalar.
! Code review.
! Rename exported nzu, nzp and related variables due to name conflict
!
! 3771 2019-02-28 12:19:33Z raasch
! rrtmg preprocessor for directives moved/added, save attribute added to temporary
! pointers to avoid compiler warnings about outlived pointer targets,
! statement added to avoid compiler warning about unused variable
!
! 3769 2019-02-28 10:16:49Z moh.hefny
! removed unused variables and subroutine radiation_radflux_gridbox
!
! 3767 2019-02-27 08:18:02Z raasch
! unused variable for file index removed from rrd-subroutines parameter list
!
! 3760 2019-02-21 18:47:35Z moh.hefny
! Bugfix: initialized simulated_time before calculating solar position
! to enable restart option with reading in SVF from file(s).
!
! 3754 2019-02-19 17:02:26Z kanani
! (resler, pavelkrc)
! Bugfixes: add further required MRT factors to read/write_svf,
! fix for aggregating view factors to eliminate local noise in reflected
! irradiance at mutually close surfaces (corners, presence of trees) in the
! angular discretization scheme.
!
! 3752 2019-02-19 09:37:22Z resler
! added read/write number of MRT factors to the respective routines
!
! 3705 2019-01-29 19:56:39Z suehring
! Make variables that are sampled in virtual measurement module public
!
! 3704 2019-01-29 19:51:41Z suehring
! Some interface calls moved to module_interface + cleanup
!
! 3667 2019-01-10 14:26:24Z schwenkel
! Modified check for rrtmg input files
!
! 3655 2019-01-07 16:51:22Z knoop
! nopointer option removed
!
! 1496 2014-12-02 17:25:50Z maronga
! Initial revision
!
!
! Description:
! ------------
!> Radiation models and interfaces:
!> constant, simple and RRTMG models, interface to external radiation model
!> Radiative Transfer Model (RTM) version 3.0 for modelling of radiation
!> interactions within urban canopy or other surface layer in complex terrain
!> Integrations of RTM with other PALM-4U modules:
!> integration with RRTMG, USM, LSM, PCM, BIO modules
!>
!> @todo move variable definitions used in radiation_init only to the subroutine
!>       as they are no longer required after initialization.
!> @todo Output of full column vertical profiles used in RRTMG
!> @todo Output of other rrtm arrays (such as volume mixing ratios)
!> @todo Optimize radiation_tendency routines
!>
!> @note Many variables have a leading dummy dimension (0:0) in order to
!>       match the assume-size shape expected by the RRTMG model.
!------------------------------------------------------------------------------!
 MODULE radiation_model_mod

    USE arrays_3d,                                                             &
        ONLY:  dzw, hyp, nc, pt, p, q, ql, u, v, w, zu, zw, exner, d_exner

    USE basic_constants_and_equations_mod,                                     &
        ONLY:  c_p, g, lv_d_cp, l_v, pi, r_d, rho_l, solar_constant, sigma_sb, &
               barometric_formula

    USE calc_mean_profile_mod,                                                 &
        ONLY:  calc_mean_profile

    USE control_parameters,                                                    &
        ONLY:  biometeorology, cloud_droplets, coupling_char,                  &
               debug_output, debug_output_timestep, debug_string,              &
               dt_3d,                                                          &
               dz, dt_spinup, end_time,                                        &
               humidity,                                                       &
               initializing_actions, io_blocks, io_group,                      &
               land_surface, large_scale_forcing,                              &
               latitude, longitude, lsf_surf,                                  &
               message_string, plant_canopy, pt_surface,                       &
               rho_surface, simulated_time, spinup_time, surface_pressure,     &
               read_svf, write_svf,                                            &
               time_since_reference_point, urban_surface, varnamelength

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point, log_point_s

    USE grid_variables,                                                        &
         ONLY:  ddx, ddy, dx, dy

    USE indices,                                                               &
        ONLY:  nnx, nny, nx, nxl, nxlg, nxr, nxrg, ny, nyn, nyng, nys, nysg,   &
               nzb, nzt, topo_top_ind

    USE, INTRINSIC :: iso_c_binding

    USE kinds

    USE bulk_cloud_model_mod,                                                  &
        ONLY:  bulk_cloud_model, microphysics_morrison, na_init, nc_const, sigma_gc

#if defined ( __netcdf )
    USE NETCDF
#endif

    USE netcdf_data_input_mod,                                                 &
        ONLY:  albedo_type_f,                                                  &
               albedo_pars_f,                                                  &
               building_type_f,                                                &
               building_surface_pars_f,                                        &
               pavement_type_f,                                                &
               vegetation_type_f,                                              &
               water_type_f,                                                   &
               char_fill,                                                      &
               char_lod,                                                       &
               check_existence,                                                &
               close_input_file,                                               &
               get_attribute,                                                  &
               get_dimension_length,                                           &
               get_variable,                                                   &
               inquire_num_variables,                                          &
               inquire_variable_names,                                         &
               input_file_dynamic,                                             &
               input_pids_dynamic,                                             &
               num_var_pids,                                                   &
               pids_id,                                                        &
               open_read_file,                                                 &
               real_1d_3d,                                                     &
               vars_pids

    USE palm_date_time_mod,                                                    &
        ONLY:  date_time_str_len, get_date_time,                               &
               hours_per_day, seconds_per_hour

    USE plant_canopy_model_mod,                                                &
        ONLY:  lad_s,                                                          &
               pcm_heating_rate,                                               &
               pcm_transpiration_rate,                                         &
               pcm_latent_rate,                                                &
               plant_canopy_transpiration,                                     &
               pcm_calc_transpiration_rate

    USE pegrid

#if defined ( __rrtmg )
    USE parrrsw,                                                               &
        ONLY:  naerec, nbndsw

    USE parrrtm,                                                               &
        ONLY:  nbndlw

    USE rrtmg_lw_init,                                                         &
        ONLY:  rrtmg_lw_ini

    USE rrtmg_sw_init,                                                         &
        ONLY:  rrtmg_sw_ini

    USE rrtmg_lw_rad,                                                          &
        ONLY:  rrtmg_lw

    USE rrtmg_sw_rad,                                                          &
        ONLY:  rrtmg_sw
#endif
    USE statistics,                                                            &
        ONLY:  hom

    USE surface_mod,                                                           &
        ONLY:  ind_pav_green, ind_veg_wall, ind_wat_win,                       &
               surf_lsm_h, surf_lsm_v, surf_type, surf_usm_h, surf_usm_v,      &
               vertical_surfaces_exist

    IMPLICIT NONE

    CHARACTER(10) :: radiation_scheme = 'clear-sky' ! 'constant', 'clear-sky', or 'rrtmg'

!
!-- Predefined Land surface classes (albedo_type) after Briegleb (1992)
    CHARACTER(37), DIMENSION(0:33), PARAMETER :: albedo_type_name = (/      &
                                   'user defined                         ', & !  0
                                   'ocean                                ', & !  1
                                   'mixed farming, tall grassland        ', & !  2
                                   'tall/medium grassland                ', & !  3
                                   'evergreen shrubland                  ', & !  4
                                   'short grassland/meadow/shrubland     ', & !  5
                                   'evergreen needleleaf forest          ', & !  6
                                   'mixed deciduous evergreen forest     ', & !  7
                                   'deciduous forest                     ', & !  8
                                   'tropical evergreen broadleaved forest', & !  9
                                   'medium/tall grassland/woodland       ', & ! 10
                                   'desert, sandy                        ', & ! 11
                                   'desert, rocky                        ', & ! 12
                                   'tundra                               ', & ! 13
                                   'land ice                             ', & ! 14
                                   'sea ice                              ', & ! 15
                                   'snow                                 ', & ! 16
                                   'bare soil                            ', & ! 17
                                   'asphalt/concrete mix                 ', & ! 18
                                   'asphalt (asphalt concrete)           ', & ! 19
                                   'concrete (Portland concrete)         ', & ! 20
                                   'sett                                 ', & ! 21
                                   'paving stones                        ', & ! 22
                                   'cobblestone                          ', & ! 23
                                   'metal                                ', & ! 24
                                   'wood                                 ', & ! 25
                                   'gravel                               ', & ! 26
                                   'fine gravel                          ', & ! 27
                                   'pebblestone                          ', & ! 28
                                   'woodchips                            ', & ! 29
                                   'tartan (sports)                      ', & ! 30
                                   'artifical turf (sports)              ', & ! 31
                                   'clay (sports)                        ', & ! 32
                                   'building (dummy)                     '  & ! 33
                                                         /)
!
!-- Indices of radiation-related input attributes in building_surface_pars
!-- (other are in urban_surface_mod)
    INTEGER(iwp), PARAMETER ::  ind_s_alb_b_wall                = 19 !< index for Broadband albedo of wall fraction
    INTEGER(iwp), PARAMETER ::  ind_s_alb_l_wall                = 20 !< index for Longwave albedo of wall fraction
    INTEGER(iwp), PARAMETER ::  ind_s_alb_s_wall                = 21 !< index for Shortwave albedo of wall fraction
    INTEGER(iwp), PARAMETER ::  ind_s_alb_b_win                 = 22 !< index for Broadband albedo of window fraction
    INTEGER(iwp), PARAMETER ::  ind_s_alb_l_win                 = 23 !< index for Longwave albedo of window fraction
    INTEGER(iwp), PARAMETER ::  ind_s_alb_s_win                 = 24 !< index for Shortwave albedo of window fraction
    INTEGER(iwp), PARAMETER ::  ind_s_alb_b_green               = 24 !< index for Broadband albedo of green fraction
    INTEGER(iwp), PARAMETER ::  ind_s_alb_l_green               = 25 !< index for Longwave albedo of green fraction
    INTEGER(iwp), PARAMETER ::  ind_s_alb_s_green               = 26 !< index for Shortwave albedo of green fraction

    INTEGER(iwp) :: albedo_type  = 9999999_iwp, &     !< Albedo surface type
                    dots_rad     = 0_iwp              !< starting index for timeseries output

    LOGICAL ::  unscheduled_radiation_calls = .TRUE., & !< flag parameter indicating whether additional calls of the radiation code are allowed
                constant_albedo = .FALSE.,            & !< flag parameter indicating whether the albedo may change depending on zenith
                force_radiation_call = .FALSE.,       & !< flag parameter for unscheduled radiation calls
                lw_radiation = .TRUE.,                & !< flag parameter indicating whether longwave radiation shall be calculated
                radiation = .FALSE.,                  & !< flag parameter indicating whether the radiation model is used
                sun_up    = .TRUE.,                   & !< flag parameter indicating whether the sun is up or down
                sw_radiation = .TRUE.,                & !< flag parameter indicating whether shortwave radiation shall be calculated
                sun_direction = .FALSE.,              & !< flag parameter indicating whether solar direction shall be calculated
                average_radiation = .FALSE.,          & !< flag to set the calculation of radiation averaging for the domain
                radiation_interactions = .FALSE.,     & !< flag to activiate RTM (TRUE only if vertical urban/land surface and trees exist)
                surface_reflections = .TRUE.,         & !< flag to switch the calculation of radiation interaction between surfaces.
                                                        !< When it switched off, only the effect of buildings and trees shadow
                                                        !< will be considered. However fewer SVFs are expected.
                radiation_interactions_on = .TRUE.      !< namelist flag to force RTM activiation regardless to vertical urban/land surface and trees

    REAL(wp) :: albedo = 9999999.9_wp,           & !< NAMELIST alpha
                albedo_lw_dif = 9999999.9_wp,    & !< NAMELIST aldif
                albedo_lw_dir = 9999999.9_wp,    & !< NAMELIST aldir
                albedo_sw_dif = 9999999.9_wp,    & !< NAMELIST asdif
                albedo_sw_dir = 9999999.9_wp,    & !< NAMELIST asdir
                decl_1,                          & !< declination coef. 1
                decl_2,                          & !< declination coef. 2
                decl_3,                          & !< declination coef. 3
                dt_radiation = 0.0_wp,           & !< radiation model timestep
                emissivity = 9999999.9_wp,       & !< NAMELIST surface emissivity
                lon = 0.0_wp,                    & !< longitude in radians
                lat = 0.0_wp,                    & !< latitude in radians
                net_radiation = 0.0_wp,          & !< net radiation at surface
                skip_time_do_radiation = 0.0_wp, & !< Radiation model is not called before this time
                sky_trans,                       & !< sky transmissivity
                time_radiation = 0.0_wp,         & !< time since last call of radiation code
                trace_fluxes_above = -1.0_wp       !< NAMELIST option for debug tracing of large radiative fluxes (W/m2;W/m3)

    INTEGER(iwp) ::  day_of_year   !< day of the current year

    REAL(wp) ::  cos_zenith        !< cosine of solar zenith angle, also z-coordinate of solar unit vector
    REAL(wp) ::  d_hours_day       !< 1 / hours-per-day
    REAL(wp) ::  d_seconds_hour    !< 1 / seconds-per-hour
    REAL(wp) ::  second_of_day     !< second of the current day
    REAL(wp) ::  sun_dir_lat       !< y-coordinate of solar unit vector
    REAL(wp) ::  sun_dir_lon       !< x-coordinate of solar unit vector

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rad_net_av       !< average of net radiation (rad_net) at surface
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rad_lw_in_xy_av  !< average of incoming longwave radiation at surface
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rad_lw_out_xy_av !< average of outgoing longwave radiation at surface
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rad_sw_in_xy_av  !< average of incoming shortwave radiation at surface
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rad_sw_out_xy_av !< average of outgoing shortwave radiation at surface

    REAL(wp), PARAMETER :: emissivity_atm_clsky = 0.8_wp       !< emissivity of the clear-sky atmosphere
!
!-- Land surface albedos for solar zenith angle of 60degree after Briegleb (1992)
!-- (broadband, longwave, shortwave ):   bb,      lw,      sw,
    REAL(wp), DIMENSION(0:2,1:33), PARAMETER :: albedo_pars = RESHAPE( (/&
                                   0.06_wp, 0.06_wp, 0.06_wp,            & !  1
                                   0.19_wp, 0.28_wp, 0.09_wp,            & !  2
                                   0.23_wp, 0.33_wp, 0.11_wp,            & !  3
                                   0.23_wp, 0.33_wp, 0.11_wp,            & !  4
                                   0.25_wp, 0.34_wp, 0.14_wp,            & !  5
                                   0.14_wp, 0.22_wp, 0.06_wp,            & !  6
                                   0.17_wp, 0.27_wp, 0.06_wp,            & !  7
                                   0.19_wp, 0.31_wp, 0.06_wp,            & !  8
                                   0.14_wp, 0.22_wp, 0.06_wp,            & !  9
                                   0.18_wp, 0.28_wp, 0.06_wp,            & ! 10
                                   0.43_wp, 0.51_wp, 0.35_wp,            & ! 11
                                   0.32_wp, 0.40_wp, 0.24_wp,            & ! 12
                                   0.19_wp, 0.27_wp, 0.10_wp,            & ! 13
                                   0.77_wp, 0.65_wp, 0.90_wp,            & ! 14
                                   0.77_wp, 0.65_wp, 0.90_wp,            & ! 15
                                   0.82_wp, 0.70_wp, 0.95_wp,            & ! 16
                                   0.08_wp, 0.08_wp, 0.08_wp,            & ! 17
                                   0.17_wp, 0.17_wp, 0.17_wp,            & ! 18
                                   0.17_wp, 0.17_wp, 0.17_wp,            & ! 19
                                   0.30_wp, 0.30_wp, 0.30_wp,            & ! 20
                                   0.17_wp, 0.17_wp, 0.17_wp,            & ! 21
                                   0.17_wp, 0.17_wp, 0.17_wp,            & ! 22
                                   0.17_wp, 0.17_wp, 0.17_wp,            & ! 23
                                   0.17_wp, 0.17_wp, 0.17_wp,            & ! 24
                                   0.17_wp, 0.17_wp, 0.17_wp,            & ! 25
                                   0.17_wp, 0.17_wp, 0.17_wp,            & ! 26
                                   0.17_wp, 0.17_wp, 0.17_wp,            & ! 27
                                   0.17_wp, 0.17_wp, 0.17_wp,            & ! 28
                                   0.17_wp, 0.17_wp, 0.17_wp,            & ! 29
                                   0.17_wp, 0.17_wp, 0.17_wp,            & ! 30
                                   0.17_wp, 0.17_wp, 0.17_wp,            & ! 31
                                   0.17_wp, 0.17_wp, 0.17_wp,            & ! 32
                                   0.17_wp, 0.17_wp, 0.17_wp             & ! 33
                                 /), (/ 3, 33 /) )

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: &
                        rad_lw_cs_hr,                  & !< longwave clear sky radiation heating rate (K/s)
                        rad_lw_cs_hr_av,               & !< average of rad_lw_cs_hr
                        rad_lw_hr,                     & !< longwave radiation heating rate (K/s)
                        rad_lw_hr_av,                  & !< average of rad_sw_hr
                        rad_lw_in,                     & !< incoming longwave radiation (W/m2)
                        rad_lw_in_av,                  & !< average of rad_lw_in
                        rad_lw_out,                    & !< outgoing longwave radiation (W/m2)
                        rad_lw_out_av,                 & !< average of rad_lw_out
                        rad_sw_cs_hr,                  & !< shortwave clear sky radiation heating rate (K/s)
                        rad_sw_cs_hr_av,               & !< average of rad_sw_cs_hr
                        rad_sw_hr,                     & !< shortwave radiation heating rate (K/s)
                        rad_sw_hr_av,                  & !< average of rad_sw_hr
                        rad_sw_in,                     & !< incoming shortwave radiation (W/m2)
                        rad_sw_in_av,                  & !< average of rad_sw_in
                        rad_sw_out,                    & !< outgoing shortwave radiation (W/m2)
                        rad_sw_out_av                    !< average of rad_sw_out


!
!-- Variables and parameters used in RRTMG only
#if defined ( __rrtmg )
    CHARACTER(LEN=12) :: rrtm_input_file = "RAD_SND_DATA" !< name of the NetCDF input file (sounding data)


!
!-- Flag parameters to be passed to RRTMG (should not be changed until ice phase in clouds is allowed)
    INTEGER(iwp), PARAMETER :: rrtm_idrv     = 1, & !< flag for longwave upward flux calculation option (0,1)
                               rrtm_inflglw  = 2, & !< flag for lw cloud optical properties (0,1,2)
                               rrtm_iceflglw = 0, & !< flag for lw ice particle specifications (0,1,2,3)
                               rrtm_liqflglw = 1, & !< flag for lw liquid droplet specifications
                               rrtm_inflgsw  = 2, & !< flag for sw cloud optical properties (0,1,2)
                               rrtm_iceflgsw = 0, & !< flag for sw ice particle specifications (0,1,2,3)
                               rrtm_liqflgsw = 1    !< flag for sw liquid droplet specifications

!
!-- The following variables should be only changed with care, as this will
!-- require further setting of some variables, which is currently not
!-- implemented (aerosols, ice phase).
    INTEGER(iwp) :: nzt_rad,           & !< upper vertical limit for radiation calculations
                    rrtm_icld = 0,     & !< cloud flag (0: clear sky column, 1: cloudy column)
                    rrtm_iaer = 0        !< aerosol option flag (0: no aerosol layers, for lw only: 6 (requires setting of rrtm_sw_ecaer), 10: one or more aerosol layers (not implemented)

    INTEGER(iwp) :: nc_stat !< local variable for storin the result of netCDF calls for error message handling

    LOGICAL :: snd_exists = .FALSE.      !< flag parameter to check whether a user-defined input files exists
    LOGICAL :: sw_exists = .FALSE.       !< flag parameter to check whether that required rrtmg sw file exists
    LOGICAL :: lw_exists = .FALSE.       !< flag parameter to check whether that required rrtmg lw file exists


    REAL(wp), PARAMETER :: mol_mass_air_d_wv = 1.607793_wp !< molecular weight dry air / water vapor

    REAL(wp), DIMENSION(:), ALLOCATABLE :: hyp_snd,     & !< hypostatic pressure from sounding data (hPa)
                                           rrtm_tsfc,   & !< dummy array for storing surface temperature
                                           t_snd          !< actual temperature from sounding data (hPa)

    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: rrtm_ccl4vmr,   & !< CCL4 volume mixing ratio (g/mol)
                                             rrtm_cfc11vmr,  & !< CFC11 volume mixing ratio (g/mol)
                                             rrtm_cfc12vmr,  & !< CFC12 volume mixing ratio (g/mol)
                                             rrtm_cfc22vmr,  & !< CFC22 volume mixing ratio (g/mol)
                                             rrtm_ch4vmr,    & !< CH4 volume mixing ratio
                                             rrtm_cicewp,    & !< in-cloud ice water path (g/m2)
                                             rrtm_cldfr,     & !< cloud fraction (0,1)
                                             rrtm_cliqwp,    & !< in-cloud liquid water path (g/m2)
                                             rrtm_co2vmr,    & !< CO2 volume mixing ratio (g/mol)
                                             rrtm_emis,      & !< surface emissivity (0-1)
                                             rrtm_h2ovmr,    & !< H2O volume mixing ratio
                                             rrtm_n2ovmr,    & !< N2O volume mixing ratio
                                             rrtm_o2vmr,     & !< O2 volume mixing ratio
                                             rrtm_o3vmr,     & !< O3 volume mixing ratio
                                             rrtm_play,      & !< pressure layers (hPa, zu-grid)
                                             rrtm_plev,      & !< pressure layers (hPa, zw-grid)
                                             rrtm_reice,     & !< cloud ice effective radius (microns)
                                             rrtm_reliq,     & !< cloud water drop effective radius (microns)
                                             rrtm_tlay,      & !< actual temperature (K, zu-grid)
                                             rrtm_tlev,      & !< actual temperature (K, zw-grid)
                                             rrtm_lwdflx,    & !< RRTM output of incoming longwave radiation flux (W/m2)
                                             rrtm_lwdflxc,   & !< RRTM output of outgoing clear sky longwave radiation flux (W/m2)
                                             rrtm_lwuflx,    & !< RRTM output of outgoing longwave radiation flux (W/m2)
                                             rrtm_lwuflxc,   & !< RRTM output of incoming clear sky longwave radiation flux (W/m2)
                                             rrtm_lwuflx_dt, & !< RRTM output of incoming clear sky longwave radiation flux (W/m2)
                                             rrtm_lwuflxc_dt,& !< RRTM output of outgoing clear sky longwave radiation flux (W/m2)
                                             rrtm_lwhr,      & !< RRTM output of longwave radiation heating rate (K/d)
                                             rrtm_lwhrc,     & !< RRTM output of incoming longwave clear sky radiation heating rate (K/d)
                                             rrtm_swdflx,    & !< RRTM output of incoming shortwave radiation flux (W/m2)
                                             rrtm_swdflxc,   & !< RRTM output of outgoing clear sky shortwave radiation flux (W/m2)
                                             rrtm_swuflx,    & !< RRTM output of outgoing shortwave radiation flux (W/m2)
                                             rrtm_swuflxc,   & !< RRTM output of incoming clear sky shortwave radiation flux (W/m2)
                                             rrtm_swhr,      & !< RRTM output of shortwave radiation heating rate (K/d)
                                             rrtm_swhrc,     & !< RRTM output of incoming shortwave clear sky radiation heating rate (K/d)
                                             rrtm_dirdflux,  & !< RRTM output of incoming direct shortwave (W/m2)
                                             rrtm_difdflux     !< RRTM output of incoming diffuse shortwave (W/m2)

    REAL(wp), DIMENSION(1) ::                rrtm_aldif,     & !< surface albedo for longwave diffuse radiation
                                             rrtm_aldir,     & !< surface albedo for longwave direct radiation
                                             rrtm_asdif,     & !< surface albedo for shortwave diffuse radiation
                                             rrtm_asdir        !< surface albedo for shortwave direct radiation

!
!-- Definition of arrays that are currently not used for calling RRTMG (due to setting of flag parameters)
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  rad_lw_cs_in,   & !< incoming clear sky longwave radiation (W/m2) (not used)
                                                rad_lw_cs_out,  & !< outgoing clear sky longwave radiation (W/m2) (not used)
                                                rad_sw_cs_in,   & !< incoming clear sky shortwave radiation (W/m2) (not used)
                                                rad_sw_cs_out,  & !< outgoing clear sky shortwave radiation (W/m2) (not used)
                                                rrtm_lw_tauaer, & !< lw aerosol optical depth
                                                rrtm_lw_taucld, & !< lw in-cloud optical depth
                                                rrtm_sw_taucld, & !< sw in-cloud optical depth
                                                rrtm_sw_ssacld, & !< sw in-cloud single scattering albedo
                                                rrtm_sw_asmcld, & !< sw in-cloud asymmetry parameter
                                                rrtm_sw_fsfcld, & !< sw in-cloud forward scattering fraction
                                                rrtm_sw_tauaer, & !< sw aerosol optical depth
                                                rrtm_sw_ssaaer, & !< sw aerosol single scattering albedo
                                                rrtm_sw_asmaer, & !< sw aerosol asymmetry parameter
                                                rrtm_sw_ecaer     !< sw aerosol optical detph at 0.55 microns (rrtm_iaer = 6 only)

#endif
!
!-- Parameters of urban and land surface models
    INTEGER(iwp)                                   ::  nz_urban                           !< number of layers of urban surface (will be calculated)
    INTEGER(iwp)                                   ::  nz_plant                           !< number of layers of plant canopy (will be calculated)
    INTEGER(iwp)                                   ::  nz_urban_b                         !< bottom layer of urban surface (will be calculated)
    INTEGER(iwp)                                   ::  nz_urban_t                         !< top layer of urban surface (will be calculated)
    INTEGER(iwp)                                   ::  nz_plant_t                         !< top layer of plant canopy (will be calculated)
!-- parameters of urban and land surface models
    INTEGER(iwp), PARAMETER                        ::  nzut_free = 3                      !< number of free layers above top of of topography
    INTEGER(iwp), PARAMETER                        ::  ndsvf = 2                          !< number of dimensions of real values in SVF
    INTEGER(iwp), PARAMETER                        ::  idsvf = 2                          !< number of dimensions of integer values in SVF
    INTEGER(iwp), PARAMETER                        ::  ndcsf = 1                          !< number of dimensions of real values in CSF
    INTEGER(iwp), PARAMETER                        ::  idcsf = 2                          !< number of dimensions of integer values in CSF
    INTEGER(iwp), PARAMETER                        ::  kdcsf = 4                          !< number of dimensions of integer values in CSF calculation array
    INTEGER(iwp), PARAMETER                        ::  id = 1                             !< position of d-index in surfl and surf
    INTEGER(iwp), PARAMETER                        ::  iz = 2                             !< position of k-index in surfl and surf
    INTEGER(iwp), PARAMETER                        ::  iy = 3                             !< position of j-index in surfl and surf
    INTEGER(iwp), PARAMETER                        ::  ix = 4                             !< position of i-index in surfl and surf
    INTEGER(iwp), PARAMETER                        ::  im = 5                             !< position of surface m-index in surfl and surf
    INTEGER(iwp), PARAMETER                        ::  nidx_surf = 5                      !< number of indices in surfl and surf

    INTEGER(iwp), PARAMETER                        ::  nsurf_type = 10                    !< number of surf types incl. phys.(land+urban) & (atm.,sky,boundary) surfaces - 1

    INTEGER(iwp), PARAMETER                        ::  iup_u    = 0                       !< 0 - index of urban upward surface (ground or roof)
    INTEGER(iwp), PARAMETER                        ::  idown_u  = 1                       !< 1 - index of urban downward surface (overhanging)
    INTEGER(iwp), PARAMETER                        ::  inorth_u = 2                       !< 2 - index of urban northward facing wall
    INTEGER(iwp), PARAMETER                        ::  isouth_u = 3                       !< 3 - index of urban southward facing wall
    INTEGER(iwp), PARAMETER                        ::  ieast_u  = 4                       !< 4 - index of urban eastward facing wall
    INTEGER(iwp), PARAMETER                        ::  iwest_u  = 5                       !< 5 - index of urban westward facing wall

    INTEGER(iwp), PARAMETER                        ::  iup_l    = 6                       !< 6 - index of land upward surface (ground or roof)
    INTEGER(iwp), PARAMETER                        ::  inorth_l = 7                       !< 7 - index of land northward facing wall
    INTEGER(iwp), PARAMETER                        ::  isouth_l = 8                       !< 8 - index of land southward facing wall
    INTEGER(iwp), PARAMETER                        ::  ieast_l  = 9                       !< 9 - index of land eastward facing wall
    INTEGER(iwp), PARAMETER                        ::  iwest_l  = 10                      !< 10- index of land westward facing wall

    INTEGER(iwp), DIMENSION(0:nsurf_type), PARAMETER ::  idir = (/0, 0,0, 0,1,-1,0,0, 0,1,-1/)   !< surface normal direction x indices
    INTEGER(iwp), DIMENSION(0:nsurf_type), PARAMETER ::  jdir = (/0, 0,1,-1,0, 0,0,1,-1,0, 0/)   !< surface normal direction y indices
    INTEGER(iwp), DIMENSION(0:nsurf_type), PARAMETER ::  kdir = (/1,-1,0, 0,0, 0,1,0, 0,0, 0/)   !< surface normal direction z indices
    REAL(wp),     DIMENSION(0:nsurf_type)          :: facearea                            !< area of single face in respective
                                                                                          !< direction (will be calc'd)


!-- indices and sizes of urban and land surface models
    INTEGER(iwp)                                   ::  startland        !< start index of block of land and roof surfaces
    INTEGER(iwp)                                   ::  endland          !< end index of block of land and roof surfaces
    INTEGER(iwp)                                   ::  nlands           !< number of land and roof surfaces in local processor
    INTEGER(iwp)                                   ::  startwall        !< start index of block of wall surfaces
    INTEGER(iwp)                                   ::  endwall          !< end index of block of wall surfaces
    INTEGER(iwp)                                   ::  nwalls           !< number of wall surfaces in local processor

!-- indices needed for RTM netcdf output subroutines
    INTEGER(iwp), PARAMETER                        :: nd = 5
    CHARACTER(LEN=6), DIMENSION(0:nd-1), PARAMETER :: dirname = (/ '_roof ', '_south', '_north', '_west ', '_east ' /)
    INTEGER(iwp), DIMENSION(0:nd-1), PARAMETER     :: dirint_u = (/ iup_u, isouth_u, inorth_u, iwest_u, ieast_u /)
    INTEGER(iwp), DIMENSION(0:nd-1), PARAMETER     :: dirint_l = (/ iup_l, isouth_l, inorth_l, iwest_l, ieast_l /)
    INTEGER(iwp), DIMENSION(0:nd-1)                :: dirstart
    INTEGER(iwp), DIMENSION(0:nd-1)                :: dirend

!-- indices and sizes of urban and land surface models
    INTEGER(iwp), DIMENSION(:,:), POINTER          ::  surfl            !< coordinates of i-th local surface in local grid - surfl[:,k] = [d, z, y, x, m]
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE,TARGET ::  surfl_linear     !< dtto (linearly allocated array)
    INTEGER(iwp), DIMENSION(:,:), POINTER          ::  surf             !< coordinates of i-th surface in grid - surf[:,k] = [d, z, y, x, m]
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE,TARGET ::  surf_linear      !< dtto (linearly allocated array)
    INTEGER(iwp)                                   ::  nsurfl           !< number of all surfaces in local processor
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE,TARGET ::  nsurfs           !< array of number of all surfaces in individual processors
    INTEGER(iwp)                                   ::  nsurf            !< global number of surfaces in index array of surfaces (nsurf = proc nsurfs)
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE,TARGET ::  surfstart        !< starts of blocks of surfaces for individual processors in array surf (indexed from 1)
                                                                        !< respective block for particular processor is surfstart[iproc+1]+1 : surfstart[iproc+1]+nsurfs[iproc+1]

!-- block variables needed for calculation of the plant canopy model inside the urban surface model
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE      ::  pct              !< top layer of the plant canopy
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE      ::  pch              !< heights of the plant canopy
    INTEGER(iwp)                                   ::  npcbl = 0        !< number of the plant canopy gridboxes in local processor
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE      ::  pcbl             !< k,j,i coordinates of l-th local plant canopy box pcbl[:,l] = [k, j, i]
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  pcbinsw          !< array of absorbed sw radiation for local plant canopy box
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  pcbinswdir       !< array of absorbed direct sw radiation for local plant canopy box
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  pcbinswdif       !< array of absorbed diffusion sw radiation for local plant canopy box
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  pcbinlw          !< array of absorbed lw radiation for local plant canopy box

!-- configuration parameters (they can be setup in PALM config)
    LOGICAL                                        ::  raytrace_mpi_rma = .TRUE.          !< use MPI RMA to access LAD and gridsurf from remote processes during raytracing
    LOGICAL                                        ::  rad_angular_discretization = .TRUE.!< whether to use fixed resolution discretization of view factors for
                                                                                          !< reflected radiation (as opposed to all mutually visible pairs)
    LOGICAL                                        ::  plant_lw_interact = .TRUE.         !< whether plant canopy interacts with LW radiation (in addition to SW)
    INTEGER(iwp)                                   ::  mrt_nlevels = 0                    !< number of vertical boxes above surface for which to calculate MRT
    LOGICAL                                        ::  mrt_skip_roof = .TRUE.             !< do not calculate MRT above roof surfaces
    LOGICAL                                        ::  mrt_include_sw = .TRUE.            !< should MRT calculation include SW radiation as well?
    INTEGER(wp)                                    ::  mrt_geom = 1                       !< method for MRT direction weights simulating a sphere or a human body
    REAL(wp), DIMENSION(2)                         ::  mrt_geom_params = (/ .12_wp, .88_wp /)   !< parameters for the selected method
    INTEGER(iwp)                                   ::  nrefsteps = 3                      !< number of reflection steps to perform
    REAL(wp), PARAMETER                            ::  ext_coef = 0.6_wp                  !< extinction coefficient (a.k.a. alpha)
    INTEGER(iwp), PARAMETER                        ::  rad_version_len = 10               !< length of identification string of rad version
    CHARACTER(rad_version_len), PARAMETER          ::  rad_version = 'RAD v. 3.0'         !< identification of version of binary svf and restart files
    INTEGER(iwp)                                   ::  raytrace_discrete_elevs = 40       !< number of discretization steps for elevation (nadir to zenith)
    INTEGER(iwp)                                   ::  raytrace_discrete_azims = 80       !< number of discretization steps for azimuth (out of 360 degrees)
    REAL(wp)                                       ::  max_raytracing_dist = -999.0_wp    !< maximum distance for raytracing (in metres)
    REAL(wp)                                       ::  min_irrf_value = 1e-6_wp           !< minimum potential irradiance factor value for raytracing
    REAL(wp), DIMENSION(1:30)                      ::  svfnorm_report_thresh = 1e21_wp    !< thresholds of SVF normalization values to report
    INTEGER(iwp)                                   ::  svfnorm_report_num                 !< number of SVF normalization thresholds to report

!-- radiation related arrays to be used in radiation_interaction routine
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          ::  rad_sw_in_dir    !< direct sw radiation
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          ::  rad_sw_in_diff   !< diffusion sw radiation
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          ::  rad_lw_in_diff   !< diffusion lw radiation

!-- parameters required for RRTMG lower boundary condition
    REAL(wp)                   :: albedo_urb      !< albedo value retuned to RRTMG boundary cond.
    REAL(wp)                   :: emissivity_urb  !< emissivity value retuned to RRTMG boundary cond.
    REAL(wp)                   :: t_rad_urb       !< temperature value retuned to RRTMG boundary cond.

!-- type for calculation of svf
    TYPE t_svf
        INTEGER(iwp)                               :: isurflt           !<
        INTEGER(iwp)                               :: isurfs            !<
        REAL(wp)                                   :: rsvf              !<
        REAL(wp)                                   :: rtransp           !<
    END TYPE

!-- type for calculation of csf
    TYPE t_csf
        INTEGER(iwp)                               :: ip                !<
        INTEGER(iwp)                               :: itx               !<
        INTEGER(iwp)                               :: ity               !<
        INTEGER(iwp)                               :: itz               !<
        INTEGER(iwp)                               :: isurfs            !< Idx of source face / -1 for sky
        REAL(wp)                                   :: rcvf              !< Canopy view factor for faces /
                                                                        !< canopy sink factor for sky (-1)
    END TYPE

!-- arrays storing the values of USM
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE      ::  svfsurf          !< svfsurf[:,isvf] = index of target and source surface for svf[isvf]
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          ::  svf              !< array of shape view factors+direct irradiation factors for local surfaces
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfins          !< array of sw radiation falling to local surface after i-th reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinl          !< array of lw radiation for local surface after i-th reflection

    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  skyvf            !< array of sky view factor for each local surface
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  skyvft           !< array of sky view factor including transparency for each local surface
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          ::  dsitrans         !< dsidir[isvfl,i] = path transmittance of i-th
                                                                        !< direction of direct solar irradiance per target surface
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          ::  dsitransc        !< dtto per plant canopy box
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          ::  dsidir           !< dsidir[:,i] = unit vector of i-th
                                                                        !< direction of direct solar irradiance
    INTEGER(iwp)                                   ::  ndsidir          !< number of apparent solar directions used
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE      ::  dsidir_rev       !< dsidir_rev[ielev,iazim] = i for dsidir or -1 if not present

    INTEGER(iwp)                                   ::  nmrtbl           !< No. of local grid boxes for which MRT is calculated
    INTEGER(iwp)                                   ::  nmrtf            !< number of MRT factors for local processor
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE      ::  mrtbl            !< coordinates of i-th local MRT box - surfl[:,i] = [z, y, x]
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE      ::  mrtfsurf         !< mrtfsurf[:,imrtf] = index of target MRT box and source surface for mrtf[imrtf]
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  mrtf             !< array of MRT factors for each local MRT box
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  mrtft            !< array of MRT factors including transparency for each local MRT box
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  mrtsky           !< array of sky view factor for each local MRT box
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  mrtskyt          !< array of sky view factor including transparency for each local MRT box
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          ::  mrtdsit          !< array of direct solar transparencies for each local MRT box
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  mrtinsw          !< mean SW radiant flux for each MRT box
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  mrtinlw          !< mean LW radiant flux for each MRT box
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  mrt              !< mean radiant temperature for each MRT box
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  mrtinsw_av       !< time average mean SW radiant flux for each MRT box
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  mrtinlw_av       !< time average mean LW radiant flux for each MRT box
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  mrt_av           !< time average mean radiant temperature for each MRT box

    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinsw         !< array of sw radiation falling to local surface including radiation from reflections
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinlw         !< array of lw radiation falling to local surface including radiation from reflections
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinswdir      !< array of direct sw radiation falling to local surface
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinswdif      !< array of diffuse sw radiation from sky and model boundary falling to local surface
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinlwdif      !< array of diffuse lw radiation from sky and model boundary falling to local surface

                                                                        !< Outward radiation is only valid for nonvirtual surfaces
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfoutsl        !< array of reflected sw radiation for local surface in i-th reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfoutll        !< array of reflected + emitted lw radiation for local surface in i-th reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfouts         !< array of reflected sw radiation for all surfaces in i-th reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfoutl         !< array of reflected + emitted lw radiation for all surfaces in i-th reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinlg         !< global array of incoming lw radiation from plant canopy
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfoutsw        !< array of total sw radiation outgoing from nonvirtual surfaces surfaces after all reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfoutlw        !< array of total lw radiation outgoing from nonvirtual surfaces surfaces after all reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfemitlwl      !< array of emitted lw radiation for local surface used to calculate effective surface temperature for radiation model

!-- block variables needed for calculation of the plant canopy model inside the urban surface model
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE      ::  csfsurf          !< csfsurf[:,icsf] = index of target surface and csf grid index for csf[icsf]
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          ::  csf              !< array of plant canopy sink fators + direct irradiation factors (transparency)
    REAL(wp), DIMENSION(:,:,:), POINTER            ::  sub_lad          !< subset of lad_s within urban surface, transformed to plain Z coordinate
#if defined( __parallel )
    REAL(wp), DIMENSION(:), POINTER                ::  sub_lad_g        !< sub_lad globalized (used to avoid MPI RMA calls in raytracing)
#endif
    REAL(wp)                                       ::  prototype_lad    !< prototype leaf area density for computing effective optical depth
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE        ::  nzterr, plantt   !< temporary global arrays for raytracing
    INTEGER(iwp)                                   ::  plantt_max

!-- arrays and variables for calculation of svf and csf
    TYPE(t_svf), DIMENSION(:), POINTER             ::  asvf             !< pointer to growing svc array
    TYPE(t_csf), DIMENSION(:), POINTER             ::  acsf             !< pointer to growing csf array
    TYPE(t_svf), DIMENSION(:), POINTER             ::  amrtf            !< pointer to growing mrtf array
    TYPE(t_svf), DIMENSION(:), ALLOCATABLE, TARGET ::  asvf1, asvf2     !< realizations of svf array
    TYPE(t_csf), DIMENSION(:), ALLOCATABLE, TARGET ::  acsf1, acsf2     !< realizations of csf array
    TYPE(t_svf), DIMENSION(:), ALLOCATABLE, TARGET ::  amrtf1, amrtf2   !< realizations of mftf array
    INTEGER(iwp)                                   ::  nsvfla           !< dimmension of array allocated for storage of svf in local processor
    INTEGER(iwp)                                   ::  ncsfla           !< dimmension of array allocated for storage of csf in local processor
    INTEGER(iwp)                                   ::  nmrtfa           !< dimmension of array allocated for storage of mrt
    INTEGER(iwp)                                   ::  msvf, mcsf, mmrtf!< mod for swapping the growing array
    INTEGER(iwp), PARAMETER                        ::  gasize = 100000_iwp  !< initial size of growing arrays
    REAL(wp), PARAMETER                            ::  grow_factor = 1.4_wp !< growth factor of growing arrays
    INTEGER(iwp)                                   ::  nsvfl            !< number of svf for local processor
    INTEGER(iwp)                                   ::  ncsfl            !< no. of csf in local processor
                                                                        !< needed only during calc_svf but must be here because it is
                                                                        !< shared between subroutines calc_svf and raytrace
    INTEGER(iwp), DIMENSION(:,:,:,:), POINTER      ::  gridsurf         !< reverse index of local surfl[d,k,j,i] (for case rad_angular_discretization)
    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE    ::  gridpcbl         !< reverse index of local pcbl[k,j,i]
    INTEGER(iwp), PARAMETER                        ::  nsurf_type_u = 6 !< number of urban surf types (used in gridsurf)

!-- temporary arrays for calculation of csf in raytracing
    INTEGER(iwp)                                   ::  maxboxesg        !< max number of boxes ray can cross in the domain
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE      ::  boxes            !< coordinates of gridboxes being crossed by ray
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  crlens           !< array of crossing lengths of ray for particular grid boxes
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE        ::  lad_ip           !< array of numbers of process where lad is stored
#if defined( __parallel )
    INTEGER(kind=MPI_ADDRESS_KIND), &
                  DIMENSION(:), ALLOCATABLE        ::  lad_disp         !< array of displaycements of lad in local array of proc lad_ip
    INTEGER(iwp)                                   ::  win_lad          !< MPI RMA window for leaf area density
    INTEGER(iwp)                                   ::  win_gridsurf     !< MPI RMA window for reverse grid surface index
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  lad_s_ray        !< array of received lad_s for appropriate gridboxes crossed by ray
#endif
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE        ::  target_surfl
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE      ::  rt2_track
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          ::  rt2_track_lad
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  rt2_track_dist
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  rt2_dist

!-- arrays for time averages
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfradnet_av    !< average of net radiation to local surface including radiation from reflections
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinsw_av      !< average of sw radiation falling to local surface including radiation from reflections
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinlw_av      !< average of lw radiation falling to local surface including radiation from reflections
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinswdir_av   !< average of direct sw radiation falling to local surface
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinswdif_av   !< average of diffuse sw radiation from sky and model boundary falling to local surface
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinlwdif_av   !< average of diffuse lw radiation from sky and model boundary falling to local surface
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinswref_av   !< average of sw radiation falling to surface from reflections
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinlwref_av   !< average of lw radiation falling to surface from reflections
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfoutsw_av     !< average of total sw radiation outgoing from nonvirtual surfaces surfaces after all reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfoutlw_av     !< average of total lw radiation outgoing from nonvirtual surfaces surfaces after all reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfins_av       !< average of array of residua of sw radiation absorbed in surface after last reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinl_av       !< average of array of residua of lw radiation absorbed in surface after last reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  pcbinlw_av       !< Average of pcbinlw
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  pcbinsw_av       !< Average of pcbinsw
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  pcbinswdir_av    !< Average of pcbinswdir
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  pcbinswdif_av    !< Average of pcbinswdif
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  pcbinswref_av    !< Average of pcbinswref


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-- Energy balance variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-- parameters of the land, roof and wall surfaces
    REAL(wp), DIMENSION(:), ALLOCATABLE            :: albedo_surf        !< albedo of the surface
    REAL(wp), DIMENSION(:), ALLOCATABLE            :: emiss_surf         !< emissivity of the wall surface
!
!-- External radiation. Depending on the given level of detail either a 1D or
!-- a 3D array will be allocated.
    TYPE( real_1d_3d ) ::  rad_lw_in_f     !< external incoming longwave radiation, from observation or model
    TYPE( real_1d_3d ) ::  rad_sw_in_f     !< external incoming shortwave radiation, from observation or model
    TYPE( real_1d_3d ) ::  rad_sw_in_dif_f !< external incoming shortwave radiation, diffuse part, from observation or model
    TYPE( real_1d_3d ) ::  time_rad_f      !< time dimension for external radiation, from observation or model

    INTERFACE radiation_check_data_output
       MODULE PROCEDURE radiation_check_data_output
    END INTERFACE radiation_check_data_output

    INTERFACE radiation_check_data_output_ts
       MODULE PROCEDURE radiation_check_data_output_ts
    END INTERFACE radiation_check_data_output_ts

    INTERFACE radiation_check_data_output_pr
       MODULE PROCEDURE radiation_check_data_output_pr
    END INTERFACE radiation_check_data_output_pr

    INTERFACE radiation_check_parameters
       MODULE PROCEDURE radiation_check_parameters
    END INTERFACE radiation_check_parameters

    INTERFACE radiation_clearsky
       MODULE PROCEDURE radiation_clearsky
    END INTERFACE radiation_clearsky

    INTERFACE radiation_constant
       MODULE PROCEDURE radiation_constant
    END INTERFACE radiation_constant

    INTERFACE radiation_control
       MODULE PROCEDURE radiation_control
    END INTERFACE radiation_control

    INTERFACE radiation_3d_data_averaging
       MODULE PROCEDURE radiation_3d_data_averaging
    END INTERFACE radiation_3d_data_averaging

    INTERFACE radiation_data_output_2d
       MODULE PROCEDURE radiation_data_output_2d
    END INTERFACE radiation_data_output_2d

    INTERFACE radiation_data_output_3d
       MODULE PROCEDURE radiation_data_output_3d
    END INTERFACE radiation_data_output_3d

    INTERFACE radiation_data_output_mask
       MODULE PROCEDURE radiation_data_output_mask
    END INTERFACE radiation_data_output_mask

    INTERFACE radiation_define_netcdf_grid
       MODULE PROCEDURE radiation_define_netcdf_grid
    END INTERFACE radiation_define_netcdf_grid

    INTERFACE radiation_header
       MODULE PROCEDURE radiation_header
    END INTERFACE radiation_header

    INTERFACE radiation_init
       MODULE PROCEDURE radiation_init
    END INTERFACE radiation_init

    INTERFACE radiation_parin
       MODULE PROCEDURE radiation_parin
    END INTERFACE radiation_parin

    INTERFACE radiation_rrtmg
       MODULE PROCEDURE radiation_rrtmg
    END INTERFACE radiation_rrtmg

#if defined( __rrtmg )
    INTERFACE radiation_tendency
       MODULE PROCEDURE radiation_tendency
       MODULE PROCEDURE radiation_tendency_ij
    END INTERFACE radiation_tendency
#endif

    INTERFACE radiation_rrd_local
       MODULE PROCEDURE radiation_rrd_local
    END INTERFACE radiation_rrd_local

    INTERFACE radiation_wrd_local
       MODULE PROCEDURE radiation_wrd_local
    END INTERFACE radiation_wrd_local

    INTERFACE radiation_interaction
       MODULE PROCEDURE radiation_interaction
    END INTERFACE radiation_interaction

    INTERFACE radiation_interaction_init
       MODULE PROCEDURE radiation_interaction_init
    END INTERFACE radiation_interaction_init

    INTERFACE radiation_presimulate_solar_pos
       MODULE PROCEDURE radiation_presimulate_solar_pos
    END INTERFACE radiation_presimulate_solar_pos

    INTERFACE radiation_calc_svf
       MODULE PROCEDURE radiation_calc_svf
    END INTERFACE radiation_calc_svf

    INTERFACE radiation_write_svf
       MODULE PROCEDURE radiation_write_svf
    END INTERFACE radiation_write_svf

    INTERFACE radiation_read_svf
       MODULE PROCEDURE radiation_read_svf
    END INTERFACE radiation_read_svf


    SAVE

    PRIVATE

!
!-- Public functions / NEEDS SORTING
    PUBLIC radiation_check_data_output, radiation_check_data_output_pr,        &
           radiation_check_data_output_ts,                                     &
           radiation_check_parameters, radiation_control,                      &
           radiation_header, radiation_init, radiation_parin,                  &
           radiation_3d_data_averaging,                                        &
           radiation_data_output_2d, radiation_data_output_3d,                 &
           radiation_define_netcdf_grid, radiation_wrd_local,                  &
           radiation_rrd_local, radiation_data_output_mask,                    &
           radiation_calc_svf, radiation_write_svf,                            &
           radiation_interaction, radiation_interaction_init,                  &
           radiation_read_svf, radiation_presimulate_solar_pos


!
!-- Public variables and constants / NEEDS SORTING
    PUBLIC albedo, albedo_type, decl_1, decl_2, decl_3, dots_rad, dt_radiation,&
           emissivity, force_radiation_call, lat, lon, mrt_geom,               &
           mrt_geom_params,                                                    &
           mrt_include_sw, mrt_nlevels, mrtbl, mrtinsw, mrtinlw, nmrtbl,       &
           rad_net_av, radiation, radiation_scheme, rad_lw_in,                 &
           rad_lw_in_av, rad_lw_out, rad_lw_out_av,                            &
           rad_lw_cs_hr, rad_lw_cs_hr_av, rad_lw_hr, rad_lw_hr_av, rad_sw_in,  &
           rad_sw_in_av, rad_sw_out, rad_sw_out_av, rad_sw_cs_hr,              &
           rad_sw_cs_hr_av, rad_sw_hr, rad_sw_hr_av, solar_constant,           &
           skip_time_do_radiation, time_radiation, unscheduled_radiation_calls,&
           cos_zenith, calc_zenith, sun_direction, sun_dir_lat, sun_dir_lon,   &
           idir, jdir, kdir, id, iz, iy, ix,                                   &
           iup_u, inorth_u, isouth_u, ieast_u, iwest_u,                        &
           iup_l, inorth_l, isouth_l, ieast_l, iwest_l,                        &
           nsurf_type, nz_urban_b, nz_urban_t, nz_urban, pch, nsurf,                 &
           idsvf, ndsvf, idcsf, ndcsf, kdcsf, pct,                             &
           radiation_interactions, startwall, startland, endland, endwall,     &
           skyvf, skyvft, radiation_interactions_on, average_radiation,        &
           rad_sw_in_diff, rad_sw_in_dir


#if defined ( __rrtmg )
    PUBLIC radiation_tendency, rrtm_aldif, rrtm_aldir, rrtm_asdif, rrtm_asdir
#endif

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This subroutine controls the calls of the radiation schemes
!------------------------------------------------------------------------------!
    SUBROUTINE radiation_control


       IMPLICIT NONE


       IF ( debug_output_timestep )  CALL debug_message( 'radiation_control', 'start' )


       SELECT CASE ( TRIM( radiation_scheme ) )

          CASE ( 'constant' )
             CALL radiation_constant

          CASE ( 'clear-sky' )
             CALL radiation_clearsky

          CASE ( 'rrtmg' )
             CALL radiation_rrtmg

          CASE ( 'external' )
!
!--          During spinup apply clear-sky model
             IF ( time_since_reference_point < 0.0_wp )  THEN
                CALL radiation_clearsky
             ELSE
                CALL radiation_external
             ENDIF

          CASE DEFAULT

       END SELECT

       IF ( debug_output_timestep )  CALL debug_message( 'radiation_control', 'end' )

    END SUBROUTINE radiation_control

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check data output for radiation model
!------------------------------------------------------------------------------!
    SUBROUTINE radiation_check_data_output( variable, unit, i, ilen, k )


       USE control_parameters,                                                 &
           ONLY: data_output, message_string

       IMPLICIT NONE

       CHARACTER (LEN=*) ::  unit          !<
       CHARACTER (LEN=*) ::  variable      !<

       INTEGER(iwp) :: i, k
       INTEGER(iwp) :: ilen
       CHARACTER(LEN=varnamelength) :: var  !< TRIM(variable)

       var = TRIM(variable)

       IF ( len(var) < 3_iwp  )  THEN
          unit = 'illegal'
          RETURN
       ENDIF

       IF ( var(1:3) /= 'rad'  .AND.  var(1:3) /= 'rtm' )  THEN
          unit = 'illegal'
          RETURN
       ENDIF

!--    first process diractional variables
       IF ( var(1:12) == 'rtm_rad_net_'  .OR.  var(1:13) == 'rtm_rad_insw_'  .OR.        &
            var(1:13) == 'rtm_rad_inlw_'  .OR.  var(1:16) == 'rtm_rad_inswdir_'  .OR.    &
            var(1:16) == 'rtm_rad_inswdif_'  .OR.  var(1:16) == 'rtm_rad_inswref_'  .OR. &
            var(1:16) == 'rtm_rad_inlwdif_'  .OR.  var(1:16) == 'rtm_rad_inlwref_'  .OR. &
            var(1:14) == 'rtm_rad_outsw_'  .OR.  var(1:14) == 'rtm_rad_outlw_'  .OR.     &
            var(1:14) == 'rtm_rad_ressw_'  .OR.  var(1:14) == 'rtm_rad_reslw_'  ) THEN
          IF ( .NOT.  radiation ) THEN
                message_string = 'output of "' // TRIM( var ) // '" require'&
                                 // 's radiation = .TRUE.'
                CALL message( 'check_parameters', 'PA0509', 1, 2, 0, 6, 0 )
          ENDIF
          unit = 'W/m2'
       ELSE IF ( var(1:7) == 'rtm_svf'  .OR.  var(1:7) == 'rtm_dif'  .OR.                &
                 var(1:9) == 'rtm_skyvf' .OR. var(1:9) == 'rtm_skyvft'  .OR.             &
                 var(1:12) == 'rtm_surfalb_'  .OR.  var(1:13) == 'rtm_surfemis_'  ) THEN
          IF ( .NOT.  radiation ) THEN
                message_string = 'output of "' // TRIM( var ) // '" require'&
                                 // 's radiation = .TRUE.'
                CALL message( 'check_parameters', 'PA0509', 1, 2, 0, 6, 0 )
          ENDIF
          unit = '1'
       ELSE
!--       non-directional variables
          SELECT CASE ( TRIM( var ) )
             CASE ( 'rad_lw_cs_hr', 'rad_lw_hr', 'rad_lw_in', 'rad_lw_out', &
                    'rad_sw_cs_hr', 'rad_sw_hr', 'rad_sw_in', 'rad_sw_out'  )
                IF (  .NOT.  radiation  .OR.  radiation_scheme /= 'rrtmg' )  THEN
                   message_string = '"output of "' // TRIM( var ) // '" requi' // &
                                    'res radiation = .TRUE. and ' //              &
                                    'radiation_scheme = "rrtmg"'
                   CALL message( 'check_parameters', 'PA0406', 1, 2, 0, 6, 0 )
                ENDIF
                unit = 'K/h'

             CASE ( 'rad_net*', 'rrtm_aldif*', 'rrtm_aldir*', 'rrtm_asdif*',      &
                    'rrtm_asdir*', 'rad_lw_in*', 'rad_lw_out*', 'rad_sw_in*',     &
                    'rad_sw_out*')
                IF ( i == 0 .AND. ilen == 0 .AND. k == 0)  THEN
                   ! Workaround for masked output (calls with i=ilen=k=0)
                   unit = 'illegal'
                   RETURN
                ENDIF
                IF ( k == 0  .OR.  data_output(i)(ilen-2:ilen) /= '_xy' )  THEN
                   message_string = 'illegal value for data_output: "' //         &
                                    TRIM( var ) // '" & only 2d-horizontal ' //   &
                                    'cross sections are allowed for this value'
                   CALL message( 'check_parameters', 'PA0111', 1, 2, 0, 6, 0 )
                ENDIF
                IF (  .NOT.  radiation  .OR.  radiation_scheme /= "rrtmg" )  THEN
                   IF ( TRIM( var ) == 'rrtm_aldif*'  .OR.                        &
                        TRIM( var ) == 'rrtm_aldir*'  .OR.                        &
                        TRIM( var ) == 'rrtm_asdif*'  .OR.                        &
                        TRIM( var ) == 'rrtm_asdir*'      )                       &
                   THEN
                      message_string = 'output of "' // TRIM( var ) // '" require'&
                                       // 's radiation = .TRUE. and radiation_sch'&
                                       // 'eme = "rrtmg"'
                      CALL message( 'check_parameters', 'PA0409', 1, 2, 0, 6, 0 )
                   ENDIF
                ENDIF

                IF ( TRIM( var ) == 'rad_net*'      ) unit = 'W/m2'
                IF ( TRIM( var ) == 'rad_lw_in*'    ) unit = 'W/m2'
                IF ( TRIM( var ) == 'rad_lw_out*'   ) unit = 'W/m2'
                IF ( TRIM( var ) == 'rad_sw_in*'    ) unit = 'W/m2'
                IF ( TRIM( var ) == 'rad_sw_out*'   ) unit = 'W/m2'
                IF ( TRIM( var ) == 'rad_sw_in'     ) unit = 'W/m2'
                IF ( TRIM( var ) == 'rrtm_aldif*'   ) unit = ''
                IF ( TRIM( var ) == 'rrtm_aldir*'   ) unit = ''
                IF ( TRIM( var ) == 'rrtm_asdif*'   ) unit = ''
                IF ( TRIM( var ) == 'rrtm_asdir*'   ) unit = ''

             CASE ( 'rtm_rad_pc_inlw', 'rtm_rad_pc_insw', 'rtm_rad_pc_inswdir', &
                    'rtm_rad_pc_inswdif', 'rtm_rad_pc_inswref')
                IF ( .NOT.  radiation ) THEN
                   message_string = 'output of "' // TRIM( var ) // '" require'&
                                    // 's radiation = .TRUE.'
                   CALL message( 'check_parameters', 'PA0509', 1, 2, 0, 6, 0 )
                ENDIF
                unit = 'W'

             CASE ( 'rtm_mrt', 'rtm_mrt_sw', 'rtm_mrt_lw'  )
                IF ( i == 0 .AND. ilen == 0 .AND. k == 0)  THEN
                   ! Workaround for masked output (calls with i=ilen=k=0)
                   unit = 'illegal'
                   RETURN
                ENDIF

                IF ( .NOT.  radiation ) THEN
                   message_string = 'output of "' // TRIM( var ) // '" require'&
                                    // 's radiation = .TRUE.'
                   CALL message( 'check_parameters', 'PA0509', 1, 2, 0, 6, 0 )
                ENDIF
                IF ( mrt_nlevels == 0 ) THEN
                   message_string = 'output of "' // TRIM( var ) // '" require'&
                                    // 's mrt_nlevels > 0'
                   CALL message( 'check_parameters', 'PA0510', 1, 2, 0, 6, 0 )
                ENDIF
                IF ( TRIM( var ) == 'rtm_mrt_sw'  .AND.  .NOT. mrt_include_sw ) THEN
                   message_string = 'output of "' // TRIM( var ) // '" require'&
                                    // 's rtm_mrt_sw = .TRUE.'
                   CALL message( 'check_parameters', 'PA0511', 1, 2, 0, 6, 0 )
                ENDIF
                IF ( TRIM( var ) == 'rtm_mrt' ) THEN
                   unit = 'K'
                ELSE
                   unit = 'W m-2'
                ENDIF

             CASE DEFAULT
                unit = 'illegal'

          END SELECT
       ENDIF

    END SUBROUTINE radiation_check_data_output


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Set module-specific timeseries units and labels
!------------------------------------------------------------------------------!
 SUBROUTINE radiation_check_data_output_ts( dots_max, dots_num )


    INTEGER(iwp),      INTENT(IN)     ::  dots_max
    INTEGER(iwp),      INTENT(INOUT)  ::  dots_num

!
!-- Next line is just to avoid compiler warning about unused variable.
    IF ( dots_max == 0 )  CONTINUE

!
!-- Temporary solution to add LSM and radiation time series to the default
!-- output
    IF ( land_surface  .OR.  radiation )  THEN
       IF ( TRIM( radiation_scheme ) == 'rrtmg' )  THEN
          dots_num = dots_num + 15
       ELSE
          dots_num = dots_num + 11
       ENDIF
    ENDIF


 END SUBROUTINE radiation_check_data_output_ts

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check data output of profiles for radiation model
!------------------------------------------------------------------------------!
    SUBROUTINE radiation_check_data_output_pr( variable, var_count, unit,      &
               dopr_unit )

       USE arrays_3d,                                                          &
           ONLY: zu

       USE control_parameters,                                                 &
           ONLY: data_output_pr, message_string

       USE indices

       USE profil_parameter

       USE statistics

       IMPLICIT NONE

       CHARACTER (LEN=*) ::  unit      !<
       CHARACTER (LEN=*) ::  variable  !<
       CHARACTER (LEN=*) ::  dopr_unit !< local value of dopr_unit

       INTEGER(iwp) ::  var_count     !<

       SELECT CASE ( TRIM( variable ) )

         CASE ( 'rad_net' )
             IF ( (  .NOT.  radiation )  .OR.  radiation_scheme == 'constant' )&
             THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) // ' is' // &
                                 'not available for radiation = .FALSE. or ' //&
                                 'radiation_scheme = "constant"'
                CALL message( 'check_parameters', 'PA0408', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 99
                dopr_unit  = 'W/m2'
                hom(:,2,99,:)  = SPREAD( zw, 2, statistic_regions+1 )
                unit = dopr_unit
             ENDIF

          CASE ( 'rad_lw_in' )
             IF ( (  .NOT.  radiation)  .OR.  radiation_scheme == 'constant' ) &
             THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) // ' is' // &
                                 'not available for radiation = .FALSE. or ' //&
                                 'radiation_scheme = "constant"'
                CALL message( 'check_parameters', 'PA0408', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 100
                dopr_unit  = 'W/m2'
                hom(:,2,100,:)  = SPREAD( zw, 2, statistic_regions+1 )
                unit = dopr_unit
             ENDIF

          CASE ( 'rad_lw_out' )
             IF ( (  .NOT. radiation )  .OR.  radiation_scheme == 'constant' ) &
             THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) // ' is' // &
                                 'not available for radiation = .FALSE. or ' //&
                                 'radiation_scheme = "constant"'
                CALL message( 'check_parameters', 'PA0408', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 101
                dopr_unit  = 'W/m2'
                hom(:,2,101,:)  = SPREAD( zw, 2, statistic_regions+1 )
                unit = dopr_unit
             ENDIF

          CASE ( 'rad_sw_in' )
             IF ( (  .NOT. radiation )  .OR.  radiation_scheme == 'constant' ) &
             THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) // ' is' // &
                                 'not available for radiation = .FALSE. or ' //&
                                 'radiation_scheme = "constant"'
                CALL message( 'check_parameters', 'PA0408', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 102
                dopr_unit  = 'W/m2'
                hom(:,2,102,:)  = SPREAD( zw, 2, statistic_regions+1 )
                unit = dopr_unit
             ENDIF

          CASE ( 'rad_sw_out')
             IF ( (  .NOT.  radiation )  .OR.  radiation_scheme == 'constant' )&
             THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) // ' is' // &
                                 'not available for radiation = .FALSE. or ' //&
                                 'radiation_scheme = "constant"'
                CALL message( 'check_parameters', 'PA0408', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 103
                dopr_unit  = 'W/m2'
                hom(:,2,103,:)  = SPREAD( zw, 2, statistic_regions+1 )
                unit = dopr_unit
             ENDIF

          CASE ( 'rad_lw_cs_hr' )
             IF ( (  .NOT.  radiation )  .OR.  radiation_scheme /= 'rrtmg' )   &
             THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) // ' is' // &
                                 'not available for radiation = .FALSE. or ' //&
                                 'radiation_scheme /= "rrtmg"'
                CALL message( 'check_parameters', 'PA0413', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 104
                dopr_unit  = 'K/h'
                hom(:,2,104,:)  = SPREAD( zu, 2, statistic_regions+1 )
                unit = dopr_unit
             ENDIF

          CASE ( 'rad_lw_hr' )
             IF ( (  .NOT.  radiation )  .OR.  radiation_scheme /= 'rrtmg' )   &
             THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) // ' is' // &
                                 'not available for radiation = .FALSE. or ' //&
                                 'radiation_scheme /= "rrtmg"'
                CALL message( 'check_parameters', 'PA0413', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 105
                dopr_unit  = 'K/h'
                hom(:,2,105,:)  = SPREAD( zu, 2, statistic_regions+1 )
                unit = dopr_unit
             ENDIF

          CASE ( 'rad_sw_cs_hr' )
             IF ( (  .NOT.  radiation )  .OR.  radiation_scheme /= 'rrtmg' )   &
             THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) // ' is' // &
                                 'not available for radiation = .FALSE. or ' //&
                                 'radiation_scheme /= "rrtmg"'
                CALL message( 'check_parameters', 'PA0413', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 106
                dopr_unit  = 'K/h'
                hom(:,2,106,:)  = SPREAD( zu, 2, statistic_regions+1 )
                unit = dopr_unit
             ENDIF

          CASE ( 'rad_sw_hr' )
             IF ( (  .NOT.  radiation )  .OR.  radiation_scheme /= 'rrtmg' )   &
             THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) // ' is' // &
                                 'not available for radiation = .FALSE. or ' //&
                                 'radiation_scheme /= "rrtmg"'
                CALL message( 'check_parameters', 'PA0413', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 107
                dopr_unit  = 'K/h'
                hom(:,2,107,:)  = SPREAD( zu, 2, statistic_regions+1 )
                unit = dopr_unit
             ENDIF


          CASE DEFAULT
             unit = 'illegal'

       END SELECT


    END SUBROUTINE radiation_check_data_output_pr


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check parameters routine for radiation model
!------------------------------------------------------------------------------!
    SUBROUTINE radiation_check_parameters

       USE control_parameters,                                                 &
           ONLY: land_surface, message_string, rotation_angle, urban_surface

       USE netcdf_data_input_mod,                                              &
           ONLY:  input_pids_static

       IMPLICIT NONE

!
!--    In case no urban-surface or land-surface model is applied, usage of
!--    a radiation model make no sense.
       IF ( .NOT. land_surface  .AND.  .NOT. urban_surface )  THEN
          message_string = 'Usage of radiation module is only allowed if ' //  &
                           'land-surface and/or urban-surface model is applied.'
          CALL message( 'check_parameters', 'PA0486', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( radiation_scheme /= 'constant'   .AND.                             &
            radiation_scheme /= 'clear-sky'  .AND.                             &
            radiation_scheme /= 'rrtmg'      .AND.                             &
            radiation_scheme /= 'external' )  THEN
          message_string = 'unknown radiation_scheme = '//                     &
                           TRIM( radiation_scheme )
          CALL message( 'check_parameters', 'PA0405', 1, 2, 0, 6, 0 )
       ELSEIF ( radiation_scheme == 'rrtmg' )  THEN
#if ! defined ( __rrtmg )
          message_string = 'radiation_scheme = "rrtmg" requires ' //           &
                           'compilation of PALM with pre-processor ' //        &
                           'directive -D__rrtmg'
          CALL message( 'check_parameters', 'PA0407', 1, 2, 0, 6, 0 )
#endif
#if defined ( __rrtmg ) && ! defined( __netcdf )
          message_string = 'radiation_scheme = "rrtmg" requires ' //           &
                           'the use of NetCDF (preprocessor directive ' //     &
                           '-D__netcdf'
          CALL message( 'check_parameters', 'PA0412', 1, 2, 0, 6, 0 )
#endif

       ENDIF
!
!--    Checks performed only if data is given via namelist only.
       IF ( .NOT. input_pids_static )  THEN
          IF ( albedo_type == 0  .AND.  albedo == 9999999.9_wp  .AND.          &
               radiation_scheme == 'clear-sky')  THEN
             message_string = 'radiation_scheme = "clear-sky" in combination'//&
                              'with albedo_type = 0 requires setting of'//     &
                              'albedo /= 9999999.9'
             CALL message( 'check_parameters', 'PA0410', 1, 2, 0, 6, 0 )
          ENDIF

          IF ( albedo_type == 0  .AND.  radiation_scheme == 'rrtmg'  .AND.     &
             ( albedo_lw_dif == 9999999.9_wp .OR. albedo_lw_dir == 9999999.9_wp&
          .OR. albedo_sw_dif == 9999999.9_wp .OR. albedo_sw_dir == 9999999.9_wp&
             ) ) THEN
             message_string = 'radiation_scheme = "rrtmg" in combination' //   &
                              'with albedo_type = 0 requires setting of ' //   &
                              'albedo_lw_dif /= 9999999.9' //                  &
                              'albedo_lw_dir /= 9999999.9' //                  &
                              'albedo_sw_dif /= 9999999.9 and' //              &
                              'albedo_sw_dir /= 9999999.9'
             CALL message( 'check_parameters', 'PA0411', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF
!
!--    Parallel rad_angular_discretization without raytrace_mpi_rma is not implemented
!--    Serial mode does not allow mpi_rma
#if defined( __parallel )
       IF ( rad_angular_discretization  .AND.  .NOT. raytrace_mpi_rma )  THEN
          message_string = 'rad_angular_discretization can only be used ' //  &
                           'together with raytrace_mpi_rma or when ' //  &
                           'no parallelization is applied.'
          CALL message( 'readiation_check_parameters', 'PA0486', 1, 2, 0, 6, 0 )
       ENDIF
#else
       IF ( raytrace_mpi_rma )  THEN
          message_string = 'raytrace_mpi_rma = .T. not allowed in serial mode'
          CALL message( 'readiation_check_parameters', 'PA0710', 1, 2, 0, 6, 0 )
       ENDIF
#endif

       IF ( cloud_droplets  .AND.   radiation_scheme == 'rrtmg'  .AND.         &
            average_radiation ) THEN
          message_string = 'average_radiation = .T. with radiation_scheme'//   &
                           '= "rrtmg" in combination cloud_droplets = .T.'//   &
                           'is not implementd'
          CALL message( 'check_parameters', 'PA0560', 1, 2, 0, 6, 0 )
       ENDIF

!
!--    Incialize svf normalization reporting histogram
       svfnorm_report_num = 1
       DO WHILE ( svfnorm_report_thresh(svfnorm_report_num) < 1e20_wp          &
                   .AND. svfnorm_report_num <= 30 )
          svfnorm_report_num = svfnorm_report_num + 1
       ENDDO
       svfnorm_report_num = svfnorm_report_num - 1
!
!--    Check for dt_radiation
       IF ( dt_radiation <= 0.0 )  THEN
          message_string = 'dt_radiation must be > 0.0'
          CALL message( 'check_parameters', 'PA0591', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Check rotation angle
       !> @todo Remove this limitation
       IF ( rotation_angle /= 0.0 )  THEN
          message_string = 'rotation of the model domain is not considered in the radiation ' //   &
                           'model.&Using rotation_angle /= 0.0 is not allowed in combination ' //  &
                           'with the radiation model at the moment!'
          CALL message( 'check_parameters', 'PA0675', 1, 2, 0, 6, 0 )
       ENDIF

    END SUBROUTINE radiation_check_parameters


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of the radiation model and Radiative Transfer Model
!------------------------------------------------------------------------------!
    SUBROUTINE radiation_init

       IMPLICIT NONE

       INTEGER(iwp) ::  i         !< running index x-direction
       INTEGER(iwp) ::  is        !< running index for input surface elements
       INTEGER(iwp) ::  ioff      !< offset in x between surface element reference grid point in atmosphere and actual surface
       INTEGER(iwp) ::  j         !< running index y-direction
       INTEGER(iwp) ::  joff      !< offset in y between surface element reference grid point in atmosphere and actual surface
       INTEGER(iwp) ::  k         !< running index z-direction
       INTEGER(iwp) ::  l         !< running index for orientation of vertical surfaces
       INTEGER(iwp) ::  m         !< running index for surface elements
       INTEGER(iwp) ::  ntime = 0 !< number of available external radiation timesteps
#if defined( __rrtmg )
       INTEGER(iwp) ::  ind_type  !< running index for subgrid-surface tiles
#endif
       LOGICAL      ::  radiation_input_root_domain !< flag indicating the existence of a dynamic input file for the root domain


       IF ( debug_output )  CALL debug_message( 'radiation_init', 'start' )
!
!--    Activate radiation_interactions according to the existence of vertical surfaces and/or trees
!      or if biometeorology output is required for flat surfaces.
!--    The namelist parameter radiation_interactions_on can override this behavior.
!--    (This check cannot be performed in check_parameters, because vertical_surfaces_exist is first set in
!--    init_surface_arrays.)
       IF ( radiation_interactions_on )  THEN
          IF ( vertical_surfaces_exist  .OR.  plant_canopy  .OR.  biometeorology )  THEN
             radiation_interactions    = .TRUE.
             average_radiation         = .TRUE.
          ELSE
             radiation_interactions_on = .FALSE.   !< reset namelist parameter: no interactions
                                                   !< calculations necessary in case of flat surface
          ENDIF
       ELSEIF ( vertical_surfaces_exist  .OR.  plant_canopy  .OR.  biometeorology )  THEN
          message_string = 'radiation_interactions_on is set to .FALSE. although '     // &
                           'vertical surfaces and/or trees or biometeorology exist '   // &
                           'is ON. The model will run without RTM (no shadows, no '    // &
                           'radiation reflections)'
          CALL message( 'init_3d_model', 'PA0348', 0, 1, 0, 6, 0 )
       ENDIF
!
!--    Precalculate some time constants
       d_hours_day    = 1.0_wp / REAL( hours_per_day, KIND = wp )
       d_seconds_hour = 1.0_wp / seconds_per_hour

!
!--    If required, initialize radiation interactions between surfaces
!--    via sky-view factors. This must be done before radiation is initialized.
       IF ( radiation_interactions )  CALL radiation_interaction_init
!
!--    Allocate array for storing the surface net radiation
       IF ( .NOT. ALLOCATED ( surf_lsm_h%rad_net )  .AND.                      &
                  surf_lsm_h%ns > 0  )   THEN
          ALLOCATE( surf_lsm_h%rad_net(1:surf_lsm_h%ns) )
          surf_lsm_h%rad_net = 0.0_wp
       ENDIF
       IF ( .NOT. ALLOCATED ( surf_usm_h%rad_net )  .AND.                      &
                  surf_usm_h%ns > 0  )  THEN
          ALLOCATE( surf_usm_h%rad_net(1:surf_usm_h%ns) )
          surf_usm_h%rad_net = 0.0_wp
       ENDIF
       DO  l = 0, 3
          IF ( .NOT. ALLOCATED ( surf_lsm_v(l)%rad_net )  .AND.                &
                     surf_lsm_v(l)%ns > 0  )  THEN
             ALLOCATE( surf_lsm_v(l)%rad_net(1:surf_lsm_v(l)%ns) )
             surf_lsm_v(l)%rad_net = 0.0_wp
          ENDIF
          IF ( .NOT. ALLOCATED ( surf_usm_v(l)%rad_net )  .AND.                &
                     surf_usm_v(l)%ns > 0  )  THEN
             ALLOCATE( surf_usm_v(l)%rad_net(1:surf_usm_v(l)%ns) )
             surf_usm_v(l)%rad_net = 0.0_wp
          ENDIF
       ENDDO


!
!--    Allocate array for storing the surface longwave (out) radiation change
       IF ( .NOT. ALLOCATED ( surf_lsm_h%rad_lw_out_change_0 )  .AND.          &
                  surf_lsm_h%ns > 0  )   THEN
          ALLOCATE( surf_lsm_h%rad_lw_out_change_0(1:surf_lsm_h%ns) )
          surf_lsm_h%rad_lw_out_change_0 = 0.0_wp
       ENDIF
       IF ( .NOT. ALLOCATED ( surf_usm_h%rad_lw_out_change_0 )  .AND.          &
                  surf_usm_h%ns > 0  )  THEN
          ALLOCATE( surf_usm_h%rad_lw_out_change_0(1:surf_usm_h%ns) )
          surf_usm_h%rad_lw_out_change_0 = 0.0_wp
       ENDIF
       DO  l = 0, 3
          IF ( .NOT. ALLOCATED ( surf_lsm_v(l)%rad_lw_out_change_0 )  .AND.    &
                     surf_lsm_v(l)%ns > 0  )  THEN
             ALLOCATE( surf_lsm_v(l)%rad_lw_out_change_0(1:surf_lsm_v(l)%ns) )
             surf_lsm_v(l)%rad_lw_out_change_0 = 0.0_wp
          ENDIF
          IF ( .NOT. ALLOCATED ( surf_usm_v(l)%rad_lw_out_change_0 )  .AND.    &
                     surf_usm_v(l)%ns > 0  )  THEN
             ALLOCATE( surf_usm_v(l)%rad_lw_out_change_0(1:surf_usm_v(l)%ns) )
             surf_usm_v(l)%rad_lw_out_change_0 = 0.0_wp
          ENDIF
       ENDDO

!
!--    Allocate surface arrays for incoming/outgoing short/longwave radiation
       IF ( .NOT. ALLOCATED ( surf_lsm_h%rad_sw_in )  .AND.                    &
                  surf_lsm_h%ns > 0  )   THEN
          ALLOCATE( surf_lsm_h%rad_sw_in(1:surf_lsm_h%ns)  )
          ALLOCATE( surf_lsm_h%rad_sw_out(1:surf_lsm_h%ns) )
          ALLOCATE( surf_lsm_h%rad_sw_dir(1:surf_lsm_h%ns) )
          ALLOCATE( surf_lsm_h%rad_sw_dif(1:surf_lsm_h%ns) )
          ALLOCATE( surf_lsm_h%rad_sw_ref(1:surf_lsm_h%ns) )
          ALLOCATE( surf_lsm_h%rad_sw_res(1:surf_lsm_h%ns) )
          ALLOCATE( surf_lsm_h%rad_lw_in(1:surf_lsm_h%ns)  )
          ALLOCATE( surf_lsm_h%rad_lw_out(1:surf_lsm_h%ns) )
          ALLOCATE( surf_lsm_h%rad_lw_dif(1:surf_lsm_h%ns) )
          ALLOCATE( surf_lsm_h%rad_lw_ref(1:surf_lsm_h%ns) )
          ALLOCATE( surf_lsm_h%rad_lw_res(1:surf_lsm_h%ns) )
          surf_lsm_h%rad_sw_in  = 0.0_wp
          surf_lsm_h%rad_sw_out = 0.0_wp
          surf_lsm_h%rad_sw_dir = 0.0_wp
          surf_lsm_h%rad_sw_dif = 0.0_wp
          surf_lsm_h%rad_sw_ref = 0.0_wp
          surf_lsm_h%rad_sw_res = 0.0_wp
          surf_lsm_h%rad_lw_in  = 0.0_wp
          surf_lsm_h%rad_lw_out = 0.0_wp
          surf_lsm_h%rad_lw_dif = 0.0_wp
          surf_lsm_h%rad_lw_ref = 0.0_wp
          surf_lsm_h%rad_lw_res = 0.0_wp
       ENDIF
       IF ( .NOT. ALLOCATED ( surf_usm_h%rad_sw_in )  .AND.                    &
                  surf_usm_h%ns > 0  )  THEN
          ALLOCATE( surf_usm_h%rad_sw_in(1:surf_usm_h%ns)  )
          ALLOCATE( surf_usm_h%rad_sw_out(1:surf_usm_h%ns) )
          ALLOCATE( surf_usm_h%rad_sw_dir(1:surf_usm_h%ns) )
          ALLOCATE( surf_usm_h%rad_sw_dif(1:surf_usm_h%ns) )
          ALLOCATE( surf_usm_h%rad_sw_ref(1:surf_usm_h%ns) )
          ALLOCATE( surf_usm_h%rad_sw_res(1:surf_usm_h%ns) )
          ALLOCATE( surf_usm_h%rad_lw_in(1:surf_usm_h%ns)  )
          ALLOCATE( surf_usm_h%rad_lw_out(1:surf_usm_h%ns) )
          ALLOCATE( surf_usm_h%rad_lw_dif(1:surf_usm_h%ns) )
          ALLOCATE( surf_usm_h%rad_lw_ref(1:surf_usm_h%ns) )
          ALLOCATE( surf_usm_h%rad_lw_res(1:surf_usm_h%ns) )
          surf_usm_h%rad_sw_in  = 0.0_wp
          surf_usm_h%rad_sw_out = 0.0_wp
          surf_usm_h%rad_sw_dir = 0.0_wp
          surf_usm_h%rad_sw_dif = 0.0_wp
          surf_usm_h%rad_sw_ref = 0.0_wp
          surf_usm_h%rad_sw_res = 0.0_wp
          surf_usm_h%rad_lw_in  = 0.0_wp
          surf_usm_h%rad_lw_out = 0.0_wp
          surf_usm_h%rad_lw_dif = 0.0_wp
          surf_usm_h%rad_lw_ref = 0.0_wp
          surf_usm_h%rad_lw_res = 0.0_wp
       ENDIF
       DO  l = 0, 3
          IF ( .NOT. ALLOCATED ( surf_lsm_v(l)%rad_sw_in )  .AND.              &
                     surf_lsm_v(l)%ns > 0  )  THEN
             ALLOCATE( surf_lsm_v(l)%rad_sw_in(1:surf_lsm_v(l)%ns)  )
             ALLOCATE( surf_lsm_v(l)%rad_sw_out(1:surf_lsm_v(l)%ns) )
             ALLOCATE( surf_lsm_v(l)%rad_sw_dir(1:surf_lsm_v(l)%ns) )
             ALLOCATE( surf_lsm_v(l)%rad_sw_dif(1:surf_lsm_v(l)%ns) )
             ALLOCATE( surf_lsm_v(l)%rad_sw_ref(1:surf_lsm_v(l)%ns) )
             ALLOCATE( surf_lsm_v(l)%rad_sw_res(1:surf_lsm_v(l)%ns) )

             ALLOCATE( surf_lsm_v(l)%rad_lw_in(1:surf_lsm_v(l)%ns)  )
             ALLOCATE( surf_lsm_v(l)%rad_lw_out(1:surf_lsm_v(l)%ns) )
             ALLOCATE( surf_lsm_v(l)%rad_lw_dif(1:surf_lsm_v(l)%ns) )
             ALLOCATE( surf_lsm_v(l)%rad_lw_ref(1:surf_lsm_v(l)%ns) )
             ALLOCATE( surf_lsm_v(l)%rad_lw_res(1:surf_lsm_v(l)%ns) )

             surf_lsm_v(l)%rad_sw_in  = 0.0_wp
             surf_lsm_v(l)%rad_sw_out = 0.0_wp
             surf_lsm_v(l)%rad_sw_dir = 0.0_wp
             surf_lsm_v(l)%rad_sw_dif = 0.0_wp
             surf_lsm_v(l)%rad_sw_ref = 0.0_wp
             surf_lsm_v(l)%rad_sw_res = 0.0_wp

             surf_lsm_v(l)%rad_lw_in  = 0.0_wp
             surf_lsm_v(l)%rad_lw_out = 0.0_wp
             surf_lsm_v(l)%rad_lw_dif = 0.0_wp
             surf_lsm_v(l)%rad_lw_ref = 0.0_wp
             surf_lsm_v(l)%rad_lw_res = 0.0_wp
          ENDIF
          IF ( .NOT. ALLOCATED ( surf_usm_v(l)%rad_sw_in )  .AND.              &
                     surf_usm_v(l)%ns > 0  )  THEN
             ALLOCATE( surf_usm_v(l)%rad_sw_in(1:surf_usm_v(l)%ns)  )
             ALLOCATE( surf_usm_v(l)%rad_sw_out(1:surf_usm_v(l)%ns) )
             ALLOCATE( surf_usm_v(l)%rad_sw_dir(1:surf_usm_v(l)%ns) )
             ALLOCATE( surf_usm_v(l)%rad_sw_dif(1:surf_usm_v(l)%ns) )
             ALLOCATE( surf_usm_v(l)%rad_sw_ref(1:surf_usm_v(l)%ns) )
             ALLOCATE( surf_usm_v(l)%rad_sw_res(1:surf_usm_v(l)%ns) )
             ALLOCATE( surf_usm_v(l)%rad_lw_in(1:surf_usm_v(l)%ns)  )
             ALLOCATE( surf_usm_v(l)%rad_lw_out(1:surf_usm_v(l)%ns) )
             ALLOCATE( surf_usm_v(l)%rad_lw_dif(1:surf_usm_v(l)%ns) )
             ALLOCATE( surf_usm_v(l)%rad_lw_ref(1:surf_usm_v(l)%ns) )
             ALLOCATE( surf_usm_v(l)%rad_lw_res(1:surf_usm_v(l)%ns) )
             surf_usm_v(l)%rad_sw_in  = 0.0_wp
             surf_usm_v(l)%rad_sw_out = 0.0_wp
             surf_usm_v(l)%rad_sw_dir = 0.0_wp
             surf_usm_v(l)%rad_sw_dif = 0.0_wp
             surf_usm_v(l)%rad_sw_ref = 0.0_wp
             surf_usm_v(l)%rad_sw_res = 0.0_wp
             surf_usm_v(l)%rad_lw_in  = 0.0_wp
             surf_usm_v(l)%rad_lw_out = 0.0_wp
             surf_usm_v(l)%rad_lw_dif = 0.0_wp
             surf_usm_v(l)%rad_lw_ref = 0.0_wp
             surf_usm_v(l)%rad_lw_res = 0.0_wp
          ENDIF
       ENDDO
!
!--    Fix net radiation in case of radiation_scheme = 'constant'
       IF ( radiation_scheme == 'constant' )  THEN
          IF ( ALLOCATED( surf_lsm_h%rad_net ) )                               &
             surf_lsm_h%rad_net    = net_radiation
          IF ( ALLOCATED( surf_usm_h%rad_net ) )                               &
             surf_usm_h%rad_net    = net_radiation
!
!--       Todo: weight with inclination angle
          DO  l = 0, 3
             IF ( ALLOCATED( surf_lsm_v(l)%rad_net ) )                         &
                surf_lsm_v(l)%rad_net = net_radiation
             IF ( ALLOCATED( surf_usm_v(l)%rad_net ) )                         &
                surf_usm_v(l)%rad_net = net_radiation
          ENDDO
!          radiation = .FALSE.
!
!--    Calculate orbital constants
       ELSE
          decl_1 = SIN(23.45_wp * pi / 180.0_wp)
          decl_2 = 2.0_wp * pi / 365.0_wp
          decl_3 = decl_2 * 81.0_wp
          lat    = latitude * pi / 180.0_wp
          lon    = longitude * pi / 180.0_wp
       ENDIF

       IF ( radiation_scheme == 'clear-sky'  .OR.                              &
            radiation_scheme == 'constant'   .OR.                              &
            radiation_scheme == 'external' )  THEN
!
!--       Allocate arrays for incoming/outgoing short/longwave radiation
          IF ( .NOT. ALLOCATED ( rad_sw_in ) )  THEN
             ALLOCATE ( rad_sw_in(0:0,nysg:nyng,nxlg:nxrg) )
             rad_sw_in = 0.0_wp
          ENDIF
          IF ( .NOT. ALLOCATED ( rad_sw_out ) )  THEN
             ALLOCATE ( rad_sw_out(0:0,nysg:nyng,nxlg:nxrg) )
             rad_sw_out = 0.0_wp
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_lw_in ) )  THEN
             ALLOCATE ( rad_lw_in(0:0,nysg:nyng,nxlg:nxrg) )
             rad_lw_in = 0.0_wp
          ENDIF
          IF ( .NOT. ALLOCATED ( rad_lw_out ) )  THEN
             ALLOCATE ( rad_lw_out(0:0,nysg:nyng,nxlg:nxrg) )
             rad_lw_out = 0.0_wp
          ENDIF

!
!--       Allocate average arrays for incoming/outgoing short/longwave radiation
          IF ( .NOT. ALLOCATED ( rad_sw_in_av ) )  THEN
             ALLOCATE ( rad_sw_in_av(0:0,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( .NOT. ALLOCATED ( rad_sw_out_av ) )  THEN
             ALLOCATE ( rad_sw_out_av(0:0,nysg:nyng,nxlg:nxrg) )
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_lw_in_av ) )  THEN
             ALLOCATE ( rad_lw_in_av(0:0,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( .NOT. ALLOCATED ( rad_lw_out_av ) )  THEN
             ALLOCATE ( rad_lw_out_av(0:0,nysg:nyng,nxlg:nxrg) )
          ENDIF
!
!--       Allocate arrays for broadband albedo, and level 1 initialization
!--       via namelist paramter, unless not already allocated.
          IF ( .NOT. ALLOCATED(surf_lsm_h%albedo) )  THEN
             ALLOCATE( surf_lsm_h%albedo(1:surf_lsm_h%ns,0:2)     )
             surf_lsm_h%albedo    = albedo
          ENDIF
          IF ( .NOT. ALLOCATED(surf_usm_h%albedo) )  THEN
             ALLOCATE( surf_usm_h%albedo(1:surf_usm_h%ns,0:2)     )
             surf_usm_h%albedo    = albedo
          ENDIF

          DO  l = 0, 3
             IF ( .NOT. ALLOCATED( surf_lsm_v(l)%albedo ) )  THEN
                ALLOCATE( surf_lsm_v(l)%albedo(1:surf_lsm_v(l)%ns,0:2) )
                surf_lsm_v(l)%albedo = albedo
             ENDIF
             IF ( .NOT. ALLOCATED( surf_usm_v(l)%albedo ) )  THEN
                ALLOCATE( surf_usm_v(l)%albedo(1:surf_usm_v(l)%ns,0:2) )
                surf_usm_v(l)%albedo = albedo
             ENDIF
          ENDDO
!
!--       Level 2 initialization of broadband albedo via given albedo_type.
!--       Only if albedo_type is non-zero. In case of urban surface and
!--       input data is read from ASCII file, albedo_type will be zero, so that
!--       albedo won't be overwritten.
          DO  m = 1, surf_lsm_h%ns
             IF ( surf_lsm_h%albedo_type(m,ind_veg_wall) /= 0 )                &
                surf_lsm_h%albedo(m,ind_veg_wall) =                            &
                           albedo_pars(0,surf_lsm_h%albedo_type(m,ind_veg_wall))
             IF ( surf_lsm_h%albedo_type(m,ind_pav_green) /= 0 )               &
                surf_lsm_h%albedo(m,ind_pav_green) =                           &
                           albedo_pars(0,surf_lsm_h%albedo_type(m,ind_pav_green))
             IF ( surf_lsm_h%albedo_type(m,ind_wat_win) /= 0 )                 &
                surf_lsm_h%albedo(m,ind_wat_win) =                             &
                           albedo_pars(0,surf_lsm_h%albedo_type(m,ind_wat_win))
          ENDDO
          DO  m = 1, surf_usm_h%ns
             IF ( surf_usm_h%albedo_type(m,ind_veg_wall) /= 0 )                &
                surf_usm_h%albedo(m,ind_veg_wall) =                            &
                           albedo_pars(0,surf_usm_h%albedo_type(m,ind_veg_wall))
             IF ( surf_usm_h%albedo_type(m,ind_pav_green) /= 0 )               &
                surf_usm_h%albedo(m,ind_pav_green) =                           &
                           albedo_pars(0,surf_usm_h%albedo_type(m,ind_pav_green))
             IF ( surf_usm_h%albedo_type(m,ind_wat_win) /= 0 )                 &
                surf_usm_h%albedo(m,ind_wat_win) =                             &
                           albedo_pars(0,surf_usm_h%albedo_type(m,ind_wat_win))
          ENDDO

          DO  l = 0, 3
             DO  m = 1, surf_lsm_v(l)%ns
                IF ( surf_lsm_v(l)%albedo_type(m,ind_veg_wall) /= 0 )          &
                   surf_lsm_v(l)%albedo(m,ind_veg_wall) =                      &
                        albedo_pars(0,surf_lsm_v(l)%albedo_type(m,ind_veg_wall))
                IF ( surf_lsm_v(l)%albedo_type(m,ind_pav_green) /= 0 )         &
                   surf_lsm_v(l)%albedo(m,ind_pav_green) =                     &
                        albedo_pars(0,surf_lsm_v(l)%albedo_type(m,ind_pav_green))
                IF ( surf_lsm_v(l)%albedo_type(m,ind_wat_win) /= 0 )           &
                   surf_lsm_v(l)%albedo(m,ind_wat_win) =                       &
                        albedo_pars(0,surf_lsm_v(l)%albedo_type(m,ind_wat_win))
             ENDDO
             DO  m = 1, surf_usm_v(l)%ns
                IF ( surf_usm_v(l)%albedo_type(m,ind_veg_wall) /= 0 )          &
                   surf_usm_v(l)%albedo(m,ind_veg_wall) =                      &
                        albedo_pars(0,surf_usm_v(l)%albedo_type(m,ind_veg_wall))
                IF ( surf_usm_v(l)%albedo_type(m,ind_pav_green) /= 0 )         &
                   surf_usm_v(l)%albedo(m,ind_pav_green) =                     &
                        albedo_pars(0,surf_usm_v(l)%albedo_type(m,ind_pav_green))
                IF ( surf_usm_v(l)%albedo_type(m,ind_wat_win) /= 0 )           &
                   surf_usm_v(l)%albedo(m,ind_wat_win) =                       &
                        albedo_pars(0,surf_usm_v(l)%albedo_type(m,ind_wat_win))
             ENDDO
          ENDDO

!
!--       Level 3 initialization at grid points where albedo type is zero.
!--       This case, albedo is taken from file. In case of constant radiation
!--       or clear sky, only broadband albedo is given.
          IF ( albedo_pars_f%from_file )  THEN
!
!--          Horizontal surfaces
             DO  m = 1, surf_lsm_h%ns
                i = surf_lsm_h%i(m)
                j = surf_lsm_h%j(m)
                IF ( albedo_pars_f%pars_xy(0,j,i) /= albedo_pars_f%fill )  THEN
                   surf_lsm_h%albedo(m,ind_veg_wall)  = albedo_pars_f%pars_xy(0,j,i)
                   surf_lsm_h%albedo(m,ind_pav_green) = albedo_pars_f%pars_xy(0,j,i)
                   surf_lsm_h%albedo(m,ind_wat_win)   = albedo_pars_f%pars_xy(0,j,i)
                ENDIF
             ENDDO
             DO  m = 1, surf_usm_h%ns
                i = surf_usm_h%i(m)
                j = surf_usm_h%j(m)
                IF ( albedo_pars_f%pars_xy(0,j,i) /= albedo_pars_f%fill )  THEN
                   surf_usm_h%albedo(m,ind_veg_wall)  = albedo_pars_f%pars_xy(0,j,i)
                   surf_usm_h%albedo(m,ind_pav_green) = albedo_pars_f%pars_xy(0,j,i)
                   surf_usm_h%albedo(m,ind_wat_win)   = albedo_pars_f%pars_xy(0,j,i)
                ENDIF
             ENDDO
!
!--          Vertical surfaces
             DO  l = 0, 3

                ioff = surf_lsm_v(l)%ioff
                joff = surf_lsm_v(l)%joff
                DO  m = 1, surf_lsm_v(l)%ns
                   i = surf_lsm_v(l)%i(m) + ioff
                   j = surf_lsm_v(l)%j(m) + joff
                   IF ( albedo_pars_f%pars_xy(0,j,i) /= albedo_pars_f%fill )  THEN
                      surf_lsm_v(l)%albedo(m,ind_veg_wall) = albedo_pars_f%pars_xy(0,j,i)
                      surf_lsm_v(l)%albedo(m,ind_pav_green) = albedo_pars_f%pars_xy(0,j,i)
                      surf_lsm_v(l)%albedo(m,ind_wat_win) = albedo_pars_f%pars_xy(0,j,i)
                   ENDIF
                ENDDO

                ioff = surf_usm_v(l)%ioff
                joff = surf_usm_v(l)%joff
                DO  m = 1, surf_usm_v(l)%ns
                   i = surf_usm_v(l)%i(m) + ioff
                   j = surf_usm_v(l)%j(m) + joff
                   IF ( albedo_pars_f%pars_xy(0,j,i) /= albedo_pars_f%fill )  THEN
                      surf_usm_v(l)%albedo(m,ind_veg_wall) = albedo_pars_f%pars_xy(0,j,i)
                      surf_usm_v(l)%albedo(m,ind_pav_green) = albedo_pars_f%pars_xy(0,j,i)
                      surf_usm_v(l)%albedo(m,ind_wat_win) = albedo_pars_f%pars_xy(0,j,i)
                   ENDIF
                ENDDO
             ENDDO

          ENDIF
!
!--       Read explicit albedo values from building surface pars. If present,
!--       they override all less specific albedo values and force a albedo_type
!--       to zero in order to take effect.
          IF ( building_surface_pars_f%from_file )  THEN
             DO  m = 1, surf_usm_h%ns
                i = surf_usm_h%i(m)
                j = surf_usm_h%j(m)
                k = surf_usm_h%k(m)
!
!--             Iterate over surfaces in column, check height and orientation
                DO  is = building_surface_pars_f%index_ji(1,j,i), &
                         building_surface_pars_f%index_ji(2,j,i)
                   IF ( building_surface_pars_f%coords(4,is) == -surf_usm_h%koff .AND. &
                        building_surface_pars_f%coords(1,is) == k )  THEN

                      IF ( building_surface_pars_f%pars(ind_s_alb_b_wall,is) /=      &
                           building_surface_pars_f%fill )  THEN
                         surf_usm_h%albedo(m,ind_veg_wall) =                         &
                                  building_surface_pars_f%pars(ind_s_alb_b_wall,is)
                         surf_usm_h%albedo_type(m,ind_veg_wall) = 0
                      ENDIF

                      IF ( building_surface_pars_f%pars(ind_s_alb_b_win,is) /=       &
                           building_surface_pars_f%fill )  THEN
                         surf_usm_h%albedo(m,ind_wat_win) =                          &
                                  building_surface_pars_f%pars(ind_s_alb_b_win,is)
                         surf_usm_h%albedo_type(m,ind_wat_win) = 0
                      ENDIF

                      IF ( building_surface_pars_f%pars(ind_s_alb_b_green,is) /=     &
                           building_surface_pars_f%fill )  THEN
                         surf_usm_h%albedo(m,ind_pav_green) =                        &
                                  building_surface_pars_f%pars(ind_s_alb_b_green,is)
                         surf_usm_h%albedo_type(m,ind_pav_green) = 0
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
!--                Iterate over surfaces in column, check height and orientation
                   DO  is = building_surface_pars_f%index_ji(1,j,i), &
                            building_surface_pars_f%index_ji(2,j,i)
                      IF ( building_surface_pars_f%coords(5,is) == -surf_usm_v(l)%joff .AND. &
                           building_surface_pars_f%coords(6,is) == -surf_usm_v(l)%ioff .AND. &
                           building_surface_pars_f%coords(1,is) == k )  THEN

                         IF ( building_surface_pars_f%pars(ind_s_alb_b_wall,is) /=      &
                              building_surface_pars_f%fill )  THEN
                            surf_usm_v(l)%albedo(m,ind_veg_wall) =                      &
                                     building_surface_pars_f%pars(ind_s_alb_b_wall,is)
                            surf_usm_v(l)%albedo_type(m,ind_veg_wall) = 0
                         ENDIF

                         IF ( building_surface_pars_f%pars(ind_s_alb_b_win,is) /=       &
                              building_surface_pars_f%fill )  THEN
                            surf_usm_v(l)%albedo(m,ind_wat_win) =                       &
                                     building_surface_pars_f%pars(ind_s_alb_b_win,is)
                            surf_usm_v(l)%albedo_type(m,ind_wat_win) = 0
                         ENDIF

                         IF ( building_surface_pars_f%pars(ind_s_alb_b_green,is) /=     &
                              building_surface_pars_f%fill )  THEN
                            surf_usm_v(l)%albedo(m,ind_pav_green) =                     &
                                     building_surface_pars_f%pars(ind_s_alb_b_green,is)
                            surf_usm_v(l)%albedo_type(m,ind_pav_green) = 0
                         ENDIF

                         EXIT ! surface was found and processed
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
!
!--    Initialization actions for RRTMG
       ELSEIF ( radiation_scheme == 'rrtmg' )  THEN
#if defined ( __rrtmg )
!
!--       Allocate albedos for short/longwave radiation, horizontal surfaces
!--       for wall/green/window (USM) or vegetation/pavement/water surfaces
!--       (LSM).
          ALLOCATE ( surf_lsm_h%aldif(1:surf_lsm_h%ns,0:2)       )
          ALLOCATE ( surf_lsm_h%aldir(1:surf_lsm_h%ns,0:2)       )
          ALLOCATE ( surf_lsm_h%asdif(1:surf_lsm_h%ns,0:2)       )
          ALLOCATE ( surf_lsm_h%asdir(1:surf_lsm_h%ns,0:2)       )
          ALLOCATE ( surf_lsm_h%rrtm_aldif(1:surf_lsm_h%ns,0:2)  )
          ALLOCATE ( surf_lsm_h%rrtm_aldir(1:surf_lsm_h%ns,0:2)  )
          ALLOCATE ( surf_lsm_h%rrtm_asdif(1:surf_lsm_h%ns,0:2)  )
          ALLOCATE ( surf_lsm_h%rrtm_asdir(1:surf_lsm_h%ns,0:2)  )

          ALLOCATE ( surf_usm_h%aldif(1:surf_usm_h%ns,0:2)       )
          ALLOCATE ( surf_usm_h%aldir(1:surf_usm_h%ns,0:2)       )
          ALLOCATE ( surf_usm_h%asdif(1:surf_usm_h%ns,0:2)       )
          ALLOCATE ( surf_usm_h%asdir(1:surf_usm_h%ns,0:2)       )
          ALLOCATE ( surf_usm_h%rrtm_aldif(1:surf_usm_h%ns,0:2)  )
          ALLOCATE ( surf_usm_h%rrtm_aldir(1:surf_usm_h%ns,0:2)  )
          ALLOCATE ( surf_usm_h%rrtm_asdif(1:surf_usm_h%ns,0:2)  )
          ALLOCATE ( surf_usm_h%rrtm_asdir(1:surf_usm_h%ns,0:2)  )

!
!--       Allocate broadband albedo (temporary for the current radiation
!--       implementations)
          IF ( .NOT. ALLOCATED(surf_lsm_h%albedo) )                            &
             ALLOCATE( surf_lsm_h%albedo(1:surf_lsm_h%ns,0:2)     )
          IF ( .NOT. ALLOCATED(surf_usm_h%albedo) )                            &
             ALLOCATE( surf_usm_h%albedo(1:surf_usm_h%ns,0:2)     )

!
!--       Allocate albedos for short/longwave radiation, vertical surfaces
          DO  l = 0, 3

             ALLOCATE ( surf_lsm_v(l)%aldif(1:surf_lsm_v(l)%ns,0:2)      )
             ALLOCATE ( surf_lsm_v(l)%aldir(1:surf_lsm_v(l)%ns,0:2)      )
             ALLOCATE ( surf_lsm_v(l)%asdif(1:surf_lsm_v(l)%ns,0:2)      )
             ALLOCATE ( surf_lsm_v(l)%asdir(1:surf_lsm_v(l)%ns,0:2)      )

             ALLOCATE ( surf_lsm_v(l)%rrtm_aldif(1:surf_lsm_v(l)%ns,0:2) )
             ALLOCATE ( surf_lsm_v(l)%rrtm_aldir(1:surf_lsm_v(l)%ns,0:2) )
             ALLOCATE ( surf_lsm_v(l)%rrtm_asdif(1:surf_lsm_v(l)%ns,0:2) )
             ALLOCATE ( surf_lsm_v(l)%rrtm_asdir(1:surf_lsm_v(l)%ns,0:2) )

             ALLOCATE ( surf_usm_v(l)%aldif(1:surf_usm_v(l)%ns,0:2)      )
             ALLOCATE ( surf_usm_v(l)%aldir(1:surf_usm_v(l)%ns,0:2)      )
             ALLOCATE ( surf_usm_v(l)%asdif(1:surf_usm_v(l)%ns,0:2)      )
             ALLOCATE ( surf_usm_v(l)%asdir(1:surf_usm_v(l)%ns,0:2)      )

             ALLOCATE ( surf_usm_v(l)%rrtm_aldif(1:surf_usm_v(l)%ns,0:2) )
             ALLOCATE ( surf_usm_v(l)%rrtm_aldir(1:surf_usm_v(l)%ns,0:2) )
             ALLOCATE ( surf_usm_v(l)%rrtm_asdif(1:surf_usm_v(l)%ns,0:2) )
             ALLOCATE ( surf_usm_v(l)%rrtm_asdir(1:surf_usm_v(l)%ns,0:2) )
!
!--          Allocate broadband albedo (temporary for the current radiation
!--          implementations)
             IF ( .NOT. ALLOCATED( surf_lsm_v(l)%albedo ) )                    &
                ALLOCATE( surf_lsm_v(l)%albedo(1:surf_lsm_v(l)%ns,0:2) )
             IF ( .NOT. ALLOCATED( surf_usm_v(l)%albedo ) )                    &
                ALLOCATE( surf_usm_v(l)%albedo(1:surf_usm_v(l)%ns,0:2) )

          ENDDO
!
!--       Level 1 initialization of spectral albedos via namelist
!--       paramters. Please note, this case all surface tiles are initialized
!--       the same.
          IF ( surf_lsm_h%ns > 0 )  THEN
             surf_lsm_h%aldif  = albedo_lw_dif
             surf_lsm_h%aldir  = albedo_lw_dir
             surf_lsm_h%asdif  = albedo_sw_dif
             surf_lsm_h%asdir  = albedo_sw_dir
             surf_lsm_h%albedo = albedo_sw_dif
          ENDIF
          IF ( surf_usm_h%ns > 0 )  THEN
             IF ( surf_usm_h%albedo_from_ascii )  THEN
                surf_usm_h%aldif  = surf_usm_h%albedo
                surf_usm_h%aldir  = surf_usm_h%albedo
                surf_usm_h%asdif  = surf_usm_h%albedo
                surf_usm_h%asdir  = surf_usm_h%albedo
             ELSE
                surf_usm_h%aldif  = albedo_lw_dif
                surf_usm_h%aldir  = albedo_lw_dir
                surf_usm_h%asdif  = albedo_sw_dif
                surf_usm_h%asdir  = albedo_sw_dir
                surf_usm_h%albedo = albedo_sw_dif
             ENDIF
          ENDIF

          DO  l = 0, 3

             IF ( surf_lsm_v(l)%ns > 0 )  THEN
                surf_lsm_v(l)%aldif  = albedo_lw_dif
                surf_lsm_v(l)%aldir  = albedo_lw_dir
                surf_lsm_v(l)%asdif  = albedo_sw_dif
                surf_lsm_v(l)%asdir  = albedo_sw_dir
                surf_lsm_v(l)%albedo = albedo_sw_dif
             ENDIF

             IF ( surf_usm_v(l)%ns > 0 )  THEN
                IF ( surf_usm_v(l)%albedo_from_ascii )  THEN
                   surf_usm_v(l)%aldif  = surf_usm_v(l)%albedo
                   surf_usm_v(l)%aldir  = surf_usm_v(l)%albedo
                   surf_usm_v(l)%asdif  = surf_usm_v(l)%albedo
                   surf_usm_v(l)%asdir  = surf_usm_v(l)%albedo
                ELSE
                   surf_usm_v(l)%aldif  = albedo_lw_dif
                   surf_usm_v(l)%aldir  = albedo_lw_dir
                   surf_usm_v(l)%asdif  = albedo_sw_dif
                   surf_usm_v(l)%asdir  = albedo_sw_dir
                ENDIF
             ENDIF
          ENDDO

!
!--       Level 2 initialization of spectral albedos via albedo_type.
!--       Please note, for natural- and urban-type surfaces, a tile approach
!--       is applied so that the resulting albedo is calculated via the weighted
!--       average of respective surface fractions.
          DO  m = 1, surf_lsm_h%ns
!
!--          Spectral albedos for vegetation/pavement/water surfaces
             DO  ind_type = 0, 2
                IF ( surf_lsm_h%albedo_type(m,ind_type) /= 0 )  THEN
                   surf_lsm_h%aldif(m,ind_type) =                              &
                               albedo_pars(1,surf_lsm_h%albedo_type(m,ind_type))
                   surf_lsm_h%asdif(m,ind_type) =                              &
                               albedo_pars(2,surf_lsm_h%albedo_type(m,ind_type))
                   surf_lsm_h%aldir(m,ind_type) =                              &
                               albedo_pars(1,surf_lsm_h%albedo_type(m,ind_type))
                   surf_lsm_h%asdir(m,ind_type) =                              &
                               albedo_pars(2,surf_lsm_h%albedo_type(m,ind_type))
                   surf_lsm_h%albedo(m,ind_type) =                             &
                               albedo_pars(0,surf_lsm_h%albedo_type(m,ind_type))
                ENDIF
             ENDDO

          ENDDO
!
!--       For urban surface only if albedo has not been already initialized
!--       in the urban-surface model via the ASCII file.
          IF ( .NOT. surf_usm_h%albedo_from_ascii )  THEN
             DO  m = 1, surf_usm_h%ns
!
!--             Spectral albedos for wall/green/window surfaces
                DO  ind_type = 0, 2
                   IF ( surf_usm_h%albedo_type(m,ind_type) /= 0 )  THEN
                      surf_usm_h%aldif(m,ind_type) =                           &
                               albedo_pars(1,surf_usm_h%albedo_type(m,ind_type))
                      surf_usm_h%asdif(m,ind_type) =                           &
                               albedo_pars(2,surf_usm_h%albedo_type(m,ind_type))
                      surf_usm_h%aldir(m,ind_type) =                           &
                               albedo_pars(1,surf_usm_h%albedo_type(m,ind_type))
                      surf_usm_h%asdir(m,ind_type) =                           &
                               albedo_pars(2,surf_usm_h%albedo_type(m,ind_type))
                      surf_usm_h%albedo(m,ind_type) =                          &
                               albedo_pars(0,surf_usm_h%albedo_type(m,ind_type))
                   ENDIF
                ENDDO

             ENDDO
          ENDIF

          DO l = 0, 3

             DO  m = 1, surf_lsm_v(l)%ns
!
!--             Spectral albedos for vegetation/pavement/water surfaces
                DO  ind_type = 0, 2
                   IF ( surf_lsm_v(l)%albedo_type(m,ind_type) /= 0 )  THEN
                      surf_lsm_v(l)%aldif(m,ind_type) =                        &
                            albedo_pars(1,surf_lsm_v(l)%albedo_type(m,ind_type))
                      surf_lsm_v(l)%asdif(m,ind_type) =                        &
                            albedo_pars(2,surf_lsm_v(l)%albedo_type(m,ind_type))
                      surf_lsm_v(l)%aldir(m,ind_type) =                        &
                            albedo_pars(1,surf_lsm_v(l)%albedo_type(m,ind_type))
                      surf_lsm_v(l)%asdir(m,ind_type) =                        &
                            albedo_pars(2,surf_lsm_v(l)%albedo_type(m,ind_type))
                      surf_lsm_v(l)%albedo(m,ind_type) =                       &
                            albedo_pars(0,surf_lsm_v(l)%albedo_type(m,ind_type))
                   ENDIF
                ENDDO
             ENDDO
!
!--          For urban surface only if albedo has not been already initialized
!--          in the urban-surface model via the ASCII file.
             IF ( .NOT. surf_usm_v(l)%albedo_from_ascii )  THEN
                DO  m = 1, surf_usm_v(l)%ns
!
!--                Spectral albedos for wall/green/window surfaces
                   DO  ind_type = 0, 2
                      IF ( surf_usm_v(l)%albedo_type(m,ind_type) /= 0 )  THEN
                         surf_usm_v(l)%aldif(m,ind_type) =                     &
                            albedo_pars(1,surf_usm_v(l)%albedo_type(m,ind_type))
                         surf_usm_v(l)%asdif(m,ind_type) =                     &
                            albedo_pars(2,surf_usm_v(l)%albedo_type(m,ind_type))
                         surf_usm_v(l)%aldir(m,ind_type) =                     &
                            albedo_pars(1,surf_usm_v(l)%albedo_type(m,ind_type))
                         surf_usm_v(l)%asdir(m,ind_type) =                     &
                            albedo_pars(2,surf_usm_v(l)%albedo_type(m,ind_type))
                         surf_usm_v(l)%albedo(m,ind_type) =                    &
                            albedo_pars(0,surf_usm_v(l)%albedo_type(m,ind_type))
                      ENDIF
                   ENDDO

                ENDDO
             ENDIF
          ENDDO
!
!--       Level 3 initialization at grid points where albedo type is zero.
!--       This case, spectral albedos are taken from file if available
          IF ( albedo_pars_f%from_file )  THEN
!
!--          Horizontal
             DO  m = 1, surf_lsm_h%ns
                i = surf_lsm_h%i(m)
                j = surf_lsm_h%j(m)
!
!--             Spectral albedos for vegetation/pavement/water surfaces
                DO  ind_type = 0, 2
                   IF ( albedo_pars_f%pars_xy(0,j,i) /= albedo_pars_f%fill )   &
                      surf_lsm_h%albedo(m,ind_type) =                          &
                                             albedo_pars_f%pars_xy(0,j,i)
                   IF ( albedo_pars_f%pars_xy(1,j,i) /= albedo_pars_f%fill )   &
                      surf_lsm_h%aldir(m,ind_type) =                           &
                                             albedo_pars_f%pars_xy(1,j,i)
                   IF ( albedo_pars_f%pars_xy(1,j,i) /= albedo_pars_f%fill )   &
                      surf_lsm_h%aldif(m,ind_type) =                           &
                                             albedo_pars_f%pars_xy(1,j,i)
                   IF ( albedo_pars_f%pars_xy(2,j,i) /= albedo_pars_f%fill )   &
                      surf_lsm_h%asdir(m,ind_type) =                           &
                                             albedo_pars_f%pars_xy(2,j,i)
                   IF ( albedo_pars_f%pars_xy(2,j,i) /= albedo_pars_f%fill )   &
                      surf_lsm_h%asdif(m,ind_type) =                           &
                                             albedo_pars_f%pars_xy(2,j,i)
                ENDDO
             ENDDO
!
!--          For urban surface only if albedo has not been already initialized
!--          in the urban-surface model via the ASCII file.
             IF ( .NOT. surf_usm_h%albedo_from_ascii )  THEN
                DO  m = 1, surf_usm_h%ns
                   i = surf_usm_h%i(m)
                   j = surf_usm_h%j(m)
!
!--                Broadband albedos for wall/green/window surfaces
                   DO  ind_type = 0, 2
                      IF ( albedo_pars_f%pars_xy(0,j,i) /= albedo_pars_f%fill )&
                         surf_usm_h%albedo(m,ind_type) =                       &
                                             albedo_pars_f%pars_xy(0,j,i)
                   ENDDO
!
!--                Spectral albedos especially for building wall surfaces
                   IF ( albedo_pars_f%pars_xy(1,j,i) /= albedo_pars_f%fill )  THEN
                      surf_usm_h%aldir(m,ind_veg_wall) =                       &
                                                albedo_pars_f%pars_xy(1,j,i)
                      surf_usm_h%aldif(m,ind_veg_wall) =                       &
                                                albedo_pars_f%pars_xy(1,j,i)
                   ENDIF
                   IF ( albedo_pars_f%pars_xy(2,j,i) /= albedo_pars_f%fill )  THEN
                      surf_usm_h%asdir(m,ind_veg_wall) =                       &
                                                albedo_pars_f%pars_xy(2,j,i)
                      surf_usm_h%asdif(m,ind_veg_wall) =                       &
                                                albedo_pars_f%pars_xy(2,j,i)
                   ENDIF
!
!--                Spectral albedos especially for building green surfaces
                   IF ( albedo_pars_f%pars_xy(3,j,i) /= albedo_pars_f%fill )  THEN
                      surf_usm_h%aldir(m,ind_pav_green) =                      &
                                                albedo_pars_f%pars_xy(3,j,i)
                      surf_usm_h%aldif(m,ind_pav_green) =                      &
                                                albedo_pars_f%pars_xy(3,j,i)
                   ENDIF
                   IF ( albedo_pars_f%pars_xy(4,j,i) /= albedo_pars_f%fill )  THEN
                      surf_usm_h%asdir(m,ind_pav_green) =                      &
                                                albedo_pars_f%pars_xy(4,j,i)
                      surf_usm_h%asdif(m,ind_pav_green) =                      &
                                                albedo_pars_f%pars_xy(4,j,i)
                   ENDIF
!
!--                Spectral albedos especially for building window surfaces
                   IF ( albedo_pars_f%pars_xy(5,j,i) /= albedo_pars_f%fill )  THEN
                      surf_usm_h%aldir(m,ind_wat_win) =                        &
                                                albedo_pars_f%pars_xy(5,j,i)
                      surf_usm_h%aldif(m,ind_wat_win) =                        &
                                                albedo_pars_f%pars_xy(5,j,i)
                   ENDIF
                   IF ( albedo_pars_f%pars_xy(6,j,i) /= albedo_pars_f%fill )  THEN
                      surf_usm_h%asdir(m,ind_wat_win) =                        &
                                                albedo_pars_f%pars_xy(6,j,i)
                      surf_usm_h%asdif(m,ind_wat_win) =                        &
                                                albedo_pars_f%pars_xy(6,j,i)
                   ENDIF

                ENDDO
             ENDIF
!
!--          Vertical
             DO  l = 0, 3
                ioff = surf_lsm_v(l)%ioff
                joff = surf_lsm_v(l)%joff

                DO  m = 1, surf_lsm_v(l)%ns
                   i = surf_lsm_v(l)%i(m)
                   j = surf_lsm_v(l)%j(m)
!
!--                Spectral albedos for vegetation/pavement/water surfaces
                   DO  ind_type = 0, 2
                      IF ( albedo_pars_f%pars_xy(0,j+joff,i+ioff) /=           &
                           albedo_pars_f%fill )                                &
                         surf_lsm_v(l)%albedo(m,ind_type) =                    &
                                       albedo_pars_f%pars_xy(0,j+joff,i+ioff)
                      IF ( albedo_pars_f%pars_xy(1,j+joff,i+ioff) /=           &
                           albedo_pars_f%fill )                                &
                         surf_lsm_v(l)%aldir(m,ind_type) =                     &
                                       albedo_pars_f%pars_xy(1,j+joff,i+ioff)
                      IF ( albedo_pars_f%pars_xy(1,j+joff,i+ioff) /=           &
                           albedo_pars_f%fill )                                &
                         surf_lsm_v(l)%aldif(m,ind_type) =                     &
                                       albedo_pars_f%pars_xy(1,j+joff,i+ioff)
                      IF ( albedo_pars_f%pars_xy(2,j+joff,i+ioff) /=           &
                           albedo_pars_f%fill )                                &
                         surf_lsm_v(l)%asdir(m,ind_type) =                     &
                                       albedo_pars_f%pars_xy(2,j+joff,i+ioff)
                      IF ( albedo_pars_f%pars_xy(2,j+joff,i+ioff) /=           &
                           albedo_pars_f%fill )                                &
                         surf_lsm_v(l)%asdif(m,ind_type) =                     &
                                       albedo_pars_f%pars_xy(2,j+joff,i+ioff)
                   ENDDO
                ENDDO
!
!--             For urban surface only if albedo has not been already initialized
!--             in the urban-surface model via the ASCII file.
                IF ( .NOT. surf_usm_v(l)%albedo_from_ascii )  THEN
                   ioff = surf_usm_v(l)%ioff
                   joff = surf_usm_v(l)%joff

                   DO  m = 1, surf_usm_v(l)%ns
                      i = surf_usm_v(l)%i(m)
                      j = surf_usm_v(l)%j(m)
!
!--                   Broadband albedos for wall/green/window surfaces
                      DO  ind_type = 0, 2
                         IF ( albedo_pars_f%pars_xy(0,j+joff,i+ioff) /=        &
                              albedo_pars_f%fill )                             &
                            surf_usm_v(l)%albedo(m,ind_type) =                 &
                                          albedo_pars_f%pars_xy(0,j+joff,i+ioff)
                      ENDDO
!
!--                   Spectral albedos especially for building wall surfaces
                      IF ( albedo_pars_f%pars_xy(1,j+joff,i+ioff) /=           &
                           albedo_pars_f%fill )  THEN
                         surf_usm_v(l)%aldir(m,ind_veg_wall) =                 &
                                         albedo_pars_f%pars_xy(1,j+joff,i+ioff)
                         surf_usm_v(l)%aldif(m,ind_veg_wall) =                 &
                                         albedo_pars_f%pars_xy(1,j+joff,i+ioff)
                      ENDIF
                      IF ( albedo_pars_f%pars_xy(2,j+joff,i+ioff) /=           &
                           albedo_pars_f%fill )  THEN
                         surf_usm_v(l)%asdir(m,ind_veg_wall) =                 &
                                         albedo_pars_f%pars_xy(2,j+joff,i+ioff)
                         surf_usm_v(l)%asdif(m,ind_veg_wall) =                 &
                                         albedo_pars_f%pars_xy(2,j+joff,i+ioff)
                      ENDIF
!
!--                   Spectral albedos especially for building green surfaces
                      IF ( albedo_pars_f%pars_xy(3,j+joff,i+ioff) /=           &
                           albedo_pars_f%fill )  THEN
                         surf_usm_v(l)%aldir(m,ind_pav_green) =                &
                                         albedo_pars_f%pars_xy(3,j+joff,i+ioff)
                         surf_usm_v(l)%aldif(m,ind_pav_green) =                &
                                         albedo_pars_f%pars_xy(3,j+joff,i+ioff)
                      ENDIF
                      IF ( albedo_pars_f%pars_xy(4,j+joff,i+ioff) /=           &
                           albedo_pars_f%fill )  THEN
                         surf_usm_v(l)%asdir(m,ind_pav_green) =                &
                                         albedo_pars_f%pars_xy(4,j+joff,i+ioff)
                         surf_usm_v(l)%asdif(m,ind_pav_green) =                &
                                         albedo_pars_f%pars_xy(4,j+joff,i+ioff)
                      ENDIF
!
!--                   Spectral albedos especially for building window surfaces
                      IF ( albedo_pars_f%pars_xy(5,j+joff,i+ioff) /=           &
                           albedo_pars_f%fill )  THEN
                         surf_usm_v(l)%aldir(m,ind_wat_win) =                  &
                                         albedo_pars_f%pars_xy(5,j+joff,i+ioff)
                         surf_usm_v(l)%aldif(m,ind_wat_win) =                  &
                                         albedo_pars_f%pars_xy(5,j+joff,i+ioff)
                      ENDIF
                      IF ( albedo_pars_f%pars_xy(6,j+joff,i+ioff) /=           &
                           albedo_pars_f%fill )  THEN
                         surf_usm_v(l)%asdir(m,ind_wat_win) =                  &
                                         albedo_pars_f%pars_xy(6,j+joff,i+ioff)
                         surf_usm_v(l)%asdif(m,ind_wat_win) =                  &
                                         albedo_pars_f%pars_xy(6,j+joff,i+ioff)
                      ENDIF
                   ENDDO
                ENDIF
             ENDDO

          ENDIF
!
!--       Read explicit albedo values from building surface pars. If present,
!--       they override all less specific albedo values and force a albedo_type
!--       to zero in order to take effect.
          IF ( building_surface_pars_f%from_file )  THEN
             DO  m = 1, surf_usm_h%ns
                i = surf_usm_h%i(m)
                j = surf_usm_h%j(m)
                k = surf_usm_h%k(m)
!
!--             Iterate over surfaces in column, check height and orientation
                DO  is = building_surface_pars_f%index_ji(1,j,i), &
                         building_surface_pars_f%index_ji(2,j,i)
                   IF ( building_surface_pars_f%coords(4,is) == -surf_usm_h%koff .AND. &
                        building_surface_pars_f%coords(1,is) == k )  THEN

                      IF ( building_surface_pars_f%pars(ind_s_alb_b_wall,is) /=      &
                           building_surface_pars_f%fill )  THEN
                         surf_usm_h%albedo(m,ind_veg_wall) =                         &
                                  building_surface_pars_f%pars(ind_s_alb_b_wall,is)
                         surf_usm_h%albedo_type(m,ind_veg_wall) = 0
                      ENDIF

                      IF ( building_surface_pars_f%pars(ind_s_alb_l_wall,is) /=      &
                           building_surface_pars_f%fill )  THEN
                         surf_usm_h%aldir(m,ind_veg_wall) =                          &
                                  building_surface_pars_f%pars(ind_s_alb_l_wall,is)
                         surf_usm_h%aldif(m,ind_veg_wall) =                          &
                                  building_surface_pars_f%pars(ind_s_alb_l_wall,is)
                         surf_usm_h%albedo_type(m,ind_veg_wall) = 0
                      ENDIF

                      IF ( building_surface_pars_f%pars(ind_s_alb_s_wall,is) /=      &
                           building_surface_pars_f%fill )  THEN
                         surf_usm_h%asdir(m,ind_veg_wall) =                          &
                                  building_surface_pars_f%pars(ind_s_alb_s_wall,is)
                         surf_usm_h%asdif(m,ind_veg_wall) =                          &
                                  building_surface_pars_f%pars(ind_s_alb_s_wall,is)
                         surf_usm_h%albedo_type(m,ind_veg_wall) = 0
                      ENDIF

                      IF ( building_surface_pars_f%pars(ind_s_alb_b_win,is) /=       &
                           building_surface_pars_f%fill )  THEN
                         surf_usm_h%albedo(m,ind_wat_win) =                          &
                                  building_surface_pars_f%pars(ind_s_alb_b_win,is)
                         surf_usm_h%albedo_type(m,ind_wat_win) = 0
                      ENDIF

                      IF ( building_surface_pars_f%pars(ind_s_alb_l_win,is) /=       &
                           building_surface_pars_f%fill )  THEN
                         surf_usm_h%aldir(m,ind_wat_win) =                           &
                                  building_surface_pars_f%pars(ind_s_alb_l_win,is)
                         surf_usm_h%aldif(m,ind_wat_win) =                           &
                                  building_surface_pars_f%pars(ind_s_alb_l_win,is)
                         surf_usm_h%albedo_type(m,ind_wat_win) = 0
                      ENDIF

                      IF ( building_surface_pars_f%pars(ind_s_alb_s_win,is) /=       &
                           building_surface_pars_f%fill )  THEN
                         surf_usm_h%asdir(m,ind_wat_win) =                           &
                                  building_surface_pars_f%pars(ind_s_alb_s_win,is)
                         surf_usm_h%asdif(m,ind_wat_win) =                           &
                                  building_surface_pars_f%pars(ind_s_alb_s_win,is)
                         surf_usm_h%albedo_type(m,ind_wat_win) = 0
                      ENDIF

                      IF ( building_surface_pars_f%pars(ind_s_alb_b_green,is) /=     &
                           building_surface_pars_f%fill )  THEN
                         surf_usm_h%albedo(m,ind_pav_green) =                        &
                                  building_surface_pars_f%pars(ind_s_alb_b_green,is)
                         surf_usm_h%albedo_type(m,ind_pav_green) = 0
                      ENDIF

                      IF ( building_surface_pars_f%pars(ind_s_alb_l_green,is) /=     &
                           building_surface_pars_f%fill )  THEN
                         surf_usm_h%aldir(m,ind_pav_green) =                         &
                                  building_surface_pars_f%pars(ind_s_alb_l_green,is)
                         surf_usm_h%aldif(m,ind_pav_green) =                         &
                                  building_surface_pars_f%pars(ind_s_alb_l_green,is)
                         surf_usm_h%albedo_type(m,ind_pav_green) = 0
                      ENDIF

                      IF ( building_surface_pars_f%pars(ind_s_alb_s_green,is) /=     &
                           building_surface_pars_f%fill )  THEN
                         surf_usm_h%asdir(m,ind_pav_green) =                         &
                                  building_surface_pars_f%pars(ind_s_alb_s_green,is)
                         surf_usm_h%asdif(m,ind_pav_green) =                         &
                                  building_surface_pars_f%pars(ind_s_alb_s_green,is)
                         surf_usm_h%albedo_type(m,ind_pav_green) = 0
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
!--                Iterate over surfaces in column, check height and orientation
                   DO  is = building_surface_pars_f%index_ji(1,j,i), &
                            building_surface_pars_f%index_ji(2,j,i)
                      IF ( building_surface_pars_f%coords(5,is) == -surf_usm_v(l)%joff .AND. &
                           building_surface_pars_f%coords(6,is) == -surf_usm_v(l)%ioff .AND. &
                           building_surface_pars_f%coords(1,is) == k )  THEN

                         IF ( building_surface_pars_f%pars(ind_s_alb_b_wall,is) /=      &
                              building_surface_pars_f%fill )  THEN
                            surf_usm_v(l)%albedo(m,ind_veg_wall) =                      &
                                     building_surface_pars_f%pars(ind_s_alb_b_wall,is)
                            surf_usm_v(l)%albedo_type(m,ind_veg_wall) = 0
                         ENDIF

                         IF ( building_surface_pars_f%pars(ind_s_alb_l_wall,is) /=      &
                              building_surface_pars_f%fill )  THEN
                            surf_usm_v(l)%aldir(m,ind_veg_wall) =                       &
                                     building_surface_pars_f%pars(ind_s_alb_l_wall,is)
                            surf_usm_v(l)%aldif(m,ind_veg_wall) =                       &
                                     building_surface_pars_f%pars(ind_s_alb_l_wall,is)
                            surf_usm_v(l)%albedo_type(m,ind_veg_wall) = 0
                         ENDIF

                         IF ( building_surface_pars_f%pars(ind_s_alb_s_wall,is) /=      &
                              building_surface_pars_f%fill )  THEN
                            surf_usm_v(l)%asdir(m,ind_veg_wall) =                       &
                                     building_surface_pars_f%pars(ind_s_alb_s_wall,is)
                            surf_usm_v(l)%asdif(m,ind_veg_wall) =                       &
                                     building_surface_pars_f%pars(ind_s_alb_s_wall,is)
                            surf_usm_v(l)%albedo_type(m,ind_veg_wall) = 0
                         ENDIF

                         IF ( building_surface_pars_f%pars(ind_s_alb_b_win,is) /=       &
                              building_surface_pars_f%fill )  THEN
                            surf_usm_v(l)%albedo(m,ind_wat_win) =                       &
                                     building_surface_pars_f%pars(ind_s_alb_b_win,is)
                            surf_usm_v(l)%albedo_type(m,ind_wat_win) = 0
                         ENDIF

                         IF ( building_surface_pars_f%pars(ind_s_alb_l_win,is) /=       &
                              building_surface_pars_f%fill )  THEN
                            surf_usm_v(l)%aldir(m,ind_wat_win) =                        &
                                     building_surface_pars_f%pars(ind_s_alb_l_win,is)
                            surf_usm_v(l)%aldif(m,ind_wat_win) =                        &
                                     building_surface_pars_f%pars(ind_s_alb_l_win,is)
                            surf_usm_v(l)%albedo_type(m,ind_wat_win) = 0
                         ENDIF

                         IF ( building_surface_pars_f%pars(ind_s_alb_s_win,is) /=       &
                              building_surface_pars_f%fill )  THEN
                            surf_usm_v(l)%asdir(m,ind_wat_win) =                        &
                                     building_surface_pars_f%pars(ind_s_alb_s_win,is)
                            surf_usm_v(l)%asdif(m,ind_wat_win) =                        &
                                     building_surface_pars_f%pars(ind_s_alb_s_win,is)
                            surf_usm_v(l)%albedo_type(m,ind_wat_win) = 0
                         ENDIF

                         IF ( building_surface_pars_f%pars(ind_s_alb_b_green,is) /=     &
                              building_surface_pars_f%fill )  THEN
                            surf_usm_v(l)%albedo(m,ind_pav_green) =                     &
                                     building_surface_pars_f%pars(ind_s_alb_b_green,is)
                            surf_usm_v(l)%albedo_type(m,ind_pav_green) = 0
                         ENDIF

                         IF ( building_surface_pars_f%pars(ind_s_alb_l_green,is) /=     &
                              building_surface_pars_f%fill )  THEN
                            surf_usm_v(l)%aldir(m,ind_pav_green) =                      &
                                     building_surface_pars_f%pars(ind_s_alb_l_green,is)
                            surf_usm_v(l)%aldif(m,ind_pav_green) =                      &
                                     building_surface_pars_f%pars(ind_s_alb_l_green,is)
                            surf_usm_v(l)%albedo_type(m,ind_pav_green) = 0
                         ENDIF

                         IF ( building_surface_pars_f%pars(ind_s_alb_s_green,is) /=     &
                              building_surface_pars_f%fill )  THEN
                            surf_usm_v(l)%asdir(m,ind_pav_green) =                      &
                                     building_surface_pars_f%pars(ind_s_alb_s_green,is)
                            surf_usm_v(l)%asdif(m,ind_pav_green) =                      &
                                     building_surface_pars_f%pars(ind_s_alb_s_green,is)
                            surf_usm_v(l)%albedo_type(m,ind_pav_green) = 0
                         ENDIF

                         EXIT ! surface was found and processed
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

!
!--       Calculate initial values of current (cosine of) the zenith angle and
!--       whether the sun is up
          CALL get_date_time( time_since_reference_point, &
                              day_of_year=day_of_year,    &
                              second_of_day=second_of_day )
          CALL calc_zenith( day_of_year, second_of_day )
!
!--       Calculate initial surface albedo for different surfaces
          IF ( .NOT. constant_albedo )  THEN
#if defined( __netcdf )
!
!--          Horizontally aligned natural and urban surfaces
             CALL calc_albedo( surf_lsm_h )
             CALL calc_albedo( surf_usm_h )
!
!--          Vertically aligned natural and urban surfaces
             DO  l = 0, 3
                CALL calc_albedo( surf_lsm_v(l) )
                CALL calc_albedo( surf_usm_v(l) )
             ENDDO
#endif
          ELSE
!
!--          Initialize sun-inclination independent spectral albedos
!--          Horizontal surfaces
             IF ( surf_lsm_h%ns > 0 )  THEN
                surf_lsm_h%rrtm_aldir = surf_lsm_h%aldir
                surf_lsm_h%rrtm_asdir = surf_lsm_h%asdir
                surf_lsm_h%rrtm_aldif = surf_lsm_h%aldif
                surf_lsm_h%rrtm_asdif = surf_lsm_h%asdif
             ENDIF
             IF ( surf_usm_h%ns > 0 )  THEN
                surf_usm_h%rrtm_aldir = surf_usm_h%aldir
                surf_usm_h%rrtm_asdir = surf_usm_h%asdir
                surf_usm_h%rrtm_aldif = surf_usm_h%aldif
                surf_usm_h%rrtm_asdif = surf_usm_h%asdif
             ENDIF
!
!--          Vertical surfaces
             DO  l = 0, 3
                IF ( surf_lsm_v(l)%ns > 0 )  THEN
                   surf_lsm_v(l)%rrtm_aldir = surf_lsm_v(l)%aldir
                   surf_lsm_v(l)%rrtm_asdir = surf_lsm_v(l)%asdir
                   surf_lsm_v(l)%rrtm_aldif = surf_lsm_v(l)%aldif
                   surf_lsm_v(l)%rrtm_asdif = surf_lsm_v(l)%asdif
                ENDIF
                IF ( surf_usm_v(l)%ns > 0 )  THEN
                   surf_usm_v(l)%rrtm_aldir = surf_usm_v(l)%aldir
                   surf_usm_v(l)%rrtm_asdir = surf_usm_v(l)%asdir
                   surf_usm_v(l)%rrtm_aldif = surf_usm_v(l)%aldif
                   surf_usm_v(l)%rrtm_asdif = surf_usm_v(l)%asdif
                ENDIF
             ENDDO

          ENDIF

!
!--       Allocate 3d arrays of radiative fluxes and heating rates
          IF ( .NOT. ALLOCATED ( rad_sw_in ) )  THEN
             ALLOCATE ( rad_sw_in(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             rad_sw_in = 0.0_wp
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_sw_in_av ) )  THEN
             ALLOCATE ( rad_sw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_sw_out ) )  THEN
             ALLOCATE ( rad_sw_out(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             rad_sw_out = 0.0_wp
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_sw_out_av ) )  THEN
             ALLOCATE ( rad_sw_out_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_sw_hr ) )  THEN
             ALLOCATE ( rad_sw_hr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             rad_sw_hr = 0.0_wp
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_sw_hr_av ) )  THEN
             ALLOCATE ( rad_sw_hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             rad_sw_hr_av = 0.0_wp
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_sw_cs_hr ) )  THEN
             ALLOCATE ( rad_sw_cs_hr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             rad_sw_cs_hr = 0.0_wp
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_sw_cs_hr_av ) )  THEN
             ALLOCATE ( rad_sw_cs_hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             rad_sw_cs_hr_av = 0.0_wp
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_lw_in ) )  THEN
             ALLOCATE ( rad_lw_in(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             rad_lw_in = 0.0_wp
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_lw_in_av ) )  THEN
             ALLOCATE ( rad_lw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_lw_out ) )  THEN
             ALLOCATE ( rad_lw_out(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
            rad_lw_out = 0.0_wp
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_lw_out_av ) )  THEN
             ALLOCATE ( rad_lw_out_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_lw_hr ) )  THEN
             ALLOCATE ( rad_lw_hr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             rad_lw_hr = 0.0_wp
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_lw_hr_av ) )  THEN
             ALLOCATE ( rad_lw_hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             rad_lw_hr_av = 0.0_wp
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_lw_cs_hr ) )  THEN
             ALLOCATE ( rad_lw_cs_hr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             rad_lw_cs_hr = 0.0_wp
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_lw_cs_hr_av ) )  THEN
             ALLOCATE ( rad_lw_cs_hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             rad_lw_cs_hr_av = 0.0_wp
          ENDIF

          ALLOCATE ( rad_sw_cs_in(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ALLOCATE ( rad_sw_cs_out(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          rad_sw_cs_in  = 0.0_wp
          rad_sw_cs_out = 0.0_wp

          ALLOCATE ( rad_lw_cs_in(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ALLOCATE ( rad_lw_cs_out(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          rad_lw_cs_in  = 0.0_wp
          rad_lw_cs_out = 0.0_wp

!
!--       Allocate 1-element array for surface temperature
!--       (RRTMG anticipates an array as passed argument).
          ALLOCATE ( rrtm_tsfc(1) )
!
!--       Allocate surface emissivity.
!--       Values will be given directly before calling rrtm_lw.
          ALLOCATE ( rrtm_emis(0:0,1:nbndlw+1) )

!
!--       Initialize RRTMG, before check if files are existent
          INQUIRE( FILE='rrtmg_lw.nc', EXIST=lw_exists )
          IF ( .NOT. lw_exists )  THEN
             message_string = 'Input file rrtmg_lw.nc' //                &
                            '&for rrtmg missing. ' // &
                            '&Please provide <jobname>_lsw file in the INPUT directory.'
             CALL message( 'radiation_init', 'PA0583', 1, 2, 0, 6, 0 )
          ENDIF
          INQUIRE( FILE='rrtmg_sw.nc', EXIST=sw_exists )
          IF ( .NOT. sw_exists )  THEN
             message_string = 'Input file rrtmg_sw.nc' //                &
                            '&for rrtmg missing. ' // &
                            '&Please provide <jobname>_rsw file in the INPUT directory.'
             CALL message( 'radiation_init', 'PA0584', 1, 2, 0, 6, 0 )
          ENDIF

          IF ( lw_radiation )  CALL rrtmg_lw_ini ( c_p )
          IF ( sw_radiation )  CALL rrtmg_sw_ini ( c_p )

!
!--       Set input files for RRTMG
          INQUIRE(FILE="RAD_SND_DATA", EXIST=snd_exists)
          IF ( .NOT. snd_exists )  THEN
             rrtm_input_file = "rrtmg_lw.nc"
          ENDIF

!
!--       Read vertical layers for RRTMG from sounding data
!--       The routine provides nzt_rad, hyp_snd(1:nzt_rad),
!--       t_snd(nzt+2:nzt_rad), rrtm_play(1:nzt_rad), rrtm_plev(1_nzt_rad+1),
!--       rrtm_tlay(nzt+2:nzt_rad), rrtm_tlev(nzt+2:nzt_rad+1)
          CALL read_sounding_data

!
!--       Read trace gas profiles from file. This routine provides
!--       the rrtm_ arrays (1:nzt_rad+1)
          CALL read_trace_gas_data
#endif
       ENDIF
!
!--    Initializaion actions exclusively required for external
!--    radiation forcing
       IF ( radiation_scheme == 'external' )  THEN
!
!--       Open the radiation input file. Note, for child domain, a dynamic
!--       input file is often not provided. In order to do not need to
!--       duplicate the dynamic input file just for the radiation input, take
!--       it from the dynamic file for the parent if not available for the
!--       child domain(s). In this case this is possible because radiation
!--       input should be the same for each model.
          INQUIRE( FILE = TRIM( input_file_dynamic ),                          &
                   EXIST = radiation_input_root_domain  )

          IF ( .NOT. input_pids_dynamic  .AND.                                 &
               .NOT. radiation_input_root_domain )  THEN
             message_string = 'In case of external radiation forcing ' //      &
                              'a dynamic input file is required. If no ' //    &
                              'dynamic input for the child domain(s) is ' //   &
                              'provided, at least one for the root domain ' // &
                              'is needed.'
             CALL message( 'radiation_init', 'PA0315', 1, 2, 0, 6, 0 )
          ENDIF
#if defined( __netcdf )
!
!--       Open dynamic input file for child domain if available, else, open
!--       dynamic input file for the root domain.
          IF ( input_pids_dynamic )  THEN
             CALL open_read_file( TRIM( input_file_dynamic ) //                &
                                  TRIM( coupling_char ),                       &
                                  pids_id )
          ELSEIF ( radiation_input_root_domain )  THEN
             CALL open_read_file( TRIM( input_file_dynamic ),                  &
                                  pids_id )
          ENDIF

          CALL inquire_num_variables( pids_id, num_var_pids )
!
!--       Allocate memory to store variable names and read them
          ALLOCATE( vars_pids(1:num_var_pids) )
          CALL inquire_variable_names( pids_id, vars_pids )
!
!--       Input time dimension.
          IF ( check_existence( vars_pids, 'time_rad' ) )  THEN
             CALL get_dimension_length( pids_id, ntime, 'time_rad' )

             ALLOCATE( time_rad_f%var1d(0:ntime-1) )
!
!--          Read variable
             CALL get_variable( pids_id, 'time_rad', time_rad_f%var1d )

             time_rad_f%from_file = .TRUE.
          ENDIF
!
!--       Input shortwave downwelling.
          IF ( check_existence( vars_pids, 'rad_sw_in' ) )  THEN
!
!--          Get _FillValue attribute
             CALL get_attribute( pids_id, char_fill, rad_sw_in_f%fill,         &
                                 .FALSE., 'rad_sw_in' )
!
!--          Get level-of-detail
             CALL get_attribute( pids_id, char_lod, rad_sw_in_f%lod,           &
                                 .FALSE., 'rad_sw_in' )
!
!--          Level-of-detail 1 - radiation depends only on time_rad
             IF ( rad_sw_in_f%lod == 1 )  THEN
                ALLOCATE( rad_sw_in_f%var1d(0:ntime-1) )
                CALL get_variable( pids_id, 'rad_sw_in', rad_sw_in_f%var1d )
                rad_sw_in_f%from_file = .TRUE.
!
!--          Level-of-detail 2 - radiation depends on time_rad, y, x
             ELSEIF ( rad_sw_in_f%lod == 2 )  THEN
                ALLOCATE( rad_sw_in_f%var3d(0:ntime-1,nys:nyn,nxl:nxr) )

                CALL get_variable( pids_id, 'rad_sw_in', rad_sw_in_f%var3d,    &
                                   nxl, nxr, nys, nyn, 0, ntime-1 )

                rad_sw_in_f%from_file = .TRUE.
             ELSE
                message_string = '"rad_sw_in" has no valid lod attribute'
                CALL message( 'radiation_init', 'PA0646', 1, 2, 0, 6, 0 )
             ENDIF
          ENDIF
!
!--       Input longwave downwelling.
          IF ( check_existence( vars_pids, 'rad_lw_in' ) )  THEN
!
!--          Get _FillValue attribute
             CALL get_attribute( pids_id, char_fill, rad_lw_in_f%fill,         &
                                 .FALSE., 'rad_lw_in' )
!
!--          Get level-of-detail
             CALL get_attribute( pids_id, char_lod, rad_lw_in_f%lod,           &
                                 .FALSE., 'rad_lw_in' )
!
!--          Level-of-detail 1 - radiation depends only on time_rad
             IF ( rad_lw_in_f%lod == 1 )  THEN
                ALLOCATE( rad_lw_in_f%var1d(0:ntime-1) )
                CALL get_variable( pids_id, 'rad_lw_in', rad_lw_in_f%var1d )
                rad_lw_in_f%from_file = .TRUE.
!
!--          Level-of-detail 2 - radiation depends on time_rad, y, x
             ELSEIF ( rad_lw_in_f%lod == 2 )  THEN
                ALLOCATE( rad_lw_in_f%var3d(0:ntime-1,nys:nyn,nxl:nxr) )

                CALL get_variable( pids_id, 'rad_lw_in', rad_lw_in_f%var3d,    &
                                   nxl, nxr, nys, nyn, 0, ntime-1 )

                rad_lw_in_f%from_file = .TRUE.
             ELSE
                message_string = '"rad_lw_in" has no valid lod attribute'
                CALL message( 'radiation_init', 'PA0646', 1, 2, 0, 6, 0 )
             ENDIF
          ENDIF
!
!--       Input shortwave downwelling, diffuse part.
          IF ( check_existence( vars_pids, 'rad_sw_in_dif' ) )  THEN
!
!--          Read _FillValue attribute
             CALL get_attribute( pids_id, char_fill, rad_sw_in_dif_f%fill,     &
                                 .FALSE., 'rad_sw_in_dif' )
!
!--          Get level-of-detail
             CALL get_attribute( pids_id, char_lod, rad_sw_in_dif_f%lod,       &
                                 .FALSE., 'rad_sw_in_dif' )
!
!--          Level-of-detail 1 - radiation depends only on time_rad
             IF ( rad_sw_in_dif_f%lod == 1 )  THEN
                ALLOCATE( rad_sw_in_dif_f%var1d(0:ntime-1) )
                CALL get_variable( pids_id, 'rad_sw_in_dif',                   &
                                   rad_sw_in_dif_f%var1d )
                rad_sw_in_dif_f%from_file = .TRUE.
!
!--          Level-of-detail 2 - radiation depends on time_rad, y, x
             ELSEIF ( rad_sw_in_dif_f%lod == 2 )  THEN
                ALLOCATE( rad_sw_in_dif_f%var3d(0:ntime-1,nys:nyn,nxl:nxr) )

                CALL get_variable( pids_id, 'rad_sw_in_dif',                   &
                                   rad_sw_in_dif_f%var3d,                      &
                                   nxl, nxr, nys, nyn, 0, ntime-1 )

                rad_sw_in_dif_f%from_file = .TRUE.
             ELSE
                message_string = '"rad_sw_in_dif" has no valid lod attribute'
                CALL message( 'radiation_init', 'PA0646', 1, 2, 0, 6, 0 )
             ENDIF
          ENDIF
!
!--       Finally, close the input file and deallocate temporary arrays
          DEALLOCATE( vars_pids )

          CALL close_input_file( pids_id )
#endif
!
!--       Make some consistency checks.
          IF ( .NOT. rad_sw_in_f%from_file  .OR.                               &
               .NOT. rad_lw_in_f%from_file )  THEN
             message_string = 'In case of external radiation forcing ' //      &
                              'both, rad_sw_in and rad_lw_in are required.'
             CALL message( 'radiation_init', 'PA0195', 1, 2, 0, 6, 0 )
          ENDIF

          IF ( .NOT. time_rad_f%from_file )  THEN
             message_string = 'In case of external radiation forcing ' //      &
                              'dimension time_rad is required.'
             CALL message( 'radiation_init', 'PA0196', 1, 2, 0, 6, 0 )
          ENDIF

          CALL get_date_time( 0.0_wp, second_of_day=second_of_day )

          IF ( end_time - spinup_time > time_rad_f%var1d(ntime-1) )  THEN
             message_string = 'External radiation forcing does not cover ' //  &
                              'the entire simulation time.'
             CALL message( 'radiation_init', 'PA0314', 1, 2, 0, 6, 0 )
          ENDIF
!
!--       Check for fill values in radiation
          IF ( ALLOCATED( rad_sw_in_f%var1d ) )  THEN
             IF ( ANY( rad_sw_in_f%var1d == rad_sw_in_f%fill ) )  THEN
                message_string = 'External radiation array "rad_sw_in" ' //    &
                                 'must not contain any fill values.'
                CALL message( 'radiation_init', 'PA0197', 1, 2, 0, 6, 0 )
             ENDIF
          ENDIF

          IF ( ALLOCATED( rad_lw_in_f%var1d ) )  THEN
             IF ( ANY( rad_lw_in_f%var1d == rad_lw_in_f%fill ) )  THEN
                message_string = 'External radiation array "rad_lw_in" ' //    &
                                 'must not contain any fill values.'
                CALL message( 'radiation_init', 'PA0198', 1, 2, 0, 6, 0 )
             ENDIF
          ENDIF

          IF ( ALLOCATED( rad_sw_in_dif_f%var1d ) )  THEN
             IF ( ANY( rad_sw_in_dif_f%var1d == rad_sw_in_dif_f%fill ) )  THEN
                message_string = 'External radiation array "rad_sw_in_dif" ' //&
                                 'must not contain any fill values.'
                CALL message( 'radiation_init', 'PA0199', 1, 2, 0, 6, 0 )
             ENDIF
          ENDIF

          IF ( ALLOCATED( rad_sw_in_f%var3d ) )  THEN
             IF ( ANY( rad_sw_in_f%var3d == rad_sw_in_f%fill ) )  THEN
                message_string = 'External radiation array "rad_sw_in" ' //    &
                                 'must not contain any fill values.'
                CALL message( 'radiation_init', 'PA0197', 1, 2, 0, 6, 0 )
             ENDIF
          ENDIF

          IF ( ALLOCATED( rad_lw_in_f%var3d ) )  THEN
             IF ( ANY( rad_lw_in_f%var3d == rad_lw_in_f%fill ) )  THEN
                message_string = 'External radiation array "rad_lw_in" ' //    &
                                 'must not contain any fill values.'
                CALL message( 'radiation_init', 'PA0198', 1, 2, 0, 6, 0 )
             ENDIF
          ENDIF

          IF ( ALLOCATED( rad_sw_in_dif_f%var3d ) )  THEN
             IF ( ANY( rad_sw_in_dif_f%var3d == rad_sw_in_dif_f%fill ) )  THEN
                message_string = 'External radiation array "rad_sw_in_dif" ' //&
                                 'must not contain any fill values.'
                CALL message( 'radiation_init', 'PA0199', 1, 2, 0, 6, 0 )
             ENDIF
          ENDIF
!
!--       Currently, 2D external radiation input is not possible in
!--       combination with topography where average radiation is used.
          IF ( ( rad_lw_in_f%lod == 2  .OR.  rad_sw_in_f%lod == 2  .OR.      &
                 rad_sw_in_dif_f%lod == 2  )  .AND. average_radiation )  THEN
             message_string = 'External radiation with lod = 2 is currently '//&
                              'not possible with average_radiation = .T..'
                CALL message( 'radiation_init', 'PA0670', 1, 2, 0, 6, 0 )
          ENDIF
!
!--       All radiation input should have the same level of detail. The sum
!--       of lods divided by the number of available radiation arrays must be
!--       1 (if all are lod = 1) or 2 (if all are lod = 2).
          IF ( REAL( MERGE( rad_lw_in_f%lod, 0, rad_lw_in_f%from_file ) +       &
                     MERGE( rad_sw_in_f%lod, 0, rad_sw_in_f%from_file ) +       &
                     MERGE( rad_sw_in_dif_f%lod, 0, rad_sw_in_dif_f%from_file ),&
                     KIND = wp ) /                                              &
                   ( MERGE( 1.0_wp, 0.0_wp, rad_lw_in_f%from_file ) +           &
                     MERGE( 1.0_wp, 0.0_wp, rad_sw_in_f%from_file ) +           &
                     MERGE( 1.0_wp, 0.0_wp, rad_sw_in_dif_f%from_file ) )       &
                     /= 1.0_wp  .AND.                                           &
               REAL( MERGE( rad_lw_in_f%lod, 0, rad_lw_in_f%from_file ) +       &
                     MERGE( rad_sw_in_f%lod, 0, rad_sw_in_f%from_file ) +       &
                     MERGE( rad_sw_in_dif_f%lod, 0, rad_sw_in_dif_f%from_file ),&
                     KIND = wp ) /                                              &
                   ( MERGE( 1.0_wp, 0.0_wp, rad_lw_in_f%from_file ) +           &
                     MERGE( 1.0_wp, 0.0_wp, rad_sw_in_f%from_file ) +           &
                     MERGE( 1.0_wp, 0.0_wp, rad_sw_in_dif_f%from_file ) )       &
                     /= 2.0_wp )  THEN
             message_string = 'External radiation input should have the same '//&
                              'lod.'
             CALL message( 'radiation_init', 'PA0673', 1, 2, 0, 6, 0 )
          ENDIF

       ENDIF
!
!--    Perform user actions if required
       CALL user_init_radiation

!
!--    Calculate radiative fluxes at model start
       SELECT CASE ( TRIM( radiation_scheme ) )

          CASE ( 'rrtmg' )
             CALL radiation_rrtmg

          CASE ( 'clear-sky' )
             CALL radiation_clearsky

          CASE ( 'constant' )
             CALL radiation_constant

          CASE ( 'external' )
!
!--          During spinup apply clear-sky model
             IF ( time_since_reference_point < 0.0_wp )  THEN
                CALL radiation_clearsky
             ELSE
                CALL radiation_external
             ENDIF

          CASE DEFAULT

       END SELECT

!
!--    Find all discretized apparent solar positions for radiation interaction.
       IF ( radiation_interactions )  CALL radiation_presimulate_solar_pos

!
!--    If required, read or calculate and write out the SVF
       IF ( radiation_interactions .AND. read_svf)  THEN
!
!--       Read sky-view factors and further required data from file
          CALL radiation_read_svf()

       ELSEIF ( radiation_interactions .AND. .NOT. read_svf)  THEN
!
!--       calculate SFV and CSF
          CALL radiation_calc_svf()
       ENDIF

       IF ( radiation_interactions .AND. write_svf)  THEN
!
!--       Write svf, csf svfsurf and csfsurf data to file
          CALL radiation_write_svf()
       ENDIF

!
!--    Adjust radiative fluxes. In case of urban and land surfaces, also
!--    call an initial interaction.
       IF ( radiation_interactions )  THEN
          CALL radiation_interaction
       ENDIF

       IF ( debug_output )  CALL debug_message( 'radiation_init', 'end' )

       RETURN !todo: remove, I don't see what we need this for here

    END SUBROUTINE radiation_init


!------------------------------------------------------------------------------!
! Description:
! ------------
!> A simple clear sky radiation model
!------------------------------------------------------------------------------!
    SUBROUTINE radiation_external

       IMPLICIT NONE

       INTEGER(iwp) ::  l   !< running index for surface orientation
       INTEGER(iwp) ::  t   !< index of current timestep
       INTEGER(iwp) ::  tm  !< index of previous timestep

       LOGICAL      ::  horizontal !< flag indicating treatment of horinzontal surfaces

       REAL(wp) ::  fac_dt               !< interpolation factor
       REAL(wp) ::  second_of_day_init   !< second of the day at model start

       TYPE(surf_type), POINTER ::  surf !< pointer on respective surface type, used to generalize routine

!
!--    Calculate current zenith angle
       CALL get_date_time( time_since_reference_point, &
                           day_of_year=day_of_year,    &
                           second_of_day=second_of_day )
       CALL calc_zenith( day_of_year, second_of_day )
!
!--    Interpolate external radiation on current timestep
       IF ( time_since_reference_point  <= 0.0_wp )  THEN
          t      = 0
          tm     = 0
          fac_dt = 0
       ELSE
          CALL get_date_time( 0.0_wp, second_of_day=second_of_day_init )
          t = 0
          DO WHILE ( time_rad_f%var1d(t) <= time_since_reference_point )
             t = t + 1
          ENDDO

          tm = MAX( t-1, 0 )

          fac_dt = ( time_since_reference_point                                &
                   - time_rad_f%var1d(tm) + dt_3d )                            &
                 / ( time_rad_f%var1d(t)  - time_rad_f%var1d(tm) )
          fac_dt = MIN( 1.0_wp, fac_dt )
       ENDIF
!
!--    Call clear-sky calculation for each surface orientation.
!--    First, horizontal surfaces
       horizontal = .TRUE.
       surf => surf_lsm_h
       CALL radiation_external_surf
       surf => surf_usm_h
       CALL radiation_external_surf
       horizontal = .FALSE.
!
!--    Vertical surfaces
       DO  l = 0, 3
          surf => surf_lsm_v(l)
          CALL radiation_external_surf
          surf => surf_usm_v(l)
          CALL radiation_external_surf
       ENDDO

       CONTAINS

          SUBROUTINE radiation_external_surf

             USE control_parameters

             IMPLICIT NONE

             INTEGER(iwp) ::  i    !< grid index along x-dimension
             INTEGER(iwp) ::  j    !< grid index along y-dimension
             INTEGER(iwp) ::  k    !< grid index along z-dimension
             INTEGER(iwp) ::  m    !< running index for surface elements

             REAL(wp) ::  lw_in     !< downwelling longwave radiation, interpolated value
             REAL(wp) ::  sw_in     !< downwelling shortwave radiation, interpolated value
             REAL(wp) ::  sw_in_dif !< downwelling diffuse shortwave radiation, interpolated value

             IF ( surf%ns < 1 )  RETURN
!
!--          level-of-detail = 1. Note, here it must be distinguished between
!--          averaged radiation and non-averaged radiation for the upwelling
!--          fluxes.
             IF ( rad_sw_in_f%lod == 1 )  THEN

                sw_in = ( 1.0_wp - fac_dt ) * rad_sw_in_f%var1d(tm)            &
                                   + fac_dt * rad_sw_in_f%var1d(t)

                lw_in = ( 1.0_wp - fac_dt ) * rad_lw_in_f%var1d(tm)            &
                                   + fac_dt * rad_lw_in_f%var1d(t)
!
!--             Limit shortwave incoming radiation to positive values, in order
!--             to overcome possible observation errors.
                sw_in = MAX( 0.0_wp, sw_in )
                sw_in = MERGE( sw_in, 0.0_wp, sun_up )

                surf%rad_sw_in = sw_in
                surf%rad_lw_in = lw_in

                IF ( average_radiation )  THEN
                   surf%rad_sw_out = albedo_urb * surf%rad_sw_in

                   surf%rad_lw_out = emissivity_urb * sigma_sb * t_rad_urb**4  &
                                  + ( 1.0_wp - emissivity_urb ) * surf%rad_lw_in

                   surf%rad_net = surf%rad_sw_in - surf%rad_sw_out             &
                                + surf%rad_lw_in - surf%rad_lw_out

                   surf%rad_lw_out_change_0 = 4.0_wp * emissivity_urb          &
                                                     * sigma_sb                &
                                                     * t_rad_urb**3
                ELSE
                   DO  m = 1, surf%ns
                      k = surf%k(m)
                      surf%rad_sw_out(m) = ( surf%frac(m,ind_veg_wall)  *      &
                                             surf%albedo(m,ind_veg_wall)       &
                                           + surf%frac(m,ind_pav_green) *      &
                                             surf%albedo(m,ind_pav_green)      &
                                           + surf%frac(m,ind_wat_win)   *      &
                                             surf%albedo(m,ind_wat_win) )      &
                                           * surf%rad_sw_in(m)

                      surf%rad_lw_out(m) = ( surf%frac(m,ind_veg_wall)  *      &
                                             surf%emissivity(m,ind_veg_wall)   &
                                           + surf%frac(m,ind_pav_green) *      &
                                             surf%emissivity(m,ind_pav_green)  &
                                           + surf%frac(m,ind_wat_win)   *      &
                                             surf%emissivity(m,ind_wat_win)    &
                                           )                                   &
                                           * sigma_sb                          &
                                           * ( surf%pt_surface(m) * exner(k) )**4

                      surf%rad_lw_out_change_0(m) =                            &
                                         ( surf%frac(m,ind_veg_wall)  *        &
                                           surf%emissivity(m,ind_veg_wall)     &
                                         + surf%frac(m,ind_pav_green) *        &
                                           surf%emissivity(im,ind_pav_green)    &
                                         + surf%frac(m,ind_wat_win)   *        &
                                           surf%emissivity(m,ind_wat_win)      &
                                         ) * 4.0_wp * sigma_sb                 &
                                         * ( surf%pt_surface(m) * exner(k) )**3
                   ENDDO

                ENDIF
!
!--             If diffuse shortwave radiation is available, store it on
!--             the respective files.
                IF ( rad_sw_in_dif_f%from_file )  THEN
                   sw_in_dif= ( 1.0_wp - fac_dt ) * rad_sw_in_dif_f%var1d(tm)  &
                                         + fac_dt * rad_sw_in_dif_f%var1d(t)

                   IF ( ALLOCATED( rad_sw_in_diff ) )  rad_sw_in_diff = sw_in_dif
                   IF ( ALLOCATED( rad_sw_in_dir  ) )  rad_sw_in_dir  = sw_in  &
                                                                    - sw_in_dif
!
!--                Diffuse longwave radiation equals the total downwelling
!--                longwave radiation
                   IF ( ALLOCATED( rad_lw_in_diff ) )  rad_lw_in_diff = lw_in
                ENDIF
!
!--          level-of-detail = 2
             ELSE

                DO  m = 1, surf%ns
                   i = surf%i(m)
                   j = surf%j(m)
                   k = surf%k(m)

                   surf%rad_sw_in(m) = ( 1.0_wp - fac_dt )                     &
                                            * rad_sw_in_f%var3d(tm,j,i)        &
                                   + fac_dt * rad_sw_in_f%var3d(t,j,i)
!
!--                Limit shortwave incoming radiation to positive values, in
!--                order to overcome possible observation errors.
                   surf%rad_sw_in(m) = MAX( 0.0_wp, surf%rad_sw_in(m) )
                   surf%rad_sw_in(m) = MERGE( surf%rad_sw_in(m), 0.0_wp, sun_up )

                   surf%rad_lw_in(m) = ( 1.0_wp - fac_dt )                     &
                                            * rad_lw_in_f%var3d(tm,j,i)        &
                                   + fac_dt * rad_lw_in_f%var3d(t,j,i)
!
!--                Weighted average according to surface fraction.
                   surf%rad_sw_out(m) = ( surf%frac(m,ind_veg_wall)  *         &
                                          surf%albedo(m,ind_veg_wall)          &
                                        + surf%frac(m,ind_pav_green) *         &
                                          surf%albedo(m,ind_pav_green)         &
                                        + surf%frac(m,ind_wat_win)   *         &
                                          surf%albedo(m,ind_wat_win) )         &
                                        * surf%rad_sw_in(m)

                   surf%rad_lw_out(m) = ( surf%frac(m,ind_veg_wall)  *         &
                                          surf%emissivity(m,ind_veg_wall)      &
                                        + surf%frac(m,ind_pav_green) *         &
                                          surf%emissivity(m,ind_pav_green)     &
                                        + surf%frac(m,ind_wat_win)   *         &
                                          surf%emissivity(m,ind_wat_win)       &
                                        )                                      &
                                        * sigma_sb                             &
                                        * ( surf%pt_surface(m) * exner(k) )**4

                   surf%rad_lw_out_change_0(m) =                               &
                                      ( surf%frac(m,ind_veg_wall)  *           &
                                        surf%emissivity(m,ind_veg_wall)        &
                                      + surf%frac(m,ind_pav_green) *           &
                                        surf%emissivity(m,ind_pav_green)       &
                                      + surf%frac(m,ind_wat_win)   *           &
                                        surf%emissivity(m,ind_wat_win)         &
                                      ) * 4.0_wp * sigma_sb                    &
                                      * ( surf%pt_surface(m) * exner(k) )**3

                   surf%rad_net(m) = surf%rad_sw_in(m) - surf%rad_sw_out(m)    &
                                   + surf%rad_lw_in(m) - surf%rad_lw_out(m)
!
!--                If diffuse shortwave radiation is available, store it on
!--                the respective files.
                   IF ( rad_sw_in_dif_f%from_file )  THEN
                      IF ( ALLOCATED( rad_sw_in_diff ) )                       &
                         rad_sw_in_diff(j,i) = ( 1.0_wp - fac_dt )             &
                                              * rad_sw_in_dif_f%var3d(tm,j,i)  &
                                     + fac_dt * rad_sw_in_dif_f%var3d(t,j,i)
!
!--                   dir = sw_in - sw_in_dif.
                      IF ( ALLOCATED( rad_sw_in_dir  ) )                       &
                         rad_sw_in_dir(j,i)  = surf%rad_sw_in(m) -             &
                                               rad_sw_in_diff(j,i)
!
!--                   Diffuse longwave radiation equals the total downwelling
!--                   longwave radiation
                      IF ( ALLOCATED( rad_lw_in_diff ) )                       &
                         rad_lw_in_diff(j,i) = surf%rad_lw_in(m)
                   ENDIF

                ENDDO

             ENDIF
!
!--          Store radiation also on 2D arrays, which are still used for
!--          direct-diffuse splitting. Note, this is only required
!--          for horizontal surfaces, which covers all x,y position.
             IF ( horizontal )  THEN
                DO  m = 1, surf%ns
                   i = surf%i(m)
                   j = surf%j(m)

                   rad_sw_in(0,j,i)  = surf%rad_sw_in(m)
                   rad_lw_in(0,j,i)  = surf%rad_lw_in(m)
                   rad_sw_out(0,j,i) = surf%rad_sw_out(m)
                   rad_lw_out(0,j,i) = surf%rad_lw_out(m)
                ENDDO
             ENDIF

          END SUBROUTINE radiation_external_surf

    END SUBROUTINE radiation_external

!------------------------------------------------------------------------------!
! Description:
! ------------
!> A simple clear sky radiation model
!------------------------------------------------------------------------------!
    SUBROUTINE radiation_clearsky

       IMPLICIT NONE

       INTEGER(iwp) ::  l         !< running index for surface orientation

       LOGICAL      ::  horizontal !< flag indicating treatment of horinzontal surfaces

       REAL(wp)     ::  pt1       !< potential temperature at first grid level or mean value at urban layer top
       REAL(wp)     ::  pt1_l     !< potential temperature at first grid level or mean value at urban layer top at local subdomain
       REAL(wp)     ::  ql1       !< liquid water mixing ratio at first grid level or mean value at urban layer top
       REAL(wp)     ::  ql1_l     !< liquid water mixing ratio at first grid level or mean value at urban layer top at local subdomain

       TYPE(surf_type), POINTER ::  surf !< pointer on respective surface type, used to generalize routine

!
!--    Calculate current zenith angle
       CALL get_date_time( time_since_reference_point, &
                           day_of_year=day_of_year,    &
                           second_of_day=second_of_day )
       CALL calc_zenith( day_of_year, second_of_day )

!
!--    Calculate sky transmissivity
       sky_trans = 0.6_wp + 0.2_wp * cos_zenith

!
!--    Calculate value of the Exner function at model surface
!
!--    In case averaged radiation is used, calculate mean temperature and
!--    liquid water mixing ratio at the urban-layer top.
       IF ( average_radiation ) THEN
          pt1   = 0.0_wp
          IF ( bulk_cloud_model  .OR.  cloud_droplets )  ql1   = 0.0_wp

          pt1_l = SUM( pt(nz_urban_t,nys:nyn,nxl:nxr) )
          IF ( bulk_cloud_model  .OR.  cloud_droplets  )  ql1_l = SUM( ql(nz_urban_t,nys:nyn,nxl:nxr) )

#if defined( __parallel )
          IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
          CALL MPI_ALLREDUCE( pt1_l, pt1, 1, MPI_REAL, MPI_SUM, comm2d, ierr )
          IF ( ierr /= 0 ) THEN
              WRITE(9,*) 'Error MPI_AllReduce1:', ierr, pt1_l, pt1
              FLUSH(9)
          ENDIF

          IF ( bulk_cloud_model  .OR.  cloud_droplets ) THEN
              CALL MPI_ALLREDUCE( ql1_l, ql1, 1, MPI_REAL, MPI_SUM, comm2d, ierr )
              IF ( ierr /= 0 ) THEN
                  WRITE(9,*) 'Error MPI_AllReduce2:', ierr, ql1_l, ql1
                  FLUSH(9)
              ENDIF
          ENDIF
#else
          pt1 = pt1_l
          IF ( bulk_cloud_model  .OR.  cloud_droplets )  ql1 = ql1_l
#endif

          IF ( bulk_cloud_model  .OR.  cloud_droplets  )  pt1 = pt1 + lv_d_cp / exner(nz_urban_t) * ql1
!
!--       Finally, divide by number of grid points
          pt1 = pt1 / REAL( ( nx + 1 ) * ( ny + 1 ), KIND=wp )
       ENDIF
!
!--    Call clear-sky calculation for each surface orientation.
!--    First, horizontal surfaces
       horizontal = .TRUE.
       surf => surf_lsm_h
       CALL radiation_clearsky_surf
       surf => surf_usm_h
       CALL radiation_clearsky_surf
       horizontal = .FALSE.
!
!--    Vertical surfaces
       DO  l = 0, 3
          surf => surf_lsm_v(l)
          CALL radiation_clearsky_surf
          surf => surf_usm_v(l)
          CALL radiation_clearsky_surf
       ENDDO

       CONTAINS

          SUBROUTINE radiation_clearsky_surf

             IMPLICIT NONE

             INTEGER(iwp) ::  i         !< index x-direction
             INTEGER(iwp) ::  j         !< index y-direction
             INTEGER(iwp) ::  k         !< index z-direction
             INTEGER(iwp) ::  m         !< running index for surface elements

             IF ( surf%ns < 1 )  RETURN

!
!--          Calculate radiation fluxes and net radiation (rad_net) assuming
!--          homogeneous urban radiation conditions.
             IF ( average_radiation ) THEN

                k = nz_urban_t

                surf%rad_sw_in  = solar_constant * sky_trans * cos_zenith
                surf%rad_sw_out = albedo_urb * surf%rad_sw_in

                surf%rad_lw_in  = emissivity_atm_clsky * sigma_sb * (pt1 * exner(k+1))**4

                surf%rad_lw_out = emissivity_urb * sigma_sb * (t_rad_urb)**4   &
                                    + (1.0_wp - emissivity_urb) * surf%rad_lw_in

                surf%rad_net = surf%rad_sw_in - surf%rad_sw_out                &
                             + surf%rad_lw_in - surf%rad_lw_out

                surf%rad_lw_out_change_0 = 4.0_wp * emissivity_urb * sigma_sb  &
                                           * (t_rad_urb)**3

!
!--          Calculate radiation fluxes and net radiation (rad_net) for each surface
!--          element.
             ELSE

                DO  m = 1, surf%ns
                   i = surf%i(m)
                   j = surf%j(m)
                   k = surf%k(m)

                   surf%rad_sw_in(m) = solar_constant * sky_trans * cos_zenith

!
!--                Weighted average according to surface fraction.
!--                ATTENTION: when radiation interactions are switched on the
!--                calculated fluxes below are not actually used as they are
!--                overwritten in radiation_interaction.
                   surf%rad_sw_out(m) = ( surf%frac(m,ind_veg_wall)  *         &
                                          surf%albedo(m,ind_veg_wall)          &
                                        + surf%frac(m,ind_pav_green) *         &
                                          surf%albedo(m,ind_pav_green)         &
                                        + surf%frac(m,ind_wat_win)   *         &
                                          surf%albedo(m,ind_wat_win) )         &
                                        * surf%rad_sw_in(m)

                   surf%rad_lw_out(m) = ( surf%frac(m,ind_veg_wall)  *         &
                                          surf%emissivity(m,ind_veg_wall)      &
                                        + surf%frac(m,ind_pav_green) *         &
                                          surf%emissivity(m,ind_pav_green)     &
                                        + surf%frac(m,ind_wat_win)   *         &
                                          surf%emissivity(m,ind_wat_win)       &
                                        )                                      &
                                        * sigma_sb                             &
                                        * ( surf%pt_surface(m) * exner(nzb) )**4

                   surf%rad_lw_out_change_0(m) =                               &
                                      ( surf%frac(m,ind_veg_wall)  *           &
                                        surf%emissivity(m,ind_veg_wall)        &
                                      + surf%frac(m,ind_pav_green) *           &
                                        surf%emissivity(m,ind_pav_green)       &
                                      + surf%frac(m,ind_wat_win)   *           &
                                        surf%emissivity(m,ind_wat_win)         &
                                      ) * 4.0_wp * sigma_sb                    &
                                      * ( surf%pt_surface(m) * exner(nzb) )** 3


                   IF ( bulk_cloud_model  .OR.  cloud_droplets  )  THEN
                      pt1 = pt(k,j,i) + lv_d_cp / exner(k) * ql(k,j,i)
                      surf%rad_lw_in(m)  = emissivity_atm_clsky * sigma_sb * (pt1 * exner(k))**4
                   ELSE
                      surf%rad_lw_in(m)  = emissivity_atm_clsky * sigma_sb * (pt(k,j,i) * exner(k))**4
                   ENDIF

                   surf%rad_net(m) = surf%rad_sw_in(m) - surf%rad_sw_out(m)    &
                                   + surf%rad_lw_in(m) - surf%rad_lw_out(m)

                ENDDO

             ENDIF

!
!--          Fill out values in radiation arrays. Note, this is only required
!--          for horizontal surfaces, which covers all x,y position.
             IF ( horizontal )  THEN
                DO  m = 1, surf%ns
                   i = surf%i(m)
                   j = surf%j(m)
                   rad_sw_in(0,j,i)  = surf%rad_sw_in(m)
                   rad_sw_out(0,j,i) = surf%rad_sw_out(m)
                   rad_lw_in(0,j,i)  = surf%rad_lw_in(m)
                   rad_lw_out(0,j,i) = surf%rad_lw_out(m)
                ENDDO
             ENDIF

          END SUBROUTINE radiation_clearsky_surf

    END SUBROUTINE radiation_clearsky


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This scheme keeps the prescribed net radiation constant during the run
!------------------------------------------------------------------------------!
    SUBROUTINE radiation_constant


       IMPLICIT NONE

       INTEGER(iwp) ::  l         !< running index for surface orientation

       LOGICAL      ::  horizontal !< flag indicating treatment of horinzontal surfaces

       REAL(wp)     ::  pt1       !< potential temperature at first grid level or mean value at urban layer top
       REAL(wp)     ::  pt1_l     !< potential temperature at first grid level or mean value at urban layer top at local subdomain
       REAL(wp)     ::  ql1       !< liquid water mixing ratio at first grid level or mean value at urban layer top
       REAL(wp)     ::  ql1_l     !< liquid water mixing ratio at first grid level or mean value at urban layer top at local subdomain

       TYPE(surf_type), POINTER ::  surf !< pointer on respective surface type, used to generalize routine

!
!--    In case averaged radiation is used, calculate mean temperature and
!--    liquid water mixing ratio at the urban-layer top.
       IF ( average_radiation ) THEN
          pt1   = 0.0_wp
          IF ( bulk_cloud_model  .OR.  cloud_droplets )  ql1   = 0.0_wp

          pt1_l = SUM( pt(nz_urban_t,nys:nyn,nxl:nxr) )
          IF ( bulk_cloud_model  .OR.  cloud_droplets )  ql1_l = SUM( ql(nz_urban_t,nys:nyn,nxl:nxr) )

#if defined( __parallel )
          IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
          CALL MPI_ALLREDUCE( pt1_l, pt1, 1, MPI_REAL, MPI_SUM, comm2d, ierr )
          IF ( ierr /= 0 ) THEN
              WRITE(9,*) 'Error MPI_AllReduce3:', ierr, pt1_l, pt1
              FLUSH(9)
          ENDIF
          IF ( bulk_cloud_model  .OR.  cloud_droplets )  THEN
             CALL MPI_ALLREDUCE( ql1_l, ql1, 1, MPI_REAL, MPI_SUM, comm2d, ierr )
             IF ( ierr /= 0 ) THEN
                 WRITE(9,*) 'Error MPI_AllReduce4:', ierr, ql1_l, ql1
                 FLUSH(9)
             ENDIF
          ENDIF
#else
          pt1 = pt1_l
          IF ( bulk_cloud_model  .OR.  cloud_droplets )  ql1 = ql1_l
#endif
          IF ( bulk_cloud_model  .OR.  cloud_droplets )  pt1 = pt1 + lv_d_cp / exner(nz_urban_t+1) * ql1
!
!--       Finally, divide by number of grid points
          pt1 = pt1 / REAL( ( nx + 1 ) * ( ny + 1 ), KIND=wp )
       ENDIF

!
!--    First, horizontal surfaces
       horizontal = .TRUE.
       surf => surf_lsm_h
       CALL radiation_constant_surf
       surf => surf_usm_h
       CALL radiation_constant_surf
       horizontal = .FALSE.
!
!--    Vertical surfaces
       DO  l = 0, 3
          surf => surf_lsm_v(l)
          CALL radiation_constant_surf
          surf => surf_usm_v(l)
          CALL radiation_constant_surf
       ENDDO

       CONTAINS

          SUBROUTINE radiation_constant_surf

             IMPLICIT NONE

             INTEGER(iwp) ::  i         !< index x-direction
             INTEGER(iwp) ::  ioff      !< offset between surface element and adjacent grid point along x
             INTEGER(iwp) ::  j         !< index y-direction
             INTEGER(iwp) ::  joff      !< offset between surface element and adjacent grid point along y
             INTEGER(iwp) ::  k         !< index z-direction
             INTEGER(iwp) ::  koff      !< offset between surface element and adjacent grid point along z
             INTEGER(iwp) ::  m         !< running index for surface elements

             IF ( surf%ns < 1 )  RETURN

!--          Calculate homogenoeus urban radiation fluxes
             IF ( average_radiation ) THEN

                surf%rad_net = net_radiation

                surf%rad_lw_in  = emissivity_atm_clsky * sigma_sb * (pt1 * exner(nz_urban_t+1))**4

                surf%rad_lw_out = emissivity_urb * sigma_sb * (t_rad_urb)**4   &
                                    + ( 1.0_wp - emissivity_urb )             & ! shouldn't be this a bulk value -- emissivity_urb?
                                    * surf%rad_lw_in

                surf%rad_lw_out_change_0 = 4.0_wp * emissivity_urb * sigma_sb  &
                                           * t_rad_urb**3

                surf%rad_sw_in = ( surf%rad_net - surf%rad_lw_in               &
                                     + surf%rad_lw_out )                       &
                                     / ( 1.0_wp - albedo_urb )

                surf%rad_sw_out =  albedo_urb * surf%rad_sw_in

!
!--          Calculate radiation fluxes for each surface element
             ELSE
!
!--             Determine index offset between surface element and adjacent
!--             atmospheric grid point
                ioff = surf%ioff
                joff = surf%joff
                koff = surf%koff

!
!--             Prescribe net radiation and estimate the remaining radiative fluxes
                DO  m = 1, surf%ns
                   i = surf%i(m)
                   j = surf%j(m)
                   k = surf%k(m)

                   surf%rad_net(m) = net_radiation

                   IF ( bulk_cloud_model  .OR.  cloud_droplets )  THEN
                      pt1 = pt(k,j,i) + lv_d_cp / exner(k) * ql(k,j,i)
                      surf%rad_lw_in(m)  = emissivity_atm_clsky * sigma_sb * (pt1 * exner(k))**4
                   ELSE
                      surf%rad_lw_in(m)  = emissivity_atm_clsky * sigma_sb *                 &
                                             ( pt(k,j,i) * exner(k) )**4
                   ENDIF

!
!--                Weighted average according to surface fraction.
                   surf%rad_lw_out(m) = ( surf%frac(m,ind_veg_wall)  *         &
                                          surf%emissivity(m,ind_veg_wall)      &
                                        + surf%frac(m,ind_pav_green) *         &
                                          surf%emissivity(m,ind_pav_green)     &
                                        + surf%frac(m,ind_wat_win)   *         &
                                          surf%emissivity(m,ind_wat_win)       &
                                        )                                      &
                                      * sigma_sb                               &
                                      * ( surf%pt_surface(m) * exner(nzb) )**4

                   surf%rad_sw_in(m) = ( surf%rad_net(m) - surf%rad_lw_in(m)   &
                                       + surf%rad_lw_out(m) )                  &
                                       / ( 1.0_wp -                            &
                                          ( surf%frac(m,ind_veg_wall)  *       &
                                            surf%albedo(m,ind_veg_wall)        &
                                         +  surf%frac(m,ind_pav_green) *       &
                                            surf%albedo(m,ind_pav_green)       &
                                         +  surf%frac(m,ind_wat_win)   *       &
                                            surf%albedo(m,ind_wat_win) )       &
                                         )

                   surf%rad_sw_out(m) = ( surf%frac(m,ind_veg_wall)  *         &
                                          surf%albedo(m,ind_veg_wall)          &
                                        + surf%frac(m,ind_pav_green) *         &
                                          surf%albedo(m,ind_pav_green)         &
                                        + surf%frac(m,ind_wat_win)   *         &
                                          surf%albedo(m,ind_wat_win) )         &
                                      * surf%rad_sw_in(m)

                ENDDO

             ENDIF

!
!--          Fill out values in radiation arrays. Note, this is only required
!--          for horizontal surfaces, which covers all x,y position.
             IF ( horizontal )  THEN
                DO  m = 1, surf%ns
                   i = surf%i(m)
                   j = surf%j(m)
                   rad_sw_in(0,j,i)  = surf%rad_sw_in(m)
                   rad_sw_out(0,j,i) = surf%rad_sw_out(m)
                   rad_lw_in(0,j,i)  = surf%rad_lw_in(m)
                   rad_lw_out(0,j,i) = surf%rad_lw_out(m)
                ENDDO
             ENDIF

          END SUBROUTINE radiation_constant_surf


    END SUBROUTINE radiation_constant

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Header output for radiation model
!------------------------------------------------------------------------------!
    SUBROUTINE radiation_header ( io )


       IMPLICIT NONE

       INTEGER(iwp), INTENT(IN) ::  io            !< Unit of the output file



!
!--    Write radiation model header
       WRITE( io, 3 )

       IF ( radiation_scheme == "constant" )  THEN
          WRITE( io, 4 ) net_radiation
       ELSEIF ( radiation_scheme == "clear-sky" )  THEN
          WRITE( io, 5 )
       ELSEIF ( radiation_scheme == "rrtmg" )  THEN
          WRITE( io, 6 )
          IF ( .NOT. lw_radiation )  WRITE( io, 10 )
          IF ( .NOT. sw_radiation )  WRITE( io, 11 )
       ELSEIF ( radiation_scheme == "external" )  THEN
          WRITE( io, 14 )
       ENDIF

       IF ( albedo_type_f%from_file  .OR.  vegetation_type_f%from_file  .OR.   &
            pavement_type_f%from_file  .OR.  water_type_f%from_file  .OR.      &
            building_type_f%from_file )  THEN
             WRITE( io, 13 )
       ELSE
          IF ( albedo_type == 0 )  THEN
             WRITE( io, 7 ) albedo
          ELSE
             WRITE( io, 8 ) TRIM( albedo_type_name(albedo_type) )
          ENDIF
       ENDIF
       IF ( constant_albedo )  THEN
          WRITE( io, 9 )
       ENDIF

       WRITE( io, 12 ) dt_radiation


 3 FORMAT (//' Radiation model information:'/                                  &
              ' ----------------------------'/)
 4 FORMAT ('    --> Using constant net radiation: net_radiation = ', F6.2,     &
           // 'W/m**2')
 5 FORMAT ('    --> Simple radiation scheme for clear sky is used (no clouds,',&
                   ' default)')
 6 FORMAT ('    --> RRTMG scheme is used')
 7 FORMAT (/'    User-specific surface albedo: albedo =', F6.3)
 8 FORMAT (/'    Albedo is set for land surface type: ', A)
 9 FORMAT (/'    --> Albedo is fixed during the run')
10 FORMAT (/'    --> Longwave radiation is disabled')
11 FORMAT (/'    --> Shortwave radiation is disabled.')
12 FORMAT  ('    Timestep: dt_radiation = ', F6.2, '  s')
13 FORMAT (/'    Albedo is set individually for each xy-location, according ', &
                 'to given surface type.')
14 FORMAT ('    --> External radiation forcing is used')


    END SUBROUTINE radiation_header


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Parin for &radiation_parameters for radiation model and RTM
!------------------------------------------------------------------------------!
    SUBROUTINE radiation_parin


       IMPLICIT NONE

       CHARACTER (LEN=80) ::  line  !< dummy string that contains the current line of the parameter file

       NAMELIST /radiation_par/   albedo, albedo_lw_dif, albedo_lw_dir,         &
                                  albedo_sw_dif, albedo_sw_dir, albedo_type,    &
                                  constant_albedo, dt_radiation, emissivity,    &
                                  lw_radiation, max_raytracing_dist,            &
                                  min_irrf_value, mrt_geom, mrt_geom_params,    &
                                  mrt_include_sw, mrt_nlevels,                  &
                                  mrt_skip_roof, net_radiation, nrefsteps,      &
                                  plant_lw_interact, rad_angular_discretization,&
                                  radiation_interactions_on, radiation_scheme,  &
                                  raytrace_discrete_azims,                      &
                                  raytrace_discrete_elevs, raytrace_mpi_rma,    &
                                  trace_fluxes_above,                           &
                                  skip_time_do_radiation, surface_reflections,  &
                                  svfnorm_report_thresh, sw_radiation,          &
                                  unscheduled_radiation_calls


       NAMELIST /radiation_parameters/ albedo, albedo_lw_dif, albedo_lw_dir,    &
                                  albedo_sw_dif, albedo_sw_dir, albedo_type,    &
                                  constant_albedo, dt_radiation, emissivity,    &
                                  lw_radiation, max_raytracing_dist,            &
                                  min_irrf_value, mrt_geom, mrt_geom_params,    &
                                  mrt_include_sw, mrt_nlevels,                  &
                                  mrt_skip_roof, net_radiation, nrefsteps,      &
                                  plant_lw_interact, rad_angular_discretization,&
                                  radiation_interactions_on, radiation_scheme,  &
                                  raytrace_discrete_azims,                      &
                                  raytrace_discrete_elevs, raytrace_mpi_rma,    &
                                  trace_fluxes_above,                           &
                                  skip_time_do_radiation, surface_reflections,  &
                                  svfnorm_report_thresh, sw_radiation,          &
                                  unscheduled_radiation_calls

       line = ' '

!
!--    Try to find radiation model namelist
       REWIND ( 11 )
       line = ' '
       DO WHILE ( INDEX( line, '&radiation_parameters' ) == 0 )
          READ ( 11, '(A)', END=12 )  line
       ENDDO
       BACKSPACE ( 11 )

!
!--    Read user-defined namelist
       READ ( 11, radiation_parameters, ERR = 10 )

!
!--    Set flag that indicates that the radiation model is switched on
       radiation = .TRUE.

       GOTO 14

 10    BACKSPACE( 11 )
       READ( 11 , '(A)') line
       CALL parin_fail_message( 'radiation_parameters', line )
!
!--    Try to find old namelist
 12    REWIND ( 11 )
       line = ' '
       DO WHILE ( INDEX( line, '&radiation_par' ) == 0 )
          READ ( 11, '(A)', END=14 )  line
       ENDDO
       BACKSPACE ( 11 )

!
!--    Read user-defined namelist
       READ ( 11, radiation_par, ERR = 13, END = 14 )

       message_string = 'namelist radiation_par is deprecated and will be ' // &
                     'removed in near future. Please use namelist ' //         &
                     'radiation_parameters instead'
       CALL message( 'radiation_parin', 'PA0487', 0, 1, 0, 6, 0 )

!
!--    Set flag that indicates that the radiation model is switched on
       radiation = .TRUE.

       IF ( .NOT.  radiation_interactions_on  .AND.  surface_reflections )  THEN
          message_string = 'surface_reflections is allowed only when '      // &
               'radiation_interactions_on is set to TRUE'
          CALL message( 'radiation_parin', 'PA0293',1, 2, 0, 6, 0 )
       ENDIF

       GOTO 14

 13    BACKSPACE( 11 )
       READ( 11 , '(A)') line
       CALL parin_fail_message( 'radiation_par', line )

 14    CONTINUE

    END SUBROUTINE radiation_parin


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Implementation of the RRTMG radiation_scheme
!------------------------------------------------------------------------------!
    SUBROUTINE radiation_rrtmg

#if defined ( __rrtmg )
       USE exchange_horiz_mod,                                                 &
           ONLY:  exchange_horiz

       USE indices,                                                            &
           ONLY:  nbgp

       USE palm_date_time_mod,                                                 &
           ONLY:  hours_per_day

       USE particle_attributes,                                                &
           ONLY:  grid_particles, number_of_particles, particles, prt_count

       IMPLICIT NONE


       INTEGER(iwp) ::  i, j, k, l, m, n !< loop indices
       INTEGER(iwp) ::  k_topo_l   !< topography top index
       INTEGER(iwp) ::  k_topo     !< topography top index

       REAL(wp)     ::  d_hours_day  !< 1 / hours-per-day
       REAL(wp)     ::  nc_rad, &    !< number concentration of cloud droplets
                        s_r2,   &    !< weighted sum over all droplets with r^2
                        s_r3         !< weighted sum over all droplets with r^3

       REAL(wp), DIMENSION(0:nzt+1) :: pt_av, q_av, ql_av
       REAL(wp), DIMENSION(0:0)     :: zenith   !< to provide indexed array
!
!--    Just dummy arguments
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: rrtm_lw_taucld_dum,          &
                                                  rrtm_lw_tauaer_dum,          &
                                                  rrtm_sw_taucld_dum,          &
                                                  rrtm_sw_ssacld_dum,          &
                                                  rrtm_sw_asmcld_dum,          &
                                                  rrtm_sw_fsfcld_dum,          &
                                                  rrtm_sw_tauaer_dum,          &
                                                  rrtm_sw_ssaaer_dum,          &
                                                  rrtm_sw_asmaer_dum,          &
                                                  rrtm_sw_ecaer_dum

!
!--    Pre-calculate parameters
       d_hours_day = 1.0_wp / REAL( hours_per_day, KIND=wp )

!
!--    Calculate current (cosine of) zenith angle and whether the sun is up
       CALL get_date_time( time_since_reference_point, &
                           day_of_year=day_of_year,    &
                           second_of_day=second_of_day )
       CALL calc_zenith( day_of_year, second_of_day )
       zenith(0) = cos_zenith
!
!--    Calculate surface albedo. In case average radiation is applied,
!--    this is not required.
#if defined( __netcdf )
       IF ( .NOT. constant_albedo )  THEN
!
!--       Horizontally aligned default, natural and urban surfaces
          CALL calc_albedo( surf_lsm_h    )
          CALL calc_albedo( surf_usm_h    )
!
!--       Vertically aligned default, natural and urban surfaces
          DO  l = 0, 3
             CALL calc_albedo( surf_lsm_v(l) )
             CALL calc_albedo( surf_usm_v(l) )
          ENDDO
       ENDIF
#endif

!
!--    Prepare input data for RRTMG

!
!--    In case of large scale forcing with surface data, calculate new pressure
!--    profile. nzt_rad might be modified by these calls and all required arrays
!--    will then be re-allocated
       IF ( large_scale_forcing  .AND.  lsf_surf )  THEN
          CALL read_sounding_data
          CALL read_trace_gas_data
       ENDIF


       IF ( average_radiation ) THEN
!
!--       Determine minimum topography top index.
          k_topo_l = MINVAL( topo_top_ind(nys:nyn,nxl:nxr,0) )
#if defined( __parallel )
          CALL MPI_ALLREDUCE( k_topo_l, k_topo, 1, MPI_INTEGER, MPI_MIN, &
                              comm2d, ierr)
#else
          k_topo = k_topo_l
#endif

          rrtm_asdir(1)  = albedo_urb
          rrtm_asdif(1)  = albedo_urb
          rrtm_aldir(1)  = albedo_urb
          rrtm_aldif(1)  = albedo_urb

          rrtm_emis = emissivity_urb
!
!--       Calculate mean pt profile.
          CALL calc_mean_profile( pt, 4 )
          pt_av = hom(:, 1, 4, 0)

          IF ( humidity )  THEN
             CALL calc_mean_profile( q, 41 )
             q_av  = hom(:, 1, 41, 0)
          ENDIF
!
!--       Prepare profiles of temperature and H2O volume mixing ratio
          rrtm_tlev(0,k_topo+1) = t_rad_urb

          IF ( bulk_cloud_model )  THEN

             CALL calc_mean_profile( ql, 54 )
             ! average ql is now in hom(:, 1, 54, 0)
             ql_av = hom(:, 1, 54, 0)

             DO k = nzb+1, nzt+1
                rrtm_tlay(0,k) = pt_av(k) * ( (hyp(k) ) / 100000._wp       &
                                 )**.286_wp + lv_d_cp * ql_av(k)
                rrtm_h2ovmr(0,k) = mol_mass_air_d_wv * (q_av(k) - ql_av(k))
             ENDDO
          ELSE
             DO k = nzb+1, nzt+1
                rrtm_tlay(0,k) = pt_av(k) * ( (hyp(k) ) / 100000._wp       &
                                 )**.286_wp
             ENDDO

             IF ( humidity )  THEN
                DO k = nzb+1, nzt+1
                   rrtm_h2ovmr(0,k) = mol_mass_air_d_wv * q_av(k)
                ENDDO
             ELSE
                rrtm_h2ovmr(0,nzb+1:nzt+1) = 0.0_wp
             ENDIF
          ENDIF

!
!--       Avoid temperature/humidity jumps at the top of the PALM domain by
!--       linear interpolation from nzt+2 to nzt+7. Jumps are induced by
!--       discrepancies between the values in the  domain and those above that
!--       are prescribed in RRTMG
          DO k = nzt+2, nzt+7
             rrtm_tlay(0,k) = rrtm_tlay(0,nzt+1)                            &
                           + ( rrtm_tlay(0,nzt+8) - rrtm_tlay(0,nzt+1) )    &
                           / ( rrtm_play(0,nzt+8) - rrtm_play(0,nzt+1) )    &
                           * ( rrtm_play(0,k) - rrtm_play(0,nzt+1) )

             rrtm_h2ovmr(0,k) = rrtm_h2ovmr(0,nzt+1)                        &
                           + ( rrtm_h2ovmr(0,nzt+8) - rrtm_h2ovmr(0,nzt+1) )&
                           / ( rrtm_play(0,nzt+8)   - rrtm_play(0,nzt+1)   )&
                           * ( rrtm_play(0,k) - rrtm_play(0,nzt+1) )

          ENDDO

!--       Linear interpolate to zw grid. Loop reaches one level further up
!--       due to the staggered grid in RRTMG
          DO k = k_topo+2, nzt+8
             rrtm_tlev(0,k)   = rrtm_tlay(0,k-1) + (rrtm_tlay(0,k) -        &
                                rrtm_tlay(0,k-1))                           &
                                / ( rrtm_play(0,k) - rrtm_play(0,k-1) )     &
                                * ( rrtm_plev(0,k) - rrtm_play(0,k-1) )
          ENDDO
!
!--       Calculate liquid water path and cloud fraction for each column.
!--       Note that LWP is required in g/m2 instead of kg/kg m.
          rrtm_cldfr  = 0.0_wp
          rrtm_reliq  = 0.0_wp
          rrtm_cliqwp = 0.0_wp
          rrtm_icld   = 0

          IF ( bulk_cloud_model )  THEN
             DO k = nzb+1, nzt+1
                rrtm_cliqwp(0,k) =  ql_av(k) * 1000._wp *                   &
                                    (rrtm_plev(0,k) - rrtm_plev(0,k+1))     &
                                    * 100._wp / g

                IF ( rrtm_cliqwp(0,k) > 0._wp )  THEN
                   rrtm_cldfr(0,k) = 1._wp
                   IF ( rrtm_icld == 0 )  rrtm_icld = 1

!
!--                Calculate cloud droplet effective radius
                   rrtm_reliq(0,k) = 1.0E6_wp * ( 3.0_wp * ql_av(k)         &
                                     * rho_surface                          &
                                     / ( 4.0_wp * pi * nc_const * rho_l )   &
                                     )**0.33333333333333_wp                 &
                                     * EXP( LOG( sigma_gc )**2 )
!
!--                Limit effective radius
                   IF ( rrtm_reliq(0,k) > 0.0_wp )  THEN
                      rrtm_reliq(0,k) = MAX(rrtm_reliq(0,k),2.5_wp)
                      rrtm_reliq(0,k) = MIN(rrtm_reliq(0,k),60.0_wp)
                   ENDIF
                ENDIF
             ENDDO
          ENDIF

!
!--       Set surface temperature
          rrtm_tsfc = t_rad_urb

          IF ( lw_radiation )  THEN
!
!--          Due to technical reasons, copy optical depth to dummy arguments
!--          which are allocated on the exact size as the rrtmg_lw is called.
!--          As one dimesion is allocated with zero size, compiler complains
!--          that rank of the array does not match that of the
!--          assumed-shaped arguments in the RRTMG library. In order to
!--          avoid this, write to dummy arguments and give pass the entire
!--          dummy array. Seems to be the only existing work-around.
             ALLOCATE( rrtm_lw_taucld_dum(1:nbndlw+1,0:0,k_topo+1:nzt_rad+1) )
             ALLOCATE( rrtm_lw_tauaer_dum(0:0,k_topo+1:nzt_rad+1,1:nbndlw+1) )

             rrtm_lw_taucld_dum =                                              &
                             rrtm_lw_taucld(1:nbndlw+1,0:0,k_topo+1:nzt_rad+1)
             rrtm_lw_tauaer_dum =                                              &
                             rrtm_lw_tauaer(0:0,k_topo+1:nzt_rad+1,1:nbndlw+1)

             CALL rrtmg_lw( 1,                                                 &
                            nzt_rad-k_topo,                                    &
                            rrtm_icld,                                         &
                            rrtm_idrv,                                         &
                            rrtm_play(:,k_topo+1:),                   &
                            rrtm_plev(:,k_topo+1:),                   &
                            rrtm_tlay(:,k_topo+1:),                   &
                            rrtm_tlev(:,k_topo+1:),                   &
                            rrtm_tsfc,                                         &
                            rrtm_h2ovmr(:,k_topo+1:),                 &
                            rrtm_o3vmr(:,k_topo+1:),                  &
                            rrtm_co2vmr(:,k_topo+1:),                 &
                            rrtm_ch4vmr(:,k_topo+1:),                 &
                            rrtm_n2ovmr(:,k_topo+1:),                 &
                            rrtm_o2vmr(:,k_topo+1:),                  &
                            rrtm_cfc11vmr(:,k_topo+1:),               &
                            rrtm_cfc12vmr(:,k_topo+1:),               &
                            rrtm_cfc22vmr(:,k_topo+1:),               &
                            rrtm_ccl4vmr(:,k_topo+1:),                &
                            rrtm_emis,                                         &
                            rrtm_inflglw,                                      &
                            rrtm_iceflglw,                                     &
                            rrtm_liqflglw,                                     &
                            rrtm_cldfr(:,k_topo+1:),                  &
                            rrtm_lw_taucld_dum,                                &
                            rrtm_cicewp(:,k_topo+1:),                 &
                            rrtm_cliqwp(:,k_topo+1:),                 &
                            rrtm_reice(:,k_topo+1:),                  &
                            rrtm_reliq(:,k_topo+1:),                  &
                            rrtm_lw_tauaer_dum,                                &
                            rrtm_lwuflx(:,k_topo:),                   &
                            rrtm_lwdflx(:,k_topo:),                   &
                            rrtm_lwhr(:,k_topo+1:),                   &
                            rrtm_lwuflxc(:,k_topo:),                  &
                            rrtm_lwdflxc(:,k_topo:),                  &
                            rrtm_lwhrc(:,k_topo+1:),                  &
                            rrtm_lwuflx_dt(:,k_topo:),                &
                            rrtm_lwuflxc_dt(:,k_topo:) )

             DEALLOCATE ( rrtm_lw_taucld_dum )
             DEALLOCATE ( rrtm_lw_tauaer_dum )
!
!--          Save fluxes
             DO k = nzb, nzt+1
                rad_lw_in(k,:,:)  = rrtm_lwdflx(0,k)
                rad_lw_out(k,:,:) = rrtm_lwuflx(0,k)
             ENDDO
             rad_lw_in_diff(:,:) = rad_lw_in(k_topo,:,:)
!
!--          Save heating rates (convert from K/d to K/h).
!--          Further, even though an aggregated radiation is computed, map
!--          signle-column profiles on top of any topography, in order to
!--          obtain correct near surface radiation heating/cooling rates.
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   k_topo_l = topo_top_ind(j,i,0)
                   DO k = k_topo_l+1, nzt+1
                      rad_lw_hr(k,j,i)     = rrtm_lwhr(0,k-k_topo_l)  * d_hours_day
                      rad_lw_cs_hr(k,j,i)  = rrtm_lwhrc(0,k-k_topo_l) * d_hours_day
                   ENDDO
                ENDDO
             ENDDO

          ENDIF

          IF ( sw_radiation .AND. sun_up )  THEN
!
!--          Due to technical reasons, copy optical depths and other
!--          to dummy arguments which are allocated on the exact size as the
!--          rrtmg_sw is called.
!--          As one dimesion is allocated with zero size, compiler complains
!--          that rank of the array does not match that of the
!--          assumed-shaped arguments in the RRTMG library. In order to
!--          avoid this, write to dummy arguments and give pass the entire
!--          dummy array. Seems to be the only existing work-around.
             ALLOCATE( rrtm_sw_taucld_dum(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1) )
             ALLOCATE( rrtm_sw_ssacld_dum(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1) )
             ALLOCATE( rrtm_sw_asmcld_dum(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1) )
             ALLOCATE( rrtm_sw_fsfcld_dum(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1) )
             ALLOCATE( rrtm_sw_tauaer_dum(0:0,k_topo+1:nzt_rad+1,1:nbndsw+1) )
             ALLOCATE( rrtm_sw_ssaaer_dum(0:0,k_topo+1:nzt_rad+1,1:nbndsw+1) )
             ALLOCATE( rrtm_sw_asmaer_dum(0:0,k_topo+1:nzt_rad+1,1:nbndsw+1) )
             ALLOCATE( rrtm_sw_ecaer_dum(0:0,k_topo+1:nzt_rad+1,1:naerec+1)  )

             rrtm_sw_taucld_dum = rrtm_sw_taucld(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1)
             rrtm_sw_ssacld_dum = rrtm_sw_ssacld(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1)
             rrtm_sw_asmcld_dum = rrtm_sw_asmcld(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1)
             rrtm_sw_fsfcld_dum = rrtm_sw_fsfcld(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1)
             rrtm_sw_tauaer_dum = rrtm_sw_tauaer(0:0,k_topo+1:nzt_rad+1,1:nbndsw+1)
             rrtm_sw_ssaaer_dum = rrtm_sw_ssaaer(0:0,k_topo+1:nzt_rad+1,1:nbndsw+1)
             rrtm_sw_asmaer_dum = rrtm_sw_asmaer(0:0,k_topo+1:nzt_rad+1,1:nbndsw+1)
             rrtm_sw_ecaer_dum  = rrtm_sw_ecaer(0:0,k_topo+1:nzt_rad+1,1:naerec+1)

             CALL rrtmg_sw( 1,                                                 &
                            nzt_rad-k_topo,                                    &
                            rrtm_icld,                                         &
                            rrtm_iaer,                                         &
                            rrtm_play(:,k_topo+1:nzt_rad+1),                   &
                            rrtm_plev(:,k_topo+1:nzt_rad+2),                   &
                            rrtm_tlay(:,k_topo+1:nzt_rad+1),                   &
                            rrtm_tlev(:,k_topo+1:nzt_rad+2),                   &
                            rrtm_tsfc,                                         &
                            rrtm_h2ovmr(:,k_topo+1:nzt_rad+1),                 &
                            rrtm_o3vmr(:,k_topo+1:nzt_rad+1),                  &
                            rrtm_co2vmr(:,k_topo+1:nzt_rad+1),                 &
                            rrtm_ch4vmr(:,k_topo+1:nzt_rad+1),                 &
                            rrtm_n2ovmr(:,k_topo+1:nzt_rad+1),                 &
                            rrtm_o2vmr(:,k_topo+1:nzt_rad+1),                  &
                            rrtm_asdir,                                        &
                            rrtm_asdif,                                        &
                            rrtm_aldir,                                        &
                            rrtm_aldif,                                        &
                            zenith,                                            &
                            0.0_wp,                                            &
                            day_of_year,                                       &
                            solar_constant,                                    &
                            rrtm_inflgsw,                                      &
                            rrtm_iceflgsw,                                     &
                            rrtm_liqflgsw,                                     &
                            rrtm_cldfr(:,k_topo+1:nzt_rad+1),                  &
                            rrtm_sw_taucld_dum,                                &
                            rrtm_sw_ssacld_dum,                                &
                            rrtm_sw_asmcld_dum,                                &
                            rrtm_sw_fsfcld_dum,                                &
                            rrtm_cicewp(:,k_topo+1:nzt_rad+1),                 &
                            rrtm_cliqwp(:,k_topo+1:nzt_rad+1),                 &
                            rrtm_reice(:,k_topo+1:nzt_rad+1),                  &
                            rrtm_reliq(:,k_topo+1:nzt_rad+1),                  &
                            rrtm_sw_tauaer_dum,                                &
                            rrtm_sw_ssaaer_dum,                                &
                            rrtm_sw_asmaer_dum,                                &
                            rrtm_sw_ecaer_dum,                                 &
                            rrtm_swuflx(:,k_topo:nzt_rad+1),                   &
                            rrtm_swdflx(:,k_topo:nzt_rad+1),                   &
                            rrtm_swhr(:,k_topo+1:nzt_rad+1),                   &
                            rrtm_swuflxc(:,k_topo:nzt_rad+1),                  &
                            rrtm_swdflxc(:,k_topo:nzt_rad+1),                  &
                            rrtm_swhrc(:,k_topo+1:nzt_rad+1),                  &
                            rrtm_dirdflux(:,k_topo:nzt_rad+1),                 &
                            rrtm_difdflux(:,k_topo:nzt_rad+1) )

             DEALLOCATE( rrtm_sw_taucld_dum )
             DEALLOCATE( rrtm_sw_ssacld_dum )
             DEALLOCATE( rrtm_sw_asmcld_dum )
             DEALLOCATE( rrtm_sw_fsfcld_dum )
             DEALLOCATE( rrtm_sw_tauaer_dum )
             DEALLOCATE( rrtm_sw_ssaaer_dum )
             DEALLOCATE( rrtm_sw_asmaer_dum )
             DEALLOCATE( rrtm_sw_ecaer_dum )

!
!--          Save radiation fluxes for the entire depth of the model domain
             DO k = nzb, nzt+1
                rad_sw_in(k,:,:)  = rrtm_swdflx(0,k)
                rad_sw_out(k,:,:) = rrtm_swuflx(0,k)
             ENDDO
!--          Save direct and diffuse SW radiation at the surface (required by RTM)
             rad_sw_in_dir(:,:) = rrtm_dirdflux(0,k_topo)
             rad_sw_in_diff(:,:) = rrtm_difdflux(0,k_topo)

!
!--          Save heating rates (convert from K/d to K/s)
             DO k = nzb+1, nzt+1
                rad_sw_hr(k,:,:)     = rrtm_swhr(0,k)  * d_hours_day
                rad_sw_cs_hr(k,:,:)  = rrtm_swhrc(0,k) * d_hours_day
             ENDDO
!
!--       Solar radiation is zero during night
          ELSE
             rad_sw_in  = 0.0_wp
             rad_sw_out = 0.0_wp
             rad_sw_in_dir(:,:) = 0.0_wp
             rad_sw_in_diff(:,:) = 0.0_wp
          ENDIF
!
!--    RRTMG is called for each (j,i) grid point separately, starting at the
!--    highest topography level. Here no RTM is used since average_radiation is false
       ELSE
!
!--       Loop over all grid points
          DO i = nxl, nxr
             DO j = nys, nyn

!
!--             Prepare profiles of temperature and H2O volume mixing ratio
                DO  m = surf_lsm_h%start_index(j,i), surf_lsm_h%end_index(j,i)
                   rrtm_tlev(0,nzb+1) = surf_lsm_h%pt_surface(m) * exner(nzb)
                ENDDO
                DO  m = surf_usm_h%start_index(j,i), surf_usm_h%end_index(j,i)
                   rrtm_tlev(0,nzb+1) = surf_usm_h%pt_surface(m) * exner(nzb)
                ENDDO


                IF ( bulk_cloud_model )  THEN
                   DO k = nzb+1, nzt+1
                      rrtm_tlay(0,k) = pt(k,j,i) * exner(k)                    &
                                        + lv_d_cp * ql(k,j,i)
                      rrtm_h2ovmr(0,k) = mol_mass_air_d_wv * (q(k,j,i) - ql(k,j,i))
                   ENDDO
                ELSEIF ( cloud_droplets )  THEN
                   DO k = nzb+1, nzt+1
                      rrtm_tlay(0,k) = pt(k,j,i) * exner(k)                     &
                                        + lv_d_cp * ql(k,j,i)
                      rrtm_h2ovmr(0,k) = mol_mass_air_d_wv * q(k,j,i)
                   ENDDO
                ELSE
                   DO k = nzb+1, nzt+1
                      rrtm_tlay(0,k) = pt(k,j,i) * exner(k)
                   ENDDO

                   IF ( humidity )  THEN
                      DO k = nzb+1, nzt+1
                         rrtm_h2ovmr(0,k) =  mol_mass_air_d_wv * q(k,j,i)
                      ENDDO
                   ELSE
                      rrtm_h2ovmr(0,nzb+1:nzt+1) = 0.0_wp
                   ENDIF
                ENDIF

!
!--             Avoid temperature/humidity jumps at the top of the LES domain by
!--             linear interpolation from nzt+2 to nzt+7
                DO k = nzt+2, nzt+7
                   rrtm_tlay(0,k) = rrtm_tlay(0,nzt+1)                         &
                                 + ( rrtm_tlay(0,nzt+8) - rrtm_tlay(0,nzt+1) ) &
                                 / ( rrtm_play(0,nzt+8) - rrtm_play(0,nzt+1) ) &
                                 * ( rrtm_play(0,k)     - rrtm_play(0,nzt+1) )

                   rrtm_h2ovmr(0,k) = rrtm_h2ovmr(0,nzt+1)                     &
                              + ( rrtm_h2ovmr(0,nzt+8) - rrtm_h2ovmr(0,nzt+1) )&
                              / ( rrtm_play(0,nzt+8)   - rrtm_play(0,nzt+1)   )&
                              * ( rrtm_play(0,k)       - rrtm_play(0,nzt+1) )

                ENDDO

!--             Linear interpolate to zw grid
                DO k = nzb+2, nzt+8
                   rrtm_tlev(0,k)   = rrtm_tlay(0,k-1) + (rrtm_tlay(0,k) -     &
                                      rrtm_tlay(0,k-1))                        &
                                      / ( rrtm_play(0,k) - rrtm_play(0,k-1) )  &
                                      * ( rrtm_plev(0,k) - rrtm_play(0,k-1) )
                ENDDO


!
!--             Calculate liquid water path and cloud fraction for each column.
!--             Note that LWP is required in g/m2 instead of kg/kg m.
                rrtm_cldfr  = 0.0_wp
                rrtm_reliq  = 0.0_wp
                rrtm_cliqwp = 0.0_wp
                rrtm_icld   = 0

                IF ( bulk_cloud_model  .OR.  cloud_droplets )  THEN
                   DO k = nzb+1, nzt+1
                      rrtm_cliqwp(0,k) =  ql(k,j,i) * 1000.0_wp *              &
                                          (rrtm_plev(0,k) - rrtm_plev(0,k+1))  &
                                          * 100.0_wp / g

                      IF ( rrtm_cliqwp(0,k) > 0.0_wp )  THEN
                         rrtm_cldfr(0,k) = 1.0_wp
                         IF ( rrtm_icld == 0 )  rrtm_icld = 1

!
!--                      Calculate cloud droplet effective radius
                         IF ( bulk_cloud_model )  THEN
!
!--                         Calculete effective droplet radius. In case of using
!--                         cloud_scheme = 'morrison' and a non reasonable number
!--                         of cloud droplets the inital aerosol number
!--                         concentration is considered.
                            IF ( microphysics_morrison )  THEN
                               IF ( nc(k,j,i) > 1.0E-20_wp )  THEN
                                  nc_rad = nc(k,j,i)
                               ELSE
                                  nc_rad = na_init
                               ENDIF
                            ELSE
                               nc_rad = nc_const
                            ENDIF

                            rrtm_reliq(0,k) = 1.0E6_wp * ( 3.0_wp * ql(k,j,i)     &
                                              * rho_surface                       &
                                              / ( 4.0_wp * pi * nc_rad * rho_l )  &
                                              )**0.33333333333333_wp              &
                                              * EXP( LOG( sigma_gc )**2 )

                         ELSEIF ( cloud_droplets )  THEN
                            number_of_particles = prt_count(k,j,i)

                            IF (number_of_particles <= 0)  CYCLE
                            particles => grid_particles(k,j,i)%particles(1:number_of_particles)
                            s_r2 = 0.0_wp
                            s_r3 = 0.0_wp

                            DO  n = 1, number_of_particles
                               IF ( particles(n)%particle_mask )  THEN
                                  s_r2 = s_r2 + particles(n)%radius**2 *       &
                                         particles(n)%weight_factor
                                  s_r3 = s_r3 + particles(n)%radius**3 *       &
                                         particles(n)%weight_factor
                               ENDIF
                            ENDDO

                            IF ( s_r2 > 0.0_wp )  rrtm_reliq(0,k) = s_r3 / s_r2

                         ENDIF

!
!--                      Limit effective radius
                         IF ( rrtm_reliq(0,k) > 0.0_wp )  THEN
                            rrtm_reliq(0,k) = MAX(rrtm_reliq(0,k),2.5_wp)
                            rrtm_reliq(0,k) = MIN(rrtm_reliq(0,k),60.0_wp)
                        ENDIF
                      ENDIF
                   ENDDO
                ENDIF

!
!--             Write surface emissivity and surface temperature at current
!--             surface element on RRTMG-shaped array.
!--             Please note, as RRTMG is a single column model, surface attributes
!--             are only obtained from horizontally aligned surfaces (for
!--             simplicity). Taking surface attributes from horizontal and
!--             vertical walls would lead to multiple solutions.
!--             Moreover, for natural- and urban-type surfaces, several surface
!--             classes can exist at a surface element next to each other.
!--             To obtain bulk parameters, apply a weighted average for these
!--             surfaces.
                DO  m = surf_lsm_h%start_index(j,i), surf_lsm_h%end_index(j,i)
                   rrtm_emis = surf_lsm_h%frac(m,ind_veg_wall)  *              &
                               surf_lsm_h%emissivity(m,ind_veg_wall)  +        &
                               surf_lsm_h%frac(m,ind_pav_green) *              &
                               surf_lsm_h%emissivity(m,ind_pav_green) +        &
                               surf_lsm_h%frac(m,ind_wat_win)   *              &
                               surf_lsm_h%emissivity(m,ind_wat_win)
                   rrtm_tsfc = surf_lsm_h%pt_surface(m) * exner(nzb)
                ENDDO
                DO  m = surf_usm_h%start_index(j,i), surf_usm_h%end_index(j,i)
                   rrtm_emis = surf_usm_h%frac(m,ind_veg_wall)  *              &
                               surf_usm_h%emissivity(m,ind_veg_wall)  +        &
                               surf_usm_h%frac(m,ind_pav_green) *              &
                               surf_usm_h%emissivity(m,ind_pav_green) +        &
                               surf_usm_h%frac(m,ind_wat_win)   *              &
                               surf_usm_h%emissivity(m,ind_wat_win)
                   rrtm_tsfc = surf_usm_h%pt_surface(m) * exner(nzb)
                ENDDO
!
!--             Obtain topography top index (lower bound of RRTMG)
                k_topo = topo_top_ind(j,i,0)

                IF ( lw_radiation )  THEN
!
!--                Due to technical reasons, copy optical depth to dummy arguments
!--                which are allocated on the exact size as the rrtmg_lw is called.
!--                As one dimesion is allocated with zero size, compiler complains
!--                that rank of the array does not match that of the
!--                assumed-shaped arguments in the RRTMG library. In order to
!--                avoid this, write to dummy arguments and give pass the entire
!--                dummy array. Seems to be the only existing work-around.
                   ALLOCATE( rrtm_lw_taucld_dum(1:nbndlw+1,0:0,k_topo+1:nzt_rad+1) )
                   ALLOCATE( rrtm_lw_tauaer_dum(0:0,k_topo+1:nzt_rad+1,1:nbndlw+1) )

                   rrtm_lw_taucld_dum =                                        &
                               rrtm_lw_taucld(1:nbndlw+1,0:0,k_topo+1:nzt_rad+1)
                   rrtm_lw_tauaer_dum =                                        &
                               rrtm_lw_tauaer(0:0,k_topo+1:nzt_rad+1,1:nbndlw+1)

                   CALL rrtmg_lw( 1,                                           &
                                  nzt_rad-k_topo,                              &
                                  rrtm_icld,                                   &
                                  rrtm_idrv,                                   &
                                  rrtm_play(:,k_topo+1:nzt_rad+1),             &
                                  rrtm_plev(:,k_topo+1:nzt_rad+2),             &
                                  rrtm_tlay(:,k_topo+1:nzt_rad+1),             &
                                  rrtm_tlev(:,k_topo+1:nzt_rad+2),             &
                                  rrtm_tsfc,                                   &
                                  rrtm_h2ovmr(:,k_topo+1:nzt_rad+1),           &
                                  rrtm_o3vmr(:,k_topo+1:nzt_rad+1),            &
                                  rrtm_co2vmr(:,k_topo+1:nzt_rad+1),           &
                                  rrtm_ch4vmr(:,k_topo+1:nzt_rad+1),           &
                                  rrtm_n2ovmr(:,k_topo+1:nzt_rad+1),           &
                                  rrtm_o2vmr(:,k_topo+1:nzt_rad+1),            &
                                  rrtm_cfc11vmr(:,k_topo+1:nzt_rad+1),         &
                                  rrtm_cfc12vmr(:,k_topo+1:nzt_rad+1),         &
                                  rrtm_cfc22vmr(:,k_topo+1:nzt_rad+1),         &
                                  rrtm_ccl4vmr(:,k_topo+1:nzt_rad+1),          &
                                  rrtm_emis,                                   &
                                  rrtm_inflglw,                                &
                                  rrtm_iceflglw,                               &
                                  rrtm_liqflglw,                               &
                                  rrtm_cldfr(:,k_topo+1:nzt_rad+1),            &
                                  rrtm_lw_taucld_dum,                          &
                                  rrtm_cicewp(:,k_topo+1:nzt_rad+1),           &
                                  rrtm_cliqwp(:,k_topo+1:nzt_rad+1),           &
                                  rrtm_reice(:,k_topo+1:nzt_rad+1),            &
                                  rrtm_reliq(:,k_topo+1:nzt_rad+1),            &
                                  rrtm_lw_tauaer_dum,                          &
                                  rrtm_lwuflx(:,k_topo:nzt_rad+1),             &
                                  rrtm_lwdflx(:,k_topo:nzt_rad+1),             &
                                  rrtm_lwhr(:,k_topo+1:nzt_rad+1),             &
                                  rrtm_lwuflxc(:,k_topo:nzt_rad+1),            &
                                  rrtm_lwdflxc(:,k_topo:nzt_rad+1),            &
                                  rrtm_lwhrc(:,k_topo+1:nzt_rad+1),            &
                                  rrtm_lwuflx_dt(:,k_topo:nzt_rad+1),          &
                                  rrtm_lwuflxc_dt(:,k_topo:nzt_rad+1) )

                   DEALLOCATE ( rrtm_lw_taucld_dum )
                   DEALLOCATE ( rrtm_lw_tauaer_dum )
!
!--                Save fluxes
                   DO k = k_topo, nzt+1
                      rad_lw_in(k,j,i)  = rrtm_lwdflx(0,k)
                      rad_lw_out(k,j,i) = rrtm_lwuflx(0,k)
                   ENDDO

!
!--                Save heating rates (convert from K/d to K/h)
                   DO k = k_topo+1, nzt+1
                      rad_lw_hr(k,j,i)     = rrtm_lwhr(0,k-k_topo)  * d_hours_day
                      rad_lw_cs_hr(k,j,i)  = rrtm_lwhrc(0,k-k_topo) * d_hours_day
                   ENDDO

!
!--                Save surface radiative fluxes and change in LW heating rate
!--                onto respective surface elements
!--                Horizontal surfaces
                   DO  m = surf_lsm_h%start_index(j,i),                        &
                           surf_lsm_h%end_index(j,i)
                      surf_lsm_h%rad_lw_in(m)           = rrtm_lwdflx(0,k_topo)
                      surf_lsm_h%rad_lw_out(m)          = rrtm_lwuflx(0,k_topo)
                      surf_lsm_h%rad_lw_out_change_0(m) = rrtm_lwuflx_dt(0,k_topo)
                   ENDDO
                   DO  m = surf_usm_h%start_index(j,i),                        &
                           surf_usm_h%end_index(j,i)
                      surf_usm_h%rad_lw_in(m)           = rrtm_lwdflx(0,k_topo)
                      surf_usm_h%rad_lw_out(m)          = rrtm_lwuflx(0,k_topo)
                      surf_usm_h%rad_lw_out_change_0(m) = rrtm_lwuflx_dt(0,k_topo)
                   ENDDO
!
!--                Vertical surfaces. Fluxes are obtain at vertical level of the
!--                respective surface element
                   DO  l = 0, 3
                      DO  m = surf_lsm_v(l)%start_index(j,i),                  &
                              surf_lsm_v(l)%end_index(j,i)
                         k                                    = surf_lsm_v(l)%k(m)
                         surf_lsm_v(l)%rad_lw_in(m)           = rrtm_lwdflx(0,k)
                         surf_lsm_v(l)%rad_lw_out(m)          = rrtm_lwuflx(0,k)
                         surf_lsm_v(l)%rad_lw_out_change_0(m) = rrtm_lwuflx_dt(0,k)
                      ENDDO
                      DO  m = surf_usm_v(l)%start_index(j,i),                  &
                              surf_usm_v(l)%end_index(j,i)
                         k                                    = surf_usm_v(l)%k(m)
                         surf_usm_v(l)%rad_lw_in(m)           = rrtm_lwdflx(0,k)
                         surf_usm_v(l)%rad_lw_out(m)          = rrtm_lwuflx(0,k)
                         surf_usm_v(l)%rad_lw_out_change_0(m) = rrtm_lwuflx_dt(0,k)
                      ENDDO
                   ENDDO

                ENDIF

                IF ( sw_radiation .AND. sun_up )  THEN
!
!--                Get albedo for direct/diffusive long/shortwave radiation at
!--                current (y,x)-location from surface variables.
!--                Only obtain it from horizontal surfaces, as RRTMG is a single
!--                column model
!--                (Please note, only one loop will entered, controlled by
!--                start-end index.)
                   DO  m = surf_lsm_h%start_index(j,i),                        &
                           surf_lsm_h%end_index(j,i)
                      rrtm_asdir(1)  = SUM( surf_lsm_h%frac(m,:) *             &
                                            surf_lsm_h%rrtm_asdir(m,:) )
                      rrtm_asdif(1)  = SUM( surf_lsm_h%frac(m,:) *             &
                                            surf_lsm_h%rrtm_asdif(m,:) )
                      rrtm_aldir(1)  = SUM( surf_lsm_h%frac(m,:) *             &
                                            surf_lsm_h%rrtm_aldir(m,:) )
                      rrtm_aldif(1)  = SUM( surf_lsm_h%frac(m,:) *             &
                                            surf_lsm_h%rrtm_aldif(m,:) )
                   ENDDO
                   DO  m = surf_usm_h%start_index(j,i),                        &
                           surf_usm_h%end_index(j,i)
                      rrtm_asdir(1)  = SUM( surf_usm_h%frac(m,:) *             &
                                            surf_usm_h%rrtm_asdir(m,:) )
                      rrtm_asdif(1)  = SUM( surf_usm_h%frac(m,:) *             &
                                            surf_usm_h%rrtm_asdif(m,:) )
                      rrtm_aldir(1)  = SUM( surf_usm_h%frac(m,:) *             &
                                            surf_usm_h%rrtm_aldir(m,:) )
                      rrtm_aldif(1)  = SUM( surf_usm_h%frac(m,:) *             &
                                            surf_usm_h%rrtm_aldif(m,:) )
                   ENDDO
!
!--                Due to technical reasons, copy optical depths and other
!--                to dummy arguments which are allocated on the exact size as the
!--                rrtmg_sw is called.
!--                As one dimesion is allocated with zero size, compiler complains
!--                that rank of the array does not match that of the
!--                assumed-shaped arguments in the RRTMG library. In order to
!--                avoid this, write to dummy arguments and give pass the entire
!--                dummy array. Seems to be the only existing work-around.
                   ALLOCATE( rrtm_sw_taucld_dum(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1) )
                   ALLOCATE( rrtm_sw_ssacld_dum(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1) )
                   ALLOCATE( rrtm_sw_asmcld_dum(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1) )
                   ALLOCATE( rrtm_sw_fsfcld_dum(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1) )
                   ALLOCATE( rrtm_sw_tauaer_dum(0:0,k_topo+1:nzt_rad+1,1:nbndsw+1) )
                   ALLOCATE( rrtm_sw_ssaaer_dum(0:0,k_topo+1:nzt_rad+1,1:nbndsw+1) )
                   ALLOCATE( rrtm_sw_asmaer_dum(0:0,k_topo+1:nzt_rad+1,1:nbndsw+1) )
                   ALLOCATE( rrtm_sw_ecaer_dum(0:0,k_topo+1:nzt_rad+1,1:naerec+1)  )

                   rrtm_sw_taucld_dum = rrtm_sw_taucld(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1)
                   rrtm_sw_ssacld_dum = rrtm_sw_ssacld(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1)
                   rrtm_sw_asmcld_dum = rrtm_sw_asmcld(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1)
                   rrtm_sw_fsfcld_dum = rrtm_sw_fsfcld(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1)
                   rrtm_sw_tauaer_dum = rrtm_sw_tauaer(0:0,k_topo+1:nzt_rad+1,1:nbndsw+1)
                   rrtm_sw_ssaaer_dum = rrtm_sw_ssaaer(0:0,k_topo+1:nzt_rad+1,1:nbndsw+1)
                   rrtm_sw_asmaer_dum = rrtm_sw_asmaer(0:0,k_topo+1:nzt_rad+1,1:nbndsw+1)
                   rrtm_sw_ecaer_dum  = rrtm_sw_ecaer(0:0,k_topo+1:nzt_rad+1,1:naerec+1)

                   CALL rrtmg_sw( 1,                                           &
                                  nzt_rad-k_topo,                              &
                                  rrtm_icld,                                   &
                                  rrtm_iaer,                                   &
                                  rrtm_play(:,k_topo+1:nzt_rad+1),             &
                                  rrtm_plev(:,k_topo+1:nzt_rad+2),             &
                                  rrtm_tlay(:,k_topo+1:nzt_rad+1),             &
                                  rrtm_tlev(:,k_topo+1:nzt_rad+2),             &
                                  rrtm_tsfc,                                   &
                                  rrtm_h2ovmr(:,k_topo+1:nzt_rad+1),           &
                                  rrtm_o3vmr(:,k_topo+1:nzt_rad+1),            &
                                  rrtm_co2vmr(:,k_topo+1:nzt_rad+1),           &
                                  rrtm_ch4vmr(:,k_topo+1:nzt_rad+1),           &
                                  rrtm_n2ovmr(:,k_topo+1:nzt_rad+1),           &
                                  rrtm_o2vmr(:,k_topo+1:nzt_rad+1),            &
                                  rrtm_asdir,                                  &
                                  rrtm_asdif,                                  &
                                  rrtm_aldir,                                  &
                                  rrtm_aldif,                                  &
                                  zenith,                                      &
                                  0.0_wp,                                      &
                                  day_of_year,                                 &
                                  solar_constant,                              &
                                  rrtm_inflgsw,                                &
                                  rrtm_iceflgsw,                               &
                                  rrtm_liqflgsw,                               &
                                  rrtm_cldfr(:,k_topo+1:nzt_rad+1),            &
                                  rrtm_sw_taucld_dum,                          &
                                  rrtm_sw_ssacld_dum,                          &
                                  rrtm_sw_asmcld_dum,                          &
                                  rrtm_sw_fsfcld_dum,                          &
                                  rrtm_cicewp(:,k_topo+1:nzt_rad+1),           &
                                  rrtm_cliqwp(:,k_topo+1:nzt_rad+1),           &
                                  rrtm_reice(:,k_topo+1:nzt_rad+1),            &
                                  rrtm_reliq(:,k_topo+1:nzt_rad+1),            &
                                  rrtm_sw_tauaer_dum,                          &
                                  rrtm_sw_ssaaer_dum,                          &
                                  rrtm_sw_asmaer_dum,                          &
                                  rrtm_sw_ecaer_dum,                           &
                                  rrtm_swuflx(:,k_topo:nzt_rad+1),             &
                                  rrtm_swdflx(:,k_topo:nzt_rad+1),             &
                                  rrtm_swhr(:,k_topo+1:nzt_rad+1),             &
                                  rrtm_swuflxc(:,k_topo:nzt_rad+1),            &
                                  rrtm_swdflxc(:,k_topo:nzt_rad+1),            &
                                  rrtm_swhrc(:,k_topo+1:nzt_rad+1),            &
                                  rrtm_dirdflux(:,k_topo:nzt_rad+1),           &
                                  rrtm_difdflux(:,k_topo:nzt_rad+1) )

                   DEALLOCATE( rrtm_sw_taucld_dum )
                   DEALLOCATE( rrtm_sw_ssacld_dum )
                   DEALLOCATE( rrtm_sw_asmcld_dum )
                   DEALLOCATE( rrtm_sw_fsfcld_dum )
                   DEALLOCATE( rrtm_sw_tauaer_dum )
                   DEALLOCATE( rrtm_sw_ssaaer_dum )
                   DEALLOCATE( rrtm_sw_asmaer_dum )
                   DEALLOCATE( rrtm_sw_ecaer_dum )
!
!--                Save fluxes
                   DO k = nzb, nzt+1
                      rad_sw_in(k,j,i)  = rrtm_swdflx(0,k)
                      rad_sw_out(k,j,i) = rrtm_swuflx(0,k)
                   ENDDO
!
!--                Save heating rates (convert from K/d to K/s)
                   DO k = nzb+1, nzt+1
                      rad_sw_hr(k,j,i)     = rrtm_swhr(0,k)  * d_hours_day
                      rad_sw_cs_hr(k,j,i)  = rrtm_swhrc(0,k) * d_hours_day
                   ENDDO

!
!--                Save surface radiative fluxes onto respective surface elements
!--                Horizontal surfaces
                   DO  m = surf_lsm_h%start_index(j,i),                        &
                           surf_lsm_h%end_index(j,i)
                      surf_lsm_h%rad_sw_in(m)     = rrtm_swdflx(0,k_topo)
                      surf_lsm_h%rad_sw_out(m)    = rrtm_swuflx(0,k_topo)
                   ENDDO
                   DO  m = surf_usm_h%start_index(j,i),                        &
                           surf_usm_h%end_index(j,i)
                      surf_usm_h%rad_sw_in(m)     = rrtm_swdflx(0,k_topo)
                      surf_usm_h%rad_sw_out(m)    = rrtm_swuflx(0,k_topo)
                   ENDDO
!
!--                Vertical surfaces. Fluxes are obtain at respective vertical
!--                level of the surface element
                   DO  l = 0, 3
                      DO  m = surf_lsm_v(l)%start_index(j,i),                  &
                              surf_lsm_v(l)%end_index(j,i)
                         k                           = surf_lsm_v(l)%k(m)
                         surf_lsm_v(l)%rad_sw_in(m)  = rrtm_swdflx(0,k)
                         surf_lsm_v(l)%rad_sw_out(m) = rrtm_swuflx(0,k)
                      ENDDO
                      DO  m = surf_usm_v(l)%start_index(j,i),                  &
                              surf_usm_v(l)%end_index(j,i)
                         k                           = surf_usm_v(l)%k(m)
                         surf_usm_v(l)%rad_sw_in(m)  = rrtm_swdflx(0,k)
                         surf_usm_v(l)%rad_sw_out(m) = rrtm_swuflx(0,k)
                      ENDDO
                   ENDDO
!
!--             Solar radiation is zero during night
                ELSE
                   rad_sw_in  = 0.0_wp
                   rad_sw_out = 0.0_wp
!--             !!!!!!!! ATTENSION !!!!!!!!!!!!!!!
!--             Surface radiative fluxes should be also set to zero here
!--                Save surface radiative fluxes onto respective surface elements
!--                Horizontal surfaces
                   DO  m = surf_lsm_h%start_index(j,i),                        &
                           surf_lsm_h%end_index(j,i)
                      surf_lsm_h%rad_sw_in(m)     = 0.0_wp
                      surf_lsm_h%rad_sw_out(m)    = 0.0_wp
                   ENDDO
                   DO  m = surf_usm_h%start_index(j,i),                        &
                           surf_usm_h%end_index(j,i)
                      surf_usm_h%rad_sw_in(m)     = 0.0_wp
                      surf_usm_h%rad_sw_out(m)    = 0.0_wp
                   ENDDO
!
!--                Vertical surfaces. Fluxes are obtain at respective vertical
!--                level of the surface element
                   DO  l = 0, 3
                      DO  m = surf_lsm_v(l)%start_index(j,i),                  &
                              surf_lsm_v(l)%end_index(j,i)
                         k                           = surf_lsm_v(l)%k(m)
                         surf_lsm_v(l)%rad_sw_in(m)  = 0.0_wp
                         surf_lsm_v(l)%rad_sw_out(m) = 0.0_wp
                      ENDDO
                      DO  m = surf_usm_v(l)%start_index(j,i),                  &
                              surf_usm_v(l)%end_index(j,i)
                         k                           = surf_usm_v(l)%k(m)
                         surf_usm_v(l)%rad_sw_in(m)  = 0.0_wp
                         surf_usm_v(l)%rad_sw_out(m) = 0.0_wp
                      ENDDO
                   ENDDO
                ENDIF

             ENDDO
          ENDDO

       ENDIF
!
!--    Finally, calculate surface net radiation for surface elements.
       IF (  .NOT.  radiation_interactions  ) THEN
!--       First, for horizontal surfaces
          DO  m = 1, surf_lsm_h%ns
             surf_lsm_h%rad_net(m) = surf_lsm_h%rad_sw_in(m)                   &
                                   - surf_lsm_h%rad_sw_out(m)                  &
                                   + surf_lsm_h%rad_lw_in(m)                   &
                                   - surf_lsm_h%rad_lw_out(m)
          ENDDO
          DO  m = 1, surf_usm_h%ns
             surf_usm_h%rad_net(m) = surf_usm_h%rad_sw_in(m)                   &
                                   - surf_usm_h%rad_sw_out(m)                  &
                                   + surf_usm_h%rad_lw_in(m)                   &
                                   - surf_usm_h%rad_lw_out(m)
          ENDDO
!
!--       Vertical surfaces.
!--       Todo: weight with azimuth and zenith angle according to their orientation!
          DO  l = 0, 3
             DO  m = 1, surf_lsm_v(l)%ns
                surf_lsm_v(l)%rad_net(m) = surf_lsm_v(l)%rad_sw_in(m)          &
                                         - surf_lsm_v(l)%rad_sw_out(m)         &
                                         + surf_lsm_v(l)%rad_lw_in(m)          &
                                         - surf_lsm_v(l)%rad_lw_out(m)
             ENDDO
             DO  m = 1, surf_usm_v(l)%ns
                surf_usm_v(l)%rad_net(m) = surf_usm_v(l)%rad_sw_in(m)          &
                                         - surf_usm_v(l)%rad_sw_out(m)         &
                                         + surf_usm_v(l)%rad_lw_in(m)          &
                                         - surf_usm_v(l)%rad_lw_out(m)
             ENDDO
          ENDDO
       ENDIF


       CALL exchange_horiz( rad_lw_in,  nbgp )
       CALL exchange_horiz( rad_lw_out, nbgp )
       CALL exchange_horiz( rad_lw_hr,    nbgp )
       CALL exchange_horiz( rad_lw_cs_hr, nbgp )

       CALL exchange_horiz( rad_sw_in,  nbgp )
       CALL exchange_horiz( rad_sw_out, nbgp )
       CALL exchange_horiz( rad_sw_hr,    nbgp )
       CALL exchange_horiz( rad_sw_cs_hr, nbgp )

#endif

    END SUBROUTINE radiation_rrtmg


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate the cosine of the zenith angle (variable is called zenith)
!------------------------------------------------------------------------------!
    SUBROUTINE calc_zenith( day_of_year, second_of_day )

       USE palm_date_time_mod,                                                 &
           ONLY:  seconds_per_day

       IMPLICIT NONE

       INTEGER(iwp), INTENT(IN) ::  day_of_year    !< day of the year

       REAL(wp)                 ::  declination    !< solar declination angle
       REAL(wp)                 ::  hour_angle     !< solar hour angle
       REAL(wp),     INTENT(IN) ::  second_of_day  !< current time of the day in UTC

!
!--    Calculate solar declination and hour angle
       declination = ASIN( decl_1 * SIN(decl_2 * REAL(day_of_year, KIND=wp) - decl_3) )
       hour_angle  = 2.0_wp * pi * ( second_of_day / seconds_per_day ) + lon - pi

!
!--    Calculate cosine of solar zenith angle
       cos_zenith = SIN(lat) * SIN(declination) + COS(lat) * COS(declination)   &
                                            * COS(hour_angle)
       cos_zenith = MAX(0.0_wp,cos_zenith)

!
!--    Calculate solar directional vector
       IF ( sun_direction )  THEN

!
!--       Direction in longitudes equals to sin(solar_azimuth) * sin(zenith)
          sun_dir_lon = -SIN(hour_angle) * COS(declination)

!
!--       Direction in latitues equals to cos(solar_azimuth) * sin(zenith)
          sun_dir_lat = SIN(declination) * COS(lat) - COS(hour_angle) &
                              * COS(declination) * SIN(lat)
       ENDIF

!
!--    Check if the sun is up (otheriwse shortwave calculations can be skipped)
       IF ( cos_zenith > 0.0_wp )  THEN
          sun_up = .TRUE.
       ELSE
          sun_up = .FALSE.
       END IF

    END SUBROUTINE calc_zenith

#if defined ( __rrtmg ) && defined ( __netcdf )
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculates surface albedo components based on Briegleb (1992) and
!> Briegleb et al. (1986)
!------------------------------------------------------------------------------!
    SUBROUTINE calc_albedo( surf )

        IMPLICIT NONE

        INTEGER(iwp)    ::  ind_type !< running index surface tiles
        INTEGER(iwp)    ::  m        !< running index surface elements

        TYPE(surf_type) ::  surf !< treated surfaces

        IF ( sun_up  .AND.  .NOT. average_radiation )  THEN

           DO  m = 1, surf%ns
!
!--           Loop over surface elements
              DO  ind_type = 0, SIZE( surf%albedo_type, 2 ) - 1

!
!--              Ocean
                 IF ( surf%albedo_type(m,ind_type) == 1 )  THEN
                    surf%rrtm_aldir(m,ind_type) = 0.026_wp /                    &
                                                ( cos_zenith**1.7_wp + 0.065_wp )&
                                     + 0.15_wp * ( cos_zenith - 0.1_wp )         &
                                               * ( cos_zenith - 0.5_wp )         &
                                               * ( cos_zenith - 1.0_wp )
                    surf%rrtm_asdir(m,ind_type) = surf%rrtm_aldir(m,ind_type)
!
!--              Snow
                 ELSEIF ( surf%albedo_type(m,ind_type) == 16 )  THEN
                    IF ( cos_zenith < 0.5_wp )  THEN
                       surf%rrtm_aldir(m,ind_type) =                           &
                                 0.5_wp * ( 1.0_wp - surf%aldif(im,ind_type) )  &
                                        * ( ( 3.0_wp / ( 1.0_wp + 4.0_wp       &
                                        * cos_zenith ) ) - 1.0_wp )
                       surf%rrtm_asdir(m,ind_type) =                           &
                                 0.5_wp * ( 1.0_wp - surf%asdif(m,ind_type) )  &
                                        * ( ( 3.0_wp / ( 1.0_wp + 4.0_wp       &
                                        * cos_zenith ) ) - 1.0_wp )

                       surf%rrtm_aldir(m,ind_type) =                           &
                                       MIN(0.98_wp, surf%rrtm_aldir(m,ind_type))
                       surf%rrtm_asdir(m,ind_type) =                           &
                                       MIN(0.98_wp, surf%rrtm_asdir(m,ind_type))
                    ELSE
                       surf%rrtm_aldir(m,ind_type) = surf%aldif(m,ind_type)
                       surf%rrtm_asdir(m,ind_type) = surf%asdif(m,ind_type)
                    ENDIF
!
!--              Sea ice
                 ELSEIF ( surf%albedo_type(m,ind_type) == 15 )  THEN
                    surf%rrtm_aldir(m,ind_type) = surf%aldif(m,ind_type)
                    surf%rrtm_asdir(m,ind_type) = surf%asdif(m,ind_type)

!
!--              Asphalt
                 ELSEIF ( surf%albedo_type(m,ind_type) == 17 )  THEN
                    surf%rrtm_aldir(m,ind_type) = surf%aldif(m,ind_type)
                    surf%rrtm_asdir(m,ind_type) = surf%asdif(m,ind_type)


!
!--              Bare soil
                 ELSEIF ( surf%albedo_type(m,ind_type) == 18 )  THEN
                    surf%rrtm_aldir(m,ind_type) = surf%aldif(m,ind_type)
                    surf%rrtm_asdir(m,ind_type) = surf%asdif(m,ind_type)

!
!--              Land surfaces
                 ELSE
                    SELECT CASE ( surf%albedo_type(m,ind_type) )

!
!--                    Surface types with strong zenith dependence
                       CASE ( 1, 2, 3, 4, 11, 12, 13 )
                          surf%rrtm_aldir(m,ind_type) =                        &
                                surf%aldif(m,ind_type) * 1.4_wp /              &
                                           ( 1.0_wp + 0.8_wp * cos_zenith )
                          surf%rrtm_asdir(m,ind_type) =                        &
                                surf%asdif(m,ind_type) * 1.4_wp /              &
                                           ( 1.0_wp + 0.8_wp * cos_zenith )
!
!--                    Surface types with weak zenith dependence
                       CASE ( 5, 6, 7, 8, 9, 10, 14 )
                          surf%rrtm_aldir(m,ind_type) =                        &
                                surf%aldif(m,ind_type) * 1.1_wp /              &
                                           ( 1.0_wp + 0.2_wp * cos_zenith )
                          surf%rrtm_asdir(m,ind_type) =                        &
                                surf%asdif(m,ind_type) * 1.1_wp /              &
                                           ( 1.0_wp + 0.2_wp * cos_zenith )

                       CASE DEFAULT

                    END SELECT
                 ENDIF
!
!--              Diffusive albedo is taken from Table 2
                 surf%rrtm_aldif(m,ind_type) = surf%aldif(m,ind_type)
                 surf%rrtm_asdif(m,ind_type) = surf%asdif(m,ind_type)
              ENDDO
           ENDDO
!
!--     Set albedo in case of average radiation
        ELSEIF ( sun_up  .AND.  average_radiation )  THEN
           surf%rrtm_asdir = albedo_urb
           surf%rrtm_asdif = albedo_urb
           surf%rrtm_aldir = albedo_urb
           surf%rrtm_aldif = albedo_urb
!
!--     Darkness
        ELSE
           surf%rrtm_aldir = 0.0_wp
           surf%rrtm_asdir = 0.0_wp
           surf%rrtm_aldif = 0.0_wp
           surf%rrtm_asdif = 0.0_wp
        ENDIF

    END SUBROUTINE calc_albedo

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Read sounding data (pressure and temperature) from RADIATION_DATA.
!------------------------------------------------------------------------------!
    SUBROUTINE read_sounding_data

       IMPLICIT NONE

       INTEGER(iwp) :: id,           & !< NetCDF id of input file
                       id_dim_zrad,  & !< pressure level id in the NetCDF file
                       id_var,       & !< NetCDF variable id
                       k,            & !< loop index
                       nz_snd,       & !< number of vertical levels in the sounding data
                       nz_snd_start, & !< start vertical index for sounding data to be used
                       nz_snd_end      !< end vertical index for souding data to be used

       REAL(wp) :: t_surface           !< actual surface temperature

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  hyp_snd_tmp, & !< temporary hydrostatic pressure profile (sounding)
                                               t_snd_tmp      !< temporary temperature profile (sounding)

!
!--    In case of updates, deallocate arrays first (sufficient to check one
!--    array as the others are automatically allocated). This is required
!--    because nzt_rad might change during the update
       IF ( ALLOCATED ( hyp_snd ) )  THEN
          DEALLOCATE( hyp_snd )
          DEALLOCATE( t_snd )
          DEALLOCATE ( rrtm_play )
          DEALLOCATE ( rrtm_plev )
          DEALLOCATE ( rrtm_tlay )
          DEALLOCATE ( rrtm_tlev )

          DEALLOCATE ( rrtm_cicewp )
          DEALLOCATE ( rrtm_cldfr )
          DEALLOCATE ( rrtm_cliqwp )
          DEALLOCATE ( rrtm_reice )
          DEALLOCATE ( rrtm_reliq )
          DEALLOCATE ( rrtm_lw_taucld )
          DEALLOCATE ( rrtm_lw_tauaer )

          DEALLOCATE ( rrtm_lwdflx  )
          DEALLOCATE ( rrtm_lwdflxc )
          DEALLOCATE ( rrtm_lwuflx  )
          DEALLOCATE ( rrtm_lwuflxc )
          DEALLOCATE ( rrtm_lwuflx_dt )
          DEALLOCATE ( rrtm_lwuflxc_dt )
          DEALLOCATE ( rrtm_lwhr  )
          DEALLOCATE ( rrtm_lwhrc )

          DEALLOCATE ( rrtm_sw_taucld )
          DEALLOCATE ( rrtm_sw_ssacld )
          DEALLOCATE ( rrtm_sw_asmcld )
          DEALLOCATE ( rrtm_sw_fsfcld )
          DEALLOCATE ( rrtm_sw_tauaer )
          DEALLOCATE ( rrtm_sw_ssaaer )
          DEALLOCATE ( rrtm_sw_asmaer )
          DEALLOCATE ( rrtm_sw_ecaer )

          DEALLOCATE ( rrtm_swdflx  )
          DEALLOCATE ( rrtm_swdflxc )
          DEALLOCATE ( rrtm_swuflx  )
          DEALLOCATE ( rrtm_swuflxc )
          DEALLOCATE ( rrtm_swhr  )
          DEALLOCATE ( rrtm_swhrc )
          DEALLOCATE ( rrtm_dirdflux )
          DEALLOCATE ( rrtm_difdflux )

       ENDIF

!
!--    Open file for reading
       nc_stat = NF90_OPEN( rrtm_input_file, NF90_NOWRITE, id )
       CALL netcdf_handle_error_rad( 'read_sounding_data', 549 )

!
!--    Inquire dimension of z axis and save in nz_snd
       nc_stat = NF90_INQ_DIMID( id, "Pressure", id_dim_zrad )
       nc_stat = NF90_INQUIRE_DIMENSION( id, id_dim_zrad, len = nz_snd )
       CALL netcdf_handle_error_rad( 'read_sounding_data', 551 )

!
! !--    Allocate temporary array for storing pressure data
       ALLOCATE( hyp_snd_tmp(1:nz_snd) )
       hyp_snd_tmp = 0.0_wp


!--    Read pressure from file
       nc_stat = NF90_INQ_VARID( id, "Pressure", id_var )
       nc_stat = NF90_GET_VAR( id, id_var, hyp_snd_tmp(:), start = (/1/),      &
                               count = (/nz_snd/) )
       CALL netcdf_handle_error_rad( 'read_sounding_data', 552 )

!
!--    Allocate temporary array for storing temperature data
       ALLOCATE( t_snd_tmp(1:nz_snd) )
       t_snd_tmp = 0.0_wp

!
!--    Read temperature from file
       nc_stat = NF90_INQ_VARID( id, "ReferenceTemperature", id_var )
       nc_stat = NF90_GET_VAR( id, id_var, t_snd_tmp(:), start = (/1/),        &
                               count = (/nz_snd/) )
       CALL netcdf_handle_error_rad( 'read_sounding_data', 553 )

!
!--    Calculate start of sounding data
       nz_snd_start = nz_snd + 1
       nz_snd_end   = nz_snd + 1

!
!--    Start filling vertical dimension at 10hPa above the model domain (hyp is
!--    in Pa, hyp_snd in hPa).
       DO  k = 1, nz_snd
          IF ( hyp_snd_tmp(k) < ( hyp(nzt+1) - 1000.0_wp) * 0.01_wp )  THEN
             nz_snd_start = k
             EXIT
          END IF
       END DO

       IF ( nz_snd_start <= nz_snd )  THEN
          nz_snd_end = nz_snd
       END IF


!
!--    Calculate of total grid points for RRTMG calculations
       nzt_rad = nzt + nz_snd_end - nz_snd_start + 1

!
!--    Save data above LES domain in hyp_snd, t_snd
       ALLOCATE( hyp_snd(nzb+1:nzt_rad) )
       ALLOCATE( t_snd(nzb+1:nzt_rad)   )
       hyp_snd = 0.0_wp
       t_snd = 0.0_wp

       hyp_snd(nzt+2:nzt_rad) = hyp_snd_tmp(nz_snd_start+1:nz_snd_end)
       t_snd(nzt+2:nzt_rad)   = t_snd_tmp(nz_snd_start+1:nz_snd_end)

       nc_stat = NF90_CLOSE( id )

!
!--    Calculate pressure levels on zu and zw grid. Sounding data is added at
!--    top of the LES domain. This routine does not consider horizontal or
!--    vertical variability of pressure and temperature
       ALLOCATE ( rrtm_play(0:0,nzb+1:nzt_rad+1)   )
       ALLOCATE ( rrtm_plev(0:0,nzb+1:nzt_rad+2)   )

       t_surface = pt_surface * exner(nzb)
       DO k = nzb+1, nzt+1
          rrtm_play(0,k) = hyp(k) * 0.01_wp
          rrtm_plev(0,k) = barometric_formula(zw(k-1),                         &
                              pt_surface * exner(nzb), &
                              surface_pressure )
       ENDDO

       DO k = nzt+2, nzt_rad
          rrtm_play(0,k) = hyp_snd(k)
          rrtm_plev(0,k) = 0.5_wp * ( rrtm_play(0,k) + rrtm_play(0,k-1) )
       ENDDO
       rrtm_plev(0,nzt_rad+1) = MAX( 0.5 * hyp_snd(nzt_rad),                   &
                                   1.5 * hyp_snd(nzt_rad)                      &
                                 - 0.5 * hyp_snd(nzt_rad-1) )
       rrtm_plev(0,nzt_rad+2)  = MIN( 1.0E-4_wp,                               &
                                      0.25_wp * rrtm_plev(0,nzt_rad+1) )

       rrtm_play(0,nzt_rad+1) = 0.5 * rrtm_plev(0,nzt_rad+1)

!
!--    Calculate temperature/humidity levels at top of the LES domain.
!--    Currently, the temperature is taken from sounding data (might lead to a
!--    temperature jump at interface. To do: Humidity is currently not
!--    calculated above the LES domain.
       ALLOCATE ( rrtm_tlay(0:0,nzb+1:nzt_rad+1)   )
       ALLOCATE ( rrtm_tlev(0:0,nzb+1:nzt_rad+2)   )

       DO k = nzt+8, nzt_rad
          rrtm_tlay(0,k)   = t_snd(k)
       ENDDO
       rrtm_tlay(0,nzt_rad+1) = 2.0_wp * rrtm_tlay(0,nzt_rad)                  &
                                - rrtm_tlay(0,nzt_rad-1)
       DO k = nzt+9, nzt_rad+1
          rrtm_tlev(0,k)   = rrtm_tlay(0,k-1) + (rrtm_tlay(0,k)                &
                             - rrtm_tlay(0,k-1))                               &
                             / ( rrtm_play(0,k) - rrtm_play(0,k-1) )           &
                             * ( rrtm_plev(0,k) - rrtm_play(0,k-1) )
       ENDDO

       rrtm_tlev(0,nzt_rad+2)   = 2.0_wp * rrtm_tlay(0,nzt_rad+1)              &
                                  - rrtm_tlev(0,nzt_rad)
!
!--    Allocate remaining RRTMG arrays
       ALLOCATE ( rrtm_cicewp(0:0,nzb+1:nzt_rad+1) )
       ALLOCATE ( rrtm_cldfr(0:0,nzb+1:nzt_rad+1) )
       ALLOCATE ( rrtm_cliqwp(0:0,nzb+1:nzt_rad+1) )
       ALLOCATE ( rrtm_reice(0:0,nzb+1:nzt_rad+1) )
       ALLOCATE ( rrtm_reliq(0:0,nzb+1:nzt_rad+1) )
       ALLOCATE ( rrtm_lw_taucld(1:nbndlw+1,0:0,nzb+1:nzt_rad+1) )
       ALLOCATE ( rrtm_lw_tauaer(0:0,nzb+1:nzt_rad+1,1:nbndlw+1) )
       ALLOCATE ( rrtm_sw_taucld(1:nbndsw+1,0:0,nzb+1:nzt_rad+1) )
       ALLOCATE ( rrtm_sw_ssacld(1:nbndsw+1,0:0,nzb+1:nzt_rad+1) )
       ALLOCATE ( rrtm_sw_asmcld(1:nbndsw+1,0:0,nzb+1:nzt_rad+1) )
       ALLOCATE ( rrtm_sw_fsfcld(1:nbndsw+1,0:0,nzb+1:nzt_rad+1) )
       ALLOCATE ( rrtm_sw_tauaer(0:0,nzb+1:nzt_rad+1,1:nbndsw+1) )
       ALLOCATE ( rrtm_sw_ssaaer(0:0,nzb+1:nzt_rad+1,1:nbndsw+1) )
       ALLOCATE ( rrtm_sw_asmaer(0:0,nzb+1:nzt_rad+1,1:nbndsw+1) )
       ALLOCATE ( rrtm_sw_ecaer(0:0,nzb+1:nzt_rad+1,1:naerec+1) )

!
!--    The ice phase is currently not considered in PALM
       rrtm_cicewp = 0.0_wp
       rrtm_reice  = 0.0_wp

!
!--    Set other parameters (move to NAMELIST parameters in the future)
       rrtm_lw_tauaer = 0.0_wp
       rrtm_lw_taucld = 0.0_wp
       rrtm_sw_taucld = 0.0_wp
       rrtm_sw_ssacld = 0.0_wp
       rrtm_sw_asmcld = 0.0_wp
       rrtm_sw_fsfcld = 0.0_wp
       rrtm_sw_tauaer = 0.0_wp
       rrtm_sw_ssaaer = 0.0_wp
       rrtm_sw_asmaer = 0.0_wp
       rrtm_sw_ecaer  = 0.0_wp


       ALLOCATE ( rrtm_swdflx(0:0,nzb:nzt_rad+1)  )
       ALLOCATE ( rrtm_swuflx(0:0,nzb:nzt_rad+1)  )
       ALLOCATE ( rrtm_swhr(0:0,nzb+1:nzt_rad+1)  )
       ALLOCATE ( rrtm_swuflxc(0:0,nzb:nzt_rad+1) )
       ALLOCATE ( rrtm_swdflxc(0:0,nzb:nzt_rad+1) )
       ALLOCATE ( rrtm_swhrc(0:0,nzb+1:nzt_rad+1) )
       ALLOCATE ( rrtm_dirdflux(0:0,nzb:nzt_rad+1) )
       ALLOCATE ( rrtm_difdflux(0:0,nzb:nzt_rad+1) )

       rrtm_swdflx  = 0.0_wp
       rrtm_swuflx  = 0.0_wp
       rrtm_swhr    = 0.0_wp
       rrtm_swuflxc = 0.0_wp
       rrtm_swdflxc = 0.0_wp
       rrtm_swhrc   = 0.0_wp
       rrtm_dirdflux = 0.0_wp
       rrtm_difdflux = 0.0_wp

       ALLOCATE ( rrtm_lwdflx(0:0,nzb:nzt_rad+1)  )
       ALLOCATE ( rrtm_lwuflx(0:0,nzb:nzt_rad+1)  )
       ALLOCATE ( rrtm_lwhr(0:0,nzb+1:nzt_rad+1)  )
       ALLOCATE ( rrtm_lwuflxc(0:0,nzb:nzt_rad+1) )
       ALLOCATE ( rrtm_lwdflxc(0:0,nzb:nzt_rad+1) )
       ALLOCATE ( rrtm_lwhrc(0:0,nzb+1:nzt_rad+1) )

       rrtm_lwdflx  = 0.0_wp
       rrtm_lwuflx  = 0.0_wp
       rrtm_lwhr    = 0.0_wp
       rrtm_lwuflxc = 0.0_wp
       rrtm_lwdflxc = 0.0_wp
       rrtm_lwhrc   = 0.0_wp

       ALLOCATE ( rrtm_lwuflx_dt(0:0,nzb:nzt_rad+1) )
       ALLOCATE ( rrtm_lwuflxc_dt(0:0,nzb:nzt_rad+1) )

       rrtm_lwuflx_dt = 0.0_wp
       rrtm_lwuflxc_dt = 0.0_wp

    END SUBROUTINE read_sounding_data


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Read trace gas data from file and convert into trace gas paths / volume
!> mixing ratios. If a user-defined input file is provided it needs to follow
!> the convections used in RRTMG (see respective netCDF files shipped with
!> RRTMG)
!------------------------------------------------------------------------------!
    SUBROUTINE read_trace_gas_data

       USE rrsw_ncpar

       IMPLICIT NONE

       INTEGER(iwp), PARAMETER :: num_trace_gases = 10 !< number of trace gases (absorbers)

       CHARACTER(LEN=5), DIMENSION(num_trace_gases), PARAMETER ::              & !< trace gas names
           trace_names = (/'O3   ', 'CO2  ', 'CH4  ', 'N2O  ', 'O2   ',        &
                           'CFC11', 'CFC12', 'CFC22', 'CCL4 ', 'H2O  '/)

       INTEGER(iwp) :: id,     & !< NetCDF id
                       k,      & !< loop index
                       m,      & !< loop index
                       n,      & !< loop index
                       nabs,   & !< number of absorbers
                       np,     & !< number of pressure levels
                       id_abs, & !< NetCDF id of the respective absorber
                       id_dim, & !< NetCDF id of asborber's dimension
                       id_var    !< NetCDf id ot the absorber

       REAL(wp) :: p_mls_l, &    !< pressure lower limit for interpolation
                   p_mls_u, &    !< pressure upper limit for interpolation
                   p_wgt_l, &    !< pressure weight lower limit for interpolation
                   p_wgt_u, &    !< pressure weight upper limit for interpolation
                   p_mls_m       !< mean pressure between upper and lower limits


       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  p_mls,          & !< pressure levels for the absorbers
                                                 rrtm_play_tmp,  & !< temporary array for pressure zu-levels
                                                 rrtm_plev_tmp,  & !< temporary array for pressure zw-levels
                                                 trace_path_tmp    !< temporary array for storing trace gas path data

       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  trace_mls,      & !< array for storing the absorber amounts
                                                 trace_mls_path, & !< array for storing trace gas path data
                                                 trace_mls_tmp     !< temporary array for storing trace gas data


!
!--    In case of updates, deallocate arrays first (sufficient to check one
!--    array as the others are automatically allocated)
       IF ( ALLOCATED ( rrtm_o3vmr ) )  THEN
          DEALLOCATE ( rrtm_o3vmr  )
          DEALLOCATE ( rrtm_co2vmr )
          DEALLOCATE ( rrtm_ch4vmr )
          DEALLOCATE ( rrtm_n2ovmr )
          DEALLOCATE ( rrtm_o2vmr  )
          DEALLOCATE ( rrtm_cfc11vmr )
          DEALLOCATE ( rrtm_cfc12vmr )
          DEALLOCATE ( rrtm_cfc22vmr )
          DEALLOCATE ( rrtm_ccl4vmr  )
          DEALLOCATE ( rrtm_h2ovmr  )
       ENDIF

!
!--    Allocate trace gas profiles
       ALLOCATE ( rrtm_o3vmr(0:0,1:nzt_rad+1)  )
       ALLOCATE ( rrtm_co2vmr(0:0,1:nzt_rad+1) )
       ALLOCATE ( rrtm_ch4vmr(0:0,1:nzt_rad+1) )
       ALLOCATE ( rrtm_n2ovmr(0:0,1:nzt_rad+1) )
       ALLOCATE ( rrtm_o2vmr(0:0,1:nzt_rad+1)  )
       ALLOCATE ( rrtm_cfc11vmr(0:0,1:nzt_rad+1) )
       ALLOCATE ( rrtm_cfc12vmr(0:0,1:nzt_rad+1) )
       ALLOCATE ( rrtm_cfc22vmr(0:0,1:nzt_rad+1) )
       ALLOCATE ( rrtm_ccl4vmr(0:0,1:nzt_rad+1)  )
       ALLOCATE ( rrtm_h2ovmr(0:0,1:nzt_rad+1)  )

!
!--    Open file for reading
       nc_stat = NF90_OPEN( rrtm_input_file, NF90_NOWRITE, id )
       CALL netcdf_handle_error_rad( 'read_trace_gas_data', 549 )
!
!--    Inquire dimension ids and dimensions
       nc_stat = NF90_INQ_DIMID( id, "Pressure", id_dim )
       CALL netcdf_handle_error_rad( 'read_trace_gas_data', 550 )
       nc_stat = NF90_INQUIRE_DIMENSION( id, id_dim, len = np)
       CALL netcdf_handle_error_rad( 'read_trace_gas_data', 550 )

       nc_stat = NF90_INQ_DIMID( id, "Absorber", id_dim )
       CALL netcdf_handle_error_rad( 'read_trace_gas_data', 550 )
       nc_stat = NF90_INQUIRE_DIMENSION( id, id_dim, len = nabs )
       CALL netcdf_handle_error_rad( 'read_trace_gas_data', 550 )


!
!--    Allocate pressure, and trace gas arrays
       ALLOCATE( p_mls(1:np) )
       ALLOCATE( trace_mls(1:num_trace_gases,1:np) )
       ALLOCATE( trace_mls_tmp(1:nabs,1:np) )


       nc_stat = NF90_INQ_VARID( id, "Pressure", id_var )
       CALL netcdf_handle_error_rad( 'read_trace_gas_data', 550 )
       nc_stat = NF90_GET_VAR( id, id_var, p_mls )
       CALL netcdf_handle_error_rad( 'read_trace_gas_data', 550 )

       nc_stat = NF90_INQ_VARID( id, "AbsorberAmountMLS", id_var )
       CALL netcdf_handle_error_rad( 'read_trace_gas_data', 550 )
       nc_stat = NF90_GET_VAR( id, id_var, trace_mls_tmp )
       CALL netcdf_handle_error_rad( 'read_trace_gas_data', 550 )


!
!--    Write absorber amounts (mls) to trace_mls
       DO n = 1, num_trace_gases
          CALL getAbsorberIndex( TRIM( trace_names(n) ), id_abs )

          trace_mls(n,1:np) = trace_mls_tmp(id_abs,1:np)

!
!--       Replace missing values by zero
          WHERE ( trace_mls(n,:) > 2.0_wp )
             trace_mls(n,:) = 0.0_wp
          END WHERE
       END DO

       DEALLOCATE ( trace_mls_tmp )

       nc_stat = NF90_CLOSE( id )
       CALL netcdf_handle_error_rad( 'read_trace_gas_data', 551 )

!
!--    Add extra pressure level for calculations of the trace gas paths
       ALLOCATE ( rrtm_play_tmp(1:nzt_rad+1) )
       ALLOCATE ( rrtm_plev_tmp(1:nzt_rad+2) )

       rrtm_play_tmp(1:nzt_rad)   = rrtm_play(0,1:nzt_rad)
       rrtm_plev_tmp(1:nzt_rad+1) = rrtm_plev(0,1:nzt_rad+1)
       rrtm_play_tmp(nzt_rad+1)   = rrtm_plev(0,nzt_rad+1) * 0.5_wp
       rrtm_plev_tmp(nzt_rad+2)   = MIN( 1.0E-4_wp, 0.25_wp                    &
                                         * rrtm_plev(0,nzt_rad+1) )

!
!--    Calculate trace gas path (zero at surface) with interpolation to the
!--    sounding levels
       ALLOCATE ( trace_mls_path(1:nzt_rad+2,1:num_trace_gases) )

       trace_mls_path(nzb+1,:) = 0.0_wp

       DO k = nzb+2, nzt_rad+2
          DO m = 1, num_trace_gases
             trace_mls_path(k,m) = trace_mls_path(k-1,m)

!
!--          When the pressure level is higher than the trace gas pressure
!--          level, assume that
             IF ( rrtm_plev_tmp(k-1) > p_mls(1) )  THEN

                trace_mls_path(k,m) = trace_mls_path(k,m) + trace_mls(m,1)     &
                                      * ( rrtm_plev_tmp(k-1)                   &
                                          - MAX( p_mls(1), rrtm_plev_tmp(k) )  &
                                        ) / g
             ENDIF

!
!--          Integrate for each sounding level from the contributing p_mls
!--          levels
             DO n = 2, np
!
!--             Limit p_mls so that it is within the model level
                p_mls_u = MIN( rrtm_plev_tmp(k-1),                             &
                          MAX( rrtm_plev_tmp(k), p_mls(n) ) )
                p_mls_l = MIN( rrtm_plev_tmp(k-1),                             &
                          MAX( rrtm_plev_tmp(k), p_mls(n-1) ) )

                IF ( p_mls_l > p_mls_u )  THEN

!
!--                Calculate weights for interpolation
                   p_mls_m = 0.5_wp * (p_mls_l + p_mls_u)
                   p_wgt_u = (p_mls(n-1) - p_mls_m) / (p_mls(n-1) - p_mls(n))
                   p_wgt_l = (p_mls_m - p_mls(n))   / (p_mls(n-1) - p_mls(n))

!
!--                Add level to trace gas path
                   trace_mls_path(k,m) = trace_mls_path(k,m)                   &
                                         +  ( p_wgt_u * trace_mls(m,n)         &
                                            + p_wgt_l * trace_mls(m,n-1) )     &
                                         * (p_mls_l - p_mls_u) / g
                ENDIF
             ENDDO

             IF ( rrtm_plev_tmp(k) < p_mls(np) )  THEN
                trace_mls_path(k,m) = trace_mls_path(k,m) + trace_mls(m,np)    &
                                      * ( MIN( rrtm_plev_tmp(k-1), p_mls(np) ) &
                                          - rrtm_plev_tmp(k)                   &
                                        ) / g
             ENDIF
          ENDDO
       ENDDO


!
!--    Prepare trace gas path profiles
       ALLOCATE ( trace_path_tmp(1:nzt_rad+1) )

       DO m = 1, num_trace_gases

          trace_path_tmp(1:nzt_rad+1) = ( trace_mls_path(2:nzt_rad+2,m)        &
                                       - trace_mls_path(1:nzt_rad+1,m) ) * g   &
                                       / ( rrtm_plev_tmp(1:nzt_rad+1)          &
                                       - rrtm_plev_tmp(2:nzt_rad+2) )

!
!--       Save trace gas paths to the respective arrays
          SELECT CASE ( TRIM( trace_names(m) ) )

             CASE ( 'O3' )

                rrtm_o3vmr(0,:) = trace_path_tmp(:)

             CASE ( 'CO2' )

                rrtm_co2vmr(0,:) = trace_path_tmp(:)

             CASE ( 'CH4' )

                rrtm_ch4vmr(0,:) = trace_path_tmp(:)

             CASE ( 'N2O' )

                rrtm_n2ovmr(0,:) = trace_path_tmp(:)

             CASE ( 'O2' )

                rrtm_o2vmr(0,:) = trace_path_tmp(:)

             CASE ( 'CFC11' )

                rrtm_cfc11vmr(0,:) = trace_path_tmp(:)

             CASE ( 'CFC12' )

                rrtm_cfc12vmr(0,:) = trace_path_tmp(:)

             CASE ( 'CFC22' )

                rrtm_cfc22vmr(0,:) = trace_path_tmp(:)

             CASE ( 'CCL4' )

                rrtm_ccl4vmr(0,:) = trace_path_tmp(:)

             CASE ( 'H2O' )

                rrtm_h2ovmr(0,:) = trace_path_tmp(:)

             CASE DEFAULT

          END SELECT

       ENDDO

       DEALLOCATE ( trace_path_tmp )
       DEALLOCATE ( trace_mls_path )
       DEALLOCATE ( rrtm_play_tmp )
       DEALLOCATE ( rrtm_plev_tmp )
       DEALLOCATE ( trace_mls )
       DEALLOCATE ( p_mls )

    END SUBROUTINE read_trace_gas_data


    SUBROUTINE netcdf_handle_error_rad( routine_name, errno )

       USE control_parameters,                                                 &
           ONLY:  message_string

       USE NETCDF

       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=6) ::  message_identifier
       CHARACTER(LEN=*) ::  routine_name

       INTEGER(iwp) ::  errno

       IF ( nc_stat /= NF90_NOERR )  THEN

          WRITE( message_identifier, '(''NC'',I4.4)' )  errno
          message_string = TRIM( NF90_STRERROR( nc_stat ) )

          CALL message( routine_name, message_identifier, 2, 2, 0, 6, 1 )

       ENDIF

    END SUBROUTINE netcdf_handle_error_rad
#endif


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate temperature tendency due to radiative cooling/heating.
!> Cache-optimized version.
!------------------------------------------------------------------------------!
#if defined( __rrtmg )
 SUBROUTINE radiation_tendency_ij ( i, j, tend )

    IMPLICIT NONE

    INTEGER(iwp) :: i, j, k !< loop indices

    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) :: tend !< pt tendency term

    IF ( radiation_scheme == 'rrtmg' )  THEN
!
!--    Calculate tendency based on heating rate
       DO k = nzb+1, nzt+1
          tend(k,j,i) = tend(k,j,i) + (rad_lw_hr(k,j,i) + rad_sw_hr(k,j,i))    &
                                         * d_exner(k) * d_seconds_hour
       ENDDO

    ENDIF

 END SUBROUTINE radiation_tendency_ij
#endif


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate temperature tendency due to radiative cooling/heating.
!> Vector-optimized version
!------------------------------------------------------------------------------!
#if defined( __rrtmg )
 SUBROUTINE radiation_tendency ( tend )

    USE indices,                                                               &
        ONLY:  nxl, nxr, nyn, nys

    IMPLICIT NONE

    INTEGER(iwp) :: i, j, k !< loop indices

    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) :: tend !< pt tendency term

    IF ( radiation_scheme == 'rrtmg' )  THEN
!
!--    Calculate tendency based on heating rate
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO k = nzb+1, nzt+1
                tend(k,j,i) = tend(k,j,i) + ( rad_lw_hr(k,j,i)                 &
                                          +  rad_sw_hr(k,j,i) ) * d_exner(k)   &
                                          * d_seconds_hour
             ENDDO
          ENDDO
       ENDDO
    ENDIF

 END SUBROUTINE radiation_tendency
#endif

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Radiative Transfer Model (RTM) version 3.0 for modelling of radiation
!> interactions within urban canopy or inside of surface layer in complex terrain.
!> This subroutine calculates interaction of the solar SW and LW radiation
!> with urban and land surfaces and updates all surface heatfluxes.
!> It also calculates interactions of SW and LW radiation with resolved
!> plant canopy and calculates the corresponding plant canopy heat fluxes.
!> The subroutine also models spatial and temporal distribution of Mean
!> Radiant Temperature (MRT). The resulting values are provided to other
!> PALM-4U modules (RRTMG, USM, LSM, PCM and BIO).
!>
!> The new version 3.0 was radically rewriten from version 1.0.
!> The most significant changes include new angular discretization scheme,
!> redesigned and significantly optimized raytracing scheme, new processes
!> included in modelling (e.g. intetrations of LW radiation with PC),
!> integrated calculation of Mean Radiant Temperature (MRT), and improved
!> and enhanced output and debug capabilities. This new version significantly
!> improves effectivity of the paralelization and the scalability of the model
!> and allows simulation of extensive domain with appropriate HPC resources.
!>
!> More info about RTM v.1.0. see:
!> Resler et al., GMD. 2017, https://doi.org/10.5194/gmd-10-3635-2017
!> Info about RTM v. 3.0 see:
!> Krc et al. 2020 (to appear in GMD),
!> Maronga et al., GMDD 2019,  https://doi.org/10.5194/gmd-2019-103
!>


!------------------------------------------------------------------------------!

 SUBROUTINE radiation_interaction

    USE control_parameters,                                                    &
        ONLY:  rotation_angle

     IMPLICIT NONE

     INTEGER(iwp)                      ::  i, j, k, kk, d, refstep, m, mm, l, ll
     INTEGER(iwp)                      ::  isurf, isurfsrc, isvf, icsf, ipcgb
     INTEGER(iwp)                      ::  imrt, imrtf
     INTEGER(iwp)                      ::  isd                !< solar direction number
     INTEGER(iwp)                      ::  pc_box_dimshift    !< transform for best accuracy
     INTEGER(iwp), DIMENSION(0:3)      ::  reorder = (/ 1, 0, 3, 2 /)

     REAL(wp), DIMENSION(3,3)          ::  mrot               !< grid rotation matrix (zyx)
     REAL(wp), DIMENSION(3,0:nsurf_type)::  vnorm             !< face direction normal vectors (zyx)
     REAL(wp), DIMENSION(3)            ::  sunorig            !< grid rotated solar direction unit vector (zyx)
     REAL(wp), DIMENSION(3)            ::  sunorig_grid       !< grid squashed solar direction unit vector (zyx)
     REAL(wp), DIMENSION(0:nsurf_type) ::  costheta           !< direct irradiance factor of solar angle
     REAL(wp), DIMENSION(nz_urban_b:nz_urban_t) ::  pchf_prep          !< precalculated factor for canopy temperature tendency
     REAL(wp)                          ::  pc_box_area, pc_abs_frac, pc_abs_eff
     REAL(wp)                          ::  asrc               !< area of source face
     REAL(wp)                          ::  pcrad              !< irradiance from plant canopy
     REAL(wp)                          ::  pabsswl  = 0.0_wp  !< total absorbed SW radiation energy in local processor (W)
     REAL(wp)                          ::  pabssw   = 0.0_wp  !< total absorbed SW radiation energy in all processors (W)
     REAL(wp)                          ::  pabslwl  = 0.0_wp  !< total absorbed LW radiation energy in local processor (W)
     REAL(wp)                          ::  pabslw   = 0.0_wp  !< total absorbed LW radiation energy in all processors (W)
     REAL(wp)                          ::  pemitlwl = 0.0_wp  !< total emitted LW radiation energy in all processors (W)
     REAL(wp)                          ::  pemitlw  = 0.0_wp  !< total emitted LW radiation energy in all processors (W)
     REAL(wp)                          ::  pinswl   = 0.0_wp  !< total received SW radiation energy in local processor (W)
     REAL(wp)                          ::  pinsw    = 0.0_wp  !< total received SW radiation energy in all processor (W)
     REAL(wp)                          ::  pinlwl   = 0.0_wp  !< total received LW radiation energy in local processor (W)
     REAL(wp)                          ::  pinlw    = 0.0_wp  !< total received LW radiation energy in all processor (W)
     REAL(wp)                          ::  emiss_sum_surfl    !< sum of emissisivity of surfaces in local processor
     REAL(wp)                          ::  emiss_sum_surf     !< sum of emissisivity of surfaces in all processor
     REAL(wp)                          ::  area_surfl         !< total area of surfaces in local processor
     REAL(wp)                          ::  area_surf          !< total area of surfaces in all processor
     REAL(wp)                          ::  area_hor           !< total horizontal area of domain in all processor
#if defined( __parallel )
     REAL(wp), DIMENSION(1:7)          ::  combine_allreduce   !< dummy array used to combine several MPI_ALLREDUCE calls
     REAL(wp), DIMENSION(1:7)          ::  combine_allreduce_l !< dummy array used to combine several MPI_ALLREDUCE calls
#endif

     IF ( debug_output_timestep )  CALL debug_message( 'radiation_interaction', 'start' )

     IF ( plant_canopy )  THEN
         pchf_prep(:) = r_d * exner(nz_urban_b:nz_urban_t)                                 &
                     / (c_p * hyp(nz_urban_b:nz_urban_t) * dx*dy*dz(1)) !< equals to 1 / (rho * c_p * Vbox * T)
     ENDIF

     sun_direction = .TRUE.
     CALL get_date_time( time_since_reference_point, &
                         day_of_year=day_of_year,    &
                         second_of_day=second_of_day )
     CALL calc_zenith( day_of_year, second_of_day ) !< required also for diffusion radiation

!
!--     prepare rotated normal vectors and irradiance factor
     vnorm(1,:) = kdir(:)
     vnorm(2,:) = jdir(:)
     vnorm(3,:) = idir(:)
     mrot(1, :) = (/ 1._wp,  0._wp,      0._wp      /)
     mrot(2, :) = (/ 0._wp,  COS(rotation_angle), SIN(rotation_angle) /)
     mrot(3, :) = (/ 0._wp, -SIN(rotation_angle), COS(rotation_angle) /)
     sunorig = (/ cos_zenith, sun_dir_lat, sun_dir_lon /)
     sunorig = MATMUL(mrot, sunorig)
     DO d = 0, nsurf_type
         costheta(d) = DOT_PRODUCT(sunorig, vnorm(:,d))
     ENDDO

     IF ( cos_zenith > 0 )  THEN
!--      now we will "squash" the sunorig vector by grid box size in
!--      each dimension, so that this new direction vector will allow us
!--      to traverse the ray path within grid coordinates directly
         sunorig_grid = (/ sunorig(1)/dz(1), sunorig(2)/dy, sunorig(3)/dx /)
!--      sunorig_grid = sunorig_grid / norm2(sunorig_grid)
         sunorig_grid = sunorig_grid / SQRT(SUM(sunorig_grid**2))

         IF ( npcbl > 0 )  THEN
!--         precompute effective box depth with prototype Leaf Area Density
            pc_box_dimshift = MAXLOC(ABS(sunorig), 1) - 1
            CALL box_absorb(CSHIFT((/dz(1),dy,dx/), pc_box_dimshift),       &
                                60, prototype_lad,                          &
                                CSHIFT(ABS(sunorig), pc_box_dimshift),      &
                                pc_box_area, pc_abs_frac)
            pc_box_area = pc_box_area * ABS(sunorig(pc_box_dimshift+1)      &
                          / sunorig(1))
            pc_abs_eff = LOG(1._wp - pc_abs_frac) / prototype_lad
         ENDIF
     ENDIF
!
!--  Split downwelling shortwave radiation into a diffuse and a direct part.
!--  Note, if radiation scheme is RRTMG or diffuse radiation is externally
!--  prescribed, this is not required. Please note, in case of external
!--  radiation, the clear-sky model is applied during spinup, so that
!--  radiation need to be split also in this case.
     IF ( radiation_scheme == 'constant'   .OR.                                &
          radiation_scheme == 'clear-sky'  .OR.                                &
          ( radiation_scheme == 'external'  .AND.                              &
            .NOT. rad_sw_in_dif_f%from_file  )  .OR.                           &
          ( radiation_scheme == 'external'  .AND.                              &
            time_since_reference_point < 0.0_wp ) )  THEN
        CALL calc_diffusion_radiation
     ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--     First pass: direct + diffuse irradiance + thermal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     surfinswdir   = 0._wp !nsurfl
     surfins       = 0._wp !nsurfl
     surfinl       = 0._wp !nsurfl
     surfoutsl(:)  = 0.0_wp !start-end
     surfoutll(:)  = 0.0_wp !start-end
     IF ( nmrtbl > 0 )  THEN
        mrtinsw(:) = 0._wp
        mrtinlw(:) = 0._wp
     ENDIF
     surfinlg(:)  = 0._wp !global


!--  Set up thermal radiation from surfaces
!--  emiss_surf is defined only for surfaces for which energy balance is calculated
!--  Workaround: reorder surface data type back on 1D array including all surfaces,
!--  which implies to reorder horizontal and vertical surfaces
!
!--  Horizontal walls
     mm = 1
     DO  i = nxl, nxr
        DO  j = nys, nyn
!
!--        urban
           DO  m = surf_usm_h%start_index(j,i), surf_usm_h%end_index(j,i)
              surfoutll(mm) = SUM ( surf_usm_h%frac(m,:) *                  &
                                    surf_usm_h%emissivity(m,:) )            &
                                  * sigma_sb                                &
                                  * surf_usm_h%pt_surface(m)**4
              albedo_surf(mm) = SUM ( surf_usm_h%frac(m,:) *                &
                                      surf_usm_h%albedo(m,:) )
              emiss_surf(mm)  = SUM ( surf_usm_h%frac(m,:) *                &
                                      surf_usm_h%emissivity(m,:) )
              mm = mm + 1
           ENDDO
!
!--        land
           DO  m = surf_lsm_h%start_index(j,i), surf_lsm_h%end_index(j,i)
              surfoutll(mm) = SUM ( surf_lsm_h%frac(m,:) *                  &
                                    surf_lsm_h%emissivity(m,:) )            &
                                  * sigma_sb                                &
                                  * surf_lsm_h%pt_surface(m)**4
              albedo_surf(mm) = SUM ( surf_lsm_h%frac(m,:) *                &
                                      surf_lsm_h%albedo(m,:) )
              emiss_surf(mm)  = SUM ( surf_lsm_h%frac(m,:) *                &
                                      surf_lsm_h%emissivity(m,:) )
              mm = mm + 1
           ENDDO
        ENDDO
     ENDDO
!
!--  Vertical walls
     DO  i = nxl, nxr
        DO  j = nys, nyn
           DO  ll = 0, 3
              l = reorder(ll)
!
!--           urban
              DO  m = surf_usm_v(l)%start_index(j,i),                       &
                      surf_usm_v(l)%end_index(j,i)
                 surfoutll(mm) = SUM ( surf_usm_v(l)%frac(m,:) *            &
                                       surf_usm_v(l)%emissivity(m,:) )      &
                                  * sigma_sb                                &
                                  * surf_usm_v(l)%pt_surface(m)**4
                 albedo_surf(mm) = SUM ( surf_usm_v(l)%frac(m,:) *          &
                                         surf_usm_v(l)%albedo(m,:) )
                 emiss_surf(mm)  = SUM ( surf_usm_v(l)%frac(m,:) *          &
                                         surf_usm_v(l)%emissivity(m,:) )
                 mm = mm + 1
              ENDDO
!
!--           land
              DO  m = surf_lsm_v(l)%start_index(j,i),                       &
                      surf_lsm_v(l)%end_index(j,i)
                 surfoutll(mm) = SUM ( surf_lsm_v(l)%frac(m,:) *            &
                                       surf_lsm_v(l)%emissivity(m,:) )      &
                                  * sigma_sb                                &
                                  * surf_lsm_v(l)%pt_surface(m)**4
                 albedo_surf(mm) = SUM ( surf_lsm_v(l)%frac(m,:) *          &
                                         surf_lsm_v(l)%albedo(m,:) )
                 emiss_surf(mm)  = SUM ( surf_lsm_v(l)%frac(m,:) *          &
                                         surf_lsm_v(l)%emissivity(m,:) )
                 mm = mm + 1
              ENDDO
           ENDDO
        ENDDO
     ENDDO

     IF ( trace_fluxes_above >= 0.0_wp )  THEN
        CALL radiation_print_debug_surf( 'surfoutll before initial pass', surfoutll )
        CALL radiation_print_debug_horz( 'rad_lw_in_diff before initial pass', rad_lw_in_diff )
        CALL radiation_print_debug_horz( 'rad_sw_in_diff before initial pass', rad_sw_in_diff )
        CALL radiation_print_debug_horz( 'rad_sw_in_dir before initial pass', rad_sw_in_dir )
     ENDIF

#if defined( __parallel )
!--     might be optimized and gather only values relevant for current processor
     CALL MPI_AllGatherv(surfoutll, nsurfl, MPI_REAL, &
                         surfoutl, nsurfs, surfstart, MPI_REAL, comm2d, ierr) !nsurf global
     IF ( ierr /= 0 ) THEN
         WRITE(9,*) 'Error MPI_AllGatherv1:', ierr, SIZE(surfoutll), nsurfl, &
                     SIZE(surfoutl), nsurfs, surfstart
         FLUSH(9)
     ENDIF
#else
     surfoutl(:) = surfoutll(:) !nsurf global
#endif

     IF ( surface_reflections)  THEN
        DO  isvf = 1, nsvfl
           isurf = svfsurf(1, isvf)
           k     = surfl(iz, isurf)
           j     = surfl(iy, isurf)
           i     = surfl(ix, isurf)
           isurfsrc = svfsurf(2, isvf)
!
!--        For surface-to-surface factors we calculate thermal radiation in 1st pass
           IF ( plant_lw_interact )  THEN
              surfinl(isurf) = surfinl(isurf) + svf(1,isvf) * svf(2,isvf) * surfoutl(isurfsrc)
           ELSE
              surfinl(isurf) = surfinl(isurf) + svf(1,isvf) * surfoutl(isurfsrc)
           ENDIF
        ENDDO
     ENDIF
!
!--  diffuse radiation using sky view factor
     DO isurf = 1, nsurfl
        j = surfl(iy, isurf)
        i = surfl(ix, isurf)
        surfinswdif(isurf) = rad_sw_in_diff(j,i) * skyvft(isurf)
        IF ( plant_lw_interact )  THEN
           surfinlwdif(isurf) = rad_lw_in_diff(j,i) * skyvft(isurf)
        ELSE
           surfinlwdif(isurf) = rad_lw_in_diff(j,i) * skyvf(isurf)
        ENDIF
     ENDDO
!
!--  MRT diffuse irradiance
     DO  imrt = 1, nmrtbl
        j = mrtbl(iy, imrt)
        i = mrtbl(ix, imrt)
        mrtinsw(imrt) = mrtskyt(imrt) * rad_sw_in_diff(j,i)
        mrtinlw(imrt) = mrtsky(imrt) * rad_lw_in_diff(j,i)
     ENDDO

     !-- direct radiation
     IF ( cos_zenith > 0 )  THEN
        !--Identify solar direction vector (discretized number) 1)
        !--
        j = FLOOR(ACOS(cos_zenith) / pi * raytrace_discrete_elevs)
        i = MODULO(NINT(ATAN2(sun_dir_lon, sun_dir_lat)               &
                        / (2.0_wp*pi) * raytrace_discrete_azims-0.5_wp, iwp), &
                   raytrace_discrete_azims)
        isd = dsidir_rev(j, i)
!-- TODO: check if isd = -1 to report that this solar position is not precalculated
        DO isurf = 1, nsurfl
           j = surfl(iy, isurf)
           i = surfl(ix, isurf)
           surfinswdir(isurf) = rad_sw_in_dir(j,i) *                        &
                costheta(surfl(id, isurf)) * dsitrans(isurf, isd) / cos_zenith
        ENDDO
!
!--     MRT direct irradiance
        DO  imrt = 1, nmrtbl
           j = mrtbl(iy, imrt)
           i = mrtbl(ix, imrt)
           mrtinsw(imrt) = mrtinsw(imrt) + mrtdsit(imrt, isd) * rad_sw_in_dir(j,i) &
                                     / cos_zenith / 4.0_wp ! normal to sphere
        ENDDO
     ENDIF
!
!--  MRT first pass thermal
     DO  imrtf = 1, nmrtf
        imrt = mrtfsurf(1, imrtf)
        isurfsrc = mrtfsurf(2, imrtf)
        mrtinlw(imrt) = mrtinlw(imrt) + mrtf(imrtf) * surfoutl(isurfsrc)
     ENDDO
!
!--  Absorption in each local plant canopy grid box from the first atmospheric
!--  pass of radiation
     IF ( npcbl > 0 )  THEN

         pcbinswdir(:) = 0.0_wp
         pcbinswdif(:) = 0.0_wp
         pcbinlw(:) = 0.0_wp

         DO icsf = 1, ncsfl
             ipcgb = csfsurf(1, icsf)
             i = pcbl(ix,ipcgb)
             j = pcbl(iy,ipcgb)
             k = pcbl(iz,ipcgb)
             isurfsrc = csfsurf(2, icsf)

             IF ( isurfsrc == -1 )  THEN
!
!--             Diffuse radiation from sky
                pcbinswdif(ipcgb) = csf(1,icsf) * rad_sw_in_diff(j,i)
!
!--             Absorbed diffuse LW radiation from sky minus emitted to sky
                IF ( plant_lw_interact )  THEN
                   pcbinlw(ipcgb) = csf(1,icsf)                                  &
                                       * (rad_lw_in_diff(j, i)                   &
                                          - sigma_sb * (pt(k,j,i)*exner(k))**4)
                ENDIF
!
!--             Direct solar radiation
                IF ( cos_zenith > 0 )  THEN
!--                Estimate directed box absorption
                   pc_abs_frac = 1.0_wp - exp(pc_abs_eff * lad_s(k,j,i))
!
!--                isd has already been established, see 1)
                   pcbinswdir(ipcgb) = rad_sw_in_dir(j, i) * pc_box_area &
                                       * pc_abs_frac * dsitransc(ipcgb, isd)
                ENDIF
             ELSE
                IF ( plant_lw_interact )  THEN
!
!--                Thermal emission from plan canopy towards respective face
                   pcrad = sigma_sb * (pt(k,j,i) * exner(k))**4 * csf(1,icsf)
                   surfinlg(isurfsrc) = surfinlg(isurfsrc) + pcrad
!
!--                Remove the flux above + absorb LW from first pass from surfaces
                   asrc = facearea(surf(id, isurfsrc))
                   pcbinlw(ipcgb) = pcbinlw(ipcgb)                      &
                                    + (csf(1,icsf) * surfoutl(isurfsrc) & ! Absorb from first pass surf emit
                                       - pcrad)                         & ! Remove emitted heatflux
                                    * asrc
                ENDIF
             ENDIF
         ENDDO

         pcbinsw(:) = pcbinswdir(:) + pcbinswdif(:)
     ENDIF

     IF ( trace_fluxes_above >= 0.0_wp )  THEN
        CALL radiation_print_debug_surf( 'surfinl after initial pass', surfinl )
        CALL radiation_print_debug_surf( 'surfinlwdif after initial pass', surfinlwdif )
        CALL radiation_print_debug_surf( 'surfinswdif after initial pass', surfinswdif )
        CALL radiation_print_debug_surf( 'surfinswdir after initial pass', surfinswdir )
        IF ( npcbl > 0 )  THEN
           CALL radiation_print_debug_pcb( 'pcbinlw after initial pass', pcbinlw )
           CALL radiation_print_debug_pcb( 'pcbinswdif after initial pass', pcbinswdif )
           CALL radiation_print_debug_pcb( 'pcbinswdir after initial pass', pcbinswdir )
        ENDIF
     ENDIF

     IF ( plant_lw_interact )  THEN
!
!--     Exchange incoming lw radiation from plant canopy
#if defined( __parallel )
        CALL MPI_Allreduce(MPI_IN_PLACE, surfinlg, nsurf, MPI_REAL, MPI_SUM, comm2d, ierr)
        IF ( ierr /= 0 )  THEN
           WRITE (9,*) 'Error MPI_Allreduce:', ierr
           FLUSH(9)
        ENDIF
        surfinl(:) = surfinl(:) + surfinlg(surfstart(myid)+1:surfstart(myid+1))
#else
        surfinl(:) = surfinl(:) + surfinlg(:)
#endif
     ENDIF

     IF ( trace_fluxes_above >= 0.0_wp )  THEN
        CALL radiation_print_debug_surf( 'surfinl after PC emiss', surfinl )
     ENDIF

     surfins = surfinswdir + surfinswdif
     surfinl = surfinl + surfinlwdif
     surfinsw = surfins
     surfinlw = surfinl
     surfoutsw = 0.0_wp
     surfoutlw = surfoutll
     surfemitlwl = surfoutll

     IF ( .NOT.  surface_reflections )  THEN
!
!--     Set nrefsteps to 0 to disable reflections
        nrefsteps = 0
        surfoutsl = albedo_surf * surfins
        surfoutll = (1.0_wp - emiss_surf) * surfinl
        surfoutsw = surfoutsw + surfoutsl
        surfoutlw = surfoutlw + surfoutll
     ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--     Next passes - reflections
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     DO refstep = 1, nrefsteps

         surfoutsl = albedo_surf * surfins
!
!--      for non-transparent surfaces, longwave albedo is 1 - emissivity
         surfoutll = (1.0_wp - emiss_surf) * surfinl

         IF ( trace_fluxes_above >= 0.0_wp )  THEN
            CALL radiation_print_debug_surf( 'surfoutll before reflective pass', surfoutll, refstep )
            CALL radiation_print_debug_surf( 'surfoutsl before reflective pass', surfoutsl, refstep )
         ENDIF

#if defined( __parallel )
         CALL MPI_AllGatherv(surfoutsl, nsurfl, MPI_REAL, &
             surfouts, nsurfs, surfstart, MPI_REAL, comm2d, ierr)
         IF ( ierr /= 0 )  THEN
             WRITE(9,*) 'Error MPI_AllGatherv2:', ierr, SIZE(surfoutsl), nsurfl, &
                        SIZE(surfouts), nsurfs, surfstart
             FLUSH(9)
         ENDIF

         CALL MPI_AllGatherv(surfoutll, nsurfl, MPI_REAL, &
             surfoutl, nsurfs, surfstart, MPI_REAL, comm2d, ierr)
         IF ( ierr /= 0 )  THEN
             WRITE(9,*) 'Error MPI_AllGatherv3:', ierr, SIZE(surfoutll), nsurfl, &
                        SIZE(surfoutl), nsurfs, surfstart
             FLUSH(9)
         ENDIF

#else
         surfouts = surfoutsl
         surfoutl = surfoutll
#endif
!
!--      Reset for the input from next reflective pass
         surfins = 0.0_wp
         surfinl = 0.00_wp
!
!--      Reflected radiation
         DO isvf = 1, nsvfl
             isurf = svfsurf(1, isvf)
             isurfsrc = svfsurf(2, isvf)
             surfins(isurf) = surfins(isurf) + svf(1,isvf) * svf(2,isvf) * surfouts(isurfsrc)
             IF ( plant_lw_interact )  THEN
                surfinl(isurf) = surfinl(isurf) + svf(1,isvf) * svf(2,isvf) * surfoutl(isurfsrc)
             ELSE
                surfinl(isurf) = surfinl(isurf) + svf(1,isvf) * surfoutl(isurfsrc)
             ENDIF
         ENDDO
!
!--      NOTE: PC absorbtion and MRT from reflected can both be done at once
!--      after all reflections if we do one more MPI_ALLGATHERV on surfout.
!--      Advantage: less local computation. Disadvantage: one more collective
!--      MPI call.
!
!--      Radiation absorbed by plant canopy
         DO  icsf = 1, ncsfl
             ipcgb = csfsurf(1, icsf)
             isurfsrc = csfsurf(2, icsf)
             IF ( isurfsrc == -1 )  CYCLE ! sky->face only in 1st pass, not here
!
!--          Calculate source surface area. If the `surf' array is removed
!--          before timestepping starts (future version), then asrc must be
!--          stored within `csf'
             asrc = facearea(surf(id, isurfsrc))
             pcbinsw(ipcgb) = pcbinsw(ipcgb) + csf(1,icsf) * surfouts(isurfsrc) * asrc
             IF ( plant_lw_interact )  THEN
                pcbinlw(ipcgb) = pcbinlw(ipcgb) + csf(1,icsf) * surfoutl(isurfsrc) * asrc
             ENDIF
         ENDDO
!
!--      MRT reflected
         DO  imrtf = 1, nmrtf
            imrt = mrtfsurf(1, imrtf)
            isurfsrc = mrtfsurf(2, imrtf)
            mrtinsw(imrt) = mrtinsw(imrt) + mrtft(imrtf) * surfouts(isurfsrc)
            mrtinlw(imrt) = mrtinlw(imrt) + mrtf(imrtf) * surfoutl(isurfsrc)
         ENDDO

         IF ( trace_fluxes_above >= 0.0_wp )  THEN
            CALL radiation_print_debug_surf( 'surfinl after reflected pass', surfinl, refstep )
            CALL radiation_print_debug_surf( 'surfins after reflected pass', surfins, refstep )
            CALL radiation_print_debug_pcb( 'pcbinlw after reflected pass', pcbinlw, refstep )
            CALL radiation_print_debug_pcb( 'pcbinsw after reflected pass', pcbinsw, refstep )
         ENDIF

         surfinsw = surfinsw  + surfins
         surfinlw = surfinlw  + surfinl
         surfoutsw = surfoutsw + surfoutsl
         surfoutlw = surfoutlw + surfoutll

     ENDDO ! refstep

!--  push heat flux absorbed by plant canopy to respective 3D arrays
     IF ( npcbl > 0 )  THEN
         pcm_heating_rate(:,:,:) = 0.0_wp
         DO ipcgb = 1, npcbl
             j = pcbl(iy, ipcgb)
             i = pcbl(ix, ipcgb)
             k = pcbl(iz, ipcgb)
!
!--          Following expression equals former kk = k - nzb_s_inner(j,i)
             kk = k - topo_top_ind(j,i,0)  !- lad arrays are defined flat
             pcm_heating_rate(kk, j, i) = (pcbinsw(ipcgb) + pcbinlw(ipcgb)) &
                 * pchf_prep(k) * pt(k, j, i) !-- = dT/dt
         ENDDO

         IF ( humidity .AND. plant_canopy_transpiration ) THEN
!--          Calculation of plant canopy transpiration rate and correspondidng latent heat rate
             pcm_transpiration_rate(:,:,:) = 0.0_wp
             pcm_latent_rate(:,:,:) = 0.0_wp
             DO ipcgb = 1, npcbl
                 i = pcbl(ix, ipcgb)
                 j = pcbl(iy, ipcgb)
                 k = pcbl(iz, ipcgb)
                 kk = k - topo_top_ind(j,i,0)  !- lad arrays are defined flat
                 CALL pcm_calc_transpiration_rate( i, j, k, kk, pcbinsw(ipcgb), pcbinlw(ipcgb), &
                                                   pcm_transpiration_rate(kk,j,i), pcm_latent_rate(kk,j,i) )
              ENDDO
         ENDIF
     ENDIF
!
!--  Calculate black body MRT (after all reflections)
     IF ( nmrtbl > 0 )  THEN
        IF ( mrt_include_sw )  THEN
           mrt(:) = ((mrtinsw(:) + mrtinlw(:)) / sigma_sb) ** 0.25_wp
        ELSE
           mrt(:) = (mrtinlw(:) / sigma_sb) ** 0.25_wp
        ENDIF
     ENDIF
!
!--     Transfer radiation arrays required for energy balance to the respective data types
     DO  i = 1, nsurfl
        m  = surfl(im,i)
!
!--     (1) Urban surfaces
!--     upward-facing
        IF ( surfl(1,i) == iup_u )  THEN
           surf_usm_h%rad_sw_in(m)  = surfinsw(i)
           surf_usm_h%rad_sw_out(m) = surfoutsw(i)
           surf_usm_h%rad_sw_dir(m) = surfinswdir(i)
           surf_usm_h%rad_sw_dif(m) = surfinswdif(i)
           surf_usm_h%rad_sw_ref(m) = surfinsw(i) - surfinswdir(i) -        &
                                      surfinswdif(i)
           surf_usm_h%rad_sw_res(m) = surfins(i)
           surf_usm_h%rad_lw_in(m)  = surfinlw(i)
           surf_usm_h%rad_lw_out(m) = surfoutlw(i)
           surf_usm_h%rad_net(m)    = surfinsw(i) - surfoutsw(i) +          &
                                      surfinlw(i) - surfoutlw(i)
           surf_usm_h%rad_net_l(m)  = surf_usm_h%rad_net(m)
           surf_usm_h%rad_lw_dif(m) = surfinlwdif(i)
           surf_usm_h%rad_lw_ref(m) = surfinlw(i) - surfinlwdif(i)
           surf_usm_h%rad_lw_res(m) = surfinl(i)
!
!--     northward-facding
        ELSEIF ( surfl(1,i) == inorth_u )  THEN
           surf_usm_v(0)%rad_sw_in(m)  = surfinsw(i)
           surf_usm_v(0)%rad_sw_out(m) = surfoutsw(i)
           surf_usm_v(0)%rad_sw_dir(m) = surfinswdir(i)
           surf_usm_v(0)%rad_sw_dif(m) = surfinswdif(i)
           surf_usm_v(0)%rad_sw_ref(m) = surfinsw(i) - surfinswdir(i) -     &
                                         surfinswdif(i)
           surf_usm_v(0)%rad_sw_res(m) = surfins(i)
           surf_usm_v(0)%rad_lw_in(m)  = surfinlw(i)
           surf_usm_v(0)%rad_lw_out(m) = surfoutlw(i)
           surf_usm_v(0)%rad_net(m)    = surfinsw(i) - surfoutsw(i) +       &
                                         surfinlw(i) - surfoutlw(i)
           surf_usm_v(0)%rad_net_l(m)  = surf_usm_v(0)%rad_net(m)
           surf_usm_v(0)%rad_lw_dif(m) = surfinlwdif(i)
           surf_usm_v(0)%rad_lw_ref(m) = surfinlw(i) - surfinlwdif(i)
           surf_usm_v(0)%rad_lw_res(m) = surfinl(i)
!
!--     southward-facding
        ELSEIF ( surfl(1,i) == isouth_u )  THEN
           surf_usm_v(1)%rad_sw_in(m)  = surfinsw(i)
           surf_usm_v(1)%rad_sw_out(m) = surfoutsw(i)
           surf_usm_v(1)%rad_sw_dir(m) = surfinswdir(i)
           surf_usm_v(1)%rad_sw_dif(m) = surfinswdif(i)
           surf_usm_v(1)%rad_sw_ref(m) = surfinsw(i) - surfinswdir(i) -     &
                                         surfinswdif(i)
           surf_usm_v(1)%rad_sw_res(m) = surfins(i)
           surf_usm_v(1)%rad_lw_in(m)  = surfinlw(i)
           surf_usm_v(1)%rad_lw_out(m) = surfoutlw(i)
           surf_usm_v(1)%rad_net(m)    = surfinsw(i) - surfoutsw(i) +       &
                                         surfinlw(i) - surfoutlw(i)
           surf_usm_v(1)%rad_net_l(m)  = surf_usm_v(1)%rad_net(m)
           surf_usm_v(1)%rad_lw_dif(m) = surfinlwdif(i)
           surf_usm_v(1)%rad_lw_ref(m) = surfinlw(i) - surfinlwdif(i)
           surf_usm_v(1)%rad_lw_res(m) = surfinl(i)
!
!--     eastward-facing
        ELSEIF ( surfl(1,i) == ieast_u )  THEN
           surf_usm_v(2)%rad_sw_in(m)  = surfinsw(i)
           surf_usm_v(2)%rad_sw_out(m) = surfoutsw(i)
           surf_usm_v(2)%rad_sw_dir(m) = surfinswdir(i)
           surf_usm_v(2)%rad_sw_dif(m) = surfinswdif(i)
           surf_usm_v(2)%rad_sw_ref(m) = surfinsw(i) - surfinswdir(i) -     &
                                         surfinswdif(i)
           surf_usm_v(2)%rad_sw_res(m) = surfins(i)
           surf_usm_v(2)%rad_lw_in(m)  = surfinlw(i)
           surf_usm_v(2)%rad_lw_out(m) = surfoutlw(i)
           surf_usm_v(2)%rad_net(m)    = surfinsw(i) - surfoutsw(i) +       &
                                         surfinlw(i) - surfoutlw(i)
           surf_usm_v(2)%rad_net_l(m)  = surf_usm_v(2)%rad_net(m)
           surf_usm_v(2)%rad_lw_dif(m) = surfinlwdif(i)
           surf_usm_v(2)%rad_lw_ref(m) = surfinlw(i) - surfinlwdif(i)
           surf_usm_v(2)%rad_lw_res(m) = surfinl(i)
!
!--     westward-facding
        ELSEIF ( surfl(1,i) == iwest_u )  THEN
           surf_usm_v(3)%rad_sw_in(m)  = surfinsw(i)
           surf_usm_v(3)%rad_sw_out(m) = surfoutsw(i)
           surf_usm_v(3)%rad_sw_dir(m) = surfinswdir(i)
           surf_usm_v(3)%rad_sw_dif(m) = surfinswdif(i)
           surf_usm_v(3)%rad_sw_ref(m) = surfinsw(i) - surfinswdir(i) -     &
                                         surfinswdif(i)
           surf_usm_v(3)%rad_sw_res(m) = surfins(i)
           surf_usm_v(3)%rad_lw_in(m)  = surfinlw(i)
           surf_usm_v(3)%rad_lw_out(m) = surfoutlw(i)
           surf_usm_v(3)%rad_net(m)    = surfinsw(i) - surfoutsw(i) +       &
                                         surfinlw(i) - surfoutlw(i)
           surf_usm_v(3)%rad_net_l(m)  = surf_usm_v(3)%rad_net(m)
           surf_usm_v(3)%rad_lw_dif(m) = surfinlwdif(i)
           surf_usm_v(3)%rad_lw_ref(m) = surfinlw(i) - surfinlwdif(i)
           surf_usm_v(3)%rad_lw_res(m) = surfinl(i)
!
!--     (2) land surfaces
!--     upward-facing
        ELSEIF ( surfl(1,i) == iup_l )  THEN
           surf_lsm_h%rad_sw_in(m)  = surfinsw(i)
           surf_lsm_h%rad_sw_out(m) = surfoutsw(i)
           surf_lsm_h%rad_sw_dir(m) = surfinswdir(i)
           surf_lsm_h%rad_sw_dif(m) = surfinswdif(i)
           surf_lsm_h%rad_sw_ref(m) = surfinsw(i) - surfinswdir(i) -        &
                                         surfinswdif(i)
           surf_lsm_h%rad_sw_res(m) = surfins(i)
           surf_lsm_h%rad_lw_in(m)  = surfinlw(i)
           surf_lsm_h%rad_lw_out(m) = surfoutlw(i)
           surf_lsm_h%rad_net(m)    = surfinsw(i) - surfoutsw(i) +          &
                                      surfinlw(i) - surfoutlw(i)
           surf_lsm_h%rad_lw_dif(m) = surfinlwdif(i)
           surf_lsm_h%rad_lw_ref(m) = surfinlw(i) - surfinlwdif(i)
           surf_lsm_h%rad_lw_res(m) = surfinl(i)
!
!--     northward-facding
        ELSEIF ( surfl(1,i) == inorth_l )  THEN
           surf_lsm_v(0)%rad_sw_in(m)  = surfinsw(i)
           surf_lsm_v(0)%rad_sw_out(m) = surfoutsw(i)
           surf_lsm_v(0)%rad_sw_dir(m) = surfinswdir(i)
           surf_lsm_v(0)%rad_sw_dif(m) = surfinswdif(i)
           surf_lsm_v(0)%rad_sw_ref(m) = surfinsw(i) - surfinswdir(i) -     &
                                         surfinswdif(i)
           surf_lsm_v(0)%rad_sw_res(m) = surfins(i)
           surf_lsm_v(0)%rad_lw_in(m)  = surfinlw(i)
           surf_lsm_v(0)%rad_lw_out(m) = surfoutlw(i)
           surf_lsm_v(0)%rad_net(m)    = surfinsw(i) - surfoutsw(i) +       &
                                         surfinlw(i) - surfoutlw(i)
           surf_lsm_v(0)%rad_lw_dif(m) = surfinlwdif(i)
           surf_lsm_v(0)%rad_lw_ref(m) = surfinlw(i) - surfinlwdif(i)
           surf_lsm_v(0)%rad_lw_res(m) = surfinl(i)
!
!--     southward-facding
        ELSEIF ( surfl(1,i) == isouth_l )  THEN
           surf_lsm_v(1)%rad_sw_in(m)  = surfinsw(i)
           surf_lsm_v(1)%rad_sw_out(m) = surfoutsw(i)
           surf_lsm_v(1)%rad_sw_dir(m) = surfinswdir(i)
           surf_lsm_v(1)%rad_sw_dif(m) = surfinswdif(i)
           surf_lsm_v(1)%rad_sw_ref(m) = surfinsw(i) - surfinswdir(i) -     &
                                         surfinswdif(i)
           surf_lsm_v(1)%rad_sw_res(m) = surfins(i)
           surf_lsm_v(1)%rad_lw_in(m)  = surfinlw(i)
           surf_lsm_v(1)%rad_lw_out(m) = surfoutlw(i)
           surf_lsm_v(1)%rad_net(m)    = surfinsw(i) - surfoutsw(i) +       &
                                         surfinlw(i) - surfoutlw(i)
           surf_lsm_v(1)%rad_lw_dif(m) = surfinlwdif(i)
           surf_lsm_v(1)%rad_lw_ref(m) = surfinlw(i) - surfinlwdif(i)
           surf_lsm_v(1)%rad_lw_res(m) = surfinl(i)
!
!--     eastward-facing
        ELSEIF ( surfl(1,i) == ieast_l )  THEN
           surf_lsm_v(2)%rad_sw_in(m)  = surfinsw(i)
           surf_lsm_v(2)%rad_sw_out(m) = surfoutsw(i)
           surf_lsm_v(2)%rad_sw_dir(m) = surfinswdir(i)
           surf_lsm_v(2)%rad_sw_dif(m) = surfinswdif(i)
           surf_lsm_v(2)%rad_sw_ref(m) = surfinsw(i) - surfinswdir(i) -     &
                                         surfinswdif(i)
           surf_lsm_v(2)%rad_sw_res(m) = surfins(i)
           surf_lsm_v(2)%rad_lw_in(m)  = surfinlw(i)
           surf_lsm_v(2)%rad_lw_out(m) = surfoutlw(i)
           surf_lsm_v(2)%rad_net(m)    = surfinsw(i) - surfoutsw(i) +       &
                                         surfinlw(i) - surfoutlw(i)
           surf_lsm_v(2)%rad_lw_dif(m) = surfinlwdif(i)
           surf_lsm_v(2)%rad_lw_ref(m) = surfinlw(i) - surfinlwdif(i)
           surf_lsm_v(2)%rad_lw_res(m) = surfinl(i)
!
!--     westward-facing
        ELSEIF ( surfl(1,i) == iwest_l )  THEN
           surf_lsm_v(3)%rad_sw_in(m)  = surfinsw(i)
           surf_lsm_v(3)%rad_sw_out(m) = surfoutsw(i)
           surf_lsm_v(3)%rad_sw_dir(m) = surfinswdir(i)
           surf_lsm_v(3)%rad_sw_dif(m) = surfinswdif(i)
           surf_lsm_v(3)%rad_sw_ref(m) = surfinsw(i) - surfinswdir(i) -     &
                                         surfinswdif(i)
           surf_lsm_v(3)%rad_sw_res(m) = surfins(i)
           surf_lsm_v(3)%rad_lw_in(m)  = surfinlw(i)
           surf_lsm_v(3)%rad_lw_out(m) = surfoutlw(i)
           surf_lsm_v(3)%rad_net(m)    = surfinsw(i) - surfoutsw(i) +       &
                                         surfinlw(i) - surfoutlw(i)
           surf_lsm_v(3)%rad_lw_dif(m) = surfinlwdif(i)
           surf_lsm_v(3)%rad_lw_ref(m) = surfinlw(i) - surfinlwdif(i)
           surf_lsm_v(3)%rad_lw_res(m) = surfinl(i)
        ENDIF

     ENDDO

     DO  m = 1, surf_usm_h%ns
        surf_usm_h%surfhf(m) = surf_usm_h%rad_sw_in(m)  +                   &
                               surf_usm_h%rad_lw_in(m)  -                   &
                               surf_usm_h%rad_sw_out(m) -                   &
                               surf_usm_h%rad_lw_out(m)
     ENDDO
     DO  m = 1, surf_lsm_h%ns
        surf_lsm_h%surfhf(m) = surf_lsm_h%rad_sw_in(m)  +                   &
                               surf_lsm_h%rad_lw_in(m)  -                   &
                               surf_lsm_h%rad_sw_out(m) -                   &
                               surf_lsm_h%rad_lw_out(m)
     ENDDO

     DO  l = 0, 3
!--     urban
        DO  m = 1, surf_usm_v(l)%ns
           surf_usm_v(l)%surfhf(m) = surf_usm_v(l)%rad_sw_in(m)  +          &
                                     surf_usm_v(l)%rad_lw_in(m)  -          &
                                     surf_usm_v(l)%rad_sw_out(m) -          &
                                     surf_usm_v(l)%rad_lw_out(m)
        ENDDO
!--     land
        DO  m = 1, surf_lsm_v(l)%ns
           surf_lsm_v(l)%surfhf(m) = surf_lsm_v(l)%rad_sw_in(m)  +          &
                                     surf_lsm_v(l)%rad_lw_in(m)  -          &
                                     surf_lsm_v(l)%rad_sw_out(m) -          &
                                     surf_lsm_v(l)%rad_lw_out(m)

        ENDDO
     ENDDO
!
!--  Calculate the average temperature, albedo, and emissivity for urban/land
!--  domain when using average_radiation in the respective radiation model

!--  calculate horizontal area
! !!! ATTENTION!!! uniform grid is assumed here
     area_hor = (nx+1) * (ny+1) * dx * dy
!
!--  absorbed/received SW & LW and emitted LW energy of all physical
!--  surfaces (land and urban) in local processor
     pinswl = 0._wp
     pinlwl = 0._wp
     pabsswl = 0._wp
     pabslwl = 0._wp
     pemitlwl = 0._wp
     emiss_sum_surfl = 0._wp
     area_surfl = 0._wp
     DO  i = 1, nsurfl
        d = surfl(id, i)
!--  received SW & LW
        pinswl = pinswl + (surfinswdir(i) + surfinswdif(i)) * facearea(d)
        pinlwl = pinlwl + surfinlwdif(i) * facearea(d)
!--   absorbed SW & LW
        pabsswl = pabsswl + (1._wp - albedo_surf(i)) *                   &
                                                surfinsw(i) * facearea(d)
        pabslwl = pabslwl + emiss_surf(i) * surfinlw(i) * facearea(d)
!--   emitted LW
        pemitlwl = pemitlwl + surfemitlwl(i) * facearea(d)
!--   emissivity and area sum
        emiss_sum_surfl = emiss_sum_surfl + emiss_surf(i) * facearea(d)
        area_surfl = area_surfl + facearea(d)
     END DO
!
!--  add the absorbed SW energy by plant canopy
     IF ( npcbl > 0 )  THEN
        pabsswl = pabsswl + SUM(pcbinsw)
        pabslwl = pabslwl + SUM(pcbinlw)
        pinswl  = pinswl + SUM(pcbinswdir) + SUM(pcbinswdif)
     ENDIF
!
!--  gather all rad flux energy in all processors. In order to reduce
!--  the number of MPI calls (to reduce latencies), combine the required
!--  quantities in one array, sum it up, and subsequently re-distribute
!--  back to the respective quantities.
#if defined( __parallel )
     combine_allreduce_l(1) = pinswl
     combine_allreduce_l(2) = pinlwl
     combine_allreduce_l(3) = pabsswl
     combine_allreduce_l(4) = pabslwl
     combine_allreduce_l(5) = pemitlwl
     combine_allreduce_l(6) = emiss_sum_surfl
     combine_allreduce_l(7) = area_surfl

     CALL MPI_ALLREDUCE( combine_allreduce_l,                                  &
                         combine_allreduce,                                    &
                         SIZE( combine_allreduce ),                            &
                         MPI_REAL,                                             &
                         MPI_SUM,                                              &
                         comm2d,                                               &
                         ierr )

     pinsw          = combine_allreduce(1)
     pinlw          = combine_allreduce(2)
     pabssw         = combine_allreduce(3)
     pabslw         = combine_allreduce(4)
     pemitlw        = combine_allreduce(5)
     emiss_sum_surf = combine_allreduce(6)
     area_surf      = combine_allreduce(7)
#else
     pinsw          = pinswl
     pinlw          = pinlwl
     pabssw         = pabsswl
     pabslw         = pabslwl
     pemitlw        = pemitlwl
     emiss_sum_surf = emiss_sum_surfl
     area_surf      = area_surfl
#endif

!--  (1) albedo
     IF ( pinsw /= 0.0_wp )  albedo_urb = ( pinsw - pabssw ) / pinsw
!--  (2) average emmsivity
     IF ( area_surf /= 0.0_wp )  emissivity_urb = emiss_sum_surf / area_surf
!
!--  Temporally comment out calculation of effective radiative temperature.
!--  See below for more explanation.
!--  (3) temperature
!--   first we calculate an effective horizontal area to account for
!--   the effect of vertical surfaces (which contributes to LW emission)
!--   We simply use the ratio of the total LW to the incoming LW flux
      area_hor = pinlw / rad_lw_in_diff(nyn,nxl)
      t_rad_urb = ( ( pemitlw - pabslw + emissivity_urb * pinlw ) / &
           (emissivity_urb * sigma_sb * area_hor) )**0.25_wp

     IF ( debug_output_timestep )  CALL debug_message( 'radiation_interaction', 'end' )


    CONTAINS

!------------------------------------------------------------------------------!
!> Calculates radiation absorbed by box with given size and LAD.
!>
!> Simulates resol**2 rays (by equally spacing a bounding horizontal square
!> conatining all possible rays that would cross the box) and calculates
!> average transparency per ray. Returns fraction of absorbed radiation flux
!> and area for which this fraction is effective.
!------------------------------------------------------------------------------!
    PURE SUBROUTINE box_absorb(boxsize, resol, dens, uvec, area, absorb)
       IMPLICIT NONE

       REAL(wp), DIMENSION(3), INTENT(in) :: &
            boxsize, &      !< z, y, x size of box in m
            uvec            !< z, y, x unit vector of incoming flux
       INTEGER(iwp), INTENT(in) :: &
            resol           !< No. of rays in x and y dimensions
       REAL(wp), INTENT(in) :: &
            dens            !< box density (e.g. Leaf Area Density)
       REAL(wp), INTENT(out) :: &
            area, &         !< horizontal area for flux absorbtion
            absorb          !< fraction of absorbed flux
       REAL(wp) :: &
            xshift, yshift, &
            xmin, xmax, ymin, ymax, &
            xorig, yorig, &
            dx1, dy1, dz1, dx2, dy2, dz2, &
            crdist, &
            transp
       INTEGER(iwp) :: &
            i, j

       xshift = uvec(3) / uvec(1) * boxsize(1)
       xmin = min(0._wp, -xshift)
       xmax = boxsize(3) + max(0._wp, -xshift)
       yshift = uvec(2) / uvec(1) * boxsize(1)
       ymin = min(0._wp, -yshift)
       ymax = boxsize(2) + max(0._wp, -yshift)

       transp = 0._wp
       DO i = 1, resol
          xorig = xmin + (xmax-xmin) * (i-.5_wp) / resol
          DO j = 1, resol
             yorig = ymin + (ymax-ymin) * (j-.5_wp) / resol

             dz1 = 0._wp
             dz2 = boxsize(1)/uvec(1)

             IF ( uvec(2) > 0._wp )  THEN
                dy1 = -yorig             / uvec(2) !< crossing with y=0
                dy2 = (boxsize(2)-yorig) / uvec(2) !< crossing with y=boxsize(2)
             ELSE !uvec(2)==0
                dy1 = -huge(1._wp)
                dy2 = huge(1._wp)
             ENDIF

             IF ( uvec(3) > 0._wp )  THEN
                dx1 = -xorig             / uvec(3) !< crossing with x=0
                dx2 = (boxsize(3)-xorig) / uvec(3) !< crossing with x=boxsize(3)
             ELSE !uvec(3)==0
                dx1 = -huge(1._wp)
                dx2 = huge(1._wp)
             ENDIF

             crdist = max(0._wp, (min(dz2, dy2, dx2) - max(dz1, dy1, dx1)))
             transp = transp + exp(-ext_coef * dens * crdist)
          ENDDO
       ENDDO
       transp = transp / resol**2
       area = (boxsize(3)+xshift)*(boxsize(2)+yshift)
       absorb = 1._wp - transp

    END SUBROUTINE box_absorb

!------------------------------------------------------------------------------!
! Description:
! ------------
!> This subroutine splits direct and diffusion dw radiation for RTM processing.
!> It sould not be called in case the radiation model already does it
!> It follows Boland, Ridley & Brown (2008)
!------------------------------------------------------------------------------!
    SUBROUTINE calc_diffusion_radiation

       USE palm_date_time_mod,                                                 &
           ONLY:  seconds_per_day

       INTEGER(iwp)        ::  i                        !< grid index x-direction
       INTEGER(iwp)        ::  j                        !< grid index y-direction
       INTEGER(iwp)        ::  days_per_year            !< days in the current year

       REAL(wp)            ::  clearnessIndex           !< clearness index
       REAL(wp)            ::  corrected_solarUp        !< corrected solar up radiation
       REAL(wp)            ::  diff_frac                !< diffusion fraction of the radiation
       REAL(wp)            ::  etr                      !< extraterestrial radiation
       REAL(wp)            ::  horizontalETR            !< horizontal extraterestrial radiation
       REAL(wp), PARAMETER ::  lowest_solarUp = 0.1_wp  !< limit the sun elevation to protect stability of the calculation
       REAL(wp)            ::  second_of_year           !< current second of the year
       REAL(wp)            ::  year_angle               !< angle

!
!--     Calculate current day and time based on the initial values and simulation time
        CALL get_date_time( time_since_reference_point,      &
                            second_of_year = second_of_year, &
                            days_per_year = days_per_year    )
        year_angle = second_of_year / ( REAL( days_per_year, KIND=wp ) * seconds_per_day ) &
                   * 2.0_wp * pi

        etr = solar_constant * (1.00011_wp +                                   &
                          0.034221_wp * cos(year_angle) +                      &
                          0.001280_wp * sin(year_angle) +                      &
                          0.000719_wp * cos(2.0_wp * year_angle) +             &
                          0.000077_wp * sin(2.0_wp * year_angle))

!--
!--     Under a very low angle, we keep extraterestrial radiation at
!--     the last small value, therefore the clearness index will be pushed
!--     towards 0 while keeping full continuity.
        IF ( cos_zenith <= lowest_solarUp )  THEN
            corrected_solarUp = lowest_solarUp
        ELSE
            corrected_solarUp = cos_zenith
        ENDIF

        horizontalETR = etr * corrected_solarUp

        DO i = nxl, nxr
            DO j = nys, nyn
                clearnessIndex = rad_sw_in(0,j,i) / horizontalETR
                diff_frac = 1.0_wp / (1.0_wp + exp(-5.0033_wp + 8.6025_wp * clearnessIndex))
                rad_sw_in_diff(j,i) = rad_sw_in(0,j,i) * diff_frac
                rad_sw_in_dir(j,i)  = rad_sw_in(0,j,i) * (1.0_wp - diff_frac)
                rad_lw_in_diff(j,i) = rad_lw_in(0,j,i)
            ENDDO
        ENDDO

    END SUBROUTINE calc_diffusion_radiation

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Print consecutive radiative extremes if requested to trace early radiation
!> interaction instabilities.
!------------------------------------------------------------------------------!
    SUBROUTINE radiation_print_debug_surf( description, values, step )

       CHARACTER (LEN=*), INTENT(in)      ::  description
       REAL(wp), DIMENSION(:), INTENT(in) ::  values
       INTEGER(iwp), INTENT(in), OPTIONAL ::  step

       CHARACTER (LEN=50)                 ::  location
       CHARACTER (LEN=1024)               ::  debug_string
       INTEGER                            ::  isurf
       REAL(wp)                           ::  x

       isurf = MAXLOC( values, DIM=1 )
       x = values(isurf)
       IF ( x < trace_fluxes_above )  RETURN

       IF ( PRESENT( step ) )  THEN
          WRITE( location, '(A," #",I0)' ) description, step
       ELSE
          location = description
       ENDIF

       WRITE( debug_string, '("Maximum of ",A50," = ",F12.1," at coords ' //  &
                            'i=",I4,", j=",I4,", k=",I4,", d=",I1,". '    //  &
                            'Alb=",F7.3,", emis=",F7.3)'                    ) &
              location, x, surfl(ix,isurf), surfl(iy,isurf),                  &
              surfl(iz,isurf), surfl(id,isurf), albedo_surf(isurf),           &
              emiss_surf(isurf)
       CALL debug_message( debug_string, 'info' )

    END SUBROUTINE

    SUBROUTINE radiation_print_debug_pcb( description, values, step )

       CHARACTER (LEN=*), INTENT(in)      ::  description
       REAL(wp), DIMENSION(:), INTENT(in) ::  values
       INTEGER(iwp), INTENT(in), OPTIONAL ::  step

       CHARACTER (LEN=50)                 ::  location
       CHARACTER (LEN=1024)               ::  debug_string
       INTEGER                            ::  ipcb
       REAL(wp)                           ::  x

       IF ( npcbl <= 0 )  RETURN
       ipcb = MAXLOC( values, DIM=1 )
       x = values(ipcb) / (dx*dy*dz(1))
       IF ( x < trace_fluxes_above )  RETURN

       IF ( PRESENT( step ) )  THEN
          WRITE( location, '(A," #",I0)' ) description, step
       ELSE
          location = description
       ENDIF

       WRITE( debug_string, '("Maximum of ",A50," = ",F12.1," at coords ' //  &
                            'i=",I4,", j=",I4,", k=",I4)'                   ) &
              location, x, pcbl(ix,ipcb), pcbl(iy,ipcb), pcbl(iz,ipcb)
       CALL debug_message( debug_string, 'info' )

    END SUBROUTINE

    SUBROUTINE radiation_print_debug_horz( description, values, step )

       CHARACTER (LEN=*), INTENT(in)        ::  description
       REAL(wp), DIMENSION(:,:), INTENT(in) ::  values
       INTEGER(iwp), INTENT(in), OPTIONAL   ::  step

       CHARACTER (LEN=50)                   ::  location
       CHARACTER (LEN=1024)                 ::  debug_string
       INTEGER, DIMENSION(2)                ::  ji
       REAL(wp)                             ::  x

       ji = MAXLOC( values )
       x = values(ji(1),ji(2))
       IF ( x < trace_fluxes_above )  RETURN

       IF ( PRESENT( step ) )  THEN
          WRITE( location, '(A," #",I0)' ) description, step
       ELSE
          location = description
       ENDIF

       WRITE( debug_string, '("Maximum of ",A50," = ",F12.1," at coords ' //  &
                            'i=",I4,", j=",I4)'                             ) &
              location, x, ji(2), ji(1)
       CALL debug_message( debug_string, 'info' )

    END SUBROUTINE

 END SUBROUTINE radiation_interaction

!------------------------------------------------------------------------------!
! Description:
! ------------
!> This subroutine initializes structures needed for Radiative Transfer
!> Model (RTM). This model calculates transformation processes of the
!> radiation inside urban and land canopy layer. The module includes also
!> the interaction of the radiation with the resolved plant canopy.
!>
!------------------------------------------------------------------------------!
    SUBROUTINE radiation_interaction_init

       USE control_parameters,                                                 &
           ONLY:  dz_stretch_level_start

       USE plant_canopy_model_mod,                                             &
           ONLY:  lad_s

       IMPLICIT NONE

       INTEGER(iwp) :: i, j, k, l, m, d
       INTEGER(iwp) :: k_topo     !< vertical index indicating topography top for given (j,i)
       INTEGER(iwp) :: nzptl, nzubl, nzutl, isurf, ipcgb, imrt
       REAL(wp)     :: mrl
#if defined( __parallel )
       INTEGER(iwp), DIMENSION(:), POINTER, SAVE ::  gridsurf_rma   !< fortran pointer, but lower bounds are 1
       TYPE(c_ptr)                               ::  gridsurf_rma_p !< allocated c pointer
       INTEGER(iwp)                              ::  minfo          !< MPI RMA window info handle
#endif

!
!--     precalculate face areas for different face directions using normal vector
        DO d = 0, nsurf_type
            facearea(d) = 1._wp
            IF ( idir(d) == 0 ) facearea(d) = facearea(d) * dx
            IF ( jdir(d) == 0 ) facearea(d) = facearea(d) * dy
            IF ( kdir(d) == 0 ) facearea(d) = facearea(d) * dz(1)
        ENDDO
!
!--    Find nz_urban_b, nz_urban_t, nz_urban via wall_flag_0 array (nzb_s_inner will be
!--    removed later). The following contruct finds the lowest / largest index
!--    for any upward-facing wall (see bit 12).
       nzubl = MINVAL( topo_top_ind(nys:nyn,nxl:nxr,0) )
       nzutl = MAXVAL( topo_top_ind(nys:nyn,nxl:nxr,0) )

       nzubl = MAX( nzubl, nzb )

       IF ( plant_canopy )  THEN
!--        allocate needed arrays
           ALLOCATE( pct(nys:nyn,nxl:nxr) )
           ALLOCATE( pch(nys:nyn,nxl:nxr) )

!--        calculate plant canopy height
           npcbl = 0
           pct   = 0
           pch   = 0
           DO i = nxl, nxr
               DO j = nys, nyn
!
!--                Find topography top index
                   k_topo = topo_top_ind(j,i,0)

                   DO k = nzt+1, 0, -1
                       IF ( lad_s(k,j,i) /= 0.0_wp )  THEN
!--                        we are at the top of the pcs
                           pct(j,i) = k + k_topo
                           pch(j,i) = k
                           npcbl = npcbl + pch(j,i)
                           EXIT
                       ENDIF
                   ENDDO
               ENDDO
           ENDDO

           nzutl = MAX( nzutl, MAXVAL( pct ) )
           nzptl = MAXVAL( pct )

           prototype_lad = MAXVAL( lad_s ) * .9_wp  !< better be *1.0 if lad is either 0 or maxval(lad) everywhere
           IF ( prototype_lad <= 0._wp ) prototype_lad = .3_wp
           !WRITE(message_string, '(a,f6.3)') 'Precomputing effective box optical ' &
           !    // 'depth using prototype leaf area density = ', prototype_lad
           !CALL message('radiation_interaction_init', 'PA0520', 0, 0, -1, 6, 0)
       ENDIF

       nzutl = MIN( nzutl + nzut_free, nzt )

#if defined( __parallel )
       CALL MPI_AllReduce(nzubl, nz_urban_b, 1, MPI_INTEGER, MPI_MIN, comm2d, ierr )
       IF ( ierr /= 0 ) THEN
           WRITE(9,*) 'Error MPI_AllReduce11:', ierr, nzubl, nz_urban_b
           FLUSH(9)
       ENDIF
       CALL MPI_AllReduce(nzutl, nz_urban_t, 1, MPI_INTEGER, MPI_MAX, comm2d, ierr )
       IF ( ierr /= 0 ) THEN
           WRITE(9,*) 'Error MPI_AllReduce12:', ierr, nzutl, nz_urban_t
           FLUSH(9)
       ENDIF
       CALL MPI_AllReduce(nzptl, nz_plant_t, 1, MPI_INTEGER, MPI_MAX, comm2d, ierr )
       IF ( ierr /= 0 ) THEN
           WRITE(9,*) 'Error MPI_AllReduce13:', ierr, nzptl, nz_plant_t
           FLUSH(9)
       ENDIF
#else
       nz_urban_b = nzubl
       nz_urban_t = nzutl
       nz_plant_t = nzptl
#endif
!
!--    Stretching (non-uniform grid spacing) is not considered in the radiation
!--    model. Therefore, vertical stretching has to be applied above the area
!--    where the parts of the radiation model which assume constant grid spacing
!--    are active. ABS (...) is required because the default value of
!--    dz_stretch_level_start is -9999999.9_wp (negative).
       IF ( ABS( dz_stretch_level_start(1) ) <= zw(nz_urban_t) ) THEN
          WRITE( message_string, * ) 'The lowest level where vertical ',       &
                                     'stretching is applied have to be ',      &
                                     'greater than ', zw(nz_urban_t)
          CALL message( 'radiation_interaction_init', 'PA0496', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    global number of urban and plant layers
       nz_urban = nz_urban_t - nz_urban_b + 1
       nz_plant = nz_plant_t - nz_urban_b + 1
!
!--    check max_raytracing_dist relative to urban surface layer height
       mrl = 2.0_wp * nz_urban * dz(1)
!--    set max_raytracing_dist to double the urban surface layer height, if not set
       IF ( max_raytracing_dist == -999.0_wp ) THEN
          max_raytracing_dist = mrl
       ENDIF
!--    check if max_raytracing_dist set too low (here we only warn the user. Other
!      option is to correct the value again to double the urban surface layer height)
       IF ( max_raytracing_dist  <  mrl ) THEN
          WRITE(message_string, '(a,f6.1)') 'Max_raytracing_dist is set less than ' // &
               'double the urban surface layer height, i.e. ', mrl
          CALL message('radiation_interaction_init', 'PA0521', 0, 0, 0, 6, 0 )
       ENDIF
!        IF ( max_raytracing_dist <= mrl ) THEN
!           IF ( max_raytracing_dist /= -999.0_wp ) THEN
! !--          max_raytracing_dist too low
!              WRITE(message_string, '(a,f6.1)') 'Max_raytracing_dist too low, ' &
!                    // 'override to value ', mrl
!              CALL message('radiation_interaction_init', 'PA0521', 0, 0, -1, 6, 0)
!           ENDIF
!           max_raytracing_dist = mrl
!        ENDIF
!
!--    allocate urban surfaces grid
!--    calc number of surfaces in local proc
       IF ( debug_output )  CALL debug_message( 'calculation of indices for surfaces', 'info' )

       nsurfl = 0
!
!--    Number of horizontal surfaces including land- and roof surfaces in both USM and LSM. Note that
!--    All horizontal surface elements are already counted in surface_mod.
       startland = 1
       nsurfl    = surf_usm_h%ns + surf_lsm_h%ns
       endland   = nsurfl
       nlands    = endland - startland + 1

!
!--    Number of vertical surfaces in both USM and LSM. Note that all vertical surface elements are
!--    already counted in surface_mod.
       startwall = nsurfl+1
       DO  i = 0,3
          nsurfl = nsurfl + surf_usm_v(i)%ns + surf_lsm_v(i)%ns
       ENDDO
       endwall = nsurfl
       nwalls  = endwall - startwall + 1
       dirstart = (/ startland, startwall, startwall, startwall, startwall /)
       dirend = (/ endland, endwall, endwall, endwall, endwall /)

!--    fill gridpcbl and pcbl
       IF ( npcbl > 0 )  THEN
           ALLOCATE( pcbl(iz:ix, 1:npcbl) )
           ALLOCATE( gridpcbl(nz_urban_b:nz_plant_t,nys:nyn,nxl:nxr) )
           pcbl = -1
           gridpcbl(:,:,:) = 0
           ipcgb = 0
           DO i = nxl, nxr
               DO j = nys, nyn
!
!--                Find topography top index
                   k_topo = topo_top_ind(j,i,0)

                   DO k = k_topo + 1, pct(j,i)
                       ipcgb = ipcgb + 1
                       gridpcbl(k,j,i) = ipcgb
                       pcbl(:,ipcgb) = (/ k, j, i /)
                   ENDDO
               ENDDO
           ENDDO
           ALLOCATE( pcbinsw( 1:npcbl ) )
           ALLOCATE( pcbinswdir( 1:npcbl ) )
           ALLOCATE( pcbinswdif( 1:npcbl ) )
           ALLOCATE( pcbinlw( 1:npcbl ) )
       ENDIF

!
!--    Fill surfl (the ordering of local surfaces given by the following
!--    cycles must not be altered, certain file input routines may depend
!--    on it).
!
!--    We allocate the array as linear and then use a two-dimensional pointer
!--    into it, because some MPI implementations crash with 2D-allocated arrays.
       ALLOCATE(surfl_linear(nidx_surf*nsurfl))
       surfl(1:nidx_surf,1:nsurfl) => surfl_linear(1:nidx_surf*nsurfl)
       isurf = 0
       IF ( rad_angular_discretization )  THEN
!
!--       Allocate and fill the reverse indexing array gridsurf
#if defined( __parallel )
!
!--       raytrace_mpi_rma is asserted

          CALL MPI_Info_create(minfo, ierr)
          IF ( ierr /= 0 ) THEN
              WRITE(9,*) 'Error MPI_Info_create1:', ierr
              FLUSH(9)
          ENDIF
          CALL MPI_Info_set(minfo, 'accumulate_ordering', 'none', ierr)
          IF ( ierr /= 0 ) THEN
              WRITE(9,*) 'Error MPI_Info_set1:', ierr
              FLUSH(9)
          ENDIF
          CALL MPI_Info_set(minfo, 'accumulate_ops', 'same_op', ierr)
          IF ( ierr /= 0 ) THEN
              WRITE(9,*) 'Error MPI_Info_set2:', ierr
              FLUSH(9)
          ENDIF
          CALL MPI_Info_set(minfo, 'same_size', 'true', ierr)
          IF ( ierr /= 0 ) THEN
              WRITE(9,*) 'Error MPI_Info_set3:', ierr
              FLUSH(9)
          ENDIF
          CALL MPI_Info_set(minfo, 'same_disp_unit', 'true', ierr)
          IF ( ierr /= 0 ) THEN
              WRITE(9,*) 'Error MPI_Info_set4:', ierr
              FLUSH(9)
          ENDIF

          CALL MPI_Win_allocate(INT(STORAGE_SIZE(1_iwp)/8*nsurf_type_u*nz_urban*nny*nnx, &
                                    kind=MPI_ADDRESS_KIND), STORAGE_SIZE(1_iwp)/8,  &
                                minfo, comm2d, gridsurf_rma_p, win_gridsurf, ierr)
          IF ( ierr /= 0 ) THEN
              WRITE(9,*) 'Error MPI_Win_allocate1:', ierr, &
                 INT(STORAGE_SIZE(1_iwp)/8*nsurf_type_u*nz_urban*nny*nnx,kind=MPI_ADDRESS_KIND), &
                 STORAGE_SIZE(1_iwp)/8, win_gridsurf
              FLUSH(9)
          ENDIF

          CALL MPI_Info_free(minfo, ierr)
          IF ( ierr /= 0 ) THEN
              WRITE(9,*) 'Error MPI_Info_free1:', ierr
              FLUSH(9)
          ENDIF

!
!--       On Intel compilers, calling c_f_pointer to transform a C pointer
!--       directly to a multi-dimensional Fotran pointer leads to strange
!--       errors on dimension boundaries. However, transforming to a 1D
!--       pointer and then redirecting a multidimensional pointer to it works
!--       fine.
          CALL c_f_pointer(gridsurf_rma_p, gridsurf_rma, (/ nsurf_type_u*nz_urban*nny*nnx /))
          gridsurf(0:nsurf_type_u-1, nz_urban_b:nz_urban_t, nys:nyn, nxl:nxr) =>                &
                     gridsurf_rma(1:nsurf_type_u*nz_urban*nny*nnx)
#else
          ALLOCATE(gridsurf(0:nsurf_type_u-1,nz_urban_b:nz_urban_t,nys:nyn,nxl:nxr) )
#endif
          gridsurf(:,:,:,:) = -999
       ENDIF

!--    add horizontal surface elements (land and urban surfaces)
!--    TODO: add urban overhanging surfaces (idown_u)
       DO i = nxl, nxr
           DO j = nys, nyn
              DO  m = surf_usm_h%start_index(j,i), surf_usm_h%end_index(j,i)
                 k = surf_usm_h%k(m)
                 isurf = isurf + 1
                 surfl(:,isurf) = (/iup_u,k,j,i,m/)
                 IF ( rad_angular_discretization ) THEN
                    gridsurf(iup_u,k,j,i) = isurf
                 ENDIF
              ENDDO

              DO  m = surf_lsm_h%start_index(j,i), surf_lsm_h%end_index(j,i)
                 k = surf_lsm_h%k(m)
                 isurf = isurf + 1
                 surfl(:,isurf) = (/iup_l,k,j,i,m/)
                 IF ( rad_angular_discretization ) THEN
                    gridsurf(iup_u,k,j,i) = isurf
                 ENDIF
              ENDDO

           ENDDO
       ENDDO

!--    add vertical surface elements (land and urban surfaces)
!--    TODO: remove the hard coding of l = 0 to l = idirection
       DO i = nxl, nxr
           DO j = nys, nyn
              l = 0
              DO  m = surf_usm_v(l)%start_index(j,i), surf_usm_v(l)%end_index(j,i)
                 k = surf_usm_v(l)%k(m)
                 isurf = isurf + 1
                 surfl(:,isurf) = (/inorth_u,k,j,i,m/)
                 IF ( rad_angular_discretization ) THEN
                    gridsurf(inorth_u,k,j,i) = isurf
                 ENDIF
              ENDDO
              DO  m = surf_lsm_v(l)%start_index(j,i), surf_lsm_v(l)%end_index(j,i)
                 k = surf_lsm_v(l)%k(m)
                 isurf = isurf + 1
                 surfl(:,isurf) = (/inorth_l,k,j,i,m/)
                 IF ( rad_angular_discretization ) THEN
                    gridsurf(inorth_u,k,j,i) = isurf
                 ENDIF
              ENDDO

              l = 1
              DO  m = surf_usm_v(l)%start_index(j,i), surf_usm_v(l)%end_index(j,i)
                 k = surf_usm_v(l)%k(m)
                 isurf = isurf + 1
                 surfl(:,isurf) = (/isouth_u,k,j,i,m/)
                 IF ( rad_angular_discretization ) THEN
                    gridsurf(isouth_u,k,j,i) = isurf
                 ENDIF
              ENDDO
              DO  m = surf_lsm_v(l)%start_index(j,i), surf_lsm_v(l)%end_index(j,i)
                 k = surf_lsm_v(l)%k(m)
                 isurf = isurf + 1
                 surfl(:,isurf) = (/isouth_l,k,j,i,m/)
                 IF ( rad_angular_discretization ) THEN
                    gridsurf(isouth_u,k,j,i) = isurf
                 ENDIF
              ENDDO

              l = 2
              DO  m = surf_usm_v(l)%start_index(j,i), surf_usm_v(l)%end_index(j,i)
                 k = surf_usm_v(l)%k(m)
                 isurf = isurf + 1
                 surfl(:,isurf) = (/ieast_u,k,j,i,m/)
                 IF ( rad_angular_discretization ) THEN
                    gridsurf(ieast_u,k,j,i) = isurf
                 ENDIF
              ENDDO
              DO  m = surf_lsm_v(l)%start_index(j,i), surf_lsm_v(l)%end_index(j,i)
                 k = surf_lsm_v(l)%k(m)
                 isurf = isurf + 1
                 surfl(:,isurf) = (/ieast_l,k,j,i,m/)
                 IF ( rad_angular_discretization ) THEN
                    gridsurf(ieast_u,k,j,i) = isurf
                 ENDIF
              ENDDO

              l = 3
              DO  m = surf_usm_v(l)%start_index(j,i), surf_usm_v(l)%end_index(j,i)
                 k = surf_usm_v(l)%k(m)
                 isurf = isurf + 1
                 surfl(:,isurf) = (/iwest_u,k,j,i,m/)
                 IF ( rad_angular_discretization ) THEN
                    gridsurf(iwest_u,k,j,i) = isurf
                 ENDIF
              ENDDO
              DO  m = surf_lsm_v(l)%start_index(j,i), surf_lsm_v(l)%end_index(j,i)
                 k = surf_lsm_v(l)%k(m)
                 isurf = isurf + 1
                 surfl(:,isurf) = (/iwest_l,k,j,i,m/)
                 IF ( rad_angular_discretization ) THEN
                    gridsurf(iwest_u,k,j,i) = isurf
                 ENDIF
              ENDDO
           ENDDO
       ENDDO
!
!--    Add local MRT boxes for specified number of levels
       nmrtbl = 0
       IF ( mrt_nlevels > 0 )  THEN
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  m = surf_usm_h%start_index(j,i), surf_usm_h%end_index(j,i)
!
!--                Skip roof if requested
                   IF ( mrt_skip_roof  .AND.  surf_usm_h%isroof_surf(m) )  CYCLE
!
!--                Cycle over specified no of levels
                   nmrtbl = nmrtbl + mrt_nlevels
                ENDDO
!
!--             Dtto for LSM
                DO  m = surf_lsm_h%start_index(j,i), surf_lsm_h%end_index(j,i)
                   nmrtbl = nmrtbl + mrt_nlevels
                ENDDO
             ENDDO
          ENDDO

          ALLOCATE( mrtbl(iz:ix,nmrtbl), mrtsky(nmrtbl), mrtskyt(nmrtbl), &
                    mrtinsw(nmrtbl), mrtinlw(nmrtbl), mrt(nmrtbl) )

          imrt = 0
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  m = surf_usm_h%start_index(j,i), surf_usm_h%end_index(j,i)
!
!--                Skip roof if requested
                   IF ( mrt_skip_roof  .AND.  surf_usm_h%isroof_surf(m) )  CYCLE
!
!--                Cycle over specified no of levels
                   l = surf_usm_h%k(m)
                   DO  k = l, l + mrt_nlevels - 1
                      imrt = imrt + 1
                      mrtbl(:,imrt) = (/k,j,i/)
                   ENDDO
                ENDDO
!
!--             Dtto for LSM
                DO  m = surf_lsm_h%start_index(j,i), surf_lsm_h%end_index(j,i)
                   l = surf_lsm_h%k(m)
                   DO  k = l, l + mrt_nlevels - 1
                      imrt = imrt + 1
                      mrtbl(:,imrt) = (/k,j,i/)
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDIF

!
!--    broadband albedo of the land, roof and wall surface
!--    for domain border and sky set artifically to 1.0
!--    what allows us to calculate heat flux leaving over
!--    side and top borders of the domain
       ALLOCATE ( albedo_surf(nsurfl) )
       albedo_surf = 1.0_wp
!
!--    Also allocate further array for emissivity with identical order of
!--    surface elements as radiation arrays.
       ALLOCATE ( emiss_surf(nsurfl)  )


!
!--    global array surf of indices of surfaces and displacement index array surfstart
       ALLOCATE(nsurfs(0:numprocs-1))

#if defined( __parallel )
       CALL MPI_Allgather(nsurfl,1,MPI_INTEGER,nsurfs,1,MPI_INTEGER,comm2d,ierr)
       IF ( ierr /= 0 ) THEN
         WRITE(9,*) 'Error MPI_AllGather1:', ierr, nsurfl, nsurfs
         FLUSH(9)
     ENDIF

#else
       nsurfs(0) = nsurfl
#endif
       ALLOCATE(surfstart(0:numprocs))
       k = 0
       DO i=0,numprocs-1
           surfstart(i) = k
           k = k+nsurfs(i)
       ENDDO
       surfstart(numprocs) = k
       nsurf = k
!
!--    We allocate the array as linear and then use a two-dimensional pointer
!--    into it, because some MPI implementations crash with 2D-allocated arrays.
       ALLOCATE(surf_linear(nidx_surf*nsurf))
       surf(1:nidx_surf,1:nsurf) => surf_linear(1:nidx_surf*nsurf)

#if defined( __parallel )
       CALL MPI_AllGatherv(surfl_linear, nsurfl*nidx_surf, MPI_INTEGER,    &
                           surf_linear, nsurfs*nidx_surf,                  &
                           surfstart(0:numprocs-1)*nidx_surf, MPI_INTEGER, &
                           comm2d, ierr)
       IF ( ierr /= 0 ) THEN
           WRITE(9,*) 'Error MPI_AllGatherv4:', ierr, SIZE(surfl_linear),    &
                      nsurfl*nidx_surf, SIZE(surf_linear), nsurfs*nidx_surf, &
                      surfstart(0:numprocs-1)*nidx_surf
           FLUSH(9)
       ENDIF
#else
       surf = surfl
#endif

!--
!--    allocation of the arrays for direct and diffusion radiation
       IF ( debug_output )  CALL debug_message( 'allocation of radiation arrays', 'info' )
!--    rad_sw_in, rad_lw_in are computed in radiation model,
!--    splitting of direct and diffusion part is done
!--    in calc_diffusion_radiation for now

       ALLOCATE( rad_sw_in_dir(nysg:nyng,nxlg:nxrg) )
       ALLOCATE( rad_sw_in_diff(nysg:nyng,nxlg:nxrg) )
       ALLOCATE( rad_lw_in_diff(nysg:nyng,nxlg:nxrg) )
       rad_sw_in_dir  = 0.0_wp
       rad_sw_in_diff = 0.0_wp
       rad_lw_in_diff = 0.0_wp

!--    allocate radiation arrays
       ALLOCATE( surfins(nsurfl) )
       ALLOCATE( surfinl(nsurfl) )
       ALLOCATE( surfinsw(nsurfl) )
       ALLOCATE( surfinlw(nsurfl) )
       ALLOCATE( surfinswdir(nsurfl) )
       ALLOCATE( surfinswdif(nsurfl) )
       ALLOCATE( surfinlwdif(nsurfl) )
       ALLOCATE( surfoutsl(nsurfl) )
       ALLOCATE( surfoutll(nsurfl) )
       ALLOCATE( surfoutsw(nsurfl) )
       ALLOCATE( surfoutlw(nsurfl) )
       ALLOCATE( surfouts(nsurf) )
       ALLOCATE( surfoutl(nsurf) )
       ALLOCATE( surfinlg(nsurf) )
       ALLOCATE( skyvf(nsurfl) )
       ALLOCATE( skyvft(nsurfl) )
       ALLOCATE( surfemitlwl(nsurfl) )

!
!--    In case of average_radiation, aggregated surface albedo and emissivity,
!--    also set initial value for t_rad_urb.
!--    For now set an arbitrary initial value.
       IF ( average_radiation )  THEN
          albedo_urb = 0.1_wp
          emissivity_urb = 0.9_wp
          t_rad_urb = pt_surface
       ENDIF

    END SUBROUTINE radiation_interaction_init

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculates shape view factors (SVF), plant sink canopy factors (PCSF),
!> sky-view factors, discretized path for direct solar radiation, MRT factors
!> and other preprocessed data needed for radiation_interaction inside RTM.
!> This subroutine is called only one at the beginning of the simulation.
!> The resulting factors can be stored to files and reused with other
!> simulations utilizing the same surface and plant canopy structure.
!------------------------------------------------------------------------------!
    SUBROUTINE radiation_calc_svf

        IMPLICIT NONE

        INTEGER(iwp)                                  :: i, j, k, d, ip, jp
        INTEGER(iwp)                                  :: isvf, ksvf, icsf, kcsf, npcsfl, isvf_surflt, imrt, imrtf, ipcgb
        INTEGER(iwp)                                  :: sd, td
        INTEGER(iwp)                                  :: iaz, izn      !< azimuth, zenith counters
        INTEGER(iwp)                                  :: naz, nzn      !< azimuth, zenith num of steps
        REAL(wp)                                      :: az0, zn0      !< starting azimuth/zenith
        REAL(wp)                                      :: azs, zns      !< azimuth/zenith cycle step
        REAL(wp)                                      :: az1, az2      !< relative azimuth of section borders
        REAL(wp)                                      :: azmid         !< ray (center) azimuth
        REAL(wp)                                      :: yxlen         !< |yxdir|
        REAL(wp), DIMENSION(2)                        :: yxdir         !< y,x *unit* vector of ray direction (in grid units)
        REAL(wp), DIMENSION(:), ALLOCATABLE           :: zdirs         !< directions in z (tangent of elevation)
        REAL(wp), DIMENSION(:), ALLOCATABLE           :: zcent         !< zenith angle centers
        REAL(wp), DIMENSION(:), ALLOCATABLE           :: zbdry         !< zenith angle boundaries
        REAL(wp), DIMENSION(:), ALLOCATABLE           :: vffrac        !< view factor fractions for individual rays
        REAL(wp), DIMENSION(:), ALLOCATABLE           :: vffrac0       !< dtto (original values)
        REAL(wp), DIMENSION(:), ALLOCATABLE           :: ztransp       !< array of transparency in z steps
        INTEGER(iwp)                                  :: lowest_free_ray !< index into zdirs
        INTEGER(iwp), DIMENSION(:), ALLOCATABLE       :: itarget       !< face indices of detected obstacles
        INTEGER(iwp)                                  :: itarg0, itarg1

        INTEGER(iwp)                                  :: udim
        REAL(wp),     DIMENSION(:), ALLOCATABLE,TARGET:: csflt_l, pcsflt_l
        REAL(wp),     DIMENSION(:,:), POINTER         :: csflt, pcsflt
        INTEGER(iwp), DIMENSION(:), ALLOCATABLE,TARGET:: kcsflt_l,kpcsflt_l
        INTEGER(iwp), DIMENSION(:,:), POINTER         :: kcsflt,kpcsflt
        INTEGER(iwp), DIMENSION(:), ALLOCATABLE       :: icsflt,dcsflt,ipcsflt,dpcsflt
        REAL(wp), DIMENSION(3)                        :: uv
        LOGICAL                                       :: visible
        REAL(wp), DIMENSION(3)                        :: sa, ta          !< real coordinates z,y,x of source and target
        REAL(wp)                                      :: difvf           !< differential view factor
        REAL(wp)                                      :: transparency, rirrf, sqdist, svfsum
        INTEGER(iwp)                                  :: isurflt, isurfs, isurflt_prev
        INTEGER(idp)                                  :: ray_skip_maxdist, ray_skip_minval !< skipped raytracing counts
        INTEGER(iwp)                                  :: max_track_len !< maximum 2d track length
#if defined( __parallel )
        INTEGER(iwp), DIMENSION(:), ALLOCATABLE,TARGET:: nzterrl_l
        INTEGER(iwp), DIMENSION(:,:), POINTER         :: nzterrl
        INTEGER(iwp)                                  :: minfo
        REAL(wp), DIMENSION(:), POINTER, SAVE         :: lad_s_rma       !< fortran 1D pointer
        TYPE(c_ptr)                                   :: lad_s_rma_p     !< allocated c pointer
        INTEGER(kind=MPI_ADDRESS_KIND)                :: size_lad_rma
#endif
!
        INTEGER(iwp), DIMENSION(0:svfnorm_report_num) :: svfnorm_counts


!--     calculation of the SVF
        CALL location_message( 'calculating view factors for radiation interaction', 'start' )

!--     initialize variables and temporary arrays for calculation of svf and csf
        nsvfl  = 0
        ncsfl  = 0
        nsvfla = gasize
        msvf   = 1
        ALLOCATE( asvf1(nsvfla) )
        asvf => asvf1
        IF ( plant_canopy )  THEN
            ncsfla = gasize
            mcsf   = 1
            ALLOCATE( acsf1(ncsfla) )
            acsf => acsf1
        ENDIF
        nmrtf = 0
        IF ( mrt_nlevels > 0 )  THEN
           nmrtfa = gasize
           mmrtf = 1
           ALLOCATE ( amrtf1(nmrtfa) )
           amrtf => amrtf1
        ENDIF
        ray_skip_maxdist = 0
        ray_skip_minval = 0

!--     initialize temporary terrain and plant canopy height arrays (global 2D array!)
        ALLOCATE( nzterr(0:(nx+1)*(ny+1)-1) )
#if defined( __parallel )
        !ALLOCATE( nzterrl(nys:nyn,nxl:nxr) )
        ALLOCATE( nzterrl_l((nyn-nys+1)*(nxr-nxl+1)) )
        nzterrl(nys:nyn,nxl:nxr) => nzterrl_l(1:(nyn-nys+1)*(nxr-nxl+1))
        nzterrl = topo_top_ind(nys:nyn,nxl:nxr,0)
        CALL MPI_AllGather( nzterrl_l, nnx*nny, MPI_INTEGER, &
                            nzterr, nnx*nny, MPI_INTEGER, comm2d, ierr )
        IF ( ierr /= 0 ) THEN
            WRITE(9,*) 'Error MPI_AllGather1:', ierr, SIZE(nzterrl_l), nnx*nny, &
                       SIZE(nzterr), nnx*nny
            FLUSH(9)
        ENDIF
        DEALLOCATE(nzterrl_l)
#else
        nzterr = RESHAPE( topo_top_ind(nys:nyn,nxl:nxr,0), (/(nx+1)*(ny+1)/) )
#endif
        IF ( plant_canopy )  THEN
            ALLOCATE( plantt(0:(nx+1)*(ny+1)-1) )
            maxboxesg = nx + ny + nz_plant + 1
            max_track_len = nx + ny + 1
!--         temporary arrays storing values for csf calculation during raytracing
            ALLOCATE( boxes(3, maxboxesg) )
            ALLOCATE( crlens(maxboxesg) )

#if defined( __parallel )
            CALL MPI_AllGather( pct, nnx*nny, MPI_INTEGER, &
                                plantt, nnx*nny, MPI_INTEGER, comm2d, ierr )
            IF ( ierr /= 0 ) THEN
                WRITE(9,*) 'Error MPI_AllGather2:', ierr, SIZE(pct), nnx*nny, &
                           SIZE(plantt), nnx*nny
                FLUSH(9)
            ENDIF

!--         temporary arrays storing values for csf calculation during raytracing
            ALLOCATE( lad_ip(maxboxesg) )
            ALLOCATE( lad_disp(maxboxesg) )

            IF ( raytrace_mpi_rma )  THEN
                ALLOCATE( lad_s_ray(maxboxesg) )

                ! set conditions for RMA communication
                CALL MPI_Info_create(minfo, ierr)
                IF ( ierr /= 0 ) THEN
                    WRITE(9,*) 'Error MPI_Info_create2:', ierr
                    FLUSH(9)
                ENDIF
                CALL MPI_Info_set(minfo, 'accumulate_ordering', 'none', ierr)
                IF ( ierr /= 0 ) THEN
                    WRITE(9,*) 'Error MPI_Info_set5:', ierr
                    FLUSH(9)
                ENDIF
                CALL MPI_Info_set(minfo, 'accumulate_ops', 'same_op', ierr)
                IF ( ierr /= 0 ) THEN
                    WRITE(9,*) 'Error MPI_Info_set6:', ierr
                    FLUSH(9)
                ENDIF
                CALL MPI_Info_set(minfo, 'same_size', 'true', ierr)
                IF ( ierr /= 0 ) THEN
                    WRITE(9,*) 'Error MPI_Info_set7:', ierr
                    FLUSH(9)
                ENDIF
                CALL MPI_Info_set(minfo, 'same_disp_unit', 'true', ierr)
                IF ( ierr /= 0 ) THEN
                    WRITE(9,*) 'Error MPI_Info_set8:', ierr
                    FLUSH(9)
                ENDIF

!--             Allocate and initialize the MPI RMA window
!--             must be in accordance with allocation of lad_s in plant_canopy_model
!--             optimization of memory should be done
!--             Argument X of function STORAGE_SIZE(X) needs arbitrary REAL(wp) value, set to 1.0_wp for now
                size_lad_rma = STORAGE_SIZE(1.0_wp)/8*nnx*nny*nz_plant
                CALL MPI_Win_allocate(size_lad_rma, STORAGE_SIZE(1.0_wp)/8, minfo, comm2d, &
                                        lad_s_rma_p, win_lad, ierr)
                IF ( ierr /= 0 ) THEN
                    WRITE(9,*) 'Error MPI_Win_allocate2:', ierr, size_lad_rma, &
                                STORAGE_SIZE(1.0_wp)/8, win_lad
                    FLUSH(9)
                ENDIF
                CALL c_f_pointer(lad_s_rma_p, lad_s_rma, (/ nz_plant*nny*nnx /))
                sub_lad(nz_urban_b:nz_plant_t, nys:nyn, nxl:nxr) => lad_s_rma(1:nz_plant*nny*nnx)
            ELSE
                ALLOCATE(sub_lad(nz_urban_b:nz_plant_t, nys:nyn, nxl:nxr))
            ENDIF
#else
            plantt = RESHAPE( pct(nys:nyn,nxl:nxr), (/(nx+1)*(ny+1)/) )
            ALLOCATE(sub_lad(nz_urban_b:nz_plant_t, nys:nyn, nxl:nxr))
#endif
            plantt_max = MAXVAL(plantt)
            ALLOCATE( rt2_track(2, max_track_len), rt2_track_lad(nz_urban_b:plantt_max, max_track_len), &
                      rt2_track_dist(0:max_track_len), rt2_dist(plantt_max-nz_urban_b+2) )

            sub_lad(:,:,:) = 0._wp
            DO i = nxl, nxr
                DO j = nys, nyn
                    k = topo_top_ind(j,i,0)

                    sub_lad(k:nz_plant_t, j, i) = lad_s(0:nz_plant_t-k, j, i)
                ENDDO
            ENDDO

#if defined( __parallel )
            IF ( raytrace_mpi_rma )  THEN
                CALL MPI_Info_free(minfo, ierr)
                IF ( ierr /= 0 ) THEN
                    WRITE(9,*) 'Error MPI_Info_free2:', ierr
                    FLUSH(9)
                ENDIF
                CALL MPI_Win_lock_all(0, win_lad, ierr)
                IF ( ierr /= 0 ) THEN
                    WRITE(9,*) 'Error MPI_Win_lock_all1:', ierr, win_lad
                    FLUSH(9)
                ENDIF

            ELSE
                ALLOCATE( sub_lad_g(0:(nx+1)*(ny+1)*nz_plant-1) )
                CALL MPI_AllGather( sub_lad, nnx*nny*nz_plant, MPI_REAL, &
                                    sub_lad_g, nnx*nny*nz_plant, MPI_REAL, comm2d, ierr )
                IF ( ierr /= 0 ) THEN
                    WRITE(9,*) 'Error MPI_AllGather3:', ierr, SIZE(sub_lad), &
                                nnx*nny*nz_plant, SIZE(sub_lad_g), nnx*nny*nz_plant
                    FLUSH(9)
                ENDIF
            ENDIF
#endif
        ENDIF

!--     prepare the MPI_Win for collecting the surface indices
!--     from the reverse index arrays gridsurf from processors of target surfaces
#if defined( __parallel )
        IF ( rad_angular_discretization )  THEN
!
!--         raytrace_mpi_rma is asserted
            CALL MPI_Win_lock_all(0, win_gridsurf, ierr)
            IF ( ierr /= 0 ) THEN
                WRITE(9,*) 'Error MPI_Win_lock_all2:', ierr, win_gridsurf
                FLUSH(9)
            ENDIF
        ENDIF
#endif


        !--Directions opposite to face normals are not even calculated,
        !--they must be preset to 0
        !--
        dsitrans(:,:) = 0._wp

        DO isurflt = 1, nsurfl
!--         determine face centers
            td = surfl(id, isurflt)
            ta = (/ REAL(surfl(iz, isurflt), wp) - 0.5_wp * kdir(td),  &
                      REAL(surfl(iy, isurflt), wp) - 0.5_wp * jdir(td),  &
                      REAL(surfl(ix, isurflt), wp) - 0.5_wp * idir(td)  /)

            !--Calculate sky view factor and raytrace DSI paths
            skyvf(isurflt) = 0._wp
            skyvft(isurflt) = 0._wp

            !--Select a proper half-sphere for 2D raytracing
            SELECT CASE ( td )
               CASE ( iup_u, iup_l )
                  az0 = 0._wp
                  naz = raytrace_discrete_azims
                  azs = 2._wp * pi / REAL(naz, wp)
                  zn0 = 0._wp
                  nzn = raytrace_discrete_elevs / 2
                  zns = pi / 2._wp / REAL(nzn, wp)
               CASE ( isouth_u, isouth_l )
                  az0 = pi / 2._wp
                  naz = raytrace_discrete_azims / 2
                  azs = pi / REAL(naz, wp)
                  zn0 = 0._wp
                  nzn = raytrace_discrete_elevs
                  zns = pi / REAL(nzn, wp)
               CASE ( inorth_u, inorth_l )
                  az0 = - pi / 2._wp
                  naz = raytrace_discrete_azims / 2
                  azs = pi / REAL(naz, wp)
                  zn0 = 0._wp
                  nzn = raytrace_discrete_elevs
                  zns = pi / REAL(nzn, wp)
               CASE ( iwest_u, iwest_l )
                  az0 = pi
                  naz = raytrace_discrete_azims / 2
                  azs = pi / REAL(naz, wp)
                  zn0 = 0._wp
                  nzn = raytrace_discrete_elevs
                  zns = pi / REAL(nzn, wp)
               CASE ( ieast_u, ieast_l )
                  az0 = 0._wp
                  naz = raytrace_discrete_azims / 2
                  azs = pi / REAL(naz, wp)
                  zn0 = 0._wp
                  nzn = raytrace_discrete_elevs
                  zns = pi / REAL(nzn, wp)
               CASE DEFAULT
                  WRITE(message_string, *) 'ERROR: the surface type ', td,     &
                                           ' is not supported for calculating',&
                                           ' SVF'
                  CALL message( 'radiation_calc_svf', 'PA0488', 1, 2, 0, 6, 0 )
            END SELECT

            ALLOCATE ( zdirs(1:nzn), zcent(1:nzn), zbdry(0:nzn), vffrac(1:nzn*naz), &
                       ztransp(1:nzn*naz), itarget(1:nzn*naz) )   !FIXME allocate itarget only
                                                                  !in case of rad_angular_discretization

            itarg0 = 1
            itarg1 = nzn
            zcent(:) = (/( zn0+(REAL(izn,wp)-.5_wp)*zns, izn=1, nzn )/)
            zbdry(:) = (/( zn0+REAL(izn,wp)*zns, izn=0, nzn )/)
            IF ( td == iup_u  .OR.  td == iup_l )  THEN
               vffrac(1:nzn) = (COS(2 * zbdry(0:nzn-1)) - COS(2 * zbdry(1:nzn))) / 2._wp / REAL(naz, wp)
!
!--            For horizontal target, vf fractions are constant per azimuth
               DO iaz = 1, naz-1
                  vffrac(iaz*nzn+1:(iaz+1)*nzn) = vffrac(1:nzn)
               ENDDO
!--            sum of whole vffrac equals 1, verified
            ENDIF
!
!--         Calculate sky-view factor and direct solar visibility using 2D raytracing
            DO iaz = 1, naz
               azmid = az0 + (REAL(iaz, wp) - .5_wp) * azs
               IF ( td /= iup_u  .AND.  td /= iup_l )  THEN
                  az2 = REAL(iaz, wp) * azs - pi/2._wp
                  az1 = az2 - azs
                  !TODO precalculate after 1st line
                  vffrac(itarg0:itarg1) = (SIN(az2) - SIN(az1))               &
                              * (zbdry(1:nzn) - zbdry(0:nzn-1)                &
                                 + SIN(zbdry(0:nzn-1))*COS(zbdry(0:nzn-1))    &
                                 - SIN(zbdry(1:nzn))*COS(zbdry(1:nzn)))       &
                              / (2._wp * pi)
!--               sum of whole vffrac equals 1, verified
               ENDIF
               yxdir(:) = (/ COS(azmid) / dy, SIN(azmid) / dx /)
               yxlen = SQRT(SUM(yxdir(:)**2))
               zdirs(:) = COS(zcent(:)) / (dz(1) * yxlen * SIN(zcent(:)))
               yxdir(:) = yxdir(:) / yxlen

               CALL raytrace_2d(ta, yxdir, nzn, zdirs,                        &
                                    surfstart(myid) + isurflt, facearea(td),  &
                                    vffrac(itarg0:itarg1), .TRUE., .TRUE.,    &
                                    .FALSE., lowest_free_ray,                 &
                                    ztransp(itarg0:itarg1),                   &
                                    itarget(itarg0:itarg1))

               skyvf(isurflt) = skyvf(isurflt) + &
                                SUM(vffrac(itarg0:itarg0+lowest_free_ray-1))
               skyvft(isurflt) = skyvft(isurflt) + &
                                 SUM(ztransp(itarg0:itarg0+lowest_free_ray-1) &
                                             * vffrac(itarg0:itarg0+lowest_free_ray-1))

!--            Save direct solar transparency
               j = MODULO(NINT(azmid/                                          &
                               (2._wp*pi)*raytrace_discrete_azims-.5_wp, iwp), &
                          raytrace_discrete_azims)

               DO k = 1, raytrace_discrete_elevs/2
                  i = dsidir_rev(k-1, j)
                  IF ( i /= -1  .AND.  k <= lowest_free_ray )  &
                     dsitrans(isurflt, i) = ztransp(itarg0+k-1)
               ENDDO

!
!--            Advance itarget indices
               itarg0 = itarg1 + 1
               itarg1 = itarg1 + nzn
            ENDDO

            IF ( rad_angular_discretization )  THEN
!--            sort itarget by face id
               CALL quicksort_itarget(itarget,vffrac,ztransp,1,nzn*naz)
!
!--            For aggregation, we need fractions multiplied by transmissivities
               ztransp(:) = vffrac(:) * ztransp(:)
!
!--            find the first valid position
               itarg0 = 1
               DO WHILE ( itarg0 <= nzn*naz )
                  IF ( itarget(itarg0) /= -1 )  EXIT
                  itarg0 = itarg0 + 1
               ENDDO

               DO  i = itarg0, nzn*naz
!
!--               For duplicate values, only sum up vf fraction value
                  IF ( i < nzn*naz )  THEN
                     IF ( itarget(i+1) == itarget(i) )  THEN
                        vffrac(i+1) = vffrac(i+1) + vffrac(i)
                        ztransp(i+1) = ztransp(i+1) + ztransp(i)
                        CYCLE
                     ENDIF
                  ENDIF
!
!--               write to the svf array
                  nsvfl = nsvfl + 1
!--               check dimmension of asvf array and enlarge it if needed
                  IF ( nsvfla < nsvfl )  THEN
                     k = CEILING(REAL(nsvfla, kind=wp) * grow_factor)
                     IF ( msvf == 0 )  THEN
                        msvf = 1
                        ALLOCATE( asvf1(k) )
                        asvf => asvf1
                        asvf1(1:nsvfla) = asvf2
                        DEALLOCATE( asvf2 )
                     ELSE
                        msvf = 0
                        ALLOCATE( asvf2(k) )
                        asvf => asvf2
                        asvf2(1:nsvfla) = asvf1
                        DEALLOCATE( asvf1 )
                     ENDIF

                     IF ( debug_output )  THEN
                        WRITE( debug_string, '(A,3I12)' ) 'Grow asvf:', nsvfl, nsvfla, k
                        CALL debug_message( debug_string, 'info' )
                     ENDIF

                     nsvfla = k
                  ENDIF
!--               write svf values into the array
                  asvf(nsvfl)%isurflt = isurflt
                  asvf(nsvfl)%isurfs = itarget(i)
                  asvf(nsvfl)%rsvf = vffrac(i)
                  asvf(nsvfl)%rtransp = ztransp(i) / vffrac(i)
               END DO

            ENDIF ! rad_angular_discretization

            DEALLOCATE ( zdirs, zcent, zbdry, vffrac, ztransp, itarget ) !FIXME itarget shall be allocated only
                                                                  !in case of rad_angular_discretization
!
!--         Following calculations only required for surface_reflections
            IF ( surface_reflections  .AND.  .NOT. rad_angular_discretization )  THEN

               DO  isurfs = 1, nsurf
                  IF ( .NOT.  surface_facing(surfl(ix, isurflt), surfl(iy, isurflt), &
                     surfl(iz, isurflt), surfl(id, isurflt), &
                     surf(ix, isurfs), surf(iy, isurfs), &
                     surf(iz, isurfs), surf(id, isurfs)) )  THEN
                     CYCLE
                  ENDIF

                  sd = surf(id, isurfs)
                  sa = (/ REAL(surf(iz, isurfs), wp) - 0.5_wp * kdir(sd),  &
                          REAL(surf(iy, isurfs), wp) - 0.5_wp * jdir(sd),  &
                          REAL(surf(ix, isurfs), wp) - 0.5_wp * idir(sd)  /)

!--               unit vector source -> target
                  uv = (/ (ta(1)-sa(1))*dz(1), (ta(2)-sa(2))*dy, (ta(3)-sa(3))*dx /)
                  sqdist = SUM(uv(:)**2)
                  uv = uv / SQRT(sqdist)

!--               reject raytracing above max distance
                  IF ( SQRT(sqdist) > max_raytracing_dist ) THEN
                     ray_skip_maxdist = ray_skip_maxdist + 1
                     CYCLE
                  ENDIF

                  difvf = dot_product((/ kdir(sd), jdir(sd), idir(sd) /), uv) & ! cosine of source normal and direction
                      * dot_product((/ kdir(td), jdir(td), idir(td) /), -uv) &  ! cosine of target normal and reverse direction
                      / (pi * sqdist) ! square of distance between centers
!
!--               irradiance factor (our unshaded shape view factor) = view factor per differential target area * source area
                  rirrf = difvf * facearea(sd)

!--               reject raytracing for potentially too small view factor values
                  IF ( rirrf < min_irrf_value ) THEN
                      ray_skip_minval = ray_skip_minval + 1
                      CYCLE
                  ENDIF

!--               raytrace + process plant canopy sinks within
                  CALL raytrace(sa, ta, isurfs, difvf, facearea(td), .TRUE., &
                                visible, transparency)

                  IF ( .NOT.  visible ) CYCLE
                 ! rsvf = rirrf * transparency

!--               write to the svf array
                  nsvfl = nsvfl + 1
!--               check dimmension of asvf array and enlarge it if needed
                  IF ( nsvfla < nsvfl )  THEN
                     k = CEILING(REAL(nsvfla, kind=wp) * grow_factor)
                     IF ( msvf == 0 )  THEN
                        msvf = 1
                        ALLOCATE( asvf1(k) )
                        asvf => asvf1
                        asvf1(1:nsvfla) = asvf2
                        DEALLOCATE( asvf2 )
                     ELSE
                        msvf = 0
                        ALLOCATE( asvf2(k) )
                        asvf => asvf2
                        asvf2(1:nsvfla) = asvf1
                        DEALLOCATE( asvf1 )
                     ENDIF

                     IF ( debug_output )  THEN
                        WRITE( debug_string, '(A,3I12)' ) 'Grow asvf:', nsvfl, nsvfla, k
                        CALL debug_message( debug_string, 'info' )
                     ENDIF

                     nsvfla = k
                  ENDIF
!--               write svf values into the array
                  asvf(nsvfl)%isurflt = isurflt
                  asvf(nsvfl)%isurfs = isurfs
                  asvf(nsvfl)%rsvf = rirrf !we postopne multiplication by transparency
                  asvf(nsvfl)%rtransp = transparency !a.k.a. Direct Irradiance Factor
               ENDDO
            ENDIF
        ENDDO

!--
!--     Raytrace to canopy boxes to fill dsitransc TODO optimize
        dsitransc(:,:) = 0._wp
        az0 = 0._wp
        naz = raytrace_discrete_azims
        azs = 2._wp * pi / REAL(naz, wp)
        zn0 = 0._wp
        nzn = raytrace_discrete_elevs / 2
        zns = pi / 2._wp / REAL(nzn, wp)
        ALLOCATE ( zdirs(1:nzn), zcent(1:nzn), vffrac(1:nzn), ztransp(1:nzn), &
               itarget(1:nzn) )
        zcent(:) = (/( zn0+(REAL(izn,wp)-.5_wp)*zns, izn=1, nzn )/)
        vffrac(:) = 0._wp

        DO  ipcgb = 1, npcbl
           ta = (/ REAL(pcbl(iz, ipcgb), wp),  &
                   REAL(pcbl(iy, ipcgb), wp),  &
                   REAL(pcbl(ix, ipcgb), wp) /)
!--        Calculate direct solar visibility using 2D raytracing
           DO  iaz = 1, naz
              azmid = az0 + (REAL(iaz, wp) - .5_wp) * azs
              yxdir(:) = (/ COS(azmid) / dy, SIN(azmid) / dx /)
              yxlen = SQRT(SUM(yxdir(:)**2))
              zdirs(:) = COS(zcent(:)) / (dz(1) * yxlen * SIN(zcent(:)))
              yxdir(:) = yxdir(:) / yxlen
              CALL raytrace_2d(ta, yxdir, nzn, zdirs,                                &
                                   -999, -999._wp, vffrac, .FALSE., .FALSE., .TRUE., &
                                   lowest_free_ray, ztransp, itarget)

!--           Save direct solar transparency
              j = MODULO(NINT(azmid/                                         &
                             (2._wp*pi)*raytrace_discrete_azims-.5_wp, iwp), &
                         raytrace_discrete_azims)
              DO  k = 1, raytrace_discrete_elevs/2
                 i = dsidir_rev(k-1, j)
                 IF ( i /= -1  .AND.  k <= lowest_free_ray ) &
                    dsitransc(ipcgb, i) = ztransp(k)
              ENDDO
           ENDDO
        ENDDO
        DEALLOCATE ( zdirs, zcent, vffrac, ztransp, itarget )
!--
!--     Raytrace to MRT boxes
        IF ( nmrtbl > 0 )  THEN
           mrtdsit(:,:) = 0._wp
           mrtsky(:) = 0._wp
           mrtskyt(:) = 0._wp
           az0 = 0._wp
           naz = raytrace_discrete_azims
           azs = 2._wp * pi / REAL(naz, wp)
           zn0 = 0._wp
           nzn = raytrace_discrete_elevs
           zns = pi / REAL(nzn, wp)
           ALLOCATE ( zdirs(1:nzn), zcent(1:nzn), zbdry(0:nzn), vffrac(1:nzn*naz), vffrac0(1:nzn), &
                      ztransp(1:nzn*naz), itarget(1:nzn*naz) )   !FIXME allocate itarget only
                                                                 !in case of rad_angular_discretization

           zcent(:) = (/( zn0+(REAL(izn,wp)-.5_wp)*zns, izn=1, nzn )/)
           zbdry(:) = (/( zn0+REAL(izn,wp)*zns, izn=0, nzn )/)
           vffrac0(:) = (COS(zbdry(0:nzn-1)) - COS(zbdry(1:nzn))) / 2._wp / REAL(naz, wp)
!
!--        Modify direction weights to simulate human body (lower weight for
!--        irradiance from zenith, higher from sides) depending on selection.
!--        For mrt_geom=0, no weighting is done (simulates spherical globe
!--        thermometer).
           SELECT CASE ( mrt_geom )

           CASE ( 1 )
              vffrac0(:) = vffrac0(:) * MAX(0._wp, SIN(zcent(:))*mrt_geom_params(2) &
                                                   + COS(zcent(:))*mrt_geom_params(1))
              vffrac0(:) = vffrac0(:) / (SUM(vffrac0) * REAL(naz, wp))

           CASE ( 2 )
              vffrac0(:) = vffrac0(:)                                          &
                           * SQRT( ( mrt_geom_params(1) * COS(zcent(:)) ) ** 2 &
                                 + ( mrt_geom_params(2) * SIN(zcent(:)) ) ** 2 )
              vffrac0(:) = vffrac0(:) / (SUM(vffrac0) * REAL(naz, wp))

           END SELECT

           DO  imrt = 1, nmrtbl
              ta = (/ REAL(mrtbl(iz, imrt), wp),  &
                      REAL(mrtbl(iy, imrt), wp),  &
                      REAL(mrtbl(ix, imrt), wp) /)
!
!--           vf fractions are constant per azimuth
              DO iaz = 0, naz-1
                 vffrac(iaz*nzn+1:(iaz+1)*nzn) = vffrac0(:)
              ENDDO
!--           sum of whole vffrac equals 1, verified
              itarg0 = 1
              itarg1 = nzn
!
!--           Calculate sky-view factor and direct solar visibility using 2D raytracing
              DO  iaz = 1, naz
                 azmid = az0 + (REAL(iaz, wp) - .5_wp) * azs
                 yxdir(:) = (/ COS(azmid) / dy, SIN(azmid) / dx /)
                 yxlen = SQRT(SUM(yxdir(:)**2))
                 zdirs(:) = COS(zcent(:)) / (dz(1) * yxlen * SIN(zcent(:)))
                 yxdir(:) = yxdir(:) / yxlen

                 CALL raytrace_2d(ta, yxdir, nzn, zdirs,                         &
                                  -999, -999._wp, vffrac(itarg0:itarg1), .TRUE., &
                                  .FALSE., .TRUE., lowest_free_ray,              &
                                  ztransp(itarg0:itarg1),                        &
                                  itarget(itarg0:itarg1))

!--              Sky view factors for MRT
                 mrtsky(imrt) = mrtsky(imrt) + &
                                  SUM(vffrac(itarg0:itarg0+lowest_free_ray-1))
                 mrtskyt(imrt) = mrtskyt(imrt) + &
                                   SUM(ztransp(itarg0:itarg0+lowest_free_ray-1) &
                                               * vffrac(itarg0:itarg0+lowest_free_ray-1))
!--              Direct solar transparency for MRT
                 j = MODULO(NINT(azmid/                                         &
                                (2._wp*pi)*raytrace_discrete_azims-.5_wp, iwp), &
                            raytrace_discrete_azims)
                 DO  k = 1, raytrace_discrete_elevs/2
                    i = dsidir_rev(k-1, j)
                    IF ( i /= -1  .AND.  k <= lowest_free_ray ) &
                       mrtdsit(imrt, i) = ztransp(itarg0+k-1)
                 ENDDO
!
!--              Advance itarget indices
                 itarg0 = itarg1 + 1
                 itarg1 = itarg1 + nzn
              ENDDO

!--           sort itarget by face id
              CALL quicksort_itarget(itarget,vffrac,ztransp,1,nzn*naz)
!
!--           find the first valid position
              itarg0 = 1
              DO WHILE ( itarg0 <= nzn*naz )
                 IF ( itarget(itarg0) /= -1 )  EXIT
                 itarg0 = itarg0 + 1
              ENDDO

              DO  i = itarg0, nzn*naz
!
!--              For duplicate values, only sum up vf fraction value
                 IF ( i < nzn*naz )  THEN
                    IF ( itarget(i+1) == itarget(i) )  THEN
                       vffrac(i+1) = vffrac(i+1) + vffrac(i)
                       CYCLE
                    ENDIF
                 ENDIF
!
!--              write to the mrtf array
                 nmrtf = nmrtf + 1
!--              check dimmension of mrtf array and enlarge it if needed
                 IF ( nmrtfa < nmrtf )  THEN
                    k = CEILING(REAL(nmrtfa, kind=wp) * grow_factor)
                    IF ( mmrtf == 0 )  THEN
                       mmrtf = 1
                       ALLOCATE( amrtf1(k) )
                       amrtf => amrtf1
                       amrtf1(1:nmrtfa) = amrtf2
                       DEALLOCATE( amrtf2 )
                    ELSE
                       mmrtf = 0
                       ALLOCATE( amrtf2(k) )
                       amrtf => amrtf2
                       amrtf2(1:nmrtfa) = amrtf1
                       DEALLOCATE( amrtf1 )
                    ENDIF

                    IF ( debug_output )  THEN
                       WRITE( debug_string, '(A,3I12)' ) 'Grow amrtf:', nmrtf, nmrtfa, k
                       CALL debug_message( debug_string, 'info' )
                    ENDIF

                    nmrtfa = k
                 ENDIF
!--              write mrtf values into the array
                 amrtf(nmrtf)%isurflt = imrt
                 amrtf(nmrtf)%isurfs = itarget(i)
                 amrtf(nmrtf)%rsvf = vffrac(i)
                 amrtf(nmrtf)%rtransp = ztransp(i)
              ENDDO ! itarg

           ENDDO ! imrt
           DEALLOCATE ( zdirs, zcent, zbdry, vffrac, vffrac0, ztransp, itarget )
!
!--        Move MRT factors to final arrays
           ALLOCATE ( mrtf(nmrtf), mrtft(nmrtf), mrtfsurf(2,nmrtf) )
           DO  imrtf = 1, nmrtf
              mrtf(imrtf) = amrtf(imrtf)%rsvf
              mrtft(imrtf) = amrtf(imrtf)%rsvf * amrtf(imrtf)%rtransp
              mrtfsurf(:,imrtf) = (/amrtf(imrtf)%isurflt, amrtf(imrtf)%isurfs /)
           ENDDO
           IF ( ALLOCATED(amrtf1) )  DEALLOCATE( amrtf1 )
           IF ( ALLOCATED(amrtf2) )  DEALLOCATE( amrtf2 )
        ENDIF ! nmrtbl > 0

        IF ( rad_angular_discretization )  THEN
#if defined( __parallel )
!--        finalize MPI_RMA communication established to get global index of the surface from grid indices
!--        flush all MPI window pending requests
           CALL MPI_Win_flush_all(win_gridsurf, ierr)
           IF ( ierr /= 0 ) THEN
               WRITE(9,*) 'Error MPI_Win_flush_all1:', ierr, win_gridsurf
               FLUSH(9)
           ENDIF
!--        unlock MPI window
           CALL MPI_Win_unlock_all(win_gridsurf, ierr)
           IF ( ierr /= 0 ) THEN
               WRITE(9,*) 'Error MPI_Win_unlock_all1:', ierr, win_gridsurf
               FLUSH(9)
           ENDIF
!--        free MPI window
           CALL MPI_Win_free(win_gridsurf, ierr)
           IF ( ierr /= 0 ) THEN
               WRITE(9,*) 'Error MPI_Win_free1:', ierr, win_gridsurf
               FLUSH(9)
           ENDIF
#else
           DEALLOCATE ( gridsurf )
#endif
        ENDIF

        IF ( debug_output )  CALL debug_message( 'waiting for completion of SVF and CSF calculation in all processes', 'info' )

!--     deallocate temporary global arrays
        DEALLOCATE(nzterr)

        IF ( plant_canopy )  THEN
!--         finalize mpi_rma communication and deallocate temporary arrays
#if defined( __parallel )
            IF ( raytrace_mpi_rma )  THEN
                CALL MPI_Win_flush_all(win_lad, ierr)
                IF ( ierr /= 0 ) THEN
                    WRITE(9,*) 'Error MPI_Win_flush_all2:', ierr, win_lad
                    FLUSH(9)
                ENDIF
!--             unlock MPI window
                CALL MPI_Win_unlock_all(win_lad, ierr)
                IF ( ierr /= 0 ) THEN
                    WRITE(9,*) 'Error MPI_Win_unlock_all2:', ierr, win_lad
                    FLUSH(9)
                ENDIF
!--             free MPI window
                CALL MPI_Win_free(win_lad, ierr)
                IF ( ierr /= 0 ) THEN
                    WRITE(9,*) 'Error MPI_Win_free2:', ierr, win_lad
                    FLUSH(9)
                ENDIF
!--             deallocate temporary arrays storing values for csf calculation during raytracing
                DEALLOCATE( lad_s_ray )
!--             sub_lad is the pointer to lad_s_rma in case of raytrace_mpi_rma
!--             and must not be deallocated here
            ELSE
                DEALLOCATE(sub_lad)
                DEALLOCATE(sub_lad_g)
            ENDIF
#else
            DEALLOCATE(sub_lad)
#endif
            DEALLOCATE( boxes )
            DEALLOCATE( crlens )
            DEALLOCATE( plantt )
            DEALLOCATE( rt2_track, rt2_track_lad, rt2_track_dist, rt2_dist )
        ENDIF

        IF ( debug_output )  CALL debug_message( 'calculation of the complete SVF array', 'info' )

        IF ( rad_angular_discretization )  THEN
           IF ( debug_output )  THEN
              WRITE( debug_string, '("Load ",I0," SVFs from the structure array to plain arrays")' ) nsvfl
              CALL debug_message( debug_string, 'info' )
           ENDIF
           ALLOCATE( svf(ndsvf,nsvfl) )
           ALLOCATE( svfsurf(idsvf,nsvfl) )

           DO isvf = 1, nsvfl
               svf(:, isvf) = (/ asvf(isvf)%rsvf, asvf(isvf)%rtransp /)
               svfsurf(:, isvf) = (/ asvf(isvf)%isurflt, asvf(isvf)%isurfs /)
           ENDDO
        ELSE
           IF ( debug_output )  CALL debug_message( 'Start SVF sort', 'info' )
!--        sort svf ( a version of quicksort )
           CALL quicksort_svf(asvf,1,nsvfl)

           !< load svf from the structure array to plain arrays
           IF ( debug_output )  THEN
              WRITE( debug_string, '("Load ",I0," SVFs from the structure array to plain arrays")' ) nsvfl
              CALL debug_message( debug_string, 'info' )
           ENDIF
           ALLOCATE( svf(ndsvf,nsvfl) )
           ALLOCATE( svfsurf(idsvf,nsvfl) )
           svfnorm_counts(:) = 0._wp
           isurflt_prev = -1
           ksvf = 1
           svfsum = 0._wp
           DO isvf = 1, nsvfl
!--            normalize svf per target face
               IF ( asvf(ksvf)%isurflt /= isurflt_prev )  THEN
                   IF ( isurflt_prev /= -1  .AND.  svfsum /= 0._wp )  THEN
                       !< update histogram of logged svf normalization values
                       i = searchsorted(svfnorm_report_thresh, svfsum / (1._wp-skyvf(isurflt_prev)))
                       svfnorm_counts(i) = svfnorm_counts(i) + 1

                       svf(1, isvf_surflt:isvf-1) = svf(1, isvf_surflt:isvf-1) / svfsum * (1._wp-skyvf(isurflt_prev))
                   ENDIF
                   isurflt_prev = asvf(ksvf)%isurflt
                   isvf_surflt = isvf
                   svfsum = asvf(ksvf)%rsvf !?? / asvf(ksvf)%rtransp
               ELSE
                   svfsum = svfsum + asvf(ksvf)%rsvf !?? / asvf(ksvf)%rtransp
               ENDIF

               svf(:, isvf) = (/ asvf(ksvf)%rsvf, asvf(ksvf)%rtransp /)
               svfsurf(:, isvf) = (/ asvf(ksvf)%isurflt, asvf(ksvf)%isurfs /)

!--            next element
               ksvf = ksvf + 1
           ENDDO

           IF ( isurflt_prev /= -1  .AND.  svfsum /= 0._wp )  THEN
               i = searchsorted(svfnorm_report_thresh, svfsum / (1._wp-skyvf(isurflt_prev)))
               svfnorm_counts(i) = svfnorm_counts(i) + 1

               svf(1, isvf_surflt:nsvfl) = svf(1, isvf_surflt:nsvfl) / svfsum * (1._wp-skyvf(isurflt_prev))
           ENDIF
           WRITE(9, *) 'SVF normalization histogram:', svfnorm_counts,  &
                       'on thresholds:', svfnorm_report_thresh(1:svfnorm_report_num), '(val < thresh <= val)'
           !TODO we should be able to deallocate skyvf, from now on we only need skyvft
        ENDIF ! rad_angular_discretization

!--     deallocate temporary asvf array
!--     DEALLOCATE(asvf) - ifort has a problem with deallocation of allocatable target
!--     via pointing pointer - we need to test original targets
        IF ( ALLOCATED(asvf1) )  THEN
            DEALLOCATE(asvf1)
        ENDIF
        IF ( ALLOCATED(asvf2) )  THEN
            DEALLOCATE(asvf2)
        ENDIF

        npcsfl = 0
        IF ( plant_canopy )  THEN

            IF ( debug_output )  CALL debug_message( 'Calculation of the complete CSF array', 'info' )
!--         sort and merge csf for the last time, keeping the array size to minimum
            CALL merge_and_grow_csf(-1)

!--         aggregate csb among processors
!--         allocate necessary arrays
            udim = max(ncsfl,1)
            ALLOCATE( csflt_l(ndcsf*udim) )
            csflt(1:ndcsf,1:udim) => csflt_l(1:ndcsf*udim)
            ALLOCATE( kcsflt_l(kdcsf*udim) )
            kcsflt(1:kdcsf,1:udim) => kcsflt_l(1:kdcsf*udim)
            ALLOCATE( icsflt(0:numprocs-1) )
            ALLOCATE( dcsflt(0:numprocs-1) )
            ALLOCATE( ipcsflt(0:numprocs-1) )
            ALLOCATE( dpcsflt(0:numprocs-1) )

!--         fill out arrays of csf values and
!--         arrays of number of elements and displacements
!--         for particular precessors
            icsflt = 0
            dcsflt = 0
            ip = -1
            j = -1
            d = 0
            DO kcsf = 1, ncsfl
                j = j+1
                IF ( acsf(kcsf)%ip /= ip )  THEN
!--                 new block of the processor
!--                 number of elements of previous block
                    IF ( ip>=0) icsflt(ip) = j
                    d = d+j
!--                 blank blocks
                    DO jp = ip+1, acsf(kcsf)%ip-1
!--                     number of elements is zero, displacement is equal to previous
                        icsflt(jp) = 0
                        dcsflt(jp) = d
                    ENDDO
!--                 the actual block
                    ip = acsf(kcsf)%ip
                    dcsflt(ip) = d
                    j = 0
                ENDIF
                csflt(1,kcsf) = acsf(kcsf)%rcvf
!--             fill out integer values of itz,ity,itx,isurfs
                kcsflt(1,kcsf) = acsf(kcsf)%itz
                kcsflt(2,kcsf) = acsf(kcsf)%ity
                kcsflt(3,kcsf) = acsf(kcsf)%itx
                kcsflt(4,kcsf) = acsf(kcsf)%isurfs
            ENDDO
!--         last blank blocks at the end of array
            j = j+1
            IF ( ip>=0 ) icsflt(ip) = j
            d = d+j
            DO jp = ip+1, numprocs-1
!--             number of elements is zero, displacement is equal to previous
                icsflt(jp) = 0
                dcsflt(jp) = d
            ENDDO

!--         deallocate temporary acsf array
!--         DEALLOCATE(acsf) - ifort has a problem with deallocation of allocatable target
!--         via pointing pointer - we need to test original targets
            IF ( ALLOCATED(acsf1) )  THEN
                DEALLOCATE(acsf1)
            ENDIF
            IF ( ALLOCATED(acsf2) )  THEN
                DEALLOCATE(acsf2)
            ENDIF

#if defined( __parallel )
!--         scatter and gather the number of elements to and from all processor
!--         and calculate displacements
            IF ( debug_output )  CALL debug_message( 'Scatter and gather the number of elements to and from all processor', 'info' )

            CALL MPI_AlltoAll(icsflt,1,MPI_INTEGER,ipcsflt,1,MPI_INTEGER,comm2d, ierr)

            IF ( ierr /= 0 ) THEN
                WRITE(9,*) 'Error MPI_AlltoAll1:', ierr, SIZE(icsflt), SIZE(ipcsflt)
                FLUSH(9)
            ENDIF

            npcsfl = SUM(ipcsflt)
            d = 0
            DO i = 0, numprocs-1
                dpcsflt(i) = d
                d = d + ipcsflt(i)
            ENDDO

!--         exchange csf fields between processors
            IF ( debug_output )  CALL debug_message( 'Exchange csf fields between processors', 'info' )
            udim = max(npcsfl,1)
            ALLOCATE( pcsflt_l(ndcsf*udim) )
            pcsflt(1:ndcsf,1:udim) => pcsflt_l(1:ndcsf*udim)
            ALLOCATE( kpcsflt_l(kdcsf*udim) )
            kpcsflt(1:kdcsf,1:udim) => kpcsflt_l(1:kdcsf*udim)
            CALL MPI_AlltoAllv(csflt_l, ndcsf*icsflt, ndcsf*dcsflt, MPI_REAL, &
                pcsflt_l, ndcsf*ipcsflt, ndcsf*dpcsflt, MPI_REAL, comm2d, ierr)
            IF ( ierr /= 0 ) THEN
                WRITE(9,*) 'Error MPI_AlltoAllv1:', ierr, SIZE(ipcsflt), ndcsf*icsflt, &
                            ndcsf*dcsflt, SIZE(pcsflt_l),ndcsf*ipcsflt, ndcsf*dpcsflt
                FLUSH(9)
            ENDIF

            CALL MPI_AlltoAllv(kcsflt_l, kdcsf*icsflt, kdcsf*dcsflt, MPI_INTEGER, &
                kpcsflt_l, kdcsf*ipcsflt, kdcsf*dpcsflt, MPI_INTEGER, comm2d, ierr)
            IF ( ierr /= 0 ) THEN
                WRITE(9,*) 'Error MPI_AlltoAllv2:', ierr, SIZE(kcsflt_l),kdcsf*icsflt, &
                           kdcsf*dcsflt, SIZE(kpcsflt_l), kdcsf*ipcsflt, kdcsf*dpcsflt
                FLUSH(9)
            ENDIF

#else
            npcsfl = ncsfl
            ALLOCATE( pcsflt(ndcsf,max(npcsfl,ndcsf)) )
            ALLOCATE( kpcsflt(kdcsf,max(npcsfl,kdcsf)) )
            pcsflt = csflt
            kpcsflt = kcsflt
#endif

!--         deallocate temporary arrays
            DEALLOCATE( csflt_l )
            DEALLOCATE( kcsflt_l )
            DEALLOCATE( icsflt )
            DEALLOCATE( dcsflt )
            DEALLOCATE( ipcsflt )
            DEALLOCATE( dpcsflt )

!--         sort csf ( a version of quicksort )
            IF ( debug_output )  CALL debug_message( 'Sort csf', 'info' )
            CALL quicksort_csf2(kpcsflt, pcsflt, 1, npcsfl)

!--         aggregate canopy sink factor records with identical box & source
!--         againg across all values from all processors
            IF ( debug_output )  CALL debug_message( 'Aggregate canopy sink factor records with identical box', 'info' )

            IF ( npcsfl > 0 )  THEN
                icsf = 1 !< reading index
                kcsf = 1 !< writing index
                DO WHILE (icsf < npcsfl)
!--                 here kpcsf(kcsf) already has values from kpcsf(icsf)
                    IF ( kpcsflt(3,icsf) == kpcsflt(3,icsf+1)  .AND.  &
                         kpcsflt(2,icsf) == kpcsflt(2,icsf+1)  .AND.  &
                         kpcsflt(1,icsf) == kpcsflt(1,icsf+1)  .AND.  &
                         kpcsflt(4,icsf) == kpcsflt(4,icsf+1) )  THEN

                        pcsflt(1,kcsf) = pcsflt(1,kcsf) + pcsflt(1,icsf+1)

!--                     advance reading index, keep writing index
                        icsf = icsf + 1
                    ELSE
!--                     not identical, just advance and copy
                        icsf = icsf + 1
                        kcsf = kcsf + 1
                        kpcsflt(:,kcsf) = kpcsflt(:,icsf)
                        pcsflt(:,kcsf) = pcsflt(:,icsf)
                    ENDIF
                ENDDO
!--             last written item is now also the last item in valid part of array
                npcsfl = kcsf
            ENDIF

            ncsfl = npcsfl
            IF ( ncsfl > 0 )  THEN
                ALLOCATE( csf(ndcsf,ncsfl) )
                ALLOCATE( csfsurf(idcsf,ncsfl) )
                DO icsf = 1, ncsfl
                    csf(:,icsf) = pcsflt(:,icsf)
                    csfsurf(1,icsf) =  gridpcbl(kpcsflt(1,icsf),kpcsflt(2,icsf),kpcsflt(3,icsf))
                    csfsurf(2,icsf) =  kpcsflt(4,icsf)
                ENDDO
            ENDIF

!--         deallocation of temporary arrays
            IF ( npcbl > 0 )  DEALLOCATE( gridpcbl )
            DEALLOCATE( pcsflt_l )
            DEALLOCATE( kpcsflt_l )
            IF ( debug_output )  THEN
               WRITE( debug_string, '("Finished aggregating ",I0," CSFs.")') ncsfl
               CALL debug_message( debug_string, 'info' )
            ENDIF

        ENDIF

#if defined( __parallel )
        CALL MPI_BARRIER( comm2d, ierr )
#endif
        CALL location_message( 'calculating view factors for radiation interaction', 'finished' )

        RETURN  !todo: remove

!        WRITE( message_string, * )  &
!            'I/O error when processing shape view factors / ',  &
!            'plant canopy sink factors / direct irradiance factors.'
!        CALL message( 'init_urban_surface', 'PA0502', 2, 2, 0, 6, 0 )

    END SUBROUTINE radiation_calc_svf


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Raytracing for detecting obstacles and calculating compound canopy sink
!> factors for RTM. (A simple obstacle detection would only need to process
!> faces in 3 dimensions without any ordering.)
!> Assumtions:
!> -----------
!> 1. The ray always originates from a face midpoint (only one coordinate equals
!>    *.5, i.e. wall) and doesn't travel parallel to the surface (that would mean
!>    shape factor=0). Therefore, the ray may never travel exactly along a face
!>    or an edge.
!> 2. From grid bottom to urban surface top the grid has to be *equidistant*
!>    within each of the dimensions, including vertical (but the resolution
!>    doesn't need to be the same in all three dimensions).
!------------------------------------------------------------------------------!
    SUBROUTINE raytrace(src, targ, isrc, difvf, atarg, create_csf, visible, transparency)
        IMPLICIT NONE

        REAL(wp), DIMENSION(3), INTENT(in)     :: src, targ    !< real coordinates z,y,x
        INTEGER(iwp), INTENT(in)               :: isrc         !< index of source face for csf
        REAL(wp), INTENT(in)                   :: difvf        !< differential view factor for csf
        REAL(wp), INTENT(in)                   :: atarg        !< target surface area for csf
        LOGICAL, INTENT(in)                    :: create_csf   !< whether to generate new CSFs during raytracing
        LOGICAL, INTENT(out)                   :: visible
        REAL(wp), INTENT(out)                  :: transparency !< along whole path
        INTEGER(iwp)                           :: i, k, d
        INTEGER(iwp)                           :: seldim       !< dimension to be incremented
        INTEGER(iwp)                           :: ncsb         !< no of written plant canopy sinkboxes
        INTEGER(iwp)                           :: maxboxes     !< max no of gridboxes visited
        REAL(wp)                               :: distance     !< euclidean along path
        REAL(wp)                               :: crlen        !< length of gridbox crossing
        REAL(wp)                               :: lastdist     !< beginning of current crossing
        REAL(wp)                               :: nextdist     !< end of current crossing
        REAL(wp)                               :: realdist     !< distance in meters per unit distance
        REAL(wp)                               :: crmid        !< midpoint of crossing
        REAL(wp)                               :: cursink      !< sink factor for current canopy box
        REAL(wp), DIMENSION(3)                 :: delta        !< path vector
        REAL(wp), DIMENSION(3)                 :: uvect        !< unit vector
        REAL(wp), DIMENSION(3)                 :: dimnextdist  !< distance for each dimension increments
        INTEGER(iwp), DIMENSION(3)             :: box          !< gridbox being crossed
        INTEGER(iwp), DIMENSION(3)             :: dimnext      !< next dimension increments along path
        INTEGER(iwp), DIMENSION(3)             :: dimdelta     !< dimension direction = +- 1
        INTEGER(iwp)                           :: px, py       !< number of processors in x and y dir before
                                                               !< the processor in the question
        INTEGER(iwp)                           :: ip           !< number of processor where gridbox reside
        INTEGER(iwp)                           :: ig           !< 1D index of gridbox in global 2D array

        REAL(wp)                               :: eps = 1E-10_wp !< epsilon for value comparison
        REAL(wp)                               :: lad_s_target   !< recieved lad_s of particular grid box

!
!--     Maximum number of gridboxes visited equals to maximum number of boundaries crossed in each dimension plus one. That's also
!--     the maximum number of plant canopy boxes written. We grow the acsf array accordingly using exponential factor.
        maxboxes = SUM(ABS(NINT(targ, iwp) - NINT(src, iwp))) + 1
        IF ( plant_canopy  .AND.  ncsfl + maxboxes > ncsfla )  THEN
!--         use this code for growing by fixed exponential increments (equivalent to case where ncsfl always increases by 1)
!--         k = CEILING(grow_factor ** real(CEILING(log(real(ncsfl + maxboxes, kind=wp)) &
!--                                                / log(grow_factor)), kind=wp))
!--         or use this code to simply always keep some extra space after growing
            k = CEILING(REAL(ncsfl + maxboxes, kind=wp) * grow_factor)

            CALL merge_and_grow_csf(k)
        ENDIF

        transparency = 1._wp
        ncsb = 0

        delta(:) = targ(:) - src(:)
        distance = SQRT(SUM(delta(:)**2))
        IF ( distance == 0._wp )  THEN
            visible = .TRUE.
            RETURN
        ENDIF
        uvect(:) = delta(:) / distance
        realdist = SQRT(SUM( (uvect(:)*(/dz(1),dy,dx/))**2 ))

        lastdist = 0._wp

!--     Since all face coordinates have values *.5 and we'd like to use
!--     integers, all these have .5 added
        DO d = 1, 3
            IF ( uvect(d) == 0._wp )  THEN
                dimnext(d) = 999999999
                dimdelta(d) = 999999999
                dimnextdist(d) = 1.0E20_wp
            ELSE IF ( uvect(d) > 0._wp )  THEN
                dimnext(d) = CEILING(src(d) + .5_wp)
                dimdelta(d) = 1
                dimnextdist(d) = (dimnext(d) - .5_wp - src(d)) / uvect(d)
            ELSE
                dimnext(d) = FLOOR(src(d) + .5_wp)
                dimdelta(d) = -1
                dimnextdist(d) = (dimnext(d) - .5_wp - src(d)) / uvect(d)
            ENDIF
        ENDDO

        DO
!--         along what dimension will the next wall crossing be?
            seldim = minloc(dimnextdist, 1)
            nextdist = dimnextdist(seldim)
            IF ( nextdist > distance ) nextdist = distance

            crlen = nextdist - lastdist
            IF ( crlen > .001_wp )  THEN
                crmid = (lastdist + nextdist) * .5_wp
                box = NINT(src(:) + uvect(:) * crmid, iwp)

!--             calculate index of the grid with global indices (box(2),box(3))
!--             in the array nzterr and plantt and id of the coresponding processor
                px = box(3)/nnx
                py = box(2)/nny
                ip = px*pdims(2)+py
                ig = ip*nnx*nny + (box(3)-px*nnx)*nny + box(2)-py*nny
                IF ( box(1) <= nzterr(ig) )  THEN
                    visible = .FALSE.
                    RETURN
                ENDIF

                IF ( plant_canopy )  THEN
                    IF ( box(1) <= plantt(ig) )  THEN
                        ncsb = ncsb + 1
                        boxes(:,ncsb) = box
                        crlens(ncsb) = crlen
#if defined( __parallel )
                        lad_ip(ncsb) = ip
                        lad_disp(ncsb) = (box(3)-px*nnx)*(nny*nz_plant) + (box(2)-py*nny)*nz_plant + box(1)-nz_urban_b
#endif
                    ENDIF
                ENDIF
            ENDIF

            IF ( ABS(distance - nextdist) < eps )  EXIT
            lastdist = nextdist
            dimnext(seldim) = dimnext(seldim) + dimdelta(seldim)
            dimnextdist(seldim) = (dimnext(seldim) - .5_wp - src(seldim)) / uvect(seldim)
        ENDDO

        IF ( plant_canopy )  THEN
#if defined( __parallel )
            IF ( raytrace_mpi_rma )  THEN
!--             send requests for lad_s to appropriate processor
                CALL cpu_log( log_point_s(77), 'rad_rma_lad', 'start' )
                DO i = 1, ncsb
                    CALL MPI_Get(lad_s_ray(i), 1, MPI_REAL, lad_ip(i), lad_disp(i), &
                                 1, MPI_REAL, win_lad, ierr)
                    IF ( ierr /= 0 )  THEN
                        WRITE(9,*) 'Error MPI_Get1:', ierr, lad_s_ray(i), &
                                   lad_ip(i), lad_disp(i), win_lad
                        FLUSH(9)
                    ENDIF
                ENDDO

!--             wait for all pending local requests complete
                CALL MPI_Win_flush_local_all(win_lad, ierr)
                IF ( ierr /= 0 )  THEN
                    WRITE(9,*) 'Error MPI_Win_flush_local_all1:', ierr, win_lad
                    FLUSH(9)
                ENDIF
                CALL cpu_log( log_point_s(77), 'rad_rma_lad', 'stop' )

            ENDIF
#endif

!--         calculate csf and transparency
            DO i = 1, ncsb
#if defined( __parallel )
                IF ( raytrace_mpi_rma )  THEN
                    lad_s_target = lad_s_ray(i)
                ELSE
                    lad_s_target = sub_lad_g(lad_ip(i)*nnx*nny*nz_plant + lad_disp(i))
                ENDIF
#else
                lad_s_target = sub_lad(boxes(1,i),boxes(2,i),boxes(3,i))
#endif
                cursink = 1._wp - exp(-ext_coef * lad_s_target * crlens(i)*realdist)

                IF ( create_csf )  THEN
!--                 write svf values into the array
                    ncsfl = ncsfl + 1
                    acsf(ncsfl)%ip = lad_ip(i)
                    acsf(ncsfl)%itx = boxes(3,i)
                    acsf(ncsfl)%ity = boxes(2,i)
                    acsf(ncsfl)%itz = boxes(1,i)
                    acsf(ncsfl)%isurfs = isrc
                    acsf(ncsfl)%rcvf = cursink*transparency*difvf*atarg
                ENDIF  !< create_csf

                transparency = transparency * (1._wp - cursink)

            ENDDO
        ENDIF

        visible = .TRUE.

    END SUBROUTINE raytrace


!------------------------------------------------------------------------------!
! Description:
! ------------
!> A new, more efficient version of ray tracing algorithm that processes a whole
!> arc instead of a single ray (new in RTM version 2.5).
!>
!> In all comments, horizon means tangent of horizon angle, i.e.
!> vertical_delta / horizontal_distance
!------------------------------------------------------------------------------!
   SUBROUTINE raytrace_2d(origin, yxdir, nrays, zdirs, iorig, aorig, vffrac,  &
                              calc_svf, create_csf, skip_1st_pcb,             &
                              lowest_free_ray, transparency, itarget)
      IMPLICIT NONE

      REAL(wp), DIMENSION(3), INTENT(IN)     ::  origin        !< z,y,x coordinates of ray origin
      REAL(wp), DIMENSION(2), INTENT(IN)     ::  yxdir         !< y,x *unit* vector of ray direction (in grid units)
      INTEGER(iwp)                           ::  nrays         !< number of rays (z directions) to raytrace
      REAL(wp), DIMENSION(nrays), INTENT(IN) ::  zdirs         !< list of z directions to raytrace (z/hdist, grid, zenith->nadir)
      INTEGER(iwp), INTENT(in)               ::  iorig         !< index of origin face for csf
      REAL(wp), INTENT(in)                   ::  aorig         !< origin face area for csf
      REAL(wp), DIMENSION(nrays), INTENT(in) ::  vffrac        !< view factor fractions of each ray for csf
      LOGICAL, INTENT(in)                    ::  calc_svf      !< whether to calculate SFV (identify obstacle surfaces)
      LOGICAL, INTENT(in)                    ::  create_csf    !< whether to create canopy sink factors
      LOGICAL, INTENT(in)                    ::  skip_1st_pcb  !< whether to skip first plant canopy box during raytracing
      INTEGER(iwp), INTENT(out)              ::  lowest_free_ray !< index into zdirs
      REAL(wp), DIMENSION(nrays), INTENT(OUT) ::  transparency !< transparencies of zdirs paths
      INTEGER(iwp), DIMENSION(nrays), INTENT(OUT) ::  itarget  !< global indices of target faces for zdirs

      INTEGER(iwp), DIMENSION(nrays)         ::  target_procs
      REAL(wp)                               ::  horizon       !< highest horizon found after raytracing (z/hdist)
      INTEGER(iwp)                           ::  i, k, l, d
      INTEGER(iwp)                           ::  seldim       !< dimension to be incremented
      REAL(wp), DIMENSION(2)                 ::  yxorigin     !< horizontal copy of origin (y,x)
      REAL(wp)                               ::  distance     !< euclidean along path
      REAL(wp)                               ::  lastdist     !< beginning of current crossing
      REAL(wp)                               ::  nextdist     !< end of current crossing
      REAL(wp)                               ::  crmid        !< midpoint of crossing
      REAL(wp)                               ::  horz_entry   !< horizon at entry to column
      REAL(wp)                               ::  horz_exit    !< horizon at exit from column
      REAL(wp)                               ::  bdydim       !< boundary for current dimension
      REAL(wp), DIMENSION(2)                 ::  crossdist    !< distances to boundary for dimensions
      REAL(wp), DIMENSION(2)                 ::  dimnextdist  !< distance for each dimension increments
      INTEGER(iwp), DIMENSION(2)             ::  column       !< grid column being crossed
      INTEGER(iwp), DIMENSION(2)             ::  dimnext      !< next dimension increments along path
      INTEGER(iwp), DIMENSION(2)             ::  dimdelta     !< dimension direction = +- 1
      INTEGER(iwp)                           ::  px, py       !< number of processors in x and y dir before
                                                              !< the processor in the question
      INTEGER(iwp)                           ::  ip           !< number of processor where gridbox reside
      INTEGER(iwp)                           ::  ig           !< 1D index of gridbox in global 2D array
      INTEGER(iwp)                           ::  maxboxes     !< max no of CSF created
      INTEGER(iwp)                           ::  nly          !< maximum  plant canopy height
      INTEGER(iwp)                           ::  ntrack

      INTEGER(iwp)                           ::  zb0
      INTEGER(iwp)                           ::  zb1
      INTEGER(iwp)                           ::  nz
      INTEGER(iwp)                           ::  iz
      INTEGER(iwp)                           ::  zsgn
      INTEGER(iwp)                           ::  lastdir      !< wall direction before hitting this column
      INTEGER(iwp), DIMENSION(2)             ::  lastcolumn

#if defined( __parallel )
      INTEGER(iwp)                           ::  lowest_lad   !< lowest column cell for which we need LAD
      INTEGER(iwp)                           ::  wcount       !< RMA window item count
      INTEGER(MPI_ADDRESS_KIND)              ::  wdisp        !< RMA window displacement
#endif

      REAL(wp)                               ::  eps = 1E-10_wp !< epsilon for value comparison
      REAL(wp)                               ::  zbottom, ztop !< urban surface boundary in real numbers
      REAL(wp)                               ::  zorig         !< z coordinate of ray column entry
      REAL(wp)                               ::  zexit         !< z coordinate of ray column exit
      REAL(wp)                               ::  qdist         !< ratio of real distance to z coord difference
      REAL(wp)                               ::  dxxyy         !< square of real horizontal distance
      REAL(wp)                               ::  curtrans      !< transparency of current PC box crossing



      yxorigin(:) = origin(2:3)
      transparency(:) = 1._wp !-- Pre-set the all rays to transparent before reducing
      horizon = -HUGE(1._wp)
      lowest_free_ray = nrays
      IF ( rad_angular_discretization  .AND.  calc_svf )  THEN
         ALLOCATE(target_surfl(nrays))
         target_surfl(:) = -1
         lastdir = -999
         lastcolumn(:) = -999
      ENDIF

!--   Determine distance to boundary (in 2D xy)
      IF ( yxdir(1) > 0._wp )  THEN
         bdydim = ny + .5_wp !< north global boundary
         crossdist(1) = (bdydim - yxorigin(1)) / yxdir(1)
      ELSEIF ( yxdir(1) == 0._wp )  THEN
         crossdist(1) = HUGE(1._wp)
      ELSE
          bdydim = -.5_wp !< south global boundary
          crossdist(1) = (bdydim - yxorigin(1)) / yxdir(1)
      ENDIF

      IF ( yxdir(2) > 0._wp )  THEN
          bdydim = nx + .5_wp !< east global boundary
          crossdist(2) = (bdydim - yxorigin(2)) / yxdir(2)
      ELSEIF ( yxdir(2) == 0._wp )  THEN
         crossdist(2) = HUGE(1._wp)
      ELSE
          bdydim = -.5_wp !< west global boundary
          crossdist(2) = (bdydim - yxorigin(2)) / yxdir(2)
      ENDIF
      distance = minval(crossdist, 1)

      IF ( plant_canopy )  THEN
         rt2_track_dist(0) = 0._wp
         rt2_track_lad(:,:) = 0._wp
         nly = plantt_max - nz_urban_b + 1
      ENDIF

      lastdist = 0._wp

!--   Since all face coordinates have values *.5 and we'd like to use
!--   integers, all these have .5 added
      DO  d = 1, 2
          IF ( yxdir(d) == 0._wp )  THEN
              dimnext(d) = HUGE(1_iwp)
              dimdelta(d) = HUGE(1_iwp)
              dimnextdist(d) = HUGE(1._wp)
          ELSE IF ( yxdir(d) > 0._wp )  THEN
              dimnext(d) = FLOOR(yxorigin(d) + .5_wp) + 1
              dimdelta(d) = 1
              dimnextdist(d) = (dimnext(d) - .5_wp - yxorigin(d)) / yxdir(d)
          ELSE
              dimnext(d) = CEILING(yxorigin(d) + .5_wp) - 1
              dimdelta(d) = -1
              dimnextdist(d) = (dimnext(d) - .5_wp - yxorigin(d)) / yxdir(d)
          ENDIF
      ENDDO

      ntrack = 0
      DO
!--      along what dimension will the next wall crossing be?
         seldim = minloc(dimnextdist, 1)
         nextdist = dimnextdist(seldim)
         IF ( nextdist > distance )  nextdist = distance

         IF ( nextdist > lastdist )  THEN
            ntrack = ntrack + 1
            crmid = (lastdist + nextdist) * .5_wp
            column = NINT(yxorigin(:) + yxdir(:) * crmid, iwp)

!--         calculate index of the grid with global indices (column(1),column(2))
!--         in the array nzterr and plantt and id of the coresponding processor
            px = column(2)/nnx
            py = column(1)/nny
            ip = px*pdims(2)+py
            ig = ip*nnx*nny + (column(2)-px*nnx)*nny + column(1)-py*nny

            IF ( lastdist == 0._wp )  THEN
               horz_entry = -HUGE(1._wp)
            ELSE
               horz_entry = (REAL(nzterr(ig), wp) + .5_wp - origin(1)) / lastdist
            ENDIF
            horz_exit = (REAL(nzterr(ig), wp) + .5_wp - origin(1)) / nextdist

            IF ( rad_angular_discretization  .AND.  calc_svf )  THEN
!
!--            Identify vertical obstacles hit by rays in current column
               DO WHILE ( lowest_free_ray > 0 )
                  IF ( zdirs(lowest_free_ray) > horz_entry )  EXIT
!
!--               This may only happen after 1st column, so lastdir and lastcolumn are valid
                  CALL request_itarget(lastdir,                                         &
                        CEILING(-0.5_wp + origin(1) + zdirs(lowest_free_ray)*lastdist), &
                        lastcolumn(1), lastcolumn(2),                                   &
                        target_surfl(lowest_free_ray), target_procs(lowest_free_ray))
                  lowest_free_ray = lowest_free_ray - 1
               ENDDO
!
!--            Identify horizontal obstacles hit by rays in current column
               DO WHILE ( lowest_free_ray > 0 )
                  IF ( zdirs(lowest_free_ray) > horz_exit )  EXIT
                  CALL request_itarget(iup_u, nzterr(ig)+1, column(1), column(2), &
                                       target_surfl(lowest_free_ray),           &
                                       target_procs(lowest_free_ray))
                  lowest_free_ray = lowest_free_ray - 1
               ENDDO
            ENDIF

            horizon = MAX(horizon, horz_entry, horz_exit)

            IF ( plant_canopy )  THEN
               rt2_track(:, ntrack) = column(:)
               rt2_track_dist(ntrack) = nextdist
            ENDIF
         ENDIF

         IF ( nextdist + eps >= distance )  EXIT

         IF ( rad_angular_discretization  .AND.  calc_svf )  THEN
!
!--         Save wall direction of coming building column (= this air column)
            IF ( seldim == 1 )  THEN
               IF ( dimdelta(seldim) == 1 )  THEN
                  lastdir = isouth_u
               ELSE
                  lastdir = inorth_u
               ENDIF
            ELSE
               IF ( dimdelta(seldim) == 1 )  THEN
                  lastdir = iwest_u
               ELSE
                  lastdir = ieast_u
               ENDIF
            ENDIF
            lastcolumn = column
         ENDIF
         lastdist = nextdist
         dimnext(seldim) = dimnext(seldim) + dimdelta(seldim)
         dimnextdist(seldim) = (dimnext(seldim) - .5_wp - yxorigin(seldim)) / yxdir(seldim)
      ENDDO

      IF ( plant_canopy )  THEN
!--      Request LAD WHERE applicable
!--
#if defined( __parallel )
         IF ( raytrace_mpi_rma )  THEN
!--         send requests for lad_s to appropriate processor
            !CALL cpu_log( log_point_s(77), 'usm_init_rma', 'start' )
            DO  i = 1, ntrack
               px = rt2_track(2,i)/nnx
               py = rt2_track(1,i)/nny
               ip = px*pdims(2)+py
               ig = ip*nnx*nny + (rt2_track(2,i)-px*nnx)*nny + rt2_track(1,i)-py*nny

               IF ( rad_angular_discretization  .AND.  calc_svf )  THEN
!
!--               For fixed view resolution, we need plant canopy even for rays
!--               to opposing surfaces
                  lowest_lad = nzterr(ig) + 1
               ELSE
!
!--               We only need LAD for rays directed above horizon (to sky)
                  lowest_lad = CEILING( -0.5_wp + origin(1) +            &
                                    MIN( horizon * rt2_track_dist(i-1),  & ! entry
                                         horizon * rt2_track_dist(i)   ) ) ! exit
               ENDIF
!
!--            Skip asking for LAD where all plant canopy is under requested level
               IF ( plantt(ig) < lowest_lad )  CYCLE

               wdisp = (rt2_track(2,i)-px*nnx)*(nny*nz_plant) + (rt2_track(1,i)-py*nny)*nz_plant + lowest_lad-nz_urban_b
               wcount = plantt(ig)-lowest_lad+1
               ! TODO send request ASAP - even during raytracing
               CALL MPI_Get(rt2_track_lad(lowest_lad:plantt(ig), i), wcount, MPI_REAL, ip,    &
                            wdisp, wcount, MPI_REAL, win_lad, ierr)
               IF ( ierr /= 0 )  THEN
                  WRITE(9,*) 'Error MPI_Get2:', ierr, rt2_track_lad(lowest_lad:plantt(ig), i), &
                             wcount, ip, wdisp, win_lad
                  FLUSH(9)
               ENDIF
            ENDDO

!--         wait for all pending local requests complete
            ! TODO WAIT selectively for each column later when needed
            CALL MPI_Win_flush_local_all(win_lad, ierr)
            IF ( ierr /= 0 )  THEN
               WRITE(9,*) 'Error MPI_Win_flush_local_all2:', ierr, win_lad
               FLUSH(9)
            ENDIF
            !CALL cpu_log( log_point_s(77), 'usm_init_rma', 'stop' )

         ELSE ! raytrace_mpi_rma = .F.
            DO  i = 1, ntrack
               px = rt2_track(2,i)/nnx
               py = rt2_track(1,i)/nny
               ip = px*pdims(2)+py
               ig = ip*nnx*nny*nz_plant + (rt2_track(2,i)-px*nnx)*(nny*nz_plant) + (rt2_track(1,i)-py*nny)*nz_plant
               rt2_track_lad(nz_urban_b:plantt_max, i) = sub_lad_g(ig:ig+nly-1)
            ENDDO
         ENDIF
#else
         DO  i = 1, ntrack
            rt2_track_lad(nz_urban_b:plantt_max, i) = sub_lad(rt2_track(1,i), rt2_track(2,i), nz_urban_b:plantt_max)
         ENDDO
#endif
      ENDIF ! plant_canopy

      IF ( rad_angular_discretization  .AND.  calc_svf )  THEN
#if defined( __parallel )
!--      wait for all gridsurf requests to complete
         CALL MPI_Win_flush_local_all(win_gridsurf, ierr)
         IF ( ierr /= 0 )  THEN
            WRITE(9,*) 'Error MPI_Win_flush_local_all3:', ierr, win_gridsurf
            FLUSH(9)
         ENDIF
#endif
!
!--      recalculate local surf indices into global ones
         DO i = 1, nrays
            IF ( target_surfl(i) == -1 )  THEN
               itarget(i) = -1
            ELSE
               itarget(i) = target_surfl(i) + surfstart(target_procs(i))
            ENDIF
         ENDDO

         DEALLOCATE( target_surfl )

      ELSE
         itarget(:) = -1
      ENDIF ! rad_angular_discretization

      IF ( plant_canopy )  THEN
!--      Skip the PCB around origin if requested (for MRT, the PCB might not be there)
!--
         IF ( skip_1st_pcb  .AND.  NINT(origin(1)) <= plantt_max )  THEN
            rt2_track_lad(NINT(origin(1), iwp), 1) = 0._wp
         ENDIF

!--      Assert that we have space allocated for CSFs
!--
         maxboxes = (ntrack + MAX(CEILING(origin(1)-.5_wp) - nz_urban_b,          &
                                  nz_urban_t - CEILING(origin(1)-.5_wp))) * nrays
         IF ( ncsfl + maxboxes > ncsfla )  THEN
!--         use this code for growing by fixed exponential increments (equivalent to case where ncsfl always increases by 1)
!--         k = CEILING(grow_factor ** real(CEILING(log(real(ncsfl + maxboxes, kind=wp)) &
!--                                                / log(grow_factor)), kind=wp))
!--         or use this code to simply always keep some extra space after growing
            k = CEILING(REAL(ncsfl + maxboxes, kind=wp) * grow_factor)
            CALL merge_and_grow_csf(k)
         ENDIF

!--      Calculate transparencies and store new CSFs
!--
         zbottom = REAL(nz_urban_b, wp) - .5_wp
         ztop = REAL(plantt_max, wp) + .5_wp

!--      Reverse direction of radiation (face->sky), only when calc_svf
!--
         IF ( calc_svf )  THEN
            DO  i = 1, ntrack ! for each column
               dxxyy = ((dy*yxdir(1))**2 + (dx*yxdir(2))**2) * (rt2_track_dist(i)-rt2_track_dist(i-1))**2
               px = rt2_track(2,i)/nnx
               py = rt2_track(1,i)/nny
               ip = px*pdims(2)+py

               DO  k = 1, nrays ! for each ray
!
!--               NOTE 6778:
!--               With traditional svf discretization, CSFs under the horizon
!--               (i.e. for surface to surface radiation)  are created in
!--               raytrace(). With rad_angular_discretization, we must create
!--               CSFs under horizon only for one direction, otherwise we would
!--               have duplicate amount of energy. Although we could choose
!--               either of the two directions (they differ only by
!--               discretization error with no bias), we choose the the backward
!--               direction, because it tends to cumulate high canopy sink
!--               factors closer to raytrace origin, i.e. it should potentially
!--               cause less moiree.
                  IF ( .NOT. rad_angular_discretization )  THEN
                     IF ( zdirs(k) <= horizon )  CYCLE
                  ENDIF

                  zorig = origin(1) + zdirs(k) * rt2_track_dist(i-1)
                  IF ( zorig <= zbottom  .OR.  zorig >= ztop )  CYCLE

                  zsgn = INT(SIGN(1._wp, zdirs(k)), iwp)
                  rt2_dist(1) = 0._wp
                  IF ( zdirs(k) == 0._wp )  THEN ! ray is exactly horizontal
                     nz = 2
                     rt2_dist(nz) = SQRT(dxxyy)
                     iz = CEILING(-.5_wp + zorig, iwp)
                  ELSE
                     zexit = MIN(MAX(origin(1) + zdirs(k) * rt2_track_dist(i), zbottom), ztop)

                     zb0 = FLOOR(  zorig * zsgn - .5_wp) + 1  ! because it must be greater than orig
                     zb1 = CEILING(zexit * zsgn - .5_wp) - 1  ! because it must be smaller than exit
                     nz = MAX(zb1 - zb0 + 3, 2)
                     rt2_dist(nz) = SQRT(((zexit-zorig)*dz(1))**2 + dxxyy)
                     qdist = rt2_dist(nz) / (zexit-zorig)
                     rt2_dist(2:nz-1) = (/( ((REAL(l, wp) + .5_wp) * zsgn - zorig) * qdist , l = zb0, zb1 )/)
                     iz = zb0 * zsgn
                  ENDIF

                  DO  l = 2, nz
                     IF ( rt2_track_lad(iz, i) > 0._wp )  THEN
                        curtrans = exp(-ext_coef * rt2_track_lad(iz, i) * (rt2_dist(l)-rt2_dist(l-1)))

                        IF ( create_csf )  THEN
                           ncsfl = ncsfl + 1
                           acsf(ncsfl)%ip = ip
                           acsf(ncsfl)%itx = rt2_track(2,i)
                           acsf(ncsfl)%ity = rt2_track(1,i)
                           acsf(ncsfl)%itz = iz
                           acsf(ncsfl)%isurfs = iorig
                           acsf(ncsfl)%rcvf = (1._wp - curtrans)*transparency(k)*vffrac(k)
                        ENDIF

                        transparency(k) = transparency(k) * curtrans
                     ENDIF
                     iz = iz + zsgn
                  ENDDO ! l = 1, nz - 1
               ENDDO ! k = 1, nrays
            ENDDO ! i = 1, ntrack

            transparency(1:lowest_free_ray) = 1._wp !-- Reset rays above horizon to transparent (see NOTE 6778)
         ENDIF

!--      Forward direction of radiation (sky->face), always
!--
         DO  i = ntrack, 1, -1 ! for each column backwards
            dxxyy = ((dy*yxdir(1))**2 + (dx*yxdir(2))**2) * (rt2_track_dist(i)-rt2_track_dist(i-1))**2
            px = rt2_track(2,i)/nnx
            py = rt2_track(1,i)/nny
            ip = px*pdims(2)+py

            DO  k = 1, nrays ! for each ray
!
!--            See NOTE 6778 above
               IF ( zdirs(k) <= horizon )  CYCLE

               zexit = origin(1) + zdirs(k) * rt2_track_dist(i-1)
               IF ( zexit <= zbottom  .OR.  zexit >= ztop )  CYCLE

               zsgn = -INT(SIGN(1._wp, zdirs(k)), iwp)
               rt2_dist(1) = 0._wp
               IF ( zdirs(k) == 0._wp )  THEN ! ray is exactly horizontal
                  nz = 2
                  rt2_dist(nz) = SQRT(dxxyy)
                  iz = NINT(zexit, iwp)
               ELSE
                  zorig = MIN(MAX(origin(1) + zdirs(k) * rt2_track_dist(i), zbottom), ztop)

                  zb0 = FLOOR(  zorig * zsgn - .5_wp) + 1  ! because it must be greater than orig
                  zb1 = CEILING(zexit * zsgn - .5_wp) - 1  ! because it must be smaller than exit
                  nz = MAX(zb1 - zb0 + 3, 2)
                  rt2_dist(nz) = SQRT(((zexit-zorig)*dz(1))**2 + dxxyy)
                  qdist = rt2_dist(nz) / (zexit-zorig)
                  rt2_dist(2:nz-1) = (/( ((REAL(l, wp) + .5_wp) * zsgn - zorig) * qdist , l = zb0, zb1 )/)
                  iz = zb0 * zsgn
               ENDIF

               DO  l = 2, nz
                  IF ( rt2_track_lad(iz, i) > 0._wp )  THEN
                     curtrans = exp(-ext_coef * rt2_track_lad(iz, i) * (rt2_dist(l)-rt2_dist(l-1)))

                     IF ( create_csf )  THEN
                        ncsfl = ncsfl + 1
                        acsf(ncsfl)%ip = ip
                        acsf(ncsfl)%itx = rt2_track(2,i)
                        acsf(ncsfl)%ity = rt2_track(1,i)
                        acsf(ncsfl)%itz = iz
                        IF ( itarget(k) /= -1 ) STOP 1 !FIXME remove after test
                        acsf(ncsfl)%isurfs = -1
                        acsf(ncsfl)%rcvf = (1._wp - curtrans)*transparency(k)*aorig*vffrac(k)
                     ENDIF  ! create_csf

                     transparency(k) = transparency(k) * curtrans
                  ENDIF
                  iz = iz + zsgn
               ENDDO ! l = 1, nz - 1
            ENDDO ! k = 1, nrays
         ENDDO ! i = 1, ntrack
      ENDIF ! plant_canopy

      IF ( .NOT. (rad_angular_discretization  .AND.  calc_svf) )  THEN
!
!--      Just update lowest_free_ray according to horizon
         DO WHILE ( lowest_free_ray > 0 )
            IF ( zdirs(lowest_free_ray) > horizon )  EXIT
            lowest_free_ray = lowest_free_ray - 1
         ENDDO
      ENDIF

   CONTAINS

      SUBROUTINE request_itarget( d, z, y, x, isurfl, iproc )

         INTEGER(iwp), INTENT(in)            ::  d, z, y, x
         INTEGER(iwp), TARGET, INTENT(out)   ::  isurfl
         INTEGER(iwp), INTENT(out)           ::  iproc
#if defined( __parallel )
         INTEGER(iwp)                        ::  px, py        !< number of processors in x and y direction
                                                               !< before the processor in the question
#endif

#if defined( __parallel )
         INTEGER(KIND=MPI_ADDRESS_KIND)      ::  target_displ  !< index of the grid in the local gridsurf array

!
!--      Calculate target processor and index in the remote local target gridsurf array
         px = x / nnx
         py = y / nny
         iproc = px * pdims(2) + py
         target_displ = ((x-px*nnx) * nny + y - py*nny ) * nz_urban * nsurf_type_u +&
                        ( z-nz_urban_b ) * nsurf_type_u + d
!
!--      Send MPI_Get request to obtain index target_surfl(i)
         CALL MPI_GET( isurfl, 1, MPI_INTEGER, iproc, target_displ,            &
                       1, MPI_INTEGER, win_gridsurf, ierr)
         IF ( ierr /= 0 )  THEN
            WRITE( 9,* ) 'Error MPI_Get3:', ierr, isurfl, iproc, target_displ, &
                         win_gridsurf
            FLUSH( 9 )
         ENDIF
#else
!--      set index target_surfl(i)
         isurfl = gridsurf(d,z,y,x)
         iproc  = 0  ! required to avoid compile error about unused variable in serial mode
#endif

      END SUBROUTINE request_itarget

   END SUBROUTINE raytrace_2d


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Calculates apparent solar positions for all timesteps and stores discretized
!> positions for RTM.
!------------------------------------------------------------------------------!
   SUBROUTINE radiation_presimulate_solar_pos

      IMPLICIT NONE

      INTEGER(iwp) ::  it, i, j                           !< loop indices

      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: dsidir_tmp !< dsidir_tmp[:,i] = unit vector of i-th
                                                          !< appreant solar direction

      ALLOCATE ( dsidir_rev(0:raytrace_discrete_elevs/2-1,                 &
                            0:raytrace_discrete_azims-1) )
      dsidir_rev(:,:) = -1
      ALLOCATE ( dsidir_tmp(3,                                             &
                     raytrace_discrete_elevs/2*raytrace_discrete_azims) )
      ndsidir = 0
      sun_direction = .TRUE.

!
!--   Process spinup time if configured
      IF ( spinup_time > 0._wp )  THEN
         DO  it = 0, CEILING(spinup_time / dt_spinup)
            CALL simulate_pos( it * dt_spinup - spinup_time )
         ENDDO
      ENDIF
!
!--   Process simulation time
      DO  it = 0, CEILING(( end_time - spinup_time ) / dt_radiation)
         CALL simulate_pos( it * dt_radiation )
      ENDDO
!
!--   Allocate global vars which depend on ndsidir
      ALLOCATE ( dsidir ( 3, ndsidir ) )
      dsidir(:,:) = dsidir_tmp(:, 1:ndsidir)
      DEALLOCATE ( dsidir_tmp )

      ALLOCATE ( dsitrans(nsurfl, ndsidir) )
      ALLOCATE ( dsitransc(npcbl, ndsidir) )
      IF ( nmrtbl > 0 )  ALLOCATE ( mrtdsit(nmrtbl, ndsidir) )

      WRITE ( message_string, * ) 'Precalculated', ndsidir, ' solar positions', &
                                  ' from', it, ' timesteps.'
      CALL message( 'radiation_presimulate_solar_pos', 'UI0013', 0, 0, 0, 6, 0 )

      CONTAINS

      !------------------------------------------------------------------------!
      ! Description:
      ! ------------
      !> Simuates a single position
      !------------------------------------------------------------------------!
      SUBROUTINE simulate_pos( time_since_reference_local )

         REAL(wp), INTENT(IN) ::  time_since_reference_local  !< local time since reference
!
!--      Update apparent solar position based on modified t_s_r_p
         CALL get_date_time( time_since_reference_local, &
                             day_of_year=day_of_year,    &
                             second_of_day=second_of_day )
         CALL calc_zenith( day_of_year, second_of_day )
         IF ( cos_zenith > 0 )  THEN
!--
!--         Identify solar direction vector (discretized number) 1)
            i = MODULO(NINT(ATAN2(sun_dir_lon, sun_dir_lat)               &
                            / (2._wp*pi) * raytrace_discrete_azims-.5_wp, iwp), &
                       raytrace_discrete_azims)
            j = FLOOR(ACOS(cos_zenith) / pi * raytrace_discrete_elevs)
            IF ( dsidir_rev(j, i) == -1 )  THEN
               ndsidir = ndsidir + 1
               dsidir_tmp(:, ndsidir) =                                              &
                     (/ COS((REAL(j,wp)+.5_wp) * pi      / raytrace_discrete_elevs), &
                        SIN((REAL(j,wp)+.5_wp) * pi      / raytrace_discrete_elevs)  &
                      * COS((REAL(i,wp)+.5_wp) * 2_wp*pi / raytrace_discrete_azims), &
                        SIN((REAL(j,wp)+.5_wp) * pi      / raytrace_discrete_elevs)  &
                      * SIN((REAL(i,wp)+.5_wp) * 2_wp*pi / raytrace_discrete_azims) /)
               dsidir_rev(j, i) = ndsidir
            ENDIF
         ENDIF
      END SUBROUTINE simulate_pos

   END SUBROUTINE radiation_presimulate_solar_pos



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Determines whether two faces are oriented towards each other in RTM. Since the
!> surfaces follow the gird box surfaces, it checks first whether the two surfaces
!> are directed in the same direction, then it checks if the two surfaces are
!> located in confronted direction but facing away from each other, e.g. <--| |-->
!------------------------------------------------------------------------------!
    PURE LOGICAL FUNCTION surface_facing(x, y, z, d, x2, y2, z2, d2)
        IMPLICIT NONE
        INTEGER(iwp),   INTENT(in)  :: x, y, z, d, x2, y2, z2, d2

        surface_facing = .FALSE.

!-- first check: are the two surfaces directed in the same direction
        IF ( (d==iup_u  .OR.  d==iup_l )                             &
             .AND. (d2==iup_u  .OR. d2==iup_l) ) RETURN
        IF ( (d==isouth_u  .OR.  d==isouth_l ) &
             .AND.  (d2==isouth_u  .OR.  d2==isouth_l) ) RETURN
        IF ( (d==inorth_u  .OR.  d==inorth_l ) &
             .AND.  (d2==inorth_u  .OR.  d2==inorth_l) ) RETURN
        IF ( (d==iwest_u  .OR.  d==iwest_l )     &
             .AND.  (d2==iwest_u  .OR.  d2==iwest_l ) ) RETURN
        IF ( (d==ieast_u  .OR.  d==ieast_l )     &
             .AND.  (d2==ieast_u  .OR.  d2==ieast_l ) ) RETURN

!-- second check: are surfaces facing away from each other
        SELECT CASE (d)
            CASE (iup_u, iup_l)                     !< upward facing surfaces
                IF ( z2 < z ) RETURN
            CASE (isouth_u, isouth_l)               !< southward facing surfaces
                IF ( y2 > y ) RETURN
            CASE (inorth_u, inorth_l)               !< northward facing surfaces
                IF ( y2 < y ) RETURN
            CASE (iwest_u, iwest_l)                 !< westward facing surfaces
                IF ( x2 > x ) RETURN
            CASE (ieast_u, ieast_l)                 !< eastward facing surfaces
                IF ( x2 < x ) RETURN
        END SELECT

        SELECT CASE (d2)
            CASE (iup_u)                            !< ground, roof
                IF ( z < z2 ) RETURN
            CASE (isouth_u, isouth_l)               !< south facing
                IF ( y > y2 ) RETURN
            CASE (inorth_u, inorth_l)               !< north facing
                IF ( y < y2 ) RETURN
            CASE (iwest_u, iwest_l)                 !< west facing
                IF ( x > x2 ) RETURN
            CASE (ieast_u, ieast_l)                 !< east facing
                IF ( x < x2 ) RETURN
            CASE (-1)
                CONTINUE
        END SELECT

        surface_facing = .TRUE.

    END FUNCTION surface_facing


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Reads svf, svfsurf, csf, csfsurf and mrt factors data from saved file.
!> This allows to skip their calculation during of RTM init phase.
!> SVF means sky view factors and CSF means canopy sink factors.
!------------------------------------------------------------------------------!
    SUBROUTINE radiation_read_svf

       IMPLICIT NONE

       CHARACTER(rad_version_len)   :: rad_version_field

       INTEGER(iwp)                 :: i
       INTEGER(iwp)                 :: ndsidir_from_file = 0
       INTEGER(iwp)                 :: npcbl_from_file = 0
       INTEGER(iwp)                 :: nsurfl_from_file = 0
       INTEGER(iwp)                 :: nmrtbl_from_file = 0


       CALL location_message( 'reading view factors for radiation interaction', 'start' )

       DO  i = 0, io_blocks-1
          IF ( i == io_group )  THEN

!
!--          numprocs_previous_run is only known in case of reading restart
!--          data. If a new initial run which reads svf data is started the
!--          following query will be skipped
             IF ( initializing_actions == 'read_restart_data' ) THEN

                IF ( numprocs_previous_run /= numprocs ) THEN
                   WRITE( message_string, * ) 'A different number of ',        &
                                              'processors between the run ',   &
                                              'that has written the svf data ',&
                                              'and the one that will read it ',&
                                              'is not allowed'
                   CALL message( 'check_open', 'PA0491', 1, 2, 0, 6, 0 )
                ENDIF

             ENDIF

!
!--          Open binary file
             CALL check_open( 88 )

!
!--          read and check version
             READ ( 88 ) rad_version_field
             IF ( TRIM(rad_version_field) /= TRIM(rad_version) )  THEN
                 WRITE( message_string, * ) 'Version of binary SVF file "',    &
                             TRIM(rad_version_field), '" does not match ',     &
                             'the version of model "', TRIM(rad_version), '"'
                 CALL message( 'radiation_read_svf', 'PA0482', 1, 2, 0, 6, 0 )
             ENDIF

!
!--          read nsvfl, ncsfl, nsurfl, nmrtf
             READ ( 88 ) nsvfl, ncsfl, nsurfl_from_file, npcbl_from_file,      &
                         ndsidir_from_file, nmrtbl_from_file, nmrtf

             IF ( nsvfl < 0  .OR.  ncsfl < 0 )  THEN
                 WRITE( message_string, * ) 'Wrong number of SVF or CSF'
                 CALL message( 'radiation_read_svf', 'PA0483', 1, 2, 0, 6, 0 )
             ELSE
                 WRITE(debug_string,*)   'Number of SVF, CSF, and nsurfl ',    &
                                         'to read', nsvfl, ncsfl,              &
                                         nsurfl_from_file
                 IF ( debug_output )  CALL debug_message( debug_string, 'info' )
             ENDIF

             IF ( nsurfl_from_file /= nsurfl )  THEN
                 WRITE( message_string, * ) 'nsurfl from SVF file does not ',  &
                                            'match calculated nsurfl from ',   &
                                            'radiation_interaction_init'
                 CALL message( 'radiation_read_svf', 'PA0490', 1, 2, 0, 6, 0 )
             ENDIF

             IF ( npcbl_from_file /= npcbl )  THEN
                 WRITE( message_string, * ) 'npcbl from SVF file does not ',   &
                                            'match calculated npcbl from ',    &
                                            'radiation_interaction_init'
                 CALL message( 'radiation_read_svf', 'PA0493', 1, 2, 0, 6, 0 )
             ENDIF

             IF ( ndsidir_from_file /= ndsidir )  THEN
                 WRITE( message_string, * ) 'ndsidir from SVF file does not ', &
                                            'match calculated ndsidir from ',  &
                                            'radiation_presimulate_solar_pos'
                 CALL message( 'radiation_read_svf', 'PA0494', 1, 2, 0, 6, 0 )
             ENDIF
             IF ( nmrtbl_from_file /= nmrtbl )  THEN
                 WRITE( message_string, * ) 'nmrtbl from SVF file does not ',  &
                                            'match calculated nmrtbl from ',   &
                                            'radiation_interaction_init'
                 CALL message( 'radiation_read_svf', 'PA0494', 1, 2, 0, 6, 0 )
             ELSE
                 WRITE(debug_string,*) 'Number of nmrtf to read ', nmrtf
                 IF ( debug_output )  CALL debug_message( debug_string, 'info' )
             ENDIF

!
!--          Arrays skyvf, skyvft, dsitrans and dsitransc are allready
!--          allocated in radiation_interaction_init and
!--          radiation_presimulate_solar_pos
             IF ( nsurfl > 0 )  THEN
                READ(88) skyvf
                READ(88) skyvft
                READ(88) dsitrans
             ENDIF

             IF ( plant_canopy  .AND.  npcbl > 0 ) THEN
                READ ( 88 )  dsitransc
             ENDIF

!
!--          The allocation of svf, svfsurf, csf, csfsurf, mrtf, mrtft, and
!--          mrtfsurf happens in routine radiation_calc_svf which is not
!--          called if the program enters radiation_read_svf. Therefore
!--          these arrays has to allocate in the following
             IF ( nsvfl > 0 )  THEN
                ALLOCATE( svf(ndsvf,nsvfl) )
                ALLOCATE( svfsurf(idsvf,nsvfl) )
                READ(88) svf
                READ(88) svfsurf
             ENDIF

             IF ( plant_canopy  .AND.  ncsfl > 0 )  THEN
                ALLOCATE( csf(ndcsf,ncsfl) )
                ALLOCATE( csfsurf(idcsf,ncsfl) )
                READ(88) csf
                READ(88) csfsurf
             ENDIF

             IF ( nmrtbl > 0 )  THEN
                READ(88) mrtsky
                READ(88) mrtskyt
                READ(88) mrtdsit
             ENDIF

             IF ( nmrtf > 0 )  THEN
                ALLOCATE ( mrtf(nmrtf) )
                ALLOCATE ( mrtft(nmrtf) )
                ALLOCATE ( mrtfsurf(2,nmrtf) )
                READ(88) mrtf
                READ(88) mrtft
                READ(88) mrtfsurf
             ENDIF

!
!--          Close binary file
             CALL close_file( 88 )

          ENDIF
#if defined( __parallel )
          CALL MPI_BARRIER( comm2d, ierr )
#endif
       ENDDO

       CALL location_message( 'reading view factors for radiation interaction', 'finished' )


    END SUBROUTINE radiation_read_svf


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine stores svf, svfsurf, csf, csfsurf and mrt data to a file.
!> The stored factors can be reused in future simulation with the same
!> geometry structure of the surfaces and resolved plant canopy.
!------------------------------------------------------------------------------!
    SUBROUTINE radiation_write_svf

       IMPLICIT NONE

       INTEGER(iwp)        :: i


       CALL location_message( 'writing view factors for radiation interaction', 'start' )

       DO  i = 0, io_blocks-1
          IF ( i == io_group )  THEN
!
!--          Open binary file
             CALL check_open( 89 )

             WRITE ( 89 )  rad_version
             WRITE ( 89 )  nsvfl, ncsfl, nsurfl, npcbl, ndsidir, nmrtbl, nmrtf
             IF ( nsurfl > 0 ) THEN
                WRITE ( 89 )  skyvf
                WRITE ( 89 )  skyvft
                WRITE ( 89 )  dsitrans
             ENDIF
             IF ( npcbl > 0 ) THEN
                WRITE ( 89 )  dsitransc
             ENDIF
             IF ( nsvfl > 0 ) THEN
                WRITE ( 89 )  svf
                WRITE ( 89 )  svfsurf
             ENDIF
             IF ( plant_canopy  .AND.  ncsfl > 0 )  THEN
                 WRITE ( 89 )  csf
                 WRITE ( 89 )  csfsurf
             ENDIF
             IF ( nmrtbl > 0 )  THEN
                WRITE ( 89 ) mrtsky
                WRITE ( 89 ) mrtskyt
                WRITE ( 89 ) mrtdsit
             ENDIF
             IF ( nmrtf > 0 )  THEN
                 WRITE ( 89 )  mrtf
                 WRITE ( 89 )  mrtft
                 WRITE ( 89 )  mrtfsurf
             ENDIF
!
!--          Close binary file
             CALL close_file( 89 )

          ENDIF
#if defined( __parallel )
          CALL MPI_BARRIER( comm2d, ierr )
#endif
       ENDDO

       CALL location_message( 'writing view factors for radiation interaction', 'finished' )


    END SUBROUTINE radiation_write_svf


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Block of auxiliary subroutines for RTM:
!> 1. quicksort and corresponding comparison
!> 2. merge_and_grow_csf for implementation of "dynamical growing"
!>    array for csf
!------------------------------------------------------------------------------!
!-- quicksort.f -*-f90-*-
!-- Author: t-nissie, adaptation J.Resler
!-- License: GPLv3
!-- Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
    RECURSIVE SUBROUTINE quicksort_itarget(itarget, vffrac, ztransp, first, last)
        IMPLICIT NONE
        INTEGER(iwp), DIMENSION(:), INTENT(INOUT)   :: itarget
        REAL(wp), DIMENSION(:), INTENT(INOUT)       :: vffrac, ztransp
        INTEGER(iwp), INTENT(IN)                    :: first, last
        INTEGER(iwp)                                :: x, t
        INTEGER(iwp)                                :: i, j
        REAL(wp)                                    :: tr

        IF ( first>=last ) RETURN
        x = itarget((first+last)/2)
        i = first
        j = last
        DO
            DO WHILE ( itarget(i) < x )
               i=i+1
            ENDDO
            DO WHILE ( x < itarget(j) )
                j=j-1
            ENDDO
            IF ( i >= j ) EXIT
            t = itarget(i);  itarget(i) = itarget(j);  itarget(j) = t
            tr = vffrac(i);  vffrac(i) = vffrac(j);  vffrac(j) = tr
            tr = ztransp(i);  ztransp(i) = ztransp(j);  ztransp(j) = tr
            i=i+1
            j=j-1
        ENDDO
        IF ( first < i-1 ) CALL quicksort_itarget(itarget, vffrac, ztransp, first, i-1)
        IF ( j+1 < last )  CALL quicksort_itarget(itarget, vffrac, ztransp, j+1, last)
    END SUBROUTINE quicksort_itarget

    PURE FUNCTION svf_lt(svf1,svf2) result (res)
      TYPE (t_svf), INTENT(in) :: svf1,svf2
      LOGICAL                  :: res
      IF ( svf1%isurflt < svf2%isurflt  .OR.    &
          (svf1%isurflt == svf2%isurflt  .AND.  svf1%isurfs < svf2%isurfs) )  THEN
          res = .TRUE.
      ELSE
          res = .FALSE.
      ENDIF
    END FUNCTION svf_lt


!-- quicksort.f -*-f90-*-
!-- Author: t-nissie, adaptation J.Resler
!-- License: GPLv3
!-- Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
    RECURSIVE SUBROUTINE quicksort_svf(svfl, first, last)
        IMPLICIT NONE
        TYPE(t_svf), DIMENSION(:), INTENT(INOUT)  :: svfl
        INTEGER(iwp), INTENT(IN)                  :: first, last
        TYPE(t_svf)                               :: x, t
        INTEGER(iwp)                              :: i, j

        IF ( first>=last ) RETURN
        x = svfl( (first+last) / 2 )
        i = first
        j = last
        DO
            DO while ( svf_lt(svfl(i),x) )
               i=i+1
            ENDDO
            DO while ( svf_lt(x,svfl(j)) )
                j=j-1
            ENDDO
            IF ( i >= j ) EXIT
            t = svfl(i);  svfl(i) = svfl(j);  svfl(j) = t
            i=i+1
            j=j-1
        ENDDO
        IF ( first < i-1 ) CALL quicksort_svf(svfl, first, i-1)
        IF ( j+1 < last )  CALL quicksort_svf(svfl, j+1, last)
    END SUBROUTINE quicksort_svf

    PURE FUNCTION csf_lt(csf1,csf2) result (res)
      TYPE (t_csf), INTENT(in) :: csf1,csf2
      LOGICAL                  :: res
      IF ( csf1%ip < csf2%ip  .OR.    &
           (csf1%ip == csf2%ip  .AND.  csf1%itx < csf2%itx)  .OR.  &
           (csf1%ip == csf2%ip  .AND.  csf1%itx == csf2%itx  .AND.  csf1%ity < csf2%ity)  .OR.  &
           (csf1%ip == csf2%ip  .AND.  csf1%itx == csf2%itx  .AND.  csf1%ity == csf2%ity  .AND.   &
            csf1%itz < csf2%itz)  .OR.  &
           (csf1%ip == csf2%ip  .AND.  csf1%itx == csf2%itx  .AND.  csf1%ity == csf2%ity  .AND.   &
            csf1%itz == csf2%itz  .AND.  csf1%isurfs < csf2%isurfs) )  THEN
          res = .TRUE.
      ELSE
          res = .FALSE.
      ENDIF
    END FUNCTION csf_lt


!-- quicksort.f -*-f90-*-
!-- Author: t-nissie, adaptation J.Resler
!-- License: GPLv3
!-- Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
    RECURSIVE SUBROUTINE quicksort_csf(csfl, first, last)
        IMPLICIT NONE
        TYPE(t_csf), DIMENSION(:), INTENT(INOUT)  :: csfl
        INTEGER(iwp), INTENT(IN)                  :: first, last
        TYPE(t_csf)                               :: x, t
        INTEGER(iwp)                              :: i, j

        IF ( first>=last ) RETURN
        x = csfl( (first+last)/2 )
        i = first
        j = last
        DO
            DO while ( csf_lt(csfl(i),x) )
                i=i+1
            ENDDO
            DO while ( csf_lt(x,csfl(j)) )
                j=j-1
            ENDDO
            IF ( i >= j ) EXIT
            t = csfl(i);  csfl(i) = csfl(j);  csfl(j) = t
            i=i+1
            j=j-1
        ENDDO
        IF ( first < i-1 ) CALL quicksort_csf(csfl, first, i-1)
        IF ( j+1 < last )  CALL quicksort_csf(csfl, j+1, last)
    END SUBROUTINE quicksort_csf


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Grows the CSF array in RTM exponentially when it is full. During that,
!> the ray canopy sink factors with common source face and target plant canopy
!> grid cell are merged together so that the size doesn't grow out of control.
!------------------------------------------------------------------------------!
    SUBROUTINE merge_and_grow_csf(newsize)
        INTEGER(iwp), INTENT(in)                :: newsize  !< new array size after grow, must be >= ncsfl
                                                            !< or -1 to shrink to minimum
        INTEGER(iwp)                            :: iread, iwrite
        TYPE(t_csf), DIMENSION(:), POINTER      :: acsfnew


        IF ( newsize == -1 )  THEN
!--         merge in-place
            acsfnew => acsf
        ELSE
!--         allocate new array
            IF ( mcsf == 0 )  THEN
                ALLOCATE( acsf1(newsize) )
                acsfnew => acsf1
            ELSE
                ALLOCATE( acsf2(newsize) )
                acsfnew => acsf2
            ENDIF
        ENDIF

        IF ( ncsfl >= 1 )  THEN
!--         sort csf in place (quicksort)
            CALL quicksort_csf(acsf,1,ncsfl)

!--         while moving to a new array, aggregate canopy sink factor records with identical box & source
            acsfnew(1) = acsf(1)
            iwrite = 1
            DO iread = 2, ncsfl
!--             here acsf(kcsf) already has values from acsf(icsf)
                IF ( acsfnew(iwrite)%itx == acsf(iread)%itx &
                         .AND.  acsfnew(iwrite)%ity == acsf(iread)%ity &
                         .AND.  acsfnew(iwrite)%itz == acsf(iread)%itz &
                         .AND.  acsfnew(iwrite)%isurfs == acsf(iread)%isurfs )  THEN

                    acsfnew(iwrite)%rcvf = acsfnew(iwrite)%rcvf + acsf(iread)%rcvf
!--                 advance reading index, keep writing index
                ELSE
!--                 not identical, just advance and copy
                    iwrite = iwrite + 1
                    acsfnew(iwrite) = acsf(iread)
                ENDIF
            ENDDO
            ncsfl = iwrite
        ENDIF

        IF ( newsize == -1 )  THEN
!--         allocate new array and copy shrinked data
            IF ( mcsf == 0 )  THEN
                ALLOCATE( acsf1(ncsfl) )
                acsf1(1:ncsfl) = acsf2(1:ncsfl)
            ELSE
                ALLOCATE( acsf2(ncsfl) )
                acsf2(1:ncsfl) = acsf1(1:ncsfl)
            ENDIF
        ENDIF

!--     deallocate old array
        IF ( mcsf == 0 )  THEN
            mcsf = 1
            acsf => acsf1
            DEALLOCATE( acsf2 )
        ELSE
            mcsf = 0
            acsf => acsf2
            DEALLOCATE( acsf1 )
        ENDIF
        ncsfla = newsize

        IF ( debug_output )  THEN
           WRITE( debug_string, '(A,2I12)' ) 'Grow acsf2:', ncsfl, ncsfla
           CALL debug_message( debug_string, 'info' )
        ENDIF

    END SUBROUTINE merge_and_grow_csf


!-- quicksort.f -*-f90-*-
!-- Author: t-nissie, adaptation J.Resler
!-- License: GPLv3
!-- Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
    RECURSIVE SUBROUTINE quicksort_csf2(kpcsflt, pcsflt, first, last)
        IMPLICIT NONE
        INTEGER(iwp), DIMENSION(:,:), INTENT(INOUT)  :: kpcsflt
        REAL(wp), DIMENSION(:,:), INTENT(INOUT)      :: pcsflt
        INTEGER(iwp), INTENT(IN)                     :: first, last
        REAL(wp), DIMENSION(ndcsf)                   :: t2
        INTEGER(iwp), DIMENSION(kdcsf)               :: x, t1
        INTEGER(iwp)                                 :: i, j

        IF ( first>=last ) RETURN
        x = kpcsflt(:, (first+last)/2 )
        i = first
        j = last
        DO
            DO while ( csf_lt2(kpcsflt(:,i),x) )
                i=i+1
            ENDDO
            DO while ( csf_lt2(x,kpcsflt(:,j)) )
                j=j-1
            ENDDO
            IF ( i >= j ) EXIT
            t1 = kpcsflt(:,i);  kpcsflt(:,i) = kpcsflt(:,j);  kpcsflt(:,j) = t1
            t2 = pcsflt(:,i);  pcsflt(:,i) = pcsflt(:,j);  pcsflt(:,j) = t2
            i=i+1
            j=j-1
        ENDDO
        IF ( first < i-1 ) CALL quicksort_csf2(kpcsflt, pcsflt, first, i-1)
        IF ( j+1 < last )  CALL quicksort_csf2(kpcsflt, pcsflt, j+1, last)
    END SUBROUTINE quicksort_csf2


    PURE FUNCTION csf_lt2(item1, item2) result(res)
        INTEGER(iwp), DIMENSION(kdcsf), INTENT(in)  :: item1, item2
        LOGICAL                                     :: res
        res = ( (item1(3) < item2(3))                                                        &
             .OR.  (item1(3) == item2(3)  .AND.  item1(2) < item2(2))                            &
             .OR.  (item1(3) == item2(3)  .AND.  item1(2) == item2(2)  .AND.  item1(1) < item2(1)) &
             .OR.  (item1(3) == item2(3)  .AND.  item1(2) == item2(2)  .AND.  item1(1) == item2(1) &
                 .AND.  item1(4) < item2(4)) )
    END FUNCTION csf_lt2

    PURE FUNCTION searchsorted(athresh, val) result(ind)
        REAL(wp), DIMENSION(:), INTENT(IN)  :: athresh
        REAL(wp), INTENT(IN)                :: val
        INTEGER(iwp)                        :: ind
        INTEGER(iwp)                        :: i

        DO i = LBOUND(athresh, 1), UBOUND(athresh, 1)
            IF ( val < athresh(i) ) THEN
                ind = i - 1
                RETURN
            ENDIF
        ENDDO
        ind = UBOUND(athresh, 1)
    END FUNCTION searchsorted


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine for averaging 3D data
!------------------------------------------------------------------------------!
SUBROUTINE radiation_3d_data_averaging( mode, variable )


    USE control_parameters

    USE indices

    USE kinds

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  mode    !<
    CHARACTER (LEN=*) :: variable !<

    LOGICAL      ::  match_lsm !< flag indicating natural-type surface
    LOGICAL      ::  match_usm !< flag indicating urban-type surface

    INTEGER(iwp) ::  i !<
    INTEGER(iwp) ::  j !<
    INTEGER(iwp) ::  k !<
    INTEGER(iwp) ::  l, m !< index of current surface element

    INTEGER(iwp)                                       :: ids, idsint_u, idsint_l, isurf
    CHARACTER(LEN=varnamelength)                       :: var

!-- find the real name of the variable
    ids = -1
    l = -1
    var = TRIM(variable)
    DO i = 0, nd-1
        k = len(TRIM(var))
        j = len(TRIM(dirname(i)))
        IF ( k-j+1 >= 1_iwp ) THEN
           IF ( TRIM(var(k-j+1:k)) == TRIM(dirname(i)) )  THEN
               ids = i
               idsint_u = dirint_u(ids)
               idsint_l = dirint_l(ids)
               var = var(:k-j)
               EXIT
           ENDIF
        ENDIF
    ENDDO
    IF ( ids == -1 )  THEN
        var = TRIM(variable)
    ENDIF

    IF ( mode == 'allocate' )  THEN

       SELECT CASE ( TRIM( var ) )
!--          block of large scale (e.g. RRTMG) radiation output variables
             CASE ( 'rad_net*' )
                IF ( .NOT. ALLOCATED( rad_net_av ) )  THEN
                   ALLOCATE( rad_net_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_net_av = 0.0_wp

             CASE ( 'rad_lw_in*' )
                IF ( .NOT. ALLOCATED( rad_lw_in_xy_av ) )  THEN
                   ALLOCATE( rad_lw_in_xy_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_lw_in_xy_av = 0.0_wp

             CASE ( 'rad_lw_out*' )
                IF ( .NOT. ALLOCATED( rad_lw_out_xy_av ) )  THEN
                   ALLOCATE( rad_lw_out_xy_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_lw_out_xy_av = 0.0_wp

             CASE ( 'rad_sw_in*' )
                IF ( .NOT. ALLOCATED( rad_sw_in_xy_av ) )  THEN
                   ALLOCATE( rad_sw_in_xy_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_sw_in_xy_av = 0.0_wp

             CASE ( 'rad_sw_out*' )
                IF ( .NOT. ALLOCATED( rad_sw_out_xy_av ) )  THEN
                   ALLOCATE( rad_sw_out_xy_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_sw_out_xy_av = 0.0_wp

             CASE ( 'rad_lw_in' )
                IF ( .NOT. ALLOCATED( rad_lw_in_av ) )  THEN
                   ALLOCATE( rad_lw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_lw_in_av = 0.0_wp

             CASE ( 'rad_lw_out' )
                IF ( .NOT. ALLOCATED( rad_lw_out_av ) )  THEN
                   ALLOCATE( rad_lw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_lw_out_av = 0.0_wp

             CASE ( 'rad_lw_cs_hr' )
                IF ( .NOT. ALLOCATED( rad_lw_cs_hr_av ) )  THEN
                   ALLOCATE( rad_lw_cs_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_lw_cs_hr_av = 0.0_wp

             CASE ( 'rad_lw_hr' )
                IF ( .NOT. ALLOCATED( rad_lw_hr_av ) )  THEN
                   ALLOCATE( rad_lw_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_lw_hr_av = 0.0_wp

             CASE ( 'rad_sw_in' )
                IF ( .NOT. ALLOCATED( rad_sw_in_av ) )  THEN
                   ALLOCATE( rad_sw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_sw_in_av = 0.0_wp

             CASE ( 'rad_sw_out' )
                IF ( .NOT. ALLOCATED( rad_sw_out_av ) )  THEN
                   ALLOCATE( rad_sw_out_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_sw_out_av = 0.0_wp

             CASE ( 'rad_sw_cs_hr' )
                IF ( .NOT. ALLOCATED( rad_sw_cs_hr_av ) )  THEN
                   ALLOCATE( rad_sw_cs_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_sw_cs_hr_av = 0.0_wp

             CASE ( 'rad_sw_hr' )
                IF ( .NOT. ALLOCATED( rad_sw_hr_av ) )  THEN
                   ALLOCATE( rad_sw_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_sw_hr_av = 0.0_wp

!--          block of RTM output variables
             CASE ( 'rtm_rad_net' )
!--              array of complete radiation balance
                 IF ( .NOT.  ALLOCATED(surfradnet_av) )  THEN
                     ALLOCATE( surfradnet_av(nsurfl) )
                     surfradnet_av = 0.0_wp
                 ENDIF

             CASE ( 'rtm_rad_insw' )
!--                 array of sw radiation falling to surface after i-th reflection
                 IF ( .NOT.  ALLOCATED(surfinsw_av) )  THEN
                     ALLOCATE( surfinsw_av(nsurfl) )
                     surfinsw_av = 0.0_wp
                 ENDIF

             CASE ( 'rtm_rad_inlw' )
!--                 array of lw radiation falling to surface after i-th reflection
                 IF ( .NOT.  ALLOCATED(surfinlw_av) )  THEN
                     ALLOCATE( surfinlw_av(nsurfl) )
                     surfinlw_av = 0.0_wp
                 ENDIF

             CASE ( 'rtm_rad_inswdir' )
!--                 array of direct sw radiation falling to surface from sun
                 IF ( .NOT.  ALLOCATED(surfinswdir_av) )  THEN
                     ALLOCATE( surfinswdir_av(nsurfl) )
                     surfinswdir_av = 0.0_wp
                 ENDIF

             CASE ( 'rtm_rad_inswdif' )
!--                 array of difusion sw radiation falling to surface from sky and borders of the domain
                 IF ( .NOT.  ALLOCATED(surfinswdif_av) )  THEN
                     ALLOCATE( surfinswdif_av(nsurfl) )
                     surfinswdif_av = 0.0_wp
                 ENDIF

             CASE ( 'rtm_rad_inswref' )
!--                 array of sw radiation falling to surface from reflections
                 IF ( .NOT.  ALLOCATED(surfinswref_av) )  THEN
                     ALLOCATE( surfinswref_av(nsurfl) )
                     surfinswref_av = 0.0_wp
                 ENDIF

             CASE ( 'rtm_rad_inlwdif' )
!--                 array of sw radiation falling to surface after i-th reflection
                IF ( .NOT.  ALLOCATED(surfinlwdif_av) )  THEN
                     ALLOCATE( surfinlwdif_av(nsurfl) )
                     surfinlwdif_av = 0.0_wp
                 ENDIF

             CASE ( 'rtm_rad_inlwref' )
!--                 array of lw radiation falling to surface from reflections
                 IF ( .NOT.  ALLOCATED(surfinlwref_av) )  THEN
                     ALLOCATE( surfinlwref_av(nsurfl) )
                     surfinlwref_av = 0.0_wp
                 ENDIF

             CASE ( 'rtm_rad_outsw' )
!--                 array of sw radiation emitted from surface after i-th reflection
                 IF ( .NOT.  ALLOCATED(surfoutsw_av) )  THEN
                     ALLOCATE( surfoutsw_av(nsurfl) )
                     surfoutsw_av = 0.0_wp
                 ENDIF

             CASE ( 'rtm_rad_outlw' )
!--                 array of lw radiation emitted from surface after i-th reflection
                 IF ( .NOT.  ALLOCATED(surfoutlw_av) )  THEN
                     ALLOCATE( surfoutlw_av(nsurfl) )
                     surfoutlw_av = 0.0_wp
                 ENDIF
             CASE ( 'rtm_rad_ressw' )
!--                 array of residua of sw radiation absorbed in surface after last reflection
                 IF ( .NOT.  ALLOCATED(surfins_av) )  THEN
                     ALLOCATE( surfins_av(nsurfl) )
                     surfins_av = 0.0_wp
                 ENDIF

             CASE ( 'rtm_rad_reslw' )
!--                 array of residua of lw radiation absorbed in surface after last reflection
                 IF ( .NOT.  ALLOCATED(surfinl_av) )  THEN
                     ALLOCATE( surfinl_av(nsurfl) )
                     surfinl_av = 0.0_wp
                 ENDIF

             CASE ( 'rtm_rad_pc_inlw' )
!--                 array of of lw radiation absorbed in plant canopy
                 IF ( .NOT.  ALLOCATED(pcbinlw_av) )  THEN
                     ALLOCATE( pcbinlw_av(1:npcbl) )
                     pcbinlw_av = 0.0_wp
                 ENDIF

             CASE ( 'rtm_rad_pc_insw' )
!--                 array of of sw radiation absorbed in plant canopy
                 IF ( .NOT.  ALLOCATED(pcbinsw_av) )  THEN
                     ALLOCATE( pcbinsw_av(1:npcbl) )
                     pcbinsw_av = 0.0_wp
                 ENDIF

             CASE ( 'rtm_rad_pc_inswdir' )
!--                 array of of direct sw radiation absorbed in plant canopy
                 IF ( .NOT.  ALLOCATED(pcbinswdir_av) )  THEN
                     ALLOCATE( pcbinswdir_av(1:npcbl) )
                     pcbinswdir_av = 0.0_wp
                 ENDIF

             CASE ( 'rtm_rad_pc_inswdif' )
!--                 array of of diffuse sw radiation absorbed in plant canopy
                 IF ( .NOT.  ALLOCATED(pcbinswdif_av) )  THEN
                     ALLOCATE( pcbinswdif_av(1:npcbl) )
                     pcbinswdif_av = 0.0_wp
                 ENDIF

             CASE ( 'rtm_rad_pc_inswref' )
!--                 array of of reflected sw radiation absorbed in plant canopy
                 IF ( .NOT.  ALLOCATED(pcbinswref_av) )  THEN
                     ALLOCATE( pcbinswref_av(1:npcbl) )
                     pcbinswref_av = 0.0_wp
                 ENDIF

             CASE ( 'rtm_mrt_sw' )
                IF ( .NOT. ALLOCATED( mrtinsw_av ) )  THEN
                   ALLOCATE( mrtinsw_av(nmrtbl) )
                ENDIF
                mrtinsw_av = 0.0_wp

             CASE ( 'rtm_mrt_lw' )
                IF ( .NOT. ALLOCATED( mrtinlw_av ) )  THEN
                   ALLOCATE( mrtinlw_av(nmrtbl) )
                ENDIF
                mrtinlw_av = 0.0_wp

             CASE ( 'rtm_mrt' )
                IF ( .NOT. ALLOCATED( mrt_av ) )  THEN
                   ALLOCATE( mrt_av(nmrtbl) )
                ENDIF
                mrt_av = 0.0_wp

          CASE DEFAULT
             CONTINUE

       END SELECT

    ELSEIF ( mode == 'sum' )  THEN

       SELECT CASE ( TRIM( var ) )
!--       block of large scale (e.g. RRTMG) radiation output variables
          CASE ( 'rad_net*' )
             IF ( ALLOCATED( rad_net_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      match_lsm = surf_lsm_h%start_index(j,i) <=               &
                                  surf_lsm_h%end_index(j,i)
                      match_usm = surf_usm_h%start_index(j,i) <=               &
                                  surf_usm_h%end_index(j,i)

                     IF ( match_lsm  .AND.  .NOT. match_usm )  THEN
                         m = surf_lsm_h%end_index(j,i)
                         rad_net_av(j,i) = rad_net_av(j,i) +                   &
                                         surf_lsm_h%rad_net(m)
                      ELSEIF ( match_usm )  THEN
                         m = surf_usm_h%end_index(j,i)
                         rad_net_av(j,i) = rad_net_av(j,i) +                   &
                                         surf_usm_h%rad_net(m)
                      ENDIF
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_lw_in*' )
             IF ( ALLOCATED( rad_lw_in_xy_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      match_lsm = surf_lsm_h%start_index(j,i) <=               &
                                  surf_lsm_h%end_index(j,i)
                      match_usm = surf_usm_h%start_index(j,i) <=               &
                                  surf_usm_h%end_index(j,i)

                      IF ( match_lsm  .AND.  .NOT. match_usm )  THEN
                         m = surf_lsm_h%end_index(j,i)
                         rad_lw_in_xy_av(j,i) = rad_lw_in_xy_av(j,i) +         &
                                         surf_lsm_h%rad_lw_in(m)
                      ELSEIF ( match_usm )  THEN
                         m = surf_usm_h%end_index(j,i)
                         rad_lw_in_xy_av(j,i) = rad_lw_in_xy_av(j,i) +         &
                                         surf_usm_h%rad_lw_in(m)
                      ENDIF
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_lw_out*' )
             IF ( ALLOCATED( rad_lw_out_xy_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      match_lsm = surf_lsm_h%start_index(j,i) <=               &
                                  surf_lsm_h%end_index(j,i)
                      match_usm = surf_usm_h%start_index(j,i) <=               &
                                  surf_usm_h%end_index(j,i)

                      IF ( match_lsm  .AND.  .NOT. match_usm )  THEN
                         m = surf_lsm_h%end_index(j,i)
                         rad_lw_out_xy_av(j,i) = rad_lw_out_xy_av(j,i) +       &
                                                 surf_lsm_h%rad_lw_out(m)
                      ELSEIF ( match_usm )  THEN
                         m = surf_usm_h%end_index(j,i)
                         rad_lw_out_xy_av(j,i) = rad_lw_out_xy_av(j,i) +       &
                                                 surf_usm_h%rad_lw_out(m)
                      ENDIF
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_sw_in*' )
             IF ( ALLOCATED( rad_sw_in_xy_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      match_lsm = surf_lsm_h%start_index(j,i) <=               &
                                  surf_lsm_h%end_index(j,i)
                      match_usm = surf_usm_h%start_index(j,i) <=               &
                                  surf_usm_h%end_index(j,i)

                      IF ( match_lsm  .AND.  .NOT. match_usm )  THEN
                         m = surf_lsm_h%end_index(j,i)
                         rad_sw_in_xy_av(j,i) = rad_sw_in_xy_av(j,i) +         &
                                                surf_lsm_h%rad_sw_in(m)
                      ELSEIF ( match_usm )  THEN
                         m = surf_usm_h%end_index(j,i)
                         rad_sw_in_xy_av(j,i) = rad_sw_in_xy_av(j,i) +         &
                                                surf_usm_h%rad_sw_in(m)
                      ENDIF
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_sw_out*' )
             IF ( ALLOCATED( rad_sw_out_xy_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      match_lsm = surf_lsm_h%start_index(j,i) <=               &
                                  surf_lsm_h%end_index(j,i)
                      match_usm = surf_usm_h%start_index(j,i) <=               &
                                  surf_usm_h%end_index(j,i)

                      IF ( match_lsm  .AND.  .NOT. match_usm )  THEN
                         m = surf_lsm_h%end_index(j,i)
                         rad_sw_out_xy_av(j,i) = rad_sw_out_xy_av(j,i) +       &
                                                 surf_lsm_h%rad_sw_out(m)
                      ELSEIF ( match_usm )  THEN
                         m = surf_usm_h%end_index(j,i)
                         rad_sw_out_xy_av(j,i) = rad_sw_out_xy_av(j,i) +       &
                                                 surf_usm_h%rad_sw_out(m)
                      ENDIF
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_lw_in' )
             IF ( ALLOCATED( rad_lw_in_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_lw_in_av(k,j,i) = rad_lw_in_av(k,j,i)             &
                                               + rad_lw_in(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_lw_out' )
             IF ( ALLOCATED( rad_lw_out_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_lw_out_av(k,j,i) = rad_lw_out_av(k,j,i)           &
                                                + rad_lw_out(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_lw_cs_hr' )
             IF ( ALLOCATED( rad_lw_cs_hr_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_lw_cs_hr_av(k,j,i) = rad_lw_cs_hr_av(k,j,i)       &
                                                  + rad_lw_cs_hr(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_lw_hr' )
             IF ( ALLOCATED( rad_lw_hr_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_lw_hr_av(k,j,i) = rad_lw_hr_av(k,j,i)             &
                                               + rad_lw_hr(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_sw_in' )
             IF ( ALLOCATED( rad_sw_in_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_sw_in_av(k,j,i) = rad_sw_in_av(k,j,i)             &
                                               + rad_sw_in(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_sw_out' )
             IF ( ALLOCATED( rad_sw_out_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_sw_out_av(k,j,i) = rad_sw_out_av(k,j,i)           &
                                                + rad_sw_out(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_sw_cs_hr' )
             IF ( ALLOCATED( rad_sw_cs_hr_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_sw_cs_hr_av(k,j,i) = rad_sw_cs_hr_av(k,j,i)       &
                                                  + rad_sw_cs_hr(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_sw_hr' )
             IF ( ALLOCATED( rad_sw_hr_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_sw_hr_av(k,j,i) = rad_sw_hr_av(k,j,i)             &
                                               + rad_sw_hr(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

!--       block of RTM output variables
          CASE ( 'rtm_rad_net' )
!--           array of complete radiation balance
              DO isurf = dirstart(ids), dirend(ids)
                 IF ( surfl(id,isurf) == idsint_u  .OR.  surfl(id,isurf) == idsint_l )  THEN
                    surfradnet_av(isurf) = surfinsw(isurf) - surfoutsw(isurf) +  surfinlw(isurf) - surfoutlw(isurf)
                 ENDIF
              ENDDO

          CASE ( 'rtm_rad_insw' )
!--           array of sw radiation falling to surface after i-th reflection
              DO isurf = dirstart(ids), dirend(ids)
                  IF ( surfl(id,isurf) == idsint_u  .OR.  surfl(id,isurf) == idsint_l )  THEN
                      surfinsw_av(isurf) = surfinsw_av(isurf) + surfinsw(isurf)
                  ENDIF
              ENDDO

          CASE ( 'rtm_rad_inlw' )
!--           array of lw radiation falling to surface after i-th reflection
              DO isurf = dirstart(ids), dirend(ids)
                  IF ( surfl(id,isurf) == idsint_u  .OR.  surfl(id,isurf) == idsint_l )  THEN
                      surfinlw_av(isurf) = surfinlw_av(isurf) + surfinlw(isurf)
                  ENDIF
              ENDDO

          CASE ( 'rtm_rad_inswdir' )
!--           array of direct sw radiation falling to surface from sun
              DO isurf = dirstart(ids), dirend(ids)
                  IF ( surfl(id,isurf) == idsint_u  .OR.  surfl(id,isurf) == idsint_l )  THEN
                      surfinswdir_av(isurf) = surfinswdir_av(isurf) + surfinswdir(isurf)
                  ENDIF
              ENDDO

          CASE ( 'rtm_rad_inswdif' )
!--           array of difusion sw radiation falling to surface from sky and borders of the domain
              DO isurf = dirstart(ids), dirend(ids)
                  IF ( surfl(id,isurf) == idsint_u  .OR.  surfl(id,isurf) == idsint_l )  THEN
                      surfinswdif_av(isurf) = surfinswdif_av(isurf) + surfinswdif(isurf)
                  ENDIF
              ENDDO

          CASE ( 'rtm_rad_inswref' )
!--           array of sw radiation falling to surface from reflections
              DO isurf = dirstart(ids), dirend(ids)
                  IF ( surfl(id,isurf) == idsint_u  .OR.  surfl(id,isurf) == idsint_l )  THEN
                      surfinswref_av(isurf) = surfinswref_av(isurf) + surfinsw(isurf) - &
                                          surfinswdir(isurf) - surfinswdif(isurf)
                  ENDIF
              ENDDO


          CASE ( 'rtm_rad_inlwdif' )
!--           array of sw radiation falling to surface after i-th reflection
              DO isurf = dirstart(ids), dirend(ids)
                  IF ( surfl(id,isurf) == idsint_u  .OR.  surfl(id,isurf) == idsint_l )  THEN
                      surfinlwdif_av(isurf) = surfinlwdif_av(isurf) + surfinlwdif(isurf)
                  ENDIF
              ENDDO
!
          CASE ( 'rtm_rad_inlwref' )
!--           array of lw radiation falling to surface from reflections
              DO isurf = dirstart(ids), dirend(ids)
                  IF ( surfl(id,isurf) == idsint_u  .OR.  surfl(id,isurf) == idsint_l )  THEN
                      surfinlwref_av(isurf) = surfinlwref_av(isurf) + &
                                          surfinlw(isurf) - surfinlwdif(isurf)
                  ENDIF
              ENDDO

          CASE ( 'rtm_rad_outsw' )
!--           array of sw radiation emitted from surface after i-th reflection
              DO isurf = dirstart(ids), dirend(ids)
                  IF ( surfl(id,isurf) == idsint_u  .OR.  surfl(id,isurf) == idsint_l )  THEN
                      surfoutsw_av(isurf) = surfoutsw_av(isurf) + surfoutsw(isurf)
                  ENDIF
              ENDDO

          CASE ( 'rtm_rad_outlw' )
!--           array of lw radiation emitted from surface after i-th reflection
              DO isurf = dirstart(ids), dirend(ids)
                  IF ( surfl(id,isurf) == idsint_u  .OR.  surfl(id,isurf) == idsint_l )  THEN
                      surfoutlw_av(isurf) = surfoutlw_av(isurf) + surfoutlw(isurf)
                  ENDIF
              ENDDO

          CASE ( 'rtm_rad_ressw' )
!--           array of residua of sw radiation absorbed in surface after last reflection
              DO isurf = dirstart(ids), dirend(ids)
                  IF ( surfl(id,isurf) == idsint_u  .OR.  surfl(id,isurf) == idsint_l )  THEN
                      surfins_av(isurf) = surfins_av(isurf) + surfins(isurf)
                  ENDIF
              ENDDO

          CASE ( 'rtm_rad_reslw' )
!--           array of residua of lw radiation absorbed in surface after last reflection
              DO isurf = dirstart(ids), dirend(ids)
                  IF ( surfl(id,isurf) == idsint_u  .OR.  surfl(id,isurf) == idsint_l )  THEN
                      surfinl_av(isurf) = surfinl_av(isurf) + surfinl(isurf)
                  ENDIF
              ENDDO

          CASE ( 'rtm_rad_pc_inlw' )
              DO l = 1, npcbl
                 pcbinlw_av(l) = pcbinlw_av(l) + pcbinlw(l)
              ENDDO

          CASE ( 'rtm_rad_pc_insw' )
              DO l = 1, npcbl
                 pcbinsw_av(l) = pcbinsw_av(l) + pcbinsw(l)
              ENDDO

          CASE ( 'rtm_rad_pc_inswdir' )
              DO l = 1, npcbl
                 pcbinswdir_av(l) = pcbinswdir_av(l) + pcbinswdir(l)
              ENDDO

          CASE ( 'rtm_rad_pc_inswdif' )
              DO l = 1, npcbl
                 pcbinswdif_av(l) = pcbinswdif_av(l) + pcbinswdif(l)
              ENDDO

          CASE ( 'rtm_rad_pc_inswref' )
              DO l = 1, npcbl
                 pcbinswref_av(l) = pcbinswref_av(l) + pcbinsw(l) - pcbinswdir(l) - pcbinswdif(l)
              ENDDO

          CASE ( 'rad_mrt_sw' )
             IF ( ALLOCATED( mrtinsw_av ) )  THEN
                mrtinsw_av(:) = mrtinsw_av(:) + mrtinsw(:)
             ENDIF

          CASE ( 'rad_mrt_lw' )
             IF ( ALLOCATED( mrtinlw_av ) )  THEN
                mrtinlw_av(:) = mrtinlw_av(:) + mrtinlw(:)
             ENDIF

          CASE ( 'rad_mrt' )
             IF ( ALLOCATED( mrt_av ) )  THEN
                mrt_av(:) = mrt_av(:) + mrt(:)
             ENDIF

          CASE DEFAULT
             CONTINUE

       END SELECT

    ELSEIF ( mode == 'average' )  THEN

       SELECT CASE ( TRIM( var ) )
!--       block of large scale (e.g. RRTMG) radiation output variables
          CASE ( 'rad_net*' )
             IF ( ALLOCATED( rad_net_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      rad_net_av(j,i) = rad_net_av(j,i)                        &
                                        / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_lw_in*' )
             IF ( ALLOCATED( rad_lw_in_xy_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      rad_lw_in_xy_av(j,i) = rad_lw_in_xy_av(j,i)              &
                                        / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_lw_out*' )
             IF ( ALLOCATED( rad_lw_out_xy_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      rad_lw_out_xy_av(j,i) = rad_lw_out_xy_av(j,i)            &
                                        / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_sw_in*' )
             IF ( ALLOCATED( rad_sw_in_xy_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      rad_sw_in_xy_av(j,i) = rad_sw_in_xy_av(j,i)              &
                                        / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_sw_out*' )
             IF ( ALLOCATED( rad_sw_out_xy_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      rad_sw_out_xy_av(j,i) = rad_sw_out_xy_av(j,i)             &
                                        / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_lw_in' )
             IF ( ALLOCATED( rad_lw_in_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_lw_in_av(k,j,i) = rad_lw_in_av(k,j,i)             &
                                               / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_lw_out' )
             IF ( ALLOCATED( rad_lw_out_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_lw_out_av(k,j,i) = rad_lw_out_av(k,j,i)           &
                                                / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_lw_cs_hr' )
             IF ( ALLOCATED( rad_lw_cs_hr_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_lw_cs_hr_av(k,j,i) = rad_lw_cs_hr_av(k,j,i)       &
                                                / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_lw_hr' )
             IF ( ALLOCATED( rad_lw_hr_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_lw_hr_av(k,j,i) = rad_lw_hr_av(k,j,i)             &
                                               / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_sw_in' )
             IF ( ALLOCATED( rad_sw_in_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_sw_in_av(k,j,i) = rad_sw_in_av(k,j,i)             &
                                               / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_sw_out' )
             IF ( ALLOCATED( rad_sw_out_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_sw_out_av(k,j,i) = rad_sw_out_av(k,j,i)           &
                                                / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_sw_cs_hr' )
             IF ( ALLOCATED( rad_sw_cs_hr_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_sw_cs_hr_av(k,j,i) = rad_sw_cs_hr_av(k,j,i)       &
                                                / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_sw_hr' )
             IF ( ALLOCATED( rad_sw_hr_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_sw_hr_av(k,j,i) = rad_sw_hr_av(k,j,i)             &
                                               / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

!--       block of RTM output variables
          CASE ( 'rtm_rad_net' )
!--           array of complete radiation balance
              DO isurf = dirstart(ids), dirend(ids)
                  IF ( surfl(id,isurf) == idsint_u  .OR.  surfl(id,isurf) == idsint_l )  THEN
                      surfradnet_av(isurf) = surfinsw_av(isurf) / REAL( average_count_3d, kind=wp )
                  ENDIF
              ENDDO

          CASE ( 'rtm_rad_insw' )
!--           array of sw radiation falling to surface after i-th reflection
              DO isurf = dirstart(ids), dirend(ids)
                  IF ( surfl(id,isurf) == idsint_u  .OR.  surfl(id,isurf) == idsint_l )  THEN
                      surfinsw_av(isurf) = surfinsw_av(isurf) / REAL( average_count_3d, kind=wp )
                  ENDIF
              ENDDO

          CASE ( 'rtm_rad_inlw' )
!--           array of lw radiation falling to surface after i-th reflection
              DO isurf = dirstart(ids), dirend(ids)
                  IF ( surfl(id,isurf) == idsint_u  .OR.  surfl(id,isurf) == idsint_l )  THEN
                      surfinlw_av(isurf) = surfinlw_av(isurf) / REAL( average_count_3d, kind=wp )
                  ENDIF
              ENDDO

          CASE ( 'rtm_rad_inswdir' )
!--           array of direct sw radiation falling to surface from sun
              DO isurf = dirstart(ids), dirend(ids)
                  IF ( surfl(id,isurf) == idsint_u  .OR.  surfl(id,isurf) == idsint_l )  THEN
                      surfinswdir_av(isurf) = surfinswdir_av(isurf) / REAL( average_count_3d, kind=wp )
                  ENDIF
              ENDDO

          CASE ( 'rtm_rad_inswdif' )
!--           array of difusion sw radiation falling to surface from sky and borders of the domain
              DO isurf = dirstart(ids), dirend(ids)
                  IF ( surfl(id,isurf) == idsint_u  .OR.  surfl(id,isurf) == idsint_l )  THEN
                      surfinswdif_av(isurf) = surfinswdif_av(isurf) / REAL( average_count_3d, kind=wp )
                  ENDIF
              ENDDO

          CASE ( 'rtm_rad_inswref' )
!--           array of sw radiation falling to surface from reflections
              DO isurf = dirstart(ids), dirend(ids)
                  IF ( surfl(id,isurf) == idsint_u  .OR.  surfl(id,isurf) == idsint_l )  THEN
                      surfinswref_av(isurf) = surfinswref_av(isurf) / REAL( average_count_3d, kind=wp )
                  ENDIF
              ENDDO

          CASE ( 'rtm_rad_inlwdif' )
!--           array of sw radiation falling to surface after i-th reflection
              DO isurf = dirstart(ids), dirend(ids)
                  IF ( surfl(id,isurf) == idsint_u  .OR.  surfl(id,isurf) == idsint_l )  THEN
                      surfinlwdif_av(isurf) = surfinlwdif_av(isurf) / REAL( average_count_3d, kind=wp )
                  ENDIF
              ENDDO

          CASE ( 'rtm_rad_inlwref' )
!--           array of lw radiation falling to surface from reflections
              DO isurf = dirstart(ids), dirend(ids)
                  IF ( surfl(id,isurf) == idsint_u  .OR.  surfl(id,isurf) == idsint_l )  THEN
                      surfinlwref_av(isurf) = surfinlwref_av(isurf) / REAL( average_count_3d, kind=wp )
                  ENDIF
              ENDDO

          CASE ( 'rtm_rad_outsw' )
!--           array of sw radiation emitted from surface after i-th reflection
              DO isurf = dirstart(ids), dirend(ids)
                  IF ( surfl(id,isurf) == idsint_u  .OR.  surfl(id,isurf) == idsint_l )  THEN
                      surfoutsw_av(isurf) = surfoutsw_av(isurf) / REAL( average_count_3d, kind=wp )
                  ENDIF
              ENDDO

          CASE ( 'rtm_rad_outlw' )
!--           array of lw radiation emitted from surface after i-th reflection
              DO isurf = dirstart(ids), dirend(ids)
                  IF ( surfl(id,isurf) == idsint_u  .OR.  surfl(id,isurf) == idsint_l )  THEN
                      surfoutlw_av(isurf) = surfoutlw_av(isurf) / REAL( average_count_3d, kind=wp )
                  ENDIF
              ENDDO

          CASE ( 'rtm_rad_ressw' )
!--           array of residua of sw radiation absorbed in surface after last reflection
              DO isurf = dirstart(ids), dirend(ids)
                  IF ( surfl(id,isurf) == idsint_u  .OR.  surfl(id,isurf) == idsint_l )  THEN
                      surfins_av(isurf) = surfins_av(isurf) / REAL( average_count_3d, kind=wp )
                  ENDIF
              ENDDO

          CASE ( 'rtm_rad_reslw' )
!--           array of residua of lw radiation absorbed in surface after last reflection
              DO isurf = dirstart(ids), dirend(ids)
                  IF ( surfl(id,isurf) == idsint_u  .OR.  surfl(id,isurf) == idsint_l )  THEN
                      surfinl_av(isurf) = surfinl_av(isurf) / REAL( average_count_3d, kind=wp )
                  ENDIF
              ENDDO

          CASE ( 'rtm_rad_pc_inlw' )
              DO l = 1, npcbl
                 pcbinlw_av(:) = pcbinlw_av(:) / REAL( average_count_3d, kind=wp )
              ENDDO

          CASE ( 'rtm_rad_pc_insw' )
              DO l = 1, npcbl
                 pcbinsw_av(:) = pcbinsw_av(:) / REAL( average_count_3d, kind=wp )
              ENDDO

          CASE ( 'rtm_rad_pc_inswdir' )
              DO l = 1, npcbl
                 pcbinswdir_av(:) = pcbinswdir_av(:) / REAL( average_count_3d, kind=wp )
              ENDDO

          CASE ( 'rtm_rad_pc_inswdif' )
              DO l = 1, npcbl
                 pcbinswdif_av(:) = pcbinswdif_av(:) / REAL( average_count_3d, kind=wp )
              ENDDO

          CASE ( 'rtm_rad_pc_inswref' )
              DO l = 1, npcbl
                 pcbinswref_av(:) = pcbinswref_av(:) / REAL( average_count_3d, kind=wp )
              ENDDO

          CASE ( 'rad_mrt_lw' )
             IF ( ALLOCATED( mrtinlw_av ) )  THEN
                mrtinlw_av(:) = mrtinlw_av(:) / REAL( average_count_3d, KIND=wp )
             ENDIF

          CASE ( 'rad_mrt' )
             IF ( ALLOCATED( mrt_av ) )  THEN
                mrt_av(:) = mrt_av(:) / REAL( average_count_3d, KIND=wp )
             ENDIF

       END SELECT

    ENDIF

END SUBROUTINE radiation_3d_data_averaging


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining appropriate grid for netcdf variables.
!> It is called out from subroutine netcdf.
!------------------------------------------------------------------------------!
SUBROUTINE radiation_define_netcdf_grid( variable, found, grid_x, grid_y, grid_z )

    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT(IN)  ::  variable    !<
    LOGICAL, INTENT(OUT)           ::  found       !<
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_x      !<
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_y      !<
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_z      !<

    CHARACTER (len=varnamelength)  :: var

    found  = .TRUE.

!
!-- Check for the grid
    var = TRIM(variable)
!-- RTM directional variables
    IF ( var(1:12) == 'rtm_rad_net_'  .OR.  var(1:13) == 'rtm_rad_insw_'  .OR.          &
         var(1:13) == 'rtm_rad_inlw_'  .OR.  var(1:16) == 'rtm_rad_inswdir_'  .OR.      &
         var(1:16) == 'rtm_rad_inswdif_'  .OR.  var(1:16) == 'rtm_rad_inswref_'  .OR.   &
         var(1:16) == 'rtm_rad_inlwdif_'  .OR.  var(1:16) == 'rtm_rad_inlwref_'  .OR.   &
         var(1:14) == 'rtm_rad_outsw_'  .OR.  var(1:14) == 'rtm_rad_outlw_'  .OR.       &
         var(1:14) == 'rtm_rad_ressw_'  .OR.  var(1:14) == 'rtm_rad_reslw_'  .OR.       &
         var == 'rtm_rad_pc_inlw'  .OR.                                                 &
         var == 'rtm_rad_pc_insw'  .OR.  var == 'rtm_rad_pc_inswdir'  .OR.              &
         var == 'rtm_rad_pc_inswdif'  .OR.  var == 'rtm_rad_pc_inswref'  .OR.           &
         var(1:7) == 'rtm_svf'  .OR.  var(1:7) == 'rtm_dif'  .OR.                       &
         var(1:9) == 'rtm_skyvf' .OR. var(1:10) == 'rtm_skyvft'  .OR.                   &
         var(1:12) == 'rtm_surfalb_'  .OR.  var(1:13) == 'rtm_surfemis_'  .OR.          &
         var == 'rtm_mrt'  .OR.  var ==  'rtm_mrt_sw'  .OR.  var == 'rtm_mrt_lw' )  THEN

         found = .TRUE.
         grid_x = 'x'
         grid_y = 'y'
         grid_z = 'zu'
    ELSE

       SELECT CASE ( TRIM( var ) )

          CASE ( 'rad_lw_cs_hr', 'rad_lw_hr', 'rad_sw_cs_hr', 'rad_sw_hr',        &
                 'rad_lw_cs_hr_xy', 'rad_lw_hr_xy', 'rad_sw_cs_hr_xy',            &
                 'rad_sw_hr_xy', 'rad_lw_cs_hr_xz', 'rad_lw_hr_xz',               &
                 'rad_sw_cs_hr_xz', 'rad_sw_hr_xz', 'rad_lw_cs_hr_yz',            &
                 'rad_lw_hr_yz', 'rad_sw_cs_hr_yz', 'rad_sw_hr_yz',               &
                 'rad_mrt', 'rad_mrt_sw', 'rad_mrt_lw' )
             grid_x = 'x'
             grid_y = 'y'
             grid_z = 'zu'

          CASE ( 'rad_lw_in', 'rad_lw_out', 'rad_sw_in', 'rad_sw_out',            &
                 'rad_lw_in_xy', 'rad_lw_out_xy', 'rad_sw_in_xy','rad_sw_out_xy', &
                 'rad_lw_in_xz', 'rad_lw_out_xz', 'rad_sw_in_xz','rad_sw_out_xz', &
                 'rad_lw_in_yz', 'rad_lw_out_yz', 'rad_sw_in_yz','rad_sw_out_yz' )
             grid_x = 'x'
             grid_y = 'y'
             grid_z = 'zw'


          CASE DEFAULT
             found  = .FALSE.
             grid_x = 'none'
             grid_y = 'none'
             grid_z = 'none'

           END SELECT
       ENDIF

    END SUBROUTINE radiation_define_netcdf_grid

!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining 2D output variables
!------------------------------------------------------------------------------!
 SUBROUTINE radiation_data_output_2d( av, variable, found, grid, mode,         &
                                      local_pf, two_d, nzb_do, nzt_do )

    USE indices

    USE kinds


    IMPLICIT NONE

    CHARACTER (LEN=*) ::  grid     !<
    CHARACTER (LEN=*) ::  mode     !<
    CHARACTER (LEN=*) ::  variable !<

    INTEGER(iwp) ::  av !<
    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  k  !<
    INTEGER(iwp) ::  m  !< index of surface element at grid point (j,i)
    INTEGER(iwp) ::  nzb_do   !<
    INTEGER(iwp) ::  nzt_do   !<

    LOGICAL      ::  found !<
    LOGICAL      ::  two_d !< flag parameter that indicates 2D variables (horizontal cross sections)

    REAL(wp) ::  fill_value = -999.0_wp    !< value for the _FillValue attribute

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf !<

    found = .TRUE.

    SELECT CASE ( TRIM( variable ) )

       CASE ( 'rad_net*_xy' )        ! 2d-array
          IF ( av == 0 ) THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
!
!--                Obtain rad_net from its respective surface type
!--                Natural-type surfaces
                   DO  m = surf_lsm_h%start_index(j,i),                        &
                           surf_lsm_h%end_index(j,i)
                      local_pf(i,j,nzb+1) = surf_lsm_h%rad_net(m)
                   ENDDO
!
!--                Urban-type surfaces
                   DO  m = surf_usm_h%start_index(j,i),                        &
                           surf_usm_h%end_index(j,i)
                      local_pf(i,j,nzb+1) = surf_usm_h%rad_net(m)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( rad_net_av ) ) THEN
                ALLOCATE( rad_net_av(nysg:nyng,nxlg:nxrg) )
                rad_net_av = REAL( fill_value, KIND = wp )
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = rad_net_av(j,i)
                ENDDO
             ENDDO
          ENDIF
          two_d = .TRUE.
          grid = 'zu1'

       CASE ( 'rad_lw_in*_xy' )        ! 2d-array
          IF ( av == 0 ) THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
!
!--                Obtain rad_net from its respective surface type
!--                Natural-type surfaces
                   DO  m = surf_lsm_h%start_index(j,i),                        &
                           surf_lsm_h%end_index(j,i)
                      local_pf(i,j,nzb+1) = surf_lsm_h%rad_lw_in(m)
                   ENDDO
!
!--                Urban-type surfaces
                   DO  m = surf_usm_h%start_index(j,i),                        &
                           surf_usm_h%end_index(j,i)
                      local_pf(i,j,nzb+1) = surf_usm_h%rad_lw_in(m)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( rad_lw_in_xy_av ) ) THEN
                ALLOCATE( rad_lw_in_xy_av(nysg:nyng,nxlg:nxrg) )
                rad_lw_in_xy_av = REAL( fill_value, KIND = wp )
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = rad_lw_in_xy_av(j,i)
                ENDDO
             ENDDO
          ENDIF
          two_d = .TRUE.
          grid = 'zu1'

       CASE ( 'rad_lw_out*_xy' )        ! 2d-array
          IF ( av == 0 ) THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
!
!--                Obtain rad_net from its respective surface type
!--                Natural-type surfaces
                   DO  m = surf_lsm_h%start_index(j,i),                        &
                           surf_lsm_h%end_index(j,i)
                      local_pf(i,j,nzb+1) = surf_lsm_h%rad_lw_out(m)
                   ENDDO
!
!--                Urban-type surfaces
                   DO  m = surf_usm_h%start_index(j,i),                        &
                           surf_usm_h%end_index(j,i)
                      local_pf(i,j,nzb+1) = surf_usm_h%rad_lw_out(m)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( rad_lw_out_xy_av ) ) THEN
                ALLOCATE( rad_lw_out_xy_av(nysg:nyng,nxlg:nxrg) )
                rad_lw_out_xy_av = REAL( fill_value, KIND = wp )
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = rad_lw_out_xy_av(j,i)
                ENDDO
             ENDDO
          ENDIF
          two_d = .TRUE.
          grid = 'zu1'

       CASE ( 'rad_sw_in*_xy' )        ! 2d-array
          IF ( av == 0 ) THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
!
!--                Obtain rad_net from its respective surface type
!--                Natural-type surfaces
                   DO  m = surf_lsm_h%start_index(j,i),                        &
                           surf_lsm_h%end_index(j,i)
                      local_pf(i,j,nzb+1) = surf_lsm_h%rad_sw_in(m)
                   ENDDO
!
!--                Urban-type surfaces
                   DO  m = surf_usm_h%start_index(j,i),                        &
                           surf_usm_h%end_index(j,i)
                      local_pf(i,j,nzb+1) = surf_usm_h%rad_sw_in(m)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( rad_sw_in_xy_av ) ) THEN
                ALLOCATE( rad_sw_in_xy_av(nysg:nyng,nxlg:nxrg) )
                rad_sw_in_xy_av = REAL( fill_value, KIND = wp )
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = rad_sw_in_xy_av(j,i)
                ENDDO
             ENDDO
          ENDIF
          two_d = .TRUE.
          grid = 'zu1'

       CASE ( 'rad_sw_out*_xy' )        ! 2d-array
          IF ( av == 0 ) THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
!
!--                Obtain rad_net from its respective surface type
!--                Natural-type surfaces
                   DO  m = surf_lsm_h%start_index(j,i),                        &
                           surf_lsm_h%end_index(j,i)
                      local_pf(i,j,nzb+1) = surf_lsm_h%rad_sw_out(m)
                   ENDDO
!
!--                Urban-type surfaces
                   DO  m = surf_usm_h%start_index(j,i),                        &
                           surf_usm_h%end_index(j,i)
                      local_pf(i,j,nzb+1) = surf_usm_h%rad_sw_out(m)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( rad_sw_out_xy_av ) ) THEN
                ALLOCATE( rad_sw_out_xy_av(nysg:nyng,nxlg:nxrg) )
                rad_sw_out_xy_av = REAL( fill_value, KIND = wp )
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = rad_sw_out_xy_av(j,i)
                ENDDO
             ENDDO
          ENDIF
          two_d = .TRUE.
          grid = 'zu1'

       CASE ( 'rad_lw_in_xy', 'rad_lw_in_xz', 'rad_lw_in_yz' )
          IF ( av == 0 ) THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = rad_lw_in(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
            IF ( .NOT. ALLOCATED( rad_lw_in_av ) ) THEN
               ALLOCATE( rad_lw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_lw_in_av = REAL( fill_value, KIND = wp )
            ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = rad_lw_in_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid = 'zu'

       CASE ( 'rad_lw_out_xy', 'rad_lw_out_xz', 'rad_lw_out_yz' )
          IF ( av == 0 ) THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = rad_lw_out(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
            IF ( .NOT. ALLOCATED( rad_lw_out_av ) ) THEN
               ALLOCATE( rad_lw_out_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_lw_out_av = REAL( fill_value, KIND = wp )
            ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = rad_lw_out_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid = 'zu'

       CASE ( 'rad_lw_cs_hr_xy', 'rad_lw_cs_hr_xz', 'rad_lw_cs_hr_yz' )
          IF ( av == 0 ) THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = rad_lw_cs_hr(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
            IF ( .NOT. ALLOCATED( rad_lw_cs_hr_av ) ) THEN
               ALLOCATE( rad_lw_cs_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_lw_cs_hr_av = REAL( fill_value, KIND = wp )
            ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = rad_lw_cs_hr_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid = 'zw'

       CASE ( 'rad_lw_hr_xy', 'rad_lw_hr_xz', 'rad_lw_hr_yz' )
          IF ( av == 0 ) THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = rad_lw_hr(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
            IF ( .NOT. ALLOCATED( rad_lw_hr_av ) ) THEN
               ALLOCATE( rad_lw_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_lw_hr_av= REAL( fill_value, KIND = wp )
            ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = rad_lw_hr_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid = 'zw'

       CASE ( 'rad_sw_in_xy', 'rad_sw_in_xz', 'rad_sw_in_yz' )
          IF ( av == 0 ) THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = rad_sw_in(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
            IF ( .NOT. ALLOCATED( rad_sw_in_av ) ) THEN
               ALLOCATE( rad_sw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_sw_in_av = REAL( fill_value, KIND = wp )
            ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = rad_sw_in_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid = 'zu'

       CASE ( 'rad_sw_out_xy', 'rad_sw_out_xz', 'rad_sw_out_yz' )
          IF ( av == 0 ) THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = rad_sw_out(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
            IF ( .NOT. ALLOCATED( rad_sw_out_av ) ) THEN
               ALLOCATE( rad_sw_out_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_sw_out_av = REAL( fill_value, KIND = wp )
            ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb, nzt+1
                      local_pf(i,j,k) = rad_sw_out_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid = 'zu'

       CASE ( 'rad_sw_cs_hr_xy', 'rad_sw_cs_hr_xz', 'rad_sw_cs_hr_yz' )
          IF ( av == 0 ) THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = rad_sw_cs_hr(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
            IF ( .NOT. ALLOCATED( rad_sw_cs_hr_av ) ) THEN
               ALLOCATE( rad_sw_cs_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_sw_cs_hr_av = REAL( fill_value, KIND = wp )
            ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = rad_sw_cs_hr_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid = 'zw'

       CASE ( 'rad_sw_hr_xy', 'rad_sw_hr_xz', 'rad_sw_hr_yz' )
          IF ( av == 0 ) THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = rad_sw_hr(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
            IF ( .NOT. ALLOCATED( rad_sw_hr_av ) ) THEN
               ALLOCATE( rad_sw_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_sw_hr_av = REAL( fill_value, KIND = wp )
            ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = rad_sw_hr_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid = 'zw'

       CASE DEFAULT
          found = .FALSE.
          grid  = 'none'

    END SELECT

 END SUBROUTINE radiation_data_output_2d


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining 3D output variables
!------------------------------------------------------------------------------!
 SUBROUTINE radiation_data_output_3d( av, variable, found, local_pf, nzb_do, nzt_do )


    USE indices

    USE kinds


    IMPLICIT NONE

    CHARACTER (LEN=*) ::  variable !<

    INTEGER(iwp) ::  av          !<
    INTEGER(iwp) ::  i, j, k, l  !<
    INTEGER(iwp) ::  nzb_do      !<
    INTEGER(iwp) ::  nzt_do      !<

    LOGICAL      ::  found       !<

    REAL(wp)     ::  fill_value = -999.0_wp    !< value for the _FillValue attribute

    REAL(sp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf !<

    CHARACTER (len=varnamelength)                   :: var, surfid
    INTEGER(iwp)                                    :: ids,idsint_u,idsint_l,isurf,isvf,isurfs,isurflt,ipcgb
    INTEGER(iwp)                                    :: is, js, ks, istat

    found = .TRUE.
    var = TRIM(variable)

!-- check if variable belongs to radiation related variables (starts with rad or rtm)
    IF ( len(var) < 3_iwp  )  THEN
       found = .FALSE.
       RETURN
    ENDIF

    IF ( var(1:3) /= 'rad'  .AND.  var(1:3) /= 'rtm' )  THEN
       found = .FALSE.
       RETURN
    ENDIF

    ids = -1
    DO i = 0, nd-1
        k = len(TRIM(var))
        j = len(TRIM(dirname(i)))
        IF ( k-j+1 >= 1_iwp ) THEN
           IF ( TRIM(var(k-j+1:k)) == TRIM(dirname(i)) )  THEN
              ids = i
              idsint_u = dirint_u(ids)
              idsint_l = dirint_l(ids)
              var = var(:k-j)
              EXIT
           ENDIF
        ENDIF
    ENDDO
    IF ( ids == -1 )  THEN
        var = TRIM(variable)
    ENDIF

    IF ( (var(1:8) == 'rtm_svf_'  .OR.  var(1:8) == 'rtm_dif_')  .AND.  len(TRIM(var)) >= 13 )  THEN
!--     svf values to particular surface
        surfid = var(9:)
        i = index(surfid,'_')
        j = index(surfid(i+1:),'_')
        READ(surfid(1:i-1),*, iostat=istat ) is
        IF ( istat == 0 )  THEN
            READ(surfid(i+1:i+j-1),*, iostat=istat ) js
        ENDIF
        IF ( istat == 0 )  THEN
            READ(surfid(i+j+1:),*, iostat=istat ) ks
        ENDIF
        IF ( istat == 0 )  THEN
            var = var(1:7)
        ENDIF
    ENDIF

    local_pf = fill_value

    SELECT CASE ( TRIM( var ) )
!--   block of large scale radiation model (e.g. RRTMG) output variables
      CASE ( 'rad_sw_in' )
         IF ( av == 0 )  THEN
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_sw_in(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            IF ( .NOT. ALLOCATED( rad_sw_in_av ) ) THEN
               ALLOCATE( rad_sw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_sw_in_av = REAL( fill_value, KIND = wp )
            ENDIF
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_sw_in_av(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

      CASE ( 'rad_sw_out' )
         IF ( av == 0 )  THEN
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_sw_out(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            IF ( .NOT. ALLOCATED( rad_sw_out_av ) ) THEN
               ALLOCATE( rad_sw_out_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_sw_out_av = REAL( fill_value, KIND = wp )
            ENDIF
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_sw_out_av(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

      CASE ( 'rad_sw_cs_hr' )
         IF ( av == 0 )  THEN
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_sw_cs_hr(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            IF ( .NOT. ALLOCATED( rad_sw_cs_hr_av ) ) THEN
               ALLOCATE( rad_sw_cs_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_sw_cs_hr_av = REAL( fill_value, KIND = wp )
            ENDIF
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_sw_cs_hr_av(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

      CASE ( 'rad_sw_hr' )
         IF ( av == 0 )  THEN
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_sw_hr(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            IF ( .NOT. ALLOCATED( rad_sw_hr_av ) ) THEN
               ALLOCATE( rad_sw_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_sw_hr_av = REAL( fill_value, KIND = wp )
            ENDIF
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_sw_hr_av(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

      CASE ( 'rad_lw_in' )
         IF ( av == 0 )  THEN
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_lw_in(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            IF ( .NOT. ALLOCATED( rad_lw_in_av ) ) THEN
               ALLOCATE( rad_lw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_lw_in_av = REAL( fill_value, KIND = wp )
            ENDIF
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_lw_in_av(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

      CASE ( 'rad_lw_out' )
         IF ( av == 0 )  THEN
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_lw_out(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            IF ( .NOT. ALLOCATED( rad_lw_out_av ) ) THEN
               ALLOCATE( rad_lw_out_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_lw_out_av = REAL( fill_value, KIND = wp )
            ENDIF
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_lw_out_av(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

      CASE ( 'rad_lw_cs_hr' )
         IF ( av == 0 )  THEN
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_lw_cs_hr(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            IF ( .NOT. ALLOCATED( rad_lw_cs_hr_av ) ) THEN
               ALLOCATE( rad_lw_cs_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_lw_cs_hr_av = REAL( fill_value, KIND = wp )
            ENDIF
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_lw_cs_hr_av(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

      CASE ( 'rad_lw_hr' )
         IF ( av == 0 )  THEN
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_lw_hr(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            IF ( .NOT. ALLOCATED( rad_lw_hr_av ) ) THEN
               ALLOCATE( rad_lw_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
              rad_lw_hr_av = REAL( fill_value, KIND = wp )
            ENDIF
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_lw_hr_av(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

      CASE ( 'rtm_rad_net' )
!--     array of complete radiation balance
         DO isurf = dirstart(ids), dirend(ids)
            IF ( surfl(id,isurf) == idsint_u  .OR.  surfl(id,isurf) == idsint_l )  THEN
               IF ( av == 0 )  THEN
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = &
                         surfinsw(isurf) - surfoutsw(isurf) +  surfinlw(isurf) - surfoutlw(isurf)
               ELSE
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfradnet_av(isurf)
               ENDIF
            ENDIF
         ENDDO

      CASE ( 'rtm_rad_insw' )
!--      array of sw radiation falling to surface after i-th reflection
         DO isurf = dirstart(ids), dirend(ids)
            IF ( surfl(id,isurf) == idsint_u  .OR.  surfl(id,isurf) == idsint_l )  THEN
               IF ( av == 0 )  THEN
                 local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfinsw(isurf)
               ELSE
                 local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfinsw_av(isurf)
               ENDIF
            ENDIF
         ENDDO

      CASE ( 'rtm_rad_inlw' )
!--      array of lw radiation falling to surface after i-th reflection
         DO isurf = dirstart(ids), dirend(ids)
            IF ( surfl(id,isurf) == idsint_u  .OR.  surfl(id,isurf) == idsint_l )  THEN
               IF ( av == 0 )  THEN
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfinlw(isurf)
               ELSE
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfinlw_av(isurf)
               ENDIF
             ENDIF
         ENDDO

      CASE ( 'rtm_rad_inswdir' )
!--      array of direct sw radiation falling to surface from sun
         DO isurf = dirstart(ids), dirend(ids)
            IF ( surfl(id,isurf) == idsint_u  .OR.  surfl(id,isurf) == idsint_l )  THEN
               IF ( av == 0 )  THEN
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfinswdir(isurf)
               ELSE
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfinswdir_av(isurf)
               ENDIF
            ENDIF
         ENDDO

      CASE ( 'rtm_rad_inswdif' )
!--      array of difusion sw radiation falling to surface from sky and borders of the domain
         DO isurf = dirstart(ids), dirend(ids)
            IF ( surfl(id,isurf) == idsint_u  .OR.  surfl(id,isurf) == idsint_l )  THEN
               IF ( av == 0 )  THEN
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfinswdif(isurf)
               ELSE
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfinswdif_av(isurf)
               ENDIF
            ENDIF
         ENDDO

      CASE ( 'rtm_rad_inswref' )
!--      array of sw radiation falling to surface from reflections
         DO isurf = dirstart(ids), dirend(ids)
            IF ( surfl(id,isurf) == idsint_u  .OR.  surfl(id,isurf) == idsint_l )  THEN
               IF ( av == 0 )  THEN
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = &
                    surfinsw(isurf) - surfinswdir(isurf) - surfinswdif(isurf)
               ELSE
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfinswref_av(isurf)
               ENDIF
            ENDIF
         ENDDO

      CASE ( 'rtm_rad_inlwdif' )
!--      array of difusion lw radiation falling to surface from sky and borders of the domain
         DO isurf = dirstart(ids), dirend(ids)
            IF ( surfl(id,isurf) == idsint_u  .OR.  surfl(id,isurf) == idsint_l )  THEN
               IF ( av == 0 )  THEN
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfinlwdif(isurf)
               ELSE
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfinlwdif_av(isurf)
               ENDIF
            ENDIF
         ENDDO

      CASE ( 'rtm_rad_inlwref' )
!--      array of lw radiation falling to surface from reflections
         DO isurf = dirstart(ids), dirend(ids)
            IF ( surfl(id,isurf) == idsint_u  .OR.  surfl(id,isurf) == idsint_l )  THEN
               IF ( av == 0 )  THEN
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfinlw(isurf) - surfinlwdif(isurf)
               ELSE
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfinlwref_av(isurf)
               ENDIF
            ENDIF
         ENDDO

      CASE ( 'rtm_rad_outsw' )
!--      array of sw radiation emitted from surface after i-th reflection
         DO isurf = dirstart(ids), dirend(ids)
            IF ( surfl(id,isurf) == idsint_u  .OR.  surfl(id,isurf) == idsint_l )  THEN
               IF ( av == 0 )  THEN
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfoutsw(isurf)
               ELSE
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfoutsw_av(isurf)
               ENDIF
            ENDIF
         ENDDO

      CASE ( 'rtm_rad_outlw' )
!--      array of lw radiation emitted from surface after i-th reflection
         DO isurf = dirstart(ids), dirend(ids)
            IF ( surfl(id,isurf) == idsint_u  .OR.  surfl(id,isurf) == idsint_l )  THEN
               IF ( av == 0 )  THEN
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfoutlw(isurf)
               ELSE
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfoutlw_av(isurf)
               ENDIF
            ENDIF
         ENDDO

      CASE ( 'rtm_rad_ressw' )
!--      average of array of residua of sw radiation absorbed in surface after last reflection
         DO isurf = dirstart(ids), dirend(ids)
            IF ( surfl(id,isurf) == idsint_u  .OR.  surfl(id,isurf) == idsint_l )  THEN
               IF ( av == 0 )  THEN
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfins(isurf)
               ELSE
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfins_av(isurf)
               ENDIF
            ENDIF
         ENDDO

      CASE ( 'rtm_rad_reslw' )
!--      average of array of residua of lw radiation absorbed in surface after last reflection
         DO isurf = dirstart(ids), dirend(ids)
            IF ( surfl(id,isurf) == idsint_u  .OR.  surfl(id,isurf) == idsint_l )  THEN
               IF ( av == 0 )  THEN
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfinl(isurf)
               ELSE
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfinl_av(isurf)
               ENDIF
            ENDIF
         ENDDO

      CASE ( 'rtm_rad_pc_inlw' )
!--      array of lw radiation absorbed by plant canopy
         DO ipcgb = 1, npcbl
            IF ( av == 0 )  THEN
               local_pf(pcbl(ix,ipcgb),pcbl(iy,ipcgb),pcbl(iz,ipcgb)) = pcbinlw(ipcgb)
            ELSE
               local_pf(pcbl(ix,ipcgb),pcbl(iy,ipcgb),pcbl(iz,ipcgb)) = pcbinlw_av(ipcgb)
            ENDIF
         ENDDO

      CASE ( 'rtm_rad_pc_insw' )
!--      array of sw radiation absorbed by plant canopy
         DO ipcgb = 1, npcbl
            IF ( av == 0 )  THEN
              local_pf(pcbl(ix,ipcgb),pcbl(iy,ipcgb),pcbl(iz,ipcgb)) = pcbinsw(ipcgb)
            ELSE
              local_pf(pcbl(ix,ipcgb),pcbl(iy,ipcgb),pcbl(iz,ipcgb)) = pcbinsw_av(ipcgb)
            ENDIF
         ENDDO

      CASE ( 'rtm_rad_pc_inswdir' )
!--      array of direct sw radiation absorbed by plant canopy
         DO ipcgb = 1, npcbl
            IF ( av == 0 )  THEN
               local_pf(pcbl(ix,ipcgb),pcbl(iy,ipcgb),pcbl(iz,ipcgb)) = pcbinswdir(ipcgb)
            ELSE
               local_pf(pcbl(ix,ipcgb),pcbl(iy,ipcgb),pcbl(iz,ipcgb)) = pcbinswdir_av(ipcgb)
            ENDIF
         ENDDO

      CASE ( 'rtm_rad_pc_inswdif' )
!--      array of diffuse sw radiation absorbed by plant canopy
         DO ipcgb = 1, npcbl
            IF ( av == 0 )  THEN
               local_pf(pcbl(ix,ipcgb),pcbl(iy,ipcgb),pcbl(iz,ipcgb)) = pcbinswdif(ipcgb)
            ELSE
               local_pf(pcbl(ix,ipcgb),pcbl(iy,ipcgb),pcbl(iz,ipcgb)) = pcbinswdif_av(ipcgb)
            ENDIF
         ENDDO

      CASE ( 'rtm_rad_pc_inswref' )
!--      array of reflected sw radiation absorbed by plant canopy
         DO ipcgb = 1, npcbl
            IF ( av == 0 )  THEN
               local_pf(pcbl(ix,ipcgb),pcbl(iy,ipcgb),pcbl(iz,ipcgb)) = &
                                    pcbinsw(ipcgb) - pcbinswdir(ipcgb) - pcbinswdif(ipcgb)
            ELSE
               local_pf(pcbl(ix,ipcgb),pcbl(iy,ipcgb),pcbl(iz,ipcgb)) = pcbinswref_av(ipcgb)
            ENDIF
         ENDDO

      CASE ( 'rtm_mrt_sw' )
         local_pf = REAL( fill_value, KIND = wp )
         IF ( av == 0 )  THEN
            DO  l = 1, nmrtbl
               local_pf(mrtbl(ix,l),mrtbl(iy,l),mrtbl(iz,l)) = mrtinsw(l)
            ENDDO
         ELSE
            IF ( ALLOCATED( mrtinsw_av ) ) THEN
               DO  l = 1, nmrtbl
                  local_pf(mrtbl(ix,l),mrtbl(iy,l),mrtbl(iz,l)) = mrtinsw_av(l)
               ENDDO
            ENDIF
         ENDIF

      CASE ( 'rtm_mrt_lw' )
         local_pf = REAL( fill_value, KIND = wp )
         IF ( av == 0 )  THEN
            DO  l = 1, nmrtbl
               local_pf(mrtbl(ix,l),mrtbl(iy,l),mrtbl(iz,l)) = mrtinlw(l)
            ENDDO
         ELSE
            IF ( ALLOCATED( mrtinlw_av ) ) THEN
               DO  l = 1, nmrtbl
                  local_pf(mrtbl(ix,l),mrtbl(iy,l),mrtbl(iz,l)) = mrtinlw_av(l)
               ENDDO
            ENDIF
         ENDIF

      CASE ( 'rtm_mrt' )
         local_pf = REAL( fill_value, KIND = wp )
         IF ( av == 0 )  THEN
            DO  l = 1, nmrtbl
               local_pf(mrtbl(ix,l),mrtbl(iy,l),mrtbl(iz,l)) = mrt(l)
            ENDDO
         ELSE
            IF ( ALLOCATED( mrt_av ) ) THEN
               DO  l = 1, nmrtbl
                  local_pf(mrtbl(ix,l),mrtbl(iy,l),mrtbl(iz,l)) = mrt_av(l)
               ENDDO
            ENDIF
         ENDIF
!
!--   block of RTM output variables
!--   variables are intended mainly for debugging and detailed analyse purposes
      CASE ( 'rtm_skyvf' )
!
!--      sky view factor
         DO isurf = dirstart(ids), dirend(ids)
            IF ( surfl(id,isurf) == idsint_u  .OR.  surfl(id,isurf) == idsint_l )  THEN
               local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = skyvf(isurf)
            ENDIF
         ENDDO

      CASE ( 'rtm_skyvft' )
!
!--      sky view factor
         DO isurf = dirstart(ids), dirend(ids)
            IF ( surfl(id,isurf) == ids )  THEN
               local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = skyvft(isurf)
            ENDIF
         ENDDO

      CASE ( 'rtm_svf', 'rtm_dif' )
!
!--      shape view factors or iradiance factors to selected surface
         IF ( TRIM(var)=='rtm_svf' )  THEN
             k = 1
         ELSE
             k = 2
         ENDIF
         DO isvf = 1, nsvfl
            isurflt = svfsurf(1, isvf)
            isurfs = svfsurf(2, isvf)

            IF ( surf(ix,isurfs) == is  .AND.  surf(iy,isurfs) == js  .AND. surf(iz,isurfs) == ks  .AND. &
                 (surf(id,isurfs) == idsint_u .OR. surfl(id,isurfs) == idsint_l ) ) THEN
!
!--            correct source surface
               local_pf(surfl(ix,isurflt),surfl(iy,isurflt),surfl(iz,isurflt)) = svf(k,isvf)
            ENDIF
         ENDDO

      CASE ( 'rtm_surfalb' )
!
!--      surface albedo
         DO isurf = dirstart(ids), dirend(ids)
            IF ( surfl(id,isurf) == idsint_u  .OR.  surfl(id,isurf) == idsint_l )  THEN
               local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = albedo_surf(isurf)
            ENDIF
         ENDDO

      CASE ( 'rtm_surfemis' )
!
!--      surface emissivity, weighted average
         DO isurf = dirstart(ids), dirend(ids)
            IF ( surfl(id,isurf) == idsint_u  .OR.  surfl(id,isurf) == idsint_l )  THEN
               local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = emiss_surf(isurf)
            ENDIF
         ENDDO

      CASE DEFAULT
         found = .FALSE.

    END SELECT


 END SUBROUTINE radiation_data_output_3d

!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining masked data output
!------------------------------------------------------------------------------!
 SUBROUTINE radiation_data_output_mask( av, variable, found, local_pf, mid )

    USE control_parameters

    USE indices

    USE kinds


    IMPLICIT NONE

    CHARACTER (LEN=*) ::  variable   !<

    CHARACTER(LEN=5) ::  grid        !< flag to distinquish between staggered grids

    INTEGER(iwp) ::  av              !<
    INTEGER(iwp) ::  i               !<
    INTEGER(iwp) ::  j               !<
    INTEGER(iwp) ::  k               !<
    INTEGER(iwp) ::  mid             !< masked output running index
    INTEGER(iwp) ::  topo_top_index  !< k index of highest horizontal surface

    LOGICAL ::  found                !< true if output array was found
    LOGICAL ::  resorted             !< true if array is resorted


    REAL(wp),                                                                  &
       DIMENSION(mask_size_l(mid,1),mask_size_l(mid,2),mask_size_l(mid,3)) ::  &
          local_pf   !<

    REAL(wp), DIMENSION(:,:,:), POINTER ::  to_be_resorted  !< points to array which needs to be resorted for output


    found    = .TRUE.
    grid     = 's'
    resorted = .FALSE.

    SELECT CASE ( TRIM( variable ) )


       CASE ( 'rad_lw_in' )
          IF ( av == 0 )  THEN
             to_be_resorted => rad_lw_in
          ELSE
             to_be_resorted => rad_lw_in_av
          ENDIF

       CASE ( 'rad_lw_out' )
          IF ( av == 0 )  THEN
             to_be_resorted => rad_lw_out
          ELSE
             to_be_resorted => rad_lw_out_av
          ENDIF

       CASE ( 'rad_lw_cs_hr' )
          IF ( av == 0 )  THEN
             to_be_resorted => rad_lw_cs_hr
          ELSE
             to_be_resorted => rad_lw_cs_hr_av
          ENDIF

       CASE ( 'rad_lw_hr' )
          IF ( av == 0 )  THEN
             to_be_resorted => rad_lw_hr
          ELSE
             to_be_resorted => rad_lw_hr_av
          ENDIF

       CASE ( 'rad_sw_in' )
          IF ( av == 0 )  THEN
             to_be_resorted => rad_sw_in
          ELSE
             to_be_resorted => rad_sw_in_av
          ENDIF

       CASE ( 'rad_sw_out' )
          IF ( av == 0 )  THEN
             to_be_resorted => rad_sw_out
          ELSE
             to_be_resorted => rad_sw_out_av
          ENDIF

       CASE ( 'rad_sw_cs_hr' )
          IF ( av == 0 )  THEN
             to_be_resorted => rad_sw_cs_hr
          ELSE
             to_be_resorted => rad_sw_cs_hr_av
          ENDIF

       CASE ( 'rad_sw_hr' )
          IF ( av == 0 )  THEN
             to_be_resorted => rad_sw_hr
          ELSE
             to_be_resorted => rad_sw_hr_av
          ENDIF

       CASE DEFAULT
          found = .FALSE.

    END SELECT

!
!-- Resort the array to be output, if not done above
    IF ( found  .AND.  .NOT. resorted )  THEN
       IF ( .NOT. mask_surface(mid) )  THEN
!
!--       Default masked output
          DO  i = 1, mask_size_l(mid,1)
             DO  j = 1, mask_size_l(mid,2)
                DO  k = 1, mask_size_l(mid,3)
                   local_pf(i,j,k) =  to_be_resorted(mask_k(mid,k), &
                                      mask_j(mid,j),mask_i(mid,i))
                ENDDO
             ENDDO
          ENDDO

       ELSE
!
!--       Terrain-following masked output
          DO  i = 1, mask_size_l(mid,1)
             DO  j = 1, mask_size_l(mid,2)
!
!--             Get k index of highest horizontal surface
                topo_top_index = topo_top_ind(mask_j(mid,j), &
                                              mask_i(mid,i),   &
                                              0 )
!
!--             Save output array
                DO  k = 1, mask_size_l(mid,3)
                   local_pf(i,j,k) = to_be_resorted(                         &
                                          MIN( topo_top_index+mask_k(mid,k), &
                                               nzt+1 ),                      &
                                          mask_j(mid,j),                     &
                                          mask_i(mid,i)                     )
                ENDDO
             ENDDO
          ENDDO

       ENDIF
    ENDIF



 END SUBROUTINE radiation_data_output_mask


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine writes local (subdomain) restart data
!------------------------------------------------------------------------------!
 SUBROUTINE radiation_wrd_local


    IMPLICIT NONE


    IF ( ALLOCATED( rad_net_av ) )  THEN
       CALL wrd_write_string( 'rad_net_av' )
       WRITE ( 14 )  rad_net_av
    ENDIF

    IF ( ALLOCATED( rad_lw_in_xy_av ) )  THEN
       CALL wrd_write_string( 'rad_lw_in_xy_av' )
       WRITE ( 14 )  rad_lw_in_xy_av
    ENDIF

    IF ( ALLOCATED( rad_lw_out_xy_av ) )  THEN
       CALL wrd_write_string( 'rad_lw_out_xy_av' )
       WRITE ( 14 )  rad_lw_out_xy_av
    ENDIF

    IF ( ALLOCATED( rad_sw_in_xy_av ) )  THEN
       CALL wrd_write_string( 'rad_sw_in_xy_av' )
       WRITE ( 14 )  rad_sw_in_xy_av
    ENDIF

    IF ( ALLOCATED( rad_sw_out_xy_av ) )  THEN
       CALL wrd_write_string( 'rad_sw_out_xy_av' )
       WRITE ( 14 )  rad_sw_out_xy_av
    ENDIF

    IF ( ALLOCATED( rad_lw_in ) )  THEN
       CALL wrd_write_string( 'rad_lw_in' )
       WRITE ( 14 )  rad_lw_in
    ENDIF

    IF ( ALLOCATED( rad_lw_in_av ) )  THEN
       CALL wrd_write_string( 'rad_lw_in_av' )
       WRITE ( 14 )  rad_lw_in_av
    ENDIF

    IF ( ALLOCATED( rad_lw_out ) )  THEN
       CALL wrd_write_string( 'rad_lw_out' )
       WRITE ( 14 )  rad_lw_out
    ENDIF

    IF ( ALLOCATED( rad_lw_out_av) )  THEN
       CALL wrd_write_string( 'rad_lw_out_av' )
       WRITE ( 14 )  rad_lw_out_av
    ENDIF

    IF ( ALLOCATED( rad_lw_cs_hr) )  THEN
       CALL wrd_write_string( 'rad_lw_cs_hr' )
       WRITE ( 14 )  rad_lw_cs_hr
    ENDIF

    IF ( ALLOCATED( rad_lw_cs_hr_av) )  THEN
       CALL wrd_write_string( 'rad_lw_cs_hr_av' )
       WRITE ( 14 )  rad_lw_cs_hr_av
    ENDIF

    IF ( ALLOCATED( rad_lw_hr) )  THEN
       CALL wrd_write_string( 'rad_lw_hr' )
       WRITE ( 14 )  rad_lw_hr
    ENDIF

    IF ( ALLOCATED( rad_lw_hr_av) )  THEN
       CALL wrd_write_string( 'rad_lw_hr_av' )
       WRITE ( 14 )  rad_lw_hr_av
    ENDIF

    IF ( ALLOCATED( rad_sw_in) )  THEN
       CALL wrd_write_string( 'rad_sw_in' )
       WRITE ( 14 )  rad_sw_in
    ENDIF

    IF ( ALLOCATED( rad_sw_in_av) )  THEN
       CALL wrd_write_string( 'rad_sw_in_av' )
       WRITE ( 14 )  rad_sw_in_av
    ENDIF

    IF ( ALLOCATED( rad_sw_out) )  THEN
       CALL wrd_write_string( 'rad_sw_out' )
       WRITE ( 14 )  rad_sw_out
    ENDIF

    IF ( ALLOCATED( rad_sw_out_av) )  THEN
       CALL wrd_write_string( 'rad_sw_out_av' )
       WRITE ( 14 )  rad_sw_out_av
    ENDIF

    IF ( ALLOCATED( rad_sw_cs_hr) )  THEN
       CALL wrd_write_string( 'rad_sw_cs_hr' )
       WRITE ( 14 )  rad_sw_cs_hr
    ENDIF

    IF ( ALLOCATED( rad_sw_cs_hr_av) )  THEN
       CALL wrd_write_string( 'rad_sw_cs_hr_av' )
       WRITE ( 14 )  rad_sw_cs_hr_av
    ENDIF

    IF ( ALLOCATED( rad_sw_hr) )  THEN
       CALL wrd_write_string( 'rad_sw_hr' )
       WRITE ( 14 )  rad_sw_hr
    ENDIF

    IF ( ALLOCATED( rad_sw_hr_av) )  THEN
       CALL wrd_write_string( 'rad_sw_hr_av' )
       WRITE ( 14 )  rad_sw_hr_av
    ENDIF


 END SUBROUTINE radiation_wrd_local

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine reads local (subdomain) restart data
!------------------------------------------------------------------------------!
 SUBROUTINE radiation_rrd_local( k, nxlf, nxlc, nxl_on_file, nxrf, nxrc,       &
                                nxr_on_file, nynf, nync, nyn_on_file, nysf,    &
                                nysc, nys_on_file, tmp_2d, tmp_3d, found )


    USE control_parameters

    USE indices

    USE kinds

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

    LOGICAL, INTENT(OUT)  :: found

    REAL(wp), DIMENSION(nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) :: tmp_2d   !<

    REAL(wp), DIMENSION(nzb:nzt+1,nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) :: tmp_3d   !<

    REAL(wp), DIMENSION(0:0,nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) :: tmp_3d2   !<


    found = .TRUE.


    SELECT CASE ( restart_string(1:length) )

       CASE ( 'rad_net_av' )
          IF ( .NOT. ALLOCATED( rad_net_av ) )  THEN
             ALLOCATE( rad_net_av(nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_2d
          rad_net_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  =           &
                        tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'rad_lw_in_xy_av' )
          IF ( .NOT. ALLOCATED( rad_lw_in_xy_av ) )  THEN
             ALLOCATE( rad_lw_in_xy_av(nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_2d
          rad_lw_in_xy_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  =      &
                        tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'rad_lw_out_xy_av' )
          IF ( .NOT. ALLOCATED( rad_lw_out_xy_av ) )  THEN
             ALLOCATE( rad_lw_out_xy_av(nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_2d
          rad_lw_out_xy_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  =     &
                        tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'rad_sw_in_xy_av' )
          IF ( .NOT. ALLOCATED( rad_sw_in_xy_av ) )  THEN
             ALLOCATE( rad_sw_in_xy_av(nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_2d
          rad_sw_in_xy_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  =      &
                        tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'rad_sw_out_xy_av' )
          IF ( .NOT. ALLOCATED( rad_sw_out_xy_av ) )  THEN
             ALLOCATE( rad_sw_out_xy_av(nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_2d
          rad_sw_out_xy_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  =     &
                        tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'rad_lw_in' )
          IF ( .NOT. ALLOCATED( rad_lw_in ) )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.                    &
                  radiation_scheme == 'constant'   .OR.                    &
                  radiation_scheme == 'external' )  THEN
                ALLOCATE( rad_lw_in(0:0,nysg:nyng,nxlg:nxrg) )
             ELSE
                ALLOCATE( rad_lw_in(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
          ENDIF
          IF ( k == 1 )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.                    &
                  radiation_scheme == 'constant'   .OR.                    &
                  radiation_scheme == 'external' )  THEN
                READ ( 13 )  tmp_3d2
                rad_lw_in(0:0,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =   &
                   tmp_3d2(0:0,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ELSE
                READ ( 13 )  tmp_3d
                rad_lw_in(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =     &
                    tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ENDIF
          ENDIF

       CASE ( 'rad_lw_in_av' )
          IF ( .NOT. ALLOCATED( rad_lw_in_av ) )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.                    &
                  radiation_scheme == 'constant'   .OR.                    &
                  radiation_scheme == 'external' )  THEN
                ALLOCATE( rad_lw_in_av(0:0,nysg:nyng,nxlg:nxrg) )
             ELSE
                ALLOCATE( rad_lw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
          ENDIF
          IF ( k == 1 )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.                    &
                  radiation_scheme == 'constant'   .OR.                    &
                  radiation_scheme == 'external' )  THEN
                READ ( 13 )  tmp_3d2
                rad_lw_in_av(0:0,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =&
                    tmp_3d2(0:0,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ELSE
                READ ( 13 )  tmp_3d
                rad_lw_in_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =  &
                    tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ENDIF
          ENDIF

       CASE ( 'rad_lw_out' )
          IF ( .NOT. ALLOCATED( rad_lw_out ) )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.                    &
                  radiation_scheme == 'constant'   .OR.                    &
                  radiation_scheme == 'external' )  THEN
                ALLOCATE( rad_lw_out(0:0,nysg:nyng,nxlg:nxrg) )
             ELSE
                ALLOCATE( rad_lw_out(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
          ENDIF
          IF ( k == 1 )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.                    &
                  radiation_scheme == 'constant'   .OR.                    &
                  radiation_scheme == 'external' )  THEN
                READ ( 13 )  tmp_3d2
                rad_lw_out(0:0,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =  &
                    tmp_3d2(0:0,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ELSE
                READ ( 13 )  tmp_3d
                rad_lw_out(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =    &
                    tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ENDIF
          ENDIF

       CASE ( 'rad_lw_out_av' )
          IF ( .NOT. ALLOCATED( rad_lw_out_av ) )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.                    &
                  radiation_scheme == 'constant'   .OR.                    &
                  radiation_scheme == 'external' )  THEN
                ALLOCATE( rad_lw_out_av(0:0,nysg:nyng,nxlg:nxrg) )
             ELSE
                ALLOCATE( rad_lw_out_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
          ENDIF
          IF ( k == 1 )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.                    &
                  radiation_scheme == 'constant'   .OR.                    &
                  radiation_scheme == 'external' )  THEN
                READ ( 13 )  tmp_3d2
                rad_lw_out_av(0:0,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) &
                   = tmp_3d2(0:0,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ELSE
                READ ( 13 )  tmp_3d
                rad_lw_out_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                    tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ENDIF
          ENDIF

       CASE ( 'rad_lw_cs_hr' )
          IF ( .NOT. ALLOCATED( rad_lw_cs_hr ) )  THEN
             ALLOCATE( rad_lw_cs_hr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
          rad_lw_cs_hr(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =        &
                  tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'rad_lw_cs_hr_av' )
          IF ( .NOT. ALLOCATED( rad_lw_cs_hr_av ) )  THEN
             ALLOCATE( rad_lw_cs_hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
          rad_lw_cs_hr_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =     &
                  tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'rad_lw_hr' )
          IF ( .NOT. ALLOCATED( rad_lw_hr ) )  THEN
             ALLOCATE( rad_lw_hr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
          rad_lw_hr(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =           &
                  tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'rad_lw_hr_av' )
          IF ( .NOT. ALLOCATED( rad_lw_hr_av ) )  THEN
             ALLOCATE( rad_lw_hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
          rad_lw_hr_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =        &
                  tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'rad_sw_in' )
          IF ( .NOT. ALLOCATED( rad_sw_in ) )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.                    &
                  radiation_scheme == 'constant'   .OR.                    &
                  radiation_scheme == 'external' )  THEN
                ALLOCATE( rad_sw_in(0:0,nysg:nyng,nxlg:nxrg) )
             ELSE
                ALLOCATE( rad_sw_in(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
          ENDIF
          IF ( k == 1 )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.                    &
                  radiation_scheme == 'constant'   .OR.                    &
                  radiation_scheme == 'external' )  THEN
                READ ( 13 )  tmp_3d2
                rad_sw_in(0:0,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =   &
                    tmp_3d2(0:0,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ELSE
                READ ( 13 )  tmp_3d
                rad_sw_in(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =     &
                    tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ENDIF
          ENDIF

       CASE ( 'rad_sw_in_av' )
          IF ( .NOT. ALLOCATED( rad_sw_in_av ) )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.                    &
                  radiation_scheme == 'constant'   .OR.                    &
                  radiation_scheme == 'external' )  THEN
                ALLOCATE( rad_sw_in_av(0:0,nysg:nyng,nxlg:nxrg) )
             ELSE
                ALLOCATE( rad_sw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
          ENDIF
          IF ( k == 1 )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.                    &
                  radiation_scheme == 'constant'   .OR.                    &
                  radiation_scheme == 'external' )  THEN
                READ ( 13 )  tmp_3d2
                rad_sw_in_av(0:0,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =&
                    tmp_3d2(0:0,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ELSE
                READ ( 13 )  tmp_3d
                rad_sw_in_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =  &
                    tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ENDIF
          ENDIF

       CASE ( 'rad_sw_out' )
          IF ( .NOT. ALLOCATED( rad_sw_out ) )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.                    &
                  radiation_scheme == 'constant'   .OR.                    &
                  radiation_scheme == 'external' )  THEN
                ALLOCATE( rad_sw_out(0:0,nysg:nyng,nxlg:nxrg) )
             ELSE
                ALLOCATE( rad_sw_out(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
          ENDIF
          IF ( k == 1 )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.                    &
                  radiation_scheme == 'constant'   .OR.                    &
                  radiation_scheme == 'external' )  THEN
                READ ( 13 )  tmp_3d2
                rad_sw_out(0:0,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =  &
                    tmp_3d2(0:0,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ELSE
                READ ( 13 )  tmp_3d
                rad_sw_out(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =    &
                    tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ENDIF
          ENDIF

       CASE ( 'rad_sw_out_av' )
          IF ( .NOT. ALLOCATED( rad_sw_out_av ) )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.                    &
                  radiation_scheme == 'constant'   .OR.                    &
                  radiation_scheme == 'external' )  THEN
                ALLOCATE( rad_sw_out_av(0:0,nysg:nyng,nxlg:nxrg) )
             ELSE
                ALLOCATE( rad_sw_out_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
          ENDIF
          IF ( k == 1 )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.                    &
                  radiation_scheme == 'constant'   .OR.                    &
                  radiation_scheme == 'external' )  THEN
                READ ( 13 )  tmp_3d2
                rad_sw_out_av(0:0,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) &
                   = tmp_3d2(0:0,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ELSE
                READ ( 13 )  tmp_3d
                rad_sw_out_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                    tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ENDIF
          ENDIF

       CASE ( 'rad_sw_cs_hr' )
          IF ( .NOT. ALLOCATED( rad_sw_cs_hr ) )  THEN
             ALLOCATE( rad_sw_cs_hr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
          rad_sw_cs_hr(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =        &
                  tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'rad_sw_cs_hr_av' )
          IF ( .NOT. ALLOCATED( rad_sw_cs_hr_av ) )  THEN
             ALLOCATE( rad_sw_cs_hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
          rad_sw_cs_hr_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =     &
                  tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'rad_sw_hr' )
          IF ( .NOT. ALLOCATED( rad_sw_hr ) )  THEN
             ALLOCATE( rad_sw_hr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
          rad_sw_hr(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =           &
                  tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'rad_sw_hr_av' )
          IF ( .NOT. ALLOCATED( rad_sw_hr_av ) )  THEN
             ALLOCATE( rad_sw_hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
          rad_lw_hr_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =        &
                  tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE DEFAULT

          found = .FALSE.

    END SELECT

 END SUBROUTINE radiation_rrd_local


 END MODULE radiation_model_mod
