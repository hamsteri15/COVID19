!> @file chemistry_model_mod.f90
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
! Copyright 2017-2019 Leibniz Universitaet Hannover
! Copyright 2017-2019 Karlsruhe Institute of Technology
! Copyright 2017-2019 Freie Universitaet Berlin
!------------------------------------------------------------------------------!
!
! Current revisions:
! -----------------
! 
! 
! Former revisions:
! -----------------
! $Id: chemistry_model_mod.f90 4457 2020-03-11 14:20:43Z raasch $
! use statement for exchange horiz added,
! bugfix for call of exchange horiz 2d
! 
! 4442 2020-03-04 19:21:13Z suehring
! Change order of dimension in surface array %frac to allow for better 
! vectorization.
!
! 4441 2020-03-04 19:20:35Z suehring
! in subroutine chem_init (ECC)
! - allows different init paths emission data for legacy
!   mode emission and on-demand mode
! in subroutine chem_init_internal (ECC)
! - reads netcdf file only when legacy mode is activated
!   (i.e., emiss_read_legacy_mode = .TRUE.)
!   otherwise file is read once at the beginning to obtain
!   header information, and emission data are extracted on
!   an on-demand basis
!
! 4372 2020-01-14 10:20:35Z banzhafs
! chem_parin : added handler for new namelist item emiss_legacy_read_mode (ECC)
! added messages
! CM0465 - legacy read mode selection
! CM0466 - legacy read mode force selection
! CM0467 - new read mode selection
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
! 4306 2019-11-25 12:04:48Z banzhafs
! Corretion for r4304 commit
!
! 4304 2019-11-25 10:43:03Z banzhafs
! Precision clean-up in drydepo_aero_zhang_vd subroutine
!
! 4292 2019-11-11 13:04:50Z banzhafs
! Bugfix for r4290
!
! 4290 2019-11-11 12:06:14Z banzhafs
! Bugfix in sedimentation resistance calculation in drydepo_aero_zhang_vd subroutine
!
! 4273 2019-10-24 13:40:54Z monakurppa
! Add logical switches nesting_chem and nesting_offline_chem (both .TRUE.
! by default)
! 
! 4272 2019-10-23 15:18:57Z schwenkel
! Further modularization of boundary conditions: moved boundary conditions to
! respective modules
!
! 4268 2019-10-17 11:29:38Z schwenkel
! Moving module specific boundary conditions from time_integration to module
! 
! 4230 2019-09-11 13:58:14Z suehring
! Bugfix, initialize mean profiles also in restart runs. Also initialize 
! array used for Runge-Kutta tendecies in restart runs.  
! 
! 4227 2019-09-10 18:04:34Z gronemeier
! implement new palm_date_time_mod
! 
! 4182 2019-08-22 15:20:23Z scharf
! Corrected "Former revisions" section
! 
! 4167 2019-08-16 11:01:48Z suehring
! Changed behaviour of masked output over surface to follow terrain and ignore 
! buildings (J.Resler, T.Gronemeier)
! 
! 
! 4166 2019-08-16 07:54:21Z resler
! Bugfix in decycling
! 
! 4115 2019-07-24 12:50:49Z suehring
! Fix faulty IF statement in decycling initialization
! 
! 4110 2019-07-22 17:05:21Z suehring
! - Decycling boundary conditions are only set at the ghost points not on the
!   prognostic grid points 
! - Allocation and initialization of special advection flags cs_advc_flags_s
!   used for chemistry. These are exclusively used for chemical species in 
!   order to distinguish from the usually-used flags which might be different
!   when decycling is applied in combination with cyclic boundary conditions. 
!   Moreover, cs_advc_flags_s considers extended zones around buildings where 
!   first-order upwind scheme is applied for the horizontal advection terms, 
!   in order to overcome high concentration peaks due to stationary numerical 
!   oscillations caused by horizontal advection discretization.
! 
! 4109 2019-07-22 17:00:34Z suehring
! Slightly revise setting of boundary conditions at horizontal walls, use 
! data-structure offset index instead of pre-calculate it for each facing
! 
! 4080 2019-07-09 18:17:37Z suehring
! Restore accidantly removed limitation to positive values
! 
! 4079 2019-07-09 18:04:41Z suehring
! Application of monotonic flux limiter for the vertical scalar advection 
! up to the topography top (only for the cache-optimized version at the 
! moment).
! 
! 4069 2019-07-01 14:05:51Z Giersch
! Masked output running index mid has been introduced as a local variable to 
! avoid runtime error (Loop variable has been modified) in time_integration
! 
! 4029 2019-06-14 14:04:35Z raasch
! nest_chemistry option removed
! 
! 4004 2019-05-24 11:32:38Z suehring
! in subroutine chem_parin check emiss_lod / mod_emis only
! when emissions_anthropogenic is activated in namelist (E.C. Chan)
! 
! 3968 2019-05-13 11:04:01Z suehring
! - added "emiss_lod" which serves the same function as "mode_emis"
!   both will be synchronized with emiss_lod having pirority over
!   mode_emis (see informational messages)
! - modified existing error message and introduced new informational messages
!    - CM0436 - now also applies to invalid LOD settings
!    - CM0463 - emiss_lod take precedence in case of conflict with mod_emis
!    - CM0464 - emiss_lod valued assigned based on mode_emis if undefined
! 
! 3930 2019-04-24 14:57:18Z forkel
! Changed chem_non_transport_physics to chem_non_advective_processes
! 
! 
! 3929 2019-04-24 12:52:08Z banzhafs
! Correct/complete module_interface introduction for chemistry model
! Add subroutine chem_exchange_horiz_bounds
! Bug fix deposition
!
! 3784 2019-03-05 14:16:20Z banzhafs
! 2D output of emission fluxes
! 
! 3784 2019-03-05 14:16:20Z banzhafs
! Bugfix, uncomment erroneous commented variable used for dry deposition. 
! Bugfix in 3D emission output. 
! 
! 3784 2019-03-05 14:16:20Z banzhafs
! Changes related to global restructuring of location messages and introduction 
! of additional debug messages
! 
! 3784 2019-03-05 14:16:20Z banzhafs
! some formatting of the deposition code
!
! 3784 2019-03-05 14:16:20Z banzhafs
! some formatting
! 
! 3784 2019-03-05 14:16:20Z banzhafs
! added cs_mech to USE chem_gasphase_mod 
! 
! 3784 2019-03-05 14:16:20Z banzhafs
! renamed get_mechanismname to get_mechanism_name
! renamed do_emiss to emissions_anthropogenic and do_depo to deposition_dry (ecc)
! 
! 3784 2019-03-05 14:16:20Z banzhafs
! Unused variables removed/taken care of
! 
!
! 3784 2019-03-05 14:16:20Z forkel
! Replaced READ from unit 10 by CALL get_mechanismname also in chem_header
!
!
! 3783 2019-03-05 13:23:50Z forkel
! Removed forgotte write statements an some unused variables (did not touch the 
! parts related to deposition)
! 
! 
! 3780 2019-03-05 11:19:45Z forkel
! Removed READ from unit 10, added CALL get_mechanismname
! 
! 
! 3767 2019-02-27 08:18:02Z raasch
! unused variable for file index removed from rrd-subroutines parameter list
! 
! 3738 2019-02-12 17:00:45Z suehring
! Clean-up debug prints
! 
! 3737 2019-02-12 16:57:06Z suehring
! Enable mesoscale offline nesting for chemistry variables as well as 
! initialization of chemistry via dynamic input file.
! 
! 3719 2019-02-06 13:10:18Z kanani
! Resolved cpu logpoint overlap with all progn.equations, moved cpu_log call
! to prognostic_equations for better overview
! 
! 3700 2019-01-26 17:03:42Z knoop
! Some interface calls moved to module_interface + cleanup
! 
! 3664 2019-01-09 14:00:37Z forkel
! Replaced misplaced location message by @todo
! 
! 
! 3654 2019-01-07 16:31:57Z suehring
! Disable misplaced location message
! 
! 3652 2019-01-07 15:29:59Z forkel
! Checks added for chemistry mechanism, parameter chem_mechanism added (basit)
! 
! 2718 2018-01-02 08:49:38Z maronga
! Initial revision
!
! 
!
!
! Authors:
! --------
! @author Renate Forkel
! @author Farah Kanani-Suehring
! @author Klaus Ketelsen
! @author Basit Khan
! @author Sabine Banzhaf
!
!
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Chemistry model for PALM-4U
!> @todo Extend chem_species type by nspec and nvar as addititional elements (RF)
!> @todo Check possibility to reduce dimension of chem_species%conc from nspec to nvar (RF)
!> @todo Adjust chem_rrd_local to CASE structure of others modules. It is not 
!>       allowed to use the chemistry model in a precursor run and additionally 
!>       not using it in a main run
!> @todo Implement turbulent inflow of chem spcs in inflow_turbulence.
!> @todo Separate boundary conditions for each chem spcs to be implemented
!> @todo Currently only total concentration are calculated. Resolved, parameterized
!>       and chemistry fluxes although partially and some completely coded but
!>       are not operational/activated in this version. bK.
!> @todo slight differences in passive scalar and chem spcs when chem reactions
!>       turned off. Need to be fixed. bK
!> @todo test nesting for chem spcs, was implemented by suehring (kanani)
!> @todo chemistry error messages
! 
!------------------------------------------------------------------------------!

 MODULE chemistry_model_mod

    USE advec_s_pw_mod,                                                                            &
         ONLY:  advec_s_pw

    USE advec_s_up_mod,                                                                            &
         ONLY:  advec_s_up

    USE advec_ws,                                                                                  &
         ONLY:  advec_s_ws, ws_init_flags_scalar

    USE diffusion_s_mod,                                                                           &
         ONLY:  diffusion_s

    USE kinds,                                                                                     &
         ONLY:  iwp, wp

    USE indices,                                                                                   &
         ONLY:  advc_flags_s,                                                                      &
                nbgp, nx, nxl, nxlg, nxr, nxrg, ny, nyn, nyng, nys, nysg, nz, nzb, nzt,            &
                wall_flags_total_0

    USE pegrid,                                                                                    &
         ONLY: myid, threads_per_task

    USE bulk_cloud_model_mod,                                                                      &
         ONLY:  bulk_cloud_model

    USE control_parameters,                                                                        &
         ONLY:  bc_lr_cyc, bc_ns_cyc,                                                              &
                bc_dirichlet_l,                                                                    &
                bc_dirichlet_n,                                                                    &
                bc_dirichlet_r,                                                                    &
                bc_dirichlet_s,                                                                    &
                bc_radiation_l,                                                                    &
                bc_radiation_n,                                                                    &
                bc_radiation_r,                                                                    &
                bc_radiation_s,                                                                    &
                debug_output,                                                                      &
                dt_3d, humidity, initializing_actions, message_string,                             &
                omega, tsc, intermediate_timestep_count, intermediate_timestep_count_max,          &
                max_pr_user,                                                                       &
                monotonic_limiter_z,                                                               &
                scalar_advec,                                                                      &
                timestep_scheme, use_prescribed_profile_data, ws_scheme_sca, air_chemistry

    USE arrays_3d,                                                                                 &
         ONLY:  exner, hyp, pt, q, ql, rdf_sc, tend, zu

    USE chem_gasphase_mod,                                                                         &
         ONLY:  atol, chem_gasphase_integrate, cs_mech, get_mechanism_name, nkppctrl,              &
         nmaxfixsteps, nphot, nreact, nspec, nvar, phot_names, rtol, spc_names, t_steps, vl_dim

    USE chem_modules

    USE chem_photolysis_mod,                                                                       &
        ONLY:  photolysis_control

    USE cpulog,                                                                                    &
        ONLY:  cpu_log, log_point_s

    USE statistics

    USE surface_mod,                                                                               &
         ONLY:  surf_def_h, surf_def_v, surf_lsm_h, surf_lsm_v, surf_usm_h, surf_usm_v

    IMPLICIT NONE

    PRIVATE
    SAVE

    REAL(kind=wp), ALLOCATABLE, DIMENSION(:,:,:,:), TARGET ::  spec_conc_1  !< pointer for swapping of timelevels for conc
    REAL(kind=wp), ALLOCATABLE, DIMENSION(:,:,:,:), TARGET ::  spec_conc_2  !< pointer for swapping of timelevels for conc
    REAL(kind=wp), ALLOCATABLE, DIMENSION(:,:,:,:), TARGET ::  spec_conc_3  !< pointer for swapping of timelevels for conc
    REAL(kind=wp), ALLOCATABLE, DIMENSION(:,:,:,:), TARGET ::  spec_conc_av !< averaged concentrations of chemical species       
    REAL(kind=wp), ALLOCATABLE, DIMENSION(:,:,:,:), TARGET ::  freq_1       !< pointer for phtolysis frequncies 
                                                                            !< (only 1 timelevel required)
    INTEGER, DIMENSION(nkppctrl)                           ::  icntrl       !< 20 integer parameters for fine tuning KPP code
                                                                            !< (e.g. solver type)
    REAL(kind=wp), DIMENSION(nkppctrl)                     ::  rcntrl       !< 20 real parameters for fine tuning of KPP code
                                                                            !< (e.g starting internal timestep of solver)
!
!-- Decycling of chemistry variables: Dirichlet BCs with cyclic is frequently not
!-- approproate for chemicals compounds since they may accumulate too much. 
!-- If no proper boundary conditions from a DYNAMIC input file are available,
!-- de-cycling applies the initial profiles at the inflow boundaries for 
!-- Dirichlet boundary conditions
    LOGICAL ::  decycle_chem_lr    = .FALSE.    !< switch for setting decycling in left-right direction
    LOGICAL ::  decycle_chem_ns    = .FALSE.    !< switch for setting decycling in south-norht direction
    CHARACTER (LEN=20), DIMENSION(4) ::  decycle_method = &
         (/'dirichlet','dirichlet','dirichlet','dirichlet'/)
                              !< Decycling method at horizontal boundaries, 
                              !< 1=left, 2=right, 3=south, 4=north
                              !< dirichlet = initial size distribution and 
                              !< chemical composition set for the ghost and 
                              !< first three layers
                              !< neumann = zero gradient

    REAL(kind=wp), PUBLIC ::  cs_time_step = 0._wp

!
!-- Parameter needed for Deposition calculation using DEPAC model (van Zanten et al., 2010)
    !
    INTEGER(iwp), PARAMETER ::  nlu_dep = 15                   !< Number of DEPAC landuse classes (lu's)
    INTEGER(iwp), PARAMETER ::  ncmp = 10                      !< Number of DEPAC gas components
    INTEGER(iwp), PARAMETER ::  nposp = 69                     !< Number of possible species for deposition
!
!-- DEPAC landuse classes as defined in LOTOS-EUROS model v2.1                              
    INTEGER(iwp) ::  ilu_grass              = 1        
    INTEGER(iwp) ::  ilu_arable             = 2        
    INTEGER(iwp) ::  ilu_permanent_crops    = 3         
    INTEGER(iwp) ::  ilu_coniferous_forest  = 4         
    INTEGER(iwp) ::  ilu_deciduous_forest   = 5         
    INTEGER(iwp) ::  ilu_water_sea          = 6        
    INTEGER(iwp) ::  ilu_urban              = 7        
    INTEGER(iwp) ::  ilu_other              = 8  
    INTEGER(iwp) ::  ilu_desert             = 9  
    INTEGER(iwp) ::  ilu_ice                = 10 
    INTEGER(iwp) ::  ilu_savanna            = 11 
    INTEGER(iwp) ::  ilu_tropical_forest    = 12 
    INTEGER(iwp) ::  ilu_water_inland       = 13 
    INTEGER(iwp) ::  ilu_mediterrean_scrub  = 14 
    INTEGER(iwp) ::  ilu_semi_natural_veg   = 15 

!
!-- NH3/SO2 ratio regimes:
    INTEGER(iwp), PARAMETER ::  iratns_low      = 1       !< low ratio NH3/SO2
    INTEGER(iwp), PARAMETER ::  iratns_high     = 2       !< high ratio NH3/SO2
    INTEGER(iwp), PARAMETER ::  iratns_very_low = 3       !< very low ratio NH3/SO2
!
!-- Default:
    INTEGER, PARAMETER ::  iratns_default = iratns_low
!
!-- Set alpha for f_light (4.57 is conversion factor from 1./(mumol m-2 s-1) to W m-2
    REAL(wp), DIMENSION(nlu_dep), PARAMETER ::  alpha   =(/ 0.009, 0.009, 0.009, 0.006, 0.006, -999., -999., 0.009, -999.,      &
         -999., 0.009, 0.006, -999., 0.009, 0.008/)*4.57
!
!-- Set temperatures per land use for f_temp
    REAL(wp), DIMENSION(nlu_dep), PARAMETER ::  tmin = (/ 12.0, 12.0,  12.0,  0.0,  0.0, -999., -999., 12.0, -999., -999.,      &
         12.0,  0.0, -999., 12.0,  8.0/)
    REAL(wp), DIMENSION(nlu_dep), PARAMETER ::  topt = (/ 26.0, 26.0,  26.0, 18.0, 20.0, -999., -999., 26.0, -999., -999.,      &
         26.0, 20.0, -999., 26.0, 24.0 /)
    REAL(wp), DIMENSION(nlu_dep), PARAMETER ::  tmax = (/ 40.0, 40.0,  40.0, 36.0, 35.0, -999., -999., 40.0, -999., -999.,      &
         40.0, 35.0, -999., 40.0, 39.0 /)
!
!-- Set f_min:
    REAL(wp), DIMENSION(nlu_dep), PARAMETER ::  f_min = (/ 0.01, 0.01, 0.01, 0.1, 0.1, -999., -999., 0.01, -999., -999., 0.01,  &
         0.1, -999., 0.01, 0.04/)

!
!-- Set maximal conductance (m/s)
!-- (R T/P) = 1/41000 mmol/m3 is given for 20 deg C to go from  mmol O3/m2/s to m/s
    REAL(wp), DIMENSION(nlu_dep), PARAMETER ::  g_max = (/ 270., 300., 300., 140., 150., -999., -999., 270., -999., -999.,      &
         270., 150., -999., 300., 422./)/41000.
!
!-- Set max, min for vapour pressure deficit vpd
    REAL(wp), DIMENSION(nlu_dep), PARAMETER ::  vpd_max = (/1.3, 0.9, 0.9, 0.5, 1.0, -999., -999., 1.3, -999., -999., 1.3,      &
         1.0, -999., 0.9, 2.8/) 
    REAL(wp), DIMENSION(nlu_dep), PARAMETER ::  vpd_min = (/3.0, 2.8, 2.8, 3.0, 3.25, -999., -999., 3.0, -999., -999., 3.0,     &
         3.25, -999., 2.8, 4.5/) 

    PUBLIC nreact
    PUBLIC nspec               !< number of gas phase chemical species including constant compound (e.g. N2)
    PUBLIC nvar                !< number of variable gas phase chemical species (nvar <= nspec)
    PUBLIC spc_names           !< names of gas phase chemical species (come from KPP) (come from KPP)
    PUBLIC spec_conc_2  
!    
!-- Interface section
    INTERFACE chem_3d_data_averaging
       MODULE PROCEDURE chem_3d_data_averaging 
    END INTERFACE chem_3d_data_averaging

    INTERFACE chem_boundary_conds
       MODULE PROCEDURE chem_boundary_conds
       MODULE PROCEDURE chem_boundary_conds_decycle
    END INTERFACE chem_boundary_conds

    INTERFACE chem_boundary_conditions
       MODULE PROCEDURE chem_boundary_conditions
    END INTERFACE chem_boundary_conditions

    INTERFACE chem_check_data_output
       MODULE PROCEDURE chem_check_data_output
    END INTERFACE chem_check_data_output

    INTERFACE chem_data_output_2d
       MODULE PROCEDURE chem_data_output_2d
    END INTERFACE chem_data_output_2d

    INTERFACE chem_data_output_3d
       MODULE PROCEDURE chem_data_output_3d
    END INTERFACE chem_data_output_3d

    INTERFACE chem_data_output_mask
       MODULE PROCEDURE chem_data_output_mask
    END INTERFACE chem_data_output_mask

    INTERFACE chem_check_data_output_pr
       MODULE PROCEDURE chem_check_data_output_pr
    END INTERFACE chem_check_data_output_pr

    INTERFACE chem_check_parameters
       MODULE PROCEDURE chem_check_parameters
    END INTERFACE chem_check_parameters

    INTERFACE chem_define_netcdf_grid
       MODULE PROCEDURE chem_define_netcdf_grid
    END INTERFACE chem_define_netcdf_grid

    INTERFACE chem_header
       MODULE PROCEDURE chem_header
    END INTERFACE chem_header

    INTERFACE chem_init_arrays
       MODULE PROCEDURE chem_init_arrays
    END INTERFACE chem_init_arrays

    INTERFACE chem_init
       MODULE PROCEDURE chem_init
    END INTERFACE chem_init

    INTERFACE chem_init_profiles
       MODULE PROCEDURE chem_init_profiles
    END INTERFACE chem_init_profiles

    INTERFACE chem_integrate
       MODULE PROCEDURE chem_integrate_ij
    END INTERFACE chem_integrate

    INTERFACE chem_parin
       MODULE PROCEDURE chem_parin
    END INTERFACE chem_parin

    INTERFACE chem_actions
       MODULE PROCEDURE chem_actions
       MODULE PROCEDURE chem_actions_ij
    END INTERFACE chem_actions

    INTERFACE chem_non_advective_processes
       MODULE PROCEDURE chem_non_advective_processes
       MODULE PROCEDURE chem_non_advective_processes_ij
    END INTERFACE chem_non_advective_processes
    
    INTERFACE chem_exchange_horiz_bounds
       MODULE PROCEDURE chem_exchange_horiz_bounds
    END INTERFACE chem_exchange_horiz_bounds    

    INTERFACE chem_prognostic_equations
       MODULE PROCEDURE chem_prognostic_equations
       MODULE PROCEDURE chem_prognostic_equations_ij
    END INTERFACE chem_prognostic_equations

    INTERFACE chem_rrd_local
       MODULE PROCEDURE chem_rrd_local
    END INTERFACE chem_rrd_local

    INTERFACE chem_statistics
       MODULE PROCEDURE chem_statistics
    END INTERFACE chem_statistics

    INTERFACE chem_swap_timelevel
       MODULE PROCEDURE chem_swap_timelevel
    END INTERFACE chem_swap_timelevel

    INTERFACE chem_wrd_local
       MODULE PROCEDURE chem_wrd_local 
    END INTERFACE chem_wrd_local

    INTERFACE chem_depo
       MODULE PROCEDURE chem_depo 
    END INTERFACE chem_depo

    INTERFACE drydepos_gas_depac
       MODULE PROCEDURE drydepos_gas_depac 
    END INTERFACE drydepos_gas_depac

    INTERFACE rc_special
       MODULE PROCEDURE rc_special 
    END INTERFACE rc_special

    INTERFACE  rc_gw
       MODULE PROCEDURE rc_gw 
    END INTERFACE rc_gw

    INTERFACE rw_so2 
       MODULE PROCEDURE rw_so2  
    END INTERFACE rw_so2

    INTERFACE rw_nh3_sutton
       MODULE PROCEDURE rw_nh3_sutton 
    END INTERFACE rw_nh3_sutton

    INTERFACE rw_constant
       MODULE PROCEDURE rw_constant 
    END INTERFACE rw_constant

    INTERFACE rc_gstom
       MODULE PROCEDURE rc_gstom 
    END INTERFACE rc_gstom

    INTERFACE rc_gstom_emb
       MODULE PROCEDURE rc_gstom_emb 
    END INTERFACE rc_gstom_emb

    INTERFACE par_dir_diff
       MODULE PROCEDURE par_dir_diff 
    END INTERFACE par_dir_diff

    INTERFACE rc_get_vpd
       MODULE PROCEDURE rc_get_vpd 
    END INTERFACE rc_get_vpd

    INTERFACE rc_gsoil_eff
       MODULE PROCEDURE rc_gsoil_eff 
    END INTERFACE rc_gsoil_eff

    INTERFACE rc_rinc
       MODULE PROCEDURE rc_rinc 
    END INTERFACE rc_rinc

    INTERFACE rc_rctot
       MODULE PROCEDURE rc_rctot  
    END INTERFACE rc_rctot

!    INTERFACE rc_comp_point_rc_eff
!       MODULE PROCEDURE rc_comp_point_rc_eff 
!    END INTERFACE rc_comp_point_rc_eff

    INTERFACE drydepo_aero_zhang_vd
       MODULE PROCEDURE drydepo_aero_zhang_vd 
    END INTERFACE drydepo_aero_zhang_vd

    INTERFACE get_rb_cell
       MODULE PROCEDURE  get_rb_cell
    END INTERFACE get_rb_cell



    PUBLIC chem_3d_data_averaging, chem_boundary_conds, chem_boundary_conditions, &
            chem_boundary_conds_decycle, chem_check_data_output,              &
         chem_check_data_output_pr, chem_check_parameters,                    &
         chem_data_output_2d, chem_data_output_3d, chem_data_output_mask,     &
         chem_define_netcdf_grid, chem_header, chem_init, chem_init_arrays,   &
         chem_init_profiles, chem_integrate, chem_parin,                      &
         chem_actions, chem_prognostic_equations, chem_rrd_local,             &
         chem_statistics, chem_swap_timelevel, chem_wrd_local, chem_depo,     &
         chem_non_advective_processes, chem_exchange_horiz_bounds

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine for averaging 3D data of chemical species. Due to the fact that 
!> the averaged chem arrays are allocated in chem_init, no if-query concerning
!> the allocation is required (in any mode). Attention: If you just specify an 
!> averaged output quantity in the _p3dr file during restarts the first output
!> includes the time between the beginning of the restart run and the first
!> output time (not necessarily the whole averaging_interval you have 
!> specified in your _p3d/_p3dr file )
!------------------------------------------------------------------------------!
 SUBROUTINE chem_3d_data_averaging( mode, variable )

    USE control_parameters

    USE exchange_horiz_mod,                                                    &
        ONLY:  exchange_horiz_2d


    CHARACTER (LEN=*) ::  mode     !< 
    CHARACTER (LEN=*) ::  variable !< 

    LOGICAL ::  match_def  !< flag indicating default-type surface
    LOGICAL ::  match_lsm  !< flag indicating natural-type surface
    LOGICAL ::  match_usm  !< flag indicating urban-type surface

    INTEGER(iwp) ::  i                  !< grid index x direction
    INTEGER(iwp) ::  j                  !< grid index y direction
    INTEGER(iwp) ::  k                  !< grid index z direction
    INTEGER(iwp) ::  m                  !< running index surface type
    INTEGER(iwp) ::  lsp               !< running index for chem spcs

    IF ( (variable(1:3) == 'kc_' .OR. variable(1:3) == 'em_')  )  THEN

       IF ( mode == 'allocate' )  THEN

          DO  lsp = 1, nspec
             IF ( TRIM( variable(1:3) ) == 'kc_' .AND. &
                  TRIM( variable(4:) ) == TRIM( chem_species(lsp)%name ) )  THEN
                chem_species(lsp)%conc_av = 0.0_wp
             ENDIF
          ENDDO

       ELSEIF ( mode == 'sum' )  THEN

          DO  lsp = 1, nspec
             IF ( TRIM( variable(1:3) ) == 'kc_' .AND. &
                  TRIM( variable(4:) ) == TRIM( chem_species(lsp)%name ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         chem_species(lsp)%conc_av(k,j,i) =                              &
                                           chem_species(lsp)%conc_av(k,j,i) +            &
                                           chem_species(lsp)%conc(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ELSEIF ( TRIM( variable(4:) ) == TRIM( 'cssws*' ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      match_def = surf_def_h(0)%start_index(j,i) <=                      &
                           surf_def_h(0)%end_index(j,i)
                      match_lsm = surf_lsm_h%start_index(j,i) <=                         &
                           surf_lsm_h%end_index(j,i)
                      match_usm = surf_usm_h%start_index(j,i) <=                         &
                           surf_usm_h%end_index(j,i)

                      IF ( match_def )  THEN
                         m = surf_def_h(0)%end_index(j,i)
                         chem_species(lsp)%cssws_av(j,i) =                               &
                              chem_species(lsp)%cssws_av(j,i) + &
                              surf_def_h(0)%cssws(lsp,m)
                      ELSEIF ( match_lsm  .AND.  .NOT. match_usm )  THEN
                         m = surf_lsm_h%end_index(j,i)
                         chem_species(lsp)%cssws_av(j,i) =                               &
                              chem_species(lsp)%cssws_av(j,i) + &
                              surf_lsm_h%cssws(lsp,m)
                      ELSEIF ( match_usm )  THEN
                         m = surf_usm_h%end_index(j,i)
                         chem_species(lsp)%cssws_av(j,i) =                               &
                              chem_species(lsp)%cssws_av(j,i) + &
                              surf_usm_h%cssws(lsp,m)
                      ENDIF
                   ENDDO
                ENDDO
             ENDIF
          ENDDO

       ELSEIF ( mode == 'average' )  THEN

          DO  lsp = 1, nspec
             IF ( TRIM( variable(1:3) ) == 'kc_' .AND. &
                  TRIM( variable(4:) ) == TRIM( chem_species(lsp)%name ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         chem_species(lsp)%conc_av(k,j,i) =                              &
                             chem_species(lsp)%conc_av(k,j,i) /                          &
                             REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO

             ELSEIF ( TRIM( variable(4:) ) == TRIM( 'cssws*' ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      chem_species(lsp)%cssws_av(j,i) =                                  &
                      chem_species(lsp)%cssws_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
                CALL exchange_horiz_2d( chem_species(lsp)%cssws_av )
             ENDIF
          ENDDO
       ENDIF

    ENDIF

 END SUBROUTINE chem_3d_data_averaging

    
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine to initialize and set all boundary conditions for chemical species
!------------------------------------------------------------------------------!
 SUBROUTINE chem_boundary_conds( mode )                                           

    USE control_parameters,                                                    &  
        ONLY:  bc_radiation_l, bc_radiation_n, bc_radiation_r, bc_radiation_s

    USE arrays_3d,                                                             &     
        ONLY:  dzu                                                

    USE surface_mod,                                                           &
        ONLY:  bc_h                                                           

    CHARACTER (LEN=*), INTENT(IN) ::  mode
    INTEGER(iwp) ::  i                            !< grid index x direction.
    INTEGER(iwp) ::  j                            !< grid index y direction.
    INTEGER(iwp) ::  k                            !< grid index z direction.
    INTEGER(iwp) ::  l                            !< running index boundary type, for up- and downward-facing walls.
    INTEGER(iwp) ::  m                            !< running index surface elements.
    INTEGER(iwp) ::  lsp                          !< running index for chem spcs.


    SELECT CASE ( TRIM( mode ) )        
       CASE ( 'init' )

          IF ( bc_cs_b == 'dirichlet' )  THEN
             ibc_cs_b = 0 
          ELSEIF ( bc_cs_b == 'neumann' )  THEN
             ibc_cs_b = 1 
          ELSE
             message_string = 'unknown boundary condition: bc_cs_b ="' // TRIM( bc_cs_b ) // '"'  
             CALL message( 'chem_boundary_conds', 'CM0429', 1, 2, 0, 6, 0 )
          ENDIF                                                                  
!
!--       Set Integer flags and check for possible erroneous settings for top
!--       boundary condition. 
          IF ( bc_cs_t == 'dirichlet' )  THEN
             ibc_cs_t = 0 
          ELSEIF ( bc_cs_t == 'neumann' )  THEN
             ibc_cs_t = 1
          ELSEIF ( bc_cs_t == 'initial_gradient' )  THEN
             ibc_cs_t = 2
          ELSEIF ( bc_cs_t == 'nested' )  THEN          
             ibc_cs_t = 3
          ELSE
             message_string = 'unknown boundary condition: bc_c_t ="' // TRIM( bc_cs_t ) // '"'      
             CALL message( 'check_parameters', 'CM0430', 1, 2, 0, 6, 0 )
          ENDIF

       CASE ( 'set_bc_bottomtop' )                   
!
!--       Boundary condtions for chemical species at horizontal walls      
          DO  lsp = 1, nspec                                                      
             IF ( ibc_cs_b == 0 )  THEN
                DO  l = 0, 1 
                   !$OMP PARALLEL DO PRIVATE( i, j, k )
                   DO  m = 1, bc_h(l)%ns
                       i = bc_h(l)%i(m)            
                       j = bc_h(l)%j(m)
                       k = bc_h(l)%k(m)
                      chem_species(lsp)%conc_p(k+bc_h(l)%koff,j,i) =           &
                                      chem_species(lsp)%conc(k+bc_h(l)%koff,j,i) 
                   ENDDO                                        
                ENDDO                                       

             ELSEIF ( ibc_cs_b == 1 )  THEN
!
!--             in boundary_conds there is som extra loop over m here for passive tracer
                DO  l = 0, 1
                   !$OMP PARALLEL DO PRIVATE( i, j, k )                                           
                   DO m = 1, bc_h(l)%ns
                      i = bc_h(l)%i(m)            
                      j = bc_h(l)%j(m)
                      k = bc_h(l)%k(m)
                      chem_species(lsp)%conc_p(k+bc_h(l)%koff,j,i) =           &
                                         chem_species(lsp)%conc_p(k,j,i)

                   ENDDO
                ENDDO
             ENDIF
       ENDDO    ! end lsp loop  
!
!--    Top boundary conditions for chemical species - Should this not be done for all species?
          IF ( ibc_cs_t == 0 )  THEN                    
             DO  lsp = 1, nspec
                chem_species(lsp)%conc_p(nzt+1,:,:) = chem_species(lsp)%conc(nzt+1,:,:)        
             ENDDO
          ELSEIF ( ibc_cs_t == 1 )  THEN
             DO  lsp = 1, nspec
                chem_species(lsp)%conc_p(nzt+1,:,:) = chem_species(lsp)%conc_p(nzt,:,:)
             ENDDO
          ELSEIF ( ibc_cs_t == 2 )  THEN
             DO  lsp = 1, nspec
                chem_species(lsp)%conc_p(nzt+1,:,:) = chem_species(lsp)%conc_p(nzt,:,:) + bc_cs_t_val(lsp) * dzu(nzt+1)
             ENDDO
          ENDIF

       CASE ( 'set_bc_lateral' )                       
!
!--             Lateral boundary conditions for chem species at inflow boundary
!--             are automatically set when chem_species concentration is
!--             initialized. The initially set value at the inflow boundary is not
!--             touched during time integration, hence, this boundary value remains
!--             at a constant value, which is the concentration that flows into the
!--             domain.                                                           
!--             Lateral boundary conditions for chem species at outflow boundary 

          IF ( bc_radiation_s )  THEN
             DO  lsp = 1, nspec
                chem_species(lsp)%conc_p(:,nys-1,:) = chem_species(lsp)%conc_p(:,nys,:)
             ENDDO
          ELSEIF ( bc_radiation_n )  THEN
             DO  lsp = 1, nspec
                chem_species(lsp)%conc_p(:,nyn+1,:) = chem_species(lsp)%conc_p(:,nyn,:)
             ENDDO
          ELSEIF ( bc_radiation_l )  THEN
             DO  lsp = 1, nspec
                chem_species(lsp)%conc_p(:,:,nxl-1) = chem_species(lsp)%conc_p(:,:,nxl)
             ENDDO
          ELSEIF ( bc_radiation_r )  THEN
             DO  lsp = 1, nspec
                chem_species(lsp)%conc_p(:,:,nxr+1) = chem_species(lsp)%conc_p(:,:,nxr)
             ENDDO
          ENDIF

    END SELECT

 END SUBROUTINE chem_boundary_conds

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine for boundary conditions
!------------------------------------------------------------------------------!
 SUBROUTINE chem_boundary_conditions

    IMPLICIT NONE

    INTEGER(iwp) ::  lsp             !<
    INTEGER(iwp) ::  lsp_usr         !<

!
!-- Top/bottom boundary conditions for chemical species
    CALL chem_boundary_conds( 'set_bc_bottomtop' )
!
!-- Lateral boundary conditions for chemical species
    CALL chem_boundary_conds( 'set_bc_lateral' )

!
!--  Boundary conditions for prognostic quantitites of other modules:
!--  Here, only decycling is carried out
     DO  lsp = 1, nvar
        lsp_usr = 1
        DO WHILE ( TRIM( cs_name( lsp_usr ) ) /= 'novalue' )
           IF ( TRIM( chem_species(lsp)%name ) == TRIM( cs_name(lsp_usr) ) )  THEN
              CALL chem_boundary_conds( chem_species(lsp)%conc_p,                          &
                                        chem_species(lsp)%conc_pr_init )
           ENDIF
           lsp_usr = lsp_usr + 1
        ENDDO
     ENDDO


 END SUBROUTINE chem_boundary_conditions

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Boundary conditions for prognostic variables in chemistry: decycling in the 
!> x-direction-
!> Decycling of chemistry variables: Dirichlet BCs with cyclic is frequently not
!> approproate for chemicals compounds since they may accumulate too much.
!> If no proper boundary conditions from a DYNAMIC input file are available,
!> de-cycling applies the initial profiles at the inflow boundaries for
!> Dirichlet boundary conditions
!------------------------------------------------------------------------------!
 SUBROUTINE chem_boundary_conds_decycle( cs_3d, cs_pr_init )

    USE control_parameters,                                                                         &
        ONLY:  nesting_offline

    INTEGER(iwp) ::  boundary  !< 
    INTEGER(iwp) ::  ee        !<
    INTEGER(iwp) ::  copied    !< 
    INTEGER(iwp) ::  i         !< 
    INTEGER(iwp) ::  j         !< 
    INTEGER(iwp) ::  k         !< 
    INTEGER(iwp) ::  ss        !<

    REAL(wp), DIMENSION(nzb:nzt+1) ::  cs_pr_init
    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  cs_3d
    REAL(wp) ::  flag          !< flag to mask topography grid points


    flag = 0.0_wp
!
!-- Skip input if forcing from a larger-scale model is applied
    IF ( nesting_offline  .AND.  nesting_offline_chem )  RETURN
!
!-- Left and right boundaries
    IF ( decycle_chem_lr  .AND.  bc_lr_cyc )  THEN

       DO  boundary = 1, 2

          IF ( decycle_method(boundary) == 'dirichlet' )  THEN
!
!--          Initial profile is copied to ghost and first three layers
             ss = 1
             ee = 0
             IF ( boundary == 1  .AND.  nxl == 0 )  THEN
                ss = nxlg
                ee = nxl-1
             ELSEIF ( boundary == 2  .AND.  nxr == nx )  THEN
                ss = nxr+1
                ee = nxrg
             ENDIF

             DO  i = ss, ee
                DO  j = nysg, nyng
                   DO  k = nzb+1, nzt
                      flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )
                      cs_3d(k,j,i) = cs_pr_init(k) * flag
                   ENDDO
                ENDDO
             ENDDO

          ELSEIF ( decycle_method(boundary) == 'neumann' )  THEN
!
!--          The value at the boundary is copied to the ghost layers to simulate
!--          an outlet with zero gradient
             ss = 1
             ee = 0
             IF ( boundary == 1  .AND.  nxl == 0 )  THEN
                ss = nxlg
                ee = nxl-1
                copied = nxl
             ELSEIF ( boundary == 2  .AND.  nxr == nx )  THEN
                ss = nxr+1
                ee = nxrg
                copied = nxr
             ENDIF

             DO  i = ss, ee
                DO  j = nysg, nyng
                   DO  k = nzb+1, nzt
                      flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )
                      cs_3d(k,j,i) = cs_3d(k,j,copied) * flag
                   ENDDO
                ENDDO
             ENDDO

          ELSE
             WRITE(message_string,*)                                           &
                  'unknown decycling method: decycle_method (', &
                  boundary, ') ="' // TRIM( decycle_method(boundary) ) // '"'
             CALL message( 'chem_boundary_conds_decycle', 'CM0431',           &
                  1, 2, 0, 6, 0 )
          ENDIF
       ENDDO
    ENDIF
!
!-- South and north boundaries
    IF ( decycle_chem_ns  .AND.  bc_ns_cyc )  THEN

       DO  boundary = 3, 4

          IF ( decycle_method(boundary) == 'dirichlet' )  THEN
!
!--          Initial profile is copied to ghost and first three layers
             ss = 1
             ee = 0
             IF ( boundary == 3  .AND.  nys == 0 )  THEN
                ss = nysg
                ee = nys-1
             ELSEIF ( boundary == 4  .AND.  nyn == ny )  THEN
                ss = nyn+1
                ee = nyng
             ENDIF

             DO  i = nxlg, nxrg
                DO  j = ss, ee
                   DO  k = nzb+1, nzt
                      flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )
                      cs_3d(k,j,i) = cs_pr_init(k) * flag
                   ENDDO
                ENDDO
             ENDDO


          ELSEIF ( decycle_method(boundary) == 'neumann' )  THEN
!
!--          The value at the boundary is copied to the ghost layers to simulate
!--          an outlet with zero gradient
             ss = 1
             ee = 0
             IF ( boundary == 3  .AND.  nys == 0 )  THEN
                ss = nysg
                ee = nys-1
                copied = nys
             ELSEIF ( boundary == 4  .AND.  nyn == ny )  THEN
                ss = nyn+1
                ee = nyng
                copied = nyn
             ENDIF

             DO  i = nxlg, nxrg
                DO  j = ss, ee
                   DO  k = nzb+1, nzt
                      flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )
                      cs_3d(k,j,i) = cs_3d(k,copied,i) * flag
                   ENDDO
                ENDDO
             ENDDO

          ELSE
             WRITE(message_string,*)                                           &
                  'unknown decycling method: decycle_method (', &
                  boundary, ') ="' // TRIM( decycle_method(boundary) ) // '"'
             CALL message( 'chem_boundary_conds_decycle', 'CM0432',           &
                  1, 2, 0, 6, 0 )
          ENDIF
       ENDDO
    ENDIF


 END SUBROUTINE chem_boundary_conds_decycle


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine for checking data output for chemical species
!------------------------------------------------------------------------------!
 SUBROUTINE chem_check_data_output( var, unit, i, ilen, k )


    CHARACTER (LEN=*) ::  unit     !<
    CHARACTER (LEN=*) ::  var      !<

    INTEGER(iwp) ::  i
    INTEGER(iwp) ::  lsp
    INTEGER(iwp) ::  ilen
    INTEGER(iwp) ::  k

    CHARACTER(LEN=16)    ::  spec_name

!
!-- Next statement is to avoid compiler warnings about unused variables
    IF ( ( i + ilen + k ) > 0  .OR.  var(1:1) == ' ' )  CONTINUE

    unit = 'illegal'

    spec_name = TRIM( var(4:) )             !< var 1:3 is 'kc_' or 'em_'.

    IF ( TRIM( var(1:3) ) == 'em_' )  THEN
       DO  lsp=1,nspec
          IF (TRIM( spec_name ) == TRIM( chem_species(lsp)%name ) )  THEN
             unit = 'mol m-2 s-1'
          ENDIF
!
!--       It is possible to plant PM10 and PM25 into the gasphase chemistry code
!--       as passive species (e.g. 'passive' in GASPHASE_PREPROC/mechanisms):
!--       set unit to micrograms per m**3 for PM10 and PM25 (PM2.5)
          IF (spec_name(1:2) == 'PM')  THEN
             unit = 'kg m-2 s-1'
          ENDIF
       ENDDO

    ELSE

       DO  lsp=1,nspec
          IF (TRIM( spec_name ) == TRIM( chem_species(lsp)%name ) )  THEN
             unit = 'ppm'
          ENDIF
!
!--            It is possible  to plant PM10 and PM25 into the gasphase chemistry code 
!--            as passive species (e.g. 'passive' in GASPHASE_PREPROC/mechanisms):
!--            set unit to kilograms per m**3 for PM10 and PM25 (PM2.5)
          IF (spec_name(1:2) == 'PM')  THEN  
            unit = 'kg m-3'
          ENDIF
       ENDDO

       DO  lsp=1,nphot
          IF (TRIM( spec_name ) == TRIM( phot_frequen(lsp)%name ) )  THEN
             unit = 'sec-1'
          ENDIF
       ENDDO
    ENDIF


    RETURN
 END SUBROUTINE chem_check_data_output


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine for checking data output of profiles for chemistry model
!------------------------------------------------------------------------------!
 SUBROUTINE chem_check_data_output_pr( variable, var_count, unit, dopr_unit )

    USE arrays_3d

    USE control_parameters,                                                    &
        ONLY:  data_output_pr, message_string

    USE profil_parameter

    USE statistics


    CHARACTER (LEN=*) ::  unit     !< 
    CHARACTER (LEN=*) ::  variable !< 
    CHARACTER (LEN=*) ::  dopr_unit
    CHARACTER (LEN=16) ::  spec_name

    INTEGER(iwp) ::  var_count, lsp  !<


    spec_name = TRIM( variable(4:) )    

       IF (  .NOT.  air_chemistry )  THEN
          message_string = 'data_output_pr = ' //                        &
          TRIM( data_output_pr(var_count) ) // ' is not imp' // &
                      'lemented for air_chemistry = .FALSE.'
          CALL message( 'chem_check_parameters', 'CM0433', 1, 2, 0, 6, 0 )              

       ELSE
          DO  lsp = 1, nspec
             IF (TRIM( spec_name ) == TRIM( chem_species(lsp)%name ) )  THEN 
                cs_pr_count = cs_pr_count+1
                cs_pr_index(cs_pr_count) = lsp
                dopr_index(var_count) = pr_palm+cs_pr_count  
                dopr_unit  = 'ppm'
                IF (spec_name(1:2) == 'PM')  THEN
                     dopr_unit = 'kg m-3'
                ENDIF
                   hom(:,2, dopr_index(var_count),:) = SPREAD( zu, 2, statistic_regions+1 )
                   unit = dopr_unit  
             ENDIF
          ENDDO
       ENDIF

 END SUBROUTINE chem_check_data_output_pr


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check parameters routine for chemistry_model_mod
!------------------------------------------------------------------------------!
 SUBROUTINE chem_check_parameters


    LOGICAL  ::  found
    INTEGER (iwp) ::  lsp_usr      !< running index for user defined chem spcs 
    INTEGER (iwp) ::  lsp          !< running index for chem spcs.
!
!-- check for chemical reactions status
    IF ( chem_gasphase_on )  THEN
       message_string = 'Chemical reactions: ON'
       CALL message( 'chem_check_parameters', 'CM0421', 0, 0, 0, 6, 0 )
    ELSEIF ( .NOT. (chem_gasphase_on) )  THEN
       message_string = 'Chemical reactions: OFF'
       CALL message( 'chem_check_parameters', 'CM0422', 0, 0, 0, 6, 0 )
    ENDIF
!
!-- check for chemistry time-step
    IF ( call_chem_at_all_substeps )  THEN
       message_string = 'Chemistry is calculated at all meteorology time-step'
       CALL message( 'chem_check_parameters', 'CM0423', 0, 0, 0, 6, 0 )
    ELSEIF ( .not. (call_chem_at_all_substeps) )  THEN
       message_string = 'Sub-time-steps are skipped for chemistry time-steps'
       CALL message( 'chem_check_parameters', 'CM0424', 0, 0, 0, 6, 0 )
    ENDIF
!
!-- check for photolysis scheme 
    IF ( (photolysis_scheme /= 'simple') .AND. (photolysis_scheme /= 'constant')  )  THEN
       message_string = 'Incorrect photolysis scheme selected, please check spelling'
       CALL message( 'chem_check_parameters', 'CM0425', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- check for  decycling of chem species
    IF ( (.NOT. any(decycle_method == 'neumann') ) .AND. (.NOT. any(decycle_method == 'dirichlet') ) )  THEN
       message_string = 'Incorrect boundary conditions. Only neumann or '   &
                // 'dirichlet &available for decycling chemical species '
       CALL message( 'chem_check_parameters', 'CM0426', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- check for chemical mechanism used 
    CALL get_mechanism_name
    IF ( chem_mechanism /= TRIM( cs_mech ) )  THEN
       message_string = 'Incorrect chemistry mechanism selected, check spelling in namelist and/or chem_gasphase_mod'
       CALL message( 'chem_check_parameters', 'CM0462', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- If nesting_chem = .F., set top boundary condition to its default value
    IF ( .NOT. nesting_chem  .AND.  ibc_cs_t == 3  )  THEN
       ibc_cs_t = 2
       bc_cs_t = 'initial_gradient'
    ENDIF
!
!-- chem_check_parameters is called before the array chem_species is allocated!
!-- temporary switch of this part of the check
!    RETURN                !bK commented

    CALL chem_init_internal
!
!-- check for initial chem species input
    lsp_usr = 1
    lsp     = 1
    DO WHILE ( cs_name (lsp_usr) /= 'novalue')
       found = .FALSE.
       DO  lsp = 1, nvar
          IF ( TRIM( cs_name (lsp_usr) ) == TRIM( chem_species(lsp)%name) )  THEN
             found = .TRUE.
             EXIT
          ENDIF 
       ENDDO
       IF ( .NOT.  found )  THEN
          message_string = 'Unused/incorrect input for initial surface value: ' //     &
                            TRIM( cs_name(lsp_usr) )
          CALL message( 'chem_check_parameters', 'CM0427', 1, 2, 0, 6, 0 )
       ENDIF
       lsp_usr = lsp_usr + 1
    ENDDO
!
!-- check for surface  emission flux chem species
    lsp_usr = 1
    lsp     = 1
    DO WHILE ( surface_csflux_name (lsp_usr) /= 'novalue')
       found = .FALSE.
       DO  lsp = 1, nvar
          IF ( TRIM( surface_csflux_name (lsp_usr) ) == TRIM( chem_species(lsp)%name ) )  THEN
             found = .TRUE.
             EXIT
          ENDIF 
       ENDDO
       IF ( .NOT.  found )  THEN
          message_string = 'Unused/incorrect input of chemical species for surface emission fluxes: '  &
                            // TRIM( surface_csflux_name(lsp_usr) )
          CALL message( 'chem_check_parameters', 'CM0428', 1, 2, 0, 6, 0 )
       ENDIF
       lsp_usr = lsp_usr + 1
    ENDDO

 END SUBROUTINE chem_check_parameters


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine defining 2D output variables for chemical species
!> @todo: Remove "mode" from argument list, not used.
!------------------------------------------------------------------------------!
 SUBROUTINE chem_data_output_2d( av, variable, found, grid, mode, local_pf,   &
                                  two_d, nzb_do, nzt_do, fill_value )


    CHARACTER (LEN=*) ::  grid       !<
    CHARACTER (LEN=*) ::  mode       !< 
    CHARACTER (LEN=*) ::  variable   !< 
    INTEGER(iwp) ::  av              !< flag to control data output of instantaneous or time-averaged data
    INTEGER(iwp) ::  nzb_do          !< lower limit of the domain (usually nzb)
    INTEGER(iwp) ::  nzt_do          !< upper limit of the domain (usually nzt+1)
    LOGICAL      ::  found           !< 
    LOGICAL      ::  two_d           !< flag parameter that indicates 2D variables (horizontal cross sections)
    REAL(wp)     ::  fill_value 
    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb:nzt+1) ::  local_pf  

!
!-- local variables.
    CHARACTER(LEN=16)    ::  spec_name
    INTEGER(iwp) ::  lsp
    INTEGER(iwp) ::  i               !< grid index along x-direction
    INTEGER(iwp) ::  j               !< grid index along y-direction
    INTEGER(iwp) ::  k               !< grid index along z-direction
    INTEGER(iwp) ::  m               !< running indices for surfaces
    INTEGER(iwp) ::  char_len        !< length of a character string
!
!-- Next statement is to avoid compiler warnings about unused variables
    IF ( mode(1:1) == ' '  .OR.  two_d )  CONTINUE

    found = .FALSE.
    char_len  = LEN_TRIM( variable )

    spec_name = TRIM( variable(4:char_len-3) )
!
!-- Output of emission values, i.e. surface fluxes cssws. 
    IF ( variable(1:3) == 'em_' )  THEN

       local_pf = 0.0_wp

       DO  lsp = 1, nvar
          IF ( TRIM( spec_name ) == TRIM( chem_species(lsp)%name) )  THEN
!
!--          No average output for now. 
             DO  m = 1, surf_lsm_h%ns
                local_pf(surf_lsm_h%i(m),surf_lsm_h%j(m),nzb+1) =              &
                              local_pf(surf_lsm_h%i(m),surf_lsm_h%j(m),nzb+1)  &
                                            + surf_lsm_h%cssws(lsp,m)
             ENDDO
             DO  m = 1, surf_usm_h%ns
                   local_pf(surf_usm_h%i(m),surf_usm_h%j(m),nzb+1) =           &
                              local_pf(surf_usm_h%i(m),surf_usm_h%j(m),nzb+1)  &
                                            + surf_usm_h%cssws(lsp,m)
             ENDDO
             grid = 'zu'
             found = .TRUE.
          ENDIF
       ENDDO

    ELSE

       DO  lsp=1,nspec
          IF (TRIM( spec_name ) == TRIM( chem_species(lsp)%name )  .AND.                           &
                ( (variable(char_len-2:) == '_xy')  .OR.                                           &
                  (variable(char_len-2:) == '_xz')  .OR.                                           &
                  (variable(char_len-2:) == '_yz') ) )  THEN             
!
!--   todo: remove or replace by "CALL message" mechanism (kanani)
!                    IF(myid == 0)  WRITE(6,*) 'Output of species ' // TRIM( variable )  //       &
!                                                             TRIM( chem_species(lsp)%name )       
             IF (av == 0)  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                           local_pf(i,j,k) = MERGE(                                                &
                                              chem_species(lsp)%conc(k,j,i),                       &
                                              REAL( fill_value, KIND = wp ),                       &
                                              BTEST( wall_flags_total_0(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO

             ELSE
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                           local_pf(i,j,k) = MERGE(                                                &
                                              chem_species(lsp)%conc_av(k,j,i),                    &
                                              REAL( fill_value, KIND = wp ),                       &
                                              BTEST( wall_flags_total_0(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
             grid = 'zu'           
             found = .TRUE.
          ENDIF
       ENDDO
    ENDIF

    RETURN

 END SUBROUTINE chem_data_output_2d      


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine defining 3D output variables for chemical species
!------------------------------------------------------------------------------!
 SUBROUTINE chem_data_output_3d( av, variable, found, local_pf, fill_value, nzb_do, nzt_do )


    USE surface_mod

    CHARACTER (LEN=*)    ::  variable     !<
    INTEGER(iwp)         ::  av           !<
    INTEGER(iwp) ::  nzb_do               !< lower limit of the data output (usually 0)
    INTEGER(iwp) ::  nzt_do               !< vertical upper limit of the data output (usually nz_do3d)

    LOGICAL      ::  found                !<

    REAL(wp)             ::  fill_value   !<
    REAL(sp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf 
!
!-- local variables
    CHARACTER(LEN=16)    ::  spec_name
    INTEGER(iwp)         ::  i
    INTEGER(iwp)         ::  j
    INTEGER(iwp)         ::  k
    INTEGER(iwp)         ::  m       !< running indices for surfaces
    INTEGER(iwp)         ::  l
    INTEGER(iwp)         ::  lsp     !< running index for chem spcs


    found = .FALSE.
    IF ( .NOT. (variable(1:3) == 'kc_' .OR. variable(1:3) == 'em_' ) )  THEN
       RETURN
    ENDIF

    spec_name = TRIM( variable(4:) )

    IF ( variable(1:3) == 'em_' )  THEN

       DO  lsp = 1, nvar   !!! cssws - nvar species, chem_species - nspec species !!!
          IF ( TRIM( spec_name ) == TRIM( chem_species(lsp)%name) )  THEN
          
             local_pf = 0.0_wp
!
!--          no average for now
             DO  m = 1, surf_usm_h%ns
                local_pf(surf_usm_h%i(m),surf_usm_h%j(m),surf_usm_h%k(m)) = &
                   local_pf(surf_usm_h%i(m),surf_usm_h%j(m),surf_usm_h%k(m)) + surf_usm_h%cssws(lsp,m)
             ENDDO
             DO  m = 1, surf_lsm_h%ns
                local_pf(surf_lsm_h%i(m),surf_lsm_h%j(m),surf_lsm_h%k(m)) = &
                  local_pf(surf_lsm_h%i(m),surf_lsm_h%j(m),surf_lsm_h%k(m)) + surf_lsm_h%cssws(lsp,m)
             ENDDO
             DO  l = 0, 3
                DO  m = 1, surf_usm_v(l)%ns
                   local_pf(surf_usm_v(l)%i(m),surf_usm_v(l)%j(m),surf_usm_v(l)%k(m)) = &
                     local_pf(surf_usm_v(l)%i(m),surf_usm_v(l)%j(m),surf_usm_v(l)%k(m)) + surf_usm_v(l)%cssws(lsp,m)
                ENDDO
                DO  m = 1, surf_lsm_v(l)%ns
                   local_pf(surf_lsm_v(l)%i(m),surf_lsm_v(l)%j(m),surf_lsm_v(l)%k(m)) = &
                      local_pf(surf_lsm_v(l)%i(m),surf_lsm_v(l)%j(m),surf_lsm_v(l)%k(m)) + surf_lsm_v(l)%cssws(lsp,m)
                ENDDO
             ENDDO
             found = .TRUE.
          ENDIF
       ENDDO
    ELSE
      DO  lsp = 1, nspec
         IF (TRIM( spec_name ) == TRIM( chem_species(lsp)%name) )  THEN
!
!--   todo: remove or replace by "CALL message" mechanism (kanani)
!              IF(myid == 0 .AND. chem_debug0 )  WRITE(6,*) 'Output of species ' // TRIM( variable )  // &
!                                                           TRIM( chem_species(lsp)%name )       
            IF (av == 0)  THEN
               DO  i = nxl, nxr
                  DO  j = nys, nyn
                     DO  k = nzb_do, nzt_do
                         local_pf(i,j,k) = MERGE(                             &
                                             chem_species(lsp)%conc(k,j,i),   &
                                             REAL( fill_value, KIND = wp ),   &
                                             BTEST( wall_flags_total_0(k,j,i), 0 ) )
                     ENDDO
                  ENDDO
               ENDDO

            ELSE

               DO  i = nxl, nxr
                  DO  j = nys, nyn
                     DO  k = nzb_do, nzt_do
                         local_pf(i,j,k) = MERGE(                             &
                                             chem_species(lsp)%conc_av(k,j,i),&
                                             REAL( fill_value, KIND = wp ),   &
                                             BTEST( wall_flags_total_0(k,j,i), 0 ) )
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF
            found = .TRUE.
         ENDIF
      ENDDO
    ENDIF

    RETURN

 END SUBROUTINE chem_data_output_3d


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine defining mask output variables for chemical species
!------------------------------------------------------------------------------!
 SUBROUTINE chem_data_output_mask( av, variable, found, local_pf, mid )


    USE control_parameters

    CHARACTER(LEN=16) ::  spec_name
    CHARACTER(LEN=*)  ::  variable    !<

    INTEGER(iwp) ::  av              !< flag to control data output of instantaneous or time-averaged data
    INTEGER(iwp) ::  lsp
    INTEGER(iwp) ::  i               !< grid index along x-direction
    INTEGER(iwp) ::  j               !< grid index along y-direction
    INTEGER(iwp) ::  k               !< grid index along z-direction
    INTEGER(iwp) ::  im              !< loop index for masked variables
    INTEGER(iwp) ::  jm              !< loop index for masked variables
    INTEGER(iwp) ::  kk              !< masked output index along z-direction
    INTEGER(iwp) ::  mid             !< masked output running index
    INTEGER(iwp) ::  ktt             !< k index of highest terrain surface
    
    LOGICAL ::  found
    
    REAL(wp), DIMENSION(mask_size_l(mid,1),mask_size_l(mid,2),mask_size_l(mid,3)) ::  &
              local_pf   !< 

    REAL(wp), PARAMETER ::  fill_value = -9999.0_wp    !< value for the _FillValue attribute

!
!-- local variables.

    spec_name = TRIM( variable(4:) )
    found = .FALSE.

    DO  lsp=1,nspec
       IF (TRIM( spec_name ) == TRIM( chem_species(lsp)%name) )  THEN             
!
!-- todo: remove or replace by "CALL message" mechanism (kanani)
!              IF(myid == 0 .AND. chem_debug0 )  WRITE(6,*) 'Output of species ' // TRIM( variable )  // &
!                                                        TRIM( chem_species(lsp)%name )       
          IF (av == 0)  THEN
             IF ( .NOT. mask_surface(mid) )  THEN

                DO  i = 1, mask_size_l(mid,1)
                   DO  j = 1, mask_size_l(mid,2) 
                      DO  k = 1, mask_size(mid,3) 
                          local_pf(i,j,k) = chem_species(lsp)%conc(  &
                                               mask_k(mid,k),        &
                                               mask_j(mid,j),        &
                                               mask_i(mid,i)      )
                      ENDDO
                   ENDDO
                ENDDO

             ELSE
!              
!--             Terrain-following masked output
                DO  i = 1, mask_size_l(mid,1)
                   DO  j = 1, mask_size_l(mid,2)
!--                   Get k index of the highest terraing surface
                      im = mask_i(mid,i)
                      jm = mask_j(mid,j)
                      ktt = MINLOC( MERGE( 1, 0, BTEST( wall_flags_total_0(:,jm,im), 5 )), DIM = 1 ) - 1
                      DO  k = 1, mask_size_l(mid,3)
                         kk = MIN( ktt+mask_k(mid,k), nzt+1 )
!--                      Set value if not in building
                         IF ( BTEST( wall_flags_total_0(kk,jm,im), 6 ) )  THEN
                            local_pf(i,j,k) = fill_value
                         ELSE
                            local_pf(i,j,k) = chem_species(lsp)%conc(kk,jm,im)
                         ENDIF
                      ENDDO
                   ENDDO
                ENDDO

             ENDIF
          ELSE
             IF ( .NOT. mask_surface(mid) )  THEN

                DO  i = 1, mask_size_l(mid,1)
                   DO  j = 1, mask_size_l(mid,2)
                      DO  k =  1, mask_size_l(mid,3)
                          local_pf(i,j,k) = chem_species(lsp)%conc_av(  &
                                               mask_k(mid,k),           &
                                               mask_j(mid,j),           &
                                               mask_i(mid,i)         )
                      ENDDO
                   ENDDO
                ENDDO

             ELSE
!              
!--             Terrain-following masked output
                DO  i = 1, mask_size_l(mid,1)
                   DO  j = 1, mask_size_l(mid,2)
!--                   Get k index of the highest terraing surface
                      im = mask_i(mid,i)
                      jm = mask_j(mid,j)
                      ktt = MINLOC( MERGE( 1, 0, BTEST( wall_flags_total_0(:,jm,im), 5 )), DIM = 1 ) - 1
                      DO  k = 1, mask_size_l(mid,3)
                         kk = MIN( ktt+mask_k(mid,k), nzt+1 )
!--                      Set value if not in building
                         IF ( BTEST( wall_flags_total_0(kk,jm,im), 6 ) )  THEN
                            local_pf(i,j,k) = fill_value
                         ELSE
                            local_pf(i,j,k) = chem_species(lsp)%conc_av(kk,jm,im)
                         ENDIF
                      ENDDO
                   ENDDO
                ENDDO

             ENDIF

          ENDIF
          found = .TRUE.
          EXIT
       ENDIF
    ENDDO

    RETURN

 END SUBROUTINE chem_data_output_mask      


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine defining appropriate grid for netcdf variables.
!> It is called out from subroutine netcdf.
!------------------------------------------------------------------------------!
 SUBROUTINE chem_define_netcdf_grid( var, found, grid_x, grid_y, grid_z )


    CHARACTER (LEN=*), INTENT(IN)  ::  var          !<
    LOGICAL, INTENT(OUT)           ::  found        !<
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_x       !<
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_y       !<
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_z       !<

    found  = .TRUE.

    IF ( var(1:3) == 'kc_' .OR. var(1:3) == 'em_' )  THEN                   !< always the same grid for chemistry variables
       grid_x = 'x'
       grid_y = 'y'
       grid_z = 'zu'                             
    ELSE
       found  = .FALSE.
       grid_x = 'none'
       grid_y = 'none'
       grid_z = 'none'
    ENDIF


 END SUBROUTINE chem_define_netcdf_grid


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine defining header output for chemistry model
!------------------------------------------------------------------------------!
 SUBROUTINE chem_header( io )


    INTEGER(iwp), INTENT(IN) ::  io            !< Unit of the output file
    INTEGER(iwp)  :: lsp                       !< running index for chem spcs
    INTEGER(iwp)  :: cs_fixed 
    CHARACTER (LEN=80)  :: docsflux_chr
    CHARACTER (LEN=80)  :: docsinit_chr
!
! Get name of chemical mechanism from chem_gasphase_mod
    CALL get_mechanism_name
!
!-- Write chemistry model  header
    WRITE( io, 1 )
!
!-- Gasphase reaction status
    IF ( chem_gasphase_on )  THEN 
       WRITE( io, 2 )
    ELSE
       WRITE( io, 3 )
    ENDIF
!
!-- Chemistry time-step
    WRITE ( io, 4 ) cs_time_step
!
!-- Emission mode info
!-- At the moment the evaluation is done with both emiss_lod and mode_emis
!-- but once salsa has been migrated to emiss_lod the .OR. mode_emis
!-- conditions can be removed (ecc 20190513)
    IF     ( (emiss_lod == 1) .OR. (mode_emis == 'DEFAULT') )        THEN
        WRITE ( io, 5 )
    ELSEIF ( (emiss_lod == 0) .OR. (mode_emis == 'PARAMETERIZED') )  THEN 
        WRITE ( io, 6 )
    ELSEIF ( (emiss_lod == 2) .OR. (mode_emis == 'PRE-PROCESSED') )  THEN 
        WRITE ( io, 7 )
    ENDIF
!
!-- Photolysis scheme info
    IF ( photolysis_scheme == "simple" )  THEN
       WRITE( io, 8 ) 
    ELSEIF (photolysis_scheme == "constant" )  THEN
       WRITE( io, 9 )
    ENDIF
!
!-- Emission flux info
    lsp = 1
    docsflux_chr ='Chemical species for surface emission flux: ' 
    DO WHILE ( surface_csflux_name(lsp) /= 'novalue' )
       docsflux_chr = TRIM( docsflux_chr ) // ' ' // TRIM( surface_csflux_name(lsp) ) // ','  
       IF ( LEN_TRIM( docsflux_chr ) >= 75 )  THEN
          WRITE ( io, 10 )  docsflux_chr
          docsflux_chr = '       '
       ENDIF
       lsp = lsp + 1
    ENDDO

    IF ( docsflux_chr /= '' )  THEN
       WRITE ( io, 10 )  docsflux_chr
    ENDIF
!
!-- initializatoin of Surface and profile chemical species

    lsp = 1
    docsinit_chr ='Chemical species for initial surface and profile emissions: ' 
    DO WHILE ( cs_name(lsp) /= 'novalue' )
       docsinit_chr = TRIM( docsinit_chr ) // ' ' // TRIM( cs_name(lsp) ) // ','  
       IF ( LEN_TRIM( docsinit_chr ) >= 75 )  THEN
        WRITE ( io, 11 )  docsinit_chr
        docsinit_chr = '       '
       ENDIF
       lsp = lsp + 1
    ENDDO

    IF ( docsinit_chr /= '' )  THEN
       WRITE ( io, 11 )  docsinit_chr
    ENDIF

    IF ( nesting_chem )  WRITE( io, 12 )  nesting_chem
    IF ( nesting_offline_chem )  WRITE( io, 13 )  nesting_offline_chem
!
!-- number of variable and fix chemical species and number of reactions
    cs_fixed = nspec - nvar

    WRITE ( io, * ) '   --> Chemical Mechanism        : ', cs_mech 
    WRITE ( io, * ) '   --> Chemical species, variable: ', nvar
    WRITE ( io, * ) '   --> Chemical species, fixed   : ', cs_fixed
    WRITE ( io, * ) '   --> Total number of reactions : ', nreact


1   FORMAT (//' Chemistry model information:'/                                  &
           ' ----------------------------'/)
2   FORMAT ('    --> Chemical reactions are turned on')
3   FORMAT ('    --> Chemical reactions are turned off')
4   FORMAT ('    --> Time-step for chemical species: ',F6.2, ' s')
5   FORMAT ('    --> Emission mode = DEFAULT ')
6   FORMAT ('    --> Emission mode = PARAMETERIZED ')
7   FORMAT ('    --> Emission mode = PRE-PROCESSED ')
8   FORMAT ('    --> Photolysis scheme used =  simple ')
9   FORMAT ('    --> Photolysis scheme used =  constant ')
10  FORMAT (/'    ',A)  
11  FORMAT (/'    ',A) 
12  FORMAT (/'   Nesting for chemistry variables: ', L1 )
13  FORMAT (/'   Offline nesting for chemistry variables: ', L1 )
!
!
 END SUBROUTINE chem_header


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine initializating chemistry_model_mod specific arrays
!------------------------------------------------------------------------------!
 SUBROUTINE chem_init_arrays
!
!-- Please use this place to allocate required arrays

 END SUBROUTINE chem_init_arrays


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine initializating chemistry_model_mod
!------------------------------------------------------------------------------!
 SUBROUTINE chem_init

!
!-- 20200203 (ECC)
!-- introduced additional interfaces for on-demand emission update

!    USE chem_emissions_mod,                                                    &
!        ONLY:  chem_emissions_init
        
    USE chem_emissions_mod,                                                    &
        ONLY:  chem_emissions_init, chem_emissions_header_init
        
    USE netcdf_data_input_mod,                                                 &
        ONLY:  init_3d


    INTEGER(iwp) ::  i !< running index x dimension
    INTEGER(iwp) ::  j !< running index y dimension
    INTEGER(iwp) ::  n !< running index for chemical species


    IF ( debug_output )  CALL debug_message( 'chem_init', 'start' )
!
!-- Next statement is to avoid compiler warning about unused variables
    IF ( ( ilu_arable + ilu_coniferous_forest + ilu_deciduous_forest + ilu_mediterrean_scrub + &
           ilu_permanent_crops + ilu_savanna + ilu_semi_natural_veg + ilu_tropical_forest +    &
           ilu_urban ) == 0 )  CONTINUE

!
!-- 20200203 (ECC)
!-- calls specific emisisons initialization subroutines
!-- for legacy mode and on-demand mode         

!    IF ( emissions_anthropogenic )  CALL chem_emissions_init

    IF  ( emissions_anthropogenic )  THEN

       IF  ( emiss_read_legacy_mode )  THEN
          CALL chem_emissions_init
       ELSE
          CALL chem_emissions_header_init
       ENDIF

    ENDIF


!
!-- Chemistry variables will be initialized if availabe from dynamic 
!-- input file. Note, it is possible to initialize only part of the chemistry
!-- variables from dynamic input. 
    IF ( INDEX( initializing_actions, 'inifor' ) /= 0 )  THEN
       DO  n = 1, nspec
          IF ( init_3d%from_file_chem(n) )  THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   chem_species(n)%conc(:,j,i) = init_3d%chem_init(:,n)
                ENDDO
             ENDDO
          ENDIF
       ENDDO
    ENDIF

    IF ( debug_output )  CALL debug_message( 'chem_init', 'end' )

 END SUBROUTINE chem_init


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine initializating chemistry_model_mod
!> internal workaround for chem_species dependency in chem_check_parameters
!------------------------------------------------------------------------------!
 SUBROUTINE chem_init_internal
 
    USE pegrid

    USE netcdf_data_input_mod,                                                 &
        ONLY:  chem_emis, chem_emis_att, input_pids_dynamic, init_3d,          &
               netcdf_data_input_chemistry_data

!
!-- Local variables
    INTEGER(iwp) ::  i                 !< running index for for horiz numerical grid points
    INTEGER(iwp) ::  j                 !< running index for for horiz numerical grid points
    INTEGER(iwp) ::  lsp               !< running index for chem spcs
    INTEGER(iwp) ::  lpr_lev           !< running index for chem spcs profile level

!
!-- 20200203 ECC
!-- reads netcdf data only under legacy mode

!    IF ( emissions_anthropogenic )  THEN
!       CALL netcdf_data_input_chemistry_data( chem_emis_att, chem_emis )
!    ENDIF

    IF ( emissions_anthropogenic )  THEN
       IF ( emiss_read_legacy_mode )  THEN
          CALL netcdf_data_input_chemistry_data( chem_emis_att, chem_emis )
       ENDIF
    ENDIF

!
!-- Allocate memory for chemical species
    ALLOCATE( chem_species(nspec) )
    ALLOCATE( spec_conc_1 (nzb:nzt+1,nysg:nyng,nxlg:nxrg,nspec) )
    ALLOCATE( spec_conc_2 (nzb:nzt+1,nysg:nyng,nxlg:nxrg,nspec) )
    ALLOCATE( spec_conc_3 (nzb:nzt+1,nysg:nyng,nxlg:nxrg,nspec) )
    ALLOCATE( spec_conc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg,nspec) ) 
    ALLOCATE( phot_frequen(nphot) ) 
    ALLOCATE( freq_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg,nphot) )
    ALLOCATE( bc_cs_t_val(nspec) )
!
!-- Initialize arrays
    spec_conc_1 (:,:,:,:) = 0.0_wp
    spec_conc_2 (:,:,:,:) = 0.0_wp
    spec_conc_3 (:,:,:,:) = 0.0_wp
    spec_conc_av(:,:,:,:) = 0.0_wp


    DO  lsp = 1, nspec
       chem_species(lsp)%name    = spc_names(lsp)

       chem_species(lsp)%conc   (nzb:nzt+1,nysg:nyng,nxlg:nxrg)       => spec_conc_1 (:,:,:,lsp)
       chem_species(lsp)%conc_p (nzb:nzt+1,nysg:nyng,nxlg:nxrg)       => spec_conc_2 (:,:,:,lsp)
       chem_species(lsp)%tconc_m(nzb:nzt+1,nysg:nyng,nxlg:nxrg)       => spec_conc_3 (:,:,:,lsp)
       chem_species(lsp)%conc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg)       => spec_conc_av(:,:,:,lsp)     

       ALLOCATE (chem_species(lsp)%cssws_av(nysg:nyng,nxlg:nxrg))                   
       chem_species(lsp)%cssws_av    = 0.0_wp
!
!--    The following block can be useful when emission module is not applied. &
!--    if emission module is applied the following block will be overwritten.
       ALLOCATE (chem_species(lsp)%flux_s_cs(nzb+1:nzt,0:threads_per_task-1))   
       ALLOCATE (chem_species(lsp)%diss_s_cs(nzb+1:nzt,0:threads_per_task-1))   
       ALLOCATE (chem_species(lsp)%flux_l_cs(nzb+1:nzt,nys:nyn,0:threads_per_task-1)) 
       ALLOCATE (chem_species(lsp)%diss_l_cs(nzb+1:nzt,nys:nyn,0:threads_per_task-1))   
       chem_species(lsp)%flux_s_cs = 0.0_wp                                     
       chem_species(lsp)%flux_l_cs = 0.0_wp                                     
       chem_species(lsp)%diss_s_cs = 0.0_wp                                      
       chem_species(lsp)%diss_l_cs = 0.0_wp                                     
!
!--   Allocate memory for initial concentration profiles
!--   (concentration values come from namelist)
!--   (@todo (FK): Because of this, chem_init is called in palm before
!--               check_parameters, since conc_pr_init is used there.
!--               We have to find another solution since chem_init should
!--               eventually be called from init_3d_model!!)
       ALLOCATE ( chem_species(lsp)%conc_pr_init(0:nz+1) )
       chem_species(lsp)%conc_pr_init(:) = 0.0_wp

    ENDDO
!
!-- Set control flags for decycling only at lateral boundary cores, within the 
!-- inner cores the decycle flag is set to .False.. Even though it does not 
!-- affect the setting of chemistry boundary conditions, this flag is used to
!-- set advection control flags appropriately. 
    decycle_chem_lr = MERGE( decycle_chem_lr, .FALSE.,                         &
                             nxl == 0  .OR.  nxr == nx )
    decycle_chem_ns = MERGE( decycle_chem_ns, .FALSE.,                         &
                             nys == 0  .OR.  nyn == ny )
!
!-- For some passive scalars decycling may be enabled. This case, the lateral 
!-- boundary conditions are non-cyclic for these scalars (chemical species
!-- and aerosols), while the other scalars may have
!-- cyclic boundary conditions. However, large gradients near the boundaries
!-- may produce stationary numerical oscillations near the lateral boundaries
!-- when a higher-order scheme is applied near these boundaries. 
!-- To get rid-off this, set-up additional flags that control the order of the 
!-- scalar advection scheme near the lateral boundaries for passive scalars 
!-- with decycling.
    IF ( scalar_advec == 'ws-scheme' )  THEN
       ALLOCATE( cs_advc_flags_s(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
!
!--    In case of decyling, set Neumann boundary conditions for wall_flags_total_0 
!--    bit 31 instead of cyclic boundary conditions.
!--    Bit 31 is used to identify extended degradation zones (please see
!--    following comment). 
!--    Note, since several also other modules like Salsa or other future 
!--    one may access this bit but may have other boundary conditions, the 
!--    original value of wall_flags_total_0 bit 31 must not be modified. Hence, 
!--    store the boundary conditions directly on cs_advc_flags_s. 
!--    cs_advc_flags_s will be later overwritten in ws_init_flags_scalar and
!--    bit 31 won't be used to control the numerical order. 
!--    Initialize with flag 31 only. 
       cs_advc_flags_s = 0
       cs_advc_flags_s = MERGE( IBSET( cs_advc_flags_s, 31 ), 0,               &
                                BTEST( wall_flags_total_0, 31 ) )

       IF ( decycle_chem_ns )  THEN
          IF ( nys == 0  )  THEN
             DO  i = 1, nbgp     
                cs_advc_flags_s(:,nys-i,:) = MERGE(                            &
                                        IBSET( cs_advc_flags_s(:,nys,:), 31 ), &
                                        IBCLR( cs_advc_flags_s(:,nys,:), 31 ), &
                                        BTEST( cs_advc_flags_s(:,nys,:), 31 )  &
                                                  )
             ENDDO
          ENDIF
          IF ( nyn == ny )  THEN
             DO  i = 1, nbgp  
                cs_advc_flags_s(:,nyn+i,:) = MERGE(                            &
                                        IBSET( cs_advc_flags_s(:,nyn,:), 31 ), &
                                        IBCLR( cs_advc_flags_s(:,nyn,:), 31 ), &
                                        BTEST( cs_advc_flags_s(:,nyn,:), 31 )  &
                                                  )
             ENDDO
          ENDIF
       ENDIF
       IF ( decycle_chem_lr )  THEN
          IF ( nxl == 0  )  THEN
             DO  i = 1, nbgp   
                cs_advc_flags_s(:,:,nxl-i) = MERGE(                            &
                                        IBSET( cs_advc_flags_s(:,:,nxl), 31 ), &
                                        IBCLR( cs_advc_flags_s(:,:,nxl), 31 ), &
                                        BTEST( cs_advc_flags_s(:,:,nxl), 31 )  &
                                                  )
             ENDDO
          ENDIF
          IF ( nxr == nx )  THEN 
             DO  i = 1, nbgp   
                cs_advc_flags_s(:,:,nxr+i) = MERGE(                            &
                                        IBSET( cs_advc_flags_s(:,:,nxr), 31 ), &
                                        IBCLR( cs_advc_flags_s(:,:,nxr), 31 ), &
                                        BTEST( cs_advc_flags_s(:,:,nxr), 31 )  &
                                                  )
             ENDDO
          ENDIF      
       ENDIF
!
!--    To initialize advection flags appropriately, pass the boundary flags. 
!--    The last argument indicates that a passive scalar is treated, where
!--    the horizontal advection terms are degraded already 2 grid points before
!--    the lateral boundary to avoid stationary oscillations at large-gradients. 
!--    Also, extended degradation zones are applied, where horizontal advection of 
!--    passive scalars is discretized by first-order scheme at all grid points 
!--    that in the vicinity of buildings (<= 3 grid points). Even though no 
!--    building is within the numerical stencil, first-order scheme is used. 
!--    At fourth and fifth grid point the order of the horizontal advection scheme
!--    is successively upgraded.
!--    These extended degradation zones are used to avoid stationary numerical 
!--    oscillations, which are responsible for high concentration maxima that may
!--    appear under shear-free stable conditions. 
       CALL ws_init_flags_scalar(                                              &
                  bc_dirichlet_l  .OR.  bc_radiation_l  .OR.  decycle_chem_lr, &
                  bc_dirichlet_n  .OR.  bc_radiation_n  .OR.  decycle_chem_ns, &
                  bc_dirichlet_r  .OR.  bc_radiation_r  .OR.  decycle_chem_lr, &
                  bc_dirichlet_s  .OR.  bc_radiation_s  .OR.  decycle_chem_ns, &
                  cs_advc_flags_s, .TRUE. )
    ENDIF
!
!-- Initial concentration of profiles is prescribed by parameters cs_profile
!-- and cs_heights in the namelist &chemistry_parameters

    CALL chem_init_profiles
!    
!-- In case there is dynamic input file, create a list of names for chemistry
!-- initial input files. Also, initialize array that indicates whether the
!-- respective variable is on file or not. 
    IF ( input_pids_dynamic )  THEN    
       ALLOCATE( init_3d%var_names_chem(1:nspec) )
       ALLOCATE( init_3d%from_file_chem(1:nspec) )
       init_3d%from_file_chem(:) = .FALSE.
       
       DO  lsp = 1, nspec
          init_3d%var_names_chem(lsp) = init_3d%init_char // TRIM( chem_species(lsp)%name )
       ENDDO
    ENDIF
!
!-- Initialize model variables
    IF ( TRIM( initializing_actions ) /= 'read_restart_data'  .AND.            &
         TRIM( initializing_actions ) /= 'cyclic_fill' )  THEN
!
!--    First model run of a possible job queue.
!--    Initial profiles of the variables must be computed.
       IF ( INDEX( initializing_actions, 'set_1d-model_profiles' ) /= 0 )  THEN
!
!--       Transfer initial profiles to the arrays of the 3D model
          DO  lsp = 1, nspec
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO lpr_lev = 1, nz + 1
                      chem_species(lsp)%conc(lpr_lev,j,i) = chem_species(lsp)%conc_pr_init(lpr_lev)
                   ENDDO
                ENDDO
             ENDDO   
          ENDDO

       ELSEIF ( INDEX(initializing_actions, 'set_constant_profiles') /= 0 )    &
       THEN

          DO  lsp = 1, nspec 
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   chem_species(lsp)%conc(:,j,i) = chem_species(lsp)%conc_pr_init    
                ENDDO
             ENDDO
          ENDDO

       ENDIF
!
!--    If required, change the surface chem spcs at the start of the 3D run
       IF ( cs_surface_initial_change(1) /= 0.0_wp )  THEN            
          DO  lsp = 1, nspec 
             chem_species(lsp)%conc(nzb,:,:) = chem_species(lsp)%conc(nzb,:,:) +  &
                                               cs_surface_initial_change(lsp)
          ENDDO
       ENDIF 

    ENDIF
!
!-- Initiale old and new time levels. Note, this has to be done also in restart runs
    DO  lsp = 1, nvar
       chem_species(lsp)%tconc_m = 0.0_wp                      
       chem_species(lsp)%conc_p  = chem_species(lsp)%conc     
    ENDDO

    DO  lsp = 1, nphot
       phot_frequen(lsp)%name = phot_names(lsp)
!
!-- todo: remove or replace by "CALL message" mechanism (kanani)
!--       IF( myid == 0 )  THEN
!--          WRITE(6,'(a,i4,3x,a)')  'Photolysis: ',lsp,TRIM( phot_names(lsp) )
!--       ENDIF 
       phot_frequen(lsp)%freq(nzb:nzt+1,nysg:nyng,nxlg:nxrg)  =>  freq_1(:,:,:,lsp)
    ENDDO

!    CALL photolysis_init   ! probably also required for restart

    RETURN

 END SUBROUTINE chem_init_internal


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine defining initial vertical profiles of chemical species (given by
!> namelist parameters chem_profiles and chem_heights)  --> which should work
!> analogue to parameters u_profile, v_profile and uv_heights)
!------------------------------------------------------------------------------!
 SUBROUTINE chem_init_profiles              

    USE chem_modules

!
!-- Local variables
    INTEGER ::  lsp        !< running index for number of species in derived data type species_def
    INTEGER ::  lsp_usr     !< running index for number of species (user defined)  in cs_names, cs_profiles etc
    INTEGER ::  lpr_lev    !< running index for profile level for each chem spcs. 
    INTEGER ::  npr_lev    !< the next available profile lev
!
!-- Parameter "cs_profile" and "cs_heights" are used to prescribe user defined initial profiles
!-- and heights. If parameter "cs_profile" is not prescribed then initial surface values
!-- "cs_surface" are used as constant initial profiles for each species. If "cs_profile" and
!-- "cs_heights" are prescribed, their values will!override the constant profile given by
!-- "cs_surface".
!     IF ( TRIM( initializing_actions ) /= 'read_restart_data' )  THEN 
    lsp_usr = 1
    DO  WHILE ( TRIM( cs_name( lsp_usr ) ) /= 'novalue' )   !'novalue' is the default
       DO  lsp = 1, nspec                                !
!
!--       create initial profile (conc_pr_init) for each chemical species
          IF ( TRIM( chem_species(lsp)%name ) == TRIM( cs_name(lsp_usr) ) )  THEN   !
             IF ( cs_profile(lsp_usr,1) == 9999999.9_wp )  THEN
!
!--            set a vertically constant profile based on the surface conc (cs_surface(lsp_usr)) of each species
                DO lpr_lev = 0, nzt+1
                   chem_species(lsp)%conc_pr_init(lpr_lev) = cs_surface(lsp_usr)
                ENDDO
             ELSE
                IF ( cs_heights(1,1) /= 0.0_wp )  THEN
                   message_string = 'The surface value of cs_heights must be 0.0'
                   CALL message( 'chem_check_parameters', 'CM0434', 1, 2, 0, 6, 0 )
                ENDIF

                use_prescribed_profile_data = .TRUE.

                npr_lev = 1
!                chem_species(lsp)%conc_pr_init(0) = 0.0_wp
                DO  lpr_lev = 1, nz+1
                   IF ( npr_lev < 100 )  THEN
                      DO  WHILE ( cs_heights(lsp_usr, npr_lev+1) <= zu(lpr_lev) )
                         npr_lev = npr_lev + 1
                         IF ( npr_lev == 100 )  THEN
                            message_string = 'number of chem spcs exceeding the limit'
                            CALL message( 'chem_check_parameters', 'CM0435', 1, 2, 0, 6, 0 )               
                            EXIT
                         ENDIF
                      ENDDO
                   ENDIF
                   IF ( npr_lev < 100  .AND.  cs_heights(lsp_usr,npr_lev+1) /= 9999999.9_wp )  THEN
                      chem_species(lsp)%conc_pr_init(lpr_lev) = cs_profile(lsp_usr, npr_lev) +       &
                           ( zu(lpr_lev) - cs_heights(lsp_usr, npr_lev) ) /                          &
                           ( cs_heights(lsp_usr, (npr_lev + 1)) - cs_heights(lsp_usr, npr_lev ) ) *  &
                           ( cs_profile(lsp_usr, (npr_lev + 1)) - cs_profile(lsp_usr, npr_lev ) )
                   ELSE
                      chem_species(lsp)%conc_pr_init(lpr_lev) = cs_profile(lsp_usr, npr_lev)
                   ENDIF
                ENDDO
             ENDIF
!
!--       If a profile is prescribed explicity using cs_profiles and cs_heights, then  
!--       chem_species(lsp)%conc_pr_init is populated with the specific "lsp" based
!--       on the cs_profiles(lsp_usr,:)  and cs_heights(lsp_usr,:). 
          ENDIF

       ENDDO

       lsp_usr = lsp_usr + 1
    ENDDO
!     ENDIF

 END SUBROUTINE chem_init_profiles

 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine to integrate chemical species in the given chemical mechanism
!------------------------------------------------------------------------------!
 SUBROUTINE chem_integrate_ij( i, j )

    USE statistics,                                                          &
        ONLY:  weight_pres

    USE control_parameters,                                                  &
        ONLY:  dt_3d, intermediate_timestep_count, time_since_reference_point


    INTEGER,INTENT(IN)       :: i
    INTEGER,INTENT(IN)       :: j
!
!--   local variables
    INTEGER(iwp) ::  lsp                                                     !< running index for chem spcs.
    INTEGER(iwp) ::  lph                                                     !< running index for photolysis frequencies
    INTEGER, DIMENSION(20)    :: istatus
    REAL(kind=wp), DIMENSION(nzb+1:nzt,nspec)                :: tmp_conc
    REAL(kind=wp), DIMENSION(nzb+1:nzt)                      :: tmp_temp
    REAL(kind=wp), DIMENSION(nzb+1:nzt)                      :: tmp_qvap
    REAL(kind=wp), DIMENSION(nzb+1:nzt,nphot)                :: tmp_phot
    REAL(kind=wp), DIMENSION(nzb+1:nzt)                      :: tmp_fact
    REAL(kind=wp), DIMENSION(nzb+1:nzt)                      :: tmp_fact_i    !< conversion factor between 
                                                                              !<    molecules cm^{-3} and ppm

    INTEGER,DIMENSION(nzb+1:nzt)                            :: nacc          !< Number of accepted steps
    INTEGER,DIMENSION(nzb+1:nzt)                            :: nrej          !< Number of rejected steps

    REAL(wp)                         ::  conv                                !< conversion factor
    REAL(wp), PARAMETER              ::  ppm2fr  = 1.0e-6_wp                 !< Conversion factor ppm to fraction
    REAL(wp), PARAMETER              ::  fr2ppm  = 1.0e6_wp                  !< Conversion factor fraction to ppm
!    REAL(wp), PARAMETER              ::  xm_air  = 28.96_wp                  !< Mole mass of dry air
!    REAL(wp), PARAMETER              ::  xm_h2o  = 18.01528_wp               !< Mole mass of water vapor
    REAL(wp), PARAMETER              ::  t_std   = 273.15_wp                 !< standard pressure (Pa)
    REAL(wp), PARAMETER              ::  p_std   = 101325.0_wp               !< standard pressure (Pa)
    REAL(wp), PARAMETER              ::  vmolcm  = 22.414e3_wp               !< Mole volume (22.414 l) in cm^3
    REAL(wp), PARAMETER              ::  xna     = 6.022e23_wp               !< Avogadro number (molecules/mol)

    REAL(wp),DIMENSION(size(rcntrl)) :: rcntrl_local

    REAL(kind=wp)  :: dt_chem                                             
!
!-- Set chem_gasphase_on to .FALSE. if you want to skip computation of gas phase chemistry
    IF (chem_gasphase_on)  THEN
       nacc = 0
       nrej = 0

       tmp_temp(:) = pt(nzb+1:nzt,j,i) * exner(nzb+1:nzt)
!
!--    convert ppm to molecules/cm**3
!--    tmp_fact = 1.e-6_wp*6.022e23_wp/(22.414_wp*1000._wp) * 273.15_wp * 
!--               hyp(nzb+1:nzt)/( 101300.0_wp * tmp_temp )  
       conv = ppm2fr * xna / vmolcm
       tmp_fact(:) = conv * t_std * hyp(nzb+1:nzt) / (tmp_temp(:) * p_std)
       tmp_fact_i = 1.0_wp/tmp_fact

       IF ( humidity )  THEN
          IF ( bulk_cloud_model )  THEN
             tmp_qvap(:) = ( q(nzb+1:nzt,j,i) - ql(nzb+1:nzt,j,i) ) *  &
                             xm_air/xm_h2o * fr2ppm * tmp_fact(:)
          ELSE
             tmp_qvap(:) = q(nzb+1:nzt,j,i) * xm_air/xm_h2o * fr2ppm * tmp_fact(:)
          ENDIF
       ELSE
          tmp_qvap(:) = 0.01 * xm_air/xm_h2o * fr2ppm * tmp_fact(:)          !< Constant value for q if water vapor is not computed
       ENDIF

       DO  lsp = 1,nspec
          tmp_conc(:,lsp) = chem_species(lsp)%conc(nzb+1:nzt,j,i) * tmp_fact(:) 
       ENDDO

       DO lph = 1,nphot
          tmp_phot(:,lph) = phot_frequen(lph)%freq(nzb+1:nzt,j,i)               
       ENDDO
!
!--    Compute length of time step
       IF ( call_chem_at_all_substeps )  THEN
          dt_chem = dt_3d * weight_pres(intermediate_timestep_count)
       ELSE
          dt_chem = dt_3d
       ENDIF

       cs_time_step = dt_chem

       IF(maxval(rcntrl) > 0.0)  THEN    ! Only if rcntrl is set
          IF( time_since_reference_point <= 2*dt_3d)  THEN
             rcntrl_local = 0
          ELSE
             rcntrl_local = rcntrl
          ENDIF
       ELSE
          rcntrl_local = 0
       END IF

       CALL chem_gasphase_integrate ( dt_chem, tmp_conc, tmp_temp, tmp_qvap, tmp_fact, tmp_phot, &
            icntrl_i = icntrl, rcntrl_i = rcntrl_local, xnacc = nacc, xnrej = nrej, istatus=istatus )

       DO  lsp = 1,nspec
          chem_species(lsp)%conc (nzb+1:nzt,j,i) = tmp_conc(:,lsp) * tmp_fact_i(:)
       ENDDO


    ENDIF

    RETURN
 END SUBROUTINE chem_integrate_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine defining parin for &chemistry_parameters for chemistry model
!------------------------------------------------------------------------------!
 SUBROUTINE chem_parin

    USE chem_modules
    USE control_parameters

    USE pegrid
    USE statistics


    CHARACTER (LEN=80) ::  line                        !< dummy string that contains the current line of the parameter file

    REAL(wp), DIMENSION(nmaxfixsteps) ::   my_steps    !< List of fixed timesteps   my_step(1) = 0.0 automatic stepping
    INTEGER(iwp) ::  i                                 !< 
    INTEGER(iwp) ::  max_pr_cs_tmp                     !< 


    NAMELIST /chemistry_parameters/  bc_cs_b,                          &
         bc_cs_t,                          &
         call_chem_at_all_substeps,        &
         chem_debug0,                      &
         chem_debug1,                      &
         chem_debug2,                      &
         chem_gasphase_on,                 &
         chem_mechanism,                   &          
         cs_heights,                       &
         cs_name,                          &
         cs_profile,                       &
         cs_surface,                       &
         cs_surface_initial_change,        &
         cs_vertical_gradient_level,       &
         daytype_mdh,                      &
         decycle_chem_lr,                  &
         decycle_chem_ns,                  &            
         decycle_method,                   &
         deposition_dry,                   &
         emissions_anthropogenic,          &
         emiss_lod,                        &
         emiss_factor_main,                &
         emiss_factor_side,                &
         emiss_read_legacy_mode,           &
         icntrl,                           &
         main_street_id,                   &
         max_street_id,                    &
         mode_emis,                        &
         my_steps,                         &
         nesting_chem,                     &
         nesting_offline_chem,             &
         rcntrl,                           &
         side_street_id,                   &
         photolysis_scheme,                &
         wall_csflux,                      &
         cs_vertical_gradient,             &
         top_csflux,                       & 
         surface_csflux,                   &
         surface_csflux_name,              &
         time_fac_type
!
!-- analogue to chem_names(nspj) we could invent chem_surfaceflux(nspj) and chem_topflux(nspj)
!-- so this way we could prescribe a specific flux value for each species
    !>  chemistry_parameters for initial profiles
    !>  cs_names = 'O3', 'NO2', 'NO', ...   to set initial profiles)
    !>  cs_heights(1,:) = 0.0, 100.0, 500.0, 2000.0, .... (height levels where concs will be prescribed for O3)
    !>  cs_heights(2,:) = 0.0, 200.0, 400.0, 1000.0, .... (same for NO2 etc.) 
    !>  cs_profiles(1,:) = 10.0, 20.0, 20.0, 30.0, .....  (chem spcs conc at height lvls chem_heights(1,:)) etc.
    !>  If the respective concentration profile should be constant with height, then use "cs_surface( number of spcs)"
    !>  then write these cs_surface values to chem_species(lsp)%conc_pr_init(:)
!
!--   Read chem namelist   

    CHARACTER(LEN=8)    :: solver_type

    icntrl    = 0
    rcntrl    = 0.0_wp
    my_steps  = 0.0_wp
    photolysis_scheme = 'simple'
    atol = 1.0_wp
    rtol = 0.01_wp
!
!--   Try to find chemistry package
    REWIND ( 11 )
    line = ' '
    DO   WHILE ( INDEX( line, '&chemistry_parameters' ) == 0 )
       READ ( 11, '(A)', END=20 )  line
    ENDDO
    BACKSPACE ( 11 )
!
!--   Read chemistry namelist
    READ ( 11, chemistry_parameters, ERR = 10, END = 20 )      
!
!--   Enable chemistry model
    air_chemistry = .TRUE.                    
    GOTO 20

 10 BACKSPACE( 11 )
    READ( 11 , '(A)') line
    CALL parin_fail_message( 'chemistry_parameters', line )

 20 CONTINUE



!
!-- synchronize emiss_lod and mod_emis only if emissions_anthropogenic
!-- is activated in the namelist.  Otherwise their values are "don't care"
    IF ( emissions_anthropogenic )  THEN

!
!--    check for emission mode for chem species

       IF ( emiss_lod < 0 )  THEN   !- if LOD not defined in namelist
          IF ( ( mode_emis /= 'PARAMETERIZED'  )    .AND.      &
               ( mode_emis /= 'DEFAULT'        )    .AND.      &
               ( mode_emis /= 'PRE-PROCESSED'  ) )  THEN
             message_string = 'Incorrect mode_emiss  option select. Please check spelling'
             CALL message( 'chem_check_parameters', 'CM0436', 1, 2, 0, 6, 0 )
          ENDIF
       ELSE
          IF ( ( emiss_lod /= 0 )    .AND.         &
               ( emiss_lod /= 1 )    .AND.         &
               ( emiss_lod /= 2 ) )  THEN
             message_string = 'Invalid value for emiss_lod (0, 1, or 2)'
             CALL message( 'chem_check_parameters', 'CM0436', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF

!
! for reference (ecc)
!    IF ( (mode_emis /= 'PARAMETERIZED')  .AND. ( mode_emis /= 'DEFAULT' ) .AND. ( mode_emis /= 'PRE-PROCESSED'  ) )  THEN
!       message_string = 'Incorrect mode_emiss  option select. Please check spelling'
!       CALL message( 'chem_check_parameters', 'CM0436', 1, 2, 0, 6, 0 )
!    ENDIF

!
!-- conflict resolution for emiss_lod and mode_emis
!-- 1) if emiss_lod is defined, have mode_emis assume same setting as emiss_lod
!-- 2) if emiss_lod it not defined, have emiss_lod assuem same setting as mode_emis
!-- this check is in place to retain backward compatibility with salsa until the
!-- code is migrated completed to emiss_lod
!-- note that 

       IF  ( emiss_lod >= 0 ) THEN

          SELECT CASE  ( emiss_lod )
             CASE (0)  !- parameterized mode
                mode_emis = 'PARAMETERIZED'
             CASE (1)  !- default mode
                mode_emis = 'DEFAULT'
             CASE (2)  !- preprocessed mode
                mode_emis = 'PRE-PROCESSED'
          END SELECT
       
          message_string = 'Synchronizing mode_emis to defined emiss_lod'               //  &
                           CHAR(10)  //  '                    '                         //  &
                           'NOTE - mode_emis will be depreciated in future releases'    //  &
                           CHAR(10)  //  '                    '                         //  &
                           'please use emiss_lod to define emission mode'
 
          CALL message ( 'parin_chem', 'CM0463', 0, 0, 0, 6, 0 )

       ELSE ! if emiss_lod is not set

          SELECT CASE ( mode_emis )
             CASE ('PARAMETERIZED')
                emiss_lod = 0
             CASE ('DEFAULT')
                emiss_lod = 1
             CASE ('PRE-PROCESSED')
                emiss_lod = 2
          END SELECT

          message_string = 'emiss_lod undefined.  Using existing mod_emiss setting'     //  &
                           CHAR(10)  //  '                    '                         //  &
                           'NOTE - mode_emis will be depreciated in future releases'    //  &
                           CHAR(10)  //  '                    '                         //  &
                           'please use emiss_lod to define emission mode'

          CALL message ( 'parin_chem', 'CM0464', 0, 0, 0, 6, 0 )
       ENDIF

!
!-- (ECC) input check for emission read mode.
!--    legacy : business as usual (everything read / set up at start of run)
!--    new    : emission based on timestamp, and for lod2 data is loaded on an hourly basis

!
!-- (ECC) handler for emiss_read_legacy_mode
!-- * emiss_read_legacy_mode is defaulted to TRUE
!-- * if emiss_read_legacy_mode is TRUE and LOD is 0 or 1,
!--       force emission_read_legacy_mode to TRUE (not yet implemented)

       IF ( emiss_read_legacy_mode )  THEN       !< notify legacy read mode

          message_string = 'Legacy emission read mode activated'            // &
                           CHAR(10)  //  '                    '             // &
                           'All emissions data will be loaded '             // &
                           'prior to start of simulation'
          CALL message ( 'parin_chem', 'CM0465', 0, 0, 0, 6, 0 )

       ELSE                                     !< if new read mode selected 

          IF ( emiss_lod < 2 )  THEN            !< check LOD compatibility

             message_string = 'New emission read mode '                    // &
                              'currently unavailable for LODs 0 and 1.'    // &
                              CHAR(10)  //  '                    '         // &
                              'Reverting to legacy emission read mode'
             CALL message ( 'parin_chem', 'CM0466', 0, 0, 0, 6, 0 )

             emiss_read_legacy_mode = .TRUE.

          ELSE                                  !< notify new read mode

             message_string = 'New emission read mode activated'           // &
                              CHAR(10)  //  '                    '         // &
                              'LOD 2 emissions will be updated on-demand ' // &           
                              'according to indicated timestamps'
             CALL message ( 'parin_chem', 'CM0467', 0, 0, 0, 6, 0 )

          ENDIF

       ENDIF ! if emiss_read_legacy_mode
 
                           
    ENDIF  ! if emissions_anthropengic


    t_steps = my_steps          
!
!--    Determine the number of user-defined profiles and append them to the
!--    standard data output (data_output_pr)
    max_pr_cs_tmp = 0 
    i = 1

    DO  WHILE ( data_output_pr(i)  /= ' '  .AND.  i <= 100 )
       IF ( TRIM( data_output_pr(i)(1:3) ) == 'kc_' )  THEN
          max_pr_cs_tmp = max_pr_cs_tmp+1
       ENDIF
       i = i +1
    ENDDO

    IF ( max_pr_cs_tmp > 0 )  THEN
       cs_pr_namelist_found = .TRUE.
       max_pr_cs = max_pr_cs_tmp
    ENDIF

    !      Set Solver Type
    IF(icntrl(3) == 0)  THEN
       solver_type = 'rodas3'           !Default
    ELSE IF(icntrl(3) == 1)  THEN
       solver_type = 'ros2'
    ELSE IF(icntrl(3) == 2)  THEN
       solver_type = 'ros3'
    ELSE IF(icntrl(3) == 3)  THEN
       solver_type = 'ro4'
    ELSE IF(icntrl(3) == 4)  THEN
       solver_type = 'rodas3'
    ELSE IF(icntrl(3) == 5)  THEN
       solver_type = 'rodas4'
    ELSE IF(icntrl(3) == 6)  THEN
       solver_type = 'Rang3'
    ELSE
       message_string = 'illegal solver type'
       CALL message( 'chem_parin', 'PA0506', 1, 2, 0, 6, 0 )
    END IF

!
!--   todo: remove or replace by "CALL message" mechanism (kanani)
!       write(text,*) 'gas_phase chemistry: solver_type = ',TRIM( solver_type )
!kk    Has to be changed to right calling sequence
!        IF(myid == 0)  THEN
!           write(9,*) ' '
!           write(9,*) 'kpp setup '
!           write(9,*) ' '
!           write(9,*) '    gas_phase chemistry: solver_type = ',TRIM( solver_type )
!           write(9,*) ' '
!           write(9,*) '    Hstart  = ',rcntrl(3)
!           write(9,*) '    FacMin  = ',rcntrl(4)
!           write(9,*) '    FacMax  = ',rcntrl(5)
!           write(9,*) ' '
!           IF(vl_dim > 1)  THEN
!              write(9,*) '    Vector mode                   vektor length = ',vl_dim
!           ELSE
!              write(9,*) '    Scalar mode'
!           ENDIF
!           write(9,*) ' '
!        END IF

    RETURN

 END SUBROUTINE chem_parin


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE chem_actions( location )


    CHARACTER (LEN=*), INTENT(IN) ::  location !< call location string

    SELECT CASE ( location )

       CASE ( 'before_prognostic_equations' )
!
!--       Chemical reactions and deposition
          IF ( chem_gasphase_on ) THEN
!
!--          If required, calculate photolysis frequencies -
!--          UNFINISHED: Why not before the intermediate timestep loop?
             IF ( intermediate_timestep_count ==  1 )  THEN
                CALL photolysis_control
             ENDIF

          ENDIF

       CASE DEFAULT
          CONTINUE

    END SELECT

    END SUBROUTINE chem_actions


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid points i,j
!------------------------------------------------------------------------------!

    SUBROUTINE chem_actions_ij( i, j, location )


    INTEGER(iwp),      INTENT(IN) ::  i         !< grid index in x-direction
    INTEGER(iwp),      INTENT(IN) ::  j         !< grid index in y-direction
    CHARACTER (LEN=*), INTENT(IN) ::  location  !< call location string
    INTEGER(iwp)  ::  dummy  !< call location string

    IF ( air_chemistry    )   dummy = i + j

    SELECT CASE ( location )

       CASE DEFAULT
          CONTINUE

    END SELECT


    END SUBROUTINE chem_actions_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE chem_non_advective_processes()


      INTEGER(iwp) ::  i  !<
      INTEGER(iwp) ::  j  !<

      !
      !-- Calculation of chemical reactions and deposition.


      IF ( intermediate_timestep_count == 1 .OR. call_chem_at_all_substeps )  THEN

         IF ( chem_gasphase_on ) THEN
            CALL cpu_log( log_point_s(19), 'chem.reactions', 'start' )
            !$OMP PARALLEL PRIVATE (i,j)
            !$OMP DO schedule(static,1)
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  CALL chem_integrate( i, j )
               ENDDO
            ENDDO
            !$OMP END PARALLEL
            CALL cpu_log( log_point_s(19), 'chem.reactions', 'stop' )
         ENDIF

         IF ( deposition_dry )  THEN
            CALL cpu_log( log_point_s(24), 'chem.deposition', 'start' )
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  CALL chem_depo( i, j )
               ENDDO
            ENDDO
            CALL cpu_log( log_point_s(24), 'chem.deposition', 'stop' )
         ENDIF

      ENDIF



    END SUBROUTINE chem_non_advective_processes


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid points i,j
!------------------------------------------------------------------------------!
 SUBROUTINE chem_non_advective_processes_ij( i, j )


   INTEGER(iwp), INTENT(IN) ::  i  !< grid index in x-direction
   INTEGER(iwp), INTENT(IN) ::  j  !< grid index in y-direction

!
!-- Calculation of chemical reactions and deposition.


   IF ( intermediate_timestep_count == 1 .OR. call_chem_at_all_substeps )  THEN

      IF ( chem_gasphase_on ) THEN
         CALL cpu_log( log_point_s(19), 'chem.reactions', 'start' )
         CALL chem_integrate( i, j )
         CALL cpu_log( log_point_s(19), 'chem.reactions', 'stop' )
      ENDIF

      IF ( deposition_dry )  THEN
         CALL cpu_log( log_point_s(24), 'chem.deposition', 'start' )
         CALL chem_depo( i, j )
         CALL cpu_log( log_point_s(24), 'chem.deposition', 'stop' )
      ENDIF

   ENDIF



 END SUBROUTINE chem_non_advective_processes_ij
 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> routine for exchange horiz of chemical quantities  
!------------------------------------------------------------------------------! 
 SUBROUTINE chem_exchange_horiz_bounds
  
    USE exchange_horiz_mod,                                                    &
        ONLY:  exchange_horiz

   INTEGER(iwp) ::  lsp       !<
   INTEGER(iwp) ::  lsp_usr   !<
 
!
!--    Loop over chemical species       
       CALL cpu_log( log_point_s(84), 'chem.exch-horiz', 'start' )
       DO  lsp = 1, nvar
          CALL exchange_horiz( chem_species(lsp)%conc, nbgp )    
          lsp_usr = 1  
          DO WHILE ( TRIM( cs_name( lsp_usr ) ) /= 'novalue' )
             IF ( TRIM(chem_species(lsp)%name) == TRIM(cs_name(lsp_usr)) )  THEN
!
!--             As chem_exchange_horiz_bounds is called at the beginning
!--             of prognostic_equations, boundary conditions are set on 
!--             %conc. 
                CALL chem_boundary_conds( chem_species(lsp)%conc,              &
                                          chem_species(lsp)%conc_pr_init ) 
             
             ENDIF
             lsp_usr = lsp_usr + 1
          ENDDO


       ENDDO
       CALL cpu_log( log_point_s(84), 'chem.exch-horiz', 'stop' )


 END SUBROUTINE chem_exchange_horiz_bounds

 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine calculating prognostic equations for chemical species 
!> (vector-optimized).
!> Routine is called separately for each chemical species over a loop from
!> prognostic_equations.
!------------------------------------------------------------------------------!
 SUBROUTINE chem_prognostic_equations()


    INTEGER ::  i   !< running index
    INTEGER ::  j   !< running index
    INTEGER ::  k   !< running index

    INTEGER(iwp) ::  ilsp   !<


    CALL cpu_log( log_point_s(25), 'chem.advec+diff+prog', 'start' )

    DO  ilsp = 1, nvar
!
!--    Tendency terms for chemical species
       tend = 0.0_wp
!
!--    Advection terms
       IF ( timestep_scheme(1:5) == 'runge' )  THEN
          IF ( ws_scheme_sca )  THEN
             CALL advec_s_ws( cs_advc_flags_s, chem_species(ilsp)%conc, 'kc',                      &
                              bc_dirichlet_l  .OR.  bc_radiation_l  .OR.  decycle_chem_lr,         &
                              bc_dirichlet_n  .OR.  bc_radiation_n  .OR.  decycle_chem_ns,         &
                              bc_dirichlet_r  .OR.  bc_radiation_r  .OR.  decycle_chem_lr,         &
                              bc_dirichlet_s  .OR.  bc_radiation_s  .OR.  decycle_chem_ns )
          ELSE
             CALL advec_s_pw( chem_species(ilsp)%conc )
          ENDIF
       ELSE
          CALL advec_s_up( chem_species(ilsp)%conc )
       ENDIF
!
!--    Diffusion terms  (the last three arguments are zero)
       CALL diffusion_s( chem_species(ilsp)%conc,                                                  &
            surf_def_h(0)%cssws(ilsp,:),                                                           &
            surf_def_h(1)%cssws(ilsp,:),                                                           &
            surf_def_h(2)%cssws(ilsp,:),                                                           &
            surf_lsm_h%cssws(ilsp,:),                                                              &
            surf_usm_h%cssws(ilsp,:),                                                              &
            surf_def_v(0)%cssws(ilsp,:),                                                           &
            surf_def_v(1)%cssws(ilsp,:),                                                           &
            surf_def_v(2)%cssws(ilsp,:),                                                           &
            surf_def_v(3)%cssws(ilsp,:),                                                           &
            surf_lsm_v(0)%cssws(ilsp,:),                                                           &
            surf_lsm_v(1)%cssws(ilsp,:),                                                           &
            surf_lsm_v(2)%cssws(ilsp,:),                                                           &
            surf_lsm_v(3)%cssws(ilsp,:),                                                           &
            surf_usm_v(0)%cssws(ilsp,:),                                                           &
            surf_usm_v(1)%cssws(ilsp,:),                                                           &
            surf_usm_v(2)%cssws(ilsp,:),                                                           &
            surf_usm_v(3)%cssws(ilsp,:) )
!
!--    Prognostic equation for chemical species
       DO  i = nxl, nxr
          DO  j = nys, nyn
             !following directive is required to vectorize on Intel19
             !DIR$ IVDEP
             DO  k = nzb+1, nzt
                chem_species(ilsp)%conc_p(k,j,i) =   chem_species(ilsp)%conc(k,j,i)                &
                     + ( dt_3d  *                                                                  &
                     (   tsc(2) * tend(k,j,i)                                                      &
                     + tsc(3) * chem_species(ilsp)%tconc_m(k,j,i)                                  &
                     )                                                                             &
                     - tsc(5) * rdf_sc(k)                                                          &
                     * ( chem_species(ilsp)%conc(k,j,i) - chem_species(ilsp)%conc_pr_init(k) )     &
                     )                                                                             &
                     * MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )

                IF ( chem_species(ilsp)%conc_p(k,j,i) < 0.0_wp )  THEN
                   chem_species(ilsp)%conc_p(k,j,i) = 0.1_wp * chem_species(ilsp)%conc(k,j,i)
                ENDIF
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
                      chem_species(ilsp)%tconc_m(k,j,i) = tend(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSEIF ( intermediate_timestep_count < &
               intermediate_timestep_count_max )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt
                      chem_species(ilsp)%tconc_m(k,j,i) = - 9.5625_wp * tend(k,j,i)                &
                           + 5.3125_wp * chem_species(ilsp)%tconc_m(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDIF

    ENDDO

    CALL cpu_log( log_point_s(25), 'chem.advec+diff+prog', 'stop' )

 END SUBROUTINE chem_prognostic_equations


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine calculating prognostic equations for chemical species 
!> (cache-optimized).
!> Routine is called separately for each chemical species over a loop from
!> prognostic_equations.
!------------------------------------------------------------------------------!
 SUBROUTINE chem_prognostic_equations_ij( i, j, i_omp_start, tn )


    INTEGER(iwp),INTENT(IN) :: i, j, i_omp_start, tn
    INTEGER(iwp) :: ilsp
!
!-- local variables

    INTEGER :: k

    DO  ilsp = 1, nvar
!
!--    Tendency-terms for chem spcs.
       tend(:,j,i) = 0.0_wp
!
!--    Advection terms
       IF ( timestep_scheme(1:5) == 'runge' )  THEN
          IF ( ws_scheme_sca )  THEN
             CALL advec_s_ws( cs_advc_flags_s,                                                     &
                              i,                                                                   &
                              j,                                                                   &
                              chem_species(ilsp)%conc,                                             &
                              'kc',                                                                &
                              chem_species(ilsp)%flux_s_cs,                                        &
                              chem_species(ilsp)%diss_s_cs,                                        &
                              chem_species(ilsp)%flux_l_cs,                                        &
                              chem_species(ilsp)%diss_l_cs,                                        &
                              i_omp_start,                                                         &
                              tn,                                                                  &
                              bc_dirichlet_l  .OR.  bc_radiation_l  .OR.  decycle_chem_lr,         &
                              bc_dirichlet_n  .OR.  bc_radiation_n  .OR.  decycle_chem_ns,         &
                              bc_dirichlet_r  .OR.  bc_radiation_r  .OR.  decycle_chem_lr,         &
                              bc_dirichlet_s  .OR.  bc_radiation_s  .OR.  decycle_chem_ns,         &
                              monotonic_limiter_z )
          ELSE
             CALL advec_s_pw( i, j, chem_species(ilsp)%conc )
          ENDIF
       ELSE
          CALL advec_s_up( i, j, chem_species(ilsp)%conc )
       ENDIF
!
!--    Diffusion terms (the last three arguments are zero)

       CALL diffusion_s( i, j, chem_species(ilsp)%conc,                                            &
            surf_def_h(0)%cssws(ilsp,:), surf_def_h(1)%cssws(ilsp,:),                              &
            surf_def_h(2)%cssws(ilsp,:),                                                           &
            surf_lsm_h%cssws(ilsp,:), surf_usm_h%cssws(ilsp,:),                                    &
            surf_def_v(0)%cssws(ilsp,:), surf_def_v(1)%cssws(ilsp,:),                              &
            surf_def_v(2)%cssws(ilsp,:), surf_def_v(3)%cssws(ilsp,:),                              &
            surf_lsm_v(0)%cssws(ilsp,:), surf_lsm_v(1)%cssws(ilsp,:),                              &
            surf_lsm_v(2)%cssws(ilsp,:), surf_lsm_v(3)%cssws(ilsp,:),                              &
            surf_usm_v(0)%cssws(ilsp,:), surf_usm_v(1)%cssws(ilsp,:),                              &
            surf_usm_v(2)%cssws(ilsp,:), surf_usm_v(3)%cssws(ilsp,:) )
!
!--    Prognostic equation for chem spcs
       DO  k = nzb+1, nzt
          chem_species(ilsp)%conc_p(k,j,i) = chem_species(ilsp)%conc(k,j,i) + ( dt_3d  *           &
               ( tsc(2) * tend(k,j,i) +                                                            &
               tsc(3) * chem_species(ilsp)%tconc_m(k,j,i) )                                        &
               - tsc(5) * rdf_sc(k)                                                                &
               * ( chem_species(ilsp)%conc(k,j,i) - chem_species(ilsp)%conc_pr_init(k) )           &
               )                                                                                   &
               * MERGE( 1.0_wp, 0.0_wp,                                                            &
               BTEST( wall_flags_total_0(k,j,i), 0 )                                               &
               )

          IF ( chem_species(ilsp)%conc_p(k,j,i) < 0.0_wp )  THEN                
             chem_species(ilsp)%conc_p(k,j,i) = 0.1_wp * chem_species(ilsp)%conc(k,j,i)    !FKS6                
          ENDIF
       ENDDO
!
!--    Calculate tendencies for the next Runge-Kutta step
       IF ( timestep_scheme(1:5) == 'runge' )  THEN
          IF ( intermediate_timestep_count == 1 )  THEN
             DO  k = nzb+1, nzt
                chem_species(ilsp)%tconc_m(k,j,i) = tend(k,j,i)
             ENDDO
          ELSEIF ( intermediate_timestep_count < &
               intermediate_timestep_count_max )  THEN
             DO  k = nzb+1, nzt
                chem_species(ilsp)%tconc_m(k,j,i) = -9.5625_wp * tend(k,j,i) +                     &
                     5.3125_wp * chem_species(ilsp)%tconc_m(k,j,i)
             ENDDO
          ENDIF
       ENDIF

    ENDDO

 END SUBROUTINE chem_prognostic_equations_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine to read restart data of chemical species
!------------------------------------------------------------------------------!
 SUBROUTINE chem_rrd_local( k, nxlf, nxlc, nxl_on_file, nxrf, nxrc,             &
                            nxr_on_file, nynf, nync, nyn_on_file, nysf, nysc,   &
                            nys_on_file, tmp_3d, found )

    USE control_parameters


    CHARACTER (LEN=20) :: spc_name_av !<   

    INTEGER(iwp) ::  lsp             !<
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

    LOGICAL, INTENT(OUT) :: found 

    REAL(wp), DIMENSION(nzb:nzt+1,nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) :: tmp_3d   !< 3D array to temp store data


    found = .FALSE.  


    IF ( ALLOCATED(chem_species) )  THEN

       DO  lsp = 1, nspec

          !< for time-averaged chemical conc.
          spc_name_av  =  TRIM( chem_species(lsp)%name )//'_av'

          IF ( restart_string(1:length) == TRIM( chem_species(lsp)%name) )    &
             THEN
             !< read data into tmp_3d
             IF ( k == 1 )  READ ( 13 )  tmp_3d  
             !< fill ..%conc in the restart run    
             chem_species(lsp)%conc(:,nysc-nbgp:nync+nbgp,                    &
                  nxlc-nbgp:nxrc+nbgp) =                  & 
                  tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             found = .TRUE.
          ELSEIF (restart_string(1:length) == spc_name_av )  THEN
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             chem_species(lsp)%conc_av(:,nysc-nbgp:nync+nbgp,                 &
                  nxlc-nbgp:nxrc+nbgp) =               &
                  tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             found = .TRUE.
          ENDIF

       ENDDO

    ENDIF


 END SUBROUTINE chem_rrd_local


!-------------------------------------------------------------------------------!
!> Description:
!> Calculation of horizontally averaged profiles
!> This routine is called for every statistic region (sr) defined by the user,
!> but at least for the region "total domain" (sr=0).
!> quantities.
!-------------------------------------------------------------------------------!
 SUBROUTINE chem_statistics( mode, sr, tn )


    USE arrays_3d

    USE statistics


    CHARACTER (LEN=*) ::  mode   !< 

    INTEGER(iwp) ::  i    !< running index on x-axis
    INTEGER(iwp) ::  j    !< running index on y-axis
    INTEGER(iwp) ::  k    !< vertical index counter
    INTEGER(iwp) ::  sr   !< statistical region
    INTEGER(iwp) ::  tn   !< thread number 
    INTEGER(iwp) ::  lpr  !< running index chem spcs

    IF ( mode == 'profiles' )  THEN
       !
!
!--    Sample on how to calculate horizontally averaged profiles of user-
!--    defined quantities. Each quantity is identified by the index
!--    "pr_palm+#" where "#" is an integer starting from 1. These
!--    user-profile-numbers must also be assigned to the respective strings
!--    given by data_output_pr_cs in routine user_check_data_output_pr.
!--    hom(:,:,:,:) =  dim-1 = vertical level, dim-2= 1: met-species,2:zu/zw, dim-3 = quantity( e.g.
!--                     w*pt*), dim-4 = statistical region.

!$OMP DO
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb, nzt+1
                DO lpr = 1, cs_pr_count

                   sums_l(k,pr_palm+max_pr_user+lpr,tn) = sums_l(k,pr_palm+max_pr_user+lpr,tn) +    &
                        chem_species(cs_pr_index(lpr))%conc(k,j,i) *       &
                        rmask(j,i,sr)  *                                   &
                        MERGE( 1.0_wp, 0.0_wp,                             &
                        BTEST( wall_flags_total_0(k,j,i), 22 ) )
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ELSEIF ( mode == 'time_series' )  THEN
!      @todo
    ENDIF

 END SUBROUTINE chem_statistics


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine for swapping of timelevels for chemical species
!> called out from subroutine swap_timelevel
!------------------------------------------------------------------------------!


 SUBROUTINE chem_swap_timelevel( level )


    INTEGER(iwp), INTENT(IN) ::  level
!
!-- local variables
    INTEGER(iwp)             ::  lsp


    IF ( level == 0 )  THEN
       DO  lsp=1, nvar                                       
          chem_species(lsp)%conc(nzb:nzt+1,nysg:nyng,nxlg:nxrg)    => spec_conc_1(:,:,:,lsp)
          chem_species(lsp)%conc_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg)  => spec_conc_2(:,:,:,lsp)
       ENDDO
    ELSE
       DO  lsp=1, nvar                                        
          chem_species(lsp)%conc(nzb:nzt+1,nysg:nyng,nxlg:nxrg)    => spec_conc_2(:,:,:,lsp)
          chem_species(lsp)%conc_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg)  => spec_conc_1(:,:,:,lsp)
       ENDDO
    ENDIF

    RETURN
 END SUBROUTINE chem_swap_timelevel


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine to write restart data for chemistry model
!------------------------------------------------------------------------------!
 SUBROUTINE chem_wrd_local


    INTEGER(iwp) ::  lsp  !< running index for chem spcs. 

    DO  lsp = 1, nspec
       CALL wrd_write_string( TRIM( chem_species(lsp)%name ) )
       WRITE ( 14 )  chem_species(lsp)%conc
       CALL wrd_write_string( TRIM( chem_species(lsp)%name )//'_av' )
       WRITE ( 14 )  chem_species(lsp)%conc_av
    ENDDO

 END SUBROUTINE chem_wrd_local


!-------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine to calculate the deposition of gases and PMs. For now deposition
!> only takes place on lsm and usm horizontal surfaces. Default surfaces are NOT
!> considered. The deposition of particles is derived following Zhang et al.,
!> 2001, gases are deposited using the DEPAC module (van Zanten et al., 2010).
!>     
!> @TODO: Consider deposition on vertical surfaces    
!> @TODO: Consider overlaying horizontal surfaces
!> @TODO: Consider resolved vegetation    
!> @TODO: Check error messages
!-------------------------------------------------------------------------------!
 SUBROUTINE chem_depo( i, j )

    USE control_parameters,                                                 &   
         ONLY:  dt_3d, intermediate_timestep_count, latitude,               &
                time_since_reference_point

    USE arrays_3d,                                                          &
         ONLY:  dzw, rho_air_zw

    USE palm_date_time_mod,                                                 &
         ONLY:  get_date_time

    USE surface_mod,                                                        &
         ONLY:  ind_pav_green, ind_veg_wall, ind_wat_win, surf_lsm_h,        &
         surf_type, surf_usm_h

    USE radiation_model_mod,                                                &
         ONLY:  cos_zenith


    INTEGER(iwp) ::  day_of_year                   !< current day of the year
    INTEGER(iwp), INTENT(IN) ::  i
    INTEGER(iwp), INTENT(IN) ::  j
    INTEGER(iwp) ::  k                             !< matching k to surface m at i,j
    INTEGER(iwp) ::  lsp                           !< running index for chem spcs.
    INTEGER(iwp) ::  luv_palm                      !< index of PALM LSM vegetation_type at current surface element
    INTEGER(iwp) ::  lup_palm                      !< index of PALM LSM pavement_type at current surface element
    INTEGER(iwp) ::  luw_palm                      !< index of PALM LSM water_type at current surface element
    INTEGER(iwp) ::  luu_palm                      !< index of PALM USM walls/roofs at current surface element
    INTEGER(iwp) ::  lug_palm                      !< index of PALM USM green walls/roofs at current surface element
    INTEGER(iwp) ::  lud_palm                      !< index of PALM USM windows at current surface element
    INTEGER(iwp) ::  luv_dep                       !< matching DEPAC LU to luv_palm
    INTEGER(iwp) ::  lup_dep                       !< matching DEPAC LU to lup_palm
    INTEGER(iwp) ::  luw_dep                       !< matching DEPAC LU to luw_palm
    INTEGER(iwp) ::  luu_dep                       !< matching DEPAC LU to luu_palm
    INTEGER(iwp) ::  lug_dep                       !< matching DEPAC LU to lug_palm
    INTEGER(iwp) ::  lud_dep                       !< matching DEPAC LU to lud_palm
    INTEGER(iwp) ::  m                             !< index for horizontal surfaces

    INTEGER(iwp) ::  pspec                         !< running index
    INTEGER(iwp) ::  i_pspec                       !< index for matching depac gas component
!
!-- Vegetation                                               !< Assign PALM classes to DEPAC land use classes
    INTEGER(iwp) ::  ind_luv_user = 0                        !<  ERROR as no class given in PALM 
    INTEGER(iwp) ::  ind_luv_b_soil = 1                      !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_luv_mixed_crops = 2                 !<  assigned to ilu_arable
    INTEGER(iwp) ::  ind_luv_s_grass = 3                     !<  assigned to ilu_grass
    INTEGER(iwp) ::  ind_luv_ev_needle_trees = 4             !<  assigned to ilu_coniferous_forest
    INTEGER(iwp) ::  ind_luv_de_needle_trees = 5             !<  assigned to ilu_coniferous_forest
    INTEGER(iwp) ::  ind_luv_ev_broad_trees = 6              !<  assigned to ilu_tropical_forest
    INTEGER(iwp) ::  ind_luv_de_broad_trees = 7              !<  assigned to ilu_deciduous_forest
    INTEGER(iwp) ::  ind_luv_t_grass = 8                     !<  assigned to ilu_grass
    INTEGER(iwp) ::  ind_luv_desert = 9                      !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_luv_tundra = 10                     !<  assigned to ilu_other
    INTEGER(iwp) ::  ind_luv_irr_crops = 11                  !<  assigned to ilu_arable
    INTEGER(iwp) ::  ind_luv_semidesert = 12                 !<  assigned to ilu_other
    INTEGER(iwp) ::  ind_luv_ice = 13                        !<  assigned to ilu_ice
    INTEGER(iwp) ::  ind_luv_marsh = 14                      !<  assigned to ilu_other
    INTEGER(iwp) ::  ind_luv_ev_shrubs = 15                  !<  assigned to ilu_mediterrean_scrub
    INTEGER(iwp) ::  ind_luv_de_shrubs = 16                  !<  assigned to ilu_mediterrean_scrub
    INTEGER(iwp) ::  ind_luv_mixed_forest = 17               !<  assigned to ilu_coniferous_forest (ave(decid+conif))
    INTEGER(iwp) ::  ind_luv_intrup_forest = 18              !<  assigned to ilu_other (ave(other+decid))
!
!-- Water
    INTEGER(iwp) ::  ind_luw_user = 0                        !<  ERROR as no class given in PALM  
    INTEGER(iwp) ::  ind_luw_lake = 1                        !<  assigned to ilu_water_inland
    INTEGER(iwp) ::  ind_luw_river = 2                       !<  assigned to ilu_water_inland
    INTEGER(iwp) ::  ind_luw_ocean = 3                       !<  assigned to ilu_water_sea
    INTEGER(iwp) ::  ind_luw_pond = 4                        !<  assigned to ilu_water_inland
    INTEGER(iwp) ::  ind_luw_fountain = 5                    !<  assigned to ilu_water_inland 
!
!-- Pavement
    INTEGER(iwp) ::  ind_lup_user = 0                        !<  ERROR as no class given in PALM
    INTEGER(iwp) ::  ind_lup_asph_conc = 1                   !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_asph = 2                        !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_conc = 3                        !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_sett = 4                        !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_pav_stones = 5                  !<  assigned to ilu_desert 
    INTEGER(iwp) ::  ind_lup_cobblest = 6                    !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_metal = 7                       !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_wood = 8                        !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_gravel = 9                      !<  assigned to ilu_desert 
    INTEGER(iwp) ::  ind_lup_f_gravel = 10                   !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_pebblest = 11                   !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_woodchips = 12                  !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_tartan = 13                     !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_art_turf = 14                   !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_clay = 15                       !<  assigned to ilu_desert
!
!-- Particle parameters according to the respective aerosol classes (PM25, PM10)
    INTEGER(iwp) ::  ind_p_size = 1     !< index for partsize in particle_pars
    INTEGER(iwp) ::  ind_p_dens = 2     !< index for rhopart in particle_pars
    INTEGER(iwp) ::  ind_p_slip = 3     !< index for slipcor in particle_pars

    INTEGER(iwp) ::  part_type          !< index for particle type (PM10 or PM25) in particle_pars

    INTEGER(iwp) ::  nwet               !< wetness indicator dor DEPAC; nwet=0 -> dry; nwet=1 -> wet; nwet=9 -> snow

    REAL(wp) ::  dt_chem                !< length of chem time step
    REAL(wp) ::  dh                     !< vertical grid size 
    REAL(wp) ::  inv_dh                 !< inverse of vertical grid size 
    REAL(wp) ::  dt_dh                  !< dt_chem/dh

    REAL(wp) ::  dens              !< density at layer k at i,j  
    REAL(wp) ::  r_aero_surf       !< aerodynamic resistance (s/m) at current surface element
    REAL(wp) ::  ustar_surf        !< ustar at current surface element
    REAL(wp) ::  z0h_surf          !< roughness length for heat at current surface element 
    REAL(wp) ::  solar_rad         !< solar radiation, direct and diffuse, at current surface element 
    REAL(wp) ::  ppm2ugm3          !< conversion factor from ppm to ug/m3
    REAL(wp) ::  rh_surf           !< relative humidity at current surface element
    REAL(wp) ::  lai               !< leaf area index at current surface element
    REAL(wp) ::  sai               !< surface area index at current surface element assumed to be lai + 1

    REAL(wp) ::  slinnfac        
    REAL(wp) ::  visc              !< Viscosity
    REAL(wp) ::  vs                !< Sedimentation velocity
    REAL(wp) ::  vd_lu             !< deposition velocity (m/s)
    REAL(wp) ::  rs                !< Sedimentaion resistance (s/m)
    REAL(wp) ::  rb                !< quasi-laminar boundary layer resistance (s/m)
    REAL(wp) ::  rc_tot            !< total canopy resistance (s/m)

    REAL(wp) ::  conc_ijk_ugm3     !< concentration at i, j, k in ug/m3
    REAL(wp) ::  diffusivity       !< diffusivity


    REAL(wp), DIMENSION(nspec) ::  bud_luv      !< budget for LSM vegetation type at current surface element
    REAL(wp), DIMENSION(nspec) ::  bud_lup      !< budget for LSM pavement type at current surface element
    REAL(wp), DIMENSION(nspec) ::  bud_luw      !< budget for LSM water type at current surface element
    REAL(wp), DIMENSION(nspec) ::  bud_luu      !< budget for USM walls/roofs at current surface element
    REAL(wp), DIMENSION(nspec) ::  bud_lug      !< budget for USM green surfaces at current surface element
    REAL(wp), DIMENSION(nspec) ::  bud_lud      !< budget for USM windows at current surface element
    REAL(wp), DIMENSION(nspec) ::  bud          !< overall budget at current surface element
    REAL(wp), DIMENSION(nspec) ::  conc_ijk     !< concentration at i,j,k
    REAL(wp), DIMENSION(nspec) ::  ccomp_tot    !< total compensation point (ug/m3), for now kept to zero for all species!
                                                  

    REAL(wp) ::  temp_tmp       !< temperatur at i,j,k
    REAL(wp) ::  ts             !< surface temperatur in degrees celsius
    REAL(wp) ::  qv_tmp         !< surface mixing ratio at current surface element
!
!-- Particle parameters (PM10 (1), PM25 (2))
!-- partsize (diameter in m), rhopart (density in kg/m3), slipcor
!-- (slip correction factor dimensionless, Seinfeld and Pandis 2006, Table 9.3)
    REAL(wp), DIMENSION(1:3,1:2), PARAMETER ::  particle_pars = RESHAPE( (/ &
         8.0e-6_wp, 1.14e3_wp, 1.016_wp, &  !<  1
         0.7e-6_wp, 1.14e3_wp, 1.082_wp &   !<  2
         /), (/ 3, 2 /) )

    LOGICAL ::  match_lsm     !< flag indicating natural-type surface
    LOGICAL ::  match_usm     !< flag indicating urban-type surface
!
!-- List of names of possible tracers
    CHARACTER(LEN=*), PARAMETER ::  pspecnames(nposp) = (/ &
         'NO2           ', &    !< NO2
         'NO            ', &    !< NO
         'O3            ', &    !< O3
         'CO            ', &    !< CO
         'form          ', &    !< FORM
         'ald           ', &    !< ALD
         'pan           ', &    !< PAN
         'mgly          ', &    !< MGLY
         'par           ', &    !< PAR
         'ole           ', &    !< OLE
         'eth           ', &    !< ETH
         'tol           ', &    !< TOL
         'cres          ', &    !< CRES
         'xyl           ', &    !< XYL
         'SO4a_f        ', &    !< SO4a_f
         'SO2           ', &    !< SO2
         'HNO2          ', &    !< HNO2
         'CH4           ', &    !< CH4
         'NH3           ', &    !< NH3
         'NO3           ', &    !< NO3
         'OH            ', &    !< OH
         'HO2           ', &    !< HO2
         'N2O5          ', &    !< N2O5
         'SO4a_c        ', &    !< SO4a_c
         'NH4a_f        ', &    !< NH4a_f
         'NO3a_f        ', &    !< NO3a_f
         'NO3a_c        ', &    !< NO3a_c
         'C2O3          ', &    !< C2O3
         'XO2           ', &    !< XO2
         'XO2N          ', &    !< XO2N
         'cro           ', &    !< CRO
         'HNO3          ', &    !< HNO3
         'H2O2          ', &    !< H2O2
         'iso           ', &    !< ISO
         'ispd          ', &    !< ISPD
         'to2           ', &    !< TO2
         'open          ', &    !< OPEN
         'terp          ', &    !< TERP
         'ec_f          ', &    !< EC_f
         'ec_c          ', &    !< EC_c
         'pom_f         ', &    !< POM_f
         'pom_c         ', &    !< POM_c
         'ppm_f         ', &    !< PPM_f
         'ppm_c         ', &    !< PPM_c
         'na_ff         ', &    !< Na_ff
         'na_f          ', &    !< Na_f
         'na_c          ', &    !< Na_c
         'na_cc         ', &    !< Na_cc
         'na_ccc        ', &    !< Na_ccc
         'dust_ff       ', &    !< dust_ff
         'dust_f        ', &    !< dust_f
         'dust_c        ', &    !< dust_c
         'dust_cc       ', &    !< dust_cc
         'dust_ccc      ', &    !< dust_ccc
         'tpm10         ', &    !< tpm10
         'tpm25         ', &    !< tpm25
         'tss           ', &    !< tss
         'tdust         ', &    !< tdust
         'tc            ', &    !< tc
         'tcg           ', &    !< tcg
         'tsoa          ', &    !< tsoa
         'tnmvoc        ', &    !< tnmvoc
         'SOxa          ', &    !< SOxa
         'NOya          ', &    !< NOya
         'NHxa          ', &    !< NHxa
         'NO2_obs       ', &    !< NO2_obs
         'tpm10_biascorr', &    !< tpm10_biascorr
         'tpm25_biascorr', &    !< tpm25_biascorr
         'O3_biascorr   ' /)    !< o3_biascorr
!
!-- tracer mole mass:
    REAL(wp), PARAMETER ::  specmolm(nposp) = (/ &
         xm_O * 2 + xm_N, &                         !< NO2
         xm_O + xm_N, &                             !< NO
         xm_O * 3, &                                !< O3
         xm_C + xm_O, &                             !< CO
         xm_H * 2 + xm_C + xm_O, &                  !< FORM
         xm_H * 3 + xm_C * 2 + xm_O, &              !< ALD
         xm_H * 3 + xm_C * 2 + xm_O * 5 + xm_N, &   !< PAN
         xm_H * 4 + xm_C * 3 + xm_O * 2, &          !< MGLY
         xm_H * 3 + xm_C, &                         !< PAR
         xm_H * 3 + xm_C * 2, &                     !< OLE
         xm_H * 4 + xm_C * 2, &                     !< ETH
         xm_H * 8 + xm_C * 7, &                     !< TOL
         xm_H * 8 + xm_C * 7 + xm_O, &              !< CRES
         xm_H * 10 + xm_C * 8, &                    !< XYL
         xm_S + xm_O * 4, &                         !< SO4a_f
         xm_S + xm_O * 2, &                         !< SO2
         xm_H + xm_O * 2 + xm_N, &                  !< HNO2
         xm_H * 4 + xm_C, &                         !< CH4
         xm_H * 3 + xm_N, &                         !< NH3
         xm_O * 3 + xm_N, &                         !< NO3
         xm_H + xm_O, &                             !< OH
         xm_H + xm_O * 2, &                         !< HO2
         xm_O * 5 + xm_N * 2, &                     !< N2O5
         xm_S + xm_O * 4, &                         !< SO4a_c
         xm_H * 4 + xm_N, &                         !< NH4a_f
         xm_O * 3 + xm_N, &                         !< NO3a_f
         xm_O * 3 + xm_N, &                         !< NO3a_c
         xm_C * 2 + xm_O * 3, &                     !< C2O3
         xm_dummy, &                                !< XO2
         xm_dummy, &                                !< XO2N
         xm_dummy, &                                !< CRO
         xm_H + xm_O * 3 + xm_N, &                  !< HNO3
         xm_H * 2 + xm_O * 2, &                     !< H2O2
         xm_H * 8 + xm_C * 5, &                     !< ISO
         xm_dummy, &                                !< ISPD
         xm_dummy, &                                !< TO2
         xm_dummy, &                                !< OPEN
         xm_H * 16 + xm_C * 10, &                   !< TERP
         xm_dummy, &                                !< EC_f
         xm_dummy, &                                !< EC_c
         xm_dummy, &                                !< POM_f
         xm_dummy, &                                !< POM_c
         xm_dummy, &                                !< PPM_f
         xm_dummy, &                                !< PPM_c
         xm_Na, &                                   !< Na_ff
         xm_Na, &                                   !< Na_f
         xm_Na, &                                   !< Na_c
         xm_Na, &                                   !< Na_cc
         xm_Na, &                                   !< Na_ccc
         xm_dummy, &                                !< dust_ff
         xm_dummy, &                                !< dust_f
         xm_dummy, &                                !< dust_c
         xm_dummy, &                                !< dust_cc
         xm_dummy, &                                !< dust_ccc
         xm_dummy, &                                !< tpm10
         xm_dummy, &                                !< tpm25
         xm_dummy, &                                !< tss
         xm_dummy, &                                !< tdust
         xm_dummy, &                                !< tc
         xm_dummy, &                                !< tcg
         xm_dummy, &                                !< tsoa
         xm_dummy, &                                !< tnmvoc
         xm_dummy, &                                !< SOxa
         xm_dummy, &                                !< NOya
         xm_dummy, &                                !< NHxa
         xm_O * 2 + xm_N, &                         !< NO2_obs
         xm_dummy, &                                !< tpm10_biascorr
         xm_dummy, &                                !< tpm25_biascorr
         xm_O * 3 /)                                !< o3_biascorr
!
!-- Get current day of the year
    CALL get_date_time( time_since_reference_point, day_of_year = day_of_year )
!
!-- Initialize surface element m
    m = 0
    k = 0
!
!-- LSM or USM surface present at i,j:
!-- Default surfaces are NOT considered for deposition
    match_lsm = surf_lsm_h%start_index(j,i) <= surf_lsm_h%end_index(j,i)
    match_usm = surf_usm_h%start_index(j,i) <= surf_usm_h%end_index(j,i)
!
!--For LSM surfaces

    IF ( match_lsm )  THEN
!
!--    Get surface element information at i,j:
       m = surf_lsm_h%start_index(j,i)
       k = surf_lsm_h%k(m)
!
!--    Get needed variables for surface element m
       ustar_surf  = surf_lsm_h%us(m)
       z0h_surf    = surf_lsm_h%z0h(m)
       r_aero_surf = surf_lsm_h%r_a(m)
       solar_rad   = surf_lsm_h%rad_sw_dir(m) + surf_lsm_h%rad_sw_dif(m)
       
       lai = surf_lsm_h%lai(m)
       sai = lai + 1
!
!--    For small grid spacing neglect R_a
       IF ( dzw(k) <= 1.0 )  THEN
          r_aero_surf = 0.0_wp
       ENDIF
!
!--    Initialize lu's
       luv_palm = 0
       luv_dep = 0
       lup_palm = 0
       lup_dep = 0
       luw_palm = 0
       luw_dep = 0
! 
!--    Initialize budgets
       bud_luv = 0.0_wp
       bud_lup = 0.0_wp
       bud_luw = 0.0_wp
!
!--    Get land use for i,j and assign to DEPAC lu
       IF ( surf_lsm_h%frac(m,ind_veg_wall) > 0 )  THEN
          luv_palm = surf_lsm_h%vegetation_type(m)
          IF ( luv_palm == ind_luv_user )  THEN
             message_string = 'No vegetation type defined. Please define vegetation type to enable deposition calculation'
             CALL message( 'chem_depo', 'CM0451', 1, 2, 0, 6, 0 )
          ELSEIF ( luv_palm == ind_luv_b_soil )  THEN
             luv_dep = 9
          ELSEIF ( luv_palm == ind_luv_mixed_crops )  THEN
             luv_dep = 2
          ELSEIF ( luv_palm == ind_luv_s_grass )  THEN
             luv_dep = 1
          ELSEIF ( luv_palm == ind_luv_ev_needle_trees )  THEN
             luv_dep = 4
          ELSEIF ( luv_palm == ind_luv_de_needle_trees )  THEN
             luv_dep = 4
          ELSEIF ( luv_palm == ind_luv_ev_broad_trees )  THEN
             luv_dep = 12
          ELSEIF ( luv_palm == ind_luv_de_broad_trees )  THEN
             luv_dep = 5
          ELSEIF ( luv_palm == ind_luv_t_grass )  THEN
             luv_dep = 1
          ELSEIF ( luv_palm == ind_luv_desert )  THEN
             luv_dep = 9
          ELSEIF ( luv_palm == ind_luv_tundra )  THEN
             luv_dep = 8
          ELSEIF ( luv_palm == ind_luv_irr_crops )  THEN
             luv_dep = 2
          ELSEIF ( luv_palm == ind_luv_semidesert )  THEN
             luv_dep = 8
          ELSEIF ( luv_palm == ind_luv_ice )  THEN
             luv_dep = 10
          ELSEIF ( luv_palm == ind_luv_marsh )  THEN
             luv_dep = 8
          ELSEIF ( luv_palm == ind_luv_ev_shrubs )  THEN
             luv_dep = 14
          ELSEIF ( luv_palm == ind_luv_de_shrubs )  THEN
             luv_dep = 14
          ELSEIF ( luv_palm == ind_luv_mixed_forest )  THEN
             luv_dep = 4
          ELSEIF ( luv_palm == ind_luv_intrup_forest )  THEN
             luv_dep = 8      
          ENDIF
       ENDIF

       IF ( surf_lsm_h%frac(m,ind_pav_green) > 0 )  THEN
          lup_palm = surf_lsm_h%pavement_type(m)
          IF ( lup_palm == ind_lup_user )  THEN
             message_string = 'No pavement type defined. Please define pavement type to enable deposition calculation'
             CALL message( 'chem_depo', 'CM0452', 1, 2, 0, 6, 0 )
          ELSEIF ( lup_palm == ind_lup_asph_conc )  THEN
             lup_dep = 9
          ELSEIF ( lup_palm == ind_lup_asph )  THEN
             lup_dep = 9
          ELSEIF ( lup_palm == ind_lup_conc )  THEN
             lup_dep = 9
          ELSEIF ( lup_palm == ind_lup_sett )  THEN
             lup_dep = 9
          ELSEIF ( lup_palm == ind_lup_pav_stones )  THEN
             lup_dep = 9
          ELSEIF ( lup_palm == ind_lup_cobblest )  THEN
             lup_dep = 9       
          ELSEIF ( lup_palm == ind_lup_metal )  THEN
             lup_dep = 9
          ELSEIF ( lup_palm == ind_lup_wood )  THEN
             lup_dep = 9    
          ELSEIF ( lup_palm == ind_lup_gravel )  THEN
             lup_dep = 9
          ELSEIF ( lup_palm == ind_lup_f_gravel )  THEN
             lup_dep = 9
          ELSEIF ( lup_palm == ind_lup_pebblest )  THEN
             lup_dep = 9
          ELSEIF ( lup_palm == ind_lup_woodchips )  THEN
             lup_dep = 9
          ELSEIF ( lup_palm == ind_lup_tartan )  THEN
             lup_dep = 9
          ELSEIF ( lup_palm == ind_lup_art_turf )  THEN
             lup_dep = 9
          ELSEIF ( lup_palm == ind_lup_clay )  THEN
             lup_dep = 9
          ENDIF
       ENDIF

       IF ( surf_lsm_h%frac(m,ind_wat_win) > 0 )  THEN
          luw_palm = surf_lsm_h%water_type(m)     
          IF ( luw_palm == ind_luw_user )  THEN
             message_string = 'No water type defined. Please define water type to enable deposition calculation'
             CALL message( 'chem_depo', 'CM0453', 1, 2, 0, 6, 0 )
          ELSEIF ( luw_palm ==  ind_luw_lake )  THEN
             luw_dep = 13
          ELSEIF ( luw_palm == ind_luw_river )  THEN
             luw_dep = 13
          ELSEIF ( luw_palm == ind_luw_ocean )  THEN
             luw_dep = 6
          ELSEIF ( luw_palm == ind_luw_pond )  THEN
             luw_dep = 13
          ELSEIF ( luw_palm == ind_luw_fountain )  THEN
             luw_dep = 13 
          ENDIF
       ENDIF
!
!--    Set wetness indicator to dry or wet for lsm vegetation or pavement
       IF ( surf_lsm_h%c_liq(m) > 0 )  THEN
          nwet = 1
       ELSE
          nwet = 0
       ENDIF
!
!--    Compute length of time step 
       IF ( call_chem_at_all_substeps )  THEN
          dt_chem = dt_3d * weight_pres(intermediate_timestep_count)
       ELSE
          dt_chem = dt_3d
       ENDIF

       dh = dzw(k)
       inv_dh = 1.0_wp / dh
       dt_dh = dt_chem/dh
!
!--    Concentration at i,j,k
       DO  lsp = 1, nspec
          conc_ijk(lsp) = chem_species(lsp)%conc(k,j,i)
       ENDDO

!--    Temperature at i,j,k
       temp_tmp = pt(k,j,i) * ( hyp(k) / 100000.0_wp )**0.286_wp

       ts       = temp_tmp - 273.15  !< in degrees celcius
!
!--    Viscosity of air
       visc = 1.496e-6 * temp_tmp**1.5 / (temp_tmp + 120.0)
!
!--    Air density at k
       dens = rho_air_zw(k)
!
!--    Calculate relative humidity from specific humidity for DEPAC
       qv_tmp = MAX(q(k,j,i),0.0_wp)
       rh_surf = relativehumidity_from_specifichumidity(qv_tmp, temp_tmp, hyp(k) )
!
!-- Check if surface fraction (vegetation, pavement or water) > 0 and calculate vd and budget
!-- for each surface fraction. Then derive overall budget taking into account the surface fractions.
!
!--    Vegetation
       IF ( surf_lsm_h%frac(m,ind_veg_wall) > 0 )  THEN

!
!--       No vegetation on bare soil, desert or ice:
          IF ( ( luv_palm == ind_luv_b_soil ) .OR. &
                ( luv_palm == ind_luv_desert ) .OR. &
                 ( luv_palm == ind_luv_ice ) ) THEN

             lai = 0.0_wp
             sai = 0.0_wp

          ENDIF
          
          slinnfac = 1.0_wp
!
!--       Get deposition velocity vd
          DO  lsp = 1, nvar
!
!--          Initialize 
             vs = 0.0_wp
             vd_lu = 0.0_wp
             rs = 0.0_wp
             rb = 0.0_wp
             rc_tot = 0.0_wp
             IF ( spc_names(lsp) == 'PM10' )  THEN
                part_type = 1
!
!--             Sedimentation velocity
                vs = slinnfac * sedimentation_velocity( particle_pars(ind_p_dens, part_type), &
                     particle_pars(ind_p_size, part_type), &
                     particle_pars(ind_p_slip, part_type), &
                     visc)

                CALL drydepo_aero_zhang_vd( vd_lu, rs, &
                     vs, &
                     particle_pars(ind_p_size, part_type), &
                     particle_pars(ind_p_slip, part_type), &
                     nwet, temp_tmp, dens, visc, &
                     luv_dep,  &
                     r_aero_surf, ustar_surf )

                bud_luv(lsp) = - conc_ijk(lsp) * &
                     (1.0_wp - exp(-vd_lu * dt_dh )) * dh


             ELSEIF ( spc_names(lsp) == 'PM25' )  THEN
                part_type = 2
!
!--             Sedimentation velocity
                vs = slinnfac * sedimentation_velocity( particle_pars(ind_p_dens, part_type), &
                     particle_pars(ind_p_size, part_type), &
                     particle_pars(ind_p_slip, part_type), &
                     visc)

                CALL drydepo_aero_zhang_vd( vd_lu, rs, &
                     vs, &
                     particle_pars(ind_p_size, part_type), &
                     particle_pars(ind_p_slip, part_type), &
                     nwet, temp_tmp, dens, visc, &
                     luv_dep , &
                     r_aero_surf, ustar_surf )

                bud_luv(lsp) = - conc_ijk(lsp) * &
                     (1.0_wp - exp(-vd_lu * dt_dh )) * dh

             ELSE !< GASES
!
!--             Read spc_name of current species for gas parameter
                IF ( ANY( pspecnames(:) == spc_names(lsp) ) )  THEN
                   i_pspec = 0
                   DO  pspec = 1, nposp
                      IF ( pspecnames(pspec) == spc_names(lsp) )  THEN
                         i_pspec = pspec
                      END IF
                   ENDDO

                ELSE
!
!--             For now species not deposited
                   CYCLE
                ENDIF
!
!--             Factor used for conversion from ppb to ug/m3 :
!--             ppb (mole tr)/(mole air)/ppb (kg tr)/(mole tr) (ug tr)/(kg tr) &
!--                 (mole air)/(kg air) (kg air)/(m3 air) (kg air(ug/m3)/ppb/(kg/mole) = / (kg/mole)
!--                 c           1e-9              xm_tracer         1e9       /       xm_air            dens
!--             thus:
!--                 c_in_ppb * xm_tracer * [ dens / xm_air ] = c_in_ugm3
!--             Use density at k:

                ppm2ugm3 =  (dens/xm_air) * 0.001_wp  !< (mole air)/m3
!
!--             Atmospheric concentration in DEPAC is requested in ug/m3:
                !   ug/m3              ppm          (ug/m3)/ppm/(kg/mole)     kg/mole
                conc_ijk_ugm3 = conc_ijk(lsp)         * ppm2ugm3 *   specmolm(i_pspec)  ! in ug/m3
!
!--             Diffusivity for DEPAC relevant gases
!--             Use default value 
                diffusivity            = 0.11e-4
!
!--             overwrite with known coefficients of diffusivity from Massman (1998)
                IF ( spc_names(lsp) == 'NO2' ) diffusivity = 0.136e-4 
                IF ( spc_names(lsp) == 'NO'  ) diffusivity = 0.199e-4
                IF ( spc_names(lsp) == 'O3'  ) diffusivity = 0.144e-4
                IF ( spc_names(lsp) == 'CO'  ) diffusivity = 0.176e-4
                IF ( spc_names(lsp) == 'SO2' ) diffusivity = 0.112e-4
                IF ( spc_names(lsp) == 'CH4' ) diffusivity = 0.191e-4
                IF ( spc_names(lsp) == 'NH3' ) diffusivity = 0.191e-4
!
!--             Get quasi-laminar boundary layer resistance rb:
                CALL get_rb_cell( (luv_dep == ilu_water_sea) .OR. (luv_dep == ilu_water_inland), &
                     z0h_surf, ustar_surf, diffusivity, &
                     rb )
!
!--             Get rc_tot
                CALL drydepos_gas_depac( spc_names(lsp), day_of_year, latitude, ts, ustar_surf, solar_rad, cos_zenith, &
                     rh_surf, lai, sai, nwet, luv_dep, 2, rc_tot, ccomp_tot(lsp), hyp(nzb), conc_ijk_ugm3, diffusivity, &
                     r_aero_surf , rb )
!
!--             Calculate budget
                IF ( rc_tot <= 0.0 )  THEN

                   bud_luv(lsp) = 0.0_wp

                ELSE

                   vd_lu = 1.0_wp / (r_aero_surf + rb + rc_tot )

                   bud_luv(lsp) = - (conc_ijk(lsp) - ccomp_tot(lsp)) * &
                        (1.0_wp - exp(-vd_lu * dt_dh )) * dh
                ENDIF

             ENDIF
          ENDDO
       ENDIF
!
!--    Pavement
       IF ( surf_lsm_h%frac(m,ind_pav_green) > 0 )  THEN
!
!--       No vegetation on pavements:
          lai = 0.0_wp
          sai = 0.0_wp

          slinnfac = 1.0_wp
!
!--       Get vd
          DO  lsp = 1, nvar
!
!--       Initialize
             vs = 0.0_wp
             vd_lu = 0.0_wp
             rs = 0.0_wp
             rb = 0.0_wp
             rc_tot = 0.0_wp
             IF ( spc_names(lsp) == 'PM10' )  THEN
                part_type = 1
!
!--             Sedimentation velocity:
                vs = slinnfac * sedimentation_velocity( particle_pars(ind_p_dens, part_type), &
                     particle_pars(ind_p_size, part_type), &
                     particle_pars(ind_p_slip, part_type), &
                     visc)

                CALL drydepo_aero_zhang_vd( vd_lu, rs, &
                     vs, &
                     particle_pars(ind_p_size, part_type), &
                     particle_pars(ind_p_slip, part_type), &
                     nwet, temp_tmp, dens, visc, &
                     lup_dep,  &
                     r_aero_surf, ustar_surf )

                bud_lup(lsp) = - conc_ijk(lsp) * &
                     (1.0_wp - exp(-vd_lu * dt_dh )) * dh


             ELSEIF ( spc_names(lsp) == 'PM25' )  THEN
                part_type = 2
!
!--             Sedimentation velocity:
                vs = slinnfac * sedimentation_velocity( particle_pars(ind_p_dens, part_type), &
                     particle_pars(ind_p_size, part_type), &
                     particle_pars(ind_p_slip, part_type), &
                     visc)

                CALL drydepo_aero_zhang_vd( vd_lu, rs, &
                     vs, &
                     particle_pars(ind_p_size, part_type), &
                     particle_pars(ind_p_slip, part_type), &
                     nwet, temp_tmp, dens, visc, &
                     lup_dep, &
                     r_aero_surf, ustar_surf )

                bud_lup(lsp) = - conc_ijk(lsp) * &
                     (1.0_wp - exp(-vd_lu * dt_dh )) * dh

             ELSE  !<GASES
!
!--             Read spc_name of current species for gas parameter

                IF ( ANY(pspecnames(:) == spc_names(lsp) ) )  THEN
                   i_pspec = 0
                   DO  pspec = 1, nposp
                      IF ( pspecnames(pspec) == spc_names(lsp) )  THEN
                         i_pspec = pspec
                      END IF
                   ENDDO

                ELSE
!
!--                For now species not deposited
                   CYCLE
                ENDIF
!
!--             Factor used for conversion from ppb to ug/m3 :
!--                 ppb (mole tr)/(mole air)/ppb (kg tr)/(mole tr) (ug tr)/(kg tr) &
!--                 (mole air)/(kg air) (kg air)/(m3 air) (kg air(ug/m3)/ppb/(kg/mole) = / (kg/mole)
!--                 c           1e-9               xm_tracer         1e9       /       xm_air            dens
!--             thus:
!--                 c_in_ppb * xm_tracer * [ dens / xm_air ] = c_in_ugm3
!--             Use density at lowest layer:

                ppm2ugm3 =  (dens/xm_air) * 0.001_wp  !< (mole air)/m3
!
!--             Atmospheric concentration in DEPAC is requested in ug/m3:
                !   ug/m3              ppm          (ug/m3)/ppm/(kg/mole)     kg/mole
                conc_ijk_ugm3 = conc_ijk(lsp)         * ppm2ugm3 *   specmolm(i_pspec)  ! in ug/m3
!
!--             Diffusivity for DEPAC relevant gases
!--             Use default value 
                diffusivity            = 0.11e-4
!
!--             overwrite with known coefficients of diffusivity from Massman (1998)
                IF ( spc_names(lsp) == 'NO2' ) diffusivity = 0.136e-4 
                IF ( spc_names(lsp) == 'NO'  ) diffusivity = 0.199e-4
                IF ( spc_names(lsp) == 'O3'  ) diffusivity = 0.144e-4
                IF ( spc_names(lsp) == 'CO'  ) diffusivity = 0.176e-4
                IF ( spc_names(lsp) == 'SO2' ) diffusivity = 0.112e-4
                IF ( spc_names(lsp) == 'CH4' ) diffusivity = 0.191e-4
                IF ( spc_names(lsp) == 'NH3' ) diffusivity = 0.191e-4
!
!--             Get quasi-laminar boundary layer resistance rb:
                CALL get_rb_cell( (lup_dep == ilu_water_sea) .OR. (lup_dep == ilu_water_inland),   &
                     z0h_surf, ustar_surf, diffusivity, rb )
!
!--             Get rc_tot
                CALL drydepos_gas_depac( spc_names(lsp), day_of_year, latitude, ts, ustar_surf,      &
                                         solar_rad, cos_zenith, rh_surf, lai, sai, nwet, lup_dep, 2,    &
                                         rc_tot, ccomp_tot(lsp), hyp(nzb), conc_ijk_ugm3, diffusivity,              &
                                         r_aero_surf , rb )
!
!--             Calculate budget
                IF ( rc_tot <= 0.0 )  THEN
                   bud_lup(lsp) = 0.0_wp
                ELSE
                   vd_lu = 1.0_wp / (r_aero_surf + rb + rc_tot )
                   bud_lup(lsp) = - (conc_ijk(lsp) - ccomp_tot(lsp)) * &
                        (1.0_wp - exp(-vd_lu * dt_dh )) * dh
                ENDIF

             ENDIF
          ENDDO
       ENDIF
!
!--    Water
       IF ( surf_lsm_h%frac(m,ind_wat_win) > 0 )  THEN
!
!--       No vegetation on water:
          lai = 0.0_wp
          sai = 0.0_wp
          slinnfac = 1.0_wp
!
!--       Get vd
          DO  lsp = 1, nvar
!
!--          Initialize
             vs = 0.0_wp
             vd_lu = 0.0_wp
             rs = 0.0_wp
             rb = 0.0_wp
             rc_tot = 0.0_wp 
             IF ( spc_names(lsp) == 'PM10' )  THEN
                part_type = 1
!
!--             Sedimentation velocity:
                vs = slinnfac * sedimentation_velocity( particle_pars(ind_p_dens, part_type), &
                     particle_pars(ind_p_size, part_type), &
                     particle_pars(ind_p_slip, part_type), &
                     visc)

                CALL drydepo_aero_zhang_vd( vd_lu, rs, &
                     vs, &
                     particle_pars(ind_p_size, part_type), &
                     particle_pars(ind_p_slip, part_type), &
                     nwet, temp_tmp, dens, visc, &
                     luw_dep, &
                     r_aero_surf, ustar_surf )

                bud_luw(lsp) = - conc_ijk(lsp) * &
                     (1.0_wp - exp(-vd_lu * dt_dh )) * dh

             ELSEIF ( spc_names(lsp) == 'PM25' )  THEN
                part_type = 2
!
!--             Sedimentation velocity:
                vs = slinnfac * sedimentation_velocity( particle_pars(ind_p_dens, part_type), &
                     particle_pars(ind_p_size, part_type), &
                     particle_pars(ind_p_slip, part_type), &
                     visc)

                CALL drydepo_aero_zhang_vd( vd_lu, rs, &
                     vs, &
                     particle_pars(ind_p_size, part_type), &
                     particle_pars(ind_p_slip, part_type), &
                     nwet, temp_tmp, dens, visc, &
                     luw_dep, &
                     r_aero_surf, ustar_surf )

                bud_luw(lsp) = - conc_ijk(lsp) * &
                     (1.0_wp - exp(-vd_lu * dt_dh )) * dh

             ELSE  !<GASES
!
!--             Read spc_name of current species for gas PARAMETER

                IF ( ANY(pspecnames(:) == spc_names(lsp) ) )  THEN
                   i_pspec = 0
                   DO  pspec = 1, nposp
                      IF ( pspecnames(pspec) == spc_names(lsp) )  THEN
                         i_pspec = pspec
                      END IF
                   ENDDO

                ELSE
!
!--                For now species not deposited
                   CYCLE
                ENDIF
!
!--             Factor used for conversion from ppb to ug/m3 :
!--                 ppb (mole tr)/(mole air)/ppb (kg tr)/(mole tr) (ug tr)/(kg tr) &
!--                 (mole air)/(kg air) (kg air)/(m3 air) (kg air(ug/m3)/ppb/(kg/mole) = / (kg/mole)
!--                 c           1e-9               xm_tracer         1e9       /       xm_air            dens
!--             thus:
!--                 c_in_ppb * xm_tracer * [ dens / xm_air ] = c_in_ugm3
!--             Use density at lowest layer:

                ppm2ugm3 = (dens/xm_air) * 0.001_wp  !< (mole air)/m3
!
!--             Atmospheric concentration in DEPAC is requested in ug/m3:
!--                 ug/m3        ppm          (ug/m3)/ppm/(kg/mole)     kg/mole
                conc_ijk_ugm3 = conc_ijk(lsp) * ppm2ugm3 *  specmolm(i_pspec)  ! in ug/m3
!
!--             Diffusivity for DEPAC relevant gases
!--             Use default value 
                diffusivity            = 0.11e-4
!
!--             overwrite with known coefficients of diffusivity from Massman (1998)
                IF ( spc_names(lsp) == 'NO2' ) diffusivity = 0.136e-4 
                IF ( spc_names(lsp) == 'NO'  ) diffusivity = 0.199e-4
                IF ( spc_names(lsp) == 'O3'  ) diffusivity = 0.144e-4
                IF ( spc_names(lsp) == 'CO'  ) diffusivity = 0.176e-4
                IF ( spc_names(lsp) == 'SO2' ) diffusivity = 0.112e-4
                IF ( spc_names(lsp) == 'CH4' ) diffusivity = 0.191e-4
                IF ( spc_names(lsp) == 'NH3' ) diffusivity = 0.191e-4
!
!--             Get quasi-laminar boundary layer resistance rb:
                CALL get_rb_cell( (luw_dep == ilu_water_sea) .OR. (luw_dep == ilu_water_inland),  &
                     z0h_surf, ustar_surf, diffusivity, rb )

!--             Get rc_tot
                CALL drydepos_gas_depac( spc_names(lsp), day_of_year, latitude, ts, ustar_surf,         & 
                                         solar_rad, cos_zenith, rh_surf, lai, sai, nwet, luw_dep, 2,    &
                                         rc_tot, ccomp_tot(lsp), hyp(nzb), conc_ijk_ugm3, diffusivity,  &
                                         r_aero_surf , rb )
!
!--             Calculate budget 
                IF ( rc_tot <= 0.0 )  THEN

                   bud_luw(lsp) = 0.0_wp

                ELSE

                   vd_lu = 1.0_wp / (r_aero_surf + rb + rc_tot )

                   bud_luw(lsp) = - (conc_ijk(lsp) - ccomp_tot(lsp)) * &
                        (1.0_wp - exp(-vd_lu * dt_dh )) * dh
                ENDIF

             ENDIF
          ENDDO
       ENDIF


       bud = 0.0_wp
!
!--    Calculate overall budget for surface m and adapt concentration
       DO  lsp = 1, nspec

          bud(lsp) = surf_lsm_h%frac(m,ind_veg_wall) * bud_luv(lsp) + &
               surf_lsm_h%frac(m,ind_pav_green) * bud_lup(lsp) + &
               surf_lsm_h%frac(m,ind_wat_win) * bud_luw(lsp)
!
!--       Compute new concentration:
          conc_ijk(lsp) = conc_ijk(lsp) + bud(lsp) * inv_dh

          chem_species(lsp)%conc(k,j,i) = MAX(0.0_wp, conc_ijk(lsp))

       ENDDO

    ENDIF
!
!-- For USM surfaces    

    IF ( match_usm )  THEN
!
!--    Get surface element information at i,j:
       m = surf_usm_h%start_index(j,i)
       k = surf_usm_h%k(m)
!
!--    Get needed variables for surface element m
       ustar_surf  = surf_usm_h%us(m)
       z0h_surf    = surf_usm_h%z0h(m)
       r_aero_surf = surf_usm_h%r_a(m)
       solar_rad   = surf_usm_h%rad_sw_dir(m) + surf_usm_h%rad_sw_dif(m)
       lai = surf_usm_h%lai(m)
       sai = lai + 1
!
!--    For small grid spacing neglect R_a
       IF ( dzw(k) <= 1.0 )  THEN
          r_aero_surf = 0.0_wp
       ENDIF
!
!--    Initialize lu's
       luu_palm = 0
       luu_dep = 0
       lug_palm = 0
       lug_dep = 0
       lud_palm = 0
       lud_dep = 0
!
!--    Initialize budgets
       bud_luu  = 0.0_wp
       bud_lug = 0.0_wp
       bud_lud = 0.0_wp
!
!--    Get land use for i,j and assign to DEPAC lu
       IF ( surf_usm_h%frac(m,ind_pav_green) > 0 )  THEN
!
!--       For green urban surfaces (e.g. green roofs
!--       assume LU short grass
          lug_palm = ind_luv_s_grass
          IF ( lug_palm == ind_luv_user )  THEN
             message_string = 'No vegetation type defined. Please define vegetation type to enable deposition calculation'
             CALL message( 'chem_depo', 'CM0454', 1, 2, 0, 6, 0 )
          ELSEIF ( lug_palm == ind_luv_b_soil )  THEN
             lug_dep = 9
          ELSEIF ( lug_palm == ind_luv_mixed_crops )  THEN
             lug_dep = 2
          ELSEIF ( lug_palm == ind_luv_s_grass )  THEN
             lug_dep = 1
          ELSEIF ( lug_palm == ind_luv_ev_needle_trees )  THEN
             lug_dep = 4
          ELSEIF ( lug_palm == ind_luv_de_needle_trees )  THEN
             lug_dep = 4
          ELSEIF ( lug_palm == ind_luv_ev_broad_trees )  THEN
             lug_dep = 12
          ELSEIF ( lug_palm == ind_luv_de_broad_trees )  THEN
             lug_dep = 5
          ELSEIF ( lug_palm == ind_luv_t_grass )  THEN
             lug_dep = 1
          ELSEIF ( lug_palm == ind_luv_desert )  THEN
             lug_dep = 9
          ELSEIF ( lug_palm == ind_luv_tundra  )  THEN
             lug_dep = 8
          ELSEIF ( lug_palm == ind_luv_irr_crops )  THEN
             lug_dep = 2
          ELSEIF ( lug_palm == ind_luv_semidesert )  THEN
             lug_dep = 8
          ELSEIF ( lug_palm == ind_luv_ice )  THEN
             lug_dep = 10
          ELSEIF ( lug_palm == ind_luv_marsh )  THEN
             lug_dep = 8
          ELSEIF ( lug_palm == ind_luv_ev_shrubs )  THEN
             lug_dep = 14
          ELSEIF ( lug_palm == ind_luv_de_shrubs  )  THEN
             lug_dep = 14
          ELSEIF ( lug_palm == ind_luv_mixed_forest )  THEN
             lug_dep = 4
          ELSEIF ( lug_palm == ind_luv_intrup_forest )  THEN
             lug_dep = 8      
          ENDIF
       ENDIF

       IF ( surf_usm_h%frac(m,ind_veg_wall) > 0 )  THEN
!
!--       For walls in USM assume concrete walls/roofs,
!--       assumed LU class desert as also assumed for
!--       pavements in LSM
          luu_palm = ind_lup_conc
          IF ( luu_palm == ind_lup_user )  THEN
             message_string = 'No pavement type defined. Please define pavement type to enable deposition calculation'
             CALL message( 'chem_depo', 'CM0455', 1, 2, 0, 6, 0 )
          ELSEIF ( luu_palm == ind_lup_asph_conc )  THEN
             luu_dep = 9
          ELSEIF ( luu_palm == ind_lup_asph )  THEN
             luu_dep = 9
          ELSEIF ( luu_palm ==  ind_lup_conc )  THEN
             luu_dep = 9
          ELSEIF ( luu_palm ==  ind_lup_sett )  THEN
             luu_dep = 9
          ELSEIF ( luu_palm == ind_lup_pav_stones )  THEN
             luu_dep = 9
          ELSEIF ( luu_palm == ind_lup_cobblest )  THEN
             luu_dep = 9       
          ELSEIF ( luu_palm == ind_lup_metal )  THEN
             luu_dep = 9
          ELSEIF ( luu_palm == ind_lup_wood )  THEN
             luu_dep = 9    
          ELSEIF ( luu_palm == ind_lup_gravel )  THEN
             luu_dep = 9
          ELSEIF ( luu_palm == ind_lup_f_gravel )  THEN
             luu_dep = 9
          ELSEIF ( luu_palm == ind_lup_pebblest )  THEN
             luu_dep = 9
          ELSEIF ( luu_palm == ind_lup_woodchips )  THEN
             luu_dep = 9
          ELSEIF ( luu_palm == ind_lup_tartan )  THEN
             luu_dep = 9
          ELSEIF ( luu_palm == ind_lup_art_turf )  THEN
             luu_dep = 9
          ELSEIF ( luu_palm == ind_lup_clay )  THEN
             luu_dep = 9
          ENDIF
       ENDIF

       IF ( surf_usm_h%frac(m,ind_wat_win) > 0 )  THEN
!
!--       For windows in USM assume metal as this is
!--       as close as we get, assumed LU class desert
!--       as also assumed for pavements in LSM
          lud_palm = ind_lup_metal     
          IF ( lud_palm == ind_lup_user )  THEN
             message_string = 'No pavement type defined. Please define pavement type to enable deposition calculation'
             CALL message( 'chem_depo', 'CM0456', 1, 2, 0, 6, 0 )
          ELSEIF ( lud_palm == ind_lup_asph_conc )  THEN
             lud_dep = 9
          ELSEIF ( lud_palm == ind_lup_asph )  THEN
             lud_dep = 9
          ELSEIF ( lud_palm ==  ind_lup_conc )  THEN
             lud_dep = 9
          ELSEIF ( lud_palm ==  ind_lup_sett )  THEN
             lud_dep = 9
          ELSEIF ( lud_palm == ind_lup_pav_stones )  THEN
             lud_dep = 9
          ELSEIF ( lud_palm == ind_lup_cobblest )  THEN
             lud_dep = 9       
          ELSEIF ( lud_palm == ind_lup_metal )  THEN
             lud_dep = 9
          ELSEIF ( lud_palm == ind_lup_wood )  THEN
             lud_dep = 9    
          ELSEIF ( lud_palm == ind_lup_gravel )  THEN
             lud_dep = 9
          ELSEIF ( lud_palm == ind_lup_f_gravel )  THEN
             lud_dep = 9
          ELSEIF ( lud_palm == ind_lup_pebblest )  THEN
             lud_dep = 9
          ELSEIF ( lud_palm == ind_lup_woodchips )  THEN
             lud_dep = 9
          ELSEIF ( lud_palm == ind_lup_tartan )  THEN
             lud_dep = 9
          ELSEIF ( lud_palm == ind_lup_art_turf )  THEN
             lud_dep = 9
          ELSEIF ( lud_palm == ind_lup_clay )  THEN
             lud_dep = 9
          ENDIF
       ENDIF
!
!--    @TODO: Activate these lines as soon as new ebsolver branch is merged:
!--    Set wetness indicator to dry or wet for usm vegetation or pavement
       !IF ( surf_usm_h%c_liq(m) > 0 )  THEN
       !   nwet = 1
       !ELSE
       nwet = 0
       !ENDIF
!
!--    Compute length of time step 
       IF ( call_chem_at_all_substeps )  THEN
          dt_chem = dt_3d * weight_pres(intermediate_timestep_count)
       ELSE
          dt_chem = dt_3d
       ENDIF

       dh = dzw(k)
       inv_dh = 1.0_wp / dh
       dt_dh = dt_chem/dh
!
!--    Concentration at i,j,k
       DO  lsp = 1, nspec
          conc_ijk(lsp) = chem_species(lsp)%conc(k,j,i)
       ENDDO
!
!--    Temperature at i,j,k
       temp_tmp = pt(k,j,i) * ( hyp(k) / 100000.0_wp )**0.286_wp

       ts       = temp_tmp - 273.15  !< in degrees celcius
!
!--    Viscosity of air
       visc = 1.496e-6 * temp_tmp**1.5 / (temp_tmp + 120.0)
!
!--    Air density at k
       dens = rho_air_zw(k)
!
!--    Calculate relative humidity from specific humidity for DEPAC
       qv_tmp = MAX(q(k,j,i),0.0_wp)
       rh_surf = relativehumidity_from_specifichumidity(qv_tmp, temp_tmp, hyp(k) )
!
!--    Check if surface fraction (vegetation, pavement or water) > 0 and calculate vd and budget
!--    for each surface fraction. Then derive overall budget taking into account the surface fractions.
!
!--    Walls/roofs
       IF ( surf_usm_h%frac(m,ind_veg_wall) > 0 )  THEN
!
!--       No vegetation on non-green walls:
          lai = 0.0_wp
          sai = 0.0_wp

          slinnfac = 1.0_wp
!
!--       Get vd
          DO  lsp = 1, nvar
!
!--          Initialize 
             vs = 0.0_wp
             vd_lu = 0.0_wp
             rs = 0.0_wp
             rb = 0.0_wp
             rc_tot = 0.0_wp
             IF (spc_names(lsp) == 'PM10' )  THEN
                part_type = 1
!
!--             Sedimentation velocity
                vs = slinnfac * sedimentation_velocity( particle_pars(ind_p_dens, part_type), &
                     particle_pars(ind_p_size, part_type), &
                     particle_pars(ind_p_slip, part_type), &
                     visc)

                CALL drydepo_aero_zhang_vd( vd_lu, rs, &
                     vs, &
                     particle_pars(ind_p_size, part_type), &
                     particle_pars(ind_p_slip, part_type), &
                     nwet, temp_tmp, dens, visc, &
                     luu_dep,  &
                     r_aero_surf, ustar_surf )

                bud_luu(lsp) = - conc_ijk(lsp) * &
                     (1.0_wp - exp(-vd_lu * dt_dh )) * dh

             ELSEIF ( spc_names(lsp) == 'PM25' )  THEN
                part_type = 2
!
!--             Sedimentation velocity
                vs = slinnfac * sedimentation_velocity( particle_pars(ind_p_dens, part_type), &
                     particle_pars(ind_p_size, part_type), &
                     particle_pars(ind_p_slip, part_type), &
                     visc)

                CALL drydepo_aero_zhang_vd( vd_lu, rs, &
                     vs, &
                     particle_pars(ind_p_size, part_type), &
                     particle_pars(ind_p_slip, part_type), &
                     nwet, temp_tmp, dens, visc, &
                     luu_dep , &
                     r_aero_surf, ustar_surf )

                bud_luu(lsp) = - conc_ijk(lsp) * &
                     (1.0_wp - exp(-vd_lu * dt_dh )) * dh

             ELSE  !< GASES
!
!--             Read spc_name of current species for gas parameter

                IF ( ANY( pspecnames(:) == spc_names(lsp) ) )  THEN
                   i_pspec = 0
                   DO  pspec = 1, nposp
                      IF ( pspecnames(pspec) == spc_names(lsp) )  THEN
                         i_pspec = pspec
                      END IF
                   ENDDO
                ELSE
!
!--                For now species not deposited
                   CYCLE
                ENDIF
!
!--             Factor used for conversion from ppb to ug/m3 :
!--                 ppb (mole tr)/(mole air)/ppb (kg tr)/(mole tr) (ug tr)/(kg tr) &
!--                 (mole air)/(kg air) (kg air)/(m3 air) (kg air(ug/m3)/ppb/(kg/mole) = / (kg/mole)
!--                 c           1e-9              xm_tracer         1e9       /       xm_air            dens
!--             thus:
!--                 c_in_ppb * xm_tracer * [ dens / xm_air ] = c_in_ugm3
!--             Use density at k:

                ppm2ugm3 =  (dens/xm_air) * 0.001_wp  !< (mole air)/m3

                !
!--             Atmospheric concentration in DEPAC is requested in ug/m3:
!--                 ug/m3              ppm          (ug/m3)/ppm/(kg/mole)     kg/mole
                conc_ijk_ugm3 = conc_ijk(lsp)         * ppm2ugm3 *   specmolm(i_pspec)  ! in ug/m3
!
!--             Diffusivity for DEPAC relevant gases
!--             Use default value 
                diffusivity            = 0.11e-4
!
!--             overwrite with known coefficients of diffusivity from Massman (1998)
                IF ( spc_names(lsp) == 'NO2' ) diffusivity = 0.136e-4 
                IF ( spc_names(lsp) == 'NO'  ) diffusivity = 0.199e-4
                IF ( spc_names(lsp) == 'O3'  ) diffusivity = 0.144e-4
                IF ( spc_names(lsp) == 'CO'  ) diffusivity = 0.176e-4
                IF ( spc_names(lsp) == 'SO2' ) diffusivity = 0.112e-4
                IF ( spc_names(lsp) == 'CH4' ) diffusivity = 0.191e-4
                IF ( spc_names(lsp) == 'NH3' ) diffusivity = 0.191e-4
!
!--             Get quasi-laminar boundary layer resistance rb:
                CALL get_rb_cell( (luu_dep == ilu_water_sea) .OR. (luu_dep == ilu_water_inland),   &
                     z0h_surf, ustar_surf, diffusivity, &
                     rb )
!
!--             Get rc_tot
                CALL drydepos_gas_depac( spc_names(lsp), day_of_year, latitude, ts, ustar_surf,          &
                                         solar_rad, cos_zenith, rh_surf, lai, sai, nwet, luu_dep, 2,     &
                                         rc_tot, ccomp_tot(lsp), hyp(nzb), conc_ijk_ugm3, diffusivity,   &
                                         r_aero_surf, rb )
!
!--             Calculate budget
                IF ( rc_tot <= 0.0 )  THEN

                   bud_luu(lsp) = 0.0_wp

                ELSE

                   vd_lu = 1.0_wp / (r_aero_surf + rb + rc_tot )

                   bud_luu(lsp) = - (conc_ijk(lsp) - ccomp_tot(lsp)) * &
                        (1.0_wp - exp(-vd_lu * dt_dh )) * dh
                ENDIF

             ENDIF
          ENDDO
       ENDIF
!
!--    Green usm surfaces
       IF ( surf_usm_h%frac(m,ind_pav_green) > 0 )  THEN

!
!--       No vegetation on bare soil, desert or ice:
          IF ( ( lug_palm == ind_luv_b_soil ) .OR. &
                ( lug_palm == ind_luv_desert ) .OR. &
                 ( lug_palm == ind_luv_ice ) ) THEN

             lai = 0.0_wp
             sai = 0.0_wp

          ENDIF

          
          slinnfac = 1.0_wp
!
!--       Get vd
          DO  lsp = 1, nvar
!
!--          Initialize
             vs = 0.0_wp
             vd_lu = 0.0_wp
             rs = 0.0_wp
             rb = 0.0_wp
             rc_tot = 0.0_wp
             IF ( spc_names(lsp) == 'PM10' )  THEN
                part_type = 1
!
!--             Sedimentation velocity
                vs = slinnfac * sedimentation_velocity( particle_pars(ind_p_dens, part_type), &
                     particle_pars(ind_p_size, part_type), &
                     particle_pars(ind_p_slip, part_type), &
                     visc)

                CALL drydepo_aero_zhang_vd( vd_lu, rs, &
                     vs, &
                     particle_pars(ind_p_size, part_type), &
                     particle_pars(ind_p_slip, part_type), &
                     nwet, temp_tmp, dens, visc, &
                     lug_dep,  &
                     r_aero_surf, ustar_surf )

                bud_lug(lsp) = - conc_ijk(lsp) * &
                     (1.0_wp - exp(-vd_lu * dt_dh )) * dh

             ELSEIF ( spc_names(lsp) == 'PM25' )  THEN
                part_type = 2
!
!--             Sedimentation velocity
                vs = slinnfac * sedimentation_velocity( particle_pars(ind_p_dens, part_type), &
                     particle_pars(ind_p_size, part_type), &
                     particle_pars(ind_p_slip, part_type), &
                     visc)

                CALL drydepo_aero_zhang_vd( vd_lu, rs, &
                     vs, &
                     particle_pars(ind_p_size, part_type), &
                     particle_pars(ind_p_slip, part_type), &
                     nwet, temp_tmp, dens, visc, &
                     lug_dep, &
                     r_aero_surf, ustar_surf )

                bud_lug(lsp) = - conc_ijk(lsp) * &
                     (1.0_wp - exp(-vd_lu * dt_dh )) * dh

             ELSE  !< GASES
!
!--             Read spc_name of current species for gas parameter

                IF ( ANY( pspecnames(:) == spc_names(lsp) ) )  THEN
                   i_pspec = 0
                   DO  pspec = 1, nposp
                      IF ( pspecnames(pspec) == spc_names(lsp) )  THEN
                         i_pspec = pspec
                      END IF
                   ENDDO
                ELSE
!
!--                For now species not deposited
                   CYCLE
                ENDIF
!
!--             Factor used for conversion from ppb to ug/m3 :
!--                 ppb (mole tr)/(mole air)/ppb (kg tr)/(mole tr) (ug tr)/(kg tr) &
!--                 (mole air)/(kg air) (kg air)/(m3 air) (kg air(ug/m3)/ppb/(kg/mole) = / (kg/mole)
!--                 c           1e-9               xm_tracer         1e9       /       xm_air            dens
!--             thus:
!--                  c_in_ppb * xm_tracer * [ dens / xm_air ] = c_in_ugm3
!--             Use density at k:

                ppm2ugm3 =  (dens/xm_air) * 0.001_wp  ! (mole air)/m3
!
!--             Atmospheric concentration in DEPAC is requested in ug/m3:
                !   ug/m3              ppm          (ug/m3)/ppm/(kg/mole)     kg/mole
                conc_ijk_ugm3 = conc_ijk(lsp)         * ppm2ugm3 *   specmolm(i_pspec)  ! in ug/m3
!
!--             Diffusivity for DEPAC relevant gases
!--             Use default value 
                diffusivity            = 0.11e-4
!
!--             overwrite with known coefficients of diffusivity from Massman (1998)
                IF ( spc_names(lsp) == 'NO2' ) diffusivity = 0.136e-4 
                IF ( spc_names(lsp) == 'NO'  ) diffusivity = 0.199e-4
                IF ( spc_names(lsp) == 'O3'  ) diffusivity = 0.144e-4
                IF ( spc_names(lsp) == 'CO'  ) diffusivity = 0.176e-4
                IF ( spc_names(lsp) == 'SO2' ) diffusivity = 0.112e-4
                IF ( spc_names(lsp) == 'CH4' ) diffusivity = 0.191e-4
                IF ( spc_names(lsp) == 'NH3' ) diffusivity = 0.191e-4
!
!--             Get quasi-laminar boundary layer resistance rb:
                CALL get_rb_cell( (lug_dep == ilu_water_sea) .OR. (lug_dep == ilu_water_inland),    &
                     z0h_surf, ustar_surf, diffusivity, &
                     rb )
!
!--             Get rc_tot
                CALL drydepos_gas_depac( spc_names(lsp), day_of_year, latitude, ts, ustar_surf,           &
                                         solar_rad, cos_zenith, rh_surf, lai, sai, nwet, lug_dep, 2,      &
                                         rc_tot, ccomp_tot(lsp), hyp(nzb), conc_ijk_ugm3, diffusivity,    &
                                         r_aero_surf , rb )
!
!--             Calculate budget
                IF ( rc_tot <= 0.0 )  THEN

                   bud_lug(lsp) = 0.0_wp

                ELSE

                   vd_lu = 1.0_wp / (r_aero_surf + rb + rc_tot )

                   bud_lug(lsp) = - (conc_ijk(lsp) - ccomp_tot(lsp)) * &
                        (1.0_wp - exp(-vd_lu * dt_dh )) * dh
                ENDIF

             ENDIF
          ENDDO
       ENDIF
!
!--    Windows
       IF ( surf_usm_h%frac(m,ind_wat_win) > 0 )  THEN
!
!--       No vegetation on windows:
          lai = 0.0_wp
          sai = 0.0_wp

          slinnfac = 1.0_wp
!
!--       Get vd
          DO  lsp = 1, nvar
!
!--          Initialize
             vs = 0.0_wp
             vd_lu = 0.0_wp
             rs = 0.0_wp
             rb = 0.0_wp
             rc_tot = 0.0_wp 
             IF ( spc_names(lsp) == 'PM10' )  THEN
                part_type = 1
!
!--             Sedimentation velocity
                vs = slinnfac * sedimentation_velocity( particle_pars(ind_p_dens, part_type), &
                     particle_pars(ind_p_size, part_type), &
                     particle_pars(ind_p_slip, part_type), &
                     visc)

                CALL drydepo_aero_zhang_vd( vd_lu, rs, vs, &
                     particle_pars(ind_p_size, part_type), &
                     particle_pars(ind_p_slip, part_type), &
                     nwet, temp_tmp, dens, visc,              &
                     lud_dep, r_aero_surf, ustar_surf )

                bud_lud(lsp) = - conc_ijk(lsp) * &
                     (1.0_wp - exp(-vd_lu * dt_dh )) * dh

             ELSEIF ( spc_names(lsp) == 'PM25' )  THEN
                part_type = 2
!
!--             Sedimentation velocity
                vs = slinnfac * sedimentation_velocity( particle_pars(ind_p_dens, part_type), &
                     particle_pars(ind_p_size, part_type), &
                     particle_pars(ind_p_slip, part_type), &
                     visc)

                CALL drydepo_aero_zhang_vd( vd_lu, rs, vs, &
                     particle_pars(ind_p_size, part_type), &
                     particle_pars(ind_p_slip, part_type), &
                     nwet, temp_tmp, dens, visc, &
                     lud_dep, &
                     r_aero_surf, ustar_surf )

                bud_lud(lsp) = - conc_ijk(lsp) * &
                     (1.0_wp - exp(-vd_lu * dt_dh )) * dh

             ELSE  !< GASES
!
!--             Read spc_name of current species for gas PARAMETER

                IF ( ANY( pspecnames(:) == spc_names(lsp) ) )  THEN
                   i_pspec = 0
                   DO  pspec = 1, nposp
                      IF ( pspecnames(pspec) == spc_names(lsp) )  THEN
                         i_pspec = pspec
                      END IF
                   ENDDO
                ELSE
!
!--                For now species not deposited
                   CYCLE
                ENDIF
!
!--             Factor used for conversion from ppb to ug/m3 :
!--                 ppb (mole tr)/(mole air)/ppb (kg tr)/(mole tr) (ug tr)/(kg tr) &
!--                 (mole air)/(kg air) (kg air)/(m3 air) (kg air(ug/m3)/ppb/(kg/mole) = / (kg/mole)
!--                 c           1e-9               xm_tracer         1e9       /       xm_air            dens
!--             thus:
!--                  c_in_ppb * xm_tracer * [ dens / xm_air ] = c_in_ugm3
!--             Use density at k:

                ppm2ugm3 =  (dens/xm_air) * 0.001_wp  ! (mole air)/m3
!
!--             Atmospheric concentration in DEPAC is requested in ug/m3:
!--                 ug/m3              ppm          (ug/m3)/ppm/(kg/mole)     kg/mole
                conc_ijk_ugm3 = conc_ijk(lsp)         * ppm2ugm3 *   specmolm(i_pspec)  ! in ug/m3
!
!--             Diffusivity for DEPAC relevant gases
!--             Use default value 
                diffusivity = 0.11e-4
!
!--             overwrite with known coefficients of diffusivity from Massman (1998)
                IF ( spc_names(lsp) == 'NO2' ) diffusivity = 0.136e-4 
                IF ( spc_names(lsp) == 'NO'  ) diffusivity = 0.199e-4
                IF ( spc_names(lsp) == 'O3'  ) diffusivity = 0.144e-4
                IF ( spc_names(lsp) == 'CO'  ) diffusivity = 0.176e-4
                IF ( spc_names(lsp) == 'SO2' ) diffusivity = 0.112e-4
                IF ( spc_names(lsp) == 'CH4' ) diffusivity = 0.191e-4
                IF ( spc_names(lsp) == 'NH3' ) diffusivity = 0.191e-4
!
!--             Get quasi-laminar boundary layer resistance rb:
                CALL get_rb_cell( (lud_dep == ilu_water_sea) .OR. (lud_dep == ilu_water_inland),   &
                     z0h_surf, ustar_surf, diffusivity, rb )
!
!--             Get rc_tot
                CALL drydepos_gas_depac( spc_names(lsp), day_of_year, latitude, ts, ustar_surf,         &
                                         solar_rad, cos_zenith, rh_surf, lai, sai, nwet, lud_dep, 2,    &
                                         rc_tot, ccomp_tot(lsp), hyp(nzb), conc_ijk_ugm3, diffusivity,  &
                                         r_aero_surf , rb )
!
!--             Calculate budget 
                IF ( rc_tot <= 0.0 )  THEN

                   bud_lud(lsp) = 0.0_wp

                ELSE

                   vd_lu = 1.0_wp / (r_aero_surf + rb + rc_tot )

                   bud_lud(lsp) = - (conc_ijk(lsp) - ccomp_tot(lsp)) * &
                        (1.0_wp - exp(-vd_lu * dt_dh )) * dh
                ENDIF

             ENDIF
          ENDDO
       ENDIF


       bud = 0.0_wp
!
!--    Calculate overall budget for surface m and adapt concentration
       DO  lsp = 1, nspec


          bud(lsp) = surf_usm_h%frac(m,ind_veg_wall) * bud_luu(lsp) + &
               surf_usm_h%frac(m,ind_pav_green) * bud_lug(lsp) + &
               surf_usm_h%frac(m,ind_wat_win) * bud_lud(lsp)
!
!--       Compute new concentration
          conc_ijk(lsp) = conc_ijk(lsp) + bud(lsp) * inv_dh

          chem_species(lsp)%conc(k,j,i) = MAX( 0.0_wp, conc_ijk(lsp) )

       ENDDO

    ENDIF


 END SUBROUTINE chem_depo


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine to compute total canopy (or surface) resistance Rc for gases
!>
!> DEPAC:
!> Code of the DEPAC routine and corresponding subroutines below from the DEPAC
!> module of the LOTOS-EUROS model (Manders et al., 2017)
!>
!> Original DEPAC routines by RIVM and TNO (2015), for Documentation see
!> van Zanten et al., 2010.
!------------------------------------------------------------------------------!
 SUBROUTINE drydepos_gas_depac( compnam, day_of_year, lat, t, ust, solar_rad, sinphi,    &
      rh, lai, sai, nwet, lu, iratns, rc_tot, ccomp_tot, p, conc_ijk_ugm3, diffusivity,  &
      ra, rb )  
!
!--   Some of depac arguments are OPTIONAL:
!--    A. compute Rc_tot without compensation points (ccomp_tot will be zero):
!--        CALL depac (compnam, day_of_year, lat, t, ust, glrad, sinphi, rh, nwet, lu, iratns, rc_tot, ccomp_tot, [smi])
!--    B. compute Rc_tot with compensation points (used for LOTOS-EUROS):
!--        CALL depac (compnam, day_of_year, lat, t, ust, glrad, sinphi, rh, nwet, lu, iratns, rc_tot, ccomp_tot, [smi], &
!--                c_ave_prev_nh3, c_ave_prev_so2, catm, gamma_soil_water)
!--
!--    C. compute effective Rc based on compensation points (used for OPS):
!--        CALL depac (compnam, day_of_year, lat, t, ust, glrad, sinphi, rh, nwet, lu, iratns, rc_tot, ccomp_tot, [smi], &
!--               c_ave_prev_nh3, c_ave_prev_so2, catm, gamma_soil_water, &
!--               ra, rb, rc_eff)
!--    X1. Extra (OPTIONAL) output variables: 
!--        CALL depac (compnam, day_of_year, lat, t, ust, glrad, sinphi, rh, nwet, lu, iratns, rc_tot, ccomp_tot, [smi], &
!--               c_ave_prev_nh3, c_ave_prev_so2, catm, gamma_soil_water, &
!--               ra, rb, rc_eff, &
!--               gw_out, gstom_out, gsoil_eff_out, cw_out, cstom_out, csoil_out, lai_out, sai_out)
!--    X2. Extra (OPTIONAL) needed for stomatal ozone flux calculation (only sunlit leaves): 
!--        CALL depac (compnam, day_of_year, lat, t, ust, glrad, sinphi, rh, nwet, lu, iratns, rc_tot, ccomp_tot, [smi], &
!--               c_ave_prev_nh3, c_ave_prev_so2, catm, gamma_soil_water, &
!--               ra, rb, rc_eff, &
!--               gw_out, gstom_out, gsoil_eff_out, cw_out, cstom_out, csoil_out, lai_out, sai_out, &
!--               calc_stom_o3flux, frac_sto_o3_lu, fac_surface_area_2_PLA)


    CHARACTER(LEN=*), INTENT(IN) ::  compnam         !< component name
                                                     !< 'HNO3','NO','NO2','O3','SO2','NH3'
    INTEGER(iwp), INTENT(IN) ::  day_of_year         !< day of year, 1 ... 365 (366)
    INTEGER(iwp), INTENT(IN) ::  nwet                !< wetness indicator; nwet=0 -> dry; nwet=1 -> wet; nwet=9 -> snow
    INTEGER(iwp), INTENT(IN) ::  lu                  !< land use type, lu = 1,...,nlu
    INTEGER(iwp), INTENT(IN) ::  iratns              !< index for NH3/SO2 ratio used for SO2:
                                                     !< iratns = 1: low NH3/SO2
                                                     !< iratns = 2: high NH3/SO2
                                                     !< iratns = 3: very low NH3/SO2
    REAL(wp), INTENT(IN) ::  lat                     !< latitude Northern hemisphere (degrees) (S. hemisphere not possible)
    REAL(wp), INTENT(IN) ::  t                       !< temperature (C)
    REAL(wp), INTENT(IN) ::  ust                     !< friction velocity (m/s)
    REAL(wp), INTENT(IN) ::  solar_rad               !< solar radiation, dirict+diffuse (W/m2)
    REAL(wp), INTENT(IN) ::  sinphi                  !< sin of solar elevation angle
    REAL(wp), INTENT(IN) ::  rh                      !< relative humidity (%)
    REAL(wp), INTENT(IN) ::  lai                     !< one-sidedleaf area index (-)
    REAL(wp), INTENT(IN) ::  sai                     !< surface area index (-) (lai + branches and stems)
    REAL(wp), INTENT(IN) ::  diffusivity             !< diffusivity
    REAL(wp), INTENT(IN) ::  p                       !< pressure (Pa)
    REAL(wp), INTENT(IN) ::  conc_ijk_ugm3           !< actual atmospheric concentration (ug/m3), in DEPAC=Catm
    REAL(wp), INTENT(IN) ::  ra                      !< aerodynamic resistance (s/m)
    REAL(wp), INTENT(IN) ::  rb                      !< boundary layer resistance (s/m)

    REAL(wp), INTENT(OUT) ::  rc_tot                 !< total canopy resistance Rc (s/m)
    REAL(wp), INTENT(OUT) ::  ccomp_tot              !< total compensation point (ug/m3)
!                                                     !< [= 0 for species that don't have a compensation
!-- Local variables:
!
!-- Component number taken from component name, paramteres matched with include files
    INTEGER(iwp) ::  icmp
!
!-- Component numbers:
    INTEGER(iwp), PARAMETER ::  icmp_o3   = 1
    INTEGER(iwp), PARAMETER ::  icmp_so2  = 2
    INTEGER(iwp), PARAMETER ::  icmp_no2  = 3
    INTEGER(iwp), PARAMETER ::  icmp_no   = 4
    INTEGER(iwp), PARAMETER ::  icmp_nh3  = 5
    INTEGER(iwp), PARAMETER ::  icmp_co   = 6
    INTEGER(iwp), PARAMETER ::  icmp_no3  = 7
    INTEGER(iwp), PARAMETER ::  icmp_hno3 = 8
    INTEGER(iwp), PARAMETER ::  icmp_n2o5 = 9
    INTEGER(iwp), PARAMETER ::  icmp_h2o2 = 10

    LOGICAL ::  ready                                !< Rc has been set:
    !< = 1 -> constant Rc
    !< = 2 -> temperature dependent Rc
!
!-- Vegetation indicators:
    LOGICAL ::  lai_present                          !< leaves are present for current land use type
    LOGICAL ::  sai_present                          !< vegetation is present for current land use type

!    REAL(wp) ::  laimax                              !< maximum leaf area index (-)
    REAL(wp) ::  gw                                  !< external leaf conductance (m/s)
    REAL(wp) ::  gstom                               !< stomatal conductance (m/s)
    REAL(wp) ::  gsoil_eff                           !< effective soil conductance (m/s)
    REAL(wp) ::  gc_tot                              !< total canopy conductance (m/s)
    REAL(wp) ::  cw                                  !< external leaf surface compensation point (ug/m3)
    REAL(wp) ::  cstom                               !< stomatal compensation point (ug/m3)
    REAL(wp) ::  csoil                               !< soil compensation point (ug/m3)
!
!-- Next statement is just to avoid compiler warning about unused variable
    IF ( day_of_year == 0  .OR.  ( conc_ijk_ugm3 + lat + ra + rb ) > 0.0_wp )  CONTINUE
!
!-- Define component number
    SELECT CASE ( TRIM( compnam ) )

    CASE ( 'O3', 'o3' )
       icmp = icmp_o3

    CASE ( 'SO2', 'so2' )
       icmp = icmp_so2

    CASE ( 'NO2', 'no2' )
       icmp = icmp_no2

    CASE ( 'NO', 'no' )
       icmp = icmp_no 

    CASE ( 'NH3', 'nh3' )
       icmp = icmp_nh3

    CASE ( 'CO', 'co' )
       icmp = icmp_co

    CASE ( 'NO3', 'no3' )
       icmp = icmp_no3

    CASE ( 'HNO3', 'hno3' )
       icmp = icmp_hno3

    CASE ( 'N2O5', 'n2o5' )
       icmp = icmp_n2o5

    CASE ( 'H2O2', 'h2o2' )
       icmp = icmp_h2o2

    CASE default
!
!--    Component not part of DEPAC --> not deposited
       RETURN

    END SELECT

!
!-- Inititalize
    gw        = 0.0_wp
    gstom     = 0.0_wp
    gsoil_eff = 0.0_wp
    gc_tot    = 0.0_wp
    cw        = 0.0_wp
    cstom     = 0.0_wp
    csoil     = 0.0_wp
!
!-- Check whether vegetation is present:
    lai_present = ( lai > 0.0 )
    sai_present = ( sai > 0.0 )
!
!-- Set Rc (i.e. rc_tot) in special cases:
    CALL rc_special( icmp, compnam, lu, t, nwet, rc_tot, ready, ccomp_tot )
!
!-- If Rc is not set:
    IF ( .NOT. ready ) then
!
!--    External conductance:
       CALL rc_gw( compnam, iratns, t, rh, nwet, sai_present, sai,gw )         
!
!--    Stomatal conductance:
       CALL rc_gstom( icmp, compnam, lu, lai_present, lai, solar_rad, sinphi, t, rh, diffusivity, gstom, p )
!
!--    Effective soil conductance:
       CALL rc_gsoil_eff( icmp, lu, sai, ust, nwet, t, gsoil_eff )
!
!--    Total canopy conductance (gc_tot) and resistance Rc (rc_tot):
       CALL rc_rctot( gstom, gsoil_eff, gw, gc_tot, rc_tot )
!
!--    Needed to include compensation point for NH3
!--    Compensation points (always returns ccomp_tot; currently ccomp_tot != 0 only for NH3):
!--    CALL rc_comp_point( compnam,lu,day_of_year,t,gw,gstom,gsoil_eff,gc_tot,&
!--          lai_present, sai_present, &
!--          ccomp_tot, &
!--          conc_ijk_ugm3=conc_ijk_ugm3,c_ave_prev_nh3=c_ave_prev_nh3, &
!--          c_ave_prev_so2=c_ave_prev_so2,gamma_soil_water=gamma_soil_water, &
!--          tsea=tsea,cw=cw,cstom=cstom,csoil=csoil )
!
!--    Effective Rc based on compensation points:
!--        IF ( present(rc_eff) ) then
!--          check on required arguments:
!--           IF ( (.not. present(conc_ijk_ugm3)) .OR. (.not. present(ra)) .OR. (.not. present(rb)) ) then
!--              stop 'output argument rc_eff requires input arguments conc_ijk_ugm3, ra and rb'
!--           END IF
!
!--       Compute rc_eff :
       !      CALL rc_comp_point_rc_eff(ccomp_tot,conc_ijk_ugm3,ra,rb,rc_tot,rc_eff)
       !   ENDIF
    ENDIF

 END SUBROUTINE drydepos_gas_depac


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine to compute total canopy resistance in special cases
!------------------------------------------------------------------------------!
 SUBROUTINE rc_special( icmp, compnam, lu, t, nwet, rc_tot, ready, ccomp_tot )

   
    CHARACTER(LEN=*), INTENT(IN)  ::  compnam     !< component name

    INTEGER(iwp), INTENT(IN)  ::  icmp            !< component index      
    INTEGER(iwp), INTENT(IN)  ::  lu              !< land use type, lu = 1,...,nlu
    INTEGER(iwp), INTENT(IN)  ::  nwet            !< wetness indicator; nwet=0 -> dry; nwet=1 -> wet; nwet=9 -> snow

    REAL(wp), INTENT(IN)  ::  t                   !< temperature (C)

    REAL(wp), INTENT(OUT) ::  rc_tot             !< total canopy resistance Rc (s/m)
    REAL(wp), INTENT(OUT) ::  ccomp_tot          !< total compensation point (ug/m3)

    LOGICAL, INTENT(OUT) ::  ready               !< Rc has been set
                                                 !< = 1 -> constant Rc
!
!-- Next line is to avoid compiler warning about unused variable
    IF ( icmp == 0 )  CONTINUE
!
!-- rc_tot is not yet set:
    ready = .FALSE.
!
!-- Default compensation point in special CASEs = 0:
    ccomp_tot = 0.0_wp

    SELECT CASE( TRIM( compnam ) )
    CASE( 'HNO3', 'N2O5', 'NO3', 'H2O2' )
!
!--    No separate resistances for HNO3; just one total canopy resistance:
       IF ( t < -5.0_wp .AND. nwet == 9 )  THEN
!
!--       T < 5 C and snow:
          rc_tot = 50.0_wp
       ELSE
!
!--       all other circumstances:
          rc_tot = 10.0_wp
       ENDIF
       ready = .TRUE.

    CASE( 'NO', 'CO' )
       IF ( lu == ilu_water_sea .OR. lu == ilu_water_inland )  THEN       ! water
          rc_tot = 2000.0_wp
          ready = .TRUE.
       ELSEIF ( nwet == 1 )  THEN  !< wet
          rc_tot = 2000.0_wp
          ready = .TRUE.
       ENDIF
    CASE( 'NO2', 'O3', 'SO2', 'NH3' )
!
!--    snow surface:
       IF ( nwet == 9 )  THEN
!
!--       To be activated when snow is implemented
          !CALL rc_snow(ipar_snow(icmp),t,rc_tot)
          ready = .TRUE.
       ENDIF
    CASE default
       message_string = 'Component '// TRIM( compnam ) // ' not supported'
       CALL message( 'rc_special', 'CM0457', 1, 2, 0, 6, 0 )
    END SELECT

 END SUBROUTINE rc_special


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine to compute external conductance
!------------------------------------------------------------------------------!
 SUBROUTINE rc_gw( compnam, iratns, t, rh, nwet, sai_present, sai, gw )

!
!-- Input/output variables:
    CHARACTER(LEN=*), INTENT(IN) ::  compnam      !< component name

    INTEGER(iwp), INTENT(IN) ::  nwet             !< wetness indicator; nwet=0 -> dry; nwet=1 -> wet; nwet=9 -> snow
    INTEGER(iwp), INTENT(IN) ::  iratns           !< index for NH3/SO2 ratio;
                                                  !< iratns = 1: low NH3/SO2
                                                  !< iratns = 2: high NH3/SO2
                                                  !< iratns = 3: very low NH3/SO2
    LOGICAL, INTENT(IN) ::  sai_present

    REAL(wp), INTENT(IN) ::  t                    !< temperature (C)
    REAL(wp), INTENT(IN) ::  rh                   !< relative humidity (%)
    REAL(wp), INTENT(IN) ::  sai                  !< one-sided leaf area index (-)

    REAL(wp), INTENT(OUT) ::  gw                  !< external leaf conductance (m/s)

    SELECT CASE( TRIM( compnam ) )

    CASE( 'NO2' )
       CALL rw_constant( 2000.0_wp, sai_present, gw )

    CASE( 'NO', 'CO' )
       CALL rw_constant( -9999.0_wp, sai_present, gw )   !< see Erisman et al, 1994 section 3.2.3

    CASE( 'O3' )
       CALL rw_constant( 2500.0_wp, sai_present, gw )

    CASE( 'SO2' )
       CALL rw_so2( t, nwet, rh, iratns, sai_present, gw )

    CASE( 'NH3' )
       CALL rw_nh3_sutton( t, rh, sai_present, gw )
!
!--    conversion from leaf resistance to canopy resistance by multiplying with sai:
       gw = sai * gw

    CASE default
       message_string = 'Component '// TRIM( compnam ) // ' not supported'
       CALL message( 'rc_gw', 'CM0458', 1, 2, 0, 6, 0 )
    END SELECT

 END SUBROUTINE rc_gw


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine to compute external leaf conductance for SO2
!------------------------------------------------------------------------------!
 SUBROUTINE rw_so2( t, nwet, rh, iratns, sai_present, gw )

!
!-- Input/output variables:
    INTEGER(iwp), INTENT(IN) ::  nwet        !< wetness indicator; nwet=0 -> dry; nwet=1 -> wet; nwet=9 -> snow
    INTEGER(iwp), INTENT(IN) ::  iratns      !< index for NH3/SO2 ratio:
                                             !< iratns = 1: low NH3/SO2
                                             !< iratns = 2: high NH3/SO2
                                             !< iratns = 3: very low NH3/SO2
    LOGICAL, INTENT(IN) ::  sai_present

    REAL(wp), INTENT(IN) ::  t               !< temperature (C)
    REAL(wp), INTENT(IN) ::  rh              !< relative humidity (%)   

    REAL(wp), INTENT(OUT) ::  gw             !< external leaf conductance (m/s)
!
!-- Local variables:
    REAL(wp) ::  rw                          !< external leaf resistance (s/m)
!
!-- Check if vegetation present:
    IF ( sai_present )  THEN

       IF ( nwet == 0 )  THEN
!
!--   ------------------------
!--         dry surface
!--   ------------------------
!--         T > -1 C
          IF ( t > -1.0_wp )  THEN
             IF ( rh < 81.3_wp )  THEN
                rw = 25000.0_wp * exp( -0.0693_wp * rh )
             ELSE
                rw = 0.58e12 * exp( -0.278_wp * rh ) + 10.0_wp
             ENDIF
          ELSE
             ! -5 C < T <= -1 C
             IF ( t > -5.0_wp )  THEN
                rw = 200.0_wp
             ELSE
                ! T <= -5 C
                rw = 500.0_wp
             ENDIF
          ENDIF
       ELSE
!
!--   ------------------------
!--         wet surface
!--   ------------------------
          rw = 10.0_wp !see Table 5, Erisman et al, 1994 Atm. Environment, 0 is impl. as 10
       ENDIF
!
!--    very low NH3/SO2 ratio:
       IF ( iratns == iratns_very_low ) rw = rw + 50.0_wp
!
!--      Conductance:
       gw = 1.0_wp / rw
    ELSE
!
!--      no vegetation:
       gw = 0.0_wp
    ENDIF

 END SUBROUTINE rw_so2


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine to compute external leaf conductance for NH3,
!>                  following Sutton & Fowler, 1993
!------------------------------------------------------------------------------!
 SUBROUTINE rw_nh3_sutton( tsurf, rh,sai_present, gw )

!
!-- Input/output variables:
    LOGICAL, INTENT(IN) ::  sai_present

    REAL(wp), INTENT(IN) ::  tsurf          !< surface temperature (C)
    REAL(wp), INTENT(IN) ::  rh             !< relative humidity (%)

    REAL(wp), INTENT(OUT) ::  gw            !< external leaf conductance (m/s)
!
!-- Local variables:
    REAL(wp) ::  rw                         !< external leaf resistance (s/m)
    REAL(wp) ::  sai_grass_haarweg          !< surface area index at experimental site Haarweg
!
!-- Fix sai_grass at value valid for Haarweg data for which gamma_w parametrization is derived
    sai_grass_haarweg = 3.5_wp
!
!-- Calculation rw:
!--                    100 - rh
!--    rw = 2.0 * exp(----------)
!--                      12

    IF ( sai_present )  THEN
!
!--    External resistance according to Sutton & Fowler, 1993
       rw = 2.0_wp * exp( ( 100.0_wp - rh ) / 12.0_wp )
       rw = sai_grass_haarweg * rw
!
!--    Frozen soil (from Depac v1):
       IF ( tsurf < 0.0_wp ) rw = 200.0_wp
!
!--    Conductance:
       gw = 1.0_wp / rw
    ELSE
       ! no vegetation:
       gw = 0.0_wp
    ENDIF

 END SUBROUTINE rw_nh3_sutton


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine to compute external leaf conductance 
!------------------------------------------------------------------------------!
 SUBROUTINE rw_constant( rw_val, sai_present, gw )

!
!-- Input/output variables:
    LOGICAL, INTENT(IN) ::  sai_present

    REAL(wp), INTENT(IN) ::  rw_val       !< constant value of Rw    

    REAL(wp), INTENT(OUT) ::  gw          !< wernal leaf conductance (m/s)
!
!-- Compute conductance:
    IF ( sai_present .AND. .NOT.missing(rw_val) )  THEN
       gw = 1.0_wp / rw_val
    ELSE
       gw = 0.0_wp
    ENDIF

 END SUBROUTINE rw_constant


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine to compute stomatal conductance 
!------------------------------------------------------------------------------!
 SUBROUTINE rc_gstom( icmp, compnam, lu, lai_present, lai, solar_rad, sinphi, t, rh, diffusivity, gstom, p ) 

!
!-- input/output variables:
    CHARACTER(LEN=*), INTENT(IN) ::  compnam       !< component name

    INTEGER(iwp), INTENT(IN) ::  icmp              !< component index
    INTEGER(iwp), INTENT(IN) ::  lu                !< land use type , lu = 1,...,nlu

    LOGICAL, INTENT(IN) ::  lai_present

    REAL(wp), INTENT(IN) ::  lai                   !< one-sided leaf area index
    REAL(wp), INTENT(IN) ::  solar_rad             !< solar radiation, dirict+diffuse (W/m2)
    REAL(wp), INTENT(IN) ::  sinphi                !< sin of solar elevation angle
    REAL(wp), INTENT(IN) ::  t                     !< temperature (C)
    REAL(wp), INTENT(IN) ::  rh                    !< relative humidity (%)
    REAL(wp), INTENT(IN) ::  diffusivity           !< diffusion coefficient of the gas involved

    REAL(wp), OPTIONAL,INTENT(IN) :: p             !< pressure (Pa)

    REAL(wp), INTENT(OUT) ::  gstom                !< stomatal conductance (m/s)
!
!-- Local variables
    REAL(wp) ::  vpd                               !< vapour pressure deficit (kPa)

    REAL(wp), PARAMETER ::  dO3 = 0.13e-4          !< diffusion coefficient of ozon (m2/s)
!
!-- Next line is to avoid compiler warning about unused variables
    IF ( icmp == 0 )  CONTINUE

    SELECT CASE( TRIM( compnam ) )

    CASE( 'NO', 'CO' )
!
!--    For no stomatal uptake is neglected:
       gstom = 0.0_wp

    CASE( 'NO2', 'O3', 'SO2', 'NH3' )
!
!--    if vegetation present:
       IF ( lai_present )  THEN

          IF ( solar_rad > 0.0_wp )  THEN
             CALL rc_get_vpd( t, rh, vpd )
             CALL rc_gstom_emb( lu, solar_rad, t, vpd, lai_present, lai, sinphi, gstom, p )
             gstom = gstom * diffusivity / dO3       !< Gstom of Emberson is derived for ozone 
          ELSE
             gstom = 0.0_wp
          ENDIF
       ELSE
!
!--       no vegetation; zero conductance (infinite resistance):
          gstom = 0.0_wp
       ENDIF

    CASE default
       message_string = 'Component '// TRIM( compnam ) // ' not supported'
       CALL message( 'rc_gstom', 'CM0459', 1, 2, 0, 6, 0 )
    END SELECT

 END SUBROUTINE rc_gstom


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine to compute stomatal conductance according to Emberson
!------------------------------------------------------------------------------!
 SUBROUTINE rc_gstom_emb( lu, solar_rad, T, vpd, lai_present, lai, sinp, Gsto, p )
!
!>  History
!>   Original code from Lotos-Euros, TNO, M. Schaap
!>   2009-08, M.C. van Zanten, Rivm
!>     Updated and extended.
!>   2009-09, Arjo Segers, TNO
!>     Limitted temperature influence to range to avoid
!>     floating point exceptions.

!> Method

!>   Code based on Emberson et al, 2000, Env. Poll., 403-413
!>   Notation conform Unified EMEP Model Description Part 1, ch 8
!
!>   In the calculation of f_light the modification of L. Zhang 2001, AE to the PARshade and PARsun
!>   parametrizations of Norman 1982 are applied
!>   f_phen and f_SWP are set to 1
!
!>   Land use types DEPAC versus Emberson (Table 5.1, EMEP model description)
!>   DEPAC                     Emberson
!>     1 = grass                 GR = grassland
!>     2 = arable land           TC = temperate crops ( lai according to RC = rootcrops)
!>     3 = permanent crops       TC = temperate crops ( lai according to RC = rootcrops)
!>     4 = coniferous forest     CF = tempareate/boREAL(wp) coniferous forest
!>     5 = deciduous forest      DF = temperate/boREAL(wp) deciduous forest
!>     6 = water                 W  = water
!>     7 = urban                 U  = urban
!>     8 = other                 GR = grassland
!>     9 = desert                DE = desert
!
!-- Emberson specific declarations
!
!-- Input/output variables:
    INTEGER(iwp), INTENT(IN) ::  lu             !< land use type, lu = 1,...,nlu

    LOGICAL, INTENT(IN) ::  lai_present

    REAL(wp), INTENT(IN) ::  solar_rad          !< solar radiation, dirict+diffuse (W/m2)
    REAL(wp), INTENT(IN) ::  t                  !< temperature (C)
    REAL(wp), INTENT(IN) ::  vpd                !< vapour pressure deficit (kPa)

    REAL(wp), INTENT(IN) ::  lai                !< one-sided leaf area index
    REAL(wp), INTENT(IN) ::  sinp               !< sin of solar elevation angle

    REAL(wp), OPTIONAL, INTENT(IN) ::  p        !< pressure (Pa)

    REAL(wp), INTENT(OUT) ::  gsto              !< stomatal conductance (m/s)
!
!-- Local variables:
    REAL(wp) ::  f_light
    REAL(wp) ::  f_phen
    REAL(wp) ::  f_temp
    REAL(wp) ::  f_vpd
    REAL(wp) ::  f_swp
    REAL(wp) ::  bt
    REAL(wp) ::  f_env
    REAL(wp) ::  pardir
    REAL(wp) ::  pardiff
    REAL(wp) ::  parshade
    REAL(wp) ::  parsun
    REAL(wp) ::  laisun
    REAL(wp) ::  laishade
    REAL(wp) ::  sinphi
    REAL(wp) ::  pres
    REAL(wp), PARAMETER ::  p_sealevel = 1.01325e05    !< Pa
!
!-- Check whether vegetation is present:
    IF ( lai_present )  THEN

       ! calculation of correction factors for stomatal conductance
       IF ( sinp <= 0.0_wp )  THEN  
          sinphi = 0.0001_wp
       ELSE
          sinphi = sinp
       END IF
!
!--    ratio between actual and sea-level pressure is used
!--    to correct for height in the computation of par;
!--    should not exceed sea-level pressure therefore ...
       IF (  present(p) )  THEN
          pres = min( p, p_sealevel )
       ELSE 
          pres = p_sealevel
       ENDIF
!
!--    direct and diffuse par, Photoactive (=visible) radiation:
       CALL par_dir_diff( solar_rad, sinphi, pres, p_sealevel, pardir, pardiff )
!
!--    par for shaded leaves (canopy averaged):
       parshade = pardiff * exp( -0.5 * lai**0.7 ) + 0.07 * pardir * ( 1.1 - 0.1 * lai ) * exp( -sinphi )     !< Norman,1982
       IF ( solar_rad > 200.0_wp .AND. lai > 2.5_wp )  THEN
          parshade = pardiff * exp( -0.5 * lai**0.8 ) + 0.07 * pardir * ( 1.1 - 0.1 * lai ) * exp( -sinphi )  !< Zhang et al., 2001
       END IF
!
!--    par for sunlit leaves (canopy averaged):
!--    alpha -> mean angle between leaves and the sun is fixed at 60 deg -> i.e. cos alpha = 0.5
       parsun = pardir * 0.5/sinphi + parshade             !< Norman, 1982
       IF ( solar_rad > 200.0_wp .AND. lai > 2.5_wp )  THEN
          parsun = pardir**0.8 * 0.5 / sinphi + parshade   !< Zhang et al., 2001
       END IF
!
!--    leaf area index for sunlit and shaded leaves:
       IF ( sinphi > 0 )  THEN
          laisun = 2 * sinphi * ( 1 - exp( -0.5 * lai / sinphi ) )
          laishade = lai - laisun
       ELSE
          laisun = 0
          laishade = lai
       END IF

       f_light = ( laisun * ( 1 - exp( -1.0_wp * alpha(lu) * parsun ) ) + &
            laishade * ( 1 - exp( -1.0_wp * alpha(lu) * parshade ) ) ) / lai 

       f_light = MAX(f_light,f_min(lu))
!
!--    temperature influence; only non-zero within range [tmin,tmax]:
       IF ( ( tmin(lu) < t ) .AND. ( t < tmax(lu) ) )  THEN
          bt = ( tmax(lu) - topt(lu) ) / ( topt(lu) - tmin(lu) )
          f_temp = ( ( t - tmin(lu) ) / ( topt(lu) - tmin(lu) ) ) * ( ( tmax(lu) - t ) / ( tmax(lu) - topt(lu) ) )**bt
       ELSE
          f_temp = 0.0_wp
       END IF
       f_temp = MAX( f_temp, f_min(lu) )
!
!--      vapour pressure deficit influence
       f_vpd = MIN( 1.0_wp, ( ( 1.0_wp - f_min(lu) ) * ( vpd_min(lu) - vpd ) / ( vpd_min(lu) - vpd_max(lu) ) + f_min(lu) ) )
       f_vpd = MAX( f_vpd, f_min(lu) )

       f_swp = 1.0_wp
!
!--      influence of phenology on stom. conductance
!--      ignored for now in DEPAC since influence of f_phen on lu classes in use is negligible.
!--      When other EMEP classes (e.g. med. broadleaf) are used f_phen might be too important to ignore
       f_phen = 1.0_wp
!
!--      evaluate total stomatal conductance
       f_env = f_temp * f_vpd * f_swp
       f_env = MAX( f_env,f_min(lu) )
       gsto = g_max(lu) * f_light * f_phen * f_env
!
!--      gstom expressed per m2 leafarea;
!--      this is converted with lai to m2 surface.
       gsto = lai * gsto    ! in m/s

    ELSE
       gsto = 0.0_wp
    ENDIF

 END SUBROUTINE rc_gstom_emb


 !-------------------------------------------------------------------
 !> par_dir_diff
 !>     Weiss, A., Norman, J.M. (1985) Partitioning solar radiation into direct and
 !>     diffuse, visible and near-infrared components. Agric. Forest Meteorol.
 !>     34, 205-213.
 !>     From a SUBROUTINE obtained from Leiming Zhang,
 !>     Meteorological Service of Canada
 !>     Leiming uses solar irradiance. This should be equal to global radiation and
 !>     Willem Asman set it to global radiation (here defined as solar radiation, dirict+diffuse)
 !>
 !>     @todo Check/connect/replace with radiation_model_mod variables    
 !-------------------------------------------------------------------
 SUBROUTINE par_dir_diff( solar_rad, sinphi, pres, pres_0, par_dir, par_diff )


    REAL(wp), INTENT(IN) ::  solar_rad       !< solar radiation, dirict+diffuse (W m-2)
    REAL(wp), INTENT(IN) ::  sinphi          !< sine of the solar elevation
    REAL(wp), INTENT(IN) ::  pres            !< actual pressure (to correct for height) (Pa)
    REAL(wp), INTENT(IN) ::  pres_0          !< pressure at sea level (Pa)

    REAL(wp), INTENT(OUT) ::  par_dir        !< par direct : visible (photoactive) direct beam radiation (W m-2)
    REAL(wp), INTENT(OUT) ::  par_diff       !< par diffuse: visible (photoactive) diffuse radiation (W m-2)


    REAL(wp) ::  sv                          !< total visible radiation
    REAL(wp) ::  fv                          !< par direct beam fraction (dimensionless)
    REAL(wp) ::  ratio                       !< ratio measured to potential solar radiation (dimensionless)
    REAL(wp) ::  rdm                         !< potential direct beam near-infrared radiation (W m-2); "potential" means clear-sky
    REAL(wp) ::  rdn                         !< potential diffuse near-infrared radiation (W m-2)
    REAL(wp) ::  rdu                         !< visible (par) direct beam radiation (W m-2)
    REAL(wp) ::  rdv                         !< potential visible (par) diffuse radiation (W m-2)
    REAL(wp) ::  rn                          !< near-infrared radiation (W m-2)
    REAL(wp) ::  rv                          !< visible radiation (W m-2)
    REAL(wp) ::  ww                          !< water absorption in the near infrared for 10 mm of precipitable water

!
!-- Calculate visible (PAR) direct beam radiation
!-- 600 W m-2 represents average amount of par (400-700 nm wavelength)
!-- at the top of the atmosphere; this is roughly 0.45*solar constant (solar constant=1320 Wm-2)
    rdu = 600.0_wp* exp( -0.185_wp * ( pres / pres_0 ) / sinphi ) * sinphi
!
!-- Calculate potential visible diffuse radiation
    rdv = 0.4_wp * ( 600.0_wp - rdu ) * sinphi
!
!-- Calculate the water absorption in the-near infrared
    ww = 1320 * 10**( -1.195_wp + 0.4459_wp * log10( 1.0_wp / sinphi ) - 0.0345_wp * ( log10( 1.0_wp / sinphi ) )**2 )
!
!-- Calculate potential direct beam near-infrared radiation
    rdm = (720.0_wp * exp(-0.06_wp * (pres / pres_0) / sinphi ) - ww ) * sinphi     !< 720 = solar constant - 600
!
!-- Calculate potential diffuse near-infrared radiation
    rdn = 0.6_wp * ( 720 - rdm - ww ) * sinphi
!
!-- Compute visible and near-infrared radiation
    rv = MAX( 0.1_wp, rdu + rdv )
    rn = MAX( 0.01_wp, rdm + rdn )
!
!-- Compute ratio between input global radiation (here defined as solar radiation, dirict+diffuse)
!-- and total radiation computed here
    ratio = MIN( 0.89_wp, solar_rad / ( rv + rn ) )
!
!-- Calculate total visible radiation
    sv = ratio * rv
!
!-- Calculate fraction of par in the direct beam
    fv = MIN( 0.99_wp, ( 0.9_wp - ratio ) / 0.7_wp )              !< help variable
    fv = MAX( 0.01_wp, rdu / rv * ( 1.0_wp - fv**0.6667_wp ) )    !< fraction of par in the direct beam
!
!-- Compute direct and diffuse parts of par
    par_dir = fv * sv
    par_diff = sv - par_dir

 END SUBROUTINE par_dir_diff

 
 !-------------------------------------------------------------------
 !> rc_get_vpd: get vapour pressure deficit (kPa)
 !-------------------------------------------------------------------
 SUBROUTINE rc_get_vpd( temp, rh, vpd )

!
!-- Input/output variables:
    REAL(wp), INTENT(IN) ::  temp    !< temperature (C)
    REAL(wp), INTENT(IN) ::  rh    !< relative humidity (%)

    REAL(wp), INTENT(OUT) ::  vpd    !< vapour pressure deficit (kPa)
!
!-- Local variables:
    REAL(wp) ::  esat
!
!-- fit parameters:
    REAL(wp), PARAMETER ::  a1 = 6.113718e-01
    REAL(wp), PARAMETER ::  a2 = 4.43839e-02
    REAL(wp), PARAMETER ::  a3 = 1.39817e-03
    REAL(wp), PARAMETER ::  a4 = 2.9295e-05
    REAL(wp), PARAMETER ::  a5 = 2.16e-07
    REAL(wp), PARAMETER ::  a6 = 3.0e-09
!
!-- esat is saturation vapour pressure (kPa) at temp(C) following Monteith(1973)
    esat = a1 + a2 * temp + a3 * temp**2 + a4 * temp**3 + a5 * temp**4 + a6 * temp**5
    vpd  = esat * ( 1 - rh / 100 )

 END SUBROUTINE rc_get_vpd


 !-------------------------------------------------------------------
 !> rc_gsoil_eff: compute effective soil conductance
 !-------------------------------------------------------------------
 SUBROUTINE rc_gsoil_eff( icmp, lu, sai, ust, nwet, t, gsoil_eff )

!
!-- Input/output variables:
    INTEGER(iwp), INTENT(IN) ::  icmp          !< component index
    INTEGER(iwp), INTENT(IN) ::  lu            !< land use type, lu = 1,..., nlu
    INTEGER(iwp), INTENT(IN) ::  nwet          !< index for wetness
                                               !< nwet = 0 -> dry; nwet = 1 -> wet; nwet = 9 -> snow
                                               !< N.B. this routine cannot be called with nwet = 9,
                                               !< nwet = 9 should be handled outside this routine.
    REAL(wp), INTENT(IN) ::  sai               !< surface area index
    REAL(wp), INTENT(IN) ::  ust               !< friction velocity (m/s)
    REAL(wp), INTENT(IN) ::  t                 !< temperature (C)
    REAL(wp), INTENT(OUT) ::  gsoil_eff        !< effective soil conductance (m/s)
!
!-- local variables:
    REAL(wp) ::  rinc                          !< in canopy resistance  (s/m)
    REAL(wp) ::  rsoil_eff                     !< effective soil resistance (s/m)
!
!-- Soil resistance (numbers matched with lu_classes and component numbers)
    !     grs    ara    crp    cnf    dec    wat    urb   oth    des    ice    sav    trf    wai    med    sem
    REAL(wp), PARAMETER ::  rsoil(nlu_dep,ncmp) = reshape( (/ &
         1000.,  200.,  200.,  200.,  200., 2000.,  400., 1000., 2000., 2000., 1000.,  200., 2000.,  200.,  400., &    !< O3
         1000., 1000., 1000., 1000., 1000.,   10., 1000., 1000., 1000.,  500., 1000., 1000.,   10., 1000., 1000., &    !< SO2
         1000., 1000., 1000., 1000., 1000., 2000., 1000., 1000., 1000., 2000., 1000., 1000., 2000., 1000., 1000., &    !< NO2
         -999., -999., -999., -999., -999., 2000., 1000., -999., 2000., 2000., -999., -999., 2000., -999., -999., &    !< NO
         100.,  100.,  100.,  100.,  100.,   10.,  100.,  100.,  100., 1000.,  100.,  100.,   10.,  100.,  100.,  &    !< NH3
         -999., -999., -999., -999., -999., 2000., 1000., -999., 2000., 2000., -999., -999., 2000., -999., -999., &    !< CO 
         -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., &    !< NO3 
         -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., &    !< HNO3
         -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., &    !< N2O5
         -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999. /),&  !< H2O2   
         (/nlu_dep,ncmp/) )
!
!-- For                                          o3    so2   no2     no    nh3     co     no3    hno3   n2o5   h2o2
    REAL(wp), PARAMETER ::  rsoil_wet(ncmp)    = (/2000., 10. , 2000., -999., 10.  , -999., -999., -999., -999., -999./)
    REAL(wp), PARAMETER ::  rsoil_frozen(ncmp) = (/2000., 500., 2000., -999., 1000., -999., -999., -999., -999., -999./)
!
!-- Compute in canopy (in crop) resistance:
    CALL rc_rinc( lu, sai, ust, rinc )
!
!-- Check for missing deposition path:
    IF ( missing(rinc) )  THEN
       rsoil_eff = -9999.0_wp
    ELSE
!
!--    Frozen soil (temperature below 0):
       IF ( t < 0.0_wp )  THEN
          IF ( missing( rsoil_frozen( icmp ) ) )  THEN
             rsoil_eff = -9999.0_wp
          ELSE
             rsoil_eff = rsoil_frozen( icmp ) + rinc
          ENDIF
       ELSE
!
!--       Non-frozen soil; dry:
          IF ( nwet == 0 )  THEN
             IF ( missing( rsoil( lu, icmp ) ) )  THEN
                rsoil_eff = -9999.0_wp
             ELSE
                rsoil_eff = rsoil( lu, icmp ) + rinc
             ENDIF
!
!--       Non-frozen soil; wet:
          ELSEIF ( nwet == 1 )  THEN
             IF ( missing( rsoil_wet( icmp ) ) )  THEN
                rsoil_eff = -9999.0_wp
             ELSE
                rsoil_eff = rsoil_wet( icmp ) + rinc
             ENDIF
          ELSE
             message_string = 'nwet can only be 0 or 1'
             CALL message( 'rc_gsoil_eff', 'CM0460', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF
    ENDIF
!
!-- Compute conductance:
    IF ( rsoil_eff > 0.0_wp )  THEN
       gsoil_eff = 1.0_wp / rsoil_eff
    ELSE
       gsoil_eff = 0.0_wp
    ENDIF

 END SUBROUTINE rc_gsoil_eff


 !-------------------------------------------------------------------
 !> rc_rinc: compute in canopy (or in crop) resistance
 !> van Pul and Jacobs, 1993, BLM
 !-------------------------------------------------------------------
 SUBROUTINE rc_rinc( lu, sai, ust, rinc )

!
!-- Input/output variables:
    INTEGER(iwp), INTENT(IN) ::  lu          !< land use class, lu = 1, ..., nlu

    REAL(wp), INTENT(IN) ::  sai             !< surface area index
    REAL(wp), INTENT(IN) ::  ust             !< friction velocity (m/s)

    REAL(wp), INTENT(OUT) ::  rinc           !< in canopy resistance (s/m)
!
!-- b = empirical constant for computation of rinc (in canopy resistance) (= 14 m-1 or -999 if not applicable)
!-- h = vegetation height (m)                     gra  ara crop con dec wat   urb   oth   des   ice   sav   trf  wai  med semi
    REAL(wp), DIMENSION(nlu_dep), PARAMETER ::  b = (/ -999, 14, 14, 14, 14, -999, -999, -999, -999, -999, -999, 14, -999,  &
         14, 14 /)
    REAL(wp), DIMENSION(nlu_dep), PARAMETER ::  h = (/ -999, 1,  1,  20, 20, -999, -999, -999, -999, -999, -999, 20, -999,  &
         1 ,  1 /)
!
!-- Compute Rinc only for arable land, perm. crops, forest; otherwise Rinc = 0:
    IF ( b(lu) > 0.0_wp )  THEN
!       !
!--    Check for u* > 0 (otherwise denominator = 0):
       IF ( ust > 0.0_wp )  THEN
          rinc = b(lu) * h(lu) * sai/ust
       ELSE
          rinc = 1000.0_wp
       ENDIF
    ELSE
       IF ( lu == ilu_grass .OR. lu == ilu_other  )  THEN
          rinc = -999.0_wp     !< no deposition path for grass, other, and semi-natural
       ELSE
          rinc = 0.0_wp        !< no in-canopy resistance
       ENDIF
    ENDIF

 END SUBROUTINE rc_rinc


 !-------------------------------------------------------------------
 !> rc_rctot: compute total canopy (or surface) resistance Rc
 !-------------------------------------------------------------------
 SUBROUTINE rc_rctot( gstom, gsoil_eff, gw, gc_tot, rc_tot )

!
!-- Input/output variables:
    REAL(wp), INTENT(IN) ::  gstom         !< stomatal conductance (s/m)
    REAL(wp), INTENT(IN) ::  gsoil_eff     !< effective soil conductance (s/m)
    REAL(wp), INTENT(IN) ::  gw            !< external leaf conductance (s/m)

    REAL(wp), INTENT(OUT) ::  gc_tot       !< total canopy conductance (m/s)
    REAL(wp), INTENT(OUT) ::  rc_tot       !< total canopy resistance Rc (s/m)
!
!-- Total conductance:
    gc_tot = gstom + gsoil_eff + gw
!
!-- Total resistance (note: gw can be negative, but no total emission allowed here):
    IF ( gc_tot <= 0.0_wp .OR. gw < 0.0_wp )  THEN
       rc_tot = -9999.0_wp
    ELSE
       rc_tot = 1.0_wp / gc_tot
    ENDIF

 END SUBROUTINE rc_rctot


 !-------------------------------------------------------------------
 !> rc_comp_point_rc_eff: calculate the effective resistance Rc
 !> based on one or more compensation points
 !-------------------------------------------------------------------
 !> NH3rc (see depac v3.6 is based on Avero workshop Marc Sutton. p. 173.
 !> Sutton 1998 AE 473-480)
 !>
 !> Documentation by Ferd Sauter, 2008; see also documentation block in header of depac subroutine.
 !> FS 2009-01-29: variable names made consistent with DEPAC 
 !> FS 2009-03-04: use total compensation point
 !>
 !> C: with total compensation point   ! D: approximation of C
 !>                                    !    with classical approach
 !>  zr --------- Catm                 !  zr --------- Catm    
 !>         |                          !         |       
 !>         Ra                         !         Ra      
 !>         |                          !         |       
 !>         Rb                         !         Rb      
 !>         |                          !         |       
 !>  z0 --------- Cc                   !  z0 --------- Cc
 !>         |                          !         |              
 !>        Rc                          !        Rc_eff          
 !>         |                          !         |              
 !>     --------- Ccomp_tot            !     --------- C=0
 !>
 !>
 !> The effective Rc is defined such that instead of using
 !>
 !>   F = -vd*[Catm - Ccomp_tot]                                    (1)
 !>
 !> we can use the 'normal' flux formula
 !>
 !>   F = -vd'*Catm,                                                (2)
 !>
 !> with vd' = 1/(Ra + Rb + Rc')                                    (3)
 !>
 !> and Rc' the effective Rc (rc_eff).
 !>                                                (Catm - Ccomp_tot)
 !> vd'*Catm = vd*(Catm - Ccomp_tot) <=> vd' = vd* ------------------
 !>                                                      Catm
 !>
 !>                                        (Catm - Ccomp_tot)
 !> 1/(Ra + Rb + Rc') = (1/Ra + Rb + Rc) * ------------------
 !>                                              Catm
 !>
 !>                                          Catm
 !> (Ra + Rb + Rc') = (Ra + Rb + Rc) * ------------------
 !>                                     (Catm - Ccomp_tot)
 !>
 !>                              Catm
 !> Rc' = (Ra + Rb + Rc) * ------------------ - Ra - Rb
 !>                        (Catm - Ccomp_tot)
 !>
 !>                        Catm                           Catm
 !> Rc' = (Ra + Rb) [------------------ - 1 ] + Rc * ------------------
 !>                  (Catm - Ccomp_tot)              (Catm - Ccomp_tot)
 !>
 !> Rc' = [(Ra + Rb)*Ccomp_tot + Rc*Catm ] / (Catm - Ccomp_tot)
 !>
 ! -------------------------------------------------------------------------------------------
! SUBROUTINE rc_comp_point_rc_eff( ccomp_tot, conc_ijk_ugm3, ra, rb, rc_tot, rc_eff )
!
!
!!-- Input/output variables:
!    REAL(wp), INTENT(IN) ::  ccomp_tot     !< total compensation point (weighed average of separate compensation points) (ug/m3)
!    REAL(wp), INTENT(IN) ::  conc_ijk_ugm3 !< atmospheric concentration (ug/m3) above Catm
!    REAL(wp), INTENT(IN) ::  ra            !< aerodynamic resistance (s/m)
!    REAL(wp), INTENT(IN) ::  rb            !< boundary layer resistance (s/m)
!    REAL(wp), INTENT(IN) ::  rc_tot        !< total canopy resistance (s/m)
!
!    REAL(wp), INTENT(OUT) ::  rc_eff       !< effective total canopy resistance (s/m)
!
!    !
!!-- Compute effective resistance:
!    IF (  ccomp_tot == 0.0_wp )  THEN
!       !
!!--    trace with no compensiation point ( or compensation point equal to zero)
!       rc_eff = rc_tot
!
!    ELSE IF ( ccomp_tot > 0.0_wp .AND. ( abs( conc_ijk_ugm3 - ccomp_tot ) < 1.e-8 ) )  THEN
!       !
!!--   surface concentration (almost) equal to atmospheric concentration
!!--    no exchange between surface and atmosphere, infinite RC --> vd=0
!       rc_eff = 9999999999.0_wp
!
!    ELSE IF ( ccomp_tot > 0.0_wp )  THEN
!       !
!!--    compensation point available, calculate effective resistance
!       rc_eff = ( ( ra + rb ) * ccomp_tot + rc_tot * conc_ijk_ugm3 ) / ( conc_ijk_ugm3 - ccomp_tot )
!
!    ELSE
!       rc_eff = -999.0_wp
!       message_string = 'This should not be possible, check ccomp_tot'
!       CALL message( 'rc_comp_point_rc_eff', 'CM0461', 1, 2, 0, 6, 0 )
!    ENDIF
!
!    RETURN
!    
! END SUBROUTINE rc_comp_point_rc_eff


 !-------------------------------------------------------------------
 !> missing: check for data that correspond with a missing deposition path
 !>          this data is represented by -999
 !-------------------------------------------------------------------
 LOGICAL function missing( x )

    REAL(wp), INTENT(IN) ::  x

!
!-- bandwidth for checking (in)equalities of floats
    REAL(wp), PARAMETER :: eps = 1.0e-5

    missing = (abs(x + 999.0_wp) <= eps)

 END function missing


 ELEMENTAL FUNCTION sedimentation_velocity( rhopart, partsize, slipcor, visc ) RESULT( vs )

!
!-- in/out

    REAL(wp), INTENT(IN) ::  rhopart                 !< particle density (kg/m3)
    REAL(wp), INTENT(IN) ::  partsize                !< particle size (m)
    REAL(wp), INTENT(IN) ::  slipcor                 !< slip correction factor (m)
    REAL(wp), INTENT(IN) ::  visc                    !< viscosity

    REAL(wp) ::  vs
!
!-- acceleration of gravity:
    REAL(wp), PARAMETER         ::  grav = 9.80665_wp   !< m/s2

!-- sedimentation velocity
    vs = rhopart * ( partsize**2 ) * grav * slipcor / ( 18.0_wp * visc )

 END FUNCTION sedimentation_velocity

 
 !------------------------------------------------------------------------
 !> Boundary-layer deposition resistance following Zhang (2001)
 !------------------------------------------------------------------------
 SUBROUTINE drydepo_aero_zhang_vd( vd, rs, vs1, partsize, slipcor, nwet, tsurf, dens1, viscos1, &
      luc, ftop_lu, ustar )

!
!-- in/out 

    INTEGER(iwp), INTENT(IN) ::  nwet        !< 1=rain, 9=snowcover
    INTEGER(iwp), INTENT(IN) ::  luc         !< DEPAC LU

    REAL(wp), INTENT(IN) ::  vs1             !< sedimentation velocity in lowest layer
    REAL(wp), INTENT(IN) ::  partsize        !< particle diameter (m)
    REAL(wp), INTENT(IN) ::  slipcor         !< slip correction factor
    REAL(wp), INTENT(IN) ::  tsurf           !< surface temperature (K)
    REAL(wp), INTENT(IN) ::  dens1           !< air density (kg/m3) in lowest layer
    REAL(wp), INTENT(IN) ::  viscos1         !< air viscosity in lowest layer
    REAL(wp), INTENT(IN) ::  ftop_lu         !< atmospheric resistnace Ra
    REAL(wp), INTENT(IN) ::  ustar           !< friction velocity u*    

    REAL(wp), INTENT(OUT) ::  vd             !< deposition velocity (m/s)
    REAL(wp), INTENT(OUT) ::  rs             !< sedimentaion resistance (s/m)
!
!-- constants 

    REAL(wp), PARAMETER ::  grav     = 9.80665_wp             !< acceleration of gravity (m/s2)

    REAL(wp), PARAMETER ::  epsilon0 = 3.0_wp
    REAL(wp), PARAMETER ::  kb       = 1.38066E-23_wp
    REAL(wp), PARAMETER ::  pi       = 3.141592654_wp      !< pi

    REAL(wp), PARAMETER :: alfa_lu(nlu_dep) = & 
         (/1.2_wp,  1.2_wp,   1.2_wp,  1.0_wp,  1.0_wp,   100.0_wp, 1.5_wp,  1.2_wp, 50.0_wp, 100.0_wp, &
               1.2_wp, 1.0_wp, 100.0_wp, 1.2_wp, 50.0_wp/)   
    REAL(wp), PARAMETER :: gamma_lu(nlu_dep) = &
         (/0.54_wp, 0.54_wp,  0.54_wp, 0.56_wp, 0.56_wp,  0.50_wp,  0.56_wp, 0.54_wp, 0.58_wp, 0.50_wp, &
              0.54_wp, 0.56_wp, 0.50_wp, 0.54_wp, 0.54_wp/)   
    REAL(wp), PARAMETER ::A_lu(nlu_dep) = &   
         (/3.0_wp,  3.0_wp,   2.0_wp,  2.0_wp,  7.0_wp, -99.0_wp, 10.0_wp, 3.0_wp, -99.0_wp, -99.0_wp,  &
              3.0_wp, 7.0_wp, -99.0_wp, 2.0_wp, -99.0_wp/)
!
!--   grass  arabl crops conif decid  water  urba  othr  desr  ice   sav  trf   wai  med   sem     
!
!-- local
    REAL(wp) ::  kinvisc
    REAL(wp) ::  diff_part
    REAL(wp) ::  schmidt
    REAL(wp) ::  stokes
    REAL(wp) ::  Ebrown
    REAL(wp) ::  Eimpac
    REAL(wp) ::  Einterc
    REAL(wp) ::  Reffic
!
!-- kinetic viscosity & diffusivity
    kinvisc = viscos1 / dens1    !< only needed at surface

    diff_part = kb * tsurf * slipcor / ( 3.0_wp * pi * viscos1 * partsize )
!
!-- Schmidt number
    schmidt = kinvisc / diff_part
!
!-- calculate collection efficiencie E
    Ebrown = Schmidt**( -gamma_lu(luc) )    !< Brownian diffusion
!
!-- determine Stokes number, interception efficiency 
!-- and sticking efficiency R (1 = no rebound)
    IF ( luc == ilu_ice .OR. nwet==9 .OR. luc == ilu_water_sea .OR. luc == ilu_water_inland )  THEN
       stokes = vs1 * ustar**2 / ( grav * kinvisc )
       Einterc = 0.0_wp
       Reffic = 1.0_wp
    ELSE IF ( luc == ilu_other .OR. luc == ilu_desert )  THEN     !<tundra of desert
       stokes = vs1 * ustar**2 / ( grav * kinvisc )
       Einterc = 0.0_wp
       Reffic = exp( -Stokes**0.5_wp )
    ELSE
       stokes = vs1 * ustar / ( grav * A_lu(luc) * 1.0E-3_wp )
       Einterc = 0.5_wp * ( partsize / (A_lu(luc) * 1.0E-3_wp ) )**2
       Reffic = exp( -Stokes**0.5_wp )
    END IF
!
!-- when surface is wet all particles do not rebound:
    IF ( nwet==1 )  Reffic = 1.0_wp
!
!-- determine impaction efficiency:
    Eimpac = ( stokes / ( alfa_lu(luc) + stokes ) )**2
!
!-- sedimentation resistance:
    rs = 1.0_wp / ( epsilon0 * MAX( 1.0E-5_wp, ustar ) * ( Ebrown + Eimpac + Einterc ) * Reffic )

!-- deposition velocity according to Seinfeld and Pandis (2006; eq 19.7):
!--  
!--              1
!--      vd = ------------------ + vs
!--           Ra + Rs + Ra*Rs*vs
!--  
!-- where: Rs = Rb (in Seinfeld and Pandis, 2006)

    vd = 1.0_wp / ( ftop_lu + rs + ftop_lu * rs * vs1) + vs1


 END SUBROUTINE drydepo_aero_zhang_vd


 !-------------------------------------------------------------------------------------
 !> Compute quasi-laminar boundary layer resistance as a function of landuse and tracer
 !> Original EMEP formulation by (Simpson et al, 2003) is used 
 !-------------------------------------------------------------------------------------
 SUBROUTINE get_rb_cell( is_water, z0h, ustar, diffusivity, rb )    

!
!-- in/out 

    LOGICAL , INTENT(IN) ::  is_water

    REAL(wp), INTENT(IN) ::  z0h                  !< roughness length for heat
    REAL(wp), INTENT(IN) ::  ustar                !< friction velocity
    REAL(wp), INTENT(IN) ::  diffusivity          !< coefficient of diffusivity

    REAL(wp), INTENT(OUT) ::  rb                  !< boundary layer resistance
!
!-- const

    REAL(wp), PARAMETER ::  thk = 0.19e-4         !< thermal diffusivity of dry air 20 C
    REAL(wp), PARAMETER ::  kappa_stab = 0.35     !< von Karman constant
!
!-- Next line is to avoid compiler warning about unused variable
    IF ( is_water  .OR.  ( z0h + kappa_stab ) > 0.0_wp )  CONTINUE
!
!-- Use Simpson et al. (2003) 
!-- @TODO: Check rb over water calculation, until then leave commented lines
!--  IF ( is_water )  THEN
!--   org: rb = 1.0_wp / (kappa_stab*MAX(0.01_wp,ustar)) * log(z0h/diffusivity*kappa_stab*MAX(0.01_wp,ustar))
!--        rb = 1.0_wp / (kappa_stab*MAX(0.1_wp,ustar)) * log(z0h/diffusivity*kappa_stab*MAX(0.1_wp,ustar))
!--  ELSE
    rb = 5.0_wp / MAX( 0.01_wp, ustar ) * ( thk / diffusivity )**0.67_wp
!--  END IF

 END SUBROUTINE get_rb_cell


 !-----------------------------------------------------------------
 !>  Compute water vapor partial pressure (e_w)
 !>  given specific humidity Q [(kg water)/(kg air)].
 !>
 !>  Use that gas law for volume V with temperature T
 !>  holds for the total mixture as well as the water part:
 !>
 !>    R T / V = p_air / n_air = p_water / n_water
 !>
 !>  thus:
 !>
 !>    p_water = p_air n_water / n_air
 !>
 !>  Use:
 !>    n_air =   m_air   /        xm_air
 !>            [kg air]  /  [(kg air)/(mole air)]
 !>  and:
 !>    n_water =  m_air * Q  /     xm_water
 !>              [kg water]  /  [(kg water)/(mole water)]
 !>  thus:
 !>    p_water = p_air Q / (xm_water/xm_air)
 !------------------------------------------------------------------

 ELEMENTAL FUNCTION watervaporpartialpressure( q, p ) RESULT( p_w )

!
!-- in/out 

    REAL(wp), INTENT(IN) ::  q                      !< specific humidity [(kg water)/(kg air)]
    REAL(wp), INTENT(IN) ::  p                      !< air pressure [Pa]

    REAL(wp) ::  p_w                                !< water vapor partial pressure [Pa]
!
!-- const 

    REAL(wp), PARAMETER  ::  eps = xm_h2o / xm_air  !< mole mass ratio ~ 0.622
!
!-- partial pressure of water vapor:
    p_w = p * q / eps

 END function watervaporpartialpressure


 !------------------------------------------------------------------    
 !>  Saturation vapor pressure.
 !>  From (Stull 1988, eq. 7.5.2d):
 !>
 !>      e_sat = p0 exp( 17.67 * (T-273.16) / (T-29.66) )     [Pa]
 !>
 !>  where:
 !>      p0 = 611.2 [Pa]   : reference pressure
 !>
 !>  Arguments:
 !>      T  [K]  : air temperature
 !>  Result:
 !>      e_sat_w  [Pa]  : saturation vapor pressure
 !>
 !>  References:
 !>      Roland B. Stull, 1988
 !>      An introduction to boundary layer meteorology.
 !-----------------------------------------------------------------

 ELEMENTAL FUNCTION saturationvaporpressure( t ) RESULT( e_sat_w )

!
!-- in/out

    REAL(wp), INTENT(IN) ::  t            !< temperature [K]

    REAL(wp) ::  e_sat_w                  !< saturation vapor pressure  [Pa]
!
!-- const 
    REAL(wp), PARAMETER ::  p0 = 611.2   !< base pressure [Pa]
!
!-- saturation vapor pressure:
    e_sat_w = p0 * exp( 17.67_wp * ( t - 273.16_wp ) / ( t - 29.66_wp ) )    !< [Pa]

 END FUNCTION saturationvaporpressure


 !------------------------------------------------------------------------
 !>  Relative humidity RH [%] is by definition:
 !>
 !>           e_w             water vapor partial pressure
 !>    Rh = -------- * 100
 !>         e_sat_w           saturation vapor pressure
 !------------------------------------------------------------------------

 ELEMENTAL FUNCTION relativehumidity_from_specifichumidity( q, t, p ) RESULT( rh )

!
!-- in/out 

    REAL(wp), INTENT(IN) ::  q    !< specific humidity [(kg water)/(kg air)]
    REAL(wp), INTENT(IN) ::  t    !< temperature [K]
    REAL(wp), INTENT(IN) ::  p    !< air pressure [Pa]

    REAL(wp) ::  rh               !< relative humidity [%]
!
!-- relative humidity:
    rh = watervaporpartialpressure( q, p ) / saturationvaporpressure( t ) * 100.0_wp

 END FUNCTION relativehumidity_from_specifichumidity

     
 END MODULE chemistry_model_mod
