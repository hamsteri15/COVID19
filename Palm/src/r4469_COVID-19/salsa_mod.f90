!> @file salsa_mod.f90
!--------------------------------------------------------------------------------!
! This file is part of PALM-4U.
!
! PALM-4U is free software: you can redistribute it and/or modify it under the 
! terms of the GNU General Public License as published by the Free Software 
! Foundation, either version 3 of the License, or (at your option) any later 
! version.
!
! PALM-4U is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
! A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with
! PALM. If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 2018-2019 University of Helsinki
! Copyright 1997-2020 Leibniz Universitaet Hannover
!--------------------------------------------------------------------------------!
!
! Current revisions:
! -----------------
! 
! 
! Former revisions:
! -----------------
! $Id: salsa_mod.f90 4457 2020-03-11 14:20:43Z raasch $
! use statement for exchange horiz added
! 
! 4442 2020-03-04 19:21:13Z suehring
! Change order of dimension in surface array %frac to allow for better 
! vectorization.
! 
! 4441 2020-03-04 19:20:35Z suehring
! Bug fixes and reformatting for the restart data and averaged data output
! - add missing arrays (averaged data output) in salsa_wrd_local and
!   salsa_rrd_local
! - set write_binary_salsa and read_restart_data_salsa to .T. by default
! - restructure the average arrays for gases and total mass concentrations of
!   chemical components: set to 4d arrays instead of separate arrays
! - add allocation checks for averaged data output arrays
! 
! 4416 2020-02-20 17:53:57Z monakurppa
! Time index error in salsa_emission_setup
! 
! 4380 2020-01-17 23:39:51Z monakurppa
! - Error in saving the surface fluxes in an array that is applied in the
!   deposition scheme
! - Corrections in the header: aerosol bin diameters and lower bin limits not
!   printed correctly
! 
! 4364 2020-01-08 02:12:31Z monakurppa
! Set time coordinate in the input data relative to origin_time rather than to
! 00:00:00 UTC
! 
! 4360 2020-01-07 11:25:50Z suehring
! Introduction of wall_flags_total_0, which currently sets bits based on static
! topography information used in wall_flags_static_0
! 
! 4342 2019-12-16 13:49:14Z Giersch
! cdc replaced by canopy_drag_coeff
! 
! 4329 2019-12-10 15:46:36Z motisi
! Renamed wall_flags_0 to wall_flags_static_0
! 
! 4315 2019-12-02 09:20:07Z monakurppa
! Add an additional check for the time dimension PIDS_SALSA in
! salsa_emission_setup and correct some error message identifiers.
! 
! 4298 2019-11-21 15:59:16Z suehring
! Bugfix, close netcdf input files after reading
! 
! 4295 2019-11-14 06:15:31Z monakurppa
! 
! 
! 4280 2019-10-29 14:34:15Z monakurppa
! Corrected a bug in boundary conditions and fac_dt in offline nesting
! 
! 4273 2019-10-24 13:40:54Z monakurppa
! - Rename nest_salsa to nesting_salsa
! - Correct some errors in boundary condition flags
! - Add a check for not trying to output gas concentrations in salsa if the
!   chemistry module is applied
! - Set the default value of nesting_salsa and nesting_offline_salsa to .TRUE.
! 
! 4272 2019-10-23 15:18:57Z schwenkel
! Further modularization of boundary conditions: moved boundary conditions to
! respective modules
!
! 4270 2019-10-23 10:46:20Z monakurppa
! - Implement offline nesting for salsa
! - Alphabetic ordering for module interfaces
! - Remove init_aerosol_type and init_gases_type from salsa_parin and define them
!   based on the initializing_actions
! - parameter definition removed from "season" and "season_z01" is added to parin
! - bugfix in application of index_hh after implementing the new
!   palm_date_time_mod
! - Reformat salsa emission data with LOD=2: size distribution given for each
!   emission category
! 
! 4268 2019-10-17 11:29:38Z schwenkel
! Moving module specific boundary conditions from time_integration to module
! 
! 4256 2019-10-07 10:08:52Z monakurppa
! Document previous changes: use global variables nx, ny and nz in salsa_header
! 
! 4227 2019-09-10 18:04:34Z gronemeier
! implement new palm_date_time_mod
! 
! 4226 2019-09-10 17:03:24Z suehring
! Netcdf input routine for dimension length renamed
! 
! 4182 2019-08-22 15:20:23Z scharf
! Corrected "Former revisions" section
! 
! 4167 2019-08-16 11:01:48Z suehring
! Changed behaviour of masked output over surface to follow terrain and ignore 
! buildings (J.Resler, T.Gronemeier)
! 
! 4131 2019-08-02 11:06:18Z monakurppa
! - Add "salsa_" before each salsa output variable
! - Add a possibility to output the number (salsa_N_UFP) and mass concentration
!   (salsa_PM0.1) of ultrafine particles, i.e. particles with a diameter smaller
!   than 100 nm
! - Implement aerosol emission mode "parameterized" which is based on the street
!   type (similar to the chemistry module).
! - Remove unnecessary nucleation subroutines.
! - Add the z-dimension for gaseous emissions to correspond the implementation
!   in the chemistry module
! 
! 4118 2019-07-25 16:11:45Z suehring
! - When Dirichlet condition is applied in decycling, the boundary conditions are
!   only set at the ghost points and not at the prognostic grid points as done
!   before
! - Rename decycle_ns/lr to decycle_salsa_ns/lr and decycle_method to
!   decycle_method_salsa
! - Allocation and initialization of special advection flags salsa_advc_flags_s
!   used for salsa. These are exclusively used for salsa variables to
!   distinguish from the usually-used flags which might be different when
!   decycling is applied in combination with cyclic boundary conditions.
!   Moreover, salsa_advc_flags_s considers extended zones around buildings where
!   the first-order upwind scheme is applied for the horizontal advection terms.
!   This is done to overcome high concentration peaks due to stationary numerical
!   oscillations caused by horizontal advection discretization.
! 
! 4117 2019-07-25 08:54:02Z monakurppa
! Pass integer flag array as well as boundary flags to WS scalar advection 
! routine
! 
! 4109 2019-07-22 17:00:34Z suehring
! Slightly revise setting of boundary conditions at horizontal walls, use 
! data-structure offset index instead of pre-calculate it for each facing
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
! 4058 2019-06-27 15:25:42Z knoop
! Bugfix: to_be_resorted was uninitialized in case of s_H2O in 3d_data_averaging
! 
! 4012 2019-05-31 15:19:05Z monakurppa
! Merge salsa branch to trunk. List of changes:
! - Error corrected in distr_update that resulted in the aerosol number size
!   distribution not converging if the concentration was nclim.
! - Added a separate output for aerosol liquid water (s_H2O)
! - aerosol processes for a size bin are now calculated only if the aerosol
!   number of concentration of that bin is > 2*nclim
! - An initialisation error in the subroutine "deposition" corrected and the
!   subroutine reformatted.
! - stuff from salsa_util_mod.f90 moved into salsa_mod.f90
! - calls for closing the netcdf input files added
! 
! 3956 2019-05-07 12:32:52Z monakurppa
! - Conceptual bug in depo_surf correct for urban and land surface model
! - Subroutine salsa_tendency_ij optimized.
! - Interfaces salsa_non_advective_processes and salsa_exchange_horiz_bounds
!   created. These are now called in module_interface.
!   salsa_exchange_horiz_bounds after calling salsa_driver only when needed
!   (i.e. every dt_salsa).
! 
! 3924 2019-04-23 09:33:06Z monakurppa
! Correct a bug introduced by the previous update.
! 
! 3899 2019-04-16 14:05:27Z monakurppa
! - remove unnecessary error / location messages
! - corrected some error message numbers
! - allocate source arrays only if emissions or dry deposition is applied.
!
! 3885 2019-04-11 11:29:34Z kanani
! Changes related to global restructuring of location messages and introduction 
! of additional debug messages
! 
! 3876 2019-04-08 18:41:49Z knoop
! Introduced salsa_actions module interface
! 
! 3871 2019-04-08 14:38:39Z knoop
! Major changes in formatting, performance and data input structure (see branch
! the history for details)
! - Time-dependent emissions enabled: lod=1 for yearly PM emissions that are
!   normalised depending on the time, and lod=2 for preprocessed emissions
!   (similar to the chemistry module).
! - Additionally, 'uniform' emissions allowed. This emission is set constant on
!   all horisontal upward facing surfaces and it is created based on parameters
!   surface_aerosol_flux, aerosol_flux_dpg/sigmag/mass_fracs_a/mass_fracs_b.
! - All emissions are now implemented as surface fluxes! No 3D sources anymore.
! - Update the emission information by calling salsa_emission_update if 
!   skip_time_do_salsa >= time_since_reference_point and
!   next_aero_emission_update <= time_since_reference_point
! - Aerosol background concentrations read from PIDS_DYNAMIC. The vertical grid
!   must match the one applied in the model.
! - Gas emissions and background concentrations can be also read in in salsa_mod
!   if the chemistry module is not applied.
! - In deposition, information on the land use type can be now imported from
!   the land use model
! - Use SI units in PARIN, i.e. n_lognorm given in #/m3 and dpg in metres.
! - Apply 100 character line limit
! - Change all variable names from capital to lowercase letter
! - Change real exponents to integer if possible. If not, precalculate the value
!   value of exponent
! - Rename in1a to start_subrange_1a, fn2a to end_subrange_1a etc.
! - Rename nbins --> nbins_aerosol, ncc_tot --> ncomponents_mass and ngast -->
!   ngases_salsa
! - Rename ibc to index_bc, idu to index_du etc.
! - Renamed loop indices b, c and sg to ib, ic and ig
! - run_salsa subroutine removed
! - Corrected a bud in salsa_driver: falsely applied ino instead of inh
! - Call salsa_tendency within salsa_prognostic_equations which is called in
!   module_interface_mod instead of prognostic_equations_mod
! - Removed tailing white spaces and unused variables
! - Change error message to start by PA instead of SA
! 
! 3833 2019-03-28 15:04:04Z forkel
! added USE chem_gasphase_mod for nvar, nspec and spc_names 
! 
! 3787 2019-03-07 08:43:54Z raasch
! unused variables removed
! 
! 3780 2019-03-05 11:19:45Z forkel
! unused variable for file index removed from rrd-subroutines parameter list
! 
! 3685 2019-01-21 01:02:11Z knoop
! Some interface calls moved to module_interface + cleanup
! 
! 3655 2019-01-07 16:51:22Z knoop
! Implementation of the PALM module interface
! 3412 2018-10-24 07:25:57Z monakurppa
!
! Authors:
! --------
! @author Mona Kurppa (University of Helsinki)
!
!
! Description:
! ------------
!> Sectional aerosol module for large scale applications SALSA
!> (Kokkola et al., 2008, ACP 8, 2469-2483). Solves the aerosol number and mass
!> concentration as well as chemical composition. Includes aerosol dynamic
!> processes: nucleation, condensation/evaporation of vapours, coagulation and
!> deposition on tree leaves, ground and roofs.
!> Implementation is based on formulations implemented in UCLALES-SALSA except
!> for deposition which is based on parametrisations by Zhang et al. (2001,
!> Atmos. Environ. 35, 549-560) or Petroff&Zhang (2010, Geosci. Model Dev. 3,
!> 753-769)
!>
!> @todo Apply information from emission_stack_height to lift emission sources
!> @todo Allow insoluble emissions
!------------------------------------------------------------------------------!
 MODULE salsa_mod

    USE basic_constants_and_equations_mod,                                                         &
        ONLY:  c_p, g, p_0, pi, r_d

    USE chem_gasphase_mod,                                                                         &
        ONLY:  nspec, nvar, spc_names

    USE chem_modules,                                                                              &
        ONLY:  call_chem_at_all_substeps, chem_gasphase_on, chem_species

    USE control_parameters,                                                                        &
        ONLY:  air_chemistry, bc_dirichlet_l, bc_dirichlet_n, bc_dirichlet_r, bc_dirichlet_s,      &
               bc_lr, bc_lr_cyc, bc_ns, bc_ns_cyc, bc_radiation_l, bc_radiation_n, bc_radiation_r, &
               bc_radiation_s, coupling_char, debug_output, dt_3d, intermediate_timestep_count,    &
               intermediate_timestep_count_max, land_surface, max_pr_salsa, message_string,        &
               monotonic_limiter_z, plant_canopy, pt_surface, salsa, scalar_advec,                 &
               surface_pressure, time_since_reference_point, timestep_scheme, tsc, urban_surface,  &
               ws_scheme_sca

    USE indices,                                                                                   &
        ONLY:  nbgp, nx, nxl, nxlg, nxr, nxrg, ny, nyn, nyng, nys, nysg, nzb, nz, nzt,             &
               wall_flags_total_0

    USE kinds

    USE netcdf_data_input_mod,                                                                     &
        ONLY:  chem_emis_att_type, chem_emis_val_type

    USE pegrid

    USE statistics,                                                                                &
        ONLY:  sums_salsa_ws_l

    IMPLICIT NONE
!
!-- SALSA constants:
!
!-- Local constants:
    INTEGER(iwp), PARAMETER ::  luc_urban = 15     !< default landuse type for urban
    INTEGER(iwp), PARAMETER ::  ngases_salsa  = 5  !< total number of gaseous tracers:
                                                   !< 1 = H2SO4, 2 = HNO3, 3 = NH3, 4 = OCNV
                                                   !< (non-volatile OC), 5 = OCSV (semi-volatile)
    INTEGER(iwp), PARAMETER ::  nmod = 7     !< number of modes for initialising the aerosol size distribution
    INTEGER(iwp), PARAMETER ::  nreg = 2     !< Number of main size subranges
    INTEGER(iwp), PARAMETER ::  maxspec = 7  !< Max. number of aerosol species


    REAL(wp), PARAMETER ::  fill_value = -9999.0_wp    !< value for the _FillValue attribute
!
!-- Universal constants
    REAL(wp), PARAMETER ::  abo    = 1.380662E-23_wp   !< Boltzmann constant (J/K)
    REAL(wp), PARAMETER ::  alv    = 2.260E+6_wp       !< latent heat for H2O vaporisation (J/kg)
    REAL(wp), PARAMETER ::  alv_d_rv  = 4896.96865_wp  !< alv / rv
    REAL(wp), PARAMETER ::  am_airmol = 4.8096E-26_wp  !< Average mass of an air molecule (Jacobson 2005, Eq.2.3)
    REAL(wp), PARAMETER ::  api6   = 0.5235988_wp      !< pi / 6
    REAL(wp), PARAMETER ::  argas  = 8.314409_wp       !< Gas constant (J/(mol K))
    REAL(wp), PARAMETER ::  argas_d_cpd = 8.281283865E-3_wp  !< argas per cpd
    REAL(wp), PARAMETER ::  avo    = 6.02214E+23_wp    !< Avogadro constant (1/mol)
    REAL(wp), PARAMETER ::  d_sa   = 5.539376964394570E-10_wp  !< diameter of condensing H2SO4 molecule (m)
    REAL(wp), PARAMETER ::  for_ppm_to_nconc =  7.243016311E+16_wp !< ppm * avo / R (K/(Pa*m3))
    REAL(wp), PARAMETER ::  epsoc  = 0.15_wp          !< water uptake of organic material
    REAL(wp), PARAMETER ::  mclim  = 1.0E-23_wp       !< mass concentration min limit (kg/m3)
    REAL(wp), PARAMETER ::  n3     = 158.79_wp        !< Number of H2SO4 molecules in 3 nm cluster if d_sa=5.54e-10m
    REAL(wp), PARAMETER ::  nclim  = 1.0_wp           !< number concentration min limit (#/m3)
    REAL(wp), PARAMETER ::  surfw0 = 0.073_wp         !< surface tension of water at 293 K (J/m2)
!
!-- Molar masses in kg/mol
    REAL(wp), PARAMETER ::  ambc     = 12.0E-3_wp     !< black carbon (BC)
    REAL(wp), PARAMETER ::  amdair   = 28.970E-3_wp   !< dry air
    REAL(wp), PARAMETER ::  amdu     = 100.0E-3_wp    !< mineral dust
    REAL(wp), PARAMETER ::  amh2o    = 18.0154E-3_wp  !< H2O
    REAL(wp), PARAMETER ::  amh2so4  = 98.06E-3_wp    !< H2SO4
    REAL(wp), PARAMETER ::  amhno3   = 63.01E-3_wp    !< HNO3
    REAL(wp), PARAMETER ::  amn2o    = 44.013E-3_wp   !< N2O
    REAL(wp), PARAMETER ::  amnh3    = 17.031E-3_wp   !< NH3
    REAL(wp), PARAMETER ::  amo2     = 31.9988E-3_wp  !< O2
    REAL(wp), PARAMETER ::  amo3     = 47.998E-3_wp   !< O3
    REAL(wp), PARAMETER ::  amoc     = 150.0E-3_wp    !< organic carbon (OC)
    REAL(wp), PARAMETER ::  amss     = 58.44E-3_wp    !< sea salt (NaCl)
!
!-- Densities in kg/m3
    REAL(wp), PARAMETER ::  arhobc     = 2000.0_wp  !< black carbon
    REAL(wp), PARAMETER ::  arhodu     = 2650.0_wp  !< mineral dust
    REAL(wp), PARAMETER ::  arhoh2o    = 1000.0_wp  !< H2O
    REAL(wp), PARAMETER ::  arhoh2so4  = 1830.0_wp  !< SO4
    REAL(wp), PARAMETER ::  arhohno3   = 1479.0_wp  !< HNO3
    REAL(wp), PARAMETER ::  arhonh3    = 1530.0_wp  !< NH3
    REAL(wp), PARAMETER ::  arhooc     = 2000.0_wp  !< organic carbon
    REAL(wp), PARAMETER ::  arhoss     = 2165.0_wp  !< sea salt (NaCl)
!
!-- Volume of molecule in m3/#
    REAL(wp), PARAMETER ::  amvh2o   = amh2o /avo / arhoh2o      !< H2O
    REAL(wp), PARAMETER ::  amvh2so4 = amh2so4 / avo / arhoh2so4 !< SO4
    REAL(wp), PARAMETER ::  amvhno3  = amhno3 / avo / arhohno3   !< HNO3
    REAL(wp), PARAMETER ::  amvnh3   = amnh3 / avo / arhonh3     !< NH3
    REAL(wp), PARAMETER ::  amvoc    = amoc / avo / arhooc       !< OC
    REAL(wp), PARAMETER ::  amvss    = amss / avo / arhoss       !< sea salt
!
!-- Constants for the dry deposition model by Petroff and Zhang (2010):
!-- obstacle characteristic dimension "L" (cm) (plane obstacle by default) and empirical constants
!-- C_B, C_IN, C_IM, beta_IM and C_IT for each land use category (15, as in Zhang et al. (2001))
    REAL(wp), DIMENSION(1:15), PARAMETER :: l_p10 = &
        (/0.15,  4.0,   0.15,  3.0,   3.0,   0.5,   3.0,   -99., 0.5,  2.0,  1.0,   -99., -99., -99., 3.0/)
    REAL(wp), DIMENSION(1:15), PARAMETER :: c_b_p10 = &
        (/0.887, 1.262, 0.887, 1.262, 1.262, 0.996, 0.996, -99., 0.7,  0.93, 0.996, -99., -99., -99., 1.262/)
    REAL(wp), DIMENSION(1:15), PARAMETER :: c_in_p10 = &
        (/0.81,  0.216, 0.81,  0.216, 0.216, 0.191, 0.162, -99., 0.7,  0.14, 0.162, -99., -99., -99., 0.216/)
    REAL(wp), DIMENSION(1:15), PARAMETER :: c_im_p10 = &
        (/0.162, 0.13,  0.162, 0.13,  0.13,  0.191, 0.081, -99., 0.191,0.086,0.081, -99., -99., -99., 0.13/)
    REAL(wp), DIMENSION(1:15), PARAMETER :: beta_im_p10 = &
        (/0.6,   0.47,  0.6,   0.47,  0.47,  0.47,  0.47,  -99., 0.6,  0.47, 0.47,  -99., -99., -99., 0.47/)
    REAL(wp), DIMENSION(1:15), PARAMETER :: c_it_p10 = &
        (/0.0,   0.056, 0.0,   0.056, 0.056, 0.042, 0.056, -99., 0.042,0.014,0.056, -99., -99., -99., 0.056/)
!
!-- Constants for the dry deposition model by Zhang et al. (2001): 
!-- empirical constants "alpha" and "gamma" and characteristic radius "A" for 
!-- each land use category (15) and season (5)
    REAL(wp), DIMENSION(1:15), PARAMETER :: alpha_z01 = &
        (/1.0, 0.6, 1.1, 0.8, 0.8, 1.2, 1.2, 50.0, 50.0, 1.3, 2.0, 50.0, 100.0, 100.0, 1.5/)
    REAL(wp), DIMENSION(1:15), PARAMETER :: gamma_z01 = &
        (/0.56, 0.58, 0.56, 0.56, 0.56, 0.54, 0.54, 0.54, 0.54, 0.54, 0.54, 0.54, 0.50, 0.50, 0.56/)
    REAL(wp), DIMENSION(1:15,1:5), PARAMETER :: A_z01 =  RESHAPE( (/& 
         2.0, 5.0, 2.0,  5.0, 5.0, 2.0, 2.0, -99., -99., 10.0, 10.0, -99., -99., -99., 10.0,&  ! SC1
         2.0, 5.0, 2.0,  5.0, 5.0, 2.0, 2.0, -99., -99., 10.0, 10.0, -99., -99., -99., 10.0,&  ! SC2
         2.0, 5.0, 5.0, 10.0, 5.0, 5.0, 5.0, -99., -99., 10.0, 10.0, -99., -99., -99., 10.0,&  ! SC3
         2.0, 5.0, 5.0, 10.0, 5.0, 5.0, 5.0, -99., -99., 10.0, 10.0, -99., -99., -99., 10.0,&  ! SC4
         2.0, 5.0, 2.0,  5.0, 5.0, 2.0, 2.0, -99., -99., 10.0, 10.0, -99., -99., -99., 10.0 &  ! SC5
                                                           /), (/ 15, 5 /) )
!-- Land use categories (based on Z01 but the same applies here also for P10):
!-- 1 = evergreen needleleaf trees,
!-- 2 = evergreen broadleaf trees,
!-- 3 = deciduous needleleaf trees,
!-- 4 = deciduous broadleaf trees,
!-- 5 = mixed broadleaf and needleleaf trees (deciduous broadleaf trees for P10),
!-- 6 = grass (short grass for P10),
!-- 7 = crops, mixed farming,
!-- 8 = desert,
!-- 9 = tundra,
!-- 10 = shrubs and interrupted woodlands (thorn shrubs for P10),
!-- 11 = wetland with plants (long grass for P10)
!-- 12 = ice cap and glacier,
!-- 13 = inland water (inland lake for P10)
!-- 14 = ocean (water for P10),
!-- 15 = urban
!
!-- SALSA variables:
    CHARACTER(LEN=20)  ::  bc_salsa_b = 'neumann'                 !< bottom boundary condition
    CHARACTER(LEN=20)  ::  bc_salsa_t = 'neumann'                 !< top boundary condition
    CHARACTER(LEN=20)  ::  depo_pcm_par = 'zhang2001'             !< or 'petroff2010'
    CHARACTER(LEN=20)  ::  depo_pcm_type = 'deciduous_broadleaf'  !< leaf type
    CHARACTER(LEN=20)  ::  depo_surf_par = 'zhang2001'            !< or 'petroff2010'
    CHARACTER(LEN=100) ::  input_file_dynamic = 'PIDS_DYNAMIC'    !< file name for dynamic input
    CHARACTER(LEN=100) ::  input_file_salsa   = 'PIDS_SALSA'      !< file name for emission data
    CHARACTER(LEN=20)  ::  salsa_emission_mode = 'no_emission'    !< 'no_emission', 'uniform',
                                                                  !< 'parameterized', 'read_from_file'

    CHARACTER(LEN=20), DIMENSION(4) ::  decycle_method_salsa =                                     &
                                                 (/'dirichlet','dirichlet','dirichlet','dirichlet'/)
                                     !< Decycling method at horizontal boundaries
                                     !< 1=left, 2=right, 3=south, 4=north
                                     !< dirichlet = initial profiles for the ghost and first 3 layers
                                     !< neumann = zero gradient

    CHARACTER(LEN=3), DIMENSION(maxspec) ::  listspec = &  !< Active aerosols
                                   (/'SO4','   ','   ','   ','   ','   ','   '/)

    INTEGER(iwp) ::  depo_pcm_par_num = 1   !< parametrisation type: 1=zhang2001, 2=petroff2010
    INTEGER(iwp) ::  depo_pcm_type_num = 0  !< index for the dry deposition type on the plant canopy
    INTEGER(iwp) ::  depo_surf_par_num = 1  !< parametrisation type: 1=zhang2001, 2=petroff2010
    INTEGER(iwp) ::  end_subrange_1a = 1    !< last index for bin subrange 1a
    INTEGER(iwp) ::  end_subrange_2a = 1    !< last index for bin subrange 2a
    INTEGER(iwp) ::  end_subrange_2b = 1    !< last index for bin subrange 2b
    INTEGER(iwp) ::  ibc_salsa_b            !< index for the bottom boundary condition
    INTEGER(iwp) ::  ibc_salsa_t            !< index for the top boundary condition
    INTEGER(iwp) ::  index_bc  = -1         !< index for black carbon (BC)
    INTEGER(iwp) ::  index_du  = -1         !< index for dust
    INTEGER(iwp) ::  index_nh  = -1         !< index for NH3
    INTEGER(iwp) ::  index_no  = -1         !< index for HNO3
    INTEGER(iwp) ::  index_oc  = -1         !< index for organic carbon (OC)
    INTEGER(iwp) ::  index_so4 = -1         !< index for SO4 or H2SO4
    INTEGER(iwp) ::  index_ss  = -1         !< index for sea salt
    INTEGER(iwp) ::  init_aerosol_type = 0  !< Initial size distribution type
                                            !< 0 = uniform (read from PARIN)
                                            !< 1 = read vertical profiles from an input file
    INTEGER(iwp) ::  init_gases_type = 0    !< Initial gas concentration type
                                            !< 0 = uniform (read from PARIN)
                                            !< 1 = read vertical profiles from an input file
    INTEGER(iwp) ::  lod_gas_emissions = 0  !< level of detail of the gaseous emission data
    INTEGER(iwp) ::  main_street_id = 0     !< lower bound of main street IDs for parameterized emission mode
    INTEGER(iwp) ::  max_street_id = 0      !< upper bound of main street IDs for parameterized emission mode
    INTEGER(iwp) ::  nbins_aerosol = 1      !< total number of size bins
    INTEGER(iwp) ::  ncc   = 1              !< number of chemical components used
    INTEGER(iwp) ::  ncomponents_mass = 1   !< total number of chemical compounds (ncc+1)
                                            !< if particle water is advected)
    INTEGER(iwp) ::  nj3 = 1                !< J3 parametrization (nucleation)
                                            !< 1 = condensational sink (Kerminen&Kulmala, 2002)
                                            !< 2 = coagulational sink (Lehtinen et al. 2007)
                                            !< 3 = coagS+self-coagulation (Anttila et al. 2010)
    INTEGER(iwp) ::  nsnucl = 0             !< Choice of the nucleation scheme:
                                            !< 0 = off
                                            !< 1 = binary nucleation
                                            !< 2 = activation type nucleation
                                            !< 3 = kinetic nucleation
                                            !< 4 = ternary nucleation
                                            !< 5 = nucleation with ORGANICs
                                            !< 6 = activation type of nucleation with H2SO4+ORG
                                            !< 7 = heteromolecular nucleation with H2SO4*ORG
                                            !< 8 = homomolecular nucleation of H2SO4
                                            !<     + heteromolecular nucleation with H2SO4*ORG
                                            !< 9 = homomolecular nucleation of H2SO4 and ORG
                                            !<     + heteromolecular nucleation with H2SO4*ORG
    INTEGER(iwp) ::  salsa_pr_count = 0     !< counter for salsa variable profiles
    INTEGER(iwp) ::  season_z01 = 1         !< For dry deposition by Zhang et al.: 1 = summer,
                                            !< 2 = autumn (no harvest yet), 3 = late autumn
                                            !< (already frost), 4 = winter, 5 = transitional spring
    INTEGER(iwp) ::  side_street_id = 0     !< lower bound of side street IDs for parameterized emission mode
    INTEGER(iwp) ::  start_subrange_1a = 1  !< start index for bin subranges: subrange 1a
    INTEGER(iwp) ::  start_subrange_2a = 1  !<                                subrange 2a
    INTEGER(iwp) ::  start_subrange_2b = 1  !<                                subrange 2b

    INTEGER(iwp), DIMENSION(nreg) ::  nbin = (/ 3, 7/)  !< Number of size bins per subrange: 1 & 2

    INTEGER(iwp), DIMENSION(ngases_salsa) ::  gas_index_chem = (/ 1, 1, 1, 1, 1/)  !< gas indices in chemistry_model_mod
                                                                                   !< 1 = H2SO4, 2 = HNO3,
                                                                                   !< 3 = NH3,   4 = OCNV, 5 = OCSV
    INTEGER(iwp), DIMENSION(ngases_salsa) ::  emission_index_chem  !< gas indices in the gas emission file
    INTEGER(iwp), DIMENSION(99) ::  salsa_pr_index  = 0            !< index for salsa profiles

    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  k_topo_top  !< vertical index of the topography top

    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE  ::  salsa_advc_flags_s !< flags used to degrade order of advection
                                                                        !< scheme for salsa variables near walls and
                                                                        !< lateral boundaries
!
!-- SALSA switches:
    LOGICAL ::  advect_particle_water   = .TRUE.   !< Advect water concentration of particles
    LOGICAL ::  decycle_salsa_lr        = .FALSE.  !< Undo cyclic boundaries: left and right
    LOGICAL ::  decycle_salsa_ns        = .FALSE.  !< Undo cyclic boundaries: north and south
    LOGICAL ::  include_emission        = .FALSE.  !< Include or not emissions
    LOGICAL ::  feedback_to_palm        = .FALSE.  !< Allow feedback due to condensation of H2O
    LOGICAL ::  nesting_salsa           = .TRUE.   !< Apply nesting for salsa
    LOGICAL ::  nesting_offline_salsa   = .TRUE.   !< Apply offline nesting for salsa
    LOGICAL ::  no_insoluble            = .FALSE.  !< Exclude insoluble chemical components
    LOGICAL ::  read_restart_data_salsa = .TRUE.   !< Read restart data for salsa
    LOGICAL ::  salsa_gases_from_chem   = .FALSE.  !< Transfer the gaseous components to SALSA
    LOGICAL ::  van_der_waals_coagc     = .FALSE.  !< Include van der Waals and viscous forces in coagulation
    LOGICAL ::  write_binary_salsa      = .TRUE.   !< read binary for salsa
!
!-- Process switches: nl* is read from the NAMELIST and is NOT changed.
!--                   ls* is the switch used and will get the value of nl*
!--                       except for special circumstances (spinup period etc.)
    LOGICAL ::  nlcoag       = .FALSE.  !< Coagulation master switch
    LOGICAL ::  lscoag       = .FALSE.  !<
    LOGICAL ::  nlcnd        = .FALSE.  !< Condensation master switch
    LOGICAL ::  lscnd        = .FALSE.  !<
    LOGICAL ::  nlcndgas     = .FALSE.  !< Condensation of precursor gases
    LOGICAL ::  lscndgas     = .FALSE.  !<
    LOGICAL ::  nlcndh2oae   = .FALSE.  !< Condensation of H2O on aerosol
    LOGICAL ::  lscndh2oae   = .FALSE.  !< particles (FALSE -> equilibrium calc.)
    LOGICAL ::  nldepo       = .FALSE.  !< Deposition master switch
    LOGICAL ::  lsdepo       = .FALSE.  !<
    LOGICAL ::  nldepo_surf  = .FALSE.  !< Deposition on vegetation master switch
    LOGICAL ::  lsdepo_surf  = .FALSE.  !<
    LOGICAL ::  nldepo_pcm   = .FALSE.  !< Deposition on walls master switch
    LOGICAL ::  lsdepo_pcm   = .FALSE.  !<
    LOGICAL ::  nldistupdate = .TRUE.   !< Size distribution update master switch
    LOGICAL ::  lsdistupdate = .FALSE.  !<
    LOGICAL ::  lspartition  = .FALSE.  !< Partition of HNO3 and NH3

    REAL(wp) ::  act_coeff = 1.0E-7_wp               !< Activation coefficient (1/s)
    REAL(wp) ::  dt_salsa  = 0.00001_wp              !< Time step of SALSA
    REAL(wp) ::  emiss_factor_main = 0.0_wp          !< relative emission factor for main streets
    REAL(wp) ::  emiss_factor_side = 0.0_wp          !< relative emission factor for side streets
    REAL(wp) ::  h2so4_init = nclim                  !< Init value for sulphuric acid gas
    REAL(wp) ::  hno3_init  = nclim                  !< Init value for nitric acid gas
    REAL(wp) ::  last_salsa_time = 0.0_wp            !< previous salsa call
    REAL(wp) ::  next_aero_emission_update = 0.0_wp  !< previous emission update
    REAL(wp) ::  next_gas_emission_update = 0.0_wp   !< previous emission update
    REAL(wp) ::  nf2a = 1.0_wp                       !< Number fraction allocated to 2a-bins
    REAL(wp) ::  nh3_init  = nclim                   !< Init value for ammonia gas
    REAL(wp) ::  ocnv_init = nclim                   !< Init value for non-volatile organic gases
    REAL(wp) ::  ocsv_init = nclim                   !< Init value for semi-volatile organic gases
    REAL(wp) ::  rhlim = 1.20_wp                     !< RH limit in %/100. Prevents unrealistical RH
    REAL(wp) ::  skip_time_do_salsa = 0.0_wp         !< Starting time of SALSA (s)
!
!-- Initial log-normal size distribution: mode diameter (dpg, metres),
!-- standard deviation (sigmag) and concentration (n_lognorm, #/m3)
    REAL(wp), DIMENSION(nmod) ::  dpg   = &
                     (/1.3E-8_wp, 5.4E-8_wp, 8.6E-7_wp, 2.0E-7_wp, 2.0E-7_wp, 2.0E-7_wp, 2.0E-7_wp/)
    REAL(wp), DIMENSION(nmod) ::  sigmag  = &
                                        (/1.8_wp, 2.16_wp, 2.21_wp, 2.0_wp, 2.0_wp, 2.0_wp, 2.0_wp/)
    REAL(wp), DIMENSION(nmod) ::  n_lognorm = &
                             (/1.04e+11_wp, 3.23E+10_wp, 5.4E+6_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp/)
!
!-- Initial mass fractions / chemical composition of the size distribution
    REAL(wp), DIMENSION(maxspec) ::  mass_fracs_a = &  !< mass fractions between
             (/1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)     !< aerosol species for A bins
    REAL(wp), DIMENSION(maxspec) ::  mass_fracs_b = &  !< mass fractions between
             (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)     !< aerosol species for B bins
    REAL(wp), DIMENSION(nreg+1) ::  reglim = &         !< Min&max diameters of size subranges
                                 (/ 3.0E-9_wp, 5.0E-8_wp, 1.0E-5_wp/)
!
!-- Initial log-normal size distribution: mode diameter (dpg, metres), standard deviation (sigmag)
!-- concentration (n_lognorm, #/m3) and mass fractions of all chemical components (listed in
!-- listspec) for both a (soluble) and b (insoluble) bins.
    REAL(wp), DIMENSION(nmod) ::  aerosol_flux_dpg   = &
                     (/1.3E-8_wp, 5.4E-8_wp, 8.6E-7_wp, 2.0E-7_wp, 2.0E-7_wp, 2.0E-7_wp, 2.0E-7_wp/)
    REAL(wp), DIMENSION(nmod) ::  aerosol_flux_sigmag  = &
                                        (/1.8_wp, 2.16_wp, 2.21_wp, 2.0_wp, 2.0_wp, 2.0_wp, 2.0_wp/)
    REAL(wp), DIMENSION(maxspec) ::  aerosol_flux_mass_fracs_a = &
                                                               (/1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    REAL(wp), DIMENSION(maxspec) ::  aerosol_flux_mass_fracs_b = &
                                                               (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    REAL(wp), DIMENSION(nmod) ::  surface_aerosol_flux = &
                                 (/1.0E+8_wp, 1.0E+9_wp, 1.0E+5_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp/)

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  bin_low_limits     !< to deliver information about
                                                               !< the lower diameters per bin
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  bc_am_t_val        !< vertical gradient of: aerosol mass
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  bc_an_t_val        !< of: aerosol number
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  bc_gt_t_val        !< salsa gases near domain top
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  gas_emission_time  !< Time array in gas emission data (s)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  nsect              !< Background number concentrations
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  massacc            !< Mass accomodation coefficients
!
!-- SALSA derived datatypes:
!
!-- Component index
    TYPE component_index
       CHARACTER(len=3), ALLOCATABLE ::  comp(:)  !< Component name
       INTEGER(iwp) ::  ncomp  !< Number of components
       INTEGER(iwp), ALLOCATABLE ::  ind(:)  !< Component index
    END TYPE component_index
!
!-- For matching LSM and USM surface types and the deposition module surface types
    TYPE match_surface
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  match_lupg  !< index for pavement / green roofs
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  match_luvw  !< index for vegetation / walls
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  match_luww  !< index for water / windows
    END TYPE match_surface
!
!-- Aerosol emission data attributes
    TYPE salsa_emission_attribute_type

       CHARACTER(LEN=25) ::   units

       CHARACTER(LEN=25), DIMENSION(:), ALLOCATABLE ::   cat_name    !<
       CHARACTER(LEN=25), DIMENSION(:), ALLOCATABLE ::   cc_name     !<
       CHARACTER(LEN=25), DIMENSION(:), ALLOCATABLE ::   unit_time   !<
       CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE ::  var_names   !<

       INTEGER(iwp) ::  lod = 0            !< level of detail
       INTEGER(iwp) ::  nbins = 10         !< number of aerosol size bins
       INTEGER(iwp) ::  ncat  = 0          !< number of emission categories
       INTEGER(iwp) ::  ncc   = 7          !< number of aerosol chemical components
       INTEGER(iwp) ::  nhoursyear = 0     !< number of hours: HOURLY mode
       INTEGER(iwp) ::  nmonthdayhour = 0  !< number of month days and hours: MDH mode
       INTEGER(iwp) ::  num_vars           !< number of variables
       INTEGER(iwp) ::  nt  = 0            !< number of time steps
       INTEGER(iwp) ::  nz  = 0            !< number of vertical levels
       INTEGER(iwp) ::  tind               !< time index for reference time in salsa emission data

       INTEGER(iwp), DIMENSION(maxspec) ::  cc_in2mod = 0   !<

       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  cat_index  !< Index of emission categories
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  cc_index   !< Index of chemical components

       REAL(wp) ::  conversion_factor  !< unit conversion factor for aerosol emissions

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  dmid         !< mean diameters of size bins (m)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rho          !< average density (kg/m3)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  time         !< time (s)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  time_factor  !< emission time factor
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  z            !< height (m)

       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  etf  !< emission time factor
       REAL(wp), DIMENSION(:,:), ALLOCATABLE :: stack_height

    END TYPE salsa_emission_attribute_type
!
!-- The default size distribution and mass composition per emission category:
!-- 1 = traffic, 2 = road dust, 3 = wood combustion, 4 = other
!-- Mass fractions: H2SO4, OC, BC, DU, SS, HNO3, NH3
    TYPE salsa_emission_mode_type

       INTEGER(iwp) ::  ndm = 3  !< number of default modes
       INTEGER(iwp) ::  ndc = 4  !< number of default categories

       CHARACTER(LEN=25), DIMENSION(1:4) ::  cat_name_table = (/'traffic exhaust', &
                                                                'road dust      ', &
                                                                'wood combustion', &
                                                                'other          '/)

       INTEGER(iwp), DIMENSION(1:4) ::  cat_input_to_model   !<

       REAL(wp), DIMENSION(1:3) ::  dpg_table = (/ 13.5E-9_wp, 1.4E-6_wp, 5.4E-8_wp/)  !<
       REAL(wp), DIMENSION(1:3) ::  ntot_table  !<
       REAL(wp), DIMENSION(1:3) ::  sigmag_table = (/ 1.6_wp, 1.4_wp, 1.7_wp /)  !<

       REAL(wp), DIMENSION(1:maxspec,1:4) ::  mass_frac_table = &  !<
          RESHAPE( (/ 0.04_wp, 0.48_wp, 0.48_wp, 0.0_wp,  0.0_wp, 0.0_wp, 0.0_wp, &
                      0.0_wp,  0.05_wp, 0.0_wp,  0.95_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                      0.0_wp,  0.5_wp,  0.5_wp,  0.0_wp,  0.0_wp, 0.0_wp, 0.0_wp, &
                      0.0_wp,  0.5_wp,  0.5_wp,  0.0_wp,  0.0_wp, 0.0_wp, 0.0_wp  &
                   /), (/maxspec,4/) )

       REAL(wp), DIMENSION(1:3,1:4) ::  pm_frac_table = & !< rel. mass
                                     RESHAPE( (/ 0.016_wp, 0.000_wp, 0.984_wp, &
                                                 0.000_wp, 1.000_wp, 0.000_wp, &
                                                 0.000_wp, 0.000_wp, 1.000_wp, &
                                                 1.000_wp, 0.000_wp, 1.000_wp  &
                                              /), (/3,4/) )

    END TYPE salsa_emission_mode_type
!
!-- Aerosol emission data values
    TYPE salsa_emission_value_type

       REAL(wp) ::  fill  !< fill value

       REAL(wp), DIMENSION(:,:), ALLOCATABLE :: mass_fracs  !< mass fractions per emis. category
       REAL(wp), DIMENSION(:,:), ALLOCATABLE :: num_fracs   !< number fractions per emis. category

       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: def_data      !< surface emission in PM
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: preproc_data  !< surface emission per category

    END TYPE salsa_emission_value_type
!
!-- Offline nesting data type
    TYPE salsa_nest_offl_type

       CHARACTER(LEN=16) ::  char_l = 'ls_forcing_left_'  !< leading substring at left boundary
       CHARACTER(LEN=17) ::  char_n = 'ls_forcing_north_' !< leading substring at north boundary
       CHARACTER(LEN=17) ::  char_r = 'ls_forcing_right_' !< leading substring at right boundary
       CHARACTER(LEN=17) ::  char_s = 'ls_forcing_south_' !< leading substring at south boundary
       CHARACTER(LEN=15) ::  char_t = 'ls_forcing_top_'   !< leading substring at top boundary

       CHARACTER(LEN=5), DIMENSION(1:ngases_salsa) ::  gas_name = (/'H2SO4','HNO3 ','NH3  ','OCNV ','OCSV '/)

       CHARACTER(LEN=25),  DIMENSION(:), ALLOCATABLE ::  cc_name    !< chemical component name
       CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE ::  var_names  !< list of variable names

       INTEGER(iwp) ::  id_dynamic  !< NetCDF id of dynamic input file
       INTEGER(iwp) ::  ncc         !< number of aerosol chemical components
       INTEGER(iwp) ::  nt          !< number of time levels in dynamic input file
       INTEGER(iwp) ::  nzu         !< number of vertical levels on scalar grid in dynamic input file
       INTEGER(iwp) ::  tind        !< time index for reference time in mesoscale-offline nesting
       INTEGER(iwp) ::  tind_p      !< time index for following time in mesoscale-offline nesting

       INTEGER(iwp), DIMENSION(maxspec) ::  cc_in2mod = 0  !< to transfer chemical composition from input to model

       LOGICAL ::  init  = .FALSE. !< flag indicating the initialisation of offline nesting

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  dmid      !< vertical profile of aerosol bin diameters
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  time      !< time in dynamic input file
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  zu_atmos  !< zu in dynamic input file

       REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  gconc_left   !< gas conc. at left boundary
       REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  gconc_north  !< gas conc. at north boundary
       REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  gconc_right  !< gas conc. at right boundary
       REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  gconc_south  !< gas conc. at south boundary
       REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  gconc_top    !< gas conc.at top boundary
       REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  mconc_left   !< aerosol mass conc. at left boundary
       REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  mconc_north  !< aerosol mass conc. at north boundary
       REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  mconc_right  !< aerosol mass conc. at right boundary
       REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  mconc_south  !< aerosol mass conc. at south boundary
       REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  mconc_top    !< aerosol mass conc. at top boundary
       REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  nconc_left   !< aerosol number conc. at left boundary
       REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  nconc_north  !< aerosol number conc. at north boundary
       REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  nconc_right  !< aerosol number conc. at right boundary
       REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  nconc_south  !< aerosol number conc. at south boundary
       REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  nconc_top    !< aerosol number conc. at top boundary

    END TYPE salsa_nest_offl_type
!
!-- Prognostic variable: Aerosol size bin information (number (#/m3) and mass (kg/m3) concentration)
!-- and the concentration of gaseous tracers (#/m3). Gas tracers are contained sequentially in
!-- dimension 4 as:
!-- 1. H2SO4, 2. HNO3, 3. NH3, 4. OCNV (non-volatile organics), 5. OCSV (semi-volatile)
    TYPE salsa_variable

       REAL(wp), DIMENSION(:), ALLOCATABLE     ::  init  !<

       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  diss_s     !<
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  flux_s     !<
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  source     !<
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  sums_ws_l  !<

       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  diss_l  !<
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  flux_l  !<

       REAL(wp), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::  conc     !<
       REAL(wp), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::  conc_p   !<
       REAL(wp), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::  tconc_m  !<

    END TYPE salsa_variable
!
!-- Datatype used to store information about the binned size distributions of aerosols
    TYPE t_section

       REAL(wp) ::  dmid     !< bin middle diameter (m)
       REAL(wp) ::  vhilim   !< bin volume at the high limit
       REAL(wp) ::  vlolim   !< bin volume at the low limit
       REAL(wp) ::  vratiohi !< volume ratio between the center and high limit
       REAL(wp) ::  vratiolo !< volume ratio between the center and low limit
       !******************************************************
       ! ^ Do NOT change the stuff above after initialization !
       !******************************************************
       REAL(wp) ::  core    !< Volume of dry particle
       REAL(wp) ::  dwet    !< Wet diameter or mean droplet diameter (m)
       REAL(wp) ::  numc    !< Number concentration of particles/droplets (#/m3)
       REAL(wp) ::  veqh2o  !< Equilibrium H2O concentration for each particle

       REAL(wp), DIMENSION(maxspec+1) ::  volc !< Volume concentrations (m^3/m^3) of aerosols +
                                               !< water. Since most of the stuff in SALSA is hard
                                               !< coded, these *have to be* in the order
                                               !< 1:SO4, 2:OC, 3:BC, 4:DU, 5:SS, 6:NO, 7:NH, 8:H2O
    END TYPE t_section

    TYPE(salsa_emission_attribute_type) ::  aero_emission_att  !< emission attributes
    TYPE(salsa_emission_value_type)     ::  aero_emission      !< emission values
    TYPE(salsa_emission_mode_type)      ::  def_modes          !< default emission modes

    TYPE(chem_emis_att_type) ::  chem_emission_att  !< chemistry emission attributes

    TYPE(chem_emis_val_type), DIMENSION(:), ALLOCATABLE ::  chem_emission  !< chemistry emissions

    TYPE(t_section), DIMENSION(:), ALLOCATABLE ::  aero  !< local aerosol properties

    TYPE(match_surface) ::  lsm_to_depo_h  !< to match the deposition module and horizontal LSM surfaces
    TYPE(match_surface) ::  usm_to_depo_h  !< to match the deposition module and horizontal USM surfaces

    TYPE(match_surface), DIMENSION(0:3) ::  lsm_to_depo_v  !< to match the deposition mod. and vertical LSM surfaces
    TYPE(match_surface), DIMENSION(0:3) ::  usm_to_depo_v  !< to match the deposition mod. and vertical USM surfaces
!
!-- SALSA variables: as x = x(k,j,i,bin).
!-- The 4th dimension contains all the size bins sequentially for each aerosol species  + water.
!
!-- Prognostic variables:
!
!-- Number concentration (#/m3)
    TYPE(salsa_variable), DIMENSION(:), ALLOCATABLE, TARGET ::  aerosol_number  !<
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET ::  nconc_1  !<
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET ::  nconc_2  !<
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET ::  nconc_3  !<
!
!-- Mass concentration (kg/m3)
    TYPE(salsa_variable), DIMENSION(:), ALLOCATABLE, TARGET ::  aerosol_mass  !<
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET ::  mconc_1  !<
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET ::  mconc_2  !<
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET ::  mconc_3  !<
!
!-- Gaseous concentrations (#/m3)
    TYPE(salsa_variable), DIMENSION(:), ALLOCATABLE, TARGET ::  salsa_gas  !<
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET ::  gconc_1  !<
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET ::  gconc_2  !<
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET ::  gconc_3  !<
!
!-- Diagnostic tracers
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  sedim_vd  !< sedimentation velocity per bin (m/s)
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  ra_dry    !< aerosol dry radius (m)

!-- Particle component index tables
    TYPE(component_index) :: prtcl  !< Contains "getIndex" which gives the index for a given aerosol
                                    !< component name: 1:SO4, 2:OC, 3:BC, 4:DU, 5:SS, 6:NO, 7:NH, 8:H2O
!
!-- Offline nesting:
    TYPE(salsa_nest_offl_type) ::  salsa_nest_offl  !< data structure for offline nesting
!
!-- Data output arrays:
!
!-- Integrated:
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  ldsa_av  !< lung-deposited surface area
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  ntot_av  !< total number concentration
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  nufp_av  !< ultrafine particles (UFP)
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  pm01_av  !< PM0.1
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  pm25_av  !< PM2.5
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  pm10_av  !< PM10
!
!-- Bin specific mass and number concentrations:
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET ::  mbins_av  !< bin mas
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET ::  nbins_av  !< bin number
!
!-- Gases:
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET ::  salsa_gases_av  !< gases
!
!-- In the particle phase:
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  s_h2o_av  !< liquid water
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET ::  s_mass_av  !< mass components

!
!-- PALM interfaces:

    INTERFACE salsa_actions
       MODULE PROCEDURE salsa_actions
       MODULE PROCEDURE salsa_actions_ij
    END INTERFACE salsa_actions

    INTERFACE salsa_3d_data_averaging
       MODULE PROCEDURE salsa_3d_data_averaging
    END INTERFACE salsa_3d_data_averaging

    INTERFACE salsa_boundary_conds
       MODULE PROCEDURE salsa_boundary_conds
       MODULE PROCEDURE salsa_boundary_conds_decycle
    END INTERFACE salsa_boundary_conds

    INTERFACE salsa_boundary_conditions
       MODULE PROCEDURE salsa_boundary_conditions
    END INTERFACE salsa_boundary_conditions

    INTERFACE salsa_check_data_output
       MODULE PROCEDURE salsa_check_data_output
    END INTERFACE salsa_check_data_output

    INTERFACE salsa_check_data_output_pr
       MODULE PROCEDURE salsa_check_data_output_pr
    END INTERFACE salsa_check_data_output_pr

    INTERFACE salsa_check_parameters
       MODULE PROCEDURE salsa_check_parameters
    END INTERFACE salsa_check_parameters

    INTERFACE salsa_data_output_2d
       MODULE PROCEDURE salsa_data_output_2d
    END INTERFACE salsa_data_output_2d

    INTERFACE salsa_data_output_3d
       MODULE PROCEDURE salsa_data_output_3d
    END INTERFACE salsa_data_output_3d

    INTERFACE salsa_data_output_mask
       MODULE PROCEDURE salsa_data_output_mask
    END INTERFACE salsa_data_output_mask

    INTERFACE salsa_define_netcdf_grid
       MODULE PROCEDURE salsa_define_netcdf_grid
    END INTERFACE salsa_define_netcdf_grid

    INTERFACE salsa_emission_update
       MODULE PROCEDURE salsa_emission_update
    END INTERFACE salsa_emission_update

    INTERFACE salsa_exchange_horiz_bounds
       MODULE PROCEDURE salsa_exchange_horiz_bounds
    END INTERFACE salsa_exchange_horiz_bounds

    INTERFACE salsa_header
       MODULE PROCEDURE salsa_header
    END INTERFACE salsa_header

    INTERFACE salsa_init
       MODULE PROCEDURE salsa_init
    END INTERFACE salsa_init

    INTERFACE salsa_init_arrays
       MODULE PROCEDURE salsa_init_arrays
    END INTERFACE salsa_init_arrays

    INTERFACE salsa_nesting_offl_bc
       MODULE PROCEDURE salsa_nesting_offl_bc
    END INTERFACE salsa_nesting_offl_bc

    INTERFACE salsa_nesting_offl_init
       MODULE PROCEDURE salsa_nesting_offl_init
    END INTERFACE salsa_nesting_offl_init

    INTERFACE salsa_nesting_offl_input
       MODULE PROCEDURE salsa_nesting_offl_input
    END INTERFACE salsa_nesting_offl_input

    INTERFACE salsa_non_advective_processes
       MODULE PROCEDURE salsa_non_advective_processes
       MODULE PROCEDURE salsa_non_advective_processes_ij
    END INTERFACE salsa_non_advective_processes

    INTERFACE salsa_parin
       MODULE PROCEDURE salsa_parin
    END INTERFACE salsa_parin

    INTERFACE salsa_prognostic_equations
       MODULE PROCEDURE salsa_prognostic_equations
       MODULE PROCEDURE salsa_prognostic_equations_ij
    END INTERFACE salsa_prognostic_equations

    INTERFACE salsa_rrd_local
       MODULE PROCEDURE salsa_rrd_local
    END INTERFACE salsa_rrd_local

    INTERFACE salsa_statistics
       MODULE PROCEDURE salsa_statistics
    END INTERFACE salsa_statistics

    INTERFACE salsa_swap_timelevel
       MODULE PROCEDURE salsa_swap_timelevel
    END INTERFACE salsa_swap_timelevel

    INTERFACE salsa_tendency
       MODULE PROCEDURE salsa_tendency
       MODULE PROCEDURE salsa_tendency_ij
    END INTERFACE salsa_tendency

    INTERFACE salsa_wrd_local
       MODULE PROCEDURE salsa_wrd_local
    END INTERFACE salsa_wrd_local


    SAVE

    PRIVATE
!
!-- Public functions:
    PUBLIC salsa_3d_data_averaging,       &
           salsa_actions,                 &
           salsa_boundary_conds,          &
           salsa_boundary_conditions,     &
           salsa_check_data_output,       &
           salsa_check_data_output_pr,    &
           salsa_check_parameters,        &
           salsa_data_output_2d,          &
           salsa_data_output_3d,          &
           salsa_data_output_mask,        &
           salsa_define_netcdf_grid,      &
           salsa_diagnostics,             &
           salsa_emission_update,         &
           salsa_exchange_horiz_bounds,   &
           salsa_header,                  &
           salsa_init,                    &
           salsa_init_arrays,             &
           salsa_nesting_offl_bc,         &
           salsa_nesting_offl_init,       &
           salsa_nesting_offl_input,      &
           salsa_non_advective_processes, &
           salsa_parin,                   &
           salsa_prognostic_equations,    &
           salsa_rrd_local,               &
           salsa_statistics,              &
           salsa_swap_timelevel,          &
           salsa_wrd_local

!
!-- Public parameters, constants and initial values
    PUBLIC bc_am_t_val,           &
           bc_an_t_val,           &
           bc_gt_t_val,           &
           ibc_salsa_b,           &
           init_aerosol_type,     &
           init_gases_type,       &
           nesting_salsa,         &
           nesting_offline_salsa, &
           salsa_gases_from_chem, &
           skip_time_do_salsa
!
!-- Public variables
    PUBLIC aerosol_mass,     &
           aerosol_number,   &
           gconc_2,          &
           mconc_2,          &
           nbins_aerosol,    &
           ncomponents_mass, &
           nconc_2,          &
           ngases_salsa,     &
           salsa_gas,        &
           salsa_nest_offl


 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Parin for &salsa_par for new modules
!------------------------------------------------------------------------------!
 SUBROUTINE salsa_parin

    USE control_parameters,                                                                        &
        ONLY:  data_output_pr

    IMPLICIT NONE

    CHARACTER(LEN=80) ::  line   !< dummy string that contains the current line of parameter file

    INTEGER(iwp) ::  i                 !< loop index
    INTEGER(iwp) ::  max_pr_salsa_tmp  !< dummy variable

    NAMELIST /salsa_parameters/      aerosol_flux_dpg,                         &
                                     aerosol_flux_mass_fracs_a,                &
                                     aerosol_flux_mass_fracs_b,                &
                                     aerosol_flux_sigmag,                      &
                                     advect_particle_water,                    &
                                     bc_salsa_b,                               &
                                     bc_salsa_t,                               &
                                     decycle_salsa_lr,                         &
                                     decycle_method_salsa,                     &
                                     decycle_salsa_ns,                         &
                                     depo_pcm_par,                             &
                                     depo_pcm_type,                            &
                                     depo_surf_par,                            &
                                     dpg,                                      &
                                     dt_salsa,                                 &
                                     emiss_factor_main,                        &
                                     emiss_factor_side,                        &
                                     feedback_to_palm,                         &
                                     h2so4_init,                               &
                                     hno3_init,                                &
                                     listspec,                                 &
                                     main_street_id,                           &
                                     mass_fracs_a,                             &
                                     mass_fracs_b,                             &
                                     max_street_id,                            &
                                     n_lognorm,                                &
                                     nbin,                                     &
                                     nesting_salsa,                            &
                                     nesting_offline_salsa,                    &
                                     nf2a,                                     &
                                     nh3_init,                                 &
                                     nj3,                                      &
                                     nlcnd,                                    &
                                     nlcndgas,                                 &
                                     nlcndh2oae,                               &
                                     nlcoag,                                   &
                                     nldepo,                                   &
                                     nldepo_pcm,                               &
                                     nldepo_surf,                              &
                                     nldistupdate,                             &
                                     nsnucl,                                   &
                                     ocnv_init,                                &
                                     ocsv_init,                                &
                                     read_restart_data_salsa,                  &
                                     reglim,                                   &
                                     salsa,                                    &
                                     salsa_emission_mode,                      &
                                     season_z01,                               &
                                     sigmag,                                   &
                                     side_street_id,                           &
                                     skip_time_do_salsa,                       &
                                     surface_aerosol_flux,                     &
                                     van_der_waals_coagc,                      &
                                     write_binary_salsa

    line = ' '
!
!-- Try to find salsa package
    REWIND ( 11 )
    line = ' '
    DO WHILE ( INDEX( line, '&salsa_parameters' ) == 0 )
       READ ( 11, '(A)', END=10 )  line
    ENDDO
    BACKSPACE ( 11 )
!
!-- Read user-defined namelist
    READ ( 11, salsa_parameters )
!
!-- Enable salsa (salsa switch in modules.f90)
    salsa = .TRUE.

 10 CONTINUE
!
!-- Update the number of output profiles
    max_pr_salsa_tmp = 0
    i = 1
    DO WHILE ( data_output_pr(i) /= ' '  .AND.  i <= 100 )
       IF ( TRIM( data_output_pr(i)(1:6) ) == 'salsa_' )  max_pr_salsa_tmp = max_pr_salsa_tmp + 1
       i = i + 1
    ENDDO
    IF ( max_pr_salsa_tmp > 0 )  max_pr_salsa = max_pr_salsa_tmp

 END SUBROUTINE salsa_parin

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check parameters routine for salsa.
!------------------------------------------------------------------------------!
 SUBROUTINE salsa_check_parameters

    USE control_parameters,                                                                        &
        ONLY:  child_domain, humidity, initializing_actions, nesting_offline

    IMPLICIT NONE

!
!-- Check that humidity is switched on
    IF ( salsa  .AND.  .NOT.  humidity )  THEN
       WRITE( message_string, * ) 'salsa = ', salsa, ' is not allowed with humidity = ', humidity
       CALL message( 'salsa_check_parameters', 'PA0594', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- For nested runs, explicitly set nesting boundary conditions.
    IF ( child_domain )  THEN
       IF ( nesting_salsa )  THEN
          bc_salsa_t = 'nested'
       ELSE
          bc_salsa_t = 'neumann'
       ENDIF
    ENDIF
!
!-- Set boundary conditions also in case the model is offline-nested in larger-scale models.
    IF ( nesting_offline )  THEN
       IF ( nesting_offline_salsa )  THEN
          bc_salsa_t = 'nesting_offline'
       ELSE
          bc_salsa_t = 'neumann'
       ENDIF
    ENDIF
!
!-- Set bottom boundary condition flag
    IF ( bc_salsa_b == 'dirichlet' )  THEN
       ibc_salsa_b = 0
    ELSEIF ( bc_salsa_b == 'neumann' )  THEN
       ibc_salsa_b = 1
    ELSE
       message_string = 'unknown boundary condition: bc_salsa_b = "' // TRIM( bc_salsa_t ) // '"'
       CALL message( 'salsa_check_parameters', 'PA0595', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Set top boundary conditions flag
    IF ( bc_salsa_t == 'dirichlet' )  THEN
       ibc_salsa_t = 0
    ELSEIF ( bc_salsa_t == 'neumann' )  THEN
       ibc_salsa_t = 1
    ELSEIF ( bc_salsa_t == 'initial_gradient' )  THEN
       ibc_salsa_t = 2
    ELSEIF ( bc_salsa_t == 'nested'  .OR.  bc_salsa_t == 'nesting_offline' )  THEN
       ibc_salsa_t = 3
    ELSE
       message_string = 'unknown boundary condition: bc_salsa_t = "' // TRIM( bc_salsa_t ) // '"'
       CALL message( 'salsa_check_parameters', 'PA0596', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Check J3 parametrisation
    IF ( nj3 < 1  .OR.  nj3 > 3 )  THEN
       message_string = 'unknown nj3 (must be 1-3)'
       CALL message( 'salsa_check_parameters', 'PA0597', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Check bottom boundary condition in case of surface emissions
    IF ( salsa_emission_mode /= 'no_emission'  .AND.  ibc_salsa_b  == 0 ) THEN
       message_string = 'salsa_emission_mode /= "no_emission" requires bc_salsa_b = "Neumann"'
       CALL message( 'salsa_check_parameters','PA0598', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Check whether emissions are applied
    IF ( salsa_emission_mode /= 'no_emission' )  include_emission = .TRUE.
!
!-- Set the initialisation type: background concentration are read from PIDS_DYNAMIC if
!-- initializing_actions = 'inifor set_constant_profiles'
    IF ( INDEX( initializing_actions, 'inifor' ) /= 0 )  THEN
       init_aerosol_type = 1
       init_gases_type = 1
    ENDIF
!
!-- If the run is not a restart run, set read_restart_data to .FALSE.
    IF ( TRIM( initializing_actions ) /= 'read_restart_data' )  THEN
       read_restart_data_salsa = .FALSE.
    ENDIF

 END SUBROUTINE salsa_check_parameters

!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining appropriate grid for netcdf variables.
!> It is called out from subroutine netcdf.
!> Same grid as for other scalars (see netcdf_interface_mod.f90)
!------------------------------------------------------------------------------!
 SUBROUTINE salsa_define_netcdf_grid( var, found, grid_x, grid_y, grid_z )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(OUT) ::  grid_x   !<
    CHARACTER(LEN=*), INTENT(OUT) ::  grid_y   !<
    CHARACTER(LEN=*), INTENT(OUT) ::  grid_z   !<
    CHARACTER(LEN=*), INTENT(IN)  ::  var      !<

    LOGICAL, INTENT(OUT) ::  found   !<

    found  = .TRUE.
!
!-- Check for the grid

    IF ( var(1:6) == 'salsa_' )  THEN  ! same grid for all salsa output variables
       grid_x = 'x'
       grid_y = 'y'
       grid_z = 'zu'
    ELSE
       found  = .FALSE.
       grid_x = 'none'
       grid_y = 'none'
       grid_z = 'none'
    ENDIF

 END SUBROUTINE salsa_define_netcdf_grid

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Header output for new module
!------------------------------------------------------------------------------!
 SUBROUTINE salsa_header( io )

    USE indices,                                                                                   &
        ONLY:  nx, ny, nz

    IMPLICIT NONE
 
    INTEGER(iwp), INTENT(IN) ::  io   !< Unit of the output file
!
!-- Write SALSA header
    WRITE( io, 1 )
    WRITE( io, 2 ) skip_time_do_salsa
    WRITE( io, 3 ) dt_salsa
    WRITE( io, 4 )  nz, ny, nx, nbins_aerosol
    IF ( advect_particle_water )  THEN
       WRITE( io, 5 )  nz, ny, nx, ncomponents_mass*nbins_aerosol, advect_particle_water
    ELSE
       WRITE( io, 5 )  nz, ny, nx, ncc*nbins_aerosol, advect_particle_water
    ENDIF
    IF ( .NOT. salsa_gases_from_chem )  THEN
       WRITE( io, 6 )  nz, ny, nx, ngases_salsa, salsa_gases_from_chem
    ENDIF
    WRITE( io, 7 )
    IF ( nsnucl > 0 )   WRITE( io, 8 ) nsnucl, nj3
    IF ( nlcoag )       WRITE( io, 9 )
    IF ( nlcnd )        WRITE( io, 10 ) nlcndgas, nlcndh2oae
    IF ( lspartition )  WRITE( io, 11 )
    IF ( nldepo )       WRITE( io, 12 ) nldepo_pcm, nldepo_surf
    WRITE( io, 13 )  reglim, nbin, ( aero(:)%vlolim / api6 )**0.33333333_wp
    WRITE( io, 25 )  aero(:)%dmid
    IF ( init_aerosol_type == 0 )  WRITE( io, 14 ) nsect
    WRITE( io, 15 ) ncc, listspec, mass_fracs_a, mass_fracs_b
    IF ( .NOT. salsa_gases_from_chem )  THEN
       WRITE( io, 16 ) ngases_salsa, h2so4_init, hno3_init, nh3_init, ocnv_init, ocsv_init
    ENDIF
    WRITE( io, 17 )  init_aerosol_type, init_gases_type
    IF ( init_aerosol_type == 0 )  THEN
       WRITE( io, 18 )  dpg, sigmag, n_lognorm
    ELSE
       WRITE( io, 19 )
    ENDIF
    IF ( nesting_salsa )  WRITE( io, 20 )  nesting_salsa
    IF ( nesting_offline_salsa )  WRITE( io, 21 )  nesting_offline_salsa
    WRITE( io, 22 ) salsa_emission_mode
    IF ( salsa_emission_mode == 'uniform' )  THEN
       WRITE( io, 23 ) surface_aerosol_flux, aerosol_flux_dpg, aerosol_flux_sigmag,                &
                       aerosol_flux_mass_fracs_a
    ENDIF
    IF ( SUM( aerosol_flux_mass_fracs_b ) > 0.0_wp  .OR. salsa_emission_mode == 'read_from_file' ) &
    THEN
       WRITE( io, 24 )
    ENDIF

1   FORMAT (//' SALSA information:'/                                                               &
              ' ------------------------------'/)
2   FORMAT   ('    Starts at: skip_time_do_salsa = ', F10.2, '  s')
3   FORMAT  (/'    Timestep: dt_salsa = ', F6.2, '  s')
4   FORMAT  (/'    Array shape (z,y,x,bins):'/                                                     &
              '       aerosol_number:  ', 4(I5)) 
5   FORMAT  (/'       aerosol_mass:    ', 4(I5),/                                                  &
              '       (advect_particle_water = ', L1, ')')
6   FORMAT   ('       salsa_gas: ', 4(I5),/                                                        &
              '       (salsa_gases_from_chem = ', L1, ')')
7   FORMAT  (/'    Aerosol dynamic processes included: ')
8   FORMAT  (/'       nucleation (scheme = ', I1, ' and J3 parametrization = ', I1, ')')
9   FORMAT  (/'       coagulation')
10  FORMAT  (/'       condensation (of precursor gases = ', L1, ' and water vapour = ', L1, ')' )
11  FORMAT  (/'       dissolutional growth by HNO3 and NH3')
12  FORMAT  (/'       dry deposition (on vegetation = ', L1, ' and on topography = ', L1, ')')
13  FORMAT  (/'    Aerosol bin subrange limits (in metres): ',  3(ES10.2E3), /                     &
              '    Number of size bins for each aerosol subrange: ', 2I3,/                         &
              '    Aerosol bin lower limits (in metres): ', 12(ES10.2E3))
25  FORMAT  (/'    Bin geometric mean diameters (in metres): ', 12(ES10.2E3))
14  FORMAT   ('    Initial number concentration in bins at the lowest level (#/m**3):', 9(ES10.2E3))
15  FORMAT  (/'    Number of chemical components used: ', I1,/                                     &
              '       Species: ',7(A6),/                                                           &
              '    Initial relative contribution of each species to particle volume in:',/         &
              '       a-bins: ', 7(F6.3),/                                                         &
              '       b-bins: ', 7(F6.3))
16  FORMAT  (/'    Number of gaseous tracers used: ', I1,/                                         &
              '    Initial gas concentrations:',/                                                  &
              '       H2SO4: ',ES12.4E3, ' #/m**3',/                                               &
              '       HNO3:  ',ES12.4E3, ' #/m**3',/                                               &
              '       NH3:   ',ES12.4E3, ' #/m**3',/                                               &
              '       OCNV:  ',ES12.4E3, ' #/m**3',/                                               &
              '       OCSV:  ',ES12.4E3, ' #/m**3')
17   FORMAT (/'   Initialising concentrations: ', /                                                &
              '      Aerosol size distribution: init_aerosol_type = ', I1,/                        &
              '      Gas concentrations: init_gases_type = ', I1 )
18   FORMAT ( '      Mode diametres: dpg(nmod) = ', 7(F7.3), ' (m)', /                             &
              '      Standard deviation: sigmag(nmod) = ', 7(F7.2),/                               &
              '      Number concentration: n_lognorm(nmod) = ', 7(ES12.4E3), ' (#/m3)' )
19   FORMAT (/'      Size distribution read from a file.')
20   FORMAT (/'   Nesting for salsa variables: ', L1 )
21   FORMAT (/'   Offline nesting for salsa variables: ', L1 )
22   FORMAT (/'   Emissions: salsa_emission_mode = ', A )
23   FORMAT (/'      surface_aerosol_flux = ', ES12.4E3, ' #/m**2/s', /                            &
              '      aerosol_flux_dpg     =  ', 7(F7.3), ' (m)', /                                 &
              '      aerosol_flux_sigmag  =  ', 7(F7.2), /                                         &
              '      aerosol_mass_fracs_a =  ', 7(ES12.4E3) )
24   FORMAT (/'      (currently all emissions are soluble!)')

 END SUBROUTINE salsa_header

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate SALSA arrays and define pointers if required
!------------------------------------------------------------------------------!
 SUBROUTINE salsa_init_arrays

    USE advec_ws,                                                                                  &
        ONLY: ws_init_flags_scalar

    USE surface_mod,                                                                               &
        ONLY:  surf_def_h, surf_def_v, surf_lsm_h, surf_lsm_v, surf_usm_h, surf_usm_v

    IMPLICIT NONE

    INTEGER(iwp) ::  gases_available !< Number of available gas components in the chemistry model
    INTEGER(iwp) ::  i               !< loop index for allocating
    INTEGER(iwp) ::  ii              !< index for indexing chemical components
    INTEGER(iwp) ::  l               !< loop index for allocating: surfaces
    INTEGER(iwp) ::  lsp             !< loop index for chem species in the chemistry model

    gases_available = 0
!
!-- Allocate prognostic variables (see salsa_swap_timelevel)
!
!-- Set derived indices:
!-- (This does the same as the subroutine salsa_initialize in SALSA/UCLALES-SALSA)
    start_subrange_1a = 1  ! 1st index of subrange 1a
    start_subrange_2a = start_subrange_1a + nbin(1)  ! 1st index of subrange 2a
    end_subrange_1a   = start_subrange_2a - 1        ! last index of subrange 1a
    end_subrange_2a   = end_subrange_1a + nbin(2)    ! last index of subrange 2a

!
!-- If the fraction of insoluble aerosols in subrange 2 is zero: do not allocate arrays for them
    IF ( nf2a > 0.999999_wp  .AND.  SUM( mass_fracs_b ) < 0.00001_wp )  THEN
       no_insoluble = .TRUE.
       start_subrange_2b = end_subrange_2a+1  ! 1st index of subrange 2b
       end_subrange_2b   = end_subrange_2a    ! last index of subrange 2b
    ELSE
       start_subrange_2b = start_subrange_2a + nbin(2)  ! 1st index of subrange 2b
       end_subrange_2b   = end_subrange_2a + nbin(2)    ! last index of subrange 2b
    ENDIF

    nbins_aerosol = end_subrange_2b   ! total number of aerosol size bins
!
!-- Create index tables for different aerosol components
    CALL component_index_constructor( prtcl, ncc, maxspec, listspec )

    ncomponents_mass = ncc
    IF ( advect_particle_water )  ncomponents_mass = ncc + 1  ! Add water
!
!-- Indices for chemical components used (-1 = not used)
    ii = 0
    IF ( is_used( prtcl, 'SO4' ) )  THEN
       index_so4 = get_index( prtcl,'SO4' )
       ii = ii + 1
    ENDIF
    IF ( is_used( prtcl,'OC' ) )  THEN
       index_oc = get_index(prtcl, 'OC')
       ii = ii + 1
    ENDIF
    IF ( is_used( prtcl, 'BC' ) )  THEN
       index_bc = get_index( prtcl, 'BC' )
       ii = ii + 1
    ENDIF
    IF ( is_used( prtcl, 'DU' ) )  THEN
       index_du = get_index( prtcl, 'DU' )
       ii = ii + 1
    ENDIF
    IF ( is_used( prtcl, 'SS' ) )  THEN
       index_ss = get_index( prtcl, 'SS' )
       ii = ii + 1
    ENDIF
    IF ( is_used( prtcl, 'NO' ) )  THEN
       index_no = get_index( prtcl, 'NO' )
       ii = ii + 1
    ENDIF
    IF ( is_used( prtcl, 'NH' ) )  THEN
       index_nh = get_index( prtcl, 'NH' )
       ii = ii + 1
    ENDIF
!
!-- All species must be known
    IF ( ii /= ncc )  THEN
       message_string = 'Unknown aerosol species/component(s) given in the initialization'
       CALL message( 'salsa_mod: salsa_init', 'PA0600', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Allocate:
    ALLOCATE( aero(nbins_aerosol), bc_am_t_val(nbins_aerosol*ncomponents_mass),                    &
              bc_an_t_val(nbins_aerosol), bc_gt_t_val(ngases_salsa), bin_low_limits(nbins_aerosol),&
              nsect(nbins_aerosol), massacc(nbins_aerosol) )
    ALLOCATE( k_topo_top(nysg:nyng,nxlg:nxrg) )
    IF ( nldepo ) ALLOCATE( sedim_vd(nzb:nzt+1,nysg:nyng,nxlg:nxrg,nbins_aerosol) )
    ALLOCATE( ra_dry(nzb:nzt+1,nysg:nyng,nxlg:nxrg,nbins_aerosol) )
!
!-- Initialise the sectional particle size distribution
    CALL set_sizebins
!
!-- Aerosol number concentration
    ALLOCATE( aerosol_number(nbins_aerosol) )
    ALLOCATE( nconc_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg,nbins_aerosol),                                &
              nconc_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg,nbins_aerosol),                                &
              nconc_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg,nbins_aerosol) )
    nconc_1 = 0.0_wp
    nconc_2 = 0.0_wp
    nconc_3 = 0.0_wp

    DO i = 1, nbins_aerosol
       aerosol_number(i)%conc(nzb:nzt+1,nysg:nyng,nxlg:nxrg)    => nconc_1(:,:,:,i)
       aerosol_number(i)%conc_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg)  => nconc_2(:,:,:,i)
       aerosol_number(i)%tconc_m(nzb:nzt+1,nysg:nyng,nxlg:nxrg) => nconc_3(:,:,:,i)
       ALLOCATE( aerosol_number(i)%flux_s(nzb+1:nzt,0:threads_per_task-1),                         &
                 aerosol_number(i)%diss_s(nzb+1:nzt,0:threads_per_task-1),                         &
                 aerosol_number(i)%flux_l(nzb+1:nzt,nys:nyn,0:threads_per_task-1),                 &
                 aerosol_number(i)%diss_l(nzb+1:nzt,nys:nyn,0:threads_per_task-1),                 &
                 aerosol_number(i)%init(nzb:nzt+1),                                                &
                 aerosol_number(i)%sums_ws_l(nzb:nzt+1,0:threads_per_task-1) )
       aerosol_number(i)%init = nclim
       IF ( include_emission  .OR.  ( nldepo  .AND.  nldepo_surf ) )  THEN
          ALLOCATE( aerosol_number(i)%source(nys:nyn,nxl:nxr) )
          aerosol_number(i)%source = 0.0_wp
       ENDIF
    ENDDO

!
!-- Aerosol mass concentration
    ALLOCATE( aerosol_mass(ncomponents_mass*nbins_aerosol) )
    ALLOCATE( mconc_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg,ncomponents_mass*nbins_aerosol),               &
              mconc_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg,ncomponents_mass*nbins_aerosol),               &
              mconc_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg,ncomponents_mass*nbins_aerosol) )
    mconc_1 = 0.0_wp
    mconc_2 = 0.0_wp
    mconc_3 = 0.0_wp

    DO i = 1, ncomponents_mass*nbins_aerosol
       aerosol_mass(i)%conc(nzb:nzt+1,nysg:nyng,nxlg:nxrg)    => mconc_1(:,:,:,i)
       aerosol_mass(i)%conc_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg)  => mconc_2(:,:,:,i)
       aerosol_mass(i)%tconc_m(nzb:nzt+1,nysg:nyng,nxlg:nxrg) => mconc_3(:,:,:,i)
       ALLOCATE( aerosol_mass(i)%flux_s(nzb+1:nzt,0:threads_per_task-1),                           &
                 aerosol_mass(i)%diss_s(nzb+1:nzt,0:threads_per_task-1),                           &
                 aerosol_mass(i)%flux_l(nzb+1:nzt,nys:nyn,0:threads_per_task-1),                   &
                 aerosol_mass(i)%diss_l(nzb+1:nzt,nys:nyn,0:threads_per_task-1),                   &
                 aerosol_mass(i)%init(nzb:nzt+1),                                                  &
                 aerosol_mass(i)%sums_ws_l(nzb:nzt+1,0:threads_per_task-1)  )
       aerosol_mass(i)%init = mclim
       IF ( include_emission  .OR.  ( nldepo  .AND.  nldepo_surf ) )  THEN
          ALLOCATE( aerosol_mass(i)%source(nys:nyn,nxl:nxr) )
          aerosol_mass(i)%source = 0.0_wp
       ENDIF
    ENDDO

!
!-- Surface fluxes: answs = aerosol number, amsws = aerosol mass
!
!-- Horizontal surfaces: default type
    DO  l = 0, 2   ! upward (l=0), downward (l=1) and model top (l=2)
       ALLOCATE( surf_def_h(l)%answs( 1:surf_def_h(l)%ns, nbins_aerosol ) )
       ALLOCATE( surf_def_h(l)%amsws( 1:surf_def_h(l)%ns, nbins_aerosol*ncomponents_mass ) )
       surf_def_h(l)%answs = 0.0_wp
       surf_def_h(l)%amsws = 0.0_wp
    ENDDO
!
!-- Horizontal surfaces: natural type
    ALLOCATE( surf_lsm_h%answs( 1:surf_lsm_h%ns, nbins_aerosol ) )
    ALLOCATE( surf_lsm_h%amsws( 1:surf_lsm_h%ns, nbins_aerosol*ncomponents_mass ) )
    surf_lsm_h%answs = 0.0_wp
    surf_lsm_h%amsws = 0.0_wp
!
!-- Horizontal surfaces: urban type
    ALLOCATE( surf_usm_h%answs( 1:surf_usm_h%ns, nbins_aerosol ) )
    ALLOCATE( surf_usm_h%amsws( 1:surf_usm_h%ns, nbins_aerosol*ncomponents_mass ) )
    surf_usm_h%answs = 0.0_wp
    surf_usm_h%amsws = 0.0_wp

!
!-- Vertical surfaces: northward (l=0), southward (l=1), eastward (l=2) and westward (l=3) facing
    DO  l = 0, 3
       ALLOCATE( surf_def_v(l)%answs( 1:surf_def_v(l)%ns, nbins_aerosol ) )
       surf_def_v(l)%answs = 0.0_wp
       ALLOCATE( surf_def_v(l)%amsws( 1:surf_def_v(l)%ns, nbins_aerosol*ncomponents_mass ) )
       surf_def_v(l)%amsws = 0.0_wp

       ALLOCATE( surf_lsm_v(l)%answs( 1:surf_lsm_v(l)%ns, nbins_aerosol ) )
       surf_lsm_v(l)%answs = 0.0_wp
       ALLOCATE( surf_lsm_v(l)%amsws( 1:surf_lsm_v(l)%ns, nbins_aerosol*ncomponents_mass ) )
       surf_lsm_v(l)%amsws = 0.0_wp

       ALLOCATE( surf_usm_v(l)%answs( 1:surf_usm_v(l)%ns, nbins_aerosol ) )
       surf_usm_v(l)%answs = 0.0_wp
       ALLOCATE( surf_usm_v(l)%amsws( 1:surf_usm_v(l)%ns, nbins_aerosol*ncomponents_mass ) )
       surf_usm_v(l)%amsws = 0.0_wp

    ENDDO

!
!-- Concentration of gaseous tracers (1. SO4, 2. HNO3, 3. NH3, 4. OCNV, 5. OCSV)
!-- (number concentration (#/m3) )
!
!-- If chemistry is on, read gas phase concentrations from there. Otherwise, 
!-- allocate salsa_gas array.

    IF ( air_chemistry )  THEN
       DO  lsp = 1, nvar
          SELECT CASE ( TRIM( chem_species(lsp)%name ) )
             CASE ( 'H2SO4', 'h2so4' )
                gases_available = gases_available + 1
                gas_index_chem(1) = lsp
             CASE ( 'HNO3', 'hno3' )
                gases_available = gases_available + 1
                gas_index_chem(2) = lsp
             CASE ( 'NH3', 'nh3' )
                gases_available = gases_available + 1
                gas_index_chem(3) = lsp
             CASE ( 'OCNV', 'ocnv' )
                gases_available = gases_available + 1
                gas_index_chem(4) = lsp
             CASE ( 'OCSV', 'ocsv' )
                gases_available = gases_available + 1
                gas_index_chem(5) = lsp
          END SELECT
       ENDDO

       IF ( gases_available == ngases_salsa )  THEN
          salsa_gases_from_chem = .TRUE.
       ELSE
          WRITE( message_string, * ) 'SALSA is run together with chemistry but not all gaseous '// &
                                     'components are provided by kpp (H2SO4, HNO3, NH3, OCNV, OCSV)'
       CALL message( 'check_parameters', 'PA0599', 1, 2, 0, 6, 0 )
       ENDIF

    ELSE

       ALLOCATE( salsa_gas(ngases_salsa) )
       ALLOCATE( gconc_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg,ngases_salsa),                 &
                 gconc_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg,ngases_salsa),                 &
                 gconc_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg,ngases_salsa) )
       gconc_1 = 0.0_wp
       gconc_2 = 0.0_wp
       gconc_3 = 0.0_wp

       DO i = 1, ngases_salsa
          salsa_gas(i)%conc(nzb:nzt+1,nysg:nyng,nxlg:nxrg)    => gconc_1(:,:,:,i)
          salsa_gas(i)%conc_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg)  => gconc_2(:,:,:,i)
          salsa_gas(i)%tconc_m(nzb:nzt+1,nysg:nyng,nxlg:nxrg) => gconc_3(:,:,:,i)
          ALLOCATE( salsa_gas(i)%flux_s(nzb+1:nzt,0:threads_per_task-1),       &
                    salsa_gas(i)%diss_s(nzb+1:nzt,0:threads_per_task-1),       &
                    salsa_gas(i)%flux_l(nzb+1:nzt,nys:nyn,0:threads_per_task-1),&
                    salsa_gas(i)%diss_l(nzb+1:nzt,nys:nyn,0:threads_per_task-1),&
                    salsa_gas(i)%init(nzb:nzt+1),                              &
                    salsa_gas(i)%sums_ws_l(nzb:nzt+1,0:threads_per_task-1) )
          salsa_gas(i)%init = nclim
          IF ( include_emission )  THEN
             ALLOCATE( salsa_gas(i)%source(nys:nys,nxl:nxr) )
             salsa_gas(i)%source = 0.0_wp
          ENDIF
       ENDDO
!
!--    Surface fluxes: gtsws = gaseous tracer flux
!
!--    Horizontal surfaces: default type
       DO  l = 0, 2   ! upward (l=0), downward (l=1) and model top (l=2)
          ALLOCATE( surf_def_h(l)%gtsws( 1:surf_def_h(l)%ns, ngases_salsa ) )
          surf_def_h(l)%gtsws = 0.0_wp
       ENDDO
!--    Horizontal surfaces: natural type
       ALLOCATE( surf_lsm_h%gtsws( 1:surf_lsm_h%ns, ngases_salsa ) )
       surf_lsm_h%gtsws = 0.0_wp
!--    Horizontal surfaces: urban type
       ALLOCATE( surf_usm_h%gtsws( 1:surf_usm_h%ns, ngases_salsa ) )
       surf_usm_h%gtsws = 0.0_wp
!
!--    Vertical surfaces: northward (l=0), southward (l=1), eastward (l=2) and 
!--    westward (l=3) facing
       DO  l = 0, 3
          ALLOCATE( surf_def_v(l)%gtsws( 1:surf_def_v(l)%ns, ngases_salsa ) )
          surf_def_v(l)%gtsws = 0.0_wp
          ALLOCATE( surf_lsm_v(l)%gtsws( 1:surf_lsm_v(l)%ns, ngases_salsa ) )
          surf_lsm_v(l)%gtsws = 0.0_wp
          ALLOCATE( surf_usm_v(l)%gtsws( 1:surf_usm_v(l)%ns, ngases_salsa ) )
          surf_usm_v(l)%gtsws = 0.0_wp
       ENDDO
    ENDIF

    IF ( ws_scheme_sca )  THEN

       IF ( salsa )  THEN
          ALLOCATE( sums_salsa_ws_l(nzb:nzt+1,0:threads_per_task-1) )
          sums_salsa_ws_l = 0.0_wp
       ENDIF

    ENDIF
!
!-- Set control flags for decycling only at lateral boundary cores. Within the inner cores the
!-- decycle flag is set to .FALSE.. Even though it does not affect the setting of chemistry boundary
!-- conditions, this flag is used to set advection control flags appropriately.
    decycle_salsa_lr = MERGE( decycle_salsa_lr, .FALSE., nxl == 0  .OR.  nxr == nx )
    decycle_salsa_ns = MERGE( decycle_salsa_ns, .FALSE., nys == 0  .OR.  nyn == ny )
!
!-- Decycling can be applied separately for aerosol variables, while wind and other scalars may have
!-- cyclic or nested boundary conditions. However, large gradients near the boundaries may produce
!-- stationary numerical oscillations near the lateral boundaries when a higher-order scheme is
!-- applied near these boundaries. To get rid-off this, set-up additional flags that control the
!-- order of the scalar advection scheme near the lateral boundaries for passive scalars with
!-- decycling.
    IF ( scalar_advec == 'ws-scheme' )  THEN
       ALLOCATE( salsa_advc_flags_s(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
!
!--    In case of decycling, set Neuman boundary conditions for wall_flags_total_0 bit 31 instead of
!--    cyclic boundary conditions. Bit 31 is used to identify extended degradation zones (please see
!--    the following comment). Note, since several also other modules may access this bit but may
!--    have other boundary conditions, the original value of wall_flags_total_0 bit 31 must not be
!--    modified. Hence, store the boundary conditions directly on salsa_advc_flags_s.
!--    salsa_advc_flags_s will be later overwritten in ws_init_flags_scalar and bit 31 won't be used
!--    to control the numerical order.
!--    Initialize with flag 31 only.
       salsa_advc_flags_s = 0
       salsa_advc_flags_s = MERGE( IBSET( salsa_advc_flags_s, 31 ), 0, BTEST( wall_flags_total_0, 31 ) )

       IF ( decycle_salsa_ns )  THEN
          IF ( nys == 0 )  THEN
             DO  i = 1, nbgp
                salsa_advc_flags_s(:,nys-i,:) = MERGE( IBSET( salsa_advc_flags_s(:,nys,:), 31 ),   &
                                                       IBCLR( salsa_advc_flags_s(:,nys,:), 31 ),   &
                                                       BTEST( salsa_advc_flags_s(:,nys,:), 31 ) )
             ENDDO
          ENDIF
          IF ( nyn == ny )  THEN
             DO  i = 1, nbgp
                salsa_advc_flags_s(:,nyn+i,:) = MERGE( IBSET( salsa_advc_flags_s(:,nyn,:), 31 ),   &
                                                       IBCLR( salsa_advc_flags_s(:,nyn,:), 31 ),   &
                                                       BTEST( salsa_advc_flags_s(:,nyn,:), 31 ) )
             ENDDO
          ENDIF
       ENDIF
       IF ( decycle_salsa_lr )  THEN
          IF ( nxl == 0 )  THEN
             DO  i = 1, nbgp
                salsa_advc_flags_s(:,:,nxl-i) = MERGE( IBSET( salsa_advc_flags_s(:,:,nxl), 31 ),   &
                                                       IBCLR( salsa_advc_flags_s(:,:,nxl), 31 ),   &
                                                       BTEST( salsa_advc_flags_s(:,:,nxl), 31 ) )
             ENDDO
          ENDIF
          IF ( nxr == nx )  THEN
             DO  i = 1, nbgp
                salsa_advc_flags_s(:,:,nxr+i) = MERGE( IBSET( salsa_advc_flags_s(:,:,nxr), 31 ),   &
                                                       IBCLR( salsa_advc_flags_s(:,:,nxr), 31 ),   &
                                                       BTEST( salsa_advc_flags_s(:,:,nxr), 31 ) )
             ENDDO
          ENDIF
       ENDIF
!
!--    To initialise the advection flags appropriately, pass the boundary flags to
!--    ws_init_flags_scalar. The last argument in ws_init_flags_scalar indicates that a passive
!--    scalar is being treated and the horizontal advection terms are degraded already 2 grid points
!--    before the lateral boundary. Also, extended degradation zones are applied, where
!--    horizontal advection of scalars is discretised by the first-order scheme at all grid points
!--    in the vicinity of buildings (<= 3 grid points). Even though no building is within the
!--    numerical stencil, the first-order scheme is used. At fourth and fifth grid points, the order
!--    of the horizontal advection scheme is successively upgraded.
!--    These degradations of the advection scheme are done to avoid stationary numerical
!--    oscillations, which are responsible for high concentration maxima that may appear e.g. under
!--    shear-free stable conditions.
       CALL ws_init_flags_scalar( bc_dirichlet_l  .OR.  bc_radiation_l  .OR.  decycle_salsa_lr,    &
                                  bc_dirichlet_n  .OR.  bc_radiation_n  .OR.  decycle_salsa_ns,    &
                                  bc_dirichlet_r  .OR.  bc_radiation_r  .OR.  decycle_salsa_lr,    &
                                  bc_dirichlet_s  .OR.  bc_radiation_s  .OR.  decycle_salsa_ns,    &
                                  salsa_advc_flags_s, .TRUE. )
    ENDIF


 END SUBROUTINE salsa_init_arrays

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of SALSA. Based on salsa_initialize in UCLALES-SALSA.
!> Subroutines salsa_initialize, SALSAinit and DiagInitAero in UCLALES-SALSA are
!> also merged here.
!------------------------------------------------------------------------------!
 SUBROUTINE salsa_init

    IMPLICIT NONE

    INTEGER(iwp) :: i   !<
    INTEGER(iwp) :: ib  !< loop index for aerosol number bins
    INTEGER(iwp) :: ic  !< loop index for aerosol mass bins
    INTEGER(iwp) :: ig  !< loop index for gases
    INTEGER(iwp) :: j   !<

    IF ( debug_output )  CALL debug_message( 'salsa_init', 'start' )

    bin_low_limits = 0.0_wp
    k_topo_top     = 0
    nsect          = 0.0_wp
    massacc        = 1.0_wp
!
!-- Initialise
    IF ( nldepo )  sedim_vd = 0.0_wp

    IF ( .NOT. salsa_gases_from_chem )  THEN
       IF ( .NOT. read_restart_data_salsa )  THEN
          salsa_gas(1)%conc = h2so4_init
          salsa_gas(2)%conc = hno3_init
          salsa_gas(3)%conc = nh3_init
          salsa_gas(4)%conc = ocnv_init
          salsa_gas(5)%conc = ocsv_init
       ENDIF
       DO  ig = 1, ngases_salsa
          salsa_gas(ig)%conc_p    = 0.0_wp
          salsa_gas(ig)%tconc_m   = 0.0_wp
          salsa_gas(ig)%flux_s    = 0.0_wp
          salsa_gas(ig)%diss_s    = 0.0_wp
          salsa_gas(ig)%flux_l    = 0.0_wp
          salsa_gas(ig)%diss_l    = 0.0_wp
          salsa_gas(ig)%sums_ws_l = 0.0_wp
          salsa_gas(ig)%conc_p    = salsa_gas(ig)%conc
       ENDDO
!
!--    Set initial value for gas compound tracer
       salsa_gas(1)%init = h2so4_init
       salsa_gas(2)%init = hno3_init
       salsa_gas(3)%init = nh3_init
       salsa_gas(4)%init = ocnv_init
       salsa_gas(5)%init = ocsv_init
    ENDIF
!
!-- Aerosol radius in each bin: dry and wet (m)
    ra_dry = 1.0E-10_wp
!
!-- Initialise location-dependent aerosol size distributions and chemical compositions:
    CALL aerosol_init

!-- Initalisation run of SALSA + calculate the vertical top index of the topography
    DO  i = nxl, nxr
       DO  j = nys, nyn

          k_topo_top(j,i) = MAXLOC( MERGE( 1, 0, BTEST( wall_flags_total_0(:,j,i), 12 ) ), &
                                       DIM = 1 ) - 1

          CALL salsa_driver( i, j, 1 )
          CALL salsa_diagnostics( i, j )
       ENDDO
    ENDDO

    DO  ib = 1, nbins_aerosol
       aerosol_number(ib)%conc_p    = aerosol_number(ib)%conc
       aerosol_number(ib)%tconc_m   = 0.0_wp
       aerosol_number(ib)%flux_s    = 0.0_wp
       aerosol_number(ib)%diss_s    = 0.0_wp
       aerosol_number(ib)%flux_l    = 0.0_wp
       aerosol_number(ib)%diss_l    = 0.0_wp
       aerosol_number(ib)%sums_ws_l = 0.0_wp
    ENDDO
    DO  ic = 1, ncomponents_mass*nbins_aerosol
       aerosol_mass(ic)%conc_p    = aerosol_mass(ic)%conc
       aerosol_mass(ic)%tconc_m   = 0.0_wp
       aerosol_mass(ic)%flux_s    = 0.0_wp
       aerosol_mass(ic)%diss_s    = 0.0_wp
       aerosol_mass(ic)%flux_l    = 0.0_wp
       aerosol_mass(ic)%diss_l    = 0.0_wp
       aerosol_mass(ic)%sums_ws_l = 0.0_wp
    ENDDO
!
!
!-- Initialise the deposition scheme and surface types
    IF ( nldepo )  CALL init_deposition

    IF ( include_emission )  THEN
!
!--    Read in and initialize emissions
       CALL salsa_emission_setup( .TRUE. )
       IF ( .NOT. salsa_gases_from_chem  .AND.  salsa_emission_mode == 'read_from_file' )  THEN
          CALL salsa_gas_emission_setup( .TRUE. )
       ENDIF
    ENDIF
!
!-- Partition and dissolutional growth by gaseous HNO3 and NH3
    IF ( index_no > 0  .AND.  index_nh > 0  .AND.  index_so4 > 0 )  lspartition = .TRUE.

    IF ( debug_output )  CALL debug_message( 'salsa_init', 'end' )

 END SUBROUTINE salsa_init

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initializes particle size distribution grid by calculating size bin limits 
!> and mid-size for *dry* particles in each bin. Called from salsa_initialize 
!> (only at the beginning of simulation).
!> Size distribution described using:
!>   1) moving center method (subranges 1 and 2)
!>      (Jacobson, Atmos. Env., 31, 131-144, 1997)
!>   2) fixed sectional method (subrange 3)
!> Size bins in each subrange are spaced logarithmically
!> based on given subrange size limits and bin number.
!
!> Mona changed 06/2017: Use geometric mean diameter to describe the mean 
!> particle diameter in a size bin, not the arithmeric mean which clearly
!> overestimates the total particle volume concentration. 
!
!> Coded by:
!> Hannele Korhonen (FMI) 2005
!> Harri Kokkola (FMI) 2006
!
!> Bug fixes for box model + updated for the new aerosol datatype:
!> Juha Tonttila (FMI) 2014
!------------------------------------------------------------------------------!
 SUBROUTINE set_sizebins

    IMPLICIT NONE

    INTEGER(iwp) ::  cc  !< running index
    INTEGER(iwp) ::  dd  !< running index

    REAL(wp) ::  ratio_d  !< ratio of the upper and lower diameter of subranges

    aero(:)%dwet     = 1.0E-10_wp
    aero(:)%veqh2o   = 1.0E-10_wp
    aero(:)%numc     = nclim
    aero(:)%core     = 1.0E-10_wp
    DO  cc = 1, maxspec+1    ! 1:SO4, 2:OC, 3:BC, 4:DU, 5:SS, 6:NO, 7:NH, 8:H2O
       aero(:)%volc(cc) = 0.0_wp
    ENDDO
!
!-- vlolim&vhilim: min & max *dry* volumes [fxm]
!-- dmid: bin mid *dry* diameter (m)
!-- vratiolo&vratiohi: volume ratio between the center and low/high limit
!
!-- 1) Size subrange 1:
    ratio_d = reglim(2) / reglim(1)   ! section spacing (m)
    DO  cc = start_subrange_1a, end_subrange_1a
       aero(cc)%vlolim = api6 * ( reglim(1) * ratio_d**( REAL( cc-1 ) / nbin(1) ) )**3
       aero(cc)%vhilim = api6 * ( reglim(1) * ratio_d**( REAL( cc ) / nbin(1) ) )**3
       aero(cc)%dmid = SQRT( ( aero(cc)%vhilim / api6 )**0.33333333_wp *                           &
                             ( aero(cc)%vlolim / api6 )**0.33333333_wp )
       aero(cc)%vratiohi = aero(cc)%vhilim / ( api6 * aero(cc)%dmid**3 )
       aero(cc)%vratiolo = aero(cc)%vlolim / ( api6 * aero(cc)%dmid**3 )
    ENDDO
!
!-- 2) Size subrange 2:
!-- 2.1) Sub-subrange 2a: high hygroscopicity
    ratio_d = reglim(3) / reglim(2)   ! section spacing
    DO  dd = start_subrange_2a, end_subrange_2a
       cc = dd - start_subrange_2a
       aero(dd)%vlolim = api6 * ( reglim(2) * ratio_d**( REAL( cc ) / nbin(2) ) )**3
       aero(dd)%vhilim = api6 * ( reglim(2) * ratio_d**( REAL( cc+1 ) / nbin(2) ) )**3
       aero(dd)%dmid = SQRT( ( aero(dd)%vhilim / api6 )**0.33333333_wp *                           &
                             ( aero(dd)%vlolim / api6 )**0.33333333_wp )
       aero(dd)%vratiohi = aero(dd)%vhilim / ( api6 * aero(dd)%dmid**3 )
       aero(dd)%vratiolo = aero(dd)%vlolim / ( api6 * aero(dd)%dmid**3 )
    ENDDO
!
!-- 2.2) Sub-subrange 2b: low hygroscopicity
    IF ( .NOT. no_insoluble )  THEN
       aero(start_subrange_2b:end_subrange_2b)%vlolim   = aero(start_subrange_2a:end_subrange_2a)%vlolim
       aero(start_subrange_2b:end_subrange_2b)%vhilim   = aero(start_subrange_2a:end_subrange_2a)%vhilim
       aero(start_subrange_2b:end_subrange_2b)%dmid     = aero(start_subrange_2a:end_subrange_2a)%dmid
       aero(start_subrange_2b:end_subrange_2b)%vratiohi = aero(start_subrange_2a:end_subrange_2a)%vratiohi
       aero(start_subrange_2b:end_subrange_2b)%vratiolo = aero(start_subrange_2a:end_subrange_2a)%vratiolo
    ENDIF
!
!-- Initialize the wet diameter with the bin dry diameter to avoid numerical problems later
    aero(:)%dwet = aero(:)%dmid
!
!-- Save bin limits (lower diameter) to be delivered to PALM if needed
    DO cc = 1, nbins_aerosol
       bin_low_limits(cc) = ( aero(cc)%vlolim / api6 )**0.33333333_wp
    ENDDO

 END SUBROUTINE set_sizebins

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initilize altitude-dependent aerosol size distributions and compositions.
!>
!> Mona added 06/2017: Correct the number and mass concentrations by normalizing 
!< by the given total number and mass concentration.
!>
!> Tomi Raatikainen, FMI, 29.2.2016
!------------------------------------------------------------------------------!
 SUBROUTINE aerosol_init

    USE netcdf_data_input_mod,                                                                     &
        ONLY:  check_existence, close_input_file, get_dimension_length,                            &
               get_attribute, get_variable,                                                        &
               inquire_num_variables, inquire_variable_names,                                      &
               open_read_file

    IMPLICIT NONE

    CHARACTER(LEN=25),  DIMENSION(:), ALLOCATABLE ::  cc_name    !< chemical component name
    CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE ::  var_names  !< variable names

    INTEGER(iwp) ::  ee        !< index: end
    INTEGER(iwp) ::  i         !< loop index: x-direction
    INTEGER(iwp) ::  ib        !< loop index: size bins
    INTEGER(iwp) ::  ic        !< loop index: chemical components
    INTEGER(iwp) ::  id_dyn    !< NetCDF id of PIDS_DYNAMIC_SALSA
    INTEGER(iwp) ::  ig        !< loop index: gases
    INTEGER(iwp) ::  j         !< loop index: y-direction
    INTEGER(iwp) ::  k         !< loop index: z-direction
    INTEGER(iwp) ::  lod_aero  !< level of detail of inital aerosol concentrations
    INTEGER(iwp) ::  num_vars  !< number of variables
    INTEGER(iwp) ::  pr_nbins  !< number of aerosol size bins in file
    INTEGER(iwp) ::  pr_ncc    !< number of aerosol chemical components in file
    INTEGER(iwp) ::  pr_nz     !< number of vertical grid-points in file
    INTEGER(iwp) ::  prunmode  !< running mode of SALSA
    INTEGER(iwp) ::  ss        !< index: start

    INTEGER(iwp), DIMENSION(maxspec) ::  cc_in2mod

    LOGICAL  ::  netcdf_extend = .FALSE. !< Flag: netcdf file exists

    REAL(wp) ::  flag  !< flag to mask topography grid points

    REAL(wp), DIMENSION(nbins_aerosol) ::  core   !< size of the bin mid aerosol particle

    REAL(wp), DIMENSION(0:nz+1) ::  pnf2a   !< number fraction in 2a
    REAL(wp), DIMENSION(0:nz+1) ::  pmfoc1a !< mass fraction of OC in 1a

    REAL(wp), DIMENSION(0:nz+1,nbins_aerosol)   ::  pndist  !< vertical profile of size dist. (#/m3)
    REAL(wp), DIMENSION(0:nz+1,maxspec)         ::  pmf2a   !< mass distributions in subrange 2a
    REAL(wp), DIMENSION(0:nz+1,maxspec)         ::  pmf2b   !< mass distributions in subrange 2b

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  pr_dmid  !< vertical profile of aerosol bin diameters
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  pr_z     !< z levels of profiles

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  pr_mass_fracs_a  !< mass fraction: a
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  pr_mass_fracs_b  !< and b

    cc_in2mod = 0
    prunmode = 1
!
!-- Bin mean aerosol particle volume (m3)
    core(1:nbins_aerosol) = api6 * aero(1:nbins_aerosol)%dmid**3
!
!-- Set concentrations to zero
    pndist(:,:)  = 0.0_wp
    pnf2a(:)     = nf2a
    pmf2a(:,:)   = 0.0_wp
    pmf2b(:,:)   = 0.0_wp
    pmfoc1a(:)   = 0.0_wp

    IF ( init_aerosol_type == 1 )  THEN
!
!--    Read input profiles from PIDS_DYNAMIC_SALSA
#if defined( __netcdf )
!
!--    Location-dependent size distributions and compositions.
       INQUIRE( FILE = TRIM( input_file_dynamic ) //  TRIM( coupling_char ), EXIST = netcdf_extend )
       IF ( netcdf_extend )  THEN
!
!--       Open file in read-only mode
          CALL open_read_file( TRIM( input_file_dynamic ) // TRIM( coupling_char ), id_dyn )
!
!--       At first, inquire all variable names
          CALL inquire_num_variables( id_dyn, num_vars )
!
!--       Allocate memory to store variable names
          ALLOCATE( var_names(1:num_vars) )
          CALL inquire_variable_names( id_dyn, var_names )
!
!--       Inquire vertical dimension and number of aerosol chemical components
          CALL get_dimension_length( id_dyn, pr_nz, 'z' )
          IF ( pr_nz /= nz )  THEN
             WRITE( message_string, * ) 'Number of inifor horizontal grid points does not match '//&
                                        'the number of numeric grid points.'
             CALL message( 'aerosol_init', 'PA0601', 1, 2, 0, 6, 0 )
          ENDIF
          CALL get_dimension_length( id_dyn, pr_ncc, 'composition_index' )
!
!--       Allocate memory
          ALLOCATE( pr_z(1:pr_nz), pr_mass_fracs_a(nzb:nzt+1,pr_ncc),                              &
                    pr_mass_fracs_b(nzb:nzt+1,pr_ncc) )
          pr_mass_fracs_a = 0.0_wp
          pr_mass_fracs_b = 0.0_wp
!
!--       Read vertical levels
          CALL get_variable( id_dyn, 'z', pr_z )
!
!--       Read the names of chemical components
          IF ( check_existence( var_names, 'composition_name' ) )  THEN
             CALL get_variable( id_dyn, 'composition_name', cc_name, pr_ncc )
          ELSE
             WRITE( message_string, * ) 'Missing composition_name in ' // TRIM( input_file_dynamic )
             CALL message( 'aerosol_init', 'PA0655', 1, 2, 0, 6, 0 )
          ENDIF
!
!--       Define the index of each chemical component in the model
          DO  ic = 1, pr_ncc
             SELECT CASE ( TRIM( cc_name(ic) ) )
                CASE ( 'H2SO4', 'SO4', 'h2so4', 'so4' )
                   cc_in2mod(1) = ic
                CASE ( 'OC', 'oc' )
                   cc_in2mod(2) = ic
                CASE ( 'BC', 'bc' )
                   cc_in2mod(3) = ic
                CASE ( 'DU', 'du' )
                   cc_in2mod(4) = ic
                CASE ( 'SS', 'ss' )
                   cc_in2mod(5) = ic
                CASE ( 'HNO3', 'hno3', 'NO3', 'no3', 'NO', 'no' )
                   cc_in2mod(6) = ic
                CASE ( 'NH3', 'nh3', 'NH4', 'nh4', 'NH', 'nh' )
                   cc_in2mod(7) = ic
             END SELECT
          ENDDO

          IF ( SUM( cc_in2mod ) == 0 )  THEN
             message_string = 'None of the aerosol chemical components in ' // TRIM(               &
                              input_file_dynamic ) // ' correspond to ones applied in SALSA.'
             CALL message( 'salsa_mod: aerosol_init', 'PA0602', 2, 2, 0, 6, 0 )
          ENDIF
!
!--       Vertical profiles of mass fractions of different chemical components:
          IF ( check_existence( var_names, 'init_atmosphere_mass_fracs_a' ) )  THEN
             CALL get_variable( id_dyn, 'init_atmosphere_mass_fracs_a', pr_mass_fracs_a,           &
                                0, pr_ncc-1, 0, pr_nz-1 )
          ELSE
             WRITE( message_string, * ) 'Missing init_atmosphere_mass_fracs_a in ' //              &
                                        TRIM( input_file_dynamic )
             CALL message( 'aerosol_init', 'PA0656', 1, 2, 0, 6, 0 )
          ENDIF
          CALL get_variable( id_dyn, 'init_atmosphere_mass_fracs_b', pr_mass_fracs_b,              &
                             0, pr_ncc-1, 0, pr_nz-1  )
!
!--       Match the input data with the chemical composition applied in the model
          DO  ic = 1, maxspec
             ss = cc_in2mod(ic)
             IF ( ss == 0 )  CYCLE
             pmf2a(nzb+1:nzt+1,ic) = pr_mass_fracs_a(nzb:nzt,ss)
             pmf2b(nzb+1:nzt+1,ic) = pr_mass_fracs_b(nzb:nzt,ss)
          ENDDO
!
!--       Aerosol concentrations: lod=1 (vertical profile of sectional number size distribution)
          CALL get_attribute( id_dyn, 'lod', lod_aero, .FALSE., 'init_atmosphere_aerosol' )
          IF ( lod_aero /= 1 )  THEN
             message_string = 'Currently only lod=1 accepted for init_atmosphere_aerosol'
             CALL message( 'salsa_mod: aerosol_init', 'PA0603', 2, 2, 0, 6, 0 )
          ELSE
!
!--          Bin mean diameters in the input file
             CALL get_dimension_length( id_dyn, pr_nbins, 'Dmid')
             IF ( pr_nbins /= nbins_aerosol )  THEN
                message_string = 'Number of size bins in init_atmosphere_aerosol does not match '  &
                                 // 'with that applied in the model'
                CALL message( 'salsa_mod: aerosol_init', 'PA0604', 2, 2, 0, 6, 0 )
             ENDIF

             ALLOCATE( pr_dmid(pr_nbins) )
             pr_dmid    = 0.0_wp

             CALL get_variable( id_dyn, 'Dmid', pr_dmid )
!
!--          Check whether the sectional representation conform to the one 
!--          applied in the model
             IF ( ANY( ABS( ( aero(1:nbins_aerosol)%dmid - pr_dmid ) /                             &
                              aero(1:nbins_aerosol)%dmid )  > 0.1_wp )  ) THEN
                message_string = 'Mean diameters of the aerosol size bins in ' // TRIM(            &
                                 input_file_dynamic ) // ' do not match with the sectional '//     &
                                 'representation of the model.'
                CALL message( 'salsa_mod: aerosol_init', 'PA0605', 2, 2, 0, 6, 0 )
             ENDIF
!
!--          Inital aerosol concentrations
             CALL get_variable( id_dyn, 'init_atmosphere_aerosol', pndist(nzb+1:nzt,:),            &
                                0, pr_nbins-1, 0, pr_nz-1 )
          ENDIF
!
!--       Set bottom and top boundary condition (Neumann)
          pmf2a(nzb,:)    = pmf2a(nzb+1,:)
          pmf2a(nzt+1,:)  = pmf2a(nzt,:)
          pmf2b(nzb,:)    = pmf2b(nzb+1,:)
          pmf2b(nzt+1,:)  = pmf2b(nzt,:)
          pndist(nzb,:)   = pndist(nzb+1,:)
          pndist(nzt+1,:) = pndist(nzt,:)

          IF ( index_so4 < 0 )  THEN
             pmf2a(:,1) = 0.0_wp
             pmf2b(:,1) = 0.0_wp
          ENDIF
          IF ( index_oc < 0 )  THEN
             pmf2a(:,2) = 0.0_wp
             pmf2b(:,2) = 0.0_wp
          ENDIF
          IF ( index_bc < 0 )  THEN
             pmf2a(:,3) = 0.0_wp
             pmf2b(:,3) = 0.0_wp
          ENDIF
          IF ( index_du < 0 )  THEN
             pmf2a(:,4) = 0.0_wp
             pmf2b(:,4) = 0.0_wp
          ENDIF
          IF ( index_ss < 0 )  THEN
             pmf2a(:,5) = 0.0_wp
             pmf2b(:,5) = 0.0_wp
          ENDIF
          IF ( index_no < 0 )  THEN
             pmf2a(:,6) = 0.0_wp
             pmf2b(:,6) = 0.0_wp
          ENDIF
          IF ( index_nh < 0 )  THEN
             pmf2a(:,7) = 0.0_wp
             pmf2b(:,7) = 0.0_wp
          ENDIF

          IF ( SUM( pmf2a ) < 0.00001_wp  .AND.  SUM( pmf2b ) < 0.00001_wp )  THEN
             message_string = 'Error in initialising mass fractions of chemical components. ' //   &
                              'Check that all chemical components are included in parameter file!'
             CALL message( 'salsa_mod: aerosol_init', 'PA0606', 2, 2, 0, 6, 0 ) 
          ENDIF
!
!--       Then normalise the mass fraction so that SUM = 1
          DO  k = nzb, nzt+1
             pmf2a(k,:) = pmf2a(k,:) / SUM( pmf2a(k,:) )
             IF ( SUM( pmf2b(k,:) ) > 0.0_wp )  pmf2b(k,:) = pmf2b(k,:) / SUM( pmf2b(k,:) )
          ENDDO

          DEALLOCATE( pr_z, pr_mass_fracs_a, pr_mass_fracs_b )
!
!--       Close input file
          CALL close_input_file( id_dyn )

       ELSE
          message_string = 'Input file '// TRIM( input_file_dynamic ) // TRIM( coupling_char ) //  &
                           ' for SALSA missing!'
          CALL message( 'salsa_mod: aerosol_init', 'PA0607', 1, 2, 0, 6, 0 )

       ENDIF   ! netcdf_extend

#else
       message_string = 'init_aerosol_type = 1 but preprocessor directive __netcdf is not used '// &
                        'in compiling!'
       CALL message( 'salsa_mod: aerosol_init', 'PA0608', 1, 2, 0, 6, 0 )

#endif

    ELSEIF ( init_aerosol_type == 0 )  THEN
!
!--    Mass fractions for species in a and b-bins
       IF ( index_so4 > 0 )  THEN
          pmf2a(:,1) = mass_fracs_a(index_so4)
          pmf2b(:,1) = mass_fracs_b(index_so4)
       ENDIF
       IF ( index_oc > 0 )  THEN
          pmf2a(:,2) = mass_fracs_a(index_oc)
          pmf2b(:,2) = mass_fracs_b(index_oc)
       ENDIF
       IF ( index_bc > 0 )  THEN
          pmf2a(:,3) = mass_fracs_a(index_bc)
          pmf2b(:,3) = mass_fracs_b(index_bc)
       ENDIF
       IF ( index_du > 0 )  THEN
          pmf2a(:,4) = mass_fracs_a(index_du)
          pmf2b(:,4) = mass_fracs_b(index_du)
       ENDIF
       IF ( index_ss > 0 )  THEN
          pmf2a(:,5) = mass_fracs_a(index_ss)
          pmf2b(:,5) = mass_fracs_b(index_ss)
       ENDIF
       IF ( index_no > 0 )  THEN
          pmf2a(:,6) = mass_fracs_a(index_no)
          pmf2b(:,6) = mass_fracs_b(index_no)
       ENDIF
       IF ( index_nh > 0 )  THEN
          pmf2a(:,7) = mass_fracs_a(index_nh)
          pmf2b(:,7) = mass_fracs_b(index_nh)
       ENDIF
       DO  k = nzb, nzt+1
          pmf2a(k,:) = pmf2a(k,:) / SUM( pmf2a(k,:) )
          IF ( SUM( pmf2b(k,:) ) > 0.0_wp ) pmf2b(k,:) = pmf2b(k,:) / SUM( pmf2b(k,:) )
       ENDDO

       CALL size_distribution( n_lognorm, dpg, sigmag, nsect )
!
!--    Normalize by the given total number concentration
       nsect = nsect * SUM( n_lognorm ) / SUM( nsect )
       DO  ib = start_subrange_1a, end_subrange_2b
          pndist(:,ib) = nsect(ib)
       ENDDO
    ENDIF

    IF ( init_gases_type == 1 )  THEN
!
!--    Read input profiles from PIDS_CHEM
#if defined( __netcdf )
!
!--    Location-dependent size distributions and compositions.
       INQUIRE( FILE = TRIM( input_file_dynamic ) //  TRIM( coupling_char ), EXIST = netcdf_extend )
       IF ( netcdf_extend  .AND.  .NOT. salsa_gases_from_chem )  THEN
!
!--       Open file in read-only mode
          CALL open_read_file( TRIM( input_file_dynamic ) // TRIM( coupling_char ), id_dyn )
!
!--       Inquire dimensions:
          CALL get_dimension_length( id_dyn, pr_nz, 'z' )
          IF ( pr_nz /= nz )  THEN
             WRITE( message_string, * ) 'Number of inifor horizontal grid points does not match '//&
                                        'the number of numeric grid points.'
             CALL message( 'aerosol_init', 'PA0609', 1, 2, 0, 6, 0 )
          ENDIF
!
!--       Read vertical profiles of gases:
          CALL get_variable( id_dyn, 'init_atmosphere_H2SO4', salsa_gas(1)%init(nzb+1:nzt) )
          CALL get_variable( id_dyn, 'init_atmosphere_HNO3',  salsa_gas(2)%init(nzb+1:nzt) )
          CALL get_variable( id_dyn, 'init_atmosphere_NH3',   salsa_gas(3)%init(nzb+1:nzt) )
          CALL get_variable( id_dyn, 'init_atmosphere_OCNV',  salsa_gas(4)%init(nzb+1:nzt) )
          CALL get_variable( id_dyn, 'init_atmosphere_OCSV',  salsa_gas(5)%init(nzb+1:nzt) )
!
!--       Set Neumann top and surface boundary condition for initial + initialise concentrations
          DO  ig = 1, ngases_salsa
             salsa_gas(ig)%init(nzb)   =  salsa_gas(ig)%init(nzb+1)
             salsa_gas(ig)%init(nzt+1) =  salsa_gas(ig)%init(nzt)
             IF ( .NOT. read_restart_data_salsa )  THEN
                DO  k = nzb, nzt+1
                   salsa_gas(ig)%conc(k,:,:) = salsa_gas(ig)%init(k)
                ENDDO
             ENDIF
          ENDDO
!
!--       Close input file
          CALL close_input_file( id_dyn )

       ELSEIF ( .NOT. netcdf_extend  .AND.  .NOT.  salsa_gases_from_chem )  THEN
          message_string = 'Input file '// TRIM( input_file_dynamic ) // TRIM( coupling_char ) //  &
                           ' for SALSA missing!'
          CALL message( 'salsa_mod: aerosol_init', 'PA0610', 1, 2, 0, 6, 0 )

       ENDIF   ! netcdf_extend
#else
       message_string = 'init_gases_type = 1 but preprocessor directive __netcdf is not used in '//&
                        'compiling!'
       CALL message( 'salsa_mod: aerosol_init', 'PA0611', 1, 2, 0, 6, 0 )

#endif

    ENDIF
!
!-- Both SO4 and OC are included, so use the given mass fractions
    IF ( index_oc > 0  .AND.  index_so4 > 0 )  THEN
       pmfoc1a(:) = pmf2a(:,2) / ( pmf2a(:,2) + pmf2a(:,1) )  ! Normalize
!
!-- Pure organic carbon
    ELSEIF ( index_oc > 0 )  THEN
       pmfoc1a(:) = 1.0_wp
!
!-- Pure SO4
    ELSEIF ( index_so4 > 0 )  THEN
       pmfoc1a(:) = 0.0_wp

    ELSE
       message_string = 'Either OC or SO4 must be active for aerosol region 1a!'
       CALL message( 'salsa_mod: aerosol_init', 'PA0612', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Initialize concentrations
    DO  i = nxlg, nxrg
       DO  j = nysg, nyng
          DO  k = nzb, nzt+1
!
!--          Predetermine flag to mask topography
             flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )
!
!--          a) Number concentrations
!--          Region 1:
             DO  ib = start_subrange_1a, end_subrange_1a
                IF ( .NOT. read_restart_data_salsa )  THEN
                   aerosol_number(ib)%conc(k,j,i) = pndist(k,ib) * flag
                ENDIF
                IF ( prunmode == 1 )  THEN
                   aerosol_number(ib)%init = pndist(:,ib)
                ENDIF
             ENDDO
!
!--          Region 2:
             IF ( nreg > 1 )  THEN
                DO  ib = start_subrange_2a, end_subrange_2a
                   IF ( .NOT. read_restart_data_salsa )  THEN
                      aerosol_number(ib)%conc(k,j,i) = MAX( 0.0_wp, pnf2a(k) ) * pndist(k,ib) * flag
                   ENDIF
                   IF ( prunmode == 1 )  THEN
                      aerosol_number(ib)%init = MAX( 0.0_wp, nf2a ) * pndist(:,ib)
                   ENDIF
                ENDDO
                IF ( .NOT. no_insoluble )  THEN
                   DO  ib = start_subrange_2b, end_subrange_2b
                      IF ( pnf2a(k) < 1.0_wp )  THEN
                         IF ( .NOT. read_restart_data_salsa )  THEN
                            aerosol_number(ib)%conc(k,j,i) = MAX( 0.0_wp, 1.0_wp - pnf2a(k) ) *    &
                                                             pndist(k,ib) * flag
                         ENDIF
                         IF ( prunmode == 1 )  THEN
                            aerosol_number(ib)%init = MAX( 0.0_wp, 1.0_wp - nf2a ) * pndist(:,ib)
                         ENDIF
                      ENDIF
                   ENDDO
                ENDIF
             ENDIF
!
!--          b) Aerosol mass concentrations
!--             bin subrange 1: done here separately due to the SO4/OC convention
!
!--          SO4:
             IF ( index_so4 > 0 )  THEN
                ss = ( index_so4 - 1 ) * nbins_aerosol + start_subrange_1a !< start
                ee = ( index_so4 - 1 ) * nbins_aerosol + end_subrange_1a !< end
                ib = start_subrange_1a
                DO  ic = ss, ee
                   IF ( .NOT. read_restart_data_salsa )  THEN
                      aerosol_mass(ic)%conc(k,j,i) = MAX( 0.0_wp, 1.0_wp - pmfoc1a(k) ) *          &
                                                     pndist(k,ib) * core(ib) * arhoh2so4 * flag
                   ENDIF
                   IF ( prunmode == 1 )  THEN
                      aerosol_mass(ic)%init(k) = MAX( 0.0_wp, 1.0_wp - pmfoc1a(k) ) * pndist(k,ib) &
                                                 * core(ib) * arhoh2so4
                   ENDIF
                   ib = ib+1
                ENDDO
             ENDIF
!
!--          OC:
             IF ( index_oc > 0 ) THEN
                ss = ( index_oc - 1 ) * nbins_aerosol + start_subrange_1a !< start
                ee = ( index_oc - 1 ) * nbins_aerosol + end_subrange_1a !< end
                ib = start_subrange_1a
                DO  ic = ss, ee
                   IF ( .NOT. read_restart_data_salsa )  THEN
                      aerosol_mass(ic)%conc(k,j,i) = MAX( 0.0_wp, pmfoc1a(k) ) * pndist(k,ib) *    &
                                                     core(ib) * arhooc * flag
                   ENDIF
                   IF ( prunmode == 1 )  THEN
                      aerosol_mass(ic)%init(k) = MAX( 0.0_wp, pmfoc1a(k) ) * pndist(k,ib) *        &
                                                 core(ib) * arhooc
                   ENDIF
                   ib = ib+1
                ENDDO 
             ENDIF
          ENDDO !< k

          prunmode = 3  ! Init only once

       ENDDO !< j
    ENDDO !< i

!
!-- c) Aerosol mass concentrations
!--    bin subrange 2:
    IF ( nreg > 1 ) THEN

       IF ( index_so4 > 0 ) THEN
          CALL set_aero_mass( index_so4, pmf2a(:,1), pmf2b(:,1), pnf2a, pndist, core, arhoh2so4 )
       ENDIF
       IF ( index_oc > 0 ) THEN
          CALL set_aero_mass( index_oc, pmf2a(:,2), pmf2b(:,2), pnf2a, pndist, core, arhooc )
       ENDIF
       IF ( index_bc > 0 ) THEN
          CALL set_aero_mass( index_bc, pmf2a(:,3), pmf2b(:,3), pnf2a, pndist, core, arhobc )
       ENDIF
       IF ( index_du > 0 ) THEN
          CALL set_aero_mass( index_du, pmf2a(:,4), pmf2b(:,4), pnf2a, pndist, core, arhodu )
       ENDIF
       IF ( index_ss > 0 ) THEN
          CALL set_aero_mass( index_ss, pmf2a(:,5), pmf2b(:,5), pnf2a, pndist, core, arhoss )
       ENDIF
       IF ( index_no > 0 ) THEN
          CALL set_aero_mass( index_no, pmf2a(:,6), pmf2b(:,6), pnf2a, pndist, core, arhohno3 )
       ENDIF
       IF ( index_nh > 0 ) THEN
          CALL set_aero_mass( index_nh, pmf2a(:,7), pmf2b(:,7), pnf2a, pndist, core, arhonh3 )
       ENDIF

    ENDIF

 END SUBROUTINE aerosol_init

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Create a lognormal size distribution and discretise to a sectional 
!> representation.
!------------------------------------------------------------------------------!
 SUBROUTINE size_distribution( in_ntot, in_dpg, in_sigma, psd_sect )

    IMPLICIT NONE

    INTEGER(iwp) ::  ib         !< running index: bin
    INTEGER(iwp) ::  iteration  !< running index: iteration

    REAL(wp) ::  d1         !< particle diameter (m, dummy)
    REAL(wp) ::  d2         !< particle diameter (m, dummy)
    REAL(wp) ::  delta_d    !< (d2-d1)/10
    REAL(wp) ::  deltadp    !< bin width
    REAL(wp) ::  dmidi      !< ( d1 + d2 ) / 2

    REAL(wp), DIMENSION(:), INTENT(in) ::  in_dpg    !< geometric mean diameter (m)
    REAL(wp), DIMENSION(:), INTENT(in) ::  in_ntot   !< number conc. (#/m3)
    REAL(wp), DIMENSION(:), INTENT(in) ::  in_sigma  !< standard deviation

    REAL(wp), DIMENSION(:), INTENT(inout) ::  psd_sect  !< sectional size distribution

    DO  ib = start_subrange_1a, end_subrange_2b
       psd_sect(ib) = 0.0_wp
!
!--    Particle diameter at the low limit (largest in the bin) (m)
       d1 = ( aero(ib)%vlolim / api6 )**0.33333333_wp
!
!--    Particle diameter at the high limit (smallest in the bin) (m)
       d2 = ( aero(ib)%vhilim / api6 )**0.33333333_wp
!
!--    Span of particle diameter in a bin (m)
       delta_d = 0.1_wp * ( d2 - d1 )
!
!--    Iterate:
       DO  iteration = 1, 10
          d1 = ( aero(ib)%vlolim / api6 )**0.33333333_wp + ( ib - 1) * delta_d
          d2 = d1 + delta_d
          dmidi = 0.5_wp * ( d1 + d2 )
          deltadp = LOG10( d2 / d1 )
!
!--       Size distribution
!--       in_ntot = total number, total area, or total volume concentration
!--       in_dpg = geometric-mean number, area, or volume diameter
!--       n(k) = number, area, or volume concentration in a bin
          psd_sect(ib) = psd_sect(ib) + SUM( in_ntot * deltadp / ( SQRT( 2.0_wp * pi ) *           &
                        LOG10( in_sigma ) ) * EXP( -LOG10( dmidi / in_dpg )**2.0_wp /              &
                        ( 2.0_wp * LOG10( in_sigma ) ** 2.0_wp ) ) )

       ENDDO
    ENDDO

 END SUBROUTINE size_distribution

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Sets the mass concentrations to aerosol arrays in 2a and 2b.
!>
!> Tomi Raatikainen, FMI, 29.2.2016
!------------------------------------------------------------------------------!
 SUBROUTINE set_aero_mass( ispec, pmf2a, pmf2b, pnf2a, pndist, pcore, prho )

    IMPLICIT NONE

    INTEGER(iwp) ::  ee        !< index: end
    INTEGER(iwp) ::  i         !< loop index
    INTEGER(iwp) ::  ib        !< loop index
    INTEGER(iwp) ::  ic        !< loop index
    INTEGER(iwp) ::  j         !< loop index
    INTEGER(iwp) ::  k         !< loop index
    INTEGER(iwp) ::  prunmode  !< 1 = initialise
    INTEGER(iwp) ::  ss        !< index: start

    INTEGER(iwp), INTENT(in) :: ispec  !< Aerosol species index

    REAL(wp) ::  flag   !< flag to mask topography grid points

    REAL(wp), INTENT(in) ::  prho !< Aerosol density

    REAL(wp), DIMENSION(nbins_aerosol), INTENT(in) ::  pcore !< Aerosol bin mid core volume
    REAL(wp), DIMENSION(0:nz+1), INTENT(in)        ::  pnf2a !< Number fraction for 2a
    REAL(wp), DIMENSION(0:nz+1), INTENT(in)        ::  pmf2a !< Mass distributions for a
    REAL(wp), DIMENSION(0:nz+1), INTENT(in)        ::  pmf2b !< and b bins

    REAL(wp), DIMENSION(0:nz+1,nbins_aerosol), INTENT(in) ::  pndist !< Aerosol size distribution

    prunmode = 1

    DO i = nxlg, nxrg
       DO j = nysg, nyng
          DO k = nzb, nzt+1
!
!--          Predetermine flag to mask topography
             flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) ) 
!
!--          Regime 2a:
             ss = ( ispec - 1 ) * nbins_aerosol + start_subrange_2a
             ee = ( ispec - 1 ) * nbins_aerosol + end_subrange_2a
             ib = start_subrange_2a
             DO ic = ss, ee
                IF ( .NOT. read_restart_data_salsa )  THEN
                   aerosol_mass(ic)%conc(k,j,i) = MAX( 0.0_wp, pmf2a(k) ) * pnf2a(k) * pndist(k,ib)&
                                                  * pcore(ib) * prho * flag
                ENDIF
                IF ( prunmode == 1 )  THEN
                   aerosol_mass(ic)%init(k) = MAX( 0.0_wp, pmf2a(k) ) * pnf2a(k) * pndist(k,ib) *  &
                                              pcore(ib) * prho
                ENDIF
                ib = ib + 1
             ENDDO
!
!--          Regime 2b:
             IF ( .NOT. no_insoluble )  THEN
                ss = ( ispec - 1 ) * nbins_aerosol + start_subrange_2b
                ee = ( ispec - 1 ) * nbins_aerosol + end_subrange_2b
                ib = start_subrange_2a
                DO ic = ss, ee
                   IF ( .NOT. read_restart_data_salsa )  THEN
                      aerosol_mass(ic)%conc(k,j,i) = MAX( 0.0_wp, pmf2b(k) ) * ( 1.0_wp - pnf2a(k))&
                                                     * pndist(k,ib) * pcore(ib) * prho * flag
                   ENDIF
                   IF ( prunmode == 1 )  THEN
                      aerosol_mass(ic)%init(k) = MAX( 0.0_wp, pmf2b(k) ) * ( 1.0_wp - pnf2a(k) ) * &
                                                 pndist(k,ib) * pcore(ib) * prho 
                   ENDIF
                   ib = ib + 1
                ENDDO  ! c

             ENDIF
          ENDDO   ! k

          prunmode = 3  ! Init only once

       ENDDO   ! j
    ENDDO   ! i

 END SUBROUTINE set_aero_mass

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialise the matching between surface types in LSM and deposition models.
!> Do the matching based on Zhang et al. (2001). Atmos. Environ. 35, 549-560
!> (here referred as Z01).
!------------------------------------------------------------------------------!
 SUBROUTINE init_deposition

    USE surface_mod,                                                                               &
        ONLY:  surf_lsm_h, surf_lsm_v, surf_usm_h, surf_usm_v

    IMPLICIT NONE

    INTEGER(iwp) ::  l  !< loop index for vertical surfaces

    LOGICAL :: match_lsm  !< flag to initilise LSM surfaces (if false, initialise USM surfaces)

    IF ( depo_pcm_par == 'zhang2001' )  THEN
       depo_pcm_par_num = 1
    ELSEIF ( depo_pcm_par == 'petroff2010' )  THEN
       depo_pcm_par_num = 2
    ENDIF

    IF ( depo_surf_par == 'zhang2001' )  THEN
       depo_surf_par_num = 1
    ELSEIF ( depo_surf_par == 'petroff2010' )  THEN
       depo_surf_par_num = 2
    ENDIF
!
!-- LSM: Pavement, vegetation and water
    IF ( nldepo_surf  .AND.  land_surface )  THEN
       match_lsm = .TRUE.
       ALLOCATE( lsm_to_depo_h%match_lupg(1:surf_lsm_h%ns),                                         &
                 lsm_to_depo_h%match_luvw(1:surf_lsm_h%ns),                                         &
                 lsm_to_depo_h%match_luww(1:surf_lsm_h%ns) )
       lsm_to_depo_h%match_lupg = 0
       lsm_to_depo_h%match_luvw = 0
       lsm_to_depo_h%match_luww = 0
       CALL match_sm_zhang( surf_lsm_h, lsm_to_depo_h%match_lupg, lsm_to_depo_h%match_luvw,        &
                            lsm_to_depo_h%match_luww, match_lsm )
       DO  l = 0, 3
          ALLOCATE( lsm_to_depo_v(l)%match_lupg(1:surf_lsm_v(l)%ns),                               &
                    lsm_to_depo_v(l)%match_luvw(1:surf_lsm_v(l)%ns),                               &
                    lsm_to_depo_v(l)%match_luww(1:surf_lsm_v(l)%ns) )
          lsm_to_depo_v(l)%match_lupg = 0
          lsm_to_depo_v(l)%match_luvw = 0
          lsm_to_depo_v(l)%match_luww = 0
          CALL match_sm_zhang( surf_lsm_v(l), lsm_to_depo_v(l)%match_lupg,                         &
                               lsm_to_depo_v(l)%match_luvw, lsm_to_depo_v(l)%match_luww, match_lsm )
       ENDDO
    ENDIF
!
!-- USM: Green roofs/walls, wall surfaces and windows
    IF ( nldepo_surf  .AND.  urban_surface )  THEN
       match_lsm = .FALSE.
       ALLOCATE( usm_to_depo_h%match_lupg(1:surf_usm_h%ns),                                        &
                 usm_to_depo_h%match_luvw(1:surf_usm_h%ns),                                        &
                 usm_to_depo_h%match_luww(1:surf_usm_h%ns) )
       usm_to_depo_h%match_lupg = 0
       usm_to_depo_h%match_luvw = 0
       usm_to_depo_h%match_luww = 0
       CALL match_sm_zhang( surf_usm_h, usm_to_depo_h%match_lupg, usm_to_depo_h%match_luvw,        &
                            usm_to_depo_h%match_luww, match_lsm )
       DO  l = 0, 3
          ALLOCATE( usm_to_depo_v(l)%match_lupg(1:surf_usm_v(l)%ns),                               &
                    usm_to_depo_v(l)%match_luvw(1:surf_usm_v(l)%ns),                               &
                    usm_to_depo_v(l)%match_luww(1:surf_usm_v(l)%ns) )
          usm_to_depo_v(l)%match_lupg = 0
          usm_to_depo_v(l)%match_luvw = 0
          usm_to_depo_v(l)%match_luww = 0
          CALL match_sm_zhang( surf_usm_v(l), usm_to_depo_v(l)%match_lupg,                         &
                               usm_to_depo_v(l)%match_luvw, usm_to_depo_v(l)%match_luww, match_lsm )
       ENDDO
    ENDIF

    IF ( nldepo_pcm )  THEN
       SELECT CASE ( depo_pcm_type )
          CASE ( 'evergreen_needleleaf' )
             depo_pcm_type_num = 1
          CASE ( 'evergreen_broadleaf' )
             depo_pcm_type_num = 2
          CASE ( 'deciduous_needleleaf' )
             depo_pcm_type_num = 3
          CASE ( 'deciduous_broadleaf' )
             depo_pcm_type_num = 4
          CASE DEFAULT
             message_string = 'depo_pcm_type not set correctly.'
             CALL message( 'salsa_mod: init_deposition', 'PA0613', 1, 2, 0, 6, 0 )
       END SELECT
    ENDIF

 END SUBROUTINE init_deposition

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Match the surface types in PALM and Zhang et al. 2001 deposition module
!------------------------------------------------------------------------------!
 SUBROUTINE match_sm_zhang( surf, match_pav_green, match_veg_wall, match_wat_win, match_lsm )

    USE surface_mod,                                                           &
        ONLY:  ind_pav_green, ind_veg_wall, ind_wat_win, surf_type

    IMPLICIT NONE

    INTEGER(iwp) ::  m              !< index for surface elements
    INTEGER(iwp) ::  pav_type_palm  !< pavement / green wall type in PALM
    INTEGER(iwp) ::  veg_type_palm  !< vegetation / wall type in PALM
    INTEGER(iwp) ::  wat_type_palm  !< water / window type in PALM

    INTEGER(iwp), DIMENSION(:), INTENT(inout) ::  match_pav_green  !<  matching pavement/green walls
    INTEGER(iwp), DIMENSION(:), INTENT(inout) ::  match_veg_wall   !<  matching vegetation/walls
    INTEGER(iwp), DIMENSION(:), INTENT(inout) ::  match_wat_win    !<  matching water/windows

    LOGICAL, INTENT(in) :: match_lsm  !< flag to initilise LSM surfaces (if false, initialise USM)

    TYPE(surf_type), INTENT(in) :: surf  !< respective surface type

    DO  m = 1, surf%ns
       IF ( match_lsm )  THEN
!
!--       Vegetation (LSM):
          IF ( surf%frac(m,ind_veg_wall) > 0 )  THEN
             veg_type_palm = surf%vegetation_type(m)
             SELECT CASE ( veg_type_palm )
                CASE ( 0 )
                   message_string = 'No vegetation type defined.'
                   CALL message( 'salsa_mod: init_depo_surfaces', 'PA0614', 1, 2, 0, 6, 0 )
                CASE ( 1 )  ! bare soil
                   match_veg_wall(m) = 6  ! grass in Z01
                CASE ( 2 )  ! crops, mixed farming
                   match_veg_wall(m) = 7  !  crops, mixed farming Z01
                CASE ( 3 )  ! short grass
                   match_veg_wall(m) = 6  ! grass in Z01
                CASE ( 4 )  ! evergreen needleleaf trees
                    match_veg_wall(m) = 1  ! evergreen needleleaf trees in Z01
                CASE ( 5 )  ! deciduous needleleaf trees
                   match_veg_wall(m) = 3  ! deciduous needleleaf trees in Z01
                CASE ( 6 )  ! evergreen broadleaf trees
                   match_veg_wall(m) = 2  ! evergreen broadleaf trees in Z01
                CASE ( 7 )  ! deciduous broadleaf trees
                   match_veg_wall(m) = 4  ! deciduous broadleaf trees in Z01
                CASE ( 8 )  ! tall grass
                   match_veg_wall(m) = 6  ! grass in Z01
                CASE ( 9 )  ! desert
                   match_veg_wall(m) = 8  ! desert in Z01
                CASE ( 10 )  ! tundra
                   match_veg_wall(m) = 9  ! tundra in Z01
                CASE ( 11 )  ! irrigated crops
                   match_veg_wall(m) = 7  !  crops, mixed farming Z01
                CASE ( 12 )  ! semidesert
                   match_veg_wall(m) = 8  ! desert in Z01
                CASE ( 13 )  ! ice caps and glaciers
                   match_veg_wall(m) = 12  ! ice cap and glacier in Z01
                CASE ( 14 )  ! bogs and marshes
                   match_veg_wall(m) = 11  ! wetland with plants in Z01
                CASE ( 15 )  ! evergreen shrubs
                   match_veg_wall(m) = 10  ! shrubs and interrupted woodlands in Z01
                CASE ( 16 )  ! deciduous shrubs
                   match_veg_wall(m) = 10  ! shrubs and interrupted woodlands in Z01
                CASE ( 17 )  ! mixed forest/woodland
                   match_veg_wall(m) = 5  ! mixed broadleaf and needleleaf trees in Z01
                CASE ( 18 )  ! interrupted forest
                   match_veg_wall(m) = 10  ! shrubs and interrupted woodlands in Z01
             END SELECT
          ENDIF
!
!--       Pavement (LSM):
          IF ( surf%frac(m,ind_pav_green) > 0 )  THEN
             pav_type_palm = surf%pavement_type(m)
             IF ( pav_type_palm == 0 )  THEN  ! error
                message_string = 'No pavement type defined.'
                CALL message( 'salsa_mod: match_sm_zhang', 'PA0615', 1, 2, 0, 6, 0 )
             ELSE
                match_pav_green(m) = 15  ! urban in Z01
             ENDIF
          ENDIF
!
!--       Water (LSM):
          IF ( surf%frac(m,ind_wat_win) > 0 )  THEN
             wat_type_palm = surf%water_type(m)
             IF ( wat_type_palm == 0 )  THEN  ! error
                message_string = 'No water type defined.'
                CALL message( 'salsa_mod: match_sm_zhang', 'PA0616', 1, 2, 0, 6, 0 )
             ELSEIF ( wat_type_palm == 3 )  THEN
                match_wat_win(m) = 14  ! ocean in Z01
             ELSEIF ( wat_type_palm == 1  .OR.  wat_type_palm == 2 .OR.  wat_type_palm == 4        &
                      .OR.  wat_type_palm == 5  )  THEN
                match_wat_win(m) = 13  ! inland water in Z01
             ENDIF
          ENDIF
       ELSE
!
!--       Wall surfaces (USM):
          IF ( surf%frac(m,ind_veg_wall) > 0 )  THEN
             match_veg_wall(m) = 15  ! urban in Z01
          ENDIF
!
!--       Green walls and roofs (USM):
          IF ( surf%frac(m,ind_pav_green) > 0 )  THEN
             match_pav_green(m) =  6 ! (short) grass in Z01
          ENDIF
!
!--       Windows (USM):
          IF ( surf%frac(m,ind_wat_win) > 0 )  THEN
             match_wat_win(m) = 15  ! urban in Z01
          ENDIF
       ENDIF

    ENDDO

 END SUBROUTINE match_sm_zhang

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Swapping of timelevels
!------------------------------------------------------------------------------!
 SUBROUTINE salsa_swap_timelevel( mod_count )

    IMPLICIT NONE

    INTEGER(iwp) ::  ib   !<
    INTEGER(iwp) ::  ic   !<
    INTEGER(iwp) ::  icc  !<
    INTEGER(iwp) ::  ig   !<

    INTEGER(iwp), INTENT(IN) ::  mod_count  !<

    IF ( time_since_reference_point >= skip_time_do_salsa )  THEN

       SELECT CASE ( mod_count )

          CASE ( 0 )

             DO  ib = 1, nbins_aerosol
                aerosol_number(ib)%conc(nzb:nzt+1,nysg:nyng,nxlg:nxrg)   => nconc_1(:,:,:,ib)
                aerosol_number(ib)%conc_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg) => nconc_2(:,:,:,ib)

                DO  ic = 1, ncomponents_mass
                   icc = ( ic-1 ) * nbins_aerosol + ib
                   aerosol_mass(icc)%conc(nzb:nzt+1,nysg:nyng,nxlg:nxrg)   => mconc_1(:,:,:,icc)
                   aerosol_mass(icc)%conc_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg) => mconc_2(:,:,:,icc)
                ENDDO
             ENDDO

             IF ( .NOT. salsa_gases_from_chem )  THEN
                DO  ig = 1, ngases_salsa
                   salsa_gas(ig)%conc(nzb:nzt+1,nysg:nyng,nxlg:nxrg)   => gconc_1(:,:,:,ig)
                   salsa_gas(ig)%conc_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg) => gconc_2(:,:,:,ig)
                ENDDO
             ENDIF

          CASE ( 1 )

             DO  ib = 1, nbins_aerosol
                aerosol_number(ib)%conc(nzb:nzt+1,nysg:nyng,nxlg:nxrg)   => nconc_2(:,:,:,ib)
                aerosol_number(ib)%conc_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg) => nconc_1(:,:,:,ib)
                DO  ic = 1, ncomponents_mass
                   icc = ( ic-1 ) * nbins_aerosol + ib
                   aerosol_mass(icc)%conc(nzb:nzt+1,nysg:nyng,nxlg:nxrg)   => mconc_2(:,:,:,icc)
                   aerosol_mass(icc)%conc_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg) => mconc_1(:,:,:,icc)
                ENDDO
             ENDDO

             IF ( .NOT. salsa_gases_from_chem )  THEN
                DO  ig = 1, ngases_salsa
                   salsa_gas(ig)%conc(nzb:nzt+1,nysg:nyng,nxlg:nxrg)   => gconc_2(:,:,:,ig)
                   salsa_gas(ig)%conc_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg) => gconc_1(:,:,:,ig)
                ENDDO
             ENDIF

       END SELECT

    ENDIF

 END SUBROUTINE salsa_swap_timelevel


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine reads the respective restart data.
!------------------------------------------------------------------------------!
 SUBROUTINE salsa_rrd_local( k, nxlf, nxlc, nxl_on_file, nxrf, nxrc, nxr_on_file, nynf, nync,      &
                             nyn_on_file, nysf, nysc, nys_on_file, tmp_3d, found )

    USE control_parameters,                                                                        &
        ONLY:  length, restart_string

    IMPLICIT NONE

    INTEGER(iwp) ::  ib              !<
    INTEGER(iwp) ::  ic              !<
    INTEGER(iwp) ::  ig              !<
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

    LOGICAL, INTENT(OUT)  ::  found  !<

    REAL(wp), &
       DIMENSION(nzb:nzt+1,nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) :: tmp_3d   !<

    found = .FALSE.

    IF ( read_restart_data_salsa )  THEN

       SELECT CASE ( restart_string(1:length) )

          CASE ( 'aerosol_mass' )
             DO  ic = 1, ncomponents_mass * nbins_aerosol
                IF ( k == 1 )  READ ( 13 ) tmp_3d
                aerosol_mass(ic)%conc(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                 &
                                                   tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ENDDO
             found = .TRUE.

          CASE ( 'aerosol_number' )
             DO  ib = 1, nbins_aerosol
                IF ( k == 1 )  READ ( 13 ) tmp_3d
                aerosol_number(ib)%conc(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =               &
                                                   tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ENDDO
             found = .TRUE.

          CASE( 'salsa_gases_av' )
             IF ( .NOT. ALLOCATED( salsa_gases_av ) )  THEN
                ALLOCATE( salsa_gases_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg,ngases_salsa) )
             ENDIF
             DO  ig = 1, ngases_salsa
                IF ( k == 1 )  READ ( 13 ) tmp_3d
                salsa_gases_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp,ig) =                     &
                                                   tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ENDDO
             found = .TRUE.

          CASE ( 'ldsa_av' )
             IF ( .NOT. ALLOCATED( ldsa_av ) )  ALLOCATE( ldsa_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             IF ( k == 1 )  READ ( 13 ) tmp_3d
             ldsa_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                  &
                                                   tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             found = .TRUE.

          CASE ( 'mbins_av' )
             IF ( .NOT. ALLOCATED( mbins_av ) )  THEN
                ALLOCATE( mbins_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg,nbins_aerosol) )
             ENDIF
             DO  ib = 1, nbins_aerosol
                IF ( k == 1 )  READ ( 13 ) tmp_3d
                mbins_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp,ib) =                           &
                                                   tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
                found = .TRUE.
             ENDDO

          CASE ( 'nbins_av' )
             IF ( .NOT. ALLOCATED( nbins_av ) )  THEN
                ALLOCATE( nbins_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg,nbins_aerosol) )
             ENDIF
             DO  ib = 1, nbins_aerosol
                IF ( k == 1 )  READ ( 13 ) tmp_3d
                nbins_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp,ib) =                           &
                                                   tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
                found = .TRUE.
             ENDDO

          CASE ( 'ntot_av' )
             IF ( .NOT. ALLOCATED( ntot_av ) )  ALLOCATE( ntot_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             IF ( k == 1 )  READ ( 13 ) tmp_3d
             ntot_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                  &
                                                   tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             found = .TRUE.

          CASE ( 'nufp_av' )
             IF ( .NOT. ALLOCATED( nufp_av ) )  ALLOCATE( nufp_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             IF ( k == 1 )  READ ( 13 ) tmp_3d
             nufp_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                  &
                                                   tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             found = .TRUE.

          CASE ( 'pm01_av' )
             IF ( .NOT. ALLOCATED( pm01_av ) )  ALLOCATE( pm01_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             IF ( k == 1 )  READ ( 13 ) tmp_3d
             pm01_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                  &
                                                   tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             found = .TRUE.

          CASE ( 'pm25_av' )
             IF ( .NOT. ALLOCATED( pm25_av ) )  ALLOCATE( pm25_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             IF ( k == 1 )  READ ( 13 ) tmp_3d
             pm25_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                  &
                                                   tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             found = .TRUE.

          CASE ( 'pm10_av' )
             IF ( .NOT. ALLOCATED( pm10_av ) )  ALLOCATE( pm10_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             IF ( k == 1 )  READ ( 13 ) tmp_3d
             pm10_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                  &
                                                   tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             found = .TRUE.

          CASE ( 's_mass_av' )
             IF ( .NOT. ALLOCATED( s_mass_av ) )  THEN
                ALLOCATE( s_mass_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg,ncomponents_mass) )
             ENDIF
             DO  ic = 1, ncomponents_mass
                IF ( k == 1 )  READ ( 13 ) tmp_3d
                s_mass_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp,ic) =                          &
                                                   tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ENDDO
             found = .TRUE.

          CASE ( 's_h2o_av' )
             IF ( .NOT. ALLOCATED( s_h2o_av ) )  ALLOCATE( s_h2o_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             IF ( k == 1 )  READ ( 13 ) tmp_3d
             s_h2o_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                 &
                                                   tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             found = .TRUE.

          CASE ( 'salsa_gas' )
             IF ( .NOT. salsa_gases_from_chem )  THEN
                DO  ig = 1, ngases_salsa
                   IF ( k == 1 )  READ ( 13 ) tmp_3d
                   salsa_gas(ig)%conc(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                 &
                                                   tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
                ENDDO
                found = .TRUE.
             ENDIF

          CASE DEFAULT
             found = .FALSE.

       END SELECT
    ENDIF

 END SUBROUTINE salsa_rrd_local

!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine writes the respective restart data.
!> Note that the following input variables in PARIN have to be equal between
!> restart runs:
!>    listspec, nbin, nbin2, nf2a, ncc, mass_fracs_a, mass_fracs_b
!------------------------------------------------------------------------------!
 SUBROUTINE salsa_wrd_local

    USE control_parameters,                                                                        &
        ONLY:  write_binary

    IMPLICIT NONE

    INTEGER(iwp) ::  ib   !<
    INTEGER(iwp) ::  ic   !<
    INTEGER(iwp) ::  ig  !<

    IF ( write_binary  .AND.  write_binary_salsa )  THEN

       CALL wrd_write_string( 'aerosol_mass' )
       DO  ic = 1, nbins_aerosol * ncomponents_mass
          WRITE ( 14 )  aerosol_mass(ic)%conc
       ENDDO

       CALL wrd_write_string( 'aerosol_number' )
       DO  ib = 1, nbins_aerosol
          WRITE ( 14 )  aerosol_number(ib)%conc
       ENDDO

       IF (  .NOT. salsa_gases_from_chem )  THEN

          IF ( ALLOCATED( salsa_gases_av ) )  THEN
             CALL wrd_write_string( 'salsa_gases_av' )
             DO  ig = 1, ngases_salsa
                WRITE ( 14 )  salsa_gases_av(:,:,:,ig)
             ENDDO
          ENDIF
       ENDIF

       IF ( ALLOCATED( ldsa_av ) )  THEN
       CALL wrd_write_string( 'ldsa_av' )
       WRITE ( 14 )  ldsa_av
       ENDIF

       IF ( ALLOCATED( mbins_av ) )  THEN
          CALL wrd_write_string( 'mbins_av' )
          DO  ib = 1, nbins_aerosol
             WRITE ( 14 )  mbins_av(:,:,:,ib)
          ENDDO
       ENDIF

       IF ( ALLOCATED( nbins_av ) )  THEN
          CALL wrd_write_string( 'nbins_av' )
          DO  ib = 1, nbins_aerosol
             WRITE ( 14 )  nbins_av(:,:,:,ib)
          ENDDO
       ENDIF

       IF ( ALLOCATED( ldsa_av ) )  THEN
          CALL wrd_write_string( 'ntot_av' )
          WRITE ( 14 )  ntot_av
       ENDIF

       IF ( ALLOCATED( nufp_av ) )  THEN
          CALL wrd_write_string( 'nufp_av' )
          WRITE ( 14 )  nufp_av
       ENDIF

       IF ( ALLOCATED( pm01_av ) )  THEN
          CALL wrd_write_string( 'pm01_av' )
          WRITE ( 14 )  pm01_av
       ENDIF

       IF ( ALLOCATED( pm25_av ) )  THEN
          CALL wrd_write_string( 'pm25_av' )
          WRITE ( 14 )  pm25_av
       ENDIF

       IF ( ALLOCATED( pm10_av ) )  THEN
          CALL wrd_write_string( 'pm10_av' )
          WRITE ( 14 )  pm10_av
       ENDIF

       IF ( ALLOCATED( s_mass_av ) )  THEN
          CALL wrd_write_string( 's_mass_av' )
          DO  ic = 1, ncomponents_mass
             WRITE ( 14 )  s_mass_av(:,:,:,ic)
          ENDDO
       ENDIF

       IF ( ALLOCATED( s_h2o_av ) )  THEN
          CALL wrd_write_string( 's_h2o_av' )
          WRITE ( 14 )  s_h2o_av
       ENDIF

       IF ( .NOT. salsa_gases_from_chem )  THEN
          CALL wrd_write_string( 'salsa_gas' )
          DO  ig = 1, ngases_salsa
             WRITE ( 14 )  salsa_gas(ig)%conc
          ENDDO
       ENDIF

    ENDIF

 END SUBROUTINE salsa_wrd_local

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Performs necessary unit and dimension conversion between the host model and 
!> SALSA module, and calls the main SALSA routine.
!> Partially adobted form the original SALSA boxmodel version.
!> Now takes masses in as kg/kg from LES!! Converted to m3/m3 for SALSA
!> 05/2016 Juha: This routine is still pretty much in its original shape. 
!>               It's dumb as a mule and twice as ugly, so implementation of
!>               an improved solution is necessary sooner or later.
!> Juha Tonttila, FMI, 2014
!> Jaakko Ahola, FMI, 2016
!> Only aerosol processes included, Mona Kurppa, UHel, 2017
!------------------------------------------------------------------------------!
 SUBROUTINE salsa_driver( i, j, prunmode )

    USE arrays_3d,                                                                                 &
        ONLY: pt_p, q_p, u, v, w

    USE plant_canopy_model_mod,                                                                    &
        ONLY: lad_s

    USE surface_mod,                                                                               &
        ONLY:  surf_def_h, surf_def_v, surf_lsm_h, surf_lsm_v, surf_usm_h, surf_usm_v

    IMPLICIT NONE

    INTEGER(iwp) ::  endi    !< end index
    INTEGER(iwp) ::  ib      !< loop index
    INTEGER(iwp) ::  ic      !< loop index
    INTEGER(iwp) ::  ig      !< loop index
    INTEGER(iwp) ::  k_wall  !< vertical index of topography top
    INTEGER(iwp) ::  k       !< loop index
    INTEGER(iwp) ::  l       !< loop index
    INTEGER(iwp) ::  nc_h2o  !< index of H2O in the prtcl index table
    INTEGER(iwp) ::  ss      !< loop index
    INTEGER(iwp) ::  str     !< start index
    INTEGER(iwp) ::  vc      !< default index in prtcl

    INTEGER(iwp), INTENT(in) ::  i         !< loop index
    INTEGER(iwp), INTENT(in) ::  j         !< loop index
    INTEGER(iwp), INTENT(in) ::  prunmode  !< 1: Initialization, 2: Spinup, 3: Regular runtime

    REAL(wp) ::  cw_old  !< previous H2O mixing ratio
    REAL(wp) ::  flag    !< flag to mask topography grid points
    REAL(wp) ::  in_lad  !< leaf area density (m2/m3)
    REAL(wp) ::  in_rh   !< relative humidity
    REAL(wp) ::  zgso4   !< SO4
    REAL(wp) ::  zghno3  !< HNO3
    REAL(wp) ::  zgnh3   !< NH3
    REAL(wp) ::  zgocnv  !< non-volatile OC
    REAL(wp) ::  zgocsv  !< semi-volatile OC

    REAL(wp), DIMENSION(nzb:nzt+1) ::  in_adn  !< air density (kg/m3)
    REAL(wp), DIMENSION(nzb:nzt+1) ::  in_cs   !< H2O sat. vapour conc.
    REAL(wp), DIMENSION(nzb:nzt+1) ::  in_cw   !< H2O vapour concentration
    REAL(wp), DIMENSION(nzb:nzt+1) ::  in_p    !< pressure (Pa)
    REAL(wp), DIMENSION(nzb:nzt+1) ::  in_t    !< temperature (K)
    REAL(wp), DIMENSION(nzb:nzt+1) ::  in_u    !< wind magnitude (m/s)
    REAL(wp), DIMENSION(nzb:nzt+1) ::  kvis    !< kinematic viscosity of air(m2/s)
    REAL(wp), DIMENSION(nzb:nzt+1) ::  ppm_to_nconc  !< Conversion factor from ppm to #/m3

    REAL(wp), DIMENSION(nzb:nzt+1,nbins_aerosol) ::  schmidt_num  !< particle Schmidt number
    REAL(wp), DIMENSION(nzb:nzt+1,nbins_aerosol) ::  vd           !< particle fall seed (m/s)

    TYPE(t_section), DIMENSION(nbins_aerosol) ::  lo_aero   !< additional variable for OpenMP
    TYPE(t_section), DIMENSION(nbins_aerosol) ::  aero_old  !< helper array

    aero_old(:)%numc = 0.0_wp
    in_lad           = 0.0_wp
    in_u             = 0.0_wp
    kvis             = 0.0_wp
    lo_aero          = aero
    schmidt_num      = 0.0_wp
    vd               = 0.0_wp
    zgso4            = nclim
    zghno3           = nclim
    zgnh3            = nclim
    zgocnv           = nclim
    zgocsv           = nclim
!
!-- Aerosol number is always set, but mass can be uninitialized
    DO ib = 1, nbins_aerosol
       lo_aero(ib)%volc(:)  = 0.0_wp
       aero_old(ib)%volc(:) = 0.0_wp
    ENDDO
!
!-- Set the salsa runtime config (How to make this more efficient?)
    CALL set_salsa_runtime( prunmode )
!
!-- Calculate thermodynamic quantities needed in SALSA
    CALL salsa_thrm_ij( i, j, p_ij=in_p, temp_ij=in_t, cw_ij=in_cw, cs_ij=in_cs, adn_ij=in_adn )
!
!-- Magnitude of wind: needed for deposition
    IF ( lsdepo )  THEN
       in_u(nzb+1:nzt) = SQRT( ( 0.5_wp * ( u(nzb+1:nzt,j,i) + u(nzb+1:nzt,j,i+1) ) )**2 +         &
                               ( 0.5_wp * ( v(nzb+1:nzt,j,i) + v(nzb+1:nzt,j+1,i) ) )**2 +         &
                               ( 0.5_wp * ( w(nzb:nzt-1,j,i) + w(nzb+1:nzt,j,  i) ) )**2 )
    ENDIF
!
!-- Calculate conversion factors for gas concentrations
    ppm_to_nconc(:) = for_ppm_to_nconc * in_p(:) / in_t(:)
!
!-- Determine topography-top index on scalar grid
    k_wall = k_topo_top(j,i)

    DO k = nzb+1, nzt
!
!--    Predetermine flag to mask topography
       flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )
!
!--    Wind velocity for dry depositon on vegetation
       IF ( lsdepo_pcm  .AND.  plant_canopy )  THEN
          in_lad = lad_s( MAX( k-k_wall,0 ),j,i)
       ENDIF
!
!--    For initialization and spinup, limit the RH with the parameter rhlim
       IF ( prunmode < 3 ) THEN
          in_cw(k) = MIN( in_cw(k), in_cs(k) * rhlim )
       ELSE
          in_cw(k) = in_cw(k)
       ENDIF
       cw_old = in_cw(k) !* in_adn(k)
!
!--    Set volume concentrations:
!--    Sulphate (SO4) or sulphuric acid H2SO4 
       IF ( index_so4 > 0 )  THEN
          vc = 1
          str = ( index_so4-1 ) * nbins_aerosol + 1    ! start index
          endi = index_so4 * nbins_aerosol             ! end index
          ic = 1
          DO ss = str, endi
             lo_aero(ic)%volc(vc) = aerosol_mass(ss)%conc(k,j,i) / arhoh2so4
             ic = ic+1
          ENDDO
          aero_old(1:nbins_aerosol)%volc(vc) = lo_aero(1:nbins_aerosol)%volc(vc)
       ENDIF
!
!--    Organic carbon (OC) compounds
       IF ( index_oc > 0 )  THEN
          vc = 2
          str = ( index_oc-1 ) * nbins_aerosol + 1
          endi = index_oc * nbins_aerosol
          ic = 1
          DO ss = str, endi
             lo_aero(ic)%volc(vc) = aerosol_mass(ss)%conc(k,j,i) / arhooc
             ic = ic+1
          ENDDO
          aero_old(1:nbins_aerosol)%volc(vc) = lo_aero(1:nbins_aerosol)%volc(vc)
       ENDIF
!
!--    Black carbon (BC)
       IF ( index_bc > 0 )  THEN
          vc = 3
          str = ( index_bc-1 ) * nbins_aerosol + 1 + end_subrange_1a
          endi = index_bc * nbins_aerosol
          ic = 1 + end_subrange_1a
          DO ss = str, endi
             lo_aero(ic)%volc(vc) = aerosol_mass(ss)%conc(k,j,i) / arhobc
             ic = ic+1
          ENDDO
          aero_old(1:nbins_aerosol)%volc(vc) = lo_aero(1:nbins_aerosol)%volc(vc)
       ENDIF
!
!--    Dust (DU)
       IF ( index_du > 0 )  THEN
          vc = 4
          str = ( index_du-1 ) * nbins_aerosol + 1 + end_subrange_1a
          endi = index_du * nbins_aerosol
          ic = 1 + end_subrange_1a
          DO ss = str, endi
             lo_aero(ic)%volc(vc) = aerosol_mass(ss)%conc(k,j,i) / arhodu
             ic = ic+1
          ENDDO
          aero_old(1:nbins_aerosol)%volc(vc) = lo_aero(1:nbins_aerosol)%volc(vc)
       ENDIF
!
!--    Sea salt (SS)
       IF ( index_ss > 0 )  THEN
          vc = 5
          str = ( index_ss-1 ) * nbins_aerosol + 1 + end_subrange_1a
          endi = index_ss * nbins_aerosol
          ic = 1 + end_subrange_1a
          DO ss = str, endi
             lo_aero(ic)%volc(vc) = aerosol_mass(ss)%conc(k,j,i) / arhoss
             ic = ic+1
          ENDDO
          aero_old(1:nbins_aerosol)%volc(vc) = lo_aero(1:nbins_aerosol)%volc(vc)
       ENDIF
!
!--    Nitrate (NO(3-)) or nitric acid HNO3
       IF ( index_no > 0 )  THEN
          vc = 6
          str = ( index_no-1 ) * nbins_aerosol + 1 
          endi = index_no * nbins_aerosol
          ic = 1
          DO ss = str, endi
             lo_aero(ic)%volc(vc) = aerosol_mass(ss)%conc(k,j,i) / arhohno3
             ic = ic+1
          ENDDO
          aero_old(1:nbins_aerosol)%volc(vc) = lo_aero(1:nbins_aerosol)%volc(vc)
       ENDIF
!
!--    Ammonium (NH(4+)) or ammonia NH3
       IF ( index_nh > 0 )  THEN
          vc = 7
          str = ( index_nh-1 ) * nbins_aerosol + 1
          endi = index_nh * nbins_aerosol
          ic = 1
          DO ss = str, endi
             lo_aero(ic)%volc(vc) = aerosol_mass(ss)%conc(k,j,i) / arhonh3
             ic = ic+1
          ENDDO
          aero_old(1:nbins_aerosol)%volc(vc) = lo_aero(1:nbins_aerosol)%volc(vc)
       ENDIF
!
!--    Water (always used)
       nc_h2o = get_index( prtcl,'H2O' )
       vc = 8
       str = ( nc_h2o-1 ) * nbins_aerosol + 1
       endi = nc_h2o * nbins_aerosol
       ic = 1
       IF ( advect_particle_water )  THEN
          DO ss = str, endi
             lo_aero(ic)%volc(vc) = aerosol_mass(ss)%conc(k,j,i) / arhoh2o
             ic = ic+1
          ENDDO
       ELSE
         lo_aero(1:nbins_aerosol)%volc(vc) = mclim
       ENDIF
       aero_old(1:nbins_aerosol)%volc(vc) = lo_aero(1:nbins_aerosol)%volc(vc)
!
!--    Number concentrations (numc) and particle sizes 
!--    (dwet = wet diameter, core = dry volume)
       DO  ib = 1, nbins_aerosol
          lo_aero(ib)%numc = aerosol_number(ib)%conc(k,j,i)
          aero_old(ib)%numc = lo_aero(ib)%numc
          IF ( lo_aero(ib)%numc > nclim )  THEN
             lo_aero(ib)%dwet = ( SUM( lo_aero(ib)%volc(:) ) / lo_aero(ib)%numc / api6 )**0.33333333_wp
             lo_aero(ib)%core = SUM( lo_aero(ib)%volc(1:7) ) / lo_aero(ib)%numc
          ELSE
             lo_aero(ib)%dwet = lo_aero(ib)%dmid
             lo_aero(ib)%core = api6 * ( lo_aero(ib)%dwet )**3
          ENDIF
       ENDDO
!
!--    Calculate the ambient sizes of particles by equilibrating soluble fraction of particles with
!--    water using the ZSR method.
       in_rh = in_cw(k) / in_cs(k)
       IF ( prunmode==1  .OR.  .NOT. advect_particle_water )  THEN
          CALL equilibration( in_rh, in_t(k), lo_aero, .TRUE. )
       ENDIF
!
!--    Gaseous tracer concentrations in #/m3
       IF ( salsa_gases_from_chem )  THEN
!
!--       Convert concentrations in ppm to #/m3
          zgso4  = chem_species(gas_index_chem(1))%conc(k,j,i) * ppm_to_nconc(k)
          zghno3 = chem_species(gas_index_chem(2))%conc(k,j,i) * ppm_to_nconc(k)
          zgnh3  = chem_species(gas_index_chem(3))%conc(k,j,i) * ppm_to_nconc(k)
          zgocnv = chem_species(gas_index_chem(4))%conc(k,j,i) * ppm_to_nconc(k)
          zgocsv = chem_species(gas_index_chem(5))%conc(k,j,i) * ppm_to_nconc(k)
       ELSE
          zgso4  = salsa_gas(1)%conc(k,j,i)
          zghno3 = salsa_gas(2)%conc(k,j,i)
          zgnh3  = salsa_gas(3)%conc(k,j,i)
          zgocnv = salsa_gas(4)%conc(k,j,i)
          zgocsv = salsa_gas(5)%conc(k,j,i)
       ENDIF
!
!--    Calculate aerosol processes:
!--    *********************************************************************************************
!
!--    Coagulation
       IF ( lscoag )   THEN
          CALL coagulation( lo_aero, dt_salsa, in_t(k), in_p(k) )
       ENDIF
!
!--    Condensation
       IF ( lscnd )   THEN
          CALL condensation( lo_aero, zgso4, zgocnv, zgocsv,  zghno3, zgnh3, in_cw(k), in_cs(k),   &
                             in_t(k), in_p(k), dt_salsa, prtcl )
       ENDIF
!
!--    Deposition
       IF ( lsdepo )  THEN
          CALL deposition( lo_aero, in_t(k), in_adn(k), in_u(k), in_lad, kvis(k), schmidt_num(k,:),&
                           vd(k,:) )
       ENDIF
!
!--    Size distribution bin update
       IF ( lsdistupdate )   THEN
          CALL distr_update( lo_aero )
       ENDIF
!--    *********************************************************************************************

       IF ( lsdepo ) sedim_vd(k,j,i,:) = vd(k,:)
!
!--    Calculate changes in concentrations
       DO ib = 1, nbins_aerosol
          aerosol_number(ib)%conc(k,j,i) = aerosol_number(ib)%conc(k,j,i) + ( lo_aero(ib)%numc -   &
                                           aero_old(ib)%numc ) * flag
       ENDDO

       IF ( index_so4 > 0 )  THEN
          vc = 1
          str = ( index_so4-1 ) * nbins_aerosol + 1
          endi = index_so4 * nbins_aerosol
          ic = 1
          DO ss = str, endi
             aerosol_mass(ss)%conc(k,j,i) = aerosol_mass(ss)%conc(k,j,i) + ( lo_aero(ic)%volc(vc) -&
                                            aero_old(ic)%volc(vc) ) * arhoh2so4 * flag
             ic = ic+1
          ENDDO
       ENDIF

       IF ( index_oc > 0 )  THEN
          vc = 2
          str = ( index_oc-1 ) * nbins_aerosol + 1
          endi = index_oc * nbins_aerosol
          ic = 1
          DO ss = str, endi
             aerosol_mass(ss)%conc(k,j,i) = aerosol_mass(ss)%conc(k,j,i) + ( lo_aero(ic)%volc(vc) -&
                                            aero_old(ic)%volc(vc) ) * arhooc * flag
             ic = ic+1
          ENDDO
       ENDIF

       IF ( index_bc > 0 )  THEN
          vc = 3
          str = ( index_bc-1 ) * nbins_aerosol + 1 + end_subrange_1a
          endi = index_bc * nbins_aerosol
          ic = 1 + end_subrange_1a
          DO ss = str, endi
             aerosol_mass(ss)%conc(k,j,i) = aerosol_mass(ss)%conc(k,j,i) + ( lo_aero(ic)%volc(vc) -&
                                            aero_old(ic)%volc(vc) ) * arhobc * flag
             ic = ic+1
          ENDDO
       ENDIF

       IF ( index_du > 0 )  THEN
          vc = 4
          str = ( index_du-1 ) * nbins_aerosol + 1 + end_subrange_1a
          endi = index_du * nbins_aerosol
          ic = 1 + end_subrange_1a
          DO ss = str, endi
             aerosol_mass(ss)%conc(k,j,i) = aerosol_mass(ss)%conc(k,j,i) + ( lo_aero(ic)%volc(vc) -&
                                            aero_old(ic)%volc(vc) ) * arhodu * flag
             ic = ic+1
          ENDDO
       ENDIF

       IF ( index_ss > 0 )  THEN
          vc = 5
          str = ( index_ss-1 ) * nbins_aerosol + 1 + end_subrange_1a
          endi = index_ss * nbins_aerosol
          ic = 1 + end_subrange_1a
          DO ss = str, endi
             aerosol_mass(ss)%conc(k,j,i) = aerosol_mass(ss)%conc(k,j,i) + ( lo_aero(ic)%volc(vc) -&
                                            aero_old(ic)%volc(vc) ) * arhoss * flag
             ic = ic+1
          ENDDO
       ENDIF

       IF ( index_no > 0 )  THEN
          vc = 6
          str = ( index_no-1 ) * nbins_aerosol + 1
          endi = index_no * nbins_aerosol
          ic = 1
          DO ss = str, endi
             aerosol_mass(ss)%conc(k,j,i) = aerosol_mass(ss)%conc(k,j,i) + ( lo_aero(ic)%volc(vc) -&
                                            aero_old(ic)%volc(vc) ) * arhohno3 * flag
             ic = ic+1
          ENDDO
       ENDIF

       IF ( index_nh > 0 )  THEN
          vc = 7
          str = ( index_nh-1 ) * nbins_aerosol + 1
          endi = index_nh * nbins_aerosol
          ic = 1
          DO ss = str, endi
             aerosol_mass(ss)%conc(k,j,i) = aerosol_mass(ss)%conc(k,j,i) + ( lo_aero(ic)%volc(vc) -&
                                            aero_old(ic)%volc(vc) ) * arhonh3 * flag
             ic = ic+1
          ENDDO
       ENDIF

       IF ( advect_particle_water )  THEN
          nc_h2o = get_index( prtcl,'H2O' )
          vc = 8
          str = ( nc_h2o-1 ) * nbins_aerosol + 1
          endi = nc_h2o * nbins_aerosol
          ic = 1
          DO ss = str, endi
             aerosol_mass(ss)%conc(k,j,i) = aerosol_mass(ss)%conc(k,j,i) + ( lo_aero(ic)%volc(vc) -&
                                            aero_old(ic)%volc(vc) ) * arhoh2o * flag
             ic = ic+1
          ENDDO
          IF ( prunmode == 1 )  THEN
             nc_h2o = get_index( prtcl,'H2O' )
             vc = 8
             str = ( nc_h2o-1 ) * nbins_aerosol + 1
             endi = nc_h2o * nbins_aerosol
             ic = 1
             DO ss = str, endi
                aerosol_mass(ss)%init(k) = MAX( aerosol_mass(ss)%init(k), ( lo_aero(ic)%volc(vc) - &
                                                aero_old(ic)%volc(vc) ) * arhoh2o )
                IF ( k == nzb+1 )  THEN
                   aerosol_mass(ss)%init(k-1) = aerosol_mass(ss)%init(k)
                ELSEIF ( k == nzt  )  THEN
                   aerosol_mass(ss)%init(k+1) = aerosol_mass(ss)%init(k)
                   aerosol_mass(ss)%conc(k+1,j,i) = aerosol_mass(ss)%init(k)
                ENDIF
                ic = ic+1
             ENDDO
          ENDIF
       ENDIF
!
!--    Condensation of precursor gases
       IF ( lscndgas )  THEN
          IF ( salsa_gases_from_chem )  THEN
!
!--          SO4 (or H2SO4)
             ig = gas_index_chem(1)
             chem_species(ig)%conc(k,j,i) = chem_species(ig)%conc(k,j,i) + ( zgso4 /               &
                                            ppm_to_nconc(k) - chem_species(ig)%conc(k,j,i) ) * flag
!
!--          HNO3
             ig = gas_index_chem(2)
             chem_species(ig)%conc(k,j,i) = chem_species(ig)%conc(k,j,i) + ( zghno3 /              &
                                            ppm_to_nconc(k) - chem_species(ig)%conc(k,j,i) ) * flag
!
!--          NH3
             ig = gas_index_chem(3)
             chem_species(ig)%conc(k,j,i) = chem_species(ig)%conc(k,j,i) + ( zgnh3 /               &
                                            ppm_to_nconc(k) - chem_species(ig)%conc(k,j,i) ) * flag
!
!--          non-volatile OC
             ig = gas_index_chem(4)
             chem_species(ig)%conc(k,j,i) = chem_species(ig)%conc(k,j,i) + ( zgocnv /              &
                                            ppm_to_nconc(k) - chem_species(ig)%conc(k,j,i) ) * flag
!
!--          semi-volatile OC
             ig = gas_index_chem(5)
             chem_species(ig)%conc(k,j,i) = chem_species(ig)%conc(k,j,i) + ( zgocsv /              &
                                            ppm_to_nconc(k) - chem_species(ig)%conc(k,j,i) ) * flag

          ELSE
!
!--          SO4 (or H2SO4)
             salsa_gas(1)%conc(k,j,i) = salsa_gas(1)%conc(k,j,i) + ( zgso4 -                       &
                                        salsa_gas(1)%conc(k,j,i) ) * flag
!
!--          HNO3
             salsa_gas(2)%conc(k,j,i) = salsa_gas(2)%conc(k,j,i) + ( zghno3 -                      &
                                        salsa_gas(2)%conc(k,j,i) ) * flag
!
!--          NH3
             salsa_gas(3)%conc(k,j,i) = salsa_gas(3)%conc(k,j,i) + ( zgnh3 -                       &
                                        salsa_gas(3)%conc(k,j,i) ) * flag
!
!--          non-volatile OC
             salsa_gas(4)%conc(k,j,i) = salsa_gas(4)%conc(k,j,i) + ( zgocnv -                      &
                                        salsa_gas(4)%conc(k,j,i) ) * flag
!
!--          semi-volatile OC
             salsa_gas(5)%conc(k,j,i) = salsa_gas(5)%conc(k,j,i) + ( zgocsv -                      &
                                        salsa_gas(5)%conc(k,j,i) ) * flag
          ENDIF
       ENDIF
!
!--    Tendency of water vapour mixing ratio is obtained from the change in RH during SALSA run.
!--    This releases heat and changes pt. Assumes no temperature change during SALSA run.
!--    q = r / (1+r), Euler method for integration
!
       IF ( feedback_to_palm )  THEN
          q_p(k,j,i) = q_p(k,j,i) + 1.0_wp / ( in_cw(k) * in_adn(k) + 1.0_wp )**2 *                &
                       ( in_cw(k) - cw_old ) * in_adn(k) * flag
          pt_p(k,j,i) = pt_p(k,j,i) + alv / c_p * ( in_cw(k) - cw_old ) * in_adn(k) / ( in_cw(k) / &
                        in_adn(k) + 1.0_wp )**2 * pt_p(k,j,i) / in_t(k) * flag
       ENDIF

    ENDDO   ! k

!
!-- Set surfaces and wall fluxes due to deposition
    IF ( lsdepo  .AND.  lsdepo_surf  .AND.  prunmode == 3 )  THEN 
       IF ( .NOT. land_surface  .AND.  .NOT. urban_surface )  THEN
          CALL depo_surf( i, j, surf_def_h(0), vd, schmidt_num, kvis, in_u, .TRUE. )
          DO  l = 0, 3
             CALL depo_surf( i, j, surf_def_v(l), vd, schmidt_num, kvis, in_u, .FALSE. )
          ENDDO
       ELSE
          CALL depo_surf( i, j, surf_usm_h, vd, schmidt_num, kvis, in_u, .TRUE., usm_to_depo_h )
          DO  l = 0, 3
             CALL depo_surf( i, j, surf_usm_v(l), vd, schmidt_num, kvis, in_u, .FALSE.,            &
                             usm_to_depo_v(l) )
          ENDDO
          CALL depo_surf( i, j, surf_lsm_h, vd, schmidt_num, kvis, in_u, .TRUE., lsm_to_depo_h )
          DO  l = 0, 3
             CALL depo_surf( i, j, surf_lsm_v(l), vd, schmidt_num, kvis, in_u, .FALSE.,            &
                             lsm_to_depo_v(l) )
          ENDDO
       ENDIF
    ENDIF

    IF ( prunmode < 3 )  THEN
       !$OMP MASTER
       aero = lo_aero
       !$OMP END MASTER
    END IF

 END SUBROUTINE salsa_driver

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Set logical switches according to the salsa_parameters options.
!> Juha Tonttila, FMI, 2014
!> Only aerosol processes included, Mona Kurppa, UHel, 2017
!------------------------------------------------------------------------------!
 SUBROUTINE set_salsa_runtime( prunmode )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(in) ::  prunmode

    SELECT CASE(prunmode)

       CASE(1) !< Initialization
          lscoag       = .FALSE.
          lscnd        = .FALSE.
          lscndgas     = .FALSE.
          lscndh2oae   = .FALSE.
          lsdepo       = .FALSE.
          lsdepo_pcm   = .FALSE.
          lsdepo_surf  = .FALSE.
          lsdistupdate = .TRUE.
          lspartition  = .FALSE.

       CASE(2)  !< Spinup period
          lscoag      = ( .FALSE. .AND. nlcoag   )
          lscnd       = ( .TRUE.  .AND. nlcnd    )
          lscndgas    = ( .TRUE.  .AND. nlcndgas )
          lscndh2oae  = ( .TRUE.  .AND. nlcndh2oae )

       CASE(3)  !< Run
          lscoag       = nlcoag
          lscnd        = nlcnd
          lscndgas     = nlcndgas
          lscndh2oae   = nlcndh2oae
          lsdepo       = nldepo
          lsdepo_pcm   = nldepo_pcm
          lsdepo_surf  = nldepo_surf
          lsdistupdate = nldistupdate
    END SELECT


 END SUBROUTINE set_salsa_runtime
 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculates the absolute temperature (using hydrostatic pressure), saturation 
!> vapour pressure and mixing ratio over water, relative humidity and air 
!> density needed in the SALSA model.
!> NOTE, no saturation adjustment takes place -> the resulting water vapour 
!> mixing ratio can be supersaturated, allowing the microphysical calculations
!> in SALSA.
!
!> Juha Tonttila, FMI, 2014 (original SALSAthrm)
!> Mona Kurppa, UHel, 2017 (adjustment for PALM and only aerosol processes)
!------------------------------------------------------------------------------!
 SUBROUTINE salsa_thrm_ij( i, j, p_ij, temp_ij, cw_ij, cs_ij, adn_ij )

    USE arrays_3d,                                                                                 &
        ONLY: pt, q, zu

    USE basic_constants_and_equations_mod,                                                         &
        ONLY:  barometric_formula, exner_function, ideal_gas_law_rho, magnus

    IMPLICIT NONE

    INTEGER(iwp), INTENT(in) ::  i  !<
    INTEGER(iwp), INTENT(in) ::  j  !<

    REAL(wp) ::  t_surface  !< absolute surface temperature (K)

    REAL(wp), DIMENSION(nzb:nzt+1) ::  e_s  !< saturation vapour pressure over water (Pa)

    REAL(wp), DIMENSION(:), INTENT(inout) ::  adn_ij   !< air density (kg/m3)
    REAL(wp), DIMENSION(:), INTENT(inout) ::  p_ij     !< air pressure (Pa)
    REAL(wp), DIMENSION(:), INTENT(inout) ::  temp_ij  !< air temperature (K)

    REAL(wp), DIMENSION(:), INTENT(inout), OPTIONAL ::  cw_ij  !< water vapour concentration (kg/m3)
    REAL(wp), DIMENSION(:), INTENT(inout), OPTIONAL ::  cs_ij  !< saturation water vap. conc.(kg/m3)
!
!-- Pressure p_ijk (Pa) = hydrostatic pressure
    t_surface = pt_surface * exner_function( surface_pressure * 100.0_wp )
    p_ij(:) = barometric_formula( zu, t_surface, surface_pressure * 100.0_wp )
!
!-- Absolute ambient temperature (K)
    temp_ij(:) = pt(:,j,i) * exner_function( p_ij(:) )
!
!-- Air density
    adn_ij(:) = ideal_gas_law_rho( p_ij(:), temp_ij(:) )
!
!-- Water vapour concentration r_v (kg/m3)
    IF ( PRESENT( cw_ij ) )  THEN
       cw_ij(:) = ( q(:,j,i) / ( 1.0_wp - q(:,j,i) ) ) * adn_ij(:)
    ENDIF
!
!-- Saturation mixing ratio r_s (kg/kg) from vapour pressure at temp (Pa)
    IF ( PRESENT( cs_ij ) )  THEN
       e_s(:) = 611.0_wp * EXP( alv_d_rv * ( 3.6609E-3_wp - 1.0_wp /           &
                temp_ij(:) ) )! magnus( temp_ij(:) )
       cs_ij(:) = ( 0.622_wp * e_s / ( p_ij(:) - e_s(:) ) ) * adn_ij(:)
    ENDIF

 END SUBROUTINE salsa_thrm_ij

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculates ambient sizes of particles by equilibrating soluble fraction of 
!> particles with water using the ZSR method (Stokes and Robinson, 1966).
!> Method:
!> Following chemical components are assumed water-soluble
!> - (ammonium) sulphate (100%)
!> - sea salt (100 %)
!> - organic carbon (epsoc * 100%) 
!> Exact thermodynamic considerations neglected.
!> - If particles contain no sea salt, calculation according to sulphate 
!>   properties
!> - If contain sea salt but no sulphate, calculation according to sea salt 
!>   properties
!> - If contain both sulphate and sea salt -> the molar fraction of these 
!>   compounds determines which one of them is used as the basis of calculation.
!> If sulphate and sea salt coexist in a particle, it is assumed that the Cl is 
!> replaced by sulphate; thus only either sulphate + organics or sea salt + 
!> organics is included in the calculation of soluble fraction.
!> Molality parameterizations taken from Table 1 of Tang: Thermodynamic and 
!> optical properties of mixed-salt aerosols of atmospheric importance, 
!> J. Geophys. Res., 102 (D2), 1883-1893 (1997)
!
!> Coded by:
!> Hannele Korhonen (FMI) 2005 
!> Harri Kokkola (FMI) 2006
!> Matti Niskanen(FMI) 2012
!> Anton Laakso  (FMI) 2013
!> Modified for the new aerosol datatype, Juha Tonttila (FMI) 2014
!
!> fxm: should sea salt form a solid particle when prh is very low (even though 
!> it could be mixed with e.g. sulphate)?
!> fxm: crashes if no sulphate or sea salt
!> fxm: do we really need to consider Kelvin effect for subrange 2
!------------------------------------------------------------------------------!
 SUBROUTINE equilibration( prh, ptemp, paero, init )

    IMPLICIT NONE

    INTEGER(iwp) :: ib      !< loop index
    INTEGER(iwp) :: counti  !< loop index

    LOGICAL, INTENT(in) ::  init   !< TRUE: Initialization, FALSE: Normal runtime: update water
                                   !< content only for 1a

    REAL(wp) ::  zaw      !< water activity [0-1]
    REAL(wp) ::  zcore    !< Volume of dry particle
    REAL(wp) ::  zdold    !< Old diameter
    REAL(wp) ::  zdwet    !< Wet diameter or mean droplet diameter
    REAL(wp) ::  zke      !< Kelvin term in the Köhler equation
    REAL(wp) ::  zlwc     !< liquid water content [kg/m3-air]
    REAL(wp) ::  zrh      !< Relative humidity

    REAL(wp), DIMENSION(maxspec) ::  zbinmol  !< binary molality of each components (mol/kg)
    REAL(wp), DIMENSION(maxspec) ::  zvpart   !< volume of chem. compounds in one particle

    REAL(wp), INTENT(in) ::  prh    !< relative humidity [0-1]
    REAL(wp), INTENT(in) ::  ptemp  !< temperature (K)

    TYPE(t_section), DIMENSION(nbins_aerosol), INTENT(inout) ::  paero  !< aerosol properties

    zaw       = 0.0_wp
    zlwc      = 0.0_wp
!
!-- Relative humidity:
    zrh = prh
    zrh = MAX( zrh, 0.05_wp )
    zrh = MIN( zrh, 0.98_wp)
!
!-- 1) Regime 1: sulphate and partly water-soluble OC. Done for every CALL
    DO  ib = start_subrange_1a, end_subrange_1a   ! size bin

       zbinmol = 0.0_wp
       zdold   = 1.0_wp
       zke     = 1.02_wp

       IF ( paero(ib)%numc > nclim )  THEN
!
!--       Volume in one particle
          zvpart = 0.0_wp
          zvpart(1:2) = paero(ib)%volc(1:2) / paero(ib)%numc
          zvpart(6:7) = paero(ib)%volc(6:7) / paero(ib)%numc
!
!--       Total volume and wet diameter of one dry particle 
          zcore = SUM( zvpart(1:2) )
          zdwet = paero(ib)%dwet

          counti = 0
          DO  WHILE ( ABS( zdwet / zdold - 1.0_wp ) > 1.0E-2_wp )

             zdold = MAX( zdwet, 1.0E-20_wp )
             zaw = MAX( 1.0E-3_wp, zrh / zke ) ! To avoid underflow
!
!--          Binary molalities (mol/kg): 
!--          Sulphate
             zbinmol(1) = 1.1065495E+2_wp - 3.6759197E+2_wp * zaw + 5.0462934E+2_wp * zaw**2 -     &
                          3.1543839E+2_wp * zaw**3 + 6.770824E+1_wp  * zaw**4
!--          Organic carbon
             zbinmol(2) = 1.0_wp / ( zaw * amh2o ) - 1.0_wp / amh2o
!--          Nitric acid
             zbinmol(6) = 2.306844303E+1_wp - 3.563608869E+1_wp * zaw - 6.210577919E+1_wp * zaw**2 &
                          + 5.510176187E+2_wp * zaw**3 - 1.460055286E+3_wp * zaw**4                &
                          + 1.894467542E+3_wp * zaw**5 - 1.220611402E+3_wp * zaw**6                &
                          + 3.098597737E+2_wp * zaw**7
!
!--          Calculate the liquid water content (kg/m3-air) using ZSR (see e.g. Eq. 10.98 in
!--          Seinfeld and Pandis (2006))
             zlwc = ( paero(ib)%volc(1) * ( arhoh2so4 / amh2so4 ) ) / zbinmol(1) +                 &
                    epsoc * paero(ib)%volc(2) * ( arhooc / amoc ) / zbinmol(2) +                   &
                    ( paero(ib)%volc(6) * ( arhohno3/amhno3 ) ) / zbinmol(6)
!
!--          Particle wet diameter (m)
             zdwet = ( zlwc / paero(ib)%numc / arhoh2o / api6 + ( SUM( zvpart(6:7) ) / api6 ) +    &
                       zcore / api6 )**0.33333333_wp
!
!--          Kelvin effect (Eq. 10.85 in in Seinfeld and Pandis (2006)). Avoid 
!--          overflow.
             zke = EXP( MIN( 50.0_wp, 4.0_wp * surfw0 * amvh2so4 / ( abo * ptemp *  zdwet ) ) )

             counti = counti + 1
             IF ( counti > 1000 )  THEN 
                message_string = 'Subrange 1: no convergence!'
                CALL message( 'salsa_mod: equilibration', 'PA0617', 1, 2, 0, 6, 0 )
             ENDIF
          ENDDO
!
!--       Instead of lwc, use the volume concentration of water from now on 
!--       (easy to convert...)
          paero(ib)%volc(8) = zlwc / arhoh2o
!
!--       If this is initialization, update the core and wet diameter
          IF ( init )  THEN
             paero(ib)%dwet = zdwet
             paero(ib)%core = zcore
          ENDIF

       ELSE
!--       If initialization
!--       1.2) empty bins given bin average values
          IF ( init )  THEN
             paero(ib)%dwet = paero(ib)%dmid
             paero(ib)%core = api6 * paero(ib)%dmid**3
          ENDIF

       ENDIF

    ENDDO  ! ib
!
!-- 2) Regime 2a: sulphate, OC, BC and sea salt
!--    This is done only for initialization call, otherwise the water contents
!--    are computed via condensation
    IF ( init )  THEN
       DO  ib = start_subrange_2a, end_subrange_2b
!
!--       Initialize
          zke     = 1.02_wp
          zbinmol = 0.0_wp
          zdold   = 1.0_wp
!
!--       1) Particle properties calculated for non-empty bins
          IF ( paero(ib)%numc > nclim )  THEN
!
!--          Volume in one particle [fxm]
             zvpart = 0.0_wp
             zvpart(1:7) = paero(ib)%volc(1:7) / paero(ib)%numc
!
!--          Total volume and wet diameter of one dry particle [fxm]
             zcore = SUM( zvpart(1:5) )
             zdwet = paero(ib)%dwet

             counti = 0
             DO  WHILE ( ABS( zdwet / zdold - 1.0_wp ) > 1.0E-12_wp )

                zdold = MAX( zdwet, 1.0E-20_wp )
                zaw = zrh / zke
!
!--             Binary molalities (mol/kg):
!--             Sulphate
                zbinmol(1) = 1.1065495E+2_wp - 3.6759197E+2_wp * zaw + 5.0462934E+2_wp * zaw**2 -  &
                             3.1543839E+2_wp * zaw**3 + 6.770824E+1_wp  * zaw**4
!--             Organic carbon
                zbinmol(2) = 1.0_wp / ( zaw * amh2o ) - 1.0_wp / amh2o
!--             Nitric acid
                zbinmol(6) = 2.306844303E+1_wp          - 3.563608869E+1_wp * zaw -                &
                             6.210577919E+1_wp * zaw**2 + 5.510176187E+2_wp * zaw**3 -             &
                             1.460055286E+3_wp * zaw**4 + 1.894467542E+3_wp * zaw**5 -             &
                             1.220611402E+3_wp * zaw**6 + 3.098597737E+2_wp * zaw**7 
!--             Sea salt (natrium chloride)
                zbinmol(5) = 5.875248E+1_wp - 1.8781997E+2_wp * zaw + 2.7211377E+2_wp * zaw**2 -   &
                             1.8458287E+2_wp * zaw**3 + 4.153689E+1_wp  * zaw**4
!
!--             Calculate the liquid water content (kg/m3-air)
                zlwc = ( paero(ib)%volc(1) * ( arhoh2so4 / amh2so4 ) ) / zbinmol(1) +              &
                       epsoc * ( paero(ib)%volc(2) * ( arhooc / amoc ) ) / zbinmol(2) +            &
                       ( paero(ib)%volc(6) * ( arhohno3 / amhno3 ) ) / zbinmol(6) +                &
                       ( paero(ib)%volc(5) * ( arhoss / amss ) ) / zbinmol(5)

!--             Particle wet radius (m)
                zdwet = ( zlwc / paero(ib)%numc / arhoh2o / api6 + ( SUM( zvpart(6:7) ) / api6 )  + &
                           zcore / api6 )**0.33333333_wp
!
!--             Kelvin effect (Eq. 10.85 in Seinfeld and Pandis (2006))
                zke = EXP( MIN( 50.0_wp, 4.0_wp * surfw0 * amvh2so4 / ( abo * zdwet * ptemp ) ) )

                counti = counti + 1
                IF ( counti > 1000 )  THEN 
                   message_string = 'Subrange 2: no convergence!'
                CALL message( 'salsa_mod: equilibration', 'PA0618', 1, 2, 0, 6, 0 )
                ENDIF
             ENDDO
!
!--          Liquid water content; instead of LWC use the volume concentration
             paero(ib)%volc(8) = zlwc / arhoh2o
             paero(ib)%dwet    = zdwet
             paero(ib)%core    = zcore

          ELSE
!--          2.2) empty bins given bin average values
             paero(ib)%dwet = paero(ib)%dmid
             paero(ib)%core = api6 * paero(ib)%dmid**3
          ENDIF

       ENDDO   ! ib
    ENDIF

 END SUBROUTINE equilibration

!------------------------------------------------------------------------------!
!> Description:
!> ------------
!> Calculation of the settling velocity vc (m/s) per aerosol size bin and 
!> deposition on plant canopy (lsdepo_pcm).
!
!> Deposition is based on either the scheme presented in:
!> Zhang et al. (2001), Atmos. Environ. 35, 549-560 (includes collection due to
!> Brownian diffusion, impaction, interception and sedimentation; hereafter ZO1)
!> OR
!> Petroff & Zhang (2010), Geosci. Model Dev. 3, 753-769 (includes also
!> collection due to turbulent impaction, hereafter P10)
!
!> Equation numbers refer to equation in Jacobson (2005): Fundamentals of
!> Atmospheric Modeling, 2nd Edition.
!
!> Subroutine follows closely sedim_SALSA in UCLALES-SALSA written by Juha
!> Tonttila (KIT/FMI) and Zubair Maalick (UEF).
!> Rewritten to PALM by Mona Kurppa (UH), 2017.
!
!> Call for grid point i,j,k
!------------------------------------------------------------------------------!

 SUBROUTINE deposition( paero, tk, adn, mag_u, lad, kvis, schmidt_num, vc )

    USE plant_canopy_model_mod,                                                &
        ONLY:  canopy_drag_coeff

    IMPLICIT NONE

    INTEGER(iwp) ::  ib   !< loop index
    INTEGER(iwp) ::  ic   !< loop index

    REAL(wp) ::  alpha             !< parameter, Table 3 in Z01
    REAL(wp) ::  avis              !< molecular viscocity of air (kg/(m*s))
    REAL(wp) ::  beta_im           !< parameter for turbulent impaction
    REAL(wp) ::  c_brownian_diff   !< coefficient for Brownian diffusion
    REAL(wp) ::  c_impaction       !< coefficient for inertial impaction
    REAL(wp) ::  c_interception    !< coefficient for interception
    REAL(wp) ::  c_turb_impaction  !< coefficient for turbulent impaction
    REAL(wp) ::  depo              !< deposition velocity (m/s)
    REAL(wp) ::  gamma             !< parameter, Table 3 in Z01
    REAL(wp) ::  lambda            !< molecular mean free path (m)
    REAL(wp) ::  mdiff             !< particle diffusivity coefficient
    REAL(wp) ::  par_a             !< parameter A for the characteristic radius of collectors,
                                   !< Table 3 in Z01
    REAL(wp) ::  par_l             !< obstacle characteristic dimension in P10
    REAL(wp) ::  pdn               !< particle density (kg/m3)
    REAL(wp) ::  ustar             !< friction velocity (m/s)
    REAL(wp) ::  va                !< thermal speed of an air molecule (m/s)

    REAL(wp), INTENT(in) ::  adn    !< air density (kg/m3)
    REAL(wp), INTENT(in) ::  lad    !< leaf area density (m2/m3)
    REAL(wp), INTENT(in) ::  mag_u  !< wind velocity (m/s)
    REAL(wp), INTENT(in) ::  tk     !< abs.temperature (K)

    REAL(wp), INTENT(inout) ::  kvis   !< kinematic viscosity of air (m2/s)

    REAL(wp), DIMENSION(nbins_aerosol) ::  beta   !< Cunningham slip-flow correction factor
    REAL(wp), DIMENSION(nbins_aerosol) ::  Kn     !< Knudsen number
    REAL(wp), DIMENSION(nbins_aerosol) ::  zdwet  !< wet diameter (m)

    REAL(wp), DIMENSION(:), INTENT(inout) ::  schmidt_num  !< particle Schmidt number
    REAL(wp), DIMENSION(:), INTENT(inout) ::  vc  !< critical fall speed i.e. settling velocity of
                                                  !< an aerosol particle (m/s)

    TYPE(t_section), DIMENSION(nbins_aerosol), INTENT(inout) ::  paero  !< aerosol properties
!
!-- Initialise
    depo  = 0.0_wp
    pdn   = 1500.0_wp    ! default value
    ustar = 0.0_wp
!
!-- Molecular viscosity of air (Eq. 4.54)
    avis = 1.8325E-5_wp * ( 416.16_wp / ( tk + 120.0_wp ) ) * ( tk / 296.16_wp )**1.5_wp
!
!-- Kinematic viscosity (Eq. 4.55)
    kvis =  avis / adn
!
!-- Thermal velocity of an air molecule (Eq. 15.32)
    va = SQRT( 8.0_wp * abo * tk / ( pi * am_airmol ) )
!
!-- Mean free path (m) (Eq. 15.24)
    lambda = 2.0_wp * avis / ( adn * va )
!
!-- Particle wet diameter (m)
    zdwet = paero(:)%dwet
!
!-- Knudsen number (Eq. 15.23)
    Kn = MAX( 1.0E-2_wp, lambda / ( zdwet * 0.5_wp ) ) ! To avoid underflow
!
!-- Cunningham slip-flow correction (Eq. 15.30)
    beta = 1.0_wp + Kn * ( 1.249_wp + 0.42_wp * EXP( -0.87_wp / Kn ) )
!
!-- Critical fall speed i.e. settling velocity  (Eq. 20.4)
    vc = MIN( 1.0_wp, zdwet**2 * ( pdn - adn ) * g * beta / ( 18.0_wp * avis ) )
!
!-- Deposition on vegetation
    IF ( lsdepo_pcm  .AND.  plant_canopy  .AND.  lad > 0.0_wp )  THEN
!
!--    Parameters for the land use category 'deciduous broadleaf trees'(Table 3)
       alpha   = alpha_z01(depo_pcm_type_num)
       gamma   = gamma_z01(depo_pcm_type_num)
       par_a   = A_z01(depo_pcm_type_num, season_z01) * 1.0E-3_wp
!
!--    Deposition efficiencies from Table 1. Constants from Table 2.
       par_l            = l_p10(depo_pcm_type_num) * 0.01_wp
       c_brownian_diff  = c_b_p10(depo_pcm_type_num)
       c_interception   = c_in_p10(depo_pcm_type_num)
       c_impaction      = c_im_p10(depo_pcm_type_num)
       beta_im          = beta_im_p10(depo_pcm_type_num)
       c_turb_impaction = c_it_p10(depo_pcm_type_num)

       DO  ib = 1, nbins_aerosol

          IF ( paero(ib)%numc < ( 2.0_wp * nclim ) )  CYCLE

!--       Particle diffusivity coefficient (Eq. 15.29)
          mdiff = ( abo * tk * beta(ib) ) / ( 3.0_wp * pi * avis * zdwet(ib) )
!
!--       Particle Schmidt number (Eq. 15.36)
          schmidt_num(ib) = kvis / mdiff
!
!--       Friction velocity for deposition on vegetation. Calculated following Prandtl (1925):
          ustar = SQRT( canopy_drag_coeff ) * mag_u
          SELECT CASE ( depo_pcm_par_num )

             CASE ( 1 )   ! Zhang et al. (2001)
                CALL depo_vel_Z01( vc(ib), ustar, schmidt_num(ib), paero(ib)%dwet, alpha,  gamma,  &
                                   par_a, depo )
             CASE ( 2 )   ! Petroff & Zhang (2010)
                CALL depo_vel_P10( vc(ib), mag_u, ustar, kvis, schmidt_num(ib), paero(ib)%dwet,    &
                                   par_l, c_brownian_diff, c_interception, c_impaction, beta_im,   &
                                   c_turb_impaction, depo )
          END SELECT
!
!--       Calculate the change in concentrations
          paero(ib)%numc = paero(ib)%numc - depo * lad * paero(ib)%numc * dt_salsa
          DO  ic = 1, maxspec+1
             paero(ib)%volc(ic) = paero(ib)%volc(ic) - depo * lad * paero(ib)%volc(ic) * dt_salsa
          ENDDO
       ENDDO

    ENDIF

 END SUBROUTINE deposition

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate deposition velocity (m/s) based on Zhan et al. (2001, case 1).
!------------------------------------------------------------------------------!

 SUBROUTINE depo_vel_Z01( vc, ustar, schmidt_num, diameter, alpha, gamma, par_a, depo )

    IMPLICIT NONE

    REAL(wp) ::  rs                !< overall quasi-laminar resistance for particles
    REAL(wp) ::  stokes_num        !< Stokes number for smooth or bluff surfaces

    REAL(wp), INTENT(in) ::  alpha        !< parameter, Table 3 in Z01
    REAL(wp), INTENT(in) ::  gamma        !< parameter, Table 3 in Z01
    REAL(wp), INTENT(in) ::  par_a        !< parameter A for the characteristic diameter of
                                          !< collectors, Table 3 in Z01
    REAL(wp), INTENT(in) ::  diameter     !< particle diameter
    REAL(wp), INTENT(in) ::  schmidt_num  !< particle Schmidt number
    REAL(wp), INTENT(in) ::  ustar        !< friction velocity (m/s)
    REAL(wp), INTENT(in) ::  vc           !< terminal velocity (m/s)

    REAL(wp), INTENT(inout)  ::  depo     !< deposition efficiency (m/s)

    IF ( par_a > 0.0_wp )  THEN
!
!--    Initialise
       rs = 0.0_wp
!
!--    Stokes number for vegetated surfaces (Seinfeld & Pandis (2006): Eq.19.24)
       stokes_num = vc * ustar / ( g * par_a )
!
!--    The overall quasi-laminar resistance for particles (Zhang et al., Eq. 5)
       rs = MAX( EPSILON( 1.0_wp ), ( 3.0_wp * ustar * EXP( -stokes_num**0.5_wp ) *                &
                 ( schmidt_num**( -gamma ) + ( stokes_num / ( alpha + stokes_num ) )**2 +          &
                 0.5_wp * ( diameter / par_a )**2 ) ) )

       depo = rs + vc

    ELSE
       depo = 0.0_wp
    ENDIF

 END SUBROUTINE depo_vel_Z01

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate deposition velocity (m/s) based on Petroff & Zhang (2010, case 2).
!------------------------------------------------------------------------------!

 SUBROUTINE depo_vel_P10( vc, mag_u, ustar, kvis_a, schmidt_num, diameter, par_l, c_brownian_diff, &
                          c_interception, c_impaction, beta_im, c_turb_impaction, depo )

    IMPLICIT NONE

    REAL(wp) ::  stokes_num        !< Stokes number for smooth or bluff surfaces
    REAL(wp) ::  tau_plus          !< dimensionless particle relaxation time
    REAL(wp) ::  v_bd              !< deposition velocity due to Brownian diffusion
    REAL(wp) ::  v_im              !< deposition velocity due to impaction
    REAL(wp) ::  v_in              !< deposition velocity due to interception
    REAL(wp) ::  v_it              !< deposition velocity due to turbulent impaction

    REAL(wp), INTENT(in) ::  beta_im           !< parameter for turbulent impaction
    REAL(wp), INTENT(in) ::  c_brownian_diff   !< coefficient for Brownian diffusion
    REAL(wp), INTENT(in) ::  c_impaction       !< coefficient for inertial impaction
    REAL(wp), INTENT(in) ::  c_interception    !< coefficient for interception
    REAL(wp), INTENT(in) ::  c_turb_impaction  !< coefficient for turbulent impaction
    REAL(wp), INTENT(in) ::  kvis_a       !< kinematic viscosity of air (m2/s)
    REAL(wp), INTENT(in) ::  mag_u        !< wind velocity (m/s)
    REAL(wp), INTENT(in) ::  par_l        !< obstacle characteristic dimension in P10
    REAL(wp), INTENT(in) ::  diameter       !< particle diameter
    REAL(wp), INTENT(in) ::  schmidt_num  !< particle Schmidt number
    REAL(wp), INTENT(in) ::  ustar        !< friction velocity (m/s)
    REAL(wp), INTENT(in) ::  vc           !< terminal velocity (m/s)

    REAL(wp), INTENT(inout)  ::  depo     !< deposition efficiency (m/s)

    IF ( par_l > 0.0_wp )  THEN
!
!--    Initialise
       tau_plus = 0.0_wp
       v_bd     = 0.0_wp
       v_im     = 0.0_wp
       v_in     = 0.0_wp
       v_it     = 0.0_wp
!
!--    Stokes number for vegetated surfaces (Seinfeld & Pandis (2006): Eq.19.24)
       stokes_num = vc * ustar / ( g * par_l )
!
!--    Non-dimensional relexation time of the particle on top of canopy
       tau_plus = vc * ustar**2 / ( kvis_a * g )
!
!--    Brownian diffusion
       v_bd = mag_u * c_brownian_diff * schmidt_num**( -0.66666666_wp ) *                          &
              ( mag_u * par_l / kvis_a )**( -0.5_wp )
!
!--    Interception
       v_in = mag_u * c_interception * diameter / par_l *                                          &
              ( 2.0_wp + LOG( 2.0_wp * par_l / diameter ) )
!
!--    Impaction: Petroff (2009) Eq. 18
       v_im = mag_u * c_impaction * ( stokes_num / ( stokes_num + beta_im ) )**2
!
!--    Turbulent impaction
       IF ( tau_plus < 20.0_wp )  THEN
          v_it = 2.5E-3_wp * c_turb_impaction * tau_plus**2
       ELSE
          v_it = c_turb_impaction
       ENDIF

       depo = ( v_bd + v_in + v_im + v_it + vc )

    ELSE
       depo = 0.0_wp
    ENDIF

 END SUBROUTINE depo_vel_P10

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate the dry deposition on horizontal and vertical surfaces. Implement
!> as a surface flux.
!> @todo aerodynamic resistance ignored for now (not important for 
!        high-resolution simulations)
!------------------------------------------------------------------------------!
 SUBROUTINE depo_surf( i, j, surf, vc, schmidt_num, kvis, mag_u, norm, match_array )

    USE arrays_3d,                                                                                 &
        ONLY: rho_air_zw

    USE surface_mod,                                                                               &
        ONLY:  ind_pav_green, ind_veg_wall, ind_wat_win, surf_type

    IMPLICIT NONE

    INTEGER(iwp) ::  ib      !< loop index
    INTEGER(iwp) ::  ic      !< loop index
    INTEGER(iwp) ::  icc     !< additional loop index
    INTEGER(iwp) ::  k       !< loop index
    INTEGER(iwp) ::  m       !< loop index
    INTEGER(iwp) ::  surf_e  !< End index of surface elements at (j,i)-gridpoint
    INTEGER(iwp) ::  surf_s  !< Start index of surface elements at (j,i)-gridpoint

    INTEGER(iwp), INTENT(in) ::  i  !< loop index
    INTEGER(iwp), INTENT(in) ::  j  !< loop index

    LOGICAL, INTENT(in) ::  norm   !< to normalise or not

    REAL(wp) ::  alpha             !< parameter, Table 3 in Z01
    REAL(wp) ::  beta_im           !< parameter for turbulent impaction
    REAL(wp) ::  c_brownian_diff   !< coefficient for Brownian diffusion
    REAL(wp) ::  c_impaction       !< coefficient for inertial impaction
    REAL(wp) ::  c_interception    !< coefficient for interception
    REAL(wp) ::  c_turb_impaction  !< coefficient for turbulent impaction
    REAL(wp) ::  gamma             !< parameter, Table 3 in Z01
    REAL(wp) ::  norm_fac          !< normalisation factor (usually air density)
    REAL(wp) ::  par_a             !< parameter A for the characteristic radius of collectors,
                                   !< Table 3 in Z01
    REAL(wp) ::  par_l             !< obstacle characteristic dimension in P10
    REAL(wp) ::  rs                !< the overall quasi-laminar resistance for particles
    REAL(wp) ::  tau_plus          !< dimensionless particle relaxation time
    REAL(wp) ::  v_bd              !< deposition velocity due to Brownian diffusion
    REAL(wp) ::  v_im              !< deposition velocity due to impaction
    REAL(wp) ::  v_in              !< deposition velocity due to interception
    REAL(wp) ::  v_it              !< deposition velocity due to turbulent impaction

    REAL(wp), DIMENSION(nbins_aerosol) ::  depo      !< deposition efficiency
    REAL(wp), DIMENSION(nbins_aerosol) ::  depo_sum  !< sum of deposition efficiencies

    REAL(wp), DIMENSION(:), INTENT(in) ::  kvis   !< kinematic viscosity of air (m2/s)
    REAL(wp), DIMENSION(:), INTENT(in) ::  mag_u  !< wind velocity (m/s)

    REAL(wp), DIMENSION(:,:), INTENT(in) ::  schmidt_num   !< particle Schmidt number
    REAL(wp), DIMENSION(:,:), INTENT(in) ::  vc            !< terminal velocity (m/s)

    TYPE(match_surface), INTENT(in), OPTIONAL ::  match_array  !< match the deposition module and
                                                               !< LSM/USM surfaces
    TYPE(surf_type), INTENT(inout) :: surf                     !< respective surface type
!
!-- Initialise
    depo     = 0.0_wp
    depo_sum = 0.0_wp
    rs       = 0.0_wp
    surf_s   = surf%start_index(j,i)
    surf_e   = surf%end_index(j,i)
    tau_plus = 0.0_wp
    v_bd     = 0.0_wp
    v_im     = 0.0_wp
    v_in     = 0.0_wp
    v_it     = 0.0_wp
!
!-- Model parameters for the land use category. If LSM or USM is applied, import
!-- characteristics. Otherwise, apply surface type "urban".
    alpha   = alpha_z01(luc_urban)
    gamma   = gamma_z01(luc_urban)
    par_a   = A_z01(luc_urban, season_z01) * 1.0E-3_wp

    par_l            = l_p10(luc_urban) * 0.01_wp
    c_brownian_diff  = c_b_p10(luc_urban)
    c_interception   = c_in_p10(luc_urban)
    c_impaction      = c_im_p10(luc_urban)
    beta_im          = beta_im_p10(luc_urban)
    c_turb_impaction = c_it_p10(luc_urban)


    IF ( PRESENT( match_array ) )  THEN  ! land or urban surface model

       DO  m = surf_s, surf_e

          k = surf%k(m)
          norm_fac = 1.0_wp

          IF ( norm )  norm_fac = rho_air_zw(k)  ! normalise vertical fluxes by air density

          IF ( match_array%match_lupg(m) > 0 )  THEN
             alpha = alpha_z01( match_array%match_lupg(m) )
             gamma = gamma_z01( match_array%match_lupg(m) )
             par_a = A_z01( match_array%match_lupg(m), season_z01 ) * 1.0E-3_wp

             beta_im          = beta_im_p10( match_array%match_lupg(m) )
             c_brownian_diff  = c_b_p10( match_array%match_lupg(m) )
             c_impaction      = c_im_p10( match_array%match_lupg(m) )
             c_interception   = c_in_p10( match_array%match_lupg(m) )
             c_turb_impaction = c_it_p10( match_array%match_lupg(m) )
             par_l            = l_p10( match_array%match_lupg(m) ) * 0.01_wp

             DO  ib = 1, nbins_aerosol
                IF ( aerosol_number(ib)%conc(k,j,i) < ( 2.0_wp * nclim )  .OR.                     &
                     schmidt_num(k+1,ib) < 1.0_wp )  CYCLE

                SELECT CASE ( depo_surf_par_num )

                   CASE ( 1 )
                      CALL depo_vel_Z01( vc(k+1,ib), surf%us(m), schmidt_num(k+1,ib),              &
                                         ra_dry(k,j,i,ib), alpha, gamma, par_a, depo(ib) )
                   CASE ( 2 )
                      CALL depo_vel_P10( vc(k+1,ib), mag_u(k+1), surf%us(m), kvis(k+1),            &
                                         schmidt_num(k+1,ib), ra_dry(k,j,i,ib), par_l,             &
                                         c_brownian_diff, c_interception, c_impaction, beta_im,    &
                                         c_turb_impaction, depo(ib) )
                END SELECT
             ENDDO
             depo_sum = depo_sum + surf%frac(m,ind_pav_green) * depo
          ENDIF

          IF ( match_array%match_luvw(m) > 0 )  THEN
             alpha = alpha_z01( match_array%match_luvw(m) )
             gamma = gamma_z01( match_array%match_luvw(m) )
             par_a = A_z01( match_array%match_luvw(m), season_z01 ) * 1.0E-3_wp

             beta_im          = beta_im_p10( match_array%match_luvw(m) )
             c_brownian_diff  = c_b_p10( match_array%match_luvw(m) )
             c_impaction      = c_im_p10( match_array%match_luvw(m) )
             c_interception   = c_in_p10( match_array%match_luvw(m) )
             c_turb_impaction = c_it_p10( match_array%match_luvw(m) )
             par_l            = l_p10( match_array%match_luvw(m) ) * 0.01_wp

             DO  ib = 1, nbins_aerosol
                IF ( aerosol_number(ib)%conc(k,j,i) < ( 2.0_wp * nclim )  .OR.                     &
                     schmidt_num(k+1,ib) < 1.0_wp )  CYCLE

                SELECT CASE ( depo_surf_par_num )

                   CASE ( 1 )
                      CALL depo_vel_Z01( vc(k+1,ib), surf%us(m), schmidt_num(k+1,ib),              &
                                         ra_dry(k,j,i,ib), alpha, gamma, par_a, depo(ib) )
                   CASE ( 2 )
                      CALL depo_vel_P10( vc(k+1,ib), mag_u(k+1), surf%us(m), kvis(k+1),            &
                                         schmidt_num(k+1,ib), ra_dry(k,j,i,ib), par_l,             &
                                         c_brownian_diff, c_interception, c_impaction, beta_im,    &
                                         c_turb_impaction, depo(ib) )
                END SELECT
             ENDDO
             depo_sum = depo_sum + surf%frac(m,ind_veg_wall) * depo
          ENDIF

          IF ( match_array%match_luww(m) > 0 )  THEN
             alpha = alpha_z01( match_array%match_luww(m) )
             gamma = gamma_z01( match_array%match_luww(m) )
             par_a = A_z01( match_array%match_luww(m), season_z01 ) * 1.0E-3_wp

             beta_im          = beta_im_p10( match_array%match_luww(m) )
             c_brownian_diff  = c_b_p10( match_array%match_luww(m) )
             c_impaction      = c_im_p10( match_array%match_luww(m) )
             c_interception   = c_in_p10( match_array%match_luww(m) )
             c_turb_impaction = c_it_p10( match_array%match_luww(m) )
             par_l            = l_p10( match_array%match_luww(m) ) * 0.01_wp

             DO  ib = 1, nbins_aerosol
                IF ( aerosol_number(ib)%conc(k,j,i) < ( 2.0_wp * nclim )  .OR.                     &
                     schmidt_num(k+1,ib) < 1.0_wp )  CYCLE

                SELECT CASE ( depo_surf_par_num )

                   CASE ( 1 )
                      CALL depo_vel_Z01( vc(k+1,ib), surf%us(m), schmidt_num(k+1,ib),              &
                                         ra_dry(k,j,i,ib), alpha, gamma, par_a, depo(ib) )
                   CASE ( 2 )
                      CALL depo_vel_P10( vc(k+1,ib), mag_u(k+1), surf%us(m), kvis(k+1),            &
                                         schmidt_num(k+1,ib), ra_dry(k,j,i,ib), par_l,             &
                                         c_brownian_diff, c_interception, c_impaction, beta_im,    &
                                         c_turb_impaction, depo(ib) )
                END SELECT
             ENDDO
             depo_sum = depo_sum + surf%frac(m,ind_wat_win) * depo
          ENDIF

          DO  ib = 1, nbins_aerosol
             IF ( aerosol_number(ib)%conc(k,j,i) < ( 2.0_wp * nclim ) )  CYCLE
!
!--          Calculate changes in surface fluxes due to dry deposition
             IF ( include_emission )  THEN
                surf%answs(m,ib) = aerosol_number(ib)%source(j,i) - MAX( 0.0_wp,                   &
                                   depo_sum(ib) * norm_fac * aerosol_number(ib)%conc(k,j,i) )
                DO  ic = 1, ncomponents_mass
                   icc = ( ic - 1 ) * nbins_aerosol + ib
                   surf%amsws(m,icc) = aerosol_mass(icc)%source(j,i) - MAX( 0.0_wp,                &
                                       depo_sum(ib) *  norm_fac * aerosol_mass(icc)%conc(k,j,i) )
                ENDDO  ! ic
             ELSE
                surf%answs(m,ib) = -depo_sum(ib) * norm_fac * aerosol_number(ib)%conc(k,j,i)
                DO  ic = 1, ncomponents_mass
                   icc = ( ic - 1 ) * nbins_aerosol + ib
                   surf%amsws(m,icc) = -depo_sum(ib) *  norm_fac * aerosol_mass(icc)%conc(k,j,i)
                ENDDO  ! ic
             ENDIF
          ENDDO  ! ib

       ENDDO

    ELSE  ! default surfaces

       DO  m = surf_s, surf_e

          k = surf%k(m)
          norm_fac = 1.0_wp

          IF ( norm )  norm_fac = rho_air_zw(k)  ! normalise vertical fluxes by air density

          DO  ib = 1, nbins_aerosol
             IF ( aerosol_number(ib)%conc(k,j,i) < ( 2.0_wp * nclim )  .OR.                        &
                  schmidt_num(k+1,ib) < 1.0_wp )  CYCLE

             SELECT CASE ( depo_surf_par_num )

                CASE ( 1 )
                   CALL depo_vel_Z01( vc(k+1,ib), surf%us(m), schmidt_num(k+1,ib),                 &
                                      ra_dry(k,j,i,ib), alpha, gamma, par_a, depo(ib) )
                CASE ( 2 )
                   CALL depo_vel_P10( vc(k+1,ib), mag_u(k+1), surf%us(m), kvis(k+1),               &
                                      schmidt_num(k+1,ib), ra_dry(k,j,i,ib), par_l,                &
                                      c_brownian_diff, c_interception, c_impaction, beta_im,       &
                                      c_turb_impaction, depo(ib) )
             END SELECT
!
!--          Calculate changes in surface fluxes due to dry deposition
             IF ( include_emission )  THEN
                surf%answs(m,ib) = aerosol_number(ib)%source(j,i) - MAX( 0.0_wp,                   &
                                   depo(ib) * norm_fac * aerosol_number(ib)%conc(k,j,i) )
                DO  ic = 1, ncomponents_mass
                   icc = ( ic - 1 ) * nbins_aerosol + ib
                   surf%amsws(m,icc) = aerosol_mass(icc)%source(j,i) - MAX( 0.0_wp,                &
                                       depo(ib) *  norm_fac * aerosol_mass(icc)%conc(k,j,i) )
                ENDDO  ! ic
             ELSE
                surf%answs(m,ib) = -depo(ib) * norm_fac * aerosol_number(ib)%conc(k,j,i)
                DO  ic = 1, ncomponents_mass
                   icc = ( ic - 1 ) * nbins_aerosol + ib
                   surf%amsws(m,icc) = -depo(ib) *  norm_fac * aerosol_mass(icc)%conc(k,j,i)
                ENDDO  ! ic
             ENDIF
          ENDDO  ! ib
       ENDDO

    ENDIF

 END SUBROUTINE depo_surf

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculates particle loss and change in size distribution due to (Brownian)
!> coagulation. Only for particles with dwet < 30 micrometres. 
!
!> Method:
!> Semi-implicit, non-iterative method: (Jacobson, 1994)
!> Volume concentrations of the smaller colliding particles added to the bin of
!> the larger colliding particles. Start from first bin and use the updated 
!> number and volume for calculation of following bins. NB! Our bin numbering
!> does not follow particle size in subrange 2.
!
!> Schematic for bin numbers in different subranges:
!>             1                            2
!>    +-------------------------------------------+
!>  a | 1 | 2 | 3 || 4 | 5 | 6 | 7 |  8 |  9 | 10||
!>  b |           ||11 |12 |13 |14 | 15 | 16 | 17||
!>    +-------------------------------------------+
!
!> Exact coagulation coefficients for each pressure level are scaled according
!> to current particle wet size (linear scaling).
!> Bins are organized in terms of the dry size of the condensation nucleus, 
!> while coagulation kernell is calculated with the actual hydrometeor
!> size.
!
!> Called from salsa_driver
!> fxm: Process selection should be made smarter - now just lots of IFs inside
!>      loops
!
!> Coded by:
!> Hannele Korhonen (FMI) 2005
!> Harri Kokkola (FMI) 2006
!> Tommi Bergman (FMI) 2012
!> Matti Niskanen(FMI) 2012
!> Anton Laakso  (FMI) 2013
!> Juha Tonttila (FMI) 2014
!------------------------------------------------------------------------------!
 SUBROUTINE coagulation( paero, ptstep, ptemp, ppres )

    IMPLICIT NONE

    INTEGER(iwp) ::  index_2a !< corresponding bin in subrange 2a
    INTEGER(iwp) ::  index_2b !< corresponding bin in subrange 2b
    INTEGER(iwp) ::  ib       !< loop index
    INTEGER(iwp) ::  ll       !< loop index
    INTEGER(iwp) ::  mm       !< loop index
    INTEGER(iwp) ::  nn       !< loop index

    REAL(wp) ::  pressi          !< pressure
    REAL(wp) ::  temppi          !< temperature
    REAL(wp) ::  zdpart_mm       !< diameter of particle (m)
    REAL(wp) ::  zdpart_nn       !< diameter of particle (m)
    REAL(wp) ::  zminusterm      !< coagulation loss in a bin (1/s)

    REAL(wp), INTENT(in) ::  ppres  !< ambient pressure (Pa)
    REAL(wp), INTENT(in) ::  ptemp  !< ambient temperature (K)
    REAL(wp), INTENT(in) ::  ptstep !< time step (s)

    REAL(wp), DIMENSION(nbins_aerosol) ::  zmpart     !< approximate mass of particles (kg)
    REAL(wp), DIMENSION(maxspec+1)     ::  zplusterm  !< coagulation gain in a bin (for each
                                                      !< chemical compound)
    REAL(wp), DIMENSION(nbins_aerosol,nbins_aerosol) ::  zcc  !< updated coagulation coeff. (m3/s)

    TYPE(t_section), DIMENSION(nbins_aerosol), INTENT(inout) ::  paero  !< Aerosol properties

    zdpart_mm = 0.0_wp
    zdpart_nn = 0.0_wp
!
!-- 1) Coagulation to coarse mode calculated in a simplified way:
!--    CoagSink ~ Dp in continuum subrange --> 'effective' number conc. of coarse particles

!-- 2) Updating coagulation coefficients 
!
!-- Aerosol mass (kg). Density of 1500 kg/m3 assumed
    zmpart(1:end_subrange_2b) = api6 * ( MIN( paero(1:end_subrange_2b)%dwet, 30.0E-6_wp )**3 )     &
                                * 1500.0_wp
    temppi = ptemp
    pressi = ppres
    zcc    = 0.0_wp
!
!-- Aero-aero coagulation
    DO  mm = 1, end_subrange_2b   ! smaller colliding particle
       IF ( paero(mm)%numc < ( 2.0_wp * nclim ) )  CYCLE
       DO  nn = mm, end_subrange_2b   ! larger colliding particle
          IF ( paero(nn)%numc < ( 2.0_wp * nclim ) )  CYCLE

          zdpart_mm = MIN( paero(mm)%dwet, 30.0E-6_wp )     ! Limit to 30 um
          zdpart_nn = MIN( paero(nn)%dwet, 30.0E-6_wp )     ! Limit to 30 um
!
!--       Coagulation coefficient of particles (m3/s)
          zcc(mm,nn) = coagc( zdpart_mm, zdpart_nn, zmpart(mm), zmpart(nn), temppi, pressi )
          zcc(nn,mm) = zcc(mm,nn)
       ENDDO
    ENDDO

!
!-- 3) New particle and volume concentrations after coagulation:
!--    Calculated according to Jacobson (2005) eq. 15.9
!
!-- Aerosols in subrange 1a:
    DO  ib = start_subrange_1a, end_subrange_1a
       IF ( paero(ib)%numc < ( 2.0_wp * nclim ) )  CYCLE
       zminusterm   = 0.0_wp
       zplusterm(:) = 0.0_wp
!
!--    Particles lost by coagulation with larger aerosols
       DO  ll = ib+1, end_subrange_2b
          zminusterm = zminusterm + zcc(ib,ll) * paero(ll)%numc
       ENDDO
!
!--    Coagulation gain in a bin: change in volume conc. (cm3/cm3):
       DO ll = start_subrange_1a, ib - 1
          zplusterm(1:2) = zplusterm(1:2) + zcc(ll,ib) * paero(ll)%volc(1:2)
          zplusterm(6:7) = zplusterm(6:7) + zcc(ll,ib) * paero(ll)%volc(6:7)
          zplusterm(8)   = zplusterm(8)   + zcc(ll,ib) * paero(ll)%volc(8)
       ENDDO
!
!--    Volume and number concentrations after coagulation update [fxm]
       paero(ib)%volc(1:2) = ( paero(ib)%volc(1:2) + ptstep * zplusterm(1:2) * paero(ib)%numc ) /  &
                            ( 1.0_wp + ptstep * zminusterm )
       paero(ib)%volc(6:8) = ( paero(ib)%volc(6:8) + ptstep * zplusterm(6:8) * paero(ib)%numc ) /  &
                            ( 1.0_wp + ptstep * zminusterm )
       paero(ib)%numc = paero(ib)%numc / ( 1.0_wp + ptstep * zminusterm + 0.5_wp * ptstep *        &
                        zcc(ib,ib) * paero(ib)%numc )
    ENDDO
!
!-- Aerosols in subrange 2a:
    DO  ib = start_subrange_2a, end_subrange_2a
       IF ( paero(ib)%numc < ( 2.0_wp * nclim ) )  CYCLE
       zminusterm   = 0.0_wp
       zplusterm(:) = 0.0_wp
!
!--    Find corresponding size bin in subrange 2b
       index_2b = ib - start_subrange_2a + start_subrange_2b
!
!--    Particles lost by larger particles in 2a
       DO  ll = ib+1, end_subrange_2a
          zminusterm = zminusterm + zcc(ib,ll) * paero(ll)%numc
       ENDDO
!
!--    Particles lost by larger particles in 2b
       IF ( .NOT. no_insoluble )  THEN
          DO  ll = index_2b+1, end_subrange_2b
             zminusterm = zminusterm + zcc(ib,ll) * paero(ll)%numc
          ENDDO
       ENDIF
!
!--    Particle volume gained from smaller particles in subranges 1, 2a and 2b
       DO  ll = start_subrange_1a, ib-1
          zplusterm(1:2) = zplusterm(1:2) + zcc(ll,ib) * paero(ll)%volc(1:2)
          zplusterm(6:8) = zplusterm(6:8) + zcc(ll,ib) * paero(ll)%volc(6:8)
       ENDDO
!
!--    Particle volume gained from smaller particles in 2a
!--    (Note, for components not included in the previous loop!)
       DO  ll = start_subrange_2a, ib-1
          zplusterm(3:5) = zplusterm(3:5) + zcc(ll,ib)*paero(ll)%volc(3:5)
       ENDDO
!
!--    Particle volume gained from smaller (and equal) particles in 2b
       IF ( .NOT. no_insoluble )  THEN
          DO  ll = start_subrange_2b, index_2b
             zplusterm(1:8) = zplusterm(1:8) + zcc(ll,ib) * paero(ll)%volc(1:8)
          ENDDO
       ENDIF
!
!--    Volume and number concentrations after coagulation update [fxm]
       paero(ib)%volc(1:8) = ( paero(ib)%volc(1:8) + ptstep * zplusterm(1:8) * paero(ib)%numc ) /  &
                            ( 1.0_wp + ptstep * zminusterm )
       paero(ib)%numc = paero(ib)%numc / ( 1.0_wp + ptstep * zminusterm + 0.5_wp * ptstep *        &
                        zcc(ib,ib) * paero(ib)%numc )
    ENDDO
!
!-- Aerosols in subrange 2b:
    IF ( .NOT. no_insoluble )  THEN
       DO  ib = start_subrange_2b, end_subrange_2b
          IF ( paero(ib)%numc < ( 2.0_wp * nclim ) )  CYCLE
          zminusterm   = 0.0_wp
          zplusterm(:) = 0.0_wp
!
!--       Find corresponding size bin in subsubrange 2a
          index_2a = ib - start_subrange_2b + start_subrange_2a
!
!--       Particles lost to larger particles in subranges 2b
          DO  ll = ib + 1, end_subrange_2b
             zminusterm = zminusterm + zcc(ib,ll) * paero(ll)%numc
          ENDDO
!
!--       Particles lost to larger and equal particles in 2a
          DO  ll = index_2a, end_subrange_2a
             zminusterm = zminusterm + zcc(ib,ll) * paero(ll)%numc
          ENDDO
!
!--       Particle volume gained from smaller particles in subranges 1 & 2a
          DO  ll = start_subrange_1a, index_2a - 1
             zplusterm(1:8) = zplusterm(1:8) + zcc(ll,ib) * paero(ll)%volc(1:8)
          ENDDO
!
!--       Particle volume gained from smaller particles in 2b
          DO  ll = start_subrange_2b, ib - 1
             zplusterm(1:8) = zplusterm(1:8) + zcc(ll,ib) * paero(ll)%volc(1:8)
          ENDDO
!
!--       Volume and number concentrations after coagulation update [fxm]
          paero(ib)%volc(1:8) = ( paero(ib)%volc(1:8) + ptstep * zplusterm(1:8) * paero(ib)%numc ) &
                                / ( 1.0_wp + ptstep * zminusterm )
          paero(ib)%numc = paero(ib)%numc / ( 1.0_wp + ptstep * zminusterm + 0.5_wp * ptstep *     &
                           zcc(ib,ib) * paero(ib)%numc )
       ENDDO
    ENDIF

 END SUBROUTINE coagulation

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculation of coagulation coefficients. Extended version of the function
!> originally found in mo_salsa_init.
!
!> J. Tonttila, FMI, 05/2014
!------------------------------------------------------------------------------!
 REAL(wp) FUNCTION coagc( diam1, diam2, mass1, mass2, temp, pres )

    IMPLICIT NONE

    REAL(wp) ::  fmdist  !< distance of flux matching (m)
    REAL(wp) ::  knud_p  !< particle Knudsen number
    REAL(wp) ::  mdiam   !< mean diameter of colliding particles (m)
    REAL(wp) ::  mfp     !< mean free path of air molecules (m)
    REAL(wp) ::  visc    !< viscosity of air (kg/(m s))

    REAL(wp), INTENT(in) ::  diam1  !< diameter of colliding particle 1 (m)
    REAL(wp), INTENT(in) ::  diam2  !< diameter of colliding particle 2 (m)
    REAL(wp), INTENT(in) ::  mass1  !< mass of colliding particle 1 (kg)
    REAL(wp), INTENT(in) ::  mass2  !< mass of colliding particle 2 (kg)
    REAL(wp), INTENT(in) ::  pres   !< ambient pressure (Pa?) [fxm]
    REAL(wp), INTENT(in) ::  temp   !< ambient temperature (K)

    REAL(wp), DIMENSION (2) ::  beta    !< Cunningham correction factor
    REAL(wp), DIMENSION (2) ::  dfpart  !< particle diffusion coefficient (m2/s)
    REAL(wp), DIMENSION (2) ::  diam    !< diameters of particles (m)
    REAL(wp), DIMENSION (2) ::  flux    !< flux in continuum and free molec. regime (m/s)
    REAL(wp), DIMENSION (2) ::  knud    !< particle Knudsen number
    REAL(wp), DIMENSION (2) ::  mpart   !< masses of particles (kg)
    REAL(wp), DIMENSION (2) ::  mtvel   !< particle mean thermal velocity (m/s)
    REAL(wp), DIMENSION (2) ::  omega   !< particle mean free path
    REAL(wp), DIMENSION (2) ::  tva     !< temporary variable (m)
!
!-- Initialisation
    coagc   = 0.0_wp
!
!-- 1) Initializing particle and ambient air variables
    diam  = (/ diam1, diam2 /) !< particle diameters (m)
    mpart = (/ mass1, mass2 /) !< particle masses (kg)
!
!-- Viscosity of air (kg/(m s))
    visc = ( 7.44523E-3_wp * temp ** 1.5_wp ) / ( 5093.0_wp * ( temp + 110.4_wp ) )
!
!-- Mean free path of air (m)
    mfp = ( 1.656E-10_wp * temp + 1.828E-8_wp ) * ( p_0 + 1325.0_wp ) / pres
!
!-- 2) Slip correction factor for small particles 
    knud = 2.0_wp * EXP( LOG(mfp) - LOG(diam) )! Knudsen number for air (15.23)
!
!-- Cunningham correction factor (Allen and Raabe, Aerosol Sci. Tech. 4, 269)
    beta = 1.0_wp + knud * ( 1.142_wp + 0.558_wp * EXP( -0.999_wp / knud ) )
!
!-- 3) Particle properties
!-- Diffusion coefficient (m2/s) (Jacobson (2005) eq. 15.29)
    dfpart = beta * abo * temp / ( 3.0_wp * pi * visc * diam )
!
!-- Mean thermal velocity (m/s) (Jacobson (2005) eq. 15.32)
    mtvel = SQRT( ( 8.0_wp * abo * temp ) / ( pi * mpart ) )
!
!-- Particle mean free path (m) (Jacobson (2005) eq. 15.34 )
    omega = 8.0_wp * dfpart / ( pi * mtvel )
!
!-- Mean diameter (m)
    mdiam = 0.5_wp * ( diam(1) + diam(2) )
!
!-- 4) Calculation of fluxes (Brownian collision kernels) and flux matching
!-- following Jacobson (2005):
!
!-- Flux in continuum regime (m3/s) (eq. 15.28)
    flux(1) = 4.0_wp * pi * mdiam * ( dfpart(1) + dfpart(2) )
!
!-- Flux in free molec. regime (m3/s) (eq. 15.31)
    flux(2) = pi * SQRT( ( mtvel(1)**2 ) + ( mtvel(2)**2 ) ) * ( mdiam**2 )
!
!-- temporary variables (m) to calculate flux matching distance (m)
    tva(1) = ( ( mdiam + omega(1) )**3 - ( mdiam**2 + omega(1)**2 ) * SQRT( ( mdiam**2 +           &
               omega(1)**2 ) ) ) / ( 3.0_wp * mdiam * omega(1) ) - mdiam
    tva(2) = ( ( mdiam + omega(2) )**3 - ( mdiam**2 + omega(2)**2 ) * SQRT( ( mdiam**2 +           &
               omega(2)**2 ) ) ) / ( 3.0_wp * mdiam * omega(2) ) - mdiam
!
!-- Flux matching distance (m): the mean distance from the centre of a sphere reached by particles
!-- that leave sphere's surface and travel a distance of particle mean free path (eq. 15.34)
    fmdist = SQRT( tva(1)**2 + tva(2)**2 )
!
!-- 5) Coagulation coefficient = coalescence efficiency * collision kernel (m3/s) (eq. 15.33).
!--    Here assumed coalescence efficiency 1!!
    coagc = flux(1) / ( mdiam / ( mdiam + fmdist) + flux(1) / flux(2) )
!
!-- Corrected collision kernel (Karl et al., 2016 (ACP)): Include van der Waals and viscous forces
    IF ( van_der_waals_coagc )  THEN
       knud_p = SQRT( omega(1)**2 + omega(2)**2 ) / mdiam
       IF ( knud_p >= 0.1_wp  .AND.  knud_p <= 10.0_wp )  THEN
          coagc = coagc * ( 2.0_wp + 0.4_wp * LOG( knud_p ) )
       ELSE
          coagc = coagc * 3.0_wp
       ENDIF
    ENDIF

 END FUNCTION coagc

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculates the change in particle volume and gas phase 
!> concentrations due to nucleation, condensation and dissolutional growth.
!
!> Sulphuric acid and organic vapour: only condensation and no evaporation.
!
!> New gas and aerosol phase concentrations calculated according to Jacobson 
!> (1997): Numerical techniques to solve condensational and dissolutional growth
!> equations when growth is coupled to reversible reactions, Aerosol Sci. Tech.,
!> 27, pp 491-498.
!
!> Following parameterization has been used:
!> Molecular diffusion coefficient of condensing vapour (m2/s) 
!> (Reid et al. (1987): Properties of gases and liquids, McGraw-Hill, New York.)
!> D = {1.d-7*sqrt(1/M_air + 1/M_gas)*T^1.75} / &
!      {p_atm/p_stand * (d_air^(1/3) + d_gas^(1/3))^2 }
!> M_air = 28.965 : molar mass of air (g/mol)
!> d_air = 19.70  : diffusion volume of air
!> M_h2so4 = 98.08 : molar mass of h2so4 (g/mol)
!> d_h2so4 = 51.96  : diffusion volume of h2so4 
!
!> Called from main aerosol model
!> For equations, see Jacobson, Fundamentals of Atmospheric Modeling, 2nd Edition (2005)
!
!> Coded by:
!> Hannele Korhonen (FMI) 2005
!> Harri Kokkola (FMI) 2006
!> Juha Tonttila (FMI) 2014
!> Rewritten to PALM by Mona Kurppa (UHel) 2017
!------------------------------------------------------------------------------!
 SUBROUTINE condensation( paero, pc_sa, pc_ocnv, pcocsv, pchno3, pc_nh3, pcw, pcs, ptemp, ppres,   &
                          ptstep, prtcl )

    IMPLICIT NONE

    INTEGER(iwp) ::  ss      !< start index
    INTEGER(iwp) ::  ee      !< end index

    REAL(wp) ::  zcs_ocnv    !< condensation sink of nonvolatile organics (1/s)
    REAL(wp) ::  zcs_ocsv    !< condensation sink of semivolatile organics (1/s)
    REAL(wp) ::  zcs_su      !< condensation sink of sulfate (1/s)
    REAL(wp) ::  zcs_tot     !< total condensation sink (1/s) (gases)
    REAL(wp) ::  zcvap_new1  !< vapour concentration after time step (#/m3): sulphuric acid
    REAL(wp) ::  zcvap_new2  !< nonvolatile organics
    REAL(wp) ::  zcvap_new3  !< semivolatile organics
    REAL(wp) ::  zdfvap      !< air diffusion coefficient (m2/s)
    REAL(wp) ::  zdvap1      !< change in vapour concentration (#/m3): sulphuric acid
    REAL(wp) ::  zdvap2      !< nonvolatile organics
    REAL(wp) ::  zdvap3      !< semivolatile organics
    REAL(wp) ::  zmfp        !< mean free path of condensing vapour (m)
    REAL(wp) ::  zrh         !< Relative humidity [0-1]
    REAL(wp) ::  zvisc       !< viscosity of air (kg/(m s))
    REAL(wp) ::  zn_vs_c     !< ratio of nucleation of all mass transfer in the smallest bin
    REAL(wp) ::  zxocnv      !< ratio of organic vapour in 3nm particles
    REAL(wp) ::  zxsa        !< Ratio in 3nm particles: sulphuric acid

    REAL(wp), INTENT(in) ::  ppres   !< ambient pressure (Pa)
    REAL(wp), INTENT(in) ::  pcs     !< Water vapour saturation concentration (kg/m3)
    REAL(wp), INTENT(in) ::  ptemp   !< ambient temperature (K)
    REAL(wp), INTENT(in) ::  ptstep  !< timestep (s)

    REAL(wp), INTENT(inout) ::  pchno3   !< Gas concentrations (#/m3): nitric acid HNO3
    REAL(wp), INTENT(inout) ::  pc_nh3   !< ammonia NH3
    REAL(wp), INTENT(inout) ::  pc_ocnv  !< non-volatile organics
    REAL(wp), INTENT(inout) ::  pcocsv   !< semi-volatile organics
    REAL(wp), INTENT(inout) ::  pc_sa    !< sulphuric acid H2SO4
    REAL(wp), INTENT(inout) ::  pcw      !< Water vapor concentration (kg/m3)

    REAL(wp), DIMENSION(nbins_aerosol)       ::  zbeta          !< transitional correction factor
    REAL(wp), DIMENSION(nbins_aerosol)       ::  zcolrate       !< collision rate (1/s)
    REAL(wp), DIMENSION(nbins_aerosol)       ::  zcolrate_ocnv  !< collision rate of OCNV (1/s)
    REAL(wp), DIMENSION(start_subrange_1a+1) ::  zdfpart        !< particle diffusion coef. (m2/s)
    REAL(wp), DIMENSION(nbins_aerosol)       ::  zdvoloc        !< change of organics volume
    REAL(wp), DIMENSION(nbins_aerosol)       ::  zdvolsa        !< change of sulphate volume
    REAL(wp), DIMENSION(2)                   ::  zj3n3          !< Formation massrate of molecules
                                                                !< in nucleation, (molec/m3s),
                                                                !< 1: H2SO4 and 2: organic vapor
    REAL(wp), DIMENSION(nbins_aerosol)       ::  zknud          !< particle Knudsen number

    TYPE(component_index), INTENT(in) :: prtcl  !< Keeps track which substances are used

    TYPE(t_section), DIMENSION(nbins_aerosol), INTENT(inout) ::  paero  !< Aerosol properties

    zj3n3  = 0.0_wp
    zrh    = pcw / pcs
    zxocnv = 0.0_wp
    zxsa   = 0.0_wp
!
!-- Nucleation
    IF ( nsnucl > 0 )  THEN
       CALL nucleation( paero, ptemp, zrh, ppres, pc_sa, pc_ocnv, pc_nh3, ptstep, zj3n3, zxsa,     &
                        zxocnv )
    ENDIF
!
!-- Condensation on pre-existing particles
    IF ( lscndgas )  THEN
!
!--    Initialise:
       zdvolsa = 0.0_wp
       zdvoloc = 0.0_wp
       zcolrate = 0.0_wp
!
!--    1) Properties of air and condensing gases:
!--    Viscosity of air (kg/(m s)) (Eq. 4.54 in Jabonson (2005))
       zvisc = ( 7.44523E-3_wp * ptemp ** 1.5_wp ) / ( 5093.0_wp * ( ptemp + 110.4_wp ) )
!
!--    Diffusion coefficient of air (m2/s)
       zdfvap = 5.1111E-10_wp * ptemp ** 1.75_wp * ( p_0 + 1325.0_wp ) / ppres
!
!--    Mean free path (m): same for H2SO4 and organic compounds
       zmfp = 3.0_wp * zdfvap * SQRT( pi * amh2so4 / ( 8.0_wp * argas * ptemp ) )
!
!--    2) Transition regime correction factor zbeta for particles (Fuchs and Sutugin (1971)):
!--       Size of condensing molecule considered only for nucleation mode (3 - 20 nm).
!
!--    Particle Knudsen number: condensation of gases on aerosols
       ss = start_subrange_1a
       ee = start_subrange_1a+1
       zknud(ss:ee) = 2.0_wp * zmfp / ( paero(ss:ee)%dwet + d_sa )
       ss = start_subrange_1a+2
       ee = end_subrange_2b
       zknud(ss:ee) = 2.0_wp * zmfp / paero(ss:ee)%dwet
!
!--    Transitional correction factor: aerosol + gas (the semi-empirical Fuchs- Sutugin 
!--    interpolation function (Fuchs and Sutugin, 1971))
       zbeta = ( zknud + 1.0_wp ) / ( 0.377_wp * zknud + 1.0_wp + 4.0_wp / ( 3.0_wp * massacc ) *  &
               ( zknud + zknud ** 2 ) )
!
!--    3) Collision rate of molecules to particles
!--       Particle diffusion coefficient considered only for nucleation mode (3 - 20 nm)
!
!--    Particle diffusion coefficient (m2/s) (e.g. Eq. 15.29 in Jacobson (2005))
       zdfpart = abo * ptemp * zbeta(start_subrange_1a:start_subrange_1a+1) / ( 3.0_wp * pi * zvisc&
                 * paero(start_subrange_1a:start_subrange_1a+1)%dwet)
!
!--    Collision rate (mass-transfer coefficient): gases on aerosols (1/s) (Eq. 16.64 in 
!--    Jacobson (2005))
       ss = start_subrange_1a
       ee = start_subrange_1a+1
       zcolrate(ss:ee) = MERGE( 2.0_wp * pi * ( paero(ss:ee)%dwet + d_sa ) * ( zdfvap + zdfpart ) *&
                               zbeta(ss:ee) * paero(ss:ee)%numc, 0.0_wp, paero(ss:ee)%numc > nclim )
       ss = start_subrange_1a+2
       ee = end_subrange_2b
       zcolrate(ss:ee) = MERGE( 2.0_wp * pi * paero(ss:ee)%dwet * zdfvap * zbeta(ss:ee) *          &
                                paero(ss:ee)%numc, 0.0_wp, paero(ss:ee)%numc > nclim )
!
!-- 4) Condensation sink (1/s)
       zcs_tot = SUM( zcolrate )   ! total sink
!
!--    5) Changes in gas-phase concentrations and particle volume
!
!--    5.1) Organic vapours 
!
!--    5.1.1) Non-volatile organic compound: condenses onto all bins
       IF ( pc_ocnv > 1.0E+10_wp  .AND.  zcs_tot > 1.0E-30_wp  .AND. index_oc > 0 )  &
       THEN
!--       Ratio of nucleation vs. condensation rates in the smallest bin
          zn_vs_c = 0.0_wp
          IF ( zj3n3(2) > 1.0_wp )  THEN
             zn_vs_c = ( zj3n3(2) ) / ( zj3n3(2) + pc_ocnv * zcolrate(start_subrange_1a) )
          ENDIF
!
!--       Collision rate in the smallest bin, including nucleation and condensation (see
!--       Jacobson (2005), eq. (16.73) )
          zcolrate_ocnv = zcolrate
          zcolrate_ocnv(start_subrange_1a) = zcolrate_ocnv(start_subrange_1a) + zj3n3(2) / pc_ocnv
!
!--       Total sink for organic vapor
          zcs_ocnv = zcs_tot + zj3n3(2) / pc_ocnv
!
!--       New gas phase concentration (#/m3)
          zcvap_new2 = pc_ocnv / ( 1.0_wp + ptstep * zcs_ocnv )
!
!--       Change in gas concentration (#/m3)
          zdvap2 = pc_ocnv - zcvap_new2
!
!--       Updated vapour concentration (#/m3)
          pc_ocnv = zcvap_new2
!
!--       Volume change of particles (m3(OC)/m3(air))
          zdvoloc = zcolrate_ocnv(start_subrange_1a:end_subrange_2b) / zcs_ocnv * amvoc * zdvap2
!
!--       Change of volume due to condensation in 1a-2b
          paero(start_subrange_1a:end_subrange_2b)%volc(2) =                                       &
                                          paero(start_subrange_1a:end_subrange_2b)%volc(2) + zdvoloc
!
!--       Change of number concentration in the smallest bin caused by nucleation (Jacobson (2005),
!--       eq. (16.75)). If zxocnv = 0, then the chosen nucleation mechanism doesn't take into
!--       account the non-volatile organic vapors and thus the paero doesn't have to be updated.
          IF ( zxocnv > 0.0_wp )  THEN
             paero(start_subrange_1a)%numc = paero(start_subrange_1a)%numc + zn_vs_c *             &
                                             zdvoloc(start_subrange_1a) / amvoc / ( n3 * zxocnv )
          ENDIF
       ENDIF
!
!--    5.1.2) Semivolatile organic compound: all bins except subrange 1
       zcs_ocsv = SUM( zcolrate(start_subrange_2a:end_subrange_2b) ) !< sink for semi-volatile org.
       IF ( pcocsv > 1.0E+10_wp  .AND.  zcs_ocsv > 1.0E-30  .AND. is_used( prtcl,'OC') )  THEN
!
!--       New gas phase concentration (#/m3)
          zcvap_new3 = pcocsv / ( 1.0_wp + ptstep * zcs_ocsv )
!
!--       Change in gas concentration (#/m3)
          zdvap3 = pcocsv - zcvap_new3 
!
!--       Updated gas concentration (#/m3)
          pcocsv = zcvap_new3
!
!--       Volume change of particles (m3(OC)/m3(air))
          ss = start_subrange_2a
          ee = end_subrange_2b
          zdvoloc(ss:ee) = zdvoloc(ss:ee) + zcolrate(ss:ee) / zcs_ocsv * amvoc * zdvap3
!
!--       Change of volume due to condensation in 1a-2b
          paero(start_subrange_1a:end_subrange_2b)%volc(2) =                                       &
                                          paero(start_subrange_1a:end_subrange_2b)%volc(2) + zdvoloc
       ENDIF
!
!--    5.2) Sulphate: condensed on all bins
       IF ( pc_sa > 1.0E+10_wp  .AND.  zcs_tot > 1.0E-30_wp  .AND.  index_so4 > 0 )  THEN
!
!--    Ratio of mass transfer between nucleation and condensation
          zn_vs_c = 0.0_wp
          IF ( zj3n3(1) > 1.0_wp )  THEN
             zn_vs_c = ( zj3n3(1) ) / ( zj3n3(1) + pc_sa * zcolrate(start_subrange_1a) )
          ENDIF
!
!--       Collision rate in the smallest bin, including nucleation and condensation (see
!--       Jacobson (2005), eq. (16.73))
          zcolrate(start_subrange_1a) = zcolrate(start_subrange_1a) + zj3n3(1) / pc_sa
!
!--       Total sink for sulfate (1/s)
          zcs_su = zcs_tot + zj3n3(1) / pc_sa
!
!--       Sulphuric acid:
!--       New gas phase concentration (#/m3)
          zcvap_new1 = pc_sa / ( 1.0_wp + ptstep * zcs_su )
!
!--       Change in gas concentration (#/m3)
          zdvap1 = pc_sa - zcvap_new1
!
!--       Updating vapour concentration (#/m3)
          pc_sa = zcvap_new1
!
!--       Volume change of particles (m3(SO4)/m3(air)) by condensation
          zdvolsa = zcolrate(start_subrange_1a:end_subrange_2b) / zcs_su * amvh2so4 * zdvap1
!
!--       Change of volume concentration of sulphate in aerosol [fxm]
          paero(start_subrange_1a:end_subrange_2b)%volc(1) =                                       &
                                          paero(start_subrange_1a:end_subrange_2b)%volc(1) + zdvolsa
!
!--       Change of number concentration in the smallest bin caused by nucleation
!--       (Jacobson (2005), equation (16.75))
          IF ( zxsa > 0.0_wp )  THEN
             paero(start_subrange_1a)%numc = paero(start_subrange_1a)%numc + zn_vs_c *             &
                                             zdvolsa(start_subrange_1a) / amvh2so4 / ( n3 * zxsa)
          ENDIF
       ENDIF
!
!--    Partitioning of H2O, HNO3, and NH3: Dissolutional growth
       IF ( lspartition  .AND.  ( pchno3 > 1.0E+10_wp  .OR.  pc_nh3 > 1.0E+10_wp ) )  THEN
          CALL gpparthno3( ppres, ptemp, paero, pchno3, pc_nh3, pcw, pcs, zbeta, ptstep )
       ENDIF
    ENDIF
!
!-- Condensation of water vapour
    IF ( lscndh2oae )  THEN 
       CALL gpparth2o( paero, ptemp, ppres, pcs, pcw, ptstep )
    ENDIF

 END SUBROUTINE condensation

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculates the particle number and volume increase, and gas-phase
!> concentration decrease due to nucleation subsequent growth to detectable size 
!> of 3 nm.
!
!> Method:
!> When the formed clusters grow by condensation (possibly also by self-
!> coagulation), their number is reduced due to scavenging to pre-existing
!> particles. Thus, the apparent nucleation rate at 3 nm is significantly lower
!> than the real nucleation rate (at ~1 nm).
!
!> Calculation of the formation rate of detectable particles at 3 nm (i.e. J3):
!> nj3 = 1: Kerminen, V.-M. and Kulmala, M. (2002), J. Aerosol Sci.,33, 609-622.
!> nj3 = 2: Lehtinen et al. (2007), J. Aerosol Sci., 38(9), 988-994.
!> nj3 = 3: Anttila et al. (2010), J. Aerosol Sci., 41(7), 621-636.
!
!> c = aerosol of critical radius (1 nm)
!> x = aerosol with radius 3 nm
!> 2 = wet or mean droplet
!
!> Called from subroutine condensation (in module salsa_dynamics_mod.f90)
!
!> Calls one of the following subroutines:
!>  - binnucl
!>  - ternucl
!>  - kinnucl
!>  - actnucl
!
!> fxm: currently only sulphuric acid grows particles from 1 to 3 nm
!>  (if asked from Markku, this is terribly wrong!!!)
!
!> Coded by:
!> Hannele Korhonen (FMI) 2005
!> Harri Kokkola (FMI) 2006
!> Matti Niskanen(FMI) 2012
!> Anton Laakso  (FMI) 2013
!------------------------------------------------------------------------------!

 SUBROUTINE nucleation( paero, ptemp, prh, ppres, pc_sa, pc_ocnv, pc_nh3, ptstep, pj3n3, pxsa,     &
                        pxocnv )

    IMPLICIT NONE

    INTEGER(iwp) ::  iteration

    REAL(wp) ::  zc_h2so4     !< H2SO4 conc. (#/cm3) !UNITS!
    REAL(wp) ::  zc_org       !< organic vapour conc. (#/cm3)
    REAL(wp) ::  zcc_c        !< Cunningham correct factor for c = critical (1nm)
    REAL(wp) ::  zcc_x        !< Cunningham correct factor for x = 3nm
    REAL(wp) ::  zcoags_c     !< coagulation sink (1/s) for c = critical (1nm)
    REAL(wp) ::  zcoags_x     !< coagulation sink (1/s) for x = 3nm
    REAL(wp) ::  zcoagstot    !< total particle losses due to coagulation, including condensation
                              !< and self-coagulation
    REAL(wp) ::  zcocnv_local !< organic vapour conc. (#/m3)
    REAL(wp) ::  zcsink       !< condensational sink (#/m2)
    REAL(wp) ::  zcsa_local   !< H2SO4 conc. (#/m3)
    REAL(wp) ::  zcv_c        !< mean relative thermal velocity (m/s) for c = critical (1nm)
    REAL(wp) ::  zcv_x        !< mean relative thermal velocity (m/s) for x = 3nm
    REAL(wp) ::  zdcrit       !< diameter of critical cluster (m)
    REAL(wp) ::  zdelta_vap   !< change of H2SO4 and organic vapour concentration (#/m3)
    REAL(wp) ::  zdfvap       !< air diffusion coefficient (m2/s)
    REAL(wp) ::  zdmean       !< mean diameter of existing particles (m) 
    REAL(wp) ::  zeta         !< constant: proportional to ratio of CS/GR (m)
                              !< (condensation sink / growth rate)
    REAL(wp) ::  zgamma       !< proportionality factor ((nm2*m2)/h)
    REAL(wp) ::  z_gr_clust   !< growth rate of formed clusters (nm/h)
    REAL(wp) ::  z_gr_tot     !< total growth rate
    REAL(wp) ::  zj3          !< number conc. of formed 3nm particles (#/m3)
    REAL(wp) ::  zjnuc        !< nucleation rate at ~1nm (#/m3s)
    REAL(wp) ::  z_k_eff      !< effective cogulation coefficient for freshly nucleated particles
    REAL(wp) ::  zknud_c      !< Knudsen number for c = critical (1nm)
    REAL(wp) ::  zknud_x      !< Knudsen number for x = 3nm
    REAL(wp) ::  zkocnv       !< lever: zkocnv=1 --> organic compounds involved in nucleation
    REAL(wp) ::  zksa         !< lever: zksa=1 --> H2SO4 involved in nucleation
    REAL(wp) ::  zlambda      !< parameter for adjusting the growth rate due to self-coagulation
    REAL(wp) ::  zm_c         !< particle mass (kg) for c = critical (1nm)
    REAL(wp) ::  zm_para      !< Parameter m for calculating the coagulation sink (Eq. 5&6 in 
                              !< Lehtinen et al. 2007)
    REAL(wp) ::  zm_x         !< particle mass (kg) for x = 3nm
    REAL(wp) ::  zmfp         !< mean free path of condesing vapour(m)
    REAL(wp) ::  zmixnh3      !< ammonia mixing ratio (ppt)
    REAL(wp) ::  zmyy         !< gas dynamic viscosity (N*s/m2)
    REAL(wp) ::  z_n_nuc      !< number of clusters/particles at the size range d1-dx (#/m3)
    REAL(wp) ::  znoc         !< number of organic molecules in critical cluster
    REAL(wp) ::  znsa         !< number of H2SO4 molecules in critical cluster

    REAL(wp), INTENT(in) ::  pc_nh3   !< ammonia concentration (#/m3) 
    REAL(wp), INTENT(in) ::  pc_ocnv  !< conc. of non-volatile OC (#/m3)
    REAL(wp), INTENT(in) ::  pc_sa    !< sulphuric acid conc. (#/m3)
    REAL(wp), INTENT(in) ::  ppres    !< ambient air pressure (Pa)
    REAL(wp), INTENT(in) ::  prh      !< ambient rel. humidity [0-1]
    REAL(wp), INTENT(in) ::  ptemp    !< ambient temperature (K)
    REAL(wp), INTENT(in) ::  ptstep   !< time step (s) of SALSA

    REAL(wp), INTENT(inout) ::  pj3n3(2) !< formation mass rate of molecules (molec/m3s) for
                                         !< 1: H2SO4 and 2: organic vapour

    REAL(wp), INTENT(out) ::  pxocnv  !< ratio of non-volatile organic vapours in 3 nm particles
    REAL(wp), INTENT(out) ::  pxsa    !< ratio of H2SO4 in 3 nm aerosol particles

    REAL(wp), DIMENSION(nbins_aerosol) ::  zbeta       !< transitional correction factor
    REAL(wp), DIMENSION(nbins_aerosol) ::  zcc_2       !< Cunningham correct factor:2
    REAL(wp), DIMENSION(nbins_aerosol) ::  zcv_2       !< mean relative thermal velocity (m/s): 2
    REAL(wp), DIMENSION(nbins_aerosol) ::  zcv_c2      !< average velocity after coagulation: c & 2
    REAL(wp), DIMENSION(nbins_aerosol) ::  zcv_x2      !< average velocity after coagulation: x & 2
    REAL(wp), DIMENSION(nbins_aerosol) ::  zdc_2       !< particle diffusion coefficient (m2/s): 2
    REAL(wp), DIMENSION(nbins_aerosol) ::  zdc_c       !< particle diffusion coefficient (m2/s): c
    REAL(wp), DIMENSION(nbins_aerosol) ::  zdc_c2      !< sum of diffusion coef. for c and 2
    REAL(wp), DIMENSION(nbins_aerosol) ::  zdc_x       !< particle diffusion coefficient (m2/s): x
    REAL(wp), DIMENSION(nbins_aerosol) ::  zdc_x2      !< sum of diffusion coef. for: x & 2
    REAL(wp), DIMENSION(nbins_aerosol) ::  zgamma_f_2  !< zgamma_f for calculating zomega
    REAL(wp), DIMENSION(nbins_aerosol) ::  zgamma_f_c  !< zgamma_f for calculating zomega
    REAL(wp), DIMENSION(nbins_aerosol) ::  zgamma_f_x  !< zgamma_f for calculating zomega
    REAL(wp), DIMENSION(nbins_aerosol) ::  z_k_c2      !< coagulation coef. in the continuum
                                                       !< regime: c & 2
    REAL(wp), DIMENSION(nbins_aerosol) ::  z_k_x2      !< coagulation coef. in the continuum
                                                       !< regime: x & 2
    REAL(wp), DIMENSION(nbins_aerosol) ::  zknud       !< particle Knudsen number
    REAL(wp), DIMENSION(nbins_aerosol) ::  zknud_2     !< particle Knudsen number: 2
    REAL(wp), DIMENSION(nbins_aerosol) ::  zm_2        !< particle mass (kg): 2
    REAL(wp), DIMENSION(nbins_aerosol) ::  zomega_2c   !< zomega (m) for calculating zsigma: c & 2
    REAL(wp), DIMENSION(nbins_aerosol) ::  zomega_2x   !< zomega (m) for calculating zsigma: x & 2
    REAL(wp), DIMENSION(nbins_aerosol) ::  zomega_c    !< zomega (m) for calculating zsigma: c
    REAL(wp), DIMENSION(nbins_aerosol) ::  zomega_x    !< zomega (m) for calculating zsigma: x
    REAL(wp), DIMENSION(nbins_aerosol) ::  z_r_c2      !< sum of the radii: c & 2
    REAL(wp), DIMENSION(nbins_aerosol) ::  z_r_x2      !< sum of the radii: x & 2
    REAL(wp), DIMENSION(nbins_aerosol) ::  zsigma_c2   !<
    REAL(wp), DIMENSION(nbins_aerosol) ::  zsigma_x2   !<

    TYPE(t_section), DIMENSION(nbins_aerosol), INTENT(inout) ::  paero  !< aerosol properties
!
!-- 1) Nucleation rate (zjnuc) and diameter of critical cluster (zdcrit)
    zjnuc  = 0.0_wp
    znsa   = 0.0_wp
    znoc   = 0.0_wp
    zdcrit = 0.0_wp
    zksa   = 0.0_wp
    zkocnv = 0.0_wp

    zc_h2so4 = pc_sa * 1.0E-6_wp   ! sulphuric acid conc. to #/cm3
    zc_org   = pc_ocnv * 1.0E-6_wp   ! conc. of non-volatile OC to #/cm3
    zmixnh3  = pc_nh3 * ptemp * argas / ( ppres * avo )

    SELECT CASE ( nsnucl )
!
!--    Binary H2SO4-H2O nucleation
       CASE(1)

          CALL binnucl( zc_h2so4, ptemp, prh, zjnuc, znsa, znoc, zdcrit,  zksa, zkocnv )
!
!--    Activation type nucleation (See Riipinen et al. (2007), Atmos. Chem. Phys., 7(8), 1899-1914)
       CASE(2)
!
!--       Nucleation rate (#/(m3 s))
          zc_h2so4  = MAX( zc_h2so4, 1.0E4_wp  )
          zc_h2so4  = MIN( zc_h2so4, 1.0E11_wp )
          zjnuc = act_coeff * pc_sa  ! (#/(m3 s))
!
!--       Organic compounds not involved when kinetic nucleation is assumed.
          zdcrit  = 7.9375E-10_wp   ! (m)
          zkocnv  = 0.0_wp
          zksa    = 1.0_wp
          znoc    = 0.0_wp
          znsa    = 2.0_wp
!
!--    Kinetically limited nucleation of (NH4)HSO4 clusters
!--    (See Sihto et al. (2006), Atmos. Chem. Phys., 6(12), 4079-4091.)
       CASE(3)
!
!--       Nucleation rate = coagcoeff*zpcsa**2 (#/(m3 s))
          zc_h2so4  = MAX( zc_h2so4, 1.0E4_wp  )
          zc_h2so4  = MIN( zc_h2so4, 1.0E11_wp )
          zjnuc = 5.0E-13_wp * zc_h2so4**2.0_wp * 1.0E+6_wp
!
!--       Organic compounds not involved when kinetic nucleation is assumed.
          zdcrit  = 7.9375E-10_wp   ! (m)
          zkocnv  = 0.0_wp
          zksa    = 1.0_wp
          znoc    = 0.0_wp
          znsa    = 2.0_wp
!
!--    Ternary H2SO4-H2O-NH3 nucleation
       CASE(4)

          CALL ternucl( zc_h2so4, zmixnh3, ptemp, prh, zjnuc, znsa, znoc, zdcrit, zksa, zkocnv )
!
!--    Organic nucleation, J~[ORG] or J~[ORG]**2
!--    (See Paasonen et al. (2010), Atmos. Chem. Phys., 10, 11223-11242.)
       CASE(5)
!
!--       Homomolecular nuleation rate
          zjnuc = 1.3E-7_wp * pc_ocnv   ! (1/s) (Paasonen et al. Table 4: median a_org)
!
!--       H2SO4 not involved when pure organic nucleation is assumed.
          zdcrit  = 1.5E-9  ! (m)
          zkocnv  = 1.0_wp
          zksa    = 0.0_wp
          znoc    = 1.0_wp
          znsa    = 0.0_wp
!
!--    Sum of H2SO4 and organic activation type nucleation, J~[H2SO4]+[ORG]
!--    (See Paasonen et al. (2010), Atmos. Chem. Phys., 10, 11223-11242)
       CASE(6)
!
!--       Nucleation rate  (#/m3/s)
          zjnuc = 6.1E-7_wp * pc_sa + 0.39E-7_wp * pc_ocnv   ! (Paasonen et al. Table 3.)
!
!--       Both organic compounds and H2SO4 are involved when sumnucleation is assumed.
          zdcrit  = 1.5E-9_wp   ! (m)
          zkocnv  = 1.0_wp
          zksa    = 1.0_wp
          znoc    = 1.0_wp
          znsa    = 1.0_wp
!
!--    Heteromolecular nucleation, J~[H2SO4]*[ORG] 
!--    (See Paasonen et al. (2010), Atmos. Chem. Phys., 10, 11223-11242.)
       CASE(7)
!
!--       Nucleation rate (#/m3/s)
          zjnuc = 4.1E-14_wp * pc_sa * pc_ocnv * 1.0E6_wp   ! (Paasonen et al. Table 4: median)
!
!--       Both organic compounds and H2SO4 are involved when heteromolecular nucleation is assumed
          zdcrit  = 1.5E-9_wp   ! (m)
          zkocnv  = 1.0_wp
          zksa    = 1.0_wp
          znoc    = 1.0_wp
          znsa    = 1.0_wp
!
!--    Homomolecular nucleation of H2SO4 and heteromolecular nucleation of H2SO4 and organic vapour,
!--    J~[H2SO4]**2 + [H2SO4]*[ORG] (EUCAARI project)
!--    (See Paasonen et al. (2010), Atmos. Chem. Phys., 10, 11223-11242)
       CASE(8)
!
!--       Nucleation rate (#/m3/s)
          zjnuc = ( 1.1E-14_wp * zc_h2so4**2 + 3.2E-14_wp * zc_h2so4 * zc_org ) * 1.0E+6_wp
!
!--       Both organic compounds and H2SO4 are involved when SAnucleation is assumed
          zdcrit  = 1.5E-9_wp   ! (m)
          zkocnv  = 1.0_wp
          zksa    = 1.0_wp
          znoc    = 1.0_wp
          znsa    = 3.0_wp
!
!--    Homomolecular nucleation of H2SO4 and organic vapour and heteromolecular nucleation of H2SO4
!--    and organic vapour, J~[H2SO4]**2 + [H2SO4]*[ORG]+[ORG]**2 (EUCAARI project)
       CASE(9)
!
!--       Nucleation rate (#/m3/s)
          zjnuc = ( 1.4E-14_wp * zc_h2so4**2 + 2.6E-14_wp * zc_h2so4 * zc_org + 0.037E-14_wp *     &
                    zc_org**2 ) * 1.0E+6_wp
!
!--       Both organic compounds and H2SO4 are involved when SAORGnucleation is assumed
          zdcrit  = 1.5E-9_wp   ! (m)
          zkocnv  = 1.0_wp
          zksa    = 1.0_wp
          znoc    = 3.0_wp
          znsa    = 3.0_wp

    END SELECT

    zcsa_local = pc_sa
    zcocnv_local = pc_ocnv
!
!-- 2) Change of particle and gas concentrations due to nucleation
!
!-- 2.1) Check that there is enough H2SO4 and organic vapour to produce the nucleation
    IF ( nsnucl <= 4 )  THEN 
!
!--    If the chosen nucleation scheme is 1-4, nucleation occurs only due to H2SO4. All of the total
!--    vapour concentration that is taking part to the nucleation is there for sulphuric acid
!--    (sa = H2SO4) and non-volatile organic vapour is zero.
       pxsa   = 1.0_wp   ! ratio of sulphuric acid in 3nm particles
       pxocnv = 0.0_wp   ! ratio of non-volatile origanic vapour
                                ! in 3nm particles
    ELSEIF ( nsnucl > 4 )  THEN
!
!--    If the chosen nucleation scheme is 5-9, nucleation occurs due to organic vapour or the
!--    combination of organic vapour and H2SO4. The number of needed molecules depends on the chosen
!--    nucleation type and it has an effect also on the minimum ratio of the molecules present.
       IF ( pc_sa * znsa + pc_ocnv * znoc < 1.E-14_wp )  THEN
          pxsa   = 0.0_wp
          pxocnv = 0.0_wp
       ELSE
          pxsa   = pc_sa * znsa / ( pc_sa * znsa + pc_ocnv * znoc ) 
          pxocnv = pc_ocnv * znoc / ( pc_sa * znsa + pc_ocnv * znoc )
       ENDIF
    ENDIF
!
!-- The change in total vapour concentration is the sum of the concentrations of the vapours taking
!-- part to the nucleation (depends on the chosen nucleation scheme)
    zdelta_vap = MIN( zjnuc * ( znoc + znsa ), ( pc_ocnv * zkocnv + pc_sa * zksa ) / ptstep )
!
!-- Nucleation rate J at ~1nm (#/m3s)
    zjnuc = zdelta_vap / ( znoc + znsa )
!
!-- H2SO4 concentration after nucleation (#/m3)
    zcsa_local = MAX( 1.0_wp, pc_sa - zdelta_vap * pxsa )
!
!-- Non-volative organic vapour concentration after nucleation (#/m3)
    zcocnv_local = MAX( 1.0_wp, pc_ocnv - zdelta_vap * pxocnv )
!
!-- 2.2) Formation rate of 3 nm particles (Kerminen & Kulmala, 2002)
!
!-- Growth rate by H2SO4 and organic vapour (nm/h, Eq. 21)
    z_gr_clust = 2.3623E-15_wp * SQRT( ptemp ) * ( zcsa_local + zcocnv_local )
!
!-- 2.2.2) Condensational sink of pre-existing particle population
!
!-- Diffusion coefficient (m2/s)
    zdfvap = 5.1111E-10_wp * ptemp**1.75_wp * ( p_0 + 1325.0_wp ) / ppres
!
!-- Mean free path of condensing vapour (m) (Jacobson (2005), Eq. 15.25 and 16.29)
    zmfp = 3.0_wp * zdfvap * SQRT( pi * amh2so4 / ( 8.0_wp * argas * ptemp ) )
!
!-- Knudsen number
    zknud = 2.0_wp * zmfp / ( paero(:)%dwet + d_sa )
!
!-- Transitional regime correction factor (zbeta) according to Fuchs and Sutugin (1971) (Eq. 4 in
!-- Kerminen and Kulmala, 2002)
    zbeta = ( zknud + 1.0_wp) / ( 0.377_wp * zknud + 1.0_wp + 4.0_wp / ( 3.0_wp * massacc ) *      &
            ( zknud + zknud**2 ) )
!
!-- Condensational sink (#/m2, Eq. 3)
    zcsink = SUM( paero(:)%dwet * zbeta * paero(:)%numc )
!
!-- 2.2.3) Parameterised formation rate of detectable 3 nm particles (i.e. J3)
    IF ( nj3 == 1 )  THEN   ! Kerminen and Kulmala (2002)
!
!--    Constants needed for the parameterisation: dapp = 3 nm and dens_nuc = 1830 kg/m3
       IF ( zcsink < 1.0E-30_wp )  THEN
          zeta = 0._dp
       ELSE
!
!--       Mean diameter of backgroud population (nm)
          zdmean = 1.0_wp / SUM( paero(:)%numc ) * SUM( paero(:)%numc * paero(:)%dwet ) * 1.0E+9_wp
!
!--       Proportionality factor (nm2*m2/h) (Eq. 22)
          zgamma = 0.23_wp * ( zdcrit * 1.0E+9_wp )**0.2_wp * ( zdmean / 150.0_wp )**0.048_wp *    &
                   ( ptemp / 293.0_wp )**( -0.75_wp ) * ( arhoh2so4 / 1000.0_wp )**( -0.33_wp )
!
!--       Factor eta (nm, Eq. 11)
          zeta = MIN( zgamma * zcsink / z_gr_clust, zdcrit * 1.0E11_wp )
       ENDIF
!
!--    Number conc. of clusters surviving to 3 nm in a time step (#/m3, Eq.14)
       zj3 = zjnuc * EXP( MIN( 0.0_wp, zeta / 3.0_wp - zeta / ( zdcrit * 1.0E9_wp ) ) )

    ELSEIF ( nj3 > 1 )  THEN   ! Lehtinen et al. (2007) or Anttila et al. (2010)
!
!--    Defining the parameter m (zm_para) for calculating the coagulation sink onto background
!--    particles (Eq. 5&6 in Lehtinen et al. 2007). The growth is investigated between
!--    [d1,reglim(1)] = [zdcrit,3nm] and m = LOG( CoagS_dx / CoagX_zdcrit ) / LOG( reglim / zdcrit )
!--    (Lehtinen et al. 2007, Eq. 6).
!--    The steps for the coagulation sink for reglim = 3nm and zdcrit ~= 1nm are explained in
!--    Kulmala et al. (2001). The particles of diameter zdcrit ~1.14 nm  and reglim = 3nm are both
!--    in turn the "number 1" variables (Kulmala et al. 2001).
!--    c = critical (1nm), x = 3nm, 2 = wet or mean droplet
!
!--    Sum of the radii, R12 = R1 + R2 (m) of two particles 1 and 2
       z_r_c2 = zdcrit / 2.0_wp + paero(:)%dwet / 2.0_wp
       z_r_x2 = reglim(1) / 2.0_wp + paero(:)%dwet / 2.0_wp
!
!--    Particle mass (kg) (comes only from H2SO4)
       zm_c = 4.0_wp / 3.0_wp * pi * ( zdcrit / 2.0_wp )**3 * arhoh2so4
       zm_x = 4.0_wp / 3.0_wp * pi * ( reglim(1) / 2.0_wp )**3 * arhoh2so4
       zm_2 = 4.0_wp / 3.0_wp * pi * ( 0.5_wp * paero(:)%dwet )**3 * arhoh2so4
!
!--    Mean relative thermal velocity between the particles (m/s)
       zcv_c = SQRT( 8.0_wp * abo * ptemp / ( pi * zm_c ) )
       zcv_x = SQRT( 8.0_wp * abo * ptemp / ( pi * zm_x ) )
       zcv_2 = SQRT( 8.0_wp * abo * ptemp / ( pi * zm_2 ) )
!
!--    Average velocity after coagulation
       zcv_c2(:) = SQRT( zcv_c**2 + zcv_2**2 )
       zcv_x2(:) = SQRT( zcv_x**2 + zcv_2**2 )
!
!--    Knudsen number (zmfp = mean free path of condensing vapour)
       zknud_c = 2.0_wp * zmfp / zdcrit
       zknud_x = 2.0_wp * zmfp / reglim(1)
       zknud_2(:) = MAX( 0.0_wp, 2.0_wp * zmfp / paero(:)%dwet )
!
!--    Cunningham correction factors (Allen and Raabe, 1985)
       zcc_c    = 1.0_wp + zknud_c    * ( 1.142_wp + 0.558_wp * EXP( -0.999_wp / zknud_c ) )
       zcc_x    = 1.0_wp + zknud_x    * ( 1.142_wp + 0.558_wp * EXP( -0.999_wp / zknud_x ) )
       zcc_2(:) = 1.0_wp + zknud_2(:) * ( 1.142_wp + 0.558_wp * EXP( -0.999_wp / zknud_2(:) ) )
!
!--    Gas dynamic viscosity (N*s/m2). Here, viscocity(air @20C) = 1.81e-5_dp N/m2 *s (Hinds, p. 25)
       zmyy = 1.81E-5_wp * ( ptemp / 293.0_wp )**0.74_wp
!
!--    Particle diffusion coefficient (m2/s) (continuum regime)
       zdc_c(:) = abo * ptemp * zcc_c    / ( 3.0_wp * pi * zmyy * zdcrit )
       zdc_x(:) = abo * ptemp * zcc_x    / ( 3.0_wp * pi * zmyy * reglim(1) )
       zdc_2(:) = abo * ptemp * zcc_2(:) / ( 3.0_wp * pi * zmyy * paero(:)%dwet )
!
!--    D12 = D1+D2 (Seinfield and Pandis, 2nd ed. Eq. 13.38)
       zdc_c2 = zdc_c + zdc_2
       zdc_x2 = zdc_x + zdc_2
!
!--    zgamma_f = 8*D/pi/zcv (m) for calculating zomega (Fuchs, 1964)
       zgamma_f_c = 8.0_wp * zdc_c / pi / zcv_c
       zgamma_f_x = 8.0_wp * zdc_x / pi / zcv_x
       zgamma_f_2 = 8.0_wp * zdc_2 / pi / zcv_2
!
!--    zomega (m) for calculating zsigma
       zomega_c = ( ( z_r_c2 + zgamma_f_c )**3 - ( z_r_c2 ** 2 + zgamma_f_c )**1.5_wp ) /          &
                  ( 3.0_wp * z_r_c2 * zgamma_f_c ) - z_r_c2
       zomega_x = ( ( z_r_x2 + zgamma_f_x )**3 - ( z_r_x2**2 + zgamma_f_x )** 1.5_wp ) /           &
                  ( 3.0_wp * z_r_x2 * zgamma_f_x ) - z_r_x2
       zomega_2c = ( ( z_r_c2 + zgamma_f_2 )**3 - ( z_r_c2**2 + zgamma_f_2 )**1.5_wp ) /           &
                   ( 3.0_wp * z_r_c2 * zgamma_f_2 ) - z_r_c2
       zomega_2x = ( ( z_r_x2 + zgamma_f_2 )**3 - ( z_r_x2**2 + zgamma_f_2 )**1.5_wp ) /           &
                   ( 3.0_wp * z_r_x2 * zgamma_f_2 ) - z_r_x2 
!
!--    The distance (m) at which the two fluxes are matched (condensation and coagulation sinks)
       zsigma_c2 = SQRT( zomega_c**2 + zomega_2c**2 )
       zsigma_x2 = SQRT( zomega_x**2 + zomega_2x**2 )
!
!--    Coagulation coefficient in the continuum regime (m*m2/s, Eq. 17 in Kulmala et al., 2001)
       z_k_c2 = 4.0_wp * pi * z_r_c2 * zdc_c2 / ( z_r_c2 / ( z_r_c2 + zsigma_c2 ) +                &
               4.0_wp * zdc_c2 / ( zcv_c2 * z_r_c2 ) )
       z_k_x2 = 4.0_wp * pi * z_r_x2 * zdc_x2 / ( z_r_x2 / ( z_r_x2 + zsigma_x2 ) +                &
               4.0_wp * zdc_x2 / ( zcv_x2 * z_r_x2 ) )
!
!--    Coagulation sink (1/s, Eq. 16 in Kulmala et al., 2001)
       zcoags_c = MAX( 1.0E-20_wp, SUM( z_k_c2 * paero(:)%numc ) )
       zcoags_x = MAX( 1.0E-20_wp, SUM( z_k_x2 * paero(:)%numc ) )
!
!--    Parameter m for calculating the coagulation sink onto background particles (Eq. 5&6 in 
!--    Lehtinen et al. 2007)
       zm_para = LOG( zcoags_x / zcoags_c ) / LOG( reglim(1) / zdcrit )
!
!--    Parameter gamma for calculating the formation rate J of particles having
!--    a diameter zdcrit < d < reglim(1) (Anttila et al. 2010, eq. 5 or Lehtinen et al.,2007, eq. 7)
       zgamma = ( ( ( reglim(1) / zdcrit )**( zm_para + 1.0_wp ) ) - 1.0_wp ) / ( zm_para + 1.0_wp )

       IF ( nj3 == 2 )  THEN   ! Lehtinen et al. (2007): coagulation sink
!
!--       Formation rate J before iteration (#/m3s)
          zj3 = zjnuc * EXP( MIN( 0.0_wp, -zgamma * zdcrit * zcoags_c / ( z_gr_clust * 1.0E-9_wp / &
                60.0_wp**2 ) ) )

       ELSEIF ( nj3 == 3 )  THEN  ! Anttila et al. (2010): coagulation sink and self-coag.
!
!--       If air is polluted, the self-coagulation becomes important. Self-coagulation of small 
!--       particles < 3 nm.
!
!--       "Effective" coagulation coefficient between freshly-nucleated particles:
          z_k_eff = 5.0E-16_wp   ! m3/s
!
!--       zlambda parameter for "adjusting" the growth rate due to the self-coagulation
          zlambda = 6.0_wp

          IF ( reglim(1) >= 10.0E-9_wp )  THEN   ! for particles >10 nm:
             z_k_eff   = 5.0E-17_wp
             zlambda = 3.0_wp
          ENDIF
!
!--       Initial values for coagulation sink and growth rate  (m/s)
          zcoagstot = zcoags_c
          z_gr_tot = z_gr_clust * 1.0E-9_wp / 60.0_wp**2
!
!--       Number of clusters/particles at the size range [d1,dx] (#/m3):
          z_n_nuc = zjnuc / zcoagstot !< Initial guess
!
!--       Coagulation sink and growth rate due to self-coagulation:
          DO  iteration = 1, 5
             zcoagstot = zcoags_c + z_k_eff * z_n_nuc * 1.0E-6_wp   ! (1/s, Anttila et al., eq. 1)
             z_gr_tot = z_gr_clust * 2.77777777E-7_wp +  1.5708E-6_wp * zlambda * zdcrit**3 *      &
                      ( z_n_nuc * 1.0E-6_wp ) * zcv_c * avo * 2.77777777E-7_wp ! (Eq. 3)
             zeta = - zcoagstot / ( ( zm_para + 1.0_wp ) * z_gr_tot * ( zdcrit**zm_para ) ) ! (Eq.7b)
!
!--          Calculate Eq. 7a (Taylor series for the number of particles between [d1,dx])
             z_n_nuc =  z_n_nuc_tayl( zdcrit, reglim(1), zm_para, zjnuc, zeta, z_gr_tot )
          ENDDO
!
!--       Calculate the final values with new z_n_nuc:
          zcoagstot = zcoags_c + z_k_eff * z_n_nuc * 1.0E-6_wp   ! (1/s)
          z_gr_tot = z_gr_clust * 1.0E-9_wp / 3600.0_wp + 1.5708E-6_wp *  zlambda * zdcrit**3 *    &
                   ( z_n_nuc * 1.0E-6_wp ) * zcv_c * avo * 1.0E-9_wp / 3600.0_wp !< (m/s)
          zj3 = zjnuc * EXP( MIN( 0.0_wp, -zgamma * zdcrit * zcoagstot / z_gr_tot ) ) ! (#/m3s, Eq.5a)

       ENDIF
    ENDIF
!
!-- If J3 very small (< 1 #/cm3), neglect particle formation. In real atmosphere this would mean
!-- that clusters form but coagulate to pre-existing particles who gain sulphate. Since
!-- CoagS ~ CS (4piD*CS'), we do *not* update H2SO4 concentration here but let condensation take
!-- care of it. Formation mass rate of molecules (molec/m3s) for 1: H2SO4 and 2: organic vapour
    pj3n3(1) = zj3 * n3 * pxsa
    pj3n3(2) = zj3 * n3 * pxocnv

 END SUBROUTINE nucleation

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate the nucleation rate and the size of critical clusters assuming
!> binary nucleation.
!> Parametrisation according to Vehkamaki et al. (2002), J. Geophys. Res.,
!> 107(D22), 4622. Called from subroutine nucleation.
!------------------------------------------------------------------------------!
 SUBROUTINE binnucl( pc_sa, ptemp, prh, pnuc_rate, pn_crit_sa, pn_crit_ocnv, pd_crit, pk_sa,       &
                     pk_ocnv )

    IMPLICIT NONE

    REAL(wp) ::  za      !<
    REAL(wp) ::  zb      !<
    REAL(wp) ::  zc      !<
    REAL(wp) ::  zcoll   !<
    REAL(wp) ::  zlogsa  !<  LOG( zpcsa )
    REAL(wp) ::  zlogrh  !<  LOG( zrh )
    REAL(wp) ::  zm1     !<
    REAL(wp) ::  zm2     !<
    REAL(wp) ::  zma     !<
    REAL(wp) ::  zmw     !<
    REAL(wp) ::  zntot   !< number of molecules in critical cluster
    REAL(wp) ::  zpcsa   !< sulfuric acid concentration
    REAL(wp) ::  zrh     !< relative humidity
    REAL(wp) ::  zroo    !<
    REAL(wp) ::  zt      !< temperature
    REAL(wp) ::  zv1     !<
    REAL(wp) ::  zv2     !<
    REAL(wp) ::  zx      !< mole fraction of sulphate in critical cluster
    REAL(wp) ::  zxmass  !<

    REAL(wp), INTENT(in) ::   pc_sa   !< H2SO4 conc. (#/cm3)
    REAL(wp), INTENT(in) ::   prh     !< relative humidity [0-1
    REAL(wp), INTENT(in) ::   ptemp   !< ambient temperature (K)

    REAL(wp), INTENT(out) ::  pnuc_rate     !< nucleation rate (#/(m3 s))
    REAL(wp), INTENT(out) ::  pn_crit_sa    !< number of H2SO4 molecules in cluster (#)
    REAL(wp), INTENT(out) ::  pn_crit_ocnv  !< number of organic molecules in cluster (#)
    REAL(wp), INTENT(out) ::  pd_crit       !< diameter of critical cluster (m)
    REAL(wp), INTENT(out) ::  pk_sa         !< Lever: if pk_sa = 1, H2SO4 is involved in nucleation.
    REAL(wp), INTENT(out) ::  pk_ocnv       !< Lever: if pk_ocnv = 1, organic compounds are involved

    pnuc_rate = 0.0_wp
    pd_crit   = 1.0E-9_wp
!
!-- 1) Checking that we are in the validity range of the parameterization
    zpcsa  = MAX( pc_sa, 1.0E4_wp  )
    zpcsa  = MIN( zpcsa, 1.0E11_wp )
    zrh    = MAX( prh,   0.0001_wp )
    zrh    = MIN( zrh,   1.0_wp    )
    zt     = MAX( ptemp, 190.15_wp )
    zt     = MIN( zt,    300.15_wp )

    zlogsa = LOG( zpcsa )
    zlogrh   = LOG( prh )
!
!-- 2) Mole fraction of sulphate in a critical cluster (Eq. 11)
    zx = 0.7409967177282139_wp                  - 0.002663785665140117_wp * zt +                   &
         0.002010478847383187_wp * zlogrh       - 0.0001832894131464668_wp* zt * zlogrh +          &
         0.001574072538464286_wp * zlogrh**2    - 0.00001790589121766952_wp * zt * zlogrh**2 +     &
         0.0001844027436573778_wp * zlogrh**3   - 1.503452308794887E-6_wp * zt * zlogrh**3 -       &
         0.003499978417957668_wp * zlogsa     + 0.0000504021689382576_wp * zt * zlogsa
!
!-- 3) Nucleation rate (Eq. 12)
    pnuc_rate = 0.1430901615568665_wp + 2.219563673425199_wp * zt -                                &
                0.02739106114964264_wp * zt**2 + 0.00007228107239317088_wp * zt**3 +               &
                5.91822263375044_wp / zx + 0.1174886643003278_wp * zlogrh +                        &
                0.4625315047693772_wp * zt * zlogrh - 0.01180591129059253_wp * zt**2 * zlogrh +    &
                0.0000404196487152575_wp * zt**3 * zlogrh +                                        &
                ( 15.79628615047088_wp * zlogrh ) / zx - 0.215553951893509_wp * zlogrh**2 -        &
                0.0810269192332194_wp * zt * zlogrh**2 +                                           &
                0.001435808434184642_wp * zt**2 * zlogrh**2 -                                      &
                4.775796947178588E-6_wp * zt**3 * zlogrh**2 -                                      &
                ( 2.912974063702185_wp * zlogrh**2 ) / zx - 3.588557942822751_wp * zlogrh**3 +     &
                0.04950795302831703_wp * zt * zlogrh**3 -                                          &
                0.0002138195118737068_wp * zt**2 * zlogrh**3 +                                     &
                3.108005107949533E-7_wp * zt**3 * zlogrh**3 -                                      &
                ( 0.02933332747098296_wp * zlogrh**3 ) / zx + 1.145983818561277_wp * zlogsa -      &
                0.6007956227856778_wp * zt * zlogsa + 0.00864244733283759_wp * zt**2 * zlogsa -    &
                0.00002289467254710888_wp * zt**3 * zlogsa -                                       &
                ( 8.44984513869014_wp * zlogsa ) / zx + 2.158548369286559_wp * zlogrh * zlogsa +   &
                0.0808121412840917_wp * zt * zlogrh * zlogsa -                                     &
                0.0004073815255395214_wp * zt**2 * zlogrh * zlogsa -                               &
                4.019572560156515E-7_wp * zt**3 * zlogrh * zlogsa +                                &
                ( 0.7213255852557236_wp * zlogrh * zlogsa ) / zx +                                 &
                1.62409850488771_wp * zlogrh**2 * zlogsa -                                         &
                0.01601062035325362_wp * zt * zlogrh**2 * zlogsa +                                 &
                0.00003771238979714162_wp*zt**2* zlogrh**2 * zlogsa +                              &
                3.217942606371182E-8_wp * zt**3 * zlogrh**2 * zlogsa -                             &
                ( 0.01132550810022116_wp * zlogrh**2 * zlogsa ) / zx +                             &
                9.71681713056504_wp * zlogsa**2 - 0.1150478558347306_wp * zt * zlogsa**2 +         &
                0.0001570982486038294_wp * zt**2 * zlogsa**2 +                                     &
                4.009144680125015E-7_wp * zt**3 * zlogsa**2 +                                      &
                ( 0.7118597859976135_wp * zlogsa**2 ) / zx -                                       &
                1.056105824379897_wp * zlogrh * zlogsa**2 +                                        &
                0.00903377584628419_wp * zt * zlogrh * zlogsa**2 -                                 &
                0.00001984167387090606_wp * zt**2 * zlogrh * zlogsa**2 +                           &
                2.460478196482179E-8_wp * zt**3 * zlogrh * zlogsa**2 -                             &
                ( 0.05790872906645181_wp * zlogrh * zlogsa**2 ) / zx -                             &
                0.1487119673397459_wp * zlogsa**3 + 0.002835082097822667_wp * zt * zlogsa**3 -     &
                9.24618825471694E-6_wp * zt**2 * zlogsa**3 +                                       &
                5.004267665960894E-9_wp * zt**3 * zlogsa**3 -                                      &
                ( 0.01270805101481648_wp * zlogsa**3 ) / zx
!
!-- Nucleation rate in #/(cm3 s)
    pnuc_rate = EXP( pnuc_rate ) 
!
!-- Check the validity of parameterization
    IF ( pnuc_rate < 1.0E-7_wp )  THEN 
       pnuc_rate = 0.0_wp
       pd_crit   = 1.0E-9_wp
    ENDIF
!
!-- 4) Total number of molecules in the critical cluster (Eq. 13)
    zntot = - 0.002954125078716302_wp - 0.0976834264241286_wp * zt +                               &
              0.001024847927067835_wp * zt**2 - 2.186459697726116E-6_wp * zt**3 -                  &
              0.1017165718716887_wp / zx - 0.002050640345231486_wp * zlogrh -                      &
              0.007585041382707174_wp * zt * zlogrh + 0.0001926539658089536_wp * zt**2 * zlogrh -  &
              6.70429719683894E-7_wp * zt**3 * zlogrh - ( 0.2557744774673163_wp * zlogrh ) / zx +  &
              0.003223076552477191_wp * zlogrh**2 + 0.000852636632240633_wp * zt * zlogrh**2 -     &
              0.00001547571354871789_wp * zt**2 * zlogrh**2 +                                      &
              5.666608424980593E-8_wp * zt**3 * zlogrh**2 +                                        &
              ( 0.03384437400744206_wp * zlogrh**2 ) / zx +                                        &
              0.04743226764572505_wp * zlogrh**3 - 0.0006251042204583412_wp * zt * zlogrh**3 +     &
              2.650663328519478E-6_wp * zt**2 * zlogrh**3 -                                        &
              3.674710848763778E-9_wp * zt**3 * zlogrh**3 -                                        &
              ( 0.0002672510825259393_wp * zlogrh**3 ) / zx - 0.01252108546759328_wp * zlogsa +    &
              0.005806550506277202_wp * zt * zlogsa - 0.0001016735312443444_wp * zt**2 * zlogsa +  &
              2.881946187214505E-7_wp * zt**3 * zlogsa + ( 0.0942243379396279_wp * zlogsa ) / zx - &
              0.0385459592773097_wp * zlogrh * zlogsa -                                            &
              0.0006723156277391984_wp * zt * zlogrh * zlogsa  +                                   &
              2.602884877659698E-6_wp * zt**2 * zlogrh * zlogsa +                                  &
              1.194163699688297E-8_wp * zt**3 * zlogrh * zlogsa -                                  &
              ( 0.00851515345806281_wp * zlogrh * zlogsa ) / zx -                                  &
              0.01837488495738111_wp * zlogrh**2 * zlogsa +                                        &
              0.0001720723574407498_wp * zt * zlogrh**2 * zlogsa -                                 &
              3.717657974086814E-7_wp * zt**2 * zlogrh**2 * zlogsa -                               &
              5.148746022615196E-10_wp * zt**3 * zlogrh**2 * zlogsa +                              &
              ( 0.0002686602132926594_wp * zlogrh**2 * zlogsa ) / zx -                             &
              0.06199739728812199_wp * zlogsa**2 + 0.000906958053583576_wp * zt * zlogsa**2 -      &
              9.11727926129757E-7_wp * zt**2 * zlogsa**2 -                                         &
              5.367963396508457E-9_wp * zt**3 * zlogsa**2 -                                        &
              ( 0.007742343393937707_wp * zlogsa**2 ) / zx +                                       &
              0.0121827103101659_wp * zlogrh * zlogsa**2 -                                         &
              0.0001066499571188091_wp * zt * zlogrh * zlogsa**2 +                                 &
              2.534598655067518E-7_wp * zt**2 * zlogrh * zlogsa**2 -                               &
              3.635186504599571E-10_wp * zt**3 * zlogrh * zlogsa**2 +                              &
              ( 0.0006100650851863252_wp * zlogrh * zlogsa **2 ) / zx +                            &
              0.0003201836700403512_wp * zlogsa**3 - 0.0000174761713262546_wp * zt * zlogsa**3 +   &
              6.065037668052182E-8_wp * zt**2 * zlogsa**3 -                                        &
              1.421771723004557E-11_wp * zt**3 * zlogsa**3 +                                       &
              ( 0.0001357509859501723_wp * zlogsa**3 ) / zx
    zntot = EXP( zntot )  ! in #
!
!-- 5) Size of the critical cluster pd_crit (m) (diameter) (Eq. 14)
    pn_crit_sa = zx * zntot
    pd_crit = 2.0E-9_wp * EXP( -1.6524245_wp + 0.42316402_wp * zx + 0.33466487_wp * LOG( zntot ) )
!
!-- 6) Organic compounds not involved when binary nucleation is assumed
    pn_crit_ocnv = 0.0_wp   ! number of organic molecules 
    pk_sa        = 1.0_wp   ! if = 1, H2SO4 involved in nucleation
    pk_ocnv      = 0.0_wp   ! if = 1, organic compounds involved
!
!-- Set nucleation rate to collision rate
    IF ( pn_crit_sa < 4.0_wp ) THEN
!
!--    Volumes of the colliding objects
       zma    = 96.0_wp   ! molar mass of SO4 in g/mol
       zmw    = 18.0_wp   ! molar mass of water in g/mol
       zxmass = 1.0_wp    ! mass fraction of H2SO4
       za = 0.7681724_wp + zxmass * ( 2.1847140_wp + zxmass *                                      &
                                      ( 7.1630022_wp + zxmass *                                    &
                                        ( -44.31447_wp + zxmass *                                  &
                                          ( 88.75606 + zxmass *                                    &
                                            ( -75.73729_wp + zxmass * 23.43228_wp ) ) ) ) )
       zb = 1.808225E-3_wp + zxmass * ( -9.294656E-3_wp + zxmass *                                 &
                                        ( -0.03742148_wp + zxmass *                                &
                                          ( 0.2565321_wp + zxmass *                                &
                                            ( -0.5362872_wp + zxmass *                             &
                                              ( 0.4857736 - zxmass * 0.1629592_wp ) ) ) ) )
       zc = - 3.478524E-6_wp + zxmass * ( 1.335867E-5_wp + zxmass *                                &
                                          ( 5.195706E-5_wp + zxmass *                              &
                                            ( -3.717636E-4_wp + zxmass *                           &
                                              ( 7.990811E-4_wp + zxmass *                          &
                                                ( -7.458060E-4_wp + zxmass * 2.58139E-4_wp ) ) ) ) )
!
!--    Density for the sulphuric acid solution (Eq. 10 in Vehkamaki)
       zroo = ( za + zt * ( zb + zc * zt ) ) * 1.0E+3_wp   ! (kg/m^3
       zm1  = 0.098_wp   ! molar mass of H2SO4 in kg/mol
       zm2  = zm1
       zv1  = zm1 / avo / zroo   ! volume
       zv2  = zv1
!
!--    Collision rate
       zcoll =  zpcsa * zpcsa * ( 3.0_wp * pi / 4.0_wp )**0.16666666_wp *                          &
                SQRT( 6.0_wp * argas * zt / zm1 + 6.0_wp * argas * zt / zm2 ) *                    &
                ( zv1**0.33333333_wp + zv2**0.33333333_wp )**2 * 1.0E+6_wp    ! m3 -> cm3
       zcoll = MIN( zcoll, 1.0E+10_wp )
       pnuc_rate  = zcoll   ! (#/(cm3 s))

    ELSE
       pnuc_rate  = MIN( pnuc_rate, 1.0E+10_wp )
    ENDIF
    pnuc_rate = pnuc_rate * 1.0E+6_wp   ! (#/(m3 s))

 END SUBROUTINE binnucl
 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate the nucleation rate and the size of critical clusters assuming
!> ternary nucleation. Parametrisation according to:
!> Napari et al. (2002), J. Chem. Phys., 116, 4221-4227 and
!> Napari et al. (2002), J. Geophys. Res., 107(D19), AAC 6-1-ACC 6-6.
!------------------------------------------------------------------------------!
 SUBROUTINE ternucl( pc_sa, pc_nh3, ptemp, prh, pnuc_rate, pn_crit_sa, pn_crit_ocnv, pd_crit,      &
                     pk_sa, pk_ocnv )

    IMPLICIT NONE

    REAL(wp) ::  zlnj     !< logarithm of nucleation rate
    REAL(wp) ::  zlognh3  !< LOG( pc_nh3 )
    REAL(wp) ::  zlogrh   !< LOG( prh )
    REAL(wp) ::  zlogsa   !< LOG( pc_sa )

    REAL(wp), INTENT(in) ::   pc_nh3  !< ammonia mixing ratio (ppt)
    REAL(wp), INTENT(in) ::   pc_sa   !< H2SO4 conc. (#/cm3)
    REAL(wp), INTENT(in) ::   prh     !< relative humidity [0-1]
    REAL(wp), INTENT(in) ::   ptemp   !< ambient temperature (K)

    REAL(wp), INTENT(out) ::  pd_crit  !< diameter of critical cluster (m)
    REAL(wp), INTENT(out) ::  pk_ocnv  !< if pk_ocnv = 1, organic compounds participate in nucleation
    REAL(wp), INTENT(out) ::  pk_sa    !< if pk_sa = 1, H2SO4 participate in nucleation
    REAL(wp), INTENT(out) ::  pn_crit_ocnv  !< number of organic molecules in cluster (#)
    REAL(wp), INTENT(out) ::  pn_crit_sa    !< number of H2SO4 molecules in cluster (#)
    REAL(wp), INTENT(out) ::  pnuc_rate     !< nucleation rate (#/(m3 s))
!
!-- 1) Checking that we are in the validity range of the parameterization.
!--    Validity of parameterization : DO NOT REMOVE!
    IF ( ptemp < 240.0_wp  .OR.  ptemp > 300.0_wp )  THEN
       message_string = 'Invalid input value: ptemp'
       CALL message( 'salsa_mod: ternucl', 'PA0689', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( prh < 0.05_wp  .OR.  prh > 0.95_wp )  THEN
       message_string = 'Invalid input value: prh'
       CALL message( 'salsa_mod: ternucl', 'PA0649', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( pc_sa < 1.0E+4_wp  .OR.  pc_sa > 1.0E+9_wp )  THEN
       message_string = 'Invalid input value: pc_sa'
       CALL message( 'salsa_mod: ternucl', 'PA0650', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( pc_nh3 < 0.1_wp  .OR.  pc_nh3 > 100.0_wp )  THEN
       message_string = 'Invalid input value: pc_nh3'
       CALL message( 'salsa_mod: ternucl', 'PA0651', 1, 2, 0, 6, 0 )
    ENDIF

    zlognh3 = LOG( pc_nh3 )
    zlogrh  = LOG( prh )
    zlogsa  = LOG( pc_sa )
!
!-- 2) Nucleation rate (Eq. 7 in Napari et al., 2002: Parameterization of 
!--    ternary nucleation of sulfuric acid - ammonia - water.
    zlnj = - 84.7551114741543_wp + 0.3117595133628944_wp * prh +                                   &
           1.640089605712946_wp * prh * ptemp - 0.003438516933381083_wp * prh * ptemp**2 -         &
           0.00001097530402419113_wp * prh * ptemp**3 - 0.3552967070274677_wp / zlogsa -           &
           ( 0.06651397829765026_wp * prh ) / zlogsa - ( 33.84493989762471_wp * ptemp ) / zlogsa - &
           ( 7.823815852128623_wp * prh * ptemp ) / zlogsa +                                       &
           ( 0.3453602302090915_wp * ptemp**2 ) / zlogsa +                                         &
           ( 0.01229375748100015_wp * prh * ptemp**2 ) / zlogsa -                                  &
           ( 0.000824007160514956_wp *ptemp**3 ) / zlogsa +                                        &
           ( 0.00006185539100670249_wp * prh * ptemp**3 ) / zlogsa +                               &
           3.137345238574998_wp * zlogsa + 3.680240980277051_wp * prh * zlogsa -                   &
           0.7728606202085936_wp * ptemp * zlogsa - 0.204098217156962_wp * prh * ptemp * zlogsa +  &
           0.005612037586790018_wp * ptemp**2 * zlogsa +                                           &
           0.001062588391907444_wp * prh * ptemp**2 * zlogsa -                                     &
           9.74575691760229E-6_wp * ptemp**3 * zlogsa -                                            &
           1.265595265137352E-6_wp * prh * ptemp**3 * zlogsa + 19.03593713032114_wp * zlogsa**2 -  &
           0.1709570721236754_wp * ptemp * zlogsa**2 +                                             &
           0.000479808018162089_wp * ptemp**2 * zlogsa**2 -                                        &
           4.146989369117246E-7_wp * ptemp**3 * zlogsa**2 + 1.076046750412183_wp * zlognh3 +       &
           0.6587399318567337_wp * prh * zlognh3 + 1.48932164750748_wp * ptemp * zlognh3 +         &
           0.1905424394695381_wp * prh * ptemp * zlognh3 -                                         &
           0.007960522921316015_wp * ptemp**2 * zlognh3 -                                          &
           0.001657184248661241_wp * prh * ptemp**2 * zlognh3 +                                    &
           7.612287245047392E-6_wp * ptemp**3 * zlognh3 +                                          &
           3.417436525881869E-6_wp * prh * ptemp**3 * zlognh3 +                                    &
           ( 0.1655358260404061_wp * zlognh3 ) / zlogsa +                                          &
           ( 0.05301667612522116_wp * prh * zlognh3 ) / zlogsa +                                   &
           ( 3.26622914116752_wp * ptemp * zlognh3 ) / zlogsa -                                    &
           ( 1.988145079742164_wp * prh * ptemp * zlognh3 ) / zlogsa -                             &
           ( 0.04897027401984064_wp * ptemp**2 * zlognh3 ) / zlogsa +                              &
           ( 0.01578269253599732_wp * prh * ptemp**2 * zlognh3 ) / zlogsa +                        &
           ( 0.0001469672236351303_wp * ptemp**3 * zlognh3 ) / zlogsa -                            &
           ( 0.00002935642836387197_wp * prh * ptemp**3 *zlognh3 ) / zlogsa +                      &
           6.526451177887659_wp * zlogsa * zlognh3 -                                               &
           0.2580021816722099_wp * ptemp * zlogsa * zlognh3 +                                      &
           0.001434563104474292_wp * ptemp**2 * zlogsa * zlognh3 -                                 &
           2.020361939304473E-6_wp * ptemp**3 * zlogsa * zlognh3 -                                 &
           0.160335824596627_wp * zlogsa**2 * zlognh3 +                                            &
           0.00889880721460806_wp * ptemp * zlogsa**2 * zlognh3 -                                  &
           0.00005395139051155007_wp * ptemp**2 * zlogsa**2 * zlognh3 +                            &
           8.39521718689596E-8_wp * ptemp**3 * zlogsa**2 * zlognh3 +                               &
           6.091597586754857_wp * zlognh3**2 + 8.5786763679309_wp * prh * zlognh3**2 -             &
           1.253783854872055_wp * ptemp * zlognh3**2 -                                             &
           0.1123577232346848_wp * prh * ptemp * zlognh3**2 +                                      &
           0.00939835595219825_wp * ptemp**2 * zlognh3**2 +                                        &
           0.0004726256283031513_wp * prh * ptemp**2 * zlognh3**2 -                                &
           0.00001749269360523252_wp * ptemp**3 * zlognh3**2 -                                     &
           6.483647863710339E-7_wp * prh * ptemp**3 * zlognh3**2 +                                 &
           ( 0.7284285726576598_wp * zlognh3**2 ) / zlogsa +                                       &
           ( 3.647355600846383_wp * ptemp * zlognh3**2 ) / zlogsa -                                &
           ( 0.02742195276078021_wp * ptemp**2 * zlognh3**2 ) / zlogsa +                           &
           ( 0.00004934777934047135_wp * ptemp**3 * zlognh3**2 ) / zlogsa +                        &
           41.30162491567873_wp * zlogsa * zlognh3**2 -                                            &
           0.357520416800604_wp * ptemp * zlogsa * zlognh3**2 +                                    &
           0.000904383005178356_wp * ptemp**2 * zlogsa * zlognh3**2 -                              &
           5.737876676408978E-7_wp * ptemp**3 * zlogsa * zlognh3**2 -                              &
           2.327363918851818_wp * zlogsa**2 * zlognh3**2 +                                         &
           0.02346464261919324_wp * ptemp * zlogsa**2 * zlognh3**2 -                               &
           0.000076518969516405_wp * ptemp**2 * zlogsa**2 * zlognh3**2 +                           &
           8.04589834836395E-8_wp * ptemp**3 * zlogsa**2 * zlognh3**2 -                            &
           0.02007379204248076_wp * zlogrh - 0.7521152446208771_wp * ptemp * zlogrh +              &
           0.005258130151226247_wp * ptemp**2 * zlogrh -                                           &
           8.98037634284419E-6_wp * ptemp**3 * zlogrh +                                            &
           ( 0.05993213079516759_wp * zlogrh ) / zlogsa +                                          &
           ( 5.964746463184173_wp * ptemp * zlogrh ) / zlogsa -                                    &
           ( 0.03624322255690942_wp * ptemp**2 * zlogrh ) / zlogsa +                               &
           ( 0.00004933369382462509_wp * ptemp**3 * zlogrh ) / zlogsa -                            &
           0.7327310805365114_wp * zlognh3 * zlogrh -                                              &
           0.01841792282958795_wp * ptemp * zlognh3 * zlogrh +                                     &
           0.0001471855981005184_wp * ptemp**2 * zlognh3 * zlogrh -                                &
           2.377113195631848E-7_wp * ptemp**3 * zlognh3 * zlogrh
    pnuc_rate = EXP( zlnj )   ! (#/(cm3 s))
!
!-- Check validity of parametrization
    IF ( pnuc_rate < 1.0E-5_wp )  THEN
       pnuc_rate = 0.0_wp
       pd_crit   = 1.0E-9_wp
    ELSEIF ( pnuc_rate > 1.0E6_wp )  THEN
       message_string = 'Invalid output value: nucleation rate > 10^6 1/cm3s'
       CALL message( 'salsa_mod: ternucl', 'PA0623', 1, 2, 0, 6, 0 )
    ENDIF
    pnuc_rate = pnuc_rate * 1.0E6_wp   ! (#/(m3 s))
!
!-- 3) Number of H2SO4 molecules in a critical cluster (Eq. 9)
    pn_crit_sa = 38.16448247950508_wp + 0.7741058259731187_wp * zlnj +                             &
                 0.002988789927230632_wp * zlnj**2 - 0.3576046920535017_wp * ptemp -               &
                 0.003663583011953248_wp * zlnj * ptemp + 0.000855300153372776_wp * ptemp**2
!
!-- Kinetic limit: at least 2 H2SO4 molecules in a cluster
    pn_crit_sa = MAX( pn_crit_sa, 2.0E0_wp )
!
!-- 4) Size of the critical cluster in nm (Eq. 12)
    pd_crit = 0.1410271086638381_wp - 0.001226253898894878_wp * zlnj -                             &
              7.822111731550752E-6_wp * zlnj**2 - 0.001567273351921166_wp * ptemp -                &
              0.00003075996088273962_wp * zlnj * ptemp + 0.00001083754117202233_wp * ptemp**2
    pd_crit = pd_crit * 2.0E-9_wp   ! Diameter in m
!
!-- 5) Organic compounds not involved when ternary nucleation assumed
    pn_crit_ocnv = 0.0_wp
    pk_sa   = 1.0_wp
    pk_ocnv = 0.0_wp

 END SUBROUTINE ternucl

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Function z_n_nuc_tayl is connected to the calculation of self-coagualtion of 
!> small particles. It calculates number of the particles in the size range
!> [zdcrit,dx] using Taylor-expansion (please note that the expansion is not 
!> valid for certain rational numbers, e.g. -4/3 and -3/2)
!------------------------------------------------------------------------------!
 FUNCTION z_n_nuc_tayl( d1, dx, zm_para, zjnuc_t, zeta, z_gr_tot )

    IMPLICIT NONE

    INTEGER(iwp) ::  i !< running index

    REAL(wp) ::  d1            !< lower diameter limit
    REAL(wp) ::  dx            !< upper diameter limit
    REAL(wp) ::  zjnuc_t       !< initial nucleation rate (1/s)
    REAL(wp) ::  zeta          !< ratio of CS/GR (m) (condensation sink / growth rate)
    REAL(wp) ::  term1         !<
    REAL(wp) ::  term2         !<
    REAL(wp) ::  term3         !<
    REAL(wp) ::  term4         !<
    REAL(wp) ::  term5         !<
    REAL(wp) ::  z_n_nuc_tayl  !< final nucleation rate (1/s)
    REAL(wp) ::  z_gr_tot      !< total growth rate (nm/h)
    REAL(wp) ::  zm_para       !< m parameter in Lehtinen et al. (2007), Eq. 6

    z_n_nuc_tayl = 0.0_wp

    DO  i = 0, 29
       IF ( i == 0  .OR.  i == 1 )  THEN
          term1 = 1.0_wp
       ELSE
          term1 = term1 * REAL( i, SELECTED_REAL_KIND(12,307) )
       END IF
       term2 = ( REAL( i, SELECTED_REAL_KIND(12,307) ) * ( zm_para + 1.0_wp ) + 1.0_wp ) * term1
       term3 = zeta**i
       term4 = term3 / term2
       term5 = REAL( i, SELECTED_REAL_KIND(12,307) ) * ( zm_para + 1.0_wp ) + 1.0_wp
       z_n_nuc_tayl = z_n_nuc_tayl + term4 * ( dx**term5 - d1**term5 )
    ENDDO
    z_n_nuc_tayl = z_n_nuc_tayl * zjnuc_t * EXP( -zeta * ( d1**( zm_para + 1 ) ) ) / z_gr_tot

 END FUNCTION z_n_nuc_tayl

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculates the condensation of water vapour on aerosol particles. Follows the 
!> analytical predictor method by Jacobson (2005).
!> For equations, see Jacobson (2005), Fundamentals of atmospheric modelling 
!> (2nd edition).
!------------------------------------------------------------------------------!
 SUBROUTINE gpparth2o( paero, ptemp, ppres, pcs, pcw, ptstep )

    IMPLICIT NONE

    INTEGER(iwp) ::  ib   !< loop index
    INTEGER(iwp) ::  nstr !<

    REAL(wp) ::  adt        !< internal timestep in this subroutine
    REAL(wp) ::  rhoair     !< air density (kg/m3)
    REAL(wp) ::  ttot       !< total time (s)
    REAL(wp) ::  zact       !< Water activity
    REAL(wp) ::  zaelwc1    !< Current aerosol water content (kg/m3)
    REAL(wp) ::  zaelwc2    !< New aerosol water content after equilibrium calculation (kg/m3)
    REAL(wp) ::  zbeta      !< Transitional correction factor
    REAL(wp) ::  zcwc       !< Current water vapour mole concentration in aerosols (mol/m3)
    REAL(wp) ::  zcwint     !< Current and new water vapour mole concentrations (mol/m3)
    REAL(wp) ::  zcwn       !< New water vapour mole concentration (mol/m3)
    REAL(wp) ::  zcwtot     !< Total water mole concentration (mol/m3)
    REAL(wp) ::  zdfh2o     !< molecular diffusion coefficient (cm2/s) for water
    REAL(wp) ::  zhlp1      !< intermediate variable to calculate the mass transfer coefficient
    REAL(wp) ::  zhlp2      !< intermediate variable to calculate the mass transfer coefficient
    REAL(wp) ::  zhlp3      !< intermediate variable to calculate the mass transfer coefficient
    REAL(wp) ::  zknud      !< Knudsen number
    REAL(wp) ::  zmfph2o    !< mean free path of H2O gas molecule
    REAL(wp) ::  zrh        !< relative humidity [0-1]
    REAL(wp) ::  zthcond    !< thermal conductivity of air (W/m/K)

    REAL(wp), DIMENSION(nbins_aerosol) ::  zcwcae     !< Current water mole concentrations
    REAL(wp), DIMENSION(nbins_aerosol) ::  zcwintae   !< Current and new aerosol water mole concentration
    REAL(wp), DIMENSION(nbins_aerosol) ::  zcwnae     !< New water mole concentration in aerosols
    REAL(wp), DIMENSION(nbins_aerosol) ::  zcwsurfae  !< Surface mole concentration
    REAL(wp), DIMENSION(nbins_aerosol) ::  zkelvin    !< Kelvin effect
    REAL(wp), DIMENSION(nbins_aerosol) ::  zmtae      !< Mass transfer coefficients
    REAL(wp), DIMENSION(nbins_aerosol) ::  zwsatae    !< Water saturation ratio above aerosols

    REAL(wp), INTENT(in) ::  ppres   !< Air pressure (Pa)
    REAL(wp), INTENT(in) ::  pcs     !< Water vapour saturation concentration (kg/m3)
    REAL(wp), INTENT(in) ::  ptemp   !< Ambient temperature (K)
    REAL(wp), INTENT(in) ::  ptstep  !< timestep (s)

    REAL(wp), INTENT(inout) ::  pcw  !< Water vapour concentration (kg/m3)

    TYPE(t_section), DIMENSION(nbins_aerosol), INTENT(inout) ::  paero  !< Aerosol properties
!
!-- Relative humidity [0-1]
    zrh = pcw / pcs
!
!-- Calculate the condensation only for 2a/2b aerosol bins
    nstr = start_subrange_2a
!
!-- Save the current aerosol water content, 8 in paero is H2O
    zaelwc1 = SUM( paero(start_subrange_1a:end_subrange_2b)%volc(8) ) * arhoh2o
!
!-- Equilibration:
    IF ( advect_particle_water )  THEN
       IF ( zrh < 0.98_wp  .OR.  .NOT. lscndh2oae )  THEN
          CALL equilibration( zrh, ptemp, paero, .TRUE. )
       ELSE
          CALL equilibration( zrh, ptemp, paero, .FALSE. )
       ENDIF
    ENDIF
!
!-- The new aerosol water content after equilibrium calculation
    zaelwc2 = SUM( paero(start_subrange_1a:end_subrange_2b)%volc(8) ) * arhoh2o
!
!-- New water vapour mixing ratio (kg/m3)
    pcw = pcw - ( zaelwc2 - zaelwc1 ) * ppres * amdair / ( argas * ptemp )
!
!-- Initialise variables
    zcwsurfae(:) = 0.0_wp
    zhlp1        = 0.0_wp
    zhlp2        = 0.0_wp
    zhlp3        = 0.0_wp
    zmtae(:)     = 0.0_wp
    zwsatae(:)   = 0.0_wp
!
!-- Air:
!-- Density (kg/m3)
    rhoair = amdair * ppres / ( argas * ptemp )
!
!-- Thermal conductivity of air
    zthcond = 0.023807_wp + 7.1128E-5_wp * ( ptemp - 273.16_wp )
!
!-- Water vapour:
!-- Molecular diffusion coefficient (cm2/s) (eq.16.17)
    zdfh2o = ( 5.0_wp / ( 16.0_wp * avo * rhoair * 1.0E-3_wp * 3.11E-8_wp**2 ) ) * SQRT( argas *   &
               1.0E+7_wp * ptemp * amdair * 1.0E+3_wp * ( amh2o + amdair ) * 1.0E+3_wp /           &
               ( pi * amh2o * 2.0E+3_wp ) )
    zdfh2o = zdfh2o * 1.0E-4   ! Unit change to m^2/s
!
!-- Mean free path (eq. 15.25 & 16.29)
    zmfph2o = 3.0_wp * zdfh2o * SQRT( pi * amh2o / ( 8.0_wp * argas * ptemp ) )
!
!-- Kelvin effect (eq. 16.33)
    zkelvin(:) = EXP( 4.0_wp * surfw0 * amh2o / ( argas * ptemp * arhoh2o * paero(:)%dwet) )

    DO  ib = 1, nbins_aerosol
       IF ( paero(ib)%numc > nclim  .AND.  zrh > 0.98_wp )  THEN
!
!--       Water activity
          zact = acth2o( paero(ib) )
!
!--       Saturation mole concentration over flat surface. Limit the super-
!--       saturation to max 1.01 for the mass transfer. Experimental!
          zcwsurfae(ib) = MAX( pcs, pcw / 1.01_wp ) * rhoair / amh2o
!
!--       Equilibrium saturation ratio
          zwsatae(ib) = zact * zkelvin(ib)
!
!--       Knudsen number (eq. 16.20)
          zknud = 2.0_wp * zmfph2o / paero(ib)%dwet
!
!--       Transitional correction factor (Fuks & Sutugin, 1971)
          zbeta = ( zknud + 1.0_wp ) / ( 0.377_wp * zknud + 1.0_wp + 4.0_wp /                      &
                  ( 3.0_wp * massacc(ib) ) * ( zknud + zknud**2 ) )
!
!--       Mass transfer of H2O: Eq. 16.64 but here D^eff =  zdfh2o * zbeta
          zhlp1 = paero(ib)%numc * 2.0_wp * pi * paero(ib)%dwet * zdfh2o * zbeta
!
!--       1st term on the left side of the denominator in eq. 16.55
          zhlp2 = amh2o * zdfh2o * alv * zwsatae(ib) * zcwsurfae(ib) / ( zthcond * ptemp )
!
!--       2nd term on the left side of the denominator in eq. 16.55
          zhlp3 = ( ( alv * amh2o ) / ( argas * ptemp ) ) - 1.0_wp
!
!--       Full eq. 16.64: Mass transfer coefficient (1/s)
          zmtae(ib) = zhlp1 / ( zhlp2 * zhlp3 + 1.0_wp )
       ENDIF
    ENDDO
!
!-- Current mole concentrations of water
    zcwc        = pcw * rhoair / amh2o   ! as vapour
    zcwcae(:)   = paero(:)%volc(8) * arhoh2o / amh2o   ! in aerosols
    zcwtot      = zcwc + SUM( zcwcae )   ! total water concentration
    zcwnae(:)   = 0.0_wp
    zcwintae(:) = zcwcae(:)
!
!-- Substepping loop
    zcwint = 0.0_wp
    ttot   = 0.0_wp
    DO  WHILE ( ttot < ptstep )
       adt = 2.0E-2_wp   ! internal timestep
!
!--    New vapour concentration: (eq. 16.71)
       zhlp1 = zcwc + adt * ( SUM( zmtae(nstr:nbins_aerosol) * zwsatae(nstr:nbins_aerosol) *       &
                                   zcwsurfae(nstr:nbins_aerosol) ) )   ! numerator
       zhlp2 = 1.0_wp + adt * ( SUM( zmtae(nstr:nbins_aerosol) ) )   ! denomin.
       zcwint = zhlp1 / zhlp2   ! new vapour concentration
       zcwint = MIN( zcwint, zcwtot )
       IF ( ANY( paero(:)%numc > nclim )  .AND. zrh > 0.98_wp )  THEN
          DO  ib = nstr, nbins_aerosol
             zcwintae(ib) = zcwcae(ib) + MIN( MAX( adt * zmtae(ib) * ( zcwint - zwsatae(ib) *      &
                                                   zcwsurfae(ib) ), -0.02_wp * zcwcae(ib) ),       &
                                            0.05_wp * zcwcae(ib) )
             zwsatae(ib) = acth2o( paero(ib), zcwintae(ib) ) * zkelvin(ib)
          ENDDO
       ENDIF
       zcwintae(nstr:nbins_aerosol) = MAX( zcwintae(nstr:nbins_aerosol), 0.0_wp )
!
!--    Update vapour concentration for consistency
       zcwint = zcwtot - SUM( zcwintae(1:nbins_aerosol) )
!
!--    Update "old" values for next cycle
       zcwcae = zcwintae

       ttot = ttot + adt

    ENDDO   ! ADT

    zcwn      = zcwint
    zcwnae(:) = zcwintae(:)
    pcw       = zcwn * amh2o / rhoair
    paero(:)%volc(8) = MAX( 0.0_wp, zcwnae(:) * amh2o / arhoh2o )

 END SUBROUTINE gpparth2o

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculates the activity coefficient of liquid water 
!------------------------------------------------------------------------------!
 REAL(wp) FUNCTION acth2o( ppart, pcw )

    IMPLICIT NONE

    REAL(wp) ::  zns  !< molar concentration of solutes (mol/m3)
    REAL(wp) ::  znw  !< molar concentration of water (mol/m3)

    REAL(wp), INTENT(in), OPTIONAL ::  pcw !< molar concentration of water (mol/m3)

    TYPE(t_section), INTENT(in) ::  ppart !< Aerosol properties of a bin

    zns = ( 3.0_wp * ( ppart%volc(1) * arhoh2so4 / amh2so4 ) + ( ppart%volc(2) * arhooc / amoc ) + &
            2.0_wp * ( ppart%volc(5) * arhoss / amss ) + ( ppart%volc(6) * arhohno3 / amhno3 ) +   &
            ( ppart%volc(7) * arhonh3 / amnh3 ) )

    IF ( PRESENT(pcw) ) THEN
       znw = pcw
    ELSE
       znw = ppart%volc(8) * arhoh2o / amh2o
    ENDIF
!
!-- Activity = partial pressure of water vapour / sat. vapour pressure of water over a liquid surface
!--          = molality * activity coefficient (Jacobson, 2005: eq. 17.20-21)
!-- Assume activity coefficient of 1 for water
    acth2o = MAX( 0.1_wp, znw / MAX( EPSILON( 1.0_wp ),( znw + zns ) ) )

 END FUNCTION acth2o

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculates the dissolutional growth of particles (i.e. gas transfers to a 
!> particle surface and dissolves in liquid water on the surface). Treated here 
!> as a non-equilibrium (time-dependent) process. Gases: HNO3 and NH3 
!> (Chapter 17.14 in Jacobson, 2005).
!
!> Called from subroutine condensation. 
!> Coded by:
!> Harri Kokkola (FMI)
!------------------------------------------------------------------------------! 
 SUBROUTINE gpparthno3( ppres, ptemp, paero, pghno3, pgnh3, pcw, pcs, pbeta, ptstep )

    IMPLICIT NONE

    INTEGER(iwp) ::  ib  !< loop index

    REAL(wp) ::  adt          !< timestep
    REAL(wp) ::  zc_nh3_c     !< Current NH3 gas concentration
    REAL(wp) ::  zc_nh3_int   !< Intermediate NH3 gas concentration
    REAL(wp) ::  zc_nh3_n     !< New NH3 gas concentration
    REAL(wp) ::  zc_nh3_tot   !< Total NH3 concentration
    REAL(wp) ::  zc_hno3_c    !< Current HNO3 gas concentration
    REAL(wp) ::  zc_hno3_int  !< Intermediate HNO3 gas concentration
    REAL(wp) ::  zc_hno3_n    !< New HNO3 gas concentration
    REAL(wp) ::  zc_hno3_tot  !< Total HNO3 concentration
    REAL(wp) ::  zdfvap       !< Diffusion coefficient for vapors
    REAL(wp) ::  zhlp1        !< intermediate variable
    REAL(wp) ::  zhlp2        !< intermediate variable
    REAL(wp) ::  zrh          !< relative humidity

    REAL(wp), INTENT(in) ::  ppres      !< ambient pressure (Pa)
    REAL(wp), INTENT(in) ::  pcs        !< water vapour saturation
                                        !< concentration (kg/m3)
    REAL(wp), INTENT(in) ::  ptemp      !< ambient temperature (K)
    REAL(wp), INTENT(in) ::  ptstep     !< time step (s)

    REAL(wp), INTENT(inout) ::  pghno3  !< nitric acid concentration (#/m3)
    REAL(wp), INTENT(inout) ::  pgnh3   !< ammonia conc. (#/m3)
    REAL(wp), INTENT(inout) ::  pcw     !< water vapour concentration (kg/m3)

    REAL(wp), DIMENSION(nbins_aerosol) ::  zac_hno3_ae     !< Activity coefficients for HNO3
    REAL(wp), DIMENSION(nbins_aerosol) ::  zac_hhso4_ae    !< Activity coefficients for HHSO4
    REAL(wp), DIMENSION(nbins_aerosol) ::  zac_nh3_ae      !< Activity coefficients for NH3
    REAL(wp), DIMENSION(nbins_aerosol) ::  zac_nh4hso2_ae  !< Activity coefficients for NH4HSO2
    REAL(wp), DIMENSION(nbins_aerosol) ::  zcg_hno3_eq_ae  !< Equilibrium gas concentration: HNO3
    REAL(wp), DIMENSION(nbins_aerosol) ::  zcg_nh3_eq_ae   !< Equilibrium gas concentration: NH3
    REAL(wp), DIMENSION(nbins_aerosol) ::  zc_hno3_int_ae  !< Intermediate HNO3 aerosol concentration
    REAL(wp), DIMENSION(nbins_aerosol) ::  zc_hno3_c_ae    !< Current HNO3 in aerosols
    REAL(wp), DIMENSION(nbins_aerosol) ::  zc_hno3_n_ae    !< New HNO3 in aerosols
    REAL(wp), DIMENSION(nbins_aerosol) ::  zc_nh3_int_ae   !< Intermediate NH3 aerosol concentration
    REAL(wp), DIMENSION(nbins_aerosol) ::  zc_nh3_c_ae     !< Current NH3 in aerosols
    REAL(wp), DIMENSION(nbins_aerosol) ::  zc_nh3_n_ae     !< New NH3 in aerosols
    REAL(wp), DIMENSION(nbins_aerosol) ::  zkel_hno3_ae    !< Kelvin effect for HNO3
    REAL(wp), DIMENSION(nbins_aerosol) ::  zkel_nh3_ae     !< Kelvin effects for NH3
    REAL(wp), DIMENSION(nbins_aerosol) ::  zmt_hno3_ae     !< Mass transfer coefficients for HNO3
    REAL(wp), DIMENSION(nbins_aerosol) ::  zmt_nh3_ae      !< Mass transfer coefficients for NH3
    REAL(wp), DIMENSION(nbins_aerosol) ::  zsat_hno3_ae    !< HNO3 saturation ratio over a surface
    REAL(wp), DIMENSION(nbins_aerosol) ::  zsat_nh3_ae     !< NH3 saturation ratio over a surface

    REAL(wp), DIMENSION(nbins_aerosol,maxspec) ::  zion_mols   !< Ion molalities from pdfite aerosols

    REAL(wp), DIMENSION(nbins_aerosol), INTENT(in) ::  pbeta !< transitional correction factor for

    TYPE(t_section), DIMENSION(nbins_aerosol), INTENT(inout) ::  paero !< Aerosol properties
!
!-- Initialise:
    adt            = ptstep
    zac_hhso4_ae   = 0.0_wp
    zac_nh3_ae     = 0.0_wp
    zac_nh4hso2_ae = 0.0_wp
    zac_hno3_ae    = 0.0_wp
    zcg_nh3_eq_ae  = 0.0_wp
    zcg_hno3_eq_ae = 0.0_wp
    zion_mols      = 0.0_wp
    zsat_nh3_ae    = 1.0_wp
    zsat_hno3_ae   = 1.0_wp
!
!-- Diffusion coefficient (m2/s)
    zdfvap = 5.1111E-10_wp * ptemp**1.75_wp * ( p_0 + 1325.0_wp ) / ppres
!
!-- Kelvin effects (Jacobson (2005), eq. 16.33)
    zkel_hno3_ae(1:nbins_aerosol) = EXP( 4.0_wp * surfw0 * amvhno3 /                               &
                                    ( abo * ptemp * paero(1:nbins_aerosol)%dwet ) )
    zkel_nh3_ae(1:nbins_aerosol) = EXP( 4.0_wp * surfw0 * amvnh3 /                                 &
                                   ( abo * ptemp * paero(1:nbins_aerosol)%dwet ) )
!
!-- Current vapour mole concentrations (mol/m3)
    zc_hno3_c = pghno3 / avo  ! HNO3
    zc_nh3_c = pgnh3 / avo   ! NH3
!
!-- Current particle mole concentrations (mol/m3)
    zc_hno3_c_ae(1:nbins_aerosol) = paero(1:nbins_aerosol)%volc(6) * arhohno3 / amhno3
    zc_nh3_c_ae(1:nbins_aerosol) = paero(1:nbins_aerosol)%volc(7) * arhonh3 / amnh3
!
!-- Total mole concentrations: gas and particle phase
    zc_hno3_tot = zc_hno3_c + SUM( zc_hno3_c_ae(1:nbins_aerosol) )
    zc_nh3_tot = zc_nh3_c + SUM( zc_nh3_c_ae(1:nbins_aerosol) )
!
!-- Relative humidity [0-1]
    zrh = pcw / pcs
!
!-- Mass transfer coefficients (Jacobson, Eq. 16.64)
    zmt_hno3_ae(:) = 2.0_wp * pi * paero(:)%dwet * zdfvap * paero(:)%numc * pbeta(:)
    zmt_nh3_ae(:)  = 2.0_wp * pi * paero(:)%dwet * zdfvap * paero(:)%numc * pbeta(:)

!
!-- Get the equilibrium concentrations above aerosols
    CALL nitrate_ammonium_equilibrium( zrh, ptemp, paero, zcg_hno3_eq_ae, zcg_nh3_eq_ae,           &
                                       zac_hno3_ae, zac_nh3_ae, zac_nh4hso2_ae, zac_hhso4_ae,      &
                                       zion_mols )
!
!-- Calculate NH3 and HNO3 saturation ratios for aerosols
    CALL nitrate_ammonium_saturation( ptemp, paero, zac_hno3_ae, zac_nh4hso2_ae, zac_hhso4_ae,     &
                                      zcg_hno3_eq_ae, zc_hno3_c_ae, zc_nh3_c_ae, zkel_hno3_ae,     &
                                      zkel_nh3_ae, zsat_hno3_ae, zsat_nh3_ae )
!
!-- Intermediate gas concentrations of HNO3 and NH3
    zhlp1 = SUM( zc_hno3_c_ae(:) / ( 1.0_wp + adt * zmt_hno3_ae(:) * zsat_hno3_ae(:) ) )
    zhlp2 = SUM( zmt_hno3_ae(:) / ( 1.0_wp + adt * zmt_hno3_ae(:) * zsat_hno3_ae(:) ) )
    zc_hno3_int = ( zc_hno3_tot - zhlp1 ) / ( 1.0_wp + adt * zhlp2 )

    zhlp1 = SUM( zc_nh3_c_ae(:) / ( 1.0_wp + adt * zmt_nh3_ae(:) * zsat_nh3_ae(:) ) )
    zhlp2 = SUM( zmt_nh3_ae(:) / ( 1.0_wp + adt * zmt_nh3_ae(:) * zsat_nh3_ae(:) ) )
    zc_nh3_int = ( zc_nh3_tot - zhlp1 )/( 1.0_wp + adt * zhlp2 )

    zc_hno3_int = MIN( zc_hno3_int, zc_hno3_tot )
    zc_nh3_int = MIN( zc_nh3_int, zc_nh3_tot )
!
!-- Calculate the new concentration on aerosol particles
    zc_hno3_int_ae = zc_hno3_c_ae
    zc_nh3_int_ae = zc_nh3_c_ae
    DO  ib = 1, nbins_aerosol
       zc_hno3_int_ae(ib) = ( zc_hno3_c_ae(ib) + adt * zmt_hno3_ae(ib) * zc_hno3_int ) /           &
                            ( 1.0_wp + adt * zmt_hno3_ae(ib) * zsat_hno3_ae(ib) )
       zc_nh3_int_ae(ib) = ( zc_nh3_c_ae(ib) + adt * zmt_nh3_ae(ib) * zc_nh3_int ) /               &
                           ( 1.0_wp + adt * zmt_nh3_ae(ib) * zsat_nh3_ae(ib) )
    ENDDO

    zc_hno3_int_ae(:) = MAX( zc_hno3_int_ae(:), 0.0_wp )
    zc_nh3_int_ae(:) = MAX( zc_nh3_int_ae(:), 0.0_wp )
!
!-- Final molar gas concentration and molar particle concentration of HNO3
    zc_hno3_n   = zc_hno3_int
    zc_hno3_n_ae = zc_hno3_int_ae
!
!-- Final molar gas concentration and molar particle concentration of NH3
    zc_nh3_n   = zc_nh3_int
    zc_nh3_n_ae = zc_nh3_int_ae
!
!-- Model timestep reached - update the gas concentrations
    pghno3 = zc_hno3_n * avo
    pgnh3  = zc_nh3_n * avo
!
!-- Update the particle concentrations
    DO  ib = start_subrange_1a, end_subrange_2b
       paero(ib)%volc(6) = zc_hno3_n_ae(ib) * amhno3 / arhohno3
       paero(ib)%volc(7) = zc_nh3_n_ae(ib) * amnh3 / arhonh3
    ENDDO

 END SUBROUTINE gpparthno3
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate the equilibrium concentrations above aerosols (reference?)
!------------------------------------------------------------------------------!
 SUBROUTINE nitrate_ammonium_equilibrium( prh, ptemp, ppart, pcg_hno3_eq, pcg_nh3_eq, pgamma_hno3, &
                                          pgamma_nh4, pgamma_nh4hso2, pgamma_hhso4, pmols )

    IMPLICIT NONE

    INTEGER(iwp) ::  ib  !< loop index: aerosol bins

    REAL(wp) ::  zhlp         !< intermediate variable
    REAL(wp) ::  zp_hcl       !< Equilibrium vapor pressures (Pa) of HCl
    REAL(wp) ::  zp_hno3      !< Equilibrium vapor pressures (Pa) of HNO3
    REAL(wp) ::  zp_nh3       !< Equilibrium vapor pressures (Pa) of NH3
    REAL(wp) ::  zwatertotal  !< Total water in particles (mol/m3)

    REAL(wp), INTENT(in) ::  prh    !< relative humidity
    REAL(wp), INTENT(in) ::  ptemp  !< ambient temperature (K)

    REAL(wp), DIMENSION(maxspec) ::  zgammas  !< Activity coefficients
    REAL(wp), DIMENSION(maxspec) ::  zions    !< molar concentration of ion (mol/m3)

    REAL(wp), DIMENSION(nbins_aerosol), INTENT(inout) ::  pcg_nh3_eq      !< equilibrium molar
                                                                          !< concentration: of NH3
    REAL(wp), DIMENSION(nbins_aerosol), INTENT(inout) ::  pcg_hno3_eq     !< of HNO3
    REAL(wp), DIMENSION(nbins_aerosol), INTENT(inout) ::  pgamma_hhso4    !< activity coeff. of HHSO4
    REAL(wp), DIMENSION(nbins_aerosol), INTENT(inout) ::  pgamma_nh4      !< activity coeff. of NH3
    REAL(wp), DIMENSION(nbins_aerosol), INTENT(inout) ::  pgamma_nh4hso2  !< activity coeff. of NH4HSO2
    REAL(wp), DIMENSION(nbins_aerosol), INTENT(inout) ::  pgamma_hno3     !< activity coeff. of HNO3

    REAL(wp), DIMENSION(nbins_aerosol,maxspec), INTENT(inout) ::  pmols  !< Ion molalities

    TYPE(t_section), DIMENSION(nbins_aerosol), INTENT(inout) ::  ppart  !< Aerosol properties

    zgammas     = 0.0_wp
    zhlp        = 0.0_wp
    zions       = 0.0_wp
    zp_hcl      = 0.0_wp
    zp_hno3     = 0.0_wp
    zp_nh3      = 0.0_wp
    zwatertotal = 0.0_wp

    DO  ib = 1, nbins_aerosol

       IF ( ppart(ib)%numc < nclim )  CYCLE
!
!--    Ion molar concentrations: 2*H2SO4 + CL + NO3 - Na - NH4
       zhlp = 2.0_wp * ppart(ib)%volc(1) * arhoh2so4 / amh2so4 + ppart(ib)%volc(5) * arhoss / amss &
              + ppart(ib)%volc(6) * arhohno3 / amhno3 - ppart(ib)%volc(5) * arhoss / amss -        &
              ppart(ib)%volc(7) * arhonh3 / amnh3

       zions(1) = zhlp                                   ! H+
       zions(2) = ppart(ib)%volc(7) * arhonh3 / amnh3     ! NH4+
       zions(3) = ppart(ib)%volc(5) * arhoss / amss       ! Na+
       zions(4) = ppart(ib)%volc(1) * arhoh2so4 / amh2so4 ! SO4(2-)
       zions(5) = 0.0_wp                                 ! HSO4-
       zions(6) = ppart(ib)%volc(6) * arhohno3 / amhno3   ! NO3-
       zions(7) = ppart(ib)%volc(5) * arhoss / amss       ! Cl-

       zwatertotal = ppart(ib)%volc(8) * arhoh2o / amh2o
       IF ( zwatertotal > 1.0E-30_wp )  THEN
          CALL inorganic_pdfite( prh, ptemp, zions, zwatertotal, zp_hno3, zp_hcl, zp_nh3, zgammas, &
                                 pmols(ib,:) )
       ENDIF
!
!--    Activity coefficients
       pgamma_hno3(ib)    = zgammas(1)  ! HNO3
       pgamma_nh4(ib)     = zgammas(3)  ! NH3
       pgamma_nh4hso2(ib) = zgammas(6)  ! NH4HSO2
       pgamma_hhso4(ib)   = zgammas(7)  ! HHSO4
!
!--    Equilibrium molar concentrations (mol/m3) from equlibrium pressures (Pa)
       pcg_hno3_eq(ib) = zp_hno3 / ( argas * ptemp )
       pcg_nh3_eq(ib) = zp_nh3 / ( argas * ptemp )

    ENDDO

  END SUBROUTINE nitrate_ammonium_equilibrium

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate saturation ratios of NH4 and HNO3 for aerosols
!------------------------------------------------------------------------------!
 SUBROUTINE nitrate_ammonium_saturation( ptemp, ppart, pachno3, pacnh4hso2, pachhso4, pchno3eq,    &
                                         pchno3, pc_nh3, pkelhno3, pkelnh3, psathno3, psatnh3 )

    IMPLICIT NONE

    INTEGER(iwp) :: ib   !< running index for aerosol bins

    REAL(wp) ::  k_ll_h2o   !< equilibrium constants of equilibrium reactions:
                            !< H2O(aq) <--> H+ + OH- (mol/kg)
    REAL(wp) ::  k_ll_nh3   !< NH3(aq) + H2O(aq) <--> NH4+ + OH- (mol/kg)
    REAL(wp) ::  k_gl_nh3   !< NH3(g) <--> NH3(aq) (mol/kg/atm)
    REAL(wp) ::  k_gl_hno3  !< HNO3(g) <--> H+ + NO3- (mol2/kg2/atm)
    REAL(wp) ::  zmol_no3   !< molality of NO3- (mol/kg)
    REAL(wp) ::  zmol_h     !< molality of H+ (mol/kg)
    REAL(wp) ::  zmol_so4   !< molality of SO4(2-) (mol/kg)
    REAL(wp) ::  zmol_cl    !< molality of Cl- (mol/kg)
    REAL(wp) ::  zmol_nh4   !< molality of NH4+ (mol/kg)
    REAL(wp) ::  zmol_na    !< molality of Na+ (mol/kg)
    REAL(wp) ::  zhlp1      !< intermediate variable
    REAL(wp) ::  zhlp2      !< intermediate variable
    REAL(wp) ::  zhlp3      !< intermediate variable
    REAL(wp) ::  zxi        !< particle mole concentration ratio: (NH3+SS)/H2SO4
    REAL(wp) ::  zt0        !< reference temp

    REAL(wp), PARAMETER ::  a1 = -22.52_wp     !<
    REAL(wp), PARAMETER ::  a2 = -1.50_wp      !<
    REAL(wp), PARAMETER ::  a3 = 13.79_wp      !<
    REAL(wp), PARAMETER ::  a4 = 29.17_wp      !<
    REAL(wp), PARAMETER ::  b1 = 26.92_wp      !<
    REAL(wp), PARAMETER ::  b2 = 26.92_wp      !<
    REAL(wp), PARAMETER ::  b3 = -5.39_wp      !<
    REAL(wp), PARAMETER ::  b4 = 16.84_wp      !<
    REAL(wp), PARAMETER ::  K01 = 1.01E-14_wp  !<
    REAL(wp), PARAMETER ::  K02 = 1.81E-5_wp   !<
    REAL(wp), PARAMETER ::  K03 = 57.64_wp     !<
    REAL(wp), PARAMETER ::  K04 = 2.51E+6_wp   !<

    REAL(wp), INTENT(in) ::  ptemp  !< ambient temperature (K)

    REAL(wp), DIMENSION(nbins_aerosol), INTENT(in) ::  pachhso4    !< activity coeff. of HHSO4
    REAL(wp), DIMENSION(nbins_aerosol), INTENT(in) ::  pacnh4hso2  !< activity coeff. of NH4HSO2
    REAL(wp), DIMENSION(nbins_aerosol), INTENT(in) ::  pachno3     !< activity coeff. of HNO3
    REAL(wp), DIMENSION(nbins_aerosol), INTENT(in) ::  pchno3eq    !< eq. surface concentration: HNO3
    REAL(wp), DIMENSION(nbins_aerosol), INTENT(in) ::  pchno3      !< current particle mole
                                                                   !< concentration of HNO3 (mol/m3)
    REAL(wp), DIMENSION(nbins_aerosol), INTENT(in) ::  pc_nh3      !< of NH3 (mol/m3)
    REAL(wp), DIMENSION(nbins_aerosol), INTENT(in) ::  pkelhno3    !< Kelvin effect for HNO3
    REAL(wp), DIMENSION(nbins_aerosol), INTENT(in) ::  pkelnh3     !< Kelvin effect for NH3

    REAL(wp), DIMENSION(nbins_aerosol), INTENT(out) ::  psathno3 !< saturation ratio of HNO3
    REAL(wp), DIMENSION(nbins_aerosol), INTENT(out) ::  psatnh3  !< saturation ratio of NH3

    TYPE(t_section), DIMENSION(nbins_aerosol), INTENT(inout) ::  ppart  !< Aerosol properties

    zmol_cl  = 0.0_wp
    zmol_h   = 0.0_wp
    zmol_na  = 0.0_wp
    zmol_nh4 = 0.0_wp
    zmol_no3 = 0.0_wp
    zmol_so4 = 0.0_wp
    zt0      = 298.15_wp
    zxi      = 0.0_wp
!
!-- Calculates equlibrium rate constants based on Table B.7 in Jacobson (2005):
!-- K^ll_H20, K^ll_NH3, K^gl_NH3, K^gl_HNO3
    zhlp1 = zt0 / ptemp
    zhlp2 = zhlp1 - 1.0_wp
    zhlp3 = 1.0_wp + LOG( zhlp1 ) - zhlp1

    k_ll_h2o  = K01 * EXP( a1 * zhlp2 + b1 * zhlp3 )
    k_ll_nh3  = K02 * EXP( a2 * zhlp2 + b2 * zhlp3 )
    k_gl_nh3  = K03 * EXP( a3 * zhlp2 + b3 * zhlp3 )
    k_gl_hno3 = K04 * EXP( a4 * zhlp2 + b4 * zhlp3 )

    DO  ib = 1, nbins_aerosol

       IF ( ppart(ib)%numc > nclim  .AND.  ppart(ib)%volc(8) > 1.0E-30_wp  )  THEN
!
!--       Molality of H+ and NO3-
          zhlp1 = pc_nh3(ib) * amnh3 + ppart(ib)%volc(1) * arhoh2so4 + ppart(ib)%volc(2) * arhooc  &
                  + ppart(ib)%volc(5) * arhoss + ppart(ib)%volc(8) * arhoh2o
          zmol_no3 = pchno3(ib) / zhlp1  !< mol/kg
!
!--       Particle mole concentration ratio: (NH3+SS)/H2SO4
          zxi = ( pc_nh3(ib) + ppart(ib)%volc(5) * arhoss / amss ) / ( ppart(ib)%volc(1) *         &
                  arhoh2so4 / amh2so4 )

          IF ( zxi <= 2.0_wp )  THEN
!
!--          Molality of SO4(2-)
             zhlp1 = pc_nh3(ib) * amnh3 + pchno3(ib) * amhno3 + ppart(ib)%volc(2) * arhooc +       &
                     ppart(ib)%volc(5) * arhoss + ppart(ib)%volc(8) * arhoh2o
             zmol_so4 = ( ppart(ib)%volc(1) * arhoh2so4 / amh2so4 ) / zhlp1
!
!--          Molality of Cl-
             zhlp1 = pc_nh3(ib) * amnh3 + pchno3(ib) * amhno3 + ppart(ib)%volc(2) * arhooc +       &
                     ppart(ib)%volc(1) * arhoh2so4 + ppart(ib)%volc(8) * arhoh2o
             zmol_cl = ( ppart(ib)%volc(5) * arhoss / amss ) / zhlp1
!
!--          Molality of NH4+
             zhlp1 =  pchno3(ib) * amhno3 + ppart(ib)%volc(1) * arhoh2so4 + ppart(ib)%volc(2) *    &
                      arhooc + ppart(ib)%volc(5) * arhoss + ppart(ib)%volc(8) * arhoh2o
             zmol_nh4 = pc_nh3(ib) / zhlp1
!
!--          Molality of Na+
             zmol_na = zmol_cl
!
!--          Molality of H+
             zmol_h = 2.0_wp * zmol_so4 + zmol_no3 + zmol_cl - ( zmol_nh4 + zmol_na )

          ELSE

             zhlp2 = pkelhno3(ib) * zmol_no3 * pachno3(ib)**2

             IF ( zhlp2 > 1.0E-30_wp )  THEN
                zmol_h = k_gl_hno3 * pchno3eq(ib) / zhlp2 ! Eq. 17.38
             ELSE
                zmol_h = 0.0_wp
             ENDIF

          ENDIF

          zhlp1 = ppart(ib)%volc(8) * arhoh2o * argas * ptemp * k_gl_hno3
!
!--       Saturation ratio for NH3 and for HNO3
          IF ( zmol_h > 0.0_wp )  THEN
             zhlp2 = pkelnh3(ib) / ( zhlp1 * zmol_h )
             zhlp3 = k_ll_h2o / ( k_ll_nh3 + k_gl_nh3 )
             psatnh3(ib) = zhlp2 * ( ( pacnh4hso2(ib) / pachhso4(ib) )**2 ) * zhlp3
             psathno3(ib) = ( pkelhno3(ib) * zmol_h * pachno3(ib)**2 ) / zhlp1
          ELSE
             psatnh3(ib) = 1.0_wp
             psathno3(ib) = 1.0_wp
          ENDIF
       ELSE
          psatnh3(ib) = 1.0_wp
          psathno3(ib) = 1.0_wp
       ENDIF

    ENDDO

  END SUBROUTINE nitrate_ammonium_saturation

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Prototype module for calculating the water content of a mixed inorganic/
!> organic particle + equilibrium water vapour pressure above the solution
!> (HNO3, HCL, NH3 and representative organic compounds. Efficient calculation
!> of the partitioning of species between gas and aerosol. Based in a chamber
!> study.
!
!> Written by Dave Topping. Pure organic component properties predicted by Mark
!> Barley based on VOCs predicted in MCM simulations performed by Mike Jenkin.
!> Delivered by Gordon McFiggans as Deliverable D22 from WP1.4 in the EU FP6
!> EUCAARI Integrated Project.
!
!> REFERENCES
!> Clegg et al. (1998) A Thermodynamic Model of the System H+-NH4+-Na+-SO42- -NO3--Cl--H2O at
!>    298.15 K, J. Phys. Chem., 102A, 2155-2171.
!> Clegg et al. (2001) Thermodynamic modelling of aqueous aerosols containing electrolytes and
!>    dissolved organic compounds. Journal of Aerosol Science 2001;32(6):713-738.
!> Topping et al. (2005a) A curved multi-component aerosol hygroscopicity model framework: Part 1 -
!>    Inorganic compounds. Atmospheric Chemistry and Physics 2005;5:1205-1222.
!> Topping et al. (2005b) A curved multi-component aerosol hygroscopicity model framework: Part 2 -
!>    Including organic compounds. Atmospheric Chemistry and Physics 2005;5:1223-1242.
!> Wagman et al. (1982). The NBS tables of chemical thermodynamic properties: selected values for
!>    inorganic and C₁ and C₂ organic substances in SI units (book)
!> Zaveri et al. (2005). A new method for multicomponent activity coefficients of electrolytes in
!>    aqueous atmospheric aerosols, JGR, 110, D02201, 2005.
!
!> Queries concerning the use of this code through Gordon McFiggans,
!> g.mcfiggans@manchester.ac.uk,
!> Ownership: D. Topping, Centre for Atmospheric Sciences, University of
!> Manchester, 2007
!
!> Rewritten to PALM by Mona Kurppa, UHel, 2017
!------------------------------------------------------------------------------!
 SUBROUTINE inorganic_pdfite( rh, temp, ions, water_total, press_hno3, press_hcl, press_nh3,       &
                              gamma_out, mols_out )

    IMPLICIT NONE

    INTEGER(iwp) ::  binary_case
    INTEGER(iwp) ::  full_complexity

    REAL(wp) ::  a                         !< auxiliary variable
    REAL(wp) ::  act_product               !< ionic activity coef. product:
                                           !< = (gamma_h2so4**3d0) / gamma_hhso4**2d0)
    REAL(wp) ::  ammonium_chloride         !<
    REAL(wp) ::  ammonium_chloride_eq_frac !<
    REAL(wp) ::  ammonium_nitrate          !<
    REAL(wp) ::  ammonium_nitrate_eq_frac  !<
    REAL(wp) ::  ammonium_sulphate         !<
    REAL(wp) ::  ammonium_sulphate_eq_frac !<
    REAL(wp) ::  b                         !< auxiliary variable
    REAL(wp) ::  binary_h2so4              !< binary H2SO4 activity coeff.
    REAL(wp) ::  binary_hcl                !< binary HCL activity coeff.
    REAL(wp) ::  binary_hhso4              !< binary HHSO4 activity coeff.
    REAL(wp) ::  binary_hno3               !< binary HNO3 activity coeff.
    REAL(wp) ::  binary_nh4hso4            !< binary NH4HSO4 activity coeff.
    REAL(wp) ::  c                         !< auxiliary variable
    REAL(wp) ::  charge_sum                !< sum of ionic charges
    REAL(wp) ::  gamma_h2so4               !< activity coefficient
    REAL(wp) ::  gamma_hcl                 !< activity coefficient
    REAL(wp) ::  gamma_hhso4               !< activity coeffient
    REAL(wp) ::  gamma_hno3                !< activity coefficient
    REAL(wp) ::  gamma_nh3                 !< activity coefficient
    REAL(wp) ::  gamma_nh4hso4             !< activity coefficient
    REAL(wp) ::  h_out                     !<
    REAL(wp) ::  h_real                    !< new hydrogen ion conc.
    REAL(wp) ::  h2so4_hcl                 !< contribution of H2SO4
    REAL(wp) ::  h2so4_hno3                !< contribution of H2SO4
    REAL(wp) ::  h2so4_nh3                 !< contribution of H2SO4
    REAL(wp) ::  h2so4_nh4hso4             !< contribution of H2SO4
    REAL(wp) ::  hcl_h2so4                 !< contribution of HCL
    REAL(wp) ::  hcl_hhso4                 !< contribution of HCL
    REAL(wp) ::  hcl_hno3                  !< contribution of HCL
    REAL(wp) ::  hcl_nh4hso4               !< contribution of HCL
    REAL(wp) ::  henrys_temp_dep           !< temperature dependence of Henry's Law
    REAL(wp) ::  hno3_h2so4                !< contribution of HNO3
    REAL(wp) ::  hno3_hcl                  !< contribution of HNO3
    REAL(wp) ::  hno3_hhso4                !< contribution of HNO3
    REAL(wp) ::  hno3_nh3                  !< contribution of HNO3
    REAL(wp) ::  hno3_nh4hso4              !< contribution of HNO3
    REAL(wp) ::  hso4_out                  !<
    REAL(wp) ::  hso4_real                 !< new bisulphate ion conc.
    REAL(wp) ::  hydrochloric_acid         !<
    REAL(wp) ::  hydrochloric_acid_eq_frac !<
    REAL(wp) ::  k_h                       !< equilibrium constant for H+
    REAL(wp) ::  k_hcl                     !< equilibrium constant of HCL
    REAL(wp) ::  k_hno3                    !< equilibrium constant of HNO3
    REAL(wp) ::  k_nh4                     !< equilibrium constant for NH4+
    REAL(wp) ::  k_h2o                     !< equil. const. for water_surface
    REAL(wp) ::  ln_h2so4_act              !< gamma_h2so4 = EXP(ln_h2so4_act)
    REAL(wp) ::  ln_HCL_act                !< gamma_hcl = EXP( ln_HCL_act )
    REAL(wp) ::  ln_hhso4_act              !< gamma_hhso4 = EXP(ln_hhso4_act)
    REAL(wp) ::  ln_hno3_act               !< gamma_hno3 = EXP( ln_hno3_act )
    REAL(wp) ::  ln_nh4hso4_act            !< gamma_nh4hso4 = EXP( ln_nh4hso4_act )
    REAL(wp) ::  molality_ratio_nh3        !< molality ratio of NH3 (NH4+ and H+)
    REAL(wp) ::  na2so4_h2so4              !< contribution of Na2SO4
    REAL(wp) ::  na2so4_hcl                !< contribution of Na2SO4
    REAL(wp) ::  na2so4_hhso4              !< contribution of Na2SO4
    REAL(wp) ::  na2so4_hno3               !< contribution of Na2SO4
    REAL(wp) ::  na2so4_nh3                !< contribution of Na2SO4
    REAL(wp) ::  na2so4_nh4hso4            !< contribution of Na2SO4
    REAL(wp) ::  nacl_h2so4                !< contribution of NaCl
    REAL(wp) ::  nacl_hcl                  !< contribution of NaCl
    REAL(wp) ::  nacl_hhso4                !< contribution of NaCl
    REAL(wp) ::  nacl_hno3                 !< contribution of NaCl
    REAL(wp) ::  nacl_nh3                  !< contribution of NaCl
    REAL(wp) ::  nacl_nh4hso4              !< contribution of NaCl
    REAL(wp) ::  nano3_h2so4               !< contribution of NaNO3
    REAL(wp) ::  nano3_hcl                 !< contribution of NaNO3
    REAL(wp) ::  nano3_hhso4               !< contribution of NaNO3
    REAL(wp) ::  nano3_hno3                !< contribution of NaNO3
    REAL(wp) ::  nano3_nh3                 !< contribution of NaNO3
    REAL(wp) ::  nano3_nh4hso4             !< contribution of NaNO3
    REAL(wp) ::  nh42so4_h2so4             !< contribution of NH42SO4
    REAL(wp) ::  nh42so4_hcl               !< contribution of NH42SO4
    REAL(wp) ::  nh42so4_hhso4             !< contribution of NH42SO4
    REAL(wp) ::  nh42so4_hno3              !< contribution of NH42SO4
    REAL(wp) ::  nh42so4_nh3               !< contribution of NH42SO4
    REAL(wp) ::  nh42so4_nh4hso4           !< contribution of NH42SO4
    REAL(wp) ::  nh4cl_h2so4               !< contribution of NH4Cl
    REAL(wp) ::  nh4cl_hcl                 !< contribution of NH4Cl
    REAL(wp) ::  nh4cl_hhso4               !< contribution of NH4Cl
    REAL(wp) ::  nh4cl_hno3                !< contribution of NH4Cl
    REAL(wp) ::  nh4cl_nh3                 !< contribution of NH4Cl
    REAL(wp) ::  nh4cl_nh4hso4             !< contribution of NH4Cl
    REAL(wp) ::  nh4no3_h2so4              !< contribution of NH4NO3
    REAL(wp) ::  nh4no3_hcl                !< contribution of NH4NO3
    REAL(wp) ::  nh4no3_hhso4              !< contribution of NH4NO3
    REAL(wp) ::  nh4no3_hno3               !< contribution of NH4NO3
    REAL(wp) ::  nh4no3_nh3                !< contribution of NH4NO3
    REAL(wp) ::  nh4no3_nh4hso4            !< contribution of NH4NO3
    REAL(wp) ::  nitric_acid               !<
    REAL(wp) ::  nitric_acid_eq_frac       !< Equivalent fractions
    REAL(wp) ::  press_hcl                 !< partial pressure of HCL
    REAL(wp) ::  press_hno3                !< partial pressure of HNO3
    REAL(wp) ::  press_nh3                 !< partial pressure of NH3
    REAL(wp) ::  rh                        !< relative humidity [0-1]
    REAL(wp) ::  root1                     !< auxiliary variable
    REAL(wp) ::  root2                     !< auxiliary variable
    REAL(wp) ::  so4_out                   !<
    REAL(wp) ::  so4_real                  !< new sulpate ion concentration
    REAL(wp) ::  sodium_chloride           !<
    REAL(wp) ::  sodium_chloride_eq_frac   !<
    REAL(wp) ::  sodium_nitrate            !<
    REAL(wp) ::  sodium_nitrate_eq_frac    !<
    REAL(wp) ::  sodium_sulphate           !<
    REAL(wp) ::  sodium_sulphate_eq_frac   !<
    REAL(wp) ::  solutes                   !<
    REAL(wp) ::  sulphuric_acid            !<
    REAL(wp) ::  sulphuric_acid_eq_frac    !<
    REAL(wp) ::  temp                      !< temperature
    REAL(wp) ::  water_total               !<

    REAL(wp), DIMENSION(:) ::  gamma_out !< Activity coefficient for calculating the non-ideal
                                         !< dissociation constants
                                         !< 1: HNO3, 2: HCL, 3: NH4+/H+ (NH3), 4: HHSO4**2/H2SO4,
                                         !< 5: H2SO4**3/HHSO4**2, 6: NH4HSO2, 7: HHSO4
    REAL(wp), DIMENSION(:) ::  ions      !< ion molarities (mol/m3): 1: H+, 2: NH4+, 3: Na+, 
                                         !< 4: SO4(2-), 5: HSO4-, 6: NO3-, 7: Cl-
    REAL(wp), DIMENSION(7) ::  ions_mol  !< ion molalities (mol/kg): 1: H+, 2: NH4+, 3: Na+, 
                                         !< 4: SO4(2-), 5: HSO4-, 6: NO3-, 7: Cl-
    REAL(wp), DIMENSION(:) ::  mols_out  !< ion molality output (mol/kg): 1: H+, 2: NH4+, 3: Na+, 
                                         !< 4: SO4(2-), 5: HSO4-, 6: NO3-, 7: Cl-
!
!-- Value initialisation
    binary_h2so4    = 0.0_wp
    binary_hcl      = 0.0_wp
    binary_hhso4    = 0.0_wp
    binary_hno3     = 0.0_wp
    binary_nh4hso4  = 0.0_wp
    henrys_temp_dep = ( 1.0_wp / temp - 0.0033557_wp ) ! 1/T - 1/298 K
    hcl_hno3        = 1.0_wp
    h2so4_hno3      = 1.0_wp
    nh42so4_hno3    = 1.0_wp
    nh4no3_hno3     = 1.0_wp
    nh4cl_hno3      = 1.0_wp
    na2so4_hno3     = 1.0_wp
    nano3_hno3      = 1.0_wp
    nacl_hno3       = 1.0_wp
    hno3_hcl        = 1.0_wp
    h2so4_hcl       = 1.0_wp
    nh42so4_hcl     = 1.0_wp
    nh4no3_hcl      = 1.0_wp
    nh4cl_hcl       = 1.0_wp
    na2so4_hcl      = 1.0_wp
    nano3_hcl       = 1.0_wp
    nacl_hcl        = 1.0_wp
    hno3_nh3        = 1.0_wp
    h2so4_nh3       = 1.0_wp
    nh42so4_nh3     = 1.0_wp
    nh4no3_nh3      = 1.0_wp
    nh4cl_nh3       = 1.0_wp
    na2so4_nh3      = 1.0_wp
    nano3_nh3       = 1.0_wp
    nacl_nh3        = 1.0_wp
    hno3_hhso4      = 1.0_wp
    hcl_hhso4       = 1.0_wp
    nh42so4_hhso4   = 1.0_wp
    nh4no3_hhso4    = 1.0_wp
    nh4cl_hhso4     = 1.0_wp
    na2so4_hhso4    = 1.0_wp
    nano3_hhso4     = 1.0_wp
    nacl_hhso4      = 1.0_wp
    hno3_h2so4      = 1.0_wp
    hcl_h2so4       = 1.0_wp
    nh42so4_h2so4   = 1.0_wp
    nh4no3_h2so4    = 1.0_wp
    nh4cl_h2so4     = 1.0_wp
    na2so4_h2so4    = 1.0_wp
    nano3_h2so4     = 1.0_wp
    nacl_h2so4      = 1.0_wp
!
!-- New NH3 variables
    hno3_nh4hso4    = 1.0_wp
    hcl_nh4hso4     = 1.0_wp
    h2so4_nh4hso4   = 1.0_wp
    nh42so4_nh4hso4 = 1.0_wp
    nh4no3_nh4hso4  = 1.0_wp
    nh4cl_nh4hso4   = 1.0_wp
    na2so4_nh4hso4  = 1.0_wp
    nano3_nh4hso4   = 1.0_wp
    nacl_nh4hso4    = 1.0_wp
!
!-- Juha Tonttila added
    mols_out   = 0.0_wp
    press_hno3 = 0.0_wp  !< Initialising vapour pressures over the
    press_hcl  = 0.0_wp  !< multicomponent particle
    press_nh3  = 0.0_wp
    gamma_out  = 1.0_wp  !< i.e. don't alter the ideal mixing ratios if there's nothing there.
!
!-- 1) - COMPOSITION DEFINITIONS
!
!-- a) Inorganic ion pairing:
!-- In order to calculate the water content, which is also used in calculating vapour pressures, one
!-- needs to pair the anions and cations for use in the ZSR mixing rule. The equation provided by
!-- Clegg et al. (2001) is used for ion pairing. The solutes chosen comprise of 9 inorganic salts
!-- and acids which provide a pairing between each anion and cation: (NH4)2SO4, NH4NO3, NH4Cl,
!-- Na2SO4, NaNO3, NaCl, H2SO4, HNO3, HCL. The organic compound is treated as a seperate solute.
!-- Ions: 1: H+, 2: NH4+, 3: Na+, 4: SO4(2-), 5: HSO4-, 6: NO3-, 7: Cl-
!
    charge_sum = ions(1) + ions(2) + ions(3) + 2.0_wp * ions(4) + ions(5) + ions(6) + ions(7)
    nitric_acid       = ( 2.0_wp * ions(1) * ions(6) ) / charge_sum
    hydrochloric_acid = ( 2.0_wp * ions(1) * ions(7) ) / charge_sum
    sulphuric_acid    = ( 2.0_wp * ions(1) * ions(4) ) / charge_sum
    ammonium_sulphate = ( 2.0_wp * ions(2) * ions(4) ) / charge_sum
    ammonium_nitrate  = ( 2.0_wp * ions(2) * ions(6) ) / charge_sum
    ammonium_chloride = ( 2.0_wp * ions(2) * ions(7) ) / charge_sum
    sodium_sulphate   = ( 2.0_wp * ions(3) * ions(4) ) / charge_sum
    sodium_nitrate    = ( 2.0_wp * ions(3) * ions(6) ) / charge_sum
    sodium_chloride   = ( 2.0_wp * ions(3) * ions(7) ) / charge_sum
    solutes = 0.0_wp
    solutes = 3.0_wp * sulphuric_acid    + 2.0_wp * hydrochloric_acid + 2.0_wp * nitric_acid +     &
              3.0_wp * ammonium_sulphate + 2.0_wp * ammonium_nitrate + 2.0_wp * ammonium_chloride +&
              3.0_wp * sodium_sulphate   + 2.0_wp * sodium_nitrate   + 2.0_wp * sodium_chloride
!
!-- b) Inorganic equivalent fractions:
!-- These values are calculated so that activity coefficients can be expressed by a linear additive
!-- rule, thus allowing more efficient calculations and future expansion (see more detailed
!-- description below)
    nitric_acid_eq_frac       = 2.0_wp * nitric_acid / solutes
    hydrochloric_acid_eq_frac = 2.0_wp * hydrochloric_acid / solutes
    sulphuric_acid_eq_frac    = 3.0_wp * sulphuric_acid / solutes
    ammonium_sulphate_eq_frac = 3.0_wp * ammonium_sulphate / solutes
    ammonium_nitrate_eq_frac  = 2.0_wp * ammonium_nitrate / solutes
    ammonium_chloride_eq_frac = 2.0_wp * ammonium_chloride / solutes
    sodium_sulphate_eq_frac   = 3.0_wp * sodium_sulphate / solutes
    sodium_nitrate_eq_frac    = 2.0_wp * sodium_nitrate / solutes
    sodium_chloride_eq_frac   = 2.0_wp * sodium_chloride / solutes
!
!-- Inorganic ion molalities
    ions_mol(1) = ions(1) / ( water_total * 18.01528E-3_wp )   ! H+
    ions_mol(2) = ions(2) / ( water_total * 18.01528E-3_wp )   ! NH4+
    ions_mol(3) = ions(3) / ( water_total * 18.01528E-3_wp )   ! Na+
    ions_mol(4) = ions(4) / ( water_total * 18.01528E-3_wp )   ! SO4(2-)
    ions_mol(5) = ions(5) / ( water_total * 18.01528E-3_wp )   ! HSO4(2-)
    ions_mol(6) = ions(6) / ( water_total * 18.01528E-3_wp )   !  NO3-
    ions_mol(7) = ions(7) / ( water_total * 18.01528E-3_wp )   ! Cl-

!-- ***
!-- At this point we may need to introduce a method for prescribing H+ when there is no 'real' value
!-- for H+..i.e. in the sulphate poor domain. This will give a value for solve quadratic proposed by
!-- Zaveri et al. 2005
!
!-- 2) - WATER CALCULATION
!
!-- a) The water content is calculated using the ZSR rule with solute concentrations calculated
!-- using 1a above. Whilst the usual approximation of ZSR relies on binary data consisting of 5th or
!-- higher order polynomials, in this code 4 different RH regimes are used, each housing cubic
!-- equations for the water associated with each solute listed above. Binary water contents for
!-- inorganic components were calculated using AIM online (Clegg et al 1998). The water associated
!-- with the organic compound is calculated assuming ideality and that aw = RH.
!
!-- b) Molality of each inorganic ion and organic solute (initial input) is calculated for use in
!-- vapour pressure calculation.
!
!-- 3) - BISULPHATE ION DISSOCIATION CALCULATION
!
!-- The dissociation of the bisulphate ion is calculated explicitly. A solution to the equilibrium
!-- equation between the bisulphate ion, hydrogen ion and sulphate ion is found using tabulated
!-- equilibrium constants (referenced). It is necessary to calculate the activity coefficients of
!-- HHSO4 and H2SO4 in a non-iterative manner. These are calculated using the same format as
!-- described in 4) below, where both activity coefficients were fit to the output from ADDEM
!-- (Topping et al 2005a,b) covering an extensive composition space, providing the activity
!-- coefficients and bisulphate ion dissociation as a function of equivalent mole fractions and
!-- relative humidity.
!
!-- NOTE: the flags "binary_case" and "full_complexity" are not used in this prototype. They are
!-- used for simplification of the fit expressions when using limited composition regions. This
!-- section of code calculates the bisulphate ion concentration.
!
    IF ( ions(1) > 0.0_wp .AND. ions(4) > 0.0_wp ) THEN
!
!--    HHSO4:
       binary_case = 1
       IF ( rh > 0.1_wp  .AND.  rh < 0.9_wp )  THEN
          binary_hhso4 = -4.9521_wp * rh**3 + 9.2881_wp * rh**2 - 10.777_wp * rh + 6.0534_wp
       ELSEIF ( rh >= 0.9_wp  .AND.  rh < 0.955_wp )  THEN
          binary_hhso4 = -6.3777_wp * rh + 5.962_wp
       ELSEIF ( rh >= 0.955_wp  .AND.  rh < 0.99_wp )  THEN
          binary_hhso4 = 2367.2_wp * rh**3 - 6849.7_wp * rh**2 + 6600.9_wp * rh - 2118.7_wp
       ELSEIF ( rh >= 0.99_wp  .AND.  rh < 0.9999_wp )  THEN
          binary_hhso4 = 3E-7_wp * rh**5 - 2E-5_wp * rh**4 + 0.0004_wp * rh**3 - 0.0035_wp * rh**2 &
                         + 0.0123_wp * rh - 0.3025_wp
       ENDIF

       IF ( nitric_acid > 0.0_wp )  THEN
          hno3_hhso4 = -4.2204_wp * rh**4 + 12.193_wp * rh**3 - 12.481_wp * rh**2 + 6.459_wp * rh  &
                       - 1.9004_wp
       ENDIF

       IF ( hydrochloric_acid > 0.0_wp )  THEN
          hcl_hhso4 = -54.845_wp * rh**7 + 209.54_wp * rh**6 - 336.59_wp * rh**5 + 294.21_wp *     &
                      rh**4 - 150.07_wp * rh**3 + 43.767_wp * rh**2 - 6.5495_wp * rh + 0.60048_wp
       ENDIF

       IF ( ammonium_sulphate > 0.0_wp )  THEN
          nh42so4_hhso4 = 16.768_wp * rh**3 - 28.75_wp * rh**2 + 20.011_wp * rh - 8.3206_wp
       ENDIF

       IF ( ammonium_nitrate > 0.0_wp )  THEN
          nh4no3_hhso4 = -17.184_wp * rh**4 + 56.834_wp * rh**3 - 65.765_wp * rh**2 +              &
                         35.321_wp * rh - 9.252_wp
       ENDIF

       IF (ammonium_chloride > 0.0_wp )  THEN
          IF ( rh < 0.2_wp .AND. rh >= 0.1_wp )  THEN
             nh4cl_hhso4 = 3.2809_wp * rh - 2.0637_wp
          ELSEIF ( rh >= 0.2_wp .AND. rh < 0.99_wp )  THEN
             nh4cl_hhso4 = -1.2981_wp * rh**3 + 4.7461_wp * rh**2 - 2.3269_wp * rh - 1.1259_wp
          ENDIF
       ENDIF

       IF ( sodium_sulphate > 0.0_wp )  THEN
          na2so4_hhso4 = 118.87_wp * rh**6 - 358.63_wp * rh**5 + 435.85_wp * rh**4 - 272.88_wp *   &
                         rh**3 + 94.411_wp * rh**2 - 18.21_wp * rh + 0.45935_wp
       ENDIF

       IF ( sodium_nitrate > 0.0_wp )  THEN
          IF ( rh < 0.2_wp  .AND.  rh >= 0.1_wp )  THEN
             nano3_hhso4 = 4.8456_wp * rh - 2.5773_wp
          ELSEIF ( rh >= 0.2_wp  .AND.  rh < 0.99_wp )  THEN
             nano3_hhso4 = 0.5964_wp * rh**3 - 0.38967_wp * rh**2 + 1.7918_wp * rh - 1.9691_wp
          ENDIF
       ENDIF

       IF ( sodium_chloride > 0.0_wp )  THEN
          IF ( rh < 0.2_wp )  THEN
             nacl_hhso4 = 0.51995_wp * rh - 1.3981_wp
          ELSEIF ( rh >= 0.2_wp  .AND.  rh < 0.99_wp )  THEN
             nacl_hhso4 = 1.6539_wp * rh - 1.6101_wp
          ENDIF
       ENDIF

       ln_hhso4_act = binary_hhso4 + nitric_acid_eq_frac * hno3_hhso4 +                            &
                      hydrochloric_acid_eq_frac * hcl_hhso4 +                                      &
                      ammonium_sulphate_eq_frac * nh42so4_hhso4 +                                  &
                      ammonium_nitrate_eq_frac  * nh4no3_hhso4 +                                   &
                      ammonium_chloride_eq_frac * nh4cl_hhso4 +                                    &
                      sodium_sulphate_eq_frac   * na2so4_hhso4 +                                   &
                      sodium_nitrate_eq_frac * nano3_hhso4 + sodium_chloride_eq_frac   * nacl_hhso4

       gamma_hhso4 = EXP( ln_hhso4_act )   ! molal activity coefficient of HHSO4

!--    H2SO4 (sulphuric acid):
       IF ( rh >= 0.1_wp  .AND.  rh < 0.9_wp )  THEN
          binary_h2so4 = 2.4493_wp * rh**2 - 6.2326_wp * rh + 2.1763_wp
       ELSEIF ( rh >= 0.9_wp  .AND.  rh < 0.98 )  THEN
          binary_h2so4 = 914.68_wp * rh**3 - 2502.3_wp * rh**2 + 2281.9_wp * rh - 695.11_wp
       ELSEIF ( rh >= 0.98  .AND.  rh < 0.9999 )  THEN
          binary_h2so4 = 3.0E-8_wp * rh**4 - 5E-6_wp * rh**3 + 0.0003_wp * rh**2 - 0.0022_wp *     &
                         rh - 1.1305_wp
       ENDIF

       IF ( nitric_acid > 0.0_wp )  THEN
          hno3_h2so4 = - 16.382_wp * rh**5 + 46.677_wp * rh**4 - 54.149_wp * rh**3 + 34.36_wp *    &
                         rh**2 - 12.54_wp * rh + 2.1368_wp
       ENDIF

       IF ( hydrochloric_acid > 0.0_wp )  THEN
          hcl_h2so4 = - 14.409_wp * rh**5 + 42.804_wp * rh**4 - 47.24_wp * rh**3 + 24.668_wp *     &
                        rh**2 - 5.8015_wp * rh + 0.084627_wp
       ENDIF

       IF ( ammonium_sulphate > 0.0_wp )  THEN
          nh42so4_h2so4 = 66.71_wp * rh**5 - 187.5_wp * rh**4 + 210.57_wp * rh**3 - 121.04_wp *    &
                          rh**2 + 39.182_wp * rh - 8.0606_wp
       ENDIF

       IF ( ammonium_nitrate > 0.0_wp )  THEN
          nh4no3_h2so4 = - 22.532_wp * rh**4 + 66.615_wp * rh**3 - 74.647_wp * rh**2 + 37.638_wp * &
                         rh - 6.9711_wp
       ENDIF

       IF ( ammonium_chloride > 0.0_wp )  THEN
          IF ( rh >= 0.1_wp  .AND.  rh < 0.2_wp )  THEN
             nh4cl_h2so4 = - 0.32089_wp * rh + 0.57738_wp
          ELSEIF ( rh >= 0.2_wp  .AND.  rh < 0.9_wp )  THEN
             nh4cl_h2so4 = 18.089_wp * rh**5 - 51.083_wp * rh**4 + 50.32_wp * rh**3 - 17.012_wp *  &
                           rh**2 - 0.93435_wp * rh + 1.0548_wp
          ELSEIF ( rh >= 0.9_wp  .AND.  rh < 0.99_wp )  THEN
             nh4cl_h2so4 = - 1.5749_wp * rh + 1.7002_wp
          ENDIF
       ENDIF

       IF ( sodium_sulphate > 0.0_wp )  THEN
          na2so4_h2so4 = 29.843_wp * rh**4 - 69.417_wp * rh**3 + 61.507_wp * rh**2 - 29.874_wp *   &
                         rh + 7.7556_wp
       ENDIF

       IF ( sodium_nitrate > 0.0_wp )  THEN
          nano3_h2so4 = - 122.37_wp * rh**6 + 427.43_wp * rh**5 - 604.68_wp * rh**4 + 443.08_wp *  &
                        rh**3 - 178.61_wp * rh**2 + 37.242_wp * rh - 1.9564_wp
       ENDIF

       IF ( sodium_chloride > 0.0_wp )  THEN
          nacl_h2so4 = - 40.288_wp * rh**5 + 115.61_wp * rh**4 - 129.99_wp * rh**3 + 72.652_wp *   &
                       rh**2 - 22.124_wp * rh + 4.2676_wp
       ENDIF

       ln_h2so4_act = binary_h2so4 + nitric_acid_eq_frac * hno3_h2so4 +                            &
                      hydrochloric_acid_eq_frac * hcl_h2so4 +                                      &
                      ammonium_sulphate_eq_frac * nh42so4_h2so4 +                                  &
                      ammonium_nitrate_eq_frac  * nh4no3_h2so4 +                                   &
                      ammonium_chloride_eq_frac * nh4cl_h2so4 +                                    &
                      sodium_sulphate_eq_frac * na2so4_h2so4 +                                     &
                      sodium_nitrate_eq_frac * nano3_h2so4 + sodium_chloride_eq_frac * nacl_h2so4

       gamma_h2so4 = EXP( ln_h2so4_act )    ! molal activity coefficient 
!
!--    Export activity coefficients
       IF ( gamma_h2so4 > 1.0E-10_wp )  THEN
          gamma_out(4) = gamma_hhso4**2 / gamma_h2so4
       ENDIF
       IF ( gamma_hhso4 > 1.0E-10_wp )  THEN
          gamma_out(5) = gamma_h2so4**3 / gamma_hhso4**2
       ENDIF
!
!--    Ionic activity coefficient product
       act_product = gamma_h2so4**3 / gamma_hhso4**2
!
!--    Solve the quadratic equation (i.e. x in ax**2 + bx + c = 0)
       a = 1.0_wp
       b = -1.0_wp * ( ions(4) + ions(1) + ( ( water_total * 18.0E-3_wp ) /                        &
           ( 99.0_wp * act_product ) ) )
       c = ions(4) * ions(1)
       root1 = ( ( -1.0_wp * b ) + ( ( ( b**2 ) - 4.0_wp * a * c )**0.5_wp ) ) / ( 2.0_wp * a )
       root2 = ( ( -1.0_wp * b ) - ( ( ( b**2 ) - 4.0_wp * a * c) **0.5_wp ) ) / ( 2.0_wp * a )

       IF ( root1 > ions(1)  .OR.  root1 < 0.0_wp )  THEN
          root1 = 0.0_wp
       ENDIF

       IF ( root2 > ions(1)  .OR.  root2 < 0.0_wp )  THEN
          root2 = 0.0_wp
       ENDIF
!
!--    Calculate the new hydrogen ion, bisulphate ion and sulphate ion 
!--    concentration
       h_real    = ions(1)
       so4_real  = ions(4)
       hso4_real = MAX( root1, root2 )
       h_real   = ions(1) - hso4_real
       so4_real = ions(4) - hso4_real
!
!--    Recalculate ion molalities
       ions_mol(1) = h_real    / ( water_total * 18.01528E-3_wp )   ! H+
       ions_mol(4) = so4_real  / ( water_total * 18.01528E-3_wp )   ! SO4(2-)
       ions_mol(5) = hso4_real / ( water_total * 18.01528E-3_wp )   ! HSO4(2-)

       h_out    = h_real
       hso4_out = hso4_real
       so4_out  = so4_real

    ELSE
       h_out    = ions(1)
       hso4_out = 0.0_wp
       so4_out  = ions(4)
    ENDIF

!
!-- 4) ACTIVITY COEFFICIENTS -for vapour pressures of HNO3,HCL and NH3
!
!-- This section evaluates activity coefficients and vapour pressures using the water content 
!-- calculated above) for each inorganic condensing species: a - HNO3, b - NH3, c - HCL.
!-- The following procedure is used: Zaveri et al (2005) found that one could express the variation 
!-- of activity coefficients linearly in log-space if equivalent mole fractions were used.
!-- So, by a taylor series expansion LOG( activity coefficient ) =
!--    LOG( binary activity coefficient at a given RH ) +
!--    (equivalent mole fraction compound A) *
!--    ('interaction' parameter between A and condensing species) +
!--    equivalent mole fraction compound B) *
!--    ('interaction' parameter between B and condensing species).
!-- Here, the interaction parameters have been fit to ADDEM by searching the whole compositon space
!-- and fit usign the Levenberg-Marquardt non-linear least squares algorithm.
!
!-- They are given as a function of RH and vary with complexity ranging from linear to 5th order
!-- polynomial expressions, the binary activity coefficients were calculated using AIM online.
!-- NOTE: for NH3, no binary activity coefficient was used and the data were fit to the ratio of the
!-- activity coefficients for the ammonium and hydrogen ions. Once the activity coefficients are
!-- obtained the vapour pressure can be easily calculated using tabulated equilibrium constants
!-- (referenced). This procedure differs from that of Zaveri et al (2005) in that it is not assumed
!-- one can carry behaviour from binary mixtures in multicomponent systems. To this end we have fit
!-- the 'interaction' parameters explicitly to a general inorganic equilibrium model
!-- (ADDEM - Topping et al. 2005a,b). Such parameters take into account bisulphate ion dissociation
!-- and water content. This also allows us to consider one regime for all composition space, rather
!-- than defining sulphate rich and sulphate poor regimes.
!-- NOTE: The flags "binary_case" and "full_complexity" are not used in this prototype. They are
!-- used for simplification of the fit expressions when using limited composition regions.
!
!-- a) - ACTIVITY COEFF/VAPOUR PRESSURE - HNO3
    IF ( ions(1) > 0.0_wp  .AND.  ions(6) > 0.0_wp )  THEN
       binary_case = 1
       IF ( rh > 0.1_wp  .AND.  rh < 0.98_wp )  THEN
          IF ( binary_case == 1 )  THEN
             binary_hno3 = 1.8514_wp * rh**3 - 4.6991_wp * rh**2 + 1.5514_wp * rh + 0.90236_wp
          ELSEIF ( binary_case == 2 )  THEN
             binary_hno3 = - 1.1751_wp * ( rh**2 ) - 0.53794_wp * rh + 1.2808_wp
          ENDIF
       ELSEIF ( rh >= 0.98_wp  .AND.  rh < 0.9999_wp )  THEN
          binary_hno3 = 1244.69635941351_wp * rh**3 - 2613.93941099991_wp * rh**2 +                &
                        1525.0684974546_wp * rh -155.946764059316_wp
       ENDIF
!
!--    Contributions from other solutes
       full_complexity = 1
       IF ( hydrochloric_acid > 0.0_wp )  THEN   ! HCL
          IF ( full_complexity == 1  .OR.  rh < 0.4_wp )  THEN
             hcl_hno3 = 16.051_wp * rh**4 - 44.357_wp * rh**3 + 45.141_wp * rh**2 - 21.638_wp *    &
                        rh + 4.8182_wp
          ELSEIF ( full_complexity == 0  .AND.  rh > 0.4_wp )  THEN
             hcl_hno3 = - 1.5833_wp * rh + 1.5569_wp
          ENDIF
       ENDIF

       IF ( sulphuric_acid > 0.0_wp )  THEN   ! H2SO4
          IF ( full_complexity == 1  .OR.  rh < 0.4_wp )  THEN
             h2so4_hno3 = - 3.0849_wp * rh**3 + 5.9609_wp * rh**2 - 4.468_wp * rh + 1.5658_wp
          ELSEIF ( full_complexity == 0  .AND.  rh > 0.4_wp )  THEN
             h2so4_hno3 = - 0.93473_wp * rh + 0.9363_wp
          ENDIF
       ENDIF

       IF ( ammonium_sulphate > 0.0_wp )  THEN   ! NH42SO4
          nh42so4_hno3 = 16.821_wp * rh**3 - 28.391_wp * rh**2 + 18.133_wp * rh - 6.7356_wp
       ENDIF

       IF ( ammonium_nitrate > 0.0_wp )  THEN   ! NH4NO3
          nh4no3_hno3 = 11.01_wp * rh**3 - 21.578_wp * rh**2 + 14.808_wp * rh - 4.2593_wp
       ENDIF

       IF ( ammonium_chloride > 0.0_wp )  THEN   ! NH4Cl
          IF ( full_complexity == 1  .OR.  rh <= 0.4_wp )  THEN
             nh4cl_hno3 = - 1.176_wp * rh**3 + 5.0828_wp * rh**2 - 3.8792_wp * rh - 0.05518_wp
          ELSEIF ( full_complexity == 0  .AND.  rh > 0.4_wp )  THEN
             nh4cl_hno3 = 2.6219_wp * rh**2 - 2.2609_wp * rh - 0.38436_wp
          ENDIF
       ENDIF

       IF ( sodium_sulphate > 0.0_wp )  THEN   ! Na2SO4
          na2so4_hno3 = 35.504_wp * rh**4 - 80.101_wp * rh**3 + 67.326_wp * rh**2 - 28.461_wp *    &
                        rh + 5.6016_wp
       ENDIF

       IF ( sodium_nitrate > 0.0_wp )  THEN   ! NaNO3
          IF ( full_complexity == 1 .OR. rh <= 0.4_wp ) THEN
             nano3_hno3 = 23.659_wp * rh**5 - 66.917_wp * rh**4 + 74.686_wp * rh**3 - 40.795_wp *  &
                          rh**2 + 10.831_wp * rh - 1.4701_wp
          ELSEIF ( full_complexity == 0  .AND.  rh > 0.4_wp )  THEN
             nano3_hno3 = 14.749_wp * rh**4 - 35.237_wp * rh**3 + 31.196_wp * rh**2 - 12.076_wp *  &
                          rh + 1.3605_wp
          ENDIF
       ENDIF

       IF ( sodium_chloride > 0.0_wp )  THEN   ! NaCl
          IF ( full_complexity == 1  .OR.  rh <= 0.4_wp )  THEN
             nacl_hno3 = 13.682_wp * rh**4 - 35.122_wp * rh**3 + 33.397_wp * rh**2 - 14.586_wp *   &
                         rh + 2.6276_wp
          ELSEIF ( full_complexity == 0  .AND.  rh > 0.4_wp )  THEN
             nacl_hno3 = 1.1882_wp * rh**3 - 1.1037_wp * rh**2 - 0.7642_wp * rh + 0.6671_wp
          ENDIF
       ENDIF

       ln_hno3_act = binary_hno3 + hydrochloric_acid_eq_frac * hcl_hno3 +                          &
                     sulphuric_acid_eq_frac    * h2so4_hno3 +                                      &
                     ammonium_sulphate_eq_frac * nh42so4_hno3 +                                    &
                     ammonium_nitrate_eq_frac  * nh4no3_hno3 +                                     &
                     ammonium_chloride_eq_frac * nh4cl_hno3 +                                      &
                     sodium_sulphate_eq_frac * na2so4_hno3 +                                       &
                     sodium_nitrate_eq_frac * nano3_hno3 + sodium_chloride_eq_frac   * nacl_hno3

       gamma_hno3   = EXP( ln_hno3_act )   ! Molal activity coefficient of HNO3
       gamma_out(1) = gamma_hno3
!
!--    Partial pressure calculation
!--    k_hno3 = 2.51 * ( 10**6 )
!--    k_hno3 = 2.628145923d6 !< calculated by AIM online (Clegg et al 1998) after Chameides (1984)
       k_hno3     = 2.6E6_wp * EXP( 8700.0_wp * henrys_temp_dep )
       press_hno3 = ( ions_mol(1) * ions_mol(6) * ( gamma_hno3**2 ) ) / k_hno3
    ENDIF
!
!-- b) - ACTIVITY COEFF/VAPOUR PRESSURE - NH3
!-- Follow the two solute approach of Zaveri et al. (2005)
    IF ( ions(2) > 0.0_wp  .AND.  ions_mol(1) > 0.0_wp )  THEN
!
!--    NH4HSO4:
       binary_nh4hso4 = 56.907_wp * rh**6 - 155.32_wp * rh**5 + 142.94_wp * rh**4 - 32.298_wp *    &
                        rh**3 - 27.936_wp * rh**2 + 19.502_wp * rh - 4.2618_wp
       IF ( nitric_acid > 0.0_wp)  THEN   ! HNO3
          hno3_nh4hso4 = 104.8369_wp * rh**8 - 288.8923_wp * rh**7 + 129.3445_wp * rh**6 +         &
                         373.0471_wp * rh**5 - 571.0385_wp * rh**4 + 326.3528_wp * rh**3 -         &
                         74.169_wp * rh**2 - 2.4999_wp * rh + 3.17_wp
       ENDIF

       IF ( hydrochloric_acid > 0.0_wp)  THEN   ! HCL
          hcl_nh4hso4 = - 7.9133_wp * rh**8 + 126.6648_wp * rh**7 - 460.7425_wp * rh**6 +          &
                         731.606_wp * rh**5 - 582.7467_wp * rh**4 + 216.7197_wp * rh**3 -          &
                         11.3934_wp * rh**2 - 17.7728_wp  * rh + 5.75_wp
       ENDIF

       IF ( sulphuric_acid > 0.0_wp)  THEN   ! H2SO4
          h2so4_nh4hso4 = 195.981_wp * rh**8 - 779.2067_wp * rh**7 + 1226.3647_wp * rh**6 -        &
                         964.0261_wp * rh**5 + 391.7911_wp * rh**4 - 84.1409_wp  * rh**3 +         &
                          20.0602_wp * rh**2 - 10.2663_wp  * rh + 3.5817_wp
       ENDIF

       IF ( ammonium_sulphate > 0.0_wp)  THEN   ! NH42SO4
          nh42so4_nh4hso4 = 617.777_wp * rh**8 -  2547.427_wp * rh**7 + 4361.6009_wp * rh**6 -     &
                           4003.162_wp * rh**5 + 2117.8281_wp * rh**4 - 640.0678_wp * rh**3 +      &
                            98.0902_wp * rh**2 -    2.2615_wp * rh - 2.3811_wp
       ENDIF

       IF ( ammonium_nitrate > 0.0_wp)  THEN   ! NH4NO3
          nh4no3_nh4hso4 = - 104.4504_wp * rh**8 + 539.5921_wp * rh**7 - 1157.0498_wp * rh**6 +    &
                            1322.4507_wp * rh**5 - 852.2475_wp * rh**4 + 298.3734_wp * rh**3 -     &
                              47.0309_wp * rh**2 +    1.297_wp * rh - 0.8029_wp
       ENDIF

       IF ( ammonium_chloride > 0.0_wp)  THEN   ! NH4Cl
          nh4cl_nh4hso4 = 258.1792_wp * rh**8 - 1019.3777_wp * rh**7 + 1592.8918_wp * rh**6 -      &
                         1221.0726_wp * rh**5 +  442.2548_wp * rh**4 -   43.6278_wp * rh**3 -      &
                            7.5282_wp * rh**2 -    3.8459_wp * rh + 2.2728_wp
       ENDIF

       IF ( sodium_sulphate > 0.0_wp)  THEN   ! Na2SO4
          na2so4_nh4hso4 = 225.4238_wp * rh**8 - 732.4113_wp * rh**7 + 843.7291_wp * rh**6 -       &
                           322.7328_wp * rh**5 -  88.6252_wp * rh**4 +  72.4434_wp * rh**3 +       &
                            22.9252_wp * rh**2 -  25.3954_wp * rh + 4.6971_wp
       ENDIF

       IF ( sodium_nitrate > 0.0_wp)  THEN   ! NaNO3
          nano3_nh4hso4 = 96.1348_wp * rh**8 - 341.6738_wp * rh**7 + 406.5314_wp * rh**6 -         &
                          98.5777_wp * rh**5 - 172.8286_wp * rh**4 + 149.3151_wp * rh**3 -         &
                          38.9998_wp * rh**2 -   0.2251_wp * rh + 0.4953_wp
       ENDIF

       IF ( sodium_chloride > 0.0_wp)  THEN   ! NaCl
          nacl_nh4hso4 = 91.7856_wp * rh**8 - 316.6773_wp * rh**7 + 358.2703_wp * rh**6 -          &
                         68.9142_wp * rh**5 - 156.5031_wp * rh**4 + 116.9592_wp * rh**3 -          &
                         22.5271_wp * rh**2 - 3.7716_wp * rh + 1.56_wp
       ENDIF

       ln_nh4hso4_act = binary_nh4hso4 + nitric_acid_eq_frac * hno3_nh4hso4 +                      &
                        hydrochloric_acid_eq_frac * hcl_nh4hso4 +                                  &
                        sulphuric_acid_eq_frac * h2so4_nh4hso4 +                                   &
                        ammonium_sulphate_eq_frac * nh42so4_nh4hso4 +                              &
                        ammonium_nitrate_eq_frac * nh4no3_nh4hso4 +                                &
                        ammonium_chloride_eq_frac * nh4cl_nh4hso4 +                                &
                        sodium_sulphate_eq_frac * na2so4_nh4hso4 +                                 &
                        sodium_nitrate_eq_frac * nano3_nh4hso4 +                                   &
                        sodium_chloride_eq_frac * nacl_nh4hso4

       gamma_nh4hso4 = EXP( ln_nh4hso4_act ) ! molal act. coefficient of NH4HSO4
!
!--    Molal activity coefficient of NO3-
       gamma_out(6)  = gamma_nh4hso4
!
!--    Molal activity coefficient of NH4+
       gamma_nh3     = gamma_nh4hso4**2 / gamma_hhso4**2
       gamma_out(3)  = gamma_nh3
!
!--    This actually represents the ratio of the ammonium to hydrogen ion activity coefficients
!--    (see Zaveri paper) - multiply this by the ratio of the ammonium to hydrogen ion molality and
!--    the ratio of appropriate equilibrium constants
!
!--    Equilibrium constants
!--    k_h = 57.64d0    ! Zaveri et al. (2005)
       k_h = 5.8E1_wp * EXP( 4085.0_wp * henrys_temp_dep )   ! after Chameides (1984)
!--    k_nh4 = 1.81E-5_wp    ! Zaveri et al. (2005)
       k_nh4 = 1.7E-5_wp * EXP( -4325.0_wp * henrys_temp_dep )   ! Chameides (1984)
!--    k_h2o = 1.01E-14_wp    ! Zaveri et al (2005)
       k_h2o = 1.E-14_wp * EXP( -6716.0_wp * henrys_temp_dep )   ! Chameides (1984)
!
       molality_ratio_nh3 = ions_mol(2) / ions_mol(1)
!
!--    Partial pressure calculation
       press_nh3 = molality_ratio_nh3 * gamma_nh3 * ( k_h2o / ( k_h * k_nh4 ) )

    ENDIF
!
!-- c) - ACTIVITY COEFF/VAPOUR PRESSURE - HCL
    IF ( ions(1) > 0.0_wp  .AND.  ions(7) > 0.0_wp )  THEN
       binary_case = 1
       IF ( rh > 0.1_wp  .AND.  rh < 0.98 )  THEN
          IF ( binary_case == 1 )  THEN
             binary_hcl = - 5.0179_wp * rh**3 + 9.8816_wp * rh**2 - 10.789_wp * rh + 5.4737_wp
          ELSEIF ( binary_case == 2 )  THEN
             binary_hcl = - 4.6221_wp * rh + 4.2633_wp
          ENDIF
       ELSEIF ( rh >= 0.98_wp  .AND.  rh < 0.9999_wp )  THEN
          binary_hcl = 775.6111008626_wp * rh**3 - 2146.01320888771_wp * rh**2 +                   &
                       1969.01979670259_wp *  rh - 598.878230033926_wp
       ENDIF
    ENDIF

    IF ( nitric_acid > 0.0_wp )  THEN   ! HNO3
       IF ( full_complexity == 1  .OR.  rh <= 0.4_wp )  THEN
          hno3_hcl = 9.6256_wp * rh**4 - 26.507_wp * rh**3 + 27.622_wp * rh**2 - 12.958_wp * rh +  &
                     2.2193_wp
       ELSEIF ( full_complexity == 0  .AND.  rh > 0.4_wp )  THEN
          hno3_hcl = 1.3242_wp * rh**2 - 1.8827_wp * rh + 0.55706_wp
       ENDIF
    ENDIF

    IF ( sulphuric_acid > 0.0_wp )  THEN   ! H2SO4
       IF ( full_complexity == 1  .OR.  rh <= 0.4 )  THEN
          h2so4_hcl = 1.4406_wp * rh**3 - 2.7132_wp * rh**2 + 1.014_wp * rh + 0.25226_wp
       ELSEIF ( full_complexity == 0 .AND. rh > 0.4_wp ) THEN
          h2so4_hcl = 0.30993_wp * rh**2 - 0.99171_wp * rh + 0.66913_wp
       ENDIF
    ENDIF

    IF ( ammonium_sulphate > 0.0_wp )  THEN   ! NH42SO4
       nh42so4_hcl = 22.071_wp * rh**3 - 40.678_wp * rh**2 + 27.893_wp * rh - 9.4338_wp
    ENDIF

    IF ( ammonium_nitrate > 0.0_wp )  THEN   ! NH4NO3
       nh4no3_hcl = 19.935_wp * rh**3 - 42.335_wp * rh**2 + 31.275_wp * rh - 8.8675_wp
    ENDIF

    IF ( ammonium_chloride > 0.0_wp )  THEN   ! NH4Cl
       IF ( full_complexity == 1  .OR.  rh <= 0.4_wp )  THEN
          nh4cl_hcl = 2.8048_wp * rh**3 - 4.3182_wp * rh**2 + 3.1971_wp * rh - 1.6824_wp
       ELSEIF ( full_complexity == 0  .AND.  rh > 0.4_wp )  THEN
          nh4cl_hcl = 1.2304_wp * rh**2 - 0.18262_wp * rh - 1.0643_wp
       ENDIF
    ENDIF

    IF ( sodium_sulphate > 0.0_wp )  THEN   ! Na2SO4
       na2so4_hcl = 36.104_wp * rh**4 - 78.658_wp * rh**3 + 63.441_wp * rh**2 - 26.727_wp * rh +   &
                    5.7007_wp
    ENDIF

    IF ( sodium_nitrate > 0.0_wp )  THEN   ! NaNO3
       IF ( full_complexity == 1  .OR.  rh <= 0.4_wp )  THEN
          nano3_hcl = 54.471_wp * rh**5 - 159.42_wp * rh**4 + 180.25_wp * rh**3 - 98.176_wp * rh**2&
                      + 25.309_wp * rh - 2.4275_wp
       ELSEIF ( full_complexity == 0  .AND.  rh > 0.4_wp )  THEN
          nano3_hcl = 21.632_wp * rh**4 - 53.088_wp * rh**3 + 47.285_wp * rh**2 - 18.519_wp * rh   &
                      + 2.6846_wp
       ENDIF
    ENDIF

    IF ( sodium_chloride > 0.0_wp )  THEN   ! NaCl
       IF ( full_complexity == 1  .OR.  rh <= 0.4_wp )  THEN
          nacl_hcl = 5.4138_wp * rh**4 - 12.079_wp * rh**3 + 9.627_wp * rh**2 - 3.3164_wp * rh +   &
                     0.35224_wp
       ELSEIF ( full_complexity == 0  .AND.  rh > 0.4_wp )  THEN
          nacl_hcl = 2.432_wp * rh**3 - 4.3453_wp * rh**2 + 2.3834_wp * rh - 0.4762_wp
       ENDIF
    ENDIF

    ln_HCL_act = binary_hcl + nitric_acid_eq_frac * hno3_hcl + sulphuric_acid_eq_frac * h2so4_hcl +&
                 ammonium_sulphate_eq_frac * nh42so4_hcl + ammonium_nitrate_eq_frac * nh4no3_hcl + &
                 ammonium_chloride_eq_frac * nh4cl_hcl + sodium_sulphate_eq_frac * na2so4_hcl +    &
                 sodium_nitrate_eq_frac    * nano3_hcl + sodium_chloride_eq_frac   * nacl_hcl

     gamma_hcl    = EXP( ln_HCL_act )   ! Molal activity coefficient
     gamma_out(2) = gamma_hcl
!
!--  Equilibrium constant after Wagman et al. (1982) (and NIST database)
     k_hcl = 2E6_wp * EXP( 9000.0_wp * henrys_temp_dep )

     press_hcl = ( ions_mol(1) * ions_mol(7) * gamma_hcl**2 ) / k_hcl
!
!-- 5) Ion molility output
    mols_out = ions_mol

 END SUBROUTINE inorganic_pdfite

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Update the particle size distribution. Put particles into corrects bins.
!>
!> Moving-centre method assumed, i.e. particles are allowed to grow to their 
!> exact size as long as they are not crossing the fixed diameter bin limits. 
!> If the particles in a size bin cross the lower or upper diameter limit, they 
!> are all moved to the adjacent diameter bin and their volume is averaged with
!> the particles in the new bin, which then get a new diameter. 
!
!> Moving-centre method minimises numerical diffusion. 
!------------------------------------------------------------------------------!
 SUBROUTINE distr_update( paero )

    IMPLICIT NONE

    INTEGER(iwp) ::  ib      !< loop index
    INTEGER(iwp) ::  mm      !< loop index
    INTEGER(iwp) ::  counti  !< number of while loops

    LOGICAL  ::  within_bins !< logical (particle belongs to the bin?)

    REAL(wp) ::  znfrac  !< number fraction to be moved to the larger bin
    REAL(wp) ::  zvfrac  !< volume fraction to be moved to the larger bin
    REAL(wp) ::  zvexc   !< Volume in the grown bin which exceeds the bin upper limit
    REAL(wp) ::  zvihi   !< particle volume at the high end of the bin
    REAL(wp) ::  zvilo   !< particle volume at the low end of the bin
    REAL(wp) ::  zvpart  !< particle volume (m3)
    REAL(wp) ::  zvrat   !< volume ratio of a size bin

    real(wp), dimension(nbins_aerosol) ::  dummy

    TYPE(t_section), DIMENSION(nbins_aerosol), INTENT(inout) ::  paero !< aerosol properties

    zvpart      = 0.0_wp
    zvfrac      = 0.0_wp
    within_bins = .FALSE.

    dummy = paero(:)%numc
!
!-- Check if the volume of the bin is within bin limits after update
    counti = 0
    DO  WHILE ( .NOT. within_bins )
       within_bins = .TRUE.
!
!--    Loop from larger to smaller size bins
       DO  ib = end_subrange_2b-1, start_subrange_1a, -1
          mm = 0
          IF ( paero(ib)%numc > nclim )  THEN
             zvpart = 0.0_wp
             zvfrac = 0.0_wp

             IF ( ib == end_subrange_2a )  CYCLE
!
!--          Dry volume
             zvpart = SUM( paero(ib)%volc(1:7) ) / paero(ib)%numc
!
!--          Smallest bin cannot decrease
             IF ( paero(ib)%vlolim > zvpart  .AND.  ib == start_subrange_1a ) CYCLE
!
!--          Decreasing bins
             IF ( paero(ib)%vlolim > zvpart )  THEN
                mm = ib - 1
                IF ( ib == start_subrange_2b )  mm = end_subrange_1a    ! 2b goes to 1a

                paero(mm)%numc = paero(mm)%numc + paero(ib)%numc
                paero(ib)%numc = 0.0_wp
                paero(mm)%volc(:) = paero(mm)%volc(:) + paero(ib)%volc(:)
                paero(ib)%volc(:) = 0.0_wp
                CYCLE
             ENDIF
!
!--          If size bin has not grown, cycle.
!--          Changed by Mona: compare to the arithmetic mean volume, as done originally. Now 
!--          particle volume is derived from the geometric mean diameter, not arithmetic (see 
!--          SUBROUTINE set_sizebins).
             IF ( zvpart <= api6 * ( ( aero(ib)%vhilim + aero(ib)%vlolim ) / ( 2.0_wp * api6 ) ) ) &
             CYCLE
!
!--          Avoid precision problems
             IF ( ABS( zvpart - api6 * paero(ib)%dmid**3 ) < 1.0E-35_wp )  CYCLE
!
!--          Volume ratio of the size bin
             zvrat = paero(ib)%vhilim / paero(ib)%vlolim
!
!--          Particle volume at the low end of the bin
             zvilo = 2.0_wp * zvpart / ( 1.0_wp + zvrat )
!
!--          Particle volume at the high end of the bin
             zvihi = zvrat * zvilo
!
!--          Volume in the grown bin which exceeds the bin upper limit
             zvexc = 0.5_wp * ( zvihi + paero(ib)%vhilim )
!
!--          Number fraction to be moved to the larger bin
             znfrac = MIN( 1.0_wp, ( zvihi - paero(ib)%vhilim) / ( zvihi - zvilo ) )
!
!--          Volume fraction to be moved to the larger bin
             zvfrac = MIN( 0.99_wp, znfrac * zvexc / zvpart )
             IF ( zvfrac < 0.0_wp )  THEN
                message_string = 'Error: zvfrac < 0'
                CALL message( 'salsa_mod: distr_update', 'PA0624', 1, 2, 0, 6, 0 )
             ENDIF
!
!--          Update bin
             mm = ib + 1
!
!--          Volume (cm3/cm3)
             paero(mm)%volc(:) = paero(mm)%volc(:) + znfrac * paero(ib)%numc * zvexc *             &
                                 paero(ib)%volc(:) / SUM( paero(ib)%volc(1:7) )
             paero(ib)%volc(:) = paero(ib)%volc(:) - znfrac * paero(ib)%numc * zvexc *             &
                                 paero(ib)%volc(:) / SUM( paero(ib)%volc(1:7) )

!--          Number concentration (#/m3)
             paero(mm)%numc = paero(mm)%numc + znfrac * paero(ib)%numc
             paero(ib)%numc = paero(ib)%numc * ( 1.0_wp - znfrac )

          ENDIF     ! nclim

          IF ( paero(ib)%numc > nclim )   THEN
             zvpart = SUM( paero(ib)%volc(1:7) ) / paero(ib)%numc  ! Note: dry volume!
             within_bins = ( paero(ib)%vlolim < zvpart  .AND. zvpart < paero(ib)%vhilim )
          ENDIF

       ENDDO ! - ib

       counti = counti + 1
       IF ( counti > 100 )  THEN
          message_string = 'Error: Aerosol bin update not converged'
          CALL message( 'salsa_mod: distr_update', 'PA0625', 1, 2, 0, 6, 0 )
       ENDIF

    ENDDO ! - within bins

 END SUBROUTINE distr_update

!------------------------------------------------------------------------------!
! Description:
! ------------
!> salsa_diagnostics: Update properties for the current timestep:
!>
!> Juha Tonttila, FMI, 2014
!> Tomi Raatikainen, FMI, 2016
!------------------------------------------------------------------------------!
 SUBROUTINE salsa_diagnostics( i, j )

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point_s

    IMPLICIT NONE

    INTEGER(iwp) ::  ib   !<
    INTEGER(iwp) ::  ic   !<
    INTEGER(iwp) ::  icc  !<
    INTEGER(iwp) ::  ig   !<
    INTEGER(iwp) ::  k    !<

    INTEGER(iwp), INTENT(in) ::  i  !<
    INTEGER(iwp), INTENT(in) ::  j  !<

    REAL(wp), DIMENSION(nzb:nzt+1) ::  flag          !< flag to mask topography
    REAL(wp), DIMENSION(nzb:nzt+1) ::  flag_zddry    !< flag to mask zddry
    REAL(wp), DIMENSION(nzb:nzt+1) ::  in_adn        !< air density (kg/m3)
    REAL(wp), DIMENSION(nzb:nzt+1) ::  in_p          !< pressure
    REAL(wp), DIMENSION(nzb:nzt+1) ::  in_t          !< temperature (K)
    REAL(wp), DIMENSION(nzb:nzt+1) ::  mcsum         !< sum of mass concentration
    REAL(wp), DIMENSION(nzb:nzt+1) ::  ppm_to_nconc  !< Conversion factor: ppm to #/m3
    REAL(wp), DIMENSION(nzb:nzt+1) ::  zddry         !< particle dry diameter
    REAL(wp), DIMENSION(nzb:nzt+1) ::  zvol          !< particle volume

    flag_zddry   = 0.0_wp
    in_adn       = 0.0_wp
    in_p         = 0.0_wp
    in_t         = 0.0_wp
    ppm_to_nconc = 1.0_wp
    zddry        = 0.0_wp
    zvol         = 0.0_wp

    !$OMP MASTER
    CALL cpu_log( log_point_s(94), 'salsa diagnostics ', 'start' )
    !$OMP END MASTER

!
!-- Calculate thermodynamic quantities needed in SALSA
    CALL salsa_thrm_ij( i, j, p_ij=in_p, temp_ij=in_t, adn_ij=in_adn )
!
!-- Calculate conversion factors for gas concentrations
    ppm_to_nconc = for_ppm_to_nconc * in_p / in_t
!
!-- Predetermine flag to mask topography
    flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(:,j,i), 0 ) )

    DO  ib = 1, nbins_aerosol   ! aerosol size bins
!
!--    Remove negative values
       aerosol_number(ib)%conc(:,j,i) = MAX( nclim, aerosol_number(ib)%conc(:,j,i) ) * flag
!
!--    Calculate total mass concentration per bin
       mcsum = 0.0_wp
       DO  ic = 1, ncomponents_mass
          icc = ( ic - 1 ) * nbins_aerosol + ib
          mcsum = mcsum + aerosol_mass(icc)%conc(:,j,i) * flag
          aerosol_mass(icc)%conc(:,j,i) = MAX( mclim, aerosol_mass(icc)%conc(:,j,i) ) * flag
       ENDDO
!
!--    Check that number and mass concentration match qualitatively
       IF ( ANY( aerosol_number(ib)%conc(:,j,i) > nclim  .AND. mcsum <= 0.0_wp ) )  THEN
          DO  k = nzb+1, nzt
             IF ( aerosol_number(ib)%conc(k,j,i) >= nclim  .AND. mcsum(k) <= 0.0_wp )  THEN
                aerosol_number(ib)%conc(k,j,i) = nclim * flag(k)
                DO  ic = 1, ncomponents_mass
                   icc = ( ic - 1 ) * nbins_aerosol + ib
                   aerosol_mass(icc)%conc(k,j,i) = mclim * flag(k)
                ENDDO
             ENDIF
          ENDDO
       ENDIF
!
!--    Update aerosol particle radius
       CALL bin_mixrat( 'dry', ib, i, j, zvol )
       zvol = zvol / arhoh2so4    ! Why on sulphate?
!
!--    Particles smaller then 0.1 nm diameter are set to zero
       zddry = ( zvol / MAX( nclim, aerosol_number(ib)%conc(:,j,i) ) / api6 )**0.33333333_wp
       flag_zddry = MERGE( 1.0_wp, 0.0_wp, ( zddry < 1.0E-10_wp  .AND.                             &
                           aerosol_number(ib)%conc(:,j,i) > nclim ) )
!
!--    Volatile species to the gas phase
       IF ( index_so4 > 0 .AND. lscndgas )  THEN
          ic = ( index_so4 - 1 ) * nbins_aerosol + ib
          IF ( salsa_gases_from_chem )  THEN
             ig = gas_index_chem(1)
             chem_species(ig)%conc(:,j,i) = chem_species(ig)%conc(:,j,i) +                         &
                                            aerosol_mass(ic)%conc(:,j,i) * avo * flag_zddry /      &
                                            ( amh2so4 * ppm_to_nconc ) * flag
          ELSE
             salsa_gas(1)%conc(:,j,i) = salsa_gas(1)%conc(:,j,i) + aerosol_mass(ic)%conc(:,j,i) /  &
                                        amh2so4 * avo * flag_zddry * flag
          ENDIF
       ENDIF
       IF ( index_oc > 0  .AND.  lscndgas )  THEN
          ic = ( index_oc - 1 ) * nbins_aerosol + ib
          IF ( salsa_gases_from_chem )  THEN
             ig = gas_index_chem(5)
             chem_species(ig)%conc(:,j,i) = chem_species(ig)%conc(:,j,i) +                         &
                                            aerosol_mass(ic)%conc(:,j,i) * avo * flag_zddry /      &
                                            ( amoc * ppm_to_nconc ) * flag
          ELSE
             salsa_gas(5)%conc(:,j,i) = salsa_gas(5)%conc(:,j,i) + aerosol_mass(ic)%conc(:,j,i) /  &
                                        amoc * avo * flag_zddry * flag
          ENDIF
       ENDIF
       IF ( index_no > 0  .AND.  lscndgas )  THEN
          ic = ( index_no - 1 ) * nbins_aerosol + ib
          IF ( salsa_gases_from_chem )  THEN
             ig = gas_index_chem(2)
             chem_species(ig)%conc(:,j,i) = chem_species(ig)%conc(:,j,i) +                         &
                                            aerosol_mass(ic)%conc(:,j,i) * avo * flag_zddry /      &
                                            ( amhno3 * ppm_to_nconc ) *flag
          ELSE
             salsa_gas(2)%conc(:,j,i) = salsa_gas(2)%conc(:,j,i) + aerosol_mass(ic)%conc(:,j,i) /  &
                                        amhno3 * avo * flag_zddry * flag
          ENDIF
       ENDIF
       IF ( index_nh > 0  .AND.  lscndgas )  THEN
          ic = ( index_nh - 1 ) * nbins_aerosol + ib
          IF ( salsa_gases_from_chem )  THEN
             ig = gas_index_chem(3)
             chem_species(ig)%conc(:,j,i) = chem_species(ig)%conc(:,j,i) +                         &
                                            aerosol_mass(ic)%conc(:,j,i) * avo * flag_zddry /      &
                                            ( amnh3 * ppm_to_nconc ) *flag
          ELSE
             salsa_gas(3)%conc(:,j,i) = salsa_gas(3)%conc(:,j,i) + aerosol_mass(ic)%conc(:,j,i) /  &
                                        amnh3 * avo * flag_zddry *flag
          ENDIF
       ENDIF
!
!--    Mass and number to zero (insoluble species and water are lost)
       DO  ic = 1, ncomponents_mass
          icc = ( ic - 1 ) * nbins_aerosol + ib
          aerosol_mass(icc)%conc(:,j,i) = MERGE( mclim * flag, aerosol_mass(icc)%conc(:,j,i),      &
                                                 flag_zddry > 0.0_wp )
       ENDDO
       aerosol_number(ib)%conc(:,j,i) = MERGE( nclim * flag, aerosol_number(ib)%conc(:,j,i),       &
                                               flag_zddry > 0.0_wp )
       ra_dry(:,j,i,ib) = MAX( 1.0E-10_wp, 0.5_wp * zddry )

    ENDDO
    IF ( .NOT. salsa_gases_from_chem )  THEN
       DO  ig = 1, ngases_salsa
          salsa_gas(ig)%conc(:,j,i) = MAX( nclim, salsa_gas(ig)%conc(:,j,i) ) * flag
       ENDDO
    ENDIF

   !$OMP MASTER
    CALL cpu_log( log_point_s(94), 'salsa diagnostics ', 'stop' )
   !$OMP END MASTER

 END SUBROUTINE salsa_diagnostics


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!------------------------------------------------------------------------------!
 SUBROUTINE salsa_actions( location )


    CHARACTER (LEN=*), INTENT(IN) ::  location !< call location string

    SELECT CASE ( location )

       CASE ( 'before_timestep' )

          IF ( ws_scheme_sca )  sums_salsa_ws_l = 0.0_wp

       CASE DEFAULT
          CONTINUE

    END SELECT

 END SUBROUTINE salsa_actions


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid points i,j
!------------------------------------------------------------------------------!

 SUBROUTINE salsa_actions_ij( i, j, location )


    INTEGER(iwp),      INTENT(IN) ::  i         !< grid index in x-direction
    INTEGER(iwp),      INTENT(IN) ::  j         !< grid index in y-direction
    CHARACTER (LEN=*), INTENT(IN) ::  location  !< call location string
    INTEGER(iwp)  ::  dummy  !< call location string

    IF ( salsa    )   dummy = i + j

    SELECT CASE ( location )

       CASE ( 'before_timestep' )

          IF ( ws_scheme_sca )  sums_salsa_ws_l = 0.0_wp

       CASE DEFAULT
          CONTINUE

    END SELECT


 END SUBROUTINE salsa_actions_ij

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!------------------------------------------------------------------------------!
 SUBROUTINE salsa_non_advective_processes

    USE cpulog,                                                                                    &
        ONLY:  cpu_log, log_point_s

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<

    IF ( time_since_reference_point >= skip_time_do_salsa )  THEN
       IF ( ( time_since_reference_point - last_salsa_time ) >= dt_salsa )  THEN
!
!--       Calculate aerosol dynamic processes. salsa_driver can be run with a longer time step.
          CALL cpu_log( log_point_s(90), 'salsa processes ', 'start' )
          DO  i = nxl, nxr
             DO  j = nys, nyn
                CALL salsa_diagnostics( i, j )
                CALL salsa_driver( i, j, 3 )
                CALL salsa_diagnostics( i, j )
             ENDDO
          ENDDO
          CALL cpu_log( log_point_s(90), 'salsa processes ', 'stop' )
       ENDIF
    ENDIF

 END SUBROUTINE salsa_non_advective_processes


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid points i,j
!------------------------------------------------------------------------------!
 SUBROUTINE salsa_non_advective_processes_ij( i, j )

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point_s

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  i  !< grid index in x-direction
    INTEGER(iwp), INTENT(IN) ::  j  !< grid index in y-direction

    IF ( time_since_reference_point >= skip_time_do_salsa )  THEN
       IF ( ( time_since_reference_point - last_salsa_time ) >= dt_salsa )  THEN
!
!--       Calculate aerosol dynamic processes. salsa_driver can be run with a longer time step.
          CALL cpu_log( log_point_s(90), 'salsa processes ', 'start' )
          CALL salsa_diagnostics( i, j )
          CALL salsa_driver( i, j, 3 )
          CALL salsa_diagnostics( i, j )
          CALL cpu_log( log_point_s(90), 'salsa processes ', 'stop' )
       ENDIF
    ENDIF

 END SUBROUTINE salsa_non_advective_processes_ij

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Routine for exchange horiz of salsa variables.
!------------------------------------------------------------------------------!
 SUBROUTINE salsa_exchange_horiz_bounds

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point_s

    USE exchange_horiz_mod,                                                    &
        ONLY:  exchange_horiz

    IMPLICIT NONE

    INTEGER(iwp) ::  ib   !<
    INTEGER(iwp) ::  ic   !<
    INTEGER(iwp) ::  icc  !<
    INTEGER(iwp) ::  ig   !<

    IF ( time_since_reference_point >= skip_time_do_salsa )  THEN
       IF ( ( time_since_reference_point - last_salsa_time ) >= dt_salsa )  THEN

          CALL cpu_log( log_point_s(91), 'salsa exch-horiz ', 'start' )
!
!--       Exchange ghost points and decycle if needed.
          DO  ib = 1, nbins_aerosol
             CALL exchange_horiz( aerosol_number(ib)%conc, nbgp )
             CALL salsa_boundary_conds( aerosol_number(ib)%conc, aerosol_number(ib)%init )
             DO  ic = 1, ncomponents_mass
                icc = ( ic - 1 ) * nbins_aerosol + ib
                CALL exchange_horiz( aerosol_mass(icc)%conc, nbgp )
                CALL salsa_boundary_conds( aerosol_mass(icc)%conc, aerosol_mass(icc)%init )
             ENDDO
          ENDDO
          IF ( .NOT. salsa_gases_from_chem )  THEN
             DO  ig = 1, ngases_salsa
                CALL exchange_horiz( salsa_gas(ig)%conc, nbgp )
                CALL salsa_boundary_conds( salsa_gas(ig)%conc, salsa_gas(ig)%init )
             ENDDO
          ENDIF
          CALL cpu_log( log_point_s(91), 'salsa exch-horiz ', 'stop' )
!
!--       Update last_salsa_time
          last_salsa_time = time_since_reference_point
       ENDIF
    ENDIF

 END SUBROUTINE salsa_exchange_horiz_bounds

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate the prognostic equation for aerosol number and mass, and gas
!> concentrations. Cache-optimized.
!------------------------------------------------------------------------------!
 SUBROUTINE salsa_prognostic_equations_ij( i, j, i_omp_start, tn )

    IMPLICIT NONE

    INTEGER(iwp) ::  i            !<
    INTEGER(iwp) ::  i_omp_start  !<
    INTEGER(iwp) ::  ib           !< loop index for aerosol number bin OR gas index
    INTEGER(iwp) ::  ic           !< loop index for aerosol mass bin
    INTEGER(iwp) ::  icc          !< (c-1)*nbins_aerosol+b
    INTEGER(iwp) ::  ig           !< loop index for salsa gases
    INTEGER(iwp) ::  j            !<
    INTEGER(iwp) ::  tn           !<

    IF ( time_since_reference_point >= skip_time_do_salsa )  THEN
!
!--    Aerosol number
       DO  ib = 1, nbins_aerosol
!kk          sums_salsa_ws_l = aerosol_number(ib)%sums_ws_l
          CALL salsa_tendency( 'aerosol_number', aerosol_number(ib)%conc_p, aerosol_number(ib)%conc,&
                               aerosol_number(ib)%tconc_m, i, j, i_omp_start, tn, ib, ib,          &
                               aerosol_number(ib)%flux_s, aerosol_number(ib)%diss_s,               &
                               aerosol_number(ib)%flux_l, aerosol_number(ib)%diss_l,               &
                               aerosol_number(ib)%init, .TRUE. )
!kk          aerosol_number(ib)%sums_ws_l = sums_salsa_ws_l
!
!--       Aerosol mass
          DO  ic = 1, ncomponents_mass
             icc = ( ic - 1 ) * nbins_aerosol + ib
!kk             sums_salsa_ws_l = aerosol_mass(icc)%sums_ws_l
             CALL salsa_tendency( 'aerosol_mass', aerosol_mass(icc)%conc_p, aerosol_mass(icc)%conc,&
                                  aerosol_mass(icc)%tconc_m, i, j, i_omp_start, tn, ib, ic,        &
                                  aerosol_mass(icc)%flux_s, aerosol_mass(icc)%diss_s,              &
                                  aerosol_mass(icc)%flux_l, aerosol_mass(icc)%diss_l,              &
                                  aerosol_mass(icc)%init, .TRUE. )
!kk             aerosol_mass(icc)%sums_ws_l = sums_salsa_ws_l

          ENDDO  ! ic
       ENDDO  ! ib
!
!--    Gases
       IF ( .NOT. salsa_gases_from_chem )  THEN

          DO  ig = 1, ngases_salsa
!kk             sums_salsa_ws_l = salsa_gas(ig)%sums_ws_l
             CALL salsa_tendency( 'salsa_gas', salsa_gas(ig)%conc_p, salsa_gas(ig)%conc,           &
                                  salsa_gas(ig)%tconc_m, i, j, i_omp_start, tn, ig, ig,            &
                                  salsa_gas(ig)%flux_s, salsa_gas(ig)%diss_s, salsa_gas(ig)%flux_l,&
                                  salsa_gas(ig)%diss_l, salsa_gas(ig)%init, .FALSE. )
!kk             salsa_gas(ig)%sums_ws_l = sums_salsa_ws_l

          ENDDO  ! ig

       ENDIF

    ENDIF

 END SUBROUTINE salsa_prognostic_equations_ij
!
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate the prognostic equation for aerosol number and mass, and gas
!> concentrations. For vector machines.
!------------------------------------------------------------------------------!
 SUBROUTINE salsa_prognostic_equations()

    IMPLICIT NONE

    INTEGER(iwp) ::  ib           !< loop index for aerosol number bin OR gas index
    INTEGER(iwp) ::  ic           !< loop index for aerosol mass bin
    INTEGER(iwp) ::  icc          !< (c-1)*nbins_aerosol+b
    INTEGER(iwp) ::  ig           !< loop index for salsa gases

    IF ( time_since_reference_point >= skip_time_do_salsa )  THEN
!
!--    Aerosol number
       DO  ib = 1, nbins_aerosol
          sums_salsa_ws_l = aerosol_number(ib)%sums_ws_l
          CALL salsa_tendency( 'aerosol_number', aerosol_number(ib)%conc_p, aerosol_number(ib)%conc,&
                               aerosol_number(ib)%tconc_m, ib, ib, aerosol_number(ib)%init, .TRUE. )
          aerosol_number(ib)%sums_ws_l = sums_salsa_ws_l
!
!--       Aerosol mass
          DO  ic = 1, ncomponents_mass
             icc = ( ic - 1 ) * nbins_aerosol + ib
             sums_salsa_ws_l = aerosol_mass(icc)%sums_ws_l
             CALL salsa_tendency( 'aerosol_mass', aerosol_mass(icc)%conc_p, aerosol_mass(icc)%conc,&
                                  aerosol_mass(icc)%tconc_m, ib, ic, aerosol_mass(icc)%init, .TRUE. )
             aerosol_mass(icc)%sums_ws_l = sums_salsa_ws_l

          ENDDO  ! ic
       ENDDO  ! ib
!
!--    Gases
       IF ( .NOT. salsa_gases_from_chem )  THEN

          DO  ig = 1, ngases_salsa
             sums_salsa_ws_l = salsa_gas(ig)%sums_ws_l
             CALL salsa_tendency( 'salsa_gas', salsa_gas(ig)%conc_p, salsa_gas(ig)%conc,           &
                                  salsa_gas(ig)%tconc_m, ig, ig, salsa_gas(ig)%init, .FALSE. )
             salsa_gas(ig)%sums_ws_l = sums_salsa_ws_l

          ENDDO  ! ig

       ENDIF

    ENDIF

 END SUBROUTINE salsa_prognostic_equations
!
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Tendencies for aerosol number and mass and gas concentrations.
!> Cache-optimized.
!------------------------------------------------------------------------------!
 SUBROUTINE salsa_tendency_ij( id, rs_p, rs, trs_m, i, j, i_omp_start, tn, ib, ic, flux_s, diss_s, &
                               flux_l, diss_l, rs_init, do_sedimentation )

    USE advec_ws,                                                                                  &
        ONLY:  advec_s_ws

    USE advec_s_pw_mod,                                                                            &
        ONLY:  advec_s_pw

    USE advec_s_up_mod,                                                                            &
        ONLY:  advec_s_up

    USE arrays_3d,                                                                                 &
        ONLY:  ddzu, rdf_sc, tend

    USE diffusion_s_mod,                                                                           &
        ONLY:  diffusion_s

    USE indices,                                                                                   &
        ONLY:  wall_flags_total_0

    USE surface_mod,                                                                               &
        ONLY :  surf_def_h, surf_def_v, surf_lsm_h, surf_lsm_v, surf_usm_h, surf_usm_v

    IMPLICIT NONE

    CHARACTER(LEN = *) ::  id  !<

    INTEGER(iwp) ::  i            !<
    INTEGER(iwp) ::  i_omp_start  !<
    INTEGER(iwp) ::  ib           !< loop index for aerosol number bin OR gas index
    INTEGER(iwp) ::  ic           !< loop index for aerosol mass bin
    INTEGER(iwp) ::  icc          !< (c-1)*nbins_aerosol+b
    INTEGER(iwp) ::  j            !<
    INTEGER(iwp) ::  k            !<
    INTEGER(iwp) ::  tn           !<

    LOGICAL ::  do_sedimentation  !<

    REAL(wp), DIMENSION(nzb:nzt+1) ::  rs_init  !<

    REAL(wp), DIMENSION(nzb+1:nzt,0:threads_per_task-1) ::  diss_s  !<
    REAL(wp), DIMENSION(nzb+1:nzt,0:threads_per_task-1) ::  flux_s  !<

    REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn,0:threads_per_task-1) ::  diss_l  !<
    REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn,0:threads_per_task-1) ::  flux_l  !<

    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  rs_p    !<
    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  rs      !<
    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  trs_m   !<

    icc = ( ic - 1 ) * nbins_aerosol + ib
!
!-- Tendency-terms for reactive scalar
    tend(:,j,i) = 0.0_wp
!
!-- Advection terms
    IF ( timestep_scheme(1:5) == 'runge' )  THEN
       IF ( ws_scheme_sca )  THEN
          CALL advec_s_ws( salsa_advc_flags_s, i, j, rs, id, flux_s, diss_s, flux_l, diss_l,       &
                           i_omp_start, tn, bc_dirichlet_l  .OR.  bc_radiation_l,                  &
                           bc_dirichlet_n  .OR.  bc_radiation_n,                                   &
                           bc_dirichlet_r  .OR.  bc_radiation_r,                                   &
                           bc_dirichlet_s  .OR.  bc_radiation_s, monotonic_limiter_z )
       ELSE
          CALL advec_s_pw( i, j, rs )
       ENDIF
    ELSE
       CALL advec_s_up( i, j, rs )
    ENDIF
!
!-- Diffusion terms
    SELECT CASE ( id )
       CASE ( 'aerosol_number' )
          CALL diffusion_s( i, j, rs, surf_def_h(0)%answs(:,ib),                                   &
                                      surf_def_h(1)%answs(:,ib), surf_def_h(2)%answs(:,ib),        &
                                      surf_lsm_h%answs(:,ib),    surf_usm_h%answs(:,ib),           &
                                      surf_def_v(0)%answs(:,ib), surf_def_v(1)%answs(:,ib),        &
                                      surf_def_v(2)%answs(:,ib), surf_def_v(3)%answs(:,ib),        &
                                      surf_lsm_v(0)%answs(:,ib), surf_lsm_v(1)%answs(:,ib),        &
                                      surf_lsm_v(2)%answs(:,ib), surf_lsm_v(3)%answs(:,ib),        &
                                      surf_usm_v(0)%answs(:,ib), surf_usm_v(1)%answs(:,ib),        &
                                      surf_usm_v(2)%answs(:,ib), surf_usm_v(3)%answs(:,ib) )
       CASE ( 'aerosol_mass' )
          CALL diffusion_s( i, j, rs, surf_def_h(0)%amsws(:,icc),                                  &
                                      surf_def_h(1)%amsws(:,icc), surf_def_h(2)%amsws(:,icc),      &
                                      surf_lsm_h%amsws(:,icc),    surf_usm_h%amsws(:,icc),         &
                                      surf_def_v(0)%amsws(:,icc), surf_def_v(1)%amsws(:,icc),      &
                                      surf_def_v(2)%amsws(:,icc), surf_def_v(3)%amsws(:,icc),      &
                                      surf_lsm_v(0)%amsws(:,icc), surf_lsm_v(1)%amsws(:,icc),      &
                                      surf_lsm_v(2)%amsws(:,icc), surf_lsm_v(3)%amsws(:,icc),      &
                                      surf_usm_v(0)%amsws(:,icc), surf_usm_v(1)%amsws(:,icc),      &
                                      surf_usm_v(2)%amsws(:,icc), surf_usm_v(3)%amsws(:,icc) )
       CASE ( 'salsa_gas' )
          CALL diffusion_s( i, j, rs, surf_def_h(0)%gtsws(:,ib),                                   &
                                      surf_def_h(1)%gtsws(:,ib), surf_def_h(2)%gtsws(:,ib),        &
                                      surf_lsm_h%gtsws(:,ib), surf_usm_h%gtsws(:,ib),              &
                                      surf_def_v(0)%gtsws(:,ib), surf_def_v(1)%gtsws(:,ib),        &
                                      surf_def_v(2)%gtsws(:,ib), surf_def_v(3)%gtsws(:,ib),        &
                                      surf_lsm_v(0)%gtsws(:,ib), surf_lsm_v(1)%gtsws(:,ib),        &
                                      surf_lsm_v(2)%gtsws(:,ib), surf_lsm_v(3)%gtsws(:,ib),        &
                                      surf_usm_v(0)%gtsws(:,ib), surf_usm_v(1)%gtsws(:,ib),        &
                                      surf_usm_v(2)%gtsws(:,ib), surf_usm_v(3)%gtsws(:,ib) )
    END SELECT
!
!-- Sedimentation and prognostic equation for aerosol number and mass
    IF ( lsdepo  .AND.  do_sedimentation )  THEN
!DIR$ IVDEP
       DO  k = nzb+1, nzt
          tend(k,j,i) = tend(k,j,i) - MAX( 0.0_wp, ( rs(k+1,j,i) * sedim_vd(k+1,j,i,ib) -          &
                                                     rs(k,j,i) * sedim_vd(k,j,i,ib) ) * ddzu(k) )  &
                                    * MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k-1,j,i), 0 ) )
          rs_p(k,j,i) = rs(k,j,i) + ( dt_3d * ( tsc(2) * tend(k,j,i) + tsc(3) * trs_m(k,j,i) )     &
                                      - tsc(5) * rdf_sc(k) * ( rs(k,j,i) - rs_init(k) ) )          &
                                  * MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )
          IF ( rs_p(k,j,i) < 0.0_wp )  rs_p(k,j,i) = 0.1_wp * rs(k,j,i)
       ENDDO
    ELSE
!
!--    Prognostic equation
!DIR$ IVDEP
       DO  k = nzb+1, nzt
          rs_p(k,j,i) = rs(k,j,i) + ( dt_3d * ( tsc(2) * tend(k,j,i) + tsc(3) * trs_m(k,j,i) )     &
                                                - tsc(5) * rdf_sc(k) * ( rs(k,j,i) - rs_init(k) ) )&
                                  * MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )
          IF ( rs_p(k,j,i) < 0.0_wp )  rs_p(k,j,i) = 0.1_wp * rs(k,j,i)
       ENDDO
    ENDIF
!
!-- Calculate tendencies for the next Runge-Kutta step
    IF ( timestep_scheme(1:5) == 'runge' )  THEN
       IF ( intermediate_timestep_count == 1 )  THEN
          DO  k = nzb+1, nzt
             trs_m(k,j,i) = tend(k,j,i)
          ENDDO
       ELSEIF ( intermediate_timestep_count < intermediate_timestep_count_max )  THEN
          DO  k = nzb+1, nzt
             trs_m(k,j,i) = -9.5625_wp * tend(k,j,i) + 5.3125_wp * trs_m(k,j,i)
          ENDDO
       ENDIF
    ENDIF

 END SUBROUTINE salsa_tendency_ij
!
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate the tendencies for aerosol number and mass concentrations.
!> For vector machines.
!------------------------------------------------------------------------------!
 SUBROUTINE salsa_tendency( id, rs_p, rs, trs_m, ib, ic, rs_init, do_sedimentation )

    USE advec_ws,                                                                                  &
        ONLY:  advec_s_ws
    USE advec_s_pw_mod,                                                                            &
        ONLY:  advec_s_pw
    USE advec_s_up_mod,                                                                            &
        ONLY:  advec_s_up
    USE arrays_3d,                                                                                 &
        ONLY:  ddzu, rdf_sc, tend
    USE diffusion_s_mod,                                                                           &
        ONLY:  diffusion_s
    USE indices,                                                                                   &
        ONLY:  wall_flags_total_0
    USE surface_mod,                                                                               &
        ONLY :  surf_def_h, surf_def_v, surf_lsm_h, surf_lsm_v, surf_usm_h, surf_usm_v

    IMPLICIT NONE

    CHARACTER(LEN = *) ::  id

    INTEGER(iwp) ::  ib           !< loop index for aerosol number bin OR gas index
    INTEGER(iwp) ::  ic           !< loop index for aerosol mass bin
    INTEGER(iwp) ::  icc  !< (c-1)*nbins_aerosol+b
    INTEGER(iwp) ::  i    !<
    INTEGER(iwp) ::  j    !<
    INTEGER(iwp) ::  k    !<

    LOGICAL ::  do_sedimentation  !<

    REAL(wp), DIMENSION(nzb:nzt+1) ::  rs_init !<

    REAL(wp), DIMENSION(:,:,:), POINTER ::  rs_p    !<
    REAL(wp), DIMENSION(:,:,:), POINTER ::  rs      !<
    REAL(wp), DIMENSION(:,:,:), POINTER ::  trs_m   !<

    icc = ( ic - 1 ) * nbins_aerosol + ib
!
!-- Tendency-terms for reactive scalar
    tend = 0.0_wp
!
!-- Advection terms
    IF ( timestep_scheme(1:5) == 'runge' )  THEN
       IF ( ws_scheme_sca )  THEN
          CALL advec_s_ws( salsa_advc_flags_s, rs, id, bc_dirichlet_l  .OR.  bc_radiation_l,       &
                           bc_dirichlet_n  .OR.  bc_radiation_n,                                   &
                           bc_dirichlet_r  .OR.  bc_radiation_r,                                   &
                           bc_dirichlet_s  .OR.  bc_radiation_s )
       ELSE
          CALL advec_s_pw( rs )
       ENDIF
    ELSE
       CALL advec_s_up( rs )
    ENDIF
!
!-- Diffusion terms
    SELECT CASE ( id )
       CASE ( 'aerosol_number' )
          CALL diffusion_s( rs, surf_def_h(0)%answs(:,ib),                                         &
                                surf_def_h(1)%answs(:,ib), surf_def_h(2)%answs(:,ib),              &
                                surf_lsm_h%answs(:,ib),    surf_usm_h%answs(:,ib),                 &
                                surf_def_v(0)%answs(:,ib), surf_def_v(1)%answs(:,ib),              &
                                surf_def_v(2)%answs(:,ib), surf_def_v(3)%answs(:,ib),              &
                                surf_lsm_v(0)%answs(:,ib), surf_lsm_v(1)%answs(:,ib),              &
                                surf_lsm_v(2)%answs(:,ib), surf_lsm_v(3)%answs(:,ib),              &
                                surf_usm_v(0)%answs(:,ib), surf_usm_v(1)%answs(:,ib),              &
                                surf_usm_v(2)%answs(:,ib), surf_usm_v(3)%answs(:,ib) )
       CASE ( 'aerosol_mass' )
          CALL diffusion_s( rs, surf_def_h(0)%amsws(:,icc),                                        &
                                surf_def_h(1)%amsws(:,icc), surf_def_h(2)%amsws(:,icc),            &
                                surf_lsm_h%amsws(:,icc),    surf_usm_h%amsws(:,icc),               &
                                surf_def_v(0)%amsws(:,icc), surf_def_v(1)%amsws(:,icc),            &
                                surf_def_v(2)%amsws(:,icc), surf_def_v(3)%amsws(:,icc),            &
                                surf_lsm_v(0)%amsws(:,icc), surf_lsm_v(1)%amsws(:,icc),            &
                                surf_lsm_v(2)%amsws(:,icc), surf_lsm_v(3)%amsws(:,icc),            &
                                surf_usm_v(0)%amsws(:,icc), surf_usm_v(1)%amsws(:,icc),            &
                                surf_usm_v(2)%amsws(:,icc), surf_usm_v(3)%amsws(:,icc) )
       CASE ( 'salsa_gas' )
          CALL diffusion_s( rs, surf_def_h(0)%gtsws(:,ib),                                         &
                                surf_def_h(1)%gtsws(:,ib), surf_def_h(2)%gtsws(:,ib),              &
                                surf_lsm_h%gtsws(:,ib),    surf_usm_h%gtsws(:,ib),                 &
                                surf_def_v(0)%gtsws(:,ib), surf_def_v(1)%gtsws(:,ib),              &
                                surf_def_v(2)%gtsws(:,ib), surf_def_v(3)%gtsws(:,ib),              &
                                surf_lsm_v(0)%gtsws(:,ib), surf_lsm_v(1)%gtsws(:,ib),              &
                                surf_lsm_v(2)%gtsws(:,ib), surf_lsm_v(3)%gtsws(:,ib),              &
                                surf_usm_v(0)%gtsws(:,ib), surf_usm_v(1)%gtsws(:,ib),              &
                                surf_usm_v(2)%gtsws(:,ib), surf_usm_v(3)%gtsws(:,ib) )
    END SELECT
!
!-- Prognostic equation for a scalar
    DO  i = nxl, nxr
       DO  j = nys, nyn
!
!--       Sedimentation for aerosol number and mass
          IF ( lsdepo  .AND.  do_sedimentation )  THEN
             tend(nzb+1:nzt,j,i) = tend(nzb+1:nzt,j,i) - MAX( 0.0_wp, ( rs(nzb+2:nzt+1,j,i) *      &
                                   sedim_vd(nzb+2:nzt+1,j,i,ib) - rs(nzb+1:nzt,j,i) *              &
                                   sedim_vd(nzb+1:nzt,j,i,ib) ) * ddzu(nzb+1:nzt) ) *              &
                                   MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(nzb:nzt-1,j,i), 0 ) )
          ENDIF
          DO  k = nzb+1, nzt
             rs_p(k,j,i) = rs(k,j,i) +  ( dt_3d  * ( tsc(2) * tend(k,j,i) + tsc(3) * trs_m(k,j,i) )&
                                                  - tsc(5) * rdf_sc(k) * ( rs(k,j,i) - rs_init(k) )&
                                        ) * MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )
             IF ( rs_p(k,j,i) < 0.0_wp )  rs_p(k,j,i) = 0.1_wp * rs(k,j,i)
          ENDDO
       ENDDO
    ENDDO
!
!-- Calculate tendencies for the next Runge-Kutta step
    IF ( timestep_scheme(1:5) == 'runge' )  THEN
       IF ( intermediate_timestep_count == 1 )  THEN
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   trs_m(k,j,i) = tend(k,j,i)
                ENDDO
             ENDDO
          ENDDO
       ELSEIF ( intermediate_timestep_count < intermediate_timestep_count_max )  THEN
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   trs_m(k,j,i) =  -9.5625_wp * tend(k,j,i) + 5.3125_wp * trs_m(k,j,i)
                ENDDO
             ENDDO
          ENDDO
       ENDIF
    ENDIF

 END SUBROUTINE salsa_tendency


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Boundary conditions for prognostic variables in SALSA from module interface
!------------------------------------------------------------------------------!
 SUBROUTINE salsa_boundary_conditions

    IMPLICIT NONE

    INTEGER(iwp) ::  ib              !< index for aerosol size bins
    INTEGER(iwp) ::  ic              !< index for aerosol mass bins
    INTEGER(iwp) ::  icc             !< additional index for aerosol mass bins
    INTEGER(iwp) ::  ig              !< index for salsa gases


!
!-- moved from boundary_conds
    CALL salsa_boundary_conds
!
!-- Boundary conditions for prognostic quantitites of other modules:
!-- Here, only decycling is carried out
    IF ( time_since_reference_point >= skip_time_do_salsa )  THEN

       DO  ib = 1, nbins_aerosol
          CALL salsa_boundary_conds( aerosol_number(ib)%conc_p, aerosol_number(ib)%init )
          DO  ic = 1, ncomponents_mass
             icc = ( ic - 1 ) * nbins_aerosol + ib
             CALL salsa_boundary_conds( aerosol_mass(icc)%conc_p, aerosol_mass(icc)%init )
          ENDDO
       ENDDO
       IF ( .NOT. salsa_gases_from_chem )  THEN
          DO  ig = 1, ngases_salsa
             CALL salsa_boundary_conds( salsa_gas(ig)%conc_p, salsa_gas(ig)%init )
          ENDDO
       ENDIF

    ENDIF

 END SUBROUTINE salsa_boundary_conditions

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Boundary conditions for prognostic variables in SALSA
!------------------------------------------------------------------------------!
 SUBROUTINE salsa_boundary_conds

    USE arrays_3d,                                                                                 &
        ONLY:  dzu

    USE surface_mod,                                                                               &
        ONLY :  bc_h

    IMPLICIT NONE

    INTEGER(iwp) ::  i    !< grid index x direction
    INTEGER(iwp) ::  ib   !< index for aerosol size bins
    INTEGER(iwp) ::  ic   !< index for chemical compounds in aerosols
    INTEGER(iwp) ::  icc  !< additional index for chemical compounds in aerosols
    INTEGER(iwp) ::  ig   !< idex for gaseous compounds
    INTEGER(iwp) ::  j    !< grid index y direction
    INTEGER(iwp) ::  k    !< grid index y direction
    INTEGER(iwp) ::  l    !< running index boundary type, for up- and downward-facing walls
    INTEGER(iwp) ::  m    !< running index surface elements

    IF ( time_since_reference_point >= skip_time_do_salsa )  THEN
!
!--    Surface conditions:
       IF ( ibc_salsa_b == 0 )  THEN   ! Dirichlet
!
!--       Run loop over all non-natural and natural walls. Note, in wall-datatype the k coordinate 
!--       belongs to the atmospheric grid point, therefore, set s_p at k-1
          DO  l = 0, 1
             !$OMP PARALLEL PRIVATE( ib, ic, icc, ig, i, j, k )
             !$OMP DO
             DO  m = 1, bc_h(l)%ns

                i = bc_h(l)%i(m)
                j = bc_h(l)%j(m)
                k = bc_h(l)%k(m)

                DO  ib = 1, nbins_aerosol
                   aerosol_number(ib)%conc_p(k+bc_h(l)%koff,j,i) =             &
                                    aerosol_number(ib)%conc(k+bc_h(l)%koff,j,i)
                   DO  ic = 1, ncomponents_mass
                      icc = ( ic - 1 ) * nbins_aerosol + ib
                      aerosol_mass(icc)%conc_p(k+bc_h(l)%koff,j,i) =           &
                                    aerosol_mass(icc)%conc(k+bc_h(l)%koff,j,i)
                   ENDDO
                ENDDO
                IF ( .NOT. salsa_gases_from_chem )  THEN
                   DO  ig = 1, ngases_salsa
                      salsa_gas(ig)%conc_p(k+bc_h(l)%koff,j,i) =               &
                                    salsa_gas(ig)%conc(k+bc_h(l)%koff,j,i)
                   ENDDO
                ENDIF

             ENDDO
             !$OMP END PARALLEL

          ENDDO

       ELSE   ! Neumann

          DO l = 0, 1
             !$OMP PARALLEL PRIVATE( ib, ic, icc, ig, i, j, k )
             !$OMP DO
             DO  m = 1, bc_h(l)%ns

                i = bc_h(l)%i(m)
                j = bc_h(l)%j(m)
                k = bc_h(l)%k(m)

                DO  ib = 1, nbins_aerosol
                   aerosol_number(ib)%conc_p(k+bc_h(l)%koff,j,i) =             &
                                               aerosol_number(ib)%conc_p(k,j,i)
                   DO  ic = 1, ncomponents_mass
                      icc = ( ic - 1 ) * nbins_aerosol + ib
                      aerosol_mass(icc)%conc_p(k+bc_h(l)%koff,j,i) =           &
                                               aerosol_mass(icc)%conc_p(k,j,i)
                   ENDDO
                ENDDO
                IF ( .NOT. salsa_gases_from_chem ) THEN
                   DO  ig = 1, ngases_salsa
                      salsa_gas(ig)%conc_p(k+bc_h(l)%koff,j,i) =               &
                                               salsa_gas(ig)%conc_p(k,j,i)
                   ENDDO
                ENDIF

             ENDDO
             !$OMP END PARALLEL
          ENDDO

       ENDIF
!
!--   Top boundary conditions:
       IF ( ibc_salsa_t == 0 )  THEN   ! Dirichlet

          DO  ib = 1, nbins_aerosol
             aerosol_number(ib)%conc_p(nzt+1,:,:) = aerosol_number(ib)%conc(nzt+1,:,:)
             DO  ic = 1, ncomponents_mass
                icc = ( ic - 1 ) * nbins_aerosol + ib
                aerosol_mass(icc)%conc_p(nzt+1,:,:) = aerosol_mass(icc)%conc(nzt+1,:,:)
             ENDDO
          ENDDO
          IF ( .NOT. salsa_gases_from_chem )  THEN
             DO  ig = 1, ngases_salsa
                salsa_gas(ig)%conc_p(nzt+1,:,:) = salsa_gas(ig)%conc(nzt+1,:,:)
             ENDDO
          ENDIF

       ELSEIF ( ibc_salsa_t == 1 )  THEN   ! Neumann

          DO  ib = 1, nbins_aerosol
             aerosol_number(ib)%conc_p(nzt+1,:,:) = aerosol_number(ib)%conc_p(nzt,:,:)
             DO  ic = 1, ncomponents_mass
                icc = ( ic - 1 ) * nbins_aerosol + ib
                aerosol_mass(icc)%conc_p(nzt+1,:,:) = aerosol_mass(icc)%conc_p(nzt,:,:)
             ENDDO
          ENDDO
          IF ( .NOT. salsa_gases_from_chem )  THEN
             DO  ig = 1, ngases_salsa
                salsa_gas(ig)%conc_p(nzt+1,:,:) = salsa_gas(ig)%conc_p(nzt,:,:)
             ENDDO
          ENDIF

       ELSEIF ( ibc_salsa_t == 2 )  THEN   ! Initial gradient

          DO  ib = 1, nbins_aerosol
             aerosol_number(ib)%conc_p(nzt+1,:,:) = aerosol_number(ib)%conc_p(nzt,:,:) +           &
                                                    bc_an_t_val(ib) * dzu(nzt+1)
             DO  ic = 1, ncomponents_mass
                icc = ( ic - 1 ) * nbins_aerosol + ib
                aerosol_mass(icc)%conc_p(nzt+1,:,:) = aerosol_mass(icc)%conc_p(nzt,:,:) +          &
                                                      bc_am_t_val(icc) * dzu(nzt+1)
             ENDDO
          ENDDO
          IF ( .NOT. salsa_gases_from_chem )  THEN
             DO  ig = 1, ngases_salsa
                salsa_gas(ig)%conc_p(nzt+1,:,:) = salsa_gas(ig)%conc_p(nzt,:,:) +                  &
                                                  bc_gt_t_val(ig) * dzu(nzt+1)
             ENDDO
          ENDIF

       ENDIF
!
!--    Lateral boundary conditions at the outflow
       IF ( bc_radiation_s )  THEN
          DO  ib = 1, nbins_aerosol
             aerosol_number(ib)%conc_p(:,nys-1,:) = aerosol_number(ib)%conc_p(:,nys,:)
             DO  ic = 1, ncomponents_mass
                icc = ( ic - 1 ) * nbins_aerosol + ib
                aerosol_mass(icc)%conc_p(:,nys-1,:) = aerosol_mass(icc)%conc_p(:,nys,:)
             ENDDO
          ENDDO
          IF ( .NOT. salsa_gases_from_chem )  THEN
             DO  ig = 1, ngases_salsa
                salsa_gas(ig)%conc_p(:,nys-1,:) = salsa_gas(ig)%conc_p(:,nys,:)
             ENDDO
          ENDIF

       ELSEIF ( bc_radiation_n )  THEN
          DO  ib = 1, nbins_aerosol
             aerosol_number(ib)%conc_p(:,nyn+1,:) = aerosol_number(ib)%conc_p(:,nyn,:)
             DO  ic = 1, ncomponents_mass
                icc = ( ic - 1 ) * nbins_aerosol + ib
                aerosol_mass(icc)%conc_p(:,nyn+1,:) = aerosol_mass(icc)%conc_p(:,nyn,:)
             ENDDO
          ENDDO
          IF ( .NOT. salsa_gases_from_chem )  THEN
             DO  ig = 1, ngases_salsa
                salsa_gas(ig)%conc_p(:,nyn+1,:) = salsa_gas(ig)%conc_p(:,nyn,:)
             ENDDO
          ENDIF

       ELSEIF ( bc_radiation_l )  THEN
          DO  ib = 1, nbins_aerosol
             aerosol_number(ib)%conc_p(:,:,nxl-1) = aerosol_number(ib)%conc_p(:,:,nxl)
             DO  ic = 1, ncomponents_mass
                icc = ( ic - 1 ) * nbins_aerosol + ib
                aerosol_mass(icc)%conc_p(:,:,nxl-1) = aerosol_mass(icc)%conc_p(:,:,nxl)
             ENDDO
          ENDDO
          IF ( .NOT. salsa_gases_from_chem )  THEN
             DO  ig = 1, ngases_salsa
                salsa_gas(ig)%conc_p(:,:,nxl-1) = salsa_gas(ig)%conc_p(:,:,nxl)
             ENDDO
          ENDIF

       ELSEIF ( bc_radiation_r )  THEN
          DO  ib = 1, nbins_aerosol
             aerosol_number(ib)%conc_p(:,:,nxr+1) = aerosol_number(ib)%conc_p(:,:,nxr)
             DO  ic = 1, ncomponents_mass
                icc = ( ic - 1 ) * nbins_aerosol + ib
                aerosol_mass(icc)%conc_p(:,:,nxr+1) = aerosol_mass(icc)%conc_p(:,:,nxr)
             ENDDO
          ENDDO
          IF ( .NOT. salsa_gases_from_chem )  THEN
             DO  ig = 1, ngases_salsa
                salsa_gas(ig)%conc_p(:,:,nxr+1) = salsa_gas(ig)%conc_p(:,:,nxr)
             ENDDO
          ENDIF

       ENDIF

    ENDIF

 END SUBROUTINE salsa_boundary_conds

!------------------------------------------------------------------------------!
! Description:
! ------------
! Undoing of the previously done cyclic boundary conditions.
!------------------------------------------------------------------------------!
 SUBROUTINE salsa_boundary_conds_decycle ( sq, sq_init )

    USE control_parameters,                                                                        &
        ONLY:  nesting_offline

    IMPLICIT NONE

    INTEGER(iwp) ::  boundary  !<
    INTEGER(iwp) ::  ee        !<
    INTEGER(iwp) ::  copied    !<
    INTEGER(iwp) ::  i         !<
    INTEGER(iwp) ::  j         !<
    INTEGER(iwp) ::  k         !<
    INTEGER(iwp) ::  ss        !<

    REAL(wp) ::  flag  !< flag to mask topography grid points

    REAL(wp), DIMENSION(nzb:nzt+1) ::  sq_init  !< initial concentration profile

    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  sq  !< concentration array

    flag = 0.0_wp
!
!-- Skip input if forcing from a larger-scale models is applied.
    IF ( nesting_offline  .AND.  nesting_offline_salsa )  RETURN
!
!-- Left and right boundaries
    IF ( decycle_salsa_lr  .AND.  ( bc_lr_cyc  .OR. bc_lr == 'nested' ) )  THEN

       DO  boundary = 1, 2

          IF ( decycle_method_salsa(boundary) == 'dirichlet' )  THEN
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
                      sq(k,j,i) = sq_init(k) * flag
                   ENDDO
                ENDDO
             ENDDO

          ELSEIF ( decycle_method_salsa(boundary) == 'neumann' )  THEN
!
!--          The value at the boundary is copied to the ghost layers to simulate an outlet with
!--          zero gradient
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
                      sq(k,j,i) = sq(k,j,copied) * flag
                   ENDDO
                ENDDO
             ENDDO

          ELSE
             WRITE(message_string,*) 'unknown decycling method: decycle_method_salsa (', boundary, &
                                     ') ="' // TRIM( decycle_method_salsa(boundary) ) // '"'
             CALL message( 'salsa_boundary_conds_decycle', 'PA0626', 1, 2, 0, 6, 0 )
          ENDIF
       ENDDO
    ENDIF

!
!-- South and north boundaries
     IF ( decycle_salsa_ns  .AND.  ( bc_ns_cyc  .OR. bc_ns == 'nested' ) )  THEN

       DO  boundary = 3, 4

          IF ( decycle_method_salsa(boundary) == 'dirichlet' )  THEN
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
                      sq(k,j,i) = sq_init(k) * flag
                   ENDDO
                ENDDO
             ENDDO

          ELSEIF ( decycle_method_salsa(boundary) == 'neumann' )  THEN
!
!--          The value at the boundary is copied to the ghost layers to simulate an outlet with 
!--          zero gradient
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
                      sq(k,j,i) = sq(k,copied,i) * flag
                   ENDDO
                ENDDO
             ENDDO

          ELSE
             WRITE(message_string,*) 'unknown decycling method: decycle_method_salsa (', boundary, &
                                     ') ="' // TRIM( decycle_method_salsa(boundary) ) // '"'
             CALL message( 'salsa_boundary_conds_decycle', 'PA0627', 1, 2, 0, 6, 0 )
          ENDIF
       ENDDO
    ENDIF

 END SUBROUTINE salsa_boundary_conds_decycle

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculates the total dry or wet mass concentration for individual bins
!> Juha Tonttila (FMI) 2015
!> Tomi Raatikainen (FMI) 2016
!------------------------------------------------------------------------------!
 SUBROUTINE bin_mixrat( itype, ibin, i, j, mconc )

    IMPLICIT NONE

    CHARACTER(len=*), INTENT(in) ::  itype  !< 'dry' or 'wet'

    INTEGER(iwp) ::  ic                 !< loop index for mass bin number
    INTEGER(iwp) ::  iend               !< end index: include water or not

    INTEGER(iwp), INTENT(in) ::  ibin   !< index of the chemical component
    INTEGER(iwp), INTENT(in) ::  i      !< loop index for x-direction
    INTEGER(iwp), INTENT(in) ::  j      !< loop index for y-direction

    REAL(wp), DIMENSION(:), INTENT(out) ::  mconc  !< total dry or wet mass concentration

!-- Number of components 
    IF ( itype == 'dry' )  THEN
       iend = prtcl%ncomp - 1 
    ELSE IF ( itype == 'wet' )  THEN
       iend = prtcl%ncomp
    ELSE
       message_string = 'Error in itype!'
       CALL message( 'bin_mixrat', 'PA0628', 2, 2, 0, 6, 0 )
    ENDIF

    mconc = 0.0_wp

    DO  ic = ibin, iend*nbins_aerosol+ibin, nbins_aerosol !< every nbins'th element
       mconc = mconc + aerosol_mass(ic)%conc(:,j,i)
    ENDDO

 END SUBROUTINE bin_mixrat

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Sets surface fluxes
!------------------------------------------------------------------------------!
 SUBROUTINE salsa_emission_update

    IMPLICIT NONE

    IF ( include_emission )  THEN

       IF ( time_since_reference_point >= skip_time_do_salsa  )  THEN

          IF ( next_aero_emission_update <=                                                        &
               MAX( time_since_reference_point, 0.0_wp ) )  THEN
             CALL salsa_emission_setup( .FALSE. )
          ENDIF

          IF ( next_gas_emission_update <=                                                         &
               MAX( time_since_reference_point, 0.0_wp ) )  THEN
             IF ( salsa_emission_mode == 'read_from_file'  .AND.  .NOT. salsa_gases_from_chem )    &
             THEN
                CALL salsa_gas_emission_setup( .FALSE. )
             ENDIF
          ENDIF

       ENDIF
    ENDIF

 END SUBROUTINE salsa_emission_update

!------------------------------------------------------------------------------!
!> Description:
!> ------------
!> Define aerosol fluxes: constant or read from a from file
!> @todo - Emission stack height is not used yet. For default mode, emissions
!>         are assumed to occur on upward facing horizontal surfaces.
!------------------------------------------------------------------------------!
 SUBROUTINE salsa_emission_setup( init )

    USE control_parameters,                                                                        &
        ONLY:  end_time, spinup_time

    USE netcdf_data_input_mod,                                                                     &
        ONLY:  check_existence, close_input_file, get_attribute, get_variable,                     &
               inquire_num_variables, inquire_variable_names,                                      &
               get_dimension_length, open_read_file, street_type_f

    USE palm_date_time_mod,                                                                        &
        ONLY:  days_per_week, get_date_time, hours_per_day, months_per_year, seconds_per_hour

    USE surface_mod,                                                                               &
        ONLY:  surf_def_h, surf_lsm_h, surf_usm_h

    IMPLICIT NONE

    CHARACTER(LEN=80) ::  daytype = 'workday'  !< default day type
    CHARACTER(LEN=25) ::  in_name              !< name of a gas in the input file
    CHARACTER(LEN=25) ::  mod_name             !< name in the input file

    INTEGER(iwp) ::  day_of_month   !< day of the month
    INTEGER(iwp) ::  day_of_week    !< day of the week
    INTEGER(iwp) ::  day_of_year    !< day of the year
    INTEGER(iwp) ::  hour_of_day    !< hour of the day
    INTEGER(iwp) ::  i              !< loop index
    INTEGER(iwp) ::  ib             !< loop index: aerosol number bins
    INTEGER(iwp) ::  ic             !< loop index: aerosol chemical components
    INTEGER(iwp) ::  id_salsa       !< NetCDF id of aerosol emission input file
    INTEGER(iwp) ::  in             !< loop index: emission category
    INTEGER(iwp) ::  index_dd       !< index day
    INTEGER(iwp) ::  index_hh       !< index hour
    INTEGER(iwp) ::  index_mm       !< index month
    INTEGER(iwp) ::  inn            !< loop index
    INTEGER(iwp) ::  j              !< loop index
    INTEGER(iwp) ::  month_of_year  !< month of the year
    INTEGER(iwp) ::  ss             !< loop index

    INTEGER(iwp), DIMENSION(maxspec) ::  cc_i2m   !<

    LOGICAL  ::  netcdf_extend = .FALSE.  !< NetCDF input file exists

    LOGICAL, INTENT(in) ::  init  !< if .TRUE. --> initialisation call

    REAL(wp) ::  second_of_day  !< second of the day

    REAL(wp), DIMENSION(24) ::  par_emis_time_factor =  & !< time factors for the parameterized mode
                                                      (/ 0.009, 0.004, 0.004, 0.009, 0.029, 0.039, &
                                                         0.056, 0.053, 0.051, 0.051, 0.052, 0.055, &
                                                         0.059, 0.061, 0.064, 0.067, 0.069, 0.069, &
                                                         0.049, 0.039, 0.039, 0.029, 0.024, 0.019 /)

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  nsect_emission  !< sectional number emission

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  source_array  !< temporary source array

!
!-- Define emissions:
    SELECT CASE ( salsa_emission_mode )

       CASE ( 'uniform', 'parameterized' )

          IF ( init )  THEN  ! Do only once
!
!-           Form a sectional size distribution for the emissions
             ALLOCATE( nsect_emission(1:nbins_aerosol),                                            &
                       source_array(nys:nyn,nxl:nxr,1:nbins_aerosol) )
!
!--          Precalculate a size distribution for the emission based on the mean diameter, standard
!--          deviation and number concentration per each log-normal mode
             CALL size_distribution( surface_aerosol_flux, aerosol_flux_dpg, aerosol_flux_sigmag,  &
                                     nsect_emission )
             IF ( salsa_emission_mode == 'uniform' )  THEN
                DO  ib = 1, nbins_aerosol
                   source_array(:,:,ib) = nsect_emission(ib)
                ENDDO
             ELSE
!
!--             Get a time factor for the specific hour
                IF ( .NOT.  ALLOCATED( aero_emission_att%time_factor ) )                           &
                   ALLOCATE( aero_emission_att%time_factor(1) )
                CALL get_date_time( MAX( time_since_reference_point, 0.0_wp ), hour=hour_of_day )
                index_hh = hour_of_day
                aero_emission_att%time_factor(1) = par_emis_time_factor(index_hh+1)

                IF ( street_type_f%from_file )  THEN
                   DO  i = nxl, nxr
                      DO  j = nys, nyn
                         IF ( street_type_f%var(j,i) >= main_street_id  .AND.                      &
                              street_type_f%var(j,i) < max_street_id )  THEN
                            source_array(j,i,:) = nsect_emission(:) * emiss_factor_main *          &
                                                  aero_emission_att%time_factor(1)
                         ELSEIF ( street_type_f%var(j,i) >= side_street_id  .AND.                  &
                                  street_type_f%var(j,i) < main_street_id )  THEN
                            source_array(j,i,:) = nsect_emission(:) * emiss_factor_side *          &
                                                  aero_emission_att%time_factor(1)
                         ENDIF
                      ENDDO
                   ENDDO
                ELSE
                   WRITE( message_string, * ) 'salsa_emission_mode = "parameterized" but the '//  &
                                              'street_type data is missing.'
                   CALL message( 'salsa_emission_setup', 'PA0695', 1, 2, 0, 6, 0 )
                ENDIF
             ENDIF
!
!--          Check which chemical components are used
             cc_i2m = 0
             IF ( index_so4 > 0 ) cc_i2m(1) = index_so4
             IF ( index_oc > 0 )  cc_i2m(2) = index_oc
             IF ( index_bc > 0 )  cc_i2m(3) = index_bc
             IF ( index_du > 0 )  cc_i2m(4) = index_du
             IF ( index_ss > 0 )  cc_i2m(5) = index_ss
             IF ( index_no > 0 )  cc_i2m(6) = index_no
             IF ( index_nh > 0 )  cc_i2m(7) = index_nh
!
!--          Normalise mass fractions so that their sum is 1
             aerosol_flux_mass_fracs_a = aerosol_flux_mass_fracs_a /                               &
                                         SUM( aerosol_flux_mass_fracs_a(1:ncc ) )
             IF ( salsa_emission_mode ==  'uniform' )  THEN
!
!--             Set uniform fluxes of default horizontal surfaces
                CALL set_flux( surf_def_h(0), cc_i2m, aerosol_flux_mass_fracs_a, source_array )
             ELSE
!
!--             Set fluxes normalised based on the street type on land surfaces
                CALL set_flux( surf_lsm_h, cc_i2m, aerosol_flux_mass_fracs_a, source_array )
             ENDIF

             DEALLOCATE( nsect_emission, source_array )
          ENDIF

       CASE ( 'read_from_file' )
!
!--       Reset surface fluxes
          surf_def_h(0)%answs = 0.0_wp
          surf_def_h(0)%amsws = 0.0_wp
          surf_lsm_h%answs = 0.0_wp
          surf_lsm_h%amsws = 0.0_wp
          surf_usm_h%answs = 0.0_wp
          surf_usm_h%amsws = 0.0_wp

!
!--       Reset source arrays:
          DO  ib = 1, nbins_aerosol
             aerosol_number(ib)%source = 0.0_wp
          ENDDO

          DO  ic = 1, ncomponents_mass * nbins_aerosol
             aerosol_mass(ic)%source = 0.0_wp
          ENDDO

#if defined( __netcdf )
!
!--       Check existence of PIDS_SALSA file
          INQUIRE( FILE = TRIM( input_file_salsa ) // TRIM( coupling_char ), EXIST = netcdf_extend )
          IF ( .NOT. netcdf_extend )  THEN
             message_string = 'Input file '// TRIM( input_file_salsa ) //  TRIM( coupling_char )&
                              // ' missing!'
             CALL message( 'salsa_emission_setup', 'PA0629', 1, 2, 0, 6, 0 )
          ENDIF
!
!--       Open file in read-only mode
          CALL open_read_file( TRIM( input_file_salsa ) // TRIM( coupling_char ), id_salsa )

          IF ( init )  THEN
!
!--          Variable names
             CALL inquire_num_variables( id_salsa, aero_emission_att%num_vars )
             ALLOCATE( aero_emission_att%var_names(1:aero_emission_att%num_vars) )
             CALL inquire_variable_names( id_salsa, aero_emission_att%var_names )
!
!--          Read the index and name of chemical components
             CALL get_dimension_length( id_salsa, aero_emission_att%ncc, 'composition_index' )
             ALLOCATE( aero_emission_att%cc_index(1:aero_emission_att%ncc) )
             CALL get_variable( id_salsa, 'composition_index', aero_emission_att%cc_index )

             IF ( check_existence( aero_emission_att%var_names, 'composition_name' ) )  THEN
                CALL get_variable( id_salsa, 'composition_name', aero_emission_att%cc_name,        &
                                   aero_emission_att%ncc )
             ELSE
                message_string = 'Missing composition_name in ' // TRIM( input_file_salsa )
                CALL message( 'salsa_emission_setup', 'PA0657', 1, 2, 0, 6, 0 )
             ENDIF
!
!--          Find the corresponding chemical components in the model
             aero_emission_att%cc_in2mod = 0
             DO  ic = 1, aero_emission_att%ncc
                in_name = aero_emission_att%cc_name(ic)
                SELECT CASE ( TRIM( in_name ) )
                   CASE ( 'H2SO4', 'h2so4', 'SO4', 'so4' )
                      aero_emission_att%cc_in2mod(1) = ic
                   CASE ( 'OC', 'oc', 'organics' )
                      aero_emission_att%cc_in2mod(2) = ic
                   CASE ( 'BC', 'bc' )
                      aero_emission_att%cc_in2mod(3) = ic
                   CASE ( 'DU', 'du' )
                      aero_emission_att%cc_in2mod(4) = ic
                   CASE ( 'SS', 'ss' )
                      aero_emission_att%cc_in2mod(5) = ic
                   CASE ( 'HNO3', 'hno3', 'NO', 'no', 'NO3', 'no3' )
                      aero_emission_att%cc_in2mod(6) = ic
                   CASE ( 'NH3', 'nh3', 'NH', 'nh', 'NH4', 'nh4' )
                      aero_emission_att%cc_in2mod(7) = ic
                END SELECT

             ENDDO

             IF ( SUM( aero_emission_att%cc_in2mod ) == 0 )  THEN
                message_string = 'None of the aerosol chemical components in ' // TRIM(            &
                                 input_file_salsa ) // ' correspond to the ones applied in SALSA.'
                CALL message( 'salsa_emission_setup', 'PA0630', 1, 2, 0, 6, 0 )
             ENDIF
!
!--          Get number of emission categories
             CALL get_dimension_length( id_salsa, aero_emission_att%ncat, 'ncat' )
!
!--          Get the chemical composition (i.e. mass fraction of different species) in aerosols
             IF ( check_existence( aero_emission_att%var_names, 'emission_mass_fracs' ) )  THEN
                ALLOCATE( aero_emission%mass_fracs(1:aero_emission_att%ncat,                       &
                                                   1:aero_emission_att%ncc) )
                CALL get_variable( id_salsa, 'emission_mass_fracs', aero_emission%mass_fracs,      &
                                   0, aero_emission_att%ncc-1, 0, aero_emission_att%ncat-1 )
             ELSE
                message_string = 'Missing emission_mass_fracs in ' //  TRIM( input_file_salsa )
                CALL message( 'salsa_emission_setup', 'PA0694', 1, 2, 0, 6, 0 )
             ENDIF
!
!--          If the chemical component is not activated, set its mass fraction to 0 to avoid
!--          inbalance between number and mass flux
             cc_i2m = aero_emission_att%cc_in2mod
             IF ( index_so4 < 0  .AND.  cc_i2m(1) > 0 )                                            &
                aero_emission%mass_fracs(:,cc_i2m(1)) = 0.0_wp
             IF ( index_oc  < 0  .AND.  cc_i2m(2) > 0 )                                            &
                aero_emission%mass_fracs(:,cc_i2m(2)) = 0.0_wp
             IF ( index_bc  < 0  .AND.  cc_i2m(3) > 0 )                                            &
                aero_emission%mass_fracs(:,cc_i2m(3)) = 0.0_wp
             IF ( index_du  < 0  .AND.  cc_i2m(4) > 0 )                                            &
                aero_emission%mass_fracs(:,cc_i2m(4)) = 0.0_wp
             IF ( index_ss  < 0  .AND.  cc_i2m(5) > 0 )                                            &
                aero_emission%mass_fracs(:,cc_i2m(5)) = 0.0_wp
             IF ( index_no  < 0  .AND.  cc_i2m(6) > 0 )                                            &
                aero_emission%mass_fracs(:,cc_i2m(6)) = 0.0_wp
             IF ( index_nh  < 0  .AND.  cc_i2m(7) > 0 )                                            &
                aero_emission%mass_fracs(:,cc_i2m(7)) = 0.0_wp
!
!--          Then normalise the mass fraction so that SUM = 1
             DO  in = 1, aero_emission_att%ncat
                aero_emission%mass_fracs(in,:) = aero_emission%mass_fracs(in,:) /                  &
                                                 SUM( aero_emission%mass_fracs(in,:) )
             ENDDO
!
!--          Inquire the fill value
             CALL get_attribute( id_salsa, '_FillValue', aero_emission%fill, .FALSE.,              &
                                 'aerosol_emission_values' )
!
!--          Inquire units of emissions
             CALL get_attribute( id_salsa, 'units', aero_emission_att%units, .FALSE.,              &
                                 'aerosol_emission_values' )
!
!--          Inquire the level of detail (lod)
             CALL get_attribute( id_salsa, 'lod', aero_emission_att%lod, .FALSE.,                  &
                                 'aerosol_emission_values' )

!
!--          Read different emission information depending on the level of detail of emissions:

!
!--          Default mode:
             IF ( aero_emission_att%lod == 1 )  THEN
!
!--             Unit conversion factor: convert to SI units (kg/m2/s)
                IF ( aero_emission_att%units == 'kg/m2/yr' )  THEN
                   aero_emission_att%conversion_factor = 1.0_wp / 3600.0_wp
                ELSEIF ( aero_emission_att%units == 'g/m2/yr' )  THEN
                   aero_emission_att%conversion_factor = 0.001_wp / 3600.0_wp
                ELSE
                   message_string = 'unknown unit for aerosol emissions: ' //                      &
                                    TRIM( aero_emission_att%units ) // ' (lod1)'
                   CALL message( 'salsa_emission_setup','PA0631', 1, 2, 0, 6, 0 )
                ENDIF
!
!--             Allocate emission arrays
                ALLOCATE( aero_emission_att%cat_index(1:aero_emission_att%ncat),                   &
                          aero_emission_att%rho(1:aero_emission_att%ncat),                         &
                          aero_emission_att%time_factor(1:aero_emission_att%ncat) )
!
!--             Get emission category names and indices
                IF ( check_existence( aero_emission_att%var_names, 'emission_category_name' ) )  THEN
                   CALL get_variable( id_salsa, 'emission_category_name',                          &
                                      aero_emission_att%cat_name,  aero_emission_att%ncat )
                ELSE
                   message_string = 'Missing emission_category_name in ' // TRIM( input_file_salsa )
                   CALL message( 'salsa_emission_setup', 'PA0658', 1, 2, 0, 6, 0 )
                ENDIF
                CALL get_variable( id_salsa, 'emission_category_index', aero_emission_att%cat_index )
!
!--             Find corresponding emission categories
                DO  in = 1, aero_emission_att%ncat
                   in_name = aero_emission_att%cat_name(in)
                   DO  ss = 1, def_modes%ndc
                      mod_name = def_modes%cat_name_table(ss)
                      IF ( TRIM( in_name(1:4) ) == TRIM( mod_name(1:4 ) ) )  THEN
                         def_modes%cat_input_to_model(ss) = in
                      ENDIF
                   ENDDO
                ENDDO

                IF ( SUM( def_modes%cat_input_to_model ) == 0 )  THEN
                   message_string = 'None of the emission categories in ' //  TRIM(                &
                                    input_file_salsa ) // ' match with the ones in the model.'
                   CALL message( 'salsa_emission_setup', 'PA0632', 1, 2, 0, 6, 0 )
                ENDIF
!
!--             Emission time factors: Find check whether emission time factors are given for each
!--             hour of year OR based on month, day and hour
!
!--             For each hour of year:
                IF ( check_existence( aero_emission_att%var_names, 'nhoursyear' ) )  THEN
                   CALL get_dimension_length( id_salsa, aero_emission_att%nhoursyear, 'nhoursyear' )
                   ALLOCATE( aero_emission_att%etf(1:aero_emission_att%ncat,                       &
                                                   1:aero_emission_att%nhoursyear) )
                   CALL get_variable( id_salsa, 'emission_time_factors', aero_emission_att%etf,    &
                                    0, aero_emission_att%nhoursyear-1, 0, aero_emission_att%ncat-1 )
!
!--             Based on the month, day and hour:
                ELSEIF ( check_existence( aero_emission_att%var_names, 'nmonthdayhour' ) )  THEN
                   CALL get_dimension_length( id_salsa, aero_emission_att%nmonthdayhour,           &
                                              'nmonthdayhour' )
                   ALLOCATE( aero_emission_att%etf(1:aero_emission_att%ncat,                       &
                                                   1:aero_emission_att%nmonthdayhour) )
                   CALL get_variable( id_salsa, 'emission_time_factors', aero_emission_att%etf,    &
                                 0, aero_emission_att%nmonthdayhour-1, 0, aero_emission_att%ncat-1 )
                ELSE
                   message_string = 'emission_time_factors should be given for each nhoursyear ' //&
                                    'OR nmonthdayhour'
                   CALL message( 'salsa_emission_setup','PA0633', 1, 2, 0, 6, 0 )
                ENDIF
!
!--             Next emission update
                CALL get_date_time( time_since_reference_point, second_of_day=second_of_day )
                next_aero_emission_update = MOD( second_of_day, seconds_per_hour ) !- seconds_per_hour
!
!--             Calculate average mass density (kg/m3)
                aero_emission_att%rho = 0.0_wp

                IF ( cc_i2m(1) /= 0 )  aero_emission_att%rho = aero_emission_att%rho +  arhoh2so4 *&
                                                               aero_emission%mass_fracs(:,cc_i2m(1))
                IF ( cc_i2m(2) /= 0 )  aero_emission_att%rho = aero_emission_att%rho + arhooc *    &
                                                               aero_emission%mass_fracs(:,cc_i2m(2))
                IF ( cc_i2m(3) /= 0 )  aero_emission_att%rho = aero_emission_att%rho + arhobc *    &
                                                               aero_emission%mass_fracs(:,cc_i2m(3))
                IF ( cc_i2m(4) /= 0 )  aero_emission_att%rho = aero_emission_att%rho + arhodu *    &
                                                               aero_emission%mass_fracs(:,cc_i2m(4))
                IF ( cc_i2m(5) /= 0 )  aero_emission_att%rho = aero_emission_att%rho + arhoss *    &
                                                               aero_emission%mass_fracs(:,cc_i2m(5))
                IF ( cc_i2m(6) /= 0 )  aero_emission_att%rho = aero_emission_att%rho + arhohno3 *  &
                                                               aero_emission%mass_fracs(:,cc_i2m(6))
                IF ( cc_i2m(7) /= 0 )  aero_emission_att%rho = aero_emission_att%rho + arhonh3 *   &
                                                               aero_emission%mass_fracs(:,cc_i2m(7))
!
!--             Allocate and read surface emission data (in total PM, get_variable_3d_real)
                ALLOCATE( aero_emission%def_data(nys:nyn,nxl:nxr,1:aero_emission_att%ncat) )
                CALL get_variable( id_salsa, 'aerosol_emission_values', aero_emission%def_data,    &
                                   0, aero_emission_att%ncat-1, nxl, nxr, nys, nyn )

!
!--          Pre-processed mode
             ELSEIF ( aero_emission_att%lod == 2 )  THEN
!
!--             Unit conversion factor: convert to SI units (#/m2/s)
                IF ( aero_emission_att%units == '#/m2/s' )  THEN
                   aero_emission_att%conversion_factor = 1.0_wp
                ELSE
                   message_string = 'unknown unit for aerosol emissions: ' //                      &
                                    TRIM( aero_emission_att%units )
                   CALL message( 'salsa_emission_setup','PA0634', 1, 2, 0, 6, 0 )
                ENDIF
!
!--             Number of aerosol size bins in the emission data
                CALL get_dimension_length( id_salsa, aero_emission_att%nbins, 'Dmid' )
                IF ( aero_emission_att%nbins /= nbins_aerosol )  THEN
                   message_string = 'The number of size bins in aerosol input data does not ' //   &
                                    'correspond to the model set-up'
                   CALL message( 'salsa_emission_setup','PA0635', 1, 2, 0, 6, 0 )
                ENDIF
!
!--             Number of time steps in the emission data
                CALL get_dimension_length( id_salsa, aero_emission_att%nt, 'time')
!
!--             Allocate bin diameters, time and mass fraction array
                ALLOCATE( aero_emission_att%dmid(1:nbins_aerosol),                                 &
                          aero_emission_att%time(0:aero_emission_att%nt-1),                        &
                          aero_emission%num_fracs(1:aero_emission_att%ncat,1:nbins_aerosol) )
!
!--             Read mean diameters
                CALL get_variable( id_salsa, 'Dmid', aero_emission_att%dmid )
!
!--             Check whether the sectional representation of the aerosol size distribution conform
!--             to the one applied in the model
                IF ( ANY( ABS( ( aero(1:nbins_aerosol)%dmid - aero_emission_att%dmid ) /           &
                               aero(1:nbins_aerosol)%dmid ) > 0.1_wp )  )  THEN
                   message_string = 'Mean diameters of size bins in ' // TRIM( input_file_salsa )  &
                                    // ' do not match with the ones in the model.'
                   CALL message( 'salsa_emission_setup','PA0636', 1, 2, 0, 6, 0 )
                ENDIF
!
!--             Read time stamps:
                IF ( check_existence( aero_emission_att%var_names, 'time' ) )  THEN
                   CALL get_variable( id_salsa, 'time', aero_emission_att%time )
                ELSE
                   message_string = 'Missing time in ' //  TRIM( input_file_salsa )
                   CALL message( 'salsa_emission_setup', 'PA0660', 1, 2, 0, 6, 0 )
                ENDIF
!
!--             Check if the provided data covers the entire simulation. The spinup time is added
!--             to the end_time, this must be considered here.
                IF ( end_time - spinup_time > aero_emission_att%time(aero_emission_att%nt-1) )  THEN
                   message_string = 'end_time of the simulation exceeds the time dimension in ' // &
                                    'the salsa input file.'
                   CALL message( 'salsa_emission_setup', 'PA0692', 1, 2, 0, 6, 0 ) 
                ENDIF
!
!--             Read emission number fractions per category
                IF ( check_existence( aero_emission_att%var_names, 'emission_number_fracs' ) )  THEN
                   CALL get_variable( id_salsa, 'emission_number_fracs', aero_emission%num_fracs,  &
                                      0, nbins_aerosol-1, 0, aero_emission_att%ncat-1 )
                ELSE
                   message_string = 'Missing emission_number_fracs in ' //  TRIM( input_file_salsa )
                   CALL message( 'salsa_emission_setup', 'PA0694', 1, 2, 0, 6, 0 )
                ENDIF

             ELSE
                message_string = 'Unknown lod for aerosol_emission_values.'
                CALL message( 'salsa_emission_setup','PA0637', 1, 2, 0, 6, 0 )

             ENDIF  ! lod

          ENDIF  ! init
!
!--       Define and set current emission values:
!
!--       Default type emissions (aerosol emission given as total mass emission per year):
          IF ( aero_emission_att%lod == 1 )  THEN
!
!--          Emission time factors for each emission category at current time step
             IF ( aero_emission_att%nhoursyear > aero_emission_att%nmonthdayhour )  THEN
!
!--             Get the index of the current hour
                CALL get_date_time( MAX( 0.0_wp, time_since_reference_point ),                     &
                                    day_of_year=day_of_year, hour=hour_of_day )
                index_hh = ( day_of_year - 1_iwp ) * hours_per_day + hour_of_day
                aero_emission_att%time_factor = aero_emission_att%etf(:,index_hh+1)

             ELSEIF ( aero_emission_att%nhoursyear < aero_emission_att%nmonthdayhour )  THEN
!
!--             Get the index of current hour (index_hh) (TODO: Now "workday" is always assumed.
!--             Needs to be calculated.)
                CALL get_date_time( MAX( 0.0_wp, time_since_reference_point ), month=month_of_year,&
                                    day=day_of_month, hour=hour_of_day, day_of_week=day_of_week )
                index_mm = month_of_year
                index_dd = months_per_year + day_of_week
                SELECT CASE(TRIM(daytype))

                   CASE ("workday")
                      index_hh = months_per_year + days_per_week + hour_of_day

                   CASE ("weekend")
                      index_hh = months_per_year + days_per_week + hours_per_day + hour_of_day

                   CASE ("holiday")
                      index_hh = months_per_year + days_per_week + 2*hours_per_day + hour_of_day

                END SELECT
                aero_emission_att%time_factor = aero_emission_att%etf(:,index_mm) *                &
                                                aero_emission_att%etf(:,index_dd) *                &
                                                aero_emission_att%etf(:,index_hh+1)
             ENDIF

!
!--          Create a sectional number size distribution for emissions
             ALLOCATE( nsect_emission(1:nbins_aerosol),source_array(nys:nyn,nxl:nxr,1:nbins_aerosol) )
             DO  in = 1, aero_emission_att%ncat

                inn = def_modes%cat_input_to_model(in)
!
!--             Calculate the number concentration (1/m3) of a log-normal size distribution 
!--             following Jacobson (2005): Eq 13.25.
                def_modes%ntot_table = 6.0_wp * def_modes%pm_frac_table(:,inn) / ( pi *            &
                                       ( def_modes%dpg_table )**3 *  EXP( 4.5_wp *                 &
                                       LOG( def_modes%sigmag_table )**2 ) )
!
!--             Sectional size distibution (1/m3) from a log-normal one
                CALL size_distribution( def_modes%ntot_table, def_modes%dpg_table,                 &
                                        def_modes%sigmag_table, nsect_emission )

                source_array = 0.0_wp
                DO  ib = 1, nbins_aerosol
                   source_array(:,:,ib) = aero_emission%def_data(:,:,in) *                         &
                                          aero_emission_att%conversion_factor /                    &
                                          aero_emission_att%rho(in) * nsect_emission(ib) *         &
                                          aero_emission_att%time_factor(in)
                ENDDO
!
!--             Set surface fluxes of aerosol number and mass on horizontal surfaces. Set fluxes
!--             only for either default, land or urban surface.
                IF ( .NOT. land_surface  .AND.  .NOT. urban_surface )  THEN
                   CALL set_flux( surf_def_h(0), aero_emission_att%cc_in2mod,                      &
                                  aero_emission%mass_fracs(in,:), source_array )
                ELSE
                   CALL set_flux( surf_usm_h, aero_emission_att%cc_in2mod,                         &
                                  aero_emission%mass_fracs(in,:), source_array )
                   CALL set_flux( surf_lsm_h, aero_emission_att%cc_in2mod,                         &
                                  aero_emission%mass_fracs(in,:), source_array )
                ENDIF
             ENDDO
!
!--          The next emission update is again after one hour
             next_aero_emission_update = next_aero_emission_update + 3600.0_wp


             DEALLOCATE( nsect_emission, source_array )
!
!--       Pre-processed:
          ELSEIF ( aero_emission_att%lod == 2 )  THEN
!
!--          Obtain time index for current point in time.
             aero_emission_att%tind = MINLOC( ABS( aero_emission_att%time -                        &
                                                   MAX( time_since_reference_point, 0.0_wp ) ),    &
                                              DIM = 1 ) - 1
!
!--          Allocate the data input array always before reading in the data and deallocate after
             ALLOCATE( aero_emission%preproc_data(nys:nyn,nxl:nxr,1:aero_emission_att%ncat),       &
                       source_array(nys:nyn,nxl:nxr,1:nbins_aerosol) )
!
!--          Read in the next time step (get_variable_4d_to_3d_real)
             CALL get_variable( id_salsa, 'aerosol_emission_values', aero_emission%preproc_data,   &
                                aero_emission_att%tind, 0, aero_emission_att%ncat-1,               &
                                nxl, nxr, nys, nyn )
!
!--          Calculate the sources per category and set surface fluxes
             source_array = 0.0_wp
             DO  in = 1, aero_emission_att%ncat
                DO  ib = 1, nbins_aerosol
                   source_array(:,:,ib) = aero_emission%preproc_data(:,:,in) *                     &
                                          aero_emission%num_fracs(in,ib)
                ENDDO
!
!--             Set fluxes only for either default, land and urban surface.
                IF ( .NOT. land_surface  .AND.  .NOT. urban_surface )  THEN
                   CALL set_flux( surf_def_h(0), aero_emission_att%cc_in2mod,                      &
                                  aero_emission%mass_fracs(in,:), source_array )
                ELSE
                   CALL set_flux( surf_usm_h, aero_emission_att%cc_in2mod,                         &
                                  aero_emission%mass_fracs(in,:), source_array )
                   CALL set_flux( surf_lsm_h, aero_emission_att%cc_in2mod,                         &
                                  aero_emission%mass_fracs(in,:), source_array )
                ENDIF
             ENDDO
!
!--          Determine the next emission update
             next_aero_emission_update = aero_emission_att%time(aero_emission_att%tind+2)

             DEALLOCATE( aero_emission%preproc_data, source_array )

          ENDIF
!
!--       Close input file
          CALL close_input_file( id_salsa )
#else
          message_string = 'salsa_emission_mode = "read_from_file", but preprocessor directive ' //&
                           ' __netcdf is not used in compiling!'
          CALL message( 'salsa_emission_setup', 'PA0638', 1, 2, 0, 6, 0 )

#endif
       CASE DEFAULT
          message_string = 'unknown salsa_emission_mode: ' // TRIM( salsa_emission_mode )
          CALL message( 'salsa_emission_setup', 'PA0639', 1, 2, 0, 6, 0 )

    END SELECT

    CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Sets the aerosol flux to aerosol arrays in 2a and 2b.
!------------------------------------------------------------------------------!
    SUBROUTINE set_flux( surface, cc_i_mod, mass_fracs, source_array )

       USE arrays_3d,                                                                              &
           ONLY:  rho_air_zw

       USE surface_mod,                                                                            &
           ONLY:  surf_type

       IMPLICIT NONE

       INTEGER(iwp) ::  i   !< loop index
       INTEGER(iwp) ::  ib  !< loop index
       INTEGER(iwp) ::  ic  !< loop index
       INTEGER(iwp) ::  j   !< loop index
       INTEGER(iwp) ::  k   !< loop index
       INTEGER(iwp) ::  m   !< running index for surface elements

       INTEGER(iwp), DIMENSION(:) ::  cc_i_mod   !< index of chemical component in the input data

       REAL(wp) ::  so4_oc  !< mass fraction between SO4 and OC in 1a

       REAL(wp), DIMENSION(:), INTENT(in) ::  mass_fracs  !< mass fractions of chemical components

       REAL(wp), DIMENSION(nys:nyn,nxl:nxr,1:nbins_aerosol), INTENT(inout) ::  source_array  !<

       TYPE(surf_type), INTENT(inout) :: surface  !< respective surface type

       so4_oc = 0.0_wp

       DO  m = 1, surface%ns
!
!--       Get indices of respective grid point
          i = surface%i(m)
          j = surface%j(m)
          k = surface%k(m)

          DO  ib = 1, nbins_aerosol
             IF ( source_array(j,i,ib) < nclim )  THEN
                source_array(j,i,ib) = 0.0_wp
             ENDIF
!
!--          Set mass fluxes.  First bins include only SO4 and/or OC.
             IF ( ib <= end_subrange_1a )  THEN
!
!--             Both sulphate and organic carbon
                IF ( index_so4 > 0  .AND.  index_oc > 0 )  THEN

                   ic = ( index_so4 - 1 ) * nbins_aerosol + ib
                   so4_oc = mass_fracs(cc_i_mod(1)) / ( mass_fracs(cc_i_mod(1)) +                  &
                                                        mass_fracs(cc_i_mod(2)) )
                   surface%amsws(m,ic) = surface%amsws(m,ic) + so4_oc * source_array(j,i,ib)       &
                                         * api6 * aero(ib)%dmid**3 * arhoh2so4 * rho_air_zw(k-1)
                   aerosol_mass(ic)%source(j,i) = aerosol_mass(ic)%source(j,i) + surface%amsws(m,ic)

                   ic = ( index_oc - 1 ) * nbins_aerosol + ib
                   surface%amsws(m,ic) = surface%amsws(m,ic) + ( 1-so4_oc ) * source_array(j,i,ib) &
                                         * api6 * aero(ib)%dmid**3 * arhooc * rho_air_zw(k-1)
                   aerosol_mass(ic)%source(j,i) = aerosol_mass(ic)%source(j,i) + surface%amsws(m,ic)
!
!--             Only sulphates
                ELSEIF ( index_so4 > 0  .AND.  index_oc < 0 )  THEN
                   ic = ( index_so4 - 1 ) * nbins_aerosol + ib
                   surface%amsws(m,ic) = surface%amsws(m,ic) + source_array(j,i,ib) * api6 *       &
                                         aero(ib)%dmid**3 * arhoh2so4 * rho_air_zw(k-1)
                   aerosol_mass(ic)%source(j,i) = aerosol_mass(ic)%source(j,i) + surface%amsws(m,ic)
!
!--             Only organic carbon
                ELSEIF ( index_so4 < 0  .AND.  index_oc > 0 )  THEN
                   ic = ( index_oc - 1 ) * nbins_aerosol + ib
                   surface%amsws(m,ic) = surface%amsws(m,ic) + source_array(j,i,ib) * api6 *       &
                                         aero(ib)%dmid**3 * arhooc * rho_air_zw(k-1)
                   aerosol_mass(ic)%source(j,i) = aerosol_mass(ic)%source(j,i) + surface%amsws(m,ic)
                ENDIF

             ELSE
!
!--             Sulphate
                IF ( index_so4 > 0 )  THEN
                   ic = cc_i_mod(1)
                   CALL set_mass_flux( surface, m, ib, index_so4, mass_fracs(ic), arhoh2so4,       &
                                       source_array(j,i,ib) )
                ENDIF
!
!--             Organic carbon
                IF ( index_oc > 0 )  THEN
                   ic = cc_i_mod(2)
                   CALL set_mass_flux( surface, m, ib, index_oc, mass_fracs(ic),arhooc,            &
                                       source_array(j,i,ib) )
                ENDIF
!
!--             Black carbon
                IF ( index_bc > 0 )  THEN
                   ic = cc_i_mod(3)
                   CALL set_mass_flux( surface, m, ib, index_bc, mass_fracs(ic), arhobc,           &
                                       source_array(j,i,ib) )
                ENDIF
!
!--             Dust
                IF ( index_du > 0 )  THEN
                   ic = cc_i_mod(4)
                   CALL set_mass_flux( surface, m, ib, index_du, mass_fracs(ic), arhodu,           &
                                       source_array(j,i,ib) )
                ENDIF
!
!--             Sea salt
                IF ( index_ss > 0 )  THEN
                   ic = cc_i_mod(5)
                   CALL set_mass_flux( surface, m, ib, index_ss, mass_fracs(ic), arhoss,           &
                                       source_array(j,i,ib) )
                ENDIF
!
!--             Nitric acid
                IF ( index_no > 0 )  THEN
                    ic = cc_i_mod(6)
                   CALL set_mass_flux( surface, m, ib, index_no, mass_fracs(ic), arhohno3,         &
                                       source_array(j,i,ib) )
                ENDIF
!
!--             Ammonia
                IF ( index_nh > 0 )  THEN
                    ic = cc_i_mod(7)
                   CALL set_mass_flux( surface, m, ib, index_nh, mass_fracs(ic), arhonh3,          &
                                       source_array(j,i,ib) )
                ENDIF

             ENDIF
!
!--          Save number fluxes in the end
             surface%answs(m,ib) = surface%answs(m,ib) + source_array(j,i,ib) * rho_air_zw(k-1)
             aerosol_number(ib)%source(j,i) = surface%answs(m,ib)

          ENDDO  ! ib
       ENDDO  ! m

    END SUBROUTINE set_flux

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Sets the mass emissions to aerosol arrays in 2a and 2b.
!------------------------------------------------------------------------------!
    SUBROUTINE set_mass_flux( surface, surf_num, ib, ispec, mass_frac, prho, nsource )

       USE arrays_3d,                                                                              &
           ONLY:  rho_air_zw

       USE surface_mod,                                                                            &
           ONLY:  surf_type

       IMPLICIT NONE

       INTEGER(iwp) ::  i   !< loop index
       INTEGER(iwp) ::  j   !< loop index
       INTEGER(iwp) ::  k   !< loop index
       INTEGER(iwp) ::  ic  !< loop index

       INTEGER(iwp), INTENT(in) :: ib        !< Aerosol size bin index
       INTEGER(iwp), INTENT(in) :: ispec     !< Aerosol species index
       INTEGER(iwp), INTENT(in) :: surf_num  !< index surface elements

       REAL(wp), INTENT(in) ::  mass_frac    !< mass fraction of a chemical compound in all bins
       REAL(wp), INTENT(in) ::  nsource      !< number source (#/m2/s)
       REAL(wp), INTENT(in) ::  prho         !< Aerosol density

       TYPE(surf_type), INTENT(inout) ::  surface  !< respective surface type
!
!--    Get indices of respective grid point
       i = surface%i(surf_num)
       j = surface%j(surf_num)
       k = surface%k(surf_num)
!
!--    Subrange 2a:
       ic = ( ispec - 1 ) * nbins_aerosol + ib
       surface%amsws(surf_num,ic) = surface%amsws(surf_num,ic) + mass_frac * nsource *             &
                                    aero(ib)%core * prho * rho_air_zw(k-1)
       aerosol_mass(ic)%source(j,i) = surface%amsws(surf_num,ic)

    END SUBROUTINE set_mass_flux

 END SUBROUTINE salsa_emission_setup

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Sets the gaseous fluxes
!------------------------------------------------------------------------------!
 SUBROUTINE salsa_gas_emission_setup( init )

    USE netcdf_data_input_mod,                                                                     &
        ONLY:  check_existence, close_input_file, get_attribute, get_variable,                     &
               inquire_num_variables, inquire_variable_names,                                      &
               get_dimension_length, open_read_file

    USE palm_date_time_mod,                                                                        &
        ONLY:  days_per_week, get_date_time, hours_per_day, months_per_year, seconds_per_hour

    USE surface_mod,                                                                               &
        ONLY:  surf_def_h, surf_lsm_h, surf_usm_h

    IMPLICIT NONE

    CHARACTER(LEN=80) ::  daytype = 'workday'  !< default day type
    CHARACTER(LEN=25) ::  in_name              !< name of a gas in the input file

    CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE ::  var_names   !<  variable names in input data


    INTEGER(iwp) ::  day_of_month   !< day of the month
    INTEGER(iwp) ::  day_of_week    !< day of the week
    INTEGER(iwp) ::  day_of_year    !< day of the year
    INTEGER(iwp) ::  hour_of_day    !< hour of the day
    INTEGER(iwp) ::  id_chem        !< NetCDF id of chemistry emission file
    INTEGER(iwp) ::  i              !< loop index
    INTEGER(iwp) ::  ig             !< loop index
    INTEGER(iwp) ::  in             !< running index for emission categories
    INTEGER(iwp) ::  index_dd       !< index day
    INTEGER(iwp) ::  index_hh       !< index hour
    INTEGER(iwp) ::  index_mm       !< index month
    INTEGER(iwp) ::  j              !< loop index
    INTEGER(iwp) ::  month_of_year  !< month of the year
    INTEGER(iwp) ::  num_vars       !< number of variables

    LOGICAL  ::  netcdf_extend = .FALSE.  !< NetCDF input file exists

    LOGICAL, INTENT(in) ::  init          !< if .TRUE. --> initialisation call

    REAL(wp) ::  second_of_day    !< second of the day

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  time_factor  !< emission time factor

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  dum_var_3d  !<

    REAL(wp), DIMENSION(:,:,:,:,:), ALLOCATABLE ::  dum_var_5d  !<

!
!-- Reset surface fluxes
    surf_def_h(0)%gtsws = 0.0_wp
    surf_lsm_h%gtsws = 0.0_wp
    surf_usm_h%gtsws = 0.0_wp

#if defined( __netcdf )
!
!-- Check existence of PIDS_CHEM file
    INQUIRE( FILE = 'PIDS_CHEM' // TRIM( coupling_char ), EXIST = netcdf_extend )
    IF ( .NOT. netcdf_extend )  THEN
       message_string = 'Input file PIDS_CHEM' //  TRIM( coupling_char ) // ' missing!'
       CALL message( 'salsa_gas_emission_setup', 'PA0640', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Open file in read-only mode
    CALL open_read_file( 'PIDS_CHEM' // TRIM( coupling_char ), id_chem )

    IF ( init )  THEN
!
!--    Read the index and name of chemical components
       CALL get_dimension_length( id_chem, chem_emission_att%n_emiss_species, 'nspecies' )
       ALLOCATE( chem_emission_att%species_index(1:chem_emission_att%n_emiss_species) )
       CALL get_variable( id_chem, 'emission_index', chem_emission_att%species_index )
       CALL get_variable( id_chem, 'emission_name', chem_emission_att%species_name,                &
                          chem_emission_att%n_emiss_species )
!
!--    Allocate emission data
       ALLOCATE( chem_emission(1:chem_emission_att%n_emiss_species) )
!
!--    Find the corresponding indices in the model
       emission_index_chem = 0
       DO  ig = 1, chem_emission_att%n_emiss_species
          in_name = chem_emission_att%species_name(ig)
          SELECT CASE ( TRIM( in_name ) )
             CASE ( 'H2SO4', 'h2so4' )
                emission_index_chem(1) = ig
             CASE ( 'HNO3', 'hno3' )
                emission_index_chem(2) = ig
             CASE ( 'NH3', 'nh3' )
                emission_index_chem(3) = ig
             CASE ( 'OCNV', 'ocnv' )
                emission_index_chem(4) = ig
             CASE ( 'OCSV', 'ocsv' )
                emission_index_chem(5) = ig
          END SELECT
       ENDDO
!
!--    Inquire the fill value
       CALL get_attribute( id_chem, '_FillValue', aero_emission%fill, .FALSE., 'emission_values' )
!
!--    Inquire units of emissions
       CALL get_attribute( id_chem, 'units', chem_emission_att%units, .FALSE., 'emission_values' )
!
!--    Inquire the level of detail (lod)
       CALL get_attribute( id_chem, 'lod', lod_gas_emissions, .FALSE., 'emission_values' )
!
!--    Variable names
       CALL inquire_num_variables( id_chem, num_vars )
       ALLOCATE( var_names(1:num_vars) )
       CALL inquire_variable_names( id_chem, var_names )
!
!--    Default mode: as total emissions per year
       IF ( lod_gas_emissions == 1 )  THEN

!
!--       Get number of emission categories and allocate emission arrays
          CALL get_dimension_length( id_chem, chem_emission_att%ncat, 'ncat' )
          ALLOCATE( chem_emission_att%cat_index(1:chem_emission_att%ncat),                         &
                    time_factor(1:chem_emission_att%ncat) )
!
!--       Get emission category names and indices
          CALL get_variable( id_chem, 'emission_category_name', chem_emission_att%cat_name,        &
                             chem_emission_att%ncat)
          CALL get_variable( id_chem, 'emission_category_index', chem_emission_att%cat_index )
!
!--       Emission time factors: Find check whether emission time factors are given for each hour
!--       of year OR based on month, day and hour
!
!--       For each hour of year:
          IF ( check_existence( var_names, 'nhoursyear' ) )  THEN
             CALL get_dimension_length( id_chem, chem_emission_att%nhoursyear, 'nhoursyear' )
             ALLOCATE( chem_emission_att%hourly_emis_time_factor(1:chem_emission_att%ncat,         &
                                                                 1:chem_emission_att%nhoursyear) )
             CALL get_variable( id_chem, 'emission_time_factors',                                  &
                                chem_emission_att%hourly_emis_time_factor,                         &
                                0, chem_emission_att%nhoursyear-1, 0, chem_emission_att%ncat-1 )
!
!--       Based on the month, day and hour:
          ELSEIF ( check_existence( var_names, 'nmonthdayhour' ) )  THEN
             CALL get_dimension_length( id_chem, chem_emission_att%nmonthdayhour, 'nmonthdayhour' )
             ALLOCATE( chem_emission_att%mdh_emis_time_factor(1:chem_emission_att%ncat,            &
                                                              1:chem_emission_att%nmonthdayhour) )
             CALL get_variable( id_chem, 'emission_time_factors',                                  &
                                chem_emission_att%mdh_emis_time_factor,                            &
                                0, chem_emission_att%nmonthdayhour-1, 0, chem_emission_att%ncat-1 )
          ELSE
             message_string = 'emission_time_factors should be given for each nhoursyear OR ' //   &
                              'nmonthdayhour'
             CALL message( 'salsa_gas_emission_setup','PA0641', 1, 2, 0, 6, 0 )
          ENDIF
!
!--       Next emission update
          CALL get_date_time( time_since_reference_point, second_of_day=second_of_day )
          next_gas_emission_update = MOD( second_of_day, seconds_per_hour ) !- seconds_per_hour
!
!--       Allocate and read surface emission data (in total PM) (NOTE that "preprocessed" input data
!--       array is applied now here)
          ALLOCATE( dum_var_5d(1,nys:nyn,nxl:nxr,1:chem_emission_att%n_emiss_species,              &
                               1:chem_emission_att%ncat) )
          CALL get_variable( id_chem, 'emission_values', dum_var_5d, 0, chem_emission_att%ncat-1,  &
                             0, chem_emission_att%n_emiss_species-1, nxl, nxr, nys, nyn, 0, 0 )
          DO  ig = 1, chem_emission_att%n_emiss_species
             ALLOCATE( chem_emission(ig)%default_emission_data(nys:nyn,nxl:nxr,                    &
                                                               1:chem_emission_att%ncat) )
             DO  in = 1, chem_emission_att%ncat
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      chem_emission(ig)%default_emission_data(j,i,in) = dum_var_5d(1,j,i,ig,in)
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
          DEALLOCATE( dum_var_5d )
!
!--    Pre-processed mode:
       ELSEIF ( lod_gas_emissions == 2 )  THEN
!
!--       Number of time steps in the emission data
          CALL get_dimension_length( id_chem, chem_emission_att%dt_emission, 'time' )
!
!--       Allocate and read time
          ALLOCATE( gas_emission_time(1:chem_emission_att%dt_emission) )
          CALL get_variable( id_chem, 'time', gas_emission_time )
       ELSE
          message_string = 'Unknown lod for emission_values.'
          CALL message( 'salsa_gas_emission_setup','PA0642', 1, 2, 0, 6, 0 )
       ENDIF  ! lod

    ENDIF  ! init
!
!-- Define and set current emission values:

    IF ( lod_gas_emissions == 1 )  THEN
!
!--    Emission time factors for each emission category at current time step
       IF ( chem_emission_att%nhoursyear > chem_emission_att%nmonthdayhour )  THEN
!
!--       Get the index of the current hour
          CALL get_date_time( time_since_reference_point, &
                              day_of_year=day_of_year, hour=hour_of_day )
          index_hh = ( day_of_year - 1_iwp ) * hours_per_day + hour_of_day
          IF ( .NOT. ALLOCATED( time_factor ) )  ALLOCATE( time_factor(1:chem_emission_att%ncat) )
          time_factor = 0.0_wp
          time_factor = chem_emission_att%hourly_emis_time_factor(:,index_hh+1)

       ELSEIF ( chem_emission_att%nhoursyear < chem_emission_att%nmonthdayhour )  THEN
!
!--       Get the index of current hour (index_hh) (TODO: Now "workday" is always assumed.
!--       Needs to be calculated.)
          CALL get_date_time( time_since_reference_point, &
                              month=month_of_year,        &
                              day=day_of_month,           &
                              hour=hour_of_day,           &
                              day_of_week=day_of_week     )
          index_mm = month_of_year
          index_dd = months_per_year + day_of_week
          SELECT CASE( TRIM( daytype ) )

             CASE ("workday")
                index_hh = months_per_year + days_per_week + hour_of_day

             CASE ("weekend")
                index_hh = months_per_year + days_per_week + hours_per_day + hour_of_day

             CASE ("holiday")
                index_hh = months_per_year + days_per_week + 2*hours_per_day + hour_of_day

          END SELECT
          time_factor = chem_emission_att%mdh_emis_time_factor(:,index_mm) *                       &
                        chem_emission_att%mdh_emis_time_factor(:,index_dd) *                       &
                        chem_emission_att%mdh_emis_time_factor(:,index_hh+1)
       ENDIF
!
!--    Set gas emissions for each emission category
       ALLOCATE( dum_var_3d(nys:nyn,nxl:nxr,1:chem_emission_att%n_emiss_species) )

       DO  in = 1, chem_emission_att%ncat
          DO  ig = 1, chem_emission_att%n_emiss_species
             dum_var_3d(:,:,ig) = chem_emission(ig)%default_emission_data(:,:,in)
          ENDDO
!
!--       Set surface fluxes only for either default, land or urban surface
          IF ( .NOT. land_surface  .AND.  .NOT. urban_surface )  THEN
             CALL set_gas_flux( surf_def_h(0), emission_index_chem, chem_emission_att%units,    &
                                dum_var_3d, time_factor(in) )
          ELSE
             CALL set_gas_flux( surf_usm_h, emission_index_chem, chem_emission_att%units,       &
                                dum_var_3d, time_factor(in) )
             CALL set_gas_flux( surf_lsm_h, emission_index_chem, chem_emission_att%units,       &
                                dum_var_3d, time_factor(in) )
          ENDIF
       ENDDO
       DEALLOCATE( dum_var_3d )
!
!--    The next emission update is again after one hour
       next_gas_emission_update = next_gas_emission_update + 3600.0_wp

    ELSEIF ( lod_gas_emissions == 2 )  THEN
!
!--    Obtain time index for current point in time.
       chem_emission_att%i_hour = MINLOC( ABS( gas_emission_time -                                 &
                                          MAX( time_since_reference_point, 0.0_wp ) ), DIM = 1 ) - 1
!
!--    Allocate the data input array always before reading in the data and deallocate after (NOTE 
!--    that "preprocessed" input data array is applied now here)
       ALLOCATE( dum_var_5d(1,1,nys:nyn,nxl:nxr,1:chem_emission_att%n_emiss_species) )
!
!--    Read in the next time step
       CALL get_variable( id_chem, 'emission_values', dum_var_5d,                                  &
                          0, chem_emission_att%n_emiss_species-1, nxl, nxr, nys, nyn, 0, 0,        &
                          chem_emission_att%i_hour, chem_emission_att%i_hour )
!
!--    Set surface fluxes only for either default, land or urban surface
       IF ( .NOT. land_surface  .AND.  .NOT. urban_surface )  THEN
          CALL set_gas_flux( surf_def_h(0), emission_index_chem, chem_emission_att%units,          &
                             dum_var_5d(1,1,:,:,:) )
       ELSE
          CALL set_gas_flux( surf_usm_h, emission_index_chem, chem_emission_att%units,             &
                             dum_var_5d(1,1,:,:,:) )
          CALL set_gas_flux( surf_lsm_h, emission_index_chem, chem_emission_att%units,             &
                             dum_var_5d(1,1,:,:,:) )
       ENDIF
       DEALLOCATE ( dum_var_5d )
!
!--    Determine the next emission update
       next_gas_emission_update = gas_emission_time(chem_emission_att%i_hour+2)

    ENDIF
!
!-- Close input file
    CALL close_input_file( id_chem )

#else
    message_string = 'salsa_emission_mode = "read_from_file", but preprocessor directive ' //   &
                     ' __netcdf is not used in compiling!'
    CALL message( 'salsa_gas_emission_setup', 'PA0643', 1, 2, 0, 6, 0 )

#endif

    CONTAINS
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Set gas fluxes for selected type of surfaces
!------------------------------------------------------------------------------!
    SUBROUTINE set_gas_flux( surface, cc_i_mod, unit, source_array, time_fac )

       USE arrays_3d,                                                                              &
           ONLY: dzw, hyp, pt, rho_air_zw

       USE grid_variables,                                                                         &
           ONLY:  dx, dy

       USE surface_mod,                                                                            &
           ONLY:  surf_type

       IMPLICIT NONE

       CHARACTER(LEN=*), INTENT(in) ::  unit  !< flux unit in the input file

       INTEGER(iwp) ::  ig  !< running index for gases
       INTEGER(iwp) ::  i   !< loop index
       INTEGER(iwp) ::  j   !< loop index
       INTEGER(iwp) ::  k   !< loop index
       INTEGER(iwp) ::  m   !< running index for surface elements

       INTEGER(iwp), DIMENSION(:) ::  cc_i_mod   !< index of different gases in the input data

       LOGICAL ::  use_time_fac  !< .TRUE. is time_fac present

       REAL(wp), OPTIONAL ::  time_fac  !< emission time factor

       REAL(wp), DIMENSION(ngases_salsa) ::  conv     !< unit conversion factor

       REAL(wp), DIMENSION(nys:nyn,nxl:nxr,1:chem_emission_att%n_emiss_species), INTENT(in) ::  source_array  !<

       TYPE(surf_type), INTENT(inout) :: surface  !< respective surface type

       conv = 1.0_wp
       use_time_fac = PRESENT( time_fac )

       DO  m = 1, surface%ns
!
!--       Get indices of respective grid point
          i = surface%i(m)
          j = surface%j(m)
          k = surface%k(m)
!
!--       Unit conversion factor: convert to SI units (#/m2/s)
          SELECT CASE ( TRIM( unit ) )
             CASE ( 'kg/m2/yr' )
                conv(1) = avo / ( amh2so4 * 3600.0_wp )
                conv(2) = avo / ( amhno3 * 3600.0_wp )
                conv(3) = avo / ( amnh3 * 3600.0_wp )
                conv(4) = avo / ( amoc * 3600.0_wp )
                conv(5) = avo / ( amoc * 3600.0_wp )
             CASE ( 'g/m2/yr' )
                conv(1) = avo / ( amh2so4 * 3.6E+6_wp )
                conv(2) = avo / ( amhno3 * 3.6E+6_wp )
                conv(3) = avo / ( amnh3 * 3.6E+6_wp )
                conv(4) = avo / ( amoc * 3.6E+6_wp )
                conv(5) = avo / ( amoc * 3.6E+6_wp )
             CASE ( 'g/m2/s' )
                conv(1) = avo / ( amh2so4 * 1000.0_wp )
                conv(2) = avo / ( amhno3 * 1000.0_wp )
                conv(3) = avo / ( amnh3 * 1000.0_wp )
                conv(4) = avo / ( amoc * 1000.0_wp )
                conv(5) = avo / ( amoc * 1000.0_wp )
             CASE ( '#/m2/s' )
                conv = 1.0_wp
             CASE ( 'ppm/m2/s' )
                conv = for_ppm_to_nconc * hyp(k) / pt(k,j,i) * ( 1.0E5_wp / hyp(k) )**0.286_wp *   &
                       dx * dy * dzw(k)
             CASE ( 'mumol/m2/s' )
                conv = 1.0E-6_wp * avo
             CASE DEFAULT
                message_string = 'unknown unit for gas emissions: ' // TRIM( chem_emission_att%units )
                CALL message( 'set_gas_flux','PA0644', 1, 2, 0, 6, 0 )

          END SELECT

          DO  ig = 1, ngases_salsa
             IF ( use_time_fac )  THEN
                surface%gtsws(m,ig) = surface%gtsws(m,ig) + rho_air_zw(k-1) * conv(ig) * time_fac  &
                                      * MAX( 0.0_wp, source_array(j,i,cc_i_mod(ig) ) )
             ELSE
                surface%gtsws(m,ig) = surface%gtsws(m,ig) + rho_air_zw(k-1) * conv(ig)             &
                                      * MAX( 0.0_wp, source_array(j,i,cc_i_mod(ig) ) )
             ENDIF
          ENDDO  ! ig

       ENDDO  ! m

    END SUBROUTINE set_gas_flux

 END SUBROUTINE salsa_gas_emission_setup

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check data output for salsa.
!------------------------------------------------------------------------------!
 SUBROUTINE salsa_check_data_output( var, unit )

    IMPLICIT NONE

    CHARACTER(LEN=*) ::  unit     !<
    CHARACTER(LEN=*) ::  var      !<

    INTEGER(iwp) ::  char_to_int   !< for converting character to integer

    IF ( var(1:6) /= 'salsa_' )  THEN
       unit = 'illegal'
       RETURN
    ENDIF
!
!-- Treat bin-specific outputs separately
    IF ( var(7:11) ==  'N_bin' )  THEN
       READ( var(12:),* ) char_to_int
       IF ( char_to_int >= 1  .AND. char_to_int <= SUM( nbin ) )  THEN
          unit = '#/m3'
       ELSE
          unit = 'illegal'
          RETURN
       ENDIF

    ELSEIF ( var(7:11) ==  'm_bin' )  THEN
       READ( var(12:),* ) char_to_int
       IF ( char_to_int >= 1  .AND. char_to_int <= SUM( nbin ) )  THEN
          unit = 'kg/m3'
       ELSE
          unit = 'illegal'
          RETURN
       ENDIF

    ELSEIF ( var(7:11) == 's_H2O' )  THEN
       IF ( .NOT. advect_particle_water )  THEN
          message_string = 'to output s_H2O/s_H2O_av requires that advect_particle_water = .T.'
          CALL message( 'check_parameters', 'PA0707', 1, 2, 0, 6, 0 )
       ENDIF

    ELSE
       SELECT CASE ( TRIM( var(7:) ) )

          CASE ( 'g_H2SO4', 'g_HNO3', 'g_NH3', 'g_OCNV',  'g_OCSV' )
             IF (  air_chemistry )  THEN
                message_string = 'gases are imported from the chemistry module and thus output '// &
                                 'of "' // TRIM( var ) // '" is not allowed'
                CALL message( 'check_parameters', 'PA0653', 1, 2, 0, 6, 0 )
             ENDIF
             unit = '#/m3'

          CASE ( 'LDSA' )
             unit = 'mum2/cm3'

          CASE ( 'PM0.1', 'PM2.5', 'PM10', 's_BC', 's_DU', 's_H2O', 's_NH', 's_NO', 's_OC',        &
                 's_SO4', 's_SS' )
             unit = 'kg/m3'

          CASE ( 'N_UFP', 'Ntot' )
             unit = '#/m3'

          CASE DEFAULT
             unit = 'illegal'

       END SELECT
    ENDIF

 END SUBROUTINE salsa_check_data_output

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check profile data output for salsa. Currently only for diagnostic variables
!> Ntot, N_UFP, PM0.1, PM2.5, PM10 and LDSA
!------------------------------------------------------------------------------!
 SUBROUTINE salsa_check_data_output_pr( var, var_count, unit, dopr_unit )

    USE arrays_3d,                                                                                 &
        ONLY: zu

    USE profil_parameter,                                                                          &
        ONLY:  dopr_index

    USE statistics,                                                                                &
        ONLY:  hom, pr_palm, statistic_regions

    IMPLICIT NONE

    CHARACTER(LEN=*) ::  dopr_unit  !<
    CHARACTER(LEN=*) ::  unit       !<
    CHARACTER(LEN=*) ::  var        !<

    INTEGER(iwp) ::  var_count     !<

    IF ( var(1:6) /= 'salsa_' )  THEN
       unit = 'illegal'
       RETURN
    ENDIF

    SELECT CASE ( TRIM( var(7:) ) )

       CASE( 'LDSA' )
          salsa_pr_count = salsa_pr_count + 1
          salsa_pr_index(salsa_pr_count) = 1
          dopr_index(var_count) = pr_palm + salsa_pr_count
          dopr_unit = 'mum2/cm3'
          unit = dopr_unit
          hom(:,2,dopr_index(var_count),:) = SPREAD( zu, 2, statistic_regions+1 )

       CASE( 'N_UFP' )
          salsa_pr_count = salsa_pr_count + 1
          salsa_pr_index(salsa_pr_count) = 2
          dopr_index(var_count) = pr_palm + salsa_pr_count
          dopr_unit = '#/m3'
          unit = dopr_unit
          hom(:,2,dopr_index(var_count),:) = SPREAD( zu, 2, statistic_regions+1 )

       CASE( 'Ntot' )
          salsa_pr_count = salsa_pr_count + 1
          salsa_pr_index(salsa_pr_count) = 3
          dopr_index(var_count) = pr_palm + salsa_pr_count
          dopr_unit = '#/m3'
          unit = dopr_unit
          hom(:,2,dopr_index(var_count),:) = SPREAD( zu, 2, statistic_regions+1 )

       CASE( 'PM0.1' )
          salsa_pr_count = salsa_pr_count + 1
          salsa_pr_index(salsa_pr_count) = 4
          dopr_index(var_count) = pr_palm + salsa_pr_count
          dopr_unit = 'kg/m3'
          unit = dopr_unit
          hom(:,2,dopr_index(var_count),:) = SPREAD( zu, 2, statistic_regions+1 )

       CASE( 'PM2.5' )
          salsa_pr_count = salsa_pr_count + 1
          salsa_pr_index(salsa_pr_count) = 5
          dopr_index(var_count) = pr_palm + salsa_pr_count
          dopr_unit = 'kg/m3'
          unit = dopr_unit
          hom(:,2,dopr_index(var_count),:) = SPREAD( zu, 2, statistic_regions+1 )

       CASE( 'PM10' )
          salsa_pr_count = salsa_pr_count + 1
          salsa_pr_index(salsa_pr_count) = 6
          dopr_index(var_count) = pr_palm + salsa_pr_count
          dopr_unit = 'kg/m3'
          unit = dopr_unit
          hom(:,2,dopr_index(var_count),:) = SPREAD( zu, 2, statistic_regions+1 )

       CASE DEFAULT
          unit = 'illegal'

    END SELECT


 END SUBROUTINE salsa_check_data_output_pr

!-------------------------------------------------------------------------------!
!> Description:
!> Calculation of horizontally averaged profiles for salsa.
!-------------------------------------------------------------------------------!
 SUBROUTINE salsa_statistics( mode, sr, tn )

    USE control_parameters,                                                                        &
        ONLY:  max_pr_user

    USE chem_modules,                                                                              &
        ONLY:  max_pr_cs

    USE statistics,                                                                                &
        ONLY:  pr_palm, rmask, sums_l

    IMPLICIT NONE

    CHARACTER(LEN=*) ::  mode  !<

    INTEGER(iwp) ::  i    !< loop index
    INTEGER(iwp) ::  ib   !< loop index
    INTEGER(iwp) ::  ic   !< loop index
    INTEGER(iwp) ::  ii   !< loop index
    INTEGER(iwp) ::  ind  !< index in the statistical output
    INTEGER(iwp) ::  j    !< loop index
    INTEGER(iwp) ::  k    !< loop index
    INTEGER(iwp) ::  sr   !< statistical region
    INTEGER(iwp) ::  tn   !< thread number

    REAL(wp) ::  df        !< For calculating LDSA: fraction of particles depositing in the alveolar
                           !< (or tracheobronchial) region of the lung. Depends on the particle size
    REAL(wp) ::  mean_d    !< Particle diameter in micrometres
    REAL(wp) ::  temp_bin  !< temporary variable

    IF ( mode == 'profiles' )  THEN
       !$OMP DO
       DO  ii = 1, salsa_pr_count

          ind = pr_palm + max_pr_user + max_pr_cs + ii

          SELECT CASE( salsa_pr_index(ii) )

             CASE( 1 )  ! LDSA
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         temp_bin = 0.0_wp
                         DO  ib = 1, nbins_aerosol
   !
   !--                      Diameter in micrometres
                            mean_d = 1.0E+6_wp * ra_dry(k,j,i,ib) * 2.0_wp
   !
   !--                      Deposition factor: alveolar
                            df = ( 0.01555_wp / mean_d ) * ( EXP( -0.416_wp * ( LOG( mean_d ) +    &
                                   2.84_wp )**2 ) + 19.11_wp * EXP( -0.482_wp * ( LOG( mean_d ) -  &
                                   1.362_wp )**2 ) )
   !
   !--                      Lung-deposited surface area LDSA (units mum2/cm3)
                            temp_bin = temp_bin + pi * mean_d**2 * df * 1.0E-6_wp *                &
                                       aerosol_number(ib)%conc(k,j,i)
                         ENDDO
                         sums_l(k,ind,tn) = sums_l(k,ind,tn) + temp_bin * rmask(j,i,sr)  *         &
                                           MERGE( 1.0_wp, 0.0_wp,                                  &
                                           BTEST( wall_flags_total_0(k,j,i), 22 ) )
                      ENDDO
                   ENDDO
                ENDDO

             CASE( 2 )  ! N_UFP
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         temp_bin = 0.0_wp
                         DO  ib = 1, nbins_aerosol
                            IF ( 2.0_wp * ra_dry(k,j,i,ib) <= 0.1E-6_wp )                          &
                               temp_bin = temp_bin + aerosol_number(ib)%conc(k,j,i)
                         ENDDO
                         sums_l(k,ind,tn) = sums_l(k,ind,tn) + temp_bin * rmask(j,i,sr)  *         &
                                           MERGE( 1.0_wp, 0.0_wp,                                  &
                                           BTEST( wall_flags_total_0(k,j,i), 22 ) )
                      ENDDO
                   ENDDO
                ENDDO

             CASE( 3 )  ! Ntot
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         temp_bin = 0.0_wp
                         DO  ib = 1, nbins_aerosol
                            temp_bin = temp_bin + aerosol_number(ib)%conc(k,j,i)
                         ENDDO
                         sums_l(k,ind,tn) = sums_l(k,ind,tn) + temp_bin * rmask(j,i,sr)  *         &
                                           MERGE( 1.0_wp, 0.0_wp,                                  &
                                           BTEST( wall_flags_total_0(k,j,i), 22 ) )
                      ENDDO
                   ENDDO
                ENDDO

             CASE( 4 )  ! PM0.1
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         temp_bin = 0.0_wp
                         DO  ib = 1, nbins_aerosol
                            IF ( 2.0_wp * ra_dry(k,j,i,ib) <= 0.1E-6_wp )  THEN
                               DO  ic = ib, nbins_aerosol * ncc, nbins_aerosol
                                  temp_bin = temp_bin + aerosol_mass(ic)%conc(k,j,i)
                               ENDDO
                            ENDIF
                         ENDDO
                         sums_l(k,ind,tn) = sums_l(k,ind,tn) + temp_bin * rmask(j,i,sr)  *         &
                                           MERGE( 1.0_wp, 0.0_wp,                                  &
                                           BTEST( wall_flags_total_0(k,j,i), 22 ) )
                      ENDDO
                   ENDDO
                ENDDO

             CASE( 5 )  ! PM2.5
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         temp_bin = 0.0_wp
                         DO  ib = 1, nbins_aerosol
                            IF ( 2.0_wp * ra_dry(k,j,i,ib) <= 2.5E-6_wp )  THEN
                               DO  ic = ib, nbins_aerosol * ncc, nbins_aerosol
                                  temp_bin = temp_bin + aerosol_mass(ic)%conc(k,j,i)
                               ENDDO
                            ENDIF
                         ENDDO
                         sums_l(k,ind,tn) = sums_l(k,ind,tn) + temp_bin * rmask(j,i,sr)  *         &
                                           MERGE( 1.0_wp, 0.0_wp,                                  &
                                           BTEST( wall_flags_total_0(k,j,i), 22 ) )
                      ENDDO
                   ENDDO
                ENDDO

             CASE( 6 )  ! PM10
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         temp_bin = 0.0_wp
                         DO  ib = 1, nbins_aerosol
                            IF ( 2.0_wp * ra_dry(k,j,i,ib) <= 10.0E-6_wp )  THEN
                               DO  ic = ib, nbins_aerosol * ncc, nbins_aerosol
                                  temp_bin = temp_bin + aerosol_mass(ic)%conc(k,j,i)
                               ENDDO
                            ENDIF
                         ENDDO
                         sums_l(k,ind,tn) = sums_l(k,ind,tn) + temp_bin * rmask(j,i,sr)  *         &
                                           MERGE( 1.0_wp, 0.0_wp,                                  &
                                           BTEST( wall_flags_total_0(k,j,i), 22 ) )
                      ENDDO
                   ENDDO
                ENDDO

          END SELECT
       ENDDO

    ELSEIF ( mode == 'time_series' )  THEN
!
!--    TODO
    ENDIF

 END SUBROUTINE salsa_statistics


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine for averaging 3D data
!------------------------------------------------------------------------------!
 SUBROUTINE salsa_3d_data_averaging( mode, variable )

    USE control_parameters,                                                                        &
        ONLY:  average_count_3d

    IMPLICIT NONE

    CHARACTER(LEN=*)  ::  mode       !<
    CHARACTER(LEN=10) ::  vari       !<
    CHARACTER(LEN=*)  ::  variable   !<

    INTEGER(iwp) ::  char_to_int  !< for converting character to integer
    INTEGER(iwp) ::  found_index  !<
    INTEGER(iwp) ::  i            !<
    INTEGER(iwp) ::  ib           !<
    INTEGER(iwp) ::  ic           !<
    INTEGER(iwp) ::  j            !<
    INTEGER(iwp) ::  k            !<

    REAL(wp) ::  df       !< For calculating LDSA: fraction of particles depositing in the alveolar
                          !< (or tracheobronchial) region of the lung. Depends on the particle size
    REAL(wp) ::  mean_d   !< Particle diameter in micrometres
    REAL(wp) ::  temp_bin !< temporary variable

    REAL(wp), DIMENSION(:,:,:), POINTER ::  to_be_resorted  !< points to selected output variable

    temp_bin = 0.0_wp

    IF ( mode == 'allocate' )  THEN

       IF ( variable(7:11) ==  'N_bin' )  THEN
          IF ( .NOT. ALLOCATED( nbins_av ) )  THEN
             ALLOCATE( nbins_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg,nbins_aerosol) )
          ENDIF
          nbins_av = 0.0_wp

       ELSEIF ( variable(7:11) ==  'm_bin' )  THEN
          IF ( .NOT. ALLOCATED( mbins_av ) )  THEN
             ALLOCATE( mbins_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg,nbins_aerosol) )
          ENDIF
          mbins_av = 0.0_wp

       ELSE

          SELECT CASE ( TRIM( variable(7:) ) )

             CASE ( 'g_H2SO4', 'g_HNO3', 'g_NH3', 'g_OCNV', 'g_OCSV' )
                IF ( .NOT. ALLOCATED( salsa_gases_av ) )  THEN
                   ALLOCATE( salsa_gases_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg,ngases_salsa) )
                ENDIF
                salsa_gases_av = 0.0_wp

             CASE ( 'LDSA' )
                IF ( .NOT. ALLOCATED( ldsa_av ) )  THEN
                   ALLOCATE( ldsa_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                ldsa_av = 0.0_wp

             CASE ( 'N_UFP' )
                IF ( .NOT. ALLOCATED( nufp_av ) )  THEN
                   ALLOCATE( nufp_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                nufp_av = 0.0_wp

             CASE ( 'Ntot' )
                IF ( .NOT. ALLOCATED( ntot_av ) )  THEN
                   ALLOCATE( ntot_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                ntot_av = 0.0_wp

             CASE ( 'PM0.1' )
                IF ( .NOT. ALLOCATED( pm01_av ) )  THEN
                   ALLOCATE( pm01_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                pm01_av = 0.0_wp

             CASE ( 'PM2.5' )
                IF ( .NOT. ALLOCATED( pm25_av ) )  THEN
                   ALLOCATE( pm25_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                pm25_av = 0.0_wp

             CASE ( 'PM10' )
                IF ( .NOT. ALLOCATED( pm10_av ) )  THEN
                   ALLOCATE( pm10_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                pm10_av = 0.0_wp

             CASE ( 's_BC', 's_DU', 's_NH', 's_NO', 's_OC', 's_SO4', 's_SS' )
                IF ( .NOT. ALLOCATED( s_mass_av ) )  THEN
                   ALLOCATE( s_mass_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg,ncomponents_mass+1) )
                ENDIF
                s_mass_av = 0.0_wp

             CASE ( 's_H2O' )
                IF ( .NOT. ALLOCATED( s_h2o_av ) )  THEN
                   ALLOCATE( s_h2o_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                s_h2o_av = 0.0_wp

             CASE DEFAULT
                CONTINUE

          END SELECT

       ENDIF

    ELSEIF ( mode == 'sum' )  THEN

       IF ( variable(7:11) ==  'N_bin' )  THEN
          IF ( ALLOCATED( nbins_av ) )  THEN
             READ( variable(12:),* ) char_to_int
             IF ( char_to_int >= 1  .AND. char_to_int <= SUM( nbin ) )  THEN
                ib = char_to_int
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         nbins_av(k,j,i,ib) = nbins_av(k,j,i,ib) + aerosol_number(ib)%conc(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ENDIF

       ELSEIF ( variable(7:11) ==  'm_bin' )  THEN
          IF ( ALLOCATED( mbins_av ) )  THEN
             READ( variable(12:),* ) char_to_int
             IF ( char_to_int >= 1  .AND. char_to_int <= SUM( nbin ) )  THEN
                ib = char_to_int
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         temp_bin = 0.0_wp
                         DO  ic = ib, nbins_aerosol * ncomponents_mass, nbins_aerosol
                            temp_bin = temp_bin + aerosol_mass(ic)%conc(k,j,i)
                         ENDDO
                         mbins_av(k,j,i,ib) = mbins_av(k,j,i,ib) + temp_bin
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ENDIF
       ELSE

          SELECT CASE ( TRIM( variable(7:) ) )

             CASE ( 'g_H2SO4', 'g_HNO3', 'g_NH3', 'g_OCNV', 'g_OCSV' )
                IF ( ALLOCATED( salsa_gases_av ) )  THEN

                   vari = TRIM( variable(9:) )  ! remove salsa_g_ from beginning

                   SELECT CASE( vari )

                      CASE( 'H2SO4' )
                         found_index = 1
                      CASE( 'HNO3' )
                         found_index = 2
                      CASE( 'NH3' )
                         found_index = 3
                      CASE( 'OCNV' )
                         found_index = 4
                      CASE( 'OCSV' )
                         found_index = 5

                   END SELECT

                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            salsa_gases_av(k,j,i,found_index) = salsa_gases_av(k,j,i,found_index)  &
                                                                + salsa_gas(found_index)%conc(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'LDSA' )
                IF ( ALLOCATED( ldsa_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            temp_bin = 0.0_wp
                            DO  ib = 1, nbins_aerosol
!
!--                            Diameter in micrometres
                               mean_d = 1.0E+6_wp * ra_dry(k,j,i,ib) * 2.0_wp
!
!--                            Deposition factor: alveolar (use ra_dry)
                               df = ( 0.01555_wp / mean_d ) * ( EXP( -0.416_wp * ( LOG( mean_d ) + &
                                      2.84_wp )**2 ) + 19.11_wp * EXP( -0.482_wp * ( LOG( mean_d ) &
                                      - 1.362_wp )**2 ) )
!
!--                            Lung-deposited surface area LDSA (units mum2/cm3)
                               temp_bin = temp_bin + pi * mean_d**2 * df * 1.0E-6_wp *             &
                                          aerosol_number(ib)%conc(k,j,i)
                            ENDDO
                            ldsa_av(k,j,i) = ldsa_av(k,j,i) + temp_bin
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'N_UFP' )
                IF ( ALLOCATED( nufp_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            temp_bin = 0.0_wp
                            DO  ib = 1, nbins_aerosol
                               IF ( 2.0_wp * ra_dry(k,j,i,ib) <= 0.1E-6_wp )  THEN
                                  temp_bin = temp_bin + aerosol_number(ib)%conc(k,j,i)
                               ENDIF
                            ENDDO
                            nufp_av(k,j,i) = nufp_av(k,j,i) + temp_bin
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'Ntot' )
               IF ( ALLOCATED( ntot_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            DO  ib = 1, nbins_aerosol
                               ntot_av(k,j,i) = ntot_av(k,j,i) + aerosol_number(ib)%conc(k,j,i)
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'PM0.1' )
                IF ( ALLOCATED( pm01_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            temp_bin = 0.0_wp
                            DO  ib = 1, nbins_aerosol
                               IF ( 2.0_wp * ra_dry(k,j,i,ib) <= 0.1E-6_wp )  THEN
                                  DO  ic = ib, nbins_aerosol * ncc, nbins_aerosol
                                     temp_bin = temp_bin + aerosol_mass(ic)%conc(k,j,i)
                                  ENDDO
                               ENDIF
                            ENDDO
                            pm01_av(k,j,i) = pm01_av(k,j,i) + temp_bin
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'PM2.5' )
                IF ( ALLOCATED( pm25_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            temp_bin = 0.0_wp
                            DO  ib = 1, nbins_aerosol
                               IF ( 2.0_wp * ra_dry(k,j,i,ib) <= 2.5E-6_wp )  THEN
                                  DO  ic = ib, nbins_aerosol * ncc, nbins_aerosol
                                     temp_bin = temp_bin + aerosol_mass(ic)%conc(k,j,i)
                                  ENDDO
                               ENDIF
                            ENDDO
                            pm25_av(k,j,i) = pm25_av(k,j,i) + temp_bin
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'PM10' )
                IF ( ALLOCATED( pm10_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            temp_bin = 0.0_wp
                            DO  ib = 1, nbins_aerosol
                               IF ( 2.0_wp * ra_dry(k,j,i,ib) <= 10.0E-6_wp )  THEN
                                  DO  ic = ib, nbins_aerosol * ncc, nbins_aerosol
                                     temp_bin = temp_bin + aerosol_mass(ic)%conc(k,j,i)
                                  ENDDO
                               ENDIF
                            ENDDO
                            pm10_av(k,j,i) = pm10_av(k,j,i) + temp_bin
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 's_BC', 's_DU', 's_NH', 's_NO', 's_OC', 's_SO4', 's_SS' )
                IF ( ALLOCATED( s_mass_av ) )  THEN
                   IF ( is_used( prtcl, TRIM( variable(9:) ) ) )  THEN  ! 9: remove salsa_s_
                      found_index = get_index( prtcl, TRIM( variable(9:) ) )
                      DO  i = nxlg, nxrg
                         DO  j = nysg, nyng
                            DO  k = nzb, nzt+1
                               DO  ic = ( found_index-1 ) * nbins_aerosol + 1, found_index * nbins_aerosol
                                  s_mass_av(k,j,i,found_index) = s_mass_av(k,j,i,found_index) +    &
                                                                 aerosol_mass(ic)%conc(k,j,i)
                               ENDDO
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDIF
                ENDIF

             CASE ( 's_H2O' )
                IF ( ALLOCATED( s_H2O_av ) )  THEN
                   found_index = get_index( prtcl,'H2O' )
                   to_be_resorted => s_h2o_av
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            DO  ic = ( found_index-1 ) * nbins_aerosol + 1, found_index * nbins_aerosol
                               s_h2o_av(k,j,i) = s_h2o_av(k,j,i) + aerosol_mass(ic)%conc(k,j,i)
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE DEFAULT
                CONTINUE

          END SELECT

       ENDIF

    ELSEIF ( mode == 'average' )  THEN

       IF ( variable(7:11) ==  'N_bin' )  THEN
          IF ( ALLOCATED( nbins_av ) )  THEN
             READ( variable(12:),* ) char_to_int
             IF ( char_to_int >= 1  .AND. char_to_int <= SUM( nbin ) )  THEN
                ib = char_to_int
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         nbins_av(k,j,i,ib) = nbins_av(k,j,i,ib) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ENDIF

       ELSEIF ( variable(7:11) ==  'm_bin' )  THEN
          IF ( ALLOCATED( mbins_av ) )  THEN
             READ( variable(12:),* ) char_to_int
             IF ( char_to_int >= 1  .AND. char_to_int <= SUM( nbin ) )  THEN
                ib = char_to_int
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         mbins_av(k,j,i,ib) = mbins_av(k,j,i,ib) / REAL( average_count_3d, KIND=wp)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ENDIF
       ELSE

          SELECT CASE ( TRIM( variable(7:) ) )

             CASE ( 'g_H2SO4', 'g_HNO3', 'g_NH3', 'g_OCNV', 'g_OCSV' )
                IF ( ALLOCATED( salsa_gases_av ) )  THEN
                   IF ( TRIM( variable(9:) ) == 'H2SO4' )  THEN  ! 9: remove salsa_g_ from beginning
                      found_index = 1
                   ELSEIF ( TRIM( variable(9:) ) == 'HNO3' )  THEN
                      found_index = 2
                   ELSEIF ( TRIM( variable(9:) ) == 'NH3' )  THEN
                      found_index = 3
                   ELSEIF ( TRIM( variable(9:) ) == 'OCNV' )  THEN
                      found_index = 4
                   ELSEIF ( TRIM( variable(9:) ) == 'OCSV' )  THEN
                      found_index = 5
                   ENDIF
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            salsa_gases_av(k,j,i,found_index) = salsa_gases_av(k,j,i,found_index)  &
                                                                / REAL( average_count_3d, KIND=wp )
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'LDSA' )
                IF ( ALLOCATED( ldsa_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            ldsa_av(k,j,i) = ldsa_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'N_UFP' )
                IF ( ALLOCATED( nufp_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            nufp_av(k,j,i) = nufp_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'Ntot' )
                IF ( ALLOCATED( ntot_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            ntot_av(k,j,i) = ntot_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'PM0.1' )
                IF ( ALLOCATED( pm01_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            pm01_av(k,j,i) = pm01_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'PM2.5' )
                IF ( ALLOCATED( pm25_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            pm25_av(k,j,i) = pm25_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'PM10' )
                IF ( ALLOCATED( pm10_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            pm10_av(k,j,i) = pm10_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 's_BC', 's_DU', 's_NH', 's_NO', 's_OC', 's_SO4', 's_SS' )
                IF ( ALLOCATED( s_mass_av ) )  THEN
                   IF ( is_used( prtcl, TRIM( variable(9:) ) ) )  THEN  ! 9: remove salsa_s_
                      found_index = get_index( prtcl, TRIM( variable(9:) ) )
                      DO  i = nxlg, nxrg
                         DO  j = nysg, nyng
                            DO  k = nzb, nzt+1
                               s_mass_av(k,j,i,found_index) = s_mass_av(k,j,i,found_index) /       &
                                                              REAL( average_count_3d, KIND=wp )
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDIF
                ENDIF

             CASE ( 's_H2O' )
                to_be_resorted => s_h2o_av
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         to_be_resorted(k,j,i) = to_be_resorted(k,j,i) /                           &
                                                 REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO

          END SELECT

       ENDIF
    ENDIF

 END SUBROUTINE salsa_3d_data_averaging


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining 2D output variables
!------------------------------------------------------------------------------!
 SUBROUTINE salsa_data_output_2d( av, variable, found, grid, mode, local_pf, two_d, nzb_do, nzt_do )

    USE indices

    USE kinds


    IMPLICIT NONE

    CHARACTER(LEN=*) ::  grid       !<
    CHARACTER(LEN=*) ::  mode       !<
    CHARACTER(LEN=*) ::  variable   !<
    CHARACTER(LEN=5) ::  vari       !<  trimmed format of variable

    INTEGER(iwp) ::  av           !<
    INTEGER(iwp) ::  char_to_int  !< for converting character to integer
    INTEGER(iwp) ::  found_index  !< index of a chemical compound
    INTEGER(iwp) ::  i            !<
    INTEGER(iwp) ::  ib           !< running index: size bins
    INTEGER(iwp) ::  ic           !< running index: mass bins
    INTEGER(iwp) ::  j            !<
    INTEGER(iwp) ::  k            !<
    INTEGER(iwp) ::  nzb_do       !<
    INTEGER(iwp) ::  nzt_do       !<

    LOGICAL ::  found  !<
    LOGICAL ::  two_d  !< flag parameter to indicate 2D variables (horizontal cross sections)

    REAL(wp) ::  df                       !< For calculating LDSA: fraction of particles
                                          !< depositing in the alveolar (or tracheobronchial)
                                          !< region of the lung. Depends on the particle size
    REAL(wp) ::  mean_d                   !< Particle diameter in micrometres
    REAL(wp) ::  temp_bin                 !< temporary array for calculating output variables

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf  !< output
!
!-- Next statement is to avoid compiler warning about unused variable. May be removed in future.
    IF ( two_d )  CONTINUE

    found = .TRUE.
    temp_bin  = 0.0_wp

    IF ( variable(7:11)  == 'N_bin' )  THEN

       READ( variable( 12:LEN( TRIM( variable ) ) - 3 ), * ) char_to_int
       IF ( char_to_int >= 1  .AND. char_to_int <= SUM( nbin ) )  THEN

          ib = char_to_int
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = MERGE( aerosol_number(ib)%conc(k,j,i), REAL( fill_value,   &
                                               KIND = wp ), BTEST( wall_flags_total_0(k,j,i), 0 ) )
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( nbins_av ) )  THEN
                ALLOCATE( nbins_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg,nbins_aerosol) )
                nbins_av = REAL( fill_value, KIND = wp )
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = MERGE( nbins_av(k,j,i,ib), REAL( fill_value, KIND = wp ),  &
                                               BTEST( wall_flags_total_0(k,j,i), 0 ) )
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid = 'zu'
       ENDIF

    ELSEIF ( variable(7:11)  == 'm_bin' )  THEN

       READ( variable( 12:LEN( TRIM( variable ) ) - 3 ), * ) char_to_int
       IF ( char_to_int >= 1  .AND. char_to_int <= SUM( nbin ) )  THEN

          ib = char_to_int
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      temp_bin = 0.0_wp
                      DO  ic = ib, ncomponents_mass * nbins_aerosol, nbins_aerosol
                         temp_bin = temp_bin + aerosol_mass(ic)%conc(k,j,i)
                      ENDDO
                      local_pf(i,j,k) = MERGE( temp_bin, REAL( fill_value, KIND = wp ),            &
                                               BTEST( wall_flags_total_0(k,j,i), 0 ) )
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( mbins_av ) )  THEN
                ALLOCATE( mbins_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg,nbins_aerosol) )
                mbins_av = REAL( fill_value, KIND = wp )
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = MERGE( mbins_av(k,j,i,ib), REAL( fill_value, KIND = wp ),  &
                                               BTEST( wall_flags_total_0(k,j,i), 0 ) )
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid = 'zu'
       ENDIF

    ELSE

       SELECT CASE ( TRIM( variable( 7:LEN( TRIM( variable ) ) - 3 ) ) )  ! cut out _xy, _xz or _yz

          CASE ( 'g_H2SO4', 'g_HNO3', 'g_NH3', 'g_OCNV', 'g_OCSV' )
             vari = TRIM( variable( 9:LEN( TRIM( variable ) ) - 3 ) )  ! 9: remove salsa_g_
             IF ( vari == 'H2SO4')  found_index = 1
             IF ( vari == 'HNO3')   found_index = 2
             IF ( vari == 'NH3')    found_index = 3
             IF ( vari == 'OCNV')   found_index = 4
             IF ( vari == 'OCSV')   found_index = 5
             IF ( av == 0 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         local_pf(i,j,k) = MERGE( salsa_gas(found_index)%conc(k,j,i),              &
                                                  REAL( fill_value,  KIND = wp ),                  &
                                                  BTEST( wall_flags_total_0(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
                IF ( .NOT. ALLOCATED( salsa_gases_av ) )  THEN
                   ALLOCATE( salsa_gases_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg,ngases_salsa) )
                   salsa_gases_av = REAL( fill_value, KIND = wp )
                ENDIF
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         local_pf(i,j,k) = MERGE( salsa_gases_av(k,j,i,found_index),               &
                                                  REAL( fill_value, KIND = wp ),                   &
                                                  BTEST( wall_flags_total_0(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

             IF ( mode == 'xy' )  grid = 'zu'

          CASE ( 'LDSA' )
             IF ( av == 0 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         temp_bin = 0.0_wp
                         DO  ib = 1, nbins_aerosol
!
!--                         Diameter in micrometres
                            mean_d = 1.0E+6_wp * ra_dry(k,j,i,ib) * 2.0_wp 
!
!--                         Deposition factor: alveolar
                            df = ( 0.01555_wp / mean_d ) * ( EXP( -0.416_wp * ( LOG( mean_d ) +    &
                                   2.84_wp )**2 ) + 19.11_wp * EXP( -0.482_wp * ( LOG( mean_d ) -  &
                                   1.362_wp )**2 ) )
!
!--                         Lung-deposited surface area LDSA (units mum2/cm3)
                            temp_bin = temp_bin + pi * mean_d**2 * df * 1.0E-6_wp *                &
                                       aerosol_number(ib)%conc(k,j,i)
                         ENDDO

                         local_pf(i,j,k) = MERGE( temp_bin, REAL( fill_value, KIND = wp ),         &
                                                  BTEST( wall_flags_total_0(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
                IF ( .NOT. ALLOCATED( ldsa_av ) )  THEN
                   ALLOCATE( ldsa_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ldsa_av = REAL( fill_value, KIND = wp )
                ENDIF
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         local_pf(i,j,k) = MERGE( ldsa_av(k,j,i), REAL( fill_value, KIND = wp ),   &
                                                  BTEST( wall_flags_total_0(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

             IF ( mode == 'xy' )  grid = 'zu'

          CASE ( 'N_UFP' )

             IF ( av == 0 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         temp_bin = 0.0_wp
                         DO  ib = 1, nbins_aerosol
                            IF ( 2.0_wp * ra_dry(k,j,i,ib) <= 0.1E-6_wp )  THEN
                               temp_bin = temp_bin + aerosol_number(ib)%conc(k,j,i)
                            ENDIF
                         ENDDO
                         local_pf(i,j,k) = MERGE( temp_bin, REAL( fill_value, KIND = wp ),         &
                                                  BTEST( wall_flags_total_0(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
                IF ( .NOT. ALLOCATED( nufp_av ) )  THEN
                   ALLOCATE( nufp_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   nufp_av = REAL( fill_value, KIND = wp )
                ENDIF
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         local_pf(i,j,k) = MERGE( nufp_av(k,j,i), REAL( fill_value, KIND = wp ),   &
                                                  BTEST( wall_flags_total_0(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

             IF ( mode == 'xy' )  grid = 'zu'

          CASE ( 'Ntot' )

             IF ( av == 0 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         temp_bin = 0.0_wp
                         DO  ib = 1, nbins_aerosol
                            temp_bin = temp_bin + aerosol_number(ib)%conc(k,j,i)
                         ENDDO
                         local_pf(i,j,k) = MERGE( temp_bin, REAL( fill_value, KIND = wp ),         &
                                                  BTEST( wall_flags_total_0(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
                IF ( .NOT. ALLOCATED( ntot_av ) )  THEN
                   ALLOCATE( ntot_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ntot_av = REAL( fill_value, KIND = wp )
                ENDIF
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         local_pf(i,j,k) = MERGE( ntot_av(k,j,i), REAL( fill_value, KIND = wp ),   &
                                                  BTEST( wall_flags_total_0(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

             IF ( mode == 'xy' )  grid = 'zu'

          CASE ( 'PM0.1' )
             IF ( av == 0 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         temp_bin = 0.0_wp
                         DO  ib = 1, nbins_aerosol
                            IF ( 2.0_wp * ra_dry(k,j,i,ib) <= 0.1E-6_wp )  THEN
                               DO  ic = ib, nbins_aerosol * ncc, nbins_aerosol
                                  temp_bin = temp_bin + aerosol_mass(ic)%conc(k,j,i)
                               ENDDO
                            ENDIF
                         ENDDO
                         local_pf(i,j,k) = MERGE( temp_bin, REAL( fill_value, KIND = wp ),         &
                                                  BTEST( wall_flags_total_0(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
                IF ( .NOT. ALLOCATED( pm01_av ) )  THEN
                   ALLOCATE( pm01_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   pm01_av = REAL( fill_value, KIND = wp )
                ENDIF
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         local_pf(i,j,k) = MERGE( pm01_av(k,j,i), REAL( fill_value, KIND = wp ),   &
                                                  BTEST( wall_flags_total_0(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

             IF ( mode == 'xy' )  grid = 'zu'

          CASE ( 'PM2.5' )
             IF ( av == 0 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         temp_bin = 0.0_wp
                         DO  ib = 1, nbins_aerosol
                            IF ( 2.0_wp * ra_dry(k,j,i,ib) <= 2.5E-6_wp )  THEN
                               DO  ic = ib, nbins_aerosol * ncc, nbins_aerosol
                                  temp_bin = temp_bin + aerosol_mass(ic)%conc(k,j,i)
                               ENDDO
                            ENDIF
                         ENDDO
                         local_pf(i,j,k) = MERGE( temp_bin, REAL( fill_value, KIND = wp ),         &
                                                  BTEST( wall_flags_total_0(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
                IF ( .NOT. ALLOCATED( pm25_av ) )  THEN
                   ALLOCATE( pm25_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   pm25_av = REAL( fill_value, KIND = wp )
                ENDIF
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         local_pf(i,j,k) = MERGE( pm25_av(k,j,i), REAL( fill_value, KIND = wp ),   &
                                                  BTEST( wall_flags_total_0(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

             IF ( mode == 'xy' )  grid = 'zu'

          CASE ( 'PM10' )
             IF ( av == 0 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         temp_bin = 0.0_wp
                         DO  ib = 1, nbins_aerosol
                            IF ( 2.0_wp * ra_dry(k,j,i,ib) <= 10.0E-6_wp )  THEN
                               DO  ic = ib, nbins_aerosol * ncc, nbins_aerosol
                                  temp_bin = temp_bin + aerosol_mass(ic)%conc(k,j,i)
                               ENDDO
                            ENDIF
                         ENDDO
                         local_pf(i,j,k) = MERGE( temp_bin,  REAL( fill_value, KIND = wp ),        &
                                                  BTEST( wall_flags_total_0(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
                IF ( .NOT. ALLOCATED( pm10_av ) )  THEN
                   ALLOCATE( pm10_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   pm10_av = REAL( fill_value, KIND = wp )
                ENDIF
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         local_pf(i,j,k) = MERGE( pm10_av(k,j,i), REAL( fill_value, KIND = wp ),   &
                                                  BTEST( wall_flags_total_0(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

             IF ( mode == 'xy' )  grid = 'zu'

          CASE ( 's_BC', 's_DU', 's_NH', 's_NO', 's_OC', 's_SO4', 's_SS' )
             vari = TRIM( variable( 9:LEN( TRIM( variable ) ) - 3 ) )  ! 9: remove salsa_s_
             IF ( is_used( prtcl, vari ) )  THEN
                found_index = get_index( prtcl, vari )
                IF ( av == 0 )  THEN
                   DO  i = nxl, nxr
                      DO  j = nys, nyn
                         DO  k = nzb_do, nzt_do
                            temp_bin = 0.0_wp
                            DO  ic = ( found_index-1 ) * nbins_aerosol+1, found_index * nbins_aerosol
                               temp_bin = temp_bin + aerosol_mass(ic)%conc(k,j,i)
                            ENDDO
                            local_pf(i,j,k) = MERGE( temp_bin, REAL( fill_value, KIND = wp ),      &
                                                     BTEST( wall_flags_total_0(k,j,i), 0 ) )
                         ENDDO
                      ENDDO
                   ENDDO
                ELSE
                   IF ( .NOT. ALLOCATED( s_mass_av ) )  THEN
                      ALLOCATE( s_mass_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg,ncomponents_mass) )
                      s_mass_av = REAL( fill_value, KIND = wp )
                   ENDIF
                   DO  i = nxl, nxr
                      DO  j = nys, nyn
                         DO  k = nzb_do, nzt_do
                            local_pf(i,j,k) = MERGE( s_mass_av(k,j,i,found_index),                 &
                                                     REAL( fill_value, KIND = wp ),                &
                                                     BTEST( wall_flags_total_0(k,j,i), 0 ) )
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF
             ELSE
                local_pf = fill_value
             ENDIF

             IF ( mode == 'xy' )  grid = 'zu'

          CASE ( 's_H2O' )
             found_index = get_index( prtcl, 'H2O' )
             IF ( av == 0 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         temp_bin = 0.0_wp
                         DO  ic = ( found_index-1 ) * nbins_aerosol+1, found_index * nbins_aerosol
                            temp_bin = temp_bin + aerosol_mass(ic)%conc(k,j,i)
                         ENDDO
                         local_pf(i,j,k) = MERGE( temp_bin, REAL( fill_value, KIND = wp ),         &
                                                  BTEST( wall_flags_total_0(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
     !           to_be_resorted => s_h2o_av
                IF ( .NOT. ALLOCATED( s_h2o_av ) )  THEN
                   ALLOCATE( s_h2o_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   s_h2o_av = REAL( fill_value, KIND = wp )
                ENDIF
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         local_pf(i,j,k) = MERGE( s_h2o_av(k,j,i), REAL( fill_value, KIND = wp ),  &
                                                  BTEST( wall_flags_total_0(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

             IF ( mode == 'xy' )  grid = 'zu'

          CASE DEFAULT
             found = .FALSE.
             grid  = 'none'

       END SELECT

    ENDIF

 END SUBROUTINE salsa_data_output_2d

!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining 3D output variables
!------------------------------------------------------------------------------!
 SUBROUTINE salsa_data_output_3d( av, variable, found, local_pf, nzb_do, nzt_do )

    USE indices

    USE kinds


    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(in) ::  variable   !<

    INTEGER(iwp) ::  av           !<
    INTEGER(iwp) ::  char_to_int  !< for converting character to integer
    INTEGER(iwp) ::  found_index  !< index of a chemical compound
    INTEGER(iwp) ::  ib           !< running index: size bins
    INTEGER(iwp) ::  ic           !< running index: mass bins
    INTEGER(iwp) ::  i            !<
    INTEGER(iwp) ::  j            !<
    INTEGER(iwp) ::  k            !<
    INTEGER(iwp) ::  nzb_do       !<
    INTEGER(iwp) ::  nzt_do       !<

    LOGICAL ::  found      !<

    REAL(wp) ::  df                       !< For calculating LDSA: fraction of particles
                                          !< depositing in the alveolar (or tracheobronchial)
                                          !< region of the lung. Depends on the particle size
    REAL(wp) ::  mean_d                   !< Particle diameter in micrometres
    REAL(wp) ::  temp_bin                 !< temporary array for calculating output variables

    REAL(sp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf  !< local

    found     = .TRUE.
    temp_bin  = 0.0_wp

    IF ( variable(7:11) == 'N_bin' )  THEN
       READ( variable(12:),* ) char_to_int
       IF ( char_to_int >= 1  .AND. char_to_int <= SUM( nbin ) )  THEN

          ib = char_to_int
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = MERGE( aerosol_number(ib)%conc(k,j,i), REAL( fill_value,   &
                                               KIND = wp ), BTEST( wall_flags_total_0(k,j,i), 0 ) )
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( nbins_av ) )  THEN
                ALLOCATE( nbins_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg,nbins_aerosol) )
                nbins_av = REAL( fill_value, KIND = wp )
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = MERGE( nbins_av(k,j,i,ib), REAL( fill_value, KIND = wp ),  &
                                               BTEST( wall_flags_total_0(k,j,i), 0 ) ) 
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDIF

    ELSEIF ( variable(7:11) == 'm_bin' )  THEN
       READ( variable(12:),* ) char_to_int
       IF ( char_to_int >= 1  .AND. char_to_int <= SUM( nbin ) )  THEN

          ib = char_to_int
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      temp_bin = 0.0_wp
                      DO  ic = ib, ncomponents_mass * nbins_aerosol, nbins_aerosol
                         temp_bin = temp_bin + aerosol_mass(ic)%conc(k,j,i)
                      ENDDO
                      local_pf(i,j,k) = MERGE( temp_bin, REAL( fill_value, KIND = wp ),            &
                                               BTEST( wall_flags_total_0(k,j,i), 0 ) )
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( mbins_av ) )  THEN
                ALLOCATE( mbins_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg,nbins_aerosol) )
                mbins_av = REAL( fill_value, KIND = wp )
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = MERGE( mbins_av(k,j,i,ib), REAL( fill_value, KIND = wp ),  &
                                               BTEST( wall_flags_total_0(k,j,i), 0 ) ) 
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDIF

    ELSE
       SELECT CASE ( TRIM( variable(7:) ) )

          CASE ( 'g_H2SO4', 'g_HNO3', 'g_NH3', 'g_OCNV',  'g_OCSV' )
             IF ( TRIM( variable(7:) ) == 'g_H2SO4')  found_index = 1
             IF ( TRIM( variable(7:) ) == 'g_HNO3')   found_index = 2
             IF ( TRIM( variable(7:) ) == 'g_NH3')    found_index = 3
             IF ( TRIM( variable(7:) ) == 'g_OCNV')   found_index = 4
             IF ( TRIM( variable(7:) ) == 'g_OCSV')   found_index = 5

             IF ( av == 0 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         local_pf(i,j,k) = MERGE( salsa_gas(found_index)%conc(k,j,i),              &
                                                  REAL( fill_value, KIND = wp ),                   &
                                                  BTEST( wall_flags_total_0(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
                IF ( .NOT. ALLOCATED( salsa_gases_av ) )  THEN
                   ALLOCATE( salsa_gases_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg,ngases_salsa) )
                   salsa_gases_av = REAL( fill_value, KIND = wp )
                ENDIF
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
!                          local_pf(i,j,k) = MERGE( to_be_resorted(k,j,i), REAL( fill_value,         &
!                                                KIND = wp ), BTEST( wall_flags_total_0(k,j,i), 0 ) )
                         local_pf(i,j,k) = MERGE( salsa_gases_av(k,j,i,found_index),               &
                                                  REAL( fill_value, KIND = wp ),                   &
                                                  BTEST( wall_flags_total_0(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'LDSA' )
             IF ( av == 0 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         temp_bin = 0.0_wp
                         DO  ib = 1, nbins_aerosol
!
!--                         Diameter in micrometres
                            mean_d = 1.0E+6_wp * ra_dry(k,j,i,ib) * 2.0_wp
!
!--                         Deposition factor: alveolar
                            df = ( 0.01555_wp / mean_d ) * ( EXP( -0.416_wp * ( LOG( mean_d ) +    &
                                   2.84_wp )**2 ) + 19.11_wp * EXP( -0.482_wp * ( LOG( mean_d ) -  &
                                   1.362_wp )**2 ) )
!
!--                         Lung-deposited surface area LDSA (units mum2/cm3)
                            temp_bin = temp_bin + pi * mean_d**2 * df * 1.0E-6_wp *                &
                                       aerosol_number(ib)%conc(k,j,i)
                         ENDDO
                         local_pf(i,j,k) = MERGE( temp_bin, REAL( fill_value, KIND = wp ),         &
                                                  BTEST( wall_flags_total_0(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
                IF ( .NOT. ALLOCATED( ldsa_av ) )  THEN
                   ALLOCATE( ldsa_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ldsa_av = REAL( fill_value, KIND = wp )
                ENDIF
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         local_pf(i,j,k) = MERGE( ldsa_av(k,j,i), REAL( fill_value, KIND = wp ),   &
                                                  BTEST( wall_flags_total_0(k,j,i), 0 ) ) 
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'N_UFP' )
             IF ( av == 0 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         temp_bin = 0.0_wp
                         DO  ib = 1, nbins_aerosol
                            IF ( 2.0_wp * ra_dry(k,j,i,ib) <= 0.1E-6_wp )  THEN
                               temp_bin = temp_bin + aerosol_number(ib)%conc(k,j,i)
                            ENDIF
                         ENDDO
                         local_pf(i,j,k) = MERGE( temp_bin, REAL( fill_value, KIND = wp ),         &
                                                  BTEST( wall_flags_total_0(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
                IF ( .NOT. ALLOCATED( nufp_av ) )  THEN
                   ALLOCATE( nufp_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   nufp_av = REAL( fill_value, KIND = wp )
                ENDIF
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         local_pf(i,j,k) = MERGE( nufp_av(k,j,i), REAL( fill_value, KIND = wp ),   &
                                                  BTEST( wall_flags_total_0(k,j,i), 0 ) ) 
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'Ntot' )
             IF ( av == 0 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         temp_bin = 0.0_wp
                         DO  ib = 1, nbins_aerosol
                            temp_bin = temp_bin + aerosol_number(ib)%conc(k,j,i)
                         ENDDO
                         local_pf(i,j,k) = MERGE( temp_bin, REAL( fill_value, KIND = wp ),         &
                                                  BTEST( wall_flags_total_0(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
                IF ( .NOT. ALLOCATED( ntot_av ) )  THEN
                   ALLOCATE( ntot_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ntot_av = REAL( fill_value, KIND = wp )
                ENDIF
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         local_pf(i,j,k) = MERGE( ntot_av(k,j,i), REAL( fill_value, KIND = wp ),   &
                                                  BTEST( wall_flags_total_0(k,j,i), 0 ) ) 
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'PM0.1' )
             IF ( av == 0 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         temp_bin = 0.0_wp
                         DO  ib = 1, nbins_aerosol
                            IF ( 2.0_wp * ra_dry(k,j,i,ib) <= 0.1E-6_wp )  THEN
                               DO  ic = ib, nbins_aerosol * ncc, nbins_aerosol
                                  temp_bin = temp_bin + aerosol_mass(ic)%conc(k,j,i)
                               ENDDO
                            ENDIF
                         ENDDO
                         local_pf(i,j,k) = MERGE( temp_bin, REAL( fill_value, KIND = wp ),         &
                                                  BTEST( wall_flags_total_0(k,j,i), 0 ) ) 
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
                IF ( .NOT. ALLOCATED( pm01_av ) )  THEN
                   ALLOCATE( pm01_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   pm01_av = REAL( fill_value, KIND = wp )
                ENDIF
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         local_pf(i,j,k) = MERGE( pm01_av(k,j,i), REAL( fill_value, KIND = wp ),   &
                                                  BTEST( wall_flags_total_0(k,j,i), 0 ) ) 
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'PM2.5' )
             IF ( av == 0 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         temp_bin = 0.0_wp
                         DO  ib = 1, nbins_aerosol
                            IF ( 2.0_wp * ra_dry(k,j,i,ib) <= 2.5E-6_wp )  THEN
                               DO  ic = ib, nbins_aerosol * ncc, nbins_aerosol
                                  temp_bin = temp_bin + aerosol_mass(ic)%conc(k,j,i)
                               ENDDO
                            ENDIF
                         ENDDO
                         local_pf(i,j,k) = MERGE( temp_bin, REAL( fill_value, KIND = wp ),         &
                                                  BTEST( wall_flags_total_0(k,j,i), 0 ) ) 
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
                IF ( .NOT. ALLOCATED( pm25_av ) )  THEN
                   ALLOCATE( pm25_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   pm25_av = REAL( fill_value, KIND = wp )
                ENDIF
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         local_pf(i,j,k) = MERGE( pm25_av(k,j,i), REAL( fill_value, KIND = wp ),   &
                                                  BTEST( wall_flags_total_0(k,j,i), 0 ) ) 
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'PM10' )
             IF ( av == 0 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         temp_bin = 0.0_wp
                         DO  ib = 1, nbins_aerosol
                            IF ( 2.0_wp * ra_dry(k,j,i,ib) <= 10.0E-6_wp )  THEN
                               DO  ic = ib, nbins_aerosol * ncc, nbins_aerosol
                                  temp_bin = temp_bin + aerosol_mass(ic)%conc(k,j,i)
                               ENDDO
                            ENDIF
                         ENDDO
                         local_pf(i,j,k) = MERGE( temp_bin, REAL( fill_value, KIND = wp ),         &
                                                  BTEST( wall_flags_total_0(k,j,i), 0 ) ) 
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
                IF ( .NOT. ALLOCATED( pm10_av ) )  THEN
                   ALLOCATE( pm10_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   pm10_av = REAL( fill_value, KIND = wp )
                ENDIF
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         local_pf(i,j,k) = MERGE( pm10_av(k,j,i), REAL( fill_value, KIND = wp ),   &
                                                  BTEST( wall_flags_total_0(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 's_BC', 's_DU', 's_NH', 's_NO', 's_OC', 's_SO4', 's_SS' )
             IF ( is_used( prtcl, TRIM( variable(9:) ) ) )  THEN  ! 9: remove salsa_s_
                found_index = get_index( prtcl, TRIM( variable(9:) ) )
                IF ( av == 0 )  THEN
                   DO  i = nxl, nxr
                      DO  j = nys, nyn
                         DO  k = nzb_do, nzt_do
                            temp_bin = 0.0_wp
                            DO  ic = ( found_index-1 ) * nbins_aerosol + 1, found_index * nbins_aerosol
                               temp_bin = temp_bin + aerosol_mass(ic)%conc(k,j,i)
                            ENDDO
                            local_pf(i,j,k) = MERGE( temp_bin, REAL( fill_value, KIND = wp ),      &
                                                     BTEST( wall_flags_total_0(k,j,i), 0 ) ) 
                         ENDDO
                      ENDDO
                   ENDDO
                ELSE
                   IF ( .NOT. ALLOCATED( s_mass_av ) )  THEN
                      ALLOCATE( s_mass_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg,ncomponents_mass) )
                      s_mass_av = REAL( fill_value, KIND = wp )
                   ENDIF
                   DO  i = nxl, nxr
                      DO  j = nys, nyn
                         DO  k = nzb_do, nzt_do
                            local_pf(i,j,k) = MERGE( s_mass_av(k,j,i,found_index),                 &
                                                     REAL( fill_value, KIND = wp ),                &
                                                     BTEST( wall_flags_total_0(k,j,i), 0 ) )
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF
             ENDIF

          CASE ( 's_H2O' )
             found_index = get_index( prtcl, 'H2O' )
             IF ( av == 0 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         temp_bin = 0.0_wp
                         DO  ic = ( found_index-1 ) * nbins_aerosol + 1, found_index * nbins_aerosol
                            temp_bin = temp_bin + aerosol_mass(ic)%conc(k,j,i)
                         ENDDO
                         local_pf(i,j,k) = MERGE( temp_bin, REAL( fill_value, KIND = wp ),         &
                                                  BTEST( wall_flags_total_0(k,j,i), 0 ) ) 
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
                IF ( .NOT. ALLOCATED( s_h2o_av ) )  THEN
                   ALLOCATE( s_h2o_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   s_h2o_av = REAL( fill_value, KIND = wp )
                ENDIF
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         local_pf(i,j,k) = MERGE( s_h2o_av(k,j,i), REAL( fill_value, KIND = wp ),  &
                                                  BTEST( wall_flags_total_0(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE DEFAULT
             found = .FALSE.

       END SELECT
    ENDIF

 END SUBROUTINE salsa_data_output_3d

!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining mask output variables
!------------------------------------------------------------------------------!
 SUBROUTINE salsa_data_output_mask( av, variable, found, local_pf, mid )

    USE arrays_3d,                                                                                 &
        ONLY:  tend

    USE control_parameters,                                                                        &
        ONLY:  mask_i, mask_j, mask_k, mask_size_l, mask_surface, nz_do3d

    IMPLICIT NONE

    CHARACTER(LEN=5) ::  grid      !< flag to distinquish between staggered grid
    CHARACTER(LEN=*) ::  variable  !<
    CHARACTER(LEN=7) ::  vari      !< trimmed format of variable

    INTEGER(iwp) ::  av             !<
    INTEGER(iwp) ::  char_to_int    !< for converting character to integer
    INTEGER(iwp) ::  found_index    !< index of a chemical compound
    INTEGER(iwp) ::  ib             !< loop index for aerosol size number bins
    INTEGER(iwp) ::  ic             !< loop index for chemical components
    INTEGER(iwp) ::  i              !< loop index in x-direction
    INTEGER(iwp) ::  j              !< loop index in y-direction
    INTEGER(iwp) ::  k              !< loop index in z-direction
    INTEGER(iwp) ::  im             !< loop index for masked variables
    INTEGER(iwp) ::  jm             !< loop index for masked variables
    INTEGER(iwp) ::  kk             !< loop index for masked output in z-direction
    INTEGER(iwp) ::  mid            !< masked output running index
    INTEGER(iwp) ::  ktt            !< k index of highest terrain surface

    LOGICAL ::  found      !<
    LOGICAL ::  resorted   !<

    REAL(wp) ::  df        !< For calculating LDSA: fraction of particles depositing in the alveolar
                           !< (or tracheobronchial) region of the lung. Depends on the particle size
    REAL(wp) ::  mean_d    !< Particle diameter in micrometres
    REAL(wp) ::  temp_bin  !< temporary array for calculating output variables

    REAL(wp), DIMENSION(mask_size_l(mid,1),mask_size_l(mid,2),mask_size_l(mid,3)) ::  local_pf   !<

    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg), TARGET ::  temp_array  !< temporary array

    REAL(wp), DIMENSION(:,:,:), POINTER ::  to_be_resorted  !< pointer

    found      = .TRUE.
    resorted   = .FALSE.
    grid       = 's'
    tend       = 0.0_wp
    temp_array = 0.0_wp
    temp_bin   = 0.0_wp

    IF ( variable(7:11) == 'N_bin' )  THEN
       READ( variable(12:),* ) char_to_int
       IF ( char_to_int >= 1  .AND. char_to_int <= SUM( nbin ) )  THEN
          ib = char_to_int
          IF ( av == 0 )  THEN
             IF ( .NOT. mask_surface(mid) )  THEN
                DO  i = 1, mask_size_l(mid,1)
                   DO  j = 1, mask_size_l(mid,2)
                      DO  k = 1, mask_size_l(mid,3)
                         local_pf(i,j,k) = aerosol_number(ib)%conc( mask_k(mid,k), mask_j(mid,j),  &
                                                                    mask_i(mid,i) )
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
                DO  i = 1, mask_size_l(mid,1)
                   DO  j = 1, mask_size_l(mid,2)
!
!--                   Get k index of the highest terraing surface
                      im = mask_i(mid,i)
                      jm = mask_j(mid,j)
                      ktt = MINLOC( MERGE( 1, 0, BTEST( wall_flags_total_0(:,jm,im), 5 )), &
                                                    DIM = 1 ) - 1
                      DO  k = 1, mask_size_l(mid,3)
                         kk = MIN( ktt+mask_k(mid,k), nzt+1 )
!
!--                      Set value if not in building
                         IF ( BTEST( wall_flags_total_0(kk,jm,im), 6 ) )  THEN
                            local_pf(i,j,k) = fill_value
                         ELSE
                            local_pf(i,j,k) = aerosol_number(ib)%conc(kk,jm,im)
                         ENDIF
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
             resorted = .TRUE.
          ELSE
             temp_array = nbins_av(:,:,:,ib)
             to_be_resorted => temp_array
          ENDIF
       ENDIF

    ELSEIF ( variable(7:11) == 'm_bin' )  THEN

       READ( variable(12:),* ) char_to_int
       IF ( char_to_int >= 1  .AND. char_to_int <= SUM( nbin ) )  THEN

          ib = char_to_int
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb, nz_do3d
                      temp_bin = 0.0_wp
                      DO  ic = ib, ncomponents_mass * nbins_aerosol, nbins_aerosol
                         temp_bin = temp_bin + aerosol_mass(ic)%conc(k,j,i)
                      ENDDO
                      tend(k,j,i) = temp_bin
                   ENDDO
                ENDDO
             ENDDO
             IF ( .NOT. mask_surface(mid) )  THEN
                DO  i = 1, mask_size_l(mid,1)
                   DO  j = 1, mask_size_l(mid,2)
                      DO  k = 1, mask_size_l(mid,3)
                         local_pf(i,j,k) = tend( mask_k(mid,k),  mask_j(mid,j), mask_i(mid,i) )
                      ENDDO
                   ENDDO
                ENDDO
             ELSE 
                DO  i = 1, mask_size_l(mid,1)
                   DO  j = 1, mask_size_l(mid,2)
!
!--                   Get k index of the highest terraing surface
                      im = mask_i(mid,i)
                      jm = mask_j(mid,j)
                      ktt = MINLOC( MERGE( 1, 0, BTEST( wall_flags_total_0(:,jm,im), 5 )), &
                                                    DIM = 1 ) - 1
                      DO  k = 1, mask_size_l(mid,3)
                         kk = MIN( ktt+mask_k(mid,k), nzt+1 )
!
!--                      Set value if not in building
                         IF ( BTEST( wall_flags_total_0(kk,jm,im), 6 ) )  THEN
                            local_pf(i,j,k) = fill_value
                         ELSE
                            local_pf(i,j,k) = tend(kk,jm,im)
                         ENDIF
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
             resorted = .TRUE.
          ELSE
             temp_array = mbins_av(:,:,:,ib)
             to_be_resorted => temp_array
          ENDIF
       ENDIF

    ELSE
       SELECT CASE ( TRIM( variable(7:) ) )

          CASE ( 'g_H2SO4', 'g_HNO3', 'g_NH3', 'g_OCNV',  'g_OCSV' )
             vari = TRIM( variable(7:) )
             IF ( av == 0 )  THEN
                IF ( vari == 'g_H2SO4')  to_be_resorted => salsa_gas(1)%conc
                IF ( vari == 'g_HNO3')   to_be_resorted => salsa_gas(2)%conc
                IF ( vari == 'g_NH3')    to_be_resorted => salsa_gas(3)%conc
                IF ( vari == 'g_OCNV')   to_be_resorted => salsa_gas(4)%conc
                IF ( vari == 'g_OCSV')   to_be_resorted => salsa_gas(5)%conc
             ELSE
                IF ( vari == 'g_H2SO4') temp_array = salsa_gases_av(:,:,:,1)
                IF ( vari == 'g_HNO3')  temp_array = salsa_gases_av(:,:,:,2)
                IF ( vari == 'g_NH3')   temp_array = salsa_gases_av(:,:,:,3)
                IF ( vari == 'g_OCNV')  temp_array = salsa_gases_av(:,:,:,4)
                IF ( vari == 'g_OCSV')  temp_array = salsa_gases_av(:,:,:,5)
                to_be_resorted => temp_array
             ENDIF

          CASE ( 'LDSA' )
             IF ( av == 0 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nz_do3d
                         temp_bin = 0.0_wp
                         DO  ib = 1, nbins_aerosol
!
!--                         Diameter in micrometres
                            mean_d = 1.0E+6_wp * ra_dry(k,j,i,ib) * 2.0_wp
!
!--                         Deposition factor: alveolar
                            df = ( 0.01555_wp / mean_d ) * ( EXP( -0.416_wp * ( LOG( mean_d ) +    &
                                   2.84_wp )**2 ) + 19.11_wp * EXP( -0.482_wp * ( LOG( mean_d ) -  &
                                   1.362_wp )**2 ) )
!
!--                         Lung-deposited surface area LDSA (units mum2/cm3)
                            temp_bin = temp_bin + pi * mean_d**2 * df * 1.0E-6_wp *                &
                                       aerosol_number(ib)%conc(k,j,i)
                         ENDDO
                         tend(k,j,i) = temp_bin
                      ENDDO
                   ENDDO
                ENDDO
                IF ( .NOT. mask_surface(mid) )  THEN
                   DO  i = 1, mask_size_l(mid,1)
                      DO  j = 1, mask_size_l(mid,2)
                         DO  k = 1, mask_size_l(mid,3)
                            local_pf(i,j,k) = tend( mask_k(mid,k),  mask_j(mid,j), mask_i(mid,i) )
                         ENDDO
                      ENDDO
                   ENDDO
                ELSE 
                   DO  i = 1, mask_size_l(mid,1)
                      DO  j = 1, mask_size_l(mid,2)
!
!--                      Get k index of the highest terraing surface
                         im = mask_i(mid,i)
                         jm = mask_j(mid,j)
                         ktt = MINLOC( MERGE( 1, 0, BTEST( wall_flags_total_0(:,jm,im), 5 )), &
                                                       DIM = 1 ) - 1
                         DO  k = 1, mask_size_l(mid,3)
                            kk = MIN( ktt+mask_k(mid,k), nzt+1 )
!
!--                         Set value if not in building
                            IF ( BTEST( wall_flags_total_0(kk,jm,im), 6 ) )  THEN
                               local_pf(i,j,k) = fill_value
                            ELSE
                               local_pf(i,j,k) = tend(kk,jm,im)
                            ENDIF
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF
                resorted = .TRUE.
             ELSE
                to_be_resorted => ldsa_av
             ENDIF

          CASE ( 'N_UFP' )
             IF ( av == 0 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nz_do3d
                         temp_bin = 0.0_wp
                         DO  ib = 1, nbins_aerosol
                            IF ( 2.0_wp * ra_dry(k,j,i,ib) <= 0.1E-6_wp )  THEN
                               temp_bin = temp_bin + aerosol_number(ib)%conc(k,j,i)
                            ENDIF
                         ENDDO
                         tend(k,j,i) = temp_bin
                      ENDDO
                   ENDDO
                ENDDO
                IF ( .NOT. mask_surface(mid) )  THEN
                   DO  i = 1, mask_size_l(mid,1)
                      DO  j = 1, mask_size_l(mid,2)
                         DO  k = 1, mask_size_l(mid,3)
                            local_pf(i,j,k) = tend( mask_k(mid,k),  mask_j(mid,j), mask_i(mid,i) )
                         ENDDO
                      ENDDO
                   ENDDO
                ELSE
                   DO  i = 1, mask_size_l(mid,1)
                      DO  j = 1, mask_size_l(mid,2)
!
!--                      Get k index of the highest terraing surface
                         im = mask_i(mid,i)
                         jm = mask_j(mid,j)
                         ktt = MINLOC( MERGE( 1, 0, BTEST( wall_flags_total_0(:,jm,im), 5 )), &
                                                       DIM = 1 ) - 1
                         DO  k = 1, mask_size_l(mid,3)
                            kk = MIN( ktt+mask_k(mid,k), nzt+1 )
!
!--                         Set value if not in building
                            IF ( BTEST( wall_flags_total_0(kk,jm,im), 6 ) )  THEN
                               local_pf(i,j,k) = fill_value
                            ELSE
                               local_pf(i,j,k) = tend(kk,jm,im)
                            ENDIF
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF
                resorted = .TRUE.
             ELSE
                to_be_resorted => nufp_av
             ENDIF

          CASE ( 'Ntot' )
             IF ( av == 0 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nz_do3d
                         temp_bin = 0.0_wp
                         DO  ib = 1, nbins_aerosol
                            temp_bin = temp_bin + aerosol_number(ib)%conc(k,j,i)
                         ENDDO
                         tend(k,j,i) = temp_bin
                      ENDDO
                   ENDDO
                ENDDO 
                IF ( .NOT. mask_surface(mid) )  THEN
                   DO  i = 1, mask_size_l(mid,1)
                      DO  j = 1, mask_size_l(mid,2)
                         DO  k = 1, mask_size_l(mid,3)
                            local_pf(i,j,k) = tend( mask_k(mid,k),  mask_j(mid,j), mask_i(mid,i) )
                         ENDDO
                      ENDDO
                   ENDDO
                ELSE 
                   DO  i = 1, mask_size_l(mid,1)
                      DO  j = 1, mask_size_l(mid,2)
!
!--                      Get k index of the highest terraing surface
                         im = mask_i(mid,i)
                         jm = mask_j(mid,j)
                         ktt = MINLOC( MERGE( 1, 0, BTEST( wall_flags_total_0(:,jm,im), 5 )), &
                                                       DIM = 1 ) - 1
                         DO  k = 1, mask_size_l(mid,3)
                            kk = MIN( ktt+mask_k(mid,k), nzt+1 )
!
!--                         Set value if not in building
                            IF ( BTEST( wall_flags_total_0(kk,jm,im), 6 ) )  THEN
                               local_pf(i,j,k) = fill_value
                            ELSE
                               local_pf(i,j,k) = tend(kk,jm,im)
                            ENDIF
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF
                resorted = .TRUE.
             ELSE
                to_be_resorted => ntot_av
             ENDIF

          CASE ( 'PM0.1' )
             IF ( av == 0 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nz_do3d
                         temp_bin = 0.0_wp
                         DO  ib = 1, nbins_aerosol
                            IF ( 2.0_wp * ra_dry(k,j,i,ib) <= 0.1E-6_wp )  THEN
                               DO  ic = ib, nbins_aerosol * ncc, nbins_aerosol
                                  temp_bin = temp_bin + aerosol_mass(ic)%conc(k,j,i)
                               ENDDO
                            ENDIF
                         ENDDO
                         tend(k,j,i) = temp_bin
                      ENDDO
                   ENDDO
                ENDDO 
                IF ( .NOT. mask_surface(mid) )  THEN
                   DO  i = 1, mask_size_l(mid,1)
                      DO  j = 1, mask_size_l(mid,2)
                         DO  k = 1, mask_size_l(mid,3)
                            local_pf(i,j,k) = tend( mask_k(mid,k),  mask_j(mid,j), mask_i(mid,i) )
                         ENDDO
                      ENDDO
                   ENDDO
                ELSE 
                   DO  i = 1, mask_size_l(mid,1)
                      DO  j = 1, mask_size_l(mid,2)
!
!--                      Get k index of the highest terraing surface
                         im = mask_i(mid,i)
                         jm = mask_j(mid,j)
                         ktt = MINLOC( MERGE( 1, 0, BTEST( wall_flags_total_0(:,jm,im), 5 )), &
                                                       DIM = 1 ) - 1
                         DO  k = 1, mask_size_l(mid,3)
                            kk = MIN( ktt+mask_k(mid,k), nzt+1 )
!
!--                         Set value if not in building
                            IF ( BTEST( wall_flags_total_0(kk,jm,im), 6 ) )  THEN
                               local_pf(i,j,k) = fill_value
                            ELSE
                               local_pf(i,j,k) = tend(kk,jm,im)
                            ENDIF
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF
                resorted = .TRUE.
             ELSE
                to_be_resorted => pm01_av
             ENDIF

          CASE ( 'PM2.5' )
             IF ( av == 0 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nz_do3d
                         temp_bin = 0.0_wp
                         DO  ib = 1, nbins_aerosol
                            IF ( 2.0_wp * ra_dry(k,j,i,ib) <= 2.5E-6_wp )  THEN
                               DO  ic = ib, nbins_aerosol * ncc, nbins_aerosol
                                  temp_bin = temp_bin + aerosol_mass(ic)%conc(k,j,i)
                               ENDDO
                            ENDIF
                         ENDDO
                         tend(k,j,i) = temp_bin
                      ENDDO
                   ENDDO
                ENDDO
                IF ( .NOT. mask_surface(mid) )  THEN
                   DO  i = 1, mask_size_l(mid,1)
                      DO  j = 1, mask_size_l(mid,2)
                         DO  k = 1, mask_size_l(mid,3)
                            local_pf(i,j,k) = tend( mask_k(mid,k),  mask_j(mid,j), mask_i(mid,i) )
                         ENDDO
                      ENDDO
                   ENDDO
                ELSE 
                   DO  i = 1, mask_size_l(mid,1)
                      DO  j = 1, mask_size_l(mid,2)
!
!--                      Get k index of the highest terraing surface
                         im = mask_i(mid,i)
                         jm = mask_j(mid,j)
                         ktt = MINLOC( MERGE( 1, 0, BTEST( wall_flags_total_0(:,jm,im), 5 )), &
                                                       DIM = 1 ) - 1
                         DO  k = 1, mask_size_l(mid,3)
                            kk = MIN( ktt+mask_k(mid,k), nzt+1 )
!
!--                         Set value if not in building
                            IF ( BTEST( wall_flags_total_0(kk,jm,im), 6 ) )  THEN
                               local_pf(i,j,k) = fill_value
                            ELSE
                               local_pf(i,j,k) = tend(kk,jm,im)
                            ENDIF
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF
                resorted = .TRUE.
             ELSE
                to_be_resorted => pm25_av
             ENDIF

          CASE ( 'PM10' )
             IF ( av == 0 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nz_do3d
                         temp_bin = 0.0_wp
                         DO  ib = 1, nbins_aerosol
                            IF ( 2.0_wp * ra_dry(k,j,i,ib) <= 10.0E-6_wp )  THEN
                               DO  ic = ib, nbins_aerosol * ncc, nbins_aerosol
                                  temp_bin = temp_bin + aerosol_mass(ic)%conc(k,j,i)
                               ENDDO
                            ENDIF
                         ENDDO
                         tend(k,j,i) = temp_bin
                      ENDDO
                   ENDDO
                ENDDO 
                IF ( .NOT. mask_surface(mid) )  THEN
                   DO  i = 1, mask_size_l(mid,1)
                      DO  j = 1, mask_size_l(mid,2)
                         DO  k = 1, mask_size_l(mid,3)
                            local_pf(i,j,k) = tend( mask_k(mid,k),  mask_j(mid,j), mask_i(mid,i) )
                         ENDDO
                      ENDDO
                   ENDDO
                ELSE 
                   DO  i = 1, mask_size_l(mid,1)
                      DO  j = 1, mask_size_l(mid,2)
!
!--                      Get k index of the highest terraing surface
                         im = mask_i(mid,i)
                         jm = mask_j(mid,j)
                         ktt = MINLOC( MERGE( 1, 0, BTEST( wall_flags_total_0(:,jm,im), 5 )), &
                                                       DIM = 1 ) - 1
                         DO  k = 1, mask_size_l(mid,3)
                            kk = MIN( ktt+mask_k(mid,k), nzt+1 )
!
!--                         Set value if not in building
                            IF ( BTEST( wall_flags_total_0(kk,jm,im), 6 ) )  THEN
                               local_pf(i,j,k) = fill_value
                            ELSE
                               local_pf(i,j,k) = tend(kk,jm,im)
                            ENDIF
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF
                resorted = .TRUE.
             ELSE
                to_be_resorted => pm10_av
             ENDIF

          CASE ( 's_BC', 's_DU', 's_NH', 's_NO', 's_OC', 's_SO4', 's_SS' )
             IF ( is_used( prtcl, TRIM( variable(9:) ) ) )  THEN
                found_index = get_index( prtcl, TRIM( variable(9:) ) )
                IF ( av == 0 )  THEN
                   DO  i = nxl, nxr
                      DO  j = nys, nyn
                         DO  k = nzb, nz_do3d
                            temp_bin = 0.0_wp
                            DO  ic = ( found_index-1 ) * nbins_aerosol + 1, found_index * nbins_aerosol
                               temp_bin = temp_bin + aerosol_mass(ic)%conc(k,j,i)
                            ENDDO
                            tend(k,j,i) = temp_bin
                         ENDDO
                      ENDDO
                   ENDDO
                   IF ( .NOT. mask_surface(mid) )  THEN
                      DO  i = 1, mask_size_l(mid,1)
                         DO  j = 1, mask_size_l(mid,2)
                            DO  k = 1, mask_size_l(mid,3)
                               local_pf(i,j,k) = tend( mask_k(mid,k), mask_j(mid,j), mask_i(mid,i) )
                            ENDDO
                         ENDDO
                      ENDDO
                   ELSE
                      DO  i = 1, mask_size_l(mid,1)
                         DO  j = 1, mask_size_l(mid,2)
   !
   !--                      Get k index of the highest terraing surface
                            im = mask_i(mid,i)
                            jm = mask_j(mid,j)
                            ktt = MINLOC( MERGE( 1, 0, BTEST( wall_flags_total_0(:,jm,im), 5 )), &
                                                          DIM = 1 ) - 1
                            DO  k = 1, mask_size_l(mid,3)
                               kk = MIN( ktt+mask_k(mid,k), nzt+1 )
   !
   !--                         Set value if not in building
                               IF ( BTEST( wall_flags_total_0(kk,jm,im), 6 ) )  THEN
                                  local_pf(i,j,k) = fill_value
                               ELSE
                                  local_pf(i,j,k) = tend(kk,jm,im)
                               ENDIF
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDIF
                   resorted = .TRUE.
                ELSE
                   temp_array = s_mass_av(:,:,:,found_index)
                   to_be_resorted => temp_array
                ENDIF
             ELSE
                local_pf = fill_value
             ENDIF

          CASE ( 's_H2O' )
             IF ( av == 0 )  THEN
                found_index = get_index( prtcl, 'H2O' )
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nz_do3d
                         temp_bin = 0.0_wp
                         DO  ic = ( found_index-1 ) * nbins_aerosol + 1, found_index * nbins_aerosol
                            temp_bin = temp_bin + aerosol_mass(ic)%conc(k,j,i)
                         ENDDO
                         tend(k,j,i) = temp_bin
                      ENDDO
                   ENDDO
                ENDDO
                IF ( .NOT. mask_surface(mid) )  THEN
                   DO  i = 1, mask_size_l(mid,1)
                      DO  j = 1, mask_size_l(mid,2)
                         DO  k = 1, mask_size_l(mid,3)
                            local_pf(i,j,k) = tend( mask_k(mid,k), mask_j(mid,j), mask_i(mid,i) )
                         ENDDO
                      ENDDO
                   ENDDO
                ELSE
                   DO  i = 1, mask_size_l(mid,1)
                      DO  j = 1, mask_size_l(mid,2)
!
!--                      Get k index of the highest terraing surface
                         im = mask_i(mid,i)
                         jm = mask_j(mid,j)
                         ktt = MINLOC( MERGE( 1, 0, BTEST( wall_flags_total_0(:,jm,im), 5 )), &
                                          DIM = 1 ) - 1
                         DO  k = 1, mask_size_l(mid,3)
                            kk = MIN( ktt+mask_k(mid,k), nzt+1 )
!
!--                         Set value if not in building
                            IF ( BTEST( wall_flags_total_0(kk,jm,im), 6 ) )  THEN
                               local_pf(i,j,k) = fill_value
                            ELSE
                               local_pf(i,j,k) =  tend(kk,jm,im)
                            ENDIF
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF
                resorted = .TRUE.
             ELSE
                to_be_resorted => s_h2o_av
             ENDIF

          CASE DEFAULT
             found = .FALSE.

       END SELECT
    ENDIF

    IF ( found  .AND.  .NOT. resorted )  THEN
       IF ( .NOT. mask_surface(mid) )  THEN
!
!--       Default masked output
          DO  i = 1, mask_size_l(mid,1)
             DO  j = 1, mask_size_l(mid,2)
                DO  k = 1, mask_size_l(mid,3)
                   local_pf(i,j,k) = to_be_resorted( mask_k(mid,k), mask_j(mid,j), mask_i(mid,i) )
                ENDDO
             ENDDO
          ENDDO
       ELSE
!
!--       Terrain-following masked output
          DO  i = 1, mask_size_l(mid,1)
             DO  j = 1, mask_size_l(mid,2)
!
!--             Get k index of the highest terraing surface
                im = mask_i(mid,i)
                jm = mask_j(mid,j)
                ktt = MINLOC( MERGE( 1, 0, BTEST( wall_flags_total_0(:,jm,im), 5 )), &
                                 DIM = 1 ) - 1
                DO  k = 1, mask_size_l(mid,3)
                   kk = MIN( ktt+mask_k(mid,k), nzt+1 )
!
!--                Set value if not in building
                   IF ( BTEST( wall_flags_total_0(kk,jm,im), 6 ) )  THEN
                      local_pf(i,j,k) = fill_value
                   ELSE
                      local_pf(i,j,k) = to_be_resorted(kk,jm,im)
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDIF
    ENDIF

 END SUBROUTINE salsa_data_output_mask

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Creates index tables for different (aerosol) components 
!------------------------------------------------------------------------------!
 SUBROUTINE component_index_constructor( self, ncomp, nlist, listcomp )

    IMPLICIT NONE

    INTEGER(iwp) ::  ii  !<
    INTEGER(iwp) ::  jj  !<

    INTEGER(iwp), INTENT(in) ::  nlist ! < Maximum number of components

    INTEGER(iwp), INTENT(inout) ::  ncomp  !< Number of components

    CHARACTER(LEN=3), INTENT(in) ::  listcomp(nlist)  !< List cof component names

    TYPE(component_index), INTENT(inout) ::  self  !< Object containing the indices of different
                                                   !< aerosol components

    ncomp = 0

    DO WHILE ( listcomp(ncomp+1) /= '  ' .AND. ncomp < nlist )
       ncomp = ncomp + 1
    ENDDO

    self%ncomp = ncomp
    ALLOCATE( self%ind(ncomp), self%comp(ncomp) )

    DO  ii = 1, ncomp
       self%ind(ii) = ii
    ENDDO

    jj = 1
    DO  ii = 1, nlist
       IF ( listcomp(ii) == '') CYCLE
       self%comp(jj) = listcomp(ii)
       jj = jj + 1
    ENDDO

 END SUBROUTINE component_index_constructor

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Gives the index of a component in the component list
!------------------------------------------------------------------------------!
 INTEGER FUNCTION get_index( self, incomp )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(in) ::  incomp !< Component name

    INTEGER(iwp) ::  ii  !< index

    TYPE(component_index), INTENT(in) ::  self  !< Object containing the indices of different
                                                !< aerosol components
    IF ( ANY( self%comp == incomp ) )  THEN
       ii = 1
       DO WHILE ( (self%comp(ii) /= incomp) )
          ii = ii + 1
       ENDDO
       get_index = ii
    ELSEIF ( incomp == 'H2O' )  THEN
       get_index = self%ncomp + 1
    ELSE
       WRITE( message_string, * ) 'Incorrect component name given!'
       CALL message( 'get_index', 'PA0591', 1, 2, 0, 6, 0 )
    ENDIF

 END FUNCTION get_index

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Tells if the (aerosol) component is being used in the simulation
!------------------------------------------------------------------------------!
 LOGICAL FUNCTION is_used( self, icomp )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(in) ::  icomp !< Component name

    TYPE(component_index), INTENT(in) ::  self  !< Object containing the indices of different
                                                !< aerosol components

    IF ( ANY(self%comp == icomp) ) THEN
       is_used = .TRUE.
    ELSE
       is_used = .FALSE.
    ENDIF

 END FUNCTION

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Set the lateral and top boundary conditions in case the PALM domain is
!> nested offline in a mesoscale model. Further, average boundary data and
!> determine mean profiles, further used for correct damping in the sponge
!> layer.
!------------------------------------------------------------------------------!
 SUBROUTINE salsa_nesting_offl_bc

    USE control_parameters,                                                                        &
        ONLY:  bc_dirichlet_l, bc_dirichlet_n, bc_dirichlet_r, bc_dirichlet_s, dt_3d,              &
               time_since_reference_point

    USE exchange_horiz_mod,                                                    &
        ONLY:  exchange_horiz

    USE indices,                                                                                   &
        ONLY:  nbgp, nxl, nxr, nyn, nys, nzb, nzt

    IMPLICIT NONE

    INTEGER(iwp) ::  i    !< running index x-direction
    INTEGER(iwp) ::  ib   !< running index for aerosol number bins
    INTEGER(iwp) ::  ic   !< running index for aerosol mass bins
    INTEGER(iwp) ::  icc  !< running index for aerosol mass bins
    INTEGER(iwp) ::  ig   !< running index for gaseous species
    INTEGER(iwp) ::  j    !< running index y-direction
    INTEGER(iwp) ::  k    !< running index z-direction

    REAL(wp) ::  fac_dt  !< interpolation factor

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ref_mconc    !< reference profile for aerosol mass
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ref_mconc_l  !< reference profile for aerosol mass: subdomain
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ref_nconc    !< reference profile for aerosol number
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ref_nconc_l  !< reference profile for aerosol_number: subdomain
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ref_gconc    !< reference profile for gases
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ref_gconc_l  !< reference profile for gases: subdomain

!
!-- Skip input if no forcing from larger-scale models is applied.
    IF ( .NOT. nesting_offline_salsa )  RETURN
!
!-- Allocate temporary arrays to compute salsa mean profiles
    ALLOCATE( ref_gconc(nzb:nzt+1,1:ngases_salsa), ref_gconc_l(nzb:nzt+1,1:ngases_salsa),          &
              ref_mconc(nzb:nzt+1,1:nbins_aerosol*ncomponents_mass),                               &
              ref_mconc_l(nzb:nzt+1,1:nbins_aerosol*ncomponents_mass),                             &
              ref_nconc(nzb:nzt+1,1:nbins_aerosol), ref_nconc_l(nzb:nzt+1,1:nbins_aerosol) )
    ref_gconc   = 0.0_wp
    ref_gconc_l = 0.0_wp
    ref_mconc   = 0.0_wp
    ref_mconc_l = 0.0_wp
    ref_nconc   = 0.0_wp
    ref_nconc_l = 0.0_wp

!
!-- Determine interpolation factor and limit it to 1. This is because t+dt can slightly exceed
!-- time(tind_p) before boundary data is updated again.
    fac_dt = ( time_since_reference_point - salsa_nest_offl%time(salsa_nest_offl%tind) + dt_3d ) / &
             ( salsa_nest_offl%time(salsa_nest_offl%tind_p) -                                      &
               salsa_nest_offl%time(salsa_nest_offl%tind) )
    fac_dt = MIN( 1.0_wp, fac_dt )

    IF ( bc_dirichlet_l )  THEN
       DO  ib = 1, nbins_aerosol
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                aerosol_number(ib)%conc(k,j,-1) = ( 1.0_wp - fac_dt ) *                            &
                                                  salsa_nest_offl%nconc_left(0,k,j,ib) + fac_dt *  &
                                                  salsa_nest_offl%nconc_left(1,k,j,ib)
             ENDDO
             ref_nconc_l(nzb+1:nzt,ib) = ref_nconc_l(nzb+1:nzt,ib) +                               &
                                         aerosol_number(ib)%conc(nzb+1:nzt,j,-1)
          ENDDO
          DO  ic = 1, ncomponents_mass
             icc = ( ic-1 ) * nbins_aerosol + ib
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   aerosol_mass(icc)%conc(k,j,-1) = ( 1.0_wp - fac_dt ) *                          &
                                                    salsa_nest_offl%mconc_left(0,k,j,icc) + fac_dt &
                                                    * salsa_nest_offl%mconc_left(1,k,j,icc)
                ENDDO
                ref_mconc_l(nzb+1:nzt,icc) = ref_mconc_l(nzb+1:nzt,icc) +                          &
                                             aerosol_mass(icc)%conc(nzb+1:nzt,j,-1)
             ENDDO
          ENDDO
       ENDDO
       IF ( .NOT. salsa_gases_from_chem )  THEN
          DO  ig = 1, ngases_salsa
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   salsa_gas(ig)%conc(k,j,-1) = ( 1.0_wp - fac_dt ) *                              &
                                                salsa_nest_offl%gconc_left(0,k,j,ig) + fac_dt *    &
                                                salsa_nest_offl%gconc_left(1,k,j,ig)
                ENDDO
                ref_gconc_l(nzb+1:nzt,ig) = ref_gconc_l(nzb+1:nzt,ig) +                            &
                                            salsa_gas(ig)%conc(nzb+1:nzt,j,-1)
             ENDDO
          ENDDO
       ENDIF
    ENDIF

    IF ( bc_dirichlet_r )  THEN
       DO  ib = 1, nbins_aerosol
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                aerosol_number(ib)%conc(k,j,nxr+1) = ( 1.0_wp - fac_dt ) *                         &
                                                  salsa_nest_offl%nconc_right(0,k,j,ib) + fac_dt * &
                                                  salsa_nest_offl%nconc_right(1,k,j,ib)
             ENDDO
             ref_nconc_l(nzb+1:nzt,ib) = ref_nconc_l(nzb+1:nzt,ib) +                               &
                                         aerosol_number(ib)%conc(nzb+1:nzt,j,nxr+1)
          ENDDO
          DO  ic = 1, ncomponents_mass
             icc = ( ic-1 ) * nbins_aerosol + ib
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   aerosol_mass(icc)%conc(k,j,nxr+1) = ( 1.0_wp - fac_dt ) *                       &
                                                    salsa_nest_offl%mconc_right(0,k,j,icc) + fac_dt&
                                                    * salsa_nest_offl%mconc_right(1,k,j,icc)
                ENDDO
                ref_mconc_l(nzb+1:nzt,icc) = ref_mconc_l(nzb+1:nzt,icc) +                          &
                                             aerosol_mass(icc)%conc(nzb+1:nzt,j,nxr+1)
             ENDDO
          ENDDO
       ENDDO
       IF ( .NOT. salsa_gases_from_chem )  THEN
          DO  ig = 1, ngases_salsa
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   salsa_gas(ig)%conc(k,j,nxr+1) = ( 1.0_wp - fac_dt ) *                           &
                                                   salsa_nest_offl%gconc_right(0,k,j,ig) + fac_dt *&
                                                   salsa_nest_offl%gconc_right(1,k,j,ig)
                ENDDO
                ref_gconc_l(nzb+1:nzt,ig) = ref_gconc_l(nzb+1:nzt,ig) +                            &
                                            salsa_gas(ig)%conc(nzb+1:nzt,j,nxr+1)
             ENDDO
          ENDDO
       ENDIF
    ENDIF

    IF ( bc_dirichlet_n )  THEN
       DO  ib = 1, nbins_aerosol
          DO  i = nxl, nxr
             DO  k = nzb+1, nzt
                aerosol_number(ib)%conc(k,nyn+1,i) = ( 1.0_wp - fac_dt ) *                         &
                                                  salsa_nest_offl%nconc_north(0,k,i,ib) + fac_dt * &
                                                  salsa_nest_offl%nconc_north(1,k,i,ib)
             ENDDO
             ref_nconc_l(nzb+1:nzt,ib) = ref_nconc_l(nzb+1:nzt,ib) +                               &
                                         aerosol_number(ib)%conc(nzb+1:nzt,nyn+1,i)
          ENDDO
          DO  ic = 1, ncomponents_mass
             icc = ( ic-1 ) * nbins_aerosol + ib
             DO  i = nxl, nxr
                DO  k = nzb+1, nzt
                   aerosol_mass(icc)%conc(k,nyn+1,i) = ( 1.0_wp - fac_dt ) *                       &
                                                    salsa_nest_offl%mconc_north(0,k,i,icc) + fac_dt&
                                                    * salsa_nest_offl%mconc_north(1,k,i,icc)
                ENDDO
                ref_mconc_l(nzb+1:nzt,icc) = ref_mconc_l(nzb+1:nzt,icc) +                          &
                                             aerosol_mass(icc)%conc(nzb+1:nzt,nyn+1,i)
             ENDDO
          ENDDO
       ENDDO
       IF ( .NOT. salsa_gases_from_chem )  THEN
          DO  ig = 1, ngases_salsa
             DO  i = nxl, nxr
                DO  k = nzb+1, nzt
                   salsa_gas(ig)%conc(k,nyn+1,i) = ( 1.0_wp - fac_dt ) *                           &
                                                   salsa_nest_offl%gconc_north(0,k,i,ig) + fac_dt *&
                                                   salsa_nest_offl%gconc_north(1,k,i,ig)
                ENDDO
                ref_gconc_l(nzb+1:nzt,ig) = ref_gconc_l(nzb+1:nzt,ig) +                            &
                                            salsa_gas(ig)%conc(nzb+1:nzt,nyn+1,i)
             ENDDO
          ENDDO
       ENDIF
    ENDIF

    IF ( bc_dirichlet_s )  THEN
       DO  ib = 1, nbins_aerosol
          DO  i = nxl, nxr
             DO  k = nzb+1, nzt
                aerosol_number(ib)%conc(k,-1,i) = ( 1.0_wp - fac_dt ) *                            &
                                                  salsa_nest_offl%nconc_south(0,k,i,ib) + fac_dt * &
                                                  salsa_nest_offl%nconc_south(1,k,i,ib)
             ENDDO
             ref_nconc_l(nzb+1:nzt,ib) = ref_nconc_l(nzb+1:nzt,ib) +                               &
                                         aerosol_number(ib)%conc(nzb+1:nzt,-1,i)
          ENDDO
          DO  ic = 1, ncomponents_mass
             icc = ( ic-1 ) * nbins_aerosol + ib
             DO  i = nxl, nxr
                DO  k = nzb+1, nzt
                   aerosol_mass(icc)%conc(k,-1,i) = ( 1.0_wp - fac_dt ) *                          &
                                                    salsa_nest_offl%mconc_south(0,k,i,icc) + fac_dt&
                                                    * salsa_nest_offl%mconc_south(1,k,i,icc)
                ENDDO
                ref_mconc_l(nzb+1:nzt,icc) = ref_mconc_l(nzb+1:nzt,icc) +                          &
                                             aerosol_mass(icc)%conc(nzb+1:nzt,-1,i)
             ENDDO
          ENDDO
       ENDDO
       IF ( .NOT. salsa_gases_from_chem )  THEN
          DO  ig = 1, ngases_salsa
             DO  i = nxl, nxr
                DO  k = nzb+1, nzt
                   salsa_gas(ig)%conc(k,-1,i) = ( 1.0_wp - fac_dt ) *                              &
                                                salsa_nest_offl%gconc_south(0,k,i,ig) + fac_dt *   &
                                                salsa_nest_offl%gconc_south(1,k,i,ig)
                ENDDO
                ref_gconc_l(nzb+1:nzt,ig) = ref_gconc_l(nzb+1:nzt,ig) +                            &
                                            salsa_gas(ig)%conc(nzb+1:nzt,-1,i)
             ENDDO
          ENDDO
       ENDIF
    ENDIF
!
!-- Top boundary
    DO  ib = 1, nbins_aerosol
       DO  i = nxl, nxr
          DO  j = nys, nyn
             aerosol_number(ib)%conc(nzt+1,j,i) = ( 1.0_wp - fac_dt ) *                            &
                                                  salsa_nest_offl%nconc_top(0,j,i,ib) + fac_dt *   &
                                                  salsa_nest_offl%nconc_top(1,j,i,ib)
             ref_nconc_l(nzt+1,ib) = ref_nconc_l(nzt+1,ib) + aerosol_number(ib)%conc(nzt+1,j,i)
          ENDDO
       ENDDO
       DO  ic = 1, ncomponents_mass
          icc = ( ic-1 ) * nbins_aerosol + ib
          DO  i = nxl, nxr
             DO  j = nys, nyn
                aerosol_mass(icc)%conc(nzt+1,j,i) = ( 1.0_wp - fac_dt ) *                          &
                                                    salsa_nest_offl%mconc_top(0,j,i,icc) + fac_dt *&
                                                    salsa_nest_offl%mconc_top(1,j,i,icc)
                ref_mconc_l(nzt+1,icc) = ref_mconc_l(nzt+1,icc) + aerosol_mass(icc)%conc(nzt+1,j,i)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    IF ( .NOT. salsa_gases_from_chem )  THEN
       DO  ig = 1, ngases_salsa
          DO  i = nxl, nxr
             DO  j = nys, nyn
                salsa_gas(ig)%conc(nzt+1,j,i) = ( 1.0_wp - fac_dt ) *                              &
                                                salsa_nest_offl%gconc_top(0,j,i,ig) + fac_dt *     &
                                                salsa_nest_offl%gconc_top(1,j,i,ig)
                ref_gconc_l(nzt+1,ig) = ref_gconc_l(nzt+1,ig) + salsa_gas(ig)%conc(nzt+1,j,i)
             ENDDO
          ENDDO
       ENDDO
    ENDIF
!
!-- Do local exchange
    DO  ib = 1, nbins_aerosol
       CALL exchange_horiz( aerosol_number(ib)%conc, nbgp )
       DO  ic = 1, ncomponents_mass
          icc = ( ic-1 ) * nbins_aerosol + ib
          CALL exchange_horiz( aerosol_mass(icc)%conc, nbgp )
       ENDDO
    ENDDO
    IF ( .NOT. salsa_gases_from_chem )  THEN
       DO  ig = 1, ngases_salsa
          CALL exchange_horiz( salsa_gas(ig)%conc, nbgp )
       ENDDO
    ENDIF
!
!-- In case of Rayleigh damping, where the initial profiles are still used, update these profiles
!-- from the averaged boundary data. But first, average these data.
#if defined( __parallel )
    IF ( .NOT. salsa_gases_from_chem )                                                             &
       CALL MPI_ALLREDUCE( ref_gconc_l, ref_gconc, ( nzt+1-nzb+1 ) * SIZE( ref_gconc(nzb,:) ),     &
                           MPI_REAL, MPI_SUM, comm2d, ierr )
    CALL MPI_ALLREDUCE( ref_mconc_l, ref_mconc, ( nzt+1-nzb+1 ) * SIZE( ref_mconc(nzb,:) ),        &
                        MPI_REAL, MPI_SUM, comm2d, ierr )
    CALL MPI_ALLREDUCE( ref_nconc_l, ref_nconc, ( nzt+1-nzb+1 ) * SIZE( ref_nconc(nzb,:) ),        &
                        MPI_REAL, MPI_SUM, comm2d, ierr )
#else
    IF ( .NOT. salsa_gases_from_chem )  ref_gconc = ref_gconc_l
    ref_mconc = ref_mconc_l
    ref_nconc = ref_nconc_l
#endif
!
!-- Average data. Note, reference profiles up to nzt are derived from lateral boundaries, at the
!-- model top it is derived from the top boundary. Thus, number of input data is different from
!-- nzb:nzt compared to nzt+1.
!-- Derived from lateral boundaries.
    IF ( .NOT. salsa_gases_from_chem )                                                             &
       ref_gconc(nzb:nzt,:) = ref_gconc(nzb:nzt,:) / REAL( 2.0_wp * ( ny + 1 + nx + 1 ), KIND = wp )
    ref_mconc(nzb:nzt,:) = ref_mconc(nzb:nzt,:) / REAL( 2.0_wp * ( ny + 1 + nx + 1 ), KIND = wp )
    ref_nconc(nzb:nzt,:) = ref_nconc(nzb:nzt,:) / REAL( 2.0_wp * ( ny + 1 + nx + 1 ), KIND = wp )
!
!-- Derived from top boundary
    IF ( .NOT. salsa_gases_from_chem )                                                             &
       ref_gconc(nzt+1,:) = ref_gconc(nzt+1,:) / REAL( ( ny + 1 ) * ( nx + 1 ), KIND = wp )
    ref_mconc(nzt+1,:) = ref_mconc(nzt+1,:) / REAL( ( ny + 1 ) * ( nx + 1 ), KIND = wp )
    ref_nconc(nzt+1,:) = ref_nconc(nzt+1,:) / REAL( ( ny + 1 ) * ( nx + 1 ), KIND = wp )
!
!-- Write onto init profiles, which are used for damping. Also set lower boundary condition.
    DO  ib = 1, nbins_aerosol
       aerosol_number(ib)%init(:)   = ref_nconc(:,ib)
       aerosol_number(ib)%init(nzb) = aerosol_number(ib)%init(nzb+1)
       DO  ic = 1, ncomponents_mass
          icc = ( ic-1 ) * nbins_aerosol + ib
          aerosol_mass(icc)%init(:)   = ref_mconc(:,icc)
          aerosol_mass(icc)%init(nzb) = aerosol_mass(icc)%init(nzb+1)
       ENDDO
    ENDDO
    IF ( .NOT. salsa_gases_from_chem )  THEN
       DO  ig = 1, ngases_salsa
          salsa_gas(ig)%init(:)   = ref_gconc(:,ig)
          salsa_gas(ig)%init(nzb) = salsa_gas(ig)%init(nzb+1)
       ENDDO
    ENDIF

    DEALLOCATE( ref_gconc, ref_gconc_l, ref_mconc, ref_mconc_l, ref_nconc, ref_nconc_l )

 END SUBROUTINE salsa_nesting_offl_bc

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate arrays used to read boundary data from NetCDF file and initialize
!> boundary data.
!------------------------------------------------------------------------------!
 SUBROUTINE salsa_nesting_offl_init

    USE control_parameters,                                                                        &
        ONLY:  end_time, initializing_actions, spinup_time

    USE palm_date_time_mod,                                                                        &
        ONLY:  get_date_time

    IMPLICIT NONE

    INTEGER(iwp) ::  ib          !< running index for aerosol number bins
    INTEGER(iwp) ::  ic          !< running index for aerosol mass bins
    INTEGER(iwp) ::  icc         !< additional running index for aerosol mass bins
    INTEGER(iwp) ::  ig          !< running index for gaseous species
    INTEGER(iwp) ::  nmass_bins  !< number of aerosol mass bins

    nmass_bins = nbins_aerosol * ncomponents_mass
!
!-- Allocate arrays for reading boundary values. Arrays will incorporate 2 time levels in order to
!-- interpolate in between.
    IF ( nesting_offline_salsa )  THEN
       IF ( bc_dirichlet_l )  THEN
          ALLOCATE( salsa_nest_offl%nconc_left(0:1,nzb+1:nzt,nys:nyn,1:nbins_aerosol) )
          ALLOCATE( salsa_nest_offl%mconc_left(0:1,nzb+1:nzt,nys:nyn,1:nmass_bins) )
       ENDIF
       IF ( bc_dirichlet_r )  THEN
          ALLOCATE( salsa_nest_offl%nconc_right(0:1,nzb+1:nzt,nys:nyn,1:nbins_aerosol) )
          ALLOCATE( salsa_nest_offl%mconc_right(0:1,nzb+1:nzt,nys:nyn,1:nmass_bins) )
       ENDIF
       IF ( bc_dirichlet_n )  THEN
          ALLOCATE( salsa_nest_offl%nconc_north(0:1,nzb+1:nzt,nxl:nxr,1:nbins_aerosol) )
          ALLOCATE( salsa_nest_offl%mconc_north(0:1,nzb+1:nzt,nxl:nxr,1:nmass_bins) )
       ENDIF
       IF ( bc_dirichlet_s )  THEN
          ALLOCATE( salsa_nest_offl%nconc_south(0:1,nzb+1:nzt,nxl:nxr,1:nbins_aerosol) )
          ALLOCATE( salsa_nest_offl%mconc_south(0:1,nzb+1:nzt,nxl:nxr,1:nmass_bins) )
       ENDIF
       ALLOCATE( salsa_nest_offl%nconc_top(0:1,nys:nyn,nxl:nxr,1:nbins_aerosol) )
       ALLOCATE( salsa_nest_offl%mconc_top(0:1,nys:nyn,nxl:nxr,1:nmass_bins) )

       IF ( .NOT. salsa_gases_from_chem )  THEN
          IF ( bc_dirichlet_l )  THEN
             ALLOCATE( salsa_nest_offl%gconc_left(0:1,nzb+1:nzt,nys:nyn,1:ngases_salsa) )
          ENDIF
          IF ( bc_dirichlet_r )  THEN
             ALLOCATE( salsa_nest_offl%gconc_right(0:1,nzb+1:nzt,nys:nyn,1:ngases_salsa) )
          ENDIF
          IF ( bc_dirichlet_n )  THEN
             ALLOCATE( salsa_nest_offl%gconc_north(0:1,nzb+1:nzt,nxl:nxr,1:ngases_salsa) )
          ENDIF
          IF ( bc_dirichlet_s )  THEN
             ALLOCATE( salsa_nest_offl%gconc_south(0:1,nzb+1:nzt,nxl:nxr,1:ngases_salsa) )
          ENDIF
          ALLOCATE( salsa_nest_offl%gconc_top(0:1,nys:nyn,nxl:nxr,1:ngases_salsa) )
       ENDIF

!
!--    Read data at lateral and top boundaries from a larger-scale model
       CALL salsa_nesting_offl_input
!
!--    Check if sufficient time steps are provided to cover the entire simulation. Note, dynamic
!--    input is only required for the 3D simulation, not for the soil/wall spinup. However, as the
!--    spinup time is added to the end_time, this must be considered here.
       IF ( end_time - spinup_time > salsa_nest_offl%time(salsa_nest_offl%nt-1) )  THEN
          message_string = 'end_time of the simulation exceeds the time dimension in the dynamic'//&
                           ' input file.'
          CALL message( 'salsa_nesting_offl_init', 'PA0690', 1, 2, 0, 6, 0 ) 
       ENDIF

       IF ( salsa_nest_offl%time(0) /= 0.0_wp )  THEN
          message_string = 'Offline nesting: time dimension must start at 0.0.'
          CALL message( 'salsa_nesting_offl_init', 'PA0691', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Initialize boundary data. Please note, do not initialize boundaries in case of restart runs.
       IF ( TRIM( initializing_actions ) /= 'read_restart_data'  .AND.  read_restart_data_salsa )  &
       THEN
          IF ( bc_dirichlet_l )  THEN
             DO  ib = 1, nbins_aerosol
                aerosol_number(ib)%conc(nzb+1:nzt,nys:nyn,-1) =                                    &
                                                 salsa_nest_offl%nconc_left(0,nzb+1:nzt,nys:nyn,ib)
                DO  ic = 1, ncomponents_mass
                   icc = ( ic - 1 ) * nbins_aerosol + ib
                   aerosol_mass(icc)%conc(nzb+1:nzt,nys:nyn,-1) =                                  &
                                                 salsa_nest_offl%mconc_left(0,nzb+1:nzt,nys:nyn,icc)
                ENDDO
             ENDDO
             DO  ig = 1, ngases_salsa
                salsa_gas(ig)%conc(nzb+1:nzt,nys:nyn,-1) =                                         &
                                                 salsa_nest_offl%gconc_left(0,nzb+1:nzt,nys:nyn,ig)
             ENDDO
          ENDIF
          IF ( bc_dirichlet_r )  THEN
             DO  ib = 1, nbins_aerosol
                aerosol_number(ib)%conc(nzb+1:nzt,nys:nyn,nxr+1) =                                 &
                                                salsa_nest_offl%nconc_right(0,nzb+1:nzt,nys:nyn,ib)
                DO  ic = 1, ncomponents_mass
                   icc = ( ic - 1 ) * nbins_aerosol + ib
                   aerosol_mass(icc)%conc(nzb+1:nzt,nys:nyn,nxr+1) =                               &
                                                salsa_nest_offl%mconc_right(0,nzb+1:nzt,nys:nyn,icc)
                ENDDO
             ENDDO
             DO  ig = 1, ngases_salsa
                salsa_gas(ig)%conc(nzb+1:nzt,nys:nyn,nxr+1) =                                      &
                                                 salsa_nest_offl%gconc_right(0,nzb+1:nzt,nys:nyn,ig)
             ENDDO
          ENDIF
          IF ( bc_dirichlet_n )  THEN
             DO  ib = 1, nbins_aerosol
                aerosol_number(ib)%conc(nzb+1:nzt,nyn+1,nxl:nxr) =                                 &
                                                salsa_nest_offl%nconc_north(0,nzb+1:nzt,nxl:nxr,ib)
                DO  ic = 1, ncomponents_mass
                   icc = ( ic - 1 ) * nbins_aerosol + ib
                   aerosol_mass(icc)%conc(nzb+1:nzt,nyn+1,nxl:nxr) =                               &
                                                salsa_nest_offl%mconc_north(0,nzb+1:nzt,nxl:nxr,icc)
                ENDDO
             ENDDO
             DO  ig = 1, ngases_salsa
                salsa_gas(ig)%conc(nzb+1:nzt,nyn+1,nxl:nxr) =                                      &
                                                 salsa_nest_offl%gconc_north(0,nzb+1:nzt,nxl:nxr,ig)
             ENDDO
          ENDIF
          IF ( bc_dirichlet_s )  THEN
             DO  ib = 1, nbins_aerosol
                aerosol_number(ib)%conc(nzb+1:nzt,-1,nxl:nxr) =                                    &
                                                salsa_nest_offl%nconc_south(0,nzb+1:nzt,nxl:nxr,ib)
                DO  ic = 1, ncomponents_mass
                   icc = ( ic - 1 ) * nbins_aerosol + ib
                   aerosol_mass(icc)%conc(nzb+1:nzt,-1,nxl:nxr) =                                  &
                                                salsa_nest_offl%mconc_south(0,nzb+1:nzt,nxl:nxr,icc)
                ENDDO
             ENDDO
             DO  ig = 1, ngases_salsa
                salsa_gas(ig)%conc(nzb+1:nzt,-1,nxl:nxr) =                                         &
                                                 salsa_nest_offl%gconc_south(0,nzb+1:nzt,nxl:nxr,ig)
             ENDDO
          ENDIF
       ENDIF
    ENDIF

 END SUBROUTINE salsa_nesting_offl_init

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Set the lateral and top boundary conditions in case the PALM domain is
!> nested offline in a mesoscale model. Further, average boundary data and
!> determine mean profiles, further used for correct damping in the sponge
!> layer.
!------------------------------------------------------------------------------!
 SUBROUTINE salsa_nesting_offl_input

    USE netcdf_data_input_mod,                                                                     &
        ONLY:  check_existence, close_input_file, get_attribute, get_variable,                     &
               inquire_num_variables, inquire_variable_names,                                      &
               get_dimension_length, open_read_file

    IMPLICIT NONE

    CHARACTER(LEN=25) ::  vname  !< variable name

    INTEGER(iwp) ::  ic        !< running index for aerosol chemical components
    INTEGER(iwp) ::  ig        !< running index for gases
    INTEGER(iwp) ::  num_vars  !< number of variables in netcdf input file

!
!-- Skip input if no forcing from larger-scale models is applied.
    IF ( .NOT. nesting_offline_salsa )  RETURN
!
!-- Initialise
    IF ( .NOT. salsa_nest_offl%init )  THEN

#if defined ( __netcdf )
!
!--    Open file in read-only mode
       CALL open_read_file( TRIM( input_file_dynamic ) // TRIM( coupling_char ),                   &
                            salsa_nest_offl%id_dynamic )
!
!--    At first, inquire all variable names.
       CALL inquire_num_variables( salsa_nest_offl%id_dynamic, num_vars )
!
!--    Allocate memory to store variable names.
       ALLOCATE( salsa_nest_offl%var_names(1:num_vars) )
       CALL inquire_variable_names( salsa_nest_offl%id_dynamic, salsa_nest_offl%var_names )
!
!--    Read time dimension, allocate memory and finally read time array
       CALL get_dimension_length( salsa_nest_offl%id_dynamic, salsa_nest_offl%nt,&
                                                    'time' )

       IF ( check_existence( salsa_nest_offl%var_names, 'time' ) )  THEN
          ALLOCATE( salsa_nest_offl%time(0:salsa_nest_offl%nt-1) )
          CALL get_variable( salsa_nest_offl%id_dynamic, 'time', salsa_nest_offl%time )
       ENDIF
!
!--    Read the vertical dimension
       CALL get_dimension_length( salsa_nest_offl%id_dynamic, salsa_nest_offl%nzu, 'z' )
       ALLOCATE( salsa_nest_offl%zu_atmos(1:salsa_nest_offl%nzu) )
       CALL get_variable( salsa_nest_offl%id_dynamic, 'z', salsa_nest_offl%zu_atmos )
!
!--    Read the number of aerosol chemical components
       CALL get_dimension_length( salsa_nest_offl%id_dynamic, salsa_nest_offl%ncc,                 &
                                  'composition_index' )
!
!--    Read the names of aerosol chemical components
       CALL get_variable( salsa_nest_offl%id_dynamic, 'composition_name', salsa_nest_offl%cc_name, &
                          salsa_nest_offl%ncc )
!
!--    Define the index of each chemical component in the model
       DO  ic = 1, salsa_nest_offl%ncc
          SELECT CASE ( TRIM( salsa_nest_offl%cc_name(ic) ) )
             CASE ( 'H2SO4', 'SO4', 'h2so4', 'so4' )
                salsa_nest_offl%cc_in2mod(1) = ic
             CASE ( 'OC', 'oc' )
                salsa_nest_offl%cc_in2mod(2) = ic
             CASE ( 'BC', 'bc' )
                salsa_nest_offl%cc_in2mod(3) = ic
             CASE ( 'DU', 'du' )
                salsa_nest_offl%cc_in2mod(4) = ic
             CASE ( 'SS', 'ss' )
                salsa_nest_offl%cc_in2mod(5) = ic
             CASE ( 'HNO3', 'hno3', 'NO3', 'no3', 'NO', 'no' )
                salsa_nest_offl%cc_in2mod(6) = ic
             CASE ( 'NH3', 'nh3', 'NH4', 'nh4', 'NH', 'nh' )
                salsa_nest_offl%cc_in2mod(7) = ic
          END SELECT
       ENDDO
       IF ( SUM( salsa_nest_offl%cc_in2mod ) == 0 )  THEN
          message_string = 'None of the aerosol chemical components in ' //                        &
                           TRIM( input_file_dynamic ) // ' correspond to ones applied in SALSA.'
          CALL message( 'salsa_mod: salsa_nesting_offl_input', 'PA0693', 2, 2, 0, 6, 0 )
       ENDIF
       
       CALL close_input_file( salsa_nest_offl%id_dynamic )
#endif
    ENDIF
!
!-- Check if dynamic driver data input is required.
    IF ( salsa_nest_offl%time(salsa_nest_offl%tind_p) <= MAX( time_since_reference_point, 0.0_wp)  &
         .OR.  .NOT.  salsa_nest_offl%init )  THEN
       CONTINUE
!
!-- Return otherwise
    ELSE
       RETURN
    ENDIF
!
!-- Obtain time index for current point in time.
    salsa_nest_offl%tind = MINLOC( ABS( salsa_nest_offl%time -                                     &
                                   MAX( time_since_reference_point, 0.0_wp ) ), DIM = 1 ) - 1
    salsa_nest_offl%tind_p = salsa_nest_offl%tind + 1
!
!-- Open file in read-only mode
#if defined ( __netcdf )

    CALL open_read_file( TRIM( input_file_dynamic ) // TRIM( coupling_char ),                      &
                         salsa_nest_offl%id_dynamic )
!
!-- Read data at the western boundary
    CALL get_variable( salsa_nest_offl%id_dynamic, 'ls_forcing_left_aerosol',                      &
                       salsa_nest_offl%nconc_left,                                                 &
                       MERGE( 0, 1, bc_dirichlet_l ), MERGE( nbins_aerosol-1, 0, bc_dirichlet_l ), &
                       MERGE( nys, 1, bc_dirichlet_l ), MERGE( nyn, 0, bc_dirichlet_l ),           &
                       MERGE( nzb, 1, bc_dirichlet_l ), MERGE( nzt-1, 0, bc_dirichlet_l ),         &
                       MERGE( salsa_nest_offl%tind,   1, bc_dirichlet_l ),                         &
                       MERGE( salsa_nest_offl%tind_p, 0, bc_dirichlet_l  ) )
    IF ( bc_dirichlet_l )  THEN
       salsa_nest_offl%nconc_left = MAX( nclim, salsa_nest_offl%nconc_left )
       CALL nesting_offl_aero_mass( salsa_nest_offl%tind, salsa_nest_offl%tind_p, nzb+1, nzt, nys, &
                                    nyn, 'ls_forcing_left_mass_fracs_a', 1 )
    ENDIF
    IF ( .NOT. salsa_gases_from_chem )  THEN
       DO  ig = 1, ngases_salsa
          vname = salsa_nest_offl%char_l // salsa_nest_offl%gas_name(ig)
          CALL get_variable( salsa_nest_offl%id_dynamic, TRIM( vname ),                            &
                             salsa_nest_offl%gconc_left(:,:,:,ig),                                 &
                             MERGE( nys, 1, bc_dirichlet_l ), MERGE( nyn, 0, bc_dirichlet_l ),     &
                             MERGE( nzb, 1, bc_dirichlet_l ), MERGE( nzt-1, 0, bc_dirichlet_l ),   &
                             MERGE( salsa_nest_offl%tind,   1, bc_dirichlet_l ),                   &
                             MERGE( salsa_nest_offl%tind_p, 0, bc_dirichlet_l ) )
          IF ( bc_dirichlet_l )  salsa_nest_offl%gconc_left(:,:,:,ig) =                            &
                                                  MAX( nclim, salsa_nest_offl%gconc_left(:,:,:,ig) )
       ENDDO
    ENDIF
!
!-- Read data at the eastern boundary
    CALL get_variable( salsa_nest_offl%id_dynamic, 'ls_forcing_right_aerosol',                     &
                       salsa_nest_offl%nconc_right,                                                &
                       MERGE( 0, 1, bc_dirichlet_r ), MERGE( nbins_aerosol-1, 0, bc_dirichlet_r ), &
                       MERGE( nys, 1, bc_dirichlet_r ), MERGE( nyn, 0, bc_dirichlet_r ),           &
                       MERGE( nzb, 1, bc_dirichlet_r ), MERGE( nzt-1, 0, bc_dirichlet_r ),         &
                       MERGE( salsa_nest_offl%tind,   1, bc_dirichlet_r ),                         &
                       MERGE( salsa_nest_offl%tind_p, 0, bc_dirichlet_r ) )
    IF ( bc_dirichlet_r )  THEN
       salsa_nest_offl%nconc_right = MAX( nclim, salsa_nest_offl%nconc_right )
       CALL nesting_offl_aero_mass( salsa_nest_offl%tind, salsa_nest_offl%tind_p, nzb+1, nzt, nys, &
                                    nyn, 'ls_forcing_right_mass_fracs_a', 2 )
    ENDIF
    IF ( .NOT. salsa_gases_from_chem )  THEN
       DO  ig = 1, ngases_salsa
          vname = salsa_nest_offl%char_r // salsa_nest_offl%gas_name(ig)
          CALL get_variable( salsa_nest_offl%id_dynamic, TRIM( vname ),                            &
                             salsa_nest_offl%gconc_right(:,:,:,ig),                                &
                             MERGE( nys, 1, bc_dirichlet_r ), MERGE( nyn, 0, bc_dirichlet_r ),     &
                             MERGE( nzb, 1, bc_dirichlet_r ), MERGE( nzt-1, 0, bc_dirichlet_r ),   &
                             MERGE( salsa_nest_offl%tind,   1, bc_dirichlet_r ),                   &
                             MERGE( salsa_nest_offl%tind_p, 0, bc_dirichlet_r ) )
          IF ( bc_dirichlet_r )  salsa_nest_offl%gconc_right(:,:,:,ig) =                           &
                                                 MAX( nclim, salsa_nest_offl%gconc_right(:,:,:,ig) )
       ENDDO
    ENDIF
!
!-- Read data at the northern boundary
    CALL get_variable( salsa_nest_offl%id_dynamic, 'ls_forcing_north_aerosol',                     &
                       salsa_nest_offl%nconc_north,                                                &
                       MERGE( 0, 1, bc_dirichlet_n ), MERGE( nbins_aerosol-1, 0, bc_dirichlet_n ), &
                       MERGE( nxl, 1, bc_dirichlet_n ), MERGE( nxr, 0, bc_dirichlet_n ),           &
                       MERGE( nzb, 1, bc_dirichlet_n ), MERGE( nzt-1, 0, bc_dirichlet_n ),         &
                       MERGE( salsa_nest_offl%tind,   1, bc_dirichlet_n ),                         &
                       MERGE( salsa_nest_offl%tind_p, 0, bc_dirichlet_n ) )
    IF ( bc_dirichlet_n )  THEN
       salsa_nest_offl%nconc_north = MAX( nclim, salsa_nest_offl%nconc_north )
       CALL nesting_offl_aero_mass( salsa_nest_offl%tind, salsa_nest_offl%tind_p, nzb+1, nzt, nxl, &
                                    nxr, 'ls_forcing_north_mass_fracs_a', 3 )
    ENDIF
    IF ( .NOT. salsa_gases_from_chem )  THEN
       DO  ig = 1, ngases_salsa
          vname = salsa_nest_offl%char_n // salsa_nest_offl%gas_name(ig)
          CALL get_variable( salsa_nest_offl%id_dynamic, TRIM( vname ),                            &
                             salsa_nest_offl%gconc_north(:,:,:,ig),                                &
                             MERGE( nxl, 1, bc_dirichlet_n ), MERGE( nxr, 0, bc_dirichlet_n ),     &
                             MERGE( nzb, 1, bc_dirichlet_n ), MERGE( nzt-1, 0, bc_dirichlet_n ),   &
                             MERGE( salsa_nest_offl%tind,   1, bc_dirichlet_n ),                   &
                             MERGE( salsa_nest_offl%tind_p, 0, bc_dirichlet_n ) )
          IF ( bc_dirichlet_n )  salsa_nest_offl%gconc_north(:,:,:,ig) =                           &
                                                 MAX( nclim, salsa_nest_offl%gconc_north(:,:,:,ig) )
       ENDDO
    ENDIF
!
!-- Read data at the southern boundary
    CALL get_variable( salsa_nest_offl%id_dynamic, 'ls_forcing_south_aerosol',                     &
                       salsa_nest_offl%nconc_south,                                                &
                       MERGE( 0, 1, bc_dirichlet_s ), MERGE( nbins_aerosol-1, 0, bc_dirichlet_s ), &
                       MERGE( nxl, 1, bc_dirichlet_s ), MERGE( nxr, 0, bc_dirichlet_s ),           &
                       MERGE( nzb, 1, bc_dirichlet_s ), MERGE( nzt-1, 0, bc_dirichlet_s ),         &
                       MERGE( salsa_nest_offl%tind,   1, bc_dirichlet_s ),                         &
                       MERGE( salsa_nest_offl%tind_p, 0, bc_dirichlet_s ) )
    IF ( bc_dirichlet_s )  THEN
       salsa_nest_offl%nconc_south = MAX( nclim, salsa_nest_offl%nconc_south )
       CALL nesting_offl_aero_mass( salsa_nest_offl%tind, salsa_nest_offl%tind_p, nzb+1, nzt, nxl, &
                                    nxr, 'ls_forcing_south_mass_fracs_a', 4 )
    ENDIF
    IF ( .NOT. salsa_gases_from_chem )  THEN
       DO  ig = 1, ngases_salsa
          vname = salsa_nest_offl%char_s // salsa_nest_offl%gas_name(ig)
          CALL get_variable( salsa_nest_offl%id_dynamic, TRIM( vname ),                            &
                             salsa_nest_offl%gconc_south(:,:,:,ig),                                &
                             MERGE( nxl, 1, bc_dirichlet_s ), MERGE( nxr, 0, bc_dirichlet_s ),     &
                             MERGE( nzb, 1, bc_dirichlet_s ), MERGE( nzt-1, 0, bc_dirichlet_s ),   &
                             MERGE( salsa_nest_offl%tind,   1, bc_dirichlet_s ),                   &
                             MERGE( salsa_nest_offl%tind_p, 0, bc_dirichlet_s ) )
          IF ( bc_dirichlet_s )  salsa_nest_offl%gconc_south(:,:,:,ig) =                           &
                                                 MAX( nclim, salsa_nest_offl%gconc_south(:,:,:,ig) )
       ENDDO
    ENDIF
!
!-- Read data at the top boundary
    CALL get_variable( salsa_nest_offl%id_dynamic, 'ls_forcing_top_aerosol',                       &
                       salsa_nest_offl%nconc_top(0:1,nys:nyn,nxl:nxr,1:nbins_aerosol),             &
                       0, nbins_aerosol-1, nxl, nxr, nys, nyn, salsa_nest_offl%tind,               &
                       salsa_nest_offl%tind_p )
    salsa_nest_offl%nconc_top = MAX( nclim, salsa_nest_offl%nconc_top )
    CALL nesting_offl_aero_mass( salsa_nest_offl%tind, salsa_nest_offl%tind_p, nys, nyn, nxl, nxr, &
                                 'ls_forcing_top_mass_fracs_a', 5 )
    IF ( .NOT. salsa_gases_from_chem )  THEN
       DO  ig = 1, ngases_salsa
          vname = salsa_nest_offl%char_t // salsa_nest_offl%gas_name(ig)
          CALL get_variable( salsa_nest_offl%id_dynamic, TRIM( vname ),                            &
                             salsa_nest_offl%gconc_top(:,:,:,ig), nxl, nxr, nys, nyn,              &
                             salsa_nest_offl%tind, salsa_nest_offl%tind_p )
          salsa_nest_offl%gconc_top(:,:,:,ig) = MAX( nclim, salsa_nest_offl%gconc_top(:,:,:,ig) )
       ENDDO
    ENDIF
!
!-- Close input file
    CALL close_input_file( salsa_nest_offl%id_dynamic )

#endif
!
!-- Set control flag to indicate that initialization is already done
    salsa_nest_offl%init = .TRUE.

 END SUBROUTINE salsa_nesting_offl_input

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Sets the mass concentrations to aerosol arrays in 2a and 2b.
!------------------------------------------------------------------------------!
 SUBROUTINE nesting_offl_aero_mass( ts, te, ks, ke, is, ie, varname_a, ibound )

    USE netcdf_data_input_mod,                                                                     &
        ONLY:  get_variable

    IMPLICIT NONE

    CHARACTER(LEN=25) ::  varname_b  !< name for bins b

    CHARACTER(LEN=*), INTENT(in) ::  varname_a  !< name for bins a

    INTEGER(iwp) ::  ee                !< loop index: end
    INTEGER(iwp) ::  i                 !< loop index
    INTEGER(iwp) ::  ib                !< loop index
    INTEGER(iwp) ::  ic                !< loop index
    INTEGER(iwp) ::  k                 !< loop index
    INTEGER(iwp) ::  ss                !< loop index: start
    INTEGER(iwp) ::  t                 !< loop index
    INTEGER(iwp) ::  type_so4_oc = -1  !<

    INTEGER(iwp), INTENT(in) ::  ibound  !< index: 1=left, 2=right, 3=north, 4=south, 5=top
    INTEGER(iwp), INTENT(in) ::  ie      !< loop index
    INTEGER(iwp), INTENT(in) ::  is      !< loop index
    INTEGER(iwp), INTENT(in) ::  ks      !< loop index
    INTEGER(iwp), INTENT(in) ::  ke      !< loop index
    INTEGER(iwp), INTENT(in) ::  ts      !< loop index
    INTEGER(iwp), INTENT(in) ::  te      !< loop index

    INTEGER(iwp), DIMENSION(maxspec) ::  cc_i2m   !<

    REAL(wp) ::  pmf1a !< mass fraction in 1a

    REAL(wp), DIMENSION(nbins_aerosol) ::  core   !< size of the bin mid aerosol particle

    REAL(wp), DIMENSION(0:1,ks:ke,is:ie,1:nbins_aerosol) ::  to_nconc                   !<
    REAL(wp), DIMENSION(0:1,ks:ke,is:ie,1:nbins_aerosol*ncomponents_mass) ::  to_mconc  !<

    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  mf2a !< Mass distributions for a
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  mf2b !< and b bins

!
!-- Variable name for insoluble mass fraction
    varname_b = varname_a(1:LEN( TRIM( varname_a ) ) - 1 ) // 'b'
!
!-- Bin mean aerosol particle volume (m3)
    core(1:nbins_aerosol) = api6 * aero(1:nbins_aerosol)%dmid**3
!
!-- Allocate and read mass fraction arrays
    ALLOCATE( mf2a(0:1,ks:ke,is:ie,1:salsa_nest_offl%ncc),                                         &
              mf2b(0:1,ks:ke,is:ie,1:salsa_nest_offl%ncc) )
    IF ( ibound == 5 )  THEN
       CALL get_variable( salsa_nest_offl%id_dynamic, varname_a,                                   &
                          mf2a(0:1,ks:ke,is:ie,1:salsa_nest_offl%ncc), 0, salsa_nest_offl%ncc-1,   &
                          is, ie, ks, ke, ts, te )
    ELSE
       CALL get_variable( salsa_nest_offl%id_dynamic, varname_a,                                   &
                          mf2a(0:1,ks:ke,is:ie,1:salsa_nest_offl%ncc), 0, salsa_nest_offl%ncc-1,   &
                          is, ie, ks-1, ke-1, ts, te )
    ENDIF
!
!-- If the chemical component is not activated, set its mass fraction to 0 to avoid mass inbalance
    cc_i2m = salsa_nest_offl%cc_in2mod
    IF ( index_so4 < 0  .AND. cc_i2m(1) > 0 )  mf2a(:,:,:,cc_i2m(1)) = 0.0_wp
    IF ( index_oc < 0   .AND. cc_i2m(2) > 0 )  mf2a(:,:,:,cc_i2m(2)) = 0.0_wp
    IF ( index_bc < 0   .AND. cc_i2m(3) > 0 )  mf2a(:,:,:,cc_i2m(3)) = 0.0_wp
    IF ( index_du < 0   .AND. cc_i2m(4) > 0 )  mf2a(:,:,:,cc_i2m(4)) = 0.0_wp
    IF ( index_ss < 0   .AND. cc_i2m(5) > 0 )  mf2a(:,:,:,cc_i2m(5)) = 0.0_wp
    IF ( index_no < 0   .AND. cc_i2m(6) > 0 )  mf2a(:,:,:,cc_i2m(6)) = 0.0_wp
    IF ( index_nh < 0   .AND. cc_i2m(7) > 0 )  mf2a(:,:,:,cc_i2m(7)) = 0.0_wp
    mf2b = 0.0_wp
!
!-- Initialise variable type_so4_oc to indicate whether SO4 and/OC is included in mass fraction data
    IF ( ( cc_i2m(1) > 0  .AND.  index_so4 > 0 )  .AND. ( cc_i2m(2) > 0  .AND.  index_oc > 0 ) )   &
    THEN
       type_so4_oc = 1
    ELSEIF ( cc_i2m(1) > 0  .AND.  index_so4 > 0 )  THEN
       type_so4_oc = 2
    ELSEIF ( cc_i2m(2) > 0  .AND.  index_oc > 0 )  THEN
       type_so4_oc = 3
    ENDIF

    SELECT CASE ( ibound )
       CASE( 1 )
          to_nconc = salsa_nest_offl%nconc_left
          to_mconc = salsa_nest_offl%mconc_left
       CASE( 2 )
          to_nconc = salsa_nest_offl%nconc_right
          to_mconc = salsa_nest_offl%mconc_right
       CASE( 3 )
          to_nconc = salsa_nest_offl%nconc_north
          to_mconc = salsa_nest_offl%mconc_north
       CASE( 4 )
          to_nconc = salsa_nest_offl%nconc_south
          to_mconc = salsa_nest_offl%mconc_south
       CASE( 5 )
          to_nconc = salsa_nest_offl%nconc_top
          to_mconc = salsa_nest_offl%mconc_top
    END SELECT
!
!-- Set mass concentrations:
!
!-- Regime 1:
    SELECT CASE ( type_so4_oc )
       CASE ( 1 )  ! Both SO4 and OC given

          ss = ( index_so4 - 1 ) * nbins_aerosol + start_subrange_1a  ! start
          ee = ( index_so4 - 1 ) * nbins_aerosol + end_subrange_1a    ! end
          ib = start_subrange_1a
          DO  ic = ss, ee
             DO i = is, ie
                DO k = ks, ke
                   DO t = 0, 1
                      pmf1a = mf2a(t,k,i,cc_i2m(1)) / ( mf2a(t,k,i,cc_i2m(1)) + mf2a(t,k,i,cc_i2m(2)) )
                      to_mconc(t,k,i,ic) = pmf1a * to_nconc(t,k,i,ib) * core(ib) * arhoh2so4
                   ENDDO
                ENDDO
             ENDDO
             ib = ib + 1
          ENDDO
          ss = ( index_oc - 1 ) * nbins_aerosol + start_subrange_1a ! start
          ee = ( index_oc - 1 ) * nbins_aerosol + end_subrange_1a   ! end
          ib = start_subrange_1a
          DO  ic = ss, ee
             DO i = is, ie
                DO k = ks, ke
                   DO t = 0, 1
                      pmf1a = mf2a(t,k,i,cc_i2m(2)) / ( mf2a(t,k,i,cc_i2m(1)) + mf2a(t,k,i,cc_i2m(2)) )
                      to_mconc(t,k,i,ic) = pmf1a * to_nconc(t,k,i,ib) * core(ib) * arhooc
                   ENDDO
                ENDDO
             ENDDO
             ib = ib + 1
          ENDDO
       CASE ( 2 )  ! Only SO4
          ss = ( index_so4 - 1 ) * nbins_aerosol + start_subrange_1a  ! start
          ee = ( index_so4 - 1 ) * nbins_aerosol + end_subrange_1a    ! end
          ib = start_subrange_1a
          DO  ic = ss, ee
             DO i = is, ie
                DO k = ks, ke
                   DO t = 0, 1
                      to_mconc(t,k,i,ic) = to_nconc(t,k,i,ib) * core(ib) * arhoh2so4
                   ENDDO
                ENDDO
             ENDDO
             ib = ib + 1
          ENDDO
       CASE ( 3 )  ! Only OC
          ss = ( index_oc - 1 ) * nbins_aerosol + start_subrange_1a ! start
          ee = ( index_oc - 1 ) * nbins_aerosol + end_subrange_1a   ! end
          ib = start_subrange_1a
          DO  ic = ss, ee
             DO i = is, ie
                DO k = ks, ke
                   DO t = 0, 1
                      to_mconc(t,k,i,ic) = to_nconc(t,k,i,ib) * core(ib) * arhooc
                   ENDDO
                ENDDO
             ENDDO
             ib = ib + 1
          ENDDO
    END SELECT
!
!-- Regimes 2a and 2b:
    IF ( index_so4 > 0 ) THEN
       CALL set_nest_mass( index_so4, 1, arhoh2so4 )
    ENDIF
    IF ( index_oc > 0 ) THEN
       CALL set_nest_mass( index_oc, 2, arhooc )
    ENDIF
    IF ( index_bc > 0 ) THEN
       CALL set_nest_mass( index_bc, 3, arhobc )
    ENDIF
    IF ( index_du > 0 ) THEN
       CALL set_nest_mass( index_du, 4, arhodu )
    ENDIF
    IF ( index_ss > 0 ) THEN
       CALL set_nest_mass( index_ss, 5, arhoss )
    ENDIF
    IF ( index_no > 0 ) THEN
       CALL set_nest_mass( index_no, 6, arhohno3 )
    ENDIF
    IF ( index_nh > 0 ) THEN
       CALL set_nest_mass( index_nh, 7, arhonh3 )
    ENDIF

    DEALLOCATE( mf2a, mf2b )

    SELECT CASE ( ibound )
       CASE( 1 )
          salsa_nest_offl%mconc_left = to_mconc
       CASE( 2 )
          salsa_nest_offl%mconc_right = to_mconc
       CASE( 3 )
          salsa_nest_offl%mconc_north = to_mconc
       CASE( 4 )
          salsa_nest_offl%mconc_south = to_mconc
       CASE( 5 )
          salsa_nest_offl%mconc_top = to_mconc
    END SELECT

    CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Set nesting boundaries for aerosol mass.
!------------------------------------------------------------------------------!
    SUBROUTINE set_nest_mass( ispec, ispec_def, prho )

       IMPLICIT NONE

       INTEGER(iwp) ::  ic   !< chemical component index: default
       INTEGER(iwp) ::  icc  !< loop index: mass bin

       INTEGER(iwp), INTENT(in) ::  ispec      !< aerosol species index
       INTEGER(iwp), INTENT(in) ::  ispec_def  !< default aerosol species index

       REAL(wp), INTENT(in) ::  prho !< aerosol density
!
!--    Define the index of the chemical component in the input data
       ic = salsa_nest_offl%cc_in2mod(ispec_def)

       DO i = is, ie
          DO k = ks, ke
             DO t = 0, 1
!
!--             Regime 2a:
                ss = ( ispec - 1 ) * nbins_aerosol + start_subrange_2a
                ee = ( ispec - 1 ) * nbins_aerosol + end_subrange_2a
                ib = start_subrange_2a
                DO icc = ss, ee
                   to_mconc(t,k,i,icc) = MAX( 0.0_wp, mf2a(t,k,i,ic) / SUM( mf2a(t,k,i,:) ) ) *    &
                                         to_nconc(t,k,i,ib) * core(ib) * prho
                   ib = ib + 1
                ENDDO
!
!--             Regime 2b:
                IF ( .NOT. no_insoluble )  THEN
!
!--                 TODO!
                    mf2b(t,k,i,ic) = mf2b(t,k,i,ic)
                ENDIF
             ENDDO   ! k

          ENDDO   ! j
       ENDDO   ! i

    END SUBROUTINE set_nest_mass

 END SUBROUTINE nesting_offl_aero_mass


 END MODULE salsa_mod
