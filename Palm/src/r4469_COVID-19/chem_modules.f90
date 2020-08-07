!> @file chem_modules.f90
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
! Copyright 2018-2019 Leibniz Universitaet Hannover
! Copyright 2018-2019 Karlsruhe Institute of Technology
! Copyright 2018-2019 Freie Universitaet Berlin
!------------------------------------------------------------------------------!
!
! Current revisions:
! -----------------
!
! 
! Former revisions:
! -----------------
! $Id: chem_modules.f90 4372 2020-01-14 10:20:35Z banzhafs $
! added namelist flag 'emiss_read_legacy_mode' to allow concurrent
! functioning of new emission read mode under development (ECC)
!
! 4273 2019-10-24 13:40:54Z monakurppa 
! Add logical switches nesting_chem and nesting_offline_chem (both .TRUE.
! by default)
!
! 4182 2019-08-22 15:20:23Z scharf
! Corrected "Former revisions" section
! 
! 4110 2019-07-22 17:05:21Z suehring
! +cs_advc_flags_s
! 
! 4109 2019-07-22 17:00:34Z suehring
! - introduced namelist item chem_modules@emiss_lod as future
! - replacement to chem_modules@mode_emis.  Currently keeping both
!   for backward compatibility.  chem_modules@mode_emis will be
!   depreciated upon migration of all dependent modules (e.g., salsa)
!   to chem_modules@emiss_lod
!
! (ecc) 20190513 replaced nspec_out with n_matched_vars
! 
! 3877 2019-04-08 19:09:16Z knoop
! Formatting, clean-up, clarified/corrected comments
! 
! 3833 2019-03-28 15:04:04Z forkel
! removed USE chem_gasphase_mod
! 
! 3827 2019-03-27 17:20:32Z forkel
! some formatting  and reordering (ecc) 
! 
! 3820 2019-03-27 11:53:41Z forkel
! renamed do_emis to emissions_anthropogenic, removed USE statistics, variables sorted by type 
!
! 
! 3780 2019-03-05 11:19:45Z forkel
! added cs_mech
! 
! 3652 2019-01-07 15:29:59Z forkel
! parameter chem_mechanism added (basit)
! 
! 3282 2018-09-27 10:49:12Z basit
! Initial revision
!
! Authors:
! --------
! @author Farah Kanani-Suehring
! @author Basit Khan
! @author Sabine Banzhaf
! @author Emmanuele Russo
! @author Edward C. Chan
!
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Definition of global PALM-4U chemistry variables
!------------------------------------------------------------------------------!
!
 MODULE chem_modules

    USE kinds

    IMPLICIT NONE

    CHARACTER (LEN=20) ::  bc_cs_b        = 'dirichlet'         !< namelist parameter: surface boundary condition for concentration
    CHARACTER (LEN=20) ::  bc_cs_t        = 'initial_gradient'  !< namelist parameter: top boudary condition for concentration
    CHARACTER (LEN=30) ::  chem_mechanism = 'phstatp'           !< namelist parameter: chemistry mechanism
    CHARACTER (LEN=80) ::  daytype_mdh    = 'workday'           !< namelist parameter: type of day - workday, weekend, holiday
    CHARACTER (LEN=80) ::  mode_emis      = 'PARAMETERIZED'     !< namelist parameter: mode of chemistry emissions - DEFAULT, EXPERT, PARAMETERIZED
    CHARACTER (LEN=80) ::  time_fac_type  = 'MDH'               !< namelist parameter: type of time treatment in the mode_emis DEFAULT - HOUR, MDH
    CHARACTER (LEN=10) ::  photolysis_scheme                    !< 'constant',
                                                                !< 'simple' (Simple parameterisation from MCM, Saunders et al., 2003, Atmos. Chem. Phys., 3, 161-180
                                                                !< 'fastj'  (Wild et al., 2000, J. Atmos. Chem., 37, 245-282) STILL NOT IMPLEMENTED

    CHARACTER (LEN=11), DIMENSION(99) ::  cs_name             = 'novalue'  !< namelist parameter: names of species with given fluxes (see csflux)
    CHARACTER (LEN=11), DIMENSION(99) ::  cs_profile_name     = 'novalue'  !< namelist parameter: tbc...???
    CHARACTER (LEN=11), DIMENSION(99) ::  data_output_pr_cs   = 'novalue'  !< namelist parameter: tbc...???
    CHARACTER (LEN=11), DIMENSION(99) ::  surface_csflux_name = 'novalue'  !< namelist parameter: tbc...???

    INTEGER(iwp) ::  cs_pr_count                           = 0      !< counter for chemical species profiles
    INTEGER(iwp) ::  cs_vertical_gradient_level_ind(99,10) = -9999  !< grid index values of cs_vertical_gradient_level
    INTEGER(iwp) ::  emiss_lod                             = -1     !< namelist parameter: chem emission LOD (same as mode_emis)
                                                                    !< -1 = unassigned, 0 = parameterized, 1 = default, 2 = pre-processed
    INTEGER(iwp) ::  ibc_cs_b                                       !< integer flag for bc_cs_b
    INTEGER(iwp) ::  ibc_cs_t                                       !< integer flag for bc_cs_t
    INTEGER(iwp) ::  main_street_id                        = 0      !< namelist parameter: lower bound of main street IDs (OpenStreetMaps) for PARAMETERIZED mode
    INTEGER(iwp) ::  max_pr_cs                             = 0      !< 
    INTEGER(iwp) ::  max_street_id                         = 0      !< namelist parameter: upper bound of main street IDs (OpenStreetMaps) for PARAMETERIZED mode      
    INTEGER(iwp) ::  n_matched_vars                                 !< number of matched emissions variables
    INTEGER(iwp) ::  side_street_id                        = 0      !< namelist parameter: lower bound of side street IDs (OpenStreetMaps) for PARAMETERIZED mode

    INTEGER(iwp), DIMENSION(99) ::  cs_pr_index  = 0      !< index for chemical species profiles
    INTEGER(iwp), DIMENSION(:)  ::  match_spec_nox(1:2)   !< results of matching the input and model's NOx
    INTEGER(iwp), DIMENSION(:)  ::  match_spec_pm(1:3)    !< results of matching the input and model's PMs
    INTEGER(iwp), DIMENSION(:)  ::  match_spec_sox(1:2)   !< results of matching the input and model's SOx! 

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  match_spec_input      !< index of input chem species for matching routine 
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  match_spec_model      !< index of model chem species for matching routine
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  match_spec_voc_input  !< index of VOC input components matching the model's VOCs
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  match_spec_voc_model  !< index of VOC model species matching the input VOCs comp.
    
    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE  ::  cs_advc_flags_s !< flags used to degrade order of advection scheme for 
                                                                     !< chemical species near walls and lateral boundaries

    LOGICAL ::  constant_top_csflux(99)   = .TRUE.   !< internal flag, set to .FALSE. if no top_csflux is prescribed
    LOGICAL ::  constant_csflux(99)       = .TRUE.   !< internal flag, set to .FALSE. if no surface_csflux is prescribed
    LOGICAL ::  call_chem_at_all_substeps = .FALSE.  !< namelist parameter: ....??? 
    LOGICAL ::  chem_debug0               = .FALSE.  !< namelist parameter: flag for minimum print output
    LOGICAL ::  chem_debug1               = .FALSE.  !< namelist parameter: flag for print output
    LOGICAL ::  chem_debug2               = .FALSE.  !< namelist parameter: flag for further print output
    LOGICAL ::  chem_gasphase_on          = .TRUE.   !< namelist parameter: flag to switch off chemical reactions
    LOGICAL ::  cs_pr_namelist_found      = .FALSE.  !< ...???
    LOGICAL ::  deposition_dry            = .FALSE.  !< namelist parameter: flag for activation of deposition calculation
    LOGICAL ::  emissions_anthropogenic   = .FALSE.  !< namelist parameter: flag for turning on anthropogenic emissions
    LOGICAL ::  emission_output_required  = .TRUE.   !< internal flag for requiring emission outputs
    LOGICAL ::  emiss_read_legacy_mode    = .TRUE.   !< namelist parameter: flag to read emission data using legacy mode
    LOGICAL ::  nesting_chem              = .TRUE.   !< apply self-nesting for the chemistry model
    LOGICAL ::  nesting_offline_chem      = .TRUE.   !< apply offline nesting for the chemistry model

    REAL(wp) ::  cs_surface_initial_change(99)     = 0.0_wp        !< namelist parameter: ...???
    REAL(wp) ::  cs_vertical_gradient(99,10)       = 0.0_wp        !< namelist parameter: ...???
    REAL(wp) ::  cs_vertical_gradient_level(99,10) = -999999.9_wp  !< namelist parameter: ...???
    REAL(wp) ::  emiss_factor_main(99)             = -9999.0_wp    !< namelist parameter: ...???
    REAL(wp) ::  emiss_factor_side(99)             = -9999.0_wp    !< namelist parameter: ...???
    REAL(wp) ::  surface_csflux(99)                = 0.0_wp        !< namelist parameter: ...???
    REAL(wp) ::  top_csflux(99)                    = 0.0_wp        !< namelist parameter: ...???
    REAL(wp) ::  wall_csflux(99,0:5)               = 0.0_wp        !< namelist parameter: ...???

    REAL(wp), DIMENSION(99)     ::  cs_surface = 0.0_wp        !< namelist parameter: chem species concentration at surface
    REAL(wp), DIMENSION(99,100) ::  cs_heights = 9999999.9_wp  !< namelist parameter: height levels for initial chem species concentrations
    REAL(wp), DIMENSION(99,100) ::  cs_profile = 9999999.9_wp  !< namelist parameter: chem species concentration values at cs_heights levels

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  bc_cs_t_val  !< vertical gradient of chemical species near domain top
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  css          !< scaling parameter for chem species

    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  emis_distribution  !< emissions final values (main module output) ???
                                  
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET ::  cs_1  !< pointer for swapping of timelevels for respective quantity
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET ::  cs_2  !< pointer for swapping of timelevels for respective quantity
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET ::  cs_3  !< pointer for swapping of timelevels for respective quantity

    REAL(wp), DIMENSION(:,:,:), POINTER ::  cs     !< pointer: sgs chem spcs  ???
    REAL(wp), DIMENSION(:,:,:), POINTER ::  cs_p   !< pointer: prognostic value of sgs chem spcs ???
    REAL(wp), DIMENSION(:,:,:), POINTER ::  tcs_m  !< pointer: to tcs array (temp)

    REAL, PARAMETER ::  xm_air   =   28.964e-3             !< air      molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_C     =   12.01115e-3           !< C        molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_Ca    =   40.07800e-3           !< Ca       molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_Cd    =  112.41000e-3           !< Cd       molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_Cl    =   35.45300e-3           !< Cl       molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_dummy = 1000.0e-3               !< dummy    molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_F     =   18.99840e-3           !< F        molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_H     =    1.00790e-3           !< H        molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_K     =   39.09800e-3           !< K        molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_Mg    =   24.30500e-3           !< Mg       molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_N     =   14.00670e-3           !< N        molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_Na    =   22.98977e-3           !< Na       molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_O     =   15.99940e-3           !< O        molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_Pb    =  207.20000e-3           !< Pb       molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_Pb210 =  210.00000e-3           !< Pb (210) molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_Rn222 =  222.00000e-3           !< Rn (222) molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_S     =   32.06400e-3           !< S        molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_CO2   = xm_C + xm_O * 2         !< CO2      molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_h2o   = xm_H * 2 + xm_O         !< H2O      molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_HNO3  = xm_H + xm_N + xm_O * 3  !< HNO3     molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_o3    = xm_O * 3                !< O3       molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_N2O5  = xm_N * 2 + xm_O * 5     !< N2O5     molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_NH4   = xm_N + xm_H * 4         !< NH4      molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_NO3   = xm_N + xm_O * 3         !< NO3      molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_SO4   = xm_S + xm_O * 4         !< SO4      molecular weight (kg/mol)
!
!-  Define chemical variables within chem_species
    TYPE species_def
       CHARACTER(LEN=15)                            ::  name         !< name of chemical species
       CHARACTER(LEN=15)                            ::  unit         !< unit (ppm for gases, kg m^-3 for aerosol tracers)
       REAL(kind=wp), POINTER, DIMENSION(:,:,:)     ::  conc         !< concentrations of trace gases
       REAL(kind=wp), POINTER, DIMENSION(:,:,:)     ::  conc_av      !< averaged concentrations
       REAL(kind=wp), POINTER, DIMENSION(:,:,:)     ::  conc_p       !< conc at prognostic time level
       REAL(kind=wp), POINTER, DIMENSION(:,:,:)     ::  tconc_m      !< weighted tendency of conc for previous sub-timestep (Runge-Kutta)
       REAL(kind=wp), ALLOCATABLE, DIMENSION(:,:)   ::  cssws_av     !< averaged fluxes of trace gases at surface
       REAL(kind=wp), ALLOCATABLE, DIMENSION(:,:)   ::  flux_s_cs    !< 6th-order advective flux at south face of grid box of chemical species (='cs')
       REAL(kind=wp), ALLOCATABLE, DIMENSION(:,:)   ::  diss_s_cs    !< artificial numerical dissipation flux at south face of grid box of chemical species
       REAL(kind=wp), ALLOCATABLE, DIMENSION(:,:,:) ::  flux_l_cs    !< 6th-order advective flux at left face of grid box of chemical species (='cs')
       REAL(kind=wp), ALLOCATABLE, DIMENSION(:,:,:) ::  diss_l_cs    !< artificial numerical dissipation flux at left face of grid box of chemical species
       REAL(kind=wp), ALLOCATABLE, DIMENSION(:)     ::  conc_pr_init !< initial profile of chemical species
    END TYPE species_def
!
!-- Define photolysis frequencies in phot_frequen
    TYPE photols_def
       CHARACTER(LEN=15)                            :: name          !< name of pgotolysis frequency
       CHARACTER(LEN=15)                            :: unit          !< unit (1/s)
       REAL(kind=wp), POINTER, DIMENSION(:,:,:)     :: freq          !< photolysis frequency
    END TYPE photols_def


    TYPE(species_def), ALLOCATABLE, DIMENSION(:), TARGET ::  chem_species
    TYPE(photols_def), ALLOCATABLE, DIMENSION(:), TARGET ::  phot_frequen

    SAVE

 END MODULE chem_modules
