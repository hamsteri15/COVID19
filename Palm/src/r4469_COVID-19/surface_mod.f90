!> @file surface_mod.f90
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
!
!------------------------------------------------------------------------------!
!
! Current revisions:
! ------------------
! 
! 
! Former revisions:
! -----------------
! $Id: surface_mod.f90 4366 2020-01-09 08:12:43Z raasch $
! workaround implemented to avoid vectorization bug on NEC Aurora
! 
! 4360 2020-01-07 11:25:50Z suehring
! Fix also remaining message calls.
! 
! 4354 2019-12-19 16:10:18Z suehring
! Bugfix in message call and specify error number
! 
! 4346 2019-12-18 11:55:56Z motisi
! Introduction of wall_flags_total_0, which currently sets bits based on static
! topography information used in wall_flags_static_0
! 
! 4331 2019-12-10 18:25:02Z suehring
! -pt_2m - array is moved to diagnostic_output_quantities
! 
! 4329 2019-12-10 15:46:36Z motisi
! Renamed wall_flags_0 to wall_flags_static_0
! 
! 4245 2019-09-30 08:40:37Z pavelkrc
! Corrected "Former revisions" section
! 
! 4168 2019-08-16 13:50:17Z suehring
! Remove functions get_topography_top_index. These are now replaced by 
! precalculated arrays because of too much CPU-time consumption
! 
! 4159 2019-08-15 13:31:35Z suehring
! Surface classification revised and adjusted to changes in init_grid
! 
! 4156 2019-08-14 09:18:14Z schwenkel
! Bugfix in case of cloud microphysics morrison
! 
! 4150 2019-08-08 20:00:47Z suehring
! Generic routine to initialize single surface properties added
! 
! 4104 2019-07-17 17:08:20Z suehring
! Bugfix, initialization of index space for boundary data structure accidantly
! run over ghost points, causing a segmentation fault.
! 
! 3943 2019-05-02 09:50:41Z maronga
! - Revise initialization of the boundary data structure
! - Add new data structure to set boundary conditions at vertical walls
! 
! 3943 2019-05-02 09:50:41Z maronga
! Removed qsws_eb as it is no longer needed.
! 
! 3933 2019-04-25 12:33:20Z kanani
! Add (de)allocation of pt_2m,
! bugfix: initialize pt_2m
! 
! 3833 2019-03-28 15:04:04Z forkel
! added USE chem_gasphase_mod (chem_modules will not transport nvar and nspec anymore)
! 
! 3772 2019-02-28 15:51:57Z suehring
! small change in the todo's
! 
! 3767 2019-02-27 08:18:02Z raasch
! unused variables removed from rrd-subroutine parameter list
! 
! 3761 2019-02-25 15:31:42Z raasch
! OpenACC directives added to avoid compiler warnings about unused variables,
! unused variable removed
! 
! 3745 2019-02-15 18:57:56Z suehring
! +waste_heat
! 
! 3744 2019-02-15 18:38:58Z suehring
! OpenACC port for SPEC
! 
! 2233 2017-05-30 18:08:54Z suehring
! Initial revision
!
!
! Description:
! ------------
!> Surface module defines derived data structures to treat surface-
!> bounded grid cells. Three different types of surfaces are defined: 
!> default surfaces, natural surfaces, and urban surfaces. The module
!> encompasses the allocation and initialization of surface arrays, and handles 
!> reading and writing restart data. 
!> In addition, a further derived data structure is defined, in order to set
!> boundary conditions at surfaces.  
!> @todo For the moment, downward-facing surfaces are only classified as  
!>        default type 
!> @todo Clean up urban-surface variables (some of them are not used any more) 
!> @todo Revise initialization of surface fluxes (especially for chemistry) 
!> @todo Get rid-off deallocation routines in restarts
!------------------------------------------------------------------------------!
 MODULE surface_mod

    USE arrays_3d,                                                             &
        ONLY:  heatflux_input_conversion, momentumflux_input_conversion,       &
               rho_air, rho_air_zw, zu, zw, waterflux_input_conversion 

    USE chem_gasphase_mod,                                                     &
        ONLY:  nvar, spc_names

    USE chem_modules

    USE control_parameters

    USE indices,                                                               &
        ONLY:  nxl, nxlg, nxr, nxrg, nys, nysg, nyn, nyng, nzb, nzt,           &
               wall_flags_total_0

    USE grid_variables,                                                        &
        ONLY:  dx, dy

    USE kinds

    USE model_1d_mod,                                                          &
        ONLY:  rif1d, us1d, usws1d, vsws1d


    IMPLICIT NONE

!
!-- Data type used to identify grid-points where horizontal boundary conditions 
!-- are applied 
    TYPE bc_type
       INTEGER(iwp) :: ioff !< offset value in x-direction, used to determine index of surface element
       INTEGER(iwp) :: joff !< offset value in y-direction, used to determine index of surface element
       INTEGER(iwp) :: koff !< offset value in z-direction, used to determine index of surface element
       INTEGER(iwp) :: ns   !< number of surface elements on the PE

       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  i !< x-index linking to the PALM 3D-grid  
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  j !< y-index linking to the PALM 3D-grid    
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  k !< z-index linking to the PALM 3D-grid   

       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE :: start_index !< start index within surface data type for given (j,i) 
       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE :: end_index   !< end index within surface data type for given (j,i)  

    END TYPE bc_type
!
!-- Data type used to identify and treat surface-bounded grid points  
    TYPE surf_type

       LOGICAL ::  albedo_from_ascii = .FALSE. !< flag indicating that albedo for urban surfaces is input via ASCII format (just for a workaround)
    
       INTEGER(iwp) :: ioff                                !< offset value in x-direction, used to determine index of surface element
       INTEGER(iwp) :: joff                                !< offset value in y-direction, used to determine index of surface element
       INTEGER(iwp) :: koff                                !< offset value in z-direction, used to determine index of surface element
       INTEGER(iwp) :: ns                                  !< number of surface elements on the PE

       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  i       !< x-index linking to the PALM 3D-grid  
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  j       !< y-index linking to the PALM 3D-grid    
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  k       !< z-index linking to the PALM 3D-grid       

       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  facing  !< Bit indicating surface orientation 
     
       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE :: start_index !< Start index within surface data type for given (j,i) 
       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE :: end_index   !< End index within surface data type for given (j,i)  

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  z_mo      !< surface-layer height 
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  uvw_abs   !< absolute surface-parallel velocity
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  us        !< friction velocity
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ts        !< scaling parameter temerature
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qs        !< scaling parameter humidity
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ss        !< scaling parameter passive scalar
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qcs       !< scaling parameter qc
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ncs       !< scaling parameter nc
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qrs       !< scaling parameter qr
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  nrs       !< scaling parameter nr

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ol        !< Obukhov length
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rib       !< Richardson bulk number

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  z0        !< roughness length for momentum
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  z0h       !< roughness length for heat
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  z0q       !< roughness length for humidity

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  pt1       !< potential temperature at first grid level
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qv1       !< mixing ratio at first grid level
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  vpt1      !< virtual potential temperature at first grid level
       
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  css     !< scaling parameter chemical species
!
!--    Define arrays for surface fluxes
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  usws      !< vertical momentum flux for u-component at horizontal surfaces
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  vsws      !< vertical momentum flux for v-component at horizontal surfaces 

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  shf       !< surface flux sensible heat
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qsws      !< surface flux latent heat 
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ssws      !< surface flux passive scalar
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qcsws     !< surface flux qc
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ncsws     !< surface flux nc
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qrsws     !< surface flux qr
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  nrsws     !< surface flux nr
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  sasws     !< surface flux salinity
!--    Added for SALSA:
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  answs   !< surface flux aerosol number: dim 1: flux, dim 2: bin
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  amsws   !< surface flux aerosol mass:   dim 1: flux, dim 2: bin
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  gtsws   !< surface flux gesous tracers: dim 1: flux, dim 2: gas
       
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  cssws   !< surface flux chemical species
!
!--    Required for horizontal walls in production_e
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  u_0       !< virtual velocity component (see production_e_init for further explanation)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  v_0       !< virtual velocity component (see production_e_init for further explanation)

       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  mom_flux_uv  !< momentum flux usvs and vsus at vertical surfaces (used in diffusion_u and diffusion_v)
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  mom_flux_w   !< momentum flux wsus and wsvs at vertical surfaces (used in diffusion_w)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  mom_flux_tke !< momentum flux usvs, vsus, wsus, wsvs at vertical surfaces at grid center (used in production_e)
!
!--    Variables required for LSM as well as for USM
       CHARACTER(LEN=40), DIMENSION(:), ALLOCATABLE ::  building_type_name    !< building type name at surface element
       CHARACTER(LEN=40), DIMENSION(:), ALLOCATABLE ::  pavement_type_name    !< pavement type name at surface element
       CHARACTER(LEN=40), DIMENSION(:), ALLOCATABLE ::  vegetation_type_name  !< water type at name surface element
       CHARACTER(LEN=40), DIMENSION(:), ALLOCATABLE ::  water_type_name       !< water type at name surface element
       
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  nzt_pavement     !< top index for pavement in soil
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  building_type    !< building type at surface element
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  pavement_type    !< pavement type at surface element
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  vegetation_type  !< vegetation type at surface element
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  water_type       !< water type at surface element
       
       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  albedo_type   !< albedo type, for each fraction (wall,green,window or vegetation,pavement water)

       LOGICAL, DIMENSION(:), ALLOCATABLE  ::  building_surface    !< flag parameter indicating that the surface element is covered by buildings (no LSM actions, not implemented yet)
       LOGICAL, DIMENSION(:), ALLOCATABLE  ::  building_covered    !< flag indicating that buildings are on top of orography, only used for vertical surfaces in LSM
       LOGICAL, DIMENSION(:), ALLOCATABLE  ::  pavement_surface    !< flag parameter for pavements
       LOGICAL, DIMENSION(:), ALLOCATABLE  ::  water_surface       !< flag parameter for water surfaces
       LOGICAL, DIMENSION(:), ALLOCATABLE  ::  vegetation_surface  !< flag parameter for natural land surfaces
       
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  albedo            !< broadband albedo for each surface fraction (LSM: vegetation, water, pavement; USM: wall, green, window)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  emissivity        !< emissivity of the surface, for each fraction  (LSM: vegetation, water, pavement; USM: wall, green, window)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  frac              !< relative surface fraction (LSM: vegetation, water, pavement; USM: wall, green, window)

       REAL(wp), DIMENSION(:,:), ALLOCATABLE   ::  aldif           !< albedo for longwave diffusive radiation, solar angle of 60 degrees
       REAL(wp), DIMENSION(:,:), ALLOCATABLE   ::  aldir           !< albedo for longwave direct radiation, solar angle of 60 degrees
       REAL(wp), DIMENSION(:,:), ALLOCATABLE   ::  asdif           !< albedo for shortwave diffusive radiation, solar angle of 60 degrees
       REAL(wp), DIMENSION(:,:), ALLOCATABLE   ::  asdir           !< albedo for shortwave direct radiation, solar angle of 60 degrees
       REAL(wp), DIMENSION(:,:), ALLOCATABLE   ::  rrtm_aldif      !< albedo for longwave diffusive radiation, solar angle of 60 degrees
       REAL(wp), DIMENSION(:,:), ALLOCATABLE   ::  rrtm_aldir      !< albedo for longwave direct radiation, solar angle of 60 degrees
       REAL(wp), DIMENSION(:,:), ALLOCATABLE   ::  rrtm_asdif      !< albedo for shortwave diffusive radiation, solar angle of 60 degrees
       REAL(wp), DIMENSION(:,:), ALLOCATABLE   ::  rrtm_asdir      !< albedo for shortwave direct radiation, solar angle of 60 degrees

       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  q_surface         !< skin-surface mixing ratio
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  pt_surface        !< skin-surface temperature
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  vpt_surface       !< skin-surface virtual temperature
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  rad_net           !< net radiation 
       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  rad_net_l         !< net radiation, used in USM
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  lambda_h          !< heat conductivity of soil/ wall (W/m/K) 
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  lambda_h_green    !< heat conductivity of green soil (W/m/K) 
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  lambda_h_window   !< heat conductivity of windows (W/m/K) 
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  lambda_h_def      !< default heat conductivity of soil (W/m/K)   

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_lw_in           !< incoming longwave radiation
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_lw_out          !< emitted longwave radiation
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_lw_dif          !< incoming longwave radiation from sky
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_lw_ref          !< incoming longwave radiation from reflection
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_lw_res          !< resedual longwave radiation in surface after last reflection step
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_in           !< incoming shortwave radiation
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_out          !< emitted shortwave radiation
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_dir          !< direct incoming shortwave radiation
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_dif          !< diffuse incoming shortwave radiation
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_ref          !< incoming shortwave radiation from reflection
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_sw_res          !< resedual shortwave radiation in surface after last reflection step

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  c_liq               !< liquid water coverage (of vegetated area)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  c_veg               !< vegetation coverage
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  f_sw_in             !< fraction of absorbed shortwave radiation by the surface layer (not implemented yet)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ghf                 !< ground heat flux
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  g_d                 !< coefficient for dependence of r_canopy on water vapour pressure deficit
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  lai                 !< leaf area index
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  lambda_surface_u    !< coupling between surface and soil (depends on vegetation type) (W/m2/K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  lambda_surface_s    !< coupling between surface and soil (depends on vegetation type) (W/m2/K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qsws_liq            !< surface flux of latent heat (liquid water portion)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qsws_soil           !< surface flux of latent heat (soil portion)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qsws_veg            !< surface flux of latent heat (vegetation portion)

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  r_a                 !< aerodynamic resistance 
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  r_a_green           !< aerodynamic resistance at green fraction
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  r_a_window          !< aerodynamic resistance at window fraction
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  r_canopy            !< canopy resistance
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  r_soil              !< soil resistance
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  r_soil_min          !< minimum soil resistance
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  r_s                 !< total surface resistance (combination of r_soil and r_canopy)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  r_canopy_min        !< minimum canopy (stomatal) resistance
       
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  pt_10cm             !< near surface air potential temperature at distance 10 cm from the surface (K)
       
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  alpha_vg          !< coef. of Van Genuchten
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  lambda_w          !< hydraulic diffusivity of soil (?)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  gamma_w           !< hydraulic conductivity of soil (W/m/K)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  gamma_w_sat       !< hydraulic conductivity at saturation
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  l_vg              !< coef. of Van Genuchten
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  m_fc              !< soil moisture at field capacity (m3/m3)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  m_res             !< residual soil moisture
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  m_sat             !< saturation soil moisture (m3/m3)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  m_wilt            !< soil moisture at permanent wilting point (m3/m3)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  n_vg              !< coef. Van Genuchten  
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rho_c_total_def   !< default volumetric heat capacity of the (soil) layer (J/m3/K)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rho_c_total       !< volumetric heat capacity of the actual soil matrix (J/m3/K)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  root_fr           !< root fraction within the soil layers
       
!--    Indoor model variables
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  waste_heat          !< waste heat 
!
!--    Urban surface variables
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  surface_types   !< array of types of wall parameters

       LOGICAL, DIMENSION(:), ALLOCATABLE  ::  isroof_surf          !< flag indicating roof surfaces
       LOGICAL, DIMENSION(:), ALLOCATABLE  ::  ground_level         !< flag indicating ground floor level surfaces 

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  target_temp_summer  !< indoor target temperature summer
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  target_temp_winter  !< indoor target temperature summer        

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  c_surface           !< heat capacity of the wall surface skin (J/m2/K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  c_surface_green     !< heat capacity of the green surface skin (J/m2/K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  c_surface_window    !< heat capacity of the window surface skin (J/m2/K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  green_type_roof     !< type of the green roof
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  lambda_surf         !< heat conductivity between air and surface (W/m2/K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  lambda_surf_green   !< heat conductivity between air and green surface (W/m2/K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  lambda_surf_window  !< heat conductivity between air and window surface (W/m2/K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  thickness_wall      !< thickness of the wall, roof and soil layers
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  thickness_green     !< thickness of the green wall, roof and soil layers
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  thickness_window    !< thickness of the window wall, roof and soil layers
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  transmissivity      !< transmissivity of windows

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfoutsl           !< reflected shortwave radiation for local surface in i-th reflection
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfoutll           !< reflected + emitted longwave radiation for local surface in i-th reflection
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfhf              !< total radiation flux incoming to minus outgoing from local surface

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  tt_surface_wall_m   !< surface temperature tendency (K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  tt_surface_window_m !< window surface temperature tendency (K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  tt_surface_green_m  !< green surface temperature tendency (K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  wshf                !< kinematic wall heat flux of sensible heat (actually no longer needed)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  wshf_eb             !< wall heat flux of sensible heat in wall normal direction

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  wghf_eb             !< wall ground heat flux
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  wghf_eb_window      !< window ground heat flux
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  wghf_eb_green       !< green ground heat flux
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  iwghf_eb            !< indoor wall ground heat flux
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  iwghf_eb_window     !< indoor window ground heat flux

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_lw_out_change_0

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfinsw            !< shortwave radiation falling to local surface including radiation from reflections
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfoutsw           !< total shortwave radiation outgoing from nonvirtual surfaces surfaces after all reflection
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfinlw            !< longwave radiation falling to local surface including radiation from reflections
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfoutlw           !< total longwave radiation outgoing from nonvirtual surfaces surfaces after all reflection
       
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  n_vg_green          !< vangenuchten parameters
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  alpha_vg_green      !< vangenuchten parameters
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  l_vg_green          !< vangenuchten parameters


       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rho_c_wall        !< volumetric heat capacity of the material ( J m-3 K-1 ) (= 2.19E6)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  dz_wall           !< wall grid spacing (center-center)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ddz_wall          !< 1/dz_wall
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  dz_wall_stag      !< wall grid spacing (edge-edge)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ddz_wall_stag     !< 1/dz_wall_stag
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  tt_wall_m         !< t_wall prognostic array 
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  zw                !< wall layer depths (m)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rho_c_window      !< volumetric heat capacity of the window material ( J m-3 K-1 ) (= 2.19E6)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  dz_window         !< window grid spacing (center-center)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ddz_window        !< 1/dz_window
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  dz_window_stag    !< window grid spacing (edge-edge)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ddz_window_stag   !< 1/dz_window_stag
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  tt_window_m       !< t_window prognostic array 
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  zw_window         !< window layer depths (m)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rho_c_green       !< volumetric heat capacity of the green material ( J m-3 K-1 ) (= 2.19E6)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rho_c_total_green !< volumetric heat capacity of the moist green material ( J m-3 K-1 ) (= 2.19E6)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  dz_green          !< green grid spacing (center-center)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ddz_green         !< 1/dz_green
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  dz_green_stag     !< green grid spacing (edge-edge)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ddz_green_stag    !< 1/dz_green_stag
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  tt_green_m        !< t_green prognostic array 
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  zw_green          !< green layer depths (m)

       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  gamma_w_green_sat !< hydraulic conductivity
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  lambda_w_green
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  gamma_w_green     !< hydraulic conductivity
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  tswc_h_m


!-- arrays for time averages
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  rad_net_av       !< average of rad_net_l
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfinsw_av      !< average of sw radiation falling to local surface including radiation from reflections
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfinlw_av      !< average of lw radiation falling to local surface including radiation from reflections
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfinswdir_av   !< average of direct sw radiation falling to local surface
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfinswdif_av   !< average of diffuse sw radiation from sky and model boundary falling to local surface
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfinlwdif_av   !< average of diffuse lw radiation from sky and model boundary falling to local surface
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfinswref_av   !< average of sw radiation falling to surface from reflections
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfinlwref_av   !< average of lw radiation falling to surface from reflections
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfoutsw_av     !< average of total sw radiation outgoing from nonvirtual surfaces surfaces after all reflection
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfoutlw_av     !< average of total lw radiation outgoing from nonvirtual surfaces surfaces after all reflection
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfins_av       !< average of array of residua of sw radiation absorbed in surface after last reflection
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfinl_av       !< average of array of residua of lw radiation absorbed in surface after last reflection
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfhf_av        !< average of total radiation flux incoming to minus outgoing from local surface  
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  wghf_eb_av       !< average of wghf_eb
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  wghf_eb_window_av  !< average of wghf_eb window
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  wghf_eb_green_av   !< average of wghf_eb window
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  iwghf_eb_av        !< indoor average of wghf_eb
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  iwghf_eb_window_av !< indoor average of wghf_eb window
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  wshf_eb_av       !< average of wshf_eb
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qsws_av           !< average of qsws
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qsws_veg_av       !< average of qsws_veg_eb
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  qsws_liq_av       !< average of qsws_liq_eb
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  t_surf_wall_av        !< average of wall surface temperature (K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  t_surf_av        !< average of wall surface temperature (K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  t_surf_window_av !< average of window surface temperature (K)
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  t_surf_green_av  !< average of green wall surface temperature (K)

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  pt_10cm_av       !< average of theta_10cm (K)
       
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  t_wall_av      !< Average of t_wall
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  t_window_av    !< Average of t_window
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  t_green_av     !< Average of t_green
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  swc_av         !< Average of swc

    END TYPE surf_type

    TYPE (bc_type), DIMENSION(0:1)           ::  bc_h        !< boundary condition data type, horizontal upward- and downward facing surfaces
    TYPE (bc_type), DIMENSION(0:3)           ::  bc_v        !< boundary condition data type, vertical surfaces

    TYPE (surf_type), DIMENSION(0:2), TARGET ::  surf_def_h  !< horizontal default surfaces (Up, Down, and Top)
    TYPE (surf_type), DIMENSION(0:3), TARGET ::  surf_def_v  !< vertical default surfaces (North, South, East, West)
    TYPE (surf_type)                , TARGET ::  surf_lsm_h  !< horizontal natural land surfaces, so far only upward-facing
    TYPE (surf_type), DIMENSION(0:3), TARGET ::  surf_lsm_v  !< vertical land surfaces (North, South, East, West)
    TYPE (surf_type)                , TARGET ::  surf_usm_h  !< horizontal urban surfaces, so far only upward-facing
    TYPE (surf_type), DIMENSION(0:3), TARGET ::  surf_usm_v  !< vertical urban surfaces (North, South, East, West)

    INTEGER(iwp), PARAMETER ::  ind_veg_wall  = 0            !< index for vegetation / wall-surface fraction, used for access of albedo, emissivity, etc., for each surface type   
    INTEGER(iwp), PARAMETER ::  ind_pav_green = 1            !< index for pavement / green-wall surface fraction, used for access of albedo, emissivity, etc., for each surface type
    INTEGER(iwp), PARAMETER ::  ind_wat_win   = 2            !< index for water / window-surface fraction, used for access of albedo, emissivity, etc., for each surface type

    INTEGER(iwp) ::  ns_h_on_file(0:2)                       !< total number of horizontal surfaces with the same facing, required for writing restart data 
    INTEGER(iwp) ::  ns_v_on_file(0:3)                       !< total number of vertical surfaces with the same facing, required for writing restart data 

    LOGICAL ::  vertical_surfaces_exist = .FALSE.   !< flag indicating that there are vertical urban/land surfaces 
                                                    !< in the domain (required to activiate RTM)

    LOGICAL ::  surf_bulk_cloud_model = .FALSE.        !< use cloud microphysics
    LOGICAL ::  surf_microphysics_morrison = .FALSE.   !< use 2-moment Morrison (add. prog. eq. for nc and qc)
    LOGICAL ::  surf_microphysics_seifert = .FALSE.    !< use 2-moment Seifert and Beheng scheme


    SAVE

    PRIVATE
    
    INTERFACE init_bc
       MODULE PROCEDURE init_bc
    END INTERFACE init_bc

    INTERFACE init_single_surface_properties
       MODULE PROCEDURE init_single_surface_properties
    END INTERFACE init_single_surface_properties
    
    INTERFACE init_surfaces
       MODULE PROCEDURE init_surfaces
    END INTERFACE init_surfaces

    INTERFACE init_surface_arrays
       MODULE PROCEDURE init_surface_arrays
    END INTERFACE init_surface_arrays

    INTERFACE surface_rrd_local
       MODULE PROCEDURE surface_rrd_local
    END INTERFACE surface_rrd_local

    INTERFACE surface_wrd_local
       MODULE PROCEDURE surface_wrd_local
    END INTERFACE surface_wrd_local

    INTERFACE surface_last_actions
       MODULE PROCEDURE surface_last_actions
    END INTERFACE surface_last_actions

    INTERFACE surface_restore_elements
       MODULE PROCEDURE surface_restore_elements_1d
       MODULE PROCEDURE surface_restore_elements_2d
    END INTERFACE surface_restore_elements

#if defined( _OPENACC )
    INTERFACE enter_surface_arrays
       MODULE PROCEDURE enter_surface_arrays
    END INTERFACE

    INTERFACE exit_surface_arrays
       MODULE PROCEDURE exit_surface_arrays
    END INTERFACE
#endif

!
!-- Public variables
    PUBLIC bc_h, bc_v, ind_pav_green, ind_veg_wall, ind_wat_win, ns_h_on_file, ns_v_on_file,       &
           surf_def_h, surf_def_v, surf_lsm_h, surf_lsm_v, surf_usm_h, surf_usm_v, surf_type,      &
           vertical_surfaces_exist, surf_bulk_cloud_model, surf_microphysics_morrison,             &
           surf_microphysics_seifert
!
!-- Public subroutines and functions
    PUBLIC init_bc,                                                                                &
           init_single_surface_properties,                                                         &
           init_surfaces,                                                                          &
           init_surface_arrays,                                                                    &
           surface_last_actions,                                                                   &
           surface_rrd_local,                                                                      &
           surface_restore_elements,                                                               &
           surface_wrd_local

#if defined( _OPENACC )
    PUBLIC enter_surface_arrays,                                                                   &
           exit_surface_arrays
#endif

 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize data type for setting boundary conditions at horizontal and  
!> vertical surfaces. 
!------------------------------------------------------------------------------!
    SUBROUTINE init_bc

       IMPLICIT NONE

       INTEGER(iwp) ::  i         !< loop index along x-direction
       INTEGER(iwp) ::  j         !< loop index along y-direction
       INTEGER(iwp) ::  k         !< loop index along y-direction
       INTEGER(iwp) ::  l         !< running index for differently aligned surfaces

       INTEGER(iwp), DIMENSION(0:1) ::  num_h         !< number of horizontal surfaces on subdomain
       INTEGER(iwp), DIMENSION(0:1) ::  num_h_kji     !< number of horizontal surfaces at (j,i)-grid point
       INTEGER(iwp), DIMENSION(0:1) ::  start_index_h !< local start index of horizontal surface elements
       
       INTEGER(iwp), DIMENSION(0:3) ::  num_v         !< number of vertical surfaces on subdomain
       INTEGER(iwp), DIMENSION(0:3) ::  num_v_kji     !< number of vertical surfaces at (j,i)-grid point
       INTEGER(iwp), DIMENSION(0:3) ::  start_index_v !< local start index of vertical surface elements
!
!--    Set offset indices, i.e. index difference between surface element and 
!--    surface-bounded grid point.
!--    Horizontal surfaces - no horizontal offsets
       bc_h(:)%ioff = 0
       bc_h(:)%joff = 0
!
!--    Horizontal surfaces, upward facing (0) and downward facing (1)
       bc_h(0)%koff   = -1
       bc_h(1)%koff   = 1
!
!--    Vertical surfaces - no vertical offset
       bc_v(0:3)%koff = 0
!
!--    North- and southward facing - no offset in x
       bc_v(0:1)%ioff = 0
!
!--    Northward facing offset in y
       bc_v(0)%joff = -1
!
!--    Southward facing offset in y
       bc_v(1)%joff = 1
!
!--    East- and westward facing - no offset in y
       bc_v(2:3)%joff = 0
!
!--    Eastward facing offset in x
       bc_v(2)%ioff = -1
!
!--    Westward facing offset in y
       bc_v(3)%ioff = 1
!
!--    Initialize data structure for horizontal surfaces, i.e. count the number
!--    of surface elements, allocate and initialize the respective index arrays, 
!--    and set the respective start and end indices at each (j,i)-location. 
!--    The index space is defined also over the ghost points, so that e.g. 
!--    boundary conditions for diagnostic quanitities can be set on ghost  
!--    points so that no exchange is required any more. 
       DO  l = 0, 1
!
!--       Count the number of upward- and downward-facing surfaces on subdomain
          num_h(l) = 0
          DO  i = nxlg, nxrg
             DO  j = nysg, nyng
                DO  k = nzb+1, nzt
!         
!--                Check if current gridpoint belongs to the atmosphere
                   IF ( BTEST( wall_flags_total_0(k,j,i), 0 ) )  THEN
                      IF ( .NOT. BTEST( wall_flags_total_0(k+bc_h(l)%koff,     &
                                                     j+bc_h(l)%joff,           &
                                                     i+bc_h(l)%ioff), 0 ) )    &
                         num_h(l) = num_h(l) + 1
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
!         
!--       Save the number of horizontal surface elements
          bc_h(l)%ns = num_h(l)
!
!--       ALLOCATE arrays for horizontal surfaces
          ALLOCATE( bc_h(l)%i(1:bc_h(l)%ns) )
          ALLOCATE( bc_h(l)%j(1:bc_h(l)%ns) )
          ALLOCATE( bc_h(l)%k(1:bc_h(l)%ns) )
          ALLOCATE( bc_h(l)%start_index(nysg:nyng,nxlg:nxrg) )
          ALLOCATE( bc_h(l)%end_index(nysg:nyng,nxlg:nxrg)   )
          bc_h(l)%start_index = 1
          bc_h(l)%end_index   = 0
          
          num_h(l)         = 1
          start_index_h(l) = 1
          DO  i = nxlg, nxrg
             DO  j = nysg, nyng
          
                num_h_kji(l) = 0
                DO  k = nzb+1, nzt
!         
!--                Check if current gridpoint belongs to the atmosphere
                   IF ( BTEST( wall_flags_total_0(k,j,i), 0 ) )  THEN
!         
!--                   Upward-facing
                      IF ( .NOT. BTEST( wall_flags_total_0(k+bc_h(l)%koff,     &
                                                     j+bc_h(l)%joff,           &
                                                     i+bc_h(l)%ioff), 0 )      &
                         )  THEN
                         bc_h(l)%i(num_h(l)) = i
                         bc_h(l)%j(num_h(l)) = j
                         bc_h(l)%k(num_h(l)) = k
                         num_h_kji(l)        = num_h_kji(l) + 1
                         num_h(l)            = num_h(l) + 1
                      ENDIF
                   ENDIF
                ENDDO
                bc_h(l)%start_index(j,i) = start_index_h(l)
                bc_h(l)%end_index(j,i)   = bc_h(l)%start_index(j,i) +          &
                                           num_h_kji(l) - 1
                start_index_h(l)         = bc_h(l)%end_index(j,i) + 1
             ENDDO
          ENDDO
       ENDDO

!
!--    Initialize data structure for vertical surfaces, i.e. count the number
!--    of surface elements, allocate and initialize the respective index arrays, 
!--    and set the respective start and end indices at each (j,i)-location.
       DO  l = 0, 3
!
!--       Count the number of upward- and downward-facing surfaces on subdomain
          num_v(l) = 0
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
!         
!--                Check if current gridpoint belongs to the atmosphere
                   IF ( BTEST( wall_flags_total_0(k,j,i), 0 ) )  THEN
                      IF ( .NOT. BTEST( wall_flags_total_0(k+bc_v(l)%koff,     &
                                                     j+bc_v(l)%joff,           &
                                                     i+bc_v(l)%ioff), 0 ) )    &
                         num_v(l) = num_v(l) + 1
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
!         
!--       Save the number of horizontal surface elements
          bc_v(l)%ns = num_v(l)
!
!--       ALLOCATE arrays for horizontal surfaces. In contrast to the
!--       horizontal surfaces, the index space is not defined over the 
!--       ghost points. 
          ALLOCATE( bc_v(l)%i(1:bc_v(l)%ns) )
          ALLOCATE( bc_v(l)%j(1:bc_v(l)%ns) )
          ALLOCATE( bc_v(l)%k(1:bc_v(l)%ns) )
          ALLOCATE( bc_v(l)%start_index(nys:nyn,nxl:nxr) )
          ALLOCATE( bc_v(l)%end_index(nys:nyn,nxl:nxr)   )
          bc_v(l)%start_index = 1
          bc_v(l)%end_index   = 0
          
          num_v(l)         = 1
          start_index_v(l) = 1
          DO  i = nxl, nxr
             DO  j = nys, nyn
          
                num_v_kji(l) = 0
                DO  k = nzb+1, nzt
!         
!--                Check if current gridpoint belongs to the atmosphere
                   IF ( BTEST( wall_flags_total_0(k,j,i), 0 ) )  THEN
!         
!--                   Upward-facing
                      IF ( .NOT. BTEST( wall_flags_total_0(k+bc_v(l)%koff,     &
                                                     j+bc_v(l)%joff,           &
                                                     i+bc_v(l)%ioff), 0 )      &
                         )  THEN
                         bc_v(l)%i(num_v(l)) = i
                         bc_v(l)%j(num_v(l)) = j
                         bc_v(l)%k(num_v(l)) = k
                         num_v_kji(l)        = num_v_kji(l) + 1
                         num_v(l)            = num_v(l) + 1
                      ENDIF
                   ENDIF
                ENDDO
                bc_v(l)%start_index(j,i) = start_index_v(l)
                bc_v(l)%end_index(j,i)   = bc_v(l)%start_index(j,i) +          &
                                           num_v_kji(l) - 1
                start_index_v(l)         = bc_v(l)%end_index(j,i) + 1
             ENDDO
          ENDDO
       ENDDO


    END SUBROUTINE init_bc


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize horizontal and vertical surfaces. Counts the number of default-,
!> natural and urban surfaces and allocates memory, respectively. 
!------------------------------------------------------------------------------!
    SUBROUTINE init_surface_arrays


       USE pegrid


       IMPLICIT NONE

       INTEGER(iwp)                 ::  i         !< running index x-direction 
       INTEGER(iwp)                 ::  j         !< running index y-direction
       INTEGER(iwp)                 ::  k         !< running index z-direction
       INTEGER(iwp)                 ::  l         !< index variable for surface facing
       INTEGER(iwp)                 ::  num_lsm_h !< number of horizontally-aligned natural surfaces 
       INTEGER(iwp)                 ::  num_usm_h !< number of horizontally-aligned urban surfaces 

       INTEGER(iwp), DIMENSION(0:2) ::  num_def_h !< number of horizontally-aligned default surfaces 
       INTEGER(iwp), DIMENSION(0:3) ::  num_def_v !< number of vertically-aligned default surfaces 
       INTEGER(iwp), DIMENSION(0:3) ::  num_lsm_v !< number of vertically-aligned natural surfaces 
       INTEGER(iwp), DIMENSION(0:3) ::  num_usm_v !< number of vertically-aligned urban surfaces 

       INTEGER(iwp)              ::  num_surf_v_l !< number of vertically-aligned local urban/land surfaces
       INTEGER(iwp)              ::  num_surf_v   !< number of vertically-aligned total urban/land surfaces

       LOGICAL ::  building                       !< flag indicating building grid point
       LOGICAL ::  terrain                        !< flag indicating natural terrain grid point
       LOGICAL ::  unresolved_building            !< flag indicating a grid point where actually a building is 
                                                  !< defined but not resolved by the vertical grid 

       num_def_h = 0
       num_def_v = 0
       num_lsm_h = 0
       num_lsm_v = 0
       num_usm_h = 0
       num_usm_v = 0
!
!--    Surfaces are classified according to the input data read from static
!--    input file. If no input file is present, all surfaces are classified
!--    either as natural, urban, or default, depending on the setting of 
!--    land_surface and urban_surface. To control this, use the control
!--    flag topo_no_distinct
!
!--    Count number of horizontal surfaces on local domain 
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Check if current gridpoint belongs to the atmosphere
                IF ( BTEST( wall_flags_total_0(k,j,i), 0 ) )  THEN
!
!--                Check if grid point adjoins to any upward-facing horizontal
!--                surface, e.g. the Earth surface, plane roofs, or ceilings.

                   IF ( .NOT. BTEST( wall_flags_total_0(k-1,j,i), 0 ) )  THEN
!
!--                   Determine flags indicating a terrain surface, a building
!--                   surface, 
                      terrain  = BTEST( wall_flags_total_0(k-1,j,i), 5 )  .OR.       &
                                 topo_no_distinct
                      building = BTEST( wall_flags_total_0(k-1,j,i), 6 )  .OR.       &
                                 topo_no_distinct
!
!--                   unresolved_building indicates a surface with equal height 
!--                   as terrain but with a non-grid resolved building on top.
!--                   These surfaces will be flagged as urban surfaces.
                      unresolved_building = BTEST( wall_flags_total_0(k-1,j,i), 5 )  &
                                     .AND.  BTEST( wall_flags_total_0(k-1,j,i), 6 )
!
!--                   Land-surface type
                      IF ( land_surface  .AND.  terrain  .AND.                 &
                           .NOT. unresolved_building )  THEN
                         num_lsm_h    = num_lsm_h    + 1 
!
!--                   Urban surface tpye
                      ELSEIF ( urban_surface  .AND.  building )  THEN
                         num_usm_h    = num_usm_h    + 1 
!
!--                   Default-surface type
                      ELSEIF ( .NOT. land_surface    .AND.                     &
                               .NOT. urban_surface )  THEN
                               
                         num_def_h(0) = num_def_h(0) + 1
!
!--                   Unclassifified surface-grid point. Give error message.
                      ELSE 
                         WRITE( message_string, * )                           &
                                          'Unclassified upward-facing ' //    &
                                          'surface element at '//             &
                                          'grid point (k,j,i) = ', k, j, i
                         CALL message( 'surface_mod', 'PA0698', 1, 2, myid, 6, 0 )
                      ENDIF

                   ENDIF
!
!--                Check for top-fluxes
                   IF ( k == nzt  .AND.  use_top_fluxes )  THEN
                      num_def_h(2) = num_def_h(2) + 1
!
!--                Check for any other downward-facing surface. So far only for 
!--                default surface type.
                   ELSEIF ( .NOT. BTEST( wall_flags_total_0(k+1,j,i), 0 ) )  THEN
                      num_def_h(1) = num_def_h(1) + 1
                   ENDIF 

                ENDIF
             ENDDO
          ENDDO
       ENDDO
!
!--    Count number of vertical surfaces on local domain 
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                IF ( BTEST( wall_flags_total_0(k,j,i), 0 ) )  THEN
!
!--                Northward-facing
                   IF ( .NOT. BTEST( wall_flags_total_0(k,j-1,i), 0 ) )  THEN
!
!--                   Determine flags indicating terrain or building

                      terrain  = BTEST( wall_flags_total_0(k,j-1,i), 5 )  .OR.       &
                                 topo_no_distinct
                      building = BTEST( wall_flags_total_0(k,j-1,i), 6 )   .OR.      &
                                 topo_no_distinct

                      unresolved_building = BTEST( wall_flags_total_0(k,j-1,i), 5 )  &
                                     .AND.  BTEST( wall_flags_total_0(k,j-1,i), 6 )
                                     
                      IF (  land_surface  .AND.  terrain  .AND.                &
                           .NOT. unresolved_building )  THEN
                         num_lsm_v(0) = num_lsm_v(0) + 1 
                      ELSEIF ( urban_surface  .AND.  building )  THEN
                         num_usm_v(0) = num_usm_v(0) + 1 
!
!--                   Default-surface type
                      ELSEIF ( .NOT. land_surface    .AND.                     &
                               .NOT. urban_surface )  THEN
                         num_def_v(0) = num_def_v(0) + 1 
!
!--                   Unclassifified surface-grid point. Give error message.
                      ELSE 
                         WRITE( message_string, * )                            &
                                          'Unclassified northward-facing ' //  &
                                          'surface element at '//              &
                                          'grid point (k,j,i) = ', k, j, i
                         CALL message( 'surface_mod', 'PA0698', 1, 2, myid, 6, 0 )

                      ENDIF
                   ENDIF
!
!--                Southward-facing
                   IF ( .NOT. BTEST( wall_flags_total_0(k,j+1,i), 0 ) )  THEN
!
!--                   Determine flags indicating terrain or building
                      terrain  = BTEST( wall_flags_total_0(k,j+1,i), 5 )  .OR.       &
                                 topo_no_distinct
                      building = BTEST( wall_flags_total_0(k,j+1,i), 6 )  .OR.       &
                                 topo_no_distinct
                                 
                      unresolved_building = BTEST( wall_flags_total_0(k,j+1,i), 5 )  &
                                     .AND.  BTEST( wall_flags_total_0(k,j+1,i), 6 ) 
                                
                      IF (  land_surface  .AND.  terrain  .AND.                &
                           .NOT. unresolved_building )  THEN
                         num_lsm_v(1) = num_lsm_v(1) + 1 
                      ELSEIF ( urban_surface  .AND.  building )  THEN
                         num_usm_v(1) = num_usm_v(1) + 1 
!
!--                   Default-surface type
                      ELSEIF ( .NOT. land_surface    .AND.                     &
                               .NOT. urban_surface )  THEN
                         num_def_v(1) = num_def_v(1) + 1 
!
!--                   Unclassifified surface-grid point. Give error message.
                      ELSE 
                         WRITE( message_string, * )                            &
                                          'Unclassified southward-facing ' //  &
                                          'surface element at '//              &
                                          'grid point (k,j,i) = ', k, j, i
                         CALL message( 'surface_mod', 'PA0698', 1, 2, myid, 6, 0 )

                      ENDIF
                   ENDIF
!
!--                Eastward-facing
                   IF ( .NOT. BTEST( wall_flags_total_0(k,j,i-1), 0 ) )  THEN
!
!--                   Determine flags indicating terrain or building
                      terrain  = BTEST( wall_flags_total_0(k,j,i-1), 5 )  .OR.       &
                                 topo_no_distinct
                      building = BTEST( wall_flags_total_0(k,j,i-1), 6 )  .OR.       &
                                 topo_no_distinct
                                 
                      unresolved_building = BTEST( wall_flags_total_0(k,j,i-1), 5 )  &
                                     .AND.  BTEST( wall_flags_total_0(k,j,i-1), 6 )
                                     
                      IF (  land_surface  .AND.  terrain  .AND.                &
                           .NOT. unresolved_building )  THEN
                         num_lsm_v(2) = num_lsm_v(2) + 1 
                      ELSEIF ( urban_surface  .AND.  building )  THEN
                         num_usm_v(2) = num_usm_v(2) + 1 
!
!--                   Default-surface type
                      ELSEIF ( .NOT. land_surface    .AND.                     &
                               .NOT. urban_surface )  THEN
                         num_def_v(2) = num_def_v(2) + 1 
!
!--                   Unclassifified surface-grid point. Give error message.
                      ELSE 
                         WRITE( message_string, * )                            &
                                          'Unclassified eastward-facing ' //   &
                                          'surface element at '//              &
                                          'grid point (k,j,i) = ', k, j, i
                         CALL message( 'surface_mod', 'PA0698', 1, 2, myid, 6, 0 )

                      ENDIF
                   ENDIF
!
!--                Westward-facing
                   IF ( .NOT. BTEST( wall_flags_total_0(k,j,i+1), 0 ) )  THEN
!
!--                   Determine flags indicating terrain or building
                      terrain  = BTEST( wall_flags_total_0(k,j,i+1), 5 )  .OR.       &
                                 topo_no_distinct
                      building = BTEST( wall_flags_total_0(k,j,i+1), 6 )  .OR.       &
                                 topo_no_distinct
                                 
                      unresolved_building = BTEST( wall_flags_total_0(k,j,i+1), 5 )  &
                                     .AND.  BTEST( wall_flags_total_0(k,j,i+1), 6 )
                                 
                      IF (  land_surface  .AND.  terrain  .AND.                &
                           .NOT. unresolved_building )  THEN
                         num_lsm_v(3) = num_lsm_v(3) + 1 
                      ELSEIF ( urban_surface  .AND.  building )  THEN
                         num_usm_v(3) = num_usm_v(3) + 1 
!
!--                   Default-surface type
                      ELSEIF ( .NOT. land_surface    .AND.                     &
                               .NOT. urban_surface )  THEN
                         num_def_v(3) = num_def_v(3) + 1 
!
!--                   Unclassifified surface-grid point. Give error message.
                      ELSE 
                         WRITE( message_string, * )                            &
                                          'Unclassified westward-facing ' //   &
                                          'surface element at '//              &
                                          'grid point (k,j,i) = ', k, j, i
                         CALL message( 'surface_mod', 'PA0698', 1, 2, myid, 6, 0 )

                      ENDIF
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO

!
!--    Store number of surfaces per core.
!--    Horizontal surface, default type, upward facing
       surf_def_h(0)%ns = num_def_h(0)
!
!--    Horizontal surface, default type, downward facing
       surf_def_h(1)%ns = num_def_h(1)
!
!--    Horizontal surface, default type, top downward facing
       surf_def_h(2)%ns = num_def_h(2)
!
!--    Horizontal surface, natural type, so far only upward-facing
       surf_lsm_h%ns    = num_lsm_h  
!
!--    Horizontal surface, urban type, so far only upward-facing
       surf_usm_h%ns    = num_usm_h    
!
!--    Vertical surface, default type, northward facing
       surf_def_v(0)%ns = num_def_v(0)
!
!--    Vertical surface, default type, southward facing
       surf_def_v(1)%ns = num_def_v(1)
!
!--    Vertical surface, default type, eastward facing
       surf_def_v(2)%ns = num_def_v(2)
!
!--    Vertical surface, default type, westward facing
       surf_def_v(3)%ns = num_def_v(3)
!
!--    Vertical surface, natural type, northward facing
       surf_lsm_v(0)%ns = num_lsm_v(0)
!
!--    Vertical surface, natural type, southward facing
       surf_lsm_v(1)%ns = num_lsm_v(1)
!
!--    Vertical surface, natural type, eastward facing
       surf_lsm_v(2)%ns = num_lsm_v(2)
!
!--    Vertical surface, natural type, westward facing
       surf_lsm_v(3)%ns = num_lsm_v(3)
!
!--    Vertical surface, urban type, northward facing
       surf_usm_v(0)%ns = num_usm_v(0)
!
!--    Vertical surface, urban type, southward facing
       surf_usm_v(1)%ns = num_usm_v(1)
!
!--    Vertical surface, urban type, eastward facing
       surf_usm_v(2)%ns = num_usm_v(2)
!
!--    Vertical surface, urban type, westward facing
       surf_usm_v(3)%ns = num_usm_v(3)
!
!--    Allocate required attributes for horizontal surfaces - default type. 
!--    Upward-facing (l=0) and downward-facing (l=1).
       DO  l = 0, 1
          CALL allocate_surface_attributes_h ( surf_def_h(l), nys, nyn, nxl, nxr )
       ENDDO
!
!--    Allocate required attributes for model top
       CALL allocate_surface_attributes_h_top ( surf_def_h(2), nys, nyn, nxl, nxr )
!
!--    Allocate required attributes for horizontal surfaces - natural type. 
       CALL allocate_surface_attributes_h ( surf_lsm_h, nys, nyn, nxl, nxr )
!
!--    Allocate required attributes for horizontal surfaces - urban type. 
       CALL allocate_surface_attributes_h ( surf_usm_h, nys, nyn, nxl, nxr )

!
!--    Allocate required attributes for vertical surfaces. 
!--    Northward-facing (l=0), southward-facing (l=1), eastward-facing (l=2)
!--    and westward-facing (l=3).
!--    Default type.
       DO  l = 0, 3
          CALL allocate_surface_attributes_v ( surf_def_v(l),                  &
                                               nys, nyn, nxl, nxr )
       ENDDO
!
!--    Natural type
       DO  l = 0, 3
          CALL allocate_surface_attributes_v ( surf_lsm_v(l),                  &
                                               nys, nyn, nxl, nxr )
       ENDDO
!
!--    Urban type
       DO  l = 0, 3
          CALL allocate_surface_attributes_v ( surf_usm_v(l),                  &
                                               nys, nyn, nxl, nxr )
       ENDDO
!
!--    Set the flag for the existence of vertical urban/land surfaces
       num_surf_v_l = 0
       DO  l = 0, 3
          num_surf_v_l = num_surf_v_l + surf_usm_v(l)%ns + surf_lsm_v(l)%ns
       ENDDO

#if defined( __parallel )
       CALL MPI_ALLREDUCE( num_surf_v_l, num_surf_v, 1, MPI_INTEGER,           &
                           MPI_SUM, comm2d, ierr)
#else
       num_surf_v = num_surf_v_l
#endif
       IF ( num_surf_v > 0 )  vertical_surfaces_exist = .TRUE.
        

    END SUBROUTINE init_surface_arrays


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Enter horizontal and vertical surfaces.
!------------------------------------------------------------------------------!
#if defined( _OPENACC )
    SUBROUTINE enter_surface_arrays

       IMPLICIT NONE

       INTEGER(iwp) ::  l     !<
       
       !$ACC ENTER DATA &
       !$ACC COPYIN(surf_def_h(0:2)) &
       !$ACC COPYIN(surf_def_v(0:3)) &
       !$ACC COPYIN(surf_lsm_h) &
       !$ACC COPYIN(surf_lsm_v(0:3)) &
       !$ACC COPYIN(surf_usm_h) &
       !$ACC COPYIN(surf_usm_v(0:3))

       ! Copy data in surf_def_h(0:2)
       DO  l = 0, 1
          CALL enter_surface_attributes_h(surf_def_h(l))
       ENDDO
       CALL enter_surface_attributes_h_top(surf_def_h(2))
       ! Copy data in surf_def_v(0:3)
       DO  l = 0, 3
          CALL enter_surface_attributes_v(surf_def_v(l))
       ENDDO
       ! Copy data in surf_lsm_h
       CALL enter_surface_attributes_h(surf_lsm_h)
       ! Copy data in surf_lsm_v(0:3)
       DO  l = 0, 3
          CALL enter_surface_attributes_v(surf_lsm_v(l))
       ENDDO
       ! Copy data in surf_usm_h
       CALL enter_surface_attributes_h(surf_usm_h)
       ! Copy data in surf_usm_v(0:3)
       DO  l = 0, 3
          CALL enter_surface_attributes_v(surf_usm_v(l))
       ENDDO

    END SUBROUTINE enter_surface_arrays
#endif

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Exit horizontal and vertical surfaces.
!------------------------------------------------------------------------------!
#if defined( _OPENACC )
    SUBROUTINE exit_surface_arrays

       IMPLICIT NONE

       INTEGER(iwp) ::  l     !<
       
       ! Delete data in surf_def_h(0:2)
       DO  l = 0, 1
          CALL exit_surface_attributes_h(surf_def_h(l))
       ENDDO
       CALL exit_surface_attributes_h(surf_def_h(2))
       ! Delete data in surf_def_v(0:3)
       DO  l = 0, 3
          CALL exit_surface_attributes_v(surf_def_v(l))
       ENDDO
       ! Delete data in surf_lsm_h
       CALL exit_surface_attributes_h(surf_lsm_h)
       ! Delete data in surf_lsm_v(0:3)
       DO  l = 0, 3
          CALL exit_surface_attributes_v(surf_lsm_v(l))
       ENDDO
       ! Delete data in surf_usm_h
       CALL exit_surface_attributes_h(surf_usm_h)
       ! Delete data in surf_usm_v(0:3)
       DO  l = 0, 3
          CALL exit_surface_attributes_v(surf_usm_v(l))
       ENDDO

       !$ACC EXIT DATA &
       !$ACC DELETE(surf_def_h(0:2)) &
       !$ACC DELETE(surf_def_v(0:3)) &
       !$ACC DELETE(surf_lsm_h) &
       !$ACC DELETE(surf_lsm_v(0:3)) &
       !$ACC DELETE(surf_usm_h) &
       !$ACC DELETE(surf_usm_v(0:3))

    END SUBROUTINE exit_surface_arrays
#endif

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Deallocating memory for upward and downward-facing horizontal surface types, 
!> except for top fluxes. 
!------------------------------------------------------------------------------!
    SUBROUTINE deallocate_surface_attributes_h( surfaces )

       IMPLICIT NONE


       TYPE(surf_type) ::  surfaces  !< respective surface type


       DEALLOCATE ( surfaces%start_index )
       DEALLOCATE ( surfaces%end_index )
!
!--    Indices to locate surface element
       DEALLOCATE ( surfaces%i )
       DEALLOCATE ( surfaces%j )
       DEALLOCATE ( surfaces%k )
!
!--    Surface-layer height
       DEALLOCATE ( surfaces%z_mo )
!
!--    Surface orientation
       DEALLOCATE ( surfaces%facing )
!
!--    Surface-parallel wind velocity
       DEALLOCATE ( surfaces%uvw_abs )
!
!--    Roughness
       DEALLOCATE ( surfaces%z0 )
       DEALLOCATE ( surfaces%z0h )
       DEALLOCATE ( surfaces%z0q )
!
!--    Friction velocity
       DEALLOCATE ( surfaces%us )
!
!--    Stability parameter
       DEALLOCATE ( surfaces%ol )
!
!--    Bulk Richardson number
       DEALLOCATE ( surfaces%rib )
!
!--    Vertical momentum fluxes of u and v
       DEALLOCATE ( surfaces%usws )  
       DEALLOCATE ( surfaces%vsws )  
!
!--    Required in production_e
       IF ( .NOT. constant_diffusion )  THEN    
          DEALLOCATE ( surfaces%u_0 )  
          DEALLOCATE ( surfaces%v_0 )
       ENDIF 
!
!--    Characteristic temperature and surface flux of sensible heat
       DEALLOCATE ( surfaces%ts )    
       DEALLOCATE ( surfaces%shf )
!
!--    surface temperature
       DEALLOCATE ( surfaces%pt_surface ) 
!
!--    Characteristic humidity and surface flux of latent heat
       IF ( humidity )  THEN          
          DEALLOCATE ( surfaces%qs ) 
          DEALLOCATE ( surfaces%qsws )  
          DEALLOCATE ( surfaces%q_surface   )
          DEALLOCATE ( surfaces%vpt_surface )
       ENDIF 
!
!--    Characteristic scalar and surface flux of scalar
       IF ( passive_scalar )  THEN
          DEALLOCATE ( surfaces%ss )   
          DEALLOCATE ( surfaces%ssws ) 
       ENDIF 
!
!--    Scaling parameter (cs*) and surface flux of chemical species
       IF ( air_chemistry )  THEN
          DEALLOCATE ( surfaces%css )   
          DEALLOCATE ( surfaces%cssws ) 
       ENDIF 
!
!--    Arrays for storing potential temperature and
!--    mixing ratio at first grid level
       DEALLOCATE ( surfaces%pt1 )
       DEALLOCATE ( surfaces%qv1 )
       DEALLOCATE ( surfaces%vpt1 )
       
!
!--       
       IF ( surf_bulk_cloud_model .AND. surf_microphysics_morrison)  THEN
          DEALLOCATE ( surfaces%qcs )
          DEALLOCATE ( surfaces%ncs )
          DEALLOCATE ( surfaces%qcsws )
          DEALLOCATE ( surfaces%ncsws )
       ENDIF
!
!--       
       IF ( surf_bulk_cloud_model .AND. surf_microphysics_seifert)  THEN
          DEALLOCATE ( surfaces%qrs )
          DEALLOCATE ( surfaces%nrs )
          DEALLOCATE ( surfaces%qrsws )
          DEALLOCATE ( surfaces%nrsws )
       ENDIF
!
!--    Salinity surface flux
       IF ( ocean_mode )  DEALLOCATE ( surfaces%sasws )

    END SUBROUTINE deallocate_surface_attributes_h


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocating memory for upward and downward-facing horizontal surface types, 
!> except for top fluxes. 
!------------------------------------------------------------------------------!
    SUBROUTINE allocate_surface_attributes_h( surfaces,                        &
                                              nys_l, nyn_l, nxl_l, nxr_l )

       IMPLICIT NONE

       INTEGER(iwp) ::  nyn_l  !< north bound of local 2d array start/end_index, is equal to nyn, except for restart-array
       INTEGER(iwp) ::  nys_l  !< south bound of local 2d array start/end_index, is equal to nyn, except for restart-array
       INTEGER(iwp) ::  nxl_l  !< west bound of local 2d array start/end_index, is equal to nyn, except for restart-array
       INTEGER(iwp) ::  nxr_l  !< east bound of local 2d array start/end_index, is equal to nyn, except for restart-array

       TYPE(surf_type) ::  surfaces  !< respective surface type

!
!--    Allocate arrays for start and end index of horizontal surface type 
!--    for each (j,i)-grid point. This is required e.g. in diffion_x, which is
!--    called for each (j,i). In order to find the location where the 
!--    respective flux is store within the surface-type, start- and end- 
!--    index are stored for each (j,i). For example, each (j,i) can have
!--    several entries where fluxes for horizontal surfaces might be stored,
!--    e.g. for overhanging structures where several upward-facing surfaces
!--    might exist for given (j,i).
!--    If no surface of respective type exist at current (j,i), set indicies
!--    such that loop in diffusion routines will not be entered.
       ALLOCATE ( surfaces%start_index(nys_l:nyn_l,nxl_l:nxr_l) )
       ALLOCATE ( surfaces%end_index(nys_l:nyn_l,nxl_l:nxr_l)   )
       surfaces%start_index = 0
       surfaces%end_index   = -1
!
!--    Indices to locate surface element
       ALLOCATE ( surfaces%i(1:surfaces%ns)  )
       ALLOCATE ( surfaces%j(1:surfaces%ns)  )
       ALLOCATE ( surfaces%k(1:surfaces%ns)  )
!
!--    Surface-layer height
       ALLOCATE ( surfaces%z_mo(1:surfaces%ns) )
!
!--    Surface orientation
       ALLOCATE ( surfaces%facing(1:surfaces%ns) )
!
!--    Surface-parallel wind velocity
       ALLOCATE ( surfaces%uvw_abs(1:surfaces%ns) )
!
!--    Roughness
       ALLOCATE ( surfaces%z0(1:surfaces%ns)  )
       ALLOCATE ( surfaces%z0h(1:surfaces%ns) )
       ALLOCATE ( surfaces%z0q(1:surfaces%ns) )
!
!--    Friction velocity
       ALLOCATE ( surfaces%us(1:surfaces%ns) )
!
!--    Stability parameter
       ALLOCATE ( surfaces%ol(1:surfaces%ns) )
!
!--    Bulk Richardson number
       ALLOCATE ( surfaces%rib(1:surfaces%ns) )
!
!--    Vertical momentum fluxes of u and v
       ALLOCATE ( surfaces%usws(1:surfaces%ns) )  
       ALLOCATE ( surfaces%vsws(1:surfaces%ns) )  
!
!--    Required in production_e
       IF ( .NOT. constant_diffusion )  THEN    
          ALLOCATE ( surfaces%u_0(1:surfaces%ns) )  
          ALLOCATE ( surfaces%v_0(1:surfaces%ns) )
       ENDIF 
!
!--    Characteristic temperature and surface flux of sensible heat
       ALLOCATE ( surfaces%ts(1:surfaces%ns)  )    
       ALLOCATE ( surfaces%shf(1:surfaces%ns) )
!
!--    Surface temperature
       ALLOCATE ( surfaces%pt_surface(1:surfaces%ns) ) 
!
!--    Characteristic humidity, surface flux of latent heat, and surface virtual potential temperature
       IF ( humidity )  THEN
          ALLOCATE ( surfaces%qs(1:surfaces%ns)   ) 
          ALLOCATE ( surfaces%qsws(1:surfaces%ns) )      
          ALLOCATE ( surfaces%q_surface(1:surfaces%ns)   ) 
          ALLOCATE ( surfaces%vpt_surface(1:surfaces%ns) ) 
       ENDIF 

!
!--    Characteristic scalar and surface flux of scalar
       IF ( passive_scalar )  THEN
          ALLOCATE ( surfaces%ss(1:surfaces%ns)   )   
          ALLOCATE ( surfaces%ssws(1:surfaces%ns) ) 
       ENDIF 
!
!--    Scaling parameter (cs*) and surface flux of chemical species
       IF ( air_chemistry )  THEN
          ALLOCATE ( surfaces%css(1:nvar,1:surfaces%ns)   )   
          ALLOCATE ( surfaces%cssws(1:nvar,1:surfaces%ns) ) 
       ENDIF 
!
!--    Arrays for storing potential temperature and
!--    mixing ratio at first grid level
       ALLOCATE ( surfaces%pt1(1:surfaces%ns) )
       ALLOCATE ( surfaces%qv1(1:surfaces%ns) )
       ALLOCATE ( surfaces%vpt1(1:surfaces%ns) )
!
!--       
       IF ( surf_bulk_cloud_model .AND. surf_microphysics_morrison)  THEN
          ALLOCATE ( surfaces%qcs(1:surfaces%ns)   )
          ALLOCATE ( surfaces%ncs(1:surfaces%ns)   )
          ALLOCATE ( surfaces%qcsws(1:surfaces%ns) )
          ALLOCATE ( surfaces%ncsws(1:surfaces%ns) )
       ENDIF
!
!--       
       IF ( surf_bulk_cloud_model .AND. surf_microphysics_seifert)  THEN
          ALLOCATE ( surfaces%qrs(1:surfaces%ns)   )
          ALLOCATE ( surfaces%nrs(1:surfaces%ns)   )
          ALLOCATE ( surfaces%qrsws(1:surfaces%ns) )
          ALLOCATE ( surfaces%nrsws(1:surfaces%ns) )
       ENDIF
!
!--    Salinity surface flux
       IF ( ocean_mode )  ALLOCATE ( surfaces%sasws(1:surfaces%ns) )

    END SUBROUTINE allocate_surface_attributes_h


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Exit memory for upward and downward-facing horizontal surface types,
!> except for top fluxes.
!------------------------------------------------------------------------------!
#if defined( _OPENACC )
    SUBROUTINE exit_surface_attributes_h( surfaces )

       IMPLICIT NONE
   
       TYPE(surf_type) ::  surfaces  !< respective surface type
   
       !$ACC EXIT DATA &
       !$ACC DELETE(surfaces%start_index(nys:nyn,nxl:nxr)) &
       !$ACC DELETE(surfaces%end_index(nys:nyn,nxl:nxr)) &
       !$ACC DELETE(surfaces%i(1:surfaces%ns)) &
       !$ACC DELETE(surfaces%j(1:surfaces%ns)) &
       !$ACC DELETE(surfaces%k(1:surfaces%ns)) &
       !$ACC DELETE(surfaces%z_mo(1:surfaces%ns)) &
       !$ACC DELETE(surfaces%uvw_abs(1:surfaces%ns)) &
       !$ACC DELETE(surfaces%z0(1:surfaces%ns)) &
       !$ACC COPYOUT(surfaces%us(1:surfaces%ns)) &
       !$ACC COPYOUT(surfaces%ol(1:surfaces%ns)) &
       !$ACC DELETE(surfaces%rib(1:surfaces%ns)) &
       !$ACC COPYOUT(surfaces%usws(1:surfaces%ns)) &
       !$ACC COPYOUT(surfaces%vsws(1:surfaces%ns)) &
       !$ACC COPYOUT(surfaces%ts(1:surfaces%ns)) &
       !$ACC COPYOUT(surfaces%shf(1:surfaces%ns)) &
       !$ACC DELETE(surfaces%pt_surface(1:surfaces%ns)) &
       !$ACC DELETE(surfaces%pt1(1:surfaces%ns)) &
       !$ACC DELETE(surfaces%qv1(1:surfaces%ns))

       IF ( .NOT. constant_diffusion )  THEN
          !$ACC EXIT DATA &
          !$ACC DELETE(surfaces%u_0(1:surfaces%ns)) &
          !$ACC DELETE(surfaces%v_0(1:surfaces%ns))
       ENDIF
   
    END SUBROUTINE exit_surface_attributes_h
#endif

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Enter memory for upward and downward-facing horizontal surface types,
!> except for top fluxes.
!------------------------------------------------------------------------------!
#if defined( _OPENACC )
    SUBROUTINE enter_surface_attributes_h( surfaces )

       IMPLICIT NONE

       TYPE(surf_type) ::  surfaces  !< respective surface type

       !$ACC ENTER DATA &
       !$ACC COPYIN(surfaces%start_index(nys:nyn,nxl:nxr)) &
       !$ACC COPYIN(surfaces%end_index(nys:nyn,nxl:nxr)) &
       !$ACC COPYIN(surfaces%i(1:surfaces%ns)) &
       !$ACC COPYIN(surfaces%j(1:surfaces%ns)) &
       !$ACC COPYIN(surfaces%k(1:surfaces%ns)) &
       !$ACC COPYIN(surfaces%z_mo(1:surfaces%ns)) &
       !$ACC COPYIN(surfaces%uvw_abs(1:surfaces%ns)) &
       !$ACC COPYIN(surfaces%z0(1:surfaces%ns)) &
       !$ACC COPYIN(surfaces%us(1:surfaces%ns)) &
       !$ACC COPYIN(surfaces%ol(1:surfaces%ns)) &
       !$ACC COPYIN(surfaces%rib(1:surfaces%ns)) &
       !$ACC COPYIN(surfaces%usws(1:surfaces%ns)) &
       !$ACC COPYIN(surfaces%vsws(1:surfaces%ns)) &
       !$ACC COPYIN(surfaces%ts(1:surfaces%ns)) &
       !$ACC COPYIN(surfaces%shf(1:surfaces%ns)) &
       !$ACC COPYIN(surfaces%pt1(1:surfaces%ns)) &
       !$ACC COPYIN(surfaces%qv1(1:surfaces%ns)) &
       !$ACC COPYIN(surfaces%pt_surface(1:surfaces%ns))

       IF ( .NOT. constant_diffusion )  THEN
          !$ACC ENTER DATA &
          !$ACC COPYIN(surfaces%u_0(1:surfaces%ns)) &
          !$ACC COPYIN(surfaces%v_0(1:surfaces%ns))
       ENDIF

    END SUBROUTINE enter_surface_attributes_h
#endif

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Deallocating memory for model-top fluxes  
!------------------------------------------------------------------------------!
    SUBROUTINE deallocate_surface_attributes_h_top( surfaces )

       IMPLICIT NONE


       TYPE(surf_type) ::  surfaces !< respective surface type

       DEALLOCATE ( surfaces%start_index )
       DEALLOCATE ( surfaces%end_index )
!
!--    Indices to locate surface (model-top) element
       DEALLOCATE ( surfaces%i )
       DEALLOCATE ( surfaces%j )
       DEALLOCATE ( surfaces%k )

       IF ( .NOT. constant_diffusion )  THEN    
          DEALLOCATE ( surfaces%u_0 )  
          DEALLOCATE ( surfaces%v_0 )
       ENDIF 
!
!--    Vertical momentum fluxes of u and v
       DEALLOCATE ( surfaces%usws )  
       DEALLOCATE ( surfaces%vsws )  
!
!--    Sensible heat flux
       DEALLOCATE ( surfaces%shf )
!
!--    Latent heat flux
       IF ( humidity .OR. coupling_mode == 'ocean_to_atmosphere')  THEN
          DEALLOCATE ( surfaces%qsws )      
       ENDIF 
!
!--    Scalar flux
       IF ( passive_scalar )  THEN
          DEALLOCATE ( surfaces%ssws ) 
       ENDIF 
!
!--    Chemical species flux
       IF ( air_chemistry )  THEN
          DEALLOCATE ( surfaces%cssws ) 
       ENDIF 
!
!--       
       IF ( surf_bulk_cloud_model .AND. surf_microphysics_morrison)  THEN
          DEALLOCATE ( surfaces%qcsws )
          DEALLOCATE ( surfaces%ncsws )
       ENDIF
!
!--       
       IF ( surf_bulk_cloud_model .AND. surf_microphysics_seifert)  THEN
          DEALLOCATE ( surfaces%qrsws )
          DEALLOCATE ( surfaces%nrsws )
       ENDIF
!
!--    Salinity flux
       IF ( ocean_mode )  DEALLOCATE ( surfaces%sasws )

    END SUBROUTINE deallocate_surface_attributes_h_top


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocating memory for model-top fluxes  
!------------------------------------------------------------------------------!
    SUBROUTINE allocate_surface_attributes_h_top( surfaces,                    &
                                                  nys_l, nyn_l, nxl_l, nxr_l )

       IMPLICIT NONE

       INTEGER(iwp) ::  nyn_l  !< north bound of local 2d array start/end_index, is equal to nyn, except for restart-array
       INTEGER(iwp) ::  nys_l  !< south bound of local 2d array start/end_index, is equal to nyn, except for restart-array
       INTEGER(iwp) ::  nxl_l  !< west bound of local 2d array start/end_index, is equal to nyn, except for restart-array
       INTEGER(iwp) ::  nxr_l  !< east bound of local 2d array start/end_index, is equal to nyn, except for restart-array

       TYPE(surf_type) ::  surfaces !< respective surface type

       ALLOCATE ( surfaces%start_index(nys_l:nyn_l,nxl_l:nxr_l) )
       ALLOCATE ( surfaces%end_index(nys_l:nyn_l,nxl_l:nxr_l)   )
       surfaces%start_index = 0
       surfaces%end_index   = -1
!
!--    Indices to locate surface (model-top) element
       ALLOCATE ( surfaces%i(1:surfaces%ns)  )
       ALLOCATE ( surfaces%j(1:surfaces%ns)  )
       ALLOCATE ( surfaces%k(1:surfaces%ns)  )

       IF ( .NOT. constant_diffusion )  THEN    
          ALLOCATE ( surfaces%u_0(1:surfaces%ns) )  
          ALLOCATE ( surfaces%v_0(1:surfaces%ns) )
       ENDIF 
!
!--    Vertical momentum fluxes of u and v
       ALLOCATE ( surfaces%usws(1:surfaces%ns) )  
       ALLOCATE ( surfaces%vsws(1:surfaces%ns) )  
!
!--    Sensible heat flux
       ALLOCATE ( surfaces%shf(1:surfaces%ns) )
!
!--    Latent heat flux
       IF ( humidity .OR. coupling_mode == 'ocean_to_atmosphere')  THEN
          ALLOCATE ( surfaces%qsws(1:surfaces%ns) )      
       ENDIF 
!
!--    Scalar flux
       IF ( passive_scalar )  THEN
          ALLOCATE ( surfaces%ssws(1:surfaces%ns) ) 
       ENDIF 
!
!--    Chemical species flux
       IF ( air_chemistry )  THEN
          ALLOCATE ( surfaces%cssws(1:nvar,1:surfaces%ns) ) 
       ENDIF 
!
!--       
       IF ( surf_bulk_cloud_model .AND. surf_microphysics_morrison)  THEN
          ALLOCATE ( surfaces%qcsws(1:surfaces%ns) )
          ALLOCATE ( surfaces%ncsws(1:surfaces%ns) )
       ENDIF
!
!--       
       IF ( surf_bulk_cloud_model .AND. surf_microphysics_seifert)  THEN
          ALLOCATE ( surfaces%qrsws(1:surfaces%ns) )
          ALLOCATE ( surfaces%nrsws(1:surfaces%ns) )
       ENDIF
!
!--    Salinity flux
       IF ( ocean_mode )  ALLOCATE ( surfaces%sasws(1:surfaces%ns) )

    END SUBROUTINE allocate_surface_attributes_h_top


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Exit memory for model-top fluxes.
!------------------------------------------------------------------------------!
#if defined( _OPENACC )
    SUBROUTINE exit_surface_attributes_h_top( surfaces )

       IMPLICIT NONE
   
       TYPE(surf_type) ::  surfaces  !< respective surface type
   
       !$ACC EXIT DATA &
       !$ACC DELETE(surfaces%start_index(nys:nyn,nxl:nxr)) &
       !$ACC DELETE(surfaces%end_index(nys:nyn,nxl:nxr)) &
       !$ACC DELETE(surfaces%i(1:surfaces%ns)) &
       !$ACC DELETE(surfaces%j(1:surfaces%ns)) &
       !$ACC DELETE(surfaces%k(1:surfaces%ns)) &
       !$ACC DELETE(surfaces%usws(1:surfaces%ns)) &
       !$ACC DELETE(surfaces%vsws(1:surfaces%ns)) &
       !$ACC DELETE(surfaces%shf(1:surfaces%ns))

       IF ( .NOT. constant_diffusion )  THEN
          !$ACC EXIT DATA &
          !$ACC DELETE(surfaces%u_0(1:surfaces%ns)) &
          !$ACC DELETE(surfaces%v_0(1:surfaces%ns))
       ENDIF
   
    END SUBROUTINE exit_surface_attributes_h_top
#endif

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Enter memory for model-top fluxes.
!------------------------------------------------------------------------------!
#if defined( _OPENACC )
    SUBROUTINE enter_surface_attributes_h_top( surfaces )

       IMPLICIT NONE

       TYPE(surf_type) ::  surfaces  !< respective surface type

       !$ACC ENTER DATA &
       !$ACC COPYIN(surfaces%start_index(nys:nyn,nxl:nxr)) &
       !$ACC COPYIN(surfaces%end_index(nys:nyn,nxl:nxr)) &
       !$ACC COPYIN(surfaces%i(1:surfaces%ns)) &
       !$ACC COPYIN(surfaces%j(1:surfaces%ns)) &
       !$ACC COPYIN(surfaces%k(1:surfaces%ns)) &
       !$ACC COPYIN(surfaces%usws(1:surfaces%ns)) &
       !$ACC COPYIN(surfaces%vsws(1:surfaces%ns)) &
       !$ACC COPYIN(surfaces%shf(1:surfaces%ns))

       IF ( .NOT. constant_diffusion )  THEN
          !$ACC ENTER DATA &
          !$ACC COPYIN(surfaces%u_0(1:surfaces%ns)) &
          !$ACC COPYIN(surfaces%v_0(1:surfaces%ns))
       ENDIF

    END SUBROUTINE enter_surface_attributes_h_top
#endif

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Deallocating memory for vertical surface types. 
!------------------------------------------------------------------------------!
    SUBROUTINE deallocate_surface_attributes_v( surfaces )

       IMPLICIT NONE


       TYPE(surf_type) ::  surfaces !< respective surface type

!
!--    Allocate arrays for start and end index of vertical surface type 
!--    for each (j,i)-grid point. This is required in diffion_x, which is
!--    called for each (j,i). In order to find the location where the 
!--    respective flux is store within the surface-type, start- and end- 
!--    index are stored for each (j,i). For example, each (j,i) can have
!--    several entries where fluxes for vertical surfaces might be stored.  
!--    In the flat case, where no vertical walls exit, set indicies such 
!--    that loop in diffusion routines will not be entered. 
       DEALLOCATE ( surfaces%start_index )
       DEALLOCATE ( surfaces%end_index )
!
!--    Indices to locate surface element.
       DEALLOCATE ( surfaces%i )
       DEALLOCATE ( surfaces%j )
       DEALLOCATE ( surfaces%k )
!
!--    Surface-layer height
       DEALLOCATE ( surfaces%z_mo )
!
!--    Surface orientation
       DEALLOCATE ( surfaces%facing )
!
!--    Surface parallel wind velocity
       DEALLOCATE ( surfaces%uvw_abs )
!
!--    Roughness
       DEALLOCATE ( surfaces%z0 )
       DEALLOCATE ( surfaces%z0h )
       DEALLOCATE ( surfaces%z0q )

!
!--    Friction velocity
       DEALLOCATE ( surfaces%us )
!
!--    Allocate Obukhov length and bulk Richardson number. Actually, at
!--    vertical surfaces these are only required for natural surfaces.  
!--    for natural land surfaces
       DEALLOCATE( surfaces%ol ) 
       DEALLOCATE( surfaces%rib ) 
!
!--    Allocate arrays for surface momentum fluxes for u and v. For u at north- 
!--    and south-facing surfaces, for v at east- and west-facing surfaces.
       DEALLOCATE ( surfaces%mom_flux_uv )
!
!--    Allocate array for surface momentum flux for w - wsus and wsvs
       DEALLOCATE ( surfaces%mom_flux_w ) 
!
!--    Allocate array for surface momentum flux for subgrid-scale tke wsus and 
!--    wsvs; first index usvs or vsws, second index for wsus or wsvs, depending
!--    on surface.
       DEALLOCATE ( surfaces%mom_flux_tke )  
!
!--    Characteristic temperature and surface flux of sensible heat
       DEALLOCATE ( surfaces%ts )    
       DEALLOCATE ( surfaces%shf )
!
!--    surface temperature
       DEALLOCATE ( surfaces%pt_surface ) 
!
!--    Characteristic humidity and surface flux of latent heat
       IF ( humidity )  THEN
          DEALLOCATE ( surfaces%qs ) 
          DEALLOCATE ( surfaces%qsws )  
          DEALLOCATE ( surfaces%q_surface   )
          DEALLOCATE ( surfaces%vpt_surface )
       ENDIF 
!
!--    Characteristic scalar and surface flux of scalar
       IF ( passive_scalar )  THEN
          DEALLOCATE ( surfaces%ss )   
          DEALLOCATE ( surfaces%ssws ) 
       ENDIF
!
!--    Scaling parameter (cs*) and surface flux of chemical species
       IF ( air_chemistry )  THEN
             DEALLOCATE ( surfaces%css )   
             DEALLOCATE ( surfaces%cssws ) 
       ENDIF 
!
!--    Arrays for storing potential temperature and
!--    mixing ratio at first grid level
       DEALLOCATE ( surfaces%pt1 )
       DEALLOCATE ( surfaces%qv1 )
       DEALLOCATE ( surfaces%vpt1 )

       IF ( surf_bulk_cloud_model .AND. surf_microphysics_morrison)  THEN
          DEALLOCATE ( surfaces%qcs )
          DEALLOCATE ( surfaces%ncs )
          DEALLOCATE ( surfaces%qcsws )
          DEALLOCATE ( surfaces%ncsws )
       ENDIF

       IF ( surf_bulk_cloud_model .AND. surf_microphysics_seifert)  THEN
          DEALLOCATE ( surfaces%qrs )
          DEALLOCATE ( surfaces%nrs )
          DEALLOCATE ( surfaces%qrsws )
          DEALLOCATE ( surfaces%nrsws )
       ENDIF
!
!--    Salinity surface flux
       IF ( ocean_mode )  DEALLOCATE ( surfaces%sasws )

    END SUBROUTINE deallocate_surface_attributes_v


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocating memory for vertical surface types. 
!------------------------------------------------------------------------------!
    SUBROUTINE allocate_surface_attributes_v( surfaces,                        &
                                              nys_l, nyn_l, nxl_l, nxr_l )

       IMPLICIT NONE

       INTEGER(iwp) ::  nyn_l  !< north bound of local 2d array start/end_index, is equal to nyn, except for restart-array
       INTEGER(iwp) ::  nys_l  !< south bound of local 2d array start/end_index, is equal to nyn, except for restart-array
       INTEGER(iwp) ::  nxl_l  !< west bound of local 2d array start/end_index, is equal to nyn, except for restart-array
       INTEGER(iwp) ::  nxr_l  !< east bound of local 2d array start/end_index, is equal to nyn, except for restart-array

       TYPE(surf_type) ::  surfaces !< respective surface type

!
!--    Allocate arrays for start and end index of vertical surface type 
!--    for each (j,i)-grid point. This is required in diffion_x, which is
!--    called for each (j,i). In order to find the location where the 
!--    respective flux is store within the surface-type, start- and end- 
!--    index are stored for each (j,i). For example, each (j,i) can have
!--    several entries where fluxes for vertical surfaces might be stored.  
!--    In the flat case, where no vertical walls exit, set indicies such 
!--    that loop in diffusion routines will not be entered. 
       ALLOCATE ( surfaces%start_index(nys_l:nyn_l,nxl_l:nxr_l) )
       ALLOCATE ( surfaces%end_index(nys_l:nyn_l,nxl_l:nxr_l)   )
       surfaces%start_index = 0
       surfaces%end_index   = -1
!
!--    Indices to locate surface element.
       ALLOCATE ( surfaces%i(1:surfaces%ns) )
       ALLOCATE ( surfaces%j(1:surfaces%ns) )
       ALLOCATE ( surfaces%k(1:surfaces%ns) )
!
!--    Surface-layer height
       ALLOCATE ( surfaces%z_mo(1:surfaces%ns) )
!
!--    Surface orientation
       ALLOCATE ( surfaces%facing(1:surfaces%ns) )
!
!--    Surface parallel wind velocity
       ALLOCATE ( surfaces%uvw_abs(1:surfaces%ns) )
!
!--    Roughness
       ALLOCATE ( surfaces%z0(1:surfaces%ns)  )
       ALLOCATE ( surfaces%z0h(1:surfaces%ns) )
       ALLOCATE ( surfaces%z0q(1:surfaces%ns) )

!
!--    Friction velocity
       ALLOCATE ( surfaces%us(1:surfaces%ns) )
!
!--    Allocate Obukhov length and bulk Richardson number. Actually, at
!--    vertical surfaces these are only required for natural surfaces.  
!--    for natural land surfaces
       ALLOCATE( surfaces%ol(1:surfaces%ns)  ) 
       ALLOCATE( surfaces%rib(1:surfaces%ns) ) 
!
!--    Allocate arrays for surface momentum fluxes for u and v. For u at north- 
!--    and south-facing surfaces, for v at east- and west-facing surfaces.
       ALLOCATE ( surfaces%mom_flux_uv(1:surfaces%ns) )
!
!--    Allocate array for surface momentum flux for w - wsus and wsvs
       ALLOCATE ( surfaces%mom_flux_w(1:surfaces%ns) ) 
!
!--    Allocate array for surface momentum flux for subgrid-scale tke wsus and 
!--    wsvs; first index usvs or vsws, second index for wsus or wsvs, depending
!--    on surface.
       ALLOCATE ( surfaces%mom_flux_tke(0:1,1:surfaces%ns) )  
!
!--    Characteristic temperature and surface flux of sensible heat
       ALLOCATE ( surfaces%ts(1:surfaces%ns)  )    
       ALLOCATE ( surfaces%shf(1:surfaces%ns) )
!
!--    surface temperature
       ALLOCATE ( surfaces%pt_surface(1:surfaces%ns) ) 
!
!--    Characteristic humidity and surface flux of latent heat
       IF ( humidity )  THEN
          ALLOCATE ( surfaces%qs(1:surfaces%ns)          ) 
          ALLOCATE ( surfaces%qsws(1:surfaces%ns)        )   
          ALLOCATE ( surfaces%q_surface(1:surfaces%ns)   )
          ALLOCATE ( surfaces%vpt_surface(1:surfaces%ns) )            
       ENDIF 
!
!--    Characteristic scalar and surface flux of scalar
       IF ( passive_scalar )  THEN
          ALLOCATE ( surfaces%ss(1:surfaces%ns)   )   
          ALLOCATE ( surfaces%ssws(1:surfaces%ns) ) 
       ENDIF
!
!--    Scaling parameter (cs*) and surface flux of chemical species
       IF ( air_chemistry )  THEN
             ALLOCATE ( surfaces%css(1:nvar,1:surfaces%ns)   )   
             ALLOCATE ( surfaces%cssws(1:nvar,1:surfaces%ns) ) 
       ENDIF 
!
!--    Arrays for storing potential temperature and
!--    mixing ratio at first grid level
       ALLOCATE ( surfaces%pt1(1:surfaces%ns) )
       ALLOCATE ( surfaces%qv1(1:surfaces%ns) )
       ALLOCATE ( surfaces%vpt1(1:surfaces%ns) )

       IF ( surf_bulk_cloud_model .AND. surf_microphysics_morrison)  THEN
          ALLOCATE ( surfaces%qcs(1:surfaces%ns)   )
          ALLOCATE ( surfaces%ncs(1:surfaces%ns)   )
          ALLOCATE ( surfaces%qcsws(1:surfaces%ns) )
          ALLOCATE ( surfaces%ncsws(1:surfaces%ns) )
       ENDIF

       IF ( surf_bulk_cloud_model .AND. surf_microphysics_seifert)  THEN
          ALLOCATE ( surfaces%qrs(1:surfaces%ns)   )
          ALLOCATE ( surfaces%nrs(1:surfaces%ns)   )
          ALLOCATE ( surfaces%qrsws(1:surfaces%ns) )
          ALLOCATE ( surfaces%nrsws(1:surfaces%ns) )
       ENDIF
!
!--    Salinity surface flux
       IF ( ocean_mode )  ALLOCATE ( surfaces%sasws(1:surfaces%ns) )

    END SUBROUTINE allocate_surface_attributes_v


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Exit memory for vertical surface types. 
!------------------------------------------------------------------------------!
#if defined( _OPENACC )
    SUBROUTINE exit_surface_attributes_v( surfaces )

       IMPLICIT NONE

       TYPE(surf_type) ::  surfaces  !< respective surface type

       !$ACC EXIT DATA &
       !$ACC DELETE(surfaces%start_index(nys:nyn,nxl:nxr)) &
       !$ACC DELETE(surfaces%end_index(nys:nyn,nxl:nxr)) &
       !$ACC DELETE(surfaces%i(1:surfaces%ns)) &
       !$ACC DELETE(surfaces%j(1:surfaces%ns)) &
       !$ACC DELETE(surfaces%k(1:surfaces%ns)) &
       !$ACC DELETE(surfaces%uvw_abs(1:surfaces%ns)) &
       !$ACC DELETE(surfaces%z0(1:surfaces%ns)) &
       !$ACC DELETE(surfaces%rib(1:surfaces%ns)) &
       !$ACC DELETE(surfaces%mom_flux_uv(1:surfaces%ns)) &
       !$ACC DELETE(surfaces%mom_flux_w(1:surfaces%ns)) &
       !$ACC DELETE(surfaces%mom_flux_tke(0:1,1:surfaces%ns)) &
       !$ACC DELETE(surfaces%ts(1:surfaces%ns)) &
       !$ACC DELETE(surfaces%shf(1:surfaces%ns)) &
       !$ACC DELETE(surfaces%pt1(1:surfaces%ns)) &
       !$ACC DELETE(surfaces%qv1(1:surfaces%ns))

    END SUBROUTINE exit_surface_attributes_v
#endif

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Enter memory for vertical surface types. 
!------------------------------------------------------------------------------!
#if defined( _OPENACC )
    SUBROUTINE enter_surface_attributes_v( surfaces )
   
       IMPLICIT NONE
   
       TYPE(surf_type) ::  surfaces  !< respective surface type
   
       !$ACC ENTER DATA &
       !$ACC COPYIN(surfaces%start_index(nys:nyn,nxl:nxr)) &
       !$ACC COPYIN(surfaces%end_index(nys:nyn,nxl:nxr)) &
       !$ACC COPYIN(surfaces%i(1:surfaces%ns)) &
       !$ACC COPYIN(surfaces%j(1:surfaces%ns)) &
       !$ACC COPYIN(surfaces%k(1:surfaces%ns)) &
       !$ACC COPYIN(surfaces%uvw_abs(1:surfaces%ns)) &
       !$ACC COPYIN(surfaces%z0(1:surfaces%ns)) &
       !$ACC COPYIN(surfaces%rib(1:surfaces%ns)) &
       !$ACC COPYIN(surfaces%mom_flux_uv(1:surfaces%ns)) &
       !$ACC COPYIN(surfaces%mom_flux_w(1:surfaces%ns)) &
       !$ACC COPYIN(surfaces%mom_flux_tke(0:1,1:surfaces%ns)) &
       !$ACC COPYIN(surfaces%ts(1:surfaces%ns)) &
       !$ACC COPYIN(surfaces%shf(1:surfaces%ns)) &
       !$ACC COPYIN(surfaces%pt1(1:surfaces%ns)) &
       !$ACC COPYIN(surfaces%qv1(1:surfaces%ns))
   
    END SUBROUTINE enter_surface_attributes_v
#endif

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize surface elements, i.e. set initial values for surface fluxes, 
!> friction velocity, calcuation of start/end indices, etc. . 
!> Please note, further initialization concerning 
!> special surface characteristics, e.g. soil- and vegatation type, 
!> building type, etc., is done in the land-surface and urban-surface module, 
!> respectively.  
!------------------------------------------------------------------------------!
    SUBROUTINE init_surfaces

       IMPLICIT NONE

       INTEGER(iwp) ::  i         !< running index x-direction
       INTEGER(iwp) ::  j         !< running index y-direction
       INTEGER(iwp) ::  k         !< running index z-direction

       INTEGER(iwp)                 ::  start_index_lsm_h !< dummy to determing local start index in surface type for given (j,i), for horizontal natural surfaces
       INTEGER(iwp)                 ::  start_index_usm_h !< dummy to determing local start index in surface type for given (j,i), for horizontal urban surfaces

       INTEGER(iwp)                 ::  num_lsm_h     !< current number of horizontal surface element, natural type
       INTEGER(iwp)                 ::  num_lsm_h_kji !< dummy to determing local end index in surface type for given (j,i), for for horizonal natural surfaces
       INTEGER(iwp)                 ::  num_usm_h     !< current number of horizontal surface element, urban type
       INTEGER(iwp)                 ::  num_usm_h_kji !< dummy to determing local end index in surface type for given (j,i), for for horizonal urban surfaces

       INTEGER(iwp), DIMENSION(0:2) ::  num_def_h     !< current number of horizontal surface element, default type
       INTEGER(iwp), DIMENSION(0:2) ::  num_def_h_kji !< dummy to determing local end index in surface type for given (j,i), for horizonal default surfaces
       INTEGER(iwp), DIMENSION(0:2) ::  start_index_def_h !< dummy to determing local start index in surface type for given (j,i), for horizontal default surfaces
     
       INTEGER(iwp), DIMENSION(0:3) ::  num_def_v     !< current number of vertical surface element, default type
       INTEGER(iwp), DIMENSION(0:3) ::  num_def_v_kji !< dummy to determing local end index in surface type for given (j,i), for vertical default surfaces
       INTEGER(iwp), DIMENSION(0:3) ::  num_lsm_v     !< current number of vertical surface element, natural type
       INTEGER(iwp), DIMENSION(0:3) ::  num_lsm_v_kji !< dummy to determing local end index in surface type for given (j,i), for vertical natural surfaces
       INTEGER(iwp), DIMENSION(0:3) ::  num_usm_v     !< current number of vertical surface element, urban type
       INTEGER(iwp), DIMENSION(0:3) ::  num_usm_v_kji !< dummy to determing local end index in surface type for given (j,i), for vertical urban surfaces

       INTEGER(iwp), DIMENSION(0:3) ::  start_index_def_v !< dummy to determing local start index in surface type for given (j,i), for vertical default surfaces
       INTEGER(iwp), DIMENSION(0:3) ::  start_index_lsm_v !< dummy to determing local start index in surface type for given (j,i), for vertical natural surfaces
       INTEGER(iwp), DIMENSION(0:3) ::  start_index_usm_v !< dummy to determing local start index in surface type for given (j,i), for vertical urban surfaces

       LOGICAL ::  building            !< flag indicating building grid point
       LOGICAL ::  terrain             !< flag indicating natural terrain grid point
       LOGICAL ::  unresolved_building !< flag indicating a grid point where actually a building is defined but not resolved by the vertical grid 
!
!--    Set offset indices, i.e. index difference between surface element and 
!--    surface-bounded grid point.
!--    Upward facing - no horizontal offsets
       surf_def_h(0:2)%ioff = 0
       surf_def_h(0:2)%joff = 0

       surf_lsm_h%ioff = 0
       surf_lsm_h%joff = 0

       surf_usm_h%ioff = 0
       surf_usm_h%joff = 0
!
!--    Upward facing vertical offsets
       surf_def_h(0)%koff   = -1
       surf_lsm_h%koff      = -1
       surf_usm_h%koff      = -1
!
!--    Downward facing vertical offset
       surf_def_h(1:2)%koff = 1
!
!--    Vertical surfaces - no vertical offset
       surf_def_v(0:3)%koff = 0
       surf_lsm_v(0:3)%koff = 0
       surf_usm_v(0:3)%koff = 0
!
!--    North- and southward facing - no offset in x
       surf_def_v(0:1)%ioff = 0
       surf_lsm_v(0:1)%ioff = 0
       surf_usm_v(0:1)%ioff = 0
!
!--    Northward facing offset in y
       surf_def_v(0)%joff = -1
       surf_lsm_v(0)%joff = -1
       surf_usm_v(0)%joff = -1
!
!--    Southward facing offset in y
       surf_def_v(1)%joff = 1
       surf_lsm_v(1)%joff = 1
       surf_usm_v(1)%joff = 1

!
!--    East- and westward facing - no offset in y
       surf_def_v(2:3)%joff = 0
       surf_lsm_v(2:3)%joff = 0
       surf_usm_v(2:3)%joff = 0
!
!--    Eastward facing offset in x
       surf_def_v(2)%ioff = -1
       surf_lsm_v(2)%ioff = -1
       surf_usm_v(2)%ioff = -1
!
!--    Westward facing offset in y
       surf_def_v(3)%ioff = 1
       surf_lsm_v(3)%ioff = 1
       surf_usm_v(3)%ioff = 1

!
!--    Initialize surface attributes, store indicies, surfaces orientation, etc., 
       num_def_h(0:2) = 1
       num_def_v(0:3) = 1

       num_lsm_h      = 1
       num_lsm_v(0:3) = 1

       num_usm_h      = 1
       num_usm_v(0:3) = 1

       start_index_def_h(0:2) = 1
       start_index_def_v(0:3) = 1

       start_index_lsm_h      = 1
       start_index_lsm_v(0:3) = 1

       start_index_usm_h      = 1
       start_index_usm_v(0:3) = 1

       DO  i = nxl, nxr
          DO  j = nys, nyn

             num_def_h_kji = 0
             num_def_v_kji = 0
             num_lsm_h_kji = 0
             num_lsm_v_kji = 0
             num_usm_h_kji = 0
             num_usm_v_kji = 0

             DO  k = nzb+1, nzt
!
!--             Check if current gridpoint belongs to the atmosphere
                IF ( BTEST( wall_flags_total_0(k,j,i), 0 ) )  THEN
!
!--                Upward-facing surface. Distinguish between differet surface types.
!--                To do, think about method to flag natural and non-natural 
!--                surfaces. 
                   IF ( .NOT. BTEST( wall_flags_total_0(k-1,j,i), 0 ) )  THEN 
!
!--                   Determine flags indicating terrain or building
                      terrain  = BTEST( wall_flags_total_0(k-1,j,i), 5 )  .OR.       &
                                 topo_no_distinct
                      building = BTEST( wall_flags_total_0(k-1,j,i), 6 )  .OR.       &
                                 topo_no_distinct
                                 
!
!--                   unresolved_building indicates a surface with equal height 
!--                   as terrain but with a non-grid resolved building on top.
!--                   These surfaces will be flagged as urban surfaces.
                      unresolved_building = BTEST( wall_flags_total_0(k-1,j,i), 5 )  &
                                     .AND.  BTEST( wall_flags_total_0(k-1,j,i), 6 )
!
!--                   Natural surface type          
                      IF ( land_surface  .AND.  terrain  .AND.                 &
                           .NOT. unresolved_building )  THEN
                         CALL initialize_horizontal_surfaces( k, j, i,         &
                                                              surf_lsm_h,      &
                                                              num_lsm_h,       &
                                                              num_lsm_h_kji,   &
                                                              .TRUE., .FALSE. )  
!
!--                   Urban surface tpye
                      ELSEIF ( urban_surface  .AND.  building )  THEN
                         CALL initialize_horizontal_surfaces( k, j, i,         &
                                                              surf_usm_h,      &
                                                              num_usm_h,       &
                                                              num_usm_h_kji,   &
                                                              .TRUE., .FALSE. )  
!
!--                   Default surface type
                      ELSE
                         CALL initialize_horizontal_surfaces( k, j, i,         &
                                                              surf_def_h(0),   &
                                                              num_def_h(0),    &
                                                              num_def_h_kji(0),&
                                                              .TRUE., .FALSE. )  
                      ENDIF
                   ENDIF  
!
!--                downward-facing surface, first, model top. Please note, 
!--                for the moment, downward-facing surfaces are always of 
!--                default type
                   IF ( k == nzt  .AND.  use_top_fluxes )  THEN
                      CALL initialize_top( k, j, i, surf_def_h(2),             &
                                           num_def_h(2), num_def_h_kji(2) )
!
!--                Check for any other downward-facing surface. So far only for 
!--                default surface type.
                   ELSEIF ( .NOT. BTEST( wall_flags_total_0(k+1,j,i), 0 ) )  THEN
                      CALL initialize_horizontal_surfaces( k, j, i,            &
                                                           surf_def_h(1),      &
                                                           num_def_h(1),       &
                                                           num_def_h_kji(1),   &
                                                           .FALSE., .TRUE. )    
                   ENDIF 
!
!--                Check for vertical walls and, if required, initialize it.
!                  Start with northward-facing surface.
                   IF ( .NOT. BTEST( wall_flags_total_0(k,j-1,i), 0 ) )  THEN
!
!--                   Determine flags indicating terrain or building
                      terrain  = BTEST( wall_flags_total_0(k,j-1,i), 5 )  .OR.       &
                                 topo_no_distinct
                      building = BTEST( wall_flags_total_0(k,j-1,i), 6 )  .OR.       &
                                 topo_no_distinct

                      unresolved_building = BTEST( wall_flags_total_0(k,j-1,i), 5 )  &
                                     .AND.  BTEST( wall_flags_total_0(k,j-1,i), 6 )
                                     
                      IF ( land_surface  .AND.  terrain  .AND.                 &
                           .NOT. unresolved_building )  THEN
                         CALL initialize_vertical_surfaces( k, j, i,           &
                                                            surf_lsm_v(0),     &
                                                            num_lsm_v(0),      &
                                                            num_lsm_v_kji(0),  &
                                                            .FALSE., .FALSE.,  &             
                                                            .FALSE., .TRUE. )
                      ELSEIF ( urban_surface  .AND.  building )  THEN
                         CALL initialize_vertical_surfaces( k, j, i,           &
                                                            surf_usm_v(0),     &
                                                            num_usm_v(0),      &
                                                            num_usm_v_kji(0),  &
                                                            .FALSE., .FALSE.,  &             
                                                            .FALSE., .TRUE. )
                      ELSE
                         CALL initialize_vertical_surfaces( k, j, i,           &
                                                            surf_def_v(0),     &
                                                            num_def_v(0),      &
                                                            num_def_v_kji(0),  &
                                                            .FALSE., .FALSE.,  &             
                                                            .FALSE., .TRUE. ) 
                      ENDIF
                   ENDIF
!
!--                southward-facing surface
                   IF ( .NOT. BTEST( wall_flags_total_0(k,j+1,i), 0 ) )  THEN
!
!--                   Determine flags indicating terrain or building
                      terrain  = BTEST( wall_flags_total_0(k,j+1,i), 5 )  .OR.       &
                                 topo_no_distinct
                      building = BTEST( wall_flags_total_0(k,j+1,i), 6 )  .OR.       &
                                 topo_no_distinct
                                 
                      unresolved_building = BTEST( wall_flags_total_0(k,j+1,i), 5 )  &
                                     .AND.  BTEST( wall_flags_total_0(k,j+1,i), 6 )
                                     
                      IF ( land_surface  .AND.  terrain  .AND.                 &
                           .NOT. unresolved_building )  THEN
                         CALL initialize_vertical_surfaces( k, j, i,           &
                                                            surf_lsm_v(1),     &
                                                            num_lsm_v(1),      &
                                                            num_lsm_v_kji(1),  &
                                                            .FALSE., .FALSE.,  &
                                                            .TRUE., .FALSE. )
                      ELSEIF ( urban_surface  .AND.  building )  THEN
                         CALL initialize_vertical_surfaces( k, j, i,           &
                                                            surf_usm_v(1),     &
                                                            num_usm_v(1),      &
                                                            num_usm_v_kji(1),  &
                                                            .FALSE., .FALSE.,  &
                                                            .TRUE., .FALSE. )
                      ELSE
                         CALL initialize_vertical_surfaces( k, j, i,           &
                                                            surf_def_v(1),     &
                                                            num_def_v(1),      &
                                                            num_def_v_kji(1),  &
                                                            .FALSE., .FALSE.,  &
                                                            .TRUE., .FALSE. ) 
                      ENDIF
                   ENDIF
!
!--                eastward-facing surface
                   IF ( .NOT. BTEST( wall_flags_total_0(k,j,i-1), 0 ) )  THEN
!
!--                   Determine flags indicating terrain or building
                      terrain  = BTEST( wall_flags_total_0(k,j,i-1), 5 )  .OR.       &
                                 topo_no_distinct
                      building = BTEST( wall_flags_total_0(k,j,i-1), 6 )  .OR.       &
                                 topo_no_distinct
                                 
                      unresolved_building = BTEST( wall_flags_total_0(k,j,i-1), 5 )  &
                                     .AND.  BTEST( wall_flags_total_0(k,j,i-1), 6 )
                                 
                      IF ( land_surface  .AND.  terrain  .AND.                 &
                           .NOT. unresolved_building )  THEN
                         CALL initialize_vertical_surfaces( k, j, i,           &
                                                            surf_lsm_v(2),     &
                                                            num_lsm_v(2),      &
                                                            num_lsm_v_kji(2),  &
                                                            .TRUE., .FALSE.,   &
                                                            .FALSE., .FALSE. )
                      ELSEIF ( urban_surface  .AND.  building )  THEN
                         CALL initialize_vertical_surfaces( k, j, i,           &
                                                            surf_usm_v(2),     &
                                                            num_usm_v(2),      &
                                                            num_usm_v_kji(2),  &
                                                            .TRUE., .FALSE.,   &
                                                            .FALSE., .FALSE. )
                      ELSE
                         CALL initialize_vertical_surfaces( k, j, i,           &
                                                            surf_def_v(2),     &
                                                            num_def_v(2),      &
                                                            num_def_v_kji(2),  &
                                                            .TRUE., .FALSE.,   &
                                                            .FALSE., .FALSE. ) 
                      ENDIF
                   ENDIF 
!   
!--                westward-facing surface
                   IF ( .NOT. BTEST( wall_flags_total_0(k,j,i+1), 0 ) )  THEN
!
!--                   Determine flags indicating terrain or building
                      terrain  = BTEST( wall_flags_total_0(k,j,i+1), 5 )  .OR.       &
                                 topo_no_distinct
                      building = BTEST( wall_flags_total_0(k,j,i+1), 6 )  .OR.       &
                                 topo_no_distinct
                                 
                      unresolved_building = BTEST( wall_flags_total_0(k,j,i+1), 5 )  &
                                     .AND.  BTEST( wall_flags_total_0(k,j,i+1), 6 )
                                  
                      IF ( land_surface  .AND.  terrain  .AND.                 &
                           .NOT. unresolved_building )  THEN
                         CALL initialize_vertical_surfaces( k, j, i,           &
                                                            surf_lsm_v(3),     &
                                                            num_lsm_v(3),      &
                                                            num_lsm_v_kji(3),  &
                                                           .FALSE., .TRUE.,    &
                                                           .FALSE., .FALSE. )
                      ELSEIF ( urban_surface  .AND.  building )  THEN
                         CALL initialize_vertical_surfaces( k, j, i,           &
                                                            surf_usm_v(3),     &
                                                            num_usm_v(3),      &
                                                            num_usm_v_kji(3),  &
                                                           .FALSE., .TRUE.,    &
                                                           .FALSE., .FALSE. )
                      ELSE
                         CALL initialize_vertical_surfaces( k, j, i,           &
                                                            surf_def_v(3),     &
                                                            num_def_v(3),      &
                                                            num_def_v_kji(3),  &
                                                           .FALSE., .TRUE.,    &
                                                           .FALSE., .FALSE. ) 
                      ENDIF
                   ENDIF
                ENDIF

 
             ENDDO
!
!--          Determine start- and end-index at grid point (j,i). Also, for 
!--          horizontal surfaces more than 1 horizontal surface element can 
!--          exist at grid point (j,i) if overhanging structures are present.
!--          Upward-facing surfaces
             surf_def_h(0)%start_index(j,i) = start_index_def_h(0)
             surf_def_h(0)%end_index(j,i)   = surf_def_h(0)%start_index(j,i) + &
                                                 num_def_h_kji(0) - 1
             start_index_def_h(0)           = surf_def_h(0)%end_index(j,i) + 1
!
!--          ATTENTION:
!--          workaround to prevent vectorization bug on NEC Aurora
             IF ( start_index_def_h(0) < -99999 )  THEN
                PRINT*, 'i=', i, ' j=',j, ' s=',surf_def_h(0)%start_index(j,i),                    &
                        ' e=', surf_def_h(0)%end_index(j,i)
             ENDIF
!
!--          Downward-facing surfaces, except model top
             surf_def_h(1)%start_index(j,i) = start_index_def_h(1)                                                 
             surf_def_h(1)%end_index(j,i)   = surf_def_h(1)%start_index(j,i) + &
                                                 num_def_h_kji(1) - 1
             start_index_def_h(1)           = surf_def_h(1)%end_index(j,i) + 1
!
!--          Downward-facing surfaces -- model top fluxes
             surf_def_h(2)%start_index(j,i) = start_index_def_h(2)                                                 
             surf_def_h(2)%end_index(j,i)   = surf_def_h(2)%start_index(j,i) + &
                                                 num_def_h_kji(2) - 1
             start_index_def_h(2)           = surf_def_h(2)%end_index(j,i) + 1
!
!--          Horizontal natural land surfaces
             surf_lsm_h%start_index(j,i)    = start_index_lsm_h
             surf_lsm_h%end_index(j,i)      = surf_lsm_h%start_index(j,i) +    &
                                                 num_lsm_h_kji - 1
             start_index_lsm_h              = surf_lsm_h%end_index(j,i) + 1
!
!--          Horizontal urban surfaces
             surf_usm_h%start_index(j,i)    = start_index_usm_h
             surf_usm_h%end_index(j,i)      = surf_usm_h%start_index(j,i) +    &
                                                 num_usm_h_kji - 1
             start_index_usm_h              = surf_usm_h%end_index(j,i) + 1

!
!--          Vertical surfaces - Default type
             surf_def_v(0)%start_index(j,i) = start_index_def_v(0)
             surf_def_v(1)%start_index(j,i) = start_index_def_v(1)
             surf_def_v(2)%start_index(j,i) = start_index_def_v(2)
             surf_def_v(3)%start_index(j,i) = start_index_def_v(3)
             surf_def_v(0)%end_index(j,i)   = start_index_def_v(0) +           & 
                                              num_def_v_kji(0) - 1
             surf_def_v(1)%end_index(j,i)   = start_index_def_v(1) +           &
                                              num_def_v_kji(1) - 1
             surf_def_v(2)%end_index(j,i)   = start_index_def_v(2) +           &
                                              num_def_v_kji(2) - 1
             surf_def_v(3)%end_index(j,i)   = start_index_def_v(3) +           &
                                              num_def_v_kji(3) - 1
             start_index_def_v(0)           = surf_def_v(0)%end_index(j,i) + 1
             start_index_def_v(1)           = surf_def_v(1)%end_index(j,i) + 1
             start_index_def_v(2)           = surf_def_v(2)%end_index(j,i) + 1
             start_index_def_v(3)           = surf_def_v(3)%end_index(j,i) + 1
!
!--          Natural type
             surf_lsm_v(0)%start_index(j,i) = start_index_lsm_v(0)
             surf_lsm_v(1)%start_index(j,i) = start_index_lsm_v(1)
             surf_lsm_v(2)%start_index(j,i) = start_index_lsm_v(2)
             surf_lsm_v(3)%start_index(j,i) = start_index_lsm_v(3)
             surf_lsm_v(0)%end_index(j,i)   = start_index_lsm_v(0) +           & 
                                              num_lsm_v_kji(0) - 1
             surf_lsm_v(1)%end_index(j,i)   = start_index_lsm_v(1) +           &
                                              num_lsm_v_kji(1) - 1
             surf_lsm_v(2)%end_index(j,i)   = start_index_lsm_v(2) +           &
                                              num_lsm_v_kji(2) - 1
             surf_lsm_v(3)%end_index(j,i)   = start_index_lsm_v(3) +           &
                                              num_lsm_v_kji(3) - 1
             start_index_lsm_v(0)           = surf_lsm_v(0)%end_index(j,i) + 1
             start_index_lsm_v(1)           = surf_lsm_v(1)%end_index(j,i) + 1
             start_index_lsm_v(2)           = surf_lsm_v(2)%end_index(j,i) + 1
             start_index_lsm_v(3)           = surf_lsm_v(3)%end_index(j,i) + 1
!
!--          Urban type
             surf_usm_v(0)%start_index(j,i) = start_index_usm_v(0)
             surf_usm_v(1)%start_index(j,i) = start_index_usm_v(1)
             surf_usm_v(2)%start_index(j,i) = start_index_usm_v(2)
             surf_usm_v(3)%start_index(j,i) = start_index_usm_v(3)
             surf_usm_v(0)%end_index(j,i)   = start_index_usm_v(0) +           & 
                                              num_usm_v_kji(0) - 1
             surf_usm_v(1)%end_index(j,i)   = start_index_usm_v(1) +           &
                                              num_usm_v_kji(1) - 1
             surf_usm_v(2)%end_index(j,i)   = start_index_usm_v(2) +           &
                                              num_usm_v_kji(2) - 1
             surf_usm_v(3)%end_index(j,i)   = start_index_usm_v(3) +           &
                                              num_usm_v_kji(3) - 1
             start_index_usm_v(0)           = surf_usm_v(0)%end_index(j,i) + 1
             start_index_usm_v(1)           = surf_usm_v(1)%end_index(j,i) + 1
             start_index_usm_v(2)           = surf_usm_v(2)%end_index(j,i) + 1
             start_index_usm_v(3)           = surf_usm_v(3)%end_index(j,i) + 1


          ENDDO
       ENDDO

       CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize horizontal surface elements, upward- and downward-facing. 
!> Note, horizontal surface type alsw comprises model-top fluxes, which are,
!> initialized in a different routine. 
!------------------------------------------------------------------------------!
          SUBROUTINE initialize_horizontal_surfaces( k, j, i, surf, num_h,     &
                                                     num_h_kji, upward_facing, &
                                                     downward_facing )       

             IMPLICIT NONE 

             INTEGER(iwp)  ::  i                !< running index x-direction
             INTEGER(iwp)  ::  j                !< running index y-direction
             INTEGER(iwp)  ::  k                !< running index z-direction
             INTEGER(iwp)  ::  num_h            !< current number of surface element
             INTEGER(iwp)  ::  num_h_kji        !< dummy increment
             INTEGER(iwp)  ::  lsp              !< running index chemical species
             INTEGER(iwp)  ::  lsp_pr           !< running index chemical species??

             LOGICAL       ::  upward_facing    !< flag indicating upward-facing surface
             LOGICAL       ::  downward_facing  !< flag indicating downward-facing surface

             TYPE( surf_type ) :: surf          !< respective surface type

!
!--          Store indices of respective surface element
             surf%i(num_h) = i
             surf%j(num_h) = j
             surf%k(num_h) = k
!
!--          Surface orientation, bit 0 is set to 1 for upward-facing surfaces, 
!--          bit 1 is for downward-facing surfaces.
             IF ( upward_facing   )  surf%facing(num_h) = IBSET( surf%facing(num_h), 0 )
             IF ( downward_facing )  surf%facing(num_h) = IBSET( surf%facing(num_h), 1 )
!
!--          Initialize surface-layer height
             IF ( upward_facing )  THEN
                surf%z_mo(num_h)  = zu(k) - zw(k-1)
             ELSE
                surf%z_mo(num_h)  = zw(k) - zu(k)
             ENDIF
 
             surf%z0(num_h)    = roughness_length
             surf%z0h(num_h)   = z0h_factor * roughness_length
             surf%z0q(num_h)   = z0h_factor * roughness_length          
!
!--          Initialization in case of 1D pre-cursor run
             IF ( INDEX( initializing_actions, 'set_1d-model_profiles' ) /= 0 )&
             THEN
                IF ( .NOT. constant_diffusion )  THEN
                   IF ( constant_flux_layer )  THEN
                      surf%ol(num_h)   = surf%z_mo(num_h) /                    &
                                            ( rif1d(nzb+1) + 1.0E-20_wp )
                      surf%us(num_h)   = us1d
                      surf%usws(num_h) = usws1d
                      surf%vsws(num_h) = vsws1d
                   ELSE
                      surf%ol(num_h)   = surf%z_mo(num_h) / zeta_min
                      surf%us(num_h)   = 0.0_wp
                      surf%usws(num_h) = 0.0_wp
                      surf%vsws(num_h) = 0.0_wp
                   ENDIF
                ELSE
                   surf%ol(num_h)   = surf%z_mo(num_h) / zeta_min
                   surf%us(num_h)   = 0.0_wp
                   surf%usws(num_h) = 0.0_wp
                   surf%vsws(num_h) = 0.0_wp
                ENDIF
!
!--          Initialization in all other cases
             ELSE

                surf%ol(num_h)   = surf%z_mo(num_h) / zeta_min
!
!--             Very small number is required for calculation of Obukhov length 
!--             at first timestep     
                surf%us(num_h)    = 1E-30_wp 
                surf%usws(num_h)  = 0.0_wp
                surf%vsws(num_h)  = 0.0_wp
        
             ENDIF

             surf%rib(num_h)     = 0.0_wp 
             surf%uvw_abs(num_h) = 0.0_wp

             IF ( .NOT. constant_diffusion )  THEN    
                surf%u_0(num_h)     = 0.0_wp  
                surf%v_0(num_h)     = 0.0_wp
             ENDIF 

             surf%ts(num_h)   = 0.0_wp
!
!--          Set initial value for surface temperature
             surf%pt_surface(num_h) = pt_surface
             
             IF ( humidity )  THEN
                surf%qs(num_h)   = 0.0_wp
                IF ( surf_bulk_cloud_model .AND. surf_microphysics_morrison)  THEN
                   surf%qcs(num_h) = 0.0_wp
                   surf%ncs(num_h) = 0.0_wp
   
                   surf%qcsws(num_h) = 0.0_wp
                   surf%ncsws(num_h) = 0.0_wp

                ENDIF
                IF ( surf_bulk_cloud_model .AND. surf_microphysics_seifert)  THEN
                   surf%qrs(num_h) = 0.0_wp
                   surf%nrs(num_h) = 0.0_wp
   
                   surf%qrsws(num_h) = 0.0_wp
                   surf%nrsws(num_h) = 0.0_wp

                   surf%pt1(num_h)  = 0.0_wp
                   surf%qv1(num_h)  = 0.0_wp
                   surf%vpt1(num_h) = 0.0_wp
                   
                ENDIF
                
                surf%q_surface(num_h)   = q_surface
                surf%vpt_surface(num_h) = surf%pt_surface(num_h) *             &
                                   ( 1.0_wp + 0.61_wp * surf%q_surface(num_h) )
             ENDIF

             IF ( passive_scalar )  surf%ss(num_h) = 0.0_wp

             DO  lsp = 1, nvar
                IF ( air_chemistry )  surf%css(lsp,num_h)   = 0.0_wp
!
!--             Ensure that fluxes of compounds which are not specified in 
!--             namelist are all zero --> kanani: revise
                IF ( air_chemistry )  surf%cssws(lsp,num_h) = 0.0_wp
             ENDDO
!
!--          Inititalize surface fluxes of sensible and latent heat, as well as
!--          passive scalar
             IF ( use_surface_fluxes )  THEN

                IF ( upward_facing )  THEN
                   IF ( constant_heatflux )  THEN
!   
!--                   Initialize surface heatflux. However, skip this for now if 
!--                   if random_heatflux is set. This case, shf is initialized later.
                      IF ( .NOT. random_heatflux )  THEN
                         surf%shf(num_h) = surface_heatflux *                  &
                                                 heatflux_input_conversion(k-1)
!
!--                      Check if surface heat flux might be replaced by 
!--                      prescribed wall heatflux
                         IF ( k-1 /= 0 )  THEN
                            surf%shf(num_h) = wall_heatflux(0) *               &
                                                 heatflux_input_conversion(k-1)
                         ENDIF
                      ENDIF
                   ELSE
                      surf%shf(num_h) = 0.0_wp
                   ENDIF
!
!--             Set heat-flux at downward-facing surfaces
                ELSE
                   surf%shf(num_h) = wall_heatflux(5) *                        &
                                             heatflux_input_conversion(k)
                ENDIF

                IF ( humidity )  THEN
                   IF ( upward_facing )  THEN
                      IF ( constant_waterflux )  THEN
                         surf%qsws(num_h) = surface_waterflux *                &
                                                 waterflux_input_conversion(k-1)
                         IF ( k-1 /= 0 )  THEN
                            surf%qsws(num_h) = wall_humidityflux(0) *          &
                                                 waterflux_input_conversion(k-1)
                         ENDIF
                      ELSE
                         surf%qsws(num_h) = 0.0_wp
                      ENDIF
                   ELSE
                      surf%qsws(num_h) = wall_humidityflux(5) *                &
                                                waterflux_input_conversion(k)
                   ENDIF
                ENDIF

                IF ( passive_scalar )  THEN
                   IF ( upward_facing )  THEN
                      IF ( constant_scalarflux )  THEN
                         surf%ssws(num_h) = surface_scalarflux  * rho_air_zw(k-1)

                         IF ( k-1 /= 0 )                                       &
                            surf%ssws(num_h) = wall_scalarflux(0) *            &
                                               rho_air_zw(k-1)

                      ELSE
                         surf%ssws(num_h) = 0.0_wp
                      ENDIF
                   ELSE
                      surf%ssws(num_h) = wall_scalarflux(5) * rho_air_zw(k)
                   ENDIF
                ENDIF

                IF ( air_chemistry )  THEN
                   lsp_pr = 1
                   DO  WHILE ( TRIM( surface_csflux_name( lsp_pr ) ) /= 'novalue' )   !<'novalue' is the default
                      DO  lsp = 1, nvar
!
!--                      Assign surface flux for each variable species
                         IF ( TRIM( spc_names(lsp) ) == TRIM( surface_csflux_name(lsp_pr) ) )  THEN   
                            IF ( upward_facing )  THEN
                               IF ( constant_csflux(lsp_pr) )  THEN
                                  surf%cssws(lsp,num_h) =                      &
                                                       surface_csflux(lsp_pr) *&
                                                       rho_air_zw(k-1)

                                  IF ( k-1 /= 0 )                              &
                                     surf%cssws(lsp,num_h) =                   &
                                                       wall_csflux(lsp,0) *    &
                                                       rho_air_zw(k-1) 
                               ELSE
                                  surf%cssws(lsp,num_h) = 0.0_wp
                               ENDIF
                            ELSE
                               surf%cssws(lsp,num_h) = wall_csflux(lsp,5) *    &
                                                       rho_air_zw(k)
                            ENDIF
                         ENDIF
                      ENDDO
                      lsp_pr = lsp_pr + 1
                   ENDDO
                ENDIF

                IF ( ocean_mode )  THEN
                   IF ( upward_facing )  THEN 
                      surf%sasws(num_h) = bottom_salinityflux * rho_air_zw(k-1)
                   ELSE
                      surf%sasws(num_h) = 0.0_wp
                   ENDIF
                ENDIF
             ENDIF
!
!--          Increment surface indices
             num_h     = num_h + 1
             num_h_kji = num_h_kji + 1      


          END SUBROUTINE initialize_horizontal_surfaces
       

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize model-top fluxes. Currently, only the heatflux and salinity flux 
!> can be prescribed, latent flux is zero in this case!
!------------------------------------------------------------------------------!
          SUBROUTINE initialize_top( k, j, i, surf, num_h, num_h_kji )       

             IMPLICIT NONE 

             INTEGER(iwp)  ::  i                !< running index x-direction
             INTEGER(iwp)  ::  j                !< running index y-direction
             INTEGER(iwp)  ::  k                !< running index z-direction
             INTEGER(iwp)  ::  num_h            !< current number of surface element
             INTEGER(iwp)  ::  num_h_kji        !< dummy increment
             INTEGER(iwp)  ::  lsp              !< running index for chemical species

             TYPE( surf_type ) :: surf          !< respective surface type
!
!--          Store indices of respective surface element
             surf%i(num_h) = i
             surf%j(num_h) = j
             surf%k(num_h) = k
!
!--          Initialize top heat flux
             IF ( constant_top_heatflux )                                      &
                surf%shf(num_h) = top_heatflux * heatflux_input_conversion(nzt+1)
!
!--          Initialization in case of a coupled model run
             IF ( coupling_mode == 'ocean_to_atmosphere' )  THEN
                surf%shf(num_h) = 0.0_wp
                surf%qsws(num_h) = 0.0_wp
             ENDIF
!
!--          Prescribe latent heat flux at the top      
             IF ( humidity )  THEN
                surf%qsws(num_h) = 0.0_wp
                IF ( surf_bulk_cloud_model  .AND.  surf_microphysics_morrison ) THEN
                   surf%ncsws(num_h) = 0.0_wp
                   surf%qcsws(num_h) = 0.0_wp
                ENDIF
                IF ( surf_bulk_cloud_model  .AND.  surf_microphysics_seifert ) THEN
                   surf%nrsws(num_h) = 0.0_wp
                   surf%qrsws(num_h) = 0.0_wp
                ENDIF
             ENDIF
!
!--          Prescribe top scalar flux
             IF ( passive_scalar .AND. constant_top_scalarflux )               &
                surf%ssws(num_h) = top_scalarflux * rho_air_zw(nzt+1)
!
!--          Prescribe top chemical species' flux
             DO  lsp = 1, nvar
                IF ( air_chemistry  .AND.  constant_top_csflux(lsp) )  THEN 
                   surf%cssws(lsp,num_h) = top_csflux(lsp) * rho_air_zw(nzt+1)
                ENDIF
             ENDDO
!
!--          Prescribe top salinity flux
             IF ( ocean_mode .AND. constant_top_salinityflux)                  &
                surf%sasws(num_h) = top_salinityflux * rho_air_zw(nzt+1)
!
!--          Top momentum fluxes
             IF ( constant_top_momentumflux )  THEN
                surf%usws(num_h) = top_momentumflux_u *                        &
                                            momentumflux_input_conversion(nzt+1)
                surf%vsws(num_h) = top_momentumflux_v *                        &
                                            momentumflux_input_conversion(nzt+1)
             ENDIF
!
!--          Increment surface indices
             num_h     = num_h + 1
             num_h_kji = num_h_kji + 1      


          END SUBROUTINE initialize_top


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize vertical surface elements. 
!------------------------------------------------------------------------------!
          SUBROUTINE initialize_vertical_surfaces( k, j, i, surf, num_v,       &
                                                   num_v_kji, east_facing,     &
                                                   west_facing, south_facing,  &
                                                   north_facing )       

             IMPLICIT NONE 

             INTEGER(iwp)  ::  component       !< index of wall_fluxes_ array for respective orientation 
             INTEGER(iwp)  ::  i               !< running index x-direction
             INTEGER(iwp)  ::  j               !< running index x-direction
             INTEGER(iwp)  ::  k               !< running index x-direction
             INTEGER(iwp)  ::  num_v           !< current number of surface element
             INTEGER(iwp)  ::  num_v_kji       !< current number of surface element at (j,i)
             INTEGER(iwp)  ::  lsp             !< running index for chemical species


             LOGICAL       ::  east_facing     !< flag indicating east-facing surfaces
             LOGICAL       ::  north_facing    !< flag indicating north-facing surfaces
             LOGICAL       ::  south_facing    !< flag indicating south-facing surfaces
             LOGICAL       ::  west_facing     !< flag indicating west-facing surfaces

             TYPE( surf_type ) :: surf         !< respective surface type

!
!--          Store indices of respective wall element
             surf%i(num_v)   = i
             surf%j(num_v)   = j
             surf%k(num_v)   = k
!
!--          Initialize surface-layer height, or more precisely, distance to surface
             IF ( north_facing  .OR.  south_facing )  THEN
                surf%z_mo(num_v)  = 0.5_wp * dy
             ELSE
                surf%z_mo(num_v)  = 0.5_wp * dx
             ENDIF

             surf%facing(num_v)  = 0
!
!--          Surface orientation. Moreover, set component id to map wall_heatflux, 
!--          etc., on surface type (further below)
             IF ( north_facing )  THEN
                surf%facing(num_v) = 5 !IBSET( surf%facing(num_v), 0 )  
                component          = 4
             ENDIF

             IF ( south_facing )  THEN
                surf%facing(num_v) = 6 !IBSET( surf%facing(num_v), 1 ) 
                component          = 3
             ENDIF

             IF ( east_facing )  THEN
                surf%facing(num_v) = 7 !IBSET( surf%facing(num_v), 2 )
                component          = 2
             ENDIF

             IF ( west_facing )  THEN
                surf%facing(num_v) = 8 !IBSET( surf%facing(num_v), 3 ) 
                component          = 1
             ENDIF

 
             surf%z0(num_v)  = roughness_length
             surf%z0h(num_v) = z0h_factor * roughness_length
             surf%z0q(num_v) = z0h_factor * roughness_length

             surf%us(num_v)  = 0.0_wp
!
!--          If required, initialize Obukhov length
             IF ( ALLOCATED( surf%ol ) )                                       &
                surf%ol(num_v) = surf%z_mo(num_v) / zeta_min

             surf%uvw_abs(num_v)   = 0.0_wp

             surf%mom_flux_uv(num_v) = 0.0_wp
             surf%mom_flux_w(num_v)  = 0.0_wp
             surf%mom_flux_tke(0:1,num_v) = 0.0_wp

             surf%ts(num_v)    = 0.0_wp
             surf%shf(num_v)   = wall_heatflux(component)
!
!--          Set initial value for surface temperature
             surf%pt_surface(num_v) = pt_surface

             IF ( humidity )  THEN
                surf%qs(num_v)   = 0.0_wp
                surf%qsws(num_v) = wall_humidityflux(component)
!
!--             Following wall fluxes are assumed to be zero 
                IF ( surf_bulk_cloud_model .AND. surf_microphysics_morrison)  THEN
                   surf%qcs(num_v) = 0.0_wp
                   surf%ncs(num_v) = 0.0_wp
   
                   surf%qcsws(num_v) = 0.0_wp
                   surf%ncsws(num_v) = 0.0_wp
                ENDIF
                IF ( surf_bulk_cloud_model .AND. surf_microphysics_seifert)  THEN
                   surf%qrs(num_v) = 0.0_wp
                   surf%nrs(num_v) = 0.0_wp
   
                   surf%qrsws(num_v) = 0.0_wp
                   surf%nrsws(num_v) = 0.0_wp
                ENDIF
             ENDIF

             IF ( passive_scalar )  THEN
                surf%ss(num_v)   = 0.0_wp
                surf%ssws(num_v) = wall_scalarflux(component)
             ENDIF

             IF ( air_chemistry )  THEN        
                DO  lsp = 1, nvar
                   surf%css(lsp,num_v)   = 0.0_wp
                   surf%cssws(lsp,num_v) = wall_csflux(lsp,component)
                ENDDO
             ENDIF

!
!--          So far, salinityflux at vertical surfaces is simply zero 
!--          at the moment  
             IF ( ocean_mode )  surf%sasws(num_v) = wall_salinityflux(component)
!
!--          Increment wall indices
             num_v                 = num_v + 1
             num_v_kji             = num_v_kji + 1

          END SUBROUTINE initialize_vertical_surfaces

    END SUBROUTINE init_surfaces

! Description:
! ------------
!> Initialize single surface properties from 2D input arrays
!------------------------------------------------------------------------------!
    SUBROUTINE init_single_surface_properties( var_surf, var_2d,               &
                                               ns, fill_value,                 &
                                               index_space_i,                  &
                                               index_space_j                   &
                                             )
    
       INTEGER(iwp) ::  m  !< running index over surface elements
       INTEGER(iwp) ::  ns !< number of surface elements in var_surf

       INTEGER(iwp), DIMENSION(1:ns) ::  index_space_i !< grid indices along x direction where surface properties should be defined
       INTEGER(iwp), DIMENSION(1:ns) ::  index_space_j !< grid indices along y direction where surface properties should be defined 
       
       REAL(wp) ::  fill_value !< fill value in var_2d
       
       REAL(wp), DIMENSION(1:ns) ::  var_surf !< 1D surface variable that should be initialized
       REAL(wp), DIMENSION(nys:nyn,nxl:nxr) ::  var_2d !< input variable

       DO  m = 1, ns
          IF ( var_2d(index_space_j(m),index_space_i(m)) /= fill_value )  THEN
             var_surf(m) = var_2d(index_space_j(m),index_space_i(m))
          ENDIF
       ENDDO
       
    END SUBROUTINE init_single_surface_properties

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Gathers all surface elements with the same facing (but possibly different 
!> type) onto a surface type, and writes binary data into restart files. 
!------------------------------------------------------------------------------!
    SUBROUTINE surface_wrd_local


       IMPLICIT NONE

       CHARACTER(LEN=1)             ::  dum           !< dummy string to create output-variable name

       INTEGER(iwp)                 ::  i             !< running index x-direction
       INTEGER(iwp)                 ::  j             !< running index y-direction
       INTEGER(iwp)                 ::  l             !< index surface type orientation
       INTEGER(iwp)                 ::  lsp           !< running index chemical species
       INTEGER(iwp)                 ::  m             !< running index for surface elements on individual surface array
       INTEGER(iwp), DIMENSION(0:2) ::  start_index_h !< start index for horizontal surface elements on gathered surface array
       INTEGER(iwp), DIMENSION(0:3) ::  mm            !< running index for surface elements on gathered surface array
       INTEGER(iwp), DIMENSION(0:3) ::  start_index_v !< start index for vertical surface elements on gathered surface array

       TYPE(surf_type), DIMENSION(0:2) ::  surf_h     !< gathered horizontal surfaces, contains all surface types
       TYPE(surf_type), DIMENSION(0:3) ::  surf_v     !< gathered vertical surfaces, contains all surface types

!
!--    Determine total number of horizontal and vertical surface elements before
!--    writing var_list
       CALL surface_last_actions
!
!--    Count number of grid points with same facing and allocate attributes respectively
!--    Horizontal upward facing
       surf_h(0)%ns = ns_h_on_file(0)
       CALL allocate_surface_attributes_h( surf_h(0), nys, nyn, nxl, nxr )
!
!--    Horizontal downward facing
       surf_h(1)%ns = ns_h_on_file(1)
       CALL allocate_surface_attributes_h( surf_h(1), nys, nyn, nxl, nxr )
!
!--    Model top
       surf_h(2)%ns = ns_h_on_file(2)
       CALL allocate_surface_attributes_h_top( surf_h(2), nys, nyn, nxl, nxr )
!
!--    Vertical surfaces
       DO  l = 0, 3
          surf_v(l)%ns = ns_v_on_file(l)
          CALL allocate_surface_attributes_v( surf_v(l),                       &
                                              nys, nyn, nxl, nxr )
       ENDDO
!
!--    In the following, gather data from surfaces elements with the same 
!--    facing (but possibly differt type) on 1 data-type array.
       mm(0:2) = 1
       DO  l = 0, 2
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  m = surf_def_h(l)%start_index(j,i),                        &
                        surf_def_h(l)%end_index(j,i)
                   IF ( ALLOCATED( surf_def_h(l)%us ) )                        &
                      surf_h(l)%us(mm(l))      = surf_def_h(l)%us(m)
                   IF ( ALLOCATED( surf_def_h(l)%ts ) )                        &
                      surf_h(l)%ts(mm(l))      = surf_def_h(l)%ts(m)
                   IF ( ALLOCATED( surf_def_h(l)%qs ) )                        &
                      surf_h(l)%qs(mm(l))      = surf_def_h(l)%qs(m)
                   IF ( ALLOCATED( surf_def_h(l)%ss ) )                        &
                      surf_h(l)%ss(mm(l))      = surf_def_h(l)%ss(m)
                   IF ( ALLOCATED( surf_def_h(l)%qcs ) )                       &
                      surf_h(l)%qcs(mm(l))     = surf_def_h(l)%qcs(m)
                   IF ( ALLOCATED( surf_def_h(l)%ncs ) )                       &
                      surf_h(l)%ncs(mm(l))     = surf_def_h(l)%ncs(m)
                   IF ( ALLOCATED( surf_def_h(l)%qrs ) )                       &
                      surf_h(l)%qrs(mm(l))     = surf_def_h(l)%qrs(m)
                   IF ( ALLOCATED( surf_def_h(l)%nrs ) )                       &
                      surf_h(l)%nrs(mm(l))     = surf_def_h(l)%nrs(m)
                   IF ( ALLOCATED( surf_def_h(l)%ol ) )                        &
                      surf_h(l)%ol(mm(l))      = surf_def_h(l)%ol(m)
                   IF ( ALLOCATED( surf_def_h(l)%rib ) )                       &
                      surf_h(l)%rib(mm(l))     = surf_def_h(l)%rib(m)
                   IF ( ALLOCATED( surf_def_h(l)%pt_surface ) )                &
                      surf_h(l)%pt_surface(mm(l)) = surf_def_h(l)%pt_surface(m)
                   IF ( ALLOCATED( surf_def_h(l)%q_surface ) )                 &
                      surf_h(l)%q_surface(mm(l)) = surf_def_h(l)%q_surface(m) 
                   IF ( ALLOCATED( surf_def_h(l)%vpt_surface ) )               &
                      surf_h(l)%vpt_surface(mm(l)) = surf_def_h(l)%vpt_surface(m)                      
                   IF ( ALLOCATED( surf_def_h(l)%usws ) )                      &
                      surf_h(l)%usws(mm(l))    = surf_def_h(l)%usws(m)
                   IF ( ALLOCATED( surf_def_h(l)%vsws ) )                      &
                      surf_h(l)%vsws(mm(l))    = surf_def_h(l)%vsws(m)
                   IF ( ALLOCATED( surf_def_h(l)%shf ) )                       &
                      surf_h(l)%shf(mm(l))     = surf_def_h(l)%shf(m)
                   IF ( ALLOCATED( surf_def_h(l)%qsws ) )                      &
                      surf_h(l)%qsws(mm(l))    = surf_def_h(l)%qsws(m)
                   IF ( ALLOCATED( surf_def_h(l)%ssws ) )                      &
                      surf_h(l)%ssws(mm(l))    = surf_def_h(l)%ssws(m)
                   IF ( ALLOCATED( surf_def_h(l)%css ) )  THEN 
                      DO  lsp = 1,nvar
                         surf_h(l)%css(lsp,mm(l)) = surf_def_h(l)%css(lsp,m)
                      ENDDO
                   ENDIF
                   IF ( ALLOCATED( surf_def_h(l)%cssws ) )  THEN 
                      DO  lsp = 1,nvar
                         surf_h(l)%cssws(lsp,mm(l)) = surf_def_h(l)%cssws(lsp,m)
                      ENDDO
                   ENDIF
                   IF ( ALLOCATED( surf_def_h(l)%qcsws ) )                     &
                      surf_h(l)%qcsws(mm(l))   = surf_def_h(l)%qcsws(m)
                   IF ( ALLOCATED( surf_def_h(l)%qrsws ) )                     &
                      surf_h(l)%qrsws(mm(l))   = surf_def_h(l)%qrsws(m)
                   IF ( ALLOCATED( surf_def_h(l)%ncsws ) )                     &
                      surf_h(l)%ncsws(mm(l))   = surf_def_h(l)%ncsws(m)
                   IF ( ALLOCATED( surf_def_h(l)%nrsws ) )                     &
                      surf_h(l)%nrsws(mm(l))   = surf_def_h(l)%nrsws(m)
                   IF ( ALLOCATED( surf_def_h(l)%sasws ) )                     &
                      surf_h(l)%sasws(mm(l))   = surf_def_h(l)%sasws(m)
                
                   mm(l) = mm(l) + 1
                ENDDO

                IF ( l == 0 )  THEN
                   DO  m = surf_lsm_h%start_index(j,i),                        &
                           surf_lsm_h%end_index(j,i)
                      IF ( ALLOCATED( surf_lsm_h%us ) )                        &
                         surf_h(0)%us(mm(0))      = surf_lsm_h%us(m)
                      IF ( ALLOCATED( surf_lsm_h%ts ) )                        &
                         surf_h(0)%ts(mm(0))      = surf_lsm_h%ts(m)
                      IF ( ALLOCATED( surf_lsm_h%qs ) )                        &
                         surf_h(0)%qs(mm(0))      = surf_lsm_h%qs(m)
                      IF ( ALLOCATED( surf_lsm_h%ss ) )                        &
                         surf_h(0)%ss(mm(0))      = surf_lsm_h%ss(m)
                      IF ( ALLOCATED( surf_lsm_h%qcs ) )                       &
                         surf_h(0)%qcs(mm(0))     = surf_lsm_h%qcs(m)
                      IF ( ALLOCATED( surf_lsm_h%ncs ) )                       &
                         surf_h(0)%ncs(mm(0))     = surf_lsm_h%ncs(m)
                      IF ( ALLOCATED( surf_lsm_h%qrs ) )                       &
                         surf_h(0)%qrs(mm(0))     = surf_lsm_h%qrs(m)
                      IF ( ALLOCATED( surf_lsm_h%nrs ) )                       &
                         surf_h(0)%nrs(mm(0))     = surf_lsm_h%nrs(m)
                      IF ( ALLOCATED( surf_lsm_h%ol ) )                        &
                         surf_h(0)%ol(mm(0))      = surf_lsm_h%ol(m)
                      IF ( ALLOCATED( surf_lsm_h%rib ) )                       &
                         surf_h(0)%rib(mm(0))     = surf_lsm_h%rib(m)
                      IF ( ALLOCATED( surf_lsm_h%pt_surface ) )                &
                         surf_h(l)%pt_surface(mm(l)) = surf_lsm_h%pt_surface(m)
                      IF ( ALLOCATED( surf_def_h(l)%q_surface ) )              &
                         surf_h(l)%q_surface(mm(l)) = surf_lsm_h%q_surface(m)
                      IF ( ALLOCATED( surf_def_h(l)%vpt_surface ) )            &
                         surf_h(l)%vpt_surface(mm(l)) = surf_lsm_h%vpt_surface(m)
                      IF ( ALLOCATED( surf_lsm_h%usws ) )                      &
                         surf_h(0)%usws(mm(0))    = surf_lsm_h%usws(m)
                      IF ( ALLOCATED( surf_lsm_h%vsws ) )                      &
                         surf_h(0)%vsws(mm(0))    = surf_lsm_h%vsws(m)
                      IF ( ALLOCATED( surf_lsm_h%shf ) )                       &
                         surf_h(0)%shf(mm(0))     = surf_lsm_h%shf(m)
                      IF ( ALLOCATED( surf_lsm_h%qsws ) )                      &
                         surf_h(0)%qsws(mm(0))    = surf_lsm_h%qsws(m)
                      IF ( ALLOCATED( surf_lsm_h%ssws ) )                      &
                         surf_h(0)%ssws(mm(0))    = surf_lsm_h%ssws(m)
                      IF ( ALLOCATED( surf_lsm_h%css ) )  THEN                  
                         DO  lsp = 1, nvar
                            surf_h(0)%css(lsp,mm(0)) = surf_lsm_h%css(lsp,m)
                         ENDDO
                      ENDIF
                      IF ( ALLOCATED( surf_lsm_h%cssws ) )  THEN
                         DO  lsp = 1, nvar
                            surf_h(0)%cssws(lsp,mm(0)) = surf_lsm_h%cssws(lsp,m)
                         ENDDO 
                      ENDIF
                      IF ( ALLOCATED( surf_lsm_h%qcsws ) )                     &
                         surf_h(0)%qcsws(mm(0))   = surf_lsm_h%qcsws(m)
                      IF ( ALLOCATED( surf_lsm_h%qrsws ) )                     &
                         surf_h(0)%qrsws(mm(0))   = surf_lsm_h%qrsws(m)
                      IF ( ALLOCATED( surf_lsm_h%ncsws ) )                     &
                         surf_h(0)%ncsws(mm(0))   = surf_lsm_h%ncsws(m)
                      IF ( ALLOCATED( surf_lsm_h%nrsws ) )                     &
                         surf_h(0)%nrsws(mm(0))   = surf_lsm_h%nrsws(m)
                      IF ( ALLOCATED( surf_lsm_h%sasws ) )                     &
                        surf_h(0)%sasws(mm(0))   = surf_lsm_h%sasws(m)
                
                      mm(0) = mm(0) + 1
             
                   ENDDO

                   DO  m = surf_usm_h%start_index(j,i),                        &
                           surf_usm_h%end_index(j,i)
                      IF ( ALLOCATED( surf_usm_h%us ) )                        &
                         surf_h(0)%us(mm(0))      = surf_usm_h%us(m)
                      IF ( ALLOCATED( surf_usm_h%ts ) )                        &
                         surf_h(0)%ts(mm(0))      = surf_usm_h%ts(m)
                      IF ( ALLOCATED( surf_usm_h%qs ) )                        &
                         surf_h(0)%qs(mm(0))      = surf_usm_h%qs(m)
                      IF ( ALLOCATED( surf_usm_h%ss ) )                        &
                         surf_h(0)%ss(mm(0))      = surf_usm_h%ss(m)
                      IF ( ALLOCATED( surf_usm_h%qcs ) )                       &
                         surf_h(0)%qcs(mm(0))     = surf_usm_h%qcs(m)
                      IF ( ALLOCATED( surf_usm_h%ncs ) )                       &
                         surf_h(0)%ncs(mm(0))     = surf_usm_h%ncs(m)
                      IF ( ALLOCATED( surf_usm_h%qrs ) )                       &
                         surf_h(0)%qrs(mm(0))     = surf_usm_h%qrs(m)
                      IF ( ALLOCATED( surf_usm_h%nrs ) )                       &
                         surf_h(0)%nrs(mm(0))     = surf_usm_h%nrs(m)
                      IF ( ALLOCATED( surf_usm_h%ol ) )                        &
                         surf_h(0)%ol(mm(0))      = surf_usm_h%ol(m)
                      IF ( ALLOCATED( surf_usm_h%rib ) )                       &
                         surf_h(0)%rib(mm(0))     = surf_usm_h%rib(m)
                      IF ( ALLOCATED( surf_usm_h%pt_surface ) )                &
                         surf_h(l)%pt_surface(mm(l)) = surf_usm_h%pt_surface(m)
                       IF ( ALLOCATED( surf_usm_h%q_surface ) )                &
                         surf_h(l)%q_surface(mm(l)) = surf_usm_h%q_surface(m)
                      IF ( ALLOCATED( surf_usm_h%vpt_surface ) )               &
                         surf_h(l)%vpt_surface(mm(l)) = surf_usm_h%vpt_surface(m)
                      IF ( ALLOCATED( surf_usm_h%usws ) )                      &
                         surf_h(0)%usws(mm(0))    = surf_usm_h%usws(m)
                      IF ( ALLOCATED( surf_usm_h%vsws ) )                      &
                         surf_h(0)%vsws(mm(0))    = surf_usm_h%vsws(m)
                      IF ( ALLOCATED( surf_usm_h%shf ) )                       &
                         surf_h(0)%shf(mm(0))     = surf_usm_h%shf(m)
                      IF ( ALLOCATED( surf_usm_h%qsws ) )                      &
                         surf_h(0)%qsws(mm(0))    = surf_usm_h%qsws(m)
                      IF ( ALLOCATED( surf_usm_h%ssws ) )                      &
                         surf_h(0)%ssws(mm(0))    = surf_usm_h%ssws(m)
                      IF ( ALLOCATED( surf_usm_h%css ) )  THEN             
                         DO lsp = 1, nvar
                            surf_h(0)%css(lsp,mm(0)) = surf_usm_h%css(lsp,m)
                         ENDDO
                      ENDIF
                      IF ( ALLOCATED( surf_usm_h%cssws ) )  THEN             
                         DO lsp = 1, nvar
                            surf_h(0)%cssws(lsp,mm(0)) = surf_usm_h%cssws(lsp,m)
                         ENDDO
                      ENDIF
                      IF ( ALLOCATED( surf_usm_h%qcsws ) )                     &
                         surf_h(0)%qcsws(mm(0))   = surf_usm_h%qcsws(m)
                      IF ( ALLOCATED( surf_usm_h%qrsws ) )                     &
                         surf_h(0)%qrsws(mm(0))   = surf_usm_h%qrsws(m)
                      IF ( ALLOCATED( surf_usm_h%ncsws ) )                     &
                         surf_h(0)%ncsws(mm(0))   = surf_usm_h%ncsws(m)
                      IF ( ALLOCATED( surf_usm_h%nrsws ) )                     &
                         surf_h(0)%nrsws(mm(0))   = surf_usm_h%nrsws(m)
                      IF ( ALLOCATED( surf_usm_h%sasws ) )                     &
                        surf_h(0)%sasws(mm(0))   = surf_usm_h%sasws(m)
                
                      mm(0) = mm(0) + 1
             
                   ENDDO


                ENDIF

             ENDDO

          ENDDO
!
!--       Recalculate start- and end indices for gathered surface type.
          start_index_h(l) = 1                                        
          DO  i = nxl, nxr
             DO  j = nys, nyn

                surf_h(l)%start_index(j,i) = start_index_h(l)
                surf_h(l)%end_index(j,i)   = surf_h(l)%start_index(j,i) - 1

                DO  m = surf_def_h(l)%start_index(j,i),                        &
                        surf_def_h(l)%end_index(j,i)
                   surf_h(l)%end_index(j,i) = surf_h(l)%end_index(j,i) + 1
                ENDDO
                IF ( l == 0 )  THEN
                   DO  m = surf_lsm_h%start_index(j,i),                        &
                           surf_lsm_h%end_index(j,i)
                      surf_h(l)%end_index(j,i) = surf_h(l)%end_index(j,i) + 1
                   ENDDO
                   DO  m = surf_usm_h%start_index(j,i),                        &
                           surf_usm_h%end_index(j,i)
                      surf_h(l)%end_index(j,i) = surf_h(l)%end_index(j,i) + 1
                   ENDDO
                ENDIF

                start_index_h(l) = surf_h(l)%end_index(j,i) + 1

             ENDDO
          ENDDO
       ENDDO
!
!--    Treat vertically orientated surface. Again, gather data from different
!--    surfaces types but identical orientation (e.g. northward-facing) onto
!--    one surface type which is output afterwards. 
       mm(0:3) = 1
       DO  l = 0, 3
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  m = surf_def_v(l)%start_index(j,i),                        &
                        surf_def_v(l)%end_index(j,i)
                   IF ( ALLOCATED( surf_def_v(l)%us ) )                        &
                      surf_v(l)%us(mm(l))      = surf_def_v(l)%us(m)
                   IF ( ALLOCATED( surf_def_v(l)%ts ) )                        &
                      surf_v(l)%ts(mm(l))      = surf_def_v(l)%ts(m)
                   IF ( ALLOCATED( surf_def_v(l)%qs ) )                        &
                      surf_v(l)%qs(mm(l))      = surf_def_v(l)%qs(m)
                   IF ( ALLOCATED( surf_def_v(l)%ss ) )                        &
                      surf_v(l)%ss(mm(l))      = surf_def_v(l)%ss(m)
                   IF ( ALLOCATED( surf_def_v(l)%qcs ) )                       &
                      surf_v(l)%qcs(mm(l))     = surf_def_v(l)%qcs(m)
                   IF ( ALLOCATED( surf_def_v(l)%ncs ) )                       &
                      surf_v(l)%ncs(mm(l))     = surf_def_v(l)%ncs(m)
                   IF ( ALLOCATED( surf_def_v(l)%qrs ) )                       &
                      surf_v(l)%qrs(mm(l))     = surf_def_v(l)%qrs(m)
                   IF ( ALLOCATED( surf_def_v(l)%nrs ) )                       &
                      surf_v(l)%nrs(mm(l))     = surf_def_v(l)%nrs(m)
                   IF ( ALLOCATED( surf_def_v(l)%ol ) )                        &
                      surf_v(l)%ol(mm(l))      = surf_def_v(l)%ol(m)
                   IF ( ALLOCATED( surf_def_v(l)%rib ) )                       &
                      surf_v(l)%rib(mm(l))     = surf_def_v(l)%rib(m)
                   IF ( ALLOCATED( surf_def_v(l)%pt_surface ) )                &
                      surf_v(l)%pt_surface(mm(l)) = surf_def_v(l)%pt_surface(m)
                   IF ( ALLOCATED( surf_def_v(l)%q_surface ) )                 &
                      surf_v(l)%q_surface(mm(l)) = surf_def_v(l)%q_surface(m)
                   IF ( ALLOCATED( surf_def_v(l)%vpt_surface ) )               &
                      surf_v(l)%vpt_surface(mm(l)) = surf_def_v(l)%vpt_surface(m)
                   IF ( ALLOCATED( surf_def_v(l)%shf ) )                       &
                      surf_v(l)%shf(mm(l))     = surf_def_v(l)%shf(m)
                   IF ( ALLOCATED( surf_def_v(l)%qsws ) )                      &
                      surf_v(l)%qsws(mm(l))    = surf_def_v(l)%qsws(m)
                   IF ( ALLOCATED( surf_def_v(l)%ssws ) )                      &
                      surf_v(l)%ssws(mm(l))    = surf_def_v(l)%ssws(m)
                   IF ( ALLOCATED( surf_def_v(l)%css ) )  THEN               
                      DO  lsp = 1, nvar
                         surf_v(l)%css(lsp,mm(l)) = surf_def_v(l)%css(lsp,m)
                      ENDDO
                   ENDIF
                   IF ( ALLOCATED( surf_def_v(l)%cssws ) )  THEN               
                      DO  lsp = 1, nvar
                         surf_v(l)%cssws(lsp,mm(l)) = surf_def_v(l)%cssws(lsp,m)
                      ENDDO
                   ENDIF
                   IF ( ALLOCATED( surf_def_v(l)%qcsws ) )                     &
                      surf_v(l)%qcsws(mm(l))   = surf_def_v(l)%qcsws(m)
                   IF ( ALLOCATED( surf_def_v(l)%qrsws ) )                     &
                      surf_v(l)%qrsws(mm(l))   = surf_def_v(l)%qrsws(m)
                   IF ( ALLOCATED( surf_def_v(l)%ncsws ) )                     &
                      surf_v(l)%ncsws(mm(l))   = surf_def_v(l)%ncsws(m)
                   IF ( ALLOCATED( surf_def_v(l)%nrsws ) )                     &
                      surf_v(l)%nrsws(mm(l))   = surf_def_v(l)%nrsws(m)
                   IF ( ALLOCATED( surf_def_v(l)%sasws ) )                     &
                      surf_v(l)%sasws(mm(l))   = surf_def_v(l)%sasws(m)
                   IF ( ALLOCATED( surf_def_v(l)%mom_flux_uv) )                &
                      surf_v(l)%mom_flux_uv(mm(l))  = surf_def_v(l)%mom_flux_uv(m)
                   IF ( ALLOCATED( surf_def_v(l)%mom_flux_w) )                 &
                      surf_v(l)%mom_flux_w(mm(l))   = surf_def_v(l)%mom_flux_w(m)
                   IF ( ALLOCATED( surf_def_v(l)%mom_flux_tke) )               &
                      surf_v(l)%mom_flux_tke(0:1,mm(l)) = surf_def_v(l)%mom_flux_tke(0:1,m)
                
                   mm(l) = mm(l) + 1
                ENDDO

                DO  m = surf_lsm_v(l)%start_index(j,i),                        &
                        surf_lsm_v(l)%end_index(j,i)
                   IF ( ALLOCATED( surf_lsm_v(l)%us ) )                        &
                      surf_v(l)%us(mm(l))      = surf_lsm_v(l)%us(m)
                   IF ( ALLOCATED( surf_lsm_v(l)%ts ) )                        &
                      surf_v(l)%ts(mm(l))      = surf_lsm_v(l)%ts(m)
                   IF ( ALLOCATED( surf_lsm_v(l)%qs ) )                        &
                      surf_v(l)%qs(mm(l))      = surf_lsm_v(l)%qs(m)
                   IF ( ALLOCATED( surf_lsm_v(l)%ss ) )                        &
                      surf_v(l)%ss(mm(l))      = surf_lsm_v(l)%ss(m)
                   IF ( ALLOCATED( surf_lsm_v(l)%qcs ) )                       &
                      surf_v(l)%qcs(mm(l))     = surf_lsm_v(l)%qcs(m)
                   IF ( ALLOCATED( surf_lsm_v(l)%ncs ) )                       &
                      surf_v(l)%ncs(mm(l))     = surf_lsm_v(l)%ncs(m)
                   IF ( ALLOCATED( surf_lsm_v(l)%qrs ) )                       &
                      surf_v(l)%qrs(mm(l))     = surf_lsm_v(l)%qrs(m)
                   IF ( ALLOCATED( surf_lsm_v(l)%nrs ) )                       &
                      surf_v(l)%nrs(mm(l))     = surf_lsm_v(l)%nrs(m)
                   IF ( ALLOCATED( surf_lsm_v(l)%ol ) )                        &
                      surf_v(l)%ol(mm(l))      = surf_lsm_v(l)%ol(m)
                   IF ( ALLOCATED( surf_lsm_v(l)%rib ) )                       &
                      surf_v(l)%rib(mm(l))     = surf_lsm_v(l)%rib(m)
                   IF ( ALLOCATED( surf_lsm_v(l)%pt_surface ) )                &
                      surf_v(l)%pt_surface(mm(l)) = surf_lsm_v(l)%pt_surface(m)
                   IF ( ALLOCATED( surf_lsm_v(l)%q_surface ) )                 &
                      surf_v(l)%q_surface(mm(l)) = surf_lsm_v(l)%q_surface(m)
                   IF ( ALLOCATED( surf_lsm_v(l)%vpt_surface ) )               &
                      surf_v(l)%vpt_surface(mm(l)) = surf_lsm_v(l)%vpt_surface(m)
                   IF ( ALLOCATED( surf_lsm_v(l)%usws ) )                      &
                      surf_v(l)%usws(mm(l))    = surf_lsm_v(l)%usws(m)
                   IF ( ALLOCATED( surf_lsm_v(l)%vsws ) )                      &
                      surf_v(l)%vsws(mm(l))    = surf_lsm_v(l)%vsws(m)
                   IF ( ALLOCATED( surf_lsm_v(l)%shf ) )                       &
                      surf_v(l)%shf(mm(l))     = surf_lsm_v(l)%shf(m)
                   IF ( ALLOCATED( surf_lsm_v(l)%qsws ) )                      &
                      surf_v(l)%qsws(mm(l))    = surf_lsm_v(l)%qsws(m)
                   IF ( ALLOCATED( surf_lsm_v(l)%ssws ) )                      &
                      surf_v(l)%ssws(mm(l))    = surf_lsm_v(l)%ssws(m)
                   IF ( ALLOCATED( surf_lsm_v(l)%css ) )  THEN              
                      DO  lsp = 1, nvar
                         surf_v(l)%css(lsp,mm(l)) = surf_lsm_v(l)%css(lsp,m)
                      ENDDO
                   ENDIF
                   IF ( ALLOCATED( surf_lsm_v(l)%cssws ) )  THEN              
                      DO  lsp = 1, nvar
                         surf_v(l)%cssws(lsp,mm(l)) = surf_lsm_v(l)%cssws(lsp,m)
                      ENDDO
                   ENDIF
                   IF ( ALLOCATED( surf_lsm_v(l)%qcsws ) )                     &
                      surf_v(l)%qcsws(mm(l))   = surf_lsm_v(l)%qcsws(m)
                   IF ( ALLOCATED( surf_lsm_v(l)%qrsws ) )                     &
                      surf_v(l)%qrsws(mm(l))   = surf_lsm_v(l)%qrsws(m)
                   IF ( ALLOCATED( surf_lsm_v(l)%ncsws ) )                     &
                      surf_v(l)%ncsws(mm(l))   = surf_lsm_v(l)%ncsws(m)
                   IF ( ALLOCATED( surf_lsm_v(l)%nrsws ) )                     &
                      surf_v(l)%nrsws(mm(l))   = surf_lsm_v(l)%nrsws(m)
                   IF ( ALLOCATED( surf_lsm_v(l)%sasws ) )                     &
                      surf_v(l)%sasws(mm(l))   = surf_lsm_v(l)%sasws(m)
                   IF ( ALLOCATED( surf_lsm_v(l)%mom_flux_uv) )                &
                      surf_v(l)%mom_flux_uv(mm(l))  = surf_lsm_v(l)%mom_flux_uv(m)
                   IF ( ALLOCATED( surf_lsm_v(l)%mom_flux_w) )                 &
                      surf_v(l)%mom_flux_w(mm(l))   = surf_lsm_v(l)%mom_flux_w(m)
                   IF ( ALLOCATED( surf_lsm_v(l)%mom_flux_tke) )               &
                      surf_v(l)%mom_flux_tke(0:1,mm(l)) = surf_lsm_v(l)%mom_flux_tke(0:1,m)
                
                   mm(l) = mm(l) + 1
                ENDDO

                DO  m = surf_usm_v(l)%start_index(j,i),                        &
                        surf_usm_v(l)%end_index(j,i)
                   IF ( ALLOCATED( surf_usm_v(l)%us ) )                        &
                      surf_v(l)%us(mm(l))      = surf_usm_v(l)%us(m)
                   IF ( ALLOCATED( surf_usm_v(l)%ts ) )                        &
                      surf_v(l)%ts(mm(l))      = surf_usm_v(l)%ts(m)
                   IF ( ALLOCATED( surf_usm_v(l)%qs ) )                        &
                      surf_v(l)%qs(mm(l))      = surf_usm_v(l)%qs(m)
                   IF ( ALLOCATED( surf_usm_v(l)%ss ) )                        &
                      surf_v(l)%ss(mm(l))      = surf_usm_v(l)%ss(m)
                   IF ( ALLOCATED( surf_usm_v(l)%qcs ) )                       &
                      surf_v(l)%qcs(mm(l))     = surf_usm_v(l)%qcs(m)
                   IF ( ALLOCATED( surf_usm_v(l)%ncs ) )                       &
                      surf_v(l)%ncs(mm(l))     = surf_usm_v(l)%ncs(m)
                   IF ( ALLOCATED( surf_usm_v(l)%qrs ) )                       &
                      surf_v(l)%qrs(mm(l))     = surf_usm_v(l)%qrs(m)
                   IF ( ALLOCATED( surf_usm_v(l)%nrs ) )                       &
                      surf_v(l)%nrs(mm(l))     = surf_usm_v(l)%nrs(m)
                   IF ( ALLOCATED( surf_usm_v(l)%ol ) )                        &
                      surf_v(l)%ol(mm(l))      = surf_usm_v(l)%ol(m)
                   IF ( ALLOCATED( surf_usm_v(l)%rib ) )                       &
                      surf_v(l)%rib(mm(l))     = surf_usm_v(l)%rib(m)
                   IF ( ALLOCATED( surf_usm_v(l)%pt_surface ) )                &
                      surf_v(l)%pt_surface(mm(l)) = surf_usm_v(l)%pt_surface(m)
                   IF ( ALLOCATED( surf_usm_v(l)%q_surface ) )                 &
                      surf_v(l)%q_surface(mm(l)) = surf_usm_v(l)%q_surface(m)
                   IF ( ALLOCATED( surf_usm_v(l)%vpt_surface ) )               &
                      surf_v(l)%vpt_surface(mm(l)) = surf_usm_v(l)%vpt_surface(m)
                   IF ( ALLOCATED( surf_usm_v(l)%usws ) )                      &
                      surf_v(l)%usws(mm(l))    = surf_usm_v(l)%usws(m)
                   IF ( ALLOCATED( surf_usm_v(l)%vsws ) )                      &
                      surf_v(l)%vsws(mm(l))    = surf_usm_v(l)%vsws(m)
                   IF ( ALLOCATED( surf_usm_v(l)%shf ) )                       &
                      surf_v(l)%shf(mm(l))     = surf_usm_v(l)%shf(m)
                   IF ( ALLOCATED( surf_usm_v(l)%qsws ) )                      &
                      surf_v(l)%qsws(mm(l))    = surf_usm_v(l)%qsws(m)
                   IF ( ALLOCATED( surf_usm_v(l)%ssws ) )                      &
                      surf_v(l)%ssws(mm(l))    = surf_usm_v(l)%ssws(m)
                   IF ( ALLOCATED( surf_usm_v(l)%css ) )  THEN              
                      DO  lsp = 1, nvar
                         surf_v(l)%css(lsp,mm(l)) = surf_usm_v(l)%css(lsp,m)
                      ENDDO
                   ENDIF
                   IF ( ALLOCATED( surf_usm_v(l)%cssws ) )  THEN              
                      DO  lsp = 1, nvar
                         surf_v(l)%cssws(lsp,mm(l)) = surf_usm_v(l)%cssws(lsp,m)
                      ENDDO
                   ENDIF
                   IF ( ALLOCATED( surf_usm_v(l)%qcsws ) )                     &
                      surf_v(l)%qcsws(mm(l))   = surf_usm_v(l)%qcsws(m)
                   IF ( ALLOCATED( surf_usm_v(l)%qrsws ) )                     &
                      surf_v(l)%qrsws(mm(l))   = surf_usm_v(l)%qrsws(m)
                   IF ( ALLOCATED( surf_usm_v(l)%ncsws ) )                     &
                      surf_v(l)%ncsws(mm(l))   = surf_usm_v(l)%ncsws(m)
                   IF ( ALLOCATED( surf_usm_v(l)%nrsws ) )                     &
                      surf_v(l)%nrsws(mm(l))   = surf_usm_v(l)%nrsws(m)
                   IF ( ALLOCATED( surf_usm_v(l)%sasws ) )                     &
                      surf_v(l)%sasws(mm(l))   = surf_usm_v(l)%sasws(m)
                   IF ( ALLOCATED( surf_usm_v(l)%mom_flux_uv) )                &
                      surf_v(l)%mom_flux_uv(mm(l))  = surf_usm_v(l)%mom_flux_uv(m)
                   IF ( ALLOCATED( surf_usm_v(l)%mom_flux_w) )                 &
                      surf_v(l)%mom_flux_w(mm(l))   = surf_usm_v(l)%mom_flux_w(m)
                   IF ( ALLOCATED( surf_usm_v(l)%mom_flux_tke) )               &
                      surf_v(l)%mom_flux_tke(0:1,mm(l)) = surf_usm_v(l)%mom_flux_tke(0:1,m)
                
                   mm(l) = mm(l) + 1
                ENDDO
             
             ENDDO
          ENDDO
!
!--       Recalculate start- and end-indices for gathered surface type
          start_index_v(l) = 1                                        
          DO  i = nxl, nxr
             DO  j = nys, nyn

                surf_v(l)%start_index(j,i) = start_index_v(l)
                surf_v(l)%end_index(j,i)   = surf_v(l)%start_index(j,i) -1

                DO  m = surf_def_v(l)%start_index(j,i),                        &
                        surf_def_v(l)%end_index(j,i)
                   surf_v(l)%end_index(j,i) = surf_v(l)%end_index(j,i) + 1
                ENDDO
                DO  m = surf_lsm_v(l)%start_index(j,i),                        &
                        surf_lsm_v(l)%end_index(j,i)
                   surf_v(l)%end_index(j,i) = surf_v(l)%end_index(j,i) + 1
                ENDDO
                DO  m = surf_usm_v(l)%start_index(j,i),                        &
                        surf_usm_v(l)%end_index(j,i)
                   surf_v(l)%end_index(j,i) = surf_v(l)%end_index(j,i) + 1
                ENDDO

                start_index_v(l) = surf_v(l)%end_index(j,i) + 1
             ENDDO
          ENDDO

       ENDDO
!
!--    Output strings for the total number of upward / downward-facing surfaces 
!--    on subdomain.
       CALL wrd_write_string( 'ns_h_on_file' )
       WRITE ( 14 )  ns_h_on_file
!
!--    Output strings for the total number of north/south/east/westward-facing surfaces 
!--    on subdomain.
       CALL wrd_write_string( 'ns_v_on_file' )
       WRITE ( 14 )  ns_v_on_file

!
!--    Write required restart data.
!--    Start with horizontal surfaces (upward-, downward-facing, and model top).
!--    Always start with %start_index followed by %end_index
       DO  l = 0, 2
          WRITE( dum, '(I1)')  l

          CALL wrd_write_string( 'surf_h(' // dum // ')%start_index' )
          WRITE ( 14 )  surf_h(l)%start_index

          CALL wrd_write_string( 'surf_h(' // dum // ')%end_index' )
          WRITE ( 14 )  surf_h(l)%end_index

          IF ( ALLOCATED ( surf_h(l)%us ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%us' ) 
             WRITE ( 14 )  surf_h(l)%us
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%ts ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%ts' ) 
             WRITE ( 14 )  surf_h(l)%ts
          ENDIF
          
          IF ( ALLOCATED ( surf_h(l)%qs ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%qs' )  
             WRITE ( 14 )  surf_h(l)%qs
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%ss ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%ss' )  
             WRITE ( 14 )  surf_h(l)%ss
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%qcs ) )  THEN  
             CALL wrd_write_string( 'surf_h(' // dum // ')%qcs' )
             WRITE ( 14 )  surf_h(l)%qcs
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%ncs ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%ncs' ) 
             WRITE ( 14 )  surf_h(l)%ncs
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%qrs ) )  THEN  
             CALL wrd_write_string( 'surf_h(' // dum // ')%qrs' )
             WRITE ( 14 )  surf_h(l)%qrs
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%nrs ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%nrs' )  
             WRITE ( 14 )  surf_h(l)%nrs
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%ol ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%ol' )  
             WRITE ( 14 )  surf_h(l)%ol
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%rib ) )  THEN
            CALL wrd_write_string( 'surf_h(' // dum // ')%rib' ) 
             WRITE ( 14 )  surf_h(l)%rib
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%pt_surface ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%pt_surface' ) 
             WRITE ( 14 )  surf_h(l)%pt_surface
          ENDIF
          
          IF ( ALLOCATED ( surf_h(l)%q_surface ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%q_surface' ) 
             WRITE ( 14 )  surf_h(l)%q_surface
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%vpt_surface ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%vpt_surface' ) 
             WRITE ( 14 )  surf_h(l)%vpt_surface
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%usws ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%usws' ) 
             WRITE ( 14 )  surf_h(l)%usws
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%vsws ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%vsws' )  
             WRITE ( 14 )  surf_h(l)%vsws
          ENDIF
          
          IF ( ALLOCATED ( surf_h(l)%shf ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%shf' ) 
             WRITE ( 14 )  surf_h(l)%shf
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%qsws ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%qsws' )  
             WRITE ( 14 )  surf_h(l)%qsws
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%ssws ) )  THEN 
             CALL wrd_write_string( 'surf_h(' // dum // ')%ssws' )  
             WRITE ( 14 )  surf_h(l)%ssws
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%css ) )  THEN 
             CALL wrd_write_string( 'surf_h(' // dum // ')%css' )
             WRITE ( 14 )  surf_h(l)%css
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%cssws ) )  THEN 
             CALL wrd_write_string( 'surf_h(' // dum // ')%cssws' )
             WRITE ( 14 )  surf_h(l)%cssws
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%qcsws ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%qcsws' )  
             WRITE ( 14 )  surf_h(l)%qcsws
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%ncsws ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%ncsws' )  
             WRITE ( 14 )  surf_h(l)%ncsws
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%qrsws ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%qrsws' )  
             WRITE ( 14 )  surf_h(l)%qrsws
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%nrsws ) )  THEN
             CALL wrd_write_string( 'surf_h(' // dum // ')%nrsws' )  
             WRITE ( 14 )  surf_h(l)%nrsws
          ENDIF

          IF ( ALLOCATED ( surf_h(l)%sasws ) )  THEN 
             CALL wrd_write_string( 'surf_h(' // dum // ')%sasws' ) 
             WRITE ( 14 )  surf_h(l)%sasws
          ENDIF
  
       ENDDO
!
!--    Write vertical surfaces.
!--    Always start with %start_index followed by %end_index.
       DO  l = 0, 3
          WRITE( dum, '(I1)')  l

          CALL wrd_write_string( 'surf_v(' // dum // ')%start_index' )
          WRITE ( 14 )  surf_v(l)%start_index

          CALL wrd_write_string( 'surf_v(' // dum // ')%end_index' )
          WRITE ( 14 )   surf_v(l)%end_index

          IF ( ALLOCATED ( surf_v(l)%us ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%us' )  
             WRITE ( 14 )  surf_v(l)%us
          ENDIF  

          IF ( ALLOCATED ( surf_v(l)%ts ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%ts' ) 
             WRITE ( 14 )  surf_v(l)%ts
          ENDIF

          IF ( ALLOCATED ( surf_v(l)%qs ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%qs' )  
             WRITE ( 14 )  surf_v(l)%qs
          ENDIF  

          IF ( ALLOCATED ( surf_v(l)%ss ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%ss' ) 
             WRITE ( 14 )  surf_v(l)%ss
          ENDIF  
          
          IF ( ALLOCATED ( surf_v(l)%qcs ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%qcs' ) 
             WRITE ( 14 )  surf_v(l)%qcs
          ENDIF

          IF ( ALLOCATED ( surf_v(l)%ncs ) )  THEN 
             CALL wrd_write_string( 'surf_v(' // dum // ')%ncs' )
             WRITE ( 14 )  surf_v(l)%ncs
          ENDIF

          IF ( ALLOCATED ( surf_v(l)%qrs ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%qrs' )  
             WRITE ( 14 )  surf_v(l)%qrs
          ENDIF

          IF ( ALLOCATED ( surf_v(l)%nrs ) )  THEN 
             CALL wrd_write_string( 'surf_v(' // dum // ')%nrs' ) 
             WRITE ( 14 )  surf_v(l)%nrs
          ENDIF  

          IF ( ALLOCATED ( surf_v(l)%ol ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%ol' )  
             WRITE ( 14 )  surf_v(l)%ol
          ENDIF

          IF ( ALLOCATED ( surf_v(l)%rib ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%rib' )
             WRITE ( 14 )  surf_v(l)%rib
          ENDIF

          IF ( ALLOCATED ( surf_v(l)%pt_surface ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%pt_surface' )
             WRITE ( 14 )  surf_v(l)%pt_surface
          ENDIF 
          
          IF ( ALLOCATED ( surf_v(l)%q_surface ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%q_surface' )
             WRITE ( 14 )  surf_v(l)%q_surface
          ENDIF

          IF ( ALLOCATED ( surf_v(l)%vpt_surface ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%vpt_surface' )
             WRITE ( 14 )  surf_v(l)%vpt_surface
          ENDIF 

          IF ( ALLOCATED ( surf_v(l)%shf ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%shf' )  
             WRITE ( 14 )  surf_v(l)%shf
          ENDIF

          IF ( ALLOCATED ( surf_v(l)%qsws ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%qsws' )  
             WRITE ( 14 )  surf_v(l)%qsws
          ENDIF

          IF ( ALLOCATED ( surf_v(l)%ssws ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%ssws' )  
             WRITE ( 14 )  surf_v(l)%ssws
          ENDIF

          IF ( ALLOCATED ( surf_v(l)%css ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%css' )  
             WRITE ( 14 )  surf_v(l)%css
          ENDIF

          IF ( ALLOCATED ( surf_v(l)%cssws ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%cssws' ) 
             WRITE ( 14 )  surf_v(l)%cssws
          ENDIF  

          IF ( ALLOCATED ( surf_v(l)%qcsws ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%qcsws' ) 
             WRITE ( 14 )  surf_v(l)%qcsws
          ENDIF  

          IF ( ALLOCATED ( surf_v(l)%ncsws ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%ncsws' ) 
             WRITE ( 14 )  surf_v(l)%ncsws
          ENDIF

          IF ( ALLOCATED ( surf_v(l)%qrsws ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%qrsws' )  
             WRITE ( 14 )  surf_v(l)%qrsws
          ENDIF  

          IF ( ALLOCATED ( surf_v(l)%nrsws ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%nrsws' )  
             WRITE ( 14 )  surf_v(l)%nrsws
          ENDIF 

          IF ( ALLOCATED ( surf_v(l)%sasws ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%sasws' ) 
             WRITE ( 14 )  surf_v(l)%sasws
          ENDIF

          IF ( ALLOCATED ( surf_v(l)%mom_flux_uv ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%mom_uv' )  
             WRITE ( 14 )  surf_v(l)%mom_flux_uv
          ENDIF

          IF ( ALLOCATED ( surf_v(l)%mom_flux_w ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%mom_w' ) 
             WRITE ( 14 )  surf_v(l)%mom_flux_w
          ENDIF

          IF ( ALLOCATED ( surf_v(l)%mom_flux_tke ) )  THEN
             CALL wrd_write_string( 'surf_v(' // dum // ')%mom_tke' )  
             WRITE ( 14 )  surf_v(l)%mom_flux_tke
          ENDIF 
         
       ENDDO


    END SUBROUTINE surface_wrd_local


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads surface-related restart data. Please note, restart data for a certain
!> surface orientation (e.g. horizontal upward-facing) is stored in one 
!> array, even if surface elements may belong to different surface types 
!> natural or urban for example). Surface elements are redistributed into its
!> respective surface types within this routine. This allows e.g. changing the 
!> surface type after reading the restart data, which might be required in case 
!> of cyclic_fill mode. 
!------------------------------------------------------------------------------!
    SUBROUTINE surface_rrd_local( kk, nxlf, nxlc, nxl_on_file, nxrf,           &
                                  nxr_on_file, nynf, nyn_on_file, nysf,        &
                                  nysc, nys_on_file, found )


       IMPLICIT NONE

       INTEGER(iwp)       ::  i           !< running index along x-direction, refers to former domain size
       INTEGER(iwp)       ::  ic          !< running index along x-direction, refers to current domain size
       INTEGER(iwp)       ::  j           !< running index along y-direction, refers to former domain size
       INTEGER(iwp)       ::  jc          !< running index along y-direction, refers to former domain size
       INTEGER(iwp)       ::  m           !< running index for surface elements, refers to gathered array encompassing all surface types
       INTEGER(iwp)       ::  mm          !< running index for surface elements, refers to individual surface types
       INTEGER(iwp)       ::  kk          !< running index over previous input files covering current local domain
       INTEGER(iwp)       ::  nxlc        !< index of left boundary on current subdomain
       INTEGER(iwp)       ::  nxlf        !< index of left boundary on former subdomain 
       INTEGER(iwp)       ::  nxl_on_file !< index of left boundary on former local domain 
       INTEGER(iwp)       ::  nxrf        !< index of right boundary on former subdomain
       INTEGER(iwp)       ::  nxr_on_file !< index of right boundary on former local domain  
       INTEGER(iwp)       ::  nynf        !< index of north boundary on former subdomain
       INTEGER(iwp)       ::  nyn_on_file !< index of norht boundary on former local domain  
       INTEGER(iwp)       ::  nysc        !< index of south boundary on current subdomain 
       INTEGER(iwp)       ::  nysf        !< index of south boundary on former subdomain
       INTEGER(iwp)       ::  nys_on_file !< index of south boundary on former local domain  

       INTEGER(iwp), SAVE ::  l           !< index variable for surface type

       LOGICAL            ::  surf_match_def     !< flag indicating that surface element is of default type
       LOGICAL            ::  surf_match_lsm     !< flag indicating that surface element is of natural type
       LOGICAL            ::  surf_match_usm     !< flag indicating that surface element is of urban type

       LOGICAL, INTENT(OUT) ::  found

       LOGICAL, SAVE ::  horizontal_surface !< flag indicating horizontal surfaces
       LOGICAL, SAVE ::  vertical_surface   !< flag indicating vertical surfaces

       TYPE(surf_type), DIMENSION(0:2), SAVE ::  surf_h !< horizontal surface type on file
       TYPE(surf_type), DIMENSION(0:3), SAVE ::  surf_v !< vertical surface type on file


       found              = .TRUE.

       SELECT CASE ( restart_string(1:length) )
!
!--       Read the number of horizontally orientated surface elements and 
!--       allocate arrays
          CASE ( 'ns_h_on_file' )
             IF ( kk == 1 )  THEN
                READ ( 13 )  ns_h_on_file

                IF ( ALLOCATED( surf_h(0)%start_index ) )                      &
                   CALL deallocate_surface_attributes_h( surf_h(0) )           
                IF ( ALLOCATED( surf_h(1)%start_index ) )                      &
                   CALL deallocate_surface_attributes_h( surf_h(1) )           
                IF ( ALLOCATED( surf_h(2)%start_index ) )                      &
                   CALL deallocate_surface_attributes_h_top( surf_h(2) )       
!
!--             Allocate memory for number of surface elements on file. 
!--             Please note, these number is not necessarily the same as 
!--             the final number of surface elements on local domain,
!--             which is the case if processor topology changes during
!--             restart runs.  
!--             Horizontal upward facing
                surf_h(0)%ns = ns_h_on_file(0)
                CALL allocate_surface_attributes_h( surf_h(0),                 &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file, nxr_on_file )
!
!--             Horizontal downward facing
                surf_h(1)%ns = ns_h_on_file(1)
                CALL allocate_surface_attributes_h( surf_h(1),                 &
                                        nys_on_file, nyn_on_file,              &
                                        nxl_on_file, nxr_on_file )
!
!--             Model top
                surf_h(2)%ns = ns_h_on_file(2)
                CALL allocate_surface_attributes_h_top( surf_h(2),             &
                                            nys_on_file, nyn_on_file,          &
                                            nxl_on_file, nxr_on_file )

!
!--             Initial setting of flags for horizontal and vertical surfaces,
!--             will be set after start- and end-indices are read. 
                horizontal_surface = .FALSE.
                vertical_surface   = .FALSE.

             ENDIF   
!
!--       Read the number of vertically orientated surface elements and 
!--       allocate arrays
          CASE ( 'ns_v_on_file' )
             IF ( kk == 1 ) THEN
                READ ( 13 )  ns_v_on_file

                DO  l = 0, 3
                   IF ( ALLOCATED( surf_v(l)%start_index ) )                   &
                      CALL deallocate_surface_attributes_v( surf_v(l) )
                ENDDO

                DO  l = 0, 3
                   surf_v(l)%ns = ns_v_on_file(l)
                   CALL allocate_surface_attributes_v( surf_v(l),              &
                                           nys_on_file, nyn_on_file,           &
                                           nxl_on_file, nxr_on_file )
               ENDDO

             ENDIF
!
!--       Read start and end indices of surface elements at each (ji)-gridpoint
          CASE ( 'surf_h(0)%start_index' )
             IF ( kk == 1 )                                                    &
                READ ( 13 )  surf_h(0)%start_index
             l = 0
          CASE ( 'surf_h(0)%end_index' )   
             IF ( kk == 1 )                                                    &
                READ ( 13 )  surf_h(0)%end_index
             horizontal_surface = .TRUE.
             vertical_surface   = .FALSE.
!
!--       Read specific attributes
          CASE ( 'surf_h(0)%us' )         
             IF ( ALLOCATED( surf_h(0)%us )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_h(0)%us
          CASE ( 'surf_h(0)%ts' )         
             IF ( ALLOCATED( surf_h(0)%ts )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_h(0)%ts
          CASE ( 'surf_h(0)%qs' )         
             IF ( ALLOCATED( surf_h(0)%qs )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_h(0)%qs
          CASE ( 'surf_h(0)%ss' )         
             IF ( ALLOCATED( surf_h(0)%ss )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_h(0)%ss
          CASE ( 'surf_h(0)%qcs' )         
             IF ( ALLOCATED( surf_h(0)%qcs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(0)%qcs
          CASE ( 'surf_h(0)%ncs' )         
             IF ( ALLOCATED( surf_h(0)%ncs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(0)%ncs
          CASE ( 'surf_h(0)%qrs' )         
             IF ( ALLOCATED( surf_h(0)%qrs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(0)%qrs
          CASE ( 'surf_h(0)%nrs' )         
             IF ( ALLOCATED( surf_h(0)%nrs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(0)%nrs
          CASE ( 'surf_h(0)%ol' )         
             IF ( ALLOCATED( surf_h(0)%ol )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_h(0)%ol
          CASE ( 'surf_h(0)%rib' )         
             IF ( ALLOCATED( surf_h(0)%rib )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(0)%rib
          CASE ( 'surf_h(0)%pt_surface' )         
             IF ( ALLOCATED( surf_h(0)%pt_surface )  .AND.  kk == 1 )          &
                READ ( 13 )  surf_h(0)%pt_surface
          CASE ( 'surf_h(0)%q_surface' )         
             IF ( ALLOCATED( surf_h(0)%q_surface )  .AND.  kk == 1 )           &
                READ ( 13 )  surf_h(0)%q_surface
          CASE ( 'surf_h(0)%vpt_surface' )         
             IF ( ALLOCATED( surf_h(0)%vpt_surface )  .AND.  kk == 1 )         &
                READ ( 13 )  surf_h(0)%vpt_surface
          CASE ( 'surf_h(0)%usws' )         
             IF ( ALLOCATED( surf_h(0)%usws )  .AND.  kk == 1 )                &
                READ ( 13 )  surf_h(0)%usws
          CASE ( 'surf_h(0)%vsws' )         
             IF ( ALLOCATED( surf_h(0)%vsws )  .AND.  kk == 1 )                &
                READ ( 13 )  surf_h(0)%vsws
          CASE ( 'surf_h(0)%shf' )         
             IF ( ALLOCATED( surf_h(0)%shf )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(0)%shf
          CASE ( 'surf_h(0)%qsws' )         
             IF ( ALLOCATED( surf_h(0)%qsws )  .AND.  kk == 1 )                &
                READ ( 13 )  surf_h(0)%qsws
          CASE ( 'surf_h(0)%ssws' )         
             IF ( ALLOCATED( surf_h(0)%ssws )  .AND.  kk == 1 )                &
                READ ( 13 )  surf_h(0)%ssws
          CASE ( 'surf_h(0)%css' )
             IF ( ALLOCATED( surf_h(0)%css )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(0)%css
          CASE ( 'surf_h(0)%cssws' )         
             IF ( ALLOCATED( surf_h(0)%cssws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_h(0)%cssws
          CASE ( 'surf_h(0)%qcsws' )         
             IF ( ALLOCATED( surf_h(0)%qcsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_h(0)%qcsws
          CASE ( 'surf_h(0)%ncsws' )         
             IF ( ALLOCATED( surf_h(0)%ncsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_h(0)%ncsws
          CASE ( 'surf_h(0)%qrsws' )         
             IF ( ALLOCATED( surf_h(0)%qrsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_h(0)%qrsws
          CASE ( 'surf_h(0)%nrsws' )         
             IF ( ALLOCATED( surf_h(0)%nrsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_h(0)%nrsws
          CASE ( 'surf_h(0)%sasws' )         
             IF ( ALLOCATED( surf_h(0)%sasws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_h(0)%sasws

          CASE ( 'surf_h(1)%start_index' )   
             IF ( kk == 1 )                                                    &
                READ ( 13 )  surf_h(1)%start_index
             l = 1
          CASE ( 'surf_h(1)%end_index' )   
             IF ( kk == 1 )                                                    &
                READ ( 13 )  surf_h(1)%end_index
          CASE ( 'surf_h(1)%us' )         
             IF ( ALLOCATED( surf_h(1)%us )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_h(1)%us
          CASE ( 'surf_h(1)%ts' )         
             IF ( ALLOCATED( surf_h(1)%ts )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_h(1)%ts
          CASE ( 'surf_h(1)%qs' )         
             IF ( ALLOCATED( surf_h(1)%qs )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_h(1)%qs
          CASE ( 'surf_h(1)%ss' )         
             IF ( ALLOCATED( surf_h(1)%ss )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_h(1)%ss
          CASE ( 'surf_h(1)%qcs' )         
             IF ( ALLOCATED( surf_h(1)%qcs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(1)%qcs
          CASE ( 'surf_h(1)%ncs' )         
             IF ( ALLOCATED( surf_h(1)%ncs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(1)%ncs
          CASE ( 'surf_h(1)%qrs' )         
             IF ( ALLOCATED( surf_h(1)%qrs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(1)%qrs
          CASE ( 'surf_h(1)%nrs' )         
             IF ( ALLOCATED( surf_h(1)%nrs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(1)%nrs
          CASE ( 'surf_h(1)%ol' )         
             IF ( ALLOCATED( surf_h(1)%ol )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_h(1)%ol
          CASE ( 'surf_h(1)%rib' )         
             IF ( ALLOCATED( surf_h(1)%rib )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(1)%rib
          CASE ( 'surf_h(1)%pt_surface' )         
             IF ( ALLOCATED( surf_h(1)%pt_surface )  .AND.  kk == 1 )          &
                READ ( 13 )  surf_h(1)%pt_surface
          CASE ( 'surf_h(1)%q_surface' )         
             IF ( ALLOCATED( surf_h(1)%q_surface )  .AND.  kk == 1 )           &
                READ ( 13 )  surf_h(1)%q_surface
          CASE ( 'surf_h(1)%vpt_surface' )         
             IF ( ALLOCATED( surf_h(1)%vpt_surface )  .AND.  kk == 1 )         &
                READ ( 13 )  surf_h(1)%vpt_surface
          CASE ( 'surf_h(1)%usws' )         
             IF ( ALLOCATED( surf_h(1)%usws )  .AND.  kk == 1 )                &
                READ ( 13 )  surf_h(1)%usws
          CASE ( 'surf_h(1)%vsws' )         
             IF ( ALLOCATED( surf_h(1)%vsws )  .AND.  kk == 1 )                &
                READ ( 13 )  surf_h(1)%vsws
          CASE ( 'surf_h(1)%shf' )         
             IF ( ALLOCATED( surf_h(1)%shf )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(1)%shf
          CASE ( 'surf_h(1)%qsws' )         
             IF ( ALLOCATED( surf_h(1)%qsws )  .AND.  kk == 1 )                &
                READ ( 13 )  surf_h(1)%qsws
          CASE ( 'surf_h(1)%ssws' )         
             IF ( ALLOCATED( surf_h(1)%ssws )  .AND.  kk == 1 )                &
                READ ( 13 )  surf_h(1)%ssws
          CASE ( 'surf_h(1)%css' )
             IF ( ALLOCATED( surf_h(1)%css )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(1)%css
          CASE ( 'surf_h(1)%cssws' )         
             IF ( ALLOCATED( surf_h(1)%cssws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_h(1)%cssws
          CASE ( 'surf_h(1)%qcsws' )         
             IF ( ALLOCATED( surf_h(1)%qcsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_h(1)%qcsws
          CASE ( 'surf_h(1)%ncsws' )         
             IF ( ALLOCATED( surf_h(1)%ncsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_h(1)%ncsws
          CASE ( 'surf_h(1)%qrsws' )         
             IF ( ALLOCATED( surf_h(1)%qrsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_h(1)%qrsws
          CASE ( 'surf_h(1)%nrsws' )         
             IF ( ALLOCATED( surf_h(1)%nrsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_h(1)%nrsws
          CASE ( 'surf_h(1)%sasws' )         
             IF ( ALLOCATED( surf_h(1)%sasws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_h(1)%sasws

          CASE ( 'surf_h(2)%start_index' )   
             IF ( kk == 1 )                                                    &
                READ ( 13 )  surf_h(2)%start_index
             l = 2
          CASE ( 'surf_h(2)%end_index' )   
             IF ( kk == 1 )                                                    &
                READ ( 13 )  surf_h(2)%end_index
          CASE ( 'surf_h(2)%us' )         
             IF ( ALLOCATED( surf_h(2)%us )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_h(2)%us
          CASE ( 'surf_h(2)%ts' )         
             IF ( ALLOCATED( surf_h(2)%ts )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_h(2)%ts
          CASE ( 'surf_h(2)%qs' )        
             IF ( ALLOCATED( surf_h(2)%qs )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_h(2)%qs
          CASE ( 'surf_h(2)%ss' )         
             IF ( ALLOCATED( surf_h(2)%ss )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_h(2)%ss
          CASE ( 'surf_h(2)%qcs' )         
             IF ( ALLOCATED( surf_h(2)%qcs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(2)%qcs
          CASE ( 'surf_h(2)%ncs' )         
             IF ( ALLOCATED( surf_h(2)%ncs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(2)%ncs
          CASE ( 'surf_h(2)%qrs' )         
             IF ( ALLOCATED( surf_h(2)%qrs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(2)%qrs
          CASE ( 'surf_h(2)%nrs' )         
             IF ( ALLOCATED( surf_h(2)%nrs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(2)%nrs
          CASE ( 'surf_h(2)%ol' )         
             IF ( ALLOCATED( surf_h(2)%ol )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_h(2)%ol
          CASE ( 'surf_h(2)%rib' )         
             IF ( ALLOCATED( surf_h(2)%rib )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(2)%rib
          CASE ( 'surf_h(2)%pt_surface' )         
             IF ( ALLOCATED( surf_h(2)%pt_surface )  .AND.  kk == 1 )          &
                READ ( 13 )  surf_h(2)%pt_surface
          CASE ( 'surf_h(2)%q_surface' )         
             IF ( ALLOCATED( surf_h(2)%q_surface )  .AND.  kk == 1 )           &
                READ ( 13 )  surf_h(2)%q_surface
          CASE ( 'surf_h(2)%vpt_surface' )         
             IF ( ALLOCATED( surf_h(2)%vpt_surface )  .AND.  kk == 1 )         &
                READ ( 13 )  surf_h(2)%vpt_surface
          CASE ( 'surf_h(2)%usws' )         
             IF ( ALLOCATED( surf_h(2)%usws )  .AND.  kk == 1 )                &
                READ ( 13 )  surf_h(2)%usws
          CASE ( 'surf_h(2)%vsws' )         
             IF ( ALLOCATED( surf_h(2)%vsws )  .AND.  kk == 1 )                &
                READ ( 13 )  surf_h(2)%vsws
          CASE ( 'surf_h(2)%shf' )         
             IF ( ALLOCATED( surf_h(2)%shf )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(2)%shf
          CASE ( 'surf_h(2)%qsws' )         
             IF ( ALLOCATED( surf_h(2)%qsws )  .AND.  kk == 1 )                &
                READ ( 13 )  surf_h(2)%qsws
          CASE ( 'surf_h(2)%ssws' )         
             IF ( ALLOCATED( surf_h(2)%ssws )  .AND.  kk == 1 )                &
                READ ( 13 )  surf_h(2)%ssws
          CASE ( 'surf_h(2)%css' )
             IF ( ALLOCATED( surf_h(2)%css )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_h(2)%css
          CASE ( 'surf_h(2)%cssws' )         
             IF ( ALLOCATED( surf_h(2)%cssws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_h(2)%cssws
          CASE ( 'surf_h(2)%qcsws' )         
             IF ( ALLOCATED( surf_h(2)%qcsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_h(2)%qcsws
          CASE ( 'surf_h(2)%ncsws' )         
             IF ( ALLOCATED( surf_h(2)%ncsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_h(2)%ncsws
          CASE ( 'surf_h(2)%qrsws' )         
             IF ( ALLOCATED( surf_h(2)%qrsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_h(2)%qrsws
          CASE ( 'surf_h(2)%nrsws' )         
             IF ( ALLOCATED( surf_h(2)%nrsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_h(2)%nrsws
          CASE ( 'surf_h(2)%sasws' )         
             IF ( ALLOCATED( surf_h(2)%sasws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_h(2)%sasws

          CASE ( 'surf_v(0)%start_index' )   
             IF ( kk == 1 )                                                    &
                READ ( 13 )  surf_v(0)%start_index
             l = 0
             horizontal_surface = .FALSE.
             vertical_surface   = .TRUE.
          CASE ( 'surf_v(0)%end_index' )   
             IF ( kk == 1 )                                                    &
                READ ( 13 )  surf_v(0)%end_index
          CASE ( 'surf_v(0)%us' )         
             IF ( ALLOCATED( surf_v(0)%us )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(0)%us
          CASE ( 'surf_v(0)%ts' )         
             IF ( ALLOCATED( surf_v(0)%ts )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(0)%ts
          CASE ( 'surf_v(0)%qs' )         
             IF ( ALLOCATED( surf_v(0)%qs )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(0)%qs
          CASE ( 'surf_v(0)%ss' )         
             IF ( ALLOCATED( surf_v(0)%ss )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(0)%ss
          CASE ( 'surf_v(0)%qcs' )         
             IF ( ALLOCATED( surf_v(0)%qcs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(0)%qcs
          CASE ( 'surf_v(0)%ncs' )         
             IF ( ALLOCATED( surf_v(0)%ncs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(0)%ncs
          CASE ( 'surf_v(0)%qrs' )         
             IF ( ALLOCATED( surf_v(0)%qrs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(0)%qrs
          CASE ( 'surf_v(0)%nrs' )         
             IF ( ALLOCATED( surf_v(0)%nrs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(0)%nrs
          CASE ( 'surf_v(0)%ol' )         
             IF ( ALLOCATED( surf_v(0)%ol )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(0)%ol
          CASE ( 'surf_v(0)%rib' )         
             IF ( ALLOCATED( surf_v(0)%rib )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(0)%rib
          CASE ( 'surf_v(0)%pt_surface' )         
             IF ( ALLOCATED( surf_v(0)%pt_surface )  .AND.  kk == 1 )          &
                READ ( 13 )  surf_v(0)%pt_surface
          CASE ( 'surf_v(0)%q_surface' )         
             IF ( ALLOCATED( surf_v(0)%q_surface )  .AND.  kk == 1 )           &
                READ ( 13 )  surf_v(0)%q_surface
          CASE ( 'surf_v(0)%vpt_surface' )         
             IF ( ALLOCATED( surf_v(0)%vpt_surface )  .AND.  kk == 1 )         &
                READ ( 13 )  surf_v(0)%vpt_surface
          CASE ( 'surf_v(0)%shf' )         
             IF ( ALLOCATED( surf_v(0)%shf )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(0)%shf
          CASE ( 'surf_v(0)%qsws' )         
             IF ( ALLOCATED( surf_v(0)%qsws )  .AND.  kk == 1 )                &
                READ ( 13 )  surf_v(0)%qsws
          CASE ( 'surf_v(0)%ssws' )         
             IF ( ALLOCATED( surf_v(0)%ssws )  .AND.  kk == 1 )                &
                READ ( 13 )  surf_v(0)%ssws
          CASE ( 'surf_v(0)%css' ) 
             IF ( ALLOCATED( surf_v(0)%css )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(0)%css
          CASE ( 'surf_v(0)%cssws' )         
             IF ( ALLOCATED( surf_v(0)%cssws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(0)%cssws
          CASE ( 'surf_v(0)%qcsws' )         
             IF ( ALLOCATED( surf_v(0)%qcsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(0)%qcsws
          CASE ( 'surf_v(0)%ncsws' )         
             IF ( ALLOCATED( surf_v(0)%ncsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(0)%ncsws
          CASE ( 'surf_v(0)%qrsws' )         
             IF ( ALLOCATED( surf_v(0)%qrsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(0)%qrsws
          CASE ( 'surf_v(0)%nrsws' )         
             IF ( ALLOCATED( surf_v(0)%nrsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(0)%nrsws
          CASE ( 'surf_v(0)%sasws' )         
             IF ( ALLOCATED( surf_v(0)%sasws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(0)%sasws
          CASE ( 'surf_v(0)%mom_uv' )         
             IF ( ALLOCATED( surf_v(0)%mom_flux_uv )  .AND.  kk == 1 )         &
                READ ( 13 )  surf_v(0)%mom_flux_uv
          CASE ( 'surf_v(0)%mom_w' )         
             IF ( ALLOCATED( surf_v(0)%mom_flux_w )  .AND.  kk == 1 )          &
                READ ( 13 )  surf_v(0)%mom_flux_w
          CASE ( 'surf_v(0)%mom_tke' )         
             IF ( ALLOCATED( surf_v(0)%mom_flux_tke )  .AND.  kk == 1 )        &
                READ ( 13 )  surf_v(0)%mom_flux_tke

          CASE ( 'surf_v(1)%start_index' )   
             IF ( kk == 1 )                                                    &
                READ ( 13 )  surf_v(1)%start_index
             l = 1
          CASE ( 'surf_v(1)%end_index' )   
             IF ( kk == 1 )                                                    &
                READ ( 13 )  surf_v(1)%end_index
          CASE ( 'surf_v(1)%us' )         
             IF ( ALLOCATED( surf_v(1)%us )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(1)%us
          CASE ( 'surf_v(1)%ts' )         
             IF ( ALLOCATED( surf_v(1)%ts )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(1)%ts
          CASE ( 'surf_v(1)%qs' )         
             IF ( ALLOCATED( surf_v(1)%qs )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(1)%qs
          CASE ( 'surf_v(1)%ss' )         
             IF ( ALLOCATED( surf_v(1)%ss )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(1)%ss
          CASE ( 'surf_v(1)%qcs' )         
             IF ( ALLOCATED( surf_v(1)%qcs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(1)%qcs
          CASE ( 'surf_v(1)%ncs' )         
             IF ( ALLOCATED( surf_v(1)%ncs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(1)%ncs
          CASE ( 'surf_v(1)%qrs' )         
             IF ( ALLOCATED( surf_v(1)%qrs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(1)%qrs
          CASE ( 'surf_v(1)%nrs' )         
             IF ( ALLOCATED( surf_v(1)%nrs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(1)%nrs
          CASE ( 'surf_v(1)%ol' )         
             IF ( ALLOCATED( surf_v(1)%ol )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(1)%ol
          CASE ( 'surf_v(1)%rib' )         
             IF ( ALLOCATED( surf_v(1)%rib )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(1)%rib
          CASE ( 'surf_v(1)%pt_surface' )         
             IF ( ALLOCATED( surf_v(1)%pt_surface )  .AND.  kk == 1 )          &
                READ ( 13 )  surf_v(1)%pt_surface
          CASE ( 'surf_v(1)%q_surface' )         
             IF ( ALLOCATED( surf_v(1)%q_surface )  .AND.  kk == 1 )           &
                READ ( 13 )  surf_v(1)%q_surface
          CASE ( 'surf_v(1)%vpt_surface' )         
             IF ( ALLOCATED( surf_v(1)%vpt_surface )  .AND.  kk == 1 )         &
                READ ( 13 )  surf_v(1)%vpt_surface
          CASE ( 'surf_v(1)%shf' )         
             IF ( ALLOCATED( surf_v(1)%shf )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(1)%shf
          CASE ( 'surf_v(1)%qsws' )         
             IF ( ALLOCATED( surf_v(1)%qsws )  .AND.  kk == 1 )                &
                READ ( 13 )  surf_v(1)%qsws
          CASE ( 'surf_v(1)%ssws' )         
             IF ( ALLOCATED( surf_v(1)%ssws )  .AND.  kk == 1 )                &
                READ ( 13 )  surf_v(1)%ssws
          CASE ( 'surf_v(1)%css' ) 
             IF ( ALLOCATED( surf_v(1)%css )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(1)%css
          CASE ( 'surf_v(1)%cssws' )         
             IF ( ALLOCATED( surf_v(1)%cssws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(1)%cssws
          CASE ( 'surf_v(1)%qcsws' )         
             IF ( ALLOCATED( surf_v(1)%qcsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(1)%qcsws
          CASE ( 'surf_v(1)%ncsws' )         
             IF ( ALLOCATED( surf_v(1)%ncsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(1)%ncsws
          CASE ( 'surf_v(1)%qrsws' )         
             IF ( ALLOCATED( surf_v(1)%qrsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(1)%qrsws
          CASE ( 'surf_v(1)%nrsws' )         
             IF ( ALLOCATED( surf_v(1)%nrsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(1)%nrsws
          CASE ( 'surf_v(1)%sasws' )         
             IF ( ALLOCATED( surf_v(1)%sasws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(1)%sasws
          CASE ( 'surf_v(1)%mom_uv' )         
             IF ( ALLOCATED( surf_v(1)%mom_flux_uv )  .AND.  kk == 1 )         &
                READ ( 13 )  surf_v(1)%mom_flux_uv
          CASE ( 'surf_v(1)%mom_w' )         
             IF ( ALLOCATED( surf_v(1)%mom_flux_w )  .AND.  kk == 1 )          &
                READ ( 13 )  surf_v(1)%mom_flux_w
          CASE ( 'surf_v(1)%mom_tke' )         
             IF ( ALLOCATED( surf_v(1)%mom_flux_tke )  .AND.  kk == 1 )        &
                READ ( 13 )  surf_v(1)%mom_flux_tke

          CASE ( 'surf_v(2)%start_index' )   
             IF ( kk == 1 )                                                    &
                READ ( 13 )  surf_v(2)%start_index
             l = 2
          CASE ( 'surf_v(2)%end_index' )   
             IF ( kk == 1 )                                                    &
                READ ( 13 )  surf_v(2)%end_index
          CASE ( 'surf_v(2)%us' )         
             IF ( ALLOCATED( surf_v(2)%us )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(2)%us
          CASE ( 'surf_v(2)%ts' )         
             IF ( ALLOCATED( surf_v(2)%ts )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(2)%ts
          CASE ( 'surf_v(2)%qs' )         
             IF ( ALLOCATED( surf_v(2)%qs )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(2)%qs
          CASE ( 'surf_v(2)%ss' )         
             IF ( ALLOCATED( surf_v(2)%ss )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(2)%ss
          CASE ( 'surf_v(2)%qcs' )         
             IF ( ALLOCATED( surf_v(2)%qcs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(2)%qcs
          CASE ( 'surf_v(2)%ncs' )         
             IF ( ALLOCATED( surf_v(2)%ncs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(2)%ncs
          CASE ( 'surf_v(2)%qrs' )         
             IF ( ALLOCATED( surf_v(2)%qrs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(2)%qrs
          CASE ( 'surf_v(2)%nrs' )         
             IF ( ALLOCATED( surf_v(2)%nrs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(2)%nrs
          CASE ( 'surf_v(2)%ol' )         
             IF ( ALLOCATED( surf_v(2)%ol )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(2)%ol
          CASE ( 'surf_v(2)%rib' )         
             IF ( ALLOCATED( surf_v(2)%rib )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(2)%rib
          CASE ( 'surf_v(2)%pt_surface' )         
             IF ( ALLOCATED( surf_v(2)%pt_surface )  .AND.  kk == 1 )          &
                READ ( 13 )  surf_v(2)%pt_surface
          CASE ( 'surf_v(2)%q_surface' )         
             IF ( ALLOCATED( surf_v(2)%q_surface )  .AND.  kk == 1 )           &
                READ ( 13 )  surf_v(2)%q_surface
          CASE ( 'surf_v(2)%vpt_surface' )         
             IF ( ALLOCATED( surf_v(2)%vpt_surface )  .AND.  kk == 1 )         &
                READ ( 13 )  surf_v(2)%vpt_surface
          CASE ( 'surf_v(2)%shf' )         
             IF ( ALLOCATED( surf_v(2)%shf )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(2)%shf
          CASE ( 'surf_v(2)%qsws' )         
             IF ( ALLOCATED( surf_v(2)%qsws )  .AND.  kk == 1 )                &
                READ ( 13 )  surf_v(2)%qsws
          CASE ( 'surf_v(2)%ssws' )         
             IF ( ALLOCATED( surf_v(2)%ssws )  .AND.  kk == 1 )                &
                READ ( 13 )  surf_v(2)%ssws
          CASE ( 'surf_v(2)%css' ) 
             IF ( ALLOCATED( surf_v(2)%css )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(2)%css
          CASE ( 'surf_v(2)%cssws' )         
             IF ( ALLOCATED( surf_v(2)%cssws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(2)%cssws
          CASE ( 'surf_v(2)%qcsws' )         
             IF ( ALLOCATED( surf_v(2)%qcsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(2)%qcsws
          CASE ( 'surf_v(2)%ncsws' )         
             IF ( ALLOCATED( surf_v(2)%ncsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(2)%ncsws
          CASE ( 'surf_v(2)%qrsws' )         
             IF ( ALLOCATED( surf_v(2)%qrsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(2)%qrsws
          CASE ( 'surf_v(2)%nrsws' )         
             IF ( ALLOCATED( surf_v(2)%nrsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(2)%nrsws
          CASE ( 'surf_v(2)%sasws' )         
             IF ( ALLOCATED( surf_v(2)%sasws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(2)%sasws
          CASE ( 'surf_v(2)%mom_uv' )         
             IF ( ALLOCATED( surf_v(2)%mom_flux_uv )  .AND.  kk == 1 )         &
                READ ( 13 )  surf_v(2)%mom_flux_uv
          CASE ( 'surf_v(2)%mom_w' )         
             IF ( ALLOCATED( surf_v(2)%mom_flux_w )  .AND.  kk == 1 )          &
                READ ( 13 )  surf_v(2)%mom_flux_w
          CASE ( 'surf_v(2)%mom_tke' )         
             IF ( ALLOCATED( surf_v(2)%mom_flux_tke )  .AND.  kk == 1 )        &
                READ ( 13 )  surf_v(2)%mom_flux_tke

          CASE ( 'surf_v(3)%start_index' )   
             IF ( kk == 1 )                                                    &
                READ ( 13 )  surf_v(3)%start_index
             l = 3
          CASE ( 'surf_v(3)%end_index' )   
             IF ( kk == 1 )                                                    &
                READ ( 13 )  surf_v(3)%end_index
          CASE ( 'surf_v(3)%us' )         
             IF ( ALLOCATED( surf_v(3)%us )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(3)%us
          CASE ( 'surf_v(3)%ts' )         
             IF ( ALLOCATED( surf_v(3)%ts )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(3)%ts
          CASE ( 'surf_v(3)%qs' )        
             IF ( ALLOCATED( surf_v(3)%qs )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(3)%qs
          CASE ( 'surf_v(3)%ss' )         
             IF ( ALLOCATED( surf_v(3)%ss )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(3)%ss
          CASE ( 'surf_v(3)%qcs' )         
             IF ( ALLOCATED( surf_v(3)%qcs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(3)%qcs
          CASE ( 'surf_v(3)%ncs' )         
             IF ( ALLOCATED( surf_v(3)%ncs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(3)%ncs
          CASE ( 'surf_v(3)%qrs' )         
             IF ( ALLOCATED( surf_v(3)%qrs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(3)%qrs
          CASE ( 'surf_v(3)%nrs' )         
             IF ( ALLOCATED( surf_v(3)%nrs )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(3)%nrs
          CASE ( 'surf_v(3)%ol' )         
             IF ( ALLOCATED( surf_v(3)%ol )  .AND.  kk == 1 )                  &
                READ ( 13 )  surf_v(3)%ol
          CASE ( 'surf_v(3)%rib' )         
             IF ( ALLOCATED( surf_v(3)%rib )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(3)%rib
          CASE ( 'surf_v(3)%pt_surface' )         
             IF ( ALLOCATED( surf_v(3)%pt_surface )  .AND.  kk == 1 )          &
                READ ( 13 )  surf_v(3)%pt_surface
          CASE ( 'surf_v(3)%q_surface' )         
             IF ( ALLOCATED( surf_v(3)%q_surface )  .AND.  kk == 1 )           &
                READ ( 13 )  surf_v(3)%q_surface
          CASE ( 'surf_v(3)%vpt_surface' )         
             IF ( ALLOCATED( surf_v(3)%vpt_surface )  .AND.  kk == 1 )         &
                READ ( 13 )  surf_v(3)%vpt_surface
          CASE ( 'surf_v(3)%shf' )         
             IF ( ALLOCATED( surf_v(3)%shf )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(3)%shf
          CASE ( 'surf_v(3)%qsws' )         
             IF ( ALLOCATED( surf_v(3)%qsws )  .AND.  kk == 1 )                & 
                READ ( 13 )  surf_v(3)%qsws
          CASE ( 'surf_v(3)%ssws' )         
             IF ( ALLOCATED( surf_v(3)%ssws )  .AND.  kk == 1 )                &
                READ ( 13 )  surf_v(3)%ssws
          CASE ( 'surf_v(3)%css' ) 
             IF ( ALLOCATED( surf_v(3)%css )  .AND.  kk == 1 )                 &
                READ ( 13 )  surf_v(3)%css
          CASE ( 'surf_v(3)%cssws' )         
             IF ( ALLOCATED( surf_v(3)%cssws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(3)%cssws
          CASE ( 'surf_v(3)%qcsws' )         
             IF ( ALLOCATED( surf_v(3)%qcsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(3)%qcsws
          CASE ( 'surf_v(3)%ncsws' )         
             IF ( ALLOCATED( surf_v(3)%ncsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(3)%ncsws
          CASE ( 'surf_v(3)%qrsws' )         
             IF ( ALLOCATED( surf_v(3)%qrsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(3)%qrsws
          CASE ( 'surf_v(3)%nrsws' )         
             IF ( ALLOCATED( surf_v(3)%nrsws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(3)%nrsws
          CASE ( 'surf_v(3)%sasws' )         
             IF ( ALLOCATED( surf_v(3)%sasws )  .AND.  kk == 1 )               &
                READ ( 13 )  surf_v(3)%sasws
          CASE ( 'surf_v(3)%mom_uv' )         
             IF ( ALLOCATED( surf_v(3)%mom_flux_uv )  .AND.  kk == 1 )         &
                READ ( 13 )  surf_v(3)%mom_flux_uv
          CASE ( 'surf_v(3)%mom_w' )         
             IF ( ALLOCATED( surf_v(3)%mom_flux_w )  .AND.  kk == 1 )          &
                READ ( 13 )  surf_v(3)%mom_flux_w
          CASE ( 'surf_v(3)%mom_tke' )         
             IF ( ALLOCATED( surf_v(3)%mom_flux_tke )  .AND.  kk == 1 )        &
                READ ( 13 )  surf_v(3)%mom_flux_tke

          CASE DEFAULT

                found = .FALSE.

       END SELECT
!
!--    Redistribute surface elements on its respective type. Start with 
!--    horizontally orientated surfaces.
       IF ( horizontal_surface  .AND.                                          &
            .NOT. INDEX( restart_string(1:length), '%start_index' ) /= 0 )     &
       THEN
       
          ic = nxlc
          DO  i = nxlf, nxrf
             jc = nysc
             DO  j = nysf, nynf
!
!--             Determine type of surface element, i.e. default, natural, 
!--             urban, at current grid point. 
                surf_match_def  = surf_def_h(l)%end_index(jc,ic) >=            &
                                  surf_def_h(l)%start_index(jc,ic)
                surf_match_lsm  = ( surf_lsm_h%end_index(jc,ic)  >=            &
                                    surf_lsm_h%start_index(jc,ic) )            &
                             .AND.  l == 0 
                surf_match_usm  = ( surf_usm_h%end_index(jc,ic)  >=            &
                                    surf_usm_h%start_index(jc,ic) )            &
                             .AND.  l == 0 
!
!--             Write restart data onto default-type surfaces if required.
                IF ( surf_match_def )  THEN
!
!--                Set the start index for the local surface element 
                   mm = surf_def_h(l)%start_index(jc,ic)
!
!--                For index pair (j,i) on file loop from start to end index,
!--                and in case the local surface element mm is smaller than
!--                the local end index, assign the respective surface data
!--                to this element. 
                   DO  m = surf_h(l)%start_index(j,i),                         &
                           surf_h(l)%end_index(j,i)
                      IF ( surf_def_h(l)%end_index(jc,ic) >= mm )              &
                         CALL restore_surface_elements( surf_def_h(l),         &
                                                        mm, surf_h(l), m )
                      mm = mm + 1
                   ENDDO
                ENDIF
!
!--             Same for natural-type surfaces. Please note, it is implicitly
!--             assumed that natural surface elements are below urban 
!--             urban surface elements if there are several horizontal surfaces
!--             at (j,i). An example would be bridges. 
                IF ( surf_match_lsm )  THEN
                   mm = surf_lsm_h%start_index(jc,ic)
                   DO  m = surf_h(l)%start_index(j,i),                         &
                           surf_h(l)%end_index(j,i)
                      IF ( surf_lsm_h%end_index(jc,ic) >= mm )                 &
                         CALL restore_surface_elements( surf_lsm_h,            &
                                                        mm, surf_h(l), m )
                      mm = mm + 1
                   ENDDO
                ENDIF
!
!--             Same for urban-type surfaces
                IF ( surf_match_usm )  THEN
                   mm = surf_usm_h%start_index(jc,ic)
                   DO  m = surf_h(l)%start_index(j,i),                         &
                           surf_h(l)%end_index(j,i)
                      IF ( surf_usm_h%end_index(jc,ic) >= mm )                 &
                         CALL restore_surface_elements( surf_usm_h,            &
                                                        mm, surf_h(l), m )
                      mm = mm + 1
                   ENDDO
                ENDIF

                jc = jc + 1
             ENDDO
             ic = ic + 1
          ENDDO
       ELSEIF ( vertical_surface  .AND.                                        &
            .NOT. INDEX( restart_string(1:length), '%start_index' ) /= 0 )     &
       THEN
          ic = nxlc
          DO  i = nxlf, nxrf
             jc = nysc
             DO  j = nysf, nynf
!
!--             Determine type of surface element, i.e. default, natural, 
!--             urban, at current grid point. 
                surf_match_def  = surf_def_v(l)%end_index(jc,ic) >=            &
                                  surf_def_v(l)%start_index(jc,ic)
                surf_match_lsm  = surf_lsm_v(l)%end_index(jc,ic) >=            &
                                  surf_lsm_v(l)%start_index(jc,ic)
                surf_match_usm  = surf_usm_v(l)%end_index(jc,ic) >=            &
                                  surf_usm_v(l)%start_index(jc,ic)
!
!--             Write restart data onto default-type surfaces if required.               
                IF ( surf_match_def )  THEN
!
!--                Set the start index for the local surface element 
                   mm = surf_def_v(l)%start_index(jc,ic)
!
!--                For index pair (j,i) on file loop from start to end index,
!--                and in case the local surface element mm is smaller than
!--                the local end index, assign the respective surface data
!--                to this element. 
                   DO  m = surf_v(l)%start_index(j,i),                         &
                           surf_v(l)%end_index(j,i)
                      IF ( surf_def_v(l)%end_index(jc,ic) >= mm )              &
                         CALL restore_surface_elements( surf_def_v(l),         &
                                                        mm, surf_v(l), m )
                      mm = mm + 1
                   ENDDO
                ENDIF
!
!--             Same for natural-type surfaces. Please note, it is implicitly
!--             assumed that natural surface elements are below urban 
!--             urban surface elements if there are several vertical surfaces
!--             at (j,i). An example a terrain elevations with a building on 
!--             top. So far, initialization of urban surfaces below natural
!--             surfaces on the same (j,i) is not possible, so that this case
!--             cannot occur. 
                IF ( surf_match_lsm )  THEN
                   mm = surf_lsm_v(l)%start_index(jc,ic)
                   DO  m = surf_v(l)%start_index(j,i),                         &
                           surf_v(l)%end_index(j,i)
                      IF ( surf_lsm_v(l)%end_index(jc,ic) >= mm )              &
                         CALL restore_surface_elements( surf_lsm_v(l),         &
                                                        mm, surf_v(l), m )
                      mm = mm + 1
                   ENDDO
                ENDIF

                IF ( surf_match_usm )  THEN
                   mm = surf_usm_v(l)%start_index(jc,ic)
                   DO  m = surf_v(l)%start_index(j,i),                         &
                           surf_v(l)%end_index(j,i)
                      IF ( surf_usm_v(l)%end_index(jc,ic) >= mm )              &
                         CALL restore_surface_elements( surf_usm_v(l),         &
                                                        mm, surf_v(l), m )
                      mm = mm + 1
                   ENDDO
                ENDIF

                jc = jc + 1
             ENDDO
             ic = ic + 1
          ENDDO
       ENDIF


    CONTAINS
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Restores surface elements back on its respective type.
!------------------------------------------------------------------------------!
          SUBROUTINE restore_surface_elements( surf_target, m_target,          &
                                               surf_file,   m_file )

             IMPLICIT NONE

             INTEGER(iwp)      ::  m_file      !< respective surface-element index of current surface array 
             INTEGER(iwp)      ::  m_target    !< respecitve surface-element index of surface array on file
             INTEGER(iwp)      ::  lsp         !< running index chemical species

             TYPE( surf_type ) ::  surf_target !< target surface type
             TYPE( surf_type ) ::  surf_file   !< surface type on file


             IF ( INDEX( restart_string(1:length), '%us' ) /= 0 )  THEN 
                IF ( ALLOCATED( surf_target%us )  .AND.                        &
                     ALLOCATED( surf_file%us   ) )                             & 
                   surf_target%us(m_target) = surf_file%us(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%ol' ) /= 0 )  THEN 
                IF ( ALLOCATED( surf_target%ol )  .AND.                        &
                     ALLOCATED( surf_file%ol   ) )                             & 
                   surf_target%ol(m_target) = surf_file%ol(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%pt_surface' ) /= 0 )  THEN 
                IF ( ALLOCATED( surf_target%pt_surface )  .AND.                &
                     ALLOCATED( surf_file%pt_surface   ) )                     & 
                   surf_target%pt_surface(m_target) = surf_file%pt_surface(m_file)
             ENDIF
             
             IF ( INDEX( restart_string(1:length), '%q_surface' ) /= 0 )  THEN 
                IF ( ALLOCATED( surf_target%q_surface )  .AND.                 &
                     ALLOCATED( surf_file%q_surface   ) )                      & 
                   surf_target%q_surface(m_target) = surf_file%q_surface(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%vpt_surface' ) /= 0 )  THEN 
                IF ( ALLOCATED( surf_target%vpt_surface )  .AND.               &
                     ALLOCATED( surf_file%vpt_surface   ) )                    & 
                   surf_target%vpt_surface(m_target) = surf_file%vpt_surface(m_file)
             ENDIF
             
             IF ( INDEX( restart_string(1:length), '%usws' ) /= 0 )  THEN 
                IF ( ALLOCATED( surf_target%usws )  .AND.                      &
                     ALLOCATED( surf_file%usws   ) )                           & 
                   surf_target%usws(m_target) = surf_file%usws(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%vsws' ) /= 0 )  THEN 
                IF ( ALLOCATED( surf_target%vsws )  .AND.                      &
                     ALLOCATED( surf_file%vsws   ) )                           & 
                   surf_target%vsws(m_target) = surf_file%vsws(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%ts' ) /= 0 )  THEN 
                IF ( ALLOCATED( surf_target%ts )  .AND.                        &
                     ALLOCATED( surf_file%ts   ) )                             & 
                   surf_target%ts(m_target) = surf_file%ts(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%shf' ) /= 0 )  THEN 
                IF ( ALLOCATED( surf_target%shf )  .AND.                       &
                     ALLOCATED( surf_file%shf   ) )                            & 
                   surf_target%shf(m_target) = surf_file%shf(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%qs' ) /= 0 )  THEN 
                IF ( ALLOCATED( surf_target%qs )  .AND.                        &
                     ALLOCATED( surf_file%qs   ) )                             & 
                   surf_target%qs(m_target) = surf_file%qs(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%qsws' ) /= 0 )  THEN 
                IF ( ALLOCATED( surf_target%qsws )  .AND.                      &
                     ALLOCATED( surf_file%qsws   ) )                           & 
                   surf_target%qsws(m_target) = surf_file%qsws(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%ss' ) /= 0 )  THEN 
                IF ( ALLOCATED( surf_target%ss )  .AND.                        &
                     ALLOCATED( surf_file%ss   ) )                             & 
                   surf_target%ss(m_target) = surf_file%ss(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%ssws' ) /= 0 )  THEN 
                IF ( ALLOCATED( surf_target%ssws )  .AND.                      &
                     ALLOCATED( surf_file%ssws   ) )                           & 
                   surf_target%ssws(m_target) = surf_file%ssws(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%css' ) /= 0 )  THEN 
                IF ( ALLOCATED( surf_target%css )  .AND.                     &
                     ALLOCATED( surf_file%css   ) )  THEN                  
                   DO  lsp = 1, nvar
                      surf_target%css(lsp,m_target) = surf_file%css(lsp,m_file)
                   ENDDO
                ENDIF
             ENDIF
             IF ( INDEX( restart_string(1:length), '%cssws' ) /= 0 )  THEN 
                IF ( ALLOCATED( surf_target%cssws )  .AND.                     &
                     ALLOCATED( surf_file%cssws   ) )  THEN                  
                   DO  lsp = 1, nvar
                      surf_target%cssws(lsp,m_target) = surf_file%cssws(lsp,m_file)
                   ENDDO
                ENDIF
             ENDIF

             IF ( INDEX( restart_string(1:length), '%qcs' ) /= 0 )  THEN 
                IF ( ALLOCATED( surf_target%qcs )  .AND.                       &
                     ALLOCATED( surf_file%qcs   ) )                            & 
                  surf_target%qcs(m_target) = surf_file%qcs(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%qcsws' ) /= 0 )  THEN 
                IF ( ALLOCATED( surf_target%qcsws )  .AND.                     &
                     ALLOCATED( surf_file%qcsws   ) )                          & 
                   surf_target%qcsws(m_target) = surf_file%qcsws(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%ncs' ) /= 0 )  THEN 
                IF ( ALLOCATED( surf_target%ncs )  .AND.                       &
                     ALLOCATED( surf_file%ncs   ) )                            & 
                   surf_target%ncs(m_target) = surf_file%ncs(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%ncsws' ) /= 0 )  THEN 
                IF ( ALLOCATED( surf_target%ncsws )  .AND.                     &
                     ALLOCATED( surf_file%ncsws   ) )                          & 
                   surf_target%ncsws(m_target) = surf_file%ncsws(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%qrs' ) /= 0 )  THEN 
                IF ( ALLOCATED( surf_target%qrs )  .AND.                       &
                     ALLOCATED( surf_file%qrs   ) )                            & 
                  surf_target%qrs(m_target) = surf_file%qrs(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%qrsws' ) /= 0 )  THEN 
                IF ( ALLOCATED( surf_target%qrsws )  .AND.                     &
                     ALLOCATED( surf_file%qrsws   ) )                          & 
                   surf_target%qrsws(m_target) = surf_file%qrsws(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%nrs' ) /= 0 )  THEN 
                IF ( ALLOCATED( surf_target%nrs )  .AND.                       &
                     ALLOCATED( surf_file%nrs   ) )                            & 
                   surf_target%nrs(m_target) = surf_file%nrs(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%nrsws' ) /= 0 )  THEN 
                IF ( ALLOCATED( surf_target%nrsws )  .AND.                     &
                     ALLOCATED( surf_file%nrsws   ) )                          & 
                   surf_target%nrsws(m_target) = surf_file%nrsws(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%sasws' ) /= 0 )  THEN 
                IF ( ALLOCATED( surf_target%sasws )  .AND.                     &
                     ALLOCATED( surf_file%sasws   ) )                          & 
                   surf_target%sasws(m_target) = surf_file%sasws(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%mom_uv' ) /= 0 )  THEN 
                IF ( ALLOCATED( surf_target%mom_flux_uv )  .AND.               &
                     ALLOCATED( surf_file%mom_flux_uv   ) )                    & 
                   surf_target%mom_flux_uv(m_target) =                         &
                                           surf_file%mom_flux_uv(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%mom_w' ) /= 0 )  THEN 
                IF ( ALLOCATED( surf_target%mom_flux_w )  .AND.                &
                     ALLOCATED( surf_file%mom_flux_w   ) )                     & 
                   surf_target%mom_flux_w(m_target) =                          &
                                           surf_file%mom_flux_w(m_file)
             ENDIF

             IF ( INDEX( restart_string(1:length), '%mom_tke' ) /= 0 )  THEN 
                IF ( ALLOCATED( surf_target%mom_flux_tke )  .AND.              &
                     ALLOCATED( surf_file%mom_flux_tke   ) )                   & 
                   surf_target%mom_flux_tke(0:1,m_target) =                    &
                                           surf_file%mom_flux_tke(0:1,m_file)
             ENDIF


          END SUBROUTINE restore_surface_elements


    END SUBROUTINE surface_rrd_local

 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Counts the number of surface elements with the same facing, required for 
!> reading and writing restart data.
!------------------------------------------------------------------------------!
    SUBROUTINE surface_last_actions

       IMPLICIT NONE
!
!--    Horizontal surfaces
       ns_h_on_file(0) = surf_def_h(0)%ns + surf_lsm_h%ns + surf_usm_h%ns
       ns_h_on_file(1) = surf_def_h(1)%ns
       ns_h_on_file(2) = surf_def_h(2)%ns
!
!--    Vertical surfaces
       ns_v_on_file(0) = surf_def_v(0)%ns + surf_lsm_v(0)%ns + surf_usm_v(0)%ns
       ns_v_on_file(1) = surf_def_v(1)%ns + surf_lsm_v(1)%ns + surf_usm_v(1)%ns
       ns_v_on_file(2) = surf_def_v(2)%ns + surf_lsm_v(2)%ns + surf_usm_v(2)%ns
       ns_v_on_file(3) = surf_def_v(3)%ns + surf_lsm_v(3)%ns + surf_usm_v(3)%ns

    END SUBROUTINE surface_last_actions

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Routine maps surface data read from file after restart - 1D arrays. 
!------------------------------------------------------------------------------!
    SUBROUTINE surface_restore_elements_1d( surf_target, surf_file,            &
                                            start_index_c,                     &
                                            start_index_on_file,               &
                                            end_index_on_file,                 &
                                            nxlc, nysc, nxlf, nxrf, nysf, nynf,&
                                            nys_on_file, nyn_on_file,          &
                                            nxl_on_file,nxr_on_file )

       IMPLICIT NONE
    
       INTEGER(iwp) ::  i         !< running index along x-direction, refers to former domain size
       INTEGER(iwp) ::  ic        !< running index along x-direction, refers to current domain size
       INTEGER(iwp) ::  j         !< running index along y-direction, refers to former domain size
       INTEGER(iwp) ::  jc        !< running index along y-direction, refers to former domain size        
       INTEGER(iwp) ::  m         !< surface-element index on file
       INTEGER(iwp) ::  mm        !< surface-element index on current subdomain
       INTEGER(iwp) ::  nxlc      !< index of left boundary on current subdomain
       INTEGER(iwp) ::  nxlf      !< index of left boundary on former subdomain 
       INTEGER(iwp) ::  nxrf      !< index of right boundary on former subdomain
       INTEGER(iwp) ::  nysc      !< index of north boundary on current subdomain
       INTEGER(iwp) ::  nynf      !< index of north boundary on former subdomain
       INTEGER(iwp) ::  nysf      !< index of south boundary on former subdomain

       INTEGER(iwp) ::  nxl_on_file !< leftmost index on file
       INTEGER(iwp) ::  nxr_on_file !< rightmost index on file
       INTEGER(iwp) ::  nyn_on_file !< northmost index on file
       INTEGER(iwp) ::  nys_on_file !< southmost index on file

       INTEGER(iwp), DIMENSION(nys:nyn,nxl:nxr) ::  start_index_c             
       INTEGER(iwp), DIMENSION(nys_on_file:nyn_on_file,nxl_on_file:nxr_on_file) :: & 
                            start_index_on_file   !< start index of surface elements on file
       INTEGER(iwp), DIMENSION(nys_on_file:nyn_on_file,nxl_on_file:nxr_on_file) :: & 
                            end_index_on_file     !< end index of surface elements on file
       
       REAL(wp), DIMENSION(:) ::  surf_target !< target surface type
       REAL(wp), DIMENSION(:) ::  surf_file   !< surface type on file

       ic = nxlc
       DO  i = nxlf, nxrf
          jc = nysc
          DO  j = nysf, nynf

             mm = start_index_c(jc,ic)
             DO  m = start_index_on_file(j,i), end_index_on_file(j,i)
                surf_target(mm) = surf_file(m)
                mm = mm + 1
             ENDDO

             jc = jc + 1
          ENDDO
          ic = ic + 1
       ENDDO


    END SUBROUTINE surface_restore_elements_1d
   
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Routine maps surface data read from file after restart - 2D arrays
!------------------------------------------------------------------------------! 
    SUBROUTINE surface_restore_elements_2d( surf_target, surf_file,            &
                                            start_index_c,                     &
                                            start_index_on_file,               &
                                            end_index_on_file,                 &
                                            nxlc, nysc, nxlf, nxrf, nysf, nynf,&
                                            nys_on_file, nyn_on_file,          &
                                            nxl_on_file,nxr_on_file )

       IMPLICIT NONE
    
       INTEGER(iwp) ::  i         !< running index along x-direction, refers to former domain size
       INTEGER(iwp) ::  ic        !< running index along x-direction, refers to current domain size
       INTEGER(iwp) ::  j         !< running index along y-direction, refers to former domain size
       INTEGER(iwp) ::  jc        !< running index along y-direction, refers to former domain size        
       INTEGER(iwp) ::  m         !< surface-element index on file
       INTEGER(iwp) ::  mm        !< surface-element index on current subdomain
       INTEGER(iwp) ::  nxlc      !< index of left boundary on current subdomain
       INTEGER(iwp) ::  nxlf      !< index of left boundary on former subdomain 
       INTEGER(iwp) ::  nxrf      !< index of right boundary on former subdomain
       INTEGER(iwp) ::  nysc      !< index of north boundary on current subdomain
       INTEGER(iwp) ::  nynf      !< index of north boundary on former subdomain
       INTEGER(iwp) ::  nysf      !< index of south boundary on former subdomain

       INTEGER(iwp) ::  nxl_on_file !< leftmost index on file
       INTEGER(iwp) ::  nxr_on_file !< rightmost index on file
       INTEGER(iwp) ::  nyn_on_file !< northmost index on file
       INTEGER(iwp) ::  nys_on_file !< southmost index on file

       INTEGER(iwp), DIMENSION(nys:nyn,nxl:nxr) ::  start_index_c !< start index of surface type
       INTEGER(iwp), DIMENSION(nys_on_file:nyn_on_file,nxl_on_file:nxr_on_file) :: & 
                            start_index_on_file   !< start index of surface elements on file
       INTEGER(iwp), DIMENSION(nys_on_file:nyn_on_file,nxl_on_file:nxr_on_file) :: & 
                            end_index_on_file     !< end index of surface elements on file
       
       REAL(wp), DIMENSION(:,:) ::  surf_target !< target surface type
       REAL(wp), DIMENSION(:,:) ::  surf_file   !< surface type on file
       
       ic = nxlc
       DO  i = nxlf, nxrf
          jc = nysc
          DO  j = nysf, nynf
             mm = start_index_c(jc,ic)
             DO  m = start_index_on_file(j,i), end_index_on_file(j,i)
                surf_target(:,mm) = surf_file(:,m)
                mm = mm + 1
             ENDDO

             jc = jc + 1
          ENDDO
          ic = ic + 1
       ENDDO

    END SUBROUTINE surface_restore_elements_2d


 END MODULE surface_mod
