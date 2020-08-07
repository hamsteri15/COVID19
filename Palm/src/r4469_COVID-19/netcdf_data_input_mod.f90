!> @file netcdf_data_input_mod.f90
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
! $Id: netcdf_data_input_mod.f90 4457 2020-03-11 14:20:43Z raasch $
! use statement for exchange horiz added,
! bugfixes for calls of exchange horiz 2d
! 
! 4435 2020-03-03 10:38:41Z raasch
! temporary bugfix to avoid compile problems with older NetCDFD libraries on IMUK machines
! 
! 4434 2020-03-03 10:02:18Z oliver.maas
! added optional netcdf data input for wtm array input parameters
! 
! 4404 2020-02-12 17:01:53Z suehring
! Fix misplaced preprocessor directives.
! 
! 4401 2020-02-11 16:19:09Z suehring
! Define a default list of coordinate reference system variables used when
! no static driver input is available
! 
! 4400 2020-02-10 20:32:41Z suehring
! - Routine to inquire default fill values added
! - netcdf_data_input_att and netcdf_data_input_var routines removed
! 
! 4392 2020-01-31 16:14:57Z pavelkrc
! (resler) Decrease length of reading buffer (fix problem of ifort/icc compilers)
! 
! 4389 2020-01-29 08:22:42Z raasch
! Error messages refined for reading ASCII topo file, also reading of topo file revised so that
! statement labels and goto statements are not required any more
! 
! 4388 2020-01-28 16:36:55Z raasch
! bugfix for error messages while reading ASCII topo file
! 
! 4387 2020-01-28 11:44:20Z banzhafs
! Added subroutine get_variable_string_generic ( )
! and added to interface get_variable to circumvent
! unknown application-specific restrictions
! in existing function get_variable_string ( ),
! which is retained for backward compatibility (ECC)
!
! 4370 2020-01-10 14:00:44Z raasch
! collective read switched off on NEC Aurora to avoid hang situations
!
! 4362 2020-01-07 17:15:02Z suehring
! Input of plant canopy variables from static driver moved to plant-canopy 
! model
! 
! 4360 2020-01-07 11:25:50Z suehring
! Correct single message calls, local checks must be given by the respective
! mpi rank.
! 
! 4346 2019-12-18 11:55:56Z motisi
! Introduction of wall_flags_total_0, which currently sets bits based on static
! topography information used in wall_flags_static_0
! 
! 4329 2019-12-10 15:46:36Z motisi
! Renamed wall_flags_0 to wall_flags_static_0
! 
! 4321 2019-12-04 10:26:38Z pavelkrc
! Further revise check for surface fractions
! 
! 4313 2019-11-27 14:07:00Z suehring
! Checks for surface fractions revised
! 
! 4312 2019-11-27 14:06:25Z suehring
! Open input files with read-only attribute instead of write attribute.
! 
! 4280 2019-10-29 14:34:15Z monakurppa
! Remove id_emis flags from get_variable_4d_to_3d_real and
! get_variable_5d_to_4d_real
! 
! 4258 2019-10-07 13:29:08Z suehring
! - Migrate input of soil temperature and moisture to land-surface model.
! - Remove interpolate routines and move the only required subroutine to 
!   land-surface model.
! 
! 4247 2019-09-30 10:18:24Z pavelkrc
! Add reading and processing of building_surface_pars
! 
! 4226 2019-09-10 17:03:24Z suehring
! - Netcdf input routine for dimension length renamed
! - Move offline-nesting-specific checks to nesting_offl_mod
! - Module-specific input of boundary data for offline nesting moved to 
!   nesting_offl_mod
! - Define module specific data type for offline nesting in nesting_offl_mod
! 
! 4190 2019-08-27 15:42:37Z suehring
! type real_1d changed to real_1d_3d
! 
! 4186 2019-08-23 16:06:14Z suehring
! Minor formatting adjustments
! 
! 4182 2019-08-22 15:20:23Z scharf
! Corrected "Former revisions" section
! 
! 4178 2019-08-21 11:13:06Z suehring
! Implement input of external radiation forcing. Therefore, provide public
! subroutines and variables. 
! 
! 4150 2019-08-08 20:00:47Z suehring
! Some variables are given the public attribute, in order to call netcdf input
! from single routines
! 
! 4125 2019-07-29 13:31:44Z suehring
! To enable netcdf-parallel access for lateral boundary data (dynamic input), 
! zero number of elements are passed to the respective get_variable routine 
! for non-boundary cores. 
! 
! 4100 2019-07-17 08:11:29Z forkel
! Made check for input_pids_dynamic and 'inifor' more general
! 
! 4012 2019-05-31 15:19:05Z monakurppa
! 
! 3994 2019-05-22 18:08:09Z suehring
! Remove single location message
! 
! 3976 2019-05-15 11:02:34Z hellstea
! Remove unused variables from last commit
! 
! 3969 2019-05-13 12:14:33Z suehring
! - clean-up index notations for emission_values to eliminate magic numbers
! - introduce temporary variable dum_var_5d as well as subroutines
!   get_var_5d_real and get_var_5d_real_dynamic
! - remove emission-specific code in generic get_variable routines
! - in subroutine netcdf_data_input_chemistry_data change netCDF LOD 1 
!   (default) emission_values to the following index order:
!   z, y, x, species, category
! - in subroutine netcdf_data_input_chemistry_data
!   changed netCDF LOD 2 pre-processed emission_values to the following index 
!   order: time, z, y, x, species
! - in type chem_emis_att_type replace nspec with n_emiss_species
!   but retained nspec for backward compatibility with salsa_mod. (E.C. Chan)
! 
! 3961 2019-05-08 16:12:31Z suehring
! Revise checks for building IDs and types
! 
! 3943 2019-05-02 09:50:41Z maronga
! Temporarily disabled some (faulty) checks for static driver.
! 
! 3942 2019-04-30 13:08:30Z kanani
! Fix: increase LEN of all NetCDF attribute values (caused crash in 
! netcdf_create_global_atts due to insufficient length)
! 
! 3941 2019-04-30 09:48:33Z suehring
! Move check for grid dimension to an earlier point in time when first array
! is read.
! Improve checks for building types / IDs with respect to 2D/3D buildings.
! 
! 3885 2019-04-11 11:29:34Z kanani
! Changes related to global restructuring of location messages and introduction 
! of additional debug messages
! 
! 3864 2019-04-05 09:01:56Z monakurppa
! get_variable_4d_to_3d_real modified to enable read in data of type
! data(t,y,x,n) one timestep at a time + some routines made public
! 
! 3855 2019-04-03 10:00:59Z suehring
! Typo removed
! 
! 3854 2019-04-02 16:59:33Z suehring
! Bugfix in one of the checks. Typo removed.
! 
! 3744 2019-02-15 18:38:58Z suehring
! Enable mesoscale offline nesting for chemistry variables as well as 
! initialization of chemistry via dynamic input file.
! 
! 3705 2019-01-29 19:56:39Z suehring
! Interface for attribute input of 8-bit and 32-bit integer 
! 
! 3704 2019-01-29 19:51:41Z suehring
! unused variables removed
! 
! 2696 2017-12-14 17:12:51Z kanani
! Initial revision (suehring)
!
! Authors:
! --------
! @author Matthias Suehring
! @author Edward C. Chan
! @author Emanuele Russo
!
! Description:
! ------------
!> Modulue contains routines to input data according to Palm input data
!> standart using dynamic and static input files.
!> @todo - Chemistry: revise reading of netcdf file and ajdust formatting
!>         according to standard!!! (ecc/done)
!> @todo - Order input alphabetically
!> @todo - Revise error messages and error numbers
!> @todo - Input of missing quantities (chemical species, emission rates)
!> @todo - Defninition and input of still missing variable attributes
!>         (ecc/what are they?)
!> @todo - Input of initial geostrophic wind profiles with cyclic conditions.
!> @todo - remove z dimension from default_emission_data nad preproc_emission_data
!          and correpsonding subroutines get_var_5d_real and get_var_5d_dynamic (ecc)
!> @todo - decpreciate chem_emis_att_type@nspec (ecc)
!> @todo - depreciate subroutines get_variable_4d_to_3d_real and
!>         get_variable_5d_to_4d_real (ecc)
!> @todo - introduce useful debug_message(s)
!------------------------------------------------------------------------------!
 MODULE netcdf_data_input_mod

    USE control_parameters,                                                    &
        ONLY:  coupling_char, io_blocks, io_group

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point_s

    USE indices,                                                               &
        ONLY:  nbgp

    USE kinds

#if defined ( __netcdf )
    USE NETCDF
#endif

    USE pegrid

    USE surface_mod,                                                           &
        ONLY:  ind_pav_green, ind_veg_wall, ind_wat_win
!
!-- Define type for dimensions.
    TYPE dims_xy
       INTEGER(iwp) :: nx                             !< dimension length in x
       INTEGER(iwp) :: ny                             !< dimension length in y
       INTEGER(iwp) :: nz                             !< dimension length in z
       REAL(wp), DIMENSION(:), ALLOCATABLE :: x       !< dimension array in x
       REAL(wp), DIMENSION(:), ALLOCATABLE :: y       !< dimension array in y
       REAL(wp), DIMENSION(:), ALLOCATABLE :: z       !< dimension array in z
    END TYPE dims_xy
    TYPE init_type

       CHARACTER(LEN=16) ::  init_char = 'init_atmosphere_'          !< leading substring for init variables 
       CHARACTER(LEN=23) ::  origin_time = '2000-01-01 00:00:00 +00' !< reference time of input data
       CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE ::  var_names_chem !< list of chemistry variable names that can potentially be on file

       INTEGER(iwp) ::  lod_msoil !< level of detail - soil moisture
       INTEGER(iwp) ::  lod_pt    !< level of detail - pt
       INTEGER(iwp) ::  lod_q     !< level of detail - q
       INTEGER(iwp) ::  lod_tsoil !< level of detail - soil temperature
       INTEGER(iwp) ::  lod_u     !< level of detail - u-component
       INTEGER(iwp) ::  lod_v     !< level of detail - v-component
       INTEGER(iwp) ::  lod_w     !< level of detail - w-component
       INTEGER(iwp) ::  nx        !< number of scalar grid points along x in dynamic input file
       INTEGER(iwp) ::  nxu       !< number of u grid points along x in dynamic input file
       INTEGER(iwp) ::  ny        !< number of scalar grid points along y in dynamic input file
       INTEGER(iwp) ::  nyv       !< number of v grid points along y in dynamic input file
       INTEGER(iwp) ::  nzs       !< number of vertical soil levels in dynamic input file
       INTEGER(iwp) ::  nzu       !< number of vertical levels on scalar grid in dynamic input file
       INTEGER(iwp) ::  nzw       !< number of vertical levels on w grid in dynamic input file
       
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  lod_chem !< level of detail - chemistry variables

       LOGICAL ::  from_file_msoil  = .FALSE. !< flag indicating whether soil moisture is already initialized from file
       LOGICAL ::  from_file_pt     = .FALSE. !< flag indicating whether pt is already initialized from file
       LOGICAL ::  from_file_q      = .FALSE. !< flag indicating whether q is already initialized from file
       LOGICAL ::  from_file_tsoil  = .FALSE. !< flag indicating whether soil temperature is already initialized from file
       LOGICAL ::  from_file_u      = .FALSE. !< flag indicating whether u is already initialized from file
       LOGICAL ::  from_file_ug     = .FALSE. !< flag indicating whether ug is already initialized from file
       LOGICAL ::  from_file_v      = .FALSE. !< flag indicating whether v is already initialized from file
       LOGICAL ::  from_file_vg     = .FALSE. !< flag indicating whether ug is already initialized from file
       LOGICAL ::  from_file_w      = .FALSE. !< flag indicating whether w is already initialized from file
       
       LOGICAL, DIMENSION(:), ALLOCATABLE ::  from_file_chem !< flag indicating whether chemistry variable is read from file

       REAL(wp) ::  fill_msoil              !< fill value for soil moisture
       REAL(wp) ::  fill_pt                 !< fill value for pt
       REAL(wp) ::  fill_q                  !< fill value for q
       REAL(wp) ::  fill_tsoil              !< fill value for soil temperature
       REAL(wp) ::  fill_u                  !< fill value for u
       REAL(wp) ::  fill_v                  !< fill value for v
       REAL(wp) ::  fill_w                  !< fill value for w
       REAL(wp) ::  latitude = 0.0_wp       !< latitude of the lower left corner
       REAL(wp) ::  longitude = 0.0_wp      !< longitude of the lower left corner
       REAL(wp) ::  origin_x = 500000.0_wp  !< UTM easting of the lower left corner
       REAL(wp) ::  origin_y = 0.0_wp       !< UTM northing of the lower left corner
       REAL(wp) ::  origin_z = 0.0_wp       !< reference height of input data
       REAL(wp) ::  rotation_angle = 0.0_wp !< rotation angle of input data

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  fill_chem    !< fill value - chemistry variables
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  msoil_1d     !< initial vertical profile of soil moisture
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  pt_init      !< initial vertical profile of pt
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  q_init       !< initial vertical profile of q
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  tsoil_1d     !< initial vertical profile of soil temperature
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  u_init       !< initial vertical profile of u
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  ug_init      !< initial vertical profile of ug
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  v_init       !< initial vertical profile of v
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  vg_init      !< initial vertical profile of ug
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  w_init       !< initial vertical profile of w
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  z_soil       !< vertical levels in soil in dynamic input file, used for interpolation
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  zu_atmos     !< vertical levels at scalar grid in dynamic input file, used for interpolation
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  zw_atmos     !< vertical levels at w grid in dynamic input file, used for interpolation
       
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  chem_init  !< initial vertical profiles of chemistry variables


       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  msoil_3d !< initial 3d soil moisture provide by Inifor and interpolated onto soil grid
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  tsoil_3d !< initial 3d soil temperature provide by Inifor and interpolated onto soil grid

    END TYPE init_type

!-- Data type for the general information of chemistry emissions, do not dependent on the particular chemical species
    TYPE chem_emis_att_type 

       !-DIMENSIONS
       
       INTEGER(iwp)                                 :: nspec=0            !< no of chem species provided in emission_values
       INTEGER(iwp)                                 :: n_emiss_species=0  !< no of chem species provided in emission_values
                                                                          !< same function as nspec, which will be depreciated (ecc)
                                                                                 
       INTEGER(iwp)                                 :: ncat=0             !< number of emission categories
       INTEGER(iwp)                                 :: nvoc=0             !< number of VOC components
       INTEGER(iwp)                                 :: npm=0              !< number of PM components
       INTEGER(iwp)                                 :: nnox=2             !< number of NOx components: NO and NO2
       INTEGER(iwp)                                 :: nsox=2             !< number of SOX components: SO and SO4
       INTEGER(iwp)                                 :: nhoursyear         !< number of hours of a specific year in the HOURLY mode
                                                                          !< of the default mode
       INTEGER(iwp)                                 :: nmonthdayhour      !< number of month days and hours in the MDH mode 
                                                                          !< of the default mode
       INTEGER(iwp)                                 :: dt_emission        !< Number of emissions timesteps for one year 
                                                                          !< in the pre-processed emissions case
       !-- 1d emission input variables
       CHARACTER (LEN=25),ALLOCATABLE, DIMENSION(:) :: pm_name       !< Names of PM components
       CHARACTER (LEN=25),ALLOCATABLE, DIMENSION(:) :: cat_name      !< Emission category names
       CHARACTER (LEN=25),ALLOCATABLE, DIMENSION(:) :: species_name  !< Names of emission chemical species
       CHARACTER (LEN=25),ALLOCATABLE, DIMENSION(:) :: voc_name      !< Names of VOCs components
       CHARACTER (LEN=25)                           :: units         !< Units

       INTEGER(iwp)                                 :: i_hour         !< indices for assigning emission values at different timesteps
       INTEGER(iwp),ALLOCATABLE, DIMENSION(:)       :: cat_index      !< Indices for emission categories
       INTEGER(iwp),ALLOCATABLE, DIMENSION(:)       :: species_index  !< Indices for emission chem species

       REAL(wp),ALLOCATABLE, DIMENSION(:)           :: xm             !< Molecular masses of emission chem species

       !-- 2d emission input variables 
       REAL(wp),ALLOCATABLE, DIMENSION(:,:)         :: hourly_emis_time_factor  !< Time factors for HOURLY emissions (DEFAULT mode)
       REAL(wp),ALLOCATABLE, DIMENSION(:,:)         :: mdh_emis_time_factor     !< Time factors for MDH emissions (DEFAULT mode)
       REAL(wp),ALLOCATABLE, DIMENSION(:,:)         :: nox_comp                 !< Composition of NO and NO2 
       REAL(wp),ALLOCATABLE, DIMENSION(:,:)         :: sox_comp                 !< Composition of SO2 and SO4
       REAL(wp),ALLOCATABLE, DIMENSION(:,:)         :: voc_comp                 !< Composition of VOC components (not fixed)

       !-- 3d emission input variables
       REAL(wp),ALLOCATABLE, DIMENSION(:,:,:)       :: pm_comp                  !< Composition of PM components (not fixed) 
  
    END TYPE chem_emis_att_type


!-- Data type for the values of chemistry emissions
    TYPE chem_emis_val_type 

       !REAL(wp),ALLOCATABLE, DIMENSION(:,:)     :: stack_height           !< stack height (ecc / to be implemented)
       REAL(wp),ALLOCATABLE, DIMENSION(:,:,:)    :: default_emission_data  !< Emission input values for LOD1 (DEFAULT mode)
       REAL(wp),ALLOCATABLE, DIMENSION(:,:,:,:)  :: preproc_emission_data  !< Emission input values for LOD2 (PRE-PROCESSED mode)

    END TYPE chem_emis_val_type

!
!-- Define data structures for different input data types.
!-- 8-bit Integer 2D
    TYPE int_2d_8bit
       INTEGER(KIND=1) ::  fill = -127                      !< fill value
       INTEGER(KIND=1), DIMENSION(:,:), ALLOCATABLE ::  var !< respective variable

       LOGICAL ::  from_file = .FALSE.  !< flag indicating whether an input variable is available and read from file or default values are used
    END TYPE int_2d_8bit
!
!-- 8-bit Integer 3D
    TYPE int_3d_8bit
       INTEGER(KIND=1) ::  fill = -127                           !< fill value
       INTEGER(KIND=1), DIMENSION(:,:,:), ALLOCATABLE ::  var_3d !< respective variable

       LOGICAL ::  from_file = .FALSE.  !< flag indicating whether an input variable is available and read from file or default values are used
    END TYPE int_3d_8bit
!
!-- 32-bit Integer 2D
    TYPE int_2d_32bit
       INTEGER(iwp) ::  fill = -9999                      !< fill value
       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  var  !< respective variable

       LOGICAL ::  from_file = .FALSE. !< flag indicating whether an input variable is available and read from file or default values are used
    END TYPE int_2d_32bit
!
!-- Define data type to read 1D or 3D real variables. 
    TYPE real_1d_3d
       LOGICAL ::  from_file = .FALSE.  !< flag indicating whether an input variable is available and read from file or default values are used

       INTEGER(iwp) ::  lod = -1        !< level-of-detail
       
       REAL(wp) ::  fill = -9999.9_wp                  !< fill value
       
       REAL(wp), DIMENSION(:),     ALLOCATABLE ::  var1d     !< respective 1D variable
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  var3d     !< respective 3D variable
    END TYPE real_1d_3d   
!
!-- Define data type to read 2D real variables
    TYPE real_2d
       LOGICAL ::  from_file = .FALSE.  !< flag indicating whether an input variable is available and read from file or default values are used

       INTEGER(iwp) ::  lod             !< level-of-detail
       
       REAL(wp) ::  fill = -9999.9_wp                !< fill value
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  var !< respective variable
    END TYPE real_2d

!
!-- Define data type to read 3D real variables
    TYPE real_3d
       LOGICAL ::  from_file = .FALSE.  !< flag indicating whether an input variable is available and read from file or default values are used

       INTEGER(iwp) ::  nz   !< number of grid points along vertical dimension

       REAL(wp) ::  fill = -9999.9_wp                  !< fill value
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  var !< respective variable
    END TYPE real_3d
!
!-- Define data structure where the dimension and type of the input depends
!-- on the given level of detail.
!-- For buildings, the input is either 2D float, or 3d byte.
    TYPE build_in
       INTEGER(iwp)    ::  lod = 1                               !< level of detail
       INTEGER(KIND=1) ::  fill2 = -127                          !< fill value for lod = 2
       INTEGER(iwp)    ::  nz                                    !< number of vertical layers in file
       INTEGER(KIND=1), DIMENSION(:,:,:), ALLOCATABLE ::  var_3d !< 3d variable (lod = 2)

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  z                 !< vertical coordinate for 3D building, used for consistency check

       LOGICAL ::  from_file = .FALSE.  !< flag indicating whether an input variable is available and read from file or default values are used

       REAL(wp)                              ::  fill1 = -9999.9_wp !< fill values for lod = 1
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  var_2d             !< 2d variable (lod = 1)
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  oro_max            !< terraing height under particular buildings
    END TYPE build_in

!
!-- For soil_type, the input is either 2D or 3D one-byte integer.
    TYPE soil_in
       INTEGER(iwp)                                   ::  lod = 1      !< level of detail
       INTEGER(KIND=1)                                ::  fill = -127  !< fill value for lod = 2
       INTEGER(KIND=1), DIMENSION(:,:), ALLOCATABLE   ::  var_2d       !< 2d variable (lod = 1)
       INTEGER(KIND=1), DIMENSION(:,:,:), ALLOCATABLE ::  var_3d       !< 3d variable (lod = 2)

       LOGICAL ::  from_file = .FALSE.  !< flag indicating whether an input variable is available and read from file or default values are used
    END TYPE soil_in

!
!-- Define data type for fractions between surface types
    TYPE fracs
       INTEGER(iwp)                            ::  nf             !< total number of fractions
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  nfracs         !< dimension array for fraction

       LOGICAL ::  from_file = .FALSE. !< flag indicating whether an input variable is available and read from file or default values are used

       REAL(wp)                                ::  fill = -9999.9_wp !< fill value
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  frac              !< respective fraction between different surface types
    END TYPE fracs
!
!-- Data type for parameter lists, Depending on the given level of detail,
!-- the input is 3D or 4D
    TYPE pars
       INTEGER(iwp)                            ::  lod = 1         !< level of detail
       INTEGER(iwp)                            ::  np              !< total number of parameters
       INTEGER(iwp)                            ::  nz              !< vertical dimension - number of soil layers
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  layers          !< dimension array for soil layers
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  pars            !< dimension array for parameters

       LOGICAL ::  from_file = .FALSE.  !< flag indicating whether an input variable is available and read from file or default values are used

       REAL(wp)                                  ::  fill = -9999.9_wp !< fill value
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE   ::  pars_xy           !< respective parameters, level of detail = 1
       REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  pars_xyz          !< respective parameters, level of detail = 2
    END TYPE pars
!
!-- Data type for surface parameter lists
    TYPE pars_surf
       INTEGER(iwp)                                ::  np          !< total number of parameters
       INTEGER(iwp)                                ::  nsurf       !< number of local surfaces
       INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  index_ji    !< index for beginning and end of surfaces at (j,i)
       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE   ::  coords      !< (k,j,i,norm_z,norm_y,norm_x)
                                                                   !< k,j,i:                surface position
                                                                   !< norm_z,norm_y,norm_x: surface normal vector

       LOGICAL ::  from_file = .FALSE.  !< flag indicating whether an input variable is available and read from file or default values are used

       REAL(wp)                              ::  fill = -9999.9_wp !< fill value
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  pars              !< respective parameters per surface
    END TYPE pars_surf
!
!-- Define type for global file attributes
!-- Please refer to the PALM data standard for a detailed description of each
!-- attribute.
    TYPE global_atts_type
       CHARACTER(LEN=200) ::  acronym = ' '                      !< acronym of institution
       CHARACTER(LEN=7)   ::  acronym_char = 'acronym'           !< name of attribute
       CHARACTER(LEN=200) ::  author  = ' '                      !< first name, last name, email adress
       CHARACTER(LEN=6)   ::  author_char = 'author'             !< name of attribute
       CHARACTER(LEN=200) ::  campaign = 'PALM-4U'               !< name of campaign
       CHARACTER(LEN=8)   ::  campaign_char = 'campaign'         !< name of attribute
       CHARACTER(LEN=200) ::  comment = ' '                      !< comment to data
       CHARACTER(LEN=7)   ::  comment_char = 'comment'           !< name of attribute
       CHARACTER(LEN=200) ::  contact_person = ' '               !< first name, last name, email adress
       CHARACTER(LEN=14)  ::  contact_person_char = 'contact_person'  !< name of attribute
       CHARACTER(LEN=200) ::  conventions = 'CF-1.7'             !< netCDF convention
       CHARACTER(LEN=11)  ::  conventions_char = 'Conventions'   !< name of attribute
       CHARACTER(LEN=23 ) ::  creation_time = ' '                !< creation time of data set
       CHARACTER(LEN=13)  ::  creation_time_char = 'creation_time'  !< name of attribute
       CHARACTER(LEN=200) ::  data_content = ' '                 !< content of data set
       CHARACTER(LEN=12)  ::  data_content_char = 'data_content' !< name of attribute
       CHARACTER(LEN=200) ::  dependencies = ' '                 !< dependencies of data set
       CHARACTER(LEN=12)  ::  dependencies_char = 'dependencies' !< name of attribute
       CHARACTER(LEN=200) ::  history = ' '                      !< information about data processing
       CHARACTER(LEN=7)   ::  history_char = 'history'           !< name of attribute
       CHARACTER(LEN=200) ::  institution = ' '                  !< name of responsible institution
       CHARACTER(LEN=11)  ::  institution_char = 'institution'   !< name of attribute
       CHARACTER(LEN=200) ::  keywords = ' '                     !< keywords of data set
       CHARACTER(LEN=8)   ::  keywords_char = 'keywords'         !< name of attribute
       CHARACTER(LEN=200) ::  licence = ' '                      !< licence of data set
       CHARACTER(LEN=7)   ::  licence_char = 'licence'           !< name of attribute
       CHARACTER(LEN=200) ::  location = ' '                     !< place which refers to data set
       CHARACTER(LEN=8)   ::  location_char = 'location'         !< name of attribute
       CHARACTER(LEN=10)  ::  origin_lat_char = 'origin_lat'     !< name of attribute
       CHARACTER(LEN=10)  ::  origin_lon_char = 'origin_lon'     !< name of attribute
       CHARACTER(LEN=23 ) ::  origin_time = '2000-01-01 00:00:00 +00'  !< reference time
       CHARACTER(LEN=11)  ::  origin_time_char = 'origin_time'   !< name of attribute
       CHARACTER(LEN=8)   ::  origin_x_char = 'origin_x'         !< name of attribute
       CHARACTER(LEN=8)   ::  origin_y_char = 'origin_y'         !< name of attribute
       CHARACTER(LEN=8)   ::  origin_z_char = 'origin_z'         !< name of attribute
       CHARACTER(LEN=12)  ::  palm_version_char = 'palm_version' !< name of attribute
       CHARACTER(LEN=200) ::  references = ' '                   !< literature referring to data set
       CHARACTER(LEN=10)  ::  references_char = 'references'     !< name of attribute
       CHARACTER(LEN=14)  ::  rotation_angle_char = 'rotation_angle'  !< name of attribute
       CHARACTER(LEN=200) ::  site = ' '                         !< name of model domain
       CHARACTER(LEN=4)   ::  site_char = 'site'                 !< name of attribute
       CHARACTER(LEN=200) ::  source = ' '                       !< source of data set
       CHARACTER(LEN=6)   ::  source_char = 'source'             !< name of attribute
       CHARACTER(LEN=200) ::  title = ' '                        !< title of data set
       CHARACTER(LEN=5)   ::  title_char = 'title'               !< name of attribute
       CHARACTER(LEN=7)   ::  version_char = 'version'           !< name of attribute

       INTEGER(iwp) ::  version              !< version of data set

       REAL(wp) ::  fillvalue = -9999.0      !< default fill value
       REAL(wp) ::  origin_lat               !< latitude of lower left corner
       REAL(wp) ::  origin_lon               !< longitude of lower left corner
       REAL(wp) ::  origin_x                 !< easting (UTM coordinate) of lower left corner
       REAL(wp) ::  origin_y                 !< northing (UTM coordinate) of lower left corner
       REAL(wp) ::  origin_z                 !< reference height
       REAL(wp) ::  palm_version             !< PALM version of data set
       REAL(wp) ::  rotation_angle           !< rotation angle of coordinate system of data set
    END TYPE global_atts_type
!
!-- Define type for coordinate reference system (crs)
    TYPE crs_type
       CHARACTER(LEN=200) ::  epsg_code = 'EPSG:25831'                   !< EPSG code
       CHARACTER(LEN=200) ::  grid_mapping_name = 'transverse_mercator'  !< name of grid mapping
       CHARACTER(LEN=200) ::  long_name = 'coordinate reference system'  !< name of variable crs
       CHARACTER(LEN=200) ::  units = 'm'                                !< unit of crs

       REAL(wp) ::  false_easting = 500000.0_wp                  !< false easting
       REAL(wp) ::  false_northing = 0.0_wp                      !< false northing
       REAL(wp) ::  inverse_flattening = 298.257223563_wp        !< 1/f (default for WGS84)
       REAL(wp) ::  latitude_of_projection_origin = 0.0_wp       !< latitude of projection origin
       REAL(wp) ::  longitude_of_central_meridian = 3.0_wp       !< longitude of central meridian of UTM zone (default: zone 31)
       REAL(wp) ::  longitude_of_prime_meridian = 0.0_wp         !< longitude of prime meridian
       REAL(wp) ::  scale_factor_at_central_meridian = 0.9996_wp !< scale factor of UTM coordinates
       REAL(wp) ::  semi_major_axis = 6378137.0_wp               !< length of semi major axis (default for WGS84)
    END TYPE crs_type

!
!-- Define variables
    TYPE(crs_type)   ::  coord_ref_sys  !< coordinate reference system

    TYPE(dims_xy)    ::  dim_static     !< data structure for x, y-dimension in static input file 

    TYPE(init_type) ::  init_3d    !< data structure for the initialization of the 3D flow and soil fields
    TYPE(init_type) ::  init_model !< data structure for the initialization of the model

!
!-- Define 2D variables of type NC_BYTE
    TYPE(int_2d_8bit)  ::  albedo_type_f     !< input variable for albedo type
    TYPE(int_2d_8bit)  ::  building_type_f   !< input variable for building type
    TYPE(int_2d_8bit)  ::  pavement_type_f   !< input variable for pavenment type
    TYPE(int_2d_8bit)  ::  street_crossing_f !< input variable for water type
    TYPE(int_2d_8bit)  ::  street_type_f     !< input variable for water type
    TYPE(int_2d_8bit)  ::  vegetation_type_f !< input variable for vegetation type
    TYPE(int_2d_8bit)  ::  water_type_f      !< input variable for water type
!
!-- Define 3D variables of type NC_BYTE
    TYPE(int_3d_8bit)  ::  building_obstruction_f    !< input variable for building obstruction
    TYPE(int_3d_8bit)  ::  building_obstruction_full !< input variable for building obstruction
!
!-- Define 2D variables of type NC_INT
    TYPE(int_2d_32bit) ::  building_id_f     !< input variable for building ID
!
!-- Define 2D variables of type NC_FLOAT
    TYPE(real_2d) ::  terrain_height_f       !< input variable for terrain height
    TYPE(real_2d) ::  uvem_irradiance_f      !< input variable for uvem irradiance lookup table
    TYPE(real_2d) ::  uvem_integration_f     !< input variable for uvem integration
!
!-- Define 3D variables of type NC_FLOAT
    TYPE(real_3d) ::  root_area_density_lsm_f !< input variable for root area density - parametrized vegetation
    TYPE(real_3d) ::  uvem_radiance_f         !< input variable for uvem radiance lookup table
    TYPE(real_3d) ::  uvem_projarea_f         !< input variable for uvem projection area lookup table
!
!-- Define input variable for buildings
    TYPE(build_in) ::  buildings_f           !< input variable for buildings
!
!-- Define input variables for soil_type
    TYPE(soil_in)  ::  soil_type_f           !< input variable for soil type

    TYPE(fracs) ::  surface_fraction_f       !< input variable for surface fraction

    TYPE(pars)  ::  albedo_pars_f              !< input variable for albedo parameters
    TYPE(pars)  ::  building_pars_f            !< input variable for building parameters
    TYPE(pars)  ::  pavement_pars_f            !< input variable for pavement parameters
    TYPE(pars)  ::  pavement_subsurface_pars_f !< input variable for pavement parameters
    TYPE(pars)  ::  soil_pars_f                !< input variable for soil parameters
    TYPE(pars)  ::  vegetation_pars_f          !< input variable for vegetation parameters
    TYPE(pars)  ::  water_pars_f               !< input variable for water parameters

    TYPE(pars_surf)  ::  building_surface_pars_f  !< input variable for building surface parameters

    TYPE(chem_emis_att_type)                             ::  chem_emis_att    !< Input Information of Chemistry Emission Data from netcdf  
    TYPE(chem_emis_val_type), ALLOCATABLE, DIMENSION(:)  ::  chem_emis        !< Input Chemistry Emission Data from netcdf  

    CHARACTER(LEN=3)  ::  char_lod  = 'lod'         !< name of level-of-detail attribute in NetCDF file

    CHARACTER(LEN=10) ::  char_fill = '_FillValue'        !< name of fill value attribute in NetCDF file

    CHARACTER(LEN=100) ::  input_file_static  = 'PIDS_STATIC'  !< Name of file which comprises static input data
    CHARACTER(LEN=100) ::  input_file_dynamic = 'PIDS_DYNAMIC' !< Name of file which comprises dynamic input data
    CHARACTER(LEN=100) ::  input_file_chem    = 'PIDS_CHEM'    !< Name of file which comprises chemistry input data
    CHARACTER(LEN=100) ::  input_file_uvem    = 'PIDS_UVEM'    !< Name of file which comprises static uv_exposure model input data
    CHARACTER(LEN=100) ::  input_file_vm      = 'PIDS_VM'      !< Name of file which comprises virtual measurement data
    CHARACTER(LEN=100) ::  input_file_wtm     = 'PIDS_WTM'     !< Name of file which comprises wind turbine model input data
        
    CHARACTER(LEN=25), ALLOCATABLE, DIMENSION(:) ::  string_values  !< output of string variables read from netcdf input files
    CHARACTER(LEN=50), DIMENSION(:), ALLOCATABLE ::  vars_pids      !< variable in input file

    INTEGER(iwp)                                     ::  id_emis        !< NetCDF id of input file for chemistry emissions: TBD: It has to be removed

    INTEGER(iwp) ::  nc_stat         !< return value of nf90 function call
    INTEGER(iwp) ::  num_var_pids    !< number of variables in file
    INTEGER(iwp) ::  pids_id         !< file id

    LOGICAL ::  input_pids_static  = .FALSE.   !< Flag indicating whether Palm-input-data-standard file containing static information exists
    LOGICAL ::  input_pids_dynamic = .FALSE.   !< Flag indicating whether Palm-input-data-standard file containing dynamic information exists
    LOGICAL ::  input_pids_chem    = .FALSE.   !< Flag indicating whether Palm-input-data-standard file containing chemistry information exists
    LOGICAL ::  input_pids_uvem    = .FALSE.   !< Flag indicating whether uv-expoure-model input file containing static information exists
    LOGICAL ::  input_pids_vm      = .FALSE.   !< Flag indicating whether input file for virtual measurements exist
    LOGICAL ::  input_pids_wtm     = .FALSE.   !< Flag indicating whether input file for wind turbine model exists

    LOGICAL ::  collective_read = .FALSE.      !< Enable NetCDF collective read

    REAL(wp), DIMENSION(8) ::  crs_list        !< list of coord_ref_sys values

    TYPE(global_atts_type) ::  input_file_atts !< global attributes of input file

    SAVE

    PRIVATE

    INTERFACE netcdf_data_input_check_dynamic
       MODULE PROCEDURE netcdf_data_input_check_dynamic
    END INTERFACE netcdf_data_input_check_dynamic

    INTERFACE netcdf_data_input_check_static
       MODULE PROCEDURE netcdf_data_input_check_static
    END INTERFACE netcdf_data_input_check_static

    INTERFACE netcdf_data_input_chemistry_data
       MODULE PROCEDURE netcdf_data_input_chemistry_data
    END INTERFACE netcdf_data_input_chemistry_data

    INTERFACE get_dimension_length
       MODULE PROCEDURE get_dimension_length
    END INTERFACE get_dimension_length

    INTERFACE inquire_fill_value
       MODULE PROCEDURE inquire_fill_value_int
       MODULE PROCEDURE inquire_fill_value_real
    END INTERFACE inquire_fill_value

    INTERFACE netcdf_data_input_inquire_file
       MODULE PROCEDURE netcdf_data_input_inquire_file
    END INTERFACE netcdf_data_input_inquire_file

    INTERFACE netcdf_data_input_init
       MODULE PROCEDURE netcdf_data_input_init
    END INTERFACE netcdf_data_input_init

    INTERFACE netcdf_data_input_init_3d
       MODULE PROCEDURE netcdf_data_input_init_3d
    END INTERFACE netcdf_data_input_init_3d
    
    INTERFACE netcdf_data_input_surface_data
       MODULE PROCEDURE netcdf_data_input_surface_data
    END INTERFACE netcdf_data_input_surface_data

    INTERFACE netcdf_data_input_uvem
       MODULE PROCEDURE netcdf_data_input_uvem
    END INTERFACE netcdf_data_input_uvem

    INTERFACE get_variable
       MODULE PROCEDURE get_variable_1d_char
       MODULE PROCEDURE get_variable_1d_int
       MODULE PROCEDURE get_variable_1d_real
       MODULE PROCEDURE get_variable_2d_int8
       MODULE PROCEDURE get_variable_2d_int32
       MODULE PROCEDURE get_variable_2d_real
       MODULE PROCEDURE get_variable_3d_int8
       MODULE PROCEDURE get_variable_3d_real
       MODULE PROCEDURE get_variable_3d_real_dynamic
       MODULE PROCEDURE get_variable_4d_to_3d_real
       MODULE PROCEDURE get_variable_4d_real
       MODULE PROCEDURE get_variable_5d_to_4d_real
       MODULE PROCEDURE get_variable_5d_real           ! (ecc) temp subroutine 4 reading 5D NC arrays
       MODULE PROCEDURE get_variable_5d_real_dynamic   ! 2B removed as z is out of emission_values
       MODULE PROCEDURE get_variable_string
       MODULE PROCEDURE get_variable_string_generic    ! (ecc) generic string function

    END INTERFACE get_variable

    INTERFACE get_variable_pr
       MODULE PROCEDURE get_variable_pr
    END INTERFACE get_variable_pr

    INTERFACE get_attribute
       MODULE PROCEDURE get_attribute_real
       MODULE PROCEDURE get_attribute_int8
       MODULE PROCEDURE get_attribute_int32
       MODULE PROCEDURE get_attribute_string
    END INTERFACE get_attribute

!
!-- Public data structures
    PUBLIC real_1d_3d,                                                         &
           real_2d,                                                            &
           real_3d
!
!-- Public variables
    PUBLIC albedo_pars_f, albedo_type_f, buildings_f,                          &
           building_id_f, building_pars_f, building_surface_pars_f,            &
           building_type_f,                                                    &
           char_fill,                                                          &
           char_lod,                                                           &
           chem_emis, chem_emis_att, chem_emis_att_type, chem_emis_val_type,   &
           coord_ref_sys,                                                      &
           crs_list,                                                           &
           init_3d, init_model, input_file_atts,                               &
           input_file_dynamic,                                                 &
           input_file_static,                                                  &
           input_pids_static,                                                  &
           input_pids_dynamic, input_pids_vm, input_file_vm,                   &
           input_pids_wtm, input_file_wtm,                                     &
           num_var_pids,                                                       &
           pavement_pars_f, pavement_subsurface_pars_f, pavement_type_f,       &
           pids_id,                                                            &
           root_area_density_lsm_f, soil_pars_f,                               &
           soil_type_f, street_crossing_f, street_type_f, surface_fraction_f,  &
           terrain_height_f, vegetation_pars_f, vegetation_type_f,             &
           vars_pids,                                                          &
           water_pars_f, water_type_f
!
!-- Public uv exposure variables
    PUBLIC building_obstruction_f, input_file_uvem, input_pids_uvem,           &
           netcdf_data_input_uvem,                                             &
           uvem_integration_f, uvem_irradiance_f,                              &
           uvem_projarea_f, uvem_radiance_f

!
!-- Public subroutines
    PUBLIC netcdf_data_input_check_dynamic,                                    &
           netcdf_data_input_check_static,                                     &
           netcdf_data_input_chemistry_data,                                   &
           get_dimension_length,                                               &
           netcdf_data_input_inquire_file,                                     &
           netcdf_data_input_init,                                             &
           netcdf_data_input_init_3d,                                          &
           netcdf_data_input_surface_data,                                     &
           netcdf_data_input_topo,                                             &
           get_attribute,                                                      &
           get_variable,                                                       &
           get_variable_pr,                                                    &
           open_read_file,                                                     &
           check_existence,                                                    &
           inquire_fill_value,                                                 &
           inquire_num_variables,                                              &
           inquire_variable_names,                                             &
           close_input_file


 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Inquires whether NetCDF input files according to Palm-input-data standard
!> exist. Moreover, basic checks are performed.
!------------------------------------------------------------------------------!
    SUBROUTINE netcdf_data_input_inquire_file

       USE control_parameters,                                                 &
           ONLY:  topo_no_distinct

       IMPLICIT NONE

#if defined ( __netcdf )
       INQUIRE( FILE = TRIM( input_file_static )   // TRIM( coupling_char ),   &
                EXIST = input_pids_static  )
       INQUIRE( FILE = TRIM( input_file_dynamic ) // TRIM( coupling_char ),    &
                EXIST = input_pids_dynamic )
       INQUIRE( FILE = TRIM( input_file_chem )    // TRIM( coupling_char ),    &
                EXIST = input_pids_chem )
       INQUIRE( FILE = TRIM( input_file_uvem )    // TRIM( coupling_char ),    &
                EXIST = input_pids_uvem  )
       INQUIRE( FILE = TRIM( input_file_vm )      // TRIM( coupling_char ),    &
                EXIST = input_pids_vm )
       INQUIRE( FILE = TRIM( input_file_wtm )     // TRIM( coupling_char ),    &
                EXIST = input_pids_wtm )
#endif

!
!--    As long as topography can be input via ASCII format, no distinction
!--    between building and terrain can be made. This case, classify all
!--    surfaces as default type. Same in case land-surface and urban-surface
!--    model are not applied.
       IF ( .NOT. input_pids_static )  THEN
          topo_no_distinct = .TRUE.
       ENDIF

    END SUBROUTINE netcdf_data_input_inquire_file

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads global attributes and coordinate reference system required for 
!> initialization of the model.
!------------------------------------------------------------------------------!
    SUBROUTINE netcdf_data_input_init

       IMPLICIT NONE

       INTEGER(iwp) ::  id_mod     !< NetCDF id of input file
       INTEGER(iwp) ::  var_id_crs !< NetCDF id of variable crs

!
!--    Define default list of crs attributes. This is required for coordinate
!--    transformation.
       crs_list = (/ coord_ref_sys%semi_major_axis,                            &
                     coord_ref_sys%inverse_flattening,                         &
                     coord_ref_sys%longitude_of_prime_meridian,                &
                     coord_ref_sys%longitude_of_central_meridian,              &
                     coord_ref_sys%scale_factor_at_central_meridian,           &
                     coord_ref_sys%latitude_of_projection_origin,              &
                     coord_ref_sys%false_easting,                              &
                     coord_ref_sys%false_northing /)

       IF ( .NOT. input_pids_static )  RETURN

#if defined ( __netcdf )
!
!--    Open file in read-only mode
       CALL open_read_file( TRIM( input_file_static ) //                       &
                            TRIM( coupling_char ), id_mod )
!
!--    Read global attributes
       CALL get_attribute( id_mod, input_file_atts%origin_lat_char,            &
                           input_file_atts%origin_lat, .TRUE. )

       CALL get_attribute( id_mod, input_file_atts%origin_lon_char,            &
                           input_file_atts%origin_lon, .TRUE. )

       CALL get_attribute( id_mod, input_file_atts%origin_time_char,           &
                           input_file_atts%origin_time, .TRUE. )

       CALL get_attribute( id_mod, input_file_atts%origin_x_char,              &
                           input_file_atts%origin_x, .TRUE. )

       CALL get_attribute( id_mod, input_file_atts%origin_y_char,              &
                           input_file_atts%origin_y, .TRUE. )

       CALL get_attribute( id_mod, input_file_atts%origin_z_char,              &
                           input_file_atts%origin_z, .TRUE. )

       CALL get_attribute( id_mod, input_file_atts%rotation_angle_char,        &
                           input_file_atts%rotation_angle, .TRUE. )

       CALL get_attribute( id_mod, input_file_atts%author_char,                &
                           input_file_atts%author, .TRUE., no_abort=.FALSE. )
       CALL get_attribute( id_mod, input_file_atts%contact_person_char,        &
                           input_file_atts%contact_person, .TRUE., no_abort=.FALSE. )
       CALL get_attribute( id_mod, input_file_atts%institution_char,           &
                           input_file_atts%institution,    .TRUE., no_abort=.FALSE. )
       CALL get_attribute( id_mod, input_file_atts%acronym_char,               &
                           input_file_atts%acronym,        .TRUE., no_abort=.FALSE. )

       CALL get_attribute( id_mod, input_file_atts%campaign_char,              &
                           input_file_atts%campaign, .TRUE., no_abort=.FALSE. )
       CALL get_attribute( id_mod, input_file_atts%location_char,              &
                           input_file_atts%location, .TRUE., no_abort=.FALSE. )
       CALL get_attribute( id_mod, input_file_atts%site_char,                  &
                           input_file_atts%site,     .TRUE., no_abort=.FALSE. )

       CALL get_attribute( id_mod, input_file_atts%source_char,                &
                           input_file_atts%source,     .TRUE., no_abort=.FALSE. )
       CALL get_attribute( id_mod, input_file_atts%references_char,            &
                           input_file_atts%references, .TRUE., no_abort=.FALSE. )
       CALL get_attribute( id_mod, input_file_atts%keywords_char,              &
                           input_file_atts%keywords,   .TRUE., no_abort=.FALSE. )
       CALL get_attribute( id_mod, input_file_atts%licence_char,               &
                           input_file_atts%licence,    .TRUE., no_abort=.FALSE. )
       CALL get_attribute( id_mod, input_file_atts%comment_char,               &
                           input_file_atts%comment,    .TRUE., no_abort=.FALSE. )
!
!--    Read coordinate reference system if available
       nc_stat = NF90_INQ_VARID( id_mod, 'crs', var_id_crs )
       IF ( nc_stat == NF90_NOERR )  THEN
          CALL get_attribute( id_mod, 'epsg_code',                             &
                              coord_ref_sys%epsg_code,                         &
                              .FALSE., 'crs' )
          CALL get_attribute( id_mod, 'false_easting',                         &
                              coord_ref_sys%false_easting,                     &
                              .FALSE., 'crs' )
          CALL get_attribute( id_mod, 'false_northing',                        &
                              coord_ref_sys%false_northing,                    &
                              .FALSE., 'crs' )
          CALL get_attribute( id_mod, 'grid_mapping_name',                     &
                              coord_ref_sys%grid_mapping_name,                 &
                              .FALSE., 'crs' )
          CALL get_attribute( id_mod, 'inverse_flattening',                    &
                              coord_ref_sys%inverse_flattening,                &
                              .FALSE., 'crs' )
          CALL get_attribute( id_mod, 'latitude_of_projection_origin',         &
                              coord_ref_sys%latitude_of_projection_origin,     &
                              .FALSE., 'crs' )
          CALL get_attribute( id_mod, 'long_name',                             &
                              coord_ref_sys%long_name,                         &
                              .FALSE., 'crs' )
          CALL get_attribute( id_mod, 'longitude_of_central_meridian',         &
                              coord_ref_sys%longitude_of_central_meridian,     &
                              .FALSE., 'crs' )
          CALL get_attribute( id_mod, 'longitude_of_prime_meridian',           &
                              coord_ref_sys%longitude_of_prime_meridian,       &
                              .FALSE., 'crs' )
          CALL get_attribute( id_mod, 'scale_factor_at_central_meridian',      &
                              coord_ref_sys%scale_factor_at_central_meridian,  &
                              .FALSE., 'crs' )
          CALL get_attribute( id_mod, 'semi_major_axis',                       &
                              coord_ref_sys%semi_major_axis,                   &
                              .FALSE., 'crs' )
          CALL get_attribute( id_mod, 'units',                                 &
                              coord_ref_sys%units,                             &
                              .FALSE., 'crs' )
       ELSE
!
!--       Calculate central meridian from origin_lon
          coord_ref_sys%longitude_of_central_meridian = &
             CEILING( input_file_atts%origin_lon / 6.0_wp ) * 6.0_wp - 3.0_wp
       ENDIF
!
!--    Finally, close input file
       CALL close_input_file( id_mod )
#endif
!
!--    Copy latitude, longitude, origin_z, rotation angle on init type
       init_model%latitude        = input_file_atts%origin_lat
       init_model%longitude       = input_file_atts%origin_lon
       init_model%origin_time     = input_file_atts%origin_time  
       init_model%origin_x        = input_file_atts%origin_x
       init_model%origin_y        = input_file_atts%origin_y
       init_model%origin_z        = input_file_atts%origin_z  
       init_model%rotation_angle  = input_file_atts%rotation_angle  

!
!--    Update list of crs attributes. This is required for coordinate 
!--    transformation. 
       crs_list = (/ coord_ref_sys%semi_major_axis,                            &
                     coord_ref_sys%inverse_flattening,                         &
                     coord_ref_sys%longitude_of_prime_meridian,                &
                     coord_ref_sys%longitude_of_central_meridian,              &
                     coord_ref_sys%scale_factor_at_central_meridian,           &
                     coord_ref_sys%latitude_of_projection_origin,              &
                     coord_ref_sys%false_easting,                              &
                     coord_ref_sys%false_northing /)
!
!--    In case of nested runs, each model domain might have different longitude
!--    and latitude, which would result in different Coriolis parameters and
!--    sun-zenith angles. To avoid this, longitude and latitude in each model
!--    domain will be set to the values of the root model. Please note, this
!--    synchronization is required already here. 
#if defined( __parallel )
       CALL MPI_BCAST( init_model%latitude,  1, MPI_REAL, 0,                   &
                       MPI_COMM_WORLD, ierr )
       CALL MPI_BCAST( init_model%longitude, 1, MPI_REAL, 0,                   &
                       MPI_COMM_WORLD, ierr )
#endif

    END SUBROUTINE netcdf_data_input_init


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads Chemistry NETCDF Input data, such as emission values, emission species, etc.
!------------------------------------------------------------------------------!

    SUBROUTINE netcdf_data_input_chemistry_data(emt_att,emt)

       USE chem_modules,                                       &
           ONLY:  emiss_lod, time_fac_type, surface_csflux_name

       USE control_parameters,                                 &
           ONLY:  message_string

       USE indices,                                            &
           ONLY:  nxl, nxr, nys, nyn

       IMPLICIT NONE

       TYPE(chem_emis_att_type), INTENT(INOUT)                             ::  emt_att
       TYPE(chem_emis_val_type), ALLOCATABLE, DIMENSION(:), INTENT(INOUT)  ::  emt
    
       INTEGER(iwp)  ::  i, j, k      !< generic counters
       INTEGER(iwp)  ::  ispec        !< index for number of emission species in input
       INTEGER(iwp)  ::  len_dims     !< Length of dimension
       INTEGER(iwp)  ::  num_vars     !< number of variables in netcdf input file

!
!-- dum_var_4d are designed to read in emission_values from the chemistry netCDF file.
!-- Currently the vestigial "z" dimension in emission_values makes it a 5D array,
!-- hence the corresponding dum_var_5d array.  When the "z" dimension is removed
!-- completely, dum_var_4d will be used instead
!-- (ecc 20190425) 

!       REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:)    ::  dum_var_4d  !< temp array 4 4D chem emission data
       REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:,:)  ::  dum_var_5d  !< temp array 4 5D chem emission data

!
!-- Start processing data
!
!-- Emission LOD 0 (Parameterized mode)

        IF  ( emiss_lod == 0 )  THEN

! for reference (ecc)
!       IF (TRIM(mode_emis) == "PARAMETERIZED" .OR. TRIM(mode_emis) == "parameterized") THEN 

           ispec=1
           emt_att%n_emiss_species = 0

!
!-- number of species

           DO  WHILE (TRIM( surface_csflux_name( ispec ) ) /= 'novalue' )

             emt_att%n_emiss_species = emt_att%n_emiss_species + 1
             ispec=ispec+1
!
!-- followling line retained for compatibility with salsa_mod
!-- which still uses emt_att%nspec heavily (ecc)

             emt_att%nspec = emt_att%nspec + 1

           ENDDO

!
!-- allocate emission values data type arrays

          ALLOCATE ( emt(emt_att%n_emiss_species) )

!
!-- Read EMISSION SPECIES NAMES

!
!-- allocate space for strings

          ALLOCATE (emt_att%species_name(emt_att%n_emiss_species) )
 
         DO ispec = 1, emt_att%n_emiss_species
            emt_att%species_name(ispec) = TRIM(surface_csflux_name(ispec))
         ENDDO

!
!-- LOD 1 (default mode) and LOD 2 (pre-processed mode)

       ELSE

#if defined ( __netcdf )

          IF ( .NOT. input_pids_chem )  RETURN

!
!-- first we allocate memory space for the emission species and then
!-- we differentiate between LOD 1 (default mode) and LOD 2 (pre-processed mode)

!
!-- open emission data file ( {palmcase}_chemistry )

          CALL open_read_file ( TRIM(input_file_chem) // TRIM(coupling_char), id_emis )

!
!-- inquire number of variables

          CALL inquire_num_variables ( id_emis, num_vars )

!
!-- Get General Dimension Lengths: only # species and # categories.
!-- Tther dimensions depend on the emission mode or specific components

          CALL get_dimension_length ( id_emis, emt_att%n_emiss_species, 'nspecies' )

!
!-- backward compatibility for salsa_mod (ecc)

          emt_att%nspec = emt_att%n_emiss_species

! 
!-- Allocate emission values data type arrays

          ALLOCATE ( emt(emt_att%n_emiss_species) )

!
!-- READING IN SPECIES NAMES

!
!-- Allocate memory for species names

          ALLOCATE ( emt_att%species_name(emt_att%n_emiss_species) )

!
!-- Retrieve variable name (again, should use n_emiss_strlen)

          CALL get_variable( id_emis, 'emission_name',    &
                             string_values, emt_att%n_emiss_species )
          emt_att%species_name=string_values

!
!-- dealocate string_values previously allocated in get_variable call

          IF  ( ALLOCATED(string_values) )  DEALLOCATE (string_values)

!
!-- READING IN SPECIES INDICES

!
!-- Allocate memory for species indices

          ALLOCATE ( emt_att%species_index(emt_att%n_emiss_species) )

!
!-- Retrieve variable data

          CALL get_variable( id_emis, 'emission_index', emt_att%species_index )
!
!-- Now the routine has to distinguish between chemistry emission
!-- LOD 1 (DEFAULT mode) and LOD 2 (PRE-PROCESSED mode) 

!
!-- START OF EMISSION LOD 1 (DEFAULT MODE)


          IF  ( emiss_lod == 1 )  THEN

! for reference (ecc)
!          IF (TRIM(mode_emis) == "DEFAULT" .OR. TRIM(mode_emis) == "default") THEN

! 
!-- get number of emission categories

             CALL get_dimension_length ( id_emis, emt_att%ncat, 'ncat' )

!-- READING IN EMISSION CATEGORIES INDICES

             ALLOCATE ( emt_att%cat_index(emt_att%ncat) )

!
!-- Retrieve variable data

             CALL get_variable( id_emis, 'emission_cat_index', emt_att%cat_index )


!
!-- Loop through individual species to get basic information on
!-- VOC/PM/NOX/SOX

!------------------------------------------------------------------------------
!-- NOTE - CHECK ARRAY INDICES FOR READING IN NAMES AND SPECIES
!--        IN LOD1 (DEFAULT MODE) FOR THE VARIOUS MODE SPLITS
!--        AS ALL ID_EMIS CONDITIONALS HAVE BEEN REMOVED FROM GET_VAR
!--        FUNCTIONS.  IN THEORY THIS WOULD MEAN ALL ARRAYS SHOULD BE
!--        READ FROM 0 to N-1 (C CONVENTION) AS OPPOSED TO 1 to N
!--        (FORTRAN CONVENTION).  KEEP THIS IN MIND !!
!--        (ecc 20190424)
!------------------------------------------------------------------------------
 
             DO  ispec = 1, emt_att%n_emiss_species

!
!-- VOC DATA (name and composition)

                IF  ( TRIM(emt_att%species_name(ispec)) == "VOC" .OR.                  &
                      TRIM(emt_att%species_name(ispec)) == "voc" )  THEN

!
!-- VOC name
                   CALL get_dimension_length ( id_emis, emt_att%nvoc, 'nvoc' )
                   ALLOCATE ( emt_att%voc_name(emt_att%nvoc) )
                   CALL get_variable ( id_emis,"emission_voc_name",  &
                                       string_values, emt_att%nvoc )
                   emt_att%voc_name = string_values
                   IF  ( ALLOCATED(string_values) )  DEALLOCATE (string_values)

! 
!-- VOC composition

                   ALLOCATE ( emt_att%voc_comp(emt_att%ncat,emt_att%nvoc) )
                   CALL get_variable ( id_emis, "composition_voc", emt_att%voc_comp,     &
                                       1, emt_att%ncat, 1, emt_att%nvoc )

                ENDIF  ! VOC

!
!-- PM DATA (name and composition)

                IF  ( TRIM(emt_att%species_name(ispec)) == "PM" .OR.                   &
                      TRIM(emt_att%species_name(ispec)) == "pm")  THEN

!
!-- PM name

                   CALL get_dimension_length ( id_emis, emt_att%npm, 'npm' )
                   ALLOCATE ( emt_att%pm_name(emt_att%npm) )
                   CALL get_variable ( id_emis, "pm_name", string_values, emt_att%npm )
                   emt_att%pm_name = string_values
                   IF  ( ALLOCATED(string_values) )  DEALLOCATE (string_values)      

!
!-- PM composition (PM1, PM2.5 and PM10)

                   len_dims = 3  ! PM1, PM2.5, PM10
                   ALLOCATE(emt_att%pm_comp(emt_att%ncat,emt_att%npm,len_dims))
                   CALL get_variable ( id_emis, "composition_pm", emt_att%pm_comp,       &
                                       1, emt_att%ncat, 1, emt_att%npm, 1, len_dims )

                ENDIF  ! PM

!
!-- NOX (NO and NO2)

                IF  ( TRIM(emt_att%species_name(ispec)) == "NOX" .OR.                  &
                      TRIM(emt_att%species_name(ispec)) == "nox" )  THEN

                   ALLOCATE ( emt_att%nox_comp(emt_att%ncat,emt_att%nnox) )
                   CALL get_variable ( id_emis, "composition_nox", emt_att%nox_comp,     &
                                       1, emt_att%ncat, 1, emt_att%nnox )

                ENDIF  ! NOX

!
!-- SOX (SO2 and SO4)

                IF  ( TRIM(emt_att%species_name(ispec)) == "SOX" .OR.                  &
                      TRIM(emt_att%species_name(ispec)) == "sox" )  THEN

                   ALLOCATE ( emt_att%sox_comp(emt_att%ncat,emt_att%nsox) )
                   CALL get_variable ( id_emis, "composition_sox", emt_att%sox_comp,     &
                                       1, emt_att%ncat, 1, emt_att%nsox )

                ENDIF  ! SOX

             ENDDO  ! do ispec

!
!-- EMISSION TIME SCALING FACTORS (hourly and MDH data)
 
!     
!-- HOUR   
             IF  ( TRIM(time_fac_type) == "HOUR" .OR.                        &
                   TRIM(time_fac_type) == "hour" )  THEN

                CALL get_dimension_length ( id_emis, emt_att%nhoursyear, 'nhoursyear' )
                ALLOCATE ( emt_att%hourly_emis_time_factor(emt_att%ncat,emt_att%nhoursyear) )
                CALL get_variable ( id_emis, "emission_time_factors",          &
                                    emt_att%hourly_emis_time_factor,           &
                                    1, emt_att%ncat, 1, emt_att%nhoursyear )

!
!-- MDH

             ELSE IF  ( TRIM(time_fac_type)  ==  "MDH" .OR.                  &
                        TRIM(time_fac_type)  ==  "mdh" )  THEN

                CALL get_dimension_length ( id_emis, emt_att%nmonthdayhour, 'nmonthdayhour' )
                ALLOCATE ( emt_att%mdh_emis_time_factor(emt_att%ncat,emt_att%nmonthdayhour) )
                CALL get_variable ( id_emis, "emission_time_factors",          &
                                    emt_att%mdh_emis_time_factor,              &
                                    1, emt_att%ncat, 1, emt_att%nmonthdayhour )

!
!-- ERROR (time factor undefined)

             ELSE

                message_string = 'We are in the DEFAULT chemistry emissions mode: '  //  &
                                 '     !no time-factor type specified!'              //  &
                                 'Please specify the value of time_fac_type:'        //  &
                                 '         either "MDH" or "HOUR"'                  
                CALL message( 'netcdf_data_input_chemistry_data', 'CM0200', 2, 2, 0, 6, 0 ) 
 

             ENDIF  ! time_fac_type

!
!-- read in default (LOD1) emissions from chemisty netCDF file per species

!
!-- NOTE - at the moment the data is read in per species, but in the future it would
!--        be much more sensible to read in per species per time step to reduce
!--        memory consumption and, to a lesser degree, dimensionality of data exchange
!--        (I expect this will be necessary when the problem size is large)

             DO ispec = 1, emt_att%n_emiss_species

!
!-- allocate space for species specific emission values
!-- NOTE - this array is extended by 1 cell in each horizontal direction
!--        to compensate for an apparent linear offset.  The reason of this
!--        offset is not known but it has been determined to take place beyond the
!--        scope of this module, and has little to do with index conventions.
!--        That is, setting the array horizontal limit from nx0:nx1 to 1:(nx1-nx0+1)
!--        or nx0+1:nx1+1 did not result in correct or definite behavior
!--        This must be looked at at some point by the Hannover team but for now
!--        this workaround is deemed reasonable (ecc 20190417)

                IF ( .NOT. ALLOCATED ( emt(ispec)%default_emission_data ) )  THEN
                    ALLOCATE ( emt(ispec)%default_emission_data(emt_att%ncat,nys:nyn+1,nxl:nxr+1) )
                ENDIF
!
!-- allocate dummy variable w/ index order identical to that shown in the netCDF header

                ALLOCATE ( dum_var_5d(1,nys:nyn,nxl:nxr,1,emt_att%ncat) )
!
!-- get variable.  be very careful
!-- I am using get_variable_5d_real_dynamic (note logical argument at the end)
!-- 1) use Fortran index convention (i.e., 1 to N)
!-- 2) index order must be in reverse order from above allocation order
 
                CALL get_variable ( id_emis, "emission_values", dum_var_5d, &
                                    1,            ispec, nxl+1,     nys+1,     1,                    &
                                    emt_att%ncat, 1,     nxr-nxl+1, nyn-nys+1, emt_att%dt_emission,  &
                                    .FALSE. )
!
!-- assign temp array to data structure then deallocate temp array
!-- NOTE - indices are shifted from nx0:nx1 to nx0+1:nx1+1 to offset
!--        the emission data array to counter said domain offset
!--        (ecc 20190417)

                DO k = 1, emt_att%ncat
                   DO j = nys+1, nyn+1
                      DO i = nxl+1, nxr+1
                         emt(ispec)%default_emission_data(k,j,i) = dum_var_5d(1,j-1,i-1,1,k)
                      ENDDO
                   ENDDO
                ENDDO

                DEALLOCATE ( dum_var_5d )

             ENDDO  ! ispec
!
!-- UNITS

             CALL get_attribute(id_emis,"units",emt_att%units,.FALSE.,"emission_values")

!
!-- END DEFAULT MODE


!
!-- START LOD 2 (PRE-PROCESSED MODE)

          ELSE IF  ( emiss_lod == 2 )  THEN

! for reference (ecc)
!          ELSE IF (TRIM(mode_emis) == "PRE-PROCESSED" .OR. TRIM(mode_emis) == "pre-processed") THEN

!
!-- For LOD 2 only VOC and emission data need be read

!------------------------------------------------------------------------------
!-- NOTE - CHECK ARRAY INDICES FOR READING IN NAMES AND SPECIES
!--        IN LOD2 (PRE-PROCESSED MODE) FOR THE VARIOUS MODE SPLITS
!--        AS ALL ID_EMIS CONDITIONALS HAVE BEEN REMOVED FROM GET_VAR
!--        FUNCTIONS.  IN THEORY THIS WOULD MEAN ALL ARRAYS SHOULD BE
!--        READ FROM 0 to N-1 (C CONVENTION) AS OPPOSED TO 1 to N
!--        (FORTRAN CONVENTION).  KEEP THIS IN MIND !!
!--        (ecc 20190424)
!------------------------------------------------------------------------------

             DO ispec = 1, emt_att%n_emiss_species

!
!-- VOC DATA (name and composition)

                IF  ( TRIM(emt_att%species_name(ispec)) == "VOC" .OR.                  &
                      TRIM(emt_att%species_name(ispec)) == "voc" )  THEN

!
!-- VOC name
                   CALL get_dimension_length ( id_emis, emt_att%nvoc, 'nvoc' )
                   ALLOCATE ( emt_att%voc_name(emt_att%nvoc) )
                   CALL get_variable ( id_emis, "emission_voc_name",                     &
                                       string_values, emt_att%nvoc)
                   emt_att%voc_name = string_values
                   IF  ( ALLOCATED(string_values) )  DEALLOCATE (string_values)

!
!-- VOC composition
 
                   ALLOCATE ( emt_att%voc_comp(emt_att%ncat,emt_att%nvoc) )
                   CALL get_variable ( id_emis, "composition_voc", emt_att%voc_comp,     &
                                       1, emt_att%ncat, 1, emt_att%nvoc )
                ENDIF  ! VOC
  
             ENDDO  ! ispec

!
!-- EMISSION DATA

             CALL get_dimension_length ( id_emis, emt_att%dt_emission, 'time' )   
 
!
!-- read in pre-processed (LOD2) emissions from chemisty netCDF file per species

!
!-- NOTE - at the moment the data is read in per species, but in the future it would
!--        be much more sensible to read in per species per time step to reduce
!--        memory consumption and, to a lesser degree, dimensionality of data exchange
!--        (I expect this will be necessary when the problem size is large)

             DO ispec = 1, emt_att%n_emiss_species

!
!-- allocate space for species specific emission values
!-- NOTE - this array is extended by 1 cell in each horizontal direction
!--        to compensate for an apparent linear offset.  The reason of this
!--        offset is not known but it has been determined to take place beyond the
!--        scope of this module, and has little to do with index conventions.
!--        That is, setting the array horizontal limit from nx0:nx1 to 1:(nx1-nx0+1)
!--        or nx0+1:nx1+1 did not result in correct or definite behavior
!--        This must be looked at at some point by the Hannover team but for now
!--        this workaround is deemed reasonable (ecc 20190417)

                IF ( .NOT. ALLOCATED( emt(ispec)%preproc_emission_data ) )  THEN
                   ALLOCATE( emt(ispec)%preproc_emission_data(                           &
                             emt_att%dt_emission, 1, nys:nyn+1, nxl:nxr+1) )
                ENDIF
!
!-- allocate dummy variable w/ index order identical to that shown in the netCDF header

                ALLOCATE ( dum_var_5d(emt_att%dt_emission,1,nys:nyn,nxl:nxr,1) )
!
!-- get variable.  be very careful
!-- I am using get_variable_5d_real_dynamic (note logical argument at the end)
!-- 1) use Fortran index convention (i.e., 1 to N)
!-- 2) index order must be in reverse order from above allocation order

                CALL get_variable ( id_emis, "emission_values", dum_var_5d, &
                                    ispec, nxl+1,     nys+1,     1, 1,                   &
                                    1,     nxr-nxl+1, nyn-nys+1, 1, emt_att%dt_emission, &
                                    .FALSE. )
!
!-- assign temp array to data structure then deallocate temp array
!-- NOTE - indices are shifted from nx0:nx1 to nx0+1:nx1+1 to offset
!--        the emission data array to counter said unkonwn offset
!--        (ecc 20190417)

                DO k = 1, emt_att%dt_emission
                   DO j = nys+1, nyn+1
                      DO i = nxl+1, nxr+1
                         emt(ispec)%preproc_emission_data(k,1,j,i) = dum_var_5d(k,1,j-1,i-1,1)
                      ENDDO
                   ENDDO
                ENDDO

                DEALLOCATE ( dum_var_5d )

             ENDDO  ! ispec
!
!-- UNITS

             CALL get_attribute ( id_emis, "units", emt_att%units, .FALSE. , "emission_values" )
       
          ENDIF  ! LOD1 & LOD2 (default and pre-processed mode)

          CALL close_input_file (id_emis)

#endif

       ENDIF ! LOD0 (parameterized mode)

    END SUBROUTINE netcdf_data_input_chemistry_data


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads surface classification data, such as vegetation and soil type, etc. .
!------------------------------------------------------------------------------!
    SUBROUTINE netcdf_data_input_surface_data

       USE control_parameters,                                                 &
           ONLY:  land_surface, urban_surface

       USE exchange_horiz_mod,                                                 &
           ONLY:  exchange_horiz_2d, exchange_horiz_2d_byte, exchange_horiz_2d_int

       USE indices,                                                            &
           ONLY:  nbgp, nxl, nxr, nyn, nys


       IMPLICIT NONE

       CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE ::  var_names  !< variable names in static input file

       INTEGER(iwp) ::  id_surf   !< NetCDF id of input file
       INTEGER(iwp) ::  k         !< running index along z-direction
       INTEGER(iwp) ::  k2        !< running index
       INTEGER(iwp) ::  num_vars  !< number of variables in input file
       INTEGER(iwp) ::  nz_soil   !< number of soil layers in file

!
!--    If not static input file is available, skip this routine
       IF ( .NOT. input_pids_static )  RETURN
!
!--    Measure CPU time
       CALL cpu_log( log_point_s(82), 'NetCDF input', 'start' )
!
!--    Skip the following if no land-surface or urban-surface module are
!--    applied. This case, no one of the following variables is used anyway.
       IF (  .NOT. land_surface  .AND.  .NOT. urban_surface )  RETURN

#if defined ( __netcdf )
!
!--    Open file in read-only mode
       CALL open_read_file( TRIM( input_file_static ) //                       &
                            TRIM( coupling_char ) , id_surf )
!
!--    Inquire all variable names.
!--    This will be used to check whether an optional input variable exist
!--    or not.
       CALL inquire_num_variables( id_surf, num_vars )

       ALLOCATE( var_names(1:num_vars) )
       CALL inquire_variable_names( id_surf, var_names )
!
!--    Read vegetation type and required attributes
       IF ( check_existence( var_names, 'vegetation_type' ) )  THEN
          vegetation_type_f%from_file = .TRUE.
          CALL get_attribute( id_surf, char_fill,                              &
                              vegetation_type_f%fill,                          &
                              .FALSE., 'vegetation_type' )

          ALLOCATE ( vegetation_type_f%var(nys:nyn,nxl:nxr)  )

          CALL get_variable( id_surf, 'vegetation_type',                       &
                             vegetation_type_f%var, nxl, nxr, nys, nyn )
       ELSE
          vegetation_type_f%from_file = .FALSE.
       ENDIF

!
!--    Read soil type and required attributes
       IF ( check_existence( var_names, 'soil_type' ) )  THEN
          soil_type_f%from_file = .TRUE.
!
!--       Note, lod is currently not on file; skip for the moment
!           CALL get_attribute( id_surf, char_lod,                       &
!                                      soil_type_f%lod,                  &
!                                      .FALSE., 'soil_type' )
          CALL get_attribute( id_surf, char_fill,                              &
                              soil_type_f%fill,                                &
                              .FALSE., 'soil_type' )

          IF ( soil_type_f%lod == 1 )  THEN

             ALLOCATE ( soil_type_f%var_2d(nys:nyn,nxl:nxr)  )

             CALL get_variable( id_surf, 'soil_type', soil_type_f%var_2d,      &
                                nxl, nxr, nys, nyn )

          ELSEIF ( soil_type_f%lod == 2 )  THEN
!
!--          Obtain number of soil layers from file.
             CALL get_dimension_length( id_surf, nz_soil, 'zsoil' )

             ALLOCATE ( soil_type_f%var_3d(0:nz_soil,nys:nyn,nxl:nxr) )

             CALL get_variable( id_surf, 'soil_type', soil_type_f%var_3d,      &
                                nxl, nxr, nys, nyn, 0, nz_soil )
 
          ENDIF
       ELSE
          soil_type_f%from_file = .FALSE.
       ENDIF

!
!--    Read pavement type and required attributes
       IF ( check_existence( var_names, 'pavement_type' ) )  THEN
          pavement_type_f%from_file = .TRUE.
          CALL get_attribute( id_surf, char_fill,                              &
                              pavement_type_f%fill, .FALSE.,                   &
                              'pavement_type' )

          ALLOCATE ( pavement_type_f%var(nys:nyn,nxl:nxr)  )

          CALL get_variable( id_surf, 'pavement_type', pavement_type_f%var,    &
                             nxl, nxr, nys, nyn )
       ELSE
          pavement_type_f%from_file = .FALSE.
       ENDIF

!
!--    Read water type and required attributes
       IF ( check_existence( var_names, 'water_type' ) )  THEN
          water_type_f%from_file = .TRUE.
          CALL get_attribute( id_surf, char_fill, water_type_f%fill,           &
                              .FALSE., 'water_type' )

          ALLOCATE ( water_type_f%var(nys:nyn,nxl:nxr)  )

          CALL get_variable( id_surf, 'water_type', water_type_f%var,          &
                             nxl, nxr, nys, nyn )

       ELSE
          water_type_f%from_file = .FALSE.
       ENDIF
!
!--    Read relative surface fractions of vegetation, pavement and water.
       IF ( check_existence( var_names, 'surface_fraction' ) )  THEN
          surface_fraction_f%from_file = .TRUE.
          CALL get_attribute( id_surf, char_fill,                              &
                              surface_fraction_f%fill,                         &
                              .FALSE., 'surface_fraction' )
!
!--       Inquire number of surface fractions
          CALL get_dimension_length( id_surf,                                  &
                                     surface_fraction_f%nf,                    &
                                     'nsurface_fraction' )
!
!--       Allocate dimension array and input array for surface fractions
          ALLOCATE( surface_fraction_f%nfracs(0:surface_fraction_f%nf-1) )
          ALLOCATE( surface_fraction_f%frac(0:surface_fraction_f%nf-1,         &
                                            nys:nyn,nxl:nxr) )
!
!--       Get dimension of surface fractions
          CALL get_variable( id_surf, 'nsurface_fraction',                     &
                             surface_fraction_f%nfracs )
!
!--       Read surface fractions
          CALL get_variable( id_surf, 'surface_fraction',                      &
                             surface_fraction_f%frac, nxl, nxr, nys, nyn,      &
                             0, surface_fraction_f%nf-1 )
       ELSE
          surface_fraction_f%from_file = .FALSE.
       ENDIF
!
!--    Read building parameters and related information
       IF ( check_existence( var_names, 'building_pars' ) )  THEN
          building_pars_f%from_file = .TRUE.
          CALL get_attribute( id_surf, char_fill,                              &
                              building_pars_f%fill,                            &
                              .FALSE., 'building_pars' )
!
!--       Inquire number of building parameters
          CALL get_dimension_length( id_surf,                                  &
                                      building_pars_f%np,                      &
                                      'nbuilding_pars' )
!
!--       Allocate dimension array and input array for building parameters
          ALLOCATE( building_pars_f%pars(0:building_pars_f%np-1) )
          ALLOCATE( building_pars_f%pars_xy(0:building_pars_f%np-1,            &
                                            nys:nyn,nxl:nxr) )
!
!--       Get dimension of building parameters
          CALL get_variable( id_surf, 'nbuilding_pars',                        &
                             building_pars_f%pars )
!
!--       Read building_pars
          CALL get_variable( id_surf, 'building_pars',                         &
                             building_pars_f%pars_xy, nxl, nxr, nys, nyn,      &
                             0, building_pars_f%np-1 )
       ELSE
          building_pars_f%from_file = .FALSE.
       ENDIF
!
!--    Read building surface parameters
       IF ( check_existence( var_names, 'building_surface_pars' ) )  THEN
          building_surface_pars_f%from_file = .TRUE.
          CALL get_attribute( id_surf, char_fill,                              &
                              building_surface_pars_f%fill,                    &
                              .FALSE., 'building_surface_pars' )
!
!--       Read building_surface_pars
          CALL get_variable_surf( id_surf, 'building_surface_pars', &
                                  building_surface_pars_f )
       ELSE
          building_surface_pars_f%from_file = .FALSE.
       ENDIF

!
!--    Read albedo type and required attributes
       IF ( check_existence( var_names, 'albedo_type' ) )  THEN
          albedo_type_f%from_file = .TRUE.
          CALL get_attribute( id_surf, char_fill, albedo_type_f%fill,          &
                              .FALSE.,  'albedo_type' )

          ALLOCATE ( albedo_type_f%var(nys:nyn,nxl:nxr)  )
          
          CALL get_variable( id_surf, 'albedo_type', albedo_type_f%var,        &
                             nxl, nxr, nys, nyn )
       ELSE
          albedo_type_f%from_file = .FALSE.
       ENDIF
!
!--    Read albedo parameters and related information
       IF ( check_existence( var_names, 'albedo_pars' ) )  THEN
          albedo_pars_f%from_file = .TRUE.
          CALL get_attribute( id_surf, char_fill, albedo_pars_f%fill,          &
                              .FALSE., 'albedo_pars' )
!
!--       Inquire number of albedo parameters
          CALL get_dimension_length( id_surf,                                  &
                                     albedo_pars_f%np,                         &
                                     'nalbedo_pars' )
!
!--       Allocate dimension array and input array for albedo parameters
          ALLOCATE( albedo_pars_f%pars(0:albedo_pars_f%np-1) )
          ALLOCATE( albedo_pars_f%pars_xy(0:albedo_pars_f%np-1,                &
                                          nys:nyn,nxl:nxr) )
!
!--       Get dimension of albedo parameters
          CALL get_variable( id_surf, 'nalbedo_pars', albedo_pars_f%pars )

          CALL get_variable( id_surf, 'albedo_pars', albedo_pars_f%pars_xy,    &
                             nxl, nxr, nys, nyn,                               &
                             0, albedo_pars_f%np-1 )
       ELSE
          albedo_pars_f%from_file = .FALSE.
       ENDIF

!
!--    Read pavement parameters and related information
       IF ( check_existence( var_names, 'pavement_pars' ) )  THEN
          pavement_pars_f%from_file = .TRUE.
          CALL get_attribute( id_surf, char_fill,                              &
                              pavement_pars_f%fill,                            &
                              .FALSE., 'pavement_pars' )
!
!--       Inquire number of pavement parameters
          CALL get_dimension_length( id_surf,                                  &
                                     pavement_pars_f%np,                       &
                                     'npavement_pars' )
!
!--       Allocate dimension array and input array for pavement parameters
          ALLOCATE( pavement_pars_f%pars(0:pavement_pars_f%np-1) )
          ALLOCATE( pavement_pars_f%pars_xy(0:pavement_pars_f%np-1,            &
                                            nys:nyn,nxl:nxr) )
!
!--       Get dimension of pavement parameters
          CALL get_variable( id_surf, 'npavement_pars', pavement_pars_f%pars )

          CALL get_variable( id_surf, 'pavement_pars', pavement_pars_f%pars_xy,&
                             nxl, nxr, nys, nyn,                               &
                             0, pavement_pars_f%np-1 )
       ELSE
          pavement_pars_f%from_file = .FALSE.
       ENDIF

!
!--    Read pavement subsurface parameters and related information
       IF ( check_existence( var_names, 'pavement_subsurface_pars' ) )         &
       THEN
          pavement_subsurface_pars_f%from_file = .TRUE.
          CALL get_attribute( id_surf, char_fill,                              &
                              pavement_subsurface_pars_f%fill,                 &
                              .FALSE., 'pavement_subsurface_pars' )
!
!--       Inquire number of parameters
          CALL get_dimension_length( id_surf,                                  &
                                     pavement_subsurface_pars_f%np,            &
                                     'npavement_subsurface_pars' )
!
!--       Inquire number of soil layers
          CALL get_dimension_length( id_surf,                                  &
                                     pavement_subsurface_pars_f%nz,            &
                                     'zsoil' )
!
!--       Allocate dimension array and input array for pavement parameters
          ALLOCATE( pavement_subsurface_pars_f%pars                            &
                            (0:pavement_subsurface_pars_f%np-1) )
          ALLOCATE( pavement_subsurface_pars_f%pars_xyz                        &
                            (0:pavement_subsurface_pars_f%np-1,                &
                             0:pavement_subsurface_pars_f%nz-1,                &
                             nys:nyn,nxl:nxr) )
!
!--       Get dimension of pavement parameters
          CALL get_variable( id_surf, 'npavement_subsurface_pars',             &
                             pavement_subsurface_pars_f%pars )

          CALL get_variable( id_surf, 'pavement_subsurface_pars',              &
                             pavement_subsurface_pars_f%pars_xyz,              &
                             nxl, nxr, nys, nyn,                               &
                             0, pavement_subsurface_pars_f%nz-1,               &
                             0, pavement_subsurface_pars_f%np-1 )
       ELSE
          pavement_subsurface_pars_f%from_file = .FALSE.
       ENDIF


!
!--    Read vegetation parameters and related information
       IF ( check_existence( var_names, 'vegetation_pars' ) )  THEN
          vegetation_pars_f%from_file = .TRUE.
          CALL get_attribute( id_surf, char_fill,                              &
                              vegetation_pars_f%fill,                          &
                              .FALSE.,  'vegetation_pars' )
!
!--       Inquire number of vegetation parameters
          CALL get_dimension_length( id_surf,                                  &
                                     vegetation_pars_f%np,                     &
                                     'nvegetation_pars' )
!
!--       Allocate dimension array and input array for surface fractions
          ALLOCATE( vegetation_pars_f%pars(0:vegetation_pars_f%np-1) )
          ALLOCATE( vegetation_pars_f%pars_xy(0:vegetation_pars_f%np-1,        &
                                              nys:nyn,nxl:nxr) )
!
!--       Get dimension of the parameters
          CALL get_variable( id_surf, 'nvegetation_pars',                      &
                             vegetation_pars_f%pars )

          CALL get_variable( id_surf, 'vegetation_pars',                       &
                             vegetation_pars_f%pars_xy, nxl, nxr, nys, nyn,    &
                             0, vegetation_pars_f%np-1 )
       ELSE
          vegetation_pars_f%from_file = .FALSE.
       ENDIF

!
!--    Read root parameters/distribution and related information
       IF ( check_existence( var_names, 'soil_pars' ) )  THEN
          soil_pars_f%from_file = .TRUE.
          CALL get_attribute( id_surf, char_fill,                              &
                              soil_pars_f%fill,                                &
                              .FALSE., 'soil_pars' )

          CALL get_attribute( id_surf, char_lod,                               &
                              soil_pars_f%lod,                                 &
                              .FALSE., 'soil_pars' )

!
!--       Inquire number of soil parameters
          CALL get_dimension_length( id_surf,                                  &
                                     soil_pars_f%np,                           &
                                     'nsoil_pars' )
!
!--       Read parameters array
          ALLOCATE( soil_pars_f%pars(0:soil_pars_f%np-1) )
          CALL get_variable( id_surf, 'nsoil_pars', soil_pars_f%pars )

!
!--       In case of level of detail 2, also inquire number of vertical
!--       soil layers, allocate memory and read the respective dimension
          IF ( soil_pars_f%lod == 2 )  THEN
             CALL get_dimension_length( id_surf,                               &
                                        soil_pars_f%nz,                        &
                                        'zsoil' )

             ALLOCATE( soil_pars_f%layers(0:soil_pars_f%nz-1) )
             CALL get_variable( id_surf, 'zsoil', soil_pars_f%layers )

          ENDIF

!
!--       Read soil parameters, depending on level of detail
          IF ( soil_pars_f%lod == 1 )  THEN
             ALLOCATE( soil_pars_f%pars_xy(0:soil_pars_f%np-1,                 &
                                           nys:nyn,nxl:nxr) )
                  
             CALL get_variable( id_surf, 'soil_pars', soil_pars_f%pars_xy,     &
                                nxl, nxr, nys, nyn, 0, soil_pars_f%np-1 )

          ELSEIF ( soil_pars_f%lod == 2 )  THEN
             ALLOCATE( soil_pars_f%pars_xyz(0:soil_pars_f%np-1,                &
                                            0:soil_pars_f%nz-1,                &
                                            nys:nyn,nxl:nxr) )
             CALL get_variable( id_surf, 'soil_pars',                          &
                                soil_pars_f%pars_xyz,                          &
                                nxl, nxr, nys, nyn, 0, soil_pars_f%nz-1,       &
                                0, soil_pars_f%np-1 )

          ENDIF
       ELSE
          soil_pars_f%from_file = .FALSE.
       ENDIF

!
!--    Read water parameters and related information
       IF ( check_existence( var_names, 'water_pars' ) )  THEN
          water_pars_f%from_file = .TRUE.
          CALL get_attribute( id_surf, char_fill,                              &
                              water_pars_f%fill,                               &
                              .FALSE., 'water_pars' )
!
!--       Inquire number of water parameters
          CALL get_dimension_length( id_surf,                                  &
                                     water_pars_f%np,                          &
                                     'nwater_pars' )
!
!--       Allocate dimension array and input array for water parameters
          ALLOCATE( water_pars_f%pars(0:water_pars_f%np-1) )
          ALLOCATE( water_pars_f%pars_xy(0:water_pars_f%np-1,                  &
                                         nys:nyn,nxl:nxr) )
!
!--       Get dimension of water parameters
          CALL get_variable( id_surf, 'nwater_pars', water_pars_f%pars )

          CALL get_variable( id_surf, 'water_pars', water_pars_f%pars_xy,      &
                             nxl, nxr, nys, nyn, 0, water_pars_f%np-1 )
       ELSE
          water_pars_f%from_file = .FALSE.
       ENDIF
!
!--    Read root area density - parametrized vegetation
       IF ( check_existence( var_names, 'root_area_dens_s' ) )  THEN
          root_area_density_lsm_f%from_file = .TRUE.
          CALL get_attribute( id_surf, char_fill,                              &
                              root_area_density_lsm_f%fill,                    &
                              .FALSE., 'root_area_dens_s' )
!
!--       Obtain number of soil layers from file and allocate variable
          CALL get_dimension_length( id_surf,                                  &
                                     root_area_density_lsm_f%nz,               &
                                     'zsoil' )
          ALLOCATE( root_area_density_lsm_f%var                                &
                                        (0:root_area_density_lsm_f%nz-1,       &
                                         nys:nyn,nxl:nxr) )

!
!--       Read root-area density
          CALL get_variable( id_surf, 'root_area_dens_s',                      &
                             root_area_density_lsm_f%var,                      &
                             nxl, nxr, nys, nyn,                               &
                             0, root_area_density_lsm_f%nz-1 )

       ELSE
          root_area_density_lsm_f%from_file = .FALSE.
       ENDIF
!
!--    Read street type and street crossing
       IF ( check_existence( var_names, 'street_type' ) )  THEN
          street_type_f%from_file = .TRUE.
          CALL get_attribute( id_surf, char_fill,                              &
                              street_type_f%fill, .FALSE.,                     &
                              'street_type' )

          ALLOCATE ( street_type_f%var(nys:nyn,nxl:nxr)  )
          
          CALL get_variable( id_surf, 'street_type', street_type_f%var,        &
                             nxl, nxr, nys, nyn )
       ELSE
          street_type_f%from_file = .FALSE.
       ENDIF

       IF ( check_existence( var_names, 'street_crossing' ) )  THEN
          street_crossing_f%from_file = .TRUE.
          CALL get_attribute( id_surf, char_fill,                              &
                              street_crossing_f%fill, .FALSE.,                 &
                              'street_crossing' )

          ALLOCATE ( street_crossing_f%var(nys:nyn,nxl:nxr)  )

          CALL get_variable( id_surf, 'street_crossing',                       &
                             street_crossing_f%var, nxl, nxr, nys, nyn )

       ELSE
          street_crossing_f%from_file = .FALSE.
       ENDIF
!
!--    Still missing: root_resolved and building_surface_pars.
!--    Will be implemented as soon as they are available.

!
!--    Finally, close input file
       CALL close_input_file( id_surf )
#endif
!
!--    End of CPU measurement
       CALL cpu_log( log_point_s(82), 'NetCDF input', 'stop' )
!
!--    Exchange ghost points for surface variables. Therefore, resize 
!--    variables. 
       IF ( albedo_type_f%from_file )  THEN
          CALL resize_array_2d_int8( albedo_type_f%var, nys, nyn, nxl, nxr )
          CALL exchange_horiz_2d_byte( albedo_type_f%var, nys, nyn, nxl, nxr,  &
                                       nbgp )
       ENDIF
       IF ( pavement_type_f%from_file )  THEN
          CALL resize_array_2d_int8( pavement_type_f%var, nys, nyn, nxl, nxr )
          CALL exchange_horiz_2d_byte( pavement_type_f%var, nys, nyn, nxl, nxr,&
                                       nbgp )
       ENDIF
       IF ( soil_type_f%from_file  .AND.  ALLOCATED( soil_type_f%var_2d ) )  THEN
          CALL resize_array_2d_int8( soil_type_f%var_2d, nys, nyn, nxl, nxr )
          CALL exchange_horiz_2d_byte( soil_type_f%var_2d, nys, nyn, nxl, nxr, &
                                       nbgp )
       ENDIF
       IF ( vegetation_type_f%from_file )  THEN
          CALL resize_array_2d_int8( vegetation_type_f%var, nys, nyn, nxl, nxr )
          CALL exchange_horiz_2d_byte( vegetation_type_f%var, nys, nyn, nxl,   &
                                       nxr, nbgp )
       ENDIF
       IF ( water_type_f%from_file )  THEN
          CALL resize_array_2d_int8( water_type_f%var, nys, nyn, nxl, nxr )
          CALL exchange_horiz_2d_byte( water_type_f%var, nys, nyn, nxl, nxr,   &
                                       nbgp )
       ENDIF
!
!--    Exchange ghost points for 3/4-D variables. For the sake of simplicity,
!--    loop further dimensions to use 2D exchange routines. Unfortunately this
!--    is necessary, else new MPI-data types need to be introduced just for 
!--    2 variables. 
       IF ( soil_type_f%from_file  .AND.  ALLOCATED( soil_type_f%var_3d ) )    &
       THEN
          CALL resize_array_3d_int8( soil_type_f%var_3d, 0, nz_soil,           &
                                     nys, nyn, nxl, nxr )
          DO  k = 0, nz_soil
             CALL exchange_horiz_2d_byte(                                       &
                        soil_type_f%var_3d(k,:,:), nys, nyn, nxl, nxr, nbgp )
          ENDDO
       ENDIF

       IF ( surface_fraction_f%from_file )  THEN
          CALL resize_array_3d_real( surface_fraction_f%frac,                  &
                                     0, surface_fraction_f%nf-1,               &
                                     nys, nyn, nxl, nxr )
          DO  k = 0, surface_fraction_f%nf-1
             CALL exchange_horiz_2d( surface_fraction_f%frac(k,:,:) )
          ENDDO
       ENDIF

       IF ( building_pars_f%from_file )  THEN         
          CALL resize_array_3d_real( building_pars_f%pars_xy,                  &
                                     0, building_pars_f%np-1,                  &
                                     nys, nyn, nxl, nxr )
          DO  k = 0, building_pars_f%np-1
             CALL exchange_horiz_2d( building_pars_f%pars_xy(k,:,:) )
          ENDDO
       ENDIF

       IF ( albedo_pars_f%from_file )  THEN          
          CALL resize_array_3d_real( albedo_pars_f%pars_xy,                    &
                                     0, albedo_pars_f%np-1,                    &
                                     nys, nyn, nxl, nxr )
          DO  k = 0, albedo_pars_f%np-1
             CALL exchange_horiz_2d( albedo_pars_f%pars_xy(k,:,:) )
          ENDDO
       ENDIF

       IF ( pavement_pars_f%from_file )  THEN          
          CALL resize_array_3d_real( pavement_pars_f%pars_xy,                  &
                                     0, pavement_pars_f%np-1,                  &
                                     nys, nyn, nxl, nxr )
          DO  k = 0, pavement_pars_f%np-1
             CALL exchange_horiz_2d( pavement_pars_f%pars_xy(k,:,:) )
          ENDDO
       ENDIF

       IF ( vegetation_pars_f%from_file )  THEN
          CALL resize_array_3d_real( vegetation_pars_f%pars_xy,                &
                                     0, vegetation_pars_f%np-1,                &
                                     nys, nyn, nxl, nxr )
          DO  k = 0, vegetation_pars_f%np-1
             CALL exchange_horiz_2d( vegetation_pars_f%pars_xy(k,:,:) )
          ENDDO
       ENDIF

       IF ( water_pars_f%from_file )  THEN
          CALL resize_array_3d_real( water_pars_f%pars_xy,                     &
                                     0, water_pars_f%np-1,                     &
                                     nys, nyn, nxl, nxr )
          DO  k = 0, water_pars_f%np-1
             CALL exchange_horiz_2d( water_pars_f%pars_xy(k,:,:) )
          ENDDO
       ENDIF

       IF ( root_area_density_lsm_f%from_file )  THEN
          CALL resize_array_3d_real( root_area_density_lsm_f%var,              &
                                     0, root_area_density_lsm_f%nz-1,          &
                                     nys, nyn, nxl, nxr )
          DO  k = 0, root_area_density_lsm_f%nz-1
             CALL exchange_horiz_2d( root_area_density_lsm_f%var(k,:,:) )
          ENDDO
       ENDIF

       IF ( soil_pars_f%from_file )  THEN
          IF ( soil_pars_f%lod == 1 )  THEN
          
             CALL resize_array_3d_real( soil_pars_f%pars_xy,                   &
                                        0, soil_pars_f%np-1,                   &
                                        nys, nyn, nxl, nxr )
             DO  k = 0, soil_pars_f%np-1
                CALL exchange_horiz_2d( soil_pars_f%pars_xy(k,:,:) )
             ENDDO
             
          ELSEIF ( soil_pars_f%lod == 2 )  THEN
             CALL resize_array_4d_real( soil_pars_f%pars_xyz,                  &
                                        0, soil_pars_f%np-1,                   &
                                        0, soil_pars_f%nz-1,                   &
                                        nys, nyn, nxl, nxr )

             DO  k2 = 0, soil_pars_f%nz-1
                DO  k = 0, soil_pars_f%np-1
                   CALL exchange_horiz_2d( soil_pars_f%pars_xyz(k,k2,:,:) )
                ENDDO
             ENDDO
          ENDIF
       ENDIF

       IF ( pavement_subsurface_pars_f%from_file )  THEN          
          CALL resize_array_4d_real( pavement_subsurface_pars_f%pars_xyz,      &
                                     0, pavement_subsurface_pars_f%np-1,       &
                                     0, pavement_subsurface_pars_f%nz-1,       &
                                     nys, nyn, nxl, nxr )

          DO  k2 = 0, pavement_subsurface_pars_f%nz-1
             DO  k = 0, pavement_subsurface_pars_f%np-1
                CALL exchange_horiz_2d(                                        &
                           pavement_subsurface_pars_f%pars_xyz(k,k2,:,:) )
             ENDDO
          ENDDO
       ENDIF

    END SUBROUTINE netcdf_data_input_surface_data

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads uvem lookup table information.
!------------------------------------------------------------------------------!
    SUBROUTINE netcdf_data_input_uvem
       
       USE indices,                                                            &
           ONLY:  nxl, nxr, nyn, nys

       IMPLICIT NONE

       CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE ::  var_names  !< variable names in static input file


       INTEGER(iwp) ::  id_uvem       !< NetCDF id of uvem lookup table input file
       INTEGER(iwp) ::  nli = 35      !< dimension length of lookup table in x
       INTEGER(iwp) ::  nlj =  9      !< dimension length of lookup table in y
       INTEGER(iwp) ::  nlk = 90      !< dimension length of lookup table in z
       INTEGER(iwp) ::  num_vars      !< number of variables in netcdf input file
!
!--    Input via uv exposure model lookup table input
       IF ( input_pids_uvem )  THEN

#if defined ( __netcdf )
!
!--       Open file in read-only mode
          CALL open_read_file( TRIM( input_file_uvem ) //                    &
                               TRIM( coupling_char ), id_uvem )
!
!--       At first, inquire all variable names.
!--       This will be used to check whether an input variable exist or not.
          CALL inquire_num_variables( id_uvem, num_vars )
!
!--       Allocate memory to store variable names and inquire them.
          ALLOCATE( var_names(1:num_vars) )
          CALL inquire_variable_names( id_uvem, var_names )
!
!--       uvem integration
          IF ( check_existence( var_names, 'int_factors' ) )  THEN
             uvem_integration_f%from_file = .TRUE.
!
!--          Input 2D uvem integration.
             ALLOCATE ( uvem_integration_f%var(0:nlj,0:nli)  )
             
             CALL get_variable( id_uvem, 'int_factors', uvem_integration_f%var, 0, nli, 0, nlj )
          ELSE
             uvem_integration_f%from_file = .FALSE.
          ENDIF
!
!--       uvem irradiance
          IF ( check_existence( var_names, 'irradiance' ) )  THEN
             uvem_irradiance_f%from_file = .TRUE.
!
!--          Input 2D uvem irradiance.
             ALLOCATE ( uvem_irradiance_f%var(0:nlk, 0:2)  )
             
             CALL get_variable( id_uvem, 'irradiance', uvem_irradiance_f%var, 0, 2, 0, nlk )
          ELSE
             uvem_irradiance_f%from_file = .FALSE.
          ENDIF
!
!--       uvem porjection areas
          IF ( check_existence( var_names, 'projarea' ) )  THEN
             uvem_projarea_f%from_file = .TRUE.
!
!--          Input 3D uvem projection area (human geometgry)
             ALLOCATE ( uvem_projarea_f%var(0:2,0:nlj,0:nli)  )
           
             CALL get_variable( id_uvem, 'projarea', uvem_projarea_f%var, 0, nli, 0, nlj, 0, 2 )
          ELSE
             uvem_projarea_f%from_file = .FALSE.
          ENDIF
!
!--       uvem radiance
          IF ( check_existence( var_names, 'radiance' ) )  THEN
             uvem_radiance_f%from_file = .TRUE.
!
!--          Input 3D uvem radiance
             ALLOCATE ( uvem_radiance_f%var(0:nlk,0:nlj,0:nli)  )
             
             CALL get_variable( id_uvem, 'radiance', uvem_radiance_f%var, 0, nli, 0, nlj, 0, nlk )
          ELSE
             uvem_radiance_f%from_file = .FALSE.
          ENDIF
!
!--       Read building obstruction
          IF ( check_existence( var_names, 'obstruction' ) )  THEN
             building_obstruction_full%from_file = .TRUE.
!--          Input 3D uvem building obstruction
              ALLOCATE ( building_obstruction_full%var_3d(0:44,0:2,0:2) )
              CALL get_variable( id_uvem, 'obstruction', building_obstruction_full%var_3d,0, 2, 0, 2, 0, 44 )        
          ELSE
             building_obstruction_full%from_file = .FALSE.
          ENDIF
!
          IF ( check_existence( var_names, 'obstruction' ) )  THEN
             building_obstruction_f%from_file = .TRUE.
!
!--          Input 3D uvem building obstruction
             ALLOCATE ( building_obstruction_f%var_3d(0:44,nys:nyn,nxl:nxr) )
!
             CALL get_variable( id_uvem, 'obstruction', building_obstruction_f%var_3d,      &
                                nxl, nxr, nys, nyn, 0, 44 )        
          ELSE
             building_obstruction_f%from_file = .FALSE.
          ENDIF
!
!--       Close uvem lookup table input file
          CALL close_input_file( id_uvem )
#else
          CONTINUE
#endif
       ENDIF
    END SUBROUTINE netcdf_data_input_uvem

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads orography and building information.
!------------------------------------------------------------------------------!
    SUBROUTINE netcdf_data_input_topo

       USE control_parameters,                                                 &
           ONLY:  message_string, topography

       USE exchange_horiz_mod,                                                 &
           ONLY:  exchange_horiz_2d_byte, exchange_horiz_2d_int

       USE grid_variables,                                                     &
           ONLY:  dx, dy    
           
       USE indices,                                                            &
           ONLY:  nbgp, nx, nxl, nxr, ny, nyn, nys, nzb


       IMPLICIT NONE

       CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE ::  var_names  !< variable names in static input file


       INTEGER(iwp) ::  i             !< running index along x-direction
       INTEGER(iwp) ::  ii            !< running index for IO blocks
       INTEGER(iwp) ::  id_topo       !< NetCDF id of topograhy input file
       INTEGER(iwp) ::  io_status     !< status after reading the ascii topo file
       INTEGER(iwp) ::  j             !< running index along y-direction
       INTEGER(iwp) ::  num_vars      !< number of variables in netcdf input file
       INTEGER(iwp) ::  skip_n_rows   !< counting variable to skip rows while reading topography file

       REAL(wp) ::  dum           !< dummy variable to skip columns while reading topography file
!
!--    CPU measurement
       CALL cpu_log( log_point_s(83), 'NetCDF/ASCII input topo', 'start' )

!
!--    Input via palm-input data standard
       IF ( input_pids_static )  THEN
#if defined ( __netcdf )
!
!--       Open file in read-only mode
          CALL open_read_file( TRIM( input_file_static ) //                    &
                               TRIM( coupling_char ), id_topo )
!
!--       At first, inquire all variable names.
!--       This will be used to check whether an  input variable exist
!--       or not.
          CALL inquire_num_variables( id_topo, num_vars )
!
!--       Allocate memory to store variable names and inquire them.
          ALLOCATE( var_names(1:num_vars) )
          CALL inquire_variable_names( id_topo, var_names )
!
!--       Read x, y - dimensions. Only required for consistency checks.
          CALL get_dimension_length( id_topo, dim_static%nx, 'x' )
          CALL get_dimension_length( id_topo, dim_static%ny, 'y' )
          ALLOCATE( dim_static%x(0:dim_static%nx-1) )
          ALLOCATE( dim_static%y(0:dim_static%ny-1) )
          CALL get_variable( id_topo, 'x', dim_static%x )
          CALL get_variable( id_topo, 'y', dim_static%y )
!
!--       Check whether dimension size in input file matches the model dimensions
          IF ( dim_static%nx-1 /= nx  .OR.  dim_static%ny-1 /= ny )  THEN
             message_string = 'Static input file: horizontal dimension in ' // &
                              'x- and/or y-direction ' //                      &
                              'do not match the respective model dimension'
             CALL message( 'netcdf_data_input_mod', 'PA0548', 1, 2, 0, 6, 0 )
          ENDIF
!
!--       Check if grid spacing of provided input data matches the respective
!--       grid spacing in the model. 
          IF ( ABS( dim_static%x(1) - dim_static%x(0) - dx ) > 10E-6_wp  .OR.  &
               ABS( dim_static%y(1) - dim_static%y(0) - dy ) > 10E-6_wp )  THEN
             message_string = 'Static input file: horizontal grid spacing ' // &
                              'in x- and/or y-direction ' //                   &
                              'do not match the respective model grid spacing.'
             CALL message( 'netcdf_data_input_mod', 'PA0549', 1, 2, 0, 6, 0 )
          ENDIF
!
!--       Terrain height. First, get variable-related _FillValue attribute
          IF ( check_existence( var_names, 'zt' ) )  THEN
             terrain_height_f%from_file = .TRUE.
             CALL get_attribute( id_topo, char_fill, terrain_height_f%fill,    &
                                 .FALSE., 'zt' )
!
!--          Input 2D terrain height.
             ALLOCATE ( terrain_height_f%var(nys:nyn,nxl:nxr)  )
             
             CALL get_variable( id_topo, 'zt', terrain_height_f%var,           &
                                nxl, nxr, nys, nyn )

          ELSE
             terrain_height_f%from_file = .FALSE.
          ENDIF

!
!--       Read building height. First, read its _FillValue attribute,
!--       as well as lod attribute
          buildings_f%from_file = .FALSE.
          IF ( check_existence( var_names, 'buildings_2d' ) )  THEN
             buildings_f%from_file = .TRUE.
             CALL get_attribute( id_topo, char_lod, buildings_f%lod,           &
                                 .FALSE., 'buildings_2d' )

             CALL get_attribute( id_topo, char_fill, buildings_f%fill1,        &
                                 .FALSE., 'buildings_2d' )

!
!--          Read 2D buildings
             IF ( buildings_f%lod == 1 )  THEN
                ALLOCATE ( buildings_f%var_2d(nys:nyn,nxl:nxr) )

                CALL get_variable( id_topo, 'buildings_2d',                    &
                                   buildings_f%var_2d,                         &
                                   nxl, nxr, nys, nyn )
             ELSE
                message_string = 'NetCDF attribute lod ' //                    &
                                 '(level of detail) is not set ' //            &
                                 'properly for buildings_2d.'
                CALL message( 'netcdf_data_input_mod', 'PA0540',               &
                               1, 2, 0, 6, 0 )
             ENDIF
          ENDIF
!
!--       If available, also read 3D building information. If both are
!--       available, use 3D information.
          IF ( check_existence( var_names, 'buildings_3d' ) )  THEN
             buildings_f%from_file = .TRUE.
             CALL get_attribute( id_topo, char_lod, buildings_f%lod,           &
                                 .FALSE., 'buildings_3d' )      

             CALL get_attribute( id_topo, char_fill, buildings_f%fill2,        &
                                 .FALSE., 'buildings_3d' )

             CALL get_dimension_length( id_topo, buildings_f%nz, 'z' )
!
!--          Read 3D buildings
             IF ( buildings_f%lod == 2 )  THEN
                ALLOCATE( buildings_f%z(nzb:buildings_f%nz-1) )
                CALL get_variable( id_topo, 'z', buildings_f%z )

                ALLOCATE( buildings_f%var_3d(nzb:buildings_f%nz-1,             &
                                             nys:nyn,nxl:nxr) )
                buildings_f%var_3d = 0
                
                CALL get_variable( id_topo, 'buildings_3d',                    &
                                   buildings_f%var_3d,                         &
                                   nxl, nxr, nys, nyn, 0, buildings_f%nz-1 )
             ELSE
                message_string = 'NetCDF attribute lod ' //                    &
                                 '(level of detail) is not set ' //            &
                                 'properly for buildings_3d.'
                CALL message( 'netcdf_data_input_mod', 'PA0541',               &
                               1, 2, 0, 6, 0 )
             ENDIF
          ENDIF
!
!--       Read building IDs and its FillValue attribute. Further required
!--       for mapping buildings on top of orography.
          IF ( check_existence( var_names, 'building_id' ) )  THEN
             building_id_f%from_file = .TRUE.
             CALL get_attribute( id_topo, char_fill,                           &
                                 building_id_f%fill, .FALSE.,                  &
                                 'building_id' )

             ALLOCATE ( building_id_f%var(nys:nyn,nxl:nxr) )
             
             CALL get_variable( id_topo, 'building_id', building_id_f%var,     &
                                nxl, nxr, nys, nyn )
          ELSE
             building_id_f%from_file = .FALSE.
          ENDIF
!
!--       Read building_type and required attributes.
          IF ( check_existence( var_names, 'building_type' ) )  THEN
             building_type_f%from_file = .TRUE.
             CALL get_attribute( id_topo, char_fill,                           &
                                 building_type_f%fill, .FALSE.,                &
                                 'building_type' )

             ALLOCATE ( building_type_f%var(nys:nyn,nxl:nxr) )

             CALL get_variable( id_topo, 'building_type', building_type_f%var, &
                                nxl, nxr, nys, nyn )

          ELSE
             building_type_f%from_file = .FALSE.
          ENDIF
!
!--       Close topography input file
          CALL close_input_file( id_topo )
#else
          CONTINUE
#endif
!
!--    ASCII input
       ELSEIF ( TRIM( topography ) == 'read_from_file' )  THEN
             
          DO  ii = 0, io_blocks-1
             IF ( ii == io_group )  THEN

                OPEN( 90, FILE='TOPOGRAPHY_DATA'//TRIM( coupling_char ),       &
                      STATUS='OLD', FORM='FORMATTED', IOSTAT=io_status )

                IF ( io_status > 0 )  THEN
                   message_string = 'file TOPOGRAPHY_DATA'//                      &
                                    TRIM( coupling_char )// ' does not exist'
                   CALL message( 'netcdf_data_input_mod', 'PA0208', 1, 2, 0, 6, 0 )
                ENDIF

!
!--             Read topography PE-wise. Rows are read from nyn to nys, columns
!--             are read from nxl to nxr. At first, ny-nyn rows need to be skipped.
                skip_n_rows = 0
                DO WHILE ( skip_n_rows < ny - nyn )
                   READ( 90, * )
                   skip_n_rows = skip_n_rows + 1
                ENDDO
!
!--             Read data from nyn to nys and nxl to nxr. Therefore, skip
!--             column until nxl-1 is reached
                ALLOCATE ( buildings_f%var_2d(nys:nyn,nxl:nxr) )
                DO  j = nyn, nys, -1

                   READ( 90, *, IOSTAT=io_status )                               &
                                   ( dum, i = 0, nxl-1 ),                      &
                                   ( buildings_f%var_2d(j,i), i = nxl, nxr )

                   IF ( io_status > 0 )  THEN
                      WRITE( message_string, '(A,1X,I5,1X,A)' ) 'error reading line', ny-j+1,      &
                                                   'of file TOPOGRAPHY_DATA'//TRIM( coupling_char )
                      CALL message( 'netcdf_data_input_mod', 'PA0209', 2, 2, myid, 6, 0 )
                   ELSEIF ( io_status < 0 )  THEN
                      WRITE( message_string, '(A,1X,I5)' ) 'end of line or file detected for '// &
                                  'file TOPOGRAPHY_DATA'//TRIM( coupling_char )//' at line', ny-j+1
                      CALL message( 'netcdf_data_input_mod', 'PA0704', 2, 2, myid, 6, 0 )
                   ENDIF

                ENDDO

                CLOSE( 90 )
                buildings_f%from_file = .TRUE.

             ENDIF
#if defined( __parallel )
             CALL MPI_BARRIER( comm2d, ierr )
#endif
          ENDDO

       ENDIF
!
!--    End of CPU measurement
       CALL cpu_log( log_point_s(83), 'NetCDF/ASCII input topo', 'stop' )
!
!--    Check for minimum requirement to setup building topography. If buildings
!--    are provided, also an ID and a type are required.
!--    Note, doing this check in check_parameters
!--    will be too late (data will be used for grid inititialization before).
       IF ( input_pids_static )  THEN
          IF ( buildings_f%from_file  .AND.                                    &
               .NOT. building_id_f%from_file )  THEN
             message_string = 'If building heights are prescribed in ' //      &
                              'static input file, also an ID is required.'
             CALL message( 'netcdf_data_input_mod', 'PA0542', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF
!
!--    In case no terrain height is provided by static input file, allocate
!--    array nevertheless and set terrain height to 0, which simplifies
!--    topography initialization.
       IF ( .NOT. terrain_height_f%from_file )  THEN
          ALLOCATE ( terrain_height_f%var(nys:nyn,nxl:nxr) )
          terrain_height_f%var = 0.0_wp
       ENDIF
!
!--    Finally, exchange 1 ghost point for building ID and type.
!--    In case of non-cyclic boundary conditions set Neumann conditions at the
!--    lateral boundaries.
       IF ( building_id_f%from_file )  THEN
          CALL resize_array_2d_int32( building_id_f%var, nys, nyn, nxl, nxr )
          CALL exchange_horiz_2d_int( building_id_f%var, nys, nyn, nxl, nxr,   &
                                      nbgp )
       ENDIF

       IF ( building_type_f%from_file )  THEN
          CALL resize_array_2d_int8( building_type_f%var, nys, nyn, nxl, nxr )
          CALL exchange_horiz_2d_byte( building_type_f%var, nys, nyn, nxl, nxr,   &
                                       nbgp )
       ENDIF

    END SUBROUTINE netcdf_data_input_topo

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads initialization data of u, v, w, pt, q, geostrophic wind components,
!> as well as soil moisture and soil temperature, derived from larger-scale
!> model (COSMO) by Inifor.
!------------------------------------------------------------------------------!
    SUBROUTINE netcdf_data_input_init_3d

       USE arrays_3d,                                                          &
           ONLY:  q, pt, u, v, w, zu, zw

       USE control_parameters,                                                 &
           ONLY:  air_chemistry, bc_lr_cyc, bc_ns_cyc, humidity,               &
                  message_string, neutral

       USE indices,                                                            &
           ONLY:  nx, nxl, nxlu, nxr, ny, nyn, nys, nysv, nzb, nz, nzt

       IMPLICIT NONE

       CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE ::  var_names

       LOGICAL      ::  dynamic_3d = .TRUE. !< flag indicating that 3D data is read from dynamic file
       
       INTEGER(iwp) ::  id_dynamic !< NetCDF id of dynamic input file
       INTEGER(iwp) ::  n          !< running index for chemistry variables
       INTEGER(iwp) ::  num_vars   !< number of variables in netcdf input file

       LOGICAL      ::  check_passed !< flag indicating if a check passed

!
!--    Skip routine if no input file with dynamic input data is available.
       IF ( .NOT. input_pids_dynamic )  RETURN
!
!--    Please note, Inifor is designed to provide initial data for u and v for
!--    the prognostic grid points in case of lateral Dirichlet conditions.
!--    This means that Inifor provides data from nxlu:nxr (for u) and
!--    from nysv:nyn (for v) at the left and south domain boundary, respectively.
!--    However, as work-around for the moment, PALM will run with cyclic
!--    conditions and will be initialized with data provided by Inifor
!--    boundaries in case of Dirichlet.
!--    Hence, simply set set nxlu/nysv to 1 (will be reset to its original value
!--    at the end of this routine.
       IF ( bc_lr_cyc  .AND.  nxl == 0 )  nxlu = 1
       IF ( bc_ns_cyc  .AND.  nys == 0 )  nysv = 1

!
!--    CPU measurement
       CALL cpu_log( log_point_s(85), 'NetCDF input init', 'start' )

#if defined ( __netcdf )
!
!--    Open file in read-only mode
       CALL open_read_file( TRIM( input_file_dynamic ) //                      &
                            TRIM( coupling_char ), id_dynamic )

!
!--    At first, inquire all variable names.
       CALL inquire_num_variables( id_dynamic, num_vars )
!
!--    Allocate memory to store variable names.
       ALLOCATE( var_names(1:num_vars) )
       CALL inquire_variable_names( id_dynamic, var_names )
!
!--    Read vertical dimension of scalar und w grid. 
       CALL get_dimension_length( id_dynamic, init_3d%nzu, 'z'     )
       CALL get_dimension_length( id_dynamic, init_3d%nzw, 'zw'    )
!
!--    Read also the horizontal dimensions. These are used just used fo
!--    checking the compatibility with the PALM grid before reading.
       CALL get_dimension_length( id_dynamic, init_3d%nx,  'x'  )
       CALL get_dimension_length( id_dynamic, init_3d%nxu, 'xu' )
       CALL get_dimension_length( id_dynamic, init_3d%ny,  'y'  )
       CALL get_dimension_length( id_dynamic, init_3d%nyv, 'yv' )

!
!--    Check for correct horizontal and vertical dimension. Please note,
!--    checks are performed directly here and not called from
!--    check_parameters as some varialbes are still not allocated there.
!--    Moreover, please note, u- and v-grid has 1 grid point less on
!--    Inifor grid.
       IF ( init_3d%nx-1 /= nx  .OR.  init_3d%nxu-1 /= nx - 1  .OR.            &
            init_3d%ny-1 /= ny  .OR.  init_3d%nyv-1 /= ny - 1 )  THEN
          message_string = 'Number of horizontal grid points in '//            &
                           'dynamic input file does not match ' //             &
                           'the number of numeric grid points.'
          CALL message( 'netcdf_data_input_mod', 'PA0543', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( init_3d%nzu /= nz )  THEN
          message_string = 'Number of vertical grid points in '//              &
                           'dynamic input file does not match ' //             &
                           'the number of numeric grid points.'
          CALL message( 'netcdf_data_input_mod', 'PA0543', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Read vertical dimensions. Later, these are required for eventual
!--    inter- and extrapolations of the initialization data.
       IF ( check_existence( var_names, 'z' ) )  THEN
          ALLOCATE( init_3d%zu_atmos(1:init_3d%nzu) )
          CALL get_variable( id_dynamic, 'z', init_3d%zu_atmos )
       ENDIF
       IF ( check_existence( var_names, 'zw' ) )  THEN
          ALLOCATE( init_3d%zw_atmos(1:init_3d%nzw) )
          CALL get_variable( id_dynamic, 'zw', init_3d%zw_atmos )
       ENDIF
!
!--    Check for consistency between vertical coordinates in dynamic 
!--    driver and numeric grid. 
!--    Please note, depending on compiler options both may be 
!--    equal up to a certain threshold, and differences between
!--    the numeric grid and vertical coordinate in the driver can built-
!--    up to 10E-1-10E-0 m. For this reason, the check is performed not
!--    for exactly matching values. 
       IF ( ANY( ABS( zu(1:nzt)   - init_3d%zu_atmos(1:init_3d%nzu) )    &
                      > 10E-1 )  .OR.                                    &
            ANY( ABS( zw(1:nzt-1) - init_3d%zw_atmos(1:init_3d%nzw) )    &
                      > 10E-1 ) )  THEN
          message_string = 'Vertical grid in dynamic driver does not '// &
                           'match the numeric grid.'
          CALL message( 'netcdf_data_input_mod', 'PA0543', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Read initial geostrophic wind components at 
!--    t = 0 (index 1 in file). 
       IF ( check_existence( var_names, 'ls_forcing_ug' ) )  THEN
          ALLOCATE( init_3d%ug_init(nzb:nzt+1) )
          init_3d%ug_init = 0.0_wp

          CALL get_variable_pr( id_dynamic, 'ls_forcing_ug', 1,          &
                                init_3d%ug_init(1:nzt) )
!
!--       Set top-boundary condition (Neumann)
          init_3d%ug_init(nzt+1) = init_3d%ug_init(nzt)

          init_3d%from_file_ug = .TRUE.
       ELSE
          init_3d%from_file_ug = .FALSE.
       ENDIF
       IF ( check_existence( var_names, 'ls_forcing_vg' ) )  THEN
          ALLOCATE( init_3d%vg_init(nzb:nzt+1) )
          init_3d%vg_init = 0.0_wp

          CALL get_variable_pr( id_dynamic, 'ls_forcing_vg', 1,          &
                                init_3d%vg_init(1:nzt) )
!
!--       Set top-boundary condition (Neumann)
          init_3d%vg_init(nzt+1) = init_3d%vg_init(nzt)

          init_3d%from_file_vg = .TRUE.
       ELSE
          init_3d%from_file_vg = .FALSE.
       ENDIF
!
!--    Read inital 3D data of u, v, w, pt and q,
!--    derived from COSMO model. Read PE-wise yz-slices.
!--    Please note, the u-, v- and w-component are defined on different
!--    grids with one element less in the x-, y-,
!--    and z-direction, respectively. Hence, reading is subdivided
!--    into separate loops.  
!--    Read u-component
       IF ( check_existence( var_names, 'init_atmosphere_u' ) )  THEN
!
!--       Read attributes for the fill value and level-of-detail
          CALL get_attribute( id_dynamic, char_fill, init_3d%fill_u,           &
                              .FALSE., 'init_atmosphere_u' )
          CALL get_attribute( id_dynamic, char_lod, init_3d%lod_u,             &
                              .FALSE., 'init_atmosphere_u' )
!
!--       level-of-detail 1 - read initialization profile
          IF ( init_3d%lod_u == 1 )  THEN
             ALLOCATE( init_3d%u_init(nzb:nzt+1) )
             init_3d%u_init = 0.0_wp

             CALL get_variable( id_dynamic, 'init_atmosphere_u',               &
                                init_3d%u_init(nzb+1:nzt) )
!
!--          Set top-boundary condition (Neumann)
             init_3d%u_init(nzt+1) = init_3d%u_init(nzt)
!
!--       level-of-detail 2 - read 3D initialization data
          ELSEIF ( init_3d%lod_u == 2 )  THEN
             CALL get_variable( id_dynamic, 'init_atmosphere_u',               &
                                u(nzb+1:nzt,nys:nyn,nxlu:nxr),                 &
                                nxlu, nys+1, nzb+1,                            &
                                nxr-nxlu+1, nyn-nys+1, init_3d%nzu,            &
                                dynamic_3d )
!
!--          Set value at leftmost model grid point nxl = 0. This is because
!--          Inifor provides data only from 1:nx-1 since it assumes non-cyclic
!--          conditions.
             IF ( nxl == 0 )                                                   &
                u(nzb+1:nzt,nys:nyn,nxl) = u(nzb+1:nzt,nys:nyn,nxlu)
!
!--          Set bottom and top-boundary
             u(nzb,:,:)   = u(nzb+1,:,:)
             u(nzt+1,:,:) = u(nzt,:,:)
             
          ENDIF
          init_3d%from_file_u = .TRUE.
       ELSE
          message_string = 'Missing initial data for u-component'
          CALL message( 'netcdf_data_input_mod', 'PA0544', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Read v-component
       IF ( check_existence( var_names, 'init_atmosphere_v' ) )  THEN
!
!--       Read attributes for the fill value and level-of-detail
          CALL get_attribute( id_dynamic, char_fill, init_3d%fill_v,           &
                              .FALSE., 'init_atmosphere_v' )
          CALL get_attribute( id_dynamic, char_lod, init_3d%lod_v,             &
                              .FALSE., 'init_atmosphere_v' )
!
!--       level-of-detail 1 - read initialization profile
          IF ( init_3d%lod_v == 1 )  THEN
             ALLOCATE( init_3d%v_init(nzb:nzt+1) )
             init_3d%v_init = 0.0_wp

             CALL get_variable( id_dynamic, 'init_atmosphere_v',               &
                                init_3d%v_init(nzb+1:nzt) )
!
!--          Set top-boundary condition (Neumann)
             init_3d%v_init(nzt+1) = init_3d%v_init(nzt)
!
!--       level-of-detail 2 - read 3D initialization data
          ELSEIF ( init_3d%lod_v == 2 )  THEN
          
             CALL get_variable( id_dynamic, 'init_atmosphere_v',               &
                                v(nzb+1:nzt,nysv:nyn,nxl:nxr),                 &
                                nxl+1, nysv, nzb+1,                            &
                                nxr-nxl+1, nyn-nysv+1, init_3d%nzu,            &
                                dynamic_3d )
!
!--          Set value at southmost model grid point nys = 0. This is because
!--          Inifor provides data only from 1:ny-1 since it assumes non-cyclic
!--          conditions.
             IF ( nys == 0 )                                                   &
                v(nzb+1:nzt,nys,nxl:nxr) = v(nzb+1:nzt,nysv,nxl:nxr)                                
!
!--          Set bottom and top-boundary
             v(nzb,:,:)   = v(nzb+1,:,:)
             v(nzt+1,:,:) = v(nzt,:,:)
             
          ENDIF
          init_3d%from_file_v = .TRUE.
       ELSE
          message_string = 'Missing initial data for v-component'
          CALL message( 'netcdf_data_input_mod', 'PA0544', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Read w-component
       IF ( check_existence( var_names, 'init_atmosphere_w' ) )  THEN
!
!--       Read attributes for the fill value and level-of-detail
          CALL get_attribute( id_dynamic, char_fill, init_3d%fill_w,           &
                              .FALSE., 'init_atmosphere_w' )
          CALL get_attribute( id_dynamic, char_lod, init_3d%lod_w,             &
                              .FALSE., 'init_atmosphere_w' )
!
!--       level-of-detail 1 - read initialization profile
          IF ( init_3d%lod_w == 1 )  THEN
             ALLOCATE( init_3d%w_init(nzb:nzt+1) )
             init_3d%w_init = 0.0_wp

             CALL get_variable( id_dynamic, 'init_atmosphere_w',               &
                                init_3d%w_init(nzb+1:nzt-1) )
!
!--          Set top-boundary condition (Neumann)
             init_3d%w_init(nzt:nzt+1) = init_3d%w_init(nzt-1)
!
!--       level-of-detail 2 - read 3D initialization data
          ELSEIF ( init_3d%lod_w == 2 )  THEN

             CALL get_variable( id_dynamic, 'init_atmosphere_w',                &
                                w(nzb+1:nzt-1,nys:nyn,nxl:nxr),                 &
                                nxl+1, nys+1, nzb+1,                            &
                                nxr-nxl+1, nyn-nys+1, init_3d%nzw,              &
                                dynamic_3d )
!
!--          Set bottom and top-boundary                                
             w(nzb,:,:)   = 0.0_wp  
             w(nzt,:,:)   = w(nzt-1,:,:)
             w(nzt+1,:,:) = w(nzt-1,:,:)

          ENDIF
          init_3d%from_file_w = .TRUE.
       ELSE
          message_string = 'Missing initial data for w-component'
          CALL message( 'netcdf_data_input_mod', 'PA0544', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Read potential temperature
       IF ( .NOT. neutral )  THEN
          IF ( check_existence( var_names, 'init_atmosphere_pt' ) )  THEN
!
!--          Read attributes for the fill value and level-of-detail
             CALL get_attribute( id_dynamic, char_fill, init_3d%fill_pt,       &
                                 .FALSE., 'init_atmosphere_pt' )
             CALL get_attribute( id_dynamic, char_lod, init_3d%lod_pt,         &
                                 .FALSE., 'init_atmosphere_pt' )
!
!--          level-of-detail 1 - read initialization profile
             IF ( init_3d%lod_pt == 1 )  THEN
                ALLOCATE( init_3d%pt_init(nzb:nzt+1) )

                CALL get_variable( id_dynamic, 'init_atmosphere_pt',           &
                                   init_3d%pt_init(nzb+1:nzt) )
!
!--             Set Neumann top and surface boundary condition for initial 
!--             profil
                init_3d%pt_init(nzb)   = init_3d%pt_init(nzb+1)
                init_3d%pt_init(nzt+1) = init_3d%pt_init(nzt)
!
!--          level-of-detail 2 - read 3D initialization data
             ELSEIF ( init_3d%lod_pt == 2 )  THEN

                CALL get_variable( id_dynamic, 'init_atmosphere_pt',           &
                                   pt(nzb+1:nzt,nys:nyn,nxl:nxr),              &
                                   nxl+1, nys+1, nzb+1,                        &
                                   nxr-nxl+1, nyn-nys+1, init_3d%nzu,          &
                                   dynamic_3d )
                                   
!
!--             Set bottom and top-boundary
                pt(nzb,:,:)   = pt(nzb+1,:,:)
                pt(nzt+1,:,:) = pt(nzt,:,:)             

             ENDIF
             init_3d%from_file_pt = .TRUE.
          ELSE
             message_string = 'Missing initial data for ' //                   &
                              'potential temperature'
             CALL message( 'netcdf_data_input_mod', 'PA0544', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF
!
!--    Read mixing ratio
       IF ( humidity )  THEN
          IF ( check_existence( var_names, 'init_atmosphere_qv' ) )  THEN
!
!--          Read attributes for the fill value and level-of-detail
             CALL get_attribute( id_dynamic, char_fill, init_3d%fill_q,        &
                                 .FALSE., 'init_atmosphere_qv' )
             CALL get_attribute( id_dynamic, char_lod, init_3d%lod_q,          &
                                 .FALSE., 'init_atmosphere_qv' )
!
!--          level-of-detail 1 - read initialization profile
             IF ( init_3d%lod_q == 1 )  THEN
                ALLOCATE( init_3d%q_init(nzb:nzt+1) )

                CALL get_variable( id_dynamic, 'init_atmosphere_qv',           &
                                    init_3d%q_init(nzb+1:nzt) )
!
!--             Set bottom and top boundary condition (Neumann)
                init_3d%q_init(nzb)   = init_3d%q_init(nzb+1)
                init_3d%q_init(nzt+1) = init_3d%q_init(nzt)
!
!--          level-of-detail 2 - read 3D initialization data
             ELSEIF ( init_3d%lod_q == 2 )  THEN
             
                CALL get_variable( id_dynamic, 'init_atmosphere_qv',           &
                                   q(nzb+1:nzt,nys:nyn,nxl:nxr),               &
                                   nxl+1, nys+1, nzb+1,                        &
                                   nxr-nxl+1, nyn-nys+1, init_3d%nzu,          &
                                   dynamic_3d )
                                   
!
!--             Set bottom and top-boundary
                q(nzb,:,:)   = q(nzb+1,:,:)
                q(nzt+1,:,:) = q(nzt,:,:)
                
             ENDIF
             init_3d%from_file_q = .TRUE.
          ELSE
             message_string = 'Missing initial data for ' //                   &
                              'mixing ratio'
             CALL message( 'netcdf_data_input_mod', 'PA0544', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF       
!
!--    Read chemistry variables.
!--    Please note, for the moment, only LOD=1 is allowed
       IF ( air_chemistry )  THEN
!
!--       Allocate chemistry input profiles, as well as arrays for fill values
!--       and LOD's.
          ALLOCATE( init_3d%chem_init(nzb:nzt+1,                               &
                                      1:UBOUND(init_3d%var_names_chem, 1 )) )
          ALLOCATE( init_3d%fill_chem(1:UBOUND(init_3d%var_names_chem, 1)) )    
          ALLOCATE( init_3d%lod_chem(1:UBOUND(init_3d%var_names_chem, 1))  )  
          
          DO  n = 1, UBOUND(init_3d%var_names_chem, 1)
             IF ( check_existence( var_names,                                  &
                                   TRIM( init_3d%var_names_chem(n) ) ) )  THEN
!
!--             Read attributes for the fill value and level-of-detail
                CALL get_attribute( id_dynamic, char_fill,                     &
                                    init_3d%fill_chem(n),                      &
                                    .FALSE.,                                   &
                                    TRIM( init_3d%var_names_chem(n) ) )
                CALL get_attribute( id_dynamic, char_lod,                      &
                                    init_3d%lod_chem(n),                       &
                                    .FALSE.,                                   &
                                    TRIM( init_3d%var_names_chem(n) ) )
!
!--             Give message that only LOD=1 is allowed. 
                IF ( init_3d%lod_chem(n) /= 1 )  THEN                
                   message_string = 'For chemistry variables only LOD=1 is ' //&
                                    'allowed.'
                   CALL message( 'netcdf_data_input_mod', 'PA0586',            &
                                 1, 2, 0, 6, 0 )
                ENDIF
!
!--             level-of-detail 1 - read initialization profile
                CALL get_variable( id_dynamic,                                 &
                                   TRIM( init_3d%var_names_chem(n) ),          &
                                   init_3d%chem_init(nzb+1:nzt,n) )
!
!--             Set bottom and top boundary condition (Neumann)
                init_3d%chem_init(nzb,n)   = init_3d%chem_init(nzb+1,n)
                init_3d%chem_init(nzt+1,n) = init_3d%chem_init(nzt,n)
                
                init_3d%from_file_chem(n) = .TRUE.
             ENDIF
          ENDDO
       ENDIF
!
!--    Close input file
       CALL close_input_file( id_dynamic )
#endif
!
!--    End of CPU measurement
       CALL cpu_log( log_point_s(85), 'NetCDF input init', 'stop' )
!
!--    Finally, check if the input data has any fill values. Please note,
!--    checks depend on the LOD of the input data.
       IF ( init_3d%from_file_u )  THEN
          check_passed = .TRUE.
          IF ( init_3d%lod_u == 1 )  THEN
             IF ( ANY( init_3d%u_init(nzb+1:nzt+1) == init_3d%fill_u ) )       &
                check_passed = .FALSE.
          ELSEIF ( init_3d%lod_u == 2 )  THEN
             IF ( ANY( u(nzb+1:nzt+1,nys:nyn,nxlu:nxr) == init_3d%fill_u ) )   &
                check_passed = .FALSE.
          ENDIF
          IF ( .NOT. check_passed )  THEN
             message_string = 'NetCDF input for init_atmosphere_u must ' //    &
                              'not contain any _FillValues'
             CALL message( 'netcdf_data_input_mod', 'PA0545', 2, 2, 0, 6, 0 )
          ENDIF
       ENDIF

       IF ( init_3d%from_file_v )  THEN
          check_passed = .TRUE.
          IF ( init_3d%lod_v == 1 )  THEN
             IF ( ANY( init_3d%v_init(nzb+1:nzt+1) == init_3d%fill_v ) )       &
                check_passed = .FALSE.
          ELSEIF ( init_3d%lod_v == 2 )  THEN
             IF ( ANY( v(nzb+1:nzt+1,nysv:nyn,nxl:nxr) == init_3d%fill_v ) )   &
                check_passed = .FALSE.
          ENDIF
          IF ( .NOT. check_passed )  THEN
             message_string = 'NetCDF input for init_atmosphere_v must ' //    &
                              'not contain any _FillValues'
             CALL message( 'netcdf_data_input_mod', 'PA0545', 2, 2, 0, 6, 0 )
          ENDIF
       ENDIF

       IF ( init_3d%from_file_w )  THEN
          check_passed = .TRUE.
          IF ( init_3d%lod_w == 1 )  THEN
             IF ( ANY( init_3d%w_init(nzb+1:nzt) == init_3d%fill_w ) )         &
                check_passed = .FALSE.
          ELSEIF ( init_3d%lod_w == 2 )  THEN
             IF ( ANY( w(nzb+1:nzt,nys:nyn,nxl:nxr) == init_3d%fill_w ) )      &
                check_passed = .FALSE.
          ENDIF
          IF ( .NOT. check_passed )  THEN
             message_string = 'NetCDF input for init_atmosphere_w must ' //    &
                              'not contain any _FillValues'
             CALL message( 'netcdf_data_input_mod', 'PA0545', 2, 2, 0, 6, 0 )
          ENDIF
       ENDIF

       IF ( init_3d%from_file_pt )  THEN
          check_passed = .TRUE.
          IF ( init_3d%lod_pt == 1 )  THEN
             IF ( ANY( init_3d%pt_init(nzb+1:nzt+1) == init_3d%fill_pt ) )     &
                check_passed = .FALSE.
          ELSEIF ( init_3d%lod_pt == 2 )  THEN
             IF ( ANY( pt(nzb+1:nzt+1,nys:nyn,nxl:nxr) == init_3d%fill_pt ) )  &
                check_passed = .FALSE.
          ENDIF
          IF ( .NOT. check_passed )  THEN
             message_string = 'NetCDF input for init_atmosphere_pt must ' //   &
                              'not contain any _FillValues'
             CALL message( 'netcdf_data_input_mod', 'PA0545', 2, 2, 0, 6, 0 )
          ENDIF
       ENDIF

       IF ( init_3d%from_file_q )  THEN
          check_passed = .TRUE.
          IF ( init_3d%lod_q == 1 )  THEN
             IF ( ANY( init_3d%q_init(nzb+1:nzt+1) == init_3d%fill_q ) )       &
                check_passed = .FALSE.
          ELSEIF ( init_3d%lod_q == 2 )  THEN
             IF ( ANY( q(nzb+1:nzt+1,nys:nyn,nxl:nxr) == init_3d%fill_q ) )    &
                check_passed = .FALSE.
          ENDIF
          IF ( .NOT. check_passed )  THEN
             message_string = 'NetCDF input for init_atmosphere_q must ' //    &
                              'not contain any _FillValues'
             CALL message( 'netcdf_data_input_mod', 'PA0545', 2, 2, 0, 6, 0 )
          ENDIF
       ENDIF
!
!--    Workaround for cyclic conditions. Please see above for further explanation.
       IF ( bc_lr_cyc  .AND.  nxl == 0 )  nxlu = nxl
       IF ( bc_ns_cyc  .AND.  nys == 0 )  nysv = nys

    END SUBROUTINE netcdf_data_input_init_3d

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Checks input file for consistency and minimum requirements.
!------------------------------------------------------------------------------!
    SUBROUTINE netcdf_data_input_check_dynamic

       USE control_parameters,                                                 &
           ONLY:  initializing_actions, message_string

       IMPLICIT NONE
!
!--    Dynamic input file must also be present if initialization via inifor is
!--    prescribed.
       IF ( .NOT. input_pids_dynamic  .AND.                                    &
            INDEX( initializing_actions, 'inifor' ) /= 0 )  THEN
          message_string = 'initializing_actions = inifor requires dynamic ' //&
                           'input file ' // TRIM( input_file_dynamic ) //      &
                           TRIM( coupling_char )
          CALL message( 'netcdf_data_input_mod', 'PA0547', 1, 2, 0, 6, 0 )
       ENDIF

    END SUBROUTINE netcdf_data_input_check_dynamic

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Checks input file for consistency and minimum requirements.
!------------------------------------------------------------------------------!
    SUBROUTINE netcdf_data_input_check_static

       USE arrays_3d,                                                          &
           ONLY:  zu

       USE control_parameters,                                                 &
           ONLY:  land_surface, message_string, urban_surface

       USE indices,                                                            &
           ONLY:  nxl, nxr, nyn, nys, wall_flags_total_0

       IMPLICIT NONE

       INTEGER(iwp) ::  i      !< loop index along x-direction
       INTEGER(iwp) ::  j      !< loop index along y-direction
       INTEGER(iwp) ::  n_surf !< number of different surface types at given location

       LOGICAL      ::  check_passed !< flag indicating if a check passed

!
!--    Return if no static input file is available
       IF ( .NOT. input_pids_static )  RETURN
!
!--    Check for correct dimension of surface_fractions, should run from 0-2.
       IF ( surface_fraction_f%from_file )  THEN
          IF ( surface_fraction_f%nf-1 > 2 )  THEN
             message_string = 'nsurface_fraction must not be larger than 3.' 
             CALL message( 'netcdf_data_input_mod', 'PA0580', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF
!
!--    Check orography for fill-values. For the moment, give an error message.
!--    More advanced methods, e.g. a nearest neighbor algorithm as used in GIS
!--    systems might be implemented later.
!--    Please note, if no terrain height is provided, it is set to 0.
       IF ( ANY( terrain_height_f%var == terrain_height_f%fill ) )  THEN
          message_string = 'NetCDF variable zt is not ' //                     &
                           'allowed to have missing data'
          CALL message( 'netcdf_data_input_mod', 'PA0550', 2, 2, myid, 6, 0 )
       ENDIF
!
!--    Check for negative terrain heights
       IF ( ANY( terrain_height_f%var < 0.0_wp ) )  THEN
          message_string = 'NetCDF variable zt is not ' //                     &
                           'allowed to have negative values'
          CALL message( 'netcdf_data_input_mod', 'PA0551', 2, 2, myid, 6, 0 )
       ENDIF
!
!--    If 3D buildings are read, check if building information is consistent
!--    to numeric grid.
       IF ( buildings_f%from_file )  THEN
          IF ( buildings_f%lod == 2 )  THEN
             IF ( buildings_f%nz > SIZE( zu ) )  THEN
                message_string = 'Reading 3D building data - too much ' //     &
                                 'data points along the vertical coordinate.'
                CALL message( 'netcdf_data_input_mod', 'PA0552', 2, 2, 0, 6, 0 )
             ENDIF

             IF ( ANY( ABS( buildings_f%z(0:buildings_f%nz-1) -                &
                       zu(0:buildings_f%nz-1) ) > 1E-6_wp ) )  THEN
                message_string = 'Reading 3D building data - vertical ' //     &
                                 'coordinate do not match numeric grid.'
                CALL message( 'netcdf_data_input_mod', 'PA0553', 2, 2, 0, 6, 0 )
             ENDIF
          ENDIF
       ENDIF

!
!--    Skip further checks concerning buildings and natural surface properties
!--    if no urban surface and land surface model are applied.
       IF (  .NOT. land_surface  .AND.  .NOT. urban_surface )  RETURN
!
!--    Check for minimum requirement of surface-classification data in case
!--    static input file is used.
       IF ( ( .NOT. vegetation_type_f%from_file  .OR.                          &
              .NOT. pavement_type_f%from_file    .OR.                          &
              .NOT. water_type_f%from_file       .OR.                          &
              .NOT. soil_type_f%from_file             ) .OR.                   &
             ( urban_surface  .AND.  .NOT. building_type_f%from_file ) )  THEN
          message_string = 'Minimum requirement for surface classification ' //&
                           'is not fulfilled. At least ' //                    &
                           'vegetation_type, pavement_type, ' //               &
                           'soil_type and water_type are '//                   &
                           'required. If urban-surface model is applied, ' //  &
                           'also building_type is required'
          CALL message( 'netcdf_data_input_mod', 'PA0554', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Check for general availability of input variables.
!--    If vegetation_type is 0 at any location, vegetation_pars as well as
!--    root_area_dens_s are required.
       IF ( vegetation_type_f%from_file )  THEN
          IF ( ANY( vegetation_type_f%var == 0 ) )  THEN
             IF ( .NOT. vegetation_pars_f%from_file )  THEN
                message_string = 'If vegetation_type = 0 at any location, ' // &
                                 'vegetation_pars is required'
                CALL message( 'netcdf_data_input_mod', 'PA0555', 2, 2, myid, 6, 0 )
             ENDIF
             IF ( .NOT. root_area_density_lsm_f%from_file )  THEN
                message_string = 'If vegetation_type = 0 at any location, ' // &
                                 'root_area_dens_s is required'
                CALL message( 'netcdf_data_input_mod', 'PA0556', 2, 2, myid, 6, 0 )
             ENDIF
          ENDIF
       ENDIF
!
!--    If soil_type is zero at any location, soil_pars is required.
       IF ( soil_type_f%from_file )  THEN
          check_passed = .TRUE.
          IF ( ALLOCATED( soil_type_f%var_2d ) )  THEN
             IF ( ANY( soil_type_f%var_2d == 0 ) )  THEN
                IF ( .NOT. soil_pars_f%from_file )  check_passed = .FALSE.
             ENDIF
          ELSE
             IF ( ANY( soil_type_f%var_3d == 0 ) )  THEN
                IF ( .NOT. soil_pars_f%from_file )  check_passed = .FALSE.
             ENDIF
          ENDIF
          IF ( .NOT. check_passed )  THEN
             message_string = 'If soil_type = 0 at any location, ' //          &
                              'soil_pars is required'
             CALL message( 'netcdf_data_input_mod', 'PA0557', 2, 2, myid, 6, 0 )
          ENDIF
       ENDIF
!
!--    Buildings require a type in case of urban-surface model.
       IF ( buildings_f%from_file  .AND.  .NOT. building_type_f%from_file  )  THEN
          message_string = 'If buildings are provided, also building_type ' // &
                           'is required'
          CALL message( 'netcdf_data_input_mod', 'PA0581', 2, 2, 0, 6, 0 )
       ENDIF
!
!--    Buildings require an ID.
       IF ( buildings_f%from_file  .AND.  .NOT. building_id_f%from_file  )  THEN
          message_string = 'If buildings are provided, also building_id ' //   &
                           'is required'
          CALL message( 'netcdf_data_input_mod', 'PA0582', 2, 2, 0, 6, 0 )
       ENDIF
!
!--    If building_type is zero at any location, building_pars is required.
       IF ( building_type_f%from_file )  THEN
          IF ( ANY( building_type_f%var == 0 ) )  THEN
             IF ( .NOT. building_pars_f%from_file )  THEN
                message_string = 'If building_type = 0 at any location, ' //   &
                                 'building_pars is required'
                CALL message( 'netcdf_data_input_mod', 'PA0558', 2, 2, myid, 6, 0 )
             ENDIF
          ENDIF
       ENDIF
!
!--    If building_type is provided, also building_id is needed (due to the 
!--    filtering algorithm).
       IF ( building_type_f%from_file  .AND.  .NOT. building_id_f%from_file )  &
       THEN
          message_string = 'If building_type is provided, also building_id '// &
                           'is required'
          CALL message( 'netcdf_data_input_mod', 'PA0519', 2, 2, 0, 6, 0 )
       ENDIF       
!
!--    If albedo_type is zero at any location, albedo_pars is required.
       IF ( albedo_type_f%from_file )  THEN
          IF ( ANY( albedo_type_f%var == 0 ) )  THEN
             IF ( .NOT. albedo_pars_f%from_file )  THEN
                message_string = 'If albedo_type = 0 at any location, ' //     &
                                 'albedo_pars is required'
                CALL message( 'netcdf_data_input_mod', 'PA0559', 2, 2, myid, 6, 0 )
             ENDIF
          ENDIF
       ENDIF
!
!--    If pavement_type is zero at any location, pavement_pars is required.
       IF ( pavement_type_f%from_file )  THEN
          IF ( ANY( pavement_type_f%var == 0 ) )  THEN
             IF ( .NOT. pavement_pars_f%from_file )  THEN
                message_string = 'If pavement_type = 0 at any location, ' //   &
                                 'pavement_pars is required'
                CALL message( 'netcdf_data_input_mod', 'PA0560', 2, 2, myid, 6, 0 )
             ENDIF
          ENDIF
       ENDIF
!
!--    If pavement_type is zero at any location, also pavement_subsurface_pars
!--    is required.
       IF ( pavement_type_f%from_file )  THEN
          IF ( ANY( pavement_type_f%var == 0 ) )  THEN
             IF ( .NOT. pavement_subsurface_pars_f%from_file )  THEN
                message_string = 'If pavement_type = 0 at any location, ' //   &
                                 'pavement_subsurface_pars is required'
                CALL message( 'netcdf_data_input_mod', 'PA0561', 2, 2, myid, 6, 0 )
             ENDIF
          ENDIF
       ENDIF
!
!--    If water_type is zero at any location, water_pars is required.
       IF ( water_type_f%from_file )  THEN
          IF ( ANY( water_type_f%var == 0 ) )  THEN
             IF ( .NOT. water_pars_f%from_file )  THEN
                message_string = 'If water_type = 0 at any location, ' //      &
                                 'water_pars is required'
                CALL message( 'netcdf_data_input_mod', 'PA0562', 2, 2,myid, 6, 0 )
             ENDIF
          ENDIF
       ENDIF
!
!--    Check for local consistency of the input data.
       DO  i = nxl, nxr
          DO  j = nys, nyn
!
!--          For each (y,x)-location at least one of the parameters
!--          vegetation_type, pavement_type, building_type, or water_type
!--          must be set to a non­missing value.
             IF ( land_surface  .AND.  .NOT. urban_surface )  THEN
                IF ( vegetation_type_f%var(j,i) == vegetation_type_f%fill  .AND.&
                     pavement_type_f%var(j,i)   == pavement_type_f%fill    .AND.&
                     water_type_f%var(j,i)      == water_type_f%fill )  THEN
                   WRITE( message_string, * )                                  &
                                    'At least one of the parameters '//        &
                                    'vegetation_type, pavement_type, '     //  &
                                    'or water_type must be set '//             &
                                    'to a non-missing value. Grid point: ', j, i
                   CALL message( 'netcdf_data_input_mod', 'PA0563', 2, 2, myid, 6, 0 )
                ENDIF
             ELSEIF ( land_surface  .AND.  urban_surface )  THEN
                IF ( vegetation_type_f%var(j,i) == vegetation_type_f%fill  .AND.&
                     pavement_type_f%var(j,i)   == pavement_type_f%fill    .AND.&
                     building_type_f%var(j,i)   == building_type_f%fill    .AND.&
                     water_type_f%var(j,i)      == water_type_f%fill )  THEN
                   WRITE( message_string, * )                                  &
                                 'At least one of the parameters '//           &
                                 'vegetation_type, pavement_type, '  //        &
                                 'building_type, or water_type must be set '// &
                                 'to a non-missing value. Grid point: ', j, i
                   CALL message( 'netcdf_data_input_mod', 'PA0563', 2, 2, myid, 6, 0 )
                ENDIF
             ENDIF
                
!
!--          Note that a soil_type is required for each location (y,x) where
!--          either vegetation_type or pavement_type is a non­missing value.
             IF ( ( vegetation_type_f%var(j,i) /= vegetation_type_f%fill  .OR. &
                    pavement_type_f%var(j,i)   /= pavement_type_f%fill ) )  THEN
                check_passed = .TRUE.
                IF ( ALLOCATED( soil_type_f%var_2d ) )  THEN
                   IF ( soil_type_f%var_2d(j,i) == soil_type_f%fill )          &
                      check_passed = .FALSE.
                ELSE
                   IF ( ANY( soil_type_f%var_3d(:,j,i) == soil_type_f%fill) )  &
                      check_passed = .FALSE.
                ENDIF

                IF ( .NOT. check_passed )  THEN
                   message_string = 'soil_type is required for each '//        &
                                 'location (y,x) where vegetation_type or ' // &
                                 'pavement_type is a non-missing value.'
                   CALL message( 'netcdf_data_input_mod', 'PA0564',            &
                                  2, 2, myid, 6, 0 )
                ENDIF
             ENDIF 
!
!--          Check for consistency of given types. At the moment, only one 
!--          of vegetation, pavement, or water-type can be set. This is 
!--          because no tile approach is yet implemented in the land-surface
!--          model. Later, when this is possible, surface fraction need to be 
!--          given and the sum must not  be larger than 1. Please note, in case
!--          more than one type is given at a pixel, an error message will be 
!--          given. 
             n_surf = 0
             IF ( vegetation_type_f%var(j,i) /= vegetation_type_f%fill )       &
                n_surf = n_surf + 1
             IF ( water_type_f%var(j,i)      /= water_type_f%fill )            &
                n_surf = n_surf + 1
             IF ( pavement_type_f%var(j,i)   /= pavement_type_f%fill )         &
                n_surf = n_surf + 1

             IF ( n_surf > 1 )  THEN
                WRITE( message_string, * )                                     &
                                 'More than one surface type (vegetation, '//  &
                                 'pavement, water) is given at a location. '// &
                                 'Please note, this is not possible at ' //    &
                                 'the moment as no tile approach has been ' // &
                                 'yet implemented. (i,j) = ', i, j
                CALL message( 'netcdf_data_input_mod', 'PA0565',               &
                               2, 2, myid, 6, 0 )

!                 IF ( .NOT. surface_fraction_f%from_file )  THEN
!                    message_string = 'More than one surface type (vegetation '//&
!                                  'pavement, water) is given at a location. '// &
!                                  'Please note, this is not possible at ' //    &
!                                  'the moment as no tile approach is yet ' //   &
!                                  'implemented.'
!                    message_string = 'If more than one surface type is ' //     &
!                                  'given at a location, surface_fraction ' //   &
!                                  'must be provided.'
!                    CALL message( 'netcdf_data_input_mod', 'PA0565',            &
!                                   2, 2, myid, 6, 0 )
!                 ELSEIF ( ANY ( surface_fraction_f%frac(:,j,i) ==               &
!                                surface_fraction_f%fill ) )  THEN
!                    message_string = 'If more than one surface type is ' //     &
!                                  'given at a location, surface_fraction ' //   &
!                                  'must be provided.'
!                    CALL message( 'netcdf_data_input_mod', 'PA0565',            &
!                                   2, 2, myid, 6, 0 )
!                 ENDIF
             ENDIF
!
!--          Check for further mismatches. e.g. relative fractions exceed 1 or
!--          vegetation_type is set but surface vegetation fraction is zero, 
!--          etc..
             IF ( surface_fraction_f%from_file )  THEN
!
!--             If surface fractions is given, also check that only one type 
!--             is given.
                IF ( SUM( MERGE( 1, 0, surface_fraction_f%frac(:,j,i) /= 0.0_wp&
                                .AND.  surface_fraction_f%frac(:,j,i) /=       &
                                       surface_fraction_f%fill  ) ) > 1 )  THEN
                   WRITE( message_string, * )                                  &
                                    'surface_fraction is given for more ' //   &
                                    'than one type. ' //                       &
                                    'Please note, this is not possible at ' // &
                                    'the moment as no tile approach has '//    &
                                    'yet been implemented. (i, j) = ', i, j
                   CALL message( 'netcdf_data_input_mod', 'PA0676',            &
                                  2, 2, myid, 6, 0 )
                ENDIF
!
!--             Sum of relative fractions must be 1. Note, attributed to type 
!--             conversions due to reading, the sum of surface fractions 
!--             might be not exactly 1. Hence, the sum is check with a
!--             tolerance. Later, in the land-surface model, the relative 
!--             fractions are normalized to one. Actually, surface fractions
!--             shall be _FillValue at building grid points, however, in order
!--             to relax this requirement and allow that surface-fraction can
!--             also be zero at these grid points, only perform this check
!--             at locations where some vegetation, pavement or water is defined.
                IF ( vegetation_type_f%var(j,i) /= vegetation_type_f%fill  .OR.&
                     pavement_type_f%var(j,i)   /= pavement_type_f%fill    .OR.&
                     water_type_f%var(j,i)      /= water_type_f%fill )  THEN
                   IF ( SUM ( surface_fraction_f%frac(0:2,j,i) ) >             &
                        1.0_wp + 1E-8_wp  .OR.                                 &
                        SUM ( surface_fraction_f%frac(0:2,j,i) ) <             &
                        1.0_wp - 1E-8_wp )  THEN
                      WRITE( message_string, * )                               &
                                    'The sum of all land-surface fractions ' //&
                                    'must equal 1. (i, j) = ', i, j
                      CALL message( 'netcdf_data_input_mod', 'PA0566',         &
                                     2, 2, myid, 6, 0 )
                   ENDIF
                ENDIF
!
!--             Relative fraction for a type must not be zero at locations where 
!--             this type is set. 
                IF (                                                           &
                  ( vegetation_type_f%var(j,i) /= vegetation_type_f%fill  .AND.&
                 ( surface_fraction_f%frac(ind_veg_wall,j,i) == 0.0_wp .OR.    &
                   surface_fraction_f%frac(ind_veg_wall,j,i) ==                &
                                                     surface_fraction_f%fill ) &
                  )  .OR.                                                      &
                  ( pavement_type_f%var(j,i) /= pavement_type_f%fill     .AND. &
                 ( surface_fraction_f%frac(ind_pav_green,j,i) == 0.0_wp .OR.   &
                   surface_fraction_f%frac(ind_pav_green,j,i) ==               &
                                                     surface_fraction_f%fill ) &
                  )  .OR.                                                      &
                  ( water_type_f%var(j,i) /= water_type_f%fill           .AND. &
                 ( surface_fraction_f%frac(ind_wat_win,j,i) == 0.0_wp .OR.     &
                   surface_fraction_f%frac(ind_wat_win,j,i) ==                 &
                                                     surface_fraction_f%fill ) &
                  ) )  THEN
                   WRITE( message_string, * ) 'Mismatch in setting of '     // &
                             'surface_fraction. Vegetation-, pavement-, or '// &
                             'water surface is given at (i,j) = ( ', i, j,     &
                             ' ), but surface fraction is 0 for the given type.'
                   CALL message( 'netcdf_data_input_mod', 'PA0567',            &
                                  2, 2, myid, 6, 0 )
                ENDIF
!
!--             Relative fraction for a type must not contain non-zero values
!--             if this type is not set. 
                IF (                                                           &
                  ( vegetation_type_f%var(j,i) == vegetation_type_f%fill  .AND.&
                 ( surface_fraction_f%frac(ind_veg_wall,j,i) /= 0.0_wp .AND.   &
                   surface_fraction_f%frac(ind_veg_wall,j,i) /=                &
                                                     surface_fraction_f%fill ) &
                  )  .OR.                                                      &
                  ( pavement_type_f%var(j,i) == pavement_type_f%fill     .AND. &
                 ( surface_fraction_f%frac(ind_pav_green,j,i) /= 0.0_wp .AND.  &
                   surface_fraction_f%frac(ind_pav_green,j,i) /=               &
                                                     surface_fraction_f%fill ) &
                  )  .OR.                                                      &
                  ( water_type_f%var(j,i) == water_type_f%fill           .AND. &
                 ( surface_fraction_f%frac(ind_wat_win,j,i) /= 0.0_wp .AND.    &
                   surface_fraction_f%frac(ind_wat_win,j,i) /=                 &
                                                     surface_fraction_f%fill ) &
                  ) )  THEN
                   WRITE( message_string, * ) 'Mismatch in setting of '     // &
                             'surface_fraction. Vegetation-, pavement-, or '// &
                             'water surface is not given at (i,j) = ( ', i, j, &
                             ' ), but surface fraction is not 0 for the ' //   &
                             'given type.'
                   CALL message( 'netcdf_data_input_mod', 'PA0568',            &
                                  2, 2, myid, 6, 0 )
                ENDIF
             ENDIF
!
!--          Check vegetation_pars. If vegetation_type is 0, all parameters
!--          need to be set, otherwise, single parameters set by
!--          vegetation_type can be overwritten.
             IF ( vegetation_type_f%from_file )  THEN
                IF ( vegetation_type_f%var(j,i) == 0 )  THEN
                   IF ( ANY( vegetation_pars_f%pars_xy(:,j,i) ==               &
                             vegetation_pars_f%fill ) )  THEN
                      message_string = 'If vegetation_type(y,x) = 0, all '  // &
                                       'parameters of vegetation_pars at '//   &
                                       'this location must be set.'
                      CALL message( 'netcdf_data_input_mod', 'PA0569',         &
                                     2, 2, myid, 6, 0 )
                   ENDIF
                ENDIF
             ENDIF
!
!--          Check root distribution. If vegetation_type is 0, all levels must
!--          be set.
             IF ( vegetation_type_f%from_file )  THEN
                IF ( vegetation_type_f%var(j,i) == 0 )  THEN
                   IF ( ANY( root_area_density_lsm_f%var(:,j,i) ==             &
                             root_area_density_lsm_f%fill ) )  THEN
                      message_string = 'If vegetation_type(y,x) = 0, all ' //  &
                                       'levels of root_area_dens_s ' //        &
                                       'must be set at this location.'
                      CALL message( 'netcdf_data_input_mod', 'PA0570',         &
                                     2, 2, myid, 6, 0 )
                   ENDIF
                ENDIF
             ENDIF
!
!--          Check soil parameters. If soil_type is 0, all parameters
!--          must be set.
             IF ( soil_type_f%from_file )  THEN
                check_passed = .TRUE.
                IF ( ALLOCATED( soil_type_f%var_2d ) )  THEN
                   IF ( soil_type_f%var_2d(j,i) == 0 )  THEN
                      IF ( ANY( soil_pars_f%pars_xy(:,j,i) ==                  &
                                soil_pars_f%fill ) )  check_passed = .FALSE.
                   ENDIF
                ELSE
                   IF ( ANY( soil_type_f%var_3d(:,j,i) == 0 ) )  THEN
                      IF ( ANY( soil_pars_f%pars_xy(:,j,i) ==                  &
                                soil_pars_f%fill ) )  check_passed = .FALSE.
                   ENDIF
                ENDIF
                IF ( .NOT. check_passed )  THEN
                   message_string = 'If soil_type(y,x) = 0, all levels of '  //&
                                    'soil_pars at this location must be set.'
                   CALL message( 'netcdf_data_input_mod', 'PA0571',            &
                                  2, 2, myid, 6, 0 )
                ENDIF
             ENDIF

!
!--          Check building parameters. If building_type is 0, all parameters
!--          must be set.
             IF ( building_type_f%from_file )  THEN
                IF ( building_type_f%var(j,i) == 0 )  THEN
                   IF ( ANY( building_pars_f%pars_xy(:,j,i) ==                 &
                             building_pars_f%fill ) )  THEN
                      message_string = 'If building_type(y,x) = 0, all ' //    &
                                       'parameters of building_pars at this '//&
                                       'location must be set.'
                      CALL message( 'netcdf_data_input_mod', 'PA0572',         &
                                     2, 2, myid, 6, 0 )
                   ENDIF
                ENDIF
             ENDIF
!
!--          Check if building_type is set at each building and vice versa. 
!--          Please note, buildings are already processed and filtered.
!--          For this reason, consistency checks are based on wall_flags_total_0 
!--          rather than buildings_f (buildings are represented by bit 6 in 
!--          wall_flags_total_0). 
             IF ( building_type_f%from_file  .AND.  buildings_f%from_file )  THEN
                IF ( ANY( BTEST ( wall_flags_total_0(:,j,i), 6 ) )  .AND.      &
                     building_type_f%var(j,i) == building_type_f%fill  .OR.    &
               .NOT. ANY( BTEST ( wall_flags_total_0(:,j,i), 6 ) )  .AND.      &
                     building_type_f%var(j,i) /= building_type_f%fill )  THEN
                   WRITE( message_string, * ) 'Each location where a ' //      &
                                   'building is set requires a type ' //       &
                                   '( and vice versa ) in case the ' //        &
                                   'urban-surface model is applied. ' //       &
                                   'i, j = ', i, j
                   CALL message( 'netcdf_data_input_mod', 'PA0573',            &
                                  2, 2, myid, 6, 0 )
                ENDIF
             ENDIF
!
!--          Check if at each location where a building is present also an ID
!--          is set and vice versa.
             IF ( buildings_f%from_file )  THEN
                IF ( ANY( BTEST ( wall_flags_total_0(:,j,i), 6 ) )  .AND.     &
                     building_id_f%var(j,i) == building_id_f%fill  .OR.       &
               .NOT. ANY( BTEST ( wall_flags_total_0(:,j,i), 6 ) )  .AND.     &
                     building_id_f%var(j,i) /= building_id_f%fill )  THEN
                   WRITE( message_string, * ) 'Each location where a ' //     &
                                   'building is set requires an ID ' //       &
                                   '( and vice versa ). i, j = ', i, j
                   CALL message( 'netcdf_data_input_mod', 'PA0574',           &
                                  2, 2, myid, 6, 0 )
                ENDIF
             ENDIF
!
!--          Check if building ID is set where a bulding is defined.
             IF ( buildings_f%from_file )  THEN
                IF ( ANY( BTEST ( wall_flags_total_0(:,j,i), 6 ) )  .AND.     &
                     building_id_f%var(j,i) == building_id_f%fill )  THEN
                   WRITE( message_string, * ) 'Each building grid point '//   &
                                              'requires an ID.', i, j
                   CALL message( 'netcdf_data_input_mod', 'PA0575',           &
                                  2, 2, myid, 6, 0 )
                ENDIF
             ENDIF
!
!--          Check albedo parameters. If albedo_type is 0, all parameters
!--          must be set.
             IF ( albedo_type_f%from_file )  THEN
                IF ( albedo_type_f%var(j,i) == 0 )  THEN
                   IF ( ANY( albedo_pars_f%pars_xy(:,j,i) ==                   &
                             albedo_pars_f%fill ) )  THEN
                      message_string = 'If albedo_type(y,x) = 0, all ' //      &
                                       'parameters of albedo_pars at this ' // &
                                       'location must be set.'
                      CALL message( 'netcdf_data_input_mod', 'PA0576',         &
                                     2, 2, myid, 6, 0 )
                   ENDIF
                ENDIF
             ENDIF

!
!--          Check pavement parameters. If pavement_type is 0, all parameters
!--          of pavement_pars must be set at this location.
             IF ( pavement_type_f%from_file )  THEN
                IF ( pavement_type_f%var(j,i) == 0 )  THEN
                   IF ( ANY( pavement_pars_f%pars_xy(:,j,i) ==                 &
                             pavement_pars_f%fill ) )  THEN
                      message_string = 'If pavement_type(y,x) = 0, all ' //    &
                                       'parameters of pavement_pars at this '//&
                                       'location must be set.'
                      CALL message( 'netcdf_data_input_mod', 'PA0577',         &
                                     2, 2, myid, 6, 0 )
                   ENDIF
                ENDIF
             ENDIF
!
!--          Check pavement-subsurface parameters. If pavement_type is 0,
!--          all parameters of pavement_subsurface_pars must be set  at this
!--          location.
             IF ( pavement_type_f%from_file )  THEN
                IF ( pavement_type_f%var(j,i) == 0 )  THEN
                   IF ( ANY( pavement_subsurface_pars_f%pars_xyz(:,:,j,i) ==   &
                             pavement_subsurface_pars_f%fill ) )  THEN
                      message_string = 'If pavement_type(y,x) = 0, all ' //    &
                                       'parameters of '                  //    &
                                       'pavement_subsurface_pars at this '//   &
                                       'location must be set.'
                      CALL message( 'netcdf_data_input_mod', 'PA0578',         &
                                     2, 2, myid, 6, 0 )
                   ENDIF
                ENDIF
             ENDIF

!
!--          Check water parameters. If water_type is 0, all parameters
!--          must be set  at this location.
             IF ( water_type_f%from_file )  THEN
                IF ( water_type_f%var(j,i) == 0 )  THEN
                   IF ( ANY( water_pars_f%pars_xy(:,j,i) ==                    &
                             water_pars_f%fill ) )  THEN
                      message_string = 'If water_type(y,x) = 0, all ' //       &
                                       'parameters of water_pars at this ' //  &
                                       'location must be set.'
                      CALL message( 'netcdf_data_input_mod', 'PA0579',         &
                                     2, 2, myid, 6, 0 )
                   ENDIF
                ENDIF
             ENDIF

          ENDDO
       ENDDO

    END SUBROUTINE netcdf_data_input_check_static

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Resize 8-bit 2D Integer array: (nys:nyn,nxl:nxr) -> (nysg:nyng,nxlg:nxrg)
!------------------------------------------------------------------------------!
    SUBROUTINE resize_array_2d_int8( var, js, je, is, ie )
    
       IMPLICIT NONE

       INTEGER(iwp) ::  je !< upper index bound along y direction
       INTEGER(iwp) ::  js !< lower index bound along y direction
       INTEGER(iwp) ::  ie !< upper index bound along x direction
       INTEGER(iwp) ::  is !< lower index bound along x direction
       
       INTEGER(KIND=1), DIMENSION(:,:), ALLOCATABLE ::  var     !< treated variable
       INTEGER(KIND=1), DIMENSION(:,:), ALLOCATABLE ::  var_tmp !< temporary copy
!
!--    Allocate temporary variable
       ALLOCATE( var_tmp(js-nbgp:je+nbgp,is-nbgp:ie+nbgp) )
!
!--    Temporary copy of the variable
       var_tmp(js:je,is:ie) = var(js:je,is:ie)
!
!--    Resize the array
       DEALLOCATE( var )
       ALLOCATE( var(js-nbgp:je+nbgp,is-nbgp:ie+nbgp) )
!
!--    Transfer temporary copy back to original array
       var(js:je,is:ie) = var_tmp(js:je,is:ie)

    END SUBROUTINE resize_array_2d_int8
    
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Resize 32-bit 2D Integer array: (nys:nyn,nxl:nxr) -> (nysg:nyng,nxlg:nxrg)
!------------------------------------------------------------------------------!
    SUBROUTINE resize_array_2d_int32( var, js, je, is, ie )

       IMPLICIT NONE
       
       INTEGER(iwp) ::  je !< upper index bound along y direction
       INTEGER(iwp) ::  js !< lower index bound along y direction
       INTEGER(iwp) ::  ie !< upper index bound along x direction
       INTEGER(iwp) ::  is !< lower index bound along x direction

       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  var     !< treated variable
       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  var_tmp !< temporary copy
!
!--    Allocate temporary variable
       ALLOCATE( var_tmp(js-nbgp:je+nbgp,is-nbgp:ie+nbgp) )
!
!--    Temporary copy of the variable
       var_tmp(js:je,is:ie) = var(js:je,is:ie)
!
!--    Resize the array
       DEALLOCATE( var )
       ALLOCATE( var(js-nbgp:je+nbgp,is-nbgp:ie+nbgp) )
!
!--    Transfer temporary copy back to original array
       var(js:je,is:ie) = var_tmp(js:je,is:ie)

    END SUBROUTINE resize_array_2d_int32
    
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Resize 8-bit 3D Integer array: (:,nys:nyn,nxl:nxr) -> (:,nysg:nyng,nxlg:nxrg)
!------------------------------------------------------------------------------!
    SUBROUTINE resize_array_3d_int8( var, ks, ke, js, je, is, ie )

       IMPLICIT NONE

       INTEGER(iwp) ::  je !< upper index bound along y direction
       INTEGER(iwp) ::  js !< lower index bound along y direction
       INTEGER(iwp) ::  ie !< upper index bound along x direction
       INTEGER(iwp) ::  is !< lower index bound along x direction
       INTEGER(iwp) ::  ke !< upper bound of treated array in z-direction  
       INTEGER(iwp) ::  ks !< lower bound of treated array in z-direction  
       
       INTEGER(KIND=1), DIMENSION(:,:,:), ALLOCATABLE ::  var     !< treated variable
       INTEGER(KIND=1), DIMENSION(:,:,:), ALLOCATABLE ::  var_tmp !< temporary copy
!
!--    Allocate temporary variable
       ALLOCATE( var_tmp(ks:ke,js-nbgp:je+nbgp,is-nbgp:ie+nbgp) )
!
!--    Temporary copy of the variable
       var_tmp(ks:ke,js:je,is:ie) = var(ks:ke,js:je,is:ie)
!
!--    Resize the array
       DEALLOCATE( var )
       ALLOCATE( var(ks:ke,js-nbgp:je+nbgp,is-nbgp:ie+nbgp) )
!
!--    Transfer temporary copy back to original array
       var(ks:ke,js:je,is:ie) = var_tmp(ks:ke,js:je,is:ie)

    END SUBROUTINE resize_array_3d_int8
    
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Resize 3D Real array: (:,nys:nyn,nxl:nxr) -> (:,nysg:nyng,nxlg:nxrg)
!------------------------------------------------------------------------------!
    SUBROUTINE resize_array_3d_real( var, ks, ke, js, je, is, ie )

       IMPLICIT NONE

       INTEGER(iwp) ::  je !< upper index bound along y direction
       INTEGER(iwp) ::  js !< lower index bound along y direction
       INTEGER(iwp) ::  ie !< upper index bound along x direction
       INTEGER(iwp) ::  is !< lower index bound along x direction
       INTEGER(iwp) ::  ke !< upper bound of treated array in z-direction  
       INTEGER(iwp) ::  ks !< lower bound of treated array in z-direction  
       
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  var     !< treated variable
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  var_tmp !< temporary copy
!
!--    Allocate temporary variable
       ALLOCATE( var_tmp(ks:ke,js-nbgp:je+nbgp,is-nbgp:ie+nbgp) )
!
!--    Temporary copy of the variable
       var_tmp(ks:ke,js:je,is:ie) = var(ks:ke,js:je,is:ie)
!
!--    Resize the array
       DEALLOCATE( var )
       ALLOCATE( var(ks:ke,js-nbgp:je+nbgp,is-nbgp:ie+nbgp) )
!
!--    Transfer temporary copy back to original array
       var(ks:ke,js:je,is:ie) = var_tmp(ks:ke,js:je,is:ie)

    END SUBROUTINE resize_array_3d_real
    
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Resize 4D Real array: (:,:,nys:nyn,nxl:nxr) -> (:,nysg:nyng,nxlg:nxrg)
!------------------------------------------------------------------------------!
    SUBROUTINE resize_array_4d_real( var, k1s, k1e, k2s, k2e, js, je, is, ie )

       IMPLICIT NONE
       
       INTEGER(iwp) ::  je  !< upper index bound along y direction
       INTEGER(iwp) ::  js  !< lower index bound along y direction
       INTEGER(iwp) ::  ie  !< upper index bound along x direction
       INTEGER(iwp) ::  is  !< lower index bound along x direction
       INTEGER(iwp) ::  k1e !< upper bound of treated array in z-direction  
       INTEGER(iwp) ::  k1s !< lower bound of treated array in z-direction
       INTEGER(iwp) ::  k2e !< upper bound of treated array along parameter space  
       INTEGER(iwp) ::  k2s !< lower bound of treated array along parameter space  
       
       REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  var     !< treated variable
       REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  var_tmp !< temporary copy
!
!--    Allocate temporary variable
       ALLOCATE( var_tmp(k1s:k1e,k2s:k2e,js-nbgp:je+nbgp,is-nbgp:ie+nbgp) )
!
!--    Temporary copy of the variable
       var_tmp(k1s:k1e,k2s:k2e,js:je,is:ie) = var(k1s:k1e,k2s:k2e,js:je,is:ie)
!
!--    Resize the array
       DEALLOCATE( var )
       ALLOCATE( var(k1s:k1e,k2s:k2e,js-nbgp:je+nbgp,is-nbgp:ie+nbgp) )
!
!--    Transfer temporary copy back to original array
       var(k1s:k1e,k2s:k2e,js:je,is:ie) = var_tmp(k1s:k1e,k2s:k2e,js:je,is:ie)

    END SUBROUTINE resize_array_4d_real

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Checks if a given variables is on file
!------------------------------------------------------------------------------!
    FUNCTION check_existence( vars_in_file, var_name )

       IMPLICIT NONE

       CHARACTER(LEN=*) ::  var_name                   !< variable to be checked
       CHARACTER(LEN=*), DIMENSION(:) ::  vars_in_file !< list of variables in file

       INTEGER(iwp) ::  i                              !< loop variable

       LOGICAL ::  check_existence                     !< flag indicating whether a variable exist or not - actual return value

       i = 1
       check_existence = .FALSE.
       DO  WHILE ( i <= SIZE( vars_in_file ) )
          check_existence = TRIM( vars_in_file(i) ) == TRIM( var_name )  .OR.  &
                            check_existence
          i = i + 1
       ENDDO

       RETURN

    END FUNCTION check_existence


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Closes an existing netCDF file.
!------------------------------------------------------------------------------!
    SUBROUTINE close_input_file( id )

       USE pegrid

       IMPLICIT NONE

       INTEGER(iwp), INTENT(INOUT)        ::  id        !< file id

#if defined( __netcdf )
       nc_stat = NF90_CLOSE( id )
       CALL handle_error( 'close', 540 )
#endif
    END SUBROUTINE close_input_file

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Opens an existing netCDF file for reading only and returns its id.
!------------------------------------------------------------------------------!
    SUBROUTINE open_read_file( filename, id )

       USE pegrid

       IMPLICIT NONE

       CHARACTER (LEN=*), INTENT(IN) ::  filename  !< filename
       INTEGER(iwp), INTENT(INOUT)   ::  id        !< file id

#if defined( __netcdf )

#if defined( __netcdf4_parallel )
!
!--    If __netcdf4_parallel is defined, parrallel NetCDF will be used.
       nc_stat = NF90_OPEN( filename, IOR( NF90_NOWRITE, NF90_MPIIO ), id,     &
                            COMM = comm2d, INFO = MPI_INFO_NULL )
!
!--    In case the previous open call fails, check for possible Netcdf 3 file,
!--    and open it. However, this case, disable parallel access. 
       IF( nc_stat /= NF90_NOERR )  THEN
          nc_stat = NF90_OPEN( filename, NF90_NOWRITE, id )
          collective_read = .FALSE.
       ELSE
#if defined( __nec )
          collective_read = .FALSE.   ! collective read causes hang situations on NEC Aurora
#else
          collective_read = .TRUE.
#endif
       ENDIF
#else
!
!--    All MPI processes open the file and read it (but not in parallel).
       nc_stat = NF90_OPEN( filename, NF90_NOWRITE, id )
#endif

       CALL handle_error( 'open_read_file', 539 )

#endif
    END SUBROUTINE open_read_file

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads global or variable-related attributes of type INTEGER (32-bit)
!------------------------------------------------------------------------------!
     SUBROUTINE get_attribute_int32( id, attribute_name, value, global,        &
                                     variable_name )

       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=*)            ::  attribute_name   !< attribute name
       CHARACTER(LEN=*), OPTIONAL  ::  variable_name    !< variable name

       INTEGER(iwp), INTENT(IN)    ::  id               !< file id
       INTEGER(iwp)                ::  id_var           !< variable id
       INTEGER(iwp), INTENT(INOUT) ::  value            !< read value

       LOGICAL, INTENT(IN) ::  global                   !< flag indicating global attribute
#if defined( __netcdf )

!
!--    Read global attribute
       IF ( global )  THEN
          nc_stat = NF90_GET_ATT( id, NF90_GLOBAL, TRIM( attribute_name ), value )
          CALL handle_error( 'get_attribute_int32 global', 522, attribute_name )
!
!--    Read attributes referring to a single variable. Therefore, first inquire
!--    variable id
       ELSE
          nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
          CALL handle_error( 'get_attribute_int32', 522, attribute_name )
          nc_stat = NF90_GET_ATT( id, id_var, TRIM( attribute_name ), value )
          CALL handle_error( 'get_attribute_int32', 522, attribute_name )
       ENDIF
#endif
    END SUBROUTINE get_attribute_int32

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads global or variable-related attributes of type INTEGER (8-bit)
!------------------------------------------------------------------------------!
     SUBROUTINE get_attribute_int8( id, attribute_name, value, global,         &
                                    variable_name )

       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=*)            ::  attribute_name   !< attribute name
       CHARACTER(LEN=*), OPTIONAL  ::  variable_name    !< variable name

       INTEGER(iwp), INTENT(IN)    ::  id               !< file id
       INTEGER(iwp)                ::  id_var           !< variable id
       INTEGER(KIND=1), INTENT(INOUT) ::  value         !< read value

       LOGICAL, INTENT(IN) ::  global                   !< flag indicating global attribute
#if defined( __netcdf )

!
!--    Read global attribute
       IF ( global )  THEN
          nc_stat = NF90_GET_ATT( id, NF90_GLOBAL, TRIM( attribute_name ), value )
          CALL handle_error( 'get_attribute_int8 global', 523, attribute_name )
!
!--    Read attributes referring to a single variable. Therefore, first inquire
!--    variable id
       ELSE
          nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
          CALL handle_error( 'get_attribute_int8', 523, attribute_name )
          nc_stat = NF90_GET_ATT( id, id_var, TRIM( attribute_name ), value )
          CALL handle_error( 'get_attribute_int8', 523, attribute_name )
       ENDIF
#endif
    END SUBROUTINE get_attribute_int8

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads global or variable-related attributes of type REAL
!------------------------------------------------------------------------------!
     SUBROUTINE get_attribute_real( id, attribute_name, value, global,         &
                                    variable_name )

       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=*)            ::  attribute_name   !< attribute name
       CHARACTER(LEN=*), OPTIONAL  ::  variable_name    !< variable name

       INTEGER(iwp), INTENT(IN)    ::  id               !< file id
       INTEGER(iwp)                ::  id_var           !< variable id

       LOGICAL, INTENT(IN) ::  global                   !< flag indicating global attribute

       REAL(wp), INTENT(INOUT)     ::  value            !< read value
#if defined( __netcdf )


!
!-- Read global attribute
       IF ( global )  THEN
          nc_stat = NF90_GET_ATT( id, NF90_GLOBAL, TRIM( attribute_name ), value )
          CALL handle_error( 'get_attribute_real global', 524, attribute_name )
!
!-- Read attributes referring to a single variable. Therefore, first inquire
!-- variable id
       ELSE
          nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
          CALL handle_error( 'get_attribute_real', 524, attribute_name )
          nc_stat = NF90_GET_ATT( id, id_var, TRIM( attribute_name ), value )
          CALL handle_error( 'get_attribute_real', 524, attribute_name )
       ENDIF
#endif
    END SUBROUTINE get_attribute_real

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads global or variable-related attributes of type CHARACTER
!> Remark: reading attributes of type NF_STRING return an error code 56 -
!> Attempt to convert between text & numbers.
!------------------------------------------------------------------------------!
     SUBROUTINE get_attribute_string( id, attribute_name, value, global,       &
                                      variable_name, no_abort )

       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=*)                ::  attribute_name   !< attribute name
       CHARACTER(LEN=*), OPTIONAL      ::  variable_name    !< variable name
       CHARACTER(LEN=*), INTENT(INOUT) ::  value            !< read value

       INTEGER(iwp), INTENT(IN)    ::  id               !< file id
       INTEGER(iwp)                ::  id_var           !< variable id

       LOGICAL ::  check_error                          !< flag indicating if handle_error shall be checked
       LOGICAL, INTENT(IN) ::  global                   !< flag indicating global attribute
       LOGICAL, INTENT(IN), OPTIONAL ::  no_abort       !< flag indicating if errors should be checked
#if defined( __netcdf )

       IF ( PRESENT( no_abort ) )  THEN
          check_error = no_abort
       ELSE
          check_error = .TRUE.
       ENDIF
!
!--    Read global attribute
       IF ( global )  THEN
          nc_stat = NF90_GET_ATT( id, NF90_GLOBAL, TRIM( attribute_name ), value )
          IF ( check_error)  CALL handle_error( 'get_attribute_string global', 525, attribute_name )
!
!--    Read attributes referring to a single variable. Therefore, first inquire
!--    variable id
       ELSE
          nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
          IF ( check_error)  CALL handle_error( 'get_attribute_string', 525, attribute_name )

          nc_stat = NF90_GET_ATT( id, id_var, TRIM( attribute_name ), value )
          IF ( check_error)  CALL handle_error( 'get_attribute_string',525, attribute_name )

       ENDIF
#endif
    END SUBROUTINE get_attribute_string



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Get dimension array for a given dimension
!------------------------------------------------------------------------------!
     SUBROUTINE get_dimension_length( id, dim_len, variable_name )
       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=*)            ::  variable_name    !< dimension name
       CHARACTER(LEN=100)          ::  dum              !< dummy variable to receive return character

       INTEGER(iwp)                ::  dim_len          !< dimension size
       INTEGER(iwp), INTENT(IN)    ::  id               !< file id
       INTEGER(iwp)                ::  id_dim           !< dimension id

#if defined( __netcdf )
!
!--    First, inquire dimension ID
       nc_stat = NF90_INQ_DIMID( id, TRIM( variable_name ), id_dim )
       CALL handle_error( 'get_dimension_length', 526, variable_name )
!
!--    Inquire dimension length
       nc_stat = NF90_INQUIRE_DIMENSION( id, id_dim, dum, LEN = dim_len )
       CALL handle_error( 'get_dimension_length', 526, variable_name )

#endif
    END SUBROUTINE get_dimension_length

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Routine for reading-in a character string from the chem emissions netcdf 
!> input file.  
!------------------------------------------------------------------------------!

    SUBROUTINE get_variable_string( id, variable_name, var_string, names_number )

       USE indices
       USE pegrid

       IMPLICIT NONE

       CHARACTER (LEN=25), ALLOCATABLE, DIMENSION(:), INTENT(INOUT)  :: var_string

       CHARACTER(LEN=*)                                              :: variable_name          !> variable name 

       CHARACTER (LEN=1), ALLOCATABLE, DIMENSION(:,:)                :: tmp_var_string         !> variable to be read


       INTEGER(iwp), INTENT(IN)                                      :: id                     !> file id

       INTEGER(iwp), INTENT(IN)                                      :: names_number           !> number of names

       INTEGER(iwp)                                                  :: id_var                 !> variable id

       INTEGER(iwp)                                                  :: i,j                    !> index to go through the length of the dimensions

       INTEGER(iwp)                                                  :: max_string_length=25   !> this is both the maximum length of a name, but also  
                                                                                            ! the number of the components of the first dimensions
                                                                                            ! (rows)
#if defined( __netcdf )

       ALLOCATE(tmp_var_string(max_string_length,names_number))

       ALLOCATE(var_string(names_number))

    !-- Inquire variable id
       nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )


    !-- Get variable
    !-- Start cycle over the emission species
       DO i = 1, names_number 
       !-- read the first letter of each component
          nc_stat = NF90_GET_VAR( id, id_var, var_string(i), start = (/ 1,i /), &
                                 count = (/ 1,1 /) )
          CALL handle_error( 'get_variable_string', 701 )

       !-- Start cycle over charachters
          DO j = 1, max_string_length
                       
          !-- read the rest of the components of the name
             nc_stat = NF90_GET_VAR( id, id_var, tmp_var_string(j,i), start = (/ j,i /),&
                                     count = (/ 1,1 /) )
             CALL handle_error( 'get_variable_string', 702 )

             IF ( iachar(tmp_var_string(j,i) ) == 0 ) THEN
                  tmp_var_string(j,i)=''
             ENDIF

             IF ( j>1 ) THEN
             !-- Concatenate first letter of the name and the others
                var_string(i)=TRIM(var_string(i)) // TRIM(tmp_var_string(j,i))

             ENDIF 
          ENDDO 
       ENDDO

#endif
    END SUBROUTINE get_variable_string


!
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Generalized routine for reading strings from a netCDF variable
!> to replace existing get_variable_string ( )
!>
!> Improvements:
!>   - Expanded string size limit from 25 to 512
!>   - No longer removes spaces between text magically (this seems to have
!>     been aimed at a very specific application, but I don't know what)
!>   - Simplified implementation
!>
!> Caveats:
!>   - Somehow I could not get the subroutine to work with str_array(:,:)
!>     so I reverted to a hard-coded str_array(:,512), hopefully large enough
!>     for most general applications.  This also means the character variable
!>     used for str_array must be of size (:,512)
!>     (ecc 20200128)    
!------------------------------------------------------------------------------!

 SUBROUTINE get_variable_string_generic ( id, var_name, str_array, num_str, str_len )

    IMPLICIT NONE

    CHARACTER(LEN=*),                INTENT(IN)    :: var_name       !> netCDF variable name
    CHARACTER(LEN=512), ALLOCATABLE, INTENT(INOUT) :: str_array(:)   !> output string array

    INTEGER(iwp)              :: buf_len   !> string buffer size
    INTEGER(iwp)              :: id_var    !> netCDF variable ID
    INTEGER(iwp)              :: k         !> generic counter

    INTEGER(iwp), INTENT(IN)  :: id        !> netCDF file ID
    INTEGER(iwp), INTENT(IN)  :: num_str   !> number of string elements in array
    INTEGER(iwp), INTENT(IN)  :: str_len   !> size of each string element

#if defined( __netcdf )

!
!-- set buffer length to up to hard-coded string size

    buf_len = MIN( ABS(str_len), 512 )

!
!-- allocate necessary memories for string variables

    ALLOCATE(str_array(num_str))
!
!-- get variable id

    nc_stat = NF90_INQ_VARID( id, TRIM(var_name), id_var )
!
!-- extract string variables

    DO k = 1, num_str
       str_array(k) = ''
       nc_stat = NF90_GET_VAR( id, id_var, str_array(k),  &
                      start = (/ 1, k /), count = (/ buf_len, 1 /)  )
       CALL handle_error ( 'get_variable_string_generic', 702 )
    ENDDO

#endif

 END SUBROUTINE get_variable_string_generic


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads a character variable in a 1D array
!------------------------------------------------------------------------------!
     SUBROUTINE get_variable_1d_char( id, variable_name, var )

       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=*)            ::  variable_name          !< variable name
       CHARACTER(LEN=*), DIMENSION(:), INTENT(INOUT) ::  var  !< variable to be read

       INTEGER(iwp)                ::  i                !< running index over variable dimension
       INTEGER(iwp), INTENT(IN)    ::  id               !< file id
       INTEGER(iwp)                ::  id_var           !< dimension id
       
       INTEGER(iwp), DIMENSION(2)  ::  dimid            !< dimension IDs
       INTEGER(iwp), DIMENSION(2)  ::  dimsize          !< dimension size

#if defined( __netcdf )

!
!--    First, inquire variable ID
       nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
       CALL handle_error( 'get_variable_1d_int', 527, variable_name )
!
!--    Inquire dimension IDs
       nc_stat = NF90_INQUIRE_VARIABLE( id, id_var, dimids = dimid(1:2) )
       CALL handle_error( 'get_variable_1d_char', 527, variable_name )
!
!--    Read dimesnion length
       nc_stat = NF90_INQUIRE_DIMENSION( id, dimid(1), LEN = dimsize(1) )
       nc_stat = NF90_INQUIRE_DIMENSION( id, dimid(2), LEN = dimsize(2) )
       
!
!--    Read character array. Note, each element is read individually, in order
!--    to better separate single strings. 
       DO  i = 1, dimsize(2)
          nc_stat = NF90_GET_VAR( id, id_var, var(i),                          &
                                  start = (/ 1, i /),                          &
                                  count = (/ dimsize(1), 1 /) )
          CALL handle_error( 'get_variable_1d_char', 527, variable_name )
       ENDDO     
                          
#endif
    END SUBROUTINE get_variable_1d_char

    
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads a 1D integer variable from file.
!------------------------------------------------------------------------------!
     SUBROUTINE get_variable_1d_int( id, variable_name, var )

       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=*)            ::  variable_name    !< variable name

       INTEGER(iwp), INTENT(IN)    ::  id               !< file id
       INTEGER(iwp)                ::  id_var           !< dimension id

       INTEGER(iwp), DIMENSION(:), INTENT(INOUT) ::  var  !< variable to be read
#if defined( __netcdf )

!
!--    First, inquire variable ID
       nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
       CALL handle_error( 'get_variable_1d_int', 527, variable_name )
!
!--    Inquire dimension length
       nc_stat = NF90_GET_VAR( id, id_var, var )
       CALL handle_error( 'get_variable_1d_int', 527, variable_name )

#endif
    END SUBROUTINE get_variable_1d_int

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads a 1D float variable from file.
!------------------------------------------------------------------------------!
     SUBROUTINE get_variable_1d_real( id, variable_name, var )

       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=*)            ::  variable_name    !< variable name

       INTEGER(iwp), INTENT(IN)    ::  id               !< file id
       INTEGER(iwp)                ::  id_var           !< dimension id

       REAL(wp), DIMENSION(:), INTENT(INOUT) ::  var    !< variable to be read
#if defined( __netcdf )

!
!--    First, inquire variable ID
       nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
       CALL handle_error( 'get_variable_1d_real', 528, variable_name )
!
!--    Inquire dimension length
       nc_stat = NF90_GET_VAR( id, id_var, var )
       CALL handle_error( 'get_variable_1d_real', 528, variable_name )

#endif
    END SUBROUTINE get_variable_1d_real


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads a time-dependent 1D float variable from file.
!------------------------------------------------------------------------------!
    SUBROUTINE get_variable_pr( id, variable_name, t, var )

       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=*)                      ::  variable_name    !< variable name

       INTEGER(iwp), INTENT(IN)              ::  id               !< file id
       INTEGER(iwp), DIMENSION(1:2)          ::  id_dim           !< dimension ids
       INTEGER(iwp)                          ::  id_var           !< dimension id
       INTEGER(iwp)                          ::  n_file           !< number of data-points in file along z dimension
       INTEGER(iwp), INTENT(IN)              ::  t                !< timestep number

       REAL(wp), DIMENSION(:), INTENT(INOUT) ::  var  !< variable to be read

#if defined( __netcdf )
!
!--    First, inquire variable ID
       nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
!
!--    Inquire dimension size of vertical dimension
       nc_stat = NF90_INQUIRE_VARIABLE( id, id_var, DIMIDS = id_dim )
       nc_stat = NF90_INQUIRE_DIMENSION( id, id_dim(1), LEN = n_file )
!
!--    Read variable.
       nc_stat = NF90_GET_VAR( id, id_var, var,                                &
                               start = (/ 1,      t     /),                    &
                               count = (/ n_file, 1     /) )
       CALL handle_error( 'get_variable_pr', 529, variable_name )

#endif
    END SUBROUTINE get_variable_pr


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads a per-surface pars variable from file. Because all surfaces are stored
!> as flat 1-D array, each PE has to scan the data and find the surface indices
!> belonging to its subdomain. During this scan, it also builds a necessary
!> (j,i) index.
!------------------------------------------------------------------------------!
    SUBROUTINE get_variable_surf( id, variable_name, surf )

       USE pegrid

       USE indices,                                            &
           ONLY:  nxl, nxr, nys, nyn

       USE control_parameters,                                 &
           ONLY: dz, message_string

       USE grid_variables,                                     &
           ONLY: dx, dy

       USE basic_constants_and_equations_mod,                  &
           ONLY: pi

       IMPLICIT NONE

       INTEGER(iwp), PARAMETER                   ::  nsurf_pars_read = 2**15 !< read buffer size (value > 10^15 makes problem with ifort)

       CHARACTER(LEN=*)                          ::  variable_name !< variable name

       INTEGER(iwp), DIMENSION(6)                ::  coords        !< integer coordinates of surface
       INTEGER(iwp)                              ::  i, j
       INTEGER(iwp)                              ::  isurf         !< netcdf surface index
       INTEGER(iwp)                              ::  is            !< local surface index
       INTEGER(iwp), INTENT(IN)                  ::  id            !< file id
       INTEGER(iwp), DIMENSION(2)                ::  id_dim        !< dimension ids
       INTEGER(iwp)                              ::  id_var        !< variable id
       INTEGER(iwp)                              ::  id_zs         !< zs variable id
       INTEGER(iwp)                              ::  id_ys         !< ys variable id
       INTEGER(iwp)                              ::  id_xs         !< xs variable id
       INTEGER(iwp)                              ::  id_zenith     !< zeith variable id
       INTEGER(iwp)                              ::  id_azimuth    !< azimuth variable id
       INTEGER(iwp)                              ::  is0, isc      !< read surface start and count
       INTEGER(iwp)                              ::  nsurf         !< total number of surfaces in file
       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  nsurf_ji      !< numbers of surfaces by coords

       TYPE(pars_surf)                           ::  surf          !< parameters variable to be loaded
       REAL(wp), DIMENSION(:,:), ALLOCATABLE     ::  pars_read     !< read buffer
       REAL(wp), DIMENSION(:), ALLOCATABLE       ::  zs, ys, xs    !< read buffer for zs(s), ys, xs
       REAL(wp), DIMENSION(:), ALLOCATABLE       ::  zenith        !< read buffer for zenith(s)
       REAL(wp), DIMENSION(:), ALLOCATABLE       ::  azimuth       !< read buffer for azimuth(s)
       REAL(wp)                                  ::  oro_max_l     !< maximum terrain height under building

#if defined( __netcdf )
!
!--    First, inquire variable ID
       nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
       nc_stat = NF90_INQ_VARID( id, 'zs',                  id_zs )
       nc_stat = NF90_INQ_VARID( id, 'ys',                  id_ys )
       nc_stat = NF90_INQ_VARID( id, 'xs',                  id_xs )
       nc_stat = NF90_INQ_VARID( id, 'zenith',              id_zenith )
       nc_stat = NF90_INQ_VARID( id, 'azimuth',             id_azimuth )
!
!--    Inquire dimension sizes
       nc_stat = NF90_INQUIRE_VARIABLE( id, id_var, DIMIDS = id_dim )
       nc_stat = NF90_INQUIRE_DIMENSION( id, id_dim(1), LEN = nsurf )
       nc_stat = NF90_INQUIRE_DIMENSION( id, id_dim(2), LEN = surf%np )

       ALLOCATE ( pars_read( nsurf_pars_read, surf%np ),        &
                  zs(nsurf_pars_read), ys(nsurf_pars_read),     &
                  xs(nsurf_pars_read), zenith(nsurf_pars_read), &
                  azimuth(nsurf_pars_read),                     &
                  nsurf_ji(nys:nyn, nxl:nxr) )

       nsurf_ji(:,:) = 0
!
!--    Scan surface coordinates, count local
       is0 = 1
       DO
          isc = MIN(nsurf_pars_read, nsurf - is0 + 1)
          IF ( isc <= 0 )  EXIT

          nc_stat = NF90_GET_VAR( id, id_ys, ys,     &
                                  start = (/ is0 /), &
                                  count = (/ isc /) )
          nc_stat = NF90_GET_VAR( id, id_xs, xs,     &
                                  start = (/ is0 /), &
                                  count = (/ isc /) )
          nc_stat = NF90_GET_VAR( id, id_zenith, zenith,      &
                                  start = (/ is0 /), &
                                  count = (/ isc /) )
          nc_stat = NF90_GET_VAR( id, id_azimuth, azimuth,    &
                                  start = (/ is0 /), &
                                  count = (/ isc /) )
          CALL handle_error( 'get_variable_surf', 682, 'azimuth' )
          
          DO  isurf = 1, isc
!
!--          Parse coordinates, detect if belongs to subdomain
             coords = transform_coords( xs(isurf), ys(isurf),         &
                                        zenith(isurf), azimuth(isurf) )
             IF ( coords(2) < nys  .OR.  coords(2) > nyn  .OR.  &
                  coords(3) < nxl  .OR.  coords(3) > nxr )  CYCLE

             nsurf_ji(coords(2), coords(3)) = nsurf_ji(coords(2), coords(3)) + 1
          ENDDO

          is0 = is0 + isc
       ENDDO
!
!--    Populate reverse index from surface counts
       ALLOCATE ( surf%index_ji( 2, nys:nyn, nxl:nxr ) )

       isurf = 1
       DO  j = nys, nyn
          DO  i = nxl, nxr
             surf%index_ji(:,j,i) = (/ isurf, isurf + nsurf_ji(j,i) - 1 /)
             isurf = isurf + nsurf_ji(j,i)
          ENDDO
       ENDDO

       surf%nsurf = isurf - 1
       ALLOCATE( surf%pars( 0:surf%np-1, surf%nsurf ), &
                 surf%coords( 6, surf%nsurf ) )
!
!--    Scan surfaces again, saving pars into allocated structures
       nsurf_ji(:,:) = 0
       is0 = 1
       DO
          isc = MIN(nsurf_pars_read, nsurf - is0 + 1)
          IF ( isc <= 0 )  EXIT

          nc_stat = NF90_GET_VAR( id, id_var, pars_read(1:isc, 1:surf%np), &
                                  start = (/ is0, 1       /),              &
                                  count = (/ isc, surf%np /) )
          CALL handle_error( 'get_variable_surf', 683, variable_name )

          nc_stat = NF90_GET_VAR( id, id_zs, zs,                           &
                                  start = (/ is0 /),                       &
                                  count = (/ isc /) )
          nc_stat = NF90_GET_VAR( id, id_ys, ys,                           &
                                  start = (/ is0 /),                       &
                                  count = (/ isc /) )
          nc_stat = NF90_GET_VAR( id, id_xs, xs,                           &
                                  start = (/ is0 /),                       &
                                  count = (/ isc /) )
          nc_stat = NF90_GET_VAR( id, id_zenith, zenith,                   &
                                  start = (/ is0 /),                       &
                                  count = (/ isc /) )
          nc_stat = NF90_GET_VAR( id, id_azimuth, azimuth,                 &
                                  start = (/ is0 /),                       &
                                  count = (/ isc /) )
          
          DO  isurf = 1, isc
!
!--          Parse coordinates, detect if belongs to subdomain
             coords = transform_coords( xs(isurf), ys(isurf),         &
                                        zenith(isurf), azimuth(isurf) )
             IF ( coords(2) < nys  .OR.  coords(2) > nyn  .OR.  &
                  coords(3) < nxl  .OR.  coords(3) > nxr )  CYCLE
!
!--          Determine maximum terrain under building (base z-coordinate). Using
!--          normal vector to locate building inner coordinates.
             oro_max_l = buildings_f%oro_max(coords(2)-coords(5), coords(3)-coords(6))
             IF  ( oro_max_l == buildings_f%fill1 )  THEN
                WRITE( message_string, * ) 'Found building surface on '   // &
                   'non-building coordinates (xs, ys, zenith, azimuth): ',   &
                   xs(isurf), ys(isurf), zenith(isurf), azimuth(isurf)
                CALL message( 'get_variable_surf', 'PA0684', 2, 2, myid, 6, 0 ) 
             ENDIF
!
!--          Urban layer has no stretching, therefore using dz(1) instead of linear
!--          searching through zu/zw
             coords(1) = NINT((zs(isurf) + oro_max_l) / dz(1) +     &
                              0.5_wp + 0.5_wp * coords(4), KIND=iwp)
!
!--          Save surface entry
             is = surf%index_ji(1, coords(2), coords(3)) + nsurf_ji(coords(2), coords(3))
             surf%pars(:,is) = pars_read(isurf,:)
             surf%coords(:,is) = coords(:)

             nsurf_ji(coords(2), coords(3)) = nsurf_ji(coords(2), coords(3)) + 1
          ENDDO

          is0 = is0 + isc
       ENDDO

       DEALLOCATE( pars_read, zs, ys, xs, zenith, azimuth, nsurf_ji )

    CONTAINS

       PURE FUNCTION transform_coords( x, y, zenith, azimuth )

          REAL(wp), INTENT(in)       ::  x, y    !< surface centre coordinates in metres from origin
          REAL(wp), INTENT(in)       ::  zenith  !< surface normal zenith angle in degrees
          REAL(wp), INTENT(in)       ::  azimuth !< surface normal azimuth angle in degrees

          INTEGER(iwp), DIMENSION(6) ::  transform_coords !< (k,j,i,norm_z,norm_y,norm_x)

          transform_coords(4) = NINT(COS(zenith*pi/180._wp), KIND=iwp)
          IF ( transform_coords(4) == 0 )  THEN
             transform_coords(5) = NINT(COS(azimuth*pi/180._wp), KIND=iwp)
             transform_coords(6) = NINT(SIN(azimuth*pi/180._wp), KIND=iwp)
          ELSE
             transform_coords(5) = 0._wp
             transform_coords(6) = 0._wp
          ENDIF

          transform_coords(1) = -999._wp ! not calculated here
          transform_coords(2) = NINT(y/dy - 0.5_wp + 0.5_wp*transform_coords(5), KIND=iwp)
          transform_coords(3) = NINT(x/dx - 0.5_wp + 0.5_wp*transform_coords(6), KIND=iwp)

       END FUNCTION transform_coords

#endif
    END SUBROUTINE get_variable_surf


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads a 2D REAL variable from a file. Reading is done processor-wise,
!> i.e. each core reads its own domain in slices along x.
!------------------------------------------------------------------------------!
    SUBROUTINE get_variable_2d_real( id, variable_name, var, is, ie, js, je )

       USE indices
       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=*)              ::  variable_name   !< variable name

       INTEGER(iwp)                  ::  i               !< running index along x direction
       INTEGER(iwp)                  ::  ie              !< start index for subdomain input along x direction
       INTEGER(iwp)                  ::  is              !< end index for subdomain input along x direction
       INTEGER(iwp), INTENT(IN)      ::  id              !< file id
       INTEGER(iwp)                  ::  id_var          !< variable id
       INTEGER(iwp)                  ::  j               !< running index along y direction
       INTEGER(iwp)                  ::  je              !< start index for subdomain input along y direction
       INTEGER(iwp)                  ::  js              !< end index for subdomain input along y direction
       
       REAL(wp), DIMENSION(:,:), ALLOCATABLE   ::  tmp   !< temporary variable to read data from file according
                                                         !< to its reverse memory access
       REAL(wp), DIMENSION(:,:), INTENT(INOUT) ::  var   !< variable to be read
#if defined( __netcdf )
!
!--    Inquire variable id
       nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
!
!--    Check for collective read-operation and set respective NetCDF flags if
!--    required. 
       IF ( collective_read )  THEN
#if defined( __netcdf4_parallel )
          nc_stat = NF90_VAR_PAR_ACCESS (id, id_var, NF90_COLLECTIVE)
#endif
       ENDIF

!
!-- Allocate temporary variable according to memory access on file. 
       ALLOCATE( tmp(is:ie,js:je) )
!
!-- Get variable
       nc_stat = NF90_GET_VAR( id, id_var, tmp,            &
                      start = (/ is+1,      js+1 /),       &
                      count = (/ ie-is + 1, je-js+1 /) )   
          CALL handle_error( 'get_variable_2d_real', 530, variable_name )
!
!-- Resort data. Please note, dimension subscripts of var all start at 1. 
          DO  i = is, ie 
             DO  j = js, je 
                var(j-js+1,i-is+1) = tmp(i,j)
             ENDDO
          ENDDO
       
          DEALLOCATE( tmp )

#endif
    END SUBROUTINE get_variable_2d_real

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads a 2D 32-bit INTEGER variable from file. Reading is done processor-wise,
!> i.e. each core reads its own domain in slices along x.
!------------------------------------------------------------------------------!
    SUBROUTINE get_variable_2d_int32( id, variable_name, var, is, ie, js, je )

       USE indices
       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=*)              ::  variable_name   !< variable name

       INTEGER(iwp)                  ::  i               !< running index along x direction
       INTEGER(iwp)                  ::  ie              !< start index for subdomain input along x direction
       INTEGER(iwp)                  ::  is              !< end index for subdomain input along x direction
       INTEGER(iwp), INTENT(IN)      ::  id              !< file id
       INTEGER(iwp)                  ::  id_var          !< variable id
       INTEGER(iwp)                  ::  j               !< running index along y direction
       INTEGER(iwp)                  ::  je              !< start index for subdomain input along y direction
       INTEGER(iwp)                  ::  js              !< end index for subdomain input along y direction
       
       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE   ::  tmp  !< temporary variable to read data from file according
                                                            !< to its reverse memory access
       INTEGER(iwp), DIMENSION(:,:), INTENT(INOUT) ::  var  !< variable to be read
#if defined( __netcdf )
!
!--    Inquire variable id
       nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
!
!--    Check for collective read-operation and set respective NetCDF flags if
!--    required. 
       IF ( collective_read )  THEN
#if defined( __netcdf4_parallel )       
          nc_stat = NF90_VAR_PAR_ACCESS (id, id_var, NF90_COLLECTIVE)
#endif
       ENDIF
!
!--    Allocate temporary variable according to memory access on file. 
       ALLOCATE( tmp(is:ie,js:je) )
!
!--    Get variable
       nc_stat = NF90_GET_VAR( id, id_var, tmp,                                &
                               start = (/ is+1,      js+1 /),                  &
                               count = (/ ie-is + 1, je-js+1 /) )    
                               
       CALL handle_error( 'get_variable_2d_int32', 531, variable_name )                             
!
!--    Resort data. Please note, dimension subscripts of var all start at 1. 
       DO  i = is, ie 
          DO  j = js, je 
             var(j-js+1,i-is+1) = tmp(i,j)
          ENDDO
       ENDDO
       
       DEALLOCATE( tmp )

#endif
    END SUBROUTINE get_variable_2d_int32

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads a 2D 8-bit INTEGER variable from file. Reading is done processor-wise,
!> i.e. each core reads its own domain in slices along x.
!------------------------------------------------------------------------------!
    SUBROUTINE get_variable_2d_int8( id, variable_name, var, is, ie, js, je )

       USE indices
       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=*)              ::  variable_name   !< variable name

       INTEGER(iwp)                  ::  i               !< running index along x direction
       INTEGER(iwp)                  ::  ie              !< start index for subdomain input along x direction
       INTEGER(iwp)                  ::  is              !< end index for subdomain input along x direction
       INTEGER(iwp), INTENT(IN)      ::  id              !< file id
       INTEGER(iwp)                  ::  id_var          !< variable id
       INTEGER(iwp)                  ::  j               !< running index along y direction
       INTEGER(iwp)                  ::  je              !< start index for subdomain input along y direction
       INTEGER(iwp)                  ::  js              !< end index for subdomain input along y direction
       
       INTEGER(KIND=1), DIMENSION(:,:), ALLOCATABLE   ::  tmp  !< temporary variable to read data from file according
                                                               !< to its reverse memory access
       INTEGER(KIND=1), DIMENSION(:,:), INTENT(INOUT) ::  var  !< variable to be read
#if defined( __netcdf )
!
!--    Inquire variable id
       nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
!
!--    Check for collective read-operation and set respective NetCDF flags if
!--    required. 
       IF ( collective_read )  THEN
#if defined( __netcdf4_parallel )        
          nc_stat = NF90_VAR_PAR_ACCESS (id, id_var, NF90_COLLECTIVE)
#endif          
       ENDIF
!
!--    Allocate temporary variable according to memory access on file. 
       ALLOCATE( tmp(is:ie,js:je) )
!
!--    Get variable
       nc_stat = NF90_GET_VAR( id, id_var, tmp,                                &
                               start = (/ is+1,      js+1 /),                  &
                               count = (/ ie-is + 1, je-js+1 /) )   
                               
       CALL handle_error( 'get_variable_2d_int8', 532, variable_name )
!
!--    Resort data. Please note, dimension subscripts of var all start at 1. 
       DO  i = is, ie 
          DO  j = js, je 
             var(j-js+1,i-is+1) = tmp(i,j)
          ENDDO
       ENDDO
       
       DEALLOCATE( tmp )

#endif
    END SUBROUTINE get_variable_2d_int8


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads a 3D 8-bit INTEGER variable from file.
!------------------------------------------------------------------------------!
    SUBROUTINE get_variable_3d_int8( id, variable_name, var, is, ie, js, je,   &
                                     ks, ke )

       USE indices
       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=*)              ::  variable_name   !< variable name

       INTEGER(iwp)                  ::  i               !< index along x direction
       INTEGER(iwp)                  ::  ie              !< start index for subdomain input along x direction
       INTEGER(iwp)                  ::  is              !< end index for subdomain input along x direction
       INTEGER(iwp), INTENT(IN)      ::  id              !< file id
       INTEGER(iwp)                  ::  id_var          !< variable id
       INTEGER(iwp)                  ::  j               !< index along y direction
       INTEGER(iwp)                  ::  je              !< start index for subdomain input along y direction
       INTEGER(iwp)                  ::  js              !< end index for subdomain input along y direction
       INTEGER(iwp)                  ::  k               !< index along any 3rd dimension
       INTEGER(iwp)                  ::  ke              !< start index of 3rd dimension
       INTEGER(iwp)                  ::  ks              !< end index of 3rd dimension

       INTEGER(KIND=1), DIMENSION(:,:,:), ALLOCATABLE   ::  tmp  !< temporary variable to read data from file according 
                                                                 !< to its reverse memory access

       INTEGER(KIND=1), DIMENSION(:,:,:), INTENT(INOUT) ::  var  !< variable to be read
#if defined( __netcdf )

!
!--    Inquire variable id
       nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
!
!--    Check for collective read-operation and set respective NetCDF flags if
!--    required. 
       IF ( collective_read )  THEN
#if defined( __netcdf4_parallel )
          nc_stat = NF90_VAR_PAR_ACCESS (id, id_var, NF90_COLLECTIVE)
#endif
       ENDIF
!
!--    Allocate temporary variable according to memory access on file. 
       ALLOCATE( tmp(is:ie,js:je,ks:ke) )
!
!--    Get variable
       nc_stat = NF90_GET_VAR( id, id_var, tmp,                                &
                               start = (/ is+1,    js+1,    ks+1 /),           &
                               count = (/ ie-is+1, je-js+1, ke-ks+1 /) )

       CALL handle_error( 'get_variable_3d_int8', 533, variable_name )
!
!--    Resort data. Please note, dimension subscripts of var all start at 1. 
       DO  i = is, ie 
          DO  j = js, je
             DO  k = ks, ke
                var(k-ks+1,j-js+1,i-is+1) = tmp(i,j,k)
             ENDDO
          ENDDO
       ENDDO

       DEALLOCATE( tmp )

#endif
    END SUBROUTINE get_variable_3d_int8


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads a 3D float variable from file.
!------------------------------------------------------------------------------!
    SUBROUTINE get_variable_3d_real( id, variable_name, var, is, ie, js, je,   &
                                     ks, ke )

       USE indices
       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=*)              ::  variable_name   !< variable name

       INTEGER(iwp)                  ::  i               !< index along x direction
       INTEGER(iwp)                  ::  ie              !< start index for subdomain input along x direction
       INTEGER(iwp)                  ::  is              !< end index for subdomain input along x direction
       INTEGER(iwp), INTENT(IN)      ::  id              !< file id
       INTEGER(iwp)                  ::  id_var          !< variable id
       INTEGER(iwp)                  ::  j               !< index along y direction
       INTEGER(iwp)                  ::  je              !< start index for subdomain input along y direction
       INTEGER(iwp)                  ::  js              !< end index for subdomain input along y direction
       INTEGER(iwp)                  ::  k               !< index along any 3rd dimension
       INTEGER(iwp)                  ::  ke              !< start index of 3rd dimension
       INTEGER(iwp)                  ::  ks              !< end index of 3rd dimension

       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE   ::  tmp !< temporary variable to read data from file according 
                                                         !< to its reverse memory access

       REAL(wp), DIMENSION(:,:,:), INTENT(INOUT) ::  var !< variable to be read
#if defined( __netcdf )

!
!--    Inquire variable id
       nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )  
!
!--    Check for collective read-operation and set respective NetCDF flags if
!--    required. 
       IF ( collective_read )  THEN
#if defined( __netcdf4_parallel )
          nc_stat = NF90_VAR_PAR_ACCESS (id, id_var, NF90_COLLECTIVE)
#endif
       ENDIF
!
!--    Allocate temporary variable according to memory access on file. 
       ALLOCATE( tmp(is:ie,js:je,ks:ke) )
!
!--    Get variable
       nc_stat = NF90_GET_VAR( id, id_var, tmp,                                &
                               start = (/ is+1,    js+1,    ks+1 /),           &
                               count = (/ ie-is+1, je-js+1, ke-ks+1 /) )   

       CALL handle_error( 'get_variable_3d_real', 534, variable_name )
!
!--    Resort data. Please note, dimension subscripts of var all start at 1. 
       DO  i = is, ie 
          DO  j = js, je
             DO  k = ks, ke
                var(k-ks+1,j-js+1,i-is+1) = tmp(i,j,k)
             ENDDO
          ENDDO
       ENDDO

       DEALLOCATE( tmp )

#endif
    END SUBROUTINE get_variable_3d_real

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads a 4D float variable from file. 
!------------------------------------------------------------------------------!
    SUBROUTINE get_variable_4d_real( id, variable_name, var, is, ie, js, je,   &
                                     k1s, k1e, k2s, k2e )

       USE indices
       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=*)              ::  variable_name   !< variable name

       INTEGER(iwp)                  ::  i               !< index along x direction
       INTEGER(iwp)                  ::  ie              !< start index for subdomain input along x direction
       INTEGER(iwp)                  ::  is              !< end index for subdomain input along x direction
       INTEGER(iwp), INTENT(IN)      ::  id              !< file id
       INTEGER(iwp)                  ::  id_var          !< variable id
       INTEGER(iwp)                  ::  j               !< index along y direction
       INTEGER(iwp)                  ::  je              !< start index for subdomain input along y direction
       INTEGER(iwp)                  ::  js              !< end index for subdomain input along y direction
       INTEGER(iwp)                  ::  k1              !< index along 3rd direction
       INTEGER(iwp)                  ::  k1e             !< start index for 3rd dimension
       INTEGER(iwp)                  ::  k1s             !< end index for 3rd dimension
       INTEGER(iwp)                  ::  k2              !< index along 4th direction
       INTEGER(iwp)                  ::  k2e             !< start index for 4th dimension
       INTEGER(iwp)                  ::  k2s             !< end index for 4th dimension

       REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE   ::  tmp  !< temporary variable to read data from file according
                                                            !< to its reverse memory access
       REAL(wp), DIMENSION(:,:,:,:), INTENT(INOUT) ::  var  !< variable to be read
#if defined( __netcdf )

!
!--    Inquire variable id
       nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
!
!--    Check for collective read-operation and set respective NetCDF flags if
!--    required. 
       IF ( collective_read )  THEN
#if defined( __netcdf4_parallel )       
          nc_stat = NF90_VAR_PAR_ACCESS (id, id_var, NF90_COLLECTIVE)
#endif
       ENDIF

!
!-- Allocate temporary variable according to memory access on file. 
       ALLOCATE( tmp(is:ie,js:je,k1s:k1e,k2s:k2e) )
!
!-- Get variable
       nc_stat = NF90_GET_VAR( id, id_var, tmp,                                &
                      start = (/ is+1,    js+1,    k1s+1,     k2s+1 /),        &
                      count = (/ ie-is+1, je-js+1, k1e-k1s+1, k2e-k2s+1 /) )

          CALL handle_error( 'get_variable_4d_real', 535, variable_name )
!
!-- Resort data. Please note, dimension subscripts of var all start at 1. 
       DO  i = is, ie 
          DO  j = js, je
             DO  k1 = k1s, k1e
                DO  k2 = k2s, k2e
                   var(k2-k2s+1,k1-k1s+1,j-js+1,i-is+1) = tmp(i,j,k1,k2)
                ENDDO
             ENDDO
          ENDDO
       ENDDO

       DEALLOCATE( tmp )

#endif

    END SUBROUTINE get_variable_4d_real

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads a 4D float variable from file and store it to a 3-d variable.
!------------------------------------------------------------------------------!
    SUBROUTINE get_variable_4d_to_3d_real( id, variable_name, var, ns, is, ie, js, je,   &
                                           ks, ke )

       USE indices
       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=*)              ::  variable_name   !< variable name

       INTEGER(iwp)                  ::  i               !< index along x direction
       INTEGER(iwp)                  ::  ie              !< end index for subdomain input along x direction
       INTEGER(iwp)                  ::  is              !< start index for subdomain input along x direction
       INTEGER(iwp), INTENT(IN)      ::  id              !< file id
       INTEGER(iwp)                  ::  id_var          !< variable id
       INTEGER(iwp)                  ::  j               !< index along y direction
       INTEGER(iwp)                  ::  je              !< end index for subdomain input along y direction
       INTEGER(iwp)                  ::  js              !< start index for subdomain input along y direction
       INTEGER(iwp)                  ::  k               !< index along any 4th dimension
       INTEGER(iwp)                  ::  ke              !< end index of 4th dimension
       INTEGER(iwp)                  ::  ks              !< start index of 4th dimension
       INTEGER(iwp)                  ::  ns              !< start index for subdomain input along n dimension

       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE   ::  tmp !< temporary variable to read data from file according 
                                                         !< to its reverse memory access

       REAL(wp), DIMENSION(:,:,:), INTENT(INOUT) ::  var  !< variable where the read data have to be stored: 
                                                          !< one dimension is reduced in the process
#if defined( __netcdf )

!
!--    Inquire variable id
       nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
!
!--    Check for collective read-operation and set respective NetCDF flags if
!--    required. 
       IF ( collective_read )  THEN
          nc_stat = NF90_VAR_PAR_ACCESS (id, id_var, NF90_COLLECTIVE)
       ENDIF
!
!--    Allocate temporary variable according to memory access on file. 
       ALLOCATE( tmp(is:ie,js:je,ks:ke) )
!
!--    Get variable
       nc_stat = NF90_GET_VAR( id, id_var, tmp,                                &
                               start = (/ is+1,    js+1,    ks+1,   ns+1 /),   &
                               count = (/ ie-is+1, je-js+1, ke-ks+1, 1   /) )

       CALL handle_error( 'get_variable_4d_to_3d_real', 536, variable_name )
!
!--    Resort data. Please note, dimension subscripts of var all start at 1.
       DO  i = is, ie
          DO  j = js, je
             DO  k = ks, ke
                var(k-ks+1,j-js+1,i-is+1) = tmp(i,j,k)
             ENDDO
          ENDDO
       ENDDO

      DEALLOCATE( tmp )

#endif
    END SUBROUTINE get_variable_4d_to_3d_real

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads a 3D float variables from dynamic driver, such as time-dependent xy-, 
!> xz- or yz-boundary data as well as 3D initialization data. Please note, 
!> the passed arguments are start indices and number of elements in each 
!> dimension, which is in contrast to the other 3d versions where start- and 
!> end indices are passed. The different handling of 3D dynamic variables is 
!> due to its asymmetry for the u- and v component.
!------------------------------------------------------------------------------!
    SUBROUTINE get_variable_3d_real_dynamic( id, variable_name, var,           &
                                             i1s, i2s, i3s,                    &
                                             count_1, count_2, count_3,        &
                                             par_access )
                                
       USE indices
       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=*)              ::  variable_name   !< variable name

       LOGICAL                       ::  par_access      !< additional flag indicating whether parallel read operations should be performed or not 
       
       INTEGER(iwp)                  ::  count_1         !< number of elements to be read along 1st dimension (with respect to file)
       INTEGER(iwp)                  ::  count_2         !< number of elements to be read along 2nd dimension (with respect to file)
       INTEGER(iwp)                  ::  count_3         !< number of elements to be read along 3rd dimension (with respect to file)
       INTEGER(iwp)                  ::  i1              !< running index along 1st dimension on file
       INTEGER(iwp)                  ::  i1s             !< start index for subdomain input along 1st dimension (with respect to file)
       INTEGER(iwp)                  ::  i2              !< running index along 2nd dimension on file       
       INTEGER(iwp)                  ::  i2s             !< start index for subdomain input along 2nd dimension (with respect to file)
       INTEGER(iwp)                  ::  i3              !< running index along 3rd dimension on file 
       INTEGER(iwp)                  ::  i3s             !< start index of 3rd dimension, in dynamic file this is either time (2D boundary) or z (3D)
       INTEGER(iwp), INTENT(IN)      ::  id              !< file id
       INTEGER(iwp)                  ::  id_var          !< variable id
       INTEGER(iwp)                  ::  lb1             !< lower bound of 1st dimension (with respect to file)
       INTEGER(iwp)                  ::  lb2             !< lower bound of 2nd dimension (with respect to file)
       INTEGER(iwp)                  ::  lb3             !< lower bound of 3rd dimension (with respect to file)
       INTEGER(iwp)                  ::  ub1             !< upper bound of 1st dimension (with respect to file)
       INTEGER(iwp)                  ::  ub2             !< upper bound of 2nd dimension (with respect to file)
       INTEGER(iwp)                  ::  ub3             !< upper bound of 3rd dimension (with respect to file)

       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE   ::  tmp !< temporary variable to read data from file according
                                                         !< to its reverse memory access
       
       REAL(wp), DIMENSION(:,:,:), INTENT(INOUT) ::  var !< input variable
       
#if defined( __netcdf )
!
!--    Inquire variable id.
       nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
!
!--    Check for collective read-operation and set respective NetCDF flags if
!--    required. 
!--    Please note, in contrast to the other input routines where each PEs 
!--    reads its subdomain data, dynamic input data not by all PEs, only
!--    by those which encompass lateral model boundaries. Hence, collective
!--    read operations are only enabled for top-boundary data.
       IF ( collective_read  .AND.  par_access )  THEN
#if defined( __netcdf4_parallel )       
          nc_stat = NF90_VAR_PAR_ACCESS (id, id_var, NF90_COLLECTIVE)
#endif
       ENDIF    
!
!--    Allocate temporary variable according to memory access on file. 
!--    Therefore, determine dimension bounds of input array. 
       lb1 = LBOUND(var,3)
       ub1 = UBOUND(var,3)
       lb2 = LBOUND(var,2)
       ub2 = UBOUND(var,2)
       lb3 = LBOUND(var,1)
       ub3 = UBOUND(var,1)
       ALLOCATE( tmp(lb1:ub1,lb2:ub2,lb3:ub3) )
!
!--    Get variable
       nc_stat = NF90_GET_VAR( id, id_var, tmp,                                &
                               start = (/ i1s,     i2s,     i3s /),            &
                               count = (/ count_1, count_2, count_3 /) )

       CALL handle_error( 'get_variable_3d_real_dynamic', 537, variable_name )
!
!--    Resort data. Please note, dimension subscripts of var all start at 1. 
       DO  i3 = lb3, ub3
          DO i2 = lb2, ub2
             DO  i1 = lb1, ub1
                var(i3,i2,i1) = tmp(i1,i2,i3)
             ENDDO
          ENDDO
       ENDDO
       
       DEALLOCATE( tmp )       
#endif
    END SUBROUTINE get_variable_3d_real_dynamic

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads a 5D float variable from file and store it to a 4-d variable.
!------------------------------------------------------------------------------!
    SUBROUTINE get_variable_5d_to_4d_real( id, variable_name, var,              &
                                           ns, ts, te, is, ie, js, je, ks, ke )

       USE indices
       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=*)              ::  variable_name   !< variable name

       INTEGER(iwp)                  ::  ns              !< start index for subdomain input along n dimension: 
                                                         !< ns coincides here with ne, since, we select only one 
                                                         !< value along the 1st dimension n

       INTEGER(iwp)                  ::  t               !< index along t direction
       INTEGER(iwp)                  ::  te              !< end index for subdomain input along t direction
       INTEGER(iwp)                  ::  ts              !< start index for subdomain input along t direction

       INTEGER(iwp)                  ::  i               !< index along x direction
       INTEGER(iwp)                  ::  ie              !< end index for subdomain input along x direction
       INTEGER(iwp)                  ::  is              !< start index for subdomain input along x direction
       INTEGER(iwp), INTENT(IN)      ::  id              !< file id
       INTEGER(iwp)                  ::  id_var          !< variable id
       INTEGER(iwp)                  ::  j               !< index along y direction
       INTEGER(iwp)                  ::  je              !< end index for subdomain input along y direction
       INTEGER(iwp)                  ::  js              !< start index for subdomain input along y direction
       INTEGER(iwp)                  ::  k               !< index along any 5th dimension
       INTEGER(iwp)                  ::  ke              !< end index of 5th dimension
       INTEGER(iwp)                  ::  ks              !< start index of 5th dimension

       REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE   ::  tmp !< temporary variable to read data from file according 
                                                           ! to its reverse memory access
       REAL(wp), DIMENSION(:,:,:,:), INTENT(INOUT) ::  var !< variable to be read
#if defined( __netcdf )
!
!--    Inquire variable id
       nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
!
!--    Check for collective read-operation and set respective NetCDF flags if
!--    required.
       IF ( collective_read )  THEN
          nc_stat = NF90_VAR_PAR_ACCESS (id, id_var, NF90_COLLECTIVE)
       ENDIF
!
!--    Allocate temporary variable according to memory access on file. 
       ALLOCATE( tmp(ks:ke,js:je,is:is,ts:te) )
!
!--    Get variable
       nc_stat = NF90_GET_VAR( id, id_var, tmp,                                &
                               start = (/ ks+1, js+1, is+1, ts+1, ns /),       &
                               count = (/ ke-ks+1, je-js+1, ie-is+1, te-ts+1, 1 /) )

       CALL handle_error( 'get_variable_5d_to_4d_real', 538, variable_name )
!
!--    Resort data. Please note, dimension subscripts of var all start at 1.

       DO  t = ts, te 
          DO  i = is, ie 
             DO  j = js, je
                DO  k = ks, ke
                   var(t-ts+1,i-is+1,j-js+1,k-ks+1) = tmp(k,j,i,t)
                ENDDO
             ENDDO
          ENDDO
       ENDDO 

       DEALLOCATE( tmp )
#endif
    END SUBROUTINE get_variable_5d_to_4d_real

    
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads a 5D float variable from file. 
!> NOTE - This subroutine is used specific for reading NC variable
!>        emission_values having a "z" dimension.  Said dimension
!>        is to be removed in the future and this subroutine shall
!>        be depreciated accordingly (ecc 20190418)
!------------------------------------------------------------------------------!
    SUBROUTINE get_variable_5d_real( id, variable_name, var, is, ie, js, je,   &
                                     k1s, k1e, k2s, k2e, k3s, k3e )

       USE indices
       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=*)          ::  variable_name  !< variable name

       INTEGER(iwp)              :: i       !< i index
       INTEGER(iwp)              :: ie      !< i index start
       INTEGER(iwp)              :: is      !< i index end
       INTEGER(iwp)              :: id_var  !< netCDF variable ID (varid)
       INTEGER(iwp)              :: j       !< j index
       INTEGER(iwp)              :: je      !< j index start
       INTEGER(iwp)              :: js      !< j index end
       INTEGER(iwp)              :: k1      !< k1 index
       INTEGER(iwp)              :: k1e     !< k1 index start
       INTEGER(iwp)              :: k1s     !< k1 index end
       INTEGER(iwp)              :: k2      !< k2 index
       INTEGER(iwp)              :: k2e     !< k2 index start
       INTEGER(iwp)              :: k2s     !< k2 index end
       INTEGER(iwp)              :: k3      !< k3 index
       INTEGER(iwp)              :: k3e     !< k3 index start
       INTEGER(iwp)              :: k3s     !< k3 index end
       INTEGER(iwp), INTENT(IN)  :: id      !< netCDF file ID (ncid)

       REAL(wp), DIMENSION(:,:,:,:,:), ALLOCATABLE    :: tmp  !< temp array to read data from file
       REAL(wp), DIMENSION(:,:,:,:,:), INTENT(INOUT)  :: var  !< variable to be read 

#if defined( __netcdf )

!
!-- Inquire variable id

       nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )

!
!-- Check for collective read-operation and set respective NetCDF flags if required.
 
       IF ( collective_read )  THEN

#if defined( __netcdf4_parallel )       
          nc_stat = NF90_VAR_PAR_ACCESS (id, id_var, NF90_COLLECTIVE)
#endif

       ENDIF

!
!-- Allocate temporary variable according to memory access on file. 

       ALLOCATE( tmp(is:ie,js:je,k1s:k1e,k2s:k2e,k3s:k3e) )

!
!-- Get variable from file

       nc_stat = NF90_GET_VAR ( id, id_var, tmp,                                         &
                      start = (/ is+1,    js+1,    k1s+1,     k2s+1,     k3s+1 /),       &
                      count = (/ ie-is+1, je-js+1, k1e-k1s+1, k2e-k2s+1, k3e-k3s+1 /) )

       CALL handle_error( 'get_variable_5d_real', 535, variable_name )

!
!-- Resort (reverse index order) and standardize (from 1 to N) output array

       DO  i = is, ie 
          DO  j = js, je
             DO  k1 = k1s, k1e
                DO  k2 = k2s, k2e
                   DO k3 = k3s, k3e
                      var(k3-k3s+1,k2-k2s+1,k1-k1s+1,j-js+1,i-is+1) = tmp(i,j,k1,k2,k3)
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO

       DEALLOCATE( tmp )

#endif

    END SUBROUTINE get_variable_5d_real


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads a 5D float variables from dynamic driver, such as time-dependent xy-, 
!> xz- or yz-boundary data as well as 5D initialization data. Please note, 
!> the passed arguments are start indices and number of elements in each 
!> dimension, which is in contrast to the other 3d versions where start- and 
!> end indices are passed. The different handling of 5D dynamic variables is 
!> due to its asymmetry for the u- and v component.
!> NOTE(1) - This subroutine is more flexible than get_variable_xd_real as it
!>           provides much better control over starting and count indices
!>           (ecc 20190418)
!> NOTE(2) - This subroutine is used specific for reading NC variable
!>           emission_values having a "z" dimension.  Said dimension
!>           is to be removed in the future and this subroutine shall
!>           be depreciated accordingly (ecc 20190418)
!------------------------------------------------------------------------------!

    SUBROUTINE get_variable_5d_real_dynamic( id, variable_name, var,                       &
                                             i1s, i2s, i3s, i4s, i5s,                      &
                                             count_1, count_2, count_3, count_4, count_5,  &
                                             par_access )

       USE indices
       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=*)          ::  variable_name  !< variable name

       LOGICAL                   ::  par_access     !< additional flag indicating parallel read

       INTEGER(iwp)              ::  count_1  !< # elements read in dimension 1 wrt file
       INTEGER(iwp)              ::  count_2  !< # elements read in dimension 2 wrt file
       INTEGER(iwp)              ::  count_3  !< # elements read in dimension 3 wrt file
       INTEGER(iwp)              ::  count_4  !< # elements read in dimension 4 wrt file
       INTEGER(iwp)              ::  count_5  !< # elements read in dimension 5 wrt file
       INTEGER(iwp)              ::  i1       !< index for dimension 1 on file
       INTEGER(iwp)              ::  i1s      !< starting index for dimension 1 hyperslab
       INTEGER(iwp)              ::  i2       !< index for dimension 2 on file
       INTEGER(iwp)              ::  i2s      !< starting index for dimension 2 hyperslab
       INTEGER(iwp)              ::  i3       !< index for dimension 3 on file
       INTEGER(iwp)              ::  i3s      !< starting index for dimension 3 hyperslab
       INTEGER(iwp)              ::  i4       !< index for dimension 4 on file
       INTEGER(iwp)              ::  i4s      !< starting index for dimension 4 hyperslab
       INTEGER(iwp)              ::  i5       !< index for dimension 5 on file
       INTEGER(iwp)              ::  i5s      !< starting index for dimension 5 hyperslab
       INTEGER(iwp)              ::  id_var   !< netCDF variable id (varid)
       INTEGER(iwp)              ::  lb1      !< lower bound of dimension 1 wrt file
       INTEGER(iwp)              ::  lb2      !< lower bound of dimension 2 wrt file
       INTEGER(iwp)              ::  lb3      !< lower bound of dimension 3 wrt file
       INTEGER(iwp)              ::  lb4      !< lower bound of dimension 4 wrt file
       INTEGER(iwp)              ::  lb5      !< lower bound of dimension 5 wrt file
       INTEGER(iwp)              ::  ub1      !< upper bound of dimension 1 wrt file
       INTEGER(iwp)              ::  ub2      !< upper bound of dimension 2 wrt file
       INTEGER(iwp)              ::  ub3      !< upper bound of dimension 3 wrt file
       INTEGER(iwp)              ::  ub4      !< upper bound of dimension 4 wrt file
       INTEGER(iwp)              ::  ub5      !< upper bound of dimension 5 wrt file
       INTEGER(iwp), INTENT(IN)  ::  id       !< netCDF file id (ncid)

       REAL(wp), DIMENSION(:,:,:,:,:), ALLOCATABLE    ::  tmp  !< temporary variable to read data
                                                               !< from file according is reverse
                                                               !< array index order
       REAL(wp), DIMENSION(:,:,:,:,:), INTENT(INOUT)  ::  var  !< input variable
       
#if defined( __netcdf )

!
!-- Inquire variable id.

       nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )

!
!-- Check for collective read-operation and set respective NetCDF flags if required. 
!-- Please note, in contrast to the other input routines where each PEs 
!-- reads its subdomain data, dynamic input data not by all PEs, only
!-- by those which encompass lateral model boundaries. Hence, collective
!-- read operations are only enabled for top-boundary data.

       IF ( collective_read  .AND.  par_access )  THEN

#if defined( __netcdf4_parallel )       
          nc_stat = NF90_VAR_PAR_ACCESS (id, id_var, NF90_COLLECTIVE)
#endif

       ENDIF

!
!-- Allocate temporary variable according to memory access on file. 
!-- Therefore, determine dimension bounds of input array. 

       lb1 = LBOUND(var,5)
       ub1 = UBOUND(var,5)
       lb2 = LBOUND(var,4)
       ub2 = UBOUND(var,4)
       lb3 = LBOUND(var,3)
       ub3 = UBOUND(var,3)
       lb4 = LBOUND(var,2)
       ub4 = UBOUND(var,2)
       lb5 = LBOUND(var,1)
       ub5 = UBOUND(var,1)
       ALLOCATE ( tmp(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4,lb5:ub5) )

!
!-- Get variable

       nc_stat = NF90_GET_VAR(  id, id_var, tmp,                                         &
                      start = (/ i1s,     i2s,     i3s,     i4s,     i5s     /),         &
                      count = (/ count_1, count_2, count_3, count_4, count_5 /) )

       CALL handle_error( 'get_variable_3d_real_dynamic', 537, variable_name )

!
!-- Assign temp array to output.  Note reverse index order 

       DO  i5 = lb5, ub5
          DO  i4 = lb4, ub4
             DO  i3 = lb3, ub3
                DO i2 = lb2, ub2
                   DO  i1 = lb1, ub1
                      var(i5,i4,i3,i2,i1) = tmp(i1,i2,i3,i4,i5)
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO

       DEALLOCATE( tmp )

#endif

    END SUBROUTINE get_variable_5d_real_dynamic


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Inquires the number of variables in a file
!------------------------------------------------------------------------------!
    SUBROUTINE inquire_num_variables( id, num_vars )

       USE indices
       USE pegrid

       IMPLICIT NONE

       INTEGER(iwp), INTENT(IN)      ::  id              !< file id
       INTEGER(iwp), INTENT(INOUT)   ::  num_vars        !< number of variables in a file
#if defined( __netcdf )

       nc_stat = NF90_INQUIRE( id, NVARIABLES = num_vars )
       CALL handle_error( 'inquire_num_variables', 539 )

#endif
    END SUBROUTINE inquire_num_variables


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Inquires the variable names belonging to a file.
!------------------------------------------------------------------------------!
    SUBROUTINE inquire_variable_names( id, var_names )

       USE indices
       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=*), DIMENSION(:), INTENT(INOUT) ::  var_names   !< return variable - variable names
       INTEGER(iwp)                                  ::  i           !< loop variable
       INTEGER(iwp), INTENT(IN)                      ::  id          !< file id
       INTEGER(iwp)                                  ::  num_vars    !< number of variables (unused return parameter)
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE       ::  varids      !< dummy array to strore variable ids temporarily
#if defined( __netcdf )

       ALLOCATE( varids(1:SIZE(var_names)) )
       nc_stat = NF90_INQ_VARIDS( id, NVARS = num_vars, VARIDS = varids )
       CALL handle_error( 'inquire_variable_names', 540 )

       DO  i = 1, SIZE(var_names)
          nc_stat = NF90_INQUIRE_VARIABLE( id, varids(i), NAME = var_names(i) )
          CALL handle_error( 'inquire_variable_names', 540 )
       ENDDO

       DEALLOCATE( varids )
#endif
    END SUBROUTINE inquire_variable_names

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Inquires the _FillValue settings of an integer variable.
!------------------------------------------------------------------------------!
    SUBROUTINE inquire_fill_value_int( id, var_name, nofill, fill_value )

       CHARACTER(LEN=*), INTENT(IN) ::  var_name    !< variable name

       INTEGER(iwp), INTENT(IN)  ::  id          !< file id
       INTEGER(iwp)              ::  nofill      !< flag indicating whether fill values are set or not
       INTEGER(iwp)              ::  fill_value  !< fill value
       INTEGER(iwp)              ::  id_var      !< netCDF variable id (varid)

#if defined( __netcdf )
       nc_stat = NF90_INQ_VARID( id, TRIM( var_name ), id_var )
       nc_stat = NF90_INQ_VAR_FILL(id, id_var, no_fill, fill_value )
#endif
!
!--    Further line is just to avoid compiler warnings. nofill might be used
!--    in future.
       IF ( nofill == 0  .OR.  nofill /= 0 )  CONTINUE

    END SUBROUTINE inquire_fill_value_int

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Inquires the _FillValue settings of a real variable.
!------------------------------------------------------------------------------!
    SUBROUTINE inquire_fill_value_real( id, var_name, nofill, fill_value )

       CHARACTER(LEN=*), INTENT(IN) ::  var_name    !< variable name

       INTEGER(iwp), INTENT(IN)  ::  id          !< file id
       INTEGER(iwp)              ::  nofill      !< flag indicating whether fill values are set or not
       INTEGER(iwp)              ::  id_var      !< netCDF variable id (varid)

#if defined( __imuk_old )
       INTEGER(iwp)              ::  fill_value_int  !< fill value workaround
#endif
       REAL(wp), INTENT(OUT)     ::  fill_value  !< fill value

#if defined( __netcdf )
       nc_stat = NF90_INQ_VARID( id, TRIM( var_name ), id_var )
#if defined( __imuk_old )
       nc_stat = NF90_INQ_VAR_FILL(id, id_var, no_fill, fill_value_int )
       fill_value = fill_value_int
#else
       nc_stat = NF90_INQ_VAR_FILL(id, id_var, no_fill, fill_value )
#endif
#endif
!
!--    Further line is just to avoid compiler warnings. nofill might be used
!--    in future.
       IF ( nofill == 0  .OR.  nofill /= 0 )  CONTINUE

    END SUBROUTINE inquire_fill_value_real

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Prints out a text message corresponding to the current status.
!------------------------------------------------------------------------------!
    SUBROUTINE handle_error( routine_name, errno, name )

       USE control_parameters,                                                 &
           ONLY:  message_string

       IMPLICIT NONE

       CHARACTER(LEN=6) ::  message_identifier !< string for the error number
       CHARACTER(LEN=*) ::  routine_name       !< routine name where the error happened
       CHARACTER(LEN=*), OPTIONAL ::  name     !< name of variable where reading failed

       INTEGER(iwp) ::  errno
#if defined( __netcdf )
       
       IF ( nc_stat /= NF90_NOERR )  THEN

          WRITE( message_identifier, '(''NC'',I4.4)' )  errno
          
          IF ( PRESENT( name ) )  THEN
             message_string = "Problem reading attribute/variable - " //       &
                              TRIM(name) // ": " //                            &
                              TRIM( NF90_STRERROR( nc_stat ) )
          ELSE
             message_string = TRIM( NF90_STRERROR( nc_stat ) )
          ENDIF

          CALL message( routine_name, message_identifier, 2, 2, myid, 6, 1 )

       ENDIF

#endif
    END SUBROUTINE handle_error


 END MODULE netcdf_data_input_mod
