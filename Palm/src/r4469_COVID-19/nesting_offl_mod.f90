!> @file nesting_offl_mod.f90
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
! $Id: nesting_offl_mod.f90 4457 2020-03-11 14:20:43Z raasch $
! use statement for exchange horiz added
! 
! 4360 2020-01-07 11:25:50Z suehring
! Bugfix, time coordinate is relative to origin_time rather than to 00:00:00 
! UTC.
! 
! 4346 2019-12-18 11:55:56Z motisi
! Introduction of wall_flags_total_0, which currently sets bits based on static
! topography information used in wall_flags_static_0
! 
! 4329 2019-12-10 15:46:36Z motisi
! Renamed wall_flags_0 to wall_flags_static_0
! 
! 4286 2019-10-30 16:01:14Z resler
! Fix wrong checks of time from dynamic driver in nesting_offl_mod
! 
! 4273 2019-10-24 13:40:54Z monakurppa
! Add a logical switch nesting_offline_chem
! 
! 4270 2019-10-23 10:46:20Z monakurppa
! Implement offline nesting for salsa variables.
! 
! 4231 2019-09-12 11:22:00Z suehring
! Bugfix in array deallocation
! 
! 4230 2019-09-11 13:58:14Z suehring
! Update mean chemistry profiles. These are also used for rayleigh damping. 
! 
! 4227 2019-09-10 18:04:34Z gronemeier
! implement new palm_date_time_mod
!
! - Data input moved into nesting_offl_mod
! - check rephrased
! - time variable is now relative to time_utc_init
! - Define module specific data type for offline nesting in nesting_offl_mod
! 
! 4182 2019-08-22 15:20:23Z scharf
! Corrected "Former revisions" section
! 
! 4169 2019-08-19 13:54:35Z suehring
! Additional check added.
! 
! 4168 2019-08-16 13:50:17Z suehring
! Replace function get_topography_top_index by topo_top_ind
! 
! 4125 2019-07-29 13:31:44Z suehring
! In order to enable netcdf parallel access, allocate dummy arrays for the 
! lateral boundary data on cores that actually do not belong to these 
! boundaries. 
! 
! 4079 2019-07-09 18:04:41Z suehring
! - Set boundary condition for w at nzt+1 at the lateral boundaries, even 
!   though these won't enter the numerical solution. However, due to the mass
!   conservation these values might some up to very large values which will
!   occur in the run-control file
! - Bugfix in offline nesting of chemical species
! - Do not set Neumann conditions for TKE and passive scalar
! 
! 4022 2019-06-12 11:52:39Z suehring
! Detection of boundary-layer depth in stable boundary layer on basis of 
! boundary data improved 
! Routine for boundary-layer depth calculation renamed and made public
! 
! 3987 2019-05-22 09:52:13Z kanani
! Introduce alternative switch for debug output during timestepping
! 
! 3964 2019-05-09 09:48:32Z suehring
! Ensure that veloctiy term in calculation of bulk Richardson number does not
! become zero
! 
! 3937 2019-04-29 15:09:07Z suehring
! Set boundary conditon on upper-left and upper-south grid point for the u- and
! v-component, respectively.
! 
! 3891 2019-04-12 17:52:01Z suehring
! Bugfix, do not overwrite lateral and top boundary data in case of restart 
! runs. 
! 
! 3885 2019-04-11 11:29:34Z kanani
! Changes related to global restructuring of location messages and introduction 
! of additional debug messages
! 
! 
! Do local data exchange for chemistry variables only when boundary data is  
! coming from dynamic file
! 
! 3737 2019-02-12 16:57:06Z suehring
! Introduce mesoscale nesting for chemical species
! 
! 3705 2019-01-29 19:56:39Z suehring
! Formatting adjustments
! 
! 3704 2019-01-29 19:51:41Z suehring
! Check implemented for offline nesting in child domain 
! 
! Initial Revision: 
! - separate offline nesting from large_scale_nudging_mod
! - revise offline nesting, adjust for usage of synthetic turbulence generator
! - adjust Rayleigh damping depending on the time-depending atmospheric
!   conditions
! 
! 
! Description:
! ------------
!> Offline nesting in larger-scale models. Boundary conditions for the simulation
!> are read from NetCDF file and are prescribed onto the respective arrays.
!> Further, a mass-flux correction is performed to maintain the mass balance. 
!--------------------------------------------------------------------------------!
 MODULE nesting_offl_mod

    USE arrays_3d,                                                             &
        ONLY:  dzw,                                                            &
               e,                                                              &
               diss,                                                           &
               pt,                                                             &
               pt_init,                                                        &
               q,                                                              &
               q_init,                                                         &
               rdf,                                                            &
               rdf_sc,                                                         &
               s,                                                              &
               u,                                                              &
               u_init,                                                         &
               ug,                                                             &
               v,                                                              &
               v_init,                                                         &
               vg,                                                             &
               w,                                                              &
               zu,                                                             &
               zw

    USE basic_constants_and_equations_mod,                                     &
           ONLY:  g,                                                           &
                  pi

    USE chem_modules,                                                          &
        ONLY:  chem_species, nesting_offline_chem

    USE control_parameters,                                                    &
        ONLY:  air_chemistry,                                                  & 
               bc_dirichlet_l,                                                 & 
               bc_dirichlet_n,                                                 &
               bc_dirichlet_r,                                                 &
               bc_dirichlet_s,                                                 &
               coupling_char,                                                  &
               dt_3d,                                                          &
               dz,                                                             &
               constant_diffusion,                                             &
               child_domain,                                                   &
               debug_output_timestep,                                          &
               end_time,                                                       &
               humidity,                                                       &
               initializing_actions,                                           &
               message_string,                                                 &
               nesting_offline,                                                &
               neutral,                                                        &
               passive_scalar,                                                 &
               rans_mode,                                                      &
               rans_tke_e,                                                     &
               rayleigh_damping_factor,                                        & 
               rayleigh_damping_height,                                        &
               salsa,                                                          &
               spinup_time,                                                    &
               time_since_reference_point,                                     &
               volume_flow

    USE cpulog,                                                                &
        ONLY:  cpu_log,                                                        &
               log_point,                                                      &
               log_point_s

    USE grid_variables

    USE indices,                                                               &
        ONLY:  nbgp, nx, nxl, nxlg, nxlu, nxr, nxrg, ny, nys,                  &
               nysv, nysg, nyn, nyng, nzb, nz, nzt,                            &
               topo_top_ind,                                                   &
               wall_flags_total_0

    USE kinds

    USE netcdf_data_input_mod,                                                 &
        ONLY:  check_existence,                                                &
               close_input_file,                                               &
               get_dimension_length,                                           &
               get_variable,                                                   &
               get_variable_pr,                                                &
               input_pids_dynamic,                                             &
               inquire_num_variables,                                          &
               inquire_variable_names,                                         &
               input_file_dynamic,                                             &
               num_var_pids,                                                   &
               open_read_file,                                                 &
               pids_id

    USE pegrid

    USE salsa_mod,                                                             &
        ONLY:  salsa_nesting_offl_bc,                                          &
               salsa_nesting_offl_init,                                        &
               salsa_nesting_offl_input

    IMPLICIT NONE

!
!-- Define data type for nesting in larger-scale models like COSMO.
!-- Data type comprises u, v, w, pt, and q at lateral and top boundaries.
    TYPE nest_offl_type

       CHARACTER(LEN=16) ::  char_l = 'ls_forcing_left_'  !< leading substring for variables at left boundary 
       CHARACTER(LEN=17) ::  char_n = 'ls_forcing_north_' !< leading substring for variables at north boundary  
       CHARACTER(LEN=17) ::  char_r = 'ls_forcing_right_' !< leading substring for variables at right boundary  
       CHARACTER(LEN=17) ::  char_s = 'ls_forcing_south_' !< leading substring for variables at south boundary 
       CHARACTER(LEN=15) ::  char_t = 'ls_forcing_top_'   !< leading substring for variables at top boundary 

       CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE ::  var_names         !< list of variable in dynamic input file
       CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE ::  var_names_chem_l  !< names of mesoscale nested chemistry variables at left boundary
       CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE ::  var_names_chem_n  !< names of mesoscale nested chemistry variables at north boundary
       CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE ::  var_names_chem_r  !< names of mesoscale nested chemistry variables at right boundary
       CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE ::  var_names_chem_s  !< names of mesoscale nested chemistry variables at south boundary
       CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE ::  var_names_chem_t  !< names of mesoscale nested chemistry variables at top boundary

       INTEGER(iwp) ::  nt     !< number of time levels in dynamic input file
       INTEGER(iwp) ::  nzu    !< number of vertical levels on scalar grid in dynamic input file
       INTEGER(iwp) ::  nzw    !< number of vertical levels on w grid in dynamic input file
       INTEGER(iwp) ::  tind   !< time index for reference time in mesoscale-offline nesting
       INTEGER(iwp) ::  tind_p !< time index for following time in mesoscale-offline nesting

       LOGICAL      ::  init         = .FALSE. !< flag indicating that offline nesting is already initialized

       LOGICAL, DIMENSION(:), ALLOCATABLE ::  chem_from_file_l !< flags inidicating whether left boundary data for chemistry is in dynamic input file  
       LOGICAL, DIMENSION(:), ALLOCATABLE ::  chem_from_file_n !< flags inidicating whether north boundary data for chemistry is in dynamic input file 
       LOGICAL, DIMENSION(:), ALLOCATABLE ::  chem_from_file_r !< flags inidicating whether right boundary data for chemistry is in dynamic input file 
       LOGICAL, DIMENSION(:), ALLOCATABLE ::  chem_from_file_s !< flags inidicating whether south boundary data for chemistry is in dynamic input file 
       LOGICAL, DIMENSION(:), ALLOCATABLE ::  chem_from_file_t !< flags inidicating whether top boundary data for chemistry is in dynamic input file 

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  surface_pressure !< time dependent surface pressure
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  time             !< time levels in dynamic input file
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  zu_atmos         !< vertical levels at scalar grid in dynamic input file
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  zw_atmos         !< vertical levels at w grid in dynamic input file

       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ug         !< domain-averaged geostrophic component
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  vg         !< domain-averaged geostrophic component

       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  u_left   !< u-component at left boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  v_left   !< v-component at left boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  w_left   !< w-component at left boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  q_left   !< mixing ratio at left boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  pt_left  !< potentital temperautre at left boundary

       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  u_north  !< u-component at north boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  v_north  !< v-component at north boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  w_north  !< w-component at north boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  q_north  !< mixing ratio at north boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  pt_north !< potentital temperautre at north boundary

       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  u_right  !< u-component at right boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  v_right  !< v-component at right boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  w_right  !< w-component at right boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  q_right  !< mixing ratio at right boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  pt_right !< potentital temperautre at right boundary

       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  u_south  !< u-component at south boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  v_south  !< v-component at south boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  w_south  !< w-component at south boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  q_south  !< mixing ratio at south boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  pt_south !< potentital temperautre at south boundary

       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  u_top    !< u-component at top boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  v_top    !< v-component at top boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  w_top    !< w-component at top boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  q_top    !< mixing ratio at top boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  pt_top   !< potentital temperautre at top boundary
       
       REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  chem_left   !< chemical species at left boundary
       REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  chem_north  !< chemical species at left boundary
       REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  chem_right  !< chemical species at left boundary
       REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  chem_south  !< chemical species at left boundary
       REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  chem_top    !< chemical species at left boundary

    END TYPE nest_offl_type

    REAL(wp) ::  fac_dt              !< interpolation factor
    REAL(wp) ::  zi_ribulk = 0.0_wp  !< boundary-layer depth according to bulk Richardson criterion, i.e. the height where Ri_bulk exceeds the critical
                                     !< bulk Richardson number of 0.2

    TYPE(nest_offl_type) ::  nest_offl  !< data structure for data input at lateral and top boundaries (provided by Inifor) 
    
    SAVE
    PRIVATE
!
!-- Public subroutines
    PUBLIC nesting_offl_bc,                                                    &
           nesting_offl_calc_zi,                                               &
           nesting_offl_check_parameters,                                      &
           nesting_offl_geostrophic_wind,                                      &
           nesting_offl_header,                                                &
           nesting_offl_init,                                                  &
           nesting_offl_input,                                                 &
           nesting_offl_interpolation_factor,                                  &
           nesting_offl_mass_conservation,                                     &
           nesting_offl_parin 
!
!-- Public variables
    PUBLIC zi_ribulk   

    INTERFACE nesting_offl_bc
       MODULE PROCEDURE nesting_offl_bc
    END INTERFACE nesting_offl_bc
    
    INTERFACE nesting_offl_calc_zi
       MODULE PROCEDURE nesting_offl_calc_zi
    END INTERFACE nesting_offl_calc_zi
    
    INTERFACE nesting_offl_check_parameters
       MODULE PROCEDURE nesting_offl_check_parameters
    END INTERFACE nesting_offl_check_parameters

    INTERFACE nesting_offl_geostrophic_wind
       MODULE PROCEDURE nesting_offl_geostrophic_wind
    END INTERFACE nesting_offl_geostrophic_wind
    
    INTERFACE nesting_offl_header
       MODULE PROCEDURE nesting_offl_header
    END INTERFACE nesting_offl_header
    
    INTERFACE nesting_offl_init
       MODULE PROCEDURE nesting_offl_init
    END INTERFACE nesting_offl_init

    INTERFACE nesting_offl_input
       MODULE PROCEDURE nesting_offl_input
    END INTERFACE nesting_offl_input

    INTERFACE nesting_offl_interpolation_factor
       MODULE PROCEDURE nesting_offl_interpolation_factor
    END INTERFACE nesting_offl_interpolation_factor
           
    INTERFACE nesting_offl_mass_conservation
       MODULE PROCEDURE nesting_offl_mass_conservation
    END INTERFACE nesting_offl_mass_conservation
    
    INTERFACE nesting_offl_parin
       MODULE PROCEDURE nesting_offl_parin
    END INTERFACE nesting_offl_parin

 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads data at lateral and top boundaries derived from larger-scale model.
!------------------------------------------------------------------------------!
    SUBROUTINE nesting_offl_input

       INTEGER(iwp) ::  n   !< running index for chemistry variables
       INTEGER(iwp) ::  t   !< running index time dimension

!
!--    Initialize INIFOR forcing in first call. 
       IF ( .NOT. nest_offl%init )  THEN
#if defined ( __netcdf )
!
!--       Open file in read-only mode
          CALL open_read_file( TRIM( input_file_dynamic ) //                   &
                               TRIM( coupling_char ), pids_id )
!
!--       At first, inquire all variable names.
          CALL inquire_num_variables( pids_id, num_var_pids )
!
!--       Allocate memory to store variable names.
          ALLOCATE( nest_offl%var_names(1:num_var_pids) )
          CALL inquire_variable_names( pids_id, nest_offl%var_names )
!
!--       Read time dimension, allocate memory and finally read time array
          CALL get_dimension_length( pids_id, nest_offl%nt, 'time' )

          IF ( check_existence( nest_offl%var_names, 'time' ) )  THEN
             ALLOCATE( nest_offl%time(0:nest_offl%nt-1) )
             CALL get_variable( pids_id, 'time', nest_offl%time )
          ENDIF
!
!--       Read vertical dimension of scalar und w grid
          CALL get_dimension_length( pids_id, nest_offl%nzu, 'z' )
          CALL get_dimension_length( pids_id, nest_offl%nzw, 'zw' )

          IF ( check_existence( nest_offl%var_names, 'z' ) )  THEN
             ALLOCATE( nest_offl%zu_atmos(1:nest_offl%nzu) )
             CALL get_variable( pids_id, 'z', nest_offl%zu_atmos )
          ENDIF
          IF ( check_existence( nest_offl%var_names, 'zw' ) )  THEN
             ALLOCATE( nest_offl%zw_atmos(1:nest_offl%nzw) )
             CALL get_variable( pids_id, 'zw', nest_offl%zw_atmos )
          ENDIF
!
!--       Read surface pressure
          IF ( check_existence( nest_offl%var_names,                           &
                                'surface_forcing_surface_pressure' ) )  THEN
             ALLOCATE( nest_offl%surface_pressure(0:nest_offl%nt-1) )
             CALL get_variable( pids_id,                                       &
                                'surface_forcing_surface_pressure',            &
                                nest_offl%surface_pressure )
          ENDIF
!
!--       Close input file
          CALL close_input_file( pids_id )
#endif
       ENDIF
!
!--    Check if dynamic driver data input is required.
       IF ( nest_offl%time(nest_offl%tind_p) <=                                &
            MAX( time_since_reference_point, 0.0_wp)  .OR.                     &
            .NOT.  nest_offl%init )  THEN
          CONTINUE
!
!--    Return otherwise
       ELSE
          RETURN
       ENDIF
!
!--    CPU measurement
       CALL cpu_log( log_point_s(86), 'NetCDF input forcing', 'start' )

!
!--    Obtain time index for current point in time. Note, the time coordinate
!--    in the input file is always relative to the initial time in UTC, i.e.
!--    the time coordinate always starts at 0.0 even if the initial UTC is e.g.
!--    7200.0. Further, since time_since_reference_point is negativ here when 
!--    spinup is applied, use MAX function to obtain correct time index. 
       nest_offl%tind = MINLOC( ABS( nest_offl%time -                          &
                                     MAX( time_since_reference_point, 0.0_wp)  &
                                   ), DIM = 1 ) - 1
       nest_offl%tind_p = nest_offl%tind + 1
!
!--    Open file in read-only mode
#if defined ( __netcdf )
       CALL open_read_file( TRIM( input_file_dynamic ) //                      &
                            TRIM( coupling_char ), pids_id )
!
!--    Read geostrophic wind components
       DO  t = nest_offl%tind, nest_offl%tind_p
          CALL get_variable_pr( pids_id, 'ls_forcing_ug', t+1,                 &
                                nest_offl%ug(t-nest_offl%tind,nzb+1:nzt) )
          CALL get_variable_pr( pids_id, 'ls_forcing_vg', t+1,                 &
                                nest_offl%vg(t-nest_offl%tind,nzb+1:nzt) )
       ENDDO
!
!--    Read data at lateral and top boundaries. Please note, at left and
!--    right domain boundary, yz-layers are read for u, v, w, pt and q.
!--    For the v-component, the data starts at nysv, while for the other
!--    quantities the data starts at nys. This is equivalent at the north
!--    and south domain boundary for the u-component. 
!--    Note, lateral data is also accessed by parallel IO, which is the reason
!--    why different arguments are passed depending on the boundary control
!--    flags. Cores that do not belong to the respective boundary just make
!--    a dummy read with count = 0, just in order to participate the collective
!--    operation. 
!--    Read data for western boundary    
       CALL get_variable( pids_id, 'ls_forcing_left_u',                        &
                          nest_offl%u_left,                                    & ! array to be read
                          MERGE( nys+1, 1, bc_dirichlet_l),                    & ! start index y direction
                          MERGE( nzb+1, 1, bc_dirichlet_l),                    & ! start index z direction
                          MERGE( nest_offl%tind+1, 1, bc_dirichlet_l),         & ! start index time dimension
                          MERGE( nyn-nys+1, 0, bc_dirichlet_l),                & ! number of elements along y
                          MERGE( nest_offl%nzu, 0, bc_dirichlet_l),            & ! number of elements alogn z
                          MERGE( 2, 0, bc_dirichlet_l),                        & ! number of time steps (2 or 0)
                          .TRUE. )                                               ! parallel IO when compiled accordingly
      
       CALL get_variable( pids_id, 'ls_forcing_left_v',                        &
                          nest_offl%v_left,                                    &
                          MERGE( nysv, 1, bc_dirichlet_l),                     &
                          MERGE( nzb+1, 1, bc_dirichlet_l),                    &
                          MERGE( nest_offl%tind+1, 1, bc_dirichlet_l),         &
                          MERGE( nyn-nysv+1, 0, bc_dirichlet_l),               &
                          MERGE( nest_offl%nzu, 0, bc_dirichlet_l),            &
                          MERGE( 2, 0, bc_dirichlet_l),                        &
                          .TRUE. )                                        

       CALL get_variable( pids_id, 'ls_forcing_left_w',                        &
                          nest_offl%w_left,                                    &
                          MERGE( nys+1, 1, bc_dirichlet_l),                    &
                          MERGE( nzb+1, 1, bc_dirichlet_l),                    &
                          MERGE( nest_offl%tind+1, 1, bc_dirichlet_l),         &
                          MERGE( nyn-nys+1, 0, bc_dirichlet_l),                &
                          MERGE( nest_offl%nzw, 0, bc_dirichlet_l),            &
                          MERGE( 2, 0, bc_dirichlet_l),                        &
                          .TRUE. )   

       IF ( .NOT. neutral )  THEN
          CALL get_variable( pids_id, 'ls_forcing_left_pt',                    &
                             nest_offl%pt_left,                                &
                             MERGE( nys+1, 1, bc_dirichlet_l),                 &
                             MERGE( nzb+1, 1, bc_dirichlet_l),                 &
                             MERGE( nest_offl%tind+1, 1, bc_dirichlet_l),      &
                             MERGE( nyn-nys+1, 0, bc_dirichlet_l),             &
                             MERGE( nest_offl%nzu, 0, bc_dirichlet_l),         &
                             MERGE( 2, 0, bc_dirichlet_l),                     &
                             .TRUE. )
       ENDIF

       IF ( humidity )  THEN
          CALL get_variable( pids_id, 'ls_forcing_left_qv',                    &
                             nest_offl%q_left,                                 &
                             MERGE( nys+1, 1, bc_dirichlet_l),                 &
                             MERGE( nzb+1, 1, bc_dirichlet_l),                 &
                             MERGE( nest_offl%tind+1, 1, bc_dirichlet_l),      &
                             MERGE( nyn-nys+1, 0, bc_dirichlet_l),             &
                             MERGE( nest_offl%nzu, 0, bc_dirichlet_l),         &
                             MERGE( 2, 0, bc_dirichlet_l),                     &
                             .TRUE. )
       ENDIF

       IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
          DO  n = 1, UBOUND(nest_offl%var_names_chem_l, 1)
             IF ( check_existence( nest_offl%var_names,                        &
                                   nest_offl%var_names_chem_l(n) ) )  THEN
                CALL get_variable( pids_id,                                    &
                           TRIM( nest_offl%var_names_chem_l(n) ),              &
                           nest_offl%chem_left(:,:,:,n),                       &
                           MERGE( nys+1, 1, bc_dirichlet_l),                   &
                           MERGE( nzb+1, 1, bc_dirichlet_l),                   &
                           MERGE( nest_offl%tind+1, 1, bc_dirichlet_l),        &
                           MERGE( nyn-nys+1, 0, bc_dirichlet_l),               &
                           MERGE( nest_offl%nzu, 0, bc_dirichlet_l),           &
                           MERGE( 2, 0, bc_dirichlet_l),                       &
                           .TRUE. )
                nest_offl%chem_from_file_l(n) = .TRUE.
             ENDIF
          ENDDO
       ENDIF
!
!--    Read data for eastern boundary    
       CALL get_variable( pids_id, 'ls_forcing_right_u',                       &
                          nest_offl%u_right,                                   &
                          MERGE( nys+1, 1, bc_dirichlet_r),                    &
                          MERGE( nzb+1, 1, bc_dirichlet_r),                    &
                          MERGE( nest_offl%tind+1, 1, bc_dirichlet_r),         &
                          MERGE( nyn-nys+1, 0, bc_dirichlet_r),                &
                          MERGE( nest_offl%nzu, 0, bc_dirichlet_r),            &
                          MERGE( 2, 0, bc_dirichlet_r),                        &
                          .TRUE. )                                             
      
       CALL get_variable( pids_id, 'ls_forcing_right_v',                       &
                          nest_offl%v_right,                                   &
                          MERGE( nysv, 1, bc_dirichlet_r),                     &
                          MERGE( nzb+1, 1, bc_dirichlet_r),                    &
                          MERGE( nest_offl%tind+1, 1, bc_dirichlet_r),         &
                          MERGE( nyn-nysv+1, 0, bc_dirichlet_r),               &
                          MERGE( nest_offl%nzu, 0, bc_dirichlet_r),            &
                          MERGE( 2, 0, bc_dirichlet_r),                        &
                          .TRUE. )                                             

       CALL get_variable( pids_id, 'ls_forcing_right_w',                       &
                          nest_offl%w_right,                                   &
                          MERGE( nys+1, 1, bc_dirichlet_r),                    &
                          MERGE( nzb+1, 1, bc_dirichlet_r),                    &
                          MERGE( nest_offl%tind+1, 1, bc_dirichlet_r),         &
                          MERGE( nyn-nys+1, 0, bc_dirichlet_r),                &
                          MERGE( nest_offl%nzw, 0, bc_dirichlet_r),            &
                          MERGE( 2, 0, bc_dirichlet_r),                        &
                          .TRUE. )   

       IF ( .NOT. neutral )  THEN
          CALL get_variable( pids_id, 'ls_forcing_right_pt',                   &
                             nest_offl%pt_right,                               &
                             MERGE( nys+1, 1, bc_dirichlet_r),                 &
                             MERGE( nzb+1, 1, bc_dirichlet_r),                 &
                             MERGE( nest_offl%tind+1, 1, bc_dirichlet_r),      &
                             MERGE( nyn-nys+1, 0, bc_dirichlet_r),             &
                             MERGE( nest_offl%nzu, 0, bc_dirichlet_r),         &
                             MERGE( 2, 0, bc_dirichlet_r),                     &
                             .TRUE. )
       ENDIF

       IF ( humidity )  THEN
          CALL get_variable( pids_id, 'ls_forcing_right_qv',                   &
                             nest_offl%q_right,                                &
                             MERGE( nys+1, 1, bc_dirichlet_r),                 &
                             MERGE( nzb+1, 1, bc_dirichlet_r),                 &
                             MERGE( nest_offl%tind+1, 1, bc_dirichlet_r),      &
                             MERGE( nyn-nys+1, 0, bc_dirichlet_r),             &
                             MERGE( nest_offl%nzu, 0, bc_dirichlet_r),         &
                             MERGE( 2, 0, bc_dirichlet_r),                     &
                             .TRUE. )
       ENDIF

       IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
          DO  n = 1, UBOUND(nest_offl%var_names_chem_r, 1)
             IF ( check_existence( nest_offl%var_names,                        &
                                   nest_offl%var_names_chem_r(n) ) )  THEN
                CALL get_variable( pids_id,                                    &
                           TRIM( nest_offl%var_names_chem_r(n) ),              &
                           nest_offl%chem_right(:,:,:,n),                      &
                           MERGE( nys+1, 1, bc_dirichlet_r),                   &
                           MERGE( nzb+1, 1, bc_dirichlet_r),                   &
                           MERGE( nest_offl%tind+1, 1, bc_dirichlet_r),        &
                           MERGE( nyn-nys+1, 0, bc_dirichlet_r),               &
                           MERGE( nest_offl%nzu, 0, bc_dirichlet_r),           &
                           MERGE( 2, 0, bc_dirichlet_r),                       &
                           .TRUE. )
                nest_offl%chem_from_file_r(n) = .TRUE.
             ENDIF
          ENDDO
       ENDIF
!
!--    Read data for northern boundary
       CALL get_variable( pids_id, 'ls_forcing_north_u',                       & ! array to be read
                          nest_offl%u_north,                                   & ! start index x direction
                          MERGE( nxlu, 1, bc_dirichlet_n ),                    & ! start index z direction
                          MERGE( nzb+1, 1, bc_dirichlet_n ),                   & ! start index time dimension
                          MERGE( nest_offl%tind+1, 1, bc_dirichlet_n ),        & ! number of elements along x
                          MERGE( nxr-nxlu+1, 0, bc_dirichlet_n ),              & ! number of elements alogn z 
                          MERGE( nest_offl%nzu, 0, bc_dirichlet_n ),           & ! number of time steps (2 or 0)
                          MERGE( 2, 0, bc_dirichlet_n ),                       & ! parallel IO when compiled accordingly
                          .TRUE. )

       CALL get_variable( pids_id, 'ls_forcing_north_v',                       & ! array to be read
                          nest_offl%v_north,                                   & ! start index x direction
                          MERGE( nxl+1, 1, bc_dirichlet_n ),                   & ! start index z direction
                          MERGE( nzb+1, 1, bc_dirichlet_n ),                   & ! start index time dimension
                          MERGE( nest_offl%tind+1, 1, bc_dirichlet_n ),        & ! number of elements along x
                          MERGE( nxr-nxl+1, 0, bc_dirichlet_n ),               & ! number of elements alogn z 
                          MERGE( nest_offl%nzu, 0, bc_dirichlet_n ),           & ! number of time steps (2 or 0)
                          MERGE( 2, 0, bc_dirichlet_n ),                       & ! parallel IO when compiled accordingly
                          .TRUE. )

       CALL get_variable( pids_id, 'ls_forcing_north_w',                       & ! array to be read
                          nest_offl%w_north,                                   & ! start index x direction
                          MERGE( nxl+1, 1, bc_dirichlet_n ),                   & ! start index z direction
                          MERGE( nzb+1, 1, bc_dirichlet_n ),                   & ! start index time dimension
                          MERGE( nest_offl%tind+1, 1, bc_dirichlet_n ),        & ! number of elements along x
                          MERGE( nxr-nxl+1, 0, bc_dirichlet_n ),               & ! number of elements alogn z 
                          MERGE( nest_offl%nzw, 0, bc_dirichlet_n ),           & ! number of time steps (2 or 0)
                          MERGE( 2, 0, bc_dirichlet_n ),                       & ! parallel IO when compiled accordingly
                          .TRUE. )

       IF ( .NOT. neutral )  THEN
          CALL get_variable( pids_id, 'ls_forcing_north_pt',                   & ! array to be read
                             nest_offl%pt_north,                               & ! start index x direction
                             MERGE( nxl+1, 1, bc_dirichlet_n ),                & ! start index z direction
                             MERGE( nzb+1, 1, bc_dirichlet_n ),                & ! start index time dimension
                             MERGE( nest_offl%tind+1, 1, bc_dirichlet_n ),     & ! number of elements along x
                             MERGE( nxr-nxl+1, 0, bc_dirichlet_n ),            & ! number of elements alogn z 
                             MERGE( nest_offl%nzu, 0, bc_dirichlet_n ),        & ! number of time steps (2 or 0)
                             MERGE( 2, 0, bc_dirichlet_n ),                    & ! parallel IO when compiled accordingly
                             .TRUE. )
       ENDIF
       IF ( humidity )  THEN
          CALL get_variable( pids_id, 'ls_forcing_north_qv',                   & ! array to be read
                             nest_offl%q_north,                                & ! start index x direction
                             MERGE( nxl+1, 1, bc_dirichlet_n ),                & ! start index z direction
                             MERGE( nzb+1, 1, bc_dirichlet_n ),                & ! start index time dimension
                             MERGE( nest_offl%tind+1, 1, bc_dirichlet_n ),     & ! number of elements along x
                             MERGE( nxr-nxl+1, 0, bc_dirichlet_n ),            & ! number of elements alogn z 
                             MERGE( nest_offl%nzu, 0, bc_dirichlet_n ),        & ! number of time steps (2 or 0)
                             MERGE( 2, 0, bc_dirichlet_n ),                    & ! parallel IO when compiled accordingly
                             .TRUE. )
       ENDIF

       IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
          DO  n = 1, UBOUND(nest_offl%var_names_chem_n, 1)
             IF ( check_existence( nest_offl%var_names,                        &
                                   nest_offl%var_names_chem_n(n) ) )  THEN
                CALL get_variable( pids_id,                                    &
                           TRIM( nest_offl%var_names_chem_n(n) ),              &
                           nest_offl%chem_north(:,:,:,n),                      &
                           MERGE( nxl+1, 1, bc_dirichlet_n ),                  &
                           MERGE( nzb+1, 1, bc_dirichlet_n ),                  &
                           MERGE( nest_offl%tind+1, 1, bc_dirichlet_n ),       &
                           MERGE( nxr-nxl+1, 0, bc_dirichlet_n ),              &
                           MERGE( nest_offl%nzu, 0, bc_dirichlet_n ),          &
                           MERGE( 2, 0, bc_dirichlet_n ),                      &
                           .TRUE. )
                nest_offl%chem_from_file_n(n) = .TRUE.
             ENDIF
          ENDDO
       ENDIF
!
!--    Read data for southern boundary
       CALL get_variable( pids_id, 'ls_forcing_south_u',                       & ! array to be read
                          nest_offl%u_south,                                   & ! start index x direction
                          MERGE( nxlu, 1, bc_dirichlet_s ),                    & ! start index z direction
                          MERGE( nzb+1, 1, bc_dirichlet_s ),                   & ! start index time dimension
                          MERGE( nest_offl%tind+1, 1, bc_dirichlet_s ),        & ! number of elements along x
                          MERGE( nxr-nxlu+1, 0, bc_dirichlet_s ),              & ! number of elements alogn z 
                          MERGE( nest_offl%nzu, 0, bc_dirichlet_s ),           & ! number of time steps (2 or 0)
                          MERGE( 2, 0, bc_dirichlet_s ),                       & ! parallel IO when compiled accordingly
                          .TRUE. )

       CALL get_variable( pids_id, 'ls_forcing_south_v',                       & ! array to be read
                          nest_offl%v_south,                                   & ! start index x direction
                          MERGE( nxl+1, 1, bc_dirichlet_s ),                   & ! start index z direction
                          MERGE( nzb+1, 1, bc_dirichlet_s ),                   & ! start index time dimension
                          MERGE( nest_offl%tind+1, 1, bc_dirichlet_s ),        & ! number of elements along x
                          MERGE( nxr-nxl+1, 0, bc_dirichlet_s ),               & ! number of elements alogn z 
                          MERGE( nest_offl%nzu, 0, bc_dirichlet_s ),           & ! number of time steps (2 or 0)
                          MERGE( 2, 0, bc_dirichlet_s ),                       & ! parallel IO when compiled accordingly
                          .TRUE. )

       CALL get_variable( pids_id, 'ls_forcing_south_w',                       & ! array to be read
                          nest_offl%w_south,                                   & ! start index x direction
                          MERGE( nxl+1, 1, bc_dirichlet_s ),                   & ! start index z direction
                          MERGE( nzb+1, 1, bc_dirichlet_s ),                   & ! start index time dimension
                          MERGE( nest_offl%tind+1, 1, bc_dirichlet_s ),        & ! number of elements along x
                          MERGE( nxr-nxl+1, 0, bc_dirichlet_s ),               & ! number of elements alogn z 
                          MERGE( nest_offl%nzw, 0, bc_dirichlet_s ),           & ! number of time steps (2 or 0)
                          MERGE( 2, 0, bc_dirichlet_s ),                       & ! parallel IO when compiled accordingly
                          .TRUE. )

       IF ( .NOT. neutral )  THEN
          CALL get_variable( pids_id, 'ls_forcing_south_pt',                   & ! array to be read
                             nest_offl%pt_south,                               & ! start index x direction
                             MERGE( nxl+1, 1, bc_dirichlet_s ),                & ! start index z direction
                             MERGE( nzb+1, 1, bc_dirichlet_s ),                & ! start index time dimension
                             MERGE( nest_offl%tind+1, 1, bc_dirichlet_s ),     & ! number of elements along x
                             MERGE( nxr-nxl+1, 0, bc_dirichlet_s ),            & ! number of elements alogn z 
                             MERGE( nest_offl%nzu, 0, bc_dirichlet_s ),        & ! number of time steps (2 or 0)
                             MERGE( 2, 0, bc_dirichlet_s ),                    & ! parallel IO when compiled accordingly
                             .TRUE. )
       ENDIF
       IF ( humidity )  THEN
          CALL get_variable( pids_id, 'ls_forcing_south_qv',                   & ! array to be read
                             nest_offl%q_south,                                & ! start index x direction
                             MERGE( nxl+1, 1, bc_dirichlet_s ),                & ! start index z direction
                             MERGE( nzb+1, 1, bc_dirichlet_s ),                & ! start index time dimension
                             MERGE( nest_offl%tind+1, 1, bc_dirichlet_s ),     & ! number of elements along x
                             MERGE( nxr-nxl+1, 0, bc_dirichlet_s ),            & ! number of elements alogn z 
                             MERGE( nest_offl%nzu, 0, bc_dirichlet_s ),        & ! number of time steps (2 or 0)
                             MERGE( 2, 0, bc_dirichlet_s ),                    & ! parallel IO when compiled accordingly
                             .TRUE. )
       ENDIF

       IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
          DO  n = 1, UBOUND(nest_offl%var_names_chem_s, 1)
             IF ( check_existence( nest_offl%var_names,                        &
                                   nest_offl%var_names_chem_s(n) ) )  THEN
                CALL get_variable( pids_id,                                    &
                           TRIM( nest_offl%var_names_chem_s(n) ),              &
                           nest_offl%chem_south(:,:,:,n),                      &
                           MERGE( nxl+1, 1, bc_dirichlet_s ),                  &
                           MERGE( nzb+1, 1, bc_dirichlet_s ),                  &
                           MERGE( nest_offl%tind+1, 1, bc_dirichlet_s ),       &
                           MERGE( nxr-nxl+1, 0, bc_dirichlet_s ),              &
                           MERGE( nest_offl%nzu, 0, bc_dirichlet_s ),          &
                           MERGE( 2, 0, bc_dirichlet_s ),                      &
                           .TRUE. )
                nest_offl%chem_from_file_s(n) = .TRUE.
             ENDIF
          ENDDO
       ENDIF
!
!--    Top boundary
       CALL get_variable( pids_id, 'ls_forcing_top_u',                         &
                             nest_offl%u_top(0:1,nys:nyn,nxlu:nxr),            &
                             nxlu, nys+1, nest_offl%tind+1,                    &
                             nxr-nxlu+1, nyn-nys+1, 2, .TRUE. )

       CALL get_variable( pids_id, 'ls_forcing_top_v',                         &
                             nest_offl%v_top(0:1,nysv:nyn,nxl:nxr),            &
                             nxl+1, nysv, nest_offl%tind+1,                    &
                             nxr-nxl+1, nyn-nysv+1, 2, .TRUE. )

       CALL get_variable( pids_id, 'ls_forcing_top_w',                         &
                             nest_offl%w_top(0:1,nys:nyn,nxl:nxr),             &
                             nxl+1, nys+1, nest_offl%tind+1,                   &
                             nxr-nxl+1, nyn-nys+1, 2, .TRUE. )

       IF ( .NOT. neutral )  THEN
          CALL get_variable( pids_id, 'ls_forcing_top_pt',                     &
                                nest_offl%pt_top(0:1,nys:nyn,nxl:nxr),         &
                                nxl+1, nys+1, nest_offl%tind+1,                &
                                nxr-nxl+1, nyn-nys+1, 2, .TRUE. )
       ENDIF
       IF ( humidity )  THEN
          CALL get_variable( pids_id, 'ls_forcing_top_qv',                     &
                                nest_offl%q_top(0:1,nys:nyn,nxl:nxr),          &
                                nxl+1, nys+1, nest_offl%tind+1,                &
                                nxr-nxl+1, nyn-nys+1, 2, .TRUE. )
       ENDIF

       IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
          DO  n = 1, UBOUND(nest_offl%var_names_chem_t, 1)
             IF ( check_existence( nest_offl%var_names,                        &
                                   nest_offl%var_names_chem_t(n) ) )  THEN
                CALL get_variable( pids_id,                                    &
                              TRIM( nest_offl%var_names_chem_t(n) ),           &
                              nest_offl%chem_top(0:1,nys:nyn,nxl:nxr,n),       &
                              nxl+1, nys+1, nest_offl%tind+1,                  &
                              nxr-nxl+1, nyn-nys+1, 2, .TRUE. )
                nest_offl%chem_from_file_t(n) = .TRUE.
             ENDIF
          ENDDO
       ENDIF

!
!--    Close input file
       CALL close_input_file( pids_id )
#endif
!
!--    Set control flag to indicate that boundary data has been initially
!--    input.
       nest_offl%init = .TRUE.
!
!--    Call offline nesting for salsa
       IF ( salsa )  CALL salsa_nesting_offl_input
!
!--    End of CPU measurement
       CALL cpu_log( log_point_s(86), 'NetCDF input forcing', 'stop' )

    END SUBROUTINE nesting_offl_input


!------------------------------------------------------------------------------!
! Description:
! ------------
!> In this subroutine a constant mass within the model domain is guaranteed. 
!> Larger-scale models may be based on a compressible equation system, which is
!> not consistent with PALMs incompressible equation system. In order to avoid
!> a decrease or increase of mass during the simulation, non-divergent flow
!> through the lateral and top boundaries is compensated by the vertical wind 
!> component at the top boundary. 
!------------------------------------------------------------------------------!
    SUBROUTINE nesting_offl_mass_conservation

       INTEGER(iwp) ::  i !< grid index in x-direction
       INTEGER(iwp) ::  j !< grid index in y-direction
       INTEGER(iwp) ::  k !< grid index in z-direction

       REAL(wp) ::  d_area_t                        !< inverse of the total area of the horizontal model domain
       REAL(wp) ::  w_correct                       !< vertical velocity increment required to compensate non-divergent flow through the boundaries
       REAL(wp), DIMENSION(1:3) ::  volume_flow_l   !< local volume flow


       IF ( debug_output_timestep )  CALL debug_message( 'nesting_offl_mass_conservation', 'start' )

       CALL  cpu_log( log_point(58), 'offline nesting', 'start' )
       
       volume_flow   = 0.0_wp
       volume_flow_l = 0.0_wp

       d_area_t = 1.0_wp / ( ( nx + 1 ) * dx * ( ny + 1 ) * dy )

       IF ( bc_dirichlet_l )  THEN
          i = nxl
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                volume_flow_l(1) = volume_flow_l(1) + u(k,j,i) * dzw(k) * dy   &
                                     * MERGE( 1.0_wp, 0.0_wp,                  &
                                       BTEST( wall_flags_total_0(k,j,i), 1 ) )
             ENDDO
          ENDDO
       ENDIF
       IF ( bc_dirichlet_r )  THEN
          i = nxr+1
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                volume_flow_l(1) = volume_flow_l(1) - u(k,j,i) * dzw(k) * dy   &
                                     * MERGE( 1.0_wp, 0.0_wp,                  &
                                       BTEST( wall_flags_total_0(k,j,i), 1 ) )
             ENDDO
          ENDDO
       ENDIF
       IF ( bc_dirichlet_s )  THEN
          j = nys
          DO  i = nxl, nxr
             DO  k = nzb+1, nzt
                volume_flow_l(2) = volume_flow_l(2) + v(k,j,i) * dzw(k) * dx   &
                                     * MERGE( 1.0_wp, 0.0_wp,                  &
                                       BTEST( wall_flags_total_0(k,j,i), 2 ) )
             ENDDO
          ENDDO
       ENDIF
       IF ( bc_dirichlet_n )  THEN
          j = nyn+1
          DO  i = nxl, nxr
             DO  k = nzb+1, nzt
                volume_flow_l(2) = volume_flow_l(2) - v(k,j,i) * dzw(k) * dx   &
                                     * MERGE( 1.0_wp, 0.0_wp,                  &
                                       BTEST( wall_flags_total_0(k,j,i), 2 ) )
             ENDDO
          ENDDO
       ENDIF
!
!--    Top boundary
       k = nzt
       DO  i = nxl, nxr
          DO  j = nys, nyn
             volume_flow_l(3) = volume_flow_l(3) - w(k,j,i) * dx * dy
          ENDDO
       ENDDO

#if defined( __parallel )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( volume_flow_l, volume_flow, 3, MPI_REAL, MPI_SUM,   &
                           comm2d, ierr )
#else
       volume_flow = volume_flow_l
#endif

       w_correct = SUM( volume_flow ) * d_area_t

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzt, nzt + 1
                w(k,j,i) = w(k,j,i) + w_correct
             ENDDO
          ENDDO
       ENDDO
       
       CALL  cpu_log( log_point(58), 'offline nesting', 'stop' )

       IF ( debug_output_timestep )  CALL debug_message( 'nesting_offl_mass_conservation', 'end' )

    END SUBROUTINE nesting_offl_mass_conservation


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Set the lateral and top boundary conditions in case the PALM domain is 
!> nested offline in a mesoscale model. Further, average boundary data and
!> determine mean profiles, further used for correct damping in the sponge
!> layer. 
!------------------------------------------------------------------------------!
    SUBROUTINE nesting_offl_bc                     

       USE exchange_horiz_mod,                                                    &
           ONLY:  exchange_horiz

       INTEGER(iwp) ::  i !< running index x-direction
       INTEGER(iwp) ::  j !< running index y-direction
       INTEGER(iwp) ::  k !< running index z-direction
       INTEGER(iwp) ::  n !< running index for chemical species
       
       REAL(wp), DIMENSION(nzb:nzt+1) ::  pt_ref   !< reference profile for potential temperature
       REAL(wp), DIMENSION(nzb:nzt+1) ::  pt_ref_l !< reference profile for potential temperature on subdomain
       REAL(wp), DIMENSION(nzb:nzt+1) ::  q_ref    !< reference profile for mixing ratio
       REAL(wp), DIMENSION(nzb:nzt+1) ::  q_ref_l  !< reference profile for mixing ratio on subdomain
       REAL(wp), DIMENSION(nzb:nzt+1) ::  u_ref    !< reference profile for u-component
       REAL(wp), DIMENSION(nzb:nzt+1) ::  u_ref_l  !< reference profile for u-component on subdomain
       REAL(wp), DIMENSION(nzb:nzt+1) ::  v_ref    !< reference profile for v-component
       REAL(wp), DIMENSION(nzb:nzt+1) ::  v_ref_l  !< reference profile for v-component on subdomain

       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ref_chem   !< reference profile for chemical species
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ref_chem_l !< reference profile for chemical species on subdomain

       IF ( debug_output_timestep )  CALL debug_message( 'nesting_offl_bc', 'start' )

       CALL  cpu_log( log_point(58), 'offline nesting', 'start' )      
!
!--    Initialize mean profiles, derived from boundary data, to zero
       pt_ref   = 0.0_wp
       q_ref    = 0.0_wp
       u_ref    = 0.0_wp
       v_ref    = 0.0_wp

       pt_ref_l = 0.0_wp
       q_ref_l  = 0.0_wp
       u_ref_l  = 0.0_wp
       v_ref_l  = 0.0_wp
!
!--    If required, allocate temporary arrays to compute chemistry mean profiles
       IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
          ALLOCATE( ref_chem(nzb:nzt+1,1:UBOUND( chem_species, 1 ) )   )
          ALLOCATE( ref_chem_l(nzb:nzt+1,1:UBOUND( chem_species, 1 ) ) )
          ref_chem   = 0.0_wp
          ref_chem_l = 0.0_wp
       ENDIF
!
!--    Set boundary conditions of u-, v-, w-component, as well as q, and pt.
!--    Note, boundary values at the left boundary: i=-1 (v,w,pt,q) and 
!--    i=0 (u), at the right boundary: i=nxr+1 (all), at the south boundary:
!--    j=-1 (u,w,pt,q) and j=0 (v), at the north boundary: j=nyn+1 (all).
!--    Please note, at the left (for u) and south (for v) boundary, values 
!--    for u and v are set also at i/j=-1, since these values are used in 
!--    boundary_conditions() to restore prognostic values.
!--    Further, sum up data to calculate mean profiles from boundary data, 
!--    used for Rayleigh damping. 
       IF ( bc_dirichlet_l )  THEN

          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                u(k,j,0) = interpolate_in_time( nest_offl%u_left(0,k,j),       &
                                                nest_offl%u_left(1,k,j),       &
                                                fac_dt ) *                     &
                             MERGE( 1.0_wp, 0.0_wp,                            &
                                    BTEST( wall_flags_total_0(k,j,0), 1 ) )
                u(k,j,-1) = u(k,j,0)
             ENDDO
             u_ref_l(nzb+1:nzt) = u_ref_l(nzb+1:nzt) + u(nzb+1:nzt,j,0)
          ENDDO

          DO  j = nys, nyn
             DO  k = nzb+1, nzt-1
                w(k,j,-1) = interpolate_in_time( nest_offl%w_left(0,k,j),      &
                                                 nest_offl%w_left(1,k,j),      &
                                                 fac_dt ) *                    &
                            MERGE( 1.0_wp, 0.0_wp,                             &
                                   BTEST( wall_flags_total_0(k,j,-1), 3 ) )
             ENDDO
             w(nzt,j,-1) = w(nzt-1,j,-1)
          ENDDO

          DO  j = nysv, nyn
             DO  k = nzb+1, nzt
                v(k,j,-1) = interpolate_in_time( nest_offl%v_left(0,k,j),      &
                                                 nest_offl%v_left(1,k,j),      &
                                                 fac_dt ) *                    &
                               MERGE( 1.0_wp, 0.0_wp,                          &
                                      BTEST( wall_flags_total_0(k,j,-1), 2 ) )
             ENDDO
             v_ref_l(nzb+1:nzt) = v_ref_l(nzb+1:nzt) + v(nzb+1:nzt,j,-1)
          ENDDO

          IF ( .NOT. neutral )  THEN
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   pt(k,j,-1) = interpolate_in_time( nest_offl%pt_left(0,k,j), &
                                                     nest_offl%pt_left(1,k,j), &
                                                     fac_dt )
 
                ENDDO
                pt_ref_l(nzb+1:nzt) = pt_ref_l(nzb+1:nzt) + pt(nzb+1:nzt,j,-1)
             ENDDO
          ENDIF

          IF ( humidity )  THEN
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   q(k,j,-1) = interpolate_in_time( nest_offl%q_left(0,k,j),   &
                                                    nest_offl%q_left(1,k,j),   &
                                                    fac_dt )
 
                ENDDO
                q_ref_l(nzb+1:nzt) = q_ref_l(nzb+1:nzt) + q(nzb+1:nzt,j,-1)
             ENDDO
          ENDIF

          IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
             DO  n = 1, UBOUND( chem_species, 1 )
                IF ( nest_offl%chem_from_file_l(n) )  THEN
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         chem_species(n)%conc(k,j,-1) = interpolate_in_time(   &
                                                  nest_offl%chem_left(0,k,j,n),&
                                                  nest_offl%chem_left(1,k,j,n),&
                                                  fac_dt                   )
                      ENDDO
                      ref_chem_l(nzb+1:nzt,n) = ref_chem_l(nzb+1:nzt,n)        &
                                         + chem_species(n)%conc(nzb+1:nzt,j,-1)
                   ENDDO
                ENDIF
             ENDDO
          ENDIF

       ENDIF

       IF ( bc_dirichlet_r )  THEN

          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                u(k,j,nxr+1) = interpolate_in_time( nest_offl%u_right(0,k,j),  &
                                                    nest_offl%u_right(1,k,j),  &
                                                    fac_dt ) *                 &
                             MERGE( 1.0_wp, 0.0_wp,                            &
                                    BTEST( wall_flags_total_0(k,j,nxr+1), 1 ) )
             ENDDO
             u_ref_l(nzb+1:nzt) = u_ref_l(nzb+1:nzt) + u(nzb+1:nzt,j,nxr+1)
          ENDDO
          DO  j = nys, nyn
             DO  k = nzb+1, nzt-1
                w(k,j,nxr+1) = interpolate_in_time( nest_offl%w_right(0,k,j),  &
                                                    nest_offl%w_right(1,k,j),  &
                                                    fac_dt ) *                 &
                             MERGE( 1.0_wp, 0.0_wp,                            &
                                    BTEST( wall_flags_total_0(k,j,nxr+1), 3 ) )
             ENDDO
             w(nzt,j,nxr+1) = w(nzt-1,j,nxr+1)
          ENDDO

          DO  j = nysv, nyn
             DO  k = nzb+1, nzt
                v(k,j,nxr+1) = interpolate_in_time( nest_offl%v_right(0,k,j),  &
                                                    nest_offl%v_right(1,k,j),  &
                                                    fac_dt ) *                 &
                             MERGE( 1.0_wp, 0.0_wp,                            &
                                    BTEST( wall_flags_total_0(k,j,nxr+1), 2 ) )
             ENDDO
             v_ref_l(nzb+1:nzt) = v_ref_l(nzb+1:nzt) + v(nzb+1:nzt,j,nxr+1)
          ENDDO

          IF ( .NOT. neutral )  THEN
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   pt(k,j,nxr+1) = interpolate_in_time(                        &
                                                  nest_offl%pt_right(0,k,j),   &
                                                  nest_offl%pt_right(1,k,j),   &
                                                  fac_dt )
                ENDDO
                pt_ref_l(nzb+1:nzt) = pt_ref_l(nzb+1:nzt) + pt(nzb+1:nzt,j,nxr+1)
             ENDDO
          ENDIF

          IF ( humidity )  THEN
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   q(k,j,nxr+1) = interpolate_in_time(                         &
                                                  nest_offl%q_right(0,k,j),    &
                                                  nest_offl%q_right(1,k,j),    &
                                                  fac_dt ) 
 
                ENDDO
                q_ref_l(nzb+1:nzt) = q_ref_l(nzb+1:nzt) + q(nzb+1:nzt,j,nxr+1)
             ENDDO
          ENDIF

          IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
             DO  n = 1, UBOUND( chem_species, 1 )
                IF ( nest_offl%chem_from_file_r(n) )  THEN
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         chem_species(n)%conc(k,j,nxr+1) = interpolate_in_time(&
                                                 nest_offl%chem_right(0,k,j,n),&
                                                 nest_offl%chem_right(1,k,j,n),&
                                                 fac_dt                       )
                      ENDDO
                      ref_chem_l(nzb+1:nzt,n) = ref_chem_l(nzb+1:nzt,n)        &
                                       + chem_species(n)%conc(nzb+1:nzt,j,nxr+1)
                   ENDDO
                ENDIF
             ENDDO
          ENDIF

       ENDIF

       IF ( bc_dirichlet_s )  THEN

          DO  i = nxl, nxr
             DO  k = nzb+1, nzt
                v(k,0,i) = interpolate_in_time( nest_offl%v_south(0,k,i),      &
                                                nest_offl%v_south(1,k,i),      &
                                                fac_dt ) *                     &
                           MERGE( 1.0_wp, 0.0_wp,                              &
                                  BTEST( wall_flags_total_0(k,0,i), 2 ) )
                v(k,-1,i) = v(k,0,i) 
             ENDDO
             v_ref_l(nzb+1:nzt) = v_ref_l(nzb+1:nzt) + v(nzb+1:nzt,0,i)
          ENDDO

          DO  i = nxl, nxr
             DO  k = nzb+1, nzt-1
                w(k,-1,i) = interpolate_in_time( nest_offl%w_south(0,k,i),     &
                                                 nest_offl%w_south(1,k,i),     &
                                                 fac_dt ) *                    &
                           MERGE( 1.0_wp, 0.0_wp,                              &
                                  BTEST( wall_flags_total_0(k,-1,i), 3 ) )
             ENDDO
             w(nzt,-1,i) = w(nzt-1,-1,i)
          ENDDO

          DO  i = nxlu, nxr
             DO  k = nzb+1, nzt
                u(k,-1,i) = interpolate_in_time( nest_offl%u_south(0,k,i),     &
                                                 nest_offl%u_south(1,k,i),     &
                                                 fac_dt ) *                    &
                           MERGE( 1.0_wp, 0.0_wp,                              &
                                  BTEST( wall_flags_total_0(k,-1,i), 1 ) )
             ENDDO
             u_ref_l(nzb+1:nzt) = u_ref_l(nzb+1:nzt) + u(nzb+1:nzt,-1,i)
          ENDDO

          IF ( .NOT. neutral )  THEN
             DO  i = nxl, nxr
                DO  k = nzb+1, nzt
                   pt(k,-1,i) = interpolate_in_time(                           &
                                                 nest_offl%pt_south(0,k,i),    &
                                                 nest_offl%pt_south(1,k,i),    &
                                                 fac_dt )
 
                ENDDO
                pt_ref_l(nzb+1:nzt) = pt_ref_l(nzb+1:nzt) + pt(nzb+1:nzt,-1,i)
             ENDDO
          ENDIF

          IF ( humidity )  THEN
             DO  i = nxl, nxr
                DO  k = nzb+1, nzt
                   q(k,-1,i) = interpolate_in_time(                            &
                                                 nest_offl%q_south(0,k,i),     &
                                                 nest_offl%q_south(1,k,i),     &
                                                 fac_dt )
 
                ENDDO
                q_ref_l(nzb+1:nzt) = q_ref_l(nzb+1:nzt) + q(nzb+1:nzt,-1,i)
             ENDDO
          ENDIF

          IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
             DO  n = 1, UBOUND( chem_species, 1 )
                IF ( nest_offl%chem_from_file_s(n) )  THEN
                   DO  i = nxl, nxr
                      DO  k = nzb+1, nzt
                         chem_species(n)%conc(k,-1,i) = interpolate_in_time(   &
                                                 nest_offl%chem_south(0,k,i,n),&
                                                 nest_offl%chem_south(1,k,i,n),&
                                                 fac_dt                    )
                      ENDDO
                      ref_chem_l(nzb+1:nzt,n) = ref_chem_l(nzb+1:nzt,n)        &
                                       + chem_species(n)%conc(nzb+1:nzt,-1,i)
                   ENDDO
                ENDIF
             ENDDO
          ENDIF

       ENDIF

       IF ( bc_dirichlet_n )  THEN

          DO  i = nxl, nxr
             DO  k = nzb+1, nzt
                v(k,nyn+1,i) = interpolate_in_time( nest_offl%v_north(0,k,i),  &
                                                    nest_offl%v_north(1,k,i),  &
                                                    fac_dt ) *                 &
                               MERGE( 1.0_wp, 0.0_wp,                          &
                                    BTEST( wall_flags_total_0(k,nyn+1,i), 2 ) )
             ENDDO
             v_ref_l(nzb+1:nzt) = v_ref_l(nzb+1:nzt) + v(nzb+1:nzt,nyn+1,i)
          ENDDO
          DO  i = nxl, nxr
             DO  k = nzb+1, nzt-1
                w(k,nyn+1,i) = interpolate_in_time( nest_offl%w_north(0,k,i),  &
                                                    nest_offl%w_north(1,k,i),  &
                                                    fac_dt ) *                 &
                               MERGE( 1.0_wp, 0.0_wp,                          &
                                    BTEST( wall_flags_total_0(k,nyn+1,i), 3 ) )
             ENDDO
             w(nzt,nyn+1,i) = w(nzt-1,nyn+1,i)
          ENDDO

          DO  i = nxlu, nxr
             DO  k = nzb+1, nzt
                u(k,nyn+1,i) = interpolate_in_time( nest_offl%u_north(0,k,i),  &
                                                    nest_offl%u_north(1,k,i),  &
                                                    fac_dt ) *                 &
                               MERGE( 1.0_wp, 0.0_wp,                          &
                                    BTEST( wall_flags_total_0(k,nyn+1,i), 1 ) )

             ENDDO
             u_ref_l(nzb+1:nzt) = u_ref_l(nzb+1:nzt) + u(nzb+1:nzt,nyn+1,i)
          ENDDO

          IF ( .NOT. neutral )  THEN
             DO  i = nxl, nxr
                DO  k = nzb+1, nzt
                   pt(k,nyn+1,i) = interpolate_in_time(                        &
                                                    nest_offl%pt_north(0,k,i), &
                                                    nest_offl%pt_north(1,k,i), &
                                                    fac_dt )
 
                ENDDO
                pt_ref_l(nzb+1:nzt) = pt_ref_l(nzb+1:nzt) + pt(nzb+1:nzt,nyn+1,i)
             ENDDO
          ENDIF

          IF ( humidity )  THEN
             DO  i = nxl, nxr
                DO  k = nzb+1, nzt
                   q(k,nyn+1,i) = interpolate_in_time(                         &
                                                    nest_offl%q_north(0,k,i),  &
                                                    nest_offl%q_north(1,k,i),  &
                                                    fac_dt )
 
                ENDDO
                q_ref_l(nzb+1:nzt) = q_ref_l(nzb+1:nzt) + q(nzb+1:nzt,nyn+1,i)
             ENDDO
          ENDIF

          IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
             DO  n = 1, UBOUND( chem_species, 1 )
                IF ( nest_offl%chem_from_file_n(n) )  THEN
                   DO  i = nxl, nxr
                      DO  k = nzb+1, nzt
                         chem_species(n)%conc(k,nyn+1,i) = interpolate_in_time(&
                                                 nest_offl%chem_north(0,k,i,n),&
                                                 nest_offl%chem_north(1,k,i,n),&
                                                 fac_dt                       )
                      ENDDO
                      ref_chem_l(nzb+1:nzt,n) = ref_chem_l(nzb+1:nzt,n)        &
                                       + chem_species(n)%conc(nzb+1:nzt,nyn+1,i)
                   ENDDO
                ENDIF
             ENDDO
          ENDIF

       ENDIF
!
!--    Top boundary
       DO  i = nxlu, nxr
          DO  j = nys, nyn
             u(nzt+1,j,i) = interpolate_in_time( nest_offl%u_top(0,j,i),       &
                                                 nest_offl%u_top(1,j,i),       &
                                                 fac_dt ) *                    &
                           MERGE( 1.0_wp, 0.0_wp,                              &
                                  BTEST( wall_flags_total_0(nzt+1,j,i), 1 ) )
             u_ref_l(nzt+1) = u_ref_l(nzt+1) + u(nzt+1,j,i)
          ENDDO
       ENDDO
!
!--    For left boundary set boundary condition for u-component also at top 
!--    grid point. 
!--    Note, this has no effect on the numeric solution, only for data output.
       IF ( bc_dirichlet_l )  u(nzt+1,:,nxl) = u(nzt+1,:,nxlu)

       DO  i = nxl, nxr
          DO  j = nysv, nyn
             v(nzt+1,j,i) = interpolate_in_time( nest_offl%v_top(0,j,i),       &
                                                 nest_offl%v_top(1,j,i),       &
                                                 fac_dt ) *                    &
                           MERGE( 1.0_wp, 0.0_wp,                              &
                                  BTEST( wall_flags_total_0(nzt+1,j,i), 2 ) )
             v_ref_l(nzt+1) = v_ref_l(nzt+1) + v(nzt+1,j,i)
          ENDDO
       ENDDO
!
!--    For south boundary set boundary condition for v-component also at top 
!--    grid point. 
!--    Note, this has no effect on the numeric solution, only for data output.
       IF ( bc_dirichlet_s )  v(nzt+1,nys,:) = v(nzt+1,nysv,:)

       DO  i = nxl, nxr
          DO  j = nys, nyn
             w(nzt,j,i) = interpolate_in_time( nest_offl%w_top(0,j,i),         &
                                               nest_offl%w_top(1,j,i),         &
                                               fac_dt ) *                      &
                           MERGE( 1.0_wp, 0.0_wp,                              &
                                  BTEST( wall_flags_total_0(nzt,j,i), 3 ) )
             w(nzt+1,j,i) = w(nzt,j,i)
          ENDDO
       ENDDO


       IF ( .NOT. neutral )  THEN
          DO  i = nxl, nxr
             DO  j = nys, nyn
                pt(nzt+1,j,i) = interpolate_in_time( nest_offl%pt_top(0,j,i),  &
                                                     nest_offl%pt_top(1,j,i),  &
                                                     fac_dt )
                pt_ref_l(nzt+1) = pt_ref_l(nzt+1) + pt(nzt+1,j,i)
             ENDDO
          ENDDO
       ENDIF

       IF ( humidity )  THEN
          DO  i = nxl, nxr
             DO  j = nys, nyn
                q(nzt+1,j,i) = interpolate_in_time( nest_offl%q_top(0,j,i),    &
                                                    nest_offl%q_top(1,j,i),    &
                                                    fac_dt )
                q_ref_l(nzt+1) = q_ref_l(nzt+1) + q(nzt+1,j,i)
             ENDDO
          ENDDO
       ENDIF

       IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
          DO  n = 1, UBOUND( chem_species, 1 )
             IF ( nest_offl%chem_from_file_t(n) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      chem_species(n)%conc(nzt+1,j,i) = interpolate_in_time(   &
                                              nest_offl%chem_top(0,j,i,n),     &
                                              nest_offl%chem_top(1,j,i,n),     &
                                              fac_dt                       )
                      ref_chem_l(nzt+1,n) = ref_chem_l(nzt+1,n) +              &
                                            chem_species(n)%conc(nzt+1,j,i)
                   ENDDO
                ENDDO
             ENDIF
          ENDDO
       ENDIF
!
!--    Moreover, set Neumann boundary condition for subgrid-scale TKE, 
!--    passive scalar, dissipation, and chemical species if required
       IF ( rans_mode  .AND.  rans_tke_e )  THEN
          IF (  bc_dirichlet_l )  diss(:,:,nxl-1) = diss(:,:,nxl)
          IF (  bc_dirichlet_r )  diss(:,:,nxr+1) = diss(:,:,nxr)
          IF (  bc_dirichlet_s )  diss(:,nys-1,:) = diss(:,nys,:)
          IF (  bc_dirichlet_n )  diss(:,nyn+1,:) = diss(:,nyn,:)
       ENDIF
!        IF ( .NOT. constant_diffusion )  THEN
!           IF (  bc_dirichlet_l )  e(:,:,nxl-1) = e(:,:,nxl)
!           IF (  bc_dirichlet_r )  e(:,:,nxr+1) = e(:,:,nxr)
!           IF (  bc_dirichlet_s )  e(:,nys-1,:) = e(:,nys,:)
!           IF (  bc_dirichlet_n )  e(:,nyn+1,:) = e(:,nyn,:)
!           e(nzt+1,:,:) = e(nzt,:,:)
!        ENDIF
!        IF ( passive_scalar )  THEN
!           IF (  bc_dirichlet_l )  s(:,:,nxl-1) = s(:,:,nxl)
!           IF (  bc_dirichlet_r )  s(:,:,nxr+1) = s(:,:,nxr)
!           IF (  bc_dirichlet_s )  s(:,nys-1,:) = s(:,nys,:)
!           IF (  bc_dirichlet_n )  s(:,nyn+1,:) = s(:,nyn,:)
!        ENDIF

       CALL exchange_horiz( u, nbgp )
       CALL exchange_horiz( v, nbgp )
       CALL exchange_horiz( w, nbgp )
       IF ( .NOT. neutral )  CALL exchange_horiz( pt, nbgp )
       IF ( humidity      )  CALL exchange_horiz( q,  nbgp )
       IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
          DO  n = 1, UBOUND( chem_species, 1 )
!
!--          Do local exchange only when necessary, i.e. when data is coming 
!--          from dynamic file.
             IF ( nest_offl%chem_from_file_t(n) )                              &
                CALL exchange_horiz( chem_species(n)%conc, nbgp )
          ENDDO
       ENDIF
!
!--    Set top boundary condition at all horizontal grid points, also at the
!--    lateral boundary grid points. 
       w(nzt+1,:,:) = w(nzt,:,:)       
!
!--    Offline nesting for salsa
       IF ( salsa )  CALL salsa_nesting_offl_bc
!
!--    In case of Rayleigh damping, where the profiles u_init, v_init
!--    q_init and pt_init are still used, update these profiles from the
!--    averaged boundary data. 
!--    But first, average these data.
#if defined( __parallel )
       CALL MPI_ALLREDUCE( u_ref_l, u_ref, nzt+1-nzb+1, MPI_REAL, MPI_SUM,     &
                           comm2d, ierr )
       CALL MPI_ALLREDUCE( v_ref_l, v_ref, nzt+1-nzb+1, MPI_REAL, MPI_SUM,     &
                           comm2d, ierr )
       IF ( humidity )  THEN
          CALL MPI_ALLREDUCE( q_ref_l, q_ref, nzt+1-nzb+1, MPI_REAL, MPI_SUM,  &
                              comm2d, ierr )
       ENDIF
       IF ( .NOT. neutral )  THEN
          CALL MPI_ALLREDUCE( pt_ref_l, pt_ref, nzt+1-nzb+1, MPI_REAL, MPI_SUM,&
                              comm2d, ierr )
       ENDIF
       IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
          CALL MPI_ALLREDUCE( ref_chem_l, ref_chem,                            &
                              ( nzt+1-nzb+1 ) * SIZE( ref_chem(nzb,:) ),       &
                              MPI_REAL, MPI_SUM, comm2d, ierr )
       ENDIF
#else
       u_ref  = u_ref_l
       v_ref  = v_ref_l
       IF ( humidity )       q_ref    = q_ref_l
       IF ( .NOT. neutral )  pt_ref   = pt_ref_l
       IF ( air_chemistry  .AND.  nesting_offline_chem )  ref_chem = ref_chem_l
#endif
!
!--    Average data. Note, reference profiles up to nzt are derived from lateral
!--    boundaries, at the model top it is derived from the top boundary. Thus,
!--    number of input data is different from nzb:nzt compared to nzt+1. 
!--    Derived from lateral boundaries. 
       u_ref(nzb:nzt) = u_ref(nzb:nzt) / REAL( 2.0_wp * ( ny + 1 + nx     ),   &
                                               KIND = wp )  
       v_ref(nzb:nzt) = v_ref(nzb:nzt) / REAL( 2.0_wp * ( ny   + nx + 1   ),   &
                                               KIND = wp )
       IF ( humidity )                                                         &
          q_ref(nzb:nzt) = q_ref(nzb:nzt)   / REAL( 2.0_wp *                   &
                                                          ( ny + 1 + nx + 1 ), &
                                                    KIND = wp )
       IF ( .NOT. neutral )                                                    &
          pt_ref(nzb:nzt) = pt_ref(nzb:nzt) / REAL( 2.0_wp *                   &
                                                          ( ny + 1 + nx + 1 ), &
                                              KIND = wp )
       IF ( air_chemistry  .AND.  nesting_offline_chem )                       &
          ref_chem(nzb:nzt,:) = ref_chem(nzb:nzt,:) / REAL( 2.0_wp *           &
                                                          ( ny + 1 + nx + 1 ), &
                                                            KIND = wp )
!
!--    Derived from top boundary.    
       u_ref(nzt+1) = u_ref(nzt+1) / REAL( ( ny + 1 ) * ( nx     ), KIND = wp ) 
       v_ref(nzt+1) = v_ref(nzt+1) / REAL( ( ny     ) * ( nx + 1 ), KIND = wp )
       IF ( humidity )                                                         &
          q_ref(nzt+1) = q_ref(nzt+1)   / REAL( ( ny + 1 ) * ( nx + 1 ),       &
                                                KIND = wp )
       IF ( .NOT. neutral )                                                    &
          pt_ref(nzt+1) = pt_ref(nzt+1) / REAL( ( ny + 1 ) * ( nx + 1 ),       &
                                                KIND = wp )
       IF ( air_chemistry  .AND.  nesting_offline_chem )                       &
          ref_chem(nzt+1,:) = ref_chem(nzt+1,:) /                              &
                              REAL( ( ny + 1 ) * ( nx + 1 ),KIND = wp )
!
!--    Write onto init profiles, which are used for damping. Also set lower
!--    boundary condition for scalars (not required for u and v as these are 
!--    zero at k=nzb.
       u_init = u_ref
       v_init = v_ref
       IF ( humidity      )  THEN 
          q_init      = q_ref
          q_init(nzb) = q_init(nzb+1)
       ENDIF
       IF ( .NOT. neutral )  THEN
          pt_init      = pt_ref
          pt_init(nzb) = pt_init(nzb+1)
       ENDIF

       IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
          DO  n = 1, UBOUND( chem_species, 1 )
             IF ( nest_offl%chem_from_file_t(n) )  THEN 
                chem_species(n)%conc_pr_init(:) = ref_chem(:,n)
                chem_species(n)%conc_pr_init(nzb) =                            &
                                            chem_species(n)%conc_pr_init(nzb+1)
             ENDIF
          ENDDO
       ENDIF

       IF ( ALLOCATED( ref_chem   ) )  DEALLOCATE( ref_chem   )
       IF ( ALLOCATED( ref_chem_l ) )  DEALLOCATE( ref_chem_l )
!
!--    Further, adjust Rayleigh damping height in case of time-changing conditions.
!--    Therefore, calculate boundary-layer depth first. 
       CALL nesting_offl_calc_zi
       CALL adjust_sponge_layer 

       CALL  cpu_log( log_point(58), 'offline nesting', 'stop' )

       IF ( debug_output_timestep )  CALL debug_message( 'nesting_offl_bc', 'end' )


    END SUBROUTINE nesting_offl_bc

!------------------------------------------------------------------------------!
! Description:
!------------------------------------------------------------------------------!
!>  Update of the geostrophic wind components.
!>  @todo: update geostrophic wind also in the child domains (should be done
!>         in the nesting.
!------------------------------------------------------------------------------!
    SUBROUTINE nesting_offl_geostrophic_wind

       INTEGER(iwp) ::  k
!
!--    Update geostrophic wind components from dynamic input file. 
       DO  k = nzb+1, nzt
          ug(k) = interpolate_in_time( nest_offl%ug(0,k), nest_offl%ug(1,k),   &
                                       fac_dt )
          vg(k) = interpolate_in_time( nest_offl%vg(0,k), nest_offl%vg(1,k),   &
                                       fac_dt )
       ENDDO
       ug(nzt+1) = ug(nzt)
       vg(nzt+1) = vg(nzt)

    END SUBROUTINE nesting_offl_geostrophic_wind

!------------------------------------------------------------------------------!
! Description:
!------------------------------------------------------------------------------!
!>  Determine the interpolation constant for time interpolation. The 
!>  calculation is separated from the nesting_offl_bc and 
!>  nesting_offl_geostrophic_wind in order to be independent on the order 
!>  of calls. 
!------------------------------------------------------------------------------!
    SUBROUTINE nesting_offl_interpolation_factor
!
!--    Determine interpolation factor and limit it to 1. This is because
!--    t+dt can slightly exceed time(tind_p) before boundary data is updated 
!--    again. 
       fac_dt = ( time_since_reference_point                                   &
                - nest_offl%time(nest_offl%tind) + dt_3d ) /                   &
           ( nest_offl%time(nest_offl%tind_p) - nest_offl%time(nest_offl%tind) )

       fac_dt = MIN( 1.0_wp, fac_dt )

    END SUBROUTINE nesting_offl_interpolation_factor

!------------------------------------------------------------------------------!
! Description:
!------------------------------------------------------------------------------!
!> Calculates the boundary-layer depth from the boundary data, according to 
!> bulk-Richardson criterion. 
!------------------------------------------------------------------------------!
    SUBROUTINE nesting_offl_calc_zi

       INTEGER(iwp) :: i                            !< loop index in x-direction
       INTEGER(iwp) :: j                            !< loop index in y-direction
       INTEGER(iwp) :: k                            !< loop index in z-direction
       INTEGER(iwp) :: k_max_loc                    !< index of maximum wind speed along z-direction
       INTEGER(iwp) :: k_surface                    !< topography top index in z-direction
       INTEGER(iwp) :: num_boundary_gp_non_cyclic   !< number of non-cyclic boundaries, used for averaging ABL depth
       INTEGER(iwp) :: num_boundary_gp_non_cyclic_l !< number of non-cyclic boundaries, used for averaging ABL depth
    
       REAL(wp) ::  ri_bulk                 !< bulk Richardson number
       REAL(wp) ::  ri_bulk_crit = 0.25_wp  !< critical bulk Richardson number
       REAL(wp) ::  topo_max                !< maximum topography level in model domain
       REAL(wp) ::  topo_max_l              !< maximum topography level in subdomain
       REAL(wp) ::  vpt_surface             !< near-surface virtual potential temperature
       REAL(wp) ::  zi_l                    !< mean boundary-layer depth on subdomain   
       REAL(wp) ::  zi_local                !< local boundary-layer depth  
    
       REAL(wp), DIMENSION(nzb:nzt+1) ::  vpt_col !< vertical profile of virtual potential temperature at (j,i)-grid point
       REAL(wp), DIMENSION(nzb:nzt+1) ::  uv_abs  !< vertical profile of horizontal wind speed at (j,i)-grid point

       
!
!--    Calculate mean boundary-layer height from boundary data.
!--    Start with the left and right boundaries. 
       zi_l      = 0.0_wp
       num_boundary_gp_non_cyclic_l = 0
       IF ( bc_dirichlet_l  .OR.  bc_dirichlet_r )  THEN
!
!--       Sum-up and store number of boundary grid points used for averaging
!--       ABL depth
          num_boundary_gp_non_cyclic_l = num_boundary_gp_non_cyclic_l +        &
                                         nxr - nxl + 1
!
!--       Determine index along x. Please note, index indicates boundary
!--       grid point for scalars. 
          i = MERGE( -1, nxr + 1, bc_dirichlet_l )
    
          DO  j = nys, nyn
!
!--          Determine topography top index at current (j,i) index
             k_surface = topo_top_ind(j,i,0)
!
!--          Pre-compute surface virtual temperature. Therefore, use 2nd 
!--          prognostic level according to Heinze et al. (2017).
             IF ( humidity )  THEN
                vpt_surface = pt(k_surface+2,j,i) *                            &
                            ( 1.0_wp + 0.61_wp * q(k_surface+2,j,i) )
                vpt_col     = pt(:,j,i) * ( 1.0_wp + 0.61_wp * q(:,j,i) )
             ELSE
                vpt_surface = pt(k_surface+2,j,i)
                vpt_col     = pt(:,j,i)
             ENDIF
!
!--          Calculate local boundary layer height from bulk Richardson number,
!--          i.e. the height where the bulk Richardson number exceeds its
!--          critical value of 0.25 (according to Heinze et al., 2017).
!--          Note, no interpolation of u- and v-component is made, as both 
!--          are mainly mean inflow profiles with very small spatial variation. 
!--          Add a safety factor in case the velocity term becomes zero. This
!--          may happen if overhanging 3D structures are directly located at
!--          the boundary, where velocity inside the building is zero 
!--          (k_surface is the index of the lowest upward-facing surface).
             uv_abs(:) = SQRT( MERGE( u(:,j,i+1), u(:,j,i),                   &
                                      bc_dirichlet_l )**2 +                   &
                               v(:,j,i)**2 )
!
!--          Determine index of the maximum wind speed                               
             k_max_loc = MAXLOC( uv_abs(:), DIM = 1 ) - 1

             zi_local = 0.0_wp
             DO  k = k_surface+1, nzt
                ri_bulk = zu(k) * g / vpt_surface *                            &
                          ( vpt_col(k) - vpt_surface ) /                       &
                          ( uv_abs(k) + 1E-5_wp ) 
!
!--             Check if critical Richardson number is exceeded. Further, check
!--             if there is a maxium in the wind profile in order to detect also
!--             ABL heights in the stable boundary layer. 
                IF ( zi_local == 0.0_wp  .AND.                                 &
                     ( ri_bulk > ri_bulk_crit  .OR.  k == k_max_loc ) )        &
                   zi_local = zu(k)           
             ENDDO
!
!--          Assure that the minimum local boundary-layer depth is at least at 
!--          the second vertical grid level.
             zi_l = zi_l + MAX( zi_local, zu(k_surface+2) )   
             
          ENDDO
       
       ENDIF
!
!--    Do the same at the north and south boundaries.
       IF ( bc_dirichlet_s  .OR.  bc_dirichlet_n )  THEN
       
          num_boundary_gp_non_cyclic_l = num_boundary_gp_non_cyclic_l +        &
                                         nxr - nxl + 1
       
          j = MERGE( -1, nyn + 1, bc_dirichlet_s )
       
          DO  i = nxl, nxr
             k_surface = topo_top_ind(j,i,0)
 
             IF ( humidity )  THEN
                vpt_surface = pt(k_surface+2,j,i) *                            &
                            ( 1.0_wp + 0.61_wp * q(k_surface+2,j,i) )
                vpt_col     = pt(:,j,i) * ( 1.0_wp + 0.61_wp * q(:,j,i) )
             ELSE
                vpt_surface = pt(k_surface+2,j,i)
                vpt_col  = pt(:,j,i)
             ENDIF
          
             uv_abs(:) = SQRT( u(:,j,i)**2 +                                   &
                               MERGE( v(:,j+1,i), v(:,j,i),                    &
                               bc_dirichlet_s )**2 )
!
!--          Determine index of the maximum wind speed
             k_max_loc = MAXLOC( uv_abs(:), DIM = 1 ) - 1
          
             zi_local = 0.0_wp
             DO  k = k_surface+1, nzt                
                ri_bulk = zu(k) * g / vpt_surface *                            &
                       ( vpt_col(k) - vpt_surface ) /                          &
                       ( uv_abs(k) + 1E-5_wp ) 
!
!--             Check if critical Richardson number is exceeded. Further, check
!--             if there is a maxium in the wind profile in order to detect also
!--             ABL heights in the stable boundary layer.                        
                IF ( zi_local == 0.0_wp  .AND.                                 &
                     ( ri_bulk > ri_bulk_crit  .OR.  k == k_max_loc ) )        &
                   zi_local = zu(k)    
             ENDDO
             zi_l = zi_l + MAX( zi_local, zu(k_surface+2) )   
          
          ENDDO
          
       ENDIF
    
#if defined( __parallel )
       CALL MPI_ALLREDUCE( zi_l, zi_ribulk, 1, MPI_REAL, MPI_SUM,              &
                           comm2d, ierr )
       CALL MPI_ALLREDUCE( num_boundary_gp_non_cyclic_l,                       &
                           num_boundary_gp_non_cyclic,                         &
                           1, MPI_INTEGER, MPI_SUM, comm2d, ierr )
#else
       zi_ribulk = zi_l
       num_boundary_gp_non_cyclic = num_boundary_gp_non_cyclic_l
#endif
       zi_ribulk = zi_ribulk / REAL( num_boundary_gp_non_cyclic, KIND = wp )
!
!--    Finally, check if boundary layer depth is not below the any topography. 
!--    zi_ribulk will be used to adjust rayleigh damping height, i.e. the 
!--    lower level of the sponge layer, as well as to adjust the synthetic
!--    turbulence generator accordingly. If Rayleigh damping would be applied
!--    near buildings, etc., this would spoil the simulation results. 
       topo_max_l = zw(MAXVAL( topo_top_ind(nys:nyn,nxl:nxr,0) ))
       
#if defined( __parallel )
       CALL MPI_ALLREDUCE( topo_max_l, topo_max, 1, MPI_REAL, MPI_MAX,         &
                           comm2d, ierr )
#else
       topo_max     = topo_max_l
#endif
!        zi_ribulk = MAX( zi_ribulk, topo_max )
       
    END SUBROUTINE nesting_offl_calc_zi
    
    
!------------------------------------------------------------------------------!
! Description:
!------------------------------------------------------------------------------!
!> Adjust the height where the rayleigh damping starts, i.e. the lower level
!> of the sponge layer. 
!------------------------------------------------------------------------------!
    SUBROUTINE adjust_sponge_layer

       INTEGER(iwp) :: k   !< loop index in z-direction
    
       REAL(wp) ::  rdh    !< updated Rayleigh damping height
 
    
       IF ( rayleigh_damping_height > 0.0_wp  .AND.                            &
            rayleigh_damping_factor > 0.0_wp )  THEN
!
!--       Update Rayleigh-damping height and re-calculate height-depending 
!--       damping coefficients. 
!--       Assure that rayleigh damping starts well above the boundary layer.   
          rdh = MIN( MAX( zi_ribulk * 1.3_wp, 10.0_wp * dz(1) ),               & 
                     0.8_wp * zu(nzt), rayleigh_damping_height )
!
!--       Update Rayleigh damping factor
          DO  k = nzb+1, nzt
             IF ( zu(k) >= rdh )  THEN
                rdf(k) = rayleigh_damping_factor *                             &
                                          ( SIN( pi * 0.5_wp * ( zu(k) - rdh ) &
                                        / ( zu(nzt) - rdh ) )                  &
                                          )**2
             ENDIF
          ENDDO
          rdf_sc = rdf
       
       ENDIF

    END SUBROUTINE adjust_sponge_layer
    
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Performs consistency checks
!------------------------------------------------------------------------------!
    SUBROUTINE nesting_offl_check_parameters 
!
!--    Check if offline nesting is applied in nested child domain.
       IF ( nesting_offline  .AND.  child_domain )  THEN
          message_string = 'Offline nesting is only applicable in root model.'
          CALL message( 'offline_nesting_check_parameters', 'PA0622', 1, 2, 0, 6, 0 )       
       ENDIF

    END SUBROUTINE nesting_offl_check_parameters
    
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads the parameter list nesting_offl_parameters
!------------------------------------------------------------------------------!
    SUBROUTINE nesting_offl_parin 
       
       CHARACTER (LEN=80) ::  line   !< dummy string that contains the current line of the parameter file


       NAMELIST /nesting_offl_parameters/   nesting_offline

       line = ' '

!
!--    Try to find stg package
       REWIND ( 11 )
       line = ' '
       DO WHILE ( INDEX( line, '&nesting_offl_parameters' ) == 0 )
          READ ( 11, '(A)', END=20 )  line
       ENDDO
       BACKSPACE ( 11 )

!
!--    Read namelist
       READ ( 11, nesting_offl_parameters, ERR = 10, END = 20 )

       GOTO 20

 10    BACKSPACE( 11 )
       READ( 11 , '(A)') line
       CALL parin_fail_message( 'nesting_offl_parameters', line )

 20    CONTINUE


    END SUBROUTINE nesting_offl_parin

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Writes information about offline nesting into HEADER file
!------------------------------------------------------------------------------!
    SUBROUTINE nesting_offl_header ( io )

       INTEGER(iwp), INTENT(IN) ::  io !< Unit of the output file

       WRITE ( io, 1 )
       IF ( nesting_offline )  THEN
          WRITE ( io, 3 )
       ELSE
          WRITE ( io, 2 )
       ENDIF

1 FORMAT (//' Offline nesting in COSMO model:'/                                &
              ' -------------------------------'/)
2 FORMAT (' --> No offlince nesting is used (default) ')
3 FORMAT (' --> Offlince nesting is used. Boundary data is read from dynamic input file ')

    END SUBROUTINE nesting_offl_header 

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate arrays used to read boundary data from NetCDF file and initialize
!> boundary data. 
!------------------------------------------------------------------------------!
    SUBROUTINE nesting_offl_init
           
       INTEGER(iwp) ::  n !< running index for chemical species

!--    Allocate arrays for geostrophic wind components. Arrays will 
!--    incorporate 2 time levels in order to interpolate in between. 
       ALLOCATE( nest_offl%ug(0:1,1:nzt) )
       ALLOCATE( nest_offl%vg(0:1,1:nzt) )
!
!--    Allocate arrays for reading left/right boundary values. Arrays will 
!--    incorporate 2  time levels in order to interpolate in between. If the core has 
!--    no boundary, allocate a dummy array, in order to enable netcdf parallel
!--    access. Dummy arrays will be allocated with dimension length zero.
       IF ( bc_dirichlet_l )  THEN
          ALLOCATE( nest_offl%u_left(0:1,nzb+1:nzt,nys:nyn)  )
          ALLOCATE( nest_offl%v_left(0:1,nzb+1:nzt,nysv:nyn) )
          ALLOCATE( nest_offl%w_left(0:1,nzb+1:nzt-1,nys:nyn) )
          IF ( humidity )       ALLOCATE( nest_offl%q_left(0:1,nzb+1:nzt,nys:nyn)  )
          IF ( .NOT. neutral )  ALLOCATE( nest_offl%pt_left(0:1,nzb+1:nzt,nys:nyn) )
          IF ( air_chemistry  .AND.  nesting_offline_chem )                                        &
             ALLOCATE( nest_offl%chem_left(0:1,nzb+1:nzt,nys:nyn,1:UBOUND( chem_species, 1 )) )
       ELSE
          ALLOCATE( nest_offl%u_left(1:1,1:1,1:1)  )
          ALLOCATE( nest_offl%v_left(1:1,1:1,1:1)  )
          ALLOCATE( nest_offl%w_left(1:1,1:1,1:1)  )
          IF ( humidity )       ALLOCATE( nest_offl%q_left(1:1,1:1,1:1)  )
          IF ( .NOT. neutral )  ALLOCATE( nest_offl%pt_left(1:1,1:1,1:1)  )
          IF ( air_chemistry  .AND.  nesting_offline_chem )                                        &
             ALLOCATE( nest_offl%chem_left(1:1,1:1,1:1,1:UBOUND( chem_species, 1 )) )
       ENDIF
       IF ( bc_dirichlet_r )  THEN
          ALLOCATE( nest_offl%u_right(0:1,nzb+1:nzt,nys:nyn)  )
          ALLOCATE( nest_offl%v_right(0:1,nzb+1:nzt,nysv:nyn) )
          ALLOCATE( nest_offl%w_right(0:1,nzb+1:nzt-1,nys:nyn) )
          IF ( humidity )       ALLOCATE( nest_offl%q_right(0:1,nzb+1:nzt,nys:nyn)  )
          IF ( .NOT. neutral )  ALLOCATE( nest_offl%pt_right(0:1,nzb+1:nzt,nys:nyn) )
          IF ( air_chemistry  .AND.  nesting_offline_chem )                                        &
             ALLOCATE( nest_offl%chem_right(0:1,nzb+1:nzt,nys:nyn,1:UBOUND( chem_species, 1 )) )
       ELSE
          ALLOCATE( nest_offl%u_right(1:1,1:1,1:1)  )
          ALLOCATE( nest_offl%v_right(1:1,1:1,1:1)  )
          ALLOCATE( nest_offl%w_right(1:1,1:1,1:1)  )
          IF ( humidity )       ALLOCATE( nest_offl%q_right(1:1,1:1,1:1)  )
          IF ( .NOT. neutral )  ALLOCATE( nest_offl%pt_right(1:1,1:1,1:1)  )
          IF ( air_chemistry  .AND.  nesting_offline_chem )                                        &
             ALLOCATE( nest_offl%chem_right(1:1,1:1,1:1,1:UBOUND( chem_species, 1 )) )
       ENDIF
!
!--    Allocate arrays for reading north/south boundary values. Arrays will 
!--    incorporate 2  time levels in order to interpolate in between. If the core has 
!--    no boundary, allocate a dummy array, in order to enable netcdf parallel
!--    access. Dummy arrays will be allocated with dimension length zero.
       IF ( bc_dirichlet_n )  THEN
          ALLOCATE( nest_offl%u_north(0:1,nzb+1:nzt,nxlu:nxr) )
          ALLOCATE( nest_offl%v_north(0:1,nzb+1:nzt,nxl:nxr)  )
          ALLOCATE( nest_offl%w_north(0:1,nzb+1:nzt-1,nxl:nxr) )
          IF ( humidity )       ALLOCATE( nest_offl%q_north(0:1,nzb+1:nzt,nxl:nxr)  )
          IF ( .NOT. neutral )  ALLOCATE( nest_offl%pt_north(0:1,nzb+1:nzt,nxl:nxr) )
          IF ( air_chemistry  .AND.  nesting_offline_chem )                                        &
             ALLOCATE( nest_offl%chem_north(0:1,nzb+1:nzt,nxl:nxr,1:UBOUND( chem_species, 1 )) )
       ELSE
          ALLOCATE( nest_offl%u_north(1:1,1:1,1:1)  )
          ALLOCATE( nest_offl%v_north(1:1,1:1,1:1)  )
          ALLOCATE( nest_offl%w_north(1:1,1:1,1:1)  )
          IF ( humidity )       ALLOCATE( nest_offl%q_north(1:1,1:1,1:1)  )
          IF ( .NOT. neutral )  ALLOCATE( nest_offl%pt_north(1:1,1:1,1:1)  )
          IF ( air_chemistry  .AND.  nesting_offline_chem )                                        &
             ALLOCATE( nest_offl%chem_north(1:1,1:1,1:1,1:UBOUND( chem_species, 1 )) )
       ENDIF
       IF ( bc_dirichlet_s )  THEN
          ALLOCATE( nest_offl%u_south(0:1,nzb+1:nzt,nxlu:nxr) )
          ALLOCATE( nest_offl%v_south(0:1,nzb+1:nzt,nxl:nxr)  )
          ALLOCATE( nest_offl%w_south(0:1,nzb+1:nzt-1,nxl:nxr)    )
          IF ( humidity )       ALLOCATE( nest_offl%q_south(0:1,nzb+1:nzt,nxl:nxr)  )
          IF ( .NOT. neutral )  ALLOCATE( nest_offl%pt_south(0:1,nzb+1:nzt,nxl:nxr) )
          IF ( air_chemistry  .AND.  nesting_offline_chem )                                        &
             ALLOCATE( nest_offl%chem_south(0:1,nzb+1:nzt,nxl:nxr,1:UBOUND( chem_species, 1 )) )
       ELSE
          ALLOCATE( nest_offl%u_south(1:1,1:1,1:1)  )
          ALLOCATE( nest_offl%v_south(1:1,1:1,1:1)  )
          ALLOCATE( nest_offl%w_south(1:1,1:1,1:1)  )
          IF ( humidity )       ALLOCATE( nest_offl%q_south(1:1,1:1,1:1)  )
          IF ( .NOT. neutral )  ALLOCATE( nest_offl%pt_south(1:1,1:1,1:1)  )
          IF ( air_chemistry  .AND.  nesting_offline_chem )                                        &
             ALLOCATE( nest_offl%chem_south(1:1,1:1,1:1,1:UBOUND( chem_species, 1 )) )
       ENDIF
!
!--    Allocate arrays for reading data at the top boundary. In contrast to the
!--    lateral boundaries, every core reads these data so that no dummy
!--    arrays need to be allocated.
       ALLOCATE( nest_offl%u_top(0:1,nys:nyn,nxlu:nxr) )
       ALLOCATE( nest_offl%v_top(0:1,nysv:nyn,nxl:nxr) )
       ALLOCATE( nest_offl%w_top(0:1,nys:nyn,nxl:nxr)  )
       IF ( humidity )       ALLOCATE( nest_offl%q_top(0:1,nys:nyn,nxl:nxr)  )
       IF ( .NOT. neutral )  ALLOCATE( nest_offl%pt_top(0:1,nys:nyn,nxl:nxr) )
       IF ( air_chemistry  .AND.  nesting_offline_chem )                                           &
          ALLOCATE( nest_offl%chem_top(0:1,nys:nyn,nxl:nxr,1:UBOUND( chem_species, 1 )) )
!
!--    For chemical species, create the names of the variables. This is necessary
!--    to identify the respective variable and write it onto the correct array 
!--    in the chem_species datatype.
       IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
          ALLOCATE( nest_offl%chem_from_file_l(1:UBOUND( chem_species, 1 )) )
          ALLOCATE( nest_offl%chem_from_file_n(1:UBOUND( chem_species, 1 )) )
          ALLOCATE( nest_offl%chem_from_file_r(1:UBOUND( chem_species, 1 )) )
          ALLOCATE( nest_offl%chem_from_file_s(1:UBOUND( chem_species, 1 )) )
          ALLOCATE( nest_offl%chem_from_file_t(1:UBOUND( chem_species, 1 )) )

          ALLOCATE( nest_offl%var_names_chem_l(1:UBOUND( chem_species, 1 )) )
          ALLOCATE( nest_offl%var_names_chem_n(1:UBOUND( chem_species, 1 )) )
          ALLOCATE( nest_offl%var_names_chem_r(1:UBOUND( chem_species, 1 )) )
          ALLOCATE( nest_offl%var_names_chem_s(1:UBOUND( chem_species, 1 )) )
          ALLOCATE( nest_offl%var_names_chem_t(1:UBOUND( chem_species, 1 )) )
!
!--       Initialize flags that indicate whether the variable is on file or
!--       not. Please note, this is only necessary for chemistry variables. 
          nest_offl%chem_from_file_l(:) = .FALSE.
          nest_offl%chem_from_file_n(:) = .FALSE.
          nest_offl%chem_from_file_r(:) = .FALSE.
          nest_offl%chem_from_file_s(:) = .FALSE.
          nest_offl%chem_from_file_t(:) = .FALSE.

          DO  n = 1, UBOUND( chem_species, 1 )
             nest_offl%var_names_chem_l(n) = nest_offl%char_l //               &
                                                  TRIM(chem_species(n)%name)
             nest_offl%var_names_chem_n(n) = nest_offl%char_n //               &
                                                  TRIM(chem_species(n)%name)
             nest_offl%var_names_chem_r(n) = nest_offl%char_r //               &
                                                  TRIM(chem_species(n)%name)
             nest_offl%var_names_chem_s(n) = nest_offl%char_s //               &
                                                  TRIM(chem_species(n)%name)
             nest_offl%var_names_chem_t(n) = nest_offl%char_t //               &
                                                  TRIM(chem_species(n)%name)
          ENDDO
       ENDIF
!
!--    Offline nesting for salsa
       IF ( salsa )  CALL salsa_nesting_offl_init
!
!--    Before initial data input is initiated, check if dynamic input file is 
!--    present.
       IF ( .NOT. input_pids_dynamic )  THEN
          message_string = 'nesting_offline = .TRUE. requires dynamic '  //    &
                            'input file ' //                                   &
                            TRIM( input_file_dynamic ) // TRIM( coupling_char )
          CALL message( 'nesting_offl_init', 'PA0546', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Read COSMO data at lateral and top boundaries
       CALL nesting_offl_input
!
!--    Check if sufficient time steps are provided to cover the entire 
!--    simulation. Note, dynamic input is only required for the 3D simulation,
!--    not for the soil/wall spinup. However, as the spinup time is added
!--    to the end_time, this must be considered here. 
       IF ( end_time - spinup_time > nest_offl%time(nest_offl%nt-1) )  THEN
          message_string = 'end_time of the simulation exceeds the ' //        &
                           'time dimension in the dynamic input file.'
          CALL message( 'nesting_offl_init', 'PA0183', 1, 2, 0, 6, 0 ) 
       ENDIF
!
!--    Initialize boundary data. Please note, do not initialize boundaries in 
!--    case of restart runs. This case the boundaries are already initialized
!--    and the boundary data from file would be on the wrong time level.  
       IF ( TRIM( initializing_actions ) /= 'read_restart_data' )  THEN
          IF ( bc_dirichlet_l )  THEN
             u(nzb+1:nzt,nys:nyn,0)    = nest_offl%u_left(0,nzb+1:nzt,nys:nyn)
             v(nzb+1:nzt,nysv:nyn,-1)  = nest_offl%v_left(0,nzb+1:nzt,nysv:nyn) 
             w(nzb+1:nzt-1,nys:nyn,-1) = nest_offl%w_left(0,nzb+1:nzt-1,nys:nyn)
             IF ( .NOT. neutral )  pt(nzb+1:nzt,nys:nyn,-1) =                  &
                                      nest_offl%pt_left(0,nzb+1:nzt,nys:nyn)
             IF ( humidity      )  q(nzb+1:nzt,nys:nyn,-1)  =                  &
                                      nest_offl%q_left(0,nzb+1:nzt,nys:nyn)
             IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
                DO  n = 1, UBOUND( chem_species, 1 )
                   IF( nest_offl%chem_from_file_l(n) )  THEN
                      chem_species(n)%conc(nzb+1:nzt,nys:nyn,-1) =             &
                                      nest_offl%chem_left(0,nzb+1:nzt,nys:nyn,n)
                   ENDIF
                ENDDO
             ENDIF
          ENDIF
          IF ( bc_dirichlet_r )  THEN
             u(nzb+1:nzt,nys:nyn,nxr+1)   = nest_offl%u_right(0,nzb+1:nzt,nys:nyn) 
             v(nzb+1:nzt,nysv:nyn,nxr+1)  = nest_offl%v_right(0,nzb+1:nzt,nysv:nyn)
             w(nzb+1:nzt-1,nys:nyn,nxr+1) = nest_offl%w_right(0,nzb+1:nzt-1,nys:nyn)
             IF ( .NOT. neutral )  pt(nzb+1:nzt,nys:nyn,nxr+1) =               &
                                      nest_offl%pt_right(0,nzb+1:nzt,nys:nyn)
             IF ( humidity      )  q(nzb+1:nzt,nys:nyn,nxr+1)  =               &
                                      nest_offl%q_right(0,nzb+1:nzt,nys:nyn)
             IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
                DO  n = 1, UBOUND( chem_species, 1 )
                   IF( nest_offl%chem_from_file_r(n) )  THEN
                      chem_species(n)%conc(nzb+1:nzt,nys:nyn,nxr+1) =          &
                                      nest_offl%chem_right(0,nzb+1:nzt,nys:nyn,n)
                   ENDIF
                ENDDO
             ENDIF
          ENDIF
          IF ( bc_dirichlet_s )  THEN
             u(nzb+1:nzt,-1,nxlu:nxr)  = nest_offl%u_south(0,nzb+1:nzt,nxlu:nxr)
             v(nzb+1:nzt,0,nxl:nxr)    = nest_offl%v_south(0,nzb+1:nzt,nxl:nxr) 
             w(nzb+1:nzt-1,-1,nxl:nxr) = nest_offl%w_south(0,nzb+1:nzt-1,nxl:nxr) 
             IF ( .NOT. neutral )  pt(nzb+1:nzt,-1,nxl:nxr) =                  &
                                      nest_offl%pt_south(0,nzb+1:nzt,nxl:nxr)
             IF ( humidity      )  q(nzb+1:nzt,-1,nxl:nxr)  =                  &
                                      nest_offl%q_south(0,nzb+1:nzt,nxl:nxr) 
             IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
                DO  n = 1, UBOUND( chem_species, 1 )
                   IF( nest_offl%chem_from_file_s(n) )  THEN
                      chem_species(n)%conc(nzb+1:nzt,-1,nxl:nxr) =             &
                                      nest_offl%chem_south(0,nzb+1:nzt,nxl:nxr,n)
                   ENDIF
                ENDDO
             ENDIF
          ENDIF
          IF ( bc_dirichlet_n )  THEN
             u(nzb+1:nzt,nyn+1,nxlu:nxr)  = nest_offl%u_north(0,nzb+1:nzt,nxlu:nxr)
             v(nzb+1:nzt,nyn+1,nxl:nxr)   = nest_offl%v_north(0,nzb+1:nzt,nxl:nxr) 
             w(nzb+1:nzt-1,nyn+1,nxl:nxr) = nest_offl%w_north(0,nzb+1:nzt-1,nxl:nxr)
             IF ( .NOT. neutral )  pt(nzb+1:nzt,nyn+1,nxl:nxr) =               &
                                      nest_offl%pt_north(0,nzb+1:nzt,nxl:nxr) 
             IF ( humidity      )  q(nzb+1:nzt,nyn+1,nxl:nxr)  =               &
                                      nest_offl%q_north(0,nzb+1:nzt,nxl:nxr)
             IF ( air_chemistry  .AND.  nesting_offline_chem )  THEN
                DO  n = 1, UBOUND( chem_species, 1 )
                   IF( nest_offl%chem_from_file_n(n) )  THEN
                      chem_species(n)%conc(nzb+1:nzt,nyn+1,nxl:nxr) =          &
                                      nest_offl%chem_north(0,nzb+1:nzt,nxl:nxr,n)
                   ENDIF
                ENDDO
             ENDIF
          ENDIF
!
!--       Initialize geostrophic wind components. Actually this is already done in
!--       init_3d_model when initializing_action = 'inifor', however, in speical 
!--       case of user-defined initialization this will be done here again, in 
!--       order to have a consistent initialization. 
          ug(nzb+1:nzt) = nest_offl%ug(0,nzb+1:nzt)
          vg(nzb+1:nzt) = nest_offl%vg(0,nzb+1:nzt)
!
!--       Set bottom and top boundary condition for geostrophic wind components
          ug(nzt+1) = ug(nzt)
          vg(nzt+1) = vg(nzt)
          ug(nzb)   = ug(nzb+1)
          vg(nzb)   = vg(nzb+1)
       ENDIF
!
!--    After boundary data is initialized, mask topography at the 
!--    boundaries for the velocity components.
       u = MERGE( u, 0.0_wp, BTEST( wall_flags_total_0, 1 ) )
       v = MERGE( v, 0.0_wp, BTEST( wall_flags_total_0, 2 ) )
       w = MERGE( w, 0.0_wp, BTEST( wall_flags_total_0, 3 ) )
!
!--    Initial calculation of the boundary layer depth from the prescribed 
!--    boundary data. This is requiered for initialize the synthetic turbulence
!--    generator correctly. 
       CALL nesting_offl_calc_zi
       
!
!--    After boundary data is initialized, ensure mass conservation. Not 
!--    necessary in restart runs.
       IF ( TRIM( initializing_actions ) /= 'read_restart_data' )  THEN
          CALL nesting_offl_mass_conservation
       ENDIF

    END SUBROUTINE nesting_offl_init
    
!------------------------------------------------------------------------------!
! Description:
!------------------------------------------------------------------------------!
!> Interpolation function, used to interpolate boundary data in time. 
!------------------------------------------------------------------------------!
    FUNCTION interpolate_in_time( var_t1, var_t2, fac  ) 

       REAL(wp)            :: interpolate_in_time !< time-interpolated boundary value
       REAL(wp)            :: var_t1              !< boundary value at t1
       REAL(wp)            :: var_t2              !< boundary value at t2
       REAL(wp)            :: fac                 !< interpolation factor

       interpolate_in_time = ( 1.0_wp - fac ) * var_t1 + fac * var_t2      

    END FUNCTION interpolate_in_time



 END MODULE nesting_offl_mod
