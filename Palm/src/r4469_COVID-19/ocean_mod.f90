!> @file ocean_mod.f90
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
! Copyright 2017-2018 Leibniz Universitaet Hannover
!--------------------------------------------------------------------------------!
!
! Current revisions:
! -----------------
! 
! 
! Former revisions:
! -----------------
! $Id: ocean_mod.f90 4370 2020-01-10 14:00:44Z raasch $
! vector directives added to force vectorization on Intel19 compiler
! 
! 4346 2019-12-18 11:55:56Z motisi
! Introduction of wall_flags_total_0, which currently sets bits based on static
! topography information used in wall_flags_static_0
! 
! 4329 2019-12-10 15:46:36Z motisi
! Renamed wall_flags_0 to wall_flags_static_0
! 
! 4272 2019-10-23 15:18:57Z schwenkel
! Further modularization of boundary conditions: moved boundary conditions to
! respective modules
!
! 4196 2019-08-29 11:02:06Z gronemeier
! Consider rotation of model domain for calculating the Stokes drift
! 
! 4182 2019-08-22 15:20:23Z scharf
! Corrected "Former revisions" section
! 
! 4110 2019-07-22 17:05:21Z suehring
! Pass integer flag array as well as boundary flags to WS scalar advection 
! routine
! 
! 4109 2019-07-22 17:00:34Z suehring
! implemented ocean_actions
! 
! 3767 2019-02-27 08:18:02Z raasch
! unused variable for file index and tmp_2d removed from rrd-subroutine parameter
! list
! 
! 3719 2019-02-06 13:10:18Z kanani
! Changed log_point to log_point_s, otherwise this overlaps with 
! 'all progn.equations' cpu measurement.
! 
! 3684 2019-01-20 20:20:58Z knoop
! nopointer option removed
! 3294 2018-10-01 02:37:10Z raasch
! initial revision
!
!
! Authors:
! --------
! @author Siegfried Raasch
!
! Description:
! ------------
!> This module contains relevant code for PALM's ocean mode, e.g. the
!> prognostic equation for salinity, the equation of state for seawater,
!> the Craik Leibovich force (Stokes force), and wave breaking effects
!------------------------------------------------------------------------------!
 MODULE ocean_mod
 

    USE arrays_3d,                                                             &
        ONLY:  prho, prho_1, rho_ocean, rho_1, sa, sa_init, sa_1, sa_2, sa_3,  &
               sa_p, tsa_m, flux_l_sa, flux_s_sa, diss_l_sa, diss_s_sa

    USE control_parameters,                                                    &
        ONLY:  atmos_ocean_sign, bottom_salinityflux,                          &
               constant_top_salinityflux, ocean_mode, top_salinityflux,        &
               wall_salinityflux, loop_optimization, ws_scheme_sca

    USE kinds

    USE pegrid,                                                                &
        ONLY:  threads_per_task

    USE statistics,                                                            &
        ONLY:  sums_wssas_ws_l

    USE indices,                                                               &
        ONLY:  advc_flags_s, nxl, nxr, nyn, nys, nzb, nzt, wall_flags_total_0

    USE surface_mod,                                                           &
        ONLY:  bc_h, surf_def_v, surf_def_h, surf_lsm_h, surf_lsm_v,           &
               surf_usm_h, surf_usm_v

    IMPLICIT NONE

    CHARACTER (LEN=20) ::  bc_sa_t = 'neumann'  !< namelist parameter

    INTEGER(iwp) ::  ibc_sa_t   !< integer flag for bc_sa_t
    INTEGER(iwp) ::  iran_ocean = -1234567  !< random number used for wave breaking effects

    INTEGER(iwp) ::  sa_vertical_gradient_level_ind(10) = -9999  !< grid index values of sa_vertical_gradient_level(s)

    LOGICAL ::  salinity = .TRUE.             !< switch for using salinity
    LOGICAL ::  stokes_force = .FALSE.        !< switch to switch on the Stokes force
    LOGICAL ::  wave_breaking = .FALSE.       !< switch to switch on wave breaking effects
    LOGICAL ::  surface_cooling_switched_off = .FALSE.  !< variable to check if surface heat flux has been switched off

    REAL(wp) ::  alpha_wave_breaking = 3.0_wp !< coefficient for wave breaking generated turbulence from Noh et al. (2004), JPO
    REAL(wp) ::  prho_reference               !< reference state of potential density at ocean surface
    REAL(wp) ::  sa_surface = 35.0_wp         !< surface salinity, namelist parameter
    REAL(wp) ::  sa_vertical_gradient(10) = 0.0_wp               !< namelist parameter
    REAL(wp) ::  sa_vertical_gradient_level(10) = -999999.9_wp   !< namelist parameter
    REAL(wp) ::  stokes_waveheight = 0.0_wp  !< wave height assumed for Stokes drift velocity
    REAL(wp) ::  stokes_wavelength = 0.0_wp  !< wavelength assumed for Stokes drift velocity
    REAL(wp) ::  surface_cooling_spinup_time = 999999.9_wp  !< time after which surface heat flux is switched off
    REAL(wp) ::  timescale_wave_breaking     !< time scale of random forcing
    REAL(wp) ::  u_star_wave_breaking        !< to store the absolute value of friction velocity at the ocean surface

    REAL(wp), DIMENSION(12), PARAMETER ::  nom =                               &
                          (/ 9.9984085444849347D2,   7.3471625860981584D0,     &
                            -5.3211231792841769D-2,  3.6492439109814549D-4,    &
                             2.5880571023991390D0,  -6.7168282786692354D-3,    &
                             1.9203202055760151D-3,  1.1798263740430364D-2,    &
                             9.8920219266399117D-8,  4.6996642771754730D-6,    &
                            -2.5862187075154352D-8, -3.2921414007960662D-12 /)
                          !< constants used in equation of state for seawater

    REAL(wp), DIMENSION(13), PARAMETER ::  den =                               &
                          (/ 1.0D0,                  7.2815210113327091D-3,    &
                            -4.4787265461983921D-5,  3.3851002965802430D-7,    &
                             1.3651202389758572D-10, 1.7632126669040377D-3,    &
                            -8.8066583251206474D-6, -1.8832689434804897D-10,   &
                             5.7463776745432097D-6,  1.4716275472242334D-9,    &
                             6.7103246285651894D-6, -2.4461698007024582D-17,   &
                            -9.1534417604289062D-18 /)
                          !< constants used in equation of state for seawater

    SAVE

    PUBLIC ::  bc_sa_t, ibc_sa_t, prho_reference, sa_surface,                  &
               sa_vertical_gradient, sa_vertical_gradient_level,               &
               sa_vertical_gradient_level_ind, stokes_force, wave_breaking


    INTERFACE eqn_state_seawater
       MODULE PROCEDURE eqn_state_seawater
       MODULE PROCEDURE eqn_state_seawater_ij
    END INTERFACE eqn_state_seawater

    INTERFACE eqn_state_seawater_func
       MODULE PROCEDURE eqn_state_seawater_func
    END INTERFACE eqn_state_seawater_func

    INTERFACE ocean_check_parameters
       MODULE PROCEDURE ocean_check_parameters
    END INTERFACE ocean_check_parameters

    INTERFACE ocean_check_data_output
       MODULE PROCEDURE ocean_check_data_output
    END INTERFACE ocean_check_data_output

    INTERFACE ocean_check_data_output_pr
       MODULE PROCEDURE ocean_check_data_output_pr
    END INTERFACE ocean_check_data_output_pr

    INTERFACE ocean_define_netcdf_grid
       MODULE PROCEDURE ocean_define_netcdf_grid
    END INTERFACE ocean_define_netcdf_grid

    INTERFACE ocean_data_output_2d
       MODULE PROCEDURE ocean_data_output_2d
    END INTERFACE ocean_data_output_2d

    INTERFACE ocean_data_output_3d
       MODULE PROCEDURE ocean_data_output_3d
    END INTERFACE ocean_data_output_3d

    INTERFACE ocean_header
       MODULE PROCEDURE ocean_header
    END INTERFACE ocean_header

    INTERFACE ocean_init
       MODULE PROCEDURE ocean_init
    END INTERFACE ocean_init

    INTERFACE ocean_init_arrays
       MODULE PROCEDURE ocean_init_arrays
    END INTERFACE ocean_init_arrays

    INTERFACE ocean_parin
       MODULE PROCEDURE ocean_parin
    END INTERFACE ocean_parin

    INTERFACE ocean_actions
       MODULE PROCEDURE ocean_actions
       MODULE PROCEDURE ocean_actions_ij
    END INTERFACE ocean_actions

    INTERFACE ocean_prognostic_equations
       MODULE PROCEDURE ocean_prognostic_equations
       MODULE PROCEDURE ocean_prognostic_equations_ij
    END INTERFACE ocean_prognostic_equations

    INTERFACE ocean_boundary_conditions
       MODULE PROCEDURE ocean_boundary_conditions
    END INTERFACE ocean_boundary_conditions

    INTERFACE ocean_swap_timelevel
       MODULE PROCEDURE ocean_swap_timelevel
    END INTERFACE ocean_swap_timelevel

    INTERFACE ocean_rrd_global
       MODULE PROCEDURE ocean_rrd_global
    END INTERFACE ocean_rrd_global

    INTERFACE ocean_rrd_local
       MODULE PROCEDURE ocean_rrd_local
    END INTERFACE ocean_rrd_local

    INTERFACE ocean_wrd_global
       MODULE PROCEDURE ocean_wrd_global
    END INTERFACE ocean_wrd_global

    INTERFACE ocean_wrd_local
       MODULE PROCEDURE ocean_wrd_local
    END INTERFACE ocean_wrd_local

    INTERFACE ocean_3d_data_averaging
       MODULE PROCEDURE ocean_3d_data_averaging
    END INTERFACE ocean_3d_data_averaging

    INTERFACE stokes_drift_terms
       MODULE PROCEDURE stokes_drift_terms
       MODULE PROCEDURE stokes_drift_terms_ij
    END INTERFACE stokes_drift_terms

    INTERFACE wave_breaking_term
       MODULE PROCEDURE wave_breaking_term
       MODULE PROCEDURE wave_breaking_term_ij
    END INTERFACE wave_breaking_term

    PRIVATE
!
!-- Add INTERFACES that must be available to other modules (alphabetical order)
    PUBLIC eqn_state_seawater, ocean_actions, ocean_check_data_output,         &
           ocean_check_data_output_pr, ocean_check_parameters,                 &
           ocean_data_output_2d, ocean_data_output_3d,                         &
           ocean_define_netcdf_grid, ocean_header, ocean_init,                 &
           ocean_init_arrays, ocean_parin, ocean_prognostic_equations,         &
           ocean_swap_timelevel, ocean_rrd_global, ocean_rrd_local,            &
           ocean_wrd_global, ocean_wrd_local, ocean_3d_data_averaging,         &
           ocean_boundary_conditions, stokes_drift_terms, wave_breaking_term


 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Equation of state for seawater as a function of potential temperature,
!> salinity, and pressure.
!> For coefficients see Jackett et al., 2006: J. Atm. Ocean Tech.
!> eqn_state_seawater calculates the potential density referred at hyp(0).
!> eqn_state_seawater_func calculates density.
!> TODO: so far, routine is not adjusted to use topography
!------------------------------------------------------------------------------!
 SUBROUTINE eqn_state_seawater

    USE arrays_3d,                                                             &
        ONLY:  hyp, prho, pt_p, rho_ocean, sa_p
    USE indices,                                                               &
        ONLY:  nxl, nxr, nyn, nys, nzb, nzt

    IMPLICIT NONE

    INTEGER(iwp) ::  i       !< running index x direction
    INTEGER(iwp) ::  j       !< running index y direction
    INTEGER(iwp) ::  k       !< running index z direction
    INTEGER(iwp) ::  m       !< running index surface elements

    REAL(wp) ::  pden   !< temporary scalar user for calculating density
    REAL(wp) ::  pnom   !< temporary scalar user for calculating density
    REAL(wp) ::  p1     !< temporary scalar user for calculating density
    REAL(wp) ::  p2     !< temporary scalar user for calculating density
    REAL(wp) ::  p3     !< temporary scalar user for calculating density
    REAL(wp) ::  pt1    !< temporary scalar user for calculating density
    REAL(wp) ::  pt2    !< temporary scalar user for calculating density
    REAL(wp) ::  pt3    !< temporary scalar user for calculating density
    REAL(wp) ::  pt4    !< temporary scalar user for calculating density
    REAL(wp) ::  sa1    !< temporary scalar user for calculating density
    REAL(wp) ::  sa15   !< temporary scalar user for calculating density
    REAL(wp) ::  sa2    !< temporary scalar user for calculating density


    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt
!
!--          Pressure is needed in dbar
             p1 = hyp(k) * 1E-4_wp
             p2 = p1 * p1
             p3 = p2 * p1

!
!--          Temperature needed in degree Celsius
             pt1 = pt_p(k,j,i) - 273.15_wp
             pt2 = pt1 * pt1
             pt3 = pt1 * pt2
             pt4 = pt2 * pt2

             sa1  = sa_p(k,j,i)
             sa15 = sa1 * SQRT( sa1 )
             sa2  = sa1 * sa1

             pnom = nom(1)           + nom(2)*pt1     + nom(3)*pt2     +       &
                    nom(4)*pt3       + nom(5)*sa1     + nom(6)*sa1*pt1 +       &
                    nom(7)*sa2

             pden = den(1)           + den(2)*pt1     + den(3)*pt2     +       &
                    den(4)*pt3       + den(5)*pt4     + den(6)*sa1     +       &
                    den(7)*sa1*pt1   + den(8)*sa1*pt3 + den(9)*sa15    +       &
                    den(10)*sa15*pt2
!
!--          Potential density (without pressure terms)
             prho(k,j,i) = pnom / pden

             pnom = pnom +             nom(8)*p1      + nom(9)*p1*pt2  +       &
                    nom(10)*p1*sa1   + nom(11)*p2     + nom(12)*p2*pt2

             pden = pden +             den(11)*p1     + den(12)*p2*pt3 +       &
                    den(13)*p3*pt1

!
!--          In-situ density
             rho_ocean(k,j,i) = pnom / pden

          ENDDO

!
!--       Neumann conditions are assumed at top boundary and currently also at
!--       bottom boundary (see comment lines below)
          prho(nzt+1,j,i)      = prho(nzt,j,i)
          rho_ocean(nzt+1,j,i) = rho_ocean(nzt,j,i)

       ENDDO
    ENDDO
!
!-- Neumann conditions at up/downward-facing surfaces
    !$OMP PARALLEL DO PRIVATE( i, j, k )
    DO  m = 1, bc_h(0)%ns
       i = bc_h(0)%i(m)
       j = bc_h(0)%j(m)
       k = bc_h(0)%k(m)
       prho(k-1,j,i)      = prho(k,j,i)
       rho_ocean(k-1,j,i) = rho_ocean(k,j,i)
    ENDDO
!
!-- Downward facing surfaces
    !$OMP PARALLEL DO PRIVATE( i, j, k )
    DO  m = 1, bc_h(1)%ns
       i = bc_h(1)%i(m)
       j = bc_h(1)%j(m)
       k = bc_h(1)%k(m)
       prho(k+1,j,i)      = prho(k,j,i)
       rho_ocean(k+1,j,i) = rho_ocean(k,j,i)
    ENDDO

 END SUBROUTINE eqn_state_seawater


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Same as above, but for grid point i,j
!------------------------------------------------------------------------------!
 SUBROUTINE eqn_state_seawater_ij( i, j )

    USE arrays_3d,                                                             &
        ONLY:  hyp, prho, pt_p, rho_ocean, sa_p

    USE indices,                                                               &
        ONLY:  nzb, nzt

    IMPLICIT NONE

    INTEGER(iwp) ::  i       !< running index x direction
    INTEGER(iwp) ::  j       !< running index y direction
    INTEGER(iwp) ::  k       !< running index z direction
    INTEGER(iwp) ::  m       !< running index surface elements
    INTEGER(iwp) ::  surf_e  !< end index of surface elements at (j,i)-gridpoint
    INTEGER(iwp) ::  surf_s  !< start index of surface elements at (j,i)-gridpoint

    REAL(wp) ::  pden   !< temporary scalar user for calculating density
    REAL(wp) ::  pnom   !< temporary scalar user for calculating density
    REAL(wp) ::  p1     !< temporary scalar user for calculating density
    REAL(wp) ::  p2     !< temporary scalar user for calculating density
    REAL(wp) ::  p3     !< temporary scalar user for calculating density
    REAL(wp) ::  pt1    !< temporary scalar user for calculating density
    REAL(wp) ::  pt2    !< temporary scalar user for calculating density
    REAL(wp) ::  pt3    !< temporary scalar user for calculating density
    REAL(wp) ::  pt4    !< temporary scalar user for calculating density
    REAL(wp) ::  sa1    !< temporary scalar user for calculating density
    REAL(wp) ::  sa15   !< temporary scalar user for calculating density
    REAL(wp) ::  sa2    !< temporary scalar user for calculating density

    DO  k = nzb+1, nzt
!
!--    Pressure is needed in dbar
       p1 = hyp(k) * 1E-4_wp
       p2 = p1 * p1
       p3 = p2 * p1
!
!--    Temperature needed in degree Celsius
       pt1 = pt_p(k,j,i) - 273.15_wp
       pt2 = pt1 * pt1
       pt3 = pt1 * pt2
       pt4 = pt2 * pt2

       sa1  = sa_p(k,j,i)
       sa15 = sa1 * SQRT( sa1 )
       sa2  = sa1 * sa1

       pnom = nom(1)           + nom(2)*pt1     + nom(3)*pt2     +             &
              nom(4)*pt3       + nom(5)*sa1     + nom(6)*sa1*pt1 + nom(7)*sa2

       pden = den(1)           + den(2)*pt1     + den(3)*pt2     +             &
              den(4)*pt3       + den(5)*pt4     + den(6)*sa1     +             &
              den(7)*sa1*pt1   + den(8)*sa1*pt3 + den(9)*sa15    +             &
              den(10)*sa15*pt2
!
!--    Potential density (without pressure terms)
       prho(k,j,i) = pnom / pden

       pnom = pnom +             nom(8)*p1      + nom(9)*p1*pt2  +             &
              nom(10)*p1*sa1   + nom(11)*p2     + nom(12)*p2*pt2
       pden = pden +             den(11)*p1     + den(12)*p2*pt3 +             &
              den(13)*p3*pt1

!
!--    In-situ density
       rho_ocean(k,j,i) = pnom / pden

    ENDDO
!
!-- Neumann conditions at up/downward-facing walls
    surf_s = bc_h(0)%start_index(j,i)
    surf_e = bc_h(0)%end_index(j,i)
    DO  m = surf_s, surf_e
       k = bc_h(0)%k(m)
       prho(k-1,j,i)      = prho(k,j,i)
       rho_ocean(k-1,j,i) = rho_ocean(k,j,i)
    ENDDO
!
!-- Downward facing surfaces
    surf_s = bc_h(1)%start_index(j,i)
    surf_e = bc_h(1)%end_index(j,i)
    DO  m = surf_s, surf_e
       k = bc_h(1)%k(m)
       prho(k+1,j,i)      = prho(k,j,i)
       rho_ocean(k+1,j,i) = rho_ocean(k,j,i)
    ENDDO
!
!-- Neumann condition are assumed at top boundary
    prho(nzt+1,j,i)      = prho(nzt,j,i)
    rho_ocean(nzt+1,j,i) = rho_ocean(nzt,j,i)

 END SUBROUTINE eqn_state_seawater_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Equation of state for seawater as a function
!------------------------------------------------------------------------------!
 REAL(wp) FUNCTION eqn_state_seawater_func( p, pt, sa )

    IMPLICIT NONE

    REAL(wp) ::  p      !< temporary scalar user for calculating density
    REAL(wp) ::  p1     !< temporary scalar user for calculating density
    REAL(wp) ::  p2     !< temporary scalar user for calculating density
    REAL(wp) ::  p3     !< temporary scalar user for calculating density
    REAL(wp) ::  pt     !< temporary scalar user for calculating density
    REAL(wp) ::  pt1    !< temporary scalar user for calculating density
    REAL(wp) ::  pt2    !< temporary scalar user for calculating density
    REAL(wp) ::  pt3    !< temporary scalar user for calculating density
    REAL(wp) ::  pt4    !< temporary scalar user for calculating density
    REAL(wp) ::  sa     !< temporary scalar user for calculating density
    REAL(wp) ::  sa15   !< temporary scalar user for calculating density
    REAL(wp) ::  sa2    !< temporary scalar user for calculating density

!
!-- Pressure is needed in dbar
    p1 = p  * 1.0E-4_wp
    p2 = p1 * p1
    p3 = p2 * p1

!
!-- Temperature needed in degree Celsius
    pt1 = pt - 273.15_wp
    pt2 = pt1 * pt1
    pt3 = pt1 * pt2
    pt4 = pt2 * pt2

    sa15 = sa * SQRT( sa )
    sa2  = sa * sa


    eqn_state_seawater_func =                                                  &
         ( nom(1)        + nom(2)*pt1       + nom(3)*pt2    + nom(4)*pt3     + &
           nom(5)*sa     + nom(6)*sa*pt1    + nom(7)*sa2    + nom(8)*p1      + &
           nom(9)*p1*pt2 + nom(10)*p1*sa    + nom(11)*p2    + nom(12)*p2*pt2   &
         ) /                                                                   &
         ( den(1)        + den(2)*pt1       + den(3)*pt2    + den(4)*pt3     + &
           den(5)*pt4    + den(6)*sa        + den(7)*sa*pt1 + den(8)*sa*pt3  + &
           den(9)*sa15   + den(10)*sa15*pt2 + den(11)*p1    + den(12)*p2*pt3 + &
           den(13)*p3*pt1                                                      &
         )


 END FUNCTION eqn_state_seawater_func


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads the ocean parameters namelist
!------------------------------------------------------------------------------!
 SUBROUTINE ocean_parin

    IMPLICIT NONE

    CHARACTER (LEN=80) ::  line  !< dummy string that contains the current line of the parameter file


    NAMELIST /ocean_parameters/  bc_sa_t, bottom_salinityflux, salinity,       &
             sa_surface, sa_vertical_gradient, sa_vertical_gradient_level,     &
             stokes_waveheight, stokes_wavelength, surface_cooling_spinup_time,&
             top_salinityflux, wall_salinityflux, wave_breaking

!
!-- Try to find the namelist
    REWIND ( 11 )
    line = ' '
    DO WHILE ( INDEX( line, '&ocean_parameters' ) == 0 )
       READ ( 11, '(A)', END=12 )  line
    ENDDO
    BACKSPACE ( 11 )

!
!-- Read namelist
    READ ( 11, ocean_parameters, ERR = 10 )
!
!-- Set switch that enables PALM's ocean mode
    ocean_mode = .TRUE.

    GOTO 12

 10 BACKSPACE( 11 )
    READ( 11 , '(A)') line
    CALL parin_fail_message( 'ocean_parameters', line )

 12 CONTINUE

 END SUBROUTINE ocean_parin

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check parameters routine for the ocean mode
!------------------------------------------------------------------------------!
 SUBROUTINE ocean_check_parameters

    USE control_parameters,                                                    &
        ONLY:  coupling_char, coupling_mode, initializing_actions,             &
               message_string, use_top_fluxes

    USE pmc_interface,                                                         &
        ONLY:  nested_run

    IMPLICIT NONE


!
!-- Check for invalid combinations
    IF ( nested_run )  THEN
       message_string = 'ocean mode not allowed for nesting'
       CALL message( 'ocean_check_parameters', 'PA0510', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( TRIM( initializing_actions ) == 'cyclic_fill' )  THEN
       message_string = 'ocean mode does not allow cyclic-fill initialization'
       CALL message( 'ocean_check_parameters', 'PA0511', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Check ocean setting
    IF ( TRIM( coupling_mode ) == 'uncoupled'  .AND.                           &
         TRIM( coupling_char ) == '_O' .AND.                                   &
         .NOT. ocean_mode )  THEN

!
!--    Check whether an (uncoupled) atmospheric run has been declared as an
!--    ocean run (this setting is done via palmrun-option -y)
       message_string = 'ocean mode does not allow coupling_char = "' //       &
                        TRIM( coupling_char ) // '" set by palmrun-option "-y"'
       CALL message( 'ocean_check_parameters', 'PA0317', 1, 2, 0, 6, 0 )

    ENDIF

!
!-- Ocean version must use flux boundary conditions at the top
    IF ( .NOT. use_top_fluxes )  THEN
       message_string = 'use_top_fluxes must be .TRUE. in ocean mode'
       CALL message( 'ocean_check_parameters', 'PA0042', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Boundary conditions for salinity
    IF ( bc_sa_t == 'dirichlet' )  THEN
       ibc_sa_t = 0
    ELSEIF ( bc_sa_t == 'neumann' )  THEN
       ibc_sa_t = 1
    ELSE
       message_string = 'unknown boundary condition: bc_sa_t = "' //           &
                        TRIM( bc_sa_t ) // '"'
       CALL message( 'ocean_check_parameters', 'PA0068', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( top_salinityflux == 9999999.9_wp )  constant_top_salinityflux = .FALSE.

    IF ( .NOT. salinity )  THEN
       IF ( ( bottom_salinityflux /= 0.0_wp  .AND.                             &
              bottom_salinityflux /= 9999999.9_wp )  .OR.                      &
            ( top_salinityflux /= 0.0_wp     .AND.                             &
              top_salinityflux /= 9999999.9_wp ) )                             &
       THEN
          message_string = 'salinityflux must not be set for ocean run ' //    &
                           'without salinity'
          CALL message( 'ocean_check_parameters', 'PA0509', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

    IF ( ibc_sa_t == 1  .AND.  top_salinityflux == 9999999.9_wp )  THEN
       message_string = 'boundary condition: bc_sa_t = "' //                   &
                        TRIM( bc_sa_t ) // '" requires to set top_salinityflux'
       CALL message( 'ocean_check_parameters', 'PA0069', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- A fixed salinity at the top implies Dirichlet boundary condition for
!-- salinity. In this case specification of a constant salinity flux is
!-- forbidden.
    IF ( ibc_sa_t == 0  .AND.  constant_top_salinityflux  .AND.                &
         top_salinityflux /= 0.0_wp )  THEN
       message_string = 'boundary condition: bc_sa_t = "' //                   &
                        TRIM( bc_sa_t ) // '" is not allowed with ' //         &
                        'top_salinityflux /= 0.0'
       CALL message( 'ocean_check_parameters', 'PA0070', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Check if Stokes force is to be used
    IF ( stokes_waveheight /= 0.0_wp  .AND.  stokes_wavelength /= 0.0_wp )  THEN
       stokes_force = .TRUE.
    ELSE
       IF ( ( stokes_waveheight <= 0.0_wp .AND. stokes_wavelength > 0.0_wp ) &
            .OR.                                                               &
            ( stokes_waveheight > 0.0_wp .AND. stokes_wavelength <= 0.0_wp ) &
            .OR.                                                               &
            ( stokes_waveheight < 0.0_wp .AND. stokes_wavelength < 0.0_wp  ) ) &
       THEN
          message_string = 'wrong settings for stokes_wavelength and/or ' //   &
                           'stokes_waveheight'
          CALL message( 'ocean_check_parameters', 'PA0460', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

 END SUBROUTINE ocean_check_parameters


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check data output.
!------------------------------------------------------------------------------!
 SUBROUTINE ocean_check_data_output( var, unit )
 
    IMPLICIT NONE

    CHARACTER (LEN=*) ::  unit     !< unit of output variable
    CHARACTER (LEN=*) ::  var      !< name of output variable


    SELECT CASE ( TRIM( var ) )

       CASE ( 'rho_sea_water' )
          unit = 'kg/m3'

       CASE ( 'sa' )
          unit = 'psu'

       CASE DEFAULT
          unit = 'illegal'

    END SELECT

 END SUBROUTINE ocean_check_data_output


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check data output of profiles
!------------------------------------------------------------------------------!
 SUBROUTINE ocean_check_data_output_pr( variable, var_count, unit, dopr_unit )

    USE arrays_3d,                                                             &
        ONLY:  zu, zw

    USE control_parameters,                                                    &
        ONLY:  data_output_pr

    USE indices

    USE profil_parameter

    USE statistics

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  unit      !<
    CHARACTER (LEN=*) ::  variable  !<
    CHARACTER (LEN=*) ::  dopr_unit !< local value of dopr_unit

    INTEGER(iwp) ::  var_count     !<

    SELECT CASE ( TRIM( variable ) )

       CASE ( 'prho' )
          dopr_index(var_count) = 71
          dopr_unit             = 'kg/m3'
          hom(:,2,71,:) = SPREAD( zu, 2, statistic_regions+1 )
          unit = dopr_unit

       CASE ( 'rho_sea_water' )
          dopr_index(var_count) = 64
          dopr_unit             = 'kg/m3'
          hom(:,2,64,:) = SPREAD( zu, 2, statistic_regions+1 )
          IF ( data_output_pr(var_count)(1:1) == '#' )  THEN
             dopr_initial_index(var_count) = 77
             hom(:,2,77,:)             = SPREAD( zu, 2, statistic_regions+1 )
             hom(nzb,2,77,:)           = 0.0_wp    ! because zu(nzb) is negative
             data_output_pr(var_count) = data_output_pr(var_count)(2:)
          ENDIF
          unit = dopr_unit

       CASE ( 'sa', '#sa' )
          dopr_index(var_count) = 23
          dopr_unit             = 'psu'
          hom(:,2,23,:) = SPREAD( zu, 2, statistic_regions+1 )
          IF ( data_output_pr(var_count)(1:1) == '#' )  THEN
             dopr_initial_index(var_count) = 26
             hom(:,2,26,:)             = SPREAD( zu, 2, statistic_regions+1 )
             hom(nzb,2,26,:)           = 0.0_wp    ! because zu(nzb) is negative
             data_output_pr(var_count) = data_output_pr(var_count)(2:)
          ENDIF
          unit = dopr_unit

       CASE ( 'w"sa"' )
          dopr_index(var_count) = 65
          dopr_unit             = 'psu m/s'
          hom(:,2,65,:) = SPREAD( zw, 2, statistic_regions+1 )
          unit = dopr_unit

       CASE ( 'w*sa*' )
          dopr_index(var_count) = 66
          dopr_unit             = 'psu m/s'
          hom(:,2,66,:) = SPREAD( zw, 2, statistic_regions+1 )
          unit = dopr_unit

       CASE ( 'wsa' )
          dopr_index(var_count) = 67
          dopr_unit             = 'psu m/s'
          hom(:,2,67,:) = SPREAD( zw, 2, statistic_regions+1 )
          unit = dopr_unit

       CASE DEFAULT
          unit = 'illegal'

    END SELECT


 END SUBROUTINE ocean_check_data_output_pr


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Define appropriate grid for netcdf variables.
!> It is called out from subroutine netcdf.
!------------------------------------------------------------------------------!
 SUBROUTINE ocean_define_netcdf_grid( var, found, grid_x, grid_y, grid_z )
    
    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT(OUT) ::  grid_x   !< x grid of output variable
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_y   !< y grid of output variable
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_z   !< z grid of output variable
    CHARACTER (LEN=*), INTENT(IN)  ::  var      !< name of output variable

    LOGICAL, INTENT(OUT) ::  found   !< flag if output variable is found

    found  = .TRUE.

!
!-- Check for the grid
    SELECT CASE ( TRIM( var ) )

       CASE ( 'rho_sea_water', 'rho_sea_water_xy', 'rho_sea_water_xz', &
              'rho_sea_water_yz', 'sa', 'sa_xy', 'sa_xz', 'sa_yz' )
          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zu'

       CASE DEFAULT
          found  = .FALSE.
          grid_x = 'none'
          grid_y = 'none'
          grid_z = 'none'

    END SELECT

 END SUBROUTINE ocean_define_netcdf_grid


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Average 3D data.
!------------------------------------------------------------------------------!
 SUBROUTINE ocean_3d_data_averaging( mode, variable )
 

    USE arrays_3d,                                                             &
        ONLY:  rho_ocean, sa

    USE averaging,                                                             &
        ONLY:  rho_ocean_av, sa_av

    USE control_parameters,                                                    &
        ONLY:  average_count_3d

    USE indices,                                                               &
        ONLY:  nxlg, nxrg, nyng, nysg, nzb, nzt

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  mode       !< flag defining mode 'allocate', 'sum' or 'average'
    CHARACTER (LEN=*) ::  variable   !< name of variable

    INTEGER(iwp) ::  i   !< loop index
    INTEGER(iwp) ::  j   !< loop index
    INTEGER(iwp) ::  k   !< loop index

    IF ( mode == 'allocate' )  THEN

       SELECT CASE ( TRIM( variable ) )

          CASE ( 'rho_sea_water' )
             IF ( .NOT. ALLOCATED( rho_ocean_av ) )  THEN
                ALLOCATE( rho_ocean_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             rho_ocean_av = 0.0_wp

          CASE ( 'sa' )
             IF ( .NOT. ALLOCATED( sa_av ) )  THEN
                ALLOCATE( sa_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             sa_av = 0.0_wp

          CASE DEFAULT
             CONTINUE

       END SELECT

    ELSEIF ( mode == 'sum' )  THEN

       SELECT CASE ( TRIM( variable ) )

          CASE ( 'rho_sea_water' )
             IF ( ALLOCATED( rho_ocean_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rho_ocean_av(k,j,i) = rho_ocean_av(k,j,i) +           &
                                               rho_ocean(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'sa' )
             IF ( ALLOCATED( sa_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         sa_av(k,j,i) = sa_av(k,j,i) + sa(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE DEFAULT
             CONTINUE

       END SELECT

    ELSEIF ( mode == 'average' )  THEN

       SELECT CASE ( TRIM( variable ) )

          CASE ( 'rho_sea_water' )
             IF ( ALLOCATED( rho_ocean_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rho_ocean_av(k,j,i) = rho_ocean_av(k,j,i) /           &
                                               REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'sa' )
             IF ( ALLOCATED( sa_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         sa_av(k,j,i) = sa_av(k,j,i) /                         &
                                        REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

       END SELECT

    ENDIF

 END SUBROUTINE ocean_3d_data_averaging


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Define 2D output variables.
!------------------------------------------------------------------------------!
 SUBROUTINE ocean_data_output_2d( av, variable, found, grid, mode, local_pf,   &
                                  nzb_do, nzt_do )
 
    USE arrays_3d,                                                             &
        ONLY:  rho_ocean, sa

    USE averaging,                                                             &
        ONLY:  rho_ocean_av, sa_av

    USE indices,                                                               &
        ONLY: nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg, nzb, nzt,            &
              wall_flags_total_0

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  grid       !< name of vertical grid
    CHARACTER (LEN=*) ::  mode       !< either 'xy', 'xz' or 'yz'
    CHARACTER (LEN=*) ::  variable   !< name of variable

    INTEGER(iwp) ::  av        !< flag for (non-)average output
    INTEGER(iwp) ::  flag_nr   !< number of masking flag
    INTEGER(iwp) ::  i         !< loop index
    INTEGER(iwp) ::  j         !< loop index
    INTEGER(iwp) ::  k         !< loop index
    INTEGER(iwp) ::  nzb_do    !< vertical output index (bottom)
    INTEGER(iwp) ::  nzt_do    !< vertical output index (top)

    LOGICAL ::  found   !< flag if output variable is found
    LOGICAL ::  resorted  !< flag if output is already resorted

    REAL(wp) ::  fill_value = -999.0_wp  !< value for the _FillValue attribute

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf !< local
       !< array to which output data is resorted to

    REAL(wp), DIMENSION(:,:,:), POINTER ::  to_be_resorted  !< points to selected output variable
    
    found = .TRUE.
    resorted = .FALSE.
!
!-- Set masking flag for topography for not resorted arrays
    flag_nr = 0

    SELECT CASE ( TRIM( variable ) )

       CASE ( 'rho_sea_water_xy', 'rho_sea_water_xz', 'rho_sea_water_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => rho_ocean
          ELSE
             IF ( .NOT. ALLOCATED( rho_ocean_av ) ) THEN
                ALLOCATE( rho_ocean_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                rho_ocean_av = REAL( fill_value, KIND = wp )
             ENDIF
             to_be_resorted => rho_ocean_av
          ENDIF

       CASE ( 'sa_xy', 'sa_xz', 'sa_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => sa
          ELSE
             IF ( .NOT. ALLOCATED( sa_av ) ) THEN
                ALLOCATE( sa_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                sa_av = REAL( fill_value, KIND = wp )
             ENDIF
             to_be_resorted => sa_av
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
                local_pf(i,j,k) = MERGE( to_be_resorted(k,j,i),                &
                                   REAL( fill_value, KIND = wp ),              &
                                   BTEST( wall_flags_total_0(k,j,i), flag_nr ) )
             ENDDO
          ENDDO
       ENDDO
       resorted = .TRUE.
    ENDIF
 
 END SUBROUTINE ocean_data_output_2d

 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Define 3D output variables.
!------------------------------------------------------------------------------!
 SUBROUTINE ocean_data_output_3d( av, variable, found, local_pf, nzb_do, nzt_do )
 

    USE arrays_3d,                                                             &
        ONLY:  rho_ocean, sa

    USE averaging,                                                             &
        ONLY:  rho_ocean_av, sa_av

    USE indices,                                                               &
        ONLY:  nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg, nzb, nzt,           &
               wall_flags_total_0

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  variable   !< name of variable

    INTEGER(iwp) ::  av        !< flag for (non-)average output
    INTEGER(iwp) ::  flag_nr   !< number of masking flag
    INTEGER(iwp) ::  i         !< loop index
    INTEGER(iwp) ::  j         !< loop index
    INTEGER(iwp) ::  k         !< loop index
    INTEGER(iwp) ::  nzb_do    !< lower limit of the data output (usually 0)
    INTEGER(iwp) ::  nzt_do    !< vertical upper limit of the data output (usually nz_do3d)

    LOGICAL ::  found     !< flag if output variable is found
    LOGICAL ::  resorted  !< flag if output is already resorted

    REAL(wp) ::  fill_value = -999.0_wp  !< value for the _FillValue attribute

    REAL(sp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf   !< local
                                  !< array to which output data is resorted to

    REAL(wp), DIMENSION(:,:,:), POINTER ::  to_be_resorted  !< points to selected output variable

    found = .TRUE.
    resorted = .FALSE.
!
!-- Set masking flag for topography for not resorted arrays
    flag_nr = 0

    SELECT CASE ( TRIM( variable ) )

       CASE ( 'rho_sea_water' )
          IF ( av == 0 )  THEN
             to_be_resorted => rho_ocean
          ELSE
             IF ( .NOT. ALLOCATED( rho_ocean_av ) ) THEN
                ALLOCATE( rho_ocean_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                rho_ocean_av = REAL( fill_value, KIND = wp )
             ENDIF
             to_be_resorted => rho_ocean_av
          ENDIF

       CASE ( 'sa' )
          IF ( av == 0 )  THEN
             to_be_resorted => sa
          ELSE
             IF ( .NOT. ALLOCATED( sa_av ) ) THEN
                ALLOCATE( sa_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                sa_av = REAL( fill_value, KIND = wp )
             ENDIF
             to_be_resorted => sa_av
          ENDIF

       CASE DEFAULT
          found = .FALSE.

    END SELECT


    IF ( found  .AND.  .NOT. resorted )  THEN
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb_do, nzt_do
                local_pf(i,j,k) = MERGE( to_be_resorted(k,j,i),                &
                                   REAL( fill_value, KIND = wp ),              &
                                   BTEST( wall_flags_total_0(k,j,i), flag_nr ) )
             ENDDO
          ENDDO
       ENDDO
       resorted = .TRUE.
    ENDIF

 END SUBROUTINE ocean_data_output_3d


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Header output for ocean parameters
!------------------------------------------------------------------------------!
 SUBROUTINE ocean_header( io )


    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  io   !< Unit of the output file

!
!-- Write ocean header
    WRITE( io, 1 )
    IF ( stokes_force  )  WRITE( io, 2 )  stokes_waveheight, stokes_wavelength
    IF ( wave_breaking )  THEN
       WRITE( io, 3 )  alpha_wave_breaking, timescale_wave_breaking
    ENDIF
    IF ( .NOT. salinity )  WRITE( io, 4 )
    IF ( surface_cooling_spinup_time /= 999999.9_wp )  THEN
       WRITE( io, 5 )  surface_cooling_spinup_time
    ENDIF

1   FORMAT (//' Ocean settings:'/                                              &
              ' ------------------------------------------'/)
2   FORMAT ('    --> Craik-Leibovich vortex force and Stokes drift switched',  &
                     ' on'/                                                    &
            '        waveheight: ',F4.1,' m   wavelength: ',F6.1,' m')
3   FORMAT ('    --> wave breaking generated turbulence switched on'/          &
            '        alpha:    ',F4.1/                                         &
            '        timescale:',F5.1,' s')
4   FORMAT ('    --> prognostic salinity equation is switched off' )
5   FORMAT ('    --> surface heat flux is switched off after ',F8.1,' s')

 END SUBROUTINE ocean_header


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate arrays and assign pointers.
!------------------------------------------------------------------------------!
 SUBROUTINE ocean_init_arrays

    USE indices,                                                               &
        ONLY:  nxlg, nxrg, nyn, nyng, nys, nysg, nzb, nzt

    IMPLICIT NONE

    ALLOCATE( prho_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                           &
              rho_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                            &
              sa_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                             &
              sa_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

    IF (  salinity )  THEN
       ALLOCATE( sa_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ENDIF

    IF ( ws_scheme_sca )  THEN
       ALLOCATE( sums_wssas_ws_l(nzb:nzt+1,0:threads_per_task-1) )
       sums_wssas_ws_l = 0.0_wp
    ENDIF

    IF ( loop_optimization /= 'vector' )  THEN

       IF ( ws_scheme_sca )  THEN
          ALLOCATE( flux_s_sa(nzb+1:nzt,0:threads_per_task-1) )
          ALLOCATE( diss_s_sa(nzb+1:nzt,0:threads_per_task-1) )
          ALLOCATE( flux_l_sa(nzb+1:nzt,nys:nyn,0:threads_per_task-1) )
          ALLOCATE( diss_l_sa(nzb+1:nzt,nys:nyn,0:threads_per_task-1) )
       ENDIF

    ENDIF

    prho => prho_1
    rho_ocean  => rho_1  ! routines calc_mean_profile and diffusion_e require
                         ! density to be a pointer

!
!-- Initial assignment of pointers
    IF ( salinity )  THEN
       sa => sa_1;  sa_p => sa_2;  tsa_m => sa_3
    ELSE
       sa => sa_1;  sa_p => sa_1;  tsa_m => sa_3
    ENDIF

 END SUBROUTINE ocean_init_arrays


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of quantities needed for the ocean mode
!------------------------------------------------------------------------------!
 SUBROUTINE ocean_init


    USE arrays_3d,                                                             &
        ONLY:  dzu, dzw, hyp, pt_init, ref_state, u_stokes_zu, u_stokes_zw,    &
               v_stokes_zu, v_stokes_zw, zu, zw

    USE basic_constants_and_equations_mod,                                     &
        ONLY:  g

    USE basic_constants_and_equations_mod,                                     &
        ONLY:  pi

    USE control_parameters,                                                    &
        ONLY:  initializing_actions, molecular_viscosity, rho_surface,         &
               rho_reference, surface_pressure, top_momentumflux_u,            &
               top_momentumflux_v, use_single_reference_value

    USE indices,                                                               &
        ONLY:  nxl, nxlg, nxrg, nyng, nys, nysg, nzb, nzt

    USE kinds

    USE pegrid,                                                                &
        ONLY:  myid

    USE statistics,                                                            &
        ONLY:  hom, statistic_regions

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !< loop index
    INTEGER(iwp) ::  j  !< loop index
    INTEGER(iwp) ::  k  !< loop index
    INTEGER(iwp) ::  n  !< loop index

    REAL(wp) ::  alpha !< angle of surface stress
    REAL(wp) ::  dum   !< dummy argument
    REAL(wp) ::  pt_l  !< local scalar for pt used in equation of state function
    REAL(wp) ::  sa_l  !< local scalar for sa used in equation of state function
    REAL(wp) ::  velocity_amplitude  !< local scalar for amplitude of Stokes drift velocity
    REAL(wp) ::  x     !< temporary variable to store surface stress along x
    REAL(wp) ::  y     !< temporary variable to store surface stress along y

    REAL(wp), DIMENSION(nzb:nzt+1) ::  rho_ocean_init  !< local array for initial density

    ALLOCATE( hyp(nzb:nzt+1) )


!
!-- In case of no restart run, calculate the inital salinity profilevcusing the
!-- given salinity gradients
    IF ( TRIM( initializing_actions ) /= 'read_restart_data' )  THEN

       sa_init = sa_surface
!
!--    Last arguments gives back the gradient at top level to be used as
!--    possible Neumann boundary condition. This is not realized for the ocean
!--    mode, therefore a dummy argument is used.
       IF ( salinity )  THEN
          CALL init_vertical_profiles( sa_vertical_gradient_level_ind,          &
                                       sa_vertical_gradient_level,              &
                                       sa_vertical_gradient, sa_init,           &
                                       sa_surface, dum )
       ENDIF
    ENDIF

!
!-- Initialize required 3d-arrays
    IF ( TRIM( initializing_actions ) /= 'read_restart_data'  .AND.            &
         TRIM( initializing_actions ) /= 'cyclic_fill' )  THEN
!
!--    Initialization via computed 1D-model profiles
       IF ( INDEX( initializing_actions, 'set_constant_profiles' ) /= 0 )  THEN

          DO  i = nxlg, nxrg
             DO  j = nysg, nyng
                sa(:,j,i) = sa_init
             ENDDO
          ENDDO

       ENDIF
!
!--    Store initial profiles for output purposes etc.
!--    Store initial salinity profile
       hom(:,1,26,:)  = SPREAD( sa(:,nys,nxl), 2, statistic_regions+1 )
!
!--    Initialize old and new time levels.
       tsa_m = 0.0_wp
       sa_p  = sa

    ELSEIF ( TRIM( initializing_actions ) == 'read_restart_data' )  THEN

!
!--    Initialize new time levels (only done in order to set boundary values
!--    including ghost points)
       sa_p = sa
!
!--    Allthough tendency arrays are set in prognostic_equations, they have
!--    have to be predefined here because they are used (but multiplied with 0)
!--    there before they are set.
       tsa_m = 0.0_wp

    ENDIF

!
!-- Set water density near the ocean surface
    rho_surface = 1027.62_wp

!
!-- Set kinematic viscosity to sea water at 20C.
!-- This changes the default value that is given for air!
    molecular_viscosity = 1.05E-6_wp

!
!-- Change sign of buoyancy/stability terms because density gradient is used
!-- instead of the potential temperature gradient to calculate the buoyancy
    atmos_ocean_sign = -1.0_wp

!
!-- Calculate initial vertical profile of hydrostatic pressure (in Pa)
!-- and the reference density (used later in buoyancy term)
!-- First step: Calculate pressure using reference density
    hyp(nzt+1) = surface_pressure * 100.0_wp
    hyp(nzt)   = hyp(nzt+1) + rho_surface * g * 0.5_wp * dzu(nzt+1)
    rho_ocean_init(nzt)   = rho_surface
    rho_ocean_init(nzt+1) = rho_surface  ! only required for output

    DO  k = nzt-1, 1, -1
       hyp(k) = hyp(k+1) + rho_surface * g * dzu(k)
    ENDDO
    hyp(0) = hyp(1) + rho_surface * g * dzu(1)

!
!-- Second step: Iteratively calculate in situ density (based on presssure)
!-- and pressure (based on in situ density)
    DO  n = 1, 5

       rho_reference = rho_surface * 0.5_wp * dzu(nzt+1)

       DO  k = nzt, 0, -1

          sa_l = 0.5_wp * ( sa_init(k) + sa_init(k+1) )
          pt_l = 0.5_wp * ( pt_init(k) + pt_init(k+1) )

          rho_ocean_init(k) = eqn_state_seawater_func( hyp(k), pt_l, sa_l )

          rho_reference = rho_reference + rho_ocean_init(k) * dzu(k+1)

       ENDDO

       rho_reference = rho_reference / ( zw(nzt) - zu(nzb) )

       hyp(nzt) = hyp(nzt+1) + rho_surface * g * 0.5_wp * dzu(nzt+1)
       DO  k = nzt-1, 0, -1
          hyp(k) = hyp(k+1) + g * 0.5_wp * ( rho_ocean_init(k)                 &
                                           + rho_ocean_init(k+1) ) * dzu(k+1)
       ENDDO

    ENDDO

!
!-- Calculate the reference potential density
    prho_reference = 0.0_wp
    DO  k = 0, nzt

       sa_l = 0.5_wp * ( sa_init(k) + sa_init(k+1) )
       pt_l = 0.5_wp * ( pt_init(k) + pt_init(k+1) )

       prho_reference = prho_reference + dzu(k+1) * &
                        eqn_state_seawater_func( 0.0_wp, pt_l, sa_l )

    ENDDO

    prho_reference = prho_reference / ( zu(nzt) - zu(nzb) )

!
!-- Calculate the 3d array of initial in situ and potential density,
!-- based on the initial temperature and salinity profile
    CALL eqn_state_seawater

!
!-- Store initial density profile
    hom(:,1,77,:)  = SPREAD( rho_ocean_init(:), 2, statistic_regions+1 )

!
!-- Set the reference state to be used in the buoyancy terms
    IF ( use_single_reference_value )  THEN
       ref_state(:) = rho_reference
    ELSE
       ref_state(:) = rho_ocean_init(:)
    ENDIF

!
!-- Calculate the Stokes drift velocity profile
    IF ( stokes_force )  THEN

!
!--    First, calculate angle of surface stress
       x = -top_momentumflux_u
       y = -top_momentumflux_v
       IF ( x == 0.0_wp )  THEN
          IF ( y > 0.0_wp )  THEN
             alpha = pi / 2.0_wp
          ELSEIF ( y < 0.0_wp )  THEN
             alpha = 3.0_wp * pi / 2.0_wp
          ENDIF
       ELSE
          IF ( x < 0.0_wp )  THEN
             alpha = ATAN( y / x ) + pi
          ELSE
             IF ( y < 0.0_wp )  THEN
                alpha = ATAN( y / x ) + 2.0_wp * pi
             ELSE
                alpha = ATAN( y / x )
             ENDIF
          ENDIF
       ENDIF

       velocity_amplitude = ( pi * stokes_waveheight / stokes_wavelength )**2 *&
                            SQRT( g * stokes_wavelength / ( 2.0_wp * pi ) )

       DO  k = nzb, nzt
          u_stokes_zu(k) = velocity_amplitude * COS( alpha ) *                 &
                           EXP( 4.0_wp * pi * zu(k) / stokes_wavelength )
          u_stokes_zw(k) = velocity_amplitude * COS( alpha ) *                 &
                           EXP( 4.0_wp * pi * zw(k) / stokes_wavelength )
          v_stokes_zu(k) = velocity_amplitude * SIN( alpha ) *                 &
                           EXP( 4.0_wp * pi * zu(k) / stokes_wavelength )
          v_stokes_zw(k) = velocity_amplitude * SIN( alpha ) *                 &
                           EXP( 4.0_wp * pi * zw(k) / stokes_wavelength )
       ENDDO
       u_stokes_zu(nzt+1) = u_stokes_zw(nzt) ! because zu(nzt+1) changes the sign
       u_stokes_zw(nzt+1) = u_stokes_zw(nzt) ! because zw(nzt+1) changes the sign
       v_stokes_zu(nzt+1) = v_stokes_zw(nzt) ! because zu(nzt+1) changes the sign
       v_stokes_zw(nzt+1) = v_stokes_zw(nzt) ! because zw(nzt+1) changes the sign

    ENDIF

!
!-- Wave breaking effects
    IF ( wave_breaking )  THEN
!
!--    Calculate friction velocity at ocean surface
       u_star_wave_breaking = SQRT( SQRT( top_momentumflux_u**2 +              &
                                          top_momentumflux_v**2 ) )
!
!--    Set the time scale of random forcing. The vertical grid spacing at the
!--    ocean surface is assumed as the length scale of turbulence.
!--    Formula follows Noh et al. (2004), JPO
       timescale_wave_breaking = 0.1_wp * dzw(nzt) / alpha_wave_breaking /     &
                                 u_star_wave_breaking
!
!--    Set random number seeds differently on the processor cores in order to
!--    create different random number sequences
       iran_ocean = iran_ocean + myid
    ENDIF

 END SUBROUTINE ocean_init


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!------------------------------------------------------------------------------!
 SUBROUTINE ocean_actions( location )


    CHARACTER (LEN=*), INTENT(IN) ::  location !< call location string

    SELECT CASE ( location )

       CASE ( 'before_timestep' )

          IF ( ws_scheme_sca )  THEN
             sums_wssas_ws_l = 0.0_wp
          ENDIF

       CASE DEFAULT
          CONTINUE

    END SELECT

 END SUBROUTINE ocean_actions


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid points i,j
!------------------------------------------------------------------------------!
 SUBROUTINE ocean_actions_ij( i, j, location )


    INTEGER(iwp),      INTENT(IN) ::  i         !< grid index in x-direction
    INTEGER(iwp),      INTENT(IN) ::  j         !< grid index in y-direction
    CHARACTER (LEN=*), INTENT(IN) ::  location  !< call location string
    INTEGER(iwp)  ::  dummy  !< call location string

    IF ( ocean_mode )   dummy = i + j

    SELECT CASE ( location )

       CASE ( 'before_timestep' )

          IF ( ws_scheme_sca )  THEN
             sums_wssas_ws_l = 0.0_wp
          ENDIF

       CASE DEFAULT
          CONTINUE

    END SELECT

 END SUBROUTINE ocean_actions_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Prognostic equation for salinity.
!> Vector-optimized version
!------------------------------------------------------------------------------!
 SUBROUTINE ocean_prognostic_equations

    USE advec_s_bc_mod,                                                        &
        ONLY:  advec_s_bc

    USE advec_s_pw_mod,                                                        &
        ONLY:  advec_s_pw

    USE advec_s_up_mod,                                                        &
        ONLY:  advec_s_up

    USE advec_ws,                                                              &
        ONLY:  advec_s_ws

    USE arrays_3d,                                                             &
        ONLY:  rdf_sc, tend, tsa_m

    USE control_parameters,                                                    &
        ONLY:  bc_dirichlet_l,                                                 &
               bc_dirichlet_n,                                                 &
               bc_dirichlet_r,                                                 &
               bc_dirichlet_s,                                                 &
               bc_radiation_l,                                                 &
               bc_radiation_n,                                                 &
               bc_radiation_r,                                                 &
               bc_radiation_s,                                                 &
               dt_3d, intermediate_timestep_count,                             &
               intermediate_timestep_count_max, scalar_advec, simulated_time,  &
               timestep_scheme, tsc, ws_scheme_sca

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point_s

    USE diffusion_s_mod,                                                       &
        ONLY:  diffusion_s

    IMPLICIT NONE

    INTEGER(iwp) ::  i       !< loop index
    INTEGER(iwp) ::  j       !< loop index
    INTEGER(iwp) ::  k       !< loop index

    REAL(wp)     ::  sbt     !< weighting factor for sub-time step

!
!-- Switch of the surface heat flux, if requested
    IF ( surface_cooling_spinup_time /= 999999.9_wp )  THEN
       IF ( .NOT. surface_cooling_switched_off  .AND.                          &
            simulated_time >= surface_cooling_spinup_time )  THEN

          surf_def_h(2)%shf = 0.0_wp
          surface_cooling_switched_off = .TRUE.

       ENDIF
    ENDIF

!
!-- Compute prognostic equations for the ocean mode
!-- First, start with salinity
    IF ( salinity )  THEN

       CALL cpu_log( log_point_s(20), 'sa-equation', 'start' )

!
!--    sa-tendency terms with communication
       sbt = tsc(2)
       IF ( scalar_advec == 'bc-scheme' )  THEN

          IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--          Bott-Chlond scheme always uses Euler time step. Thus:
             sbt = 1.0_wp
          ENDIF
          tend = 0.0_wp
          CALL advec_s_bc( sa, 'sa' )

       ENDIF

!
!--    sa-tendency terms with no communication
       IF ( scalar_advec /= 'bc-scheme' )  THEN
          tend = 0.0_wp
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( ws_scheme_sca )  THEN
                CALL advec_s_ws( advc_flags_s, sa, 'sa',                       &
                                 bc_dirichlet_l  .OR.  bc_radiation_l,         &
                                 bc_dirichlet_n  .OR.  bc_radiation_n,         &
                                 bc_dirichlet_r  .OR.  bc_radiation_r,         &
                                 bc_dirichlet_s  .OR.  bc_radiation_s )
             ELSE
                CALL advec_s_pw( sa )
             ENDIF
          ELSE
             CALL advec_s_up( sa )
          ENDIF
       ENDIF

       CALL diffusion_s( sa,                                                   &
                         surf_def_h(0)%sasws, surf_def_h(1)%sasws,             &
                         surf_def_h(2)%sasws,                                  &
                         surf_lsm_h%sasws,    surf_usm_h%sasws,                &
                         surf_def_v(0)%sasws, surf_def_v(1)%sasws,             &
                         surf_def_v(2)%sasws, surf_def_v(3)%sasws,             &
                         surf_lsm_v(0)%sasws, surf_lsm_v(1)%sasws,             &
                         surf_lsm_v(2)%sasws, surf_lsm_v(3)%sasws,             &
                         surf_usm_v(0)%sasws, surf_usm_v(1)%sasws,             &
                         surf_usm_v(2)%sasws, surf_usm_v(3)%sasws )

!       CALL user_actions( 'sa-tendency' ) ToDo: find general solution for dependency between modules

!
!--    Prognostic equation for salinity
       DO  i = nxl, nxr
          DO  j = nys, nyn
             !following directive is required to vectorize on Intel19
             !DIR$ IVDEP
             DO  k = nzb+1, nzt
                sa_p(k,j,i) = sa(k,j,i) + ( dt_3d * ( sbt * tend(k,j,i) +            &
                                                      tsc(3) * tsa_m(k,j,i) )        &
                                                  - tsc(5) * rdf_sc(k) *             &
                                                    ( sa(k,j,i) - sa_init(k) )       &
                                          )                                          &
                                            * MERGE( 1.0_wp, 0.0_wp,                 &
                                               BTEST( wall_flags_total_0(k,j,i), 0 ) &
                                                   )
                IF ( sa_p(k,j,i) < 0.0_wp )  sa_p(k,j,i) = 0.1_wp * sa(k,j,i)
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
                      tsa_m(k,j,i) = tend(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSEIF ( intermediate_timestep_count < intermediate_timestep_count_max ) &
          THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt
                      tsa_m(k,j,i) =   -9.5625_wp * tend(k,j,i) +              &
                                        5.3125_wp * tsa_m(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDIF

       CALL cpu_log( log_point_s(20), 'sa-equation', 'stop' )

    ENDIF

!
!-- Calculate density by the equation of state for seawater
    CALL cpu_log( log_point_s(21), 'eqns-seawater', 'start' )
    CALL eqn_state_seawater
    CALL cpu_log( log_point_s(21), 'eqns-seawater', 'stop' )

 END SUBROUTINE ocean_prognostic_equations


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Prognostic equations for ocean mode (so far, salinity only)
!> Cache-optimized version
!------------------------------------------------------------------------------!
 SUBROUTINE ocean_prognostic_equations_ij( i, j, i_omp_start, tn )

    USE advec_s_pw_mod,                                                        &
        ONLY:  advec_s_pw

    USE advec_s_up_mod,                                                        &
        ONLY:  advec_s_up

    USE advec_ws,                                                              &
        ONLY:  advec_s_ws

    USE arrays_3d,                                                             &
        ONLY:  diss_l_sa, diss_s_sa, flux_l_sa, flux_s_sa, rdf_sc, tend, tsa_m

    USE control_parameters,                                                    &
        ONLY:  bc_dirichlet_l,                                                 &
               bc_dirichlet_n,                                                 &
               bc_dirichlet_r,                                                 &
               bc_dirichlet_s,                                                 &
               bc_radiation_l,                                                 &
               bc_radiation_n,                                                 &
               bc_radiation_r,                                                 &
               bc_radiation_s,                                                 &
               dt_3d, intermediate_timestep_count,                             &
               intermediate_timestep_count_max, simulated_time,                &
               timestep_scheme, tsc, ws_scheme_sca

    USE diffusion_s_mod,                                                       &
        ONLY:  diffusion_s

    IMPLICIT NONE

    INTEGER(iwp) ::  i             !< loop index x direction
    INTEGER(iwp) ::  i_omp_start   !< first loop index of i-loop in calling    &
                                   !< routine prognostic_equations
    INTEGER(iwp) ::  j             !< loop index y direction
    INTEGER(iwp) ::  k             !< loop index z direction
    INTEGER(iwp) ::  tn            !< task number of openmp task


!
!-- Switch of the surface heat flux, if requested
    IF ( surface_cooling_spinup_time /= 999999.9_wp )  THEN
       IF ( .NOT. surface_cooling_switched_off  .AND.                          &
            simulated_time >= surface_cooling_spinup_time )  THEN

          surf_def_h(2)%shf = 0.0_wp
          surface_cooling_switched_off = .TRUE.

       ENDIF
    ENDIF

!
!-- Compute prognostic equations for the ocean mode
!-- First, start with tendency-terms for salinity
    IF ( salinity )  THEN

       tend(:,j,i) = 0.0_wp
       IF ( timestep_scheme(1:5) == 'runge' ) &
       THEN
          IF ( ws_scheme_sca )  THEN
             CALL advec_s_ws( advc_flags_s,                                    &
                              i, j, sa, 'sa', flux_s_sa,  diss_s_sa, flux_l_sa,&
                              diss_l_sa, i_omp_start, tn,                      &
                              bc_dirichlet_l  .OR.  bc_radiation_l,            &
                              bc_dirichlet_n  .OR.  bc_radiation_n,            &
                              bc_dirichlet_r  .OR.  bc_radiation_r,            &
                              bc_dirichlet_s  .OR.  bc_radiation_s )
          ELSE
             CALL advec_s_pw( i, j, sa )
          ENDIF
       ELSE
          CALL advec_s_up( i, j, sa )
       ENDIF
       CALL diffusion_s( i, j, sa,                                             &
                         surf_def_h(0)%sasws, surf_def_h(1)%sasws,             &
                         surf_def_h(2)%sasws,                                  &
                         surf_lsm_h%sasws,    surf_usm_h%sasws,                &
                         surf_def_v(0)%sasws, surf_def_v(1)%sasws,             &
                         surf_def_v(2)%sasws, surf_def_v(3)%sasws,             &
                         surf_lsm_v(0)%sasws, surf_lsm_v(1)%sasws,             &
                         surf_lsm_v(2)%sasws, surf_lsm_v(3)%sasws,             &
                         surf_usm_v(0)%sasws, surf_usm_v(1)%sasws,             &
                         surf_usm_v(2)%sasws, surf_usm_v(3)%sasws )

!       CALL user_actions( i, j, 'sa-tendency' ) ToDo: find general solution for dependency between modules

!
!--    Prognostic equation for salinity
       DO  k = nzb+1, nzt

          sa_p(k,j,i) = sa(k,j,i) + ( dt_3d *                                  &
                                              ( tsc(2) * tend(k,j,i) +         &
                                                tsc(3) * tsa_m(k,j,i) )        &
                                    - tsc(5) * rdf_sc(k)                       &
                                             * ( sa(k,j,i) - sa_init(k) )      &
                                    ) * MERGE( 1.0_wp, 0.0_wp,                 &
                                         BTEST( wall_flags_total_0(k,j,i), 0 ) )

          IF ( sa_p(k,j,i) < 0.0_wp )  sa_p(k,j,i) = 0.1_wp * sa(k,j,i)

       ENDDO

!
!--    Calculate tendencies for the next Runge-Kutta step
       IF ( timestep_scheme(1:5) == 'runge' )  THEN
          IF ( intermediate_timestep_count == 1 )  THEN
             DO  k = nzb+1, nzt
                tsa_m(k,j,i) = tend(k,j,i)
             ENDDO
          ELSEIF ( intermediate_timestep_count < intermediate_timestep_count_max ) &
          THEN
             DO  k = nzb+1, nzt
                tsa_m(k,j,i) =   -9.5625_wp * tend(k,j,i) +                    &
                                  5.3125_wp * tsa_m(k,j,i)
             ENDDO
          ENDIF
       ENDIF

    ENDIF

!
!-- Calculate density by the equation of state for seawater
    CALL eqn_state_seawater( i, j )

 END SUBROUTINE ocean_prognostic_equations_ij

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Boundary conditions for ocean model
!------------------------------------------------------------------------------!
 SUBROUTINE ocean_boundary_conditions

    IMPLICIT NONE

    INTEGER(iwp) ::  i                            !< grid index x direction.
    INTEGER(iwp) ::  j                            !< grid index y direction.
    INTEGER(iwp) ::  k                            !< grid index z direction.
    INTEGER(iwp) ::  l                            !< running index boundary type, for up- and downward-facing walls.
    INTEGER(iwp) ::  m                            !< running index surface elements.

!
!--    Boundary conditions for salinity
!--    Bottom boundary: Neumann condition because salinity flux is always
!--    given.
       DO  l = 0, 1
          !$OMP PARALLEL DO PRIVATE( i, j, k )
          DO  m = 1, bc_h(l)%ns
             i = bc_h(l)%i(m)
             j = bc_h(l)%j(m)
             k = bc_h(l)%k(m)
             sa_p(k+bc_h(l)%koff,j,i) = sa_p(k,j,i)
          ENDDO
       ENDDO
!
!--    Top boundary: Dirichlet or Neumann
       IF ( ibc_sa_t == 0 )  THEN
           sa_p(nzt+1,:,:) = sa(nzt+1,:,:)
       ELSEIF ( ibc_sa_t == 1 )  THEN
           sa_p(nzt+1,:,:) = sa_p(nzt,:,:)
       ENDIF

 END SUBROUTINE ocean_boundary_conditions

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Swapping of timelevels.
!------------------------------------------------------------------------------!
 SUBROUTINE ocean_swap_timelevel( mod_count )

    IMPLICIT NONE

    INTEGER, INTENT(IN) ::  mod_count  !< flag defining where pointers point to


    SELECT CASE ( mod_count )

       CASE ( 0 )
          IF ( salinity )  THEN
             sa => sa_1;    sa_p => sa_2
          ENDIF

       CASE ( 1 )
          IF ( salinity )  THEN
             sa => sa_2;    sa_p => sa_1
          ENDIF

    END SELECT

 END SUBROUTINE ocean_swap_timelevel


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine reads the respective restart data for the ocean module.
!------------------------------------------------------------------------------!
 SUBROUTINE ocean_rrd_global( found )


    USE control_parameters,                                                    &
        ONLY: length, restart_string


    IMPLICIT NONE

    LOGICAL, INTENT(OUT)  ::  found


    found = .TRUE.

    SELECT CASE ( restart_string(1:length) )

       CASE ( 'bc_sa_t' )
          READ ( 13 )  bc_sa_t

       CASE ( 'bottom_salinityflux' )
          READ ( 13 )  bottom_salinityflux

       CASE ( 'salinity' )
          READ ( 13 )  salinity

       CASE ( 'sa_init' )
          READ ( 13 )  sa_init

       CASE ( 'sa_surface' )
          READ ( 13 )  sa_surface

       CASE ( 'sa_vertical_gradient' )
          READ ( 13 )  sa_vertical_gradient

       CASE ( 'sa_vertical_gradient_level' )
          READ ( 13 )  sa_vertical_gradient_level

       CASE ( 'stokes_waveheight' )
          READ ( 13 )  stokes_waveheight

       CASE ( 'stokes_wavelength' )
          READ ( 13 )  stokes_wavelength

       CASE ( 'surface_cooling_spinup_time' )
          READ ( 13 )  surface_cooling_spinup_time

       CASE ( 'top_salinityflux' )
          READ ( 13 )  top_salinityflux

       CASE ( 'wall_salinityflux' )
          READ ( 13 )  wall_salinityflux

       CASE ( 'wave_breaking' )
          READ ( 13 )  wave_breaking

       CASE DEFAULT

          found = .FALSE.

    END SELECT

 END SUBROUTINE ocean_rrd_global


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine reads the respective restart data for the ocean module.
!------------------------------------------------------------------------------!
 SUBROUTINE ocean_rrd_local( k, nxlf, nxlc, nxl_on_file, nxrf, nxrc,           &
                             nxr_on_file, nynf, nync, nyn_on_file, nysf,       &
                             nysc, nys_on_file, tmp_3d, found )

    USE averaging,                                                             &
        ONLY:  rho_ocean_av, sa_av

    USE control_parameters,                                                    &
        ONLY:  length, restart_string

    USE indices,                                                               &
        ONLY:  nbgp, nxlg, nxrg, nyng, nysg, nzb, nzt

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

    REAL(wp), DIMENSION(nzb:nzt+1,nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) :: tmp_3d   !<


    found = .TRUE.

    SELECT CASE ( restart_string(1:length) )

       CASE ( 'rho_ocean_av' )
          IF ( .NOT. ALLOCATED( rho_ocean_av ) )  THEN
             ALLOCATE( rho_ocean_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
          rho_ocean_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =            &
                              tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'sa' )
          IF ( k == 1 )  READ ( 13 )  tmp_3d
          sa(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                      &
                              tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'sa_av' )
          IF ( .NOT. ALLOCATED( sa_av ) )  THEN
              ALLOCATE( sa_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
          sa_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                   &
                              tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE DEFAULT
          found = .FALSE.

    END SELECT

 END SUBROUTINE ocean_rrd_local


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine writes the respective restart data for the ocean module.
!------------------------------------------------------------------------------!
 SUBROUTINE ocean_wrd_global


    IMPLICIT NONE

    CALL wrd_write_string( 'bc_sa_t' )
    WRITE ( 14 )  bc_sa_t

    CALL wrd_write_string( 'bottom_salinityflux' )
    WRITE ( 14 )  bottom_salinityflux

    CALL wrd_write_string( 'salinity' )
    WRITE ( 14 )  salinity

    CALL wrd_write_string( 'sa_init' )
    WRITE ( 14 )  sa_init

    CALL wrd_write_string( 'sa_surface' )
    WRITE ( 14 )  sa_surface

    CALL wrd_write_string( 'sa_vertical_gradient' )
    WRITE ( 14 )  sa_vertical_gradient

    CALL wrd_write_string( 'sa_vertical_gradient_level' )
    WRITE ( 14 )  sa_vertical_gradient_level

    CALL wrd_write_string( 'stokes_waveheight' )
    WRITE ( 14 )  stokes_waveheight

    CALL wrd_write_string( 'stokes_wavelength' )
    WRITE ( 14 )  stokes_wavelength

    CALL wrd_write_string( 'surface_cooling_spinup_time' )
    WRITE ( 14 )  surface_cooling_spinup_time

    CALL wrd_write_string( 'top_salinityflux' )
    WRITE ( 14 )  top_salinityflux

    CALL wrd_write_string( 'wall_salinityflux' )
    WRITE ( 14 )  wall_salinityflux

    CALL wrd_write_string( 'wave_breaking' )
    WRITE ( 14 )  wave_breaking

 END SUBROUTINE ocean_wrd_global


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine writes the respective restart data for the ocean module.
!------------------------------------------------------------------------------!
 SUBROUTINE ocean_wrd_local

    USE averaging,                                                             &
        ONLY:  rho_ocean_av, sa_av

    IMPLICIT NONE

    IF ( ALLOCATED( rho_ocean_av ) )  THEN
       CALL wrd_write_string( 'rho_ocean_av' )
       WRITE ( 14 )  rho_ocean_av
    ENDIF

    CALL wrd_write_string( 'sa' )
    WRITE ( 14 )  sa

    IF ( ALLOCATED( sa_av ) )  THEN
       CALL wrd_write_string( 'sa_av' )
       WRITE ( 14 )  sa_av
    ENDIF

 END SUBROUTINE ocean_wrd_local


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine calculates the Craik Leibovich vortex force and the additional
!> effect of the Stokes drift on the Coriolis force
!> Call for all gridpoints.
!------------------------------------------------------------------------------!
 SUBROUTINE stokes_drift_terms( component )

    USE arrays_3d,                                                             &
        ONLY:  ddzu, u, u_stokes_zu, u_stokes_zw, v, v_stokes_zu,              &
               v_stokes_zw, w, tend

    USE basic_constants_and_equations_mod,                                     &
        ONLY:  pi

    USE control_parameters,                                                    &
        ONLY:  f, fs, message_string, rotation_angle

    USE grid_variables,                                                        &
        ONLY:  ddx, ddy

    USE indices,                                                               &
        ONLY:  nxl, nxr, nys, nysv, nyn, nzb, nzt

    IMPLICIT NONE

    INTEGER(iwp) ::  component      !< component of momentum equation
    INTEGER(iwp) ::  i              !< loop index along x
    INTEGER(iwp) ::  j              !< loop index along y
    INTEGER(iwp) ::  k              !< loop index along z

    REAL(wp)     ::  cos_rot_angle  !< cosine of model rotation angle
    REAL(wp)     ::  sin_rot_angle  !< sine of model rotation angle

!
!-- Compute Stokes terms for the respective velocity components
    SELECT CASE ( component )

!
!--    u-component
       CASE ( 1 )
          DO  i = nxl, nxr
             DO  j = nysv, nyn
                DO  k = nzb+1, nzt
                   tend(k,j,i) = tend(k,j,i) + v_stokes_zu(k) * (              &
                                   0.5 * ( v(k,j+1,i) - v(k,j+1,i-1)           &
                                         + v(k,j,i)   - v(k,j,i-1)   ) * ddx   &
                                 - 0.5 * ( u(k,j+1,i) - u(k,j-1,i) )   * ddy   &
                                                                )              &
                                 + f * v_stokes_zu(k)
                ENDDO
             ENDDO
          ENDDO

!
!--    v-component
       CASE ( 2 )
          DO  i = nxl, nxr
             DO  j = nysv, nyn
                DO  k = nzb+1, nzt
                   tend(k,j,i) = tend(k,j,i) - u_stokes_zu(k) * (              &
                                   0.5 * ( v(k,j,i+1) - v(k,j,i-1) )   * ddx   &
                                 - 0.5 * ( u(k,j,i) - u(k,j-1,i)               &
                                         + u(k,j,i+1) - u(k,j-1,i+1) ) * ddy   &
                                                                )              &
                                 - f * u_stokes_zu(k)
                ENDDO
             ENDDO
          ENDDO

!
!--    w-component
       CASE ( 3 )

!
!--       Precalculate cosine and sine of rotation angle
          cos_rot_angle = COS( rotation_angle * pi / 180.0_wp )
          sin_rot_angle = SIN( rotation_angle * pi / 180.0_wp )

          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   tend(k,j,i) = tend(k,j,i) + u_stokes_zw(k) * (              &
                                             0.5 * ( u(k+1,j,i) - u(k,j,i)     &
                                                   + u(k+1,j,i+1) - u(k,j,i+1) &
                                                   ) * ddzu(k+1)               &
                                           - 0.5 * ( w(k,j,i+1) - w(k,j,i-1)   &
                                                   ) * ddx      )              &
                                             - v_stokes_zw(k) * (              &
                                             0.5 * ( w(k,j+1,i) - w(k,j-1,i)   &
                                                   ) * ddy                     &
                                           - 0.5 * ( v(k+1,j,i) - v(k,j,i)     &
                                                   + v(k+1,j+1,i) - v(k,j+1,i) &
                                                   ) * ddzu(k)  )              &
                                           + fs * (                            &
                                               sin_rot_angle * v_stokes_zw(k)  &
                                             + cos_rot_angle * u_stokes_zw(k)  &
                                                  )
                ENDDO
             ENDDO
          ENDDO

       CASE DEFAULT
          WRITE( message_string, * ) 'wrong component of Stokes force: ',      &
                                     component
          CALL message( 'stokes_drift_terms', 'PA0091', 1, 2, 0, 6, 0 )

    END SELECT

 END SUBROUTINE stokes_drift_terms


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine calculates the Craik Leibovich vortex force and the additional
!> effect of the Stokes drift on the Coriolis force
!> Call for gridpoints i,j.
!------------------------------------------------------------------------------!

 SUBROUTINE stokes_drift_terms_ij( i, j, component )

    USE arrays_3d,                                                             &
        ONLY:  ddzu, u, u_stokes_zu, u_stokes_zw, v, v_stokes_zu,              &
               v_stokes_zw, w, tend

    USE basic_constants_and_equations_mod,                                     &
        ONLY:  pi

    USE control_parameters,                                                    &
        ONLY:  f, fs, message_string, rotation_angle

    USE grid_variables,                                                        &
        ONLY:  ddx, ddy

    USE indices,                                                               &
        ONLY:  nzb, nzt

    IMPLICIT NONE

    INTEGER(iwp) ::  component      !< component of momentum equation
    INTEGER(iwp) ::  i              !< loop index along x
    INTEGER(iwp) ::  j              !< loop index along y
    INTEGER(iwp) ::  k              !< loop incex along z

    REAL(wp)     ::  cos_rot_angle  !< cosine of model rotation angle
    REAL(wp)     ::  sin_rot_angle  !< sine of model rotation angle

!
!-- Compute Stokes terms for the respective velocity components
    SELECT CASE ( component )

!
!--    u-component
       CASE ( 1 )
          DO  k = nzb+1, nzt
             tend(k,j,i) = tend(k,j,i) + v_stokes_zu(k) * (                    &
                                     0.5 * ( v(k,j+1,i) - v(k,j+1,i-1)         &
                                           + v(k,j,i)   - v(k,j,i-1)   ) * ddx &
                                   - 0.5 * ( u(k,j+1,i) - u(k,j-1,i) )   * ddy &
                                                          )                    &
                                       + f * v_stokes_zu(k)
          ENDDO
!
!--    v-component
       CASE ( 2 )
          DO  k = nzb+1, nzt
             tend(k,j,i) = tend(k,j,i) - u_stokes_zu(k) * (                    &
                                     0.5 * ( v(k,j,i+1) - v(k,j,i-1) )   * ddx &
                                   - 0.5 * ( u(k,j,i) - u(k,j-1,i)             &
                                           + u(k,j,i+1) - u(k,j-1,i+1) ) * ddy &
                                                          )                    &
                                       - f * u_stokes_zu(k)
          ENDDO

!
!--    w-component
       CASE ( 3 )

!
!--       Precalculate cosine and sine of rotation angle
          cos_rot_angle = COS( rotation_angle * pi / 180.0_wp )
          sin_rot_angle = SIN( rotation_angle * pi / 180.0_wp )

          DO  k = nzb+1, nzt
             tend(k,j,i) = tend(k,j,i) + u_stokes_zw(k) * (              &
                                     0.5 * ( u(k+1,j,i) - u(k,j,i)     &
                                                   + u(k+1,j,i+1) - u(k,j,i+1) &
                                                   ) * ddzu(k+1)               &
                                           - 0.5 * ( w(k,j,i+1) - w(k,j,i-1)   &
                                                   ) * ddx )                   &
                                       - v_stokes_zw(k) * (                    &
                                             0.5 * ( w(k,j+1,i) - w(k,j-1,i)   &
                                                   ) * ddy                     &
                                           - 0.5 * ( v(k+1,j,i) - v(k,j,i)     &
                                                   + v(k+1,j+1,i) - v(k,j+1,i) &
                                                   ) * ddzu(k)  )              &
                                       + fs * ( sin_rot_angle * v_stokes_zw(k) &
                                              + cos_rot_angle * u_stokes_zw(k) &
                                              )
          ENDDO

       CASE DEFAULT
          WRITE( message_string, * ) ' wrong component: ', component
          CALL message( 'stokes_drift_terms', 'PA0091', 1, 2, 0, 6, 0 )

    END SELECT

 END SUBROUTINE stokes_drift_terms_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine calculates turbulence generated by wave breaking near the ocean
!> surface, following a parameterization given in Noh et al. (2004), JPO
!> Call for all gridpoints.
!> TODO: so far, this routine only works if the model time step has about the
!>       same value as the time scale of wave breaking!
!------------------------------------------------------------------------------!
 SUBROUTINE wave_breaking_term( component )

    USE arrays_3d,                                                             &
        ONLY:  u_p, v_p

    USE control_parameters,                                                    &
        ONLY:  dt_3d, message_string

    USE indices,                                                               &
        ONLY:  nxl, nxlu, nxr, nys, nysv, nyn, nzt

    IMPLICIT NONE

    INTEGER(iwp) ::  component  !< component of momentum equation
    INTEGER(iwp) ::  i          !< loop index along x
    INTEGER(iwp) ::  j          !< loop index along y

    REAL(wp) ::  random_gauss  !< function that creates a random number with a
                               !< Gaussian distribution


!
!-- Compute wave breaking terms for the respective velocity components.
!-- Velocities are directly manipulated, since this is not a real force
    SELECT CASE ( component )

!
!--    u-component
       CASE ( 1 )
          DO  i = nxlu, nxr
             DO  j = nys, nyn
                u_p(nzt,j,i) = u_p(nzt,j,i) +                                  &
                               ( random_gauss( iran_ocean, 1.0_wp ) - 1.0_wp ) &
                               * alpha_wave_breaking * u_star_wave_breaking    &
                               / timescale_wave_breaking * dt_3d
             ENDDO
          ENDDO
!
!--    v-component
       CASE ( 2 )
          DO  i = nxl, nxr
             DO  j = nysv, nyn
                v_p(nzt,j,i) = v_p(nzt,j,i) +                                  &
                               ( random_gauss( iran_ocean, 1.0_wp ) - 1.0_wp ) &
                               * alpha_wave_breaking * u_star_wave_breaking    &
                               / timescale_wave_breaking * dt_3d
             ENDDO
          ENDDO

       CASE DEFAULT
          WRITE( message_string, * ) 'wrong component of wave breaking: ',     &
                                     component
          CALL message( 'stokes_drift_terms', 'PA0466', 1, 2, 0, 6, 0 )

    END SELECT

 END SUBROUTINE wave_breaking_term


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine calculates turbulence generated by wave breaking near the ocean
!> surface, following a parameterization given in Noh et al. (2004), JPO
!> Call for gridpoint i,j.
!> TODO: so far, this routine only works if the model time step has about the
!>       same value as the time scale of wave breaking!
!------------------------------------------------------------------------------!
 SUBROUTINE wave_breaking_term_ij( i, j, component )

    USE arrays_3d,                                                             &
        ONLY:  u_p, v_p

    USE control_parameters,                                                    &
        ONLY:  dt_3d, message_string

    USE indices,                                                               &
        ONLY:  nzt

    IMPLICIT NONE

    INTEGER(iwp) ::  component  !< component of momentum equation
    INTEGER(iwp) ::  i          !< loop index along x
    INTEGER(iwp) ::  j          !< loop index along y

    REAL(wp) ::  random_gauss  !< function that creates a random number with a
                               !< Gaussian distribution

!
!-- Compute wave breaking terms for the respective velocity components
    SELECT CASE ( component )

!
!--    u-/v-component
       CASE ( 1 )
          u_p(nzt,j,i) = u_p(nzt,j,i) +                                        &
                         ( random_gauss( iran_ocean, 1.0_wp ) - 1.0_wp )       &
                         * alpha_wave_breaking * u_star_wave_breaking          &
                         / timescale_wave_breaking * dt_3d

       CASE ( 2 )
          v_p(nzt,j,i) = v_p(nzt,j,i) +                                        &
                         ( random_gauss( iran_ocean, 1.0_wp ) - 1.0_wp )       &
                         * alpha_wave_breaking * u_star_wave_breaking          &
                         / timescale_wave_breaking * dt_3d

       CASE DEFAULT
          WRITE( message_string, * ) 'wrong component of wave breaking: ',     &
                                     component
          CALL message( 'stokes_drift_terms', 'PA0466', 1, 2, 0, 6, 0 )

    END SELECT

 END SUBROUTINE wave_breaking_term_ij


 END MODULE ocean_mod
