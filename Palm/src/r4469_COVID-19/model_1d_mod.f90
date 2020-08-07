!> @file model_1d_mod.f90
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
! $Id: model_1d_mod.f90 4449 2020-03-09 14:43:16Z suehring $
! Set intermediate_timestep_count back to zero after 1D-model integration.
! This is required e.g. for initial calls of calc_mean_profile.
!
! 4360 2020-01-07 11:25:50Z suehring
! Corrected "Former revisions" section
!
! 3655 2019-01-07 16:51:22Z knoop
! Modularization of all bulk cloud physics code components
!
! Revision 1.1  1998/03/09 16:22:10  raasch
! Initial revision
!
!
! Description:
! ------------
!> 1D-model to initialize the 3D-arrays.
!> The temperature profile is set as steady and a corresponding steady solution
!> of the wind profile is being computed.
!> All subroutines required can be found within this file.
!>
!> @todo harmonize code with new surface_layer_fluxes module
!> @bug 1D model crashes when using small grid spacings in the order of 1 m
!> @fixme option "as_in_3d_model" seems to be an inappropriate option because
!>        the 1D model uses different turbulence closure approaches at least if
!>        the 3D model is set to LES-mode.
!------------------------------------------------------------------------------!
 MODULE model_1d_mod

    USE arrays_3d,                                                             &
        ONLY:  dd2zu, ddzu, ddzw, dzu, dzw, pt_init, q_init, ug, u_init,       &
               vg, v_init, zu

    USE basic_constants_and_equations_mod,                                     &
        ONLY:  g, kappa, pi

    USE control_parameters,                                                    &
        ONLY:  constant_diffusion, constant_flux_layer, dissipation_1d, f,     &
               humidity, ibc_e_b, intermediate_timestep_count,                 &
               intermediate_timestep_count_max, km_constant,                   &
               message_string, mixing_length_1d, prandtl_number,               &
               roughness_length, run_description_header, simulated_time_chr,   &
               timestep_scheme, tsc, z0h_factor

    USE indices,                                                               &
        ONLY:  nzb, nzb_diff, nzt

    USE kinds

    USE pegrid,                                                                &
        ONLY:  myid


    IMPLICIT NONE

    INTEGER(iwp) ::  current_timestep_number_1d = 0  !< current timestep number (1d-model)
    INTEGER(iwp) ::  damp_level_ind_1d               !< lower grid index of damping layer (1d-model)

    LOGICAL ::  run_control_header_1d = .FALSE.  !< flag for output of run control header (1d-model)
    LOGICAL ::  stop_dt_1d = .FALSE.             !< termination flag, used in case of too small timestep (1d-model)

    REAL(wp) ::  alpha_buoyancy                !< model constant according to Koblitz (2013)
    REAL(wp) ::  c_0 = 0.416179145_wp          !< = 0.03^0.25; model constant according to Koblitz (2013)
    REAL(wp) ::  c_1 = 1.52_wp                 !< model constant according to Koblitz (2013)
    REAL(wp) ::  c_2 = 1.83_wp                 !< model constant according to Koblitz (2013)
    REAL(wp) ::  c_3                           !< model constant
    REAL(wp) ::  c_mu                          !< model constant
    REAL(wp) ::  damp_level_1d = -1.0_wp       !< namelist parameter
    REAL(wp) ::  dt_1d = 60.0_wp               !< dynamic timestep (1d-model)
    REAL(wp) ::  dt_max_1d = 300.0_wp          !< timestep limit (1d-model)
    REAL(wp) ::  dt_pr_1d = 9999999.9_wp       !< namelist parameter
    REAL(wp) ::  dt_run_control_1d = 60.0_wp   !< namelist parameter
    REAL(wp) ::  end_time_1d = 864000.0_wp     !< namelist parameter
    REAL(wp) ::  lambda                        !< maximum mixing length
    REAL(wp) ::  qs1d                          !< characteristic humidity scale (1d-model)
    REAL(wp) ::  simulated_time_1d = 0.0_wp    !< updated simulated time (1d-model)
    REAL(wp) ::  sig_diss = 2.95_wp            !< model constant according to Koblitz (2013)
    REAL(wp) ::  sig_e = 2.95_wp               !< model constant according to Koblitz (2013)
    REAL(wp) ::  time_pr_1d = 0.0_wp           !< updated simulated time for profile output (1d-model)
    REAL(wp) ::  time_run_control_1d = 0.0_wp  !< updated simulated time for run-control output (1d-model)
    REAL(wp) ::  ts1d                          !< characteristic temperature scale (1d-model)
    REAL(wp) ::  us1d                          !< friction velocity (1d-model)
    REAL(wp) ::  usws1d                        !< u-component of the momentum flux (1d-model)
    REAL(wp) ::  vsws1d                        !< v-component of the momentum flux (1d-model)
    REAL(wp) ::  z01d                          !< roughness length for momentum (1d-model)
    REAL(wp) ::  z0h1d                         !< roughness length for scalars (1d-model)

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  diss1d   !< tke dissipation rate (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  diss1d_p !< prognostic value of tke dissipation rate (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  e1d      !< tke (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  e1d_p    !< prognostic value of tke (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  kh1d     !< turbulent diffusion coefficient for heat (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  km1d     !< turbulent diffusion coefficient for momentum (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  l1d      !< mixing length for turbulent diffusion coefficients (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  l1d_init !< initial mixing length (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  l1d_diss !< mixing length for dissipation (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  rif1d    !< Richardson flux number (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  te_diss  !< tendency of diss (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  te_dissm !< weighted tendency of diss for previous sub-timestep (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  te_e     !< tendency of e (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  te_em    !< weighted tendency of e for previous sub-timestep (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  te_u     !< tendency of u (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  te_um    !< weighted tendency of u for previous sub-timestep (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  te_v     !< tendency of v (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  te_vm    !< weighted tendency of v for previous sub-timestep (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  u1d      !< u-velocity component (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  u1d_p    !< prognostic value of u-velocity component (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  v1d      !< v-velocity component (1d-model)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  v1d_p    !< prognostic value of v-velocity component (1d-model)

!
!-- Initialize 1D model
    INTERFACE init_1d_model
       MODULE PROCEDURE init_1d_model
    END INTERFACE init_1d_model

!
!-- Print profiles
    INTERFACE print_1d_model
       MODULE PROCEDURE print_1d_model
    END INTERFACE print_1d_model

!
!-- Print run control information
    INTERFACE run_control_1d
       MODULE PROCEDURE run_control_1d
    END INTERFACE run_control_1d

!
!-- Main procedure
    INTERFACE time_integration_1d
       MODULE PROCEDURE time_integration_1d
    END INTERFACE time_integration_1d

!
!-- Calculate time step
    INTERFACE timestep_1d
       MODULE PROCEDURE timestep_1d
    END INTERFACE timestep_1d

    SAVE

    PRIVATE
!
!-- Public interfaces
    PUBLIC  init_1d_model

!
!-- Public variables
    PUBLIC  damp_level_1d, damp_level_ind_1d, diss1d, dt_pr_1d,                &
            dt_run_control_1d, e1d, end_time_1d, kh1d, km1d, l1d, rif1d, u1d,  &
            us1d, usws1d, v1d, vsws1d


    CONTAINS

 SUBROUTINE init_1d_model

    USE grid_variables,                                                        &
        ONLY:  dx, dy

    IMPLICIT NONE

    CHARACTER (LEN=9) ::  time_to_string  !< function to transform time from real to character string

    INTEGER(iwp) ::  k  !< loop index

!
!-- Allocate required 1D-arrays
    ALLOCATE( diss1d(nzb:nzt+1), diss1d_p(nzb:nzt+1),                          &
              e1d(nzb:nzt+1), e1d_p(nzb:nzt+1), kh1d(nzb:nzt+1),               &
              km1d(nzb:nzt+1), l1d(nzb:nzt+1), l1d_init(nzb:nzt+1),            &
              l1d_diss(nzb:nzt+1), rif1d(nzb:nzt+1), te_diss(nzb:nzt+1),       &
              te_dissm(nzb:nzt+1), te_e(nzb:nzt+1),                            &
              te_em(nzb:nzt+1), te_u(nzb:nzt+1), te_um(nzb:nzt+1),             &
              te_v(nzb:nzt+1), te_vm(nzb:nzt+1), u1d(nzb:nzt+1),               &
              u1d_p(nzb:nzt+1),  v1d(nzb:nzt+1), v1d_p(nzb:nzt+1) )

!
!-- Initialize arrays
    IF ( constant_diffusion )  THEN
       km1d = km_constant
       kh1d = km_constant / prandtl_number
    ELSE
       diss1d = 0.0_wp; diss1d_p = 0.0_wp
       e1d = 0.0_wp; e1d_p = 0.0_wp
       kh1d = 0.0_wp; km1d = 0.0_wp
       rif1d = 0.0_wp
!
!--    Compute the mixing length
       l1d_init(nzb) = 0.0_wp

       IF ( TRIM( mixing_length_1d ) == 'blackadar' )  THEN
!
!--       Blackadar mixing length
          IF ( f /= 0.0_wp )  THEN
             lambda = 2.7E-4_wp * SQRT( ug(nzt+1)**2 + vg(nzt+1)**2 ) /        &
                               ABS( f ) + 1E-10_wp
          ELSE
             lambda = 30.0_wp
          ENDIF

          DO  k = nzb+1, nzt+1
             l1d_init(k) = kappa * zu(k) / ( 1.0_wp + kappa * zu(k) / lambda )
          ENDDO

       ELSEIF ( TRIM( mixing_length_1d ) == 'as_in_3d_model' )  THEN
!
!--       Use the same mixing length as in 3D model (LES-mode)
          !> @todo rename (delete?) this option
          !>  As the mixing length is different between RANS and LES mode, it
          !>  must be distinguished here between these modes. For RANS mode,
          !>  the mixing length is calculated accoding to Blackadar, which is
          !>  the other option at this point.
          !>  Maybe delete this option entirely (not appropriate in LES case)
          !>  2018-03-20, gronemeier
          DO  k = nzb+1, nzt
             l1d_init(k)  = ( dx * dy * dzw(k) )**0.33333333333333_wp
          ENDDO
          l1d_init(nzt+1) = l1d_init(nzt)

       ENDIF
    ENDIF
    l1d      = l1d_init
    l1d_diss = l1d_init
    u1d      = u_init
    u1d_p    = u_init
    v1d      = v_init
    v1d_p    = v_init

!
!-- Set initial horizontal velocities at the lowest grid levels to a very small
!-- value in order to avoid too small time steps caused by the diffusion limit
!-- in the initial phase of a run (at k=1, dz/2 occurs in the limiting formula!)
    u1d(0:1)   = 0.1_wp
    u1d_p(0:1) = 0.1_wp
    v1d(0:1)   = 0.1_wp
    v1d_p(0:1) = 0.1_wp

!
!-- For u*, theta* and the momentum fluxes plausible values are set
    IF ( constant_flux_layer )  THEN
       us1d = 0.1_wp   ! without initial friction the flow would not change
    ELSE
       diss1d(nzb+1) = 0.001_wp
       e1d(nzb+1)  = 1.0_wp
       km1d(nzb+1) = 1.0_wp
       us1d = 0.0_wp
    ENDIF
    ts1d = 0.0_wp
    usws1d = 0.0_wp
    vsws1d = 0.0_wp
    z01d  = roughness_length
    z0h1d = z0h_factor * z01d
    IF ( humidity )  qs1d = 0.0_wp

!
!-- Tendencies must be preset in order to avoid runtime errors
    te_diss  = 0.0_wp
    te_dissm = 0.0_wp
    te_e  = 0.0_wp
    te_em = 0.0_wp
    te_um = 0.0_wp
    te_vm = 0.0_wp

!
!-- Set model constant
    IF ( dissipation_1d == 'as_in_3d_model' )  c_0 = 0.1_wp
    c_mu = c_0**4

!
!-- Set start time in hh:mm:ss - format
    simulated_time_chr = time_to_string( simulated_time_1d )

!
!-- Integrate the 1D-model equations using the Runge-Kutta scheme
    CALL time_integration_1d


 END SUBROUTINE init_1d_model



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Runge-Kutta time differencing scheme for the 1D-model.
!------------------------------------------------------------------------------!

 SUBROUTINE time_integration_1d

    IMPLICIT NONE

    CHARACTER (LEN=9) ::  time_to_string  !< function to transform time from real to character string

    INTEGER(iwp) ::  k  !< loop index

    REAL(wp) ::  a            !< auxiliary variable
    REAL(wp) ::  b            !< auxiliary variable
    REAL(wp) ::  dpt_dz       !< vertical temperature gradient
    REAL(wp) ::  flux         !< vertical temperature gradient
    REAL(wp) ::  kmzm         !< Km(z-dz/2)
    REAL(wp) ::  kmzp         !< Km(z+dz/2)
    REAL(wp) ::  l_stable     !< mixing length for stable case
    REAL(wp) ::  pt_0         !< reference temperature
    REAL(wp) ::  uv_total     !< horizontal wind speed

!
!-- Determine the time step at the start of a 1D-simulation and
!-- determine and printout quantities used for run control
    dt_1d = 0.01_wp
    CALL run_control_1d

!
!-- Start of time loop
    DO  WHILE ( simulated_time_1d < end_time_1d  .AND.  .NOT. stop_dt_1d )

!
!--    Depending on the timestep scheme, carry out one or more intermediate
!--    timesteps

       intermediate_timestep_count = 0
       DO  WHILE ( intermediate_timestep_count < &
                   intermediate_timestep_count_max )

          intermediate_timestep_count = intermediate_timestep_count + 1

          CALL timestep_scheme_steering

!
!--       Compute all tendency terms. If a constant-flux layer is simulated,
!--       k starts at nzb+2.
          DO  k = nzb_diff, nzt

             kmzm = 0.5_wp * ( km1d(k-1) + km1d(k) )
             kmzp = 0.5_wp * ( km1d(k) + km1d(k+1) )
!
!--          u-component
             te_u(k) =  f * ( v1d(k) - vg(k) ) + ( &
                              kmzp * ( u1d(k+1) - u1d(k) ) * ddzu(k+1) &
                            - kmzm * ( u1d(k) - u1d(k-1) ) * ddzu(k)   &
                                                 ) * ddzw(k)
!
!--          v-component
             te_v(k) = -f * ( u1d(k) - ug(k) ) + (                     &
                              kmzp * ( v1d(k+1) - v1d(k) ) * ddzu(k+1) &
                            - kmzm * ( v1d(k) - v1d(k-1) ) * ddzu(k)   &
                                                 ) * ddzw(k)
          ENDDO
          IF ( .NOT. constant_diffusion )  THEN
             DO  k = nzb_diff, nzt
!
!--             TKE and dissipation rate
                kmzm = 0.5_wp * ( km1d(k-1) + km1d(k) )
                kmzp = 0.5_wp * ( km1d(k) + km1d(k+1) )
                IF ( .NOT. humidity )  THEN
                   pt_0 = pt_init(k)
                   flux =  ( pt_init(k+1)-pt_init(k-1) ) * dd2zu(k)
                ELSE
                   pt_0 = pt_init(k) * ( 1.0_wp + 0.61_wp * q_init(k) )
                   flux = ( ( pt_init(k+1) - pt_init(k-1) ) +                  &
                            0.61_wp * ( pt_init(k+1) * q_init(k+1) -           &
                                        pt_init(k-1) * q_init(k-1)   )         &
                          ) * dd2zu(k)
                ENDIF

!
!--             Calculate dissipation rate if no prognostic equation is used for
!--             dissipation rate
                IF ( dissipation_1d == 'detering' )  THEN
                   diss1d(k) = c_0**3 * e1d(k) * SQRT( e1d(k) ) / l1d_diss(k)
                ELSEIF ( dissipation_1d == 'as_in_3d_model' )  THEN
                   diss1d(k) = ( 0.19_wp + 0.74_wp * l1d_diss(k) / l1d_init(k) &
                               ) * e1d(k) * SQRT( e1d(k) ) / l1d_diss(k)
                ENDIF
!
!--             TKE
                te_e(k) = km1d(k) * ( ( ( u1d(k+1) - u1d(k-1) ) * dd2zu(k) )**2&
                                    + ( ( v1d(k+1) - v1d(k-1) ) * dd2zu(k) )**2&
                                    )                                          &
                                    - g / pt_0 * kh1d(k) * flux                &
                                    +            (                             &
                                     kmzp * ( e1d(k+1) - e1d(k) ) * ddzu(k+1)  &
                                   - kmzm * ( e1d(k) - e1d(k-1) ) * ddzu(k)    &
                                                 ) * ddzw(k) / sig_e           &
                                   - diss1d(k)

                IF ( dissipation_1d == 'prognostic' )  THEN
!
!--                dissipation rate
                   IF ( rif1d(k) >= 0.0_wp )  THEN
                      alpha_buoyancy = 1.0_wp - l1d(k) / lambda
                   ELSE
                      alpha_buoyancy = 1.0_wp - ( 1.0_wp + ( c_2 - 1.0_wp )    &
                                                         / ( c_2 - c_1    ) )  &
                                              * l1d(k) / lambda
                   ENDIF
                   c_3 = ( c_1 - c_2 ) * alpha_buoyancy + 1.0_wp
                   te_diss(k) = ( km1d(k) *                                    &
                                  ( ( ( u1d(k+1) - u1d(k-1) ) * dd2zu(k) )**2  &
                                  + ( ( v1d(k+1) - v1d(k-1) ) * dd2zu(k) )**2  &
                                  ) * ( c_1 + (c_2 - c_1) * l1d(k) / lambda )  &
                                  - g / pt_0 * kh1d(k) * flux * c_3            &
                                  - c_2 * diss1d(k)                            &
                                ) * diss1d(k) / ( e1d(k) + 1.0E-20_wp )        &
                                + (   kmzp * ( diss1d(k+1) - diss1d(k) )       &
                                           * ddzu(k+1)                         &
                                    - kmzm * ( diss1d(k) - diss1d(k-1) )       &
                                           * ddzu(k)                           &
                                  ) * ddzw(k) / sig_diss

                ENDIF

             ENDDO
          ENDIF

!
!--       Tendency terms at the top of the constant-flux layer.
!--       Finite differences of the momentum fluxes are computed using half the
!--       normal grid length (2.0*ddzw(k)) for the sake of enhanced accuracy
          IF ( constant_flux_layer )  THEN

             k = nzb+1
             kmzm = 0.5_wp * ( km1d(k-1) + km1d(k) )
             kmzp = 0.5_wp * ( km1d(k) + km1d(k+1) )
             IF ( .NOT. humidity )  THEN
                pt_0 = pt_init(k)
                flux =  ( pt_init(k+1)-pt_init(k-1) ) * dd2zu(k)
             ELSE
                pt_0 = pt_init(k) * ( 1.0_wp + 0.61_wp * q_init(k) )
                flux = ( ( pt_init(k+1) - pt_init(k-1) ) +                     &
                         0.61_wp * ( pt_init(k+1) * q_init(k+1) -              &
                                     pt_init(k-1) * q_init(k-1)   )            &
                       ) * dd2zu(k)
             ENDIF

!
!--          Calculate dissipation rate if no prognostic equation is used for
!--          dissipation rate
             IF ( dissipation_1d == 'detering' )  THEN
                diss1d(k) = c_0**3 * e1d(k) * SQRT( e1d(k) ) / l1d_diss(k)
             ELSEIF ( dissipation_1d == 'as_in_3d_model' )  THEN
                diss1d(k) = ( 0.19_wp + 0.74_wp * l1d_diss(k) / l1d_init(k) )  &
                            * e1d(k) * SQRT( e1d(k) ) / l1d_diss(k)
             ENDIF

!
!--          u-component
             te_u(k) = f * ( v1d(k) - vg(k) ) + (                              &
                       kmzp * ( u1d(k+1) - u1d(k) ) * ddzu(k+1) + usws1d       &
                                                ) * 2.0_wp * ddzw(k)
!
!--          v-component
             te_v(k) = -f * ( u1d(k) - ug(k) ) + (                             &
                       kmzp * ( v1d(k+1) - v1d(k) ) * ddzu(k+1) + vsws1d       &
                                                 ) * 2.0_wp * ddzw(k)
!
!--          TKE
             IF ( .NOT. dissipation_1d == 'prognostic' )  THEN
                !> @query why integrate over 2dz
                !>   Why is it allowed to integrate over two delta-z for e
                !>   while for u and v it is not?
                !>   2018-04-23, gronemeier
                te_e(k) = km1d(k) * ( ( ( u1d(k+1) - u1d(k-1) ) * dd2zu(k) )**2&
                                    + ( ( v1d(k+1) - v1d(k-1) ) * dd2zu(k) )**2&
                                    )                                          &
                                    - g / pt_0 * kh1d(k) * flux                &
                                    +           (                              &
                                     kmzp * ( e1d(k+1) - e1d(k) ) * ddzu(k+1)  &
                                   - kmzm * ( e1d(k) - e1d(k-1) ) * ddzu(k)    &
                                                 ) * ddzw(k) / sig_e           &
                                   - diss1d(k)
             ENDIF

          ENDIF

!
!--       Prognostic equations for all 1D variables
          DO  k = nzb+1, nzt

             u1d_p(k) = u1d(k) + dt_1d * ( tsc(2) * te_u(k) + &
                                           tsc(3) * te_um(k) )
             v1d_p(k) = v1d(k) + dt_1d * ( tsc(2) * te_v(k) + &
                                           tsc(3) * te_vm(k) )

          ENDDO
          IF ( .NOT. constant_diffusion )  THEN

             DO  k = nzb+1, nzt
                e1d_p(k) = e1d(k) + dt_1d * ( tsc(2) * te_e(k) + &
                                              tsc(3) * te_em(k) )
             ENDDO

!
!--          Eliminate negative TKE values, which can result from the
!--          integration due to numerical inaccuracies. In such cases the TKE
!--          value is reduced to 10 percent of its old value.
             WHERE ( e1d_p < 0.0_wp )  e1d_p = 0.1_wp * e1d

             IF ( dissipation_1d == 'prognostic' )  THEN
                DO  k = nzb+1, nzt
                   diss1d_p(k) = diss1d(k) + dt_1d * ( tsc(2) * te_diss(k) + &
                                                       tsc(3) * te_dissm(k) )
                ENDDO
                WHERE ( diss1d_p < 0.0_wp )  diss1d_p = 0.1_wp * diss1d
             ENDIF
          ENDIF

!
!--       Calculate tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' ) THEN
             IF ( intermediate_timestep_count == 1 )  THEN

                DO  k = nzb+1, nzt
                   te_um(k) = te_u(k)
                   te_vm(k) = te_v(k)
                ENDDO

                IF ( .NOT. constant_diffusion )  THEN
                   DO k = nzb+1, nzt
                      te_em(k) = te_e(k)
                   ENDDO
                   IF ( dissipation_1d == 'prognostic' )  THEN
                      DO k = nzb+1, nzt
                         te_dissm(k) = te_diss(k)
                      ENDDO
                   ENDIF
                ENDIF

             ELSEIF ( intermediate_timestep_count < &
                         intermediate_timestep_count_max )  THEN

                DO  k = nzb+1, nzt
                   te_um(k) = -9.5625_wp * te_u(k) + 5.3125_wp * te_um(k)
                   te_vm(k) = -9.5625_wp * te_v(k) + 5.3125_wp * te_vm(k)
                ENDDO

                IF ( .NOT. constant_diffusion )  THEN
                   DO k = nzb+1, nzt
                      te_em(k) = -9.5625_wp * te_e(k) + 5.3125_wp * te_em(k)
                   ENDDO
                   IF ( dissipation_1d == 'prognostic' )  THEN
                      DO k = nzb+1, nzt
                         te_dissm(k) = -9.5625_wp * te_diss(k)  &
                                     +  5.3125_wp * te_dissm(k)
                      ENDDO
                   ENDIF
                ENDIF

             ENDIF
          ENDIF

!
!--       Boundary conditions for the prognostic variables.
!--       At the top boundary (nzt+1) u, v, e, and diss keep their initial
!--       values (ug(nzt+1), vg(nzt+1), 0, 0).
!--       At the bottom boundary, Dirichlet condition is used for u and v (0)
!--       and Neumann condition for e and diss (e(nzb)=e(nzb+1)).
          u1d_p(nzb) = 0.0_wp
          v1d_p(nzb) = 0.0_wp

!
!--       Swap the time levels in preparation for the next time step.
          u1d  = u1d_p
          v1d  = v1d_p
          IF ( .NOT. constant_diffusion )  THEN
             e1d  = e1d_p
             IF ( dissipation_1d == 'prognostic' )  THEN
                diss1d = diss1d_p
             ENDIF
          ENDIF

!
!--       Compute diffusion quantities
          IF ( .NOT. constant_diffusion )  THEN

!
!--          First compute the vertical fluxes in the constant-flux layer
             IF ( constant_flux_layer )  THEN
!
!--             Compute theta* using Rif numbers of the previous time step
                IF ( rif1d(nzb+1) >= 0.0_wp )  THEN
!
!--                Stable stratification
                   ts1d = kappa * ( pt_init(nzb+1) - pt_init(nzb) ) /          &
                          ( LOG( zu(nzb+1) / z0h1d ) + 5.0_wp * rif1d(nzb+1) * &
                                          ( zu(nzb+1) - z0h1d ) / zu(nzb+1)    &
                          )
                ELSE
!
!--                Unstable stratification
                   a = SQRT( 1.0_wp - 16.0_wp * rif1d(nzb+1) )
                   b = SQRT( 1.0_wp - 16.0_wp * rif1d(nzb+1) /                 &
                       zu(nzb+1) * z0h1d )

                   ts1d = kappa * ( pt_init(nzb+1) - pt_init(nzb) ) /          &
                          LOG( (a-1.0_wp) / (a+1.0_wp) *                       &
                               (b+1.0_wp) / (b-1.0_wp) )
                ENDIF

             ENDIF    ! constant_flux_layer
             !> @todo combine if clauses
             !>   The previous and following if clauses can be combined into a
             !>   single clause
             !>   2018-04-23, gronemeier
!
!--          Compute the Richardson-flux numbers,
!--          first at the top of the constant-flux layer using u* of the
!--          previous time step (+1E-30, if u* = 0), then in the remaining area.
!--          There the rif-numbers of the previous time step are used.

             IF ( constant_flux_layer )  THEN
                IF ( .NOT. humidity )  THEN
                   pt_0 = pt_init(nzb+1)
                   flux = ts1d
                ELSE
                   pt_0 = pt_init(nzb+1) * ( 1.0_wp + 0.61_wp * q_init(nzb+1) )
                   flux = ts1d + 0.61_wp * pt_init(k) * qs1d
                ENDIF
                rif1d(nzb+1) = zu(nzb+1) * kappa * g * flux / &
                               ( pt_0 * ( us1d**2 + 1E-30_wp ) )
             ENDIF

             DO  k = nzb_diff, nzt
                IF ( .NOT. humidity )  THEN
                   pt_0 = pt_init(k)
                   flux = ( pt_init(k+1) - pt_init(k-1) ) * dd2zu(k)
                ELSE
                   pt_0 = pt_init(k) * ( 1.0_wp + 0.61_wp * q_init(k) )
                   flux = ( ( pt_init(k+1) - pt_init(k-1) )                    &
                            + 0.61_wp                                          &
                            * (   pt_init(k+1) * q_init(k+1)                   &
                                - pt_init(k-1) * q_init(k-1) )                 &
                          ) * dd2zu(k)
                ENDIF
                IF ( rif1d(k) >= 0.0_wp )  THEN
                   rif1d(k) = g / pt_0 * flux /                                &
                              (  ( ( u1d(k+1) - u1d(k-1) ) * dd2zu(k) )**2     &
                               + ( ( v1d(k+1) - v1d(k-1) ) * dd2zu(k) )**2     &
                               + 1E-30_wp                                      &
                              )
                ELSE
                   rif1d(k) = g / pt_0 * flux /                                &
                              (  ( ( u1d(k+1) - u1d(k-1) ) * dd2zu(k) )**2     &
                               + ( ( v1d(k+1) - v1d(k-1) ) * dd2zu(k) )**2     &
                               + 1E-30_wp                                      &
                              ) * ( 1.0_wp - 16.0_wp * rif1d(k) )**0.25_wp
                ENDIF
             ENDDO
!
!--          Richardson-numbers must remain restricted to a realistic value
!--          range. It is exceeded excessively for very small velocities
!--          (u,v --> 0).
             WHERE ( rif1d < -5.0_wp )  rif1d = -5.0_wp
             WHERE ( rif1d > 1.0_wp )  rif1d = 1.0_wp

!
!--          Compute u* from the absolute velocity value
             IF ( constant_flux_layer )  THEN
                uv_total = SQRT( u1d(nzb+1)**2 + v1d(nzb+1)**2 )

                IF ( rif1d(nzb+1) >= 0.0_wp )  THEN
!
!--                Stable stratification
                   us1d = kappa * uv_total / (                                 &
                             LOG( zu(nzb+1) / z01d ) + 5.0_wp * rif1d(nzb+1) * &
                                              ( zu(nzb+1) - z01d ) / zu(nzb+1) &
                                             )
                ELSE
!
!--                Unstable stratification
                   a = 1.0_wp / SQRT( SQRT( 1.0_wp - 16.0_wp * rif1d(nzb+1) ) )
                   b = 1.0_wp / SQRT( SQRT( 1.0_wp - 16.0_wp * rif1d(nzb+1) /  &
                                                     zu(nzb+1) * z01d ) )
                   us1d = kappa * uv_total / (                                 &
                              LOG( (1.0_wp+b) / (1.0_wp-b) * (1.0_wp-a) /      &
                                   (1.0_wp+a) ) +                              &
                              2.0_wp * ( ATAN( b ) - ATAN( a ) )               &
                                             )
                ENDIF

!
!--             Compute the momentum fluxes for the diffusion terms
                usws1d  = - u1d(nzb+1) / uv_total * us1d**2
                vsws1d  = - v1d(nzb+1) / uv_total * us1d**2

!
!--             Boundary condition for the turbulent kinetic energy and
!--             dissipation rate at the top of the constant-flux layer.
!--             Additional Neumann condition de/dz = 0 at nzb is set to ensure
!--             compatibility with the 3D model.
                IF ( ibc_e_b == 2 )  THEN
                   e1d(nzb+1) = ( us1d / c_0 )**2
                ENDIF
                IF ( dissipation_1d == 'prognostic' )  THEN
                   e1d(nzb+1) = ( us1d / c_0 )**2
                   diss1d(nzb+1) = us1d**3 / ( kappa * zu(nzb+1) )
                   diss1d(nzb) = diss1d(nzb+1)
                ENDIF
                e1d(nzb) = e1d(nzb+1)

                IF ( humidity ) THEN
!
!--                Compute q*
                   IF ( rif1d(nzb+1) >= 0.0_wp )  THEN
!
!--                   Stable stratification
                      qs1d = kappa * ( q_init(nzb+1) - q_init(nzb) ) /         &
                          ( LOG( zu(nzb+1) / z0h1d ) + 5.0_wp * rif1d(nzb+1) * &
                                          ( zu(nzb+1) - z0h1d ) / zu(nzb+1)    &
                          )
                   ELSE
!
!--                   Unstable stratification
                      a = SQRT( 1.0_wp - 16.0_wp * rif1d(nzb+1) )
                      b = SQRT( 1.0_wp - 16.0_wp * rif1d(nzb+1) /              &
                                         zu(nzb+1) * z0h1d )
                      qs1d = kappa * ( q_init(nzb+1) - q_init(nzb) ) /         &
                             LOG( (a-1.0_wp) / (a+1.0_wp) *                    &
                                  (b+1.0_wp) / (b-1.0_wp) )
                   ENDIF
                ELSE
                   qs1d = 0.0_wp
                ENDIF

             ENDIF   !  constant_flux_layer

!
!--          Compute the diabatic mixing length. The unstable stratification
!--          must not be considered for l1d (km1d) as it is already considered
!--          in the dissipation of TKE via l1d_diss. Otherwise, km1d would be
!--          too large.
             IF ( dissipation_1d /= 'prognostic' )  THEN
                IF ( mixing_length_1d == 'blackadar' )  THEN
                   DO  k = nzb+1, nzt
                      IF ( rif1d(k) >= 0.0_wp )  THEN
                         l1d(k) = l1d_init(k) / ( 1.0_wp + 5.0_wp * rif1d(k) )
                         l1d_diss(k) = l1d(k)
                      ELSE
                         l1d(k) = l1d_init(k)
                         l1d_diss(k) = l1d_init(k) *                           &
                                       SQRT( 1.0_wp - 16.0_wp * rif1d(k) )
                      ENDIF
                   ENDDO
                ELSEIF ( mixing_length_1d == 'as_in_3d_model' )  THEN
                   DO  k = nzb+1, nzt
                      dpt_dz = ( pt_init(k+1) - pt_init(k-1) ) * dd2zu(k)
                      IF ( dpt_dz > 0.0_wp )  THEN
                         l_stable = 0.76_wp * SQRT( e1d(k) )                   &
                                  / SQRT( g / pt_init(k) * dpt_dz ) + 1E-5_wp
                      ELSE
                         l_stable = l1d_init(k)
                      ENDIF
                      l1d(k) = MIN( l1d_init(k), l_stable )
                      l1d_diss(k) = l1d(k)
                   ENDDO
                ENDIF
             ELSE
                DO  k = nzb+1, nzt
                   l1d(k) = c_0**3 * e1d(k) * SQRT( e1d(k) )                   &
                          / ( diss1d(k) + 1.0E-30_wp )
                ENDDO
             ENDIF

!
!--          Compute the diffusion coefficients for momentum via the
!--          corresponding Prandtl-layer relationship and according to
!--          Prandtl-Kolmogorov, respectively
             IF ( constant_flux_layer )  THEN
                IF ( rif1d(nzb+1) >= 0.0_wp )  THEN
                   km1d(nzb+1) = us1d * kappa * zu(nzb+1) /                    &
                                 ( 1.0_wp + 5.0_wp * rif1d(nzb+1) )
                ELSE
                   km1d(nzb+1) = us1d * kappa * zu(nzb+1) *                    &
                                 ( 1.0_wp - 16.0_wp * rif1d(nzb+1) )**0.25_wp
                ENDIF
             ENDIF

             IF ( dissipation_1d == 'prognostic' )  THEN
                DO  k = nzb_diff, nzt
                   km1d(k) = c_mu * e1d(k)**2 / ( diss1d(k) + 1.0E-30_wp )
                ENDDO
             ELSE
                DO  k = nzb_diff, nzt
                   km1d(k) = c_0 * SQRT( e1d(k) ) * l1d(k)
                ENDDO
             ENDIF

!
!--          Add damping layer
             DO  k = damp_level_ind_1d+1, nzt+1
                km1d(k) = 1.1_wp * km1d(k-1)
                km1d(k) = MIN( km1d(k), 10.0_wp )
             ENDDO

!
!--          Compute the diffusion coefficient for heat via the relationship
!--          kh = phim / phih * km
             DO  k = nzb+1, nzt
                IF ( rif1d(k) >= 0.0_wp )  THEN
                   kh1d(k) = km1d(k)
                ELSE
                   kh1d(k) = km1d(k) * ( 1.0_wp - 16.0_wp * rif1d(k) )**0.25_wp
                ENDIF
             ENDDO

          ENDIF   ! .NOT. constant_diffusion

       ENDDO   ! intermediate step loop

!
!--    Increment simulated time and output times
       current_timestep_number_1d = current_timestep_number_1d + 1
       simulated_time_1d          = simulated_time_1d + dt_1d
       simulated_time_chr         = time_to_string( simulated_time_1d )
       time_pr_1d                 = time_pr_1d          + dt_1d
       time_run_control_1d        = time_run_control_1d + dt_1d

!
!--    Determine and print out quantities for run control
       IF ( time_run_control_1d >= dt_run_control_1d )  THEN
          CALL run_control_1d
          time_run_control_1d = time_run_control_1d - dt_run_control_1d
       ENDIF

!
!--    Profile output on file
       IF ( time_pr_1d >= dt_pr_1d )  THEN
          CALL print_1d_model
          time_pr_1d = time_pr_1d - dt_pr_1d
       ENDIF

!
!--    Determine size of next time step
       CALL timestep_1d

    ENDDO   ! time loop
!
!-- Set intermediate_timestep_count back to zero. This is required e.g. for
!-- initial calls of calc_mean_profile.
    intermediate_timestep_count = 0

 END SUBROUTINE time_integration_1d


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute and print out quantities for run control of the 1D model.
!------------------------------------------------------------------------------!

 SUBROUTINE run_control_1d


    IMPLICIT NONE

    INTEGER(iwp) ::  k     !< loop index

    REAL(wp) ::  alpha     !< angle of wind vector at top of constant-flux layer
    REAL(wp) ::  energy    !< kinetic energy
    REAL(wp) ::  umax      !< maximum of u
    REAL(wp) ::  uv_total  !< horizontal wind speed
    REAL(wp) ::  vmax      !< maximum of v

!
!-- Output
    IF ( myid == 0 )  THEN
!
!--    If necessary, write header
       IF ( .NOT. run_control_header_1d )  THEN
          CALL check_open( 15 )
          WRITE ( 15, 100 )
          run_control_header_1d = .TRUE.
       ENDIF

!
!--    Compute control quantities
!--    grid level nzp is excluded due to mirror boundary condition
       umax = 0.0_wp; vmax = 0.0_wp; energy = 0.0_wp
       DO  k = nzb+1, nzt+1
          umax = MAX( ABS( umax ), ABS( u1d(k) ) )
          vmax = MAX( ABS( vmax ), ABS( v1d(k) ) )
          energy = energy + 0.5_wp * ( u1d(k)**2 + v1d(k)**2 )
       ENDDO
       energy = energy / REAL( nzt - nzb + 1, KIND=wp )

       uv_total = SQRT( u1d(nzb+1)**2 + v1d(nzb+1)**2 )
       IF ( ABS( v1d(nzb+1) ) < 1.0E-5_wp )  THEN
          alpha = ACOS( SIGN( 1.0_wp , u1d(nzb+1) ) )
       ELSE
          alpha = ACOS( u1d(nzb+1) / uv_total )
          IF ( v1d(nzb+1) <= 0.0_wp )  alpha = 2.0_wp * pi - alpha
       ENDIF
       alpha = alpha / ( 2.0_wp * pi ) * 360.0_wp

       WRITE ( 15, 101 )  current_timestep_number_1d, simulated_time_chr, &
                          dt_1d, umax, vmax, us1d, alpha, energy
!
!--    Write buffer contents to disc immediately
       FLUSH( 15 )

    ENDIF

!
!-- formats
100 FORMAT (///'1D run control output:'/ &
              &'------------------------------'// &
           &'ITER.   HH:MM:SS    DT      UMAX   VMAX    U*   ALPHA   ENERG.'/ &
           &'-------------------------------------------------------------')
101 FORMAT (I7,1X,A9,1X,F6.2,2X,F6.2,1X,F6.2,1X,F6.3,2X,F5.1,2X,F7.2)


 END SUBROUTINE run_control_1d



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the time step w.r.t. the diffusion criterion
!------------------------------------------------------------------------------!

 SUBROUTINE timestep_1d

    IMPLICIT NONE

    INTEGER(iwp) ::  k    !< loop index

    REAL(wp) ::  dt_diff  !< time step accorind to diffusion criterion
    REAL(wp) ::  dt_old   !< previous time step
    REAL(wp) ::  fac      !< factor of criterion
    REAL(wp) ::  value    !< auxiliary variable

!
!-- Save previous time step
    dt_old = dt_1d

!
!-- Compute the currently feasible time step according to the diffusion
!-- criterion. At nzb+1 the half grid length is used.
    fac = 0.125
    dt_diff = dt_max_1d
    DO  k = nzb+2, nzt
       value   = fac * dzu(k) * dzu(k) / ( km1d(k) + 1E-20_wp )
       dt_diff = MIN( value, dt_diff )
    ENDDO
    value   = fac * zu(nzb+1) * zu(nzb+1) / ( km1d(nzb+1) + 1E-20_wp )
    dt_1d = MIN( value, dt_diff )

!
!-- Limit the new time step to a maximum of 10 times the previous time step
    dt_1d = MIN( dt_old * 10.0_wp, dt_1d )

!
!-- Set flag when the time step becomes too small
    IF ( dt_1d < ( 1.0E-15_wp * dt_max_1d ) )  THEN
       stop_dt_1d = .TRUE.

       WRITE( message_string, * ) 'timestep has exceeded the lower limit&',    &
                                  'dt_1d = ',dt_1d,' s   simulation stopped!'
       CALL message( 'timestep_1d', 'PA0192', 1, 2, 0, 6, 0 )

    ENDIF

 END SUBROUTINE timestep_1d



!------------------------------------------------------------------------------!
! Description:
! ------------
!> List output of profiles from the 1D-model
!------------------------------------------------------------------------------!

 SUBROUTINE print_1d_model

    IMPLICIT NONE

    INTEGER(iwp) ::  k  !< loop parameter

    LOGICAL, SAVE :: write_first = .TRUE. !< flag for writing header


    IF ( myid == 0 )  THEN
!
!--    Open list output file for profiles from the 1D-model
       CALL check_open( 17 )

!
!--    Write Header
       IF ( write_first )  THEN
          WRITE ( 17, 100 )  TRIM( run_description_header )
          write_first = .FALSE.
       ENDIF

!
!--    Write the values
       WRITE ( 17, 104 )  TRIM( simulated_time_chr )
       WRITE ( 17, 101 )
       WRITE ( 17, 102 )
       WRITE ( 17, 101 )
       DO  k = nzt+1, nzb, -1
          WRITE ( 17, 103)  k, zu(k), u1d(k), v1d(k), pt_init(k), e1d(k), &
                            rif1d(k), km1d(k), kh1d(k), l1d(k), diss1d(k)
       ENDDO
       WRITE ( 17, 101 )
       WRITE ( 17, 102 )
       WRITE ( 17, 101 )

!
!--    Write buffer contents to disc immediately
       FLUSH( 17 )

    ENDIF

!
!-- Formats
100 FORMAT ('# ',A/'#',10('-')/'# 1d-model profiles')
104 FORMAT (//'# Time: ',A)
101 FORMAT ('#',111('-'))
102 FORMAT ('#  k     zu      u          v          pt         e          ',   &
            'rif        Km         Kh         l          diss   ')
103 FORMAT (1X,I4,1X,F7.1,9(1X,E10.3))


 END SUBROUTINE print_1d_model


 END MODULE
