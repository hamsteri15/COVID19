!> @file prognostic_equations.f90
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
! $Id: prognostic_equations.f90 4370 2020-01-10 14:00:44Z raasch $
! vector directives added to force vectorization on Intel19 compiler
! 
! 4360 2020-01-07 11:25:50Z suehring
! Introduction of wall_flags_total_0, which currently sets bits based on static
! topography information used in wall_flags_static_0
! 
! 4329 2019-12-10 15:46:36Z motisi
! Renamed wall_flags_0 to wall_flags_static_0
! 
! 4182 2019-08-22 15:20:23Z scharf
! Corrected "Former revisions" section
! 
! 4110 2019-07-22 17:05:21Z suehring
! pass integer flag array to WS scalar advection routine which is now necessary
! as the flags may differ for scalars, e.g. pt can be cyclic while chemical
! species may be non-cyclic. Further, pass boundary flags. 
! 
! 4109 2019-07-22 17:00:34Z suehring
! Application of monotonic flux limiter for the vertical scalar advection 
! up to the topography top (only for the cache-optimized version at the 
! moment). Please note, at the moment the limiter is only applied for passive
! scalars.
! 
! 4048 2019-06-21 21:00:21Z knoop
! Moved tcm_prognostic_equations to module_interface
! 
! 3987 2019-05-22 09:52:13Z kanani
! Introduce alternative switch for debug output during timestepping
! 
! 3956 2019-05-07 12:32:52Z monakurppa
! Removed salsa calls.
! 
! 3931 2019-04-24 16:34:28Z schwenkel
! Correct/complete module_interface introduction for chemistry model
! 
! 3899 2019-04-16 14:05:27Z monakurppa
! Corrections in the OpenMP version of salsa
!
! 3887 2019 -04-12 08:47:41Z schwenkel
! Implicit Bugfix for chemistry model, loop for non_transport_physics over 
! ghost points is avoided. Instead introducing module_interface_exchange_horiz.
! 
! 3885 2019-04-11 11:29:34Z kanani
! Changes related to global restructuring of location messages and introduction 
! of additional debug messages
! 
! 3881 2019-04-10 09:31:22Z suehring
! Bugfix in OpenMP directive
! 
! 3880 2019-04-08 21:43:02Z knoop
! Moved wtm_tendencies to module_interface_actions
! 
! 3874 2019-04-08 16:53:48Z knoop
! Added non_transport_physics module interfaces and moved bcm code into it
! 
! 3872 2019-04-08 15:03:06Z knoop
! Moving prognostic equations of bcm into bulk_cloud_model_mod
! 
! 3864 2019-04-05 09:01:56Z monakurppa
! Modifications made for salsa:
! - salsa_prognostic_equations moved to salsa_mod (and the call to 
!   module_interface_mod)
! - Renamed nbins --> nbins_aerosol, ncc_tot --> ncomponents_mass and 
!   ngast --> ngases_salsa and loop indices b, c and sg to ib, ic and ig
! 
! 3840 2019-03-29 10:35:52Z knoop
! added USE chem_gasphase_mod for nspec, nspec and spc_names
! 
! 3820 2019-03-27 11:53:41Z forkel
! renamed do_depo to deposition_dry (ecc)
! 
! 3797 2019-03-15 11:15:38Z forkel
! Call chem_integegrate in OpenMP loop   (ketelsen)
! 
! 
! 3771 2019-02-28 12:19:33Z raasch
! preprocessor directivs fro rrtmg added
! 
! 3761 2019-02-25 15:31:42Z raasch
! unused variable removed
! 
! 3719 2019-02-06 13:10:18Z kanani
! Cleaned up chemistry cpu measurements
! 
! 3684 2019-01-20 20:20:58Z knoop
! OpenACC port for SPEC
!
! Revision 1.1  2000/04/13 14:56:27  schroeter
! Initial revision
!
!
! Description:
! ------------
!> Solving the prognostic equations.
!------------------------------------------------------------------------------!
 MODULE prognostic_equations_mod

    USE advec_s_bc_mod,                                                        &
        ONLY:  advec_s_bc

    USE advec_s_pw_mod,                                                        &
        ONLY:  advec_s_pw

    USE advec_s_up_mod,                                                        &
        ONLY:  advec_s_up

    USE advec_u_pw_mod,                                                        &
        ONLY:  advec_u_pw

    USE advec_u_up_mod,                                                        &
        ONLY:  advec_u_up

    USE advec_v_pw_mod,                                                        &
        ONLY:  advec_v_pw

    USE advec_v_up_mod,                                                        &
        ONLY:  advec_v_up

    USE advec_w_pw_mod,                                                        &
        ONLY:  advec_w_pw

    USE advec_w_up_mod,                                                        &
        ONLY:  advec_w_up

    USE advec_ws,                                                              &
        ONLY:  advec_s_ws, advec_u_ws, advec_v_ws, advec_w_ws

    USE arrays_3d,                                                             &
        ONLY:  diss_l_e, diss_l_pt, diss_l_q,                                  &
               diss_l_s, diss_l_sa, diss_s_e,                                  &
               diss_s_pt, diss_s_q, diss_s_s, diss_s_sa,                       &
               e, e_p, flux_s_e, flux_s_pt, flux_s_q,                          &
               flux_s_s, flux_s_sa, flux_l_e,                                  &
               flux_l_pt, flux_l_q, flux_l_s,                                  &
               flux_l_sa, pt, ptdf_x, ptdf_y, pt_init,                         &
               pt_p, prho, q, q_init, q_p, rdf, rdf_sc,                        &
               ref_state, rho_ocean, s, s_init, s_p, tend, te_m,               &
               tpt_m, tq_m, ts_m, tu_m, tv_m, tw_m, u,                         &
               ug, u_init, u_p, v, vg, vpt, v_init, v_p, w, w_p

    USE buoyancy_mod,                                                          &
        ONLY:  buoyancy

    USE control_parameters,                                                    &
        ONLY:  bc_dirichlet_l,                                                 &
               bc_dirichlet_n,                                                 &
               bc_dirichlet_r,                                                 &
               bc_dirichlet_s,                                                 &
               bc_radiation_l,                                                 &
               bc_radiation_n,                                                 &
               bc_radiation_r,                                                 &
               bc_radiation_s,                                                 &
               constant_diffusion,                                             &
               debug_output_timestep,                                          &
               dp_external, dp_level_ind_b, dp_smooth_factor, dpdxy, dt_3d,    &
               humidity, intermediate_timestep_count,                          &
               intermediate_timestep_count_max, large_scale_forcing,           &
               large_scale_subsidence,                                         &
               monotonic_limiter_z,                                            &
               neutral, nudging,                                               &
               ocean_mode, passive_scalar, plant_canopy, pt_reference,         &
               scalar_advec, scalar_advec, simulated_time, sloping_surface,    &
               timestep_scheme, tsc, use_subsidence_tendencies,                &
               use_upstream_for_tke, wind_turbine, ws_scheme_mom,              & 
               ws_scheme_sca, urban_surface, land_surface,                     &
               time_since_reference_point, salsa

    USE coriolis_mod,                                                          &
        ONLY:  coriolis

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point, log_point_s

    USE diffusion_s_mod,                                                       &
        ONLY:  diffusion_s

    USE diffusion_u_mod,                                                       &
        ONLY:  diffusion_u

    USE diffusion_v_mod,                                                       &
        ONLY:  diffusion_v

    USE diffusion_w_mod,                                                       &
        ONLY:  diffusion_w

    USE indices,                                                               &
        ONLY:  advc_flags_s,                                                   &
               nbgp, nxl, nxlg, nxlu, nxr, nxrg, nyn, nyng, nys, nysg, nysv,   &
               nzb, nzt, wall_flags_total_0

    USE kinds

    USE lsf_nudging_mod,                                                       &
        ONLY:  ls_advec, nudge

    USE module_interface,                                                      &
        ONLY:  module_interface_actions, &
               module_interface_non_advective_processes, &
               module_interface_exchange_horiz, &
               module_interface_prognostic_equations

    USE ocean_mod,                                                             &
        ONLY:  stokes_drift_terms, stokes_force,   &
               wave_breaking, wave_breaking_term

    USE plant_canopy_model_mod,                                                &
        ONLY:  cthf, pcm_tendency

#if defined( __rrtmg )
    USE radiation_model_mod,                                                   &
        ONLY:  radiation, radiation_tendency,                                  &
               skip_time_do_radiation
#endif

    USE statistics,                                                            &
        ONLY:  hom

    USE subsidence_mod,                                                        &
        ONLY:  subsidence

    USE surface_mod,                                                           &
        ONLY :  surf_def_h, surf_def_v, surf_lsm_h, surf_lsm_v, surf_usm_h,    &
                surf_usm_v

    IMPLICIT NONE

    PRIVATE
    PUBLIC prognostic_equations_cache, prognostic_equations_vector

    INTERFACE prognostic_equations_cache
       MODULE PROCEDURE prognostic_equations_cache
    END INTERFACE prognostic_equations_cache

    INTERFACE prognostic_equations_vector
       MODULE PROCEDURE prognostic_equations_vector
    END INTERFACE prognostic_equations_vector


 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Version with one optimized loop over all equations. It is only allowed to
!> be called for the Wicker and Skamarock or Piascek-Williams advection scheme.
!>
!> Here the calls of most subroutines are embedded in two DO loops over i and j,
!> so communication between CPUs is not allowed (does not make sense) within
!> these loops.
!>
!> (Optimized to avoid cache missings, i.e. for Power4/5-architectures.)
!------------------------------------------------------------------------------!

 SUBROUTINE prognostic_equations_cache


    INTEGER(iwp) ::  i                   !<
    INTEGER(iwp) ::  i_omp_start         !<
    INTEGER(iwp) ::  j                   !<
    INTEGER(iwp) ::  k                   !<
!$  INTEGER(iwp) ::  omp_get_thread_num  !<
    INTEGER(iwp) ::  tn = 0              !<

    LOGICAL      ::  loop_start          !<


    IF ( debug_output_timestep )  CALL debug_message( 'prognostic_equations_cache', 'start' )
!
!-- Time measurement can only be performed for the whole set of equations
    CALL cpu_log( log_point(32), 'all progn.equations', 'start' )

    !$OMP PARALLEL PRIVATE (i,j)
    !$OMP DO
    DO  i = nxl, nxr
       DO  j = nys, nyn
!
!--       Calculate non advective processes for all other modules
          CALL module_interface_non_advective_processes( i, j )
       ENDDO
    ENDDO
!
!-- Module Inferface for exchange horiz after non_advective_processes but before
!-- advection. Therefore, non_advective_processes must not run for ghost points.
    !$OMP END PARALLEL
    CALL module_interface_exchange_horiz()
!
!-- Loop over all prognostic equations
    !$OMP PARALLEL PRIVATE (i,i_omp_start,j,k,loop_start,tn)

    !$  tn = omp_get_thread_num()
    loop_start = .TRUE.

    !$OMP DO
    DO  i = nxl, nxr

!
!--    Store the first loop index. It differs for each thread and is required
!--    later in advec_ws
       IF ( loop_start )  THEN
          loop_start  = .FALSE.
          i_omp_start = i
       ENDIF

       DO  j = nys, nyn
!
!--       Tendency terms for u-velocity component. Please note, in case of 
!--       non-cyclic boundary conditions the grid point i=0 is excluded from
!--       the prognostic equations for the u-component.
          IF ( i >= nxlu )  THEN

             tend(:,j,i) = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( ws_scheme_mom )  THEN
                   CALL advec_u_ws( i, j, i_omp_start, tn )
                ELSE
                   CALL advec_u_pw( i, j )
                ENDIF
             ELSE
                CALL advec_u_up( i, j )
             ENDIF
             CALL diffusion_u( i, j )
             CALL coriolis( i, j, 1 )
             IF ( sloping_surface  .AND.  .NOT. neutral )  THEN
                CALL buoyancy( i, j, pt, 1 )
             ENDIF

!
!--          Drag by plant canopy
             IF ( plant_canopy )  CALL pcm_tendency( i, j, 1 )

!
!--          External pressure gradient
             IF ( dp_external )  THEN
                DO  k = dp_level_ind_b+1, nzt
                   tend(k,j,i) = tend(k,j,i) - dpdxy(1) * dp_smooth_factor(k)
                ENDDO
             ENDIF

!
!--          Nudging
             IF ( nudging )  CALL nudge( i, j, simulated_time, 'u' )

!
!--          Effect of Stokes drift (in ocean mode only)
             IF ( stokes_force )  CALL stokes_drift_terms( i, j, 1 )

             CALL module_interface_actions( i, j, 'u-tendency' )
!
!--          Prognostic equation for u-velocity component
             DO  k = nzb+1, nzt

                u_p(k,j,i) = u(k,j,i) + ( dt_3d *                               &
                                            ( tsc(2) * tend(k,j,i) +            &
                                              tsc(3) * tu_m(k,j,i) )            &
                                            - tsc(5) * rdf(k)                   &
                                                     * ( u(k,j,i) - u_init(k) ) &
                                        ) * MERGE( 1.0_wp, 0.0_wp,              &
                                           BTEST( wall_flags_total_0(k,j,i), 1 )&
                                                 )
             ENDDO

!
!--          Add turbulence generated by wave breaking (in ocean mode only)
             IF ( wave_breaking  .AND.                                         &
               intermediate_timestep_count == intermediate_timestep_count_max )&
             THEN
                CALL wave_breaking_term( i, j, 1 )
             ENDIF

!
!--          Calculate tendencies for the next Runge-Kutta step
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( intermediate_timestep_count == 1 )  THEN
                   DO  k = nzb+1, nzt
                      tu_m(k,j,i) = tend(k,j,i)
                   ENDDO
                ELSEIF ( intermediate_timestep_count < &
                         intermediate_timestep_count_max )  THEN
                   DO  k = nzb+1, nzt
                      tu_m(k,j,i) =   -9.5625_wp * tend(k,j,i)                &
                                     + 5.3125_wp * tu_m(k,j,i)
                   ENDDO
                ENDIF
             ENDIF

          ENDIF
!
!--       Tendency terms for v-velocity component. Please note, in case of 
!--       non-cyclic boundary conditions the grid point j=0 is excluded from
!--       the prognostic equations for the v-component. !--       
          IF ( j >= nysv )  THEN

             tend(:,j,i) = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( ws_scheme_mom )  THEN
                    CALL advec_v_ws( i, j, i_omp_start, tn )
                ELSE
                    CALL advec_v_pw( i, j )
                ENDIF
             ELSE
                CALL advec_v_up( i, j )
             ENDIF
             CALL diffusion_v( i, j )
             CALL coriolis( i, j, 2 )

!
!--          Drag by plant canopy
             IF ( plant_canopy )  CALL pcm_tendency( i, j, 2 )

!
!--          External pressure gradient
             IF ( dp_external )  THEN
                DO  k = dp_level_ind_b+1, nzt
                   tend(k,j,i) = tend(k,j,i) - dpdxy(2) * dp_smooth_factor(k)
                ENDDO
             ENDIF

!
!--          Nudging
             IF ( nudging )  CALL nudge( i, j, simulated_time, 'v' )

!
!--          Effect of Stokes drift (in ocean mode only)
             IF ( stokes_force )  CALL stokes_drift_terms( i, j, 2 )

             CALL module_interface_actions( i, j, 'v-tendency' )
!
!--          Prognostic equation for v-velocity component
             DO  k = nzb+1, nzt
                v_p(k,j,i) = v(k,j,i) + ( dt_3d *                              &
                                            ( tsc(2) * tend(k,j,i) +           &
                                              tsc(3) * tv_m(k,j,i) )           &
                                            - tsc(5) * rdf(k)                  &
                                                     * ( v(k,j,i) - v_init(k) )&
                                        ) * MERGE( 1.0_wp, 0.0_wp,             &
                                          BTEST( wall_flags_total_0(k,j,i), 2 )&
                                                 )
             ENDDO

!
!--          Add turbulence generated by wave breaking (in ocean mode only)
             IF ( wave_breaking  .AND.                                         &
               intermediate_timestep_count == intermediate_timestep_count_max )&
             THEN
                CALL wave_breaking_term( i, j, 2 )
             ENDIF

!
!--          Calculate tendencies for the next Runge-Kutta step
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( intermediate_timestep_count == 1 )  THEN
                   DO  k = nzb+1, nzt
                      tv_m(k,j,i) = tend(k,j,i)
                   ENDDO
                ELSEIF ( intermediate_timestep_count < &
                         intermediate_timestep_count_max )  THEN
                   DO  k = nzb+1, nzt
                      tv_m(k,j,i) =   -9.5625_wp * tend(k,j,i)                 &
                                     + 5.3125_wp * tv_m(k,j,i)
                   ENDDO
                ENDIF
             ENDIF

          ENDIF

!
!--       Tendency terms for w-velocity component
          tend(:,j,i) = 0.0_wp
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( ws_scheme_mom )  THEN
                CALL advec_w_ws( i, j, i_omp_start, tn )
             ELSE
                CALL advec_w_pw( i, j )
             END IF
          ELSE
             CALL advec_w_up( i, j )
          ENDIF
          CALL diffusion_w( i, j )
          CALL coriolis( i, j, 3 )

          IF ( .NOT. neutral )  THEN
             IF ( ocean_mode )  THEN
                CALL buoyancy( i, j, rho_ocean, 3 )
             ELSE
                IF ( .NOT. humidity )  THEN
                   CALL buoyancy( i, j, pt, 3 )
                ELSE
                   CALL buoyancy( i, j, vpt, 3 )
                ENDIF
             ENDIF
          ENDIF

!
!--       Drag by plant canopy
          IF ( plant_canopy )  CALL pcm_tendency( i, j, 3 )

!
!--       Effect of Stokes drift (in ocean mode only)
          IF ( stokes_force )  CALL stokes_drift_terms( i, j, 3 )

          CALL module_interface_actions( i, j, 'w-tendency' )
!
!--       Prognostic equation for w-velocity component
          DO  k = nzb+1, nzt-1
             w_p(k,j,i) = w(k,j,i) + ( dt_3d *                                 &
                                             ( tsc(2) * tend(k,j,i) +          &
                                               tsc(3) * tw_m(k,j,i) )          &
                                             - tsc(5) * rdf(k) * w(k,j,i)      &
                                     ) * MERGE( 1.0_wp, 0.0_wp,                &
                                          BTEST( wall_flags_total_0(k,j,i), 3 )&
                                              )
          ENDDO

!
!--       Calculate tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                DO  k = nzb+1, nzt-1
                   tw_m(k,j,i) = tend(k,j,i)
                ENDDO
             ELSEIF ( intermediate_timestep_count < &
                      intermediate_timestep_count_max )  THEN
                DO  k = nzb+1, nzt-1
                   tw_m(k,j,i) =   -9.5625_wp * tend(k,j,i)                    &
                                  + 5.3125_wp * tw_m(k,j,i)
                ENDDO
             ENDIF
          ENDIF

!
!--       If required, compute prognostic equation for potential temperature
          IF ( .NOT. neutral )  THEN
!
!--          Tendency terms for potential temperature
             tend(:,j,i) = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                   IF ( ws_scheme_sca )  THEN
                      CALL advec_s_ws( advc_flags_s,                           &
                                       i, j, pt, 'pt', flux_s_pt, diss_s_pt,   &
                                       flux_l_pt, diss_l_pt, i_omp_start, tn,  &
                                       bc_dirichlet_l  .OR.  bc_radiation_l,   &
                                       bc_dirichlet_n  .OR.  bc_radiation_n,   &
                                       bc_dirichlet_r  .OR.  bc_radiation_r,   &
                                       bc_dirichlet_s  .OR.  bc_radiation_s )
                   ELSE
                      CALL advec_s_pw( i, j, pt )
                   ENDIF
             ELSE
                CALL advec_s_up( i, j, pt )
             ENDIF
             CALL diffusion_s( i, j, pt,                                       &
                               surf_def_h(0)%shf, surf_def_h(1)%shf,           &
                               surf_def_h(2)%shf,                              &
                               surf_lsm_h%shf,    surf_usm_h%shf,              &
                               surf_def_v(0)%shf, surf_def_v(1)%shf,           &
                               surf_def_v(2)%shf, surf_def_v(3)%shf,           &
                               surf_lsm_v(0)%shf, surf_lsm_v(1)%shf,           &
                               surf_lsm_v(2)%shf, surf_lsm_v(3)%shf,           &
                               surf_usm_v(0)%shf, surf_usm_v(1)%shf,           &
                               surf_usm_v(2)%shf, surf_usm_v(3)%shf )

!
!--          Consideration of heat sources within the plant canopy
             IF ( plant_canopy  .AND.                                          &
                (cthf /= 0.0_wp  .OR. urban_surface  .OR.  land_surface) )  THEN
                CALL pcm_tendency( i, j, 4 )
             ENDIF

!
!--          Large scale advection
             IF ( large_scale_forcing )  THEN
                CALL ls_advec( i, j, simulated_time, 'pt' )
             ENDIF

!
!--          Nudging
             IF ( nudging )  CALL nudge( i, j, simulated_time, 'pt' )

!
!--          If required, compute effect of large-scale subsidence/ascent
             IF ( large_scale_subsidence  .AND.                                &
                  .NOT. use_subsidence_tendencies )  THEN
                CALL subsidence( i, j, tend, pt, pt_init, 2 )
             ENDIF

#if defined( __rrtmg )
!
!--          If required, add tendency due to radiative heating/cooling
             IF ( radiation  .AND.                                             &
                  simulated_time > skip_time_do_radiation )  THEN
                CALL radiation_tendency ( i, j, tend )
             ENDIF
#endif

             CALL module_interface_actions( i, j, 'pt-tendency' )
!
!--          Prognostic equation for potential temperature
             DO  k = nzb+1, nzt
                pt_p(k,j,i) = pt(k,j,i) + ( dt_3d *                            &
                                                  ( tsc(2) * tend(k,j,i) +     &
                                                    tsc(3) * tpt_m(k,j,i) )    &
                                                  - tsc(5)                     &
                                                  * ( pt(k,j,i) - pt_init(k) ) &
                                                  * ( rdf_sc(k) + ptdf_x(i)    &
                                                                + ptdf_y(j) )  &
                                          )                                    &
                                       * MERGE( 1.0_wp, 0.0_wp,                &
                                          BTEST( wall_flags_total_0(k,j,i), 0 )&
                                              )
             ENDDO

!
!--          Calculate tendencies for the next Runge-Kutta step
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( intermediate_timestep_count == 1 )  THEN
                   DO  k = nzb+1, nzt
                      tpt_m(k,j,i) = tend(k,j,i)
                   ENDDO
                ELSEIF ( intermediate_timestep_count < &
                         intermediate_timestep_count_max )  THEN
                   DO  k = nzb+1, nzt
                      tpt_m(k,j,i) =   -9.5625_wp * tend(k,j,i) +              &
                                        5.3125_wp * tpt_m(k,j,i)
                   ENDDO
                ENDIF
             ENDIF

          ENDIF

!
!--       If required, compute prognostic equation for total water content
          IF ( humidity )  THEN

!
!--          Tendency-terms for total water content / scalar
             tend(:,j,i) = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' )                            &
             THEN
                IF ( ws_scheme_sca )  THEN
                   CALL advec_s_ws( advc_flags_s,                              &
                                    i, j, q, 'q', flux_s_q,                    &
                                    diss_s_q, flux_l_q, diss_l_q,              &
                                    i_omp_start, tn,                           &
                                    bc_dirichlet_l  .OR.  bc_radiation_l,      &
                                    bc_dirichlet_n  .OR.  bc_radiation_n,      &
                                    bc_dirichlet_r  .OR.  bc_radiation_r,      &
                                    bc_dirichlet_s  .OR.  bc_radiation_s )
                ELSE
                   CALL advec_s_pw( i, j, q )
                ENDIF
             ELSE
                CALL advec_s_up( i, j, q )
             ENDIF
             CALL diffusion_s( i, j, q,                                        &
                               surf_def_h(0)%qsws, surf_def_h(1)%qsws,         &
                               surf_def_h(2)%qsws,                             &
                               surf_lsm_h%qsws,    surf_usm_h%qsws,            &
                               surf_def_v(0)%qsws, surf_def_v(1)%qsws,         &
                               surf_def_v(2)%qsws, surf_def_v(3)%qsws,         &
                               surf_lsm_v(0)%qsws, surf_lsm_v(1)%qsws,         &
                               surf_lsm_v(2)%qsws, surf_lsm_v(3)%qsws,         &
                               surf_usm_v(0)%qsws, surf_usm_v(1)%qsws,         &
                               surf_usm_v(2)%qsws, surf_usm_v(3)%qsws )

!
!--          Sink or source of humidity due to canopy elements
             IF ( plant_canopy )  CALL pcm_tendency( i, j, 5 )

!
!--          Large scale advection
             IF ( large_scale_forcing )  THEN
                CALL ls_advec( i, j, simulated_time, 'q' )
             ENDIF

!
!--          Nudging
             IF ( nudging )  CALL nudge( i, j, simulated_time, 'q' )

!
!--          If required compute influence of large-scale subsidence/ascent
             IF ( large_scale_subsidence  .AND.                                &
                  .NOT. use_subsidence_tendencies )  THEN
                CALL subsidence( i, j, tend, q, q_init, 3 )
             ENDIF

             CALL module_interface_actions( i, j, 'q-tendency' )

!
!--          Prognostic equation for total water content / scalar
             DO  k = nzb+1, nzt
                q_p(k,j,i) = q(k,j,i) + ( dt_3d *                              &
                                                ( tsc(2) * tend(k,j,i) +       &
                                                  tsc(3) * tq_m(k,j,i) )       &
                                                - tsc(5) * rdf_sc(k) *         &
                                                      ( q(k,j,i) - q_init(k) ) &
                                        )                                      &
                                       * MERGE( 1.0_wp, 0.0_wp,                &
                                          BTEST( wall_flags_total_0(k,j,i), 0 )&
                                              )                
                IF ( q_p(k,j,i) < 0.0_wp )  q_p(k,j,i) = 0.1_wp * q(k,j,i)
             ENDDO

!
!--          Calculate tendencies for the next Runge-Kutta step
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( intermediate_timestep_count == 1 )  THEN
                   DO  k = nzb+1, nzt
                      tq_m(k,j,i) = tend(k,j,i)
                   ENDDO
                ELSEIF ( intermediate_timestep_count < &
                         intermediate_timestep_count_max )  THEN
                   DO  k = nzb+1, nzt
                      tq_m(k,j,i) =   -9.5625_wp * tend(k,j,i) +               &
                                       5.3125_wp * tq_m(k,j,i)
                   ENDDO
                ENDIF
             ENDIF

          ENDIF

!
!--       If required, compute prognostic equation for scalar
          IF ( passive_scalar )  THEN
!
!--          Tendency-terms for total water content / scalar
             tend(:,j,i) = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' )                            &
             THEN
                IF ( ws_scheme_sca )  THEN
!
!--                For scalar advection apply monotonic flux limiter near 
!--                topography.
                   CALL advec_s_ws( advc_flags_s,                              &
                                    i, j, s, 's', flux_s_s,                    &
                                    diss_s_s, flux_l_s, diss_l_s, i_omp_start, &
                                    tn,                                        &
                                    bc_dirichlet_l  .OR.  bc_radiation_l,      &
                                    bc_dirichlet_n  .OR.  bc_radiation_n,      &
                                    bc_dirichlet_r  .OR.  bc_radiation_r,      &
                                    bc_dirichlet_s  .OR.  bc_radiation_s,      &
                                    monotonic_limiter_z )
                ELSE
                   CALL advec_s_pw( i, j, s )
                ENDIF
             ELSE
                CALL advec_s_up( i, j, s )
             ENDIF
             CALL diffusion_s( i, j, s,                                        &
                               surf_def_h(0)%ssws, surf_def_h(1)%ssws,         &
                               surf_def_h(2)%ssws,                             &
                               surf_lsm_h%ssws,    surf_usm_h%ssws,            &
                               surf_def_v(0)%ssws, surf_def_v(1)%ssws,         &
                               surf_def_v(2)%ssws, surf_def_v(3)%ssws,         &
                               surf_lsm_v(0)%ssws, surf_lsm_v(1)%ssws,         &
                               surf_lsm_v(2)%ssws, surf_lsm_v(3)%ssws,         &
                               surf_usm_v(0)%ssws, surf_usm_v(1)%ssws,         &
                               surf_usm_v(2)%ssws, surf_usm_v(3)%ssws )

!
!--          Sink or source of scalar concentration due to canopy elements
             IF ( plant_canopy )  CALL pcm_tendency( i, j, 7 )

!
!--          Large scale advection, still need to be extended for scalars
!              IF ( large_scale_forcing )  THEN
!                 CALL ls_advec( i, j, simulated_time, 's' )
!              ENDIF

!
!--          Nudging, still need to be extended for scalars
!              IF ( nudging )  CALL nudge( i, j, simulated_time, 's' )

!
!--          If required compute influence of large-scale subsidence/ascent.
!--          Note, the last argument is of no meaning in this case, as it is
!--          only used in conjunction with large_scale_forcing, which is to
!--          date not implemented for scalars.
             IF ( large_scale_subsidence  .AND.                                &
                  .NOT. use_subsidence_tendencies  .AND.                       &
                  .NOT. large_scale_forcing )  THEN
                CALL subsidence( i, j, tend, s, s_init, 3 )
             ENDIF

             CALL module_interface_actions( i, j, 's-tendency' )

!
!--          Prognostic equation for scalar
             DO  k = nzb+1, nzt
                s_p(k,j,i) = s(k,j,i) + (  dt_3d *                             &
                                            ( tsc(2) * tend(k,j,i) +           &
                                              tsc(3) * ts_m(k,j,i) )           &
                                            - tsc(5) * rdf_sc(k)               &
                                                     * ( s(k,j,i) - s_init(k) )&
                                        )                                      &
                                       * MERGE( 1.0_wp, 0.0_wp,                &
                                          BTEST( wall_flags_total_0(k,j,i), 0 )&
                                              )
                IF ( s_p(k,j,i) < 0.0_wp )  s_p(k,j,i) = 0.1_wp * s(k,j,i)
             ENDDO

!
!--          Calculate tendencies for the next Runge-Kutta step
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( intermediate_timestep_count == 1 )  THEN
                   DO  k = nzb+1, nzt
                      ts_m(k,j,i) = tend(k,j,i)
                   ENDDO
                ELSEIF ( intermediate_timestep_count < &
                         intermediate_timestep_count_max )  THEN
                   DO  k = nzb+1, nzt
                      ts_m(k,j,i) =   -9.5625_wp * tend(k,j,i) +               &
                                       5.3125_wp * ts_m(k,j,i)
                   ENDDO
                ENDIF
             ENDIF

          ENDIF
!
!--       Calculate prognostic equations for all other modules
          CALL module_interface_prognostic_equations( i, j, i_omp_start, tn )

       ENDDO  ! loop over j
    ENDDO  ! loop over i
    !$OMP END PARALLEL


    CALL cpu_log( log_point(32), 'all progn.equations', 'stop' )

    IF ( debug_output_timestep )  CALL debug_message( 'prognostic_equations_cache', 'end' )

 END SUBROUTINE prognostic_equations_cache


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Version for vector machines
!------------------------------------------------------------------------------!

 SUBROUTINE prognostic_equations_vector


    INTEGER(iwp) ::  i     !<
    INTEGER(iwp) ::  j     !<
    INTEGER(iwp) ::  k     !<

    REAL(wp)     ::  sbt  !<


    IF ( debug_output_timestep )  CALL debug_message( 'prognostic_equations_vector', 'start' )
!
!-- Calculate non advective processes for all other modules
    CALL module_interface_non_advective_processes
!
!-- Module Inferface for exchange horiz after non_advective_processes but before
!-- advection. Therefore, non_advective_processes must not run for ghost points.
    CALL module_interface_exchange_horiz()
!
!-- u-velocity component
    CALL cpu_log( log_point(5), 'u-equation', 'start' )

    !$ACC KERNELS PRESENT(tend)
    tend = 0.0_wp
    !$ACC END KERNELS
    IF ( timestep_scheme(1:5) == 'runge' )  THEN
       IF ( ws_scheme_mom )  THEN
          CALL advec_u_ws
       ELSE
          CALL advec_u_pw
       ENDIF
    ELSE
       CALL advec_u_up
    ENDIF
    CALL diffusion_u
    CALL coriolis( 1 )
    IF ( sloping_surface  .AND.  .NOT. neutral )  THEN
       CALL buoyancy( pt, 1 )
    ENDIF

!
!-- Drag by plant canopy
    IF ( plant_canopy )  CALL pcm_tendency( 1 )

!
!-- External pressure gradient
    IF ( dp_external )  THEN
       DO  i = nxlu, nxr
          DO  j = nys, nyn
             DO  k = dp_level_ind_b+1, nzt
                tend(k,j,i) = tend(k,j,i) - dpdxy(1) * dp_smooth_factor(k)
             ENDDO
          ENDDO
       ENDDO
    ENDIF

!
!-- Nudging
    IF ( nudging )  CALL nudge( simulated_time, 'u' )

!
!-- Effect of Stokes drift (in ocean mode only)
    IF ( stokes_force )  CALL stokes_drift_terms( 1 )

    CALL module_interface_actions( 'u-tendency' )

!
!-- Prognostic equation for u-velocity component
    !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i, j, k) &
    !$ACC PRESENT(u, tend, tu_m, u_init, rdf, wall_flags_total_0) &
    !$ACC PRESENT(tsc(2:5)) &
    !$ACC PRESENT(u_p)
    DO  i = nxlu, nxr
       DO  j = nys, nyn
          !following directive is required to vectorize on Intel19
          !DIR$ IVDEP
          DO  k = nzb+1, nzt
             u_p(k,j,i) = u(k,j,i) + ( dt_3d * ( tsc(2) * tend(k,j,i) +          &
                                                 tsc(3) * tu_m(k,j,i) )          &
                                               - tsc(5) * rdf(k) *               &
                                                        ( u(k,j,i) - u_init(k) ) &
                                     ) * MERGE( 1.0_wp, 0.0_wp,                  &
                                           BTEST( wall_flags_total_0(k,j,i), 1 ) &
                                              )
          ENDDO
       ENDDO
    ENDDO

!
!-- Add turbulence generated by wave breaking (in ocean mode only)
    IF ( wave_breaking  .AND.                                                  &
         intermediate_timestep_count == intermediate_timestep_count_max )      &
    THEN
       CALL wave_breaking_term( 1 )
    ENDIF

!
!-- Calculate tendencies for the next Runge-Kutta step
    IF ( timestep_scheme(1:5) == 'runge' )  THEN
       IF ( intermediate_timestep_count == 1 )  THEN
          !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i, j, k) &
          !$ACC PRESENT(tend, tu_m)
          DO  i = nxlu, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   tu_m(k,j,i) = tend(k,j,i)
                ENDDO
             ENDDO
          ENDDO
       ELSEIF ( intermediate_timestep_count < &
                intermediate_timestep_count_max )  THEN
          !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i, j, k) &
          !$ACC PRESENT(tend, tu_m)
          DO  i = nxlu, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   tu_m(k,j,i) =    -9.5625_wp * tend(k,j,i)                   &
                                   + 5.3125_wp * tu_m(k,j,i)
                ENDDO
             ENDDO
          ENDDO
       ENDIF
    ENDIF

    CALL cpu_log( log_point(5), 'u-equation', 'stop' )

!
!-- v-velocity component
    CALL cpu_log( log_point(6), 'v-equation', 'start' )

    !$ACC KERNELS PRESENT(tend)
    tend = 0.0_wp
    !$ACC END KERNELS
    IF ( timestep_scheme(1:5) == 'runge' )  THEN
       IF ( ws_scheme_mom )  THEN
          CALL advec_v_ws
       ELSE
          CALL advec_v_pw
       END IF
    ELSE
       CALL advec_v_up
    ENDIF
    CALL diffusion_v
    CALL coriolis( 2 )

!
!-- Drag by plant canopy
    IF ( plant_canopy )  CALL pcm_tendency( 2 )

!
!-- External pressure gradient
    IF ( dp_external )  THEN
       DO  i = nxl, nxr
          DO  j = nysv, nyn
             DO  k = dp_level_ind_b+1, nzt
                tend(k,j,i) = tend(k,j,i) - dpdxy(2) * dp_smooth_factor(k)
             ENDDO
          ENDDO
       ENDDO
    ENDIF

!
!-- Nudging
    IF ( nudging )  CALL nudge( simulated_time, 'v' )

!
!-- Effect of Stokes drift (in ocean mode only)
    IF ( stokes_force )  CALL stokes_drift_terms( 2 )

    CALL module_interface_actions( 'v-tendency' )

!
!-- Prognostic equation for v-velocity component
    !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i, j, k) &
    !$ACC PRESENT(v, tend, tv_m, v_init, rdf, wall_flags_total_0) &
    !$ACC PRESENT(tsc(2:5)) &
    !$ACC PRESENT(v_p)
    DO  i = nxl, nxr
       DO  j = nysv, nyn
          !following directive is required to vectorize on Intel19
          !DIR$ IVDEP
          DO  k = nzb+1, nzt
             v_p(k,j,i) = v(k,j,i) + ( dt_3d * ( tsc(2) * tend(k,j,i) +        &
                                                 tsc(3) * tv_m(k,j,i) )        &
                                               - tsc(5) * rdf(k) *             &
                                                      ( v(k,j,i) - v_init(k) ) &
                                     ) * MERGE( 1.0_wp, 0.0_wp,                &
                                          BTEST( wall_flags_total_0(k,j,i), 2 )&
                                              )
          ENDDO
       ENDDO
    ENDDO

!
!-- Add turbulence generated by wave breaking (in ocean mode only)
    IF ( wave_breaking  .AND.                                                  &
         intermediate_timestep_count == intermediate_timestep_count_max )      &
    THEN
       CALL wave_breaking_term( 2 )
    ENDIF

!
!-- Calculate tendencies for the next Runge-Kutta step
    IF ( timestep_scheme(1:5) == 'runge' )  THEN
       IF ( intermediate_timestep_count == 1 )  THEN
          !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i, j, k) &
          !$ACC PRESENT(tend, tv_m)
          DO  i = nxl, nxr
             DO  j = nysv, nyn
                DO  k = nzb+1, nzt
                   tv_m(k,j,i) = tend(k,j,i)
                ENDDO
             ENDDO
          ENDDO
       ELSEIF ( intermediate_timestep_count < &
                intermediate_timestep_count_max )  THEN
          !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i, j, k) &
          !$ACC PRESENT(tend, tv_m)
          DO  i = nxl, nxr
             DO  j = nysv, nyn
                DO  k = nzb+1, nzt
                   tv_m(k,j,i) =   -9.5625_wp * tend(k,j,i)                    &
                                  + 5.3125_wp * tv_m(k,j,i)
                ENDDO
             ENDDO
          ENDDO
       ENDIF
    ENDIF

    CALL cpu_log( log_point(6), 'v-equation', 'stop' )

!
!-- w-velocity component
    CALL cpu_log( log_point(7), 'w-equation', 'start' )

    !$ACC KERNELS PRESENT(tend)
    tend = 0.0_wp
    !$ACC END KERNELS
    IF ( timestep_scheme(1:5) == 'runge' )  THEN
       IF ( ws_scheme_mom )  THEN
          CALL advec_w_ws
       ELSE
          CALL advec_w_pw
       ENDIF
    ELSE
       CALL advec_w_up
    ENDIF
    CALL diffusion_w
    CALL coriolis( 3 )

    IF ( .NOT. neutral )  THEN
       IF ( ocean_mode )  THEN
          CALL buoyancy( rho_ocean, 3 )
       ELSE
          IF ( .NOT. humidity )  THEN
             CALL buoyancy( pt, 3 )
          ELSE
             CALL buoyancy( vpt, 3 )
          ENDIF
       ENDIF
    ENDIF

!
!-- Drag by plant canopy
    IF ( plant_canopy )  CALL pcm_tendency( 3 )

!
!-- Effect of Stokes drift (in ocean mode only)
    IF ( stokes_force )  CALL stokes_drift_terms( 3 )

    CALL module_interface_actions( 'w-tendency' )

!
!-- Prognostic equation for w-velocity component
    !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i, j, k) &
    !$ACC PRESENT(w, tend, tw_m, v_init, rdf, wall_flags_total_0) &
    !$ACC PRESENT(tsc(2:5)) &
    !$ACC PRESENT(w_p)
    DO  i = nxl, nxr
       DO  j = nys, nyn
          !following directive is required to vectorize on Intel19
          !DIR$ IVDEP
          DO  k = nzb+1, nzt-1
             w_p(k,j,i) = w(k,j,i) + ( dt_3d * ( tsc(2) * tend(k,j,i) +        &
                                                 tsc(3) * tw_m(k,j,i) )        &
                                               - tsc(5) * rdf(k) * w(k,j,i)    &
                                     ) * MERGE( 1.0_wp, 0.0_wp,                &
                                          BTEST( wall_flags_total_0(k,j,i), 3 )&
                                              )
          ENDDO
       ENDDO
    ENDDO

!
!-- Calculate tendencies for the next Runge-Kutta step
    IF ( timestep_scheme(1:5) == 'runge' )  THEN
       IF ( intermediate_timestep_count == 1 )  THEN
          !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i, j, k) &
          !$ACC PRESENT(tend, tw_m)
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt-1
                   tw_m(k,j,i) = tend(k,j,i)
                ENDDO
             ENDDO
          ENDDO
       ELSEIF ( intermediate_timestep_count < &
                intermediate_timestep_count_max )  THEN
          !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i, j, k) &
          !$ACC PRESENT(tend, tw_m)
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt-1
                   tw_m(k,j,i) =   -9.5625_wp * tend(k,j,i)                    &
                                  + 5.3125_wp * tw_m(k,j,i)
                ENDDO
             ENDDO
          ENDDO
       ENDIF
    ENDIF

    CALL cpu_log( log_point(7), 'w-equation', 'stop' )


!
!-- If required, compute prognostic equation for potential temperature
    IF ( .NOT. neutral )  THEN

       CALL cpu_log( log_point(13), 'pt-equation', 'start' )

!
!--    pt-tendency terms with communication
       sbt = tsc(2)
       IF ( scalar_advec == 'bc-scheme' )  THEN

          IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--          Bott-Chlond scheme always uses Euler time step. Thus:
             sbt = 1.0_wp
          ENDIF
          tend = 0.0_wp
          CALL advec_s_bc( pt, 'pt' )

       ENDIF

!
!--    pt-tendency terms with no communication
       IF ( scalar_advec /= 'bc-scheme' )  THEN
          !$ACC KERNELS PRESENT(tend)
          tend = 0.0_wp
          !$ACC END KERNELS
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( ws_scheme_sca )  THEN
                CALL advec_s_ws( advc_flags_s, pt, 'pt',                       &
                                 bc_dirichlet_l  .OR.  bc_radiation_l,         &
                                 bc_dirichlet_n  .OR.  bc_radiation_n,         &
                                 bc_dirichlet_r  .OR.  bc_radiation_r,         &
                                 bc_dirichlet_s  .OR.  bc_radiation_s )
             ELSE
                CALL advec_s_pw( pt )
             ENDIF
          ELSE
             CALL advec_s_up( pt )
          ENDIF
       ENDIF

       CALL diffusion_s( pt,                                                   &
                         surf_def_h(0)%shf, surf_def_h(1)%shf,                 &
                         surf_def_h(2)%shf,                                    &
                         surf_lsm_h%shf,    surf_usm_h%shf,                    &
                         surf_def_v(0)%shf, surf_def_v(1)%shf,                 &
                         surf_def_v(2)%shf, surf_def_v(3)%shf,                 &
                         surf_lsm_v(0)%shf, surf_lsm_v(1)%shf,                 &
                         surf_lsm_v(2)%shf, surf_lsm_v(3)%shf,                 &
                         surf_usm_v(0)%shf, surf_usm_v(1)%shf,                 &
                         surf_usm_v(2)%shf, surf_usm_v(3)%shf )

!
!--    Consideration of heat sources within the plant canopy
       IF ( plant_canopy  .AND.                                          &
            (cthf /= 0.0_wp  .OR. urban_surface  .OR.  land_surface) )  THEN
          CALL pcm_tendency( 4 )
       ENDIF

!
!--    Large scale advection
       IF ( large_scale_forcing )  THEN
          CALL ls_advec( simulated_time, 'pt' )
       ENDIF

!
!--    Nudging
       IF ( nudging )  CALL nudge( simulated_time, 'pt' )

!
!--    If required compute influence of large-scale subsidence/ascent
       IF ( large_scale_subsidence  .AND.                                      &
            .NOT. use_subsidence_tendencies )  THEN
          CALL subsidence( tend, pt, pt_init, 2 )
       ENDIF

#if defined( __rrtmg )
!
!--    If required, add tendency due to radiative heating/cooling
       IF ( radiation  .AND.                                                   &
            simulated_time > skip_time_do_radiation )  THEN
            CALL radiation_tendency ( tend )
       ENDIF
#endif

       CALL module_interface_actions( 'pt-tendency' )

!
!--    Prognostic equation for potential temperature
       !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i, j, k) &
       !$ACC PRESENT(pt, tend, tpt_m, wall_flags_total_0) &
       !$ACC PRESENT(pt_init, rdf_sc, ptdf_x, ptdf_y) &
       !$ACC PRESENT(tsc(3:5)) &
       !$ACC PRESENT(pt_p)
       DO  i = nxl, nxr
          DO  j = nys, nyn
             !following directive is required to vectorize on Intel19
             !DIR$ IVDEP
             DO  k = nzb+1, nzt
                pt_p(k,j,i) = pt(k,j,i) + ( dt_3d * ( sbt * tend(k,j,i) +      &
                                                   tsc(3) * tpt_m(k,j,i) )     &
                                                 - tsc(5) *                    &
                                                   ( pt(k,j,i) - pt_init(k) ) *&
                                          ( rdf_sc(k) + ptdf_x(i) + ptdf_y(j) )&
                                          )                                    &
                                   * MERGE( 1.0_wp, 0.0_wp,                    &
                                       BTEST( wall_flags_total_0(k,j,i), 0 )   &
                                          )
             ENDDO
          ENDDO
       ENDDO
!
!--    Calculate tendencies for the next Runge-Kutta step
       IF ( timestep_scheme(1:5) == 'runge' )  THEN
          IF ( intermediate_timestep_count == 1 )  THEN
             !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i, j, k) &
             !$ACC PRESENT(tend, tpt_m)
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt
                      tpt_m(k,j,i) = tend(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSEIF ( intermediate_timestep_count < &
                   intermediate_timestep_count_max )  THEN
             !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i, j, k) &
             !$ACC PRESENT(tend, tpt_m)
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt
                      tpt_m(k,j,i) =   -9.5625_wp * tend(k,j,i) +              &
                                        5.3125_wp * tpt_m(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDIF

       CALL cpu_log( log_point(13), 'pt-equation', 'stop' )

    ENDIF

!
!-- If required, compute prognostic equation for total water content
    IF ( humidity )  THEN

       CALL cpu_log( log_point(29), 'q-equation', 'start' )

!
!--    Scalar/q-tendency terms with communication
       sbt = tsc(2)
       IF ( scalar_advec == 'bc-scheme' )  THEN

          IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--          Bott-Chlond scheme always uses Euler time step. Thus:
             sbt = 1.0_wp
          ENDIF
          tend = 0.0_wp
          CALL advec_s_bc( q, 'q' )

       ENDIF

!
!--    Scalar/q-tendency terms with no communication
       IF ( scalar_advec /= 'bc-scheme' )  THEN
          tend = 0.0_wp
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( ws_scheme_sca )  THEN
                CALL advec_s_ws( advc_flags_s, q, 'q',                         &
                                 bc_dirichlet_l  .OR.  bc_radiation_l,         &
                                 bc_dirichlet_n  .OR.  bc_radiation_n,         &
                                 bc_dirichlet_r  .OR.  bc_radiation_r,         &
                                 bc_dirichlet_s  .OR.  bc_radiation_s )
             ELSE
                CALL advec_s_pw( q )
             ENDIF
          ELSE
             CALL advec_s_up( q )
          ENDIF
       ENDIF

       CALL diffusion_s( q,                                                    &
                         surf_def_h(0)%qsws, surf_def_h(1)%qsws,               &
                         surf_def_h(2)%qsws,                                   &
                         surf_lsm_h%qsws,    surf_usm_h%qsws,                  &
                         surf_def_v(0)%qsws, surf_def_v(1)%qsws,               &
                         surf_def_v(2)%qsws, surf_def_v(3)%qsws,               &
                         surf_lsm_v(0)%qsws, surf_lsm_v(1)%qsws,               &
                         surf_lsm_v(2)%qsws, surf_lsm_v(3)%qsws,               &
                         surf_usm_v(0)%qsws, surf_usm_v(1)%qsws,               &
                         surf_usm_v(2)%qsws, surf_usm_v(3)%qsws )

!
!--    Sink or source of humidity due to canopy elements
       IF ( plant_canopy ) CALL pcm_tendency( 5 )

!
!--    Large scale advection
       IF ( large_scale_forcing )  THEN
          CALL ls_advec( simulated_time, 'q' )
       ENDIF

!
!--    Nudging
       IF ( nudging )  CALL nudge( simulated_time, 'q' )

!
!--    If required compute influence of large-scale subsidence/ascent
       IF ( large_scale_subsidence  .AND.                                      &
            .NOT. use_subsidence_tendencies )  THEN
         CALL subsidence( tend, q, q_init, 3 )
       ENDIF

       CALL module_interface_actions( 'q-tendency' )

!
!--    Prognostic equation for total water content
       DO  i = nxl, nxr
          DO  j = nys, nyn
             !following directive is required to vectorize on Intel19
             !DIR$ IVDEP
             DO  k = nzb+1, nzt
                q_p(k,j,i) = q(k,j,i) + ( dt_3d * ( sbt * tend(k,j,i) +        &
                                                 tsc(3) * tq_m(k,j,i) )        &
                                               - tsc(5) * rdf_sc(k) *          &
                                                      ( q(k,j,i) - q_init(k) ) &
                                        ) * MERGE( 1.0_wp, 0.0_wp,             &
                                         BTEST( wall_flags_total_0(k,j,i), 0 ) &
                                                 )
                IF ( q_p(k,j,i) < 0.0_wp )  q_p(k,j,i) = 0.1_wp * q(k,j,i)
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
                      tq_m(k,j,i) = tend(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSEIF ( intermediate_timestep_count < &
                   intermediate_timestep_count_max )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt
                      tq_m(k,j,i) =   -9.5625_wp * tend(k,j,i)                 &
                                     + 5.3125_wp * tq_m(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDIF

       CALL cpu_log( log_point(29), 'q-equation', 'stop' )

    ENDIF
!
!-- If required, compute prognostic equation for scalar
    IF ( passive_scalar )  THEN

       CALL cpu_log( log_point(66), 's-equation', 'start' )

!
!--    Scalar/q-tendency terms with communication
       sbt = tsc(2)
       IF ( scalar_advec == 'bc-scheme' )  THEN

          IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--          Bott-Chlond scheme always uses Euler time step. Thus:
             sbt = 1.0_wp
          ENDIF
          tend = 0.0_wp
          CALL advec_s_bc( s, 's' )

       ENDIF

!
!--    Scalar/q-tendency terms with no communication
       IF ( scalar_advec /= 'bc-scheme' )  THEN
          tend = 0.0_wp
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( ws_scheme_sca )  THEN
                CALL advec_s_ws( advc_flags_s, s, 's',                         &
                                 bc_dirichlet_l  .OR.  bc_radiation_l,         &
                                 bc_dirichlet_n  .OR.  bc_radiation_n,         &
                                 bc_dirichlet_r  .OR.  bc_radiation_r,         &
                                 bc_dirichlet_s  .OR.  bc_radiation_s )
             ELSE
                CALL advec_s_pw( s )
             ENDIF
          ELSE
             CALL advec_s_up( s )
          ENDIF
       ENDIF

       CALL diffusion_s( s,                                                    &
                         surf_def_h(0)%ssws, surf_def_h(1)%ssws,               &
                         surf_def_h(2)%ssws,                                   &
                         surf_lsm_h%ssws,    surf_usm_h%ssws,                  &
                         surf_def_v(0)%ssws, surf_def_v(1)%ssws,               &
                         surf_def_v(2)%ssws, surf_def_v(3)%ssws,               &
                         surf_lsm_v(0)%ssws, surf_lsm_v(1)%ssws,               &
                         surf_lsm_v(2)%ssws, surf_lsm_v(3)%ssws,               &
                         surf_usm_v(0)%ssws, surf_usm_v(1)%ssws,               &
                         surf_usm_v(2)%ssws, surf_usm_v(3)%ssws )

!
!--    Sink or source of humidity due to canopy elements
       IF ( plant_canopy ) CALL pcm_tendency( 7 )

!
!--    Large scale advection. Not implemented for scalars so far.
!        IF ( large_scale_forcing )  THEN
!           CALL ls_advec( simulated_time, 'q' )
!        ENDIF

!
!--    Nudging. Not implemented for scalars so far.
!        IF ( nudging )  CALL nudge( simulated_time, 'q' )

!
!--    If required compute influence of large-scale subsidence/ascent.
!--    Not implemented for scalars so far.
       IF ( large_scale_subsidence  .AND.                                      &
            .NOT. use_subsidence_tendencies  .AND.                             &
            .NOT. large_scale_forcing )  THEN
         CALL subsidence( tend, s, s_init, 3 )
       ENDIF

       CALL module_interface_actions( 's-tendency' )

!
!--    Prognostic equation for total water content
       DO  i = nxl, nxr
          DO  j = nys, nyn
             !following directive is required to vectorize on Intel19
             !DIR$ IVDEP
             DO  k = nzb+1, nzt
                s_p(k,j,i) = s(k,j,i) + ( dt_3d * ( sbt * tend(k,j,i) +        &
                                                 tsc(3) * ts_m(k,j,i) )        &
                                               - tsc(5) * rdf_sc(k) *          &
                                                    ( s(k,j,i) - s_init(k) )   &
                                        )                                      &
                                   * MERGE( 1.0_wp, 0.0_wp,                    &
                                       BTEST( wall_flags_total_0(k,j,i), 0 )   &
                                          )
                IF ( s_p(k,j,i) < 0.0_wp )  s_p(k,j,i) = 0.1_wp * s(k,j,i)
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
                      ts_m(k,j,i) = tend(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSEIF ( intermediate_timestep_count < &
                   intermediate_timestep_count_max )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt
                      ts_m(k,j,i) =   -9.5625_wp * tend(k,j,i)                 &
                                     + 5.3125_wp * ts_m(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDIF

       CALL cpu_log( log_point(66), 's-equation', 'stop' )

    ENDIF
!
!-- Calculate prognostic equations for all other modules
    CALL module_interface_prognostic_equations()

    IF ( debug_output_timestep )  CALL debug_message( 'prognostic_equations_vector', 'end' )

 END SUBROUTINE prognostic_equations_vector


 END MODULE prognostic_equations_mod
