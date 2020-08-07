!> @file flow_statistics.f90
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
! $Id: flow_statistics.f90 4464 2020-03-17 11:08:46Z Giersch $
! Reset last change (r4463)
! 
! 4463 2020-03-17 09:27:36Z Giersch
! Calculate horizontally averaged profiles of all velocity components at the
! same place
!
! 4444 2020-03-05 15:59:50Z raasch
! bugfix: cpp-directives for serial mode added
! 
! 4442 2020-03-04 19:21:13Z suehring
! Change order of dimension in surface array %frac to allow for better 
! vectorization.
! 
! 4441 2020-03-04 19:20:35Z suehring
! Introduction of wall_flags_total_0, which currently sets bits based on static
! topography information used in wall_flags_static_0
! 
! 4329 2019-12-10 15:46:36Z motisi
! Renamed wall_flags_0 to wall_flags_static_0
! 
! 4182 2019-08-22 15:20:23Z scharf
! Corrected "Former revisions" section
! 
! 4131 2019-08-02 11:06:18Z monakurppa
! Allow profile output for salsa variables.
! 
! 4039 2019-06-18 10:32:41Z suehring
! Correct conversion to kinematic scalar fluxes in case of pw-scheme and 
! statistic regions
! 
! 3828 2019-03-27 19:36:23Z raasch
! unused variables removed
! 
! 3676 2019-01-16 15:07:05Z knoop
! Bugfix, terminate OMP Parallel block
!
! Revision 1.1  1997/08/11 06:15:17  raasch
! Initial revision
!
!
! Description:
! ------------
!> Compute average profiles and further average flow quantities for the different
!> user-defined (sub-)regions. The region indexed 0 is the total model domain.
!>
!> @note For simplicity, nzb_s_inner and nzb_diff_s_inner are being used as a
!>       lower vertical index for k-loops for all variables, although strictly 
!>       speaking the k-loops would have to be split up according to the staggered 
!>       grid. However, this implies no error since staggered velocity components 
!>       are zero at the walls and inside buildings.
!------------------------------------------------------------------------------!
 SUBROUTINE flow_statistics


    USE arrays_3d,                                                             &
        ONLY:  ddzu, ddzw, e, heatflux_output_conversion, hyp, km, kh,         &
               momentumflux_output_conversion, nc, nr, p, prho, prr, pt, q,    &
               qc, ql, qr, rho_air, rho_air_zw, rho_ocean, s,                  &
               sa, u, ug, v, vg, vpt, w, w_subs, waterflux_output_conversion,  &
               zw, d_exner

    USE basic_constants_and_equations_mod,                                     &
        ONLY:  g, lv_d_cp

    USE bulk_cloud_model_mod,                                                  &
        ONLY: bulk_cloud_model, microphysics_morrison, microphysics_seifert

    USE chem_modules,                                                          &
        ONLY:  max_pr_cs

    USE control_parameters,                                                    &
        ONLY:   air_chemistry, average_count_pr, cloud_droplets, do_sum,       &
                dt_3d, humidity, initializing_actions, land_surface,           &
                large_scale_forcing, large_scale_subsidence, max_pr_user,      &
                message_string, neutral, ocean_mode, passive_scalar,           &
                simulated_time, simulated_time_at_begin,                       &
                use_subsidence_tendencies, use_surface_fluxes, use_top_fluxes, &
                ws_scheme_mom, ws_scheme_sca, salsa, max_pr_salsa

    USE cpulog,                                                                &
        ONLY:   cpu_log, log_point

    USE grid_variables,                                                        &
        ONLY:   ddx, ddy
        
    USE indices,                                                               &
        ONLY:   ngp_2dh, ngp_2dh_s_inner, ngp_3d, ngp_3d_inner, nxl, nxr, nyn, &
                nys, nzb, nzt, topo_min_level, wall_flags_total_0

#if defined( __parallel )
    USE indices,                                                               &
        ONLY:  ngp_sums, ngp_sums_ls
#endif
        
    USE kinds
    
    USE land_surface_model_mod,                                                &
        ONLY:   m_soil_h, nzb_soil, nzt_soil, t_soil_h

    USE lsf_nudging_mod,                                                       &
        ONLY:   td_lsa_lpt, td_lsa_q, td_sub_lpt, td_sub_q, time_vert

    USE module_interface,                                                      &
        ONLY:  module_interface_statistics

    USE netcdf_interface,                                                      &
        ONLY:  dots_rad, dots_soil, dots_max

    USE pegrid

    USE radiation_model_mod,                                                   &
        ONLY:  radiation, radiation_scheme,                                    &
               rad_lw_in, rad_lw_out, rad_lw_cs_hr, rad_lw_hr,                 &
               rad_sw_in, rad_sw_out, rad_sw_cs_hr, rad_sw_hr

    USE statistics

    USE surface_mod,                                                           &
        ONLY :  surf_def_h, surf_lsm_h, surf_usm_h


    IMPLICIT NONE

    INTEGER(iwp) ::  i                   !<
    INTEGER(iwp) ::  j                   !<
    INTEGER(iwp) ::  k                   !<
    INTEGER(iwp) ::  ki                  !< 
    INTEGER(iwp) ::  k_surface_level     !<
    INTEGER(iwp) ::  m                   !< loop variable over all horizontal wall elements 
    INTEGER(iwp) ::  l                   !< loop variable over surface facing -- up- or downward-facing
    INTEGER(iwp) ::  nt                  !<
!$  INTEGER(iwp) ::  omp_get_thread_num  !<
    INTEGER(iwp) ::  sr                  !<
    INTEGER(iwp) ::  tn                  !<

    LOGICAL ::  first  !<
    
    REAL(wp) ::  dptdz_threshold  !< 
    REAL(wp) ::  fac              !<
    REAL(wp) ::  flag             !<
    REAL(wp) ::  height           !<
    REAL(wp) ::  pts              !<
    REAL(wp) ::  sums_l_etot      !<
    REAL(wp) ::  ust              !<
    REAL(wp) ::  ust2             !<
    REAL(wp) ::  u2               !<
    REAL(wp) ::  vst              !<
    REAL(wp) ::  vst2             !<
    REAL(wp) ::  v2               !<
    REAL(wp) ::  w2               !<
    
    REAL(wp) ::  dptdz(nzb+1:nzt+1)    !<
    REAL(wp) ::  sums_ll(nzb:nzt+1,2)  !<

    CALL cpu_log( log_point(10), 'flow_statistics', 'start' )


!
!-- To be on the safe side, check whether flow_statistics has already been
!-- called once after the current time step
    IF ( flow_statistics_called )  THEN

       message_string = 'flow_statistics is called two times within one ' // &
                        'timestep'
       CALL message( 'flow_statistics', 'PA0190', 1, 2, 0, 6, 0 )

    ENDIF

!
!-- Compute statistics for each (sub-)region
    DO  sr = 0, statistic_regions

!
!--    Initialize (local) summation array
       sums_l = 0.0_wp
#ifdef _OPENACC
       !$ACC KERNELS PRESENT(sums_l)
       sums_l = 0.0_wp
       !$ACC END KERNELS
#endif

!
!--    Store sums that have been computed in other subroutines in summation
!--    array
       sums_l(:,11,:) = sums_l_l(:,sr,:)      ! mixing length from diffusivities
!--    WARNING: next line still has to be adjusted for OpenMP 
       sums_l(:,21,0) = sums_wsts_bc_l(:,sr) *                                 &
                        heatflux_output_conversion  ! heat flux from advec_s_bc
       sums_l(nzb+9,pr_palm,0)  = sums_divold_l(sr)  ! old divergence from pres
       sums_l(nzb+10,pr_palm,0) = sums_divnew_l(sr)  ! new divergence from pres

!
!--    When calcuating horizontally-averaged total (resolved- plus subgrid-
!--    scale) vertical fluxes and velocity variances by using commonly-
!--    applied Reynolds-based methods ( e.g. <w'pt'> = (w-<w>)*(pt-<pt>) )
!--    in combination with the 5th order advection scheme, pronounced 
!--    artificial kinks could be observed in the vertical profiles near the 
!--    surface. Please note: these kinks were not related to the model truth, 
!--    i.e. these kinks are just related to an evaluation problem.   
!--    In order avoid these kinks, vertical fluxes and horizontal as well 
!--    vertical velocity variances are calculated directly within the advection
!--    routines, according to the numerical discretization, to evaluate the 
!--    statistical quantities as they will appear within the prognostic 
!--    equations.
!--    Copy the turbulent quantities, evaluated in the advection routines to 
!--    the local array sums_l() for further computations.
       IF ( ws_scheme_mom .AND. sr == 0 )  THEN

!
!--       According to the Neumann bc for the horizontal velocity components,
!--       the corresponding fluxes has to satisfiy the same bc.
          IF ( ocean_mode )  THEN
             sums_us2_ws_l(nzt+1,:) = sums_us2_ws_l(nzt,:)
             sums_vs2_ws_l(nzt+1,:) = sums_vs2_ws_l(nzt,:)
          ENDIF

          DO  i = 0, threads_per_task-1
!
!--          Swap the turbulent quantities evaluated in advec_ws.
             sums_l(:,13,i) = sums_wsus_ws_l(:,i)                              &
                              * momentumflux_output_conversion ! w*u*
             sums_l(:,15,i) = sums_wsvs_ws_l(:,i)                              &
                              * momentumflux_output_conversion ! w*v*
             sums_l(:,30,i) = sums_us2_ws_l(:,i)        ! u*2 
             sums_l(:,31,i) = sums_vs2_ws_l(:,i)        ! v*2 
             sums_l(:,32,i) = sums_ws2_ws_l(:,i)        ! w*2 
             sums_l(:,34,i) = sums_l(:,34,i) + 0.5_wp *                        & 
                              ( sums_us2_ws_l(:,i) + sums_vs2_ws_l(:,i) +      &
                                sums_ws2_ws_l(:,i) )    ! e*
          ENDDO

       ENDIF

       IF ( ws_scheme_sca .AND. sr == 0 )  THEN

          DO  i = 0, threads_per_task-1
             sums_l(:,17,i)                        = sums_wspts_ws_l(:,i)      &
                                           * heatflux_output_conversion  ! w*pt*
             IF ( ocean_mode     ) sums_l(:,66,i)  = sums_wssas_ws_l(:,i) ! w*sa*
             IF ( humidity       ) sums_l(:,49,i)  = sums_wsqs_ws_l(:,i)       &
                                           * waterflux_output_conversion  ! w*q*
             IF ( passive_scalar ) sums_l(:,114,i) = sums_wsss_ws_l(:,i)  ! w*s*
          ENDDO

       ENDIF
! 
!--    Horizontally averaged profiles of horizontal velocities and temperature.
!--    They must have been computed before, because they are already required
!--    for other horizontal averages.
       tn = 0
       !$OMP PARALLEL PRIVATE( i, j, k, tn, flag )
       !$ tn = omp_get_thread_num()
       !$OMP DO
       !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i, j, k, flag) &
       !$ACC PRESENT(wall_flags_total_0, u, v, pt, rmask, sums_l)
       DO  i = nxl, nxr
          DO  j =  nys, nyn
             DO  k = nzb, nzt+1
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 22 ) )
                !$ACC ATOMIC
                sums_l(k,1,tn)  = sums_l(k,1,tn)  + u(k,j,i)  * rmask(j,i,sr)  &
                                                              * flag
                !$ACC ATOMIC
                sums_l(k,2,tn)  = sums_l(k,2,tn)  + v(k,j,i)  * rmask(j,i,sr)  &
                                                              * flag
                !$ACC ATOMIC
                sums_l(k,4,tn)  = sums_l(k,4,tn)  + pt(k,j,i) * rmask(j,i,sr)  &
                                                              * flag
             ENDDO
          ENDDO
       ENDDO
       !$ACC UPDATE HOST(sums_l(:,1,tn), sums_l(:,2,tn), sums_l(:,4,tn))

!
!--    Horizontally averaged profile of salinity
       IF ( ocean_mode )  THEN
          !$OMP DO
          DO  i = nxl, nxr
             DO  j =  nys, nyn
                DO  k = nzb, nzt+1
                   sums_l(k,23,tn)  = sums_l(k,23,tn) + sa(k,j,i)              &
                                    * rmask(j,i,sr)                            &
                                    * MERGE( 1.0_wp, 0.0_wp,                   &
                                             BTEST( wall_flags_total_0(k,j,i), 22 ) )
                ENDDO
             ENDDO
          ENDDO
       ENDIF

!
!--    Horizontally averaged profiles of virtual potential temperature,
!--    total water content, water vapor mixing ratio and liquid water potential
!--    temperature
       IF ( humidity )  THEN
          !$OMP DO
          DO  i = nxl, nxr
             DO  j =  nys, nyn
                DO  k = nzb, nzt+1
                   flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 22 ) )
                   sums_l(k,44,tn)  = sums_l(k,44,tn) +                        &
                                      vpt(k,j,i) * rmask(j,i,sr) * flag
                   sums_l(k,41,tn)  = sums_l(k,41,tn) +                        &
                                      q(k,j,i) * rmask(j,i,sr)   * flag
                ENDDO
             ENDDO
          ENDDO
          IF ( bulk_cloud_model )  THEN
             !$OMP DO
             DO  i = nxl, nxr
                DO  j =  nys, nyn
                   DO  k = nzb, nzt+1
                      flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 22 ) )
                      sums_l(k,42,tn) = sums_l(k,42,tn) +                      &
                                      ( q(k,j,i) - ql(k,j,i) ) * rmask(j,i,sr) &
                                                               * flag
                      sums_l(k,43,tn) = sums_l(k,43,tn) + (                    &
                                      pt(k,j,i) + lv_d_cp * d_exner(k) * ql(k,j,i) &
                                                          ) * rmask(j,i,sr)    &
                                                            * flag
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDIF

!
!--    Horizontally averaged profiles of passive scalar
       IF ( passive_scalar )  THEN
          !$OMP DO
          DO  i = nxl, nxr
             DO  j =  nys, nyn
                DO  k = nzb, nzt+1
                   sums_l(k,115,tn)  = sums_l(k,115,tn) + s(k,j,i)             &
                                    * rmask(j,i,sr)                            &
                                    * MERGE( 1.0_wp, 0.0_wp,                   &
                                             BTEST( wall_flags_total_0(k,j,i), 22 ) )
                ENDDO
             ENDDO
          ENDDO
       ENDIF
       !$OMP END PARALLEL
!
!--    Summation of thread sums
       IF ( threads_per_task > 1 )  THEN
          DO  i = 1, threads_per_task-1
             sums_l(:,1,0) = sums_l(:,1,0) + sums_l(:,1,i)
             sums_l(:,2,0) = sums_l(:,2,0) + sums_l(:,2,i)
             sums_l(:,4,0) = sums_l(:,4,0) + sums_l(:,4,i)
             IF ( ocean_mode )  THEN
                sums_l(:,23,0) = sums_l(:,23,0) + sums_l(:,23,i)
             ENDIF
             IF ( humidity )  THEN
                sums_l(:,41,0) = sums_l(:,41,0) + sums_l(:,41,i)
                sums_l(:,44,0) = sums_l(:,44,0) + sums_l(:,44,i)
                IF ( bulk_cloud_model )  THEN
                   sums_l(:,42,0) = sums_l(:,42,0) + sums_l(:,42,i)
                   sums_l(:,43,0) = sums_l(:,43,0) + sums_l(:,43,i)
                ENDIF
             ENDIF
             IF ( passive_scalar )  THEN
                sums_l(:,115,0) = sums_l(:,115,0) + sums_l(:,115,i)
             ENDIF
          ENDDO
       ENDIF

#if defined( __parallel )
!
!--    Compute total sum from local sums
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( sums_l(nzb,1,0), sums(nzb,1), nzt+2-nzb, MPI_REAL,  &
                           MPI_SUM, comm2d, ierr )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( sums_l(nzb,2,0), sums(nzb,2), nzt+2-nzb, MPI_REAL,  &
                           MPI_SUM, comm2d, ierr )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( sums_l(nzb,4,0), sums(nzb,4), nzt+2-nzb, MPI_REAL,  &
                           MPI_SUM, comm2d, ierr )
       IF ( ocean_mode )  THEN
          IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
          CALL MPI_ALLREDUCE( sums_l(nzb,23,0), sums(nzb,23), nzt+2-nzb,       &
                              MPI_REAL, MPI_SUM, comm2d, ierr )
       ENDIF
       IF ( humidity ) THEN
          IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
          CALL MPI_ALLREDUCE( sums_l(nzb,44,0), sums(nzb,44), nzt+2-nzb,       &
                              MPI_REAL, MPI_SUM, comm2d, ierr )
          IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
          CALL MPI_ALLREDUCE( sums_l(nzb,41,0), sums(nzb,41), nzt+2-nzb,       &
                              MPI_REAL, MPI_SUM, comm2d, ierr )
          IF ( bulk_cloud_model ) THEN
             IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
             CALL MPI_ALLREDUCE( sums_l(nzb,42,0), sums(nzb,42), nzt+2-nzb,    &
                                 MPI_REAL, MPI_SUM, comm2d, ierr )
             IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
             CALL MPI_ALLREDUCE( sums_l(nzb,43,0), sums(nzb,43), nzt+2-nzb,    &
                                 MPI_REAL, MPI_SUM, comm2d, ierr )
          ENDIF
       ENDIF

       IF ( passive_scalar )  THEN
          IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
          CALL MPI_ALLREDUCE( sums_l(nzb,115,0), sums(nzb,115), nzt+2-nzb,       &
                              MPI_REAL, MPI_SUM, comm2d, ierr )
       ENDIF
#else
       sums(:,1) = sums_l(:,1,0)
       sums(:,2) = sums_l(:,2,0)
       sums(:,4) = sums_l(:,4,0)
       IF ( ocean_mode )  sums(:,23) = sums_l(:,23,0)
       IF ( humidity ) THEN
          sums(:,44) = sums_l(:,44,0)
          sums(:,41) = sums_l(:,41,0)
          IF ( bulk_cloud_model ) THEN
             sums(:,42) = sums_l(:,42,0)
             sums(:,43) = sums_l(:,43,0)
          ENDIF
       ENDIF
       IF ( passive_scalar )  sums(:,115) = sums_l(:,115,0)
#endif

!
!--    Final values are obtained by division by the total number of grid points 
!--    used for summation. After that store profiles.
       sums(:,1) = sums(:,1) / ngp_2dh(sr)
       sums(:,2) = sums(:,2) / ngp_2dh(sr)
       sums(:,4) = sums(:,4) / ngp_2dh_s_inner(:,sr)
       hom(:,1,1,sr) = sums(:,1)             ! u
       hom(:,1,2,sr) = sums(:,2)             ! v
       hom(:,1,4,sr) = sums(:,4)             ! pt
       !$ACC UPDATE DEVICE(hom(:,1,1,sr), hom(:,1,2,sr), hom(:,1,4,sr))


!
!--    Salinity
       IF ( ocean_mode )  THEN
          sums(:,23) = sums(:,23) / ngp_2dh_s_inner(:,sr)
          hom(:,1,23,sr) = sums(:,23)             ! sa
       ENDIF

!
!--    Humidity and cloud parameters
       IF ( humidity ) THEN
          sums(:,44) = sums(:,44) / ngp_2dh_s_inner(:,sr)
          sums(:,41) = sums(:,41) / ngp_2dh_s_inner(:,sr)
          hom(:,1,44,sr) = sums(:,44)             ! vpt
          hom(:,1,41,sr) = sums(:,41)             ! qv (q)
          IF ( bulk_cloud_model ) THEN
             sums(:,42) = sums(:,42) / ngp_2dh_s_inner(:,sr)
             sums(:,43) = sums(:,43) / ngp_2dh_s_inner(:,sr)
             hom(:,1,42,sr) = sums(:,42)             ! qv
             hom(:,1,43,sr) = sums(:,43)             ! pt
          ENDIF
       ENDIF

!
!--    Passive scalar
       IF ( passive_scalar )  hom(:,1,115,sr) = sums(:,115) /                  &
            ngp_2dh_s_inner(:,sr)                    ! s 

!
!--    Horizontally averaged profiles of the remaining prognostic variables,
!--    variances, the total and the perturbation energy (single values in last
!--    column of sums_l) and some diagnostic quantities.
!--    NOTE: for simplicity, nzb_s_inner is used below, although strictly 
!--    ----  speaking the following k-loop would have to be split up and 
!--          rearranged according to the staggered grid.
!--          However, this implies no error since staggered velocity components
!--          are zero at the walls and inside buildings.
       tn = 0
       !$OMP PARALLEL PRIVATE( i, j, k, pts, sums_ll,                          &
       !$OMP                   sums_l_etot, tn, ust, ust2, u2, vst, vst2, v2,  &
       !$OMP                   w2, flag, m, ki, l )
       !$ tn = omp_get_thread_num()
       !$OMP DO
       !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(i, j, k, m) &
       !$ACC PRIVATE(sums_l_etot, flag) &
       !$ACC PRESENT(wall_flags_total_0, rmask, momentumflux_output_conversion) &
       !$ACC PRESENT(hom(:,1,4,sr)) &
       !$ACC PRESENT(e, u, v, w, km, kh, p, pt) &
       !$ACC PRESENT(surf_def_h(0), surf_lsm_h, surf_usm_h) &
       !$ACC PRESENT(sums_l)
       DO  i = nxl, nxr
          DO  j =  nys, nyn
             sums_l_etot = 0.0_wp
             DO  k = nzb, nzt+1
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 22 ) )
!
!--             Prognostic and diagnostic variables
                !$ACC ATOMIC
                sums_l(k,3,tn)  = sums_l(k,3,tn)  + w(k,j,i)  * rmask(j,i,sr)  &
                                                              * flag
                !$ACC ATOMIC
                sums_l(k,8,tn)  = sums_l(k,8,tn)  + e(k,j,i)  * rmask(j,i,sr)  &
                                                              * flag
                !$ACC ATOMIC
                sums_l(k,9,tn)  = sums_l(k,9,tn)  + km(k,j,i) * rmask(j,i,sr)  &
                                                              * flag
                !$ACC ATOMIC
                sums_l(k,10,tn) = sums_l(k,10,tn) + kh(k,j,i) * rmask(j,i,sr)  &
                                                              * flag
                !$ACC ATOMIC
                sums_l(k,40,tn) = sums_l(k,40,tn) + ( p(k,j,i)                 &
                                         / momentumflux_output_conversion(k) ) &
                                                              * flag

                !$ACC ATOMIC
                sums_l(k,33,tn) = sums_l(k,33,tn) + &
                                  ( pt(k,j,i)-hom(k,1,4,sr) )**2 * rmask(j,i,sr)&
                                                                 * flag
#ifndef _OPENACC
                IF ( humidity )  THEN
                   sums_l(k,70,tn) = sums_l(k,70,tn) + &
                                  ( q(k,j,i)-hom(k,1,41,sr) )**2 * rmask(j,i,sr)&
                                                                 * flag
                ENDIF
                IF ( passive_scalar )  THEN
                   sums_l(k,116,tn) = sums_l(k,116,tn) + &
                                  ( s(k,j,i)-hom(k,1,115,sr) )**2 * rmask(j,i,sr)&
                                                                  * flag
                ENDIF
#endif
!
!--             Higher moments
!--             (Computation of the skewness of w further below)
                !$ACC ATOMIC
                sums_l(k,38,tn) = sums_l(k,38,tn) + w(k,j,i)**3 * rmask(j,i,sr) &
                                                                * flag

                sums_l_etot  = sums_l_etot + &
                                        0.5_wp * ( u(k,j,i)**2 + v(k,j,i)**2 +  &
                                        w(k,j,i)**2 )            * rmask(j,i,sr)&
                                                                 * flag
             ENDDO
!
!--          Total and perturbation energy for the total domain (being
!--          collected in the last column of sums_l). Summation of these
!--          quantities is seperated from the previous loop in order to
!--          allow vectorization of that loop.
             !$ACC ATOMIC
             sums_l(nzb+4,pr_palm,tn) = sums_l(nzb+4,pr_palm,tn) + sums_l_etot
!
!--          2D-arrays (being collected in the last column of sums_l)
             IF ( surf_def_h(0)%end_index(j,i) >=                              &
                  surf_def_h(0)%start_index(j,i) )  THEN
                m = surf_def_h(0)%start_index(j,i)
                !$ACC ATOMIC
                sums_l(nzb,pr_palm,tn)   = sums_l(nzb,pr_palm,tn) +            &
                                        surf_def_h(0)%us(m)   * rmask(j,i,sr)
                !$ACC ATOMIC
                sums_l(nzb+1,pr_palm,tn) = sums_l(nzb+1,pr_palm,tn) +          &
                                        surf_def_h(0)%usws(m) * rmask(j,i,sr)
                !$ACC ATOMIC
                sums_l(nzb+2,pr_palm,tn) = sums_l(nzb+2,pr_palm,tn) +          &
                                        surf_def_h(0)%vsws(m) * rmask(j,i,sr)
                !$ACC ATOMIC
                sums_l(nzb+3,pr_palm,tn) = sums_l(nzb+3,pr_palm,tn) +          &
                                        surf_def_h(0)%ts(m)   * rmask(j,i,sr)
#ifndef _OPENACC
                IF ( humidity )  THEN
                   sums_l(nzb+12,pr_palm,tn) = sums_l(nzb+12,pr_palm,tn) +     &
                                            surf_def_h(0)%qs(m)   * rmask(j,i,sr)
                ENDIF
                IF ( passive_scalar )  THEN
                   sums_l(nzb+13,pr_palm,tn) = sums_l(nzb+13,pr_palm,tn) +     &
                                            surf_def_h(0)%ss(m)   * rmask(j,i,sr)
                ENDIF
#endif
!
!--             Summation of surface temperature.
                !$ACC ATOMIC
                sums_l(nzb+14,pr_palm,tn) = sums_l(nzb+14,pr_palm,tn)   +      &
                                            surf_def_h(0)%pt_surface(m) *      &
                                            rmask(j,i,sr)
             ENDIF
             IF ( surf_lsm_h%end_index(j,i) >= surf_lsm_h%start_index(j,i) )  THEN
                m = surf_lsm_h%start_index(j,i)
                !$ACC ATOMIC
                sums_l(nzb,pr_palm,tn)   = sums_l(nzb,pr_palm,tn) +            &
                                        surf_lsm_h%us(m)   * rmask(j,i,sr)
                !$ACC ATOMIC
                sums_l(nzb+1,pr_palm,tn) = sums_l(nzb+1,pr_palm,tn) +          &
                                        surf_lsm_h%usws(m) * rmask(j,i,sr)
                !$ACC ATOMIC
                sums_l(nzb+2,pr_palm,tn) = sums_l(nzb+2,pr_palm,tn) +          &
                                        surf_lsm_h%vsws(m) * rmask(j,i,sr)
                !$ACC ATOMIC
                sums_l(nzb+3,pr_palm,tn) = sums_l(nzb+3,pr_palm,tn) +          &
                                        surf_lsm_h%ts(m)   * rmask(j,i,sr)
#ifndef _OPENACC
                IF ( humidity )  THEN
                   sums_l(nzb+12,pr_palm,tn) = sums_l(nzb+12,pr_palm,tn) +     &
                                            surf_lsm_h%qs(m)   * rmask(j,i,sr)
                ENDIF
                IF ( passive_scalar )  THEN
                   sums_l(nzb+13,pr_palm,tn) = sums_l(nzb+13,pr_palm,tn) +     &
                                            surf_lsm_h%ss(m)   * rmask(j,i,sr)
                ENDIF
#endif
!
!--             Summation of surface temperature.
                !$ACC ATOMIC
                sums_l(nzb+14,pr_palm,tn) = sums_l(nzb+14,pr_palm,tn)   +      &
                                            surf_lsm_h%pt_surface(m)    *      &
                                            rmask(j,i,sr)
             ENDIF
             IF ( surf_usm_h%end_index(j,i) >= surf_usm_h%start_index(j,i) )  THEN
                m = surf_usm_h%start_index(j,i)
                !$ACC ATOMIC
                sums_l(nzb,pr_palm,tn)   = sums_l(nzb,pr_palm,tn) +            &
                                        surf_usm_h%us(m)   * rmask(j,i,sr)
                !$ACC ATOMIC
                sums_l(nzb+1,pr_palm,tn) = sums_l(nzb+1,pr_palm,tn) +          &
                                        surf_usm_h%usws(m) * rmask(j,i,sr)
                !$ACC ATOMIC
                sums_l(nzb+2,pr_palm,tn) = sums_l(nzb+2,pr_palm,tn) +          &
                                        surf_usm_h%vsws(m) * rmask(j,i,sr)
                !$ACC ATOMIC
                sums_l(nzb+3,pr_palm,tn) = sums_l(nzb+3,pr_palm,tn) +          &
                                        surf_usm_h%ts(m)   * rmask(j,i,sr)
#ifndef _OPENACC
                IF ( humidity )  THEN
                   sums_l(nzb+12,pr_palm,tn) = sums_l(nzb+12,pr_palm,tn) +     &
                                            surf_usm_h%qs(m)   * rmask(j,i,sr)
                ENDIF
                IF ( passive_scalar )  THEN
                   sums_l(nzb+13,pr_palm,tn) = sums_l(nzb+13,pr_palm,tn) +     &
                                            surf_usm_h%ss(m)   * rmask(j,i,sr)
                ENDIF
#endif
!
!--             Summation of surface temperature.
                !$ACC ATOMIC
                sums_l(nzb+14,pr_palm,tn) = sums_l(nzb+14,pr_palm,tn)   +      &
                                            surf_usm_h%pt_surface(m)    *      &
                                            rmask(j,i,sr)
             ENDIF
          ENDDO
       ENDDO
       !$ACC UPDATE &
       !$ACC HOST(sums_l(:,3,tn), sums_l(:,8,tn), sums_l(:,9,tn)) &
       !$ACC HOST(sums_l(:,10,tn), sums_l(:,40,tn), sums_l(:,33,tn)) &
       !$ACC HOST(sums_l(:,38,tn)) &
       !$ACC HOST(sums_l(nzb:nzb+4,pr_palm,tn), sums_l(nzb+14:nzb+14,pr_palm,tn))

!
!--    Computation of statistics when ws-scheme is not used. Else these 
!--    quantities are evaluated in the advection routines.
       IF ( .NOT. ws_scheme_mom .OR. sr /= 0 .OR. simulated_time == 0.0_wp )   &
       THEN
          !$OMP DO
          DO  i = nxl, nxr
             DO  j =  nys, nyn
                DO  k = nzb, nzt+1
                   flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 22 ) )

                   u2   = u(k,j,i)**2
                   v2   = v(k,j,i)**2
                   w2   = w(k,j,i)**2
                   ust2 = ( u(k,j,i) - hom(k,1,1,sr) )**2
                   vst2 = ( v(k,j,i) - hom(k,1,2,sr) )**2

                   sums_l(k,30,tn) = sums_l(k,30,tn) + ust2 * rmask(j,i,sr)    &
                                                            * flag
                   sums_l(k,31,tn) = sums_l(k,31,tn) + vst2 * rmask(j,i,sr)    &
                                                            * flag
                   sums_l(k,32,tn) = sums_l(k,32,tn) + w2   * rmask(j,i,sr)    &
                                                            * flag
!
!--                Perturbation energy

                   sums_l(k,34,tn) = sums_l(k,34,tn) + 0.5_wp *                &
                                  ( ust2 + vst2 + w2 )      * rmask(j,i,sr)    &
                                                            * flag
                ENDDO
             ENDDO
          ENDDO
       ENDIF
!
!--    Computaion of domain-averaged perturbation energy. Please note, 
!--    to prevent that perturbation energy is larger (even if only slightly) 
!--    than the total kinetic energy, calculation is based on deviations from
!--    the horizontal mean, instead of spatial descretization of the advection
!--    term. 
       !$OMP DO
       !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i, j, k, flag, w2, ust2, vst2) &
       !$ACC PRESENT(wall_flags_total_0, u, v, w, rmask, hom(:,1,1:2,sr)) &
       !$ACC PRESENT(sums_l)
       DO  i = nxl, nxr
          DO  j =  nys, nyn
             DO  k = nzb, nzt+1
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 22 ) )

                w2   = w(k,j,i)**2
                ust2 = ( u(k,j,i) - hom(k,1,1,sr) )**2
                vst2 = ( v(k,j,i) - hom(k,1,2,sr) )**2
                w2   = w(k,j,i)**2

                !$ACC ATOMIC
                sums_l(nzb+5,pr_palm,tn) = sums_l(nzb+5,pr_palm,tn)            &
                                 + 0.5_wp * ( ust2 + vst2 + w2 )               &
                                 * rmask(j,i,sr)                               &
                                 * flag
             ENDDO
          ENDDO
       ENDDO
       !$ACC UPDATE HOST(sums_l(nzb+5:nzb+5,pr_palm,tn))

!
!--    Horizontally averaged profiles of the vertical fluxes

       !$OMP DO
       !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(i, j, k, l, m) &
       !$ACC PRIVATE(ki, flag, ust, vst, pts) &
       !$ACC PRESENT(kh, km, u, v, w, pt) &
       !$ACC PRESENT(wall_flags_total_0, rmask, ddzu, rho_air_zw, hom(:,1,1:4,sr)) &
       !$ACC PRESENT(heatflux_output_conversion, momentumflux_output_conversion) &
       !$ACC PRESENT(surf_def_h(0:2), surf_lsm_h, surf_usm_h) &
       !$ACC PRESENT(sums_l)
       DO  i = nxl, nxr
          DO  j = nys, nyn
!
!--          Subgridscale fluxes (without Prandtl layer from k=nzb, 
!--          oterwise from k=nzb+1)
!--          NOTE: for simplicity, nzb_diff_s_inner is used below, although
!--          ----  strictly speaking the following k-loop would have to be
!--                split up according to the staggered grid.
!--                However, this implies no error since staggered velocity 
!--                components are zero at the walls and inside buildings.
!--          Flag 23 is used to mask surface fluxes as well as model-top fluxes,
!--          which are added further below. 
             DO  k = nzb, nzt
                flag = MERGE( 1.0_wp, 0.0_wp,                                  &
                              BTEST( wall_flags_total_0(k,j,i), 23 ) ) *       &
                       MERGE( 1.0_wp, 0.0_wp,                                  &
                              BTEST( wall_flags_total_0(k,j,i), 9  ) )
!
!--             Momentum flux w"u"
                !$ACC ATOMIC
                sums_l(k,12,tn) = sums_l(k,12,tn) - 0.25_wp * (                &
                               km(k,j,i)+km(k+1,j,i)+km(k,j,i-1)+km(k+1,j,i-1) &
                                                           ) * (               &
                                   ( u(k+1,j,i) - u(k,j,i)   ) * ddzu(k+1)     &
                                 + ( w(k,j,i)   - w(k,j,i-1) ) * ddx           &
                                                           ) * rmask(j,i,sr)   &
                                         * rho_air_zw(k)                       &
                                         * momentumflux_output_conversion(k)   &
                                         * flag
!
!--             Momentum flux w"v"
                !$ACC ATOMIC
                sums_l(k,14,tn) = sums_l(k,14,tn) - 0.25_wp * (                &
                               km(k,j,i)+km(k+1,j,i)+km(k,j-1,i)+km(k+1,j-1,i) &
                                                           ) * (               &
                                   ( v(k+1,j,i) - v(k,j,i)   ) * ddzu(k+1)     &
                                 + ( w(k,j,i)   - w(k,j-1,i) ) * ddy           &
                                                           ) * rmask(j,i,sr)   &
                                         * rho_air_zw(k)                       &
                                         * momentumflux_output_conversion(k)   &
                                         * flag
!
!--             Heat flux w"pt"
                !$ACC ATOMIC
                sums_l(k,16,tn) = sums_l(k,16,tn)                              &
                                         - 0.5_wp * ( kh(k,j,i) + kh(k+1,j,i) )&
                                               * ( pt(k+1,j,i) - pt(k,j,i) )   &
                                               * rho_air_zw(k)                 &
                                               * heatflux_output_conversion(k) &
                                               * ddzu(k+1) * rmask(j,i,sr)     &
                                               * flag

!
!--             Salinity flux w"sa"
#ifndef _OPENACC
                IF ( ocean_mode )  THEN
                   sums_l(k,65,tn) = sums_l(k,65,tn)                           &
                                         - 0.5_wp * ( kh(k,j,i) + kh(k+1,j,i) )&
                                               * ( sa(k+1,j,i) - sa(k,j,i) )   &
                                               * ddzu(k+1) * rmask(j,i,sr)     &
                                               * flag
                ENDIF

!
!--             Buoyancy flux, water flux (humidity flux) w"q"
                IF ( humidity ) THEN
                   sums_l(k,45,tn) = sums_l(k,45,tn)                           &
                                         - 0.5_wp * ( kh(k,j,i) + kh(k+1,j,i) )&
                                               * ( vpt(k+1,j,i) - vpt(k,j,i) ) &
                                               * rho_air_zw(k)                 &
                                               * heatflux_output_conversion(k) &
                                               * ddzu(k+1) * rmask(j,i,sr) * flag
                   sums_l(k,48,tn) = sums_l(k,48,tn)                           &
                                         - 0.5_wp * ( kh(k,j,i) + kh(k+1,j,i) )&
                                               * ( q(k+1,j,i) - q(k,j,i) )     &
                                               * rho_air_zw(k)                 &
                                               * waterflux_output_conversion(k)&
                                               * ddzu(k+1) * rmask(j,i,sr) * flag

                   IF ( bulk_cloud_model ) THEN
                      sums_l(k,51,tn) = sums_l(k,51,tn)                        &
                                         - 0.5_wp * ( kh(k,j,i) + kh(k+1,j,i) )&
                                               * ( ( q(k+1,j,i) - ql(k+1,j,i) )&
                                                - ( q(k,j,i) - ql(k,j,i) ) )   &
                                               * rho_air_zw(k)                 &
                                               * waterflux_output_conversion(k)&
                                               * ddzu(k+1) * rmask(j,i,sr) * flag
                   ENDIF
                ENDIF

!
!--             Passive scalar flux
                IF ( passive_scalar )  THEN
                   sums_l(k,117,tn) = sums_l(k,117,tn)                         &
                                         - 0.5_wp * ( kh(k,j,i) + kh(k+1,j,i) )&
                                                  * ( s(k+1,j,i) - s(k,j,i) )  &
                                                  * ddzu(k+1) * rmask(j,i,sr)  &
                                                              * flag
                ENDIF
#endif

             ENDDO

!
!--          Subgridscale fluxes in the Prandtl layer
             IF ( use_surface_fluxes )  THEN
                DO  l = 0, 1
                   ! The original code using MERGE doesn't work with the PGI
                   ! compiler when running on the GPU.
                   ! This is submitted as a compiler Bug in PGI ticket TPR#26718
                   ! ki = MERGE( -1, 0, l == 0 )
                   ki = -1 + l
                   IF ( surf_def_h(l)%ns >= 1 )  THEN
                      DO  m = surf_def_h(l)%start_index(j,i),                  &
                              surf_def_h(l)%end_index(j,i)
                         k = surf_def_h(l)%k(m)

                         !$ACC ATOMIC
                         sums_l(k+ki,12,tn) = sums_l(k+ki,12,tn) + &
                                    momentumflux_output_conversion(k+ki) * &
                                    surf_def_h(l)%usws(m) * rmask(j,i,sr)     ! w"u"
                         !$ACC ATOMIC
                         sums_l(k+ki,14,tn) = sums_l(k+ki,14,tn) + &
                                    momentumflux_output_conversion(k+ki) * &
                                    surf_def_h(l)%vsws(m) * rmask(j,i,sr)     ! w"v"
                         !$ACC ATOMIC
                         sums_l(k+ki,16,tn) = sums_l(k+ki,16,tn) + &
                                    heatflux_output_conversion(k+ki) * &
                                    surf_def_h(l)%shf(m)  * rmask(j,i,sr)     ! w"pt"
#if 0
                         sums_l(k+ki,58,tn) = sums_l(k+ki,58,tn) + &
                                    0.0_wp * rmask(j,i,sr)        ! u"pt"
                         sums_l(k+ki,61,tn) = sums_l(k+ki,61,tn) + &
                                    0.0_wp * rmask(j,i,sr)        ! v"pt"
#endif
#ifndef _OPENACC
                         IF ( ocean_mode )  THEN
                            sums_l(k+ki,65,tn) = sums_l(k+ki,65,tn) + &
                                       surf_def_h(l)%sasws(m) * rmask(j,i,sr)  ! w"sa"
                         ENDIF
                         IF ( humidity )  THEN
                            sums_l(k+ki,48,tn) = sums_l(k+ki,48,tn) +                     &
                                       waterflux_output_conversion(k+ki) *      &
                                       surf_def_h(l)%qsws(m) * rmask(j,i,sr)  ! w"q" (w"qv")
                            sums_l(k+ki,45,tn) = sums_l(k+ki,45,tn) + (                   &
                                       ( 1.0_wp + 0.61_wp * q(k+ki,j,i) ) *     &
                                       surf_def_h(l)%shf(m) + 0.61_wp * pt(k+ki,j,i) *      &
                                                  surf_def_h(l)%qsws(m) )                  &
                                       * heatflux_output_conversion(k+ki)
                            IF ( cloud_droplets )  THEN
                               sums_l(k+ki,45,tn) = sums_l(k+ki,45,tn) + (                &
                                         ( 1.0_wp + 0.61_wp * q(k+ki,j,i) -     &
                                           ql(k+ki,j,i) ) * surf_def_h(l)%shf(m) +          &
                                           0.61_wp * pt(k+ki,j,i) * surf_def_h(l)%qsws(m) ) &
                                          * heatflux_output_conversion(k+ki)
                            ENDIF
                            IF ( bulk_cloud_model )  THEN
!
!--                            Formula does not work if ql(k+ki) /= 0.0
                               sums_l(k+ki,51,tn) = sums_l(k+ki,51,tn) +                  &
                                          waterflux_output_conversion(k+ki) *   &
                                          surf_def_h(l)%qsws(m) * rmask(j,i,sr) ! w"q" (w"qv")
                            ENDIF
                         ENDIF
                         IF ( passive_scalar )  THEN
                            sums_l(k+ki,117,tn) = sums_l(k+ki,117,tn) +                     &
                                        surf_def_h(l)%ssws(m) * rmask(j,i,sr) ! w"s"
                         ENDIF
#endif

                      ENDDO

                   ENDIF
                ENDDO
                IF ( surf_lsm_h%end_index(j,i) >=                              &
                     surf_lsm_h%start_index(j,i) )  THEN
                   m = surf_lsm_h%start_index(j,i)
                   !$ACC ATOMIC
                   sums_l(nzb,12,tn) = sums_l(nzb,12,tn) + &
                                    momentumflux_output_conversion(nzb) * &
                                    surf_lsm_h%usws(m) * rmask(j,i,sr)     ! w"u"
                   !$ACC ATOMIC
                   sums_l(nzb,14,tn) = sums_l(nzb,14,tn) + &
                                    momentumflux_output_conversion(nzb) * &
                                    surf_lsm_h%vsws(m) * rmask(j,i,sr)     ! w"v"
                   !$ACC ATOMIC
                   sums_l(nzb,16,tn) = sums_l(nzb,16,tn) + &
                                    heatflux_output_conversion(nzb) * &
                                    surf_lsm_h%shf(m)  * rmask(j,i,sr)     ! w"pt"
#if 0
                   sums_l(nzb,58,tn) = sums_l(nzb,58,tn) + &
                                    0.0_wp * rmask(j,i,sr)        ! u"pt"
                   sums_l(nzb,61,tn) = sums_l(nzb,61,tn) + &
                                    0.0_wp * rmask(j,i,sr)        ! v"pt"
#endif
#ifndef _OPENACC
                   IF ( ocean_mode )  THEN
                      sums_l(nzb,65,tn) = sums_l(nzb,65,tn) + &
                                       surf_lsm_h%sasws(m) * rmask(j,i,sr)  ! w"sa"
                   ENDIF
                   IF ( humidity )  THEN
                      sums_l(nzb,48,tn) = sums_l(nzb,48,tn) +                     &
                                       waterflux_output_conversion(nzb) *      &
                                       surf_lsm_h%qsws(m) * rmask(j,i,sr)  ! w"q" (w"qv")
                      sums_l(nzb,45,tn) = sums_l(nzb,45,tn) + (                   &
                                       ( 1.0_wp + 0.61_wp * q(nzb,j,i) ) *     &
                                       surf_lsm_h%shf(m) + 0.61_wp * pt(nzb,j,i) *      &
                                                  surf_lsm_h%qsws(m) )                  &
                                       * heatflux_output_conversion(nzb)
                      IF ( cloud_droplets )  THEN
                         sums_l(nzb,45,tn) = sums_l(nzb,45,tn) + (                &
                                         ( 1.0_wp + 0.61_wp * q(nzb,j,i) -     &
                                           ql(nzb,j,i) ) * surf_lsm_h%shf(m) +          &
                                           0.61_wp * pt(nzb,j,i) * surf_lsm_h%qsws(m) ) &
                                          * heatflux_output_conversion(nzb)
                      ENDIF
                      IF ( bulk_cloud_model )  THEN
!
!--                      Formula does not work if ql(nzb) /= 0.0
                         sums_l(nzb,51,tn) = sums_l(nzb,51,tn) +                  &
                                          waterflux_output_conversion(nzb) *   &
                                          surf_lsm_h%qsws(m) * rmask(j,i,sr) ! w"q" (w"qv")
                      ENDIF
                   ENDIF
                   IF ( passive_scalar )  THEN
                      sums_l(nzb,117,tn) = sums_l(nzb,117,tn) +                     &
                                        surf_lsm_h%ssws(m) * rmask(j,i,sr) ! w"s"
                   ENDIF
#endif

                ENDIF
                IF ( surf_usm_h%end_index(j,i) >=                              &
                     surf_usm_h%start_index(j,i) )  THEN
                   m = surf_usm_h%start_index(j,i)
                   !$ACC ATOMIC
                   sums_l(nzb,12,tn) = sums_l(nzb,12,tn) + &
                                    momentumflux_output_conversion(nzb) * &
                                    surf_usm_h%usws(m) * rmask(j,i,sr)     ! w"u"
                   !$ACC ATOMIC
                   sums_l(nzb,14,tn) = sums_l(nzb,14,tn) + &
                                    momentumflux_output_conversion(nzb) * &
                                    surf_usm_h%vsws(m) * rmask(j,i,sr)     ! w"v"
                   !$ACC ATOMIC
                   sums_l(nzb,16,tn) = sums_l(nzb,16,tn) + &
                                    heatflux_output_conversion(nzb) * &
                                    surf_usm_h%shf(m)  * rmask(j,i,sr)     ! w"pt"
#if 0
                   sums_l(nzb,58,tn) = sums_l(nzb,58,tn) + &
                                    0.0_wp * rmask(j,i,sr)        ! u"pt"
                   sums_l(nzb,61,tn) = sums_l(nzb,61,tn) + &
                                    0.0_wp * rmask(j,i,sr)        ! v"pt"
#endif
#ifndef _OPENACC
                   IF ( ocean_mode )  THEN
                      sums_l(nzb,65,tn) = sums_l(nzb,65,tn) + &
                                       surf_usm_h%sasws(m) * rmask(j,i,sr)  ! w"sa"
                   ENDIF
                   IF ( humidity )  THEN
                      sums_l(nzb,48,tn) = sums_l(nzb,48,tn) +                     &
                                       waterflux_output_conversion(nzb) *      &
                                       surf_usm_h%qsws(m) * rmask(j,i,sr)  ! w"q" (w"qv")
                      sums_l(nzb,45,tn) = sums_l(nzb,45,tn) + (                   &
                                       ( 1.0_wp + 0.61_wp * q(nzb,j,i) ) *     &
                                       surf_usm_h%shf(m) + 0.61_wp * pt(nzb,j,i) *      &
                                                  surf_usm_h%qsws(m) )                  &
                                       * heatflux_output_conversion(nzb)
                      IF ( cloud_droplets )  THEN
                         sums_l(nzb,45,tn) = sums_l(nzb,45,tn) + (                &
                                         ( 1.0_wp + 0.61_wp * q(nzb,j,i) -     &
                                           ql(nzb,j,i) ) * surf_usm_h%shf(m) +          &
                                           0.61_wp * pt(nzb,j,i) * surf_usm_h%qsws(m) ) &
                                          * heatflux_output_conversion(nzb)
                      ENDIF
                      IF ( bulk_cloud_model )  THEN
!
!--                      Formula does not work if ql(nzb) /= 0.0
                         sums_l(nzb,51,tn) = sums_l(nzb,51,tn) +                  &
                                          waterflux_output_conversion(nzb) *   &
                                          surf_usm_h%qsws(m) * rmask(j,i,sr) ! w"q" (w"qv")
                      ENDIF
                   ENDIF
                   IF ( passive_scalar )  THEN
                      sums_l(nzb,117,tn) = sums_l(nzb,117,tn) +                     &
                                        surf_usm_h%ssws(m) * rmask(j,i,sr) ! w"s"
                   ENDIF
#endif

                ENDIF

             ENDIF

#ifndef _OPENACC
             IF ( .NOT. neutral )  THEN
                IF ( surf_def_h(0)%end_index(j,i) >=                           &
                     surf_def_h(0)%start_index(j,i) )  THEN
                   m = surf_def_h(0)%start_index(j,i)
                   sums_l(nzb,112,tn) = sums_l(nzb,112,tn) +                   &
                                        surf_def_h(0)%ol(m)  * rmask(j,i,sr) ! L
                ENDIF
                IF ( surf_lsm_h%end_index(j,i) >=                              &
                     surf_lsm_h%start_index(j,i) )  THEN
                   m = surf_lsm_h%start_index(j,i)
                   sums_l(nzb,112,tn) = sums_l(nzb,112,tn) +                   &
                                        surf_lsm_h%ol(m)  * rmask(j,i,sr) ! L
                ENDIF
                IF ( surf_usm_h%end_index(j,i) >=                              &
                     surf_usm_h%start_index(j,i) )  THEN
                   m = surf_usm_h%start_index(j,i)
                   sums_l(nzb,112,tn) = sums_l(nzb,112,tn) +                   &
                                        surf_usm_h%ol(m)  * rmask(j,i,sr) ! L
                ENDIF
             ENDIF

             IF ( radiation )  THEN
                IF ( surf_def_h(0)%end_index(j,i) >=                           &
                     surf_def_h(0)%start_index(j,i) )  THEN
                   m = surf_def_h(0)%start_index(j,i)
                   sums_l(nzb,99,tn)  = sums_l(nzb,99,tn)  +                   &
                                        surf_def_h(0)%rad_net(m) * rmask(j,i,sr)
                   sums_l(nzb,100,tn) = sums_l(nzb,100,tn)  +                  &
                                        surf_def_h(0)%rad_lw_in(m) * rmask(j,i,sr)
                   sums_l(nzb,101,tn) = sums_l(nzb,101,tn)  +                  &
                                        surf_def_h(0)%rad_lw_out(m) * rmask(j,i,sr)
                   sums_l(nzb,102,tn) = sums_l(nzb,102,tn)  +                  &
                                        surf_def_h(0)%rad_sw_in(m) * rmask(j,i,sr)
                   sums_l(nzb,103,tn) = sums_l(nzb,103,tn)  +                  &
                                        surf_def_h(0)%rad_sw_out(m) * rmask(j,i,sr)
                ENDIF
                IF ( surf_lsm_h%end_index(j,i) >=                              &
                     surf_lsm_h%start_index(j,i) )  THEN
                   m = surf_lsm_h%start_index(j,i)
                   sums_l(nzb,99,tn)  = sums_l(nzb,99,tn)  +                   &
                                        surf_lsm_h%rad_net(m) * rmask(j,i,sr)
                   sums_l(nzb,100,tn) = sums_l(nzb,100,tn)  +                  &
                                        surf_lsm_h%rad_lw_in(m) * rmask(j,i,sr)
                   sums_l(nzb,101,tn) = sums_l(nzb,101,tn)  +                  &
                                        surf_lsm_h%rad_lw_out(m) * rmask(j,i,sr)
                   sums_l(nzb,102,tn) = sums_l(nzb,102,tn)  +                  &
                                        surf_lsm_h%rad_sw_in(m) * rmask(j,i,sr)
                   sums_l(nzb,103,tn) = sums_l(nzb,103,tn)  +                  &
                                        surf_lsm_h%rad_sw_out(m) * rmask(j,i,sr)
                ENDIF
                IF ( surf_usm_h%end_index(j,i) >=                              &
                     surf_usm_h%start_index(j,i) )  THEN
                   m = surf_usm_h%start_index(j,i)
                   sums_l(nzb,99,tn)  = sums_l(nzb,99,tn)  +                   &
                                        surf_usm_h%rad_net(m) * rmask(j,i,sr)
                   sums_l(nzb,100,tn) = sums_l(nzb,100,tn)  +                  &
                                        surf_usm_h%rad_lw_in(m) * rmask(j,i,sr)
                   sums_l(nzb,101,tn) = sums_l(nzb,101,tn)  +                  &
                                        surf_usm_h%rad_lw_out(m) * rmask(j,i,sr)
                   sums_l(nzb,102,tn) = sums_l(nzb,102,tn)  +                  &
                                        surf_usm_h%rad_sw_in(m) * rmask(j,i,sr)
                   sums_l(nzb,103,tn) = sums_l(nzb,103,tn)  +                  &
                                        surf_usm_h%rad_sw_out(m) * rmask(j,i,sr)
                ENDIF

#if defined ( __rrtmg )
                IF ( radiation_scheme == 'rrtmg' )  THEN

                   IF ( surf_def_h(0)%end_index(j,i) >=                        &
                        surf_def_h(0)%start_index(j,i) )  THEN
                      m = surf_def_h(0)%start_index(j,i)
                      sums_l(nzb,108,tn)  = sums_l(nzb,108,tn)  +              &
                                   surf_def_h(0)%rrtm_aldif(m,0) * rmask(j,i,sr)
                      sums_l(nzb,109,tn) = sums_l(nzb,109,tn)  +               &
                                   surf_def_h(0)%rrtm_aldir(m,0) * rmask(j,i,sr)
                      sums_l(nzb,110,tn) = sums_l(nzb,110,tn)  +               &
                                   surf_def_h(0)%rrtm_asdif(m,0) * rmask(j,i,sr)
                      sums_l(nzb,111,tn) = sums_l(nzb,111,tn)  +               &
                                   surf_def_h(0)%rrtm_asdir(m,0) * rmask(j,i,sr)
                   ENDIF
                   IF ( surf_lsm_h%end_index(j,i) >=                           &
                        surf_lsm_h%start_index(j,i) )  THEN
                      m = surf_lsm_h%start_index(j,i)
                      sums_l(nzb,108,tn)  = sums_l(nzb,108,tn)  +              &
                               SUM( surf_lsm_h%frac(m,:) *                     &
                                    surf_lsm_h%rrtm_aldif(m,:) ) * rmask(j,i,sr)
                      sums_l(nzb,109,tn) = sums_l(nzb,109,tn)  +               &
                               SUM( surf_lsm_h%frac(m,:) *                     &
                                    surf_lsm_h%rrtm_aldir(m,:) ) * rmask(j,i,sr)
                      sums_l(nzb,110,tn) = sums_l(nzb,110,tn)  +               &
                               SUM( surf_lsm_h%frac(m,:) *                     &
                                    surf_lsm_h%rrtm_asdif(m,:) ) * rmask(j,i,sr)
                      sums_l(nzb,111,tn) = sums_l(nzb,111,tn)  +               &
                               SUM( surf_lsm_h%frac(m,:) *                     &
                                    surf_lsm_h%rrtm_asdir(m,:) ) * rmask(j,i,sr)
                   ENDIF
                   IF ( surf_usm_h%end_index(j,i) >=                           &
                        surf_usm_h%start_index(j,i) )  THEN
                      m = surf_usm_h%start_index(j,i)
                      sums_l(nzb,108,tn)  = sums_l(nzb,108,tn)  +              &
                               SUM( surf_usm_h%frac(m,:) *                     &
                                    surf_usm_h%rrtm_aldif(m,:) ) * rmask(j,i,sr)
                      sums_l(nzb,109,tn) = sums_l(nzb,109,tn)  +               &
                               SUM( surf_usm_h%frac(m,:) *                     &
                                    surf_usm_h%rrtm_aldir(m,:) ) * rmask(j,i,sr)
                      sums_l(nzb,110,tn) = sums_l(nzb,110,tn)  +               &
                               SUM( surf_usm_h%frac(m,:) *                     &
                                    surf_usm_h%rrtm_asdif(m,:) ) * rmask(j,i,sr)
                      sums_l(nzb,111,tn) = sums_l(nzb,111,tn)  +               &
                               SUM( surf_usm_h%frac(m,:) *                     &
                                    surf_usm_h%rrtm_asdir(m,:) ) * rmask(j,i,sr)
                   ENDIF

                ENDIF
#endif
             ENDIF
#endif
!
!--          Subgridscale fluxes at the top surface
             IF ( use_top_fluxes )  THEN
                m = surf_def_h(2)%start_index(j,i)
                !$ACC ATOMIC
                sums_l(nzt,12,tn) = sums_l(nzt,12,tn) + &
                                    momentumflux_output_conversion(nzt) * &
                                    surf_def_h(2)%usws(m) * rmask(j,i,sr)    ! w"u"
                !$ACC ATOMIC
                sums_l(nzt+1,12,tn) = sums_l(nzt+1,12,tn) + &
                                    momentumflux_output_conversion(nzt+1) * &
                                    surf_def_h(2)%usws(m) * rmask(j,i,sr)    ! w"u"
                !$ACC ATOMIC
                sums_l(nzt,14,tn) = sums_l(nzt,14,tn) + &
                                    momentumflux_output_conversion(nzt) * &
                                    surf_def_h(2)%vsws(m) * rmask(j,i,sr)    ! w"v"
                !$ACC ATOMIC
                sums_l(nzt+1,14,tn) = sums_l(nzt+1,14,tn) + &
                                    momentumflux_output_conversion(nzt+1) * &
                                    surf_def_h(2)%vsws(m) * rmask(j,i,sr)    ! w"v"
                !$ACC ATOMIC
                sums_l(nzt,16,tn) = sums_l(nzt,16,tn) + &
                                    heatflux_output_conversion(nzt) * &
                                    surf_def_h(2)%shf(m)  * rmask(j,i,sr)   ! w"pt"
                !$ACC ATOMIC
                sums_l(nzt+1,16,tn) = sums_l(nzt+1,16,tn) + &
                                    heatflux_output_conversion(nzt+1) * &
                                    surf_def_h(2)%shf(m)  * rmask(j,i,sr)   ! w"pt"
#if 0
                sums_l(nzt:nzt+1,58,tn) = sums_l(nzt:nzt+1,58,tn) + &
                                    0.0_wp * rmask(j,i,sr)        ! u"pt"
                sums_l(nzt:nzt+1,61,tn) = sums_l(nzt:nzt+1,61,tn) + &
                                    0.0_wp * rmask(j,i,sr)        ! v"pt" 
#endif
#ifndef _OPENACC
                IF ( ocean_mode )  THEN
                   sums_l(nzt,65,tn) = sums_l(nzt,65,tn) + &
                                       surf_def_h(2)%sasws(m) * rmask(j,i,sr)  ! w"sa"
                ENDIF
                IF ( humidity )  THEN
                   sums_l(nzt,48,tn) = sums_l(nzt,48,tn) +                     &
                                       waterflux_output_conversion(nzt) *      &
                                       surf_def_h(2)%qsws(m) * rmask(j,i,sr) ! w"q" (w"qv")
                   sums_l(nzt,45,tn) = sums_l(nzt,45,tn) + (                   &
                                       ( 1.0_wp + 0.61_wp * q(nzt,j,i) ) *     &
                                       surf_def_h(2)%shf(m) +                  &
                                       0.61_wp * pt(nzt,j,i) *    &
                                       surf_def_h(2)%qsws(m) )      &
                                       * heatflux_output_conversion(nzt)
                   IF ( cloud_droplets )  THEN
                      sums_l(nzt,45,tn) = sums_l(nzt,45,tn) + (                &
                                          ( 1.0_wp + 0.61_wp * q(nzt,j,i) -    &
                                            ql(nzt,j,i) ) *                    &
                                            surf_def_h(2)%shf(m) +             &
                                           0.61_wp * pt(nzt,j,i) *             &
                                           surf_def_h(2)%qsws(m) )&
                                           * heatflux_output_conversion(nzt)
                   ENDIF
                   IF ( bulk_cloud_model )  THEN
!
!--                   Formula does not work if ql(nzb) /= 0.0
                      sums_l(nzt,51,tn) = sums_l(nzt,51,tn) + &   ! w"q" (w"qv")
                                          waterflux_output_conversion(nzt) *   &
                                          surf_def_h(2)%qsws(m) * rmask(j,i,sr)
                   ENDIF
                ENDIF
                IF ( passive_scalar )  THEN
                   sums_l(nzt,117,tn) = sums_l(nzt,117,tn) + &
                                        surf_def_h(2)%ssws(m) * rmask(j,i,sr) ! w"s"
                ENDIF
#endif
             ENDIF

!
!--          Resolved fluxes (can be computed for all horizontal points)
!--          NOTE: for simplicity, nzb_s_inner is used below, although strictly 
!--          ----  speaking the following k-loop would have to be split up and 
!--                rearranged according to the staggered grid.
             DO  k = nzb, nzt
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 22 ) )
                ust = 0.5_wp * ( u(k,j,i)   - hom(k,1,1,sr) +                  &
                                 u(k+1,j,i) - hom(k+1,1,1,sr) )
                vst = 0.5_wp * ( v(k,j,i)   - hom(k,1,2,sr) +                  &
                                 v(k+1,j,i) - hom(k+1,1,2,sr) )
                pts = 0.5_wp * ( pt(k,j,i)   - hom(k,1,4,sr) +                 &
                                 pt(k+1,j,i) - hom(k+1,1,4,sr) )

!--             Higher moments
                !$ACC ATOMIC
                sums_l(k,35,tn) = sums_l(k,35,tn) + pts * w(k,j,i)**2 *        &
                                                    rmask(j,i,sr) * flag
                !$ACC ATOMIC
                sums_l(k,36,tn) = sums_l(k,36,tn) + pts**2 * w(k,j,i) *        &
                                                    rmask(j,i,sr) * flag

!
!--             Salinity flux and density (density does not belong to here,
!--             but so far there is no other suitable place to calculate)
#ifndef _OPENACC
                IF ( ocean_mode )  THEN
                   IF( .NOT. ws_scheme_sca .OR. sr /= 0 )  THEN
                      pts = 0.5_wp * ( sa(k,j,i)   - hom(k,1,23,sr) +          &
                                       sa(k+1,j,i) - hom(k+1,1,23,sr) )
                      sums_l(k,66,tn) = sums_l(k,66,tn) + pts * w(k,j,i) *     &
                                        rmask(j,i,sr) * flag
                   ENDIF
                   sums_l(k,64,tn) = sums_l(k,64,tn) + rho_ocean(k,j,i) *      &
                                                       rmask(j,i,sr) * flag
                   sums_l(k,71,tn) = sums_l(k,71,tn) + prho(k,j,i) *           &
                                                       rmask(j,i,sr) * flag
                ENDIF

!
!--             Buoyancy flux, water flux, humidity flux, liquid water
!--             content, rain drop concentration and rain water content
                IF ( humidity )  THEN
                   IF ( bulk_cloud_model .OR. cloud_droplets )  THEN
                      pts = 0.5_wp * ( vpt(k,j,i)   - hom(k,1,44,sr) +         &
                                    vpt(k+1,j,i) - hom(k+1,1,44,sr) )
                      sums_l(k,46,tn) = sums_l(k,46,tn) + pts * w(k,j,i) *     &
                                               rho_air_zw(k) *                 &
                                               heatflux_output_conversion(k) * &
                                                          rmask(j,i,sr) * flag
                      sums_l(k,54,tn) = sums_l(k,54,tn) + ql(k,j,i) * rmask(j,i,sr) &
                                                                    * flag

                      IF ( .NOT. cloud_droplets )  THEN
                         pts = 0.5_wp *                                        &
                              ( ( q(k,j,i) - ql(k,j,i) ) -                     &
                              hom(k,1,42,sr) +                                 &
                              ( q(k+1,j,i) - ql(k+1,j,i) ) -                   &
                              hom(k+1,1,42,sr) )
                         sums_l(k,52,tn) = sums_l(k,52,tn) + pts * w(k,j,i) *  &
                                             rho_air_zw(k) *                   &
                                             waterflux_output_conversion(k) *  &
                                                             rmask(j,i,sr)  *  &
                                                             flag
                         sums_l(k,75,tn) = sums_l(k,75,tn) + qc(k,j,i) *       &
                                                             rmask(j,i,sr) *   &
                                                             flag
                         sums_l(k,76,tn) = sums_l(k,76,tn) + prr(k,j,i) *      &
                                                             rmask(j,i,sr) *   &
                                                             flag
                         IF ( microphysics_morrison )  THEN
                            sums_l(k,123,tn) = sums_l(k,123,tn) + nc(k,j,i) *  &
                                                                rmask(j,i,sr) *&
                                                                flag
                         ENDIF
                         IF ( microphysics_seifert )  THEN
                            sums_l(k,73,tn) = sums_l(k,73,tn) + nr(k,j,i) *    &
                                                                rmask(j,i,sr) *&
                                                                flag 
                            sums_l(k,74,tn) = sums_l(k,74,tn) + qr(k,j,i) *    &
                                                                rmask(j,i,sr) *&
                                                                flag
                         ENDIF
                      ENDIF

                   ELSE
                      IF( .NOT. ws_scheme_sca .OR. sr /= 0 )  THEN
                         pts = 0.5_wp * ( vpt(k,j,i)   - hom(k,1,44,sr) +      &
                                          vpt(k+1,j,i) - hom(k+1,1,44,sr) )
                         sums_l(k,46,tn) = sums_l(k,46,tn) + pts * w(k,j,i) *  &
                                              rho_air_zw(k) *                  &
                                              heatflux_output_conversion(k) *  &
                                                             rmask(j,i,sr)  *  &
                                                             flag
                      ELSE IF ( ws_scheme_sca .AND. sr == 0 )  THEN
                         sums_l(k,46,tn) = ( ( 1.0_wp + 0.61_wp *              & 
                                               hom(k,1,41,sr) ) *              &
                                             sums_l(k,17,tn) +                 &
                                             0.61_wp * hom(k,1,4,sr) *         &
                                             sums_l(k,49,tn)                   &
                                           ) * heatflux_output_conversion(k) * &
                                               flag
                      END IF
                   END IF
                ENDIF
!
!--             Passive scalar flux
                IF ( passive_scalar .AND. ( .NOT. ws_scheme_sca                &
                     .OR. sr /= 0 ) )  THEN
                   pts = 0.5_wp * ( s(k,j,i)   - hom(k,1,115,sr) +             &
                                    s(k+1,j,i) - hom(k+1,1,115,sr) )
                   sums_l(k,114,tn) = sums_l(k,114,tn) + pts * w(k,j,i) *      &
                                                       rmask(j,i,sr) * flag
                ENDIF
#endif

!
!--             Energy flux w*e*
!--             has to be adjusted
                !$ACC ATOMIC
                sums_l(k,37,tn) = sums_l(k,37,tn) + w(k,j,i) * 0.5_wp *        &
                                             ( ust**2 + vst**2 + w(k,j,i)**2 ) &
                                           * rho_air_zw(k)                     &
                                           * momentumflux_output_conversion(k) &
                                           * rmask(j,i,sr) * flag
             ENDDO
          ENDDO
       ENDDO
       !$OMP END PARALLEL

       !$ACC UPDATE &
       !$ACC HOST(sums_l(:,12,tn), sums_l(:,14,tn), sums_l(:,16,tn)) &
       !$ACC HOST(sums_l(:,35,tn), sums_l(:,36,tn), sums_l(:,37,tn))
!
!--    Treat land-surface quantities according to new wall model structure. 
       IF ( land_surface )  THEN
          tn = 0
          !$OMP PARALLEL PRIVATE( i, j, m, tn )
          !$ tn = omp_get_thread_num()
          !$OMP DO
          DO  m = 1, surf_lsm_h%ns
             i = surf_lsm_h%i(m)
             j = surf_lsm_h%j(m)
        
             IF ( i >= nxl  .AND.  i <= nxr  .AND.                             &
                  j >= nys  .AND.  j <= nyn )  THEN 
                sums_l(nzb,93,tn)  = sums_l(nzb,93,tn) + surf_lsm_h%ghf(m)
                sums_l(nzb,94,tn)  = sums_l(nzb,94,tn) + surf_lsm_h%qsws_liq(m)
                sums_l(nzb,95,tn)  = sums_l(nzb,95,tn) + surf_lsm_h%qsws_soil(m)
                sums_l(nzb,96,tn)  = sums_l(nzb,96,tn) + surf_lsm_h%qsws_veg(m)
                sums_l(nzb,97,tn)  = sums_l(nzb,97,tn) + surf_lsm_h%r_a(m)
                sums_l(nzb,98,tn) = sums_l(nzb,98,tn)+ surf_lsm_h%r_s(m)
             ENDIF
          ENDDO
          !$OMP END PARALLEL

          tn = 0
          !$OMP PARALLEL PRIVATE( i, j, k, m, tn )
          !$ tn = omp_get_thread_num()
          !$OMP DO
          DO  m = 1, surf_lsm_h%ns

             i = surf_lsm_h%i(m)            
             j = surf_lsm_h%j(m)

             IF ( i >= nxl  .AND.  i <= nxr  .AND.                             &
                  j >= nys  .AND.  j <= nyn )  THEN 

                DO  k = nzb_soil, nzt_soil
                   sums_l(k,89,tn)  = sums_l(k,89,tn)  + t_soil_h%var_2d(k,m)  &
                                      * rmask(j,i,sr)
                   sums_l(k,91,tn)  = sums_l(k,91,tn)  + m_soil_h%var_2d(k,m)  &
                                      * rmask(j,i,sr)
                ENDDO
             ENDIF
          ENDDO
          !$OMP END PARALLEL
       ENDIF
!
!--    For speed optimization fluxes which have been computed in part directly
!--    inside the WS advection routines are treated seperatly
!--    Momentum fluxes first:

       tn = 0
       !$OMP PARALLEL PRIVATE( i, j, k, tn, flag )
       !$ tn = omp_get_thread_num()
       IF ( .NOT. ws_scheme_mom .OR. sr /= 0  )  THEN
          !$OMP DO
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb, nzt
!
!--                Flag 23 is used to mask surface fluxes as well as model-top 
!--                fluxes, which are added further below. 
                   flag = MERGE( 1.0_wp, 0.0_wp,                               &
                                 BTEST( wall_flags_total_0(k,j,i), 23 ) ) *    &
                          MERGE( 1.0_wp, 0.0_wp,                               &
                                 BTEST( wall_flags_total_0(k,j,i), 9  ) )

                   ust = 0.5_wp * ( u(k,j,i)   - hom(k,1,1,sr) +               &
                                    u(k+1,j,i) - hom(k+1,1,1,sr) )
                   vst = 0.5_wp * ( v(k,j,i)   - hom(k,1,2,sr) +               &
                                    v(k+1,j,i) - hom(k+1,1,2,sr) )
!
!--                Momentum flux w*u*
                   sums_l(k,13,tn) = sums_l(k,13,tn) + 0.5_wp *                &
                                                     ( w(k,j,i-1) + w(k,j,i) ) &
                                           * rho_air_zw(k)                     &
                                           * momentumflux_output_conversion(k) &
                                                     * ust * rmask(j,i,sr)     &
                                                           * flag
!
!--                Momentum flux w*v*
                   sums_l(k,15,tn) = sums_l(k,15,tn) + 0.5_wp *                &
                                                     ( w(k,j-1,i) + w(k,j,i) ) &
                                           * rho_air_zw(k)                     &
                                           * momentumflux_output_conversion(k) &
                                                     * vst * rmask(j,i,sr)     &
                                                           * flag
                ENDDO
             ENDDO
          ENDDO

       ENDIF
       IF ( .NOT. ws_scheme_sca .OR. sr /= 0 )  THEN
          !$OMP DO
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb, nzt
                   flag = MERGE( 1.0_wp, 0.0_wp,                               &
                                 BTEST( wall_flags_total_0(k,j,i), 23 ) ) *    &
                          MERGE( 1.0_wp, 0.0_wp,                               &
                                 BTEST( wall_flags_total_0(k,j,i), 9  ) )
!
!--                Vertical heat flux
                   sums_l(k,17,tn) = sums_l(k,17,tn) + 0.5_wp *                &
                           ( pt(k,j,i)   - hom(k,1,4,sr) +                     &
                             pt(k+1,j,i) - hom(k+1,1,4,sr) )                   &
                           * rho_air_zw(k)                                     &
                           * heatflux_output_conversion(k)                     &
                           * w(k,j,i) * rmask(j,i,sr) * flag
                   IF ( humidity )  THEN
                      pts = 0.5_wp * ( q(k,j,i)   - hom(k,1,41,sr) +           &
                                      q(k+1,j,i) - hom(k+1,1,41,sr) )
                      sums_l(k,49,tn) = sums_l(k,49,tn) + pts * w(k,j,i) *     &
                                       rho_air_zw(k) *                         &
                                       waterflux_output_conversion(k) *        &
                                       rmask(j,i,sr) * flag
                   ENDIF
                   IF ( passive_scalar )  THEN
                      pts = 0.5_wp * ( s(k,j,i)   - hom(k,1,115,sr) +           &
                                      s(k+1,j,i) - hom(k+1,1,115,sr) )
                      sums_l(k,114,tn) = sums_l(k,114,tn) + pts * w(k,j,i) *    &
                                        rmask(j,i,sr) * flag
                   ENDIF
                ENDDO
             ENDDO
          ENDDO

       ENDIF

!
!--    Density at top follows Neumann condition
       IF ( ocean_mode )  THEN
          sums_l(nzt+1,64,tn) = sums_l(nzt,64,tn)
          sums_l(nzt+1,71,tn) = sums_l(nzt,71,tn)
       ENDIF

!
!--    Divergence of vertical flux of resolved scale energy and pressure
!--    fluctuations as well as flux of pressure fluctuation itself (68).
!--    First calculate the products, then the divergence.
!--    Calculation is time consuming. Do it only, if profiles shall be plotted.
       IF ( hom(nzb+1,2,55,0) /= 0.0_wp  .OR.  hom(nzb+1,2,68,0) /= 0.0_wp )   &
       THEN
          sums_ll = 0.0_wp  ! local array

          !$OMP DO
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )

                   sums_ll(k,1) = sums_ll(k,1) + 0.5_wp * w(k,j,i) * (         &
                  ( 0.25_wp * ( u(k,j,i)+u(k+1,j,i)+u(k,j,i+1)+u(k+1,j,i+1) )  &
                            - 0.5_wp * ( hom(k,1,1,sr) + hom(k+1,1,1,sr) ) )**2&
                + ( 0.25_wp * ( v(k,j,i)+v(k+1,j,i)+v(k,j+1,i)+v(k+1,j+1,i) )  &
                            - 0.5_wp * ( hom(k,1,2,sr) + hom(k+1,1,2,sr) ) )**2&
                + w(k,j,i)**2                                        ) * flag

                   sums_ll(k,2) = sums_ll(k,2) + 0.5_wp * w(k,j,i)             &
                                       * ( ( p(k,j,i) + p(k+1,j,i) )           &
                                         / momentumflux_output_conversion(k) ) &
                                       * flag

                ENDDO
             ENDDO
          ENDDO
          sums_ll(0,1)     = 0.0_wp    ! because w is zero at the bottom
          sums_ll(nzt+1,1) = 0.0_wp
          sums_ll(0,2)     = 0.0_wp
          sums_ll(nzt+1,2) = 0.0_wp

          DO  k = nzb+1, nzt
             sums_l(k,55,tn) = ( sums_ll(k,1) - sums_ll(k-1,1) ) * ddzw(k)
             sums_l(k,56,tn) = ( sums_ll(k,2) - sums_ll(k-1,2) ) * ddzw(k)
             sums_l(k,68,tn) = sums_ll(k,2)
          ENDDO
          sums_l(nzb,55,tn) = sums_l(nzb+1,55,tn)
          sums_l(nzb,56,tn) = sums_l(nzb+1,56,tn)
          sums_l(nzb,68,tn) = 0.0_wp    ! because w* = 0 at nzb

       ENDIF

!
!--    Divergence of vertical flux of SGS TKE and the flux itself (69)
       IF ( hom(nzb+1,2,57,0) /= 0.0_wp  .OR.  hom(nzb+1,2,69,0) /= 0.0_wp )   &
       THEN
          !$OMP DO
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt

                   flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )

                   sums_l(k,57,tn) = sums_l(k,57,tn) - 0.5_wp * (              &
                   (km(k,j,i)+km(k+1,j,i)) * (e(k+1,j,i)-e(k,j,i)) * ddzu(k+1) &
                 - (km(k-1,j,i)+km(k,j,i)) * (e(k,j,i)-e(k-1,j,i)) * ddzu(k)   &
                                                                ) * ddzw(k)    &
                                                                  * flag

                   sums_l(k,69,tn) = sums_l(k,69,tn) - 0.5_wp * (              &
                   (km(k,j,i)+km(k+1,j,i)) * (e(k+1,j,i)-e(k,j,i)) * ddzu(k+1) &
                                                                )  * flag

                ENDDO
             ENDDO
          ENDDO
          sums_l(nzb,57,tn) = sums_l(nzb+1,57,tn)
          sums_l(nzb,69,tn) = sums_l(nzb+1,69,tn)

       ENDIF

!
!--    Horizontal heat fluxes (subgrid, resolved, total). 
!--    Do it only, if profiles shall be plotted. 
       IF ( hom(nzb+1,2,58,0) /= 0.0_wp ) THEN

          !$OMP DO
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )
!
!--                Subgrid horizontal heat fluxes u"pt", v"pt"
                   sums_l(k,58,tn) = sums_l(k,58,tn) - 0.5_wp *                &
                                                   ( kh(k,j,i) + kh(k,j,i-1) ) &
                                                 * ( pt(k,j,i-1) - pt(k,j,i) ) &
                                               * rho_air_zw(k)                 &
                                               * heatflux_output_conversion(k) &
                                                 * ddx * rmask(j,i,sr) * flag
                   sums_l(k,61,tn) = sums_l(k,61,tn) - 0.5_wp *                &
                                                   ( kh(k,j,i) + kh(k,j-1,i) ) &
                                                 * ( pt(k,j-1,i) - pt(k,j,i) ) &
                                               * rho_air_zw(k)                 &
                                               * heatflux_output_conversion(k) &
                                                 * ddy * rmask(j,i,sr) * flag
!
!--                Resolved horizontal heat fluxes u*pt*, v*pt*
                   sums_l(k,59,tn) = sums_l(k,59,tn) +                         &
                                                  ( u(k,j,i) - hom(k,1,1,sr) ) &
                                    * 0.5_wp * ( pt(k,j,i-1) - hom(k,1,4,sr) + &
                                                 pt(k,j,i)   - hom(k,1,4,sr) ) &
                                               * heatflux_output_conversion(k) &
                                               * flag
                   pts = 0.5_wp * ( pt(k,j-1,i) - hom(k,1,4,sr) +              &
                                    pt(k,j,i)   - hom(k,1,4,sr) )
                   sums_l(k,62,tn) = sums_l(k,62,tn) +                         &
                                                  ( v(k,j,i) - hom(k,1,2,sr) ) &
                                    * 0.5_wp * ( pt(k,j-1,i) - hom(k,1,4,sr) + &
                                                 pt(k,j,i)   - hom(k,1,4,sr) ) &
                                               * heatflux_output_conversion(k) &
                                               * flag
                ENDDO
             ENDDO
          ENDDO
!
!--       Fluxes at the surface must be zero (e.g. due to the Prandtl-layer)
          sums_l(nzb,58,tn) = 0.0_wp
          sums_l(nzb,59,tn) = 0.0_wp
          sums_l(nzb,60,tn) = 0.0_wp
          sums_l(nzb,61,tn) = 0.0_wp
          sums_l(nzb,62,tn) = 0.0_wp
          sums_l(nzb,63,tn) = 0.0_wp

       ENDIF
       !$OMP END PARALLEL

!
!--    Collect current large scale advection and subsidence tendencies for
!--    data output
       IF ( large_scale_forcing  .AND.  ( simulated_time > 0.0_wp ) )  THEN
!
!--       Interpolation in time of LSF_DATA 
          nt = 1
          DO WHILE ( simulated_time - dt_3d > time_vert(nt) )
             nt = nt + 1
          ENDDO
          IF ( simulated_time - dt_3d /= time_vert(nt) )  THEN
            nt = nt - 1
          ENDIF

          fac = ( simulated_time - dt_3d - time_vert(nt) )                     &
                / ( time_vert(nt+1)-time_vert(nt) )


          DO  k = nzb, nzt
             sums_ls_l(k,0) = td_lsa_lpt(k,nt)                                 &
                              + fac * ( td_lsa_lpt(k,nt+1) - td_lsa_lpt(k,nt) )
             sums_ls_l(k,1) = td_lsa_q(k,nt)                                   &
                              + fac * ( td_lsa_q(k,nt+1) - td_lsa_q(k,nt) )
          ENDDO

          sums_ls_l(nzt+1,0) = sums_ls_l(nzt,0)
          sums_ls_l(nzt+1,1) = sums_ls_l(nzt,1)

          IF ( large_scale_subsidence .AND. use_subsidence_tendencies )  THEN

             DO  k = nzb, nzt
                sums_ls_l(k,2) = td_sub_lpt(k,nt) + fac *                      &
                                 ( td_sub_lpt(k,nt+1) - td_sub_lpt(k,nt) )
                sums_ls_l(k,3) = td_sub_q(k,nt) + fac *                        &
                                 ( td_sub_q(k,nt+1) - td_sub_q(k,nt) )
             ENDDO

             sums_ls_l(nzt+1,2) = sums_ls_l(nzt,2)
             sums_ls_l(nzt+1,3) = sums_ls_l(nzt,3)

          ENDIF

       ENDIF

       tn = 0
       !$OMP PARALLEL PRIVATE( i, j, k, tn )
       !$ tn = omp_get_thread_num()       
       IF ( radiation .AND. radiation_scheme == 'rrtmg' )  THEN
          !$OMP DO
          DO  i = nxl, nxr
             DO  j =  nys, nyn
                DO  k = nzb+1, nzt+1
                   flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )

                   sums_l(k,100,tn)  = sums_l(k,100,tn)  + rad_lw_in(k,j,i)    &
                                       * rmask(j,i,sr) * flag
                   sums_l(k,101,tn)  = sums_l(k,101,tn)  + rad_lw_out(k,j,i)   &
                                       * rmask(j,i,sr) * flag
                   sums_l(k,102,tn)  = sums_l(k,102,tn)  + rad_sw_in(k,j,i)    &
                                       * rmask(j,i,sr) * flag
                   sums_l(k,103,tn)  = sums_l(k,103,tn)  + rad_sw_out(k,j,i)   &
                                       * rmask(j,i,sr) * flag
                   sums_l(k,104,tn)  = sums_l(k,104,tn)  + rad_lw_cs_hr(k,j,i) &
                                       * rmask(j,i,sr) * flag
                   sums_l(k,105,tn)  = sums_l(k,105,tn)  + rad_lw_hr(k,j,i)    &
                                       * rmask(j,i,sr) * flag
                   sums_l(k,106,tn)  = sums_l(k,106,tn)  + rad_sw_cs_hr(k,j,i) &
                                       * rmask(j,i,sr) * flag
                   sums_l(k,107,tn)  = sums_l(k,107,tn)  + rad_sw_hr(k,j,i)    &
                                       * rmask(j,i,sr) * flag
                ENDDO
             ENDDO
          ENDDO
       ENDIF

!
!--    Calculate the profiles for all other modules
       CALL module_interface_statistics( 'profiles', sr, tn, dots_max )
       !$OMP END PARALLEL

!
!--    Summation of thread sums
       IF ( threads_per_task > 1 )  THEN
          DO  i = 1, threads_per_task-1
             sums_l(:,3,0)          = sums_l(:,3,0) + sums_l(:,3,i)
             sums_l(:,4:40,0)       = sums_l(:,4:40,0) + sums_l(:,4:40,i)
             sums_l(:,45:pr_palm,0) = sums_l(:,45:pr_palm,0) + &
                                      sums_l(:,45:pr_palm,i)
             IF ( max_pr_user > 0 )  THEN
                sums_l(:,pr_palm+1:pr_palm+max_pr_user,0) = &
                                   sums_l(:,pr_palm+1:pr_palm+max_pr_user,0) + &
                                   sums_l(:,pr_palm+1:pr_palm+max_pr_user,i)
             ENDIF

             IF ( air_chemistry )  THEN
                IF ( max_pr_cs > 0 )  THEN                                 
                     sums_l(:,pr_palm+max_pr_user+1:pr_palm + max_pr_user+ max_pr_cs,0) =          &
                               sums_l(:,pr_palm+max_pr_user+1:pr_palm + max_pr_user+max_pr_cs,0) + &
                               sums_l(:,pr_palm+max_pr_user+1:pr_palm + max_pr_user+max_pr_cs,i)

                ENDIF
             ENDIF
             IF ( salsa )  THEN
                IF ( max_pr_cs > 0 )  THEN
                   sums_l(:,pr_palm+max_pr_user+max_pr_cs+1:pr_palm+max_pr_user+max_pr_cs+max_pr_salsa,0) =    &
                      sums_l(:,pr_palm+max_pr_user+max_pr_cs+1:pr_palm+max_pr_user+max_pr_cs+max_pr_salsa,0) + &
                      sums_l(:,pr_palm+max_pr_user+max_pr_cs+1:pr_palm+max_pr_user+max_pr_cs+max_pr_salsa,i)

                ENDIF
             ENDIF
          ENDDO
       ENDIF

#if defined( __parallel )

!
!--    Compute total sum from local sums
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( sums_l(nzb,1,0), sums(nzb,1), ngp_sums, MPI_REAL,   &
                           MPI_SUM, comm2d, ierr )
       IF ( large_scale_forcing )  THEN
          CALL MPI_ALLREDUCE( sums_ls_l(nzb,2), sums(nzb,83), ngp_sums_ls,     &
                              MPI_REAL, MPI_SUM, comm2d, ierr )
       ENDIF

       IF ( air_chemistry  .AND.  max_pr_cs > 0 )  THEN
          IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
             DO  i = 1, max_pr_cs
                CALL MPI_ALLREDUCE( sums_l(nzb,pr_palm+max_pr_user+i,0),       &
                                    sums(nzb,pr_palm+max_pr_user+i),           &
                                    nzt+2-nzb, MPI_REAL, MPI_SUM, comm2d, ierr )
             ENDDO
       ENDIF

       IF ( salsa  .AND.  max_pr_salsa > 0 )  THEN
          IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
             DO  i = 1, max_pr_salsa
                CALL MPI_ALLREDUCE( sums_l(nzb,pr_palm+max_pr_user+max_pr_cs+i,0),                 &
                                    sums(nzb,pr_palm+max_pr_user+max_pr_user+i),                   &
                                    nzt+2-nzb, MPI_REAL, MPI_SUM, comm2d, ierr )
             ENDDO
       ENDIF

#else
       sums = sums_l(:,:,0)
       IF ( large_scale_forcing )  THEN
          sums(:,81:88) = sums_ls_l
       ENDIF
#endif

!
!--    Final values are obtained by division by the total number of grid points 
!--    used for summation. After that store profiles.
!--    Check, if statistical regions do contain at least one grid point at the
!--    respective k-level, otherwise division by zero will lead to undefined
!--    values, which may cause e.g. problems with NetCDF output
!--    Profiles:
       DO  k = nzb, nzt+1
          sums(k,3)             = sums(k,3)             / ngp_2dh(sr)
          sums(k,12:22)         = sums(k,12:22)         / ngp_2dh(sr)
          sums(k,30:32)         = sums(k,30:32)         / ngp_2dh(sr)
          sums(k,35:39)         = sums(k,35:39)         / ngp_2dh(sr)
          sums(k,45:53)         = sums(k,45:53)         / ngp_2dh(sr)
          sums(k,55:63)         = sums(k,55:63)         / ngp_2dh(sr)
          sums(k,81:88)         = sums(k,81:88)         / ngp_2dh(sr)
          sums(k,89:112)        = sums(k,89:112)        / ngp_2dh(sr)
          sums(k,114)           = sums(k,114)           / ngp_2dh(sr)
          sums(k,117)           = sums(k,117)           / ngp_2dh(sr)
          IF ( ngp_2dh_s_inner(k,sr) /= 0 )  THEN
             sums(k,8:11)          = sums(k,8:11)          / ngp_2dh_s_inner(k,sr)
             sums(k,23:29)         = sums(k,23:29)         / ngp_2dh_s_inner(k,sr)
             sums(k,33:34)         = sums(k,33:34)         / ngp_2dh_s_inner(k,sr)
             sums(k,40)            = sums(k,40)            / ngp_2dh_s_inner(k,sr)
             sums(k,54)            = sums(k,54)            / ngp_2dh_s_inner(k,sr)
             sums(k,64)            = sums(k,64)            / ngp_2dh_s_inner(k,sr)
             sums(k,70:80)         = sums(k,70:80)         / ngp_2dh_s_inner(k,sr)
             sums(k,116)           = sums(k,116)           / ngp_2dh_s_inner(k,sr)
             sums(k,118:pr_palm-2) = sums(k,118:pr_palm-2) / ngp_2dh_s_inner(k,sr)
             sums(k,123)           = sums(k,123) * ngp_2dh_s_inner(k,sr)  / ngp_2dh(sr)
          ENDIF
       ENDDO

!--    u* and so on
!--    As sums(nzb:nzb+3,pr_palm) are full 2D arrays (us, usws, vsws, ts) whose
!--    size is always ( nx + 1 ) * ( ny + 1 ), defined at the first grid layer
!--    above the topography, they are being divided by ngp_2dh(sr)
       sums(nzb:nzb+3,pr_palm)    = sums(nzb:nzb+3,pr_palm)    / &
                                    ngp_2dh(sr)
       sums(nzb+12,pr_palm)       = sums(nzb+12,pr_palm)       / &    ! qs
                                    ngp_2dh(sr)
       sums(nzb+13,pr_palm)       = sums(nzb+13,pr_palm)       / &    ! ss
                                    ngp_2dh(sr)
       sums(nzb+14,pr_palm)       = sums(nzb+14,pr_palm)       / &    ! surface temperature
                                    ngp_2dh(sr)
!--    eges, e*
       sums(nzb+4:nzb+5,pr_palm)  = sums(nzb+4:nzb+5,pr_palm)  / &
                                    ngp_3d(sr)
!--    Old and new divergence
       sums(nzb+9:nzb+10,pr_palm) = sums(nzb+9:nzb+10,pr_palm) / &
                                    ngp_3d_inner(sr)

!--    User-defined profiles
       IF ( max_pr_user > 0 )  THEN
          DO  k = nzb, nzt+1
             sums(k,pr_palm+1:pr_palm+max_pr_user) = &
                                    sums(k,pr_palm+1:pr_palm+max_pr_user) / &
                                    ngp_2dh_s_inner(k,sr)
          ENDDO
       ENDIF

       IF ( air_chemistry ) THEN
          IF ( max_pr_cs > 0 )  THEN                  
             DO k = nzb, nzt+1
                sums(k, pr_palm+1:pr_palm+max_pr_user+max_pr_cs) = &
                                 sums(k, pr_palm+1:pr_palm+max_pr_user+max_pr_cs) / &
                                 ngp_2dh_s_inner(k,sr)
             ENDDO
          ENDIF  
       ENDIF

       IF ( salsa ) THEN
          IF ( max_pr_salsa > 0 )  THEN
             DO k = nzb, nzt+1
                sums(k,pr_palm+max_pr_user+max_pr_cs+1:pr_palm+max_pr_user+max_pr_cs+max_pr_salsa) = &
                  sums(k,pr_palm+max_pr_user+max_pr_cs+1:pr_palm+max_pr_user+max_pr_cs+max_pr_salsa) &
                  / ngp_2dh_s_inner(k,sr)
             ENDDO
          ENDIF  
       ENDIF

!
!--    Collect horizontal average in hom.
!--    Compute deduced averages (e.g. total heat flux)
       hom(:,1,3,sr)  = sums(:,3)      ! w
       hom(:,1,8,sr)  = sums(:,8)      ! e     profiles 5-7 are initial profiles
       hom(:,1,9,sr)  = sums(:,9)      ! km
       hom(:,1,10,sr) = sums(:,10)     ! kh
       hom(:,1,11,sr) = sums(:,11)     ! l
       hom(:,1,12,sr) = sums(:,12)     ! w"u"
       hom(:,1,13,sr) = sums(:,13)     ! w*u*
       hom(:,1,14,sr) = sums(:,14)     ! w"v"
       hom(:,1,15,sr) = sums(:,15)     ! w*v*
       hom(:,1,16,sr) = sums(:,16)     ! w"pt"
       hom(:,1,17,sr) = sums(:,17)     ! w*pt*
       hom(:,1,18,sr) = sums(:,16) + sums(:,17)    ! wpt
       hom(:,1,19,sr) = sums(:,12) + sums(:,13)    ! wu
       hom(:,1,20,sr) = sums(:,14) + sums(:,15)    ! wv
       hom(:,1,21,sr) = sums(:,21)     ! w*pt*BC
       hom(:,1,22,sr) = sums(:,16) + sums(:,21)    ! wptBC
                                       ! profile 24 is initial profile (sa)
                                       ! profiles 25-29 left empty for initial 
                                       ! profiles
       hom(:,1,30,sr) = sums(:,30)     ! u*2
       hom(:,1,31,sr) = sums(:,31)     ! v*2
       hom(:,1,32,sr) = sums(:,32)     ! w*2
       hom(:,1,33,sr) = sums(:,33)     ! pt*2
       hom(:,1,34,sr) = sums(:,34)     ! e*
       hom(:,1,35,sr) = sums(:,35)     ! w*2pt*
       hom(:,1,36,sr) = sums(:,36)     ! w*pt*2
       hom(:,1,37,sr) = sums(:,37)     ! w*e*
       hom(:,1,38,sr) = sums(:,38)     ! w*3
       hom(:,1,39,sr) = sums(:,38) / ( abs( sums(:,32) ) + 1E-20_wp )**1.5_wp   ! Sw
       hom(:,1,40,sr) = sums(:,40)     ! p
       hom(:,1,45,sr) = sums(:,45)     ! w"vpt"
       hom(:,1,46,sr) = sums(:,46)     ! w*vpt*       
       hom(:,1,47,sr) = sums(:,45) + sums(:,46)    ! wvpt
       hom(:,1,48,sr) = sums(:,48)     ! w"q" (w"qv")
       hom(:,1,49,sr) = sums(:,49)     ! w*q* (w*qv*)
       hom(:,1,50,sr) = sums(:,48) + sums(:,49)    ! wq (wqv)
       hom(:,1,51,sr) = sums(:,51)     ! w"qv"
       hom(:,1,52,sr) = sums(:,52)     ! w*qv*       
       hom(:,1,53,sr) = sums(:,52) + sums(:,51)    ! wq (wqv)
       hom(:,1,54,sr) = sums(:,54)     ! ql
       hom(:,1,55,sr) = sums(:,55)     ! w*u*u*/dz
       hom(:,1,56,sr) = sums(:,56)     ! w*p*/dz
       hom(:,1,57,sr) = sums(:,57)     ! ( w"e + w"p"/rho_ocean )/dz
       hom(:,1,58,sr) = sums(:,58)     ! u"pt"
       hom(:,1,59,sr) = sums(:,59)     ! u*pt*
       hom(:,1,60,sr) = sums(:,58) + sums(:,59)    ! upt_t
       hom(:,1,61,sr) = sums(:,61)     ! v"pt"
       hom(:,1,62,sr) = sums(:,62)     ! v*pt*
       hom(:,1,63,sr) = sums(:,61) + sums(:,62)    ! vpt_t
       hom(:,1,64,sr) = sums(:,64)     ! rho_ocean
       hom(:,1,65,sr) = sums(:,65)     ! w"sa"
       hom(:,1,66,sr) = sums(:,66)     ! w*sa*
       hom(:,1,67,sr) = sums(:,65) + sums(:,66)    ! wsa
       hom(:,1,68,sr) = sums(:,68)     ! w*p*
       hom(:,1,69,sr) = sums(:,69)     ! w"e + w"p"/rho_ocean
       hom(:,1,70,sr) = sums(:,70)     ! q*2
       hom(:,1,71,sr) = sums(:,71)     ! prho
       hom(:,1,72,sr) = hyp * 1E-2_wp  ! hyp in hPa
       hom(:,1,123,sr) = sums(:,123)   ! nc
       hom(:,1,73,sr) = sums(:,73)     ! nr
       hom(:,1,74,sr) = sums(:,74)     ! qr
       hom(:,1,75,sr) = sums(:,75)     ! qc
       hom(:,1,76,sr) = sums(:,76)     ! prr (precipitation rate)
                                       ! 77 is initial density profile
       hom(:,1,78,sr) = ug             ! ug
       hom(:,1,79,sr) = vg             ! vg
       hom(:,1,80,sr) = w_subs         ! w_subs

       IF ( large_scale_forcing )  THEN
          hom(:,1,81,sr) = sums_ls_l(:,0)          ! td_lsa_lpt
          hom(:,1,82,sr) = sums_ls_l(:,1)          ! td_lsa_q
          IF ( use_subsidence_tendencies )  THEN
             hom(:,1,83,sr) = sums_ls_l(:,2)       ! td_sub_lpt
             hom(:,1,84,sr) = sums_ls_l(:,3)       ! td_sub_q
          ELSE
             hom(:,1,83,sr) = sums(:,83)           ! td_sub_lpt
             hom(:,1,84,sr) = sums(:,84)           ! td_sub_q
          ENDIF
          hom(:,1,85,sr) = sums(:,85)              ! td_nud_lpt
          hom(:,1,86,sr) = sums(:,86)              ! td_nud_q
          hom(:,1,87,sr) = sums(:,87)              ! td_nud_u
          hom(:,1,88,sr) = sums(:,88)              ! td_nud_v
       ENDIF

       IF ( land_surface )  THEN
          hom(:,1,89,sr) = sums(:,89)              ! t_soil
                                                   ! 90 is initial t_soil profile
          hom(:,1,91,sr) = sums(:,91)              ! m_soil
                                                   ! 92 is initial m_soil profile
          hom(:,1,93,sr)  = sums(:,93)             ! ghf
          hom(:,1,94,sr)  = sums(:,94)             ! qsws_liq
          hom(:,1,95,sr)  = sums(:,95)             ! qsws_soil
          hom(:,1,96,sr)  = sums(:,96)             ! qsws_veg
          hom(:,1,97,sr)  = sums(:,97)             ! r_a
          hom(:,1,98,sr) = sums(:,98)              ! r_s

       ENDIF

       IF ( radiation )  THEN
          hom(:,1,99,sr) = sums(:,99)            ! rad_net
          hom(:,1,100,sr) = sums(:,100)            ! rad_lw_in
          hom(:,1,101,sr) = sums(:,101)            ! rad_lw_out
          hom(:,1,102,sr) = sums(:,102)            ! rad_sw_in
          hom(:,1,103,sr) = sums(:,103)            ! rad_sw_out

          IF ( radiation_scheme == 'rrtmg' )  THEN
             hom(:,1,104,sr) = sums(:,104)            ! rad_lw_cs_hr
             hom(:,1,105,sr) = sums(:,105)            ! rad_lw_hr
             hom(:,1,106,sr) = sums(:,106)            ! rad_sw_cs_hr
             hom(:,1,107,sr) = sums(:,107)            ! rad_sw_hr

             hom(:,1,108,sr) = sums(:,108)            ! rrtm_aldif
             hom(:,1,109,sr) = sums(:,109)            ! rrtm_aldir
             hom(:,1,110,sr) = sums(:,110)            ! rrtm_asdif
             hom(:,1,111,sr) = sums(:,111)            ! rrtm_asdir
          ENDIF
       ENDIF

       hom(:,1,112,sr) = sums(:,112)            !: L

       IF ( passive_scalar )  THEN
          hom(:,1,117,sr) = sums(:,117)     ! w"s"
          hom(:,1,114,sr) = sums(:,114)     ! w*s*
          hom(:,1,118,sr) = sums(:,117) + sums(:,114)    ! ws 
          hom(:,1,116,sr) = sums(:,116)     ! s*2
       ENDIF

       hom(:,1,119,sr) = rho_air       ! rho_air in Kg/m^3
       hom(:,1,120,sr) = rho_air_zw    ! rho_air_zw in Kg/m^3

       hom(:,1,pr_palm,sr) =   sums(:,pr_palm)
                                       ! u*, w'u', w'v', t* (in last profile)

       IF ( max_pr_user > 0 )  THEN    ! user-defined profiles
          hom(:,1,pr_palm+1:pr_palm+max_pr_user,sr) = &
                               sums(:,pr_palm+1:pr_palm+max_pr_user)
       ENDIF

       IF ( air_chemistry )  THEN
          IF ( max_pr_cs > 0 )  THEN    ! chem_spcs profiles      
             hom(:, 1, pr_palm+max_pr_user+1:pr_palm + max_pr_user+max_pr_cs, sr) = &
                               sums(:, pr_palm+max_pr_user+1:pr_palm+max_pr_user+max_pr_cs)
          ENDIF
       ENDIF

       IF ( salsa )  THEN
          IF ( max_pr_salsa > 0 )  THEN    ! salsa profiles
             hom(:,1,pr_palm+max_pr_user+max_pr_cs+1:pr_palm+max_pr_user+max_pr_cs+max_pr_salsa, sr) = &
                  sums(:,pr_palm+max_pr_user+max_pr_cs+1:pr_palm+max_pr_user+max_pr_cs+max_pr_salsa)
          ENDIF
       ENDIF
!
!--    Determine the boundary layer height using two different schemes.
!--    First scheme: Starting from the Earth's (Ocean's) surface, look for the
!--    first relative minimum (maximum) of the total heat flux.
!--    The corresponding height is assumed as the boundary layer height, if it
!--    is less than 1.5 times the height where the heat flux becomes negative
!--    (positive) for the first time. Attention: the resolved vertical sensible 
!--    heat flux (hom(:,1,17,sr) = w*pt*) is not known at the beginning because
!--    the calculation happens in advec_s_ws which is called after 
!--    flow_statistics. Therefore z_i is directly taken from restart data at 
!--    the beginning of restart runs.  
       IF ( TRIM( initializing_actions ) /= 'read_restart_data' .OR.           &
            simulated_time_at_begin /= simulated_time ) THEN

          z_i(1) = 0.0_wp
          first = .TRUE.

          IF ( ocean_mode )  THEN
             DO  k = nzt, nzb+1, -1
                IF ( first  .AND.  hom(k,1,18,sr) < -1.0E-8_wp )  THEN
                   first = .FALSE.
                   height = zw(k)
                ENDIF
                IF ( hom(k,1,18,sr) < -1.0E-8_wp  .AND.                        &
                     hom(k-1,1,18,sr) > hom(k,1,18,sr) )  THEN
                   IF ( zw(k) < 1.5_wp * height )  THEN
                      z_i(1) = zw(k)
                   ELSE
                      z_i(1) = height
                   ENDIF
                   EXIT
                ENDIF
             ENDDO
          ELSE
             DO  k = nzb, nzt-1
                IF ( first  .AND.  hom(k,1,18,sr) < -1.0E-8_wp )  THEN
                   first = .FALSE.
                   height = zw(k)
                ENDIF
                IF ( hom(k,1,18,sr) < -1.0E-8_wp  .AND.                        &
                     hom(k+1,1,18,sr) > hom(k,1,18,sr) )  THEN
                   IF ( zw(k) < 1.5_wp * height )  THEN
                      z_i(1) = zw(k)
                   ELSE
                      z_i(1) = height
                   ENDIF
                   EXIT
                ENDIF
             ENDDO
          ENDIF

!
!--       Second scheme: Gradient scheme from Sullivan et al. (1998), modified 
!--       by Uhlenbrock(2006). The boundary layer height is the height with the 
!--       maximal local temperature gradient: starting from the second (the
!--       last but one) vertical gridpoint, the local gradient must be at least 
!--       0.2K/100m and greater than the next four gradients.
!--       WARNING: The threshold value of 0.2K/100m must be adjusted for the
!--       ocean case! 
          z_i(2) = 0.0_wp
          DO  k = nzb+1, nzt+1
             dptdz(k) = ( hom(k,1,4,sr) - hom(k-1,1,4,sr) ) * ddzu(k)
          ENDDO
          dptdz_threshold = 0.2_wp / 100.0_wp

          IF ( ocean_mode )  THEN
             DO  k = nzt+1, nzb+5, -1
                IF ( dptdz(k) > dptdz_threshold  .AND.                         &
                     dptdz(k) > dptdz(k-1)  .AND.  dptdz(k) > dptdz(k-2)  .AND.&
                     dptdz(k) > dptdz(k-3)  .AND.  dptdz(k) > dptdz(k-4) )  THEN
                   z_i(2) = zw(k-1)
                   EXIT
                ENDIF
             ENDDO
          ELSE
             DO  k = nzb+1, nzt-3
                IF ( dptdz(k) > dptdz_threshold  .AND.                         &
                     dptdz(k) > dptdz(k+1)  .AND.  dptdz(k) > dptdz(k+2)  .AND.&
                     dptdz(k) > dptdz(k+3)  .AND.  dptdz(k) > dptdz(k+4) )  THEN
                   z_i(2) = zw(k-1)
                   EXIT
                ENDIF
             ENDDO
          ENDIF

       ENDIF

       hom(nzb+6,1,pr_palm,sr) = z_i(1)
       hom(nzb+7,1,pr_palm,sr) = z_i(2)

!
!--    Determine vertical index which is nearest to the mean surface level
!--    height of the respective statistic region
       DO  k = nzb, nzt
          IF ( zw(k) >= mean_surface_level_height(sr) )  THEN
             k_surface_level = k
             EXIT
          ENDIF
       ENDDO

!
!--    Computation of both the characteristic vertical velocity and
!--    the characteristic convective boundary layer temperature.
!--    The inversion height entering into the equation is defined with respect
!--    to the mean surface level height of the respective statistic region.
!--    The horizontal average at surface level index + 1 is input for the
!--    average temperature.
       IF ( hom(k_surface_level,1,18,sr) > 1.0E-8_wp  .AND.  z_i(1) /= 0.0_wp )&
       THEN
          hom(nzb+8,1,pr_palm,sr) =                                            &
             ( g / hom(k_surface_level+1,1,4,sr) *                             &
             ( hom(k_surface_level,1,18,sr) /                                  &
             ( heatflux_output_conversion(nzb) * rho_air(nzb) ) )              &
             * ABS( z_i(1) - mean_surface_level_height(sr) ) )**0.333333333_wp
       ELSE
          hom(nzb+8,1,pr_palm,sr)  = 0.0_wp
       ENDIF

!
!--    Collect the time series quantities. Please note, timeseries quantities 
!--    which are collected from horizontally averaged profiles, e.g. wpt 
!--    or pt(zp), are treated specially. In case of elevated model surfaces, 
!--    index nzb+1 might be within topography and data will be zero. Therefore,
!--    take value for the first atmosphere index, which is topo_min_level+1. 
       ts_value(1,sr) = hom(nzb+4,1,pr_palm,sr)        ! E
       ts_value(2,sr) = hom(nzb+5,1,pr_palm,sr)        ! E*
       ts_value(3,sr) = dt_3d
       ts_value(4,sr) = hom(nzb,1,pr_palm,sr)          ! u*
       ts_value(5,sr) = hom(nzb+3,1,pr_palm,sr)        ! th*
       ts_value(6,sr) = u_max
       ts_value(7,sr) = v_max
       ts_value(8,sr) = w_max
       ts_value(9,sr) = hom(nzb+10,1,pr_palm,sr)       ! new divergence
       ts_value(10,sr) = hom(nzb+9,1,pr_palm,sr)       ! old Divergence
       ts_value(11,sr) = hom(nzb+6,1,pr_palm,sr)       ! z_i(1)
       ts_value(12,sr) = hom(nzb+7,1,pr_palm,sr)       ! z_i(2)
       ts_value(13,sr) = hom(nzb+8,1,pr_palm,sr)       ! w*
       ts_value(14,sr) = hom(nzb,1,16,sr)              ! w'pt'   at k=0
       ts_value(15,sr) = hom(topo_min_level+1,1,16,sr) ! w'pt'   at k=1
       ts_value(16,sr) = hom(topo_min_level+1,1,18,sr) ! wpt     at k=1
       ts_value(17,sr) = hom(nzb+14,1,pr_palm,sr)      ! pt(0)
       ts_value(18,sr) = hom(topo_min_level+1,1,4,sr)  ! pt(zp)
       ts_value(19,sr) = hom(nzb+1,1,pr_palm,sr)       ! u'w'    at k=0
       ts_value(20,sr) = hom(nzb+2,1,pr_palm,sr)       ! v'w'    at k=0
       ts_value(21,sr) = hom(nzb,1,48,sr)              ! w"q"    at k=0

       IF ( .NOT. neutral )  THEN
          ts_value(22,sr) = hom(nzb,1,112,sr)          ! L
       ELSE
          ts_value(22,sr) = 1.0E10_wp
       ENDIF

       ts_value(23,sr) = hom(nzb+12,1,pr_palm,sr)   ! q*

       IF ( passive_scalar )  THEN
          ts_value(24,sr) = hom(nzb+13,1,117,sr)       ! w"s" ( to do ! )
          ts_value(25,sr) = hom(nzb+13,1,pr_palm,sr)   ! s*
       ENDIF

!
!--    Collect land surface model timeseries
       IF ( land_surface )  THEN
          ts_value(dots_soil  ,sr) = hom(nzb,1,93,sr)           ! ghf
          ts_value(dots_soil+1,sr) = hom(nzb,1,94,sr)           ! qsws_liq
          ts_value(dots_soil+2,sr) = hom(nzb,1,95,sr)           ! qsws_soil
          ts_value(dots_soil+3,sr) = hom(nzb,1,96,sr)           ! qsws_veg
          ts_value(dots_soil+4,sr) = hom(nzb,1,97,sr)           ! r_a
          ts_value(dots_soil+5,sr) = hom(nzb,1,98,sr)           ! r_s
       ENDIF
!
!--    Collect radiation model timeseries
       IF ( radiation )  THEN
          ts_value(dots_rad,sr)   = hom(nzb,1,99,sr)           ! rad_net
          ts_value(dots_rad+1,sr) = hom(nzb,1,100,sr)          ! rad_lw_in
          ts_value(dots_rad+2,sr) = hom(nzb,1,101,sr)          ! rad_lw_out
          ts_value(dots_rad+3,sr) = hom(nzb,1,102,sr)          ! rad_sw_in
          ts_value(dots_rad+4,sr) = hom(nzb,1,103,sr)          ! rad_sw_out

          IF ( radiation_scheme == 'rrtmg' )  THEN
             ts_value(dots_rad+5,sr) = hom(nzb,1,108,sr)          ! rrtm_aldif
             ts_value(dots_rad+6,sr) = hom(nzb,1,109,sr)          ! rrtm_aldir
             ts_value(dots_rad+7,sr) = hom(nzb,1,110,sr)          ! rrtm_asdif
             ts_value(dots_rad+8,sr) = hom(nzb,1,111,sr)          ! rrtm_asdir
          ENDIF

       ENDIF

!
!--    Calculate additional statistics provided by other modules
       CALL module_interface_statistics( 'time_series', sr, 0, dots_max )

    ENDDO    ! loop of the subregions

!
!-- If required, sum up horizontal averages for subsequent time averaging.
!-- Do not sum, if flow statistics is called before the first initial time step.
    IF ( do_sum  .AND.  simulated_time /= 0.0_wp )  THEN
       IF ( average_count_pr == 0 )  hom_sum = 0.0_wp
       hom_sum = hom_sum + hom(:,1,:,:)
       average_count_pr = average_count_pr + 1
       do_sum = .FALSE.
    ENDIF

!
!-- Set flag for other UPs (e.g. output routines, but also buoyancy).
!-- This flag is reset after each time step in time_integration.
    flow_statistics_called = .TRUE.

    CALL cpu_log( log_point(10), 'flow_statistics', 'stop' )


 END SUBROUTINE flow_statistics
