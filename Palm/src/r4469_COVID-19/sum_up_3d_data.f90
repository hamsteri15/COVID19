!> @file sum_up_3d_data.f90
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
! $Id: sum_up_3d_data.f90 4442 2020-03-04 19:21:13Z suehring $
! Change order of dimension in surface array %frac to allow for better 
! vectorization.
! 
! 4441 2020-03-04 19:20:35Z suehring
! Move 2-m potential temperature output to diagnostic_output_quantities
! 
! 4182 2019-08-22 15:20:23Z scharf
! Corrected "Former revisions" section
! 
! 4048 2019-06-21 21:00:21Z knoop
! Moved tcm_3d_data_averaging to module_interface
! 
! 4039 2019-06-18 10:32:41Z suehring
! Modularize diagnostic output
! 
! 3994 2019-05-22 18:08:09Z suehring
! output of turbulence intensity added
! 
! 3943 2019-05-02 09:50:41Z maronga
! Added output of qsws_av for green roofs.
! 
! 3933 2019-04-25 12:33:20Z kanani
! Formatting
! 
! 3773 2019-03-01 08:56:57Z maronga
! Added output of theta_2m*_xy_av
! 
! 3761 2019-02-25 15:31:42Z raasch
! unused variables removed
! 
! 3655 2019-01-07 16:51:22Z knoop
! Implementation of the PALM module interface
!
! Revision 1.1  2006/02/23 12:55:23  raasch
! Initial revision
!
!
! Description:
! ------------
!> Sum-up the values of 3d-arrays. The real averaging is later done in routine
!> average_3d_data. 
!------------------------------------------------------------------------------!
 SUBROUTINE sum_up_3d_data
 

    USE arrays_3d,                                                             &
        ONLY:  dzw, d_exner, e, heatflux_output_conversion, p,                 &
               pt, q, ql, ql_c, ql_v, s, u, v, vpt, w,                         &
               waterflux_output_conversion

    USE averaging,                                                             &
        ONLY:  e_av, ghf_av, lpt_av, lwp_av, ol_av, p_av, pc_av, pr_av, pt_av, &
               q_av, ql_av, ql_c_av, ql_v_av, ql_vp_av, qsws_av,               &
               qv_av, r_a_av, s_av, shf_av, ssws_av, ts_av, tsurf_av, u_av,    &
               us_av, v_av, vpt_av, w_av, z0_av, z0h_av, z0q_av

    USE basic_constants_and_equations_mod,                                     &
        ONLY:  c_p, lv_d_cp, l_v

    USE bulk_cloud_model_mod,                                                  &
        ONLY:  bulk_cloud_model

    USE control_parameters,                                                    &
        ONLY:  average_count_3d, doav, doav_n, rho_surface, urban_surface,     &
               varnamelength

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point

    USE indices,                                                               &
        ONLY:  nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg, nzb, nzt 

    USE kinds

    USE module_interface,                                                      &
        ONLY:  module_interface_3d_data_averaging

    USE particle_attributes,                                                   &
        ONLY:  grid_particles, number_of_particles, particles, prt_count

    USE surface_mod,                                                           &
        ONLY:  ind_pav_green, ind_veg_wall, ind_wat_win,                       &
               surf_def_h, surf_lsm_h, surf_usm_h

    USE urban_surface_mod,                                                     &
        ONLY:  usm_3d_data_averaging


    IMPLICIT NONE

    LOGICAL      ::  match_def !< flag indicating default-type surface
    LOGICAL      ::  match_lsm !< flag indicating natural-type surface
    LOGICAL      ::  match_usm !< flag indicating urban-type surface
    
    INTEGER(iwp) ::  i   !< grid index x direction
    INTEGER(iwp) ::  ii  !< running index
    INTEGER(iwp) ::  j   !< grid index y direction
    INTEGER(iwp) ::  k   !< grid index x direction
    INTEGER(iwp) ::  m   !< running index over surfacle elements
    INTEGER(iwp) ::  n   !< running index over number of particles per grid box

    REAL(wp)     ::  mean_r !< mean-particle radius witin grid box
    REAL(wp)     ::  s_r2   !< mean-particle radius witin grid box to the power of two
    REAL(wp)     ::  s_r3   !< mean-particle radius witin grid box to the power of three

    CHARACTER (LEN=varnamelength) ::  trimvar  !< TRIM of output-variable string


    CALL cpu_log (log_point(34),'sum_up_3d_data','start')

!
!-- Allocate and initialize the summation arrays if called for the very first
!-- time or the first time after average_3d_data has been called 
!-- (some or all of the arrays may have been already allocated
!-- in rrd_local)
    IF ( average_count_3d == 0 )  THEN

       DO  ii = 1, doav_n

          trimvar = TRIM( doav(ii) )

          SELECT CASE ( trimvar )

             CASE ( 'ghf*' )
                IF ( .NOT. ALLOCATED( ghf_av ) )  THEN
                   ALLOCATE( ghf_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                ghf_av = 0.0_wp

             CASE ( 'e' )
                IF ( .NOT. ALLOCATED( e_av ) )  THEN
                   ALLOCATE( e_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                e_av = 0.0_wp

             CASE ( 'thetal' )
                IF ( .NOT. ALLOCATED( lpt_av ) )  THEN
                   ALLOCATE( lpt_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                lpt_av = 0.0_wp

             CASE ( 'lwp*' )
                IF ( .NOT. ALLOCATED( lwp_av ) )  THEN
                   ALLOCATE( lwp_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                lwp_av = 0.0_wp

             CASE ( 'ol*' )
                IF ( .NOT. ALLOCATED( ol_av ) )  THEN
                   ALLOCATE( ol_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                ol_av = 0.0_wp

             CASE ( 'p' )
                IF ( .NOT. ALLOCATED( p_av ) )  THEN
                   ALLOCATE( p_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                p_av = 0.0_wp

             CASE ( 'pc' )
                IF ( .NOT. ALLOCATED( pc_av ) )  THEN
                   ALLOCATE( pc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                pc_av = 0.0_wp

             CASE ( 'pr' )
                IF ( .NOT. ALLOCATED( pr_av ) )  THEN
                   ALLOCATE( pr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                pr_av = 0.0_wp

             CASE ( 'theta' )
                IF ( .NOT. ALLOCATED( pt_av ) )  THEN
                   ALLOCATE( pt_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                pt_av = 0.0_wp

             CASE ( 'q' )
                IF ( .NOT. ALLOCATED( q_av ) )  THEN
                   ALLOCATE( q_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                q_av = 0.0_wp

             CASE ( 'ql' )
                IF ( .NOT. ALLOCATED( ql_av ) )  THEN
                   ALLOCATE( ql_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                ql_av = 0.0_wp

             CASE ( 'ql_c' )
                IF ( .NOT. ALLOCATED( ql_c_av ) )  THEN
                   ALLOCATE( ql_c_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                ql_c_av = 0.0_wp

             CASE ( 'ql_v' )
                IF ( .NOT. ALLOCATED( ql_v_av ) )  THEN
                   ALLOCATE( ql_v_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                ql_v_av = 0.0_wp

             CASE ( 'ql_vp' )
                IF ( .NOT. ALLOCATED( ql_vp_av ) )  THEN
                   ALLOCATE( ql_vp_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                ql_vp_av = 0.0_wp

             CASE ( 'qsws*' )
                IF ( .NOT. ALLOCATED( qsws_av ) )  THEN
                   ALLOCATE( qsws_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                qsws_av = 0.0_wp

             CASE ( 'qv' )
                IF ( .NOT. ALLOCATED( qv_av ) )  THEN
                   ALLOCATE( qv_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                qv_av = 0.0_wp

             CASE ( 'r_a*' )
                IF ( .NOT. ALLOCATED( r_a_av ) )  THEN
                   ALLOCATE( r_a_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                r_a_av = 0.0_wp

             CASE ( 's' )
                IF ( .NOT. ALLOCATED( s_av ) )  THEN
                   ALLOCATE( s_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                s_av = 0.0_wp

             CASE ( 'shf*' )
                IF ( .NOT. ALLOCATED( shf_av ) )  THEN
                   ALLOCATE( shf_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                shf_av = 0.0_wp
                
             CASE ( 'ssws*' )
                IF ( .NOT. ALLOCATED( ssws_av ) )  THEN
                   ALLOCATE( ssws_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                ssws_av = 0.0_wp                

             CASE ( 't*' )
                IF ( .NOT. ALLOCATED( ts_av ) )  THEN
                   ALLOCATE( ts_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                ts_av = 0.0_wp

             CASE ( 'tsurf*' )
                IF ( .NOT. ALLOCATED( tsurf_av ) )  THEN
                   ALLOCATE( tsurf_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                tsurf_av = 0.0_wp

             CASE ( 'u' )
                IF ( .NOT. ALLOCATED( u_av ) )  THEN
                   ALLOCATE( u_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                u_av = 0.0_wp

             CASE ( 'us*' )
                IF ( .NOT. ALLOCATED( us_av ) )  THEN
                   ALLOCATE( us_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                us_av = 0.0_wp

             CASE ( 'v' )
                IF ( .NOT. ALLOCATED( v_av ) )  THEN
                   ALLOCATE( v_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                v_av = 0.0_wp

             CASE ( 'thetav' )
                IF ( .NOT. ALLOCATED( vpt_av ) )  THEN
                   ALLOCATE( vpt_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                vpt_av = 0.0_wp

             CASE ( 'w' )
                IF ( .NOT. ALLOCATED( w_av ) )  THEN
                   ALLOCATE( w_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                w_av = 0.0_wp

             CASE ( 'z0*' )
                IF ( .NOT. ALLOCATED( z0_av ) )  THEN
                   ALLOCATE( z0_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                z0_av = 0.0_wp

             CASE ( 'z0h*' )
                IF ( .NOT. ALLOCATED( z0h_av ) )  THEN
                   ALLOCATE( z0h_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                z0h_av = 0.0_wp

             CASE ( 'z0q*' )
                IF ( .NOT. ALLOCATED( z0q_av ) )  THEN
                   ALLOCATE( z0q_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                z0q_av = 0.0_wp


             CASE DEFAULT

!
!--             Allocating and initializing data arrays for all other modules
                CALL module_interface_3d_data_averaging( 'allocate', trimvar )


          END SELECT

       ENDDO

    ENDIF

!
!-- Loop of all variables to be averaged.
    DO  ii = 1, doav_n

       trimvar = TRIM( doav(ii) )
!
!--    Store the array chosen on the temporary array.
       SELECT CASE ( trimvar )

          CASE ( 'ghf*' )
             IF ( ALLOCATED( ghf_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
!
!--                   Check whether grid point is a natural- or urban-type 
!--                   surface.
                      match_lsm = surf_lsm_h%start_index(j,i) <=               &
                                  surf_lsm_h%end_index(j,i)
                      match_usm = surf_usm_h%start_index(j,i) <=               &
                                  surf_usm_h%end_index(j,i)
!
!--                   In order to avoid double-counting of surface properties,
!--                   always assume that natural-type surfaces are below urban-
!--                   type surfaces, e.g. in case of bridges.
!--                   Further, take only the last suface element, i.e. the
!--                   uppermost surface which would be visible from above
                      IF ( match_lsm  .AND.  .NOT. match_usm )  THEN
                         m = surf_lsm_h%end_index(j,i)
                         ghf_av(j,i) = ghf_av(j,i) +                           &
                                         surf_lsm_h%ghf(m)
                      ELSEIF ( match_usm )  THEN
                         m = surf_usm_h%end_index(j,i)
                         ghf_av(j,i) = ghf_av(j,i) +                           &
                                         surf_usm_h%frac(m,ind_veg_wall)  *    &
                                         surf_usm_h%wghf_eb(m)        +        &
                                         surf_usm_h%frac(m,ind_pav_green) *    &
                                         surf_usm_h%wghf_eb_green(m)  +        &
                                         surf_usm_h%frac(m,ind_wat_win)   *    &
                                         surf_usm_h%wghf_eb_window(m)
                      ENDIF
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'e' )
             IF ( ALLOCATED( e_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         e_av(k,j,i) = e_av(k,j,i) + e(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'thetal' )
             IF ( ALLOCATED( lpt_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         lpt_av(k,j,i) = lpt_av(k,j,i) + pt(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'lwp*' )
             IF ( ALLOCATED( lwp_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      lwp_av(j,i) = lwp_av(j,i) + SUM( ql(nzb:nzt,j,i)            &
                                                  * dzw(1:nzt+1) ) * rho_surface
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'ol*' )
             IF ( ALLOCATED( ol_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      match_def = surf_def_h(0)%start_index(j,i) <=            &
                                  surf_def_h(0)%end_index(j,i)
                      match_lsm = surf_lsm_h%start_index(j,i) <=               &
                                  surf_lsm_h%end_index(j,i)
                      match_usm = surf_usm_h%start_index(j,i) <=               &
                                  surf_usm_h%end_index(j,i)

                      IF ( match_def )  THEN
                         m = surf_def_h(0)%end_index(j,i)
                         ol_av(j,i) = ol_av(j,i) +                             &
                                         surf_def_h(0)%ol(m)
                      ELSEIF ( match_lsm  .AND.  .NOT. match_usm )  THEN
                         m = surf_lsm_h%end_index(j,i)
                         ol_av(j,i) = ol_av(j,i) +                             &
                                         surf_lsm_h%ol(m)
                      ELSEIF ( match_usm )  THEN
                         m = surf_usm_h%end_index(j,i)
                         ol_av(j,i) = ol_av(j,i) +                             &
                                         surf_usm_h%ol(m)
                      ENDIF
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'p' )
             IF ( ALLOCATED( p_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         p_av(k,j,i) = p_av(k,j,i) + p(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'pc' )
             IF ( ALLOCATED( pc_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         pc_av(k,j,i) = pc_av(k,j,i) + prt_count(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'pr' )
             IF ( ALLOCATED( pr_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         number_of_particles = prt_count(k,j,i)
                         IF ( number_of_particles <= 0 )  CYCLE
                         particles =>                                          &
                         grid_particles(k,j,i)%particles(1:number_of_particles)
                         s_r2 = 0.0_wp
                         s_r3 = 0.0_wp

                         DO  n = 1, number_of_particles
                            IF ( particles(n)%particle_mask )  THEN
                               s_r2 = s_r2 + particles(n)%radius**2 *          &
                                   particles(n)%weight_factor
                               s_r3 = s_r3 + particles(n)%radius**3 *          &
                                   particles(n)%weight_factor
                            ENDIF
                         ENDDO

                         IF ( s_r2 > 0.0_wp )  THEN
                            mean_r = s_r3 / s_r2
                         ELSE
                            mean_r = 0.0_wp
                         ENDIF
                         pr_av(k,j,i) = pr_av(k,j,i) + mean_r
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'theta' )
             IF ( ALLOCATED( pt_av ) ) THEN
                IF ( .NOT. bulk_cloud_model ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                            pt_av(k,j,i) = pt_av(k,j,i) + pt(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
                ELSE
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                            pt_av(k,j,i) = pt_av(k,j,i) + pt(k,j,i) + lv_d_cp * &
                                                          d_exner(k) * ql(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF
             ENDIF

          CASE ( 'q' )
             IF ( ALLOCATED( q_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         q_av(k,j,i) = q_av(k,j,i) + q(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'ql' )
             IF ( ALLOCATED( ql_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         ql_av(k,j,i) = ql_av(k,j,i) + ql(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'ql_c' )
             IF ( ALLOCATED( ql_c_av ) ) THEN 
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         ql_c_av(k,j,i) = ql_c_av(k,j,i) + ql_c(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'ql_v' )
             IF ( ALLOCATED( ql_v_av ) ) THEN 
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         ql_v_av(k,j,i) = ql_v_av(k,j,i) + ql_v(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'ql_vp' )
             IF ( ALLOCATED( ql_vp_av ) ) THEN 
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         number_of_particles = prt_count(k,j,i)
                         IF ( number_of_particles <= 0 )  CYCLE
                         particles =>                                          & 
                         grid_particles(k,j,i)%particles(1:number_of_particles)
                         DO  n = 1, number_of_particles
                            IF ( particles(n)%particle_mask )  THEN
                               ql_vp_av(k,j,i) = ql_vp_av(k,j,i) + &
                                                 particles(n)%weight_factor /  &
                                                 number_of_particles
                            ENDIF
                         ENDDO
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'qsws*' )
!
!--          In case of default surfaces, clean-up flux by density.
!--          In case of land- and urban-surfaces, convert fluxes into
!--          dynamic units.
!--          Question (maronga): are the .NOT. statements really required?
             IF ( ALLOCATED( qsws_av ) ) THEN 
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      match_def = surf_def_h(0)%start_index(j,i) <=            &
                                  surf_def_h(0)%end_index(j,i)
                      match_lsm = surf_lsm_h%start_index(j,i) <=               &
                                  surf_lsm_h%end_index(j,i)
                      match_usm = surf_usm_h%start_index(j,i) <=               &
                                  surf_usm_h%end_index(j,i)

                      IF ( match_def )  THEN
                         m = surf_def_h(0)%end_index(j,i)
                         qsws_av(j,i) = qsws_av(j,i) +                         &
                                         surf_def_h(0)%qsws(m) *               &
                                         waterflux_output_conversion(nzb)
                      ELSEIF ( match_lsm  .AND.  .NOT. match_usm )  THEN
                         m = surf_lsm_h%end_index(j,i)
                         qsws_av(j,i) = qsws_av(j,i) +                         &
                                         surf_lsm_h%qsws(m) * l_v
                      ELSEIF ( match_usm  .AND.  .NOT. match_lsm )  THEN
                         m = surf_usm_h%end_index(j,i)
                         qsws_av(j,i) = qsws_av(j,i) +                         &
                                         surf_usm_h%qsws(m) * l_v
                      ENDIF
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'qv' )
             IF ( ALLOCATED( qv_av ) ) THEN 
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         qv_av(k,j,i) = qv_av(k,j,i) + q(k,j,i) - ql(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'r_a*' )
             IF ( ALLOCATED( r_a_av ) ) THEN 
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      match_lsm = surf_lsm_h%start_index(j,i) <=               &
                                  surf_lsm_h%end_index(j,i)
                      match_usm = surf_usm_h%start_index(j,i) <=               &
                                  surf_usm_h%end_index(j,i)

                      IF ( match_lsm  .AND.  .NOT. match_usm )  THEN
                         m = surf_lsm_h%end_index(j,i)
                         r_a_av(j,i) = r_a_av(j,i) +                           &
                                         surf_lsm_h%r_a(m)
                      ELSEIF ( match_usm )  THEN
                         m = surf_usm_h%end_index(j,i)
                         r_a_av(j,i) = r_a_av(j,i) +                           &
                                         surf_usm_h%frac(m,ind_veg_wall)  *    &
                                         surf_usm_h%r_a(m)       +             & 
                                         surf_usm_h%frac(m,ind_pav_green) *    &
                                         surf_usm_h%r_a_green(m) +             & 
                                         surf_usm_h%frac(m,ind_wat_win)   *    &
                                         surf_usm_h%r_a_window(m)
                      ENDIF
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 's' )
             IF ( ALLOCATED( s_av ) ) THEN 
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         s_av(k,j,i) = s_av(k,j,i) + s(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'shf*' )
!
!--          In case of default surfaces, clean-up flux by density.
!--          In case of land- and urban-surfaces, convert fluxes into
!--          dynamic units.
             IF ( ALLOCATED( shf_av ) ) THEN 
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      match_def = surf_def_h(0)%start_index(j,i) <=            &
                                  surf_def_h(0)%end_index(j,i)
                      match_lsm = surf_lsm_h%start_index(j,i) <=               &
                                  surf_lsm_h%end_index(j,i)
                      match_usm = surf_usm_h%start_index(j,i) <=               &
                                  surf_usm_h%end_index(j,i)

                      IF ( match_def )  THEN
                         m = surf_def_h(0)%end_index(j,i)
                         shf_av(j,i) = shf_av(j,i) +                           &
                                         surf_def_h(0)%shf(m)  *               &
                                         heatflux_output_conversion(nzb)
                      ELSEIF ( match_lsm  .AND.  .NOT. match_usm )  THEN
                         m = surf_lsm_h%end_index(j,i)
                         shf_av(j,i) = shf_av(j,i) +                           &
                                         surf_lsm_h%shf(m) * c_p
                      ELSEIF ( match_usm )  THEN
                         m = surf_usm_h%end_index(j,i)
                         shf_av(j,i) = shf_av(j,i) +                           &
                                         surf_usm_h%shf(m) * c_p
                      ENDIF
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'ssws*' )
             IF ( ALLOCATED( ssws_av ) ) THEN 
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      match_def = surf_def_h(0)%start_index(j,i) <=            &
                                  surf_def_h(0)%end_index(j,i)
                      match_lsm = surf_lsm_h%start_index(j,i) <=               &
                                  surf_lsm_h%end_index(j,i)
                      match_usm = surf_usm_h%start_index(j,i) <=               &
                                  surf_usm_h%end_index(j,i)

                      IF ( match_def )  THEN
                         m = surf_def_h(0)%end_index(j,i)
                         ssws_av(j,i) = ssws_av(j,i) +                         &
                                         surf_def_h(0)%ssws(m)
                      ELSEIF ( match_lsm  .AND.  .NOT. match_usm )  THEN
                         m = surf_lsm_h%end_index(j,i)
                         ssws_av(j,i) = ssws_av(j,i) +                         &
                                         surf_lsm_h%ssws(m)
                      ELSEIF ( match_usm )  THEN
                         m = surf_usm_h%end_index(j,i)
                         ssws_av(j,i) = ssws_av(j,i) +                         &
                                         surf_usm_h%ssws(m)
                      ENDIF
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 't*' )
             IF ( ALLOCATED( ts_av ) ) THEN 
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      match_def = surf_def_h(0)%start_index(j,i) <=            &
                                  surf_def_h(0)%end_index(j,i)
                      match_lsm = surf_lsm_h%start_index(j,i) <=               &
                                  surf_lsm_h%end_index(j,i)
                      match_usm = surf_usm_h%start_index(j,i) <=               &
                                  surf_usm_h%end_index(j,i)

                      IF ( match_def )  THEN
                         m = surf_def_h(0)%end_index(j,i)
                         ts_av(j,i) = ts_av(j,i) +                             &
                                         surf_def_h(0)%ts(m)
                      ELSEIF ( match_lsm  .AND.  .NOT. match_usm )  THEN
                         m = surf_lsm_h%end_index(j,i)
                         ts_av(j,i) = ts_av(j,i) +                             &
                                         surf_lsm_h%ts(m)
                      ELSEIF ( match_usm )  THEN
                         m = surf_usm_h%end_index(j,i)
                         ts_av(j,i) = ts_av(j,i) +                             &
                                         surf_usm_h%ts(m)
                      ENDIF
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'tsurf*' )
             IF ( ALLOCATED( tsurf_av ) ) THEN    
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      match_def = surf_def_h(0)%start_index(j,i) <=            &
                                  surf_def_h(0)%end_index(j,i)
                      match_lsm = surf_lsm_h%start_index(j,i) <=               &
                                  surf_lsm_h%end_index(j,i)
                      match_usm = surf_usm_h%start_index(j,i) <=               &
                                  surf_usm_h%end_index(j,i)

                      IF ( match_def )  THEN
                         m = surf_def_h(0)%end_index(j,i)
                         tsurf_av(j,i) = tsurf_av(j,i) +                       &
                                         surf_def_h(0)%pt_surface(m)
                      ELSEIF ( match_lsm  .AND.  .NOT. match_usm )  THEN
                         m = surf_lsm_h%end_index(j,i)
                         tsurf_av(j,i) = tsurf_av(j,i) +                       &
                                         surf_lsm_h%pt_surface(m)
                      ELSEIF ( match_usm )  THEN
                         m = surf_usm_h%end_index(j,i)
                         tsurf_av(j,i) = tsurf_av(j,i) +                       &
                                         surf_usm_h%pt_surface(m)
                      ENDIF
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'u' )
             IF ( ALLOCATED( u_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         u_av(k,j,i) = u_av(k,j,i) + u(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'us*' )
             IF ( ALLOCATED( us_av ) ) THEN   
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      match_def = surf_def_h(0)%start_index(j,i) <=            &
                                  surf_def_h(0)%end_index(j,i)
                      match_lsm = surf_lsm_h%start_index(j,i) <=               &
                                  surf_lsm_h%end_index(j,i)
                      match_usm = surf_usm_h%start_index(j,i) <=               &
                                  surf_usm_h%end_index(j,i)

                      IF ( match_def )  THEN
                         m = surf_def_h(0)%end_index(j,i)
                         us_av(j,i) = us_av(j,i) +                             &
                                         surf_def_h(0)%us(m)
                      ELSEIF ( match_lsm  .AND.  .NOT. match_usm )  THEN
                         m = surf_lsm_h%end_index(j,i)
                         us_av(j,i) = us_av(j,i) +                             &
                                         surf_lsm_h%us(m)
                      ELSEIF ( match_usm )  THEN
                         m = surf_usm_h%end_index(j,i)
                         us_av(j,i) = us_av(j,i) +                             &
                                         surf_usm_h%us(m)
                      ENDIF
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'v' )
             IF ( ALLOCATED( v_av ) ) THEN 
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         v_av(k,j,i) = v_av(k,j,i) + v(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'thetav' )
             IF ( ALLOCATED( vpt_av ) ) THEN 
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         vpt_av(k,j,i) = vpt_av(k,j,i) + vpt(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'w' )
             IF ( ALLOCATED( w_av ) ) THEN 
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         w_av(k,j,i) = w_av(k,j,i) + w(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'z0*' )
             IF ( ALLOCATED( z0_av ) ) THEN 
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      match_def = surf_def_h(0)%start_index(j,i) <=            &
                                  surf_def_h(0)%end_index(j,i)
                      match_lsm = surf_lsm_h%start_index(j,i) <=               &
                                  surf_lsm_h%end_index(j,i)
                      match_usm = surf_usm_h%start_index(j,i) <=               &
                                  surf_usm_h%end_index(j,i)

                      IF ( match_def )  THEN
                         m = surf_def_h(0)%end_index(j,i)
                         z0_av(j,i) = z0_av(j,i) +                             &
                                         surf_def_h(0)%z0(m)
                      ELSEIF ( match_lsm  .AND.  .NOT. match_usm )  THEN
                         m = surf_lsm_h%end_index(j,i)
                         z0_av(j,i) = z0_av(j,i) +                             &
                                         surf_lsm_h%z0(m)
                      ELSEIF ( match_usm )  THEN
                         m = surf_usm_h%end_index(j,i)
                         z0_av(j,i) = z0_av(j,i) +                             &
                                         surf_usm_h%z0(m)
                      ENDIF
                   ENDDO
                ENDDO   
             ENDIF

          CASE ( 'z0h*' )
             IF ( ALLOCATED( z0h_av ) ) THEN 
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      match_def = surf_def_h(0)%start_index(j,i) <=            &
                                  surf_def_h(0)%end_index(j,i)
                      match_lsm = surf_lsm_h%start_index(j,i) <=               &
                                  surf_lsm_h%end_index(j,i)
                      match_usm = surf_usm_h%start_index(j,i) <=               &
                                  surf_usm_h%end_index(j,i)

                      IF ( match_def )  THEN
                         m = surf_def_h(0)%end_index(j,i)
                         z0h_av(j,i) = z0h_av(j,i) +                           &
                                         surf_def_h(0)%z0h(m)
                      ELSEIF ( match_lsm  .AND.  .NOT. match_usm )  THEN
                         m = surf_lsm_h%end_index(j,i)
                         z0h_av(j,i) = z0h_av(j,i) +                           &
                                         surf_lsm_h%z0h(m)
                      ELSEIF ( match_usm )  THEN
                         m = surf_usm_h%end_index(j,i)
                         z0h_av(j,i) = z0h_av(j,i) +                           &
                                         surf_usm_h%z0h(m)
                      ENDIF
                   ENDDO
                ENDDO
             ENDIF
   
          CASE ( 'z0q*' )
             IF ( ALLOCATED( z0q_av ) ) THEN 
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      match_def = surf_def_h(0)%start_index(j,i) <=            &
                                  surf_def_h(0)%end_index(j,i)
                      match_lsm = surf_lsm_h%start_index(j,i) <=               &
                                  surf_lsm_h%end_index(j,i)
                      match_usm = surf_usm_h%start_index(j,i) <=               &
                                  surf_usm_h%end_index(j,i)

                      IF ( match_def )  THEN
                         m = surf_def_h(0)%end_index(j,i)
                         z0q_av(j,i) = z0q_av(j,i) +                           &
                                         surf_def_h(0)%z0q(m)
                      ELSEIF ( match_lsm  .AND.  .NOT. match_usm )  THEN
                         m = surf_lsm_h%end_index(j,i)
                         z0q_av(j,i) = z0q_av(j,i) +                           &
                                         surf_lsm_h%z0q(m)
                      ELSEIF ( match_usm )  THEN
                         m = surf_usm_h%end_index(j,i)
                         z0q_av(j,i) = z0q_av(j,i) +                           &
                                         surf_usm_h%z0q(m)
                      ENDIF
                   ENDDO
                ENDDO
             ENDIF

          CASE DEFAULT

!--          In case of urban surface variables it should be always checked
!--          if respective arrays are allocated, at least in case of a restart
!--          run, as averaged usm arrays are not read from file at the moment.
             IF ( urban_surface )  THEN
                CALL usm_3d_data_averaging( 'allocate', trimvar )
             ENDIF

!
!--          Summing up data from all other modules
             CALL module_interface_3d_data_averaging( 'sum', trimvar )


       END SELECT

    ENDDO

    CALL cpu_log( log_point(34), 'sum_up_3d_data', 'stop' )


 END SUBROUTINE sum_up_3d_data
