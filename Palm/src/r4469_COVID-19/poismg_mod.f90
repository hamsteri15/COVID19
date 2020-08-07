!> @file poismg.f90
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
! $Id: poismg_mod.f90 4457 2020-03-11 14:20:43Z raasch $
! use statement for exchange horiz added
! 
! 4432 2020-02-28 07:43:21Z raasch
! bugfix for previous revision (vector directive was changed by mistake)
! 
! 4429 2020-02-27 15:24:30Z raasch
! statement added to avoid compile error due to unused dummy argument
! bugfix: cpp-directives added for serial mode
! 
! 4360 2020-01-07 11:25:50Z suehring
! Corrected "Former revisions" section
! 
! 3725 2019-02-07 10:11:02Z raasch
! unused subroutine removed
! 
! 3655 2019-01-07 16:51:22Z knoop
! unnecessary check eliminated
!
! Following optimisations have been made:
! - vectorisation (for Intel-CPUs) of the red-black algorithm by resorting
!   array elements with even and odd indices
! - explicit boundary conditions for building walls removed (solver is
!   running through the buildings
! - reduced data transfer in case of ghost point exchange, because only
!   "red" or "black" data points need to be exchanged. This is not applied
!   for coarser grid levels, since for then the transfer time is latency bound
!
!
! Description:
! ------------
!> Solves the Poisson equation for the perturbation pressure with a multigrid
!> V- or W-Cycle scheme.
!>
!> This multigrid method was originally developed for PALM by Joerg Uhlenbrock,
!> September 2000 - July 2001. It has been optimised for speed by Klaus
!> Ketelsen in November 2014.
!> 
!> @attention Loop unrolling and cache optimization in SOR-Red/Black method
!>            still does not give the expected speedup! 
!>
!> @todo Further work required.
!------------------------------------------------------------------------------!
 MODULE poismg_mod
 
    USE control_parameters,                                                    &
        ONLY:  bc_dirichlet_l, bc_dirichlet_n, bc_dirichlet_r,                 &
               bc_dirichlet_s, bc_radiation_l, bc_radiation_n, bc_radiation_r, &
               bc_radiation_s, grid_level, nesting_offline

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point_s

    USE exchange_horiz_mod,                                                    &
        ONLY:  exchange_horiz

    USE kinds

    USE pegrid

    PRIVATE

    INTEGER, SAVE                             ::  ind_even_odd    !< border index between even and odd k index
    INTEGER, DIMENSION(:), SAVE, ALLOCATABLE  ::  even_odd_level  !< stores ind_even_odd for all MG levels

    REAL(wp), DIMENSION(:,:), SAVE, ALLOCATABLE ::  f1_mg_b, f2_mg_b, f3_mg_b  !< blocked version of f1_mg ...

    INTERFACE poismg
       MODULE PROCEDURE poismg
    END INTERFACE poismg

    INTERFACE sort_k_to_even_odd_blocks
       MODULE PROCEDURE sort_k_to_even_odd_blocks
!       MODULE PROCEDURE sort_k_to_even_odd_blocks_int
       MODULE PROCEDURE sort_k_to_even_odd_blocks_1d
    END INTERFACE sort_k_to_even_odd_blocks

    PUBLIC poismg

 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Solves the Poisson equation for the perturbation pressure with a multigrid
!> V- or W-Cycle scheme.
!------------------------------------------------------------------------------!
    SUBROUTINE poismg( r )

       USE arrays_3d,                                                          &
           ONLY:  d, p_loc

       USE control_parameters,                                                 &
           ONLY:  bc_lr_cyc, bc_ns_cyc, gathered_size, grid_level,             &
                  grid_level_count, ibc_p_t,                                   &
                  maximum_grid_level, message_string, mgcycles, mg_cycles,     &
                  mg_switch_to_pe0_level, residual_limit, subdomain_size

       USE cpulog,                                                             &
           ONLY:  cpu_log, log_point_s

       USE indices,                                                            &
           ONLY:  nxl, nxlg, nxl_mg, nxr, nxrg, nxr_mg, nys, nysg, nys_mg, nyn,&
                  nyng, nyn_mg, nzb, nzt, nzt_mg

       IMPLICIT NONE

       REAL(wp) ::  maxerror          !<
       REAL(wp) ::  maximum_mgcycles  !<
       REAL(wp) ::  residual_norm     !<

       REAL(wp), DIMENSION(nzb:nzt+1,nys-1:nyn+1,nxl-1:nxr+1) ::  r  !<

       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  p3  !<


       CALL cpu_log( log_point_s(29), 'poismg', 'start' )
!
!--    Initialize arrays and variables used in this subroutine

!--    If the number of grid points of the gathered grid, which is collected
!--    on PE0, is larger than the number of grid points of an PE, than array
!--    p3 will be enlarged.
       IF ( gathered_size > subdomain_size )  THEN
          ALLOCATE( p3(nzb:nzt_mg(mg_switch_to_pe0_level)+1,nys_mg(            &
                    mg_switch_to_pe0_level)-1:nyn_mg(mg_switch_to_pe0_level)+1,&
                    nxl_mg(mg_switch_to_pe0_level)-1:nxr_mg(                   &
                    mg_switch_to_pe0_level)+1) )
       ELSE
          ALLOCATE ( p3(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       ENDIF

       p3 = 0.0_wp

 
!
!--    Ghost boundaries have to be added to divergence array.
!--    Exchange routine needs to know the grid level!
       grid_level = maximum_grid_level
       CALL exchange_horiz( d, 1)
!
!--    Set bottom and top boundary conditions
       d(nzb,:,:) = d(nzb+1,:,:)
       IF ( ibc_p_t == 1 )  d(nzt+1,:,: ) = d(nzt,:,:)
!
!--    Set lateral boundary conditions in non-cyclic case
       IF ( .NOT. bc_lr_cyc )  THEN
          IF ( bc_dirichlet_l  .OR.  bc_radiation_l )                          &
             d(:,:,nxl-1) = d(:,:,nxl)
          IF ( bc_dirichlet_r  .OR.  bc_radiation_r )                          &
             d(:,:,nxr+1) = d(:,:,nxr)
       ENDIF
       IF ( .NOT. bc_ns_cyc )  THEN
          IF ( bc_dirichlet_n  .OR.  bc_radiation_n )                          &
             d(:,nyn+1,:) = d(:,nyn,:)
          IF ( bc_dirichlet_s  .OR.  bc_radiation_s )                          &
             d(:,nys-1,:) = d(:,nys,:)
       ENDIF
!
!--    Initiation of the multigrid scheme. Does n cycles until the
!--    residual is smaller than the given limit. The accuracy of the solution
!--    of the poisson equation will increase with the number of cycles.
!--    If the number of cycles is preset by the user, this number will be
!--    carried out regardless of the accuracy.
       grid_level_count =  0
       mgcycles         =  0
       IF ( mg_cycles == -1 )  THEN
          maximum_mgcycles = 0
          residual_norm    = 1.0_wp
       ELSE
          maximum_mgcycles = mg_cycles
          residual_norm    = 0.0_wp
       ENDIF

!
!--    Initial settings for sorting k-dimension from sequential order (alternate
!--    even/odd) into blocks of even and odd or vice versa
       CALL init_even_odd_blocks

!
!--    Sort input arrays in even/odd blocks along k-dimension
       CALL sort_k_to_even_odd_blocks( d, grid_level )
       CALL sort_k_to_even_odd_blocks( p_loc, grid_level )

!
!--    The complete multigrid cycles are running in block mode, i.e. over
!--    seperate data blocks of even and odd indices
       DO WHILE ( residual_norm > residual_limit  .OR. &
                  mgcycles < maximum_mgcycles )
 
          CALL next_mg_level( d, p_loc, p3, r)

!
!--       Calculate the residual if the user has not preset the number of
!--       cycles to be performed
          IF ( maximum_mgcycles == 0 )  THEN
             CALL resid( d, p_loc, r )
             maxerror = SUM( r(nzb+1:nzt,nys:nyn,nxl:nxr)**2 )

#if defined( __parallel )
             IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
                CALL MPI_ALLREDUCE( maxerror, residual_norm, 1, MPI_REAL,      &
                                    MPI_SUM, comm2d, ierr)
#else
                residual_norm = maxerror
#endif
             residual_norm = SQRT( residual_norm )
          ENDIF

          mgcycles = mgcycles + 1

!
!--       If the user has not limited the number of cycles, stop the run in case
!--       of insufficient convergence
          IF ( mgcycles > 1000  .AND.  mg_cycles == -1 )  THEN
             message_string = 'no sufficient convergence within 1000 cycles'
             CALL message( 'poismg', 'PA0283', 1, 2, 0, 6, 0 )
          ENDIF

       ENDDO

       DEALLOCATE( p3 )
!
!--    Result has to be sorted back from even/odd blocks to sequential order
       CALL sort_k_to_sequential( p_loc )
!
!--    Unset the grid level. Variable is used to determine the MPI datatypes for
!--    ghost point exchange
       grid_level = 0

       CALL cpu_log( log_point_s(29), 'poismg', 'stop' )

    END SUBROUTINE poismg


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Computes the residual of the perturbation pressure.
!------------------------------------------------------------------------------!
    SUBROUTINE resid( f_mg, p_mg, r )


       USE arrays_3d,                                                          &
           ONLY:  rho_air_mg

       USE control_parameters,                                                 &
           ONLY:  bc_lr_cyc, bc_ns_cyc, ibc_p_b, ibc_p_t
       USE grid_variables,                                                     &
           ONLY:  ddx2_mg, ddy2_mg

       USE indices,                                                            &
           ONLY:  nxl_mg, nxr_mg, nys_mg, nyn_mg, nzb, nzt_mg

       IMPLICIT NONE

       INTEGER(iwp) ::  i        !< index variable along x
       INTEGER(iwp) ::  j        !< index variable along y
       INTEGER(iwp) ::  k        !< index variable along z
       INTEGER(iwp) ::  l        !< index indicating grid level
       INTEGER(iwp) ::  km1      !< index variable along z dimension (k-1)
       INTEGER(iwp) ::  kp1      !< index variable along z dimension (k+1)

       REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,                           &
                       nys_mg(grid_level)-1:nyn_mg(grid_level)+1,              &
                       nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) ::  f_mg  !< velocity divergence
       REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,                           &
                       nys_mg(grid_level)-1:nyn_mg(grid_level)+1,              &
                       nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) ::  p_mg  !< perturbation pressure
       REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,                           &
                       nys_mg(grid_level)-1:nyn_mg(grid_level)+1,              &
                       nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) ::  r     !< residuum of perturbation pressure

!
!--    Calculate the residual
       l = grid_level

       CALL cpu_log( log_point_s(53), 'resid', 'start' )
       !$OMP PARALLEL PRIVATE (i,j,k,km1,kp1)
       !$OMP DO
       DO  i = nxl_mg(l), nxr_mg(l)
          DO  j = nys_mg(l), nyn_mg(l)
                !DIR$ IVDEP
             DO k = ind_even_odd+1, nzt_mg(l)
                km1 = k-ind_even_odd-1
                kp1 = k-ind_even_odd
                r(k,j,i) = f_mg(k,j,i)                                         &
                      - rho_air_mg(k,l) * ddx2_mg(l) *                         &
                      ( p_mg(k,j,i+1) +  p_mg(k,j,i-1)  )                      &
                      - rho_air_mg(k,l) * ddy2_mg(l) *                         &
                      ( p_mg(k,j+1,i) + p_mg(k,j-1,i)  )                       &
                      - f2_mg_b(k,l) * p_mg(kp1,j,i)                           &
                      - f3_mg_b(k,l) * p_mg(km1,j,i)                           &
                      + f1_mg_b(k,l) * p_mg(k,j,i)
             ENDDO
             !DIR$ IVDEP
             DO k = nzb+1, ind_even_odd
                km1 = k+ind_even_odd
                kp1 = k+ind_even_odd+1
                r(k,j,i) = f_mg(k,j,i)                                         &
                      - rho_air_mg(k,l) * ddx2_mg(l) *                         &
                      ( p_mg(k,j,i+1) +  p_mg(k,j,i-1)  )                      &
                      - rho_air_mg(k,l) * ddy2_mg(l) *                         &
                      ( p_mg(k,j+1,i) + p_mg(k,j-1,i)  )                       &
                      - f2_mg_b(k,l) * p_mg(kp1,j,i)                           &
                      - f3_mg_b(k,l) * p_mg(km1,j,i)                           &
                      + f1_mg_b(k,l) * p_mg(k,j,i)
             ENDDO
          ENDDO
       ENDDO
       !$OMP END PARALLEL
!
!--    Horizontal boundary conditions
       CALL exchange_horiz( r, 1)

       IF ( .NOT. bc_lr_cyc )  THEN
          IF ( bc_dirichlet_l  .OR.  bc_radiation_l )  THEN
             r(:,:,nxl_mg(l)-1) = r(:,:,nxl_mg(l))
          ENDIF
          IF ( bc_dirichlet_r  .OR.  bc_radiation_r )  THEN
             r(:,:,nxr_mg(l)+1) = r(:,:,nxr_mg(l))
          ENDIF
       ENDIF

       IF ( .NOT. bc_ns_cyc )  THEN
          IF ( bc_dirichlet_n  .OR.  bc_radiation_n )  THEN
             r(:,nyn_mg(l)+1,:) = r(:,nyn_mg(l),:)
          ENDIF
          IF ( bc_dirichlet_s  .OR.  bc_radiation_s )  THEN
             r(:,nys_mg(l)-1,:) = r(:,nys_mg(l),:)
          ENDIF
       ENDIF

!
!--    Boundary conditions at bottom and top of the domain. Points may be within
!--    buildings, but that doesn't matter.
       IF ( ibc_p_b == 1 )  THEN
!
!--       equivalent to r(nzb,:,: ) = r(nzb+1,:,:)
          r(nzb,:,: ) = r(ind_even_odd+1,:,:)
       ELSE
          r(nzb,:,: ) = 0.0_wp
       ENDIF

       IF ( ibc_p_t == 1 )  THEN
!
!--       equivalent to r(nzt_mg(l)+1,:,: ) = r(nzt_mg(l),:,:)
          r(nzt_mg(l)+1,:,: ) = r(ind_even_odd,:,:)
       ELSE
          r(nzt_mg(l)+1,:,: ) = 0.0_wp
       ENDIF

       CALL cpu_log( log_point_s(53), 'resid', 'stop' )

    END SUBROUTINE resid


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Interpolates the residual on the next coarser grid with "full weighting"
!> scheme
!------------------------------------------------------------------------------!
    SUBROUTINE restrict( f_mg, r )


       USE control_parameters,                                                 &
           ONLY:  bc_lr_cyc, bc_ns_cyc, ibc_p_b, ibc_p_t

       USE indices,                                                            &
           ONLY:  nxl_mg, nxr_mg, nys_mg, nyn_mg, nzb, nzt_mg

       IMPLICIT NONE

       INTEGER(iwp) ::  i    !< index variable along x on finer grid
       INTEGER(iwp) ::  ic   !< index variable along x on coarser grid
       INTEGER(iwp) ::  j    !< index variable along y on finer grid
       INTEGER(iwp) ::  jc   !< index variable along y on coarser grid
       INTEGER(iwp) ::  k    !< index variable along z on finer grid
       INTEGER(iwp) ::  kc   !< index variable along z on coarser grid
       INTEGER(iwp) ::  l    !< index indicating finer grid level
       INTEGER(iwp) ::  km1  !< index variable along z dimension (k-1 on finer level) 
       INTEGER(iwp) ::  kp1  !< index variable along z dimension (k+1 on finer level)


       REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,                           &
                           nys_mg(grid_level)-1:nyn_mg(grid_level)+1,          &
                           nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) ::       &
                                         f_mg  !< Residual on coarser grid level

       REAL(wp), DIMENSION(nzb:nzt_mg(grid_level+1)+1,                         &
                           nys_mg(grid_level+1)-1:nyn_mg(grid_level+1)+1,      &
                           nxl_mg(grid_level+1)-1:nxr_mg(grid_level+1)+1) ::   &
                                         r !< Residual on finer grid level

!
!--    Interpolate the residual
       l = grid_level

       CALL cpu_log( log_point_s(54), 'restrict', 'start' )
!
!--    No wall treatment
       !$OMP PARALLEL PRIVATE (i,j,k,ic,jc,kc,km1,kp1)
       !$OMP DO SCHEDULE( STATIC )
       DO  ic = nxl_mg(l), nxr_mg(l)
          i = 2*ic
          DO  jc = nys_mg(l), nyn_mg(l)
!
!--          Calculation for the first point along k
             j  = 2*jc
!
!--          Calculation for the other points along k
             !DIR$ IVDEP
             DO k = ind_even_odd+1, nzt_mg(l+1)    ! Fine grid at this point
                km1 = k-ind_even_odd-1
                kp1 = k-ind_even_odd
                kc  = k-ind_even_odd               ! Coarse grid index

                f_mg(kc,jc,ic) = 1.0_wp / 64.0_wp * (                      &
                               8.0_wp * r(k,j,i)                            &
                             + 4.0_wp * ( r(k,j,i-1)     + r(k,j,i+1)     + &
                                          r(k,j+1,i)     + r(k,j-1,i)     ) &
                             + 2.0_wp * ( r(k,j-1,i-1)   + r(k,j+1,i-1)   + &
                                          r(k,j-1,i+1)   + r(k,j+1,i+1)   ) &
                             + 4.0_wp * r(km1,j,i)                          &
                             + 2.0_wp * ( r(km1,j,i-1)   + r(km1,j,i+1)   + &
                                          r(km1,j+1,i)   + r(km1,j-1,i)   ) &
                             +          ( r(km1,j-1,i-1) + r(km1,j+1,i-1) + &
                                          r(km1,j-1,i+1) + r(km1,j+1,i+1) ) &
                             + 4.0_wp * r(kp1,j,i)                          &
                             + 2.0_wp * ( r(kp1,j,i-1)   + r(kp1,j,i+1)   + &
                                          r(kp1,j+1,i)   + r(kp1,j-1,i)   ) &
                             +          ( r(kp1,j-1,i-1) + r(kp1,j+1,i-1) + &
                                          r(kp1,j-1,i+1) + r(kp1,j+1,i+1) ) &
                                        )
             ENDDO
          ENDDO
       ENDDO
       !$OMP ENDDO
       !$OMP END PARALLEL

!
!--    Ghost point exchange
       CALL exchange_horiz( f_mg, 1)
!
!--    Horizontal boundary conditions
       IF ( .NOT. bc_lr_cyc )  THEN
          IF ( bc_dirichlet_l  .OR.  bc_radiation_l )  THEN
             f_mg(:,:,nxl_mg(l)-1) = f_mg(:,:,nxl_mg(l))
          ENDIF
          IF ( bc_dirichlet_r  .OR.  bc_radiation_r )  THEN
             f_mg(:,:,nxr_mg(l)+1) = f_mg(:,:,nxr_mg(l))
          ENDIF
       ENDIF

       IF ( .NOT. bc_ns_cyc )  THEN
          IF ( bc_dirichlet_n  .OR.  bc_radiation_n )  THEN
             f_mg(:,nyn_mg(l)+1,:) = f_mg(:,nyn_mg(l),:)
          ENDIF
          IF ( bc_dirichlet_s  .OR.  bc_radiation_s )  THEN
             f_mg(:,nys_mg(l)-1,:) = f_mg(:,nys_mg(l),:)
          ENDIF
       ENDIF

!
!--    Boundary conditions at bottom and top of the domain.
!--    These points are not handled by the above loop. Points may be within
!--    buildings, but that doesn't matter. Remark: f_mg is ordered sequentielly
!--    after interpolation on coarse grid (is ordered in odd-even blocks further
!--    below). 
       IF ( ibc_p_b == 1 )  THEN
          f_mg(nzb,:,: ) = f_mg(nzb+1,:,:)
       ELSE
          f_mg(nzb,:,: ) = 0.0_wp
       ENDIF

       IF ( ibc_p_t == 1 )  THEN
          f_mg(nzt_mg(l)+1,:,: ) = f_mg(nzt_mg(l),:,:)
       ELSE
          f_mg(nzt_mg(l)+1,:,: ) = 0.0_wp
       ENDIF

       CALL cpu_log( log_point_s(54), 'restrict', 'stop' )
!
!--    Since residual is in sequential order after interpolation, an additional 
!--    sorting in odd-even blocks along z dimension is required at this point.
       CALL sort_k_to_even_odd_blocks( f_mg , l)

    END SUBROUTINE restrict


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Interpolates the correction of the perturbation pressure 
!> to the next finer grid.
!------------------------------------------------------------------------------!
    SUBROUTINE prolong( p, temp )


       USE control_parameters,                                                 &
           ONLY:  bc_lr_cyc, bc_ns_cyc, ibc_p_b, ibc_p_t
       USE indices,                                                            &
           ONLY:  nxl_mg, nxr_mg, nys_mg, nyn_mg, nzb, nzt_mg

       IMPLICIT NONE

       INTEGER(iwp) ::  i   !< index variable along x on coarser grid level
       INTEGER(iwp) ::  j   !< index variable along y on coarser grid level
       INTEGER(iwp) ::  k   !< index variable along z on coarser grid level
       INTEGER(iwp) ::  l   !< index indicating finer grid level
       INTEGER(iwp) ::  kp1 !< index variable along z
       INTEGER(iwp) ::  ke  !< index for prolog even
       INTEGER(iwp) ::  ko  !< index for prolog odd

       REAL(wp), DIMENSION(nzb:nzt_mg(grid_level-1)+1,                         &
                           nys_mg(grid_level-1)-1:nyn_mg(grid_level-1)+1,      &
                           nxl_mg(grid_level-1)-1:nxr_mg(grid_level-1)+1 ) ::  &
                               p     !< perturbation pressure on coarser grid level

       REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,                           &
                           nys_mg(grid_level)-1:nyn_mg(grid_level)+1,          &
                           nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) ::       &
                               temp  !< perturbation pressure on finer grid level


       CALL cpu_log( log_point_s(55), 'prolong', 'start' )

!
!--    First, store elements of the coarser grid on the next finer grid
       l = grid_level
       ind_even_odd = even_odd_level(grid_level-1)

       !$OMP PARALLEL PRIVATE (i,j,k,kp1,ke,ko)
       !$OMP DO
       DO  i = nxl_mg(l-1), nxr_mg(l-1)
          DO  j = nys_mg(l-1), nyn_mg(l-1)

             !DIR$ IVDEP
             DO k = ind_even_odd+1, nzt_mg(l-1)
                kp1 = k - ind_even_odd
                ke  = 2 * ( k-ind_even_odd - 1 ) + 1
                ko  = 2 * k - 1
!
!--             Points of the coarse grid are directly stored on the next finer
!--             grid
                temp(ko,2*j,2*i)   = p(k,j,i)
!
!--             Points between two coarse-grid points
                temp(ko,2*j,2*i+1) = 0.5_wp * ( p(k,j,i) + p(k,j,i+1) )
                temp(ko,2*j+1,2*i) = 0.5_wp * ( p(k,j,i) + p(k,j+1,i) )
                temp(ke,2*j,2*i)   = 0.5_wp * ( p(k,j,i) + p(kp1,j,i) )
!
!--             Points in the center of the planes stretched by four points
!--             of the coarse grid cube
                temp(ko,2*j+1,2*i+1) = 0.25_wp * ( p(k,j,i)   + p(k,j,i+1) +   &
                                                   p(k,j+1,i) + p(k,j+1,i+1) )
                temp(ke,2*j,2*i+1)   = 0.25_wp * ( p(k,j,i)   + p(k,j,i+1) +   &
                                                   p(kp1,j,i) + p(kp1,j,i+1) )
                temp(ke,2*j+1,2*i)   = 0.25_wp * ( p(k,j,i)   + p(k,j+1,i) +   &
                                                   p(kp1,j,i) + p(kp1,j+1,i) )
!
!--             Points in the middle of coarse grid cube
                temp(ke,2*j+1,2*i+1) = 0.125_wp *                              &
                                               ( p(k,j,i)     + p(k,j,i+1)   + &
                                                 p(k,j+1,i)   + p(k,j+1,i+1) + &
                                                 p(kp1,j,i)   + p(kp1,j,i+1) + &
                                                 p(kp1,j+1,i) + p(kp1,j+1,i+1) )

             ENDDO

             !DIR$ IVDEP
             DO k = nzb+1, ind_even_odd
                kp1 = k + ind_even_odd + 1
                ke  = 2 * k
                ko  = 2 * ( k + ind_even_odd )
!
!--             Points of the coarse grid are directly stored on the next finer
!--             grid
                temp(ko,2*j,2*i)   = p(k,j,i)
!
!--             Points between two coarse-grid points
                temp(ko,2*j,2*i+1) = 0.5_wp * ( p(k,j,i) + p(k,j,i+1) )
                temp(ko,2*j+1,2*i) = 0.5_wp * ( p(k,j,i) + p(k,j+1,i) )
                temp(ke,2*j,2*i)   = 0.5_wp * ( p(k,j,i) + p(kp1,j,i) )
!
!--             Points in the center of the planes stretched by four points
!--             of the coarse grid cube
                temp(ko,2*j+1,2*i+1) = 0.25_wp * ( p(k,j,i)   + p(k,j,i+1) +   &
                                                   p(k,j+1,i) + p(k,j+1,i+1) )
                temp(ke,2*j,2*i+1)   = 0.25_wp * ( p(k,j,i)   + p(k,j,i+1) +   &
                                                   p(kp1,j,i) + p(kp1,j,i+1) )
                temp(ke,2*j+1,2*i)   = 0.25_wp * ( p(k,j,i)   + p(k,j+1,i) +   &
                                                   p(kp1,j,i) + p(kp1,j+1,i) )
!
!--             Points in the middle of coarse grid cube
                temp(ke,2*j+1,2*i+1) = 0.125_wp *                              &
                                               ( p(k,j,i)     + p(k,j,i+1)   + &
                                                 p(k,j+1,i)   + p(k,j+1,i+1) + &
                                                 p(kp1,j,i)   + p(kp1,j,i+1) + &
                                                 p(kp1,j+1,i) + p(kp1,j+1,i+1) )

             ENDDO

          ENDDO
       ENDDO
       !$OMP END PARALLEL

       ind_even_odd = even_odd_level(grid_level)
!
!--    Horizontal boundary conditions
       CALL exchange_horiz( temp, 1)

       IF ( .NOT. bc_lr_cyc )  THEN
          IF ( bc_dirichlet_l  .OR.  bc_radiation_l )  THEN
             temp(:,:,nxl_mg(l)-1) = temp(:,:,nxl_mg(l))
          ENDIF
          IF ( bc_dirichlet_r  .OR.  bc_radiation_r )  THEN
             temp(:,:,nxr_mg(l)+1) = temp(:,:,nxr_mg(l))
          ENDIF
       ENDIF

       IF ( .NOT. bc_ns_cyc )  THEN
          IF ( bc_dirichlet_n  .OR.  bc_radiation_n )  THEN
             temp(:,nyn_mg(l)+1,:) = temp(:,nyn_mg(l),:)
          ENDIF
          IF ( bc_dirichlet_s  .OR.  bc_radiation_s )  THEN
             temp(:,nys_mg(l)-1,:) = temp(:,nys_mg(l),:)
          ENDIF
       ENDIF

!
!--    Bottom and top boundary conditions
       IF ( ibc_p_b == 1 )  THEN
!
!--       equivalent to temp(nzb,:,: ) = temp(nzb+1,:,:)
          temp(nzb,:,: ) = temp(ind_even_odd+1,:,:)
       ELSE
          temp(nzb,:,: ) = 0.0_wp
       ENDIF

       IF ( ibc_p_t == 1 )  THEN
!
!--       equivalent to temp(nzt_mg(l)+1,:,: ) = temp(nzt_mg(l),:,:)
          temp(nzt_mg(l)+1,:,: ) = temp(ind_even_odd,:,:)
       ELSE
          temp(nzt_mg(l)+1,:,: ) = 0.0_wp
       ENDIF

       CALL cpu_log( log_point_s(55), 'prolong', 'stop' )

    END SUBROUTINE prolong


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Relaxation method for the multigrid scheme. A Gauss-Seidel iteration with
!> 3D-Red-Black decomposition (GS-RB) is used.
!------------------------------------------------------------------------------!
    SUBROUTINE redblack( f_mg, p_mg )


       USE arrays_3d,                                                          &
           ONLY:  rho_air_mg

       USE control_parameters,                                                 &
           ONLY:  bc_lr_cyc, bc_ns_cyc, ibc_p_b, ibc_p_t, ngsrb 

       USE grid_variables,                                                     &
           ONLY:  ddx2_mg, ddy2_mg

       USE indices,                                                            &
           ONLY:  nxl_mg, nxr_mg, nys_mg, nyn_mg, nzb, nzt_mg

       IMPLICIT NONE

       INTEGER(iwp) :: color    !< grid point color, either red or black
       INTEGER(iwp) :: i        !< index variable along x
       INTEGER(iwp) :: ic       !< index variable along x
       INTEGER(iwp) :: j        !< index variable along y
       INTEGER(iwp) :: jc       !< index variable along y
       INTEGER(iwp) :: jj       !< index variable along y
       INTEGER(iwp) :: k        !< index variable along z
       INTEGER(iwp) :: l        !< grid level
       INTEGER(iwp) :: n        !< loop variable Gauß-Seidel iterations
       INTEGER(iwp) :: km1      !< index variable (k-1)
       INTEGER(iwp) :: kp1      !< index variable (k+1)

       LOGICAL      :: unroll   !< flag indicating whether loop unrolling is possible

       REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,                           &
                           nys_mg(grid_level)-1:nyn_mg(grid_level)+1,          &
                           nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) ::       &
                                      f_mg  !< residual of perturbation pressure 
       REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,                           &
                           nys_mg(grid_level)-1:nyn_mg(grid_level)+1,          &
                           nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) ::       &
                                      p_mg  !< perturbation pressure

       l = grid_level

       unroll = ( MOD( nyn_mg(l)-nys_mg(l)+1, 4 ) == 0  .AND.                  &
                  MOD( nxr_mg(l)-nxl_mg(l)+1, 2 ) == 0 )

       DO  n = 1, ngsrb
       
          DO  color = 1, 2

             IF ( .NOT. unroll )  THEN

                CALL cpu_log( log_point_s(36), 'redblack_no_unroll_f', 'start' )
!
!--             Without unrolling of loops, no cache optimization
                !$OMP PARALLEL PRIVATE (i,j,k,km1,kp1)
                !$OMP DO
                DO  i = nxl_mg(l), nxr_mg(l), 2
                   DO  j = nys_mg(l) + 2 - color, nyn_mg(l), 2
                      !DIR$ IVDEP
                      DO  k = ind_even_odd+1, nzt_mg(l)
                         km1 = k-ind_even_odd-1
                         kp1 = k-ind_even_odd
                         p_mg(k,j,i) = 1.0_wp / f1_mg_b(k,l) * (               &
                                 rho_air_mg(k,l) * ddx2_mg(l) *                &
                               ( p_mg(k,j,i+1) + p_mg(k,j,i-1) )               &
                               + rho_air_mg(k,l) * ddy2_mg(l) *                &
                               ( p_mg(k,j+1,i) + p_mg(k,j-1,i) )               &
                               + f2_mg_b(k,l) * p_mg(kp1,j,i)                  &
                               + f3_mg_b(k,l) * p_mg(km1,j,i)                  &
                               - f_mg(k,j,i)                   )
                      ENDDO
                   ENDDO
                ENDDO
   
                !$OMP DO
                DO  i = nxl_mg(l)+1, nxr_mg(l), 2
                   DO  j = nys_mg(l) + (color-1), nyn_mg(l), 2
                       !DIR$ IVDEP
                       DO  k = ind_even_odd+1, nzt_mg(l)
                         km1 = k-ind_even_odd-1
                         kp1 = k-ind_even_odd
                         p_mg(k,j,i) = 1.0_wp / f1_mg_b(k,l) * (               &
                                 rho_air_mg(k,l) * ddx2_mg(l) *                &
                               ( p_mg(k,j,i+1) + p_mg(k,j,i-1) )               &
                               + rho_air_mg(k,l) * ddy2_mg(l) *                &
                               ( p_mg(k,j+1,i) + p_mg(k,j-1,i) )               &
                               + f2_mg_b(k,l) * p_mg(kp1,j,i)                  &
                               + f3_mg_b(k,l) * p_mg(km1,j,i)                  &
                               - f_mg(k,j,i)                   )
                      ENDDO
                   ENDDO
                ENDDO
  
                !$OMP DO
                DO  i = nxl_mg(l), nxr_mg(l), 2
                   DO  j = nys_mg(l) + (color-1), nyn_mg(l), 2
                      !DIR$ IVDEP
                      DO  k = nzb+1, ind_even_odd
                         km1 = k+ind_even_odd
                         kp1 = k+ind_even_odd+1
                         p_mg(k,j,i) = 1.0_wp / f1_mg_b(k,l) * (               &
                                 rho_air_mg(k,l) * ddx2_mg(l) *                &
                               ( p_mg(k,j,i+1) + p_mg(k,j,i-1) )               &
                               + rho_air_mg(k,l) * ddy2_mg(l) *                &
                               ( p_mg(k,j+1,i) + p_mg(k,j-1,i) )               &
                               + f2_mg_b(k,l) * p_mg(kp1,j,i)                  &
                               + f3_mg_b(k,l) * p_mg(km1,j,i)                  &
                               - f_mg(k,j,i)                   )
                      ENDDO
                   ENDDO
                ENDDO

                !$OMP DO
                DO  i = nxl_mg(l)+1, nxr_mg(l), 2
                   DO  j = nys_mg(l) + 2 - color, nyn_mg(l), 2
                      !DIR$ IVDEP
                      DO  k = nzb+1, ind_even_odd
                         km1 = k+ind_even_odd
                         kp1 = k+ind_even_odd+1
                         p_mg(k,j,i) = 1.0_wp / f1_mg_b(k,l) * (               &
                                 rho_air_mg(k,l) * ddx2_mg(l) *                &
                               ( p_mg(k,j,i+1) + p_mg(k,j,i-1) )               &
                               + rho_air_mg(k,l) * ddy2_mg(l) *                &
                               ( p_mg(k,j+1,i) + p_mg(k,j-1,i) )               &
                               + f2_mg_b(k,l) * p_mg(kp1,j,i)                  &
                               + f3_mg_b(k,l) * p_mg(km1,j,i)                  &
                               - f_mg(k,j,i)                   )
                      ENDDO
                   ENDDO
                ENDDO
                !$OMP END PARALLEL

                CALL cpu_log( log_point_s(36), 'redblack_no_unroll_f', 'stop' )

             ELSE
!
!--              Loop unrolling along y, only one i loop for better cache use
                CALL cpu_log( log_point_s(38), 'redblack_unroll_f', 'start' )

                !$OMP PARALLEL PRIVATE (i,j,k,ic,jc,km1,kp1,jj)
                !$OMP DO
                DO  ic = nxl_mg(l), nxr_mg(l), 2
                   DO  jc = nys_mg(l), nyn_mg(l), 4
                      i  = ic
                      jj = jc+2-color
                      !DIR$ IVDEP
                      DO  k = ind_even_odd+1, nzt_mg(l)
                         km1 = k-ind_even_odd-1
                         kp1 = k-ind_even_odd
                         j   = jj
                         p_mg(k,j,i) = 1.0_wp / f1_mg_b(k,l) * (               &
                                 rho_air_mg(k,l) * ddx2_mg(l) *                &
                               ( p_mg(k,j,i+1) + p_mg(k,j,i-1) )               &
                               + rho_air_mg(k,l) * ddy2_mg(l) *                &
                               ( p_mg(k,j+1,i) + p_mg(k,j-1,i) )               &
                               + f2_mg_b(k,l) * p_mg(kp1,j,i)                  &
                               + f3_mg_b(k,l) * p_mg(km1,j,i)                  &
                               - f_mg(k,j,i)                   )
                         j = jj+2
                         p_mg(k,j,i) = 1.0_wp / f1_mg_b(k,l) * (               &
                                 rho_air_mg(k,l) * ddx2_mg(l) *                &
                               ( p_mg(k,j,i+1) + p_mg(k,j,i-1) )               &
                               + rho_air_mg(k,l) * ddy2_mg(l) *                &
                               ( p_mg(k,j+1,i) + p_mg(k,j-1,i) )               &
                               + f2_mg_b(k,l) * p_mg(kp1,j,i)                  &
                               + f3_mg_b(k,l) * p_mg(km1,j,i)                  &
                               - f_mg(k,j,i)                   )
                      ENDDO

                      i  = ic+1
                      jj = jc+color-1
                      !DIR$ IVDEP
                      DO  k = ind_even_odd+1, nzt_mg(l)
                         km1 = k-ind_even_odd-1
                         kp1 = k-ind_even_odd
                         j   = jj
                         p_mg(k,j,i) = 1.0_wp / f1_mg_b(k,l) * (               &
                                 rho_air_mg(k,l) * ddx2_mg(l) *                &
                               ( p_mg(k,j,i+1) + p_mg(k,j,i-1) )               &
                               + rho_air_mg(k,l) * ddy2_mg(l) *                &
                               ( p_mg(k,j+1,i) + p_mg(k,j-1,i) )               &
                               + f2_mg_b(k,l) * p_mg(kp1,j,i)                  &
                               + f3_mg_b(k,l) * p_mg(km1,j,i)                  &
                               - f_mg(k,j,i)                   )
                         j = jj+2
                         p_mg(k,j,i) = 1.0_wp / f1_mg_b(k,l) * (               &
                                 rho_air_mg(k,l) * ddx2_mg(l) *                &
                               ( p_mg(k,j,i+1) + p_mg(k,j,i-1) )               &
                               + rho_air_mg(k,l) * ddy2_mg(l) *                &
                               ( p_mg(k,j+1,i) + p_mg(k,j-1,i) )               &
                               + f2_mg_b(k,l) * p_mg(kp1,j,i)                  &
                               + f3_mg_b(k,l) * p_mg(km1,j,i)                  &
                               - f_mg(k,j,i)                   )
                      ENDDO

                      i  = ic
                      jj = jc+color-1
                      !DIR$ IVDEP
                      DO  k = nzb+1, ind_even_odd
                         km1 = k+ind_even_odd
                         kp1 = k+ind_even_odd+1
                         j   = jj
                         p_mg(k,j,i) = 1.0_wp / f1_mg_b(k,l) * (               &
                                 rho_air_mg(k,l) * ddx2_mg(l) *                &
                               ( p_mg(k,j,i+1) + p_mg(k,j,i-1) )               &
                               + rho_air_mg(k,l) * ddy2_mg(l) *                &
                               ( p_mg(k,j+1,i) + p_mg(k,j-1,i) )               &
                               + f2_mg_b(k,l) * p_mg(kp1,j,i)                  &
                               + f3_mg_b(k,l) * p_mg(km1,j,i)                  &
                               - f_mg(k,j,i)                   )
                         j = jj+2
                         p_mg(k,j,i) = 1.0_wp / f1_mg_b(k,l) * (               &
                                 rho_air_mg(k,l) * ddx2_mg(l) *                &
                               ( p_mg(k,j,i+1) + p_mg(k,j,i-1) )               &
                               + rho_air_mg(k,l) * ddy2_mg(l) *                &
                               ( p_mg(k,j+1,i) + p_mg(k,j-1,i) )               &
                               + f2_mg_b(k,l) * p_mg(kp1,j,i)                  &
                               + f3_mg_b(k,l) * p_mg(km1,j,i)                  &
                               - f_mg(k,j,i)                   )
                      ENDDO

                      i  = ic+1
                      jj = jc+2-color
                      !DIR$ IVDEP
                      DO  k = nzb+1, ind_even_odd
                         km1 = k+ind_even_odd
                         kp1 = k+ind_even_odd+1
                         j   = jj
                         p_mg(k,j,i) = 1.0_wp / f1_mg_b(k,l) * (               &
                                 rho_air_mg(k,l) * ddx2_mg(l) *                &
                               ( p_mg(k,j,i+1) + p_mg(k,j,i-1) )               &
                               + rho_air_mg(k,l) * ddy2_mg(l) *                &
                               ( p_mg(k,j+1,i) + p_mg(k,j-1,i) )               &
                               + f2_mg_b(k,l) * p_mg(kp1,j,i)                  &
                               + f3_mg_b(k,l) * p_mg(km1,j,i)                  &
                               - f_mg(k,j,i)                   )
                         j = jj+2
                         p_mg(k,j,i) = 1.0_wp / f1_mg_b(k,l) * (               &
                                 rho_air_mg(k,l) * ddx2_mg(l) *                &
                               ( p_mg(k,j,i+1) + p_mg(k,j,i-1) )               &
                               + rho_air_mg(k,l) * ddy2_mg(l) *                &
                               ( p_mg(k,j+1,i) + p_mg(k,j-1,i) )               &
                               + f2_mg_b(k,l) * p_mg(kp1,j,i)                  &
                               + f3_mg_b(k,l) * p_mg(km1,j,i)                  &
                               - f_mg(k,j,i)                   )
                      ENDDO

                   ENDDO
                ENDDO
                !$OMP END PARALLEL

                CALL cpu_log( log_point_s(38), 'redblack_unroll_f', 'stop' )

             ENDIF

!
!--          Horizontal boundary conditions
             CALL special_exchange_horiz( p_mg, color )

             IF ( .NOT. bc_lr_cyc )  THEN
                IF ( bc_dirichlet_l  .OR.  bc_radiation_l )  THEN
                   p_mg(:,:,nxl_mg(l)-1) = p_mg(:,:,nxl_mg(l))
                ENDIF
                IF ( bc_dirichlet_r  .OR.  bc_radiation_r )  THEN
                   p_mg(:,:,nxr_mg(l)+1) = p_mg(:,:,nxr_mg(l))
                ENDIF
             ENDIF

             IF ( .NOT. bc_ns_cyc )  THEN
                IF ( bc_dirichlet_n  .OR.  bc_radiation_n )  THEN
                   p_mg(:,nyn_mg(l)+1,:) = p_mg(:,nyn_mg(l),:)
                ENDIF
                IF ( bc_dirichlet_s  .OR.  bc_radiation_s )  THEN
                   p_mg(:,nys_mg(l)-1,:) = p_mg(:,nys_mg(l),:)
                ENDIF
             ENDIF

!
!--          Bottom and top boundary conditions
             IF ( ibc_p_b == 1 )  THEN
!
!--             equivalent to p_mg(nzb,:,: ) = p_mg(nzb+1,:,:) 
                p_mg(nzb,:,: ) = p_mg(ind_even_odd+1,:,:)
             ELSE
                p_mg(nzb,:,: ) = 0.0_wp
             ENDIF

             IF ( ibc_p_t == 1 )  THEN
!
!--             equivalent to p_mg(nzt_mg(l)+1,:,: ) = p_mg(nzt_mg(l),:,:)
                p_mg(nzt_mg(l)+1,:,: ) = p_mg(ind_even_odd,:,:)
             ELSE
                p_mg(nzt_mg(l)+1,:,: ) = 0.0_wp
             ENDIF

          ENDDO

       ENDDO

    END SUBROUTINE redblack


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Sort k-Dimension from sequential into blocks of even and odd.
!> This is required to vectorize the red-black subroutine.
!> Version for 3D-REAL arrays
!------------------------------------------------------------------------------!
    SUBROUTINE sort_k_to_even_odd_blocks( p_mg , glevel )


       USE indices,                                                            &
           ONLY:  nxl_mg, nxr_mg, nys_mg, nyn_mg, nzb, nzt_mg

       IMPLICIT NONE

       INTEGER(iwp), INTENT(IN) ::  glevel  !< grid level 

       REAL(wp), DIMENSION(nzb:nzt_mg(glevel)+1,                               &
                           nys_mg(glevel)-1:nyn_mg(glevel)+1,                  &
                           nxl_mg(glevel)-1:nxr_mg(glevel)+1) ::               &
                                      p_mg  !< array to be sorted
!
!--    Local variables
       INTEGER(iwp) :: i        !< index variable along x
       INTEGER(iwp) :: j        !< index variable along y
       INTEGER(iwp) :: k        !< index variable along z
       INTEGER(iwp) :: l        !< grid level
       INTEGER(iwp) :: ind      !< index variable along z
       REAL(wp), DIMENSION(nzb:nzt_mg(glevel)+1) ::  tmp  !< odd-even sorted temporary array


       CALL cpu_log( log_point_s(52), 'sort_k_to_even_odd', 'start' )

       l = glevel
       ind_even_odd = even_odd_level(l)

       !$OMP PARALLEL PRIVATE (i,j,k,ind,tmp)
       !$OMP DO
       DO  i = nxl_mg(l)-1, nxr_mg(l)+1
          DO  j = nys_mg(l)-1, nyn_mg(l)+1

!
!--          Sort the data with even k index
             ind = nzb-1
             DO  k = nzb, nzt_mg(l), 2
                ind = ind + 1
                tmp(ind) = p_mg(k,j,i)
             ENDDO
!
!--          Sort the data with odd k index
             DO  k = nzb+1, nzt_mg(l)+1, 2
                ind = ind + 1
                tmp(ind) = p_mg(k,j,i)
             ENDDO

             p_mg(:,j,i) = tmp

          ENDDO
       ENDDO
       !$OMP END PARALLEL

       CALL cpu_log( log_point_s(52), 'sort_k_to_even_odd', 'stop' )

    END SUBROUTINE sort_k_to_even_odd_blocks


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Sort k-Dimension from sequential into blocks of even and odd.
!> This is required to vectorize the red-black subroutine.
!> Version for 1D-REAL arrays
!------------------------------------------------------------------------------!
    SUBROUTINE sort_k_to_even_odd_blocks_1d( f_mg, f_mg_b, glevel )


       USE indices,                                                            &
           ONLY:  nzb, nzt_mg

       IMPLICIT NONE

       INTEGER(iwp), INTENT(IN) ::  glevel  !< grid level

       REAL(wp), DIMENSION(nzb+1:nzt_mg(glevel)) ::  f_mg    !< 1D input array
       REAL(wp), DIMENSION(nzb:nzt_mg(glevel)+1) ::  f_mg_b  !< 1D output array

!
!--    Local variables
       INTEGER(iwp) :: ind   !< index variable along z
       INTEGER(iwp) :: k     !< index variable along z 


       ind = nzb - 1
!
!--    Sort the data with even k index
       DO  k = nzb, nzt_mg(glevel), 2
          ind = ind + 1
          IF ( k >= nzb+1  .AND.  k <= nzt_mg(glevel) )  THEN
             f_mg_b(ind) = f_mg(k)
          ENDIF
       ENDDO
!
!--    Sort the data with odd k index
       DO  k = nzb+1, nzt_mg(glevel)+1, 2
          ind = ind + 1
          IF( k >= nzb+1  .AND.  k <= nzt_mg(glevel) )  THEN
             f_mg_b(ind) = f_mg(k)
          ENDIF
       ENDDO

    END SUBROUTINE sort_k_to_even_odd_blocks_1d


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Sort k-Dimension from sequential into blocks of even and odd.
!> This is required to vectorize the red-black subroutine.
!> Version for 2D-INTEGER arrays
!------------------------------------------------------------------------------!
!    SUBROUTINE sort_k_to_even_odd_blocks_int( i_mg , glevel )
!
!
!       USE indices,                                                            &
!           ONLY:  nxl_mg, nxr_mg, nys_mg, nyn_mg, nzb, nzt_mg
!
!       IMPLICIT NONE
!
!       INTEGER(iwp), INTENT(IN) ::  glevel  !< grid level
!
!       INTEGER(iwp), DIMENSION(nzb:nzt_mg(glevel)+1,                           &
!                               nys_mg(glevel)-1:nyn_mg(glevel)+1,              &
!                               nxl_mg(glevel)-1:nxr_mg(glevel)+1) ::           &
!                                    i_mg    !< array to be sorted
!
!--    Local variables
!       INTEGER(iwp) :: i        !< index variabel along x
!       INTEGER(iwp) :: j        !< index variable along y
!       INTEGER(iwp) :: k        !< index variable along z
!       INTEGER(iwp) :: l        !< grid level
!       INTEGER(iwp) :: ind      !< index variable along z
!       INTEGER(iwp),DIMENSION(nzb:nzt_mg(glevel)+1) ::  tmp  !< temporary odd-even sorted array
!
!
!       CALL cpu_log( log_point_s(52), 'sort_k_to_even_odd', 'start' )
!
!       l = glevel
!       ind_even_odd = even_odd_level(l)
!
!       DO  i = nxl_mg(l)-1, nxr_mg(l)+1
!          DO  j = nys_mg(l)-1, nyn_mg(l)+1
!
!
!--          Sort the data with even k index
!             ind = nzb-1
!             DO  k = nzb, nzt_mg(l), 2
!                ind = ind + 1
!                tmp(ind) = i_mg(k,j,i)
!             ENDDO
!
!--          Sort the data with odd k index
!             DO  k = nzb+1, nzt_mg(l)+1, 2
!                ind = ind + 1
!                tmp(ind) = i_mg(k,j,i)
!             ENDDO
!
!             i_mg(:,j,i) = tmp
!
!          ENDDO
!       ENDDO
!
!       CALL cpu_log( log_point_s(52), 'sort_k_to_even_odd', 'stop' )
!
!    END SUBROUTINE sort_k_to_even_odd_blocks_int


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Sort k-dimension from blocks of even and odd into sequential
!------------------------------------------------------------------------------!
    SUBROUTINE sort_k_to_sequential( p_mg )


       USE control_parameters,                                                 &
           ONLY:  grid_level

       USE indices,                                                            &
           ONLY:  nxl_mg, nxr_mg, nys_mg, nyn_mg, nzb, nzt_mg

       IMPLICIT NONE

       REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,                           &
                           nys_mg(grid_level)-1:nyn_mg(grid_level)+1,          &
                           nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) ::       &
                                     p_mg  !< array to be sorted
!
!--    Local variables
       INTEGER(iwp) :: i        !< index variable along x
       INTEGER(iwp) :: j        !< index variable along y
       INTEGER(iwp) :: k        !< index variable along z
       INTEGER(iwp) :: l        !< grid level
       INTEGER(iwp) :: ind      !< index variable along z

       REAL(wp),DIMENSION(nzb:nzt_mg(grid_level)+1) ::  tmp


       l = grid_level

       !$OMP PARALLEL PRIVATE (i,j,k,ind,tmp)
       !$OMP DO
       DO  i = nxl_mg(l)-1, nxr_mg(l)+1
          DO  j = nys_mg(l)-1, nyn_mg(l)+1

             ind = nzb - 1
             tmp = p_mg(:,j,i)
             DO  k = nzb, nzt_mg(l), 2
                ind = ind + 1
                p_mg(k,j,i) = tmp(ind)
             ENDDO

             DO  k = nzb+1, nzt_mg(l)+1, 2
                ind = ind + 1
                p_mg(k,j,i) = tmp(ind)
             ENDDO
          ENDDO
       ENDDO
       !$OMP END PARALLEL

    END SUBROUTINE sort_k_to_sequential


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Gather subdomain data from all PEs.
!------------------------------------------------------------------------------!
#if defined( __parallel )
    SUBROUTINE mg_gather( f2, f2_sub )

       USE control_parameters,                                                 &
           ONLY:  grid_level

       USE cpulog,                                                             &
           ONLY:  cpu_log, log_point_s

       USE indices,                                                            &
           ONLY:  mg_loc_ind, nxl_mg, nxr_mg, nys_mg, nyn_mg, nzb, nzt_mg

       IMPLICIT NONE

       INTEGER(iwp) ::  i       !<
       INTEGER(iwp) ::  il      !<
       INTEGER(iwp) ::  ir      !<
       INTEGER(iwp) ::  j       !<
       INTEGER(iwp) ::  jn      !<
       INTEGER(iwp) ::  js      !<
       INTEGER(iwp) ::  k       !<
       INTEGER(iwp) ::  nwords  !<

       REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,                           &
                       nys_mg(grid_level)-1:nyn_mg(grid_level)+1,              &
                       nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) ::  f2    !<
       REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,                           &
                       nys_mg(grid_level)-1:nyn_mg(grid_level)+1,              &
                       nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) ::  f2_l  !<

       REAL(wp), DIMENSION(nzb:mg_loc_ind(5,myid)+1,                           &
                           mg_loc_ind(3,myid)-1:mg_loc_ind(4,myid)+1,          &
                           mg_loc_ind(1,myid)-1:mg_loc_ind(2,myid)+1) ::  f2_sub  !<


       CALL cpu_log( log_point_s(34), 'mg_gather', 'start' )

       f2_l = 0.0_wp

!
!--    Store the local subdomain array on the total array
       js = mg_loc_ind(3,myid)
       IF ( south_border_pe )  js = js - 1
       jn = mg_loc_ind(4,myid)
       IF ( north_border_pe )  jn = jn + 1
       il = mg_loc_ind(1,myid)
       IF ( left_border_pe )   il = il - 1
       ir = mg_loc_ind(2,myid)
       IF ( right_border_pe )  ir = ir + 1
       DO  i = il, ir
          DO  j = js, jn
             DO  k = nzb, nzt_mg(grid_level)+1
                f2_l(k,j,i) = f2_sub(k,j,i)
             ENDDO
          ENDDO
       ENDDO

!
!--    Find out the number of array elements of the total array
       nwords = SIZE( f2 )

!
!--    Gather subdomain data from all PEs
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( f2_l(nzb,nys_mg(grid_level)-1,nxl_mg(grid_level)-1), &
                           f2(nzb,nys_mg(grid_level)-1,nxl_mg(grid_level)-1),   &
                           nwords, MPI_REAL, MPI_SUM, comm2d, ierr )

       CALL cpu_log( log_point_s(34), 'mg_gather', 'stop' )
    
    END SUBROUTINE mg_gather
#endif


!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo It might be possible to improve the speed of this routine by using
!>       non-blocking communication
!------------------------------------------------------------------------------!
#if defined( __parallel )
    SUBROUTINE mg_scatter( p2, p2_sub )

       USE control_parameters,                                                 &
           ONLY:  grid_level

       USE cpulog,                                                             &
           ONLY:  cpu_log, log_point_s

       USE indices,                                                            &
           ONLY:  mg_loc_ind, nxl_mg, nxr_mg, nys_mg, nyn_mg, nzb, nzt_mg

       IMPLICIT NONE

       REAL(wp), DIMENSION(nzb:nzt_mg(grid_level-1)+1,                         &
                           nys_mg(grid_level-1)-1:nyn_mg(grid_level-1)+1,      &
                           nxl_mg(grid_level-1)-1:nxr_mg(grid_level-1)+1) ::  p2  !<

       REAL(wp), DIMENSION(nzb:mg_loc_ind(5,myid)+1,                           &
                           mg_loc_ind(3,myid)-1:mg_loc_ind(4,myid)+1,          &
                           mg_loc_ind(1,myid)-1:mg_loc_ind(2,myid)+1) ::  p2_sub  !<


       CALL cpu_log( log_point_s(35), 'mg_scatter', 'start' )

       p2_sub = p2(:,mg_loc_ind(3,myid)-1:mg_loc_ind(4,myid)+1, &
                     mg_loc_ind(1,myid)-1:mg_loc_ind(2,myid)+1)

       CALL cpu_log( log_point_s(35), 'mg_scatter', 'stop' )
    
    END SUBROUTINE mg_scatter
#endif

!------------------------------------------------------------------------------!
! Description:
! ------------
!> This is where the multigrid technique takes place. V- and W- Cycle are
!> implemented and steered by the parameter "gamma". Parameter "nue" determines
!> the convergence of the multigrid iterative solution. There are nue times
!> RB-GS iterations. It should be set to "1" or "2", considering the time effort
!> one would like to invest. Last choice shows a very good converging factor,
!> but leads to an increase in computing time.
!------------------------------------------------------------------------------!
    RECURSIVE SUBROUTINE next_mg_level( f_mg, p_mg, p3, r )

       USE control_parameters,                                                 &
           ONLY:  bc_lr_dirrad, bc_lr_raddir, bc_ns_dirrad, bc_ns_raddir,      &
                  child_domain, gamma_mg, grid_level_count, maximum_grid_level,&
                  mg_switch_to_pe0_level, mg_switch_to_pe0, ngsrb

       USE indices,                                                            &
           ONLY:  mg_loc_ind, nxl, nxl_mg, nxr, nxr_mg, nys, nys_mg, nyn,      &
                  nyn_mg, nzb, nzt, nzt_mg

       IMPLICIT NONE

       INTEGER(iwp) ::  i            !< index variable along x
       INTEGER(iwp) ::  j            !< index variable along y
       INTEGER(iwp) ::  k            !< index variable along z
       INTEGER(iwp) ::  nxl_mg_save  !<
       INTEGER(iwp) ::  nxr_mg_save  !<
       INTEGER(iwp) ::  nyn_mg_save  !<
       INTEGER(iwp) ::  nys_mg_save  !<
       INTEGER(iwp) ::  nzt_mg_save  !<

       REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,                           &
                           nys_mg(grid_level)-1:nyn_mg(grid_level)+1,          &
                           nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) :: f_mg  !<
       REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,                           &
                           nys_mg(grid_level)-1:nyn_mg(grid_level)+1,          &
                           nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) :: p_mg  !<
       REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,                           &
                           nys_mg(grid_level)-1:nyn_mg(grid_level)+1,          &
                           nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) :: p3    !<
       REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,                           &
                           nys_mg(grid_level)-1:nyn_mg(grid_level)+1,          &
                           nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) :: r     !<

       REAL(wp), DIMENSION(nzb:nzt_mg(grid_level-1)+1,                         &
                           nys_mg(grid_level-1)-1:nyn_mg(grid_level-1)+1,      &
                           nxl_mg(grid_level-1)-1:nxr_mg(grid_level-1)+1) ::  f2  !<
       REAL(wp), DIMENSION(nzb:nzt_mg(grid_level-1)+1,                         &
                           nys_mg(grid_level-1)-1:nyn_mg(grid_level-1)+1,      &
                           nxl_mg(grid_level-1)-1:nxr_mg(grid_level-1)+1) ::  p2  !<

       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  f2_sub  !<

#if defined( __parallel )
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  p2_sub  !<
#endif

!
!--    Restriction to the coarsest grid
    10 IF ( grid_level == 1 )  THEN

!
!--       Solution on the coarsest grid. Double the number of Gauss-Seidel
!--       iterations in order to get a more accurate solution.
          ngsrb = 2 * ngsrb

          ind_even_odd = even_odd_level(grid_level)

          CALL redblack( f_mg, p_mg )

          ngsrb = ngsrb / 2


       ELSEIF ( grid_level /= 1 )  THEN

          grid_level_count(grid_level) = grid_level_count(grid_level) + 1

!
!--       Solution on the actual grid level
          ind_even_odd = even_odd_level(grid_level)

          CALL redblack( f_mg, p_mg )

!
!--       Determination of the actual residual
          CALL resid( f_mg, p_mg, r )

!--       Restriction of the residual (finer grid values!) to the next coarser
!--       grid. Therefore, the grid level has to be decremented now. nxl..nzt have
!--       to be set to the coarse grid values, because these variables are needed
!--       for the exchange of ghost points in routine exchange_horiz
          grid_level = grid_level - 1

          nxl = nxl_mg(grid_level)
          nys = nys_mg(grid_level)
          nxr = nxr_mg(grid_level)
          nyn = nyn_mg(grid_level)
          nzt = nzt_mg(grid_level)

          IF ( grid_level == mg_switch_to_pe0_level )  THEN

!
!--          From this level on, calculations are done on PE0 only.
!--          First, carry out restriction on the subdomain.
!--          Therefore, indices of the level have to be changed to subdomain
!--          values in between (otherwise, the restrict routine would expect
!--          the gathered array)

             nxl_mg_save = nxl_mg(grid_level)
             nxr_mg_save = nxr_mg(grid_level)
             nys_mg_save = nys_mg(grid_level)
             nyn_mg_save = nyn_mg(grid_level)
             nzt_mg_save = nzt_mg(grid_level)
             nxl_mg(grid_level) = mg_loc_ind(1,myid)
             nxr_mg(grid_level) = mg_loc_ind(2,myid)
             nys_mg(grid_level) = mg_loc_ind(3,myid)
             nyn_mg(grid_level) = mg_loc_ind(4,myid)
             nzt_mg(grid_level) = mg_loc_ind(5,myid)
             nxl = mg_loc_ind(1,myid)
             nxr = mg_loc_ind(2,myid)
             nys = mg_loc_ind(3,myid)
             nyn = mg_loc_ind(4,myid)
             nzt = mg_loc_ind(5,myid)

             ALLOCATE( f2_sub(nzb:nzt_mg(grid_level)+1,                    &
                              nys_mg(grid_level)-1:nyn_mg(grid_level)+1,   &
                              nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) )

             CALL restrict( f2_sub, r )

!
!--          Restore the correct indices of this level
             nxl_mg(grid_level) = nxl_mg_save
             nxr_mg(grid_level) = nxr_mg_save
             nys_mg(grid_level) = nys_mg_save
             nyn_mg(grid_level) = nyn_mg_save
             nzt_mg(grid_level) = nzt_mg_save
             nxl = nxl_mg(grid_level)
             nxr = nxr_mg(grid_level)
             nys = nys_mg(grid_level)
             nyn = nyn_mg(grid_level)
             nzt = nzt_mg(grid_level)
!
!--          Gather all arrays from the subdomains on PE0
#if defined( __parallel )
             CALL mg_gather( f2, f2_sub )
#endif

!
!--          Set switch for routine exchange_horiz, that no ghostpoint exchange
!--          has to be carried out from now on
             mg_switch_to_pe0 = .TRUE.

!
!--          In case of non-cyclic lateral boundary conditions, both in- and
!--          outflow conditions have to be used on all PEs after the switch,
!--          because then they have the total domain.
             IF ( bc_lr_dirrad )  THEN
                bc_dirichlet_l  = .TRUE.
                bc_dirichlet_r  = .FALSE.
                bc_radiation_l = .FALSE.
                bc_radiation_r = .TRUE.
             ELSEIF ( bc_lr_raddir )  THEN
                bc_dirichlet_l  = .FALSE.
                bc_dirichlet_r  = .TRUE.
                bc_radiation_l = .TRUE.
                bc_radiation_r = .FALSE.
             ELSEIF ( child_domain  .OR.  nesting_offline )  THEN
                bc_dirichlet_l = .TRUE.
                bc_dirichlet_r = .TRUE.
             ENDIF

             IF ( bc_ns_dirrad )  THEN
                bc_dirichlet_n  = .TRUE.
                bc_dirichlet_s  = .FALSE.
                bc_radiation_n = .FALSE.
                bc_radiation_s = .TRUE.
             ELSEIF ( bc_ns_raddir )  THEN
                bc_dirichlet_n  = .FALSE.
                bc_dirichlet_s  = .TRUE.
                bc_radiation_n = .TRUE.
                bc_radiation_s = .FALSE.
             ELSEIF ( child_domain  .OR.  nesting_offline)  THEN
                bc_dirichlet_s = .TRUE.
                bc_dirichlet_n = .TRUE.
             ENDIF

             DEALLOCATE( f2_sub )

          ELSE

             CALL restrict( f2, r )

             ind_even_odd = even_odd_level(grid_level)  ! must be after restrict

          ENDIF

          p2 = 0.0_wp

!
!--       Repeat the same procedure till the coarsest grid is reached
          CALL next_mg_level( f2, p2, p3, r )

       ENDIF

!
!--    Now follows the prolongation
       IF ( grid_level >= 2 )  THEN

!
!--       Prolongation of the new residual. The values are transferred
!--       from the coarse to the next finer grid.
          IF ( grid_level == mg_switch_to_pe0_level+1 )  THEN

#if defined( __parallel )
!
!--          At this level, the new residual first has to be scattered from
!--          PE0 to the other PEs
             ALLOCATE( p2_sub(nzb:mg_loc_ind(5,myid)+1,             &
                       mg_loc_ind(3,myid)-1:mg_loc_ind(4,myid)+1,   &
                       mg_loc_ind(1,myid)-1:mg_loc_ind(2,myid)+1) )

             CALL mg_scatter( p2, p2_sub )

!
!--          Therefore, indices of the previous level have to be changed to
!--          subdomain values in between (otherwise, the prolong routine would
!--          expect the gathered array)
             nxl_mg_save = nxl_mg(grid_level-1)
             nxr_mg_save = nxr_mg(grid_level-1)
             nys_mg_save = nys_mg(grid_level-1)
             nyn_mg_save = nyn_mg(grid_level-1)
             nzt_mg_save = nzt_mg(grid_level-1)
             nxl_mg(grid_level-1) = mg_loc_ind(1,myid)
             nxr_mg(grid_level-1) = mg_loc_ind(2,myid)
             nys_mg(grid_level-1) = mg_loc_ind(3,myid)
             nyn_mg(grid_level-1) = mg_loc_ind(4,myid)
             nzt_mg(grid_level-1) = mg_loc_ind(5,myid)

!
!--          Set switch for routine exchange_horiz, that ghostpoint exchange
!--          has to be carried again out from now on
             mg_switch_to_pe0 = .FALSE.

!
!--          For non-cyclic lateral boundary conditions and in case of nesting, 
!--          restore the in-/outflow conditions. 
             bc_dirichlet_l = .FALSE.;  bc_dirichlet_r = .FALSE.
             bc_dirichlet_n = .FALSE.;  bc_dirichlet_s = .FALSE.
             bc_radiation_l = .FALSE.;  bc_radiation_r = .FALSE.
             bc_radiation_n = .FALSE.;  bc_radiation_s = .FALSE.

             IF ( pleft == MPI_PROC_NULL )  THEN
                IF ( bc_lr_dirrad  .OR.  child_domain  .OR.  nesting_offline ) &
                THEN
                   bc_dirichlet_l = .TRUE.
                ELSEIF ( bc_lr_raddir )  THEN
                   bc_radiation_l = .TRUE.
                ENDIF
             ENDIF

             IF ( pright == MPI_PROC_NULL )  THEN
                IF ( bc_lr_dirrad )  THEN
                   bc_radiation_r = .TRUE.
                ELSEIF ( bc_lr_raddir  .OR.  child_domain  .OR.                &
                         nesting_offline )  THEN
                   bc_dirichlet_r = .TRUE.
                ENDIF
             ENDIF

             IF ( psouth == MPI_PROC_NULL )  THEN
                IF ( bc_ns_dirrad )  THEN
                   bc_radiation_s = .TRUE.
                ELSEIF ( bc_ns_raddir  .OR.  child_domain  .OR.                &
                         nesting_offline )  THEN
                   bc_dirichlet_s = .TRUE.
                ENDIF
             ENDIF

             IF ( pnorth == MPI_PROC_NULL )  THEN
                IF ( bc_ns_dirrad  .OR.  child_domain  .OR.  nesting_offline ) &
                THEN
                   bc_dirichlet_n = .TRUE.
                ELSEIF ( bc_ns_raddir )  THEN
                   bc_radiation_n = .TRUE.
                ENDIF
             ENDIF

             CALL prolong( p2_sub, p3 )

!
!--          Restore the correct indices of the previous level
             nxl_mg(grid_level-1) = nxl_mg_save
             nxr_mg(grid_level-1) = nxr_mg_save
             nys_mg(grid_level-1) = nys_mg_save
             nyn_mg(grid_level-1) = nyn_mg_save
             nzt_mg(grid_level-1) = nzt_mg_save

             DEALLOCATE( p2_sub )
#endif

          ELSE

             CALL prolong( p2, p3 )

          ENDIF

!
!--       Computation of the new pressure correction. Therefore,
!--       values from prior grids are added up automatically stage by stage.
          DO  i = nxl_mg(grid_level)-1, nxr_mg(grid_level)+1
             DO  j = nys_mg(grid_level)-1, nyn_mg(grid_level)+1
                DO  k = nzb, nzt_mg(grid_level)+1
                   p_mg(k,j,i) = p_mg(k,j,i) + p3(k,j,i)
                ENDDO
             ENDDO
          ENDDO

!
!--       Relaxation of the new solution
          CALL redblack( f_mg, p_mg )

       ENDIF


!
!--    The following few lines serve the steering of the multigrid scheme
       IF ( grid_level == maximum_grid_level )  THEN

          GOTO 20

       ELSEIF ( grid_level /= maximum_grid_level  .AND.  grid_level /= 1  .AND. &
                grid_level_count(grid_level) /= gamma_mg )  THEN

          GOTO 10

       ENDIF

!
!--    Reset counter for the next call of poismg
       grid_level_count(grid_level) = 0

!
!--    Continue with the next finer level. nxl..nzt have to be
!--    set to the finer grid values, because these variables are needed for the
!--    exchange of ghost points in routine exchange_horiz
       grid_level = grid_level + 1
       ind_even_odd = even_odd_level(grid_level)

       nxl = nxl_mg(grid_level)
       nxr = nxr_mg(grid_level)
       nys = nys_mg(grid_level)
       nyn = nyn_mg(grid_level)
       nzt = nzt_mg(grid_level)

    20 CONTINUE

    END SUBROUTINE next_mg_level


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initial settings for sorting k-dimension from sequential order (alternate
!> even/odd) into blocks of even and odd or vice versa
!------------------------------------------------------------------------------!
    SUBROUTINE init_even_odd_blocks


       USE arrays_3d,                                                          &
           ONLY:  f1_mg, f2_mg, f3_mg

       USE control_parameters,                                                 &
           ONLY:  grid_level, maximum_grid_level

       USE indices,                                                            &
           ONLY:  nzb, nzt, nzt_mg

       USE indices,                                                            &
           ONLY:  nzb, nzt_mg

       IMPLICIT NONE
!
!--    Local variables
       INTEGER(iwp) ::  i     !<  
       INTEGER(iwp) ::  l     !<

       LOGICAL, SAVE ::  lfirst = .TRUE.


       IF ( .NOT. lfirst )  RETURN

       ALLOCATE( even_odd_level(maximum_grid_level) )

       ALLOCATE( f1_mg_b(nzb:nzt+1,maximum_grid_level),                        &
                 f2_mg_b(nzb:nzt+1,maximum_grid_level),                        &
                 f3_mg_b(nzb:nzt+1,maximum_grid_level) )

!
!--    Set border index between the even and odd block
       DO  i = maximum_grid_level, 1, -1
          even_odd_level(i) = nzt_mg(i) / 2
       ENDDO

!
!--    Sort grid coefficients used in red/black scheme and for calculating the
!--    residual to block (even/odd) structure
       DO  l = maximum_grid_level, 1 , -1
          CALL sort_k_to_even_odd_blocks( f1_mg(nzb+1:nzt_mg(grid_level),l),   &
                                          f1_mg_b(nzb:nzt_mg(grid_level)+1,l), &
                                          l )
          CALL sort_k_to_even_odd_blocks( f2_mg(nzb+1:nzt_mg(grid_level),l),   &
                                          f2_mg_b(nzb:nzt_mg(grid_level)+1,l), &
                                          l )
          CALL sort_k_to_even_odd_blocks( f3_mg(nzb+1:nzt_mg(grid_level),l),   &
                                          f3_mg_b(nzb:nzt_mg(grid_level)+1,l), &
                                          l )
       ENDDO

       lfirst = .FALSE.

     END SUBROUTINE init_even_odd_blocks


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Special exchange_horiz subroutine for use in redblack. Transfers only
!> "red" or "black" data points.
!------------------------------------------------------------------------------!
     SUBROUTINE special_exchange_horiz( p_mg, color )


       USE control_parameters,                                                 &
           ONLY:  grid_level

#if defined( __parallel )
       USE control_parameters,                                                 &
           ONLY:  mg_switch_to_pe0_level, synchronous_exchange
#endif

       USE indices,                                                            &
           ONLY:  nxl_mg, nxr_mg, nys_mg, nyn_mg, nzb, nzt_mg

#if defined( __parallel )
       USE indices,                                                            &
           ONLY:  nxl, nxr, nys, nyn, nzt
#endif

       IMPLICIT NONE

       REAL(wp), DIMENSION(nzb:nzt_mg(grid_level)+1,                           &
                           nys_mg(grid_level)-1:nyn_mg(grid_level)+1,          &
                           nxl_mg(grid_level)-1:nxr_mg(grid_level)+1) ::       &
                                    p_mg   !< treated array

       INTEGER(iwp), INTENT(IN) ::  color  !< flag for grid point type (red or black)

#if defined ( __parallel )
!
!--    Local variables
       INTEGER(iwp) ::  i        !< index variable along x
       INTEGER(iwp) ::  i1       !< index variable along x on coarse level 
       INTEGER(iwp) ::  i2       !< index variable along x on coarse level 

       INTEGER(iwp) ::  j        !< index variable along y
       INTEGER(iwp) ::  j1       !< index variable along y on coarse level 
       INTEGER(iwp) ::  j2       !< index variable along y on coarse level 
       INTEGER(iwp) ::  k        !< index variable along z
       INTEGER(iwp) ::  l        !< short for grid level
       INTEGER(iwp) ::  jys      !< index for lower local PE boundary along y
       INTEGER(iwp) ::  jyn      !< index for upper local PE boundary along y
       INTEGER(iwp) ::  ixl      !< index for lower local PE boundary along x
       INTEGER(iwp) ::  ixr      !< index for upper local PE boundary along x 

       LOGICAL      ::  synchronous_exchange_save    !< dummy to reset synchronous_exchange to prescribed value

       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  temp  !< temporary array on next coarser grid level

       synchronous_exchange_save   = synchronous_exchange
       synchronous_exchange        = .FALSE.

       l = grid_level

       ind_even_odd = even_odd_level(grid_level)

!
!--    Restricted transfer only on finer levels with enough data.
!--    Restricted transfer is not possible for levels smaller or equal to 
!--    'switch to PE0 levels', since array bounds does not fit. Moreover, 
!--    it is not possible for the coarsest grid level, since the dimensions
!--    of temp are not defined. For such cases, normal exchange_horiz is called.
       IF ( l > 1 .AND. l > mg_switch_to_pe0_level + 1 .AND.                   &
            ( ngp_xz(grid_level) >= 900 .OR. ngp_yz(grid_level) >= 900 ) )  THEN

          jys = nys_mg(grid_level-1)
          jyn = nyn_mg(grid_level-1)
          ixl = nxl_mg(grid_level-1)
          ixr = nxr_mg(grid_level-1)
          ALLOCATE( temp(nzb:nzt_mg(l-1)+1,jys-1:jyn+1,ixl-1:ixr+1) )
!
!--       Handling the even k Values
!--       Collecting data for the north - south exchange
!--       Since only every second value has to be transfered, data are stored
!--       on the next coarser grid level, because the arrays on that level
!--       have just the required size
          i1 = nxl_mg(grid_level-1)
          i2 = nxl_mg(grid_level-1)

          DO  i = nxl_mg(l), nxr_mg(l), 2
             DO  j = nys_mg(l) + 2 - color, nyn_mg(l), 2

                IF ( j == nys_mg(l) )  THEN
                   !DIR$ IVDEP
                   DO  k = ind_even_odd+1, nzt_mg(l)
                      temp(k-ind_even_odd,jys,i1) = p_mg(k,j,i)
                   ENDDO
                   i1 = i1 + 1

                ENDIF

                IF ( j == nyn_mg(l) )  THEN
                   !DIR$ IVDEP
                   DO  k = ind_even_odd+1, nzt_mg(l)
                      temp(k-ind_even_odd,jyn,i2) = p_mg(k,j,i)
                   ENDDO
                   i2 = i2 + 1

                ENDIF

             ENDDO
          ENDDO

          DO  i = nxl_mg(l)+1, nxr_mg(l), 2
             DO  j = nys_mg(l) + (color-1), nyn_mg(l), 2

                IF ( j == nys_mg(l) )  THEN
                   !DIR$ IVDEP
                   DO  k = ind_even_odd+1, nzt_mg(l)
                      temp(k-ind_even_odd,jys,i1) = p_mg(k,j,i)
                   ENDDO
                   i1 = i1 + 1

                ENDIF

                IF ( j == nyn_mg(l) )  THEN
                   !DIR$ IVDEP
                   DO  k = ind_even_odd+1, nzt_mg(l)
                      temp(k-ind_even_odd,jyn,i2) = p_mg(k,j,i)
                   ENDDO
                   i2 = i2 + 1

                ENDIF

             ENDDO
          ENDDO

          grid_level = grid_level-1

          nxl = nxl_mg(grid_level)
          nys = nys_mg(grid_level)
          nxr = nxr_mg(grid_level)
          nyn = nyn_mg(grid_level)
          nzt = nzt_mg(grid_level)

          send_receive = 'ns'
          CALL exchange_horiz( temp, 1 )

          grid_level = grid_level+1

          i1 = nxl_mg(grid_level-1)
          i2 = nxl_mg(grid_level-1)

          DO  i = nxl_mg(l), nxr_mg(l), 2
             DO  j = nys_mg(l) + 2 - color, nyn_mg(l), 2

                IF ( j == nys_mg(l) )  THEN
                   !DIR$ IVDEP
                   DO  k = ind_even_odd+1, nzt_mg(l)
                      p_mg(k,nyn_mg(l)+1,i) = temp(k-ind_even_odd,jyn+1,i1)
                   ENDDO
                   i1 = i1 + 1

                ENDIF

                IF ( j == nyn_mg(l) )  THEN
                   !DIR$ IVDEP
                   DO  k = ind_even_odd+1, nzt_mg(l)
                      p_mg(k,nys_mg(l)-1,i) = temp(k-ind_even_odd,jys-1,i2)
                   ENDDO
                   i2 = i2 + 1

                ENDIF

             ENDDO
          ENDDO

          DO  i = nxl_mg(l)+1, nxr_mg(l), 2
             DO  j = nys_mg(l) + (color-1), nyn_mg(l), 2

                IF ( j == nys_mg(l) )  THEN
                   !DIR$ IVDEP
                   DO  k = ind_even_odd+1, nzt_mg(l)
                      p_mg(k,nyn_mg(l)+1,i) = temp(k-ind_even_odd,jyn+1,i1)
                   ENDDO
                   i1 = i1 + 1

                ENDIF

                IF ( j == nyn_mg(l) )  THEN
                   !DIR$ IVDEP
                   DO  k = ind_even_odd+1, nzt_mg(l)
                      p_mg(k,nys_mg(l)-1,i) = temp(k-ind_even_odd,jys-1,i2)
                   ENDDO
                   i2 = i2 + 1

                ENDIF

             ENDDO
          ENDDO

!
!--       Collecting data for the left - right exchange
!--       Since only every second value has to be transfered, data are stored
!--       on the next coarser grid level, because the arrays on that level
!--       have just the required size
          j1 = nys_mg(grid_level-1)
          j2 = nys_mg(grid_level-1)

          DO j = nys_mg(l) + 2 - color, nyn_mg(l), 2
             DO  i = nxl_mg(l), nxr_mg(l), 2

                IF ( i == nxl_mg(l) )  THEN
                   !DIR$ IVDEP
                   DO  k = ind_even_odd+1, nzt_mg(l)
                      temp(k-ind_even_odd,j1,ixl) = p_mg(k,j,i)
                   ENDDO
                   j1 = j1 + 1

                ENDIF

                IF ( i == nxr_mg(l) )  THEN
                   !DIR$ IVDEP
                   DO  k = ind_even_odd+1, nzt_mg(l)
                      temp(k-ind_even_odd,j2,ixr) = p_mg(k,j,i)
                   ENDDO
                   j2 = j2 + 1

                ENDIF

             ENDDO
          ENDDO

          DO j = nys_mg(l) + (color-1), nyn_mg(l), 2
             DO  i = nxl_mg(l)+1, nxr_mg(l), 2

                IF ( i == nxl_mg(l) )  THEN
                   !DIR$ IVDEP
                   DO  k = ind_even_odd+1, nzt_mg(l)
                      temp(k-ind_even_odd,j1,ixl) = p_mg(k,j,i)
                   ENDDO
                   j1 = j1 + 1

                ENDIF

                IF ( i == nxr_mg(l) )  THEN
                   !DIR$ IVDEP
                   DO  k = ind_even_odd+1, nzt_mg(l)
                      temp(k-ind_even_odd,j2,ixr) = p_mg(k,j,i)
                   ENDDO
                   j2 = j2 + 1

                ENDIF

             ENDDO
          ENDDO

          grid_level = grid_level-1
          send_receive = 'lr'

          CALL exchange_horiz( temp, 1 )

          grid_level = grid_level+1

          j1 = nys_mg(grid_level-1)
          j2 = nys_mg(grid_level-1)

          DO j = nys_mg(l) + 2 - color, nyn_mg(l), 2
             DO  i = nxl_mg(l), nxr_mg(l), 2

                IF ( i == nxl_mg(l) )  THEN
                   !DIR$ IVDEP
                   DO  k = ind_even_odd+1, nzt_mg(l)
                      p_mg(k,j,nxr_mg(l)+1)  = temp(k-ind_even_odd,j1,ixr+1)
                   ENDDO
                   j1 = j1 + 1

                ENDIF

                IF ( i == nxr_mg(l) )  THEN
                   !DIR$ IVDEP
                   DO  k = ind_even_odd+1, nzt_mg(l)
                      p_mg(k,j,nxl_mg(l)-1)  = temp(k-ind_even_odd,j2,ixl-1)
                   ENDDO
                   j2 = j2 + 1

                ENDIF

             ENDDO
          ENDDO

          DO j = nys_mg(l) + (color-1), nyn_mg(l), 2
             DO  i = nxl_mg(l)+1, nxr_mg(l), 2

                IF ( i == nxl_mg(l) )  THEN
                   !DIR$ IVDEP
                   DO  k = ind_even_odd+1, nzt_mg(l)
                      p_mg(k,j,nxr_mg(l)+1)  = temp(k-ind_even_odd,j1,ixr+1)
                   ENDDO
                   j1 = j1 + 1

                ENDIF

                IF ( i == nxr_mg(l) )  THEN
                   !DIR$ IVDEP
                   DO  k = ind_even_odd+1, nzt_mg(l)
                      p_mg(k,j,nxl_mg(l)-1)  = temp(k-ind_even_odd,j2,ixl-1)
                   ENDDO
                   j2 = j2 + 1

                ENDIF

             ENDDO
          ENDDO

!
!--       Now handling the even k values
!--       Collecting data for the north - south exchange
!--       Since only every second value has to be transfered, data are stored
!--       on the next coarser grid level, because the arrays on that level
!--       have just the required size
          i1 = nxl_mg(grid_level-1)
          i2 = nxl_mg(grid_level-1)

          DO  i = nxl_mg(l), nxr_mg(l), 2
             DO  j = nys_mg(l) + (color-1), nyn_mg(l), 2

                IF ( j == nys_mg(l) )  THEN
                   !DIR$ IVDEP
                   DO  k = nzb+1, ind_even_odd
                      temp(k,jys,i1) = p_mg(k,j,i)
                   ENDDO
                   i1 = i1 + 1

                ENDIF

                IF ( j == nyn_mg(l) )  THEN
                   !DIR$ IVDEP
                   DO  k = nzb+1, ind_even_odd
                      temp(k,jyn,i2) = p_mg(k,j,i)
                   ENDDO
                   i2 = i2 + 1

                ENDIF

             ENDDO
          ENDDO

          DO  i = nxl_mg(l)+1, nxr_mg(l), 2
             DO j = nys_mg(l) + 2 - color, nyn_mg(l), 2

                IF ( j == nys_mg(l) )  THEN
                   !DIR$ IVDEP
                   DO  k = nzb+1, ind_even_odd
                      temp(k,jys,i1) = p_mg(k,j,i)
                   ENDDO
                   i1 = i1 + 1

                ENDIF

                IF ( j == nyn_mg(l) )  THEN
                   !DIR$ IVDEP
                   DO  k = nzb+1, ind_even_odd
                      temp(k,jyn,i2) = p_mg(k,j,i)
                   ENDDO
                   i2 = i2 + 1

                ENDIF

             ENDDO
          ENDDO

          grid_level = grid_level-1

          send_receive = 'ns'
          CALL exchange_horiz( temp, 1 )

          grid_level = grid_level+1

          i1 = nxl_mg(grid_level-1)
          i2 = nxl_mg(grid_level-1)

          DO  i = nxl_mg(l), nxr_mg(l), 2
             DO  j = nys_mg(l) + (color-1), nyn_mg(l), 2

                IF ( j == nys_mg(l) )  THEN
                   !DIR$ IVDEP
                   DO  k = nzb+1, ind_even_odd
                      p_mg(k,nyn_mg(l)+1,i) = temp(k,jyn+1,i1)
                   ENDDO
                   i1 = i1 + 1

                ENDIF

                IF ( j == nyn_mg(l) )  THEN
                   !DIR$ IVDEP
                   DO  k = nzb+1, ind_even_odd
                      p_mg(k,nys_mg(l)-1,i) = temp(k,jys-1,i2)
                   ENDDO
                   i2 = i2 + 1

                ENDIF

             ENDDO
          ENDDO

          DO  i = nxl_mg(l)+1, nxr_mg(l), 2
             DO j = nys_mg(l) + 2 - color, nyn_mg(l), 2

                IF ( j == nys_mg(l) )  THEN
                   !DIR$ IVDEP
                   DO  k = nzb+1, ind_even_odd
                      p_mg(k,nyn_mg(l)+1,i) = temp(k,jyn+1,i1)
                   ENDDO
                   i1 = i1 + 1

                ENDIF

                IF ( j == nyn_mg(l) )  THEN
                   !DIR$ IVDEP
                   DO  k = nzb+1, ind_even_odd
                      p_mg(k,nys_mg(l)-1,i) = temp(k,jys-1,i2)
                   ENDDO
                   i2 = i2 + 1

                ENDIF

             ENDDO
          ENDDO

          j1 = nys_mg(grid_level-1)
          j2 = nys_mg(grid_level-1)

          DO  i = nxl_mg(l), nxr_mg(l), 2
             DO  j = nys_mg(l) + (color-1), nyn_mg(l), 2

                IF ( i == nxl_mg(l) )  THEN
                   !DIR$ IVDEP
                   DO  k = nzb+1, ind_even_odd
                      temp(k,j1,ixl) = p_mg(k,j,i)
                   ENDDO
                   j1 = j1 + 1

                ENDIF

                IF ( i == nxr_mg(l) )  THEN
                   !DIR$ IVDEP
                   DO  k = nzb+1, ind_even_odd
                      temp(k,j2,ixr) = p_mg(k,j,i)
                   ENDDO
                   j2 = j2 + 1

                ENDIF

             ENDDO
          ENDDO

          DO  i = nxl_mg(l)+1, nxr_mg(l), 2
             DO j = nys_mg(l) + 2 - color, nyn_mg(l), 2

                IF ( i == nxl_mg(l) )  THEN
                   !DIR$ IVDEP
                   DO  k = nzb+1, ind_even_odd
                      temp(k,j1,ixl) = p_mg(k,j,i)
                   ENDDO
                   j1 = j1 + 1

                ENDIF

                IF ( i == nxr_mg(l) )  THEN
                   !DIR$ IVDEP
                   DO  k = nzb+1, ind_even_odd
                      temp(k,j2,ixr) = p_mg(k,j,i)
                   ENDDO
                   j2 = j2 + 1

                ENDIF

             ENDDO
          ENDDO

          grid_level = grid_level-1

          send_receive = 'lr'
          CALL exchange_horiz( temp, 1 )

          grid_level = grid_level+1

          nxl = nxl_mg(grid_level)
          nys = nys_mg(grid_level)
          nxr = nxr_mg(grid_level)
          nyn = nyn_mg(grid_level)
          nzt = nzt_mg(grid_level)

          j1 = nys_mg(grid_level-1)
          j2 = nys_mg(grid_level-1)

          DO  i = nxl_mg(l), nxr_mg(l), 2
             DO  j = nys_mg(l) + (color-1), nyn_mg(l), 2

                IF ( i == nxl_mg(l) )  THEN
                   !DIR$ IVDEP
                   DO  k = nzb+1, ind_even_odd
                      p_mg(k,j,nxr_mg(l)+1)  = temp(k,j1,ixr+1)
                   ENDDO
                   j1 = j1 + 1

                ENDIF

                IF ( i == nxr_mg(l) )  THEN
                   !DIR$ IVDEP
                   DO  k = nzb+1, ind_even_odd
                      p_mg(k,j,nxl_mg(l)-1)  = temp(k,j2,ixl-1)
                   ENDDO
                   j2 = j2 + 1

                ENDIF

             ENDDO
          ENDDO

          DO  i = nxl_mg(l)+1, nxr_mg(l), 2
             DO j = nys_mg(l) + 2 - color, nyn_mg(l), 2

                IF ( i == nxl_mg(l) )  THEN
                   !DIR$ IVDEP
                   DO  k = nzb+1, ind_even_odd
                      p_mg(k,j,nxr_mg(l)+1)  = temp(k,j1,ixr+1)
                   ENDDO
                   j1 = j1 + 1

                ENDIF

                IF ( i == nxr_mg(l) )  THEN
                   !DIR$ IVDEP
                   DO  k = nzb+1, ind_even_odd
                      p_mg(k,j,nxl_mg(l)-1)  = temp(k,j2,ixl-1)
                   ENDDO
                   j2 = j2 + 1

                ENDIF

             ENDDO
          ENDDO

          DEALLOCATE( temp )

       ELSE

!
!--       Standard horizontal ghost boundary exchange for small coarse grid
!--       levels, where the transfer time is latency bound
          CALL exchange_horiz( p_mg, 1 )

       ENDIF

!
!--    Reset values to default PALM setup
       synchronous_exchange   = synchronous_exchange_save
       send_receive = 'al'
#else

!
!--    Next line is to avoid compile error due to unused dummy argument
       IF ( color == 1234567 )  RETURN
!
!--    Standard horizontal ghost boundary exchange for small coarse grid
!--    levels, where the transfer time is latency bound
       CALL exchange_horiz( p_mg, 1 )
#endif

    END SUBROUTINE special_exchange_horiz

 END MODULE poismg_mod
