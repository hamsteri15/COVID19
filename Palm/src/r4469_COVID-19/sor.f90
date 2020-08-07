!> @file sor.f90
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
! $Id: sor.f90 4457 2020-03-11 14:20:43Z raasch $
! use statement for exchange horiz added
! 
! 4360 2020-01-07 11:25:50Z suehring
! Corrected "Former revisions" section
! 
! 3655 2019-01-07 16:51:22Z knoop
! Rename variables in mesoscale-offline nesting mode
!
! Revision 1.1  1997/08/11 06:25:56  raasch
! Initial revision
!
!
! Description:
! ------------
!> Solve the Poisson-equation with the SOR-Red/Black-scheme.
!------------------------------------------------------------------------------!
 SUBROUTINE sor( d, ddzu, ddzw, p )

    USE arrays_3d,                                                             &
        ONLY:  rho_air, rho_air_zw

    USE control_parameters,                                                    &
        ONLY:  bc_dirichlet_l, bc_dirichlet_n, bc_dirichlet_r,                 &
               bc_dirichlet_s, bc_lr_cyc, bc_ns_cyc, bc_radiation_l,           &
               bc_radiation_n, bc_radiation_r, bc_radiation_s, ibc_p_b,        &
               ibc_p_t, n_sor, omega_sor

    USE exchange_horiz_mod,                                                    &
        ONLY:  exchange_horiz

    USE grid_variables,                                                        &
        ONLY:  ddx2, ddy2

    USE indices,                                                               &
        ONLY:  nbgp, nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg, nz, nzb, nzt

    USE kinds

    IMPLICIT NONE

    INTEGER(iwp) ::  i              !<
    INTEGER(iwp) ::  j              !<
    INTEGER(iwp) ::  k              !<
    INTEGER(iwp) ::  n              !<
    INTEGER(iwp) ::  nxl1           !<
    INTEGER(iwp) ::  nxl2           !<
    INTEGER(iwp) ::  nys1           !<
    INTEGER(iwp) ::  nys2           !<

    REAL(wp)     ::  ddzu(1:nz+1)   !<
    REAL(wp)     ::  ddzw(1:nzt+1)  !<

    REAL(wp)     ::  d(nzb+1:nzt,nys:nyn,nxl:nxr)      !<
    REAL(wp)     ::  p(nzb:nzt+1,nysg:nyng,nxlg:nxrg)  !<

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  f1         !<
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  f2         !<
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  f3         !<

    ALLOCATE( f1(1:nz), f2(1:nz), f3(1:nz) )

!
!-- Compute pre-factors.
    DO  k = 1, nz
         f2(k) = ddzu(k+1) * ddzw(k) * rho_air_zw(k)
         f3(k) = ddzu(k)   * ddzw(k) * rho_air_zw(k-1)
         f1(k) = 2.0_wp * ( ddx2 + ddy2 ) * rho_air(k) + f2(k) + f3(k)
    ENDDO

!
!-- Limits for RED- and BLACK-part.
    IF ( MOD( nxl , 2 ) == 0 )  THEN
       nxl1 = nxl
       nxl2 = nxl + 1
    ELSE
       nxl1 = nxl + 1
       nxl2 = nxl
    ENDIF
    IF ( MOD( nys , 2 ) == 0 )  THEN
       nys1 = nys
       nys2 = nys + 1
    ELSE
       nys1 = nys + 1
       nys2 = nys
    ENDIF

    DO  n = 1, n_sor

!
!--    RED-part
       DO  i = nxl1, nxr, 2
          DO  j = nys2, nyn, 2
             DO  k = nzb+1, nzt
                p(k,j,i) = p(k,j,i) + omega_sor / f1(k) * (            &
                           rho_air(k) * ddx2 * ( p(k,j,i+1) + p(k,j,i-1) ) +   &
                           rho_air(k) * ddy2 * ( p(k,j+1,i) + p(k,j-1,i) ) +   &
                           f2(k) * p(k+1,j,i)                              +   &
                           f3(k) * p(k-1,j,i)                              -   &
                           d(k,j,i)                                        -   &
                           f1(k) * p(k,j,i)           )
             ENDDO
          ENDDO
       ENDDO

       DO  i = nxl2, nxr, 2
          DO  j = nys1, nyn, 2
             DO  k = nzb+1, nzt
                p(k,j,i) = p(k,j,i) + omega_sor / f1(k) * (                    &
                           rho_air(k) * ddx2 * ( p(k,j,i+1) + p(k,j,i-1) ) +   &
                           rho_air(k) * ddy2 * ( p(k,j+1,i) + p(k,j-1,i) ) +   &
                           f2(k) * p(k+1,j,i)                              +   &
                           f3(k) * p(k-1,j,i)                              -   &
                           d(k,j,i)                                        -   &
                           f1(k) * p(k,j,i)           )
             ENDDO
          ENDDO
       ENDDO

!
!--    Exchange of boundary values for p.
       CALL exchange_horiz( p, nbgp )

!
!--    Horizontal (Neumann) boundary conditions in case of non-cyclic boundaries
       IF ( .NOT. bc_lr_cyc )  THEN
          IF ( bc_dirichlet_l  .OR.  bc_radiation_l )  p(:,:,nxl-1) = p(:,:,nxl)
          IF ( bc_dirichlet_r  .OR.  bc_radiation_r )  p(:,:,nxr+1) = p(:,:,nxr)
       ENDIF
       IF ( .NOT. bc_ns_cyc )  THEN
          IF ( bc_dirichlet_n  .OR.  bc_radiation_n )  p(:,nyn+1,:) = p(:,nyn,:)
          IF ( bc_dirichlet_s  .OR.  bc_radiation_s )  p(:,nys-1,:) = p(:,nys,:)
       ENDIF

!
!--    BLACK-part
       DO  i = nxl1, nxr, 2
          DO  j = nys1, nyn, 2
             DO  k = nzb+1, nzt
                p(k,j,i) = p(k,j,i) + omega_sor / f1(k) * (            &
                           rho_air(k) * ddx2 * ( p(k,j,i+1) + p(k,j,i-1) ) +   &
                           rho_air(k) * ddy2 * ( p(k,j+1,i) + p(k,j-1,i) ) +   &
                           f2(k) * p(k+1,j,i)                              +   &
                           f3(k) * p(k-1,j,i)                              -   &
                           d(k,j,i)                                        -   &
                           f1(k) * p(k,j,i)           )
             ENDDO
          ENDDO
       ENDDO

       DO  i = nxl2, nxr, 2
          DO  j = nys2, nyn, 2
             DO  k = nzb+1, nzt
                p(k,j,i) = p(k,j,i) + omega_sor / f1(k) * (            &
                           rho_air(k) * ddx2 * ( p(k,j,i+1) + p(k,j,i-1) ) +   &
                           rho_air(k) * ddy2 * ( p(k,j+1,i) + p(k,j-1,i) ) +   &
                           f2(k) * p(k+1,j,i)                              +   &
                           f3(k) * p(k-1,j,i)                              -   &
                           d(k,j,i)                                        -   &
                           f1(k) * p(k,j,i)           )
             ENDDO
          ENDDO
       ENDDO

!
!--    Exchange of boundary values for p.
       CALL exchange_horiz( p, nbgp )

!
!--    Boundary conditions top/bottom.
!--    Bottom boundary
       IF ( ibc_p_b == 1 )  THEN       !       Neumann
          p(nzb,:,:) = p(nzb+1,:,:)
       ELSE                            !       Dirichlet
          p(nzb,:,:) = 0.0_wp
       ENDIF

!
!--    Top boundary
       IF ( ibc_p_t == 1 )  THEN                 !  Neumann
          p(nzt+1,:,:) = p(nzt,:,:)
       ELSE                      !  Dirichlet
          p(nzt+1,:,:) = 0.0_wp
       ENDIF

!
!--    Horizontal (Neumann) boundary conditions in case of non-cyclic boundaries
       IF ( .NOT. bc_lr_cyc )  THEN
          IF ( bc_dirichlet_l  .OR.  bc_radiation_l )  p(:,:,nxl-1) = p(:,:,nxl)
          IF ( bc_dirichlet_r  .OR.  bc_radiation_r )  p(:,:,nxr+1) = p(:,:,nxr)
       ENDIF
       IF ( .NOT. bc_ns_cyc )  THEN
          IF ( bc_dirichlet_n  .OR.  bc_radiation_n )  p(:,nyn+1,:) = p(:,nyn,:)
          IF ( bc_dirichlet_s  .OR.  bc_radiation_s )  p(:,nys-1,:) = p(:,nys,:)
       ENDIF


    ENDDO

    DEALLOCATE( f1, f2, f3 )

 END SUBROUTINE sor
