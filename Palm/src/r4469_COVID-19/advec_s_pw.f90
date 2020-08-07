!> @file advec_s_pw.f90
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
! $Id: advec_s_pw.f90 4360 2020-01-07 11:25:50Z suehring $
! Corrected "Former revisions" section
! 
! 3927 2019-04-23 13:24:29Z raasch
! pointer attribute removed from scalar 3d-array for performance reasons
! 
! 3665 2019-01-10 08:28:24Z raasch
! unused variables removed
! 
! 3655 2019-01-07 16:51:22Z knoop
! nopointer option removed
!
! Revision 1.1  1997/08/29 08:54:20  raasch
! Initial revision
!
!
! Description:
! ------------
!> Advection term for scalar variables using the Piacsek and Williams scheme
!> (form C3). Contrary to PW itself, for reasons of accuracy their scheme is
!> slightly modified as follows: the values of those scalars that are used for
!> the computation of the flux divergence are reduced by the value of the
!> relevant scalar at the location where the difference is computed (sk(k,j,i)).
!> NOTE: at the first grid point above the surface computation still takes place!
!------------------------------------------------------------------------------!
 MODULE advec_s_pw_mod
 

    PRIVATE
    PUBLIC advec_s_pw

    INTERFACE advec_s_pw
       MODULE PROCEDURE advec_s_pw
       MODULE PROCEDURE advec_s_pw_ij
    END INTERFACE
 
 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE advec_s_pw( sk )

       USE arrays_3d,                                                          &
           ONLY:  dd2zu, tend, u, u_stokes_zu, v, v_stokes_zu, w

       USE control_parameters,                                                 &
           ONLY:  u_gtrans, v_gtrans

       USE grid_variables,                                                     &
           ONLY:  ddx, ddy

       USE indices,                                                            &
           ONLY:  nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg, nzb, nzt

       USE kinds


       IMPLICIT NONE

       INTEGER(iwp) ::  i !< grid index along x-direction
       INTEGER(iwp) ::  j !< grid index along y-direction
       INTEGER(iwp) ::  k !< grid index along z-direction

       REAL(wp)     ::  gu  !< local additional advective velocity
       REAL(wp)     ::  gv  !< local additional advective velocity

       REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  sk
 

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt

!
!--             Galilean transformation + Stokes drift velocity (ocean only)
                gu = u_gtrans - u_stokes_zu(k)
                gv = v_gtrans - v_stokes_zu(k)

                tend(k,j,i) = tend(k,j,i) +                                    &
                                     ( -0.5_wp * ( ( u(k,j,i+1) - gu       ) * &
                                                   ( sk(k,j,i+1) - sk(k,j,i) ) &
                                                 - ( u(k,j,i)   - gu       ) * &
                                                   ( sk(k,j,i-1) - sk(k,j,i) ) &
                                                 ) * ddx                       &
                                       -0.5_wp * ( ( v(k,j+1,i) - gv       ) * &
                                                   ( sk(k,j+1,i) - sk(k,j,i) ) &
                                                 - ( v(k,j,i)   - gv       ) * &
                                                   ( sk(k,j-1,i) - sk(k,j,i) ) &
                                                 ) * ddy                       &
                                       -         (   w(k,j,i)                * &
                                                   ( sk(k+1,j,i) - sk(k,j,i) ) &
                                                 -   w(k-1,j,i)              * &
                                                   ( sk(k-1,j,i) - sk(k,j,i) ) &
                                                 ) * dd2zu(k)                  &
                                      )
             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE advec_s_pw


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE advec_s_pw_ij( i, j, sk )

       USE arrays_3d,                                                          &
           ONLY:  dd2zu, tend, u, u_stokes_zu, v, v_stokes_zu, w

       USE control_parameters,                                                 &
           ONLY:  u_gtrans, v_gtrans

       USE grid_variables,                                                     &
           ONLY:  ddx, ddy

       USE indices,                                                            &
           ONLY:  nxlg, nxrg, nyng, nysg, nzb, nzt

       USE kinds


       IMPLICIT NONE

       INTEGER(iwp) ::  i !< grid index along x-direction
       INTEGER(iwp) ::  j !< grid index along y-direction
       INTEGER(iwp) ::  k !< grid index along z-direction

       REAL(wp)     ::  gu  !< local additional advective velocity
       REAL(wp)     ::  gv  !< local additional advective velocity

       REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  sk


       DO  k = nzb+1, nzt

!
!--       Galilean transformation + Stokes drift velocity (ocean only)
          gu = u_gtrans - u_stokes_zu(k)
          gv = v_gtrans - v_stokes_zu(k)

          tend(k,j,i) = tend(k,j,i) +                                          &
                                    ( -0.5_wp * ( ( u(k,j,i+1) - gu        ) * &
                                                  ( sk(k,j,i+1) - sk(k,j,i) )  &
                                                - ( u(k,j,i)   - gu        ) * &
                                                  ( sk(k,j,i-1) - sk(k,j,i) )  &
                                                ) * ddx                        &
                                      -0.5_wp * ( ( v(k,j+1,i) - gv       ) *  &
                                                  ( sk(k,j+1,i) - sk(k,j,i) )  &
                                                - ( v(k,j,i)   - gv       ) *  &
                                                  ( sk(k,j-1,i) - sk(k,j,i) )  &
                                                ) * ddy                        &
                                      -         (   w(k,j,i)                *  &
                                                  ( sk(k+1,j,i) - sk(k,j,i) )  &
                                                -   w(k-1,j,i)              *  &
                                                  ( sk(k-1,j,i) - sk(k,j,i) )  &
                                                ) * dd2zu(k)                   &
                                    )
       ENDDO

    END SUBROUTINE advec_s_pw_ij

 END MODULE advec_s_pw_mod
