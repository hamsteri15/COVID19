!> @file advec_w_pw.f90
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
! $Id: advec_w_pw.f90 4360 2020-01-07 11:25:50Z suehring $
! Corrected "Former revisions" section
! 
! 3655 2019-01-07 16:51:22Z knoop
! variables documented
!
! Revision 1.1  1997/08/11 06:10:29  raasch
! Initial revision
!
!
! Description:
! ------------
!> Advection term for w velocity-component using Piacsek and Williams.
!> Vertical advection at the first grid point above the surface is done with
!> normal centred differences, because otherwise no information from the surface
!> would be communicated upwards due to w=0 at k=nzb.
!------------------------------------------------------------------------------!
 MODULE advec_w_pw_mod
 

    PRIVATE
    PUBLIC advec_w_pw

    INTERFACE advec_w_pw
       MODULE PROCEDURE advec_w_pw
       MODULE PROCEDURE advec_w_pw_ij
    END INTERFACE advec_w_pw
 
 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE advec_w_pw

       USE arrays_3d,                                                          &
           ONLY:  ddzu, tend, u, v, w

       USE control_parameters,                                                 &
           ONLY:  u_gtrans, v_gtrans

       USE grid_variables,                                                     &
           ONLY:  ddx, ddy

       USE indices,                                                            &
           ONLY:  nxl, nxr, nyn, nys, nzb, nzt

       USE kinds


       IMPLICIT NONE

       INTEGER(iwp) ::  i !< grid index along x-direction
       INTEGER(iwp) ::  j !< grid index along y-direction
       INTEGER(iwp) ::  k !< grid index along z-direction
       
       REAL(wp)    ::  gu !< Galilei-transformation velocity along x
       REAL(wp)    ::  gv !< Galilei-transformation velocity along y

 
       gu = 2.0_wp * u_gtrans
       gv = 2.0_wp * v_gtrans
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                tend(k,j,i) = tend(k,j,i) - 0.25_wp * (                        &
                         ( w(k,j,i+1) * ( u(k+1,j,i+1) + u(k,j,i+1) - gu )     &
                         - w(k,j,i-1) * ( u(k+1,j,i) + u(k,j,i) - gu ) ) * ddx &
                       + ( w(k,j+1,i) * ( v(k+1,j+1,i) + v(k,j+1,i) - gv )     &
                         - w(k,j-1,i) * ( v(k+1,j,i) + v(k,j,i) - gv ) ) * ddy &
                       + ( w(k+1,j,i) * ( w(k+1,j,i) + w(k,j,i) )              &
                         - w(k-1,j,i) * ( w(k,j,i) + w(k-1,j,i) ) )            &
                                                                  * ddzu(k+1)  &
                                                      )
             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE advec_w_pw


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE advec_w_pw_ij( i, j )

       USE arrays_3d,                                                          &
           ONLY:  ddzu, tend, u, v, w

       USE control_parameters,                                                 &
           ONLY:  u_gtrans, v_gtrans

       USE grid_variables,                                                     &
           ONLY:  ddx, ddy

       USE indices,                                                            &
           ONLY:  nzb, nzt

       USE kinds


       IMPLICIT NONE

       INTEGER(iwp) ::  i !< grid index along x-direction
       INTEGER(iwp) ::  j !< grid index along y-direction
       INTEGER(iwp) ::  k !< grid index along z-direction
       
       REAL(wp)    ::  gu !< Galilei-transformation velocity along x
       REAL(wp)    ::  gv !< Galilei-transformation velocity along y

       gu = 2.0_wp * u_gtrans
       gv = 2.0_wp * v_gtrans
       DO  k = nzb+1, nzt
          tend(k,j,i) = tend(k,j,i) - 0.25_wp * (                              &
                         ( w(k,j,i+1) * ( u(k+1,j,i+1) + u(k,j,i+1) - gu )     &
                         - w(k,j,i-1) * ( u(k+1,j,i) + u(k,j,i) - gu ) ) * ddx &
                       + ( w(k,j+1,i) * ( v(k+1,j+1,i) + v(k,j+1,i) - gv )     &
                         - w(k,j-1,i) * ( v(k+1,j,i) + v(k,j,i) - gv ) ) * ddy &
                       + ( w(k+1,j,i) * ( w(k+1,j,i) + w(k,j,i) )              &
                         - w(k-1,j,i) * ( w(k,j,i) + w(k-1,j,i) ) )            &
                                                                  * ddzu(k+1)  &
                                                )
       ENDDO
    END SUBROUTINE advec_w_pw_ij

 END MODULE advec_w_pw_mod
