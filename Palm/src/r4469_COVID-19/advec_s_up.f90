!> @file advec_s_up.f90
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
! $Id: advec_s_up.f90 4360 2020-01-07 11:25:50Z suehring $
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
! Revision 1.1  1997/08/29 08:54:33  raasch
! Initial revision
!
!
! Description:
! ------------
!> Advection term for scalar quantities using the Upstream scheme.
!> NOTE: vertical advection at k=1 still has wrong grid spacing for w>0!
!>       The same problem occurs for all topography boundaries!
!------------------------------------------------------------------------------!
 MODULE advec_s_up_mod
 

    PRIVATE
    PUBLIC advec_s_up

    INTERFACE advec_s_up
       MODULE PROCEDURE advec_s_up
       MODULE PROCEDURE advec_s_up_ij
    END INTERFACE advec_s_up

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE advec_s_up( sk )

       USE arrays_3d,                                                          &
           ONLY:  ddzu, tend, u, v, w

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

       REAL(wp) ::  ukomp !< advection velocity along x-direction 
       REAL(wp) ::  vkomp !< advection velocity along y-direction 
       REAL(wp) ::  wkomp !< advection velocity along z-direction 

       REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  sk !< treated scalar


       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             x-direction
                ukomp = 0.5_wp * ( u(k,j,i) + u(k,j,i+1) ) - u_gtrans
                IF ( ukomp > 0.0_wp )  THEN
                   tend(k,j,i) = tend(k,j,i) - ukomp *                         &
                                         ( sk(k,j,i) - sk(k,j,i-1) ) * ddx
                ELSE
                   tend(k,j,i) = tend(k,j,i) - ukomp *                         &
                                         ( sk(k,j,i+1) - sk(k,j,i) ) * ddx
                ENDIF
!
!--             y-direction
                vkomp = 0.5_wp * ( v(k,j,i) + v(k,j+1,i) ) - v_gtrans
                IF ( vkomp > 0.0_wp )  THEN
                   tend(k,j,i) = tend(k,j,i) - vkomp *                         &
                                         ( sk(k,j,i) - sk(k,j-1,i) ) * ddy
                ELSE
                   tend(k,j,i) = tend(k,j,i) - vkomp *                         &
                                         ( sk(k,j+1,i) - sk(k,j,i) ) * ddy
                ENDIF
!
!--             z-direction
                wkomp = 0.5_wp * ( w(k,j,i) + w(k-1,j,i) )
                IF ( wkomp > 0.0_wp )  THEN
                   tend(k,j,i) = tend(k,j,i) - wkomp *                         &
                                         ( sk(k,j,i) - sk(k-1,j,i) ) * ddzu(k)
                ELSE
                   tend(k,j,i) = tend(k,j,i) - wkomp *                         &
                                         ( sk(k+1,j,i)-sk(k,j,i) ) * ddzu(k+1)
                ENDIF

             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE advec_s_up


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE advec_s_up_ij( i, j, sk )

       USE arrays_3d,                                                          &
           ONLY:  ddzu, tend, u, v, w

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

       REAL(wp) ::  ukomp !< advection velocity along x-direction 
       REAL(wp) ::  vkomp !< advection velocity along y-direction 
       REAL(wp) ::  wkomp !< advection velocity along z-direction 
       
       REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  sk !< treated scalar


       DO  k = nzb+1, nzt
!
!--       x-direction
          ukomp = 0.5_wp * ( u(k,j,i) + u(k,j,i+1) ) - u_gtrans
          IF ( ukomp > 0.0_wp )  THEN
             tend(k,j,i) = tend(k,j,i) - ukomp *                               &
                                         ( sk(k,j,i) - sk(k,j,i-1) ) * ddx
          ELSE
             tend(k,j,i) = tend(k,j,i) - ukomp *                               &
                                         ( sk(k,j,i+1) - sk(k,j,i) ) * ddx
          ENDIF
!
!--       y-direction
          vkomp = 0.5_wp * ( v(k,j,i) + v(k,j+1,i) ) - v_gtrans
          IF ( vkomp > 0.0_wp )  THEN
             tend(k,j,i) = tend(k,j,i) - vkomp *                               &
                                         ( sk(k,j,i) - sk(k,j-1,i) ) * ddy
          ELSE
             tend(k,j,i) = tend(k,j,i) - vkomp *                               &
                                         ( sk(k,j+1,i) - sk(k,j,i) ) * ddy
          ENDIF
!
!--       z-direction
          wkomp = 0.5_wp * ( w(k,j,i) + w(k-1,j,i) )
          IF ( wkomp > 0.0_wp )  THEN
             tend(k,j,i) = tend(k,j,i) - wkomp *                               &
                                         ( sk(k,j,i) - sk(k-1,j,i) ) * ddzu(k)
          ELSE
             tend(k,j,i) = tend(k,j,i) - wkomp *                               &
                                         ( sk(k+1,j,i)-sk(k,j,i) ) * ddzu(k+1)
          ENDIF

       ENDDO

    END SUBROUTINE advec_s_up_ij

 END MODULE advec_s_up_mod
