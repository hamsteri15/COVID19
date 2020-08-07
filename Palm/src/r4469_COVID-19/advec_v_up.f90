!> @file advec_v_up.f90
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
! $Id: advec_v_up.f90 4360 2020-01-07 11:25:50Z suehring $
! Corrected "Former revisions" section
! 
! 3655 2019-01-07 16:51:22Z knoop
! variables documented
!
! Revision 1.1  1997/08/29 08:56:05  raasch
! Initial revision
!
!
! Description:
! ------------
!> Advection term for the v velocity-component using upstream scheme.
!> NOTE: vertical advection at k=1 still has wrong grid spacing for w>0!
!>       The same problem occurs for all topography boundaries!
!------------------------------------------------------------------------------!
 MODULE advec_v_up_mod
 

    PRIVATE
    PUBLIC advec_v_up

    INTERFACE advec_v_up
       MODULE PROCEDURE advec_v_up
       MODULE PROCEDURE advec_v_up_ij
    END INTERFACE advec_v_up

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE advec_v_up

       USE arrays_3d,                                                          &
           ONLY:  ddzu, tend, u, v, w

       USE control_parameters,                                                 &
           ONLY:  u_gtrans, v_gtrans

       USE grid_variables,                                                     &
           ONLY:  ddx, ddy

       USE indices,                                                            &
           ONLY:  nxl, nxr, nyn, nysv, nzb, nzt

       USE kinds


       IMPLICIT NONE

       INTEGER(iwp) ::  i !< grid index along x-direction
       INTEGER(iwp) ::  j !< grid index along y-direction
       INTEGER(iwp) ::  k !< grid index along z-direction

       REAL(wp) ::  ukomp !< advection velocity along x-direction 
       REAL(wp) ::  vkomp !< advection velocity along y-direction 
       REAL(wp) ::  wkomp !< advection velocity along z-direction       


       DO  i = nxl, nxr
          DO  j = nysv, nyn
             DO  k = nzb+1, nzt
!
!--             x-direction
                ukomp = 0.25_wp * ( u(k,j,i)   + u(k,j-1,i) +                  &
                                 u(k,j,i+1) + u(k,j-1,i+1) ) - u_gtrans
                IF ( ukomp > 0.0_wp )  THEN
                   tend(k,j,i) = tend(k,j,i) - ukomp *                         &
                                         ( v(k,j,i) - v(k,j,i-1) ) * ddx
                ELSE
                   tend(k,j,i) = tend(k,j,i) - ukomp *                         &
                                         ( v(k,j,i+1) - v(k,j,i) ) * ddx
                ENDIF
!
!--             y-direction
                vkomp = v(k,j,i) - v_gtrans
                IF ( vkomp > 0.0_wp )  THEN
                   tend(k,j,i) = tend(k,j,i) - vkomp *                         &
                                         ( v(k,j,i) - v(k,j-1,i) ) * ddy
                ELSE
                   tend(k,j,i) = tend(k,j,i) - vkomp *                         &
                                         ( v(k,j+1,i) - v(k,j,i) ) * ddy
                ENDIF
!
!--             z-direction
                wkomp = 0.25_wp * ( w(k,j,i)  + w(k-1,j,i) +                   &
                                 w(k,j-1,i) + w(k-1,j-1,i) )
                IF ( wkomp > 0.0_wp )  THEN
                   tend(k,j,i) = tend(k,j,i) - wkomp *                         &
                                         ( v(k,j,i) - v(k-1,j,i) ) * ddzu(k)
                ELSE
                   tend(k,j,i) = tend(k,j,i) - wkomp *                         &
                                         ( v(k+1,j,i) - v(k,j,i) ) * ddzu(k+1)
                ENDIF

             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE advec_v_up


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE advec_v_up_ij( i, j )

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

       REAL(wp) ::  ukomp !< advection velocity along x-direction 
       REAL(wp) ::  vkomp !< advection velocity along y-direction 
       REAL(wp) ::  wkomp !< advection velocity along z-direction 


       DO  k = nzb+1, nzt
!
!--       x-direction
          ukomp = 0.25_wp * ( u(k,j,i) + u(k,j-1,i) + u(k,j,i+1) + u(k,j-1,i+1) &
                         ) - u_gtrans
          IF ( ukomp > 0.0_wp )  THEN
             tend(k,j,i) = tend(k,j,i) - ukomp *                               &
                                         ( v(k,j,i) - v(k,j,i-1) ) * ddx
          ELSE
             tend(k,j,i) = tend(k,j,i) - ukomp *                               &
                                         ( v(k,j,i+1) - v(k,j,i) ) * ddx
          ENDIF
!
!--       y-direction
          vkomp = v(k,j,i) - v_gtrans
          IF ( vkomp > 0.0_wp )  THEN
             tend(k,j,i) = tend(k,j,i) - vkomp *                               &
                                         ( v(k,j,i) - v(k,j-1,i) ) * ddy
          ELSE
             tend(k,j,i) = tend(k,j,i) - vkomp *                               &
                                         ( v(k,j+1,i) - v(k,j,i) ) * ddy
          ENDIF
!
!--       z-direction
          wkomp = 0.25_wp * ( w(k,j,i) + w(k-1,j,i) + w(k,j-1,i) + w(k-1,j-1,i) )
          IF ( wkomp > 0.0_wp )  THEN
             tend(k,j,i) = tend(k,j,i) - wkomp *                               &
                                         ( v(k,j,i) - v(k-1,j,i) ) * ddzu(k)
          ELSE
             tend(k,j,i) = tend(k,j,i) - wkomp *                               &
                                         ( v(k+1,j,i) - v(k,j,i) ) * ddzu(k+1)
          ENDIF

       ENDDO

    END SUBROUTINE advec_v_up_ij

 END MODULE advec_v_up_mod
