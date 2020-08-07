!> @file coriolis.f90
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
! $Id: coriolis.f90 4360 2020-01-07 11:25:50Z suehring $
! Introduction of wall_flags_total_0, which currently sets bits based on static
! topography information used in wall_flags_static_0
! 
! 4329 2019-12-10 15:46:36Z motisi
! Renamed wall_flags_0 to wall_flags_static_0
! 
! 4196 2019-08-29 11:02:06Z gronemeier
! Consider rotation of model domain
! 
! 4182 2019-08-22 15:20:23Z scharf
! Corrected "Former revisions" section
! 
! 3655 2019-01-07 16:51:22Z knoop
! OpenACC port for SPEC
!
! Revision 1.1  1997/08/29 08:57:38  raasch
! Initial revision
!
!
! Description:
! ------------
!> Computation of all Coriolis terms in the equations of motion.
!> 
!> @note In this routine the topography is masked, even though this 
!>       is again done in prognostic_equations. However, omitting the masking 
!>       here lead to slightly different results. Reason unknown. 
!------------------------------------------------------------------------------!
 MODULE coriolis_mod
 

    PRIVATE
    PUBLIC coriolis

    INTERFACE coriolis
       MODULE PROCEDURE coriolis
       MODULE PROCEDURE coriolis_ij
    END INTERFACE coriolis

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE coriolis( component )

       USE arrays_3d,                                                          &
           ONLY:  tend, u, ug, v, vg, w 

       USE basic_constants_and_equations_mod,                                  &
           ONLY:  pi

       USE control_parameters,                                                 &
           ONLY:  f, fs, message_string, rotation_angle
           
       USE indices,                                                            &
           ONLY:  nxl, nxlu, nxr, nyn, nys, nysv, nzb, nzt,                    &
                  wall_flags_total_0
                   
       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  component      !< component of momentum equation
       INTEGER(iwp) ::  i              !< running index x direction 
       INTEGER(iwp) ::  j              !< running index y direction 
       INTEGER(iwp) ::  k              !< running index z direction 

       REAL(wp)     ::  cos_rot_angle  !< cosine of model rotation angle
       REAL(wp)     ::  flag           !< flag to mask topography
       REAL(wp)     ::  sin_rot_angle  !< sine of model rotation angle

!
!--    Precalculate cosine and sine of rotation angle
       cos_rot_angle = COS( rotation_angle * pi / 180.0_wp )
       sin_rot_angle = SIN( rotation_angle * pi / 180.0_wp )

!
!--    Compute Coriolis terms for the three velocity components
       SELECT CASE ( component )

!
!--       u-component
          CASE ( 1 )
             !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i, j, k, flag) &
             !$ACC PRESENT(wall_flags_total_0) &
             !$ACC PRESENT(v, w, vg) &
             !$ACC PRESENT(tend)
             DO  i = nxlu, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt
!
!--                   Predetermine flag to mask topography
                      flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 1 ) )

                      tend(k,j,i) = tend(k,j,i) + flag *                                           &
                            ( f                                                                    &
                              * ( 0.25_wp * ( v(k,j,i-1) + v(k,j,i) + v(k,j+1,i-1) + v(k,j+1,i) )  &
                                - vg(k) )                                                          &
                            - fs * cos_rot_angle                                                   &
                              * 0.25_wp * ( w(k-1,j,i-1) + w(k-1,j,i) + w(k,j,i-1) + w(k,j,i) )    &
                            )
                   ENDDO
                ENDDO
             ENDDO

!
!--       v-component
          CASE ( 2 )
             !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i, j, k, flag) &
             !$ACC PRESENT(wall_flags_total_0) &
             !$ACC PRESENT(u, w, ug) &
             !$ACC PRESENT(tend)
             DO  i = nxl, nxr
                DO  j = nysv, nyn
                   DO  k = nzb+1, nzt
!
!--                   Predetermine flag to mask topography
                      flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 2 ) )

                      tend(k,j,i) = tend(k,j,i) - flag *                                           &
                            ( f                                                                    &
                              * ( 0.25_wp * ( u(k,j-1,i) + u(k,j,i) + u(k,j-1,i+1) + u(k,j,i+1) )  &
                                - ug(k) )                                                          &
                            + fs * sin_rot_angle                                                   &
                              * 0.25_wp * ( w(k,j,i) + w(k-1,j,i) + w(k,j-1,i) + w(k-1,j-1,i) )    &
                            )
                   ENDDO
                ENDDO
             ENDDO

!
!--       w-component
          CASE ( 3 )
             !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i, j, k, flag) &
             !$ACC PRESENT(wall_flags_total_0) &
             !$ACC PRESENT(u, v) &
             !$ACC PRESENT(tend)
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt
!
!--                   Predetermine flag to mask topography
                      flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 3 ) )

                      tend(k,j,i) = tend(k,j,i)                                                 &
                                  + fs * 0.25_wp * flag                                         &
                                    * ( cos_rot_angle                                           &
                                        * ( u(k,j,i) + u(k+1,j,i) + u(k,j,i+1) + u(k+1,j,i+1) ) &
                                      + sin_rot_angle                                           &
                                        * ( v(k,j,i) + v(k+1,j,i) + v(k,j+1,i) + v(k+1,j+1,i) ) &
                                      )
                   ENDDO
                ENDDO
             ENDDO

          CASE DEFAULT

             WRITE( message_string, * ) ' wrong component: ', component
             CALL message( 'coriolis', 'PA0173', 1, 2, 0, 6, 0 )

       END SELECT

    END SUBROUTINE coriolis


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE coriolis_ij( i, j, component )

       USE arrays_3d,                                                          &
           ONLY:  tend, u, ug, v, vg, w 

       USE basic_constants_and_equations_mod,                                  &
           ONLY:  pi

       USE control_parameters,                                                 &
           ONLY:  f, fs, message_string, rotation_angle
           
       USE indices,                                                            &
           ONLY:  nzb, nzt, wall_flags_total_0
           
       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  component  !< component of momentum equation
       INTEGER(iwp) ::  i          !< running index x direction 
       INTEGER(iwp) ::  j          !< running index y direction 
       INTEGER(iwp) ::  k          !< running index z direction 

       REAL(wp)     ::  cos_rot_angle  !< cosine of model rotation angle
       REAL(wp)     ::  flag           !< flag to mask topography
       REAL(wp)     ::  sin_rot_angle  !< sine of model rotation angle

!
!--    Precalculate cosine and sine of rotation angle
       cos_rot_angle = COS( rotation_angle * pi / 180.0_wp )
       sin_rot_angle = SIN( rotation_angle * pi / 180.0_wp )

!
!--    Compute Coriolis terms for the three velocity components
       SELECT CASE ( component )

!
!--       u-component
          CASE ( 1 )
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 1 ) )

                tend(k,j,i) = tend(k,j,i) + flag *                                                 &
                            ( f                                                                    &
                              * ( 0.25_wp * ( v(k,j,i-1) + v(k,j,i) + v(k,j+1,i-1) + v(k,j+1,i) )  &
                                - vg(k) )                                                          &
                            - fs * cos_rot_angle                                                   &
                              * 0.25_wp * ( w(k-1,j,i-1) + w(k-1,j,i) + w(k,j,i-1) + w(k,j,i) )    &
                            )
             ENDDO

!
!--       v-component
          CASE ( 2 )
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 2 ) )

                tend(k,j,i) = tend(k,j,i) - flag *                                                 &
                            ( f                                                                    &
                              * ( 0.25_wp * ( u(k,j-1,i) + u(k,j,i) + u(k,j-1,i+1) + u(k,j,i+1) )  &
                                - ug(k) )                                                          &
                            + fs * sin_rot_angle                                                   &
                              * 0.25_wp * ( w(k,j,i) + w(k-1,j,i) + w(k,j-1,i) + w(k-1,j-1,i) )    &
                            )
             ENDDO

!
!--       w-component
          CASE ( 3 )
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 3 ) )

                tend(k,j,i) = tend(k,j,i)                                                 &
                            + fs * 0.25_wp * flag                                         &
                              * ( cos_rot_angle                                           &
                                  * ( u(k,j,i) + u(k+1,j,i) + u(k,j,i+1) + u(k+1,j,i+1) ) &
                                + sin_rot_angle                                           &
                                  * ( v(k,j,i) + v(k+1,j,i) + v(k,j+1,i) + v(k+1,j+1,i) ) &
                                )
             ENDDO

          CASE DEFAULT

             WRITE( message_string, * ) ' wrong component: ', component
             CALL message( 'coriolis', 'PA0173', 1, 2, 0, 6, 0 )

       END SELECT

    END SUBROUTINE coriolis_ij

 END MODULE coriolis_mod
