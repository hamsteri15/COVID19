!> @file init_pt_anomaly.f90
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
! $Id: init_pt_anomaly.f90 4457 2020-03-11 14:20:43Z raasch $
! use statement for exchange horiz added
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
! 3655 2019-01-07 16:51:22Z knoop
! Added topography flags
!
! Revision 1.1  1997/08/29 08:58:56  raasch
! Initial revision
!
!
! Description:
! ------------
!> Impose a temperature perturbation for an advection test.
!------------------------------------------------------------------------------!
 SUBROUTINE init_pt_anomaly


    USE arrays_3d,                                                             &
        ONLY:  pt, zu

    USE control_parameters

    USE exchange_horiz_mod,                                                    &
        ONLY:  exchange_horiz

    USE grid_variables,                                                        &
        ONLY:  dx, dy

    USE indices,                                                               &
        ONLY:  nbgp, nx, nxl, nxr, ny, nyn, nys, nzb, nzt, wall_flags_total_0

    USE kinds

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !< grid index along x
    INTEGER(iwp) ::  ic !< center index along x 
    INTEGER(iwp) ::  j  !< grid index along y
    INTEGER(iwp) ::  jc !< center index along y
    INTEGER(iwp) ::  k  !< grid index along z
    INTEGER(iwp) ::  kc !< center index along z

    REAL(wp)     ::  amount                               !< amount of temperature perturbation
    REAL(wp)     ::  bubble_center_y                      !< center of bubble in y
    REAL(wp)     ::  bubble_center_z = 170.0              !< center of bubble in z
    REAL(wp)     ::  bubble_sigma_y = 300.0               !< width of bubble in y
    REAL(wp)     ::  bubble_sigma_z = 150.0               !< width of bubble in z
    REAL(wp)     ::  flag                                 !< flag to mask topography grid points
    REAL(wp)     ::  initial_temperature_difference = 0.4 !< temperature perturbation for bubble in K
    REAL(wp)     ::  radius                               !< radius of pt anomaly
    REAL(wp)     ::  rc                                   !< radius of pt anomaly
    REAL(wp)     ::  x                                    !< x dimension of pt anomaly
    REAL(wp)     ::  y                                    !< y dimension of pt anomaly
    REAL(wp)     ::  z                                    !< z dimension of pt anomaly


!
!-- Defaults: radius rc, strength z,
!--           position of center: ic, jc, kc
    rc =  10.0_wp * dx
    ic =  ( nx+1 ) / 2
    jc =  ic
    kc =  nzt / 2

    IF ( INDEX( initializing_actions, 'initialize_ptanom' ) /= 0 )  THEN
!
!--    Compute the perturbation.
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )

                x = ( i - ic ) * dx
                y = ( j - jc ) * dy
                z = ABS( zu(k) - zu(kc) )
                radius = SQRT( x**2 + y**2 + z**2 )
                IF ( radius <= rc )  THEN
                   amount = 5.0_wp * EXP( -( radius * 0.001_wp / 2.0_wp )**2 )
                ELSE
                   amount = 0.0_wp
                ENDIF

                pt(k,j,i) = pt(k,j,i) + amount * flag

             ENDDO
          ENDDO
       ENDDO

!
!-- Initialize warm air bubble close to surface and homogenous elegonated 
!-- along x-Axis
    ELSEIF ( INDEX( initializing_actions, 'initialize_bubble' ) /= 0 )  THEN
!
!--    Calculate y-center of model domain
       bubble_center_y = ( ny + 1.0 ) * dy / 2.0

!
!--    Compute perturbation for potential temperaure
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 0 ) )

                pt(k,j,i) = pt(k,j,i) +                                        &
                               EXP( -0.5 * ( (j* dy  - bubble_center_y) /      &
                                                       bubble_sigma_y )**2) *  &
                               EXP( -0.5 * ( (zu(k)  - bubble_center_z) /      &
                                                       bubble_sigma_z)**2) *   &
                               initial_temperature_difference * flag
             ENDDO
          ENDDO
       ENDDO
    ENDIF

!
!-- Exchange of boundary values for temperature
    CALL exchange_horiz( pt, nbgp )


 END SUBROUTINE init_pt_anomaly
