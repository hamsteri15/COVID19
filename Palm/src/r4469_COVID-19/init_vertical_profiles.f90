!> @file ocean_mod.f90
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
! Copyright 2017-2018 Leibniz Universitaet Hannover
!--------------------------------------------------------------------------------!
!
! Current revisions:
! -----------------
! 
! 
! Former revisions:
! -----------------
! $Id: init_vertical_profiles.f90 3294 2018-10-01 02:37:10Z raasch $
! split from check_parameters as separate file to avoid circular dependency
! with ocean_mod
!
! 
!
!
! Authors:
! --------
! @author Siegfried Raasch
!
! Description:
! ------------
!> Inititalizes the vertical profiles of scalar quantities.
!------------------------------------------------------------------------------!
 SUBROUTINE init_vertical_profiles( vertical_gradient_level_ind,               &
                                    vertical_gradient_level,                   &
                                    vertical_gradient, initial_profile,        &
                                    surface_value, bc_top_gradient )

    USE arrays_3d,                                                             &
        ONLY:  dzu, zu

    USE control_parameters,                                                    &
        ONLY:  ocean_mode

    USE indices,                                                               &
        ONLY:  nz, nzt

    USE kinds

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !< loop counter
    INTEGER(iwp) ::  k  !< loop counter
    INTEGER(iwp), DIMENSION(1:10) ::  vertical_gradient_level_ind  !< vertical grid indices for gradient levels

    REAL(wp)     ::  bc_top_gradient  !< model top gradient
    REAL(wp)     ::  gradient         !< vertica gradient of the respective quantity
    REAL(wp)     ::  surface_value    !< surface value of the respecitve quantity

    REAL(wp), DIMENSION(0:nz+1) ::  initial_profile          !< initialisation profile
    REAL(wp), DIMENSION(1:10)   ::  vertical_gradient        !< given vertical gradient
    REAL(wp), DIMENSION(1:10)   ::  vertical_gradient_level  !< given vertical gradient level

    i = 1
    gradient = 0.0_wp

    IF ( .NOT. ocean_mode )  THEN

       vertical_gradient_level_ind(1) = 0
       DO  k = 1, nzt+1
          IF ( i < 11 )  THEN
             IF ( vertical_gradient_level(i) < zu(k)  .AND.            &
                  vertical_gradient_level(i) >= 0.0_wp )  THEN
                gradient = vertical_gradient(i) / 100.0_wp
                vertical_gradient_level_ind(i) = k - 1
                i = i + 1
             ENDIF
          ENDIF
          IF ( gradient /= 0.0_wp )  THEN
             IF ( k /= 1 )  THEN
                initial_profile(k) = initial_profile(k-1) + dzu(k) * gradient
             ELSE
                initial_profile(k) = initial_profile(k-1) + dzu(k) * gradient
             ENDIF
          ELSE
             initial_profile(k) = initial_profile(k-1)
          ENDIF
!
!--       Avoid negative values of scalars
          IF ( initial_profile(k) < 0.0_wp )  THEN
             initial_profile(k) = 0.0_wp
          ENDIF
       ENDDO

    ELSE

!
!--    In ocean mode, profiles are constructed starting from the ocean surface,
!--    which is at the top of the model domain
       vertical_gradient_level_ind(1) = nzt+1
       DO  k = nzt, 0, -1
          IF ( i < 11 )  THEN
             IF ( vertical_gradient_level(i) > zu(k)  .AND.            &
                  vertical_gradient_level(i) <= 0.0_wp )  THEN
                gradient = vertical_gradient(i) / 100.0_wp
                vertical_gradient_level_ind(i) = k + 1
                i = i + 1
             ENDIF
          ENDIF
          IF ( gradient /= 0.0_wp )  THEN
             IF ( k /= nzt )  THEN
                initial_profile(k) = initial_profile(k+1) - dzu(k+1) * gradient
             ELSE
                initial_profile(k)   = surface_value - 0.5_wp * dzu(k+1) *     &
                                                       gradient
                initial_profile(k+1) = surface_value + 0.5_wp * dzu(k+1) *     &
                                                       gradient
             ENDIF
          ELSE
             initial_profile(k) = initial_profile(k+1)
          ENDIF
!
!--       Avoid negative values of scalars
          IF ( initial_profile(k) < 0.0_wp )  THEN
             initial_profile(k) = 0.0_wp
          ENDIF
       ENDDO

    ENDIF

!
!-- In case of no given gradients, choose zero gradient conditions
    IF ( vertical_gradient_level(1) == -999999.9_wp )  THEN
       vertical_gradient_level(1) = 0.0_wp
    ENDIF
!
!-- Store gradient at the top boundary for possible Neumann boundary condition
    bc_top_gradient  = ( initial_profile(nzt+1) - initial_profile(nzt) ) /     &
                       dzu(nzt+1)

 END SUBROUTINE init_vertical_profiles
