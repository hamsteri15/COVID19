!> @file user_init_urban_surface.f90
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
! Copyright 2015 Czech Technical University in Prague
! Copyright 1997-2020 Leibniz Universitaet Hannover
!--------------------------------------------------------------------------------!
!
! Current revisions:
! -----------------
! 
! 
! Former revisions:
! -----------------
! $Id: user_init_urban_surface.f90 4360 2020-01-07 11:25:50Z suehring $
! Corrected "Former revisions" section
! 
! 3768 2019-02-27 14:35:58Z raasch
! unused variables commented out to avoid compiler warnings
! 
! 3655 2019-01-07 16:51:22Z knoop
! Corrected "Former revisions" section
! 
! 2007 2016-08-24 15:47:17Z kanani
! Initial revision
!
!
! Description:
! ------------
!> Execution of user-defined actions to initiate the urban surface model
!------------------------------------------------------------------------------!
 SUBROUTINE user_init_urban_surface

    USE arrays_3d
    
!    USE control_parameters,                                                    &
!        ONLY:  urban_surface
    
    USE indices
    
    USE kinds

    USE urban_surface_mod

    USE surface_mod    
    
    USE user

    IMPLICIT NONE

!    INTEGER(iwp) ::  i  !< grid index
!    INTEGER(iwp) ::  j  !< grid index
!    INTEGER(iwp) ::  m  !< running index on 1D wall-type grid

!
!-- Here the user-defined urban surface initialization actions follow.
!-- Example: set roughness length at urban surface
!     DO  m = 1, surf_usm_h%ns
!        surf_usm_h%z0(m) = 0.1_wp
!     ENDDO



 END SUBROUTINE user_init_urban_surface

