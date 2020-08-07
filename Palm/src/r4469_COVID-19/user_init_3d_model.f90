!> @file user_init_3d_model.f90
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
! $Id: user_init_3d_model.f90 4360 2020-01-07 11:25:50Z suehring $
! Corrected "Former revisions" section
! 
! 3768 2019-02-27 14:35:58Z raasch
! variables commented out to avoid compiler warnings about unused variables
! 
! 3655 2019-01-07 16:51:22Z knoop
! Corrected "Former revisions" section
!
! 211 2008-11-11 04:46:24Z raasch
! Former file user_interface.f90 split into one file per subroutine
!
! Description:
! ------------
!> Allows the complete initialization of the 3d model.
!>
!> @attention The user is responsible to set at least all those quantities which
!>            are normally set within init_3d_model!
!------------------------------------------------------------------------------!
 SUBROUTINE user_init_3d_model
 

    USE arrays_3d
    
    USE control_parameters
    
    USE indices
    
    USE kinds

    USE surface_mod
    
    USE user

    IMPLICIT NONE

!    INTEGER(iwp) ::  l !< running index surface orientation
!    INTEGER(iwp) ::  m !< running index surface elements

!
!-- Initialization of surface-related quantities.
!-- The following example shows required initialization of surface quantitites
!-- at default-type upward-facing surfaces.  
!   DO  m = 1, surf_def_h(0)%ns
!      surf_def_h(0)%ol(m)   = ...    ! Obukhov length
!      surf_def_h(0)%us(m  ) = ...    ! friction velocity
!      surf_def_h(0)%usws(m) = ...    ! vertical momentum flux, u-component
!      surf_def_h(0)%vsws(m) = ...    ! vertical momentum flux, v-component
!      surf_def_h(0)%z0(m)   = ...    ! roughness length for momentum
!      IF ( .NOT. neutral )  THEN
!         surf_def_h(0)%ts(m)   = ... ! scaling parameter
!         surf_def_h(0)%shf(m)  = ... ! surface sensible heat flux
!         surf_def_h(0)%z0h(m)  = ... ! roughness length for heat
!      ENDIF
!      IF ( humditiy )  THEN
!         surf_def_h(0)%qs(m)   = ... ! scaling parameter
!         surf_def_h(0)%qsws(m) = ... ! surface latent heat flux
!         surf_def_h(0)%z0q(m)  = ... ! roughness length for moisture
!      ENDIF
!      IF ( passive_scalar )  THEN
!         surf_def_h(0)%ss(m)   = ... ! scaling parameter
!         surf_def_h(0)%ssws(m) = ... ! surface latent heat flux
!      ENDIF
!   ENDDO 
!
!-- Same for natural and urban type surfaces
!   DO  m = 1, surf_lsm_h%ns
!      ...
!   ENDDO 
!   DO  m = 1, surf_usm_h%ns
!      ...
!   ENDDO
!
!-- Also care for vertically aligned surfaces (default-, natural-, and 
!-- urban-type).
!   DO  l = 0, 3
!      DO  m = 1, surf_def_v(l)%ns
!         ...
!      ENDDO
!      DO  m = 1, surf_lsm_v(l)%ns
!         ...
!      ENDDO
!      DO  m = 1, surf_usm_v(l)%ns
!         ...
!      ENDDO
!   ENDDO
!
!
!-- In the following, initialize 3D quantities, e.g. u, v, w, pt, etc..

 END SUBROUTINE user_init_3d_model

