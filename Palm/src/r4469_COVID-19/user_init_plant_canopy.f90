!> @file user_init_plant_canopy.f90
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
! $Id: user_init_plant_canopy.f90 4360 2020-01-07 11:25:50Z suehring $
! Renamed canopy_mode 'block' to 'homogeneous'
! 
! 4182 2019-08-22 15:20:23Z scharf
! Corrected "Former revisions" section
! 
! 3768 2019-02-27 14:35:58Z raasch
! unused variables commented out to avoid compiler warnings
! 
! 3655 2019-01-07 16:51:22Z knoop
! Corrected "Former revisions" section
!
! 211 2008-11-11 04:46:24Z raasch
! Former file user_interface.f90 split into one file per subroutine
!
! Description:
! ------------
!> Initialisation of the leaf area density array (for scalar grid points) and 
!> the array of the canopy drag coefficient, if the user has not chosen any 
!> of the default cases 
!------------------------------------------------------------------------------!
 SUBROUTINE user_init_plant_canopy
 

    USE arrays_3d
    
    USE control_parameters
    
    USE indices
    
    USE kinds

    USE plant_canopy_model_mod
    
    USE user

    IMPLICIT NONE

!    INTEGER(iwp) :: i   !< running index
!    INTEGER(iwp) :: j   !< running index

!
!-- Here the user-defined grid initializing actions follow:

!
!-- Set the 3D-array lad_s for user defined canopies
    SELECT CASE ( TRIM( canopy_mode ) )

       CASE ( 'homogeneous' )
!
!--       Not allowed here since this is the standard case used in init_3d_model.

       CASE ( 'user_defined_canopy_1' )
!
!--       Here the user can define his own forest topography. 
!--       The following lines contain an example, where the plant canopy extends
!--       only over the second half of the model domain along x.
!--       Attention: DO-loops have to include the ghost points (nxlg-nxrg, 
!--       nysg-nyng), because no exchange of ghost point information is intended,
!--       in order to minimize communication between CPUs.
!          DO  i = nxlg, nxrg
!             IF ( i >= INT( ( nx+1 ) / 2 ) ) THEN
!                DO  j = nysg, nyng
!                   lad_s(:,j,i) = lad(:)
!                ENDDO
!             ELSE
!                lad_s(:,:,i) = 0.0_wp
!             ENDIF
!          ENDDO 
!
!--       After definition, please
!--       remove the following three lines!
          message_string = 'canopy_mode "' // canopy_mode // &
                           '" not available yet'
          CALL message( 'user_init_plant_canopy', 'UI0007', 0, 1, 0, 6, 0 )
          
       CASE DEFAULT
!
!--       The DEFAULT case is reached if the parameter canopy_mode contains a
!--       wrong character string that is neither recognized in init_3d_model nor
!--       here in user_init_plant_canopy.
          message_string = 'unknown canopy_mode "' // canopy_mode // '"'
          CALL message( 'user_init_plant_canopy', 'UI0008', 1, 2, 0, 6, 0 )

    END SELECT


 END SUBROUTINE user_init_plant_canopy

