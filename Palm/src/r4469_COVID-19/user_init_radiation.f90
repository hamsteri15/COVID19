!> @file user_init_radiation.f90
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
! $Id: user_init_radiation.f90 4360 2020-01-07 11:25:50Z suehring $
! Corrected "Former revisions" section
! 
! 3768 2019-02-27 14:35:58Z raasch
! unused variables commented out to avoid compiler warnings
! 
! 3655 2019-01-07 16:51:22Z knoop
! Corrected "Former revisions" section
! 
! 1585 2015-04-30 07:05:52Z maronga
! Initial revision
! 
! Description:
! ------------
!> Execution of user-defined actions to initiate the radiation model
!------------------------------------------------------------------------------!
 SUBROUTINE user_init_radiation
 

    USE arrays_3d
    
    USE control_parameters
    
    USE indices
    
    USE kinds

    USE radiation_model_mod
    
    USE user

    IMPLICIT NONE

!    INTEGER(iwp) :: i   !< running index
!    INTEGER(iwp) :: j   !< running index

!
!-- Here the user-defined radiation initialization actions follow:



 END SUBROUTINE user_init_radiation

