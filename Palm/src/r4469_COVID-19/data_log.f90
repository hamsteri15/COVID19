!> @file data_log.f90
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
! $Id: data_log.f90 4360 2020-01-07 11:25:50Z suehring $
! Corrected "Former revisions" section
! 
! 3725 2019-02-07 10:11:02Z raasch
! preprocessor directives removed to avoid compiler warnings
! 
! 3655 2019-01-07 16:51:22Z knoop
! Corrected "Former revisions" section
!
! Revision 1.1  2006/02/23 10:09:29  raasch
! Initial revision
!
!
! Description:
! ------------
!> Complete logging of data
!------------------------------------------------------------------------------!
 SUBROUTINE data_log( array, i1, i2, j1, j2, k1, k2 )
 
    USE control_parameters,                                                    &
        ONLY:  log_message, simulated_time
        
    USE kinds
        
    USE pegrid

    IMPLICIT NONE

    INTEGER(iwp) ::  i1  !< 
    INTEGER(iwp) ::  i2  !< 
    INTEGER(iwp) ::  j1  !< 
    INTEGER(iwp) ::  j2  !< 
    INTEGER(iwp) ::  k1  !< 
    INTEGER(iwp) ::  k2  !< 

    REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2) ::  array  !< 


!
!-- Open the file for data logging
    CALL check_open( 20 )

!
!-- Write the message string
    WRITE ( 20 )  log_message

!
!-- Write the simulated time and the array indices
    WRITE ( 20 )  simulated_time, i1, i2, j1, j2, k1, k2

!
!-- Write the array
    WRITE ( 20 )  array

 END SUBROUTINE data_log



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Complete logging of data for 2d arrays
!------------------------------------------------------------------------------!
 
 SUBROUTINE data_log_2d( array, i1, i2, j1, j2)

    USE control_parameters,                                                    &
        ONLY:  log_message, simulated_time

    USE kinds
            
    USE pegrid

    IMPLICIT NONE

    INTEGER(iwp) ::  i1  !< 
    INTEGER(iwp) ::  i2  !< 
    INTEGER(iwp) ::  j1  !< 
    INTEGER(iwp) ::  j2  !< 

    REAL(wp), DIMENSION(i1:i2,j1:j2) ::  array  !< 


!
!-- Open the file for data logging
    CALL check_open( 20 )

!
!-- Write the message string
    WRITE ( 20 )  log_message

!
!-- Write the simulated time and the array indices
    WRITE ( 20 )  simulated_time, i1, i2, j1, j2

!
!-- Write the array
    WRITE ( 20 )  array

 END SUBROUTINE data_log_2d



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Complete logging of data for 2d integer arrays
!------------------------------------------------------------------------------!
 
 SUBROUTINE data_log_2d_int( array, i1, i2, j1, j2)

    USE control_parameters,                                                    &
        ONLY:  log_message, simulated_time

    USE kinds
            
    USE pegrid

    IMPLICIT NONE

    INTEGER(iwp) ::  i1  !< 
    INTEGER(iwp) ::  i2  !< 
    INTEGER(iwp) ::  j1  !< 
    INTEGER(iwp) ::  j2  !< 

    INTEGER(iwp), DIMENSION(i1:i2,j1:j2) ::  array  !< 


!
!-- Open the file for data logging
    CALL check_open( 20 )

!
!-- Write the message string
    WRITE ( 20 )  log_message

!
!-- Write the simulated time and the array indices
    WRITE ( 20 )  simulated_time, i1, i2, j1, j2

!
!-- Write the array
    WRITE ( 20 )  array

 END SUBROUTINE data_log_2d_int
