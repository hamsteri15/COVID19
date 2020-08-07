!> @file time_to_string.f90
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
! $Id: time_to_string.f90 4360 2020-01-07 11:25:50Z suehring $
! Corrected "Former revisions" section
! 
! 3655 2019-01-07 16:51:22Z knoop
! Corrected "Former revisions" section
!
! Revision 1.1  1997/08/11 06:26:08  raasch
! Initial revision
!
!
! Description:
! ------------
!> Transforming the time from real to character-string hh:mm:ss
!------------------------------------------------------------------------------!
 FUNCTION time_to_string( time )
 

    USE kinds

    IMPLICIT NONE

    CHARACTER (LEN=9) ::  time_to_string !< 

    INTEGER(iwp)      ::  hours   !< 
    INTEGER(iwp)      ::  minutes !< 
    INTEGER(iwp)      ::  seconds !< 

    REAL(wp)          ::  rest_time !< 
    REAL(wp)          ::  time      !< 

!
!-- Calculate the number of hours, minutes, and seconds
    hours     = INT( time / 3600.0_wp )
    rest_time = time - hours * 3600_wp
    minutes   = INT( rest_time / 60.0_wp )
    seconds   = rest_time - minutes * 60

!
!-- Build the string
    IF ( hours < 100 )  THEN
       WRITE (time_to_string,'(I2.2,'':'',I2.2,'':'',I2.2)')  hours, minutes, &
                                                              seconds
    ELSE
       WRITE (time_to_string,'(I3.3,'':'',I2.2,'':'',I2.2)')  hours, minutes, &
                                                              seconds
    ENDIF

 END FUNCTION time_to_string
