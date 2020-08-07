!> @file local_tremain_ini.f90
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
! $Id: local_tremain_ini.f90 4360 2020-01-07 11:25:50Z suehring $
! Corrected "Former revisions" section
! 
! 3655 2019-01-07 16:51:22Z knoop
! Corrected "Former revisions" section
!
! Revision 1.1  1998/03/18 20:15:05  raasch
! Initial revision
!
!
! Description:
! ------------
!> Initialization of CPU-time measurements for different operating systems
!------------------------------------------------------------------------------!
 SUBROUTINE local_tremain_ini
 
      
    USE cpulog,                                                                &
        ONLY:  initial_wallclock_time
        
    USE kinds

    IMPLICIT NONE

    INTEGER(idp)     ::  count      !<
    INTEGER(idp)     ::  count_rate !<

!
!-- Get initial wall clock time
    CALL SYSTEM_CLOCK( count, count_rate )
    initial_wallclock_time = REAL( count, KIND=wp ) / REAL( count_rate, KIND=wp )

 END SUBROUTINE local_tremain_ini
