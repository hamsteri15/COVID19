!> @posix_calls_from_fortran.f90
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
! $Id: posix_calls_from_fortran.f90 2696 2017-12-14 17:12:51Z kanani $
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! add variable description
! 
! 1986 2016-08-10 14:07:17Z gronemeier
! Initial revision
! 
! Description:
! ------------
!> Collection of POSIX-command calls for Fortran
!------------------------------------------------------------------------------!
 MODULE posix_calls_from_fortran

    USE, INTRINSIC ::  iso_c_binding, ONLY: c_int

    IMPLICIT none

    PRIVATE


    INTERFACE
!
!--    Sleep function from C library
       FUNCTION fsleep( seconds )  BIND( C, NAME='sleep' )
          IMPORT
          INTEGER(c_int) ::  fsleep
          INTEGER(c_int), INTENT(IN), VALUE ::  seconds
       END FUNCTION fsleep

    END INTERFACE

    INTERFACE fortran_sleep
       MODULE PROCEDURE fortran_sleep
    END INTERFACE fortran_sleep

    PUBLIC fortran_sleep


 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Wait a specified amount of seconds
!------------------------------------------------------------------------------!
 SUBROUTINE fortran_sleep( seconds )

    INTEGER, INTENT(IN) ::  seconds             !< seconds to wait

    INTEGER(c_int)      ::  seconds_in_c        !< same as seconds
    INTEGER(c_int)      ::  sleep_return_value  !< returned value to sleep

    seconds_in_c = seconds

    sleep_return_value = fsleep( seconds_in_c )

 END SUBROUTINE fortran_sleep

 END MODULE posix_calls_from_fortran
