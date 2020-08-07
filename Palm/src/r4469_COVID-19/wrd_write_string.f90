!> @file wrd_write_string.f90
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
! $Id: wrd_write_string.f90 4360 2020-01-07 11:25:50Z suehring $
! Initial revision
!
!
! Description:
! ------------
!> Calculates length of string and write the respective value together with the
!> string into binary file(s) for restart runs
!------------------------------------------------------------------------------!
 
 SUBROUTINE wrd_write_string( string )


    USE kinds


    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT (IN) ::  string  !< for storing strings in case of writing/reading restart data

    INTEGER(iwp) ::  length  !< integer that specifies the length of a string in case of writing/reading restart
                             !< data


    length = LEN_TRIM( string )

    WRITE ( 14 )  length
    WRITE ( 14 )  string(1:length)


 END SUBROUTINE wrd_write_string
