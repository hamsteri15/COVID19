!> @file random_function_mod.f90
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
! $Id: random_function_mod.f90 4360 2020-01-07 11:25:50Z suehring $
! Corrected "Former revisions" section
! 
! 4144 2019-08-06 09:11:47Z raasch
! relational operators .EQ., .NE., etc. replaced by ==, /=, etc.
! 
! 3802 2019-03-17 13:33:42Z raasch
! type conversion added to avoid compiler warning about constant integer
! division truncation
! 
! 3655 2019-01-07 16:51:22Z knoop
! Corrected "Former revisions" section
!
! Revision 1.1  1998/02/04 16:09:45  raasch
! Initial revision
!
!
! Description:
! ------------
!> Random number generator, produces numbers equally distributed in interval [0,1]
!> This routine is taken from the "numerical recipes"
!------------------------------------------------------------------------------!
 MODULE random_function_mod
 

    USE kinds

    IMPLICIT NONE

    PRIVATE

    PUBLIC random_function, random_function_ini

    INTEGER(iwp), PUBLIC, SAVE ::  random_iv(32)  !<
    INTEGER(iwp), PUBLIC, SAVE ::  random_iy      !<

    INTERFACE random_function_ini
       MODULE PROCEDURE random_function_ini
    END INTERFACE random_function_ini

    INTERFACE random_function
       MODULE PROCEDURE random_function
    END INTERFACE random_function

 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE random_function_ini

       IMPLICIT NONE

       random_iv = 0
       random_iy = 0

    END SUBROUTINE random_function_ini

    FUNCTION random_function( idum )


       IMPLICIT NONE

       INTEGER(iwp) ::  ia               !<
       INTEGER(iwp) ::  idum             !<
       INTEGER(iwp) ::  im               !<
       INTEGER(iwp) ::  iq               !<
       INTEGER(iwp) ::  ir               !<
       INTEGER(iwp) ::  ndiv             !<
       INTEGER(iwp) ::  ntab             !<

       INTEGER(iwp) ::  j                !<
       INTEGER(iwp) ::  k                !<

       REAL(wp)     ::  am               !<
       REAL(wp)     ::  eps              !<
       REAL(wp)     ::  random_function  !<
       REAL(wp)     ::  rnmx             !<

       PARAMETER ( ia=16807, im=2147483647, am=1.0_wp/im, iq=127773, ir=2836, &
                   ntab=32, ndiv=1+INT(REAL(im-1)/ntab), eps=1.2e-7_wp, rnmx=1.0_wp-eps )

       IF ( idum <= 0  .OR.  random_iy == 0 )  THEN
          idum = max (-idum,1)
          DO  j = ntab+8,1,-1
             k    = idum / iq
             idum = ia * ( idum - k * iq ) - ir * k
             IF ( idum < 0 )  idum = idum + im
             IF ( j <= ntab )  random_iv(j) = idum
          ENDDO
          random_iy = random_iv(1)
       ENDIF

       k    = idum / iq
       idum = ia * ( idum - k * iq ) - ir * k
       IF ( idum < 0 )  idum = idum + im
       j            = 1 + random_iy / ndiv
       random_iy    = random_iv(j)
       random_iv(j) = idum
       random_function  = MIN( am * random_iy , rnmx )

    END FUNCTION random_function

 END MODULE random_function_mod
