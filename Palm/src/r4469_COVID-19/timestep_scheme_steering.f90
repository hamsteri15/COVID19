!> @file timestep_scheme_steering.f90
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
! ------------------
! 
! 
! Former revisions:
! -----------------
! $Id: timestep_scheme_steering.f90 4360 2020-01-07 11:25:50Z suehring $
! Corrected "Former revisions" section
! 
! 3655 2019-01-07 16:51:22Z knoop
! OpenACC port for SPEC
!
! Revision 1.1  2004/01/28 15:34:47  raasch
! Initial revision
!
!
! Description:
! ------------
!> Depending on the timestep scheme set the steering factors for the prognostic
!> equations.
!------------------------------------------------------------------------------!
 SUBROUTINE timestep_scheme_steering
 

    USE control_parameters,                                                    &
        ONLY:  intermediate_timestep_count, timestep_scheme, tsc

    USE kinds

    IMPLICIT NONE


    IF ( timestep_scheme(1:5) == 'runge' )  THEN
!
!--    Runge-Kutta schemes (here the factors depend on the respective
!--    intermediate step)
       IF ( timestep_scheme == 'runge-kutta-2' )  THEN
          IF ( intermediate_timestep_count == 1 )  THEN
             tsc(1:5) = (/ 1.0_wp, 1.0_wp,  0.0_wp, 0.0_wp, 0.0_wp /)
          ELSE
             tsc(1:5) = (/ 1.0_wp, 0.5_wp, -0.5_wp, 0.0_wp, 1.0_wp /)
          ENDIF
       ELSE
          IF ( intermediate_timestep_count == 1 )  THEN
             tsc(1:5) = (/ 1.0_wp,  1.0_wp /  3.0_wp,           0.0_wp, 0.0_wp, 0.0_wp /)
          ELSEIF ( intermediate_timestep_count == 2 )  THEN
             tsc(1:5) = (/ 1.0_wp, 15.0_wp / 16.0_wp, -25.0_wp/48.0_wp, 0.0_wp, 0.0_wp /)
          ELSE
             tsc(1:5) = (/ 1.0_wp,  8.0_wp / 15.0_wp,   1.0_wp/15.0_wp, 0.0_wp, 1.0_wp /)
          ENDIF          
       ENDIF

       !$ACC UPDATE DEVICE(tsc(1:5))

    ELSEIF ( timestep_scheme == 'euler' )  THEN
!
!--    Euler scheme
       tsc(1:5) = (/ 1.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp /)

    ENDIF


 END SUBROUTINE timestep_scheme_steering
