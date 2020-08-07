!> @file user_lpm_advec.f90
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
! $Id: user_lpm_advec.f90 4360 2020-01-07 11:25:50Z suehring $
! Corrected "Former revisions" section
! 
! 3768 2019-02-27 14:35:58Z raasch
! variables commented + statement added to avoid compiler warnings about unused variables
! 
! 3655 2019-01-07 16:51:22Z knoop
! Corrected "Former revisions" section
!
! 211 2008-11-11 04:46:24Z raasch
! Former file user_interface.f90 split into one file per subroutine
!
! Description:
! ------------
!> Modification of initial particles by the user.
!------------------------------------------------------------------------------!
 SUBROUTINE user_lpm_advec( ip, jp, kp )
 

    USE kinds
    
    USE particle_attributes
    
    USE user

    IMPLICIT NONE

    INTEGER(iwp) ::  ip   !< index of particle grid box, x-direction
    INTEGER(iwp) ::  jp   !< index of particle grid box, y-direction
    INTEGER(iwp) ::  kp   !< index of particle grid box, z-direction
!    INTEGER(iwp) ::  n    !< particle index
!    INTEGER(iwp) ::  nb   !< index of sub-box particles are sorted in

!    INTEGER(iwp), DIMENSION(0:7)  ::  start_index !< start particle index for current sub-box
!    INTEGER(iwp), DIMENSION(0:7)  ::  end_index   !< start particle index for current sub-box

!
!-- Next line is to avoid compiler warning about unused variables. Please remove.
    IF ( ip == 0  .OR.  jp == 0  .OR.  kp == 0 )  CONTINUE

!
!-- Here the user-defined actions follow
!   number_of_particles = prt_count(kp,jp,ip)
!   particles => grid_particles(kp,jp,ip)%particles(1:number_of_particles)
!
!   start_index = grid_particles(kp,jp,ip)%start_index
!   end_index   = grid_particles(kp,jp,ip)%end_index
!
!   IF ( number_of_particles <= 0 )  CYCLE
!   DO  n = 1, number_of_particles
!   DO  nb = 0, 7
!      DO  n = start_index(nb), end_index(nb)
!         particles(n)%xxx = 
!      ENDDO
!   ENDDO

 END SUBROUTINE user_lpm_advec

