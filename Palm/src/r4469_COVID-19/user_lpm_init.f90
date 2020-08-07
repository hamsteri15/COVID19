!> @file user_lpm_init.f90
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
! $Id: user_lpm_init.f90 4360 2020-01-07 11:25:50Z suehring $
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
!> Modification of initial particles by the user.
!------------------------------------------------------------------------------!
 SUBROUTINE user_lpm_init
 
    USE indices

    USE kinds
    
    USE particle_attributes
    
    USE user

    IMPLICIT NONE

!    INTEGER(iwp) ::  ip   !<
!    INTEGER(iwp) ::  jp   !<
!    INTEGER(iwp) ::  kp   !<
!    INTEGER(iwp) ::  n    !<

!
!-- Here the user-defined actions follow
!     DO  ip = nxl, nxr
!        DO  jp = nys, nyn
!           DO  kp = nzb+1, nzt
!              number_of_particles = prt_count(kp,jp,ip)
!              particles => grid_particles(kp,jp,ip)%particles(1:number_of_particles)
!              IF ( number_of_particles <= 0 )  CYCLE
!              DO  n = 1, number_of_particles
!
!              ENDDO
!           ENDDO
!        ENDDO
!     ENDDO

 END SUBROUTINE user_lpm_init

