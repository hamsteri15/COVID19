!> @file init_advec.f90
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
! $Id: init_advec.f90 4360 2020-01-07 11:25:50Z suehring $
! Corrected "Former revisions" section
! 
! 3655 2019-01-07 16:51:22Z knoop
! Corrected "Former revisions" section
!
! Revision 1.1  1999/02/05 09:07:38  raasch
! Initial revision
!
!
! Description:
! ------------
!> Initialize constant coefficients and parameters for certain advection schemes.
!------------------------------------------------------------------------------!
 SUBROUTINE init_advec
 

    USE advection,                                                             &
        ONLY:  aex, bex, dex, eex
        
    USE kinds
    
    USE control_parameters,                                                    &
        ONLY:  scalar_advec

    IMPLICIT NONE

    INTEGER(iwp) ::  i          !<
    INTEGER(iwp) ::  intervals  !< 
    INTEGER(iwp) ::  j          !<
    
    REAL(wp) :: delt   !<
    REAL(wp) :: dn     !<
    REAL(wp) :: dnneu  !<
    REAL(wp) :: ex1    !<
    REAL(wp) :: ex2    !<
    REAL(wp) :: ex3    !<
    REAL(wp) :: ex4    !<
    REAL(wp) :: ex5    !<
    REAL(wp) :: ex6    !<
    REAL(wp) :: sterm  !<


    IF ( scalar_advec == 'bc-scheme' )  THEN

!
!--    Compute exponential coefficients for the Bott-Chlond scheme
       intervals = 1000
       ALLOCATE( aex(intervals), bex(intervals), dex(intervals), eex(intervals) )

       delt  = 1.0_wp / REAL( intervals, KIND=wp )
       sterm = delt * 0.5_wp

       DO  i = 1, intervals

          IF ( sterm > 0.5_wp )  THEN
             dn = -5.0_wp
          ELSE
             dn = 5.0_wp
          ENDIF

          DO  j = 1, 15
             ex1 = dn * EXP( -dn ) - EXP( 0.5_wp * dn ) + EXP( -0.5_wp * dn )
             ex2 = EXP( dn ) - EXP( -dn )
             ex3 = EXP( -dn ) * ( 1.0_wp - dn ) - 0.5_wp * EXP(  0.5_wp * dn ) &
                                                - 0.5_wp * EXP( -0.5_wp * dn )
             ex4 = EXP( dn ) + EXP( -dn )
             ex5 = dn * sterm + ex1 / ex2
             ex6 = sterm + ( ex3 * ex2 - ex4 * ex1 ) / ( ex2 * ex2 )
             dnneu = dn - ex5 / ex6
             dn  = dnneu
          ENDDO

          IF ( sterm < 0.5_wp )  dn = MAX(  2.95E-2_wp, dn )
          IF ( sterm > 0.5_wp )  dn = MIN( -2.95E-2_wp, dn )
          ex1 = EXP( -dn )
          ex2 = EXP( dn ) - ex1
          aex(i) = -ex1 / ex2
          bex(i) = 1.0_wp / ex2
          dex(i) = dn
          eex(i) = EXP( dex(i) * 0.5_wp )
          sterm = sterm + delt

       ENDDO

    ENDIF

 END SUBROUTINE init_advec
