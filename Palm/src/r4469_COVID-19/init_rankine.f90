!> @file init_rankine.f90
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
! $Id: init_rankine.f90 4457 2020-03-11 14:20:43Z raasch $
! use statement for exchange horiz added
! 
! 4360 2020-01-07 11:25:50Z suehring
! Corrected "Former revisions" section
! 
! 3655 2019-01-07 16:51:22Z knoop
! Modularization of all bulk cloud physics code components
!
! Revision 1.1  1997/08/11 06:18:43  raasch
! Initial revision
!
!
! Description:
! ------------
!> Initialize a (nondivergent) Rankine eddy with a vertical axis in order to test 
!> the advection terms and the pressure solver.
!------------------------------------------------------------------------------!
 SUBROUTINE init_rankine
 

    USE arrays_3d,                                                             &
        ONLY:  pt, pt_init, u, u_init, v, v_init

    USE control_parameters,                                                    &
        ONLY:  initializing_actions, n_sor, nsor, nsor_ini   

    USE basic_constants_and_equations_mod,                                     &
        ONLY:  pi

    USE exchange_horiz_mod,                                                    &
        ONLY:  exchange_horiz

    USE grid_variables,                                                        &
        ONLY:  dx, dy 

    USE indices,                                                               &
        ONLY:  nbgp, nx, nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg, nzb, nzt     
                
    USE kinds

    IMPLICIT NONE

    INTEGER(iwp) ::  i   !<
    INTEGER(iwp) ::  ic  !<
    INTEGER(iwp) ::  j   !<
    INTEGER(iwp) ::  jc  !<
    INTEGER(iwp) ::  k   !<
    INTEGER(iwp) ::  kc1 !<
    INTEGER(iwp) ::  kc2 !<
    
    REAL(wp)     ::  alpha  !<
    REAL(wp)     ::  betrag !<
    REAL(wp)     ::  radius !<
    REAL(wp)     ::  rc     !<
    REAL(wp)     ::  uw     !<
    REAL(wp)     ::  vw     !<
    REAL(wp)     ::  x      !<
    REAL(wp)     ::  y      !<

!
!-- Default: eddy radius rc, eddy strength z,
!--          position of eddy centre: ic, jc, kc1, kc2
    rc  =  4.0_wp * dx
    ic  =  ( nx+1 ) / 2
    jc  =  ic
    kc1 = nzb
    kc2 = nzt+1

!
!-- Reset initial profiles to constant profiles
    IF ( INDEX(initializing_actions, 'set_constant_profiles') /= 0 )  THEN
       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
             pt(:,j,i) = pt_init
             u(:,j,i)  = u_init
             v(:,j,i)  = v_init
          ENDDO
       ENDDO
    ENDIF

!
!-- Compute the u-component.
    DO  i = nxl, nxr
       DO  j = nys, nyn
          x = ( i - ic - 0.5_wp ) * dx
          y = ( j - jc          ) * dy
          radius = SQRT( x**2 + y**2 )
          IF ( radius <= 2.0_wp * rc )  THEN
             betrag = radius / ( 2.0_wp * rc ) * 0.08_wp
          ELSEIF ( radius > 2.0_wp * rc  .AND.  radius < 8.0_wp * rc )  THEN
             betrag = 0.08_wp * EXP( -( radius - 2.0_wp * rc ) / 2.0_wp )
          ELSE
             betrag = 0.0_wp
          ENDIF

          IF ( x == 0.0_wp )  THEN
             IF ( y > 0.0_wp )  THEN
                alpha = pi / 2.0_wp
             ELSEIF ( y < 0.0_wp )  THEN
                alpha = 3.0_wp * pi / 2.0_wp
             ENDIF
          ELSE
             IF ( x < 0.0_wp )  THEN
                alpha = ATAN( y / x ) + pi
             ELSE
                IF ( y < 0.0_wp )  THEN
                   alpha = ATAN( y / x ) + 2.0_wp * pi
                ELSE
                   alpha = ATAN( y / x )
                ENDIF
             ENDIF
          ENDIF

          uw = -SIN( alpha ) * betrag

          DO  k = kc1, kc2
             u(k,j,i) = u(k,j,i) + uw
          ENDDO
       ENDDO
    ENDDO

!
!-- Compute the v-component.
    DO  i = nxl, nxr
       DO  j = nys, nyn
          x = ( i - ic          ) * dx
          y = ( j - jc - 0.5_wp ) * dy
          radius = SQRT( x**2 + y**2 )
          IF ( radius <= 2.0_wp * rc )  THEN
             betrag = radius / ( 2.0_wp * rc ) * 0.08_wp
          ELSEIF ( radius > 2.0_wp * rc  .AND.  radius < 8.0_wp * rc )  THEN
             betrag = 0.08_wp * EXP( -( radius - 2.0_wp * rc ) / 2.0_wp )
          ELSE
             betrag = 0.0_wp
          ENDIF

          IF ( x == 0.0_wp )  THEN
             IF ( y > 0.0_wp )  THEN
                alpha = pi / 2.0_wp
             ELSEIF ( y < 0.0_wp )  THEN
                alpha = 3.0_wp * pi / 2.0_wp
             ENDIF
          ELSE
             IF ( x < 0.0_wp )  THEN
                alpha = ATAN( y / x ) + pi
             ELSE
                IF ( y < 0.0_wp )  THEN
                   alpha = ATAN( y / x ) + 2.0_wp * pi
                ELSE
                   alpha = ATAN( y / x )
                ENDIF
             ENDIF
          ENDIF

          vw = COS( alpha ) * betrag

          DO  k = kc1, kc2
             v(k,j,i) = v(k,j,i) + vw
          ENDDO
       ENDDO
    ENDDO

!
!-- Exchange of boundary values for the velocities.
    CALL exchange_horiz( u, nbgp)
    CALL exchange_horiz( v, nbgp )
!
!-- Make velocity field nondivergent.
    n_sor = nsor_ini
    CALL pres
    n_sor = nsor

 END SUBROUTINE init_rankine
