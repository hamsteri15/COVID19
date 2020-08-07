!> @file user_init_grid.f90
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
! $Id: user_init_grid.f90 4360 2020-01-07 11:25:50Z suehring $
! Corrected "Former revisions" section
! 
! 3768 2019-02-27 14:35:58Z raasch
! variables commented + statement added to avoid compiler warnings about unused variables
! 
! 3655 2019-01-07 16:51:22Z knoop
! dz was replaced by dz(1)
!
! 217 2008-12-09 18:00:48Z letzel
! +topography_grid_convention
! Former file user_interface.f90 split into one file per subroutine
!
! Description:
! ------------
!> Execution of user-defined grid initializing actions
!------------------------------------------------------------------------------!
 SUBROUTINE user_init_grid( topo_3d )
 

    USE control_parameters
    
    USE indices
    
    USE kinds
    
    USE user

    IMPLICIT NONE

!    INTEGER(iwp)                                           ::  k_topo      !< topography top index
    INTEGER(iwp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  topo_3d     !< 3D topography field

!    REAL(wp) ::  h_topo !< user-defined topography height

!
!-- Next line is to avoid compiler warning about unused variables. Please remove.
    IF ( topo_3d(nzb,nysg,nxlg) == 0 )  CONTINUE


!
!-- Here the user-defined grid initializing actions follow:

!
!-- Set the index array nzb_local for non-flat topography.
!-- Here consistency checks concerning domain size and periodicity are necessary
    SELECT CASE ( TRIM( topography ) )

       CASE ( 'flat', 'single_building', 'single_street_canyon', 'tunnel' )
!
!--       Not allowed here since these are the standard cases used in init_grid.

       CASE ( 'user_defined_topography_1' )
!
!--       Here the user can define his own topography.
!--       After definition, please remove the following three lines!
          message_string = 'topography "' // topography // '" not available yet'
          CALL message( 'user_init_grid', 'UI0005', 1, 2, 0, 6, 0 )
!
!--       The user is allowed to set surface-mounted as well as non-surface
!--       mounted topography (e.g. overhanging structures). For both, use 
!--       3D array topo_3d and set bit 0. The convention is: bit is zero inside 
!--       topography, bit is 1 for atmospheric grid point. 
!--       The following example shows how to prescribe sine-like topography 
!--       along x-direction with amplitude of 10 * dz(1) and wavelength 10 * dy.
!           DO  i = nxlg, nxrg
!              h_topo = 10.0_wp * dz(1) * (SIN(3.14_wp*0.5_wp)*i*dx / ( 5.0_wp * dy ) )**2
! 
!              k_topo = MINLOC( ABS( zw - h_topo ), 1 ) - 1
! 
!              topo_3d(k_topo+1:nzt+1,:,i) =                                     &
!                                          IBSET( topo_3d(k_topo+1:nzt+1,:,i), 0 ) 
!           ENDDO 
! 
!           CALL exchange_horiz_int( topo_3d, nys, nyn, nxl, nxr, nzt, nbgp )

       CASE DEFAULT
!
!--       The DEFAULT case is reached if the parameter topography contains a
!--       wrong character string that is neither recognized in init_grid nor 
!--       here in user_init_grid.
          message_string = 'unknown topography "' // topography // '"'
          CALL message( 'user_init_grid', 'UI0006', 1, 2, 0, 6, 0 )

    END SELECT




 END SUBROUTINE user_init_grid

