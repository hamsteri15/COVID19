!> @file init_slope.f90
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
! $Id: init_slope.f90 4360 2020-01-07 11:25:50Z suehring $
! Corrected "Former revisions" section
! 
! 3655 2019-01-07 16:51:22Z knoop
! Modularization of all bulk cloud physics code components
!
! Revision 1.1  2000/04/27 07:06:24  raasch
! Initial revision
!
!
! Description:
! ------------
!> Initialization of the temperature field and other variables used in case
!> of a sloping surface.
!> @note when a sloping surface is used, only one constant temperature
!>       gradient is allowed!
!------------------------------------------------------------------------------!
 SUBROUTINE init_slope
 

    USE arrays_3d,                                                             &
        ONLY:  pt, pt_init, pt_slope_ref, zu
        
    USE basic_constants_and_equations_mod,                                     &
        ONLY:  pi
                    
    USE control_parameters,                                                    &
        ONLY:  alpha_surface, initializing_actions, pt_slope_offset,           &
               pt_surface, pt_vertical_gradient, sin_alpha_surface
        
    USE grid_variables,                                                        &
        ONLY:  dx
        
    USE indices,                                                               &
        ONLY:  ngp_2dh, nx, nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg, nzb, nzt
        
    USE kinds

    USE pegrid


    IMPLICIT NONE

    INTEGER(iwp) ::  i        !<
    INTEGER(iwp) ::  j        !<
    INTEGER(iwp) ::  k        !<
    
    REAL(wp)     ::  alpha    !<
    REAL(wp)     ::  height   !<
    REAL(wp)     ::  pt_value !<
    REAL(wp)     ::  radius   !<
    
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  pt_init_local !<

!
!-- Calculate reference temperature field needed for computing buoyancy
    ALLOCATE( pt_slope_ref(nzb:nzt+1,nxlg:nxrg) )

    DO  i = nxlg, nxrg
       DO  k = nzb, nzt+1

!
!--       Compute height of grid-point relative to lower left corner of
!--       the total domain.
!--       First compute the distance between the actual grid point and the
!--       lower left corner as well as the angle between the line connecting
!--       these points and the bottom of the model.
          IF ( k /= nzb )  THEN
             radius = SQRT( ( i * dx )**2 + zu(k)**2 )
             height = zu(k)
          ELSE
             radius = SQRT( ( i * dx )**2 )
             height = 0.0_wp
          ENDIF
          IF ( radius /= 0.0_wp )  THEN
             alpha = ASIN( height / radius )
          ELSE
             alpha = 0.0_wp
          ENDIF
!
!--       Compute temperatures in the rotated coordinate system
          alpha    = alpha + alpha_surface / 180.0_wp * pi
          pt_value = pt_surface + radius * SIN( alpha ) * &
                                  pt_vertical_gradient(1) / 100.0_wp
          pt_slope_ref(k,i) = pt_value
       ENDDO                
    ENDDO

!
!-- Temperature difference between left and right boundary of the total domain,
!-- used for the cyclic boundary in x-direction
    pt_slope_offset = (nx+1) * dx * sin_alpha_surface * &
                      pt_vertical_gradient(1) / 100.0_wp


!
!-- Following action must only be executed for initial runs
    IF ( TRIM( initializing_actions ) /= 'read_restart_data' )  THEN
!
!--    Set initial temperature equal to the reference temperature field
       DO  j = nysg, nyng
          pt(:,j,:) = pt_slope_ref
       ENDDO

!
!--    Recompute the mean initial temperature profile (mean along x-direction of
!--    the rotated coordinate system)
       ALLOCATE( pt_init_local(nzb:nzt+1) )
       pt_init_local = 0.0_wp
       DO  i = nxl, nxr
          DO  j =  nys, nyn
             DO  k = nzb, nzt+1
                pt_init_local(k) = pt_init_local(k) + pt(k,j,i)
             ENDDO
          ENDDO
       ENDDO

#if defined( __parallel )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( pt_init_local, pt_init, nzt+2-nzb, MPI_REAL, &
                            MPI_SUM, comm2d, ierr )
#else
       pt_init = pt_init_local
#endif

       pt_init = pt_init / ngp_2dh(0)
       DEALLOCATE( pt_init_local )

    ENDIF

 END SUBROUTINE init_slope
