!> @file init_grid.f90
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
! $Id: init_grid.f90 4457 2020-03-11 14:20:43Z raasch $
! use statement for exchange horiz added,
! bugfix for call of exchange horiz 2d
! 
! 4444 2020-03-05 15:59:50Z raasch
! bugfix: cpp-directives for serial mode added
! 
! 4414 2020-02-19 20:16:04Z suehring
! - Remove deprecated topography arrays nzb_s_inner, nzb_u_inner, etc.
! - Move initialization of boundary conditions and multigrid into an extra
!   module interface.
! 
! 4386 2020-01-27 15:07:30Z Giersch
! Allocation statements, comments, naming of variables revised and _wp added to 
! real type values
! 
! 4360 2020-01-07 11:25:50Z suehring
! Revise error messages for generic tunnel setup.
! 
! 4346 2019-12-18 11:55:56Z motisi
! Introduction of wall_flags_total_0, which currently sets bits based on static
! topography information used in wall_flags_static_0
! 
! 4340 2019-12-16 08:17:03Z Giersch
! Topography closed channel flow with symmetric boundaries implemented
! 
! 4329 2019-12-10 15:46:36Z motisi
! Renamed wall_flags_0 to wall_flags_static_0
! 
! 4328 2019-12-09 18:53:04Z suehring
! Minor change in nzb_max computation. Commentation added.
! 
! 4314 2019-11-29 10:29:20Z suehring
! Set additional topography flag 4 to mark topography grid points emerged
! from the filtering process.
! 
! 4294 2019-11-13 18:34:16Z suehring
! Bugfix, always set bit 5 and 6 of wall_flags, indicating terrain- and 
! building surfaces in all  cases, in order to enable terrain-following output 
! also when no land- or urban-surface model is applied. 
! 
! 4265 2019-10-15 16:16:24Z suehring
! Bugfix for last commit, exchange oro_max variable only when it is allocated 
! (not necessarily the case when topography is input from ASCII file).
! 
! 4245 2019-09-30 08:40:37Z pavelkrc
! Store oro_max (building z-offset) in 2D for building surfaces
! 
! 4189 2019-08-26 16:19:38Z suehring
! - Add check for proper setting of namelist parameter topography
! - Set flag to indicate land surfaces in case no topography is provided
! 
! 4182 2019-08-22 15:20:23Z scharf
! Corrected "Former revisions" section
! 
! 4168 2019-08-16 13:50:17Z suehring
! Pre-calculate topography top index and store it on an array (replaces former
! functions get_topography_top_index)
! 
! 4159 2019-08-15 13:31:35Z suehring
! Revision of topography processing. This was not consistent between 2D and 3D
! buildings. 
! 
! 4144 2019-08-06 09:11:47Z raasch
! relational operators .EQ., .NE., etc. replaced by ==, /=, etc.
! 
! 4115 2019-07-24 12:50:49Z suehring
! Bugfix in setting near-surface flag 24, inidicating wall-bounded grid points 
! 
! 4110 2019-07-22 17:05:21Z suehring
! - Separate initialization of advection flags for momentum and scalars.
! - Change subroutine interface for ws_init_flags_scalar to pass boundary flags
! 
! 4109 2019-07-22 17:00:34Z suehring
! Fix bad commit
! 
! 3926 2019-04-23 12:56:42Z suehring
! Minor bugfix in building mapping when all building IDs in the model domain
! are missing
! 
! 3857 2019-04-03 13:00:16Z knoop
! In projection of non-building 3D objects onto numerical grid remove 
! dependency on building_type
! 
! 3763 2019-02-25 17:33:49Z suehring
! Replace work-around for ghost point exchange of 1-byte arrays with specific 
! routine as already done in other routines
! 
! 3761 2019-02-25 15:31:42Z raasch
! unused variables removed
! 
! 3661 2019-01-08 18:22:50Z suehring
! Remove setting of nzb_max to nzt at non-cyclic boundary PEs, instead, 
! order degradation of advection scheme is handeled directly in advec_ws
! 
! 3655 2019-01-07 16:51:22Z knoop
! Comment added
!
! Revision 1.1  1997/08/11 06:17:45  raasch
! Initial revision (Testversion)
!
!
! Description:
! -----------------------------------------------------------------------------!
!> Creating grid depending constants
!> @todo: Rearrange topo flag list
!> @todo: reference 3D buildings on top of orography is not tested and may need
!>        further improvement for steep slopes
!> @todo: Use more advanced setting of building type at filled holes 
!------------------------------------------------------------------------------!
 SUBROUTINE init_grid

    USE arrays_3d,                                                             &
        ONLY:  dd2zu, ddzu, ddzu_pres, ddzw, dzu, dzw, x, xu, y, yv, zu, zw

    USE control_parameters,                                                    &
        ONLY:  constant_flux_layer, dz, dz_max, dz_stretch_factor,             &
               dz_stretch_factor_array, dz_stretch_level, dz_stretch_level_end,&
               dz_stretch_level_end_index, dz_stretch_level_start_index,       &
               dz_stretch_level_start, ibc_uv_b, message_string,               &
               number_stretch_level_end,                                       &
               number_stretch_level_start,                                     &
               ocean_mode,                                                     &
               psolver,                                                        &
               symmetry_flag,                                                  &
               topography,                                                     &
               use_surface_fluxes

    USE grid_variables,                                                        &
        ONLY:  ddx, ddx2, ddy, ddy2, dx, dx2, dy, dy2, zu_s_inner, zw_w_inner

    USE indices,                                                               &
        ONLY:  nbgp,                                                           &
               nx,                                                             &
               nxl,                                                            &
               nxlg,                                                           &
               nxr,                                                            &
               nxrg,                                                           &
               ny,                                                             &
               nyn,                                                            &
               nyng,                                                           &
               nys,                                                            &
               nysg,                                                           &
               nz,                                                             &
               nzb,                                                            &
               nzb_diff,                                                       &
               nzb_max,                                                        &
               nzt,                                                            &
               topo_top_ind,                                                   &
               topo_min_level

    USE kinds

    USE pegrid

#if defined( __parallel )
    USE vertical_nesting_mod,                                                  &
        ONLY:  vnested, vnest_init_grid
#endif

    IMPLICIT NONE

    INTEGER(iwp) ::  i             !< index variable along x 
    INTEGER(iwp) ::  j             !< index variable along y
    INTEGER(iwp) ::  k             !< index variable along z
    INTEGER(iwp) ::  k_top         !< topography top index on local PE
    INTEGER(iwp) ::  n             !< loop variable for stretching
    INTEGER(iwp) ::  number_dz     !< number of user-specified dz values       
    INTEGER(iwp) ::  nzb_local_max !< vertical grid index of maximum topography height
    INTEGER(iwp) ::  nzb_local_min !< vertical grid index of minimum topography height

    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  topo !< input array for 3D topography and dummy array for setting "outer"-flags

    REAL(wp) ::  dz_level_end  !< distance between calculated height level for u/v-grid and user-specified end level for stretching
    REAL(wp) ::  dz_stretched  !< stretched vertical grid spacing
    
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  min_dz_stretch_level_end !< Array that contains all minimum heights where the stretching can end


!
!-- Calculation of horizontal array bounds including ghost layers
    nxlg = nxl - nbgp
    nxrg = nxr + nbgp
    nysg = nys - nbgp
    nyng = nyn + nbgp

!
!-- Allocate grid arrays
    ALLOCATE( x(0:nx) )
    ALLOCATE( xu(0:nx) )
    
    DO i = 0, nx
       xu(i) = i * dx
       x(i)  = i * dx + 0.5_wp * dx
    ENDDO

    ALLOCATE( y(0:ny) )
    ALLOCATE( yv(0:ny) )
    
    DO j = 0, ny
       yv(j) = j * dy
       y(j)  = j * dy + 0.5_wp * dy
    ENDDO

    ALLOCATE( ddzu(1:nzt+1) )
    ALLOCATE( ddzw(1:nzt+1) )
    ALLOCATE( dd2zu(1:nzt) )
    ALLOCATE( dzu(1:nzt+1) )
    ALLOCATE( dzw(1:nzt+1) )
    ALLOCATE( zu(nzb:nzt+1) )
    ALLOCATE( zw(nzb:nzt+1) )

!
!-- For constructing an appropriate grid, the vertical grid spacing dz has to
!-- be specified with a non-negative value in the parameter file
    IF ( dz(1) == -1.0_wp )  THEN
       message_string = 'missing dz'
       CALL message( 'init_grid', 'PA0200', 1, 2, 0, 6, 0 ) 
    ELSEIF ( dz(1) <= 0.0_wp )  THEN
       WRITE( message_string, * ) 'dz=',dz(1),' <= 0.0'
       CALL message( 'init_grid', 'PA0201', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Initialize dz_stretch_level_start with the value of dz_stretch_level
!-- if it was set by the user
    IF ( dz_stretch_level /= -9999999.9_wp ) THEN
       dz_stretch_level_start(1) = dz_stretch_level
    ENDIF
       
!
!-- Determine number of dz values and stretching levels specified by the
!-- user to allow right controlling of the stretching mechanism and to
!-- perform error checks. The additional requirement that dz /= dz_max
!-- for counting number of user-specified dz values is necessary. Otherwise 
!-- restarts would abort if the old stretching mechanism with dz_stretch_level
!-- is used (Attention: The user is not allowed to specify a dz value equal
!-- to the default of dz_max = 999.0). 
    number_dz = COUNT( dz /= -1.0_wp .AND. dz /= dz_max)
    number_stretch_level_start = COUNT( dz_stretch_level_start /=              &
                                       -9999999.9_wp )
    number_stretch_level_end = COUNT( dz_stretch_level_end /=                  &
                                      9999999.9_wp )

!
!-- The number of specified end levels +1 has to be the same as the number 
!-- of specified dz values
    IF ( number_dz /= number_stretch_level_end + 1 ) THEN
       WRITE( message_string, * ) 'The number of values for dz = ',            &
                                   number_dz, 'has to be the same as& ',       &
                                   'the number of values for ',                &
                                   'dz_stretch_level_end + 1 = ',              &
                                   number_stretch_level_end+1
          CALL message( 'init_grid', 'PA0156', 1, 2, 0, 6, 0 )
    ENDIF
    
!
!-- The number of specified start levels has to be the same or one less than 
!-- the number of specified dz values
    IF ( number_dz /= number_stretch_level_start + 1 .AND.                     &
         number_dz /= number_stretch_level_start ) THEN
       WRITE( message_string, * ) 'The number of values for dz = ',            &
                                   number_dz, 'has to be the same as or one ', &
                                   'more than& the number of values for ',     &
                                   'dz_stretch_level_start = ',                &
                                   number_stretch_level_start
          CALL message( 'init_grid', 'PA0211', 1, 2, 0, 6, 0 )
    ENDIF
    
!-- The number of specified start levels has to be the same or one more than 
!-- the number of specified end levels
    IF ( number_stretch_level_start /= number_stretch_level_end + 1 .AND.      &
         number_stretch_level_start /= number_stretch_level_end ) THEN
       WRITE( message_string, * ) 'The number of values for ',                 &
                                  'dz_stretch_level_start = ',                 &
                                   dz_stretch_level_start, 'has to be the ',   &
                                   'same or one more than& the number of ',    &
                                   'values for dz_stretch_level_end = ',       &
                                   number_stretch_level_end
          CALL message( 'init_grid', 'PA0216', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Initialize dz for the free atmosphere with the value of dz_max
    IF ( dz(number_stretch_level_start+1) == -1.0_wp .AND.                     &
         number_stretch_level_start /= 0 ) THEN 
       dz(number_stretch_level_start+1) = dz_max
    ENDIF
       
!
!-- Initialize the stretching factor if (infinitely) stretching in the free 
!-- atmosphere is desired (dz_stretch_level_end was not specified for the 
!-- free atmosphere)
    IF ( number_stretch_level_start == number_stretch_level_end + 1 ) THEN 
       dz_stretch_factor_array(number_stretch_level_start) =                   &
       dz_stretch_factor
    ENDIF
    
!
!-- Allocation of arrays for stretching
    ALLOCATE( min_dz_stretch_level_end(number_stretch_level_start) )

!
!-- Define the vertical grid levels. Start with atmosphere branch
    IF ( .NOT. ocean_mode )  THEN
    
!
!--    The stretching region has to be large enough to allow for a smooth
!--    transition between two different grid spacings. The number 4 is an
!--    empirical value
       DO n = 1, number_stretch_level_start
          min_dz_stretch_level_end(n) = dz_stretch_level_start(n) +            &
                                        4 * MAX( dz(n),dz(n+1) )
       ENDDO

       IF ( ANY( min_dz_stretch_level_end(1:number_stretch_level_start) >      &
                 dz_stretch_level_end(1:number_stretch_level_start) ) ) THEN
             message_string= 'Each dz_stretch_level_end has to be larger ' //  &
                             'than its corresponding value for &' //           &
                             'dz_stretch_level_start + 4*MAX(dz(n),dz(n+1)) '//&
                             'to allow for smooth grid stretching'
             CALL message( 'init_grid', 'PA0224', 1, 2, 0, 6, 0 )
       ENDIF
       
!
!--    Stretching must not be applied within the surface layer 
!--    (first two grid points). For the default case dz_stretch_level_start 
!--    is negative. Therefore the absolut value is checked here.
       IF ( ANY( ABS( dz_stretch_level_start ) <= dz(1) * 1.5_wp ) ) THEN
          WRITE( message_string, * ) 'Each dz_stretch_level_start has to be ',&
                                     'larger than ', dz(1) * 1.5
             CALL message( 'init_grid', 'PA0226', 1, 2, 0, 6, 0 )
       ENDIF

!
!--    The stretching has to start and end on a grid level. Therefore 
!--    user-specified values are mapped to the next lowest level. The  
!--    calculation of the first level is realized differently just because of 
!--    historical reasons (the advanced/new stretching mechanism was realized  
!--    in a way that results don't change if the old parameters 
!--    dz_stretch_level, dz_stretch_factor and dz_max are used)
       IF ( number_stretch_level_start /= 0 ) THEN
          dz_stretch_level_start(1) = INT( (dz_stretch_level_start(1) -        &
                                            dz(1)/2.0) / dz(1) )               &
                                      * dz(1) + dz(1)/2.0
       ENDIF
       
       IF ( number_stretch_level_start > 1 ) THEN
          DO n = 2, number_stretch_level_start
             dz_stretch_level_start(n) = INT( dz_stretch_level_start(n) /      &
                                              dz(n) ) * dz(n)
          ENDDO
       ENDIF
       
       IF ( number_stretch_level_end /= 0 ) THEN
          DO n = 1, number_stretch_level_end
             dz_stretch_level_end(n) = INT( dz_stretch_level_end(n) /          &
                                            dz(n+1) ) * dz(n+1)
          ENDDO
       ENDIF

!
!--    Determine stretching factor if necessary
       IF ( number_stretch_level_end >= 1 ) THEN 
          CALL calculate_stretching_factor( number_stretch_level_end )
       ENDIF

!
!--    Grid for atmosphere with surface at z=0 (k=0, w-grid).
!--    First compute the u- and v-levels. In case of dirichlet bc for u and v
!--    the first u/v- and w-level (k=0) are defined at same height (z=0). 
!--    The second u-level (k=1) corresponds to the top of the 
!--    surface layer. In case of symmetric boundaries (closed channel flow), 
!--    the first grid point is always at z=0.
       IF ( ibc_uv_b == 0 .OR. ibc_uv_b == 2 .OR.                              & 
            topography == 'closed_channel' ) THEN
          zu(0) = 0.0_wp
       ELSE
          zu(0) = - dz(1) * 0.5_wp
       ENDIF
          
       zu(1) =   dz(1) * 0.5_wp
       
!
!--    Determine u and v height levels considering the possibility of grid
!--    stretching in several heights.
       n = 1
       dz_stretch_level_start_index = nzt+1
       dz_stretch_level_end_index = nzt+1
       dz_stretched = dz(1)

!--    The default value of dz_stretch_level_start is negative, thus the first
!--    condition is true even if no stretching shall be applied. Hence, the 
!--    second condition is also necessary.
       DO  k = 2, nzt+1-symmetry_flag
          IF ( dz_stretch_level_start(n) <= zu(k-1) .AND.                      &
               dz_stretch_level_start(n) /= -9999999.9_wp ) THEN
             dz_stretched = dz_stretched * dz_stretch_factor_array(n)
             
             IF ( dz(n) > dz(n+1) ) THEN
                dz_stretched = MAX( dz_stretched, dz(n+1) ) !Restrict dz_stretched to the user-specified (higher) dz
             ELSE
                dz_stretched = MIN( dz_stretched, dz(n+1) ) !Restrict dz_stretched to the user-specified (lower) dz
             ENDIF
             
             IF ( dz_stretch_level_start_index(n) == nzt+1 )                         &
             dz_stretch_level_start_index(n) = k-1
             
          ENDIF
          
          zu(k) = zu(k-1) + dz_stretched
          
!
!--       Make sure that the stretching ends exactly at dz_stretch_level_end 
          dz_level_end = ABS( zu(k) - dz_stretch_level_end(n) ) 
          
          IF ( dz_level_end  < dz(n+1)/3.0 ) THEN
             zu(k) = dz_stretch_level_end(n)
             dz_stretched = dz(n+1)
             dz_stretch_level_end_index(n) = k
             n = n + 1             
          ENDIF
       ENDDO
       
!
!--    If a closed channel flow is simulated, make sure that grid structure is  
!--    the same for both bottom and top boundary. (Hint: Using a different dz 
!--    at the bottom and at the top makes no sense due to symmetric boundaries
!--    where dz should be equal. Therefore, different dz at the bottom and top  
!--    causes an abort (see check_parameters).)
       IF ( topography == 'closed_channel' ) THEN
          zu(nzt+1) = zu(nzt) + dz(1) * 0.5_wp
       ENDIF

!
!--    Compute the w-levels. They are always staggered half-way between the 
!--    corresponding u-levels. In case of dirichlet bc for u and v at the 
!--    ground the first u- and w-level (k=0) are defined at same height (z=0). 
!--    Per default, the top w-level is extrapolated linearly. In case of 
!--    a closed channel flow, zu(nzt+1) and zw(nzt) must be set explicitely.
!--    (Hint: Using a different dz at the bottom and at the top makes no sense
!--    due to symmetric boundaries where dz should be equal. Therefore, 
!--    different dz at the bottom and top causes an abort (see 
!--    check_parameters).)
       zw(0) = 0.0_wp
       DO  k = 1, nzt-symmetry_flag
          zw(k) = ( zu(k) + zu(k+1) ) * 0.5_wp
       ENDDO
       IF ( topography == 'closed_channel' ) THEN
          zw(nzt)   = zw(nzt-1) + dz(1)
          zw(nzt+1) = zw(nzt) + dz(1)
       ELSE
          zw(nzt+1) = zw(nzt) + 2.0_wp * ( zu(nzt+1) - zw(nzt) )
       ENDIF

    ELSE !ocean branch

!
!--    The stretching region has to be large enough to allow for a smooth
!--    transition between two different grid spacings. The number 4 is an
!--    empirical value
       DO n = 1, number_stretch_level_start
          min_dz_stretch_level_end(n) = dz_stretch_level_start(n) -            &
                                        4 * MAX( dz(n),dz(n+1) )
       ENDDO
       
       IF ( ANY( min_dz_stretch_level_end (1:number_stretch_level_start) <     &
                 dz_stretch_level_end(1:number_stretch_level_start) ) ) THEN
             message_string= 'Each dz_stretch_level_end has to be less ' //   &
                             'than its corresponding value for &' //           &
                             'dz_stretch_level_start - 4*MAX(dz(n),dz(n+1)) '//&
                             'to allow for smooth grid stretching'
             CALL message( 'init_grid', 'PA0224', 1, 2, 0, 6, 0 )
       ENDIF
       
!
!--    Stretching must not be applied close to the surface (last two grid 
!--    points). For the default case dz_stretch_level_start is negative. 
       IF ( ANY( dz_stretch_level_start >= - dz(1) * 1.5_wp ) ) THEN
          WRITE( message_string, * ) 'Each dz_stretch_level_start has to be ',&
                                     'less than ', -dz(1) * 1.5
             CALL message( 'init_grid', 'PA0226', 1, 2, 0, 6, 0 )
       ENDIF

!
!--    The stretching has to start and end on a grid level. Therefore 
!--    user-specified values are mapped to the next highest level. The  
!--    calculation of the first level is realized differently just because of 
!--    historical reasons (the advanced/new stretching mechanism was realized  
!--    in a way that results don't change if the old parameters 
!--    dz_stretch_level, dz_stretch_factor and dz_max are used)
       IF ( number_stretch_level_start /= 0 ) THEN
          dz_stretch_level_start(1) = INT( (dz_stretch_level_start(1) +        &
                                            dz(1)/2.0) / dz(1) )               &
                                      * dz(1) - dz(1)/2.0
       ENDIF
       
       IF ( number_stretch_level_start > 1 ) THEN
          DO n = 2, number_stretch_level_start
             dz_stretch_level_start(n) = INT( dz_stretch_level_start(n) /      &
                                              dz(n) ) * dz(n)
          ENDDO
       ENDIF
       
       IF ( number_stretch_level_end /= 0 ) THEN
          DO n = 1, number_stretch_level_end
             dz_stretch_level_end(n) = INT( dz_stretch_level_end(n) /          &
                                            dz(n+1) ) * dz(n+1)
          ENDDO
       ENDIF
       
!
!--    Determine stretching factor if necessary
       IF ( number_stretch_level_end >= 1 ) THEN 
          CALL calculate_stretching_factor( number_stretch_level_end )
       ENDIF

!
!--    Grid for ocean with free water surface is at k=nzt (w-grid).
!--    In case of neumann bc at the ground the first first u-level (k=0) lies
!--    below the first w-level (k=0). In case of dirichlet bc the first u- and
!--    w-level are defined at same height, but staggered from the second level.
!--    The second u-level (k=1) corresponds to the top of the surface layer.
!--    z values are negative starting from z=0 (surface)
       zu(nzt+1) =   dz(1) * 0.5_wp
       zu(nzt)   = - dz(1) * 0.5_wp

!
!--    Determine u and v height levels considering the possibility of grid
!--    stretching in several heights.
       n = 1
       dz_stretch_level_start_index = 0
       dz_stretch_level_end_index = 0
       dz_stretched = dz(1)

       DO  k = nzt-1, 0, -1
          
          IF ( dz_stretch_level_start(n) >= zu(k+1) ) THEN
             dz_stretched = dz_stretched * dz_stretch_factor_array(n)

             IF ( dz(n) > dz(n+1) ) THEN
                dz_stretched = MAX( dz_stretched, dz(n+1) ) !Restrict dz_stretched to the user-specified (higher) dz
             ELSE
                dz_stretched = MIN( dz_stretched, dz(n+1) ) !Restrict dz_stretched to the user-specified (lower) dz
             ENDIF
             
             IF ( dz_stretch_level_start_index(n) == 0 )                             &
             dz_stretch_level_start_index(n) = k+1
             
          ENDIF
          
          zu(k) = zu(k+1) - dz_stretched
          
!
!--       Make sure that the stretching ends exactly at dz_stretch_level_end 
          dz_level_end = ABS( zu(k) - dz_stretch_level_end(n) ) 
          
          IF ( dz_level_end  < dz(n+1)/3.0 ) THEN
             zu(k) = dz_stretch_level_end(n)
             dz_stretched = dz(n+1)
             dz_stretch_level_end_index(n) = k
             n = n + 1             
          ENDIF
       ENDDO
       
!
!--    Compute the w-levels. They are always staggered half-way between the 
!--    corresponding u-levels, except in case of dirichlet bc for u and v
!--    at the ground. In this case the first u- and w-level are defined at 
!--    same height. The top w-level (nzt+1) is not used but set for 
!--    consistency, since w and all scalar variables are defined up tp nzt+1.
       zw(nzt+1) = dz(1)
       zw(nzt)   = 0.0_wp
       DO  k = 0, nzt
          zw(k) = ( zu(k) + zu(k+1) ) * 0.5_wp
       ENDDO

!
!--    In case of dirichlet bc for u and v the first u- and w-level are defined
!--    at same height.
       IF ( ibc_uv_b == 0 ) THEN
          zu(0) = zw(0)
       ENDIF

    ENDIF !End of defining the vertical grid levels

!
!-- Compute grid lengths.
    DO  k = 1, nzt+1
       dzu(k)  = zu(k) - zu(k-1)
       ddzu(k) = 1.0_wp / dzu(k)
       dzw(k)  = zw(k) - zw(k-1)
       ddzw(k) = 1.0_wp / dzw(k)
    ENDDO

    DO  k = 1, nzt
       dd2zu(k) = 1.0_wp / ( dzu(k) + dzu(k+1) )
    ENDDO
    
!    
!-- The FFT- SOR-pressure solvers assume grid spacings of a staggered grid
!-- everywhere. For the actual grid, the grid spacing at the lowest level
!-- is only dz/2, but should be dz. Therefore, an additional array
!-- containing with appropriate grid information is created for these
!-- solvers.
    IF ( psolver(1:9) /= 'multigrid' )  THEN
       ALLOCATE( ddzu_pres(1:nzt+1) )
       ddzu_pres = ddzu
       ddzu_pres(1) = ddzu_pres(2)  ! change for lowest level
    ENDIF

!
!-- Compute the reciprocal values of the horizontal grid lengths.
    ddx = 1.0_wp / dx
    ddy = 1.0_wp / dy
    dx2 = dx * dx
    dy2 = dy * dy
    ddx2 = 1.0_wp / dx2
    ddy2 = 1.0_wp / dy2

!
!-- Allocate 3D array to set topography
    ALLOCATE( topo(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    topo = 0
!
!-- Initialize topography by generic topography or read topography from file.  
    CALL init_topo( topo )
!
!-- Set flags to mask topography on the grid. 
    CALL set_topo_flags( topo )

!
!-- Determine the maximum level of topography. It is used for
!-- steering the degradation of order of the applied advection scheme, 
!-- as well in the lpm.
    k_top = 0
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb, nzt + 1
             k_top = MAX( k_top, MERGE( k, 0, .NOT. BTEST( topo(k,j,i), 0 ) ) )
          ENDDO
       ENDDO
    ENDDO
#if defined( __parallel )
    CALL MPI_ALLREDUCE( k_top, nzb_max, 1, MPI_INTEGER,                        &
                        MPI_MAX, comm2d, ierr )
#else
    nzb_max = k_top
#endif
!
!-- Increment nzb_max by 1 in order to allow for proper diverengence correction.
!-- Further, in case topography extents up to the model top, limit to nzt.
    nzb_max = MIN( nzb_max+1, nzt ) 
!
!-- Determine minimum index of topography. Usually, this will be nzb. In case
!-- there is elevated topography, however, the lowest topography will be higher. 
!-- This index is e.g. used to calculate mean first-grid point atmosphere 
!-- temperature, surface pressure and density, etc. .
    topo_min_level   = 0
#if defined( __parallel )
    CALL MPI_ALLREDUCE( MINVAL( topo_top_ind(nys:nyn,nxl:nxr,0) ),             &
                        topo_min_level, 1, MPI_INTEGER, MPI_MIN, comm2d, ierr )
#else
    topo_min_level = MINVAL( topo_top_ind(nys:nyn,nxl:nxr,0) )
#endif

!
!-- Check topography for consistency with model domain. Therefore, use
!-- maximum and minium topography-top indices. Note, minimum topography top
!-- index is already calculated.  
    IF ( TRIM( topography ) /= 'flat' )  THEN
#if defined( __parallel )
       CALL MPI_ALLREDUCE( MAXVAL( topo_top_ind(nys:nyn,nxl:nxr,0) ),          &
                           nzb_local_max, 1, MPI_INTEGER, MPI_MAX, comm2d, ierr )               
#else
       nzb_local_max = MAXVAL( topo_top_ind(nys:nyn,nxl:nxr,0) )
#endif
       nzb_local_min = topo_min_level
!
!--    Consistency checks
       IF ( nzb_local_min < 0  .OR.  nzb_local_max  > nz + 1 )  THEN
          WRITE( message_string, * ) 'nzb_local values are outside the',       &
                                ' model domain',                               &
                                '&MINVAL( nzb_local ) = ', nzb_local_min,      &
                                '&MAXVAL( nzb_local ) = ', nzb_local_max
          CALL message( 'init_grid', 'PA0210', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF
!
!-- Define vertical gridpoint from (or to) which on the usual finite difference
!-- form (which does not use surface fluxes) is applied 
    IF ( constant_flux_layer  .OR.  use_surface_fluxes )  THEN
       nzb_diff = nzb + 2
    ELSE
       nzb_diff = nzb + 1
    ENDIF

    IF ( TRIM( topography ) /= 'flat' )  THEN
!
!--    Allocate and set the arrays containing the topography height (for output
!--    reasons only).
       IF ( nxr == nx  .AND.  nyn /= ny )  THEN
          ALLOCATE( zu_s_inner(nxl:nxr+1,nys:nyn),                             &
                    zw_w_inner(nxl:nxr+1,nys:nyn) )
       ELSEIF ( nxr /= nx  .AND.  nyn == ny )  THEN
          ALLOCATE( zu_s_inner(nxl:nxr,nys:nyn+1),                             &
                    zw_w_inner(nxl:nxr,nys:nyn+1) )
       ELSEIF ( nxr == nx  .AND.  nyn == ny )  THEN
          ALLOCATE( zu_s_inner(nxl:nxr+1,nys:nyn+1),                           &
                    zw_w_inner(nxl:nxr+1,nys:nyn+1) )
       ELSE
          ALLOCATE( zu_s_inner(nxl:nxr,nys:nyn),                               &
                    zw_w_inner(nxl:nxr,nys:nyn) )
       ENDIF

       zu_s_inner   = 0.0_wp
       zw_w_inner   = 0.0_wp
!
!--    Determine local topography height on scalar and w-grid. Note, setting
!--    lateral boundary values is not necessary, realized via wall_flags_static_0
!--    array. Further, please note that loop bounds are different from
!--    nxl to nxr and nys to nyn on south and right model boundary, hence,
!--    use intrinsic lbound and ubound functions to infer array bounds.
       DO  i = LBOUND(zu_s_inner, 1), UBOUND(zu_s_inner, 1)
          DO  j = LBOUND(zu_s_inner, 2), UBOUND(zu_s_inner, 2)
!
!--          Topography height on scalar grid. Therefore, determine index of
!--          upward-facing surface element on scalar grid.
             zu_s_inner(i,j) = zu(topo_top_ind(j,i,0))
!
!--          Topography height on w grid. Therefore, determine index of
!--          upward-facing surface element on w grid.
             zw_w_inner(i,j) = zw(topo_top_ind(j,i,3))
          ENDDO
       ENDDO
    ENDIF

#if defined( __parallel )
!
!-- Vertical nesting: communicate vertical grid level arrays between fine and
!-- coarse grid
    IF ( vnested )  CALL vnest_init_grid
#endif

 END SUBROUTINE init_grid


! Description:
! -----------------------------------------------------------------------------!
!> Calculation of the stretching factor through an iterative method. Ideas were 
!> taken from the paper "Regional stretched grid generation and its application
!> to the NCAR RegCM (1999)". Normally, no analytic solution exists because the
!> system of equations has two variables (r,l) but four requirements 
!> (l=integer, r=[0,88;1,2], Eq(6), Eq(5) starting from index j=1) which
!> results into an overdetermined system. 
!------------------------------------------------------------------------------!
 SUBROUTINE calculate_stretching_factor( number_end )
 
    USE control_parameters,                                                    &
        ONLY:  dz, dz_stretch_factor_array,                 &
               dz_stretch_level_end, dz_stretch_level_start, message_string
  
    USE kinds
    
    IMPLICIT NONE
    
    INTEGER(iwp) ::  iterations  !< number of iterations until stretch_factor_lower/upper_limit is reached  
    INTEGER(iwp) ::  l_rounded   !< after l_rounded grid levels dz(n) is strechted to dz(n+1) with stretch_factor_2 
    INTEGER(iwp) ::  n           !< loop variable for stretching
    
    INTEGER(iwp), INTENT(IN) ::  number_end !< number of user-specified end levels for stretching
        
    REAL(wp) ::  delta_l               !< absolute difference between l and l_rounded
    REAL(wp) ::  delta_stretch_factor  !< absolute difference between stretch_factor_1 and stretch_factor_2
    REAL(wp) ::  delta_total_new       !< sum of delta_l and delta_stretch_factor for the next iteration (should be as small as possible) 
    REAL(wp) ::  delta_total_old       !< sum of delta_l and delta_stretch_factor for the last iteration 
    REAL(wp) ::  distance              !< distance between dz_stretch_level_start and dz_stretch_level_end (stretching region)
    REAL(wp) ::  l                     !< value that fulfil Eq. (5) in the paper mentioned above together with stretch_factor_1 exactly
    REAL(wp) ::  numerator             !< numerator of the quotient
    REAL(wp) ::  stretch_factor_1      !< stretching factor that fulfil Eq. (5) togehter with l exactly
    REAL(wp) ::  stretch_factor_2      !< stretching factor that fulfil Eq. (6) togehter with l_rounded exactly
    
    REAL(wp) ::  dz_stretch_factor_array_2(9) = 1.08_wp  !< Array that contains all stretch_factor_2 that belongs to stretch_factor_1 
    
    REAL(wp), PARAMETER ::  stretch_factor_interval = 1.0E-06_wp  !< interval for sampling possible stretching factors
    REAL(wp), PARAMETER ::  stretch_factor_lower_limit = 0.88_wp  !< lowest possible stretching factor
    REAL(wp), PARAMETER ::  stretch_factor_upper_limit = 1.12_wp  !< highest possible stretching factor
 
 
    l = 0
    DO  n = 1, number_end
    
       iterations = 1
       stretch_factor_1 = 1.0_wp 
       stretch_factor_2 = 1.0_wp
       delta_total_old = 1.0_wp
       
!
!--    First branch for stretching from rough to fine
       IF ( dz(n) > dz(n+1) ) THEN
          DO WHILE ( stretch_factor_1 >= stretch_factor_lower_limit ) 
             
             stretch_factor_1 = 1.0_wp - iterations * stretch_factor_interval
             distance = ABS( dz_stretch_level_end(n) -                         &
                        dz_stretch_level_start(n) )   
             numerator = distance*stretch_factor_1/dz(n) +                     &
                         stretch_factor_1 - distance/dz(n)
             
             IF ( numerator > 0.0_wp ) THEN
                l = LOG( numerator ) / LOG( stretch_factor_1 ) - 1.0_wp
                l_rounded = NINT( l )
                delta_l = ABS( l_rounded - l ) / l
             ENDIF
             
             stretch_factor_2 = EXP( LOG( dz(n+1)/dz(n) ) / (l_rounded) )
             
             delta_stretch_factor = ABS( stretch_factor_1 -                    &
                                         stretch_factor_2 ) /                  &
                                    stretch_factor_2
             
             delta_total_new = delta_l + delta_stretch_factor

!
!--          stretch_factor_1 is taken to guarantee that the stretching
!--          procedure ends as close as possible to dz_stretch_level_end.
!--          stretch_factor_2 would guarantee that the stretched dz(n) is
!--          equal to dz(n+1) after l_rounded grid levels.
             IF (delta_total_new < delta_total_old) THEN
                dz_stretch_factor_array(n) = stretch_factor_1
                dz_stretch_factor_array_2(n) = stretch_factor_2
                delta_total_old = delta_total_new
             ENDIF
             
             iterations = iterations + 1
            
          ENDDO

!
!--    Second branch for stretching from fine to rough
       ELSEIF ( dz(n) < dz(n+1) ) THEN 
          DO WHILE ( stretch_factor_1 <= stretch_factor_upper_limit )
                     
             stretch_factor_1 = 1.0_wp + iterations * stretch_factor_interval
             distance = ABS( dz_stretch_level_end(n) -                         &
                        dz_stretch_level_start(n) ) 
             numerator = distance*stretch_factor_1/dz(n) +                     &
                         stretch_factor_1 - distance/dz(n)
             
             l = LOG( numerator ) / LOG( stretch_factor_1 ) - 1.0_wp
             l_rounded = NINT( l )
             delta_l = ABS( l_rounded - l ) / l
             
             stretch_factor_2 = EXP( LOG( dz(n+1)/dz(n) ) / (l_rounded) )

             delta_stretch_factor = ABS( stretch_factor_1 -                    &
                                        stretch_factor_2 ) /                   &
                                        stretch_factor_2
             
             delta_total_new = delta_l + delta_stretch_factor
             
!
!--          stretch_factor_1 is taken to guarantee that the stretching
!--          procedure ends as close as possible to dz_stretch_level_end.
!--          stretch_factor_2 would guarantee that the stretched dz(n) is
!--          equal to dz(n+1) after l_rounded grid levels.
             IF (delta_total_new < delta_total_old) THEN
                dz_stretch_factor_array(n) = stretch_factor_1
                dz_stretch_factor_array_2(n) = stretch_factor_2
                delta_total_old = delta_total_new
             ENDIF
             
             iterations = iterations + 1
          ENDDO
          
       ELSE
          message_string= 'Two adjacent values of dz must be different'
          CALL message( 'init_grid', 'PA0228', 1, 2, 0, 6, 0 )
          
       ENDIF

!
!--    Check if also the second stretching factor fits into the allowed
!--    interval. If not, print a warning for the user.
       IF ( dz_stretch_factor_array_2(n) < stretch_factor_lower_limit .OR.     & 
            dz_stretch_factor_array_2(n) > stretch_factor_upper_limit ) THEN
          WRITE( message_string, * ) 'stretch_factor_2 = ',                    &
                                     dz_stretch_factor_array_2(n), ' which is',&
                                     ' responsible for exactly reaching& dz =',&
                                      dz(n+1), 'after a specific amount of',   & 
                                     ' grid levels& exceeds the upper',        &
                                     ' limit =', stretch_factor_upper_limit,   &
                                     ' &or lower limit = ',                    &
                                     stretch_factor_lower_limit
          CALL message( 'init_grid', 'PA0499', 0, 1, 0, 6, 0 )
            
       ENDIF
    ENDDO
        
 END SUBROUTINE calculate_stretching_factor
 
 
! Description:
! -----------------------------------------------------------------------------!
!> Set temporary topography flags and reference buildings on top of underlying
!> orography. 
!------------------------------------------------------------------------------!
 SUBROUTINE process_topography( topo_3d )

    USE arrays_3d,                                                             &
        ONLY:  zu, zw

    USE control_parameters,                                                    &
        ONLY:  bc_lr_cyc, bc_ns_cyc, message_string, ocean_mode

    USE exchange_horiz_mod,                                                    &
        ONLY:  exchange_horiz_int, exchange_horiz_2d

    USE indices,                                                               &
        ONLY:  nbgp, nx, nxl, nxlg, nxr, nxrg, ny, nyn, nyng, nys, nysg, nzb,  &
               nzt

    USE netcdf_data_input_mod,                                                 &
        ONLY:  buildings_f, building_id_f, building_type_f, input_pids_static, &
               terrain_height_f

    USE kinds

    USE pegrid

    IMPLICIT NONE

    INTEGER(iwp) ::  i                !< running index along x-direction
    INTEGER(iwp) ::  j                !< running index along y-direction
    INTEGER(iwp) ::  k                !< running index along z-direction with respect to numeric grid
    INTEGER(iwp) ::  k2               !< running index along z-direction with respect to netcdf grid
    INTEGER(iwp) ::  nr               !< index variable indication maximum terrain height for respective building ID
    INTEGER(iwp) ::  num_build        !< counter for number of buildings
    INTEGER(iwp) ::  topo_top_index   !< orography top index, used to map 3D buildings onto terrain

#if defined( __parallel )
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  displace_dum        !< displacements of start addresses, used for MPI_ALLGATHERV
#endif
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  build_ids           !< building IDs on entire model domain
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  build_ids_final     !< building IDs on entire model domain, multiple occurences are sorted out 
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  build_ids_final_tmp !< temporary array used for resizing
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  build_ids_l         !< building IDs on local subdomain
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  build_ids_l_tmp     !< temporary array used to resize array of building IDs

    INTEGER(iwp), DIMENSION(0:numprocs-1) ::  num_buildings     !< number of buildings with different ID on entire model domain
    INTEGER(iwp), DIMENSION(0:numprocs-1) ::  num_buildings_l   !< number of buildings with different ID on local subdomain

    INTEGER(iwp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  topo_3d !< input array for 3D topography and dummy array for setting "outer"-flags

    REAL(wp)                            ::  ocean_offset        !< offset to consider inverse vertical coordinate at topography definition
    REAL(wp)                            ::  oro_min = 0.0_wp    !< minimum terrain height in entire model domain, used to reference terrain to zero
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  oro_max             !< maximum terrain height occupied by an building with certain id
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  oro_max_l           !< maximum terrain height occupied by an building with certain id, on local subdomain

!
!-- Reference lowest terrain height to zero. In case the minimum terrain height
!-- is non-zero, all grid points of the lower vertical grid levels might be 
!-- entirely below the surface, meaning a waste of computational resources.
!-- In order to avoid this, remove the lowest terrain height. Please note,
!-- in case of a nested run, the global minimum from all parent and childs 
!-- need to be remove to avoid steep edges at the child-domain boundaries.
    IF ( input_pids_static )  THEN
    
#if defined( __parallel ) 
       CALL MPI_ALLREDUCE( MINVAL( terrain_height_f%var ), oro_min, 1,         &
                           MPI_REAL, MPI_MIN, MPI_COMM_WORLD, ierr )
#else
       oro_min = MINVAL( terrain_height_f%var )
#endif
       terrain_height_f%var = terrain_height_f%var - oro_min
!                           
!--    Give an informative message that terrain height is referenced to zero    
       IF ( oro_min > 0.0_wp )  THEN
          WRITE( message_string, * ) 'Terrain height was internally shifted '//&
                          'downwards by ', oro_min, 'meter(s) to save ' //     &
                          'computational resources.'
          CALL message( 'init_grid', 'PA0505', 0, 0, 0, 6, 0 )
       ENDIF
    ENDIF    
    
!
!-- In the following, buildings and orography are further preprocessed 
!-- before they are mapped on the LES grid.
!-- Buildings are mapped on top of the orography by maintaining the roof 
!-- shape of the building. This can be achieved by referencing building on 
!-- top of the maximum terrain height within the area occupied by the 
!-- respective building. As buildings and terrain height are defined PE-wise,
!-- parallelization of this referencing is required (a building can be 
!-- distributed between different PEs).  
!-- In a first step, determine the number of buildings with different 
!-- building id on each PE. In a next step, all building ids are gathered
!-- into one array which is present to all PEs. For each building ID, 
!-- the maximum terrain height occupied by the respective building is 
!-- computed and distributed to each PE.  
!-- Finally, for each building id and its respective reference orography, 
!-- builidings are mapped on top.   
!-- 
!-- First, pre-set topography flags, bit 1 indicates orography, bit 2 
!-- buildings 
!-- classify the respective surfaces.
    topo_3d          = IBSET( topo_3d, 0 )
    topo_3d(nzb,:,:) = IBCLR( topo_3d(nzb,:,:), 0 )
!
!-- In order to map topography on PALM grid also in case of ocean simulations,
!-- pre-calculate an offset value.
    ocean_offset = MERGE( zw(0), 0.0_wp, ocean_mode )
!
!-- Reference buildings on top of orography. This is not necessary
!-- if topography is read from ASCII file as no distinction between buildings
!-- and terrain height can be made. Moreover, this is also not necessary if
!-- urban-surface and land-surface model are used at the same time. 
    IF ( input_pids_static )  THEN

       IF ( buildings_f%from_file )  THEN 
          num_buildings_l = 0
          num_buildings   = 0
!
!--       Allocate at least one element for building ids and give it an inital
!--       negative value that will be overwritten later. This, however, is 
!--       necessary in case there all IDs in the model domain are fill values. 
          ALLOCATE( build_ids_l(1) )
          build_ids_l = -1 
          DO  i = nxl, nxr
             DO  j = nys, nyn
                IF ( building_id_f%var(j,i) /= building_id_f%fill )  THEN
                   IF ( num_buildings_l(myid) > 0 )  THEN
                      IF ( ANY( building_id_f%var(j,i) ==  build_ids_l ) )   &
                      THEN
                         CYCLE
                      ELSE
                         num_buildings_l(myid) = num_buildings_l(myid) + 1
!
!--                   Resize array with different local building ids
                      ALLOCATE( build_ids_l_tmp(1:SIZE(build_ids_l)) )
                      build_ids_l_tmp = build_ids_l
                      DEALLOCATE( build_ids_l )
                      ALLOCATE( build_ids_l(1:num_buildings_l(myid)) )
                      build_ids_l(1:num_buildings_l(myid)-1) =                 &
                                  build_ids_l_tmp(1:num_buildings_l(myid)-1)
                      build_ids_l(num_buildings_l(myid)) = building_id_f%var(j,i)
                      DEALLOCATE( build_ids_l_tmp )
                   ENDIF
!
!--                First occuring building id on PE 
                   ELSE 
                      num_buildings_l(myid) = num_buildings_l(myid) + 1
                      build_ids_l(1) = building_id_f%var(j,i)
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
!
!--       Determine number of different building ids for the entire domain 
#if defined( __parallel ) 
          CALL MPI_ALLREDUCE( num_buildings_l, num_buildings, numprocs,              &
                              MPI_INTEGER, MPI_SUM, comm2d, ierr ) 
#else
          num_buildings = num_buildings_l
#endif
!
!--       Gather all buildings ids on each PEs. 
!--       First, allocate array encompassing all building ids in model domain.  
          ALLOCATE( build_ids(1:SUM(num_buildings)) )
#if defined( __parallel ) 
!
!--       Allocate array for displacements. 
!--       As each PE may has a different number of buildings, so that
!--       the block sizes send by each PE may not be equal. Hence, 
!--       information about the respective displacement is required, indicating 
!--       the respective adress where each MPI-task writes into the receive 
!--       buffer array  
          ALLOCATE( displace_dum(0:numprocs-1) )
          displace_dum(0) = 0
          DO i = 1, numprocs-1
             displace_dum(i) = displace_dum(i-1) + num_buildings(i-1)
          ENDDO

          CALL MPI_ALLGATHERV( build_ids_l(1:num_buildings_l(myid)),                 &
                               num_buildings(myid),                                  &
                               MPI_INTEGER,                                          &
                               build_ids,                                            &
                               num_buildings,                                        &
                               displace_dum,                                         & 
                               MPI_INTEGER,                                          &
                               comm2d, ierr )   

          DEALLOCATE( displace_dum )

#else
          build_ids = build_ids_l
#endif

!
!--       Note, in parallel mode building ids can occure mutliple times, as 
!--       each PE has send its own ids. Therefore, sort out building ids which 
!--       appear more than one time. 
          num_build = 0
          DO  nr = 1, SIZE(build_ids)

             IF ( ALLOCATED(build_ids_final) )  THEN
                IF ( ANY( build_ids(nr) == build_ids_final ) )  THEN
                   CYCLE
                ELSE
                   num_build = num_build + 1
!
!--                Resize
                   ALLOCATE( build_ids_final_tmp(1:num_build) )
                   build_ids_final_tmp(1:num_build-1) = build_ids_final(1:num_build-1)
                   DEALLOCATE( build_ids_final )
                   ALLOCATE( build_ids_final(1:num_build) )
                   build_ids_final(1:num_build-1) = build_ids_final_tmp(1:num_build-1)
                   build_ids_final(num_build) = build_ids(nr)
                   DEALLOCATE( build_ids_final_tmp )
                ENDIF             
             ELSE
                num_build = num_build + 1
                ALLOCATE( build_ids_final(1:num_build) )
                build_ids_final(num_build) = build_ids(nr)
             ENDIF
          ENDDO

!
!--       Determine maximumum terrain height occupied by the respective 
!--       building and temporalily store on oro_max 
          ALLOCATE( oro_max_l(1:SIZE(build_ids_final)) )
          ALLOCATE( oro_max(1:SIZE(build_ids_final))   )
          oro_max_l = 0.0_wp

          DO  nr = 1, SIZE(build_ids_final)
             oro_max_l(nr) = MAXVAL(                                           &
                              MERGE( terrain_height_f%var(nys:nyn,nxl:nxr),    &
                                     0.0_wp,                                   &
                                     building_id_f%var(nys:nyn,nxl:nxr) ==     &
                                     build_ids_final(nr) ) )
          ENDDO
   
#if defined( __parallel )    
          IF ( SIZE(build_ids_final) >= 1 ) THEN
             CALL MPI_ALLREDUCE( oro_max_l, oro_max, SIZE( oro_max ), MPI_REAL,&
                                 MPI_MAX, comm2d, ierr ) 
          ENDIF
#else
          oro_max = oro_max_l
#endif
!
!--       Finally, determine discrete grid height of maximum orography occupied
!--       by a building. Use all-or-nothing approach, i.e. if terrain
!--       exceeds the scalar level the grid box is fully terrain and the 
!--       maximum terrain is set to the zw level. 
!--       terrain or 
          oro_max_l = 0.0
          DO  nr = 1, SIZE(build_ids_final)
             DO  k = nzb, nzt
                IF ( zu(k) - ocean_offset <= oro_max(nr) )                     &
                   oro_max_l(nr) = zw(k) - ocean_offset
             ENDDO
             oro_max(nr) = oro_max_l(nr)
          ENDDO
       ENDIF
!
!--    Allocate array for storing terrain height under buildings
       IF ( buildings_f%from_file )  THEN
          ALLOCATE( buildings_f%oro_max(nysg:nyng,nxlg:nxrg) )
          buildings_f%oro_max = buildings_f%fill1
       END IF
!
!--    Map orography as well as buildings onto grid. 
       DO  i = nxl, nxr
          DO  j = nys, nyn
             topo_top_index = 0
!
!--          Obtain index in global building_id array
             IF ( buildings_f%from_file )  THEN
                IF ( building_id_f%var(j,i) /= building_id_f%fill )  THEN
!
!--                Determine index where maximum terrain height occupied by 
!--                the respective building height is stored.
                   nr = MINLOC( ABS( build_ids_final -                         &
                                     building_id_f%var(j,i) ), DIM = 1 )
!
!--                Save grid-indexed oro_max
                   buildings_f%oro_max(j,i) = oro_max(nr)
                ENDIF
             ENDIF
             DO  k = nzb, nzt
!
!--             In a first step, if grid point is below or equal the given 
!--             terrain height, grid point is flagged to be of type natural. 
!--             Please note, in case there is also a building which is lower
!--             than the vertical grid spacing, initialization of surface
!--             attributes will not be correct as given surface information
!--             will not be in accordance to the classified grid points. 
!--             Hence, in this case, also a building flag.
                IF ( zu(k) - ocean_offset <= terrain_height_f%var(j,i) )  THEN
                   topo_3d(k,j,i) = IBCLR( topo_3d(k,j,i), 0 )
                   topo_3d(k,j,i) = IBSET( topo_3d(k,j,i), 1 )
                   topo_top_index = k ! topo_top_index + 1
                ENDIF
!
!--             Set building grid points. Here, only consider 2D buildings. 
!--             3D buildings require separate treatment. 
                IF ( buildings_f%from_file  .AND.  buildings_f%lod == 1 )  THEN
!
!--                Fill-up the terrain to the level of maximum orography
!--                within the building-covered area.
                   IF ( building_id_f%var(j,i) /= building_id_f%fill )  THEN
!
!--                   Note, oro_max is always on zw level                   
                      IF ( zu(k) - ocean_offset < oro_max(nr) )  THEN
                         topo_3d(k,j,i) = IBCLR( topo_3d(k,j,i), 0 )
                         topo_3d(k,j,i) = IBSET( topo_3d(k,j,i), 1 )
                      ELSEIF ( zu(k) - ocean_offset <=                         &
                               oro_max(nr) + buildings_f%var_2d(j,i) )  THEN
                         topo_3d(k,j,i) = IBCLR( topo_3d(k,j,i), 0 )
                         topo_3d(k,j,i) = IBSET( topo_3d(k,j,i), 2 )
                      ENDIF
                   ENDIF
                ENDIF
             ENDDO
!
!--          Special treatment for non grid-resolved buildings. This case,
!--          the uppermost terrain grid point is flagged as building as well
!--          well, even though no building exists at all. However, the 
!--          surface element will be identified as urban-surface and the 
!--          input data provided by the drivers is consistent to the surface
!--          classification. Else, all non grid-resolved buildings would vanish
!--          and identified as terrain grid points, which, however, won't be
!--          consistent with the input data. 
             IF ( buildings_f%from_file  .AND.  buildings_f%lod == 1 )  THEN
                IF ( building_id_f%var(j,i) /= building_id_f%fill )  THEN
                   DO  k = nzb, nzt
                      IF( zw(k) - ocean_offset == oro_max(nr) )  THEN
                         IF ( buildings_f%var_2d(j,i) <= zu(k+1) - zw(k) )  THEN
                            topo_3d(k,j,i) = IBSET( topo_3d(k,j,i), 2 )
                         ENDIF
                      ENDIF
                   ENDDO
                ENDIF
             ENDIF
!
!--          Map 3D buildings onto terrain height.  
!--          In case of any slopes, map building on top of maximum terrain 
!--          height covered by the building. In other words, extend 
!--          building down to the respective local terrain-surface height. 
             IF ( buildings_f%from_file  .AND.  buildings_f%lod == 2 )  THEN
                IF ( building_id_f%var(j,i) /= building_id_f%fill )  THEN
!
!--                Extend building down to the terrain surface, i.e. fill-up
!--                surface irregularities below a building. Note, oro_max
!--                is already a discrete height according to the all-or-nothing
!--                approach, i.e. grid box is either topography or atmosphere, 
!--                terrain top is defined at upper bound of the grid box.
!--                Hence, check for zw in this case. 
!--                Note, do this only for buildings which are surface mounted, 
!--                i.e. building types 1-6. Below bridges, which are represented
!--                exclusively by building type 7, terrain shape should be
!--                maintained. 
                   IF ( building_type_f%from_file )  THEN
                      IF ( building_type_f%var(j,i) /= 7 )  THEN
                         DO k = topo_top_index + 1, nzt + 1      
                            IF ( zu(k) - ocean_offset <= oro_max(nr) )  THEN
                               topo_3d(k,j,i) = IBCLR( topo_3d(k,j,i), 0 )
                               topo_3d(k,j,i) = IBSET( topo_3d(k,j,i), 1 )
                            ENDIF
                         ENDDO        
!                     
!--                      After surface irregularities are smoothen, determine 
!--                      lower start index where building starts. 
                         DO  k = nzb, nzt
                            IF ( zu(k) - ocean_offset <= oro_max(nr) )         &
                               topo_top_index = k
                         ENDDO
                      ENDIF
                   ENDIF
!
!--                Finally, map building on top.
                   k2 = 0
                   DO k = topo_top_index, nzt + 1
                      IF ( k2 <= buildings_f%nz-1 )  THEN
                         IF ( buildings_f%var_3d(k2,j,i) == 1 )  THEN
                            topo_3d(k,j,i) = IBCLR( topo_3d(k,j,i), 0 )
                            topo_3d(k,j,i) = IBSET( topo_3d(k,j,i), 2 )
                         ENDIF
                      ENDIF
                      k2 = k2 + 1
                   ENDDO
                ENDIF
             ENDIF
          ENDDO
       ENDDO
!
!--    Horizontal exchange the oro_max array, which is required to for
!--    initialization of building-surface properties.
       IF ( ALLOCATED( buildings_f%oro_max ) )  THEN
          CALL exchange_horiz_2d( buildings_f%oro_max(:,:) )
       ENDIF
!
!--    Deallocate temporary arrays required for processing and reading data
       IF ( ALLOCATED( oro_max         ) )  DEALLOCATE( oro_max         )
       IF ( ALLOCATED( oro_max_l       ) )  DEALLOCATE( oro_max_l       )
       IF ( ALLOCATED( build_ids_final ) )  DEALLOCATE( build_ids_final )
!
!-- Topography input via ASCII format. 
    ELSE
       ocean_offset     = MERGE( zw(0), 0.0_wp, ocean_mode )
!
!--    Initialize topography bit 0 (indicates obstacle) everywhere to zero
!--    and clear all grid points at nzb, where alway a surface is defined.
!--    Further, set also bit 1 (indicates terrain) at nzb, which is further
!--    used for masked data output and further processing. Note, in the 
!--    ASCII case no distinction is made between buildings and terrain, 
!--    so that setting of bit 1 and 2 at the same time has no effect. 
       topo_3d          = IBSET( topo_3d, 0 )
       topo_3d(nzb,:,:) = IBCLR( topo_3d(nzb,:,:), 0 )
       topo_3d(nzb,:,:) = IBSET( topo_3d(nzb,:,:), 1 )
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb, nzt
!
!--             Flag topography for all grid points which are below
!--             the local topography height. 
!--             Note, each topography is flagged as building. 
                IF ( zu(k) - ocean_offset <= buildings_f%var_2d(j,i) )  THEN
                   topo_3d(k,j,i) = IBCLR( topo_3d(k,j,i), 0 )
                   topo_3d(k,j,i) = IBSET( topo_3d(k,j,i), 2 ) !indicates building
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    CALL exchange_horiz_int( topo_3d, nys, nyn, nxl, nxr, nzt, nbgp )

    IF ( .NOT. bc_ns_cyc )  THEN
       IF ( nys == 0  )  topo_3d(:,-1,:)   = topo_3d(:,0,:)
       IF ( nyn == ny )  topo_3d(:,ny+1,:) = topo_3d(:,ny,:)
    ENDIF

    IF ( .NOT. bc_lr_cyc )  THEN
       IF ( nxl == 0  )  topo_3d(:,:,-1)   = topo_3d(:,:,0)
       IF ( nxr == nx )  topo_3d(:,:,nx+1) = topo_3d(:,:,nx)          
    ENDIF

 END SUBROUTINE process_topography


! Description:
! -----------------------------------------------------------------------------!
!> Filter topography, i.e. fill holes resolved by only one grid point.  
!> Such holes are suspected to lead to velocity blow-ups as continuity 
!> equation on discrete grid cannot be fulfilled in such case.
!------------------------------------------------------------------------------!
 SUBROUTINE filter_topography( topo_3d )

    USE control_parameters,                                                    &
        ONLY:  bc_lr_cyc, bc_ns_cyc, message_string

    USE exchange_horiz_mod,                                                    &
        ONLY:  exchange_horiz_int, exchange_horiz_2d_byte, exchange_horiz_2d_int

    USE indices,                                                               &
        ONLY:  nbgp, nx, nxl, nxlg, nxr, nxrg, ny, nyn, nyng, nys, nysg, nzb, nzt

    USE netcdf_data_input_mod,                                                 &
        ONLY:  building_id_f, building_type_f 

    USE  pegrid

    IMPLICIT NONE

    LOGICAL      ::  filled = .FALSE. !< flag indicating if holes were filled

    INTEGER(iwp) ::  i          !< running index along x-direction
    INTEGER(iwp) ::  j          !< running index along y-direction
    INTEGER(iwp) ::  k          !< running index along z-direction
    INTEGER(iwp) ::  num_hole   !< number of holes (in topography) resolved by only one grid point 
    INTEGER(iwp) ::  num_hole_l !< number of holes (in topography) resolved by only one grid point on local PE     
    INTEGER(iwp) ::  num_wall   !< number of surrounding vertical walls for a single grid point

    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE            ::  topo_tmp          !< temporary 3D-topography used to fill holes
    INTEGER(iwp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  topo_3d           !< 3D-topography array merging buildings and orography
!
!-- Before checking for holes, set lateral boundary conditions for 
!-- topography. After hole-filling, boundary conditions must be set again.
!-- Several iterations are performed, in order to fill holes which might 
!-- emerge by the filling-algorithm itself.
    ALLOCATE( topo_tmp(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    topo_tmp = 0

    num_hole = 99999
    DO WHILE ( num_hole > 0 )       

       num_hole = 0    
       CALL exchange_horiz_int( topo_3d, nys, nyn, nxl, nxr, nzt, nbgp )
!
!--    Exchange also building ID and type. Note, building_type is an one-byte 
!--    variable.
       IF ( building_id_f%from_file )                                          &
          CALL exchange_horiz_2d_int( building_id_f%var, nys, nyn, nxl, nxr, nbgp )
       IF ( building_type_f%from_file )                                        &
          CALL exchange_horiz_2d_byte( building_type_f%var, nys, nyn, nxl, nxr, nbgp )

       topo_tmp = topo_3d
!
!--    In case of non-cyclic lateral boundaries, assume lateral boundary to be 
!--    a solid wall. Thus, intermediate spaces of one grid point between 
!--    boundary and some topographic structure will be filled.           
       IF ( .NOT. bc_ns_cyc )  THEN
          IF ( nys == 0  )  topo_tmp(:,-1,:)   = IBCLR( topo_tmp(:,0,:),  0 )
          IF ( nyn == ny )  topo_tmp(:,ny+1,:) = IBCLR( topo_tmp(:,ny,:), 0 )
       ENDIF

       IF ( .NOT. bc_lr_cyc )  THEN
          IF ( nxl == 0  )  topo_tmp(:,:,-1)   = IBCLR( topo_tmp(:,:,0),  0 )
          IF ( nxr == nx )  topo_tmp(:,:,nx+1) = IBCLR( topo_tmp(:,:,nx), 0 )          
       ENDIF

       num_hole_l = 0
       DO i = nxl, nxr
          DO j = nys, nyn
             DO  k = nzb+1, nzt
                IF ( BTEST( topo_tmp(k,j,i), 0 ) )  THEN
                   num_wall = 0
                   IF ( .NOT. BTEST( topo_tmp(k,j-1,i), 0 ) )                  &
                      num_wall = num_wall + 1
                   IF ( .NOT. BTEST( topo_tmp(k,j+1,i), 0 ) )                  &
                      num_wall = num_wall + 1
                   IF ( .NOT. BTEST( topo_tmp(k,j,i-1), 0 ) )                  &
                      num_wall = num_wall + 1
                   IF ( .NOT. BTEST( topo_tmp(k,j,i+1), 0 ) )                  &
                      num_wall = num_wall + 1
                   IF ( .NOT. BTEST( topo_tmp(k-1,j,i), 0 ) )                  &
                      num_wall = num_wall + 1   
                   IF ( .NOT. BTEST( topo_tmp(k+1,j,i), 0 ) )                  &
                      num_wall = num_wall + 1

                   IF ( num_wall >= 4 )  THEN
                      num_hole_l     = num_hole_l + 1
!
!--                   Clear flag 0 and set special flag ( bit 4) to indicate 
!--                   that new topography point is a result of filtering process.
                      topo_3d(k,j,i) = IBCLR( topo_3d(k,j,i), 0 )
                      topo_3d(k,j,i) = IBSET( topo_3d(k,j,i), 4 )
!
!--                   If filled grid point is occupied by a building, classify
!--                   it as building grid point.
                      IF ( building_type_f%from_file )  THEN
                         IF ( building_type_f%var(j,i)   /=                    &  
                              building_type_f%fill            .OR.             &       
                              building_type_f%var(j+1,i) /=                    &  
                              building_type_f%fill            .OR.             &                
                              building_type_f%var(j-1,i) /=                    &                
                              building_type_f%fill            .OR.             &                
                              building_type_f%var(j,i+1) /=                    &                
                              building_type_f%fill            .OR.             &                
                              building_type_f%var(j,i-1) /=                    &                
                              building_type_f%fill )  THEN
!
!--                         Set flag indicating building surfaces
                            topo_3d(k,j,i) = IBSET( topo_3d(k,j,i), 2 )
!
!--                         Set building_type and ID at this position if not 
!--                         already set. This is required for proper 
!--                         initialization of urban-surface energy balance 
!--                         solver.
                            IF ( building_type_f%var(j,i) ==                   &
                                 building_type_f%fill )  THEN

                               IF ( building_type_f%var(j+1,i) /=              &
                                    building_type_f%fill )  THEN
                                  building_type_f%var(j,i) =                   &
                                                    building_type_f%var(j+1,i)
                                  building_id_f%var(j,i) =                     &
                                                    building_id_f%var(j+1,i)
                               ELSEIF ( building_type_f%var(j-1,i) /=          &
                                        building_type_f%fill )  THEN
                                  building_type_f%var(j,i) =                   &
                                                    building_type_f%var(j-1,i)
                                  building_id_f%var(j,i) =                     &
                                                    building_id_f%var(j-1,i)
                               ELSEIF ( building_type_f%var(j,i+1) /=          &
                                        building_type_f%fill )  THEN
                                  building_type_f%var(j,i) =                   &
                                                    building_type_f%var(j,i+1)
                                  building_id_f%var(j,i) =                     &
                                                    building_id_f%var(j,i+1)
                               ELSEIF ( building_type_f%var(j,i-1) /=          &
                                        building_type_f%fill )  THEN
                                  building_type_f%var(j,i) =                   &
                                                    building_type_f%var(j,i-1)
                                  building_id_f%var(j,i) =                     &
                                                    building_id_f%var(j,i-1)
                               ENDIF
                            ENDIF
                         ENDIF
                      ENDIF
!
!--                   If filled grid point is already classified as building
!--                   everything is fine, else classify this grid point as
!--                   natural type grid point. This case, values for the 
!--                   surface type are already set.
                      IF ( .NOT. BTEST( topo_3d(k,j,i), 2 ) )  THEN
                         topo_3d(k,j,i) = IBSET( topo_3d(k,j,i), 1 )
                      ENDIF
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO
!
!--    Count the total number of holes, required for informative message.
#if defined( __parallel )
       CALL MPI_ALLREDUCE( num_hole_l, num_hole, 1, MPI_INTEGER, MPI_SUM,      &
                           comm2d, ierr )
#else
       num_hole = num_hole_l
#endif    
       IF ( num_hole > 0  .AND.  .NOT. filled )  filled = .TRUE.

    ENDDO
!
!-- Create an informative message if any holes were filled.
    IF ( filled )  THEN
       WRITE( message_string, * ) 'Topography was filtered, i.e. holes ' //    &
                                  'resolved by only one grid point '     //    &
                                  'were filled during initialization.'
       CALL message( 'init_grid', 'PA0430', 0, 0, 0, 6, 0 )
    ENDIF

    DEALLOCATE( topo_tmp )
!
!-- Finally, exchange topo_3d array again and if necessary set Neumann boundary
!-- condition in case of non-cyclic lateral boundaries. 
    CALL exchange_horiz_int( topo_3d, nys, nyn, nxl, nxr, nzt, nbgp )

    IF ( .NOT. bc_ns_cyc )  THEN
       IF ( nys == 0  )  topo_3d(:,-1,:)   = topo_3d(:,0,:)
       IF ( nyn == ny )  topo_3d(:,ny+1,:) = topo_3d(:,ny,:)
    ENDIF

    IF ( .NOT. bc_lr_cyc )  THEN
       IF ( nxl == 0  )  topo_3d(:,:,-1)   = topo_3d(:,:,0)
       IF ( nxr == nx )  topo_3d(:,:,nx+1) = topo_3d(:,:,nx)          
    ENDIF
!
!-- Exchange building ID and type. Note, building_type is an one-byte variable.
    IF ( building_id_f%from_file )                                             &
       CALL exchange_horiz_2d_int( building_id_f%var, nys, nyn, nxl, nxr, nbgp )
    IF ( building_type_f%from_file )                                           &
       CALL exchange_horiz_2d_byte( building_type_f%var, nys, nyn, nxl, nxr, nbgp )

 END SUBROUTINE filter_topography


! Description:
! -----------------------------------------------------------------------------!
!> Reads topography information from file or sets generic topography. Moreover, 
!> all topography-relevant topography arrays are initialized, and grid flags
!> are set.  
!------------------------------------------------------------------------------!
 SUBROUTINE init_topo( topo )

    USE arrays_3d,                                                             &
        ONLY:  zw
        
    USE control_parameters,                                                    &
        ONLY:  bc_lr_cyc, bc_ns_cyc, building_height, building_length_x,       &
               building_length_y, building_wall_left, building_wall_south,     &
               canyon_height, canyon_wall_left, canyon_wall_south,             &
               canyon_width_x, canyon_width_y, dp_level_ind_b, dz,             &
               message_string, topography, topography_grid_convention,         &
               tunnel_height, tunnel_length, tunnel_width_x, tunnel_width_y,   &
               tunnel_wall_depth
         
    USE exchange_horiz_mod,                                                    &
        ONLY:  exchange_horiz_int

    USE grid_variables,                                                        &
        ONLY:  dx, dy
        
    USE indices,                                                               &
        ONLY:  nbgp, nx, nxl, nxlg, nxr, nxrg, ny, nyn, nyng, nys, nysg, nz,   &
               nzb, nzt
    
    USE kinds
    
    USE netcdf_data_input_mod,                                                 &
        ONLY:  buildings_f, terrain_height_f 

    USE pegrid

    IMPLICIT NONE

    INTEGER(iwp) ::  bh                !< temporary vertical index of building height
    INTEGER(iwp) ::  ngp_bx            !< grid point number of building size along x
    INTEGER(iwp) ::  ngp_by            !< grid point number of building size along y
    INTEGER(iwp) ::  index_left_bwall  !< index for left building wall
    INTEGER(iwp) ::  index_right_bwall !< index for right building wall
    INTEGER(iwp) ::  index_north_bwall !< index for north building wall
    INTEGER(iwp) ::  index_south_bwall !< index for south building wall
    INTEGER(iwp) ::  ch                !< temporary vertical index for canyon height
    INTEGER(iwp) ::  ngp_cx            !< grid point number of canyon size along x
    INTEGER(iwp) ::  ngp_cy            !< grid point number of canyon size along y
    INTEGER(iwp) ::  index_left_cwall  !< index for left canyon wall
    INTEGER(iwp) ::  index_right_cwall !< index for right canyon wall
    INTEGER(iwp) ::  index_north_cwall !< index for north canyon wall
    INTEGER(iwp) ::  index_south_cwall !< index for south canyon wall
    INTEGER(iwp) ::  i                 !< index variable along x
    INTEGER(iwp) ::  j                 !< index variable along y
    INTEGER(iwp) ::  k                 !< index variable along z
    INTEGER(iwp) ::  hv_in             !< heavyside function to model inner tunnel surface 
    INTEGER(iwp) ::  hv_out            !< heavyside function to model outer tunnel surface 
    INTEGER(iwp) ::  txe_out           !< end position of outer tunnel wall in x
    INTEGER(iwp) ::  txs_out           !< start position of outer tunnel wall in x
    INTEGER(iwp) ::  tye_out           !< end position of outer tunnel wall in y
    INTEGER(iwp) ::  tys_out           !< start position of outer tunnel wall in y
    INTEGER(iwp) ::  txe_in            !< end position of inner tunnel wall in x
    INTEGER(iwp) ::  txs_in            !< start position of inner tunnel wall in x
    INTEGER(iwp) ::  tye_in            !< end position of inner tunnel wall in y
    INTEGER(iwp) ::  tys_in            !< start position of inner tunnel wall in y
    INTEGER(iwp) ::  td                !< tunnel wall depth
    INTEGER(iwp) ::  th                !< height of outer tunnel wall

    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  nzb_local         !< index for topography top at cell-center
    INTEGER(iwp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  topo !< input array for 3D topography and dummy array for setting "outer"-flags
!
!-- Check for correct setting of the namelist parameter topography. If 
!-- topography information is read from file but topography = 'flat', 
!-- initialization does not work properly.
    IF ( ( buildings_f%from_file  .OR.  terrain_height_f%from_file )  .AND.    &
           TRIM( topography ) /= 'read_from_file' )  THEN
       message_string =  'If topography information is provided (via ' //      &
                         'Netcdf or ASCII input) topography = '        //      &
                         '"read_from_file" is required.'
       CALL message( 'init_grid', 'PA0437', 1, 2, 0, 6, 0 )      
    ENDIF
!
!-- Set outer and inner index arrays for non-flat topography.
!-- Here consistency checks concerning domain size and periodicity are
!-- necessary.
!-- Within this SELECT CASE structure only nzb_local is initialized
!-- individually depending on the chosen topography type, all other index 
!-- arrays are initialized further below.
    SELECT CASE ( TRIM( topography ) )

       CASE ( 'flat' )
!   
!--       Initialilize 3D topography array, used later for initializing flags
          topo(nzb+1:nzt+1,:,:) = IBSET( topo(nzb+1:nzt+1,:,:), 0 )
          
       CASE ( 'closed_channel' )
!   
!--       Initialilize 3D topography array, used later for initializing flags
          topo(nzb+1:nzt,:,:) = IBSET( topo(nzb+1:nzt,:,:), 0 ) 

       CASE ( 'single_building' )
!
!--       Single rectangular building, by default centered in the middle of the
!--       total domain
          ngp_bx = NINT( building_length_x / dx )
          ngp_by = NINT( building_length_y / dy )
          bh  = MINLOC( ABS( zw - building_height ), 1 ) - 1
          IF ( ABS( zw(bh)   - building_height ) == &
               ABS( zw(bh+1) - building_height )    )  bh = bh + 1
          IF ( building_wall_left == 9999999.9_wp )  THEN
             building_wall_left = ( nx + 1 - ngp_bx ) / 2 * dx
          ENDIF
          index_left_bwall = NINT( building_wall_left / dx )
          index_right_bwall = index_left_bwall + ngp_bx

          IF ( building_wall_south == 9999999.9_wp )  THEN
              building_wall_south = ( ny + 1 - ngp_by ) / 2 * dy
          ENDIF
          index_south_bwall = NINT( building_wall_south / dy )
          index_north_bwall = index_south_bwall + ngp_by

!
!--       Building size has to meet some requirements
          IF ( ( index_left_bwall < 1 ) .OR. ( index_right_bwall > nx-1 ) .OR. &
               ( index_right_bwall < index_left_bwall+3 ) .OR.                 &
               ( index_south_bwall < 1 ) .OR. ( index_north_bwall > ny-1 ) .OR.&
               ( index_north_bwall < index_south_bwall+3 ) )  THEN
             WRITE( message_string, * ) 'inconsistent building parameters:',   &
                                      '&index_left_bwall=', index_left_bwall,  &
                                      'index_right_bwall=', index_right_bwall, &
                                      'index_south_bwall=', index_south_bwall, &
                                      'index_north_bwall=', index_north_bwall, &
                                      'nx=', nx, 'ny=', ny
             CALL message( 'init_grid', 'PA0203', 1, 2, 0, 6, 0 )
          ENDIF

          ALLOCATE( nzb_local(nysg:nyng,nxlg:nxrg) )
          nzb_local = 0
!
!--       Define the building. 
          IF ( index_left_bwall <= nxr  .AND.  index_right_bwall >= nxl  .AND. &
               index_south_bwall <= nyn  .AND.  index_north_bwall >= nys )     & 
             nzb_local(MAX(nys,index_south_bwall):MIN(nyn,index_north_bwall),  &
                       MAX(nxl,index_left_bwall):MIN(nxr,index_right_bwall)) = bh
!
!--       Set bit array on basis of nzb_local
          DO  i = nxl, nxr
             DO  j = nys, nyn
                topo(nzb_local(j,i)+1:nzt+1,j,i) =                             &
                                 IBSET( topo(nzb_local(j,i)+1:nzt+1,j,i), 0 ) 
             ENDDO
          ENDDO
        
          DEALLOCATE( nzb_local )

          CALL exchange_horiz_int( topo, nys, nyn, nxl, nxr, nzt, nbgp )
!
!--       Set boundary conditions also for flags. Can be interpreted as Neumann
!--       boundary conditions for topography. 
          IF ( .NOT. bc_ns_cyc )  THEN
             IF ( nys == 0  )  THEN
                DO  i = 1, nbgp     
                   topo(:,nys-i,:)   = topo(:,nys,:)
                ENDDO
             ENDIF
             IF ( nyn == ny )  THEN
                DO  i = 1, nbgp  
                   topo(:,nyn+i,:) = topo(:,nyn,:)
                ENDDO
             ENDIF
          ENDIF
          IF ( .NOT. bc_lr_cyc )  THEN
             IF ( nxl == 0  )  THEN
                DO  i = 1, nbgp   
                   topo(:,:,nxl-i)   = topo(:,:,nxl)
                ENDDO
             ENDIF
             IF ( nxr == nx )  THEN 
                DO  i = 1, nbgp   
                   topo(:,:,nxr+i) = topo(:,:,nxr)     
                ENDDO
             ENDIF      
          ENDIF

       CASE ( 'single_street_canyon' )
!
!--       Single quasi-2D street canyon of infinite length in x or y direction.
!--       The canyon is centered in the other direction by default.
          IF ( canyon_width_x /= 9999999.9_wp )  THEN
!
!--          Street canyon in y direction
             ngp_cx = NINT( canyon_width_x / dx )
             IF ( canyon_wall_left == 9999999.9_wp )  THEN
                canyon_wall_left = ( nx + 1 - ngp_cx ) / 2 * dx
             ENDIF
             index_left_cwall= NINT( canyon_wall_left / dx )
             index_right_cwall= index_left_cwall+ ngp_cx
          ELSEIF ( canyon_width_y /= 9999999.9_wp )  THEN
!
!--          Street canyon in x direction
             ngp_cy = NINT( canyon_width_y / dy )
             IF ( canyon_wall_south == 9999999.9_wp )  THEN
                canyon_wall_south = ( ny + 1 - ngp_cy ) / 2 * dy
             ENDIF
             index_south_cwall = NINT( canyon_wall_south / dy )
             index_north_cwall = index_south_cwall + ngp_cy
      
          ELSE
             
             message_string = 'no street canyon width given'
             CALL message( 'init_grid', 'PA0204', 1, 2, 0, 6, 0 )
  
          ENDIF

          ch  = MINLOC( ABS( zw - canyon_height ), 1 ) - 1
          IF ( ABS( zw(ch)   - canyon_height ) == &
               ABS( zw(ch+1) - canyon_height )    )  ch = ch + 1
          dp_level_ind_b = ch
!
!--       Street canyon size has to meet some requirements
          IF ( canyon_width_x /= 9999999.9_wp )  THEN
             IF ( ( index_left_cwall< 1 ) .OR. ( index_right_cwall> nx-1 ) .OR.&
                  ( ngp_cx < 3 ) .OR. ( ch < 3 ) )  THEN
                WRITE( message_string, * ) 'inconsistent canyon parameters:',  &
                                           '&index_left_cwall=', index_left_cwall, &
                                           ' index_right_cwall=', index_right_cwall, &
                                           ' ngp_cx=', ngp_cx,                 &
                                           ' ch=', ch, ' nx=', nx, ' ny=', ny
                CALL message( 'init_grid', 'PA0205', 1, 2, 0, 6, 0 ) 
             ENDIF
          ELSEIF ( canyon_width_y /= 9999999.9_wp )  THEN
             IF ( ( index_south_cwall < 1 ) .OR.                               &
                  ( index_north_cwall > ny-1 ) .OR. ( ngp_cy < 3 ) .OR.        &
                  ( ch < 3 ) )  THEN
                WRITE( message_string, * ) 'inconsistent canyon parameters:',  &
                                           '&index_south_cwall=', index_south_cwall, & 
                                           ' index_north_cwall=', index_north_cwall, &
                                           ' ngp_cy=', ngp_cy,                 &
                                           ' ch=', ch, ' nx=', nx, ' ny=', ny
                CALL message( 'init_grid', 'PA0206', 1, 2, 0, 6, 0 ) 
             ENDIF
          ENDIF
          IF ( canyon_width_x /= 9999999.9_wp .AND.                            &                 
               canyon_width_y /= 9999999.9_wp )  THEN
             message_string = 'inconsistent canyon parameters:' //             &   
                              '&street canyon can only be oriented' //         &
                              ' either in x- or in y-direction'
             CALL message( 'init_grid', 'PA0207', 1, 2, 0, 6, 0 )
          ENDIF

          ALLOCATE( nzb_local(nysg:nyng,nxlg:nxrg) )
          nzb_local = ch
          IF ( canyon_width_x /= 9999999.9_wp )  THEN
             IF ( index_left_cwall<= nxr  .AND.  index_right_cwall>= nxl )     &
                nzb_local(:,MAX(nxl,index_left_cwall+1):MIN(nxr,index_right_cwall-1)) = 0
          ELSEIF ( canyon_width_y /= 9999999.9_wp )  THEN
             IF ( index_south_cwall <= nyn  .AND.  index_north_cwall >= nys )  &          
                nzb_local(MAX(nys,index_south_cwall+1):MIN(nyn,index_north_cwall-1),:) = 0
          ENDIF
!
!--       Set bit array on basis of nzb_local
          DO  i = nxl, nxr
             DO  j = nys, nyn
                topo(nzb_local(j,i)+1:nzt+1,j,i) =                             &
                                 IBSET( topo(nzb_local(j,i)+1:nzt+1,j,i), 0 ) 
             ENDDO
          ENDDO
          DEALLOCATE( nzb_local )

          CALL exchange_horiz_int( topo, nys, nyn, nxl, nxr, nzt, nbgp )
!
!--       Set boundary conditions also for flags. Can be interpreted as Neumann
!--       boundary conditions for topography. 
          IF ( .NOT. bc_ns_cyc )  THEN
             IF ( nys == 0  )  THEN
                DO  i = 1, nbgp     
                   topo(:,nys-i,:)   = topo(:,nys,:)
                ENDDO
             ENDIF
             IF ( nyn == ny )  THEN
                DO  i = 1, nbgp  
                   topo(:,nyn+i,:) = topo(:,nyn,:)
                ENDDO
             ENDIF
          ENDIF
          IF ( .NOT. bc_lr_cyc )  THEN
             IF ( nxl == 0  )  THEN
                DO  i = 1, nbgp   
                   topo(:,:,nxl-i)   = topo(:,:,nxl)
                ENDDO
             ENDIF
             IF ( nxr == nx )  THEN 
                DO  i = 1, nbgp   
                   topo(:,:,nxr+i) = topo(:,:,nxr)     
                ENDDO
             ENDIF      
          ENDIF

       CASE ( 'tunnel' )

!
!--       Tunnel height
          IF ( tunnel_height == 9999999.9_wp )  THEN
             th = zw( INT( 0.2 * nz) )
          ELSE
             th = tunnel_height
          ENDIF
!
!--       Tunnel-wall depth
          IF ( tunnel_wall_depth == 9999999.9_wp )  THEN  
             td = MAX ( dx, dy, dz(1) )
          ELSE
             td = tunnel_wall_depth
          ENDIF
!
!--       Check for tunnel width
          IF ( tunnel_width_x == 9999999.9_wp  .AND.                           &
               tunnel_width_y == 9999999.9_wp  )  THEN
             message_string = 'No tunnel width is given. '
             CALL message( 'init_grid', 'PA0280', 1, 2, 0, 6, 0 )
          ENDIF
          IF ( tunnel_width_x /= 9999999.9_wp  .AND.                           &
               tunnel_width_y /= 9999999.9_wp  )  THEN
             message_string = 'Inconsistent tunnel parameters:' //             &   
                              'tunnel can only be oriented' //                 &
                              'either in x- or in y-direction.'
             CALL message( 'init_grid', 'PA0281', 1, 2, 0, 6, 0 )
          ENDIF
!
!--       Check for too small tunnel width in x- and y-direction
          IF ( tunnel_width_x /= 9999999.9_wp  .AND.                           &   
               tunnel_width_x - 2.0_wp * td <= 2.0_wp * dx )  THEN
             message_string = 'tunnel_width_x too small'
             CALL message( 'init_grid', 'PA0175', 1, 2, 0, 6, 0 )
          ENDIF
          IF ( tunnel_width_y /= 9999999.9_wp  .AND.                           &
               tunnel_width_y - 2.0_wp * td <= 2.0_wp * dy )  THEN
             message_string = 'tunnel_width_y too small'
             CALL message( 'init_grid', 'PA0455', 1, 2, 0, 6, 0 )
          ENDIF
!
!--       Check for too large tunnel width. 
!--       Tunnel axis along y.
          IF ( tunnel_width_x /= 9999999.9_wp )  THEN
             IF ( tunnel_width_x > ( nx + 1 ) * dx )  THEN
                message_string = 'tunnel_width_x too large'
                CALL message( 'init_grid', 'PA0282', 1, 2, 0, 6, 0 )
             ENDIF

             txs_out = INT( ( nx + 1 ) * 0.5_wp * dx - tunnel_width_x * 0.5_wp )
             txe_out = INT( ( nx + 1 ) * 0.5_wp * dx + tunnel_width_x * 0.5_wp )
             txs_in  = INT( ( nx + 1 ) * 0.5_wp * dx -                         &
                                      ( tunnel_width_x * 0.5_wp - td ) )
             txe_in  = INT( ( nx + 1 ) * 0.5_wp * dx +                         &
                                   ( tunnel_width_x * 0.5_wp - td ) )

             tys_out = INT( ( ny + 1 ) * 0.5_wp * dy - tunnel_length * 0.5_wp )
             tye_out = INT( ( ny + 1 ) * 0.5_wp * dy + tunnel_length * 0.5_wp )
             tys_in  = tys_out
             tye_in  = tye_out
          ENDIF
!
!--       Tunnel axis along x.
          IF ( tunnel_width_y /= 9999999.9_wp )  THEN
             IF ( tunnel_width_y > ( ny + 1 ) * dy )  THEN
                message_string = 'tunnel_width_y too large'
                CALL message( 'init_grid', 'PA0456', 1, 2, 0, 6, 0 )
             ENDIF

             txs_out = INT( ( nx + 1 ) * 0.5_wp * dx - tunnel_length * 0.5_wp )
             txe_out = INT( ( nx + 1 ) * 0.5_wp * dx + tunnel_length * 0.5_wp )
             txs_in  = txs_out
             txe_in  = txe_out

             tys_out = INT( ( ny + 1 ) * 0.5_wp * dy - tunnel_width_y * 0.5_wp )
             tye_out = INT( ( ny + 1 ) * 0.5_wp * dy + tunnel_width_y * 0.5_wp )
             tys_in  = INT( ( ny + 1 ) * 0.5_wp * dy -                         &
                                        ( tunnel_width_y * 0.5_wp - td ) )
             tye_in  = INT( ( ny + 1 ) * 0.5_wp * dy +                         &
                                     ( tunnel_width_y * 0.5_wp - td ) )
          ENDIF

          topo = 0
          DO  i = nxl, nxr
             DO  j = nys, nyn
!
!--             Use heaviside function to model outer tunnel surface
                hv_out = th * 0.5_wp *                                         &
                              ( ( SIGN( 1.0_wp, i * dx - txs_out ) + 1.0_wp )  &
                              - ( SIGN( 1.0_wp, i * dx - txe_out ) + 1.0_wp ) )

                hv_out = hv_out * 0.5_wp *                                     &
                            ( ( SIGN( 1.0_wp, j * dy - tys_out ) + 1.0_wp )    &
                            - ( SIGN( 1.0_wp, j * dy - tye_out ) + 1.0_wp ) )
!    
!--             Use heaviside function to model inner tunnel surface
                hv_in  = ( th - td ) * 0.5_wp *                                &
                                ( ( SIGN( 1.0_wp, i * dx - txs_in ) + 1.0_wp ) &
                                - ( SIGN( 1.0_wp, i * dx - txe_in ) + 1.0_wp ) )

                hv_in = hv_in * 0.5_wp *                                       &
                                ( ( SIGN( 1.0_wp, j * dy - tys_in ) + 1.0_wp ) &
                                - ( SIGN( 1.0_wp, j * dy - tye_in ) + 1.0_wp ) )
!
!--             Set flags at x-y-positions without any tunnel surface
                IF ( hv_out - hv_in == 0.0_wp )  THEN
                   topo(nzb+1:nzt+1,j,i) = IBSET( topo(nzb+1:nzt+1,j,i), 0 )
!
!--             Set flags at x-y-positions with tunnel surfaces
                ELSE
                   DO  k = nzb + 1, nzt + 1
!
!--                   Inner tunnel
                      IF ( hv_out - hv_in == th )  THEN
                         IF ( zw(k) <= hv_out )  THEN
                            topo(k,j,i) = IBCLR( topo(k,j,i), 0 )
                         ELSE
                            topo(k,j,i) = IBSET( topo(k,j,i), 0 )
                         ENDIF
                      ENDIF
!
!--                   Lateral tunnel walls
                      IF ( hv_out - hv_in == td )  THEN
                         IF ( zw(k) <= hv_in )  THEN 
                            topo(k,j,i) = IBSET( topo(k,j,i), 0 )
                         ELSEIF ( zw(k) > hv_in  .AND.  zw(k) <= hv_out )  THEN 
                            topo(k,j,i) = IBCLR( topo(k,j,i), 0 )
                         ELSEIF ( zw(k) > hv_out )  THEN 
                            topo(k,j,i) = IBSET( topo(k,j,i), 0 )
                         ENDIF
                      ENDIF
                   ENDDO
                ENDIF
             ENDDO
          ENDDO

          CALL exchange_horiz_int( topo, nys, nyn, nxl, nxr, nzt, nbgp )
!
!--       Set boundary conditions also for flags. Can be interpreted as Neumann
!--       boundary conditions for topography. 
          IF ( .NOT. bc_ns_cyc )  THEN
             IF ( nys == 0  )  THEN
                DO  i = 1, nbgp     
                   topo(:,nys-i,:)   = topo(:,nys,:)
                ENDDO
             ENDIF
             IF ( nyn == ny )  THEN
                DO  i = 1, nbgp  
                   topo(:,nyn+i,:) = topo(:,nyn,:)
                ENDDO
             ENDIF
          ENDIF
          IF ( .NOT. bc_lr_cyc )  THEN
             IF ( nxl == 0  )  THEN
                DO  i = 1, nbgp   
                   topo(:,:,nxl-i)   = topo(:,:,nxl)
                ENDDO
             ENDIF
             IF ( nxr == nx )  THEN 
                DO  i = 1, nbgp   
                   topo(:,:,nxr+i) = topo(:,:,nxr)     
                ENDDO
             ENDIF      
          ENDIF

       CASE ( 'read_from_file' )
!
!--       Note, topography information have been already read.  
!--       If required, further process topography, i.e. reference buildings on
!--       top of orography and set temporary 3D topography array, which is 
!--       used later to set grid flags. Calling of this rouinte is also 
!--       required in case of ASCII input, even though no distinction between 
!--       terrain- and building height is made in this case.  
          CALL process_topography( topo )
!
!--       Filter holes resolved by only one grid-point
          CALL filter_topography( topo )
!
!--       Exchange ghost-points, as well as add cyclic or Neumann boundary 
!--       conditions.
          CALL exchange_horiz_int( topo, nys, nyn, nxl, nxr, nzt, nbgp )
!
!--       Set lateral boundary conditions for topography on all ghost layers
          IF ( .NOT. bc_ns_cyc )  THEN
             IF ( nys == 0  )  THEN
                DO  i = 1, nbgp         
                   topo(:,nys-i,:) = topo(:,nys,:)
                ENDDO
             ENDIF
             IF ( nyn == ny )  THEN
                DO  i = 1, nbgp         
                   topo(:,nyn+i,:) = topo(:,nyn,:)
                ENDDO
             ENDIF
          ENDIF

          IF ( .NOT. bc_lr_cyc )  THEN
             IF ( nxl == 0  )  THEN
                DO  i = 1, nbgp 
                   topo(:,:,nxl-i) = topo(:,:,nxl)
                ENDDO
             ENDIF
             IF ( nxr == nx )  THEN
                DO  i = 1, nbgp 
                   topo(:,:,nxr+i) = topo(:,:,nxr)
                ENDDO
             ENDIF
          ENDIF


       CASE DEFAULT
!   
!--       The DEFAULT case is reached either if the parameter topography
!--       contains a wrong character string or if the user has defined a special
!--       case in the user interface. There, the subroutine user_init_grid 
!--       checks which of these two conditions applies.
          CALL user_init_grid( topo )
          CALL filter_topography( topo )

    END SELECT
!
!-- Consistency checks and index array initialization are only required for
!-- non-flat topography.
    IF ( TRIM( topography ) /= 'flat' )  THEN
!
!--    In case of non-flat topography, check whether the convention how to 
!--    define the topography grid has been set correctly, or whether the default
!--    is applicable. If this is not possible, abort.
       IF ( TRIM( topography_grid_convention ) == ' ' )  THEN
          IF ( TRIM( topography ) /= 'closed_channel' .AND.                    &
               TRIM( topography ) /= 'single_building' .AND.                   &
               TRIM( topography ) /= 'single_street_canyon' .AND.              &
               TRIM( topography ) /= 'tunnel'  .AND.                           &
               TRIM( topography ) /= 'read_from_file')  THEN
!--          The default value is not applicable here, because it is only valid
!--          for the four standard cases 'single_building', 
!--          'single_street_canyon', 'tunnel' and 'read_from_file'
!--          defined in init_grid.
             WRITE( message_string, * )                                        &
               'The value for "topography_grid_convention" ',                  &
               'is not set. Its default value is & only valid for ',           &
               '"topography" = ''single_building'', ''tunnel'' ',              &
               '''single_street_canyon'', ''closed_channel'' & or ',           &
               '''read_from_file''.',                                          &
               '& Choose ''cell_edge'' or ''cell_center''.'
             CALL message( 'init_grid', 'PA0239', 1, 2, 0, 6, 0 )
          ELSE
!--          The default value is applicable here.
!--          Set convention according to topography.
             IF ( TRIM( topography ) == 'single_building' .OR.                 &
                  TRIM( topography ) == 'single_street_canyon' )  THEN
                topography_grid_convention = 'cell_edge'
             ELSEIF ( TRIM( topography ) == 'read_from_file'  .OR.             &
                      TRIM( topography ) == 'tunnel')  THEN
                topography_grid_convention = 'cell_center'
             ENDIF
          ENDIF
       ELSEIF ( TRIM( topography_grid_convention ) /= 'cell_edge' .AND.        &
                TRIM( topography_grid_convention ) /= 'cell_center' )  THEN
          WRITE( message_string, * )                                           &
            'The value for "topography_grid_convention" is ',                  &
            'not recognized.& Choose ''cell_edge'' or ''cell_center''.'
          CALL message( 'init_grid', 'PA0240', 1, 2, 0, 6, 0 )
       ENDIF


       IF ( topography_grid_convention == 'cell_edge' )  THEN
! 
!--       The array nzb_local as defined using the 'cell_edge' convention 
!--       describes the actual total size of topography which is defined at the 
!--       cell edges where u=0 on the topography walls in x-direction and v=0 
!--       on the topography walls in y-direction. However, PALM uses individual
!--       arrays nzb_u|v|w|s_inner|outer that are based on nzb_s_inner.
!--       Therefore, the extent of topography in nzb_local is now reduced by 
!--       1dx at the E topography walls and by 1dy at the N topography walls 
!--       to form the basis for nzb_s_inner. 
!--       Note, the reverse memory access (i-j instead of j-i) is absolutely
!--       required at this point.
          DO  j = nys+1, nyn+1
             DO  i = nxl-1, nxr
                DO  k = nzb, nzt+1
                   IF ( BTEST( topo(k,j,i), 0 )  .OR.                          &
                        BTEST( topo(k,j,i+1), 0 ) )                            &
                       topo(k,j,i) = IBSET( topo(k,j,i), 0 )
                ENDDO
             ENDDO
          ENDDO      
          CALL exchange_horiz_int( topo, nys, nyn, nxl, nxr, nzt, nbgp )

          DO  i = nxl, nxr+1
             DO  j = nys-1, nyn
                DO  k = nzb, nzt+1
                   IF ( BTEST( topo(k,j,i), 0 )  .OR.                          &
                        BTEST( topo(k,j+1,i), 0 ) )                            &
                      topo(k,j,i) = IBSET( topo(k,j,i), 0 )
                ENDDO
             ENDDO
          ENDDO  
          CALL exchange_horiz_int( topo, nys, nyn, nxl, nxr, nzt, nbgp )
   
       ENDIF
    ENDIF


 END SUBROUTINE init_topo

 SUBROUTINE set_topo_flags(topo)

    USE control_parameters,                                                    &
        ONLY:  bc_lr_cyc, bc_ns_cyc, constant_flux_layer,                      &
               scalar_advec, topography, use_surface_fluxes, use_top_fluxes

    USE exchange_horiz_mod,                                                    &
        ONLY:  exchange_horiz_int

    USE indices,                                                               &
        ONLY:  nbgp, nx, nxl, nxlg, nxr, nxrg, ny, nyn, nyng, nys, nysg, nzb,  &
               nzt, topo_top_ind, wall_flags_static_0, wall_flags_total_0

    USE kinds

    IMPLICIT NONE

    INTEGER(iwp) ::  i             !< index variable along x
    INTEGER(iwp) ::  ibit          !< integer bit position of topgraphy masking array
    INTEGER(iwp) ::  j             !< index variable along y
    INTEGER(iwp) ::  k             !< index variable along z

    INTEGER(iwp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  topo !< input array for 3D topography and dummy array for setting "outer"-flags

    ALLOCATE( wall_flags_static_0(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    wall_flags_static_0 = 0
!
!-- Set-up topography flags. First, set flags only for s, u, v and w-grid.
!-- Further special flags will be set in following loops. 
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb, nzt+1
!
!--          scalar grid
             IF ( BTEST( topo(k,j,i), 0 ) )                                    &
                wall_flags_static_0(k,j,i) = IBSET( wall_flags_static_0(k,j,i), 0 )
!
!--          u grid
             IF ( BTEST( topo(k,j,i),   0 )  .AND.                             &
                  BTEST( topo(k,j,i-1), 0 ) )                                  &
                wall_flags_static_0(k,j,i) = IBSET( wall_flags_static_0(k,j,i), 1 )
!
!--          v grid
             IF ( BTEST( topo(k,j,i),   0 )  .AND.                             &
                  BTEST( topo(k,j-1,i), 0 ) )                                  &
                 wall_flags_static_0(k,j,i) = IBSET( wall_flags_static_0(k,j,i), 2 )

          ENDDO

          DO k = nzb, nzt
!
!--          w grid
             IF ( BTEST( topo(k,j,i),   0 )  .AND.                             &
                  BTEST( topo(k+1,j,i), 0 ) )                                  &
                wall_flags_static_0(k,j,i) = IBSET( wall_flags_static_0(k,j,i), 3 )
          ENDDO
          
          IF ( topography /= 'closed_channel' ) THEN
             wall_flags_static_0(nzt+1,j,i) = IBSET( wall_flags_static_0(nzt+1,j,i), 3 )
          ENDIF

       ENDDO
    ENDDO

    CALL exchange_horiz_int( wall_flags_static_0, nys, nyn, nxl, nxr, nzt, nbgp )

!
!-- Set outer array for scalars to mask near-surface grid points. Note, on 
!-- basis of flag 24 futher flags will be derived which are used to control
!-- production of subgrid TKE production near walls. 
    
    ALLOCATE( wall_flags_total_0(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    wall_flags_total_0 = 0
                                    
    DO i = nxl, nxr
       DO j = nys, nyn
          DO k = nzb, nzt+1
             wall_flags_total_0(k,j,i) = IOR( wall_flags_total_0(k,j,i), wall_flags_static_0(k,j,i) )
          ENDDO
       ENDDO
    ENDDO
    
    CALL exchange_horiz_int( wall_flags_total_0, nys, nyn, nxl, nxr, nzt, nbgp )
    
    DO i = nxl, nxr
       DO j = nys, nyn
          DO k = nzb, nzt+1
             IF ( BTEST( wall_flags_total_0(k,j-1,i), 0 )    .AND.                   &
                  BTEST( wall_flags_total_0(k,j+1,i), 0 )    .AND.                   &
                  BTEST( wall_flags_total_0(k,j,i-1), 0 )    .AND.                   &
                  BTEST( wall_flags_total_0(k,j,i+1), 0 )    .AND.                   &
                  BTEST( wall_flags_total_0(k,j-1,i-1), 0 )  .AND.                   &
                  BTEST( wall_flags_total_0(k,j+1,i-1), 0 )  .AND.                   &
                  BTEST( wall_flags_total_0(k,j-1,i+1), 0 )  .AND.                   &
                  BTEST( wall_flags_total_0(k,j+1,i+1), 0 ) )                        &
                wall_flags_total_0(k,j,i) = IBSET( wall_flags_total_0(k,j,i), 24 )
          ENDDO
       ENDDO
    ENDDO
!
!-- Set further special flags
    DO i = nxl, nxr
       DO j = nys, nyn
          DO k = nzb, nzt+1
!
!--          scalar grid, former nzb_diff_s_inner.
!--          Note, use this flag also to mask topography in diffusion_u and 
!--          diffusion_v along the vertical direction. In case of 
!--          use_surface_fluxes, fluxes are calculated via MOST, else, simple
!--          gradient approach is applied. Please note, in case of u- and v-
!--          diffuison, a small error is made at edges (on the east side for u,
!--          at the north side for v), since topography on scalar grid point
!--          is used instead of topography on u/v-grid. As number of topography grid
!--          points on uv-grid is different than s-grid, different number of 
!--          surface elements would be required. In order to avoid this, 
!--          treat edges (u(k,j,i+1)) simply by a gradient approach, i.e. these
!--          points are not masked within diffusion_u. Tests had shown that the
!--          effect on the flow is negligible. 
             IF ( constant_flux_layer  .OR.  use_surface_fluxes )  THEN
                IF ( BTEST( wall_flags_total_0(k,j,i), 0 ) )                         &
                   wall_flags_total_0(k,j,i) = IBSET( wall_flags_total_0(k,j,i), 8 )
             ELSE
                wall_flags_total_0(k,j,i) = IBSET( wall_flags_total_0(k,j,i), 8 )
             ENDIF

          ENDDO
!
!--       Special flag to control vertical diffusion at model top - former 
!--       nzt_diff
          wall_flags_total_0(:,j,i) = IBSET( wall_flags_total_0(:,j,i), 9 )
          IF ( use_top_fluxes )                                                &
             wall_flags_total_0(nzt+1,j,i) = IBCLR( wall_flags_total_0(nzt+1,j,i), 9 )


          DO k = nzb+1, nzt
!
!--          Special flag on u grid, former nzb_u_inner + 1, required   
!--          for disturb_field and initialization. Do not disturb directly at
!--          topography, as well as initialize u with zero one grid point outside
!--          of topography.
             IF ( BTEST( wall_flags_total_0(k-1,j,i), 1 )  .AND.                     &
                  BTEST( wall_flags_total_0(k,j,i),   1 )  .AND.                     &
                  BTEST( wall_flags_total_0(k+1,j,i), 1 ) )                          &
                wall_flags_total_0(k,j,i) = IBSET( wall_flags_total_0(k,j,i), 20 )
!
!--          Special flag on v grid, former nzb_v_inner + 1, required   
!--          for disturb_field and initialization. Do not disturb directly at
!--          topography, as well as initialize v with zero one grid point outside
!--          of topography.
             IF ( BTEST( wall_flags_total_0(k-1,j,i), 2 )  .AND.                     &
                  BTEST( wall_flags_total_0(k,j,i),   2 )  .AND.                     &
                  BTEST( wall_flags_total_0(k+1,j,i), 2 ) )                          &
                wall_flags_total_0(k,j,i) = IBSET( wall_flags_total_0(k,j,i), 21 )
!
!--          Special flag on scalar grid, former nzb_s_inner+1. Used for 
!--          lpm_sgs_tke
             IF ( BTEST( wall_flags_total_0(k,j,i),   0 )  .AND.                     &
                  BTEST( wall_flags_total_0(k-1,j,i), 0 )  .AND.                     &
                  BTEST( wall_flags_total_0(k+1,j,i), 0 ) )                          &
                wall_flags_total_0(k,j,i) = IBSET( wall_flags_total_0(k,j,i), 25 )
!
!--          Special flag on scalar grid, nzb_diff_s_outer - 1, required in 
!--          in production_e
             IF ( constant_flux_layer  .OR.  use_surface_fluxes )  THEN
                IF ( BTEST( wall_flags_total_0(k,j,i),   24 )  .AND.                 &
                     BTEST( wall_flags_total_0(k-1,j,i), 24 )  .AND.                 &
                     BTEST( wall_flags_total_0(k+1,j,i), 0 ) )                       &
                   wall_flags_total_0(k,j,i) = IBSET( wall_flags_total_0(k,j,i), 29 )
             ELSE
                IF ( BTEST( wall_flags_total_0(k,j,i), 0 ) )                         &
                   wall_flags_total_0(k,j,i) = IBSET( wall_flags_total_0(k,j,i), 29 )
             ENDIF
!
!--          Special flag on scalar grid, nzb_diff_s_outer - 1, required in 
!--          in production_e
             IF ( constant_flux_layer  .OR.  use_surface_fluxes )  THEN
                IF ( BTEST( wall_flags_total_0(k,j,i),   0 )  .AND.                  &
                     BTEST( wall_flags_total_0(k-1,j,i), 0 )  .AND.                  &
                     BTEST( wall_flags_total_0(k+1,j,i), 0 ) )                       &
                   wall_flags_total_0(k,j,i) = IBSET( wall_flags_total_0(k,j,i), 30 )
             ELSE
                IF ( BTEST( wall_flags_total_0(k,j,i), 0 ) )                         &
                   wall_flags_total_0(k,j,i) = IBSET( wall_flags_total_0(k,j,i), 30 )
             ENDIF
          ENDDO
!
!--       Flags indicating downward facing walls
          DO k = nzb+1, nzt+1
!
!--          Scalar grid
             IF ( BTEST( wall_flags_total_0(k-1,j,i), 0 )  .AND.                     &
            .NOT. BTEST( wall_flags_total_0(k,j,i), 0   ) )                          & 
                 wall_flags_total_0(k,j,i) = IBSET( wall_flags_total_0(k,j,i), 13 ) 
!
!--          Downward facing wall on u grid
             IF ( BTEST( wall_flags_total_0(k-1,j,i), 1 )  .AND.                     &
            .NOT. BTEST( wall_flags_total_0(k,j,i), 1   ) )                          & 
                wall_flags_total_0(k,j,i) = IBSET( wall_flags_total_0(k,j,i), 15 )
!
!--          Downward facing wall on v grid
             IF ( BTEST( wall_flags_total_0(k-1,j,i), 2 )  .AND.                     &
            .NOT. BTEST( wall_flags_total_0(k,j,i), 2   ) )                          & 
                wall_flags_total_0(k,j,i) = IBSET( wall_flags_total_0(k,j,i), 17 )
!
!--          Downward facing wall on w grid
             IF ( BTEST( wall_flags_total_0(k-1,j,i), 3 )  .AND.                     &
            .NOT. BTEST( wall_flags_total_0(k,j,i), 3 ) )                            & 
                wall_flags_total_0(k,j,i) = IBSET( wall_flags_total_0(k,j,i), 19 )
          ENDDO
!
!--       Flags indicating upward facing walls
          DO k = nzb, nzt
!
!--          Upward facing wall on scalar grid
             IF ( .NOT. BTEST( wall_flags_total_0(k,j,i),   0 )  .AND.               &
                        BTEST( wall_flags_total_0(k+1,j,i), 0 ) )                    & 
                wall_flags_total_0(k,j,i) = IBSET( wall_flags_total_0(k,j,i), 12 )
!
!--          Upward facing wall on u grid
             IF ( .NOT. BTEST( wall_flags_total_0(k,j,i),   1 )  .AND.               &
                        BTEST( wall_flags_total_0(k+1,j,i), 1 ) )                    & 
                wall_flags_total_0(k,j,i) = IBSET( wall_flags_total_0(k,j,i), 14 )

!   
!--          Upward facing wall on v grid
             IF ( .NOT. BTEST( wall_flags_total_0(k,j,i),   2 )  .AND.               &
                        BTEST( wall_flags_total_0(k+1,j,i), 2 ) )                    & 
                wall_flags_total_0(k,j,i) = IBSET( wall_flags_total_0(k,j,i), 16 )
   
!
!--          Upward facing wall on w grid
             IF ( .NOT. BTEST( wall_flags_total_0(k,j,i),   3 )  .AND.               &
                        BTEST( wall_flags_total_0(k+1,j,i), 3 ) )                    & 
                wall_flags_total_0(k,j,i) = IBSET( wall_flags_total_0(k,j,i), 18 )
!
!--          Special flag on scalar grid, former nzb_s_inner
             IF ( BTEST( wall_flags_total_0(k,j,i), 0 )  .OR.                        &
                  BTEST( wall_flags_total_0(k,j,i), 12 ) .OR.                        &
                  BTEST( wall_flags_total_0(k,j,i), 13 ) )                           &
                   wall_flags_total_0(k,j,i) = IBSET( wall_flags_total_0(k,j,i), 22 )
!
!--          Special flag on scalar grid, nzb_diff_s_inner - 1, required for 
!--          flow_statistics
             IF ( constant_flux_layer  .OR.  use_surface_fluxes )  THEN
                IF ( BTEST( wall_flags_total_0(k,j,i),   0 )  .AND.                  &
                     BTEST( wall_flags_total_0(k+1,j,i), 0 ) )                       &
                  wall_flags_total_0(k,j,i) = IBSET( wall_flags_total_0(k,j,i), 23 )
             ELSE
                IF ( BTEST( wall_flags_total_0(k,j,i), 22 ) )                        &
                   wall_flags_total_0(k,j,i) = IBSET( wall_flags_total_0(k,j,i), 23 )
             ENDIF
   

          ENDDO
          wall_flags_total_0(nzt+1,j,i) = IBSET( wall_flags_total_0(nzt+1,j,i), 22 )
          wall_flags_total_0(nzt+1,j,i) = IBSET( wall_flags_total_0(nzt+1,j,i), 23 )
!
!--       Set flags indicating that topography is close by in horizontal 
!--       direction, i.e. flags that infold the topography. These will be used
!--       to set advection flags for passive scalars, where due to large 
!--       gradients near buildings stationary numerical oscillations can produce 
!--       unrealistically high concentrations. This is only necessary if 
!--       WS-scheme is applied for scalar advection. Note, these flags will be 
!--       only used for passive scalars such as chemical species or aerosols.
          IF ( scalar_advec == 'ws-scheme' )  THEN
             DO k = nzb, nzt
                IF ( BTEST( wall_flags_total_0(k,j,i), 0 )  .AND. (                  &
                     ANY( .NOT. BTEST( wall_flags_total_0(k,j-3:j+3,i-1), 0 ) )  .OR.&
                     ANY( .NOT. BTEST( wall_flags_total_0(k,j-3:j+3,i-2), 0 ) )  .OR.&
                     ANY( .NOT. BTEST( wall_flags_total_0(k,j-3:j+3,i-3), 0 ) )  .OR.&
                     ANY( .NOT. BTEST( wall_flags_total_0(k,j-3:j+3,i+1), 0 ) )  .OR.&
                     ANY( .NOT. BTEST( wall_flags_total_0(k,j-3:j+3,i+2), 0 ) )  .OR.&
                     ANY( .NOT. BTEST( wall_flags_total_0(k,j-3:j+3,i+3), 0 ) )  .OR.&
                     ANY( .NOT. BTEST( wall_flags_total_0(k,j-1,i-3:i+3), 0 ) )  .OR.&
                     ANY( .NOT. BTEST( wall_flags_total_0(k,j-2,i-3:i+3), 0 ) )  .OR.&
                     ANY( .NOT. BTEST( wall_flags_total_0(k,j-3,i-3:i+3), 0 ) )  .OR.&
                     ANY( .NOT. BTEST( wall_flags_total_0(k,j+1,i-3:i+3), 0 ) )  .OR.&
                     ANY( .NOT. BTEST( wall_flags_total_0(k,j+2,i-3:i+3), 0 ) )  .OR.&
                     ANY( .NOT. BTEST( wall_flags_total_0(k,j+3,i-3:i+3), 0 ) )      &
                                                            ) )                      &
                   wall_flags_total_0(k,j,i) = IBSET( wall_flags_total_0(k,j,i), 31 )
                     
             ENDDO
          ENDIF
       ENDDO
    ENDDO
!
!-- Finally, set identification flags indicating natural terrain or buildings.
!-- Natural terrain grid points. Information on the type of the surface is 
!-- stored in bit 1 of 3D Integer array topo. However, this bit is only set 
!-- when topography is read from file. In order to run the land-surface model
!-- also without topography information, set bit 1 explicitely in this case.
!-- 
!-- Natural terrain grid points
!-- If no topography is initialized, the land-surface is at k = nzb.
    IF ( TRIM( topography ) /= 'read_from_file' )  THEN
       wall_flags_static_0(nzb,:,:) = IBSET( wall_flags_static_0(nzb,:,:), 5 )
    ELSE
       DO i = nxl, nxr
          DO j = nys, nyn
             DO k = nzb, nzt+1
!         
!--             Natural terrain grid point
                IF ( BTEST( topo(k,j,i), 1 ) )                                 &
                   wall_flags_static_0(k,j,i) = IBSET( wall_flags_static_0(k,j,i), 5 )
             ENDDO
          ENDDO
       ENDDO
    ENDIF
!
!-- Building grid points.
    DO i = nxl, nxr
       DO j = nys, nyn
          DO k = nzb, nzt+1
             IF ( BTEST( topo(k,j,i), 2 ) )                                    &
                wall_flags_static_0(k,j,i) = IBSET( wall_flags_static_0(k,j,i), 6 )
          ENDDO
       ENDDO
    ENDDO
!
!-- Set flag 4, indicating new topography grid points due to filtering.
    DO i = nxl, nxr
       DO j = nys, nyn
          DO k = nzb, nzt+1
             IF ( BTEST( topo(k,j,i), 4 ) )                                    &
                wall_flags_static_0(k,j,i) = IBSET( wall_flags_static_0(k,j,i), 4 )
          ENDDO
       ENDDO
    ENDDO
    
    CALL exchange_horiz_int( wall_flags_static_0, nys, nyn, nxl, nxr, nzt, nbgp )
    
    DO i = nxl, nxr
       DO j = nys, nyn
          DO k = nzb, nzt+1
             wall_flags_total_0(k,j,i) = IOR( wall_flags_total_0(k,j,i), wall_flags_static_0(k,j,i) )
          ENDDO
       ENDDO
    ENDDO
!
!-- Exchange ghost points for wall flags
    CALL exchange_horiz_int( wall_flags_total_0, nys, nyn, nxl, nxr, nzt, nbgp )
!
!-- Set boundary conditions also for flags. Can be interpreted as Neumann
!-- boundary conditions for topography. 
    IF ( .NOT. bc_ns_cyc )  THEN
       IF ( nys == 0  )  THEN
          DO  i = 1, nbgp     
             wall_flags_total_0(:,nys-i,:)   = wall_flags_total_0(:,nys,:)
          ENDDO
       ENDIF
       IF ( nyn == ny )  THEN
          DO  i = 1, nbgp  
             wall_flags_total_0(:,nyn+i,:) = wall_flags_total_0(:,nyn,:)
          ENDDO
       ENDIF
    ENDIF
    IF ( .NOT. bc_lr_cyc )  THEN
       IF ( nxl == 0  )  THEN
          DO  i = 1, nbgp   
             wall_flags_total_0(:,:,nxl-i)   = wall_flags_total_0(:,:,nxl)
          ENDDO
       ENDIF
       IF ( nxr == nx )  THEN 
          DO  i = 1, nbgp   
             wall_flags_total_0(:,:,nxr+i) = wall_flags_total_0(:,:,nxr)     
          ENDDO
       ENDIF      
    ENDIF
!
!-- Pre-calculate topography top indices (former get_topography_top_index 
!-- function)
    ALLOCATE( topo_top_ind(nysg:nyng,nxlg:nxrg,0:4) )
!
!-- Uppermost topography index on scalar grid
    ibit = 12
    topo_top_ind(:,:,0) = MAXLOC(                                              &
                                  MERGE( 1, 0,                                 &
                                    BTEST( wall_flags_total_0(:,:,:), ibit )   &
                                       ), DIM = 1                              &
                                ) - 1 
!
!-- Uppermost topography index on u grid 
    ibit = 14
    topo_top_ind(:,:,1) = MAXLOC(                                              &
                                  MERGE( 1, 0,                                 &
                                    BTEST( wall_flags_total_0(:,:,:), ibit )   &
                                       ), DIM = 1                              &
                                ) - 1 
!
!-- Uppermost topography index on v grid 
    ibit = 16
    topo_top_ind(:,:,2) = MAXLOC(                                              &
                                  MERGE( 1, 0,                                 &
                                    BTEST( wall_flags_total_0(:,:,:), ibit )   &
                                       ), DIM = 1                              &
                                ) - 1 
!
!-- Uppermost topography index on w grid
    ibit = 18
    topo_top_ind(:,:,3) = MAXLOC(                                              &
                                  MERGE( 1, 0,                                 &
                                    BTEST( wall_flags_total_0(:,:,:), ibit )   &
                                       ), DIM = 1                              &
                                ) - 1 
!
!-- Uppermost topography index on scalar outer grid
    ibit = 24
    topo_top_ind(:,:,4) = MAXLOC(                                              &
                                  MERGE( 1, 0,                                 &
                                    BTEST( wall_flags_total_0(:,:,:), ibit )   &
                                       ), DIM = 1                              &
                                ) - 1

 END SUBROUTINE set_topo_flags




