!> @file surface_coupler.f90
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
! $Id: surface_coupler.f90 4429 2020-02-27 15:24:30Z raasch $
! bugfix: preprocessor directives rearranged for serial mode
! 
! 4360 2020-01-07 11:25:50Z suehring
! Corrected "Former revisions" section
! 
! 3655 2019-01-07 16:51:22Z knoop
! Modularization of all bulk cloud physics code components
!
! 109 2007-08-28 15:26:47Z letzel
! Initial revision
!
! Description:
! ------------
!> Data exchange at the interface between coupled models
!------------------------------------------------------------------------------!
 SUBROUTINE surface_coupler
#if defined( __parallel )
 

    USE arrays_3d,                                                             &
        ONLY:  pt, rho_ocean, sa, total_2d_a, total_2d_o, u, v

    USE basic_constants_and_equations_mod,                                     &
        ONLY:  c_p, l_v

    USE control_parameters,                                                    &
        ONLY:  coupling_mode, coupling_mode_remote, coupling_topology,         &
               humidity, humidity_remote, land_surface, message_string,        &
               terminate_coupled, terminate_coupled_remote,                    &
               time_since_reference_point, urban_surface

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point

    USE indices,                                                               &
        ONLY:  nbgp, nx, nxl, nxlg, nxr, nxrg, nx_a, nx_o, ny, nyn, nyng, nys, &
               nysg, ny_a, ny_o, nzt

    USE kinds

    USE pegrid

    USE surface_mod,                                                           &
        ONLY :  surf_def_h, surf_lsm_h, surf_type, surf_usm_h

    IMPLICIT NONE

    INTEGER(iwp) ::  i                                    !< index variable x-direction 
    INTEGER(iwp) ::  j                                    !< index variable y-direction 
    INTEGER(iwp) ::  m                                    !< running index for surface elements

    REAL(wp)    ::  cpw = 4218.0_wp                       !< heat capacity of water at constant pressure
    REAL(wp)    ::  time_since_reference_point_rem        !< 
    REAL(wp)    ::  total_2d(-nbgp:ny+nbgp,-nbgp:nx+nbgp) !< 

    REAL(wp), DIMENSION(nysg:nyng,nxlg:nxrg) ::  surface_flux !< dummy array for surface fluxes on 2D grid


    CALL cpu_log( log_point(39), 'surface_coupler', 'start' )



!
!-- In case of model termination initiated by the remote model 
!-- (terminate_coupled_remote > 0), initiate termination of the local model. 
!-- The rest of the coupler must then be skipped because it would cause an MPI 
!-- intercomminucation hang.
!-- If necessary, the coupler will be called at the beginning of the next 
!-- restart run.

    IF ( coupling_topology == 0 ) THEN
       CALL MPI_SENDRECV( terminate_coupled,        1, MPI_INTEGER, target_id, &
                          0,                                                   &
                          terminate_coupled_remote, 1, MPI_INTEGER, target_id, &
                          0, comm_inter, status, ierr )
    ELSE
       IF ( myid == 0) THEN
          CALL MPI_SENDRECV( terminate_coupled,        1, MPI_INTEGER, &
                             target_id, 0,                             &
                             terminate_coupled_remote, 1, MPI_INTEGER, & 
                             target_id, 0,                             &
                             comm_inter, status, ierr )
       ENDIF
       CALL MPI_BCAST( terminate_coupled_remote, 1, MPI_INTEGER, 0, comm2d, &
                       ierr )

       ALLOCATE( total_2d_a(-nbgp:ny_a+nbgp,-nbgp:nx_a+nbgp),       &
                 total_2d_o(-nbgp:ny_o+nbgp,-nbgp:nx_o+nbgp) )

    ENDIF

    IF ( terminate_coupled_remote > 0 )  THEN
       WRITE( message_string, * ) 'remote model "',                            &
                                  TRIM( coupling_mode_remote ),                &
                                  '" terminated',                              &
                                  '&with terminate_coupled_remote = ',         &
                                  terminate_coupled_remote,                    &
                                  '&local model  "', TRIM( coupling_mode ),    &
                                  '" has',                                     &
                                  '&terminate_coupled = ',                     &
                                   terminate_coupled
       CALL message( 'surface_coupler', 'PA0310', 1, 2, 0, 6, 0 )
       RETURN
    ENDIF
 

!
!-- Exchange the current simulated time between the models,
!-- currently just for total_2d
    IF ( coupling_topology == 0 ) THEN
   
       CALL MPI_SEND( time_since_reference_point, 1, MPI_REAL, target_id, 11, &
                      comm_inter, ierr )
       CALL MPI_RECV( time_since_reference_point_rem, 1, MPI_REAL, target_id, &
                      11, comm_inter, status, ierr )
    ELSE

       IF ( myid == 0 ) THEN

          CALL MPI_SEND( time_since_reference_point, 1, MPI_REAL, target_id, &
                         11, comm_inter, ierr )
          CALL MPI_RECV( time_since_reference_point_rem, 1, MPI_REAL,        &
                         target_id, 11, comm_inter, status, ierr )

       ENDIF

       CALL MPI_BCAST( time_since_reference_point_rem, 1, MPI_REAL, 0, comm2d, &
                       ierr )

    ENDIF

!
!-- Exchange the interface data
    IF ( coupling_mode == 'atmosphere_to_ocean' )  THEN
    
!
!--    Horizontal grid size and number of processors is equal in ocean and
!--    atmosphere
       IF ( coupling_topology == 0 )  THEN

!
!--       Send heat flux at bottom surface to the ocean. First, transfer from 
!--       1D surface type to 2D grid.
          CALL transfer_1D_to_2D_equal( surf_def_h(0)%shf, surf_lsm_h%shf,     &
                                        surf_usm_h%shf )
          CALL MPI_SEND( surface_flux(nysg,nxlg), ngp_xy, MPI_REAL, target_id, &
                         12, comm_inter, ierr )
!
!--       Send humidity flux at bottom surface to the ocean. First, transfer
!--       from 1D surface type to 2D grid.
          CALL transfer_1D_to_2D_equal( surf_def_h(0)%qsws, surf_lsm_h%qsws,   &
                                        surf_usm_h%qsws )
          IF ( humidity )  THEN
             CALL MPI_SEND( surface_flux(nysg,nxlg), ngp_xy, MPI_REAL,         &
                            target_id, 13, comm_inter, ierr )
          ENDIF
!
!--       Receive temperature at the bottom surface from the ocean
          CALL MPI_RECV( pt(0,nysg,nxlg), 1, type_xy, target_id, 14,           &
                         comm_inter, status, ierr )
!
!--       Send the momentum flux (u) at bottom surface to the ocean. First, 
!--       transfer from 1D surface type to 2D grid.
          CALL transfer_1D_to_2D_equal( surf_def_h(0)%usws, surf_lsm_h%usws,   &
                                        surf_usm_h%usws )
          CALL MPI_SEND( surface_flux(nysg,nxlg), ngp_xy, MPI_REAL, target_id, &
                         15, comm_inter, ierr )
!
!--       Send the momentum flux (v) at bottom surface to the ocean. First, 
!--       transfer from 1D surface type to 2D grid.
          CALL transfer_1D_to_2D_equal( surf_def_h(0)%vsws, surf_lsm_h%vsws,   &
                                        surf_usm_h%vsws )
          CALL MPI_SEND( surface_flux(nysg,nxlg), ngp_xy, MPI_REAL, target_id, &
                         16, comm_inter, ierr )
!
!--       Receive u at the bottom surface from the ocean
          CALL MPI_RECV( u(0,nysg,nxlg), 1, type_xy, target_id, 17,            &
                         comm_inter, status, ierr )
!
!--       Receive v at the bottom surface from the ocean
          CALL MPI_RECV( v(0,nysg,nxlg), 1, type_xy, target_id, 18,            &
                         comm_inter, status, ierr )
!
!--    Horizontal grid size or number of processors differs between 
!--    ocean and atmosphere
       ELSE
      
!
!--       Send heat flux at bottom surface to the ocean
          total_2d_a = 0.0_wp
          total_2d   = 0.0_wp
!
!--       Transfer from 1D surface type to 2D grid.
          CALL transfer_1D_to_2D_unequal( surf_def_h(0)%shf, surf_lsm_h%shf,   &
                                          surf_usm_h%shf )

          CALL MPI_REDUCE( total_2d, total_2d_a, ngp_a, MPI_REAL, MPI_SUM, 0,  &
                           comm2d, ierr )
          CALL interpolate_to_ocean( 12 )   
!
!--       Send humidity flux at bottom surface to the ocean
          IF ( humidity )  THEN
             total_2d_a = 0.0_wp
             total_2d   = 0.0_wp
!
!--          Transfer from 1D surface type to 2D grid.
             CALL transfer_1D_to_2D_unequal( surf_def_h(0)%qsws,              &
                                             surf_lsm_h%qsws,                 &
                                             surf_usm_h%qsws )

             CALL MPI_REDUCE( total_2d, total_2d_a, ngp_a, MPI_REAL, MPI_SUM, &
                              0, comm2d, ierr )
             CALL interpolate_to_ocean( 13 )
          ENDIF
!
!--       Receive temperature at the bottom surface from the ocean
          IF ( myid == 0 )  THEN
             CALL MPI_RECV( total_2d_a(-nbgp,-nbgp), ngp_a, MPI_REAL,          &
                            target_id, 14, comm_inter, status, ierr )   
          ENDIF
          CALL MPI_BARRIER( comm2d, ierr )
          CALL MPI_BCAST( total_2d_a(-nbgp,-nbgp), ngp_a, MPI_REAL, 0, comm2d, &
                          ierr )
          pt(0,nysg:nyng,nxlg:nxrg) = total_2d_a(nysg:nyng,nxlg:nxrg)
!
!--       Send momentum flux (u) at bottom surface to the ocean
          total_2d_a = 0.0_wp 
          total_2d   = 0.0_wp
!
!--       Transfer from 1D surface type to 2D grid.
          CALL transfer_1D_to_2D_unequal( surf_def_h(0)%usws, surf_lsm_h%usws, &
                                          surf_usm_h%usws )
          CALL MPI_REDUCE( total_2d, total_2d_a, ngp_a, MPI_REAL, MPI_SUM, 0, &
                           comm2d, ierr )
          CALL interpolate_to_ocean( 15 )
!
!--       Send momentum flux (v) at bottom surface to the ocean
          total_2d_a = 0.0_wp
          total_2d   = 0.0_wp
!
!--       Transfer from 1D surface type to 2D grid.
          CALL transfer_1D_to_2D_unequal( surf_def_h(0)%usws, surf_lsm_h%usws, &
                                          surf_usm_h%usws )
          CALL MPI_REDUCE( total_2d, total_2d_a, ngp_a, MPI_REAL, MPI_SUM, 0, &
                           comm2d, ierr )
          CALL interpolate_to_ocean( 16 )
!
!--       Receive u at the bottom surface from the ocean
          IF ( myid == 0 )  THEN
             CALL MPI_RECV( total_2d_a(-nbgp,-nbgp), ngp_a, MPI_REAL, &
                            target_id, 17, comm_inter, status, ierr )
          ENDIF
          CALL MPI_BARRIER( comm2d, ierr )
          CALL MPI_BCAST( total_2d_a(-nbgp,-nbgp), ngp_a, MPI_REAL, 0, comm2d, &
                          ierr )
          u(0,nysg:nyng,nxlg:nxrg) = total_2d_a(nysg:nyng,nxlg:nxrg)
!
!--       Receive v at the bottom surface from the ocean
          IF ( myid == 0 )  THEN
             CALL MPI_RECV( total_2d_a(-nbgp,-nbgp), ngp_a, MPI_REAL, &
                            target_id, 18, comm_inter, status, ierr )
          ENDIF
          CALL MPI_BARRIER( comm2d, ierr )
          CALL MPI_BCAST( total_2d_a(-nbgp,-nbgp), ngp_a, MPI_REAL, 0, comm2d, &
                          ierr )
          v(0,nysg:nyng,nxlg:nxrg) = total_2d_a(nysg:nyng,nxlg:nxrg)

       ENDIF

    ELSEIF ( coupling_mode == 'ocean_to_atmosphere' )  THEN

!
!--    Horizontal grid size and number of processors is equal
!--    in ocean and atmosphere
       IF ( coupling_topology == 0 ) THEN
!
!--       Receive heat flux at the sea surface (top) from the atmosphere
          CALL MPI_RECV( surface_flux(nysg,nxlg), ngp_xy, MPI_REAL, target_id, 12, &
                         comm_inter, status, ierr )
          CALL transfer_2D_to_1D_equal( surf_def_h(2)%shf )
!
!--       Receive humidity flux from the atmosphere (bottom)
!--       and add it to the heat flux at the sea surface (top)...
          IF ( humidity_remote )  THEN
             CALL MPI_RECV( surface_flux(nysg,nxlg), ngp_xy, MPI_REAL, &
                            target_id, 13, comm_inter, status, ierr )
             CALL transfer_2D_to_1D_equal( surf_def_h(2)%qsws )
          ENDIF
!
!--       Send sea surface temperature to the atmosphere model
          CALL MPI_SEND( pt(nzt,nysg,nxlg), 1, type_xy, target_id, 14, &
                         comm_inter, ierr )
!
!--       Receive momentum flux (u) at the sea surface (top) from the atmosphere
          CALL MPI_RECV( surface_flux(nysg,nxlg), ngp_xy, MPI_REAL, target_id, 15, &
                         comm_inter, status, ierr )
          CALL transfer_2D_to_1D_equal( surf_def_h(2)%usws )
!
!--       Receive momentum flux (v) at the sea surface (top) from the atmosphere
          CALL MPI_RECV( surface_flux(nysg,nxlg), ngp_xy, MPI_REAL, target_id, 16, &
                         comm_inter, status, ierr )
          CALL transfer_2D_to_1D_equal( surf_def_h(2)%vsws )
!
!--       Send u to the atmosphere
          CALL MPI_SEND( u(nzt,nysg,nxlg), 1, type_xy, target_id, 17, &
                         comm_inter, ierr )
!
!--       Send v to the atmosphere
          CALL MPI_SEND( v(nzt,nysg,nxlg), 1, type_xy, target_id, 18, &
                         comm_inter, ierr )
!
!--    Horizontal gridsize or number of processors differs between 
!--    ocean and atmosphere
       ELSE
!
!--       Receive heat flux at the sea surface (top) from the atmosphere
          IF ( myid == 0 )  THEN
             CALL MPI_RECV( total_2d_o(-nbgp,-nbgp), ngp_o, MPI_REAL, &
                            target_id, 12, comm_inter, status, ierr )
          ENDIF
          CALL MPI_BARRIER( comm2d, ierr )
          CALL MPI_BCAST( total_2d_o(-nbgp,-nbgp), ngp_o, MPI_REAL, 0, comm2d, &
                          ierr )
          CALL transfer_2D_to_1D_unequal( surf_def_h(2)%shf )
!
!--       Receive humidity flux at the sea surface (top) from the atmosphere
          IF ( humidity_remote )  THEN
             IF ( myid == 0 )  THEN
                CALL MPI_RECV( total_2d_o(-nbgp,-nbgp), ngp_o, MPI_REAL, &
                               target_id, 13, comm_inter, status, ierr )
             ENDIF
             CALL MPI_BARRIER( comm2d, ierr )
             CALL MPI_BCAST( total_2d_o(-nbgp,-nbgp), ngp_o, MPI_REAL, 0, &
                             comm2d, ierr)
             CALL transfer_2D_to_1D_unequal( surf_def_h(2)%qsws )
          ENDIF
!
!--       Send surface temperature to atmosphere
          total_2d_o = 0.0_wp
          total_2d   = 0.0_wp
          total_2d(nys:nyn,nxl:nxr) = pt(nzt,nys:nyn,nxl:nxr)

          CALL MPI_REDUCE( total_2d, total_2d_o, ngp_o, MPI_REAL, MPI_SUM, 0, &
                           comm2d, ierr) 
          CALL interpolate_to_atmos( 14 )
!
!--       Receive momentum flux (u) at the sea surface (top) from the atmosphere
          IF ( myid == 0 )  THEN
             CALL MPI_RECV( total_2d_o(-nbgp,-nbgp), ngp_o, MPI_REAL, &
                            target_id, 15, comm_inter, status, ierr )
          ENDIF
          CALL MPI_BARRIER( comm2d, ierr )
          CALL MPI_BCAST( total_2d_o(-nbgp,-nbgp), ngp_o, MPI_REAL, &
                          0, comm2d, ierr )
          CALL transfer_2D_to_1D_unequal( surf_def_h(2)%usws )
!
!--       Receive momentum flux (v) at the sea surface (top) from the atmosphere
          IF ( myid == 0 )  THEN
             CALL MPI_RECV( total_2d_o(-nbgp,-nbgp), ngp_o, MPI_REAL, &
                            target_id, 16, comm_inter, status, ierr )
          ENDIF
          CALL MPI_BARRIER( comm2d, ierr )
          CALL MPI_BCAST( total_2d_o(-nbgp,-nbgp), ngp_o, MPI_REAL, 0, comm2d, &
                          ierr )
          CALL transfer_2D_to_1D_unequal( surf_def_h(2)%vsws )
!
!--       Send u to atmosphere
          total_2d_o = 0.0_wp  
          total_2d   = 0.0_wp
          total_2d(nys:nyn,nxl:nxr) = u(nzt,nys:nyn,nxl:nxr)
          CALL MPI_REDUCE( total_2d, total_2d_o, ngp_o, MPI_REAL, MPI_SUM, 0, &
                           comm2d, ierr )
          CALL interpolate_to_atmos( 17 )
!
!--       Send v to atmosphere
          total_2d_o = 0.0_wp
          total_2d   = 0.0_wp
          total_2d(nys:nyn,nxl:nxr) = v(nzt,nys:nyn,nxl:nxr)
          CALL MPI_REDUCE( total_2d, total_2d_o, ngp_o, MPI_REAL, MPI_SUM, 0, &
                           comm2d, ierr )
          CALL interpolate_to_atmos( 18 )
        
       ENDIF

!
!--    Conversions of fluxes received from atmosphere
       IF ( humidity_remote )  THEN
!
!--       Here top heat flux is still the sum of atmospheric bottom heat fluxes,
!--       * latent heat of vaporization in m2/s2, or 540 cal/g, or 40.65 kJ/mol
!--       /(rho_atm(=1.0)*c_p)
          DO  m = 1, surf_def_h(2)%ns
             i = surf_def_h(2)%i(m)
             j = surf_def_h(2)%j(m)
             
             surf_def_h(2)%shf(m) = surf_def_h(2)%shf(m) +                     &
                                    surf_def_h(2)%qsws(m) * l_v / c_p
!
!--          ...and convert it to a salinity flux at the sea surface (top)
!--          following Steinhorn (1991), JPO 21, pp. 1681-1683:
!--          S'w' = -S * evaporation / ( rho_water * ( 1 - S ) )
             surf_def_h(2)%sasws(m) = -1.0_wp * sa(nzt,j,i) * 0.001_wp *       &
                                      surf_def_h(2)%qsws(m) /                  &
                                    ( rho_ocean(nzt,j,i) *                     &
                                      ( 1.0_wp - sa(nzt,j,i) * 0.001_wp )      &
                                    )
          ENDDO
       ENDIF

!
!--    Adjust the kinematic heat flux with respect to ocean density
!--    (constants are the specific heat capacities for air and water), as well 
!--    as momentum fluxes 
       DO  m = 1, surf_def_h(2)%ns
          i = surf_def_h(2)%i(m)
          j = surf_def_h(2)%j(m)
          surf_def_h(2)%shf(m) = surf_def_h(2)%shf(m) / rho_ocean(nzt,j,i) *   &
                                 c_p / cpw

          surf_def_h(2)%usws(m) = surf_def_h(2)%usws(m) / rho_ocean(nzt,j,i)
          surf_def_h(2)%vsws(m) = surf_def_h(2)%vsws(m) / rho_ocean(nzt,j,i)
       ENDDO

    ENDIF

    IF ( coupling_topology == 1 )  THEN
       DEALLOCATE( total_2d_o, total_2d_a )
    ENDIF

    CALL cpu_log( log_point(39), 'surface_coupler', 'stop' )


     CONTAINS 

!       Description:
!------------------------------------------------------------------------------!
!>      Data transfer from 1D surface-data type to 2D dummy array for equal 
!>      grids in atmosphere and ocean.
!------------------------------------------------------------------------------!
        SUBROUTINE transfer_1D_to_2D_equal( def_1d, lsm_1d, usm_1d )

           IMPLICIT NONE

            INTEGER(iwp) ::  i   !< running index x
            INTEGER(iwp) ::  j   !< running index y
            INTEGER(iwp) ::  m   !< running index surface type

            REAL(wp), DIMENSION(1:surf_def_h(0)%ns) ::  def_1d !< 1D surface flux, default surfaces
            REAL(wp), DIMENSION(1:surf_lsm_h%ns)    ::  lsm_1d !< 1D surface flux, natural surfaces
            REAL(wp), DIMENSION(1:surf_usm_h%ns)    ::  usm_1d !< 1D surface flux, urban surfaces
!
!--         Transfer surface flux at default surfaces to 2D grid 
            DO  m = 1, surf_def_h(0)%ns
               i = surf_def_h(0)%i(m)
               j = surf_def_h(0)%j(m)
               surface_flux(j,i) = def_1d(m)
            ENDDO
!
!--         Transfer surface flux at natural surfaces to 2D grid 
            IF ( land_surface )  THEN 
               DO  m = 1, SIZE(lsm_1d)
                  i = surf_lsm_h%i(m)
                  j = surf_lsm_h%j(m)
                  surface_flux(j,i) = lsm_1d(m)
               ENDDO
            ENDIF
!
!--         Transfer surface flux at natural surfaces to 2D grid 
            IF ( urban_surface )  THEN 
               DO  m = 1, SIZE(usm_1d)
                  i = surf_usm_h%i(m)
                  j = surf_usm_h%j(m)
                  surface_flux(j,i) = usm_1d(m)
               ENDDO
            ENDIF

        END SUBROUTINE transfer_1D_to_2D_equal

!       Description:
!------------------------------------------------------------------------------!
!>      Data transfer from 2D array for equal grids onto 1D surface-data type
!>      array.
!------------------------------------------------------------------------------!
        SUBROUTINE transfer_2D_to_1D_equal( def_1d )

           IMPLICIT NONE

            INTEGER(iwp) ::  i   !< running index x
            INTEGER(iwp) ::  j   !< running index y
            INTEGER(iwp) ::  m   !< running index surface type

            REAL(wp), DIMENSION(1:surf_def_h(2)%ns) ::  def_1d !< 1D surface flux, default surfaces
!
!--         Transfer surface flux to 1D surface type, only for default surfaces 
            DO  m = 1, surf_def_h(2)%ns
               i = surf_def_h(2)%i(m)
               j = surf_def_h(2)%j(m)
               def_1d(m) = surface_flux(j,i)
            ENDDO

        END SUBROUTINE transfer_2D_to_1D_equal

!       Description:
!------------------------------------------------------------------------------!
!>      Data transfer from 1D surface-data type to 2D dummy array from unequal 
!>      grids in atmosphere and ocean.
!------------------------------------------------------------------------------!
        SUBROUTINE transfer_1D_to_2D_unequal( def_1d, lsm_1d, usm_1d )

           IMPLICIT NONE

            INTEGER(iwp) ::  i   !< running index x
            INTEGER(iwp) ::  j   !< running index y
            INTEGER(iwp) ::  m   !< running index surface type

            REAL(wp), DIMENSION(1:surf_def_h(0)%ns) ::  def_1d !< 1D surface flux, default surfaces
            REAL(wp), DIMENSION(1:surf_lsm_h%ns)    ::  lsm_1d !< 1D surface flux, natural surfaces
            REAL(wp), DIMENSION(1:surf_usm_h%ns)    ::  usm_1d !< 1D surface flux, urban surfaces
!
!--         Transfer surface flux at default surfaces to 2D grid. Transfer no
!--         ghost-grid points since total_2d is a global array.
            DO  m = 1, SIZE(def_1d)
               i = surf_def_h(0)%i(m)
               j = surf_def_h(0)%j(m)

               IF ( i >= nxl  .AND.  i <= nxr  .AND.                           &
                    j >= nys  .AND.  j <= nyn )  THEN
                  total_2d(j,i) = def_1d(m)
               ENDIF
            ENDDO
!
!--         Transfer surface flux at natural surfaces to 2D grid 
            IF ( land_surface )  THEN 
               DO  m = 1, SIZE(lsm_1d)
                  i = surf_lsm_h%i(m)
                  j = surf_lsm_h%j(m)

                  IF ( i >= nxl  .AND.  i <= nxr  .AND.                        &
                       j >= nys  .AND.  j <= nyn )  THEN
                     total_2d(j,i) = lsm_1d(m)
                  ENDIF
               ENDDO
            ENDIF
!
!--         Transfer surface flux at natural surfaces to 2D grid 
            IF ( urban_surface )  THEN 
               DO  m = 1, SIZE(usm_1d)
                  i = surf_usm_h%i(m)
                  j = surf_usm_h%j(m)

                  IF ( i >= nxl  .AND.  i <= nxr  .AND.                        &
                       j >= nys  .AND.  j <= nyn )  THEN
                     total_2d(j,i) = usm_1d(m)
                  ENDIF
               ENDDO
            ENDIF

        END SUBROUTINE transfer_1D_to_2D_unequal

!       Description:
!------------------------------------------------------------------------------!
!>      Data transfer from 2D dummy array from unequal grids to 1D surface-data 
!>      type.
!------------------------------------------------------------------------------!
        SUBROUTINE transfer_2D_to_1D_unequal( def_1d )

           IMPLICIT NONE

            INTEGER(iwp) ::  i   !< running index x
            INTEGER(iwp) ::  j   !< running index y
            INTEGER(iwp) ::  m   !< running index surface type

            REAL(wp), DIMENSION(1:surf_def_h(2)%ns) ::  def_1d !< 1D surface flux, default surfaces
!
!--         Transfer 2D surface flux to default surfaces data type. Transfer no
!--         ghost-grid points since total_2d is a global array.
            DO  m = 1, SIZE(def_1d)
               i = surf_def_h(2)%i(m)
               j = surf_def_h(2)%j(m)

               IF ( i >= nxl  .AND.  i <= nxr  .AND.                           &
                    j >= nys  .AND.  j <= nyn )  THEN
                  def_1d(m) = total_2d_o(j,i)
               ENDIF
            ENDDO


        END SUBROUTINE transfer_2D_to_1D_unequal

#endif
  END SUBROUTINE surface_coupler



!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
#if defined( __parallel )

  SUBROUTINE interpolate_to_atmos( tag )

    USE arrays_3d,                                                             &
        ONLY:  total_2d_a, total_2d_o

    USE indices,                                                               &
        ONLY:  nbgp, nx, nx_a, nx_o, ny, ny_a, ny_o

    USE kinds

    USE pegrid

    IMPLICIT NONE

    INTEGER(iwp) ::  dnx  !< 
    INTEGER(iwp) ::  dnx2 !< 
    INTEGER(iwp) ::  dny  !< 
    INTEGER(iwp) ::  dny2 !< 
    INTEGER(iwp) ::  i    !< 
    INTEGER(iwp) ::  ii   !< 
    INTEGER(iwp) ::  j    !< 
    INTEGER(iwp) ::  jj   !< 

    INTEGER(iwp), intent(in) ::  tag !< 

    CALL MPI_BARRIER( comm2d, ierr )

    IF ( myid == 0 )  THEN
!
!--    Cyclic boundary conditions for the total 2D-grid
       total_2d_o(-nbgp:-1,:) = total_2d_o(ny+1-nbgp:ny,:)
       total_2d_o(:,-nbgp:-1) = total_2d_o(:,nx+1-nbgp:nx)

       total_2d_o(ny+1:ny+nbgp,:) = total_2d_o(0:nbgp-1,:)
       total_2d_o(:,nx+1:nx+nbgp) = total_2d_o(:,0:nbgp-1)

!
!--    Number of gridpoints of the fine grid within one mesh of the coarse grid
       dnx = (nx_o+1) / (nx_a+1) 
       dny = (ny_o+1) / (ny_a+1) 

!
!--    Distance for interpolation around coarse grid points within the fine
!--    grid (note: 2*dnx2 must not be equal with dnx)
       dnx2 = 2 * ( dnx / 2 )
       dny2 = 2 * ( dny / 2 )

       total_2d_a = 0.0_wp
!
!--    Interpolation from ocean-grid-layer to atmosphere-grid-layer
       DO  j = 0, ny_a
          DO  i = 0, nx_a 
             DO  jj = 0, dny2
                DO  ii = 0, dnx2
                   total_2d_a(j,i) = total_2d_a(j,i) &
                                     + total_2d_o(j*dny+jj,i*dnx+ii)
                ENDDO
             ENDDO
             total_2d_a(j,i) = total_2d_a(j,i) / ( ( dnx2 + 1 ) * ( dny2 + 1 ) )
          ENDDO
       ENDDO
! 
!--    Cyclic boundary conditions for atmosphere grid
       total_2d_a(-nbgp:-1,:) = total_2d_a(ny_a+1-nbgp:ny_a,:)
       total_2d_a(:,-nbgp:-1) = total_2d_a(:,nx_a+1-nbgp:nx_a)
       
       total_2d_a(ny_a+1:ny_a+nbgp,:) = total_2d_a(0:nbgp-1,:)
       total_2d_a(:,nx_a+1:nx_a+nbgp) = total_2d_a(:,0:nbgp-1)
!
!--    Transfer of the atmosphere-grid-layer to the atmosphere
       CALL MPI_SEND( total_2d_a(-nbgp,-nbgp), ngp_a, MPI_REAL, target_id, &
                      tag, comm_inter, ierr )

    ENDIF

    CALL MPI_BARRIER( comm2d, ierr )

  END SUBROUTINE interpolate_to_atmos

#endif


!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
#if defined( __parallel )

  SUBROUTINE interpolate_to_ocean( tag )

    USE arrays_3d,                                                             &
        ONLY:  total_2d_a, total_2d_o

    USE indices,                                                               &
        ONLY:  nbgp, nx, nx_a, nx_o, ny, ny_a, ny_o

    USE kinds

    USE pegrid

    IMPLICIT NONE

    INTEGER(iwp)             ::  dnx !< 
    INTEGER(iwp)             ::  dny !< 
    INTEGER(iwp)             ::  i   !< 
    INTEGER(iwp)             ::  ii  !< 
    INTEGER(iwp)             ::  j   !< 
    INTEGER(iwp)             ::  jj  !< 
    INTEGER(iwp), intent(in) ::  tag !< 

    REAL(wp)                 ::  fl  !< 
    REAL(wp)                 ::  fr  !< 
    REAL(wp)                 ::  myl !< 
    REAL(wp)                 ::  myr !< 

    CALL MPI_BARRIER( comm2d, ierr )

    IF ( myid == 0 )  THEN   

!
!--    Number of gridpoints of the fine grid within one mesh of the coarse grid
       dnx = ( nx_o + 1 ) / ( nx_a + 1 ) 
       dny = ( ny_o + 1 ) / ( ny_a + 1 ) 

! 
!--    Cyclic boundary conditions for atmosphere grid
       total_2d_a(-nbgp:-1,:) = total_2d_a(ny+1-nbgp:ny,:)
       total_2d_a(:,-nbgp:-1) = total_2d_a(:,nx+1-nbgp:nx)
       
       total_2d_a(ny+1:ny+nbgp,:) = total_2d_a(0:nbgp-1,:)
       total_2d_a(:,nx+1:nx+nbgp) = total_2d_a(:,0:nbgp-1)
!
!--    Bilinear Interpolation from atmosphere grid-layer to ocean grid-layer
       DO  j = 0, ny
          DO  i = 0, nx
             myl = ( total_2d_a(j+1,i)   - total_2d_a(j,i)   ) / dny
             myr = ( total_2d_a(j+1,i+1) - total_2d_a(j,i+1) ) / dny
             DO  jj = 0, dny-1
                fl = myl*jj + total_2d_a(j,i)  
                fr = myr*jj + total_2d_a(j,i+1)  
                DO  ii = 0, dnx-1
                   total_2d_o(j*dny+jj,i*dnx+ii) = ( fr - fl ) / dnx * ii + fl
                ENDDO
             ENDDO
          ENDDO
       ENDDO
! 
!--    Cyclic boundary conditions for ocean grid
       total_2d_o(-nbgp:-1,:) = total_2d_o(ny_o+1-nbgp:ny_o,:)
       total_2d_o(:,-nbgp:-1) = total_2d_o(:,nx_o+1-nbgp:nx_o)

       total_2d_o(ny_o+1:ny_o+nbgp,:) = total_2d_o(0:nbgp-1,:)
       total_2d_o(:,nx_o+1:nx_o+nbgp) = total_2d_o(:,0:nbgp-1)

       CALL MPI_SEND( total_2d_o(-nbgp,-nbgp), ngp_o, MPI_REAL, &
                      target_id, tag, comm_inter, ierr )

    ENDIF

    CALL MPI_BARRIER( comm2d, ierr )  

  END SUBROUTINE interpolate_to_ocean

#endif
