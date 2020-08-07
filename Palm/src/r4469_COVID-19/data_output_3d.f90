!> @file data_output_3d.f90
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
! $Id: data_output_3d.f90 4457 2020-03-11 14:20:43Z raasch $
! use statement for exchange horiz added
! 
! 4444 2020-03-05 15:59:50Z raasch
! bugfix: cpp-directives for serial mode added
! 
! 4360 2020-01-07 11:25:50Z suehring
! Introduction of wall_flags_total_0, which currently sets bits based on static
! topography information used in wall_flags_static_0
! 
! 4329 2019-12-10 15:46:36Z motisi
! Renamed wall_flags_0 to wall_flags_static_0
! 
! 4182 2019-08-22 15:20:23Z scharf
! Corrected "Former revisions" section
! 
! 4162 2019-08-16 05:54:29Z raasch
! bugfix for r4155
! 
! 4155 2019-08-14 06:25:18Z raasch
! bugfix for 3d-output in serial mode (ghost points must not be written)
! 
! 4127 2019-07-30 14:47:10Z suehring
! Adjustment for top boundary index for plant-canopy model outputs 
! (merge from branch resler)
! 
! 4048 2019-06-21 21:00:21Z knoop
! Moved tcm_data_output_3d to module_interface
! 
! 4039 2019-06-18 10:32:41Z suehring
! modularize diagnostic output
! 
! 3994 2019-05-22 18:08:09Z suehring
! output of turbulence intensity added
! 
! 3987 2019-05-22 09:52:13Z kanani
! Introduce alternative switch for debug output during timestepping
! 
! 3885 2019-04-11 11:29:34Z kanani
! Changes related to global restructuring of location messages and introduction 
! of additional debug messages
! 
! 3814 2019-03-26 08:40:31Z pavelkrc
! unused variables removed
! 
! 3655 2019-01-07 16:51:22Z knoop
! Bugfix: use time_since_reference_point instead of simulated_time (relevant
! when using wall/soil spinup)
!
! Revision 1.1  1997/09/03 06:29:36  raasch
! Initial revision
!
!
! Description:
! ------------
!> Output of the 3D-arrays in netCDF and/or AVS format.
!------------------------------------------------------------------------------!
 SUBROUTINE data_output_3d( av )
 

    USE arrays_3d,                                                             &
        ONLY:  d_exner, e, p, pt, q, ql, ql_c, ql_v, s, tend, u, v, vpt, w

    USE averaging

    USE basic_constants_and_equations_mod,                                     &
        ONLY:  lv_d_cp

    USE bulk_cloud_model_mod,                                                  &
        ONLY:  bulk_cloud_model

    USE control_parameters,                                                    &
        ONLY:  debug_output_timestep, do3d, do3d_no, do3d_time_count,          &
               land_surface, message_string, ntdim_3d, nz_do3d, plant_canopy,  &
               psolver, time_since_reference_point, urban_surface,             &
               varnamelength

#if defined( __parallel )
    USE control_parameters,                                                    &
        ONLY:  io_blocks, io_group
#endif

    USE cpulog,                                                                &
        ONLY:  log_point, cpu_log

    USE exchange_horiz_mod,                                                    &
        ONLY:  exchange_horiz

    USE indices,                                                               &
        ONLY:  nbgp, nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg, nzb, nzt,     &
               wall_flags_total_0

#if ! defined( __parallel )
    USE indices,                                                               &
        ONLY:  nx, ny
#endif

    USE kinds

    USE land_surface_model_mod,                                                &
        ONLY: lsm_data_output_3d, nzb_soil, nzt_soil

    USE module_interface,                                                      &
        ONLY:  module_interface_data_output_3d

#if defined( __netcdf )
    USE NETCDF
#endif

    USE netcdf_interface,                                                      &
        ONLY:  fill_value, id_set_3d, id_var_do3d, id_var_time_3d, nc_stat,    &
               netcdf_data_format, netcdf_handle_error

    USE particle_attributes,                                                   &
        ONLY:  grid_particles, number_of_particles, particles,                 &
               particle_advection_start, prt_count

    USE pegrid

    USE plant_canopy_model_mod,                                                &
        ONLY:  pch_index

    USE radiation_model_mod,                                                   &
        ONLY:  nz_urban_b, nz_urban_t

    USE urban_surface_mod,                                                     &
        ONLY:  usm_data_output_3d


    IMPLICIT NONE

    INTEGER(iwp) ::  av        !< flag for (non-)average output
    INTEGER(iwp) ::  flag_nr   !< number of masking flag
    INTEGER(iwp) ::  i         !< loop index
    INTEGER(iwp) ::  ivar      !< variable index
    INTEGER(iwp) ::  j         !< loop index
    INTEGER(iwp) ::  k         !< loop index
    INTEGER(iwp) ::  n         !< loop index
    INTEGER(iwp) ::  nzb_do    !< vertical lower limit for data output
    INTEGER(iwp) ::  nzt_do    !< vertical upper limit for data output

    LOGICAL      ::  found     !< true if output variable was found
    LOGICAL      ::  resorted  !< true if variable is resorted

    REAL(wp)     ::  mean_r    !< mean particle radius
    REAL(wp)     ::  s_r2      !< sum( particle-radius**2 )
    REAL(wp)     ::  s_r3      !< sum( particle-radius**3 )

    REAL(sp), DIMENSION(:,:,:), ALLOCATABLE ::  local_pf  !< output array

    REAL(wp), DIMENSION(:,:,:), POINTER ::  to_be_resorted  !< pointer to array which shall be output

    CHARACTER (LEN=varnamelength) ::  trimvar  !< TRIM of output-variable string

!
!-- Return, if nothing to output
    IF ( do3d_no(av) == 0 )  RETURN

    IF ( debug_output_timestep )  CALL debug_message( 'data_output_3d', 'start' )

    CALL cpu_log (log_point(14),'data_output_3d','start')

!
!-- Open output file.
!-- For classic or 64bit netCDF output on more than one PE, each PE opens its
!-- own file and writes the data of its subdomain in binary format. After the
!-- run, these files are combined to one NetCDF file by combine_plot_fields.
!-- For netCDF4/HDF5 output, data is written in parallel into one file.
    IF ( netcdf_data_format < 5 )  THEN
#if defined( __parallel )
       CALL check_open( 30 )
#endif
       IF ( myid == 0 )  CALL check_open( 106+av*10 )
    ELSE
       CALL check_open( 106+av*10 )
    ENDIF

!
!-- For parallel netcdf output the time axis must be limited. Return, if this
!-- limit is exceeded. This could be the case, if the simulated time exceeds 
!-- the given end time by the length of the given output interval.
    IF ( netcdf_data_format > 4 )  THEN
       IF ( do3d_time_count(av) + 1 > ntdim_3d(av) )  THEN
          WRITE ( message_string, * ) 'Output of 3d data is not given at t=',               &
                                      time_since_reference_point, 's because the maximum ', & 
                                      'number of output time levels is ',                   &
                                      'exceeded.'
          CALL message( 'data_output_3d', 'PA0387', 0, 1, 0, 6, 0 )
          CALL cpu_log( log_point(14), 'data_output_3d', 'stop' )
          RETURN
       ENDIF
    ENDIF

!
!-- Update the netCDF time axis
!-- In case of parallel output, this is only done by PE0 to increase the
!-- performance.
#if defined( __netcdf )
    do3d_time_count(av) = do3d_time_count(av) + 1
    IF ( myid == 0 )  THEN
       nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_time_3d(av),           &
                               (/ time_since_reference_point /),            &
                               start = (/ do3d_time_count(av) /),           &
                               count = (/ 1 /) )
       CALL netcdf_handle_error( 'data_output_3d', 376 )
    ENDIF
#endif

!
!-- Loop over all variables to be written.
    ivar = 1

    DO  WHILE ( do3d(av,ivar)(1:1) /= ' ' )

!
!--    Initiate found flag and resorting flag
       found = .FALSE.
       resorted = .FALSE.
       trimvar = TRIM( do3d(av,ivar) )

!
!--    Temporary solution to account for data output within the new urban
!--    surface model (urban_surface_mod.f90), see also SELECT CASE ( trimvar ).
!--    Store the array chosen on the temporary array.
       nzb_do   = nzb
!
!--    Set top index for 3D output. Note in case of plant-canopy model
!--    these index is determined by pch_index. 
       IF ( plant_canopy  .AND.  trimvar(1:4) == 'pcm_' )  THEN
          nzt_do   = pch_index
       ELSE
          nzt_do   = nz_do3d
       ENDIF

!
!--    Allocate a temporary array with the desired output dimensions.
       ALLOCATE( local_pf(nxl:nxr,nys:nyn,nzb_do:nzt_do) )
!
!--    Before each output, set array local_pf to fill value
       local_pf = fill_value
!
!--    Set masking flag for topography for not resorted arrays
       flag_nr = 0

       SELECT CASE ( trimvar )

          CASE ( 'e' )
             IF ( av == 0 )  THEN
                to_be_resorted => e
             ELSE
                IF ( .NOT. ALLOCATED( e_av ) ) THEN
                   ALLOCATE( e_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   e_av = REAL( fill_value, KIND = wp )
                ENDIF
                to_be_resorted => e_av
             ENDIF

          CASE ( 'thetal' )
             IF ( av == 0 )  THEN
                to_be_resorted => pt
             ELSE
                IF ( .NOT. ALLOCATED( lpt_av ) ) THEN
                   ALLOCATE( lpt_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   lpt_av = REAL( fill_value, KIND = wp )
                ENDIF
                to_be_resorted => lpt_av
             ENDIF

          CASE ( 'p' )
             IF ( av == 0 )  THEN
                IF ( psolver /= 'sor' )  CALL exchange_horiz( p, nbgp )
                to_be_resorted => p
             ELSE
                IF ( .NOT. ALLOCATED( p_av ) ) THEN
                   ALLOCATE( p_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   p_av = REAL( fill_value, KIND = wp )
                ENDIF
                IF ( psolver /= 'sor' )  CALL exchange_horiz( p_av, nbgp )
                to_be_resorted => p_av
             ENDIF

          CASE ( 'pc' )  ! particle concentration (requires ghostpoint exchange)
             IF ( av == 0 )  THEN
                IF ( time_since_reference_point >= particle_advection_start )  THEN
                   tend = prt_count
                ELSE
                   tend = 0.0_wp
                ENDIF
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         local_pf(i,j,k) = tend(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
                resorted = .TRUE.
             ELSE
                IF ( .NOT. ALLOCATED( pc_av ) ) THEN
                   ALLOCATE( pc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   pc_av = REAL( fill_value, KIND = wp )
                ENDIF
                to_be_resorted => pc_av
             ENDIF

          CASE ( 'pr' )  ! mean particle radius (effective radius)
             IF ( av == 0 )  THEN
                IF ( time_since_reference_point >= particle_advection_start )  THEN
                   DO  i = nxl, nxr
                      DO  j = nys, nyn
                         DO  k = nzb_do, nzt_do
                            number_of_particles = prt_count(k,j,i)
                            IF (number_of_particles <= 0)  CYCLE
                            particles => grid_particles(k,j,i)%particles(1:number_of_particles)
                            s_r2 = 0.0_wp
                            s_r3 = 0.0_wp
                            DO  n = 1, number_of_particles
                               IF ( particles(n)%particle_mask )  THEN
                                  s_r2 = s_r2 + particles(n)%radius**2 * &
                                         particles(n)%weight_factor
                                  s_r3 = s_r3 + particles(n)%radius**3 * &
                                         particles(n)%weight_factor
                               ENDIF
                            ENDDO
                            IF ( s_r2 > 0.0_wp )  THEN
                               mean_r = s_r3 / s_r2
                            ELSE
                               mean_r = 0.0_wp
                            ENDIF
                            tend(k,j,i) = mean_r
                         ENDDO
                      ENDDO
                   ENDDO
                ELSE
                   tend = 0.0_wp
                ENDIF
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         local_pf(i,j,k) = tend(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
                resorted = .TRUE.
             ELSE
                IF ( .NOT. ALLOCATED( pr_av ) ) THEN
                   ALLOCATE( pr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   pr_av = REAL( fill_value, KIND = wp )
                ENDIF
                to_be_resorted => pr_av
             ENDIF

          CASE ( 'theta' )
             IF ( av == 0 )  THEN
                IF ( .NOT. bulk_cloud_model ) THEN
                   to_be_resorted => pt
                ELSE
                   DO  i = nxl, nxr
                      DO  j = nys, nyn
                         DO  k = nzb_do, nzt_do
                            local_pf(i,j,k) = pt(k,j,i) + lv_d_cp *            &
                                                          d_exner(k) *         &
                                                          ql(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
                   resorted = .TRUE.
                ENDIF
             ELSE
                IF ( .NOT. ALLOCATED( pt_av ) ) THEN
                   ALLOCATE( pt_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   pt_av = REAL( fill_value, KIND = wp )
                ENDIF
                to_be_resorted => pt_av
             ENDIF

          CASE ( 'q' )
             IF ( av == 0 )  THEN
                to_be_resorted => q
             ELSE
                IF ( .NOT. ALLOCATED( q_av ) ) THEN
                   ALLOCATE( q_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   q_av = REAL( fill_value, KIND = wp )
                ENDIF
                to_be_resorted => q_av
             ENDIF

          CASE ( 'ql' )
             IF ( av == 0 )  THEN
                to_be_resorted => ql
             ELSE
                IF ( .NOT. ALLOCATED( ql_av ) ) THEN
                   ALLOCATE( ql_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ql_av = REAL( fill_value, KIND = wp )
                ENDIF
                to_be_resorted => ql_av
             ENDIF

          CASE ( 'ql_c' )
             IF ( av == 0 )  THEN
                to_be_resorted => ql_c
             ELSE
                IF ( .NOT. ALLOCATED( ql_c_av ) ) THEN
                   ALLOCATE( ql_c_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ql_c_av = REAL( fill_value, KIND = wp )
                ENDIF
                to_be_resorted => ql_c_av
             ENDIF

          CASE ( 'ql_v' )
             IF ( av == 0 )  THEN
                to_be_resorted => ql_v
             ELSE
                IF ( .NOT. ALLOCATED( ql_v_av ) ) THEN
                   ALLOCATE( ql_v_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ql_v_av = REAL( fill_value, KIND = wp )
                ENDIF
                to_be_resorted => ql_v_av
             ENDIF

          CASE ( 'ql_vp' )
             IF ( av == 0 )  THEN
                IF ( time_since_reference_point >= particle_advection_start )  THEN
                   DO  i = nxl, nxr
                      DO  j = nys, nyn
                         DO  k = nzb_do, nzt_do
                            number_of_particles = prt_count(k,j,i)
                            IF (number_of_particles <= 0)  CYCLE
                            particles => grid_particles(k,j,i)%particles(1:number_of_particles)
                            DO  n = 1, number_of_particles
                               IF ( particles(n)%particle_mask )  THEN
                                  tend(k,j,i) =  tend(k,j,i) +                 &
                                                 particles(n)%weight_factor /  &
                                                 prt_count(k,j,i)
                               ENDIF
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDDO
                ELSE
                   tend = 0.0_wp
                ENDIF
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         local_pf(i,j,k) = tend(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
                resorted = .TRUE.
             ELSE
                IF ( .NOT. ALLOCATED( ql_vp_av ) ) THEN
                   ALLOCATE( ql_vp_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ql_vp_av = REAL( fill_value, KIND = wp )
                ENDIF
                to_be_resorted => ql_vp_av
             ENDIF

          CASE ( 'qv' )
             IF ( av == 0 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         local_pf(i,j,k) = q(k,j,i) - ql(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
                resorted = .TRUE.
             ELSE
                IF ( .NOT. ALLOCATED( qv_av ) ) THEN
                   ALLOCATE( qv_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   qv_av = REAL( fill_value, KIND = wp )
                ENDIF
                to_be_resorted => qv_av
             ENDIF

          CASE ( 's' )
             IF ( av == 0 )  THEN
                to_be_resorted => s
             ELSE
                IF ( .NOT. ALLOCATED( s_av ) ) THEN
                   ALLOCATE( s_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   s_av = REAL( fill_value, KIND = wp )
                ENDIF
                to_be_resorted => s_av
             ENDIF

          CASE ( 'u' )
             flag_nr = 1
             IF ( av == 0 )  THEN
                to_be_resorted => u
             ELSE
                IF ( .NOT. ALLOCATED( u_av ) ) THEN
                   ALLOCATE( u_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   u_av = REAL( fill_value, KIND = wp )
                ENDIF
                to_be_resorted => u_av
             ENDIF

          CASE ( 'v' )
             flag_nr = 2
             IF ( av == 0 )  THEN
                to_be_resorted => v
             ELSE
                IF ( .NOT. ALLOCATED( v_av ) ) THEN
                   ALLOCATE( v_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   v_av = REAL( fill_value, KIND = wp )
                ENDIF
                to_be_resorted => v_av
             ENDIF

          CASE ( 'thetav' )
             IF ( av == 0 )  THEN
                to_be_resorted => vpt
             ELSE
                IF ( .NOT. ALLOCATED( vpt_av ) ) THEN
                   ALLOCATE( vpt_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   vpt_av = REAL( fill_value, KIND = wp )
                ENDIF
                to_be_resorted => vpt_av
             ENDIF

          CASE ( 'w' )
             flag_nr = 3
             IF ( av == 0 )  THEN
                to_be_resorted => w
             ELSE
                IF ( .NOT. ALLOCATED( w_av ) ) THEN
                   ALLOCATE( w_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   w_av = REAL( fill_value, KIND = wp )
                ENDIF
                to_be_resorted => w_av
             ENDIF

          CASE DEFAULT
!
!--          Quantities of other modules
             IF ( .NOT. found )  THEN
                CALL module_interface_data_output_3d(                          &
                        av, trimvar, found, local_pf,                          &
                        fill_value, resorted, nzb_do, nzt_do                   &
                     )
             ENDIF

!
!--          Temporary workaround: ToDo: refactor local_pf allocation
             IF ( .NOT. found  .AND.  urban_surface  .AND.  trimvar(1:4) == 'usm_' )  THEN
!
!--             For urban model quantities, it is required to re-allocate local_pf
                nzb_do = nz_urban_b
                nzt_do = nz_urban_t

                DEALLOCATE ( local_pf )
                ALLOCATE( local_pf(nxl:nxr,nys:nyn,nzb_do:nzt_do) )
                local_pf = fill_value

                CALL usm_data_output_3d( av, trimvar, found, local_pf,         &
                                         nzb_do, nzt_do )
                resorted = .TRUE.

!
!--             If no soil model variable was found, re-allocate local_pf
                IF ( .NOT. found )  THEN
                   nzb_do = nzb
                   nzt_do = nz_do3d

                   DEALLOCATE ( local_pf )
                   ALLOCATE( local_pf(nxl:nxr,nys:nyn,nzb_do:nzt_do) )
                ENDIF

             ENDIF

!
!--          Temporary workaround: ToDo: refactor local_pf allocation
             IF ( .NOT. found  .AND.  land_surface )  THEN
!
!--             For soil model quantities, it is required to re-allocate local_pf
                nzb_do = nzb_soil
                nzt_do = nzt_soil

                DEALLOCATE ( local_pf )
                ALLOCATE( local_pf(nxl:nxr,nys:nyn,nzb_do:nzt_do) )
                local_pf = fill_value

                CALL lsm_data_output_3d( av, trimvar, found, local_pf )
                resorted = .TRUE.

!
!--             If no soil model variable was found, re-allocate local_pf
                IF ( .NOT. found )  THEN
                   nzb_do = nzb
                   nzt_do = nz_do3d

                   DEALLOCATE ( local_pf )
                   ALLOCATE( local_pf(nxl:nxr,nys:nyn,nzb_do:nzt_do) )
                ENDIF

             ENDIF

             IF ( .NOT. found )  THEN
                message_string =  'no output available for: ' //               &
                                  TRIM( do3d(av,ivar) )
                CALL message( 'data_output_3d', 'PA0182', 0, 0, 0, 6, 0 )
             ENDIF

       END SELECT

!
!--    Resort the array to be output, if not done above
       IF ( .NOT. resorted )  THEN
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb_do, nzt_do
                   local_pf(i,j,k) = MERGE(                                    &
                                      to_be_resorted(k,j,i),                   &
                                      REAL( fill_value, KIND = wp ),           &
                                      BTEST( wall_flags_total_0(k,j,i), flag_nr ) )
                ENDDO
             ENDDO
          ENDDO
       ENDIF

!
!--    Output of the 3D-array
#if defined( __parallel )
       IF ( netcdf_data_format < 5 )  THEN
!
!--       Non-parallel netCDF output. Data is output in parallel in
!--       FORTRAN binary format here, and later collected into one file by
!--       combine_plot_fields
          IF ( myid == 0 )  THEN
             WRITE ( 30 )  time_since_reference_point,                   &
                           do3d_time_count(av), av
          ENDIF
          DO  i = 0, io_blocks-1
             IF ( i == io_group )  THEN
                WRITE ( 30 )  nxl, nxr, nys, nyn, nzb_do, nzt_do
                WRITE ( 30 )  local_pf(:,:,nzb_do:nzt_do)
             ENDIF

             CALL MPI_BARRIER( comm2d, ierr )

          ENDDO

       ELSE
#if defined( __netcdf )
!
!--       Parallel output in netCDF4/HDF5 format.
!          IF ( nxr == nx  .AND.  nyn /= ny )  THEN
!             nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_do3d(av,ivar),  &
!                               local_pf(nxl:nxr+1,nys:nyn,nzb_do:nzt_do),    &
!                start = (/ nxl+1, nys+1, nzb_do+1, do3d_time_count(av) /),  &
!                count = (/ nxr-nxl+2, nyn-nys+1, nzt_do-nzb_do+1, 1 /) )
!          ELSEIF ( nxr /= nx  .AND.  nyn == ny )  THEN
!             nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_do3d(av,ivar),  &
!                               local_pf(nxl:nxr,nys:nyn+1,nzb_do:nzt_do),    &
!                start = (/ nxl+1, nys+1, nzb_do+1, do3d_time_count(av) /),  &
!                count = (/ nxr-nxl+1, nyn-nys+2, nzt_do-nzb_do+1, 1 /) )
!          ELSEIF ( nxr == nx  .AND.  nyn == ny )  THEN
!             nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_do3d(av,ivar),  &
!                             local_pf(nxl:nxr+1,nys:nyn+1,nzb_do:nzt_do  ),  &
!                start = (/ nxl+1, nys+1, nzb_do+1, do3d_time_count(av) /),  &
!                count = (/ nxr-nxl+2, nyn-nys+2, nzt_do-nzb_do+1, 1 /) )
!          ELSE
             nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_do3d(av,ivar),  &
                                 local_pf(nxl:nxr,nys:nyn,nzb_do:nzt_do),    &
                start = (/ nxl+1, nys+1, nzb_do+1, do3d_time_count(av) /),  &
                count = (/ nxr-nxl+1, nyn-nys+1, nzt_do-nzb_do+1, 1 /) )
!          ENDIF
          CALL netcdf_handle_error( 'data_output_3d', 386 )
#endif
       ENDIF
#else
#if defined( __netcdf )
       nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_do3d(av,ivar),        &
                         local_pf(nxl:nxr,nys:nyn,nzb_do:nzt_do),          &
                         start = (/ 1, 1, 1, do3d_time_count(av) /),       &
                         count = (/ nx+1, ny+1, nzt_do-nzb_do+1, 1 /) )
       CALL netcdf_handle_error( 'data_output_3d', 446 )
#endif
#endif

       ivar = ivar + 1

!
!--    Deallocate temporary array
       DEALLOCATE ( local_pf )

    ENDDO

    CALL cpu_log( log_point(14), 'data_output_3d', 'stop' )

    IF ( debug_output_timestep )  CALL debug_message( 'data_output_3d', 'end' )


 END SUBROUTINE data_output_3d
