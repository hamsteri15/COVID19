!> @file surface_data_output_mod.f90
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
! $Id: surface_data_output_mod.f90 4444 2020-03-05 15:59:50Z raasch $
! bugfix: cpp-directives for serial mode added
! 
! 4360 2020-01-07 11:25:50Z suehring
! Fix wrongly declared nc_stat variable in surface_data_output_mod
! 
! 4205 2019-08-30 13:25:00Z suehring
! - Correct x,y-coordinates of vertical surfaces in netcdf output
! - Change definition of azimuth angle, reference is north 0 degree
! - zenith angle is always defined, also for vertical surfaces where it is 
!   90 degree, while azimuth angle is only defined for vertical surfaces, not
!   for horizontal ones
! 
! 4182 2019-08-22 15:20:23Z scharf
! Corrected "Former revisions" section
! 
! 4129 2019-07-31 12:56:07Z gronemeier
! - bugfix: corrected loop over horizontal default surfaces
! - change default setting of to_vtk and to_netcdf
!
! 4029 2019-06-14 14:04:35Z raasch
! netcdf variable NF90_NOFILL is used as argument instead of "1" in call to NF90_DEF_VAR_FILL
!
! 3881 2019-04-10 09:31:22Z suehring
! Check for zero output timestep (not allowed in parallel NetCDF output mode)
!
! 3817 2019-03-26 13:53:57Z suehring
! Correct output coordinates of vertical surface elements
!
! 3766 2019-02-26 16:23:41Z raasch
! bugfix in surface_data_output_rrd_local (variable k removed)
!
! 3762 2019-02-25 16:54:16Z suehring
! Remove unused variables and add preprocessor directives for variables that
! are used only when netcdf4 is defined
!
! 3745 2019-02-15 18:57:56Z suehring
! Output of waste_heat and innermost wall flux from indoor model
!
! 3744 2019-02-15 18:38:58Z suehring
! Add azimuth and zenith to output file; set long-name attributes;
! clean-up coding layout
!
! 3735 2019-02-12 09:52:40Z suehring
! - Split initialization into initialization of arrays and further initialization
!   in order to enable reading of restart data.
! - Consider restarts in surface data averaging.
! - Correct error message numbers
!
! 3731 2019-02-11 13:06:27Z suehring
! Bugfix: add cpp options
!
! 3727 2019-02-08 14:52:10Z gronemeier
! Enable NetCDF output for surface data (suehring, gronemeier)
!
! 3691 2019-01-23 09:57:04Z suehring
! Add output of surface-parallel flow speed
!
! 3648 2019-01-02 16:35:46Z suehring
! Rename module and subroutines
! 3420 2018-10-24 17:30:08Z gronemeier
! Initial implementation from Klaus Ketelsen and Matthias Suehring
!
!
! Authors:
! --------
! @author Klaus Ketelsen, Matthias Suehring, Tobias Gronemeier
!
! Description:
! ------------
!> Generate output for surface data.
!>
!> @todo Create namelist file for post-processing tool.
!------------------------------------------------------------------------------!

MODULE surface_data_output_mod

   USE kinds

   USE arrays_3d,                                                              &
       ONLY:  zu, zw

   USE control_parameters,                                                     &
       ONLY:  coupling_char, data_output_during_spinup, end_time,              &
              message_string, run_description_header, simulated_time_at_begin, &
              spinup_time, surface_output

   USE grid_variables,                                                         &
       ONLY: dx,dy

   USE indices,                                                                &
       ONLY: nxl, nxr, nys, nyn, nzb, nzt

#if defined( __netcdf )
   USE NETCDF
#endif

   USE netcdf_data_input_mod,                                                  &
       ONLY:  init_model

   USE netcdf_interface,                                                       &
       ONLY:  nc_stat, netcdf_create_att, netcdf_create_dim,                   &
              netcdf_create_file, netcdf_create_global_atts,                   &
              netcdf_create_var, netcdf_data_format, netcdf_handle_error

   USE pegrid

   USE surface_mod,                                                            &
       ONLY:  surf_def_h, surf_def_v, surf_lsm_h, surf_lsm_v,                  &
              surf_usm_h, surf_usm_v

   IMPLICIT NONE

   TYPE surf_out                      !< data structure which contains all surfaces elements of all types on subdomain

      INTEGER(iwp) ::  ns             !< number of surface elements on subdomain
      INTEGER(iwp) ::  ns_total       !< total number of surface elements
      INTEGER(iwp) ::  npoints        !< number of points / vertices which define a surface element (on subdomain)
      INTEGER(iwp) ::  npoints_total  !< total number of points / vertices which define a surface element

      INTEGER(iwp), DIMENSION(:), ALLOCATABLE   ::  s        !< coordinate for NetCDF output, number of the surface element

      REAL(wp) ::  fillvalue = -9999.0_wp !< fillvalue for surface elements which are not defined

      REAL(wp), DIMENSION(:), ALLOCATABLE   ::  azimuth  !< azimuth orientation coordinate for NetCDF output
      REAL(wp), DIMENSION(:), ALLOCATABLE   ::  es_utm   !< E-UTM coordinate for NetCDF output
      REAL(wp), DIMENSION(:), ALLOCATABLE   ::  ns_utm   !< E-UTM coordinate for NetCDF output
      REAL(wp), DIMENSION(:), ALLOCATABLE   ::  xs       !< x-coordinate for NetCDF output
      REAL(wp), DIMENSION(:), ALLOCATABLE   ::  ys       !< y-coordinate for NetCDF output
      REAL(wp), DIMENSION(:), ALLOCATABLE   ::  zs       !< z-coordinate for NetCDF output
      REAL(wp), DIMENSION(:), ALLOCATABLE   ::  zenith   !< zenith orientation coordinate for NetCDF output
      REAL(wp), DIMENSION(:), ALLOCATABLE   ::  var_out  !< output variables
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  var_av   !< variables used for averaging
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  points   !< points  / vertices of a surface element
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  polygons !< polygon data of a surface element
   END TYPE surf_out

   CHARACTER(LEN=100), DIMENSION(300)       ::  data_output_surf = ' '  !< namelist variable which describes the output variables
   CHARACTER(LEN=100), DIMENSION(0:1,300)   ::  dosurf = ' '            !< internal variable which describes the output variables
                                                                        !<  and separates averaged from non-averaged output
   CHARACTER(LEN=100), DIMENSION(0:1,300)   ::  dosurf_unit = ' '       !< internal variable which holds the unit of the given output variable

   INTEGER(iwp) ::  average_count_surf = 0   !< number of ensemble members used for averaging
   INTEGER(iwp) ::  dosurf_no(0:1) = 0       !< number of surface output quantities
#if defined( __netcdf4_parallel )
   INTEGER(iwp) ::  oldmode                  !< save old set-fill-mode of netcdf file (not needed, but required for routine call)

   INTEGER(iwp), DIMENSION(0:1) ::  dosurf_time_count = 0 !< count of output time steps
   INTEGER(iwp), DIMENSION(0:1) ::  id_dim_s_surf         !< netcdf ID for dimension s
   INTEGER(iwp), DIMENSION(0:1) ::  id_dim_time_surf      !< netcdf ID for dimension time
   INTEGER(iwp), DIMENSION(0:1) ::  id_set_surf           !< netcdf ID for file
   INTEGER(iwp), DIMENSION(0:1) ::  id_var_azimuth_surf   !< netcdf ID for variable azimuth
   INTEGER(iwp), DIMENSION(0:1) ::  id_var_etum_surf      !< netcdf ID for variable Es_UTM
   INTEGER(iwp), DIMENSION(0:1) ::  id_var_nutm_surf      !< netcdf ID for variable Ns_UTM
   INTEGER(iwp), DIMENSION(0:1) ::  id_var_time_surf      !< netcdf ID for variable time
   INTEGER(iwp), DIMENSION(0:1) ::  id_var_s_surf         !< netcdf ID for variable s
   INTEGER(iwp), DIMENSION(0:1) ::  id_var_xs_surf        !< netcdf ID for variable xs
   INTEGER(iwp), DIMENSION(0:1) ::  id_var_ys_surf        !< netcdf ID for variable ys
   INTEGER(iwp), DIMENSION(0:1) ::  id_var_zenith_surf    !< netcdf ID for variable zenith
   INTEGER(iwp), DIMENSION(0:1) ::  id_var_zs_surf        !< netcdf ID for variable zs
   INTEGER(iwp), DIMENSION(0:1) ::  ntdim_surf            !< number of output time steps

   INTEGER(iwp), DIMENSION(0:1,300) ::  id_var_dosurf     !< netcdf ID for output variables
#endif

   LOGICAL :: first_output(0:1) = .FALSE.                 !< true if first output was already called
   LOGICAL :: to_netcdf = .FALSE.                         !< flag indicating parallel NetCDF output
   LOGICAL :: to_vtk = .FALSE.                            !< flag indicating binary surface-data output that can be further
                                                          !< processed to VTK format

   REAL(wp) ::  averaging_interval_surf  = 9999999.9_wp   !< averaging interval
   REAL(wp) ::  dt_dosurf = 9999999.9_wp                  !< time interval for instantaneous data output
   REAL(wp) ::  dt_dosurf_av = 9999999.9_wp               !< time interval for averaged data output
   REAL(wp) ::  skip_time_dosurf = 0.0_wp                 !< skip time for instantaneous data output
   REAL(wp) ::  skip_time_dosurf_av = 0.0_wp              !< skip time for averaged data output
   REAL(wp) ::  time_dosurf = 0.0_wp                      !< internal counter variable to check for instantaneous data output
   REAL(wp) ::  time_dosurf_av = 0.0_wp                   !< internal counter variable to check for averaged data output

   TYPE(surf_out) ::  surfaces      !< variable which contains all required output information

   SAVE

   PRIVATE

   INTERFACE  surface_data_output
      MODULE PROCEDURE surface_data_output
   END INTERFACE  surface_data_output

   INTERFACE  surface_data_output_averaging
      MODULE PROCEDURE surface_data_output_averaging
   END INTERFACE  surface_data_output_averaging

   INTERFACE  surface_data_output_check_parameters
      MODULE PROCEDURE surface_data_output_check_parameters
   END INTERFACE  surface_data_output_check_parameters

   INTERFACE  surface_data_output_init
      MODULE PROCEDURE surface_data_output_init
   END INTERFACE  surface_data_output_init

   INTERFACE  surface_data_output_init_arrays
      MODULE PROCEDURE surface_data_output_init_arrays
   END INTERFACE  surface_data_output_init_arrays

   INTERFACE  surface_data_output_last_action
      MODULE PROCEDURE surface_data_output_last_action
   END INTERFACE  surface_data_output_last_action

   INTERFACE  surface_data_output_parin
      MODULE PROCEDURE surface_data_output_parin
   END INTERFACE  surface_data_output_parin

   INTERFACE  surface_data_output_rrd_global
      MODULE PROCEDURE surface_data_output_rrd_global
   END INTERFACE  surface_data_output_rrd_global

   INTERFACE  surface_data_output_rrd_local
      MODULE PROCEDURE surface_data_output_rrd_local
   END INTERFACE  surface_data_output_rrd_local

   INTERFACE  surface_data_output_wrd_global
      MODULE PROCEDURE surface_data_output_wrd_global
   END INTERFACE  surface_data_output_wrd_global

   INTERFACE  surface_data_output_wrd_local
      MODULE PROCEDURE surface_data_output_wrd_local
   END INTERFACE  surface_data_output_wrd_local

!
!--Public subroutines
   PUBLIC surface_data_output, surface_data_output_averaging,                  &
          surface_data_output_check_parameters, surface_data_output_init,      &
          surface_data_output_init_arrays, surface_data_output_last_action,    &
          surface_data_output_parin, surface_data_output_rrd_global,           &
          surface_data_output_rrd_local, surface_data_output_wrd_local,        &
          surface_data_output_wrd_global
!
!--Public variables
   PUBLIC average_count_surf, averaging_interval_surf, dt_dosurf, dt_dosurf_av,&
          skip_time_dosurf, skip_time_dosurf_av, time_dosurf, time_dosurf_av

 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine counts the number of surfaces on each core and allocates
!> arrays.
!------------------------------------------------------------------------------!
   SUBROUTINE surface_data_output_init_arrays

      IMPLICIT NONE

!
!--   Determine the number of surface elements on subdomain
      surfaces%ns = surf_def_h(0)%ns + surf_lsm_h%ns + surf_usm_h%ns           & !horizontal upward-facing
                  + surf_def_h(1)%ns                                           & !horizontal downard-facing
                  + surf_def_v(0)%ns + surf_lsm_v(0)%ns + surf_usm_v(0)%ns     & !northward-facing
                  + surf_def_v(1)%ns + surf_lsm_v(1)%ns + surf_usm_v(1)%ns     & !southward-facing
                  + surf_def_v(2)%ns + surf_lsm_v(2)%ns + surf_usm_v(2)%ns     & !westward-facing
                  + surf_def_v(3)%ns + surf_lsm_v(3)%ns + surf_usm_v(3)%ns       !eastward-facing
!
!--    Determine the total number of surfaces in the model domain
#if defined( __parallel )
       CALL MPI_ALLREDUCE( surfaces%ns, surfaces%ns_total, 1,                  &
                           MPI_INTEGER, MPI_SUM, comm2d, ierr )
#else
       surfaces%ns_total = surfaces%ns
#endif
!
!--   Allocate output variable and set to _FillValue attribute
      ALLOCATE ( surfaces%var_out(1:surfaces%ns) )
      surfaces%var_out = surfaces%fillvalue
!
!--   If there is an output of time average output variables, allocate the
!--   required array.
      IF ( dosurf_no(1) > 0 )  THEN
         ALLOCATE ( surfaces%var_av(1:surfaces%ns,1:dosurf_no(1)) )
         surfaces%var_av = 0.0_wp
      ENDIF

   END SUBROUTINE surface_data_output_init_arrays


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization surface-data output data structure: calculation of vertices
!> and polygon data for the surface elements, allocation of required arrays.
!------------------------------------------------------------------------------!
   SUBROUTINE surface_data_output_init

      IMPLICIT NONE

#if defined( __netcdf4_parallel )
      CHARACTER (LEN=100)  :: filename            !< name of output file
      CHARACTER (LEN=80)   ::  time_average_text  !< string written to file attribute time_avg
      CHARACTER (LEN=4000) ::  var_list           !< list of variables written to NetCDF file

      INTEGER(iwp) ::  av                !< flag for averaged (=1) and non-averaged (=0) data
#endif
      INTEGER(iwp) ::  i                 !< grid index in x-direction, also running variable for counting non-average data output
      INTEGER(iwp) ::  j                 !< grid index in y-direction, also running variable for counting average data output
      INTEGER(iwp) ::  k                 !< grid index in z-direction
      INTEGER(iwp) ::  l                 !< running index for surface-element orientation
      INTEGER(iwp) ::  m                 !< running index for surface elements
      INTEGER(iwp) ::  mm                !< local counting variable for surface elements
      INTEGER(iwp) ::  npg               !< counter variable for all surface elements ( or polygons )
      INTEGER(iwp) ::  point_index_count !< local counter variable for point index
      INTEGER(iwp) ::  start_count       !< local start counter for the surface index

      INTEGER(iwp), DIMENSION(0:numprocs-1) :: num_points_on_pe   !< array which contains the number of points on all mpi ranks
      INTEGER(iwp), DIMENSION(0:numprocs-1) :: num_surfaces_on_pe !< array which contains the number of surfaces on all mpi ranks
      INTEGER(iwp), ALLOCATABLE, DIMENSION(:,:,:) ::  point_index !< dummy array used to check where the reference points for surface polygons are located

      REAL(wp) ::  az    !< azimuth angle, indicated the vertical orientation of a surface element
      REAL(wp) ::  off_x !< grid offset in x-direction between the stored grid index and the actual wall
      REAL(wp) ::  off_y !< grid offset in y-direction between the stored grid index and the actual wall
#if defined( __netcdf4_parallel )
      REAL(wp), DIMENSION(:), ALLOCATABLE ::  netcdf_data_1d  !< dummy array to output 1D data into netcdf file
#endif

!
!--   If output to VTK format is enabled, initialize point and polygon data.
!--   In a first step, count the number of points which are defining
!--   the surfaces and the polygons.
      IF ( to_vtk )  THEN
         ALLOCATE( point_index(nzb:nzt+1,nys-1:nyn+1,nxl-1:nxr+1) )
         point_index = -1
!
!--      Horizontal default surfaces
         surfaces%npoints = 0
         DO  l = 0, 1
            DO  m = 1, surf_def_h(0)%ns
!
!--            Determine the indices of the respective grid cell inside the topography
               i = surf_def_h(0)%i(m) + surf_def_h(0)%ioff
               j = surf_def_h(0)%j(m) + surf_def_h(0)%joff
               k = surf_def_h(0)%k(m) + surf_def_h(0)%koff
!
!--            Check if the vertices that define the surface element are already
!--            defined, if not, increment the counter.
               IF ( point_index(k,j,i) < 0 )  THEN
                  surfaces%npoints   = surfaces%npoints + 1
                  point_index(k,j,i) = surfaces%npoints - 1
               ENDIF
               IF ( point_index(k,j,i+1) < 0 )  THEN
                  surfaces%npoints     = surfaces%npoints + 1
                  point_index(k,j,i+1) = surfaces%npoints - 1
               ENDIF
               IF ( point_index(k,j+1,i+1) < 0 )  THEN
                  surfaces%npoints       = surfaces%npoints + 1
                  point_index(k,j+1,i+1) = surfaces%npoints - 1
               ENDIF
               IF ( point_index(k,j+1,i) < 0 )  THEN
                  surfaces%npoints     = surfaces%npoints + 1
                  point_index(k,j+1,i) = surfaces%npoints - 1
               ENDIF
            ENDDO
         ENDDO
         DO  m = 1, surf_lsm_h%ns
            i = surf_lsm_h%i(m) + surf_lsm_h%ioff
            j = surf_lsm_h%j(m) + surf_lsm_h%joff
            k = surf_lsm_h%k(m) + surf_lsm_h%koff

            IF ( point_index(k,j,i) < 0 )  THEN
               surfaces%npoints   = surfaces%npoints + 1
               point_index(k,j,i) = surfaces%npoints - 1
            ENDIF
            IF ( point_index(k,j,i+1) < 0 )  THEN
               surfaces%npoints     = surfaces%npoints + 1
               point_index(k,j,i+1) = surfaces%npoints - 1
            ENDIF
            IF ( point_index(k,j+1,i+1) < 0 )  THEN
               surfaces%npoints       = surfaces%npoints + 1
               point_index(k,j+1,i+1) = surfaces%npoints - 1
            ENDIF
            IF ( point_index(k,j+1,i) < 0 )  THEN
               surfaces%npoints     = surfaces%npoints + 1
               point_index(k,j+1,i) = surfaces%npoints - 1
            ENDIF
         ENDDO
         DO  m = 1, surf_usm_h%ns
            i = surf_usm_h%i(m) + surf_usm_h%ioff
            j = surf_usm_h%j(m) + surf_usm_h%joff
            k = surf_usm_h%k(m) + surf_usm_h%koff

            IF ( point_index(k,j,i) < 0 )  THEN
               surfaces%npoints   = surfaces%npoints + 1
               point_index(k,j,i) = surfaces%npoints - 1
            ENDIF
            IF ( point_index(k,j,i+1) < 0 )  THEN
               surfaces%npoints     = surfaces%npoints + 1
               point_index(k,j,i+1) = surfaces%npoints - 1
            ENDIF
            IF ( point_index(k,j+1,i+1) < 0 )  THEN
               surfaces%npoints       = surfaces%npoints + 1
               point_index(k,j+1,i+1) = surfaces%npoints - 1
            ENDIF
            IF ( point_index(k,j+1,i) < 0 )  THEN
               surfaces%npoints     = surfaces%npoints + 1
               point_index(k,j+1,i) = surfaces%npoints - 1
            ENDIF
         ENDDO
!
!--      Vertical surfaces
         DO  l = 0, 3
            DO  m = 1, surf_def_v(l)%ns
!
!--            Determine the indices of the respective grid cell inside the
!--            topography. Please note, j-index for north-facing surfaces
!--            ( l==0 ) is identical to the reference j-index outside the grid
!--            box. Equivalent for east-facing surfaces and i-index.
               i = surf_def_v(l)%i(m) + MERGE( surf_def_v(l)%ioff, 0, l == 3 )
               j = surf_def_v(l)%j(m) + MERGE( surf_def_v(l)%joff, 0, l == 1 )
               k = surf_def_v(l)%k(m) + surf_def_v(l)%koff
!
!--            Lower left /front vertex
               IF ( point_index(k,j,i) < 0 )  THEN
                  surfaces%npoints   = surfaces%npoints + 1
                  point_index(k,j,i) = surfaces%npoints - 1
               ENDIF
!
!--            Upper / lower right index for north- and south-facing surfaces
               IF ( l == 0  .OR.  l == 1 )  THEN
                  IF ( point_index(k,j,i+1) < 0 )  THEN
                     surfaces%npoints     = surfaces%npoints + 1
                     point_index(k,j,i+1) = surfaces%npoints - 1
                  ENDIF
                  IF ( point_index(k+1,j,i+1) < 0 )  THEN
                     surfaces%npoints       = surfaces%npoints + 1
                     point_index(k+1,j,i+1) = surfaces%npoints - 1
                  ENDIF
!
!--            Upper / lower front index for east- and west-facing surfaces
               ELSEIF ( l == 2  .OR.  l == 3 )  THEN
                  IF ( point_index(k,j+1,i) < 0 )  THEN
                     surfaces%npoints     = surfaces%npoints + 1
                     point_index(k,j+1,i) = surfaces%npoints - 1
                  ENDIF
                  IF ( point_index(k+1,j+1,i) < 0 )  THEN
                     surfaces%npoints       = surfaces%npoints + 1
                     point_index(k+1,j+1,i) = surfaces%npoints - 1
                  ENDIF
               ENDIF
!
!--            Upper left / front vertex
               IF ( point_index(k+1,j,i) < 0 )  THEN
                  surfaces%npoints     = surfaces%npoints + 1
                  point_index(k+1,j,i) = surfaces%npoints - 1
               ENDIF
            ENDDO
            DO  m = 1, surf_lsm_v(l)%ns
               i = surf_lsm_v(l)%i(m) + MERGE( surf_lsm_v(l)%ioff, 0, l == 3 )
               j = surf_lsm_v(l)%j(m) + MERGE( surf_lsm_v(l)%joff, 0, l == 1 )
               k = surf_lsm_v(l)%k(m) + surf_lsm_v(l)%koff
!
!--            Lower left /front vertex
               IF ( point_index(k,j,i) < 0 )  THEN
                  surfaces%npoints   = surfaces%npoints + 1
                  point_index(k,j,i) = surfaces%npoints - 1
               ENDIF
!
!--            Upper / lower right index for north- and south-facing surfaces
               IF ( l == 0  .OR.  l == 1 )  THEN
                  IF ( point_index(k,j,i+1) < 0 )  THEN
                     surfaces%npoints     = surfaces%npoints + 1
                     point_index(k,j,i+1) = surfaces%npoints - 1
                  ENDIF
                  IF ( point_index(k+1,j,i+1) < 0 )  THEN
                     surfaces%npoints       = surfaces%npoints + 1
                     point_index(k+1,j,i+1) = surfaces%npoints - 1
                  ENDIF
!
!--            Upper / lower front index for east- and west-facing surfaces
               ELSEIF ( l == 2  .OR.  l == 3 )  THEN
                  IF ( point_index(k,j+1,i) < 0 )  THEN
                     surfaces%npoints     = surfaces%npoints + 1
                     point_index(k,j+1,i) = surfaces%npoints - 1
                  ENDIF
                  IF ( point_index(k+1,j+1,i) < 0 )  THEN
                     surfaces%npoints       = surfaces%npoints + 1
                     point_index(k+1,j+1,i) = surfaces%npoints - 1
                  ENDIF
               ENDIF
!
!--            Upper left / front vertex
               IF ( point_index(k+1,j,i) < 0 )  THEN
                  surfaces%npoints     = surfaces%npoints + 1
                  point_index(k+1,j,i) = surfaces%npoints - 1
               ENDIF
            ENDDO

            DO  m = 1, surf_usm_v(l)%ns
               i = surf_usm_v(l)%i(m) + MERGE( surf_usm_v(l)%ioff, 0, l == 3 )
               j = surf_usm_v(l)%j(m) + MERGE( surf_usm_v(l)%joff, 0, l == 1 )
               k = surf_usm_v(l)%k(m) + surf_usm_v(l)%koff
!
!--            Lower left /front vertex
               IF ( point_index(k,j,i) < 0 )  THEN
                  surfaces%npoints   = surfaces%npoints + 1
                  point_index(k,j,i) = surfaces%npoints - 1
               ENDIF
!
!--            Upper / lower right index for north- and south-facing surfaces
               IF ( l == 0  .OR.  l == 1 )  THEN
                  IF ( point_index(k,j,i+1) < 0 )  THEN
                     surfaces%npoints     = surfaces%npoints + 1
                     point_index(k,j,i+1) = surfaces%npoints - 1
                  ENDIF
                  IF ( point_index(k+1,j,i+1) < 0 )  THEN
                     surfaces%npoints       = surfaces%npoints + 1
                     point_index(k+1,j,i+1) = surfaces%npoints - 1
                  ENDIF
!
!--            Upper / lower front index for east- and west-facing surfaces
               ELSEIF ( l == 2  .OR.  l == 3 )  THEN
                  IF ( point_index(k,j+1,i) < 0 )  THEN
                     surfaces%npoints     = surfaces%npoints + 1
                     point_index(k,j+1,i) = surfaces%npoints - 1
                  ENDIF
                  IF ( point_index(k+1,j+1,i) < 0 )  THEN
                     surfaces%npoints       = surfaces%npoints + 1
                     point_index(k+1,j+1,i) = surfaces%npoints - 1
                  ENDIF
               ENDIF
!
!--            Upper left / front vertex
               IF ( point_index(k+1,j,i) < 0 )  THEN
                  surfaces%npoints     = surfaces%npoints + 1
                  point_index(k+1,j,i) = surfaces%npoints - 1
               ENDIF
            ENDDO

         ENDDO
!
!--      Allocate the number of points and polygons. Note, the number of
!--      polygons is identical to the number of surfaces elements, whereas the
!--      number of points (vertices), which define the polygons, can be larger.
         ALLOCATE( surfaces%points(3,1:surfaces%npoints) )
         ALLOCATE( surfaces%polygons(5,1:surfaces%ns)    )
!
!--      Note, PARAVIEW expects consecutively ordered points, in order to
!--      unambiguously identify surfaces.
!--      Hence, all PEs should know where they start counting, depending on the
!--      number of points on the other PE's with lower MPI rank.
#if defined( __parallel )
         CALL MPI_ALLGATHER( surfaces%npoints, 1, MPI_INTEGER,                 &
                             num_points_on_pe, 1, MPI_INTEGER, comm2d, ierr  )
#else
         num_points_on_pe = surfaces%npoints
#endif

!
!--      After the number of vertices is counted, repeat the loops and define
!--      the vertices. Start with the horizontal default surfaces.
!--      First, however, determine the offset where couting of points should be
!--      started, which is the sum of points of all PE's with lower MPI rank.
         i                 = 0
         point_index_count = 0
         DO WHILE ( i < myid  .AND.  i <= SIZE( num_points_on_pe ) )
            point_index_count = point_index_count + num_points_on_pe(i)
            i                 = i + 1
         ENDDO

         surfaces%npoints = 0
         point_index      = -1
         npg              = 0

         DO  l = 0, 1
            DO  m = 1, surf_def_h(l)%ns
!
!--            Determine the indices of the respective grid cell inside the
!--            topography.
               i = surf_def_h(l)%i(m) + surf_def_h(l)%ioff
               j = surf_def_h(l)%j(m) + surf_def_h(l)%joff
               k = surf_def_h(l)%k(m) + surf_def_h(l)%koff
!
!--            Check if the vertices that define the surface element are
!--            already defined, if not, increment the counter.
               IF ( point_index(k,j,i) < 0 )  THEN
                  surfaces%npoints   = surfaces%npoints + 1
                  point_index(k,j,i) = point_index_count
                  point_index_count  = point_index_count + 1
                  surfaces%points(1,surfaces%npoints) = ( i - 0.5_wp ) * dx
                  surfaces%points(2,surfaces%npoints) = ( j - 0.5_wp ) * dy
                  surfaces%points(3,surfaces%npoints) = zw(k)
               ENDIF
               IF ( point_index(k,j,i+1) < 0 )  THEN
                  surfaces%npoints     = surfaces%npoints + 1
                  point_index(k,j,i+1) = point_index_count
                  point_index_count    = point_index_count + 1
                  surfaces%points(1,surfaces%npoints) = ( i + 1 - 0.5_wp ) * dx
                  surfaces%points(2,surfaces%npoints) = ( j     - 0.5_wp ) * dy
                  surfaces%points(3,surfaces%npoints) = zw(k)
               ENDIF
               IF ( point_index(k,j+1,i+1) < 0 )  THEN
                  surfaces%npoints       = surfaces%npoints + 1
                  point_index(k,j+1,i+1) = point_index_count
                  point_index_count      = point_index_count + 1
                  surfaces%points(1,surfaces%npoints) = ( i + 1 - 0.5_wp ) * dx
                  surfaces%points(2,surfaces%npoints) = ( j + 1 - 0.5_wp ) * dy
                  surfaces%points(3,surfaces%npoints) = zw(k)
               ENDIF
               IF ( point_index(k,j+1,i) < 0 )  THEN
                  surfaces%npoints     = surfaces%npoints + 1
                  point_index(k,j+1,i) = point_index_count
                  point_index_count    = point_index_count + 1
                  surfaces%points(1,surfaces%npoints) = ( i     - 0.5_wp ) * dx
                  surfaces%points(2,surfaces%npoints) = ( j + 1 - 0.5_wp ) * dy
                  surfaces%points(3,surfaces%npoints) = zw(k)
               ENDIF

               npg                        = npg + 1
               surfaces%polygons(1,npg)   = 4
               surfaces%polygons(2,npg)   = point_index(k,j,i)
               surfaces%polygons(3,npg)   = point_index(k,j,i+1)
               surfaces%polygons(4,npg)   = point_index(k,j+1,i+1)
               surfaces%polygons(5,npg)   = point_index(k,j+1,i)
            ENDDO
         ENDDO
         DO  m = 1, surf_lsm_h%ns
            i = surf_lsm_h%i(m) + surf_lsm_h%ioff
            j = surf_lsm_h%j(m) + surf_lsm_h%joff
            k = surf_lsm_h%k(m) + surf_lsm_h%koff
            IF ( point_index(k,j,i) < 0 )  THEN
               surfaces%npoints   = surfaces%npoints + 1
               point_index(k,j,i) = point_index_count
               point_index_count  = point_index_count + 1
               surfaces%points(1,surfaces%npoints) = ( i - 0.5_wp ) * dx
               surfaces%points(2,surfaces%npoints) = ( j - 0.5_wp ) * dy
               surfaces%points(3,surfaces%npoints) = zw(k)
            ENDIF
            IF ( point_index(k,j,i+1) < 0 )  THEN
               surfaces%npoints     = surfaces%npoints + 1
               point_index(k,j,i+1) = point_index_count
               point_index_count    = point_index_count + 1
               surfaces%points(1,surfaces%npoints) = ( i + 1 - 0.5_wp ) * dx
               surfaces%points(2,surfaces%npoints) = ( j     - 0.5_wp ) * dy
               surfaces%points(3,surfaces%npoints) = zw(k)
            ENDIF
            IF ( point_index(k,j+1,i+1) < 0 )  THEN
               surfaces%npoints       = surfaces%npoints + 1
               point_index(k,j+1,i+1) = point_index_count
               point_index_count      = point_index_count + 1
               surfaces%points(1,surfaces%npoints) = ( i + 1 - 0.5_wp ) * dx
               surfaces%points(2,surfaces%npoints) = ( j + 1 - 0.5_wp ) * dy
               surfaces%points(3,surfaces%npoints) = zw(k)
            ENDIF
            IF ( point_index(k,j+1,i) < 0 )  THEN
               surfaces%npoints     = surfaces%npoints + 1
               point_index(k,j+1,i) = point_index_count
               point_index_count    = point_index_count + 1
               surfaces%points(1,surfaces%npoints) = ( i     - 0.5_wp ) * dx
               surfaces%points(2,surfaces%npoints) = ( j + 1 - 0.5_wp ) * dy
               surfaces%points(3,surfaces%npoints) = zw(k)
            ENDIF

            npg                        = npg + 1
            surfaces%polygons(1,npg)   = 4
            surfaces%polygons(2,npg)   = point_index(k,j,i)
            surfaces%polygons(3,npg)   = point_index(k,j,i+1)
            surfaces%polygons(4,npg)   = point_index(k,j+1,i+1)
            surfaces%polygons(5,npg)   = point_index(k,j+1,i)
         ENDDO

         DO  m = 1, surf_usm_h%ns
            i = surf_usm_h%i(m) + surf_usm_h%ioff
            j = surf_usm_h%j(m) + surf_usm_h%joff
            k = surf_usm_h%k(m) + surf_usm_h%koff

            IF ( point_index(k,j,i) < 0 )  THEN
               surfaces%npoints   = surfaces%npoints + 1
               point_index(k,j,i) = point_index_count
               point_index_count  = point_index_count + 1
               surfaces%points(1,surfaces%npoints) = ( i - 0.5_wp ) * dx
               surfaces%points(2,surfaces%npoints) = ( j - 0.5_wp ) * dy
               surfaces%points(3,surfaces%npoints) = zw(k)
            ENDIF
            IF ( point_index(k,j,i+1) < 0 )  THEN
               surfaces%npoints     = surfaces%npoints + 1
               point_index(k,j,i+1) = point_index_count
               point_index_count    = point_index_count + 1
               surfaces%points(1,surfaces%npoints) = ( i + 1 - 0.5_wp ) * dx
               surfaces%points(2,surfaces%npoints) = ( j     - 0.5_wp ) * dy
               surfaces%points(3,surfaces%npoints) = zw(k)
            ENDIF
            IF ( point_index(k,j+1,i+1) < 0 )  THEN
               surfaces%npoints       = surfaces%npoints + 1
               point_index(k,j+1,i+1) = point_index_count
               point_index_count      = point_index_count + 1
               surfaces%points(1,surfaces%npoints) = ( i + 1 - 0.5_wp ) * dx
               surfaces%points(2,surfaces%npoints) = ( j + 1 - 0.5_wp ) * dy
               surfaces%points(3,surfaces%npoints) = zw(k)
            ENDIF
            IF ( point_index(k,j+1,i) < 0 )  THEN
               surfaces%npoints     = surfaces%npoints + 1
               point_index(k,j+1,i) = point_index_count
               point_index_count    = point_index_count + 1
               surfaces%points(1,surfaces%npoints) = ( i     - 0.5_wp ) * dx
               surfaces%points(2,surfaces%npoints) = ( j + 1 - 0.5_wp ) * dy
               surfaces%points(3,surfaces%npoints) = zw(k)
            ENDIF

            npg                        = npg + 1
            surfaces%polygons(1,npg)   = 4
            surfaces%polygons(2,npg)   = point_index(k,j,i)
            surfaces%polygons(3,npg)   = point_index(k,j,i+1)
            surfaces%polygons(4,npg)   = point_index(k,j+1,i+1)
            surfaces%polygons(5,npg)   = point_index(k,j+1,i)
         ENDDO

         DO  l = 0, 3
            DO  m = 1, surf_def_v(l)%ns
!
!--            Determine the indices of the respective grid cell inside the
!--            topography.
!--            NOTE, j-index for north-facing surfaces ( l==0 ) is
!--            identical to the reference j-index outside the grid box.
!--            Equivalent for east-facing surfaces and i-index.
               i = surf_def_v(l)%i(m) + MERGE( surf_def_v(l)%ioff, 0, l == 3 )
               j = surf_def_v(l)%j(m) + MERGE( surf_def_v(l)%joff, 0, l == 1 )
               k = surf_def_v(l)%k(m) + surf_def_v(l)%koff
!
!--            Lower left /front vertex
               IF ( point_index(k,j,i) < 0 )  THEN
                  surfaces%npoints   = surfaces%npoints + 1
                  point_index(k,j,i) = point_index_count
                  point_index_count  = point_index_count + 1
                  surfaces%points(1,surfaces%npoints) = ( i - 0.5_wp ) * dx
                  surfaces%points(2,surfaces%npoints) = ( j - 0.5_wp ) * dy
                  surfaces%points(3,surfaces%npoints) = zw(k-1)
               ENDIF
!
!--            Upper / lower right index for north- and south-facing surfaces
               IF ( l == 0  .OR.  l == 1 )  THEN
                  IF ( point_index(k,j,i+1) < 0 )  THEN
                     surfaces%npoints     = surfaces%npoints + 1
                     point_index(k,j,i+1) = point_index_count
                     point_index_count    = point_index_count + 1
                     surfaces%points(1,surfaces%npoints) = ( i + 1 - 0.5_wp ) * dx
                     surfaces%points(2,surfaces%npoints) = ( j     - 0.5_wp ) * dy
                     surfaces%points(3,surfaces%npoints) = zw(k-1)
                  ENDIF
                  IF ( point_index(k+1,j,i+1) < 0 )  THEN
                     surfaces%npoints       = surfaces%npoints + 1
                     point_index(k+1,j,i+1) = point_index_count
                     point_index_count      = point_index_count + 1
                     surfaces%points(1,surfaces%npoints) = ( i + 1 - 0.5_wp ) * dx
                     surfaces%points(2,surfaces%npoints) = ( j     - 0.5_wp ) * dy
                     surfaces%points(3,surfaces%npoints) = zw(k)
                  ENDIF
!
!--            Upper / lower front index for east- and west-facing surfaces
               ELSEIF ( l == 2  .OR.  l == 3 )  THEN
                  IF ( point_index(k,j+1,i) < 0 )  THEN
                     surfaces%npoints     = surfaces%npoints + 1
                     point_index(k,j+1,i) = point_index_count
                     point_index_count    = point_index_count + 1
                     surfaces%points(1,surfaces%npoints) = ( i     - 0.5_wp ) * dx
                     surfaces%points(2,surfaces%npoints) = ( j + 1 - 0.5_wp ) * dy
                     surfaces%points(3,surfaces%npoints) = zw(k-1)
                  ENDIF
                  IF ( point_index(k+1,j+1,i) < 0 )  THEN
                     surfaces%npoints       = surfaces%npoints + 1
                     point_index(k+1,j+1,i) = point_index_count
                     point_index_count      = point_index_count + 1
                     surfaces%points(1,surfaces%npoints) = ( i     - 0.5_wp ) * dx
                     surfaces%points(2,surfaces%npoints) = ( j + 1 - 0.5_wp ) * dy
                     surfaces%points(3,surfaces%npoints) = zw(k)
                  ENDIF
               ENDIF
!
!--            Upper left / front vertex
               IF ( point_index(k+1,j,i) < 0 )  THEN
                  surfaces%npoints     = surfaces%npoints + 1
                  point_index(k+1,j,i) = point_index_count
                  point_index_count    = point_index_count + 1
                  surfaces%points(1,surfaces%npoints) = ( i - 0.5_wp ) * dx
                  surfaces%points(2,surfaces%npoints) = ( j - 0.5_wp ) * dy
                  surfaces%points(3,surfaces%npoints) = zw(k)
               ENDIF

               npg = npg + 1
               IF ( l == 0  .OR.  l == 1 )  THEN
                  surfaces%polygons(1,npg)   = 4
                  surfaces%polygons(2,npg)   = point_index(k,j,i)
                  surfaces%polygons(3,npg)   = point_index(k,j,i+1)
                  surfaces%polygons(4,npg)   = point_index(k+1,j,i+1)
                  surfaces%polygons(5,npg)   = point_index(k+1,j,i)
               ELSE
                  surfaces%polygons(1,npg)   = 4
                  surfaces%polygons(2,npg)   = point_index(k,j,i)
                  surfaces%polygons(3,npg)   = point_index(k,j+1,i)
                  surfaces%polygons(4,npg)   = point_index(k+1,j+1,i)
                  surfaces%polygons(5,npg)   = point_index(k+1,j,i)
               ENDIF

            ENDDO

            DO  m = 1, surf_lsm_v(l)%ns
               i = surf_lsm_v(l)%i(m) + MERGE( surf_lsm_v(l)%ioff, 0, l == 3 )
               j = surf_lsm_v(l)%j(m) + MERGE( surf_lsm_v(l)%joff, 0, l == 1 )
               k = surf_lsm_v(l)%k(m) + surf_lsm_v(l)%koff
!
!--            Lower left /front vertex
               IF ( point_index(k,j,i) < 0 )  THEN
                  surfaces%npoints   = surfaces%npoints + 1
                  point_index(k,j,i) = point_index_count
                  point_index_count  = point_index_count + 1
                  surfaces%points(1,surfaces%npoints) = ( i - 0.5_wp ) * dx
                  surfaces%points(2,surfaces%npoints) = ( j - 0.5_wp ) * dy
                  surfaces%points(3,surfaces%npoints) = zw(k-1)
               ENDIF
!
!--            Upper / lower right index for north- and south-facing surfaces
               IF ( l == 0  .OR.  l == 1 )  THEN
                  IF ( point_index(k,j,i+1) < 0 )  THEN
                     surfaces%npoints     = surfaces%npoints + 1
                     point_index(k,j,i+1) = point_index_count
                     point_index_count    = point_index_count + 1
                     surfaces%points(1,surfaces%npoints) = ( i + 1 - 0.5_wp ) * dx
                     surfaces%points(2,surfaces%npoints) = ( j     - 0.5_wp ) * dy
                     surfaces%points(3,surfaces%npoints) = zw(k-1)
                  ENDIF
                  IF ( point_index(k+1,j,i+1) < 0 )  THEN
                     surfaces%npoints       = surfaces%npoints + 1
                     point_index(k+1,j,i+1) = point_index_count
                     point_index_count      = point_index_count + 1
                     surfaces%points(1,surfaces%npoints) = ( i + 1 - 0.5_wp ) * dx
                     surfaces%points(2,surfaces%npoints) = ( j     - 0.5_wp ) * dy
                     surfaces%points(3,surfaces%npoints) = zw(k)
                  ENDIF
!
!--            Upper / lower front index for east- and west-facing surfaces
               ELSEIF ( l == 2  .OR.  l == 3 )  THEN
                  IF ( point_index(k,j+1,i) < 0 )  THEN
                     surfaces%npoints     = surfaces%npoints + 1
                     point_index(k,j+1,i) = point_index_count
                     point_index_count    = point_index_count + 1
                     surfaces%points(1,surfaces%npoints) = ( i     - 0.5_wp ) * dx
                     surfaces%points(2,surfaces%npoints) = ( j + 1 - 0.5_wp ) * dy
                     surfaces%points(3,surfaces%npoints) = zw(k-1)
                  ENDIF
                  IF ( point_index(k+1,j+1,i) < 0 )  THEN
                     surfaces%npoints       = surfaces%npoints + 1
                     point_index(k+1,j+1,i) = point_index_count
                     point_index_count      = point_index_count + 1
                     surfaces%points(1,surfaces%npoints) = ( i     - 0.5_wp ) * dx
                     surfaces%points(2,surfaces%npoints) = ( j + 1 - 0.5_wp ) * dy
                     surfaces%points(3,surfaces%npoints) = zw(k)
                  ENDIF
               ENDIF
!
!--            Upper left / front vertex
               IF ( point_index(k+1,j,i) < 0 )  THEN
                  surfaces%npoints     = surfaces%npoints + 1
                  point_index(k+1,j,i) = point_index_count
                  point_index_count    = point_index_count + 1
                  surfaces%points(1,surfaces%npoints) = ( i - 0.5_wp ) * dx
                  surfaces%points(2,surfaces%npoints) = ( j - 0.5_wp ) * dy
                  surfaces%points(3,surfaces%npoints) = zw(k)
               ENDIF

               npg = npg + 1
               IF ( l == 0  .OR.  l == 1 )  THEN
                  surfaces%polygons(1,npg)   = 4
                  surfaces%polygons(2,npg)   = point_index(k,j,i)
                  surfaces%polygons(3,npg)   = point_index(k,j,i+1)
                  surfaces%polygons(4,npg)   = point_index(k+1,j,i+1)
                  surfaces%polygons(5,npg)   = point_index(k+1,j,i)
               ELSE
                  surfaces%polygons(1,npg)   = 4
                  surfaces%polygons(2,npg)   = point_index(k,j,i)
                  surfaces%polygons(3,npg)   = point_index(k,j+1,i)
                  surfaces%polygons(4,npg)   = point_index(k+1,j+1,i)
                  surfaces%polygons(5,npg)   = point_index(k+1,j,i)
               ENDIF
            ENDDO
            DO  m = 1, surf_usm_v(l)%ns
               i = surf_usm_v(l)%i(m) + MERGE( surf_usm_v(l)%ioff, 0, l == 3 )
               j = surf_usm_v(l)%j(m) + MERGE( surf_usm_v(l)%joff, 0, l == 1 )
               k = surf_usm_v(l)%k(m) + surf_usm_v(l)%koff
!
!--            Lower left /front vertex
               IF ( point_index(k,j,i) < 0 )  THEN
                  surfaces%npoints   = surfaces%npoints + 1
                  point_index(k,j,i) = point_index_count
                  point_index_count  = point_index_count + 1
                  surfaces%points(1,surfaces%npoints) = ( i - 0.5_wp ) * dx
                  surfaces%points(2,surfaces%npoints) = ( j - 0.5_wp ) * dy
                  surfaces%points(3,surfaces%npoints) = zw(k-1)
               ENDIF
!
!--            Upper / lower right index for north- and south-facing surfaces
               IF ( l == 0  .OR.  l == 1 )  THEN
                  IF ( point_index(k,j,i+1) < 0 )  THEN
                     surfaces%npoints     = surfaces%npoints + 1
                     point_index(k,j,i+1) = point_index_count
                     point_index_count    = point_index_count + 1
                     surfaces%points(1,surfaces%npoints) = ( i + 1 - 0.5_wp ) * dx
                     surfaces%points(2,surfaces%npoints) = ( j     - 0.5_wp ) * dy
                     surfaces%points(3,surfaces%npoints) = zw(k-1)
                  ENDIF
                  IF ( point_index(k+1,j,i+1) < 0 )  THEN
                     surfaces%npoints       = surfaces%npoints + 1
                     point_index(k+1,j,i+1) = point_index_count
                     point_index_count      = point_index_count + 1
                     surfaces%points(1,surfaces%npoints) = ( i + 1 - 0.5_wp ) * dx
                     surfaces%points(2,surfaces%npoints) = ( j     - 0.5_wp ) * dy
                     surfaces%points(3,surfaces%npoints) = zw(k)
                  ENDIF
!
!--            Upper / lower front index for east- and west-facing surfaces
               ELSEIF ( l == 2  .OR.  l == 3 )  THEN
                  IF ( point_index(k,j+1,i) < 0 )  THEN
                     surfaces%npoints     = surfaces%npoints + 1
                     point_index(k,j+1,i) = point_index_count
                     point_index_count    = point_index_count + 1
                     surfaces%points(1,surfaces%npoints) = ( i     - 0.5_wp ) * dx
                     surfaces%points(2,surfaces%npoints) = ( j + 1 - 0.5_wp ) * dy
                     surfaces%points(3,surfaces%npoints) = zw(k-1)
                  ENDIF
                  IF ( point_index(k+1,j+1,i) < 0 )  THEN
                     surfaces%npoints       = surfaces%npoints + 1
                     point_index(k+1,j+1,i) = point_index_count
                     point_index_count      = point_index_count + 1
                     surfaces%points(1,surfaces%npoints) = ( i     - 0.5_wp ) * dx
                     surfaces%points(2,surfaces%npoints) = ( j + 1 - 0.5_wp ) * dy
                     surfaces%points(3,surfaces%npoints) = zw(k)
                  ENDIF
               ENDIF
!
!--            Upper left / front vertex
               IF ( point_index(k+1,j,i) < 0 )  THEN
                  surfaces%npoints     = surfaces%npoints + 1
                  point_index(k+1,j,i) = point_index_count
                  point_index_count    = point_index_count + 1
                  surfaces%points(1,surfaces%npoints) = ( i - 0.5_wp ) * dx
                  surfaces%points(2,surfaces%npoints) = ( j - 0.5_wp ) * dy
                  surfaces%points(3,surfaces%npoints) = zw(k)
               ENDIF

               npg = npg + 1
               IF ( l == 0  .OR.  l == 1 )  THEN
                  surfaces%polygons(1,npg)   = 4
                  surfaces%polygons(2,npg)   = point_index(k,j,i)
                  surfaces%polygons(3,npg)   = point_index(k,j,i+1)
                  surfaces%polygons(4,npg)   = point_index(k+1,j,i+1)
                  surfaces%polygons(5,npg)   = point_index(k+1,j,i)
               ELSE
                  surfaces%polygons(1,npg)   = 4
                  surfaces%polygons(2,npg)   = point_index(k,j,i)
                  surfaces%polygons(3,npg)   = point_index(k,j+1,i)
                  surfaces%polygons(4,npg)   = point_index(k+1,j+1,i)
                  surfaces%polygons(5,npg)   = point_index(k+1,j,i)
               ENDIF
            ENDDO

         ENDDO
!
!--      Deallocate temporary dummy variable
         DEALLOCATE ( point_index )
!
!--      Sum-up total number of vertices on domain. This
!--      will be needed for post-processing.
         surfaces%npoints_total = 0
#if defined( __parallel )
          CALL MPI_ALLREDUCE( surfaces%npoints, surfaces%npoints_total, 1,     &
                              MPI_INTEGER, MPI_SUM, comm2d, ierr )
#else
          surfaces%npoints_total = surfaces%npoints
#endif
       ENDIF
!
!--    If output to netcdf is enabled, set-up the coordinate arrays that
!--    unambiguously describe the position and orientation of each surface
!--    element.
       IF ( to_netcdf )  THEN
!
!--       Allocate local coordinate arrays
          ALLOCATE( surfaces%s(1:surfaces%ns)       )
          ALLOCATE( surfaces%xs(1:surfaces%ns)      )
          ALLOCATE( surfaces%ys(1:surfaces%ns)      )
          ALLOCATE( surfaces%zs(1:surfaces%ns)      )
          ALLOCATE( surfaces%azimuth(1:surfaces%ns) )
          ALLOCATE( surfaces%zenith(1:surfaces%ns)  )
          ALLOCATE( surfaces%es_utm(1:surfaces%ns)  )
          ALLOCATE( surfaces%ns_utm(1:surfaces%ns)  )
!
!--       Gather the number of surface on each processor, in order to number
!--       the surface elements in ascending order with respect to the total
!--       number of surfaces in the domain.
#if defined( __parallel )
          CALL MPI_ALLGATHER( surfaces%ns, 1, MPI_INTEGER,                     &
                              num_surfaces_on_pe, 1, MPI_INTEGER, comm2d, ierr  )
#else
          num_surfaces_on_pe = surfaces%ns
#endif
!
!--       First, however, determine the offset where couting of the surfaces
!--       should start (the sum of surfaces on all PE's with lower MPI rank).
          i           = 0
          start_count = 1
          DO WHILE ( i < myid  .AND.  i <= SIZE( num_surfaces_on_pe ) )
             start_count = start_count + num_surfaces_on_pe(i)
             i           = i + 1
          ENDDO
!
!--       Set coordinate arrays. For horizontal surfaces, azimuth
!--       angles are not defined (fill value). Zenith angle is 0 (180) for
!--       upward (downward)-facing surfaces.
          i  = start_count
          mm = 1
          DO  m = 1, surf_def_h(0)%ns
             surfaces%s(mm)       = i
             surfaces%xs(mm)      = ( surf_def_h(0)%i(m) + 0.5_wp ) * dx
             surfaces%ys(mm)      = ( surf_def_h(0)%j(m) + 0.5_wp ) * dy
             surfaces%zs(mm)      = zw(surf_def_h(0)%k(m)+surf_def_h(0)%koff)
             surfaces%azimuth(mm) = surfaces%fillvalue
             surfaces%zenith(mm)  = 0.0
             i                    = i + 1
             mm                   = mm + 1
          ENDDO
          DO  m = 1, surf_lsm_h%ns
             surfaces%s(mm)       = i
             surfaces%xs(mm)      = ( surf_lsm_h%i(m) + 0.5_wp ) * dx
             surfaces%ys(mm)      = ( surf_lsm_h%j(m) + 0.5_wp ) * dy
             surfaces%zs(mm)      = zw(surf_lsm_h%k(m)+surf_lsm_h%koff)
             surfaces%azimuth(mm) = surfaces%fillvalue
             surfaces%zenith(mm)  = 0.0
             i                    = i + 1
             mm                   = mm + 1
          ENDDO
          DO  m = 1, surf_usm_h%ns
             surfaces%s(mm)       = i
             surfaces%xs(mm)      = ( surf_usm_h%i(m) + 0.5_wp ) * dx
             surfaces%ys(mm)      = ( surf_usm_h%j(m) + 0.5_wp ) * dy
             surfaces%zs(mm)      = zw(surf_usm_h%k(m)+surf_usm_h%koff)
             surfaces%azimuth(mm) = surfaces%fillvalue
             surfaces%zenith(mm)  = 0.0
             i                    = i + 1
             mm                   = mm + 1
          ENDDO
          DO  m = 1, surf_def_h(1)%ns
             surfaces%s(mm)       = i
             surfaces%xs(mm)      = ( surf_def_h(1)%i(m) + 0.5_wp ) * dx
             surfaces%ys(mm)      = ( surf_def_h(1)%j(m) + 0.5_wp ) * dy
             surfaces%zs(mm)      = zw(surf_def_h(1)%k(m)+surf_def_h(1)%koff)
             surfaces%azimuth(mm) = surfaces%fillvalue
             surfaces%zenith(mm)  = 180.0
             i                    = i + 1
             mm                   = mm + 1
          ENDDO
!
!--       For vertical surfaces, zenith angles are not defined (fill value).
!--       Azimuth angle: northward (0), eastward (90), southward (180), 
!--       westward (270). 
!--       Note, for vertical surfaces, zenith angles are 90.0_wp. 
          DO  l = 0, 3
             IF ( l == 0 )  THEN
                az    = 0.0_wp
                off_x = 0.5_wp
                off_y = 0.0_wp
             ELSEIF ( l == 1 )  THEN
                az    = 180.0_wp
                off_x = 0.5_wp
                off_y = 1.0_wp
             ELSEIF ( l == 2 )  THEN
                az    = 90.0_wp
                off_x = 0.0_wp
                off_y = 0.5_wp
             ELSEIF ( l == 3 )  THEN
                az    = 270.0_wp
                off_x = 1.0_wp
                off_y = 0.5_wp
             ENDIF

             DO  m = 1, surf_def_v(l)%ns
                surfaces%s(mm)       = i
                surfaces%xs(mm)      = ( surf_def_v(l)%i(m) + off_x ) * dx
                surfaces%ys(mm)      = ( surf_def_v(l)%j(m) + off_y ) * dy
                surfaces%zs(mm)      = zu(surf_def_v(l)%k(m))
                surfaces%azimuth(mm) = az
                surfaces%zenith(mm)  = 90.0_wp
                i                    = i + 1
                mm                   = mm + 1
             ENDDO
             DO  m = 1, surf_lsm_v(l)%ns
                surfaces%s(mm)       = i
                surfaces%xs(mm)      = ( surf_lsm_v(l)%i(m) + off_x ) * dx
                surfaces%ys(mm)      = ( surf_lsm_v(l)%j(m) + off_y ) * dy
                surfaces%zs(mm)      = zu(surf_lsm_v(l)%k(m))
                surfaces%azimuth(mm) = az
                surfaces%zenith(mm)  = 90.0_wp
                i                    = i + 1
                mm                   = mm + 1
             ENDDO
             DO  m = 1, surf_usm_v(l)%ns
                surfaces%s(mm)       = i
                surfaces%xs(mm)      = ( surf_usm_v(l)%i(m) + off_x ) * dx
                surfaces%ys(mm)      = ( surf_usm_v(l)%j(m) + off_y ) * dy
                surfaces%zs(mm)      = zu(surf_usm_v(l)%k(m))
                surfaces%azimuth(mm) = az
                surfaces%zenith(mm)  = 90.0_wp
                i                    = i + 1
                mm                   = mm + 1
             ENDDO
          ENDDO
!
!--       Finally, define UTM coordinates, which are the x/y-coordinates
!--       plus the origin (lower-left coordinate of the model domain).
          surfaces%es_utm = surfaces%xs + init_model%origin_x
          surfaces%ns_utm = surfaces%ys + init_model%origin_y
!
!--       Initialize NetCDF data output. Please note, local start position for
!--       the surface elements in the NetCDF file is surfaces%s(1), while
!--       the number of surfaces on the subdomain is given by surfaces%ns.
#if defined( __netcdf4_parallel )

!
!--       Calculate number of time steps to be output
          ntdim_surf(0) = dosurf_time_count(0) + CEILING(                      &
                        ( end_time - MAX(                                      &
                            MERGE( skip_time_dosurf,                           &
                                   skip_time_dosurf + spinup_time,             &
                                   data_output_during_spinup ),                &
                            simulated_time_at_begin )                          &
                        ) / dt_dosurf )

          ntdim_surf(1) = dosurf_time_count(1) + CEILING(                      &
                        ( end_time - MAX(                                      &
                            MERGE( skip_time_dosurf_av,                        &
                                   skip_time_dosurf_av + spinup_time,          &
                                   data_output_during_spinup ),                &
                            simulated_time_at_begin )                          &
                        ) / dt_dosurf_av )

!
!--       Create NetCDF4 files for parallel writing
          DO  av = 0, 1
!
!--          If there is no instantaneous data (av=0) or averaged data (av=1)
!--          requested, do not create the corresponding NetCDF file
             IF ( dosurf_no(av) == 0 ) CYCLE

             IF ( av == 0 )  THEN
                filename = 'SURFACE_DATA_NETCDF' // TRIM( coupling_char )
             ELSE
                filename = 'SURFACE_DATA_AV_NETCDF' // TRIM( coupling_char )
             ENDIF
!
!--          Open file using netCDF4/HDF5 format, parallel
             nc_stat = NF90_CREATE( TRIM(filename),                            &
                                    IOR( NF90_NOCLOBBER,                       &
                                    IOR( NF90_NETCDF4, NF90_MPIIO ) ),         &
                                    id_set_surf(av),                           &
                                    COMM = comm2d, INFO = MPI_INFO_NULL )
             CALL netcdf_handle_error( 'surface_data_output_mod', 5550 )

             !- Write some global attributes
             IF ( av == 0 )  THEN
                CALL netcdf_create_global_atts( id_set_surf(av),               &
                                                'surface-data',                &
                                                TRIM( run_description_header ),&
                                                5551 )
                time_average_text = ' '
             ELSE
                CALL netcdf_create_global_atts( id_set_surf(av),               &
                                                'surface-data_av',             &
                                                TRIM( run_description_header ),&
                                                5552 )
                WRITE ( time_average_text,'(F7.1,'' s avg'')' )  &
                   averaging_interval_surf
                nc_stat = NF90_PUT_ATT( id_set_surf(av), NF90_GLOBAL,          &
                                        'time_avg',                            &
                                        TRIM( time_average_text ) )
                CALL netcdf_handle_error( 'surface_data_output_mod', 5553 )
             ENDIF


!
!--          Define time coordinate for surface data.
!--          For parallel output the time dimension has to be limited
!--          (ntdim_surf), otherwise the performance drops significantly.
             CALL netcdf_create_dim( id_set_surf(av), 'time', ntdim_surf(av),  &
                                     id_dim_time_surf(av), 5554 )

             CALL netcdf_create_var( id_set_surf(av),                          &
                                     (/ id_dim_time_surf(av) /),               &
                                     'time', NF90_DOUBLE,                      &
                                     id_var_time_surf(av),                     &
                                     'seconds since '//                        &
                                     TRIM(init_model%origin_time),             &
                                     'time', 5555, 5555, 5555 )

             CALL netcdf_create_att( id_set_surf(av), id_var_time_surf(av),    &
                                     'standard_name', 'time', 5556)

             CALL netcdf_create_att( id_set_surf(av), id_var_time_surf(av),    &
                                     'axis', 'T', 5557)
!
!--          Define spatial dimensions and coordinates:
!--          Define index of surface element
             CALL netcdf_create_dim( id_set_surf(av), 's', surfaces%ns_total,  &
                                     id_dim_s_surf(av), 5558 )
             CALL netcdf_create_var( id_set_surf(av), (/ id_dim_s_surf(av) /), &
                                     's', NF90_DOUBLE, id_var_s_surf(av),      &
                                     '1', 'number of surface element',         &
                                     5559, 5559, 5559 )
!
!--          Define x coordinate
             CALL netcdf_create_var( id_set_surf(av), (/ id_dim_s_surf(av) /), &
                                     'xs', NF90_DOUBLE, id_var_xs_surf(av),    &
                                     'meters',                                 &
                                     'distance to origin in x-direction',      &
                                     5561, 5561, 5561 )
!
!--           Define y coordinate
             CALL netcdf_create_var( id_set_surf(av), (/ id_dim_s_surf(av) /), &
                                     'ys', NF90_DOUBLE, id_var_ys_surf(av),    &
                                     'meters',                                 &
                                     'distance to origin in y-direction',      &
                                     5562, 5562, 5562 )
!
!--          Define z coordinate
             CALL netcdf_create_var( id_set_surf(av), (/ id_dim_s_surf(av) /), &
                                     'zs', NF90_DOUBLE, id_var_zs_surf(av),    &
                                     'meters', 'height', 5560, 5560, 5560 )
             CALL netcdf_create_att( id_set_surf(av), id_var_zs_surf(av),      &
                                     'standard_name', 'height', 5583 )

!
!--          Define UTM coordinates
             CALL netcdf_create_var( id_set_surf(av), (/ id_dim_s_surf(av) /), &
                                     'Es_UTM', NF90_DOUBLE,                    &
                                     id_var_etum_surf(av),                     &
                                     'meters', '', 5563, 5563, 5563 )
             CALL netcdf_create_var( id_set_surf(av), (/ id_dim_s_surf(av) /), &
                                     'Ns_UTM', NF90_DOUBLE,                    &
                                     id_var_nutm_surf(av),                     &
                                     'meters', '', 5564, 5564, 5564 )

!
!--          Define angles
             CALL netcdf_create_var( id_set_surf(av), (/ id_dim_s_surf(av) /), &
                                     'azimuth', NF90_DOUBLE,                   &
                                     id_var_azimuth_surf(av),                  &
                                     'degree', 'azimuth angle',                &
                                     5577, 5578, 5579,                         &
                                     fill = .TRUE. )
             CALL netcdf_create_att( id_set_surf(av), id_var_azimuth_surf(av), &
                                     'standard_name', 'surface_azimuth_angle', &
                                     5584 )

             CALL netcdf_create_var( id_set_surf(av), (/ id_dim_s_surf(av) /), &
                                     'zenith', NF90_DOUBLE,                    &
                                     id_var_zenith_surf(av),                   &
                                     'degree', '', 5580, 5581, 5582,           &
                                     fill = .TRUE. )
!
!--          Define the variables
             var_list = ';'
             i = 1

             DO WHILE ( dosurf(av,i)(1:1) /= ' ' )

                CALL netcdf_create_var( id_set_surf(av),(/ id_dim_s_surf(av),  &
                                        id_dim_time_surf(av) /), dosurf(av,i), &
                                        NF90_REAL4, id_var_dosurf(av,i),       &
                                        dosurf_unit(av,i), dosurf(av,i), 5565, &
                                        5565, 5565, .TRUE. )
!
!--                Set no fill for every variable to increase performance.
                nc_stat = NF90_DEF_VAR_FILL( id_set_surf(av),                  &
                                             id_var_dosurf(av,i),              &
                                             NF90_NOFILL, 0 )
                CALL netcdf_handle_error( 'surface_data_output_init', 5566 )
!
!--                Set collective io operations for parallel io
                nc_stat = NF90_VAR_PAR_ACCESS( id_set_surf(av),                &
                                               id_var_dosurf(av,i),            &
                                               NF90_COLLECTIVE )
                CALL netcdf_handle_error( 'surface_data_output_init', 5567 )
                var_list = TRIM( var_list ) // TRIM( dosurf(av,i) ) // ';'

                i = i + 1

             ENDDO
!
!--          Write the list of variables as global attribute (this is used by
!--          restart runs and by combine_plot_fields)
             nc_stat = NF90_PUT_ATT( id_set_surf(av), NF90_GLOBAL, 'VAR_LIST', &
                                     var_list )
             CALL netcdf_handle_error( 'surface_data_output_init', 5568 )

!
!--          Set general no fill, otherwise the performance drops significantly
!--          for parallel output.
             nc_stat = NF90_SET_FILL( id_set_surf(av), NF90_NOFILL, oldmode )
             CALL netcdf_handle_error( 'surface_data_output_init', 5569 )

!
!--          Leave netCDF define mode
             nc_stat = NF90_ENDDEF( id_set_surf(av) )
             CALL netcdf_handle_error( 'surface_data_output_init', 5570 )

!
!--          These data are only written by PE0 for parallel output to increase
!--          the performance.
             IF ( myid == 0 )  THEN
!
!--             Write data for surface indices
                ALLOCATE( netcdf_data_1d(1:surfaces%ns_total) )

                DO  i = 1, surfaces%ns_total
                   netcdf_data_1d(i) = i
                ENDDO

                nc_stat = NF90_PUT_VAR( id_set_surf(av), id_var_s_surf(av),    &
                                        netcdf_data_1d, start = (/ 1 /),       &
                                        count = (/ surfaces%ns_total /) )
                CALL netcdf_handle_error( 'surface_data_output_init', 5571 )

                DEALLOCATE( netcdf_data_1d )

             ENDIF

!
!--          Write surface positions to file
             nc_stat = NF90_PUT_VAR( id_set_surf(av), id_var_xs_surf(av),      &
                                     surfaces%xs, start = (/ surfaces%s(1) /), &
                                     count = (/ surfaces%ns /) )
             CALL netcdf_handle_error( 'surface_data_output_init', 5572 )

             nc_stat = NF90_PUT_VAR( id_set_surf(av), id_var_ys_surf(av),      &
                                     surfaces%ys, start = (/ surfaces%s(1) /), &
                                     count = (/ surfaces%ns /) )
             CALL netcdf_handle_error( 'surface_data_output_init', 5573 )

             nc_stat = NF90_PUT_VAR( id_set_surf(av), id_var_zs_surf(av),      &
                                     surfaces%zs, start = (/ surfaces%s(1) /), &
                                     count = (/ surfaces%ns /) )
             CALL netcdf_handle_error( 'surface_data_output_init', 5574 )

             nc_stat = NF90_PUT_VAR( id_set_surf(av), id_var_etum_surf(av),    &
                                     surfaces%es_utm,                          &
                                     start = (/ surfaces%s(1) /),              &
                                     count = (/ surfaces%ns /) )
             CALL netcdf_handle_error( 'surface_data_output_init', 5575 )

             nc_stat = NF90_PUT_VAR( id_set_surf(av), id_var_nutm_surf(av),    &
                                     surfaces%ns_utm,                          &
                                     start = (/ surfaces%s(1) /),              &
                                     count = (/ surfaces%ns /) )
             CALL netcdf_handle_error( 'surface_data_output_init', 5576 )

             nc_stat = NF90_PUT_VAR( id_set_surf(av), id_var_azimuth_surf(av), &
                                     surfaces%azimuth,                         &
                                     start = (/ surfaces%s(1) /),              &
                                     count = (/ surfaces%ns /) )
             CALL netcdf_handle_error( 'surface_data_output_init', 5585 )

             nc_stat = NF90_PUT_VAR( id_set_surf(av), id_var_zenith_surf(av),  &
                                     surfaces%zenith,                          &
                                     start = (/ surfaces%s(1) /),              &
                                     count = (/ surfaces%ns /) )
             CALL netcdf_handle_error( 'surface_data_output_init', 5586 )

          ENDDO
#endif

       ENDIF

   END SUBROUTINE surface_data_output_init

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Routine for controlling the data output. Surface data is collected from
!> different types of surfaces (default, natural, urban) and different
!> orientation and written to one 1D-output array. Further, NetCDF routines
!> are called to write the surface data in the respective NetCDF files.
!------------------------------------------------------------------------------!
   SUBROUTINE surface_data_output( av )

      USE control_parameters,                                                  &
          ONLY:  io_blocks, io_group, time_since_reference_point

#if defined( __parallel )
      USE pegrid,                                                              &
          ONLY:  comm2d, ierr
#endif


      IMPLICIT NONE

      CHARACTER(LEN=100) ::  trimvar = ' ' !< dummy for single output variable

      INTEGER(iwp) ::  av     !< id indicating average or non-average data output
      INTEGER(iwp) ::  i      !< loop index
      INTEGER(iwp) ::  n_out  !< counter variables for surface output

!
!--   Return, if nothing to output
      IF ( dosurf_no(av) == 0 )  RETURN
!
!--   In case of VTK output, check if binary files are open and write coordinates.
      IF ( to_vtk )  THEN

         CALL check_open( 25+av )

         IF ( .NOT. first_output(av) )  THEN
            DO  i = 0, io_blocks-1
               IF ( i == io_group )  THEN
                  WRITE ( 25+av )  surfaces%npoints
                  WRITE ( 25+av )  surfaces%npoints_total
                  WRITE ( 25+av )  surfaces%ns
                  WRITE ( 25+av )  surfaces%ns_total
                  WRITE ( 25+av )  surfaces%points
                  WRITE ( 25+av )  surfaces%polygons
               ENDIF
#if defined( __parallel )
               CALL MPI_BARRIER( comm2d, ierr )
#endif
               first_output(av) = .TRUE.
            ENDDO
         ENDIF
      ENDIF
!
!--   In case of NetCDF output, check if enough time steps are available in file
!--   and update time axis.
      IF ( to_netcdf )  THEN
#if defined( __netcdf4_parallel )
         IF ( dosurf_time_count(av) + 1 > ntdim_surf(av) )  THEN
            WRITE ( message_string, * )                               &
               'Output of surface data is not given at t=',           &
               time_since_reference_point, 's because the maximum ',  &
               'number of output time levels is exceeded.'
            CALL message( 'surface_data_output', 'PA0539', 0, 1, 0, 6, 0 )

            RETURN

         ENDIF
!
!--      Update the netCDF time axis
!--      In case of parallel output, this is only done by PE0 to increase the
!--      performance.
         dosurf_time_count(av) = dosurf_time_count(av) + 1
         IF ( myid == 0 )  THEN
            nc_stat = NF90_PUT_VAR( id_set_surf(av), id_var_time_surf(av),  &
                                    (/ time_since_reference_point /),       &
                                    start = (/ dosurf_time_count(av) /),    &
                                    count = (/ 1 /) )
            CALL netcdf_handle_error( 'surface_data_output', 6666 )
         ENDIF
#endif
      ENDIF

!
!--   Cycle through output quantities and write them to file.
      n_out = 0
      DO  WHILE ( dosurf(av,n_out+1)(1:1) /= ' ' )

         n_out   = n_out + 1
         trimvar = TRIM( dosurf(av,n_out) )
!
!--      Set the output array to the _FillValue in case it is not
!--      defined for each type of surface.
         surfaces%var_out = surfaces%fillvalue
         SELECT CASE ( trimvar )

            CASE ( 'us' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%us,          &
                                               surf_def_h(1)%us,               &
                                               surf_lsm_h%us,                  &
                                               surf_usm_h%us,                  &
                                               surf_def_v(0)%us,               &
                                               surf_lsm_v(0)%us,               &
                                               surf_usm_v(0)%us,               &
                                               surf_def_v(1)%us,               &
                                               surf_lsm_v(1)%us,               &
                                               surf_usm_v(1)%us,               &
                                               surf_def_v(2)%us,               &
                                               surf_lsm_v(2)%us,               &
                                               surf_usm_v(2)%us,               &
                                               surf_def_v(3)%us,               &
                                               surf_lsm_v(3)%us,               &
                                               surf_usm_v(3)%us )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'ts' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%ts,          &
                                               surf_def_h(1)%ts,               &
                                               surf_lsm_h%ts,                  &
                                               surf_usm_h%ts,                  &
                                               surf_def_v(0)%ts,               &
                                               surf_lsm_v(0)%ts,               &
                                               surf_usm_v(0)%ts,               &
                                               surf_def_v(1)%ts,               &
                                               surf_lsm_v(1)%ts,               &
                                               surf_usm_v(1)%ts,               &
                                               surf_def_v(2)%ts,               &
                                               surf_lsm_v(2)%ts,               &
                                               surf_usm_v(2)%ts,               &
                                               surf_def_v(3)%ts,               &
                                               surf_lsm_v(3)%ts,               &
                                               surf_usm_v(3)%ts )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'qs' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%qs,          &
                                               surf_def_h(1)%qs,               &
                                               surf_lsm_h%qs,                  &
                                               surf_usm_h%qs,                  &
                                               surf_def_v(0)%qs,               &
                                               surf_lsm_v(0)%qs,               &
                                               surf_usm_v(0)%qs,               &
                                               surf_def_v(1)%qs,               &
                                               surf_lsm_v(1)%qs,               &
                                               surf_usm_v(1)%qs,               &
                                               surf_def_v(2)%qs,               &
                                               surf_lsm_v(2)%qs,               &
                                               surf_usm_v(2)%qs,               &
                                               surf_def_v(3)%qs,               &
                                               surf_lsm_v(3)%qs,               &
                                               surf_usm_v(3)%qs )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'ss' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%ss,          &
                                               surf_def_h(1)%ss,               &
                                               surf_lsm_h%ss,                  &
                                               surf_usm_h%ss,                  &
                                               surf_def_v(0)%ss,               &
                                               surf_lsm_v(0)%ss,               &
                                               surf_usm_v(0)%ss,               &
                                               surf_def_v(1)%ss,               &
                                               surf_lsm_v(1)%ss,               &
                                               surf_usm_v(1)%ss,               &
                                               surf_def_v(2)%ss,               &
                                               surf_lsm_v(2)%ss,               &
                                               surf_usm_v(2)%ss,               &
                                               surf_def_v(3)%ss,               &
                                               surf_lsm_v(3)%ss,               &
                                               surf_usm_v(3)%ss )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'qcs' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%qcs,         &
                                               surf_def_h(1)%qcs,              &
                                               surf_lsm_h%qcs,                 &
                                               surf_usm_h%qcs,                 &
                                               surf_def_v(0)%qcs,              &
                                               surf_lsm_v(0)%qcs,              &
                                               surf_usm_v(0)%qcs,              &
                                               surf_def_v(1)%qcs,              &
                                               surf_lsm_v(1)%qcs,              &
                                               surf_usm_v(1)%qcs,              &
                                               surf_def_v(2)%qcs,              &
                                               surf_lsm_v(2)%qcs,              &
                                               surf_usm_v(2)%qcs,              &
                                               surf_def_v(3)%qcs,              &
                                               surf_lsm_v(3)%qcs,              &
                                               surf_usm_v(3)%qcs )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'ncs' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%ncs,         &
                                               surf_def_h(1)%ncs,              &
                                               surf_lsm_h%ncs,                 &
                                               surf_usm_h%ncs,                 &
                                               surf_def_v(0)%ncs,              &
                                               surf_lsm_v(0)%ncs,              &
                                               surf_usm_v(0)%ncs,              &
                                               surf_def_v(1)%ncs,              &
                                               surf_lsm_v(1)%ncs,              &
                                               surf_usm_v(1)%ncs,              &
                                               surf_def_v(2)%ncs,              &
                                               surf_lsm_v(2)%ncs,              &
                                               surf_usm_v(2)%ncs,              &
                                               surf_def_v(3)%ncs,              &
                                               surf_lsm_v(3)%ncs,              &
                                               surf_usm_v(3)%ncs )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'qrs' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%qrs,         &
                                               surf_def_h(1)%qrs,              &
                                               surf_lsm_h%qrs,                 &
                                               surf_usm_h%qrs,                 &
                                               surf_def_v(0)%qrs,              &
                                               surf_lsm_v(0)%qrs,              &
                                               surf_usm_v(0)%qrs,              &
                                               surf_def_v(1)%qrs,              &
                                               surf_lsm_v(1)%qrs,              &
                                               surf_usm_v(1)%qrs,              &
                                               surf_def_v(2)%qrs,              &
                                               surf_lsm_v(2)%qrs,              &
                                               surf_usm_v(2)%qrs,              &
                                               surf_def_v(3)%qrs,              &
                                               surf_lsm_v(3)%qrs,              &
                                               surf_usm_v(3)%qrs )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'nrs' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%nrs,         &
                                               surf_def_h(1)%nrs,              &
                                               surf_lsm_h%nrs,                 &
                                               surf_usm_h%nrs,                 &
                                               surf_def_v(0)%nrs,              &
                                               surf_lsm_v(0)%nrs,              &
                                               surf_usm_v(0)%nrs,              &
                                               surf_def_v(1)%nrs,              &
                                               surf_lsm_v(1)%nrs,              &
                                               surf_usm_v(1)%nrs,              &
                                               surf_def_v(2)%nrs,              &
                                               surf_lsm_v(2)%nrs,              &
                                               surf_usm_v(2)%nrs,              &
                                               surf_def_v(3)%nrs,              &
                                               surf_lsm_v(3)%nrs,              &
                                               surf_usm_v(3)%nrs )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'ol' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%ol,          &
                                               surf_def_h(1)%ol,               &
                                               surf_lsm_h%ol,                  &
                                               surf_usm_h%ol,                  &
                                               surf_def_v(0)%ol,               &
                                               surf_lsm_v(0)%ol,               &
                                               surf_usm_v(0)%ol,               &
                                               surf_def_v(1)%ol,               &
                                               surf_lsm_v(1)%ol,               &
                                               surf_usm_v(1)%ol,               &
                                               surf_def_v(2)%ol,               &
                                               surf_lsm_v(2)%ol,               &
                                               surf_usm_v(2)%ol,               &
                                               surf_def_v(3)%ol,               &
                                               surf_lsm_v(3)%ol,               &
                                               surf_usm_v(3)%ol )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'z0' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%z0,          &
                                               surf_def_h(1)%z0,               &
                                               surf_lsm_h%z0,                  &
                                               surf_usm_h%z0,                  &
                                               surf_def_v(0)%z0,               &
                                               surf_lsm_v(0)%z0,               &
                                               surf_usm_v(0)%z0,               &
                                               surf_def_v(1)%z0,               &
                                               surf_lsm_v(1)%z0,               &
                                               surf_usm_v(1)%z0,               &
                                               surf_def_v(2)%z0,               &
                                               surf_lsm_v(2)%z0,               &
                                               surf_usm_v(2)%z0,               &
                                               surf_def_v(3)%z0,               &
                                               surf_lsm_v(3)%z0,               &
                                               surf_usm_v(3)%z0 )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'z0h' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%z0h,         &
                                               surf_def_h(1)%z0h,              &
                                               surf_lsm_h%z0h,                 &
                                               surf_usm_h%z0h,                 &
                                               surf_def_v(0)%z0h,              &
                                               surf_lsm_v(0)%z0h,              &
                                               surf_usm_v(0)%z0h,              &
                                               surf_def_v(1)%z0h,              &
                                               surf_lsm_v(1)%z0h,              &
                                               surf_usm_v(1)%z0h,              &
                                               surf_def_v(2)%z0h,              &
                                               surf_lsm_v(2)%z0h,              &
                                               surf_usm_v(2)%z0h,              &
                                               surf_def_v(3)%z0h,              &
                                               surf_lsm_v(3)%z0h,              &
                                               surf_usm_v(3)%z0h )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'z0q' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%z0q,         &
                                               surf_def_h(1)%z0q,              &
                                               surf_lsm_h%z0q,                 &
                                               surf_usm_h%z0q,                 &
                                               surf_def_v(0)%z0q,              &
                                               surf_lsm_v(0)%z0q,              &
                                               surf_usm_v(0)%z0q,              &
                                               surf_def_v(1)%z0q,              &
                                               surf_lsm_v(1)%z0q,              &
                                               surf_usm_v(1)%z0q,              &
                                               surf_def_v(2)%z0q,              &
                                               surf_lsm_v(2)%z0q,              &
                                               surf_usm_v(2)%z0q,              &
                                               surf_def_v(3)%z0q,              &
                                               surf_lsm_v(3)%z0q,              &
                                               surf_usm_v(3)%z0q )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'theta1' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%pt1,         &
                                               surf_def_h(1)%pt1,              &
                                               surf_lsm_h%pt1,                 &
                                               surf_usm_h%pt1,                 &
                                               surf_def_v(0)%pt1,              &
                                               surf_lsm_v(0)%pt1,              &
                                               surf_usm_v(0)%pt1,              &
                                               surf_def_v(1)%pt1,              &
                                               surf_lsm_v(1)%pt1,              &
                                               surf_usm_v(1)%pt1,              &
                                               surf_def_v(2)%pt1,              &
                                               surf_lsm_v(2)%pt1,              &
                                               surf_usm_v(2)%pt1,              &
                                               surf_def_v(3)%pt1,              &
                                               surf_lsm_v(3)%pt1,              &
                                               surf_usm_v(3)%pt1 )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'qv1' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%qv1,         &
                                               surf_def_h(1)%qv1,              &
                                               surf_lsm_h%qv1,                 &
                                               surf_usm_h%qv1,                 &
                                               surf_def_v(0)%qv1,              &
                                               surf_lsm_v(0)%qv1,              &
                                               surf_usm_v(0)%qv1,              &
                                               surf_def_v(1)%qv1,              &
                                               surf_lsm_v(1)%qv1,              &
                                               surf_usm_v(1)%qv1,              &
                                               surf_def_v(2)%qv1,              &
                                               surf_lsm_v(2)%qv1,              &
                                               surf_usm_v(2)%qv1,              &
                                               surf_def_v(3)%qv1,              &
                                               surf_lsm_v(3)%qv1,              &
                                               surf_usm_v(3)%qv1 )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'thetav1' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%vpt1,        &
                                               surf_def_h(1)%vpt1,             &
                                               surf_lsm_h%vpt1,                &
                                               surf_usm_h%vpt1,                &
                                               surf_def_v(0)%vpt1,             &
                                               surf_lsm_v(0)%vpt1,             &
                                               surf_usm_v(0)%vpt1,             &
                                               surf_def_v(1)%vpt1,             &
                                               surf_lsm_v(1)%vpt1,             &
                                               surf_usm_v(1)%vpt1,             &
                                               surf_def_v(2)%vpt1,             &
                                               surf_lsm_v(2)%vpt1,             &
                                               surf_usm_v(2)%vpt1,             &
                                               surf_def_v(3)%vpt1,             &
                                               surf_lsm_v(3)%vpt1,             &
                                               surf_usm_v(3)%vpt1 )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'usws' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%usws,        &
                                               surf_def_h(1)%usws,             &
                                               surf_lsm_h%usws,                &
                                               surf_usm_h%usws,                &
                                               surf_def_v(0)%usws,             &
                                               surf_lsm_v(0)%usws,             &
                                               surf_usm_v(0)%usws,             &
                                               surf_def_v(1)%usws,             &
                                               surf_lsm_v(1)%usws,             &
                                               surf_usm_v(1)%usws,             &
                                               surf_def_v(2)%usws,             &
                                               surf_lsm_v(2)%usws,             &
                                               surf_usm_v(2)%usws,             &
                                               surf_def_v(3)%usws,             &
                                               surf_lsm_v(3)%usws,             &
                                               surf_usm_v(3)%usws )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'vsws' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%vsws,        &
                                               surf_def_h(1)%vsws,             &
                                               surf_lsm_h%vsws,                &
                                               surf_usm_h%vsws,                &
                                               surf_def_v(0)%vsws,             &
                                               surf_lsm_v(0)%vsws,             &
                                               surf_usm_v(0)%vsws,             &
                                               surf_def_v(1)%vsws,             &
                                               surf_lsm_v(1)%vsws,             &
                                               surf_usm_v(1)%vsws,             &
                                               surf_def_v(2)%vsws,             &
                                               surf_lsm_v(2)%vsws,             &
                                               surf_usm_v(2)%vsws,             &
                                               surf_def_v(3)%vsws,             &
                                               surf_lsm_v(3)%vsws,             &
                                               surf_usm_v(3)%vsws )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'shf' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%shf,         &
                                               surf_def_h(1)%shf,              &
                                               surf_lsm_h%shf,                 &
                                               surf_usm_h%shf,                 &
                                               surf_def_v(0)%shf,              &
                                               surf_lsm_v(0)%shf,              &
                                               surf_usm_v(0)%shf,              &
                                               surf_def_v(1)%shf,              &
                                               surf_lsm_v(1)%shf,              &
                                               surf_usm_v(1)%shf,              &
                                               surf_def_v(2)%shf,              &
                                               surf_lsm_v(2)%shf,              &
                                               surf_usm_v(2)%shf,              &
                                               surf_def_v(3)%shf,              &
                                               surf_lsm_v(3)%shf,              &
                                               surf_usm_v(3)%shf )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp
               ENDIF

            CASE ( 'qsws' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%qsws,        &
                                               surf_def_h(1)%qsws,             &
                                               surf_lsm_h%qsws,                &
                                               surf_usm_h%qsws,                &
                                               surf_def_v(0)%qsws,             &
                                               surf_lsm_v(0)%qsws,             &
                                               surf_usm_v(0)%qsws,             &
                                               surf_def_v(1)%qsws,             &
                                               surf_lsm_v(1)%qsws,             &
                                               surf_usm_v(1)%qsws,             &
                                               surf_def_v(2)%qsws,             &
                                               surf_lsm_v(2)%qsws,             &
                                               surf_usm_v(2)%qsws,             &
                                               surf_def_v(3)%qsws,             &
                                               surf_lsm_v(3)%qsws,             &
                                               surf_usm_v(3)%qsws )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'ssws' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%ssws,        &
                                               surf_def_h(1)%ssws,             &
                                               surf_lsm_h%ssws,                &
                                               surf_usm_h%ssws,                &
                                               surf_def_v(0)%ssws,             &
                                               surf_lsm_v(0)%ssws,             &
                                               surf_usm_v(0)%ssws,             &
                                               surf_def_v(1)%ssws,             &
                                               surf_lsm_v(1)%ssws,             &
                                               surf_usm_v(1)%ssws,             &
                                               surf_def_v(2)%ssws,             &
                                               surf_lsm_v(2)%ssws,             &
                                               surf_usm_v(2)%ssws,             &
                                               surf_def_v(3)%ssws,             &
                                               surf_lsm_v(3)%ssws,             &
                                               surf_usm_v(3)%ssws )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'qcsws' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%qcsws,       &
                                               surf_def_h(1)%qcsws,            &
                                               surf_lsm_h%qcsws,               &
                                               surf_usm_h%qcsws,               &
                                               surf_def_v(0)%qcsws,            &
                                               surf_lsm_v(0)%qcsws,            &
                                               surf_usm_v(0)%qcsws,            &
                                               surf_def_v(1)%qcsws,            &
                                               surf_lsm_v(1)%qcsws,            &
                                               surf_usm_v(1)%qcsws,            &
                                               surf_def_v(2)%qcsws,            &
                                               surf_lsm_v(2)%qcsws,            &
                                               surf_usm_v(2)%qcsws,            &
                                               surf_def_v(3)%qcsws,            &
                                               surf_lsm_v(3)%qcsws,            &
                                               surf_usm_v(3)%qcsws )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'ncsws' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%ncsws,       &
                                               surf_def_h(1)%ncsws,            &
                                               surf_lsm_h%ncsws,               &
                                               surf_usm_h%ncsws,               &
                                               surf_def_v(0)%ncsws,            &
                                               surf_lsm_v(0)%ncsws,            &
                                               surf_usm_v(0)%ncsws,            &
                                               surf_def_v(1)%ncsws,            &
                                               surf_lsm_v(1)%ncsws,            &
                                               surf_usm_v(1)%ncsws,            &
                                               surf_def_v(2)%ncsws,            &
                                               surf_lsm_v(2)%ncsws,            &
                                               surf_usm_v(2)%ncsws,            &
                                               surf_def_v(3)%ncsws,            &
                                               surf_lsm_v(3)%ncsws,            &
                                               surf_usm_v(3)%ncsws )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'qrsws' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%qrsws,       &
                                               surf_def_h(1)%qrsws,            &
                                               surf_lsm_h%qrsws,               &
                                               surf_usm_h%qrsws,               &
                                               surf_def_v(0)%qrsws,            &
                                               surf_lsm_v(0)%qrsws,            &
                                               surf_usm_v(0)%qrsws,            &
                                               surf_def_v(1)%qrsws,            &
                                               surf_lsm_v(1)%qrsws,            &
                                               surf_usm_v(1)%qrsws,            &
                                               surf_def_v(2)%qrsws,            &
                                               surf_lsm_v(2)%qrsws,            &
                                               surf_usm_v(2)%qrsws,            &
                                               surf_def_v(3)%qrsws,            &
                                               surf_lsm_v(3)%qrsws,            &
                                               surf_usm_v(3)%qrsws )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'nrsws' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%nrsws,       &
                                               surf_def_h(1)%nrsws,            &
                                               surf_lsm_h%nrsws,               &
                                               surf_usm_h%nrsws,               &
                                               surf_def_v(0)%nrsws,            &
                                               surf_lsm_v(0)%nrsws,            &
                                               surf_usm_v(0)%nrsws,            &
                                               surf_def_v(1)%nrsws,            &
                                               surf_lsm_v(1)%nrsws,            &
                                               surf_usm_v(1)%nrsws,            &
                                               surf_def_v(2)%nrsws,            &
                                               surf_lsm_v(2)%nrsws,            &
                                               surf_usm_v(2)%nrsws,            &
                                               surf_def_v(3)%nrsws,            &
                                               surf_lsm_v(3)%nrsws,            &
                                               surf_usm_v(3)%nrsws )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'sasws' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%sasws,       &
                                               surf_def_h(1)%sasws,            &
                                               surf_lsm_h%sasws,               &
                                               surf_usm_h%sasws,               &
                                               surf_def_v(0)%sasws,            &
                                               surf_lsm_v(0)%sasws,            &
                                               surf_usm_v(0)%sasws,            &
                                               surf_def_v(1)%sasws,            &
                                               surf_lsm_v(1)%sasws,            &
                                               surf_usm_v(1)%sasws,            &
                                               surf_def_v(2)%sasws,            &
                                               surf_lsm_v(2)%sasws,            &
                                               surf_usm_v(2)%sasws,            &
                                               surf_def_v(3)%sasws,            &
                                               surf_lsm_v(3)%sasws,            &
                                               surf_usm_v(3)%sasws )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'q_surface' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%q_surface,   &
                                               surf_def_h(1)%q_surface,        &
                                               surf_lsm_h%q_surface,           &
                                               surf_usm_h%q_surface,           &
                                               surf_def_v(0)%q_surface,        &
                                               surf_lsm_v(0)%q_surface,        &
                                               surf_usm_v(0)%q_surface,        &
                                               surf_def_v(1)%q_surface,        &
                                               surf_lsm_v(1)%q_surface,        &
                                               surf_usm_v(1)%q_surface,        &
                                               surf_def_v(2)%q_surface,        &
                                               surf_lsm_v(2)%q_surface,        &
                                               surf_usm_v(2)%q_surface,        &
                                               surf_def_v(3)%q_surface,        &
                                               surf_lsm_v(3)%q_surface,        &
                                               surf_usm_v(3)%q_surface )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'theta_surface' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%pt_surface,  &
                                               surf_def_h(1)%pt_surface,       &
                                               surf_lsm_h%pt_surface,          &
                                               surf_usm_h%pt_surface,          &
                                               surf_def_v(0)%pt_surface,       &
                                               surf_lsm_v(0)%pt_surface,       &
                                               surf_usm_v(0)%pt_surface,       &
                                               surf_def_v(1)%pt_surface,       &
                                               surf_lsm_v(1)%pt_surface,       &
                                               surf_usm_v(1)%pt_surface,       &
                                               surf_def_v(2)%pt_surface,       &
                                               surf_lsm_v(2)%pt_surface,       &
                                               surf_usm_v(2)%pt_surface,       &
                                               surf_def_v(3)%pt_surface,       &
                                               surf_lsm_v(3)%pt_surface,       &
                                               surf_usm_v(3)%pt_surface )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'thetav_surface' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%vpt_surface, &
                                               surf_def_h(1)%vpt_surface,      &
                                               surf_lsm_h%vpt_surface,         &
                                               surf_usm_h%vpt_surface,         &
                                               surf_def_v(0)%vpt_surface,      &
                                               surf_lsm_v(0)%vpt_surface,      &
                                               surf_usm_v(0)%vpt_surface,      &
                                               surf_def_v(1)%vpt_surface,      &
                                               surf_lsm_v(1)%vpt_surface,      &
                                               surf_usm_v(1)%vpt_surface,      &
                                               surf_def_v(2)%vpt_surface,      &
                                               surf_lsm_v(2)%vpt_surface,      &
                                               surf_usm_v(2)%vpt_surface,      &
                                               surf_def_v(3)%vpt_surface,      &
                                               surf_lsm_v(3)%vpt_surface,      &
                                               surf_usm_v(3)%vpt_surface)
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'rad_net' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%rad_net,     &
                                               surf_def_h(1)%rad_net,          &
                                               surf_lsm_h%rad_net,             &
                                               surf_usm_h%rad_net,             &
                                               surf_def_v(0)%rad_net,          &
                                               surf_lsm_v(0)%rad_net,          &
                                               surf_usm_v(0)%rad_net,          &
                                               surf_def_v(1)%rad_net,          &
                                               surf_lsm_v(1)%rad_net,          &
                                               surf_usm_v(1)%rad_net,          &
                                               surf_def_v(2)%rad_net,          &
                                               surf_lsm_v(2)%rad_net,          &
                                               surf_usm_v(2)%rad_net,          &
                                               surf_def_v(3)%rad_net,          &
                                               surf_lsm_v(3)%rad_net,          &
                                               surf_usm_v(3)%rad_net )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'rad_lw_in' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%rad_lw_in,   &
                                               surf_def_h(1)%rad_lw_in,        &
                                               surf_lsm_h%rad_lw_in,           &
                                               surf_usm_h%rad_lw_in,           &
                                               surf_def_v(0)%rad_lw_in,        &
                                               surf_lsm_v(0)%rad_lw_in,        &
                                               surf_usm_v(0)%rad_lw_in,        &
                                               surf_def_v(1)%rad_lw_in,        &
                                               surf_lsm_v(1)%rad_lw_in,        &
                                               surf_usm_v(1)%rad_lw_in,        &
                                               surf_def_v(2)%rad_lw_in,        &
                                               surf_lsm_v(2)%rad_lw_in,        &
                                               surf_usm_v(2)%rad_lw_in,        &
                                               surf_def_v(3)%rad_lw_in,        &
                                               surf_lsm_v(3)%rad_lw_in,        &
                                               surf_usm_v(3)%rad_lw_in )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'rad_lw_out' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%rad_lw_out,  &
                                               surf_def_h(1)%rad_lw_out,       &
                                               surf_lsm_h%rad_lw_out,          &
                                               surf_usm_h%rad_lw_out,          &
                                               surf_def_v(0)%rad_lw_out,       &
                                               surf_lsm_v(0)%rad_lw_out,       &
                                               surf_usm_v(0)%rad_lw_out,       &
                                               surf_def_v(1)%rad_lw_out,       &
                                               surf_lsm_v(1)%rad_lw_out,       &
                                               surf_usm_v(1)%rad_lw_out,       &
                                               surf_def_v(2)%rad_lw_out,       &
                                               surf_lsm_v(2)%rad_lw_out,       &
                                               surf_usm_v(2)%rad_lw_out,       &
                                               surf_def_v(3)%rad_lw_out,       &
                                               surf_lsm_v(3)%rad_lw_out,       &
                                               surf_usm_v(3)%rad_lw_out )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'rad_sw_in' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%rad_sw_in,   &
                                               surf_def_h(1)%rad_sw_in,        &
                                               surf_lsm_h%rad_sw_in,           &
                                               surf_usm_h%rad_sw_in,           &
                                               surf_def_v(0)%rad_sw_in,        &
                                               surf_lsm_v(0)%rad_sw_in,        &
                                               surf_usm_v(0)%rad_sw_in,        &
                                               surf_def_v(1)%rad_sw_in,        &
                                               surf_lsm_v(1)%rad_sw_in,        &
                                               surf_usm_v(1)%rad_sw_in,        &
                                               surf_def_v(2)%rad_sw_in,        &
                                               surf_lsm_v(2)%rad_sw_in,        &
                                               surf_usm_v(2)%rad_sw_in,        &
                                               surf_def_v(3)%rad_sw_in,        &
                                               surf_lsm_v(3)%rad_sw_in,        &
                                               surf_usm_v(3)%rad_sw_in )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'rad_sw_out' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%rad_sw_out,  &
                                               surf_def_h(1)%rad_sw_out,       &
                                               surf_lsm_h%rad_sw_out,          &
                                               surf_usm_h%rad_sw_out,          &
                                               surf_def_v(0)%rad_sw_out,       &
                                               surf_lsm_v(0)%rad_sw_out,       &
                                               surf_usm_v(0)%rad_sw_out,       &
                                               surf_def_v(1)%rad_sw_out,       &
                                               surf_lsm_v(1)%rad_sw_out,       &
                                               surf_usm_v(1)%rad_sw_out,       &
                                               surf_def_v(2)%rad_sw_out,       &
                                               surf_lsm_v(2)%rad_sw_out,       &
                                               surf_usm_v(2)%rad_sw_out,       &
                                               surf_def_v(3)%rad_sw_out,       &
                                               surf_lsm_v(3)%rad_sw_out,       &
                                               surf_usm_v(3)%rad_sw_out )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'ghf' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%ghf,         &
                                               surf_def_h(1)%ghf,              &
                                               surf_lsm_h%ghf,                 &
                                               surf_usm_h%ghf,                 &
                                               surf_def_v(0)%ghf,              &
                                               surf_lsm_v(0)%ghf,              &
                                               surf_usm_v(0)%ghf,              &
                                               surf_def_v(1)%ghf,              &
                                               surf_lsm_v(1)%ghf,              &
                                               surf_usm_v(1)%ghf,              &
                                               surf_def_v(2)%ghf,              &
                                               surf_lsm_v(2)%ghf,              &
                                               surf_usm_v(2)%ghf,              &
                                               surf_def_v(3)%ghf,              &
                                               surf_lsm_v(3)%ghf,              &
                                               surf_usm_v(3)%ghf )
                                                                        ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'r_a' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%r_a,         &
                                               surf_def_h(1)%r_a,              &
                                               surf_lsm_h%r_a,                 &
                                               surf_usm_h%r_a,                 &
                                               surf_def_v(0)%r_a,              &
                                               surf_lsm_v(0)%r_a,              &
                                               surf_usm_v(0)%r_a,              &
                                               surf_def_v(1)%r_a,              &
                                               surf_lsm_v(1)%r_a,              &
                                               surf_usm_v(1)%r_a,              &
                                               surf_def_v(2)%r_a,              &
                                               surf_lsm_v(2)%r_a,              &
                                               surf_usm_v(2)%r_a,              &
                                               surf_def_v(3)%r_a,              &
                                               surf_lsm_v(3)%r_a,              &
                                               surf_usm_v(3)%r_a )
                                                                      ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'r_soil' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%r_soil,      &
                                               surf_def_h(1)%r_soil,           &
                                               surf_lsm_h%r_soil,              &
                                               surf_usm_h%r_soil,              &
                                               surf_def_v(0)%r_soil,           &
                                               surf_lsm_v(0)%r_soil,           &
                                               surf_usm_v(0)%r_soil,           &
                                               surf_def_v(1)%r_soil,           &
                                               surf_lsm_v(1)%r_soil,           &
                                               surf_usm_v(1)%r_soil,           &
                                               surf_def_v(2)%r_soil,           &
                                               surf_lsm_v(2)%r_soil,           &
                                               surf_usm_v(2)%r_soil,           &
                                               surf_def_v(3)%r_soil,           &
                                               surf_lsm_v(3)%r_soil,           &
                                               surf_usm_v(3)%r_soil )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'r_canopy' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%r_canopy,    &
                                               surf_def_h(1)%r_canopy,         &
                                               surf_lsm_h%r_canopy,            &
                                               surf_usm_h%r_canopy,            &
                                               surf_def_v(0)%r_canopy,         &
                                               surf_lsm_v(0)%r_canopy,         &
                                               surf_usm_v(0)%r_canopy,         &
                                               surf_def_v(1)%r_canopy,         &
                                               surf_lsm_v(1)%r_canopy,         &
                                               surf_usm_v(1)%r_canopy,         &
                                               surf_def_v(2)%r_canopy,         &
                                               surf_lsm_v(2)%r_canopy,         &
                                               surf_usm_v(2)%r_canopy,         &
                                               surf_def_v(3)%r_canopy,         &
                                               surf_lsm_v(3)%r_canopy,         &
                                               surf_usm_v(3)%r_canopy )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'r_s' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%r_s,         &
                                               surf_def_h(1)%r_s,              &
                                               surf_lsm_h%r_s,                 &
                                               surf_usm_h%r_s,                 &
                                               surf_def_v(0)%r_s,              &
                                               surf_lsm_v(0)%r_s,              &
                                               surf_usm_v(0)%r_s,              &
                                               surf_def_v(1)%r_s,              &
                                               surf_lsm_v(1)%r_s,              &
                                               surf_usm_v(1)%r_s,              &
                                               surf_def_v(2)%r_s,              &
                                               surf_lsm_v(2)%r_s,              &
                                               surf_usm_v(2)%r_s,              &
                                               surf_def_v(3)%r_s,              &
                                               surf_lsm_v(3)%r_s,              &
                                               surf_usm_v(3)%r_s )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'rad_sw_dir' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%rad_sw_dir,  &
                                               surf_def_h(1)%rad_sw_dir,       &
                                               surf_lsm_h%rad_sw_dir,          &
                                               surf_usm_h%rad_sw_dir,          &
                                               surf_def_v(0)%rad_sw_dir,       &
                                               surf_lsm_v(0)%rad_sw_dir,       &
                                               surf_usm_v(0)%rad_sw_dir,       &
                                               surf_def_v(1)%rad_sw_dir,       &
                                               surf_lsm_v(1)%rad_sw_dir,       &
                                               surf_usm_v(1)%rad_sw_dir,       &
                                               surf_def_v(2)%rad_sw_dir,       &
                                               surf_lsm_v(2)%rad_sw_dir,       &
                                               surf_usm_v(2)%rad_sw_dir,       &
                                               surf_def_v(3)%rad_sw_dir,       &
                                               surf_lsm_v(3)%rad_sw_dir,       &
                                               surf_usm_v(3)%rad_sw_dir )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'rad_sw_dif' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%rad_sw_dif,  &
                                               surf_def_h(1)%rad_sw_dif,       &
                                               surf_lsm_h%rad_sw_dif,          &
                                               surf_usm_h%rad_sw_dif,          &
                                               surf_def_v(0)%rad_sw_dif,       &
                                               surf_lsm_v(0)%rad_sw_dif,       &
                                               surf_usm_v(0)%rad_sw_dif,       &
                                               surf_def_v(1)%rad_sw_dif,       &
                                               surf_lsm_v(1)%rad_sw_dif,       &
                                               surf_usm_v(1)%rad_sw_dif,       &
                                               surf_def_v(2)%rad_sw_dif,       &
                                               surf_lsm_v(2)%rad_sw_dif,       &
                                               surf_usm_v(2)%rad_sw_dif,       &
                                               surf_def_v(3)%rad_sw_dif,       &
                                               surf_lsm_v(3)%rad_sw_dif,       &
                                               surf_usm_v(3)%rad_sw_dif )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'rad_sw_ref' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%rad_sw_ref,  &
                                               surf_def_h(1)%rad_sw_ref,       &
                                               surf_lsm_h%rad_sw_ref,          &
                                               surf_usm_h%rad_sw_ref,          &
                                               surf_def_v(0)%rad_sw_ref,       &
                                               surf_lsm_v(0)%rad_sw_ref,       &
                                               surf_usm_v(0)%rad_sw_ref,       &
                                               surf_def_v(1)%rad_sw_ref,       &
                                               surf_lsm_v(1)%rad_sw_ref,       &
                                               surf_usm_v(1)%rad_sw_ref,       &
                                               surf_def_v(2)%rad_sw_ref,       &
                                               surf_lsm_v(2)%rad_sw_ref,       &
                                               surf_usm_v(2)%rad_sw_ref,       &
                                               surf_def_v(3)%rad_sw_ref,       &
                                               surf_lsm_v(3)%rad_sw_ref,       &
                                               surf_usm_v(3)%rad_sw_ref )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'rad_sw_res' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%rad_sw_res,  &
                                               surf_def_h(1)%rad_sw_res,       &
                                               surf_lsm_h%rad_sw_res,          &
                                               surf_usm_h%rad_sw_res,          &
                                               surf_def_v(0)%rad_sw_res,       &
                                               surf_lsm_v(0)%rad_sw_res,       &
                                               surf_usm_v(0)%rad_sw_res,       &
                                               surf_def_v(1)%rad_sw_res,       &
                                               surf_lsm_v(1)%rad_sw_res,       &
                                               surf_usm_v(1)%rad_sw_res,       &
                                               surf_def_v(2)%rad_sw_res,       &
                                               surf_lsm_v(2)%rad_sw_res,       &
                                               surf_usm_v(2)%rad_sw_res,       &
                                               surf_def_v(3)%rad_sw_res,       &
                                               surf_lsm_v(3)%rad_sw_res,       &
                                               surf_usm_v(3)%rad_sw_res )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'rad_lw_dif' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%rad_lw_dif,  &
                                               surf_def_h(1)%rad_lw_dif,       &
                                               surf_lsm_h%rad_lw_dif,          &
                                               surf_usm_h%rad_lw_dif,          &
                                               surf_def_v(0)%rad_lw_dif,       &
                                               surf_lsm_v(0)%rad_lw_dif,       &
                                               surf_usm_v(0)%rad_lw_dif,       &
                                               surf_def_v(1)%rad_lw_dif,       &
                                               surf_lsm_v(1)%rad_lw_dif,       &
                                               surf_usm_v(1)%rad_lw_dif,       &
                                               surf_def_v(2)%rad_lw_dif,       &
                                               surf_lsm_v(2)%rad_lw_dif,       &
                                               surf_usm_v(2)%rad_lw_dif,       &
                                               surf_def_v(3)%rad_lw_dif,       &
                                               surf_lsm_v(3)%rad_lw_dif,       &
                                               surf_usm_v(3)%rad_lw_dif )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'rad_lw_ref' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%rad_lw_ref,  &
                                               surf_def_h(1)%rad_lw_ref,       &
                                               surf_lsm_h%rad_lw_ref,          &
                                               surf_usm_h%rad_lw_ref,          &
                                               surf_def_v(0)%rad_lw_ref,       &
                                               surf_lsm_v(0)%rad_lw_ref,       &
                                               surf_usm_v(0)%rad_lw_ref,       &
                                               surf_def_v(1)%rad_lw_ref,       &
                                               surf_lsm_v(1)%rad_lw_ref,       &
                                               surf_usm_v(1)%rad_lw_ref,       &
                                               surf_def_v(2)%rad_lw_ref,       &
                                               surf_lsm_v(2)%rad_lw_ref,       &
                                               surf_usm_v(2)%rad_lw_ref,       &
                                               surf_def_v(3)%rad_lw_ref,       &
                                               surf_lsm_v(3)%rad_lw_ref,       &
                                               surf_usm_v(3)%rad_lw_ref )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'rad_lw_res' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%rad_lw_res,  &
                                               surf_def_h(1)%rad_lw_res,       &
                                               surf_lsm_h%rad_lw_res,          &
                                               surf_usm_h%rad_lw_res,          &
                                               surf_def_v(0)%rad_lw_res,       &
                                               surf_lsm_v(0)%rad_lw_res,       &
                                               surf_usm_v(0)%rad_lw_res,       &
                                               surf_def_v(1)%rad_lw_res,       &
                                               surf_lsm_v(1)%rad_lw_res,       &
                                               surf_usm_v(1)%rad_lw_res,       &
                                               surf_def_v(2)%rad_lw_res,       &
                                               surf_lsm_v(2)%rad_lw_res,       &
                                               surf_usm_v(2)%rad_lw_res,       &
                                               surf_def_v(3)%rad_lw_res,       &
                                               surf_lsm_v(3)%rad_lw_res,       &
                                               surf_usm_v(3)%rad_lw_res )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF

            CASE ( 'uvw1' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%uvw_abs,     &
                                               surf_def_h(1)%uvw_abs,          &
                                               surf_lsm_h%uvw_abs,             &
                                               surf_usm_h%uvw_abs,             &
                                               surf_def_v(0)%uvw_abs,          &
                                               surf_lsm_v(0)%uvw_abs,          &
                                               surf_usm_v(0)%uvw_abs,          &
                                               surf_def_v(1)%uvw_abs,          &
                                               surf_lsm_v(1)%uvw_abs,          &
                                               surf_usm_v(1)%uvw_abs,          &
                                               surf_def_v(2)%uvw_abs,          &
                                               surf_lsm_v(2)%uvw_abs,          &
                                               surf_usm_v(2)%uvw_abs,          &
                                               surf_def_v(3)%uvw_abs,          &
                                               surf_lsm_v(3)%uvw_abs,          &
                                               surf_usm_v(3)%uvw_abs )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF
!
!--         Waste heat from indoor model
            CASE ( 'waste_heat' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%waste_heat,  &
                                               surf_def_h(1)%waste_heat,       &
                                               surf_lsm_h%waste_heat,          &
                                               surf_usm_h%waste_heat,          &
                                               surf_def_v(0)%waste_heat,       &
                                               surf_lsm_v(0)%waste_heat,       &
                                               surf_usm_v(0)%waste_heat,       &
                                               surf_def_v(1)%waste_heat,       &
                                               surf_lsm_v(1)%waste_heat,       &
                                               surf_usm_v(1)%waste_heat,       &
                                               surf_def_v(2)%waste_heat,       &
                                               surf_lsm_v(2)%waste_heat,       &
                                               surf_usm_v(2)%waste_heat,       &
                                               surf_def_v(3)%waste_heat,       &
                                               surf_lsm_v(3)%waste_heat,       &
                                               surf_usm_v(3)%waste_heat )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF
!
!--         Innermost building wall flux from indoor model
            CASE ( 'im_hf' )
!
!--            Output of instantaneous data
               IF ( av == 0 )  THEN
                  CALL surface_data_output_collect( surf_def_h(0)%iwghf_eb,    &
                                               surf_def_h(1)%iwghf_eb,         &
                                               surf_lsm_h%iwghf_eb,            &
                                               surf_usm_h%iwghf_eb,            &
                                               surf_def_v(0)%iwghf_eb,         &
                                               surf_lsm_v(0)%iwghf_eb,         &
                                               surf_usm_v(0)%iwghf_eb,         &
                                               surf_def_v(1)%iwghf_eb,         &
                                               surf_lsm_v(1)%iwghf_eb,         &
                                               surf_usm_v(1)%iwghf_eb,         &
                                               surf_def_v(2)%iwghf_eb,         &
                                               surf_lsm_v(2)%iwghf_eb,         &
                                               surf_usm_v(2)%iwghf_eb,         &
                                               surf_def_v(3)%iwghf_eb,         &
                                               surf_lsm_v(3)%iwghf_eb,         &
                                               surf_usm_v(3)%iwghf_eb )
               ELSE
!
!--               Output of averaged data
                  surfaces%var_out(:) = surfaces%var_av(:,n_out) /             &
                                        REAL( average_count_surf, KIND=wp )
                  surfaces%var_av(:,n_out) = 0.0_wp

               ENDIF
!
!--            Add further variables:
!--            'css', 'cssws', 'qsws_liq', 'qsws_soil', 'qsws_veg'

         END SELECT
!
!--      Write to binary file:
!--      - surfaces%points ( 3, 1-npoints )
!--      - surfaces%polygons ( 5, 1-ns )
!--      - surfaces%var_out ( 1-ns, time )
!--      - Dimension: 1-nsurfaces, 1-npoints - can be ordered consecutively
!--      - Distinguish between average and non-average data
         IF ( to_vtk )  THEN
            DO  i = 0, io_blocks-1
               IF ( i == io_group )  THEN
                  WRITE ( 25+av )  LEN_TRIM( 'time' )
                  WRITE ( 25+av )  'time'
                  WRITE ( 25+av )  time_since_reference_point
                  WRITE ( 25+av )  LEN_TRIM( trimvar )
                  WRITE ( 25+av )  TRIM( trimvar )
                  WRITE ( 25+av )  surfaces%var_out
               ENDIF
#if defined( __parallel )
               CALL MPI_BARRIER( comm2d, ierr )
#endif
            ENDDO
         ENDIF

         IF ( to_netcdf )  THEN
#if defined( __netcdf4_parallel )
!
!--         Write output array to file
            nc_stat = NF90_PUT_VAR( id_set_surf(av), id_var_dosurf(av,n_out),  &
                                    surfaces%var_out,                          &
                                    start = (/ surfaces%s(1),                  &
                                               dosurf_time_count(av) /),       &
                                    count = (/ surfaces%ns, 1 /) )
            CALL netcdf_handle_error( 'surface_data_output', 6667 )
#endif
         ENDIF

      ENDDO

!
!--   If averaged output was written to NetCDF file, set the counter to zero
      IF ( av == 1 )  average_count_surf = 0

   END SUBROUTINE surface_data_output

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Routine for controlling the data averaging.
!------------------------------------------------------------------------------!
   SUBROUTINE surface_data_output_averaging

      IMPLICIT NONE

      CHARACTER(LEN=100) ::  trimvar !< dummy variable for current output variable

      INTEGER(iwp) ::  n_out  !< counter variables for surface output

      n_out = 0
      DO  WHILE ( dosurf(1,n_out+1)(1:1) /= ' ' )

         n_out   = n_out + 1
         trimvar = TRIM( dosurf(1,n_out) )

         SELECT CASE ( trimvar )

            CASE ( 'us' )
               CALL surface_data_output_sum_up( surf_def_h(0)%us,              &
                                           surf_def_h(1)%us,                   &
                                           surf_lsm_h%us,                      &
                                           surf_usm_h%us,                      &
                                           surf_def_v(0)%us,                   &
                                           surf_lsm_v(0)%us,                   &
                                           surf_usm_v(0)%us,                   &
                                           surf_def_v(1)%us,                   &
                                           surf_lsm_v(1)%us,                   &
                                           surf_usm_v(1)%us,                   &
                                           surf_def_v(2)%us,                   &
                                           surf_lsm_v(2)%us,                   &
                                           surf_usm_v(2)%us,                   &
                                           surf_def_v(3)%us,                   &
                                           surf_lsm_v(3)%us,                   &
                                           surf_usm_v(3)%us, n_out )

            CASE ( 'ts' )
               CALL surface_data_output_sum_up( surf_def_h(0)%ts,              &
                                           surf_def_h(1)%ts,                   &
                                           surf_lsm_h%ts,                      &
                                           surf_usm_h%ts,                      &
                                           surf_def_v(0)%ts,                   &
                                           surf_lsm_v(0)%ts,                   &
                                           surf_usm_v(0)%ts,                   &
                                           surf_def_v(1)%ts,                   &
                                           surf_lsm_v(1)%ts,                   &
                                           surf_usm_v(1)%ts,                   &
                                           surf_def_v(2)%ts,                   &
                                           surf_lsm_v(2)%ts,                   &
                                           surf_usm_v(2)%ts,                   &
                                           surf_def_v(3)%ts,                   &
                                           surf_lsm_v(3)%ts,                   &
                                           surf_usm_v(3)%ts, n_out )

            CASE ( 'qs' )
               CALL surface_data_output_sum_up( surf_def_h(0)%qs,              &
                                           surf_def_h(1)%qs,                   &
                                           surf_lsm_h%qs,                      &
                                           surf_usm_h%qs,                      &
                                           surf_def_v(0)%qs,                   &
                                           surf_lsm_v(0)%qs,                   &
                                           surf_usm_v(0)%qs,                   &
                                           surf_def_v(1)%qs,                   &
                                           surf_lsm_v(1)%qs,                   &
                                           surf_usm_v(1)%qs,                   &
                                           surf_def_v(2)%qs,                   &
                                           surf_lsm_v(2)%qs,                   &
                                           surf_usm_v(2)%qs,                   &
                                           surf_def_v(3)%qs,                   &
                                           surf_lsm_v(3)%qs,                   &
                                           surf_usm_v(3)%qs, n_out )

            CASE ( 'ss' )
               CALL surface_data_output_sum_up( surf_def_h(0)%ss,              &
                                           surf_def_h(1)%ss,                   &
                                           surf_lsm_h%ss,                      &
                                           surf_usm_h%ss,                      &
                                           surf_def_v(0)%ss,                   &
                                           surf_lsm_v(0)%ss,                   &
                                           surf_usm_v(0)%ss,                   &
                                           surf_def_v(1)%ss,                   &
                                           surf_lsm_v(1)%ss,                   &
                                           surf_usm_v(1)%ss,                   &
                                           surf_def_v(2)%ss,                   &
                                           surf_lsm_v(2)%ss,                   &
                                           surf_usm_v(2)%ss,                   &
                                           surf_def_v(3)%ss,                   &
                                           surf_lsm_v(3)%ss,                   &
                                           surf_usm_v(3)%ss, n_out )

            CASE ( 'qcs' )
               CALL surface_data_output_sum_up( surf_def_h(0)%qcs,             &
                                           surf_def_h(1)%qcs,                  &
                                           surf_lsm_h%qcs,                     &
                                           surf_usm_h%qcs,                     &
                                           surf_def_v(0)%qcs,                  &
                                           surf_lsm_v(0)%qcs,                  &
                                           surf_usm_v(0)%qcs,                  &
                                           surf_def_v(1)%qcs,                  &
                                           surf_lsm_v(1)%qcs,                  &
                                           surf_usm_v(1)%qcs,                  &
                                           surf_def_v(2)%qcs,                  &
                                           surf_lsm_v(2)%qcs,                  &
                                           surf_usm_v(2)%qcs,                  &
                                           surf_def_v(3)%qcs,                  &
                                           surf_lsm_v(3)%qcs,                  &
                                           surf_usm_v(3)%qcs, n_out )

            CASE ( 'ncs' )
               CALL surface_data_output_sum_up( surf_def_h(0)%ncs,             &
                                           surf_def_h(1)%ncs,                  &
                                           surf_lsm_h%ncs,                     &
                                           surf_usm_h%ncs,                     &
                                           surf_def_v(0)%ncs,                  &
                                           surf_lsm_v(0)%ncs,                  &
                                           surf_usm_v(0)%ncs,                  &
                                           surf_def_v(1)%ncs,                  &
                                           surf_lsm_v(1)%ncs,                  &
                                           surf_usm_v(1)%ncs,                  &
                                           surf_def_v(2)%ncs,                  &
                                           surf_lsm_v(2)%ncs,                  &
                                           surf_usm_v(2)%ncs,                  &
                                           surf_def_v(3)%ncs,                  &
                                           surf_lsm_v(3)%ncs,                  &
                                           surf_usm_v(3)%ncs, n_out )

            CASE ( 'qrs' )
               CALL surface_data_output_sum_up( surf_def_h(0)%qrs,             &
                                           surf_def_h(1)%qrs,                  &
                                           surf_lsm_h%qrs,                     &
                                           surf_usm_h%qrs,                     &
                                           surf_def_v(0)%qrs,                  &
                                           surf_lsm_v(0)%qrs,                  &
                                           surf_usm_v(0)%qrs,                  &
                                           surf_def_v(1)%qrs,                  &
                                           surf_lsm_v(1)%qrs,                  &
                                           surf_usm_v(1)%qrs,                  &
                                           surf_def_v(2)%qrs,                  &
                                           surf_lsm_v(2)%qrs,                  &
                                           surf_usm_v(2)%qrs,                  &
                                           surf_def_v(3)%qrs,                  &
                                           surf_lsm_v(3)%qrs,                  &
                                           surf_usm_v(3)%qrs, n_out )

            CASE ( 'nrs' )
               CALL surface_data_output_sum_up( surf_def_h(0)%nrs,             &
                                           surf_def_h(1)%nrs,                  &
                                           surf_lsm_h%nrs,                     &
                                           surf_usm_h%nrs,                     &
                                           surf_def_v(0)%nrs,                  &
                                           surf_lsm_v(0)%nrs,                  &
                                           surf_usm_v(0)%nrs,                  &
                                           surf_def_v(1)%nrs,                  &
                                           surf_lsm_v(1)%nrs,                  &
                                           surf_usm_v(1)%nrs,                  &
                                           surf_def_v(2)%nrs,                  &
                                           surf_lsm_v(2)%nrs,                  &
                                           surf_usm_v(2)%nrs,                  &
                                           surf_def_v(3)%nrs,                  &
                                           surf_lsm_v(3)%nrs,                  &
                                           surf_usm_v(3)%nrs, n_out )

            CASE ( 'ol' )
               CALL surface_data_output_sum_up( surf_def_h(0)%ol,              &
                                           surf_def_h(1)%ol,                   &
                                           surf_lsm_h%ol,                      &
                                           surf_usm_h%ol,                      &
                                           surf_def_v(0)%ol,                   &
                                           surf_lsm_v(0)%ol,                   &
                                           surf_usm_v(0)%ol,                   &
                                           surf_def_v(1)%ol,                   &
                                           surf_lsm_v(1)%ol,                   &
                                           surf_usm_v(1)%ol,                   &
                                           surf_def_v(2)%ol,                   &
                                           surf_lsm_v(2)%ol,                   &
                                           surf_usm_v(2)%ol,                   &
                                           surf_def_v(3)%ol,                   &
                                           surf_lsm_v(3)%ol,                   &
                                           surf_usm_v(3)%ol, n_out )

            CASE ( 'z0' )
               CALL surface_data_output_sum_up( surf_def_h(0)%z0,              &
                                           surf_def_h(1)%z0,                   &
                                           surf_lsm_h%z0,                      &
                                           surf_usm_h%z0,                      &
                                           surf_def_v(0)%z0,                   &
                                           surf_lsm_v(0)%z0,                   &
                                           surf_usm_v(0)%z0,                   &
                                           surf_def_v(1)%z0,                   &
                                           surf_lsm_v(1)%z0,                   &
                                           surf_usm_v(1)%z0,                   &
                                           surf_def_v(2)%z0,                   &
                                           surf_lsm_v(2)%z0,                   &
                                           surf_usm_v(2)%z0,                   &
                                           surf_def_v(3)%z0,                   &
                                           surf_lsm_v(3)%z0,                   &
                                           surf_usm_v(3)%z0, n_out )

            CASE ( 'z0h' )
               CALL surface_data_output_sum_up( surf_def_h(0)%z0h,             &
                                           surf_def_h(1)%z0h,                  &
                                           surf_lsm_h%z0h,                     &
                                           surf_usm_h%z0h,                     &
                                           surf_def_v(0)%z0h,                  &
                                           surf_lsm_v(0)%z0h,                  &
                                           surf_usm_v(0)%z0h,                  &
                                           surf_def_v(1)%z0h,                  &
                                           surf_lsm_v(1)%z0h,                  &
                                           surf_usm_v(1)%z0h,                  &
                                           surf_def_v(2)%z0h,                  &
                                           surf_lsm_v(2)%z0h,                  &
                                           surf_usm_v(2)%z0h,                  &
                                           surf_def_v(3)%z0h,                  &
                                           surf_lsm_v(3)%z0h,                  &
                                           surf_usm_v(3)%z0h, n_out )

            CASE ( 'z0q' )
               CALL surface_data_output_sum_up( surf_def_h(0)%z0q,             &
                                           surf_def_h(1)%z0q,                  &
                                           surf_lsm_h%z0q,                     &
                                           surf_usm_h%z0q,                     &
                                           surf_def_v(0)%z0q,                  &
                                           surf_lsm_v(0)%z0q,                  &
                                           surf_usm_v(0)%z0q,                  &
                                           surf_def_v(1)%z0q,                  &
                                           surf_lsm_v(1)%z0q,                  &
                                           surf_usm_v(1)%z0q,                  &
                                           surf_def_v(2)%z0q,                  &
                                           surf_lsm_v(2)%z0q,                  &
                                           surf_usm_v(2)%z0q,                  &
                                           surf_def_v(3)%z0q,                  &
                                           surf_lsm_v(3)%z0q,                  &
                                           surf_usm_v(3)%z0q, n_out )

            CASE ( 'theta1' )
               CALL surface_data_output_sum_up( surf_def_h(0)%pt1,             &
                                           surf_def_h(1)%pt1,                  &
                                           surf_lsm_h%pt1,                     &
                                           surf_usm_h%pt1,                     &
                                           surf_def_v(0)%pt1,                  &
                                           surf_lsm_v(0)%pt1,                  &
                                           surf_usm_v(0)%pt1,                  &
                                           surf_def_v(1)%pt1,                  &
                                           surf_lsm_v(1)%pt1,                  &
                                           surf_usm_v(1)%pt1,                  &
                                           surf_def_v(2)%pt1,                  &
                                           surf_lsm_v(2)%pt1,                  &
                                           surf_usm_v(2)%pt1,                  &
                                           surf_def_v(3)%pt1,                  &
                                           surf_lsm_v(3)%pt1,                  &
                                           surf_usm_v(3)%pt1, n_out )

            CASE ( 'qv1' )
               CALL surface_data_output_sum_up( surf_def_h(0)%qv1,             &
                                           surf_def_h(1)%qv1,                  &
                                           surf_lsm_h%qv1,                     &
                                           surf_usm_h%qv1,                     &
                                           surf_def_v(0)%qv1,                  &
                                           surf_lsm_v(0)%qv1,                  &
                                           surf_usm_v(0)%qv1,                  &
                                           surf_def_v(1)%qv1,                  &
                                           surf_lsm_v(1)%qv1,                  &
                                           surf_usm_v(1)%qv1,                  &
                                           surf_def_v(2)%qv1,                  &
                                           surf_lsm_v(2)%qv1,                  &
                                           surf_usm_v(2)%qv1,                  &
                                           surf_def_v(3)%qv1,                  &
                                           surf_lsm_v(3)%qv1,                  &
                                           surf_usm_v(3)%qv1, n_out )

            CASE ( 'thetav1' )
               CALL surface_data_output_sum_up( surf_def_h(0)%vpt1,            &
                                           surf_def_h(1)%vpt1,                 &
                                           surf_lsm_h%vpt1,                    &
                                           surf_usm_h%vpt1,                    &
                                           surf_def_v(0)%vpt1,                 &
                                           surf_lsm_v(0)%vpt1,                 &
                                           surf_usm_v(0)%vpt1,                 &
                                           surf_def_v(1)%vpt1,                 &
                                           surf_lsm_v(1)%vpt1,                 &
                                           surf_usm_v(1)%vpt1,                 &
                                           surf_def_v(2)%vpt1,                 &
                                           surf_lsm_v(2)%vpt1,                 &
                                           surf_usm_v(2)%vpt1,                 &
                                           surf_def_v(3)%vpt1,                 &
                                           surf_lsm_v(3)%vpt1,                 &
                                           surf_usm_v(3)%vpt1, n_out )

            CASE ( 'usws' )
               CALL surface_data_output_sum_up( surf_def_h(0)%usws,            &
                                           surf_def_h(1)%usws,                 &
                                           surf_lsm_h%usws,                    &
                                           surf_usm_h%usws,                    &
                                           surf_def_v(0)%usws,                 &
                                           surf_lsm_v(0)%usws,                 &
                                           surf_usm_v(0)%usws,                 &
                                           surf_def_v(1)%usws,                 &
                                           surf_lsm_v(1)%usws,                 &
                                           surf_usm_v(1)%usws,                 &
                                           surf_def_v(2)%usws,                 &
                                           surf_lsm_v(2)%usws,                 &
                                           surf_usm_v(2)%usws,                 &
                                           surf_def_v(3)%usws,                 &
                                           surf_lsm_v(3)%usws,                 &
                                           surf_usm_v(3)%usws, n_out )

            CASE ( 'vsws' )
               CALL surface_data_output_sum_up( surf_def_h(0)%vsws,            &
                                           surf_def_h(1)%vsws,                 &
                                           surf_lsm_h%vsws,                    &
                                           surf_usm_h%vsws,                    &
                                           surf_def_v(0)%vsws,                 &
                                           surf_lsm_v(0)%vsws,                 &
                                           surf_usm_v(0)%vsws,                 &
                                           surf_def_v(1)%vsws,                 &
                                           surf_lsm_v(1)%vsws,                 &
                                           surf_usm_v(1)%vsws,                 &
                                           surf_def_v(2)%vsws,                 &
                                           surf_lsm_v(2)%vsws,                 &
                                           surf_usm_v(2)%vsws,                 &
                                           surf_def_v(3)%vsws,                 &
                                           surf_lsm_v(3)%vsws,                 &
                                           surf_usm_v(3)%vsws, n_out )

            CASE ( 'shf' )
               CALL surface_data_output_sum_up( surf_def_h(0)%shf,             &
                                           surf_def_h(1)%shf,                  &
                                           surf_lsm_h%shf,                     &
                                           surf_usm_h%shf,                     &
                                           surf_def_v(0)%shf,                  &
                                           surf_lsm_v(0)%shf,                  &
                                           surf_usm_v(0)%shf,                  &
                                           surf_def_v(1)%shf,                  &
                                           surf_lsm_v(1)%shf,                  &
                                           surf_usm_v(1)%shf,                  &
                                           surf_def_v(2)%shf,                  &
                                           surf_lsm_v(2)%shf,                  &
                                           surf_usm_v(2)%shf,                  &
                                           surf_def_v(3)%shf,                  &
                                           surf_lsm_v(3)%shf,                  &
                                           surf_usm_v(3)%shf, n_out )

            CASE ( 'qsws' )
               CALL surface_data_output_sum_up( surf_def_h(0)%qsws,            &
                                           surf_def_h(1)%qsws,                 &
                                           surf_lsm_h%qsws,                    &
                                           surf_usm_h%qsws,                    &
                                           surf_def_v(0)%qsws,                 &
                                           surf_lsm_v(0)%qsws,                 &
                                           surf_usm_v(0)%qsws,                 &
                                           surf_def_v(1)%qsws,                 &
                                           surf_lsm_v(1)%qsws,                 &
                                           surf_usm_v(1)%qsws,                 &
                                           surf_def_v(2)%qsws,                 &
                                           surf_lsm_v(2)%qsws,                 &
                                           surf_usm_v(2)%qsws,                 &
                                           surf_def_v(3)%qsws,                 &
                                           surf_lsm_v(3)%qsws,                 &
                                           surf_usm_v(3)%qsws, n_out )

            CASE ( 'ssws' )
               CALL surface_data_output_sum_up( surf_def_h(0)%ssws,            &
                                           surf_def_h(1)%ssws,                 &
                                           surf_lsm_h%ssws,                    &
                                           surf_usm_h%ssws,                    &
                                           surf_def_v(0)%ssws,                 &
                                           surf_lsm_v(0)%ssws,                 &
                                           surf_usm_v(0)%ssws,                 &
                                           surf_def_v(1)%ssws,                 &
                                           surf_lsm_v(1)%ssws,                 &
                                           surf_usm_v(1)%ssws,                 &
                                           surf_def_v(2)%ssws,                 &
                                           surf_lsm_v(2)%ssws,                 &
                                           surf_usm_v(2)%ssws,                 &
                                           surf_def_v(3)%ssws,                 &
                                           surf_lsm_v(3)%ssws,                 &
                                           surf_usm_v(3)%ssws, n_out )

            CASE ( 'qcsws' )
               CALL surface_data_output_sum_up( surf_def_h(0)%qcsws,           &
                                           surf_def_h(1)%qcsws,                &
                                           surf_lsm_h%qcsws,                   &
                                           surf_usm_h%qcsws,                   &
                                           surf_def_v(0)%qcsws,                &
                                           surf_lsm_v(0)%qcsws,                &
                                           surf_usm_v(0)%qcsws,                &
                                           surf_def_v(1)%qcsws,                &
                                           surf_lsm_v(1)%qcsws,                &
                                           surf_usm_v(1)%qcsws,                &
                                           surf_def_v(2)%qcsws,                &
                                           surf_lsm_v(2)%qcsws,                &
                                           surf_usm_v(2)%qcsws,                &
                                           surf_def_v(3)%qcsws,                &
                                           surf_lsm_v(3)%qcsws,                &
                                           surf_usm_v(3)%qcsws, n_out )

            CASE ( 'ncsws' )
               CALL surface_data_output_sum_up( surf_def_h(0)%ncsws,           &
                                           surf_def_h(1)%ncsws,                &
                                           surf_lsm_h%ncsws,                   &
                                           surf_usm_h%ncsws,                   &
                                           surf_def_v(0)%ncsws,                &
                                           surf_lsm_v(0)%ncsws,                &
                                           surf_usm_v(0)%ncsws,                &
                                           surf_def_v(1)%ncsws,                &
                                           surf_lsm_v(1)%ncsws,                &
                                           surf_usm_v(1)%ncsws,                &
                                           surf_def_v(2)%ncsws,                &
                                           surf_lsm_v(2)%ncsws,                &
                                           surf_usm_v(2)%ncsws,                &
                                           surf_def_v(3)%ncsws,                &
                                           surf_lsm_v(3)%ncsws,                &
                                           surf_usm_v(3)%ncsws, n_out )

            CASE ( 'qrsws' )
               CALL surface_data_output_sum_up( surf_def_h(0)%qrsws,           &
                                           surf_def_h(1)%qrsws,                &
                                           surf_lsm_h%qrsws,                   &
                                           surf_usm_h%qrsws,                   &
                                           surf_def_v(0)%qrsws,                &
                                           surf_lsm_v(0)%qrsws,                &
                                           surf_usm_v(0)%qrsws,                &
                                           surf_def_v(1)%qrsws,                &
                                           surf_lsm_v(1)%qrsws,                &
                                           surf_usm_v(1)%qrsws,                &
                                           surf_def_v(2)%qrsws,                &
                                           surf_lsm_v(2)%qrsws,                &
                                           surf_usm_v(2)%qrsws,                &
                                           surf_def_v(3)%qrsws,                &
                                           surf_lsm_v(3)%qrsws,                &
                                           surf_usm_v(3)%qrsws, n_out )

            CASE ( 'nrsws' )
               CALL surface_data_output_sum_up( surf_def_h(0)%nrsws,           &
                                           surf_def_h(1)%nrsws,                &
                                           surf_lsm_h%nrsws,                   &
                                           surf_usm_h%nrsws,                   &
                                           surf_def_v(0)%nrsws,                &
                                           surf_lsm_v(0)%nrsws,                &
                                           surf_usm_v(0)%nrsws,                &
                                           surf_def_v(1)%nrsws,                &
                                           surf_lsm_v(1)%nrsws,                &
                                           surf_usm_v(1)%nrsws,                &
                                           surf_def_v(2)%nrsws,                &
                                           surf_lsm_v(2)%nrsws,                &
                                           surf_usm_v(2)%nrsws,                &
                                           surf_def_v(3)%nrsws,                &
                                           surf_lsm_v(3)%nrsws,                &
                                           surf_usm_v(3)%nrsws, n_out )

            CASE ( 'sasws' )
               CALL surface_data_output_sum_up( surf_def_h(0)%sasws,           &
                                           surf_def_h(1)%sasws,                &
                                           surf_lsm_h%sasws,                   &
                                           surf_usm_h%sasws,                   &
                                           surf_def_v(0)%sasws,                &
                                           surf_lsm_v(0)%sasws,                &
                                           surf_usm_v(0)%sasws,                &
                                           surf_def_v(1)%sasws,                &
                                           surf_lsm_v(1)%sasws,                &
                                           surf_usm_v(1)%sasws,                &
                                           surf_def_v(2)%sasws,                &
                                           surf_lsm_v(2)%sasws,                &
                                           surf_usm_v(2)%sasws,                &
                                           surf_def_v(3)%sasws,                &
                                           surf_lsm_v(3)%sasws,                &
                                           surf_usm_v(3)%sasws, n_out )

            CASE ( 'q_surface' )
               CALL surface_data_output_sum_up( surf_def_h(0)%q_surface,       &
                                           surf_def_h(1)%q_surface,            &
                                           surf_lsm_h%q_surface,               &
                                           surf_usm_h%q_surface,               &
                                           surf_def_v(0)%q_surface,            &
                                           surf_lsm_v(0)%q_surface,            &
                                           surf_usm_v(0)%q_surface,            &
                                           surf_def_v(1)%q_surface,            &
                                           surf_lsm_v(1)%q_surface,            &
                                           surf_usm_v(1)%q_surface,            &
                                           surf_def_v(2)%q_surface,            &
                                           surf_lsm_v(2)%q_surface,            &
                                           surf_usm_v(2)%q_surface,            &
                                           surf_def_v(3)%q_surface,            &
                                           surf_lsm_v(3)%q_surface,            &
                                           surf_usm_v(3)%q_surface, n_out )


            CASE ( 'theta_surface' )
               CALL surface_data_output_sum_up( surf_def_h(0)%pt_surface,      &
                                           surf_def_h(1)%pt_surface,           &
                                           surf_lsm_h%pt_surface,              &
                                           surf_usm_h%pt_surface,              &
                                           surf_def_v(0)%pt_surface,           &
                                           surf_lsm_v(0)%pt_surface,           &
                                           surf_usm_v(0)%pt_surface,           &
                                           surf_def_v(1)%pt_surface,           &
                                           surf_lsm_v(1)%pt_surface,           &
                                           surf_usm_v(1)%pt_surface,           &
                                           surf_def_v(2)%pt_surface,           &
                                           surf_lsm_v(2)%pt_surface,           &
                                           surf_usm_v(2)%pt_surface,           &
                                           surf_def_v(3)%pt_surface,           &
                                           surf_lsm_v(3)%pt_surface,           &
                                           surf_usm_v(3)%pt_surface, n_out )

            CASE ( 'thetav_surface' )
               CALL surface_data_output_sum_up( surf_def_h(0)%vpt_surface,     &
                                           surf_def_h(1)%vpt_surface,          &
                                           surf_lsm_h%vpt_surface,             &
                                           surf_usm_h%vpt_surface,             &
                                           surf_def_v(0)%vpt_surface,          &
                                           surf_lsm_v(0)%vpt_surface,          &
                                           surf_usm_v(0)%vpt_surface,          &
                                           surf_def_v(1)%vpt_surface,          &
                                           surf_lsm_v(1)%vpt_surface,          &
                                           surf_usm_v(1)%vpt_surface,          &
                                           surf_def_v(2)%vpt_surface,          &
                                           surf_lsm_v(2)%vpt_surface,          &
                                           surf_usm_v(2)%vpt_surface,          &
                                           surf_def_v(3)%vpt_surface,          &
                                           surf_lsm_v(3)%vpt_surface,          &
                                           surf_usm_v(3)%vpt_surface, n_out )

            CASE ( 'rad_net' )
               CALL surface_data_output_sum_up( surf_def_h(0)%rad_net,         &
                                           surf_def_h(1)%rad_net,              &
                                           surf_lsm_h%rad_net,                 &
                                           surf_usm_h%rad_net,                 &
                                           surf_def_v(0)%rad_net,              &
                                           surf_lsm_v(0)%rad_net,              &
                                           surf_usm_v(0)%rad_net,              &
                                           surf_def_v(1)%rad_net,              &
                                           surf_lsm_v(1)%rad_net,              &
                                           surf_usm_v(1)%rad_net,              &
                                           surf_def_v(2)%rad_net,              &
                                           surf_lsm_v(2)%rad_net,              &
                                           surf_usm_v(2)%rad_net,              &
                                           surf_def_v(3)%rad_net,              &
                                           surf_lsm_v(3)%rad_net,              &
                                           surf_usm_v(3)%rad_net, n_out )

            CASE ( 'rad_lw_in' )
               CALL surface_data_output_sum_up( surf_def_h(0)%rad_lw_in,       &
                                           surf_def_h(1)%rad_lw_in,            &
                                           surf_lsm_h%rad_lw_in,               &
                                           surf_usm_h%rad_lw_in,               &
                                           surf_def_v(0)%rad_lw_in,            &
                                           surf_lsm_v(0)%rad_lw_in,            &
                                           surf_usm_v(0)%rad_lw_in,            &
                                           surf_def_v(1)%rad_lw_in,            &
                                           surf_lsm_v(1)%rad_lw_in,            &
                                           surf_usm_v(1)%rad_lw_in,            &
                                           surf_def_v(2)%rad_lw_in,            &
                                           surf_lsm_v(2)%rad_lw_in,            &
                                           surf_usm_v(2)%rad_lw_in,            &
                                           surf_def_v(3)%rad_lw_in,            &
                                           surf_lsm_v(3)%rad_lw_in,            &
                                           surf_usm_v(3)%rad_lw_in, n_out )

            CASE ( 'rad_lw_out' )
               CALL surface_data_output_sum_up( surf_def_h(0)%rad_lw_out,      &
                                           surf_def_h(1)%rad_lw_out,           &
                                           surf_lsm_h%rad_lw_out,              &
                                           surf_usm_h%rad_lw_out,              &
                                           surf_def_v(0)%rad_lw_out,           &
                                           surf_lsm_v(0)%rad_lw_out,           &
                                           surf_usm_v(0)%rad_lw_out,           &
                                           surf_def_v(1)%rad_lw_out,           &
                                           surf_lsm_v(1)%rad_lw_out,           &
                                           surf_usm_v(1)%rad_lw_out,           &
                                           surf_def_v(2)%rad_lw_out,           &
                                           surf_lsm_v(2)%rad_lw_out,           &
                                           surf_usm_v(2)%rad_lw_out,           &
                                           surf_def_v(3)%rad_lw_out,           &
                                           surf_lsm_v(3)%rad_lw_out,           &
                                           surf_usm_v(3)%rad_lw_out, n_out )

            CASE ( 'rad_sw_in' )
               CALL surface_data_output_sum_up( surf_def_h(0)%rad_sw_in,       &
                                           surf_def_h(1)%rad_sw_in,            &
                                           surf_lsm_h%rad_sw_in,               &
                                           surf_usm_h%rad_sw_in,               &
                                           surf_def_v(0)%rad_sw_in,            &
                                           surf_lsm_v(0)%rad_sw_in,            &
                                           surf_usm_v(0)%rad_sw_in,            &
                                           surf_def_v(1)%rad_sw_in,            &
                                           surf_lsm_v(1)%rad_sw_in,            &
                                           surf_usm_v(1)%rad_sw_in,            &
                                           surf_def_v(2)%rad_sw_in,            &
                                           surf_lsm_v(2)%rad_sw_in,            &
                                           surf_usm_v(2)%rad_sw_in,            &
                                           surf_def_v(3)%rad_sw_in,            &
                                           surf_lsm_v(3)%rad_sw_in,            &
                                           surf_usm_v(3)%rad_sw_in, n_out )

            CASE ( 'rad_sw_out' )
               CALL surface_data_output_sum_up( surf_def_h(0)%rad_sw_out,      &
                                           surf_def_h(1)%rad_sw_out,           &
                                           surf_lsm_h%rad_sw_out,              &
                                           surf_usm_h%rad_sw_out,              &
                                           surf_def_v(0)%rad_sw_out,           &
                                           surf_lsm_v(0)%rad_sw_out,           &
                                           surf_usm_v(0)%rad_sw_out,           &
                                           surf_def_v(1)%rad_sw_out,           &
                                           surf_lsm_v(1)%rad_sw_out,           &
                                           surf_usm_v(1)%rad_sw_out,           &
                                           surf_def_v(2)%rad_sw_out,           &
                                           surf_lsm_v(2)%rad_sw_out,           &
                                           surf_usm_v(2)%rad_sw_out,           &
                                           surf_def_v(3)%rad_sw_out,           &
                                           surf_lsm_v(3)%rad_sw_out,           &
                                           surf_usm_v(3)%rad_sw_out, n_out )

            CASE ( 'ghf' )
               CALL surface_data_output_sum_up( surf_def_h(0)%ghf,             &
                                           surf_def_h(1)%ghf,                  &
                                           surf_lsm_h%ghf,                     &
                                           surf_usm_h%ghf,                     &
                                           surf_def_v(0)%ghf,                  &
                                           surf_lsm_v(0)%ghf,                  &
                                           surf_usm_v(0)%ghf,                  &
                                           surf_def_v(1)%ghf,                  &
                                           surf_lsm_v(1)%ghf,                  &
                                           surf_usm_v(1)%ghf,                  &
                                           surf_def_v(2)%ghf,                  &
                                           surf_lsm_v(2)%ghf,                  &
                                           surf_usm_v(2)%ghf,                  &
                                           surf_def_v(3)%ghf,                  &
                                           surf_lsm_v(3)%ghf,                  &
                                           surf_usm_v(3)%ghf, n_out )

            CASE ( 'r_a' )
               CALL surface_data_output_sum_up( surf_def_h(0)%r_a,             &
                                           surf_def_h(1)%r_a,                  &
                                           surf_lsm_h%r_a,                     &
                                           surf_usm_h%r_a,                     &
                                           surf_def_v(0)%r_a,                  &
                                           surf_lsm_v(0)%r_a,                  &
                                           surf_usm_v(0)%r_a,                  &
                                           surf_def_v(1)%r_a,                  &
                                           surf_lsm_v(1)%r_a,                  &
                                           surf_usm_v(1)%r_a,                  &
                                           surf_def_v(2)%r_a,                  &
                                           surf_lsm_v(2)%r_a,                  &
                                           surf_usm_v(2)%r_a,                  &
                                           surf_def_v(3)%r_a,                  &
                                           surf_lsm_v(3)%r_a,                  &
                                           surf_usm_v(3)%r_a, n_out )

            CASE ( 'r_soil' )
               CALL surface_data_output_sum_up( surf_def_h(0)%r_soil,          &
                                           surf_def_h(1)%r_soil,               &
                                           surf_lsm_h%r_soil,                  &
                                           surf_usm_h%r_soil,                  &
                                           surf_def_v(0)%r_soil,               &
                                           surf_lsm_v(0)%r_soil,               &
                                           surf_usm_v(0)%r_soil,               &
                                           surf_def_v(1)%r_soil,               &
                                           surf_lsm_v(1)%r_soil,               &
                                           surf_usm_v(1)%r_soil,               &
                                           surf_def_v(2)%r_soil,               &
                                           surf_lsm_v(2)%r_soil,               &
                                           surf_usm_v(2)%r_soil,               &
                                           surf_def_v(3)%r_soil,               &
                                           surf_lsm_v(3)%r_soil,               &
                                           surf_usm_v(3)%r_soil, n_out )

            CASE ( 'r_canopy' )
               CALL surface_data_output_sum_up( surf_def_h(0)%r_canopy,        &
                                           surf_def_h(1)%r_canopy,             &
                                           surf_lsm_h%r_canopy,                &
                                           surf_usm_h%r_canopy,                &
                                           surf_def_v(0)%r_canopy,             &
                                           surf_lsm_v(0)%r_canopy,             &
                                           surf_usm_v(0)%r_canopy,             &
                                           surf_def_v(1)%r_canopy,             &
                                           surf_lsm_v(1)%r_canopy,             &
                                           surf_usm_v(1)%r_canopy,             &
                                           surf_def_v(2)%r_canopy,             &
                                           surf_lsm_v(2)%r_canopy,             &
                                           surf_usm_v(2)%r_canopy,             &
                                           surf_def_v(3)%r_canopy,             &
                                           surf_lsm_v(3)%r_canopy,             &
                                           surf_usm_v(3)%r_canopy, n_out )

            CASE ( 'r_s' )
               CALL surface_data_output_sum_up( surf_def_h(0)%r_s,             &
                                           surf_def_h(1)%r_s,                  &
                                           surf_lsm_h%r_s,                     &
                                           surf_usm_h%r_s,                     &
                                           surf_def_v(0)%r_s,                  &
                                           surf_lsm_v(0)%r_s,                  &
                                           surf_usm_v(0)%r_s,                  &
                                           surf_def_v(1)%r_s,                  &
                                           surf_lsm_v(1)%r_s,                  &
                                           surf_usm_v(1)%r_s,                  &
                                           surf_def_v(2)%r_s,                  &
                                           surf_lsm_v(2)%r_s,                  &
                                           surf_usm_v(2)%r_s,                  &
                                           surf_def_v(3)%r_s,                  &
                                           surf_lsm_v(3)%r_s,                  &
                                           surf_usm_v(3)%r_s, n_out )


            CASE ( 'rad_sw_dir' )
               CALL surface_data_output_sum_up( surf_def_h(0)%rad_sw_dir,      &
                                           surf_def_h(1)%rad_sw_dir,           &
                                           surf_lsm_h%rad_sw_dir,              &
                                           surf_usm_h%rad_sw_dir,              &
                                           surf_def_v(0)%rad_sw_dir,           &
                                           surf_lsm_v(0)%rad_sw_dir,           &
                                           surf_usm_v(0)%rad_sw_dir,           &
                                           surf_def_v(1)%rad_sw_dir,           &
                                           surf_lsm_v(1)%rad_sw_dir,           &
                                           surf_usm_v(1)%rad_sw_dir,           &
                                           surf_def_v(2)%rad_sw_dir,           &
                                           surf_lsm_v(2)%rad_sw_dir,           &
                                           surf_usm_v(2)%rad_sw_dir,           &
                                           surf_def_v(3)%rad_sw_dir,           &
                                           surf_lsm_v(3)%rad_sw_dir,           &
                                           surf_usm_v(3)%rad_sw_dir, n_out )
            CASE ( 'rad_sw_dif' )
               CALL surface_data_output_sum_up( surf_def_h(0)%rad_sw_dif,      &
                                           surf_def_h(1)%rad_sw_dif,           &
                                           surf_lsm_h%rad_sw_dif,              &
                                           surf_usm_h%rad_sw_dif,              &
                                           surf_def_v(0)%rad_sw_dif,           &
                                           surf_lsm_v(0)%rad_sw_dif,           &
                                           surf_usm_v(0)%rad_sw_dif,           &
                                           surf_def_v(1)%rad_sw_dif,           &
                                           surf_lsm_v(1)%rad_sw_dif,           &
                                           surf_usm_v(1)%rad_sw_dif,           &
                                           surf_def_v(2)%rad_sw_dif,           &
                                           surf_lsm_v(2)%rad_sw_dif,           &
                                           surf_usm_v(2)%rad_sw_dif,           &
                                           surf_def_v(3)%rad_sw_dif,           &
                                           surf_lsm_v(3)%rad_sw_dif,           &
                                           surf_usm_v(3)%rad_sw_dif, n_out )

            CASE ( 'rad_sw_ref' )
               CALL surface_data_output_sum_up( surf_def_h(0)%rad_sw_ref,      &
                                           surf_def_h(1)%rad_sw_ref,           &
                                           surf_lsm_h%rad_sw_ref,              &
                                           surf_usm_h%rad_sw_ref,              &
                                           surf_def_v(0)%rad_sw_ref,           &
                                           surf_lsm_v(0)%rad_sw_ref,           &
                                           surf_usm_v(0)%rad_sw_ref,           &
                                           surf_def_v(1)%rad_sw_ref,           &
                                           surf_lsm_v(1)%rad_sw_ref,           &
                                           surf_usm_v(1)%rad_sw_ref,           &
                                           surf_def_v(2)%rad_sw_ref,           &
                                           surf_lsm_v(2)%rad_sw_ref,           &
                                           surf_usm_v(2)%rad_sw_ref,           &
                                           surf_def_v(3)%rad_sw_ref,           &
                                           surf_lsm_v(3)%rad_sw_ref,           &
                                           surf_usm_v(3)%rad_sw_ref, n_out )

            CASE ( 'rad_sw_res' )
               CALL surface_data_output_sum_up( surf_def_h(0)%rad_sw_res,      &
                                           surf_def_h(1)%rad_sw_res,           &
                                           surf_lsm_h%rad_sw_res,              &
                                           surf_usm_h%rad_sw_res,              &
                                           surf_def_v(0)%rad_sw_res,           &
                                           surf_lsm_v(0)%rad_sw_res,           &
                                           surf_usm_v(0)%rad_sw_res,           &
                                           surf_def_v(1)%rad_sw_res,           &
                                           surf_lsm_v(1)%rad_sw_res,           &
                                           surf_usm_v(1)%rad_sw_res,           &
                                           surf_def_v(2)%rad_sw_res,           &
                                           surf_lsm_v(2)%rad_sw_res,           &
                                           surf_usm_v(2)%rad_sw_res,           &
                                           surf_def_v(3)%rad_sw_res,           &
                                           surf_lsm_v(3)%rad_sw_res,           &
                                           surf_usm_v(3)%rad_sw_res, n_out )

            CASE ( 'rad_lw_dif' )
               CALL surface_data_output_sum_up( surf_def_h(0)%rad_lw_dif,      &
                                           surf_def_h(1)%rad_lw_dif,           &
                                           surf_lsm_h%rad_lw_dif,              &
                                           surf_usm_h%rad_lw_dif,              &
                                           surf_def_v(0)%rad_lw_dif,           &
                                           surf_lsm_v(0)%rad_lw_dif,           &
                                           surf_usm_v(0)%rad_lw_dif,           &
                                           surf_def_v(1)%rad_lw_dif,           &
                                           surf_lsm_v(1)%rad_lw_dif,           &
                                           surf_usm_v(1)%rad_lw_dif,           &
                                           surf_def_v(2)%rad_lw_dif,           &
                                           surf_lsm_v(2)%rad_lw_dif,           &
                                           surf_usm_v(2)%rad_lw_dif,           &
                                           surf_def_v(3)%rad_lw_dif,           &
                                           surf_lsm_v(3)%rad_lw_dif,           &
                                           surf_usm_v(3)%rad_lw_dif, n_out )

            CASE ( 'rad_lw_ref' )
               CALL surface_data_output_sum_up( surf_def_h(0)%rad_lw_ref,      &
                                           surf_def_h(1)%rad_lw_ref,           &
                                           surf_lsm_h%rad_lw_ref,              &
                                           surf_usm_h%rad_lw_ref,              &
                                           surf_def_v(0)%rad_lw_ref,           &
                                           surf_lsm_v(0)%rad_lw_ref,           &
                                           surf_usm_v(0)%rad_lw_ref,           &
                                           surf_def_v(1)%rad_lw_ref,           &
                                           surf_lsm_v(1)%rad_lw_ref,           &
                                           surf_usm_v(1)%rad_lw_ref,           &
                                           surf_def_v(2)%rad_lw_ref,           &
                                           surf_lsm_v(2)%rad_lw_ref,           &
                                           surf_usm_v(2)%rad_lw_ref,           &
                                           surf_def_v(3)%rad_lw_ref,           &
                                           surf_lsm_v(3)%rad_lw_ref,           &
                                           surf_usm_v(3)%rad_lw_ref, n_out )

            CASE ( 'rad_lw_res' )
               CALL surface_data_output_sum_up( surf_def_h(0)%rad_lw_res,      &
                                           surf_def_h(1)%rad_lw_res,           &
                                           surf_lsm_h%rad_lw_res,              &
                                           surf_usm_h%rad_lw_res,              &
                                           surf_def_v(0)%rad_lw_res,           &
                                           surf_lsm_v(0)%rad_lw_res,           &
                                           surf_usm_v(0)%rad_lw_res,           &
                                           surf_def_v(1)%rad_lw_res,           &
                                           surf_lsm_v(1)%rad_lw_res,           &
                                           surf_usm_v(1)%rad_lw_res,           &
                                           surf_def_v(2)%rad_lw_res,           &
                                           surf_lsm_v(2)%rad_lw_res,           &
                                           surf_usm_v(2)%rad_lw_res,           &
                                           surf_def_v(3)%rad_lw_res,           &
                                           surf_lsm_v(3)%rad_lw_res,           &
                                           surf_usm_v(3)%rad_lw_res, n_out )

            CASE ( 'uvw1' )
               CALL surface_data_output_sum_up( surf_def_h(0)%uvw_abs,         &
                                           surf_def_h(1)%uvw_abs,              &
                                           surf_lsm_h%uvw_abs,                 &
                                           surf_usm_h%uvw_abs,                 &
                                           surf_def_v(0)%uvw_abs,              &
                                           surf_lsm_v(0)%uvw_abs,              &
                                           surf_usm_v(0)%uvw_abs,              &
                                           surf_def_v(1)%uvw_abs,              &
                                           surf_lsm_v(1)%uvw_abs,              &
                                           surf_usm_v(1)%uvw_abs,              &
                                           surf_def_v(2)%uvw_abs,              &
                                           surf_lsm_v(2)%uvw_abs,              &
                                           surf_usm_v(2)%uvw_abs,              &
                                           surf_def_v(3)%uvw_abs,              &
                                           surf_lsm_v(3)%uvw_abs,              &
                                           surf_usm_v(3)%uvw_abs, n_out )

            CASE ( 'waste_heat' )
               CALL surface_data_output_sum_up( surf_def_h(0)%waste_heat,      &
                                           surf_def_h(1)%waste_heat,           &
                                           surf_lsm_h%waste_heat,              &
                                           surf_usm_h%waste_heat,              &
                                           surf_def_v(0)%waste_heat,           &
                                           surf_lsm_v(0)%waste_heat,           &
                                           surf_usm_v(0)%waste_heat,           &
                                           surf_def_v(1)%waste_heat,           &
                                           surf_lsm_v(1)%waste_heat,           &
                                           surf_usm_v(1)%waste_heat,           &
                                           surf_def_v(2)%waste_heat,           &
                                           surf_lsm_v(2)%waste_heat,           &
                                           surf_usm_v(2)%waste_heat,           &
                                           surf_def_v(3)%waste_heat,           &
                                           surf_lsm_v(3)%waste_heat,           &
                                           surf_usm_v(3)%waste_heat, n_out )

            CASE ( 'im_hf' )
               CALL surface_data_output_sum_up( surf_def_h(0)%iwghf_eb,        &
                                           surf_def_h(1)%iwghf_eb,             &
                                           surf_lsm_h%iwghf_eb,                &
                                           surf_usm_h%iwghf_eb,                &
                                           surf_def_v(0)%iwghf_eb,             &
                                           surf_lsm_v(0)%iwghf_eb,             &
                                           surf_usm_v(0)%iwghf_eb,             &
                                           surf_def_v(1)%iwghf_eb,             &
                                           surf_lsm_v(1)%iwghf_eb,             &
                                           surf_usm_v(1)%iwghf_eb,             &
                                           surf_def_v(2)%iwghf_eb,             &
                                           surf_lsm_v(2)%iwghf_eb,             &
                                           surf_usm_v(2)%iwghf_eb,             &
                                           surf_def_v(3)%iwghf_eb,             &
                                           surf_lsm_v(3)%iwghf_eb,             &
                                           surf_usm_v(3)%iwghf_eb, n_out )

         END SELECT
      ENDDO


   END SUBROUTINE surface_data_output_averaging

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Sum-up the surface data for average output variables.
!------------------------------------------------------------------------------!
   SUBROUTINE surface_data_output_sum_up( var_def_h0, var_def_h1,              &
                                     var_lsm_h,  var_usm_h,                    &
                                     var_def_v0, var_lsm_v0, var_usm_v0,       &
                                     var_def_v1, var_lsm_v1, var_usm_v1,       &
                                     var_def_v2, var_lsm_v2, var_usm_v2,       &
                                     var_def_v3, var_lsm_v3, var_usm_v3, n_out )

      IMPLICIT NONE

      INTEGER(iwp) ::  m          !< running index for surface elements
      INTEGER(iwp) ::  n_out      !< index for output variable
      INTEGER(iwp) ::  n_surf     !< running index for surface elements

      REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(IN) ::  var_def_h0 !< output variable at upward-facing default-type surfaces
      REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(IN) ::  var_def_h1 !< output variable at downward-facing default-type surfaces
      REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(IN) ::  var_lsm_h  !< output variable at upward-facing natural-type surfaces
      REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(IN) ::  var_usm_h  !< output variable at upward-facing urban-type surfaces
      REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(IN) ::  var_def_v0 !< output variable at northward-facing default-type surfaces
      REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(IN) ::  var_def_v1 !< output variable at southward-facing default-type surfaces
      REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(IN) ::  var_def_v2 !< output variable at eastward-facing default-type surfaces
      REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(IN) ::  var_def_v3 !< output variable at westward-facing default-type surfaces
      REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(IN) ::  var_lsm_v0 !< output variable at northward-facing natural-type surfaces
      REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(IN) ::  var_lsm_v1 !< output variable at southward-facing natural-type surfaces
      REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(IN) ::  var_lsm_v2 !< output variable at eastward-facing natural-type surfaces
      REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(IN) ::  var_lsm_v3 !< output variable at westward-facing natural-type surfaces
      REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(IN) ::  var_usm_v0 !< output variable at northward-facing urban-type surfaces
      REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(IN) ::  var_usm_v1 !< output variable at southward-facing urban-type surfaces
      REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(IN) ::  var_usm_v2 !< output variable at eastward-facing urban-type surfaces
      REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(IN) ::  var_usm_v3 !< output variable at westward-facing urban-type surfaces

!
!--   Set counter variable to zero before the variable is written to
!--   the output array.
      n_surf = 0

!
!--   Write the horizontal surfaces.
!--   Before each the variable is written to the output data structure, first
!--   check if the variable for the respective surface type is defined.
!--   If a variable is not defined, skip the block and increment the counter
!--   variable by the number of surface elements of this type. Usually this
!--   is zere, however, there might be the situation that e.g. urban surfaces
!--   are defined but the respective variable is not allocated for this surface
!--   type. To write the data on the exact position, increment the counter.
      IF ( ALLOCATED( var_def_h0 ) )  THEN
         DO  m = 1, surf_def_h(0)%ns
            n_surf                        = n_surf + 1
            surfaces%var_av(n_surf,n_out) = surfaces%var_av(n_surf,n_out)      &
                                          + var_def_h0(m)
         ENDDO
      ELSE
         n_surf = n_surf + surf_def_h(0)%ns
      ENDIF
      IF ( ALLOCATED( var_def_h1 ) )  THEN
         DO  m = 1, surf_def_h(1)%ns
            n_surf                   = n_surf + 1
            surfaces%var_av(n_surf,n_out) = surfaces%var_av(n_surf,n_out)      &
                                          + var_def_h1(m)
         ENDDO
      ELSE
         n_surf = n_surf + surf_def_h(1)%ns
      ENDIF
      IF ( ALLOCATED( var_lsm_h ) )  THEN
         DO  m = 1, surf_lsm_h%ns
            n_surf                        = n_surf + 1
            surfaces%var_av(n_surf,n_out) = surfaces%var_av(n_surf,n_out)      &
                                          + var_lsm_h(m)
         ENDDO
      ELSE
         n_surf = n_surf + surf_lsm_h%ns
      ENDIF
      IF ( ALLOCATED( var_usm_h ) )  THEN
         DO  m = 1, surf_usm_h%ns
            n_surf                        = n_surf + 1
            surfaces%var_av(n_surf,n_out) = surfaces%var_av(n_surf,n_out)      &
                                          + var_usm_h(m)
         ENDDO
      ELSE
         n_surf = n_surf + surf_usm_h%ns
      ENDIF
!
!--   Write northward-facing
      IF ( ALLOCATED( var_def_v0 ) )  THEN
         DO  m = 1, surf_def_v(0)%ns
            n_surf                        = n_surf + 1
            surfaces%var_av(n_surf,n_out) = surfaces%var_av(n_surf,n_out)      &
                                          + var_def_v0(m)
         ENDDO
      ELSE
         n_surf = n_surf + surf_def_v(0)%ns
      ENDIF
      IF ( ALLOCATED( var_lsm_v0 ) )  THEN
         DO  m = 1, surf_lsm_v(0)%ns
            n_surf                        = n_surf + 1
            surfaces%var_av(n_surf,n_out) = surfaces%var_av(n_surf,n_out)      &
                                          + var_lsm_v0(m)
         ENDDO
      ELSE
         n_surf = n_surf + surf_lsm_v(0)%ns
      ENDIF
      IF ( ALLOCATED( var_usm_v0 ) )  THEN
         DO  m = 1, surf_usm_v(0)%ns
            n_surf                        = n_surf + 1
            surfaces%var_av(n_surf,n_out) = surfaces%var_av(n_surf,n_out)      &
                                          + var_usm_v0(m)
         ENDDO
      ELSE
         n_surf = n_surf + surf_usm_v(0)%ns
      ENDIF
!
!--   Write southward-facing
      IF ( ALLOCATED( var_def_v1 ) )  THEN
         DO  m = 1, surf_def_v(1)%ns
            n_surf                        = n_surf + 1
            surfaces%var_av(n_surf,n_out) = surfaces%var_av(n_surf,n_out)      &
                                          + var_def_v1(m)
         ENDDO
      ELSE
         n_surf = n_surf + surf_def_v(1)%ns
      ENDIF
      IF ( ALLOCATED( var_lsm_v1 ) )  THEN
         DO  m = 1, surf_lsm_v(1)%ns
            n_surf                        = n_surf + 1
            surfaces%var_av(n_surf,n_out) = surfaces%var_av(n_surf,n_out)      &
                                          + var_lsm_v1(m)
         ENDDO
      ELSE
         n_surf = n_surf + surf_lsm_v(1)%ns
      ENDIF
      IF ( ALLOCATED( var_usm_v1 ) )  THEN
         DO  m = 1, surf_usm_v(1)%ns
            n_surf                        = n_surf + 1
            surfaces%var_av(n_surf,n_out) = surfaces%var_av(n_surf,n_out)      &
                                          + var_usm_v1(m)
         ENDDO
      ELSE
         n_surf = n_surf + surf_usm_v(1)%ns
      ENDIF
!
!--   Write eastward-facing
      IF ( ALLOCATED( var_def_v2 ) )  THEN
         DO  m = 1, surf_def_v(2)%ns
            n_surf                        = n_surf + 1
            surfaces%var_av(n_surf,n_out) = surfaces%var_av(n_surf,n_out)      &
                                          + var_def_v2(m)
         ENDDO
      ELSE
         n_surf = n_surf + surf_def_v(2)%ns
      ENDIF
      IF ( ALLOCATED( var_lsm_v2 ) )  THEN
         DO  m = 1, surf_lsm_v(2)%ns
            n_surf                        = n_surf + 1
            surfaces%var_av(n_surf,n_out) = surfaces%var_av(n_surf,n_out)      &
                                          + var_lsm_v2(m)
         ENDDO
      ELSE
         n_surf = n_surf + surf_lsm_v(2)%ns
      ENDIF
      IF ( ALLOCATED( var_usm_v2 ) )  THEN
         DO  m = 1, surf_usm_v(2)%ns
            n_surf                        = n_surf + 1
            surfaces%var_av(n_surf,n_out) = surfaces%var_av(n_surf,n_out)      &
                                          + var_usm_v2(m)
         ENDDO
      ELSE
         n_surf = n_surf + surf_usm_v(2)%ns
      ENDIF
!
!--   Write westward-facing
      IF ( ALLOCATED( var_def_v3 ) )  THEN
         DO  m = 1, surf_def_v(3)%ns
            n_surf                        = n_surf + 1
            surfaces%var_av(n_surf,n_out) = surfaces%var_av(n_surf,n_out)      &
                                          + var_def_v3(m)
         ENDDO
      ELSE
         n_surf = n_surf + surf_def_v(3)%ns
      ENDIF
      IF ( ALLOCATED( var_lsm_v3 ) )  THEN
         DO  m = 1, surf_lsm_v(3)%ns
            n_surf                        = n_surf + 1
            surfaces%var_av(n_surf,n_out) = surfaces%var_av(n_surf,n_out)      &
                                          + var_lsm_v3(m)
         ENDDO
      ELSE
         n_surf = n_surf + surf_lsm_v(3)%ns
      ENDIF
      IF ( ALLOCATED( var_usm_v3 ) )  THEN
         DO  m = 1, surf_usm_v(3)%ns
            n_surf                        = n_surf + 1
            surfaces%var_av(n_surf,n_out) = surfaces%var_av(n_surf,n_out)      &
                                          + var_usm_v3(m)
         ENDDO
      ELSE
         n_surf = n_surf + surf_usm_v(3)%ns
      ENDIF

   END SUBROUTINE surface_data_output_sum_up

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Collect the surface data from different types and different orientation.
!------------------------------------------------------------------------------!
   SUBROUTINE surface_data_output_collect( var_def_h0, var_def_h1,             &
                                      var_lsm_h,  var_usm_h,                   &
                                      var_def_v0, var_lsm_v0, var_usm_v0,      &
                                      var_def_v1, var_lsm_v1, var_usm_v1,      &
                                      var_def_v2, var_lsm_v2, var_usm_v2,      &
                                      var_def_v3, var_lsm_v3, var_usm_v3 )

      IMPLICIT NONE

      INTEGER(iwp) ::  m      !< running index for surface elements
      INTEGER(iwp) ::  n_surf !< running index for surface elements

      REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(IN) ::  var_def_h0 !< output variable at upward-facing default-type surfaces
      REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(IN) ::  var_def_h1 !< output variable at downward-facing default-type surfaces
      REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(IN) ::  var_lsm_h  !< output variable at upward-facing natural-type surfaces
      REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(IN) ::  var_usm_h  !< output variable at upward-facing urban-type surfaces
      REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(IN) ::  var_def_v0 !< output variable at northward-facing default-type surfaces
      REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(IN) ::  var_def_v1 !< output variable at southward-facing default-type surfaces
      REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(IN) ::  var_def_v2 !< output variable at eastward-facing default-type surfaces
      REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(IN) ::  var_def_v3 !< output variable at westward-facing default-type surfaces
      REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(IN) ::  var_lsm_v0 !< output variable at northward-facing natural-type surfaces
      REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(IN) ::  var_lsm_v1 !< output variable at southward-facing natural-type surfaces
      REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(IN) ::  var_lsm_v2 !< output variable at eastward-facing natural-type surfaces
      REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(IN) ::  var_lsm_v3 !< output variable at westward-facing natural-type surfaces
      REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(IN) ::  var_usm_v0 !< output variable at northward-facing urban-type surfaces
      REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(IN) ::  var_usm_v1 !< output variable at southward-facing urban-type surfaces
      REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(IN) ::  var_usm_v2 !< output variable at eastward-facing urban-type surfaces
      REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(IN) ::  var_usm_v3 !< output variable at westward-facing urban-type surfaces

!
!--   Set counter variable to zero before the variable is written to
!--   the output array.
      n_surf = 0

!
!--   Write the horizontal surfaces.
!--   Before each the variable is written to the output data structure, first
!--   check if the variable for the respective surface type is defined.
!--   If a variable is not defined, skip the block and increment the counter
!--   variable by the number of surface elements of this type. Usually this
!--   is zere, however, there might be the situation that e.g. urban surfaces
!--   are defined but the respective variable is not allocated for this surface
!--   type. To write the data on the exact position, increment the counter.
      IF ( ALLOCATED( var_def_h0 ) )  THEN
         DO  m = 1, surf_def_h(0)%ns
            n_surf                   = n_surf + 1
            surfaces%var_out(n_surf) = var_def_h0(m)
         ENDDO
      ELSE
         n_surf = n_surf + surf_def_h(0)%ns
      ENDIF
      IF ( ALLOCATED( var_def_h1 ) )  THEN
         DO  m = 1, surf_def_h(1)%ns
            n_surf                   = n_surf + 1
            surfaces%var_out(n_surf) = var_def_h1(m)
         ENDDO
      ELSE
         n_surf = n_surf + surf_def_h(1)%ns
      ENDIF
      IF ( ALLOCATED( var_lsm_h ) )  THEN
         DO  m = 1, surf_lsm_h%ns
            n_surf                   = n_surf + 1
            surfaces%var_out(n_surf) = var_lsm_h(m)
         ENDDO
      ELSE
         n_surf = n_surf + surf_lsm_h%ns
      ENDIF
      IF ( ALLOCATED( var_usm_h ) )  THEN
         DO  m = 1, surf_usm_h%ns
            n_surf                   = n_surf + 1
            surfaces%var_out(n_surf) = var_usm_h(m)
         ENDDO
      ELSE
         n_surf = n_surf + surf_usm_h%ns
      ENDIF
!
!--   Write northward-facing
      IF ( ALLOCATED( var_def_v0 ) )  THEN
         DO  m = 1, surf_def_v(0)%ns
            n_surf                   = n_surf + 1
            surfaces%var_out(n_surf) = var_def_v0(m)
         ENDDO
      ELSE
         n_surf = n_surf + surf_def_v(0)%ns
      ENDIF
      IF ( ALLOCATED( var_lsm_v0 ) )  THEN
         DO  m = 1, surf_lsm_v(0)%ns
            n_surf                   = n_surf + 1
            surfaces%var_out(n_surf) = var_lsm_v0(m)
         ENDDO
      ELSE
         n_surf = n_surf + surf_lsm_v(0)%ns
      ENDIF
      IF ( ALLOCATED( var_usm_v0 ) )  THEN
         DO  m = 1, surf_usm_v(0)%ns
            n_surf                   = n_surf + 1
            surfaces%var_out(n_surf) = var_usm_v0(m)
         ENDDO
      ELSE
         n_surf = n_surf + surf_usm_v(0)%ns
      ENDIF
!
!--   Write southward-facing
      IF ( ALLOCATED( var_def_v1 ) )  THEN
         DO  m = 1, surf_def_v(1)%ns
            n_surf                   = n_surf + 1
            surfaces%var_out(n_surf) = var_def_v1(m)
         ENDDO
      ELSE
         n_surf = n_surf + surf_def_v(1)%ns
      ENDIF
      IF ( ALLOCATED( var_lsm_v1 ) )  THEN
         DO  m = 1, surf_lsm_v(1)%ns
            n_surf                   = n_surf + 1
            surfaces%var_out(n_surf) = var_lsm_v1(m)
         ENDDO
      ELSE
         n_surf = n_surf + surf_lsm_v(1)%ns
      ENDIF
      IF ( ALLOCATED( var_usm_v1 ) )  THEN
         DO  m = 1, surf_usm_v(1)%ns
            n_surf                   = n_surf + 1
            surfaces%var_out(n_surf) = var_usm_v1(m)
         ENDDO
      ELSE
         n_surf = n_surf + surf_usm_v(1)%ns
      ENDIF
!
!--   Write eastward-facing
      IF ( ALLOCATED( var_def_v2 ) )  THEN
         DO  m = 1, surf_def_v(2)%ns
            n_surf                   = n_surf + 1
            surfaces%var_out(n_surf) = var_def_v2(m)
         ENDDO
      ELSE
         n_surf = n_surf + surf_def_v(2)%ns
      ENDIF
      IF ( ALLOCATED( var_lsm_v2 ) )  THEN
         DO  m = 1, surf_lsm_v(2)%ns
            n_surf                   = n_surf + 1
            surfaces%var_out(n_surf) = var_lsm_v2(m)
         ENDDO
      ELSE
         n_surf = n_surf + surf_lsm_v(2)%ns
      ENDIF
      IF ( ALLOCATED( var_usm_v2 ) )  THEN
         DO  m = 1, surf_usm_v(2)%ns
            n_surf                   = n_surf + 1
            surfaces%var_out(n_surf) = var_usm_v2(m)
         ENDDO
      ELSE
         n_surf = n_surf + surf_usm_v(2)%ns
      ENDIF
!
!--   Write westward-facing
      IF ( ALLOCATED( var_def_v3 ) )  THEN
         DO  m = 1, surf_def_v(3)%ns
            n_surf                   = n_surf + 1
            surfaces%var_out(n_surf) = var_def_v3(m)
         ENDDO
      ELSE
         n_surf = n_surf + surf_def_v(3)%ns
      ENDIF
      IF ( ALLOCATED( var_lsm_v3 ) )  THEN
         DO  m = 1, surf_lsm_v(3)%ns
            n_surf                   = n_surf + 1
            surfaces%var_out(n_surf) = var_lsm_v3(m)
         ENDDO
      ELSE
         n_surf = n_surf + surf_lsm_v(3)%ns
      ENDIF
      IF ( ALLOCATED( var_usm_v3 ) )  THEN
         DO  m = 1, surf_usm_v(3)%ns
            n_surf                   = n_surf + 1
            surfaces%var_out(n_surf) = var_usm_v3(m)
         ENDDO
      ELSE
         n_surf = n_surf + surf_usm_v(3)%ns
      ENDIF

   END SUBROUTINE surface_data_output_collect

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Parin for output of surface parameters
!------------------------------------------------------------------------------!
    SUBROUTINE surface_data_output_parin

       IMPLICIT NONE

       CHARACTER (LEN=80) ::  line  !< dummy string that contains the current line of the parameter file


       NAMELIST /surface_data_output_parameters/                               &
                                  averaging_interval_surf, data_output_surf,   &
                                  dt_dosurf, dt_dosurf_av,                     &
                                  skip_time_dosurf, skip_time_dosurf_av,       &
                                  to_netcdf, to_vtk

       line = ' '

!
!--    Try to find the namelist
       REWIND ( 11 )
       line = ' '
       DO WHILE ( INDEX( line, '&surface_data_output_parameters' ) == 0 )
          READ ( 11, '(A)', END=14 )  line
       ENDDO
       BACKSPACE ( 11 )

!
!--    Read namelist
       READ ( 11, surface_data_output_parameters, ERR = 10 )
!
!--    Set flag that indicates that surface data output is switched on
       surface_output = .TRUE.
       GOTO 14

 10    BACKSPACE( 11 )
       READ( 11 , '(A)') line
       CALL parin_fail_message( 'surface_data_output_parameters', line )

 14    CONTINUE


    END SUBROUTINE surface_data_output_parin


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check the input parameters for consistency. Further pre-process the given
!> output variables, i.e. separate them into average and non-average output
!> variables and map them onto internal output array.
!------------------------------------------------------------------------------!
    SUBROUTINE surface_data_output_check_parameters

       USE control_parameters,                                                 &
           ONLY:  averaging_interval, dt_data_output, indoor_model,            &
                  initializing_actions, message_string

       USE pegrid,                                                             &
           ONLY:  numprocs_previous_run

       IMPLICIT NONE

       CHARACTER(LEN=100) ::  trimvar !< dummy for single output variable
       CHARACTER(LEN=100) ::  unit    !< dummy for unit of output variable

       INTEGER(iwp) ::  av    !< id indicating average or non-average data output
       INTEGER(iwp) ::  ilen  !< string length
       INTEGER(iwp) ::  n_out !< running index for number of output variables
!
!--    Check if any output file type is selected
       IF ( .NOT. to_vtk  .AND.  .NOT. to_netcdf )  THEN
          WRITE( message_string, * )                                     &
             'no output file type selected for surface-data output!&' // &
             'Set at least either "to_vtk" or "to_netcdf" to .TRUE.'
          CALL message( 'surface_data_output_check_parameters',          &
                        'PA0662', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Check the average interval
       IF ( averaging_interval_surf == 9999999.9_wp )  THEN
          averaging_interval_surf = averaging_interval
       ENDIF
!
!--    Set the default data-output interval dt_data_output if necessary
       IF ( dt_dosurf    == 9999999.9_wp )  dt_dosurf    = dt_data_output
       IF ( dt_dosurf_av == 9999999.9_wp )  dt_dosurf_av = dt_data_output

       IF ( averaging_interval_surf > dt_dosurf_av )  THEN
          WRITE( message_string, * )  'averaging_interval_surf = ',            &
                averaging_interval_surf, ' must be <= dt_dosurf_av = ',        &
                dt_dosurf_av
          CALL message( 'surface_data_output_check_parameters',                &
                        'PA0536', 1, 2, 0, 6, 0 )
       ENDIF

#if ! defined( __netcdf4_parallel )
!
!--    Surface output via NetCDF requires parallel NetCDF
       IF ( to_netcdf )  THEN
          message_string = 'to_netcdf = .True. requires parallel NetCDF'
          CALL message( 'surface_data_output_check_parameters',                &
                        'PA0116', 1, 2, 0, 6, 0 )
       ENDIF
#endif
!
!--    In case of parallel NetCDF output the output timestep must not be zero.
!--    This is because the number of requiered output timesteps is
!--    pre-calculated, which is not possible with zero output timestep.
       IF ( netcdf_data_format > 4 )  THEN
          IF ( dt_dosurf == 0.0_wp )  THEN
             message_string = 'dt_dosurf = 0.0 while using a ' //             &
                              'variable timestep and parallel netCDF4 ' //    &
                              'is not allowed.'
             CALL message( 'surface_data_output_check_parameters', 'PA0081',  &
                           1, 2, 0, 6, 0 )
          ENDIF

          IF ( dt_dosurf_av == 0.0_wp )  THEN
             message_string = 'dt_dosurf_av = 0.0 while using a ' //          &
                              'variable timestep and parallel netCDF4 ' //    &
                              'is not allowed.'
             CALL message( 'surface_data_output_check_parameters', 'PA0081',  &
                           1, 2, 0, 6, 0 )
          ENDIF
       ENDIF

!
!--   In case of restart runs, check it the number of cores has been changed.
!--   With surface output this is not allowed.
      IF ( TRIM( initializing_actions ) == 'read_restart_data'  .AND.          &
           numprocs_previous_run /= numprocs ) THEN
         message_string = 'The number of cores has been changed between ' //   &
                          'restart runs. This is not allowed when surface ' // &
                          'data output is used.'
          CALL message( 'surface_data_output_check_parameters',                &
                        'PA0585', 1, 2, 0, 6, 0 )
      ENDIF
!
!--   Count number of output variables and separate output strings for
!--   average and non-average output variables.
      n_out = 0
      DO WHILE ( data_output_surf(n_out+1)(1:1) /= ' ' )

         n_out = n_out + 1
         ilen = LEN_TRIM( data_output_surf(n_out) )
         trimvar = TRIM( data_output_surf(n_out) )

!
!--      Check for data averaging
         av = 0
         IF ( ilen > 3 )  THEN
            IF ( data_output_surf(n_out)(ilen-2:ilen) == '_av' )  THEN
               trimvar = data_output_surf(n_out)(1:ilen-3)
               av      = 1
            ENDIF
         ENDIF

         dosurf_no(av) = dosurf_no(av) + 1
         dosurf(av,dosurf_no(av)) = TRIM( trimvar )

!
!--      Check if all output variables are known and assign a unit
         unit = 'not set'
         SELECT CASE ( TRIM( trimvar ) )

            CASE ( 'css', 'cssws', 'qsws_liq', 'qsws_soil', 'qsws_veg' )
               message_string = TRIM( trimvar ) //                             &
                             ' is not yet implemented in the surface output'
               CALL message( 'surface_data_output_check_parameters',           &
                             'PA0537', 1, 2, 0, 6, 0 )

            CASE ( 'us', 'uvw1' )
               unit = 'm/s'

            CASE ( 'ss', 'qcs', 'ncs', 'qrs', 'nrs' )
               unit = '1'

            CASE ( 'z0', 'z0h', 'z0q', 'ol' )
               unit = 'm'

            CASE ( 'ts', 'theta1', 'thetav1', 'theta_surface', 'thetav_surface' )
               unit = 'K'

            CASE ( 'usws', 'vsws' )
               unit = 'm2/s2'

            CASE ( 'qcsws', 'ncsws', 'qrsws', 'nrsws', 'sasws' )

            CASE ( 'shf' )
               unit = 'K m/s'

            CASE ( 'qsws' )
               unit = 'kg/kg m/s'

            CASE ( 'ssws' )
               unit = 'kg/m2/s'

            CASE ( 'qs', 'q_surface', 'qv1' )
               unit = 'kg/kg'

            CASE ( 'rad_net' )
               unit = 'W/m2'

            CASE ( 'rad_lw_in', 'rad_lw_out', 'rad_lw_dif', 'rad_lw_ref',     &
                   'rad_lw_res' )
               unit = 'W/m2'

            CASE ( 'rad_sw_in', 'rad_sw_out', 'rad_sw_dif', 'rad_sw_ref',     &
                   'rad_sw_res', 'rad_sw_dir' )
               unit = 'W/m2'

            CASE ( 'ghf' )
               unit = 'W/m2'

            CASE ( 'r_a', 'r_canopy', 'r_soil', 'r_s' )
               unit = 's/m'

            CASE ( 'waste_heat', 'im_hf' )
               IF ( .NOT. indoor_model )  THEN
                  message_string = TRIM( trimvar ) //                          &
                             ' requires the indoor model'
               CALL message( 'surface_data_output_check_parameters',           &
                             'PA0588', 1, 2, 0, 6, 0 )
               ENDIF

               unit = 'W/m2'

            CASE DEFAULT
               message_string = TRIM( trimvar ) //                             &
                             ' is not part of the surface output'
               CALL message( 'surface_data_output_check_parameters',           &
                             'PA0538', 1, 2, 0, 6, 0 )
         END SELECT

         dosurf_unit(av,dosurf_no(av)) = unit

       ENDDO

    END SUBROUTINE surface_data_output_check_parameters


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Last action.
!------------------------------------------------------------------------------!
    SUBROUTINE surface_data_output_last_action( av )

      USE control_parameters,                                                  &
          ONLY:  io_blocks, io_group

#if defined( __parallel )
      USE pegrid,                                                              &
          ONLY:  comm2d, ierr
#endif

      IMPLICIT NONE

      INTEGER(iwp) ::  av     !< id indicating average or non-average data output
      INTEGER(iwp) ::  i      !< loop index

!
!--   Return, if nothing to output
      IF ( dosurf_no(av) == 0 )  RETURN
!
!--   If output to VTK files is enabled, check if files are open and write
!--   an end-of-file statement.
      IF ( to_vtk )  THEN
         CALL check_open( 25+av )
!
!--      Write time coordinate
         DO  i = 0, io_blocks-1
            IF ( i == io_group )  THEN
               WRITE ( 25+av )  LEN_TRIM( 'END' )
               WRITE ( 25+av )  'END'
            ENDIF
#if defined( __parallel )
            CALL MPI_BARRIER( comm2d, ierr )
#endif
         ENDDO
      ENDIF

    END SUBROUTINE surface_data_output_last_action

!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine reads globally used restart data.
!------------------------------------------------------------------------------!
    SUBROUTINE surface_data_output_rrd_global( found )


       USE control_parameters,                                                 &
           ONLY: length, restart_string

       IMPLICIT NONE

       LOGICAL, INTENT(OUT)  ::  found !< flag indicating if variable was found

       found = .TRUE.

       SELECT CASE ( restart_string(1:length) )

          CASE ( 'average_count_surf' )
             READ ( 13 )  average_count_surf

          CASE DEFAULT

             found = .FALSE.

       END SELECT


    END SUBROUTINE surface_data_output_rrd_global

!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine reads the respective restart data.
!------------------------------------------------------------------------------!
    SUBROUTINE surface_data_output_rrd_local( found )


       USE control_parameters,                                                 &
           ONLY: length, restart_string

       IMPLICIT NONE

       LOGICAL, INTENT(OUT)  ::  found

!
!--    Here the reading of user-defined restart data follows:
!--    Sample for user-defined output
       found = .TRUE.

       SELECT CASE ( restart_string(1:length) )

          CASE ( 'surfaces%var_av' )
             READ ( 13 )  surfaces%var_av

          CASE DEFAULT

             found = .FALSE.

          END SELECT


    END SUBROUTINE surface_data_output_rrd_local

!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine writes the respective restart data.
!------------------------------------------------------------------------------!
    SUBROUTINE surface_data_output_wrd_global

       IMPLICIT NONE

       CALL wrd_write_string( 'average_count_surf' )
       WRITE ( 14 )  average_count_surf

    END SUBROUTINE surface_data_output_wrd_global

!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine writes restart data which individual on each PE
!------------------------------------------------------------------------------!
    SUBROUTINE surface_data_output_wrd_local

       IMPLICIT NONE

         IF ( ALLOCATED( surfaces%var_av ) )  THEN
            CALL wrd_write_string( 'surfaces%var_av' )
            WRITE ( 14 )  surfaces%var_av
         ENDIF


    END SUBROUTINE surface_data_output_wrd_local


END MODULE surface_data_output_mod
