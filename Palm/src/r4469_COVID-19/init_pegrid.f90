!> @file init_pegrid.f90
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
! $Id: init_pegrid.f90 4461 2020-03-12 16:51:59Z raasch $
! communicator configurations for four virtual pe grids defined
! 
! 4444 2020-03-05 15:59:50Z raasch
! bugfix: cpp-directives for serial mode added
! 
! 4360 2020-01-07 11:25:50Z suehring
! changed message PA0467
! 
! 4264 2019-10-15 16:00:23Z scharf
! corrected error message string
! 
! 4241 2019-09-27 06:32:47Z raasch
! Check added to ensure that subdomain grid has at least the size as given by the number
! of ghost points
! 
! 4182 2019-08-22 15:20:23Z scharf
! Corrected "Former revisions" section
! 
! 4045 2019-06-21 10:58:47Z raasch
! bugfix: kind attribute added to nint function to allow for large integers which may appear in
! case of default recycling width and small grid spacings
! 
! 3999 2019-05-23 16:09:37Z suehring
! Spend 3 ghost points also in case of pw-scheme when nesting is applied
! 
! 3897 2019-04-15 11:51:14Z suehring
! Minor revision of multigrid check; give warning instead of an abort.
! 
! 3890 2019-04-12 15:59:20Z suehring
! Check if grid coarsening is possible on subdomain, in order to avoid that 
! multigrid approach effectively reduces to a Gauss-Seidel scheme.
! 
! 3885 2019-04-11 11:29:34Z kanani
! Changes related to global restructuring of location messages and introduction 
! of additional debug messages
! 
! 3884 2019-04-10 13:31:55Z Giersch
! id_recycling is only calculated in case of tubulent inflow
! 
! 3761 2019-02-25 15:31:42Z raasch
! unused variable removed
! 
! 3655 2019-01-07 16:51:22Z knoop
! variables documented
!
! Revision 1.1  1997/07/24 11:15:09  raasch
! Initial revision
!
!
! Description:
! ------------
!> Determination of the virtual processor topology (if not prescribed by the
!> user)and computation of the grid point number and array bounds of the local
!> domains.
!> @todo: remove MPI-data types for 2D exchange on coarse multigrid level (not
!>        used any more) 
!------------------------------------------------------------------------------!
 SUBROUTINE init_pegrid
 

    USE control_parameters,                                                    &
        ONLY:  bc_dirichlet_l, bc_dirichlet_n, bc_dirichlet_r, bc_dirichlet_s, &
               bc_lr, bc_ns, bc_radiation_l, bc_radiation_n, bc_radiation_r,   &
               bc_radiation_s, &
               grid_level, grid_level_count, maximum_grid_level,               &
               message_string, mg_switch_to_pe0_level,         &
               psolver


#if defined( __parallel )
    USE control_parameters,                                                    &
        ONLY:  coupling_mode, coupling_topology, gathered_size, momentum_advec, &
               outflow_source_plane, recycling_width, scalar_advec, subdomain_size, &
               turbulent_inflow, turbulent_outflow, y_shift

    USE grid_variables,                                                        &
        ONLY:  dx
#endif
        
    USE indices,                                                               &
        ONLY:  nnx, nny, nnz, nx, nxl, nxl_mg,   &
               nxlu, nxr, nxr_mg, ny, nyn, nyn_mg, nys, nys_mg,    &
               nysv, nz, nzb, nzt, nzt_mg, wall_flags_1, wall_flags_2,         &
               wall_flags_3, wall_flags_4, wall_flags_5, wall_flags_6,         &
               wall_flags_7, wall_flags_8, wall_flags_9, wall_flags_10

#if defined( __parallel )
    USE indices,                                                               &
        ONLY:  mg_loc_ind, nbgp, nx_a, nx_o, ny_a, ny_o
#endif

    USE kinds
      
    USE pegrid
    
#if defined( __parallel )
    USE pmc_interface,                                                         &
        ONLY:  nested_run

    USE spectra_mod,                                                           &
        ONLY:  calculate_spectra

    USE synthetic_turbulence_generator_mod,                                    &
        ONLY:  id_stg_left, id_stg_north, id_stg_right, id_stg_south,          &
               use_syn_turb_gen
#endif

    USE transpose_indices,                                                     &
        ONLY:  nxl_y, nxl_z, nxr_y, nxr_z, nyn_x, nyn_z, nys_x,&
               nys_z, nzb_x, nzb_y, nzt_x, nzt_y

#if defined( __parallel )
    USE transpose_indices,                                                     &
        ONLY:  nxl_yd, nxr_yd, nzb_yd, nzt_yd

    USE vertical_nesting_mod,                                                  &
        ONLY:  vnested, vnest_init_pegrid_domain, vnest_init_pegrid_rank
#endif

    IMPLICIT NONE

    INTEGER(iwp) ::  i                        !< running index over number of processors or number of multigrid level
#if defined( __parallel )
    INTEGER(iwp) ::  id_inflow_l              !< ID indicating processors located at the left inflow boundary
    INTEGER(iwp) ::  id_outflow_l             !< local value of id_outflow
    INTEGER(iwp) ::  id_outflow_source_l      !< local value of id_outflow_source
    INTEGER(iwp) ::  id_recycling_l           !< ID indicating processors located at the recycling plane
    INTEGER(iwp) ::  id_stg_left_l            !< left lateral boundary local core id in case of turbulence generator  
    INTEGER(iwp) ::  id_stg_north_l           !< north lateral boundary local core id in case of turbulence generator  
    INTEGER(iwp) ::  id_stg_right_l           !< right lateral boundary local core id in case of turbulence generator  
    INTEGER(iwp) ::  id_stg_south_l           !< south lateral boundary local core id in case of turbulence generator  
    INTEGER(iwp) ::  ind(5)                   !< array containing the subdomain bounds
#endif
    INTEGER(iwp) ::  j                        !< running index, used for various loops
    INTEGER(iwp) ::  k                        !< number of vertical grid points in different multigrid level
    INTEGER(iwp) ::  maximum_grid_level_l     !< maximum number of grid level without switching to PE 0
    INTEGER(iwp) ::  mg_levels_x              !< maximum number of grid level allowed along x-direction
    INTEGER(iwp) ::  mg_levels_y              !< maximum number of grid level allowed along y-direction
    INTEGER(iwp) ::  mg_levels_z              !< maximum number of grid level allowed along z-direction
    INTEGER(iwp) ::  mg_switch_to_pe0_level_l !< maximum number of grid level with switching to PE 0
#if defined( __parallel )
    INTEGER(iwp) ::  nnx_y                    !< quotient of number of grid points along x-direction and number of PEs used along y-direction
    INTEGER(iwp) ::  nny_x                    !< quotient of number of grid points along y-direction and number of PEs used along x-direction
    INTEGER(iwp) ::  nny_z                    !< quotient of number of grid points along y-direction and number of PEs used along x-direction
    INTEGER(iwp) ::  nnz_x                    !< quotient of number of grid points along z-direction and number of PEs used along x-direction
    INTEGER(iwp) ::  nnz_y                    !< quotient of number of grid points along z-direction and number of PEs used along x-direction
    INTEGER(iwp) ::  numproc_sqr              !< square root of the number of processors
#endif
    INTEGER(iwp) ::  nxl_l                    !< lower index bound along x-direction on subdomain and different multigrid level
    INTEGER(iwp) ::  nxr_l                    !< upper index bound along x-direction on subdomain and different multigrid level
    INTEGER(iwp) ::  nyn_l                    !< lower index bound along y-direction on subdomain and different multigrid level
    INTEGER(iwp) ::  nys_l                    !< upper index bound along y-direction on subdomain and different multigrid level
#if defined( __parallel )
    INTEGER(iwp) ::  nzb_l                    !< lower index bound along z-direction on subdomain and different multigrid level
#endif
    INTEGER(iwp) ::  nzt_l                    !< upper index bound along z-direction on subdomain and different multigrid level
!$  INTEGER(iwp) ::  omp_get_num_threads      !< number of OpenMP threads

#if defined( __parallel )
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  ind_all !< dummy array containing index bounds on subdomain, used for gathering
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  nxlf    !< lower index bound allong x-direction for every PE
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  nxrf    !< upper index bound allong x-direction for every PE
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  nynf    !< lower index bound allong y-direction for every PE
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  nysf    !< lower index bound allong y-direction for every PE

    INTEGER(iwp), DIMENSION(2) ::  pdims_remote         !< number of PEs used for coupled model (only in atmospher-ocean coupling)
    INTEGER(iwp)               ::  lcoord(2)            !< PE coordinates of left neighbor along x and y
    INTEGER(iwp)               ::  rcoord(2)            !< PE coordinates of right neighbor along x and y
#endif

!
!-- Get the number of OpenMP threads
    !$OMP PARALLEL
!$  threads_per_task = omp_get_num_threads()
    !$OMP END PARALLEL


#if defined( __parallel )

    CALL location_message( 'creating virtual PE grids + MPI derived data types', 'start' )

!
!-- Determine the processor topology or check it, if prescribed by the user
    IF ( npex == -1  .AND.  npey == -1 )  THEN

!
!--    Automatic determination of the topology
       numproc_sqr = SQRT( REAL( numprocs, KIND=wp ) )
       pdims(1)    = MAX( numproc_sqr , 1 )
       DO  WHILE ( MOD( numprocs , pdims(1) ) /= 0 )
          pdims(1) = pdims(1) - 1
       ENDDO
       pdims(2) = numprocs / pdims(1)

    ELSEIF ( npex /= -1  .AND.  npey /= -1 )  THEN

!
!--    Prescribed by user. Number of processors on the prescribed topology
!--    must be equal to the number of PEs available to the job
       IF ( ( npex * npey ) /= numprocs )  THEN
          WRITE( message_string, * ) 'number of PEs of the prescribed ',       &
              'topology (', npex*npey,') does not match & the number of ',     &
              'PEs available to the job (', numprocs, ')'
          CALL message( 'init_pegrid', 'PA0221', 1, 2, 0, 6, 0 )
       ENDIF
       pdims(1) = npex
       pdims(2) = npey

    ELSE
!
!--    If the processor topology is prescribed by the user, the number of
!--    PEs must be given in both directions
       message_string = 'if the processor topology is prescribed by th' //     &
                'e user & both values of "npex" and "npey" must be given' //   &
                ' in the &NAMELIST-parameter file'
       CALL message( 'init_pegrid', 'PA0222', 1, 2, 0, 6, 0 )

    ENDIF

!
!-- Create four default MPI communicators for the 2d virtual PE grid. One of them will be used
!-- as the main communicator for this run, while others might be used for specific quantities like
!-- aerosol, chemical species, or passive scalars), if their horizontal boundary conditions shall
!-- be different from those of the other quantities (e.g. non-cyclic conditions for aerosols, and
!-- cyclic conditions for all others).
    DO  i = 1, 4

       IF ( i == 1 )  cyclic = (/  .TRUE., .TRUE.  /)   ! cyclic along x and y
       IF ( i == 2 )  cyclic = (/  .TRUE., .FALSE. /)   ! cyclic along x
       IF ( i == 3 )  cyclic = (/ .FALSE., .TRUE.  /)   ! cyllic along y
       IF ( i == 4 )  cyclic = (/ .FALSE., .FALSE. /)   ! non-cyclic

       CALL MPI_CART_CREATE( comm_palm, ndim, pdims, cyclic, reorder,                              &
                             communicator_configurations(i)%mpi_communicator, ierr )

       CALL MPI_CART_SHIFT( communicator_configurations(i)%mpi_communicator, 0, 1,                 &
                            communicator_configurations(i)%pleft,                                  &
                            communicator_configurations(i)%pright, ierr )

       CALL MPI_CART_SHIFT( communicator_configurations(i)%mpi_communicator, 1, 1,                 &
                            communicator_configurations(i)%psouth,                                 &
                            communicator_configurations(i)%pnorth, ierr )

    ENDDO

!
!-- If necessary, set horizontal boundary conditions to non-cyclic
    IF ( bc_lr /= 'cyclic' )  cyclic(1) = .FALSE.
    IF ( bc_ns /= 'cyclic' )  cyclic(2) = .FALSE.


!
!-- Set the main communicator (virtual pe grid) for this run
    IF ( bc_lr == 'cyclic'  .AND.  bc_ns == 'cyclic' )  i = 1
    IF ( bc_lr == 'cyclic'  .AND.  bc_ns /= 'cyclic' )  i = 2
    IF ( bc_lr /= 'cyclic'  .AND.  bc_ns == 'cyclic' )  i = 3
    IF ( bc_lr /= 'cyclic'  .AND.  bc_ns /= 'cyclic' )  i = 4

    comm2d = communicator_configurations(i)%mpi_communicator
    pleft  = communicator_configurations(i)%pleft
    pright = communicator_configurations(i)%pright
    psouth = communicator_configurations(i)%psouth
    pnorth = communicator_configurations(i)%pnorth

!
!-- Set rank and coordinates of the main communicator
    CALL MPI_COMM_RANK( comm2d, myid, ierr )
    WRITE (myid_char,'(''_'',I6.6)')  myid

    CALL MPI_CART_COORDS( comm2d, myid, ndim, pcoord, ierr )

!
!-- In case of cyclic boundary conditions, a y-shift at the boundaries in
!-- x-direction can be introduced via parameter y_shift. The shift is done
!-- by modifying the processor grid in such a way that processors located
!-- at the x-boundary communicate across it to processors with y-coordinate
!-- shifted by y_shift relative to their own. This feature can not be used
!-- in combination with an fft pressure solver. It has been implemented to
!-- counter the effect of streak structures in case of cyclic boundary
!-- conditions. For a description of these see Munters
!-- (2016; dx.doi.org/10.1063/1.4941912)
!--
!-- Get coordinates of left and right neighbor on PE grid
    IF ( y_shift /= 0 ) THEN
       IF ( bc_lr == 'cyclic' ) THEN
          IF ( TRIM( psolver ) /= 'multigrid' .AND.                            &
                TRIM( psolver ) /= 'multigrid_noopt')                          &
          THEN
             message_string = 'y_shift /= 0 requires a multigrid pressure solver '
             CALL message( 'check_parameters', 'PA0468', 1, 2, 0, 6, 0 )
          ENDIF

          CALL MPI_CART_COORDS( comm2d, pright, ndim, rcoord, ierr )
          CALL MPI_CART_COORDS( comm2d, pleft, ndim, lcoord, ierr )

!   
!--       If the x(y)-coordinate of the right (left) neighbor is smaller (greater)
!--       than that of the calling process, then the calling process is located on
!--       the right (left) boundary of the processor grid. In that case,
!--       the y-coordinate of that neighbor is increased (decreased) by y_shift.
!--       The rank of the process with that coordinate is then inquired and the
!--       neighbor rank for MPI_SENDRECV, pright (pleft) is set to it.
!--       In this way, the calling process receives a new right (left) neighbor
!--       for all future MPI_SENDRECV calls. That neighbor has a y-coordinate
!--       of y+(-)y_shift, where y is the original right (left) neighbor's
!--       y-coordinate. The modulo-operation ensures that if the neighbor's
!--       y-coordinate exceeds the grid-boundary, it will be relocated to
!--       the opposite part of the grid cyclicly.
          IF ( rcoord(1) < pcoord(1) ) THEN
             rcoord(2) = MODULO( rcoord(2) + y_shift, pdims(2) )
             CALL MPI_CART_RANK( comm2d, rcoord, pright, ierr )
          ENDIF

          IF ( lcoord(1) > pcoord(1) ) THEN
             lcoord(2) = MODULO( lcoord(2) - y_shift, pdims(2) )
             CALL MPI_CART_RANK( comm2d, lcoord, pleft, ierr )
          ENDIF
         
       ELSE
!
!--       y-shift for non-cyclic boundary conditions is only implemented 
!--       for the turbulence recycling method in inflow_turbulence.f90
          IF ( .NOT. turbulent_inflow )  THEN
             message_string = 'y_shift /= 0 is only allowed for cyclic ' //    &
                              'boundary conditions in both directions '  //    &
                              'or with turbulent_inflow == .TRUE.'
             CALL message( 'check_parameters', 'PA0467', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF
    ENDIF

!
!-- Vertical nesting: store four lists that identify partner ranks to exchange 
!-- data
    IF ( vnested )  CALL vnest_init_pegrid_rank

!
!-- Determine sub-topologies for transpositions
!-- Transposition from z to x:
    remain_dims(1) = .TRUE.
    remain_dims(2) = .FALSE.
    CALL MPI_CART_SUB( comm2d, remain_dims, comm1dx, ierr )
    CALL MPI_COMM_RANK( comm1dx, myidx, ierr )
!
!-- Transposition from x to y
    remain_dims(1) = .FALSE.
    remain_dims(2) = .TRUE.
    CALL MPI_CART_SUB( comm2d, remain_dims, comm1dy, ierr )
    CALL MPI_COMM_RANK( comm1dy, myidy, ierr )


!
!-- Calculate array bounds along x-direction for every PE.
    ALLOCATE( nxlf(0:pdims(1)-1), nxrf(0:pdims(1)-1), nynf(0:pdims(2)-1),      &
              nysf(0:pdims(2)-1) )

    IF ( MOD( nx+1 , pdims(1) ) /= 0 )  THEN
       WRITE( message_string, * ) 'x-direction: gridpoint number (',nx+1,') ', &
                                  'is not an& integral multiple of the number',&
                                  ' of processors (',pdims(1),')'
       CALL message( 'init_pegrid', 'PA0225', 1, 2, 0, 6, 0 )
    ELSE
       nnx  = ( nx + 1 ) / pdims(1)
    ENDIF

!
!-- Left and right array bounds, number of gridpoints
    DO  i = 0, pdims(1)-1
       nxlf(i)   = i * nnx
       nxrf(i)   = ( i + 1 ) * nnx - 1
    ENDDO

!
!-- Calculate array bounds in y-direction for every PE.
    IF ( MOD( ny+1 , pdims(2) ) /= 0 )  THEN
       WRITE( message_string, * ) 'y-direction: gridpoint number (',ny+1,') ', &
                                  'is not an& integral multiple of the number',&
                                  ' of processors (',pdims(2),')'
       CALL message( 'init_pegrid', 'PA0227', 1, 2, 0, 6, 0 )
    ELSE
       nny  = ( ny + 1 ) / pdims(2)
    ENDIF

!
!-- South and north array bounds
    DO  j = 0, pdims(2)-1
       nysf(j)   = j * nny
       nynf(j)   = ( j + 1 ) * nny - 1
    ENDDO

!
!-- Local array bounds of the respective PEs
    nxl = nxlf(pcoord(1))
    nxr = nxrf(pcoord(1))
    nys = nysf(pcoord(2))
    nyn = nynf(pcoord(2))
    nzb = 0
    nzt = nz
    nnz = nz

!
!-- Set switches to define if the PE is situated at the border of the virtual
!-- processor grid
    IF ( nxl == 0 )   left_border_pe  = .TRUE.
    IF ( nxr == nx )  right_border_pe = .TRUE.
    IF ( nys == 0 )   south_border_pe = .TRUE.
    IF ( nyn == ny )  north_border_pe = .TRUE.

!
!-- Calculate array bounds and gridpoint numbers for the transposed arrays
!-- (needed in the pressure solver)
!-- For the transposed arrays, cyclic boundaries as well as top and bottom
!-- boundaries are omitted, because they are obstructive to the transposition

!
!-- 1. transposition  z --> x
!-- This transposition is not neccessary in case of a 1d-decomposition along x
    IF ( psolver == 'poisfft'  .OR.  calculate_spectra )  THEN

       IF ( pdims(2) /= 1 )  THEN
          IF ( MOD( nz , pdims(1) ) /= 0 )  THEN
             WRITE( message_string, * ) 'transposition z --> x:& ',              &
                                        'nz=',nz,' is not an integral multiple ',&
                                        'of pdims(1)=',pdims(1)
             CALL message( 'init_pegrid', 'PA0230', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF

       nys_x = nys
       nyn_x = nyn
       nny_x = nny
       nnz_x = nz / pdims(1)
       nzb_x = 1 + myidx * nnz_x
       nzt_x = ( myidx + 1 ) * nnz_x
       sendrecvcount_zx = nnx * nny * nnz_x

    ENDIF


    IF ( psolver == 'poisfft' )  THEN
!
!--    2. transposition  x --> y
       IF ( MOD( nx+1 , pdims(2) ) /= 0 )  THEN
          WRITE( message_string, * ) 'transposition x --> y:& ',              &
                                     'nx+1=',nx+1,' is not an integral ',     &
                                     'multiple of pdims(2)=',pdims(2)
          CALL message( 'init_pegrid', 'PA0231', 1, 2, 0, 6, 0 )
       ENDIF

       nnz_y = nnz_x
       nzb_y = nzb_x
       nzt_y = nzt_x
       nnx_y = (nx+1) / pdims(2)
       nxl_y = myidy * nnx_y
       nxr_y = ( myidy + 1 ) * nnx_y - 1
       sendrecvcount_xy = nnx_y * nny_x * nnz_y
!
!--    3. transposition  y --> z
!--    (ELSE:  x --> y  in case of 1D-decomposition along x)
       nxl_z = nxl_y
       nxr_z = nxr_y
       nny_z = (ny+1) / pdims(1)
       nys_z = myidx * nny_z
       nyn_z = ( myidx + 1 ) * nny_z - 1
       sendrecvcount_yz = nnx_y * nny_z * nnz_y

       IF ( pdims(2) /= 1 )  THEN
!
!--       y --> z
!--       This transposition is not neccessary in case of a 1d-decomposition
!--       along x, except that the uptream-spline method is switched on
          IF ( MOD( ny+1 , pdims(1) ) /= 0 )  THEN
             WRITE( message_string, * ) 'transposition y --> z:& ',            &
                                        'ny+1=',ny+1,' is not an integral ',   &
                                        'multiple of pdims(1)=',pdims(1)
             CALL message( 'init_pegrid', 'PA0232', 1, 2, 0, 6, 0 )
          ENDIF

       ELSE
!
!--       x --> y
!--       This condition must be fulfilled for a 1D-decomposition along x
          IF ( MOD( ny+1 , pdims(1) ) /= 0 )  THEN
             WRITE( message_string, * ) 'transposition x --> y:& ',            &
                                        'ny+1=',ny+1,' is not an integral ',   &
                                        'multiple of pdims(1)=',pdims(1)
             CALL message( 'init_pegrid', 'PA0233', 1, 2, 0, 6, 0 )
          ENDIF

       ENDIF

    ENDIF

!
!-- Indices for direct transpositions z --> y (used for calculating spectra)
    IF ( calculate_spectra )  THEN
       IF ( MOD( nz, pdims(2) ) /= 0 )  THEN
          WRITE( message_string, * ) 'direct transposition z --> y (needed ',  &
                                     'for spectra):& nz=',nz,' is not an ',    &
                                     'integral multiple of pdims(2)=',pdims(2)
          CALL message( 'init_pegrid', 'PA0234', 1, 2, 0, 6, 0 )
       ELSE
          nxl_yd = nxl
          nxr_yd = nxr
          nzb_yd = 1 + myidy * ( nz / pdims(2) )
          nzt_yd = ( myidy + 1 ) * ( nz / pdims(2) )
          sendrecvcount_zyd = nnx * nny * ( nz / pdims(2) )
       ENDIF
    ENDIF

    IF ( psolver == 'poisfft'  .OR.  calculate_spectra )  THEN
!
!--    Indices for direct transpositions y --> x 
!--    (they are only possible in case of a 1d-decomposition along x)
       IF ( pdims(2) == 1 )  THEN
          nny_x = nny / pdims(1)
          nys_x = myid * nny_x
          nyn_x = ( myid + 1 ) * nny_x - 1
          nzb_x = 1
          nzt_x = nz
          sendrecvcount_xy = nnx * nny_x * nz
       ENDIF

    ENDIF

    IF ( psolver == 'poisfft' )  THEN
!
!--    Indices for direct transpositions x --> y 
!--    (they are only possible in case of a 1d-decomposition along y)
       IF ( pdims(1) == 1 )  THEN
          nnx_y = nnx / pdims(2)
          nxl_y = myid * nnx_y
          nxr_y = ( myid + 1 ) * nnx_y - 1
          nzb_y = 1
          nzt_y = nz
          sendrecvcount_xy = nnx_y * nny * nz
       ENDIF

    ENDIF

!
!-- Arrays for storing the array bounds are needed any more
    DEALLOCATE( nxlf , nxrf , nynf , nysf )


!
!-- Collect index bounds from other PEs (to be written to restart file later)
    ALLOCATE( hor_index_bounds(4,0:numprocs-1) )

    IF ( myid == 0 )  THEN

       hor_index_bounds(1,0) = nxl
       hor_index_bounds(2,0) = nxr
       hor_index_bounds(3,0) = nys
       hor_index_bounds(4,0) = nyn

!
!--    Receive data from all other PEs
       DO  i = 1, numprocs-1
          CALL MPI_RECV( ibuf, 4, MPI_INTEGER, i, MPI_ANY_TAG, comm2d, status, &
                         ierr )
          hor_index_bounds(:,i) = ibuf(1:4)
       ENDDO

    ELSE
!
!--    Send index bounds to PE0
       ibuf(1) = nxl
       ibuf(2) = nxr
       ibuf(3) = nys
       ibuf(4) = nyn
       CALL MPI_SEND( ibuf, 4, MPI_INTEGER, 0, myid, comm2d, ierr )

    ENDIF


#if defined( __print )
!
!-- Control output
    IF ( myid == 0 )  THEN
       PRINT*, '*** processor topology ***'
       PRINT*, ' '
       PRINT*, 'myid   pcoord    left right  south north  idx idy   nxl: nxr',&
               &'   nys: nyn'
       PRINT*, '------------------------------------------------------------',&
               &'-----------'
       WRITE (*,1000)  0, pcoord(1), pcoord(2), pleft, pright, psouth, pnorth, &
                       myidx, myidy, nxl, nxr, nys, nyn
1000   FORMAT (I4,2X,'(',I3,',',I3,')',3X,I4,2X,I4,3X,I4,2X,I4,2X,I3,1X,I3, &
               2(2X,I4,':',I4))

!
!--    Receive data from the other PEs
       DO  i = 1,numprocs-1
          CALL MPI_RECV( ibuf, 12, MPI_INTEGER, i, MPI_ANY_TAG, comm2d, status, &
                         ierr )
          WRITE (*,1000)  i, ( ibuf(j) , j = 1,12 )
       ENDDO
    ELSE

!
!--    Send data to PE0
       ibuf(1) = pcoord(1); ibuf(2) = pcoord(2); ibuf(3) = pleft
       ibuf(4) = pright; ibuf(5) = psouth; ibuf(6) = pnorth; ibuf(7) = myidx
       ibuf(8) = myidy; ibuf(9) = nxl; ibuf(10) = nxr; ibuf(11) = nys
       ibuf(12) = nyn
       CALL MPI_SEND( ibuf, 12, MPI_INTEGER, 0, myid, comm2d, ierr )       
    ENDIF
#endif

! 
!-- Determine the number of ghost point layers
    IF ( scalar_advec   == 'ws-scheme'  .OR.                                   &
         momentum_advec == 'ws-scheme'  .OR.  nested_run )  THEN
       nbgp = 3
    ELSE
       nbgp = 1
    ENDIF 

!
!-- Check that the number of computational grid points is not smaller than the number of
!-- ghost points.
    IF ( nnx < nbgp )  THEN
       WRITE( message_string, * ) 'number of subdomain grid points along x (', nnx, ') is smaller',&
                                  'than the number of ghost points (', nbgp, ')'
       CALL message( 'init_pegrid', 'PA0682', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( nny < nbgp )  THEN
       WRITE( message_string, * ) 'number of subdomain grid points along y (', nny, ') is smaller',&
                                  'than the number of ghost points (', nbgp, ')'
       CALL message( 'init_pegrid', 'PA0683', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Create a new MPI derived datatype for the exchange of surface (xy) data,
!-- which is needed for coupled atmosphere-ocean runs.
!-- First, calculate number of grid points of an xy-plane.
    ngp_xy  = ( nxr - nxl + 1 + 2 * nbgp ) * ( nyn - nys + 1 + 2 * nbgp )
    CALL MPI_TYPE_VECTOR( ngp_xy, 1, nzt-nzb+2, MPI_REAL, type_xy, ierr )
    CALL MPI_TYPE_COMMIT( type_xy, ierr )

    IF ( TRIM( coupling_mode ) /= 'uncoupled' .AND. .NOT. vnested )  THEN
   
!
!--    Pass the number of grid points of the atmosphere model to
!--    the ocean model and vice versa
       IF ( coupling_mode == 'atmosphere_to_ocean' )  THEN

          nx_a = nx
          ny_a = ny

          IF ( myid == 0 )  THEN

             CALL MPI_SEND( nx_a, 1, MPI_INTEGER, numprocs, 1, comm_inter,  &
                            ierr )
             CALL MPI_SEND( ny_a, 1, MPI_INTEGER, numprocs, 2, comm_inter,  &
                            ierr )
             CALL MPI_SEND( pdims, 2, MPI_INTEGER, numprocs, 3, comm_inter, &
                            ierr )
             CALL MPI_RECV( nx_o, 1, MPI_INTEGER, numprocs, 4, comm_inter,  &
                            status, ierr )
             CALL MPI_RECV( ny_o, 1, MPI_INTEGER, numprocs, 5, comm_inter,  &
                            status, ierr )
             CALL MPI_RECV( pdims_remote, 2, MPI_INTEGER, numprocs, 6,      &
                            comm_inter, status, ierr )
          ENDIF

          CALL MPI_BCAST( nx_o, 1, MPI_INTEGER, 0, comm2d, ierr )
          CALL MPI_BCAST( ny_o, 1, MPI_INTEGER, 0, comm2d, ierr ) 
          CALL MPI_BCAST( pdims_remote, 2, MPI_INTEGER, 0, comm2d, ierr )
       
       ELSEIF ( coupling_mode == 'ocean_to_atmosphere' )  THEN

          nx_o = nx
          ny_o = ny 

          IF ( myid == 0 ) THEN

             CALL MPI_RECV( nx_a, 1, MPI_INTEGER, 0, 1, comm_inter, status, &
                            ierr )
             CALL MPI_RECV( ny_a, 1, MPI_INTEGER, 0, 2, comm_inter, status, &
                            ierr )
             CALL MPI_RECV( pdims_remote, 2, MPI_INTEGER, 0, 3, comm_inter, &
                            status, ierr )
             CALL MPI_SEND( nx_o, 1, MPI_INTEGER, 0, 4, comm_inter, ierr )
             CALL MPI_SEND( ny_o, 1, MPI_INTEGER, 0, 5, comm_inter, ierr )
             CALL MPI_SEND( pdims, 2, MPI_INTEGER, 0, 6, comm_inter, ierr )
          ENDIF

          CALL MPI_BCAST( nx_a, 1, MPI_INTEGER, 0, comm2d, ierr)
          CALL MPI_BCAST( ny_a, 1, MPI_INTEGER, 0, comm2d, ierr) 
          CALL MPI_BCAST( pdims_remote, 2, MPI_INTEGER, 0, comm2d, ierr) 

       ENDIF
  
       ngp_a = ( nx_a+1 + 2 * nbgp ) * ( ny_a+1 + 2 * nbgp )
       ngp_o = ( nx_o+1 + 2 * nbgp ) * ( ny_o+1 + 2 * nbgp )

!
!--    Determine if the horizontal grid and the number of PEs in ocean and
!--    atmosphere is same or not
       IF ( nx_o == nx_a  .AND.  ny_o == ny_a  .AND.  &
            pdims(1) == pdims_remote(1) .AND. pdims(2) == pdims_remote(2) ) &
       THEN
          coupling_topology = 0
       ELSE
          coupling_topology = 1
       ENDIF 

!
!--    Determine the target PEs for the exchange between ocean and
!--    atmosphere (comm2d)
       IF ( coupling_topology == 0 )  THEN
!
!--       In case of identical topologies, every atmosphere PE has exactly one
!--       ocean PE counterpart and vice versa
          IF ( TRIM( coupling_mode ) == 'atmosphere_to_ocean' ) THEN
             target_id = myid + numprocs
          ELSE
             target_id = myid 
          ENDIF

       ELSE
!
!--       In case of nonequivalent topology in ocean and atmosphere only for
!--       PE0 in ocean and PE0 in atmosphere a target_id is needed, since
!--       data echxchange between ocean and atmosphere will be done only
!--       between these PEs.    
          IF ( myid == 0 )  THEN

             IF ( TRIM( coupling_mode ) == 'atmosphere_to_ocean' )  THEN
                target_id = numprocs 
             ELSE
                target_id = 0
             ENDIF

          ENDIF

       ENDIF

    ENDIF

!
!-- Store partner grid point co-ordinates as lists.
!-- Create custom MPI vector datatypes for contiguous data transfer
    IF ( vnested )  CALL vnest_init_pegrid_domain

#else

!
!-- Array bounds when running on a single PE (respectively a non-parallel
!-- machine)
    nxl = 0
    nxr = nx
    nnx = nxr - nxl + 1
    nys = 0
    nyn = ny
    nny = nyn - nys + 1
    nzb = 0
    nzt = nz
    nnz = nz

    ALLOCATE( hor_index_bounds(4,0:0) )
    hor_index_bounds(1,0) = nxl
    hor_index_bounds(2,0) = nxr
    hor_index_bounds(3,0) = nys
    hor_index_bounds(4,0) = nyn

!
!-- Array bounds for the pressure solver (in the parallel code, these bounds
!-- are the ones for the transposed arrays)
    nys_x = nys
    nyn_x = nyn
    nzb_x = nzb + 1
    nzt_x = nzt

    nxl_y = nxl
    nxr_y = nxr
    nzb_y = nzb + 1
    nzt_y = nzt

    nxl_z = nxl
    nxr_z = nxr
    nys_z = nys
    nyn_z = nyn

#endif

!
!-- Calculate number of grid levels necessary for the multigrid poisson solver
!-- as well as the gridpoint indices on each level
    IF ( psolver(1:9) == 'multigrid' )  THEN

!
!--    First calculate number of possible grid levels for the subdomains
       mg_levels_x = 1
       mg_levels_y = 1
       mg_levels_z = 1

       i = nnx
       DO WHILE ( MOD( i, 2 ) == 0  .AND.  i /= 2 )
          i = i / 2
          mg_levels_x = mg_levels_x + 1
       ENDDO

       j = nny
       DO WHILE ( MOD( j, 2 ) == 0  .AND.  j /= 2 )
          j = j / 2
          mg_levels_y = mg_levels_y + 1
       ENDDO

       k = nz    ! do not use nnz because it might be > nz due to transposition
                 ! requirements
       DO WHILE ( MOD( k, 2 ) == 0  .AND.  k /= 2 )
          k = k / 2
          mg_levels_z = mg_levels_z + 1
       ENDDO
!
!--    The optimized MG-solver does not allow odd values for nz at the coarsest
!--    grid level
       IF ( TRIM( psolver ) /= 'multigrid_noopt' )  THEN
          IF ( MOD( k, 2 ) /= 0 )  mg_levels_z = mg_levels_z - 1
!
!--       An odd value of nz does not work. The finest level must have an even
!--       value.
          IF (  mg_levels_z == 0 )  THEN
             message_string = 'optimized multigrid method requires nz to be even'
             CALL message( 'init_pegrid', 'PA0495', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF

       maximum_grid_level = MIN( mg_levels_x, mg_levels_y, mg_levels_z )
!       
!--    Check if subdomain sizes prevents any coarsening. 
!--    This case, the maximum number of grid levels is 1, i.e. effectively 
!--    a Gauss-Seidel scheme is applied rather than a multigrid approach.
!--    Give a warning in this case.
       IF ( maximum_grid_level == 1  .AND.  mg_switch_to_pe0_level == -1 )  THEN
          message_string = 'No grid coarsening possible, multigrid ' //        &
                           'approach effectively reduces to a Gauss-Seidel ' //&
                           'scheme.'
          
          CALL message( 'poismg', 'PA0648', 0, 1, 0, 6, 0 )
       ENDIF

!
!--    Find out, if the total domain allows more levels. These additional
!--    levels are identically processed on all PEs.
       IF ( numprocs > 1  .AND.  mg_switch_to_pe0_level /= -1 )  THEN

          IF ( mg_levels_z > MIN( mg_levels_x, mg_levels_y ) )  THEN

             mg_switch_to_pe0_level_l = maximum_grid_level

             mg_levels_x = 1
             mg_levels_y = 1

             i = nx+1
             DO WHILE ( MOD( i, 2 ) == 0  .AND.  i /= 2 )
                i = i / 2
                mg_levels_x = mg_levels_x + 1
             ENDDO

             j = ny+1
             DO WHILE ( MOD( j, 2 ) == 0  .AND.  j /= 2 )
                j = j / 2
                mg_levels_y = mg_levels_y + 1
             ENDDO

             maximum_grid_level_l = MIN( mg_levels_x, mg_levels_y, mg_levels_z )

             IF ( maximum_grid_level_l > mg_switch_to_pe0_level_l )  THEN
                mg_switch_to_pe0_level_l = maximum_grid_level_l - &
                                           mg_switch_to_pe0_level_l + 1
             ELSE
                mg_switch_to_pe0_level_l = 0
             ENDIF

          ELSE

             mg_switch_to_pe0_level_l = 0
             maximum_grid_level_l = maximum_grid_level

          ENDIF

!
!--       Use switch level calculated above only if it is not pre-defined
!--       by user
          IF ( mg_switch_to_pe0_level == 0 )  THEN
             IF ( mg_switch_to_pe0_level_l /= 0 )  THEN
                mg_switch_to_pe0_level = mg_switch_to_pe0_level_l
                maximum_grid_level     = maximum_grid_level_l
             ENDIF

          ELSE
!
!--          Check pre-defined value and reset to default, if neccessary
             IF ( mg_switch_to_pe0_level < mg_switch_to_pe0_level_l  .OR.      &
                  mg_switch_to_pe0_level >= maximum_grid_level_l )  THEN
                message_string = 'mg_switch_to_pe0_level ' //                  &
                                 'out of range and reset to 0'
                CALL message( 'init_pegrid', 'PA0235', 0, 1, 0, 6, 0 )
                mg_switch_to_pe0_level = 0
             ELSE
!
!--             Use the largest number of possible levels anyway and recalculate
!--             the switch level to this largest number of possible values
                maximum_grid_level = maximum_grid_level_l

             ENDIF

          ENDIF

       ENDIF

       ALLOCATE( grid_level_count(maximum_grid_level),                       &
                 nxl_mg(0:maximum_grid_level), nxr_mg(0:maximum_grid_level), &
                 nyn_mg(0:maximum_grid_level), nys_mg(0:maximum_grid_level), &
                 nzt_mg(0:maximum_grid_level) )

       grid_level_count = 0
!
!--    Index zero required as dummy due to definition of arrays f2 and p2 in
!--    recursive subroutine next_mg_level
       nxl_mg(0) = 0; nxr_mg(0) = 0; nyn_mg(0) = 0; nys_mg(0) = 0; nzt_mg(0) = 0

       nxl_l = nxl; nxr_l = nxr; nys_l = nys; nyn_l = nyn; nzt_l = nzt

       DO  i = maximum_grid_level, 1 , -1

          IF ( i == mg_switch_to_pe0_level )  THEN
#if defined( __parallel )
!
!--          Save the grid size of the subdomain at the switch level, because
!--          it is needed in poismg.
             ind(1) = nxl_l; ind(2) = nxr_l
             ind(3) = nys_l; ind(4) = nyn_l
             ind(5) = nzt_l
             ALLOCATE( ind_all(5*numprocs), mg_loc_ind(5,0:numprocs-1) )
             CALL MPI_ALLGATHER( ind, 5, MPI_INTEGER, ind_all, 5, &
                                 MPI_INTEGER, comm2d, ierr )
             DO  j = 0, numprocs-1
                DO  k = 1, 5
                   mg_loc_ind(k,j) = ind_all(k+j*5)
                ENDDO
             ENDDO
             DEALLOCATE( ind_all )
!
!--          Calculate the grid size of the total domain
             nxr_l = ( nxr_l-nxl_l+1 ) * pdims(1) - 1
             nxl_l = 0
             nyn_l = ( nyn_l-nys_l+1 ) * pdims(2) - 1
             nys_l = 0
!
!--          The size of this gathered array must not be larger than the
!--          array tend, which is used in the multigrid scheme as a temporary
!--          array. Therefore the subdomain size of an PE is calculated and 
!--          the size of the gathered grid. These values are used in  
!--          routines pres and poismg
             subdomain_size = ( nxr - nxl + 2 * nbgp + 1 ) * &
                              ( nyn - nys + 2 * nbgp + 1 ) * ( nzt - nzb + 2 )
             gathered_size  = ( nxr_l - nxl_l + 3 ) * ( nyn_l - nys_l + 3 ) *  &
                              ( nzt_l - nzb + 2 )

#else
             message_string = 'multigrid gather/scatter impossible ' //        &
                          'in non parallel mode'
             CALL message( 'init_pegrid', 'PA0237', 1, 2, 0, 6, 0 )
#endif
          ENDIF

          nxl_mg(i) = nxl_l
          nxr_mg(i) = nxr_l
          nys_mg(i) = nys_l
          nyn_mg(i) = nyn_l
          nzt_mg(i) = nzt_l

          nxl_l = nxl_l / 2 
          nxr_l = nxr_l / 2
          nys_l = nys_l / 2 
          nyn_l = nyn_l / 2 
          nzt_l = nzt_l / 2 
          
       ENDDO

!
!--    Temporary problem: Currently calculation of maxerror in routine poismg crashes
!--    if grid data are collected on PE0 already on the finest grid level.
!--    To be solved later.
       IF ( maximum_grid_level == mg_switch_to_pe0_level )  THEN
          message_string = 'grid coarsening on subdomain level cannot be performed'
          CALL message( 'poismg', 'PA0236', 1, 2, 0, 6, 0 )
       ENDIF

    ELSE

       maximum_grid_level = 0

    ENDIF

!
!-- Default level 0 tells exchange_horiz that all ghost planes have to be
!-- exchanged. grid_level is adjusted in poismg, where only one ghost plane
!-- is required.
    grid_level = 0

#if defined( __parallel )
!
!-- Gridpoint number for the exchange of ghost points (y-line for 2D-arrays)
    ngp_y  = nyn - nys + 1 + 2 * nbgp

!
!-- Define new MPI derived datatypes for the exchange of ghost points in
!-- x- and y-direction for 2D-arrays (line)
    CALL MPI_TYPE_VECTOR( nxr-nxl+1+2*nbgp, nbgp, ngp_y, MPI_REAL, type_x,     &
                          ierr )
    CALL MPI_TYPE_COMMIT( type_x, ierr )

    CALL MPI_TYPE_VECTOR( nbgp, ngp_y, ngp_y, MPI_REAL, type_y, ierr )
    CALL MPI_TYPE_COMMIT( type_y, ierr )
!
!-- Define new MPI derived datatypes for the exchange of ghost points in
!-- x- and y-direction for 2D-INTEGER arrays (line) - on normal grid.
!-- Define types for 32-bit and 8-bit Integer. The 8-bit Integer are only 
!-- required on normal grid, while 32-bit Integer may be also required on 
!-- coarser grid level in case of multigrid solver.
!
!-- 8-bit Integer
    CALL MPI_TYPE_VECTOR( nxr-nxl+1+2*nbgp, nbgp, ngp_y, MPI_BYTE,             &
                          type_x_byte, ierr )
    CALL MPI_TYPE_COMMIT( type_x_byte, ierr )

    CALL MPI_TYPE_VECTOR( nbgp, ngp_y, ngp_y, MPI_BYTE,                        &
                          type_y_byte, ierr )
    CALL MPI_TYPE_COMMIT( type_y_byte, ierr )
!
!-- 32-bit Integer
    ALLOCATE( type_x_int(0:maximum_grid_level),                                &
              type_y_int(0:maximum_grid_level) )
              
    CALL MPI_TYPE_VECTOR( nxr-nxl+1+2*nbgp, nbgp, ngp_y, MPI_INTEGER,          &
                          type_x_int(0), ierr )
    CALL MPI_TYPE_COMMIT( type_x_int(0), ierr )

    CALL MPI_TYPE_VECTOR( nbgp, ngp_y, ngp_y, MPI_INTEGER, type_y_int(0), ierr )
    CALL MPI_TYPE_COMMIT( type_y_int(0), ierr )
!
!-- Calculate gridpoint numbers for the exchange of ghost points along x
!-- (yz-plane for 3D-arrays) and define MPI derived data type(s) for the
!-- exchange of ghost points in y-direction (xz-plane).
!-- Do these calculations for the model grid and (if necessary) also
!-- for the coarser grid levels used in the multigrid method
    ALLOCATE ( ngp_xz(0:maximum_grid_level),                                   &
               ngp_xz_int(0:maximum_grid_level),                               &
               ngp_yz(0:maximum_grid_level),                                   &
               ngp_yz_int(0:maximum_grid_level),                               &
               type_xz(0:maximum_grid_level),                                  &
               type_xz_int(0:maximum_grid_level),                              &
               type_yz(0:maximum_grid_level),                                  &
               type_yz_int(0:maximum_grid_level) )

    nxl_l = nxl; nxr_l = nxr; nys_l = nys; nyn_l = nyn; nzb_l = nzb; nzt_l = nzt

!
!-- Discern between the model grid, which needs nbgp ghost points and
!-- grid levels for the multigrid scheme. In the latter case only one
!-- ghost point is necessary.
!-- First definition of MPI-datatypes for exchange of ghost layers on normal 
!-- grid. The following loop is needed for data exchange in poismg.f90.
!
!-- Determine number of grid points of yz-layer for exchange
    ngp_yz(0) = (nzt - nzb + 2) * (nyn - nys + 1 + 2 * nbgp)

!
!-- Define an MPI-datatype for the exchange of left/right boundaries.
!-- Although data are contiguous in physical memory (which does not
!-- necessarily require an MPI-derived datatype), the data exchange between
!-- left and right PE's using the MPI-derived type is 10% faster than without.
    CALL MPI_TYPE_VECTOR( nxr-nxl+1+2*nbgp, nbgp*(nzt-nzb+2), ngp_yz(0), &
                          MPI_REAL, type_xz(0), ierr )
    CALL MPI_TYPE_COMMIT( type_xz(0), ierr )

    CALL MPI_TYPE_VECTOR( nbgp, ngp_yz(0), ngp_yz(0), MPI_REAL, type_yz(0), &
                          ierr ) 
    CALL MPI_TYPE_COMMIT( type_yz(0), ierr )

!
!-- Define data types for exchange of 3D Integer arrays.
    ngp_yz_int(0) = (nzt - nzb + 2) * (nyn - nys + 1 + 2 * nbgp)

    CALL MPI_TYPE_VECTOR( nxr-nxl+1+2*nbgp, nbgp*(nzt-nzb+2), ngp_yz_int(0),   &
                          MPI_INTEGER, type_xz_int(0), ierr )
    CALL MPI_TYPE_COMMIT( type_xz_int(0), ierr )

    CALL MPI_TYPE_VECTOR( nbgp, ngp_yz_int(0), ngp_yz_int(0), MPI_INTEGER,     &
                          type_yz_int(0), ierr )
    CALL MPI_TYPE_COMMIT( type_yz_int(0), ierr )

!
!-- Definition of MPI-datatypes for multigrid method (coarser level grids)
    IF ( psolver(1:9) == 'multigrid' )  THEN
!    
!--    Definition of MPI-datatyoe as above, but only 1 ghost level is used
       DO  i = maximum_grid_level, 1 , -1
!
!--       For 3D-exchange on different multigrid level, one ghost point for 
!--       REAL arrays, two ghost points for INTEGER arrays
          ngp_xz(i) = (nzt_l - nzb_l + 2) * (nxr_l - nxl_l + 3)
          ngp_yz(i) = (nzt_l - nzb_l + 2) * (nyn_l - nys_l + 3)

          ngp_xz_int(i) = (nzt_l - nzb_l + 2) * (nxr_l - nxl_l + 3)
          ngp_yz_int(i) = (nzt_l - nzb_l + 2) * (nyn_l - nys_l + 3)
!
!--       MPI data type for REAL arrays, for xz-layers
          CALL MPI_TYPE_VECTOR( nxr_l-nxl_l+3, nzt_l-nzb_l+2, ngp_yz(i),       &
                                MPI_REAL, type_xz(i), ierr )
          CALL MPI_TYPE_COMMIT( type_xz(i), ierr )

!
!--       MPI data type for INTEGER arrays, for xz-layers
          CALL MPI_TYPE_VECTOR( nxr_l-nxl_l+3, nzt_l-nzb_l+2, ngp_yz_int(i),   &
                                MPI_INTEGER, type_xz_int(i), ierr )
          CALL MPI_TYPE_COMMIT( type_xz_int(i), ierr )

!
!--       MPI data type for REAL arrays, for yz-layers
          CALL MPI_TYPE_VECTOR( 1, ngp_yz(i), ngp_yz(i), MPI_REAL, type_yz(i), &
                                ierr )
          CALL MPI_TYPE_COMMIT( type_yz(i), ierr )
!
!--       MPI data type for INTEGER arrays, for yz-layers
          CALL MPI_TYPE_VECTOR( 1, ngp_yz_int(i), ngp_yz_int(i), MPI_INTEGER,  &
                                type_yz_int(i), ierr )
          CALL MPI_TYPE_COMMIT( type_yz_int(i), ierr )


!--       For 2D-exchange of INTEGER arrays on coarser grid level, where 2 ghost
!--       points need to be exchanged. Only required for 32-bit Integer arrays.
          CALL MPI_TYPE_VECTOR( nxr_l-nxl_l+5, 2, nyn_l-nys_l+5, MPI_INTEGER,  &
                                type_x_int(i), ierr )
          CALL MPI_TYPE_COMMIT( type_x_int(i), ierr )


          CALL MPI_TYPE_VECTOR( 2, nyn_l-nys_l+5, nyn_l-nys_l+5, MPI_INTEGER,  &
                                type_y_int(i), ierr )
          CALL MPI_TYPE_COMMIT( type_y_int(i), ierr )

          nxl_l = nxl_l / 2
          nxr_l = nxr_l / 2
          nys_l = nys_l / 2
          nyn_l = nyn_l / 2
          nzt_l = nzt_l / 2

       ENDDO

    ENDIF

#endif

#if defined( __parallel )
!
!-- Setting of flags for inflow/outflow/nesting conditions.
    IF ( pleft == MPI_PROC_NULL )  THEN
       IF ( bc_lr == 'dirichlet/radiation'  .OR.  bc_lr == 'nested'  .OR.      &
            bc_lr == 'nesting_offline' )  THEN
          bc_dirichlet_l  = .TRUE.
       ELSEIF ( bc_lr == 'radiation/dirichlet' )  THEN
          bc_radiation_l = .TRUE.
       ENDIF
    ENDIF
 
    IF ( pright == MPI_PROC_NULL )  THEN
       IF ( bc_lr == 'dirichlet/radiation' )  THEN
          bc_radiation_r = .TRUE.
       ELSEIF ( bc_lr == 'radiation/dirichlet'  .OR.  bc_lr == 'nested'  .OR.  &
                bc_lr == 'nesting_offline' )  THEN
          bc_dirichlet_r  = .TRUE.
       ENDIF
    ENDIF

    IF ( psouth == MPI_PROC_NULL )  THEN
       IF ( bc_ns == 'dirichlet/radiation' )  THEN
          bc_radiation_s = .TRUE.
       ELSEIF ( bc_ns == 'radiation/dirichlet'  .OR.  bc_ns == 'nested'  .OR.  &
                bc_ns == 'nesting_offline' )  THEN
          bc_dirichlet_s  = .TRUE.
       ENDIF
    ENDIF

    IF ( pnorth == MPI_PROC_NULL )  THEN
       IF ( bc_ns == 'dirichlet/radiation'  .OR.  bc_ns == 'nested'  .OR.      &
            bc_ns == 'nesting_offline' )  THEN
          bc_dirichlet_n  = .TRUE.
       ELSEIF ( bc_ns == 'radiation/dirichlet' )  THEN
          bc_radiation_n = .TRUE.
       ENDIF
    ENDIF
!
!-- In case of synthetic turbulence geneartor determine ids. 
!-- Please note, if no forcing or nesting is applied, the generator is applied
!-- only at the left lateral boundary.
    IF ( use_syn_turb_gen )  THEN
       IF ( bc_dirichlet_l )  THEN
          id_stg_left_l = myidx
       ELSE
          id_stg_left_l = 0
       ENDIF
       IF ( bc_dirichlet_r )  THEN
          id_stg_right_l = myidx
       ELSE
          id_stg_right_l = 0
       ENDIF
       IF ( bc_dirichlet_s )  THEN
          id_stg_south_l = myidy
       ELSE
          id_stg_south_l = 0
       ENDIF
       IF ( bc_dirichlet_n )  THEN
          id_stg_north_l = myidy
       ELSE
          id_stg_north_l = 0
       ENDIF

       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( id_stg_left_l, id_stg_left,   1, MPI_INTEGER,       &
                           MPI_SUM, comm1dx, ierr )

       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( id_stg_right_l, id_stg_right, 1, MPI_INTEGER,       &
                           MPI_SUM, comm1dx, ierr )

       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( id_stg_south_l, id_stg_south, 1, MPI_INTEGER,       &
                           MPI_SUM, comm1dy, ierr )

       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( id_stg_north_l, id_stg_north, 1, MPI_INTEGER,       &
                           MPI_SUM, comm1dy, ierr )

    ENDIF 
 
!
!-- Broadcast the id of the inflow PE
    IF ( bc_dirichlet_l )  THEN
       id_inflow_l = myidx
    ELSE
       id_inflow_l = 0
    ENDIF
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( id_inflow_l, id_inflow, 1, MPI_INTEGER, MPI_SUM, &
                        comm1dx, ierr )

!
!-- Broadcast the id of the recycling plane
!-- WARNING: needs to be adjusted in case of inflows other than from left side!
    IF ( turbulent_inflow ) THEN
    
       IF ( NINT( recycling_width / dx, KIND=idp ) >= nxl  .AND.                                   &
            NINT( recycling_width / dx, KIND=idp ) <= nxr )  THEN
          id_recycling_l = myidx
       ELSE
          id_recycling_l = 0
       ENDIF
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( id_recycling_l, id_recycling, 1, MPI_INTEGER, MPI_SUM, &
                           comm1dx, ierr )
                           
    ENDIF

!
!-- Broadcast the id of the outflow PE and outflow-source plane
    IF ( turbulent_outflow )  THEN

       IF ( bc_radiation_r )  THEN
          id_outflow_l = myidx
       ELSE
          id_outflow_l = 0
       ENDIF
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( id_outflow_l, id_outflow, 1, MPI_INTEGER, MPI_SUM, &
                           comm1dx, ierr )

       IF ( NINT( outflow_source_plane / dx ) >= nxl  .AND. &
            NINT( outflow_source_plane / dx ) <= nxr )  THEN
          id_outflow_source_l = myidx
       ELSE
          id_outflow_source_l = 0
       ENDIF
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( id_outflow_source_l, id_outflow_source, 1, &
                           MPI_INTEGER, MPI_SUM, comm1dx, ierr )

    ENDIF

    CALL location_message( 'creating virtual PE grids + MPI derived data types', 'finished' )

#else
    IF ( bc_lr == 'dirichlet/radiation' )  THEN
       bc_dirichlet_l = .TRUE.
       bc_radiation_r = .TRUE.
    ELSEIF ( bc_lr == 'radiation/dirichlet' )  THEN
       bc_radiation_l = .TRUE.
       bc_dirichlet_r = .TRUE.
    ENDIF

    IF ( bc_ns == 'dirichlet/radiation' )  THEN
       bc_dirichlet_n = .TRUE.
       bc_radiation_s = .TRUE.
    ELSEIF ( bc_ns == 'radiation/dirichlet' )  THEN
       bc_radiation_n = .TRUE.
       bc_dirichlet_s = .TRUE.
    ENDIF
#endif

!
!-- At the inflow or outflow, u or v, respectively, have to be calculated for
!-- one more grid point.
    IF ( bc_dirichlet_l  .OR.  bc_radiation_l )  THEN
       nxlu = nxl + 1
    ELSE
       nxlu = nxl
    ENDIF
    IF ( bc_dirichlet_s  .OR.  bc_radiation_s )  THEN
       nysv = nys + 1
    ELSE
       nysv = nys
    ENDIF

!
!-- Allocate wall flag arrays used in the multigrid solver
    IF ( psolver(1:9) == 'multigrid' )  THEN

       DO  i = maximum_grid_level, 1, -1

           SELECT CASE ( i )

              CASE ( 1 )
                 ALLOCATE( wall_flags_1(nzb:nzt_mg(i)+1,         &
                                        nys_mg(i)-1:nyn_mg(i)+1, &
                                        nxl_mg(i)-1:nxr_mg(i)+1) )

              CASE ( 2 )
                 ALLOCATE( wall_flags_2(nzb:nzt_mg(i)+1,         &
                                        nys_mg(i)-1:nyn_mg(i)+1, &
                                        nxl_mg(i)-1:nxr_mg(i)+1) )

              CASE ( 3 )
                 ALLOCATE( wall_flags_3(nzb:nzt_mg(i)+1,         &
                                        nys_mg(i)-1:nyn_mg(i)+1, &
                                        nxl_mg(i)-1:nxr_mg(i)+1) )

              CASE ( 4 )
                 ALLOCATE( wall_flags_4(nzb:nzt_mg(i)+1,         &
                                        nys_mg(i)-1:nyn_mg(i)+1, &
                                        nxl_mg(i)-1:nxr_mg(i)+1) )

              CASE ( 5 )
                 ALLOCATE( wall_flags_5(nzb:nzt_mg(i)+1,         &
                                        nys_mg(i)-1:nyn_mg(i)+1, &
                                        nxl_mg(i)-1:nxr_mg(i)+1) )

              CASE ( 6 )
                 ALLOCATE( wall_flags_6(nzb:nzt_mg(i)+1,         &
                                        nys_mg(i)-1:nyn_mg(i)+1, &
                                        nxl_mg(i)-1:nxr_mg(i)+1) )

              CASE ( 7 )
                 ALLOCATE( wall_flags_7(nzb:nzt_mg(i)+1,         &
                                        nys_mg(i)-1:nyn_mg(i)+1, &
                                        nxl_mg(i)-1:nxr_mg(i)+1) )

              CASE ( 8 )
                 ALLOCATE( wall_flags_8(nzb:nzt_mg(i)+1,         &
                                        nys_mg(i)-1:nyn_mg(i)+1, &
                                        nxl_mg(i)-1:nxr_mg(i)+1) )

              CASE ( 9 )
                 ALLOCATE( wall_flags_9(nzb:nzt_mg(i)+1,         &
                                        nys_mg(i)-1:nyn_mg(i)+1, &
                                        nxl_mg(i)-1:nxr_mg(i)+1) )

              CASE ( 10 )
                 ALLOCATE( wall_flags_10(nzb:nzt_mg(i)+1,        &
                                        nys_mg(i)-1:nyn_mg(i)+1, &
                                        nxl_mg(i)-1:nxr_mg(i)+1) )

              CASE DEFAULT
                 message_string = 'more than 10 multigrid levels'
                 CALL message( 'init_pegrid', 'PA0238', 1, 2, 0, 6, 0 )

          END SELECT

       ENDDO

    ENDIF

 END SUBROUTINE init_pegrid
