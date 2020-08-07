!> @file vertical_nesting_mod.f90
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
! Copyright 2017-2018 Karlsruhe Institute of Technology
!------------------------------------------------------------------------------!
!
! Current revisions:
! -----------------
! 
! 
! Former revisions:
! -----------------
! $Id: vertical_nesting_mod.f90 4457 2020-03-11 14:20:43Z raasch $
! use statement for exchange horiz added
! 
! 4444 2020-03-05 15:59:50Z raasch
! bugfix: cpp-directives for serial mode added
! 
! 4360 2020-01-07 11:25:50Z suehring
! Corrected "Former revisions" section
! 
! 4102 2019-07-17 16:00:03Z suehring
! - Slightly revise setting of boundary conditions at horizontal walls, use 
!   data-structure offset index instead of pre-calculate it for each facing
! 
! 4101 2019-07-17 15:14:26Z gronemeier
! remove old_dt
! 
! 3802 2019-03-17 13:33:42Z raasch
! unused subroutines commented out
! 
! 3655 2019-01-07 16:51:22Z knoop
! unused variables removed
! 
! 2365 2017-08-21 14:59:59Z kanani
! Initial revision (SadiqHuq)
! 
! 
! Description:
! ------------
!> Module for vertical nesting.
!>
!> Definition of parameters and variables for vertical nesting
!> The horizontal extent of the parent (Coarse Grid) and the child (Fine Grid)
!> have to be identical. The vertical extent of the FG should be smaller than CG.
!> Only integer nesting ratio supported. Odd nesting ratio preferred
!> The code follows MPI-1 standards. The available processors are split into
!> two groups using MPI_COMM_SPLIT. Exchange of data from CG to FG is called
!> interpolation. FG initialization by interpolation is done once at the start.
!> FG boundary conditions are set by interpolated at every timestep. 
!> Exchange of data from CG to FG is called anterpolation, the two-way interaction
!> occurs at every timestep.
!> vnest_start_time set in PARIN allows delayed start of the coupling
!> after spin-up of the CG
!>
!> @todo Replace dz(1) appropriatly to account for grid stretching
!> @todo Ensure that code can be compiled for serial and parallel mode. Please
!>       check the placement of the directive "__parallel".
!> @todo Add descriptions for all declared variables/parameters, one declaration
!>       statement per variable
!> @todo Add a descriptive header above each subroutine (see land_surface_model)
!> @todo FORTRAN language statements must not be used as names for variables 
!>       (e.g. if). Please rename it.
!> @todo Revise code according to PALM Coding Standard
!------------------------------------------------------------------------------!
 MODULE vertical_nesting_mod

    USE exchange_horiz_mod,                                                    &
        ONLY:  exchange_horiz, exchange_horiz_2d

    USE kinds

    IMPLICIT NONE

    LOGICAL                                   ::  vnested = .FALSE.            !> set to true if palmrun
                                                                               !> provides specific information via stdin
    LOGICAL                                   ::  vnest_init = .FALSE.         !> set to true when FG is initialized
    REAL(wp)                                  ::  vnest_start_time = 9999999.9 !> simulated time when FG should be 
                                                                               !> initialized. Should be
                                                                               !> identical in PARIN & PARIN_N

#if defined( __parallel )

    INTEGER(iwp),DIMENSION(3,2)               :: bdims = 0        !> sub-domain grid topology of current PE
    INTEGER(iwp),DIMENSION(3,2)               :: bdims_rem = 0    !> sub-domain grid topology of partner PE
    INTEGER(iwp)                              :: cg_nprocs        !> no. of PE in CG. Set by palmrun -Y
    INTEGER(iwp)                              :: fg_nprocs        !> no. of PE in FG. Set by palmrun -Y
    INTEGER(iwp)                              :: TYPE_VNEST_BC    !> derived contiguous data type for interpolation
    INTEGER(iwp)                              :: TYPE_VNEST_ANTER !> derived contiguous data type for anterpolation
    INTEGER(iwp),DIMENSION(:,:,:),ALLOCATABLE :: c2f_dims_cg      !> One CG PE sends data to multiple FG PEs
                                                                  !> list of grid-topology of partners
    INTEGER(iwp),DIMENSION(:,:,:),ALLOCATABLE :: f2c_dims_cg      !> One CG PE receives data from multiple FG PEs
                                                                  !> list of grid-topology of partners
    INTEGER(iwp),DIMENSION(:),ALLOCATABLE     :: c2f_dims_fg      !> One FG PE sends data to multiple CG PE
                                                                  !> list of grid-topology of partner
    INTEGER(iwp),DIMENSION(:),ALLOCATABLE     :: f2c_dims_fg      !> One FG PE sends data to only one CG PE
                                                                  !> list of grid-topology of partner

    INTEGER(iwp),DIMENSION(:,:),ALLOCATABLE   :: f_rnk_lst        !> list storing rank of FG PE denoted by pdims
    INTEGER(iwp),DIMENSION(:,:),ALLOCATABLE   :: c_rnk_lst        !> list storing rank of CG PE denoted by pdims
    INTEGER(iwp),DIMENSION(3)                 :: cfratio          !> Nesting ratio in x,y and z-directions

    INTEGER(iwp)                              :: nxc              !> no. of CG grid points in x-direction
    INTEGER(iwp)                              :: nxf              !> no. of FG grid points in x-direction
    INTEGER(iwp)                              :: nyc              !> no. of CG grid points in y-direction
    INTEGER(iwp)                              :: nyf              !> no. of FG grid points in y-direction
    INTEGER(iwp)                              :: nzc              !> no. of CG grid points in z-direction
    INTEGER(iwp)                              :: nzf              !> no. of FG grid points in z-direction
    INTEGER(iwp)                              :: ngp_c            !> no. of CG grid points in one vertical level
    INTEGER(iwp)                              :: ngp_f            !> no. of FG grid points in one vertical level

    INTEGER(iwp)                              :: n_cell_c         !> total no. of CG grid points in a PE
    INTEGER(iwp),DIMENSION(2)                 :: pdims_partner    !> processor topology of partner PE
    INTEGER(iwp)                              :: target_idex      !> temporary variable
    INTEGER(iwp),DIMENSION(2)                 :: offset           !> temporary variable
    INTEGER(iwp),DIMENSION(2)                 :: map_coord        !> temporary variable

    REAL(wp)                                  :: dxc              !> CG grid pacing in x-direction
    REAL(wp)                                  :: dyc              !> FG grid pacing in x-direction
    REAL(wp)                                  :: dxf              !> CG grid pacing in y-direction
    REAL(wp)                                  :: dyf              !> FG grid pacing in y-direction
    REAL(wp)                                  :: dzc              !> CG grid pacing in z-direction
    REAL(wp)                                  :: dzf              !> FG grid pacing in z-direction
    REAL(wp)                                  :: dtc              !> dt calculated for CG
    REAL(wp)                                  :: dtf              !> dt calculated for FG

    REAL(wp), DIMENSION(:),    ALLOCATABLE    :: zuc              !> CG vertical u-levels
    REAL(wp), DIMENSION(:),    ALLOCATABLE    :: zuf              !> FG vertical u-levels
    REAL(wp), DIMENSION(:),    ALLOCATABLE    :: zwc              !> CG vertical w-levels
    REAL(wp), DIMENSION(:),    ALLOCATABLE    :: zwf              !> FG vertical w-levels
    REAL(wp), DIMENSION(:,:,:), POINTER       :: interpol3d       !> pointers to simplify function calls
    REAL(wp), DIMENSION(:,:,:), POINTER       :: anterpol3d       !> pointers to simplify function calls


    REAL(wp),DIMENSION(:,:,:), ALLOCATABLE    :: work3d           !> temporary array for exchange of 3D data
    REAL(wp),DIMENSION(:,:),   ALLOCATABLE    :: work2dusws       !> temporary array for exchange of 2D data
    REAL(wp),DIMENSION(:,:),   ALLOCATABLE    :: work2dvsws       !> temporary array for exchange of 2D data
    REAL(wp),DIMENSION(:,:),   ALLOCATABLE    :: work2dts         !> temporary array for exchange of 2D data
    REAL(wp),DIMENSION(:,:),   ALLOCATABLE    :: work2dus         !> temporary array for exchange of 2D data

    SAVE

!-- Public functions
    PUBLIC vnest_init_fine, vnest_boundary_conds, vnest_anterpolate,          &
           vnest_boundary_conds_khkm, vnest_anterpolate_e,                    &
           vnest_init_pegrid_rank, vnest_init_pegrid_domain, vnest_init_grid, &
           vnest_timestep_sync, vnest_deallocate

!-- Public constants and variables
    PUBLIC vnested, vnest_init, vnest_start_time

    PRIVATE bdims, bdims_rem,                                                 &
            work3d, work2dusws, work2dvsws, work2dts, work2dus,               &
            dxc, dyc, dxf, dyf, dzc, dzf, dtc, dtf,                           &
            zuc, zuf, zwc, zwf, interpol3d, anterpol3d,                       &
            cg_nprocs, fg_nprocs,                                             &
            c2f_dims_cg, c2f_dims_fg, f2c_dims_cg, f2c_dims_fg,               &
            f_rnk_lst, c_rnk_lst, cfratio, pdims_partner,                     &
            nxc, nxf, nyc, nyf, nzc, nzf,                                     &
            ngp_c, ngp_f, target_idex, n_cell_c,                              &
            offset, map_coord, TYPE_VNEST_BC, TYPE_VNEST_ANTER                

    INTERFACE vnest_anterpolate
       MODULE PROCEDURE vnest_anterpolate
    END INTERFACE vnest_anterpolate

    INTERFACE vnest_anterpolate_e
       MODULE PROCEDURE vnest_anterpolate_e
    END INTERFACE vnest_anterpolate_e

    INTERFACE vnest_boundary_conds
       MODULE PROCEDURE vnest_boundary_conds
    END INTERFACE vnest_boundary_conds

    INTERFACE vnest_boundary_conds_khkm
       MODULE PROCEDURE vnest_boundary_conds_khkm
    END INTERFACE vnest_boundary_conds_khkm

    INTERFACE vnest_check_parameters
       MODULE PROCEDURE vnest_check_parameters
    END INTERFACE vnest_check_parameters

    INTERFACE vnest_deallocate
       MODULE PROCEDURE vnest_deallocate
    END INTERFACE vnest_deallocate

    INTERFACE vnest_init_fine
       MODULE PROCEDURE vnest_init_fine
    END INTERFACE vnest_init_fine

    INTERFACE vnest_init_grid
       MODULE PROCEDURE vnest_init_grid
    END INTERFACE vnest_init_grid

    INTERFACE vnest_init_pegrid_domain
       MODULE PROCEDURE vnest_init_pegrid_domain
    END INTERFACE vnest_init_pegrid_domain

    INTERFACE vnest_init_pegrid_rank
       MODULE PROCEDURE vnest_init_pegrid_rank
    END INTERFACE vnest_init_pegrid_rank

    INTERFACE vnest_timestep_sync
       MODULE PROCEDURE vnest_timestep_sync
    END INTERFACE vnest_timestep_sync

 CONTAINS
    
    
    
    SUBROUTINE vnest_init_fine
#if defined( __parallel )
   
        !--------------------------------------------------------------------------------!
        ! Description:
        ! ------------
        ! At the specified vnest_start_time initialize the Fine Grid based on the coarse
        ! grid values
        !------------------------------------------------------------------------------!
   
   
        USE arrays_3d
        USE control_parameters
        USE grid_variables
        USE indices
        USE interfaces
        USE pegrid
        USE turbulence_closure_mod,                                            &
            ONLY :  tcm_diffusivities
        
   
        IMPLICIT NONE
   
        REAL(wp)                              :: time_since_reference_point_rem
        INTEGER(iwp)                          :: i
        INTEGER(iwp)                          :: j
        INTEGER(iwp)                          :: iif
        INTEGER(iwp)                          :: jjf
        INTEGER(iwp)                          :: kkf
   
   
        if (myid ==0 ) print *, ' TIME TO INIT FINE from COARSE', simulated_time
   
        !
        !-- In case of model termination initiated by the remote model
        !-- (terminate_coupled_remote > 0), initiate termination of the local model.
        !-- The rest of the coupler must then be skipped because it would cause an MPI
        !-- intercomminucation hang.
        !-- If necessary, the coupler will be called at the beginning of the next
        !-- restart run.
   
        IF ( myid == 0) THEN
            CALL MPI_SENDRECV( terminate_coupled,        1, MPI_INTEGER,       &
                target_id, 0,                                                  &
                terminate_coupled_remote, 1, MPI_INTEGER,                      &
                target_id, 0,                                                  &
                comm_inter, status, ierr )
        ENDIF
        CALL MPI_BCAST( terminate_coupled_remote, 1, MPI_INTEGER, 0, comm2d, &
            ierr )
   
        IF ( terminate_coupled_remote > 0 )  THEN
            WRITE( message_string, * ) 'remote model "',                       &
                TRIM( coupling_mode_remote ),                                  &
                '" terminated',                                                &
                '&with terminate_coupled_remote = ',                           &
                terminate_coupled_remote,                                      &
                '&local model  "', TRIM( coupling_mode ),                      &
                '" has',                                                       &
                '&terminate_coupled = ',                                       &
                terminate_coupled
            CALL message( 'vnest_init_fine', 'PA0310', 1, 2, 0, 6, 0 )
            RETURN
        ENDIF
   
   
        !
        !-- Exchange the current simulated time between the models,
        !-- currently just for total_2ding
        IF ( myid == 0 ) THEN
   
            CALL MPI_SEND( time_since_reference_point, 1, MPI_REAL, target_id, &
                11, comm_inter, ierr )
            CALL MPI_RECV( time_since_reference_point_rem, 1, MPI_REAL,        &
                target_id, 11, comm_inter, status, ierr )
   
        ENDIF
   
        CALL MPI_BCAST( time_since_reference_point_rem, 1, MPI_REAL, 0, comm2d, &
            ierr )
   
   
        IF ( coupling_mode == 'vnested_crse' )  THEN
!-- Send data to fine grid for initialization
   
            offset(1) = ( pdims_partner(1) / pdims(1) ) * pcoord(1)
            offset(2) = ( pdims_partner(2) / pdims(2) ) * pcoord(2)
   
            do j = 0,   ( pdims_partner(2) / pdims(2) ) - 1
                do i = 0,   ( pdims_partner(1) / pdims(1) ) - 1
                    map_coord(1) = i+offset(1)
                    map_coord(2) = j+offset(2)
   
                    target_idex = f_rnk_lst(map_coord(1),map_coord(2)) + numprocs
   
                    CALL MPI_RECV( bdims_rem,   6, MPI_INTEGER, target_idex, 10, &
                        comm_inter,status, ierr )
   
                    bdims (1,1) = bdims_rem (1,1) / cfratio(1)
                    bdims (1,2) = bdims_rem (1,2) / cfratio(1)
                    bdims (2,1) = bdims_rem (2,1) / cfratio(2)
                    bdims (2,2) = bdims_rem (2,2) / cfratio(2)
                    bdims (3,1) = bdims_rem (3,1)
                    bdims (3,2) = bdims_rem (3,2) / cfratio(3)
   
   
                    CALL MPI_SEND( bdims, 6, MPI_INTEGER, target_idex, 9, &
                        comm_inter, ierr )
   
   
                    n_cell_c = (bdims(1,2)-bdims(1,1)+3) * &
                        (bdims(2,2)-bdims(2,1)+3) * &
                        (bdims(3,2)-bdims(3,1)+3)
   
                    CALL MPI_SEND( u( bdims(3,1):bdims(3,2)+2, &
                        bdims(2,1)-1:bdims(2,2)+1, &
                        bdims(1,1)-1:bdims(1,2)+1),&
                        n_cell_c, MPI_REAL, target_idex,  &
                        101, comm_inter, ierr)
   
                    CALL MPI_SEND( v( bdims(3,1):bdims(3,2)+2, &
                        bdims(2,1)-1:bdims(2,2)+1, &
                        bdims(1,1)-1:bdims(1,2)+1),&
                        n_cell_c, MPI_REAL, target_idex,  &
                        102, comm_inter, ierr)
   
                    CALL MPI_SEND( w( bdims(3,1):bdims(3,2)+2, &
                        bdims(2,1)-1:bdims(2,2)+1, &
                        bdims(1,1)-1:bdims(1,2)+1),&
                        n_cell_c, MPI_REAL, target_idex,  &
                        103, comm_inter, ierr)
   
                    CALL MPI_SEND( pt(bdims(3,1):bdims(3,2)+2, &
                        bdims(2,1)-1:bdims(2,2)+1, &
                        bdims(1,1)-1:bdims(1,2)+1),&
                        n_cell_c, MPI_REAL, target_idex,  &
                        105, comm_inter, ierr)
   
            IF ( humidity )  THEN
                    CALL MPI_SEND( q(bdims(3,1):bdims(3,2)+2, &
                        bdims(2,1)-1:bdims(2,2)+1, &
                        bdims(1,1)-1:bdims(1,2)+1),&
                        n_cell_c, MPI_REAL, target_idex,  &
                        116, comm_inter, ierr)
            ENDIF
 
                     CALL MPI_SEND( e( bdims(3,1):bdims(3,2)+2, &
                        bdims(2,1)-1:bdims(2,2)+1, &
                        bdims(1,1)-1:bdims(1,2)+1),&
                        n_cell_c, MPI_REAL, target_idex,  &
                        104, comm_inter, ierr)
   
                   CALL MPI_SEND(kh( bdims(3,1):bdims(3,2)+2, &
                        bdims(2,1)-1:bdims(2,2)+1, &
                        bdims(1,1)-1:bdims(1,2)+1),&
                        n_cell_c, MPI_REAL, target_idex,  &
                        106, comm_inter, ierr)
   
                    CALL MPI_SEND(km( bdims(3,1):bdims(3,2)+2, &
                        bdims(2,1)-1:bdims(2,2)+1, &
                        bdims(1,1)-1:bdims(1,2)+1),&
                        n_cell_c, MPI_REAL, target_idex,  &
                        107, comm_inter, ierr)
   
!-- Send Surface fluxes
            IF ( use_surface_fluxes )  THEN
   
                   n_cell_c = (bdims(1,2)-bdims(1,1)+3) * &
                        (bdims(2,2)-bdims(2,1)+3)
   
!
!--     shf and z0 for CG / FG need to initialized in input file or user_code
!--     TODO
!--     initialization of usws, vsws, ts and us not vital to vnest FG
!--     variables are not compatible with the new surface layer module
!
!                 CALL MPI_SEND(surf_def_h(0)%usws( bdims(2,1)-1:bdims(2,2)+1, &
!                     bdims(1,1)-1:bdims(1,2)+1),&
!                     n_cell_c, MPI_REAL, target_idex,  &
!                     110, comm_inter, ierr   )
!
!                 CALL MPI_SEND(surf_def_h(0)%vsws( bdims(2,1)-1:bdims(2,2)+1, &
!                     bdims(1,1)-1:bdims(1,2)+1),&
!                     n_cell_c, MPI_REAL, target_idex,  &
!                     111, comm_inter, ierr   )
!
!                 CALL MPI_SEND(ts  ( bdims(2,1)-1:bdims(2,2)+1, &
!                     bdims(1,1)-1:bdims(1,2)+1),&
!                     n_cell_c, MPI_REAL, target_idex,  &
!                     112, comm_inter, ierr   )
!   
!                 CALL MPI_SEND(us  ( bdims(2,1)-1:bdims(2,2)+1, &
!                     bdims(1,1)-1:bdims(1,2)+1),&
!                     n_cell_c, MPI_REAL, target_idex,  &
!                     113, comm_inter, ierr   )
!   
            ENDIF
   
   
   
   
                end do
            end do
   
        ELSEIF ( coupling_mode == 'vnested_fine' )  THEN
!-- Receive data from coarse grid for initialization
   
            offset(1) = pcoord(1) / ( pdims(1)/pdims_partner(1) )
            offset(2) = pcoord(2) / ( pdims(2)/pdims_partner(2) )
            map_coord(1) = offset(1)
            map_coord(2) = offset(2)
            target_idex = c_rnk_lst(map_coord(1),map_coord(2))
   
            bdims (1,1) = nxl
            bdims (1,2) = nxr
            bdims (2,1) = nys
            bdims (2,2) = nyn
            bdims (3,1) = nzb
            bdims (3,2) = nzt
   
            CALL MPI_SEND( bdims,       6, MPI_INTEGER, target_idex, 10, &
                comm_inter, ierr )
   
            CALL MPI_RECV( bdims_rem,   6, MPI_INTEGER, target_idex, 9, &
                comm_inter,status, ierr )
   
            n_cell_c = (bdims_rem(1,2)-bdims_rem(1,1)+3) * &
                (bdims_rem(2,2)-bdims_rem(2,1)+3) * &
                (bdims_rem(3,2)-bdims_rem(3,1)+3)
   
            ALLOCATE( work3d ( bdims_rem(3,1)  :bdims_rem(3,2)+2, &
                bdims_rem(2,1)-1:bdims_rem(2,2)+1, &
                bdims_rem(1,1)-1:bdims_rem(1,2)+1))
   
   
            CALL MPI_RECV( work3d,n_cell_c, MPI_REAL, target_idex, 101, &
                comm_inter,status, ierr )
            interpol3d => u
            call interpolate_to_fine_u
   
            CALL MPI_RECV( work3d,n_cell_c, MPI_REAL, target_idex, 102, &
                comm_inter,status, ierr )
            interpol3d => v
            call interpolate_to_fine_v
   
            CALL MPI_RECV( work3d,n_cell_c, MPI_REAL, target_idex, 103, &
                comm_inter,status, ierr )
            interpol3d => w
            call interpolate_to_fine_w
   
            CALL MPI_RECV( work3d,n_cell_c, MPI_REAL, target_idex, 105, &
                comm_inter,status, ierr )
            interpol3d => pt
            call interpolate_to_fine_s
   
            IF ( humidity )  THEN
            CALL MPI_RECV( work3d,n_cell_c, MPI_REAL, target_idex, 116, &
                comm_inter,status, ierr )
            interpol3d => q
            call interpolate_to_fine_s
            ENDIF

            CALL MPI_RECV( work3d,n_cell_c, MPI_REAL, target_idex, 104, &
                comm_inter,status, ierr )
            interpol3d => e
            call interpolate_to_fine_s
   
!-- kh,km no target attribute, use of pointer not possible
            CALL MPI_RECV( work3d,n_cell_c, MPI_REAL, target_idex, 106, &
                comm_inter,status, ierr )
            call interpolate_to_fine_kh
   
            CALL MPI_RECV( work3d,n_cell_c, MPI_REAL, target_idex, 107, &
                comm_inter,status, ierr )
            call interpolate_to_fine_km
   
            DEALLOCATE(   work3d       )
            NULLIFY   (   interpol3d   )
   
!-- Recv Surface Fluxes    
            IF ( use_surface_fluxes )  THEN
            n_cell_c = (bdims_rem(1,2)-bdims_rem(1,1)+3) * &
                (bdims_rem(2,2)-bdims_rem(2,1)+3)
   
            ALLOCATE( work2dusws ( bdims_rem(2,1)-1:bdims_rem(2,2)+1, &
                bdims_rem(1,1)-1:bdims_rem(1,2)+1) )
            ALLOCATE( work2dvsws ( bdims_rem(2,1)-1:bdims_rem(2,2)+1, &
                bdims_rem(1,1)-1:bdims_rem(1,2)+1) )
            ALLOCATE( work2dts   ( bdims_rem(2,1)-1:bdims_rem(2,2)+1, &
                bdims_rem(1,1)-1:bdims_rem(1,2)+1) )
            ALLOCATE( work2dus   ( bdims_rem(2,1)-1:bdims_rem(2,2)+1, &
                bdims_rem(1,1)-1:bdims_rem(1,2)+1) )
   
!
!--     shf and z0 for CG / FG need to initialized in input file or user_code
!--     TODO
!--     initialization of usws, vsws, ts and us not vital to vnest FG
!--     variables are not compatible with the new surface layer module
! 
!          CALL MPI_RECV( work2dusws,n_cell_c, MPI_REAL, target_idex, 110, &
!              comm_inter,status, ierr )
! 
!          CALL MPI_RECV( work2dvsws,n_cell_c, MPI_REAL, target_idex, 111, &
!              comm_inter,status, ierr )
! 
!          CALL MPI_RECV( work2dts  ,n_cell_c, MPI_REAL, target_idex, 112, &
!              comm_inter,status, ierr )
!   
!          CALL MPI_RECV( work2dus  ,n_cell_c, MPI_REAL, target_idex, 113, &
!              comm_inter,status, ierr )
!   
!           CALL interpolate_to_fine_flux ( 108 )
   
            DEALLOCATE(   work2dusws   )
            DEALLOCATE(   work2dvsws   )
            DEALLOCATE(   work2dts     )
            DEALLOCATE(   work2dus     )
          ENDIF
   
          IF ( .NOT. constant_diffusion ) THEN
             DO kkf =  bdims(3,1)+1,bdims(3,2)+1
                 DO jjf =  bdims(2,1),bdims(2,2)
                     DO iif =  bdims(1,1),bdims(1,2)
   
                         IF ( e(kkf,jjf,iif) < 0.0 ) THEN
                              e(kkf,jjf,iif) = 1E-15_wp
                         END IF
   
                     END DO
                 END DO
             END DO
          ENDIF
   
           w(nzt+1,:,:) = w(nzt,:,:)
          
           CALL exchange_horiz( u, nbgp )
           CALL exchange_horiz( v, nbgp )
           CALL exchange_horiz( w, nbgp )
           CALL exchange_horiz( pt, nbgp )
           IF ( .NOT. constant_diffusion )  CALL exchange_horiz( e, nbgp )
           IF ( humidity )  CALL exchange_horiz( q, nbgp )
   
           !
           !--         Velocity boundary conditions at the bottom boundary
           IF ( ibc_uv_b == 0 ) THEN
               u(nzb,:,:) = 0.0_wp
               v(nzb,:,:) = 0.0_wp
           ELSE
               u(nzb,:,:) = u(nzb+1,:,:)
               v(nzb,:,:) = v(nzb+1,:,:)
           END IF
           
           
           w(nzb,:,:) = 0.0_wp

!
!-- Temperature boundary conditions at the bottom boundary
           IF ( ibc_pt_b /= 0 ) THEN
               pt(nzb,:,:) = pt(nzb+1,:,:)
           END IF
           
           !
           !-- Bottom boundary condition for the turbulent kinetic energy
           !-- Generally a Neumann condition with de/dz=0 is assumed
           IF ( .NOT. constant_diffusion ) THEN
               e(nzb,:,:) = e(nzb+1,:,:)
           END IF
           
           !
           !-- Bottom boundary condition for turbulent diffusion coefficients
           km(nzb,:,:) = km(nzb+1,:,:)
           kh(nzb,:,:) = kh(nzb+1,:,:)
           
           !diffusivities required
           IF ( .NOT. humidity ) THEN
               CALL tcm_diffusivities( pt, pt_reference )
           ELSE
               CALL tcm_diffusivities( vpt, pt_reference )
           ENDIF
   

!
!-- Reset Fine Grid top Boundary Condition
!-- At the top of the FG, the scalars always follow Dirichlet condition

            ibc_pt_t = 0

!-- Initialize old time levels
            pt_p = pt; u_p = u; v_p = v; w_p = w
            IF ( .NOT. constant_diffusion ) e_p = e
            IF ( humidity ) THEN
               ibc_q_t = 0
               q_p     = q
            ENDIF
   
        ENDIF
   

        if (myid==0) print *, '** Fine Initalized ** simulated_time:', simulated_time

    CONTAINS
   
       SUBROUTINE interpolate_to_fine_w
      
           USE arrays_3d
           USE control_parameters
           USE grid_variables
           USE indices
           USE pegrid
           
      
           IMPLICIT NONE
      
           INTEGER(iwp)                       :: i
           INTEGER(iwp)                       :: j
           INTEGER(iwp)                       :: k
           INTEGER(iwp)                       :: iif
           INTEGER(iwp)                       :: jjf
           INTEGER(iwp)                       :: kkf
           INTEGER(iwp)                       :: nzbottom
           INTEGER(iwp)                       :: nztop
           INTEGER(iwp)                       :: bottomx
           INTEGER(iwp)                       :: bottomy
           INTEGER(iwp)                       :: bottomz
           INTEGER(iwp)                       :: topx
           INTEGER(iwp)                       :: topy
           INTEGER(iwp)                       :: topz
           REAL(wp)                           :: eps
           REAL(wp)                           :: alpha
           REAL(wp)                           :: eminus
           REAL(wp)                           :: edot
           REAL(wp)                           :: eplus
           REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: wprs
           REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: wprf
      
      
           nzbottom = bdims_rem (3,1)
           nztop  = bdims_rem (3,2)
      
           ALLOCATE( wprf(nzbottom:nztop, bdims_rem(2,1)-1: bdims_rem(2,2)+1,nxl:nxr) )
           ALLOCATE( wprs(nzbottom:nztop,nys:nyn,nxl:nxr) )
      
      
           !
           !-- Initialisation of the velocity component w
           !
           !-- Interpolation in x-direction
           DO k = nzbottom, nztop
               DO j = bdims_rem(2,1)-1, bdims_rem(2,2)+1
                   DO i = bdims_rem(1,1),bdims_rem(1,2)
      
                       bottomx = (nxf+1)/(nxc+1) * i
                       topx    = (nxf+1)/(nxc+1) * (i+1) - 1
      
                       DO iif = bottomx, topx
      
                           eps    = ( iif * dxf + 0.5 * dxf - i * dxc - 0.5 * dxc ) / dxc
                           alpha  = ( ( dxf / dxc )**2.0 - 1.0 ) / 24.0
                           eminus = eps * ( eps - 1.0 ) / 2.0 + alpha
                           edot  = ( 1.0 - eps**2.0 ) - 2.0 * alpha
                           eplus  = eps * ( eps + 1.0 ) / 2.0 + alpha
                      
                           wprf(k,j,iif) = eminus * work3d(k,j,i-1) &
                               + edot  * work3d(k,j,i)   &
                               + eplus  * work3d(k,j,i+1)
                       END DO
      
                   END DO
               END DO
           END DO
      
           !
           !-- Interpolation in y-direction (quadratic, Clark and Farley)
           DO k = nzbottom, nztop
               DO j = bdims_rem(2,1), bdims_rem(2,2)
                
                   bottomy = (nyf+1)/(nyc+1) * j
                   topy    = (nyf+1)/(nyc+1) * (j+1) - 1
      
                   DO iif = nxl, nxr
                       DO jjf = bottomy, topy
                      
                           eps    = ( jjf * dyf + 0.5 * dyf - j * dyc - 0.5 * dyc ) / dyc
                           alpha  = ( ( dyf / dyc )**2.0 - 1.0 ) / 24.0
                           eminus = eps * ( eps - 1.0 ) / 2.0 + alpha
                           edot  = ( 1.0 - eps**2.0 ) - 2.0 * alpha
                           eplus  = eps * ( eps + 1.0 ) / 2.0 + alpha
      
                           wprs(k,jjf,iif) = eminus * wprf(k,j-1,iif) &
                               + edot  * wprf(k,j,iif)   &
                               + eplus  * wprf(k,j+1,iif)
      
                       END DO
                   END DO
      
               END DO
           END DO
      
           !
           !-- Interpolation in z-direction (linear)
      
           DO k = nzbottom, nztop-1
      
               bottomz = (dzc/dzf) * k
               topz    = (dzc/dzf) * (k+1) - 1
      
               DO jjf = nys, nyn
                   DO iif = nxl, nxr
                       DO kkf = bottomz, topz
      
                           w(kkf,jjf,iif) = wprs(k,jjf,iif) + ( zwf(kkf) - zwc(k) ) &
                               * ( wprs(k+1,jjf,iif) - wprs(k,jjf,iif) ) / dzc
      
                       END DO
                   END DO
               END DO
      
           END DO
      
           DO jjf = nys, nyn
               DO iif = nxl, nxr
      
                   w(nzt,jjf,iif) = wprs(nztop,jjf,iif)
       
               END DO
           END DO
           !
           !    w(nzb:nzt+1,nys:nyn,nxl:nxr) = 0
      
           DEALLOCATE( wprf, wprs )
      
       END SUBROUTINE interpolate_to_fine_w
      
       SUBROUTINE interpolate_to_fine_u
      
      
           USE arrays_3d
           USE control_parameters
           USE grid_variables
           USE indices
           USE pegrid
           
      
           IMPLICIT NONE
      
           INTEGER(iwp)                       :: i
           INTEGER(iwp)                       :: j
           INTEGER(iwp)                       :: k
           INTEGER(iwp)                       :: iif
           INTEGER(iwp)                       :: jjf
           INTEGER(iwp)                       :: kkf
           INTEGER(iwp)                       :: nzbottom
           INTEGER(iwp)                       :: nztop
           INTEGER(iwp)                       :: bottomx
           INTEGER(iwp)                       :: bottomy
           INTEGER(iwp)                       :: bottomz
           INTEGER(iwp)                       :: topx
           INTEGER(iwp)                       :: topy
           INTEGER(iwp)                       :: topz
           REAL(wp)                           :: eps
           REAL(wp)                           :: alpha
           REAL(wp)                           :: eminus
           REAL(wp)                           :: edot
           REAL(wp)                           :: eplus
           REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: uprf
           REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: uprs
      
      
      
           nzbottom = bdims_rem (3,1)
           nztop  = bdims_rem (3,2)
      
           ALLOCATE( uprf(nzbottom:nztop+2,nys:nyn,bdims_rem(1,1)-1:bdims_rem(1,2)+1) )
           ALLOCATE( uprs(nzb+1:nzt+1,nys:nyn,bdims_rem(1,1)-1:bdims_rem(1,2)+1) )
      
           !
           !-- Initialisation of the velocity component uf
      
           !
           !-- Interpolation in y-direction (quadratic, Clark and Farley)
      
           DO k = nzbottom, nztop+2
               DO j = bdims_rem(2,1), bdims_rem(2,2)
      
                   bottomy = (nyf+1)/(nyc+1) * j
                   topy    = (nyf+1)/(nyc+1) * (j+1) - 1
      
                   DO i = bdims_rem(1,1)-1, bdims_rem(1,2)+1
                       DO jjf = bottomy, topy
      
                           eps    = ( jjf * dyf + 0.5 * dyf - j * dyc - 0.5 * dyc ) / dyc
                           alpha  = ( ( dyf / dyc )**2.0 - 1.0 ) / 24.0
                           eminus = eps * ( eps - 1.0 ) / 2.0 + alpha
                           edot  = ( 1.0 - eps**2.0 ) - 2.0 * alpha
                           eplus  = eps * ( eps + 1.0 ) / 2.0 + alpha
      
                           uprf(k,jjf,i) = eminus * work3d(k,j-1,i) &
                               + edot  * work3d(k,j,i)   &
                               + eplus  * work3d(k,j+1,i)
      
                       END DO
                   END DO
      
               END DO
           END DO
      
           !
           !-- Interpolation in z-direction (quadratic, Clark and Farley)
      
           DO k = nzbottom+1, nztop
      
               bottomz = (dzc/dzf) * (k-1) + 1
               topz    = (dzc/dzf) * k
      
               DO jjf = nys, nyn
                   DO i = bdims_rem(1,1)-1, bdims_rem(1,2)+1
                       DO kkf = bottomz, topz
                      
                           eps    = ( zuf(kkf) - zuc(k) ) / dzc
                           alpha  = ( ( dzf / dzc )**2.0 - 1.0 ) / 24.0
                           eminus = eps * ( eps - 1.0 ) / 2.0 + alpha
                           edot  = ( 1.0 - eps**2.0 ) - 2.0 * alpha
                           eplus  = eps * ( eps + 1.0 ) / 2.0 + alpha
      
                           uprs(kkf,jjf,i) = eminus * uprf(k-1,jjf,i) &
                               + edot  * uprf(k,jjf,i)   &
                               + eplus  * uprf(k+1,jjf,i)
      
                       END DO
                   END DO
               END DO
      
           END DO
      
           DO jjf = nys, nyn
               DO i = bdims_rem(1,1)-1, bdims_rem(1,2)+1
      
                   eps    = ( zuf(nzt+1) - zuc(nztop+1) ) / dzc
                   alpha  = ( ( dzf / dzc )**2.0 - 1.0 ) / 24.0
                   eminus = eps * ( eps - 1.0 ) / 2.0 + alpha
                   edot  = ( 1.0 - eps**2.0 ) - 2.0 * alpha
                   eplus  = eps * ( eps + 1.0 ) / 2.0 + alpha
      
                   uprs(nzt+1,jjf,i)  = eminus * uprf(nztop,jjf,i)   &
                       + edot  * uprf(nztop+1,jjf,i) &
                       + eplus  * uprf(nztop+2,jjf,i)
      
               END DO
           END DO
      
           !
           !-- Interpolation in x-direction (linear)
      
           DO kkf = nzb+1, nzt+1
               DO jjf = nys, nyn
                   DO i = bdims_rem(1,1), bdims_rem(1,2)
      
                       bottomx = (nxf+1)/(nxc+1) * i
                       topx    = (nxf+1)/(nxc+1) * (i+1) - 1
      
                       DO iif = bottomx, topx
                           u(kkf,jjf,iif)  = uprs(kkf,jjf,i) + ( iif * dxf - i * dxc ) &
                               * ( uprs(kkf,jjf,i+1) - uprs(kkf,jjf,i) ) / dxc
                       END DO
      
                   END DO
               END DO
           END DO
           !
           !-- Determination of uf at the bottom boundary
      
      
      
           DEALLOCATE( uprf, uprs )
      
       END SUBROUTINE interpolate_to_fine_u
      
      
       SUBROUTINE interpolate_to_fine_v
      
      
           USE arrays_3d
           USE control_parameters
           USE grid_variables
           USE indices
           USE pegrid
           
      
           IMPLICIT NONE

           INTEGER(iwp)                       :: i
           INTEGER(iwp)                       :: j
           INTEGER(iwp)                       :: k
           INTEGER(iwp)                       :: iif
           INTEGER(iwp)                       :: jjf
           INTEGER(iwp)                       :: kkf
           INTEGER(iwp)                       :: nzbottom
           INTEGER(iwp)                       :: nztop
           INTEGER(iwp)                       :: bottomx
           INTEGER(iwp)                       :: bottomy
           INTEGER(iwp)                       :: bottomz
           INTEGER(iwp)                       :: topx
           INTEGER(iwp)                       :: topy
           INTEGER(iwp)                       :: topz
           REAL(wp)                           :: eps
           REAL(wp)                           :: alpha
           REAL(wp)                           :: eminus
           REAL(wp)                           :: edot
           REAL(wp)                           :: eplus
           REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: vprs  
           REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: vprf
          
      
           nzbottom = bdims_rem (3,1)
           nztop  = bdims_rem (3,2)
      
           ALLOCATE( vprf(nzbottom:nztop+2,bdims_rem(2,1)-1:bdims_rem(2,2)+1,nxl:nxr) )
           ALLOCATE( vprs(nzb+1:nzt+1, bdims_rem(2,1)-1:bdims_rem(2,2)+1,nxl:nxr) )
           !
           !-- Initialisation of the velocity component vf
      
           !
           !-- Interpolation in x-direction (quadratic, Clark and Farley)
      
           DO k = nzbottom, nztop+2
               DO j = bdims_rem(2,1)-1, bdims_rem(2,2)+1
                   DO i = bdims_rem(1,1), bdims_rem(1,2)
      
                       bottomx = (nxf+1)/(nxc+1) * i
                       topx    = (nxf+1)/(nxc+1) * (i+1) - 1
      
                       DO iif = bottomx, topx
      
                           eps    = ( iif * dxf + 0.5 * dxf - i * dxc - 0.5 * dxc ) / dxc
                           alpha  = ( ( dxf / dxc )**2.0 - 1.0 ) / 24.0
                           eminus = eps * ( eps - 1.0 ) / 2.0 + alpha
                           edot  = ( 1.0 - eps**2.0 ) - 2.0 * alpha
                           eplus  = eps * ( eps + 1.0 ) / 2.0 + alpha
      
                           vprf(k,j,iif) = eminus * work3d(k,j,i-1) &
                               + edot  * work3d(k,j,i)   &
                               + eplus  * work3d(k,j,i+1)
      
                       END DO
      
                   END DO
               END DO
           END DO
      
           !
           !-- Interpolation in z-direction (quadratic, Clark and Farley)
      
           DO k = nzbottom+1, nztop
             
               bottomz = (dzc/dzf) * (k-1) + 1
               topz    = (dzc/dzf) * k
      
               DO j = bdims_rem(2,1)-1, bdims_rem(2,2)+1
                   DO iif = nxl, nxr
                       DO kkf = bottomz, topz
      
                           eps    = ( zuf(kkf) - zuc(k) ) / dzc
                           alpha  = ( ( dzf / dzc )**2.0 - 1.0 ) / 24.0
                           eminus = eps * ( eps - 1.0 ) / 2.0 + alpha
                           edot  = ( 1.0 - eps**2.0 ) - 2.0 * alpha
                           eplus  = eps * ( eps + 1.0 ) / 2.0 + alpha
      
                           vprs(kkf,j,iif) = eminus * vprf(k-1,j,iif) &
                               + edot  * vprf(k,j,iif)   &
                               + eplus  * vprf(k+1,j,iif)
      
                       END DO
                   END DO
               END DO
      
           END DO
             
           DO j = bdims_rem(2,1)-1, bdims_rem(2,2)+1
               DO iif = nxl, nxr
      
                   eps    = ( zuf(nzt+1) - zuc(nztop+1) ) / dzc
                   alpha  = ( ( dzf / dzc )**2.0 - 1.0 ) / 24.0
                   eminus = eps * ( eps - 1.0 ) / 2.0 + alpha
                   edot  = ( 1.0 - eps**2.0 ) - 2.0 * alpha
                   eplus  = eps * ( eps + 1.0 ) / 2.0 + alpha
      
                   vprs(nzt+1,j,iif)  = eminus * vprf(nztop,j,iif)   &
                       + edot  * vprf(nztop+1,j,iif) &
                       + eplus  * vprf(nztop+2,j,iif)
      
               END DO
           END DO
      
           !
           !-- Interpolation in y-direction (linear)
      
           DO kkf = nzb+1, nzt+1
               DO j = bdims_rem(2,1), bdims_rem(2,2)
      
                   bottomy = (nyf+1)/(nyc+1) * j
                   topy    = (nyf+1)/(nyc+1) * (j+1) - 1
      
                   DO iif = nxl, nxr
                       DO jjf = bottomy, topy
                           v (kkf,jjf,iif) = vprs(kkf,j,iif) + ( jjf * dyf - j * dyc ) &
                               * ( vprs(kkf,j+1,iif) - vprs(kkf,j,iif) ) / dyc
                       END DO
                   END DO
        
               END DO
           END DO
      
           !
           !-- Determination of vf at the bottom boundary
      
      
           DEALLOCATE( vprf, vprs )
      
       END SUBROUTINE interpolate_to_fine_v
      
      
       SUBROUTINE interpolate_to_fine_s
      
      
           USE arrays_3d
           USE control_parameters
           USE grid_variables
           USE indices
           USE pegrid
           
      
           IMPLICIT NONE

           INTEGER(iwp)                       :: i
           INTEGER(iwp)                       :: j
           INTEGER(iwp)                       :: k
           INTEGER(iwp)                       :: iif
           INTEGER(iwp)                       :: jjf
           INTEGER(iwp)                       :: kkf
           INTEGER(iwp)                       :: nzbottom
           INTEGER(iwp)                       :: nztop
           INTEGER(iwp)                       :: bottomx
           INTEGER(iwp)                       :: bottomy
           INTEGER(iwp)                       :: bottomz
           INTEGER(iwp)                       :: topx
           INTEGER(iwp)                       :: topy
           INTEGER(iwp)                       :: topz
           REAL(wp)                           :: eps
           REAL(wp)                           :: alpha
           REAL(wp)                           :: eminus
           REAL(wp)                           :: edot
           REAL(wp)                           :: eplus
           REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: ptprs
           REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: ptprf
      
      
           nzbottom = bdims_rem (3,1)
           nztop  = bdims_rem (3,2)
      
           ALLOCATE( ptprf(nzbottom:nztop+2,bdims_rem(2,1)-1:bdims_rem(2,2)+1,nxl:nxr) )
           ALLOCATE( ptprs(nzbottom:nztop+2,nys:nyn,nxl:nxr) )
      
           !
           !-- Initialisation of scalar variables
      
           !
           !-- Interpolation in x-direction (quadratic, Clark and Farley)
      
           DO k = nzbottom, nztop+2
               DO j = bdims_rem(2,1)-1, bdims_rem(2,2)+1
                   DO i = bdims_rem(1,1), bdims_rem(1,2)
      
                       bottomx = (nxf+1)/(nxc+1) * i
                       topx    = (nxf+1)/(nxc+1) * (i+1) - 1
      
                       DO iif = bottomx, topx
      
                           eps    = ( iif * dxf + 0.5 * dxf - i * dxc - 0.5 * dxc ) / dxc
                           alpha  = ( ( dxf / dxc )**2.0 - 1.0 ) / 24.0
                           eminus = eps * ( eps - 1.0 ) / 2.0 + alpha
                           edot  = ( 1.0 - eps**2.0 ) - 2.0 * alpha
                           eplus  = eps * ( eps + 1.0 ) / 2.0 + alpha
      
                           ptprf(k,j,iif) = eminus * work3d(k,j,i-1) &
                               + edot  * work3d(k,j,i)   &
                               + eplus  * work3d(k,j,i+1)
                       END DO
      
                   END DO
               END DO
           END DO
      
           !
           !-- Interpolation in y-direction (quadratic, Clark and Farley)
      
           DO k = nzbottom, nztop+2
               DO j = bdims_rem(2,1), bdims_rem(2,2)
        
                   bottomy = (nyf+1)/(nyc+1) * j
                   topy    = (nyf+1)/(nyc+1) * (j+1) - 1
      
                   DO iif = nxl, nxr
                       DO jjf = bottomy, topy
      
                           eps    = ( jjf * dyf + 0.5 * dyf - j * dyc - 0.5 * dyc ) / dyc
                           alpha  = ( ( dyf / dyc )**2.0 - 1.0 ) / 24.0
                           eminus = eps * ( eps - 1.0 ) / 2.0 + alpha
                           edot  = ( 1.0 - eps**2.0 ) - 2.0 * alpha
                           eplus  = eps * ( eps + 1.0 ) / 2.0 + alpha
      
                           ptprs(k,jjf,iif) = eminus * ptprf(k,j-1,iif) &
                               + edot  * ptprf(k,j,iif)   &
                               + eplus  * ptprf(k,j+1,iif)
      
                       END DO
                   END DO
      
               END DO
           END DO
      
           !
           !-- Interpolation in z-direction (quadratic, Clark and Farley)
      
           DO k = nzbottom+1, nztop
            
               bottomz = (dzc/dzf) * (k-1) + 1
               topz    = (dzc/dzf) * k
      
               DO jjf = nys, nyn
                   DO iif = nxl, nxr
                       DO kkf = bottomz, topz
                       
                           eps    = ( zuf(kkf) - zuc(k) ) / dzc
                           alpha  = ( ( dzf / dzc )**2.0 - 1.0 ) / 24.0
                           eminus = eps * ( eps - 1.0 ) / 2.0 + alpha
                           edot  = ( 1.0 - eps**2.0 ) - 2.0 * alpha
                           eplus  = eps * ( eps + 1.0 ) / 2.0 + alpha
       
                           interpol3d(kkf,jjf,iif) = eminus * ptprs(k-1,jjf,iif) &
                               + edot  * ptprs(k,jjf,iif)   &
                               + eplus  * ptprs(k+1,jjf,iif)
      
                       END DO
                   END DO
               END DO
      
           END DO
               
           DO jjf = nys, nyn
               DO iif = nxl, nxr
      
                   eps    = ( zuf(nzt+1) - zuc(nztop+1) ) / dzc
                   alpha  = ( ( dzf / dzc )**2.0 - 1.0 ) / 24.0
                   eminus = eps * ( eps - 1.0 ) / 2.0 + alpha
                   edot  = ( 1.0 - eps**2.0 ) - 2.0 * alpha
                   eplus  = eps * ( eps + 1.0 ) / 2.0 + alpha
      
                   interpol3d(nzt+1,jjf,iif) = eminus * ptprs(nztop,jjf,iif)   &
                       + edot  * ptprs(nztop+1,jjf,iif) &
                       + eplus  * ptprs(nztop+2,jjf,iif)
      
               END DO
           END DO
      
      
           DEALLOCATE( ptprf, ptprs ) 
      
       END SUBROUTINE interpolate_to_fine_s
      
      
       SUBROUTINE interpolate_to_fine_kh
      
      
           USE arrays_3d
           USE control_parameters
           USE grid_variables
           USE indices
           USE pegrid
           
      
           IMPLICIT NONE

           INTEGER(iwp)                       :: i
           INTEGER(iwp)                       :: j
           INTEGER(iwp)                       :: k
           INTEGER(iwp)                       :: iif
           INTEGER(iwp)                       :: jjf
           INTEGER(iwp)                       :: kkf
           INTEGER(iwp)                       :: nzbottom
           INTEGER(iwp)                       :: nztop
           INTEGER(iwp)                       :: bottomx
           INTEGER(iwp)                       :: bottomy
           INTEGER(iwp)                       :: bottomz
           INTEGER(iwp)                       :: topx
           INTEGER(iwp)                       :: topy
           INTEGER(iwp)                       :: topz
           REAL(wp)                           :: eps
           REAL(wp)                           :: alpha
           REAL(wp)                           :: eminus
           REAL(wp)                           :: edot
           REAL(wp)                           :: eplus
           REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: ptprs
           REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: ptprf
      
      
           nzbottom = bdims_rem (3,1)
           nztop  = bdims_rem (3,2)
           !   nztop  = blk_dim_rem (3,2)+1
      
      
           ALLOCATE( ptprf(nzbottom:nztop+2,bdims_rem(2,1)-1:bdims_rem(2,2)+1,nxl:nxr) )
           ALLOCATE( ptprs(nzbottom:nztop+2,nys:nyn,nxl:nxr) )
      
      
           !
           !-- Initialisation of scalar variables
      
           !
           !-- Interpolation in x-direction (quadratic, Clark and Farley)
      
           DO k = nzbottom, nztop+2
               DO j = bdims_rem(2,1)-1, bdims_rem(2,2)+1
                   DO i = bdims_rem(1,1), bdims_rem(1,2)
      
                       bottomx = (nxf+1)/(nxc+1) * i
                       topx    = (nxf+1)/(nxc+1) * (i+1) - 1
      
                       DO iif = bottomx, topx
      
                           eps    = ( iif * dxf + 0.5 * dxf - i * dxc - 0.5 * dxc ) / dxc
                           alpha  = ( ( dxf / dxc )**2.0 - 1.0 ) / 24.0
                           eminus = eps * ( eps - 1.0 ) / 2.0 + alpha
                           edot  = ( 1.0 - eps**2.0 ) - 2.0 * alpha
                           eplus  = eps * ( eps + 1.0 ) / 2.0 + alpha
      
                           ptprf(k,j,iif) = eminus * work3d(k,j,i-1) &
                               + edot  * work3d(k,j,i)   &
                               + eplus  * work3d(k,j,i+1)
                       END DO
      
                   END DO
               END DO
           END DO
      
           !
           !-- Interpolation in y-direction (quadratic, Clark and Farley)
      
           DO k = nzbottom, nztop+2
               DO j = bdims_rem(2,1), bdims_rem(2,2)
        
                   bottomy = (nyf+1)/(nyc+1) * j
                   topy    = (nyf+1)/(nyc+1) * (j+1) - 1
      
                   DO iif = nxl, nxr
                       DO jjf = bottomy, topy
      
                           eps    = ( jjf * dyf + 0.5 * dyf - j * dyc - 0.5 * dyc ) / dyc
                           alpha  = ( ( dyf / dyc )**2.0 - 1.0 ) / 24.0
                           eminus = eps * ( eps - 1.0 ) / 2.0 + alpha
                           edot  = ( 1.0 - eps**2.0 ) - 2.0 * alpha
                           eplus  = eps * ( eps + 1.0 ) / 2.0 + alpha
      
                           ptprs(k,jjf,iif) = eminus * ptprf(k,j-1,iif) &
                               + edot  * ptprf(k,j,iif)   &
                               + eplus  * ptprf(k,j+1,iif)
      
                       END DO
                   END DO
      
               END DO
           END DO
      
           !
           !-- Interpolation in z-direction (quadratic, Clark and Farley)
      
           DO k = nzbottom+1, nztop
            
               bottomz = (dzc/dzf) * (k-1) + 1
               topz    = (dzc/dzf) * k
      
               DO jjf = nys, nyn
                   DO iif = nxl, nxr
                       DO kkf = bottomz, topz
                       
                           eps    = ( zuf(kkf) - zuc(k) ) / dzc
                           alpha  = ( ( dzf / dzc )**2.0 - 1.0 ) / 24.0
                           eminus = eps * ( eps - 1.0 ) / 2.0 + alpha
                           edot  = ( 1.0 - eps**2.0 ) - 2.0 * alpha
                           eplus  = eps * ( eps + 1.0 ) / 2.0 + alpha
       
                           kh(kkf,jjf,iif) = eminus * ptprs(k-1,jjf,iif) &
                               + edot  * ptprs(k,jjf,iif)   &
                               + eplus  * ptprs(k+1,jjf,iif)
      
                       END DO
                   END DO
               END DO
      
           END DO
               
           DO jjf = nys, nyn
               DO iif = nxl, nxr
      
                   eps    = ( zuf(nzt+1) - zuc(nztop+1) ) / dzc
                   alpha  = ( ( dzf / dzc )**2.0 - 1.0 ) / 24.0
                   eminus = eps * ( eps - 1.0 ) / 2.0 + alpha
                   edot  = ( 1.0 - eps**2.0 ) - 2.0 * alpha
                   eplus  = eps * ( eps + 1.0 ) / 2.0 + alpha
      
                   kh(nzt+1,jjf,iif) = eminus * ptprs(nztop,jjf,iif)   &
                       + edot  * ptprs(nztop+1,jjf,iif) &
                       + eplus  * ptprs(nztop+2,jjf,iif)
      
               END DO
           END DO
      
      
           DEALLOCATE( ptprf, ptprs ) 
      
       END SUBROUTINE interpolate_to_fine_kh
      
       SUBROUTINE interpolate_to_fine_km
      
      
           USE arrays_3d
           USE control_parameters
           USE grid_variables
           USE indices
           USE pegrid
           
      
           IMPLICIT NONE

           INTEGER(iwp)                       :: i
           INTEGER(iwp)                       :: j
           INTEGER(iwp)                       :: k
           INTEGER(iwp)                       :: iif
           INTEGER(iwp)                       :: jjf
           INTEGER(iwp)                       :: kkf
           INTEGER(iwp)                       :: nzbottom
           INTEGER(iwp)                       :: nztop
           INTEGER(iwp)                       :: bottomx
           INTEGER(iwp)                       :: bottomy
           INTEGER(iwp)                       :: bottomz
           INTEGER(iwp)                       :: topx
           INTEGER(iwp)                       :: topy
           INTEGER(iwp)                       :: topz
           REAL(wp)                           :: eps
           REAL(wp)                           :: alpha
           REAL(wp)                           :: eminus
           REAL(wp)                           :: edot
           REAL(wp)                           :: eplus
           REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: ptprs
           REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: ptprf
      
      
           nzbottom = bdims_rem (3,1)
           nztop  = bdims_rem (3,2)
           !   nztop  = blk_dim_rem (3,2)+1
      
      
           ALLOCATE( ptprf(nzbottom:nztop+2,bdims_rem(2,1)-1:bdims_rem(2,2)+1,nxl:nxr) )
           ALLOCATE( ptprs(nzbottom:nztop+2,nys:nyn,nxl:nxr) )
      
      
           !
           !-- Initialisation of scalar variables
      
           !
           !-- Interpolation in x-direction (quadratic, Clark and Farley)
      
           DO k = nzbottom, nztop+2
               DO j = bdims_rem(2,1)-1, bdims_rem(2,2)+1
                   DO i = bdims_rem(1,1), bdims_rem(1,2)
      
                       bottomx = (nxf+1)/(nxc+1) * i
                       topx    = (nxf+1)/(nxc+1) * (i+1) - 1
      
                       DO iif = bottomx, topx
      
                           eps    = ( iif * dxf + 0.5 * dxf - i * dxc - 0.5 * dxc ) / dxc
                           alpha  = ( ( dxf / dxc )**2.0 - 1.0 ) / 24.0
                           eminus = eps * ( eps - 1.0 ) / 2.0 + alpha
                           edot  = ( 1.0 - eps**2.0 ) - 2.0 * alpha
                           eplus  = eps * ( eps + 1.0 ) / 2.0 + alpha
      
                           ptprf(k,j,iif) = eminus * work3d(k,j,i-1) &
                               + edot  * work3d(k,j,i)   &
                               + eplus  * work3d(k,j,i+1)
                       END DO
      
                   END DO
               END DO
           END DO
      
           !
           !-- Interpolation in y-direction (quadratic, Clark and Farley)
      
           DO k = nzbottom, nztop+2
               DO j = bdims_rem(2,1), bdims_rem(2,2)
        
                   bottomy = (nyf+1)/(nyc+1) * j
                   topy    = (nyf+1)/(nyc+1) * (j+1) - 1
      
                   DO iif = nxl, nxr
                       DO jjf = bottomy, topy
      
                           eps    = ( jjf * dyf + 0.5 * dyf - j * dyc - 0.5 * dyc ) / dyc
                           alpha  = ( ( dyf / dyc )**2.0 - 1.0 ) / 24.0
                           eminus = eps * ( eps - 1.0 ) / 2.0 + alpha
                           edot  = ( 1.0 - eps**2.0 ) - 2.0 * alpha
                           eplus  = eps * ( eps + 1.0 ) / 2.0 + alpha
      
                           ptprs(k,jjf,iif) = eminus * ptprf(k,j-1,iif) &
                               + edot  * ptprf(k,j,iif)   &
                               + eplus  * ptprf(k,j+1,iif)
      
                       END DO
                   END DO
      
               END DO
           END DO
      
           !
           !-- Interpolation in z-direction (quadratic, Clark and Farley)
      
           DO k = nzbottom+1, nztop
            
               bottomz = (dzc/dzf) * (k-1) + 1
               topz    = (dzc/dzf) * k
      
               DO jjf = nys, nyn
                   DO iif = nxl, nxr
                       DO kkf = bottomz, topz
                       
                           eps    = ( zuf(kkf) - zuc(k) ) / dzc
                           alpha  = ( ( dzf / dzc )**2.0 - 1.0 ) / 24.0
                           eminus = eps * ( eps - 1.0 ) / 2.0 + alpha
                           edot  = ( 1.0 - eps**2.0 ) - 2.0 * alpha
                           eplus  = eps * ( eps + 1.0 ) / 2.0 + alpha
       
                           km(kkf,jjf,iif) = eminus * ptprs(k-1,jjf,iif) &
                               + edot  * ptprs(k,jjf,iif)   &
                               + eplus  * ptprs(k+1,jjf,iif)
      
                       END DO
                   END DO
               END DO
      
           END DO
               
           DO jjf = nys, nyn
               DO iif = nxl, nxr
      
                   eps    = ( zuf(nzt+1) - zuc(nztop+1) ) / dzc
                   alpha  = ( ( dzf / dzc )**2.0 - 1.0 ) / 24.0
                   eminus = eps * ( eps - 1.0 ) / 2.0 + alpha
                   edot  = ( 1.0 - eps**2.0 ) - 2.0 * alpha
                   eplus  = eps * ( eps + 1.0 ) / 2.0 + alpha
      
                   km(nzt+1,jjf,iif)  = eminus * ptprs(nztop,jjf,iif)   &
                       + edot  * ptprs(nztop+1,jjf,iif) &
                       + eplus  * ptprs(nztop+2,jjf,iif)
      
               END DO
           END DO
      
      
           DEALLOCATE( ptprf, ptprs ) 
      
       END SUBROUTINE interpolate_to_fine_km
      
      
      
      
!       SUBROUTINE interpolate_to_fine_flux
!      
!      
!           USE arrays_3d
!           USE control_parameters
!           USE grid_variables
!           USE indices
!           USE pegrid
!           
!      
!           IMPLICIT NONE
!
!           INTEGER(iwp)                       :: i
!           INTEGER(iwp)                       :: j
!           INTEGER(iwp)                       :: iif
!           INTEGER(iwp)                       :: jjf
!           INTEGER(iwp)                       :: bottomx
!           INTEGER(iwp)                       :: bottomy
!           INTEGER(iwp)                       :: topx
!           INTEGER(iwp)                       :: topy
!           REAL(wp)                           :: eps
!           REAL(wp)                           :: alpha
!           REAL(wp)                           :: eminus
!           REAL(wp)                           :: edot
!           REAL(wp)                           :: eplus
!           REAL(wp), DIMENSION(:,:), ALLOCATABLE   :: uswspr
!           REAL(wp), DIMENSION(:,:), ALLOCATABLE   :: vswspr
!           REAL(wp), DIMENSION(:,:), ALLOCATABLE   :: tspr 
!           REAL(wp), DIMENSION(:,:), ALLOCATABLE   :: uspr
!      
!           ALLOCATE( uswspr(bdims_rem(2,1)-1:bdims_rem(2,2)+1,nxl:nxr) )
!           ALLOCATE( vswspr(bdims_rem(2,1)-1:bdims_rem(2,2)+1,nxl:nxr) )
!           ALLOCATE( tspr  (bdims_rem(2,1)-1:bdims_rem(2,2)+1,nxl:nxr) )
!           ALLOCATE( uspr  (bdims_rem(2,1)-1:bdims_rem(2,2)+1,nxl:nxr) )
!      
!           !
!           !-- Initialisation of scalar variables (2D)
!      
!           !
!           !-- Interpolation in x-direction (quadratic, Clark and Farley)
!      
!           DO j = bdims_rem(2,1)-1, bdims_rem(2,2)+1
!               DO i = bdims_rem(1,1), bdims_rem(1,2)
!              
!                   bottomx = (nxf+1)/(nxc+1) * i
!                   topx    = (nxf+1)/(nxc+1) * (i+1) - 1
!      
!                   DO iif = bottomx, topx
!      
!                       eps    = ( iif * dxf + 0.5 * dxf - i * dxc - 0.5 * dxc ) / dxc
!                       alpha  = ( ( dxf / dxc )**2.0 - 1.0 ) / 24.0
!                       eminus = eps * ( eps - 1.0 ) / 2.0 + alpha
!                       edot  = ( 1.0 - eps**2.0 ) - 2.0 * alpha
!                       eplus  = eps * ( eps + 1.0 ) / 2.0 + alpha
!      
!                       uswspr(j,iif) = eminus * work2dusws(j,i-1) &
!                           + edot  * work2dusws(j,i)   &
!                           + eplus  * work2dusws(j,i+1)
!      
!                       vswspr(j,iif) = eminus * work2dvsws(j,i-1) &
!                           + edot  * work2dvsws(j,i)   &
!                           + eplus  * work2dvsws(j,i+1)
!      
!                       tspr(j,iif)   = eminus * work2dts(j,i-1) &
!                           + edot  * work2dts(j,i)   &
!                           + eplus  * work2dts(j,i+1)
!      
!                       uspr(j,iif)   = eminus * work2dus(j,i-1) &
!                           + edot  * work2dus(j,i)   &
!                           + eplus  * work2dus(j,i+1)
!      
!                   END DO
!      
!               END DO
!           END DO
!      
!           !
!           !-- Interpolation in y-direction (quadratic, Clark and Farley)
!      
!           DO j = bdims_rem(2,1), bdims_rem(2,2)
!             
!               bottomy = (nyf+1)/(nyc+1) * j
!               topy    = (nyf+1)/(nyc+1) * (j+1) - 1
!      
!               DO iif = nxl, nxr
!                   DO jjf = bottomy, topy
!      
!                       eps    = ( jjf * dyf + 0.5 * dyf - j * dyc - 0.5 * dyc ) / dyc
!                       alpha  = ( ( dyf / dyc )**2.0 - 1.0 ) / 24.0
!                       eminus = eps * ( eps - 1.0 ) / 2.0 + alpha
!                       edot  = ( 1.0 - eps**2.0 ) - 2.0 * alpha
!                       eplus  = eps * ( eps + 1.0 ) / 2.0 + alpha
! 
!!
!!--   TODO
!--    variables are not compatible with the new surface layer module
!   
!                    surf_def_h(0)%usws(jjf,iif) = eminus * uswspr(j-1,if) &
!                        + edot  * uswspr(j,iif)   &
!                        + eplus  * uswspr(j+1,iif)
!   
!                    surf_def_h(0)%vsws(jjf,iif) = eminus * vswspr(j-1,if) &
!                        + edot  * vswspr(j,iif)   &
!                        + eplus  * vswspr(j+1,iif)
!   
!                    ts(jjf,iif)   = eminus * tspr(j-1,if) &
!                        + edot  * tspr(j,iif)   &
!                        + eplus  * tspr(j+1,iif)
!
!                    us(jjf,iif)   = eminus * uspr(j-1,if) &
!                        + edot  * uspr(j,iif)   &
!                        + eplus  * uspr(j+1,iif)
!      
!                   END DO
!               END DO
!      
!           END DO
!      
!      
!           DEALLOCATE( uswspr, vswspr )
!           DEALLOCATE( tspr, uspr )
!      
!      
!       END SUBROUTINE interpolate_to_fine_flux
   
   
#endif       
    END SUBROUTINE vnest_init_fine
   
    SUBROUTINE vnest_boundary_conds
#if defined( __parallel )
        !------------------------------------------------------------------------------!
        ! Description:
        ! ------------
        ! Boundary conditions for the prognostic quantities.
        ! One additional bottom boundary condition is applied for the TKE (=(u*)**2)
        ! in prandtl_fluxes. The cyclic lateral boundary conditions are implicitly
        ! handled in routine exchange_horiz. Pressure boundary conditions are
        ! explicitly set in routines pres, poisfft, poismg and sor.
        !------------------------------------------------------------------------------!
    
        USE arrays_3d
        USE control_parameters
        USE grid_variables
        USE indices
        USE pegrid
        
    
        IMPLICIT NONE
    
        INTEGER(iwp)                          :: i
        INTEGER(iwp)                          :: j
        INTEGER(iwp)                          :: iif
        INTEGER(iwp)                          :: jjf
    
    
        !
        !-- vnest: top boundary conditions
    
        IF ( coupling_mode == 'vnested_crse' )  THEN
        !-- Send data to fine grid for TOP BC
    
            offset(1) = ( pdims_partner(1) / pdims(1) ) * pcoord(1)
            offset(2) = ( pdims_partner(2) / pdims(2) ) * pcoord(2)
    
            do j = 0,   ( pdims_partner(2) / pdims(2) ) - 1
                do i = 0,   ( pdims_partner(1) / pdims(1) ) - 1
                    map_coord(1) = i+offset(1)
                    map_coord(2) = j+offset(2)
    
                    target_idex = f_rnk_lst(map_coord(1),map_coord(2)) + numprocs
      
                    bdims (1,1) =  c2f_dims_cg (0,map_coord(1),map_coord(2))
                    bdims (1,2) =  c2f_dims_cg (1,map_coord(1),map_coord(2))
                    bdims (2,1) =  c2f_dims_cg (2,map_coord(1),map_coord(2))
                    bdims (2,2) =  c2f_dims_cg (3,map_coord(1),map_coord(2))
                    bdims (3,1) =  c2f_dims_cg (4,map_coord(1),map_coord(2))
                    bdims (3,2) =  c2f_dims_cg (5,map_coord(1),map_coord(2))
  
                    n_cell_c =      ( (bdims(1,2)-bdims(1,1)) + 3 ) * &
                        ( (bdims(2,2)-bdims(2,1)) + 3 ) * &
                        ( (bdims(3,2)-bdims(3,1)) + 1 )
    
                    CALL MPI_SEND(u (bdims(3,1), bdims(2,1)-1, bdims(1,1)-1), &
                        1, TYPE_VNEST_BC, target_idex,    &
                        201,    comm_inter, ierr)
    
                    CALL MPI_SEND(v(bdims(3,1), bdims(2,1)-1, bdims(1,1)-1),&
                        1, TYPE_VNEST_BC, target_idex,    &
                        202,    comm_inter, ierr)

                   CALL MPI_SEND(w(bdims(3,1), bdims(2,1)-1, bdims(1,1)-1),&
                        1, TYPE_VNEST_BC, target_idex,    &
                        203,    comm_inter, ierr)
    
                    CALL MPI_SEND(pt(bdims(3,1), bdims(2,1)-1, bdims(1,1)-1),&
                        1, TYPE_VNEST_BC, target_idex,    &
                        205,    comm_inter, ierr)
    
                    IF ( humidity )  THEN 
                    CALL MPI_SEND(q(bdims(3,1), bdims(2,1)-1, bdims(1,1)-1),&
                        1, TYPE_VNEST_BC, target_idex,    &
                        209,    comm_inter, ierr)
                    ENDIF 
 
                end do
            end do
    
        ELSEIF ( coupling_mode == 'vnested_fine' )  THEN
        !-- Receive data from coarse grid for TOP BC
    
            offset(1) = pcoord(1) / ( pdims(1)/pdims_partner(1) )
            offset(2) = pcoord(2) / ( pdims(2)/pdims_partner(2) )
            map_coord(1) = offset(1)
            map_coord(2) = offset(2)
            target_idex = c_rnk_lst(map_coord(1),map_coord(2))
    
            bdims_rem (1,1) = c2f_dims_fg(0)
            bdims_rem (1,2) = c2f_dims_fg(1)
            bdims_rem (2,1) = c2f_dims_fg(2)
            bdims_rem (2,2) = c2f_dims_fg(3)
            bdims_rem (3,1) = c2f_dims_fg(4)
            bdims_rem (3,2) = c2f_dims_fg(5)
 
            n_cell_c =                                    &
                ( (bdims_rem(1,2)-bdims_rem(1,1)) + 3 ) * &
                ( (bdims_rem(2,2)-bdims_rem(2,1)) + 3 ) * &
                ( (bdims_rem(3,2)-bdims_rem(3,1)) + 1 )
    
            ALLOCATE( work3d  (                    &
                bdims_rem(3,1)  :bdims_rem(3,2)  , &
                bdims_rem(2,1)-1:bdims_rem(2,2)+1, &
                bdims_rem(1,1)-1:bdims_rem(1,2)+1))
    
    
            CALL MPI_RECV( work3d ,n_cell_c, MPI_REAL, target_idex, 201, &
                comm_inter,status, ierr )
            interpol3d => u
            call vnest_set_topbc_u
    
            CALL MPI_RECV( work3d ,n_cell_c, MPI_REAL, target_idex, 202, &
                comm_inter,status, ierr )
            interpol3d => v
            call vnest_set_topbc_v
    
            CALL MPI_RECV( work3d ,n_cell_c, MPI_REAL, target_idex, 203, &
                comm_inter,status, ierr )
            interpol3d => w
            call vnest_set_topbc_w
    
            CALL MPI_RECV( work3d,n_cell_c, MPI_REAL, target_idex, 205, &
                comm_inter,status, ierr )
            interpol3d => pt
            call vnest_set_topbc_s
    
            IF ( humidity )  THEN
               CALL MPI_RECV( work3d,n_cell_c, MPI_REAL, target_idex, 209, &
                        comm_inter,status, ierr )
                interpol3d => q
                call vnest_set_topbc_s
    
                CALL exchange_horiz_2d(q (nzt+1,:,:) )
            ENDIF
    
!-- TKE Neumann BC for FG top
            DO jjf = nys, nyn
                DO iif = nxl, nxr
                   e(nzt+1,jjf,iif) = e(nzt,jjf,iif)
                END DO
            END DO

!
!-- w level nzt+1 does not impact results. Only to avoid jumps while
!-- plotting profiles
            w(nzt+1,:,:) = w(nzt,:,:)

            CALL exchange_horiz_2d(u (nzt+1,:,:) )
            CALL exchange_horiz_2d(v (nzt+1,:,:) )
            CALL exchange_horiz_2d(pt(nzt+1,:,:) )
            CALL exchange_horiz_2d(e (nzt+1,:,:) )
            CALL exchange_horiz_2d(w (nzt+1,:,:) )
            CALL exchange_horiz_2d(w (nzt  ,:,:) )
    
            NULLIFY    (   interpol3d  )
            DEALLOCATE (   work3d      )
    
        ENDIF
    
    
    CONTAINS
    
       SUBROUTINE vnest_set_topbc_w
       
       
           USE arrays_3d
           USE control_parameters
           USE grid_variables
           USE indices
           USE pegrid
           
       
           IMPLICIT NONE

           INTEGER(iwp)                       :: i
           INTEGER(iwp)                       :: j
           INTEGER(iwp)                       :: iif
           INTEGER(iwp)                       :: jjf
           INTEGER(iwp)                       :: bottomx
           INTEGER(iwp)                       :: bottomy
           INTEGER(iwp)                       :: topx
           INTEGER(iwp)                       :: topy
           REAL(wp)                           :: eps
           REAL(wp)                           :: alpha
           REAL(wp)                           :: eminus
           REAL(wp)                           :: edot
           REAL(wp)                           :: eplus
           REAL(wp), DIMENSION(:,:), ALLOCATABLE :: wprf
       
           
           ALLOCATE( wprf(bdims_rem(2,1)-1:bdims_rem(2,2)+1,nxl:nxr) )
       
           !
           !-- Determination of a boundary condition for the vertical velocity component w:
           !-- In this case only interpolation in x- and y- direction is necessary, as the
           !-- boundary w-node of the fine grid coincides with a w-node in the coarse grid.
           !-- For both interpolations the scheme of Clark and Farley is used.
       
           !
           !-- Interpolation in x-direction
       
           DO j = bdims_rem(2,1)-1, bdims_rem(2,2)+1
       
               DO i = bdims_rem(1,1), bdims_rem(1,2)
       
                   bottomx = (nxf+1)/(nxc+1) * i
                   topx    = (nxf+1)/(nxc+1) * (i+1) - 1
       
                   DO iif = bottomx, topx
       
                       eps = (iif * dxf + 0.5 * dxf - i * dxc - 0.5 * dxc) / dxc
                       alpha = ( (dxf/dxc)**2.0 - 1.0) / 24.0
                       eminus = eps * ( eps - 1.0 ) / 2.0 + alpha
                       edot  = ( 1.0 - eps**2.0 ) - 2.0 * alpha
                       eplus  = eps * ( eps + 1.0 ) / 2.0 + alpha
                       wprf(j,iif) = eminus * work3d(bdims_rem(3,1),j,i-1) &
                           + edot  * work3d(bdims_rem(3,1),j,i)   &
                           + eplus  * work3d(bdims_rem(3,1),j,i+1)
       
                   END DO
       
               END DO
       
           END DO
              
           !
           !-- Interpolation in y-direction
       
           DO j = bdims_rem(2,1), bdims_rem(2,2)
       
               bottomy = (nyf+1)/(nyc+1) * j
               topy    = (nyf+1)/(nyc+1) * (j+1) - 1
       
               DO iif = nxl, nxr
       
                   DO jjf = bottomy, topy
       
                       eps = (jjf * dyf + 0.5 * dyf - j * dyc - 0.5 * dyc) / dyc
       
                       alpha = ( (dyf/dyc)**2.0 - 1.0) / 24.0
       
                       eminus = eps * ( eps - 1.0 ) / 2.0 + alpha
       
                       edot = ( 1.0 - eps**2.0 ) - 2.0 * alpha
       
                       eplus = eps * ( eps + 1.0 ) / 2.0 + alpha
       
                       w(nzt,jjf,iif) = eminus * wprf(j-1,iif) &
                           + edot  * wprf(j,iif)   &
                           + eplus  * wprf(j+1,iif)
       
                   END DO
       
               END DO
       
           END DO

           DEALLOCATE( wprf )
       
       END SUBROUTINE vnest_set_topbc_w
       
       
       SUBROUTINE vnest_set_topbc_u
       
       
           USE arrays_3d
           USE control_parameters
           USE grid_variables
           USE indices
           USE pegrid
           
       
           IMPLICIT NONE

           INTEGER(iwp)                       :: i
           INTEGER(iwp)                       :: j
           INTEGER(iwp)                       :: k
           INTEGER(iwp)                       :: iif
           INTEGER(iwp)                       :: jjf
           INTEGER(iwp)                       :: bottomx
           INTEGER(iwp)                       :: bottomy
           INTEGER(iwp)                       :: topx
           INTEGER(iwp)                       :: topy
           REAL(wp)                           :: eps
           REAL(wp)                           :: alpha
           REAL(wp)                           :: eminus
           REAL(wp)                           :: edot
           REAL(wp)                           :: eplus
           REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: uprf
           REAL(wp), DIMENSION(:,:), ALLOCATABLE   :: uprs
       
           ALLOCATE( uprf(bdims_rem(3,1):bdims_rem(3,2),nys:nyn,bdims_rem(1,1)-1:bdims_rem(1,2)+1) )
           ALLOCATE( uprs(nys:nyn,bdims_rem(1,1)-1:bdims_rem(1,2)+1) )
       
       
           !
           !-- Interpolation in y-direction
       
           DO k = bdims_rem(3,1), bdims_rem(3,2)
               DO j = bdims_rem(2,1), bdims_rem(2,2)
       
                   bottomy = (nyf+1)/(nyc+1) * j
                   topy    = (nyf+1)/(nyc+1) * (j+1) - 1
       
                   DO i = bdims_rem(1,1)-1, bdims_rem(1,2)+1
                       DO jjf = bottomy, topy
       
                           eps = (jjf * dyf + 0.5 * dyf - j * dyc - 0.5 * dyc) / dyc
                           alpha = ( (dyf/dyc)**2.0 - 1.0) / 24.0
                           eminus = eps * ( eps - 1.0 ) / 2.0 + alpha
                           edot  = ( 1.0 - eps**2.0 ) - 2.0 * alpha
                           eplus  = eps * ( eps + 1.0 ) / 2.0 + alpha
       
                           uprf(k,jjf,i) = eminus * work3d(k,j-1,i) &
                               + edot  * work3d(k,j,i)   &
                               + eplus  * work3d(k,j+1,i)
                       END DO
                   END DO
       
               END DO
           END DO
       
           !
           !-- Interpolation in z-direction
       
           DO jjf = nys, nyn
               DO i = bdims_rem(1,1)-1, bdims_rem(1,2)+1
                   eps = ( zuf(nzt+1) - zuc(bdims_rem(3,1)+1) ) / dzc
                   alpha = ( (dzf/dzc)**2.0 - 1.0) / 24.0
                   eminus = eps * ( eps - 1.0 ) / 2.0 + alpha
                   edot  = ( 1.0 - eps**2.0 ) - 2.0 * alpha
                   eplus  = eps * ( eps + 1.0 ) / 2.0 + alpha
                   uprs(jjf,i) = eminus * uprf(bdims_rem(3,1),jjf,i)   &
                       + edot  * uprf(bdims_rem(3,1)+1,jjf,i) &
                       + eplus  * uprf(bdims_rem(3,1)+2,jjf,i)
               END DO
           END DO
       
           !
           !-- Interpolation in x-direction
       
           DO jjf = nys, nyn
               DO i = bdims_rem(1,1), bdims_rem(1,2)
       
                   bottomx = (nxf+1)/(nxc+1) * i
                   topx    = (nxf+1)/(nxc+1) * (i+1) - 1
       
                   DO iif = bottomx, topx
                       u(nzt+1,jjf,iif) = uprs(jjf,i) + ( iif * dxf - i * dxc ) * ( uprs(jjf,i+1) - uprs(jjf,i) ) / dxc
                   END DO
       
               END DO
           END DO           
       
       
       
           DEALLOCATE ( uprf, uprs )
       
       END SUBROUTINE vnest_set_topbc_u
       
       
       SUBROUTINE vnest_set_topbc_v
       
       
           USE arrays_3d
           USE control_parameters
           USE grid_variables
           USE indices
           USE pegrid
           
       
           IMPLICIT NONE

           INTEGER(iwp)                       :: i
           INTEGER(iwp)                       :: j
           INTEGER(iwp)                       :: k
           INTEGER(iwp)                       :: iif
           INTEGER(iwp)                       :: jjf
           INTEGER(iwp)                       :: bottomx
           INTEGER(iwp)                       :: bottomy
           INTEGER(iwp)                       :: topx
           INTEGER(iwp)                       :: topy
           REAL(wp)                           :: eps
           REAL(wp)                           :: alpha
           REAL(wp)                           :: eminus
           REAL(wp)                           :: edot
           REAL(wp)                           :: eplus
           REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: vprf
           REAL(wp), DIMENSION(:,:), ALLOCATABLE   :: vprs
         
           
       
           ALLOCATE( vprf(bdims_rem(3,1):bdims_rem(3,2),bdims_rem(2,1)-1:bdims_rem(2,2)+1,nxl:nxr) )
           ALLOCATE( vprs(bdims_rem(2,1)-1:bdims_rem(2,2)+1,nxl:nxr) )
           !
           !-- Determination of a boundary condition for the horizontal velocity component v:
           !-- Interpolation in x- and z-direction is carried out by using the scheme,
           !-- which was derived by Clark and Farley (1984). In y-direction a
           !-- linear interpolation is carried out.
       
           !
           !-- Interpolation in x-direction
       
           DO k = bdims_rem(3,1), bdims_rem(3,2)
               DO j = bdims_rem(2,1)-1, bdims_rem(2,2)+1
                   DO i = bdims_rem(1,1), bdims_rem(1,2)
       
                       bottomx = (nxf+1)/(nxc+1) * i
                       topx    = (nxf+1)/(nxc+1) * (i+1) - 1
       
                       DO iif = bottomx, topx
       
                           eps = (iif * dxf + 0.5 * dxf - i * dxc - 0.5 * dxc) / dxc
                           alpha = ( (dxf/dxc)**2.0 - 1.0) / 24.0
                           eminus = eps * ( eps - 1.0 ) / 2.0 + alpha
                           edot  = ( 1.0 - eps**2.0 ) - 2.0 * alpha
                           eplus  = eps * ( eps + 1.0 ) / 2.0 + alpha
                           vprf(k,j,iif) = eminus * work3d(k,j,i-1) &
                               + edot  * work3d(k,j,i)   &
                               + eplus  * work3d(k,j,i+1)
                       END DO
       
                   END DO
               END DO
           END DO
       
           !
           !-- Interpolation in z-direction
       
           DO j = bdims_rem(2,1)-1, bdims_rem(2,2)+1
               DO iif = nxl, nxr
       
                   eps = ( zuf(nzt+1) - zuc(bdims_rem(3,1)+1) ) / dzc
                   alpha = ( (dzf/dzc)**2.0 - 1.0) / 24.0
                   eminus = eps * ( eps - 1.0 ) / 2.0 + alpha
                   edot  = ( 1.0 - eps**2.0 ) - 2.0 * alpha
                   eplus  = eps * ( eps + 1.0 ) / 2.0 + alpha
                   vprs(j,iif) = eminus * vprf(bdims_rem(3,1),j,iif)   &
                       + edot  * vprf(bdims_rem(3,1)+1,j,iif) &
                       + eplus  * vprf(bdims_rem(3,1)+2,j,iif)
       
               END DO
           END DO
       
           !
           !-- Interpolation in y-direction
       
           DO j = bdims_rem(2,1), bdims_rem(2,2)
               DO iif = nxl, nxr
       
                   bottomy = (nyf+1)/(nyc+1) * j
                   topy    = (nyf+1)/(nyc+1) * (j+1) - 1
       
                   DO jjf = bottomy, topy
       
                       v(nzt+1,jjf,iif) = vprs(j,iif) + ( jjf * dyf - j * dyc ) * ( vprs(j+1,iif) - vprs(j,iif) ) / dyc
       
                   END DO
               END DO
           END DO      
       
       
           DEALLOCATE ( vprf, vprs)
       
       
       
       END SUBROUTINE vnest_set_topbc_v
       
       
       SUBROUTINE vnest_set_topbc_s
       
       
           USE arrays_3d
           USE control_parameters
           USE grid_variables
           USE indices
           USE pegrid
           
       
           IMPLICIT NONE

           INTEGER(iwp)                       :: i
           INTEGER(iwp)                       :: j
           INTEGER(iwp)                       :: k
           INTEGER(iwp)                       :: iif
           INTEGER(iwp)                       :: jjf
           INTEGER(iwp)                       :: bottomx
           INTEGER(iwp)                       :: bottomy
           INTEGER(iwp)                       :: topx
           INTEGER(iwp)                       :: topy
           REAL(wp)                           :: eps
           REAL(wp)                           :: alpha
           REAL(wp)                           :: eminus
           REAL(wp)                           :: edot
           REAL(wp)                           :: eplus
           REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: ptprf
           REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: ptprs
       
           
       
           ALLOCATE( ptprf(bdims_rem(3,1):bdims_rem(3,2),bdims_rem(2,1)-1:bdims_rem(2,2)+1,nxl:nxr) )
           ALLOCATE( ptprs(bdims_rem(3,1):bdims_rem(3,2),nys:nyn,nxl:nxr) )
       
           !
           !-- Determination of a boundary condition for the potential temperature pt:
           !-- The scheme derived by Clark and Farley can be used in all three dimensions.
       
           !
           !-- Interpolation in x-direction
       
           DO k = bdims_rem(3,1), bdims_rem(3,2)
       
               DO j = bdims_rem(2,1)-1, bdims_rem(2,2)+1
       
                   DO i = bdims_rem(1,1), bdims_rem(1,2)
       
                       bottomx = (nxf+1)/(nxc+1) * i
                       topx    = (nxf+1)/(nxc+1) *(i+1) - 1
       
                       DO iif = bottomx, topx
       
                           eps = (iif * dxf + 0.5 * dxf - i * dxc - 0.5 * dxc) / dxc
       
                           alpha = ( (dxf/dxc)**2.0 - 1.0) / 24.0
       
                           eminus = eps * (eps - 1.0 ) / 2.0 + alpha
       
                           edot = ( 1.0 - eps**2.0 ) - 2.0 * alpha
       
                           eplus = eps * ( eps + 1.0 ) / 2.0 + alpha
       
                           ptprf(k,j,iif) = eminus * work3d(k,j,i-1) &
                               + edot  * work3d(k,j,i)   &
                               + eplus  * work3d(k,j,i+1)
                       END DO
       
                   END DO
       
               END DO
       
           END DO
       
           !
           !-- Interpolation in y-direction
       
           DO k = bdims_rem(3,1), bdims_rem(3,2)
       
               DO j = bdims_rem(2,1), bdims_rem(2,2)
       
                   bottomy = (nyf+1)/(nyc+1) * j
                   topy    = (nyf+1)/(nyc+1) * (j+1) - 1
       
                   DO iif = nxl, nxr
       
                       DO jjf = bottomy, topy
       
                           eps = (jjf * dyf + 0.5 * dyf - j * dyc - 0.5 * dyc) / dyc
       
                           alpha = ( (dyf/dyc)**2.0 - 1.0) / 24.0
                   
                           eminus = eps * (eps - 1.0) / 2.0 + alpha
       
                           edot  = ( 1.0 - eps**2.0 ) - 2.0 * alpha
       
                           eplus  = eps * ( eps + 1.0 ) / 2.0 + alpha
       
                           ptprs(k,jjf,iif) = eminus * ptprf(k,j-1,iif) &
                               + edot  * ptprf(k,j,iif)   &
                               + eplus  * ptprf(k,j+1,iif)
                       END DO
       
                   END DO
       
               END DO
       
           END DO
       
           !
           !-- Interpolation in z-direction
       
           DO jjf = nys, nyn
               DO iif = nxl, nxr
       
                   eps = ( zuf(nzt+1) - zuc(bdims_rem(3,1)+1) ) / dzc
       
                   alpha = ( (dzf/dzc)**2.0 - 1.0) / 24.0
       
                   eminus = eps * ( eps - 1.0 ) / 2.0 + alpha
       
                   edot  = ( 1.0 - eps**2.0 ) - 2.0 * alpha
       
                   eplus  = eps * ( eps + 1.0 ) / 2.0 + alpha
       
                   interpol3d (nzt+1,jjf,iif) = eminus * ptprs(bdims_rem(3,1),jjf,iif)   &
                       + edot  * ptprs(bdims_rem(3,1)+1,jjf,iif) &
                       + eplus  * ptprs(bdims_rem(3,1)+2,jjf,iif)
       
               END DO
           END DO          
           
           DEALLOCATE ( ptprf, ptprs )
       
       
       
       END SUBROUTINE vnest_set_topbc_s
#endif
    END SUBROUTINE vnest_boundary_conds
    
 
    SUBROUTINE vnest_boundary_conds_khkm
#if defined( __parallel )
    
        !--------------------------------------------------------------------------------!
        ! Description:
        ! ------------
        ! Boundary conditions for the prognostic quantities.
        ! One additional bottom boundary condition is applied for the TKE (=(u*)**2)
        ! in prandtl_fluxes. The cyclic lateral boundary conditions are implicitly
        ! handled in routine exchange_horiz. Pressure boundary conditions are
        ! explicitly set in routines pres, poisfft, poismg and sor.
        !------------------------------------------------------------------------------!
    
        USE arrays_3d
        USE control_parameters
        USE grid_variables
        USE indices
        USE pegrid
        
    
        IMPLICIT NONE
    
        INTEGER(iwp)                          :: i
        INTEGER(iwp)                          :: j
        INTEGER(iwp)                          :: iif
        INTEGER(iwp)                          :: jjf
    
    
        IF ( coupling_mode == 'vnested_crse' )  THEN
        ! Send data to fine grid for TOP BC
    
            offset(1) = ( pdims_partner(1) / pdims(1) ) * pcoord(1)
            offset(2) = ( pdims_partner(2) / pdims(2) ) * pcoord(2)
    
            do j = 0,   ( pdims_partner(2) / pdims(2) ) - 1
                do i = 0,   ( pdims_partner(1) / pdims(1) ) - 1
                    map_coord(1) = i+offset(1)
                    map_coord(2) = j+offset(2)
    
                    target_idex = f_rnk_lst(map_coord(1),map_coord(2)) + numprocs
    
                    CALL MPI_RECV( bdims_rem,   6, MPI_INTEGER, target_idex, 10, &
                        comm_inter,status, ierr )
    
                    bdims (1,1) = bdims_rem (1,1) / cfratio(1)
                    bdims (1,2) = bdims_rem (1,2) / cfratio(1)
                    bdims (2,1) = bdims_rem (2,1) / cfratio(2)
                    bdims (2,2) = bdims_rem (2,2) / cfratio(2)
                    bdims (3,1) = bdims_rem (3,2) / cfratio(3)
                    bdims (3,2) = bdims     (3,1) + 2
    
                    CALL MPI_SEND( bdims,       6, MPI_INTEGER, target_idex, 9, &
                        comm_inter, ierr )
    
    
                    n_cell_c =      ( (bdims(1,2)-bdims(1,1)) + 3 ) * &
                        ( (bdims(2,2)-bdims(2,1)) + 3 ) * &
                        ( (bdims(3,2)-bdims(3,1)) + 1 )
    
                    CALL MPI_SEND(kh(bdims(3,1)  :bdims(3,2)  , &
                        bdims(2,1)-1:bdims(2,2)+1, &
                        bdims(1,1)-1:bdims(1,2)+1),&
                        n_cell_c, MPI_REAL, target_idex,    &
                        207, comm_inter, ierr)
    
                    CALL MPI_SEND(km(bdims(3,1)  :bdims(3,2)  , &
                        bdims(2,1)-1:bdims(2,2)+1, &
                        bdims(1,1)-1:bdims(1,2)+1),&
                        n_cell_c, MPI_REAL, target_idex,    &
                        208, comm_inter, ierr)
    
    
    
                end do
            end do
    
        ELSEIF ( coupling_mode == 'vnested_fine' )  THEN
                      ! Receive data from coarse grid for TOP BC
    
            offset(1) = pcoord(1) / ( pdims(1)/pdims_partner(1) )
            offset(2) = pcoord(2) / ( pdims(2)/pdims_partner(2) )
            map_coord(1) = offset(1)
            map_coord(2) = offset(2)
            target_idex = c_rnk_lst(map_coord(1),map_coord(2))
    
            bdims (1,1) = nxl
            bdims (1,2) = nxr
            bdims (2,1) = nys
            bdims (2,2) = nyn
            bdims (3,1) = nzb
            bdims (3,2) = nzt
    
            CALL MPI_SEND( bdims,       6, MPI_INTEGER, target_idex, 10, &
                comm_inter, ierr )
    
            CALL MPI_RECV( bdims_rem,   6, MPI_INTEGER, target_idex, 9, &
                comm_inter,status, ierr )
    
            n_cell_c =       ( (bdims_rem(1,2)-bdims_rem(1,1)) + 3 ) * &
                ( (bdims_rem(2,2)-bdims_rem(2,1)) + 3 ) * &
                ( (bdims_rem(3,2)-bdims_rem(3,1)) + 1 )
    
            ALLOCATE( work3d  ( bdims_rem(3,1)  :bdims_rem(3,2)  , &
                bdims_rem(2,1)-1:bdims_rem(2,2)+1, &
                bdims_rem(1,1)-1:bdims_rem(1,2)+1))
    
    
            CALL MPI_RECV( work3d,n_cell_c, MPI_REAL, target_idex, 207, &
                comm_inter,status, ierr )
    
        ! Neumann BC for FG kh
        DO jjf = nys, nyn
            DO iif = nxl, nxr
               kh(nzt+1,jjf,iif) = kh(nzt,jjf,iif)
            END DO
        END DO
    
        CALL MPI_RECV( work3d,n_cell_c, MPI_REAL, target_idex, 208, &
             comm_inter,status, ierr )
    
        ! Neumann BC for FG kh
        DO jjf = nys, nyn
            DO iif = nxl, nxr
               km(nzt+1,jjf,iif) = km(nzt,jjf,iif)
            END DO
        END DO
    
    
        !
        !-- The following evaluation can only be performed, if the fine grid is situated below the inversion
        !!    DO jjf = nys-1, nyn+1
        !!       DO iif = nxl-1, nxr+1
        !!
        !!          km(nzt+1,jjf,iif) = 0.1 * l_grid(nzt+1) * SQRT( e(nzt+1,jjf,iif) )
        !!          kh(nzt+1,jjf,iif) = 3.0 * km(nzt+1,jjf,iif)
        !!
        !!       END DO
        !!    END DO
    
        CALL exchange_horiz_2d(km(nzt+1,:,:) )
        CALL exchange_horiz_2d(kh(nzt+1,:,:) )
    
        DEALLOCATE (   work3d      )
    
        ENDIF
    
    
!    CONTAINS
! 
!       SUBROUTINE vnest_set_topbc_kh
!       
!       
!           USE arrays_3d
!           USE control_parameters
!           USE grid_variables
!           USE indices
!           USE pegrid
!           
!       
!           IMPLICIT NONE
!
!           INTEGER(iwp)                       :: i
!           INTEGER(iwp)                       :: j
!           INTEGER(iwp)                       :: k
!           INTEGER(iwp)                       :: iif
!           INTEGER(iwp)                       :: jjf
!           INTEGER(iwp)                       :: bottomx
!           INTEGER(iwp)                       :: bottomy
!           INTEGER(iwp)                       :: topx
!           INTEGER(iwp)                       :: topy
!           REAL(wp)                           :: eps
!           REAL(wp)                           :: alpha
!           REAL(wp)                           :: eminus
!           REAL(wp)                           :: edot
!           REAL(wp)                           :: eplus
!           REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: ptprf
!           REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: ptprs
!       
!           
!       
!           ALLOCATE( ptprf(bdims_rem(3,1):bdims_rem(3,2),bdims_rem(2,1)-1:bdims_rem(2,2)+1,nxl:nxr) )
!           ALLOCATE( ptprs(bdims_rem(3,1):bdims_rem(3,2),nys:nyn,nxl:nxr) )
!       
!           !
!           !-- Determination of a boundary condition for the potential temperature pt:
!           !-- The scheme derived by Clark and Farley can be used in all three dimensions.
!       
!           !
!           !-- Interpolation in x-direction
!       
!           DO k = bdims_rem(3,1), bdims_rem(3,2)
!       
!               DO j = bdims_rem(2,1)-1, bdims_rem(2,2)+1
!       
!                   DO i = bdims_rem(1,1), bdims_rem(1,2)
!       
!                       bottomx = (nxf+1)/(nxc+1) * i
!                       topx    = (nxf+1)/(nxc+1) *(i+1) - 1
!       
!                       DO iif = bottomx, topx
!       
!                           eps = (iif * dxf + 0.5 * dxf - i * dxc - 0.5 * dxc) / dxc
!       
!                           alpha = ( (dxf/dxc)**2.0 - 1.0) / 24.0
!       
!                           eminus = eps * (eps - 1.0 ) / 2.0 + alpha
!       
!                           edot = ( 1.0 - eps**2.0 ) - 2.0 * alpha
!       
!                           eplus = eps * ( eps + 1.0 ) / 2.0 + alpha
!       
!                           ptprf(k,j,iif) = eminus * work3d(k,j,i-1) &
!                               + edot  * work3d(k,j,i)   &
!                               + eplus  * work3d(k,j,i+1)
!                       END DO
!       
!                   END DO
!       
!               END DO
!       
!           END DO
!       
!           !
!           !-- Interpolation in y-direction
!       
!           DO k = bdims_rem(3,1), bdims_rem(3,2)
!       
!               DO j = bdims_rem(2,1), bdims_rem(2,2)
!       
!                   bottomy = (nyf+1)/(nyc+1) * j
!                   topy    = (nyf+1)/(nyc+1) * (j+1) - 1
!       
!                   DO iif = nxl, nxr
!       
!                       DO jjf = bottomy, topy
!       
!                           eps = (jjf * dyf + 0.5 * dyf - j * dyc - 0.5 * dyc) / dyc
!       
!                           alpha = ( (dyf/dyc)**2.0 - 1.0) / 24.0
!                   
!                           eminus = eps * (eps - 1.0) / 2.0 + alpha
!       
!                           edot  = ( 1.0 - eps**2.0 ) - 2.0 * alpha
!       
!                           eplus  = eps * ( eps + 1.0 ) / 2.0 + alpha
!       
!                           ptprs(k,jjf,iif) = eminus * ptprf(k,j-1,iif) &
!                               + edot  * ptprf(k,j,iif)   &
!                               + eplus  * ptprf(k,j+1,iif)
!                       END DO
!       
!                   END DO
!       
!               END DO
!       
!           END DO
!       
!           !
!           !-- Interpolation in z-direction
!       
!           DO jjf = nys, nyn
!               DO iif = nxl, nxr
!       
!                   eps = ( zuf(nzt+1) - zuc(bdims_rem(3,1)+1) ) / dzc
!       
!                   alpha = ( (dzf/dzc)**2.0 - 1.0) / 24.0
!       
!                   eminus = eps * ( eps - 1.0 ) / 2.0 + alpha
!       
!                   edot  = ( 1.0 - eps**2.0 ) - 2.0 * alpha
!       
!                   eplus  = eps * ( eps + 1.0 ) / 2.0 + alpha
!       
!                   kh (nzt+1,jjf,iif) = eminus * ptprs(bdims_rem(3,1),jjf,iif)   &
!                       + edot  * ptprs(bdims_rem(3,1)+1,jjf,iif) &
!                       + eplus  * ptprs(bdims_rem(3,1)+2,jjf,iif)
!       
!               END DO
!           END DO          
!           
!           DEALLOCATE ( ptprf, ptprs )
!       
!       
!       
!       END SUBROUTINE vnest_set_topbc_kh
       
!       SUBROUTINE vnest_set_topbc_km
!       
!       
!           USE arrays_3d
!           USE control_parameters
!           USE grid_variables
!           USE indices
!           USE pegrid
!           
!       
!           IMPLICIT NONE
!
!           INTEGER(iwp)                       :: i
!           INTEGER(iwp)                       :: j
!           INTEGER(iwp)                       :: k
!           INTEGER(iwp)                       :: iif
!           INTEGER(iwp)                       :: jjf
!           INTEGER(iwp)                       :: bottomx
!           INTEGER(iwp)                       :: bottomy
!           INTEGER(iwp)                       :: topx
!           INTEGER(iwp)                       :: topy
!           REAL(wp)                           :: eps
!           REAL(wp)                           :: alpha
!           REAL(wp)                           :: eminus
!           REAL(wp)                           :: edot
!           REAL(wp)                           :: eplus
!           REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: ptprf
!           REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: ptprs
!       
!           
!       
!           ALLOCATE( ptprf(bdims_rem(3,1):bdims_rem(3,2),bdims_rem(2,1)-1:bdims_rem(2,2)+1,nxl:nxr) )
!           ALLOCATE( ptprs(bdims_rem(3,1):bdims_rem(3,2),nys:nyn,nxl:nxr) )
!       
!           !
!           !-- Determination of a boundary condition for the potential temperature pt:
!           !-- The scheme derived by Clark and Farley can be used in all three dimensions.
!       
!           !
!           !-- Interpolation in x-direction
!       
!           DO k = bdims_rem(3,1), bdims_rem(3,2)
!       
!               DO j = bdims_rem(2,1)-1, bdims_rem(2,2)+1
!       
!                   DO i = bdims_rem(1,1), bdims_rem(1,2)
!       
!                       bottomx = (nxf+1)/(nxc+1) * i
!                       topx    = (nxf+1)/(nxc+1) *(i+1) - 1
!       
!                       DO iif = bottomx, topx
!       
!                           eps = (iif * dxf + 0.5 * dxf - i * dxc - 0.5 * dxc) / dxc
!       
!                           alpha = ( (dxf/dxc)**2.0 - 1.0) / 24.0
!       
!                           eminus = eps * (eps - 1.0 ) / 2.0 + alpha
!       
!                           edot = ( 1.0 - eps**2.0 ) - 2.0 * alpha
!       
!                           eplus = eps * ( eps + 1.0 ) / 2.0 + alpha
!       
!                           ptprf(k,j,iif) = eminus * work3d(k,j,i-1) &
!                               + edot  * work3d(k,j,i)   &
!                               + eplus  * work3d(k,j,i+1)
!                       END DO
!       
!                   END DO
!       
!               END DO
!       
!           END DO
!       
!           !
!           !-- Interpolation in y-direction
!       
!           DO k = bdims_rem(3,1), bdims_rem(3,2)
!       
!               DO j = bdims_rem(2,1), bdims_rem(2,2)
!       
!                   bottomy = (nyf+1)/(nyc+1) * j
!                   topy    = (nyf+1)/(nyc+1) * (j+1) - 1
!       
!                   DO iif = nxl, nxr
!       
!                       DO jjf = bottomy, topy
!       
!                           eps = (jjf * dyf + 0.5 * dyf - j * dyc - 0.5 * dyc) / dyc
!       
!                           alpha = ( (dyf/dyc)**2.0 - 1.0) / 24.0
!                   
!                           eminus = eps * (eps - 1.0) / 2.0 + alpha
!       
!                           edot  = ( 1.0 - eps**2.0 ) - 2.0 * alpha
!       
!                           eplus  = eps * ( eps + 1.0 ) / 2.0 + alpha
!       
!                           ptprs(k,jjf,iif) = eminus * ptprf(k,j-1,iif) &
!                               + edot  * ptprf(k,j,iif)   &
!                               + eplus  * ptprf(k,j+1,iif)
!                       END DO
!       
!                   END DO
!       
!               END DO
!       
!           END DO
!       
!           !
!           !-- Interpolation in z-direction
!       
!           DO jjf = nys, nyn
!               DO iif = nxl, nxr
!       
!                   eps = ( zuf(nzt+1) - zuc(bdims_rem(3,1)+1) ) / dzc
!       
!                   alpha = ( (dzf/dzc)**2.0 - 1.0) / 24.0
!       
!                   eminus = eps * ( eps - 1.0 ) / 2.0 + alpha
!       
!                   edot  = ( 1.0 - eps**2.0 ) - 2.0 * alpha
!       
!                   eplus  = eps * ( eps + 1.0 ) / 2.0 + alpha
!       
!                   km (nzt+1,jjf,iif) = eminus * ptprs(bdims_rem(3,1),jjf,iif)   &
!                       + edot  * ptprs(bdims_rem(3,1)+1,jjf,iif) &
!                       + eplus  * ptprs(bdims_rem(3,1)+2,jjf,iif)
!       
!               END DO
!           END DO          
!           
!           DEALLOCATE ( ptprf, ptprs )
!       
!       
!       
!       END SUBROUTINE vnest_set_topbc_km
    
 
#endif
    END SUBROUTINE vnest_boundary_conds_khkm
   
   
   
    SUBROUTINE vnest_anterpolate

#if defined( __parallel )
   
        !--------------------------------------------------------------------------------!
        ! Description:
        ! ------------
        ! Anterpolate data from fine grid to coarse grid.
        !------------------------------------------------------------------------------!
   
        USE arrays_3d
        USE control_parameters
        USE grid_variables
        USE indices
        USE interfaces
        USE pegrid
        USE surface_mod,                                                           &
            ONLY :  bc_h
        
   
        IMPLICIT NONE
   
        REAL(wp)                              ::  time_since_reference_point_rem
        INTEGER(iwp)                          ::  i
        INTEGER(iwp)                          ::  j
        INTEGER(iwp)                          ::  k 
        INTEGER(iwp)                          ::  l  !< running index boundary type, for up- and downward-facing walls
        INTEGER(iwp)                          ::  m  !< running index surface elements
   
   
   
        !
        !-- In case of model termination initiated by the remote model
        !-- (terminate_coupled_remote > 0), initiate termination of the local model.
        !-- The rest of the coupler must then be skipped because it would cause an MPI
        !-- intercomminucation hang.
        !-- If necessary, the coupler will be called at the beginning of the next
        !-- restart run.
   
        IF ( myid == 0) THEN
            CALL MPI_SENDRECV( terminate_coupled,        1, MPI_INTEGER, &
                target_id, 0,                             &
                terminate_coupled_remote, 1, MPI_INTEGER, &
                target_id, 0,                             &
                comm_inter, status, ierr )
        ENDIF
        CALL MPI_BCAST( terminate_coupled_remote, 1, MPI_INTEGER, 0, comm2d, &
            ierr )
   
        IF ( terminate_coupled_remote > 0 )  THEN
            WRITE( message_string, * ) 'remote model "',                       &
                TRIM( coupling_mode_remote ),                                  &
                '" terminated',                                                &
                '&with terminate_coupled_remote = ',                           &
                terminate_coupled_remote,                                      &
                '&local model  "', TRIM( coupling_mode ),                      &
                '" has',                                                       &
                '&terminate_coupled = ',                                       &
                terminate_coupled
            CALL message( 'vnest_anterpolate', 'PA0310', 1, 2, 0, 6, 0 )
            RETURN
        ENDIF
   
   
        !
        !-- Exchange the current simulated time between the models
   
        IF ( myid == 0 ) THEN
   
            CALL MPI_SEND( time_since_reference_point, 1, MPI_REAL, target_id, &
                11, comm_inter, ierr )
            CALL MPI_RECV( time_since_reference_point_rem, 1, MPI_REAL,        &
                target_id, 11, comm_inter, status, ierr )
   
        ENDIF
   
        CALL MPI_BCAST( time_since_reference_point_rem, 1, MPI_REAL, 0, comm2d, &
            ierr )
   
        IF ( coupling_mode == 'vnested_crse' )  THEN
                      ! Receive data from fine grid for anterpolation
   
            offset(1) = ( pdims_partner(1) / pdims(1) ) * pcoord(1)
            offset(2) = ( pdims_partner(2) / pdims(2) ) * pcoord(2)
   
            do j = 0,   ( pdims_partner(2) / pdims(2) ) - 1
                do i = 0,   ( pdims_partner(1) / pdims(1) ) - 1
                    map_coord(1) = i+offset(1)
                    map_coord(2) = j+offset(2)
   
                    target_idex = f_rnk_lst(map_coord(1),map_coord(2)) + numprocs
   
                    CALL MPI_RECV( bdims_rem,   6, MPI_INTEGER, target_idex, 10, &
                        comm_inter,status, ierr )
   
                    bdims (1,1) = bdims_rem (1,1) / cfratio(1)
                    bdims (1,2) = bdims_rem (1,2) / cfratio(1)
                    bdims (2,1) = bdims_rem (2,1) / cfratio(2)
                    bdims (2,2) = bdims_rem (2,2) / cfratio(2)
                    bdims (3,1) = bdims_rem (3,1)
                    bdims (3,2) = bdims_rem (3,2) / cfratio(3)
   
                    CALL MPI_SEND( bdims, 6, MPI_INTEGER, target_idex, 9, &
                        comm_inter, ierr )
   
                    n_cell_c =                      &
                        (bdims(1,2)-bdims(1,1)+1) * &
                        (bdims(2,2)-bdims(2,1)+1) * &
                        (bdims(3,2)-bdims(3,1)+0)
   
                    CALL MPI_RECV( u(            &
                        bdims(3,1)+1:bdims(3,2), &
                        bdims(2,1)  :bdims(2,2), &
                        bdims(1,1)  :bdims(1,2)),&
                        n_cell_c, MPI_REAL, target_idex, 101, &
                        comm_inter,status, ierr )
   
                    CALL MPI_RECV( v(            &
                        bdims(3,1)+1:bdims(3,2), &
                        bdims(2,1)  :bdims(2,2), &
                        bdims(1,1)  :bdims(1,2)),&
                        n_cell_c, MPI_REAL, target_idex, 102, &
                        comm_inter,status, ierr )
   
                    CALL MPI_RECV(pt(            &
                        bdims(3,1)+1:bdims(3,2), &
                        bdims(2,1)  :bdims(2,2), &
                        bdims(1,1)  :bdims(1,2)),&
                        n_cell_c, MPI_REAL, target_idex, 105, &
                        comm_inter,status, ierr )
   
                 IF ( humidity ) THEN
                    CALL MPI_RECV(q(            &
                        bdims(3,1)+1:bdims(3,2), &
                        bdims(2,1)  :bdims(2,2), &
                        bdims(1,1)  :bdims(1,2)),&
                        n_cell_c, MPI_REAL, target_idex, 106, &
                        comm_inter,status, ierr )
                 ENDIF
   
                    CALL MPI_RECV( w(              &
                        bdims(3,1)  :bdims(3,2)-1, &
                        bdims(2,1)  :bdims(2,2),   &
                        bdims(1,1)  :bdims(1,2)),  &
                        n_cell_c, MPI_REAL, target_idex, 103, &
                        comm_inter,status, ierr )
   
                end do
            end do
   
   
   
            !
            !-- Boundary conditions for the velocity components u and v
   
   
            IF ( ibc_uv_b == 0 ) THEN
                u(nzb,:,:) = 0.0_wp 
                v(nzb,:,:) = 0.0_wp 
            ELSE
                u(nzb,:,:) = u(nzb+1,:,:)
                v(nzb,:,:) = v(nzb+1,:,:)
            END IF
            !
            !-- Boundary conditions for the velocity components w

            w(nzb,:,:) = 0.0_wp
   
!
!-- Temperature at bottom boundary.
!-- Neumann, zero-gradient
            IF ( ibc_pt_b == 1 )  THEN 
               DO  l = 0, 1 
                  DO  m = 1, bc_h(l)%ns
                     i = bc_h(l)%i(m)     
                     j = bc_h(l)%j(m)
                     k = bc_h(l)%k(m)
                     pt(k+bc_h(l)%koff,j,i) = pt(k,j,i)
                  ENDDO
               ENDDO
            ENDIF
   
   
            CALL exchange_horiz( u, nbgp )
            CALL exchange_horiz( v, nbgp )
            CALL exchange_horiz( w, nbgp )
            CALL exchange_horiz( pt, nbgp )
   
   
        ELSEIF ( coupling_mode == 'vnested_fine' )  THEN
                      ! Send data to coarse grid for anterpolation
   
            offset(1) = pcoord(1) / ( pdims(1)/pdims_partner(1) )
            offset(2) = pcoord(2) / ( pdims(2)/pdims_partner(2) )
            map_coord(1) = offset(1)
            map_coord(2) = offset(2)
            target_idex = c_rnk_lst(map_coord(1),map_coord(2))
   
!-- Limit anterpolation level to nzt - z nesting ratio (a pseudo-buffer layer)
            bdims (1,1) = nxl
            bdims (1,2) = nxr
            bdims (2,1) = nys
            bdims (2,2) = nyn
            bdims (3,1) = nzb
            bdims (3,2) = nzt-cfratio(3) 
   
            CALL MPI_SEND( bdims, 6, MPI_INTEGER, target_idex, 10, &
                comm_inter, ierr )
   
            CALL MPI_RECV( bdims_rem,   6, MPI_INTEGER, target_idex, 9, &
                comm_inter,status, ierr )
   
   
            ALLOCATE( work3d (                   &
                bdims_rem(3,1)+1:bdims_rem(3,2), &
                bdims_rem(2,1)  :bdims_rem(2,2), &
                bdims_rem(1,1)  :bdims_rem(1,2)))
   
   
            anterpol3d => u
      
            CALL anterpolate_to_crse_u
            CALL MPI_SEND( work3d, 1, TYPE_VNEST_ANTER, target_idex,  &
                101, comm_inter, ierr)
   
            anterpol3d => v
   
            CALL anterpolate_to_crse_v
            CALL MPI_SEND( work3d, 1, TYPE_VNEST_ANTER, target_idex,  &
                102, comm_inter, ierr)
   
            anterpol3d => pt
   
            CALL anterpolate_to_crse_s
            CALL MPI_SEND( work3d, 1, TYPE_VNEST_ANTER, target_idex,  &
                105, comm_inter, ierr)
   
   
          IF ( humidity ) THEN
   
            anterpol3d => q
   
            CALL anterpolate_to_crse_s
            CALL MPI_SEND( work3d, 1, TYPE_VNEST_ANTER, target_idex,  &
                106, comm_inter, ierr)
          ENDIF
   
   
            DEALLOCATE(   work3d                    )
            ALLOCATE(   work3d ( bdims_rem(3,1)  :bdims_rem(3,2)-1, &
                                 bdims_rem(2,1)  :bdims_rem(2,2), &
                                 bdims_rem(1,1)  :bdims_rem(1,2)))
            anterpol3d => w
            CALL anterpolate_to_crse_w
            CALL MPI_SEND( work3d, 1, TYPE_VNEST_ANTER, target_idex,  &
                103, comm_inter, ierr)
   
            NULLIFY   (   anterpol3d                )
            DEALLOCATE(   work3d                    )
   
        ENDIF
   
   

    CONTAINS
       SUBROUTINE anterpolate_to_crse_u
      
      
           USE arrays_3d
           USE control_parameters
           USE grid_variables
           USE indices
           USE pegrid
           
      
           IMPLICIT NONE

           INTEGER(iwp)                       :: i
           INTEGER(iwp)                       :: j
           INTEGER(iwp)                       :: k
           INTEGER(iwp)                       :: iif
           INTEGER(iwp)                       :: jjf
           INTEGER(iwp)                       :: kkf
           INTEGER(iwp)                       :: bottomy
           INTEGER(iwp)                       :: bottomz
           INTEGER(iwp)                       :: topy
           INTEGER(iwp)                       :: topz
           REAL(wp)                           :: aweight
      
           !
           !-- Anterpolation of the velocity components u
           !-- only values in yz-planes that coincide in the fine and
           !-- the coarse grid are considered
      
           DO k = bdims_rem(3,1)+1, bdims_rem(3,2)
      
               bottomz = (dzc/dzf) * (k-1) + 1
               topz    = (dzc/dzf) * k
      
               DO j = bdims_rem(2,1),bdims_rem(2,2)
       
                   bottomy = (nyf+1) / (nyc+1) * j
                   topy    = (nyf+1) / (nyc+1) * (j+1) - 1
      
                   DO i = bdims_rem(1,1),bdims_rem(1,2)
      
                       iif = (nxf+1) / (nxc+1) * i
      
                       aweight   = 0.0
      
                       DO kkf = bottomz, topz
                           DO jjf = bottomy, topy
      
                               aweight   = aweight + anterpol3d(kkf,jjf,iif) *        &
                                   (dzf/dzc) * (dyf/dyc)
      
                           END DO
                       END DO
      
                       work3d(k,j,i)   = aweight
      
                   END DO
      
               END DO
      
           END DO
      
      
      
       END SUBROUTINE anterpolate_to_crse_u
      
      
       SUBROUTINE anterpolate_to_crse_v
      
      
           USE arrays_3d
           USE control_parameters
           USE grid_variables
           USE indices
           USE pegrid
           
      
           IMPLICIT NONE

           INTEGER(iwp)                       :: i
           INTEGER(iwp)                       :: j
           INTEGER(iwp)                       :: k
           INTEGER(iwp)                       :: iif
           INTEGER(iwp)                       :: jjf
           INTEGER(iwp)                       :: kkf
           INTEGER(iwp)                       :: bottomx
           INTEGER(iwp)                       :: bottomz
           INTEGER(iwp)                       :: topx
           INTEGER(iwp)                       :: topz
           REAL(wp)                           :: aweight

           !
           !-- Anterpolation of the velocity components v
           !-- only values in xz-planes that coincide in the fine and
           !-- the coarse grid are considered
      
           DO k = bdims_rem(3,1)+1, bdims_rem(3,2)
      
               bottomz = (dzc/dzf) * (k-1) + 1
               topz    = (dzc/dzf) * k
      
               DO j = bdims_rem(2,1), bdims_rem(2,2)
      
                   jjf = (nyf+1) / (nyc+1) * j
      
                   DO i = bdims_rem(1,1), bdims_rem(1,2)
      
                       bottomx = (nxf+1) / (nxc+1) * i
                       topx    = (nxf+1) / (nxc+1) * (i+1) - 1
      
                       aweight   = 0.0
      
                       DO kkf = bottomz, topz
                           DO iif = bottomx, topx
      
                               aweight   = aweight + anterpol3d(kkf,jjf,iif) *        &
                                   (dzf/dzc) * (dxf/dxc)
      
      
                           END DO
                       END DO
      
                       work3d(k,j,i)   = aweight
      
                   END DO
               END DO
           END DO
      
      
      
       END SUBROUTINE anterpolate_to_crse_v
      
      
       SUBROUTINE anterpolate_to_crse_w
      
      
           USE arrays_3d
           USE control_parameters
           USE grid_variables
           USE indices
           USE pegrid
           
      
           IMPLICIT NONE

           INTEGER(iwp)                       :: i
           INTEGER(iwp)                       :: j
           INTEGER(iwp)                       :: k
           INTEGER(iwp)                       :: iif
           INTEGER(iwp)                       :: jjf
           INTEGER(iwp)                       :: kkf
           INTEGER(iwp)                       :: bottomx
           INTEGER(iwp)                       :: bottomy
           INTEGER(iwp)                       :: topx
           INTEGER(iwp)                       :: topy
           REAL(wp)                           :: aweight

           !
           !-- Anterpolation of the velocity components w
           !-- only values in xy-planes that coincide in the fine and
           !-- the coarse grid are considered
      
           DO k = bdims_rem(3,1), bdims_rem(3,2)-1
      
               kkf = cfratio(3) * k
      
               DO j = bdims_rem(2,1), bdims_rem(2,2)
      
                   bottomy = (nyf+1) / (nyc+1) * j
                   topy    = (nyf+1) / (nyc+1) * (j+1) - 1
      
                   DO i = bdims_rem(1,1), bdims_rem(1,2)
      
                       bottomx = (nxf+1) / (nxc+1) * i
                       topx    = (nxf+1) / (nxc+1) * (i+1) - 1
      
                       aweight   = 0.0
      
                       DO jjf = bottomy, topy
                           DO iif = bottomx, topx
      
                               aweight   = aweight + anterpol3d (kkf,jjf,iif) *        &
                                   (dxf/dxc) * (dyf/dyc)
      
                           END DO
                       END DO
      
                       work3d(k,j,i)   = aweight
      
                   END DO
      
               END DO
      
           END DO
      
      
       END SUBROUTINE anterpolate_to_crse_w
      
      
       SUBROUTINE anterpolate_to_crse_s
      
      
           USE arrays_3d
           USE control_parameters
           USE grid_variables
           USE indices
           USE pegrid
           
      
           IMPLICIT NONE

           INTEGER(iwp)                       :: i
           INTEGER(iwp)                       :: j
           INTEGER(iwp)                       :: k
           INTEGER(iwp)                       :: iif
           INTEGER(iwp)                       :: jjf
           INTEGER(iwp)                       :: kkf
           INTEGER(iwp)                       :: bottomx
           INTEGER(iwp)                       :: bottomy
           INTEGER(iwp)                       :: bottomz
           INTEGER(iwp)                       :: topx
           INTEGER(iwp)                       :: topy
           INTEGER(iwp)                       :: topz
           REAL(wp)                           :: aweight
      
           !
           !-- Anterpolation of the potential temperature pt
           !-- all fine grid values are considered
      
           DO k = bdims_rem(3,1)+1, bdims_rem(3,2)
      
               bottomz = (dzc/dzf) * (k-1) + 1
               topz    = (dzc/dzf) * k
      
               DO j = bdims_rem(2,1), bdims_rem(2,2)
      
                   bottomy = (nyf+1) / (nyc+1) * j
                   topy    = (nyf+1) / (nyc+1) * (j+1) - 1
      
                   DO i = bdims_rem(1,1), bdims_rem(1,2)
      
                       bottomx = (nxf+1) / (nxc+1) * i
                       topx    = (nxf+1) / (nxc+1) * (i+1) - 1
      
                       aweight   = 0.0
          
                       DO kkf = bottomz, topz
                           DO jjf = bottomy, topy
                               DO iif = bottomx, topx
      
                                   aweight   = aweight  + anterpol3d(kkf,jjf,iif) *                  &
                                       (dzf/dzc) * (dyf/dyc) * (dxf/dxc)
      
                               END DO
                           END DO
                       END DO
      
                       work3d(k,j,i)   = aweight
      
                   END DO
      
               END DO
      
           END DO
      
      
       END SUBROUTINE anterpolate_to_crse_s
#endif       
    END SUBROUTINE vnest_anterpolate
   
   
   
    SUBROUTINE vnest_anterpolate_e
#if defined( __parallel )
   
        !--------------------------------------------------------------------------------!
        ! Description:
        ! ------------
        ! Anterpolate TKE from fine grid to coarse grid.
        !------------------------------------------------------------------------------!
   
        USE arrays_3d
        USE control_parameters
        USE grid_variables
        USE indices
        USE interfaces
        USE pegrid
        
   
        IMPLICIT NONE
   
        REAL(wp)                              :: time_since_reference_point_rem
        INTEGER(iwp)                          :: i
        INTEGER(iwp)                          :: j
   
        !
        !-- In case of model termination initiated by the remote model
        !-- (terminate_coupled_remote > 0), initiate termination of the local model.
        !-- The rest of the coupler must then be skipped because it would cause an MPI
        !-- intercomminucation hang.
        !-- If necessary, the coupler will be called at the beginning of the next
        !-- restart run.
   
        IF ( myid == 0) THEN
            CALL MPI_SENDRECV( terminate_coupled,        1, MPI_INTEGER,       &
                target_id, 0,                                                  &
                terminate_coupled_remote, 1, MPI_INTEGER,                      &
                target_id, 0,                                                  &
                comm_inter, status, ierr )
        ENDIF
        CALL MPI_BCAST( terminate_coupled_remote, 1, MPI_INTEGER, 0, comm2d,   &
            ierr )
   
        IF ( terminate_coupled_remote > 0 )  THEN
            WRITE( message_string, * ) 'remote model "',                       &
                TRIM( coupling_mode_remote ),                                  &
                '" terminated',                                                &
                '&with terminate_coupled_remote = ',                           &
                terminate_coupled_remote,                                      &
                '&local model  "', TRIM( coupling_mode ),                      & 
                '" has',                                                       &
                '&terminate_coupled = ',                                       &
                terminate_coupled
            CALL message( 'vnest_anterpolate_e', 'PA0310', 1, 2, 0, 6, 0 )
            RETURN
        ENDIF
   
   
        !
        !-- Exchange the current simulated time between the models
        IF ( myid == 0 ) THEN
   
            CALL MPI_SEND( time_since_reference_point, 1, MPI_REAL, target_id, &
                11, comm_inter, ierr )
            CALL MPI_RECV( time_since_reference_point_rem, 1, MPI_REAL,        &
                target_id, 11, comm_inter, status, ierr )
   
        ENDIF
   
        CALL MPI_BCAST( time_since_reference_point_rem, 1, MPI_REAL, 0, comm2d, &
            ierr )
   
        IF ( coupling_mode == 'vnested_crse' )  THEN
                      ! Receive data from fine grid for anterpolation
   
            offset(1) = ( pdims_partner(1) / pdims(1) ) * pcoord(1)
            offset(2) = ( pdims_partner(2) / pdims(2) ) * pcoord(2)
   
            do j = 0,   ( pdims_partner(2) / pdims(2) ) - 1
                do i = 0,   ( pdims_partner(1) / pdims(1) ) - 1
                    map_coord(1) = i+offset(1)
                    map_coord(2) = j+offset(2)
   
                    target_idex = f_rnk_lst(map_coord(1),map_coord(2)) + numprocs
  
                    bdims (1,1) = f2c_dims_cg (0,map_coord(1),map_coord(2)) 
                    bdims (1,2) = f2c_dims_cg (1,map_coord(1),map_coord(2)) 
                    bdims (2,1) = f2c_dims_cg (2,map_coord(1),map_coord(2)) 
                    bdims (2,2) = f2c_dims_cg (3,map_coord(1),map_coord(2)) 
                    bdims (3,1) = f2c_dims_cg (4,map_coord(1),map_coord(2))
                    bdims (3,2) = f2c_dims_cg (5,map_coord(1),map_coord(2)) 


                     n_cell_c = (bdims(1,2)-bdims(1,1)+1) * &
                         (bdims(2,2)-bdims(2,1)+1) * &
                         (bdims(3,2)-bdims(3,1)+0)


                    CALL MPI_RECV( e( bdims(3,1)+1:bdims(3,2), &
                        bdims(2,1)  :bdims(2,2), &
                        bdims(1,1)  :bdims(1,2)),&
                        n_cell_c, MPI_REAL, target_idex, 104, &
                        comm_inter,status, ierr )
                end do
            end do
   
   
            !
            !-- Boundary conditions
   
            IF ( .NOT. constant_diffusion ) THEN
                e(nzb,:,:)   = e(nzb+1,:,:)
            END IF
   
            IF ( .NOT. constant_diffusion )  CALL exchange_horiz( e, nbgp )
   
        ELSEIF ( coupling_mode == 'vnested_fine' )  THEN
                      ! Send data to coarse grid for anterpolation
   
            offset(1) = pcoord(1) / ( pdims(1)/pdims_partner(1) )
            offset(2) = pcoord(2) / ( pdims(2)/pdims_partner(2) )
            map_coord(1) = offset(1)
            map_coord(2) = offset(2)
            target_idex = c_rnk_lst(map_coord(1),map_coord(2))

            bdims_rem (1,1) =  f2c_dims_fg (0)
            bdims_rem (1,2) =  f2c_dims_fg (1)
            bdims_rem (2,1) =  f2c_dims_fg (2)
            bdims_rem (2,2) =  f2c_dims_fg (3)
            bdims_rem (3,1) =  f2c_dims_fg (4)
            bdims_rem (3,2) =  f2c_dims_fg (5)

            ALLOCATE( work3d (                   &
                bdims_rem(3,1)+1:bdims_rem(3,2), &
                bdims_rem(2,1)  :bdims_rem(2,2), &
                bdims_rem(1,1)  :bdims_rem(1,2)))
   
            anterpol3d => e
   
            CALL anterpolate_to_crse_e

            CALL MPI_SEND( work3d, 1, TYPE_VNEST_ANTER, target_idex,  &
                           104, comm_inter, ierr)

            NULLIFY   (   anterpol3d                )
            DEALLOCATE(   work3d                    )
        ENDIF
   
   
    CONTAINS 
   
   
   
   
   
       SUBROUTINE anterpolate_to_crse_e
      
      
           USE arrays_3d
           USE control_parameters
           USE grid_variables
           USE indices
           USE pegrid
           
      
           IMPLICIT NONE

           INTEGER(iwp)                       :: i
           INTEGER(iwp)                       :: j
           INTEGER(iwp)                       :: k
           INTEGER(iwp)                       :: iif
           INTEGER(iwp)                       :: jjf
           INTEGER(iwp)                       :: kkf
           INTEGER(iwp)                       :: bottomx
           INTEGER(iwp)                       :: bottomy
           INTEGER(iwp)                       :: bottomz
           INTEGER(iwp)                       :: topx
           INTEGER(iwp)                       :: topy
           INTEGER(iwp)                       :: topz
           REAL(wp)                           :: aweight_a
           REAL(wp)                           :: aweight_b
           REAL(wp)                           :: aweight_c
           REAL(wp)                           :: aweight_d
           REAL(wp)                           :: aweight_e
           REAL(wp)                           :: energ 
          
           DO k = bdims_rem(3,1)+1, bdims_rem(3,2)
      
               bottomz = (dzc/dzf) * (k-1) + 1
               topz    = (dzc/dzf) * k
             
               DO j = bdims_rem(2,1), bdims_rem(2,2)
      
                   bottomy = (nyf+1) / (nyc+1) * j
                   topy    = (nyf+1) / (nyc+1) * (j+1) - 1
      
                   DO i = bdims_rem(1,1), bdims_rem(1,2)
      
                       bottomx = (nxf+1) / (nxc+1) * i
                       topx    = (nxf+1) / (nxc+1) * (i+1) - 1
      
                       aweight_a   = 0.0
                       aweight_b   = 0.0
                       aweight_c   = 0.0
                       aweight_d   = 0.0
                       aweight_e   = 0.0
      
                       DO kkf = bottomz, topz
                           DO jjf = bottomy, topy
                               DO iif = bottomx, topx
      
                                   aweight_a = aweight_a + anterpol3d(kkf,jjf,iif)  *                     &
                                       (dzf/dzc) * (dyf/dyc) * (dxf/dxc)
                              
      
                                   energ = ( 0.5 * ( u(kkf,jjf,iif)   + u(kkf,jjf,iif+1) ) )**2.0 +          &
                                       ( 0.5 * ( v(kkf,jjf,iif)   + v(kkf,jjf+1,iif) ) )**2.0 +              &
                                       ( 0.5 * ( w(kkf-1,jjf,iif) + w(kkf,jjf,iif) ) )**2.0
      
                                   aweight_b   =   aweight_b + energ         *                         &
                                       (dzf/dzc) * (dyf/dyc) * (dxf/dxc)
       
                                   aweight_c   =   aweight_c + 0.5 * ( u(kkf,jjf,iif) + u(kkf,jjf,iif+1) ) * &
                                       (dzf/dzc) * (dyf/dyc) * (dxf/dxc)
      
                                   aweight_d   =   aweight_d + 0.5 * ( v(kkf,jjf,iif) + v(kkf,jjf+1,iif) ) * &
                                       (dzf/dzc) * (dyf/dyc) * (dxf/dxc)
       
                                   aweight_e   =   aweight_e + 0.5 * ( w(kkf-1,jjf,iif) + w(kkf,jjf,iif) ) * &
                                       (dzf/dzc) * (dyf/dyc) * (dxf/dxc)
      
      
                               END DO
                           END DO
                       END DO
      
                       work3d(k,j,i)  = aweight_a + 0.5 * ( aweight_b -  &
                           aweight_c**2.0 -                              &
                           aweight_d**2.0 -                              &
                           aweight_e**2.0                 )
      
                   END DO
      
               END DO
      
           END DO
      
       
      
       END SUBROUTINE anterpolate_to_crse_e
#endif       
    END SUBROUTINE vnest_anterpolate_e

    SUBROUTINE vnest_init_pegrid_rank
#if defined( __parallel )
! Domain decomposition and exchange of grid variables between coarse and fine
! Given processor coordinates as index f_rnk_lst(pcoord(1), pcoord(2))
!   returns the rank. A single coarse block will have to send data to multiple
!   fine blocks. In the coarse grid the pcoords of the remote block is first found and then using
!   f_rnk_lst the target_idex is identified.
! blk_dim stores the index limits of a given block. blk_dim_remote is received
! from the asscoiated nest partner.
! cf_ratio(1:3) is the ratio between fine and coarse grid: nxc/nxf, nyc/nyf and
! ceiling(dxc/dxf) 


       USE control_parameters,                                                 &
           ONLY:  coupling_mode

       USE kinds
         
       USE pegrid
 

       IMPLICIT NONE

       INTEGER(iwp)                           :: dest_rnk
       INTEGER(iwp)                           ::  i                        !<

       IF (myid == 0) THEN
           IF ( coupling_mode == 'vnested_crse') THEN
               CALL MPI_SEND( pdims, 2, MPI_INTEGER, numprocs, 33, comm_inter, &
                   ierr )
               CALL MPI_RECV( pdims_partner, 2, MPI_INTEGER, numprocs, 66,     &
                   comm_inter, status, ierr )
           ELSEIF ( coupling_mode == 'vnested_fine') THEN
               CALL MPI_RECV( pdims_partner, 2, MPI_INTEGER, 0, 33,      &
                   comm_inter, status, ierr )
               CALL MPI_SEND( pdims, 2, MPI_INTEGER, 0, 66, comm_inter, &
                   ierr )
           ENDIF
       ENDIF


       IF ( coupling_mode == 'vnested_crse') THEN
           CALL MPI_BCAST( pdims_partner, 2, MPI_INTEGER, 0, comm2d, ierr )
           ALLOCATE( c_rnk_lst( 0:(pdims(1)-1)    ,0:(pdims(2)-1) )   )
           ALLOCATE( f_rnk_lst( 0:(pdims_partner(1)-1)  ,0:(pdims_partner(2)-1) ) )
           do i=0,numprocs-1
               CALL MPI_CART_COORDS( comm2d, i, ndim, pcoord, ierr )
               call MPI_Cart_rank(comm2d, pcoord, dest_rnk, ierr)
               c_rnk_lst(pcoord(1),pcoord(2)) = dest_rnk
           end do
       ELSEIF ( coupling_mode == 'vnested_fine') THEN
           CALL MPI_BCAST( pdims_partner, 2, MPI_INTEGER, 0, comm2d, ierr )
           ALLOCATE( c_rnk_lst( 0:(pdims_partner(1)-1)  ,0:(pdims_partner(2)-1) ) )
           ALLOCATE( f_rnk_lst( 0:(pdims(1)-1)    ,0:(pdims(2)-1) )   )

           do i=0,numprocs-1
               CALL MPI_CART_COORDS( comm2d, i, ndim, pcoord, ierr )
               call MPI_Cart_rank(comm2d, pcoord, dest_rnk, ierr)
               f_rnk_lst(pcoord(1),pcoord(2)) = dest_rnk
           enddo
       ENDIF


       IF ( coupling_mode == 'vnested_crse') THEN
           if (myid == 0) then
               CALL MPI_SEND( c_rnk_lst, pdims(1)*pdims(2),     MPI_INTEGER, numprocs, 0, comm_inter, ierr )
               CALL MPI_RECV( f_rnk_lst, pdims_partner(1)*pdims_partner(2), MPI_INTEGER, numprocs, 4, comm_inter,status, ierr )
           end if
           CALL MPI_BCAST( f_rnk_lst, pdims_partner(1)*pdims_partner(2), MPI_INTEGER, 0, comm2d, ierr )
       ELSEIF ( coupling_mode == 'vnested_fine') THEN
           if (myid == 0) then
               CALL MPI_RECV( c_rnk_lst, pdims_partner(1)*pdims_partner(2),  MPI_INTEGER, 0, 0, comm_inter,status, ierr )
               CALL MPI_SEND( f_rnk_lst, pdims(1)*pdims(2),      MPI_INTEGER, 0, 4, comm_inter, ierr )
           end if
           CALL MPI_BCAST( c_rnk_lst, pdims_partner(1)*pdims_partner(2), MPI_INTEGER, 0, comm2d, ierr )
       ENDIF

       !-- Reason for MPI error unknown; solved if three lines duplicated
       CALL MPI_CART_COORDS( comm2d, myid, ndim, pcoord, ierr )
       CALL MPI_CART_SHIFT( comm2d, 0, 1, pleft, pright, ierr )
       CALL MPI_CART_SHIFT( comm2d, 1, 1, psouth, pnorth, ierr )


#endif
 
    END SUBROUTINE vnest_init_pegrid_rank


    SUBROUTINE vnest_init_pegrid_domain
#if defined( __parallel )

       USE control_parameters,                                                 &
           ONLY:  coupling_mode, coupling_topology, dz,                        &
                  dz_stretch_level_start, message_string

       USE grid_variables,                                                     &
           ONLY:  dx, dy
           
       USE indices,                                                            &
           ONLY: nbgp, nx, ny, nz, nxl, nxr, nys, nyn, nzb, nzt

       USE kinds
         
       USE pegrid
       
       IMPLICIT NONE

       INTEGER(iwp)                           :: i              !<
       INTEGER(iwp)                           :: j              !<
       INTEGER(iwp)                           :: tempx
       INTEGER(iwp)                           :: tempy
       INTEGER(iwp)                           :: TYPE_INT_YZ
       INTEGER(iwp)                           :: SIZEOFREAL
       INTEGER(iwp)                           :: MTV_X
       INTEGER(iwp)                           :: MTV_Y
       INTEGER(iwp)                           :: MTV_Z
       INTEGER(iwp)                           :: MTV_RX
       INTEGER(iwp)                           :: MTV_RY
       INTEGER(iwp)                           :: MTV_RZ
   
       !
       !--    Pass the number of grid points of the coarse model to
       !--    the nested model and vice versa
       IF ( coupling_mode ==  'vnested_crse' )  THEN

           nxc = nx
           nyc = ny
           nzc = nz
           dxc = dx
           dyc = dy
           dzc = dz(1)
           cg_nprocs = numprocs

           IF ( myid == 0 )  THEN

               CALL MPI_SEND( nxc, 1, MPI_INTEGER  , numprocs, 1,  comm_inter, &
                   ierr )
               CALL MPI_SEND( nyc, 1, MPI_INTEGER  , numprocs, 2,  comm_inter, &
                   ierr )
               CALL MPI_SEND( nzc, 1, MPI_INTEGER  , numprocs, 3,  comm_inter, &
                   ierr )
               CALL MPI_SEND( dxc, 1, MPI_REAL     , numprocs, 4,  comm_inter, &
                   ierr )
               CALL MPI_SEND( dyc, 1, MPI_REAL     , numprocs, 5,  comm_inter, &
                   ierr )
               CALL MPI_SEND( dzc, 1, MPI_REAL     , numprocs, 6,  comm_inter, &
                   ierr )
               CALL MPI_SEND( pdims, 2, MPI_INTEGER, numprocs, 7,  comm_inter, &
                   ierr )
               CALL MPI_SEND( cg_nprocs, 1, MPI_INTEGER, numprocs, 8,  comm_inter,  &
                   ierr )
               CALL MPI_RECV( nxf, 1, MPI_INTEGER,   numprocs, 21, comm_inter, &
                   status, ierr )
               CALL MPI_RECV( nyf, 1, MPI_INTEGER,   numprocs, 22, comm_inter, &
                   status, ierr )
               CALL MPI_RECV( nzf, 1, MPI_INTEGER,   numprocs, 23, comm_inter, &
                   status, ierr )
               CALL MPI_RECV( dxf, 1, MPI_REAL,      numprocs, 24, comm_inter, &
                   status, ierr )
               CALL MPI_RECV( dyf, 1, MPI_REAL,      numprocs, 25, comm_inter, &
                   status, ierr )
               CALL MPI_RECV( dzf, 1, MPI_REAL,      numprocs, 26, comm_inter, &
                   status, ierr )
               CALL MPI_RECV( pdims_partner, 2, MPI_INTEGER,                   &
                   numprocs, 27, comm_inter, status, ierr )
               CALL MPI_RECV( fg_nprocs, 1, MPI_INTEGER,                       &
                   numprocs, 28, comm_inter, status, ierr )
           ENDIF

           CALL MPI_BCAST( nxf, 1,   MPI_INTEGER, 0, comm2d, ierr )
           CALL MPI_BCAST( nyf, 1,   MPI_INTEGER, 0, comm2d, ierr )
           CALL MPI_BCAST( nzf, 1,   MPI_INTEGER, 0, comm2d, ierr )
           CALL MPI_BCAST( dxf, 1,   MPI_REAL,    0, comm2d, ierr )
           CALL MPI_BCAST( dyf, 1,   MPI_REAL,    0, comm2d, ierr )
           CALL MPI_BCAST( dzf, 1,   MPI_REAL,    0, comm2d, ierr )
           CALL MPI_BCAST( pdims_partner, 2, MPI_INTEGER, 0, comm2d, ierr )
           CALL MPI_BCAST( fg_nprocs,  1, MPI_INTEGER, 0, comm2d, ierr )
           
!
!--        Check if stretching is used within the nested domain. ABS(...) is 
!--        necessary because of the default value of -9999999.9_wp (negative)
           IF ( ABS( dz_stretch_level_start(1) ) <= (nzf+1)*dzf )  THEN       
               message_string = 'Stretching in the parent domain is '//        &
                                'only allowed above the nested domain'
               CALL message( 'vnest_init_pegrid_domain', 'PA0497', 1, 2, 0, 6, 0 )
           ENDIF

       ELSEIF ( coupling_mode ==  'vnested_fine' )  THEN

           nxf = nx
           nyf = ny
           nzf = nz
           dxf = dx
           dyf = dy
           dzf = dz(1)
           fg_nprocs = numprocs

           IF ( myid == 0 ) THEN

               CALL MPI_RECV( nxc, 1, MPI_INTEGER, 0, 1, comm_inter, status, &
                   ierr )
               CALL MPI_RECV( nyc, 1, MPI_INTEGER, 0, 2, comm_inter, status, &
                   ierr )
               CALL MPI_RECV( nzc, 1, MPI_INTEGER, 0, 3, comm_inter, status, &
                   ierr )
               CALL MPI_RECV( dxc, 1, MPI_REAL,    0, 4, comm_inter, status, &
                   ierr )
               CALL MPI_RECV( dyc, 1, MPI_REAL,    0, 5, comm_inter, status, &
                   ierr )
               CALL MPI_RECV( dzc, 1, MPI_REAL,    0, 6, comm_inter, status, &
                   ierr )
               CALL MPI_RECV( pdims_partner, 2, MPI_INTEGER, 0, 7, comm_inter, &
                   status, ierr )
               CALL MPI_RECV( cg_nprocs,  1, MPI_INTEGER, 0, 8, comm_inter, &
                   status, ierr )
               CALL MPI_SEND( nxf, 1, MPI_INTEGER, 0, 21,  comm_inter, ierr )
               CALL MPI_SEND( nyf, 1, MPI_INTEGER, 0, 22,  comm_inter, ierr )
               CALL MPI_SEND( nzf, 1, MPI_INTEGER, 0, 23, comm_inter, ierr )
               CALL MPI_SEND( dxf, 1, MPI_REAL,    0, 24, comm_inter, ierr )
               CALL MPI_SEND( dyf, 1, MPI_REAL,    0, 25, comm_inter, ierr )
               CALL MPI_SEND( dzf, 1, MPI_REAL,    0, 26, comm_inter, ierr )
               CALL MPI_SEND( pdims,2,MPI_INTEGER, 0, 27, comm_inter, ierr )
               CALL MPI_SEND( fg_nprocs,1,MPI_INTEGER, 0, 28, comm_inter, ierr )
           ENDIF

           CALL MPI_BCAST( nxc, 1, MPI_INTEGER, 0, comm2d, ierr)
           CALL MPI_BCAST( nyc, 1, MPI_INTEGER, 0, comm2d, ierr)
           CALL MPI_BCAST( nzc, 1, MPI_INTEGER, 0, comm2d, ierr)
           CALL MPI_BCAST( dxc, 1, MPI_REAL,    0, comm2d, ierr)
           CALL MPI_BCAST( dyc, 1, MPI_REAL,    0, comm2d, ierr)
           CALL MPI_BCAST( dzc, 1, MPI_REAL,    0, comm2d, ierr)
           CALL MPI_BCAST( pdims_partner, 2, MPI_INTEGER, 0, comm2d, ierr)
           CALL MPI_BCAST( cg_nprocs,  1, MPI_INTEGER, 0, comm2d, ierr)

       ENDIF
       
       ngp_c = ( nxc+1 + 2 * nbgp ) * ( nyc+1 + 2 * nbgp )
       ngp_f = ( nxf+1 + 2 * nbgp ) * ( nyf+1 + 2 * nbgp )

       IF ( coupling_mode(1:8) == 'vnested_')  coupling_topology = 1


!-- Nesting Ratio: For each coarse grid cell how many fine grid cells exist
       cfratio(1) = INT ( (nxf+1) / (nxc+1) )
       cfratio(2) = INT ( (nyf+1) / (nyc+1) )
       cfratio(3) = CEILING (  dzc    /  dzf  )

!-- target_id is used only for exhange of information like simulated_time
!-- which are then MPI_BCAST to other processors in the group
          IF ( myid == 0 )  THEN
   
              IF ( TRIM( coupling_mode ) == 'vnested_crse' )  THEN
                  target_id = numprocs
              ELSE IF ( TRIM( coupling_mode ) == 'vnested_fine' )  THEN
                  target_id = 0
              ENDIF
   
          ENDIF

!-- Store partner grid dimenstions and create MPI derived types

         IF ( coupling_mode == 'vnested_crse' )  THEN

            offset(1) = ( pdims_partner(1) / pdims(1) ) * pcoord(1)
            offset(2) = ( pdims_partner(2) / pdims(2) ) * pcoord(2)
      
            tempx =  ( pdims_partner(1) / pdims(1) ) - 1
            tempy =  ( pdims_partner(2) / pdims(2) ) - 1
            ALLOCATE( c2f_dims_cg (0:5,offset(1):tempx+offset(1),offset(2):tempy+offset(2) ) )
            ALLOCATE( f2c_dims_cg (0:5,offset(1):tempx+offset(1),offset(2):tempy+offset(2) ) )

            do j = 0,   ( pdims_partner(2) / pdims(2) ) - 1
                do i = 0,   ( pdims_partner(1) / pdims(1) ) - 1
                   map_coord(1) = i+offset(1)
                   map_coord(2) = j+offset(2)
      
                   target_idex = f_rnk_lst(map_coord(1),map_coord(2)) + numprocs

                   CALL MPI_RECV( bdims_rem,   6, MPI_INTEGER, target_idex, 10, &
                       comm_inter,status, ierr )

!-- Store the CG dimensions that correspond to the FG partner; needed for FG top BC
!-- One CG can have multiple FG partners. The 3D array is mapped by partner proc co-ord
                   c2f_dims_cg (0,map_coord(1),map_coord(2)) = bdims_rem (1,1) / cfratio(1)
                   c2f_dims_cg (1,map_coord(1),map_coord(2)) = bdims_rem (1,2) / cfratio(1)
                   c2f_dims_cg (2,map_coord(1),map_coord(2)) = bdims_rem (2,1) / cfratio(2)
                   c2f_dims_cg (3,map_coord(1),map_coord(2)) = bdims_rem (2,2) / cfratio(2)
                   c2f_dims_cg (4,map_coord(1),map_coord(2)) = bdims_rem (3,2) / cfratio(3)
                   c2f_dims_cg (5,map_coord(1),map_coord(2)) =(bdims_rem (3,2) / cfratio(3)) + 2
 
!-- Store the CG dimensions that  correspond to the FG partner; needed for anterpolation
                   f2c_dims_cg (0,map_coord(1),map_coord(2)) = bdims_rem (1,1) / cfratio(1)
                   f2c_dims_cg (1,map_coord(1),map_coord(2)) = bdims_rem (1,2) / cfratio(1)
                   f2c_dims_cg (2,map_coord(1),map_coord(2)) = bdims_rem (2,1) / cfratio(2)
                   f2c_dims_cg (3,map_coord(1),map_coord(2)) = bdims_rem (2,2) / cfratio(2)
                   f2c_dims_cg (4,map_coord(1),map_coord(2)) = bdims_rem (3,1)
                   f2c_dims_cg (5,map_coord(1),map_coord(2)) =(bdims_rem (3,2)-cfratio(3))/ cfratio(3)
   
                   CALL MPI_SEND( c2f_dims_cg (:,map_coord(1),map_coord(2)), 6, &
                      MPI_INTEGER, target_idex, 100, comm_inter, ierr )
      
                   CALL MPI_SEND( f2c_dims_cg (:,map_coord(1),map_coord(2)), 6, &
                      MPI_INTEGER, target_idex, 101, comm_inter, ierr )
      
                end do
            end do
     
!-- A derived data type to pack 3 Z-levels of CG to set FG top BC
            MTV_X = ( nxr - nxl + 1 ) + 2*nbgp
            MTV_Y = ( nyn - nys + 1 ) + 2*nbgp
            MTV_Z = nzt+1 - nzb +1
            
            MTV_RX = ( c2f_dims_cg (1,offset(1),offset(2)) - c2f_dims_cg (0,offset(1),offset(2)) ) +1+2
            MTV_RY = ( c2f_dims_cg (3,offset(1),offset(2)) - c2f_dims_cg (2,offset(1),offset(2)) ) +1+2
            MTV_RZ = ( c2f_dims_cg (5,offset(1),offset(2)) - c2f_dims_cg (4,offset(1),offset(2)) ) +1
            
            CALL MPI_TYPE_EXTENT(MPI_REAL, SIZEOFREAL, IERR)

            CALL MPI_TYPE_VECTOR ( MTV_RY, MTV_RZ, MTV_Z, MPI_REAL, TYPE_INT_YZ, IERR)
            CALL MPI_TYPE_HVECTOR( MTV_RX, 1, MTV_Z*MTV_Y*SIZEOFREAL, &
                                              TYPE_INT_YZ, TYPE_VNEST_BC, IERR)
            CALL MPI_TYPE_FREE(TYPE_INT_YZ, IERR)
            CALL MPI_TYPE_COMMIT(TYPE_VNEST_BC, IERR)   


         ELSEIF ( coupling_mode == 'vnested_fine' )  THEN

           ALLOCATE( c2f_dims_fg (0:5) )
           ALLOCATE( f2c_dims_fg (0:5) )
    
           offset(1) = pcoord(1) / ( pdims(1)/pdims_partner(1) )
           offset(2) = pcoord(2) / ( pdims(2)/pdims_partner(2) )
           map_coord(1) = offset(1)
           map_coord(2) = offset(2)
           target_idex = c_rnk_lst(map_coord(1),map_coord(2))
    
           bdims (1,1) = nxl
           bdims (1,2) = nxr
           bdims (2,1) = nys
           bdims (2,2) = nyn
           bdims (3,1) = nzb
           bdims (3,2) = nzt

           CALL MPI_SEND( bdims,    6, MPI_INTEGER, target_idex, 10, &
               comm_inter, ierr )
           
!-- Store the CG dimensions that correspond to the FG partner; needed for FG top BC
!-- One FG can have only one CG partner
           CALL MPI_RECV( c2f_dims_fg,   6, MPI_INTEGER, target_idex, 100, &
               comm_inter,status, ierr )
           
           CALL MPI_RECV( f2c_dims_fg,   6, MPI_INTEGER, target_idex, 101, &
               comm_inter,status, ierr )
 
!-- Store the CG dimensions that  correspond to the FG partner; needed for anterpolation

           n_cell_c = (f2c_dims_fg(1)-f2c_dims_fg(0)+1) * &
                      (f2c_dims_fg(3)-f2c_dims_fg(2)+1) * &
                      (f2c_dims_fg(5)-f2c_dims_fg(4)+0)

           CALL MPI_TYPE_CONTIGUOUS(n_cell_c, MPI_REAL, TYPE_VNEST_ANTER, IERR)
           CALL MPI_TYPE_COMMIT(TYPE_VNEST_ANTER, ierr) 

        ENDIF
#endif   
       END SUBROUTINE vnest_init_pegrid_domain


       SUBROUTINE vnest_init_grid
 
#if defined( __parallel )
          USE arrays_3d,                                                       &
              ONLY:  zu, zw
              
          USE control_parameters,                                              &
              ONLY:  coupling_mode, message_string, number_stretch_level_start
             
          USE indices,                                                         &
              ONLY:  nzt
         
          USE kinds
          
          USE pegrid

          IMPLICIT NONE
      
!
!--       Allocate and Exchange zuc and zuf, zwc and zwf
          IF ( coupling_mode(1:8)  == 'vnested_' )  THEN
          
             ALLOCATE( zuc(0:nzc+1), zuf(0:nzf+1) )
             ALLOCATE( zwc(0:nzc+1), zwf(0:nzf+1) )
          
             IF ( coupling_mode ==  'vnested_crse' )  THEN
                
                zuc = zu
                zwc = zw
                
                IF ( myid == 0 )  THEN
          
                   CALL MPI_SEND( zuc, nzt+2, MPI_REAL, numprocs, 41, comm_inter,  &
                                  ierr )
                   CALL MPI_RECV( zuf, nzf+2, MPI_REAL, numprocs, 42, comm_inter,  &
                                  status, ierr )
          
                   CALL MPI_SEND( zwc, nzt+2, MPI_REAL, numprocs, 43, comm_inter,  &
                                  ierr )
                   CALL MPI_RECV( zwf, nzf+2, MPI_REAL, numprocs, 44, comm_inter,  &
                                  status, ierr )
          
                ENDIF
          
                CALL MPI_BCAST( zuf,nzf+2,MPI_REAL,    0, comm2d, ierr ) 
                CALL MPI_BCAST( zwf,nzf+2,MPI_REAL,    0, comm2d, ierr ) 
          
             ELSEIF ( coupling_mode ==  'vnested_fine' )  THEN
          
!
!--             Check if stretching is used within the nested domain
                IF ( number_stretch_level_start > 0 )  THEN
                   message_string = 'Stretching in the nested domain is not '//&
                           'allowed'
                   CALL message( 'vnest_init_grid', 'PA0498', 1, 2, 0, 6, 0 )
                ENDIF
                
                zuf = zu
                zwf = zw
                
                IF ( myid == 0 )  THEN
          
                   CALL MPI_RECV( zuc,nzc+2, MPI_REAL, 0, 41, comm_inter, status, &
                                  ierr )
                   CALL MPI_SEND( zuf,nzt+2, MPI_REAL, 0, 42, comm_inter, ierr )
          
                   CALL MPI_RECV( zwc,nzc+2, MPI_REAL, 0, 43, comm_inter, status, &
                                  ierr )
                   CALL MPI_SEND( zwf,nzt+2, MPI_REAL, 0, 44, comm_inter, ierr )
                ENDIF
          
                CALL MPI_BCAST( zuc,nzc+2,MPI_REAL,  0, comm2d, ierr ) 
                CALL MPI_BCAST( zwc,nzc+2,MPI_REAL,  0, comm2d, ierr ) 
          
             ENDIF
          ENDIF

#endif
       END SUBROUTINE vnest_init_grid


       SUBROUTINE vnest_check_parameters
#if defined( __parallel )

          USE pegrid,                                                                &
              ONLY:  myid

          IMPLICIT NONE

          IF (myid==0) PRINT*, '*** vnest: check parameters not implemented yet ***'
      
#endif
       END SUBROUTINE vnest_check_parameters


       SUBROUTINE vnest_timestep_sync

#if defined( __parallel )
         USE control_parameters,                                                    &
             ONLY:  coupling_mode, dt_3d
      
         USE interfaces
      
         USE kinds

         USE pegrid

         IMPLICIT NONE

         IF ( coupling_mode == 'vnested_crse')  THEN
            dtc = dt_3d
               if (myid == 0) then
                   CALL MPI_SEND( dt_3d, 1, MPI_REAL, target_id,                     &
                        31, comm_inter, ierr )
                   CALL MPI_RECV( dtf, 1, MPI_REAL,                                  &
                        target_id, 32, comm_inter, status, ierr )

                endif
                CALL MPI_BCAST( dtf, 1,   MPI_REAL,    0, comm2d, ierr ) 
         ELSE
                dtf = dt_3d
                if (myid == 0) then
                   CALL MPI_RECV( dtc, 1, MPI_REAL,                                  &
                        target_id, 31, comm_inter, status, ierr )
                   CALL MPI_SEND( dt_3d, 1, MPI_REAL, target_id,                     &
                        32, comm_inter, ierr )

                endif
                CALL MPI_BCAST( dtc, 1,   MPI_REAL,    0, comm2d, ierr ) 
         
         ENDIF
!-- Identical timestep for coarse and fine grids
          dt_3d = MIN( dtc, dtf )
#endif
       END SUBROUTINE vnest_timestep_sync
       
       SUBROUTINE vnest_deallocate
#if defined( __parallel )
          USE control_parameters,                                                    &
              ONLY:  coupling_mode
             
          IMPLICIT NONE
          
          IF ( ALLOCATED(c_rnk_lst) ) DEALLOCATE (c_rnk_lst)
          IF ( ALLOCATED(f_rnk_lst) ) DEALLOCATE (f_rnk_lst)
          
          IF ( coupling_mode == 'vnested_crse')  THEN
             IF ( ALLOCATED (c2f_dims_cg) ) DEALLOCATE (c2f_dims_cg)
             IF ( ALLOCATED (f2c_dims_cg) ) DEALLOCATE (f2c_dims_cg)
          ELSEIF( coupling_mode == 'vnested_fine' )  THEN
             IF ( ALLOCATED (c2f_dims_fg) ) DEALLOCATE (c2f_dims_fg)
             IF ( ALLOCATED (f2c_dims_fg) ) DEALLOCATE (f2c_dims_fg)
          ENDIF
#endif 
       END SUBROUTINE vnest_deallocate

#endif
 END MODULE vertical_nesting_mod
