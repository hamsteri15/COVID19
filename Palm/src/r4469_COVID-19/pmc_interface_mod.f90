!> @file pmc_interface_mod.f90
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
! $Id: pmc_interface_mod.f90 4444 2020-03-05 15:59:50Z raasch $
! bugfix: cpp-directives and variable declarations for serial mode added
! 
! 4413 2020-02-19 15:52:19Z hellstea
! All the USE-statements within subroutines moved up to the module declaration section.
! 
! 4385 2020-01-27 08:37:37Z hellstea
! Error messages PA0425 and PA0426 made more specific
! 
! 4360 2020-01-07 11:25:50Z suehring
! Introduction of wall_flags_total_0, which currently sets bits based on static
! topography information used in wall_flags_static_0
! 
! 4329 2019-12-10 15:46:36Z motisi
! Renamed wall_flags_0 to wall_flags_static_0
! 
! 4273 2019-10-24 13:40:54Z monakurppa
! Add a logical switch nesting_chem and rename nest_salsa to nesting_salsa
! 
! 4260 2019-10-09 14:04:03Z hellstea
! Rest of the possibly round-off-error sensitive grid-line matching tests
! changed to round-off-error tolerant forms throughout the module.
! 
! 4249 2019-10-01 12:27:47Z hellstea
! Several grid-line matching tests changed to a round-off-error tolerant form
! in pmci_setup_parent, pmci_define_index_mapping and pmci_check_grid_matching.
! 
! 4182 2019-08-22 15:20:23Z scharf
! Corrected "Former revisions" section
! 
! 4168 2019-08-16 13:50:17Z suehring
! Replace function get_topography_top_index by topo_top_ind
! 
! 4029 2019-06-14 14:04:35Z raasch
! nest_chemistry switch removed
! 
! 4026 2019-06-12 16:50:15Z suehring
! Masked topography at boundary grid points in mass conservation, in order to 
! avoid that mean velocities within topography are imposed
! 
! 4011 2019-05-31 14:34:03Z hellstea
! Mass (volume) flux correction included to ensure global mass conservation for child domains.
! 
! 3987 2019-05-22 09:52:13Z kanani
! Introduce alternative switch for debug output during timestepping
! 
! 3984 2019-05-16 15:17:03Z hellstea
! Commenting improved, pmci_map_fine_to_coarse_grid renamed as pmci_map_child_grid_to_parent_grid,
! set_child_edge_coords renamed as pmci_set_child_edge_coords, some variables renamed, etc. 
! 
! 3979 2019-05-15 13:54:29Z hellstea
! Bugfix in pmc_interp_1sto_sn. This bug had effect only in case of 1-d domain
! decomposition with npex = 1. 
! 
! 3976 2019-05-15 11:02:34Z hellstea
! Child initialization also for the redundant ghost points behind the nested 
! boundaries added (2nd and 3rd ghost-point layers and corners).
! 
! 3948 2019-05-03 14:49:57Z hellstea
! Some variables renamed, a little cleaning up and some commenting improvements 
! 
! 3947 2019-05-03 07:56:44Z hellstea
! The checks included in 3946 are extended for the z-direction and moved into its
! own subroutine called from pmci_define_index_mapping.
! 
! 3946 2019-05-02 14:18:59Z hellstea
! Check added for child domains too small in terms of number of parent-grid cells so
! that anterpolation is not possible. Checks added for too wide anterpolation buffer
! for the same reason. Some minor code reformatting done.
!
! 3945 2019-05-02 11:29:27Z raasch
!
! 3932 2019-04-24 17:31:34Z suehring
! Add missing if statements for call of pmc_set_dataarray_name for TKE and
! dissipation.
!
! 3888 2019-04-12 09:18:10Z hellstea
! Variables renamed, commenting improved etc.
! 
! 3885 2019-04-11 11:29:34Z kanani
! Changes related to global restructuring of location messages and introduction 
! of additional debug messages
! 
! 3883 2019-04-10 12:51:50Z hellstea
! Checks and error messages improved and extended. All the child index bounds in the
! parent-grid index space are made module variables. Function get_number_of_childs
! renamed get_number_of_children. A number of variables renamed 
! and qite a lot of other code reshaping made all around the module. 
! 
! 3876 2019-04-08 18:41:49Z knoop
! Implemented nesting for salsa variables.
! 
! 3833 2019-03-28 15:04:04Z forkel
! replaced USE chem_modules by USE chem_gasphase_mod 
! 
! 3822 2019-03-27 13:10:23Z hellstea
! Temporary increase of the vertical dimension of the parent-grid arrays and 
! workarrc_t is cancelled as unnecessary.
! 
! 3819 2019-03-27 11:01:36Z hellstea
! Adjustable anterpolation buffer introduced on all nest boundaries, it is controlled
! by the new nesting_parameters parameter anterpolation_buffer_width.
! 
! 3804 2019-03-19 13:46:20Z hellstea
! Anterpolation domain is lowered from kct-1 to kct-3 to avoid exessive       
! kinetic energy from building up in CBL flows.
! 
! 3803 2019-03-19 13:44:40Z hellstea
! A bug fixed in lateral boundary interpolations. Dimension of val changed from  
! 5 to 3 in pmci_setup_parent and pmci_setup_child.
! 
! 3794 2019-03-15 09:36:33Z raasch
! two remaining unused variables removed
! 
! 3792 2019-03-14 16:50:07Z hellstea
! Interpolations improved. Large number of obsolete subroutines removed.
! All unused variables removed. 
! 
! 3741 2019-02-13 16:24:49Z hellstea
! Interpolations and child initialization adjusted to handle set ups with child 
! pe-subdomain dimension not integer divisible by the grid-spacing ratio in the 
! respective direction. Set ups with pe-subdomain dimension smaller than the 
! grid-spacing ratio in the respective direction are now forbidden.
! 
! 3708 2019-01-30 12:58:13Z hellstea
! Checks for parent / child grid line matching introduced.
! Interpolation of nest-boundary-tangential velocity components revised.
! 
! 3697 2019-01-24 17:16:13Z hellstea
! Bugfix: upper k-bound in the child initialization interpolation 
! pmci_interp_1sto_all corrected.
! Copying of the nest boundary values into the redundant 2nd and 3rd ghost-node 
! layers is added to the pmci_interp_1sto_*-routines.
! 
! 3681 2019-01-18 15:06:05Z hellstea
! Linear interpolations are replaced by first order interpolations. The linear 
! interpolation routines are still included but not called. In the child 
! inititialization the interpolation is also changed to 1st order and the linear
! interpolation is not kept.
! Subroutine pmci_map_fine_to_coarse_grid is rewritten.
! Several changes in pmci_init_anterp_tophat.
! Child's parent-grid arrays (uc, vc,...) are made non-overlapping on the PE-
! subdomain boundaries in order to allow grid-spacing ratios higher than nbgp. 
! Subroutine pmci_init_tkefactor is removed as unnecessary.
! 
! 3655 2019-01-07 16:51:22Z knoop
! Remove unused variable simulated_time
! 
! 1762 2016-02-25 12:31:13Z hellstea
! Initial revision by A. Hellsten
!
! Description:
! ------------
! Domain nesting interface routines. The low-level inter-domain communication   
! is conducted by the PMC-library routines.
!
! @todo Remove array_3d variables from USE statements thate not used in the 
!       routine
! @todo Data transfer of qc and nc is prepared but not activated
!------------------------------------------------------------------------------!
 MODULE pmc_interface

#if ! defined( __parallel )
!
!-- Serial mode does not allow nesting, but requires the following variables as steering
!-- quantities
    USE kinds

    IMPLICIT NONE

    PUBLIC

    CHARACTER(LEN=8), SAVE ::  nesting_mode = 'none'   !< steering parameter for 1- or 2-way nesting

    INTEGER(iwp), SAVE     ::  comm_world_nesting    !< Global nesting communicator
    INTEGER(iwp), SAVE     ::  cpl_id  = 1           !<

    LOGICAL, SAVE ::  nested_run = .FALSE.        !< general switch
    LOGICAL, SAVE ::  rans_mode_parent = .FALSE.  !< parent model mode (.F.-LES mode, .T.-RANS mode)

#else

    USE ISO_C_BINDING


    USE arrays_3d,                                                             &
        ONLY:  diss, diss_2, dzu, dzw, e, e_p, e_2, nc, nc_2, nc_p, nr, nr_2,  &
               pt, pt_2, q, q_2, qc, qc_2, qr, qr_2, s, s_2,                   &
               u, u_p, u_2, v, v_p, v_2, w, w_p, w_2, zu, zw
    
    USE chem_gasphase_mod,                                                     &
        ONLY:  nspec

    USE chem_modules,                                                          &
        ONLY:  chem_species, ibc_cs_b, nesting_chem

    USE chemistry_model_mod,                                                   &
        ONLY:  spec_conc_2
    
    USE control_parameters,                                                    &
        ONLY:  air_chemistry, bc_dirichlet_l, bc_dirichlet_n, bc_dirichlet_r,  &
               bc_dirichlet_s, child_domain,                                   &
               constant_diffusion, constant_flux_layer,                        &
               coupling_char, end_time,                                        &
               debug_output_timestep,                                          &
               dt_restart, dt_3d, dz, humidity,                                &
               ibc_pt_b, ibc_q_b, ibc_s_b, ibc_uv_b,                           &
               message_string, neutral, passive_scalar, rans_mode, rans_tke_e, &
               restart_time,                                                   &
               roughness_length, salsa, topography, volume_flow, time_restart
    
    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point_s

    USE grid_variables,                                                        &
        ONLY:  dx, dy

    USE indices,                                                               &
        ONLY:  nbgp, nx, nxl, nxlg, nxlu, nxr, nxrg, ny, nyn, nyng, nys, nysg, &
               nysv, nz, nzb, nzt, topo_top_ind, wall_flags_total_0

    USE bulk_cloud_model_mod,                                                  &
        ONLY: bulk_cloud_model, microphysics_morrison, microphysics_seifert

    USE particle_attributes,                                                   &
        ONLY:  particle_advection

    USE kinds

#if defined( __parallel )
#if !defined( __mpifh )
    USE MPI
#endif

    USE pegrid,                                                                &
        ONLY:  collective_wait, comm1dx, comm1dy, comm2d, myid, myidx, myidy,  &
               numprocs, pdims, pleft, pnorth, pright, psouth, status

    USE pmc_child,                                                             &
        ONLY:  pmc_childinit, pmc_c_clear_next_array_list,                     &
               pmc_c_getnextarray, pmc_c_get_2d_index_list, pmc_c_getbuffer,   &
               pmc_c_putbuffer, pmc_c_setind_and_allocmem,                     &
               pmc_c_set_dataarray, pmc_set_dataarray_name

    USE pmc_general,                                                           &
        ONLY:  da_namelen, pmc_max_array

    USE pmc_handle_communicator,                                               &
        ONLY:  pmc_get_model_info, pmc_init_model, pmc_is_rootmodel,           &
               pmc_no_namelist_found, pmc_parent_for_child, m_couplers

    USE pmc_mpi_wrapper,                                                       &
        ONLY:  pmc_bcast, pmc_recv_from_child, pmc_recv_from_parent,           &
               pmc_send_to_child, pmc_send_to_parent

    USE pmc_parent,                                                            &
        ONLY:  pmc_parentinit, pmc_s_clear_next_array_list, pmc_s_fillbuffer,  &
               pmc_s_getdata_from_buffer, pmc_s_getnextarray,                  &
               pmc_s_setind_and_allocmem, pmc_s_set_active_data_array,         &
               pmc_s_set_dataarray, pmc_s_set_2d_index_list

#endif

    USE salsa_mod,                                                             &
        ONLY:  aerosol_mass, aerosol_number, gconc_2, ibc_salsa_b,             &
               mconc_2, nbins_aerosol,                                         &
               ncomponents_mass, nconc_2, nesting_salsa, ngases_salsa,         &
               salsa_gas, salsa_gases_from_chem

    USE surface_mod,                                                           &
        ONLY:  bc_h, surf_def_h, surf_lsm_h, surf_usm_h

    IMPLICIT NONE

#if defined( __parallel )
#if defined( __mpifh )
    INCLUDE "mpif.h"
#endif
#endif

    PRIVATE
!
!-- Constants
    INTEGER(iwp), PARAMETER ::  child_to_parent = 2   !< Parameter for pmci_parent_datatrans indicating the direction of transfer
    INTEGER(iwp), PARAMETER ::  parent_to_child = 1   !< Parameter for pmci_parent_datatrans indicating the direction of transfer
    INTEGER(iwp), PARAMETER ::  interpolation_scheme_lrsn  = 2  !< Interpolation scheme to be used on lateral boundaries
    INTEGER(iwp), PARAMETER ::  interpolation_scheme_t     = 3  !< Interpolation scheme to be used on top boundary

    REAL(wp), PARAMETER ::  tolefac = 1.0E-6_wp                 !< Relative tolerence for grid-line matching tests and comparisons
!
!-- Coupler setup
    INTEGER(iwp), SAVE      ::  comm_world_nesting    !< Global nesting communicator
    INTEGER(iwp), SAVE      ::  cpl_id  = 1           !< 
    INTEGER(iwp), SAVE      ::  cpl_npe_total         !<
    INTEGER(iwp), SAVE      ::  cpl_parent_id         !<
    
    CHARACTER(LEN=32), SAVE ::  cpl_name              !<

!
!-- Control parameters
    INTEGER(iwp),     SAVE ::  anterpolation_buffer_width = 2       !< Boundary buffer width for anterpolation
    CHARACTER(LEN=7), SAVE ::  nesting_datatransfer_mode = 'mixed'  !< steering parameter for data-transfer mode
    CHARACTER(LEN=8), SAVE ::  nesting_mode = 'two-way'             !< steering parameter for 1- or 2-way nesting
    
    LOGICAL, SAVE ::  nested_run = .FALSE.  !< general switch
    LOGICAL, SAVE ::  rans_mode_parent = .FALSE. !< mode of parent model (.F. - LES mode, .T. - RANS mode)
!
!-- Geometry
    REAL(wp), SAVE, DIMENSION(:), ALLOCATABLE, PUBLIC ::  coord_x            !< Array for the absolute x-coordinates
    REAL(wp), SAVE, DIMENSION(:), ALLOCATABLE, PUBLIC ::  coord_y            !< Array for the absolute y-coordinates
    REAL(wp), SAVE, PUBLIC                            ::  lower_left_coord_x !< x-coordinate of the lower left corner of the domain 
    REAL(wp), SAVE, PUBLIC                            ::  lower_left_coord_y !< y-coordinate of the lower left corner of the domain 
!
!-- Children's parent-grid arrays
    INTEGER(iwp), SAVE, DIMENSION(5), PUBLIC    ::  parent_bound        !< subdomain index bounds for children's parent-grid arrays

    REAL(wp), SAVE, DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  dissc !< Parent-grid array on child domain - dissipation rate
    REAL(wp), SAVE, DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  ec    !< Parent-grid array on child domain - SGS TKE
    REAL(wp), SAVE, DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  ptc   !< Parent-grid array on child domain - potential temperature
    REAL(wp), SAVE, DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  uc    !< Parent-grid array on child domain - velocity component u
    REAL(wp), SAVE, DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  vc    !< Parent-grid array on child domain - velocity component v
    REAL(wp), SAVE, DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  wc    !< Parent-grid array on child domain - velocity component w
    REAL(wp), SAVE, DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  q_c   !< Parent-grid array on child domain - 
    REAL(wp), SAVE, DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  qcc   !< Parent-grid array on child domain - 
    REAL(wp), SAVE, DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  qrc   !< Parent-grid array on child domain - 
    REAL(wp), SAVE, DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  nrc   !< Parent-grid array on child domain - 
    REAL(wp), SAVE, DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  ncc   !< Parent-grid array on child domain - 
    REAL(wp), SAVE, DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  sc    !< Parent-grid array on child domain - 
    INTEGER(idp), SAVE, DIMENSION(:,:), ALLOCATABLE, TARGET, PUBLIC ::  nr_partc    !<
    INTEGER(idp), SAVE, DIMENSION(:,:), ALLOCATABLE, TARGET, PUBLIC ::  part_adrc   !<

    REAL(wp), SAVE, DIMENSION(:,:,:,:), ALLOCATABLE, TARGET ::  chem_spec_c  !< Parent-grid array on child domain - chemical species

    REAL(wp), SAVE, DIMENSION(:,:,:,:), ALLOCATABLE, TARGET ::  aerosol_mass_c    !< Aerosol mass
    REAL(wp), SAVE, DIMENSION(:,:,:,:), ALLOCATABLE, TARGET ::  aerosol_number_c  !< Aerosol number 
    REAL(wp), SAVE, DIMENSION(:,:,:,:), ALLOCATABLE, TARGET ::  salsa_gas_c       !< SALSA gases
!
!-- Grid-spacing ratios.
    INTEGER(iwp), SAVE ::  igsr    !< Integer grid-spacing ratio in i-direction 
    INTEGER(iwp), SAVE ::  jgsr    !< Integer grid-spacing ratio in j-direction 
    INTEGER(iwp), SAVE ::  kgsr    !< Integer grid-spacing ratio in k-direction
!
!-- Global parent-grid index bounds
    INTEGER(iwp), SAVE ::  iplg    !< Leftmost parent-grid array ip index of the whole child domain
    INTEGER(iwp), SAVE ::  iprg    !< Rightmost parent-grid array ip index of the whole child domain
    INTEGER(iwp), SAVE ::  jpsg    !< Southmost parent-grid array jp index of the whole child domain
    INTEGER(iwp), SAVE ::  jpng    !< Northmost parent-grid array jp index of the whole child domain
!
!-- Local parent-grid index bounds. Different sets of index bounds are needed for parent-grid arrays (uc, etc),
!-- for index mapping arrays (iflu, etc) and for work arrays (workarr_lr, etc). This is because these arrays
!-- have different dimensions depending on the location of the subdomain relative to boundaries and corners.
    INTEGER(iwp), SAVE ::  ipl     !< Left index limit for children's parent-grid arrays
    INTEGER(iwp), SAVE ::  ipla    !< Left index limit for allocation of index-mapping and other auxiliary arrays
    INTEGER(iwp), SAVE ::  iplw    !< Left index limit for children's parent-grid work arrays
    INTEGER(iwp), SAVE ::  ipr     !< Right index limit for children's parent-grid arrays
    INTEGER(iwp), SAVE ::  ipra    !< Right index limit for allocation of index-mapping and other auxiliary arrays
    INTEGER(iwp), SAVE ::  iprw    !< Right index limit for children's parent-grid work arrays
    INTEGER(iwp), SAVE ::  jpn     !< North index limit for children's parent-grid arrays
    INTEGER(iwp), SAVE ::  jpna    !< North index limit for allocation of index-mapping and other auxiliary arrays
    INTEGER(iwp), SAVE ::  jpnw    !< North index limit for children's parent-grid work arrays
    INTEGER(iwp), SAVE ::  jps     !< South index limit for children's parent-grid arrays
    INTEGER(iwp), SAVE ::  jpsa    !< South index limit for allocation of index-mapping and other auxiliary arrays
    INTEGER(iwp), SAVE ::  jpsw    !< South index limit for children's parent-grid work arrays
!
!-- Highest prognostic parent-grid k-indices.
    INTEGER(iwp), SAVE ::  kcto     !< Upper bound for k in anterpolation of variables other than w.
    INTEGER(iwp), SAVE ::  kctw     !< Upper bound for k in anterpolation of w.
!
!-- Child-array indices to be precomputed and stored for anterpolation.
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  iflu   !< child index indicating left bound of parent grid box on u-grid
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  ifuu   !< child index indicating right bound of parent grid box on u-grid
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  iflo   !< child index indicating left bound of parent grid box on scalar-grid
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  ifuo   !< child index indicating right bound of parent grid box on scalar-grid
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  jflv   !< child index indicating south bound of parent grid box on v-grid
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  jfuv   !< child index indicating north bound of parent grid box on v-grid
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  jflo   !< child index indicating south bound of parent grid box on scalar-grid
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  jfuo   !< child index indicating north bound of parent grid box on scalar-grid
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  kflw   !< child index indicating lower bound of parent grid box on w-grid
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  kfuw   !< child index indicating upper bound of parent grid box on w-grid
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  kflo   !< child index indicating lower bound of parent grid box on scalar-grid
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  kfuo   !< child index indicating upper bound of parent grid box on scalar-grid
!
!-- Number of child-grid nodes within anterpolation cells to be precomputed for anterpolation.
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::  ijkfc_u  !< number of child grid points contributing to a parent grid
                                                                   !< node in anterpolation, u-grid
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::  ijkfc_v  !< number of child grid points contributing to a parent grid
                                                                   !< node in anterpolation, v-grid
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::  ijkfc_w  !< number of child grid points contributing to a parent grid
                                                                   !< node in anterpolation, w-grid
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::  ijkfc_s  !< number of child grid points contributing to a parent grid
                                                                   !< node in anterpolation, scalar-grid
!    
!-- Work arrays for interpolation and user-defined type definitions for horizontal work-array exchange    
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: workarr_lr
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: workarr_sn
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: workarr_t

    INTEGER(iwp) :: workarr_lr_exchange_type
    INTEGER(iwp) :: workarr_sn_exchange_type
    INTEGER(iwp) :: workarr_t_exchange_type_x
    INTEGER(iwp) :: workarr_t_exchange_type_y
 
    INTEGER(iwp), DIMENSION(3)          ::  parent_grid_info_int    !< Array for communicating the parent-grid dimensions
                                                                    !< to its children.

    REAL(wp), DIMENSION(6)              ::  face_area               !< Surface area of each boundary face 
    REAL(wp), DIMENSION(7)              ::  parent_grid_info_real   !< Array for communicating the real-type parent-grid 
                                                                    !< parameters to its children.

    TYPE parentgrid_def
       INTEGER(iwp)                        ::  nx                 !<
       INTEGER(iwp)                        ::  ny                 !<
       INTEGER(iwp)                        ::  nz                 !<
       REAL(wp)                            ::  dx                 !<
       REAL(wp)                            ::  dy                 !<
       REAL(wp)                            ::  dz                 !<
       REAL(wp)                            ::  lower_left_coord_x !<
       REAL(wp)                            ::  lower_left_coord_y !<
       REAL(wp)                            ::  xend               !<
       REAL(wp)                            ::  yend               !<
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  coord_x            !<
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  coord_y            !<
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  dzu                !<
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  dzw                !<
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  zu                 !<
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  zw                 !<
    END TYPE parentgrid_def

    TYPE(parentgrid_def), SAVE, PUBLIC     ::  pg                 !< Parent-grid information package of type parentgrid_def
!
!-- Variables for particle coupling
    TYPE, PUBLIC :: childgrid_def
       INTEGER(iwp)                        ::  nx                   !<
       INTEGER(iwp)                        ::  ny                   !<
       INTEGER(iwp)                        ::  nz                   !<
       REAL(wp)                            ::  dx                   !<
       REAL(wp)                            ::  dy                   !<
       REAL(wp)                            ::  dz                   !<
       REAL(wp)                            ::  lx_coord, lx_coord_b !<   ! split onto separate lines
       REAL(wp)                            ::  rx_coord, rx_coord_b !<
       REAL(wp)                            ::  sy_coord, sy_coord_b !<
       REAL(wp)                            ::  ny_coord, ny_coord_b !<
       REAL(wp)                            ::  uz_coord, uz_coord_b !<
    END TYPE childgrid_def

    TYPE(childgrid_def), SAVE, ALLOCATABLE, DIMENSION(:), PUBLIC ::  childgrid  !<

    INTEGER(idp), ALLOCATABLE,DIMENSION(:,:), PUBLIC,TARGET ::  nr_part  !<
    INTEGER(idp), ALLOCATABLE,DIMENSION(:,:), PUBLIC,TARGET ::  part_adr !<

   
    INTERFACE pmci_boundary_conds
       MODULE PROCEDURE pmci_boundary_conds
    END INTERFACE pmci_boundary_conds
    
    INTERFACE pmci_check_setting_mismatches
       MODULE PROCEDURE pmci_check_setting_mismatches
    END INTERFACE

    INTERFACE pmci_child_initialize
       MODULE PROCEDURE pmci_child_initialize
    END INTERFACE

    INTERFACE pmci_synchronize
       MODULE PROCEDURE pmci_synchronize
    END INTERFACE

    INTERFACE pmci_datatrans
       MODULE PROCEDURE pmci_datatrans
    END INTERFACE pmci_datatrans

    INTERFACE pmci_ensure_nest_mass_conservation
       MODULE PROCEDURE pmci_ensure_nest_mass_conservation
    END INTERFACE pmci_ensure_nest_mass_conservation

    INTERFACE pmci_ensure_nest_mass_conservation_vertical
       MODULE PROCEDURE pmci_ensure_nest_mass_conservation_vertical
    END INTERFACE pmci_ensure_nest_mass_conservation_vertical

    INTERFACE pmci_init
       MODULE PROCEDURE pmci_init
    END INTERFACE

    INTERFACE pmci_modelconfiguration
       MODULE PROCEDURE pmci_modelconfiguration
    END INTERFACE

    INTERFACE pmci_parent_initialize
       MODULE PROCEDURE pmci_parent_initialize
    END INTERFACE

    INTERFACE get_number_of_children
       MODULE PROCEDURE get_number_of_children
    END  INTERFACE get_number_of_children

    INTERFACE get_childid
       MODULE PROCEDURE get_childid
    END  INTERFACE get_childid

    INTERFACE get_child_edges
       MODULE PROCEDURE get_child_edges
    END  INTERFACE get_child_edges

    INTERFACE get_child_gridspacing
       MODULE PROCEDURE get_child_gridspacing
    END  INTERFACE get_child_gridspacing

    INTERFACE pmci_set_swaplevel
       MODULE PROCEDURE pmci_set_swaplevel
    END INTERFACE pmci_set_swaplevel

    PUBLIC child_to_parent, comm_world_nesting, cpl_id, nested_run,                                 &
           nesting_datatransfer_mode, nesting_mode, parent_to_child, rans_mode_parent

    PUBLIC pmci_boundary_conds
    PUBLIC pmci_child_initialize
    PUBLIC pmci_datatrans
    PUBLIC pmci_init
    PUBLIC pmci_modelconfiguration
    PUBLIC pmci_parent_initialize
    PUBLIC pmci_synchronize
    PUBLIC pmci_set_swaplevel
    PUBLIC get_number_of_children, get_childid, get_child_edges, get_child_gridspacing
    PUBLIC pmci_ensure_nest_mass_conservation
    PUBLIC pmci_ensure_nest_mass_conservation_vertical
    
 CONTAINS


 SUBROUTINE pmci_init( world_comm )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(OUT) ::  world_comm   !<

#if defined( __parallel )

    INTEGER(iwp) ::  pmc_status   !<


    CALL pmc_init_model( world_comm, nesting_datatransfer_mode, nesting_mode,                       &
                         anterpolation_buffer_width, pmc_status )

    IF ( pmc_status == pmc_no_namelist_found )  THEN
!
!--    This is not a nested run
       world_comm = MPI_COMM_WORLD
       cpl_id     = 1
       cpl_name   = ""

       RETURN

    ENDIF
!
!-- Check steering parameter values
    IF ( TRIM( nesting_mode ) /= 'one-way'  .AND.                                                   &
         TRIM( nesting_mode ) /= 'two-way'  .AND.                                                   &
         TRIM( nesting_mode ) /= 'vertical' )                                                       &
    THEN
       message_string = 'illegal nesting mode: ' // TRIM( nesting_mode )
       CALL message( 'pmci_init', 'PA0417', 3, 2, 0, 6, 0 )
    ENDIF

    IF ( TRIM( nesting_datatransfer_mode ) /= 'cascade'  .AND.                                      &
         TRIM( nesting_datatransfer_mode ) /= 'mixed'    .AND.                                      &
         TRIM( nesting_datatransfer_mode ) /= 'overlap' )                                           &
    THEN
       message_string = 'illegal nesting datatransfer mode: ' // TRIM( nesting_datatransfer_mode )
       CALL message( 'pmci_init', 'PA0418', 3, 2, 0, 6, 0 )
    ENDIF
!
!-- Set the general steering switch which tells PALM that its a nested run
    nested_run = .TRUE.
!
!-- Get some variables required by the pmc-interface (and in some cases in the
!-- PALM code out of the pmci) out of the pmc-core
    CALL pmc_get_model_info( comm_world_nesting = comm_world_nesting,                               &
                             cpl_id = cpl_id, cpl_parent_id = cpl_parent_id,                        &
                             cpl_name = cpl_name, npe_total = cpl_npe_total,                        &
                             lower_left_x = lower_left_coord_x,                                     &
                             lower_left_y = lower_left_coord_y )
!
!-- Set the steering switch which tells the models that they are nested (of
!-- course the root domain is not nested)
    IF ( .NOT.  pmc_is_rootmodel() )  THEN
       child_domain = .TRUE.
       WRITE( coupling_char, '(A2,I2.2)') '_N', cpl_id
    ENDIF

!
!-- Message that communicators for nesting are initialized.
!-- Attention: myid has been set at the end of pmc_init_model in order to
!-- guarantee that only PE0 of the root domain does the output.
    CALL location_message( 'initialize model nesting', 'finished' )
!
!-- Reset myid to its default value
    myid = 0
#else
!
!-- Nesting cannot be used in serial mode. cpl_id is set to root domain (1)
!-- because no location messages would be generated otherwise.
!-- world_comm is given a dummy value to avoid compiler warnings (INTENT(OUT)
!-- must get an explicit value).
!-- Note that this branch is only to avoid compiler warnings. The actual 
!-- execution never reaches here because the call of this subroutine is 
!-- already enclosed by  #if defined( __parallel ).
    cpl_id     = 1
    nested_run = .FALSE.
    world_comm = 1
#endif

 END SUBROUTINE pmci_init



 SUBROUTINE pmci_modelconfiguration

    IMPLICIT NONE

    INTEGER(iwp) ::  ncpl   !<  number of nest domains

    
#if defined( __parallel )
    CALL location_message( 'setup the nested model configuration', 'start' )
    CALL cpu_log( log_point_s(79), 'pmci_model_config', 'start' )
!
!-- Compute absolute coordinates for all models
    CALL pmci_setup_coordinates         ! CONTAIN THIS 
!
!-- Determine the number of coupled arrays
    CALL pmci_num_arrays                ! CONTAIN THIS
!
!-- Initialize the child (must be called before pmc_setup_parent)
!-- Klaus, extend this comment to explain why it must be called before    
    CALL pmci_setup_child               ! CONTAIN THIS
!
!-- Initialize PMC parent
    CALL pmci_setup_parent              ! CONTAIN THIS
!
!-- Check for mismatches between settings of master and child variables
!-- (e.g., all children have to follow the end_time settings of the root master)
    CALL pmci_check_setting_mismatches  ! CONTAIN THIS
!
!-- Set flag file for combine_plot_fields for precessing the nest output data
    OPEN( 90, FILE='3DNESTING', FORM='FORMATTED' )
    CALL pmc_get_model_info( ncpl = ncpl )
    WRITE( 90, '(I2)' )  ncpl
    CLOSE( 90 )

    CALL cpu_log( log_point_s(79), 'pmci_model_config', 'stop' )
    CALL location_message( 'setup the nested model configuration', 'finished' )
#endif

 END SUBROUTINE pmci_modelconfiguration



 SUBROUTINE pmci_setup_parent

#if defined( __parallel )
    IMPLICIT NONE

    INTEGER(iwp) ::  child_id           !< Child id-number for the child m
    INTEGER(iwp) ::  ierr               !< MPI-error code 
    INTEGER(iwp) ::  kp                 !< Parent-grid index n the z-direction
    INTEGER(iwp) ::  lb = 1             !< Running index for aerosol size bins
    INTEGER(iwp) ::  lc = 1             !< Running index for aerosol mass bins
    INTEGER(iwp) ::  lg = 1             !< Running index for SALSA gases
    INTEGER(iwp) ::  m                  !< Loop index over all children of the current parent
    INTEGER(iwp) ::  msib               !< Loop index over all other children than m in case of siblings (parallel children)
    INTEGER(iwp) ::  n = 1              !< Running index for chemical species
    INTEGER(iwp) ::  nx_child           !< Number of child-grid points in the x-direction
    INTEGER(iwp) ::  ny_child           !< Number of child-grid points in the y-direction
    INTEGER(iwp) ::  nz_child           !< Number of child-grid points in the z-direction
    INTEGER(iwp) ::  sibling_id         !< Child id-number for the child msib (sibling of child m)
    
    INTEGER(iwp), DIMENSION(3) ::  child_grid_dim  !< Array for receiving the child-grid dimensions from the children
    
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  child_x_left   !< Minimum x-coordinate of the child domain including the ghost
                                                           !< point layers 
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  child_x_right  !< Maximum x-coordinate of the child domain including the ghost
                                                           !< point layers   
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  child_y_south  !< Minimum y-coordinate of the child domain including the ghost
                                                           !< point layers
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  child_y_north  !< Maximum y-coordinate of the child domain including the ghost
                                                           !< point layers
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  child_coord_x  !< Child domain x-coordinate array
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  child_coord_y  !< Child domain y-coordinate array
    
    REAL(wp), DIMENSION(5) ::  child_grid_info  !< Array for receiving the child-grid spacings etc from the children
    
    REAL(wp) ::  child_height         !< Height of the child domain defined on the child side as zw(nzt+1)
    REAL(wp) ::  dx_child             !< Child-grid spacing in the x-direction
    REAL(wp) ::  dy_child             !< Child-grid spacing in the y-direction
    REAL(wp) ::  dz_child             !< Child-grid spacing in the z-direction
    REAL(wp) ::  left_limit           !< Left limit for the absolute x-coordinate of the child left boundary 
    REAL(wp) ::  north_limit          !< North limit for the absolute y-coordinate of the child north boundary
    REAL(wp) ::  right_limit          !< Right limit for the absolute x-coordinate of the child right boundary
    REAL(wp) ::  south_limit          !< South limit for the absolute y-coordinate of the child south boundary 
    REAL(wp) ::  upper_right_coord_x  !< Absolute x-coordinate of the upper right corner of the child domain 
    REAL(wp) ::  upper_right_coord_y  !< Absolute y-coordinate of the upper right corner of the child domain  
    REAL(wp) ::  xez                  !< Minimum separation in the x-direction required between the child and
                                      !< parent boundaries (left or right)
    REAL(wp) ::  yez                  !< Minimum separation in the y-direction required between the child and
                                      !< parent boundaries (south or north)
    REAL(wp)     ::  tolex            !< Tolerance for grid-line matching in x-direction
    REAL(wp)     ::  toley            !< Tolerance for grid-line matching in y-direction
    REAL(wp)     ::  tolez            !< Tolerance for grid-line matching in z-direction    

    CHARACTER(LEN=32) ::  myname      !< String for variable name such as 'u'

    LOGICAL :: m_left_in_msib         !< Logical auxiliary parameter for the overlap test: true if the left border
                                      !< of the child m is within the x-range of the child msib
    LOGICAL :: m_right_in_msib        !< Logical auxiliary parameter for the overlap test: true if the right border
                                      !< of the child m is within the x-range of the child msib 
    LOGICAL :: msib_left_in_m         !< Logical auxiliary parameter for the overlap test: true if the left border
                                      !< of the child msib is within the x-range of the child m
    LOGICAL :: msib_right_in_m        !< Logical auxiliary parameter for the overlap test: true if the right border
                                      !< of the child msib is within the x-range of the child m
    LOGICAL :: m_south_in_msib        !< Logical auxiliary parameter for the overlap test: true if the south border
                                      !< of the child m is within the y-range of the child msib
    LOGICAL :: m_north_in_msib        !< Logical auxiliary parameter for the overlap test: true if the north border
                                      !< of the child m is within the y-range of the child msib 
    LOGICAL :: msib_south_in_m        !< Logical auxiliary parameter for the overlap test: true if the south border
                                      !< of the child msib is within the y-range of the child m
    LOGICAL :: msib_north_in_m        !< Logical auxiliary parameter for the overlap test: true if the north border
                                      !< of the child msib is within the y-range of the child m

!
!-- Grid-line tolerances.
    tolex = tolefac * dx
    toley = tolefac * dy
    tolez = tolefac * dz(1)    
!
!-- Initialize the current pmc parent.
    CALL pmc_parentinit
!
!-- Corners of all children of the present parent. Note that 
!-- SIZE( pmc_parent_for_child ) = 1 if we have no children.
    IF ( ( SIZE( pmc_parent_for_child ) - 1 > 0 )  .AND.  myid == 0 )  THEN 
       ALLOCATE( child_x_left(1:SIZE( pmc_parent_for_child ) - 1) )
       ALLOCATE( child_x_right(1:SIZE( pmc_parent_for_child ) - 1) )
       ALLOCATE( child_y_south(1:SIZE( pmc_parent_for_child ) - 1) )
       ALLOCATE( child_y_north(1:SIZE( pmc_parent_for_child ) - 1) )
    ENDIF
    IF ( ( SIZE( pmc_parent_for_child ) - 1 > 0 ) )  THEN
       ALLOCATE( childgrid(1:SIZE( pmc_parent_for_child ) - 1) )
    ENDIF
!
!-- Get coordinates from all children and check that the children match the parent
!-- domain and each others. Note that SIZE( pmc_parent_for_child ) = 1 
!-- if we have no children, thence the loop is not executed at all. 
    DO  m = 1, SIZE( pmc_parent_for_child ) - 1

       child_id = pmc_parent_for_child(m)

       IF ( myid == 0 )  THEN

          CALL pmc_recv_from_child( child_id, child_grid_dim,  SIZE(child_grid_dim), 0, 123, ierr )
          CALL pmc_recv_from_child( child_id, child_grid_info, SIZE(child_grid_info), 0, 124, ierr )
         
          nx_child     = child_grid_dim(1)
          ny_child     = child_grid_dim(2)
          dx_child     = child_grid_info(3)
          dy_child     = child_grid_info(4)
          dz_child     = child_grid_info(5)
          child_height = child_grid_info(1)
!
!--       Find the highest child-domain level in the parent grid for the reduced z transfer
          DO  kp = 1, nzt                 
             IF ( zw(kp) - child_height > tolez )  THEN                   
                nz_child = kp
                EXIT
             ENDIF
          ENDDO
!    
!--       Get absolute coordinates from the child
          ALLOCATE( child_coord_x(-nbgp:nx_child+nbgp) )
          ALLOCATE( child_coord_y(-nbgp:ny_child+nbgp) )
         
          CALL pmc_recv_from_child( child_id, child_coord_x, SIZE( child_coord_x ), 0, 11, ierr )
          CALL pmc_recv_from_child( child_id, child_coord_y, SIZE( child_coord_y ), 0, 12, ierr )
         
          parent_grid_info_real(1) = lower_left_coord_x
          parent_grid_info_real(2) = lower_left_coord_y
          parent_grid_info_real(3) = dx
          parent_grid_info_real(4) = dy

          upper_right_coord_x      = lower_left_coord_x + ( nx + 1 ) * dx
          upper_right_coord_y      = lower_left_coord_y + ( ny + 1 ) * dy
          parent_grid_info_real(5) = upper_right_coord_x
          parent_grid_info_real(6) = upper_right_coord_y
          parent_grid_info_real(7) = dz(1)

          parent_grid_info_int(1)  = nx
          parent_grid_info_int(2)  = ny
          parent_grid_info_int(3)  = nz_child
!
!--       Check that the child domain matches its parent domain. 
          IF ( nesting_mode == 'vertical' )  THEN
!
!--          In case of vertical nesting, the lateral boundaries must match exactly. 
             right_limit = upper_right_coord_x
             north_limit = upper_right_coord_y
             IF ( ABS( child_coord_x(nx_child+1) - right_limit ) > tolex )  THEN
                WRITE ( message_string, "(a,i2,a)" ) 'nested child (id: ',child_id,                 &
                     ') domain right edge does not match its parent right edge'
                CALL message( 'pmci_setup_parent', 'PA0425', 3, 2, 0, 6, 0 )
             ENDIF
             IF ( ABS( child_coord_y(ny_child+1) - north_limit ) > toley )  THEN
                WRITE ( message_string, "(a,i2,a)" ) 'nested child (id: ',child_id,                 &
                     ') domain north edge does not match its parent north edge'
                CALL message( 'pmci_setup_parent', 'PA0425', 3, 2, 0, 6, 0 )
             ENDIF
          ELSE       
!
!--          In case of 3-D nesting, check that the child domain is completely 
!--          inside its parent domain. 
             xez = ( nbgp + 1 ) * dx 
             yez = ( nbgp + 1 ) * dy 
             left_limit  = lower_left_coord_x + xez
             right_limit = upper_right_coord_x - xez
             south_limit = lower_left_coord_y + yez
             north_limit = upper_right_coord_y - yez
             IF ( left_limit - child_coord_x(0) > tolex )  THEN
                WRITE ( message_string, "(a,i2,a)" ) 'nested child (id: ',child_id,                 &
                     ') domain does not fit in its parent domain, left edge is either too ' //      &
                     'close or outside its parent left edge' 
                CALL message( 'pmci_setup_parent', 'PA0425', 3, 2, 0, 6, 0 )
             ENDIF
             IF ( child_coord_x(nx_child+1) - right_limit > tolex )  THEN
                WRITE ( message_string, "(a,i2,a)" ) 'nested child (id: ',child_id,                 &
                     ') domain does not fit in its parent domain, right edge is either too ' //     &
                     'close or outside its parent right edge' 
                CALL message( 'pmci_setup_parent', 'PA0425', 3, 2, 0, 6, 0 )
             ENDIF
             IF ( south_limit - child_coord_y(0) > toley )  THEN
                WRITE ( message_string, "(a,i2,a)" ) 'nested child (id: ',child_id,                 &
                     ') domain does not fit in its parent domain, south edge is either too ' //     &
                     'close or outside its parent south edge' 
                CALL message( 'pmci_setup_parent', 'PA0425', 3, 2, 0, 6, 0 )
             ENDIF
             IF ( child_coord_y(ny_child+1) - north_limit > toley )  THEN
                WRITE ( message_string, "(a,i2,a)" ) 'nested child (id: ',child_id,                 &
                     ') domain does not fit in its parent domain, north edge is either too ' //     &
                     'close or outside its parent north edge' 
                CALL message( 'pmci_setup_parent', 'PA0425', 3, 2, 0, 6, 0 )
             ENDIF
          ENDIF
!             
!--       Child domain must be lower than the parent domain such that the top ghost
!--       layer of the child grid does not exceed the parent domain top boundary.
          IF ( child_height - zw(nzt) > tolez ) THEN
             WRITE ( message_string, "(a,i2,a)" ) 'nested child (id: ',child_id,                    &
                     ') domain does not fit in its parent domain, top edge is either too ' //       &
                     'close or above its parent top edge' 
             CALL message( 'pmci_setup_parent', 'PA0425', 3, 2, 0, 6, 0 )
          ENDIF
!
!--       If parallel child domains (siblings) do exist ( m > 1 ), 
!--       check that they do not overlap.
          child_x_left(m)  = child_coord_x(-nbgp)
          child_x_right(m) = child_coord_x(nx_child+nbgp)
          child_y_south(m) = child_coord_y(-nbgp)
          child_y_north(m) = child_coord_y(ny_child+nbgp)

          IF ( nesting_mode /= 'vertical' )  THEN
!
!--          Note that the msib-loop is executed only if ( m > 1 ).  
!--          Also note that the tests have to be made both ways (m vs msib and msib vs m)
!--          in order to detect all the possible overlap situations.
             DO  msib = 1, m - 1
!
!--             Set some logical auxiliary parameters to simplify the IF-condition.                  
                m_left_in_msib  = ( child_x_left(m)  >= child_x_left(msib)  - tolex )  .AND.        &
                                  ( child_x_left(m)  <= child_x_right(msib) + tolex )
                m_right_in_msib = ( child_x_right(m) >= child_x_left(msib)  - tolex )  .AND.        &
                                  ( child_x_right(m) <= child_x_right(msib) + tolex )
                msib_left_in_m  = ( child_x_left(msib)  >= child_x_left(m)  - tolex )  .AND.        &
                                  ( child_x_left(msib)  <= child_x_right(m) + tolex )
                msib_right_in_m = ( child_x_right(msib) >= child_x_left(m)  - tolex )  .AND.        &
                                  ( child_x_right(msib) <= child_x_right(m) + tolex )
                m_south_in_msib = ( child_y_south(m) >= child_y_south(msib) - toley )  .AND.        &
                                  ( child_y_south(m) <= child_y_north(msib) + toley )
                m_north_in_msib = ( child_y_north(m) >= child_y_south(msib) - toley )  .AND.        &
                                  ( child_y_north(m) <= child_y_north(msib) + toley )
                msib_south_in_m = ( child_y_south(msib) >= child_y_south(m) - toley )  .AND.        &
                                  ( child_y_south(msib) <= child_y_north(m) + toley )
                msib_north_in_m = ( child_y_north(msib) >= child_y_south(m) - toley )  .AND.        &
                                  ( child_y_north(msib) <= child_y_north(m) + toley )
                
                IF ( ( m_left_in_msib  .OR.  m_right_in_msib  .OR.                                  &
                       msib_left_in_m  .OR.  msib_right_in_m )                                      &
                     .AND.                                                                          &
                     ( m_south_in_msib  .OR.  m_north_in_msib  .OR.                                 &
                       msib_south_in_m  .OR.  msib_north_in_m ) )  THEN
                   sibling_id = pmc_parent_for_child(msib)
                   WRITE ( message_string, "(a,i2,a,i2,a)" ) 'nested parallel child domains (ids: ',&
                        child_id, ' and ', sibling_id, ') overlap'
                   CALL message( 'pmci_setup_parent', 'PA0426', 3, 2, 0, 6, 0 )
                ENDIF

             ENDDO
          ENDIF          

          CALL pmci_set_child_edge_coords

          DEALLOCATE( child_coord_x )
          DEALLOCATE( child_coord_y )
!
!--       Send information about operating mode (LES or RANS) to child. This will be 
!--       used to control TKE nesting and setting boundary conditions properly.
          CALL pmc_send_to_child( child_id, rans_mode, 1, 0, 19, ierr )  
!
!--       Send parent grid information to child
          CALL pmc_send_to_child( child_id, parent_grid_info_real,                                  &
                                  SIZE( parent_grid_info_real ), 0, 21,                             &
                                  ierr )
          CALL pmc_send_to_child( child_id, parent_grid_info_int,  3, 0,                            &
                                  22, ierr )
!
!--       Send local grid to child
          CALL pmc_send_to_child( child_id, coord_x, nx+1+2*nbgp, 0, 24,                            &
                                  ierr )
          CALL pmc_send_to_child( child_id, coord_y, ny+1+2*nbgp, 0, 25,                            &
                                  ierr )
!
!--       Also send the dzu-, dzw-, zu- and zw-arrays here
          CALL pmc_send_to_child( child_id, dzu, nz_child + 1, 0, 26, ierr )
          CALL pmc_send_to_child( child_id, dzw, nz_child + 1, 0, 27, ierr )
          CALL pmc_send_to_child( child_id, zu,  nz_child + 2, 0, 28, ierr )
          CALL pmc_send_to_child( child_id, zw,  nz_child + 2, 0, 29, ierr )
          
       ENDIF  ! ( myid == 0 ) 

       CALL MPI_BCAST( nz_child, 1, MPI_INTEGER, 0, comm2d, ierr )

       CALL MPI_BCAST( childgrid(m), STORAGE_SIZE(childgrid(1))/8, MPI_BYTE, 0, comm2d, ierr )
!
!--    Set up the index-list which is an integer array that maps the child index space on 
!--    the parent index- and subdomain spaces.
       CALL pmci_create_index_list
!
!--    Include couple arrays into parent content.
!--    The adresses of the PALM 2D or 3D array (here parent grid) which are candidates
!--    for coupling are stored once into the pmc context. While data transfer, the array do not
!--    have to be specified again
       CALL pmc_s_clear_next_array_list
       DO WHILE ( pmc_s_getnextarray( child_id, myname ) )
          IF ( INDEX( TRIM( myname ), 'chem_' ) /= 0 )  THEN             
             CALL pmci_set_array_pointer( myname, child_id = child_id, nz_child = nz_child, n = n )
             n = n + 1  
          ELSEIF ( INDEX( TRIM( myname ), 'an_' ) /= 0 )  THEN
             CALL pmci_set_array_pointer( myname, child_id = child_id, nz_child = nz_child, n = lb )
             lb = lb + 1 
          ELSEIF ( INDEX( TRIM( myname ), 'am_' ) /= 0 )  THEN
             CALL pmci_set_array_pointer( myname, child_id = child_id, nz_child = nz_child, n = lc )
             lc = lc + 1 
          ELSEIF ( INDEX( TRIM( myname ), 'sg_' ) /= 0  .AND.  .NOT. salsa_gases_from_chem )  THEN
             CALL pmci_set_array_pointer( myname, child_id = child_id, nz_child = nz_child, n = lg )
             lg = lg + 1
          ELSE
             CALL pmci_set_array_pointer( myname, child_id = child_id, nz_child = nz_child )
          ENDIF
       ENDDO

       CALL pmc_s_setind_and_allocmem( child_id )
       
    ENDDO  ! m

    IF ( ( SIZE( pmc_parent_for_child ) - 1 > 0 ) .AND. myid == 0 )  THEN
       DEALLOCATE( child_x_left )
       DEALLOCATE( child_x_right )
       DEALLOCATE( child_y_south )
       DEALLOCATE( child_y_north )
    ENDIF

    
 CONTAINS


    SUBROUTINE pmci_create_index_list

       IMPLICIT NONE

       INTEGER(iwp) ::  ilist            !< Index-list index running over the child's parent-grid jc,ic-space
       INTEGER(iwp) ::  index_list_size  !< Dimension 2 of the array index_list
       INTEGER(iwp) ::  ierr             !< MPI error code
       INTEGER(iwp) ::  ip               !< Running parent-grid index on the child domain in the x-direction 
       INTEGER(iwp) ::  jp               !< Running parent-grid index on the child domain in the y-direction
       INTEGER(iwp) ::  n                !< Running index over child subdomains
       INTEGER(iwp) ::  nrx              !< Parent subdomain dimension in the x-direction
       INTEGER(iwp) ::  nry              !< Parent subdomain dimension in the y-direction
       INTEGER(iwp) ::  pex              !< Two-dimensional subdomain (pe) index in the x-direction
       INTEGER(iwp) ::  pey              !< Two-dimensional subdomain (pe) index in the y-direction
       INTEGER(iwp) ::  parent_pe        !< Parent subdomain index (one-dimensional)

       INTEGER(iwp), DIMENSION(2) ::  pe_indices_2d                                  !< Array for two-dimensional subdomain (pe)
                                                                                     !< indices needed for MPI_CART_RANK
       INTEGER(iwp), DIMENSION(2) ::  size_of_childs_parent_grid_bounds_all          !< Dimensions of childs_parent_grid_bounds_all
       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE  ::  childs_parent_grid_bounds_all  !< Array that contains the child's
                                                                                     !< parent-grid index bounds for all its
                                                                                     !< subdomains (pes)
       INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE  ::  index_list                     !< Array that maps the child index space on 
                                                                                     !< the parent index- and subdomain spaces
       
       IF ( myid == 0 )  THEN
          
          CALL pmc_recv_from_child( child_id, size_of_childs_parent_grid_bounds_all,                &
                                    2, 0, 40, ierr )
          ALLOCATE( childs_parent_grid_bounds_all(size_of_childs_parent_grid_bounds_all(1),         &
                                                  size_of_childs_parent_grid_bounds_all(2)) )
          CALL pmc_recv_from_child( child_id, childs_parent_grid_bounds_all,                        &
                                    SIZE( childs_parent_grid_bounds_all ), 0, 41, ierr )
!
!--       Compute size (dimension) of the index_list.
          index_list_size = 0         
          DO  n = 1, size_of_childs_parent_grid_bounds_all(2)
             index_list_size = index_list_size +                                                    &
                  ( childs_parent_grid_bounds_all(4,n) - childs_parent_grid_bounds_all(3,n) + 1 ) * &
                  ( childs_parent_grid_bounds_all(2,n) - childs_parent_grid_bounds_all(1,n) + 1 )
          ENDDO

          ALLOCATE( index_list(6,index_list_size) )

          nrx = nxr - nxl + 1
          nry = nyn - nys + 1
          ilist = 0
!
!--       Loop over all children PEs
          DO  n = 1, size_of_childs_parent_grid_bounds_all(2)           ! 
!
!--          Subspace along y required by actual child PE
             DO  jp = childs_parent_grid_bounds_all(3,n), childs_parent_grid_bounds_all(4,n)  ! jp = jps, jpn of child PE# n 
!
!--             Subspace along x required by actual child PE
                DO  ip = childs_parent_grid_bounds_all(1,n), childs_parent_grid_bounds_all(2,n)  ! ip = ipl, ipr of child PE# n 

                   pex = ip / nrx
                   pey = jp / nry
                   pe_indices_2d(1) = pex
                   pe_indices_2d(2) = pey
                   CALL MPI_CART_RANK( comm2d, pe_indices_2d, parent_pe, ierr )
                  
                   ilist = ilist + 1
!
!--                First index in parent array  ! TO_DO: Klaus, please explain better
                   index_list(1,ilist) = ip - ( pex * nrx ) + 1 + nbgp
!
!--                Second index in parent array  ! TO_DO: Klaus, please explain better
                   index_list(2,ilist) = jp - ( pey * nry ) + 1 + nbgp
!
!--                x index of child's parent grid 
                   index_list(3,ilist) = ip - childs_parent_grid_bounds_all(1,n) + 1
!
!--                y index of child's parent grid
                   index_list(4,ilist) = jp - childs_parent_grid_bounds_all(3,n) + 1
!
!--                PE number of child
                   index_list(5,ilist) = n - 1
!
!--                PE number of parent
                   index_list(6,ilist) = parent_pe

                ENDDO
             ENDDO
          ENDDO
!
!--       TO_DO: Klaus: comment what is done here
          CALL pmc_s_set_2d_index_list( child_id, index_list(:,1:ilist) )

       ELSE
!
!--       TO_DO: Klaus: comment why this dummy allocation is required
          ALLOCATE( index_list(6,1) )
          CALL pmc_s_set_2d_index_list( child_id, index_list )
       ENDIF

       DEALLOCATE(index_list)

     END SUBROUTINE pmci_create_index_list



     SUBROUTINE pmci_set_child_edge_coords
        IMPLICIT  NONE

        INTEGER(iwp) ::  nbgp_lpm = 1  !< Number of ghost-point layers used for lpm (Klaus, is this correct?)

        
        nbgp_lpm = MIN( nbgp_lpm, nbgp )

        childgrid(m)%nx = nx_child
        childgrid(m)%ny = ny_child
        childgrid(m)%nz = nz_child
        childgrid(m)%dx = dx_child
        childgrid(m)%dy = dy_child
        childgrid(m)%dz = dz_child

        childgrid(m)%lx_coord   = child_coord_x(0)
        childgrid(m)%lx_coord_b = child_coord_x(-nbgp_lpm)
        childgrid(m)%rx_coord   = child_coord_x(nx_child) + dx_child
        childgrid(m)%rx_coord_b = child_coord_x(nx_child+nbgp_lpm) + dx_child
        childgrid(m)%sy_coord   = child_coord_y(0)
        childgrid(m)%sy_coord_b = child_coord_y(-nbgp_lpm)
        childgrid(m)%ny_coord   = child_coord_y(ny_child) + dy_child
        childgrid(m)%ny_coord_b = child_coord_y(ny_child+nbgp_lpm) + dy_child
        childgrid(m)%uz_coord   = child_grid_info(2)
        childgrid(m)%uz_coord_b = child_grid_info(1)

     END SUBROUTINE pmci_set_child_edge_coords

#endif
 END SUBROUTINE pmci_setup_parent



 SUBROUTINE pmci_setup_child

#if defined( __parallel )
    IMPLICIT NONE

    INTEGER(iwp) ::  ierr                          !< MPI error code
    INTEGER(iwp) ::  lb                            !< Running index for aerosol size bins
    INTEGER(iwp) ::  lc                            !< Running index for aerosol mass bins
    INTEGER(iwp) ::  lg                            !< Running index for SALSA gases
    INTEGER(iwp) ::  n                             !< Running index for number of chemical species
    INTEGER(iwp), DIMENSION(3) ::  child_grid_dim  !< Array for sending the child-grid dimensions to parent 

    REAL(wp), DIMENSION(5) ::  child_grid_info     !< Array for sending the child-grid spacings etc to parent 
         
    CHARACTER( LEN=da_namelen ) ::  myname         !< Name of the variable to be coupled
    CHARACTER(LEN=5) ::  salsa_char                !< Name extension for the variable name in case of SALSA variable
    
!
!-- Child setup
!-- Root model does not have a parent and is not a child, therefore no child setup on root model
    IF ( .NOT. pmc_is_rootmodel() )  THEN
!
!--    KLaus, add a description here what pmc_childinit does       
       CALL pmc_childinit
!
!--    The arrays, which actually will be exchanged between child and parent
!--    are defined Here AND ONLY HERE.
!--    If a variable is removed, it only has to be removed from here.
!--    Please check, if the arrays are in the list of POSSIBLE exchange arrays
!--    in subroutines:
!--    pmci_set_array_pointer (for parent arrays)
!--    pmci_create_childs_parent_grid_arrays (for child's parent-grid arrays)
       CALL pmc_set_dataarray_name( 'parent', 'u', 'child', 'u', ierr )
       CALL pmc_set_dataarray_name( 'parent', 'v', 'child', 'v', ierr )
       CALL pmc_set_dataarray_name( 'parent', 'w', 'child', 'w', ierr )
!
!--    Set data array name for TKE. Please note, nesting of TKE is actually
!--    only done if both parent and child are in LES or in RANS mode. Due to 
!--    design of model coupler, however, data array names must be already 
!--    available at this point.
       IF ( (        rans_mode_parent  .AND.         rans_mode )  .OR.                              &
            (  .NOT. rans_mode_parent  .AND.  .NOT.  rans_mode  .AND.                               &
               .NOT. constant_diffusion ) )  THEN
          CALL pmc_set_dataarray_name( 'parent', 'e', 'child', 'e', ierr )
       ENDIF
!
!--    Nesting of dissipation rate only if both parent and child are in RANS
!--    mode and TKE-epsilon closure is applied. Please see also comment for TKE
!--    above.
       IF ( rans_mode_parent  .AND.  rans_mode  .AND.  rans_tke_e )  THEN
          CALL pmc_set_dataarray_name( 'parent', 'diss', 'child', 'diss', ierr )
       ENDIF

       IF ( .NOT. neutral )  THEN
          CALL pmc_set_dataarray_name( 'parent', 'pt' ,'child', 'pt', ierr )
       ENDIF

       IF ( humidity )  THEN

          CALL pmc_set_dataarray_name( 'parent', 'q', 'child', 'q', ierr )

          IF ( bulk_cloud_model  .AND.  microphysics_morrison )  THEN
            CALL pmc_set_dataarray_name( 'parent', 'qc', 'child', 'qc', ierr )  
            CALL pmc_set_dataarray_name( 'parent', 'nc', 'child', 'nc', ierr ) 
          ENDIF

          IF ( bulk_cloud_model  .AND.  microphysics_seifert )  THEN
             CALL pmc_set_dataarray_name( 'parent', 'qr', 'child', 'qr', ierr )
             CALL pmc_set_dataarray_name( 'parent', 'nr', 'child', 'nr', ierr ) 
          ENDIF
     
       ENDIF

       IF ( passive_scalar )  THEN
          CALL pmc_set_dataarray_name( 'parent', 's', 'child', 's', ierr )
       ENDIF

       IF ( particle_advection )  THEN
          CALL pmc_set_dataarray_name( 'parent', 'nr_part', 'child', 'nr_part', ierr )
          CALL pmc_set_dataarray_name( 'parent', 'part_adr', 'child', 'part_adr', ierr )
       ENDIF
       
       IF ( air_chemistry  .AND.  nesting_chem )  THEN
          DO n = 1, nspec
             CALL pmc_set_dataarray_name( 'parent', 'chem_' // TRIM( chem_species(n)%name ),        &
                                          'child',  'chem_' // TRIM( chem_species(n)%name ), ierr )
          ENDDO 
       ENDIF

       IF ( salsa  .AND.  nesting_salsa )  THEN
          DO  lb = 1, nbins_aerosol
             WRITE(salsa_char,'(i0)') lb
             CALL pmc_set_dataarray_name( 'parent', 'an_' // TRIM( salsa_char ),                    &
                                          'child',  'an_' // TRIM( salsa_char ), ierr )
          ENDDO
          DO  lc = 1, nbins_aerosol * ncomponents_mass
             WRITE(salsa_char,'(i0)') lc
             CALL pmc_set_dataarray_name( 'parent', 'am_' // TRIM( salsa_char ),                    &
                                          'child',  'am_' // TRIM( salsa_char ), ierr )
          ENDDO
          IF ( .NOT. salsa_gases_from_chem )  THEN
             DO  lg = 1, ngases_salsa
                WRITE(salsa_char,'(i0)') lg
                CALL pmc_set_dataarray_name( 'parent', 'sg_' // TRIM( salsa_char ),                 &
                                             'child',  'sg_' // TRIM( salsa_char ), ierr )
             ENDDO
          ENDIF
       ENDIF

       CALL pmc_set_dataarray_name( lastentry = .TRUE. )
!
!--    Send grid to parent
       child_grid_dim(1)  = nx
       child_grid_dim(2)  = ny
       child_grid_dim(3)  = nz
       child_grid_info(1) = zw(nzt+1)
       child_grid_info(2) = zw(nzt)
       child_grid_info(3) = dx
       child_grid_info(4) = dy
       child_grid_info(5) = dz(1)

       IF ( myid == 0 )  THEN

          CALL pmc_send_to_parent( child_grid_dim, SIZE( child_grid_dim ), 0, 123, ierr )
          CALL pmc_send_to_parent( child_grid_info, SIZE( child_grid_info ), 0, 124, ierr )
          CALL pmc_send_to_parent( coord_x, nx + 1 + 2 * nbgp, 0, 11, ierr )
          CALL pmc_send_to_parent( coord_y, ny + 1 + 2 * nbgp, 0, 12, ierr )

          CALL pmc_recv_from_parent( rans_mode_parent, 1, 0, 19, ierr )
!
!--       Receive parent-grid information.
          CALL pmc_recv_from_parent( parent_grid_info_real,                    &
                                     SIZE(parent_grid_info_real), 0, 21, ierr )
          CALL pmc_recv_from_parent( parent_grid_info_int,  3, 0, 22, ierr )

       ENDIF

       CALL MPI_BCAST( parent_grid_info_real, SIZE(parent_grid_info_real),     &
                       MPI_REAL, 0, comm2d, ierr )
       CALL MPI_BCAST( parent_grid_info_int, 3, MPI_INTEGER, 0, comm2d, ierr )

       pg%dx = parent_grid_info_real(3)
       pg%dy = parent_grid_info_real(4)
       pg%dz = parent_grid_info_real(7)
       pg%nx = parent_grid_info_int(1)
       pg%ny = parent_grid_info_int(2)
       pg%nz = parent_grid_info_int(3)
!
!--    Allocate 1-D arrays for parent-grid coordinates and grid-spacings in the z-direction
       ALLOCATE( pg%coord_x(-nbgp:pg%nx+nbgp) )
       ALLOCATE( pg%coord_y(-nbgp:pg%ny+nbgp) )
       ALLOCATE( pg%dzu(1:pg%nz+1) )
       ALLOCATE( pg%dzw(1:pg%nz+1) )
       ALLOCATE( pg%zu(0:pg%nz+1) )
       ALLOCATE( pg%zw(0:pg%nz+1) )
!
!--    Get parent-grid coordinates and grid-spacings in the z-direction from the parent
       IF ( myid == 0)  THEN
          CALL pmc_recv_from_parent( pg%coord_x, pg%nx+1+2*nbgp, 0, 24, ierr )
          CALL pmc_recv_from_parent( pg%coord_y, pg%ny+1+2*nbgp, 0, 25, ierr )
          CALL pmc_recv_from_parent( pg%dzu, pg%nz+1, 0, 26, ierr )
          CALL pmc_recv_from_parent( pg%dzw, pg%nz+1, 0, 27, ierr )
          CALL pmc_recv_from_parent( pg%zu, pg%nz+2, 0, 28, ierr )
          CALL pmc_recv_from_parent( pg%zw, pg%nz+2, 0, 29, ierr )
       ENDIF
!
!--    Broadcast this information
       CALL MPI_BCAST( pg%coord_x, pg%nx+1+2*nbgp, MPI_REAL, 0, comm2d, ierr )
       CALL MPI_BCAST( pg%coord_y, pg%ny+1+2*nbgp, MPI_REAL, 0, comm2d, ierr )
       CALL MPI_BCAST( pg%dzu, pg%nz+1, MPI_REAL, 0, comm2d, ierr )
       CALL MPI_BCAST( pg%dzw, pg%nz+1, MPI_REAL, 0, comm2d, ierr )
       CALL MPI_BCAST( pg%zu, pg%nz+2,  MPI_REAL, 0, comm2d, ierr )
       CALL MPI_BCAST( pg%zw, pg%nz+2,  MPI_REAL, 0, comm2d, ierr )
       CALL MPI_BCAST( rans_mode_parent, 1, MPI_LOGICAL, 0, comm2d, ierr )       
!
!--    Find the index bounds for the nest domain in the parent-grid index space
       CALL pmci_map_child_grid_to_parent_grid
!
!--    TO_DO: Klaus give a comment what is happening here
       CALL pmc_c_get_2d_index_list
!
!--    Include couple arrays into child content
!--    TO_DO: Klaus: better explain the above comment (what is child content?)
       CALL  pmc_c_clear_next_array_list

       n  = 1
       lb = 1
       lc = 1
       lg = 1

       DO  WHILE ( pmc_c_getnextarray( myname ) )
!
!--       Note that pg%nz is not the original nz of parent, but the highest
!--       parent-grid level needed for nesting.
!--       Note that in case of chemical species or SALSA variables an additional 
!--       parameter needs to be passed. The parameter is required to set the pointer 
!--       correctlyto the chemical-species or SALSA data structure. Hence, first check if 
!--       the current variable is a chemical species or a SALSA variable. If so, pass 
!--       index id of respective sub-variable (species or bin) and increment this subsequently.
          IF ( INDEX( TRIM( myname ), 'chem_' ) /= 0 )  THEN             
             CALL pmci_create_childs_parent_grid_arrays ( myname, ipl, ipr, jps, jpn, pg%nz, n )
             n = n + 1    
          ELSEIF ( INDEX( TRIM( myname ), 'an_' ) /= 0 )  THEN
             CALL pmci_create_childs_parent_grid_arrays ( myname, ipl, ipr, jps, jpn, pg%nz, lb )
             lb = lb + 1
          ELSEIF ( INDEX( TRIM( myname ), 'am_' ) /= 0 )  THEN
             CALL pmci_create_childs_parent_grid_arrays ( myname, ipl, ipr, jps, jpn, pg%nz, lc )
             lc = lc + 1
          ELSEIF ( INDEX( TRIM( myname ), 'sg_' ) /= 0  .AND.  .NOT.  salsa_gases_from_chem )  THEN
             CALL pmci_create_childs_parent_grid_arrays ( myname, ipl, ipr, jps, jpn, pg%nz, lg )
             lg = lg + 1
          ELSE
             CALL pmci_create_childs_parent_grid_arrays ( myname, ipl, ipr, jps, jpn, pg%nz )
          ENDIF
       ENDDO
       CALL pmc_c_setind_and_allocmem
!
!--    Precompute the index-mapping arrays
       CALL pmci_define_index_mapping
!
!--    Check that the child and parent grid lines do match 
       CALL pmci_check_grid_matching
!       
!--    Compute surface areas of the nest-boundary faces 
       CALL pmci_compute_face_areas
       
    ENDIF

 CONTAINS


    SUBROUTINE pmci_map_child_grid_to_parent_grid
!
!--    Determine index bounds of interpolation/anterpolation area in the parent-grid index space
       IMPLICIT NONE

       INTEGER(iwp), DIMENSION(5,numprocs) ::  parent_bound_all     !< Transfer array for parent-grid index bounds

       INTEGER(iwp), DIMENSION(4)          ::  parent_bound_global  !< Transfer array for global parent-grid index bounds
       INTEGER(iwp), DIMENSION(2)          ::  size_of_array        !< For sending the dimensions of parent_bound_all to parent

       INTEGER(iwp) ::  ip      !< Running parent-grid index in the x-direction 
       INTEGER(iwp) ::  iauxl   !< Offset between the index bound ipl and the auxiliary index bound ipla
       INTEGER(iwp) ::  iauxr   !< Offset between the index bound ipr and the auxiliary index bound ipra
       INTEGER(iwp) ::  ijaux   !< Temporary variable for receiving the index bound from the neighbouring subdomain
       INTEGER(iwp) ::  jp      !< Running parent-grid index in the y-direction 
       INTEGER(iwp) ::  jauxs   !< Offset between the index bound jps and the auxiliary index bound jpsa
       INTEGER(iwp) ::  jauxn   !< Offset between the index bound jpn and the auxiliary index bound jpna

       REAL(wp) ::  tolex       !< Tolerance for grid-line matching in x-direction    
       REAL(wp) ::  toley       !< Tolerance for grid-line matching in y-direction    
       REAL(wp) ::  xexl        !< Parent-grid array exceedance behind the left edge of the child PE subdomain
       REAL(wp) ::  xexr        !< Parent-grid array exceedance behind the right edge of the child PE subdomain 
       REAL(wp) ::  yexs        !< Parent-grid array exceedance behind the south edge of the child PE subdomain
       REAL(wp) ::  yexn        !< Parent-grid array exceedance behind the north edge of the child PE subdomain
       REAL(wp) ::  xpl         !< Requested left-edge x-coordinate of the parent-grid array domain (at the internal boundaries
                                !< the real edge may differ from this in some cases as explained in the comment block below)  
       REAL(wp) ::  xpr         !< Requested right-edge x-coordinate of the parent-grid array domain (at the internal boundaries
                                !< the real edge may differ from this in some cases as explained in the comment block below)
       REAL(wp) ::  yps         !< Requested south-edge y-coordinate of the parent-grid array domain (at the internal boundaries
                                !< the real edge may differ from this in some cases as explained in the comment block below)
       REAL(wp) ::  ypn         !< Requested south-edge y-coordinate of the parent-grid array domain (at the internal boundaries
                                !< the real edge may differ from this in some cases as explained in the comment block below)

!
!--    Determine the index limits for the child's parent-grid arrays (such as uc for example).
!--    Note that at the outer edges of the child domain (nest boundaries) these arrays exceed
!--    the boundary by two parent-grid cells. At the internal boundaries, there are no 
!--    exceedances and thus no overlaps with the neighbouring subdomain. If at least half 
!--    of the parent-grid cell is within the current child sub-domain, then it is included 
!--    in the current sub-domain's parent-grid array. Else the parent-grid cell is 
!--    included in the neighbouring subdomain's parent-grid array, or not included at all if 
!--    we are at the outer edge of the child domain. This may occur especially when a large 
!--    grid-spacing ratio is used.
!
!--    Tolerances for grid-line matching.
       tolex = tolefac * dx
       toley = tolefac * dy
!
!--    Left boundary.
!--    Extension by two parent-grid cells behind the boundary, see the comment block above.    
       IF ( bc_dirichlet_l )  THEN
          xexl  = 2.0_wp * pg%dx
          iauxl = 0
       ELSE
          xexl  = 0.0_wp
          iauxl = 1
       ENDIF
       xpl     = coord_x(nxl) - xexl
       DO  ip = 0, pg%nx
          IF ( pg%coord_x(ip) + 0.5_wp * pg%dx >= xpl - tolex )  THEN
             ipl = MAX( 0, ip )
             EXIT
          ENDIF
       ENDDO
!
!--    Right boundary.
!--    Extension by two parent-grid cells behind the boundary, see the comment block above.       
       IF ( bc_dirichlet_r )  THEN
          xexr  = 2.0_wp * pg%dx
          iauxr = 0  
       ELSE
          xexr  = 0.0_wp
          iauxr = 1  
       ENDIF
       xpr  = coord_x(nxr+1) + xexr
       DO  ip = pg%nx, 0 , -1
          IF ( pg%coord_x(ip) + 0.5_wp * pg%dx <= xpr + tolex )  THEN
             ipr = MIN( pg%nx, MAX( ipl, ip ) )
             EXIT
          ENDIF
       ENDDO
!
!--    South boundary.
!--    Extension by two parent-grid cells behind the boundary, see the comment block above.   
       IF ( bc_dirichlet_s )  THEN
          yexs  = 2.0_wp * pg%dy
          jauxs = 0  
       ELSE
          yexs  = 0.0_wp
          jauxs = 1  
       ENDIF
       yps  = coord_y(nys) - yexs
       DO  jp = 0, pg%ny
          IF ( pg%coord_y(jp) + 0.5_wp * pg%dy >= yps - toley )  THEN             
             jps = MAX( 0, jp )
             EXIT
          ENDIF
       ENDDO
!
!--    North boundary.
!--    Extension by two parent-grid cells behind the boundary, see the comment block above.  
       IF  ( bc_dirichlet_n )  THEN
          yexn  = 2.0_wp * pg%dy
          jauxn = 0
       ELSE
          yexn  = 0.0_wp
          jauxn = 1
       ENDIF
       ypn  = coord_y(nyn+1) + yexn
       DO  jp = pg%ny, 0 , -1
          IF ( pg%coord_y(jp) + 0.5_wp * pg%dy <= ypn + toley )  THEN
             jpn = MIN( pg%ny, MAX( jps, jp ) )
             EXIT
          ENDIF
       ENDDO
!
!--    Make sure that the indexing is contiguous (no gaps, no overlaps). 
!--    This is a safety measure mainly for cases with high grid-spacing 
!--    ratio and narrow child subdomains.
       IF ( pdims(1) > 1 )  THEN
          IF ( nxl == 0 )  THEN
             CALL MPI_SEND( ipr, 1, MPI_INTEGER, pright, 717, comm2d, ierr )
          ELSE IF ( nxr == nx )  THEN
             CALL MPI_RECV( ijaux, 1, MPI_INTEGER, pleft, 717, comm2d, status, ierr )
             ipl = ijaux + 1
          ELSE
             CALL MPI_SEND( ipr, 1, MPI_INTEGER, pright, 717, comm2d, ierr )
             CALL MPI_RECV( ijaux, 1, MPI_INTEGER, pleft, 717, comm2d, status, ierr ) 
             ipl = ijaux + 1
          ENDIF
       ENDIF

       IF ( pdims(2) > 1 )  THEN
          IF ( nys == 0 )  THEN
             CALL MPI_SEND( jpn, 1, MPI_INTEGER, pnorth, 719, comm2d, ierr )
          ELSE IF ( nyn == ny )  THEN
             CALL MPI_RECV( ijaux, 1, MPI_INTEGER, psouth, 719, comm2d, status, ierr )
             jps = ijaux + 1
          ELSE
             CALL MPI_SEND( jpn, 1, MPI_INTEGER, pnorth, 719, comm2d, ierr )
             CALL MPI_RECV( ijaux, 1, MPI_INTEGER, psouth, 719, comm2d, status, ierr ) 
             jps = ijaux + 1
          ENDIF
       ENDIF
          
       WRITE(9,"('pmci_map_child_grid_to_parent_grid. Parent-grid array bounds: ',4(i4,2x))")             &
            ipl, ipr, jps, jpn
       FLUSH(9)

       parent_bound(1) = ipl
       parent_bound(2) = ipr
       parent_bound(3) = jps
       parent_bound(4) = jpn
       parent_bound(5) = myid
!
!--    The following auxiliary index bounds are used for allocating index mapping and 
!--    some other auxiliary arrays.
       ipla = ipl - iauxl
       ipra = ipr + iauxr
       jpsa = jps - jauxs
       jpna = jpn + jauxn
!
!--    The index-bounds parent_bound of all subdomains of the current child domain
!--    must be sent to the parent in order for the parent to create the index list.
!--    For this reason, the parent_bound arrays are packed together in single
!--    array parent_bound_all using MPI_GATHER.       
!--    Note that MPI_Gather receives data from all processes in the rank order
!--    This fact is exploited in creating the index list in pmci_create_index_list.
       CALL MPI_GATHER( parent_bound, 5, MPI_INTEGER, parent_bound_all, 5,                          &
                        MPI_INTEGER, 0, comm2d, ierr )

       IF ( myid == 0 )  THEN
          size_of_array(1) = SIZE( parent_bound_all, 1 )
          size_of_array(2) = SIZE( parent_bound_all, 2 )
          CALL pmc_send_to_parent( size_of_array, 2, 0, 40, ierr )
          CALL pmc_send_to_parent( parent_bound_all, SIZE( parent_bound_all ), 0, 41, ierr )
!
!--       Determine the global parent-grid index bounds       
          parent_bound_global(1) = MINVAL( parent_bound_all(1,:) )
          parent_bound_global(2) = MAXVAL( parent_bound_all(2,:) )
          parent_bound_global(3) = MINVAL( parent_bound_all(3,:) )
          parent_bound_global(4) = MAXVAL( parent_bound_all(4,:) )
       ENDIF
!
!--    Broadcast the global parent-grid index bounds to all current child processes
       CALL MPI_BCAST( parent_bound_global, 4, MPI_INTEGER, 0, comm2d, ierr )
       iplg = parent_bound_global(1)
       iprg = parent_bound_global(2)
       jpsg = parent_bound_global(3)
       jpng = parent_bound_global(4)
       WRITE( 9, "('pmci_map_child_grid_to_parent_grid. Global parent-grid index bounds iplg, iprg, jpsg, jpng: ',4(i4,2x))" ) &
            iplg, iprg, jpsg, jpng
       FLUSH( 9 )
       
    END SUBROUTINE pmci_map_child_grid_to_parent_grid

     
      
    SUBROUTINE pmci_define_index_mapping
!
!--    Precomputation of the mapping of the child- and parent-grid indices.

       IMPLICIT NONE

       INTEGER(iwp) ::  i         !< Child-grid index in the x-direction
       INTEGER(iwp) ::  ii        !< Parent-grid index in the x-direction
       INTEGER(iwp) ::  istart    !<
       INTEGER(iwp) ::  ir        !<
       INTEGER(iwp) ::  iw        !< Child-grid index limited to -1 <= iw <= nx+1 for wall_flags_total_0
       INTEGER(iwp) ::  j         !< Child-grid index in the y-direction
       INTEGER(iwp) ::  jj        !< Parent-grid index in the y-direction
       INTEGER(iwp) ::  jstart    !<
       INTEGER(iwp) ::  jr        !<
       INTEGER(iwp) ::  jw        !< Child-grid index limited to -1 <= jw <= ny+1 for wall_flags_total_0
       INTEGER(iwp) ::  k         !< Child-grid index in the z-direction
       INTEGER(iwp) ::  kk        !< Parent-grid index in the z-direction
       INTEGER(iwp) ::  kstart    !<
       INTEGER(iwp) ::  kw        !< Child-grid index limited to kw <= nzt+1 for wall_flags_total_0

       REAL(wp)     ::  tolex     !< Tolerance for grid-line matching in x-direction    
       REAL(wp)     ::  toley     !< Tolerance for grid-line matching in y-direction    
       REAL(wp)     ::  tolez     !< Tolerance for grid-line matching in z-direction    

!
!--    Grid-line tolerances.
       tolex = tolefac * dx
       toley = tolefac * dy
       tolez = tolefac * dz(1)
!
!--    Allocate child-grid work arrays for interpolation.
       igsr = NINT( pg%dx / dx, iwp )
       jgsr = NINT( pg%dy / dy, iwp )
       kgsr = NINT( pg%dzw(1) / dzw(1), iwp )
       WRITE(9,"('igsr, jgsr, kgsr: ',3(i3,2x))") igsr, jgsr, kgsr
       FLUSH(9)
!       
!--    Determine index bounds for the parent-grid work arrays for
!--    interpolation and allocate them.
       CALL pmci_allocate_workarrays
!       
!--    Define the MPI-datatypes for parent-grid work array 
!--    exchange between the PE-subdomains.
       CALL pmci_create_workarray_exchange_datatypes
!
!--    First determine kcto and kctw which refer to the uppermost 
!--    parent-grid levels below the child top-boundary level.
!--    Note that these comparison tests are not round-off-error
!--    sensitive and therefore tolerance buffering is not needed here.
       kk = 0
       DO WHILE ( pg%zu(kk) <= zu(nzt) )
          kk = kk + 1
       ENDDO
       kcto = kk - 1

       kk = 0
       DO WHILE ( pg%zw(kk) <= zw(nzt-1) )
          kk = kk + 1
       ENDDO
       kctw = kk - 1

       WRITE( 9, "('kcto, kctw = ', 2(i3,2x))" ) kcto, kctw
       FLUSH( 9 )
!       
!--    In case of two-way coupling, check that the child domain is sufficiently 
!--    large in terms of the number of parent-grid cells covered. Otherwise
!--    anterpolation is not possible.
       IF ( nesting_mode == 'two-way')  THEN
          CALL pmci_check_child_domain_size
       ENDIF
       
       ALLOCATE( iflu(ipla:ipra) )
       ALLOCATE( iflo(ipla:ipra) )
       ALLOCATE( ifuu(ipla:ipra) )
       ALLOCATE( ifuo(ipla:ipra) )
       ALLOCATE( jflv(jpsa:jpna) )
       ALLOCATE( jflo(jpsa:jpna) )
       ALLOCATE( jfuv(jpsa:jpna) )
       ALLOCATE( jfuo(jpsa:jpna) )       
       ALLOCATE( kflw(0:pg%nz+1) )
       ALLOCATE( kflo(0:pg%nz+1) )
       ALLOCATE( kfuw(0:pg%nz+1) )
       ALLOCATE( kfuo(0:pg%nz+1) )
       ALLOCATE( ijkfc_u(0:pg%nz+1,jpsa:jpna,ipla:ipra) )
       ALLOCATE( ijkfc_v(0:pg%nz+1,jpsa:jpna,ipla:ipra) )
       ALLOCATE( ijkfc_w(0:pg%nz+1,jpsa:jpna,ipla:ipra) )
       ALLOCATE( ijkfc_s(0:pg%nz+1,jpsa:jpna,ipla:ipra) )

       ijkfc_u = 0
       ijkfc_v = 0
       ijkfc_w = 0
       ijkfc_s = 0
!
!--    i-indices of u for each ii-index value
       istart = nxlg
       DO  ii = ipla, ipra
!
!--       The parent and child grid lines do always match in x, hence we 
!--       use only the local k,j-child-grid plane for the anterpolation.
!--       However, icru still has to be stored separately as these index bounds
!--       are passed as arguments to the interpolation and anterpolation 
!--       subroutines.
!--       Note that this comparison test is round-off-error sensitive 
!--       and therefore tolerance buffering is needed here.
          i = istart
          DO WHILE ( pg%coord_x(ii) - coord_x(i) > tolex  .AND. i < nxrg )
             i = i + 1
          ENDDO
          iflu(ii) = MIN( MAX( i, nxlg ), nxrg )
          ifuu(ii) = iflu(ii)
          istart   = iflu(ii)
!
!--       Print out the index bounds for checking and debugging purposes
          WRITE( 9, "('pmci_define_index_mapping, ii, iflu, ifuu: ', 3(i4,2x))" )                   &
               ii, iflu(ii), ifuu(ii)
          FLUSH( 9 )
       ENDDO
       WRITE( 9, * )
!
!--    i-indices of others for each ii-index value.
!--    Note that these comparison tests are not round-off-error
!--    sensitive and therefore tolerance buffering is not needed here.
       istart = nxlg
       DO  ii = ipla, ipra
          i = istart
          DO WHILE ( ( coord_x(i) + 0.5_wp * dx < pg%coord_x(ii) )  .AND.  ( i < nxrg ) )
             i  = i + 1
          ENDDO
          iflo(ii) = MIN( MAX( i, nxlg ), nxrg )
          ir = i
          DO WHILE ( ( coord_x(ir) + 0.5_wp * dx < pg%coord_x(ii) + pg%dx )  .AND.  ( i < nxrg+1 ) )
             i  = i + 1
             ir = MIN( i, nxrg )
          ENDDO
          ifuo(ii) = MIN( MAX( i-1, iflo(ii) ), nxrg )
          istart = iflo(ii)
!
!--       Print out the index bounds for checking and debugging purposes
          WRITE( 9, "('pmci_define_index_mapping, ii, iflo, ifuo: ', 3(i4,2x))" )                   &
               ii, iflo(ii), ifuo(ii)
          FLUSH( 9 )
       ENDDO
       WRITE( 9, * )
!
!--    j-indices of v for each jj-index value
       jstart = nysg
       DO  jj = jpsa, jpna
!
!--       The parent and child grid lines do always match in y, hence we 
!--       use only the local k,i-child-grid plane for the anterpolation.
!--       However, jcnv still has to be stored separately as these index bounds 
!--       are passed as arguments to the interpolation and anterpolation 
!--       subroutines.
!--       Note that this comparison test is round-off-error sensitive 
!--       and therefore tolerance buffering is needed here.
          j = jstart
          DO WHILE ( pg%coord_y(jj) - coord_y(j) > toley  .AND. j < nyng )
             j = j + 1
          ENDDO
          jflv(jj) = MIN( MAX( j, nysg ), nyng )
          jfuv(jj) = jflv(jj)
          jstart   = jflv(jj)
!
!--       Print out the index bounds for checking and debugging purposes
          WRITE( 9, "('pmci_define_index_mapping, jj, jflv, jfuv: ', 3(i4,2x))" )                   &
               jj, jflv(jj), jfuv(jj)
          FLUSH(9)
       ENDDO
       WRITE( 9, * )
!
!--    j-indices of others for each jj-index value
!--    Note that these comparison tests are not round-off-error
!--    sensitive and therefore tolerance buffering is not needed here.
       jstart = nysg
       DO  jj = jpsa, jpna
          j = jstart
          DO WHILE ( ( coord_y(j) + 0.5_wp * dy < pg%coord_y(jj) ) .AND. ( j < nyng ) )
             j  = j + 1
          ENDDO
          jflo(jj) = MIN( MAX( j, nysg ), nyng )
          jr = j
          DO WHILE ( ( coord_y(jr) + 0.5_wp * dy < pg%coord_y(jj) + pg%dy ) .AND. ( j < nyng+1 ) )
             j  = j + 1
             jr = MIN( j, nyng )
          ENDDO
          jfuo(jj) = MIN( MAX( j-1, jflo(jj) ), nyng )
          jstart = jflo(jj)
!
!--       Print out the index bounds for checking and debugging purposes
          WRITE( 9, "('pmci_define_index_mapping, jj, jflo, jfuo: ', 3(i4,2x))" )                   &
               jj, jflo(jj), jfuo(jj)
          FLUSH( 9 )
       ENDDO
       WRITE( 9, * )
!
!--    k-indices of w for each kk-index value
!--    Note that anterpolation index limits are needed also for the top boundary 
!--    ghost cell level because they are used also in the interpolation.
       kstart  = 0
       kflw(0) = 0
       kfuw(0) = 0
       DO  kk = 1, pg%nz+1
!
!--       The parent and child grid lines do always match in z, hence we 
!--       use only the local j,i-child-grid plane for the anterpolation.
!--       However, kctw still has to be stored separately as these index bounds 
!--       are passed as arguments to the interpolation and anterpolation 
!--       subroutines.
!--       Note that this comparison test is round-off-error sensitive 
!--       and therefore tolerance buffering is needed here.
          k = kstart
          DO WHILE ( ( pg%zw(kk) - zw(k) > tolez )  .AND.  ( k < nzt+1 ) )
             k = k + 1
          ENDDO
          kflw(kk) = MIN( MAX( k, 1 ), nzt + 1 )
          kfuw(kk) = kflw(kk)
          kstart   = kflw(kk)
!
!--       Print out the index bounds for checking and debugging purposes
          WRITE( 9, "('pmci_define_index_mapping, kk, kflw, kfuw: ', 4(i4,2x), 2(e12.5,2x))" )      &
               kk, kflw(kk), kfuw(kk), nzt,  pg%zu(kk), pg%zw(kk)
          FLUSH( 9 )
       ENDDO
       WRITE( 9, * )
!
!--    k-indices of others for each kk-index value
       kstart  = 0
       kflo(0) = 0
       kfuo(0) = 0
!
!--    Note that anterpolation index limits are needed also for the top boundary 
!--    ghost cell level because they are used also in the interpolation.
!--    Note that these comparison tests are not round-off-error
!--    sensitive and therefore tolerance buffering is not needed here.
       DO  kk = 1, pg%nz+1
          k = kstart
          DO WHILE ( ( zu(k) < pg%zw(kk-1) )  .AND.  ( k <= nzt ) )
             k = k + 1
          ENDDO
          kflo(kk) = MIN( MAX( k, 1 ), nzt + 1 )
          DO WHILE ( ( zu(k) < pg%zw(kk) )  .AND.  ( k <= nzt+1 ) )
             k = k + 1
             IF ( k > nzt + 1 ) EXIT  ! This EXIT is to prevent zu(k) from flowing over.
          ENDDO
          kfuo(kk) = MIN( MAX( k-1, kflo(kk) ), nzt + 1 )
          kstart = kflo(kk)
       ENDDO
!
!--    Print out the index bounds for checking and debugging purposes
       DO  kk = 1, pg%nz+1
          WRITE( 9, "('pmci_define_index_mapping, kk, kflo, kfuo: ', 4(i4,2x), 2(e12.5,2x))" )      &
               kk, kflo(kk), kfuo(kk), nzt,  pg%zu(kk), pg%zw(kk)
          FLUSH( 9 )
       ENDDO
       WRITE( 9, * )
!
!--    Precomputation of number of child-grid nodes inside parent-grid cells.
!--    Note that ii, jj, and kk are parent-grid indices. 
!--    This information is needed in the anterpolation. 
!--    The indices for wall_flags_total_0 (kw,jw,iw) must be limited to the range 
!--    [-1,...,nx/ny/nzt+1] in order to avoid zero values on the outer ghost nodes.
       DO  ii = ipla, ipra
          DO  jj = jpsa, jpna
             DO  kk = 0, pg%nz+1
!
!--             u-component
                DO  i = iflu(ii), ifuu(ii)
                   iw = MAX( MIN( i, nx+1 ), -1 )
                   DO  j = jflo(jj), jfuo(jj)
                      jw = MAX( MIN( j, ny+1 ), -1 )
                      DO  k = kflo(kk), kfuo(kk)
                         kw = MIN( k, nzt+1 )                
                         ijkfc_u(kk,jj,ii) = ijkfc_u(kk,jj,ii)                                      &
                              + MERGE( 1, 0, BTEST( wall_flags_total_0(kw,jw,iw), 1 ) )
                      ENDDO
                   ENDDO
                ENDDO
!
!--             v-component 
                DO  i = iflo(ii), ifuo(ii)
                   iw = MAX( MIN( i, nx+1 ), -1 )
                   DO  j = jflv(jj), jfuv(jj)
                      jw = MAX( MIN( j, ny+1 ), -1 )
                      DO  k = kflo(kk), kfuo(kk)
                         kw = MIN( k, nzt+1 )                                        
                         ijkfc_v(kk,jj,ii) = ijkfc_v(kk,jj,ii)                                      &
                              + MERGE( 1, 0, BTEST( wall_flags_total_0(kw,jw,iw), 2 ) )
                      ENDDO
                   ENDDO
                ENDDO
!
!--             scalars
                DO  i = iflo(ii), ifuo(ii)
                   iw = MAX( MIN( i, nx+1 ), -1 )
                   DO  j = jflo(jj), jfuo(jj)
                      jw = MAX( MIN( j, ny+1 ), -1 )
                      DO  k = kflo(kk), kfuo(kk)
                         kw = MIN( k, nzt+1 )
                         ijkfc_s(kk,jj,ii) = ijkfc_s(kk,jj,ii)                                      &
                              + MERGE( 1, 0, BTEST( wall_flags_total_0(kw,jw,iw), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
!
!--             w-component 
                DO  i = iflo(ii), ifuo(ii)
                   iw = MAX( MIN( i, nx+1 ), -1 )
                   DO  j = jflo(jj), jfuo(jj)
                      jw = MAX( MIN( j, ny+1 ), -1 )
                      DO  k = kflw(kk), kfuw(kk)
                         kw = MIN( k, nzt+1 )
                         ijkfc_w(kk,jj,ii) = ijkfc_w(kk,jj,ii)                                      &
                              + MERGE( 1, 0, BTEST( wall_flags_total_0(kw,jw,iw), 3 ) )
                      ENDDO
                   ENDDO
                ENDDO

             ENDDO  ! kk        
          ENDDO  ! jj
       ENDDO  ! ii

    END SUBROUTINE pmci_define_index_mapping



    SUBROUTINE pmci_check_child_domain_size
!       
!--    Check if the child domain is too small in terms of number of parent-grid cells 
!--    covered so that anterpolation buffers fill the whole domain so that anterpolation
!--    not possible. Also, check that anterpolation_buffer_width is not too large to  
!--    prevent anterpolation.
       IMPLICIT NONE
      
!
!--    First x-direction
       IF ( iplg + 3 + anterpolation_buffer_width > iprg - 3 - anterpolation_buffer_width )  THEN
          IF ( iprg - iplg + 1 < 7 )  THEN
!
!--          Error
             WRITE( message_string, * ) 'child domain too narrow for anterpolation in x-direction'
             CALL message( 'pmci_check_child_domain_size', 'PA0652', 3, 2, 0, 6, 0 )
          ELSE IF ( iprg - iplg + 1 < 11 )  THEN
!                
!--          Warning
             WRITE( message_string, * ) 'anterpolation_buffer_width value too high, reset to 0'
             CALL message( 'pmci_check_child_domain_size', 'PA0653', 0, 1, 0, 6, 0 )
             anterpolation_buffer_width = 0
          ELSE
!                
!--          Informative message
             WRITE( message_string, * ) 'anterpolation_buffer_width value too high, reset to default value 2'
             CALL message( 'pmci_check_child_domain_size', 'PA0654', 0, 0, 0, 6, 0 )
             anterpolation_buffer_width = 2
          ENDIF
       ENDIF
!
!--    Then y-direction          
       IF ( jpsg + 3 + anterpolation_buffer_width > jpng - 3 - anterpolation_buffer_width )  THEN
          IF ( jpng - jpsg + 1 < 7 )  THEN
!
!--          Error
             WRITE( message_string, * ) 'child domain too narrow for anterpolation in y-direction'
             CALL message( 'pmci_check_child_domain_size', 'PA0652', 3, 2, 0, 6, 0 )
          ELSE IF ( jpng - jpsg + 1 < 11 )  THEN
!                
!--          Warning
             WRITE( message_string, * ) 'anterpolation_buffer_width value too high, reset to 0'
             CALL message( 'pmci_check_child_domain_size', 'PA0653', 0, 1, 0, 6, 0 )
             anterpolation_buffer_width = 0
          ELSE
!                
!--          Informative message
             WRITE( message_string, * ) 'anterpolation_buffer_width value too high, reset to default value 2'
             CALL message( 'pmci_check_child_domain_size', 'PA0654', 0, 0, 0, 6, 0 )
             anterpolation_buffer_width = 2
          ENDIF
       ENDIF
!
!--    Finally z-direction               
       IF ( kctw - 1 - anterpolation_buffer_width < 1 )  THEN
          IF ( kctw - 1 < 1 )  THEN
!
!--          Error
             WRITE( message_string, * ) 'child domain too shallow for anterpolation in z-direction'
             CALL message( 'pmci_check_child_domain_size', 'PA0652', 3, 2, 0, 6, 0 )
          ELSE IF ( kctw - 3 < 1 )  THEN
!                
!--          Warning
             WRITE( message_string, * ) 'anterpolation_buffer_width value too high, reset to 0'
             CALL message( 'pmci_check_child_domain_size', 'PA0653', 0, 1, 0, 6, 0 )
             anterpolation_buffer_width = 0
          ELSE
!                
!--          Informative message
             WRITE( message_string, * ) 'anterpolation_buffer_width value too high, reset to default value 2'
             CALL message( 'pmci_check_child_domain_size', 'PA0654', 0, 0, 0, 6, 0 )
             anterpolation_buffer_width = 2  
          ENDIF
       ENDIF

    END SUBROUTINE pmci_check_child_domain_size


    
    SUBROUTINE pmci_allocate_workarrays
!
!--    Allocate parent-grid work-arrays for interpolation
       IMPLICIT NONE

!
!--    Determine and store the PE-subdomain dependent index bounds
       IF ( bc_dirichlet_l )  THEN
          iplw = ipl + 1
       ELSE
          iplw = ipl - 1
       ENDIF

       IF ( bc_dirichlet_r )  THEN
          iprw = ipr - 1
       ELSE
          iprw = ipr + 1
       ENDIF

       IF ( bc_dirichlet_s )  THEN
          jpsw = jps + 1
       ELSE
          jpsw = jps - 1
       ENDIF

       IF ( bc_dirichlet_n )  THEN
          jpnw = jpn - 1
       ELSE
          jpnw = jpn + 1
       ENDIF
!
!--    Left and right boundaries.
       ALLOCATE( workarr_lr(0:pg%nz+1,jpsw:jpnw,0:2) )
!
!--    South and north boundaries.
       ALLOCATE( workarr_sn(0:pg%nz+1,0:2,iplw:iprw) )
!
!--    Top boundary.
       ALLOCATE( workarr_t(0:2,jpsw:jpnw,iplw:iprw) )

    END SUBROUTINE pmci_allocate_workarrays



    SUBROUTINE pmci_create_workarray_exchange_datatypes
!
!--    Define specific MPI types for workarr-exchange.
       IMPLICIT NONE

!
!--    For the left and right boundaries
       CALL MPI_TYPE_VECTOR( 3, pg%nz+2, (jpnw-jpsw+1)*(pg%nz+2), MPI_REAL,                         &
            workarr_lr_exchange_type, ierr )
       CALL MPI_TYPE_COMMIT( workarr_lr_exchange_type, ierr )
!
!--    For the south and north boundaries
       CALL MPI_TYPE_VECTOR( 1, 3*(pg%nz+2), 3*(pg%nz+2), MPI_REAL,                                 &
            workarr_sn_exchange_type, ierr )
       CALL MPI_TYPE_COMMIT( workarr_sn_exchange_type, ierr )
!
!--    For the top-boundary x-slices
       CALL MPI_TYPE_VECTOR( iprw-iplw+1, 3, 3*(jpnw-jpsw+1), MPI_REAL,                             &
            workarr_t_exchange_type_x, ierr )
       CALL MPI_TYPE_COMMIT( workarr_t_exchange_type_x, ierr )
!
!--    For the top-boundary y-slices
       CALL MPI_TYPE_VECTOR( 1, 3*(jpnw-jpsw+1), 3*(jpnw-jpsw+1), MPI_REAL,                         &
            workarr_t_exchange_type_y, ierr )
       CALL MPI_TYPE_COMMIT( workarr_t_exchange_type_y, ierr )
       
    END SUBROUTINE pmci_create_workarray_exchange_datatypes



    SUBROUTINE pmci_check_grid_matching 
!
!--    Check that the grid lines of child and parent do match.
!--    Also check that the child subdomain width is not smaller than 
!--    the parent grid spacing in the respective direction. 
       IMPLICIT NONE

       INTEGER(iwp) ::  non_matching_height = 0              !< Flag for non-matching child-domain height
       INTEGER(iwp) ::  non_matching_lower_left_corner = 0   !< Flag for non-matching lower left corner
       INTEGER(iwp) ::  non_matching_upper_right_corner = 0  !< Flag for non-matching upper right corner
       INTEGER(iwp) ::  non_int_gsr_x = 0                    !< Flag for non-integer grid-spacing ration in x-direction
       INTEGER(iwp) ::  non_int_gsr_y = 0                    !< Flag for non-integer grid-spacing ration in y-direction
       INTEGER(iwp) ::  non_int_gsr_z = 0                    !< Flag for non-integer grid-spacing ration in z-direction
       INTEGER(iwp) ::  too_narrow_pesd_x = 0                !< Flag for too narrow pe-subdomain in x-direction
       INTEGER(iwp) ::  too_narrow_pesd_y = 0                !< Flag for too narrow pe-subdomain in y-direction
                                                                                                                  
       REAL(wp) ::  child_ngp_x_l                            !< Number of gridpoints in child subdomain in x-direction
                                                             !< converted to REAL(wp) 
       REAL(wp) ::  child_ngp_y_l                            !< Number of gridpoints in child subdomain in y-direction
                                                             !< converted to REAL(wp)
       REAL(wp) ::  gridline_mismatch_x                      !< Mismatch between the parent and child gridlines in the x-direction
       REAL(wp) ::  gridline_mismatch_y                      !< Mismatch between the parent and child gridlines in the y-direction
       REAL(wp) ::  gsr_mismatch_x                           !< Deviation of the grid-spacing ratio from the nearest integer value, the x-direction
       REAL(wp) ::  gsr_mismatch_y                           !< Deviation of the grid-spacing ratio from the nearest integer value, the y-direction
       REAL(wp) ::  upper_right_coord_x                      !< X-coordinate of the upper right corner of the child domain
       REAL(wp) ::  upper_right_coord_y                      !< Y-coordinate of the upper right corner of the child domain
       REAL(wp) ::  tolex                                    !< Tolerance for grid-line matching in x-direction
       REAL(wp) ::  toley                                    !< Tolerance for grid-line matching in y-direction
       REAL(wp) ::  tolez                                    !< Tolerance for grid-line matching in z-direction

       
       IF ( myid == 0 )  THEN

          tolex = tolefac * dx
          toley = tolefac * dy
          tolez = tolefac * dz(1)
!
!--       First check that the child domain lower left corner matches the parent grid lines.
          gridline_mismatch_x = ABS( NINT( lower_left_coord_x / pg%dx ) * pg%dx - lower_left_coord_x )
          gridline_mismatch_y = ABS( NINT( lower_left_coord_y / pg%dy ) * pg%dy - lower_left_coord_y )
          IF ( gridline_mismatch_x > tolex ) non_matching_lower_left_corner = 1
          IF ( gridline_mismatch_y > toley ) non_matching_lower_left_corner = 1
!
!--       Then check that the child doman upper right corner matches the parent grid lines.
          upper_right_coord_x = lower_left_coord_x + ( nx + 1 ) * dx
          upper_right_coord_y = lower_left_coord_y + ( ny + 1 ) * dy
          gridline_mismatch_x = ABS( NINT( upper_right_coord_x / pg%dx ) * pg%dx - upper_right_coord_x )
          gridline_mismatch_y = ABS( NINT( upper_right_coord_y / pg%dy ) * pg%dy - upper_right_coord_y )
          IF ( gridline_mismatch_x > tolex ) non_matching_upper_right_corner = 1
          IF ( gridline_mismatch_y > toley ) non_matching_upper_right_corner = 1
!
!--       Also check that the cild domain height matches the parent grid lines.
          IF ( MOD( zw(nzt), pg%dz ) > tolez ) non_matching_height = 1
!
!--       Check that the grid-spacing ratios in each direction are integer valued.    
          gsr_mismatch_x = ABS( NINT( pg%dx / dx ) * dx - pg%dx )
          gsr_mismatch_y = ABS( NINT( pg%dy / dy ) * dy - pg%dy )
          IF ( gsr_mismatch_x > tolex )  non_int_gsr_x = 1
          IF ( gsr_mismatch_y > toley )  non_int_gsr_y = 1
!
!--       In the z-direction, all levels need to be checked separately against grid stretching  
!--       which is not allowed.
          DO  n = 0, kctw+1
             IF ( ABS( pg%zw(n) - zw(kflw(n)) ) > tolez )  non_int_gsr_z = 1
          ENDDO

          child_ngp_x_l = REAL( nxr - nxl + 1, KIND=wp )
          IF ( child_ngp_x_l / REAL( igsr, KIND=wp ) < 1.0_wp )  too_narrow_pesd_x = 1
          child_ngp_y_l = REAL( nyn - nys + 1, KIND=wp )
          IF ( child_ngp_y_l / REAL( jgsr, KIND=wp ) < 1.0_wp )  too_narrow_pesd_y = 1
         
          IF ( non_matching_height > 0 )  THEN
             WRITE( message_string, * ) 'nested child domain height must match ',                   &
                                        'its parent grid lines'
             CALL message( 'pmci_check_grid_matching', 'PA0414', 3, 2, 0, 6, 0 )
          ENDIF

          IF ( non_matching_lower_left_corner > 0 )  THEN
             WRITE( message_string, * ) 'nested child domain lower left ',                          &
                                        'corner must match its parent grid lines'
             CALL message( 'pmci_check_grid_matching', 'PA0414', 3, 2, 0, 6, 0 )
          ENDIF

          IF ( non_matching_upper_right_corner > 0 )  THEN
             WRITE( message_string, * ) 'nested child domain upper right ',                         &
                                        'corner must match its parent grid lines'
             CALL message( 'pmci_check_grid_matching', 'PA0414', 3, 2, 0, 6, 0 )
          ENDIF

          IF ( non_int_gsr_x > 0 )  THEN
             WRITE( message_string, * ) 'nesting grid-spacing ratio ( parent dx / child dx ) ',     &
                                        'must have an integer value'
             CALL message( 'pmci_check_grid_matching', 'PA0416', 3, 2, 0, 6, 0 )
          ENDIF

          IF ( non_int_gsr_y > 0 )  THEN
             WRITE( message_string, * ) 'nesting grid-spacing ratio ( parent dy / child dy ) ',     &
                                        'must have an integer value'
             CALL message( 'pmci_check_grid_matching', 'PA0416', 3, 2, 0, 6, 0 )
          ENDIF

          IF ( non_int_gsr_z > 0 )  THEN
             WRITE( message_string, * ) 'nesting grid-spacing ratio ( parent dz / child dz ) ',     &
                                        'must have an integer value for each z-level'
             CALL message( 'pmci_check_grid_matching', 'PA0416', 3, 2, 0, 6, 0 )
          ENDIF

          IF ( too_narrow_pesd_x > 0 )  THEN
            WRITE( message_string, * ) 'child subdomain width in x-direction must not be ',        &
                                        'smaller than its parent grid dx. Change the PE-grid ',     &
                                        'setting (npex, npey) to satisfy this requirement.'  
             CALL message( 'pmci_check_grid_matching', 'PA0587', 3, 2, 0, 6, 0 )
          ENDIF
 
          IF ( too_narrow_pesd_y > 0 )  THEN
             WRITE( message_string, * ) 'child subdomain width in y-direction must not be ',        &
                                        'smaller than its parent grid dy. Change the PE-grid ',     &
                                        'setting (npex, npey) to satisfy this requirement.'  
             CALL message( 'pmci_check_grid_matching', 'PA0587', 3, 2, 0, 6, 0 )
          ENDIF
                 
       ENDIF  !  ( myid == 0 )
       
    END SUBROUTINE pmci_check_grid_matching



    SUBROUTINE pmci_compute_face_areas

       IMPLICIT NONE
       REAL(wp)  :: face_area_local   !< Local (for the current pe) integral face area of the left boundary
       REAL(wp)  :: sub_sum           !< Intermediate sum in order to improve the accuracy of the summation

       INTEGER(iwp) :: i              !< Running index in the x-direction
       INTEGER(iwp) :: j              !< Running index in the y-direction
       INTEGER(iwp) :: k              !< Running index in the z-direction 
       INTEGER(iwp) :: k_wall         !< Local topography top k-index
       INTEGER(iwp) :: n              !< Running index over boundary faces

       
!
!--    Sum up the volume flow through the left boundary
       face_area(1) = 0.0_wp
       face_area_local = 0.0_wp
       IF ( bc_dirichlet_l )  THEN
          i = 0
          DO  j = nys, nyn
             sub_sum = 0.0_wp
             k_wall = topo_top_ind(j,i,1)
             DO   k = k_wall + 1, nzt
                sub_sum = sub_sum + dzw(k)
             ENDDO
             face_area_local =  face_area_local + dy * sub_sum
          ENDDO
       ENDIF
       
#if defined( __parallel )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( face_area_local, face_area(1), 1, MPI_REAL, MPI_SUM, comm2d, ierr )
#else
       face_area(1) = face_area_local
#endif
!
!--    Sum up the volume flow through the right boundary
       face_area(2) = 0.0_wp
       face_area_local = 0.0_wp
       IF ( bc_dirichlet_r )  THEN
          i = nx
          DO  j = nys, nyn
             sub_sum = 0.0_wp
             k_wall = topo_top_ind(j,i,1)
             DO   k = k_wall + 1, nzt
                sub_sum = sub_sum + dzw(k)
             ENDDO
             face_area_local =  face_area_local + dy * sub_sum
          ENDDO
       ENDIF
       
#if defined( __parallel )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( face_area_local, face_area(2), 1, MPI_REAL, MPI_SUM, comm2d, ierr )
#else
       face_area(2) = face_area_local
#endif
!
!--    Sum up the volume flow through the south boundary
       face_area(3) = 0.0_wp
       face_area_local = 0.0_wp
       IF ( bc_dirichlet_s )  THEN
          j = 0
          DO  i = nxl, nxr
             sub_sum = 0.0_wp
             k_wall = topo_top_ind(j,i,2)
             DO  k = k_wall + 1, nzt
                sub_sum = sub_sum + dzw(k)
             ENDDO
             face_area_local = face_area_local + dx * sub_sum
          ENDDO
       ENDIF
       
#if defined( __parallel )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( face_area_local, face_area(3), 1, MPI_REAL, MPI_SUM, comm2d, ierr )
#else
       face_area(3) = face_area_local
#endif
!
!--    Sum up the volume flow through the north boundary
       face_area(4) = 0.0_wp
       face_area_local = 0.0_wp
       IF ( bc_dirichlet_n )  THEN
          j = ny
          DO  i = nxl, nxr
             sub_sum = 0.0_wp
             k_wall = topo_top_ind(j,i,2)
             DO  k = k_wall + 1, nzt
                sub_sum = sub_sum + dzw(k)
             ENDDO
             face_area_local = face_area_local + dx * sub_sum
          ENDDO
       ENDIF
       
#if defined( __parallel )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( face_area_local, face_area(4), 1, MPI_REAL, MPI_SUM, comm2d, ierr )
#else
       face_area(4) = face_area_local
#endif
!
!--    The top face area does not depend on the topography at all.       
       face_area(5) = ( nx + 1 ) * ( ny + 1 ) * dx * dy
!
!--    The 6th element is used for the total area
       face_area(6) = 0.0_wp
       DO  n = 1, 5
          face_area(6) = face_area(6) + face_area(n)
       ENDDO

!       write( 9, "(6(e12.5,2x))") ( face_area(n), n = 1, 6 )
!       flush( 9 )
       
    END SUBROUTINE pmci_compute_face_areas
#endif
    
 END SUBROUTINE pmci_setup_child



 SUBROUTINE pmci_setup_coordinates

#if defined( __parallel )
    IMPLICIT NONE

    INTEGER(iwp) ::  i   !<
    INTEGER(iwp) ::  j   !<

!
!-- Create coordinate arrays.
    ALLOCATE( coord_x(-nbgp:nx+nbgp) )
    ALLOCATE( coord_y(-nbgp:ny+nbgp) )
     
    DO  i = -nbgp, nx + nbgp
       coord_x(i) = lower_left_coord_x + i * dx
    ENDDO

    DO  j = -nbgp, ny + nbgp
       coord_y(j) = lower_left_coord_y + j * dy
    ENDDO

#endif
 END SUBROUTINE pmci_setup_coordinates

!------------------------------------------------------------------------------!
! Description:
! ------------
!> In this subroutine the number of coupled arrays is determined. 
!------------------------------------------------------------------------------! 
  SUBROUTINE pmci_num_arrays  
               
#if defined( __parallel ) 
    IMPLICIT NONE
!
!-- The number of coupled arrays depends on the model settings. At least
!-- 5 arrays need to be coupled (u, v, w, e, diss).  Please note, actually
!-- e and diss (TKE and dissipation rate) are only required if RANS-RANS
!-- nesting is applied, but memory is allocated nevertheless. This is because
!-- the information whether they are needed or not is retrieved at a later
!-- point in time. In case e and diss are not needed, they are also not 
!-- exchanged between parent and child. 
    pmc_max_array = 5
!
!-- pt
    IF ( .NOT. neutral )  pmc_max_array = pmc_max_array + 1
    
    IF ( humidity )  THEN
!
!--    q
       pmc_max_array = pmc_max_array + 1
!
!--    qc, nc
       IF ( bulk_cloud_model  .AND.  microphysics_morrison )                   &
          pmc_max_array = pmc_max_array + 2
!
!--    qr, nr
       IF ( bulk_cloud_model  .AND.  microphysics_seifert )                    &
          pmc_max_array = pmc_max_array + 2
    ENDIF
!
!-- s
    IF ( passive_scalar )  pmc_max_array = pmc_max_array + 1
!
!-- nr_part, part_adr
    IF ( particle_advection )  pmc_max_array = pmc_max_array + 2
!
!-- Chemistry, depends on number of species
    IF ( air_chemistry  .AND.  nesting_chem )  pmc_max_array = pmc_max_array + nspec
!
!-- SALSA, depens on the number aerosol size bins and chemical components + 
!-- the number of default gases
    IF ( salsa  .AND.  nesting_salsa )  pmc_max_array = pmc_max_array + nbins_aerosol +            &
                                                        nbins_aerosol * ncomponents_mass
    IF ( .NOT. salsa_gases_from_chem )  pmc_max_array = pmc_max_array + ngases_salsa

#endif
    
 END SUBROUTINE pmci_num_arrays


 SUBROUTINE pmci_set_array_pointer( name, child_id, nz_child, n )
   
    IMPLICIT NONE
   
    INTEGER(iwp), INTENT(IN) ::  child_id  !<
    INTEGER(iwp), INTENT(IN) ::  nz_child  !<
    
    INTEGER(iwp), INTENT(IN), OPTIONAL ::  n          !< index of chemical species
    
    CHARACTER(LEN=*), INTENT(IN) ::  name             !<

#if defined( __parallel )     
!
!-- Local variables:       
    INTEGER(iwp) ::  ierr                             !< MPI error code

    INTEGER(idp), POINTER, DIMENSION(:,:) ::  i_2d    !<
       
    REAL(wp), POINTER, DIMENSION(:,:)   ::  p_2d      !<
    REAL(wp), POINTER, DIMENSION(:,:,:) ::  p_3d      !<
    REAL(wp), POINTER, DIMENSION(:,:,:) ::  p_3d_sec  !<
    

    NULLIFY( p_3d )
    NULLIFY( p_2d )
    NULLIFY( i_2d )
!
!-- List of array names, which can be coupled.
!-- In case of 3D please change also the second array for the pointer version
    IF ( TRIM(name) == "u"          )  p_3d => u
    IF ( TRIM(name) == "v"          )  p_3d => v
    IF ( TRIM(name) == "w"          )  p_3d => w
    IF ( TRIM(name) == "e"          )  p_3d => e
    IF ( TRIM(name) == "pt"         )  p_3d => pt
    IF ( TRIM(name) == "q"          )  p_3d => q
    IF ( TRIM(name) == "qc"         )  p_3d => qc
    IF ( TRIM(name) == "qr"         )  p_3d => qr
    IF ( TRIM(name) == "nr"         )  p_3d => nr
    IF ( TRIM(name) == "nc"         )  p_3d => nc
    IF ( TRIM(name) == "s"          )  p_3d => s
    IF ( TRIM(name) == "diss"       )  p_3d => diss   
    IF ( TRIM(name) == "nr_part"    )  i_2d => nr_part
    IF ( TRIM(name) == "part_adr"   )  i_2d => part_adr
    IF ( INDEX( TRIM(name), "chem_" ) /= 0      )  p_3d => chem_species(n)%conc
    IF ( INDEX( TRIM(name), "an_" ) /= 0  )  p_3d => aerosol_number(n)%conc
    IF ( INDEX( TRIM(name), "am_" ) /= 0 )  p_3d => aerosol_mass(n)%conc
    IF ( INDEX( TRIM(name), "sg_" ) /= 0  .AND.  .NOT. salsa_gases_from_chem ) &
       p_3d => salsa_gas(n)%conc
!
!-- Next line is just an example for a 2D array (not active for coupling!) 
!-- Please note, that z0 has to be declared as TARGET array in modules.f90.
!    IF ( TRIM(name) == "z0" )    p_2d => z0
    IF ( TRIM(name) == "u"    )  p_3d_sec => u_2
    IF ( TRIM(name) == "v"    )  p_3d_sec => v_2
    IF ( TRIM(name) == "w"    )  p_3d_sec => w_2
    IF ( TRIM(name) == "e"    )  p_3d_sec => e_2
    IF ( TRIM(name) == "pt"   )  p_3d_sec => pt_2
    IF ( TRIM(name) == "q"    )  p_3d_sec => q_2
    IF ( TRIM(name) == "qc"   )  p_3d_sec => qc_2
    IF ( TRIM(name) == "qr"   )  p_3d_sec => qr_2
    IF ( TRIM(name) == "nr"   )  p_3d_sec => nr_2
    IF ( TRIM(name) == "nc"   )  p_3d_sec => nc_2
    IF ( TRIM(name) == "s"    )  p_3d_sec => s_2
    IF ( TRIM(name) == "diss" )  p_3d_sec => diss_2
    IF ( INDEX( TRIM(name), "chem_" ) /= 0 )  p_3d_sec => spec_conc_2(:,:,:,n)
    IF ( INDEX( TRIM(name), "an_" )   /= 0 )  p_3d_sec => nconc_2(:,:,:,n)
    IF ( INDEX( TRIM(name), "am_" )   /= 0 )  p_3d_sec => mconc_2(:,:,:,n)
    IF ( INDEX( TRIM(name), "sg_" )   /= 0  .AND.  .NOT. salsa_gases_from_chem ) &
                                 p_3d_sec => gconc_2(:,:,:,n)

    IF ( ASSOCIATED( p_3d ) )  THEN
       CALL pmc_s_set_dataarray( child_id, p_3d, nz_child, nz, array_2 = p_3d_sec )
    ELSEIF  ( ASSOCIATED( p_2d ) )  THEN
       CALL pmc_s_set_dataarray( child_id, p_2d )
    ELSEIF  ( ASSOCIATED( i_2d ) )  THEN
       CALL pmc_s_set_dataarray( child_id, i_2d )
    ELSE
!
!--    Give only one message for the root domain
       IF ( pmc_is_rootmodel()  .AND.  myid == 0 )  THEN
          message_string = 'pointer for array "' // TRIM( name ) // '" can''t be associated'
          CALL message( 'pmci_set_array_pointer', 'PA0117', 3, 2, 0, 6, 0 )
       ELSE
!
!--       Avoid others to continue
          CALL MPI_BARRIER( comm2d, ierr )
       ENDIF
       
    ENDIF
    
#endif
    
 END SUBROUTINE pmci_set_array_pointer


     
 INTEGER FUNCTION get_number_of_children()

    IMPLICIT NONE

    
#if defined( __parallel )
    get_number_of_children = SIZE( pmc_parent_for_child ) - 1
#else
    get_number_of_children = 0
#endif

    RETURN

 END FUNCTION get_number_of_children


 
 INTEGER FUNCTION get_childid( id_index )

    IMPLICIT NONE

    INTEGER, INTENT(IN) ::  id_index   !<

    
#if defined( __parallel )
    get_childid = pmc_parent_for_child(id_index)
#else
    get_childid = 0
#endif

    RETURN

 END FUNCTION get_childid



 SUBROUTINE get_child_edges( m, lx_coord, lx_coord_b, rx_coord, rx_coord_b, sy_coord, sy_coord_b,   &
      ny_coord, ny_coord_b, uz_coord, uz_coord_b )
   
    IMPLICIT NONE

    INTEGER,INTENT(IN)   ::  m                     !<

    REAL(wp),INTENT(OUT) ::  lx_coord, lx_coord_b  !<
    REAL(wp),INTENT(OUT) ::  rx_coord, rx_coord_b  !<
    REAL(wp),INTENT(OUT) ::  ny_coord, ny_coord_b  !<
    REAL(wp),INTENT(OUT) ::  sy_coord, sy_coord_b  !<
    REAL(wp),INTENT(OUT) ::  uz_coord, uz_coord_b  !<

    
#if defined( __parallel )
    
    lx_coord = childgrid(m)%lx_coord
    rx_coord = childgrid(m)%rx_coord
    sy_coord = childgrid(m)%sy_coord
    ny_coord = childgrid(m)%ny_coord
    uz_coord = childgrid(m)%uz_coord
    
    lx_coord_b = childgrid(m)%lx_coord_b
    rx_coord_b = childgrid(m)%rx_coord_b
    sy_coord_b = childgrid(m)%sy_coord_b
    ny_coord_b = childgrid(m)%ny_coord_b
    uz_coord_b = childgrid(m)%uz_coord_b
    
#endif
    
 END SUBROUTINE get_child_edges



 SUBROUTINE  get_child_gridspacing( m, dx, dy,dz )

    IMPLICIT NONE
   
    INTEGER, INTENT(IN)             ::  m      !<

    REAL(wp), INTENT(OUT)           ::  dx,dy  !<

    REAL(wp), INTENT(OUT), OPTIONAL ::  dz     !<


#if defined( __parallel )
    
    dx = childgrid(m)%dx
    dy = childgrid(m)%dy
    IF ( PRESENT( dz ) )  THEN
       dz = childgrid(m)%dz
    ENDIF
    
#endif
    
 END SUBROUTINE get_child_gridspacing



 SUBROUTINE pmci_create_childs_parent_grid_arrays( name, is, ie, js, je, nzc, n  )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  ie      !<  RENAME ie, is, je, js?
    INTEGER(iwp), INTENT(IN) ::  is      !<
    INTEGER(iwp), INTENT(IN) ::  je      !<
    INTEGER(iwp), INTENT(IN) ::  js      !<
    INTEGER(iwp), INTENT(IN) ::  nzc     !<  nzc is pg%nz, but note that pg%nz is not the original nz of parent,
                                         !<  but the highest parent-grid level needed for nesting.
    INTEGER(iwp), INTENT(IN), OPTIONAL ::  n  !< number of chemical species / salsa variables
    
    CHARACTER(LEN=*), INTENT(IN) ::  name    !<

#if defined( __parallel )
!       
!-- Local variables:
    INTEGER(iwp) ::  ierr    !<
    
    INTEGER(idp), POINTER,DIMENSION(:,:)   ::  i_2d    !<
    
    REAL(wp), POINTER,DIMENSION(:,:)       ::  p_2d    !<
    REAL(wp), POINTER,DIMENSION(:,:,:)     ::  p_3d    !<
    
    NULLIFY( p_3d )
    NULLIFY( p_2d )
    NULLIFY( i_2d )
!
!-- List of array names, which can be coupled
    IF ( TRIM( name ) == "u" )  THEN
       IF ( .NOT. ALLOCATED( uc ) )  ALLOCATE( uc(0:nzc+1,js:je,is:ie) )
       p_3d => uc
    ELSEIF ( TRIM( name ) == "v" )  THEN
       IF ( .NOT. ALLOCATED( vc ) )  ALLOCATE( vc(0:nzc+1,js:je,is:ie) )
       p_3d => vc
    ELSEIF ( TRIM( name ) == "w" )  THEN
       IF ( .NOT. ALLOCATED( wc ) )  ALLOCATE( wc(0:nzc+1,js:je,is:ie) )
       p_3d => wc
    ELSEIF ( TRIM( name ) == "e" )  THEN
       IF ( .NOT. ALLOCATED( ec ) )  ALLOCATE( ec(0:nzc+1,js:je,is:ie) )
       p_3d => ec
    ELSEIF ( TRIM( name ) == "diss" )  THEN
       IF ( .NOT. ALLOCATED( dissc ) )  ALLOCATE( dissc(0:nzc+1,js:je,is:ie) )
       p_3d => dissc
    ELSEIF ( TRIM( name ) == "pt")  THEN
       IF ( .NOT. ALLOCATED( ptc ) )  ALLOCATE( ptc(0:nzc+1,js:je,is:ie) )
       p_3d => ptc
    ELSEIF ( TRIM( name ) == "q")  THEN
       IF ( .NOT. ALLOCATED( q_c ) )  ALLOCATE( q_c(0:nzc+1,js:je,is:ie) )
       p_3d => q_c
    ELSEIF ( TRIM( name ) == "qc")  THEN
       IF ( .NOT. ALLOCATED( qcc ) )  ALLOCATE( qcc(0:nzc+1,js:je,is:ie) )
       p_3d => qcc
    ELSEIF ( TRIM( name ) == "qr")  THEN
       IF ( .NOT. ALLOCATED( qrc ) )  ALLOCATE( qrc(0:nzc+1,js:je,is:ie) )
       p_3d => qrc
    ELSEIF ( TRIM( name ) == "nr")  THEN
       IF ( .NOT. ALLOCATED( nrc ) )  ALLOCATE( nrc(0:nzc+1,js:je,is:ie) )
       p_3d => nrc
    ELSEIF ( TRIM( name ) == "nc")  THEN
       IF ( .NOT. ALLOCATED( ncc ) )  ALLOCATE( ncc(0:nzc+1,js:je,is:ie) )
       p_3d => ncc
    ELSEIF ( TRIM( name ) == "s")  THEN
       IF ( .NOT. ALLOCATED( sc ) )  ALLOCATE( sc(0:nzc+1,js:je,is:ie) )
       p_3d => sc
    ELSEIF ( TRIM( name ) == "nr_part") THEN
       IF ( .NOT. ALLOCATED( nr_partc ) )  ALLOCATE( nr_partc(js:je,is:ie) )
       i_2d => nr_partc
    ELSEIF ( TRIM( name ) == "part_adr") THEN
       IF ( .NOT. ALLOCATED( part_adrc ) )  ALLOCATE( part_adrc(js:je,is:ie) )
       i_2d => part_adrc
    ELSEIF ( TRIM( name(1:5) ) == "chem_" )  THEN
       IF ( .NOT. ALLOCATED( chem_spec_c ) ) ALLOCATE( chem_spec_c(0:nzc+1,js:je,is:ie,1:nspec) )
       p_3d => chem_spec_c(:,:,:,n)
    ELSEIF ( TRIM( name(1:3) ) == "an_" )  THEN
       IF ( .NOT. ALLOCATED( aerosol_number_c ) )                              &
          ALLOCATE( aerosol_number_c(0:nzc+1,js:je,is:ie,1:nbins_aerosol) )
       p_3d => aerosol_number_c(:,:,:,n)
    ELSEIF ( TRIM( name(1:3) ) == "am_" )  THEN
       IF ( .NOT. ALLOCATED( aerosol_mass_c ) )                                &
          ALLOCATE( aerosol_mass_c(0:nzc+1,js:je,is:ie,1:(nbins_aerosol*ncomponents_mass) ) )
       p_3d => aerosol_mass_c(:,:,:,n)
    ELSEIF ( TRIM( name(1:3) ) == "sg_"  .AND.  .NOT. salsa_gases_from_chem )  &
    THEN
       IF ( .NOT. ALLOCATED( salsa_gas_c ) )                                   &
          ALLOCATE( salsa_gas_c(0:nzc+1,js:je,is:ie,1:ngases_salsa) )
       p_3d => salsa_gas_c(:,:,:,n)
    !ELSEIF (trim(name) == "z0") then
       !IF (.not.allocated(z0c))  allocate(z0c(js:je, is:ie))
       !p_2d => z0c
    ENDIF

    IF ( ASSOCIATED( p_3d ) )  THEN
       CALL pmc_c_set_dataarray( p_3d )
    ELSEIF ( ASSOCIATED( p_2d ) )  THEN
       CALL pmc_c_set_dataarray( p_2d )
    ELSEIF ( ASSOCIATED( i_2d ) )  THEN
       CALL pmc_c_set_dataarray( i_2d )
    ELSE
!
!--    Give only one message for the first child domain.
       IF ( cpl_id == 2  .AND.  myid == 0 )  THEN
          message_string = 'pointer for array "' // TRIM( name ) //            &
               '" can''t be associated'
          CALL message( 'pmci_create_childs_parent_grid_arrays', 'PA0170', 3, 2, 0, 6, 0 )
       ELSE
!
!--          Prevent others from continuing in case the abort is to come.
          CALL MPI_BARRIER( comm2d, ierr )
       ENDIF

    ENDIF

#endif
 END SUBROUTINE pmci_create_childs_parent_grid_arrays


 SUBROUTINE pmci_parent_initialize

!
!-- Send data for the children in order to let them create initial 
!-- conditions by interpolating the parent-domain fields.
#if defined( __parallel )
    IMPLICIT NONE

    INTEGER(iwp) ::  child_id    !<
    INTEGER(iwp) ::  m           !<
    REAL(wp) ::  waittime        !<


    DO  m = 1, SIZE( pmc_parent_for_child ) - 1
       child_id = pmc_parent_for_child(m)
       CALL pmc_s_fillbuffer( child_id, waittime=waittime )
    ENDDO

#endif
 END SUBROUTINE pmci_parent_initialize



 SUBROUTINE pmci_child_initialize

!
!-- Create initial conditions for the current child domain by interpolating 
!-- the parent-domain fields.
#if defined( __parallel )
    IMPLICIT NONE

    INTEGER(iwp) ::  ic         !< Child-grid index in x-direction
    INTEGER(iwp) ::  jc         !< Child-grid index in y-direction
    INTEGER(iwp) ::  kc         !< Child-grid index in z-direction
    INTEGER(iwp) ::  lb         !< Running index for aerosol size bins
    INTEGER(iwp) ::  lc         !< Running index for aerosol mass bins 
    INTEGER(iwp) ::  lg         !< Running index for salsa gases
    INTEGER(iwp) ::  n          !< Running index for chemical species
    REAL(wp) ::  waittime       !< Waiting time 

!
!-- Root model is never anyone's child
    IF ( .NOT.  pmc_is_rootmodel() )  THEN
!
!--    Get data from the parent
       CALL pmc_c_getbuffer( waittime = waittime )
!
!--    The interpolation.
       CALL pmci_interp_1sto_all ( u, uc, kcto, iflu, ifuu, jflo, jfuo, kflo, kfuo, 'u' )
       CALL pmci_interp_1sto_all ( v, vc, kcto, iflo, ifuo, jflv, jfuv, kflo, kfuo, 'v' )
       CALL pmci_interp_1sto_all ( w, wc, kctw, iflo, ifuo, jflo, jfuo, kflw, kfuw, 'w' )

       IF ( (        rans_mode_parent  .AND.         rans_mode )  .OR.                              &
            (  .NOT. rans_mode_parent  .AND.  .NOT.  rans_mode  .AND.                               &
               .NOT. constant_diffusion ) )  THEN
          CALL pmci_interp_1sto_all ( e, ec, kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, 'e' )
       ENDIF

       IF ( rans_mode_parent  .AND.  rans_mode  .AND.  rans_tke_e )  THEN
          CALL pmci_interp_1sto_all ( diss, dissc, kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, 's' )
       ENDIF

       IF ( .NOT. neutral )  THEN
          CALL pmci_interp_1sto_all ( pt, ptc, kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, 's' )
       ENDIF

       IF ( humidity )  THEN

          CALL pmci_interp_1sto_all ( q, q_c, kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, 's' )

          IF ( bulk_cloud_model  .AND.  microphysics_morrison )  THEN
             CALL pmci_interp_1sto_all ( qc, qcc, kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, 's' ) 
             CALL pmci_interp_1sto_all ( nc, ncc, kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, 's' )   
          ENDIF

          IF ( bulk_cloud_model  .AND.  microphysics_seifert )  THEN
             CALL pmci_interp_1sto_all ( qr, qrc, kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, 's' )
             CALL pmci_interp_1sto_all ( nr, nrc, kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, 's' )
          ENDIF

       ENDIF

       IF ( passive_scalar )  THEN
          CALL pmci_interp_1sto_all ( s, sc, kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, 's' )
       ENDIF

       IF ( air_chemistry  .AND.  nesting_chem )  THEN
          DO  n = 1, nspec
             CALL pmci_interp_1sto_all ( chem_species(n)%conc, chem_spec_c(:,:,:,n),                &
                                         kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, 's' )
          ENDDO
       ENDIF

       IF ( salsa  .AND.  nesting_salsa )  THEN
          DO  lb = 1, nbins_aerosol
             CALL pmci_interp_1sto_all ( aerosol_number(lb)%conc, aerosol_number_c(:,:,:,lb),       &
                                         kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, 's' )
          ENDDO
          DO  lc = 1, nbins_aerosol * ncomponents_mass
             CALL pmci_interp_1sto_all ( aerosol_mass(lc)%conc, aerosol_mass_c(:,:,:,lc),           &
                                         kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, 's' )
          ENDDO
          IF ( .NOT. salsa_gases_from_chem )  THEN
             DO  lg = 1, ngases_salsa
                CALL pmci_interp_1sto_all ( salsa_gas(lg)%conc, salsa_gas_c(:,:,:,lg),              &
                                            kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, 's' )
             ENDDO
          ENDIF
       ENDIF

       IF ( topography /= 'flat' )  THEN
!
!--       Inside buildings set velocities back to zero.
          DO  ic = nxlg, nxrg
             DO  jc = nysg, nyng
                DO  kc = nzb, nzt
                   u(kc,jc,ic)   = MERGE( u(kc,jc,ic), 0.0_wp, BTEST( wall_flags_total_0(kc,jc,ic), 1 ) )
                   v(kc,jc,ic)   = MERGE( v(kc,jc,ic), 0.0_wp, BTEST( wall_flags_total_0(kc,jc,ic), 2 ) )
                   w(kc,jc,ic)   = MERGE( w(kc,jc,ic), 0.0_wp, BTEST( wall_flags_total_0(kc,jc,ic), 3 ) )
                   u_p(kc,jc,ic) = MERGE( u_p(kc,jc,ic), 0.0_wp, BTEST( wall_flags_total_0(kc,jc,ic), 1 ) )
                   v_p(kc,jc,ic) = MERGE( v_p(kc,jc,ic), 0.0_wp, BTEST( wall_flags_total_0(kc,jc,ic), 2 ) )
                   w_p(kc,jc,ic) = MERGE( w_p(kc,jc,ic), 0.0_wp, BTEST( wall_flags_total_0(kc,jc,ic), 3 ) )
                ENDDO
             ENDDO
          ENDDO
       ENDIF
    ENDIF


 CONTAINS


    SUBROUTINE pmci_interp_1sto_all( child_array, parent_array, kct, ifl, ifu, jfl, jfu, kfl, kfu,  &
         var )
!
!--    Interpolation of the internal values for the child-domain initialization
       IMPLICIT NONE

       INTEGER(iwp), INTENT(IN) ::  kct  !< The parent-grid index in z-direction just below the boundary value node

       INTEGER(iwp), DIMENSION(ipla:ipra), INTENT(IN) ::  ifl  !<  Indicates start index of child cells belonging to certain
                                                               !<  parent cell - x direction
       INTEGER(iwp), DIMENSION(ipla:ipra), INTENT(IN) ::  ifu  !<  Indicates end index of child cells belonging to certain
                                                               !<  parent cell - x direction
       INTEGER(iwp), DIMENSION(jpsa:jpna), INTENT(IN) ::  jfl  !<  Indicates start index of child cells belonging to certain
                                                               !<  parent cell - y direction
       INTEGER(iwp), DIMENSION(jpsa:jpna), INTENT(IN) ::  jfu  !<  Indicates end index of child cells belonging to certain
                                                               !<  parent cell - y direction
       INTEGER(iwp), DIMENSION(0:pg%nz+1), INTENT(IN) ::  kfl  !<  Indicates start index of child cells belonging to certain
                                                               !<  parent cell - z direction
       INTEGER(iwp), DIMENSION(0:pg%nz+1), INTENT(IN) ::  kfu  !<  Indicates end index of child cells belonging to certain
                                                               !<  parent cell - z direction
       REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg), INTENT(INOUT) ::  child_array  !<  Child-grid array 
       REAL(wp), DIMENSION(0:pg%nz+1,jps:jpn,ipl:ipr), INTENT(IN) ::  parent_array        !<  Parent-grid array

       CHARACTER(LEN=1), INTENT(IN) ::  var  !<  Variable symbol: 'u', 'v', 'w' or 's'
!
!--    Local variables:
       INTEGER(iwp) ::  ic        !< Running child-grid index in the x-direction
       INTEGER(iwp) ::  icb       !< Index pointing to the first redundant ghost point layer behind the actual boundary
                                  !< ghost point layer in the x-direction
       INTEGER(iwp) ::  icbc      !< Index pointing to the boundary ghost point layer in the x-direction
       INTEGER(iwp) ::  icfirst   !< Leftmost child-grid index initialized by the main loops (usually icfirst == icl_init)
       INTEGER(iwp) ::  iclast    !< Rightmost child-grid index initialized by the main loops (usually iclast == icr_init)
       INTEGER(iwp) ::  icl_init  !< Left child-grid index bound for initialization in the x-direction
       INTEGER(iwp) ::  icr_init  !< Right child-grid index bound for initialization in the x-direction
       INTEGER(iwp) ::  jc        !< Running child-grid index in the y-direction
       INTEGER(iwp) ::  jcb       !< Index pointing to the first redundant ghost point layer behind the actual boundary
                                  !< ghost point layer in the y-direction
       INTEGER(iwp) ::  jcbc      !< Index pointing to the boundary ghost point layer in the y-direction
       INTEGER(iwp) ::  jcfirst   !< Southmost child-grid index initialized by the main loops (usually jcfirst == jcs_init)
       INTEGER(iwp) ::  jclast    !< Northmost child-grid index initialized by the main loops (usually jclast == jcn_init)
       INTEGER(iwp) ::  jcs_init  !< South child-grid index bound for initialization in the y-direction
       INTEGER(iwp) ::  jcn_init  !< North child-grid index bound for initialization in the y-direction
       INTEGER(iwp) ::  kc        !< Running child-grid index in the z-direction
       INTEGER(iwp) ::  ip        !< Running parent-grid index in the x-direction
       INTEGER(iwp) ::  ipl_init  !< Left parent-grid index bound for initialization in the x-direction
       INTEGER(iwp) ::  ipr_init  !< Right parent-grid index bound for initialization in the x-direction
       INTEGER(iwp) ::  jp        !< Running parent-grid index in the y-direction
       INTEGER(iwp) ::  jps_init  !< South parent-grid index bound for initialization in the y-direction
       INTEGER(iwp) ::  jpn_init  !< North parent-grid index bound for initialization in the y-direction
       INTEGER(iwp) ::  kp        !< Running parent-grid index in the z-direction


       ipl_init = ipl
       ipr_init = ipr
       jps_init = jps
       jpn_init = jpn
       icl_init = nxl
       icr_init = nxr
       jcs_init = nys
       jcn_init = nyn

       icbc = -1
       icb  = -2
       jcbc = -1
       jcb  = -2
       IF ( var == 'u' )  THEN
          icbc =  0
          icb  = -1
       ELSE IF ( var == 'v' )  THEN
          jcbc =  0
          jcb  = -1
       ENDIF
       
       IF ( nesting_mode /= 'vertical' )  THEN
          IF ( bc_dirichlet_l )  THEN
             ipl_init = ipl + 1
             icl_init = nxl - 1
!
!--          For u, nxl is a ghost node, but not for the other variables
             IF ( var == 'u' )  THEN
                ipl_init = ipl + 2
                icl_init = nxl
             ENDIF
          ENDIF
          IF ( bc_dirichlet_s )  THEN
             jps_init = jps + 1
             jcs_init = nys - 1 
!
!--          For v, nys is a ghost node, but not for the other variables
             IF ( var == 'v' )  THEN
                jps_init = jps + 2
                jcs_init = nys
             ENDIF
          ENDIF
          IF ( bc_dirichlet_r )  THEN
             ipr_init = ipr - 1
             icr_init = nxr + 1 
          ENDIF
          IF ( bc_dirichlet_n )  THEN
             jpn_init = jpn - 1
             jcn_init = nyn + 1
          ENDIF
       ENDIF      

       child_array(:,:,:) = 0.0_wp

       IF ( var == 'u' )  THEN

          icfirst = ifl(ipl_init) 
          iclast  = ifl(ipr_init+1) - 1
          jcfirst = jfl(jps_init)
          jclast  = jfu(jpn_init)
          DO  ip = ipl_init, ipr_init
             DO  jp = jps_init, jpn_init
                DO  kp = 0, kct + 1 

                   DO  ic = ifl(ip), ifl(ip+1)-1
                      DO  jc = jfl(jp), jfu(jp)
                         DO  kc = kfl(kp), MIN( kfu(kp), nzt+1 )
                            child_array(kc,jc,ic) = parent_array(kp,jp,ip)
                         ENDDO
                      ENDDO
                   ENDDO
                   
                ENDDO
             ENDDO
          ENDDO

       ELSE IF ( var == 'v' )  THEN

          icfirst = ifl(ipl_init) 
          iclast  = ifu(ipr_init)
          jcfirst = jfl(jps_init)
          jclast  = jfl(jpn_init+1) - 1
          DO  ip = ipl_init, ipr_init
             DO  jp = jps_init, jpn_init
                DO  kp = 0, kct + 1  

                   DO  ic = ifl(ip), ifu(ip)
                      DO  jc = jfl(jp), jfl(jp+1)-1
                         DO  kc = kfl(kp), MIN( kfu(kp), nzt+1 )
                            child_array(kc,jc,ic) = parent_array(kp,jp,ip) 
                         ENDDO
                      ENDDO
                   ENDDO

                ENDDO
             ENDDO
          ENDDO

       ELSE IF ( var == 'w' )  THEN

          icfirst = ifl(ipl_init) 
          iclast  = ifu(ipr_init)
          jcfirst = jfl(jps_init)
          jclast  = jfu(jpn_init)
          DO  ip = ipl_init, ipr_init
             DO  jp = jps_init, jpn_init
                DO  kp = 1, kct + 1  

                   DO  ic = ifl(ip), ifu(ip)
                      DO  jc = jfl(jp), jfu(jp)
!                         
!--                      Because the kp-loop for w starts from kp=1 instead of 0
                         child_array(nzb,jc,ic) = 0.0_wp
                         DO  kc = kfu(kp-1)+1, kfu(kp) 
                            child_array(kc,jc,ic) = parent_array(kp,jp,ip)
                         ENDDO
                      ENDDO
                   ENDDO
                   
                ENDDO
             ENDDO
          ENDDO

       ELSE   ! scalars

          icfirst = ifl(ipl_init) 
          iclast  = ifu(ipr_init)
          jcfirst = jfl(jps_init)
          jclast  = jfu(jpn_init)
          DO  ip = ipl_init, ipr_init
             DO  jp = jps_init, jpn_init
                DO  kp = 0, kct + 1
                                      
                   DO  ic = ifl(ip), ifu(ip)
                      DO  jc = jfl(jp), jfu(jp)                         
                         DO  kc = kfl(kp), MIN( kfu(kp), nzt+1 )
                            child_array(kc,jc,ic) = parent_array(kp,jp,ip)
                         ENDDO
                      ENDDO
                   ENDDO
                   
                ENDDO
             ENDDO
          ENDDO

       ENDIF  ! var  
!
!--    If the number of grid points in child subdomain in x- or y-direction
!--    (nxr - nxl + 1 and/or nyn - nys + 1) is not integer divisible by the grid spacing
!--    ratio in its direction (igsr and/or jgsr), the above loops will return with 
!--    unfilled gaps in the initial fields. These gaps, if present, are filled here.  
       IF ( icfirst > icl_init )  THEN
          DO  ic = icl_init, icfirst - 1
             child_array(:,:,ic) = child_array(:,:,icfirst)
          ENDDO
       ENDIF
       IF ( iclast < icr_init )  THEN
          DO  ic = iclast + 1, icr_init
             child_array(:,:,ic) = child_array(:,:,iclast)
          ENDDO
       ENDIF
       IF ( jcfirst > jcs_init )  THEN
          DO  jc = jcs_init, jcfirst - 1
             child_array(:,jc,:) = child_array(:,jcfirst,:)
          ENDDO
       ENDIF
       IF ( jclast < jcn_init )  THEN
          DO  jc = jclast + 1, jcn_init
             child_array(:,jc,:) = child_array(:,jclast,:)
          ENDDO
       ENDIF
!
!--    Finally, make sure that also the redundant 2nd and 3rd ghost-node layers
!--    including the corners are properly filled up.
       IF ( nys == 0 )  THEN
          DO  jc = -nbgp, jcb  ! jcb = -2 if var == v, else jcb = -1
             child_array(0:nzt+1,jc,nxlg:nxrg) = child_array(0:nzt+1,jcbc,nxlg:nxrg)
          ENDDO          
       ENDIF
       IF ( nyn == ny )  THEN
          DO  jc = ny+2, ny+nbgp
             child_array(0:nzt+1,jc,nxlg:nxrg) = child_array(0:nzt+1,ny+1,nxlg:nxrg)
          ENDDO
       ENDIF
       IF ( nxl == 0 )  THEN
          DO  ic = -nbgp, icb  ! icb = -2 if var == u, else icb = -1
             child_array(0:nzt+1,nysg:nyng,ic) = child_array(0:nzt+1,nysg:nyng,icbc)
          ENDDO          
       ENDIF
       IF ( nxr == nx )  THEN
          DO  ic = nx+2, nx+nbgp
             child_array(0:nzt+1,nysg:nyng,ic) = child_array(0:nzt+1,nysg:nyng,nx+1)
          ENDDO   
       ENDIF

    END SUBROUTINE pmci_interp_1sto_all

#endif
 END SUBROUTINE pmci_child_initialize



 SUBROUTINE pmci_check_setting_mismatches
!
!-- Check for mismatches between settings of root and child variables
!-- (e.g., all children have to follow the end_time settings of the root model).
!-- The root model overwrites variables in the other models, so these variables
!-- only need to be set once in file PARIN.

#if defined( __parallel )
    IMPLICIT NONE

    INTEGER ::  ierr                 !<  MPI error code

    REAL(wp) ::  dt_restart_root     !< 
    REAL(wp) ::  end_time_root       !<  
    REAL(wp) ::  restart_time_root   !< 
    REAL(wp) ::  time_restart_root   !<  

!
!-- Check the time to be simulated.
!-- Here, and in the following, the root process communicates the respective
!-- variable to all others, and its value will then be compared with the local
!-- values.
    IF ( pmc_is_rootmodel() )  end_time_root = end_time
    CALL MPI_BCAST( end_time_root, 1, MPI_REAL, 0, comm_world_nesting, ierr )

    IF ( .NOT. pmc_is_rootmodel() )  THEN
       IF ( end_time /= end_time_root )  THEN
          WRITE( message_string, * )  'mismatch between root model and ',                           &
               'child settings:& end_time(root) = ', end_time_root,                                 &
               '& end_time(child) = ', end_time, '& child value is set',                            &
               ' to root value'
          CALL message( 'pmci_check_setting_mismatches', 'PA0419', 0, 1, 0, 6,                      &
                        0 )
          end_time = end_time_root
       ENDIF
    ENDIF
!
!-- Same for restart time
    IF ( pmc_is_rootmodel() )  restart_time_root = restart_time
    CALL MPI_BCAST( restart_time_root, 1, MPI_REAL, 0, comm_world_nesting, ierr )

    IF ( .NOT. pmc_is_rootmodel() )  THEN
       IF ( restart_time /= restart_time_root )  THEN
          WRITE( message_string, * )  'mismatch between root model and ',      &
               'child settings: & restart_time(root) = ', restart_time_root,   &
               '& restart_time(child) = ', restart_time, '& child ',           &
               'value is set to root value'
          CALL message( 'pmci_check_setting_mismatches', 'PA0419', 0, 1, 0, 6, &
                        0 )
          restart_time = restart_time_root
       ENDIF
    ENDIF
!
!-- Same for dt_restart
    IF ( pmc_is_rootmodel() )  dt_restart_root = dt_restart
    CALL MPI_BCAST( dt_restart_root, 1, MPI_REAL, 0, comm_world_nesting, ierr )

    IF ( .NOT. pmc_is_rootmodel() )  THEN
       IF ( dt_restart /= dt_restart_root )  THEN
          WRITE( message_string, * )  'mismatch between root model and ',      &
               'child settings: & dt_restart(root) = ', dt_restart_root,       &
               '& dt_restart(child) = ', dt_restart, '& child ',               &
               'value is set to root value'
          CALL message( 'pmci_check_setting_mismatches', 'PA0419', 0, 1, 0, 6, &
                        0 )
          dt_restart = dt_restart_root
       ENDIF
    ENDIF
!
!-- Same for time_restart
    IF ( pmc_is_rootmodel() )  time_restart_root = time_restart
    CALL MPI_BCAST( time_restart_root, 1, MPI_REAL, 0, comm_world_nesting, ierr )

    IF ( .NOT. pmc_is_rootmodel() )  THEN
       IF ( time_restart /= time_restart_root )  THEN
          WRITE( message_string, * )  'mismatch between root model and ',      &
               'child settings: & time_restart(root) = ', time_restart_root,   &
               '& time_restart(child) = ', time_restart, '& child ',           &
               'value is set to root value'
          CALL message( 'pmci_check_setting_mismatches', 'PA0419', 0, 1, 0, 6, &
                        0 )
          time_restart = time_restart_root
       ENDIF
    ENDIF

#endif

 END SUBROUTINE pmci_check_setting_mismatches



 SUBROUTINE pmci_synchronize

#if defined( __parallel )
!
!-- Unify the time steps for each model and synchronize using 
!-- MPI_ALLREDUCE with the MPI_MIN operator over all processes using 
!-- the global communicator MPI_COMM_WORLD.
   
   IMPLICIT NONE

   INTEGER(iwp) ::  ierr  !<  MPI error code
   REAL(wp) ::  dtl       !<  Local time step of the current process
   REAL(wp) ::  dtg       !<  Global time step defined as the global minimum of dtl of all processes


   IF ( debug_output_timestep )  CALL debug_message( 'pmci_synchronize', 'start' )
   
   dtl = dt_3d
   CALL MPI_ALLREDUCE( dtl, dtg, 1, MPI_REAL, MPI_MIN, MPI_COMM_WORLD, ierr )
   dt_3d  = dtg

   IF ( debug_output_timestep )  CALL debug_message( 'pmci_synchronize', 'end' )

#endif
 END SUBROUTINE pmci_synchronize
                


 SUBROUTINE pmci_set_swaplevel( swaplevel )

!
!-- After each Runge-Kutta sub-timestep, alternately set buffer one or buffer
!-- two active

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  swaplevel  !< swaplevel (1 or 2) of PALM's timestep

    INTEGER(iwp) ::  child_id    !<  Child id of the child number m
    INTEGER(iwp) ::  m           !<  Loop index over all children of the current parent 

#if defined( __parallel )
    DO  m = 1, SIZE( pmc_parent_for_child ) - 1
       child_id = pmc_parent_for_child(m)
       CALL pmc_s_set_active_data_array( child_id, swaplevel )
    ENDDO
#endif
 END SUBROUTINE pmci_set_swaplevel



 SUBROUTINE pmci_datatrans( local_nesting_mode )   
!
!-- This subroutine controls the nesting according to the nestpar 
!-- parameter nesting_mode (two-way (default) or one-way) and the 
!-- order of anterpolations according to the nestpar parameter 
!-- nesting_datatransfer_mode (cascade, overlap or mixed (default)).
!-- Although nesting_mode is a variable of this model, pass it as 
!-- an argument to allow for example to force one-way initialization 
!-- phase.
!-- Note that interpolation ( parent_to_child ) must always be carried 
!-- out before anterpolation ( child_to_parent ).

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) ::  local_nesting_mode  !<  Nesting mode: 'one-way', 'two-way' or 'vertical'

#if defined( __parallel )    

    IF ( debug_output_timestep )  CALL debug_message( 'pmci_datatrans', 'start' )

    IF ( TRIM( local_nesting_mode ) == 'one-way' )  THEN

       CALL pmci_child_datatrans( parent_to_child )
       CALL pmci_parent_datatrans( parent_to_child )

    ELSE

       IF ( nesting_datatransfer_mode == 'cascade' )  THEN

          CALL pmci_child_datatrans( parent_to_child )
          CALL pmci_parent_datatrans( parent_to_child )

          CALL pmci_parent_datatrans( child_to_parent )
          CALL pmci_child_datatrans( child_to_parent )

       ELSEIF ( nesting_datatransfer_mode == 'overlap')  THEN

          CALL pmci_parent_datatrans( parent_to_child )
          CALL pmci_child_datatrans( parent_to_child )

          CALL pmci_child_datatrans( child_to_parent )
          CALL pmci_parent_datatrans( child_to_parent )

       ELSEIF ( TRIM( nesting_datatransfer_mode ) == 'mixed' )  THEN

          CALL pmci_parent_datatrans( parent_to_child )
          CALL pmci_child_datatrans( parent_to_child )

          CALL pmci_parent_datatrans( child_to_parent )
          CALL pmci_child_datatrans( child_to_parent )

       ENDIF

    ENDIF

    IF ( debug_output_timestep )  CALL debug_message( 'pmci_datatrans', 'end' )

#endif
 END SUBROUTINE pmci_datatrans



 SUBROUTINE pmci_parent_datatrans( direction )
   
    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  direction   !<  Direction of the data transfer: 'parent_to_child' or 'child_to_parent'

#if defined( __parallel )
    INTEGER(iwp) ::  child_id    !<  Child id of the child number m
    INTEGER(iwp) ::  i           !<  Parent-grid index in x-direction
    INTEGER(iwp) ::  j           !<  Parent-grid index in y-direction
    INTEGER(iwp) ::  k           !<  Parent-grid index in z-direction
    INTEGER(iwp) ::  m           !<  Loop index over all children of the current parent 


    DO  m = 1, SIZE( pmc_parent_for_child ) - 1
       child_id = pmc_parent_for_child(m)
       IF ( direction == parent_to_child )  THEN
          CALL cpu_log( log_point_s(71), 'pmc parent send', 'start' )
          CALL pmc_s_fillbuffer( child_id )
          CALL cpu_log( log_point_s(71), 'pmc parent send', 'stop' )
       ELSE
!
!--       Communication from child to parent
          CALL cpu_log( log_point_s(72), 'pmc parent recv', 'start' )
          child_id = pmc_parent_for_child(m)
          CALL pmc_s_getdata_from_buffer( child_id )
          CALL cpu_log( log_point_s(72), 'pmc parent recv', 'stop' )
!
!--       The anterpolated data is now available in u etc
          IF ( topography /= 'flat' )  THEN
!
!--          Inside buildings/topography reset velocities back to zero.
!--          Scalars (pt, q, s, km, kh, p, sa, ...) are ignored at
!--          present, maybe revise later.
!--          Resetting of e is removed as unnecessary since e is not 
!--          anterpolated, and as incorrect since it overran the default 
!--          Neumann condition (bc_e_b). 
             DO   i = nxlg, nxrg
                DO   j = nysg, nyng
                   DO  k = nzb, nzt+1
                      u(k,j,i) = MERGE( u(k,j,i), 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 1 ) )
                      v(k,j,i) = MERGE( v(k,j,i), 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 2 ) )
                      w(k,j,i) = MERGE( w(k,j,i), 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 3 ) )
!
!--                 TO_DO: zero setting of temperature within topography creates
!--                       wrong results
!                   pt(nzb:nzb_s_inner(j,i),j,i) = 0.0_wp
!                   IF ( humidity  .OR.  passive_scalar )  THEN
!                      q(nzb:nzb_s_inner(j,i),j,i) = 0.0_wp
!                   ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDIF
    ENDDO  ! m

#endif
 END SUBROUTINE pmci_parent_datatrans



 SUBROUTINE pmci_child_datatrans( direction )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  direction  !< Transfer direction: parent_to_child or child_to_parent

#if defined( __parallel )

    REAL(wp), DIMENSION(1) ::  dtl          !< Time step size


    dtl = dt_3d
    IF ( .NOT.  pmc_is_rootmodel() )  THEN

       IF ( direction == parent_to_child )  THEN
    
          CALL cpu_log( log_point_s(73), 'pmc child recv', 'start' )
          CALL pmc_c_getbuffer( )
          CALL cpu_log( log_point_s(73), 'pmc child recv', 'stop' )

          CALL cpu_log( log_point_s(75), 'pmc interpolation', 'start' )
          CALL pmci_interpolation
          CALL cpu_log( log_point_s(75), 'pmc interpolation', 'stop' )
     
       ELSE
!
!--       direction == child_to_parent
          CALL cpu_log( log_point_s(76), 'pmc anterpolation', 'start' )
          CALL pmci_anterpolation
          CALL cpu_log( log_point_s(76), 'pmc anterpolation', 'stop' )

          CALL cpu_log( log_point_s(74), 'pmc child send', 'start' )
          CALL pmc_c_putbuffer( )
          CALL cpu_log( log_point_s(74), 'pmc child send', 'stop' )

       ENDIF
    ENDIF

 CONTAINS

   
    SUBROUTINE pmci_interpolation

!
!--    A wrapper routine for all interpolation actions
      
       IMPLICIT NONE

       INTEGER(iwp) ::  ibgp       !< Index running over the nbgp boundary ghost points in i-direction
       INTEGER(iwp) ::  jbgp       !< Index running over the nbgp boundary ghost points in j-direction
       INTEGER(iwp) ::  lb         !< Running index for aerosol size bins
       INTEGER(iwp) ::  lc         !< Running index for aerosol mass bins
       INTEGER(iwp) ::  lg         !< Running index for salsa gases
       INTEGER(iwp) ::  n          !< Running index for number of chemical species
      
!
!--    In case of vertical nesting no interpolation is needed for the
!--    horizontal boundaries
       IF ( nesting_mode /= 'vertical' )  THEN
!
!--       Left border pe:
          IF ( bc_dirichlet_l )  THEN

             CALL pmci_interp_1sto_lr( u, uc, kcto, jflo, jfuo, kflo, kfuo, 'l', 'u' )
             CALL pmci_interp_1sto_lr( v, vc, kcto, jflv, jfuv, kflo, kfuo, 'l', 'v' )
             CALL pmci_interp_1sto_lr( w, wc, kctw, jflo, jfuo, kflw, kfuw, 'l', 'w' )

             IF ( (         rans_mode_parent  .AND.         rans_mode )  .OR.                       &
                  (  .NOT.  rans_mode_parent  .AND.  .NOT.  rans_mode  .AND.                        &
                     .NOT.  constant_diffusion ) )  THEN
!                CALL pmci_interp_1sto_lr( e, ec, kcto, jflo, jfuo, kflo, kfuo, 'l', 'e' )
!
!--             Interpolation of e is replaced by the Neumann condition. 
                DO  ibgp = -nbgp, -1
                   e(nzb:nzt,nys:nyn,ibgp) = e(nzb:nzt,nys:nyn,0)
                ENDDO

             ENDIF

             IF ( rans_mode_parent  .AND.  rans_mode  .AND.  rans_tke_e )  THEN
                CALL pmci_interp_1sto_lr( diss, dissc, kcto, jflo, jfuo, kflo, kfuo, 'l', 's' )
             ENDIF

             IF ( .NOT. neutral )  THEN
                CALL pmci_interp_1sto_lr( pt, ptc, kcto, jflo, jfuo, kflo, kfuo, 'l', 's' )
             ENDIF

             IF ( humidity )  THEN

                CALL pmci_interp_1sto_lr( q, q_c, kcto, jflo, jfuo, kflo, kfuo, 'l', 's' )

                IF ( bulk_cloud_model  .AND.  microphysics_morrison )  THEN
                   CALL pmci_interp_1sto_lr( qc, qcc, kcto, jflo, jfuo, kflo, kfuo, 'l', 's' )  
                   CALL pmci_interp_1sto_lr( nc, ncc, kcto, jflo, jfuo, kflo, kfuo, 'l', 's' )          
                ENDIF

                IF ( bulk_cloud_model  .AND.  microphysics_seifert )  THEN
                   CALL pmci_interp_1sto_lr( qr, qrc, kcto, jflo, jfuo, kflo, kfuo, 'l', 's' ) 
                   CALL pmci_interp_1sto_lr( nr, nrc, kcto, jflo, jfuo, kflo, kfuo, 'l', 's' )             
                ENDIF

             ENDIF

             IF ( passive_scalar )  THEN
                CALL pmci_interp_1sto_lr( s, sc, kcto, jflo, jfuo, kflo, kfuo, 'l', 's' )
             ENDIF

             IF ( air_chemistry  .AND.  nesting_chem )  THEN
                DO  n = 1, nspec
                   CALL pmci_interp_1sto_lr( chem_species(n)%conc, chem_spec_c(:,:,:,n),            &
                        kcto, jflo, jfuo, kflo, kfuo, 'l', 's' )
                ENDDO
             ENDIF

             IF ( salsa  .AND.  nesting_salsa )  THEN
                DO  lb = 1, nbins_aerosol
                   CALL pmci_interp_1sto_lr( aerosol_number(lb)%conc, aerosol_number_c(:,:,:,lb),   &
                                             kcto, jflo, jfuo, kflo, kfuo, 'l', 's')
                ENDDO
                DO  lc = 1, nbins_aerosol * ncomponents_mass
                   CALL pmci_interp_1sto_lr( aerosol_mass(lc)%conc, aerosol_mass_c(:,:,:,lc),       &
                                             kcto, jflo, jfuo, kflo, kfuo, 'l', 's')
                ENDDO
                IF ( .NOT. salsa_gases_from_chem )  THEN
                   DO  lg = 1, ngases_salsa
                      CALL pmci_interp_1sto_lr( salsa_gas(lg)%conc, salsa_gas_c(:,:,:,lg),          &
                                                kcto, jflo, jfuo, kflo, kfuo, 'l', 's')
                   ENDDO
                ENDIF
             ENDIF

          ENDIF
!
!--       Right border pe
          IF ( bc_dirichlet_r )  THEN
             
             CALL pmci_interp_1sto_lr( u, uc, kcto, jflo, jfuo, kflo, kfuo, 'r', 'u' )
             CALL pmci_interp_1sto_lr( v, vc, kcto, jflv, jfuv, kflo, kfuo, 'r', 'v' )
             CALL pmci_interp_1sto_lr( w, wc, kctw, jflo, jfuo, kflw, kfuw, 'r', 'w' )

             IF ( (         rans_mode_parent  .AND.         rans_mode )  .OR.                       &
                  (  .NOT.  rans_mode_parent  .AND.  .NOT.  rans_mode  .AND.                        &
                     .NOT.  constant_diffusion ) )  THEN
!                CALL pmci_interp_1sto_lr( e, ec, kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, 'r', 'e' )
!
!--             Interpolation of e is replaced by the Neumann condition. 
                DO  ibgp = nx+1, nx+nbgp
                   e(nzb:nzt,nys:nyn,ibgp) = e(nzb:nzt,nys:nyn,nx)
                ENDDO
             ENDIF

             IF ( rans_mode_parent  .AND.  rans_mode  .AND.  rans_tke_e )  THEN
                CALL pmci_interp_1sto_lr( diss, dissc, kcto, jflo, jfuo, kflo, kfuo, 'r', 's' )
             ENDIF

             IF (  .NOT.  neutral )  THEN
                CALL pmci_interp_1sto_lr( pt, ptc, kcto, jflo, jfuo, kflo, kfuo, 'r', 's' )
             ENDIF

             IF ( humidity )  THEN
                CALL pmci_interp_1sto_lr( q, q_c, kcto, jflo, jfuo, kflo, kfuo, 'r', 's' )

                IF ( bulk_cloud_model  .AND.  microphysics_morrison )  THEN
                   CALL pmci_interp_1sto_lr( qc, qcc, kcto, jflo, jfuo, kflo, kfuo, 'r', 's' ) 
                   CALL pmci_interp_1sto_lr( nc, ncc, kcto, jflo, jfuo, kflo, kfuo, 'r', 's' )
                ENDIF

                IF ( bulk_cloud_model  .AND.  microphysics_seifert )  THEN
                   CALL pmci_interp_1sto_lr( qr, qrc, kcto, jflo, jfuo, kflo, kfuo, 'r', 's' ) 
                   CALL pmci_interp_1sto_lr( nr, nrc, kcto, jflo, jfuo, kflo, kfuo, 'r', 's' ) 
                ENDIF

             ENDIF

             IF ( passive_scalar )  THEN
                CALL pmci_interp_1sto_lr( s, sc, kcto, jflo, jfuo, kflo, kfuo, 'r', 's' )
             ENDIF

             IF ( air_chemistry  .AND.  nesting_chem )  THEN
                DO  n = 1, nspec
                   CALL pmci_interp_1sto_lr( chem_species(n)%conc, chem_spec_c(:,:,:,n),            &
                        kcto, jflo, jfuo, kflo, kfuo, 'r', 's' )
                ENDDO
             ENDIF

             IF ( salsa  .AND.  nesting_salsa )  THEN
                DO  lb = 1, nbins_aerosol
                   CALL pmci_interp_1sto_lr( aerosol_number(lb)%conc, aerosol_number_c(:,:,:,lb),   &
                                             kcto, jflo, jfuo, kflo, kfuo, 'r', 's' )
                ENDDO
                DO  lc = 1, nbins_aerosol * ncomponents_mass
                   CALL pmci_interp_1sto_lr( aerosol_mass(lc)%conc, aerosol_mass_c(:,:,:,lc),       &
                                             kcto, jflo, jfuo, kflo, kfuo, 'r', 's' )
                ENDDO
                IF ( .NOT. salsa_gases_from_chem )  THEN
                   DO  lg = 1, ngases_salsa
                      CALL pmci_interp_1sto_lr( salsa_gas(lg)%conc, salsa_gas_c(:,:,:,lg),          &
                                                kcto, jflo, jfuo, kflo, kfuo, 'r', 's' )
                   ENDDO
                ENDIF
             ENDIF

          ENDIF
!
!--       South border pe
          IF ( bc_dirichlet_s )  THEN

             CALL pmci_interp_1sto_sn( v, vc, kcto, iflo, ifuo, kflo, kfuo, 's', 'v' )
             CALL pmci_interp_1sto_sn( w, wc, kctw, iflo, ifuo, kflw, kfuw, 's', 'w' )
             CALL pmci_interp_1sto_sn( u, uc, kcto, iflu, ifuu, kflo, kfuo, 's', 'u' )

             IF ( (         rans_mode_parent  .AND.         rans_mode )  .OR.                       &
                  (  .NOT.  rans_mode_parent  .AND.  .NOT.  rans_mode  .AND.                        &
                     .NOT.  constant_diffusion ) )  THEN
!                CALL pmci_interp_1sto_sn( e, ec, kcto, iflo, ifuo, kflo, kfuo, 's', 'e' )
!
!--             Interpolation of e is replaced by the Neumann condition. 
                DO  jbgp = -nbgp, -1
                   e(nzb:nzt,jbgp,nxl:nxr) = e(nzb:nzt,0,nxl:nxr)
                ENDDO
             ENDIF

             IF ( rans_mode_parent  .AND.  rans_mode  .AND.  rans_tke_e )  THEN
                CALL pmci_interp_1sto_sn( diss, dissc, kcto, iflo, ifuo, kflo, kfuo, 's', 's' )
             ENDIF

             IF (  .NOT.  neutral )  THEN
                CALL pmci_interp_1sto_sn( pt, ptc, kcto, iflo, ifuo, kflo, kfuo, 's', 's' )
             ENDIF

             IF ( humidity )  THEN
                CALL pmci_interp_1sto_sn( q, q_c, kcto, iflo, ifuo, kflo, kfuo, 's', 's' )

                IF ( bulk_cloud_model  .AND.  microphysics_morrison )  THEN
                   CALL pmci_interp_1sto_sn( qc, qcc, kcto, iflo, ifuo, kflo, kfuo, 's', 's' )
                   CALL pmci_interp_1sto_sn( nc, ncc, kcto, iflo, ifuo, kflo, kfuo, 's', 's' )
                ENDIF

                IF ( bulk_cloud_model  .AND.  microphysics_seifert )  THEN
                   CALL pmci_interp_1sto_sn( qr, qrc, kcto, iflo, ifuo, kflo, kfuo, 's', 's' )
                   CALL pmci_interp_1sto_sn( nr, nrc, kcto, iflo, ifuo, kflo, kfuo, 's', 's' )
                ENDIF

             ENDIF

             IF ( passive_scalar )  THEN
                CALL pmci_interp_1sto_sn( s, sc, kcto, iflo, ifuo, kflo, kfuo, 's', 's' )
             ENDIF

             IF ( air_chemistry  .AND.  nesting_chem )  THEN
                DO  n = 1, nspec
                   CALL pmci_interp_1sto_sn( chem_species(n)%conc, chem_spec_c(:,:,:,n),            &
                        kcto, iflo, ifuo, kflo, kfuo, 's', 's' )
                ENDDO
             ENDIF
             
             IF ( salsa  .AND.  nesting_salsa )  THEN
                DO  lb = 1, nbins_aerosol
                   CALL pmci_interp_1sto_sn( aerosol_number(lb)%conc, aerosol_number_c(:,:,:,lb),   &
                                             kcto, iflo, ifuo, kflo, kfuo, 's', 's' )
                ENDDO
                DO  lc = 1, nbins_aerosol * ncomponents_mass
                   CALL pmci_interp_1sto_sn( aerosol_mass(lc)%conc, aerosol_mass_c(:,:,:,lc),       &
                                             kcto, iflo, ifuo, kflo, kfuo, 's', 's' )
                ENDDO
                IF ( .NOT. salsa_gases_from_chem )  THEN
                   DO  lg = 1, ngases_salsa
                      CALL pmci_interp_1sto_sn( salsa_gas(lg)%conc, salsa_gas_c(:,:,:,lg),          &
                                                kcto, iflo, ifuo, kflo, kfuo, 's', 's' )
                   ENDDO
                ENDIF
             ENDIF
                       
          ENDIF
!
!--       North border pe
          IF ( bc_dirichlet_n )  THEN
             
             CALL pmci_interp_1sto_sn( v, vc, kcto, iflo, ifuo, kflo, kfuo, 'n', 'v' )
             CALL pmci_interp_1sto_sn( w, wc, kctw, iflo, ifuo, kflw, kfuw, 'n', 'w' )
             CALL pmci_interp_1sto_sn( u, uc, kcto, iflu, ifuu, kflo, kfuo, 'n', 'u' )

             IF ( (         rans_mode_parent  .AND.         rans_mode )  .OR.                       & 
                  (  .NOT.  rans_mode_parent  .AND.  .NOT.  rans_mode  .AND.                        &
                     .NOT.  constant_diffusion ) )  THEN
!                CALL pmci_interp_1sto_sn( e, ec, kcto, iflo, ifuo, kflo, kfuo, 'n', 'e' )
!
!--             Interpolation of e is replaced by the Neumann condition. 
                DO  jbgp = ny+1, ny+nbgp
                   e(nzb:nzt,jbgp,nxl:nxr) = e(nzb:nzt,ny,nxl:nxr)
                ENDDO
             ENDIF

             IF ( rans_mode_parent  .AND.  rans_mode  .AND.  rans_tke_e )  THEN
                CALL pmci_interp_1sto_sn( diss, dissc, kcto, iflo, ifuo, kflo, kfuo, 'n', 's' )
             ENDIF

             IF (  .NOT.  neutral )  THEN
                CALL pmci_interp_1sto_sn( pt, ptc, kcto, iflo, ifuo, kflo, kfuo, 'n', 's' )
             ENDIF

             IF ( humidity )  THEN
                CALL pmci_interp_1sto_sn( q, q_c, kcto, iflo, ifuo, kflo, kfuo, 'n', 's' )

                IF ( bulk_cloud_model  .AND.  microphysics_morrison )  THEN
                   CALL pmci_interp_1sto_sn( qc, qcc, kcto, iflo, ifuo, kflo, kfuo, 'n', 's' )
                   CALL pmci_interp_1sto_sn( nc, ncc, kcto, iflo, ifuo, kflo, kfuo, 'n', 's' )
                ENDIF

                IF ( bulk_cloud_model  .AND.  microphysics_seifert )  THEN
                   CALL pmci_interp_1sto_sn( qr, qrc, kcto, iflo, ifuo, kflo, kfuo, 'n', 's' )
                   CALL pmci_interp_1sto_sn( nr, nrc, kcto, iflo, ifuo, kflo, kfuo, 'n', 's' )
                ENDIF

             ENDIF

             IF ( passive_scalar )  THEN
                CALL pmci_interp_1sto_sn( s, sc, kcto, iflo, ifuo, kflo, kfuo, 'n', 's' )
             ENDIF

             IF ( air_chemistry  .AND.  nesting_chem )  THEN
                DO  n = 1, nspec
                   CALL pmci_interp_1sto_sn( chem_species(n)%conc, chem_spec_c(:,:,:,n),            &
                        kcto, iflo, ifuo, kflo, kfuo, 'n', 's' )
                ENDDO
             ENDIF
             
             IF ( salsa  .AND.  nesting_salsa )  THEN
                DO  lb = 1, nbins_aerosol
                   CALL pmci_interp_1sto_sn( aerosol_number(lb)%conc, aerosol_number_c(:,:,:,lb),   &
                                             kcto, iflo, ifuo, kflo, kfuo, 'n', 's' )
                ENDDO
                DO  lc = 1, nbins_aerosol * ncomponents_mass
                   CALL pmci_interp_1sto_sn( aerosol_mass(lc)%conc, aerosol_mass_c(:,:,:,lc),       &
                                             kcto, iflo, ifuo, kflo, kfuo, 'n', 's' )
                ENDDO
                IF ( .NOT. salsa_gases_from_chem )  THEN
                   DO  lg = 1, ngases_salsa
                      CALL pmci_interp_1sto_sn( salsa_gas(lg)%conc, salsa_gas_c(:,:,:,lg),          &
                                                kcto, iflo, ifuo, kflo, kfuo, 'n', 's' )
                   ENDDO
                ENDIF
             ENDIF
                          
          ENDIF

       ENDIF       ! IF ( nesting_mode /= 'vertical' )
!
!--    All PEs are top-border PEs
       CALL pmci_interp_1sto_t( w, wc, kctw, iflo, ifuo, jflo, jfuo, 'w' )
       CALL pmci_interp_1sto_t( u, uc, kcto, iflu, ifuu, jflo, jfuo, 'u' )
       CALL pmci_interp_1sto_t( v, vc, kcto, iflo, ifuo, jflv, jfuv, 'v' )

       IF ( (         rans_mode_parent  .AND.         rans_mode )  .OR.                             &
            (  .NOT.  rans_mode_parent  .AND.  .NOT.  rans_mode  .AND.                              &
               .NOT.  constant_diffusion ) )  THEN
!          CALL pmci_interp_1sto_t( e, ec, kcto, iflo, ifuo, jflo, jfuo, 'e' )
!
!--       Interpolation of e is replaced by the Neumann condition. 
          e(nzt+1,nys:nyn,nxl:nxr) = e(nzt,nys:nyn,nxl:nxr)
       ENDIF

       IF ( rans_mode_parent  .AND.  rans_mode  .AND.  rans_tke_e )  THEN
          CALL pmci_interp_1sto_t( diss, dissc, kcto, iflo, ifuo, jflo, jfuo, 's' )
       ENDIF

       IF ( .NOT. neutral )  THEN
          CALL pmci_interp_1sto_t( pt, ptc, kcto, iflo, ifuo, jflo, jfuo, 's' )
       ENDIF

       IF ( humidity )  THEN
          CALL pmci_interp_1sto_t( q, q_c, kcto, iflo, ifuo, jflo, jfuo, 's' )
          IF ( bulk_cloud_model  .AND.  microphysics_morrison )  THEN
             CALL pmci_interp_1sto_t( qc, qcc, kcto, iflo, ifuo, jflo, jfuo, 's' )
             CALL pmci_interp_1sto_t( nc, ncc, kcto, iflo, ifuo, jflo, jfuo, 's' )
          ENDIF
          IF ( bulk_cloud_model  .AND.  microphysics_seifert )  THEN
             CALL pmci_interp_1sto_t( qr, qrc, kcto, iflo, ifuo, jflo, jfuo, 's' )
             CALL pmci_interp_1sto_t( nr, nrc, kcto, iflo, ifuo, jflo, jfuo, 's' )
          ENDIF
       ENDIF

       IF ( passive_scalar )  THEN
          CALL pmci_interp_1sto_t( s, sc, kcto, iflo, ifuo, jflo, jfuo, 's' )
       ENDIF

       IF ( air_chemistry  .AND.  nesting_chem )  THEN
          DO  n = 1, nspec
             CALL pmci_interp_1sto_t( chem_species(n)%conc, chem_spec_c(:,:,:,n),                   &
                                      kcto, iflo, ifuo, jflo, jfuo, 's' )
          ENDDO
       ENDIF 
       
       IF ( salsa  .AND.  nesting_salsa )  THEN
          DO  lb = 1, nbins_aerosol
             CALL pmci_interp_1sto_t( aerosol_number(lb)%conc, aerosol_number_c(:,:,:,lb),          &
                                      kcto, iflo, ifuo, jflo, jfuo, 's' )
          ENDDO
          DO  lc = 1, nbins_aerosol * ncomponents_mass
             CALL pmci_interp_1sto_t( aerosol_mass(lc)%conc, aerosol_mass_c(:,:,:,lc),              &
                                      kcto, iflo, ifuo, jflo, jfuo, 's' )
          ENDDO
          IF ( .NOT. salsa_gases_from_chem )  THEN
             DO  lg = 1, ngases_salsa
                CALL pmci_interp_1sto_t( salsa_gas(lg)%conc, salsa_gas_c(:,:,:,lg),                 &
                                         kcto, iflo, ifuo, jflo, jfuo, 's' )
             ENDDO
          ENDIF
       ENDIF

   END SUBROUTINE pmci_interpolation



   SUBROUTINE pmci_anterpolation

!
!--   A wrapper routine for all anterpolation actions.
!--   Note that TKE is not anterpolated.
      IMPLICIT NONE
      INTEGER(iwp) ::  lb         !< Running index for aerosol size bins
      INTEGER(iwp) ::  lc         !< Running index for aerosol mass bins
      INTEGER(iwp) ::  lg         !< Running index for salsa gases
      INTEGER(iwp) ::  n          !< Running index for number of chemical species

      
      CALL pmci_anterp_tophat( u,  uc,  kcto, iflu, ifuu, jflo, jfuo, kflo, kfuo, ijkfc_u, 'u' )
      CALL pmci_anterp_tophat( v,  vc,  kcto, iflo, ifuo, jflv, jfuv, kflo, kfuo, ijkfc_v, 'v' )
      CALL pmci_anterp_tophat( w,  wc,  kctw, iflo, ifuo, jflo, jfuo, kflw, kfuw, ijkfc_w, 'w' )
!
!--   Anterpolation of TKE and dissipation rate if parent and child are in 
!--   RANS mode.
      IF ( rans_mode_parent  .AND.  rans_mode )  THEN
         CALL pmci_anterp_tophat( e, ec, kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, ijkfc_s, 'e' )
!
!--      Anterpolation of dissipation rate only if TKE-e closure is applied.
         IF ( rans_tke_e )  THEN
            CALL pmci_anterp_tophat( diss, dissc, kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo,         &
                 ijkfc_s, 'diss' )
         ENDIF

      ENDIF

      IF ( .NOT. neutral )  THEN
         CALL pmci_anterp_tophat( pt, ptc, kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, ijkfc_s, 'pt' )
      ENDIF

      IF ( humidity )  THEN

         CALL pmci_anterp_tophat( q, q_c, kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, ijkfc_s, 'q' )

         IF ( bulk_cloud_model  .AND.  microphysics_morrison )  THEN

            CALL pmci_anterp_tophat( qc, qcc, kcto, iflo, ifuo, jflo, jfuo,                         &
                                     kflo, kfuo, ijkfc_s, 'qc' )
            
            CALL pmci_anterp_tophat( nc, ncc, kcto, iflo, ifuo, jflo, jfuo,                         &
                                     kflo, kfuo, ijkfc_s, 'nc' )

         ENDIF

         IF ( bulk_cloud_model  .AND.  microphysics_seifert )  THEN

            CALL pmci_anterp_tophat( qr, qrc, kcto, iflo, ifuo, jflo, jfuo,                         &
                                     kflo, kfuo, ijkfc_s, 'qr' )

            CALL pmci_anterp_tophat( nr, nrc, kcto, iflo, ifuo, jflo, jfuo,                         &
                                     kflo, kfuo, ijkfc_s, 'nr' )

         ENDIF

      ENDIF

      IF ( passive_scalar )  THEN
         CALL pmci_anterp_tophat( s, sc, kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, ijkfc_s, 's' )
      ENDIF

      IF ( air_chemistry  .AND.  nesting_chem )  THEN
         DO  n = 1, nspec
            CALL pmci_anterp_tophat( chem_species(n)%conc, chem_spec_c(:,:,:,n),                    &
                                     kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, ijkfc_s, 's' )
         ENDDO
      ENDIF

      IF ( salsa  .AND.  nesting_salsa )  THEN
         DO  lb = 1, nbins_aerosol
            CALL pmci_anterp_tophat( aerosol_number(lb)%conc, aerosol_number_c(:,:,:,lb),           &
                                     kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, ijkfc_s, 's' )
         ENDDO
         DO  lc = 1, nbins_aerosol * ncomponents_mass
            CALL pmci_anterp_tophat( aerosol_mass(lc)%conc, aerosol_mass_c(:,:,:,lc),               &
                                     kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, ijkfc_s, 's' )
         ENDDO
         IF ( .NOT. salsa_gases_from_chem )  THEN
            DO  lg = 1, ngases_salsa
               CALL pmci_anterp_tophat( salsa_gas(lg)%conc, salsa_gas_c(:,:,:,lg),                  &
                                        kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, ijkfc_s, 's' )
            ENDDO
         ENDIF
      ENDIF

   END SUBROUTINE pmci_anterpolation



   SUBROUTINE pmci_interp_1sto_lr( child_array, parent_array, kct, jfl, jfu, kfl, kfu, edge, var )
!
!--   Interpolation of ghost-node values used as the child-domain boundary
!--   conditions. This subroutine handles the left and right boundaries. 
      IMPLICIT NONE

      INTEGER(iwp), INTENT(IN) ::  kct  !< The parent-grid index in z-direction just below the boundary value node
      
      INTEGER(iwp), DIMENSION(jpsa:jpna), INTENT(IN) ::  jfl  !< Indicates start index of child cells belonging to certain
                                                              !< parent cell - y direction
      INTEGER(iwp), DIMENSION(jpsa:jpna), INTENT(IN) ::  jfu  !< Indicates end index of child cells belonging to certain
                                                              !< parent cell - y direction
      INTEGER(iwp), DIMENSION(0:pg%nz+1), INTENT(IN) ::  kfl  !< Indicates start index of child cells belonging to certain
                                                              !< parent cell - z direction
      INTEGER(iwp), DIMENSION(0:pg%nz+1), INTENT(IN) ::  kfu  !< Indicates end index of child cells belonging to certain
                                                              !< parent cell - z direction

      REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg), INTENT(INOUT) ::  child_array   !< Child-grid array
      REAL(wp), DIMENSION(0:pg%nz+1,jps:jpn,ipl:ipr), INTENT(IN)        ::  parent_array  !< Parent-grid array

      CHARACTER(LEN=1), INTENT(IN) ::  edge                   !< Edge symbol: 'l' or 'r'
      CHARACTER(LEN=1), INTENT(IN) ::  var                    !< Variable symbol: 'u', 'v', 'w' or 's'      
!
!--   Local variables:
      INTEGER(iwp) ::  icb      !< Fixed child-grid index in the x-direction pointing to the node just behind the
                                !< boundary-value node
      INTEGER(iwp) ::  icbc     !< Fixed child-grid index in the x-direction pointing to the boundary-value nodes
      INTEGER(iwp) ::  icbgp    !< Index running over the redundant boundary ghost points in the x-direction
      INTEGER(iwp) ::  ierr     !< MPI error code
      INTEGER(iwp) ::  ipbeg    !< Parent-grid index in the x-direction pointing to the starting point of workarr_lr
                                !< in the parent-grid array
      INTEGER(iwp) ::  ipw      !< Reduced parent-grid index in the x-direction for workarr_lr pointing to
                                !< the boundary ghost node
      INTEGER(iwp) ::  ipwp     !< Reduced parent-grid index in the x-direction for workarr_lr pointing to
                                !< the first prognostic node
      INTEGER(iwp) ::  jc       !< Running child-grid index in the y-direction
      INTEGER(iwp) ::  jp       !< Running parent-grid index in the y-direction
      INTEGER(iwp) ::  kc       !< Running child-grid index in the z-direction
      INTEGER(iwp) ::  kp       !< Running parent-grid index in the z-direction      
      
      REAL(wp) ::  cb           !< Interpolation coefficient for the boundary ghost node  
      REAL(wp) ::  cp           !< Interpolation coefficient for the first prognostic node
      REAL(wp) ::  c_interp_1   !< Value interpolated to the flux point in x direction from the parent-grid data
      REAL(wp) ::  c_interp_2   !< Auxiliary value interpolated  to the flux point in x direction from the parent-grid data
! 
!--   Check which edge is to be handled
      IF ( edge == 'l' )  THEN
!
!--      For u, nxl is a ghost node, but not for the other variables
         IF ( var == 'u' )  THEN
            icbc  = nxl    
            icb   = icbc - 1
            ipw   = 2
            ipwp  = ipw        ! This is redundant when var == 'u' 
            ipbeg = ipl
         ELSE
            icbc  = nxl - 1
            icb   = icbc - 1
            ipw   = 1
            ipwp  = 2
            ipbeg = ipl
         ENDIF        
      ELSEIF ( edge == 'r' )  THEN
         IF ( var == 'u' )  THEN
            icbc  = nxr + 1 
            icb   = icbc + 1
            ipw   = 1
            ipwp  = ipw        ! This is redundant when var == 'u'            
            ipbeg = ipr - 2
         ELSE
            icbc  = nxr + 1
            icb   = icbc + 1
            ipw   = 1
            ipwp  = 0
            ipbeg = ipr - 2
         ENDIF         
      ENDIF
!
!--   Interpolation coefficients
      IF ( interpolation_scheme_lrsn == 1 )  THEN
         cb = 1.0_wp  ! 1st-order upwind 
      ELSE IF ( interpolation_scheme_lrsn == 2 )  THEN
         cb = 0.5_wp  ! 2nd-order central
      ELSE
         cb = 0.5_wp  ! 2nd-order central (default)
      ENDIF         
      cp = 1.0_wp - cb
!
!--   Substitute the necessary parent-grid data to the work array workarr_lr. 
      workarr_lr = 0.0_wp     
      IF ( pdims(2) > 1 )  THEN
         
         IF ( bc_dirichlet_s )  THEN 
            workarr_lr(0:pg%nz+1,jpsw:jpnw-1,0:2) = parent_array(0:pg%nz+1,jpsw:jpnw-1,ipbeg:ipbeg+2)
         ELSE IF ( bc_dirichlet_n )  THEN
            workarr_lr(0:pg%nz+1,jpsw+1:jpnw,0:2) = parent_array(0:pg%nz+1,jpsw+1:jpnw,ipbeg:ipbeg+2)
         ELSE
            workarr_lr(0:pg%nz+1,jpsw+1:jpnw-1,0:2)                                                 &
                 = parent_array(0:pg%nz+1,jpsw+1:jpnw-1,ipbeg:ipbeg+2)
         ENDIF
!
!--      South-north exchange if more than one PE subdomain in the y-direction.
!--      Note that in case of 3-D nesting the south (psouth == MPI_PROC_NULL) 
!--      and north (pnorth == MPI_PROC_NULL) boundaries are not exchanged 
!--      because the nest domain is not cyclic.
!--      From south to north
         CALL MPI_SENDRECV( workarr_lr(0,jpsw+1,0), 1, workarr_lr_exchange_type, psouth,  0,        &
                            workarr_lr(0,jpnw,0), 1, workarr_lr_exchange_type, pnorth,  0, comm2d,  &
                            status, ierr )                              
!                                                                             
!--      From north to south                                                  
         CALL MPI_SENDRECV( workarr_lr(0,jpnw-1,0), 1, workarr_lr_exchange_type, pnorth,  1,        &
                            workarr_lr(0,jpsw,0), 1, workarr_lr_exchange_type, psouth,  1, comm2d,  &
                            status, ierr )                               

      ELSE                                                                    
         workarr_lr(0:pg%nz+1,jpsw:jpnw,0:2) = parent_array(0:pg%nz+1,jpsw:jpnw,ipbeg:ipbeg+2)            
      ENDIF

      IF ( var == 'u' )  THEN

         DO  jp = jpsw, jpnw
            DO  kp = 0, kct 
               
               DO  jc = jfl(jp), jfu(jp)
                  DO  kc = kfl(kp), kfu(kp)
                     child_array(kc,jc,icbc) = workarr_lr(kp,jp,ipw)
                  ENDDO
               ENDDO

            ENDDO
         ENDDO

      ELSE IF ( var == 'v' )  THEN
         
         DO  jp = jpsw, jpnw-1
            DO  kp = 0, kct 
!
!--            First interpolate to the flux point
               c_interp_1 = cb * workarr_lr(kp,jp,ipw)   + cp * workarr_lr(kp,jp,ipwp)
               c_interp_2 = cb * workarr_lr(kp,jp+1,ipw) + cp * workarr_lr(kp,jp+1,ipwp)
!
!--            Use averages of the neighbouring matching grid-line values
               DO  jc = jfl(jp), jfl(jp+1)
                  child_array(kfl(kp):kfu(kp),jc,icbc) = 0.5_wp * ( c_interp_1 + c_interp_2 )
               ENDDO
!
!--            Then set the values along the matching grid-lines  
               IF  ( MOD( jfl(jp), jgsr ) == 0 )  THEN
                  child_array(kfl(kp):kfu(kp),jfl(jp),icbc) = c_interp_1
               ENDIF
            ENDDO
         ENDDO
!
!--      Finally, set the values along the last matching grid-line  
         IF ( MOD( jfl(jpnw), jgsr ) == 0 )  THEN
            DO  kp = 0, kct
               c_interp_1 = cb * workarr_lr(kp,jpnw,ipw) + cp * workarr_lr(kp,jpnw,ipwp)
               child_array(kfl(kp):kfu(kp),jfl(jpnw),icbc) = c_interp_1
            ENDDO
         ENDIF
!
!--      A gap may still remain in some cases if the subdomain size is not 
!--      divisible by the grid-spacing ratio. In such a case, fill the 
!--      gap. Note however, this operation may produce some additional 
!--      momentum conservation error.
         IF ( jfl(jpnw) < nyn )  THEN
            DO  kp = 0, kct
               DO  jc = jfl(jpnw) + 1, nyn
                  child_array(kfl(kp):kfu(kp),jc,icbc) = child_array(kfl(kp):kfu(kp),jfl(jpnw),icbc)
               ENDDO
            ENDDO
         ENDIF

      ELSE IF ( var == 'w' )  THEN

         DO  jp = jpsw, jpnw
            DO  kp = 0, kct + 1   ! It is important to go up to kct+1  
!
!--            Interpolate to the flux point
               c_interp_1 = cb * workarr_lr(kp,jp,ipw) + cp * workarr_lr(kp,jp,ipwp)
!
!--            First substitute only the matching-node values
               child_array(kfu(kp),jfl(jp):jfu(jp),icbc) = c_interp_1
                  
            ENDDO
         ENDDO

         DO  jp = jpsw, jpnw
            DO  kp = 1, kct + 1   ! It is important to go up to kct+1  
!
!--            Then fill up the nodes in between with the averages                  
               DO  kc = kfu(kp-1) + 1, kfu(kp) - 1 
                  child_array(kc,jfl(jp):jfu(jp),icbc) =                                            &
                       0.5_wp * ( child_array(kfu(kp-1),jfl(jp):jfu(jp),icbc)                       &
                       + child_array(kfu(kp),jfl(jp):jfu(jp),icbc) )
               ENDDO
                  
            ENDDO
         ENDDO

      ELSE   ! any scalar
         
         DO  jp = jpsw, jpnw
            DO  kp = 0, kct 
!
!--            Interpolate to the flux point
               c_interp_1 = cb * workarr_lr(kp,jp,ipw) + cp * workarr_lr(kp,jp,ipwp)
               DO  jc = jfl(jp), jfu(jp)
                  DO  kc = kfl(kp), kfu(kp)
                     child_array(kc,jc,icbc) = c_interp_1
                  ENDDO
               ENDDO

            ENDDO
         ENDDO

      ENDIF  ! var
!
!--   Fill up also the redundant 2nd and 3rd ghost-node layers
      IF ( edge == 'l' )  THEN
         DO  icbgp = -nbgp, icb
            child_array(0:nzt+1,nysg:nyng,icbgp) = child_array(0:nzt+1,nysg:nyng,icbc)
         ENDDO
      ELSEIF ( edge == 'r' )  THEN
         DO  icbgp = icb, nx+nbgp
            child_array(0:nzt+1,nysg:nyng,icbgp) = child_array(0:nzt+1,nysg:nyng,icbc)
         ENDDO
      ENDIF

   END SUBROUTINE pmci_interp_1sto_lr



   SUBROUTINE pmci_interp_1sto_sn( child_array, parent_array, kct, ifl, ifu, kfl, kfu, edge, var )
!
!--   Interpolation of ghost-node values used as the child-domain boundary
!--   conditions. This subroutine handles the south and north boundaries. 
      IMPLICIT NONE

      INTEGER(iwp), INTENT(IN) ::  kct  !< The parent-grid index in z-direction just below the boundary value node
      
      INTEGER(iwp), DIMENSION(ipla:ipra), INTENT(IN) ::  ifl  !< Indicates start index of child cells belonging to certain
                                                              !< parent cell - x direction
      INTEGER(iwp), DIMENSION(ipla:ipra), INTENT(IN) ::  ifu  !< Indicates end index of child cells belonging to certain
                                                              !< parent cell - x direction
      INTEGER(iwp), DIMENSION(0:pg%nz+1), INTENT(IN) ::  kfl  !< Indicates start index of child cells belonging to certain
                                                              !< parent cell - z direction
      INTEGER(iwp), DIMENSION(0:pg%nz+1), INTENT(IN) ::  kfu  !< Indicates end index of child cells belonging to certain
                                                              !< parent cell - z direction
      
      REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg), INTENT(INOUT) ::  child_array   !< Child-grid array
      REAL(wp), DIMENSION(0:pg%nz+1,jps:jpn,ipl:ipr), INTENT(IN)        ::  parent_array  !< Parent-grid array

      CHARACTER(LEN=1), INTENT(IN) ::  edge   !< Edge symbol: 's' or 'n' 
      CHARACTER(LEN=1), INTENT(IN) ::  var    !< Variable symbol: 'u', 'v', 'w' or 's'
!
!--   Local variables:      
      INTEGER(iwp) ::  ic       !< Running child-grid index in the x-direction
      INTEGER(iwp) ::  ierr     !< MPI error code
      INTEGER(iwp) ::  ip       !< Running parent-grid index in the x-direction
      INTEGER(iwp) ::  jcb      !< Fixed child-grid index in the y-direction pointing to the node just behind the
                                !< boundary-value node
      INTEGER(iwp) ::  jcbc     !< Fixed child-grid index in the y-direction pointing to the boundary-value nodes
      INTEGER(iwp) ::  jcbgp    !< Index running over the redundant boundary ghost points in y-direction
      INTEGER(iwp) ::  jpbeg    !< Parent-grid index in the y-direction pointing to the starting point of workarr_sn
                                !< in the parent-grid array
      INTEGER(iwp) ::  jpw      !< Reduced parent-grid index in the y-direction for workarr_sn pointing to
                                !< the boundary ghost node 
      INTEGER(iwp) ::  jpwp     !< Reduced parent-grid index in the y-direction for workarr_sn pointing to
                                !< the first prognostic node
      INTEGER(iwp) ::  kc       !< Running child-grid index in the z-direction      
      INTEGER(iwp) ::  kp       !< Running parent-grid index in the z-direction
      REAL(wp) ::  cb           !< Interpolation coefficient for the boundary ghost node  
      REAL(wp) ::  cp           !< Interpolation coefficient for the first prognostic node
      REAL(wp) ::  c_interp_1   !< Value interpolated to the flux point in x direction from the parent-grid data
      REAL(wp) ::  c_interp_2   !< Auxiliary value interpolated  to the flux point in x direction from the parent-grid data

      
!
!--   Check which edge is to be handled: south or north
      IF ( edge == 's' )  THEN
!
!--      For v, nys is a ghost node, but not for the other variables
         IF ( var == 'v' )  THEN
            jcbc  = nys
            jcb   = jcbc - 1
            jpw   = 2
            jpwp  = 2         ! This is redundant when var == 'v'
            jpbeg = jps
         ELSE
            jcbc  = nys - 1
            jcb   = jcbc - 1
            jpw   = 1
            jpwp  = 2
            jpbeg = jps
         ENDIF
      ELSEIF ( edge == 'n' )  THEN
         IF ( var == 'v' )  THEN
            jcbc  = nyn + 1
            jcb   = jcbc + 1
            jpw   = 1
            jpwp  = 0         ! This is redundant when var == 'v'     
            jpbeg = jpn - 2
         ELSE
            jcbc  = nyn + 1
            jcb   = jcbc + 1
            jpw   = 1
            jpwp  = 0
            jpbeg = jpn - 2
         ENDIF
      ENDIF
!
!--   Interpolation coefficients
      IF ( interpolation_scheme_lrsn == 1 )  THEN
         cb = 1.0_wp  ! 1st-order upwind 
      ELSE IF ( interpolation_scheme_lrsn == 2 )  THEN
         cb = 0.5_wp  ! 2nd-order central
      ELSE
         cb = 0.5_wp  ! 2nd-order central (default)
      ENDIF         
      cp = 1.0_wp - cb
!
!--   Substitute the necessary parent-grid data to the work array workarr_sn. 
      workarr_sn = 0.0_wp     
      IF ( pdims(1) > 1 )  THEN

         IF ( bc_dirichlet_l )  THEN
            workarr_sn(0:pg%nz+1,0:2,iplw:iprw-1) = parent_array(0:pg%nz+1,jpbeg:jpbeg+2,iplw:iprw-1)
         ELSE IF ( bc_dirichlet_r )  THEN
            workarr_sn(0:pg%nz+1,0:2,iplw+1:iprw) = parent_array(0:pg%nz+1,jpbeg:jpbeg+2,iplw+1:iprw)
         ELSE
            workarr_sn(0:pg%nz+1,0:2,iplw+1:iprw-1)                                                 &
                 = parent_array(0:pg%nz+1,jpbeg:jpbeg+2,iplw+1:iprw-1)
         ENDIF
!
!--      Left-right exchange if more than one PE subdomain in the x-direction.
!--      Note that in case of 3-D nesting the left (pleft == MPI_PROC_NULL) and 
!--      right (pright == MPI_PROC_NULL) boundaries are not exchanged because 
!--      the nest domain is not cyclic.
!--      From left to right
         CALL MPI_SENDRECV( workarr_sn(0,0,iplw+1), 1, workarr_sn_exchange_type, pleft,   0,        &
                            workarr_sn(0,0,iprw), 1, workarr_sn_exchange_type, pright, 0, comm2d,   &
                            status, ierr )
!                                                                           
!--      From right to left                                                 
         CALL MPI_SENDRECV( workarr_sn(0,0,iprw-1), 1, workarr_sn_exchange_type, pright,  1,        &
                            workarr_sn(0,0,iplw), 1, workarr_sn_exchange_type, pleft, 1, comm2d,    &
                            status, ierr )

      ELSE
         workarr_sn(0:pg%nz+1,0:2,iplw:iprw) = parent_array(0:pg%nz+1,jpbeg:jpbeg+2,iplw:iprw)
      ENDIF

      IF ( var == 'v' )  THEN

         DO  ip = iplw, iprw
            DO  kp = 0, kct 
               
               DO  ic = ifl(ip), ifu(ip)
                  DO  kc = kfl(kp), kfu(kp)
                     child_array(kc,jcbc,ic) = workarr_sn(kp,jpw,ip)
                  ENDDO
               ENDDO

            ENDDO
         ENDDO

      ELSE IF ( var == 'u' )  THEN
         
         DO  ip = iplw, iprw - 1
            DO  kp = 0, kct
!
!--            First interpolate to the flux point
               c_interp_1 = cb * workarr_sn(kp,jpw,ip)   + cp * workarr_sn(kp,jpwp,ip)
               c_interp_2 = cb * workarr_sn(kp,jpw,ip+1) + cp * workarr_sn(kp,jpwp,ip+1)
!
!--            Use averages of the neighbouring matching grid-line values
               DO  ic = ifl(ip), ifl(ip+1)
                  child_array(kfl(kp):kfu(kp),jcbc,ic) = 0.5_wp * ( c_interp_1 + c_interp_2 )
               ENDDO
!
!--            Then set the values along the matching grid-lines  
               IF ( MOD( ifl(ip), igsr ) == 0 )  THEN
                  child_array(kfl(kp):kfu(kp),jcbc,ifl(ip)) = c_interp_1
               ENDIF

            ENDDO
         ENDDO
!
!--      Finally, set the values along the last matching grid-line  
         IF ( MOD( ifl(iprw), igsr ) == 0 )  THEN
            DO  kp = 0, kct
               c_interp_1 = cb * workarr_sn(kp,jpw,iprw) + cp * workarr_sn(kp,jpwp,iprw)
               child_array(kfl(kp):kfu(kp),jcbc,ifl(iprw)) = c_interp_1
            ENDDO
         ENDIF
!
!--      A gap may still remain in some cases if the subdomain size is not 
!--      divisible by the grid-spacing ratio. In such a case, fill the 
!--      gap. Note however, this operation may produce some additional 
!--      momentum conservation error.
         IF ( ifl(iprw) < nxr )  THEN
            DO  kp = 0, kct
               DO  ic = ifl(iprw) + 1, nxr
                  child_array(kfl(kp):kfu(kp),jcbc,ic) = child_array(kfl(kp):kfu(kp),jcbc,ifl(iprw))
               ENDDO
            ENDDO
         ENDIF

      ELSE IF ( var == 'w' )  THEN

         DO  ip = iplw, iprw
            DO  kp = 0, kct + 1   ! It is important to go up to kct+1  
!
!--            Interpolate to the flux point
               c_interp_1 = cb * workarr_sn(kp,jpw,ip) + cp * workarr_sn(kp,jpwp,ip)
!
!--            First substitute only the matching-node values
               child_array(kfu(kp),jcbc,ifl(ip):ifu(ip)) = c_interp_1

            ENDDO
         ENDDO

         DO  ip = iplw, iprw
            DO  kp = 1, kct + 1   ! It is important to go up to kct + 1  
!
!--            Then fill up the nodes in between with the averages
               DO  kc = kfu(kp-1) + 1, kfu(kp) - 1 
                  child_array(kc,jcbc,ifl(ip):ifu(ip)) =                                            &
                       0.5_wp * ( child_array(kfu(kp-1),jcbc,ifl(ip):ifu(ip))                       &
                       + child_array(kfu(kp),jcbc,ifl(ip):ifu(ip)) )
               ENDDO

            ENDDO
         ENDDO

      ELSE   ! Any scalar
         
         DO  ip = iplw, iprw
            DO  kp = 0, kct 
!
!--            Interpolate to the flux point
               c_interp_1 = cb * workarr_sn(kp,jpw,ip) + cp * workarr_sn(kp,jpwp,ip)
               DO  ic = ifl(ip), ifu(ip)
                  DO  kc = kfl(kp), kfu(kp)
                     child_array(kc,jcbc,ic) = c_interp_1
                  ENDDO
               ENDDO

            ENDDO
         ENDDO

      ENDIF  ! var
!
!--   Fill up also the redundant 2nd and 3rd ghost-node layers
      IF ( edge == 's' )  THEN
         DO  jcbgp = -nbgp, jcb
            child_array(0:nzt+1,jcbgp,nxlg:nxrg) = child_array(0:nzt+1,jcbc,nxlg:nxrg)
         ENDDO
      ELSEIF ( edge == 'n' )  THEN
         DO  jcbgp = jcb, ny+nbgp
            child_array(0:nzt+1,jcbgp,nxlg:nxrg) = child_array(0:nzt+1,jcbc,nxlg:nxrg)
         ENDDO
      ENDIF

   END SUBROUTINE pmci_interp_1sto_sn



   SUBROUTINE pmci_interp_1sto_t( child_array, parent_array, kct, ifl, ifu, jfl, jfu, var )
!
!--   Interpolation of ghost-node values used as the child-domain boundary
!--   conditions. This subroutine handles the top boundary. 
      IMPLICIT NONE

      INTEGER(iwp), INTENT(IN) ::  kct  !< The parent-grid index in z-direction just below the boundary value node
      
      INTEGER(iwp), DIMENSION(ipla:ipra), INTENT(IN) ::  ifl  !< Indicates start index of child cells belonging to certain
                                                              !< parent cell - x direction
      INTEGER(iwp), DIMENSION(ipla:ipra), INTENT(IN) ::  ifu  !< Indicates end index of child cells belonging to certain
                                                              !< parent cell - x direction
      INTEGER(iwp), DIMENSION(jpsa:jpna), INTENT(IN) ::  jfl  !< Indicates start index of child cells belonging to certain
                                                              !< parent cell - y direction
      INTEGER(iwp), DIMENSION(jpsa:jpna), INTENT(IN) ::  jfu  !< Indicates end index of child cells belonging to certain
                                                              !< parent cell - y direction

      REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg), INTENT(INOUT) ::  child_array   !< Child-grid array
      REAL(wp), DIMENSION(0:pg%nz+1,jps:jpn,ipl:ipr), INTENT(IN)        ::  parent_array  !< Parent-grid array

      CHARACTER(LEN=1), INTENT(IN) ::  var                    !< Variable symbol: 'u', 'v', 'w' or 's'
!
!--   Local variables:      
      INTEGER(iwp) ::  ic          !< Running child-grid index in the x-direction
      INTEGER(iwp) ::  ierr        !< MPI error code
      INTEGER(iwp) ::  iplc        !< Lower parent-grid index limit in the x-direction for copying parent-grid
                                   !< array data to workarr_t
      INTEGER(iwp) ::  iprc        !< Upper parent-grid index limit in the x-direction for copying parent-grid
                                   !< array data to workarr_t
      INTEGER(iwp) ::  jc          !< Running child-grid index in the y-direction
      INTEGER(iwp) ::  jpsc        !< Lower parent-grid index limit in the y-direction for copying parent-grid
                                   !< array data to workarr_t
      INTEGER(iwp) ::  jpnc        !< Upper parent-grid-index limit in the y-direction for copying parent-grid
                                   !< array data to workarr_t
      INTEGER(iwp) ::  kc          !< Vertical child-grid index fixed to the boundary-value level 
      INTEGER(iwp) ::  ip          !< Running parent-grid index in the x-direction
      INTEGER(iwp) ::  jp          !< Running parent-grid index in the y-direction
      INTEGER(iwp) ::  kpw         !< Reduced parent-grid index in the z-direction for workarr_t pointing to
                                   !< the boundary ghost node
      REAL(wp)     ::  c31         !< Interpolation coefficient for the 3rd-order WS scheme
      REAL(wp)     ::  c32         !< Interpolation coefficient for the 3rd-order WS scheme
      REAL(wp)     ::  c33         !< Interpolation coefficient for the 3rd-order WS scheme
      REAL(wp)     ::  c_interp_1  !< Value interpolated to the flux point in z direction from the parent-grid data
      REAL(wp)     ::  c_interp_2  !< Auxiliary value interpolated to the flux point in z direction from the parent-grid data


      IF ( var == 'w' )  THEN
         kc = nzt
      ELSE
         kc = nzt + 1
      ENDIF
      kpw = 1
!
!--   Interpolation coefficients
      IF ( interpolation_scheme_t == 1 )  THEN
         c31 =  0.0_wp           ! 1st-order upwind
         c32 =  1.0_wp
         c33 =  0.0_wp
      ELSE IF ( interpolation_scheme_t == 2 )  THEN
         c31 =  0.5_wp           ! 2nd-order central
         c32 =  0.5_wp
         c33 =  0.0_wp
      ELSE            
         c31 =  2.0_wp / 6.0_wp  ! 3rd-order WS upwind biased (default)
         c32 =  5.0_wp / 6.0_wp
         c33 = -1.0_wp / 6.0_wp         
      ENDIF         
!
!--   Substitute the necessary parent-grid data to the work array. 
!--   Note that the dimension of workarr_t is (0:2,jpsw:jpnw,iplw:iprw),
!--   And the jc?w and ic?w-index bounds depend on the location of the PE-
!--   subdomain relative to the side boundaries. 
      iplc = iplw + 1
      iprc = iprw - 1      
      jpsc = jpsw + 1
      jpnc = jpnw - 1
      IF ( bc_dirichlet_l )  THEN
         iplc = iplw
      ENDIF
      IF ( bc_dirichlet_r )  THEN
         iprc = iprw
      ENDIF
      IF ( bc_dirichlet_s )  THEN
         jpsc = jpsw
      ENDIF
      IF ( bc_dirichlet_n )  THEN
         jpnc = jpnw
      ENDIF
      workarr_t = 0.0_wp
      workarr_t(0:2,jpsc:jpnc,iplc:iprc) = parent_array(kct:kct+2,jpsc:jpnc,iplc:iprc)
!
!--   Left-right exchange if more than one PE subdomain in the x-direction.
!--   Note that in case of 3-D nesting the left and right boundaries are 
!--   not exchanged because the nest domain is not cyclic.
      IF ( pdims(1) > 1 )  THEN
!
!--      From left to right
         CALL MPI_SENDRECV( workarr_t(0,jpsw,iplw+1), 1, workarr_t_exchange_type_y, pleft, 0,       &
                            workarr_t(0,jpsw,iprw), 1, workarr_t_exchange_type_y, pright, 0,        &
                            comm2d, status, ierr )
!                                                                              
!--      From right to left                                                    
         CALL MPI_SENDRECV( workarr_t(0,jpsw,iprw-1), 1, workarr_t_exchange_type_y, pright, 1,      &
                            workarr_t(0,jpsw,iplw), 1, workarr_t_exchange_type_y, pleft,  1,        &
                            comm2d, status, ierr )                                           
      ENDIF                                                                    
!                                                                              
!--   South-north exchange if more than one PE subdomain in the y-direction.   
!--   Note that in case of 3-D nesting the south and north boundaries are      
!--   not exchanged because the nest domain is not cyclic.                     
      IF ( pdims(2) > 1 )  THEN                                               
!                                                                              
!--      From south to north                                                   
         CALL MPI_SENDRECV( workarr_t(0,jpsw+1,iplw), 1, workarr_t_exchange_type_x, psouth, 2,      &
                            workarr_t(0,jpnw,iplw), 1, workarr_t_exchange_type_x, pnorth, 2,        &
                            comm2d, status, ierr )                                           
!                                                                              
!--      From north to south                                                   
         CALL MPI_SENDRECV( workarr_t(0,jpnw-1,iplw), 1, workarr_t_exchange_type_x, pnorth, 3,      &
                            workarr_t(0,jpsw,iplw), 1, workarr_t_exchange_type_x, psouth, 3,        &
                            comm2d, status, ierr )
      ENDIF

      IF  ( var == 'w' )  THEN
         DO  ip = iplw, iprw
            DO  jp = jpsw, jpnw
 
               DO  ic = ifl(ip), ifu(ip)
                  DO  jc = jfl(jp), jfu(jp)
                     child_array(kc,jc,ic) = workarr_t(kpw,jp,ip)
                  ENDDO
               ENDDO

            ENDDO
         ENDDO

      ELSE IF  ( var == 'u' )  THEN

         DO  ip = iplw, iprw - 1
            DO  jp = jpsw, jpnw
!
!--            First interpolate to the flux point using the 3rd-order WS scheme
               c_interp_1 = c31 * workarr_t(kpw-1,jp,ip)   + c32 * workarr_t(kpw,jp,ip)             &
                          + c33 * workarr_t(kpw+1,jp,ip)
               c_interp_2 = c31 * workarr_t(kpw-1,jp,ip+1) + c32 * workarr_t(kpw,jp,ip+1)           &
                          + c33 * workarr_t(kpw+1,jp,ip+1)
!
!--            Use averages of the neighbouring matching grid-line values
               DO  ic = ifl(ip), ifl(ip+1)
                  child_array(kc,jfl(jp):jfu(jp),ic) = 0.5_wp * ( c_interp_1 + c_interp_2 )
               ENDDO
!
!--            Then set the values along the matching grid-lines  
               IF ( MOD( ifl(ip), igsr ) == 0 )  THEN
!
!--               First interpolate to the flux point using the 3rd-order WS scheme
                  c_interp_1 = c31 * workarr_t(kpw-1,jp,ip) + c32 * workarr_t(kpw,jp,ip)            &
                             + c33 * workarr_t(kpw+1,jp,ip)                  
                  child_array(kc,jfl(jp):jfu(jp),ifl(ip)) = c_interp_1
               ENDIF

            ENDDO
         ENDDO
!
!--      Finally, set the values along the last matching grid-line  
         IF  ( MOD( ifl(iprw), igsr ) == 0 )  THEN
            DO  jp = jpsw, jpnw
!
!--            First interpolate to the flux point using the 3rd-order WS scheme
               c_interp_1 = c31 * workarr_t(kpw-1,jp,iprw) + c32 * workarr_t(kpw,jp,iprw)           &
                          + c33 * workarr_t(kpw+1,jp,iprw)
               child_array(kc,jfl(jp):jfu(jp),ifl(iprw)) = c_interp_1
            ENDDO
         ENDIF
!
!--      A gap may still remain in some cases if the subdomain size is not 
!--      divisible by the grid-spacing ratio. In such a case, fill the 
!--      gap. Note however, this operation may produce some additional 
!--      momentum conservation error.
         IF ( ifl(iprw) < nxr )  THEN
            DO  jp = jpsw, jpnw
               DO  ic = ifl(iprw) + 1, nxr
                  child_array(kc,jfl(jp):jfu(jp),ic) = child_array(kc,jfl(jp):jfu(jp),ifl(iprw))
               ENDDO
            ENDDO
         ENDIF

      ELSE IF  ( var == 'v' )  THEN

         DO  ip = iplw, iprw
            DO  jp = jpsw, jpnw-1
!
!--            First interpolate to the flux point using the 3rd-order WS scheme
               c_interp_1 = c31 * workarr_t(kpw-1,jp,ip)   + c32 * workarr_t(kpw,jp,ip)             &
                          + c33 * workarr_t(kpw+1,jp,ip)
               c_interp_2 = c31 * workarr_t(kpw-1,jp+1,ip) + c32 * workarr_t(kpw,jp+1,ip)           &
                          + c33 * workarr_t(kpw+1,jp+1,ip)
!
!--            Use averages of the neighbouring matching grid-line values
               DO  jc = jfl(jp), jfl(jp+1)          
                  child_array(kc,jc,ifl(ip):ifu(ip)) = 0.5_wp * ( c_interp_1 + c_interp_2 )
               ENDDO
!
!--            Then set the values along the matching grid-lines  
               IF ( MOD( jfl(jp), jgsr ) == 0 )  THEN
                  c_interp_1 = c31 * workarr_t(kpw-1,jp,ip) + c32 * workarr_t(kpw,jp,ip)            &
                             + c33 * workarr_t(kpw+1,jp,ip)
                  child_array(kc,jfl(jp),ifl(ip):ifu(ip)) = c_interp_1
               ENDIF
               
            ENDDO

         ENDDO
!
!--      Finally, set the values along the last matching grid-line
         IF ( MOD( jfl(jpnw), jgsr ) == 0 )  THEN
            DO  ip = iplw, iprw
!
!--            First interpolate to the flux point using the 3rd-order WS scheme
               c_interp_1 = c31 * workarr_t(kpw-1,jpnw,ip) + c32 * workarr_t(kpw,jpnw,ip)           &
                          + c33 * workarr_t(kpw+1,jpnw,ip)
               child_array(kc,jfl(jpnw),ifl(ip):ifu(ip)) = c_interp_1
            ENDDO
         ENDIF
!
!--      A gap may still remain in some cases if the subdomain size is not 
!--      divisible by the grid-spacing ratio. In such a case, fill the 
!--      gap. Note however, this operation may produce some additional 
!--      momentum conservation error.
         IF  ( jfl(jpnw) < nyn )  THEN
            DO  ip = iplw, iprw
               DO  jc = jfl(jpnw)+1, nyn
                  child_array(kc,jc,ifl(ip):ifu(ip)) = child_array(kc,jfl(jpnw),ifl(ip):ifu(ip))
               ENDDO
            ENDDO
         ENDIF

      ELSE  ! any scalar variable

         DO  ip = iplw, iprw
            DO  jp = jpsw, jpnw
!
!--            First interpolate to the flux point using the 3rd-order WS scheme
               c_interp_1 = c31 * workarr_t(kpw-1,jp,ip) + c32 * workarr_t(kpw,jp,ip)               &
                          + c33 * workarr_t(kpw+1,jp,ip)
               DO  ic = ifl(ip), ifu(ip)
                  DO  jc = jfl(jp), jfu(jp)
                     child_array(kc,jc,ic) = c_interp_1
                  ENDDO
               ENDDO

            ENDDO
         ENDDO

      ENDIF  ! var
!
!--   Just fill up the redundant second ghost-node layer in case of var == w.
      IF ( var == 'w' )  THEN
         child_array(nzt+1,:,:) = child_array(nzt,:,:)
      ENDIF

   END SUBROUTINE pmci_interp_1sto_t



   SUBROUTINE pmci_anterp_tophat( child_array, parent_array, kct, ifl, ifu, jfl, jfu, kfl, kfu,     &
                                  ijkfc, var )
!
!--   Anterpolation of internal-node values to be used as the parent-domain
!--   values. This subroutine is based on the first-order numerical
!--   integration of the child-grid values contained within the anterpolation
!--   cell (Clark & Farley, Journal of the Atmospheric Sciences 41(3), 1984).

      IMPLICIT NONE

      INTEGER(iwp), INTENT(IN) ::  kct  !< Top boundary index for anterpolation along z
      
      INTEGER(iwp), DIMENSION(0:pg%nz+1,jpsa:jpna,ipla:ipra), INTENT(IN) ::  ijkfc  !< number of child grid points contributing
                                                                                    !< to a parent grid box
      INTEGER(iwp), DIMENSION(ipla:ipra), INTENT(IN) ::  ifl  !< Indicates start index of child cells belonging to certain
                                                              !< parent cell - x direction
      INTEGER(iwp), DIMENSION(ipla:ipra), INTENT(IN) ::  ifu  !< Indicates end index of child cells belonging to certain
                                                              !< parent cell - x direction
      INTEGER(iwp), DIMENSION(jpsa:jpna), INTENT(IN) ::  jfl  !< Indicates start index of child cells belonging to certain
                                                              !< parent cell - y direction
      INTEGER(iwp), DIMENSION(jpsa:jpna), INTENT(IN) ::  jfu  !< Indicates end index of child cells belonging to certain
                                                              !< parent cell - y direction
      INTEGER(iwp), DIMENSION(0:pg%nz+1), INTENT(IN) ::  kfl  !< Indicates start index of child cells belonging to certain
                                                              !< parent cell - z direction
      INTEGER(iwp), DIMENSION(0:pg%nz+1), INTENT(IN) ::  kfu  !< Indicates end index of child cells belonging to certain
                                                              !< parent cell - z direction

      REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg), INTENT(IN) ::  child_array   !< Child-grid array
      REAL(wp), DIMENSION(0:pg%nz+1,jps:jpn,ipl:ipr), INTENT(INOUT)  ::  parent_array  !< Parent-grid array

      CHARACTER(LEN=*), INTENT(IN) ::  var                   !< Variable symbol: 'u', 'v', 'w' or 's'
!
!--   Local variables:  
      INTEGER(iwp) ::  ic              !< Running index x-direction - child grid
      INTEGER(iwp) ::  ipl_anterp      !< Left boundary index for anterpolation along x
      INTEGER(iwp) ::  ipr_anterp      !< Right boundary index for anterpolation along x
      INTEGER(iwp) ::  jc              !< Running index y-direction - child grid
      INTEGER(iwp) ::  jpn_anterp      !< North boundary index for anterpolation along y
      INTEGER(iwp) ::  jps_anterp      !< South boundary index for anterpolation along y
      INTEGER(iwp) ::  kc              !< Running index z-direction - child grid     
      INTEGER(iwp) ::  kpb_anterp = 0  !< Bottom boundary index for anterpolation along z
      INTEGER(iwp) ::  kpt_anterp      !< Top boundary index for anterpolation along z
      INTEGER(iwp) ::  ip              !< Running index x-direction - parent grid
      INTEGER(iwp) ::  jp              !< Running index y-direction - parent grid
      INTEGER(iwp) ::  kp              !< Running index z-direction - parent grid
      INTEGER(iwp) ::  var_flag        !< bit number used to flag topography on respective grid

      REAL(wp) ::  cellsum       !< sum of respective child cells belonging to parent cell 

!
!--   Define the index bounds ipl_anterp, ipr_anterp, jps_anterp and jpn_anterp.
!--   Note that kcb_anterp is simply zero and kct_anterp depends on kct which enters 
!--   here as a parameter and it is determined in pmci_define_index_mapping.
!--   Note that the grid points directly used also for interpolation (from parent to
!--   child) are always excluded from anterpolation, e.g. anterpolation is maximally
!--   only from 0:kct-1, since kct is directly used for interpolation. Similar restriction is
!--   applied to the lateral boundaries as well. An additional buffer is
!--   also applied (default value for anterpolation_buffer_width = 2) in order
!--   to avoid unphysical accumulation of kinetic energy.
      ipl_anterp = ipl
      ipr_anterp = ipr
      jps_anterp = jps
      jpn_anterp = jpn
      kpb_anterp = 0
      kpt_anterp = kct - 1 - anterpolation_buffer_width

      IF ( nesting_mode /= 'vertical' )  THEN
!
!--      Set the anterpolation buffers on the lateral boundaries
         ipl_anterp = MAX( ipl, iplg + 3 + anterpolation_buffer_width )
         ipr_anterp = MIN( ipr, iprg - 3 - anterpolation_buffer_width )
         jps_anterp = MAX( jps, jpsg + 3 + anterpolation_buffer_width )
         jpn_anterp = MIN( jpn, jpng - 3 - anterpolation_buffer_width )
         
      ENDIF
!
!--   Set masking bit for topography flags 
      IF ( var == 'u' )  THEN 
         var_flag = 1 
      ELSEIF ( var == 'v' )  THEN
         var_flag = 2 
      ELSEIF ( var == 'w' )  THEN
         var_flag = 3
      ELSE
         var_flag = 0
      ENDIF
!
!--   Note that ip, jp, and kp are parent-grid indices and ic,jc, and kc 
!--   are child-grid indices.
      DO  ip = ipl_anterp, ipr_anterp
         DO  jp = jps_anterp, jpn_anterp
!
!--         For simplicity anterpolate within buildings and under elevated
!--         terrain too
            DO  kp = kpb_anterp, kpt_anterp
               cellsum = 0.0_wp
               DO  ic = ifl(ip), ifu(ip)
                  DO  jc = jfl(jp), jfu(jp)
                     DO  kc = kfl(kp), kfu(kp)
                        cellsum = cellsum + MERGE( child_array(kc,jc,ic), 0.0_wp,                   &
                             BTEST( wall_flags_total_0(kc,jc,ic), var_flag ) )
                     ENDDO
                  ENDDO
               ENDDO
!
!--            In case all child grid points are inside topography, i.e. 
!--            ijkfc and cellsum are zero, also parent solution would have 
!--            zero values at that grid point, which may cause problems in 
!--            particular for the temperature. Therefore, in case cellsum is 
!--            zero, keep the parent solution at this point. 
               IF ( ijkfc(kp,jp,ip) /= 0 )  THEN
                  parent_array(kp,jp,ip) = cellsum / REAL( ijkfc(kp,jp,ip), KIND=wp )
               ENDIF

            ENDDO
         ENDDO
      ENDDO

   END SUBROUTINE pmci_anterp_tophat

#endif

 END SUBROUTINE pmci_child_datatrans

! Description:
! ------------
!> Set boundary conditions for the prognostic quantities after interpolation 
!> and anterpolation at upward- and downward facing surfaces.  
!> @todo: add Dirichlet boundary conditions for pot. temperature, humdidity and
!> passive scalar.
!------------------------------------------------------------------------------!
 SUBROUTINE pmci_boundary_conds

#if defined( __parallel )
    IMPLICIT NONE

    INTEGER(iwp) ::  ic  !< Index along x-direction
    INTEGER(iwp) ::  jc  !< Index along y-direction
    INTEGER(iwp) ::  kc  !< Index along z-direction
    INTEGER(iwp) ::  lb  !< Running index for aerosol size bins
    INTEGER(iwp) ::  lc  !< Running index for aerosol mass bins
    INTEGER(iwp) ::  lg  !< Running index for salsa gases
    INTEGER(iwp) ::  m   !< Running index for surface type
    INTEGER(iwp) ::  n   !< Running index for number of chemical species
    

    IF ( debug_output_timestep )  CALL debug_message( 'pmci_boundary_conds', 'start' )
!
!-- Set Dirichlet boundary conditions for horizontal velocity components
    IF ( ibc_uv_b == 0 )  THEN
!
!--    Upward-facing surfaces
       DO  m = 1, bc_h(0)%ns
          ic = bc_h(0)%i(m)            
          jc = bc_h(0)%j(m)
          kc = bc_h(0)%k(m)
          u(kc-1,jc,ic) = 0.0_wp
          v(kc-1,jc,ic) = 0.0_wp
       ENDDO
!
!--    Downward-facing surfaces
       DO  m = 1, bc_h(1)%ns
          ic = bc_h(1)%i(m)            
          jc = bc_h(1)%j(m)
          kc = bc_h(1)%k(m)
          u(kc+1,jc,ic) = 0.0_wp
          v(kc+1,jc,ic) = 0.0_wp
       ENDDO
    ENDIF
!
!-- Set Dirichlet boundary conditions for vertical velocity component
!-- Upward-facing surfaces
    DO  m = 1, bc_h(0)%ns
       ic = bc_h(0)%i(m)            
       jc = bc_h(0)%j(m)
       kc = bc_h(0)%k(m)
       w(kc-1,jc,ic) = 0.0_wp
    ENDDO
!
!-- Downward-facing surfaces
    DO  m = 1, bc_h(1)%ns
       ic = bc_h(1)%i(m)            
       jc = bc_h(1)%j(m)
       kc = bc_h(1)%k(m)
       w(kc+1,jc,ic) = 0.0_wp
    ENDDO
!
!-- Set Neumann boundary conditions for potential temperature
    IF ( .NOT. neutral )  THEN
       IF ( ibc_pt_b == 1 )  THEN
          DO  m = 1, bc_h(0)%ns
             ic = bc_h(0)%i(m)            
             jc = bc_h(0)%j(m)
             kc = bc_h(0)%k(m)
             pt(kc-1,jc,ic) = pt(kc,jc,ic)
          ENDDO
          DO  m = 1, bc_h(1)%ns
             ic = bc_h(1)%i(m)            
             jc = bc_h(1)%j(m)
             kc = bc_h(1)%k(m)
             pt(kc+1,jc,ic) = pt(kc,jc,ic)
          ENDDO   
       ENDIF
    ENDIF
!
!-- Set Neumann boundary conditions for humidity and cloud-physical quantities
    IF ( humidity )  THEN
       IF ( ibc_q_b == 1 )  THEN
          DO  m = 1, bc_h(0)%ns
             ic = bc_h(0)%i(m)            
             jc = bc_h(0)%j(m)
             kc = bc_h(0)%k(m)
             q(kc-1,jc,ic) = q(kc,jc,ic)
          ENDDO  
          DO  m = 1, bc_h(1)%ns
             ic = bc_h(1)%i(m)            
             jc = bc_h(1)%j(m)
             kc = bc_h(1)%k(m)
             q(kc+1,jc,ic) = q(kc,jc,ic)
          ENDDO  
       ENDIF
       IF ( bulk_cloud_model  .AND.  microphysics_morrison )  THEN
          DO  m = 1, bc_h(0)%ns
             ic = bc_h(0)%i(m)            
             jc = bc_h(0)%j(m)
             kc = bc_h(0)%k(m)
             nc(kc-1,jc,ic) = 0.0_wp
             qc(kc-1,jc,ic) = 0.0_wp
          ENDDO  
          DO  m = 1, bc_h(1)%ns
             ic = bc_h(1)%i(m)            
             jc = bc_h(1)%j(m)
             kc = bc_h(1)%k(m)

             nc(kc+1,jc,ic) = 0.0_wp
             qc(kc+1,jc,ic) = 0.0_wp
          ENDDO  
       ENDIF

       IF ( bulk_cloud_model  .AND.  microphysics_seifert )  THEN
          DO  m = 1, bc_h(0)%ns
             ic = bc_h(0)%i(m)            
             jc = bc_h(0)%j(m)
             kc = bc_h(0)%k(m)
             nr(kc-1,jc,ic) = 0.0_wp
             qr(kc-1,jc,ic) = 0.0_wp
          ENDDO  
          DO  m = 1, bc_h(1)%ns
             ic = bc_h(1)%i(m)            
             jc = bc_h(1)%j(m)
             kc = bc_h(1)%k(m)
             nr(kc+1,jc,ic) = 0.0_wp
             qr(kc+1,jc,ic) = 0.0_wp
          ENDDO  
       ENDIF

    ENDIF
!
!-- Set Neumann boundary conditions for passive scalar
    IF ( passive_scalar )  THEN
       IF ( ibc_s_b == 1 )  THEN
          DO  m = 1, bc_h(0)%ns
             ic = bc_h(0)%i(m)            
             jc = bc_h(0)%j(m)
             kc = bc_h(0)%k(m)
             s(kc-1,jc,ic) = s(kc,jc,ic)
          ENDDO 
          DO  m = 1, bc_h(1)%ns
             ic = bc_h(1)%i(m)            
             jc = bc_h(1)%j(m)
             kc = bc_h(1)%k(m)
             s(kc+1,jc,ic) = s(kc,jc,ic)
          ENDDO  
       ENDIF
    ENDIF
!
!-- Set Neumann boundary conditions for chemical species
    IF ( air_chemistry  .AND.  nesting_chem )  THEN
       IF ( ibc_cs_b == 1 )  THEN
          DO  n = 1, nspec
             DO  m = 1, bc_h(0)%ns
                ic = bc_h(0)%i(m)            
                jc = bc_h(0)%j(m)
                kc = bc_h(0)%k(m)
                chem_species(n)%conc(kc-1,jc,ic) = chem_species(n)%conc(kc,jc,ic)
             ENDDO 
             DO  m = 1, bc_h(1)%ns
                ic = bc_h(1)%i(m)            
                jc = bc_h(1)%j(m)
                kc = bc_h(1)%k(m)
                chem_species(n)%conc(kc+1,jc,ic) = chem_species(n)%conc(kc,jc,ic)
             ENDDO
          ENDDO
       ENDIF
    ENDIF 
!
!-- Set Neumann boundary conditions for aerosols and salsa gases
    IF ( salsa  .AND.  nesting_salsa )  THEN
       IF ( ibc_salsa_b == 1 )  THEN
          DO  m = 1, bc_h(0)%ns
             ic = bc_h(0)%i(m)
             jc = bc_h(0)%j(m)
             kc = bc_h(0)%k(m)
             DO  lb = 1, nbins_aerosol
                aerosol_number(lb)%conc(kc-1,jc,ic) = aerosol_number(lb)%conc(kc,jc,ic)
             ENDDO
             DO  lc = 1, nbins_aerosol * ncomponents_mass
                aerosol_mass(lc)%conc(kc-1,jc,ic) = aerosol_mass(lc)%conc(kc,jc,ic)
             ENDDO
             IF ( .NOT. salsa_gases_from_chem )  THEN
                DO  lg = 1, ngases_salsa
                   salsa_gas(lg)%conc(kc-1,jc,ic) = salsa_gas(lg)%conc(kc,jc,ic)
                ENDDO
             ENDIF
          ENDDO
          DO  m = 1, bc_h(1)%ns
             ic = bc_h(1)%i(m)
             jc = bc_h(1)%j(m)
             kc = bc_h(1)%k(m)
             DO  lb = 1, nbins_aerosol
                aerosol_number(lb)%conc(kc+1,jc,ic) = aerosol_number(lb)%conc(kc,jc,ic)
             ENDDO
             DO  lc = 1, nbins_aerosol * ncomponents_mass
                aerosol_mass(lc)%conc(kc+1,jc,ic) = aerosol_mass(lc)%conc(kc,jc,ic)
             ENDDO
             IF ( .NOT. salsa_gases_from_chem )  THEN
                DO  lg = 1, ngases_salsa
                   salsa_gas(lg)%conc(kc+1,jc,ic) = salsa_gas(lg)%conc(kc,jc,ic)
                ENDDO
             ENDIF
          ENDDO
       ENDIF
    ENDIF    

    IF ( debug_output_timestep )  CALL debug_message( 'pmci_boundary_conds', 'end' )

#endif
 END SUBROUTINE pmci_boundary_conds


 
 SUBROUTINE pmci_ensure_nest_mass_conservation

!
!-- Adjust the volume-flow rate through the nested boundaries so that the net volume
!-- flow through all boundaries of the current nest domain becomes zero.
    IMPLICIT NONE

    INTEGER(iwp) ::  i                        !< Running index in the x-direction
    INTEGER(iwp) ::  ierr                     !< MPI error code
    INTEGER(iwp) ::  j                        !< Running index in the y-direction
    INTEGER(iwp) ::  k                        !< Running index in the z-direction
    INTEGER(iwp) ::  n                        !< Running index over the boundary faces: l, r, s, n and t

    REAL(wp) ::  dxdy                         !< Surface area of grid cell top face
    REAL(wp) ::  innor                        !< Inner normal vector of the grid cell face
    REAL(wp) ::  sub_sum                      !< Intermediate sum for reducing the loss of signifigant digits in 2-D summations
    REAL(wp) ::  u_corr_left                  !< Correction added to the left boundary value of u
    REAL(wp) ::  u_corr_right                 !< Correction added to the right boundary value of u
    REAL(wp) ::  v_corr_south                 !< Correction added to the south boundary value of v
    REAL(wp) ::  v_corr_north                 !< Correction added to the north boundary value of v
    REAL(wp) ::  volume_flux_integral         !< Surface integral of volume flux over the domain boundaries
    REAL(wp) ::  volume_flux_local            !< Surface integral of volume flux over the subdomain boundary face
    REAL(wp) ::  w_corr_top                   !< Correction added to the top boundary value of w 

    REAL(wp), DIMENSION(5) ::  volume_flux    !< Surface integral of volume flux over each boundary face of the domain

    
!
!-- Sum up the volume flow through the left boundary
    volume_flux(1) = 0.0_wp
    volume_flux_local = 0.0_wp
    IF ( bc_dirichlet_l )  THEN
       i = 0
       innor = dy
       DO   j = nys, nyn
          sub_sum = 0.0_wp
          DO   k = nzb+1, nzt
             sub_sum = sub_sum + innor * u(k,j,i) * dzw(k)                                          &
                               * MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 1 ) )
          ENDDO
          volume_flux_local = volume_flux_local + sub_sum
       ENDDO
    ENDIF

#if defined( __parallel )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( volume_flux_local, volume_flux(1), 1, MPI_REAL, MPI_SUM, comm2d, ierr )
#else
    volume_flux(1) = volume_flux_local
#endif
!
!-- Sum up the volume flow through the right boundary
    volume_flux(2) = 0.0_wp
    volume_flux_local = 0.0_wp
    IF ( bc_dirichlet_r )  THEN
       i = nx + 1
       innor = -dy
       DO   j = nys, nyn
          sub_sum = 0.0_wp
          DO   k = nzb+1, nzt
             sub_sum = sub_sum + innor * u(k,j,i) * dzw(k)                                          &
                               * MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 1 ) )
          ENDDO
          volume_flux_local = volume_flux_local + sub_sum
       ENDDO
    ENDIF

#if defined( __parallel )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( volume_flux_local, volume_flux(2), 1, MPI_REAL, MPI_SUM, comm2d, ierr )
#else
    volume_flux(2) = volume_flux_local
#endif
!
!-- Sum up the volume flow through the south boundary
    volume_flux(3) = 0.0_wp    
    volume_flux_local = 0.0_wp
    IF ( bc_dirichlet_s )  THEN
       j = 0
       innor = dx
       DO   i = nxl, nxr
          sub_sum = 0.0_wp
          DO   k = nzb+1, nzt
             sub_sum = sub_sum + innor * v(k,j,i) * dzw(k)                                          &
                               * MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 2 ) )
          ENDDO
          volume_flux_local = volume_flux_local + sub_sum
       ENDDO
    ENDIF
    
#if defined( __parallel )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( volume_flux_local, volume_flux(3), 1, MPI_REAL, MPI_SUM, comm2d, ierr )
#else
    volume_flux(3) = volume_flux_local
#endif
!    
!-- Sum up the volume flow through the north boundary
    volume_flux(4) = 0.0_wp
    volume_flux_local = 0.0_wp
    IF ( bc_dirichlet_n )  THEN
       j = ny + 1
       innor = -dx
       DO  i = nxl, nxr
          sub_sum = 0.0_wp
          DO  k = nzb+1, nzt
             sub_sum = sub_sum + innor * v(k,j,i) * dzw(k)                                          &
                               * MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 2 ) )
          ENDDO
          volume_flux_local = volume_flux_local + sub_sum
       ENDDO
    ENDIF

#if defined( __parallel )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( volume_flux_local, volume_flux(4), 1, MPI_REAL, MPI_SUM, comm2d, ierr )
#else
    volume_flux(4) = volume_flux_local
#endif
!
!-- Sum up the volume flow through the top boundary
    volume_flux(5) = 0.0_wp
    volume_flux_local = 0.0_wp
    dxdy = dx * dy
    k = nzt
    DO  i = nxl, nxr
       sub_sum = 0.0_wp
       DO   j = nys, nyn
          sub_sum = sub_sum - w(k,j,i) * dxdy  ! Minus, because the inner unit normal vector is (0,0,-1)
       ENDDO
       volume_flux_local = volume_flux_local + sub_sum
    ENDDO

#if defined( __parallel )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( volume_flux_local, volume_flux(5), 1, MPI_REAL, MPI_SUM, comm2d, ierr )
#else
    volume_flux(5) = volume_flux_local
#endif

    volume_flux_integral = 0.0_wp
    DO  n = 1, 5
       volume_flux_integral = volume_flux_integral + volume_flux(n)
    ENDDO
!    
!-- Correction equally distributed to all nest boundaries, area_total must be used as area.
!-- Note that face_area(6) is the total area (=sum from 1 to 5)
    w_corr_top   = volume_flux_integral / face_area(6)
    u_corr_left  =-w_corr_top
    u_corr_right = w_corr_top
    v_corr_south =-w_corr_top
    v_corr_north = w_corr_top
!!
!!-- Just print out the net volume fluxes through each boundary. Only the root process prints.    
!    if ( myid == 0 )  then       
!       write( 9, "(5(e14.7,2x),4x,e14.7,4x,e12.5,4x,5(e14.7,2x))" )                                 &
!            volume_flux(1), volume_flux(2), volume_flux(3), volume_flux(4), volume_flux(5),         &
!            volume_flux_integral, c_correc,                                                         &
!            u_corr_left, u_corr_right,  v_corr_south, v_corr_north, w_corr_top
!       flush( 9 )
!    endif    
!
!-- Correct the top-boundary value of w
    DO   i = nxl, nxr
       DO   j = nys, nyn
          DO  k = nzt, nzt + 1
             w(k,j,i) = w(k,j,i) + w_corr_top
          ENDDO
       ENDDO
    ENDDO
!
!-- Correct the left-boundary value of u
    IF ( bc_dirichlet_l )  THEN
       DO  i = nxlg, nxl
          DO  j = nys, nyn
             DO  k = nzb + 1, nzt
                u(k,j,i) = u(k,j,i) + u_corr_left                              &
                     * MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 1 ) )
             ENDDO
          ENDDO
       ENDDO
    ENDIF
!
!-- Correct the right-boundary value of u
    IF ( bc_dirichlet_r )  THEN
       DO  i = nxr+1, nxrg
          DO  j = nys, nyn
             DO  k = nzb + 1, nzt
                u(k,j,i) = u(k,j,i) + u_corr_right                              &
                      * MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 1 ) )
             ENDDO
          ENDDO
       ENDDO
    ENDIF
!
!-- Correct the south-boundary value of v
    IF ( bc_dirichlet_s )  THEN
       DO  i = nxl, nxr
          DO  j = nysg, nys
             DO  k = nzb + 1, nzt
                v(k,j,i) = v(k,j,i) + v_corr_south                              &
                      * MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 2 ) )
             ENDDO
          ENDDO
       ENDDO
    ENDIF
!
!-- Correct the north-boundary value of v
    IF ( bc_dirichlet_n )  THEN
       DO  i = nxl, nxr
          DO  j = nyn+1, nyng
             DO  k = nzb + 1, nzt
                v(k,j,i) = v(k,j,i) + v_corr_north                              &
                      * MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), 2 ) )
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    
    
 END SUBROUTINE pmci_ensure_nest_mass_conservation


 
 SUBROUTINE pmci_ensure_nest_mass_conservation_vertical

!
!-- Adjust the volume-flow rate through the top boundary so that the net volume
!-- flow through all boundaries of the current nest domain becomes zero.
    IMPLICIT NONE

    INTEGER(iwp) ::  i                        !< Running index in the x-direction
    INTEGER(iwp) ::  ierr                     !< MPI error code
    INTEGER(iwp) ::  j                        !< Running index in the y-direction
    INTEGER(iwp) ::  k                        !< Running index in the z-direction

    REAL(wp) ::  dxdy                         !< Surface area of grid cell top face
    REAL(wp) ::  sub_sum                      !< Intermediate sum for reducing the loss of signifigant digits in 2-D summations
    REAL(wp) ::  top_area                     !< Top boundary face area
    REAL(wp) ::  volume_flux                  !< Surface integral of volume flux over the top boundary face
    REAL(wp) ::  volume_flux_local            !< Surface integral of volume flux over the subdomain boundary face
    REAL(wp) ::  w_corr_top                   !< Correction added to the top boundary value of w 


    top_area = face_area(5)
!
!-- Sum up the volume flow through the top boundary
    volume_flux = 0.0_wp
    volume_flux_local = 0.0_wp
    dxdy = dx * dy
    k = nzt
    DO  i = nxl, nxr
       sub_sum = 0.0_wp
       DO   j = nys, nyn
          sub_sum = sub_sum - w(k,j,i) * dxdy  ! Minus, because the inner unit normal vector is (0,0,-1)
       ENDDO
       volume_flux_local = volume_flux_local + sub_sum
    ENDDO

#if defined( __parallel )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( volume_flux_local, volume_flux, 1, MPI_REAL, MPI_SUM, comm2d, ierr )
#else
    volume_flux = volume_flux_local
#endif

    w_corr_top   = volume_flux / top_area
!!
!!-- Just print out the net volume fluxes through each boundary. Only the root process prints.    
!    if ( myid == 0 )  then       
!       write( 9, "(5(e14.7,2x),4x,e14.7,4x,e12.5,4x,5(e14.7,2x))" )                                 &
!            volume_flux(1), volume_flux(2), volume_flux(3), volume_flux(4), volume_flux(5),         &
!            volume_flux_integral, c_correc,                                                         &
!            u_corr_left, u_corr_right,  v_corr_south, v_corr_north, w_corr_top
!       flush( 9 )
!    endif    
!
!-- Correct the top-boundary value of w
    DO   i = nxl, nxr
       DO   j = nys, nyn
          DO  k = nzt, nzt + 1
             w(k,j,i) = w(k,j,i) + w_corr_top
          ENDDO
       ENDDO
    ENDDO
    
 END SUBROUTINE pmci_ensure_nest_mass_conservation_vertical

#endif
END MODULE pmc_interface
