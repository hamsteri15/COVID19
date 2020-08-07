!> @file mod_particle_attributes.f90
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
! $Id: mod_particle_attributes.f90 4360 2020-01-07 11:25:50Z suehring $
! Corrected "Former revisions" section
! 
! 4043 2019-06-18 16:59:00Z schwenkel
! Remove min_nr_particle
! 
! 4017 2019-06-06 12:16:46Z schwenkel
! interoperable C datatypes introduced in particle type to avoid compiler warnings
! 
! 3720 2019-02-06 13:19:55Z knoop
! time_prel replaced by last_particle_release_time
! 1359 2014-04-11 17:15:14Z hoffmann
! new module containing all particle related variables
! -dt_sort_particles
!
! Description:
! ------------
!> Definition of variables used to compute particle transport
!------------------------------------------------------------------------------!
MODULE particle_attributes

    USE, INTRINSIC ::  ISO_C_BINDING

    USE kinds

    INTEGER(iwp) ::  dissipation_classes = 10                     !< namelist parameter (see documentation)
    INTEGER(iwp) ::  ibc_par_b                                    !< particle bottom boundary condition dummy
    INTEGER(iwp) ::  ibc_par_lr                                   !< particle left/right boundary condition dummy
    INTEGER(iwp) ::  ibc_par_ns                                   !< particle north/south boundary condition dummy
    INTEGER(iwp) ::  ibc_par_t                                    !< particle top boundary condition dummy
    INTEGER(iwp) ::  number_of_particles = 0                      !< number of particles for each grid box (3d array is saved on prt_count)
    INTEGER(iwp) ::  number_of_particle_groups = 1                !< namelist parameter (see documentation)

    INTEGER(iwp), PARAMETER ::  max_number_of_particle_groups = 10 !< maximum allowed number of particle groups

    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  prt_count  !< 3d array of number of particles of every grid box
    
    LOGICAL ::  particle_advection = .FALSE.              !< parameter to steer the advection of particles
    LOGICAL ::  use_sgs_for_particles = .FALSE.           !< namelist parameter (see documentation)    
    LOGICAL ::  wang_kernel = .FALSE.                     !< flag for collision kernel

    REAL(wp) ::  alloc_factor = 20.0_wp                    !< namelist parameter (see documentation)
    REAL(wp) ::  particle_advection_start = 0.0_wp         !< namelist parameter (see documentation)

!
!-- WARNING: For compatibility of derived types, the BIND attribute is required, and interoperable C
!-- datatypes must be used. These type are hard wired here! So changes in working precision (wp, iwp)
!-- will not affect the particle_type!
!-- The main reason for introducing the interoperable datatypes was to avoid compiler warnings of
!-- the gfortran compiler.
!-- The BIND attribite is required because of C_F_POINTER usage in the pmc particle interface.
    TYPE, BIND(C) ::  particle_type
        REAL(C_DOUBLE) ::  aux1          !< auxiliary multi-purpose feature
        REAL(C_DOUBLE) ::  aux2          !< auxiliary multi-purpose feature
        REAL(C_DOUBLE) ::  radius        !< radius of particle
        REAL(C_DOUBLE) ::  age           !< age of particle
        REAL(C_DOUBLE) ::  age_m         !<
        REAL(C_DOUBLE) ::  dt_sum        !<
        REAL(C_DOUBLE) ::  e_m           !< interpolated sgs tke
        REAL(C_DOUBLE) ::  origin_x      !< origin x-position of particle (changed cyclic bc)
        REAL(C_DOUBLE) ::  origin_y      !< origin y-position of particle (changed cyclic bc)
        REAL(C_DOUBLE) ::  origin_z      !< origin z-position of particle (changed cyclic bc)
        REAL(C_DOUBLE) ::  rvar1         !<
        REAL(C_DOUBLE) ::  rvar2         !<
        REAL(C_DOUBLE) ::  rvar3         !<
        REAL(C_DOUBLE) ::  speed_x       !< speed of particle in x
        REAL(C_DOUBLE) ::  speed_y       !< speed of particle in y
        REAL(C_DOUBLE) ::  speed_z       !< speed of particle in z
        REAL(C_DOUBLE) ::  weight_factor !< weighting factor
        REAL(C_DOUBLE) ::  x             !< x-position
        REAL(C_DOUBLE) ::  y             !< y-position
        REAL(C_DOUBLE) ::  z             !< z-position
        INTEGER(C_INT) ::  class         !< radius class needed for collision
        INTEGER(C_INT) ::  group         !< number of particle group
        INTEGER(C_LONG_LONG) ::  id            !< particle ID (64 bit integer)
        LOGICAL(C_BOOL) ::  particle_mask !< if this parameter is set to false the particle will be deleted
        INTEGER(C_INT) ::  block_nr      !< number for sorting (removable?)
    END TYPE particle_type

    TYPE(particle_type), DIMENSION(:), POINTER ::  particles       !< Particle array for this grid cell
    TYPE(particle_type)                        ::  zero_particle   !< zero particle to avoid weird thinge

    TYPE particle_groups_type
        SEQUENCE
        REAL(wp) ::  density_ratio  !< density ratio of the fluid and the particles
        REAL(wp) ::  radius         !< radius of particle
        REAL(wp) ::  exp_arg        !< exponential term of particle inertia
        REAL(wp) ::  exp_term       !< exponential term of particle inertia
    END TYPE particle_groups_type

    TYPE(particle_groups_type), DIMENSION(max_number_of_particle_groups) ::    &
       particle_groups

    TYPE  grid_particle_def
        INTEGER(iwp), DIMENSION(0:7)               ::  start_index        !< start particle index for current block
        INTEGER(iwp), DIMENSION(0:7)               ::  end_index          !< end particle index for current block
        INTEGER(iwp)                               ::  id_counter         !< particle id counter
        LOGICAL                                    ::  time_loop_done     !< timestep loop for particle advection
        TYPE(particle_type), POINTER, DIMENSION(:) ::  particles          !< Particle array for this grid cell
    END TYPE grid_particle_def

    TYPE(grid_particle_def), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  grid_particles

    TYPE block_offset_def          !<
        INTEGER(iwp) ::  i_off     !<
        INTEGER(iwp) ::  j_off     !<
        INTEGER(iwp) ::  k_off     !<
    END TYPE block_offset_def

    TYPE(block_offset_def), DIMENSION(0:7)         ::  block_offset

    SAVE


END MODULE particle_attributes
