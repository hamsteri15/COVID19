!> @file multi_agent_system_mod.f90
!--------------------------------------------------------------------------------!
! This file is part of PALM-4U.
!
! PALM-4U is free software: you can redistribute it and/or modify it under the 
! terms of the GNU General Public License as published by the Free Software 
! Foundation, either version 3 of the License, or (at your option) any later 
! version.
!
! PALM-4U is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
! A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with
! PALM. If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 2016-2019 Leibniz Universitaet Hannover
!------------------------------------------------------------------------------!
!
! Current revisions:
! ------------------
! 
! 
! Former revisions:
! -----------------
! $Id: multi_agent_system_mod.f90 4444 2020-03-05 15:59:50Z raasch $
! bugfix: cpp-directives for serial mode added
! 
! 4346 2019-12-18 11:55:56Z motisi
! Removed wall_flags_static_0 from USE statements as it's not used within
! the module
! 
! 4329 2019-12-10 15:46:36Z motisi
! Renamed wall_flags_0 to wall_flags_static_0
! 
! 4307 2019-11-26 14:12:36Z maronga
! Activated output of iPT
! 
! 4182 2019-08-22 15:20:23Z scharf
! Corrected "Former revisions" section
! 
! 4168 2019-08-16 13:50:17Z suehring
! Replace function get_topography_top_index by topo_top_ind
! 
! 3987 2019-05-22 09:52:13Z kanani
! Introduce alternative switch for debug output during timestepping
! 
! 3885 2019-04-11 11:29:34Z kanani
! Changes related to global restructuring of location messages and introduction 
! of additional debug messages
! 
! 3876 2019-04-08 18:41:49Z knoop
! replaced nspec by nvar: only variable species should bconsidered, fixed species are not relevant 
! 
! 3766 2019-02-26 16:23:41Z raasch
! save attribute added to local targets to avoid outlive pointer target warning
! 
! 3665 2019-01-10 08:28:24Z raasch
! unused variables removed
! 
! 3159 2018-07-20 11:20:01Z sward
! Initial revision
!
! 
!
! Authors:
! --------
! @author sward
!
!
! Description:
! ------------
!> Multi Agent System for the simulation of pedestrian movement in urban
!> environments
!------------------------------------------------------------------------------!
 MODULE multi_agent_system_mod

    USE, INTRINSIC ::  ISO_C_BINDING

    USE basic_constants_and_equations_mod,                                     &
        ONLY:  pi

    USE control_parameters,                                                    &
        ONLY:  biometeorology,                                                 &
               debug_output_timestep,                                          &
               dt_3d,                                                          &
               dt_write_agent_data,                                            &
               message_string,                                                 &
               time_since_reference_point

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point, log_point_s

    USE grid_variables,                                                        &
        ONLY:  ddx, ddy, dx, dy

    USE indices,                                                               &
        ONLY:  nx, nxl, nxlg, nxr, nxrg, ny, nyn, nyng, nys, nysg, nzb,        &
               topo_top_ind
               
    USE random_function_mod,                                                   &
        ONLY:  random_function

    USE kinds

    USE pegrid

    CHARACTER(LEN=15) ::  bc_mas_lr = 'absorb'  !< left/right boundary condition
    CHARACTER(LEN=15) ::  bc_mas_ns = 'absorb'  !< north/south boundary condition

    INTEGER(iwp) ::  deleted_agents = 0                !< number of deleted agents per time step
    INTEGER(iwp) ::  dim_size_agtnum_manual = 9999999  !< namelist parameter (see documentation)
    INTEGER(iwp) ::  heap_count                        !< number of items in binary heap (for pathfinding)
    INTEGER(iwp) ::  ibc_mas_lr                        !< agent left/right boundary condition dummy
    INTEGER(iwp) ::  ibc_mas_ns                        !< agent north/south boundary condition dummy
!    INTEGER(iwp) ::  ind_pm10 = -9                     !< chemical species index of PM10
!    INTEGER(iwp) ::  ind_pm25 = -9                     !< chemical species index of PM2.5
    INTEGER(iwp) ::  iran_agent = -1234567             !< number for random generator
    INTEGER(iwp) ::  min_nr_agent = 2                  !< namelist parameter (see documentation)
#if defined( __parallel )
    INTEGER(iwp) ::  ghla_count_recv                   !< number of agents in left ghost layer
    INTEGER(iwp) ::  ghna_count_recv                   !< number of agents in north ghost layer
    INTEGER(iwp) ::  ghra_count_recv                   !< number of agents in right ghost layer
    INTEGER(iwp) ::  ghsa_count_recv                   !< number of agents in south ghost layer
    INTEGER(iwp) ::  nr_move_north                     !< number of agts to move north during exchange_horiz
    INTEGER(iwp) ::  nr_move_south                     !< number of agts to move south during exchange_horiz
#endif
    INTEGER(iwp) ::  maximum_number_of_agents = 0      !< maximum number of agents during run
    INTEGER(iwp) ::  number_of_agents = 0              !< number of agents for each grid box (3d array is saved on agt_count)
    INTEGER(iwp) ::  number_of_agent_groups = 1        !< namelist parameter (see documentation)
    INTEGER(iwp) ::  sort_count_mas = 0                !< counter for sorting agents
    INTEGER(iwp) ::  agt_path_size = 15                !< size of agent path array
    INTEGER(iwp) ::  step_dealloc_mas = 100            !< namelist parameter (see documentation)
    INTEGER(iwp) ::  total_number_of_agents            !< total number of agents in the whole model domain

    INTEGER(iwp), PARAMETER ::  NR_2_direction_move = 10000 !< parameter for agent exchange
    INTEGER(iwp), PARAMETER ::  PHASE_INIT    = 1           !< phase parameter
    INTEGER(iwp), PARAMETER ::  PHASE_RELEASE = 2           !< phase parameter

    INTEGER(iwp), PARAMETER ::  max_number_of_agent_groups = 100 !< maximum allowed number of agent groups

    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  agt_count         !< 3d array of number of agents of every grid box
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  s_measure_height  !< k-index(s-grid) for measurement
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  top_top_s         !< k-index of first s-gridpoint above topography
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  top_top_w         !< k-index of first v-gridpoint above topography
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  obstacle_flags    !< flags to identify corners and edges of topography that cannot be crossed by agents

    LOGICAL ::  deallocate_memory_mas = .TRUE.          !< namelist parameter (see documentation)
    LOGICAL ::  dt_3d_reached_mas                       !< flag: agent timestep has reached model timestep
    LOGICAL ::  dt_3d_reached_l_mas                     !< flag: agent timestep has reached model timestep
    LOGICAL ::  agents_active = .FALSE.                 !< flag for agent system
    LOGICAL ::  random_start_position_agents = .TRUE.   !< namelist parameter (see documentation)
    LOGICAL ::  read_agents_from_restartfile = .FALSE.  !< namelist parameter (see documentation)
    LOGICAL ::  agent_own_timestep = .FALSE.            !< namelist parameter (see documentation)

    LOGICAL, DIMENSION(max_number_of_agent_groups) ::  a_rand_target = .FALSE. !< namelist parameter (see documentation)

    REAL(wp) ::  agent_maximum_age = 9999999.9_wp          !< namelist parameter (see documentation)
    REAL(wp) ::  agent_substep_time = 0.0_wp               !< time measurement during one LES timestep
    REAL(wp) ::  alloc_factor_mas = 20.0_wp                !< namelist parameter (see documentation)
    REAL(wp) ::  coll_t_0 = 3.                             !< namelist parameter (see documentation)
    REAL(wp) ::  corner_gate_start = 0.5_wp                !< namelist parameter (see documentation)
    REAL(wp) ::  corner_gate_width = 1.0_wp                !< namelist parameter (see documentation)
    REAL(wp) ::  dim_size_factor_agtnum = 1.0_wp           !< namelist parameter (see documentation)
    REAL(wp) ::  d_sigma_rep_agent                         !< inverse of sigma_rep_agent
    REAL(wp) ::  d_sigma_rep_wall                          !< inverse of sigma_rep_wall
    REAL(wp) ::  d_tau_accel_agent                         !< inverse of tau_accel_agent
    REAL(wp) ::  desired_speed = 1.2_wp                    !< namelist parameter (see documentation)
    REAL(wp) ::  des_sp_sig = .2_wp                        !< namelist parameter (see documentation)
    REAL(wp) ::  dist_target_reached = 2.0_wp              !< distance at which target counts as reached
    REAL(wp) ::  dist_to_int_target = .25_wp               !< namelist parameter (see documentation)
    REAL(wp) ::  dt_agent = 0.02_wp                        !< namelist parameter (see documentation)
    REAL(wp) ::  dt_arel = 9999999.9_wp                    !< namelist parameter (see documentation)
    REAL(wp) ::  end_time_arel = 9999999.9_wp              !< namelist parameter (see documentation)
    REAL(wp) ::  force_x                                   !< dummy value for force on current agent in x-direction
    REAL(wp) ::  force_y                                   !< dummy value for force on current agent in y-direction
    REAL(wp) ::  max_dist_from_path = 0.25_wp              !< distance from current path at which a new path is calculated
    REAL(wp) ::  radius_agent = .25_wp                     !< namelist parameter (see documentation)
    REAL(wp) ::  repuls_agent = 1.5_wp                     !< namelist parameter (see documentation)
    REAL(wp) ::  repuls_wall = 7.0_wp                      !< namelist parameter (see documentation)
    REAL(wp) ::  scan_radius_agent = 3.0_wp                !< namelist parameter (see documentation)
    REAL(wp) ::  scan_radius_wall = 2.0_wp                 !< namelist parameter (see documentation)
    REAL(wp) ::  sigma_rep_agent = 0.3_wp                  !< namelist parameter (see documentation)
    REAL(wp) ::  sigma_rep_wall = 0.1_wp                   !< namelist parameter (see documentation)
    REAL(wp) ::  tau_accel_agent = 0.5_wp                  !< namelist parameter (see documentation)
    REAL(wp) ::  time_arel = 0.0_wp                        !< time for agent release
    REAL(wp) ::  time_write_agent_data = 0.0_wp            !< write agent data at current time on file
    REAL(wp) ::  v_max_agent = 1.3_wp                      !< namelist parameter (see documentation)

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  dummy_path_x  !<  dummy path (x-coordinate)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  dummy_path_y  !<  dummy path (y-coordinate)

    REAL(wp), DIMENSION(max_number_of_agent_groups) ::  adx = 9999999.9_wp  !< namelist parameter (see documentation)
    REAL(wp), DIMENSION(max_number_of_agent_groups) ::  ady = 9999999.9_wp  !< namelist parameter (see documentation)
    REAL(wp), DIMENSION(max_number_of_agent_groups) ::  asl = 9999999.9_wp  !< namelist parameter (see documentation)
    REAL(wp), DIMENSION(max_number_of_agent_groups) ::  asn = 9999999.9_wp  !< namelist parameter (see documentation)
    REAL(wp), DIMENSION(max_number_of_agent_groups) ::  asr = 9999999.9_wp  !< namelist parameter (see documentation)
    REAL(wp), DIMENSION(max_number_of_agent_groups) ::  ass = 9999999.9_wp  !< namelist parameter (see documentation)
    REAL(wp), DIMENSION(max_number_of_agent_groups) ::  at_x = 9999999.9_wp !< namelist parameter (see documentation)
    REAL(wp), DIMENSION(max_number_of_agent_groups) ::  at_y = 9999999.9_wp !< namelist parameter (see documentation)
!
!-- Type for the definition of an agent
    TYPE agent_type
        INTEGER(iwp) ::  block_nr             !< number for sorting
        INTEGER(iwp) ::  group                !< number of agent group
        INTEGER(idp) ::  id                   !< particle ID (64 bit integer)
        INTEGER(iwp) ::  path_counter         !< current target along path (path_x/y)
        LOGICAL      ::  agent_mask           !< if this parameter is set to false the agent will be deleted
        REAL(wp)     ::  age                  !< age of agent
        REAL(wp)     ::  age_m                !< age of agent
        REAL(wp)     ::  dt_sum               !< sum of agents subtimesteps
        REAL(wp)     ::  clo                  !< clothing index
        REAL(wp)     ::  energy_storage       !< energy stored by agent
        REAL(wp)     ::  clothing_temp        !< energy stored by agent
        REAL(wp)     ::  actlev               !< metabolic + work energy of the person
        REAL(wp)     ::  age_years            !< physical age of the person
        REAL(wp)     ::  weight               !< total weight of the person (kg)
        REAL(wp)     ::  height               !< height of the person (m)
        REAL(wp)     ::  work                 !< workload of the agent (W)
        INTEGER(iwp) ::  sex                  !< agents gender: 1 = male, 2 = female
        REAL(wp)     ::  force_x              !< force term x-direction
        REAL(wp)     ::  force_y              !< force term y-direction
        REAL(wp)     ::  origin_x             !< origin x-position of agent
        REAL(wp)     ::  origin_y             !< origin y-position of agent
        REAL(wp)     ::  pm10                 !< PM10 concentration at agent position
        REAL(wp)     ::  pm25                 !< PM25 concentration at agent position
        REAL(wp)     ::  speed_abs            !< absolute value of agent speed
        REAL(wp)     ::  speed_e_x            !< normalized speed of agent in x
        REAL(wp)     ::  speed_e_y            !< normalized speed of agent in y
        REAL(wp)     ::  speed_des            !< agent's desired speed
        REAL(wp)     ::  speed_x              !< speed of agent in x
        REAL(wp)     ::  speed_y              !< speed of agent in y
        REAL(wp)     ::  ipt                  !< instationary thermal index iPT (degree_C)
        REAL(wp)     ::  windspeed            !< absolute value of windspeed at agent position
        REAL(wp)     ::  x                    !< x-position
        REAL(wp)     ::  y                    !< y-position
        REAL(wp)     ::  t                    !< temperature
        REAL(wp)     ::  t_x                  !< x-position
        REAL(wp)     ::  t_y                  !< y-position
        REAL(wp), DIMENSION(0:15) ::  path_x  !< agent path to target (x)
        REAL(wp), DIMENSION(0:15) ::  path_y  !< agent path to target (y)
    END TYPE agent_type

    TYPE(agent_type), DIMENSION(:), POINTER ::  agents               !< Agent array for this grid cell
    TYPE(agent_type)                        ::  zero_agent           !< zero agent to avoid weird thing
#if defined( __parallel )
    TYPE(agent_type), DIMENSION(:), ALLOCATABLE ::  move_also_north  !< for agent exchange between PEs
    TYPE(agent_type), DIMENSION(:), ALLOCATABLE ::  move_also_south  !< for agent exchange between PEs
    TYPE(agent_type), DIMENSION(:), ALLOCATABLE ::  agt_gh_l         !< ghost layer left of pe domain
    TYPE(agent_type), DIMENSION(:), ALLOCATABLE ::  agt_gh_n         !< ghost layer north of pe domain
    TYPE(agent_type), DIMENSION(:), ALLOCATABLE ::  agt_gh_r         !< ghost layer right of pe domain
    TYPE(agent_type), DIMENSION(:), ALLOCATABLE ::  agt_gh_s         !< ghost layer south of pe domain
#endif
!
!-- Type for 2D grid on which agents are stored
    TYPE  grid_agent_def
        INTEGER(iwp), DIMENSION(0:3)            ::  start_index        !< start agent index for current block
        INTEGER(iwp), DIMENSION(0:3)            ::  end_index          !< end agent index for current block
        INTEGER(iwp)                            ::  id_counter         !< agent id counter (removeable?)
        LOGICAL                                 ::  time_loop_done     !< timestep loop for agent advection
        TYPE(agent_type), POINTER, DIMENSION(:) ::  agents             !< Particle array for this grid cell
    END TYPE grid_agent_def

    TYPE(grid_agent_def), DIMENSION(:,:), ALLOCATABLE, TARGET ::  grid_agents !< 2D grid on which agents are stored
!
!-- Item in a priority queue (binary heap)
    TYPE heap_item
       INTEGER(iwp) ::  mesh_id       !< id of the submitted mesh point
       REAL(wp)     ::  priority      !< priority of the mesh point (= distance so far + heuristic to goal)
    END TYPE heap_item

    TYPE(heap_item), DIMENSION(:), ALLOCATABLE ::  queue  !< priority queue realized as binary heap
!
!-- Type for mesh point in visibility graph
    TYPE  mesh_point
        INTEGER(iwp)                            ::  polygon_id          !< Polygon the point belongs to
        INTEGER(iwp)                            ::  vertex_id           !< Vertex in the polygon
        INTEGER(iwp)                            ::  noc                 !< number of connections
        INTEGER(iwp)                            ::  origin_id           !< ID of previous mesh point on path (A*)
        INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  connected_vertices  !< Index of connected vertices
        REAL(wp)                                ::  cost_so_far         !< Cost to reach this mesh point (A*)
        REAL(wp)                                ::  x                   !< x-coordinate
        REAL(wp)                                ::  y                   !< y-coordinate
        REAL(wp)                                ::  x_s                 !< corner shifted outward from building by 1m (x)
        REAL(wp)                                ::  y_s                 !< corner shifted outward from building by 1m (y)
        REAL(wp), DIMENSION(:), ALLOCATABLE     ::  distance_to_vertex  !< Distance to each vertex
    END TYPE mesh_point

    TYPE(mesh_point), DIMENSION(:), ALLOCATABLE ::  mesh     !< navigation mesh
    TYPE(mesh_point), DIMENSION(:), ALLOCATABLE ::  tmp_mesh !< temporary navigation mesh
!
!-- Vertex of a polygon
    TYPE  vertex_type
        LOGICAL               ::  delete  !< Flag to mark vertex for deletion
        REAL(wp)              ::  x       !< x-coordinate
        REAL(wp)              ::  y       !< y-coordinate
    END TYPE vertex_type
!
!-- Polygon containing a number of vertices
    TYPE  polygon_type
        INTEGER(iwp)                                 ::  nov       !< Number of vertices in this polygon
        TYPE(vertex_type), DIMENSION(:), ALLOCATABLE ::  vertices  !< Array of vertices
    END TYPE polygon_type

    TYPE(polygon_type), DIMENSION(:), ALLOCATABLE ::  polygons  !< Building data in polygon form

    SAVE

    PRIVATE
!
!-- Public functions
    PUBLIC mas_init, mas_last_actions, mas_parin, multi_agent_system

!
!-- Public parameters, constants and initial values
    PUBLIC agents_active

    INTERFACE mas_parin
       MODULE PROCEDURE mas_parin
    END INTERFACE mas_parin

    INTERFACE mas_init
       MODULE PROCEDURE mas_init
    END INTERFACE mas_init

    INTERFACE mas_last_actions
       MODULE PROCEDURE mas_last_actions
    END INTERFACE mas_last_actions

    INTERFACE multi_agent_system
       MODULE PROCEDURE multi_agent_system
    END INTERFACE multi_agent_system

    CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Multi Agent System:
!> executes a number of agents sub-timesteps until the model timestep is reached.
!> The agent timestep is usually smaller than the model timestep
!------------------------------------------------------------------------------!
 SUBROUTINE multi_agent_system

    USE biometeorology_mod,                                                    &
        ONLY:  bio_calc_ipt,                                                   &
               bio_calculate_mrt_grid,                                         &
               bio_get_thermal_index_input_ij


    IMPLICIT NONE

    INTEGER(iwp)       ::  i                  !< counter
    INTEGER(iwp)       ::  ie                 !< counter
    INTEGER(iwp)       ::  is                 !< counter
    INTEGER(iwp)       ::  j                  !< counter
    INTEGER(iwp)       ::  je                 !< counter
    INTEGER(iwp)       ::  js                 !< counter
    INTEGER(iwp), SAVE ::  mas_count = 0      !< counts the mas-calls
    INTEGER(iwp)                :: a     !< agent iterator
    !-- local meteorological conditions
    REAL(wp)                    :: tmrt  !< mean radiant temperature        (degree_C)
    REAL(wp)                    :: ta    !< air temperature                 (degree_C)
    REAL(wp)                    :: vp    !< vapour pressure                 (hPa)
    REAL(wp)                    :: v     !< wind speed    (local level)     (m/s)
    REAL(wp)                    :: pair  !< air pressure                    (hPa)


    LOGICAL       ::  first_loop_stride   !< flag for first loop stride of agent sub-timesteps
    LOGICAL, SAVE ::  first_call = .TRUE. !< first call of mas flag for output


    IF ( debug_output_timestep )  CALL debug_message( 'multi_agent_system', 'start' )

    CALL cpu_log( log_point(9), 'mas', 'start' )
!
!-- Initialize variables for the next (sub-) timestep, i.e., for marking
!-- those agents to be deleted after the timestep
    deleted_agents = 0
    agent_substep_time = 0.0_wp
!
!-- If necessary, release new set of agents
    IF ( time_arel >= dt_arel  .AND.  end_time_arel > time_since_reference_point )  THEN

       CALL mas_create_agent(PHASE_RELEASE)
!
!--    The MOD function allows for changes in the output interval with
!--    restart runs.
       time_arel = MOD( time_arel, MAX( dt_arel, dt_3d ) )

    ENDIF

    first_loop_stride = .TRUE.
    grid_agents(:,:)%time_loop_done = .TRUE.
!
!-- Set timestep variable
    IF ( .NOT. agent_own_timestep ) dt_agent = dt_3d
!
!-- Timestep loop for agent transport.
!-- This loop has to be repeated until the transport time of every agent
!-- (within the total domain!) has reached the LES timestep (dt_3d).
!-- Timestep scheme is Euler-forward
    DO
!
!--    Write agent data at current time on file.
       time_write_agent_data = time_write_agent_data + dt_agent
       agent_substep_time    = agent_substep_time    + dt_agent
       IF ( time_write_agent_data >= dt_write_agent_data )  THEN
#if defined( __netcdf )
          IF ( first_loop_stride ) CALL mas_get_prognostic_quantities
          CALL mas_data_output_agents ( first_call )
#else
          WRITE( message_string, * ) 'NetCDF is needed for agent output. ',    &
                                     'Set __netcdf in compiler options'
          CALL message( 'multi_agent_system', 'PA0071', 1, 2, 0, 6, 0 )
#endif
          IF(first_call) first_call = .FALSE.
          time_write_agent_data = time_write_agent_data - dt_write_agent_data
       ENDIF
!
!--    Flag is true by default, will be set to false if an agent has not yet
!--    reached the model timestep
       grid_agents(:,:)%time_loop_done = .TRUE.

!
!--    First part of agent transport:
!--    Evaluate social forces for all agents at current positions
       CALL cpu_log( log_point_s(9), 'mas_social_forces', 'start' )
       DO  i = nxl, nxr
          DO  j = nys, nyn

             number_of_agents = agt_count(j,i)
!
!--          If grid cell is empty, cycle
             IF ( number_of_agents <= 0 ) CYCLE

             agents => grid_agents(j,i)%agents(1:number_of_agents)
!
!--          Evaluation of social forces
             CALL mas_timestep_forces_call(i,j)

          ENDDO
       ENDDO
       CALL cpu_log( log_point_s(9), 'mas_social_forces', 'stop' )
!
!--    Second part of agent transport:
!--    timestep
       CALL cpu_log( log_point_s(16), 'mas_timestep', 'start' )
       DO  i = nxl, nxr
          DO  j = nys, nyn

             number_of_agents = agt_count(j,i)
!
!--          If grid cell is empty, flag must be true
             IF ( number_of_agents <= 0 )  THEN
                grid_agents(j,i)%time_loop_done = .TRUE.
                CYCLE
             ENDIF

             agents => grid_agents(j,i)%agents(1:number_of_agents)

             agents(1:number_of_agents)%agent_mask = .TRUE.
!
!--          Initialize the variable storing the total time that an agent
!--          has advanced within the timestep procedure
             IF ( first_loop_stride )  THEN
                agents(1:number_of_agents)%dt_sum = 0.0_wp
             ENDIF
!
!--          Initialize the switch used for the loop exit condition checked
!--          at the end of this loop. If at least one agent has failed to
!--          reach the LES timestep, this switch will be set false in 
!--          mas_transport.
             dt_3d_reached_l_mas = .TRUE.
!
!--          Timestep
             CALL mas_timestep
!
!--          Delete agents that have been simulated longer than allowed
             CALL mas_boundary_conds( 'max_sim_time' )
!
!--          Delete agents that have reached target area
             CALL mas_boundary_conds( 'target_area' )
!
!---         If not all agents of the actual grid cell have reached the 
!--          LES timestep, this cell has to to another loop iteration. Due to
!--          the fact that agents can move into neighboring grid cell, 
!--          these neighbor cells also have to perform another loop iteration
             IF ( .NOT. dt_3d_reached_l_mas )  THEN
                js = MAX(nys,j-1)
                je = MIN(nyn,j+1)
                is = MAX(nxl,i-1)
                ie = MIN(nxr,i+1)
                grid_agents(js:je,is:ie)%time_loop_done = .FALSE.
             ENDIF

          ENDDO
       ENDDO
       CALL cpu_log( log_point_s(16), 'mas_timestep', 'stop' )

!
!--    Find out, if all agents on every PE have completed the LES timestep
!--    and set the switch corespondingly
       dt_3d_reached_l_mas = ALL(grid_agents(:,:)%time_loop_done)
#if defined( __parallel )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( dt_3d_reached_l_mas, dt_3d_reached_mas, 1, MPI_LOGICAL, &
                           MPI_LAND, comm2d, ierr )
#else
       dt_3d_reached_mas = dt_3d_reached_l_mas
#endif

!
!--    Increment time since last release
       IF ( dt_3d_reached_mas )  time_arel = time_arel + dt_3d

!
!--    Move Agents local to PE to a different grid cell
       CALL cpu_log( log_point_s(18), 'mas_move_exch_sort', 'start' )
       CALL mas_eh_move_agent
!
!--    Horizontal boundary conditions including exchange between subdmains
       CALL mas_eh_exchange_horiz
!
!--    Pack agents (eliminate those marked for deletion),
!--    determine new number of agents
       CALL mas_ps_sort_in_subboxes
       CALL cpu_log( log_point_s(18), 'mas_move_exch_sort', 'stop' )
!
!--    Initialize variables for the next (sub-) timestep, i.e., for marking
!--    those agents to be deleted after the timestep
       deleted_agents = 0

       IF ( biometeorology )  THEN
!
!--       Fill out the MRT 2D grid from appropriate source (RTM, RRTMG,...)
          CALL bio_calculate_mrt_grid ( .FALSE. )
!
!--       Call of human thermal comfort mod (and UV exposure)
          DO  i = nxl, nxr
             DO  j = nys, nyn

                number_of_agents = agt_count(j,i)
!
!--             If grid cell gets empty, cycle
                IF ( number_of_agents <= 0 )  CYCLE

                agents => grid_agents(j,i)%agents(1:number_of_agents)
!
!--             Evaluation of social forces
!                CALL bio_dynamic( i, j )
!
!--             Determine local meteorological conditions
                CALL bio_get_thermal_index_input_ij ( .FALSE., i, j, ta, vp,  &
                                                      v, pair, tmrt )

                DO  a = 1, number_of_agents
!
!--                Calculate instationary thermal indices based on local tmrt

                   CALL bio_calc_ipt ( ta, vp, v, pair, tmrt,                 &
                                       agents(a)%dt_sum,                      &
                                       agents(a)%energy_storage,              &
                                       agents(a)%clothing_temp,               &
                                       agents(a)%clo,                         &
                                       agents(a)%actlev,                      &
                                       agents(a)%age_years,                   &
                                       agents(a)%weight,                      &
                                       agents(a)%height,                      &
                                       agents(a)%work,                        &
                                       agents(a)%sex,                         &
                                       agents(a)%ipt )
                END DO

             ENDDO
          ENDDO
       ENDIF

       IF ( dt_3d_reached_mas )  EXIT

       first_loop_stride = .FALSE.
    ENDDO   ! timestep loop

!
!-- Deallocate unused memory
    IF ( deallocate_memory_mas  .AND.  mas_count == step_dealloc_mas )  THEN
       CALL mas_eh_dealloc_agents_array
       mas_count = 0
    ELSEIF ( deallocate_memory_mas )  THEN
       mas_count = mas_count + 1
    ENDIF

    CALL cpu_log( log_point(9), 'mas', 'stop' )

    IF ( debug_output_timestep )  CALL debug_message( 'multi_agent_system', 'end' )


 END SUBROUTINE multi_agent_system

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculation of the direction vector from each agent to its current
!> intermittent target
!------------------------------------------------------------------------------!
    SUBROUTINE mas_agent_direction

       IMPLICIT NONE

       LOGICAL ::  path_flag !< true if new path must be calculated

       INTEGER(iwp) ::  n  !< loop variable over all agents in a grid box
       INTEGER(iwp) ::  pc !< agent path counter

       REAL(wp) ::  abs_dir         !< length of direction vector (for normalization)
!       REAL(wp) ::  d_curr_target   !< rounding influence expressed as x speed component
!       REAL(wp) ::  d_prev_target   !< rounding influence expressed as x speed component
       REAL(wp) ::  dir_x           !< direction of agent (x)
       REAL(wp) ::  dir_y           !< direction of agent (y)
!       REAL(wp) ::  dist_round = 3. !< distance at which agents start rounding a corner
       REAL(wp) ::  dtit            !< distance to intermittent target
!       REAL(wp) ::  round_fac  = 0.2 !< factor for rounding influence
!       REAL(wp) ::  speed_round_x   !< rounding influence expressed as x speed component
!       REAL(wp) ::  speed_round_y   !< rounding influence expressed as x speed component

!
!--    loop over all agents in the current grid box
       DO n = 1, number_of_agents
          path_flag = .FALSE.
          pc = agents(n)%path_counter
!
!--       If no path was calculated for agent yet, do it
          IF ( pc >= 999 ) THEN
             CALL mas_nav_find_path(n)
             pc = agents(n)%path_counter
!
!--       Check if new path must be calculated and if so, do it
          ELSE
!
!--          Case one: Agent has come close enough to intermittent target.
!--                    -> chose new int target and calculate rest of path if no
!--                       new intermittent targets are left
             dtit = SQRT((agents(n)%x - agents(n)%path_x(pc))**2               &
                       + (agents(n)%y - agents(n)%path_y(pc))**2)
             IF ( dtit < dist_to_int_target ) THEN
                agents(n)%path_counter = agents(n)%path_counter + 1
                pc = agents(n)%path_counter
!
!--             Path counter out of scope (each agent can store a maximum of 15
!--             intermittent targets on the way to her final target); new path
!--             must be calculated
                IF ( pc >= SIZE(agents(n)%path_x) ) THEN
                   path_flag = .TRUE.
                ENDIF
!
!--          Case two: Agent too far from path 
!--                    -> set flag for new path to be calculated
             ELSEIF ( dist_point_to_edge(agents(n)%path_x(pc-1),               &
                                         agents(n)%path_y(pc-1),               &
                                         agents(n)%path_x(pc),                 &
                                         agents(n)%path_y(pc),                 &
                                         agents(n)%x, agents(n)%y)             &
                      > max_dist_from_path )                                   &
             THEN
                path_flag = .TRUE.
             ENDIF
!
!--          If either of the above two cases was true, calculate new path and 
!--          reset 0th path point. This point (the last target the agent had)
!--          is needed for the agents rounding of corners and the calculation
!--          of her deviation from her current path
             IF ( path_flag ) THEN
                CALL mas_nav_find_path(n)
                pc = agents(n)%path_counter
             ENDIF
          ENDIF
!
!--       Normalize direction vector
          abs_dir             = 1.0d-12
          dir_x               = agents(n)%path_x(pc) - agents(n)%x
          dir_y               = agents(n)%path_y(pc) - agents(n)%y
          abs_dir             = SQRT(dir_x**2 + dir_y**2)+1.0d-12
!--         needed later for corner rounding
!           dir_x               = dir_x/abs_dir
!           dir_y               = dir_y/abs_dir
!           dir_x               = dir_x + speed_round_x
!           dir_y               = dir_y + speed_round_y
!           abs_dir             = SQRT(dir_x**2 + dir_y**2)+1.0d-12
          agents(n)%speed_e_x = dir_x/abs_dir
          agents(n)%speed_e_y = dir_y/abs_dir
       ENDDO

!
!-- corner rounding; to be added
!
!--       Calculate direction change due to rounding of corners

!           speed_round_x = 0.
!           speed_round_y = 0.
!           
!           d_curr_target = SQRT( (agents(n)%path_x(pc) - agents(n)%x)**2 +      &
!                                 (agents(n)%path_y(pc) - agents(n)%y)**2 )
!           d_prev_target = SQRT( (agents(n)%path_x(pc-1) - agents(n)%x)**2 +    &
!                                 (agents(n)%path_y(pc-1) - agents(n)%y)**2 )
! !
! !--       Agent is close to next target and that target is not the final one
!           IF ( d_curr_target < dist_round .AND. dist_round <                   &
!                           SQRT( (agents(n)%path_x(pc) - agents(n)%t_x)**2 +    &
!                                 (agents(n)%path_y(pc) - agents(n)%t_y)**2 ) )  &
!           THEN
!              speed_round_x = (agents(n)%path_x(pc+1) - agents(n)%path_x(pc)) / &
!                              ABS( agents(n)%path_x(pc)                         &
!                                 - agents(n)%path_x(pc+1)) * round_fac *        &
!                              SIN( pi/dist_round*d_curr_target )
!              speed_round_y = (agents(n)%path_y(pc+1) - agents(n)%path_y(pc)) / &
!                              ABS( agents(n)%path_y(pc)                         &
!                                 - agents(n)%path_y(pc+1)) * round_fac *        &
!                              SIN( pi/dist_round*d_curr_target )
!           ENDIF
! 
!           IF ( d_prev_target < dist_round ) THEN
!              IF ( agents(n)%path_x(pc) /= agents(n)%path_x(pc+1) ) THEN
!                 speed_round_x = speed_round_x +                                   &
!                                 (agents(n)%path_x(pc) - agents(n)%path_x(pc+1)) / &
!                                 ABS( agents(n)%path_x(pc)                         &
!                                    - agents(n)%path_x(pc+1)) * round_fac *        &
!                                 SIN( pi/dist_round*d_prev_target )
!              ENDIF
!              
!              IF ( agents(n)%path_y(pc) /= agents(n)%path_y(pc+1) ) THEN
!                 speed_round_y = speed_round_y +                                   &
!                              (agents(n)%path_y(pc) - agents(n)%path_y(pc+1)) / &
!                              ABS( agents(n)%path_y(pc)                         &
!                                 - agents(n)%path_y(pc+1)) * round_fac *        &
!                              SIN( pi/dist_round*d_prev_target )
!              ENDIF
             
!           ENDIF


    END SUBROUTINE mas_agent_direction

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Boundary conditions for maximum time, target reached and out of domain
!------------------------------------------------------------------------------!
    SUBROUTINE mas_boundary_conds( location )

       IMPLICIT NONE

       CHARACTER (LEN=*) ::  location !< Identifier

       INTEGER(iwp) ::  n   !< agent number
       INTEGER(iwp) ::  grp !< agent group

       REAL(wp) ::  dist_to_target !< distance to target

       IF ( location == 'max_sim_time' )  THEN

!
!--       Delete agents that have been simulated longer than allowed
          DO  n = 1, number_of_agents

             IF ( agents(n)%age > agent_maximum_age  .AND.                     &
                  agents(n)%agent_mask )                                       &
             THEN
                agents(n)%agent_mask  = .FALSE.
                deleted_agents = deleted_agents + 1
             ENDIF

          ENDDO
       ENDIF

       IF ( location == 'target_area' )  THEN

!
!--       Delete agents that entered target region
          DO  n = 1, number_of_agents
             grp = agents(n)%group
             dist_to_target = SQRT((agents(n)%x-at_x(grp))**2                  &
                                 + (agents(n)%y-at_y(grp))**2)
             IF ( dist_to_target < dist_target_reached ) THEN
                agents(n)%agent_mask  = .FALSE.
                deleted_agents = deleted_agents + 1
             ENDIF

          ENDDO
       ENDIF

    END SUBROUTINE mas_boundary_conds

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Release new agents at their respective sources
!------------------------------------------------------------------------------!
    SUBROUTINE mas_create_agent (phase)

       IMPLICIT  NONE

       INTEGER(iwp) ::  alloc_size  !< relative increase of allocated memory for agents
       INTEGER(iwp) ::  i           !< loop variable ( agent groups )
       INTEGER(iwp) ::  ip          !< index variable along x
       INTEGER(iwp) ::  jp          !< index variable along y
       INTEGER(iwp) ::  loop_stride !< loop variable for initialization
       INTEGER(iwp) ::  n           !< loop variable ( number of agents )
       INTEGER(iwp) ::  new_size    !< new size of allocated memory for agents
       INTEGER(iwp) ::  rn_side     !< index of agent path

       INTEGER(iwp), INTENT(IN) ::  phase       !< mode of inititialization

       INTEGER(iwp), DIMENSION(nysg:nyng,nxlg:nxrg) ::  local_count !< start address of new agent
       INTEGER(iwp), DIMENSION(nysg:nyng,nxlg:nxrg) ::  local_start !< start address of new agent

       LOGICAL ::  first_stride !< flag for initialization

       REAL(wp) ::  pos_x       !< increment for agent position in x
       REAL(wp) ::  pos_y       !< increment for agent position in y
       REAL(wp) ::  rand_contr  !< dummy argument for random position
       REAL(wp) ::  rn_side_dum !< index of agent path

       TYPE(agent_type),TARGET ::  tmp_agent !< temporary agent used for initialization

!
!--    Calculate agent positions and store agent attributes, if
!--    agent is situated on this PE
       DO  loop_stride = 1, 2
          first_stride = (loop_stride == 1)
          IF ( first_stride )   THEN
             local_count = 0           ! count number of agents
          ELSE
             local_count = agt_count   ! Start address of new agents
          ENDIF

          DO  i = 1, number_of_agent_groups

             pos_y = ass(i)

             DO WHILE ( pos_y <= asn(i) )

                IF ( pos_y >= nys * dy  .AND.                                  &
                           pos_y <  ( nyn + 1 ) * dy  )                        &
                THEN

                   pos_x = asl(i)

            xloop: DO WHILE ( pos_x <= asr(i) )

                      IF ( pos_x >= nxl * dx  .AND.                            &
                                 pos_x <  ( nxr + 1) * dx )                    &
                      THEN

                         tmp_agent%agent_mask = .TRUE.
                         tmp_agent%group         = i
                         tmp_agent%id            = 0_idp
                         tmp_agent%block_nr      = -1
                         tmp_agent%path_counter  = 999 !SIZE(tmp_agent%path_x)
                         tmp_agent%age           = 0.0_wp
                         tmp_agent%age_m         = 0.0_wp
                         tmp_agent%dt_sum        = 0.0_wp
                         tmp_agent%clo           = -999.0_wp
                         tmp_agent%energy_storage= 0.0_wp
                         tmp_agent%ipt           = 99999.0_wp
                         tmp_agent%clothing_temp = -999._wp      !< energy stored by agent (W)
                         tmp_agent%actlev        = 134.6862_wp   !< metabolic + work energy of the person
                         tmp_agent%age_years     = 35._wp        !< physical age of the person
                         tmp_agent%weight        = 75._wp        !< total weight of the person (kg)
                         tmp_agent%height        = 1.75_wp       !< height of the person (m)
                         tmp_agent%work          = 134.6862_wp   !< workload of the agent (W)
                         tmp_agent%sex           = 1             !< agents gender: 1 = male, 2 = female
                         tmp_agent%force_x       = 0.0_wp
                         tmp_agent%force_y       = 0.0_wp
                         tmp_agent%origin_x      = pos_x
                         tmp_agent%origin_y      = pos_y
                         tmp_agent%speed_abs     = 0.0_wp
                         tmp_agent%speed_e_x     = 0.0_wp
                         tmp_agent%speed_e_y     = 0.0_wp
                         tmp_agent%speed_des     = random_normal(desired_speed,&
                                                                 des_sp_sig)
                         tmp_agent%speed_x       = 0.0_wp
                         tmp_agent%speed_y       = 0.0_wp
                         tmp_agent%x             = pos_x
                         tmp_agent%y             = pos_y
                         tmp_agent%path_x        = -1.0_wp
                         tmp_agent%path_y        = -1.0_wp
                         tmp_agent%t_x           = - pi
                         tmp_agent%t_y           = - pi
!
!--                      Determine the grid indices of the agent position
                         ip = tmp_agent%x * ddx
                         jp = tmp_agent%y * ddy
!
!--                      Give each agent its target
                         IF ( a_rand_target(i) ) THEN
!
!--                         Agent shall receive random target just outside 
!--                         simulated area
                            rn_side_dum = random_function(iran_agent)
                            rn_side     = FLOOR(4.*rn_side_dum)
                            IF ( rn_side < 2 ) THEN
                               IF ( rn_side == 0 ) THEN
                                  tmp_agent%t_y = -2*dy
                               ELSE
                                  tmp_agent%t_y = (ny+3)*dy
                               ENDIF
                               tmp_agent%t_x = random_function(iran_agent) *   &
                                               (nx+1)*dx
                            ELSE
                               IF ( rn_side == 2 ) THEN
                                  tmp_agent%t_x = -2*dx
                               ELSE
                                  tmp_agent%t_x = (nx+3)*dx
                               ENDIF
                               tmp_agent%t_y = random_function(iran_agent) *   &
                                               (ny+1)*dy
                            ENDIF
!
!--                      Agent gets target of her group
                         ELSE
                            tmp_agent%t_x = at_x(i)
                            tmp_agent%t_y = at_y(i)
                         ENDIF

                         local_count(jp,ip) = local_count(jp,ip) + 1

                         IF ( .NOT. first_stride )  THEN
                            grid_agents(jp,ip)%agents(local_count(jp,ip))      &
                                              = tmp_agent
                         ENDIF

                      ENDIF

                      pos_x = pos_x + adx(i)

                   ENDDO xloop

                ENDIF

                pos_y = pos_y + ady(i)

             ENDDO

          ENDDO

!
!--       Allocate or reallocate agents array to new size
          IF ( first_stride )  THEN
             DO  ip = nxlg, nxrg
                DO  jp = nysg, nyng
                   IF ( phase == PHASE_INIT )  THEN
                      IF ( local_count(jp,ip) > 0 )  THEN
                         alloc_size = MAX( INT( local_count(jp,ip) *           &
                            ( 1.0_wp + alloc_factor_mas / 100.0_wp ) ),        &
                            min_nr_agent )
                      ELSE
                         alloc_size = min_nr_agent
                      ENDIF
                      ALLOCATE(grid_agents(jp,ip)%agents(1:alloc_size))
                      DO  n = 1, alloc_size
                         grid_agents(jp,ip)%agents(n) = zero_agent
                      ENDDO
                   ELSEIF ( phase == PHASE_RELEASE )  THEN
                      IF ( local_count(jp,ip) > 0 )  THEN
                         new_size   = local_count(jp,ip) + agt_count(jp,ip)
                         alloc_size = MAX( INT( new_size * ( 1.0_wp +          &
                            alloc_factor_mas / 100.0_wp ) ), min_nr_agent )
                         IF( alloc_size > SIZE( grid_agents(jp,ip)%agents) )   &
                         THEN
                            CALL mas_eh_realloc_agents_array(ip,jp,alloc_size)
                         ENDIF
                      ENDIF
                   ENDIF
                ENDDO
             ENDDO
          ENDIF

       ENDDO

       local_start = agt_count+1
       agt_count   = local_count

!
!--    Calculate agent IDs
       DO  ip = nxl, nxr
          DO  jp = nys, nyn
             number_of_agents = agt_count(jp,ip)
             IF ( number_of_agents <= 0 )  CYCLE
             agents => grid_agents(jp,ip)%agents(1:number_of_agents)

             DO  n = local_start(jp,ip), number_of_agents  !only new agents 

                agents(n)%id = 10000_idp**2 * grid_agents(jp,ip)%id_counter +  &
                               10000_idp * jp + ip
!
!--             Count the number of agents that have been released before
                grid_agents(jp,ip)%id_counter = grid_agents(jp,ip)%id_counter  &
                                                + 1

             ENDDO

          ENDDO
       ENDDO

!
!--    Add random fluctuation to agent positions.
       IF ( random_start_position_agents )  THEN
          DO  ip = nxl, nxr
             DO  jp = nys, nyn
                number_of_agents = agt_count(jp,ip)
                IF ( number_of_agents <= 0 )  CYCLE
                agents => grid_agents(jp,ip)%agents(1:number_of_agents)
!
!--             Move only new agents. Moreover, limit random fluctuation 
!--             in order to prevent that agents move more than one grid box, 
!--             which would lead to problems concerning agent exchange 
!--             between processors in case adx/ady are larger than dx/dy, 
!--             respectively.  
                DO  n = local_start(jp,ip), number_of_agents
                   IF ( asl(agents(n)%group) /= asr(agents(n)%group) )  THEN
                      rand_contr = ( random_function( iran_agent ) - 0.5_wp ) *&
                                     adx(agents(n)%group)
                      agents(n)%x = agents(n)%x +                        &
                              MERGE( rand_contr, SIGN( dx, rand_contr ),       &
                                     ABS( rand_contr ) < dx                    &
                                   ) 
                   ENDIF
                   IF ( ass(agents(n)%group) /= asn(agents(n)%group) )  THEN
                      rand_contr = ( random_function( iran_agent ) - 0.5_wp ) *&
                                     ady(agents(n)%group)
                      agents(n)%y = agents(n)%y +                        &
                              MERGE( rand_contr, SIGN( dy, rand_contr ),       &
                                     ABS( rand_contr ) < dy )
                   ENDIF
                ENDDO
!
!--             Delete agents that have been simulated longer than allowed
                CALL mas_boundary_conds( 'max_sim_time' )
!
!--             Delete agents that have reached target area
                CALL mas_boundary_conds( 'target_area' )

             ENDDO
          ENDDO
!
!--       Exchange agents between grid cells and processors
          CALL mas_eh_move_agent
          CALL mas_eh_exchange_horiz

       ENDIF
!
!--    In case of random_start_position_agents, delete agents identified by 
!--    mas_eh_exchange_horiz and mas_boundary_conds. Then sort agents into 
!--    blocks, which is needed for a fast interpolation of the LES fields 
!--    on the agent position.
       CALL mas_ps_sort_in_subboxes

!
!--    Determine the current number of agents
       number_of_agents = 0
       DO  ip = nxl, nxr
          DO  jp = nys, nyn
             number_of_agents = number_of_agents + agt_count(jp,ip)
          ENDDO
       ENDDO
!
!--    Calculate the number of agents of the total domain
#if defined( __parallel )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( number_of_agents, total_number_of_agents, 1, &
       MPI_INTEGER, MPI_SUM, comm2d, ierr )
#else
       total_number_of_agents = number_of_agents
#endif

       RETURN

    END SUBROUTINE mas_create_agent

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Creates flags that indicate if a gridbox contains edges or corners. These
!> flags are used for agents to check if obstacles are close to them.
!------------------------------------------------------------------------------!
    SUBROUTINE mas_create_obstacle_flags

       USE arrays_3d,                                                          &
           ONLY:  zw

       IMPLICIT NONE

       INTEGER(iwp) ::  il
       INTEGER(iwp) ::  jl

       ALLOCATE(obstacle_flags(nysg:nyng,nxlg:nxrg))

       obstacle_flags = 0

       DO il = nxlg, nxrg
          DO jl = nysg, nyng
!
!--          Exclude cyclic topography boundary
             IF ( il < 0 .OR. il > nx .OR. jl < 0 .OR. jl > ny ) CYCLE
!
!--          North edge
             IF ( jl < nyng ) THEN
                IF ( ( top_top_s(jl,il) - top_top_s(jl+1,il) ) > 1 .AND.       &
                     ( zw( top_top_w(jl,il) ) -                                &
                       zw( top_top_w(jl+1,il) ) ) > .51_wp )                   &
                THEN
                   obstacle_flags(jl,il) = IBSET( obstacle_flags(jl,il), 0 )
                ENDIF
             ENDIF
!
!--          North right corner
             IF ( jl < nyng .AND. il < nxrg ) THEN
                IF ( ( top_top_s(jl,il) - top_top_s(jl+1,il) )   > 1 .AND.     &
                     ( top_top_s(jl,il) - top_top_s(jl+1,il+1) ) > 1 .AND.     &
                     ( top_top_s(jl,il) - top_top_s(jl,il+1) )   > 1 .AND.     &
                     ( zw( top_top_w(jl,il) ) -                                &
                       zw( top_top_w(jl+1,il+1) ) ) > .51_wp )                 &
                THEN
                   obstacle_flags(jl,il) = IBSET( obstacle_flags(jl,il), 1 )
                ENDIF
             ENDIF
!
!--          Right edge
             IF ( il < nxrg ) THEN
                IF ( ( top_top_s(jl,il) - top_top_s(jl,il+1) ) > 1 .AND.       &
                     ( zw( top_top_w(jl,il) ) -                                &
                       zw( top_top_w(jl,il+1) ) ) > .51_wp )                   &
                THEN
                   obstacle_flags(jl,il) = IBSET( obstacle_flags(jl,il), 2 )
                ENDIF
             ENDIF
!
!--          South right corner
             IF ( jl > nysg .AND. il < nxrg ) THEN
                IF ( ( top_top_s(jl,il) - top_top_s(jl,il+1) )   > 1 .AND.     &
                     ( top_top_s(jl,il) - top_top_s(jl-1,il+1) ) > 1 .AND.     &
                     ( top_top_s(jl,il) - top_top_s(jl-1,il) )   > 1 .AND.     &
                     ( zw(top_top_w(jl,il)) -                                  &
                       zw( top_top_w(jl-1,il+1) ) ) > .51_wp )                 &
                THEN
                   obstacle_flags(jl,il) = IBSET( obstacle_flags(jl,il), 3 )
                ENDIF
             ENDIF
!
!--          South edge
             IF ( jl > nysg ) THEN
                IF ( ( top_top_s(jl,il) - top_top_s(jl-1,il) ) > 1 .AND.       &
                     ( zw( top_top_w(jl,il) ) -                                &
                       zw( top_top_w(jl-1,il) ) ) > .51_wp )                   &
                THEN
                   obstacle_flags(jl,il) = IBSET( obstacle_flags(jl,il), 4 )
                ENDIF
             ENDIF
!
!--          South left corner
             IF ( jl > nysg .AND. il > nxlg ) THEN
                IF ( ( top_top_s(jl,il) - top_top_s(jl-1,il) )   > 1 .AND.     &
                     ( top_top_s(jl,il) - top_top_s(jl-1,il-1) ) > 1 .AND.     &
                     ( top_top_s(jl,il) - top_top_s(jl,il-1) )   > 1 .AND.     &
                     ( zw( top_top_w(jl,il) ) -                                &
                       zw(top_top_w(jl-1,il-1) ) ) > .51_wp )                  &
                THEN
                   obstacle_flags(jl,il) = IBSET( obstacle_flags(jl,il), 5 )
                ENDIF
             ENDIF
!
!--          Left edge
             IF ( il > nxlg ) THEN
                IF ( ( top_top_s(jl,il) - top_top_s(jl,il-1) ) > 1 .AND.       &
                     ( zw(top_top_w(jl,il) ) -                                 &
                       zw(top_top_w(jl,il-1) ) ) > .51_wp )                 &
                THEN
                   obstacle_flags(jl,il) = IBSET( obstacle_flags(jl,il), 6 )
                ENDIF
             ENDIF
!
!--          North left corner
             IF ( jl < nyng .AND. il > nxlg ) THEN
                IF ( ( top_top_s(jl,il) - top_top_s(jl,il-1) )   > 1 .AND.     &
                     ( top_top_s(jl,il) - top_top_s(jl+1,il-1) ) > 1 .AND.     &
                     ( top_top_s(jl,il) - top_top_s(jl+1,il) )   > 1 .AND.     &
                     ( zw( top_top_w(jl,il) ) -                                &
                       zw( top_top_w(jl+1,il-1) ) ) > .51_wp )                 &
                THEN
                   obstacle_flags(jl,il) = IBSET( obstacle_flags(jl,il), 7 )
                ENDIF
             ENDIF

          ENDDO
       ENDDO

    END SUBROUTINE mas_create_obstacle_flags

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Write agent data in netCDF format 
!------------------------------------------------------------------------------!
    SUBROUTINE mas_data_output_agents( ftest )

       USE control_parameters,                                                 &
           ONLY:  agt_time_count, biometeorology, end_time, message_string,    &
                  multi_agent_system_end, multi_agent_system_start

       USE netcdf_interface,                                                   &
           ONLY:  nc_stat, id_set_agt, id_var_time_agt,        &
                  id_var_agt, netcdf_handle_error

       USE pegrid

#if defined( __netcdf )
       USE NETCDF
#endif
       USE mas_global_attributes,                                              &
           ONLY:  dim_size_agtnum

       IMPLICIT NONE

#if defined( __parallel )
       INTEGER(iwp) ::  agt_size !< Agent size in bytes
       INTEGER(iwp) ::  n        !< counter (number of PEs)
       INTEGER(iwp) ::  noa_rcv  !< received number of agents
#endif
       INTEGER(iwp) ::  dummy    !< dummy
       INTEGER(iwp) ::  ii       !< counter (x)
       INTEGER(iwp) ::  ip       !< counter (x)
       INTEGER(iwp) ::  jp       !< counter (y)
       INTEGER(iwp) ::  noa      !< number of agents
       INTEGER(iwp) ::  out_noa  !< number of agents for output

#if defined( __parallel )
       INTEGER(iwp), DIMENSION(0:numprocs-1) ::  noa_arr !< number of agents on each PE
#endif
!
!--    SAVE attribute required to avoid compiler warning about pointer outlive the pointer target
       TYPE(agent_type), DIMENSION(:), ALLOCATABLE, TARGET, SAVE ::  trf_agents !< all agents on current PE
#if defined( __parallel )
       TYPE(agent_type), DIMENSION(:), ALLOCATABLE, TARGET, SAVE ::  out_agents !< all agents in entire domain
#endif

       LOGICAL, INTENT (INOUT) :: ftest

       LOGICAL, SAVE :: agt_dimension_exceeded = .FALSE.

       CALL cpu_log( log_point_s(17), 'mas_data_output', 'start' )
!
!--    Get total number of agents and put all agents on one PE in one array
       noa = 0
       DO  ip = nxl, nxr
          DO  jp = nys, nyn
             noa  = noa  + agt_count(jp,ip)
          ENDDO
       ENDDO
       IF(noa > 0) THEN
          ALLOCATE(trf_agents(1:noa))
          dummy = 1
          DO  ip = nxl, nxr
             DO  jp = nys, nyn
                IF ( agt_count(jp,ip) == 0 ) CYCLE
                agents => grid_agents(jp,ip)%agents(1:agt_count(jp,ip))
                trf_agents(dummy:(dummy-1+agt_count(jp,ip))) = agents
                dummy = dummy + agt_count(jp,ip)
             ENDDO
          ENDDO
       ENDIF
#if defined( __parallel )
!
!--    Gather all agents on PE0 for output
       IF ( myid == 0 )  THEN
          noa_arr(0) = noa
!
!--       Receive data from all other PEs.
          DO  n = 1, numprocs-1
              CALL MPI_RECV( noa_arr(n), 1, MPI_INTEGER,                       &
                              n, 0, comm2d, status, ierr )
          ENDDO
       ELSE
          CALL MPI_SEND( noa, 1, MPI_INTEGER, 0, 0, comm2d, ierr )
       ENDIF
       CALL MPI_BARRIER( comm2d, ierr )
       agt_size = STORAGE_SIZE(zero_agent)/8
       IF ( myid == 0 )  THEN
!
!--       Receive data from all other PEs.
          out_noa = SUM(noa_arr)
          IF ( out_noa > 0 ) THEN
             ALLOCATE( out_agents(1:out_noa) )
             IF ( noa > 0 ) THEN
                out_agents(1:noa) = trf_agents
             ENDIF
             noa_rcv = noa
             DO n = 1, numprocs-1
                IF ( noa_arr(n) > 0 ) THEN
                   CALL MPI_RECV( out_agents(noa_rcv+1), noa_arr(n)*agt_size,  &
                                  MPI_BYTE, n, 0, comm2d, status, ierr )
                   noa_rcv = noa_rcv + noa_arr(n)
                ENDIF
             ENDDO
          ELSE
             ALLOCATE( out_agents(1:2) )
             out_agents = zero_agent
             out_noa    = 2
          ENDIF
       ELSE
          IF ( noa > 0 ) THEN
             CALL MPI_SEND( trf_agents(1), noa*agt_size, MPI_BYTE, 0, 0,       &
                                        comm2d, ierr )
          ENDIF
       ENDIF
!
!--    A barrier has to be set, because otherwise some PEs may
!--    proceed too fast so that PE0 may receive wrong data on
!--    tag 0
       CALL MPI_BARRIER( comm2d, ierr )
#endif
       IF ( myid == 0 ) THEN
#if defined( __parallel )
          agents => out_agents
#else
          agents => trf_agents
#endif

#if defined( __netcdf )
!
!--       Update maximum number of agents
          maximum_number_of_agents = MAX(maximum_number_of_agents, out_noa)
!
!--       Output in netCDF format
          IF ( ftest ) THEN
!
!--          First, define size of agent number dimension from amount of agents
!--          released, release interval, time of agent simulation and max
!--          age of agents
             dim_size_agtnum = MIN( MIN( multi_agent_system_end, end_time )    &
                                       - multi_agent_system_start,             &
                                    agent_maximum_age)

             DO ii = 1, number_of_agent_groups
                dim_size_agtnum = dim_size_agtnum                              &
                                + (FLOOR( ( asr(ii)-asl(ii) ) / adx(ii) ) + 1) &
                                * (FLOOR( ( asn(ii)-ass(ii) ) / ady(ii) ) + 1) &
                                * (FLOOR( dim_size_agtnum / dt_arel )     + 1) &
                                * dim_size_factor_agtnum
                dim_size_agtnum = MIN( dim_size_agtnum, dim_size_agtnum_manual )
             ENDDO
             CALL check_open( 118 )
          ENDIF

!
!--       Update the NetCDF time axis
          agt_time_count = agt_time_count + 1

          IF ( .NOT. agt_dimension_exceeded ) THEN
!
!--          if number of agents to be output exceeds dimension, set flag and
!--          print warning
             IF ( out_noa > dim_size_agtnum ) THEN

                agt_dimension_exceeded = .TRUE.
                WRITE(message_string,'(A,F11.1,2(A,I8))')                      &
                                'Number of agents exceeds agent dimension.' // &
                                '&Starting at time_since_reference_point = ',  &
                                time_since_reference_point,                    &
                                ' s, &data may be missing.'//                  &
                                '&Number of agents:     ', out_noa,            &
                                '&Agent dimension size: ', dim_size_agtnum

                CALL message( 'mas_data_output_agents',                        &
                              'PA0420', 0, 1, 0, 6, 0 )

             ENDIF
          ENDIF

!
!--       reduce number of output agents to dimension size, if necessary
          IF ( agt_dimension_exceeded ) THEN

             out_noa = MIN( out_noa, dim_size_agtnum )

          ENDIF

          nc_stat = NF90_PUT_VAR( id_set_agt, id_var_time_agt,                 &
                                  (/ time_since_reference_point + agent_substep_time /),            &
                                  start = (/ agt_time_count /),                &
                                  count = (/ 1 /) )
          CALL netcdf_handle_error( 'mas_data_output_agents', 1 )

!
!--       Output agent attributes

          nc_stat = NF90_PUT_VAR( id_set_agt, id_var_agt(1), agents%id,        &
                                  start = (/ 1, agt_time_count /),             &
                                  count = (/ out_noa /) )
          CALL netcdf_handle_error( 'mas_data_output_agents', 2 )
 
          nc_stat = NF90_PUT_VAR( id_set_agt, id_var_agt(2), agents%x,         &
                                  start = (/ 1, agt_time_count /),             &
                                  count = (/ out_noa /) )
          CALL netcdf_handle_error( 'mas_data_output_agents', 3 )

          nc_stat = NF90_PUT_VAR( id_set_agt, id_var_agt(3), agents%y,         &
                                  start = (/ 1, agt_time_count /),             &
                                  count = (/ out_noa /) )
          CALL netcdf_handle_error( 'mas_data_output_agents', 4 )

          nc_stat = NF90_PUT_VAR( id_set_agt, id_var_agt(4), agents%windspeed, &
                                  start = (/ 1, agt_time_count /),             &
                                  count = (/ out_noa /) )
          CALL netcdf_handle_error( 'mas_data_output_agents', 5 )

          nc_stat = NF90_PUT_VAR( id_set_agt, id_var_agt(5), agents%t,         &
                                  start = (/ 1, agt_time_count /),             &
                                  count = (/ out_noa /) )
          CALL netcdf_handle_error( 'mas_data_output_agents', 6 )

          nc_stat = NF90_PUT_VAR( id_set_agt, id_var_agt(6), agents%group,     &
                                  start = (/ 1, agt_time_count /),             &
                                  count = (/ out_noa /) )
          CALL netcdf_handle_error( 'mas_data_output_agents', 7 )


          IF ( biometeorology )  THEN
             nc_stat = NF90_PUT_VAR( id_set_agt, id_var_agt(7), agents%ipt,    &
                                     start = (/ 1, agt_time_count /),          &
                                     count = (/ out_noa /) )
             CALL netcdf_handle_error( 'mas_data_output_agents', 8 )  
          ENDIF
          

          
!           nc_stat = NF90_PUT_VAR( id_set_agt, id_var_agt(8), agents%pm10,      &
!                                   start = (/ 1, agt_time_count /),             &
!                                   count = (/ out_noa /) )
!           CALL netcdf_handle_error( 'mas_data_output_agents', 9 )
! 
!           nc_stat = NF90_PUT_VAR( id_set_agt, id_var_agt(9), agents%pm25,      &
!                                   start = (/ 1, agt_time_count /),             &
!                                   count = (/ out_noa /) )
!           CALL netcdf_handle_error( 'mas_data_output_agents', 10 )
! 
! 
!           nc_stat = NF90_PUT_VAR( id_set_agt, id_var_agt(9), agents%uv,        &
!                                   start = (/ 1, agt_time_count /),             &
!                                   count = (/ out_noa /) )
!           CALL netcdf_handle_error( 'mas_data_output_agents', 10 )

          CALL cpu_log( log_point_s(17), 'mas_data_output', 'stop' )


#endif

#if defined( __parallel )
          IF ( ALLOCATED( out_agents ) ) DEALLOCATE( out_agents )
#endif
       ELSE
          CALL cpu_log( log_point_s(17), 'mas_data_output', 'stop' )
       ENDIF

       IF ( ALLOCATED( trf_agents ) ) DEALLOCATE( trf_agents )

    END SUBROUTINE mas_data_output_agents

#if defined( __parallel )
!------------------------------------------------------------------------------!
! Description:
! ------------
!> If an agent moves from one processor to another, this subroutine moves 
!> the corresponding elements from the agent arrays of the old grid cells 
!> to the agent arrays of the new grid cells.
!------------------------------------------------------------------------------!
    SUBROUTINE mas_eh_add_agents_to_gridcell (agent_array)

       IMPLICIT NONE

       INTEGER(iwp) ::  aindex !< dummy argument for new number of agents per grid box
       INTEGER(iwp) ::  ip     !< grid index (x) of agent
       INTEGER(iwp) ::  jp     !< grid index (x) of agent
       INTEGER(iwp) ::  n      !< index variable of agent

       LOGICAL ::  pack_done !< flag to indicate that packing is done

       TYPE(agent_type), DIMENSION(:), INTENT(IN)  ::  agent_array !< new agents in a grid box
       TYPE(agent_type), DIMENSION(:), ALLOCATABLE ::  temp_ns     !< temporary agent array for reallocation

       pack_done     = .FALSE.

       DO n = 1, SIZE(agent_array)

          IF ( .NOT. agent_array(n)%agent_mask )  CYCLE

          ip = agent_array(n)%x * ddx
          jp = agent_array(n)%y * ddy

          IF ( ip >= nxl  .AND.  ip <= nxr  .AND.                              &
               jp >= nys  .AND.  jp <= nyn )  &
          THEN ! agent stays on processor
             number_of_agents = agt_count(jp,ip)
             agents => grid_agents(jp,ip)%agents(1:number_of_agents)

             aindex = agt_count(jp,ip)+1
             IF( aindex > SIZE(grid_agents(jp,ip)%agents) )  THEN
                IF ( pack_done )  THEN
                   CALL mas_eh_realloc_agents_array (ip,jp)
                ELSE
                   CALL mas_ps_pack
                   agt_count(jp,ip) = number_of_agents
                   aindex = agt_count(jp,ip)+1
                   IF ( aindex > SIZE(grid_agents(jp,ip)%agents) )  THEN
                      CALL mas_eh_realloc_agents_array (ip,jp)
                   ENDIF
                   pack_done = .TRUE.
                ENDIF
             ENDIF
             grid_agents(jp,ip)%agents(aindex) = agent_array(n)
             agt_count(jp,ip) = aindex
          ELSE
             IF ( jp <= nys - 1 )  THEN
                nr_move_south = nr_move_south+1
!
!--             Before agent information is swapped to exchange-array, check 
!--             if enough memory is allocated. If required, reallocate exchange
!--             array.
                IF ( nr_move_south > SIZE(move_also_south) )  THEN
!
!--                At first, allocate further temporary array to swap agent 
!--                information.
                   ALLOCATE( temp_ns(SIZE(move_also_south)+NR_2_direction_move))
                   temp_ns(1:nr_move_south-1) = move_also_south                &
                                                (1:nr_move_south-1)
                   DEALLOCATE( move_also_south )
                   ALLOCATE( move_also_south(SIZE(temp_ns)) )
                   move_also_south(1:nr_move_south-1) = temp_ns                &
                                                        (1:nr_move_south-1)
                   DEALLOCATE( temp_ns )

                ENDIF

                move_also_south(nr_move_south) = agent_array(n)

                IF ( jp == -1 )  THEN
!
!--                Apply boundary condition along y
                   IF ( ibc_mas_ns == 0 )  THEN
                      move_also_south(nr_move_south)%y =                       &
                                         move_also_south(nr_move_south)%y      &
                                       + ( ny + 1 ) * dy
                      move_also_south(nr_move_south)%origin_y =                &
                                       move_also_south(nr_move_south)%origin_y &
                                       + ( ny + 1 ) * dy
                   ELSEIF ( ibc_mas_ns == 1 )  THEN
!
!--                   Agent absorption
                      move_also_south(nr_move_south)%agent_mask = .FALSE.
                      deleted_agents = deleted_agents + 1

                   ENDIF
                ENDIF
             ELSEIF ( jp >= nyn+1 )  THEN
                nr_move_north = nr_move_north+1
!
!--             Before agent information is swapped to exchange-array, check 
!--             if enough memory is allocated. If required, reallocate exchange
!--             array.
                IF ( nr_move_north > SIZE(move_also_north) )  THEN
!
!--                At first, allocate further temporary array to swap agent 
!--                information.
                   ALLOCATE( temp_ns(SIZE(move_also_north)+NR_2_direction_move))
                   temp_ns(1:nr_move_north-1) =                                &
                                move_also_south(1:nr_move_north-1)
                   DEALLOCATE( move_also_north )
                   ALLOCATE( move_also_north(SIZE(temp_ns)) )
                   move_also_north(1:nr_move_north-1) =                        &
                               temp_ns(1:nr_move_north-1)
                   DEALLOCATE( temp_ns )

                ENDIF

                move_also_north(nr_move_north) = agent_array(n)
                IF ( jp == ny+1 )  THEN
!
!--                Apply boundary condition along y
                   IF ( ibc_mas_ns == 0 )  THEN

                      move_also_north(nr_move_north)%y =                       &
                         move_also_north(nr_move_north)%y                      &
                       - ( ny + 1 ) * dy
                      move_also_north(nr_move_north)%origin_y =                &
                         move_also_north(nr_move_north)%origin_y               &
                       - ( ny + 1 ) * dy
                   ELSEIF ( ibc_mas_ns == 1 )  THEN
!
!--                   Agent absorption
                      move_also_north(nr_move_north)%agent_mask = .FALSE.
                      deleted_agents = deleted_agents + 1

                   ENDIF
                ENDIF
             ENDIF
          ENDIF
       ENDDO

       RETURN

    END SUBROUTINE mas_eh_add_agents_to_gridcell
#endif


#if defined( __parallel )
!------------------------------------------------------------------------------!
! Description:
! ------------
!> After ghost layer agents have been received from neighboring PEs, this
!> subroutine sorts them into the corresponding grid cells
!------------------------------------------------------------------------------!
    SUBROUTINE mas_eh_add_ghost_agents_to_gridcell (agent_array)

       IMPLICIT NONE

       INTEGER(iwp) ::  ip     !< grid index (x) of agent
       INTEGER(iwp) ::  jp     !< grid index (x) of agent
       INTEGER(iwp) ::  n      !< index variable of agent
       INTEGER(iwp) ::  aindex !< dummy argument for new number of agents per grid box

       LOGICAL ::  pack_done !< flag to indicate that packing is done

       TYPE(agent_type), DIMENSION(:), INTENT(IN)  ::  agent_array !< new agents in a grid box

       pack_done     = .FALSE.

       DO n = 1, SIZE(agent_array)

          IF ( .NOT. agent_array(n)%agent_mask )  CYCLE

          ip = agent_array(n)%x * ddx
          jp = agent_array(n)%y * ddy

          IF ( ip < nxl  .OR.  ip > nxr  .OR.  jp < nys  .OR.  jp > nyn ) THEN
             number_of_agents = agt_count(jp,ip)
             agents => grid_agents(jp,ip)%agents(1:number_of_agents)

             aindex = agt_count(jp,ip)+1
             IF( aindex > SIZE(grid_agents(jp,ip)%agents) )  THEN
                IF ( pack_done )  THEN
                   CALL mas_eh_realloc_agents_array (ip,jp)
                ELSE
                   CALL mas_ps_pack
                   agt_count(jp,ip) = number_of_agents
                   aindex = agt_count(jp,ip)+1
                   IF ( aindex > SIZE(grid_agents(jp,ip)%agents) )  THEN
                      CALL mas_eh_realloc_agents_array (ip,jp)
                   ENDIF
                   pack_done = .TRUE.
                ENDIF
             ENDIF
             grid_agents(jp,ip)%agents(aindex) = agent_array(n)
             agt_count(jp,ip) = aindex
          ENDIF
       ENDDO
    END SUBROUTINE mas_eh_add_ghost_agents_to_gridcell
#endif

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Resizing of agent arrays
!------------------------------------------------------------------------------!
    SUBROUTINE mas_eh_dealloc_agents_array

       IMPLICIT NONE

       INTEGER(iwp) ::  i         !< grid index (x) of agent
       INTEGER(iwp) ::  j         !< grid index (y) of agent
       INTEGER(iwp) ::  old_size  !< old array size
       INTEGER(iwp) ::  new_size  !< new array size
       INTEGER(iwp) ::  noa       !< number of agents

       LOGICAL ::  dealloc  !< flag that indicates if reallocation is necessary

       TYPE(agent_type), DIMENSION(10) ::  tmp_agents_s !< temporary static agent array

       TYPE(agent_type), DIMENSION(:), ALLOCATABLE ::  tmp_agents_d !< temporary dynamic agent array

       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
!
!--          Determine number of active agents
             noa = agt_count(j,i)
!
!--          Determine allocated memory size
             old_size = SIZE( grid_agents(j,i)%agents )
!
!--          Check for large unused memory
             dealloc = ( ( noa < min_nr_agent .AND. old_size  > min_nr_agent ) &
                    .OR. ( noa > min_nr_agent .AND. old_size - noa *           &
                         ( 1.0_wp + 0.01_wp * alloc_factor_mas ) > 0.0_wp ) )
!
!--          If large unused memory was found, resize the corresponding array
             IF ( dealloc )  THEN
                IF ( noa < min_nr_agent )  THEN
                   new_size = min_nr_agent
                ELSE 
                   new_size = INT( noa * ( 1.0_wp +                            &
                                            0.01_wp * alloc_factor_mas ) )
                ENDIF

                IF ( noa <= 10 )  THEN

                   tmp_agents_s(1:noa) = grid_agents(j,i)%agents(1:noa)

                   DEALLOCATE(grid_agents(j,i)%agents)
                   ALLOCATE(grid_agents(j,i)%agents(1:new_size))

                   grid_agents(j,i)%agents(1:noa)          = tmp_agents_s(1:noa)
                   grid_agents(j,i)%agents(noa+1:new_size) = zero_agent

                ELSE

                   ALLOCATE(tmp_agents_d(noa))
                   tmp_agents_d(1:noa) = grid_agents(j,i)%agents(1:noa)

                   DEALLOCATE(grid_agents(j,i)%agents)
                   ALLOCATE(grid_agents(j,i)%agents(new_size))

                   grid_agents(j,i)%agents(1:noa)          = tmp_agents_d(1:noa)
                   grid_agents(j,i)%agents(noa+1:new_size) = zero_agent

                   DEALLOCATE(tmp_agents_d)

                ENDIF

             ENDIF
          ENDDO
       ENDDO

    END SUBROUTINE mas_eh_dealloc_agents_array

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Exchange between subdomains.
!> As soon as one agent has moved beyond the boundary of the domain, it
!> is included in the relevant transfer arrays and marked for subsequent
!> deletion on this PE.
!> First sweep for crossings in x direction. Find out first the number of
!> agents to be transferred and allocate temporary arrays needed to store
!> them.
!> For a one-dimensional decomposition along y, no transfer is necessary,
!> because the agent remains on the PE, but the agent coordinate has to
!> be adjusted.
!------------------------------------------------------------------------------!
    SUBROUTINE mas_eh_exchange_horiz

       IMPLICIT NONE

       INTEGER(iwp) ::  ip               !< index variable along x
       INTEGER(iwp) ::  jp               !< index variable along y
       INTEGER(iwp) ::  n                !< agent index variable

#if defined( __parallel )

       INTEGER(iwp) ::  i                !< grid index (x) of agent positition
       INTEGER(iwp) ::  j                !< grid index (y) of agent positition
       INTEGER(iwp) ::  par_size         !< Agent size in bytes

       INTEGER(iwp) ::  trla_count       !< number of agents send to left PE
       INTEGER(iwp) ::  trla_count_recv  !< number of agents receive from right PE
       INTEGER(iwp) ::  trna_count       !< number of agents send to north PE
       INTEGER(iwp) ::  trna_count_recv  !< number of agents receive from south PE
       INTEGER(iwp) ::  trra_count       !< number of agents send to right PE
       INTEGER(iwp) ::  trra_count_recv  !< number of agents receive from left PE
       INTEGER(iwp) ::  trsa_count       !< number of agents send to south PE
       INTEGER(iwp) ::  trsa_count_recv  !< number of agents receive from north PE

       TYPE(agent_type), DIMENSION(:), ALLOCATABLE ::  rvla  !< agents received from right PE
       TYPE(agent_type), DIMENSION(:), ALLOCATABLE ::  rvna  !< agents received from south PE
       TYPE(agent_type), DIMENSION(:), ALLOCATABLE ::  rvra  !< agents received from left PE
       TYPE(agent_type), DIMENSION(:), ALLOCATABLE ::  rvsa  !< agents received from north PE
       TYPE(agent_type), DIMENSION(:), ALLOCATABLE ::  trla  !< agents send to left PE
       TYPE(agent_type), DIMENSION(:), ALLOCATABLE ::  trna  !< agents send to north PE
       TYPE(agent_type), DIMENSION(:), ALLOCATABLE ::  trra  !< agents send to right PE
       TYPE(agent_type), DIMENSION(:), ALLOCATABLE ::  trsa  !< agents send to south PE

!
!--    Exchange between subdomains.
!--    As soon as one agent has moved beyond the boundary of the domain, it
!--    is included in the relevant transfer arrays and marked for subsequent
!--    deletion on this PE.
!--    First sweep for crossings in x direction. Find out first the number of
!--    agents to be transferred and allocate temporary arrays needed to store
!--    them.
!--    For a one-dimensional decomposition along y, no transfer is necessary,
!--    because the agent remains on the PE, but the agent coordinate has to
!--    be adjusted.
       trla_count  = 0
       trra_count  = 0

       trla_count_recv   = 0
       trra_count_recv   = 0

       IF ( pdims(1) /= 1 )  THEN
!
!--       First calculate the storage necessary for sending and receiving the data.
!--       Compute only first (nxl) and last (nxr) loop iterration.
          DO  ip = nxl, nxr, nxr - nxl
             DO  jp = nys, nyn

                number_of_agents = agt_count(jp,ip)
                IF ( number_of_agents <= 0 )  CYCLE
                agents => grid_agents(jp,ip)%agents(1:number_of_agents)
                DO  n = 1, number_of_agents
                   IF ( agents(n)%agent_mask )  THEN
                      i = agents(n)%x * ddx
!
!--                   Above calculation does not work for indices less than zero
                      IF ( agents(n)%x < 0.0_wp )  i = -1

                      IF ( i < nxl )  THEN
                         trla_count = trla_count + 1
                      ELSEIF ( i > nxr )  THEN
                         trra_count = trra_count + 1
                      ENDIF
                   ENDIF
                ENDDO

             ENDDO
          ENDDO

          IF ( trla_count  == 0 )  trla_count  = 1
          IF ( trra_count  == 0 )  trra_count  = 1

          ALLOCATE( trla(trla_count), trra(trra_count) )

          trla = zero_agent
          trra = zero_agent

          trla_count  = 0
          trra_count  = 0

       ENDIF
!
!--    Compute only first (nxl) and last (nxr) loop iterration
       DO  ip = nxl, nxr, nxr-nxl
          DO  jp = nys, nyn
             number_of_agents = agt_count(jp,ip)
             IF ( number_of_agents <= 0 ) CYCLE
             agents => grid_agents(jp,ip)%agents(1:number_of_agents)
             DO  n = 1, number_of_agents
!
!--             Only those agents that have not been marked as 'deleted' may
!--             be moved.
                IF ( agents(n)%agent_mask )  THEN

                   i = agents(n)%x * ddx
!
!--                Above calculation does not work for indices less than zero
                   IF ( agents(n)%x < 0.0_wp )  i = -1

                   IF ( i <  nxl )  THEN
                      IF ( i < 0 )  THEN
!
!--                      Apply boundary condition along x
                         IF ( ibc_mas_lr == 0 )  THEN
!
!--                         Cyclic condition
                            IF ( pdims(1) == 1 )  THEN
                               agents(n)%x        = ( nx + 1 ) * dx +          &
                                                    agents(n)%x
                               agents(n)%origin_x = ( nx + 1 ) * dx +          &
                                                    agents(n)%origin_x
                            ELSE
                               trla_count         = trla_count + 1
                               trla(trla_count)   = agents(n)
                               trla(trla_count)%x = ( nx + 1 ) * dx +          &
                                                    trla(trla_count)%x
                               trla(trla_count)%origin_x =                     &
                                                trla(trla_count)%origin_x +    &
                                                ( nx + 1 ) * dx
                               agents(n)%agent_mask  = .FALSE.
                               deleted_agents = deleted_agents + 1

                               IF ( trla(trla_count)%x >=                      &
                                        (nx + 1)* dx - 1.0E-12_wp )            &
                               THEN
                                  trla(trla_count)%x = trla(trla_count)%x -    &
                                                   1.0E-10_wp
                                  trla(trla_count)%origin_x =                  &
                                                   trla(trla_count)%origin_x - 1
                               ENDIF

                            ENDIF

                         ELSEIF ( ibc_mas_lr == 1 )  THEN
!
!--                         Agent absorption
                            agents(n)%agent_mask = .FALSE.
                            deleted_agents = deleted_agents + 1

                         ENDIF
                      ELSE
!
!--                      Store agent data in the transfer array, which will be 
!--                      send to the neighbouring PE
                         trla_count = trla_count + 1
                         trla(trla_count) = agents(n)
                         agents(n)%agent_mask = .FALSE.
                         deleted_agents = deleted_agents + 1

                      ENDIF

                   ELSEIF ( i > nxr )  THEN
                      IF ( i > nx )  THEN
!
!--                      Apply boundary condition along x
                         IF ( ibc_mas_lr == 0 )  THEN
!
!--                         Cyclic condition
                            IF ( pdims(1) == 1 )  THEN
                               agents(n)%x = agents(n)%x - ( nx + 1 ) * dx
                               agents(n)%origin_x = agents(n)%origin_x - &
                               ( nx + 1 ) * dx
                            ELSE
                               trra_count = trra_count + 1
                               trra(trra_count) = agents(n)
                               trra(trra_count)%x = trra(trra_count)%x -       &
                                                    ( nx + 1 ) * dx
                               trra(trra_count)%origin_x =                     &
                                                   trra(trra_count)%origin_x - &
                                                   ( nx + 1 ) * dx
                               agents(n)%agent_mask = .FALSE.
                               deleted_agents = deleted_agents + 1

                            ENDIF

                         ELSEIF ( ibc_mas_lr == 1 )  THEN
!
!--                         Agent absorption
                            agents(n)%agent_mask = .FALSE.
                            deleted_agents = deleted_agents + 1

                         ENDIF
                      ELSE
!
!--                      Store agent data in the transfer array, which will be send
!--                      to the neighbouring PE
                         trra_count = trra_count + 1
                         trra(trra_count) = agents(n)
                         agents(n)%agent_mask = .FALSE.
                         deleted_agents = deleted_agents + 1

                      ENDIF

                   ENDIF
                ENDIF

             ENDDO
          ENDDO
       ENDDO

!
!--    Allocate arrays required for north-south exchange, as these
!--    are used directly after agents are exchange along x-direction.
       ALLOCATE( move_also_north(1:NR_2_direction_move) )
       ALLOCATE( move_also_south(1:NR_2_direction_move) )

       nr_move_north = 0
       nr_move_south = 0
!
!--    Send left boundary, receive right boundary (but first exchange how many
!--    and chec if agent storage must be extended)
       IF ( pdims(1) /= 1 )  THEN

          CALL MPI_SENDRECV( trla_count,      1, MPI_INTEGER, pleft,  0, &
                             trra_count_recv, 1, MPI_INTEGER, pright, 0, &
                             comm2d, status, ierr )

          ALLOCATE(rvra(MAX(1,trra_count_recv)))
!
!--       This MPI_SENDRECV should work even with odd mixture on 32 and 64 Bit 
!--       variables in structure agent_type (due to the calculation of par_size)
          par_size = STORAGE_SIZE(trla(1))/8
          CALL MPI_SENDRECV( trla, max(1,trla_count)*par_size, MPI_BYTE, pleft,&
                    1, rvra, max(1,trra_count_recv)*par_size, MPI_BYTE, pright,&
                             1, comm2d, status, ierr )

          IF ( trra_count_recv > 0 ) THEN
             CALL mas_eh_add_agents_to_gridcell(rvra(1:trra_count_recv))
          ENDIF

          DEALLOCATE(rvra)

!
!--       Send right boundary, receive left boundary
          CALL MPI_SENDRECV( trra_count,      1, MPI_INTEGER, pright, 0,       &
                             trla_count_recv, 1, MPI_INTEGER, pleft,  0,       &
                             comm2d, status, ierr )

          ALLOCATE(rvla(MAX(1,trla_count_recv)))
!
!--       This MPI_SENDRECV should work even with odd mixture on 32 and 64 Bit 
!--       variables in structure agent_type (due to the calculation of par_size)
          par_size = STORAGE_SIZE(trra(1))/8
          CALL MPI_SENDRECV( trra, max(1,trra_count)*par_size, MPI_BYTE,       &
                             pright, 1, rvla,                                  &
                             max(1,trla_count_recv)*par_size, MPI_BYTE,        &
                             pleft, 1, comm2d, status, ierr )

          IF ( trla_count_recv > 0 ) THEN
             CALL mas_eh_add_agents_to_gridcell(rvla(1:trla_count_recv))
          ENDIF

          DEALLOCATE( rvla )
          DEALLOCATE( trla, trra )

       ENDIF

!
!--    Check whether agents have crossed the boundaries in y direction. Note 
!--    that this case can also apply to agents that have just been received
!--    from the adjacent right or left PE.
!--    Find out first the number of agents to be transferred and allocate
!--    temporary arrays needed to store them.
!--    For a one-dimensional decomposition along y, no transfer is necessary,
!--    because the agent remains on the PE.
       trsa_count  = nr_move_south
       trna_count  = nr_move_north

       trsa_count_recv   = 0
       trna_count_recv   = 0

       IF ( pdims(2) /= 1 )  THEN
!
!--       First calculate the storage necessary for sending and receiving the
!--       data
          DO  ip = nxl, nxr
             DO  jp = nys, nyn, nyn-nys    !compute only first (nys) and last (nyn) loop iterration
                number_of_agents = agt_count(jp,ip)
                IF ( number_of_agents <= 0 )  CYCLE
                agents => grid_agents(jp,ip)%agents(1:number_of_agents)
                DO  n = 1, number_of_agents
                   IF ( agents(n)%agent_mask )  THEN
                      j = agents(n)%y * ddy
!
!--                   Above calculation does not work for indices less than zero
                      IF ( agents(n)%y < 0.0_wp )  j = -1

                      IF ( j < nys )  THEN
                         trsa_count = trsa_count + 1
                      ELSEIF ( j > nyn )  THEN
                         trna_count = trna_count + 1
                      ENDIF
                   ENDIF
                ENDDO
             ENDDO
          ENDDO

          IF ( trsa_count  == 0 )  trsa_count  = 1
          IF ( trna_count  == 0 )  trna_count  = 1

          ALLOCATE( trsa(trsa_count), trna(trna_count) )

          trsa = zero_agent
          trna = zero_agent

          trsa_count  = nr_move_south
          trna_count  = nr_move_north

          trsa(1:nr_move_south) = move_also_south(1:nr_move_south)
          trna(1:nr_move_north) = move_also_north(1:nr_move_north)

       ENDIF

       DO  ip = nxl, nxr
          DO  jp = nys, nyn, nyn-nys ! compute only first (nys) and last (nyn) loop iterration
             number_of_agents = agt_count(jp,ip)
             IF ( number_of_agents <= 0 )  CYCLE
             agents => grid_agents(jp,ip)%agents(1:number_of_agents)
             DO  n = 1, number_of_agents
!
!--             Only those agents that have not been marked as 'deleted' may
!--             be moved.
                IF ( agents(n)%agent_mask )  THEN

                   j = agents(n)%y * ddy
!
!--                Above calculation does not work for indices less than zero
                   IF ( agents(n)%y < 0.0_wp * dy )  j = -1

                   IF ( j < nys )  THEN
                      IF ( j < 0 )  THEN
!
!--                      Apply boundary condition along y
                         IF ( ibc_mas_ns == 0 )  THEN
!
!--                         Cyclic condition
                            IF ( pdims(2) == 1 )  THEN
                               agents(n)%y = ( ny + 1 ) * dy + agents(n)%y
                               agents(n)%origin_y = ( ny + 1 ) * dy +          &
                                                     agents(n)%origin_y
                            ELSE
                               trsa_count         = trsa_count + 1
                               trsa(trsa_count)   = agents(n)
                               trsa(trsa_count)%y = ( ny + 1 ) * dy +          &
                                                 trsa(trsa_count)%y
                               trsa(trsa_count)%origin_y =                     &
                                                  trsa(trsa_count)%origin_y    &
                                                + ( ny + 1 ) * dy
                               agents(n)%agent_mask = .FALSE.
                               deleted_agents = deleted_agents + 1

                               IF ( trsa(trsa_count)%y >=                      &
                                                    (ny+1)* dy - 1.0E-12_wp )  &
                               THEN
                                  trsa(trsa_count)%y = trsa(trsa_count)%y -    &
                                                       1.0E-10_wp
                                  trsa(trsa_count)%origin_y =                  &
                                                  trsa(trsa_count)%origin_y - 1
                               ENDIF

                            ENDIF

                         ELSEIF ( ibc_mas_ns == 1 )  THEN
!
!--                         Agent absorption
                            agents(n)%agent_mask = .FALSE.
                            deleted_agents          = deleted_agents + 1

                         ENDIF
                      ELSE
!
!--                      Store agent data in the transfer array, which will 
!--                      be send to the neighbouring PE
                         trsa_count = trsa_count + 1
                         trsa(trsa_count) = agents(n)
                         agents(n)%agent_mask = .FALSE.
                         deleted_agents = deleted_agents + 1

                      ENDIF

                   ELSEIF ( j > nyn )  THEN
                      IF ( j > ny )  THEN
!
!--                      Apply boundary condition along y
                         IF ( ibc_mas_ns == 0 )  THEN
!
!--                         Cyclic condition
                            IF ( pdims(2) == 1 )  THEN
                               agents(n)%y        = agents(n)%y -              &
                                                    ( ny + 1 ) * dy
                               agents(n)%origin_y = agents(n)%origin_y -       &
                                                    ( ny + 1 ) * dy
                            ELSE
                               trna_count         = trna_count + 1
                               trna(trna_count)   = agents(n)
                               trna(trna_count)%y =                            &
                                          trna(trna_count)%y - ( ny + 1 ) * dy
                               trna(trna_count)%origin_y =                     &
                                         trna(trna_count)%origin_y -           &
                                         ( ny + 1 ) * dy
                               agents(n)%agent_mask = .FALSE.
                               deleted_agents          = deleted_agents + 1
                            ENDIF

                         ELSEIF ( ibc_mas_ns == 1 )  THEN
!
!--                         Agent absorption
                            agents(n)%agent_mask = .FALSE.
                            deleted_agents = deleted_agents + 1

                         ENDIF
                      ELSE
!
!--                      Store agent data in the transfer array, which will 
!--                      be send to the neighbouring PE
                         trna_count = trna_count + 1
                         trna(trna_count) = agents(n)
                         agents(n)%agent_mask = .FALSE.
                         deleted_agents = deleted_agents + 1

                      ENDIF

                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO

!
!--    Send front boundary, receive back boundary (but first exchange how many
!--    and chec if agent storage must be extended)
       IF ( pdims(2) /= 1 )  THEN

          CALL MPI_SENDRECV( trsa_count,      1, MPI_INTEGER, psouth, 0,       &
                             trna_count_recv, 1, MPI_INTEGER, pnorth, 0,       &
                             comm2d, status, ierr )

          ALLOCATE(rvna(MAX(1,trna_count_recv)))
!
!--       This MPI_SENDRECV should work even with odd mixture on 32 and 64 Bit 
!--       variables in structure agent_type (due to the calculation of par_size)
          par_size = STORAGE_SIZE(trsa(1))/8
          CALL MPI_SENDRECV( trsa, trsa_count*par_size, MPI_BYTE,              &
                             psouth, 1, rvna,                                  &
                             trna_count_recv*par_size, MPI_BYTE, pnorth, 1,    &
                             comm2d, status, ierr )

          IF ( trna_count_recv  > 0 ) THEN
             CALL mas_eh_add_agents_to_gridcell(rvna(1:trna_count_recv))
          ENDIF

          DEALLOCATE(rvna)

!
!--       Send back boundary, receive front boundary
          CALL MPI_SENDRECV( trna_count,      1, MPI_INTEGER, pnorth, 0,       &
                             trsa_count_recv, 1, MPI_INTEGER, psouth, 0,       &
                             comm2d, status, ierr )

          ALLOCATE(rvsa(MAX(1,trsa_count_recv)))
!
!--       This MPI_SENDRECV should work even with odd mixture on 32 and 64 Bit 
!--       variables in structure agent_type (due to the calculation of par_size)
          par_size = STORAGE_SIZE(trna(1))/8
          CALL MPI_SENDRECV( trna, trna_count*par_size, MPI_BYTE,              &
                             pnorth, 1, rvsa,                                  &
                             trsa_count_recv*par_size, MPI_BYTE, psouth, 1,    &
                             comm2d, status, ierr )

          IF ( trsa_count_recv > 0 ) THEN
             CALL mas_eh_add_agents_to_gridcell(rvsa(1:trsa_count_recv))
          ENDIF

          DEALLOCATE(rvsa)

          number_of_agents = number_of_agents + trsa_count_recv

          DEALLOCATE( trsa, trna )

       ENDIF

       DEALLOCATE( move_also_north )
       DEALLOCATE( move_also_south )

!
!--    Accumulate the number of agents transferred between the subdomains)
       CALL mas_eh_ghost_exchange

#else

       DO  ip = nxl, nxr, nxr-nxl
          DO  jp = nys, nyn
             number_of_agents = agt_count(jp,ip)
             IF ( number_of_agents <= 0 )  CYCLE
             agents => grid_agents(jp,ip)%agents(1:number_of_agents)
             DO  n = 1, number_of_agents
!
!--             Apply boundary conditions
                IF ( agents(n)%x < 0.0_wp )  THEN

                   IF ( ibc_mas_lr == 0 )  THEN
!
!--                   Cyclic boundary. Relevant coordinate has to be changed.
                      agents(n)%x = ( nx + 1 ) * dx + agents(n)%x
                      agents(n)%origin_x = ( nx + 1 ) * dx + &
                                  agents(n)%origin_x
                   ELSEIF ( ibc_mas_lr == 1 )  THEN
!
!--                   Agent absorption
                      agents(n)%agent_mask = .FALSE.
                      deleted_agents = deleted_agents + 1
                   ENDIF

                ELSEIF ( agents(n)%x >= ( nx + 1 ) * dx )  THEN

                   IF ( ibc_mas_lr == 0 )  THEN
!
!--                   Cyclic boundary. Relevant coordinate has to be changed.
                      agents(n)%x = agents(n)%x - ( nx + 1 ) * dx

                   ELSEIF ( ibc_mas_lr == 1 )  THEN
!
!--                   Agent absorption
                      agents(n)%agent_mask = .FALSE.
                      deleted_agents = deleted_agents + 1
                   ENDIF

                ENDIF
             ENDDO
          ENDDO
       ENDDO

       DO  ip = nxl, nxr
          DO  jp = nys, nyn, nyn-nys
             number_of_agents = agt_count(jp,ip)
             IF ( number_of_agents <= 0 )  CYCLE
             agents => grid_agents(jp,ip)%agents(1:number_of_agents)
             DO  n = 1, number_of_agents

                IF ( agents(n)%y < 0.0_wp )  THEN

                   IF ( ibc_mas_ns == 0 )  THEN
!
!--                   Cyclic boundary. Relevant coordinate has to be changed.
                      agents(n)%y = ( ny + 1 ) * dy + agents(n)%y
                      agents(n)%origin_y = ( ny + 1 ) * dy + &
                           agents(n)%origin_y

                   ELSEIF ( ibc_mas_ns == 1 )  THEN
!
!--                   Agent absorption
                      agents(n)%agent_mask = .FALSE.
                      deleted_agents = deleted_agents + 1
                   ENDIF

                ELSEIF ( agents(n)%y >= ( ny + 0.5_wp ) * dy )  THEN

                   IF ( ibc_mas_ns == 0 )  THEN
!
!--                   Cyclic boundary. Relevant coordinate has to be changed.
                      agents(n)%y = agents(n)%y - ( ny + 1 ) * dy

                   ELSEIF ( ibc_mas_ns == 1 )  THEN
!
!--                   Agent absorption
                      agents(n)%agent_mask = .FALSE.
                      deleted_agents = deleted_agents + 1
                   ENDIF

                ENDIF

             ENDDO
          ENDDO
       ENDDO
#endif

    END SUBROUTINE mas_eh_exchange_horiz


#if defined( __parallel )
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Sends the agents from the three gridcells closest to the 
!> north/south/left/right border of a PE to the corresponding neighbors ghost
!> layer (which is three grid boxes deep)
!------------------------------------------------------------------------------!
    SUBROUTINE mas_eh_ghost_exchange

       IMPLICIT NONE

       INTEGER(iwp) ::  ip          !< index variable along x
       INTEGER(iwp) ::  jp          !< index variable along y
       INTEGER(iwp) ::  agt_size    !< Bit size of agent datatype
       INTEGER(iwp) ::  ghla_count  !< ghost points left agent
       INTEGER(iwp) ::  ghna_count  !< ghost points north agent
       INTEGER(iwp) ::  ghra_count  !< ghost points right agent
       INTEGER(iwp) ::  ghsa_count  !< ghost points south agent

       LOGICAL ::  ghla_empty      !< ghost points left agent
       LOGICAL ::  ghla_empty_rcv  !< ghost points left agent
       LOGICAL ::  ghna_empty      !< ghost points north agent
       LOGICAL ::  ghna_empty_rcv  !< ghost points north agent
       LOGICAL ::  ghra_empty      !< ghost points right agent
       LOGICAL ::  ghra_empty_rcv  !< ghost points right agent
       LOGICAL ::  ghsa_empty      !< ghost points south agent
       LOGICAL ::  ghsa_empty_rcv  !< ghost points south agent

       TYPE(agent_type), DIMENSION(:), ALLOCATABLE ::  ghla  !< agents received from right PE
       TYPE(agent_type), DIMENSION(:), ALLOCATABLE ::  ghna  !< agents received from south PE
       TYPE(agent_type), DIMENSION(:), ALLOCATABLE ::  ghra  !< agents received from left PE
       TYPE(agent_type), DIMENSION(:), ALLOCATABLE ::  ghsa  !< agents received from north PE

       ghla_empty = .TRUE.
       ghna_empty = .TRUE.
       ghra_empty = .TRUE.
       ghsa_empty = .TRUE.
!
!--    reset ghost layer
       DO ip = nxlg, nxl-1
          DO jp = nysg, nyng
             agt_count(jp,ip) = 0
          ENDDO
       ENDDO
       DO ip = nxr+1, nxrg
          DO jp = nysg, nyng
             agt_count(jp,ip) = 0
          ENDDO
       ENDDO
       DO ip = nxl, nxr
          DO jp = nysg, nys-1
             agt_count(jp,ip) = 0
          ENDDO
       ENDDO
       DO ip = nxl, nxr
          DO jp = nyn+1, nyng
             agt_count(jp,ip) = 0
          ENDDO
       ENDDO
!
!--    Transfer of agents from left to right and vice versa
       IF ( pdims(1) /= 1 )  THEN
!
!--       Reset left and right ghost layers
          ghla_count  = 0
          ghra_count  = 0
!
!--       First calculate the storage necessary for sending 
!--       and receiving the data.
          ghla_count = SUM(agt_count(nys:nyn,nxl:nxl+2))
          ghra_count = SUM(agt_count(nys:nyn,nxr-2:nxr))
!
!--       No cyclic boundaries for agents
          IF ( nxl == 0 .OR. ghla_count == 0 ) THEN
             ghla_count = 1
          ELSE
             ghla_empty = .FALSE.
          ENDIF
          IF ( nxr == nx .OR. ghra_count == 0 ) THEN
             ghra_count = 1
          ELSE
             ghra_empty = .FALSE.
          ENDIF
          ALLOCATE( ghla(1:ghla_count), ghra(1:ghra_count) )
          ghla = zero_agent
          ghra = zero_agent
!
!--       Get all agents that will be sent left into one array
          ghla_count = 0
          IF ( nxl /= 0 ) THEN
             DO ip = nxl, nxl+2
                DO jp = nys, nyn

                   number_of_agents = agt_count(jp,ip)
                   IF ( number_of_agents <= 0 )  CYCLE
                   ghla(ghla_count+1:ghla_count+number_of_agents)              &
                       = grid_agents(jp,ip)%agents(1:number_of_agents)
                   ghla_count = ghla_count + number_of_agents

                ENDDO
             ENDDO
          ENDIF
          IF ( ghla_count == 0 )  ghla_count = 1
!
!--       Get all agents that will be sent right into one array
          ghra_count = 0
          IF ( nxr /= nx ) THEN
             DO ip = nxr-2, nxr
                DO jp = nys, nyn

                   number_of_agents = agt_count(jp,ip)
                   IF ( number_of_agents <= 0 )  CYCLE
                   ghra(ghra_count+1:ghra_count+number_of_agents)              &
                       = grid_agents(jp,ip)%agents(1:number_of_agents)
                   ghra_count = ghra_count + number_of_agents

                ENDDO
             ENDDO
          ENDIF
          IF ( ghra_count == 0 ) ghra_count = 1
!
!--       Send/receive number of agents that 
!--       will be transferred to/from left/right neighbor
          CALL MPI_SENDRECV( ghla_count,      1, MPI_INTEGER, pleft,  0,       &
                             ghra_count_recv, 1, MPI_INTEGER, pright, 0,       &
                             comm2d, status, ierr )
          ALLOCATE ( agt_gh_r(1:ghra_count_recv) )
!
!--       Send/receive number of agents that 
!--       will be transferred to/from right/left neighbor
          CALL MPI_SENDRECV( ghra_count,      1, MPI_INTEGER, pright,  0,      &
                             ghla_count_recv, 1, MPI_INTEGER, pleft,   0,      &
                             comm2d, status, ierr )
!
!--       Send/receive flag that indicates if there are actually any agents 
!--       in ghost  layer
          CALL MPI_SENDRECV( ghla_empty,     1, MPI_LOGICAL, pleft, 1,         &
                             ghra_empty_rcv, 1, MPI_LOGICAL, pright,1,         &
                             comm2d, status, ierr )
          CALL MPI_SENDRECV( ghra_empty,     1, MPI_LOGICAL, pright,1,         &
                             ghla_empty_rcv, 1, MPI_LOGICAL, pleft, 1,         &
                             comm2d, status, ierr )


          ALLOCATE ( agt_gh_l(1:ghla_count_recv) )
!
!--       Get bit size of one agent
          agt_size = STORAGE_SIZE(zero_agent)/8
!
!--       Send/receive agents to/from left/right neighbor
          CALL MPI_SENDRECV( ghla,     ghla_count      * agt_size, MPI_BYTE,   &
                                pleft, 1,                                      &
                             agt_gh_r, ghra_count_recv * agt_size, MPI_BYTE,   &
                                pright,1,                                      &
                             comm2d, status, ierr )
!
!--       Send/receive agents to/from left/right neighbor
          CALL MPI_SENDRECV( ghra,     ghra_count      * agt_size, MPI_BYTE,   &
                                pright,1,                                      &
                             agt_gh_l, ghla_count_recv * agt_size, MPI_BYTE,   &
                                pleft, 1,                                      &
                             comm2d, status, ierr )
!
!--       If agents were received, add them to the respective ghost layer cells
          IF ( .NOT. ghra_empty_rcv ) THEN
             CALL mas_eh_add_ghost_agents_to_gridcell(agt_gh_r)
          ENDIF

          IF ( .NOT. ghla_empty_rcv ) THEN
             CALL mas_eh_add_ghost_agents_to_gridcell(agt_gh_l)
          ENDIF

          DEALLOCATE( ghla, ghra, agt_gh_l, agt_gh_r )

       ENDIF

!
!--    Transfer of agents from south to north and vice versa
       IF ( pdims(2) /= 1 )  THEN
!
!--       Reset south and north ghost layers
          ghsa_count  = 0
          ghna_count  = 0
!
!--       First calculate the storage necessary for sending 
!--       and receiving the data.
          ghsa_count = SUM(agt_count(nys:nys+2,nxlg:nxrg))
          ghna_count = SUM(agt_count(nyn-2:nyn,nxlg:nxrg))
!
!--       No cyclic boundaries for agents
          IF ( nys == 0 .OR. ghsa_count == 0 ) THEN
             ghsa_count = 1
          ELSE
             ghsa_empty = .FALSE.
          ENDIF
          IF ( nyn == ny .OR. ghna_count == 0 ) THEN
             ghna_count = 1
          ELSE
             ghna_empty = .FALSE.
          ENDIF
          ALLOCATE( ghsa(1:ghsa_count), ghna(1:ghna_count) )
          ghsa = zero_agent
          ghna = zero_agent
!
!--       Get all agents that will be sent south into one array
          ghsa_count = 0
          IF ( nys /= 0 ) THEN
             DO ip = nxlg, nxrg
                DO jp = nys, nys+2

                   number_of_agents = agt_count(jp,ip)
                   IF ( number_of_agents <= 0 )  CYCLE
                   ghsa(ghsa_count+1:ghsa_count+number_of_agents)              &
                       = grid_agents(jp,ip)%agents(1:number_of_agents)
                   ghsa_count = ghsa_count + number_of_agents

                ENDDO
             ENDDO
          ENDIF
          IF ( ghsa_count == 0 )  ghsa_count = 1
!
!--       Get all agents that will be sent north into one array
          ghna_count = 0
          IF ( nyn /= ny ) THEN
             DO ip = nxlg, nxrg
                DO jp = nyn-2, nyn

                   number_of_agents = agt_count(jp,ip)
                   IF ( number_of_agents <= 0 )  CYCLE
                   ghna(ghna_count+1:ghna_count+number_of_agents)              &
                       = grid_agents(jp,ip)%agents(1:number_of_agents)
                   ghna_count = ghna_count + number_of_agents

                ENDDO
             ENDDO
          ENDIF
          IF ( ghna_count == 0 ) ghna_count = 1
!
!--       Send/receive number of agents that 
!--       will be transferred to/from south/north neighbor
          CALL MPI_SENDRECV( ghsa_count, 1, MPI_INTEGER, psouth, 0,            &
                             ghna_count_recv,   1, MPI_INTEGER, pnorth, 0,     &
                             comm2d, status, ierr )
          ALLOCATE ( agt_gh_n(1:ghna_count_recv) )
!
!--       Send/receive number of agents that 
!--       will be transferred to/from north/south neighbor
          CALL MPI_SENDRECV( ghna_count, 1, MPI_INTEGER, pnorth, 0,            &
                             ghsa_count_recv,   1, MPI_INTEGER, psouth, 0,     &
                             comm2d, status, ierr )
!
!--       Send/receive flag that indicates if there are actually any agents 
!--       in ghost  layer
          CALL MPI_SENDRECV( ghsa_empty,     1, MPI_LOGICAL, psouth, 1,        &
                             ghna_empty_rcv, 1, MPI_LOGICAL, pnorth, 1,        &
                             comm2d, status, ierr )
          CALL MPI_SENDRECV( ghna_empty,     1, MPI_LOGICAL, pnorth, 1,        &
                             ghsa_empty_rcv, 1, MPI_LOGICAL, psouth, 1,        &
                             comm2d, status, ierr )


          ALLOCATE ( agt_gh_s(1:ghsa_count_recv) )
!
!--       Get bit size of one agent
          agt_size = STORAGE_SIZE(zero_agent)/8
!
!--       Send/receive agents to/from south/north neighbor
          CALL MPI_SENDRECV( ghsa,     ghsa_count      * agt_size, MPI_BYTE,   &
                                psouth,1,                                      &
                             agt_gh_n, ghna_count_recv * agt_size, MPI_BYTE,   &
                                pnorth,1,                                      &
                             comm2d, status, ierr )
!
!--       Send/receive agents to/from south/north neighbor
          CALL MPI_SENDRECV( ghna,     ghna_count      * agt_size, MPI_BYTE,   &
                                pnorth,1,                                      &
                             agt_gh_s, ghsa_count_recv * agt_size, MPI_BYTE,   &
                                psouth,1,                                      &
                             comm2d, status, ierr )
!
!--       If agents were received, add them to the respective ghost layer cells
          IF ( .NOT. ghna_empty_rcv ) THEN
             CALL mas_eh_add_ghost_agents_to_gridcell(agt_gh_n)
          ENDIF

          IF ( .NOT. ghsa_empty_rcv ) THEN
             CALL mas_eh_add_ghost_agents_to_gridcell(agt_gh_s)
          ENDIF

          DEALLOCATE( ghna, ghsa, agt_gh_n, agt_gh_s )

       ENDIF

    END SUBROUTINE mas_eh_ghost_exchange
#endif

!------------------------------------------------------------------------------!
! Description:
! ------------
!> If an agent moves from one grid cell to another (on the current 
!> processor!), this subroutine moves the corresponding element from the
!> agent array of the old grid cell to the agent array of the new grid 
!> cell.
!------------------------------------------------------------------------------!
    SUBROUTINE mas_eh_move_agent

       IMPLICIT NONE

       INTEGER(iwp) ::  i              !< grid index (x) of agent position
       INTEGER(iwp) ::  ip             !< index variable along x
       INTEGER(iwp) ::  j              !< grid index (y) of agent position
       INTEGER(iwp) ::  jp             !< index variable along y
       INTEGER(iwp) ::  n              !< index variable for agent array
       INTEGER(iwp) ::  na_before_move !< number of agents per grid box before moving
       INTEGER(iwp) ::  aindex         !< dummy argument for number of new agent per grid box

       TYPE(agent_type), DIMENSION(:), POINTER ::  agents_before_move !< agents before moving

       DO  ip = nxl, nxr
          DO  jp = nys, nyn

                na_before_move = agt_count(jp,ip)
                IF ( na_before_move <= 0 )  CYCLE
                agents_before_move => grid_agents(jp,ip)%agents(1:na_before_move)

                DO  n = 1, na_before_move
                   i = agents_before_move(n)%x * ddx
                   j = agents_before_move(n)%y * ddy

!--                For mas_eh_exchange_horiz to work properly agents need to be
!--                moved to the outermost gridboxes of the respective processor.
!--                If the agent index is inside the processor the following
!--                lines will not change the index
                   i = MIN ( i , nxr )
                   i = MAX ( i , nxl )
                   j = MIN ( j , nyn )
                   j = MAX ( j , nys )

!
!--                Check if agent has moved to another grid cell.
                   IF ( i /= ip  .OR.  j /= jp )  THEN
!
!--                   If the agent stays on the same processor, the agent
!--                   will be added to the agent array of the new processor.
                      number_of_agents = agt_count(j,i)
                      agents => grid_agents(j,i)%agents(1:number_of_agents)

                      aindex = number_of_agents+1
                      IF (  aindex > SIZE(grid_agents(j,i)%agents)  )     &
                      THEN
                         CALL mas_eh_realloc_agents_array(i,j)
                      ENDIF

                      grid_agents(j,i)%agents(aindex) = agents_before_move(n)
                      agt_count(j,i) = aindex

                      agents_before_move(n)%agent_mask = .FALSE.
                   ENDIF
                ENDDO

          ENDDO
       ENDDO

       RETURN

    END SUBROUTINE mas_eh_move_agent

!------------------------------------------------------------------------------!
! Description:
! ------------
!> If the allocated memory for the agent array do not suffice to add arriving
!> agents from neighbour grid cells, this subrouting reallocates the 
!> agent array to assure enough memory is available. 
!------------------------------------------------------------------------------!
    SUBROUTINE mas_eh_realloc_agents_array (i,j,size_in)

       IMPLICIT NONE

       INTEGER(iwp) :: old_size  !< old array size
       INTEGER(iwp) :: new_size  !< new array size

       INTEGER(iwp), INTENT(in) ::  i  !< grid index (y)
       INTEGER(iwp), INTENT(in) ::  j  !< grid index (y)

       INTEGER(iwp), INTENT(in), OPTIONAL ::  size_in  !< size of input array

       TYPE(agent_type), DIMENSION(10) :: tmp_agents_s !< temporary static agent array

       TYPE(agent_type), DIMENSION(:), ALLOCATABLE :: tmp_agents_d !< temporary dynamic agent array

       old_size = SIZE(grid_agents(j,i)%agents)

       IF ( PRESENT(size_in) )   THEN
          new_size = size_in
       ELSE
          new_size = old_size * ( 1.0_wp + alloc_factor_mas / 100.0_wp )
       ENDIF

       new_size = MAX( new_size, min_nr_agent, old_size + 1 )

       IF ( old_size <= 10 )  THEN

          tmp_agents_s(1:old_size) = grid_agents(j,i)%agents(1:old_size)

          DEALLOCATE(grid_agents(j,i)%agents)
          ALLOCATE(grid_agents(j,i)%agents(new_size))

          grid_agents(j,i)%agents(1:old_size)         = tmp_agents_s(1:old_size)
          grid_agents(j,i)%agents(old_size+1:new_size) = zero_agent

       ELSE

          ALLOCATE(tmp_agents_d(new_size))
          tmp_agents_d(1:old_size) = grid_agents(j,i)%agents

          DEALLOCATE(grid_agents(j,i)%agents)
          ALLOCATE(grid_agents(j,i)%agents(new_size))

          grid_agents(j,i)%agents(1:old_size)         = tmp_agents_d(1:old_size)
          grid_agents(j,i)%agents(old_size+1:new_size) = zero_agent

          DEALLOCATE(tmp_agents_d)

       ENDIF
       agents => grid_agents(j,i)%agents(1:number_of_agents)

       RETURN
    END SUBROUTINE mas_eh_realloc_agents_array

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Inquires prognostic model quantities at the position of each agent and
!> stores them in that agent for later output
!------------------------------------------------------------------------------! 
    SUBROUTINE mas_get_prognostic_quantities

       USE arrays_3d,                                                          &
           ONLY:  u, v, pt, exner

       IMPLICIT NONE

       INTEGER(iwp) ::  i_offset  !<  index offset for windspeed measurement
       INTEGER(iwp) ::  il        !<  x-index
       INTEGER(iwp) ::  is        !<  subgrid box counter
       INTEGER(iwp) ::  j_offset  !<  index offset for windspeed measurement
       INTEGER(iwp) ::  jl        !<  y-index
       INTEGER(iwp) ::  kl        !<  z-index
       INTEGER(iwp) ::  nl        !<  agent counter
       INTEGER(iwp) ::  se        !<  subgrid box end index
       INTEGER(iwp) ::  si        !<  subgrid box start index

       REAL(wp) ::  u_a  !< windspeed at agent position (x)
       REAL(wp) ::  v_a  !< windspeed at agent position (y)

       DO  il = nxl, nxr
          DO  jl = nys, nyn

             number_of_agents = agt_count(jl,il)
!
!--          If grid cell is empty, cycle
             IF ( number_of_agents <= 0 ) CYCLE
             kl = s_measure_height(jl,il)

             agents => grid_agents(jl,il)%agents(1:number_of_agents)
!
!--          loop over the four subgrid boxes
             DO is = 0,3
!
!--             Set indices
                si = grid_agents(jl,il)%start_index(is)
                se = grid_agents(jl,il)%end_index(is)
                DO nl = si, se
!
!--                Calculate index offset in x-direction:
!--                Left value if wall right of grid box
!--                Right value if wall left of grid box
!--                Else the one that is closer to the agent
                   IF ( BTEST( obstacle_flags( jl, il+1 ), 6 ) ) THEN
                      i_offset = 0
                   ELSEIF ( BTEST( obstacle_flags( jl, il-1 ), 2 ) ) THEN
                      i_offset = 1
                   ELSE
                      i_offset = MERGE( 0, 1, BTEST(is,1) )
                   ENDIF
                   u_a = u( kl, jl, il + i_offset )
!
!--                Calculate index offset in y-direction:
!--                South value if wall north of grid box
!--                North value if wall south of grid box
!--                Else the one that is closer to the agent
                   IF ( BTEST( obstacle_flags( jl+1, il ), 4 ) ) THEN
                      j_offset = 0
                   ELSEIF ( BTEST( obstacle_flags( jl-1, il ), 0 ) ) THEN
                      j_offset = 1
                   ELSE
                      j_offset = MERGE( 0, 1, BTEST(is,0) )
                   ENDIF
                   v_a = v( kl, jl + j_offset, il )
!
!--                Calculate windspeed at agent postion
                   agents(nl)%windspeed = SQRT(u_a**2 + v_a**2)
!
!--                Calculate temperature at agent position
                   agents(nl)%t = pt(kl,jl,il) * exner(kl)

                ENDDO

             ENDDO

          ENDDO
       ENDDO

    END SUBROUTINE mas_get_prognostic_quantities

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Adds an item to the priority queue (binary heap) at the correct position
!------------------------------------------------------------------------------! 
    SUBROUTINE mas_heap_insert_item( id, priority )

       IMPLICIT NONE

       INTEGER(iwp) ::  cur_pos !< current position
       INTEGER(iwp) ::  id      !< mesh ID of item

       REAL(wp) ::  priority !< item priority

       TYPE(heap_item) ::  item !< heap item

       item%mesh_id  = id
       item%priority = priority
!
!--    Extend heap, if necessary
       IF ( heap_count + 1 > SIZE(queue) ) THEN
          CALL mas_heap_extend
       ENDIF
!
!--    Insert item at first unoccupied postion (highest index) of heap
       cur_pos = heap_count
       queue(cur_pos) = item
!
!--    Sort while inserted item is not at top of heap
       DO WHILE ( cur_pos /= 0 )
!
!--       If priority < its parent's priority, swap them. 
!--       Else, sorting is done.
          IF ( queue(cur_pos)%priority                                         &
              < queue(FLOOR((cur_pos)/2.))%priority )                          &
          THEN
             item = queue(cur_pos)
             queue(cur_pos) = queue(FLOOR((cur_pos)/2.))
             queue(FLOOR((cur_pos)/2.)) = item
             cur_pos = FLOOR((cur_pos)/2.)
          ELSE
             EXIT
          ENDIF
       ENDDO
!
!--    Item was added to heap, so the heap count increases
       heap_count = heap_count + 1

    END SUBROUTINE mas_heap_insert_item

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Extends the size of the priority queue (binary heap)
!------------------------------------------------------------------------------! 
    SUBROUTINE mas_heap_extend

       IMPLICIT NONE

       INTEGER(iwp) ::  soh !< size of heap

       TYPE(heap_item), DIMENSION(:), ALLOCATABLE ::  dummy_heap !< dummy heap

       soh = SIZE(queue)-1
       ALLOCATE(dummy_heap(0:soh))
       dummy_heap = queue
       DEALLOCATE(queue)
       ALLOCATE(queue(0:2*soh+1))
       queue(0:soh) = dummy_heap(0:soh)

    END SUBROUTINE mas_heap_extend

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Removes first (smallest) element from the priority queue, reorders the rest
!> and returns the ID of the removed mesh point
!------------------------------------------------------------------------------!
    SUBROUTINE mas_heap_extract_item ( id )

       IMPLICIT NONE

       INTEGER(iwp) ::  id      !< ID of item extracted item
       INTEGER(iwp) ::  child   !< child of item in heap
       INTEGER(iwp) ::  cur_pos !< current position of item in heap

       TYPE(heap_item) ::  dummy
!
!--    Get ID of mesh point with lowest priority (extracted item: top of heap)
       id = queue(0)%mesh_id
!
!--    Put last item in heap at first position
       queue(0) = queue(heap_count-1)
       cur_pos = 0
       DO
!
!--       If current item has no children, sorting is done
          IF( 2*cur_pos+1 > heap_count - 1 ) THEN
             EXIT
!
!--       If current item has only one child, check if item and its child are
!--       ordered correctly. Else, swap them.
          ELSEIF ( 2*cur_pos+2 > heap_count - 1 ) THEN
             IF ( queue(cur_pos)%priority > queue(2*cur_pos+1)%priority ) THEN
                dummy = queue(cur_pos)
                queue(cur_pos) = queue(2*cur_pos+1)
                queue(2*cur_pos+1) = dummy
                cur_pos = 2*cur_pos+1
             ELSE
                EXIT
             ENDIF
          ELSE
!
!--          determine the smaller child
             IF ( queue(2*cur_pos+1)%priority                                  &
                 >= queue(2*cur_pos+2)%priority )                              &
             THEN
                child = 2
             ELSE
                child = 1
             ENDIF
!
!--          Check if item and its smaller child are ordered falsely. If so, 
!--          swap them. Else, sorting is done.
             IF ( queue(cur_pos)%priority > queue(2*cur_pos+child )%priority ) &
             THEN
                dummy = queue(cur_pos)
                queue(cur_pos) = queue(2*cur_pos+child)
                queue(2*cur_pos+child) = dummy
                cur_pos = 2*cur_pos+child
             ELSE
                EXIT
             ENDIF
          ENDIF
       ENDDO
!
!--    Top item was removed from heap, thus, heap_cout decreases by one
       heap_count = heap_count-1

    END SUBROUTINE mas_heap_extract_item

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of Multi Agent System
!------------------------------------------------------------------------------!
    SUBROUTINE mas_init

       USE control_parameters,                                                 &
           ONLY:  coupling_char, initializing_actions, io_blocks, io_group

       USE arrays_3d,                                                          &
           ONLY:  zu, zw

       USE indices,                                                            &
           ONLY:  nzt 

       IMPLICIT NONE

       INTEGER(iwp) ::  i             !< grid cell (x)
       INTEGER(iwp) ::  ii            !< io-block counter
       INTEGER(iwp) ::  il            !< io-block counter
       INTEGER(iwp) ::  jl            !< io-block counter
       INTEGER(iwp) ::  kl            !< io-block counter
       INTEGER(iwp) ::  kdum            !< io-block counter
       INTEGER(iwp) ::  locdum            !< io-block counter
       INTEGER(iwp) ::  j             !< grid cell (y)
       INTEGER(iwp) ::  size_of_mesh  !< temporary value for read
       INTEGER(iwp) ::  size_of_pols  !< temporary value for read
       INTEGER(iwp) ::  ioerr         !< IOSTAT flag for IO-commands ( 0 = no error )

       LOGICAL ::  navigation_data_present  !< Flag: check for input file

       REAL(wp) ::  zdum  !< dummy for measurement height
       REAL(wp) ::  avg_agt_height = 1.8_wp


!
!--    Check the number of agent groups.
       IF ( number_of_agent_groups > max_number_of_agent_groups )  THEN
          WRITE( message_string, * ) 'max_number_of_agent_groups =',      &
                                     max_number_of_agent_groups ,         &
                                     '&number_of_agent_groups reset to ', &
                                     max_number_of_agent_groups
          CALL message( 'mas_init', 'PA0072', 0, 1, 0, 6, 0 )
          number_of_agent_groups = max_number_of_agent_groups
       ENDIF

!
!--    Set some parameters
       d_sigma_rep_agent = 1.0_wp/sigma_rep_agent
       d_sigma_rep_wall  = 1.0_wp/sigma_rep_wall
       d_tau_accel_agent = 1.0_wp/tau_accel_agent
       IF ( dt_agent /= 999.0_wp ) THEN
          agent_own_timestep = .TRUE.
       ENDIF

!
!--    Get index of first grid box above topography
       ALLOCATE( top_top_s(nysg:nyng,nxlg:nxrg),                               &
                 top_top_w(nysg:nyng,nxlg:nxrg),                               &
                 s_measure_height(nys:nyn,nxl:nxr) )
!
!--    Get first index above topography for scalar grid and last index in
!--    topography for z-component of wind
       DO il = nxlg, nxrg
          DO jl = nysg, nyng
             top_top_s(jl,il) = topo_top_ind(jl,il,0) + 1
             top_top_w(jl,il) = topo_top_ind(jl,il,3)
          ENDDO
       ENDDO
!
!--    Create 2D array containing the index at which measurements are done by
!--    agents. The height of this measurement is given by avg_agt_height.
       DO il = nxl, nxr
          DO jl = nys, nyn

             kdum = top_top_w(jl,il)
             zdum = zw(kdum)
             zdum = zdum + avg_agt_height
             locdum = 0
!
!--          Locate minimum distance from u-grid to measurement height (zdum)
             DO kl = 1, nzt
                IF ( ABS(zu(kl)-zdum) < ABS(zu(locdum)-zdum) ) locdum = kl
             ENDDO
             s_measure_height(jl,il) = locdum

          ENDDO
       ENDDO

       CALL mas_create_obstacle_flags

!
!--    Set default start positions, if necessary
       IF ( asl(1) == 9999999.9_wp )  asl(1) = 0.0_wp
       IF ( asr(1) == 9999999.9_wp )  asr(1) = ( nx + 1 ) * dx
       IF ( ass(1) == 9999999.9_wp )  ass(1) = 0.0_wp
       IF ( asn(1) == 9999999.9_wp )  asn(1) = ( ny + 1 ) * dy
       IF ( adx(1) == 9999999.9_wp .OR. adx(1) == 0.0_wp ) adx(1) = dx
       IF ( ady(1) == 9999999.9_wp .OR. ady(1) == 0.0_wp ) ady(1) = dy

       DO  j = 2, number_of_agent_groups
          IF ( asl(j) == 9999999.9_wp )  asl(j) = asl(j-1)
          IF ( asr(j) == 9999999.9_wp )  asr(j) = asr(j-1)
          IF ( ass(j) == 9999999.9_wp )  ass(j) = ass(j-1)
          IF ( asn(j) == 9999999.9_wp )  asn(j) = asn(j-1)
          IF ( adx(j) == 9999999.9_wp .OR. adx(j) == 0.0_wp ) adx(j) = adx(j-1)
          IF ( ady(j) == 9999999.9_wp .OR. ady(j) == 0.0_wp ) ady(j) = ady(j-1)
       ENDDO

!
!--    Check boundary condition and set internal variables
       SELECT CASE ( bc_mas_lr )

          CASE ( 'cyclic' )
             ibc_mas_lr = 0

          CASE ( 'absorb' )
             ibc_mas_lr = 1

          CASE DEFAULT
             WRITE( message_string, * ) 'unknown boundary condition ',         &
                                        'bc_mas_lr = "', TRIM( bc_mas_lr ), '"'
             CALL message( 'mas_init', 'PA0073', 1, 2, 0, 6, 0 )

       END SELECT
       SELECT CASE ( bc_mas_ns )

          CASE ( 'cyclic' )
             ibc_mas_ns = 0

          CASE ( 'absorb' )
             ibc_mas_ns = 1

          CASE DEFAULT
             WRITE( message_string, * ) 'unknown boundary condition ',         &
                                        'bc_mas_ns = "', TRIM( bc_mas_ns ), '"'
             CALL message( 'mas_init', 'PA0074', 1, 2, 0, 6, 0 )

       END SELECT

!
!--    For the first model run of a possible job chain initialize the
!--    agents, otherwise read the agent data from restart file.
       IF ( TRIM( initializing_actions ) == 'read_restart_data'                &
            .AND.  read_agents_from_restartfile )  THEN

!           CALL mas_read_restart_file

       ELSE
!
!--       Read preprocessed data of navigation mesh and building polygons
!--       for agent pathfinding
          DO ii = 0, io_blocks-1
             IF ( ii == io_group )  THEN
!
!--             Check for naviation input file and open it
                INQUIRE( FILE='NAVIGATION_DATA' // TRIM( coupling_char ), EXIST=navigation_data_present )
                IF ( .NOT. navigation_data_present )  THEN
                   message_string = 'Input file NAVIGATION_DATA' //                &
                                     TRIM( coupling_char ) // ' for MAS missing. ' // &
                                     '&Please run agent_preprocessing before the job to create it.'
                   CALL message( 'mas_init', 'PA0525', 1, 2, 0, 6, 0 )
                ENDIF
                OPEN ( 119, FILE='NAVIGATION_DATA'//TRIM( coupling_char ),     &
                            FORM='UNFORMATTED', IOSTAT=ioerr )
!
!--             Read mesh data
                READ(119) size_of_mesh
                ALLOCATE( mesh(1:size_of_mesh))
                DO i = 1, size_of_mesh
                   READ(119) mesh(i)%polygon_id, mesh(i)%vertex_id,            &
                             mesh(i)%noc, mesh(i)%origin_id,                   &
                             mesh(i)%cost_so_far, mesh(i)%x,                   &
                             mesh(i)%y, mesh(i)%x_s, mesh(i)%y_s
                   ALLOCATE( mesh(i)%connected_vertices(1:mesh(i)%noc),        &
                             mesh(i)%distance_to_vertex(1:mesh(i)%noc) )
                   DO j = 1, mesh(i)%noc
                      READ(119) mesh(i)%connected_vertices(j),                 &
                                mesh(i)%distance_to_vertex(j)
                   ENDDO
                ENDDO
!
!--             Read polygon data
                READ(119) size_of_pols
                ALLOCATE( polygons(1:size_of_pols) )
                DO i = 1, size_of_pols
                   READ(119) polygons(i)%nov
                   ALLOCATE( polygons(i)%vertices(0:polygons(i)%nov+1) )
                   DO j = 0, polygons(i)%nov+1
                      READ(119) polygons(i)%vertices(j)%delete,                &
                                polygons(i)%vertices(j)%x,                     &
                                polygons(i)%vertices(j)%y
                   ENDDO
                ENDDO
                CLOSE(119)

             ENDIF
#if defined( __parallel ) && ! defined ( __check )
             CALL MPI_BARRIER( comm2d, ierr )
#endif
          ENDDO

!
!--       Allocate agent arrays and set attributes of the initial set of
!--       agents, which can be also periodically released at later times.
          ALLOCATE( agt_count  (nysg:nyng,nxlg:nxrg),                          &
                    grid_agents(nysg:nyng,nxlg:nxrg) )
!
!--       Allocate dummy arrays for pathfinding
          ALLOCATE( dummy_path_x(0:agt_path_size),                 &
                    dummy_path_y(0:agt_path_size) )

          number_of_agents = 0
          sort_count_mas   = 0
          agt_count        = 0

!
!--       initialize counter for agent IDs 
          grid_agents%id_counter = 1

!
!--       Initialize all agents with dummy values (otherwise errors may
!--       occur within restart runs). The reason for this is still not clear
!--       and may be presumably caused by errors in the respective user-interface.
          zero_agent%agent_mask = .FALSE.
          zero_agent%block_nr      = -1
          zero_agent%group         = 0
          zero_agent%id            = 0_idp
          zero_agent%path_counter  = agt_path_size
          zero_agent%age           = 0.0_wp
          zero_agent%age_m         = 0.0_wp
          zero_agent%dt_sum        = 0.0_wp
          zero_agent%clo           = 0.0_wp
          zero_agent%energy_storage= 0.0_wp
          zero_agent%force_x       = 0.0_wp
          zero_agent%force_y       = 0.0_wp
          zero_agent%origin_x      = 0.0_wp
          zero_agent%origin_y      = 0.0_wp
          zero_agent%speed_abs     = 0.0_wp
          zero_agent%speed_e_x     = 0.0_wp
          zero_agent%speed_e_y     = 0.0_wp
          zero_agent%speed_des     = random_normal(desired_speed, des_sp_sig)
          zero_agent%speed_x       = 0.0_wp
          zero_agent%speed_y       = 0.0_wp
          zero_agent%ipt           = 0.0_wp
          zero_agent%x             = 0.0_wp
          zero_agent%y             = 0.0_wp
          zero_agent%path_x        = 0.0_wp
          zero_agent%path_y        = 0.0_wp
          zero_agent%t_x           = 0.0_wp
          zero_agent%t_y           = 0.0_wp

!
!--    Set a seed value for the random number generator to be exclusively
!--    used for the agent code. The generated random numbers should be
!--    different on the different PEs.
          iran_agent = iran_agent + myid

          CALL mas_create_agent (PHASE_INIT)

       ENDIF

!
!--    To avoid programm abort, assign agents array to the local version of 
!--    first grid cell
       number_of_agents = agt_count(nys,nxl)
       agents => grid_agents(nys,nxl)%agents(1:number_of_agents)

    END SUBROUTINE mas_init

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Output of informative message about maximum agent number
!------------------------------------------------------------------------------!
    SUBROUTINE mas_last_actions

       USE control_parameters,                                                 &
           ONLY:  message_string

       IMPLICIT NONE

       WRITE(message_string,'(A,I8,A)')                                        &
                         'The maximumn number of agents during this run was',  &
                         maximum_number_of_agents,                             &
                         '&Consider adjusting the INPUT parameter'//           &
                         '&dim_size_agtnum_manual accordingly for the next run.'

       CALL message( 'mas_data_output_agents', 'PA0457', 0, 0, 0, 6, 0 )

    END SUBROUTINE mas_last_actions

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Finds the shortest path from a start position to a target position using the
!> A*-algorithm
!------------------------------------------------------------------------------!
    SUBROUTINE mas_nav_a_star( start_x, start_y, target_x, target_y, nsteps )

       IMPLICIT NONE

       LOGICAL ::  target_reached !< flag

       INTEGER(iwp) ::  cur_node     !< current node of binary heap
       INTEGER(iwp) ::  il           !< counter (x)
       INTEGER(iwp) ::  neigh_node   !< neighbor node
       INTEGER(iwp) ::  node_counter !< binary heap node counter
       INTEGER(iwp) ::  path_ag      !< index of agent path
       INTEGER(iwp) ::  som          !< size of mesh
       INTEGER(iwp) ::  steps        !< steps along the path
       INTEGER(iwp) ::  nsteps        !< number of steps

       REAL(wp) ::  start_x      !< x-coordinate agent
       REAL(wp) ::  start_y      !< y-coordinate agent
       REAL(wp) ::  new_cost     !< updated cost to reach node
       REAL(wp) ::  new_priority !< priority of node to be added to queue
       REAL(wp) ::  rn_gate      !< random number for corner gate
       REAL(wp) ::  target_x     !< x-coordinate target
       REAL(wp) ::  target_y     !< y-coordinate target
!
!--    Coordinate Type
       TYPE coord
          REAL(wp) ::  x   !< x-coordinate
          REAL(wp) ::  x_s !< x-coordinate (shifted)
          REAL(wp) ::  y   !< y-coordinate
          REAL(wp) ::  y_s !< y-coordinate (shifted)
       END TYPE coord

       TYPE(coord), DIMENSION(:), ALLOCATABLE, TARGET ::  path     !< path array
       TYPE(coord), DIMENSION(:), ALLOCATABLE, TARGET ::  tmp_path !< temporary path for resizing

       node_counter = 0
!
!--    Create temporary navigation mesh including agent and target positions
       CALL mas_nav_create_tmp_mesh( start_x, start_y, target_x, target_y, som )
       tmp_mesh(som)%cost_so_far = 0.0_wp
!
!--    Initialize priority queue
       heap_count = 0_iwp
       ALLOCATE(queue(0:100))
       target_reached = .FALSE.
!
!--    Add starting point (agent position) to frontier (the frontier consists
!--    of all the nodes that are to be visited. The node with the smallest
!--    priority will be visited first. The priority consists of the distance
!--    from the start node to this node plus a minimal guess (direct distance)
!--    from this node to the goal). For the starting node, the priority is set
!--    to 0, as it's the only node thus far
       CALL mas_heap_insert_item(som,0.0_wp)
       cur_node = som
       DO WHILE ( heap_count > 0 )
!
!--       Step one: Pick lowest priority item from queue
          node_counter = node_counter + 1
          CALL mas_heap_extract_item(cur_node)
!
!--       Node 0 is the goal node
          IF ( cur_node == 0 ) THEN
             EXIT
          ENDIF
!
!--       Loop over all of cur_node's neighbors
          DO il = 1, tmp_mesh(cur_node)%noc
             neigh_node = tmp_mesh(cur_node)%connected_vertices(il)
!
!--          Check, if the way from the start node to this neigh_node via
!--          cur_node is shorter than the previously found shortest path to it.
!--          If so, replace said cost and add neigh_node to the frontier.
!--          cost_so_far is initialized as 1.d12 so that all found distances
!--          should be smaller.
             new_cost   = tmp_mesh(cur_node)%cost_so_far                       &
                         + tmp_mesh(cur_node)%distance_to_vertex(il)
             IF ( new_cost < tmp_mesh(neigh_node)%cost_so_far ) THEN
                tmp_mesh(neigh_node)%cost_so_far = new_cost
                tmp_mesh(neigh_node)%origin_id   = cur_node
!
!--             Priority in the queue is cost_so_far + heuristic to goal
                new_priority = new_cost                                        &
                              + heuristic(tmp_mesh(neigh_node)%x,              &
                                tmp_mesh(neigh_node)%y, tmp_mesh(0)%x,         &
                                tmp_mesh(0)%y)
                CALL mas_heap_insert_item(neigh_node,new_priority)
             ENDIF
          ENDDO
       ENDDO
!
!--    Add nodes to a path array. To do this, we must backtrack from the target
!--    node to its origin to its origin and so on until an node is reached that
!--    has no origin (%origin_id == -1). This is the starting node.
       DEALLOCATE(queue)
       cur_node = 0
       steps = 0
       ALLOCATE(path(1:100))
       DO WHILE ( cur_node /= -1 )
          steps = steps + 1
!
!--       Resize path array if necessary
          IF ( steps > SIZE(path) ) THEN
             ALLOCATE(tmp_path(1:steps-1))
             tmp_path(1:steps-1) = path(1:steps-1)
             DEALLOCATE(path)
             ALLOCATE(path(1:2*(steps-1)))
             path(1:steps-1) = tmp_path(1:steps-1)
             DEALLOCATE(tmp_path)
          ENDIF
          path(steps)%x = tmp_mesh(cur_node)%x
          path(steps)%y = tmp_mesh(cur_node)%y
          path(steps)%x_s = tmp_mesh(cur_node)%x_s
          path(steps)%y_s = tmp_mesh(cur_node)%y_s
          cur_node = tmp_mesh(cur_node)%origin_id
       ENDDO
!
!--    Add calculated intermittent targets to the path until either the
!--    target or the maximum number of intermittent targets is reached.
!--    Ignore starting point (reduce index by one), it is agent position.
       dummy_path_x = -1
       dummy_path_y = -1
       path_ag = 1
       steps = steps - 1
       nsteps = 0
       DO WHILE( steps > 0 .AND. path_ag <= agt_path_size )
!
!--       Each target point is randomly chosen along a line target along the
!--       bisector of the building corner that starts at corner_gate_start
!--       and has a width of corner_gate_width. This is to avoid clustering
!--       when opposing agent groups try to reach the same corner target.
          rn_gate = random_function(iran_agent) * corner_gate_width            &
                                                + corner_gate_start
          dummy_path_x(path_ag) = path(steps)%x + rn_gate                      &
                         * (path(steps)%x_s - path(steps)%x)
          dummy_path_y(path_ag) = path(steps)%y + rn_gate                      &
                         * (path(steps)%y_s - path(steps)%y)
          steps = steps - 1
          path_ag = path_ag + 1
          nsteps = nsteps + 1
       ENDDO
!
!--    Set current intermittent target of this agent
       DEALLOCATE(tmp_mesh, path)

    END SUBROUTINE mas_nav_a_star

!------------------------------------------------------------------------------! 
! Description:
! ------------
!> Adds a connection between two points of the navigation mesh 
!> (one-way: in_mp1 to in_mp2)
!------------------------------------------------------------------------------! 
    SUBROUTINE mas_nav_add_connection ( in_mp1, id2, in_mp2 )

       IMPLICIT NONE

       LOGICAL ::  connection_established  !< Flag to indicate if connection has already been established

       INTEGER(iwp) ::  id2  !< ID of in_mp2
       INTEGER(iwp) ::  il   !< local counter
       INTEGER(iwp) ::  noc1 !< number of connections in in_mp1

       INTEGER, DIMENSION(:), ALLOCATABLE ::  dum_cv !< dummy array for connected_vertices

       REAL(wp) ::  dist  !< Distance between the two points

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  dum_dtv

       TYPE(mesh_point) ::  in_mp1  !< mesh point that gets a new connection
       TYPE(mesh_point) ::  in_mp2  !< mesh point in_mp1 will be connected to

       connection_established = .FALSE.
!
!--    Check if connection has already been established
       noc1 = SIZE(in_mp1%connected_vertices)
       DO il = 1, in_mp1%noc
          IF ( in_mp1%connected_vertices(il) == id2 ) THEN
             connection_established = .TRUE.
             EXIT
          ENDIF
       ENDDO

       IF ( .NOT. connection_established ) THEN
!
!--       Resize arrays, if necessary
          IF ( in_mp1%noc >= noc1 ) THEN
             ALLOCATE( dum_cv(1:noc1),dum_dtv(1:noc1) )
             dum_cv  = in_mp1%connected_vertices
             dum_dtv = in_mp1%distance_to_vertex
             DEALLOCATE( in_mp1%connected_vertices, in_mp1%distance_to_vertex )
             ALLOCATE( in_mp1%connected_vertices(1:2*noc1),                    &
                       in_mp1%distance_to_vertex(1:2*noc1) )
             in_mp1%connected_vertices         = -999
             in_mp1%distance_to_vertex         = -999.
             in_mp1%connected_vertices(1:noc1) = dum_cv
             in_mp1%distance_to_vertex(1:noc1) = dum_dtv
          ENDIF

!
!--       Add connection
          in_mp1%noc = in_mp1%noc+1
          dist = SQRT( (in_mp1%x - in_mp2%x)**2 + (in_mp1%y - in_mp2%y)**2 )
          in_mp1%connected_vertices(in_mp1%noc) = id2
          in_mp1%distance_to_vertex(in_mp1%noc) = dist
       ENDIF

    END SUBROUTINE mas_nav_add_connection

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Adds a vertex (curren position of agent or target) to the existing tmp_mesh
!------------------------------------------------------------------------------! 
    SUBROUTINE mas_nav_add_vertex_to_mesh ( in_mp, in_id )

       IMPLICIT NONE

       LOGICAL ::  intersection_found !< flag

       INTEGER(iwp) ::  jl    !< mesh point counter
       INTEGER(iwp) ::  pl    !< polygon counter
       INTEGER(iwp) ::  vl    !< vertex counter
       INTEGER(iwp) ::  pid_t !< polygon id of tested mesh point
       INTEGER(iwp) ::  vid_t !< vertex id of tested mesh point
       INTEGER(iwp) ::  in_id !< vertex id of tested mesh point

       REAL(wp) ::  v1x !< x-coordinate of test vertex 1 for intersection test
       REAL(wp) ::  v1y !< y-coordinate of test vertex 1 for intersection test
       REAL(wp) ::  v2x !< x-coordinate of test vertex 2 for intersection test
       REAL(wp) ::  v2y !< y-coordinate of test vertex 2 for intersection test
       REAL(wp) ::  x   !< x-coordinate of current mesh point
       REAL(wp) ::  x_t !< x-coordinate of tested mesh point
       REAL(wp) ::  y   !< y-coordinate of current mesh point
       REAL(wp) ::  y_t !< y-coordinate of tested mesh point

       TYPE(mesh_point) ::  in_mp !< Input mesh point
!
!--
       x = in_mp%x
       y = in_mp%y
       DO jl = 0, SIZE(tmp_mesh)-2
          IF ( in_id == jl ) CYCLE
!
!--       Ignore mesh points with 0 connections
          IF ( tmp_mesh(jl)%polygon_id /= -1 ) THEN
             IF ( tmp_mesh(jl)%noc == 0 ) CYCLE
          ENDIF
          x_t = tmp_mesh(jl)%x
          y_t = tmp_mesh(jl)%y
          pid_t = tmp_mesh(jl)%polygon_id
          vid_t = tmp_mesh(jl)%vertex_id
!
!--       If the connecting line between the target and a mesh point points 
!--       into the mesh point's polygon, no connection will be 
!--       established between the two points. This is the case if the 
!--       previous (next) vertex of the polygon is right of the connecting
!--       line and the next (previous) vertex of the polygon is left of the
!--       connecting line.
          IF ( pid_t > 0 .AND. pid_t <= SIZE(polygons) ) THEN
             IF ( (((is_left(x,y,x_t,y_t,polygons(pid_t)%vertices(vid_t-1)%x,  &
                                       polygons(pid_t)%vertices(vid_t-1)%y)    &
                  .AND. is_right(x,y,x_t,y_t,                                  &
                                       polygons(pid_t)%vertices(vid_t+1)%x,    &
                                       polygons(pid_t)%vertices(vid_t+1)%y) )  &
                  .OR. (is_right(x,y,x_t,y_t,                                  &
                                       polygons(pid_t)%vertices(vid_t-1)%x,    &
                                       polygons(pid_t)%vertices(vid_t-1)%y)    &
                  .AND. is_left(x,y,x_t,y_t,                                   &
                                       polygons(pid_t)%vertices(vid_t+1)%x,    &
                                       polygons(pid_t)%vertices(vid_t+1)%y)))))&
             THEN
                CYCLE
             ENDIF
          ENDIF
!
!--       For each edge of each polygon, check if it intersects with the 
!--       potential connection. If at least one intersection is found,
!--       no connection can be made
          intersection_found = .FALSE.
          DO pl = 1, SIZE(polygons)
             DO vl = 1, polygons(pl)%nov
                v1x = polygons(pl)%vertices(vl)%x
                v1y = polygons(pl)%vertices(vl)%y
                v2x = polygons(pl)%vertices(vl+1)%x
                v2y = polygons(pl)%vertices(vl+1)%y
                intersection_found = intersect(x,y,x_t,y_t,v1x,v1y,v2x,v2y)
                IF ( intersection_found ) THEN
                   EXIT
                ENDIF
             ENDDO
             IF ( intersection_found ) EXIT
          ENDDO
          IF ( intersection_found ) CYCLE
!
!--       If neither of the above two test was true, a connection will be
!--       established between the two mesh points.
          CALL mas_nav_add_connection(in_mp,jl, tmp_mesh(jl))
          CALL mas_nav_add_connection(tmp_mesh(jl),in_id, in_mp)
       ENDDO
       CALL mas_nav_reduce_connections(in_mp)

    END SUBROUTINE mas_nav_add_vertex_to_mesh

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Creates a temporary copy of the navigation mesh to be used for pathfinding
!------------------------------------------------------------------------------!
    SUBROUTINE mas_nav_create_tmp_mesh( a_x, a_y, t_x, t_y, som )

       IMPLICIT NONE

       INTEGER(iwp) ::  som !< size of mesh
       INTEGER(iwp) ::  noc !< number of connetions
       INTEGER(iwp) ::  im  !< local mesh point counter

       REAL(wp) ::  a_x !< x-coordinate agent
       REAL(wp) ::  a_y !< y-coordinate agent
       REAL(wp) ::  t_x !< x-coordinate target
       REAL(wp) ::  t_y !< y-coordinate target
!
!--    give tmp_mesh the size of mesh
       som = SIZE(mesh)+1
       ALLOCATE(tmp_mesh(0:som))
!
!--    give the allocatable variables in tmp_mesh their respctive sizes
       DO im = 1, som-1
          noc = mesh(im)%noc
          ALLOCATE(tmp_mesh(im)%connected_vertices(1:noc))
          ALLOCATE(tmp_mesh(im)%distance_to_vertex(1:noc))
       ENDDO
!
!--    copy mesh to tmp_mesh
       tmp_mesh(1:som-1) = mesh(1:som-1)
!
!--    Add target point ...
       CALL mas_nav_init_mesh_point(tmp_mesh(0),-1_iwp,-1_iwp,t_x, t_y)
       CALL mas_nav_add_vertex_to_mesh(tmp_mesh(0),0_iwp)
!
!--    ... and start point to temp mesh
       CALL mas_nav_init_mesh_point(tmp_mesh(som),-1_iwp,-1_iwp,a_x, a_y)
       CALL mas_nav_add_vertex_to_mesh(tmp_mesh(som),som)

    END SUBROUTINE mas_nav_create_tmp_mesh
    

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Finds the shortest path from an agents' position to her target. As the
!> actual pathfinding algorithm uses the obstacle corners and then shifts them
!> outward after pathfinding, cases can uccur in which the connection between
!> these intermittent targets then intersect with obstacles. To remedy this
!> the pathfinding algorithm is then run on every two subsequent intermittent
!> targets iteratively and new intermittent targets may be added to the path
!> this way.
!------------------------------------------------------------------------------!
    SUBROUTINE mas_nav_find_path( nl )

       IMPLICIT NONE

       INTEGER(iwp) ::  nl            !< local agent counter
       INTEGER(iwp) ::  il            !< local counter
       INTEGER(iwp) ::  jl            !< local counter
       INTEGER(iwp) ::  kl            !< local counter
       INTEGER(iwp) ::  nsteps_total  !< number of steps on path
       INTEGER(iwp) ::  nsteps_dummy  !< number of steps on path
       
       REAL(wp), DIMENSION(0:30) ::  ld_path_x !< local dummy agent path to target (x)
       REAL(wp), DIMENSION(0:30) ::  ld_path_y !< local dummy agent path to target (y)
!
!--    Initialize agent path arrays
       agents(nl)%path_x    = -1
       agents(nl)%path_y    = -1
       agents(nl)%path_x(0) = agents(nl)%x
       agents(nl)%path_y(0) = agents(nl)%y
!
!--    Calculate initial path
       CALL mas_nav_a_star( agents(nl)%x,   agents(nl)%y,                      &
                            agents(nl)%t_x, agents(nl)%t_y, nsteps_total )
!
!--    Set the rest of the agent path that was just calculated
       agents(nl)%path_x(1:nsteps_total) = dummy_path_x(1:nsteps_total)
       agents(nl)%path_y(1:nsteps_total) = dummy_path_y(1:nsteps_total)
!
!--    Iterate through found path and check more intermittent targets need
!--    to be added. For this, run pathfinding between every two consecutive
!--    intermittent targets.
       DO il = 0, MIN(agt_path_size-1, nsteps_total-1)
!
!--       pathfinding between two consecutive intermittent targets
          CALL mas_nav_a_star( agents(nl)%path_x(il),   agents(nl)%path_y(il), &
                              agents(nl)%path_x(il+1), agents(nl)%path_y(il+1),&
                              nsteps_dummy )
          nsteps_dummy = nsteps_dummy - 1
!
!--       If additional intermittent targets are found, add them to the path
          IF ( nsteps_dummy > 0 ) THEN
             ld_path_x = -1
             ld_path_y = -1
             ld_path_x(il+1:il+nsteps_dummy) = dummy_path_x(1:nsteps_dummy)
             ld_path_y(il+1:il+nsteps_dummy) = dummy_path_y(1:nsteps_dummy)
             kl = 1
             DO jl = il+1,nsteps_total
               ld_path_x( il+nsteps_dummy+kl ) = agents(nl)%path_x(jl)
               ld_path_y( il+nsteps_dummy+kl ) = agents(nl)%path_y(jl)
               kl = kl + 1
               IF ( kl > agt_path_size ) EXIT
             ENDDO
             nsteps_total = MIN(nsteps_total + nsteps_dummy, agt_path_size)
             agents(nl)%path_x(il+1:nsteps_total) = ld_path_x(il+1:nsteps_total)
             agents(nl)%path_y(il+1:nsteps_total) = ld_path_y(il+1:nsteps_total)
          ENDIF

       ENDDO
!
!--    reset path counter to first intermittent target
       agents(nl)%path_counter = 1

    END SUBROUTINE mas_nav_find_path

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reduces the size of connection array to the amount of actual connections
!> after all connetions were added to a mesh point
!------------------------------------------------------------------------------!
    SUBROUTINE mas_nav_reduce_connections ( in_mp )

       IMPLICIT NONE

       INTEGER(iwp) ::  noc  !< number of connections

       INTEGER, DIMENSION(:), ALLOCATABLE ::  dum_cv   !< dummy connected_vertices

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  dum_dtv !< dummy distance_to_vertex

       TYPE(mesh_point) ::  in_mp

       noc = in_mp%noc
       ALLOCATE( dum_cv(1:noc),dum_dtv(1:noc) )
       dum_cv  = in_mp%connected_vertices(1:noc)
       dum_dtv = in_mp%distance_to_vertex(1:noc)
       DEALLOCATE( in_mp%connected_vertices, in_mp%distance_to_vertex )
       ALLOCATE( in_mp%connected_vertices(1:noc),                    &
                 in_mp%distance_to_vertex(1:noc) )
       in_mp%connected_vertices(1:noc) = dum_cv(1:noc)
       in_mp%distance_to_vertex(1:noc) = dum_dtv(1:noc)

    END SUBROUTINE mas_nav_reduce_connections

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initializes a point of the navigation mesh
!------------------------------------------------------------------------------!
    SUBROUTINE mas_nav_init_mesh_point ( in_mp, pid, vid, x, y )

       IMPLICIT NONE

       INTEGER(iwp) ::  pid !< polygon ID
       INTEGER(iwp) ::  vid !< vertex ID

       REAL(wp) ::  x !< x-coordinate
       REAL(wp) ::  y !< y-coordinate

       TYPE(mesh_point) ::  in_mp !< mesh point to be initialized

       in_mp%origin_id          = -1
       in_mp%polygon_id         = pid
       in_mp%vertex_id          = vid
       in_mp%cost_so_far        = 1.d12
       in_mp%x                  = x
       in_mp%y                  = y
       in_mp%x_s                = x
       in_mp%y_s                = y
       ALLOCATE(in_mp%connected_vertices(1:100),                               &
                in_mp%distance_to_vertex(1:100))
       in_mp%connected_vertices = -999
       in_mp%distance_to_vertex = -999.
       in_mp%noc                = 0

    END SUBROUTINE mas_nav_init_mesh_point

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reading of namlist from parin file
!------------------------------------------------------------------------------!
    SUBROUTINE mas_parin

       USE control_parameters,                                                 &
           ONLY: agent_time_unlimited, multi_agent_system_end,                 &
                 multi_agent_system_start

       IMPLICIT NONE

       CHARACTER (LEN=80) ::  line  !<

       NAMELIST /agent_parameters/  a_rand_target,                             &
                                    adx,                                       &
                                    ady,                                       &
                                    agent_maximum_age,                         &
                                    agent_time_unlimited,                      &
                                    alloc_factor_mas,                          &
                                    asl,                                       &
                                    asn,                                       &
                                    asr,                                       &
                                    ass,                                       &
                                    at_x,                                      &
                                    at_y,                                      &
                                    bc_mas_lr,                                 &
                                    bc_mas_ns,                                 &
                                    coll_t_0,                                  &
                                    corner_gate_start,                         &
                                    corner_gate_width,                         &
                                    dim_size_agtnum_manual,                    &
                                    dim_size_factor_agtnum,                    &
                                    deallocate_memory_mas,                     &
                                    dist_to_int_target,                        &
                                    dt_agent,                                  &
                                    dt_arel,                                   &
                                    dt_write_agent_data,                       &
                                    end_time_arel,                             &
                                    max_dist_from_path,                        &
                                    min_nr_agent,                              &
                                    multi_agent_system_end,                    &
                                    multi_agent_system_start,                  &
                                    number_of_agent_groups,                    &
                                    radius_agent,                              &
                                    random_start_position_agents,              &
                                    read_agents_from_restartfile,              &
                                    repuls_agent,                              &
                                    repuls_wall,                               &
                                    scan_radius_agent,                         &
                                    sigma_rep_agent,                           &
                                    sigma_rep_wall,                            &
                                    step_dealloc_mas,                          &
                                    tau_accel_agent

!
!--    Try to find agent package
       REWIND ( 11 )
       line = ' '
       DO WHILE ( INDEX( line, '&agent_parameters' ) == 0 )
          READ ( 11, '(A)', END=20 )  line
       ENDDO
       BACKSPACE ( 11 )

!
!--    Read user-defined namelist
       READ ( 11, agent_parameters, ERR = 10, END = 20 )

!
!--    Set flag that indicates that agents are switched on
       agents_active = .TRUE.
       GOTO 20

 10    BACKSPACE( 11 )
       READ( 11 , '(A)') line
       CALL parin_fail_message( 'agent_parameters', line )

 20    CONTINUE

    END SUBROUTINE mas_parin

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Routine for the whole processor
!> Sort all agents into the 4 respective subgrid boxes
!------------------------------------------------------------------------------!
    SUBROUTINE mas_ps_sort_in_subboxes

       IMPLICIT NONE

       INTEGER(iwp) ::  i           !< grid box (x)
       INTEGER(iwp) ::  ip          !< counter (x)
       INTEGER(iwp) ::  is          !< box counter
       INTEGER(iwp) ::  j           !< grid box (y)
       INTEGER(iwp) ::  jp          !< counter (y)
       INTEGER(iwp) ::  m           !< sorting index
       INTEGER(iwp) ::  n           !< agent index
       INTEGER(iwp) ::  nn          !< agent counter
       INTEGER(iwp) ::  sort_index  !< sorting index

       INTEGER(iwp), DIMENSION(0:3) ::  sort_count  !< number of agents in one subbox

       TYPE(agent_type), DIMENSION(:,:), ALLOCATABLE ::  sort_agents  !< sorted agent array

       DO  ip = nxl, nxr
          DO  jp = nys, nyn
             number_of_agents = agt_count(jp,ip)
             IF ( number_of_agents <= 0 )  CYCLE
             agents => grid_agents(jp,ip)%agents(1:number_of_agents)

             nn = 0
             sort_count = 0
             ALLOCATE( sort_agents(number_of_agents, 0:3) )

             DO  n = 1, number_of_agents
                sort_index = 0

                IF ( agents(n)%agent_mask )  THEN
                   nn = nn + 1
!
!--                Sorting agents with a binary scheme
!--                sort_index=11_2=3_10 -> agent at the left,south subgridbox
!--                sort_index=10_2=2_10 -> agent at the left,north subgridbox
!--                sort_index=01_2=1_10 -> agent at the right,south subgridbox
!--                sort_index=00_2=0_10 -> agent at the right,north subgridbox
!--                For this the center of the gridbox is calculated 
                   i = (agents(n)%x + 0.5_wp * dx) * ddx
                   j = (agents(n)%y + 0.5_wp * dy) * ddy

                   IF ( i == ip )  sort_index = sort_index + 2
                   IF ( j == jp )  sort_index = sort_index + 1

                   sort_count(sort_index) = sort_count(sort_index) + 1
                   m = sort_count(sort_index)
                   sort_agents(m,sort_index) = agents(n)
                   sort_agents(m,sort_index)%block_nr = sort_index
                ENDIF
             ENDDO

             nn = 0
             DO is = 0,3
                grid_agents(jp,ip)%start_index(is) = nn + 1
                DO n = 1,sort_count(is)
                   nn = nn + 1
                   agents(nn) = sort_agents(n,is)
                ENDDO
                grid_agents(jp,ip)%end_index(is) = nn
             ENDDO

             number_of_agents = nn
             agt_count(jp,ip) = number_of_agents
             DEALLOCATE(sort_agents)
          ENDDO
       ENDDO

    END SUBROUTINE mas_ps_sort_in_subboxes

#if defined( __parallel )
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Move all agents not marked for deletion to lowest indices (packing)
!------------------------------------------------------------------------------!
    SUBROUTINE mas_ps_pack

       IMPLICIT NONE

       INTEGER(iwp) ::  n  !< agent counter
       INTEGER(iwp) ::  nn !< number of agents
!
!--    Find out elements marked for deletion and move data from highest index
!--    values to these free indices
       nn = number_of_agents

       DO WHILE ( .NOT. agents(nn)%agent_mask )
          nn = nn-1
          IF ( nn == 0 )  EXIT
       ENDDO

       IF ( nn > 0 )  THEN
          DO  n = 1, number_of_agents
             IF ( .NOT. agents(n)%agent_mask )  THEN
                agents(n) = agents(nn)
                nn = nn - 1
                DO WHILE ( .NOT. agents(nn)%agent_mask )
                   nn = nn-1
                   IF ( n == nn )  EXIT
                ENDDO
             ENDIF
             IF ( n == nn )  EXIT
          ENDDO
       ENDIF

!
!--    The number of deleted agents has been determined in routines
!--    mas_boundary_conds, mas_droplet_collision, and mas_eh_exchange_horiz
       number_of_agents = nn

    END SUBROUTINE mas_ps_pack
#endif

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Sort agents in each sub-grid box into two groups: agents that already
!> completed the LES timestep, and agents that need further timestepping to
!> complete the LES timestep.
!------------------------------------------------------------------------------!
!    SUBROUTINE mas_ps_sort_timeloop_done
!
!       IMPLICIT NONE
!
!       INTEGER(iwp) :: end_index     !< agent end index for each sub-box
!       INTEGER(iwp) :: i             !< index of agent grid box in x-direction
!       INTEGER(iwp) :: j             !< index of agent grid box in y-direction
!       INTEGER(iwp) :: n             !< running index for number of agents
!       INTEGER(iwp) :: nb            !< index of subgrid boux
!       INTEGER(iwp) :: nf            !< indices for agents in each sub-box that already finalized their substeps
!       INTEGER(iwp) :: nnf           !< indices for agents in each sub-box that need further treatment
!       INTEGER(iwp) :: num_finalized !< number of agents in each sub-box that already finalized their substeps
!       INTEGER(iwp) :: start_index   !< agent start index for each sub-box
!
!       TYPE(agent_type), DIMENSION(:), ALLOCATABLE :: sort_agents  !< temporary agent array
!
!       DO  i = nxl, nxr
!          DO  j = nys, nyn
!
!             number_of_agents = agt_count(j,i)
!             IF ( number_of_agents <= 0 )  CYCLE
!
!             agents => grid_agents(j,i)%agents(1:number_of_agents)
!
!             DO  nb = 0, 3
!
!--             Obtain start and end index for each subgrid box
!                start_index = grid_agents(j,i)%start_index(nb)
!                end_index   = grid_agents(j,i)%end_index(nb)
!
!--             Allocate temporary array used for sorting
!                ALLOCATE( sort_agents(start_index:end_index) )
!
!--             Determine number of agents already completed the LES 
!--             timestep, and write them into a temporary array
!                nf = start_index
!                num_finalized = 0
!                DO  n = start_index, end_index
!                   IF ( dt_3d - agents(n)%dt_sum < 1E-8_wp )  THEN
!                      sort_agents(nf) = agents(n)
!                      nf              = nf + 1
!                      num_finalized   = num_finalized + 1
!                   ENDIF
!                ENDDO
!
!--             Determine number of agents that not completed the LES 
!--             timestep, and write them into a temporary array
!                nnf = nf
!                DO  n = start_index, end_index
!                   IF ( dt_3d - agents(n)%dt_sum > 1E-8_wp )  THEN
!                      sort_agents(nnf) = agents(n)
!                      nnf              = nnf + 1
!                   ENDIF
!                ENDDO
!
!--             Write back sorted agents
!                agents(start_index:end_index) =                          &
!                                        sort_agents(start_index:end_index)
!
!--             Determine updated start_index, used to masked already 
!--             completed agents. 
!                grid_agents(j,i)%start_index(nb) =                     &
!                                   grid_agents(j,i)%start_index(nb)    &
!                                 + num_finalized
!
!--             Deallocate dummy array
!                DEALLOCATE ( sort_agents )
!
!--             Finally, if number of non-completed agents is non zero 
!--             in any of the sub-boxes, set control flag appropriately. 
!                IF ( nnf > nf )                                             &
!                   grid_agents(j,i)%time_loop_done = .FALSE.
!
!             ENDDO
!          ENDDO
!       ENDDO
!
!    END SUBROUTINE mas_ps_sort_timeloop_done

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calls social forces calculations
!------------------------------------------------------------------------------!
    SUBROUTINE mas_timestep_forces_call ( ip, jp )

       IMPLICIT NONE

       INTEGER(iwp) ::  ip  !< counter, x-direction
       INTEGER(iwp) ::  jp  !< counter, y-direction
       INTEGER(iwp) ::  n   !< loop variable over all agents in a grid box

!
!--    Get direction for all agents in current grid cell
       CALL mas_agent_direction

       DO n = 1, number_of_agents

          force_x = 0.0_wp
          force_y = 0.0_wp

          CALL mas_timestep_social_forces ( 'acceleration', n, ip, jp )

          CALL mas_timestep_social_forces ( 'other_agents', n, ip, jp )

          CALL mas_timestep_social_forces ( 'walls',        n, ip, jp )
!
!--       Update forces
          agents(n)%force_x = force_x
          agents(n)%force_y = force_y
       ENDDO

    END SUBROUTINE mas_timestep_forces_call

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Euler timestep of agent transport
!------------------------------------------------------------------------------!
    SUBROUTINE mas_timestep

       IMPLICIT NONE

       INTEGER(iwp) ::  n !< loop variable over all agents in a grid box

       REAL(wp) ::  abs_v !< absolute value of velocity
       REAL(wp) ::  abs_f !< absolute value of force

       DO n = 1, number_of_agents
!
!--       Limit absolute force to a maximum to prevent unrealistic acceleration
          abs_f = SQRT((agents(n)%force_x)**2 + (agents(n)%force_y)**2)
          IF ( abs_f > 20. ) THEN
             agents(n)%force_x = agents(n)%force_x * 20. / abs_f
             agents(n)%force_y = agents(n)%force_y * 20. / abs_f
          ENDIF
!
!--       Update agent speed
          agents(n)%speed_x = agents(n)%speed_x + agents(n)%force_x * dt_agent
          agents(n)%speed_y = agents(n)%speed_y + agents(n)%force_y * dt_agent
!
!--       Reduction of agent speed to maximum agent speed
          abs_v = SQRT((agents(n)%speed_x)**2 + (agents(n)%speed_y)**2)
          IF ( abs_v > v_max_agent ) THEN
             agents(n)%speed_x = agents(n)%speed_x * v_max_agent / abs_v
             agents(n)%speed_y = agents(n)%speed_y * v_max_agent / abs_v
          ENDIF
!
!--       Update agent position
          agents(n)%x = agents(n)%x + agents(n)%speed_x * dt_agent
          agents(n)%y = agents(n)%y + agents(n)%speed_y * dt_agent
!
!--       Update absolute value of agent speed
          agents(n)%speed_abs = abs_v
!
!--       Increment the agent age and the total time that the agent
!--       has advanced within the agent timestep procedure
          agents(n)%age_m  = agents(n)%age
          agents(n)%age    = agents(n)%age    + dt_agent
          agents(n)%dt_sum = agents(n)%dt_sum + dt_agent
!
!--       Check whether there is still an agent that has not yet completed 
!--       the total LES timestep
          IF ( ( dt_3d - agents(n)%dt_sum ) > 1E-8_wp )  THEN
             dt_3d_reached_l_mas = .FALSE.
          ENDIF

       ENDDO

    END SUBROUTINE mas_timestep

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculates the Social Forces (Helbing and Molnar, 1995) that the agent
!> experiences due to acceleration towards target and repulsion by obstacles
!------------------------------------------------------------------------------!
    SUBROUTINE mas_timestep_social_forces ( mode, nl, ip, jp )

       IMPLICIT NONE

       CHARACTER (LEN=*) ::  mode  !< identifier for the mode of calculation

       INTEGER(iwp) ::  ij_dum      !< index of nearest wall
       INTEGER(iwp) ::  il          !< index variable along x 
       INTEGER(iwp) ::  ip          !< index variable along x 
       INTEGER(iwp) ::  jl          !< index variable along y
       INTEGER(iwp) ::  jp          !< index variable along y
       INTEGER(iwp) ::  nl          !< loop variable over all agents in a grid box
       INTEGER(iwp) ::  no          !< loop variable over all agents in a grid box
       INTEGER(iwp) ::  noa         !< amount of agents in a grid box 
       INTEGER(iwp) ::  sc_x_end    !< index for scan for topography/other agents
       INTEGER(iwp) ::  sc_x_start  !< index for scan for topography/other agents
       INTEGER(iwp) ::  sc_y_end    !< index for scan for topography/other agents
       INTEGER(iwp) ::  sc_y_start  !< index for scan for topography/other agents

       LOGICAL ::  corner_found  !< flag that indicates a corner has been found near agent

       REAL(wp) ::  a_pl             !< factor for collision avoidance
       REAL(wp) ::  ax_semimaj       !< semiminor axis of repulsive ellipse
       REAL(wp) ::  b_pl             !< factor for collision avoidance
       REAL(wp) ::  c_pl             !< factor for collision avoidance
       REAL(wp) ::  coll_t           !< time at which the next collision would happen
       REAL(wp) ::  d_coll_t_0       !< inverse of collision cutoff time
       REAL(wp) ::  d_pl             !< factor for collision avoidance
       REAL(wp) ::  ddum_f           !< dummy devisor collision avoidance
       REAL(wp) ::  dist             !< distance to obstacle
       REAL(wp) ::  dist_sq          !< distance to obstacle squared
       REAL(wp) ::  pos_rel_x        !< relative position of two agents (x)
       REAL(wp) ::  pos_rel_y        !< relative position of two agents (y)
       REAL(wp) ::  r_sq             !< y-position
       REAL(wp) ::  sra              !< scan radius (agents)
       REAL(wp) ::  srw              !< local variable for scan radius (walls)
       REAL(wp) ::  v_rel_x          !< relative velocity (x); collision avoidance
       REAL(wp) ::  v_rel_y          !< relative velocity (y); collision avoidance
       REAL(wp) ::  x_a              !< x-position
       REAL(wp) ::  x_wall           !< x-position of wall
       REAL(wp) ::  y_a              !< y-position
       REAL(wp) ::  y_wall           !< y-position of wall

       REAL(wp), PARAMETER ::  k_pl = 1.5  !< factor for collision avoidance

       TYPE(agent_type), DIMENSION(:), POINTER ::  l_agts !< agents that repulse current agent

!
!--    Initialization
       x_a = agents(nl)%x
       y_a = agents(nl)%y

       SELECT CASE ( TRIM( mode ) )
!
!--       Calculation of force due to agent trying to approach desired velocity
          CASE ( 'acceleration' )

             force_x = force_x + d_tau_accel_agent                             &
                          * ( agents(nl)%speed_des*agents(nl)%speed_e_x        &
                             -agents(nl)%speed_x )

             force_y = force_y + d_tau_accel_agent                             &
                          * ( agents(nl)%speed_des*agents(nl)%speed_e_y        &
                             -agents(nl)%speed_y )

!
!--       Calculation of repulsive forces by other agents in a radius around the
!--       current one
          CASE ( 'other_agents' )

             sra = scan_radius_agent
             d_coll_t_0 = 1./coll_t_0
!
!--          Find relevant gridboxes (those that could contain agents within
!--          scan radius)
             sc_x_start = FLOOR( (x_a - sra) * ddx )
             sc_x_end   = FLOOR( (x_a + sra) * ddx )
             sc_y_start = FLOOR( (y_a - sra) * ddx )
             sc_y_end   = FLOOR( (y_a + sra) * ddx )
             IF ( sc_x_start < nxlg ) sc_x_start = nxlg
             IF ( sc_x_end   > nxrg ) sc_x_end   = nxrg
             IF ( sc_y_start < nysg ) sc_y_start = nysg
             IF ( sc_y_end   > nyng ) sc_y_end   = nyng

             sra = sra**2
!
!--          Loop over all previously found relevant gridboxes
             DO il = sc_x_start, sc_x_end
                DO jl = sc_y_start, sc_y_end
                   noa = agt_count(jl,il)
                   IF ( noa <= 0 )  CYCLE
                   l_agts => grid_agents(jl,il)%agents(1:noa)
                   DO no = 1, noa
!
!--                   Skip self
                      IF ( jl == jp .AND. il == ip .AND. no == nl ) CYCLE
                      pos_rel_x = l_agts(no)%x - x_a
                      pos_rel_y = l_agts(no)%y - y_a
                      dist_sq = pos_rel_x**2 + pos_rel_y**2
                      IF ( dist_sq > sra ) CYCLE
                      r_sq    = (2*radius_agent)**2
                      v_rel_x   = agents(nl)%speed_x - l_agts(no)%speed_x
                      v_rel_y   = agents(nl)%speed_y - l_agts(no)%speed_y
!
!--                   Collision is already occuring, default to standard
!--                   social forces
                      IF ( dist_sq <= r_sq ) THEN
                         dist = SQRT(dist_sq) + 1.0d-12
                         ax_semimaj = .5_wp*SQRT( dist )

                         force_x = force_x - 0.125_wp * repuls_agent           &
                                        * d_sigma_rep_agent / ax_semimaj       &
                                        * EXP( -ax_semimaj*d_sigma_rep_agent ) &
                                        * (pos_rel_x/dist)

                         force_y = force_y - 0.125_wp * repuls_agent           &
                                        * d_sigma_rep_agent / ax_semimaj       &
                                        * EXP( -ax_semimaj*d_sigma_rep_agent ) &
                                        * (pos_rel_y/dist)
!
!--                   Currently no collision, calculate collision avoidance
!--                   force according to Karamouzas et al (2014, PRL 113,238701)
                      ELSE
!
!--                     factors
                         a_pl = v_rel_x**2 +  v_rel_y**2
                         b_pl = pos_rel_x*v_rel_x + pos_rel_y*v_rel_y
                         c_pl = dist_sq - r_sq
                         d_pl = b_pl**2 - a_pl*c_pl
!
!--                      If the two agents are moving non-parallel, calculate
!--                      collision avoidance social force
                         IF ( d_pl > 0.0_wp .AND.                              &
                            ( a_pl < -0.00001 .OR. a_pl > 0.00001 ) )          &
                         THEN

                            d_pl   = SQRT(d_pl)
                            coll_t = (b_pl - d_pl)/a_pl
                            IF ( coll_t > 0.0_wp ) THEN
!
!--                            Dummy factor
                               ddum_f = 1. / ( a_pl * coll_t**2 )              &
                                           * ( 2. / coll_t + 1.0 * d_coll_t_0 )
!
!--                            x-component of social force
                               force_x = force_x - k_pl *                      &
                                         EXP( -coll_t * d_coll_t_0 ) *         &
                                         ( v_rel_x -                           &
                                           ( b_pl * v_rel_x -                  &
                                             a_pl * pos_rel_x ) / d_pl ) *     &
                                         ddum_f
!
!--                            y-component of social force
                               force_y = force_y - k_pl *                      &
                                         EXP( -coll_t * d_coll_t_0 ) *         &
                                         ( v_rel_y -                           &
                                           ( b_pl * v_rel_y -                  &
                                             a_pl * pos_rel_y ) / d_pl ) *     &
                                         ddum_f

                            ENDIF
                         ENDIF
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO

          CASE ( 'walls' )

             srw = scan_radius_wall
             corner_found = .FALSE.
!
!--          find relevant grid boxes (those that could contain topography
!--          within radius)
             sc_x_start = (x_a - srw) * ddx
             sc_x_end   = (x_a + srw) * ddx
             sc_y_start = (y_a - srw) * ddx
             sc_y_end   = (y_a + srw) * ddx
             IF ( sc_x_start < nxlg ) sc_x_start = nxlg
             IF ( sc_x_end   > nxrg ) sc_x_end   = nxrg
             IF ( sc_y_start < nysg ) sc_y_start = nysg
             IF ( sc_y_end   > nyng ) sc_y_end   = nyng
!
!--          Find "walls" ( i.e. topography steps (up or down) higher than one
!--          grid box ) that are perpendicular to the agent within the defined
!--          search radius. Such obstacles cannot be passed and a social force
!--          to that effect is applied.
!--          Walls only apply a force perpendicular to the wall to the agent.
!--          There is therefore a search for walls directly right, left, south 
!--          and north of the agent. All other walls are ignored.
!--
!--          Check for wall left of current agent
             ij_dum = 0
             IF ( sc_x_start < ip ) THEN
                DO il = ip - 1, sc_x_start, -1
!
!--                Going left from the agent, check for a right wall
                   IF ( BTEST( obstacle_flags(jp,il), 2 ) ) THEN
!
!--                   obstacle found in grid box il, wall at right side
                      x_wall = (il+1)*dx
!
!--                   Calculate force of found wall on agent
                      CALL mas_timestep_wall_corner_force( x_a, x_wall, y_a,   &
                                                           y_a )
!
!--                   calculate new x starting index for later scan for corners
                      ij_dum = il + 1
                      EXIT
                   ENDIF
                ENDDO
             ENDIF
             IF ( ij_dum /= 0 ) sc_x_start = ij_dum 

!
!--          Check for wall right of current agent
             ij_dum = 0
             IF ( sc_x_end > ip ) THEN
                DO il = ip + 1, sc_x_end
!
!--                Going right from the agent, check for a left wall
                   IF ( BTEST( obstacle_flags(jp,il), 6 ) ) THEN
!
!--                   obstacle found in grid box il, wall at left side
                      x_wall = il*dx
!
!--                   Calculate force of found wall on agent
                      CALL mas_timestep_wall_corner_force( x_a, x_wall, y_a,   &
                                                           y_a )
!
!--                   calculate new x end index for later scan for corners
                      ij_dum = il - 1
                      EXIT
                   ENDIF
                ENDDO
             ENDIF
             IF ( ij_dum /= 0 ) sc_x_end = ij_dum 

!
!--          Check for wall south of current agent
             ij_dum = 0
             IF ( sc_y_start < jp ) THEN
                DO jl = jp - 1, sc_y_start, -1
!
!--                Going south from the agent, check for a north wall
                   IF ( BTEST( obstacle_flags(jl,ip), 0 ) ) THEN
!
!--                   obstacle found in grid box jl, wall at left side
                      y_wall = (jl+1)*dy

                      CALL mas_timestep_wall_corner_force( x_a, x_a, y_a,      &
                                                           y_wall )
!
!--                   calculate new y starting index for later scan for corners
                      ij_dum = jl + 1
                      EXIT
                   ENDIF
                ENDDO
             ENDIF
             IF ( ij_dum /= 0 ) sc_y_start = ij_dum 

!
!--          Check for wall north of current agent
             ij_dum = 0
             IF ( sc_y_end > jp ) THEN
                DO jl = jp + 1, sc_y_end 
!
!--                Going north from the agent, check for a south wall
                   IF ( BTEST( obstacle_flags(jl,ip), 4 ) ) THEN
!
!--                   obstacle found in grid box jl, wall at left side
                      y_wall = jl*dy

                      CALL mas_timestep_wall_corner_force( x_a, x_a, y_a,      &
                                                           y_wall )
!
!--                   calculate new y end index for later scan for corners
                      ij_dum = jl - 1
                   ENDIF
                ENDDO
             ENDIF
             IF ( ij_dum /= 0 ) sc_y_end = ij_dum 

!
!--          Scan for corners surrounding current agent.
!--          Only gridcells that are closer than the closest wall in each
!--          direction (n,s,r,l) are considered in the search since those
!--          further away would have a significantly smaller resulting force
!--          than the closer wall. 
             DO il = sc_x_start, sc_x_end
                DO jl = sc_y_start, sc_y_end
                   IF ( il == ip .OR. jl == jp ) CYCLE
!
!--                corners left of agent
                   IF ( il < ip ) THEN
!
!--                   south left quadrant: look for north right corner
                      IF ( jl < jp ) THEN
                         IF ( BTEST( obstacle_flags(jl,il), 1 ) ) THEN
!
!--                         calculate coordinates of the found corner
                            x_wall = (il+1)*dx
                            y_wall = (jl+1)*dy

                            CALL mas_timestep_wall_corner_force( x_a, x_wall,  &
                                                                 y_a, y_wall )

                         ENDIF
!
!--                   north left quadrant: look for south right corner
                      ELSEIF ( jl > jp ) THEN
                         IF ( BTEST( obstacle_flags(jl,il), 3 ) ) THEN
!
!--                         calculate coordinates of the corner of said gridcell
!--                         that is closest to the current agent
                            x_wall = (il+1)*dx
                            y_wall = jl*dy

                            CALL mas_timestep_wall_corner_force( x_a, x_wall,  &
                                                                 y_a, y_wall )

                         ENDIF
                      ENDIF
                   ELSEIF ( il > ip ) THEN
!
!--                   south right quadrant: look for north left corner
                      IF ( jl < jp ) THEN
                         IF ( BTEST( obstacle_flags(jl,il), 7 ) ) THEN
!
!--                         calculate coordinates of the corner of said gridcell
!--                         that is closest to the current agent
                            x_wall = il*dx
                            y_wall = (jl+1)*dy

                            CALL mas_timestep_wall_corner_force( x_a, x_wall,  &
                                                                 y_a, y_wall )

                         ENDIF
!
!--                   north right quadrant: look for south left corner
                      ELSEIF ( jl > jp ) THEN
                         IF ( BTEST( obstacle_flags(jl,il), 5 ) ) THEN
!
!--                         calculate coordinates of the corner of said gridcell
!--                         that is closest to the current agent
                            x_wall = il*dx
                            y_wall = jl*dy

                            CALL mas_timestep_wall_corner_force( x_a, x_wall,  &
                                                                 y_a, y_wall )

                         ENDIF
                      ENDIF
                   ENDIF
                ENDDO
             ENDDO

          CASE DEFAULT

       END SELECT

    END SUBROUTINE mas_timestep_social_forces

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Given a distance to the current agent, calculates the force a found corner
!> or wall exerts on that agent
!------------------------------------------------------------------------------!
    SUBROUTINE mas_timestep_wall_corner_force( xa, xw, ya, yw )

       IMPLICIT NONE

       REAL(wp) ::  dist_l     !< distance to obstacle
       REAL(wp) ::  force_d_x  !< increment of social force, x-direction
       REAL(wp) ::  force_d_y  !< increment of social force, x-direction
       REAL(wp) ::  xa         !< x-position of agent
       REAL(wp) ::  xw         !< x-position of wall
       REAL(wp) ::  ya         !< x-position of agent
       REAL(wp) ::  yw         !< y-position of wall

       force_d_x = 0.0_wp
       force_d_y = 0.0_wp
!
!--    calculate coordinates of corner relative to agent 
!--    postion and distance between corner and agent
       xw = xa - xw
       yw = ya - yw
       dist_l = SQRT( (xw)**2 + (yw)**2 )
!
!--    calculate x and y component of repulsive force 
!--    induced by previously found corner
       IF ( dist_l > 0 ) THEN
          force_d_x = repuls_wall * d_sigma_rep_wall         &
                      * EXP( -dist_l * d_sigma_rep_wall )      &
                      * xw / (dist_l)
          force_d_y = repuls_wall * d_sigma_rep_wall         &
                      * EXP( -dist_l * d_sigma_rep_wall )      &
                      * yw / (dist_l)
       ENDIF

! !--    forces that are located outside of a sight radius of 
! !--    200 degrees (-> COS(100./180.*pi) = COS(.555*pi)) of
! !--    current agent are considered to have an effect of 50%
!        IF ( force_d_x * agents(nl)%speed_e_x +               &
!             force_d_y * agents(nl)%speed_e_y <               &
!             SQRT(force_d_x**2 + force_d_y**2) *              &
!             COS( .55555555 * 3.1415 ) )                      &
!        THEN
!           force_d_x = force_d_x * .5_wp
!           force_d_y = force_d_y * .5_wp
!        ENDIF

!
!--    add force increment to total force of current agent
       force_x = force_x + force_d_x
       force_y = force_y + force_d_y

    END SUBROUTINE mas_timestep_wall_corner_force

!
!-- Calculates distance of point P to edge (A,B). If A = B, calculates
!-- point-to-point distance from A/B to P
    FUNCTION dist_point_to_edge ( a_x, a_y, b_x, b_y, p_x, p_y )

       IMPLICIT NONE

       REAL(wp)  :: ab_x                !< x-coordinate of vector from A to B
       REAL(wp)  :: ab_y                !< y-coordinate of vector from A to B
       REAL(wp)  :: ab_d                !< inverse length of vector from A to B
       REAL(wp)  :: ab_u_x              !< x-coordinate of vector with direction of ab and length 1
       REAL(wp)  :: ab_u_y              !< y-coordinate of vector with direction of ab and length 1
       REAL(wp)  :: ba_x                !< x-coordinate of vector from B to A
       REAL(wp)  :: ba_y                !< y-coordinate of vector from B to A
       REAL(wp)  :: ap_x                !< x-coordinate of vector from A to P
       REAL(wp)  :: ap_y                !< y-coordinate of vector from A to P
       REAL(wp)  :: bp_x                !< x-coordinate of vector from B to P
       REAL(wp)  :: bp_y                !< y-coordinate of vector from B to P
       REAL(wp)  :: a_x                 !< x-coordinate of point A of edge
       REAL(wp)  :: a_y                 !< y-coordinate of point A of edge
       REAL(wp)  :: b_x                 !< x-coordinate of point B of edge
       REAL(wp)  :: b_y                 !< y-coordinate of point B of edge
       REAL(wp)  :: p_x                 !< x-coordinate of point P
       REAL(wp)  :: p_y                 !< y-coordinate of point P
       REAL(wp)  :: dist_x              !< x-coordinate of point P
       REAL(wp)  :: dist_y              !< y-coordinate of point P
       REAL(wp)  :: dist_point_to_edge  !< y-coordinate of point P

       ab_x = - a_x + b_x
       ab_y = - a_y + b_y
       ba_x = - b_x + a_x 
       ba_y = - b_y + a_y 
       ap_x = - a_x + p_x
       ap_y = - a_y + p_y
       bp_x = - b_x + p_x
       bp_y = - b_y + p_y

       IF ( ab_x * ap_x + ab_y * ap_y <= 0. ) THEN
          dist_point_to_edge = SQRT((a_x - p_x)**2 + (a_y - p_y)**2)
       ELSEIF ( ba_x * bp_x + ba_y * bp_y <= 0. ) THEN
          dist_point_to_edge = SQRT((b_x - p_x)**2 + (b_y - p_y)**2)
       ELSE
          ab_d = 1./SQRT((ab_x)**2+(ab_y)**2)
          ab_u_x = ab_x*ab_d
          ab_u_y = ab_y*ab_d
          dist_x = ap_x - (ap_x*ab_u_x+ap_y*ab_u_y)*ab_u_x
          dist_y = ap_y - (ap_x*ab_u_x+ap_y*ab_u_y)*ab_u_y
          dist_point_to_edge = SQRT( dist_x**2 + dist_y**2 )
       ENDIF

    END FUNCTION dist_point_to_edge

!
!-- Returns the heuristic between points A and B (currently the straight 
!-- distance)
    FUNCTION heuristic ( ax, ay, bx, by )

       IMPLICIT NONE

       REAL(wp)  :: ax           !< x-coordinate of point A
       REAL(wp)  :: ay           !< y-coordinate of point A
       REAL(wp)  :: bx           !< x-coordinate of point B
       REAL(wp)  :: by           !< y-coordinate of point B
       REAL(wp)  :: heuristic    !< return value

       heuristic = SQRT(( ax - bx )**2 + ( ay - by )**2)

    END FUNCTION heuristic 

!
!-- Calculates if point P is left of the infinite
!-- line that contains A and B (direction: A to B) 
!-- Concept: 2D rotation of two vectors
    FUNCTION is_left ( ax, ay, bx, by, px, py )

       IMPLICIT NONE

       LOGICAL  :: is_left !< return value; TRUE if P is left of AB

       REAL(wp)  :: ax     !< x-coordinate of point A
       REAL(wp)  :: ay     !< y-coordinate of point A
       REAL(wp)  :: bx     !< x-coordinate of point B
       REAL(wp)  :: by     !< y-coordinate of point B
       REAL(wp)  :: px     !< x-coordinate of point P
       REAL(wp)  :: py     !< y-coordinate of point P

       is_left = (bx-ax)*(py-ay)-(px-ax)*(by-ay) > 0
       IF ( (ABS(ax-px) < .001 .AND. ABS(ay-py) < .001) .OR.                  &
            (ABS(bx-px) < .001 .AND. ABS(by-py) < .001) )                     &
       THEN
          is_left = .FALSE.
       ENDIF

       RETURN

    END FUNCTION is_left 

!
!-- Calculates if point P is right of the infinite
!-- line that contains A and B (direction: A to B) 
!-- Concept: 2D rotation of two vectors
    FUNCTION is_right ( ax, ay, bx, by, px, py )

       IMPLICIT NONE

       LOGICAL  :: is_right !< return value; TRUE if P is right of AB

       REAL(wp), INTENT(IN)  :: ax     !< x-coordinate of point A
       REAL(wp), INTENT(IN)  :: ay     !< y-coordinate of point A
       REAL(wp), INTENT(IN)  :: bx     !< x-coordinate of point B
       REAL(wp), INTENT(IN)  :: by     !< y-coordinate of point B
       REAL(wp), INTENT(IN)  :: px     !< x-coordinate of point P
       REAL(wp), INTENT(IN)  :: py     !< y-coordinate of point P

       is_right = (bx-ax)*(py-ay)-(px-ax)*(by-ay) < 0
       IF ( (ABS(ax-px) < .001 .AND. ABS(ay-py) < .001) .OR.                  &
            (ABS(bx-px) < .001 .AND. ABS(by-py) < .001) )                     &
       THEN
          is_right = .FALSE.
       ENDIF

       RETURN

    END FUNCTION is_right 

!
!-- Returns true if the line segments AB and PQ share an intersection
    FUNCTION intersect ( ax, ay, bx, by, px, py, qx, qy )

       IMPLICIT NONE

       LOGICAL  :: intersect !< return value; TRUE if intersection was found
       LOGICAL  :: la        !< T if a is left of PQ
       LOGICAL  :: lb        !< T if b is left of PQ
       LOGICAL  :: lp        !< T if p is left of AB
       LOGICAL  :: lq        !< T if q is left of AB
       LOGICAL  :: poss      !< flag that indicates if an intersection is still possible
       LOGICAL  :: ra        !< T if a is right of PQ
       LOGICAL  :: rb        !< T if b is right of PQ
       LOGICAL  :: rp        !< T if p is right of AB
       LOGICAL  :: rq        !< T if q is right of AB

       REAL(wp)  :: ax     !< x-coordinate of point A
       REAL(wp)  :: ay     !< y-coordinate of point A
       REAL(wp)  :: bx     !< x-coordinate of point B
       REAL(wp)  :: by     !< y-coordinate of point B
       REAL(wp)  :: px     !< x-coordinate of point P
       REAL(wp)  :: py     !< y-coordinate of point P
       REAL(wp)  :: qx     !< x-coordinate of point Q
       REAL(wp)  :: qy     !< y-coordinate of point Q

       intersect = .FALSE.
       poss      = .FALSE.
!
!--    Intersection is possible only if P and Q are on opposing sides of AB
       lp = is_left(ax,ay,bx,by,px,py)
       rq = is_right(ax,ay,bx,by,qx,qy)
       IF ( lp .AND. rq ) poss = .TRUE.
       IF ( .NOT. poss ) THEN
          lq = is_left(ax,ay,bx,by,qx,qy)
          rp = is_right(ax,ay,bx,by,px,py)
          IF ( lq .AND. rp ) poss = .TRUE.
       ENDIF
!
!--    Intersection occurs only if above test (poss) was true AND
!--    A and B are on opposing sides of PQ
       IF ( poss ) THEN
          la = is_left(px,py,qx,qy,ax,ay)
          rb = is_right(px,py,qx,qy,bx,by)
          IF ( la .AND. rb ) intersect = .TRUE.
          IF ( .NOT. intersect ) THEN
             lb = is_left(px,py,qx,qy,bx,by)
             ra = is_right(px,py,qx,qy,ax,ay)
             IF ( lb .AND. ra ) intersect = .TRUE.
          ENDIF
       ENDIF

       RETURN

    END FUNCTION intersect 

!
!-- Gives a nuber randomly distributed around an average
    FUNCTION random_normal ( avg, variation )

       IMPLICIT NONE

       REAL(wp)  :: avg            !< x-coordinate of vector from A to B
       REAL(wp)  :: variation      !< y-coordinate of vector from A to B
       REAL(wp)  :: random_normal  !< y-coordinate of vector from A to B

       REAL(wp), DIMENSION(12)  :: random_arr  !< inverse length of vector from A to B

       CALL RANDOM_NUMBER(random_arr)
       random_normal = avg + variation*(SUM(random_arr)-6.)

    END FUNCTION random_normal


 END MODULE multi_agent_system_mod
