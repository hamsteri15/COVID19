!> @file dynamics_mod.f90
!--------------------------------------------------------------------------------------------------!
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
!--------------------------------------------------------------------------------------------------!
!
! Current revisions:
! -----------------
! 
! 
! Former revisions:
! -----------------
! $Id: dynamics_mod.f90 4360 2020-01-07 11:25:50Z suehring $
! Bugfix for last commit.
! 
! 4359 2019-12-30 13:36:50Z suehring
! Refine post-initialization check for realistically inital values of mixing ratio. Give an error
! message for faulty initial values, but only a warning in a restart run. 
! 
! 4347 2019-12-18 13:18:33Z suehring
! Implement post-initialization check for realistically inital values of mixing ratio
! 
! 4281 2019-10-29 15:15:39Z schwenkel
! Moved boundary conditions in dynamics module
! 
! 4097 2019-07-15 11:59:11Z suehring
! Avoid overlong lines - limit is 132 characters per line
!
! 4047 2019-06-21 18:58:09Z knoop
! Initial introduction of the dynamics module with only dynamics_swap_timelevel implemented
!
!
! Description:
! ------------
!> This module contains the dynamics of PALM.
!--------------------------------------------------------------------------------------------------!
 MODULE dynamics_mod


    USE arrays_3d, &
        ONLY:  c_u, c_u_m, c_u_m_l, c_v, c_v_m, c_v_m_l, c_w, c_w_m, c_w_m_l,  &
               dzu, &
               exner, &
               hyp, &
               pt, pt_1, pt_2, pt_init, pt_p, &
               q, q_1, q_2, q_p, &
               s, s_1, s_2, s_p, &
               u, u_1, u_2, u_init, u_p, u_m_l, u_m_n, u_m_r, u_m_s, &
               v, v_1, v_2, v_p, v_init, v_m_l, v_m_n, v_m_r, v_m_s, &
               w, w_1, w_2, w_p, w_m_l, w_m_n, w_m_r, w_m_s

    USE basic_constants_and_equations_mod,                                                         &
        ONLY:  magnus,                                                                             &
               rd_d_rv

    USE control_parameters, &
        ONLY:  bc_dirichlet_l, &
               bc_dirichlet_s, &
               bc_radiation_l, &
               bc_radiation_n, &
               bc_radiation_r, &
               bc_radiation_s, &
               bc_pt_t_val, &
               bc_q_t_val, &
               bc_s_t_val, &
               child_domain, &
               coupling_mode, &
               dt_3d, &
               ibc_pt_b, &
               ibc_pt_t, &
               ibc_q_b, &
               ibc_q_t, &
               ibc_s_b, &
               ibc_s_t, &
               ibc_uv_b, &
               ibc_uv_t, &
               initializing_actions, &
               intermediate_timestep_count, &
               length, &
               message_string, &
               nesting_offline, &
               nudging, &
               restart_string, &
               humidity, &
               neutral, &
               passive_scalar, &
               tsc, &
               use_cmax

    USE grid_variables, &
        ONLY:  ddx, &
               ddy, &
               dx, &
               dy

    USE indices, &
        ONLY:  nbgp, &
               nx, &
               nxl, &
               nxlg, &
               nxr, &
               nxrg, &
               ny, &
               nys, &
               nysg, &
               nyn, &
               nyng, &
               nzb, &
               nzt

    USE kinds

    USE pegrid

    USE pmc_interface, &
        ONLY : nesting_mode

    USE surface_mod, &
        ONLY :  bc_h


    IMPLICIT NONE

    LOGICAL ::  dynamics_module_enabled = .FALSE.   !<

    SAVE

    PRIVATE

!
!-- Public functions
    PUBLIC &
       dynamics_parin, &
       dynamics_check_parameters, &
       dynamics_check_data_output_ts, &
       dynamics_check_data_output_pr, &
       dynamics_check_data_output, &
       dynamics_init_masks, &
       dynamics_define_netcdf_grid, &
       dynamics_init_arrays, &
       dynamics_init, &
       dynamics_init_checks, &
       dynamics_header, &
       dynamics_actions, &
       dynamics_non_advective_processes, &
       dynamics_exchange_horiz, &
       dynamics_prognostic_equations, &
       dynamics_boundary_conditions, &
       dynamics_swap_timelevel, &
       dynamics_3d_data_averaging, &
       dynamics_data_output_2d, &
       dynamics_data_output_3d, &
       dynamics_statistics, &
       dynamics_rrd_global, &
       dynamics_rrd_local, &
       dynamics_wrd_global, &
       dynamics_wrd_local, &
       dynamics_last_actions

!
!-- Public parameters, constants and initial values
    PUBLIC &
       dynamics_module_enabled

    INTERFACE dynamics_parin
       MODULE PROCEDURE dynamics_parin
    END INTERFACE dynamics_parin

    INTERFACE dynamics_check_parameters
       MODULE PROCEDURE dynamics_check_parameters
    END INTERFACE dynamics_check_parameters

    INTERFACE dynamics_check_data_output_ts
       MODULE PROCEDURE dynamics_check_data_output_ts
    END INTERFACE dynamics_check_data_output_ts

    INTERFACE dynamics_check_data_output_pr
       MODULE PROCEDURE dynamics_check_data_output_pr
    END INTERFACE dynamics_check_data_output_pr

    INTERFACE dynamics_check_data_output
       MODULE PROCEDURE dynamics_check_data_output
    END INTERFACE dynamics_check_data_output

    INTERFACE dynamics_init_masks
       MODULE PROCEDURE dynamics_init_masks
    END INTERFACE dynamics_init_masks

    INTERFACE dynamics_define_netcdf_grid
       MODULE PROCEDURE dynamics_define_netcdf_grid
    END INTERFACE dynamics_define_netcdf_grid

    INTERFACE dynamics_init_arrays
       MODULE PROCEDURE dynamics_init_arrays
    END INTERFACE dynamics_init_arrays

    INTERFACE dynamics_init
       MODULE PROCEDURE dynamics_init
    END INTERFACE dynamics_init

    INTERFACE dynamics_init_checks
       MODULE PROCEDURE dynamics_init_checks
    END INTERFACE dynamics_init_checks

    INTERFACE dynamics_header
       MODULE PROCEDURE dynamics_header
    END INTERFACE dynamics_header

    INTERFACE dynamics_actions
       MODULE PROCEDURE dynamics_actions
       MODULE PROCEDURE dynamics_actions_ij
    END INTERFACE dynamics_actions

    INTERFACE dynamics_non_advective_processes
       MODULE PROCEDURE dynamics_non_advective_processes
       MODULE PROCEDURE dynamics_non_advective_processes_ij
    END INTERFACE dynamics_non_advective_processes

    INTERFACE dynamics_exchange_horiz
       MODULE PROCEDURE dynamics_exchange_horiz
    END INTERFACE dynamics_exchange_horiz

    INTERFACE dynamics_prognostic_equations
       MODULE PROCEDURE dynamics_prognostic_equations
       MODULE PROCEDURE dynamics_prognostic_equations_ij
    END INTERFACE dynamics_prognostic_equations

    INTERFACE dynamics_boundary_conditions
       MODULE PROCEDURE dynamics_boundary_conditions
    END INTERFACE dynamics_boundary_conditions

    INTERFACE dynamics_swap_timelevel
       MODULE PROCEDURE dynamics_swap_timelevel
    END INTERFACE dynamics_swap_timelevel

    INTERFACE dynamics_3d_data_averaging
       MODULE PROCEDURE dynamics_3d_data_averaging
    END INTERFACE dynamics_3d_data_averaging

    INTERFACE dynamics_data_output_2d
       MODULE PROCEDURE dynamics_data_output_2d
    END INTERFACE dynamics_data_output_2d

    INTERFACE dynamics_data_output_3d
       MODULE PROCEDURE dynamics_data_output_3d
    END INTERFACE dynamics_data_output_3d

    INTERFACE dynamics_statistics
       MODULE PROCEDURE dynamics_statistics
    END INTERFACE dynamics_statistics

    INTERFACE dynamics_rrd_global
       MODULE PROCEDURE dynamics_rrd_global
    END INTERFACE dynamics_rrd_global

    INTERFACE dynamics_rrd_local
       MODULE PROCEDURE dynamics_rrd_local
    END INTERFACE dynamics_rrd_local

    INTERFACE dynamics_wrd_global
       MODULE PROCEDURE dynamics_wrd_global
    END INTERFACE dynamics_wrd_global

    INTERFACE dynamics_wrd_local
       MODULE PROCEDURE dynamics_wrd_local
    END INTERFACE dynamics_wrd_local

    INTERFACE dynamics_last_actions
       MODULE PROCEDURE dynamics_last_actions
    END INTERFACE dynamics_last_actions


 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific namelist
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_parin


    CHARACTER (LEN=80)  ::  line  !< dummy string that contains the current line of the parameter file

    NAMELIST /dynamics_parameters/  &
       dynamics_module_enabled

    line = ' '
!
!-- Try to find module-specific namelist
    REWIND ( 11 )
    line = ' '
    DO   WHILE ( INDEX( line, '&dynamics_parameters' ) == 0 )
       READ ( 11, '(A)', END=12 )  line
    ENDDO
    BACKSPACE ( 11 )

!-- Set default module switch to true
    dynamics_module_enabled = .TRUE.

!-- Read user-defined namelist
    READ ( 11, dynamics_parameters, ERR = 10 )

    GOTO 12

10  BACKSPACE( 11 )
    READ( 11 , '(A)') line
    CALL parin_fail_message( 'dynamics_parameters', line )

12  CONTINUE

 END SUBROUTINE dynamics_parin


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check control parameters and deduce further quantities.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_check_parameters


 END SUBROUTINE dynamics_check_parameters


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set module-specific timeseries units and labels
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_check_data_output_ts( dots_max, dots_num, dots_label, dots_unit )


    INTEGER(iwp),      INTENT(IN)     ::  dots_max
    INTEGER(iwp),      INTENT(INOUT)  ::  dots_num
    CHARACTER (LEN=*), DIMENSION(dots_max), INTENT(INOUT)  :: dots_label
    CHARACTER (LEN=*), DIMENSION(dots_max), INTENT(INOUT)  :: dots_unit

!
!-- Next line is to avoid compiler warning about unused variables. Please remove.
    IF ( dots_num == 0  .OR.  dots_label(1)(1:1) == ' '  .OR.  dots_unit(1)(1:1) == ' ' )  CONTINUE


 END SUBROUTINE dynamics_check_data_output_ts


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set the unit of module-specific profile output quantities. For those variables not recognized,
!> the parameter unit is set to "illegal", which tells the calling routine that the output variable
!> is not defined and leads to a program abort.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_check_data_output_pr( variable, var_count, unit, dopr_unit )


    CHARACTER (LEN=*) ::  unit     !<
    CHARACTER (LEN=*) ::  variable !<
    CHARACTER (LEN=*) ::  dopr_unit !< local value of dopr_unit

    INTEGER(iwp) ::  var_count     !<

!
!-- Next line is to avoid compiler warning about unused variables. Please remove.
    IF ( unit(1:1) == ' '  .OR.  dopr_unit(1:1) == ' '  .OR.  var_count == 0 )  CONTINUE

    SELECT CASE ( TRIM( variable ) )

!       CASE ( 'var_name' )

       CASE DEFAULT
          unit = 'illegal'

    END SELECT


 END SUBROUTINE dynamics_check_data_output_pr


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set the unit of module-specific output quantities. For those variables not recognized,
!> the parameter unit is set to "illegal", which tells the calling routine that the output variable
!< is not defined and leads to a program abort.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_check_data_output( variable, unit )


    CHARACTER (LEN=*) ::  unit     !<
    CHARACTER (LEN=*) ::  variable !<

    SELECT CASE ( TRIM( variable ) )

!       CASE ( 'u2' )

       CASE DEFAULT
          unit = 'illegal'

    END SELECT


 END SUBROUTINE dynamics_check_data_output


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Initialize module-specific masked output
!------------------------------------------------------------------------------!
 SUBROUTINE dynamics_init_masks( variable, unit )


    CHARACTER (LEN=*) ::  unit     !<
    CHARACTER (LEN=*) ::  variable !<


    SELECT CASE ( TRIM( variable ) )

!       CASE ( 'u2' )

       CASE DEFAULT
          unit = 'illegal'

    END SELECT


 END SUBROUTINE dynamics_init_masks


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize module-specific arrays
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_init_arrays


 END SUBROUTINE dynamics_init_arrays


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Execution of module-specific initializing actions
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_init


 END SUBROUTINE dynamics_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Perform module-specific post-initialization checks
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_init_checks

    INTEGER(iwp) ::  i !< loop index in x-direction
    INTEGER(iwp) ::  j !< loop index in y-direction
    INTEGER(iwp) ::  k !< loop index in z-direction

    LOGICAL      ::  realistic_q = .TRUE. !< flag indicating realistic mixing ratios

    REAL(wp)     ::  e_s !< saturation water vapor pressure
    REAL(wp)     ::  q_s !< saturation mixing ratio
    REAL(wp)     ::  t_l !< actual temperature

!
!-- Check for realistic initial mixing ratio. This must be in a realistic phyiscial range and must 
!-- not exceed the saturation mixing ratio by more than 2 percent. Please note, the check is 
!-- performed for each grid point (not just for a vertical profile), in order to cover also 
!-- three-dimensional initialization. Note, this check gives an error only for the initial run not
!-- for a restart run. In case there are no cloud physics considered, the mixing ratio can exceed
!-- the saturation moisture. This case a warning is given. 
    IF ( humidity  .AND.  .NOT. neutral )  THEN
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Calculate actual temperature, water vapor saturation pressure, and based on this
!--             the saturation mixing ratio.
                t_l = exner(k) * pt(k,j,i)
                e_s = magnus( t_l )
                q_s = rd_d_rv * e_s / ( hyp(k) - e_s )

                IF ( q(k,j,i) > 1.02_wp * q_s )  realistic_q = .FALSE. 
             ENDDO
          ENDDO
       ENDDO
!
!--    Since the check is performed locally, merge the logical flag from all mpi ranks, 
!--    in order to do not print the error message multiple times.
#if defined( __parallel )
       CALL MPI_ALLREDUCE( MPI_IN_PLACE, realistic_q, 1, MPI_LOGICAL, MPI_LAND, comm2d, ierr)
#endif

       IF ( .NOT. realistic_q  .AND.                                                               &
            TRIM( initializing_actions ) /= 'read_restart_data' )  THEN
          message_string = 'The initial mixing ratio exceeds the saturation mixing ratio.'
          CALL message( 'dynamic_init_checks', 'PA0697', 2, 2, 0, 6, 0 )
       ELSEIF ( .NOT. realistic_q  .AND.                                                           &
                TRIM( initializing_actions ) == 'read_restart_data' )  THEN
          message_string = 'The mixing ratio exceeds the saturation mixing ratio.'
          CALL message( 'dynamic_init_checks', 'PA0697', 0, 1, 0, 6, 0 )
       ENDIF
    ENDIF

 END SUBROUTINE dynamics_init_checks


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set the grids on which module-specific output quantities are defined. Allowed values for
!> grid_x are "x" and "xu", for grid_y "y" and "yv", and for grid_z "zu" and "zw".
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_define_netcdf_grid( variable, found, grid_x, grid_y, grid_z )


    CHARACTER (LEN=*) ::  grid_x     !<
    CHARACTER (LEN=*) ::  grid_y     !<
    CHARACTER (LEN=*) ::  grid_z     !<
    CHARACTER (LEN=*) ::  variable   !<

    LOGICAL ::  found   !<


    SELECT CASE ( TRIM( variable ) )

!       CASE ( 'u2' )

       CASE DEFAULT
          found  = .FALSE.
          grid_x = 'none'
          grid_y = 'none'
          grid_z = 'none'

    END SELECT


 END SUBROUTINE dynamics_define_netcdf_grid


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Print a header with module-specific information.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_header( io )


    INTEGER(iwp) ::  io   !<

!
!-- If no module-specific variables are read from the namelist-file, no information will be printed.
    IF ( .NOT. dynamics_module_enabled )  THEN
       WRITE ( io, 100 )
       RETURN
    ENDIF

!
!-- Printing the information.
    WRITE ( io, 110 )

!
!-- Format-descriptors
100 FORMAT (//' *** dynamic module disabled'/)
110 FORMAT (//1X,78('#')                                                       &
            //' User-defined variables and actions:'/                          &
              ' -----------------------------------'//)

 END SUBROUTINE dynamics_header


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Execute module-specific actions for all grid points
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_actions( location )


    CHARACTER (LEN=*) ::  location !<

!    INTEGER(iwp) ::  i !<
!    INTEGER(iwp) ::  j !<
!    INTEGER(iwp) ::  k !<

!
!-- Here the user-defined actions follow
!-- No calls for single grid points are allowed at locations before and
!-- after the timestep, since these calls are not within an i,j-loop
    SELECT CASE ( location )

       CASE ( 'before_timestep' )


       CASE ( 'before_prognostic_equations' )


       CASE ( 'after_integration' )


       CASE ( 'after_timestep' )


       CASE ( 'u-tendency' )


       CASE ( 'v-tendency' )


       CASE ( 'w-tendency' )


       CASE ( 'pt-tendency' )


       CASE ( 'sa-tendency' )


       CASE ( 'e-tendency' )


       CASE ( 'q-tendency' )


       CASE ( 's-tendency' )


       CASE DEFAULT
          CONTINUE

    END SELECT

 END SUBROUTINE dynamics_actions


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Execute module-specific actions for grid point i,j
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_actions_ij( i, j, location )


    CHARACTER (LEN=*) ::  location

    INTEGER(iwp) ::  i
    INTEGER(iwp) ::  j

!
!-- Here the user-defined actions follow
    SELECT CASE ( location )

       CASE ( 'u-tendency' )

!--       Next line is to avoid compiler warning about unused variables. Please remove.
          IF ( i +  j < 0 )  CONTINUE

       CASE ( 'v-tendency' )


       CASE ( 'w-tendency' )


       CASE ( 'pt-tendency' )


       CASE ( 'sa-tendency' )


       CASE ( 'e-tendency' )


       CASE ( 'q-tendency' )


       CASE ( 's-tendency' )


       CASE DEFAULT
          CONTINUE

    END SELECT

 END SUBROUTINE dynamics_actions_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute module-specific non-advective processes for all grid points
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_non_advective_processes



 END SUBROUTINE dynamics_non_advective_processes


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute module-specific non-advective processes for grid points i,j
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_non_advective_processes_ij( i, j )


    INTEGER(iwp) ::  i                 !<
    INTEGER(iwp) ::  j                 !<

!
!--    Next line is just to avoid compiler warnings about unused variables. You may remove it.
       IF ( i + j < 0 )  CONTINUE


 END SUBROUTINE dynamics_non_advective_processes_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Perform module-specific horizontal boundary exchange
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_exchange_horiz



 END SUBROUTINE dynamics_exchange_horiz


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute module-specific prognostic equations for all grid points
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_prognostic_equations



 END SUBROUTINE dynamics_prognostic_equations


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute module-specific prognostic equations for grid point i,j
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_prognostic_equations_ij( i, j, i_omp_start, tn )


    INTEGER(iwp), INTENT(IN) ::  i            !< grid index in x-direction
    INTEGER(iwp), INTENT(IN) ::  j            !< grid index in y-direction
    INTEGER(iwp), INTENT(IN) ::  i_omp_start  !< first loop index of i-loop in prognostic_equations
    INTEGER(iwp), INTENT(IN) ::  tn           !< task number of openmp task

!
!-- Next line is just to avoid compiler warnings about unused variables. You may remove it.
    IF ( i + j + i_omp_start + tn < 0 )  CONTINUE

 END SUBROUTINE dynamics_prognostic_equations_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute boundary conditions of dynamics model
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_boundary_conditions

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !< grid index x direction
    INTEGER(iwp) ::  j  !< grid index y direction
    INTEGER(iwp) ::  k  !< grid index z direction
    INTEGER(iwp) ::  l  !< running index boundary type, for up- and downward-facing walls
    INTEGER(iwp) ::  m  !< running index surface elements

    REAL(wp)    ::  c_max !< maximum phase velocity allowed by CFL criterion, used for outflow boundary condition
    REAL(wp)    ::  denom !< horizontal gradient of velocity component normal to the outflow boundary

!
!-- Bottom boundary
    IF ( ibc_uv_b == 1 )  THEN
       u_p(nzb,:,:) = u_p(nzb+1,:,:)
       v_p(nzb,:,:) = v_p(nzb+1,:,:)
    ENDIF
!
!-- Set zero vertical velocity at topography top (l=0), or bottom (l=1) in case
!-- of downward-facing surfaces.
    DO  l = 0, 1
       !$OMP PARALLEL DO PRIVATE( i, j, k )
       !$ACC PARALLEL LOOP PRIVATE(i, j, k) &
       !$ACC PRESENT(bc_h, w_p)
       DO  m = 1, bc_h(l)%ns
          i = bc_h(l)%i(m)
          j = bc_h(l)%j(m)
          k = bc_h(l)%k(m)
          w_p(k+bc_h(l)%koff,j,i) = 0.0_wp
       ENDDO
    ENDDO

!
!-- Top boundary. A nested domain ( ibc_uv_t = 3 ) does not require settings.
    IF ( ibc_uv_t == 0 )  THEN
        !$ACC KERNELS PRESENT(u_p, v_p, u_init, v_init)
        u_p(nzt+1,:,:) = u_init(nzt+1)
        v_p(nzt+1,:,:) = v_init(nzt+1)
        !$ACC END KERNELS
    ELSEIF ( ibc_uv_t == 1 )  THEN
        u_p(nzt+1,:,:) = u_p(nzt,:,:)
        v_p(nzt+1,:,:) = v_p(nzt,:,:)
    ENDIF

!
!-- Vertical nesting: Vertical velocity not zero at the top of the fine grid
    IF (  .NOT.  child_domain  .AND.  .NOT.  nesting_offline  .AND.            &
                 TRIM(coupling_mode) /= 'vnested_fine' )  THEN
       !$ACC KERNELS PRESENT(w_p)
       w_p(nzt:nzt+1,:,:) = 0.0_wp  !< nzt is not a prognostic level (but cf. pres)
       !$ACC END KERNELS
    ENDIF

!
!-- Temperature at bottom and top boundary.
!-- In case of coupled runs (ibc_pt_b = 2) the temperature is given by
!-- the sea surface temperature of the coupled ocean model.
!-- Dirichlet
    IF ( .NOT. neutral )  THEN
       IF ( ibc_pt_b == 0 )  THEN
          DO  l = 0, 1
             !$OMP PARALLEL DO PRIVATE( i, j, k )
             DO  m = 1, bc_h(l)%ns
                i = bc_h(l)%i(m)
                j = bc_h(l)%j(m)
                k = bc_h(l)%k(m)
                pt_p(k+bc_h(l)%koff,j,i) = pt(k+bc_h(l)%koff,j,i)
             ENDDO
          ENDDO
!
!--    Neumann, zero-gradient
       ELSEIF ( ibc_pt_b == 1 )  THEN
          DO  l = 0, 1
             !$OMP PARALLEL DO PRIVATE( i, j, k )
             !$ACC PARALLEL LOOP PRIVATE(i, j, k) &
             !$ACC PRESENT(bc_h, pt_p)
             DO  m = 1, bc_h(l)%ns
                i = bc_h(l)%i(m)
                j = bc_h(l)%j(m)
                k = bc_h(l)%k(m)
                pt_p(k+bc_h(l)%koff,j,i) = pt_p(k,j,i)
             ENDDO
          ENDDO
       ENDIF

!
!--    Temperature at top boundary
       IF ( ibc_pt_t == 0 )  THEN
           pt_p(nzt+1,:,:) = pt(nzt+1,:,:)
!
!--        In case of nudging adjust top boundary to pt which is
!--        read in from NUDGING-DATA
           IF ( nudging )  THEN
              pt_p(nzt+1,:,:) = pt_init(nzt+1)
           ENDIF
       ELSEIF ( ibc_pt_t == 1 )  THEN
           pt_p(nzt+1,:,:) = pt_p(nzt,:,:)
       ELSEIF ( ibc_pt_t == 2 )  THEN
           !$ACC KERNELS PRESENT(pt_p, dzu)
           pt_p(nzt+1,:,:) = pt_p(nzt,:,:) + bc_pt_t_val * dzu(nzt+1)
           !$ACC END KERNELS
       ENDIF
    ENDIF
!
!-- Boundary conditions for total water content,
!-- bottom and top boundary (see also temperature)
    IF ( humidity )  THEN
!
!--    Surface conditions for constant_humidity_flux
!--    Run loop over all non-natural and natural walls. Note, in wall-datatype
!--    the k coordinate belongs to the atmospheric grid point, therefore, set
!--    q_p at k-1
       IF ( ibc_q_b == 0 ) THEN

          DO  l = 0, 1
             !$OMP PARALLEL DO PRIVATE( i, j, k )
             DO  m = 1, bc_h(l)%ns
                i = bc_h(l)%i(m)
                j = bc_h(l)%j(m)
                k = bc_h(l)%k(m)
                q_p(k+bc_h(l)%koff,j,i) = q(k+bc_h(l)%koff,j,i)
             ENDDO
          ENDDO

       ELSE

          DO  l = 0, 1
             !$OMP PARALLEL DO PRIVATE( i, j, k )
             DO  m = 1, bc_h(l)%ns
                i = bc_h(l)%i(m)
                j = bc_h(l)%j(m)
                k = bc_h(l)%k(m)
                q_p(k+bc_h(l)%koff,j,i) = q_p(k,j,i)
             ENDDO
          ENDDO
       ENDIF
!
!--    Top boundary
       IF ( ibc_q_t == 0 ) THEN
          q_p(nzt+1,:,:) = q(nzt+1,:,:)
       ELSEIF ( ibc_q_t == 1 ) THEN
          q_p(nzt+1,:,:) = q_p(nzt,:,:) + bc_q_t_val * dzu(nzt+1)
       ENDIF
    ENDIF
!
!-- Boundary conditions for scalar,
!-- bottom and top boundary (see also temperature)
    IF ( passive_scalar )  THEN
!
!--    Surface conditions for constant_humidity_flux
!--    Run loop over all non-natural and natural walls. Note, in wall-datatype
!--    the k coordinate belongs to the atmospheric grid point, therefore, set
!--    s_p at k-1
       IF ( ibc_s_b == 0 ) THEN

          DO  l = 0, 1
             !$OMP PARALLEL DO PRIVATE( i, j, k )
             DO  m = 1, bc_h(l)%ns
                i = bc_h(l)%i(m)
                j = bc_h(l)%j(m)
                k = bc_h(l)%k(m)
                s_p(k+bc_h(l)%koff,j,i) = s(k+bc_h(l)%koff,j,i)
             ENDDO
          ENDDO

       ELSE

          DO  l = 0, 1
             !$OMP PARALLEL DO PRIVATE( i, j, k )
             DO  m = 1, bc_h(l)%ns
                i = bc_h(l)%i(m)
                j = bc_h(l)%j(m)
                k = bc_h(l)%k(m)
                s_p(k+bc_h(l)%koff,j,i) = s_p(k,j,i)
             ENDDO
          ENDDO
       ENDIF
!
!--    Top boundary condition
       IF ( ibc_s_t == 0 )  THEN
          s_p(nzt+1,:,:) = s(nzt+1,:,:)
       ELSEIF ( ibc_s_t == 1 )  THEN
          s_p(nzt+1,:,:) = s_p(nzt,:,:)
       ELSEIF ( ibc_s_t == 2 )  THEN
          s_p(nzt+1,:,:) = s_p(nzt,:,:) + bc_s_t_val * dzu(nzt+1)
       ENDIF

    ENDIF
!
!-- In case of inflow or nest boundary at the south boundary the boundary for v
!-- is at nys and in case of inflow or nest boundary at the left boundary the
!-- boundary for u is at nxl. Since in prognostic_equations (cache optimized
!-- version) these levels are handled as a prognostic level, boundary values
!-- have to be restored here.
    IF ( bc_dirichlet_s )  THEN
       v_p(:,nys,:) = v_p(:,nys-1,:)
    ELSEIF ( bc_dirichlet_l ) THEN
       u_p(:,:,nxl) = u_p(:,:,nxl-1)
    ENDIF

!
!-- The same restoration for u at i=nxl and v at j=nys as above must be made
!-- in case of nest boundaries. This must not be done in case of vertical nesting
!-- mode as in that case the lateral boundaries are actually cyclic.
!-- Lateral oundary conditions for TKE and dissipation are set
!-- in tcm_boundary_conds.
    IF ( nesting_mode /= 'vertical'  .OR.  nesting_offline )  THEN
       IF ( bc_dirichlet_s )  THEN
          v_p(:,nys,:) = v_p(:,nys-1,:)
       ENDIF
       IF ( bc_dirichlet_l )  THEN
          u_p(:,:,nxl) = u_p(:,:,nxl-1)
       ENDIF
    ENDIF

!
!-- Lateral boundary conditions for scalar quantities at the outflow.
!-- Lateral oundary conditions for TKE and dissipation are set
!-- in tcm_boundary_conds.
    IF ( bc_radiation_s )  THEN
       pt_p(:,nys-1,:)     = pt_p(:,nys,:)
       IF ( humidity )  THEN
          q_p(:,nys-1,:) = q_p(:,nys,:)
       ENDIF
       IF ( passive_scalar )  s_p(:,nys-1,:) = s_p(:,nys,:)
    ELSEIF ( bc_radiation_n )  THEN
       pt_p(:,nyn+1,:)     = pt_p(:,nyn,:)
       IF ( humidity )  THEN
          q_p(:,nyn+1,:) = q_p(:,nyn,:)
       ENDIF
       IF ( passive_scalar )  s_p(:,nyn+1,:) = s_p(:,nyn,:)
    ELSEIF ( bc_radiation_l )  THEN
       pt_p(:,:,nxl-1)     = pt_p(:,:,nxl)
       IF ( humidity )  THEN
          q_p(:,:,nxl-1) = q_p(:,:,nxl)
       ENDIF
       IF ( passive_scalar )  s_p(:,:,nxl-1) = s_p(:,:,nxl)
    ELSEIF ( bc_radiation_r )  THEN
       pt_p(:,:,nxr+1)     = pt_p(:,:,nxr)
       IF ( humidity )  THEN
          q_p(:,:,nxr+1) = q_p(:,:,nxr)
       ENDIF
       IF ( passive_scalar )  s_p(:,:,nxr+1) = s_p(:,:,nxr)
    ENDIF

!
!-- Radiation boundary conditions for the velocities at the respective outflow.
!-- The phase velocity is either assumed to the maximum phase velocity that
!-- ensures numerical stability (CFL-condition) or calculated after
!-- Orlanski(1976) and averaged along the outflow boundary.
    IF ( bc_radiation_s )  THEN

       IF ( use_cmax )  THEN
          u_p(:,-1,:) = u(:,0,:)
          v_p(:,0,:)  = v(:,1,:)
          w_p(:,-1,:) = w(:,0,:)
       ELSEIF ( .NOT. use_cmax )  THEN

          c_max = dy / dt_3d

          c_u_m_l = 0.0_wp
          c_v_m_l = 0.0_wp
          c_w_m_l = 0.0_wp

          c_u_m = 0.0_wp
          c_v_m = 0.0_wp
          c_w_m = 0.0_wp

!
!--       Calculate the phase speeds for u, v, and w, first local and then
!--       average along the outflow boundary.
          DO  k = nzb+1, nzt+1
             DO  i = nxl, nxr

                denom = u_m_s(k,0,i) - u_m_s(k,1,i)

                IF ( denom /= 0.0_wp )  THEN
                   c_u(k,i) = -c_max * ( u(k,0,i) - u_m_s(k,0,i) ) / ( denom * tsc(2) )
                   IF ( c_u(k,i) < 0.0_wp )  THEN
                      c_u(k,i) = 0.0_wp
                   ELSEIF ( c_u(k,i) > c_max )  THEN
                      c_u(k,i) = c_max
                   ENDIF
                ELSE
                   c_u(k,i) = c_max
                ENDIF

                denom = v_m_s(k,1,i) - v_m_s(k,2,i)

                IF ( denom /= 0.0_wp )  THEN
                   c_v(k,i) = -c_max * ( v(k,1,i) - v_m_s(k,1,i) ) / ( denom * tsc(2) )
                   IF ( c_v(k,i) < 0.0_wp )  THEN
                      c_v(k,i) = 0.0_wp
                   ELSEIF ( c_v(k,i) > c_max )  THEN
                      c_v(k,i) = c_max
                   ENDIF
                ELSE
                   c_v(k,i) = c_max
                ENDIF

                denom = w_m_s(k,0,i) - w_m_s(k,1,i)

                IF ( denom /= 0.0_wp )  THEN
                   c_w(k,i) = -c_max * ( w(k,0,i) - w_m_s(k,0,i) ) / ( denom * tsc(2) )
                   IF ( c_w(k,i) < 0.0_wp )  THEN
                      c_w(k,i) = 0.0_wp
                   ELSEIF ( c_w(k,i) > c_max )  THEN
                      c_w(k,i) = c_max
                   ENDIF
                ELSE
                   c_w(k,i) = c_max
                ENDIF

                c_u_m_l(k) = c_u_m_l(k) + c_u(k,i)
                c_v_m_l(k) = c_v_m_l(k) + c_v(k,i)
                c_w_m_l(k) = c_w_m_l(k) + c_w(k,i)

             ENDDO
          ENDDO

#if defined( __parallel )
          IF ( collective_wait )  CALL MPI_BARRIER( comm1dx, ierr )
          CALL MPI_ALLREDUCE( c_u_m_l(nzb+1), c_u_m(nzb+1), nzt-nzb, MPI_REAL, &
                              MPI_SUM, comm1dx, ierr )
          IF ( collective_wait )  CALL MPI_BARRIER( comm1dx, ierr )
          CALL MPI_ALLREDUCE( c_v_m_l(nzb+1), c_v_m(nzb+1), nzt-nzb, MPI_REAL, &
                              MPI_SUM, comm1dx, ierr )
          IF ( collective_wait )  CALL MPI_BARRIER( comm1dx, ierr )
          CALL MPI_ALLREDUCE( c_w_m_l(nzb+1), c_w_m(nzb+1), nzt-nzb, MPI_REAL, &
                              MPI_SUM, comm1dx, ierr )
#else
          c_u_m = c_u_m_l
          c_v_m = c_v_m_l
          c_w_m = c_w_m_l
#endif

          c_u_m = c_u_m / (nx+1)
          c_v_m = c_v_m / (nx+1)
          c_w_m = c_w_m / (nx+1)

!
!--       Save old timelevels for the next timestep
          IF ( intermediate_timestep_count == 1 )  THEN
             u_m_s(:,:,:) = u(:,0:1,:)
             v_m_s(:,:,:) = v(:,1:2,:)
             w_m_s(:,:,:) = w(:,0:1,:)
          ENDIF

!
!--       Calculate the new velocities
          DO  k = nzb+1, nzt+1
             DO  i = nxlg, nxrg
                u_p(k,-1,i) = u(k,-1,i) - dt_3d * tsc(2) * c_u_m(k) *          &
                                       ( u(k,-1,i) - u(k,0,i) ) * ddy

                v_p(k,0,i)  = v(k,0,i)  - dt_3d * tsc(2) * c_v_m(k) *          &
                                       ( v(k,0,i) - v(k,1,i) ) * ddy

                w_p(k,-1,i) = w(k,-1,i) - dt_3d * tsc(2) * c_w_m(k) *          &
                                       ( w(k,-1,i) - w(k,0,i) ) * ddy
             ENDDO
          ENDDO

!
!--       Bottom boundary at the outflow
          IF ( ibc_uv_b == 0 )  THEN
             u_p(nzb,-1,:) = 0.0_wp
             v_p(nzb,0,:)  = 0.0_wp
          ELSE
             u_p(nzb,-1,:) =  u_p(nzb+1,-1,:)
             v_p(nzb,0,:)  =  v_p(nzb+1,0,:)
          ENDIF
          w_p(nzb,-1,:) = 0.0_wp

!
!--       Top boundary at the outflow
          IF ( ibc_uv_t == 0 )  THEN
             u_p(nzt+1,-1,:) = u_init(nzt+1)
             v_p(nzt+1,0,:)  = v_init(nzt+1)
          ELSE
             u_p(nzt+1,-1,:) = u_p(nzt,-1,:)
             v_p(nzt+1,0,:)  = v_p(nzt,0,:)
          ENDIF
          w_p(nzt:nzt+1,-1,:) = 0.0_wp

       ENDIF

    ENDIF

    IF ( bc_radiation_n )  THEN

       IF ( use_cmax )  THEN
          u_p(:,ny+1,:) = u(:,ny,:)
          v_p(:,ny+1,:) = v(:,ny,:)
          w_p(:,ny+1,:) = w(:,ny,:)
       ELSEIF ( .NOT. use_cmax )  THEN

          c_max = dy / dt_3d

          c_u_m_l = 0.0_wp
          c_v_m_l = 0.0_wp
          c_w_m_l = 0.0_wp

          c_u_m = 0.0_wp
          c_v_m = 0.0_wp
          c_w_m = 0.0_wp

!
!--       Calculate the phase speeds for u, v, and w, first local and then
!--       average along the outflow boundary.
          DO  k = nzb+1, nzt+1
             DO  i = nxl, nxr

                denom = u_m_n(k,ny,i) - u_m_n(k,ny-1,i)

                IF ( denom /= 0.0_wp )  THEN
                   c_u(k,i) = -c_max * ( u(k,ny,i) - u_m_n(k,ny,i) ) / ( denom * tsc(2) )
                   IF ( c_u(k,i) < 0.0_wp )  THEN
                      c_u(k,i) = 0.0_wp
                   ELSEIF ( c_u(k,i) > c_max )  THEN
                      c_u(k,i) = c_max
                   ENDIF
                ELSE
                   c_u(k,i) = c_max
                ENDIF

                denom = v_m_n(k,ny,i) - v_m_n(k,ny-1,i)

                IF ( denom /= 0.0_wp )  THEN
                   c_v(k,i) = -c_max * ( v(k,ny,i) - v_m_n(k,ny,i) ) / ( denom * tsc(2) )
                   IF ( c_v(k,i) < 0.0_wp )  THEN
                      c_v(k,i) = 0.0_wp
                   ELSEIF ( c_v(k,i) > c_max )  THEN
                      c_v(k,i) = c_max
                   ENDIF
                ELSE
                   c_v(k,i) = c_max
                ENDIF

                denom = w_m_n(k,ny,i) - w_m_n(k,ny-1,i)

                IF ( denom /= 0.0_wp )  THEN
                   c_w(k,i) = -c_max * ( w(k,ny,i) - w_m_n(k,ny,i) ) / ( denom * tsc(2) )
                   IF ( c_w(k,i) < 0.0_wp )  THEN
                      c_w(k,i) = 0.0_wp
                   ELSEIF ( c_w(k,i) > c_max )  THEN
                      c_w(k,i) = c_max
                   ENDIF
                ELSE
                   c_w(k,i) = c_max
                ENDIF

                c_u_m_l(k) = c_u_m_l(k) + c_u(k,i)
                c_v_m_l(k) = c_v_m_l(k) + c_v(k,i)
                c_w_m_l(k) = c_w_m_l(k) + c_w(k,i)

             ENDDO
          ENDDO

#if defined( __parallel )
          IF ( collective_wait )  CALL MPI_BARRIER( comm1dx, ierr )
          CALL MPI_ALLREDUCE( c_u_m_l(nzb+1), c_u_m(nzb+1), nzt-nzb, MPI_REAL, &
                              MPI_SUM, comm1dx, ierr )
          IF ( collective_wait )  CALL MPI_BARRIER( comm1dx, ierr )
          CALL MPI_ALLREDUCE( c_v_m_l(nzb+1), c_v_m(nzb+1), nzt-nzb, MPI_REAL, &
                              MPI_SUM, comm1dx, ierr )
          IF ( collective_wait )  CALL MPI_BARRIER( comm1dx, ierr )
          CALL MPI_ALLREDUCE( c_w_m_l(nzb+1), c_w_m(nzb+1), nzt-nzb, MPI_REAL, &
                              MPI_SUM, comm1dx, ierr )
#else
          c_u_m = c_u_m_l
          c_v_m = c_v_m_l
          c_w_m = c_w_m_l
#endif

          c_u_m = c_u_m / (nx+1)
          c_v_m = c_v_m / (nx+1)
          c_w_m = c_w_m / (nx+1)

!
!--       Save old timelevels for the next timestep
          IF ( intermediate_timestep_count == 1 )  THEN
                u_m_n(:,:,:) = u(:,ny-1:ny,:)
                v_m_n(:,:,:) = v(:,ny-1:ny,:)
                w_m_n(:,:,:) = w(:,ny-1:ny,:)
          ENDIF

!
!--       Calculate the new velocities
          DO  k = nzb+1, nzt+1
             DO  i = nxlg, nxrg
                u_p(k,ny+1,i) = u(k,ny+1,i) - dt_3d * tsc(2) * c_u_m(k) *      &
                                       ( u(k,ny+1,i) - u(k,ny,i) ) * ddy

                v_p(k,ny+1,i) = v(k,ny+1,i)  - dt_3d * tsc(2) * c_v_m(k) *     &
                                       ( v(k,ny+1,i) - v(k,ny,i) ) * ddy

                w_p(k,ny+1,i) = w(k,ny+1,i) - dt_3d * tsc(2) * c_w_m(k) *      &
                                       ( w(k,ny+1,i) - w(k,ny,i) ) * ddy
             ENDDO
          ENDDO

!
!--       Bottom boundary at the outflow
          IF ( ibc_uv_b == 0 )  THEN
             u_p(nzb,ny+1,:) = 0.0_wp
             v_p(nzb,ny+1,:) = 0.0_wp
          ELSE
             u_p(nzb,ny+1,:) =  u_p(nzb+1,ny+1,:)
             v_p(nzb,ny+1,:) =  v_p(nzb+1,ny+1,:)
          ENDIF
          w_p(nzb,ny+1,:) = 0.0_wp

!
!--       Top boundary at the outflow
          IF ( ibc_uv_t == 0 )  THEN
             u_p(nzt+1,ny+1,:) = u_init(nzt+1)
             v_p(nzt+1,ny+1,:) = v_init(nzt+1)
          ELSE
             u_p(nzt+1,ny+1,:) = u_p(nzt,nyn+1,:)
             v_p(nzt+1,ny+1,:) = v_p(nzt,nyn+1,:)
          ENDIF
          w_p(nzt:nzt+1,ny+1,:) = 0.0_wp

       ENDIF

    ENDIF

    IF ( bc_radiation_l )  THEN

       IF ( use_cmax )  THEN
          u_p(:,:,0)  = u(:,:,1)
          v_p(:,:,-1) = v(:,:,0)
          w_p(:,:,-1) = w(:,:,0)
       ELSEIF ( .NOT. use_cmax )  THEN

          c_max = dx / dt_3d

          c_u_m_l = 0.0_wp
          c_v_m_l = 0.0_wp
          c_w_m_l = 0.0_wp

          c_u_m = 0.0_wp
          c_v_m = 0.0_wp
          c_w_m = 0.0_wp

!
!--       Calculate the phase speeds for u, v, and w, first local and then
!--       average along the outflow boundary.
          DO  k = nzb+1, nzt+1
             DO  j = nys, nyn

                denom = u_m_l(k,j,1) - u_m_l(k,j,2)

                IF ( denom /= 0.0_wp )  THEN
                   c_u(k,j) = -c_max * ( u(k,j,1) - u_m_l(k,j,1) ) / ( denom * tsc(2) )
                   IF ( c_u(k,j) < 0.0_wp )  THEN
                      c_u(k,j) = 0.0_wp
                   ELSEIF ( c_u(k,j) > c_max )  THEN
                      c_u(k,j) = c_max
                   ENDIF
                ELSE
                   c_u(k,j) = c_max
                ENDIF

                denom = v_m_l(k,j,0) - v_m_l(k,j,1)

                IF ( denom /= 0.0_wp )  THEN
                   c_v(k,j) = -c_max * ( v(k,j,0) - v_m_l(k,j,0) ) / ( denom * tsc(2) )
                   IF ( c_v(k,j) < 0.0_wp )  THEN
                      c_v(k,j) = 0.0_wp
                   ELSEIF ( c_v(k,j) > c_max )  THEN
                      c_v(k,j) = c_max
                   ENDIF
                ELSE
                   c_v(k,j) = c_max
                ENDIF

                denom = w_m_l(k,j,0) - w_m_l(k,j,1)

                IF ( denom /= 0.0_wp )  THEN
                   c_w(k,j) = -c_max * ( w(k,j,0) - w_m_l(k,j,0) ) / ( denom * tsc(2) )
                   IF ( c_w(k,j) < 0.0_wp )  THEN
                      c_w(k,j) = 0.0_wp
                   ELSEIF ( c_w(k,j) > c_max )  THEN
                      c_w(k,j) = c_max
                   ENDIF
                ELSE
                   c_w(k,j) = c_max
                ENDIF

                c_u_m_l(k) = c_u_m_l(k) + c_u(k,j)
                c_v_m_l(k) = c_v_m_l(k) + c_v(k,j)
                c_w_m_l(k) = c_w_m_l(k) + c_w(k,j)

             ENDDO
          ENDDO

#if defined( __parallel )
          IF ( collective_wait )  CALL MPI_BARRIER( comm1dy, ierr )
          CALL MPI_ALLREDUCE( c_u_m_l(nzb+1), c_u_m(nzb+1), nzt-nzb, MPI_REAL, &
                              MPI_SUM, comm1dy, ierr )
          IF ( collective_wait )  CALL MPI_BARRIER( comm1dy, ierr )
          CALL MPI_ALLREDUCE( c_v_m_l(nzb+1), c_v_m(nzb+1), nzt-nzb, MPI_REAL, &
                              MPI_SUM, comm1dy, ierr )
          IF ( collective_wait )  CALL MPI_BARRIER( comm1dy, ierr )
          CALL MPI_ALLREDUCE( c_w_m_l(nzb+1), c_w_m(nzb+1), nzt-nzb, MPI_REAL, &
                              MPI_SUM, comm1dy, ierr )
#else
          c_u_m = c_u_m_l
          c_v_m = c_v_m_l
          c_w_m = c_w_m_l
#endif

          c_u_m = c_u_m / (ny+1)
          c_v_m = c_v_m / (ny+1)
          c_w_m = c_w_m / (ny+1)

!
!--       Save old timelevels for the next timestep
          IF ( intermediate_timestep_count == 1 )  THEN
                u_m_l(:,:,:) = u(:,:,1:2)
                v_m_l(:,:,:) = v(:,:,0:1)
                w_m_l(:,:,:) = w(:,:,0:1)
          ENDIF

!
!--       Calculate the new velocities
          DO  k = nzb+1, nzt+1
             DO  j = nysg, nyng
                u_p(k,j,0) = u(k,j,0) - dt_3d * tsc(2) * c_u_m(k) *            &
                                       ( u(k,j,0) - u(k,j,1) ) * ddx

                v_p(k,j,-1) = v(k,j,-1) - dt_3d * tsc(2) * c_v_m(k) *          &
                                       ( v(k,j,-1) - v(k,j,0) ) * ddx

                w_p(k,j,-1) = w(k,j,-1) - dt_3d * tsc(2) * c_w_m(k) *          &
                                       ( w(k,j,-1) - w(k,j,0) ) * ddx
             ENDDO
          ENDDO

!
!--       Bottom boundary at the outflow
          IF ( ibc_uv_b == 0 )  THEN
             u_p(nzb,:,0)  = 0.0_wp
             v_p(nzb,:,-1) = 0.0_wp
          ELSE
             u_p(nzb,:,0)  =  u_p(nzb+1,:,0)
             v_p(nzb,:,-1) =  v_p(nzb+1,:,-1)
          ENDIF
          w_p(nzb,:,-1) = 0.0_wp

!
!--       Top boundary at the outflow
          IF ( ibc_uv_t == 0 )  THEN
             u_p(nzt+1,:,0)  = u_init(nzt+1)
             v_p(nzt+1,:,-1) = v_init(nzt+1)
          ELSE
             u_p(nzt+1,:,0)  = u_p(nzt,:,0)
             v_p(nzt+1,:,-1) = v_p(nzt,:,-1)
          ENDIF
          w_p(nzt:nzt+1,:,-1) = 0.0_wp

       ENDIF

    ENDIF

    IF ( bc_radiation_r )  THEN

       IF ( use_cmax )  THEN
          u_p(:,:,nx+1) = u(:,:,nx)
          v_p(:,:,nx+1) = v(:,:,nx)
          w_p(:,:,nx+1) = w(:,:,nx)
       ELSEIF ( .NOT. use_cmax )  THEN

          c_max = dx / dt_3d

          c_u_m_l = 0.0_wp
          c_v_m_l = 0.0_wp
          c_w_m_l = 0.0_wp

          c_u_m = 0.0_wp
          c_v_m = 0.0_wp
          c_w_m = 0.0_wp

!
!--       Calculate the phase speeds for u, v, and w, first local and then
!--       average along the outflow boundary.
          DO  k = nzb+1, nzt+1
             DO  j = nys, nyn

                denom = u_m_r(k,j,nx) - u_m_r(k,j,nx-1)

                IF ( denom /= 0.0_wp )  THEN
                   c_u(k,j) = -c_max * ( u(k,j,nx) - u_m_r(k,j,nx) ) / ( denom * tsc(2) )
                   IF ( c_u(k,j) < 0.0_wp )  THEN
                      c_u(k,j) = 0.0_wp
                   ELSEIF ( c_u(k,j) > c_max )  THEN
                      c_u(k,j) = c_max
                   ENDIF
                ELSE
                   c_u(k,j) = c_max
                ENDIF

                denom = v_m_r(k,j,nx) - v_m_r(k,j,nx-1)

                IF ( denom /= 0.0_wp )  THEN
                   c_v(k,j) = -c_max * ( v(k,j,nx) - v_m_r(k,j,nx) ) / ( denom * tsc(2) )
                   IF ( c_v(k,j) < 0.0_wp )  THEN
                      c_v(k,j) = 0.0_wp
                   ELSEIF ( c_v(k,j) > c_max )  THEN
                      c_v(k,j) = c_max
                   ENDIF
                ELSE
                   c_v(k,j) = c_max
                ENDIF

                denom = w_m_r(k,j,nx) - w_m_r(k,j,nx-1)

                IF ( denom /= 0.0_wp )  THEN
                   c_w(k,j) = -c_max * ( w(k,j,nx) - w_m_r(k,j,nx) ) / ( denom * tsc(2) )
                   IF ( c_w(k,j) < 0.0_wp )  THEN
                      c_w(k,j) = 0.0_wp
                   ELSEIF ( c_w(k,j) > c_max )  THEN
                      c_w(k,j) = c_max
                   ENDIF
                ELSE
                   c_w(k,j) = c_max
                ENDIF

                c_u_m_l(k) = c_u_m_l(k) + c_u(k,j)
                c_v_m_l(k) = c_v_m_l(k) + c_v(k,j)
                c_w_m_l(k) = c_w_m_l(k) + c_w(k,j)

             ENDDO
          ENDDO

#if defined( __parallel )
          IF ( collective_wait )  CALL MPI_BARRIER( comm1dy, ierr )
          CALL MPI_ALLREDUCE( c_u_m_l(nzb+1), c_u_m(nzb+1), nzt-nzb, MPI_REAL, &
                              MPI_SUM, comm1dy, ierr )
          IF ( collective_wait )  CALL MPI_BARRIER( comm1dy, ierr )
          CALL MPI_ALLREDUCE( c_v_m_l(nzb+1), c_v_m(nzb+1), nzt-nzb, MPI_REAL, &
                              MPI_SUM, comm1dy, ierr )
          IF ( collective_wait )  CALL MPI_BARRIER( comm1dy, ierr )
          CALL MPI_ALLREDUCE( c_w_m_l(nzb+1), c_w_m(nzb+1), nzt-nzb, MPI_REAL, &
                              MPI_SUM, comm1dy, ierr )
#else
          c_u_m = c_u_m_l
          c_v_m = c_v_m_l
          c_w_m = c_w_m_l
#endif

          c_u_m = c_u_m / (ny+1)
          c_v_m = c_v_m / (ny+1)
          c_w_m = c_w_m / (ny+1)

!
!--       Save old timelevels for the next timestep
          IF ( intermediate_timestep_count == 1 )  THEN
                u_m_r(:,:,:) = u(:,:,nx-1:nx)
                v_m_r(:,:,:) = v(:,:,nx-1:nx)
                w_m_r(:,:,:) = w(:,:,nx-1:nx)
          ENDIF

!
!--       Calculate the new velocities
          DO  k = nzb+1, nzt+1
             DO  j = nysg, nyng
                u_p(k,j,nx+1) = u(k,j,nx+1) - dt_3d * tsc(2) * c_u_m(k) *      &
                                       ( u(k,j,nx+1) - u(k,j,nx) ) * ddx

                v_p(k,j,nx+1) = v(k,j,nx+1) - dt_3d * tsc(2) * c_v_m(k) *      &
                                       ( v(k,j,nx+1) - v(k,j,nx) ) * ddx

                w_p(k,j,nx+1) = w(k,j,nx+1) - dt_3d * tsc(2) * c_w_m(k) *      &
                                       ( w(k,j,nx+1) - w(k,j,nx) ) * ddx
             ENDDO
          ENDDO

!
!--       Bottom boundary at the outflow
          IF ( ibc_uv_b == 0 )  THEN
             u_p(nzb,:,nx+1) = 0.0_wp
             v_p(nzb,:,nx+1) = 0.0_wp
          ELSE
             u_p(nzb,:,nx+1) =  u_p(nzb+1,:,nx+1)
             v_p(nzb,:,nx+1) =  v_p(nzb+1,:,nx+1)
          ENDIF
          w_p(nzb,:,nx+1) = 0.0_wp

!
!--       Top boundary at the outflow
          IF ( ibc_uv_t == 0 )  THEN
             u_p(nzt+1,:,nx+1) = u_init(nzt+1)
             v_p(nzt+1,:,nx+1) = v_init(nzt+1)
          ELSE
             u_p(nzt+1,:,nx+1) = u_p(nzt,:,nx+1)
             v_p(nzt+1,:,nx+1) = v_p(nzt,:,nx+1)
          ENDIF
          w_p(nzt:nzt+1,:,nx+1) = 0.0_wp

       ENDIF

    ENDIF

 END SUBROUTINE dynamics_boundary_conditions
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Swap timelevels of module-specific array pointers
!------------------------------------------------------------------------------!
 SUBROUTINE dynamics_swap_timelevel ( mod_count )


    INTEGER, INTENT(IN) :: mod_count


    SELECT CASE ( mod_count )

       CASE ( 0 )

          u  => u_1;   u_p  => u_2
          v  => v_1;   v_p  => v_2
          w  => w_1;   w_p  => w_2
          IF ( .NOT. neutral )  THEN
             pt => pt_1;  pt_p => pt_2
          ENDIF
          IF ( humidity )  THEN
             q => q_1;    q_p => q_2
          ENDIF
          IF ( passive_scalar )  THEN
             s => s_1;    s_p => s_2
          ENDIF

       CASE ( 1 )

          u  => u_2;   u_p  => u_1
          v  => v_2;   v_p  => v_1
          w  => w_2;   w_p  => w_1
          IF ( .NOT. neutral )  THEN
             pt => pt_2;  pt_p => pt_1
          ENDIF
          IF ( humidity )  THEN
             q => q_2;    q_p => q_1
          ENDIF
          IF ( passive_scalar )  THEN
             s => s_2;    s_p => s_1
          ENDIF

    END SELECT

 END SUBROUTINE dynamics_swap_timelevel


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Sum up and time-average module-specific output quantities
!> as well as allocate the array necessary for storing the average.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_3d_data_averaging( mode, variable )


    CHARACTER (LEN=*) ::  mode    !<
    CHARACTER (LEN=*) :: variable !<


    IF ( mode == 'allocate' )  THEN

       SELECT CASE ( TRIM( variable ) )

!          CASE ( 'u2' )

          CASE DEFAULT
             CONTINUE

       END SELECT

    ELSEIF ( mode == 'sum' )  THEN

       SELECT CASE ( TRIM( variable ) )

!          CASE ( 'u2' )

          CASE DEFAULT
             CONTINUE

       END SELECT

    ELSEIF ( mode == 'average' )  THEN

       SELECT CASE ( TRIM( variable ) )

!          CASE ( 'u2' )

       END SELECT

    ENDIF


 END SUBROUTINE dynamics_3d_data_averaging


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Resorts the module-specific output quantity with indices (k,j,i) to a
!> temporary array with indices (i,j,k) and sets the grid on which it is defined.
!> Allowed values for grid are "zu" and "zw".
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_data_output_2d( av, variable, found, grid, mode, local_pf, &
                                     two_d, nzb_do, nzt_do, fill_value )


    CHARACTER (LEN=*) ::  grid     !<
    CHARACTER (LEN=*), INTENT(IN) ::  mode       !< either 'xy', 'xz' or 'yz'
    CHARACTER (LEN=*) ::  variable !<

    INTEGER(iwp) ::  av     !< flag to control data output of instantaneous or time-averaged data
!    INTEGER(iwp) ::  i      !< grid index along x-direction
!    INTEGER(iwp) ::  j      !< grid index along y-direction
!    INTEGER(iwp) ::  k      !< grid index along z-direction
!    INTEGER(iwp) ::  m      !< running index surface elements
    INTEGER(iwp) ::  nzb_do !< lower limit of the domain (usually nzb)
    INTEGER(iwp) ::  nzt_do !< upper limit of the domain (usually nzt+1)

    LOGICAL      ::  found !<
    LOGICAL      ::  two_d !< flag parameter that indicates 2D variables (horizontal cross sections)

    REAL(wp), INTENT(IN) ::  fill_value

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf !<

!
!-- Next line is just to avoid compiler warnings about unused variables. You may remove it.
    IF ( two_d .AND. av + LEN( mode ) + local_pf(nxl,nys,nzb_do) + fill_value < 0.0 )  CONTINUE

    found = .TRUE.

    SELECT CASE ( TRIM( variable ) )

!       CASE ( 'u2_xy', 'u2_xz', 'u2_yz' )

       CASE DEFAULT
          found = .FALSE.
          grid  = 'none'

    END SELECT


 END SUBROUTINE dynamics_data_output_2d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Resorts the module-specific output quantity with indices (k,j,i)
!> to a temporary array with indices (i,j,k).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_data_output_3d( av, variable, found, local_pf, fill_value, nzb_do, nzt_do )


    CHARACTER (LEN=*) ::  variable !<

    INTEGER(iwp) ::  av    !<
!    INTEGER(iwp) ::  i     !<
!    INTEGER(iwp) ::  j     !<
!    INTEGER(iwp) ::  k     !<
    INTEGER(iwp) ::  nzb_do !< lower limit of the data output (usually 0)
    INTEGER(iwp) ::  nzt_do !< vertical upper limit of the data output (usually nz_do3d)

    LOGICAL      ::  found !<

    REAL(wp), INTENT(IN) ::  fill_value    !< value for the _FillValue attribute

    REAL(sp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf !<

!
!-- Next line is to avoid compiler warning about unused variables. Please remove.
    IF ( av + local_pf(nxl,nys,nzb_do) + fill_value < 0.0 )  CONTINUE


    found = .TRUE.

    SELECT CASE ( TRIM( variable ) )

!       CASE ( 'u2' )

       CASE DEFAULT
          found = .FALSE.

    END SELECT


 END SUBROUTINE dynamics_data_output_3d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculation of module-specific statistics, i.e. horizontally averaged profiles and time series.
!> This is called for every statistic region sr, but at least for the region "total domain" (sr=0).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_statistics( mode, sr, tn )


    CHARACTER (LEN=*) ::  mode   !<
!    INTEGER(iwp) ::  i    !<
!    INTEGER(iwp) ::  j    !<
!    INTEGER(iwp) ::  k    !<
    INTEGER(iwp) ::  sr   !<
    INTEGER(iwp) ::  tn   !<

!
!-- Next line is to avoid compiler warning about unused variables. Please remove.
    IF ( sr == 0  .OR.  tn == 0 )  CONTINUE

    IF ( mode == 'profiles' )  THEN

    ELSEIF ( mode == 'time_series' )  THEN

    ENDIF

 END SUBROUTINE dynamics_statistics


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific global restart data.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_rrd_global( found )


    LOGICAL, INTENT(OUT)  ::  found


    found = .TRUE.


    SELECT CASE ( restart_string(1:length) )

       CASE ( 'global_paramter' )
!          READ ( 13 )  global_parameter

       CASE DEFAULT

          found = .FALSE.

    END SELECT


 END SUBROUTINE dynamics_rrd_global


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific processor specific restart data from file(s).
!> Subdomain index limits on file are given by nxl_on_file, etc.
!> Indices nxlc, etc. indicate the range of gridpoints to be mapped from the subdomain on file (f)
!> to the subdomain of the current PE (c). They have been calculated in routine rrd_local.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_rrd_local( k, nxlf, nxlc, nxl_on_file, nxrf, nxrc, nxr_on_file, nynf, nync,   &
                                nyn_on_file, nysf, nysc, nys_on_file, tmp_2d, tmp_3d, found )


    INTEGER(iwp) ::  k               !<
    INTEGER(iwp) ::  nxlc            !<
    INTEGER(iwp) ::  nxlf            !<
    INTEGER(iwp) ::  nxl_on_file     !<
    INTEGER(iwp) ::  nxrc            !<
    INTEGER(iwp) ::  nxrf            !<
    INTEGER(iwp) ::  nxr_on_file     !<
    INTEGER(iwp) ::  nync            !<
    INTEGER(iwp) ::  nynf            !<
    INTEGER(iwp) ::  nyn_on_file     !<
    INTEGER(iwp) ::  nysc            !<
    INTEGER(iwp) ::  nysf            !<
    INTEGER(iwp) ::  nys_on_file     !<

    LOGICAL, INTENT(OUT)  ::  found

    REAL(wp), DIMENSION(nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) :: tmp_2d   !<
    REAL(wp), DIMENSION(nzb:nzt+1,nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) :: tmp_3d   !<

!
!-- Next line is to avoid compiler warning about unused variables. Please remove.
    IF ( k + nxlc + nxlf + nxrc + nxrf + nync + nynf + nysc + nysf +           &
         tmp_2d(nys_on_file,nxl_on_file) +                                     &
         tmp_3d(nzb,nys_on_file,nxl_on_file) < 0.0 )  CONTINUE
!
!-- Here the reading of user-defined restart data follows:
!-- Sample for user-defined output

    found = .TRUE.

    SELECT CASE ( restart_string(1:length) )

!       CASE ( 'u2_av' )

       CASE DEFAULT

          found = .FALSE.

    END SELECT

 END SUBROUTINE dynamics_rrd_local


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Writes global module-specific restart data into binary file(s) for restart runs.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_wrd_global


 END SUBROUTINE dynamics_wrd_global


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Writes processor specific and module-specific restart data into binary file(s) for restart runs.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_wrd_local


 END SUBROUTINE dynamics_wrd_local


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Execute module-specific actions at the very end of the program.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_last_actions


 END SUBROUTINE dynamics_last_actions

 END MODULE dynamics_mod
