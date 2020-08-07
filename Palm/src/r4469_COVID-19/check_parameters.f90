!> @file check_parameters.f90
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
! $Id: check_parameters.f90 4444 2020-03-05 15:59:50Z raasch $
! bugfix: cpp-directives for serial mode added
! 
! 4392 2020-01-31 16:14:57Z pavelkrc
! Some error numbers revised to prevent double usage
! 
! 11:55:33Z oliver.maas
! Checks for closed channel flow implemented
! 
! 11:55:33Z oliver.maas
! Move 2-m potential temperature output to diagnostic_output_quantities
! 
! 11:55:33Z oliver.maas
! removed message PA0421, concerning old parameter recycling_yshift
! 
! 11:55:33Z oliver.maas
! adjust message to the modified parameter recycling_yshift
! 
! 11:55:33Z oliver.maas
! Check if a cross section is specified if any output cross-section quantity
! is given
! 
! 11:55:33Z oliver.maas
! Overwrite rotation_angle from namelist by value from static driver
! 
! 11:55:33Z oliver.maas
! removed conversion from recycle_absolute_quantities to raq, added check and 
! error message for correct input of recycling_method_for_thermodynamic_quantities
! 
! 11:55:33Z oliver.maas
! Corrected "Former revisions" section
! 
! 11:55:33Z oliver.maas
! bugfix error message: replaced PA184 by PA0184
! 
! 11:55:33Z oliver.maas
! added conversion from recycle_absolute_quantities to raq for recycling of 
! absolute quantities and added error message PA184 for not implemented quantities
! 
! 4142 2019-08-05 12:38:31Z suehring
! Consider spinup in number of output timesteps for averaged 2D output (merge
! from branch resler).
! 
! 4069 2019-07-01 14:05:51Z Giersch
! Masked output running index mid has been introduced as a local variable to 
! avoid runtime error (Loop variable has been modified) in time_integration
! 
! 4048 2019-06-21 21:00:21Z knoop
! Moved tcm_check_data_output to module_interface
! 
! 4039 2019-06-18 10:32:41Z suehring
! Modularize diagnostic output
! 
! 4017 2019-06-06 12:16:46Z schwenkel
! output of turbulence intensity added
! 
! 3933 2019-04-25 12:33:20Z kanani
! Alphabetical resorting in CASE, condense settings for theta_2m* into one IF clause
! 
! 3885 2019-04-11 11:29:34Z kanani
! Changes related to global restructuring of location messages and introduction 
! of additional debug messages
! 
! 3766 2019-02-26 16:23:41Z raasch
! trim added to avoid truncation compiler warnings
! 
! 3761 2019-02-25 15:31:42Z raasch
! unused variables removed
! 
! 3735 2019-02-12 09:52:40Z dom_dwd_user
! Passing variable j (averaged output?) to 
! module_interface.f90:chem_check_data_output.
! 
! 3705 2019-01-29 19:56:39Z suehring
! bugfix: renamed thetav_t to vtheta_t
! 
! 3702 2019-01-28 13:19:30Z gronemeier
! most_method removed
! 
! 3655 2019-01-07 16:51:22Z knoop
! Formatting
!
! Revision 1.1  1997/08/26 06:29:23  raasch
! Initial revision
!
!
! Description:
! ------------
!> Check control parameters and deduce further quantities.
!------------------------------------------------------------------------------!
 SUBROUTINE check_parameters


    USE arrays_3d

    USE basic_constants_and_equations_mod

    USE bulk_cloud_model_mod,                                                  &
        ONLY:  bulk_cloud_model

    USE chem_modules

    USE chemistry_model_mod,                                                   &
        ONLY:  chem_boundary_conds

    USE control_parameters

    USE grid_variables

    USE kinds

    USE indices

    USE model_1d_mod,                                                          &
        ONLY:  damp_level_1d, damp_level_ind_1d

    USE module_interface,                                                      &
        ONLY:  module_interface_check_parameters,                              &
               module_interface_check_data_output_ts,                          &
               module_interface_check_data_output_pr,                          &
               module_interface_check_data_output

    USE netcdf_data_input_mod,                                                 &
        ONLY:  init_model, input_pids_static, netcdf_data_input_check_dynamic, &
               netcdf_data_input_check_static

    USE netcdf_interface,                                                      &
        ONLY:  dopr_unit, do2d_unit, do3d_unit, netcdf_data_format,            &
               netcdf_data_format_string, dots_unit, heatflux_output_unit,     &
               waterflux_output_unit, momentumflux_output_unit,                &
               dots_max, dots_num, dots_label

    USE particle_attributes,                                                   &
        ONLY:  particle_advection, use_sgs_for_particles
        
    USE pegrid

    USE pmc_interface,                                                         &
        ONLY:  cpl_id, nested_run

    USE profil_parameter

    USE statistics

    USE subsidence_mod

    USE transpose_indices

#if defined( __parallel )
    USE vertical_nesting_mod,                                                  &
        ONLY:  vnested,                                                        &
               vnest_check_parameters
#endif


    IMPLICIT NONE

    CHARACTER (LEN=varnamelength)  ::  var           !< variable name
    CHARACTER (LEN=7)   ::  unit                     !< unit of variable
    CHARACTER (LEN=8)   ::  date                     !< current date string
    CHARACTER (LEN=10)  ::  time                     !< current time string
    CHARACTER (LEN=20)  ::  ensemble_string          !< string containing number of ensemble member
    CHARACTER (LEN=15)  ::  nest_string              !< string containing id of nested domain
    CHARACTER (LEN=40)  ::  coupling_string          !< string containing type of coupling
    CHARACTER (LEN=100) ::  action                   !< flag string

    INTEGER(iwp) ::  i                               !< loop index
    INTEGER(iwp) ::  ilen                            !< string length
    INTEGER(iwp) ::  j                               !< loop index
    INTEGER(iwp) ::  k                               !< loop index
    INTEGER(iwp) ::  kk                              !< loop index
    INTEGER(iwp) ::  mid                             !< masked output running index
    INTEGER(iwp) ::  netcdf_data_format_save         !< initial value of netcdf_data_format
    INTEGER(iwp) ::  position                        !< index position of string

    LOGICAL     ::  found                            !< flag, true if output variable is already marked for averaging

    REAL(wp)    ::  gradient                         !< local gradient
#if defined( __parallel )
    REAL(wp)    ::  dt_spinup_max                    !< maximum spinup timestep in nested domains
    REAL(wp)    ::  remote = 0.0_wp                  !< MPI id of remote processor
    REAL(wp)    ::  spinup_time_max                  !< maximum spinup time in nested domains
    REAL(wp)    ::  time_to_be_simulated_from_reference_point  !< time to be simulated from reference point
#endif


    CALL location_message( 'checking parameters', 'start' )
!
!-- At first, check static and dynamic input for consistency
    CALL netcdf_data_input_check_dynamic
    CALL netcdf_data_input_check_static
!
!-- Check for overlap combinations, which are not realized yet
    IF ( transpose_compute_overlap  .AND. numprocs == 1 )  THEN
          message_string = 'transpose-compute-overlap not implemented for single PE runs'
          CALL message( 'check_parameters', 'PA0000', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Check the coupling mode
    IF ( coupling_mode /= 'uncoupled'            .AND.                         &
         coupling_mode /= 'precursor_atmos'      .AND.                         &
         coupling_mode /= 'precursor_ocean'      .AND.                         &
         coupling_mode /= 'vnested_crse'         .AND.                         &
         coupling_mode /= 'vnested_fine'         .AND.                         &
         coupling_mode /= 'atmosphere_to_ocean'  .AND.                         &
         coupling_mode /= 'ocean_to_atmosphere' )  THEN
       message_string = 'illegal coupling mode: ' // TRIM( coupling_mode )
       CALL message( 'check_parameters', 'PA0002', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Check if humidity is set to TRUE in case of the atmospheric run (for coupled runs)
    IF ( coupling_mode == 'atmosphere_to_ocean' .AND. .NOT. humidity) THEN
       message_string = ' Humidity has to be set to .T. in the _p3d file ' //  &
                        'for coupled runs between ocean and atmosphere.'
       CALL message( 'check_parameters', 'PA0476', 1, 2, 0, 6, 0 )
    ENDIF
   
!
!-- Check dt_coupling, restart_time, dt_restart, end_time, dx, dy, nx and ny
    IF ( coupling_mode /= 'uncoupled'       .AND.                              &
         coupling_mode(1:8) /= 'vnested_'   .AND.                              &
         coupling_mode /= 'precursor_atmos' .AND.                              &
         coupling_mode /= 'precursor_ocean' )  THEN

       IF ( dt_coupling == 9999999.9_wp )  THEN
          message_string = 'dt_coupling is not set but required for coup' //   &
                           'ling mode "' //  TRIM( coupling_mode ) // '"'
          CALL message( 'check_parameters', 'PA0003', 1, 2, 0, 6, 0 )
       ENDIF

#if defined( __parallel )


       IF ( myid == 0 ) THEN
          CALL MPI_SEND( dt_coupling, 1, MPI_REAL, target_id, 11, comm_inter,  &
                         ierr )
          CALL MPI_RECV( remote, 1, MPI_REAL, target_id, 11, comm_inter,       &
                         status, ierr )
       ENDIF
       CALL MPI_BCAST( remote, 1, MPI_REAL, 0, comm2d, ierr)

       IF ( dt_coupling /= remote )  THEN
          WRITE( message_string, * ) 'coupling mode "', TRIM( coupling_mode ), &
                 '": dt_coupling = ', dt_coupling, '& is not equal to ',       &
                 'dt_coupling_remote = ', remote
          CALL message( 'check_parameters', 'PA0004', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( dt_coupling <= 0.0_wp )  THEN

          IF ( myid == 0  ) THEN
             CALL MPI_SEND( dt_max, 1, MPI_REAL, target_id, 19, comm_inter, ierr )
             CALL MPI_RECV( remote, 1, MPI_REAL, target_id, 19, comm_inter,    &
                            status, ierr )
          ENDIF
          CALL MPI_BCAST( remote, 1, MPI_REAL, 0, comm2d, ierr)

          dt_coupling = MAX( dt_max, remote )
          WRITE( message_string, * ) 'coupling mode "', TRIM( coupling_mode ), &
                 '": dt_coupling <= 0.0 & is not allowed and is reset to ',    &
                 'MAX(dt_max(A,O)) = ', dt_coupling
          CALL message( 'check_parameters', 'PA0005', 0, 1, 0, 6, 0 )
       ENDIF

       IF ( myid == 0 ) THEN
          CALL MPI_SEND( restart_time, 1, MPI_REAL, target_id, 12, comm_inter, &
                         ierr )
          CALL MPI_RECV( remote, 1, MPI_REAL, target_id, 12, comm_inter,       &
                         status, ierr )
       ENDIF
       CALL MPI_BCAST( remote, 1, MPI_REAL, 0, comm2d, ierr)

       IF ( restart_time /= remote )  THEN
          WRITE( message_string, * ) 'coupling mode "', TRIM( coupling_mode ), &
                 '": restart_time = ', restart_time, '& is not equal to ',     &
                 'restart_time_remote = ', remote
          CALL message( 'check_parameters', 'PA0006', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( myid == 0 ) THEN
          CALL MPI_SEND( dt_restart, 1, MPI_REAL, target_id, 13, comm_inter,   &
                         ierr )
          CALL MPI_RECV( remote, 1, MPI_REAL, target_id, 13, comm_inter,       &
                         status, ierr )
       ENDIF
       CALL MPI_BCAST( remote, 1, MPI_REAL, 0, comm2d, ierr)

       IF ( dt_restart /= remote )  THEN
          WRITE( message_string, * ) 'coupling mode "', TRIM( coupling_mode ), &
                 '": dt_restart = ', dt_restart, '& is not equal to ',         &
                 'dt_restart_remote = ', remote
          CALL message( 'check_parameters', 'PA0007', 1, 2, 0, 6, 0 )
       ENDIF

       time_to_be_simulated_from_reference_point = end_time-coupling_start_time

       IF  ( myid == 0 ) THEN
          CALL MPI_SEND( time_to_be_simulated_from_reference_point, 1,         &
                         MPI_REAL, target_id, 14, comm_inter, ierr )
          CALL MPI_RECV( remote, 1, MPI_REAL, target_id, 14, comm_inter,       &
                         status, ierr )
       ENDIF
       CALL MPI_BCAST( remote, 1, MPI_REAL, 0, comm2d, ierr)

       IF ( time_to_be_simulated_from_reference_point /= remote )  THEN
          WRITE( message_string, * ) 'coupling mode "', TRIM( coupling_mode ), &
                 '": time_to_be_simulated_from_reference_point = ',            &
                 time_to_be_simulated_from_reference_point, '& is not equal ', &
                 'to time_to_be_simulated_from_reference_point_remote = ',     &
                 remote
          CALL message( 'check_parameters', 'PA0008', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( myid == 0 ) THEN
          CALL MPI_SEND( dx, 1, MPI_REAL, target_id, 15, comm_inter, ierr )
          CALL MPI_RECV( remote, 1, MPI_REAL, target_id, 15, comm_inter,       &
                                                             status, ierr )
       ENDIF
       CALL MPI_BCAST( remote, 1, MPI_REAL, 0, comm2d, ierr)


       IF ( coupling_mode == 'atmosphere_to_ocean') THEN

          IF ( dx < remote ) THEN
             WRITE( message_string, * ) 'coupling mode "',                     &
                   TRIM( coupling_mode ),                                      &
           '": dx in Atmosphere is not equal to or not larger than dx in ocean'
             CALL message( 'check_parameters', 'PA0009', 1, 2, 0, 6, 0 )
          ENDIF

          IF ( (nx_a+1)*dx /= (nx_o+1)*remote )  THEN
             WRITE( message_string, * ) 'coupling mode "',                     &
                    TRIM( coupling_mode ),                                     &
             '": Domain size in x-direction is not equal in ocean and atmosphere'
             CALL message( 'check_parameters', 'PA0010', 1, 2, 0, 6, 0 )
          ENDIF

       ENDIF

       IF ( myid == 0) THEN
          CALL MPI_SEND( dy, 1, MPI_REAL, target_id, 16, comm_inter, ierr )
          CALL MPI_RECV( remote, 1, MPI_REAL, target_id, 16, comm_inter,       &
                         status, ierr )
       ENDIF
       CALL MPI_BCAST( remote, 1, MPI_REAL, 0, comm2d, ierr)

       IF ( coupling_mode == 'atmosphere_to_ocean') THEN

          IF ( dy < remote )  THEN
             WRITE( message_string, * ) 'coupling mode "',                     &
                    TRIM( coupling_mode ),                                     &
                 '": dy in Atmosphere is not equal to or not larger than dy in ocean'
             CALL message( 'check_parameters', 'PA0011', 1, 2, 0, 6, 0 )
          ENDIF

          IF ( (ny_a+1)*dy /= (ny_o+1)*remote )  THEN
             WRITE( message_string, * ) 'coupling mode "',                     &
                   TRIM( coupling_mode ),                                      &
             '": Domain size in y-direction is not equal in ocean and atmosphere'
             CALL message( 'check_parameters', 'PA0012', 1, 2, 0, 6, 0 )
          ENDIF

          IF ( MOD(nx_o+1,nx_a+1) /= 0 )  THEN
             WRITE( message_string, * ) 'coupling mode "',                     &
                   TRIM( coupling_mode ),                                      &
             '": nx+1 in ocean is not divisible by nx+1 in',                   &
             ' atmosphere without remainder'
             CALL message( 'check_parameters', 'PA0339', 1, 2, 0, 6, 0 )
          ENDIF

          IF ( MOD(ny_o+1,ny_a+1) /= 0 )  THEN
             WRITE( message_string, * ) 'coupling mode "',                     &
                   TRIM( coupling_mode ),                                      &
             '": ny+1 in ocean is not divisible by ny+1 in', &
             ' atmosphere without remainder'

             CALL message( 'check_parameters', 'PA0340', 1, 2, 0, 6, 0 )
          ENDIF

       ENDIF
#else
       WRITE( message_string, * ) 'coupling requires PALM to be compiled with',&
            ' cpp-option "-D__parallel"'
       CALL message( 'check_parameters', 'PA0141', 1, 2, 0, 6, 0 )
#endif
    ENDIF

#if defined( __parallel )
!
!-- Exchange via intercommunicator
    IF ( coupling_mode == 'atmosphere_to_ocean' .AND. myid == 0 )  THEN
       CALL MPI_SEND( humidity, 1, MPI_LOGICAL, target_id, 19, comm_inter,     &
                      ierr )
    ELSEIF ( coupling_mode == 'ocean_to_atmosphere' .AND. myid == 0)  THEN
       CALL MPI_RECV( humidity_remote, 1, MPI_LOGICAL, target_id, 19,          &
                      comm_inter, status, ierr )
    ENDIF
    CALL MPI_BCAST( humidity_remote, 1, MPI_LOGICAL, 0, comm2d, ierr)

#endif

!
!-- User settings for restart times requires that "restart" has been given as
!-- file activation string. Otherwise, binary output would not be saved by
!-- palmrun.
    IF (  ( restart_time /= 9999999.9_wp  .OR.  dt_restart /= 9999999.9_wp )   &
         .AND.  .NOT. write_binary )  THEN
       WRITE( message_string, * ) 'manual restart settings requires file ',    &
                                  'activation string "restart"'
       CALL message( 'check_parameters', 'PA0001', 1, 2, 0, 6, 0 )
    ENDIF


!
!-- Generate the file header which is used as a header for most of PALM's
!-- output files
    CALL DATE_AND_TIME( date, time, run_zone )
    run_date = date(1:4)//'-'//date(5:6)//'-'//date(7:8)
    run_time = time(1:2)//':'//time(3:4)//':'//time(5:6)
    IF ( coupling_mode == 'uncoupled' )  THEN
       coupling_string = ''
    ELSEIF ( coupling_mode == 'vnested_crse' )  THEN
       coupling_string = ' nested (coarse)'
    ELSEIF ( coupling_mode == 'vnested_fine' )  THEN
       coupling_string = ' nested (fine)'
    ELSEIF ( coupling_mode == 'atmosphere_to_ocean' )  THEN
       coupling_string = ' coupled (atmosphere)'
    ELSEIF ( coupling_mode == 'ocean_to_atmosphere' )  THEN
       coupling_string = ' coupled (ocean)'
    ENDIF
    IF ( ensemble_member_nr /= 0 )  THEN
       WRITE( ensemble_string, '(2X,A,I2.2)' )  'en-no: ', ensemble_member_nr
    ELSE
       ensemble_string = ''
    ENDIF
    IF ( nested_run )  THEN
       WRITE( nest_string, '(2X,A,I2.2)' )  'nest-id: ', cpl_id
    ELSE
       nest_string = ''
    ENDIF

    WRITE ( run_description_header,                                            &
            '(A,2X,A,2X,A,A,A,I2.2,A,A,A,2X,A,A,2X,A,1X,A)' )                  &
          TRIM( version ), TRIM( revision ), 'run: ',                          &
          TRIM( run_identifier ), '.', runnr, TRIM( coupling_string ),         &
          TRIM( nest_string ), TRIM( ensemble_string), 'host: ', TRIM( host ), &
          run_date, run_time

!
!-- Check the general loop optimization method
    SELECT CASE ( TRIM( loop_optimization ) )

       CASE ( 'cache', 'vector' )
          CONTINUE

       CASE DEFAULT
          message_string = 'illegal value given for loop_optimization: "' //   &
                           TRIM( loop_optimization ) // '"'
          CALL message( 'check_parameters', 'PA0013', 1, 2, 0, 6, 0 )

    END SELECT

!
!-- Check topography setting (check for illegal parameter combinations)
    IF ( topography /= 'flat' )  THEN
       action = ' '
       IF ( scalar_advec /= 'pw-scheme' .AND. scalar_advec /= 'ws-scheme'      &
          )  THEN
          WRITE( action, '(A,A)' )  'scalar_advec = ', scalar_advec
       ENDIF
       IF ( momentum_advec /= 'pw-scheme' .AND. momentum_advec /= 'ws-scheme' )&
       THEN
          WRITE( action, '(A,A)' )  'momentum_advec = ', momentum_advec
       ENDIF
       IF ( psolver == 'sor' )  THEN
          WRITE( action, '(A,A)' )  'psolver = ', psolver
       ENDIF
       IF ( sloping_surface )  THEN
          WRITE( action, '(A)' )  'sloping surface = .TRUE.'
       ENDIF
       IF ( galilei_transformation )  THEN
          WRITE( action, '(A)' )  'galilei_transformation = .TRUE.'
       ENDIF
       IF ( cloud_droplets )  THEN
          WRITE( action, '(A)' )  'cloud_droplets = .TRUE.'
       ENDIF
       IF ( .NOT. constant_flux_layer .AND. topography /= 'closed_channel' )   &
       THEN
          WRITE( action, '(A)' )  'constant_flux_layer = .FALSE.'
       ENDIF
       IF ( action /= ' ' )  THEN
          message_string = 'The specified topography does not allow ' //       &
                           TRIM( action )
          CALL message( 'check_parameters', 'PA0014', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Check illegal/untested parameter combinations for closed channel
       If ( topography == 'closed_channel' ) THEN 
          symmetry_flag = 1
          message_string = 'Bottom and top boundary are treated equal'
          CALL message( 'check_parameters', 'PA0699', 0, 0, 0, 6, 0 )
       
          IF ( dz(1) /= dz(COUNT( dz /= -1.0_wp )) .OR.                        &
               dz_stretch_level /= -9999999.9_wp) THEN
             WRITE( message_string, * )  'dz should be equal close to the ' // &
                                         'boundaries due to symmetrical problem'
             CALL message( 'check_parameters', 'PA0700', 1, 2, 0, 6, 0 )
          ENDIF
       
          IF ( constant_flux_layer ) THEN
             WRITE( message_string, * )  'A constant flux layer is not '//     &
                                         'allowed if a closed channel '//      &
                                         'shall be used'
             CALL message( 'check_parameters', 'PA0701', 1, 2, 0, 6, 0 )
          ENDIF
       
          IF ( ocean_mode ) THEN
             WRITE( message_string, * )  'The ocean mode is not allowed if '// &
                                         'a closed channel shall be used'
             CALL message( 'check_parameters', 'PA0702', 1, 2, 0, 6, 0 )
          ENDIF
       
          IF ( momentum_advec /= 'ws-scheme' .OR.                              &
               scalar_advec /= 'ws-scheme' ) THEN
             WRITE( message_string, * )  'A closed channel require the '//     &
                                         'upwind scheme of Wicker and ' //     &
                                         'Skamarock as the advection scheme'
             CALL message( 'check_parameters', 'PA0703', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF
    ENDIF

!
!-- Check approximation
    IF ( TRIM( approximation ) /= 'boussinesq'   .AND.                         &
         TRIM( approximation ) /= 'anelastic' )  THEN
       message_string = 'unknown approximation: approximation = "' //          &
                        TRIM( approximation ) // '"'
       CALL message( 'check_parameters', 'PA0446', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Check approximation requirements
    IF ( TRIM( approximation ) == 'anelastic'   .AND.                          &
         TRIM( momentum_advec ) /= 'ws-scheme' )  THEN
       message_string = 'Anelastic approximation requires ' //                 &
                        'momentum_advec = "ws-scheme"'
       CALL message( 'check_parameters', 'PA0447', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( TRIM( approximation ) == 'anelastic'   .AND.                          &
         TRIM( psolver ) == 'multigrid' )  THEN
       message_string = 'Anelastic approximation currently only supports ' //  &
                        'psolver = "poisfft", ' //                             &
                        'psolver = "sor" and ' //                              &
                        'psolver = "multigrid_noopt"'
       CALL message( 'check_parameters', 'PA0448', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( TRIM( approximation ) == 'anelastic'   .AND.                          &
         conserve_volume_flow )  THEN
       message_string = 'Anelastic approximation is not allowed with ' //      &
                        'conserve_volume_flow = .TRUE.'
       CALL message( 'check_parameters', 'PA0449', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Check flux input mode
    IF ( TRIM( flux_input_mode ) /= 'dynamic'    .AND.                         &
         TRIM( flux_input_mode ) /= 'kinematic'  .AND.                         &
         TRIM( flux_input_mode ) /= 'approximation-specific' )  THEN
       message_string = 'unknown flux input mode: flux_input_mode = "' //      &
                        TRIM( flux_input_mode ) // '"'
       CALL message( 'check_parameters', 'PA0450', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Set flux input mode according to approximation if applicable
    IF ( TRIM( flux_input_mode ) == 'approximation-specific' )  THEN
       IF ( TRIM( approximation ) == 'anelastic' )  THEN
          flux_input_mode = 'dynamic'
       ELSEIF ( TRIM( approximation ) == 'boussinesq' )  THEN
          flux_input_mode = 'kinematic'
       ENDIF
    ENDIF

!
!-- Check flux output mode
    IF ( TRIM( flux_output_mode ) /= 'dynamic'    .AND.                        &
         TRIM( flux_output_mode ) /= 'kinematic'  .AND.                        &
         TRIM( flux_output_mode ) /= 'approximation-specific' )  THEN
       message_string = 'unknown flux output mode: flux_output_mode = "' //    &
                        TRIM( flux_output_mode ) // '"'
       CALL message( 'check_parameters', 'PA0451', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Set flux output mode according to approximation if applicable
    IF ( TRIM( flux_output_mode ) == 'approximation-specific' )  THEN
       IF ( TRIM( approximation ) == 'anelastic' )  THEN
          flux_output_mode = 'dynamic'
       ELSEIF ( TRIM( approximation ) == 'boussinesq' )  THEN
          flux_output_mode = 'kinematic'
       ENDIF
    ENDIF


!
!-- When the land- or urban-surface model is used, the flux output must be 
!-- dynamic.
    IF ( land_surface  .OR.  urban_surface )  THEN
       flux_output_mode = 'dynamic'
    ENDIF

!
!-- Set the flux output units according to flux_output_mode
    IF ( TRIM( flux_output_mode ) == 'kinematic' ) THEN
        heatflux_output_unit              = 'K m/s'
        waterflux_output_unit             = 'kg/kg m/s'
        momentumflux_output_unit          = 'm2/s2'
    ELSEIF ( TRIM( flux_output_mode ) == 'dynamic' ) THEN
        heatflux_output_unit              = 'W/m2'
        waterflux_output_unit             = 'W/m2'
        momentumflux_output_unit          = 'N/m2'
    ENDIF

!
!-- set time series output units for fluxes
    dots_unit(14:16) = TRIM( heatflux_output_unit )
    dots_unit(21)    = TRIM( waterflux_output_unit )
    dots_unit(19:20) = TRIM( momentumflux_output_unit )

!
!-- Add other module specific timeseries
    CALL module_interface_check_data_output_ts( dots_max, dots_num, dots_label, dots_unit )

!
!-- Check if maximum number of allowed timeseries is exceeded
    IF ( dots_num > dots_max )  THEN
       WRITE( message_string, * ) 'number of time series quantities exceeds',  &
                                  ' its maximum of dots_max = ', dots_max,     &
                                  '&Please increase dots_max in modules.f90.'
       CALL message( 'init_3d_model', 'PA0194', 1, 2, 0, 6, 0 )    
    ENDIF

!
!-- Check whether there are any illegal values
!-- Pressure solver:
    IF ( psolver /= 'poisfft'  .AND.  psolver /= 'sor'  .AND.                  &
         psolver /= 'multigrid'  .AND.  psolver /= 'multigrid_noopt' )  THEN
       message_string = 'unknown solver for perturbation pressure: psolver' // &
                        ' = "' // TRIM( psolver ) // '"'
       CALL message( 'check_parameters', 'PA0016', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( psolver(1:9) == 'multigrid' )  THEN
       IF ( cycle_mg == 'w' )  THEN
          gamma_mg = 2
       ELSEIF ( cycle_mg == 'v' )  THEN
          gamma_mg = 1
       ELSE
          message_string = 'unknown multigrid cycle: cycle_mg = "' //          &
                           TRIM( cycle_mg ) // '"'
          CALL message( 'check_parameters', 'PA0020', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

    IF ( fft_method /= 'singleton-algorithm'  .AND.                            &
         fft_method /= 'temperton-algorithm'  .AND.                            &
         fft_method /= 'fftw'                 .AND.                            &
         fft_method /= 'system-specific' )  THEN
       message_string = 'unknown fft-algorithm: fft_method = "' //             &
                        TRIM( fft_method ) // '"'
       CALL message( 'check_parameters', 'PA0021', 1, 2, 0, 6, 0 )
    ENDIF

    IF( momentum_advec == 'ws-scheme' .AND.                                    &
        .NOT. call_psolver_at_all_substeps  ) THEN
        message_string = 'psolver must be called at each RK3 substep when "'// &
                      TRIM(momentum_advec) // ' "is used for momentum_advec'
        CALL message( 'check_parameters', 'PA0344', 1, 2, 0, 6, 0 )
    END IF
!
!-- Advection schemes:
    IF ( momentum_advec /= 'pw-scheme'  .AND.                                  &  
         momentum_advec /= 'ws-scheme'  .AND.                                  &
         momentum_advec /= 'up-scheme' )                                       &
    THEN
       message_string = 'unknown advection scheme: momentum_advec = "' //      &
                        TRIM( momentum_advec ) // '"'
       CALL message( 'check_parameters', 'PA0022', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( ( momentum_advec == 'ws-scheme' .OR.  scalar_advec == 'ws-scheme' )   &
           .AND. ( timestep_scheme == 'euler' .OR.                             &
                   timestep_scheme == 'runge-kutta-2' ) )                      &
    THEN
       message_string = 'momentum_advec or scalar_advec = "'                   &
         // TRIM( momentum_advec ) // '" is not allowed with ' //              &
         'timestep_scheme = "' // TRIM( timestep_scheme ) // '"'
       CALL message( 'check_parameters', 'PA0023', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( scalar_advec /= 'pw-scheme'  .AND.  scalar_advec /= 'ws-scheme' .AND. &
         scalar_advec /= 'bc-scheme' .AND. scalar_advec /= 'up-scheme' )       &
    THEN
       message_string = 'unknown advection scheme: scalar_advec = "' //        &
                        TRIM( scalar_advec ) // '"'
       CALL message( 'check_parameters', 'PA0024', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( scalar_advec == 'bc-scheme'  .AND.  loop_optimization == 'cache' )    &
    THEN
       message_string = 'advection_scheme scalar_advec = "'                    &
         // TRIM( scalar_advec ) // '" not implemented for ' //                &
         'loop_optimization = "' // TRIM( loop_optimization ) // '"'
       CALL message( 'check_parameters', 'PA0026', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( use_sgs_for_particles  .AND.  .NOT. cloud_droplets  .AND.             &
         .NOT. use_upstream_for_tke  .AND.                                     &
         scalar_advec /= 'ws-scheme'                                           &
       )  THEN
       use_upstream_for_tke = .TRUE.
       message_string = 'use_upstream_for_tke is set to .TRUE. because ' //    &
                        'use_sgs_for_particles = .TRUE. '          //          &
                        'and scalar_advec /= ws-scheme'
       CALL message( 'check_parameters', 'PA0025', 0, 1, 0, 6, 0 )
    ENDIF

!
!-- Set LOGICAL switches to enhance performance
    IF ( momentum_advec == 'ws-scheme' )  ws_scheme_mom = .TRUE.
    IF ( scalar_advec   == 'ws-scheme' )  ws_scheme_sca = .TRUE.


!
!-- Timestep schemes:
    SELECT CASE ( TRIM( timestep_scheme ) )

       CASE ( 'euler' )
          intermediate_timestep_count_max = 1

       CASE ( 'runge-kutta-2' )
          intermediate_timestep_count_max = 2

       CASE ( 'runge-kutta-3' )
          intermediate_timestep_count_max = 3

       CASE DEFAULT
          message_string = 'unknown timestep scheme: timestep_scheme = "' //   &
                           TRIM( timestep_scheme ) // '"'
          CALL message( 'check_parameters', 'PA0027', 1, 2, 0, 6, 0 )

    END SELECT

    IF ( (momentum_advec /= 'pw-scheme' .AND. momentum_advec /= 'ws-scheme')   &
         .AND. timestep_scheme(1:5) == 'runge' ) THEN
       message_string = 'momentum advection scheme "' // &
                        TRIM( momentum_advec ) // '" & does not work with ' // &
                        'timestep_scheme "' // TRIM( timestep_scheme ) // '"'
       CALL message( 'check_parameters', 'PA0029', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Check for proper settings for microphysics
    IF ( bulk_cloud_model  .AND.  cloud_droplets )  THEN
       message_string = 'bulk_cloud_model = .TRUE. is not allowed with ' //    &
                        'cloud_droplets = .TRUE.'
       CALL message( 'check_parameters', 'PA0442', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Initializing actions must have been set by the user
    IF ( TRIM( initializing_actions ) == '' )  THEN
       message_string = 'no value specified for initializing_actions'
       CALL message( 'check_parameters', 'PA0149', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( TRIM( initializing_actions ) /= 'read_restart_data'  .AND.            &
         TRIM( initializing_actions ) /= 'cyclic_fill' )  THEN
!
!--    No restart run: several initialising actions are possible
       action = initializing_actions
       DO  WHILE ( TRIM( action ) /= '' )
          position = INDEX( action, ' ' )
          SELECT CASE ( action(1:position-1) )

             CASE ( 'set_constant_profiles', 'set_1d-model_profiles',          &
                    'by_user', 'initialize_vortex', 'initialize_ptanom',       &
                    'initialize_bubble', 'inifor' )
                action = action(position+1:)

             CASE DEFAULT
                message_string = 'initializing_action = "' //                  &
                                 TRIM( action ) // '" unknown or not allowed'
                CALL message( 'check_parameters', 'PA0030', 1, 2, 0, 6, 0 )

          END SELECT
       ENDDO
    ENDIF

    IF ( TRIM( initializing_actions ) == 'initialize_vortex'  .AND.            &
         conserve_volume_flow ) THEN
         message_string = 'initializing_actions = "initialize_vortex"' //      &
                        ' is not allowed with conserve_volume_flow = .T.'
       CALL message( 'check_parameters', 'PA0343', 1, 2, 0, 6, 0 )
    ENDIF


    IF ( INDEX( initializing_actions, 'set_constant_profiles' ) /= 0  .AND.    &
         INDEX( initializing_actions, 'set_1d-model_profiles' ) /= 0 )  THEN
       message_string = 'initializing_actions = "set_constant_profiles"' //    &
                        ' and "set_1d-model_profiles" are not allowed ' //     &
                        'simultaneously'
       CALL message( 'check_parameters', 'PA0031', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( INDEX( initializing_actions, 'set_constant_profiles' ) /= 0  .AND.    &
         INDEX( initializing_actions, 'by_user' ) /= 0 )  THEN
       message_string = 'initializing_actions = "set_constant_profiles"' //    &
                        ' and "by_user" are not allowed simultaneously'
       CALL message( 'check_parameters', 'PA0032', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( INDEX( initializing_actions, 'by_user' ) /= 0  .AND.                  &
         INDEX( initializing_actions, 'set_1d-model_profiles' ) /= 0 )  THEN
       message_string = 'initializing_actions = "by_user" and ' //             &
                        '"set_1d-model_profiles" are not allowed simultaneously'
       CALL message( 'check_parameters', 'PA0033', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- In case of spinup and nested run, spinup end time must be identical 
!-- in order to have synchronously running simulations. 
    IF ( nested_run )  THEN
#if defined( __parallel )
       CALL MPI_ALLREDUCE( spinup_time, spinup_time_max, 1, MPI_REAL,          &
                           MPI_MAX, MPI_COMM_WORLD, ierr )
       CALL MPI_ALLREDUCE( dt_spinup,   dt_spinup_max,   1, MPI_REAL,          &
                           MPI_MAX, MPI_COMM_WORLD, ierr )

       IF ( spinup_time /= spinup_time_max  .OR.  dt_spinup /= dt_spinup_max ) &
       THEN
          message_string = 'In case of nesting, spinup_time and ' //           &
                           'dt_spinup must be identical in all parent ' //     &
                           'and child domains.'
          CALL message( 'check_parameters', 'PA0489', 3, 2, 0, 6, 0 )
       ENDIF
#endif
    ENDIF

    IF ( bulk_cloud_model  .AND.  .NOT.  humidity )  THEN
       WRITE( message_string, * ) 'bulk_cloud_model = ', bulk_cloud_model,     &
              ' is not allowed with humidity = ', humidity
       CALL message( 'check_parameters', 'PA0034', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( humidity  .AND.  sloping_surface )  THEN
       message_string = 'humidity = .TRUE. and sloping_surface = .TRUE. ' //   &
                        'are not allowed simultaneously'
       CALL message( 'check_parameters', 'PA0036', 1, 2, 0, 6, 0 )
    ENDIF

!-- Check the module settings
    CALL module_interface_check_parameters

!
!-- In case of no restart run, check initialising parameters and calculate
!-- further quantities
    IF ( TRIM( initializing_actions ) /= 'read_restart_data' )  THEN

!
!--    Initial profiles for 1D and 3D model, respectively (u,v further below)
       pt_init = pt_surface
       IF ( humidity       )  q_init  = q_surface
       IF ( passive_scalar )  s_init  = s_surface

!--
!--    If required, compute initial profile of the geostrophic wind
!--    (component ug)
       i = 1
       gradient = 0.0_wp

       IF ( .NOT. ocean_mode )  THEN

          ug_vertical_gradient_level_ind(1) = 0
          ug(0) = ug_surface
          DO  k = 1, nzt+1
             IF ( i < 11 )  THEN
                IF ( ug_vertical_gradient_level(i) < zu(k)  .AND.              &
                     ug_vertical_gradient_level(i) >= 0.0_wp )  THEN
                   gradient = ug_vertical_gradient(i) / 100.0_wp
                   ug_vertical_gradient_level_ind(i) = k - 1
                   i = i + 1
                ENDIF
             ENDIF
             IF ( gradient /= 0.0_wp )  THEN
                IF ( k /= 1 )  THEN
                   ug(k) = ug(k-1) + dzu(k) * gradient
                ELSE
                   ug(k) = ug_surface + dzu(k) * gradient
                ENDIF
             ELSE
                ug(k) = ug(k-1)
             ENDIF
          ENDDO

       ELSE

          ug_vertical_gradient_level_ind(1) = nzt+1
          ug(nzt+1) = ug_surface
          DO  k = nzt, nzb, -1
             IF ( i < 11 )  THEN
                IF ( ug_vertical_gradient_level(i) > zu(k)  .AND.              &
                     ug_vertical_gradient_level(i) <= 0.0_wp )  THEN
                   gradient = ug_vertical_gradient(i) / 100.0_wp
                   ug_vertical_gradient_level_ind(i) = k + 1
                   i = i + 1
                ENDIF
             ENDIF
             IF ( gradient /= 0.0_wp )  THEN
                IF ( k /= nzt )  THEN
                   ug(k) = ug(k+1) - dzu(k+1) * gradient
                ELSE
                   ug(k)   = ug_surface - 0.5_wp * dzu(k+1) * gradient
                   ug(k+1) = ug_surface + 0.5_wp * dzu(k+1) * gradient
                ENDIF
             ELSE
                ug(k) = ug(k+1)
             ENDIF
          ENDDO

       ENDIF

!
!--    In case of no given gradients for ug, choose a zero gradient
       IF ( ug_vertical_gradient_level(1) == -9999999.9_wp )  THEN
          ug_vertical_gradient_level(1) = 0.0_wp
       ENDIF

!
!--
!--    If required, compute initial profile of the geostrophic wind
!--    (component vg)
       i = 1
       gradient = 0.0_wp

       IF ( .NOT. ocean_mode )  THEN

          vg_vertical_gradient_level_ind(1) = 0
          vg(0) = vg_surface
          DO  k = 1, nzt+1
             IF ( i < 11 )  THEN
                IF ( vg_vertical_gradient_level(i) < zu(k)  .AND.              &
                     vg_vertical_gradient_level(i) >= 0.0_wp )  THEN
                   gradient = vg_vertical_gradient(i) / 100.0_wp
                   vg_vertical_gradient_level_ind(i) = k - 1
                   i = i + 1
                ENDIF
             ENDIF
             IF ( gradient /= 0.0_wp )  THEN
                IF ( k /= 1 )  THEN
                   vg(k) = vg(k-1) + dzu(k) * gradient
                ELSE
                   vg(k) = vg_surface + dzu(k) * gradient
                ENDIF
             ELSE
                vg(k) = vg(k-1)
             ENDIF
          ENDDO

       ELSE

          vg_vertical_gradient_level_ind(1) = nzt+1
          vg(nzt+1) = vg_surface
          DO  k = nzt, nzb, -1
             IF ( i < 11 )  THEN
                IF ( vg_vertical_gradient_level(i) > zu(k)  .AND.              &
                     vg_vertical_gradient_level(i) <= 0.0_wp )  THEN
                   gradient = vg_vertical_gradient(i) / 100.0_wp
                   vg_vertical_gradient_level_ind(i) = k + 1
                   i = i + 1
                ENDIF
             ENDIF
             IF ( gradient /= 0.0_wp )  THEN
                IF ( k /= nzt )  THEN
                   vg(k) = vg(k+1) - dzu(k+1) * gradient
                ELSE
                   vg(k)   = vg_surface - 0.5_wp * dzu(k+1) * gradient
                   vg(k+1) = vg_surface + 0.5_wp * dzu(k+1) * gradient
                ENDIF
             ELSE
                vg(k) = vg(k+1)
             ENDIF
          ENDDO

       ENDIF

!
!--    In case of no given gradients for vg, choose a zero gradient
       IF ( vg_vertical_gradient_level(1) == -9999999.9_wp )  THEN
          vg_vertical_gradient_level(1) = 0.0_wp
       ENDIF

!
!--    Let the initial wind profiles be the calculated ug/vg profiles or
!--    interpolate them from wind profile data (if given)
       IF ( u_profile(1) == 9999999.9_wp  .AND.  v_profile(1) == 9999999.9_wp )  THEN

          u_init = ug
          v_init = vg

       ELSEIF ( u_profile(1) == 0.0_wp  .AND.  v_profile(1) == 0.0_wp )  THEN

          IF ( uv_heights(1) /= 0.0_wp )  THEN
             message_string = 'uv_heights(1) must be 0.0'
             CALL message( 'check_parameters', 'PA0345', 1, 2, 0, 6, 0 )
          ENDIF

          IF ( omega /= 0.0_wp )  THEN
             message_string = 'Coriolis force must be switched off (by setting omega=0.0)' //  &
                              ' when prescribing the forcing by u_profile and v_profile'
             CALL message( 'check_parameters', 'PA0347', 1, 2, 0, 6, 0 )
          ENDIF

          use_prescribed_profile_data = .TRUE.

          kk = 1
          u_init(0) = 0.0_wp
          v_init(0) = 0.0_wp

          DO  k = 1, nz+1

             IF ( kk < 200 )  THEN
                DO  WHILE ( uv_heights(kk+1) <= zu(k) )
                   kk = kk + 1
                   IF ( kk == 200 )  EXIT
                ENDDO
             ENDIF

             IF ( kk < 200  .AND.  uv_heights(kk+1) /= 9999999.9_wp )  THEN
                u_init(k) = u_profile(kk) + ( zu(k) - uv_heights(kk) ) /       &
                                       ( uv_heights(kk+1) - uv_heights(kk) ) * &
                                       ( u_profile(kk+1) - u_profile(kk) )
                v_init(k) = v_profile(kk) + ( zu(k) - uv_heights(kk) ) /       &
                                       ( uv_heights(kk+1) - uv_heights(kk) ) * &
                                       ( v_profile(kk+1) - v_profile(kk) )
             ELSE
                u_init(k) = u_profile(kk)
                v_init(k) = v_profile(kk)
             ENDIF

          ENDDO

       ELSE

          message_string = 'u_profile(1) and v_profile(1) must be 0.0'
          CALL message( 'check_parameters', 'PA0346', 1, 2, 0, 6, 0 )

       ENDIF

!
!--    Compute initial temperature profile using the given temperature gradients
       IF (  .NOT.  neutral )  THEN
          CALL init_vertical_profiles( pt_vertical_gradient_level_ind,          &
                                       pt_vertical_gradient_level,              &
                                       pt_vertical_gradient, pt_init,           &
                                       pt_surface, bc_pt_t_val )
       ENDIF
!
!--    Compute initial humidity profile using the given humidity gradients
       IF ( humidity )  THEN
          CALL init_vertical_profiles( q_vertical_gradient_level_ind,          &
                                       q_vertical_gradient_level,              &
                                       q_vertical_gradient, q_init,            &
                                       q_surface, bc_q_t_val )
       ENDIF
!
!--    Compute initial scalar profile using the given scalar gradients
       IF ( passive_scalar )  THEN
          CALL init_vertical_profiles( s_vertical_gradient_level_ind,          &
                                       s_vertical_gradient_level,              &
                                       s_vertical_gradient, s_init,            &
                                       s_surface, bc_s_t_val )
       ENDIF
!
!--    TODO
!--    Compute initial chemistry profile using the given chemical species gradients
!--    Russo: Is done in chem_init --> kanani: Revise

    ENDIF

!
!-- Check if the control parameter use_subsidence_tendencies is used correctly
    IF ( use_subsidence_tendencies  .AND.  .NOT.  large_scale_subsidence )  THEN
       message_string = 'The usage of use_subsidence_tendencies ' //           &
                            'requires large_scale_subsidence = .T..'
       CALL message( 'check_parameters', 'PA0396', 1, 2, 0, 6, 0 )
    ELSEIF ( use_subsidence_tendencies  .AND.  .NOT. large_scale_forcing )  THEN
       message_string = 'The usage of use_subsidence_tendencies ' //           &
                            'requires large_scale_forcing = .T..'
       CALL message( 'check_parameters', 'PA0397', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Initialize large scale subsidence if required
    If ( large_scale_subsidence )  THEN
       IF ( subs_vertical_gradient_level(1) /= -9999999.9_wp  .AND.            &
                                     .NOT.  large_scale_forcing )  THEN
          CALL init_w_subsidence
       ENDIF
!
!--    In case large_scale_forcing is used, profiles for subsidence velocity
!--    are read in from file LSF_DATA

       IF ( subs_vertical_gradient_level(1) == -9999999.9_wp  .AND.            &
            .NOT.  large_scale_forcing )  THEN
          message_string = 'There is no default large scale vertical ' //      &
                           'velocity profile set. Specify the subsidence ' //  &
                           'velocity profile via subs_vertical_gradient ' //   &
                           'and subs_vertical_gradient_level.'
          CALL message( 'check_parameters', 'PA0380', 1, 2, 0, 6, 0 )
       ENDIF
    ELSE
        IF ( subs_vertical_gradient_level(1) /= -9999999.9_wp )  THEN
           message_string = 'Enable usage of large scale subsidence by ' //    &
                            'setting large_scale_subsidence = .T..'
          CALL message( 'check_parameters', 'PA0381', 1, 2, 0, 6, 0 )
        ENDIF
    ENDIF

!
!-- Overwrite parameters from namelist if necessary and compute Coriolis parameter. 
!-- @todo - move initialization of f and fs to coriolis_mod.
    IF ( input_pids_static )  THEN
       latitude       = init_model%latitude
       longitude      = init_model%longitude
       rotation_angle = init_model%rotation_angle
    ENDIF

    f  = 2.0_wp * omega * SIN( latitude / 180.0_wp * pi )
    fs = 2.0_wp * omega * COS( latitude / 180.0_wp * pi )

!
!-- Check and set buoyancy related parameters and switches
    IF ( reference_state == 'horizontal_average' )  THEN
       CONTINUE
    ELSEIF ( reference_state == 'initial_profile' )  THEN
       use_initial_profile_as_reference = .TRUE.
    ELSEIF ( reference_state == 'single_value' )  THEN
       use_single_reference_value = .TRUE.
       IF ( pt_reference == 9999999.9_wp )  pt_reference = pt_surface
       vpt_reference = pt_reference * ( 1.0_wp + 0.61_wp * q_surface )
    ELSE
       message_string = 'illegal value for reference_state: "' //              &
                        TRIM( reference_state ) // '"'
       CALL message( 'check_parameters', 'PA0056', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- In case of a given slope, compute the relevant quantities
    IF ( alpha_surface /= 0.0_wp )  THEN
       IF ( ABS( alpha_surface ) > 90.0_wp )  THEN
          WRITE( message_string, * ) 'ABS( alpha_surface = ', alpha_surface,   &
                                     ' ) must be < 90.0'
          CALL message( 'check_parameters', 'PA0043', 1, 2, 0, 6, 0 )
       ENDIF
       sloping_surface = .TRUE.
       cos_alpha_surface = COS( alpha_surface / 180.0_wp * pi )
       sin_alpha_surface = SIN( alpha_surface / 180.0_wp * pi )
    ENDIF

!
!-- Check time step and cfl_factor
    IF ( dt /= -1.0_wp )  THEN
       IF ( dt <= 0.0_wp )  THEN
          WRITE( message_string, * ) 'dt = ', dt , ' <= 0.0'
          CALL message( 'check_parameters', 'PA0044', 1, 2, 0, 6, 0 )
       ENDIF
       dt_3d = dt
       dt_fixed = .TRUE.
    ENDIF

    IF ( cfl_factor <= 0.0_wp  .OR.  cfl_factor > 1.0_wp )  THEN
       IF ( cfl_factor == -1.0_wp )  THEN
          IF ( timestep_scheme == 'runge-kutta-2' )  THEN
             cfl_factor = 0.8_wp
          ELSEIF ( timestep_scheme == 'runge-kutta-3' )  THEN
             cfl_factor = 0.9_wp
          ELSE
             cfl_factor = 0.9_wp
          ENDIF
       ELSE
          WRITE( message_string, * ) 'cfl_factor = ', cfl_factor,              &
                 ' out of range &0.0 < cfl_factor <= 1.0 is required'
          CALL message( 'check_parameters', 'PA0045', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

!
!-- Store simulated time at begin
    simulated_time_at_begin = simulated_time

!
!-- Store reference time for coupled runs and change the coupling flag,
!-- if ...
    IF ( simulated_time == 0.0_wp )  THEN
       IF ( coupling_start_time == 0.0_wp )  THEN
          time_since_reference_point = 0.0_wp
       ELSEIF ( time_since_reference_point < 0.0_wp )  THEN
          run_coupled = .FALSE.
       ENDIF
    ENDIF

!
!-- Set wind speed in the Galilei-transformed system
    IF ( galilei_transformation )  THEN
       IF ( use_ug_for_galilei_tr                    .AND.                     &
            ug_vertical_gradient_level(1) == 0.0_wp  .AND.                     &
            ug_vertical_gradient(1) == 0.0_wp        .AND.                     &
            vg_vertical_gradient_level(1) == 0.0_wp  .AND.                     &
            vg_vertical_gradient(1) == 0.0_wp )  THEN
          u_gtrans = ug_surface * 0.6_wp
          v_gtrans = vg_surface * 0.6_wp
       ELSEIF ( use_ug_for_galilei_tr  .AND.                                   &
                ( ug_vertical_gradient_level(1) /= 0.0_wp  .OR.                &
                ug_vertical_gradient(1) /= 0.0_wp ) )  THEN
          message_string = 'baroclinity (ug) not allowed simultaneously' //    &
                           ' with galilei transformation'
          CALL message( 'check_parameters', 'PA0046', 1, 2, 0, 6, 0 )
       ELSEIF ( use_ug_for_galilei_tr  .AND.                                   &
                ( vg_vertical_gradient_level(1) /= 0.0_wp  .OR.                &
                vg_vertical_gradient(1) /= 0.0_wp ) )  THEN
          message_string = 'baroclinity (vg) not allowed simultaneously' //    &
                           ' with galilei transformation'
          CALL message( 'check_parameters', 'PA0047', 1, 2, 0, 6, 0 )
       ELSE
          message_string = 'variable translation speed used for Galilei-' //   &
             'transformation, which may cause & instabilities in stably ' //   &
             'stratified regions'
          CALL message( 'check_parameters', 'PA0048', 0, 1, 0, 6, 0 )
       ENDIF
    ENDIF

!
!-- In case of using a prandtl-layer, calculated (or prescribed) surface
!-- fluxes have to be used in the diffusion-terms
    IF ( constant_flux_layer )  use_surface_fluxes = .TRUE.

!
!-- Check boundary conditions and set internal variables:
!-- Attention: the lateral boundary conditions have been already checked in
!-- parin
!
!-- Non-cyclic lateral boundaries require the multigrid method and Piascek-
!-- Willimas or Wicker - Skamarock advection scheme. Several schemes
!-- and tools do not work with non-cyclic boundary conditions.
    IF ( bc_lr /= 'cyclic'  .OR.  bc_ns /= 'cyclic' )  THEN
       IF ( psolver(1:9) /= 'multigrid' )  THEN
          message_string = 'non-cyclic lateral boundaries do not allow ' //    &
                           'psolver = "' // TRIM( psolver ) // '"'
          CALL message( 'check_parameters', 'PA0051', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( momentum_advec /= 'pw-scheme'  .AND.                               &
            momentum_advec /= 'ws-scheme' )  THEN

          message_string = 'non-cyclic lateral boundaries do not allow ' //    &
                           'momentum_advec = "' // TRIM( momentum_advec ) // '"'
          CALL message( 'check_parameters', 'PA0052', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( scalar_advec /= 'pw-scheme'  .AND.                                 &
            scalar_advec /= 'ws-scheme' )  THEN
          message_string = 'non-cyclic lateral boundaries do not allow ' //    &
                           'scalar_advec = "' // TRIM( scalar_advec ) // '"'
          CALL message( 'check_parameters', 'PA0053', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( galilei_transformation )  THEN
          message_string = 'non-cyclic lateral boundaries do not allow ' //    &
                           'galilei_transformation = .T.'
          CALL message( 'check_parameters', 'PA0054', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

!
!-- Bottom boundary condition for the turbulent Kinetic energy
    IF ( bc_e_b == 'neumann' )  THEN
       ibc_e_b = 1
    ELSEIF ( bc_e_b == '(u*)**2+neumann' )  THEN
       ibc_e_b = 2
       IF ( .NOT. constant_flux_layer )  THEN
          bc_e_b = 'neumann'
          ibc_e_b = 1
          message_string = 'boundary condition bc_e_b changed to "' //         &
                           TRIM( bc_e_b ) // '"'
          CALL message( 'check_parameters', 'PA0057', 0, 1, 0, 6, 0 )
       ENDIF
    ELSE
       message_string = 'unknown boundary condition: bc_e_b = "' //            &
                        TRIM( bc_e_b ) // '"'
       CALL message( 'check_parameters', 'PA0058', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Boundary conditions for perturbation pressure
    IF ( bc_p_b == 'dirichlet' )  THEN
       ibc_p_b = 0
    ELSEIF ( bc_p_b == 'neumann' )  THEN
       ibc_p_b = 1
    ELSE
       message_string = 'unknown boundary condition: bc_p_b = "' //            &
                        TRIM( bc_p_b ) // '"'
       CALL message( 'check_parameters', 'PA0059', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( bc_p_t == 'dirichlet' )  THEN
       ibc_p_t = 0
!-- TO_DO: later set bc_p_t to neumann before, in case of nested domain
    ELSEIF ( bc_p_t == 'neumann' .OR. bc_p_t == 'nested' )  THEN
       ibc_p_t = 1
    ELSE
       message_string = 'unknown boundary condition: bc_p_t = "' //            &
                        TRIM( bc_p_t ) // '"'
       CALL message( 'check_parameters', 'PA0061', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Boundary conditions for potential temperature
    IF ( coupling_mode == 'atmosphere_to_ocean' )  THEN
       ibc_pt_b = 2
    ELSE
       IF ( bc_pt_b == 'dirichlet' )  THEN
          ibc_pt_b = 0
       ELSEIF ( bc_pt_b == 'neumann' )  THEN
          ibc_pt_b = 1
       ELSE
          message_string = 'unknown boundary condition: bc_pt_b = "' //        &
                           TRIM( bc_pt_b ) // '"'
          CALL message( 'check_parameters', 'PA0062', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

    IF ( bc_pt_t == 'dirichlet' )  THEN
       ibc_pt_t = 0
    ELSEIF ( bc_pt_t == 'neumann' )  THEN
       ibc_pt_t = 1
    ELSEIF ( bc_pt_t == 'initial_gradient' )  THEN
       ibc_pt_t = 2
    ELSEIF ( bc_pt_t == 'nested'  .OR.  bc_pt_t == 'nesting_offline' )  THEN
       ibc_pt_t = 3
    ELSE
       message_string = 'unknown boundary condition: bc_pt_t = "' //           &
                        TRIM( bc_pt_t ) // '"'
       CALL message( 'check_parameters', 'PA0063', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( ANY( wall_heatflux /= 0.0_wp )  .AND.                        &
         surface_heatflux == 9999999.9_wp )  THEN
       message_string = 'wall_heatflux additionally requires ' //     &
                        'setting of surface_heatflux'
       CALL message( 'check_parameters', 'PA0443', 1, 2, 0, 6, 0 )
    ENDIF

!
!   This IF clause needs revision, got too complex!!
    IF ( surface_heatflux == 9999999.9_wp  )  THEN
       constant_heatflux = .FALSE.
       IF ( large_scale_forcing  .OR.  land_surface  .OR.  urban_surface )  THEN
          IF ( ibc_pt_b == 0 )  THEN
             constant_heatflux = .FALSE.
          ELSEIF ( ibc_pt_b == 1 )  THEN
             constant_heatflux = .TRUE.
             surface_heatflux = 0.0_wp
          ENDIF
       ENDIF
    ELSE
       constant_heatflux = .TRUE.
    ENDIF

    IF ( top_heatflux     == 9999999.9_wp )  constant_top_heatflux = .FALSE.

    IF ( neutral )  THEN

       IF ( surface_heatflux /= 0.0_wp  .AND.                                  &
            surface_heatflux /= 9999999.9_wp )  THEN
          message_string = 'heatflux must not be set for pure neutral flow'
          CALL message( 'check_parameters', 'PA0351', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( top_heatflux /= 0.0_wp  .AND.  top_heatflux /= 9999999.9_wp )      &
       THEN
          message_string = 'heatflux must not be set for pure neutral flow'
          CALL message( 'check_parameters', 'PA0351', 1, 2, 0, 6, 0 )
       ENDIF

    ENDIF

    IF ( top_momentumflux_u /= 9999999.9_wp  .AND.                             &
         top_momentumflux_v /= 9999999.9_wp )  THEN
       constant_top_momentumflux = .TRUE.
    ELSEIF (  .NOT. ( top_momentumflux_u == 9999999.9_wp  .AND.                &
           top_momentumflux_v == 9999999.9_wp ) )  THEN
       message_string = 'both, top_momentumflux_u AND top_momentumflux_v ' //  &
                        'must be set'
       CALL message( 'check_parameters', 'PA0064', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- A given surface temperature implies Dirichlet boundary condition for
!-- temperature. In this case specification of a constant heat flux is
!-- forbidden.
    IF ( ibc_pt_b == 0  .AND.  constant_heatflux  .AND.                        &
         surface_heatflux /= 0.0_wp )  THEN
       message_string = 'boundary_condition: bc_pt_b = "' // TRIM( bc_pt_b ) //&
                        '& is not allowed with constant_heatflux = .TRUE.'
       CALL message( 'check_parameters', 'PA0065', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( constant_heatflux  .AND.  pt_surface_initial_change /= 0.0_wp )  THEN
       WRITE ( message_string, * )  'constant_heatflux = .TRUE. is not allo',  &
               'wed with pt_surface_initial_change (/=0) = ',                  &
               pt_surface_initial_change
       CALL message( 'check_parameters', 'PA0066', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- A given temperature at the top implies Dirichlet boundary condition for
!-- temperature. In this case specification of a constant heat flux is
!-- forbidden.
    IF ( ibc_pt_t == 0  .AND.  constant_top_heatflux  .AND.                    &
         top_heatflux /= 0.0_wp )  THEN
       message_string = 'boundary_condition: bc_pt_t = "' // TRIM( bc_pt_t ) //&
                        '" is not allowed with constant_top_heatflux = .TRUE.'
       CALL message( 'check_parameters', 'PA0067', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Set boundary conditions for total water content
    IF ( humidity )  THEN

       IF ( ANY( wall_humidityflux /= 0.0_wp )  .AND.                        &
            surface_waterflux == 9999999.9_wp )  THEN
          message_string = 'wall_humidityflux additionally requires ' //     &
                           'setting of surface_waterflux'
          CALL message( 'check_parameters', 'PA0444', 1, 2, 0, 6, 0 )
       ENDIF

       CALL set_bc_scalars( 'q', bc_q_b, bc_q_t, ibc_q_b, ibc_q_t,           &
                            'PA0071', 'PA0072' )

       IF ( surface_waterflux == 9999999.9_wp  )  THEN
          constant_waterflux = .FALSE.
          IF ( large_scale_forcing .OR. land_surface )  THEN
             IF ( ibc_q_b == 0 )  THEN
                constant_waterflux = .FALSE.
             ELSEIF ( ibc_q_b == 1 )  THEN
                constant_waterflux = .TRUE.
             ENDIF
          ENDIF
       ELSE
          constant_waterflux = .TRUE.
       ENDIF

       CALL check_bc_scalars( 'q', bc_q_b, ibc_q_b, 'PA0073', 'PA0074',        &
                              constant_waterflux, q_surface_initial_change )

    ENDIF

    IF ( passive_scalar )  THEN

       IF ( ANY( wall_scalarflux /= 0.0_wp )  .AND.                            &
            surface_scalarflux == 9999999.9_wp )  THEN
          message_string = 'wall_scalarflux additionally requires ' //         &
                           'setting of surface_scalarflux'
          CALL message( 'check_parameters', 'PA0445', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( surface_scalarflux == 9999999.9_wp )  constant_scalarflux = .FALSE.

       CALL set_bc_scalars( 's', bc_s_b, bc_s_t, ibc_s_b, ibc_s_t,             &
                            'PA0071', 'PA0072' )

       CALL check_bc_scalars( 's', bc_s_b, ibc_s_b, 'PA0073', 'PA0074',        &
                              constant_scalarflux, s_surface_initial_change )

       IF ( top_scalarflux == 9999999.9_wp )  constant_top_scalarflux = .FALSE.
!
!--    A fixed scalar concentration at the top implies Dirichlet boundary
!--    condition for scalar. Hence, in this case specification of a constant
!--    scalar flux is forbidden.
       IF ( ( ibc_s_t == 0 .OR. ibc_s_t == 2 )  .AND.  constant_top_scalarflux &
               .AND.  top_scalarflux /= 0.0_wp )  THEN
          message_string = 'boundary condition: bc_s_t = "' //                 &
                           TRIM( bc_s_t ) // '" is not allowed with ' //       &
                           'top_scalarflux /= 0.0'
          CALL message( 'check_parameters', 'PA0441', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

!
!-- Boundary conditions for chemical species
    IF ( air_chemistry )  CALL chem_boundary_conds( 'init' )

!
!-- Boundary conditions for horizontal components of wind speed
    IF ( bc_uv_b == 'dirichlet' )  THEN
       ibc_uv_b = 0
    ELSEIF ( bc_uv_b == 'neumann' )  THEN
       ibc_uv_b = 1
       IF ( constant_flux_layer )  THEN
          message_string = 'boundary condition: bc_uv_b = "' //                &
               TRIM( bc_uv_b ) // '" is not allowed with constant_flux_layer'  &
               // ' = .TRUE.'
          CALL message( 'check_parameters', 'PA0075', 1, 2, 0, 6, 0 )
       ENDIF
    ELSE
       message_string = 'unknown boundary condition: bc_uv_b = "' //           &
                        TRIM( bc_uv_b ) // '"'
       CALL message( 'check_parameters', 'PA0076', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- In case of coupled simulations u and v at the ground in atmosphere will be
!-- assigned with the u and v values of the ocean surface
    IF ( coupling_mode == 'atmosphere_to_ocean' )  THEN
       ibc_uv_b = 2
    ENDIF

    IF ( coupling_mode == 'ocean_to_atmosphere' )  THEN
       bc_uv_t = 'neumann'
       ibc_uv_t = 1
    ELSE
       IF ( bc_uv_t == 'dirichlet' .OR. bc_uv_t == 'dirichlet_0' )  THEN
          ibc_uv_t = 0
          IF ( bc_uv_t == 'dirichlet_0' )  THEN
!
!--          Velocities for the initial u,v-profiles are set zero at the top
!--          in case of dirichlet_0 conditions
             u_init(nzt+1)    = 0.0_wp
             v_init(nzt+1)    = 0.0_wp
          ENDIF
       ELSEIF ( bc_uv_t == 'neumann' )  THEN
          ibc_uv_t = 1
       ELSEIF ( bc_uv_t == 'nested'  .OR.  bc_uv_t == 'nesting_offline' )  THEN
          ibc_uv_t = 3
       ELSE
          message_string = 'unknown boundary condition: bc_uv_t = "' //        &
                           TRIM( bc_uv_t ) // '"'
          CALL message( 'check_parameters', 'PA0077', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

!
!-- Compute and check, respectively, the Rayleigh Damping parameter
    IF ( rayleigh_damping_factor == -1.0_wp )  THEN
       rayleigh_damping_factor = 0.0_wp
    ELSE
       IF ( rayleigh_damping_factor < 0.0_wp  .OR.                             &
            rayleigh_damping_factor > 1.0_wp )  THEN
          WRITE( message_string, * )  'rayleigh_damping_factor = ',            &
                              rayleigh_damping_factor, ' out of range [0.0,1.0]'
          CALL message( 'check_parameters', 'PA0078', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

    IF ( rayleigh_damping_height == -1.0_wp )  THEN
       IF (  .NOT.  ocean_mode )  THEN
          rayleigh_damping_height = 0.66666666666_wp * zu(nzt)
       ELSE
          rayleigh_damping_height = 0.66666666666_wp * zu(nzb)
       ENDIF
    ELSE
       IF (  .NOT.  ocean_mode )  THEN
          IF ( rayleigh_damping_height < 0.0_wp  .OR.                          &
               rayleigh_damping_height > zu(nzt) )  THEN
             WRITE( message_string, * )  'rayleigh_damping_height = ',         &
                   rayleigh_damping_height, ' out of range [0.0,', zu(nzt), ']'
             CALL message( 'check_parameters', 'PA0079', 1, 2, 0, 6, 0 )
          ENDIF
       ELSE
          IF ( rayleigh_damping_height > 0.0_wp  .OR.                          &
               rayleigh_damping_height < zu(nzb) )  THEN
             WRITE( message_string, * )  'rayleigh_damping_height = ',         &
                   rayleigh_damping_height, ' out of range [0.0,', zu(nzb), ']'
             CALL message( 'check_parameters', 'PA0079', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF
    ENDIF

!
!-- Check number of chosen statistic regions
    IF ( statistic_regions < 0 )  THEN
       WRITE ( message_string, * ) 'number of statistic_regions = ',           &
                   statistic_regions+1, ' is not allowed'
       CALL message( 'check_parameters', 'PA0082', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( normalizing_region > statistic_regions  .OR.                          &
         normalizing_region < 0)  THEN
       WRITE ( message_string, * ) 'normalizing_region = ',                    &
                normalizing_region, ' must be >= 0 and <= ',statistic_regions, &
                ' (value of statistic_regions)'
       CALL message( 'check_parameters', 'PA0083', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Set the default intervals for data output, if necessary
!-- NOTE: dt_dosp has already been set in spectra_parin
    IF ( dt_data_output /= 9999999.9_wp )  THEN
       IF ( dt_dopr           == 9999999.9_wp )  dt_dopr           = dt_data_output
       IF ( dt_dopts          == 9999999.9_wp )  dt_dopts          = dt_data_output
       IF ( dt_do2d_xy        == 9999999.9_wp )  dt_do2d_xy        = dt_data_output
       IF ( dt_do2d_xz        == 9999999.9_wp )  dt_do2d_xz        = dt_data_output
       IF ( dt_do2d_yz        == 9999999.9_wp )  dt_do2d_yz        = dt_data_output
       IF ( dt_do3d           == 9999999.9_wp )  dt_do3d           = dt_data_output
       IF ( dt_data_output_av == 9999999.9_wp )  dt_data_output_av = dt_data_output
       DO  mid = 1, max_masks
          IF ( dt_domask(mid) == 9999999.9_wp )  dt_domask(mid)    = dt_data_output
       ENDDO
    ENDIF

!
!-- Set the default skip time intervals for data output, if necessary
    IF ( skip_time_dopr    == 9999999.9_wp )                                   &
                                       skip_time_dopr    = skip_time_data_output
    IF ( skip_time_do2d_xy == 9999999.9_wp )                                   &
                                       skip_time_do2d_xy = skip_time_data_output
    IF ( skip_time_do2d_xz == 9999999.9_wp )                                   &
                                       skip_time_do2d_xz = skip_time_data_output
    IF ( skip_time_do2d_yz == 9999999.9_wp )                                   &
                                       skip_time_do2d_yz = skip_time_data_output
    IF ( skip_time_do3d    == 9999999.9_wp )                                   &
                                       skip_time_do3d    = skip_time_data_output
    IF ( skip_time_data_output_av == 9999999.9_wp )                            &
                                skip_time_data_output_av = skip_time_data_output
    DO  mid = 1, max_masks
       IF ( skip_time_domask(mid) == 9999999.9_wp )                            &
                                skip_time_domask(mid)    = skip_time_data_output
    ENDDO

!
!-- Check the average intervals (first for 3d-data, then for profiles)
    IF ( averaging_interval > dt_data_output_av )  THEN
       WRITE( message_string, * )  'averaging_interval = ',                    &
             averaging_interval, ' must be <= dt_data_output_av = ',           &
             dt_data_output_av
       CALL message( 'check_parameters', 'PA0085', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( averaging_interval_pr == 9999999.9_wp )  THEN
       averaging_interval_pr = averaging_interval
    ENDIF

    IF ( averaging_interval_pr > dt_dopr )  THEN
       WRITE( message_string, * )  'averaging_interval_pr = ',                 &
             averaging_interval_pr, ' must be <= dt_dopr = ', dt_dopr
       CALL message( 'check_parameters', 'PA0086', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Set the default interval for profiles entering the temporal average
    IF ( dt_averaging_input_pr == 9999999.9_wp )  THEN
       dt_averaging_input_pr = dt_averaging_input
    ENDIF

!
!-- Set the default interval for the output of timeseries to a reasonable
!-- value (tries to minimize the number of calls of flow_statistics)
    IF ( dt_dots == 9999999.9_wp )  THEN
       IF ( averaging_interval_pr == 0.0_wp )  THEN
          dt_dots = MIN( dt_run_control, dt_dopr )
       ELSE
          dt_dots = MIN( dt_run_control, dt_averaging_input_pr )
       ENDIF
    ENDIF

!
!-- Check the sample rate for averaging (first for 3d-data, then for profiles)
    IF ( dt_averaging_input > averaging_interval )  THEN
       WRITE( message_string, * )  'dt_averaging_input = ',                    &
                dt_averaging_input, ' must be <= averaging_interval = ',       &
                averaging_interval
       CALL message( 'check_parameters', 'PA0088', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( dt_averaging_input_pr > averaging_interval_pr )  THEN
       WRITE( message_string, * )  'dt_averaging_input_pr = ',                 &
                dt_averaging_input_pr, ' must be <= averaging_interval_pr = ', &
                averaging_interval_pr
       CALL message( 'check_parameters', 'PA0089', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Determine the number of output profiles and check whether they are
!-- permissible
    DO  WHILE ( data_output_pr(dopr_n+1) /= '          ' )

       dopr_n = dopr_n + 1
       i = dopr_n

!
!--    Determine internal profile number (for hom, homs)
!--    and store height levels
       SELECT CASE ( TRIM( data_output_pr(i) ) )

          CASE ( 'u', '#u' )
             dopr_index(i) = 1
             dopr_unit(i)  = 'm/s'
             hom(:,2,1,:)  = SPREAD( zu, 2, statistic_regions+1 )
             IF ( data_output_pr(i)(1:1) == '#' )  THEN
                dopr_initial_index(i) = 5
                hom(:,2,5,:)          = SPREAD( zu, 2, statistic_regions+1 )
                data_output_pr(i)     = data_output_pr(i)(2:)
             ENDIF

          CASE ( 'v', '#v' )
             dopr_index(i) = 2
             dopr_unit(i)  = 'm/s'
             hom(:,2,2,:)  = SPREAD( zu, 2, statistic_regions+1 )
             IF ( data_output_pr(i)(1:1) == '#' )  THEN
                dopr_initial_index(i) = 6
                hom(:,2,6,:)          = SPREAD( zu, 2, statistic_regions+1 )
                data_output_pr(i)     = data_output_pr(i)(2:)
             ENDIF

          CASE ( 'w' )
             dopr_index(i) = 3
             dopr_unit(i)  = 'm/s'
             hom(:,2,3,:)  = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'theta', '#theta' )
             IF ( .NOT. bulk_cloud_model ) THEN
                dopr_index(i) = 4
                dopr_unit(i)  = 'K'
                hom(:,2,4,:)  = SPREAD( zu, 2, statistic_regions+1 )
                IF ( data_output_pr(i)(1:1) == '#' )  THEN
                   dopr_initial_index(i) = 7
                   hom(:,2,7,:)          = SPREAD( zu, 2, statistic_regions+1 )
                   hom(nzb,2,7,:)        = 0.0_wp    ! because zu(nzb) is negative
                   data_output_pr(i)     = data_output_pr(i)(2:)
                ENDIF
             ELSE
                dopr_index(i) = 43
                dopr_unit(i)  = 'K'
                hom(:,2,43,:)  = SPREAD( zu, 2, statistic_regions+1 )
                IF ( data_output_pr(i)(1:1) == '#' )  THEN
                   dopr_initial_index(i) = 28
                   hom(:,2,28,:)         = SPREAD( zu, 2, statistic_regions+1 )
                   hom(nzb,2,28,:)       = 0.0_wp    ! because zu(nzb) is negative
                   data_output_pr(i)     = data_output_pr(i)(2:)
                ENDIF
             ENDIF

          CASE ( 'e', '#e' )
             dopr_index(i)  = 8
             dopr_unit(i)   = 'm2/s2'
             hom(:,2,8,:)   = SPREAD( zu, 2, statistic_regions+1 )
             hom(nzb,2,8,:) = 0.0_wp
             IF ( data_output_pr(i)(1:1) == '#' )  THEN
                dopr_initial_index(i) = 8
                hom(:,2,8,:)          = SPREAD( zu, 2, statistic_regions+1 )
                data_output_pr(i)     = data_output_pr(i)(2:)
             ENDIF

          CASE ( 'km', '#km' )
             dopr_index(i)  = 9
             dopr_unit(i)   = 'm2/s'
             hom(:,2,9,:)   = SPREAD( zu, 2, statistic_regions+1 )
             hom(nzb,2,9,:) = 0.0_wp
             IF ( data_output_pr(i)(1:1) == '#' )  THEN
                dopr_initial_index(i) = 23
                hom(:,2,23,:)         = hom(:,2,9,:)
                data_output_pr(i)     = data_output_pr(i)(2:)
             ENDIF

          CASE ( 'kh', '#kh' )
             dopr_index(i)   = 10
             dopr_unit(i)    = 'm2/s'
             hom(:,2,10,:)   = SPREAD( zu, 2, statistic_regions+1 )
             hom(nzb,2,10,:) = 0.0_wp
             IF ( data_output_pr(i)(1:1) == '#' )  THEN
                dopr_initial_index(i) = 24
                hom(:,2,24,:)         = hom(:,2,10,:)
                data_output_pr(i)     = data_output_pr(i)(2:)
             ENDIF

          CASE ( 'l', '#l' )
             dopr_index(i)   = 11
             dopr_unit(i)    = 'm'
             hom(:,2,11,:)   = SPREAD( zu, 2, statistic_regions+1 )
             hom(nzb,2,11,:) = 0.0_wp
             IF ( data_output_pr(i)(1:1) == '#' )  THEN
                dopr_initial_index(i) = 25
                hom(:,2,25,:)         = hom(:,2,11,:)
                data_output_pr(i)     = data_output_pr(i)(2:)
             ENDIF

          CASE ( 'w"u"' )
             dopr_index(i) = 12
             dopr_unit(i)  = TRIM ( momentumflux_output_unit )
             hom(:,2,12,:) = SPREAD( zw, 2, statistic_regions+1 )
             IF ( constant_flux_layer )  hom(nzb,2,12,:) = zu(1)

          CASE ( 'w*u*' )
             dopr_index(i) = 13
             dopr_unit(i)  = TRIM ( momentumflux_output_unit )
             hom(:,2,13,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'w"v"' )
             dopr_index(i) = 14
             dopr_unit(i)  = TRIM ( momentumflux_output_unit )
             hom(:,2,14,:) = SPREAD( zw, 2, statistic_regions+1 )
             IF ( constant_flux_layer )  hom(nzb,2,14,:) = zu(1)

          CASE ( 'w*v*' )
             dopr_index(i) = 15
             dopr_unit(i)  = TRIM ( momentumflux_output_unit )
             hom(:,2,15,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'w"theta"' )
             dopr_index(i) = 16
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,16,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'w*theta*' )
             dopr_index(i) = 17
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,17,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'wtheta' )
             dopr_index(i) = 18
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,18,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'wu' )
             dopr_index(i) = 19
             dopr_unit(i)  = TRIM ( momentumflux_output_unit )
             hom(:,2,19,:) = SPREAD( zw, 2, statistic_regions+1 )
             IF ( constant_flux_layer )  hom(nzb,2,19,:) = zu(1)

          CASE ( 'wv' )
             dopr_index(i) = 20
             dopr_unit(i)  = TRIM ( momentumflux_output_unit )
             hom(:,2,20,:) = SPREAD( zw, 2, statistic_regions+1 )
             IF ( constant_flux_layer )  hom(nzb,2,20,:) = zu(1)

          CASE ( 'w*theta*BC' )
             dopr_index(i) = 21
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,21,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'wthetaBC' )
             dopr_index(i) = 22
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,22,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'u*2' )
             dopr_index(i) = 30
             dopr_unit(i)  = 'm2/s2'
             hom(:,2,30,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'v*2' )
             dopr_index(i) = 31
             dopr_unit(i)  = 'm2/s2'
             hom(:,2,31,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'w*2' )
             dopr_index(i) = 32
             dopr_unit(i)  = 'm2/s2'
             hom(:,2,32,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'theta*2' )
             dopr_index(i) = 33
             dopr_unit(i)  = 'K2'
             hom(:,2,33,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'e*' )
             dopr_index(i) = 34
             dopr_unit(i)  = 'm2/s2'
             hom(:,2,34,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'w*2theta*' )
             dopr_index(i) = 35
             dopr_unit(i)  = 'K m2/s2'
             hom(:,2,35,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'w*theta*2' )
             dopr_index(i) = 36
             dopr_unit(i)  = 'K2 m/s'
             hom(:,2,36,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'w*e*' )
             dopr_index(i) = 37
             dopr_unit(i)  = 'm3/s3'
             hom(:,2,37,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'w*3' )
             dopr_index(i) = 38
             dopr_unit(i)  = 'm3/s3'
             hom(:,2,38,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'Sw' )
             dopr_index(i) = 39
             dopr_unit(i)  = 'none'
             hom(:,2,39,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'p' )
             dopr_index(i) = 40
             dopr_unit(i)  = 'Pa'
             hom(:,2,40,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'q', '#q' )
             IF ( .NOT. humidity )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for humidity = .FALSE.'
                CALL message( 'check_parameters', 'PA0092', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 41
                dopr_unit(i)  = 'kg/kg'
                hom(:,2,41,:) = SPREAD( zu, 2, statistic_regions+1 )
                IF ( data_output_pr(i)(1:1) == '#' )  THEN
                   dopr_initial_index(i) = 26
                   hom(:,2,26,:)         = SPREAD( zu, 2, statistic_regions+1 )
                   hom(nzb,2,26,:)       = 0.0_wp    ! because zu(nzb) is negative
                   data_output_pr(i)     = data_output_pr(i)(2:)
                ENDIF
             ENDIF

          CASE ( 's', '#s' )
             IF ( .NOT. passive_scalar )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for passive_scalar = .FALSE.'
                CALL message( 'check_parameters', 'PA0093', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 115
                dopr_unit(i)  = 'kg/m3'
                hom(:,2,115,:) = SPREAD( zu, 2, statistic_regions+1 )
                IF ( data_output_pr(i)(1:1) == '#' )  THEN
                   dopr_initial_index(i) = 121
                   hom(:,2,121,:)        = SPREAD( zu, 2, statistic_regions+1 )
                   hom(nzb,2,121,:)      = 0.0_wp    ! because zu(nzb) is negative
                   data_output_pr(i)     = data_output_pr(i)(2:)
                ENDIF
             ENDIF

          CASE ( 'qv', '#qv' )
             IF ( .NOT. bulk_cloud_model ) THEN
                dopr_index(i) = 41
                dopr_unit(i)  = 'kg/kg'
                hom(:,2,41,:) = SPREAD( zu, 2, statistic_regions+1 )
                IF ( data_output_pr(i)(1:1) == '#' )  THEN
                   dopr_initial_index(i) = 26
                   hom(:,2,26,:)         = SPREAD( zu, 2, statistic_regions+1 )
                   hom(nzb,2,26,:)       = 0.0_wp    ! because zu(nzb) is negative
                   data_output_pr(i)     = data_output_pr(i)(2:)
                ENDIF
             ELSE
                dopr_index(i) = 42
                dopr_unit(i)  = 'kg/kg'
                hom(:,2,42,:) = SPREAD( zu, 2, statistic_regions+1 )
                IF ( data_output_pr(i)(1:1) == '#' )  THEN
                   dopr_initial_index(i) = 27
                   hom(:,2,27,:)         = SPREAD( zu, 2, statistic_regions+1 )
                   hom(nzb,2,27,:)       = 0.0_wp   ! because zu(nzb) is negative
                   data_output_pr(i)     = data_output_pr(i)(2:)
                ENDIF
             ENDIF

          CASE ( 'thetal', '#thetal' )
             IF ( .NOT. bulk_cloud_model ) THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for bulk_cloud_model = .FALSE.'
                CALL message( 'check_parameters', 'PA0094', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 4
                dopr_unit(i)  = 'K'
                hom(:,2,4,:)  = SPREAD( zu, 2, statistic_regions+1 )
                IF ( data_output_pr(i)(1:1) == '#' )  THEN
                   dopr_initial_index(i) = 7
                   hom(:,2,7,:)          = SPREAD( zu, 2, statistic_regions+1 )
                   hom(nzb,2,7,:)        = 0.0_wp    ! because zu(nzb) is negative
                   data_output_pr(i)     = data_output_pr(i)(2:)
                ENDIF
             ENDIF

          CASE ( 'thetav', '#thetav' )
             dopr_index(i) = 44
             dopr_unit(i)  = 'K'
             hom(:,2,44,:) = SPREAD( zu, 2, statistic_regions+1 )
             IF ( data_output_pr(i)(1:1) == '#' )  THEN
                dopr_initial_index(i) = 29
                hom(:,2,29,:)         = SPREAD( zu, 2, statistic_regions+1 )
                hom(nzb,2,29,:)       = 0.0_wp    ! because zu(nzb) is negative
                data_output_pr(i)     = data_output_pr(i)(2:)
             ENDIF

          CASE ( 'w"thetav"' )
             dopr_index(i) = 45
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,45,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'w*thetav*' )
             dopr_index(i) = 46
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,46,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'wthetav' )
             dopr_index(i) = 47
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,47,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'w"q"' )
             IF (  .NOT.  humidity )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for humidity = .FALSE.'
                CALL message( 'check_parameters', 'PA0092', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 48
                dopr_unit(i)  = TRIM ( waterflux_output_unit )
                hom(:,2,48,:) = SPREAD( zw, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'w*q*' )
             IF (  .NOT.  humidity )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for humidity = .FALSE.'
                CALL message( 'check_parameters', 'PA0092', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 49
                dopr_unit(i)  = TRIM ( waterflux_output_unit )
                hom(:,2,49,:) = SPREAD( zw, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'wq' )
             IF (  .NOT.  humidity )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for humidity = .FALSE.'
                CALL message( 'check_parameters', 'PA0092', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 50
                dopr_unit(i)  = TRIM ( waterflux_output_unit )
                hom(:,2,50,:) = SPREAD( zw, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'w"s"' )
             IF (  .NOT.  passive_scalar )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for passive_scalar = .FALSE.'
                CALL message( 'check_parameters', 'PA0093', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 117
                dopr_unit(i)  = 'kg/m3 m/s'
                hom(:,2,117,:) = SPREAD( zw, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'w*s*' )
             IF (  .NOT.  passive_scalar )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for passive_scalar = .FALSE.'
                CALL message( 'check_parameters', 'PA0093', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 114
                dopr_unit(i)  = 'kg/m3 m/s'
                hom(:,2,114,:) = SPREAD( zw, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'ws' )
             IF (  .NOT.  passive_scalar )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for passive_scalar = .FALSE.'
                CALL message( 'check_parameters', 'PA0093', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 118
                dopr_unit(i)  = 'kg/m3 m/s'
                hom(:,2,118,:) = SPREAD( zw, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'w"qv"' )
             IF ( humidity  .AND.  .NOT.  bulk_cloud_model )  THEN
                dopr_index(i) = 48
                dopr_unit(i)  = TRIM ( waterflux_output_unit )
                hom(:,2,48,:) = SPREAD( zw, 2, statistic_regions+1 )
             ELSEIF ( humidity  .AND.  bulk_cloud_model )  THEN
                dopr_index(i) = 51
                dopr_unit(i)  = TRIM ( waterflux_output_unit )
                hom(:,2,51,:) = SPREAD( zw, 2, statistic_regions+1 )
             ELSE
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for bulk_cloud_model = .FALSE. ' // &
                                 'and humidity = .FALSE.'
                CALL message( 'check_parameters', 'PA0095', 1, 2, 0, 6, 0 )
             ENDIF

          CASE ( 'w*qv*' )
             IF ( humidity  .AND.  .NOT. bulk_cloud_model )  THEN
                dopr_index(i) = 49
                dopr_unit(i)  = TRIM ( waterflux_output_unit )
                hom(:,2,49,:) = SPREAD( zw, 2, statistic_regions+1 )
             ELSEIF( humidity .AND. bulk_cloud_model ) THEN
                dopr_index(i) = 52
                dopr_unit(i)  = TRIM ( waterflux_output_unit )
                hom(:,2,52,:) = SPREAD( zw, 2, statistic_regions+1 )
             ELSE
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for bulk_cloud_model = .FALSE. ' // &
                                 'and humidity = .FALSE.'
                CALL message( 'check_parameters', 'PA0095', 1, 2, 0, 6, 0 )
             ENDIF

          CASE ( 'wqv' )
             IF ( humidity  .AND.  .NOT.  bulk_cloud_model )  THEN
                dopr_index(i) = 50
                dopr_unit(i)  = TRIM ( waterflux_output_unit )
                hom(:,2,50,:) = SPREAD( zw, 2, statistic_regions+1 )
             ELSEIF ( humidity  .AND.  bulk_cloud_model )  THEN
                dopr_index(i) = 53
                dopr_unit(i)  = TRIM ( waterflux_output_unit )
                hom(:,2,53,:) = SPREAD( zw, 2, statistic_regions+1 )
             ELSE
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for bulk_cloud_model = .FALSE. ' // &
                                 'and humidity = .FALSE.'
                CALL message( 'check_parameters', 'PA0095', 1, 2, 0, 6, 0 )
             ENDIF

          CASE ( 'ql' )
             IF (  .NOT.  bulk_cloud_model  .AND.  .NOT.  cloud_droplets )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for bulk_cloud_model = .FALSE. ' // &
                                 'and cloud_droplets = .FALSE.'
                CALL message( 'check_parameters', 'PA0096', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 54
                dopr_unit(i)  = 'kg/kg'
                hom(:,2,54,:)  = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'w*u*u*:dz' )
             dopr_index(i) = 55
             dopr_unit(i)  = 'm2/s3'
             hom(:,2,55,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'w*p*:dz' )
             dopr_index(i) = 56
             dopr_unit(i)  = 'm2/s3'
             hom(:,2,56,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'w"e:dz' )
             dopr_index(i) = 57
             dopr_unit(i)  = 'm2/s3'
             hom(:,2,57,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'u"theta"' )
             dopr_index(i) = 58
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,58,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'u*theta*' )
             dopr_index(i) = 59
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,59,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'utheta_t' )
             dopr_index(i) = 60
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,60,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'v"theta"' )
             dopr_index(i) = 61
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,61,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'v*theta*' )
             dopr_index(i) = 62
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,62,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'vtheta_t' )
             dopr_index(i) = 63
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,63,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'w*p*' )
             dopr_index(i) = 68
             dopr_unit(i)  = 'm3/s3'
             hom(:,2,68,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'w"e' )
             dopr_index(i) = 69
             dopr_unit(i)  = 'm3/s3'
             hom(:,2,69,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'q*2' )
             IF (  .NOT.  humidity )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for humidity = .FALSE.'
                CALL message( 'check_parameters', 'PA0092', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 70
                dopr_unit(i)  = 'kg2/kg2'
                hom(:,2,70,:) = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'hyp' )
             dopr_index(i) = 72
             dopr_unit(i)  = 'hPa'
             hom(:,2,72,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'rho' )
             dopr_index(i)  = 119
             dopr_unit(i)   = 'kg/m3'
             hom(:,2,119,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'rho_zw' )
             dopr_index(i)  = 120
             dopr_unit(i)   = 'kg/m3'
             hom(:,2,120,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'ug' )
             dopr_index(i) = 78
             dopr_unit(i)  = 'm/s'
             hom(:,2,78,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'vg' )
             dopr_index(i) = 79
             dopr_unit(i)  = 'm/s'
             hom(:,2,79,:) = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'w_subs' )
             IF (  .NOT.  large_scale_subsidence )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for large_scale_subsidence = .FALSE.'
                CALL message( 'check_parameters', 'PA0382', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 80
                dopr_unit(i)  = 'm/s'
                hom(:,2,80,:) = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF

          CASE ( 's*2' )
             IF (  .NOT.  passive_scalar )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for passive_scalar = .FALSE.'
                CALL message( 'check_parameters', 'PA0185', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 116
                dopr_unit(i)  = 'kg2/m6'
                hom(:,2,116,:) = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF

          CASE DEFAULT
             unit = 'illegal'
!
!--          Check for other modules
             CALL module_interface_check_data_output_pr( data_output_pr(i), i, &
                                                         unit, dopr_unit(i) )

!
!--          No valid quantity found
             IF ( unit == 'illegal' )  THEN
                IF ( data_output_pr_user(1) /= ' ' )  THEN
                   message_string = 'illegal value for data_output_pr or ' //  &
                                    'data_output_pr_user = "' //               &
                                    TRIM( data_output_pr(i) ) // '"'
                   CALL message( 'check_parameters', 'PA0097', 1, 2, 0, 6, 0 )
                ELSE
                   message_string = 'illegal value for data_output_pr = "' //  &
                                    TRIM( data_output_pr(i) ) // '"'
                   CALL message( 'check_parameters', 'PA0098', 1, 2, 0, 6, 0 )
                ENDIF
             ENDIF

       END SELECT

    ENDDO


!
!-- Append user-defined data output variables to the standard data output
    IF ( data_output_user(1) /= ' ' )  THEN
       i = 1
       DO  WHILE ( data_output(i) /= ' '  .AND.  i <= 500 )
          i = i + 1
       ENDDO
       j = 1
       DO  WHILE ( data_output_user(j) /= ' '  .AND.  j <= 500 )
          IF ( i > 500 )  THEN
             message_string = 'number of output quantitities given by data' // &
                '_output and data_output_user exceeds the limit of 500'
             CALL message( 'check_parameters', 'PA0102', 1, 2, 0, 6, 0 )
          ENDIF
          data_output(i) = data_output_user(j)
          i = i + 1
          j = j + 1
       ENDDO
    ENDIF

!
!-- Check and set steering parameters for 2d/3d data output and averaging
    i   = 1
    DO  WHILE ( data_output(i) /= ' '  .AND.  i <= 500 )
!
!--    Check for data averaging
       ilen = LEN_TRIM( data_output(i) )
       j = 0                                                 ! no data averaging
       IF ( ilen > 3 )  THEN
          IF ( data_output(i)(ilen-2:ilen) == '_av' )  THEN
             j = 1                                           ! data averaging
             data_output(i) = data_output(i)(1:ilen-3)
          ENDIF
       ENDIF
!
!--    Check for cross section or volume data
       ilen = LEN_TRIM( data_output(i) )
       k = 0                                                   ! 3d data
       var = data_output(i)(1:ilen)
       IF ( ilen > 3 )  THEN
          IF ( data_output(i)(ilen-2:ilen) == '_xy'  .OR.                      &
               data_output(i)(ilen-2:ilen) == '_xz'  .OR.                      &
               data_output(i)(ilen-2:ilen) == '_yz' )  THEN
             k = 1                                             ! 2d data
             var = data_output(i)(1:ilen-3)
          ENDIF
       ENDIF

!
!--    Check for allowed value and set units
       SELECT CASE ( TRIM( var ) )

          CASE ( 'e' )
             IF ( constant_diffusion )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                                 'res constant_diffusion = .FALSE.'
                CALL message( 'check_parameters', 'PA0103', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'm2/s2'

          CASE ( 'thetal' )
             IF (  .NOT.  bulk_cloud_model )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                         'res bulk_cloud_model = .TRUE.'
                CALL message( 'check_parameters', 'PA0108', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'K'

          CASE ( 'pc', 'pr' )
             IF (  .NOT.  particle_advection )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requir' // &
                   'es a "particle_parameters"-NAMELIST in the parameter ' //  &
                   'file (PARIN)'
                CALL message( 'check_parameters', 'PA0104', 1, 2, 0, 6, 0 )
             ENDIF
             IF ( TRIM( var ) == 'pc' )  unit = 'number'
             IF ( TRIM( var ) == 'pr' )  unit = 'm'

          CASE ( 'q', 'thetav' )
             IF (  .NOT.  humidity )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                                 'res humidity = .TRUE.'
                CALL message( 'check_parameters', 'PA0105', 1, 2, 0, 6, 0 )
             ENDIF
             IF ( TRIM( var ) == 'q'   )  unit = 'kg/kg'
             IF ( TRIM( var ) == 'thetav' )  unit = 'K'

          CASE ( 'ql' )
             IF ( .NOT.  ( bulk_cloud_model  .OR.  cloud_droplets ) )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                      'res bulk_cloud_model = .TRUE. or cloud_droplets = .TRUE.'
                CALL message( 'check_parameters', 'PA0106', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'kg/kg'

          CASE ( 'ql_c', 'ql_v', 'ql_vp' )
             IF (  .NOT.  cloud_droplets )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                                 'res cloud_droplets = .TRUE.'
                CALL message( 'check_parameters', 'PA0107', 1, 2, 0, 6, 0 )
             ENDIF
             IF ( TRIM( var ) == 'ql_c'  )  unit = 'kg/kg'
             IF ( TRIM( var ) == 'ql_v'  )  unit = 'm3'
             IF ( TRIM( var ) == 'ql_vp' )  unit = 'none'

          CASE ( 'qv' )
             IF (  .NOT.  bulk_cloud_model )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                                 'res bulk_cloud_model = .TRUE.'
                CALL message( 'check_parameters', 'PA0108', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'kg/kg'

          CASE ( 's' )
             IF (  .NOT.  passive_scalar )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                                 'res passive_scalar = .TRUE.'
                CALL message( 'check_parameters', 'PA0110', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'kg/m3'

          CASE ( 'p', 'theta', 'u', 'v', 'w' )
             IF ( TRIM( var ) == 'p'  )  unit = 'Pa'
             IF ( TRIM( var ) == 'theta' )  unit = 'K'
             IF ( TRIM( var ) == 'u'  )  unit = 'm/s'
             IF ( TRIM( var ) == 'v'  )  unit = 'm/s'
             IF ( TRIM( var ) == 'w'  )  unit = 'm/s'
             CONTINUE

          CASE ( 'ghf*', 'lwp*', 'ol*', 'qsws*', 'r_a*',                       &
                 'shf*', 'ssws*', 't*', 'tsurf*', 'us*',                       &
                 'z0*', 'z0h*', 'z0q*' )
             IF ( k == 0  .OR.  data_output(i)(ilen-2:ilen) /= '_xy' )  THEN
                message_string = 'illegal value for data_output: "' //         &
                                 TRIM( var ) // '" & only 2d-horizontal ' //   &
                                 'cross sections are allowed for this value'
                CALL message( 'check_parameters', 'PA0111', 1, 2, 0, 6, 0 )
             ENDIF

             IF ( TRIM( var ) == 'lwp*'  .AND.  .NOT. bulk_cloud_model )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                                 'res bulk_cloud_model = .TRUE.'
                CALL message( 'check_parameters', 'PA0108', 1, 2, 0, 6, 0 )
             ENDIF
             IF ( TRIM( var ) == 'qsws*'  .AND.  .NOT.  humidity )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                                 'res humidity = .TRUE.'
                CALL message( 'check_parameters', 'PA0322', 1, 2, 0, 6, 0 )
             ENDIF

             IF ( TRIM( var ) == 'ghf*'  .AND.  .NOT.  land_surface )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                                 'res land_surface = .TRUE.'
                CALL message( 'check_parameters', 'PA0404', 1, 2, 0, 6, 0 )
             ENDIF

             IF ( ( TRIM( var ) == 'r_a*' .OR.  TRIM( var ) == 'ghf*' )        &
                 .AND.  .NOT.  land_surface  .AND.  .NOT.  urban_surface )     &         
             THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                                 'res land_surface = .TRUE. or ' //            &
                                 'urban_surface = .TRUE.'
                CALL message( 'check_parameters', 'PA0404', 1, 2, 0, 6, 0 )
             ENDIF
             
             IF ( TRIM( var ) == 'ssws*'  .AND.  .NOT.  passive_scalar )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                                 'res passive_scalar = .TRUE.'
                CALL message( 'check_parameters', 'PA0361', 1, 2, 0, 6, 0 )
             ENDIF

             IF ( TRIM( var ) == 'ghf*'   )  unit = 'W/m2'
             IF ( TRIM( var ) == 'lwp*'   )  unit = 'kg/m2'
             IF ( TRIM( var ) == 'ol*'    )  unit = 'm'
             IF ( TRIM( var ) == 'qsws*'  )  unit = 'kgm/kgs'
             IF ( TRIM( var ) == 'r_a*'   )  unit = 's/m'
             IF ( TRIM( var ) == 'shf*'   )  unit = 'K*m/s'
             IF ( TRIM( var ) == 'ssws*'  )  unit = 'kg/m2*s'
             IF ( TRIM( var ) == 't*'     )  unit = 'K'
             IF ( TRIM( var ) == 'tsurf*' )  unit = 'K'
             IF ( TRIM( var ) == 'us*'    )  unit = 'm/s'
             IF ( TRIM( var ) == 'z0*'    )  unit = 'm'
             IF ( TRIM( var ) == 'z0h*'   )  unit = 'm'
!
!--          Output of surface latent and sensible heat flux will be in W/m2 
!--          in case of natural- and urban-type surfaces, even if 
!--          flux_output_mode is set to kinematic units.
             IF ( land_surface  .OR.  urban_surface )  THEN
                IF ( TRIM( var ) == 'shf*'   )  unit = 'W/m2'
                IF ( TRIM( var ) == 'qsws*'  )  unit = 'W/m2'
             ENDIF

          CASE DEFAULT
!
!--          Check for other modules
             CALL module_interface_check_data_output( var, unit, i, j, ilen, k )

             IF ( unit == 'illegal' )  THEN
                IF ( data_output_user(1) /= ' ' )  THEN
                   message_string = 'illegal value for data_output or ' //     &
                         'data_output_user = "' // TRIM( data_output(i) ) // '"'
                   CALL message( 'check_parameters', 'PA0114', 1, 2, 0, 6, 0 )
                ELSE
                   message_string = 'illegal value for data_output = "' //     &
                                    TRIM( data_output(i) ) // '"'
                   CALL message( 'check_parameters', 'PA0115', 1, 2, 0, 6, 0 )
                ENDIF
             ENDIF

       END SELECT
!
!--    Set the internal steering parameters appropriately
       IF ( k == 0 )  THEN
          do3d_no(j)              = do3d_no(j) + 1
          do3d(j,do3d_no(j))      = data_output(i)
          do3d_unit(j,do3d_no(j)) = unit
       ELSE
          do2d_no(j)              = do2d_no(j) + 1
          do2d(j,do2d_no(j))      = data_output(i)
          do2d_unit(j,do2d_no(j)) = unit
          IF ( data_output(i)(ilen-2:ilen) == '_xy' )  THEN
             data_output_xy(j) = .TRUE.
          ENDIF
          IF ( data_output(i)(ilen-2:ilen) == '_xz' )  THEN
             data_output_xz(j) = .TRUE.
          ENDIF
          IF ( data_output(i)(ilen-2:ilen) == '_yz' )  THEN
             data_output_yz(j) = .TRUE.
          ENDIF
       ENDIF

       IF ( j == 1 )  THEN
!
!--       Check, if variable is already subject to averaging
          found = .FALSE.
          DO  k = 1, doav_n
             IF ( TRIM( doav(k) ) == TRIM( var ) )  found = .TRUE.
          ENDDO

          IF ( .NOT. found )  THEN
             doav_n = doav_n + 1
             doav(doav_n) = var
          ENDIF
       ENDIF

       i = i + 1
    ENDDO

!
!-- Averaged 2d or 3d output requires that an averaging interval has been set
    IF ( doav_n > 0  .AND.  averaging_interval == 0.0_wp )  THEN
       WRITE( message_string, * )  'output of averaged quantity "',            &
                                   TRIM( doav(1) ), '_av" requires to set a ', &
                                   'non-zero averaging interval'
       CALL message( 'check_parameters', 'PA0323', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Check sectional planes and store them in one shared array
    IF ( ANY( section_xy > nz + 1 ) )  THEN
       WRITE( message_string, * )  'section_xy must be <= nz + 1 = ', nz + 1
       CALL message( 'check_parameters', 'PA0319', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( ANY( section_xz > ny + 1 ) )  THEN
       WRITE( message_string, * )  'section_xz must be <= ny + 1 = ', ny + 1
       CALL message( 'check_parameters', 'PA0320', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( ANY( section_yz > nx + 1 ) )  THEN
       WRITE( message_string, * )  'section_yz must be <= nx + 1 = ', nx + 1
       CALL message( 'check_parameters', 'PA0321', 1, 2, 0, 6, 0 )
    ENDIF
    section(:,1) = section_xy
    section(:,2) = section_xz
    section(:,3) = section_yz

    IF ( ANY( data_output_xy ) .AND. .NOT. ANY( section(:,1) /= -9999 ) )  THEN
       WRITE( message_string, * )  'section_xy not defined for requested '  // &
                                   'xy-cross section output.&At least one ' // &
                                   'cross section must be given.'
       CALL message( 'check_parameters', 'PA0681', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( ANY( data_output_xz ) .AND. .NOT. ANY( section(:,2) /= -9999 ) )  THEN
       WRITE( message_string, * )  'section_xz not defined for requested '  // &
                                   'xz-cross section output.&At least one ' // &
                                   'cross section must be given.'
       CALL message( 'check_parameters', 'PA0681', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( ANY( data_output_yz ) .AND. .NOT. ANY( section(:,3) /= -9999 ) )  THEN
       WRITE( message_string, * )  'section_yz not defined for requested '  // &
                                   'yz-cross section output.&At least one ' // &
                                   'cross section must be given.'
       CALL message( 'check_parameters', 'PA0681', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Upper plot limit for 3D arrays
    IF ( nz_do3d == -9999 )  nz_do3d = nzt + 1

!
!-- Set output format string (used in header)
    SELECT CASE ( netcdf_data_format )
       CASE ( 1 )
          netcdf_data_format_string = 'netCDF classic'
       CASE ( 2 )
          netcdf_data_format_string = 'netCDF 64bit offset'
       CASE ( 3 )
          netcdf_data_format_string = 'netCDF4/HDF5'
       CASE ( 4 )
          netcdf_data_format_string = 'netCDF4/HDF5 classic'
       CASE ( 5 )
          netcdf_data_format_string = 'parallel netCDF4/HDF5'
       CASE ( 6 )
          netcdf_data_format_string = 'parallel netCDF4/HDF5 classic'

    END SELECT

!
!-- Check mask conditions
    DO mid = 1, max_masks
       IF ( data_output_masks(mid,1) /= ' '  .OR.                              &
            data_output_masks_user(mid,1) /= ' ' ) THEN
          masks = masks + 1
       ENDIF
    ENDDO

    IF ( masks < 0  .OR.  masks > max_masks )  THEN
       WRITE( message_string, * )  'illegal value: masks must be >= 0 and ',   &
            '<= ', max_masks, ' (=max_masks)'
       CALL message( 'check_parameters', 'PA0325', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( masks > 0 )  THEN
       mask_scale(1) = mask_scale_x
       mask_scale(2) = mask_scale_y
       mask_scale(3) = mask_scale_z
       IF ( ANY( mask_scale <= 0.0_wp ) )  THEN
          WRITE( message_string, * )                                           &
               'illegal value: mask_scale_x, mask_scale_y and mask_scale_z',   &
               'must be > 0.0'
          CALL message( 'check_parameters', 'PA0326', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Generate masks for masked data output
!--    Parallel netcdf output is not tested so far for masked data, hence
!--    netcdf_data_format is switched back to non-parallel output.
       netcdf_data_format_save = netcdf_data_format
       IF ( netcdf_data_format > 4 )  THEN
          IF ( netcdf_data_format == 5 ) netcdf_data_format = 3
          IF ( netcdf_data_format == 6 ) netcdf_data_format = 4
          message_string = 'netCDF file formats '//                            &
                           '5 (parallel netCDF 4) and ' //                     &
                           '6 (parallel netCDF 4 Classic model) '//            &
                           '& are currently not supported (not yet tested) ' //&
                           'for masked data. &Using respective non-parallel' //&
                           ' output for masked data.'
          CALL message( 'check_parameters', 'PA0383', 0, 0, 0, 6, 0 )
       ENDIF
       CALL init_masks
       netcdf_data_format = netcdf_data_format_save
    ENDIF

!
!-- Check the NetCDF data format
    IF ( netcdf_data_format > 2 )  THEN
#if defined( __netcdf4 )
       CONTINUE
#else
       message_string = 'netCDF: netCDF4 format requested but no ' //          &
                        'cpp-directive __netcdf4 given & switch '  //          &
                        'back to 64-bit offset format'
       CALL message( 'check_parameters', 'PA0171', 0, 1, 0, 6, 0 )
       netcdf_data_format = 2
#endif
    ENDIF
    IF ( netcdf_data_format > 4 )  THEN
#if defined( __netcdf4 ) && defined( __netcdf4_parallel )
       CONTINUE
#else
       message_string = 'netCDF: netCDF4 parallel output requested but no ' // &
                        'cpp-directive __netcdf4_parallel given & switch '   //&
                        'back to netCDF4 non-parallel output'
       CALL message( 'check_parameters', 'PA0099', 0, 1, 0, 6, 0 )
       netcdf_data_format = netcdf_data_format - 2
#endif
    ENDIF

!
!-- Calculate fixed number of output time levels for parallel netcdf output.
!-- The time dimension has to be defined as limited for parallel output,
!-- because otherwise the I/O performance drops significantly.
    IF ( netcdf_data_format > 4 )  THEN

!
!--    Check if any of the follwoing data output interval is 0.0s, which is
!--    not allowed for parallel output.
       CALL check_dt_do( dt_do3d,           'dt_do3d'           )
       CALL check_dt_do( dt_do2d_xy,        'dt_do2d_xy'        )
       CALL check_dt_do( dt_do2d_xz,        'dt_do2d_xz'        )
       CALL check_dt_do( dt_do2d_yz,        'dt_do2d_yz'        )
       CALL check_dt_do( dt_data_output_av, 'dt_data_output_av' )

!--    Set needed time levels (ntdim) to 
!--    saved time levels + to be saved time levels.
       ntdim_3d(0) = do3d_time_count(0) + CEILING(                                    &
                     ( end_time - MAX(                                                &
                         MERGE(skip_time_do3d, skip_time_do3d + spinup_time,          &
                               data_output_during_spinup ),                           &
                         simulated_time_at_begin )                                    &
                     ) / dt_do3d )
       IF ( do3d_at_begin ) ntdim_3d(0) = ntdim_3d(0) + 1

       ntdim_3d(1) = do3d_time_count(1) + CEILING(                                    &
                     ( end_time - MAX(                                                &
                         MERGE(   skip_time_data_output_av, skip_time_data_output_av  &
                                + spinup_time, data_output_during_spinup ),           &
                         simulated_time_at_begin )                                    &
                     ) / dt_data_output_av )

       ntdim_2d_xy(0) = do2d_xy_time_count(0) + CEILING(                              &
                        ( end_time - MAX(                                             &
                           MERGE(skip_time_do2d_xy, skip_time_do2d_xy + spinup_time,  &
                                 data_output_during_spinup ),                         &
                           simulated_time_at_begin )                                  &
                        ) / dt_do2d_xy )

       ntdim_2d_xz(0) = do2d_xz_time_count(0) + CEILING(                              &
                        ( end_time - MAX(                                             &
                         MERGE(skip_time_do2d_xz, skip_time_do2d_xz + spinup_time,    &
                               data_output_during_spinup ),                           &
                         simulated_time_at_begin )                                    &
                        ) / dt_do2d_xz )

       ntdim_2d_yz(0) = do2d_yz_time_count(0) + CEILING(                              &
                        ( end_time - MAX(                                             &
                         MERGE(skip_time_do2d_yz, skip_time_do2d_yz + spinup_time,    &
                               data_output_during_spinup ),                           &
                         simulated_time_at_begin )                                    &
                        ) / dt_do2d_yz )

       IF ( do2d_at_begin )  THEN
          ntdim_2d_xy(0) = ntdim_2d_xy(0) + 1
          ntdim_2d_xz(0) = ntdim_2d_xz(0) + 1
          ntdim_2d_yz(0) = ntdim_2d_yz(0) + 1
       ENDIF
!
!--    Please note, for averaged 2D data skip_time_data_output_av is the relavant 
!--    output control parameter. 
       ntdim_2d_xy(1) = do2d_xy_time_count(1) + CEILING(                              &
                     ( end_time - MAX( MERGE( skip_time_data_output_av,               &
                                              skip_time_data_output_av + spinup_time, &
                                              data_output_during_spinup ),            &
                                       simulated_time_at_begin )                      &
                     ) / dt_data_output_av )

       ntdim_2d_xz(1) = do2d_xz_time_count(1) + CEILING(                              &
                     ( end_time - MAX( MERGE( skip_time_data_output_av,               &
                                              skip_time_data_output_av + spinup_time, &
                                              data_output_during_spinup ),            &
                                       simulated_time_at_begin )                      &
                     ) / dt_data_output_av )

       ntdim_2d_yz(1) = do2d_yz_time_count(1) + CEILING(                              &
                     ( end_time - MAX( MERGE( skip_time_data_output_av,               &
                                              skip_time_data_output_av + spinup_time, &
                                              data_output_during_spinup ),            &
                                       simulated_time_at_begin )                      &
                     ) / dt_data_output_av )

    ENDIF

!
!-- Check, whether a constant diffusion coefficient shall be used
    IF ( km_constant /= -1.0_wp )  THEN
       IF ( km_constant < 0.0_wp )  THEN
          WRITE( message_string, * )  'km_constant = ', km_constant, ' < 0.0'
          CALL message( 'check_parameters', 'PA0121', 1, 2, 0, 6, 0 )
       ELSE
          IF ( prandtl_number < 0.0_wp )  THEN
             WRITE( message_string, * )  'prandtl_number = ', prandtl_number,  &
                                         ' < 0.0'
             CALL message( 'check_parameters', 'PA0122', 1, 2, 0, 6, 0 )
          ENDIF
          constant_diffusion = .TRUE.

          IF ( constant_flux_layer )  THEN
             message_string = 'constant_flux_layer is not allowed with fixed ' &
                              // 'value of km'
             CALL message( 'check_parameters', 'PA0123', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF
    ENDIF

!
!-- In case of non-cyclic lateral boundaries and a damping layer for the
!-- potential temperature, check the width of the damping layer
    IF ( bc_lr /= 'cyclic' ) THEN
       IF ( pt_damping_width < 0.0_wp  .OR.                                    &
            pt_damping_width > REAL( (nx+1) * dx ) )  THEN
          message_string = 'pt_damping_width out of range'
          CALL message( 'check_parameters', 'PA0124', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

    IF ( bc_ns /= 'cyclic' )  THEN
       IF ( pt_damping_width < 0.0_wp  .OR.                                    &
            pt_damping_width > REAL( (ny+1) * dy ) )  THEN
          message_string = 'pt_damping_width out of range'
          CALL message( 'check_parameters', 'PA0124', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

!
!-- Check value range for zeta = z/L
    IF ( zeta_min >= zeta_max )  THEN
       WRITE( message_string, * )  'zeta_min = ', zeta_min, ' must be less ',  &
                                   'than zeta_max = ', zeta_max
       CALL message( 'check_parameters', 'PA0125', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Check random generator
    IF ( (random_generator /= 'system-specific'      .AND.                     &
          random_generator /= 'random-parallel'   )  .AND.                     &
          random_generator /= 'numerical-recipes' )  THEN
       message_string = 'unknown random generator: random_generator = "' //    &
                        TRIM( random_generator ) // '"'
       CALL message( 'check_parameters', 'PA0135', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Determine upper and lower hight level indices for random perturbations
    IF ( disturbance_level_b == -9999999.9_wp )  THEN
       IF ( ocean_mode )  THEN
          disturbance_level_b     = zu((nzt*2)/3)
          disturbance_level_ind_b = ( nzt * 2 ) / 3
       ELSE
          disturbance_level_b     = zu(nzb+3)
          disturbance_level_ind_b = nzb + 3
       ENDIF
    ELSEIF ( disturbance_level_b < zu(3) )  THEN
       WRITE( message_string, * )  'disturbance_level_b = ',                   &
                           disturbance_level_b, ' must be >= ', zu(3), '(zu(3))'
       CALL message( 'check_parameters', 'PA0126', 1, 2, 0, 6, 0 )
    ELSEIF ( disturbance_level_b > zu(nzt-2) )  THEN
       WRITE( message_string, * )  'disturbance_level_b = ',                   &
                   disturbance_level_b, ' must be <= ', zu(nzt-2), '(zu(nzt-2))'
       CALL message( 'check_parameters', 'PA0127', 1, 2, 0, 6, 0 )
    ELSE
       DO  k = 3, nzt-2
          IF ( disturbance_level_b <= zu(k) )  THEN
             disturbance_level_ind_b = k
             EXIT
          ENDIF
       ENDDO
    ENDIF

    IF ( disturbance_level_t == -9999999.9_wp )  THEN
       IF ( ocean_mode )  THEN
          disturbance_level_t     = zu(nzt-3)
          disturbance_level_ind_t = nzt - 3
       ELSE
          disturbance_level_t     = zu(nzt/3)
          disturbance_level_ind_t = nzt / 3
       ENDIF
    ELSEIF ( disturbance_level_t > zu(nzt-2) )  THEN
       WRITE( message_string, * )  'disturbance_level_t = ',                   &
                   disturbance_level_t, ' must be <= ', zu(nzt-2), '(zu(nzt-2))'
       CALL message( 'check_parameters', 'PA0128', 1, 2, 0, 6, 0 )
    ELSEIF ( disturbance_level_t < disturbance_level_b )  THEN
       WRITE( message_string, * )  'disturbance_level_t = ',                   &
                   disturbance_level_t, ' must be >= disturbance_level_b = ',  &
                   disturbance_level_b
       CALL message( 'check_parameters', 'PA0129', 1, 2, 0, 6, 0 )
    ELSE
       DO  k = 3, nzt-2
          IF ( disturbance_level_t <= zu(k) )  THEN
             disturbance_level_ind_t = k
             EXIT
          ENDIF
       ENDDO
    ENDIF

!
!-- Check again whether the levels determined this way are ok.
!-- Error may occur at automatic determination and too few grid points in
!-- z-direction.
    IF ( disturbance_level_ind_t < disturbance_level_ind_b )  THEN
       WRITE( message_string, * )  'disturbance_level_ind_t = ',               &
                disturbance_level_ind_t, ' must be >= ',                       &
                'disturbance_level_ind_b = ', disturbance_level_ind_b
       CALL message( 'check_parameters', 'PA0130', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Determine the horizontal index range for random perturbations.
!-- In case of non-cyclic horizontal boundaries, no perturbations are imposed
!-- near the inflow and the perturbation area is further limited to ...(1)
!-- after the initial phase of the flow.

    IF ( bc_lr /= 'cyclic' )  THEN
       IF ( inflow_disturbance_begin == -1 )  THEN
          inflow_disturbance_begin = MIN( 10, nx/2 )
       ENDIF
       IF ( inflow_disturbance_begin < 0  .OR.  inflow_disturbance_begin > nx )&
       THEN
          message_string = 'inflow_disturbance_begin out of range'
          CALL message( 'check_parameters', 'PA0131', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( inflow_disturbance_end == -1 )  THEN
          inflow_disturbance_end = MIN( 100, 3*nx/4 )
       ENDIF
       IF ( inflow_disturbance_end < 0  .OR.  inflow_disturbance_end > nx )    &
       THEN
          message_string = 'inflow_disturbance_end out of range'
          CALL message( 'check_parameters', 'PA0132', 1, 2, 0, 6, 0 )
       ENDIF
    ELSEIF ( bc_ns /= 'cyclic' )  THEN
       IF ( inflow_disturbance_begin == -1 )  THEN
          inflow_disturbance_begin = MIN( 10, ny/2 )
       ENDIF
       IF ( inflow_disturbance_begin < 0  .OR.  inflow_disturbance_begin > ny )&
       THEN
          message_string = 'inflow_disturbance_begin out of range'
          CALL message( 'check_parameters', 'PA0131', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( inflow_disturbance_end == -1 )  THEN
          inflow_disturbance_end = MIN( 100, 3*ny/4 )
       ENDIF
       IF ( inflow_disturbance_end < 0  .OR.  inflow_disturbance_end > ny )    &
       THEN
          message_string = 'inflow_disturbance_end out of range'
          CALL message( 'check_parameters', 'PA0132', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

    IF ( random_generator == 'random-parallel' )  THEN
       dist_nxl = nxl;  dist_nxr = nxr
       dist_nys = nys;  dist_nyn = nyn
       IF ( bc_lr == 'radiation/dirichlet' )  THEN
          dist_nxr    = MIN( nx - inflow_disturbance_begin, nxr )
          dist_nxl(1) = MAX( nx - inflow_disturbance_end, nxl )
       ELSEIF ( bc_lr == 'dirichlet/radiation' )  THEN
          dist_nxl    = MAX( inflow_disturbance_begin, nxl )
          dist_nxr(1) = MIN( inflow_disturbance_end, nxr )
       ELSEIF ( bc_lr == 'nested'  .OR.  bc_lr == 'nesting_offline' )  THEN
          dist_nxl    = MAX( inflow_disturbance_begin, nxl )
          dist_nxr    = MIN( nx - inflow_disturbance_begin, nxr )
       ENDIF
       IF ( bc_ns == 'dirichlet/radiation' )  THEN
          dist_nyn    = MIN( ny - inflow_disturbance_begin, nyn )
          dist_nys(1) = MAX( ny - inflow_disturbance_end, nys )
       ELSEIF ( bc_ns == 'radiation/dirichlet' )  THEN
          dist_nys    = MAX( inflow_disturbance_begin, nys )
          dist_nyn(1) = MIN( inflow_disturbance_end, nyn )
       ELSEIF ( bc_ns == 'nested'  .OR.  bc_ns == 'nesting_offline' )  THEN
          dist_nys    = MAX( inflow_disturbance_begin, nys )
          dist_nyn    = MIN( ny - inflow_disturbance_begin, nyn )
       ENDIF
    ELSE
       dist_nxl = 0;  dist_nxr = nx
       dist_nys = 0;  dist_nyn = ny
       IF ( bc_lr == 'radiation/dirichlet' )  THEN
          dist_nxr    = nx - inflow_disturbance_begin
          dist_nxl(1) = nx - inflow_disturbance_end
       ELSEIF ( bc_lr == 'dirichlet/radiation' )  THEN
          dist_nxl    = inflow_disturbance_begin
          dist_nxr(1) = inflow_disturbance_end
       ELSEIF ( bc_lr == 'nested'  .OR.  bc_lr == 'nesting_offline' )  THEN
          dist_nxr    = nx - inflow_disturbance_begin
          dist_nxl    = inflow_disturbance_begin
       ENDIF
       IF ( bc_ns == 'dirichlet/radiation' )  THEN
          dist_nyn    = ny - inflow_disturbance_begin
          dist_nys(1) = ny - inflow_disturbance_end
       ELSEIF ( bc_ns == 'radiation/dirichlet' )  THEN
          dist_nys    = inflow_disturbance_begin
          dist_nyn(1) = inflow_disturbance_end
       ELSEIF ( bc_ns == 'nested'  .OR.  bc_ns == 'nesting_offline' )  THEN
          dist_nyn    = ny - inflow_disturbance_begin
          dist_nys    = inflow_disturbance_begin
       ENDIF
    ENDIF

!
!-- A turbulent inflow requires Dirichlet conditions at the respective inflow
!-- boundary (so far, a turbulent inflow is realized from the left side only)
    IF ( turbulent_inflow  .AND.  bc_lr /= 'dirichlet/radiation' )  THEN
       message_string = 'turbulent_inflow = .T. requires a Dirichlet ' //      &
                        'condition at the inflow boundary'
       CALL message( 'check_parameters', 'PA0133', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Turbulent inflow requires that 3d arrays have been cyclically filled with
!-- data from prerun in the first main run
    IF ( turbulent_inflow  .AND.  initializing_actions /= 'cyclic_fill'  .AND. &
         initializing_actions /= 'read_restart_data' )  THEN
       message_string = 'turbulent_inflow = .T. requires ' //                  &
                        'initializing_actions = ''cyclic_fill'' or ' //        &
                        'initializing_actions = ''read_restart_data'' '
       CALL message( 'check_parameters', 'PA0055', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- In case of turbulent inflow
    IF ( turbulent_inflow )  THEN

!
!--    Calculate the index of the recycling plane
       IF ( recycling_width <= dx  .OR.  recycling_width >= nx * dx )  THEN
          WRITE( message_string, * )  'illegal value for recycling_width: ',   &
                                      recycling_width
          CALL message( 'check_parameters', 'PA0134', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Calculate the index
       recycling_plane = recycling_width / dx
!
!--   Check for correct input of recycling method for thermodynamic quantities
       IF ( TRIM( recycling_method_for_thermodynamic_quantities ) /= 'turbulent_fluctuation' .AND. &
            TRIM( recycling_method_for_thermodynamic_quantities ) /= 'absolute_value' )  THEN
          WRITE( message_string, * )  'unknown recycling method for thermodynamic quantities: ',   &
               TRIM( recycling_method_for_thermodynamic_quantities )
          CALL message( 'check_parameters', 'PA0184', 1, 2, 0, 6, 0 )
       ENDIF

    ENDIF


    IF ( turbulent_outflow )  THEN
!
!--    Turbulent outflow requires Dirichlet conditions at the respective inflow
!--    boundary (so far, a turbulent outflow is realized at the right side only)
       IF ( bc_lr /= 'dirichlet/radiation' )  THEN
          message_string = 'turbulent_outflow = .T. requires ' //              &
                           'bc_lr = "dirichlet/radiation"'
          CALL message( 'check_parameters', 'PA0038', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    The ouflow-source plane must lay inside the model domain
       IF ( outflow_source_plane < dx  .OR.  &
            outflow_source_plane > nx * dx )  THEN
          WRITE( message_string, * )  'illegal value for outflow_source'//     &
                                      '_plane: ', outflow_source_plane
          CALL message( 'check_parameters', 'PA0145', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

!
!-- Determine damping level index for 1D model
    IF ( INDEX( initializing_actions, 'set_1d-model_profiles' ) /= 0 )  THEN
       IF ( damp_level_1d == -1.0_wp )  THEN
          damp_level_1d     = zu(nzt+1)
          damp_level_ind_1d = nzt + 1
       ELSEIF ( damp_level_1d < 0.0_wp  .OR.  damp_level_1d > zu(nzt+1) )  THEN
          WRITE( message_string, * )  'damp_level_1d = ', damp_level_1d,       &
                 ' must be >= 0.0 and <= ', zu(nzt+1), '(zu(nzt+1))'
          CALL message( 'check_parameters', 'PA0136', 1, 2, 0, 6, 0 )
       ELSE
          DO  k = 1, nzt+1
             IF ( damp_level_1d <= zu(k) )  THEN
                damp_level_ind_1d = k
                EXIT
             ENDIF
          ENDDO
       ENDIF
    ENDIF

!
!-- Check some other 1d-model parameters
    IF ( TRIM( mixing_length_1d ) /= 'as_in_3d_model'  .AND.                   &
         TRIM( mixing_length_1d ) /= 'blackadar' )  THEN
       message_string = 'mixing_length_1d = "' // TRIM( mixing_length_1d ) //  &
                        '" is unknown'
       CALL message( 'check_parameters', 'PA0137', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( TRIM( dissipation_1d ) /= 'as_in_3d_model'  .AND.                     &
         TRIM( dissipation_1d ) /= 'detering'  .AND.                           &
         TRIM( dissipation_1d ) /= 'prognostic' )  THEN
       message_string = 'dissipation_1d = "' // TRIM( dissipation_1d ) //      &
                        '" is unknown'
       CALL message( 'check_parameters', 'PA0138', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( TRIM( mixing_length_1d ) /= 'as_in_3d_model'  .AND.                   &
         TRIM( dissipation_1d ) == 'as_in_3d_model' )  THEN
       message_string = 'dissipation_1d = "' // TRIM( dissipation_1d ) //      &
                        '" requires mixing_length_1d = "as_in_3d_model"'
       CALL message( 'check_parameters', 'PA0485', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Set time for the next user defined restart (time_restart is the
!-- internal parameter for steering restart events)
    IF ( restart_time /= 9999999.9_wp )  THEN
       IF ( restart_time > time_since_reference_point )  THEN
          time_restart = restart_time
       ENDIF
    ELSE
!
!--    In case of a restart run, set internal parameter to default (no restart)
!--    if the NAMELIST-parameter restart_time is omitted
       time_restart = 9999999.9_wp
    ENDIF

!
!-- Check pressure gradient conditions
    IF ( dp_external  .AND.  conserve_volume_flow )  THEN
       WRITE( message_string, * )  'Both dp_external and conserve_volume_flo', &
            'w are .TRUE. but one of them must be .FALSE.'
       CALL message( 'check_parameters', 'PA0150', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( dp_external )  THEN
       IF ( dp_level_b < zu(nzb)  .OR.  dp_level_b > zu(nzt) )  THEN
          WRITE( message_string, * )  'dp_level_b = ', dp_level_b, ' is out ', &
               ' of range [zu(nzb), zu(nzt)]'
          CALL message( 'check_parameters', 'PA0151', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( .NOT. ANY( dpdxy /= 0.0_wp ) )  THEN
          WRITE( message_string, * )  'dp_external is .TRUE. but dpdxy is ze', &
               'ro, i.e. the external pressure gradient will not be applied'
          CALL message( 'check_parameters', 'PA0152', 0, 1, 0, 6, 0 )
       ENDIF
    ENDIF
    IF ( ANY( dpdxy /= 0.0_wp )  .AND.  .NOT.  dp_external )  THEN
       WRITE( message_string, * )  'dpdxy is nonzero but dp_external is ',     &
            '.FALSE., i.e. the external pressure gradient & will not be applied'
       CALL message( 'check_parameters', 'PA0153', 0, 1, 0, 6, 0 )
    ENDIF
    IF ( conserve_volume_flow )  THEN
       IF ( TRIM( conserve_volume_flow_mode ) == 'default' )  THEN

          conserve_volume_flow_mode = 'initial_profiles'

       ELSEIF ( TRIM( conserve_volume_flow_mode ) /= 'initial_profiles' .AND.  &
            TRIM( conserve_volume_flow_mode ) /= 'bulk_velocity' )  THEN
          WRITE( message_string, * )  'unknown conserve_volume_flow_mode: ',   &
               conserve_volume_flow_mode
          CALL message( 'check_parameters', 'PA0154', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( (bc_lr /= 'cyclic'  .OR.  bc_ns /= 'cyclic')  .AND.                &
          TRIM( conserve_volume_flow_mode ) == 'bulk_velocity' )  THEN
          WRITE( message_string, * )  'non-cyclic boundary conditions ',       &
               'require  conserve_volume_flow_mode = ''initial_profiles'''
          CALL message( 'check_parameters', 'PA0155', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF
    IF ( ( u_bulk /= 0.0_wp  .OR.  v_bulk /= 0.0_wp )  .AND.                   &
         ( .NOT. conserve_volume_flow  .OR.                                    &
         TRIM( conserve_volume_flow_mode ) /= 'bulk_velocity' ) )  THEN
       WRITE( message_string, * )  'nonzero bulk velocity requires ',          &
            'conserve_volume_flow = .T. and ',                                 &
            'conserve_volume_flow_mode = ''bulk_velocity'''
       CALL message( 'check_parameters', 'PA0157', 1, 2, 0, 6, 0 )
    ENDIF
    
!
!-- Prevent empty time records in volume, cross-section and masked data in case
!-- of non-parallel netcdf-output in restart runs
    IF ( netcdf_data_format < 5 )  THEN
       IF ( TRIM( initializing_actions ) == 'read_restart_data' )  THEN
          do3d_time_count    = 0
          do2d_xy_time_count = 0
          do2d_xz_time_count = 0
          do2d_yz_time_count = 0
          domask_time_count  = 0
       ENDIF
    ENDIF


!
!-- Check roughness length, which has to be smaller than dz/2
    IF ( ( constant_flux_layer .OR.  &
           INDEX( initializing_actions, 'set_1d-model_profiles' ) /= 0 )       &
         .AND. roughness_length >= 0.5 * dz(1) )  THEN
       message_string = 'roughness_length must be smaller than dz/2'
       CALL message( 'check_parameters', 'PA0424', 1, 2, 0, 6, 0 )
    ENDIF

#if defined( __parallel )
!
!-- Vertical nesting: check fine and coarse grid compatibility for data exchange
    IF ( vnested )  CALL vnest_check_parameters
#endif

!
!-- Check if topography is read from file in case of complex terrain simulations
    IF ( complex_terrain  .AND.  TRIM( topography ) /= 'read_from_file' )  THEN
       message_string = 'complex_terrain requires topography' //               &
                        ' = ''read_from_file'''
       CALL message( 'check_parameters', 'PA0295', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Check if vertical grid stretching is switched off in case of complex 
!-- terrain simulations
    IF ( complex_terrain  .AND.                                                &
         ANY( dz_stretch_level_start /= -9999999.9_wp ) )  THEN
       message_string = 'Vertical grid stretching is not allowed for ' //      &
                        'complex_terrain = .T.'
       CALL message( 'check_parameters', 'PA0473', 1, 2, 0, 6, 0 )
    ENDIF

    CALL location_message( 'checking parameters', 'finished' )

 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check the length of data output intervals. In case of parallel NetCDF output
!> the time levels of the output files need to be fixed. Therefore setting the
!> output interval to 0.0s (usually used to output each timestep) is not
!> possible as long as a non-fixed timestep is used.
!------------------------------------------------------------------------------!

    SUBROUTINE check_dt_do( dt_do, dt_do_name )

       IMPLICIT NONE

       CHARACTER (LEN=*), INTENT (IN) :: dt_do_name !< parin variable name

       REAL(wp), INTENT (INOUT)       :: dt_do      !< data output interval

       IF ( dt_do == 0.0_wp )  THEN
          IF ( dt_fixed )  THEN
             WRITE( message_string, '(A,F9.4,A)' )  'Output at every '  //     &
                    'timestep is wanted (' // dt_do_name // ' = 0.0).&'//      &
                    'The output interval is set to the fixed timestep dt '//   &
                    '= ', dt, 's.'
             CALL message( 'check_parameters', 'PA0060', 0, 0, 0, 6, 0 )
             dt_do = dt
          ELSE
             message_string = dt_do_name // ' = 0.0 while using a ' //         &
                              'variable timestep and parallel netCDF4 ' //     &
                              'is not allowed.'
             CALL message( 'check_parameters', 'PA0081', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF

    END SUBROUTINE check_dt_do



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Set the bottom and top boundary conditions for humidity and scalars.
!------------------------------------------------------------------------------!

    SUBROUTINE set_bc_scalars( sq, bc_b, bc_t, ibc_b, ibc_t, err_nr_b, err_nr_t )


       IMPLICIT NONE

       CHARACTER (LEN=1)   ::  sq         !< name of scalar quantity
       CHARACTER (LEN=*)   ::  bc_b       !< bottom boundary condition
       CHARACTER (LEN=*)   ::  bc_t       !< top boundary condition
       CHARACTER (LEN=*)   ::  err_nr_b   !< error number if bottom bc is unknown
       CHARACTER (LEN=*)   ::  err_nr_t   !< error number if top bc is unknown

       INTEGER(iwp)        ::  ibc_b      !< index for bottom boundary condition
       INTEGER(iwp)        ::  ibc_t      !< index for top boundary condition

!
!--    Set Integer flags and check for possilbe errorneous settings for bottom
!--    boundary condition
       IF ( bc_b == 'dirichlet' )  THEN
          ibc_b = 0
       ELSEIF ( bc_b == 'neumann' )  THEN
          ibc_b = 1
       ELSE
          message_string = 'unknown boundary condition: bc_' // TRIM( sq ) //  &
                           '_b ="' // TRIM( bc_b ) // '"'
          CALL message( 'check_parameters', err_nr_b, 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Set Integer flags and check for possilbe errorneous settings for top
!--    boundary condition
       IF ( bc_t == 'dirichlet' )  THEN
          ibc_t = 0
       ELSEIF ( bc_t == 'neumann' )  THEN
          ibc_t = 1
       ELSEIF ( bc_t == 'initial_gradient' )  THEN
          ibc_t = 2
       ELSEIF ( bc_t == 'nested'  .OR.  bc_t == 'nesting_offline' )  THEN
          ibc_t = 3
       ELSE
          message_string = 'unknown boundary condition: bc_' // TRIM( sq ) //  &
                           '_t ="' // TRIM( bc_t ) // '"'
          CALL message( 'check_parameters', err_nr_t, 1, 2, 0, 6, 0 )
       ENDIF


    END SUBROUTINE set_bc_scalars



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check for consistent settings of bottom boundary conditions for humidity
!> and scalars.
!------------------------------------------------------------------------------!

    SUBROUTINE check_bc_scalars( sq, bc_b, ibc_b,                      &
                                 err_nr_1, err_nr_2,                   &
                                 constant_flux, surface_initial_change )


       IMPLICIT NONE

       CHARACTER (LEN=1)   ::  sq                       !< name of scalar quantity
       CHARACTER (LEN=*)   ::  bc_b                     !< bottom boundary condition
       CHARACTER (LEN=*)   ::  err_nr_1                 !< error number of first error
       CHARACTER (LEN=*)   ::  err_nr_2                 !< error number of second error

       INTEGER(iwp)        ::  ibc_b                    !< index of bottom boundary condition

       LOGICAL             ::  constant_flux            !< flag for constant-flux layer

       REAL(wp)            ::  surface_initial_change   !< value of initial change at the surface

!
!--    A given surface value implies Dirichlet boundary condition for
!--    the respective quantity. In this case specification of a constant flux is
!--    forbidden. However, an exception is made for large-scale forcing as well
!--    as land-surface model.
       IF ( .NOT. land_surface  .AND.  .NOT. large_scale_forcing )  THEN
          IF ( ibc_b == 0  .AND.  constant_flux )  THEN
             message_string = 'boundary condition: bc_' // TRIM( sq ) //       &
                              '_b ' // '= "' // TRIM( bc_b ) //                &
                              '" is not allowed with prescribed surface flux'
             CALL message( 'check_parameters', err_nr_1, 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF
       IF ( constant_flux  .AND.  surface_initial_change /= 0.0_wp )  THEN
          WRITE( message_string, * )  'a prescribed surface flux is not allo', &
                 'wed with ', sq, '_surface_initial_change (/=0) = ',          &
                 surface_initial_change
          CALL message( 'check_parameters', err_nr_2, 1, 2, 0, 6, 0 )
       ENDIF


    END SUBROUTINE check_bc_scalars



 END SUBROUTINE check_parameters
