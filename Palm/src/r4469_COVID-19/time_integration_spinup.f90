!> @file time_integration_spinup.f90
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
! $Id: time_integration_spinup.f90 4457 2020-03-11 14:20:43Z raasch $
! use statement for exchange horiz added
! 
! 4444 2020-03-05 15:59:50Z raasch
! bugfix: cpp-directives for serial mode added
! 
! 4360 2020-01-07 11:25:50Z suehring
! Enable output of diagnostic quantities, e.g. 2-m temperature
! 
! 4227 2019-09-10 18:04:34Z gronemeier
! implement new palm_date_time_mod
! 
! 4223 2019-09-10 09:20:47Z gronemeier
! Corrected "Former revisions" section
! 
! 4064 2019-07-01 05:33:33Z gronemeier
! Moved call to radiation module out of intermediate time loop
! 
! 4023 2019-06-12 13:20:01Z maronga
! Time stamps are now negative in run control output
! 
! 3885 2019-04-11 11:29:34Z kanani
! Changes related to global restructuring of location messages and introduction 
! of additional debug messages
! 
! 3766 2019-02-26 16:23:41Z raasch
! unused variable removed
! 
! 3719 2019-02-06 13:10:18Z kanani
! Removed log_point(19,54,74,50,75), since they count together with same log
! points in time_integration, impossible to separate the contributions.
! Instead, the entire spinup gets an individual log_point in palm.f90
! 
! 3655 2019-01-07 16:51:22Z knoop
! Removed call to calculation of near air (10 cm) potential temperature (now in
! surface layer fluxes)
! 
! 2296 2017-06-28 07:53:56Z maronga
! Initial revision
!
!
! Description:
! ------------
!> Integration in time of the non-atmospheric model components such as land
!> surface model and urban surface model
!------------------------------------------------------------------------------!
 SUBROUTINE time_integration_spinup
 
    USE arrays_3d,                                                             &
        ONLY:  pt, pt_p, u, u_init, v, v_init

    USE control_parameters,                                                    &
        ONLY:  averaging_interval_pr, calc_soil_moisture_during_spinup,        &
               constant_diffusion, constant_flux_layer, coupling_start_time,   &
               data_output_during_spinup, dopr_n, do_sum,                      &
               dt_averaging_input_pr, dt_dopr, dt_dots, dt_do2d_xy, dt_do3d,   &
               dt_spinup, dt_3d, humidity, intermediate_timestep_count,        &
               intermediate_timestep_count_max, land_surface,                  &
               simulated_time, simulated_time_chr, skip_time_dopr,             &
               skip_time_do2d_xy, skip_time_do3d, spinup_pt_amplitude,         &
               spinup_pt_mean, spinup_time, timestep_count, time_dopr,         &
               time_dopr_av, time_dots, time_do2d_xy, time_do3d,               &
               time_run_control, time_since_reference_point, urban_surface

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point_s

    USE diagnostic_output_quantities_mod,                                      &
        ONLY:  doq_calculate

    USE exchange_horiz_mod,                                                    &
        ONLY:  exchange_horiz

    USE indices,                                                               &
        ONLY:  nbgp, nzb, nzt, nysg, nyng, nxlg, nxrg

    USE land_surface_model_mod,                                                &
        ONLY:  lsm_energy_balance, lsm_soil_model, lsm_swap_timelevel

    USE pegrid

#if defined( __parallel )
    USE pmc_interface,                                                         &
        ONLY:  nested_run
#endif

    USE kinds

    USE palm_date_time_mod,                                                    &
        ONLY:  get_date_time, seconds_per_hour

    USE radiation_model_mod,                                                   &
        ONLY:  force_radiation_call, radiation, radiation_control,             &
               radiation_interaction, radiation_interactions, time_radiation

    USE statistics,                                                            &
        ONLY:  flow_statistics_called

    USE surface_layer_fluxes_mod,                                              &
        ONLY:  surface_layer_fluxes

    USE surface_mod,                                                           &
        ONLY :  surf_lsm_h, surf_lsm_v, surf_usm_h,    &
                surf_usm_v

    USE urban_surface_mod,                                                     &
        ONLY:  usm_material_heat_model, usm_material_model,                    &
               usm_surface_energy_balance, usm_swap_timelevel,                 &
               usm_green_heat_model




    IMPLICIT NONE

    CHARACTER (LEN=9) ::  time_to_string                  !< 
  
  
    CHARACTER (LEN=1) ::  sign_chr                        !< String containing '-' or ' ' 
    CHARACTER (LEN=9) ::  time_since_reference_point_chr  !< time since reference point, i.e., negative during spinup
  
    INTEGER(iwp) ::  i !< running index
    INTEGER(iwp) ::  j !< running index
    INTEGER(iwp) ::  k !< running index
    INTEGER(iwp) ::  l !< running index
    INTEGER(iwp) ::  m !< running index

    INTEGER(iwp) :: current_timestep_number_spinup = 0  !< number if timestep during spinup
    INTEGER(iwp) :: day_of_year                         !< day of the year
  
    LOGICAL :: run_control_header_spinup = .FALSE.  !< flag parameter for steering whether the header information must be output

    REAL(wp) ::  pt_spinup      !< temporary storage of temperature
    REAL(wp) ::  dt_save        !< temporary storage for time step
    REAL(wp) ::  second_of_day  !< second of the day
                  
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  pt_save  !< temporary storage of temperature
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  u_save   !< temporary storage of u wind component
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  v_save   !< temporary storage of v wind component


!
!-- Save 3D arrays because they are to be changed for spinup purpose
    ALLOCATE( pt_save(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( u_save(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( v_save(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

    CALL exchange_horiz( pt, nbgp )   
    CALL exchange_horiz( u,  nbgp )  
    CALL exchange_horiz( v,  nbgp )  
 
    pt_save = pt
    u_save  = u
    v_save  = v

!
!-- Set the same wall-adjacent velocity to all grid points. The sign of the 
!-- original velocity field must be preserved because the surface schemes crash 
!-- otherwise. The precise reason is still unknown. A minimum velocity of 0.1
!-- m/s is used to maintain turbulent transfer at the surface.
    IF ( land_surface )  THEN
       DO  m = 1, surf_lsm_h%ns
          i   = surf_lsm_h%i(m)            
          j   = surf_lsm_h%j(m)
          k   = surf_lsm_h%k(m)
          u(k,j,i) = SIGN(1.0_wp,u_init(k)) * MAX( ABS( u_init(k) ),0.1_wp)
          v(k,j,i) = SIGN(1.0_wp,v_init(k)) * MAX( ABS( v_init(k) ),0.1_wp)
       ENDDO

       DO  l = 0, 3
          DO  m = 1, surf_lsm_v(l)%ns
             i   = surf_lsm_v(l)%i(m)            
             j   = surf_lsm_v(l)%j(m)
             k   = surf_lsm_v(l)%k(m)
             u(k,j,i) = SIGN(1.0_wp,u_init(k)) * MAX( ABS( u_init(k) ),0.1_wp)
             v(k,j,i) = SIGN(1.0_wp,v_init(k)) * MAX( ABS( v_init(k) ),0.1_wp)
          ENDDO
       ENDDO
    ENDIF

    IF ( urban_surface )  THEN
       DO  m = 1, surf_usm_h%ns
          i   = surf_usm_h%i(m)            
          j   = surf_usm_h%j(m)
          k   = surf_usm_h%k(m)
          u(k,j,i) = SIGN(1.0_wp,u_init(k)) * MAX( ABS( u_init(k) ),0.1_wp)
          v(k,j,i) = SIGN(1.0_wp,v_init(k)) * MAX( ABS( v_init(k) ),0.1_wp)
       ENDDO

       DO  l = 0, 3
          DO  m = 1, surf_usm_v(l)%ns
             i   = surf_usm_v(l)%i(m)            
             j   = surf_usm_v(l)%j(m)
             k   = surf_usm_v(l)%k(m)
             u(k,j,i) = SIGN(1.0_wp,u_init(k)) * MAX( ABS( u_init(k) ),0.1_wp)
             v(k,j,i) = SIGN(1.0_wp,v_init(k)) * MAX( ABS( v_init(k) ),0.1_wp)
          ENDDO
       ENDDO
    ENDIF

    CALL exchange_horiz( u,  nbgp )
    CALL exchange_horiz( v,  nbgp )

    dt_save = dt_3d
    dt_3d   = dt_spinup

    CALL location_message( 'wall/soil spinup time-stepping', 'start' )
!
!-- Start of the time loop
    DO  WHILE ( simulated_time < spinup_time )

       CALL cpu_log( log_point_s(15), 'timesteps spinup', 'start' )
   
!
!--    Start of intermediate step loop
       intermediate_timestep_count = 0
       DO  WHILE ( intermediate_timestep_count < &
                   intermediate_timestep_count_max )

          intermediate_timestep_count = intermediate_timestep_count + 1

!
!--       Set the steering factors for the prognostic equations which depend
!--       on the timestep scheme
          CALL timestep_scheme_steering


!
!--       Estimate a near-surface air temperature based on the position of the 
!--       sun and user input about mean temperature and amplitude. The time is
!--       shifted by one hour to simulate a lag between air temperature and
!--       incoming radiation
          CALL get_date_time( simulated_time - spinup_time - seconds_per_hour, &
                              day_of_year=day_of_year,                         &
                              second_of_day=second_of_day                      )

          pt_spinup = spinup_pt_mean + spinup_pt_amplitude                     &
                                     * solar_angle(day_of_year, second_of_day)

!
!--       Map air temperature to all grid points in the vicinity of a surface
!--       element
          IF ( land_surface )  THEN
             DO  m = 1, surf_lsm_h%ns
                i   = surf_lsm_h%i(m)            
                j   = surf_lsm_h%j(m)
                k   = surf_lsm_h%k(m)
                pt(k,j,i) = pt_spinup
             ENDDO

             DO  l = 0, 3
                DO  m = 1, surf_lsm_v(l)%ns
                   i   = surf_lsm_v(l)%i(m)            
                   j   = surf_lsm_v(l)%j(m)
                   k   = surf_lsm_v(l)%k(m)
                   pt(k,j,i) = pt_spinup
                ENDDO
             ENDDO
          ENDIF

          IF ( urban_surface )  THEN
             DO  m = 1, surf_usm_h%ns
                i   = surf_usm_h%i(m)            
                j   = surf_usm_h%j(m)
                k   = surf_usm_h%k(m)
                pt(k,j,i) = pt_spinup
                !!!!!!!!!!!!!!!!HACK!!!!!!!!!!!!!
                surf_usm_h%pt1 = pt_spinup
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ENDDO

             DO  l = 0, 3
                DO  m = 1, surf_usm_v(l)%ns
                   i   = surf_usm_v(l)%i(m)            
                   j   = surf_usm_v(l)%j(m)
                   k   = surf_usm_v(l)%k(m)
                   pt(k,j,i) = pt_spinup
                   !!!!!!!!!!!!!!!!HACK!!!!!!!!!!!!!
                   surf_usm_v(l)%pt1 = pt_spinup
                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ENDDO
             ENDDO
          ENDIF

          CALL exchange_horiz( pt,  nbgp )    


!
!--       Swap the time levels in preparation for the next time step.
          timestep_count = timestep_count + 1
     
          IF ( land_surface )  THEN
              CALL lsm_swap_timelevel ( 0 )
          ENDIF

          IF ( urban_surface )  THEN
             CALL usm_swap_timelevel ( 0 )
          ENDIF

          IF ( land_surface )  THEN
             CALL lsm_swap_timelevel ( MOD( timestep_count, 2) )
          ENDIF

          IF ( urban_surface )  THEN
             CALL usm_swap_timelevel ( MOD( timestep_count, 2) )
          ENDIF
         
!
!--       If required, compute virtual potential temperature 
          IF ( humidity )  THEN 
             CALL compute_vpt 
          ENDIF 

!
!--       Compute the diffusion quantities
          IF ( .NOT. constant_diffusion )  THEN

!
!--          First the vertical (and horizontal) fluxes in the surface 
!--          (constant flux) layer are computed
             IF ( constant_flux_layer )  THEN
                CALL surface_layer_fluxes
             ENDIF

!
!--          If required, solve the energy balance for the surface and run soil 
!--          model. Call for horizontal as well as vertical surfaces.
!--          The prognostic equation for soil moisure is switched off
             IF ( land_surface )  THEN

!
!--             Call for horizontal upward-facing surfaces
                CALL lsm_energy_balance( .TRUE., -1 )
                CALL lsm_soil_model( .TRUE., -1, calc_soil_moisture_during_spinup )
!
!--             Call for northward-facing surfaces
                CALL lsm_energy_balance( .FALSE., 0 )
                CALL lsm_soil_model( .FALSE., 0, calc_soil_moisture_during_spinup )
!
!--             Call for southward-facing surfaces
                CALL lsm_energy_balance( .FALSE., 1 )
                CALL lsm_soil_model( .FALSE., 1, calc_soil_moisture_during_spinup )
!
!--             Call for eastward-facing surfaces
                CALL lsm_energy_balance( .FALSE., 2 )
                CALL lsm_soil_model( .FALSE., 2, calc_soil_moisture_during_spinup )
!
!--             Call for westward-facing surfaces
                CALL lsm_energy_balance( .FALSE., 3 )
                CALL lsm_soil_model( .FALSE., 3, calc_soil_moisture_during_spinup )

             ENDIF

!
!--          If required, solve the energy balance for urban surfaces and run 
!--          the material heat model
             IF (urban_surface) THEN

                CALL usm_surface_energy_balance( .TRUE. )
                IF ( usm_material_model )  THEN
                   CALL usm_green_heat_model
                   CALL usm_material_heat_model( .TRUE. )
                ENDIF

             ENDIF

          ENDIF

       ENDDO   ! Intermediate step loop

!
!--    If required, calculate radiative fluxes and heating rates
       IF ( radiation )  THEN

            time_radiation = time_radiation + dt_3d

          IF ( time_radiation >= dt_3d .OR. force_radiation_call )  THEN

             IF ( .NOT. force_radiation_call )  THEN
                time_radiation = time_radiation - dt_3d
             ENDIF

             CALL radiation_control

             IF ( radiation_interactions )  THEN
                CALL radiation_interaction
             ENDIF
          ENDIF
       ENDIF

!
!--    Increase simulation time and output times
       current_timestep_number_spinup = current_timestep_number_spinup + 1
       simulated_time             = simulated_time   + dt_3d
       simulated_time_chr         = time_to_string( simulated_time )
       time_since_reference_point = simulated_time - coupling_start_time
       time_since_reference_point_chr = time_to_string( ABS(time_since_reference_point) )
       
       IF ( time_since_reference_point < 0.0_wp )  THEN
          sign_chr = '-'
       ELSE
          sign_chr = ' '
       ENDIF
      
       
       IF ( data_output_during_spinup )  THEN
          IF ( simulated_time >= skip_time_do2d_xy )  THEN
             time_do2d_xy       = time_do2d_xy     + dt_3d
          ENDIF
          IF ( simulated_time >= skip_time_do3d    )  THEN
             time_do3d          = time_do3d        + dt_3d
          ENDIF
          time_dots          = time_dots        + dt_3d
          IF ( simulated_time >= skip_time_dopr )  THEN
             time_dopr       = time_dopr        + dt_3d
          ENDIF
          time_run_control   = time_run_control + dt_3d

!
!--       Carry out statistical analysis and output at the requested output times.
!--       The MOD function is used for calculating the output time counters (like
!--       time_dopr) in order to regard a possible decrease of the output time
!--       interval in case of restart runs

!
!--       Set a flag indicating that so far no statistics have been created
!--       for this time step
          flow_statistics_called = .FALSE.

!
!--       If required, call flow_statistics for averaging in time
          IF ( averaging_interval_pr /= 0.0_wp  .AND.                          &
             ( dt_dopr - time_dopr ) <= averaging_interval_pr  .AND.           &
             simulated_time >= skip_time_dopr )  THEN
             time_dopr_av = time_dopr_av + dt_3d
             IF ( time_dopr_av >= dt_averaging_input_pr )  THEN
                do_sum = .TRUE.
                time_dopr_av = MOD( time_dopr_av,                              &
                               MAX( dt_averaging_input_pr, dt_3d ) )
             ENDIF
          ENDIF
          IF ( do_sum )  CALL flow_statistics

!
!--       Output of profiles
          IF ( time_dopr >= dt_dopr )  THEN
             IF ( dopr_n /= 0 )  CALL data_output_profiles
             time_dopr = MOD( time_dopr, MAX( dt_dopr, dt_3d ) )
             time_dopr_av = 0.0_wp    ! due to averaging (see above)
          ENDIF

!
!--       Output of time series
          IF ( time_dots >= dt_dots )  THEN
             CALL data_output_tseries
             time_dots = MOD( time_dots, MAX( dt_dots, dt_3d ) )
          ENDIF

!
!--       2d-data output (cross-sections)
          IF ( time_do2d_xy >= dt_do2d_xy )  THEN
             CALL doq_calculate
             CALL data_output_2d( 'xy', 0 )
             time_do2d_xy = MOD( time_do2d_xy, MAX( dt_do2d_xy, dt_3d ) )
          ENDIF

!
!--       3d-data output (volume data)
          IF ( time_do3d >= dt_do3d )  THEN
             CALL doq_calculate
             CALL data_output_3d( 0 )
             time_do3d = MOD( time_do3d, MAX( dt_do3d, dt_3d ) )
          ENDIF


       ENDIF

!
!--    Computation and output of run control parameters.
!--    This is also done whenever perturbations have been imposed
!        IF ( time_run_control >= dt_run_control  .OR.                           &
!             timestep_scheme(1:5) /= 'runge'  .OR.  disturbance_created )       &
!        THEN
!           CALL run_control
!           IF ( time_run_control >= dt_run_control )  THEN
!              time_run_control = MOD( time_run_control,                         &
!                                      MAX( dt_run_control, dt_3d ) )
!           ENDIF
!        ENDIF

       CALL cpu_log( log_point_s(15), 'timesteps spinup', 'stop' )


!
!--    Run control output
       IF ( myid == 0 )  THEN
!
!--       If necessary, write header
          IF ( .NOT. run_control_header_spinup )  THEN
             CALL check_open( 15 )
             WRITE ( 15, 100 )
             run_control_header_spinup = .TRUE.
          ENDIF
!
!--       Write some general information about the spinup in run control file
          WRITE ( 15, 101 )  current_timestep_number_spinup, sign_chr, time_since_reference_point_chr, dt_3d, pt_spinup
!
!--       Write buffer contents to disc immediately
          FLUSH( 15 )
       ENDIF



    ENDDO   ! time loop

!
!-- Write back saved arrays to the 3D arrays
    pt   = pt_save
    pt_p = pt_save
    u    = u_save
    v    = v_save

!
!-- Reset time step
    dt_3d = dt_save

    DEALLOCATE(pt_save)
    DEALLOCATE(u_save)
    DEALLOCATE(v_save)

#if defined( __parallel )
    IF ( nested_run )  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
#endif

    CALL location_message( 'wall/soil spinup time-stepping', 'finished' )


!
!-- Formats
100 FORMAT (///'Spinup control output:'/  &
            '---------------------------------'// &
            'ITER.   HH:MM:SS    DT   PT(z_MO)'/   &
            '---------------------------------')
101 FORMAT (I5,2X,A1,A9,1X,F6.2,3X,F6.2,2X,F6.2)

 CONTAINS

!
!-- Returns the cosine of the solar zenith angle at a given time. This routine
!-- is similar to that for calculation zenith (see radiation_model_mod.f90)
    !> @todo Load function calc_zenith of radiation model instead of
    !>       rewrite the function here.
    FUNCTION solar_angle( day_of_year, second_of_day ) 

       USE basic_constants_and_equations_mod,                                  &
           ONLY:  pi
      
       USE kinds

       USE radiation_model_mod,                                                &
           ONLY:  decl_1, decl_2, decl_3, lat, lon

       IMPLICIT NONE 


       INTEGER(iwp), INTENT(IN) ::  day_of_year  !< day of the year

       REAL(wp)             ::  declination      !< solar declination angle
       REAL(wp)             ::  hour_angle       !< solar hour angle
       REAL(wp), INTENT(IN) ::  second_of_day    !< current time of the day in UTC
       REAL(wp)             ::  solar_angle      !< cosine of the solar zenith angle
!
!--    Calculate solar declination and hour angle   
       declination = ASIN( decl_1 * SIN(decl_2 * REAL(day_of_year, KIND=wp) - decl_3) )
       hour_angle  = 2.0_wp * pi * (second_of_day / 86400.0_wp) + lon - pi

!
!--    Calculate cosine of solar zenith angle
       solar_angle = SIN(lat) * SIN(declination) + COS(lat) * COS(declination) &
                     * COS(hour_angle)

    END FUNCTION solar_angle


 END SUBROUTINE time_integration_spinup
