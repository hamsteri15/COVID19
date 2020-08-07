!> @file virtual_flights_mod.f90
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
! $Id: virtual_flight_mod.f90 4360 2020-01-07 11:25:50Z suehring $
! Corrected "Former revisions" section
! 
! 4004 2019-05-24 11:32:38Z suehring
! Allow variable start- and end locations also in return mode
! 
! 3885 2019-04-11 11:29:34Z kanani
! Changes related to global restructuring of location messages and introduction 
! of additional debug messages
! 
! 3655 2019-01-07 16:51:22Z knoop
! variables documented
! 
! 1957 2016-07-07 10:43:48Z suehring
! Initial revision
!
! Description:
! ------------
!> Module for virtual flight measurements.
!> @todo Err msg PA0438: flight can be inside topography -> extra check?
!--------------------------------------------------------------------------------!
 MODULE flight_mod
 
    USE control_parameters,                                                    &
        ONLY:  debug_output, fl_max, num_leg, num_var_fl, num_var_fl_user, virtual_flight
 
    USE kinds

    CHARACTER(LEN=6), DIMENSION(fl_max) ::  leg_mode = 'cyclic'  !< flight mode through the model domain, either 'cyclic' or 'return' 

    INTEGER(iwp) ::  l           !< index for flight leg
    INTEGER(iwp) ::  var_index   !< index for measured variable

    LOGICAL, DIMENSION(:), ALLOCATABLE  ::  cyclic_leg !< flag to identify fly mode

    REAL(wp) ::  flight_end = 9999999.9_wp  !< end time of virtual flight
    REAL(wp) ::  flight_begin = 0.0_wp      !< end time of virtual flight

    REAL(wp), DIMENSION(fl_max) ::  flight_angle = 45.0_wp   !< angle determining the horizontal flight direction
    REAL(wp), DIMENSION(fl_max) ::  flight_level = 100.0_wp  !< flight level 
    REAL(wp), DIMENSION(fl_max) ::  max_elev_change = 0.0_wp !< maximum elevation change for the respective flight leg 
    REAL(wp), DIMENSION(fl_max) ::  rate_of_climb = 0.0_wp   !< rate of climb or descent
    REAL(wp), DIMENSION(fl_max) ::  speed_agl = 25.0_wp      !< absolute horizontal flight speed above ground level (agl)
    REAL(wp), DIMENSION(fl_max) ::  x_start = 999999999.0_wp !< start x position 
    REAL(wp), DIMENSION(fl_max) ::  x_end   = 999999999.0_wp !< end x position
    REAL(wp), DIMENSION(fl_max) ::  y_start = 999999999.0_wp !< start y position
    REAL(wp), DIMENSION(fl_max) ::  y_end   = 999999999.0_wp !< end y position

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  u_agl      !< u-component of flight speed
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  v_agl      !< v-component of flight speed
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  w_agl      !< w-component of flight speed
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  x_pos      !< aircraft x-position 
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  y_pos      !< aircraft y-position 
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  z_pos      !< aircraft z-position 

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  sensor_l !< measured data on local PE
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  sensor   !< measured data

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  var_u  !< dummy array for possibly user-defined quantities

    SAVE

    PRIVATE

    INTERFACE flight_header
       MODULE PROCEDURE flight_header
    END INTERFACE flight_header
    
    INTERFACE flight_init
       MODULE PROCEDURE flight_init
    END INTERFACE flight_init

    INTERFACE flight_init_output
       MODULE PROCEDURE flight_init_output
    END INTERFACE flight_init_output

    INTERFACE flight_check_parameters
       MODULE PROCEDURE flight_check_parameters
    END INTERFACE flight_check_parameters

    INTERFACE flight_parin
       MODULE PROCEDURE flight_parin
    END INTERFACE flight_parin

    INTERFACE interpolate_xyz
       MODULE PROCEDURE interpolate_xyz
    END INTERFACE interpolate_xyz

    INTERFACE flight_measurement
       MODULE PROCEDURE flight_measurement
    END INTERFACE flight_measurement
    
    INTERFACE flight_rrd_global 
       MODULE PROCEDURE flight_rrd_global
    END INTERFACE flight_rrd_global
    
    INTERFACE flight_wrd_global 
       MODULE PROCEDURE flight_wrd_global 
    END INTERFACE flight_wrd_global 

!
!-- Private interfaces
    PRIVATE flight_check_parameters, flight_init_output, interpolate_xyz
!
!-- Public interfaces
    PUBLIC flight_init, flight_header, flight_parin, flight_measurement,       &
           flight_wrd_global, flight_rrd_global                   
!
!-- Public variables
    PUBLIC fl_max, sensor, x_pos, y_pos, z_pos

 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Header output for flight module.
!------------------------------------------------------------------------------!
    SUBROUTINE flight_header ( io )

    
       IMPLICIT NONE

       INTEGER(iwp), INTENT(IN) ::  io            !< Unit of the output file

       WRITE ( io, 1  )
       WRITE ( io, 2  )
       WRITE ( io, 3  ) num_leg
       WRITE ( io, 4  ) flight_begin
       WRITE ( io, 5  ) flight_end
       
       DO l=1, num_leg
          WRITE ( io, 6   ) l
          WRITE ( io, 7   ) speed_agl(l)
          WRITE ( io, 8   ) flight_level(l)
          WRITE ( io, 9   ) max_elev_change(l)
          WRITE ( io, 10  ) rate_of_climb(l)
          WRITE ( io, 11  ) leg_mode(l)
       ENDDO

       
     1   FORMAT (' Virtual flights:'/                                           &
               ' ----------------')
     2   FORMAT ('       Output every timestep')
     3   FORMAT ('       Number of flight legs:',    I3                )
     4   FORMAT ('       Begin of measurements:',    F10.1    , ' s'   )
     5   FORMAT ('       End of measurements:',      F10.1    , ' s'   )
     6   FORMAT ('       Leg', I3/,                                             &
                '       ------' )
     7   FORMAT ('          Flight speed            : ', F5.1, ' m/s' )
     8   FORMAT ('          Flight level            : ', F5.1, ' m'   )
     9   FORMAT ('          Maximum elevation change: ', F5.1, ' m/s' )
     10  FORMAT ('          Rate of climb / descent : ', F5.1, ' m/s' )
     11  FORMAT ('          Leg mode                : ', A/           )
 
    END SUBROUTINE flight_header 
 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads the namelist flight_par.
!------------------------------------------------------------------------------!
    SUBROUTINE flight_parin 

       USE control_parameters,                                                 &
           ONLY:  message_string
      
       IMPLICIT NONE

       CHARACTER (LEN=80) ::  line  !< dummy string that contains the current line of the parameter file

       NAMELIST /flight_par/  flight_angle, flight_end, flight_begin, leg_mode,&
                              flight_level, max_elev_change, rate_of_climb,    &
                              speed_agl, x_end, x_start, y_end, y_start
  

       NAMELIST /virtual_flight_parameters/                                    &
                              flight_angle, flight_end, flight_begin, leg_mode,&
                              flight_level, max_elev_change, rate_of_climb,    &
                              speed_agl, x_end, x_start, y_end, y_start
!
!--    Try to find the namelist flight_par
       REWIND ( 11 )
       line = ' '
       DO WHILE ( INDEX( line, '&virtual_flight_parameters' ) == 0 )
          READ ( 11, '(A)', END=12 )  line
       ENDDO
       BACKSPACE ( 11 )

!
!--    Read namelist
       READ ( 11, virtual_flight_parameters, ERR = 10 )
!
!--    Set switch that virtual flights shall be carried out
       virtual_flight = .TRUE.

       GOTO 14

 10    BACKSPACE( 11 )
       READ( 11 , '(A)') line
       CALL parin_fail_message( 'virtual_flight_parameters', line )
!
!--    Try to find the old namelist
 12    REWIND ( 11 )
       line = ' '
       DO WHILE ( INDEX( line, '&flight_par' ) == 0 )
          READ ( 11, '(A)', END=14 )  line
       ENDDO
       BACKSPACE ( 11 )

!
!--    Read namelist
       READ ( 11, flight_par, ERR = 13, END = 14 )
       
       message_string = 'namelist flight_par is deprecated and will be ' // &
                     'removed in near future.& Please use namelist ' //     &
                     'virtual_flight_parameters instead' 
       CALL message( 'flight_parin', 'PA0487', 0, 1, 0, 6, 0 )     
!
!--    Set switch that virtual flights shall be carried out
       virtual_flight = .TRUE.

       GOTO 14

 13    BACKSPACE( 11 )
       READ( 11 , '(A)') line
       CALL parin_fail_message( 'flight_par', line )

 14    CONTINUE

    END SUBROUTINE flight_parin

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Inititalization of required arrays, number of legs and flags. Moreover,
!> initialize flight speed and -direction, as well as start positions.
!------------------------------------------------------------------------------!
    SUBROUTINE flight_init
 
       USE basic_constants_and_equations_mod,                                  &
           ONLY:  pi
    
       USE control_parameters,                                                 &
           ONLY:  initializing_actions 
                  
       USE indices,                                                            &
           ONLY:  nxlg, nxrg, nysg, nyng, nzb, nzt

       IMPLICIT NONE

       REAL(wp) ::  distance  !< distance between start and end position of a flight leg


       IF ( debug_output )  CALL debug_message( 'flight_init', 'start' )
!
!--    Determine the number of flight legs
       l = 1
       DO WHILE ( x_start(l) /= 999999999.0_wp  .AND.  l <= SIZE(x_start) )
          l       = l + 1
       ENDDO
       num_leg = l-1
!
!--    Check for proper parameter settings
       CALL flight_check_parameters
!
!--    Allocate and initialize logical array for flight pattern
       ALLOCATE( cyclic_leg(1:num_leg) )
!
!--    Initialize flags for cyclic/return legs
       DO l = 1, num_leg
          cyclic_leg(l) = MERGE( .TRUE., .FALSE.,                              &
                                 TRIM( leg_mode(l) ) == 'cyclic'               &
                               )
       ENDDO
!
!--    Allocate and initialize arraxs for flight position and speed. In case of 
!--    restart runs these data are read by the routine read_flight_restart_data
!--    instead.
       IF (  TRIM( initializing_actions ) /= 'read_restart_data' )  THEN 
       
          ALLOCATE( x_pos(1:num_leg), y_pos(1:num_leg ), z_pos(1:num_leg) )
!
!--       Inititalize x-, y-, and z-positions with initial start position          
          x_pos(1:num_leg) = x_start(1:num_leg)
          y_pos(1:num_leg) = y_start(1:num_leg)
          z_pos(1:num_leg) = flight_level(1:num_leg)
!
!--       Allocate arrays for flight-speed components
          ALLOCATE( u_agl(1:num_leg),                                          &
                    v_agl(1:num_leg),                                          &
                    w_agl(1:num_leg) )
!
!--       Inititalize u-, v- and w-component. 
          DO  l = 1, num_leg
!
!--          In case of return-legs, the flight direction, i.e. the horizontal
!--          flight-speed components, are derived from the given start/end 
!--          positions.
             IF (  .NOT.  cyclic_leg(l) )  THEN
                distance = SQRT( ( x_end(l) - x_start(l) )**2                  &
                               + ( y_end(l) - y_start(l) )**2 )
                u_agl(l) = speed_agl(l) * ( x_end(l) - x_start(l) ) / distance
                v_agl(l) = speed_agl(l) * ( y_end(l) - y_start(l) ) / distance
                w_agl(l) = rate_of_climb(l)
!
!--          In case of cyclic-legs, flight direction is directly derived from 
!--          the given flight angle.
             ELSE
                u_agl(l) = speed_agl(l) * COS( flight_angle(l) * pi / 180.0_wp )
                v_agl(l) = speed_agl(l) * SIN( flight_angle(l) * pi / 180.0_wp )
                w_agl(l) = rate_of_climb(l)
             ENDIF

          ENDDO
             
       ENDIF    
!
!--    Initialized data output
       CALL flight_init_output       
!
!--    Allocate array required for user-defined quantities if necessary.
       IF ( num_var_fl_user  > 0 )                                             &
          ALLOCATE( var_u(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
!
!--    Allocate and initialize arrays containing the measured data
       ALLOCATE( sensor_l(1:num_var_fl,1:num_leg) )
       ALLOCATE( sensor(1:num_var_fl,1:num_leg)   )
       sensor_l = 0.0_wp
       sensor   = 0.0_wp

       IF ( debug_output )  CALL debug_message( 'flight_init', 'end' )

    END SUBROUTINE flight_init
    
    
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of output-variable names and units.
!------------------------------------------------------------------------------!
    SUBROUTINE flight_init_output
    
       USE control_parameters,                                                 &
          ONLY:  cloud_droplets, humidity, neutral, passive_scalar

       USE bulk_cloud_model_mod,                                               &
           ONLY:  bulk_cloud_model

       USE netcdf_interface
    
       IMPLICIT NONE

       CHARACTER(LEN=10) ::  label_leg  !< dummy argument to convert integer to string
       
       INTEGER(iwp) ::  i               !< loop variable
       INTEGER(iwp) ::  id_pt           !< identifyer for labeling
       INTEGER(iwp) ::  id_q            !< identifyer for labeling
       INTEGER(iwp) ::  id_ql           !< identifyer for labeling
       INTEGER(iwp) ::  id_s            !< identifyer for labeling        
       INTEGER(iwp) ::  id_u = 1        !< identifyer for labeling  
       INTEGER(iwp) ::  id_v = 2        !< identifyer for labeling
       INTEGER(iwp) ::  id_w = 3        !< identifyer for labeling
       INTEGER(iwp) ::  k               !< index variable
       
       LOGICAL      ::  init = .TRUE.   !< flag to distiquish calls of user_init_flight
!
!--    Define output quanities, at least three variables are measured (u,v,w)
       num_var_fl = 3
       IF ( .NOT. neutral                     )  THEN
          num_var_fl = num_var_fl + 1
          id_pt      = num_var_fl
       ENDIF
       IF ( humidity                          )  THEN
          num_var_fl = num_var_fl + 1
          id_q       = num_var_fl
       ENDIF
       IF ( bulk_cloud_model .OR. cloud_droplets )  THEN
          num_var_fl = num_var_fl + 1  
          id_ql      = num_var_fl
       ENDIF
       IF ( passive_scalar                    )  THEN
          num_var_fl = num_var_fl + 1
          id_s       = num_var_fl
       ENDIF
!
!--    Write output strings for dimensions x, y, z
       DO l=1, num_leg

          IF ( l < 10                    )  WRITE( label_leg, '(I1)')  l
          IF ( l >= 10   .AND.  l < 100  )  WRITE( label_leg, '(I2)')  l
          IF ( l >= 100  .AND.  l < 1000 )  WRITE( label_leg, '(I3)')  l

          dofl_dim_label_x(l)  = 'x_' // TRIM( label_leg )
          dofl_dim_label_y(l)  = 'y_' // TRIM( label_leg )
          dofl_dim_label_z(l)  = 'z_' // TRIM( label_leg )

       ENDDO
       
!
!--    Call user routine to initialize further variables
       CALL user_init_flight( init )
!
!--    Write output labels and units for the quanities
       k = 1
       DO l=1, num_leg
       
          IF ( l < 10                    )  WRITE( label_leg, '(I1)')  l
          IF ( l >= 10   .AND.  l < 100  )  WRITE( label_leg, '(I2)')  l
          IF ( l >= 100  .AND.  l < 1000 )  WRITE( label_leg, '(I3)')  l
          
          label_leg = 'leg_' // TRIM(label_leg) 
          DO i=1, num_var_fl

             IF ( i == id_u      )  THEN         
                dofl_label(k) = TRIM( label_leg ) // '_u'
                dofl_unit(k)  = 'm/s'
                k             = k + 1
             ELSEIF ( i == id_v  )  THEN       
                dofl_label(k) = TRIM( label_leg ) // '_v'
                dofl_unit(k)  = 'm/s'
                k             = k + 1
             ELSEIF ( i == id_w  )  THEN          
                dofl_label(k) = TRIM( label_leg ) // '_w'
                dofl_unit(k)  = 'm/s'
                k             = k + 1
             ELSEIF ( i == id_pt )  THEN       
                dofl_label(k) = TRIM( label_leg ) // '_theta'
                dofl_unit(k)  = 'K'
                k             = k + 1
             ELSEIF ( i == id_q  )  THEN       
                dofl_label(k) = TRIM( label_leg ) // '_q'
                dofl_unit(k)  = 'kg/kg'
                k             = k + 1
             ELSEIF ( i == id_ql )  THEN       
                dofl_label(k) = TRIM( label_leg ) // '_ql'
                dofl_unit(k)  = 'kg/kg'
                k             = k + 1
             ELSEIF ( i == id_s  )  THEN                          
                dofl_label(k) = TRIM( label_leg ) // '_s'
                dofl_unit(k)  = 'kg/kg'
                k             = k + 1
             ENDIF
          ENDDO
          
          DO i=1, num_var_fl_user
             CALL user_init_flight( init, k, i, label_leg )
          ENDDO
          
       ENDDO 
!
!--    Finally, set the total number of flight-output quantities.
       num_var_fl = num_var_fl + num_var_fl_user
       
    END SUBROUTINE flight_init_output

!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine calculates the current flight positions and calls the 
!> respective interpolation routine to measures the data at the current 
!> flight position.
!------------------------------------------------------------------------------!
    SUBROUTINE flight_measurement

       USE arrays_3d,                                                          &
           ONLY:  ddzu, ddzw, pt, q, ql, s, u, v, w, zu, zw

       USE control_parameters,                                                 &
           ONLY:  cloud_droplets, dt_3d, humidity, neutral,                    &
                  passive_scalar, simulated_time
                  
       USE cpulog,                                                             &
           ONLY:  cpu_log, log_point

       USE grid_variables,                                                     &
           ONLY:  ddx, ddy, dx, dy

       USE indices,                                                            &
           ONLY:  nx, nxl, nxr, ny, nys, nyn

       USE bulk_cloud_model_mod,                                               &
           ONLY:  bulk_cloud_model

       USE pegrid

       IMPLICIT NONE

       LOGICAL  ::  on_pe  !< flag to check if current flight position is on current PE

       REAL(wp) ::  x  !< distance between left edge of current grid box and flight position
       REAL(wp) ::  y  !< distance between south edge of current grid box and flight position

       INTEGER(iwp) ::  i   !< index of current grid box along x
       INTEGER(iwp) ::  j   !< index of current grid box along y
       INTEGER(iwp) ::  n   !< loop variable for number of user-defined output quantities
       
       CALL cpu_log( log_point(65), 'flight_measurement', 'start' )
!
!--    Perform flight measurement if start time is reached.
       IF ( simulated_time >= flight_begin  .AND.                              &
            simulated_time <= flight_end )  THEN

          sensor_l = 0.0_wp
          sensor   = 0.0_wp
!
!--       Loop over all flight legs
          DO l=1, num_leg
!
!--          Update location for each leg
             x_pos(l) = x_pos(l) + u_agl(l) * dt_3d 
             y_pos(l) = y_pos(l) + v_agl(l) * dt_3d 
             z_pos(l) = z_pos(l) + w_agl(l) * dt_3d
!
!--          Check if location must be modified for return legs.  
!--          Carry out horizontal reflection if required.
             IF ( .NOT. cyclic_leg(l) )  THEN

                IF ( x_start(l) <= x_end(l) )  THEN
!
!--                Outward flight, i.e. from start to end
                   IF ( u_agl(l) >= 0.0_wp  .AND.  x_pos(l) > x_end(l)      )  THEN
                      x_pos(l) = 2.0_wp * x_end(l)   - x_pos(l)
                      u_agl(l) = - u_agl(l)
!                  
!--                Return flight, i.e. from end to start
                   ELSEIF ( u_agl(l) < 0.0_wp  .AND.  x_pos(l) < x_start(l) )  THEN 
                      x_pos(l) = 2.0_wp * x_start(l) - x_pos(l)
                      u_agl(l) = - u_agl(l)
                   ENDIF
                ELSE
!
!--                Outward flight, i.e. from start to end
                   IF ( u_agl(l) < 0.0_wp  .AND.  x_pos(l) < x_end(l)      )  THEN
                      x_pos(l) = 2.0_wp * x_end(l)   - x_pos(l)
                      u_agl(l) = - u_agl(l)
!                  
!--                Return flight, i.e. from end to start
                   ELSEIF ( u_agl(l) >= 0.0_wp  .AND.  x_pos(l) > x_start(l) )  THEN 
                      x_pos(l) = 2.0_wp * x_start(l) - x_pos(l)
                      u_agl(l) = - u_agl(l)
                   ENDIF
                ENDIF
                
                IF ( y_start(l) <= y_end(l) )  THEN
!
!--                Outward flight, i.e. from start to end
                   IF ( v_agl(l) >= 0.0_wp  .AND.  y_pos(l) > y_end(l)      )  THEN
                      y_pos(l) = 2.0_wp * y_end(l)   - y_pos(l)
                      v_agl(l) = - v_agl(l)
!                  
!--                Return flight, i.e. from end to start                  
                   ELSEIF ( v_agl(l) < 0.0_wp  .AND.  y_pos(l) < y_start(l) )  THEN 
                      y_pos(l) = 2.0_wp * y_start(l) - y_pos(l)
                      v_agl(l) = - v_agl(l)
                   ENDIF
                ELSE
!
!--                Outward flight, i.e. from start to end
                   IF ( v_agl(l) < 0.0_wp  .AND.  y_pos(l) < y_end(l)      )  THEN
                      y_pos(l) = 2.0_wp * y_end(l)   - y_pos(l)
                      v_agl(l) = - v_agl(l)
!                  
!--                Return flight, i.e. from end to start                  
                   ELSEIF ( v_agl(l) >= 0.0_wp  .AND.  y_pos(l) > y_start(l) )  THEN 
                      y_pos(l) = 2.0_wp * y_start(l) - y_pos(l)
                      v_agl(l) = - v_agl(l)
                   ENDIF
                ENDIF
!
!--          Check if flight position is out of the model domain and apply
!--          cyclic conditions if required
             ELSEIF ( cyclic_leg(l) )  THEN
!
!--             Check if aircraft leaves the model domain at the right boundary
                IF ( ( flight_angle(l) >= 0.0_wp     .AND.                     &
                       flight_angle(l) <= 90.0_wp )  .OR.                      &  
                     ( flight_angle(l) >= 270.0_wp   .AND.                     &
                       flight_angle(l) <= 360.0_wp ) )  THEN
                   IF ( x_pos(l) >= ( nx + 0.5_wp ) * dx )                     &
                      x_pos(l) = x_pos(l) - ( nx + 1 ) * dx 
!
!--             Check if aircraft leaves the model domain at the left boundary
                ELSEIF ( flight_angle(l) > 90.0_wp  .AND.                      &
                         flight_angle(l) < 270.0_wp )  THEN
                   IF ( x_pos(l) < -0.5_wp * dx )                             &
                      x_pos(l) = ( nx + 1 ) * dx + x_pos(l)  
                ENDIF 
!
!--             Check if aircraft leaves the model domain at the north boundary
                IF ( flight_angle(l) >= 0.0_wp  .AND.                          &
                     flight_angle(l) <= 180.0_wp )  THEN
                   IF ( y_pos(l) >= ( ny + 0.5_wp ) * dy )                     &
                      y_pos(l) = y_pos(l) - ( ny + 1 ) * dy 
!
!--             Check if aircraft leaves the model domain at the south boundary
                ELSEIF ( flight_angle(l) > 180.0_wp  .AND.                     &
                         flight_angle(l) < 360.0_wp )  THEN
                   IF ( y_pos(l) < -0.5_wp * dy )                              &
                      y_pos(l) = ( ny + 1 ) * dy + y_pos(l) 
                ENDIF
                
             ENDIF
!
!--          Check if maximum elevation change is already reached. If required
!--          reflect vertically.
             IF ( rate_of_climb(l) /= 0.0_wp )  THEN
!
!--             First check if aircraft is too high
                IF (  w_agl(l) > 0.0_wp  .AND.                                 &
                      z_pos(l) - flight_level(l) > max_elev_change(l) )  THEN
                   z_pos(l) = 2.0_wp * ( flight_level(l) + max_elev_change(l) )&
                              - z_pos(l)
                   w_agl(l) = - w_agl(l)
!
!--             Check if aircraft is too low
                ELSEIF (  w_agl(l) < 0.0_wp  .AND.  z_pos(l) < flight_level(l) )  THEN
                   z_pos(l) = 2.0_wp * flight_level(l) - z_pos(l)
                   w_agl(l) = - w_agl(l)
                ENDIF
                
             ENDIF 
!
!--          Determine grid indices for flight position along x- and y-direction.
!--          Please note, there is a special treatment for the index 
!--          along z-direction, which is due to vertical grid stretching. 
             i = ( x_pos(l) + 0.5_wp * dx ) * ddx
             j = ( y_pos(l) + 0.5_wp * dy ) * ddy
!
!--          Check if indices are on current PE
             on_pe = ( i >= nxl  .AND.  i <= nxr  .AND.                        &
                       j >= nys  .AND.  j <= nyn )

             IF ( on_pe )  THEN

                var_index = 1
!
!--             Recalculate indices, as u is shifted by -0.5*dx.
                i =   x_pos(l) * ddx
                j = ( y_pos(l) + 0.5_wp * dy ) * ddy
!
!--             Calculate distance from left and south grid-box edges.
                x  = x_pos(l) - ( 0.5_wp - i ) * dx
                y  = y_pos(l) - j * dy
!
!--             Interpolate u-component onto current flight position.
                CALL interpolate_xyz( u, zu, ddzu, 1.0_wp, x, y, var_index, j, i )
                var_index = var_index + 1
!
!--             Recalculate indices, as v is shifted by -0.5*dy.
                i = ( x_pos(l) + 0.5_wp * dx ) * ddx
                j =   y_pos(l) * ddy

                x  = x_pos(l) - i * dx
                y  = y_pos(l) - ( 0.5_wp - j ) * dy
                CALL interpolate_xyz( v, zu, ddzu, 1.0_wp, x, y, var_index, j, i )
                var_index = var_index + 1
!
!--             Interpolate w and scalar quantities. Recalculate indices.
                i  = ( x_pos(l) + 0.5_wp * dx ) * ddx
                j  = ( y_pos(l) + 0.5_wp * dy ) * ddy
                x  = x_pos(l) - i * dx
                y  = y_pos(l) - j * dy
!
!--             Interpolate w-velocity component.
                CALL interpolate_xyz( w, zw, ddzw, 0.0_wp, x, y, var_index, j, i )
                var_index = var_index + 1
!
!--             Potential temerature
                IF ( .NOT. neutral )  THEN
                   CALL interpolate_xyz( pt, zu, ddzu, 1.0_wp, x, y, var_index, j, i )
                   var_index = var_index + 1
                ENDIF
!
!--             Humidity
                IF ( humidity )  THEN
                   CALL interpolate_xyz( q, zu, ddzu, 1.0_wp, x, y, var_index, j, i )
                   var_index = var_index + 1
                ENDIF
!
!--             Liquid water content
                IF ( bulk_cloud_model .OR. cloud_droplets )  THEN
                   CALL interpolate_xyz( ql, zu, ddzu, 1.0_wp, x, y, var_index, j, i )
                   var_index = var_index + 1
                ENDIF
!
!--             Passive scalar
                IF ( passive_scalar )  THEN
                   CALL interpolate_xyz( s, zu, ddzu, 1.0_wp, x, y, var_index, j, i )
                   var_index = var_index + 1
                ENDIF
!
!--             Treat user-defined variables if required
                DO n = 1, num_var_fl_user                
                   CALL user_flight( var_u, n )
                   CALL interpolate_xyz( var_u, zu, ddzu, 1.0_wp, x, y, var_index, j, i )
                   var_index = var_index + 1
                ENDDO
             ENDIF

          ENDDO
!
!--       Write local data on global array.
#if defined( __parallel )
          CALL MPI_ALLREDUCE(sensor_l(1,1), sensor(1,1),                       &
                             num_var_fl*num_leg, MPI_REAL, MPI_SUM,               &
                             comm2d, ierr )
#else
          sensor     = sensor_l
#endif
       ENDIF
       
       CALL cpu_log( log_point(65), 'flight_measurement', 'stop' )

    END SUBROUTINE flight_measurement

!------------------------------------------------------------------------------!
! Description:
! ------------
!> This subroutine bi-linearly interpolates the respective data onto the current
!> flight position.
!------------------------------------------------------------------------------!
    SUBROUTINE interpolate_xyz( var, z_uw, ddz_uw, fac, x, y, var_ind, j, i )

       USE control_parameters,                                                 &
           ONLY:  dz, dz_stretch_level_start   
  
      USE grid_variables,                                                      &
          ONLY:  dx, dy
    
       USE indices,                                                            &
           ONLY:  nzb, nzt, nxlg, nxrg, nysg, nyng

       IMPLICIT NONE

       INTEGER(iwp) ::  i        !< index along x
       INTEGER(iwp) ::  j        !< index along y
       INTEGER(iwp) ::  k        !< index along z
       INTEGER(iwp) ::  k1       !< dummy variable
       INTEGER(iwp) ::  var_ind  !< index variable for output quantity

       REAL(wp) ::  aa        !< dummy argument for horizontal interpolation   
       REAL(wp) ::  bb        !< dummy argument for horizontal interpolation
       REAL(wp) ::  cc        !< dummy argument for horizontal interpolation
       REAL(wp) ::  dd        !< dummy argument for horizontal interpolation
       REAL(wp) ::  gg        !< dummy argument for horizontal interpolation
       REAL(wp) ::  fac       !< flag to indentify if quantity is on zu or zw level
       REAL(wp) ::  var_int   !< horizontal interpolated variable at current position
       REAL(wp) ::  var_int_l !< horizontal interpolated variable at k-level
       REAL(wp) ::  var_int_u !< horizontal interpolated variable at (k+1)-level
       REAL(wp) ::  x         !< distance between left edge of current grid box and flight position
       REAL(wp) ::  y         !< distance between south edge of current grid box and flight position

       REAL(wp), DIMENSION(1:nzt+1)   ::  ddz_uw !< inverse vertical grid spacing
       REAL(wp), DIMENSION(nzb:nzt+1) ::  z_uw   !< height level

       REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  var !< treted quantity
!
!--    Calculate interpolation coefficients
       aa = x**2          + y**2
       bb = ( dx - x )**2 + y**2
       cc = x**2          + ( dy - y )**2
       dd = ( dx - x )**2 + ( dy - y )**2
       gg = aa + bb + cc + dd
!
!--    Obtain vertical index. Special treatment for grid index along z-direction
!--    if flight position is above the vertical grid-stretching level.
!--    fac=1 if variable is on scalar grid level, fac=0 for w-component.
       IF ( z_pos(l) < dz_stretch_level_start(1) )  THEN
          k = ( z_pos(l) + fac * 0.5_wp * dz(1) ) / dz(1)
       ELSE
!
!--       Search for k-index
          DO k1=nzb, nzt
             IF ( z_pos(l) >= z_uw(k1) .AND. z_pos(l) < z_uw(k1+1) )  THEN
                k = k1
                EXIT
             ENDIF                   
          ENDDO
       ENDIF
!
!--    (x,y)-interpolation at k-level
       var_int_l = ( ( gg - aa ) * var(k,j,i)       +                          &
                     ( gg - bb ) * var(k,j,i+1)     +                          &
                     ( gg - cc ) * var(k,j+1,i)     +                          &
                     ( gg - dd ) * var(k,j+1,i+1)                              &
                   ) / ( 3.0_wp * gg )
!
!--    (x,y)-interpolation on (k+1)-level
       var_int_u = ( ( gg - aa ) * var(k+1,j,i)     +                          &
                     ( gg - bb ) * var(k+1,j,i+1)   +                          &
                     ( gg - cc ) * var(k+1,j+1,i)   +                          &
                     ( gg - dd ) * var(k+1,j+1,i+1)                            &
                   ) / ( 3.0_wp * gg )
!
!--    z-interpolation onto current flight postion
       var_int = var_int_l                                                     &
                           + ( z_pos(l) - z_uw(k) ) * ddz_uw(k+1)              &
                           * (var_int_u - var_int_l )
!
!--    Store on local data array
       sensor_l(var_ind,l) = var_int

    END SUBROUTINE interpolate_xyz

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Perform parameter checks.
!------------------------------------------------------------------------------!
    SUBROUTINE flight_check_parameters

       USE arrays_3d,                                                          &
           ONLY:  zu
    
       USE control_parameters,                                                 &
           ONLY:  bc_lr_cyc, bc_ns_cyc, message_string

       USE grid_variables,                                                     &
           ONLY:  dx, dy   
          
       USE indices,                                                            &
           ONLY:  nx, ny, nz
           
       USE netcdf_interface,                                                   &
           ONLY:  netcdf_data_format

       IMPLICIT NONE
       
!
!--    Check if start positions are properly set.
       DO l=1, num_leg
          IF ( x_start(l) < 0.0_wp  .OR.  x_start(l) > ( nx + 1 ) * dx )  THEN
             message_string = 'Start x position is outside the model domain'
             CALL message( 'flight_check_parameters', 'PA0431', 1, 2, 0, 6, 0 )
          ENDIF 
          IF ( y_start(l) < 0.0_wp  .OR.  y_start(l) > ( ny + 1 ) * dy )  THEN
             message_string = 'Start y position is outside the model domain'
             CALL message( 'flight_check_parameters', 'PA0432', 1, 2, 0, 6, 0 )
          ENDIF 

       ENDDO
!
!--    Check for leg mode
       DO l=1, num_leg
!
!--       Check if leg mode matches the overall lateral model boundary 
!--       conditions.
          IF ( TRIM( leg_mode(l) ) == 'cyclic' )  THEN
             IF ( .NOT. bc_lr_cyc  .OR.  .NOT. bc_ns_cyc )  THEN
                message_string = 'Cyclic flight leg does not match ' //        &
                                 'lateral boundary condition'
                CALL message( 'flight_check_parameters', 'PA0433', 1, 2, 0, 6, 0 )
             ENDIF 
!
!--       Check if end-positions are inside the model domain in case of 
!..       return-legs.  
          ELSEIF ( TRIM( leg_mode(l) ) == 'return' )  THEN
             IF ( x_end(l) > ( nx + 1 ) * dx  .OR.                             &
                  y_end(l) > ( ny + 1 ) * dx )  THEN
                message_string = 'Flight leg or parts of it are outside ' //   &
                                 'the model domain'
                CALL message( 'flight_check_parameters', 'PA0434', 1, 2, 0, 6, 0 )
             ENDIF 
          ELSE
             message_string = 'Unknown flight mode'
             CALL message( 'flight_check_parameters', 'PA0435', 1, 2, 0, 6, 0 )
          ENDIF

       ENDDO          
!
!--    Check if given flight object remains inside model domain if a rate of 
!--    climb / descent is prescribed. 
       DO l=1, num_leg
          IF ( flight_level(l) + max_elev_change(l) > zu(nz) .OR.              &
               flight_level(l) + max_elev_change(l) <= 0.0_wp )  THEN
             message_string = 'Flight level is outside the model domain '
             CALL message( 'flight_check_parameters', 'PA0438', 1, 2, 0, 6, 0 )
          ENDIF 
       ENDDO       
!
!--    Check for appropriate NetCDF format. Definition of more than one 
!--    unlimited dimension is unfortunately only possible with NetCDF4/HDF5. 
       IF (  netcdf_data_format <= 2 )  THEN
          message_string = 'netcdf_data_format must be > 2'
          CALL message( 'flight_check_parameters', 'PA0439', 1, 2, 0, 6, 0 )
       ENDIF


    END SUBROUTINE flight_check_parameters
        

!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine reads the respective restart data.
!------------------------------------------------------------------------------!
    SUBROUTINE flight_rrd_global( found )  


       USE control_parameters,                                                 &
           ONLY: length, restart_string

    
       IMPLICIT NONE
       
       LOGICAL, INTENT(OUT)  ::  found !< flag indicating if a variable string is found


       found = .TRUE.


       SELECT CASE ( restart_string(1:length) )
          
          CASE ( 'u_agl' )
             IF ( .NOT. ALLOCATED( u_agl ) )  ALLOCATE( u_agl(1:num_leg) )
             READ ( 13 )  u_agl   
          CASE ( 'v_agl' )
             IF ( .NOT. ALLOCATED( v_agl ) )  ALLOCATE( v_agl(1:num_leg) )
             READ ( 13 )  v_agl
          CASE ( 'w_agl' )
             IF ( .NOT. ALLOCATED( w_agl ) )  ALLOCATE( w_agl(1:num_leg) )
             READ ( 13 )  w_agl
          CASE ( 'x_pos' )
             IF ( .NOT. ALLOCATED( x_pos ) )  ALLOCATE( x_pos(1:num_leg) )
             READ ( 13 )  x_pos
          CASE ( 'y_pos' )
             IF ( .NOT. ALLOCATED( y_pos ) )  ALLOCATE( y_pos(1:num_leg) )
             READ ( 13 )  y_pos
          CASE ( 'z_pos' )
             IF ( .NOT. ALLOCATED( z_pos ) )  ALLOCATE( z_pos(1:num_leg) )
             READ ( 13 )  z_pos

          CASE DEFAULT

             found = .FALSE.
          
       END SELECT


    END SUBROUTINE flight_rrd_global  
    
!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine writes the respective restart data.
!------------------------------------------------------------------------------!
    SUBROUTINE flight_wrd_global  


       IMPLICIT NONE
 
       CALL wrd_write_string( 'u_agl' )
       WRITE ( 14 )  u_agl

       CALL wrd_write_string( 'v_agl' )

       WRITE ( 14 )  v_agl

       CALL wrd_write_string( 'w_agl' )
       WRITE ( 14 )  w_agl

       CALL wrd_write_string( 'x_pos' )
       WRITE ( 14 )  x_pos

       CALL wrd_write_string( 'y_pos' )
       WRITE ( 14 )  y_pos

       CALL wrd_write_string( 'z_pos' )
       WRITE ( 14 )  z_pos
       
    END SUBROUTINE flight_wrd_global    
    

 END MODULE flight_mod
