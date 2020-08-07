!> @file progress_bar_mod.f90
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
! $Id: progress_bar_mod.f90 4360 2020-01-07 11:25:50Z suehring $
! Corrected "Former revisions" section
! 
! 3655 2019-01-07 16:51:22Z knoop
! Increase printed length of run identifier,
! bugfix for restarts
! 1468 2014-09-24 14:06:57Z maronga
! Added support for progress file PROGRESS which is used in case of batch jobs
!
! Description:
! ------------
!> This routine prints either a progress bar on the standard output in case of
!> interactive runs, or it prints the progress in a separate file called
!> PROGRESS.
!------------------------------------------------------------------------------!
 MODULE progress_bar
 

    USE control_parameters,                                                    &
        ONLY : end_time, initializing_actions, run_identifier, simulated_time, &
               simulated_time_at_begin, spinup_time, time_restart

    USE, INTRINSIC ::  ISO_FORTRAN_ENV,                                        &
        ONLY :  OUTPUT_UNIT

    USE kinds

    IMPLICIT NONE

    PRIVATE
    PUBLIC   progress_bar_disabled, finish_progress_bar, output_progress_bar

    CHARACTER(LEN=60) ::  bar      !< progress bar, initially filled with "_"
    CHARACTER(LEN=60) ::  crosses  !< filled with "X"

    INTEGER(iwp) ::  ilength !< length of progress bar filled with "X"

    LOGICAL ::  progress_bar_disabled = .FALSE.   !< envpar-Namelist switch

    REAL(wp) ::  time_to_be_simulated !< in sec

    LOGICAL ::  initialized = .FALSE. !< switch to determine if bar is initialized

    SAVE

 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize the progress bar/file
!------------------------------------------------------------------------------!
 
    SUBROUTINE init_progress_bar

       IMPLICIT NONE

!
!--    Calculate the time to be simulated in this job
!--    (in case of automatic restarts the calculated time will probably be
!--    larger than the time which will actually be simulated)
       IF ( time_restart /= 9999999.9_wp  .AND.  time_restart < end_time  .AND.&
            time_restart > simulated_time_at_begin )  THEN
          time_to_be_simulated = time_restart - simulated_time_at_begin
       ELSEIF ( initializing_actions == 'read_restart_data' )  THEN
          time_to_be_simulated = end_time     - simulated_time_at_begin
       ELSE
          time_to_be_simulated = end_time     - simulated_time_at_begin        &
                                              - spinup_time
       ENDIF

       IF ( progress_bar_disabled )  THEN

          CALL check_open ( 117 )
          WRITE ( 117, FMT='(A34,/)' ) run_identifier

       ELSE
          bar = '____________________________________________________________'
          crosses = 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
!
!--       Line feed on stdout to seperate the progress bar from previous messages
          WRITE ( OUTPUT_UNIT, '(1X)' )
#if defined( __intel_compiler )
!
!--       The Intel compiler does not allow to immediately flush the output buffer
!--       in case that option ADVANCE='NO' is used in the write statement.
!--       A workaround is to set a special carriage control feature and use "+" as
!--       first output character, but this non-standard and only available with the
!--       Intel compiler
          OPEN ( OUTPUT_UNIT, CARRIAGECONTROL='FORTRAN' )
#endif

       ENDIF

       initialized = .TRUE.

    END SUBROUTINE init_progress_bar


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Print progress data to standard output (interactive) or to file (batch jobs)
!------------------------------------------------------------------------------!
 
    SUBROUTINE output_progress_bar

       IMPLICIT NONE

       REAL(wp) ::  remaining_time_in_percent  !< remaining time to be simulated
                                               !< in the job
       REAL(wp) ::  remaining_time_in_percent_total !< total remaining time of
                                                    !< the job chain

       IF ( .NOT. initialized )  CALL init_progress_bar

       IF ( initializing_actions == 'read_restart_data' )  THEN
          remaining_time_in_percent =                                             &
             ( simulated_time - simulated_time_at_begin )                         &
             / time_to_be_simulated

          remaining_time_in_percent_total = simulated_time / end_time
       ELSE
          remaining_time_in_percent =                                             &
             ( simulated_time - simulated_time_at_begin - spinup_time )           &
             / time_to_be_simulated

          remaining_time_in_percent_total = ( ( simulated_time - spinup_time )    &
                                            / ( end_time       - spinup_time ) )
       ENDIF

!
!--    In batch mode, use a file (PROGRESS), otherwise use progress bar
       IF ( progress_bar_disabled )  THEN

          BACKSPACE ( 117 )
          WRITE ( 117, FMT='(F5.2,1X,F5.2)' ) remaining_time_in_percent,       &
                                              remaining_time_in_percent_total
          FLUSH( 117 )

       ELSE

!
!--       Calculate length of progress bar
          ilength = remaining_time_in_percent * 60.0_wp
          ilength = MIN( ilength, 60 )

          bar(1:ilength) = crosses(1:ilength)

#if defined( __intel_compiler )
          WRITE ( OUTPUT_UNIT, '(A,6X,''['',A,''] '',F5.1,'' left'')' )        &
                  '+', bar,                                                    &
                   MAX( 0.0_wp, ( 1.0_wp - remaining_time_in_percent ) *       &
                                  100.0_wp )
#else
          WRITE ( OUTPUT_UNIT, '(A,6X,''['',A,''] '',F5.1,'' left'')',         &
                  ADVANCE='NO' )  CHAR( 13 ), bar,                             &
                   MAX( 0.0_wp, ( 1.0_wp - remaining_time_in_percent ) *       &
                                  100.0_wp )
#endif
          FLUSH( OUTPUT_UNIT )

       ENDIF

    END SUBROUTINE output_progress_bar

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Finalization of the progress bar/file
!------------------------------------------------------------------------------!
 
    SUBROUTINE finish_progress_bar

       IMPLICIT NONE

       IF ( progress_bar_disabled )  THEN

          CALL close_file ( 117 )

       ELSE
       
#if defined( __intel_compiler )
!
!--       Reset to the default carriage control
          OPEN ( OUTPUT_UNIT, CARRIAGECONTROL='LIST' )
#endif
!
!--       Line feed when simulation has finished
          WRITE ( OUTPUT_UNIT, '(1X)' )

       ENDIF

    END SUBROUTINE finish_progress_bar


 END MODULE progress_bar
