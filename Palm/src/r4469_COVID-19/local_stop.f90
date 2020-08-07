!> @file local_stop.f90
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
! $Id: local_stop.f90 4444 2020-03-05 15:59:50Z raasch $
! bugfix: misplaced cpp-directive moved
! 
! 4360 2020-01-07 11:25:50Z suehring
! Corrected "Former revisions" section
! 
! 3655 2019-01-07 16:51:22Z knoop
! Added an empty output string to stop keywords to clean up job protocol
!
! Revision 1.1  2002/12/19 15:46:23  raasch
! Initial revision
!
!
! Description:
! ------------
!> Stop program execution
!------------------------------------------------------------------------------!
 SUBROUTINE local_stop
 
#if defined( __parallel )

    USE control_parameters,                                                    &
        ONLY:  abort_mode, coupling_mode, coupling_mode_remote, dt_restart,    &
               stop_dt, terminate_coupled, terminate_coupled_remote,           &
               terminate_run, time_restart

    USE pegrid

    USE pmc_interface,                                                         &
        ONLY:  nested_run


    IF ( coupling_mode == 'uncoupled' )  THEN
       IF ( nested_run )  THEN
!
!--       Workaround: If any of the nested model crashes, it aborts the whole
!--       run with MPI_ABORT, regardless of the reason given by abort_mode
          CALL MPI_ABORT( MPI_COMM_WORLD, 9999, ierr )
       ELSE
          IF ( abort_mode == 1 )  THEN
             CALL MPI_FINALIZE( ierr )
             STOP 1
          ELSEIF ( abort_mode == 2 )  THEN
             CALL MPI_ABORT( comm2d, 9999, ierr )
          ELSEIF ( abort_mode == 3 )  THEN
             CALL MPI_ABORT( MPI_COMM_WORLD, 9999, ierr )
          ENDIF
       ENDIF
    ELSEIF ( coupling_mode(1:8) == 'vnested_' )  THEN

       PRINT*, '+++ local_stop:'
       PRINT*, '     model "', TRIM( coupling_mode ), '" terminated'
!
!--    Abort both coarse and fine grid
       CALL MPI_ABORT( MPI_COMM_WORLD, 9999, ierr )
    ELSE

       SELECT CASE ( terminate_coupled_remote )

          CASE ( 0 )
             IF ( myid == 0 )  THEN
                PRINT*, '+++ local_stop:'
                PRINT*, '    local model "', TRIM( coupling_mode ), &
                     '" stops now'
             ENDIF
!
!--          Inform the remote model of the termination and its reason, provided
!--          the remote model has not already been informed of another 
!--          termination reason (terminate_coupled > 0) before.
             IF ( terminate_coupled == 0 )  THEN
                terminate_coupled = 1
                IF ( myid == 0 ) THEN
                   CALL MPI_SENDRECV( &
                        terminate_coupled,        1, MPI_INTEGER, target_id,  0, &
                        terminate_coupled_remote, 1, MPI_INTEGER, target_id,  0, &
                        comm_inter, status, ierr )
                ENDIF
                CALL MPI_BCAST( terminate_coupled_remote, 1, MPI_REAL, 0, comm2d, ierr)
             ENDIF
             CALL MPI_FINALIZE( ierr )
             STOP 1

          CASE ( 1 )
             IF ( myid == 0 )  THEN
                PRINT*, '+++ local_stop:'
                PRINT*, '    remote model "', TRIM( coupling_mode_remote ), &
                     '" stopped'
             ENDIF
             CALL MPI_FINALIZE( ierr )
             STOP 1

          CASE ( 2 )
             IF ( myid == 0 )  THEN
                PRINT*, '+++ local_stop:'
                PRINT*, '    remote model "', TRIM( coupling_mode_remote ), &
                     '" terminated'
                PRINT*, '    with stop_dt = .T.'
             ENDIF
             stop_dt = .TRUE.

          CASE ( 3 )
             IF ( myid == 0 )  THEN
                PRINT*, '+++ local_stop:'
                PRINT*, '    remote model "', TRIM( coupling_mode_remote ), &
                     '" terminated'
                PRINT*, '    with terminate_run = .T. (CPU-time limit)'
             ENDIF
             terminate_run = .TRUE.

          CASE ( 4 )
             IF ( myid == 0 )  THEN
                PRINT*, '+++ local_stop:'
                PRINT*, '    remote model "', TRIM( coupling_mode_remote ), &
                     '" terminated'
                PRINT*, '    with terminate_run = .T. (restart)'
             ENDIF
             terminate_run = .TRUE.
             time_restart = time_restart + dt_restart

          CASE ( 5 )
             IF ( myid == 0 )  THEN
                PRINT*, '+++ local_stop:'
                PRINT*, '    remote model "', TRIM( coupling_mode_remote ), &
                     '" terminated'
                PRINT*, '    with terminate_run = .T. (single restart)'
             ENDIF
             terminate_run = .TRUE.
             time_restart = 9999999.9_wp

       END SELECT

    ENDIF

#else

    STOP 1

#endif

 END SUBROUTINE local_stop    
