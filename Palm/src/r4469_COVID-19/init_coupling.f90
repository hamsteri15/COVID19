!> @file init_coupling.f90
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
! ------------------
! $Id: init_coupling.f90 4444 2020-03-05 15:59:50Z raasch $
! bugfix: cpp-directives for serial mode added
! 
! 4360 2020-01-07 11:25:50Z suehring
! Corrected "Former revisions" section
! 
! 3655 2019-01-07 16:51:22Z knoop
! references to mrun replaced by palmrun, and updated
!
! 222 2009-01-12 16:04:16Z letzel
! Initial revision
!
! Description:
! ------------
!> Initializing coupling via MPI-1 or MPI-2 if the coupled version of PALM is
!> called.
!------------------------------------------------------------------------------!
  SUBROUTINE init_coupling
 

    USE control_parameters,                                                    &
        ONLY:  coupling_char, coupling_mode
        
    USE kinds
    
    USE pegrid

    USE vertical_nesting_mod

    IMPLICIT NONE

!
!-- Local variables
#if defined( __parallel )
    INTEGER(iwp) ::  i            !<
    INTEGER(iwp) ::  inter_color  !<
#endif
    
    INTEGER(iwp), DIMENSION(:) ::  bc_data(0:3) = 0  !<

!
!-- Get information about the coupling mode from the environment variable
!-- which has been set by the mpiexec command.
!-- This method is currently not used because the mpiexec command is not
!-- available on some machines
!    CALL GET_ENVIRONMENT_VARIABLE( 'coupling_mode', coupling_mode, i )
!    IF ( i == 0 )  coupling_mode = 'uncoupled'
!    IF ( coupling_mode == 'ocean_to_atmosphere' )  coupling_char = '_O'

!
!-- Get information about the coupling mode from standard input (PE0 only) and
!-- distribute it to the other PEs. Distribute PEs to 2 new communicators.
!-- ATTENTION: numprocs will be reset according to the new communicators
#if defined ( __parallel )

    IF ( myid == 0 )  THEN
       READ (*,*,ERR=10,END=10)  coupling_mode, bc_data(1), bc_data(2)
10     CONTINUE
       IF ( TRIM( coupling_mode ) == 'coupled_run' )  THEN
          i = 1
       ELSEIF ( TRIM( coupling_mode ) == 'vnested_twi' )  THEN
          i = 9
       ELSE
          i = 0
       ENDIF
       bc_data(0) = i

!
!--    Check if '_O' has to be used as file extension in an uncoupled ocean
!--    run. This is required, if this run shall be continued as a coupled run.
       IF ( TRIM( coupling_mode ) == 'precursor_ocean' )  bc_data(3) = 1

    ENDIF

    CALL MPI_BCAST( bc_data(0), 4, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
    i = bc_data(0)

    IF ( i == 0 )  THEN
       coupling_mode = 'uncoupled'
!
!--    In case of a precursor ocean run, an additional flag file is created.
!--    This is necessary for data_output_2d_on_each_pe = .T.
       IF ( bc_data(3) == 1 )  THEN
          OPEN( 90, FILE='PRECURSOR_OCEAN', FORM='FORMATTED' )
          WRITE ( 90, '(''TRUE'')' )
          CLOSE ( 90 )
       ENDIF
    ELSEIF ( i == 9 )  THEN

!
!--    Set a flag to identify runs with vertical nesting
       vnested = .TRUE.
       
       comm_inter = MPI_COMM_WORLD
       
!
!--    Split the total available PE's into two groups
!--    numprocs for coarse and fine grid are read from stdin (see above, and
!--    execution command in the palmrun script, numprocs are provided via
!--    palmrun option -Y)
       IF ( myid < bc_data(1) )  THEN
          inter_color     = 0
          numprocs        = bc_data(1)
          coupling_mode   = 'vnested_crse'
       ELSE
          inter_color     = 1
          numprocs        = bc_data(2)
          coupling_mode   = 'vnested_fine'
       ENDIF
       
       CALL MPI_COMM_SPLIT( MPI_COMM_WORLD, inter_color, 0, comm_palm, ierr )
       comm2d = comm_palm
       
       OPEN( 90, FILE='VNESTING_PORT_OPENED', FORM='FORMATTED' )
       WRITE ( 90, '(''TRUE'')' )
       CLOSE ( 90 )
      
    ELSE
       comm_inter = MPI_COMM_WORLD

       IF ( myid < bc_data(1) ) THEN
          inter_color     = 0
          numprocs        = bc_data(1)
          coupling_mode   = 'atmosphere_to_ocean'
       ELSE
          inter_color     = 1
          numprocs        = bc_data(2)
          coupling_mode   = 'ocean_to_atmosphere'
       ENDIF

       CALL MPI_COMM_SPLIT( MPI_COMM_WORLD, inter_color, 0, comm_palm, ierr )
       comm2d = comm_palm

!
!--    Write a flag file for the ocean model and the other atmosphere
!--    processes.
       OPEN( 90, FILE='COUPLING_PORT_OPENED', FORM='FORMATTED' )
       WRITE ( 90, '(''TRUE'')' )
       CLOSE ( 90 )
    ENDIF
#endif

!
!-- In case of a precursor ocean run (followed by a coupled run), or a
!-- coupled atmosphere-ocean run, set the file extension for the ocean files
    IF ( TRIM( coupling_mode ) == 'ocean_to_atmosphere' .OR. bc_data(3) == 1 ) &
    THEN
       coupling_char = '_O'
    ENDIF

    IF (  TRIM( coupling_mode ) == 'vnested_fine' )  THEN
!
!-- Set file extension for vertical nesting
       coupling_char = '_NV'
    ENDIF

 END SUBROUTINE init_coupling
