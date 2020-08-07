!> @file exchange_horiz.f90
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
! $Id: exchange_horiz_mod.f90 4461 2020-03-12 16:51:59Z raasch $
! optional communicator added to exchange_horiz
! 
! 4457 2020-03-11 14:20:43Z raasch
! routine has been modularized, file exchange_horiz_2d has been merged
! 
! 4429 2020-02-27 15:24:30Z raasch
! bugfix: cpp-directives added for serial mode
! 
! 4360 2020-01-07 11:25:50Z suehring
! Corrected "Former revisions" section
! 
! 3761 2019-02-25 15:31:42Z raasch
! OpenACC directives re-formatted
! 
! 3657 2019-01-07 20:14:18Z knoop
! OpenACC port for SPEC
!
! Revision 1.1  1997/07/24 11:13:29  raasch
! Initial revision
!
!
! Description:
! ------------
!> Exchange of ghost point layers for subdomains (in parallel mode) and setting
!> of cyclic lateral boundary conditions for the total domain .
!------------------------------------------------------------------------------!
 MODULE exchange_horiz_mod

    USE kinds

    USE pegrid

    IMPLICIT NONE

    PRIVATE
    PUBLIC exchange_horiz, exchange_horiz_int, exchange_horiz_2d, exchange_horiz_2d_byte,          &
           exchange_horiz_2d_int

    INTERFACE exchange_horiz
       MODULE PROCEDURE exchange_horiz
    END INTERFACE exchange_horiz

    INTERFACE exchange_horiz_int
       MODULE PROCEDURE exchange_horiz_int
    END INTERFACE exchange_horiz_int

    INTERFACE exchange_horiz_2d
       MODULE PROCEDURE exchange_horiz_2d
    END INTERFACE exchange_horiz_2d

    INTERFACE exchange_horiz_2d_byte
       MODULE PROCEDURE exchange_horiz_2d_byte
    END INTERFACE exchange_horiz_2d_byte

    INTERFACE exchange_horiz_2d_int
       MODULE PROCEDURE exchange_horiz_2d_int
    END INTERFACE exchange_horiz_2d_int


 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Exchange of ghost point layers for subdomains (in parallel mode) and setting
!> of cyclic lateral boundary conditions for the total domain.
!> This routine is for REAL 3d-arrays.
!------------------------------------------------------------------------------!
 SUBROUTINE exchange_horiz( ar, nbgp_local, alternative_communicator)

    USE control_parameters,                                                    &
        ONLY:  bc_lr_cyc, bc_ns_cyc

#if defined( __parallel )
    USE control_parameters,                                                    &
        ONLY:  grid_level, mg_switch_to_pe0, synchronous_exchange
#endif
                
    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point_s
        
    USE indices,                                                               &
        ONLY:  nxl, nxr, nyn, nys, nzb, nzt
        

#if defined( _OPENACC )
    INTEGER(iwp) ::  i           !<
#endif

    INTEGER(iwp), OPTIONAL ::  alternative_communicator  !< alternative MPI communicator to be used
    INTEGER(iwp) ::  communicator  !< communicator that is used as argument in MPI calls
    INTEGER(iwp) ::  left_pe       !< id of left pe that is used as argument in MPI calls
    INTEGER(iwp) ::  nbgp_local    !< number of ghost point layers
    INTEGER(iwp) ::  north_pe      !< id of north pe that is used as argument in MPI calls
    INTEGER(iwp) ::  right_pe      !< id of right pe that is used as argument in MPI calls
    INTEGER(iwp) ::  south_pe      !< id of south pe that is used as argument in MPI calls
    
    REAL(wp), DIMENSION(nzb:nzt+1,nys-nbgp_local:nyn+nbgp_local,               &
                        nxl-nbgp_local:nxr+nbgp_local) ::  ar !< 3d-array for which exchange is done
                        

    CALL cpu_log( log_point_s(2), 'exchange_horiz', 'start' )

#if defined( _OPENACC )
    !$ACC UPDATE IF_PRESENT ASYNC(1) &
    !$ACC HOST(ar(:,:,nxr-nbgp_local+1:nxr)) &
    !$ACC HOST(ar(:,:,nxl:nxl+nbgp_local-1))

!
!-- Wait for first UPDATE to complete before starting the others.
    !$ACC WAIT(1) ASYNC(2)
    ! ar(:,:,nxl-nbgp_local:nxl-1) is overwritten by first part below
    ! ar(:,:,nxl:nxl+nbgp_local-1) has been transferred above
    DO i = nxl+nbgp_local, nxr-nbgp_local
       !$ACC UPDATE IF_PRESENT ASYNC(2) &
       !$ACC HOST(ar(:,nyn-nbgp_local+1:nyn,i)) &
       !$ACC HOST(ar(:,nys:nys+nbgp_local-1,i))
    ENDDO
    ! ar(:,:,nxr-nbgp_local+1:nxr) has been transferred above
    ! ar(:,:,nxr+1:nxr+nbgp_local) is overwritten by first part below

!
!-- Wait for first UPDATE to complete before starting MPI.
    !$ACC WAIT(1)
#endif

!
!-- Set the communicator to be used
    IF ( PRESENT( alternative_communicator ) )  THEN
!
!--    Alternative communicator is to be used
       communicator = communicator_configurations(alternative_communicator)%mpi_communicator
       left_pe  = communicator_configurations(alternative_communicator)%pleft
       right_pe = communicator_configurations(alternative_communicator)%pright
       south_pe = communicator_configurations(alternative_communicator)%psouth
       north_pe = communicator_configurations(alternative_communicator)%pnorth

    ELSE
!
!--    Main communicator is to be used
       communicator = comm2d
       left_pe  = pleft
       right_pe = pright
       south_pe = psouth
       north_pe = pnorth

    ENDIF

#if defined( __parallel )

!
!-- Exchange in x-direction of lateral boundaries
    IF ( pdims(1) == 1  .OR.  mg_switch_to_pe0 )  THEN
!
!--    One-dimensional decomposition along y, boundary values can be exchanged
!--    within the PE memory
       IF ( bc_lr_cyc )  THEN
          ar(:,:,nxl-nbgp_local:nxl-1) = ar(:,:,nxr-nbgp_local+1:nxr)
          ar(:,:,nxr+1:nxr+nbgp_local) = ar(:,:,nxl:nxl+nbgp_local-1)
       ENDIF

    ELSE

       IF ( synchronous_exchange )  THEN
!
!--       Send left boundary, receive right one (synchronous)
          CALL MPI_SENDRECV( ar(nzb,nys-nbgp_local,nxl),   1, type_yz(grid_level), left_pe,  0,    &
                             ar(nzb,nys-nbgp_local,nxr+1), 1, type_yz(grid_level), right_pe, 0,    &
                             communicator, status, ierr )
!
!--       Send right boundary, receive left one (synchronous)
          CALL MPI_SENDRECV( ar(nzb,nys-nbgp_local,nxr+1-nbgp_local), 1,                           &
                             type_yz(grid_level), right_pe, 1,                                     &
                             ar(nzb,nys-nbgp_local,nxl-nbgp_local), 1,                             &
                             type_yz(grid_level), left_pe,  1,                                     &
                             communicator, status, ierr )

       ELSE

!
!--       Asynchroneous exchange
          IF ( send_receive == 'lr'  .OR.  send_receive == 'al' )  THEN

             req(1:4)  = 0
             req_count = 0
!
!--          Send left boundary, receive right one (asynchronous)
             CALL MPI_ISEND( ar(nzb,nys-nbgp_local,nxl),   1, type_yz(grid_level), left_pe,        &
                             req_count, communicator, req(req_count+1), ierr )
             CALL MPI_IRECV( ar(nzb,nys-nbgp_local,nxr+1), 1, type_yz(grid_level), right_pe,       &
                             req_count, communicator, req(req_count+2), ierr )
!
!--          Send right boundary, receive left one (asynchronous)
             CALL MPI_ISEND( ar(nzb,nys-nbgp_local,nxr+1-nbgp_local), 1, type_yz(grid_level),      &
                             right_pe, req_count+1, communicator, req(req_count+3), ierr )
             CALL MPI_IRECV( ar(nzb,nys-nbgp_local,nxl-nbgp_local),   1, type_yz(grid_level),      &
                             left_pe,  req_count+1, communicator, req(req_count+4), ierr )

             CALL MPI_WAITALL( 4, req, wait_stat, ierr )

          ENDIF

       ENDIF

    ENDIF

    !$ACC UPDATE IF_PRESENT ASYNC(1) &
    !$ACC DEVICE(ar(:,:,nxl-nbgp_local:nxl-1)) &
    !$ACC DEVICE(ar(:,:,nxr+1:nxr+nbgp_local))

!
!-- Wait for UPDATES above to complete before starting MPI.
    !$ACC WAIT(2)

    IF ( pdims(2) == 1  .OR.  mg_switch_to_pe0 )  THEN
!
!--    One-dimensional decomposition along x, boundary values can be exchanged
!--    within the PE memory
       IF ( bc_ns_cyc )  THEN
          ar(:,nys-nbgp_local:nys-1,:) = ar(:,nyn-nbgp_local+1:nyn,:)
          ar(:,nyn+1:nyn+nbgp_local,:) = ar(:,nys:nys+nbgp_local-1,:)
       ENDIF

    ELSE

       IF ( synchronous_exchange )  THEN
!
!--       Send front boundary, receive rear one (synchronous)
          CALL MPI_SENDRECV( ar(nzb,nys,nxl-nbgp_local),   1, type_xz(grid_level), south_pe, 0,    &
                             ar(nzb,nyn+1,nxl-nbgp_local), 1, type_xz(grid_level), north_pe, 0,    &
                             communicator, status, ierr )
!
!--       Send rear boundary, receive front one (synchronous)
          CALL MPI_SENDRECV( ar(nzb,nyn-nbgp_local+1,nxl-nbgp_local), 1,                           &
                             type_xz(grid_level), north_pe, 1,                                     &
                             ar(nzb,nys-nbgp_local,nxl-nbgp_local),   1,                           &
                             type_xz(grid_level), south_pe, 1,                                     &
                             communicator, status, ierr )

       ELSE

!
!--       Asynchroneous exchange
          IF ( send_receive == 'ns'  .OR.  send_receive == 'al' )  THEN

             req(1:4)  = 0
             req_count = 0

!
!--          Send front boundary, receive rear one (asynchronous)
             CALL MPI_ISEND( ar(nzb,nys,nxl-nbgp_local),   1, type_xz(grid_level), south_pe,       &
                             req_count, communicator, req(req_count+1), ierr )
             CALL MPI_IRECV( ar(nzb,nyn+1,nxl-nbgp_local), 1, type_xz(grid_level), north_pe,       &
                             req_count, communicator, req(req_count+2), ierr )
!
!--          Send rear boundary, receive front one (asynchronous)
             CALL MPI_ISEND( ar(nzb,nyn-nbgp_local+1,nxl-nbgp_local), 1, type_xz(grid_level),      &
                             north_pe, req_count+1, communicator, req(req_count+3), ierr )
             CALL MPI_IRECV( ar(nzb,nys-nbgp_local,nxl-nbgp_local),   1, type_xz(grid_level),      &
                             south_pe, req_count+1, communicator, req(req_count+4), ierr )

             CALL MPI_WAITALL( 4, req, wait_stat, ierr )

          ENDIF

       ENDIF

    ENDIF

#else

!
!-- Lateral boundary conditions in the non-parallel case.
!-- Case dependent, because in GPU mode still not all arrays are on device. This
!-- workaround has to be removed later. Also, since PGI compiler 12.5 has problems
!-- with array syntax, explicit loops are used.
    IF ( bc_lr_cyc )  THEN
       ar(:,:,nxl-nbgp_local:nxl-1) = ar(:,:,nxr-nbgp_local+1:nxr)
       ar(:,:,nxr+1:nxr+nbgp_local) = ar(:,:,nxl:nxl+nbgp_local-1)
    ENDIF

    !$ACC UPDATE IF_PRESENT ASYNC(1) &
    !$ACC DEVICE(ar(:,:,nxl-nbgp_local:nxl-1)) &
    !$ACC DEVICE(ar(:,:,nxr+1:nxr+nbgp_local))

!
!-- Wait for UPDATES above to complete before starting MPI.
    !$ACC WAIT(2)

    IF ( bc_ns_cyc )  THEN
       ar(:,nys-nbgp_local:nys-1,:) = ar(:,nyn-nbgp_local+1:nyn,:)
       ar(:,nyn+1:nyn+nbgp_local,:) = ar(:,nys:nys+nbgp_local-1,:)
    ENDIF

#endif

#if defined( _OPENACC )
    DO i = nxl-nbgp_local, nxr+nbgp_local
       !$ACC UPDATE IF_PRESENT ASYNC(2) &
       !$ACC DEVICE(ar(:,nys-nbgp_local:nys-1,i)) &
       !$ACC DEVICE(ar(:,nyn+1:nyn+nbgp_local,i))
    ENDDO

!
!-- Wait for all UPDATEs to finish.
    !$ACC WAIT
#endif

    CALL cpu_log( log_point_s(2), 'exchange_horiz', 'stop' )

 END SUBROUTINE exchange_horiz


!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
 SUBROUTINE exchange_horiz_int( ar, nys_l, nyn_l, nxl_l, nxr_l, nzt_l, nbgp_local )


    USE control_parameters,                                                    &
        ONLY:  bc_lr_cyc, bc_ns_cyc

#if defined( __parallel )
    USE control_parameters,                                                    &
        ONLY:  grid_level
#endif
                       
    USE indices,                                                               &
        ONLY:  nzb

    INTEGER(iwp) ::  nxl_l       !< local index bound at current grid level, left side
    INTEGER(iwp) ::  nxr_l       !< local index bound at current grid level, right side
    INTEGER(iwp) ::  nyn_l       !< local index bound at current grid level, north side
    INTEGER(iwp) ::  nys_l       !< local index bound at current grid level, south side
    INTEGER(iwp) ::  nzt_l       !< local index bound at current grid level, top
    INTEGER(iwp) ::  nbgp_local  !< number of ghost points
    
    INTEGER(iwp), DIMENSION(nzb:nzt_l+1,nys_l-nbgp_local:nyn_l+nbgp_local,     &
                            nxl_l-nbgp_local:nxr_l+nbgp_local) ::  ar  !< treated array


#if defined( __parallel )
    IF ( pdims(1) == 1 )  THEN
!
!--    One-dimensional decomposition along y, boundary values can be exchanged
!--    within the PE memory
       IF ( bc_lr_cyc )  THEN
          ar(:,:,nxl_l-nbgp_local:nxl_l-1) = ar(:,:,nxr_l-nbgp_local+1:nxr_l)
          ar(:,:,nxr_l+1:nxr_l+nbgp_local) = ar(:,:,nxl_l:nxl_l+nbgp_local-1)
       ENDIF
    ELSE
!
!--    Send left boundary, receive right one (synchronous)
       CALL MPI_SENDRECV(                                                          &
           ar(nzb,nys_l-nbgp_local,nxl_l),   1, type_yz_int(grid_level), pleft,  0,&
           ar(nzb,nys_l-nbgp_local,nxr_l+1), 1, type_yz_int(grid_level), pright, 0,&
           comm2d, status, ierr )
!
!--    Send right boundary, receive left one (synchronous)
       CALL MPI_SENDRECV(                                                          &
           ar(nzb,nys_l-nbgp_local,nxr_l+1-nbgp_local), 1, type_yz_int(grid_level),&
           pright, 1,                                                              &
           ar(nzb,nys_l-nbgp_local,nxl_l-nbgp_local),   1, type_yz_int(grid_level),&
           pleft,  1,                                                              &
           comm2d, status, ierr )
    ENDIF


    IF ( pdims(2) == 1 )  THEN
!
!--    One-dimensional decomposition along x, boundary values can be exchanged
!--    within the PE memory
       IF ( bc_ns_cyc )  THEN
          ar(:,nys_l-nbgp_local:nys_l-1,:) = ar(:,nyn_l-nbgp_local+1:nyn_l,:)
          ar(:,nyn_l+1:nyn_l+nbgp_local,:) = ar(:,nys_l:nys_l+nbgp_local-1,:)
       ENDIF

    ELSE

!
!--    Send front boundary, receive rear one (synchronous)
       CALL MPI_SENDRECV(                                                          &
           ar(nzb,nys_l,nxl_l-nbgp_local),   1, type_xz_int(grid_level), psouth, 0,&
           ar(nzb,nyn_l+1,nxl_l-nbgp_local), 1, type_xz_int(grid_level), pnorth, 0,&
           comm2d, status, ierr )
!
!--    Send rear boundary, receive front one (synchronous)
       CALL MPI_SENDRECV( ar(nzb,nyn_l-nbgp_local+1,nxl_l-nbgp_local), 1,          &
                          type_xz_int(grid_level), pnorth, 1,                      &
                          ar(nzb,nys_l-nbgp_local,nxl_l-nbgp_local),   1,          &
                          type_xz_int(grid_level), psouth, 1,                      &
                          comm2d, status, ierr )

    ENDIF

#else

    IF ( bc_lr_cyc )  THEN
       ar(:,:,nxl_l-nbgp_local:nxl_l-1) = ar(:,:,nxr_l-nbgp_local+1:nxr_l)
       ar(:,:,nxr_l+1:nxr_l+nbgp_local) = ar(:,:,nxl_l:nxl_l+nbgp_local-1)
    ENDIF

    IF ( bc_ns_cyc )  THEN
       ar(:,nys_l-nbgp_local:nys_l-1,:) = ar(:,nyn_l-nbgp_local+1:nyn_l,:)
       ar(:,nyn_l+1:nyn_l+nbgp_local,:) = ar(:,nys_l:nys_l+nbgp_local-1,:)
    ENDIF

#endif

 END SUBROUTINE exchange_horiz_int

! Description:
! ------------
!> Exchange of lateral (ghost) boundaries (parallel computers) and cyclic
!> boundary conditions, respectively, for 2D-arrays.
!------------------------------------------------------------------------------!
 SUBROUTINE exchange_horiz_2d( ar )

    USE control_parameters,                                                    &
        ONLY :  bc_dirichlet_l, bc_dirichlet_n, bc_dirichlet_r,                &
                bc_dirichlet_s, bc_radiation_l,                                &
                bc_radiation_n, bc_radiation_r, bc_radiation_s

    USE cpulog,                                                                &
        ONLY :  cpu_log, log_point_s

    USE indices,                                                               &
        ONLY :  nbgp, nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg

#if ! defined( __parallel )
    USE control_parameters,                                                    &
        ONLY:  bc_lr_cyc, bc_ns_cyc
#endif


    INTEGER(iwp) :: i  !<

    REAL(wp) ::  ar(nysg:nyng,nxlg:nxrg)  !<


    CALL cpu_log( log_point_s(13), 'exchange_horiz_2d', 'start' )

#if defined( __parallel )

!
!-- Exchange of lateral boundary values for parallel computers
    IF ( pdims(1) == 1 )  THEN

!
!--    One-dimensional decomposition along y, boundary values can be exchanged
!--    within the PE memory
       ar(:,nxlg:nxl-1) = ar(:,nxr-nbgp+1:nxr)
       ar(:,nxr+1:nxrg) = ar(:,nxl:nxl+nbgp-1)

    ELSE
!
!--    Send left boundary, receive right one

       CALL MPI_SENDRECV( ar(nysg,nxl), 1, type_y, pleft,  0,                 &
                          ar(nysg,nxr+1), 1, type_y, pright, 0,               &
                          comm2d, status, ierr )
!
!--    Send right boundary, receive left one
       CALL MPI_SENDRECV( ar(nysg,nxr+1-nbgp), 1, type_y, pright,  1,         &
                          ar(nysg,nxlg), 1, type_y, pleft,   1,               &
                          comm2d, status, ierr )


    ENDIF

    IF ( pdims(2) == 1 )  THEN
!
!--    One-dimensional decomposition along x, boundary values can be exchanged
!--    within the PE memory
       ar(nysg:nys-1,:) = ar(nyn-nbgp+1:nyn,:)
       ar(nyn+1:nyng,:) = ar(nys:nys+nbgp-1,:)

    ELSE
!
!--    Send front boundary, receive rear one

       CALL MPI_SENDRECV( ar(nys,nxlg), 1, type_x, psouth, 0,                 &
                          ar(nyn+1,nxlg), 1, type_x, pnorth, 0,               &
                          comm2d, status, ierr )
!
!--    Send rear boundary, receive front one
       CALL MPI_SENDRECV( ar(nyn+1-nbgp,nxlg), 1, type_x, pnorth, 1,          &
                          ar(nysg,nxlg), 1, type_x, psouth, 1,                &
                          comm2d, status, ierr )

    ENDIF

#else

!
!-- Lateral boundary conditions in the non-parallel case
    IF ( bc_lr_cyc )  THEN
       ar(:,nxlg:nxl-1) = ar(:,nxr-nbgp+1:nxr)
       ar(:,nxr+1:nxrg) = ar(:,nxl:nxl+nbgp-1)
    ENDIF

    IF ( bc_ns_cyc )  THEN
       ar(nysg:nys-1,:) = ar(nyn-nbgp+1:nyn,:)
       ar(nyn+1:nyng,:) = ar(nys:nys+nbgp-1,:)
    ENDIF

#endif

!
!-- Neumann-conditions at inflow/outflow/nested boundaries
    IF ( bc_dirichlet_l  .OR.  bc_radiation_l )  THEN
       DO  i = nbgp, 1, -1
          ar(:,nxl-i) = ar(:,nxl)
       ENDDO
    ENDIF
    IF ( bc_dirichlet_r  .OR.  bc_radiation_r )  THEN
       DO  i = 1, nbgp
          ar(:,nxr+i) = ar(:,nxr)
       ENDDO
    ENDIF
    IF ( bc_dirichlet_s  .OR.  bc_radiation_s )  THEN
       DO  i = nbgp, 1, -1
          ar(nys-i,:) = ar(nys,:)
       ENDDO
    ENDIF
    IF ( bc_dirichlet_n  .OR.  bc_radiation_n )  THEN
       DO  i = 1, nbgp
          ar(nyn+i,:) = ar(nyn,:)
       ENDDO
    ENDIF

    CALL cpu_log( log_point_s(13), 'exchange_horiz_2d', 'stop' )

 END SUBROUTINE exchange_horiz_2d


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Exchange of lateral (ghost) boundaries (parallel computers) and cyclic
!> boundary conditions, respectively, for 2D 8-bit integer arrays.
!------------------------------------------------------------------------------!
 SUBROUTINE exchange_horiz_2d_byte( ar, nys_l, nyn_l, nxl_l, nxr_l, nbgp_local )


    USE control_parameters,                                                    &
        ONLY:  bc_dirichlet_l, bc_dirichlet_n, bc_dirichlet_r, bc_dirichlet_s, &
               bc_radiation_l, bc_radiation_n, bc_radiation_r, bc_radiation_s, &
               bc_radiation_l, bc_radiation_n, bc_radiation_r, bc_radiation_s

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point_s

#if ! defined( __parallel )
    USE control_parameters,                                                    &
        ONLY:  bc_lr_cyc, bc_ns_cyc
#endif

    INTEGER(iwp) ::  i           !< dummy index to zero-gradient conditions at in/outflow boundaries
    INTEGER(iwp) ::  nxl_l       !< local index bound at current grid level, left side
    INTEGER(iwp) ::  nxr_l       !< local index bound at current grid level, right side
    INTEGER(iwp) ::  nyn_l       !< local index bound at current grid level, north side
    INTEGER(iwp) ::  nys_l       !< local index bound at current grid level, south side
    INTEGER(iwp) ::  nbgp_local  !< number of ghost layers to be exchanged

    INTEGER(KIND=1), DIMENSION(nys_l-nbgp_local:nyn_l+nbgp_local,              &
                               nxl_l-nbgp_local:nxr_l+nbgp_local) ::  ar  !< treated array

    CALL cpu_log( log_point_s(13), 'exchange_horiz_2d', 'start' )

#if defined( __parallel )

!
!-- Exchange of lateral boundary values for parallel computers
    IF ( pdims(1) == 1 )  THEN

!
!--    One-dimensional decomposition along y, boundary values can be exchanged
!--    within the PE memory
       ar(:,nxl_l-nbgp_local:nxl_l-1) = ar(:,nxr_l-nbgp_local+1:nxr_l)
       ar(:,nxr_l+1:nxr_l+nbgp_local) = ar(:,nxl_l:nxl_l+nbgp_local-1)

    ELSE
!
!--    Send left boundary, receive right one
       CALL MPI_SENDRECV( ar(nys_l-nbgp_local,nxl_l),   1,                     &
                          type_y_byte, pleft,  0,                              &
                          ar(nys_l-nbgp_local,nxr_l+1), 1,                     &
                          type_y_byte, pright, 0,                              &
                          comm2d, status, ierr )
!
!--    Send right boundary, receive left one
       CALL MPI_SENDRECV( ar(nys_l-nbgp_local,nxr_l+1-nbgp_local), 1,          &
                          type_y_byte, pright, 1,                              &
                          ar(nys_l-nbgp_local,nxl_l-nbgp_local),   1,          &
                          type_y_byte, pleft,  1,                              &
                          comm2d, status, ierr )

    ENDIF

    IF ( pdims(2) == 1 )  THEN
!
!--    One-dimensional decomposition along x, boundary values can be exchanged
!--    within the PE memory
       ar(nys_l-nbgp_local:nys_l-1,:) = ar(nyn_l+1-nbgp_local:nyn_l,:)
       ar(nyn_l+1:nyn_l+nbgp_local,:) = ar(nys_l:nys_l-1+nbgp_local,:)


    ELSE
!
!--    Send front boundary, receive rear one
       CALL MPI_SENDRECV( ar(nys_l,nxl_l-nbgp_local),   1,                    &
                          type_x_byte, psouth, 0,                             &
                          ar(nyn_l+1,nxl_l-nbgp_local), 1,                    &
                          type_x_byte, pnorth, 0,                             &
                          comm2d, status, ierr )

!
!--    Send rear boundary, receive front one
       CALL MPI_SENDRECV( ar(nyn_l+1-nbgp_local,nxl_l-nbgp_local), 1,         &
                          type_x_byte, pnorth, 1,                             &
                          ar(nys_l-nbgp_local,nxl_l-nbgp_local),   1,         &
                          type_x_byte, psouth, 1,                             &
                          comm2d, status, ierr )

    ENDIF

#else

!
!-- Lateral boundary conditions in the non-parallel case
    IF ( bc_lr_cyc )  THEN
       ar(:,nxl_l-nbgp_local:nxl_l-1) = ar(:,nxr_l-nbgp_local+1:nxr_l)
       ar(:,nxr_l+1:nxr_l+nbgp_local) = ar(:,nxl_l:nxl_l+nbgp_local-1)
    ENDIF

    IF ( bc_ns_cyc )  THEN
       ar(nys_l-nbgp_local:nys_l-1,:) = ar(nyn_l+1-nbgp_local:nyn_l,:)
       ar(nyn_l+1:nyn_l+nbgp_local,:) = ar(nys_l:nys_l-1+nbgp_local,:)
    ENDIF

#endif
!
!-- Neumann-conditions at inflow/outflow/nested boundaries
    IF ( bc_dirichlet_l  .OR.  bc_radiation_l )  THEN
       DO  i = nbgp_local, 1, -1
         ar(:,nxl_l-i) = ar(:,nxl_l)
       ENDDO
    ENDIF
    IF ( bc_dirichlet_r  .OR.  bc_radiation_r  )  THEN
       DO  i = 1, nbgp_local
          ar(:,nxr_l+i) = ar(:,nxr_l)
       ENDDO
    ENDIF
    IF ( bc_dirichlet_s  .OR.  bc_radiation_s  )  THEN
       DO  i = nbgp_local, 1, -1
         ar(nys_l-i,:) = ar(nys_l,:)
       ENDDO
    ENDIF
    IF ( bc_dirichlet_n  .OR.  bc_radiation_n  )  THEN
       DO  i = 1, nbgp_local
         ar(nyn_l+i,:) = ar(nyn_l,:)
       ENDDO
    ENDIF

    CALL cpu_log( log_point_s(13), 'exchange_horiz_2d', 'stop' )

 END SUBROUTINE exchange_horiz_2d_byte


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Exchange of lateral (ghost) boundaries (parallel computers) and cyclic
!> boundary conditions, respectively, for 2D 32-bit integer arrays.
!------------------------------------------------------------------------------!
 SUBROUTINE exchange_horiz_2d_int( ar, nys_l, nyn_l, nxl_l, nxr_l, nbgp_local )


    USE control_parameters,                                                    &
        ONLY:  bc_dirichlet_l, bc_dirichlet_n, bc_dirichlet_r, bc_dirichlet_s, &
               bc_radiation_l, bc_radiation_n, bc_radiation_r, bc_radiation_s, &
               bc_radiation_l, bc_radiation_n, bc_radiation_r, bc_radiation_s

#if defined( __parallel )
    USE control_parameters,                                                    &
        ONLY:  grid_level
#endif

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point_s

#if ! defined( __parallel )
    USE control_parameters,                                                    &
        ONLY:  bc_lr_cyc, bc_ns_cyc
#endif

    INTEGER(iwp) ::  i           !< dummy index to zero-gradient conditions at in/outflow boundaries
    INTEGER(iwp) ::  nxl_l       !< local index bound at current grid level, left side
    INTEGER(iwp) ::  nxr_l       !< local index bound at current grid level, right side
    INTEGER(iwp) ::  nyn_l       !< local index bound at current grid level, north side
    INTEGER(iwp) ::  nys_l       !< local index bound at current grid level, south side
    INTEGER(iwp) ::  nbgp_local  !< number of ghost layers to be exchanged

    INTEGER(iwp), DIMENSION(nys_l-nbgp_local:nyn_l+nbgp_local,                 &
                            nxl_l-nbgp_local:nxr_l+nbgp_local) ::  ar  !< treated array

    CALL cpu_log( log_point_s(13), 'exchange_horiz_2d', 'start' )

#if defined( __parallel )

!
!-- Exchange of lateral boundary values for parallel computers
    IF ( pdims(1) == 1 )  THEN

!
!--    One-dimensional decomposition along y, boundary values can be exchanged
!--    within the PE memory
       ar(:,nxl_l-nbgp_local:nxl_l-1) = ar(:,nxr_l-nbgp_local+1:nxr_l)
       ar(:,nxr_l+1:nxr_l+nbgp_local) = ar(:,nxl_l:nxl_l+nbgp_local-1)

    ELSE
!
!--    Send left boundary, receive right one
       CALL MPI_SENDRECV( ar(nys_l-nbgp_local,nxl_l),   1,                     &
                          type_y_int(grid_level), pleft,  0,                   &
                          ar(nys_l-nbgp_local,nxr_l+1), 1,                     &
                          type_y_int(grid_level), pright, 0,                   &
                          comm2d, status, ierr )
!
!--    Send right boundary, receive left one
       CALL MPI_SENDRECV( ar(nys_l-nbgp_local,nxr_l+1-nbgp_local), 1,          &
                          type_y_int(grid_level), pright, 1,                   &
                          ar(nys_l-nbgp_local,nxl_l-nbgp_local),   1,          &
                          type_y_int(grid_level), pleft,  1,                   &
                          comm2d, status, ierr )

    ENDIF

    IF ( pdims(2) == 1 )  THEN
!
!--    One-dimensional decomposition along x, boundary values can be exchanged
!--    within the PE memory
       ar(nys_l-nbgp_local:nys_l-1,:) = ar(nyn_l+1-nbgp_local:nyn_l,:)
       ar(nyn_l+1:nyn_l+nbgp_local,:) = ar(nys_l:nys_l-1+nbgp_local,:)


    ELSE
!
!--    Send front boundary, receive rear one
       CALL MPI_SENDRECV( ar(nys_l,nxl_l-nbgp_local),   1,                    &
                          type_x_int(grid_level), psouth, 0,                  &
                          ar(nyn_l+1,nxl_l-nbgp_local), 1,                    &
                          type_x_int(grid_level), pnorth, 0,                  &
                          comm2d, status, ierr )

!
!--    Send rear boundary, receive front one
       CALL MPI_SENDRECV( ar(nyn_l+1-nbgp_local,nxl_l-nbgp_local), 1,         &
                          type_x_int(grid_level), pnorth, 1,                  &
                          ar(nys_l-nbgp_local,nxl_l-nbgp_local),   1,         &
                          type_x_int(grid_level), psouth, 1,                  &
                          comm2d, status, ierr )

    ENDIF

#else

!
!-- Lateral boundary conditions in the non-parallel case
    IF ( bc_lr_cyc )  THEN
       ar(:,nxl_l-nbgp_local:nxl_l-1) = ar(:,nxr_l-nbgp_local+1:nxr_l)
       ar(:,nxr_l+1:nxr_l+nbgp_local) = ar(:,nxl_l:nxl_l+nbgp_local-1)
    ENDIF

    IF ( bc_ns_cyc )  THEN
       ar(nys_l-nbgp_local:nys_l-1,:) = ar(nyn_l+1-nbgp_local:nyn_l,:)
       ar(nyn_l+1:nyn_l+nbgp_local,:) = ar(nys_l:nys_l-1+nbgp_local,:)
    ENDIF

#endif
!
!-- Neumann-conditions at inflow/outflow/nested boundaries
    IF ( bc_dirichlet_l  .OR.  bc_radiation_l )  THEN
       DO  i = nbgp_local, 1, -1
         ar(:,nxl_l-i) = ar(:,nxl_l)
       ENDDO
    ENDIF
    IF ( bc_dirichlet_r  .OR.  bc_radiation_r  )  THEN
       DO  i = 1, nbgp_local
          ar(:,nxr_l+i) = ar(:,nxr_l)
       ENDDO
    ENDIF
    IF ( bc_dirichlet_s  .OR.  bc_radiation_s  )  THEN
       DO  i = nbgp_local, 1, -1
         ar(nys_l-i,:) = ar(nys_l,:)
       ENDDO
    ENDIF
    IF ( bc_dirichlet_n  .OR.  bc_radiation_n  )  THEN
       DO  i = 1, nbgp_local
         ar(nyn_l+i,:) = ar(nyn_l,:)
       ENDDO
    ENDIF

    CALL cpu_log( log_point_s(13), 'exchange_horiz_2d', 'stop' )

 END SUBROUTINE exchange_horiz_2d_int


 END MODULE exchange_horiz_mod
