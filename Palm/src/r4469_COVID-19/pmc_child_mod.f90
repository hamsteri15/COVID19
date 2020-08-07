MODULE pmc_child

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
! $Id: pmc_child_mod.f90 4360 2020-01-07 11:25:50Z suehring $
! 
! 
! 4182 2019-08-22 15:20:23Z scharf
! Corrected "Former revisions" section
! 
! 3964 2019-05-09 09:48:32Z suehring
! Remove unused variable
! 
! 3963 2019-05-08 20:09:11Z suehring
! Bugfixes in initial settings of child and parent communication patterns. 
!
! 3945 2019-05-02 11:29:27Z raasch
!
! 3932 2019-04-24 17:31:34Z suehring
! typo removed
!
! 2019-02-25 15:31:42Z raasch
! statement added to avoid compiler warning
! 
! 3655 2019-01-07 16:51:22Z knoop
! explicit kind settings
!
! 1762 2016-02-25 12:31:13Z hellstea
! Initial revision by K. Ketelsen
!
! Description:
! ------------
!> Child part of Palm Model Coupler
!------------------------------------------------------------------------------!

#if defined( __parallel )

    USE, INTRINSIC ::  iso_c_binding

#if !defined( __mpifh )
    USE MPI
#endif

    USE kinds

    USE pmc_general,                                                           &
        ONLY:  arraydef, childdef, da_desclen, da_namedef, da_namelen, pedef,  &
               pmc_da_name_err,  pmc_g_setname, pmc_max_array, pmc_status_ok

    USE pmc_handle_communicator,                                               &
        ONLY:  m_model_comm, m_model_npes, m_model_rank, m_to_parent_comm

    USE pmc_mpi_wrapper,                                                       &
        ONLY:  pmc_alloc_mem, pmc_bcast, pmc_inter_bcast, pmc_time

    IMPLICIT NONE

#if defined( __mpifh )
    INCLUDE "mpif.h"
#endif

    PRIVATE
    SAVE

    TYPE(childdef), PUBLIC ::  me   !<

    INTEGER(iwp) ::  myindex = 0         !< counter and unique number for data arrays
    INTEGER(iwp) ::  next_array_in_list = 0   !<


    INTERFACE pmc_childinit
        MODULE PROCEDURE pmc_childinit
    END INTERFACE pmc_childinit

    INTERFACE pmc_c_clear_next_array_list
        MODULE PROCEDURE pmc_c_clear_next_array_list
    END INTERFACE pmc_c_clear_next_array_list

    INTERFACE pmc_c_getbuffer
        MODULE PROCEDURE pmc_c_getbuffer
    END INTERFACE pmc_c_getbuffer

    INTERFACE pmc_c_getnextarray
        MODULE PROCEDURE pmc_c_getnextarray
    END INTERFACE pmc_c_getnextarray

    INTERFACE pmc_c_get_2d_index_list
        MODULE PROCEDURE pmc_c_get_2d_index_list
    END INTERFACE pmc_c_get_2d_index_list

    INTERFACE pmc_c_putbuffer
        MODULE PROCEDURE pmc_c_putbuffer
    END INTERFACE pmc_c_putbuffer

    INTERFACE pmc_c_setind_and_allocmem
        MODULE PROCEDURE pmc_c_setind_and_allocmem
    END INTERFACE pmc_c_setind_and_allocmem

    INTERFACE pmc_c_set_dataarray
        MODULE PROCEDURE pmc_c_set_dataarray_2d
        MODULE PROCEDURE pmc_c_set_dataarray_3d
        MODULE PROCEDURE pmc_c_set_dataarray_ip2d
    END INTERFACE pmc_c_set_dataarray

    INTERFACE pmc_set_dataarray_name
        MODULE PROCEDURE pmc_set_dataarray_name
        MODULE PROCEDURE pmc_set_dataarray_name_lastentry
    END INTERFACE pmc_set_dataarray_name


    PUBLIC pmc_childinit, pmc_c_clear_next_array_list, pmc_c_getbuffer,        &
           pmc_c_getnextarray, pmc_c_putbuffer, pmc_c_setind_and_allocmem,     &
           pmc_c_set_dataarray, pmc_set_dataarray_name, pmc_c_get_2d_index_list

 CONTAINS



 SUBROUTINE pmc_childinit

     IMPLICIT NONE

     INTEGER(iwp) ::  i        !<
     INTEGER(iwp) ::  istat    !<

!
!--  Get / define the MPI environment
     me%model_comm = m_model_comm
     me%inter_comm = m_to_parent_comm

     CALL MPI_COMM_RANK( me%model_comm, me%model_rank, istat )
     CALL MPI_COMM_SIZE( me%model_comm, me%model_npes, istat )
     CALL MPI_COMM_REMOTE_SIZE( me%inter_comm, me%inter_npes, istat )
!
!--  Intra-communicator is used for MPI_GET
     CALL MPI_INTERCOMM_MERGE( me%inter_comm, .TRUE., me%intra_comm, istat )
     CALL MPI_COMM_RANK( me%intra_comm, me%intra_rank, istat )

     ALLOCATE( me%pes(me%inter_npes) )
!
!--  Allocate an array of type arraydef for all parent processes to store 
!--  information of then transfer array
     DO  i = 1, me%inter_npes
        ALLOCATE( me%pes(i)%array_list(pmc_max_array) )
     ENDDO

 END SUBROUTINE pmc_childinit



 SUBROUTINE pmc_set_dataarray_name( parentarraydesc, parentarrayname,          &
                                    childarraydesc, childarrayname, istat )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) ::  parentarrayname  !<
    CHARACTER(LEN=*), INTENT(IN) ::  parentarraydesc  !<
    CHARACTER(LEN=*), INTENT(IN) ::  childarrayname   !<
    CHARACTER(LEN=*), INTENT(IN) ::  childarraydesc   !<

    INTEGER(iwp), INTENT(OUT) ::  istat  !<
!
!-- Local variables
    TYPE(da_namedef) ::  myname  !<

    INTEGER(iwp) ::  mype  !<


    istat = pmc_status_ok
!
!-- Check length of array names
    IF ( LEN( TRIM( parentarrayname) ) > da_namelen  .OR.                      &
         LEN( TRIM( childarrayname) ) > da_namelen )  THEN
       istat = pmc_da_name_err
    ENDIF

    IF ( m_model_rank == 0 )  THEN
       myindex = myindex + 1
       myname%couple_index = myindex
       myname%parentdesc   = TRIM( parentarraydesc )
       myname%nameonparent = TRIM( parentarrayname )
       myname%childdesc    = TRIM( childarraydesc )
       myname%nameonchild  = TRIM( childarrayname )
    ENDIF

!
!-- Broadcast to all child processes
!
!-- The complete description of an transfer names array is broadcasted

    CALL pmc_bcast( myname%couple_index, 0, comm=m_model_comm )
    CALL pmc_bcast( myname%parentdesc,   0, comm=m_model_comm )
    CALL pmc_bcast( myname%nameonparent, 0, comm=m_model_comm )
    CALL pmc_bcast( myname%childdesc,    0, comm=m_model_comm )
    CALL pmc_bcast( myname%nameonchild,  0, comm=m_model_comm )
!
!-- Broadcast to all parent processes
!-- The complete description of an transfer array names is broadcasted als to all parent processe
!   Only the root PE of the broadcasts to parent using intra communicator

    IF ( m_model_rank == 0 )  THEN
        mype = MPI_ROOT
    ELSE
        mype = MPI_PROC_NULL
    ENDIF

    CALL pmc_bcast( myname%couple_index, mype, comm=m_to_parent_comm )
    CALL pmc_bcast( myname%parentdesc,   mype, comm=m_to_parent_comm )
    CALL pmc_bcast( myname%nameonparent, mype, comm=m_to_parent_comm )
    CALL pmc_bcast( myname%childdesc,    mype, comm=m_to_parent_comm )
    CALL pmc_bcast( myname%nameonchild,  mype, comm=m_to_parent_comm )

    CALL pmc_g_setname( me, myname%couple_index, myname%nameonchild )

 END SUBROUTINE pmc_set_dataarray_name



 SUBROUTINE pmc_set_dataarray_name_lastentry( lastentry )

    IMPLICIT NONE

    LOGICAL, INTENT(IN), OPTIONAL ::  lastentry  !<
!
!-- Local variables
    INTEGER ::  idum  !<
    INTEGER ::  mype  !<
    TYPE(dA_namedef) ::  myname  !<

    myname%couple_index = -1

    IF ( m_model_rank == 0 )  THEN
       mype = MPI_ROOT
    ELSE
       mype = MPI_PROC_NULL
    ENDIF

    CALL pmc_bcast( myname%couple_index, mype, comm=m_to_parent_comm )

!
!-- Next statement is just to avoid compiler warnings about unused variables
    IF ( PRESENT( lastentry ) )  idum = 1

 END SUBROUTINE pmc_set_dataarray_name_lastentry



 SUBROUTINE pmc_c_get_2d_index_list

    IMPLICIT NONE

    INTEGER(iwp) :: dummy               !<
    INTEGER(iwp) :: i, ierr, i2, j, nr  !<
    INTEGER(iwp) :: indwin              !< MPI window object
    INTEGER(iwp) :: indwin2             !< MPI window object

    INTEGER(KIND=MPI_ADDRESS_KIND) :: win_size !< Size of MPI window 1 (in bytes)
    INTEGER(KIND=MPI_ADDRESS_KIND) :: disp     !< Displacement unit (Integer = 4, floating poit = 8
    INTEGER(KIND=MPI_ADDRESS_KIND) :: winsize  !< Size of MPI window 2 (in bytes)

    INTEGER, DIMENSION(me%inter_npes*2) :: nrele  !< Number of Elements of a
                                                  !< horizontal slice
    INTEGER, DIMENSION(:), POINTER ::  myind  !<

    TYPE(pedef), POINTER ::  ape  !> Pointer to pedef structure


    win_size = STORAGE_SIZE( dummy )/8
    CALL MPI_WIN_CREATE( dummy, win_size, iwp, MPI_INFO_NULL, me%intra_comm,   &
                         indwin, ierr )
!
!-- Close window on child side and open on parent side
    CALL MPI_WIN_FENCE( 0, indwin, ierr )

!   Between the two MPI_WIN_FENCE calls, the parent can fill the RMA window

!-- Close window on parent side and open on child side

    CALL MPI_WIN_FENCE( 0, indwin, ierr )

    DO  i = 1, me%inter_npes
       disp = me%model_rank * 2
       CALL MPI_GET( nrele((i-1)*2+1), 2, MPI_INTEGER, i-1, disp, 2,           &
                     MPI_INTEGER, indwin, ierr )
    ENDDO
!
!-- MPI_GET is non-blocking -> data in nrele is not available until MPI_FENCE is
!-- called
    CALL MPI_WIN_FENCE( 0, indwin, ierr )
!
!-- Allocate memory for index array
    winsize = 0
    DO  i = 1, me%inter_npes
       ape => me%pes(i)
       i2 = ( i-1 ) * 2 + 1
       nr = nrele(i2+1)
       IF ( nr > 0 )  THEN
          ALLOCATE( ape%locind(nr) )
       ELSE
          NULLIFY( ape%locind )
       ENDIF
       winsize = MAX( INT( nr, MPI_ADDRESS_KIND ), winsize )
    ENDDO

    ALLOCATE( myind(2*winsize) )
    winsize = 1
!
!-- Local buffer used in MPI_GET can but must not be inside the MPI Window.
!-- Here, we use a dummy for the MPI window because the parent processes do 
!-- not access the RMA window via MPI_GET or MPI_PUT
    CALL MPI_WIN_CREATE( dummy, winsize, iwp, MPI_INFO_NULL, me%intra_comm,    &
                         indwin2, ierr )
!
!-- MPI_GET is non-blocking -> data in nrele is not available until MPI_FENCE is
!-- called

    CALL MPI_WIN_FENCE( 0, indwin2, ierr )

!   Between the two MPI_WIN_FENCE calls, the parent can fill the RMA window

    CALL MPI_WIN_FENCE( 0, indwin2, ierr )

    DO  i = 1, me%inter_npes
       ape => me%pes(i)
       nr = nrele(i*2)
       IF ( nr > 0 )  THEN
          disp = nrele(2*(i-1)+1)
          CALL MPI_WIN_LOCK( MPI_LOCK_SHARED , i-1, 0, indwin2, ierr )
          CALL MPI_GET( myind, 2*nr, MPI_INTEGER, i-1, disp, 2*nr,             &
                        MPI_INTEGER, indwin2, ierr )
          CALL MPI_WIN_UNLOCK( i-1, indwin2, ierr )
          DO  j = 1, nr
             ape%locind(j)%i = myind(2*j-1)
             ape%locind(j)%j = myind(2*j)
          ENDDO
          ape%nrele = nr
       ELSE
          ape%nrele = -1
       ENDIF
    ENDDO
!
!-- Don't know why, but this barrier is necessary before we can free the windows
    CALL MPI_BARRIER( me%intra_comm, ierr )

    CALL MPI_WIN_FREE( indWin,  ierr )
    CALL MPI_WIN_FREE( indwin2, ierr )
    DEALLOCATE( myind )

 END SUBROUTINE pmc_c_get_2d_index_list



 SUBROUTINE pmc_c_clear_next_array_list

    IMPLICIT NONE

    next_array_in_list = 0

 END SUBROUTINE pmc_c_clear_next_array_list



 LOGICAL FUNCTION pmc_c_getnextarray( myname )

!
!--  List handling is still required to get minimal interaction with
!--  pmc_interface
     CHARACTER(LEN=*), INTENT(OUT) ::  myname  !<
!
!-- Local variables
    TYPE(pedef), POINTER    :: ape
    TYPE(arraydef), POINTER :: ar


    next_array_in_list = next_array_in_list + 1
!
!-- Array names are the same on all child PEs, so take first process to
!-- get the name    
    ape => me%pes(1)
!
!-- Check if all arrays have been processed
    IF ( next_array_in_list > ape%nr_arrays )  THEN
       pmc_c_getnextarray = .FALSE.
       RETURN
    ENDIF

    ar => ape%array_list( next_array_in_list )

    myname = ar%name
!
!-- Return true if annother array
!-- If all array have been processed, the RETURN statement a couple of lines above is active

    pmc_c_getnextarray = .TRUE.

 END FUNCTION pmc_c_getnextarray



 SUBROUTINE pmc_c_set_dataarray_2d( array )

    IMPLICIT NONE

    REAL(wp), INTENT(IN) , DIMENSION(:,:), POINTER ::  array  !<

    INTEGER(iwp)               ::  i       !<
    INTEGER(iwp)               ::  nrdims  !<
    INTEGER(iwp), DIMENSION(4) ::  dims    !<

    TYPE(C_PTR)             ::  array_adr
    TYPE(arraydef), POINTER ::  ar
    TYPE(pedef), POINTER    ::  ape


    dims    = 1
    nrdims  = 2
    dims(1) = SIZE( array, 1 )
    dims(2) = SIZE( array, 2 )

    array_adr = C_LOC( array )

    DO  i = 1, me%inter_npes
       ape => me%pes(i)
       ar  => ape%array_list(next_array_in_list)
       ar%nrdims = nrdims
       ar%dimkey = nrdims
       ar%a_dim  = dims
       ar%data   = array_adr
    ENDDO

 END SUBROUTINE pmc_c_set_dataarray_2d

 SUBROUTINE pmc_c_set_dataarray_ip2d( array )

    IMPLICIT NONE

    INTEGER(idp), INTENT(IN) , DIMENSION(:,:), POINTER ::  array  !<

    INTEGER(iwp)               ::  i       !<
    INTEGER(iwp)               ::  nrdims  !<
    INTEGER(iwp), DIMENSION(4) ::  dims    !<

    TYPE(C_PTR)             ::  array_adr
    TYPE(arraydef), POINTER ::  ar
    TYPE(pedef), POINTER    ::  ape

    dims    = 1
    nrdims  = 2
    dims(1) = SIZE( array, 1 )
    dims(2) = SIZE( array, 2 )

    array_adr = C_LOC( array )

    DO  i = 1, me%inter_npes
       ape => me%pes(i)
       ar  => ape%array_list(next_array_in_list)
       ar%nrdims = nrdims
       ar%dimkey = 22
       ar%a_dim  = dims
       ar%data   = array_adr
    ENDDO

 END SUBROUTINE pmc_c_set_dataarray_ip2d

 SUBROUTINE pmc_c_set_dataarray_3d (array)

    IMPLICIT NONE

    REAL(wp), INTENT(IN), DIMENSION(:,:,:), POINTER ::  array  !<

    INTEGER(iwp)                ::  i
    INTEGER(iwp)                ::  nrdims
    INTEGER(iwp), DIMENSION (4) ::  dims
    
    TYPE(C_PTR)             ::  array_adr
    TYPE(pedef), POINTER    ::  ape
    TYPE(arraydef), POINTER ::  ar


    dims    = 1
    nrdims  = 3
    dims(1) = SIZE( array, 1 )
    dims(2) = SIZE( array, 2 )
    dims(3) = SIZE( array, 3 )

    array_adr = C_LOC( array )

    DO  i = 1, me%inter_npes
       ape => me%pes(i)
       ar  => ape%array_list(next_array_in_list)
       ar%nrdims = nrdims
       ar%dimkey = nrdims
       ar%a_dim  = dims
       ar%data   = array_adr
    ENDDO

 END SUBROUTINE pmc_c_set_dataarray_3d



 SUBROUTINE pmc_c_setind_and_allocmem

    IMPLICIT NONE

!
!-- Naming convention for appendices:  _pc  -> parent to child transfer
!--                                    _cp  -> child to parent transfer
!--                                    recv -> parent to child transfer
!--                                    send -> child to parent transfer
    INTEGER(iwp) ::  arlen        !<
    INTEGER(iwp) ::  myindex      !<
    INTEGER(iwp) ::  i            !<
    INTEGER(iwp) ::  ierr         !<
    INTEGER(iwp) ::  istat        !<
    INTEGER(iwp) ::  j            !<
    INTEGER(iwp) ::  lo_nr_arrays !<
    INTEGER(iwp) ::  rcount       !<
    INTEGER(iwp) ::  tag          !<
    INTEGER(iwp) ::  total_npes   !<

    INTEGER(iwp), PARAMETER ::  noindex = -1  !<

    INTEGER(idp)                   ::  bufsize  !< size of MPI data window
    INTEGER(KIND=MPI_ADDRESS_KIND) ::  winsize  !<
    
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  myindex_s
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  myindex_r

    REAL(wp), DIMENSION(:), POINTER, SAVE ::  base_array_pc  !< base array
    REAL(wp), DIMENSION(:), POINTER, SAVE ::  base_array_cp  !< base array

    TYPE(pedef), POINTER    ::  ape       !<
    TYPE(arraydef), POINTER ::  ar        !<
    Type(C_PTR)             ::  base_ptr  !<

  
    CALL MPI_COMM_SIZE (me%intra_comm, total_npes, ierr)

    lo_nr_arrays = me%pes(1)%nr_arrays

    ALLOCATE(myindex_s(lo_nr_arrays,0:total_npes-1))
    ALLOCATE(myindex_r(lo_nr_arrays,0:total_npes-1))

    myindex_s = 0

!
!-- Receive indices from child
    CALL MPI_ALLTOALL( myindex_s, lo_nr_arrays, MPI_INTEGER,                   &
                       myindex_r, lo_nr_arrays, MPI_INTEGER,                   &
                       me%intra_comm, ierr )

    myindex = 0
    bufsize = 8
!
!-- Parent to child direction.
!-- First stride: compute size and set index
    DO  i = 1, me%inter_npes
       ape => me%pes(i)
       DO  j = 1, ape%nr_arrays
          ar => ape%array_list(j)
          ar%recvindex = myindex_r(j,i-1)
!
!--       Determine max, because child buffer is allocated only once
!--       All 2D and 3d arrays use the same buffer

          IF ( ar%nrdims == 3 )  THEN
             bufsize = MAX( bufsize,                                           &
                            INT( ar%a_dim(1)*ar%a_dim(2)*ar%a_dim(3),          &
                                 MPI_ADDRESS_KIND ) )
          ELSE
             bufsize = MAX( bufsize,                                           &
                            INT( ar%a_dim(1)*ar%a_dim(2), MPI_ADDRESS_KIND ) )
          ENDIF
       ENDDO
    ENDDO

!
!-- Create RMA (one sided communication) data buffer.
!-- The buffer for MPI_GET can be PE local, i.e. it can but must not be part of
!-- the MPI RMA window
    CALL pmc_alloc_mem( base_array_pc, bufsize, base_ptr )
    me%totalbuffersize = bufsize*wp  ! total buffer size in byte
!
!-- Second stride: set buffer pointer
    DO  i = 1, me%inter_npes
       ape => me%pes(i)
       DO  j = 1, ape%nr_arrays
          ar => ape%array_list(j)
          ar%recvbuf = base_ptr
       ENDDO
    ENDDO
!
!-- Child to parent direction
    myindex = 1
    rcount  = 0
    bufsize = 8

    myindex_s = 0
    myindex_r = 0

    DO  i = 1, me%inter_npes
       ape => me%pes(i)
       tag = 300
       DO  j = 1, ape%nr_arrays
          ar => ape%array_list(j)
          IF ( ar%nrdims == 2 )  THEN
             arlen = ape%nrele
          ELSEIF( ar%nrdims == 3 )  THEN
             arlen = ape%nrele*ar%a_dim(1)
          ENDIF

          IF ( ape%nrele > 0 )  THEN
             ar%sendindex = myindex
          ELSE
             ar%sendindex = noindex
          ENDIF

          myindex_s(j,i-1) = ar%sendindex

          IF ( ape%nrele > 0 )  THEN
             ar%sendsize = arlen
             myindex     = myindex + arlen
             bufsize     = bufsize + arlen
          ENDIF

       ENDDO

    ENDDO
!
!-- Send indices to parent

    CALL MPI_ALLTOALL( myindex_s, lo_nr_arrays, MPI_INTEGER,                   &
                       myindex_r, lo_nr_arrays, MPI_INTEGER,                   &
                       me%intra_comm, ierr)

    DEALLOCATE( myindex_s )
    DEALLOCATE( myindex_r )

!
!-- Create RMA (one sided communication) window for data buffer child to parent
!-- transfer.
!-- The buffer of MPI_GET (counter part of transfer) can be PE-local, i.e. it
!-- can but must not be part of the MPI RMA window. Only one RMA window is
!-- required to prepare the data
!--        for parent -> child transfer on the parent side
!-- and
!--        for child -> parent transfer on the child side
    CALL pmc_alloc_mem( base_array_cp, bufsize )
    me%totalbuffersize = bufsize * wp  ! total buffer size in byte

    winSize = me%totalbuffersize

    CALL MPI_WIN_CREATE( base_array_cp, winsize, wp, MPI_INFO_NULL,            &
                         me%intra_comm, me%win_parent_child, ierr )
    CALL MPI_WIN_FENCE( 0, me%win_parent_child, ierr )
    CALL MPI_BARRIER( me%intra_comm, ierr )
!
!-- Second stride: set buffer pointer
    DO  i = 1, me%inter_npes
       ape => me%pes(i)
       DO  j = 1, ape%nr_arrays
          ar => ape%array_list(j)          
          IF ( ape%nrele > 0 )  THEN
             ar%sendbuf = C_LOC( base_array_cp(ar%sendindex) )
!
!--          TODO: if this is an error to be really expected, replace the
!--                following message by a meaningful standard PALM message using
!--                the message-routine
             IF ( ar%sendindex+ar%sendsize > bufsize )  THEN
                WRITE( 0,'(a,i4,4i7,1x,a)') 'Child buffer too small ', i,      &
                          ar%sendindex, ar%sendsize, ar%sendindex+ar%sendsize, &
                          bufsize, TRIM( ar%name )
                CALL MPI_ABORT( MPI_COMM_WORLD, istat, ierr )
             ENDIF
          ENDIF
       ENDDO
    ENDDO

 END SUBROUTINE pmc_c_setind_and_allocmem



 SUBROUTINE pmc_c_getbuffer( waittime, particle_transfer )

    IMPLICIT NONE

    REAL(wp), INTENT(OUT), OPTIONAL ::  waittime  !<
    LOGICAL, INTENT(IN), OPTIONAL   ::  particle_transfer  !<

    LOGICAL                        ::  lo_ptrans!<
    
    INTEGER(iwp)                        ::  ierr    !<
    INTEGER(iwp)                        ::  ij      !<
    INTEGER(iwp)                        ::  ip      !<
    INTEGER(iwp)                        ::  j       !<
    INTEGER(iwp)                        ::  myindex !<
    INTEGER(iwp)                        ::  nr      !< number of elements to get
                                                    !< from parent
    INTEGER(KIND=MPI_ADDRESS_KIND) ::  target_disp
    INTEGER,DIMENSION(1)           ::  buf_shape

    REAL(wp)                            ::  t1
    REAL(wp)                            ::  t2

    REAL(wp), POINTER, DIMENSION(:)     ::  buf
    REAL(wp), POINTER, DIMENSION(:,:)   ::  data_2d
    REAL(wp), POINTER, DIMENSION(:,:,:) ::  data_3d
    TYPE(pedef), POINTER                ::  ape
    TYPE(arraydef), POINTER             ::  ar
    INTEGER(idp), POINTER, DIMENSION(:)     ::  ibuf      !<
    INTEGER(idp), POINTER, DIMENSION(:,:)   ::  idata_2d  !<

!
!-- Synchronization of the model is done in pmci_synchronize. 
!-- Therefore the RMA window can be filled without
!-- sychronization at this point and a barrier is not necessary.

!-- In case waittime is present, the following barrier is necessary to
!-- insure the same number of barrier calls on parent and child
!-- This means, that here on child side two barriers are call successively
!-- The parent is filling its buffer between the two barrier calls

!-- Please note that waittime has to be set in pmc_s_fillbuffer AND
!-- pmc_c_getbuffer
    IF ( PRESENT( waittime ) )  THEN
       t1 = pmc_time()
       CALL MPI_BARRIER( me%intra_comm, ierr )
       t2 = pmc_time()
       waittime = t2 - t1
    ENDIF

    lo_ptrans = .FALSE.
    IF ( PRESENT( particle_transfer))    lo_ptrans = particle_transfer

!
!-- Wait for buffer is filled.
!
!-- The parent side (in pmc_s_fillbuffer) is filling the buffer in the MPI RMA window
!-- When the filling is complet, a MPI_BARRIER is called.
!-- The child is not allowd to access the parent-buffer before it is completely filled
!-- therefore the following barrier is required.

    CALL MPI_BARRIER( me%intra_comm, ierr )

    DO  ip = 1, me%inter_npes
       ape => me%pes(ip)
       DO  j = 1, ape%nr_arrays
          ar => ape%array_list(j)

          IF ( ar%dimkey == 2 .AND. .NOT.lo_ptrans)  THEN
             nr = ape%nrele
          ELSEIF ( ar%dimkey == 3 .AND. .NOT.lo_ptrans)  THEN
             nr = ape%nrele * ar%a_dim(1)
          ELSE IF ( ar%dimkey == 22 .AND. lo_ptrans)  THEN
             nr = ape%nrele
          ELSE
             CYCLE                    ! Particle array ar not transferd here
          ENDIF
          buf_shape(1) = nr
          IF ( lo_ptrans )   THEN
             CALL C_F_POINTER( ar%recvbuf, ibuf, buf_shape )
          ELSE
             CALL C_F_POINTER( ar%recvbuf, buf, buf_shape )
          ENDIF
!
!--       MPI passive target RMA
!--       One data array is fetcht from MPI RMA window on parent

          IF ( nr > 0 )  THEN
             target_disp = ar%recvindex - 1
             CALL MPI_WIN_LOCK( MPI_LOCK_SHARED , ip-1, 0,                     &
                                me%win_parent_child, ierr )
             IF ( lo_ptrans )   THEN
                CALL MPI_GET( ibuf, nr*8, MPI_BYTE, ip-1, target_disp, nr*8, MPI_BYTE,  &               !There is no MPI_INTEGER8 datatype
                                   me%win_parent_child, ierr )
             ELSE
                CALL MPI_GET( buf, nr, MPI_REAL, ip-1, target_disp, nr,        &
                              MPI_REAL, me%win_parent_child, ierr )
             ENDIF
             CALL MPI_WIN_UNLOCK( ip-1, me%win_parent_child, ierr )
          ENDIF
          myindex = 1
          IF ( ar%dimkey == 2 .AND. .NOT.lo_ptrans)  THEN

             CALL C_F_POINTER( ar%data, data_2d, ar%a_dim(1:2) )
             DO  ij = 1, ape%nrele
                data_2d(ape%locind(ij)%j,ape%locind(ij)%i) = buf(myindex)
                myindex = myindex + 1
             ENDDO

          ELSEIF ( ar%dimkey == 3 .AND. .NOT.lo_ptrans)  THEN

             CALL C_F_POINTER( ar%data, data_3d, ar%a_dim(1:3) )
             DO  ij = 1, ape%nrele
                data_3d(:,ape%locind(ij)%j,ape%locind(ij)%i) =                 &
                                              buf(myindex:myindex+ar%a_dim(1)-1)
                myindex = myindex+ar%a_dim(1)
             ENDDO

          ELSEIF ( ar%dimkey == 22 .AND. lo_ptrans)  THEN
             CALL C_F_POINTER( ar%data, idata_2d, ar%a_dim(1:2) )

             DO  ij = 1, ape%nrele
                idata_2d(ape%locind(ij)%j,ape%locind(ij)%i) = ibuf(myindex)
                myindex = myindex + 1
             ENDDO

          ENDIF
       ENDDO
    ENDDO

 END SUBROUTINE pmc_c_getbuffer



 SUBROUTINE pmc_c_putbuffer( waittime , particle_transfer )

    IMPLICIT NONE

    REAL(wp), INTENT(OUT), OPTIONAL ::  waittime  !<
    LOGICAL, INTENT(IN), OPTIONAL   ::  particle_transfer  !<

    LOGICAL ::  lo_ptrans!<
    
    INTEGER(iwp) ::  ierr         !<
    INTEGER(iwp) ::  ij           !<
    INTEGER(iwp) ::  ip           !<
    INTEGER(iwp) ::  j            !<
    INTEGER(iwp) ::  myindex      !<

    INTEGER(iwp), DIMENSION(1) ::  buf_shape    !<

    REAL(wp) ::  t1  !<
    REAL(wp) ::  t2  !<

    REAL(wp), POINTER, DIMENSION(:)         ::  buf      !<
    REAL(wp), POINTER, DIMENSION(:,:)       ::  data_2d  !<
    REAL(wp), POINTER, DIMENSION(:,:,:)     ::  data_3d  !<
    
    INTEGER(idp), POINTER, DIMENSION(:)     ::  ibuf      !<
    INTEGER(idp), POINTER, DIMENSION(:,:)   ::  idata_2d  !<

    TYPE(pedef), POINTER                    ::  ape  !<
    TYPE(arraydef), POINTER                 ::  ar   !<

!
!-- Wait for empty buffer
!-- Switch RMA epoche

    t1 = pmc_time()
    CALL MPI_BARRIER( me%intra_comm, ierr )
    t2 = pmc_time()
    IF ( PRESENT( waittime ) )  waittime = t2 - t1

    lo_ptrans = .FALSE.
    IF ( PRESENT( particle_transfer))    lo_ptrans = particle_transfer

    DO  ip = 1, me%inter_npes
       ape => me%pes(ip)
       DO  j = 1, ape%nr_arrays
          ar => aPE%array_list(j)
          myindex = 1

          IF ( ar%dimkey == 2 .AND. .NOT.lo_ptrans )  THEN

             buf_shape(1) = ape%nrele
             CALL C_F_POINTER( ar%sendbuf, buf,     buf_shape     )
             CALL C_F_POINTER( ar%data,    data_2d, ar%a_dim(1:2) )
             DO  ij = 1, ape%nrele
                buf(myindex) = data_2d(ape%locind(ij)%j,ape%locind(ij)%i)
                myindex = myindex + 1
             ENDDO

          ELSEIF ( ar%dimkey == 3 .AND. .NOT.lo_ptrans )  THEN

             buf_shape(1) = ape%nrele*ar%a_dim(1)
             CALL C_F_POINTER( ar%sendbuf, buf,     buf_shape     )
             CALL C_F_POINTER( ar%data,    data_3d, ar%a_dim(1:3) )
             DO  ij = 1, ape%nrele
                buf(myindex:myindex+ar%a_dim(1)-1) =                            &
                                    data_3d(:,ape%locind(ij)%j,ape%locind(ij)%i)
                myindex = myindex + ar%a_dim(1)
             ENDDO

          ELSE IF ( ar%dimkey == 22 .AND. lo_ptrans)  THEN

             buf_shape(1) = ape%nrele
             CALL C_F_POINTER( ar%sendbuf, ibuf,     buf_shape     )
             CALL C_F_POINTER( ar%data,    idata_2d, ar%a_dim(1:2) )

             DO  ij = 1, ape%nrele
                ibuf(myindex) = idata_2d(ape%locind(ij)%j,ape%locind(ij)%i)
                myindex = myindex + 1
             ENDDO

          ENDIF
       ENDDO
    ENDDO
!
!-- Buffer is filled
!-- Switch RMA epoche

    CALL MPI_Barrier(me%intra_comm, ierr)

 END SUBROUTINE pmc_c_putbuffer

#endif
 END MODULE pmc_child
