!> @file close_file.f90
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
! $Id: close_file.f90 4360 2020-01-07 11:25:50Z suehring $
! Corrected "Former revisions" section
! 
! 4069 2019-07-01 14:05:51Z Giersch
! Masked output running index mid has been introduced as a local variable to 
! avoid runtime error (Loop variable has been modified) in time_integration
! 
! 3655 2019-01-07 16:51:22Z knoop
! unused variables and format statements removed
!
! Revision 1.1 (close_files) 1997/08/11 06:11:18  raasch
! Initial revision
!
!
! Description:
! ------------
!> Close specified file or all open files, if "0" has been given as the
!> calling argument. In that case, execute last actions for certain unit 
!> numbers, if required.
!------------------------------------------------------------------------------!
 SUBROUTINE close_file( file_id )
 

    USE control_parameters,                                                    &
        ONLY:  max_masks, openfile
                
    USE kinds
    
#if defined( __netcdf )
    USE NETCDF
#endif

    USE netcdf_interface,                                                      &
        ONLY:  id_set_mask, id_set_pr, id_set_pts, id_set_sp,                  &
               id_set_ts, id_set_xy, id_set_xz, id_set_yz, id_set_3d,          &
               id_set_fl, nc_stat, netcdf_data_format, netcdf_handle_error
                
    USE pegrid                                           

    IMPLICIT NONE

    CHARACTER (LEN=10)  ::  datform = 'lit_endian' !< 
    CHARACTER (LEN=80)  ::  title                  !< 

    INTEGER(iwp) ::  av           !< 
    INTEGER(iwp) ::  dimx         !< 
    INTEGER(iwp) ::  dimy         !< 
    INTEGER(iwp) ::  fid          !< 
    INTEGER(iwp) ::  file_id      !< 
    INTEGER(iwp) ::  mid          !< masked output running index 
    INTEGER(iwp) ::  planz        !< 

    LOGICAL ::  checkuf = .TRUE.  !< 
    LOGICAL ::  datleg = .TRUE.   !< 
    LOGICAL ::  dbp = .FALSE.     !< 

    NAMELIST /GLOBAL/  checkuf, datform, dimx, dimy, dbp, planz,               &
                       title
    NAMELIST /RAHMEN/  datleg

!
!-- Close specified unit number (if opened) and set a flag that it has
!-- been opened one time at least
    IF ( file_id /= 0 )  THEN
       IF ( openfile(file_id)%opened )  THEN
          CLOSE ( file_id )
          openfile(file_id)%opened        = .FALSE.
          openfile(file_id)%opened_before = .TRUE.
       ENDIF
       RETURN
    ENDIF

!
!-- Close all open unit numbers
    DO  fid = 1, 200+2*max_masks

       IF ( openfile(fid)%opened .OR. openfile(fid)%opened_before )  THEN
!
!--       Last actions for certain unit numbers
          SELECT CASE ( fid )

#if defined( __netcdf )
             CASE ( 101 )

                IF ( myid == 0  .OR.  netcdf_data_format > 4 )  THEN
                   nc_stat = NF90_CLOSE( id_set_xy(0) )
                   CALL netcdf_handle_error( 'close_file', 44 )
                ENDIF

             CASE ( 102 )

                IF ( myid == 0  .OR.  netcdf_data_format > 4 )  THEN
                   nc_stat = NF90_CLOSE( id_set_xz(0) )
                   CALL netcdf_handle_error( 'close_file', 45 )
                ENDIF

             CASE ( 103 )

                IF ( myid == 0  .OR.  netcdf_data_format > 4 )  THEN
                   nc_stat = NF90_CLOSE( id_set_yz(0) )
                   CALL netcdf_handle_error( 'close_file', 46 )
                ENDIF

             CASE ( 104 )

                IF ( myid == 0 )  THEN
                   nc_stat = NF90_CLOSE( id_set_pr )
                   CALL netcdf_handle_error( 'close_file', 47 )
                ENDIF

             CASE ( 105 )

                IF ( myid == 0 )  THEN
                   nc_stat = NF90_CLOSE( id_set_ts )
                   CALL netcdf_handle_error( 'close_file', 48 )
                ENDIF

             CASE ( 106 )

                IF ( myid == 0  .OR.  netcdf_data_format > 4 )  THEN
                   nc_stat = NF90_CLOSE( id_set_3d(0) )
                   CALL netcdf_handle_error( 'close_file', 49 )
                ENDIF

             CASE ( 107 )

                IF ( myid == 0 )  THEN
                   nc_stat = NF90_CLOSE( id_set_sp )
                   CALL netcdf_handle_error( 'close_file', 50 )
                ENDIF

!
!--           Currently disabled
!             CASE ( 108 )

!                nc_stat = NF90_CLOSE( id_set_prt )
!                CALL netcdf_handle_error( 'close_file', 51 )

             CASE ( 109 ) 

                nc_stat = NF90_CLOSE( id_set_pts )
                CALL netcdf_handle_error( 'close_file', 412 )

             CASE ( 111 )

                IF ( myid == 0  .OR.  netcdf_data_format > 4 )  THEN
                   nc_stat = NF90_CLOSE( id_set_xy(1) )
                   CALL netcdf_handle_error( 'close_file', 52 )
                ENDIF

             CASE ( 112 )

                IF ( myid == 0  .OR.  netcdf_data_format > 4 )  THEN
                   nc_stat = NF90_CLOSE( id_set_xz(1) )
                   CALL netcdf_handle_error( 'close_file', 352 )
                ENDIF

             CASE ( 113 )

                IF ( myid == 0  .OR.  netcdf_data_format > 4 )  THEN
                   nc_stat = NF90_CLOSE( id_set_yz(1) )
                   CALL netcdf_handle_error( 'close_file', 353 )
                ENDIF

             CASE ( 116 )

                IF ( myid == 0  .OR.  netcdf_data_format > 4 )  THEN
                   nc_stat = NF90_CLOSE( id_set_3d(1) )
                   CALL netcdf_handle_error( 'close_file', 353 )
                ENDIF

             CASE ( 199 )

                IF ( myid == 0  .OR.  netcdf_data_format > 4 )  THEN
                   nc_stat = NF90_CLOSE( id_set_fl )
                   CALL netcdf_handle_error( 'close_file', 353 )
                ENDIF

             CASE ( 201:200+2*max_masks )

                IF ( myid == 0  .OR.  netcdf_data_format > 4 )  THEN
!
!--                decompose fid into mid and av
                   IF ( fid <= 200+max_masks )  THEN
                      mid = fid - 200
                      av = 0
                   ELSE
                      mid = fid - (200+max_masks)
                      av = 1
                   ENDIF
                   nc_stat = NF90_CLOSE( id_set_mask(mid,av) )
                   CALL netcdf_handle_error( 'close_file', 459 )

                ENDIF

#endif

          END SELECT
!
!--       Close file
          IF ( openfile(fid)%opened )  CLOSE ( fid )

       ENDIF

    ENDDO

 END SUBROUTINE close_file
