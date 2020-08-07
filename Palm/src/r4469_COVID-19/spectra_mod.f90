!> @file spectra_mod.f90
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
! $Id: spectra_mod.f90 4429 2020-02-27 15:24:30Z raasch $
! bugfix: preprocessor directives rearranged for serial mode
! 
! 4360 2020-01-07 11:25:50Z suehring
! Corrected "Former revisions" section
! 
! 3805 2019-03-20 15:26:35Z raasch
! unused variable removed
! 
! 3655 2019-01-07 16:51:22Z knoop
! Renamed output variables
!
! Revision 1.1  2001/01/05 15:08:07  raasch
! Initial revision
!
!
! Description:
! ------------
!> Calculate horizontal spectra along x and y. 
!> ATTENTION: 1d-decomposition along y still needs improvement, because in that
!>            case the gridpoint number along z still depends on the PE number
!>            because transpose_xz has to be used (and possibly also
!>            transpose_zyd needs modification).
!------------------------------------------------------------------------------!
 MODULE spectra_mod

    USE kinds

    PRIVATE

    CHARACTER (LEN=2),  DIMENSION(10) ::  spectra_direction = 'x'
    CHARACTER (LEN=10), DIMENSION(10) ::  data_output_sp  = ' '

    INTEGER(iwp) ::  average_count_sp = 0
    INTEGER(iwp) ::  dosp_time_count = 0
    INTEGER(iwp) ::  n_sp_x = 0, n_sp_y = 0

    INTEGER(iwp) ::  comp_spectra_level(100) = 999999

    LOGICAL ::  calculate_spectra   = .FALSE.  !< internal switch that spectra are calculated
    LOGICAL ::  spectra_initialized = .FALSE.  !< internal switch that spectra related quantities are initialized

    REAL(wp) ::  averaging_interval_sp = 9999999.9_wp  !< averaging interval for spectra output
    REAL(wp) ::  dt_dosp = 9999999.9_wp                !< time interval for spectra output
    REAL(wp) ::  skip_time_dosp = 9999999.9_wp         !< no output of spectra data before this interval has passed

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  var_d

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  spectrum_x, spectrum_y

    SAVE

    INTERFACE calc_spectra
       MODULE PROCEDURE calc_spectra
    END INTERFACE calc_spectra

    INTERFACE preprocess_spectra
       MODULE PROCEDURE preprocess_spectra
    END INTERFACE preprocess_spectra

#if defined( __parallel )
    INTERFACE calc_spectra_x
       MODULE PROCEDURE calc_spectra_x
    END INTERFACE calc_spectra_x

    INTERFACE calc_spectra_y
       MODULE PROCEDURE calc_spectra_y
    END INTERFACE calc_spectra_y
#endif

    INTERFACE spectra_check_parameters
       MODULE PROCEDURE spectra_check_parameters
    END INTERFACE spectra_check_parameters

    INTERFACE spectra_header
       MODULE PROCEDURE spectra_header
    END INTERFACE spectra_header

    INTERFACE spectra_init
       MODULE PROCEDURE spectra_init
    END INTERFACE spectra_init

    INTERFACE spectra_parin
       MODULE PROCEDURE spectra_parin
    END INTERFACE spectra_parin

    PUBLIC average_count_sp, averaging_interval_sp, calc_spectra,              &
           calculate_spectra, comp_spectra_level, data_output_sp,              &
           dosp_time_count, dt_dosp, n_sp_x, n_sp_y,                           &
           skip_time_dosp, spectra_check_parameters, spectra_direction,        &
           spectra_header, spectra_init, spectra_parin, spectrum_x,            &
           spectrum_y, var_d


 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Parin for &spectra_par for calculating spectra
!------------------------------------------------------------------------------!
    SUBROUTINE spectra_parin

       USE control_parameters,                                                 &
           ONLY:  dt_data_output, message_string

       IMPLICIT NONE

       CHARACTER (LEN=80) ::  line  !< dummy string that contains the current  &
                                    !< line of the parameter file

       NAMELIST /spectra_par/  averaging_interval_sp, comp_spectra_level,      &
                               data_output_sp, dt_dosp, skip_time_dosp,        &
                               spectra_direction

       NAMELIST /spectra_parameters/                                           &
                               averaging_interval_sp, comp_spectra_level,      &
                               data_output_sp, dt_dosp, skip_time_dosp,        &
                               spectra_direction
!
!--    Position the namelist-file at the beginning (it was already opened in
!--    parin), search for the namelist-group of the package and position the
!--    file at this line.
       line = ' '

!
!--    Try to find the spectra package
       REWIND ( 11 )
       line = ' '
       DO WHILE ( INDEX( line, '&spectra_parameters' ) == 0 )
          READ ( 11, '(A)', END=12 )  line
       ENDDO
       BACKSPACE ( 11 )

!
!--    Read namelist
       READ ( 11, spectra_parameters, ERR = 10 )

!
!--    Default setting of dt_dosp here (instead of check_parameters), because
!--    its current value is needed in init_pegrid
       IF ( dt_dosp == 9999999.9_wp )  dt_dosp = dt_data_output

!
!--    Set general switch that spectra shall be calculated
       calculate_spectra = .TRUE.

       GOTO 14

 10    BACKSPACE( 11 )
       READ( 11 , '(A)') line
       CALL parin_fail_message( 'spectra_parameters', line )
!
!--    Try to find the old namelist
 12    REWIND ( 11 )
       line = ' '
       DO WHILE ( INDEX( line, '&spectra_par' ) == 0 )
          READ ( 11, '(A)', END=14 )  line
       ENDDO
       BACKSPACE ( 11 )

!
!--    Read namelist
       READ ( 11, spectra_par, ERR = 13, END = 14 )

       
       message_string = 'namelist spectra_par is deprecated and will be ' // &
                     'removed in near future. Please use namelist ' //       &
                     'spectra_parameters instead'
       CALL message( 'spectra_parin', 'PA0487', 0, 1, 0, 6, 0 )
!
!--    Default setting of dt_dosp here (instead of check_parameters), because
!--    its current value is needed in init_pegrid
       IF ( dt_dosp == 9999999.9_wp )  dt_dosp = dt_data_output

!
!--    Set general switch that spectra shall be calculated
       calculate_spectra = .TRUE.
       
       GOTO 14

 13    BACKSPACE( 11 )
       READ( 11 , '(A)') line
       CALL parin_fail_message( 'spectra_par', line )
       
       
 14    CONTINUE

    END SUBROUTINE spectra_parin



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of spectra related variables
!------------------------------------------------------------------------------!
    SUBROUTINE spectra_init

       USE indices,                                                            &
           ONLY:  nx, ny, nzb, nzt

       IMPLICIT NONE

       IF ( spectra_initialized )  RETURN

       IF ( dt_dosp /= 9999999.9_wp )  THEN

          IF ( .NOT. ALLOCATED( spectrum_x ) )  THEN
             ALLOCATE( spectrum_x( 1:nx/2, 1:100, 1:10 ) )
             spectrum_x = 0.0_wp
          ENDIF

          IF ( .NOT. ALLOCATED( spectrum_y ) )  THEN
             ALLOCATE( spectrum_y( 1:ny/2, 1:100, 1:10 ) )
             spectrum_y = 0.0_wp
          ENDIF

          ALLOCATE( var_d(nzb:nzt+1) )
          var_d = 0.0_wp
       ENDIF

       spectra_initialized = .TRUE.

    END SUBROUTINE spectra_init



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check spectra related quantities
!------------------------------------------------------------------------------!
    SUBROUTINE spectra_check_parameters

       USE control_parameters,                                                 &
           ONLY:  averaging_interval, message_string, skip_time_data_output

       IMPLICIT NONE

!
!--    Check the average interval
       IF ( averaging_interval_sp == 9999999.9_wp )  THEN
          averaging_interval_sp = averaging_interval
       ENDIF

       IF ( averaging_interval_sp > dt_dosp )  THEN
          WRITE( message_string, * )  'averaging_interval_sp = ',              &
                averaging_interval_sp, ' must be <= dt_dosp = ', dt_dosp
          CALL message( 'spectra_check_parameters', 'PA0087', 1, 2, 0, 6, 0 )
       ENDIF

!
!--    Set the default skip time interval for data output, if necessary
       IF ( skip_time_dosp == 9999999.9_wp )                                   &
                                          skip_time_dosp = skip_time_data_output

    END SUBROUTINE spectra_check_parameters



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Header output for spectra
!>
!> @todo Output of netcdf data format and compression level
!------------------------------------------------------------------------------!
    SUBROUTINE spectra_header ( io )

       USE control_parameters,                                                 &
           ONLY:  dt_averaging_input_pr

!       USE netcdf_interface,                                                  &
!           ONLY:  netcdf_data_format_string, netcdf_deflate

       IMPLICIT NONE

!       CHARACTER (LEN=40) ::  output_format       !< internal string

       INTEGER(iwp) ::  i                         !< internal counter
       INTEGER(iwp), INTENT(IN) ::  io            !< Unit of the output file

!
!--    Spectra output
       IF ( dt_dosp /= 9999999.9_wp )  THEN
          WRITE ( io, 1 )

!          output_format = netcdf_data_format_string
!          IF ( netcdf_deflate == 0 )  THEN
!             WRITE ( io, 2 )  output_format
!          ELSE
!             WRITE ( io, 3 )  TRIM( output_format ), netcdf_deflate
!          ENDIF
          WRITE ( io, 2 )  'see profiles or other quantities'
          WRITE ( io, 4 )  dt_dosp
          IF ( skip_time_dosp /= 0.0_wp )  WRITE ( io, 5 )  skip_time_dosp
          WRITE ( io, 6 )  ( data_output_sp(i), i = 1,10 ),     &
                           ( spectra_direction(i), i = 1,10 ),  &
                           ( comp_spectra_level(i), i = 1,100 ), &
                           averaging_interval_sp, dt_averaging_input_pr
       ENDIF

     1 FORMAT ('    Spectra:')
     2 FORMAT ('       Output format: ',A/)
!     3 FORMAT ('       Output format: ',A, '   compressed with level: ',I1/)
     4 FORMAT ('       Output every ',F7.1,' s'/)
     5 FORMAT ('       No output during initial ',F8.2,' s')
     6 FORMAT ('       Arrays:     ', 10(A5,',')/                         &
               '       Directions: ', 10(A5,',')/                         &
               '       height levels  k = ', 20(I3,',')/                  &
               '                          ', 20(I3,',')/                  &
               '                          ', 20(I3,',')/                  &
               '                          ', 20(I3,',')/                  &
               '                          ', 19(I3,','),I3,'.'/           &
               '       Time averaged over ', F7.1, ' s,' /                &
               '       Profiles for the time averaging are taken every ', &
                    F6.1,' s')

    END SUBROUTINE spectra_header



    SUBROUTINE calc_spectra

       USE arrays_3d,                                                          &
           ONLY:  d
#if defined( __parallel )
       USE arrays_3d,                                                          &
           ONLY:  tend
#endif

       USE control_parameters,                                                 &
           ONLY:  bc_lr_cyc, bc_ns_cyc, message_string, psolver

       USE cpulog,                                                             &
           ONLY:  cpu_log, log_point

       USE fft_xy,                                                             &
           ONLY:  fft_init

       USE indices,                                                            &
           ONLY:  nxl, nxr, nyn, nys, nzb, nzt

       USE kinds

       USE pegrid,                                                             &
           ONLY:  myid
#if defined( __parallel )
       USE pegrid,                                                             &
           ONLY:  pdims
#endif

       IMPLICIT NONE

       INTEGER(iwp) ::  m  !<
       INTEGER(iwp) ::  pr !<


!
!--    Check if user gave any levels for spectra to be calculated
       IF ( comp_spectra_level(1) == 999999 )  RETURN

       CALL cpu_log( log_point(30), 'calc_spectra', 'start' )

!
!--    Initialize spectra related quantities
       CALL spectra_init

!
!--    Initialize ffts
       CALL fft_init

!
!--    Reallocate array d in required size
       IF ( psolver(1:9) == 'multigrid' )  THEN
          DEALLOCATE( d )
          ALLOCATE( d(nzb+1:nzt,nys:nyn,nxl:nxr) )
       ENDIF

       m = 1
       DO WHILE ( data_output_sp(m) /= ' '  .AND.  m <= 10 )
!
!--       Transposition from z --> x  ( y --> x in case of a 1d-decomposition
!--       along x)
          IF ( INDEX( spectra_direction(m), 'x' ) /= 0 )  THEN

!
!--          Calculation of spectra works for cyclic boundary conditions only
             IF ( .NOT. bc_lr_cyc )  THEN

                message_string = 'non-cyclic lateral boundaries along x do'//  &
                                 ' not &  allow calculation of spectra along x'
                CALL message( 'calc_spectra', 'PA0160', 1, 2, 0, 6, 0 )
             ENDIF

             CALL preprocess_spectra( m, pr )

#if defined( __parallel )
             IF ( pdims(2) /= 1 )  THEN
                CALL resort_for_zx( d, tend )
                CALL transpose_zx( tend, d )
             ELSE
                CALL transpose_yxd( d, d )
             ENDIF
             CALL calc_spectra_x( d, m )
#else
             message_string = 'sorry, calculation of spectra in non paral' //  &
                              'lel mode& is still not realized'
             CALL message( 'calc_spectra', 'PA0161', 1, 2, 0, 6, 0 )
#endif

          ENDIF

!
!--       Transposition from z --> y (d is rearranged only in case of a
!--       1d-decomposition along x)
          IF ( INDEX( spectra_direction(m), 'y' ) /= 0 )  THEN

!
!--          Calculation of spectra works for cyclic boundary conditions only
             IF ( .NOT. bc_ns_cyc )  THEN
                IF ( myid == 0 )  THEN
                   message_string = 'non-cyclic lateral boundaries along y' // &
                                    ' do not & allow calculation of spectra' //&
                                    ' along y'
                   CALL message( 'calc_spectra', 'PA0162', 1, 2, 0, 6, 0 )
                ENDIF
                CALL local_stop
             ENDIF

             CALL preprocess_spectra( m, pr )

#if defined( __parallel )
             CALL transpose_zyd( d, d )
             CALL calc_spectra_y( d, m )
#else
             message_string = 'sorry, calculation of spectra in non paral' //  &
                              'lel mode& is still not realized'
             CALL message( 'calc_spectra', 'PA0161', 1, 2, 0, 6, 0 )
#endif

          ENDIF

!
!--       Increase counter for next spectrum
          m = m + 1
          
       ENDDO

!
!--    Increase counter for averaging process in routine plot_spectra
       average_count_sp = average_count_sp + 1

       CALL cpu_log( log_point(30), 'calc_spectra', 'stop' )

    END SUBROUTINE calc_spectra


!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE preprocess_spectra( m, pr )

       USE arrays_3d,                                                          &
           ONLY:  d, pt, q, s, u, v, w

       USE indices,                                                            &
           ONLY:  ngp_2dh, nxl, nxr, nyn, nys, nzb, nzt

       USE kinds

#if defined( __parallel )
#if !defined( __mpifh )
       USE MPI
#endif

       USE pegrid,                                                             &
           ONLY:  collective_wait, comm2d, ierr
#endif

       USE statistics,                                                         &
           ONLY:  hom


       IMPLICIT NONE

#if defined( __parallel )
#if defined( __mpifh )
       INCLUDE "mpif.h"
#endif
#endif

       INTEGER(iwp) :: i  !<
       INTEGER(iwp) :: j  !<
       INTEGER(iwp) :: k  !<
       INTEGER(iwp) :: m  !<
       INTEGER(iwp) :: pr !<

       REAL(wp), DIMENSION(nzb:nzt+1) :: var_d_l

       SELECT CASE ( TRIM( data_output_sp(m) ) )
          
       CASE ( 'u' )
          pr = 1
          d(nzb+1:nzt,nys:nyn,nxl:nxr) = u(nzb+1:nzt,nys:nyn,nxl:nxr)
       
       CASE ( 'v' )
          pr = 2
          d(nzb+1:nzt,nys:nyn,nxl:nxr) = v(nzb+1:nzt,nys:nyn,nxl:nxr)
       
       CASE ( 'w' )
          pr = 3
          d(nzb+1:nzt,nys:nyn,nxl:nxr) = w(nzb+1:nzt,nys:nyn,nxl:nxr)
       
       CASE ( 'theta' )
          pr = 4
          d(nzb+1:nzt,nys:nyn,nxl:nxr) = pt(nzb+1:nzt,nys:nyn,nxl:nxr)
       
       CASE ( 'q' )
          pr = 41
          d(nzb+1:nzt,nys:nyn,nxl:nxr) = q(nzb+1:nzt,nys:nyn,nxl:nxr)
          
       CASE ( 's' )
          pr = 117
          d(nzb+1:nzt,nys:nyn,nxl:nxr) = s(nzb+1:nzt,nys:nyn,nxl:nxr)
       
       CASE DEFAULT
!
!--       The DEFAULT case is reached either if the parameter data_output_sp(m)
!--       contains a wrong character string or if the user has coded a special
!--       case in the user interface. There, the subroutine user_spectra
!--       checks which of these two conditions applies.
          CALL user_spectra( 'preprocess', m, pr )
          
       END SELECT

!
!--    Subtract horizontal mean from the array, for which spectra have to be
!--    calculated. Moreover, calculate variance of the respective quantitiy, 
!--    later used for normalizing spectra output. 
       var_d_l(:) = 0.0_wp
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                d(k,j,i)   = d(k,j,i) - hom(k,1,pr,0)
                var_d_l(k) = var_d_l(k) + d(k,j,i) * d(k,j,i)
             ENDDO
          ENDDO
       ENDDO
!
!--    Compute total variance from local variances
       var_d(:) = 0.0_wp
#if defined( __parallel ) 
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( var_d_l(0), var_d(0), nzt+1-nzb, MPI_REAL, MPI_SUM, &
                           comm2d, ierr )
#else
       var_d(:) = var_d_l(:)
#endif
       var_d(:) = var_d(:) / ngp_2dh(0)

    END SUBROUTINE preprocess_spectra


!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
#if defined( __parallel )
    SUBROUTINE calc_spectra_x( ddd, m )

       USE fft_xy,                                                             &
           ONLY:  fft_x_1d

       USE grid_variables,                                                     &
           ONLY:  dx

       USE indices,                                                            &
           ONLY:  nx, ny

       USE kinds

#if !defined( __mpifh )
       USE MPI
#endif

       USE pegrid,                                                             &
           ONLY:  comm2d, ierr, myid

       USE transpose_indices,                                                  &
           ONLY:  nyn_x, nys_x, nzb_x, nzt_x


       IMPLICIT NONE

#if defined( __mpifh )
       INCLUDE "mpif.h"
#endif

       INTEGER(iwp) ::  i         !<
       INTEGER(iwp) ::  j         !<
       INTEGER(iwp) ::  k         !<
       INTEGER(iwp) ::  m         !<
       INTEGER(iwp) ::  n         !<

       REAL(wp) ::  exponent     !<
       REAL(wp) ::  sum_spec_dum !< wavenumber-integrated spectrum
    
       REAL(wp), DIMENSION(0:nx) ::  work !<
    
       REAL(wp), DIMENSION(0:nx/2) ::  sums_spectra_l !<
    
       REAL(wp), DIMENSION(0:nx/2,100) ::  sums_spectra !<
    
       REAL(wp), DIMENSION(0:nx,nys_x:nyn_x,nzb_x:nzt_x) ::  ddd !<

!
!--    Exponent for geometric average
       exponent = 1.0_wp / ( ny + 1.0_wp )

!
!--    Loop over all levels defined by the user
       n = 1
       DO WHILE ( comp_spectra_level(n) /= 999999  .AND.  n <= 100 )

          k = comp_spectra_level(n)

!
!--       Calculate FFT only if the corresponding level is situated on this PE
          IF ( k >= nzb_x  .AND.  k <= nzt_x )  THEN
          
             DO  j = nys_x, nyn_x

                work = ddd(0:nx,j,k)
                CALL fft_x_1d( work, 'forward' )

                ddd(0,j,k) = dx * work(0)**2
                DO  i = 1, nx/2
                   ddd(i,j,k) = dx * ( work(i)**2 + work(nx+1-i)**2 )
                ENDDO

             ENDDO

!
!--          Local sum and geometric average of these spectra
!--          (WARNING: no global sum should be performed, because floating
!--          point overflow may occur)
             DO  i = 0, nx/2

                sums_spectra_l(i) = 1.0_wp
                DO  j = nys_x, nyn_x
                   sums_spectra_l(i) = sums_spectra_l(i) * ddd(i,j,k)**exponent
                ENDDO

             ENDDO
          
          ELSE

             sums_spectra_l = 1.0_wp

          ENDIF

!
!--       Global sum of spectra on PE0 (from where they are written on file)
          sums_spectra(:,n) = 0.0_wp
#if defined( __parallel )   
          CALL MPI_BARRIER( comm2d, ierr )  ! Necessary?
          CALL MPI_REDUCE( sums_spectra_l(0), sums_spectra(0,n), nx/2+1,       &
                           MPI_REAL, MPI_PROD, 0, comm2d, ierr )
#else
          sums_spectra(:,n) = sums_spectra_l
#endif
!
!--       Normalize spectra by variance
          sum_spec_dum = SUM( sums_spectra(1:nx/2,n) )

          IF ( sum_spec_dum /= 0.0_wp )  THEN
             sums_spectra(1:nx/2,n) = sums_spectra(1:nx/2,n) *                 &
                                      var_d(k) / sum_spec_dum
          ENDIF
          n = n + 1

       ENDDO
       n = n - 1

       IF ( myid == 0 )  THEN
!
!--       Sum of spectra for later averaging (see routine data_output_spectra)
          DO  i = 1, nx/2
             DO k = 1, n
                spectrum_x(i,k,m) = spectrum_x(i,k,m) + sums_spectra(i,k)
             ENDDO
          ENDDO

       ENDIF
!
!--    n_sp_x is needed by data_output_spectra_x
       n_sp_x = n

    END SUBROUTINE calc_spectra_x
#endif


!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
#if defined( __parallel )
    SUBROUTINE calc_spectra_y( ddd, m )

       USE fft_xy,                                                             &
           ONLY:  fft_y_1d

       USE grid_variables,                                                     &
           ONLY:  dy

       USE indices,                                                            &
           ONLY:  nx, ny

       USE kinds

#if !defined( __mpifh )
       USE MPI
#endif

       USE pegrid,                                                             &
           ONLY:  comm2d, ierr, myid

       USE transpose_indices,                                                  &
           ONLY:  nxl_yd, nxr_yd, nzb_yd, nzt_yd


       IMPLICIT NONE

#if defined( __mpifh )
       INCLUDE "mpif.h"
#endif

       INTEGER(iwp) ::  i         !<
       INTEGER(iwp) ::  j         !<
       INTEGER(iwp) ::  k         !<
       INTEGER(iwp) ::  m         !<
       INTEGER(iwp) ::  n         !<

       REAL(wp) ::  exponent !<
       REAL(wp) ::  sum_spec_dum !< wavenumber-integrated spectrum
    
       REAL(wp), DIMENSION(0:ny) ::  work !<
    
       REAL(wp), DIMENSION(0:ny/2) ::  sums_spectra_l !<
    
       REAL(wp), DIMENSION(0:ny/2,100) ::  sums_spectra !<
    
       REAL(wp), DIMENSION(0:ny,nxl_yd:nxr_yd,nzb_yd:nzt_yd) :: ddd !<


!
!--    Exponent for geometric average
       exponent = 1.0_wp / ( nx + 1.0_wp )

!
!--    Loop over all levels defined by the user
       n = 1
       DO WHILE ( comp_spectra_level(n) /= 999999  .AND.  n <= 100 )

          k = comp_spectra_level(n)

!
!--       Calculate FFT only if the corresponding level is situated on this PE
          IF ( k >= nzb_yd  .AND.  k <= nzt_yd )  THEN
          
             DO  i = nxl_yd, nxr_yd

                work = ddd(0:ny,i,k)
                CALL fft_y_1d( work, 'forward' )

                ddd(0,i,k) = dy * work(0)**2
                DO  j = 1, ny/2
                   ddd(j,i,k) = dy * ( work(j)**2 + work(ny+1-j)**2 )
                ENDDO

             ENDDO

!
!--          Local sum and geometric average of these spectra
!--          (WARNING: no global sum should be performed, because floating
!--          point overflow may occur)
             DO  j = 0, ny/2

                sums_spectra_l(j) = 1.0_wp
                DO  i = nxl_yd, nxr_yd
                   sums_spectra_l(j) = sums_spectra_l(j) * ddd(j,i,k)**exponent
                ENDDO

             ENDDO
          
          ELSE

             sums_spectra_l = 1.0_wp

          ENDIF

!
!--       Global sum of spectra on PE0 (from where they are written on file)
          sums_spectra(:,n) = 0.0_wp
#if defined( __parallel )   
          CALL MPI_BARRIER( comm2d, ierr )  ! Necessary?
          CALL MPI_REDUCE( sums_spectra_l(0), sums_spectra(0,n), ny/2+1,       &
                           MPI_REAL, MPI_PROD, 0, comm2d, ierr )
#else
          sums_spectra(:,n) = sums_spectra_l
#endif
!
!--       Normalize spectra by variance
          sum_spec_dum = SUM( sums_spectra(1:ny/2,n) )
          IF ( sum_spec_dum /= 0.0_wp )  THEN
             sums_spectra(1:ny/2,n) = sums_spectra(1:ny/2,n) *                 &
                                      var_d(k) / sum_spec_dum
          ENDIF
          n = n + 1

       ENDDO
       n = n - 1


       IF ( myid == 0 )  THEN
!
!--       Sum of spectra for later averaging (see routine data_output_spectra)
          DO  j = 1, ny/2
             DO k = 1, n
                spectrum_y(j,k,m) = spectrum_y(j,k,m) + sums_spectra(j,k)
             ENDDO
          ENDDO

       ENDIF

!
!--    n_sp_y is needed by data_output_spectra_y
       n_sp_y = n

    END SUBROUTINE calc_spectra_y
#endif

 END MODULE spectra_mod
