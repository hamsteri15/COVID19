!> @file print_1d.f90
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
! $Id: print_1d.f90 4360 2020-01-07 11:25:50Z suehring $
! Corrected "Former revisions" section
! 
! 3655 2019-01-07 16:51:22Z knoop
! Renamed output variables
!
! Revision 1.1  1997/09/19 07:45:22  raasch
! Initial revision
!
!
! Description:
! ------------
!> List output of 1D-profiles.
!------------------------------------------------------------------------------!
 SUBROUTINE print_1d
 

    USE arrays_3d,                                                             &
        ONLY:  zu, zw

    USE control_parameters,                                                    &
        ONLY:  run_description_header, simulated_time_chr

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point

    USE indices,                                                               &
        ONLY:  nzb, nzt

    USE kinds

    USE pegrid

    USE statistics,                                                            &
        ONLY:  flow_statistics_called, hom, region, statistic_regions

    IMPLICIT NONE


    CHARACTER (LEN=20) ::  period_chr  !<

    INTEGER(iwp) ::  k   !<
    INTEGER(iwp) ::  sr  !<


!
!-- If required, compute statistics.
    IF ( .NOT. flow_statistics_called )  CALL flow_statistics

!
!-- Flow_statistics has its own cpu-time measuring.
    CALL cpu_log( log_point(18), 'print_1d', 'start' )

    IF ( myid == 0 )  THEN
!
!--    Open file for list output of profiles.
       CALL check_open( 16 )

!
!--    Prepare header.
       period_chr = ' no time-average!'

!
!--    Output for the total domain (and each subregion, if applicable).
       DO  sr = 0, statistic_regions
!
!--       Write header.
          WRITE ( 16, 112 )
          WRITE ( 16, 100 )  TRIM( run_description_header ) // '    ' // &
                             TRIM( region( sr ) ), TRIM( period_chr ), 'uv'
          WRITE ( 16, 105 )  TRIM( simulated_time_chr )
          WRITE ( 16, 111 )

!
!--       Output of values on the scalar levels.
          WRITE ( 16, 120 )
          WRITE ( 16, 111 )
          DO  k = nzt+1, nzb, -1
             WRITE ( 16, 121)  k, zu(k), hom(k,1,1,sr),           &
                               hom(k,1,1,sr) - hom(k,1,5,sr),     &
                               hom(k,1,2,sr),                     &
                               hom(k,1,2,sr) - hom(k,1,6,sr),     &
                               hom(k,1,4,sr),                     &
                               hom(k,1,4,sr) - hom(k,1,7,sr),     &
                               hom(k,1,8,sr), hom(k,1,9,sr),      &
                               hom(k,1,10,sr), hom(k,1,11,sr), zu(k), k
          ENDDO
          WRITE ( 16, 111 )
          WRITE ( 16, 120 )
          WRITE ( 16, 111 )

!
!--       Output of values on the w-levels.
          WRITE ( 16, 112 )
          WRITE ( 16, 100 )  TRIM( run_description_header ) // '    ' // &
                             TRIM( region( sr ) ), TRIM( period_chr ), 'w'
          WRITE ( 16, 105 )  TRIM( simulated_time_chr )
          WRITE ( 16, 111 )

          WRITE ( 16, 130 )
          WRITE ( 16, 111 )
          DO  k = nzt+1, nzb, -1
             WRITE ( 16, 131)  k, zw(k), hom(k,1,16,sr),            &
                               hom(k,1,18,sr), hom(k,1,12,sr), &
                               hom(k,1,19,sr), hom(k,1,14,sr), &
                               hom(k,1,20,sr), zw(k), k
          ENDDO
          WRITE ( 16, 111 )
          WRITE ( 16, 130 )
          WRITE ( 16, 111 )

       ENDDO

    ENDIF

    CALL cpu_log( log_point(18), 'print_1d', 'stop' )

!
!-- Formats.
100 FORMAT (1X,A/1X,10('-')/ &
            ' Horizontally',A,' averaged profiles on the ',A,'-level')
105 FORMAT (' Time: ',A)
111 FORMAT (1X,131('-'))
112 FORMAT (/)
120 FORMAT ('   k     zu      u     du     v     dv     theta dtheta ', &
            ' e      Km    Kh     l      zu      k')
121 FORMAT (1X,I4,1X,F7.1,1X,F6.2,1X,F5.2,1X,F6.2,1X,F5.2,2X,F6.2,1X,F5.2, &
            1X,F7.4,1X,F5.2,1X,F5.2,1X,F6.2,1X,F7.1,2X,I4)
130 FORMAT ('   k     zw       w''theta''   wtheta     w''u''       wu       ',&
            '  w''v''       wv        zw      k')
131 FORMAT (1X,I4,1X,F7.1,6(1X,E10.3),1X,F7.1,2X,I4)


 END SUBROUTINE print_1d
