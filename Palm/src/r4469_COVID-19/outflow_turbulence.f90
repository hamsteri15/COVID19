!> @file outflow_turbulence.f90
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
! $Id: outflow_turbulence.f90 3241 2018-09-12 15:02:00Z raasch $
! Corrected "Former revisions" section
! 
! 3241 2018-09-12 15:02:00Z raasch
! unused variables removed
!
! 2050 2016-11-08 15:00:55Z gronemeier
! Initial version
! 
!
! Description:
! ------------
!> Routine based on inflow_turbulence.f90. Copies values of 3d data from a 2d
!> vertical source plane (defined by outflow_source_plane) to the outflow
!> boundary.
!------------------------------------------------------------------------------!
 SUBROUTINE outflow_turbulence

    USE arrays_3d,                                                             &
        ONLY:  e, pt, q, s, u, v, w

    USE control_parameters,                                                    &
        ONLY:  humidity, passive_scalar, outflow_source_plane

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point

    USE grid_variables,                                                        &
        ONLY:  ddx

    USE indices,                                                               &
        ONLY:  nbgp, nx, nxr, nyn, nys, nyng, nysg, nzb, nzt

    USE kinds

    USE pegrid!,                                                                &
        !ONLY:  comm1dx, id_outflow, id_outflow_source, ierr, myidx, status

    IMPLICIT NONE

    INTEGER(iwp) ::  i        !< loop index
    INTEGER(iwp) ::  j        !< loop index
    INTEGER(iwp) ::  k        !< loop index
    INTEGER(iwp) ::  l        !< loop index
    INTEGER(iwp) ::  ngp_ofv  !< number of grid points stored in outflow_val

    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,7,nbgp) ::                         &
       outflow_val            !< values to be copied to the outflow boundary


    CALL cpu_log( log_point(56), 'outflow_turbulence', 'start' )
!
!-- Get number of grid points stored in outflow_val
    ngp_ofv = ( nzt - nzb + 2 ) * ( nyn - nys + 1 + 2 * nbgp ) * 7 * nbgp
!
!-- Get position of the source plane inside the palm grid
    i = outflow_source_plane * ddx
!
!-- Use instantaneous values instead of averaged profiles
#if defined( __parallel )
    IF ( myidx == id_outflow_source )  THEN
       DO  l = 1, nbgp
          DO  j = nysg, nyng
             DO  k = nzb, nzt + 1

                outflow_val(k,j,1,l) = u(k,j,i)
                outflow_val(k,j,2,l) = v(k,j,i)
                outflow_val(k,j,3,l) = w(k,j,i)
                outflow_val(k,j,4,l) = pt(k,j,i)
                outflow_val(k,j,5,l) = e(k,j,i)
                IF ( humidity  )                                               &
                   outflow_val(k,j,6,l) = q(k,j,i)
                IF ( passive_scalar )                                          &
                   outflow_val(k,j,7,l) = s(k,j,i)

            ENDDO
          ENDDO
          i = i + 1
       ENDDO

    ENDIF
#else
    DO  l = 1, nbgp
       DO  j = nysg, nyng
          DO  k = nzb, nzt+1

             outflow_val(k,j,1,l) = u(k,j,i)
             outflow_val(k,j,2,l) = v(k,j,i)
             outflow_val(k,j,3,l) = w(k,j,i)
             outflow_val(k,j,4,l) = pt(k,j,i)
             outflow_val(k,j,5,l) = e(k,j,i)
             IF ( humidity  )                                                  &
                outflow_val(k,j,6,l) = q(k,j,i)
             IF ( passive_scalar )                                             &
                outflow_val(k,j,7,l) = s(k,j,i)

          ENDDO
       ENDDO
       i = i + 1
    ENDDO
#endif

!
!-- For parallel runs, send the values to the respective outflow PE
#if defined( __parallel )
    IF ( myidx == id_outflow_source  .AND.  myidx /= id_outflow )  THEN

       CALL MPI_SEND( outflow_val(nzb,nysg,1,1), ngp_ofv, MPI_REAL,            &
                      id_outflow, 1, comm1dx, ierr )

    ELSEIF ( myidx /= id_outflow_source  .AND.  myidx == id_outflow )  THEN

       outflow_val = 0.0_wp
       CALL MPI_RECV( outflow_val(nzb,nysg,1,1), ngp_ofv, MPI_REAL,            &
                      id_outflow_source, 1, comm1dx, status, ierr )

    ENDIF
#endif

!
!-- Copy values to the outflow
    IF ( nxr == nx )  THEN

       DO  j = nysg, nyng
          DO  k = nzb, nzt + 1

             u(k,j,nx+1:nx+nbgp)  = outflow_val(k,j,1,1:nbgp)
             v(k,j,nx+1:nx+nbgp)  = outflow_val(k,j,2,1:nbgp)
             w(k,j,nx+1:nx+nbgp)  = outflow_val(k,j,3,1:nbgp)
             pt(k,j,nx+1:nx+nbgp) = outflow_val(k,j,4,1:nbgp)
             e(k,j,nx+1:nx+nbgp)  = outflow_val(k,j,5,1:nbgp)
             e(k,j,nx+1:nx+nbgp)  = MAX( e(k,j,nx+1:nx+nbgp), 0.0_wp )

             IF ( humidity )                                                   &
                q(k,j,nx+1:nx+nbgp)  = outflow_val(k,j,6,1:nbgp)
             IF ( passive_scalar )                                             &
                s(k,j,nx+1:nx+nbgp)  = outflow_val(k,j,7,1:nbgp)

          ENDDO
       ENDDO


    ENDIF

    CALL cpu_log( log_point(56), 'outflow_turbulence', 'stop' )

 END SUBROUTINE outflow_turbulence
