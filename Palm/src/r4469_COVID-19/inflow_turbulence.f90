!> @file inflow_turbulence.f90
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
! $Id: inflow_turbulence.f90 4429 2020-02-27 15:24:30Z raasch $
! bugfix: cpp-directives added for serial mode
! 
! 4360 2020-01-07 11:25:50Z suehring
! use y_shift instead of old parameter recycling_yshift
! 
! 4297 2019-11-21 10:37:50Z oliver.maas
! changed recycling_yshift so that the y-shift can be a multiple of PE
! instead of y-shift of a half domain width
! 
! 4183 2019-08-23 07:33:16Z oliver.maas
! simplified steering of recycling of absolute values by initialization 
! parameter recycling_method_for_thermodynamic_quantities
! 
! 4182 2019-08-22 15:20:23Z scharf
! Corrected "Former revisions" section
! 
! 4172 2019-08-20 11:55:33Z oliver.maas
! added optional recycling of absolute values for pt and q
! 
! 3655 2019-01-07 16:51:22Z knoop
! Corrected "Former revisions" section
!
! Initial version (2008/03/07)
!
! Description:
! ------------
!> Imposing turbulence at the respective inflow using the turbulence
!> recycling method of Kataoka and Mizuno (2002).
!------------------------------------------------------------------------------!
 SUBROUTINE inflow_turbulence
 

    USE arrays_3d,                                                             &
        ONLY:  e, inflow_damping_factor, mean_inflow_profiles, pt, q, s, u, v, w
        
#if defined( __parallel )
    USE control_parameters,                                                    &
        ONLY:  humidity, passive_scalar, recycling_plane, y_shift,             &
               recycling_method_for_thermodynamic_quantities
#else
    USE control_parameters,                                                    &
        ONLY:  humidity, passive_scalar, recycling_plane,                      &
               recycling_method_for_thermodynamic_quantities
#endif
        
    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point
        
    USE indices,                                                               &
        ONLY:  nbgp, nxl, ny, nyn, nys, nyng, nysg, nzb, nzt
        
    USE kinds
    
    USE pegrid


    IMPLICIT NONE
    
    INTEGER(iwp) ::  i        !< loop index
    INTEGER(iwp) ::  j        !< loop index
    INTEGER(iwp) ::  k        !< loop index
    INTEGER(iwp) ::  l        !< loop index
    INTEGER(iwp) ::  ngp_ifd  !< number of grid points stored in avpr
    INTEGER(iwp) ::  ngp_pr   !< number of grid points stored in inflow_dist
#if defined( __parallel )
    INTEGER(iwp) ::  next     !< ID of receiving PE for y-shift
    INTEGER(iwp) ::  prev     !< ID of sending PE for y-shift
#endif

    REAL(wp), DIMENSION(nzb:nzt+1,7,nbgp)           ::                         &
       avpr               !< stores averaged profiles at recycling plane
    REAL(wp), DIMENSION(nzb:nzt+1,7,nbgp)           ::                         &
       avpr_l             !< auxiliary variable to calculate avpr
    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,7,nbgp) ::                         &
       inflow_dist        !< turbulence signal of vars, added at inflow boundary
#if defined( __parallel )
    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,7,nbgp) ::                         &
       local_inflow_dist  !< auxiliary variable for inflow_dist, used for y-shift
#endif
    
    CALL cpu_log( log_point(40), 'inflow_turbulence', 'start' )
    
!
!-- Carry out spanwise averaging in the recycling plane
    avpr_l = 0.0_wp
    ngp_pr = ( nzt - nzb + 2 ) * 7 * nbgp
    ngp_ifd = ngp_pr * ( nyn - nys + 1 + 2 * nbgp )

!
!-- First, local averaging within the recycling domain
    i = recycling_plane

#if defined( __parallel )
    IF ( myidx == id_recycling )  THEN
       
       DO  l = 1, nbgp
          DO  j = nys, nyn
             DO  k = nzb, nzt + 1

                avpr_l(k,1,l) = avpr_l(k,1,l) + u(k,j,i)
                avpr_l(k,2,l) = avpr_l(k,2,l) + v(k,j,i)
                avpr_l(k,3,l) = avpr_l(k,3,l) + w(k,j,i)
                avpr_l(k,4,l) = avpr_l(k,4,l) + pt(k,j,i)
                avpr_l(k,5,l) = avpr_l(k,5,l) + e(k,j,i)
                IF ( humidity )                                                &
                   avpr_l(k,6,l) = avpr_l(k,6,l) + q(k,j,i)
                IF ( passive_scalar )                                          &
                   avpr_l(k,7,l) = avpr_l(k,7,l) + s(k,j,i)

             ENDDO
          ENDDO
          i = i + 1
       ENDDO

    ENDIF
!
!-- Now, averaging over all PEs
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( avpr_l(nzb,1,1), avpr(nzb,1,1), ngp_pr, MPI_REAL,      &
                        MPI_SUM, comm2d, ierr )

#else
    DO  l = 1, nbgp
       DO  j = nys, nyn
          DO  k = nzb, nzt + 1

             avpr_l(k,1,l) = avpr_l(k,1,l) + u(k,j,i)
             avpr_l(k,2,l) = avpr_l(k,2,l) + v(k,j,i)
             avpr_l(k,3,l) = avpr_l(k,3,l) + w(k,j,i)
             avpr_l(k,4,l) = avpr_l(k,4,l) + pt(k,j,i)
             avpr_l(k,5,l) = avpr_l(k,5,l) + e(k,j,i)
             IF ( humidity )                                                   &
                avpr_l(k,6,l) = avpr_l(k,6,l) + q(k,j,i)
             IF ( passive_scalar )                                             &
                avpr_l(k,7,l) = avpr_l(k,7,l) + s(k,j,i)

          ENDDO
       ENDDO
       i = i + 1 
    ENDDO
    
    avpr = avpr_l
#endif

    avpr = avpr / ( ny + 1 )
!
!-- Calculate the disturbances at the recycling plane
!-- for recycling of absolute quantities, the disturbance is defined as the absolute value
!-- (and not as the deviation from the mean profile)
    i = recycling_plane

#if defined( __parallel )
    IF ( myidx == id_recycling )  THEN
       DO  l = 1, nbgp
          DO  j = nysg, nyng
             DO  k = nzb, nzt + 1
                inflow_dist(k,j,1,l) = u(k,j,i+1) - avpr(k,1,l)
                inflow_dist(k,j,2,l) = v(k,j,i)   - avpr(k,2,l)
                inflow_dist(k,j,3,l) = w(k,j,i)   - avpr(k,3,l)
                IF ( TRIM( recycling_method_for_thermodynamic_quantities )     &
                   == 'turbulent_fluctuation' )  THEN
                   inflow_dist(k,j,4,l) = pt(k,j,i) - avpr(k,4,l)
                ELSEIF ( TRIM( recycling_method_for_thermodynamic_quantities ) &
                   == 'absolute_value' )  THEN
                   inflow_dist(k,j,4,l) = pt(k,j,i)
                ENDIF
                inflow_dist(k,j,5,l) = e(k,j,i)   - avpr(k,5,l)
                IF ( humidity ) THEN
                   IF ( TRIM( recycling_method_for_thermodynamic_quantities )  &
                      == 'turbulent_fluctuation' )  THEN
                      inflow_dist(k,j,6,l) = q(k,j,i) - avpr(k,6,l)
                   ELSEIF ( TRIM( recycling_method_for_thermodynamic_quantities )  &
                      == 'absolute_value' )  THEN
                      inflow_dist(k,j,6,l) = q(k,j,i)
                   ENDIF
                ENDIF
                IF ( passive_scalar )                                          &
                   inflow_dist(k,j,7,l) = s(k,j,i) - avpr(k,7,l)
            ENDDO
          ENDDO
          i = i + 1
       ENDDO

    ENDIF
#else
    DO  l = 1, nbgp
       DO  j = nysg, nyng
          DO  k = nzb, nzt+1
             inflow_dist(k,j,1,l) = u(k,j,i+1) - avpr(k,1,l)
             inflow_dist(k,j,2,l) = v(k,j,i)   - avpr(k,2,l)
             inflow_dist(k,j,3,l) = w(k,j,i)   - avpr(k,3,l)
             IF ( TRIM( recycling_method_for_thermodynamic_quantities )        &
                   == 'turbulent_fluctuation' )  THEN
                inflow_dist(k,j,4,l) = pt(k,j,i) - avpr(k,4,l)
             ELSEIF ( TRIM( recycling_method_for_thermodynamic_quantities )    &
                   == 'absolute_value' )  THEN
                inflow_dist(k,j,4,l) = pt(k,j,i)
             ENDIF
             inflow_dist(k,j,5,l) = e(k,j,i)   - avpr(k,5,l)
             IF ( humidity )  THEN
                IF ( TRIM( recycling_method_for_thermodynamic_quantities )     &
                      == 'turbulent_fluctuation' )  THEN
                   inflow_dist(k,j,6,l) = q(k,j,i) - avpr(k,6,l)
                ELSEIF ( TRIM( recycling_method_for_thermodynamic_quantities ) &
                      == 'absolute_value' )  THEN
                   inflow_dist(k,j,6,l) = q(k,j,i)
                ENDIF
             ENDIF
             IF ( passive_scalar )                                             &
                inflow_dist(k,j,7,l) = s(k,j,i) - avpr(k,7,l)
              
          ENDDO
       ENDDO
       i = i + 1
    ENDDO
#endif

!
!-- For parallel runs, send the disturbances to the respective inflow PE
#if defined( __parallel )
    IF ( myidx == id_recycling  .AND.  myidx /= id_inflow )  THEN

       CALL MPI_SEND( inflow_dist(nzb,nysg,1,1), ngp_ifd, MPI_REAL,            &
                      id_inflow, 1, comm1dx, ierr )

    ELSEIF ( myidx /= id_recycling  .AND.  myidx == id_inflow )  THEN

       inflow_dist = 0.0_wp
       CALL MPI_RECV( inflow_dist(nzb,nysg,1,1), ngp_ifd, MPI_REAL,            &
                      id_recycling, 1, comm1dx, status, ierr )

    ENDIF

!
!-- y-shift for inflow_dist
!-- Shift inflow_dist in positive y direction by a number of
!-- PEs equal to y_shift
    IF ( ( y_shift /= 0 ) .AND. myidx == id_inflow ) THEN

!
!--    Calculate the ID of the PE which sends data to this PE (prev) and of the
!--    PE which receives data from this PE (next).
       prev = MODULO(myidy - y_shift , pdims(2))
       next = MODULO(myidy + y_shift , pdims(2))
       
       local_inflow_dist = 0.0_wp

       CALL MPI_SENDRECV( inflow_dist(nzb,nysg,1,1), ngp_ifd, MPI_REAL,        &
                          next, 1, local_inflow_dist(nzb,nysg,1,1), ngp_ifd,   &
                          MPI_REAL, prev, 1, comm1dy, status, ierr )

       inflow_dist = local_inflow_dist

    ENDIF

#endif

!
!-- Add the disturbance at the inflow
    IF ( nxl == 0 )  THEN

       DO  j = nysg, nyng
          DO  k = nzb, nzt + 1

             u(k,j,-nbgp+1:0) = mean_inflow_profiles(k,1) +                    &
                        inflow_dist(k,j,1,1:nbgp) * inflow_damping_factor(k)
             v(k,j,-nbgp:-1)  = mean_inflow_profiles(k,2) +                    &
                        inflow_dist(k,j,2,1:nbgp) * inflow_damping_factor(k)
             w(k,j,-nbgp:-1)  =                                                &
                        inflow_dist(k,j,3,1:nbgp) * inflow_damping_factor(k)
             IF ( TRIM( recycling_method_for_thermodynamic_quantities )        &
                   == 'turbulent_fluctuation' )  THEN
                pt(k,j,-nbgp:-1) = mean_inflow_profiles(k,4) +                 &
                inflow_dist(k,j,4,1:nbgp) * inflow_damping_factor(k)
             ELSEIF ( TRIM( recycling_method_for_thermodynamic_quantities )    &
                   == 'absolute_value' )  THEN
                pt(k,j,-nbgp:-1) = inflow_dist(k,j,4,1:nbgp)
             ENDIF
             e(k,j,-nbgp:-1)  = mean_inflow_profiles(k,5) +                    &
                        inflow_dist(k,j,5,1:nbgp) * inflow_damping_factor(k)
             e(k,j,-nbgp:-1)  = MAX( e(k,j,-nbgp:-1), 0.0_wp )
             IF ( humidity )  THEN
                IF ( TRIM( recycling_method_for_thermodynamic_quantities )     &
                      == 'turbulent_fluctuation' )  THEN
                   q(k,j,-nbgp:-1)  = mean_inflow_profiles(k,6) +              & 
                      inflow_dist(k,j,6,1:nbgp) * inflow_damping_factor(k)
                ELSEIF ( TRIM( recycling_method_for_thermodynamic_quantities ) &
                      == 'absolute_value' )  THEN
                   q(k,j,-nbgp:-1)  = inflow_dist(k,j,6,1:nbgp)
                ENDIF
             ENDIF
             IF ( passive_scalar )                                             &
                s(k,j,-nbgp:-1)  = mean_inflow_profiles(k,7) +                 &
                        inflow_dist(k,j,7,1:nbgp) * inflow_damping_factor(k)
                        
          ENDDO
       ENDDO

    ENDIF


    CALL cpu_log( log_point(40), 'inflow_turbulence', 'stop' )


 END SUBROUTINE inflow_turbulence
