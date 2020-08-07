!> @file advec_s_bc.f90
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
! $Id: advec_s_bc.f90 4429 2020-02-27 15:24:30Z raasch $
! bugfix: cpp-directives added for serial mode
! 
! 4360 2020-01-07 11:25:50Z suehring
! Corrected "Former revisions" section
! 
! 3761 2019-02-25 15:31:42Z raasch
! unused variables removed
! 
! 3655 2019-01-07 16:51:22Z knoop
! nopointer option removed
!
! Revision 1.1  1997/08/29 08:53:46  raasch
! Initial revision
!
!
! Description:
! ------------
!> Advection term for scalar quantities using the Bott-Chlond scheme.
!> Computation in individual steps for each of the three dimensions.
!> Limiting assumptions:
!> So far the scheme has been assuming equidistant grid spacing. As this is not
!> the case in the stretched portion of the z-direction, there dzw(k) is used as
!> a substitute for a constant grid length. This certainly causes incorrect
!> results; however, it is hoped that they are not too apparent for weakly 
!> stretched grids.
!> NOTE: This is a provisional, non-optimised version!
!------------------------------------------------------------------------------!
MODULE advec_s_bc_mod
 

    PRIVATE
    PUBLIC advec_s_bc

    INTERFACE advec_s_bc
       MODULE PROCEDURE advec_s_bc
    END INTERFACE advec_s_bc

 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE advec_s_bc( sk, sk_char )

       USE advection,                                                             &
           ONLY:  aex, bex, dex, eex

       USE arrays_3d,                                                             &
           ONLY:  d, ddzw, dzu, dzw, tend, u, v, w

       USE control_parameters,                                                    &
           ONLY:  dt_3d, bc_pt_t_val, bc_q_t_val, bc_s_t_val, ibc_pt_b, ibc_pt_t, &
                  ibc_q_t, ibc_s_t, message_string, pt_slope_offset,              &
                  sloping_surface, u_gtrans, v_gtrans

       USE cpulog,                                                                &
           ONLY:  cpu_log, log_point_s

       USE grid_variables,                                                        &
           ONLY:  ddx, ddy

       USE indices,                                                               &
           ONLY:  nx, nxl, nxr, nyn, nys, nzb, nzt

       USE kinds

       USE pegrid

       USE statistics,                                                            &
           ONLY:  rmask, statistic_regions, sums_wsts_bc_l


       IMPLICIT NONE

       CHARACTER (LEN=*) ::  sk_char !<

       INTEGER(iwp) ::  i         !<
       INTEGER(iwp) ::  ix        !<
       INTEGER(iwp) ::  j         !<
       INTEGER(iwp) ::  k         !<
       INTEGER(iwp) ::  sr        !<
#if defined( __parallel )
       INTEGER(iwp) ::  ngp       !<
       INTEGER(iwp) ::  type_xz_2 !<
#endif

       REAL(wp) ::  cim    !<
       REAL(wp) ::  cimf   !<
       REAL(wp) ::  cip    !<
       REAL(wp) ::  cipf   !<
       REAL(wp) ::  d_new  !<
       REAL(wp) ::  denomi !< denominator
       REAL(wp) ::  fminus !<
       REAL(wp) ::  fplus  !<
       REAL(wp) ::  f2     !<
       REAL(wp) ::  f4     !<
       REAL(wp) ::  f8     !<
       REAL(wp) ::  f12    !<
       REAL(wp) ::  f24    !<
       REAL(wp) ::  f48    !<
       REAL(wp) ::  f1920  !<
       REAL(wp) ::  im     !<
       REAL(wp) ::  ip     !<
       REAL(wp) ::  m1n    !<
       REAL(wp) ::  m1z    !<
       REAL(wp) ::  m2     !<
       REAL(wp) ::  m3     !<
       REAL(wp) ::  numera !< numerator
       REAL(wp) ::  snenn  !<
       REAL(wp) ::  sterm  !<
       REAL(wp) ::  tendcy !<
       REAL(wp) ::  t1     !<
       REAL(wp) ::  t2     !<

       REAL(wp) ::  fmax(2)   !<
       REAL(wp) ::  fmax_l(2) !<
       
       REAL(wp), DIMENSION(:,:,:), POINTER ::  sk

       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  a0   !<
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  a1   !<
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  a12  !<
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  a2   !<
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  a22  !<
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  immb !<
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  imme !<
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  impb !<
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  impe !<
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ipmb !<
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ipme !<
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ippb !<
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ippe !<
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  m1   !<
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  sw   !<
       
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  sk_p !<

!
!--    Array sk_p requires 2 extra elements for each dimension
       ALLOCATE( sk_p(nzb-2:nzt+3,nys-3:nyn+3,nxl-3:nxr+3) )
       sk_p = 0.0_wp

!
!--    Assign reciprocal values in order to avoid divisions later
       f2    = 0.5_wp
       f4    = 0.25_wp
       f8    = 0.125_wp
       f12   = 0.8333333333333333E-01_wp
       f24   = 0.4166666666666666E-01_wp
       f48   = 0.2083333333333333E-01_wp
       f1920 = 0.5208333333333333E-03_wp

!
!--    Advection in x-direction:

!
!--    Save the quantity to be advected in a local array
!--    add an enlarged boundary in x-direction
       DO  i = nxl-1, nxr+1
          DO  j = nys, nyn
             DO  k = nzb, nzt+1
                sk_p(k,j,i) = sk(k,j,i)
             ENDDO
          ENDDO
       ENDDO
#if defined( __parallel )
       ngp = 2 * ( nzt - nzb + 6 ) * ( nyn - nys + 7 )
       CALL cpu_log( log_point_s(11), 'advec_s_bc:sendrecv', 'start' )
!
!--    Send left boundary, receive right boundary
       CALL MPI_SENDRECV( sk_p(nzb-2,nys-3,nxl+1), ngp, MPI_REAL, pleft,  0,      &
                          sk_p(nzb-2,nys-3,nxr+2), ngp, MPI_REAL, pright, 0,      &
                          comm2d, status, ierr )
!
!--    Send right boundary, receive left boundary
       CALL MPI_SENDRECV( sk_p(nzb-2,nys-3,nxr-2), ngp, MPI_REAL, pright, 1,      &
                          sk_p(nzb-2,nys-3,nxl-3), ngp, MPI_REAL, pleft,  1,      &
                          comm2d, status, ierr )
       CALL cpu_log( log_point_s(11), 'advec_s_bc:sendrecv', 'pause' )
#else

!
!--    Cyclic boundary conditions
       sk_p(:,nys:nyn,nxl-3) = sk_p(:,nys:nyn,nxr-2)
       sk_p(:,nys:nyn,nxl-2) = sk_p(:,nys:nyn,nxr-1)
       sk_p(:,nys:nyn,nxr+2) = sk_p(:,nys:nyn,nxl+1)
       sk_p(:,nys:nyn,nxr+3) = sk_p(:,nys:nyn,nxl+2)
#endif

!
!--    In case of a sloping surface, the additional gridpoints in x-direction
!--    of the temperature field at the left and right boundary of the total
!--    domain must be adjusted by the temperature difference between this distance
       IF ( sloping_surface  .AND.  sk_char == 'pt' )  THEN
          IF ( nxl ==  0 )  THEN
             sk_p(:,nys:nyn,nxl-3) = sk_p(:,nys:nyn,nxl-3) - pt_slope_offset
             sk_p(:,nys:nyn,nxl-2) = sk_p(:,nys:nyn,nxl-2) - pt_slope_offset
          ENDIF
          IF ( nxr == nx )  THEN
             sk_p(:,nys:nyn,nxr+2) = sk_p(:,nys:nyn,nxr+2) + pt_slope_offset
             sk_p(:,nys:nyn,nxr+3) = sk_p(:,nys:nyn,nxr+3) + pt_slope_offset
          ENDIF
       ENDIF

!
!--    Initialise control density
       d = 0.0_wp

!
!--    Determine maxima of the first and second derivative in x-direction
       fmax_l = 0.0_wp
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                numera = ABS( sk_p(k,j,i+1) - 2.0_wp * sk_p(k,j,i) + sk_p(k,j,i-1) )
                denomi  = ABS( sk_p(k,j,i+1) - sk_p(k,j,i-1) )
                fmax_l(1) = MAX( fmax_l(1) , numera )
                fmax_l(2) = MAX( fmax_l(2) , denomi  )
             ENDDO
          ENDDO
       ENDDO
#if defined( __parallel )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( fmax_l, fmax, 2, MPI_REAL, MPI_MAX, comm2d, ierr )
#else
       fmax = fmax_l
#endif

       fmax = 0.04_wp * fmax

!
!--    Allocate temporary arrays
       ALLOCATE( a0(nzb+1:nzt,nxl-1:nxr+1),   a1(nzb+1:nzt,nxl-1:nxr+1),          &
                 a2(nzb+1:nzt,nxl-1:nxr+1),   a12(nzb+1:nzt,nxl-1:nxr+1),         &
                 a22(nzb+1:nzt,nxl-1:nxr+1),  immb(nzb+1:nzt,nxl-1:nxr+1),        &
                 imme(nzb+1:nzt,nxl-1:nxr+1), impb(nzb+1:nzt,nxl-1:nxr+1),        &
                 impe(nzb+1:nzt,nxl-1:nxr+1), ipmb(nzb+1:nzt,nxl-1:nxr+1),        &
                 ipme(nzb+1:nzt,nxl-1:nxr+1), ippb(nzb+1:nzt,nxl-1:nxr+1),        &
                 ippe(nzb+1:nzt,nxl-1:nxr+1), m1(nzb+1:nzt,nxl-2:nxr+2),          &
                 sw(nzb+1:nzt,nxl-1:nxr+1)                                        &
               )
       imme = 0.0_wp; impe = 0.0_wp; ipme = 0.0_wp; ippe = 0.0_wp

!
!--    Initialise point of time measuring of the exponential portion (this would
!--    not work if done locally within the loop)
       CALL cpu_log( log_point_s(12), 'advec_s_bc:exp', 'start' )
       CALL cpu_log( log_point_s(12), 'advec_s_bc:exp', 'pause' )

!
!--    Outer loop of all j
       DO  j = nys, nyn

!
!--       Compute polynomial coefficients
          DO  i = nxl-1, nxr+1
             DO  k = nzb+1, nzt
                a12(k,i) = 0.5_wp * ( sk_p(k,j,i+1) - sk_p(k,j,i-1) )
                a22(k,i) = 0.5_wp * ( sk_p(k,j,i+1) - 2.0_wp * sk_p(k,j,i)        &
                                                    + sk_p(k,j,i-1) )
                a0(k,i) = ( 9.0_wp * sk_p(k,j,i+2)    - 116.0_wp * sk_p(k,j,i+1)  &
                            + 2134.0_wp * sk_p(k,j,i) - 116.0_wp * sk_p(k,j,i-1)  &
                            + 9.0_wp * sk_p(k,j,i-2) ) * f1920
                a1(k,i) = ( -5.0_wp * sk_p(k,j,i+2)   + 34.0_wp * sk_p(k,j,i+1)   &
                            - 34.0_wp * sk_p(k,j,i-1) + 5.0_wp * sk_p(k,j,i-2)    &
                          ) * f48
                a2(k,i) = ( -3.0_wp * sk_p(k,j,i+2) + 36.0_wp * sk_p(k,j,i+1)     &
                            - 66.0_wp * sk_p(k,j,i) + 36.0_wp * sk_p(k,j,i-1)     &
                            - 3.0_wp * sk_p(k,j,i-2) ) * f48
             ENDDO
          ENDDO

!
!--       Fluxes using the Bott scheme
!--       *VOCL LOOP,UNROLL(2)
          DO  i = nxl, nxr
             DO  k = nzb+1, nzt
                cip  =  MAX( 0.0_wp, ( u(k,j,i+1) - u_gtrans ) * dt_3d * ddx )
                cim  = -MIN( 0.0_wp, ( u(k,j,i+1) - u_gtrans ) * dt_3d * ddx )
                cipf = 1.0_wp - 2.0_wp * cip
                cimf = 1.0_wp - 2.0_wp * cim
                ip   =   a0(k,i)   * f2  * ( 1.0_wp - cipf )                      &
                       + a1(k,i)   * f8  * ( 1.0_wp - cipf*cipf )                 &
                       + a2(k,i)   * f24 * ( 1.0_wp - cipf*cipf*cipf )
                im   =   a0(k,i+1) * f2  * ( 1.0_wp - cimf )                      &
                       - a1(k,i+1) * f8  * ( 1.0_wp - cimf*cimf )                 &
                       + a2(k,i+1) * f24 * ( 1.0_wp - cimf*cimf*cimf )
                ip   = MAX( ip, 0.0_wp )
                im   = MAX( im, 0.0_wp )
                ippb(k,i) = ip * MIN( 1.0_wp, sk_p(k,j,i)   / (ip+im+1E-15_wp) )
                impb(k,i) = im * MIN( 1.0_wp, sk_p(k,j,i+1) / (ip+im+1E-15_wp) )

                cip  =  MAX( 0.0_wp, ( u(k,j,i) - u_gtrans ) * dt_3d * ddx )
                cim  = -MIN( 0.0_wp, ( u(k,j,i) - u_gtrans ) * dt_3d * ddx )
                cipf = 1.0_wp - 2.0_wp * cip
                cimf = 1.0_wp - 2.0_wp * cim
                ip   =   a0(k,i-1) * f2  * ( 1.0_wp - cipf )                      &
                       + a1(k,i-1) * f8  * ( 1.0_wp - cipf*cipf )                 &
                       + a2(k,i-1) * f24 * ( 1.0_wp - cipf*cipf*cipf )
                im   =   a0(k,i)   * f2  * ( 1.0_wp - cimf )                      &
                       - a1(k,i)   * f8  * ( 1.0_wp - cimf*cimf )                 &
                       + a2(k,i)   * f24 * ( 1.0_wp - cimf*cimf*cimf )
                ip   = MAX( ip, 0.0_wp )
                im   = MAX( im, 0.0_wp )
                ipmb(k,i) = ip * MIN( 1.0_wp, sk_p(k,j,i-1) / (ip+im+1E-15_wp) )
                immb(k,i) = im * MIN( 1.0_wp, sk_p(k,j,i)   / (ip+im+1E-15_wp) )
             ENDDO
          ENDDO

!
!--       Compute monitor function m1
          DO  i = nxl-2, nxr+2
             DO  k = nzb+1, nzt
                m1z = ABS( sk_p(k,j,i+1) - 2.0_wp * sk_p(k,j,i) + sk_p(k,j,i-1) )
                m1n = ABS( sk_p(k,j,i+1) - sk_p(k,j,i-1) )
                IF ( m1n /= 0.0_wp  .AND.  m1n >= m1z )  THEN
                   m1(k,i) = m1z / m1n
                   IF ( m1(k,i) /= 2.0_wp  .AND.  m1n < fmax(2) )  m1(k,i) = 0.0_wp
                ELSEIF ( m1n < m1z )  THEN
                   m1(k,i) = -1.0_wp
                ELSE
                   m1(k,i) = 0.0_wp
                ENDIF
             ENDDO
          ENDDO

!
!--       Compute switch sw
          sw = 0.0_wp
          DO  i = nxl-1, nxr+1
             DO  k = nzb+1, nzt
                m2 = 2.0_wp * ABS( a1(k,i) - a12(k,i) ) /                         &
                     MAX( ABS( a1(k,i) + a12(k,i) ), 1E-35_wp )
                IF ( ABS( a1(k,i) + a12(k,i) ) < fmax(2) )  m2 = 0.0_wp

                m3 = 2.0_wp * ABS( a2(k,i) - a22(k,i) ) /                         &
                     MAX( ABS( a2(k,i) + a22(k,i) ), 1E-35_wp )
                IF ( ABS( a2(k,i) + a22(k,i) ) < fmax(1) )  m3 = 0.0_wp

                t1 = 0.35_wp
                t2 = 0.35_wp
                IF ( m1(k,i) == -1.0_wp )  t2 = 0.12_wp

!--             *VOCL STMT,IF(10)
                IF ( m1(k,i-1) == 1.0_wp .OR. m1(k,i) == 1.0_wp                   &
                     .OR. m1(k,i+1) == 1.0_wp .OR.  m2 > t2  .OR.  m3 > t2  .OR.  &
                     ( m1(k,i) > t1  .AND.  m1(k,i-1) /= -1.0_wp  .AND.           &
                       m1(k,i) /= -1.0_wp  .AND.  m1(k,i+1) /= -1.0_wp )          &
                   )  sw(k,i) = 1.0_wp
             ENDDO
          ENDDO

!
!--       Fluxes using the exponential scheme
          CALL cpu_log( log_point_s(12), 'advec_s_bc:exp', 'continue' )
          DO  i = nxl, nxr
             DO  k = nzb+1, nzt

!--             *VOCL STMT,IF(10)
                IF ( sw(k,i) == 1.0_wp )  THEN
                   snenn = sk_p(k,j,i+1) - sk_p(k,j,i-1)
                   IF ( ABS( snenn ) < 1E-9_wp )  snenn = 1E-9_wp
                   sterm = ( sk_p(k,j,i) - sk_p(k,j,i-1) ) / snenn
                   sterm = MIN( sterm, 0.9999_wp )
                   sterm = MAX( sterm, 0.0001_wp )

                   ix = INT( sterm * 1000 ) + 1

                   cip =  MAX( 0.0_wp, ( u(k,j,i+1) - u_gtrans ) * dt_3d * ddx )

                   ippe(k,i) = sk_p(k,j,i-1) * cip + snenn * (                    &
                               aex(ix) * cip + bex(ix) / dex(ix) * (              &
                               eex(ix) -                                          &
                               EXP( dex(ix)*0.5_wp * ( 1.0_wp - 2.0_wp * cip ) )  &
                                                                   )              &
                                                             )
                   IF ( sterm == 0.0001_wp )  ippe(k,i) = sk_p(k,j,i) * cip
                   IF ( sterm == 0.9999_wp )  ippe(k,i) = sk_p(k,j,i) * cip

                   snenn = sk_p(k,j,i-1) - sk_p(k,j,i+1)
                   IF ( ABS( snenn ) < 1E-9_wp )  snenn = 1E-9_wp
                   sterm = ( sk_p(k,j,i) - sk_p(k,j,i+1) ) / snenn
                   sterm = MIN( sterm, 0.9999_wp )
                   sterm = MAX( sterm, 0.0001_wp )

                   ix = INT( sterm * 1000 ) + 1

                   cim = -MIN( 0.0_wp, ( u(k,j,i) - u_gtrans ) * dt_3d * ddx )

                   imme(k,i) = sk_p(k,j,i+1) * cim + snenn * (                    &
                               aex(ix) * cim + bex(ix) / dex(ix) * (              &
                               eex(ix) -                                          &
                               EXP( dex(ix)*0.5_wp * ( 1.0_wp - 2.0_wp * cim ) )  &
                                                                   )              &
                                                             )
                   IF ( sterm == 0.0001_wp )  imme(k,i) = sk_p(k,j,i) * cim
                   IF ( sterm == 0.9999_wp )  imme(k,i) = sk_p(k,j,i) * cim
                ENDIF

!--             *VOCL STMT,IF(10)
                IF ( sw(k,i+1) == 1.0_wp )  THEN
                   snenn = sk_p(k,j,i) - sk_p(k,j,i+2)
                   IF ( ABS( snenn ) < 1E-9_wp )  snenn = 1E-9_wp
                   sterm = ( sk_p(k,j,i+1) - sk_p(k,j,i+2) ) / snenn
                   sterm = MIN( sterm, 0.9999_wp )
                   sterm = MAX( sterm, 0.0001_wp )

                   ix = INT( sterm * 1000 ) + 1

                   cim = -MIN( 0.0_wp, ( u(k,j,i+1) - u_gtrans ) * dt_3d * ddx )

                   impe(k,i) = sk_p(k,j,i+2) * cim + snenn * (                    &
                               aex(ix) * cim + bex(ix) / dex(ix) * (              &
                               eex(ix) -                                          &
                               EXP( dex(ix)*0.5_wp * ( 1.0_wp - 2.0_wp * cim ) )  &
                                                                   )              &
                                                             )
                   IF ( sterm == 0.0001_wp )  impe(k,i) = sk_p(k,j,i+1) * cim
                   IF ( sterm == 0.9999_wp )  impe(k,i) = sk_p(k,j,i+1) * cim
                ENDIF

!--             *VOCL STMT,IF(10)
                IF ( sw(k,i-1) == 1.0_wp )  THEN
                   snenn = sk_p(k,j,i) - sk_p(k,j,i-2)
                   IF ( ABS( snenn ) < 1E-9_wp )  snenn = 1E-9_wp
                   sterm = ( sk_p(k,j,i-1) - sk_p(k,j,i-2) ) / snenn
                   sterm = MIN( sterm, 0.9999_wp )
                   sterm = MAX( sterm, 0.0001_wp )

                   ix = INT( sterm * 1000 ) + 1

                   cip = MAX( 0.0_wp, ( u(k,j,i) - u_gtrans ) * dt_3d * ddx )

                   ipme(k,i) = sk_p(k,j,i-2) * cip + snenn * (                    &
                               aex(ix) * cip + bex(ix) / dex(ix) * (              &
                               eex(ix) -                                          &
                               EXP( dex(ix)*0.5_wp * ( 1.0_wp - 2.0_wp * cip ) )  &
                                                                   )              &
                                                             )
                   IF ( sterm == 0.0001_wp )  ipme(k,i) = sk_p(k,j,i-1) * cip
                   IF ( sterm == 0.9999_wp )  ipme(k,i) = sk_p(k,j,i-1) * cip
                ENDIF

             ENDDO
          ENDDO
          CALL cpu_log( log_point_s(12), 'advec_s_bc:exp', 'pause' )

!
!--       Prognostic equation
          DO  i = nxl, nxr
             DO  k = nzb+1, nzt
                fplus  = ( 1.0_wp - sw(k,i)   ) * ippb(k,i) + sw(k,i)   * ippe(k,i)  &
                       - ( 1.0_wp - sw(k,i+1) ) * impb(k,i) - sw(k,i+1) * impe(k,i)
                fminus = ( 1.0_wp - sw(k,i-1) ) * ipmb(k,i) + sw(k,i-1) * ipme(k,i)  &
                       - ( 1.0_wp - sw(k,i)   ) * immb(k,i) - sw(k,i)   * imme(k,i)
                tendcy = fplus - fminus
!
!--             Removed in order to optimize speed
!                ffmax   = MAX( ABS( fplus ), ABS( fminus ), 1E-35_wp )
!                IF ( ( ABS( tendcy ) / ffmax ) < 1E-7_wp )  tendcy = 0.0
!
!--             Density correction because of possible remaining divergences
                d_new = d(k,j,i) - ( u(k,j,i+1) - u(k,j,i) ) * dt_3d * ddx
                sk_p(k,j,i) = ( ( 1.0_wp + d(k,j,i) ) * sk_p(k,j,i) - tendcy ) /    &
                              ( 1.0_wp + d_new )
                d(k,j,i)  = d_new
             ENDDO
          ENDDO

       ENDDO   ! End of the advection in x-direction

!
!--    Deallocate temporary arrays
       DEALLOCATE( a0, a1, a2, a12, a22, immb, imme, impb, impe, ipmb, ipme,      &
                   ippb, ippe, m1, sw )


!
!--    Enlarge boundary of local array cyclically in y-direction
#if defined( __parallel )
       ngp = ( nzt - nzb + 6 ) * ( nyn - nys + 7 )
       CALL MPI_TYPE_VECTOR( nxr-nxl+7, 3*(nzt-nzb+6), ngp, MPI_REAL,             &
                             type_xz_2, ierr )
       CALL MPI_TYPE_COMMIT( type_xz_2, ierr )
!
!--    Send front boundary, receive rear boundary
       CALL cpu_log( log_point_s(11), 'advec_s_bc:sendrecv', 'continue' )
       CALL MPI_SENDRECV( sk_p(nzb-2,nys,nxl-3),   1, type_xz_2, psouth, 0,       &
                          sk_p(nzb-2,nyn+1,nxl-3), 1, type_xz_2, pnorth, 0,       &
                          comm2d, status, ierr )
!
!--    Send rear boundary, receive front boundary
       CALL MPI_SENDRECV( sk_p(nzb-2,nyn-2,nxl-3), 1, type_xz_2, pnorth, 1,       &
                          sk_p(nzb-2,nys-3,nxl-3), 1, type_xz_2, psouth, 1,       &
                          comm2d, status, ierr )
       CALL MPI_TYPE_FREE( type_xz_2, ierr )
       CALL cpu_log( log_point_s(11), 'advec_s_bc:sendrecv', 'pause' )
#else
       DO  i = nxl, nxr
          DO  k = nzb+1, nzt
             sk_p(k,nys-1,i) = sk_p(k,nyn,i)
             sk_p(k,nys-2,i) = sk_p(k,nyn-1,i)
             sk_p(k,nys-3,i) = sk_p(k,nyn-2,i)
             sk_p(k,nyn+1,i) = sk_p(k,nys,i)
             sk_p(k,nyn+2,i) = sk_p(k,nys+1,i)
             sk_p(k,nyn+3,i) = sk_p(k,nys+2,i)
          ENDDO
       ENDDO
#endif

!
!--    Determine the maxima of the first and second derivative in y-direction
       fmax_l = 0.0_wp
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                numera = ABS( sk_p(k,j+1,i) - 2.0_wp * sk_p(k,j,i) + sk_p(k,j-1,i) )
                denomi  = ABS( sk_p(k,j+1,i) - sk_p(k,j-1,i) )
                fmax_l(1) = MAX( fmax_l(1) , numera )
                fmax_l(2) = MAX( fmax_l(2) , denomi  )
             ENDDO
          ENDDO
       ENDDO
#if defined( __parallel )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( fmax_l, fmax, 2, MPI_REAL, MPI_MAX, comm2d, ierr )
#else
       fmax = fmax_l
#endif

       fmax = 0.04_wp * fmax

!
!--    Allocate temporary arrays
       ALLOCATE( a0(nzb+1:nzt,nys-1:nyn+1),   a1(nzb+1:nzt,nys-1:nyn+1),          &
                 a2(nzb+1:nzt,nys-1:nyn+1),   a12(nzb+1:nzt,nys-1:nyn+1),         &
                 a22(nzb+1:nzt,nys-1:nyn+1),  immb(nzb+1:nzt,nys-1:nyn+1),        &
                 imme(nzb+1:nzt,nys-1:nyn+1), impb(nzb+1:nzt,nys-1:nyn+1),        &
                 impe(nzb+1:nzt,nys-1:nyn+1), ipmb(nzb+1:nzt,nys-1:nyn+1),        &
                 ipme(nzb+1:nzt,nys-1:nyn+1), ippb(nzb+1:nzt,nys-1:nyn+1),        &
                 ippe(nzb+1:nzt,nys-1:nyn+1), m1(nzb+1:nzt,nys-2:nyn+2),          &
                 sw(nzb+1:nzt,nys-1:nyn+1)                                        &
               )
       imme = 0.0_wp; impe = 0.0_wp; ipme = 0.0_wp; ippe = 0.0_wp

!
!--    Outer loop of all i
       DO  i = nxl, nxr

!
!--       Compute polynomial coefficients
          DO  j = nys-1, nyn+1
             DO  k = nzb+1, nzt
                a12(k,j) = 0.5_wp * ( sk_p(k,j+1,i) - sk_p(k,j-1,i) )
                a22(k,j) = 0.5_wp * ( sk_p(k,j+1,i) - 2.0_wp * sk_p(k,j,i)        &
                                                    + sk_p(k,j-1,i) )
                a0(k,j) = ( 9.0_wp * sk_p(k,j+2,i)    - 116.0_wp * sk_p(k,j+1,i)  &
                            + 2134.0_wp * sk_p(k,j,i) - 116.0_wp * sk_p(k,j-1,i)  &
                            + 9.0_wp * sk_p(k,j-2,i) ) * f1920
                a1(k,j) = ( -5.0_wp   * sk_p(k,j+2,i) + 34.0_wp * sk_p(k,j+1,i)   &
                            - 34.0_wp * sk_p(k,j-1,i) + 5.0_wp  * sk_p(k,j-2,i)   &
                          ) * f48
                a2(k,j) = ( -3.0_wp * sk_p(k,j+2,i) + 36.0_wp * sk_p(k,j+1,i)     &
                            - 66.0_wp * sk_p(k,j,i) + 36.0_wp * sk_p(k,j-1,i)     &
                            - 3.0_wp * sk_p(k,j-2,i) ) * f48
             ENDDO
          ENDDO

!
!--       Fluxes using the Bott scheme
!--       *VOCL LOOP,UNROLL(2)
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                cip  =  MAX( 0.0_wp, ( v(k,j+1,i) - v_gtrans ) * dt_3d * ddy )
                cim  = -MIN( 0.0_wp, ( v(k,j+1,i) - v_gtrans ) * dt_3d * ddy )
                cipf = 1.0_wp - 2.0_wp * cip
                cimf = 1.0_wp - 2.0_wp * cim
                ip   =   a0(k,j)   * f2  * ( 1.0_wp - cipf )                      &
                       + a1(k,j)   * f8  * ( 1.0_wp - cipf*cipf )                 &
                       + a2(k,j)   * f24 * ( 1.0_wp - cipf*cipf*cipf )
                im   =   a0(k,j+1) * f2  * ( 1.0_wp - cimf )                      &
                       - a1(k,j+1) * f8  * ( 1.0_wp - cimf*cimf )                 &
                       + a2(k,j+1) * f24 * ( 1.0_wp - cimf*cimf*cimf )
                ip   = MAX( ip, 0.0_wp )
                im   = MAX( im, 0.0_wp )
                ippb(k,j) = ip * MIN( 1.0_wp, sk_p(k,j,i)   / (ip+im+1E-15_wp) )
                impb(k,j) = im * MIN( 1.0_wp, sk_p(k,j+1,i) / (ip+im+1E-15_wp) )

                cip  =  MAX( 0.0_wp, ( v(k,j,i) - v_gtrans ) * dt_3d * ddy )
                cim  = -MIN( 0.0_wp, ( v(k,j,i) - v_gtrans ) * dt_3d * ddy )
                cipf = 1.0_wp - 2.0_wp * cip
                cimf = 1.0_wp - 2.0_wp * cim
                ip   =   a0(k,j-1) * f2  * ( 1.0_wp - cipf )                      &
                       + a1(k,j-1) * f8  * ( 1.0_wp - cipf*cipf )                 &
                       + a2(k,j-1) * f24 * ( 1.0_wp - cipf*cipf*cipf )
                im   =   a0(k,j)   * f2  * ( 1.0_wp - cimf )                      &
                       - a1(k,j)   * f8  * ( 1.0_wp - cimf*cimf )                 &
                       + a2(k,j)   * f24 * ( 1.0_wp - cimf*cimf*cimf )
                ip   = MAX( ip, 0.0_wp )
                im   = MAX( im, 0.0_wp )
                ipmb(k,j) = ip * MIN( 1.0_wp, sk_p(k,j-1,i) / (ip+im+1E-15_wp) )
                immb(k,j) = im * MIN( 1.0_wp, sk_p(k,j,i)   / (ip+im+1E-15_wp) )
             ENDDO
          ENDDO

!
!--       Compute monitor function m1
          DO  j = nys-2, nyn+2
             DO  k = nzb+1, nzt
                m1z = ABS( sk_p(k,j+1,i) - 2.0_wp * sk_p(k,j,i) + sk_p(k,j-1,i) )
                m1n = ABS( sk_p(k,j+1,i) - sk_p(k,j-1,i) )
                IF ( m1n /= 0.0_wp  .AND.  m1n >= m1z )  THEN
                   m1(k,j) = m1z / m1n
                   IF ( m1(k,j) /= 2.0_wp  .AND.  m1n < fmax(2) )  m1(k,j) = 0.0_wp
                ELSEIF ( m1n < m1z )  THEN
                   m1(k,j) = -1.0_wp
                ELSE
                   m1(k,j) = 0.0_wp
                ENDIF
             ENDDO
          ENDDO

!
!--       Compute switch sw
          sw = 0.0_wp
          DO  j = nys-1, nyn+1
             DO  k = nzb+1, nzt
                m2 = 2.0_wp * ABS( a1(k,j) - a12(k,j) ) /                            &
                     MAX( ABS( a1(k,j) + a12(k,j) ), 1E-35_wp )
                IF ( ABS( a1(k,j) + a12(k,j) ) < fmax(2) )  m2 = 0.0_wp

                m3 = 2.0_wp * ABS( a2(k,j) - a22(k,j) ) /                            &
                     MAX( ABS( a2(k,j) + a22(k,j) ), 1E-35_wp )
                IF ( ABS( a2(k,j) + a22(k,j) ) < fmax(1) )  m3 = 0.0_wp

                t1 = 0.35_wp
                t2 = 0.35_wp
                IF ( m1(k,j) == -1.0_wp )  t2 = 0.12_wp

!--             *VOCL STMT,IF(10)
                IF ( m1(k,j-1) == 1.0_wp .OR. m1(k,j) == 1.0_wp                   &
                     .OR. m1(k,j+1) == 1.0_wp .OR.  m2 > t2  .OR.  m3 > t2  .OR.  &
                     ( m1(k,j) > t1  .AND.  m1(k,j-1) /= -1.0_wp  .AND.           &
                       m1(k,j) /= -1.0_wp  .AND.  m1(k,j+1) /= -1.0_wp )          &
                   )  sw(k,j) = 1.0_wp
             ENDDO
          ENDDO

!
!--       Fluxes using exponential scheme
          CALL cpu_log( log_point_s(12), 'advec_s_bc:exp', 'continue' )
          DO  j = nys, nyn
             DO  k = nzb+1, nzt

!--             *VOCL STMT,IF(10)
                IF ( sw(k,j) == 1.0_wp )  THEN
                   snenn = sk_p(k,j+1,i) - sk_p(k,j-1,i)
                   IF ( ABS( snenn ) < 1E-9_wp )  snenn = 1E-9_wp
                   sterm = ( sk_p(k,j,i) - sk_p(k,j-1,i) ) / snenn
                   sterm = MIN( sterm, 0.9999_wp )
                   sterm = MAX( sterm, 0.0001_wp )

                   ix = INT( sterm * 1000 ) + 1

                   cip =  MAX( 0.0_wp, ( v(k,j+1,i) - v_gtrans ) * dt_3d * ddy )

                   ippe(k,j) = sk_p(k,j-1,i) * cip + snenn * (                    &
                               aex(ix) * cip + bex(ix) / dex(ix) * (              &
                               eex(ix) -                                          &
                               EXP( dex(ix)*0.5_wp * ( 1.0_wp - 2.0_wp * cip ) )  &
                                                                   )              &
                                                             )
                   IF ( sterm == 0.0001_wp )  ippe(k,j) = sk_p(k,j,i) * cip
                   IF ( sterm == 0.9999_wp )  ippe(k,j) = sk_p(k,j,i) * cip

                   snenn = sk_p(k,j-1,i) - sk_p(k,j+1,i)
                   IF ( ABS( snenn ) < 1E-9_wp )  snenn = 1E-9_wp
                   sterm = ( sk_p(k,j,i) - sk_p(k,j+1,i) ) / snenn
                   sterm = MIN( sterm, 0.9999_wp )
                   sterm = MAX( sterm, 0.0001_wp )

                   ix = INT( sterm * 1000 ) + 1

                   cim = -MIN( 0.0_wp, ( v(k,j,i) - v_gtrans ) * dt_3d * ddy )

                   imme(k,j) = sk_p(k,j+1,i) * cim + snenn * (                    &
                               aex(ix) * cim + bex(ix) / dex(ix) * (              &
                               eex(ix) -                                          &
                               EXP( dex(ix)*0.5_wp * ( 1.0_wp - 2.0_wp * cim ) )  &
                                                                   )              &
                                                             )
                   IF ( sterm == 0.0001_wp )  imme(k,j) = sk_p(k,j,i) * cim
                   IF ( sterm == 0.9999_wp )  imme(k,j) = sk_p(k,j,i) * cim
                ENDIF

!--             *VOCL STMT,IF(10)
                IF ( sw(k,j+1) == 1.0_wp )  THEN
                   snenn = sk_p(k,j,i) - sk_p(k,j+2,i)
                   IF ( ABS( snenn ) < 1E-9_wp )  snenn = 1E-9_wp
                   sterm = ( sk_p(k,j+1,i) - sk_p(k,j+2,i) ) / snenn
                   sterm = MIN( sterm, 0.9999_wp )
                   sterm = MAX( sterm, 0.0001_wp )

                   ix = INT( sterm * 1000 ) + 1

                   cim = -MIN( 0.0_wp, ( v(k,j+1,i) - v_gtrans ) * dt_3d * ddy )

                   impe(k,j) = sk_p(k,j+2,i) * cim + snenn * (                    &
                               aex(ix) * cim + bex(ix) / dex(ix) * (              &
                               eex(ix) -                                          &
                               EXP( dex(ix)*0.5_wp * ( 1.0_wp - 2.0_wp * cim ) )  &
                                                                   )              &
                                                             )
                   IF ( sterm == 0.0001_wp )  impe(k,j) = sk_p(k,j+1,i) * cim
                   IF ( sterm == 0.9999_wp )  impe(k,j) = sk_p(k,j+1,i) * cim
                ENDIF

!--             *VOCL STMT,IF(10)
                IF ( sw(k,j-1) == 1.0_wp )  THEN
                   snenn = sk_p(k,j,i) - sk_p(k,j-2,i)
                   IF ( ABS( snenn ) < 1E-9_wp )  snenn = 1E-9_wp
                   sterm = ( sk_p(k,j-1,i) - sk_p(k,j-2,i) ) / snenn
                   sterm = MIN( sterm, 0.9999_wp )
                   sterm = MAX( sterm, 0.0001_wp )

                   ix = INT( sterm * 1000 ) + 1

                   cip = MAX( 0.0_wp, ( v(k,j,i) - v_gtrans ) * dt_3d * ddy )

                   ipme(k,j) = sk_p(k,j-2,i) * cip + snenn * (                    &
                               aex(ix) * cip + bex(ix) / dex(ix) * (              &
                               eex(ix) -                                          &
                               EXP( dex(ix)*0.5_wp * ( 1.0_wp - 2.0_wp * cip ) )  &
                                                                   )              &
                                                             )
                   IF ( sterm == 0.0001_wp )  ipme(k,j) = sk_p(k,j-1,i) * cip
                   IF ( sterm == 0.9999_wp )  ipme(k,j) = sk_p(k,j-1,i) * cip
                ENDIF

             ENDDO
          ENDDO
          CALL cpu_log( log_point_s(12), 'advec_s_bc:exp', 'pause' )

!
!--       Prognostic equation
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                fplus  = ( 1.0_wp - sw(k,j)   ) * ippb(k,j) + sw(k,j)   * ippe(k,j) &
                       - ( 1.0_wp - sw(k,j+1) ) * impb(k,j) - sw(k,j+1) * impe(k,j)
                fminus = ( 1.0_wp - sw(k,j-1) ) * ipmb(k,j) + sw(k,j-1) * ipme(k,j) &
                       - ( 1.0_wp - sw(k,j)   ) * immb(k,j) - sw(k,j)   * imme(k,j)
                tendcy = fplus - fminus
!
!--             Removed in order to optimise speed
!                ffmax   = MAX( ABS( fplus ), ABS( fminus ), 1E-35_wp )
!                 IF ( ( ABS( tendcy ) / ffmax ) < 1E-7_wp )  tendcy = 0.0
!
!--             Density correction because of possible remaining divergences
                d_new = d(k,j,i) - ( v(k,j+1,i) - v(k,j,i) ) * dt_3d * ddy
                sk_p(k,j,i) = ( ( 1.0_wp + d(k,j,i) ) * sk_p(k,j,i) - tendcy ) /  &
                              ( 1.0_wp + d_new )
                d(k,j,i)  = d_new
             ENDDO
          ENDDO

       ENDDO   ! End of the advection in y-direction
       CALL cpu_log( log_point_s(11), 'advec_s_bc:sendrecv', 'continue' )
       CALL cpu_log( log_point_s(11), 'advec_s_bc:sendrecv', 'stop' )

!
!--    Deallocate temporary arrays
       DEALLOCATE( a0, a1, a2, a12, a22, immb, imme, impb, impe, ipmb, ipme,      &
                   ippb, ippe, m1, sw )


!
!--    Initialise for the computation of heat fluxes (see below; required in
!--    UP flow_statistics)
       IF ( sk_char == 'pt' )  sums_wsts_bc_l = 0.0_wp

!
!--    Add top and bottom boundaries according to the relevant boundary conditions
       IF ( sk_char == 'pt' )  THEN

!
!--       Temperature boundary condition at the bottom boundary
          IF ( ibc_pt_b == 0 )  THEN
!
!--       Dirichlet (fixed surface temperature)
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   sk_p(nzb-1,j,i) = sk_p(nzb,j,i)
                   sk_p(nzb-2,j,i) = sk_p(nzb,j,i)
                ENDDO
             ENDDO

          ELSE
!
!--          Neumann (i.e. here zero gradient)
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   sk_p(nzb,j,i)   = sk_p(nzb+1,j,i)
                   sk_p(nzb-1,j,i) = sk_p(nzb,j,i)
                   sk_p(nzb-2,j,i) = sk_p(nzb,j,i)
                ENDDO
             ENDDO

          ENDIF

!
!--       Temperature boundary condition at the top boundary
          IF ( ibc_pt_t == 0  .OR.  ibc_pt_t == 1 )  THEN
!
!--          Dirichlet or Neumann (zero gradient)
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   sk_p(nzt+2,j,i)   = sk_p(nzt+1,j,i)
                   sk_p(nzt+3,j,i)   = sk_p(nzt+1,j,i)
                ENDDO
             ENDDO

          ELSEIF ( ibc_pt_t == 2 )  THEN
!
!--          Neumann: dzu(nzt+2:3) are not defined, dzu(nzt+1) is used instead
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   sk_p(nzt+2,j,i)   = sk_p(nzt+1,j,i) + bc_pt_t_val * dzu(nzt+1)
                   sk_p(nzt+3,j,i)   = sk_p(nzt+2,j,i) + bc_pt_t_val * dzu(nzt+1)
                ENDDO
             ENDDO

          ENDIF

       ELSEIF ( sk_char == 'sa' )  THEN

!
!--       Salinity boundary condition at the bottom boundary.
!--       So far, always Neumann (i.e. here zero gradient) is used
          DO  i = nxl, nxr
             DO  j = nys, nyn
                sk_p(nzb,j,i)   = sk_p(nzb+1,j,i)
                sk_p(nzb-1,j,i) = sk_p(nzb,j,i)
                sk_p(nzb-2,j,i) = sk_p(nzb,j,i)
             ENDDO
          ENDDO

!
!--       Salinity boundary condition at the top boundary.
!--       Dirichlet or Neumann (zero gradient)
          DO  i = nxl, nxr
             DO  j = nys, nyn
                sk_p(nzt+2,j,i)   = sk_p(nzt+1,j,i)
                sk_p(nzt+3,j,i)   = sk_p(nzt+1,j,i)
             ENDDO
          ENDDO

       ELSEIF ( sk_char == 'q' )  THEN

!
!--       Specific humidity boundary condition at the bottom boundary.
!--       Dirichlet (fixed surface humidity) or Neumann (i.e. zero gradient)
          DO  i = nxl, nxr
             DO  j = nys, nyn
                sk_p(nzb,j,i)   = sk_p(nzb+1,j,i)
                sk_p(nzb-1,j,i) = sk_p(nzb,j,i)
                sk_p(nzb-2,j,i) = sk_p(nzb,j,i)
             ENDDO
          ENDDO

!
!--       Specific humidity boundary condition at the top boundary
          IF ( ibc_q_t == 0 )  THEN
!
!--          Dirichlet
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   sk_p(nzt+2,j,i)   = sk_p(nzt+1,j,i)
                   sk_p(nzt+3,j,i)   = sk_p(nzt+1,j,i)
                ENDDO
             ENDDO

          ELSE
!
!--          Neumann: dzu(nzt+2:3) are not defined, dzu(nzt+1) is used instead
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   sk_p(nzt+2,j,i)   = sk_p(nzt+1,j,i) + bc_q_t_val * dzu(nzt+1)
                   sk_p(nzt+3,j,i)   = sk_p(nzt+2,j,i) + bc_q_t_val * dzu(nzt+1)
                ENDDO
             ENDDO

          ENDIF

       ELSEIF ( sk_char == 's' )  THEN
!
!--       Specific scalar boundary condition at the bottom boundary.
!--       Dirichlet (fixed surface humidity) or Neumann (i.e. zero gradient)
          DO  i = nxl, nxr
             DO  j = nys, nyn
                sk_p(nzb,j,i)   = sk_p(nzb+1,j,i)
                sk_p(nzb-1,j,i) = sk_p(nzb,j,i)
                sk_p(nzb-2,j,i) = sk_p(nzb,j,i)
             ENDDO
          ENDDO

!
!--       Specific scalar boundary condition at the top boundary
          IF ( ibc_s_t == 0 )  THEN
!
!--          Dirichlet
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   sk_p(nzt+2,j,i)   = sk_p(nzt+1,j,i)
                   sk_p(nzt+3,j,i)   = sk_p(nzt+1,j,i)
                ENDDO
             ENDDO

          ELSE
!
!--          Neumann: dzu(nzt+2:3) are not defined, dzu(nzt+1) is used instead
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   sk_p(nzt+2,j,i)   = sk_p(nzt+1,j,i) + bc_s_t_val * dzu(nzt+1)
                   sk_p(nzt+3,j,i)   = sk_p(nzt+2,j,i) + bc_s_t_val * dzu(nzt+1)
                ENDDO
             ENDDO

          ENDIF

       ELSEIF ( sk_char == 'qc' )  THEN

!
!--       Cloud water content boundary condition at the bottom boundary:
!--       Dirichlet (fixed surface rain water content).
          DO  i = nxl, nxr
             DO  j = nys, nyn
                sk_p(nzb,j,i)   = sk_p(nzb+1,j,i)
                sk_p(nzb-1,j,i) = sk_p(nzb,j,i)
                sk_p(nzb-2,j,i) = sk_p(nzb,j,i)
             ENDDO
          ENDDO

!
!--       Cloud water content boundary condition at the top boundary: Dirichlet
          DO  i = nxl, nxr
             DO  j = nys, nyn
                sk_p(nzt+2,j,i)   = sk_p(nzt+1,j,i)
                sk_p(nzt+3,j,i)   = sk_p(nzt+1,j,i)
             ENDDO
          ENDDO

       ELSEIF ( sk_char == 'qr' )  THEN

!
!--       Rain water content boundary condition at the bottom boundary:
!--       Dirichlet (fixed surface rain water content).
          DO  i = nxl, nxr
             DO  j = nys, nyn
                sk_p(nzb,j,i)   = sk_p(nzb+1,j,i)
                sk_p(nzb-1,j,i) = sk_p(nzb,j,i)
                sk_p(nzb-2,j,i) = sk_p(nzb,j,i)
             ENDDO
          ENDDO

!
!--       Rain water content boundary condition at the top boundary: Dirichlet
          DO  i = nxl, nxr
             DO  j = nys, nyn
                sk_p(nzt+2,j,i)   = sk_p(nzt+1,j,i)
                sk_p(nzt+3,j,i)   = sk_p(nzt+1,j,i)
             ENDDO
          ENDDO

       ELSEIF ( sk_char == 'nc' )  THEN

!
!--       Cloud drop concentration boundary condition at the bottom boundary:
!--       Dirichlet (fixed surface cloud drop concentration).
          DO  i = nxl, nxr
             DO  j = nys, nyn
                sk_p(nzb,j,i)   = sk_p(nzb+1,j,i)
                sk_p(nzb-1,j,i) = sk_p(nzb,j,i)
                sk_p(nzb-2,j,i) = sk_p(nzb,j,i)
             ENDDO
          ENDDO

!
!--       Cloud drop concentration boundary condition at the top boundary: Dirichlet
          DO  i = nxl, nxr
             DO  j = nys, nyn
                sk_p(nzt+2,j,i)   = sk_p(nzt+1,j,i)
                sk_p(nzt+3,j,i)   = sk_p(nzt+1,j,i)
             ENDDO
          ENDDO

       ELSEIF ( sk_char == 'nr' )  THEN

!
!--       Rain drop concentration boundary condition at the bottom boundary:
!--       Dirichlet (fixed surface rain drop concentration).
          DO  i = nxl, nxr
             DO  j = nys, nyn
                sk_p(nzb,j,i)   = sk_p(nzb+1,j,i)
                sk_p(nzb-1,j,i) = sk_p(nzb,j,i)
                sk_p(nzb-2,j,i) = sk_p(nzb,j,i)
             ENDDO
          ENDDO

!
!--       Rain drop concentration boundary condition at the top boundary: Dirichlet
          DO  i = nxl, nxr
             DO  j = nys, nyn
                sk_p(nzt+2,j,i)   = sk_p(nzt+1,j,i)
                sk_p(nzt+3,j,i)   = sk_p(nzt+1,j,i)
             ENDDO
          ENDDO

       ELSEIF ( sk_char == 'e' )  THEN

!
!--       TKE boundary condition at bottom and top boundary (generally Neumann)
          DO  i = nxl, nxr
             DO  j = nys, nyn
                sk_p(nzb,j,i)   = sk_p(nzb+1,j,i)
                sk_p(nzb-1,j,i) = sk_p(nzb,j,i)
                sk_p(nzb-2,j,i) = sk_p(nzb,j,i)
                sk_p(nzt+2,j,i) = sk_p(nzt+1,j,i)
                sk_p(nzt+3,j,i) = sk_p(nzt+1,j,i)
             ENDDO
          ENDDO

       ELSE

          WRITE( message_string, * ) 'no vertical boundary condi',                &
                                   'tion for variable "', sk_char, '"'
          CALL message( 'advec_s_bc', 'PA0158', 1, 2, 0, 6, 0 )
         
       ENDIF

!
!--    Determine the maxima of the first and second derivative in z-direction
       fmax_l = 0.0_wp
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb, nzt+1
                numera = ABS( sk_p(k+1,j,i) - 2.0_wp * sk_p(k,j,i) + sk_p(k-1,j,i) )
                denomi  = ABS( sk_p(k+1,j,i+1) - sk_p(k-1,j,i) )
                fmax_l(1) = MAX( fmax_l(1) , numera )
                fmax_l(2) = MAX( fmax_l(2) , denomi  )
             ENDDO
          ENDDO
       ENDDO
#if defined( __parallel )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( fmax_l, fmax, 2, MPI_REAL, MPI_MAX, comm2d, ierr )
#else
       fmax = fmax_l
#endif

       fmax = 0.04_wp * fmax

!
!--    Allocate temporary arrays
       ALLOCATE( a0(nzb:nzt+1,nys:nyn),   a1(nzb:nzt+1,nys:nyn),                  &
                 a2(nzb:nzt+1,nys:nyn),   a12(nzb:nzt+1,nys:nyn),                 &
                 a22(nzb:nzt+1,nys:nyn),  immb(nzb+1:nzt,nys:nyn),                &
                 imme(nzb+1:nzt,nys:nyn), impb(nzb+1:nzt,nys:nyn),                &
                 impe(nzb+1:nzt,nys:nyn), ipmb(nzb+1:nzt,nys:nyn),                &
                 ipme(nzb+1:nzt,nys:nyn), ippb(nzb+1:nzt,nys:nyn),                &
                 ippe(nzb+1:nzt,nys:nyn), m1(nzb-1:nzt+2,nys:nyn),                &
                 sw(nzb:nzt+1,nys:nyn)                                            &
               )
       imme = 0.0_wp; impe = 0.0_wp; ipme = 0.0_wp; ippe = 0.0_wp

!
!--    Outer loop of all i
       DO  i = nxl, nxr

!
!--       Compute polynomial coefficients
          DO  j = nys, nyn
             DO  k = nzb, nzt+1
                a12(k,j) = 0.5_wp * ( sk_p(k+1,j,i) - sk_p(k-1,j,i) )
                a22(k,j) = 0.5_wp * ( sk_p(k+1,j,i) - 2.0_wp * sk_p(k,j,i)        &
                                                    + sk_p(k-1,j,i) )
                a0(k,j) = ( 9.0_wp * sk_p(k+2,j,i)    - 116.0_wp * sk_p(k+1,j,i)  &
                            + 2134.0_wp * sk_p(k,j,i) - 116.0_wp * sk_p(k-1,j,i)  &
                            + 9.0_wp * sk_p(k-2,j,i) ) * f1920
                a1(k,j) = ( -5.0_wp   * sk_p(k+2,j,i) + 34.0_wp * sk_p(k+1,j,i)   &
                            - 34.0_wp * sk_p(k-1,j,i) + 5.0_wp  * sk_p(k-2,j,i)   &
                          ) * f48
                a2(k,j) = ( -3.0_wp * sk_p(k+2,j,i) + 36.0_wp * sk_p(k+1,j,i)     &
                            - 66.0_wp * sk_p(k,j,i) + 36.0_wp * sk_p(k-1,j,i)     &
                            - 3.0_wp * sk_p(k-2,j,i) ) * f48
             ENDDO
          ENDDO

!
!--       Fluxes using the Bott scheme
!--       *VOCL LOOP,UNROLL(2)
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                cip  =  MAX( 0.0_wp, w(k,j,i) * dt_3d * ddzw(k) )
                cim  = -MIN( 0.0_wp, w(k,j,i) * dt_3d * ddzw(k) )
                cipf = 1.0_wp - 2.0_wp * cip
                cimf = 1.0_wp - 2.0_wp * cim
                ip   =   a0(k,j)   * f2  * ( 1.0_wp - cipf )                      &
                       + a1(k,j)   * f8  * ( 1.0_wp - cipf*cipf )                 &
                       + a2(k,j)   * f24 * ( 1.0_wp - cipf*cipf*cipf )
                im   =   a0(k+1,j) * f2  * ( 1.0_wp - cimf )                      &
                       - a1(k+1,j) * f8  * ( 1.0_wp - cimf*cimf )                 &
                       + a2(k+1,j) * f24 * ( 1.0_wp - cimf*cimf*cimf )
                ip   = MAX( ip, 0.0_wp )
                im   = MAX( im, 0.0_wp )
                ippb(k,j) = ip * MIN( 1.0_wp, sk_p(k,j,i)   / (ip+im+1E-15_wp) )
                impb(k,j) = im * MIN( 1.0_wp, sk_p(k+1,j,i) / (ip+im+1E-15_wp) )

                cip  =  MAX( 0.0_wp, w(k-1,j,i) * dt_3d * ddzw(k) )
                cim  = -MIN( 0.0_wp, w(k-1,j,i) * dt_3d * ddzw(k) )
                cipf = 1.0_wp - 2.0_wp * cip
                cimf = 1.0_wp - 2.0_wp * cim
                ip   =   a0(k-1,j) * f2  * ( 1.0_wp - cipf )                      &
                       + a1(k-1,j) * f8  * ( 1.0_wp - cipf*cipf )                 &
                       + a2(k-1,j) * f24 * ( 1.0_wp - cipf*cipf*cipf )
                im   =   a0(k,j)   * f2  * ( 1.0_wp - cimf )                      &
                       - a1(k,j)   * f8  * ( 1.0_wp - cimf*cimf )                 &
                       + a2(k,j)   * f24 * ( 1.0_wp - cimf*cimf*cimf )
                ip   = MAX( ip, 0.0_wp )
                im   = MAX( im, 0.0_wp )
                ipmb(k,j) = ip * MIN( 1.0_wp, sk_p(k-1,j,i) / (ip+im+1E-15_wp) )
                immb(k,j) = im * MIN( 1.0_wp, sk_p(k,j,i)   / (ip+im+1E-15_wp) )
             ENDDO
          ENDDO

!
!--       Compute monitor function m1
          DO  j = nys, nyn
             DO  k = nzb-1, nzt+2
                m1z = ABS( sk_p(k+1,j,i) - 2.0_wp * sk_p(k,j,i) + sk_p(k-1,j,i) )
                m1n = ABS( sk_p(k+1,j,i) - sk_p(k-1,j,i) )
                IF ( m1n /= 0.0_wp  .AND.  m1n >= m1z )  THEN
                   m1(k,j) = m1z / m1n
                   IF ( m1(k,j) /= 2.0_wp  .AND.  m1n < fmax(2) )  m1(k,j) = 0.0_wp
                ELSEIF ( m1n < m1z )  THEN
                   m1(k,j) = -1.0_wp
                ELSE
                   m1(k,j) = 0.0_wp
                ENDIF
             ENDDO
          ENDDO

!
!--       Compute switch sw
          sw = 0.0_wp
          DO  j = nys, nyn
             DO  k = nzb, nzt+1
                m2 = 2.0_wp * ABS( a1(k,j) - a12(k,j) ) /                         &
                     MAX( ABS( a1(k,j) + a12(k,j) ), 1E-35_wp )
                IF ( ABS( a1(k,j) + a12(k,j) ) < fmax(2) )  m2 = 0.0_wp

                m3 = 2.0_wp * ABS( a2(k,j) - a22(k,j) ) /                         &
                     MAX( ABS( a2(k,j) + a22(k,j) ), 1E-35_wp )
                IF ( ABS( a2(k,j) + a22(k,j) ) < fmax(1) )  m3 = 0.0_wp

                t1 = 0.35_wp
                t2 = 0.35_wp
                IF ( m1(k,j) == -1.0_wp )  t2 = 0.12_wp

!--             *VOCL STMT,IF(10)
                IF ( m1(k-1,j) == 1.0_wp .OR. m1(k,j) == 1.0_wp                   &
                     .OR. m1(k+1,j) == 1.0_wp .OR.  m2 > t2  .OR.  m3 > t2  .OR.  &
                     ( m1(k,j) > t1  .AND.  m1(k-1,j) /= -1.0_wp  .AND.           &
                       m1(k,j) /= -1.0_wp  .AND.  m1(k+1,j) /= -1.0_wp )          &
                   )  sw(k,j) = 1.0_wp
             ENDDO
          ENDDO

!
!--       Fluxes using exponential scheme
          CALL cpu_log( log_point_s(12), 'advec_s_bc:exp', 'continue' )
          DO  j = nys, nyn
             DO  k = nzb+1, nzt

!--             *VOCL STMT,IF(10)
                IF ( sw(k,j) == 1.0_wp )  THEN
                   snenn = sk_p(k+1,j,i) - sk_p(k-1,j,i)
                   IF ( ABS( snenn ) < 1E-9_wp )  snenn = 1E-9_wp
                   sterm = ( sk_p(k,j,i) - sk_p(k-1,j,i) ) / snenn
                   sterm = MIN( sterm, 0.9999_wp )
                   sterm = MAX( sterm, 0.0001_wp )

                   ix = INT( sterm * 1000 ) + 1

                   cip =  MAX( 0.0_wp, w(k,j,i) * dt_3d * ddzw(k) )

                   ippe(k,j) = sk_p(k-1,j,i) * cip + snenn * (                    &
                               aex(ix) * cip + bex(ix) / dex(ix) * (              &
                               eex(ix) -                                          &
                               EXP( dex(ix)*0.5_wp * ( 1.0_wp - 2.0_wp * cip ) )  &
                                                                   )              &
                                                             )
                   IF ( sterm == 0.0001_wp )  ippe(k,j) = sk_p(k,j,i) * cip
                   IF ( sterm == 0.9999_wp )  ippe(k,j) = sk_p(k,j,i) * cip

                   snenn = sk_p(k-1,j,i) - sk_p(k+1,j,i)
                   IF ( ABS( snenn ) < 1E-9_wp )  snenn = 1E-9_wp
                   sterm = ( sk_p(k,j,i) - sk_p(k+1,j,i) ) / snenn
                   sterm = MIN( sterm, 0.9999_wp )
                   sterm = MAX( sterm, 0.0001_wp )

                   ix = INT( sterm * 1000 ) + 1

                   cim = -MIN( 0.0_wp, w(k-1,j,i) * dt_3d * ddzw(k) )

                   imme(k,j) = sk_p(k+1,j,i) * cim + snenn * (                    &
                               aex(ix) * cim + bex(ix) / dex(ix) * (              &
                               eex(ix) -                                          &
                               EXP( dex(ix)*0.5_wp * ( 1.0_wp - 2.0_wp * cim ) ) &
                                                                   )              &
                                                             )
                   IF ( sterm == 0.0001_wp )  imme(k,j) = sk_p(k,j,i) * cim
                   IF ( sterm == 0.9999_wp )  imme(k,j) = sk_p(k,j,i) * cim
                ENDIF

!--             *VOCL STMT,IF(10)
                IF ( sw(k+1,j) == 1.0_wp )  THEN
                   snenn = sk_p(k,j,i) - sk_p(k+2,j,i)
                   IF ( ABS( snenn ) < 1E-9_wp )  snenn = 1E-9_wp
                   sterm = ( sk_p(k+1,j,i) - sk_p(k+2,j,i) ) / snenn
                   sterm = MIN( sterm, 0.9999_wp )
                   sterm = MAX( sterm, 0.0001_wp )

                   ix = INT( sterm * 1000 ) + 1

                   cim = -MIN( 0.0_wp, w(k,j,i) * dt_3d * ddzw(k) )

                   impe(k,j) = sk_p(k+2,j,i) * cim + snenn * (                    &
                               aex(ix) * cim + bex(ix) / dex(ix) * (              &
                               eex(ix) -                                          &
                               EXP( dex(ix)*0.5_wp * ( 1.0_wp - 2.0_wp * cim ) )  &
                                                                   )              &
                                                             )
                   IF ( sterm == 0.0001_wp )  impe(k,j) = sk_p(k+1,j,i) * cim
                   IF ( sterm == 0.9999_wp )  impe(k,j) = sk_p(k+1,j,i) * cim
                ENDIF

!--             *VOCL STMT,IF(10)
                IF ( sw(k-1,j) == 1.0_wp )  THEN
                   snenn = sk_p(k,j,i) - sk_p(k-2,j,i)
                   IF ( ABS( snenn ) < 1E-9_wp )  snenn = 1E-9_wp
                   sterm = ( sk_p(k-1,j,i) - sk_p(k-2,j,i) ) / snenn
                   sterm = MIN( sterm, 0.9999_wp )
                   sterm = MAX( sterm, 0.0001_wp )

                   ix = INT( sterm * 1000 ) + 1

                   cip = MAX( 0.0_wp, w(k-1,j,i) * dt_3d * ddzw(k) )

                   ipme(k,j) = sk_p(k-2,j,i) * cip + snenn * (                    &
                               aex(ix) * cip + bex(ix) / dex(ix) * (              &
                               eex(ix) -                                          &
                               EXP( dex(ix)*0.5_wp * ( 1.0_wp - 2.0_wp * cip ) )  &
                                                                   )              &
                                                             )
                   IF ( sterm == 0.0001_wp )  ipme(k,j) = sk_p(k-1,j,i) * cip
                   IF ( sterm == 0.9999_wp )  ipme(k,j) = sk_p(k-1,j,i) * cip
                ENDIF

             ENDDO
          ENDDO
          CALL cpu_log( log_point_s(12), 'advec_s_bc:exp', 'pause' )

!
!--       Prognostic equation
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                fplus  = ( 1.0_wp - sw(k,j)   ) * ippb(k,j) + sw(k,j)   * ippe(k,j) &
                       - ( 1.0_wp - sw(k+1,j) ) * impb(k,j) - sw(k+1,j) * impe(k,j)
                fminus = ( 1.0_wp - sw(k-1,j) ) * ipmb(k,j) + sw(k-1,j) * ipme(k,j) &
                       - ( 1.0_wp - sw(k,j)   ) * immb(k,j) - sw(k,j)   * imme(k,j)
                tendcy = fplus - fminus
!
!--              Removed in order to optimise speed
!                ffmax   = MAX( ABS( fplus ), ABS( fminus ), 1E-35_wp )
!                IF ( ( ABS( tendcy ) / ffmax ) < 1E-7_wp )  tendcy = 0.0
!
!--             Density correction because of possible remaining divergences
                d_new = d(k,j,i) - ( w(k,j,i) - w(k-1,j,i) ) * dt_3d * ddzw(k)
                sk_p(k,j,i) = ( ( 1.0_wp + d(k,j,i) ) * sk_p(k,j,i) - tendcy ) /  &
                              ( 1.0_wp + d_new )
!
!--             Store heat flux for subsequent statistics output.
!--             array m1 is here used as temporary storage
                m1(k,j) = fplus / dt_3d * dzw(k)
             ENDDO
          ENDDO

!
!--       Sum up heat flux in order to order to obtain horizontal averages
          IF ( sk_char == 'pt' )  THEN
             DO  sr = 0, statistic_regions
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt
                      sums_wsts_bc_l(k,sr) = sums_wsts_bc_l(k,sr) +               &
                                             m1(k,j) * rmask(j,i,sr)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       ENDDO   ! End of the advection in z-direction
       CALL cpu_log( log_point_s(12), 'advec_s_bc:exp', 'continue' )
       CALL cpu_log( log_point_s(12), 'advec_s_bc:exp', 'stop' )

!
!--    Deallocate temporary arrays
       DEALLOCATE( a0, a1, a2, a12, a22, immb, imme, impb, impe, ipmb, ipme,      &
                   ippb, ippe, m1, sw )

!
!--    Store results as tendency and deallocate local array
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                tend(k,j,i) = tend(k,j,i) + ( sk_p(k,j,i) - sk(k,j,i) ) / dt_3d
             ENDDO
          ENDDO
       ENDDO

       DEALLOCATE( sk_p )

    END SUBROUTINE advec_s_bc

 END MODULE advec_s_bc_mod
