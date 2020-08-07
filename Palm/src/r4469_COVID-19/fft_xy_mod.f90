!> @file fft_xy_mod.f90
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
! $Id: fft_xy_mod.f90 4370 2020-01-10 14:00:44Z raasch $
! bugfix for Temperton-fft usage on GPU
! 
! 4366 2020-01-09 08:12:43Z raasch
! Vectorized Temperton-fft added
! 
! 4360 2020-01-07 11:25:50Z suehring
! Corrected "Former revisions" section
! 
! 4069 2019-07-01 14:05:51Z Giersch
! Code added to avoid compiler warnings
! 
! 3655 2019-01-07 16:51:22Z knoop
! OpenACC port for SPEC
!
! Revision 1.1  2002/06/11 13:00:49  raasch
! Initial revision
!
!
! Description:
! ------------
!> Fast Fourier transformation along x and y for 1d domain decomposition along x.
!> Original version: Klaus Ketelsen (May 2002)
!> @todo openmp support for vectorized Temperton fft
!------------------------------------------------------------------------------!
 MODULE fft_xy
 

    USE control_parameters,                                                    &
        ONLY:  fft_method, loop_optimization, message_string
        
    USE cuda_fft_interfaces
        
    USE indices,                                                               &
        ONLY:  nx, ny, nz

#if defined( __cuda_fft )
    USE ISO_C_BINDING
#elif defined( __fftw )
    USE, INTRINSIC ::  ISO_C_BINDING
#endif

    USE kinds
    
    USE singleton,                                                             &
        ONLY: fftn
    
    USE temperton_fft
    
    USE transpose_indices,                                                     &
        ONLY:  nxl_y, nxr_y, nyn_x, nys_x, nzb_x, nzb_y, nzt_x, nzt_y

    IMPLICIT NONE

    PRIVATE
    PUBLIC fft_x, fft_x_1d, fft_y, fft_y_1d, fft_init, fft_x_m, fft_y_m, f_vec_x, temperton_fft_vec

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE, SAVE ::  ifax_x  !<
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE, SAVE ::  ifax_y  !<

    LOGICAL, SAVE ::  init_fft = .FALSE.           !<
    LOGICAL, SAVE ::  temperton_fft_vec = .FALSE.  !<

    REAL(wp), SAVE ::  dnx      !<
    REAL(wp), SAVE ::  dny      !<
    REAL(wp), SAVE ::  sqr_dnx  !<
    REAL(wp), SAVE ::  sqr_dny  !<
    
    REAL(wp), DIMENSION(:), ALLOCATABLE, SAVE ::  trigs_x  !< 
    REAL(wp), DIMENSION(:), ALLOCATABLE, SAVE ::  trigs_y  !<

    REAL(wp), DIMENSION(:,:), ALLOCATABLE, SAVE ::  f_vec_x

#if defined( __ibm )
    INTEGER(iwp), PARAMETER ::  nau1 = 20000  !< 
    INTEGER(iwp), PARAMETER ::  nau2 = 22000  !<
!
!-- The following working arrays contain tables and have to be "save" and
!-- shared in OpenMP sense
    REAL(wp), DIMENSION(nau1), SAVE ::  aux1  !< 
    REAL(wp), DIMENSION(nau1), SAVE ::  auy1  !<
    REAL(wp), DIMENSION(nau1), SAVE ::  aux3  !< 
    REAL(wp), DIMENSION(nau1), SAVE ::  auy3  !<
    
#elif defined( __nec_fft )
    INTEGER(iwp), SAVE ::  nz1  !<
    
    REAL(wp), DIMENSION(:), ALLOCATABLE, SAVE ::  trig_xb  !<
    REAL(wp), DIMENSION(:), ALLOCATABLE, SAVE ::  trig_xf  !< 
    REAL(wp), DIMENSION(:), ALLOCATABLE, SAVE ::  trig_yb  !<
    REAL(wp), DIMENSION(:), ALLOCATABLE, SAVE ::  trig_yf  !<
    
#elif defined( __cuda_fft )
    INTEGER(C_INT), SAVE ::  plan_xf  !<
    INTEGER(C_INT), SAVE ::  plan_xi  !<
    INTEGER(C_INT), SAVE ::  plan_yf  !<
    INTEGER(C_INT), SAVE ::  plan_yi  !<

#endif

#if defined( __fftw )
    INCLUDE  'fftw3.f03'
    INTEGER(KIND=C_INT) ::  nx_c  !<
    INTEGER(KIND=C_INT) ::  ny_c  !<
    
    COMPLEX(KIND=C_DOUBLE_COMPLEX), DIMENSION(:), ALLOCATABLE, SAVE ::  x_out  !<
    COMPLEX(KIND=C_DOUBLE_COMPLEX), DIMENSION(:), ALLOCATABLE, SAVE ::         &
       y_out  !<
    
    REAL(KIND=C_DOUBLE), DIMENSION(:), ALLOCATABLE, SAVE ::                    &
       x_in   !<
    REAL(KIND=C_DOUBLE), DIMENSION(:), ALLOCATABLE, SAVE ::                    &
       y_in   !<
    !$OMP THREADPRIVATE( x_out, y_out, x_in, y_in )
    
    
    TYPE(C_PTR), SAVE ::  plan_xf, plan_xi, plan_yf, plan_yi
#endif

!
!-- Public interfaces
    INTERFACE fft_init
       MODULE PROCEDURE fft_init
    END INTERFACE fft_init

    INTERFACE fft_x
       MODULE PROCEDURE fft_x
    END INTERFACE fft_x

    INTERFACE fft_x_1d
       MODULE PROCEDURE fft_x_1d
    END INTERFACE fft_x_1d

    INTERFACE fft_y
       MODULE PROCEDURE fft_y
    END INTERFACE fft_y

    INTERFACE fft_y_1d
       MODULE PROCEDURE fft_y_1d
    END INTERFACE fft_y_1d

    INTERFACE fft_x_m
       MODULE PROCEDURE fft_x_m
    END INTERFACE fft_x_m

    INTERFACE fft_y_m
       MODULE PROCEDURE fft_y_m
    END INTERFACE fft_y_m

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE fft_init

       USE pegrid,                                                                                 &
           ONLY:  pdims

       IMPLICIT NONE

!
!--    The following temporary working arrays have to be on stack or private
!--    in OpenMP sense
#if defined( __ibm )
       REAL(wp), DIMENSION(0:nx+2) ::  workx  !<
       REAL(wp), DIMENSION(0:ny+2) ::  worky  !<
       REAL(wp), DIMENSION(nau2)   ::  aux2   !< 
       REAL(wp), DIMENSION(nau2)   ::  auy2   !< 
       REAL(wp), DIMENSION(nau2)   ::  aux4   !< 
       REAL(wp), DIMENSION(nau2)   ::  auy4   !<
#elif defined( __nec_fft )
       REAL(wp), DIMENSION(0:nx+3,nz+1)   ::  work_x  !<
       REAL(wp), DIMENSION(0:ny+3,nz+1)   ::  work_y  !<
       REAL(wp), DIMENSION(6*(nx+3),nz+1) ::  workx   !<
       REAL(wp), DIMENSION(6*(ny+3),nz+1) ::  worky   !<
#endif 

!
!--    Return, if already called
       IF ( init_fft )  THEN
          RETURN
       ELSE
          init_fft = .TRUE.
       ENDIF

#if defined( _OPENACC ) && defined( __cuda_fft )
       fft_method = 'system-specific'
#endif

!
!--    Switch to tell the Poisson-solver that the vectorized version of Temperton-fft is to be used.
       IF ( fft_method == 'temperton-algorithm'  .AND.  loop_optimization == 'vector'  .AND.       &
            pdims(1) /= 1  .AND.  pdims(2) /= 1 )  THEN
          temperton_fft_vec = .TRUE.
       ENDIF

       IF ( fft_method == 'system-specific' )  THEN

          dnx = 1.0_wp / ( nx + 1.0_wp )
          dny = 1.0_wp / ( ny + 1.0_wp )
          sqr_dnx = SQRT( dnx )
          sqr_dny = SQRT( dny )
#if defined( __ibm )
!
!--       Initialize tables for fft along x
          CALL DRCFT( 1, workx, 1, workx, 1, nx+1, 1,  1, sqr_dnx, aux1, nau1, &
                      aux2, nau2 )
          CALL DCRFT( 1, workx, 1, workx, 1, nx+1, 1, -1, sqr_dnx, aux3, nau1, &
                      aux4, nau2 )
!
!--       Initialize tables for fft along y
          CALL DRCFT( 1, worky, 1, worky, 1, ny+1, 1,  1, sqr_dny, auy1, nau1, &
                      auy2, nau2 )
          CALL DCRFT( 1, worky, 1, worky, 1, ny+1, 1, -1, sqr_dny, auy3, nau1, &
                      auy4, nau2 )
#elif defined( __nec_fft )
          message_string = 'fft method "' // TRIM( fft_method) // &
                           '" currently does not work on NEC'
          CALL message( 'fft_init', 'PA0187', 1, 2, 0, 6, 0 )

          ALLOCATE( trig_xb(2*(nx+1)), trig_xf(2*(nx+1)),                      &
                    trig_yb(2*(ny+1)), trig_yf(2*(ny+1)) )

          work_x = 0.0_wp
          work_y = 0.0_wp
          nz1  = nz + MOD( nz+1, 2 )  ! odd nz slows down fft significantly
                                      ! when using the NEC ffts

!
!--       Initialize tables for fft along x (non-vector and vector case (M))
          CALL DZFFT( 0, nx+1, sqr_dnx, work_x, work_x, trig_xf, workx, 0 )
          CALL ZDFFT( 0, nx+1, sqr_dnx, work_x, work_x, trig_xb, workx, 0 )
          CALL DZFFTM( 0, nx+1, nz1, sqr_dnx, work_x, nx+4, work_x, nx+4,      &
                       trig_xf, workx, 0 )
          CALL ZDFFTM( 0, nx+1, nz1, sqr_dnx, work_x, nx+4, work_x, nx+4,      &
                       trig_xb, workx, 0 )
!
!--       Initialize tables for fft along y (non-vector and vector case (M))
          CALL DZFFT( 0, ny+1, sqr_dny, work_y, work_y, trig_yf, worky, 0 )
          CALL ZDFFT( 0, ny+1, sqr_dny, work_y, work_y, trig_yb, worky, 0 )
          CALL DZFFTM( 0, ny+1, nz1, sqr_dny, work_y, ny+4, work_y, ny+4,      &
                       trig_yf, worky, 0 )
          CALL ZDFFTM( 0, ny+1, nz1, sqr_dny, work_y, ny+4, work_y, ny+4,      &
                       trig_yb, worky, 0 )
#elif defined( __cuda_fft )
          CALL CUFFTPLAN1D( plan_xf, nx+1, CUFFT_D2Z, (nyn_x-nys_x+1) * (nzt_x-nzb_x+1) )
          CALL CUFFTPLAN1D( plan_xi, nx+1, CUFFT_Z2D, (nyn_x-nys_x+1) * (nzt_x-nzb_x+1) )
          CALL CUFFTPLAN1D( plan_yf, ny+1, CUFFT_D2Z, (nxr_y-nxl_y+1) * (nzt_y-nzb_y+1) )
          CALL CUFFTPLAN1D( plan_yi, ny+1, CUFFT_Z2D, (nxr_y-nxl_y+1) * (nzt_y-nzb_y+1) )
#else
          message_string = 'no system-specific fft-call available'
          CALL message( 'fft_init', 'PA0188', 1, 2, 0, 6, 0 )
#endif
       ELSEIF ( fft_method == 'temperton-algorithm' )  THEN
!
!--       Temperton-algorithm
!--       Initialize tables for fft along x and y
          ALLOCATE( ifax_x(nx+1), ifax_y(ny+1), trigs_x(nx+1), trigs_y(ny+1) )

          CALL set99( trigs_x, ifax_x, nx+1 )
          CALL set99( trigs_y, ifax_y, ny+1 )

          IF ( temperton_fft_vec )  THEN
             ALLOCATE( f_vec_x((nyn_x-nys_x+1)*(nzt_x-nzb_x+1),0:nx+2) )
          ENDIF



       ELSEIF ( fft_method == 'fftw' )  THEN
!
!--       FFTW
#if defined( __fftw )
          nx_c = nx+1
          ny_c = ny+1
          !$OMP PARALLEL
          ALLOCATE( x_in(0:nx+2), y_in(0:ny+2), x_out(0:(nx+1)/2),             &
                    y_out(0:(ny+1)/2) )
          !$OMP END PARALLEL
          plan_xf = FFTW_PLAN_DFT_R2C_1D( nx_c, x_in, x_out, FFTW_ESTIMATE )
          plan_xi = FFTW_PLAN_DFT_C2R_1D( nx_c, x_out, x_in, FFTW_ESTIMATE )
          plan_yf = FFTW_PLAN_DFT_R2C_1D( ny_c, y_in, y_out, FFTW_ESTIMATE )
          plan_yi = FFTW_PLAN_DFT_C2R_1D( ny_c, y_out, y_in, FFTW_ESTIMATE )
#else
          message_string = 'preprocessor switch for fftw is missing'
          CALL message( 'fft_init', 'PA0080', 1, 2, 0, 6, 0 )
#endif

       ELSEIF ( fft_method == 'singleton-algorithm' )  THEN

          CONTINUE

       ELSE

          message_string = 'fft method "' // TRIM( fft_method) // &
                           '" not available'
          CALL message( 'fft_init', 'PA0189', 1, 2, 0, 6, 0 )
       ENDIF

    END SUBROUTINE fft_init


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Fourier-transformation along x-direction.                 
!> Version for 2D-decomposition.
!> It uses internal algorithms (Singleton or Temperton) or      
!> system-specific routines, if they are available           
!------------------------------------------------------------------------------!
 
    SUBROUTINE fft_x( ar, direction, ar_2d, ar_inv )


       IMPLICIT NONE

       CHARACTER (LEN=*) ::  direction  !<
       
       COMPLEX(wp), DIMENSION(:), ALLOCATABLE ::  cwork  !<

       INTEGER(iwp) ::  i          !< 
       INTEGER(iwp) ::  ishape(1)  !<
       INTEGER(iwp) ::  j          !<
       INTEGER(iwp) ::  k          !<
       INTEGER(iwp) ::  mm         !<

       LOGICAL ::  forward_fft !<
       
       REAL(wp), DIMENSION(0:nx+2) ::  work   !<
       REAL(wp), DIMENSION(nx+2)   ::  work1  !<
       
       REAL(wp), DIMENSION(:,:), ALLOCATABLE           ::  work_vec  !<
       REAL(wp), DIMENSION(0:nx,nys_x:nyn_x), OPTIONAL ::  ar_2d     !<

       REAL(wp), DIMENSION(nys_x:nyn_x,nzb_x:nzt_x,0:nx), OPTIONAL ::  ar_inv   !<
       REAL(wp), DIMENSION(0:nx,nys_x:nyn_x,nzb_x:nzt_x)           ::  ar       !<

#if defined( __ibm )
       REAL(wp), DIMENSION(nau2) ::  aux2  !< 
       REAL(wp), DIMENSION(nau2) ::  aux4  !<
#elif defined( __nec_fft )
       REAL(wp), DIMENSION(6*(nx+1)) ::  work2  !<
#elif defined( __cuda_fft )
       COMPLEX(dp), DIMENSION(0:(nx+1)/2,nys_x:nyn_x,nzb_x:nzt_x) ::  ar_tmp  !<
       !$ACC DECLARE CREATE(ar_tmp)
#endif

!
!--    To avoid compiler warning: Unused dummy argument ‘ar_2d’
       IF ( PRESENT( ar_2d ) )  CONTINUE

       IF ( direction == 'forward' )  THEN
          forward_fft = .TRUE.
       ELSE
          forward_fft = .FALSE.
       ENDIF

       IF ( fft_method == 'singleton-algorithm' )  THEN

!
!--       Performing the fft with singleton's software works on every system,
!--       since it is part of the model
          ALLOCATE( cwork(0:nx) )
      
          IF ( forward_fft )   then

             !$OMP PARALLEL PRIVATE ( cwork, i, ishape, j, k )
             !$OMP DO
             DO  k = nzb_x, nzt_x
                DO  j = nys_x, nyn_x

                   DO  i = 0, nx
                      cwork(i) = CMPLX( ar(i,j,k), KIND=wp )
                   ENDDO

                   ishape = SHAPE( cwork )
                   CALL FFTN( cwork, ishape )

                   DO  i = 0, (nx+1)/2
                      ar(i,j,k) = REAL( cwork(i), KIND=wp )
                   ENDDO
                   DO  i = 1, (nx+1)/2 - 1
                      ar(nx+1-i,j,k) = -AIMAG( cwork(i) )
                   ENDDO

                ENDDO
             ENDDO
             !$OMP END PARALLEL

          ELSE

             !$OMP PARALLEL PRIVATE ( cwork, i, ishape, j, k )
             !$OMP DO
             DO  k = nzb_x, nzt_x
                DO  j = nys_x, nyn_x

                   cwork(0) = CMPLX( ar(0,j,k), 0.0_wp, KIND=wp )
                   DO  i = 1, (nx+1)/2 - 1
                      cwork(i)      = CMPLX( ar(i,j,k), -ar(nx+1-i,j,k),       &
                                             KIND=wp )
                      cwork(nx+1-i) = CMPLX( ar(i,j,k),  ar(nx+1-i,j,k),       &
                                             KIND=wp )
                   ENDDO
                   cwork((nx+1)/2) = CMPLX( ar((nx+1)/2,j,k), 0.0_wp, KIND=wp )

                   ishape = SHAPE( cwork )
                   CALL FFTN( cwork, ishape, inv = .TRUE. )

                   DO  i = 0, nx
                      ar(i,j,k) = REAL( cwork(i), KIND=wp )
                   ENDDO

                ENDDO
             ENDDO
             !$OMP END PARALLEL

          ENDIF

          DEALLOCATE( cwork )

       ELSEIF ( fft_method == 'temperton-algorithm' )  THEN

!
!--       Performing the fft with Temperton's software works on every system,
!--       since it is part of the model
          IF ( forward_fft )  THEN

             IF ( .NOT. temperton_fft_vec )  THEN

                !$OMP PARALLEL PRIVATE ( work, work1, i, j, k )
                !$OMP DO
                DO  k = nzb_x, nzt_x
                   DO  j = nys_x, nyn_x

                      work(0:nx) = ar(0:nx,j,k)
                      CALL fft991cy( work, work1, trigs_x, ifax_x, 1, nx+1, nx+1, 1, -1 )

                      DO  i = 0, (nx+1)/2
                         ar(i,j,k) = work(2*i)
                      ENDDO
                      DO  i = 1, (nx+1)/2 - 1
                         ar(nx+1-i,j,k) = work(2*i+1)
                      ENDDO

                   ENDDO
                ENDDO
                !$OMP END PARALLEL

             ELSE

!
!--             Vector version of the Temperton-algorithm. Computes multiple 1-D FFT's.
                ALLOCATE( work_vec( (nyn_x-nys_x+1)*(nzt_x-nzb_x+1),nx+2) )
!
!--             f_vec_x is already set in transpose_zx
                CALL fft991cy_vec( f_vec_x, work_vec, trigs_x, ifax_x, nx+1, -1 )
                DEALLOCATE( work_vec )

                IF ( PRESENT( ar_inv ) )  THEN

                   DO  k = nzb_x, nzt_x
                      DO  j = nys_x, nyn_x
                         mm = j-nys_x+1+(k-nzb_x)*(nyn_x-nys_x+1)
                         DO  i = 0, (nx+1)/2
                            ar_inv(j,k,i) = f_vec_x(mm,2*i)
                         ENDDO
                         DO  i = 1, (nx+1)/2-1
                            ar_inv(j,k,nx+1-i) = f_vec_x(mm,2*i+1)
                         ENDDO
                      ENDDO
                   ENDDO

                ELSE

                   DO  k = nzb_x, nzt_x
                      DO  j = nys_x, nyn_x
                         mm = j-nys_x+1+(k-nzb_x)*(nyn_x-nys_x+1)
                         DO  i = 0, (nx+1)/2
                            ar(i,j,k) = f_vec_x(mm,2*i)
                         ENDDO
                         DO  i = 1, (nx+1)/2-1
                            ar(nx+1-i,j,k) = f_vec_x(mm,2*i+1)
                         ENDDO
                      ENDDO
                   ENDDO

                ENDIF

             ENDIF

          ELSE

!
!--          Backward fft
             IF ( .NOT. temperton_fft_vec )  THEN

                !$OMP PARALLEL PRIVATE ( work, work1, i, j, k )
                !$OMP DO
                DO  k = nzb_x, nzt_x
                   DO  j = nys_x, nyn_x

                      DO  i = 0, (nx+1)/2
                         work(2*i) = ar(i,j,k)
                      ENDDO
                      DO  i = 1, (nx+1)/2 - 1
                         work(2*i+1) = ar(nx+1-i,j,k)
                      ENDDO
                      work(1)    = 0.0_wp
                      work(nx+2) = 0.0_wp

                      CALL fft991cy( work, work1, trigs_x, ifax_x, 1, nx+1, nx+1, 1, 1 )
                      ar(0:nx,j,k) = work(0:nx)

                   ENDDO
                ENDDO
                !$OMP END PARALLEL

             ELSE

                IF ( PRESENT( ar_inv ) )  THEN

                   DO  k = nzb_x, nzt_x
                      DO  j = nys_x, nyn_x
                         mm = j-nys_x+1+(k-nzb_x)*(nyn_x-nys_x+1)
                         DO  i = 0, (nx+1)/2
                            f_vec_x(mm,2*i) = ar_inv(j,k,i)
                         ENDDO
                         DO  i = 1, (nx+1)/2-1
                            f_vec_x(mm,2*i+1) = ar_inv(j,k,nx+1-i)
                         ENDDO
                      ENDDO
                   ENDDO

                ELSE

                   DO  k = nzb_x, nzt_x
                      DO  j = nys_x, nyn_x
                         mm = j-nys_x+1+(k-nzb_x)*(nyn_x-nys_x+1)
                         DO  i = 0, (nx+1)/2
                            f_vec_x(mm,2*i) = ar(i,j,k)
                         ENDDO
                         DO  i = 1, (nx+1)/2-1
                            f_vec_x(mm,2*i+1) = ar(nx+1-i,j,k)
                         ENDDO
                      ENDDO
                   ENDDO

                ENDIF
                f_vec_x(:,1)    = 0.0_wp
                f_vec_x(:,nx+2) = 0.0_wp

                ALLOCATE( work_vec((nyn_x-nys_x+1)*(nzt_x-nzb_x+1),nx+2) )
                CALL fft991cy_vec( f_vec_x, work_vec, trigs_x, ifax_x, nx+1, 1 )
                DEALLOCATE( work_vec )

             ENDIF

          ENDIF

       ELSEIF ( fft_method == 'fftw' )  THEN

#if defined( __fftw )
          IF ( forward_fft )  THEN

             !$OMP PARALLEL PRIVATE ( work, i, j, k )
             !$OMP DO
             DO  k = nzb_x, nzt_x
                DO  j = nys_x, nyn_x

                   x_in(0:nx) = ar(0:nx,j,k)
                   CALL FFTW_EXECUTE_DFT_R2C( plan_xf, x_in, x_out )

                   IF ( PRESENT( ar_2d ) )  THEN

                      DO  i = 0, (nx+1)/2
                         ar_2d(i,j) = REAL( x_out(i), KIND=wp ) / ( nx+1 )
                      ENDDO
                      DO  i = 1, (nx+1)/2 - 1
                         ar_2d(nx+1-i,j) = AIMAG( x_out(i) ) / ( nx+1 )
                      ENDDO

                   ELSE

                      DO  i = 0, (nx+1)/2
                         ar(i,j,k) = REAL( x_out(i), KIND=wp ) / ( nx+1 )
                      ENDDO
                      DO  i = 1, (nx+1)/2 - 1
                         ar(nx+1-i,j,k) = AIMAG( x_out(i) ) / ( nx+1 )
                      ENDDO

                   ENDIF

                ENDDO
             ENDDO
             !$OMP END PARALLEL

          ELSE
             !$OMP PARALLEL PRIVATE ( work, i, j, k )
             !$OMP DO
             DO  k = nzb_x, nzt_x
                DO  j = nys_x, nyn_x

                   IF ( PRESENT( ar_2d ) )  THEN

                      x_out(0) = CMPLX( ar_2d(0,j), 0.0_wp, KIND=wp )
                      DO  i = 1, (nx+1)/2 - 1
                         x_out(i) = CMPLX( ar_2d(i,j), ar_2d(nx+1-i,j),        &
                                           KIND=wp )
                      ENDDO
                      x_out((nx+1)/2) = CMPLX( ar_2d((nx+1)/2,j), 0.0_wp,      &
                                               KIND=wp )

                   ELSE

                      x_out(0) = CMPLX( ar(0,j,k), 0.0_wp, KIND=wp )
                      DO  i = 1, (nx+1)/2 - 1
                         x_out(i) = CMPLX( ar(i,j,k), ar(nx+1-i,j,k), KIND=wp )
                      ENDDO
                      x_out((nx+1)/2) = CMPLX( ar((nx+1)/2,j,k), 0.0_wp,       &
                                               KIND=wp )

                   ENDIF

                   CALL FFTW_EXECUTE_DFT_C2R( plan_xi, x_out, x_in)
                   ar(0:nx,j,k) = x_in(0:nx)

                ENDDO
             ENDDO
             !$OMP END PARALLEL

          ENDIF
#endif

       ELSEIF ( fft_method == 'system-specific' )  THEN

#if defined( __ibm )
          IF ( forward_fft )  THEN

             !$OMP PARALLEL PRIVATE ( work, i, j, k )
             !$OMP DO
             DO  k = nzb_x, nzt_x
                DO  j = nys_x, nyn_x

                   CALL DRCFT( 0, ar, 1, work, 1, nx+1, 1, 1, sqr_dnx, aux1,   &
                               nau1, aux2, nau2 )

                   DO  i = 0, (nx+1)/2
                      ar(i,j,k) = work(2*i)
                   ENDDO
                   DO  i = 1, (nx+1)/2 - 1
                      ar(nx+1-i,j,k) = work(2*i+1)
                   ENDDO

                ENDDO
             ENDDO
             !$OMP END PARALLEL

          ELSE

             !$OMP PARALLEL PRIVATE ( work, i, j, k )
             !$OMP DO
             DO  k = nzb_x, nzt_x
                DO  j = nys_x, nyn_x

                   DO  i = 0, (nx+1)/2
                      work(2*i) = ar(i,j,k)
                   ENDDO
                   DO  i = 1, (nx+1)/2 - 1
                      work(2*i+1) = ar(nx+1-i,j,k)
                   ENDDO
                   work(1) = 0.0_wp
                   work(nx+2) = 0.0_wp

                   CALL DCRFT( 0, work, 1, work, 1, nx+1, 1, -1, sqr_dnx,      & 
                               aux3, nau1, aux4, nau2 )

                   DO  i = 0, nx
                      ar(i,j,k) = work(i)
                   ENDDO

                ENDDO
             ENDDO
             !$OMP END PARALLEL

          ENDIF

#elif defined( __nec_fft )

          IF ( forward_fft )  THEN

             !$OMP PARALLEL PRIVATE ( work, i, j, k )
             !$OMP DO
             DO  k = nzb_x, nzt_x
                DO  j = nys_x, nyn_x

                   work(0:nx) = ar(0:nx,j,k)

                   CALL DZFFT( 1, nx+1, sqr_dnx, work, work, trig_xf, work2, 0 )
      
                   DO  i = 0, (nx+1)/2
                      ar(i,j,k) = work(2*i)
                   ENDDO
                   DO  i = 1, (nx+1)/2 - 1
                      ar(nx+1-i,j,k) = work(2*i+1)
                   ENDDO

                ENDDO
             ENDDO
             !$END OMP PARALLEL

          ELSE

             !$OMP PARALLEL PRIVATE ( work, i, j, k )
             !$OMP DO
             DO  k = nzb_x, nzt_x
                DO  j = nys_x, nyn_x

                   DO  i = 0, (nx+1)/2
                      work(2*i) = ar(i,j,k)
                   ENDDO
                   DO  i = 1, (nx+1)/2 - 1
                      work(2*i+1) = ar(nx+1-i,j,k)
                   ENDDO
                   work(1) = 0.0_wp
                   work(nx+2) = 0.0_wp

                   CALL ZDFFT( -1, nx+1, sqr_dnx, work, work, trig_xb, work2, 0 )

                   ar(0:nx,j,k) = work(0:nx)

                ENDDO
             ENDDO
             !$OMP END PARALLEL

          ENDIF

#elif defined( __cuda_fft )

          IF ( forward_fft )  THEN

             !$ACC HOST_DATA USE_DEVICE(ar, ar_tmp)
             CALL CUFFTEXECD2Z( plan_xf, ar, ar_tmp )
             !$ACC END HOST_DATA

             !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(i,j,k) &
             !$ACC PRESENT(ar, ar_tmp)
             DO  k = nzb_x, nzt_x
                DO  j = nys_x, nyn_x

                   DO  i = 0, (nx+1)/2
                      ar(i,j,k)      = REAL( ar_tmp(i,j,k), KIND=wp )  * dnx
                   ENDDO

                   DO  i = 1, (nx+1)/2 - 1
                      ar(nx+1-i,j,k) = AIMAG( ar_tmp(i,j,k) ) * dnx
                   ENDDO

                ENDDO
             ENDDO

          ELSE

             !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(i,j,k) &
             !$ACC PRESENT(ar, ar_tmp)
             DO  k = nzb_x, nzt_x
                DO  j = nys_x, nyn_x

                   ar_tmp(0,j,k) = CMPLX( ar(0,j,k), 0.0_wp, KIND=wp )

                   DO  i = 1, (nx+1)/2 - 1
                      ar_tmp(i,j,k) = CMPLX( ar(i,j,k), ar(nx+1-i,j,k),        &
                                             KIND=wp )
                   ENDDO
                   ar_tmp((nx+1)/2,j,k) = CMPLX( ar((nx+1)/2,j,k), 0.0_wp,     &
                                                 KIND=wp )

                ENDDO
             ENDDO

             !$ACC HOST_DATA USE_DEVICE(ar, ar_tmp)
             CALL CUFFTEXECZ2D( plan_xi, ar_tmp, ar )
             !$ACC END HOST_DATA

          ENDIF

#endif

       ENDIF

    END SUBROUTINE fft_x

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Fourier-transformation along x-direction.
!> Version for 1D-decomposition.
!> It uses internal algorithms (Singleton or Temperton) or
!> system-specific routines, if they are available
!------------------------------------------------------------------------------!
 
    SUBROUTINE fft_x_1d( ar, direction )


       IMPLICIT NONE

       CHARACTER (LEN=*) ::  direction  !<
       
       INTEGER(iwp) ::  i               !<
       INTEGER(iwp) ::  ishape(1)       !<

       LOGICAL ::  forward_fft          !<

       REAL(wp), DIMENSION(0:nx)   ::  ar     !<
       REAL(wp), DIMENSION(0:nx+2) ::  work   !<
       REAL(wp), DIMENSION(nx+2)   ::  work1  !<
       
       COMPLEX(wp), DIMENSION(:), ALLOCATABLE ::  cwork  !<
       
#if defined( __ibm )
       REAL(wp), DIMENSION(nau2) ::  aux2       !<
       REAL(wp), DIMENSION(nau2) ::  aux4       !<
#elif defined( __nec_fft )
       REAL(wp), DIMENSION(6*(nx+1)) ::  work2  !<
#endif

       IF ( direction == 'forward' )  THEN
          forward_fft = .TRUE.
       ELSE
          forward_fft = .FALSE.
       ENDIF

       IF ( fft_method == 'singleton-algorithm' )  THEN

!
!--       Performing the fft with singleton's software works on every system,
!--       since it is part of the model
          ALLOCATE( cwork(0:nx) )
      
          IF ( forward_fft )   then

             DO  i = 0, nx
                cwork(i) = CMPLX( ar(i), KIND=wp )
             ENDDO
             ishape = SHAPE( cwork )
             CALL FFTN( cwork, ishape )
             DO  i = 0, (nx+1)/2
                ar(i) = REAL( cwork(i), KIND=wp )
             ENDDO
             DO  i = 1, (nx+1)/2 - 1
                ar(nx+1-i) = -AIMAG( cwork(i) )
             ENDDO

          ELSE

             cwork(0) = CMPLX( ar(0), 0.0_wp, KIND=wp )
             DO  i = 1, (nx+1)/2 - 1
                cwork(i)      = CMPLX( ar(i), -ar(nx+1-i), KIND=wp )
                cwork(nx+1-i) = CMPLX( ar(i),  ar(nx+1-i), KIND=wp )
             ENDDO
             cwork((nx+1)/2) = CMPLX( ar((nx+1)/2), 0.0_wp, KIND=wp )

             ishape = SHAPE( cwork )
             CALL FFTN( cwork, ishape, inv = .TRUE. )

             DO  i = 0, nx
                ar(i) = REAL( cwork(i), KIND=wp )
             ENDDO

          ENDIF

          DEALLOCATE( cwork )

       ELSEIF ( fft_method == 'temperton-algorithm' )  THEN

!
!--       Performing the fft with Temperton's software works on every system,
!--       since it is part of the model
          IF ( forward_fft )  THEN

             work(0:nx) = ar
             CALL fft991cy( work, work1, trigs_x, ifax_x, 1, nx+1, nx+1, 1, -1 )

             DO  i = 0, (nx+1)/2
                ar(i) = work(2*i)
             ENDDO
             DO  i = 1, (nx+1)/2 - 1
                ar(nx+1-i) = work(2*i+1)
             ENDDO

          ELSE

             DO  i = 0, (nx+1)/2
                work(2*i) = ar(i)
             ENDDO
             DO  i = 1, (nx+1)/2 - 1
                work(2*i+1) = ar(nx+1-i)
             ENDDO
             work(1)    = 0.0_wp
             work(nx+2) = 0.0_wp

             CALL fft991cy( work, work1, trigs_x, ifax_x, 1, nx+1, nx+1, 1, 1 )
             ar = work(0:nx)

          ENDIF

       ELSEIF ( fft_method == 'fftw' )  THEN

#if defined( __fftw )
          IF ( forward_fft )  THEN

             x_in(0:nx) = ar(0:nx)
             CALL FFTW_EXECUTE_DFT_R2C( plan_xf, x_in, x_out )

             DO  i = 0, (nx+1)/2
                ar(i) = REAL( x_out(i), KIND=wp ) / ( nx+1 )
             ENDDO
             DO  i = 1, (nx+1)/2 - 1
                ar(nx+1-i) = AIMAG( x_out(i) ) / ( nx+1 )
             ENDDO

         ELSE

             x_out(0) = CMPLX( ar(0), 0.0_wp, KIND=wp )
             DO  i = 1, (nx+1)/2 - 1
                x_out(i) = CMPLX( ar(i), ar(nx+1-i), KIND=wp )
             ENDDO
             x_out((nx+1)/2) = CMPLX( ar((nx+1)/2), 0.0_wp, KIND=wp )

             CALL FFTW_EXECUTE_DFT_C2R( plan_xi, x_out, x_in)
             ar(0:nx) = x_in(0:nx)

         ENDIF
#endif

       ELSEIF ( fft_method == 'system-specific' )  THEN

#if defined( __ibm )
          IF ( forward_fft )  THEN

             CALL DRCFT( 0, ar, 1, work, 1, nx+1, 1, 1, sqr_dnx, aux1, nau1,   &
                         aux2, nau2 )

             DO  i = 0, (nx+1)/2
                ar(i) = work(2*i)
             ENDDO
             DO  i = 1, (nx+1)/2 - 1
                ar(nx+1-i) = work(2*i+1)
             ENDDO

          ELSE

             DO  i = 0, (nx+1)/2
                work(2*i) = ar(i)
             ENDDO
             DO  i = 1, (nx+1)/2 - 1
                work(2*i+1) = ar(nx+1-i)
             ENDDO
             work(1) = 0.0_wp
             work(nx+2) = 0.0_wp

             CALL DCRFT( 0, work, 1, work, 1, nx+1, 1, -1, sqr_dnx, aux3, nau1, &
                         aux4, nau2 )

             DO  i = 0, nx
                ar(i) = work(i)
             ENDDO

          ENDIF
#elif defined( __nec_fft )
          IF ( forward_fft )  THEN

             work(0:nx) = ar(0:nx)

             CALL DZFFT( 1, nx+1, sqr_dnx, work, work, trig_xf, work2, 0 )
      
             DO  i = 0, (nx+1)/2
                ar(i) = work(2*i)
             ENDDO
             DO  i = 1, (nx+1)/2 - 1
                ar(nx+1-i) = work(2*i+1)
             ENDDO

          ELSE

             DO  i = 0, (nx+1)/2
                work(2*i) = ar(i)
             ENDDO
             DO  i = 1, (nx+1)/2 - 1
                work(2*i+1) = ar(nx+1-i)
             ENDDO
             work(1) = 0.0_wp
             work(nx+2) = 0.0_wp

             CALL ZDFFT( -1, nx+1, sqr_dnx, work, work, trig_xb, work2, 0 )

             ar(0:nx) = work(0:nx)

          ENDIF
#endif

       ENDIF

    END SUBROUTINE fft_x_1d

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Fourier-transformation along y-direction.
!> Version for 2D-decomposition.
!> It uses internal algorithms (Singleton or Temperton) or
!> system-specific routines, if they are available.
!> 
!> direction:  'forward' or 'backward'
!> ar, ar_tr:  3D data arrays 
!>             forward:   ar: before  ar_tr: after transformation
!>             backward:  ar_tr: before  ar: after transfosition
!> 
!> In case of non-overlapping transposition/transformation:
!> nxl_y_bound = nxl_y_l = nxl_y 
!> nxr_y_bound = nxr_y_l = nxr_y 
!> 
!> In case of overlapping transposition/transformation
!> - nxl_y_bound  and  nxr_y_bound have the original values of
!>   nxl_y, nxr_y.  ar_tr is dimensioned using these values.
!> - nxl_y_l = nxr_y_r.  ar is dimensioned with these values, so that
!>   transformation is carried out for a 2D-plane only.
!------------------------------------------------------------------------------!
 
    SUBROUTINE fft_y( ar, direction, ar_tr, nxl_y_bound, nxr_y_bound, nxl_y_l, &
                      nxr_y_l, ar_inv )


       IMPLICIT NONE

       CHARACTER (LEN=*) ::  direction  !<
       
       INTEGER(iwp) ::  i            !<
       INTEGER(iwp) ::  j            !< 
       INTEGER(iwp) ::  jshape(1)    !<
       INTEGER(iwp) ::  k            !<
       INTEGER(iwp) ::  mm           !<
       INTEGER(iwp) ::  nxl_y_bound  !<
       INTEGER(iwp) ::  nxl_y_l      !<
       INTEGER(iwp) ::  nxr_y_bound  !<
       INTEGER(iwp) ::  nxr_y_l      !<

       LOGICAL ::  forward_fft  !<

       REAL(wp), DIMENSION(0:ny+2) ::  work   !<
       REAL(wp), DIMENSION(ny+2)   ::  work1  !<
       
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  f_vec_y
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  work_vec

       REAL(wp), DIMENSION(0:ny,nxl_y_l:nxr_y_l,nzb_y:nzt_y)                   ::  ar      !<
       REAL(wp), DIMENSION(nxl_y:nxr_y,nzb_y:nzt_y,0:ny), OPTIONAL             ::  ar_inv  !<
       REAL(wp), DIMENSION(0:ny,nxl_y_bound:nxr_y_bound,nzb_y:nzt_y), OPTIONAL ::  ar_tr   !<

       COMPLEX(wp), DIMENSION(:), ALLOCATABLE ::  cwork  !<
       
#if defined( __ibm )
       REAL(wp), DIMENSION(nau2) ::  auy2  !<
       REAL(wp), DIMENSION(nau2) ::  auy4  !<
#elif defined( __nec_fft )
       REAL(wp), DIMENSION(6*(ny+1)) ::  work2  !<
#elif defined( __cuda_fft )
       COMPLEX(dp), DIMENSION(0:(ny+1)/2,nxl_y:nxr_y,nzb_y:nzt_y) ::           &
          ar_tmp  !<
       !$ACC DECLARE CREATE(ar_tmp)
#endif


       IF ( direction == 'forward' )  THEN
          forward_fft = .TRUE.
       ELSE
          forward_fft = .FALSE.
       ENDIF

       IF ( fft_method == 'singleton-algorithm' )  THEN

!
!--       Performing the fft with singleton's software works on every system,
!--       since it is part of the model
          ALLOCATE( cwork(0:ny) )

          IF ( forward_fft )   then

             !$OMP PARALLEL PRIVATE ( cwork, i, jshape, j, k )
             !$OMP DO
             DO  k = nzb_y, nzt_y
                DO  i = nxl_y_l, nxr_y_l

                   DO  j = 0, ny
                      cwork(j) = CMPLX( ar(j,i,k), KIND=wp )
                   ENDDO

                   jshape = SHAPE( cwork )
                   CALL FFTN( cwork, jshape )

                   DO  j = 0, (ny+1)/2
                      ar_tr(j,i,k) = REAL( cwork(j), KIND=wp )
                   ENDDO
                   DO  j = 1, (ny+1)/2 - 1
                      ar_tr(ny+1-j,i,k) = -AIMAG( cwork(j) )
                   ENDDO

                ENDDO
             ENDDO
             !$OMP END PARALLEL

          ELSE

             !$OMP PARALLEL PRIVATE ( cwork, i, jshape, j, k )
             !$OMP DO
             DO  k = nzb_y, nzt_y
                DO  i = nxl_y_l, nxr_y_l

                   cwork(0) = CMPLX( ar_tr(0,i,k), 0.0_wp, KIND=wp )
                   DO  j = 1, (ny+1)/2 - 1
                      cwork(j)      = CMPLX( ar_tr(j,i,k), -ar_tr(ny+1-j,i,k), &
                                             KIND=wp )
                      cwork(ny+1-j) = CMPLX( ar_tr(j,i,k),  ar_tr(ny+1-j,i,k), &
                                             KIND=wp )
                   ENDDO
                   cwork((ny+1)/2) = CMPLX( ar_tr((ny+1)/2,i,k), 0.0_wp,       &
                                            KIND=wp )

                   jshape = SHAPE( cwork )
                   CALL FFTN( cwork, jshape, inv = .TRUE. )

                   DO  j = 0, ny
                      ar(j,i,k) = REAL( cwork(j), KIND=wp )
                   ENDDO

                ENDDO
             ENDDO
             !$OMP END PARALLEL

          ENDIF

          DEALLOCATE( cwork )

       ELSEIF ( fft_method == 'temperton-algorithm' )  THEN

!
!--       Performing the fft with Temperton's software works on every system,
!--       since it is part of the model
          IF ( forward_fft )  THEN

             IF ( .NOT. temperton_fft_vec )  THEN

                !$OMP PARALLEL PRIVATE ( work, work1, i, j, k )
                !$OMP DO
                DO  k = nzb_y, nzt_y
                   DO  i = nxl_y_l, nxr_y_l

                      work(0:ny) = ar(0:ny,i,k)
                      CALL fft991cy( work, work1, trigs_y, ifax_y, 1, ny+1, ny+1, 1, -1 )

                      DO  j = 0, (ny+1)/2
                         ar_tr(j,i,k) = work(2*j)
                      ENDDO
                      DO  j = 1, (ny+1)/2 - 1
                         ar_tr(ny+1-j,i,k) = work(2*j+1)
                      ENDDO

                   ENDDO
                ENDDO
                !$OMP END PARALLEL

             ELSE
!
!--             Vector version of Temperton-fft. Computes multiple 1-D FFT's.
                ALLOCATE( f_vec_y((nxr_y_l-nxl_y_l+1)*(nzt_y-nzb_y+1),0:ny+2) )

                mm = 1
                DO  k = nzb_y, nzt_y
                   DO  i = nxl_y_l, nxr_y_l
                      f_vec_y(mm,0:nx) = ar(0:nx,i,k)
                      mm = mm+1
                   ENDDO
                ENDDO

                ALLOCATE( work_vec( (nxr_y_l-nxl_y_l+1)*(nzt_y-nzb_y+1),ny+2) )
                CALL fft991cy_vec( f_vec_y, work_vec, trigs_y, ifax_y, ny+1, -1 )
                DEALLOCATE( work_vec )

                IF( PRESENT( ar_inv ) )  THEN

                   DO  k = nzb_y, nzt_y
                      DO  i = nxl_y_l, nxr_y_l
                         mm = i-nxl_y_l+1+(k-nzb_y)*(nxr_y_l-nxl_y_l+1)
                         DO  j = 0, (ny+1)/2
                            ar_inv(i,k,j) = f_vec_y(mm,2*j)
                         ENDDO
                         DO  j = 1, (ny+1)/2 - 1
                            ar_inv(i,k,ny+1-j) = f_vec_y(mm,2*j+1)
                         ENDDO
                      ENDDO
                   ENDDO

                ELSE

                   DO  k = nzb_y, nzt_y
                      DO  i = nxl_y_l, nxr_y_l
                         mm = i-nxl_y_l+1+(k-nzb_y)*(nxr_y_l-nxl_y_l+1)
                         DO  j = 0, (ny+1)/2
                            ar(j,i,k) = f_vec_y(mm,2*j)
                         ENDDO
                         DO  j = 1, (ny+1)/2 - 1
                            ar(ny+1-j,i,k) = f_vec_y(mm,2*j+1)
                         ENDDO
                      ENDDO
                   ENDDO

                ENDIF

                DEALLOCATE( f_vec_y )

             ENDIF

          ELSE

             IF ( .NOT. temperton_fft_vec )  THEN

                !$OMP PARALLEL PRIVATE ( work, work1, i, j, k )
                !$OMP DO
                DO  k = nzb_y, nzt_y
                   DO  i = nxl_y_l, nxr_y_l

                      DO  j = 0, (ny+1)/2
                         work(2*j) = ar_tr(j,i,k)
                      ENDDO
                      DO  j = 1, (ny+1)/2 - 1
                         work(2*j+1) = ar_tr(ny+1-j,i,k)
                      ENDDO
                      work(1)    = 0.0_wp
                      work(ny+2) = 0.0_wp

                      CALL fft991cy( work, work1, trigs_y, ifax_y, 1, ny+1, ny+1, 1, 1 )
                      ar(0:ny,i,k) = work(0:ny)

                   ENDDO
                ENDDO
                !$OMP END PARALLEL

             ELSE

                ALLOCATE( f_vec_y((nxr_y_l-nxl_y_l+1)*(nzt_y-nzb_y+1),0:ny+2) )

                IF ( PRESENT( ar_inv ) )  THEN

                   DO  k = nzb_y, nzt_y
                      DO  i = nxl_y_l, nxr_y_l
                         mm = i-nxl_y_l+1+(k-nzb_y)*(nxr_y_l-nxl_y_l+1)
                         DO  j = 0, (ny+1)/2
                            f_vec_y(mm,2*j) = ar_inv(i,k,j)
                         ENDDO
                         DO  j = 1, (ny+1)/2 - 1
                            f_vec_y(mm,2*j+1) = ar_inv(i,k,ny+1-j)
                         ENDDO
                      ENDDO
                   ENDDO

                ELSE

                   DO  k = nzb_y, nzt_y
                      DO  i = nxl_y_l, nxr_y_l
                         mm = i-nxl_y_l+1+(k-nzb_y)*(nxr_y_l-nxl_y_l+1)
                         DO  j = 0, (ny+1)/2
                            f_vec_y(mm,2*j) = ar(j,i,k)
                         ENDDO
                         DO  j = 1, (ny+1)/2 - 1
                            f_vec_y(mm,2*j+1) = ar(ny+1-j,i,k)
                         ENDDO
                      ENDDO
                   ENDDO

                ENDIF

                f_vec_y(:,1)    = 0.0_wp
                f_vec_y(:,ny+2) = 0.0_wp

                ALLOCATE( work_vec((nxr_y_l-nxl_y_l+1)*(nzt_y-nzb_y+1),ny+2) )
                CALL fft991cy_vec( f_vec_y, work_vec, trigs_y, ifax_y, ny+1, 1 )
                DEALLOCATE( work_vec )

                mm = 1
                DO  k = nzb_y, nzt_y
                   DO  i = nxl_y_l, nxr_y_l
                      ar(0:ny,i,k) = f_vec_y(mm,0:ny)
                      mm = mm+1
                   ENDDO
                ENDDO

                DEALLOCATE( f_vec_y )

             ENDIF

          ENDIF

       ELSEIF ( fft_method == 'fftw' )  THEN

#if defined( __fftw )
          IF ( forward_fft )  THEN

             !$OMP PARALLEL PRIVATE ( work, i, j, k )
             !$OMP DO
             DO  k = nzb_y, nzt_y
                DO  i = nxl_y_l, nxr_y_l

                   y_in(0:ny) = ar(0:ny,i,k)
                   CALL FFTW_EXECUTE_DFT_R2C( plan_yf, y_in, y_out )

                   DO  j = 0, (ny+1)/2
                      ar_tr(j,i,k) = REAL( y_out(j), KIND=wp ) / (ny+1)
                   ENDDO
                   DO  j = 1, (ny+1)/2 - 1
                      ar_tr(ny+1-j,i,k) = AIMAG( y_out(j) ) / (ny+1)
                   ENDDO

                ENDDO
             ENDDO
             !$OMP END PARALLEL

          ELSE

             !$OMP PARALLEL PRIVATE ( work, i, j, k )
             !$OMP DO
             DO  k = nzb_y, nzt_y
                DO  i = nxl_y_l, nxr_y_l

                   y_out(0) = CMPLX( ar_tr(0,i,k), 0.0_wp, KIND=wp )
                   DO  j = 1, (ny+1)/2 - 1
                      y_out(j) = CMPLX( ar_tr(j,i,k), ar_tr(ny+1-j,i,k),       &
                                        KIND=wp )
                   ENDDO
                   y_out((ny+1)/2) = CMPLX( ar_tr((ny+1)/2,i,k), 0.0_wp,       &
                                            KIND=wp )

                   CALL FFTW_EXECUTE_DFT_C2R( plan_yi, y_out, y_in )
                   ar(0:ny,i,k) = y_in(0:ny)

                ENDDO
             ENDDO
             !$OMP END PARALLEL

          ENDIF
#endif

       ELSEIF ( fft_method == 'system-specific' )  THEN

#if defined( __ibm )
          IF ( forward_fft)  THEN

             !$OMP PARALLEL PRIVATE ( work, i, j, k )
             !$OMP DO
             DO  k = nzb_y, nzt_y
                DO  i = nxl_y_l, nxr_y_l

                   CALL DRCFT( 0, ar, 1, work, 1, ny+1, 1, 1, sqr_dny, auy1,   & 
                               nau1, auy2, nau2 )

                   DO  j = 0, (ny+1)/2
                      ar_tr(j,i,k) = work(2*j)
                   ENDDO
                   DO  j = 1, (ny+1)/2 - 1
                      ar_tr(ny+1-j,i,k) = work(2*j+1)
                   ENDDO

                ENDDO
             ENDDO
             !$OMP END PARALLEL

          ELSE

             !$OMP PARALLEL PRIVATE ( work, i, j, k )
             !$OMP DO
             DO  k = nzb_y, nzt_y
                DO  i = nxl_y_l, nxr_y_l

                   DO  j = 0, (ny+1)/2
                      work(2*j) = ar_tr(j,i,k)
                   ENDDO
                   DO  j = 1, (ny+1)/2 - 1
                      work(2*j+1) = ar_tr(ny+1-j,i,k)
                   ENDDO
                   work(1)    = 0.0_wp
                   work(ny+2) = 0.0_wp

                   CALL DCRFT( 0, work, 1, work, 1, ny+1, 1, -1, sqr_dny,      &
                               auy3, nau1, auy4, nau2 )

                   DO  j = 0, ny
                      ar(j,i,k) = work(j)
                   ENDDO

                ENDDO
             ENDDO
             !$OMP END PARALLEL

          ENDIF
#elif defined( __nec_fft )
          IF ( forward_fft )  THEN

             !$OMP PARALLEL PRIVATE ( work, i, j, k )
             !$OMP DO
             DO  k = nzb_y, nzt_y
                DO  i = nxl_y_l, nxr_y_l

                   work(0:ny) = ar(0:ny,i,k)

                   CALL DZFFT( 1, ny+1, sqr_dny, work, work, trig_yf, work2, 0 )

                   DO  j = 0, (ny+1)/2
                      ar_tr(j,i,k) = work(2*j)
                   ENDDO
                   DO  j = 1, (ny+1)/2 - 1
                      ar_tr(ny+1-j,i,k) = work(2*j+1)
                   ENDDO

                ENDDO
             ENDDO
             !$END OMP PARALLEL

          ELSE

             !$OMP PARALLEL PRIVATE ( work, i, j, k )
             !$OMP DO
             DO  k = nzb_y, nzt_y
                DO  i = nxl_y_l, nxr_y_l

                   DO  j = 0, (ny+1)/2
                      work(2*j) = ar_tr(j,i,k)
                   ENDDO
                   DO  j = 1, (ny+1)/2 - 1
                      work(2*j+1) = ar_tr(ny+1-j,i,k)
                   ENDDO
                   work(1) = 0.0_wp
                   work(ny+2) = 0.0_wp

                   CALL ZDFFT( -1, ny+1, sqr_dny, work, work, trig_yb, work2, 0 )

                   ar(0:ny,i,k) = work(0:ny)

                ENDDO
             ENDDO
             !$OMP END PARALLEL

          ENDIF
#elif defined( __cuda_fft )

          IF ( forward_fft )  THEN

             !$ACC HOST_DATA USE_DEVICE(ar, ar_tmp)
             CALL CUFFTEXECD2Z( plan_yf, ar, ar_tmp )
             !$ACC END HOST_DATA

             !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(i,j,k) &
             !$ACC PRESENT(ar, ar_tmp)
             DO  k = nzb_y, nzt_y
                DO  i = nxl_y, nxr_y

                   DO  j = 0, (ny+1)/2
                      ar(j,i,k)      = REAL( ar_tmp(j,i,k), KIND=wp )  * dny
                   ENDDO

                   DO  j = 1, (ny+1)/2 - 1
                      ar(ny+1-j,i,k) = AIMAG( ar_tmp(j,i,k) ) * dny
                   ENDDO

                ENDDO
             ENDDO

          ELSE

             !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(i,j,k) &
             !$ACC PRESENT(ar, ar_tmp)
             DO  k = nzb_y, nzt_y
                DO  i = nxl_y, nxr_y

                   ar_tmp(0,i,k) = CMPLX( ar(0,i,k), 0.0_wp, KIND=wp )

                   DO  j = 1, (ny+1)/2 - 1
                      ar_tmp(j,i,k) = CMPLX( ar(j,i,k), ar(ny+1-j,i,k),        &
                                             KIND=wp )
                   ENDDO
                   ar_tmp((ny+1)/2,i,k) = CMPLX( ar((ny+1)/2,i,k), 0.0_wp,     &
                                                 KIND=wp )

                ENDDO
             ENDDO

             !$ACC HOST_DATA USE_DEVICE(ar, ar_tmp)
             CALL CUFFTEXECZ2D( plan_yi, ar_tmp, ar )
             !$ACC END HOST_DATA

          ENDIF

#endif

       ENDIF

    END SUBROUTINE fft_y

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Fourier-transformation along y-direction.
!> Version for 1D-decomposition.
!> It uses internal algorithms (Singleton or Temperton) or
!> system-specific routines, if they are available.
!------------------------------------------------------------------------------!
 
    SUBROUTINE fft_y_1d( ar, direction )


       IMPLICIT NONE

       CHARACTER (LEN=*) ::  direction
       
       INTEGER(iwp) ::  j          !<
       INTEGER(iwp) ::  jshape(1)  !<

       LOGICAL ::  forward_fft  !<

       REAL(wp), DIMENSION(0:ny)    ::  ar     !<
       REAL(wp), DIMENSION(0:ny+2)  ::  work   !<
       REAL(wp), DIMENSION(ny+2)    ::  work1  !<
       
       COMPLEX(wp), DIMENSION(:), ALLOCATABLE ::  cwork  !<
       
#if defined( __ibm )
       REAL(wp), DIMENSION(nau2) ::  auy2  !<
       REAL(wp), DIMENSION(nau2) ::  auy4  !<
#elif defined( __nec_fft )
       REAL(wp), DIMENSION(6*(ny+1)) ::  work2  !<
#endif

       IF ( direction == 'forward' )  THEN
          forward_fft = .TRUE.
       ELSE
          forward_fft = .FALSE.
       ENDIF

       IF ( fft_method == 'singleton-algorithm' )  THEN

!
!--       Performing the fft with singleton's software works on every system,
!--       since it is part of the model
          ALLOCATE( cwork(0:ny) )

          IF ( forward_fft )  THEN

             DO  j = 0, ny
                cwork(j) = CMPLX( ar(j), KIND=wp )
             ENDDO

             jshape = SHAPE( cwork )
             CALL FFTN( cwork, jshape )

             DO  j = 0, (ny+1)/2
                ar(j) = REAL( cwork(j), KIND=wp )
             ENDDO
             DO  j = 1, (ny+1)/2 - 1
                ar(ny+1-j) = -AIMAG( cwork(j) )
             ENDDO

          ELSE

             cwork(0) = CMPLX( ar(0), 0.0_wp, KIND=wp )
             DO  j = 1, (ny+1)/2 - 1
                cwork(j)      = CMPLX( ar(j), -ar(ny+1-j), KIND=wp )
                cwork(ny+1-j) = CMPLX( ar(j),  ar(ny+1-j), KIND=wp )
             ENDDO
             cwork((ny+1)/2) = CMPLX( ar((ny+1)/2), 0.0_wp, KIND=wp )

             jshape = SHAPE( cwork )
             CALL FFTN( cwork, jshape, inv = .TRUE. )

             DO  j = 0, ny
                ar(j) = REAL( cwork(j), KIND=wp )
             ENDDO

          ENDIF

          DEALLOCATE( cwork )

       ELSEIF ( fft_method == 'temperton-algorithm' )  THEN

!
!--       Performing the fft with Temperton's software works on every system,
!--       since it is part of the model
          IF ( forward_fft )  THEN

             work(0:ny) = ar
             CALL fft991cy( work, work1, trigs_y, ifax_y, 1, ny+1, ny+1, 1, -1 )

             DO  j = 0, (ny+1)/2
                ar(j) = work(2*j)
             ENDDO
             DO  j = 1, (ny+1)/2 - 1
                ar(ny+1-j) = work(2*j+1)
             ENDDO

          ELSE

             DO  j = 0, (ny+1)/2
                work(2*j) = ar(j)
             ENDDO
             DO  j = 1, (ny+1)/2 - 1
                work(2*j+1) = ar(ny+1-j)
             ENDDO
             work(1)    = 0.0_wp
             work(ny+2) = 0.0_wp

             CALL fft991cy( work, work1, trigs_y, ifax_y, 1, ny+1, ny+1, 1, 1 )
             ar = work(0:ny)

          ENDIF

       ELSEIF ( fft_method == 'fftw' )  THEN

#if defined( __fftw )
          IF ( forward_fft )  THEN

             y_in(0:ny) = ar(0:ny)
             CALL FFTW_EXECUTE_DFT_R2C( plan_yf, y_in, y_out )

             DO  j = 0, (ny+1)/2
                ar(j) = REAL( y_out(j), KIND=wp ) / (ny+1)
             ENDDO
             DO  j = 1, (ny+1)/2 - 1
                ar(ny+1-j) = AIMAG( y_out(j) ) / (ny+1)
             ENDDO

          ELSE

             y_out(0) = CMPLX( ar(0), 0.0_wp, KIND=wp )
             DO  j = 1, (ny+1)/2 - 1
                y_out(j) = CMPLX( ar(j), ar(ny+1-j), KIND=wp )
             ENDDO
             y_out((ny+1)/2) = CMPLX( ar((ny+1)/2), 0.0_wp, KIND=wp )

             CALL FFTW_EXECUTE_DFT_C2R( plan_yi, y_out, y_in )
             ar(0:ny) = y_in(0:ny)

          ENDIF
#endif

       ELSEIF ( fft_method == 'system-specific' )  THEN

#if defined( __ibm )
          IF ( forward_fft )  THEN

             CALL DRCFT( 0, ar, 1, work, 1, ny+1, 1, 1, sqr_dny, auy1, nau1,   &
                         auy2, nau2 )

             DO  j = 0, (ny+1)/2
                ar(j) = work(2*j)
             ENDDO
             DO  j = 1, (ny+1)/2 - 1
                ar(ny+1-j) = work(2*j+1)
             ENDDO

          ELSE

             DO  j = 0, (ny+1)/2
                work(2*j) = ar(j)
             ENDDO
             DO  j = 1, (ny+1)/2 - 1
                work(2*j+1) = ar(ny+1-j)
             ENDDO
             work(1)    = 0.0_wp
             work(ny+2) = 0.0_wp

             CALL DCRFT( 0, work, 1, work, 1, ny+1, 1, -1, sqr_dny, auy3,      &
                         nau1, auy4, nau2 )

             DO  j = 0, ny
                ar(j) = work(j)
             ENDDO

          ENDIF
#elif defined( __nec_fft )
          IF ( forward_fft )  THEN

             work(0:ny) = ar(0:ny)

             CALL DZFFT( 1, ny+1, sqr_dny, work, work, trig_yf, work2, 0 )

             DO  j = 0, (ny+1)/2
                ar(j) = work(2*j)
             ENDDO
             DO  j = 1, (ny+1)/2 - 1
                ar(ny+1-j) = work(2*j+1)
             ENDDO

          ELSE

             DO  j = 0, (ny+1)/2
                work(2*j) = ar(j)
             ENDDO
             DO  j = 1, (ny+1)/2 - 1
                work(2*j+1) = ar(ny+1-j)
             ENDDO
             work(1) = 0.0_wp
             work(ny+2) = 0.0_wp

             CALL ZDFFT( -1, ny+1, sqr_dny, work, work, trig_yb, work2, 0 )

             ar(0:ny) = work(0:ny)

          ENDIF
#endif

       ENDIF

    END SUBROUTINE fft_y_1d

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Fourier-transformation along x-direction.
!> Version for 1d domain decomposition
!> using multiple 1D FFT from Math Keisan on NEC or Temperton-algorithm
!> (no singleton-algorithm on NEC because it does not vectorize)
!------------------------------------------------------------------------------!
 
    SUBROUTINE fft_x_m( ar, direction )


       IMPLICIT NONE

       CHARACTER (LEN=*) ::  direction  !<
       
       INTEGER(iwp) ::  i     !<
       INTEGER(iwp) ::  k     !<
       INTEGER(iwp) ::  siza  !<
#if defined( __nec_fft )
       INTEGER(iwp) ::  sizw
#endif

       REAL(wp), DIMENSION(0:nx,nz)       ::  ar     !<
       REAL(wp), DIMENSION(0:nx+3,nz+1)   ::  ai     !<
       REAL(wp), DIMENSION(6*(nx+4),nz+1) ::  work1  !<
       
#if defined( __nec_fft )
       COMPLEX(wp), DIMENSION(:,:), ALLOCATABLE ::  work
#endif

       IF ( fft_method == 'temperton-algorithm' )  THEN

          siza = SIZE( ai, 1 )

          IF ( direction == 'forward')  THEN

             ai(0:nx,1:nz) = ar(0:nx,1:nz)
             ai(nx+1:,:)   = 0.0_wp

             CALL fft991cy( ai, work1, trigs_x, ifax_x, 1, siza, nx+1, nz, -1 )

             DO  k = 1, nz
                DO  i = 0, (nx+1)/2
                   ar(i,k) = ai(2*i,k)
                ENDDO
                DO  i = 1, (nx+1)/2 - 1
                   ar(nx+1-i,k) = ai(2*i+1,k)
                ENDDO
             ENDDO

          ELSE

             DO  k = 1, nz
                DO  i = 0, (nx+1)/2
                   ai(2*i,k) = ar(i,k)
                ENDDO
                DO  i = 1, (nx+1)/2 - 1
                   ai(2*i+1,k) = ar(nx+1-i,k)
                ENDDO
                ai(1,k) = 0.0_wp
                ai(nx+2,k) = 0.0_wp
             ENDDO

             CALL fft991cy( ai, work1, trigs_x, ifax_x, 1, siza, nx+1, nz, 1 )

             ar(0:nx,1:nz) = ai(0:nx,1:nz)

          ENDIF

       ELSEIF ( fft_method == 'system-specific' )  THEN

#if defined( __nec_fft )
          ALLOCATE( work((nx+4)/2+1,nz+1) )
          siza = SIZE( ai, 1 )
          sizw = SIZE( work, 1 )

          IF ( direction == 'forward')  THEN

!
!--          Tables are initialized once more. This call should not be
!--          necessary, but otherwise program aborts in asymmetric case
             CALL DZFFTM( 0, nx+1, nz1, sqr_dnx, work, nx+4, work, nx+4,       &
                          trig_xf, work1, 0 )

             ai(0:nx,1:nz) = ar(0:nx,1:nz)
             IF ( nz1 > nz )  THEN
                ai(:,nz1) = 0.0_wp
             ENDIF

             CALL DZFFTM( 1, nx+1, nz1, sqr_dnx, ai, siza, work, sizw,         &
                          trig_xf, work1, 0 )

             DO  k = 1, nz
                DO  i = 0, (nx+1)/2
                   ar(i,k) = REAL( work(i+1,k), KIND=wp )
                ENDDO
                DO  i = 1, (nx+1)/2 - 1
                   ar(nx+1-i,k) = AIMAG( work(i+1,k) )
                ENDDO
             ENDDO

          ELSE

!
!--          Tables are initialized once more. This call should not be
!--          necessary, but otherwise program aborts in asymmetric case
             CALL ZDFFTM( 0, nx+1, nz1, sqr_dnx, work, nx+4, work, nx+4,       &
                          trig_xb, work1, 0 )

             IF ( nz1 > nz )  THEN
                work(:,nz1) = 0.0_wp
             ENDIF
             DO  k = 1, nz
                work(1,k) = CMPLX( ar(0,k), 0.0_wp, KIND=wp )
                DO  i = 1, (nx+1)/2 - 1
                   work(i+1,k) = CMPLX( ar(i,k), ar(nx+1-i,k), KIND=wp )
                ENDDO
                work(((nx+1)/2)+1,k) = CMPLX( ar((nx+1)/2,k), 0.0_wp, KIND=wp )
             ENDDO

             CALL ZDFFTM( -1, nx+1, nz1, sqr_dnx, work, sizw, ai, siza, &
                          trig_xb, work1, 0 )

             ar(0:nx,1:nz) = ai(0:nx,1:nz)

          ENDIF

          DEALLOCATE( work )
#endif

       ENDIF

    END SUBROUTINE fft_x_m

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Fourier-transformation along y-direction.
!> Version for 1d domain decomposition
!> using multiple 1D FFT from Math Keisan on NEC or Temperton-algorithm
!> (no singleton-algorithm on NEC because it does not vectorize)
!------------------------------------------------------------------------------!
 
    SUBROUTINE fft_y_m( ar, ny1, direction )


       IMPLICIT NONE

       CHARACTER (LEN=*) ::  direction  !<
       
       INTEGER(iwp) ::  j     !< 
       INTEGER(iwp) ::  k     !<
       INTEGER(iwp) ::  ny1   !<
       INTEGER(iwp) ::  siza  !<
#if defined( __nec_fft )
       INTEGER(iwp) ::  sizw
#endif

       REAL(wp), DIMENSION(0:ny1,nz)      ::  ar     !<
       REAL(wp), DIMENSION(0:ny+3,nz+1)   ::  ai     !<
       REAL(wp), DIMENSION(6*(ny+4),nz+1) ::  work1  !<

#if defined( __nec_fft )
       COMPLEX(wp), DIMENSION(:,:), ALLOCATABLE ::  work
#endif


       IF ( fft_method == 'temperton-algorithm' )  THEN

          siza = SIZE( ai, 1 )

          IF ( direction == 'forward')  THEN

             ai(0:ny,1:nz) = ar(0:ny,1:nz)
             ai(ny+1:,:)   = 0.0_wp

             CALL fft991cy( ai, work1, trigs_y, ifax_y, 1, siza, ny+1, nz, -1 )

             DO  k = 1, nz
                DO  j = 0, (ny+1)/2
                   ar(j,k) = ai(2*j,k)
                ENDDO
                DO  j = 1, (ny+1)/2 - 1
                   ar(ny+1-j,k) = ai(2*j+1,k)
                ENDDO
             ENDDO

          ELSE

             DO  k = 1, nz
                DO  j = 0, (ny+1)/2
                   ai(2*j,k) = ar(j,k)
                ENDDO
                DO  j = 1, (ny+1)/2 - 1
                   ai(2*j+1,k) = ar(ny+1-j,k)
                ENDDO
                ai(1,k) = 0.0_wp
                ai(ny+2,k) = 0.0_wp
             ENDDO

             CALL fft991cy( ai, work1, trigs_y, ifax_y, 1, siza, ny+1, nz, 1 )

             ar(0:ny,1:nz) = ai(0:ny,1:nz)

          ENDIF

       ELSEIF ( fft_method == 'system-specific' )  THEN

#if defined( __nec_fft )
          ALLOCATE( work((ny+4)/2+1,nz+1) )
          siza = SIZE( ai, 1 )
          sizw = SIZE( work, 1 )

          IF ( direction == 'forward')  THEN

!
!--          Tables are initialized once more. This call should not be
!--          necessary, but otherwise program aborts in asymmetric case
             CALL DZFFTM( 0, ny+1, nz1, sqr_dny, work, ny+4, work, ny+4, &
                          trig_yf, work1, 0 )

             ai(0:ny,1:nz) = ar(0:ny,1:nz)
             IF ( nz1 > nz )  THEN
                ai(:,nz1) = 0.0_wp
             ENDIF

             CALL DZFFTM( 1, ny+1, nz1, sqr_dny, ai, siza, work, sizw, &
                          trig_yf, work1, 0 )

             DO  k = 1, nz
                DO  j = 0, (ny+1)/2
                   ar(j,k) = REAL( work(j+1,k), KIND=wp )
                ENDDO
                DO  j = 1, (ny+1)/2 - 1
                   ar(ny+1-j,k) = AIMAG( work(j+1,k) )
                ENDDO
             ENDDO

          ELSE

!
!--          Tables are initialized once more. This call should not be
!--          necessary, but otherwise program aborts in asymmetric case
             CALL ZDFFTM( 0, ny+1, nz1, sqr_dny, work, ny+4, work, ny+4, &
                          trig_yb, work1, 0 )

             IF ( nz1 > nz )  THEN
                work(:,nz1) = 0.0_wp
             ENDIF
             DO  k = 1, nz
                work(1,k) = CMPLX( ar(0,k), 0.0_wp, KIND=wp )
                DO  j = 1, (ny+1)/2 - 1
                   work(j+1,k) = CMPLX( ar(j,k), ar(ny+1-j,k), KIND=wp )
                ENDDO
                work(((ny+1)/2)+1,k) = CMPLX( ar((ny+1)/2,k), 0.0_wp, KIND=wp )
             ENDDO

             CALL ZDFFTM( -1, ny+1, nz1, sqr_dny, work, sizw, ai, siza, &
                          trig_yb, work1, 0 )

             ar(0:ny,1:nz) = ai(0:ny,1:nz)

          ENDIF

          DEALLOCATE( work )
#endif

       ENDIF

    END SUBROUTINE fft_y_m


 END MODULE fft_xy
