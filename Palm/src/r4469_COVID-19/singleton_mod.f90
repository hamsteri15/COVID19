!> @file singleton_mod.f90
!------------------------------------------------------------------------------!
! Current revisions:
! -----------------
! 
! 
! Former revisions:
! -----------------
! $Id: singleton_mod.f90 4182 2019-08-22 15:20:23Z scharf $
! Corrected "Former revisions" section
! 
! 3761 2019-02-25 15:31:42Z raasch
! statement added to prevent compiler warning about unused variables
!
! Revision 1.1  2002/05/02 18:56:59  raasch
! Initial revision
!
!
! Description:
! ------------
!> Multivariate Fast Fourier Transform
!>
!> Fortran 90 Implementation of Singleton's mixed-radix algorithm,
!> RC Singleton, Stanford Research Institute, Sept. 1968.
!>
!> Adapted from fftn.c, translated from Fortran 66 to C by Mark Olesen and
!> John Beale.
!>
!> Fourier transforms can be computed either in place, using assumed size
!> arguments, or by generic function, using assumed shape arguments.
!>
!>
!> Public:
!>
!>   fftkind                              kind parameter of complex arguments
!>                                        and function results.
!>
!>   fft(array, dim, inv, stat)           generic transform function
!>    COMPLEX(fftkind), DIMENSION(:,...,:), INTENT(IN)           :: array
!>    INTEGER,          DIMENSION(:),       INTENT(IN),  OPTIONAL:: dim
!>    LOGICAL,                              INTENT(IN),  OPTIONAL:: inv
!>    INTEGER,                              INTENT(OUT), OPTIONAL:: stat
!>
!>   fftn(array, shape, dim, inv, stat)   in place transform subroutine
!>    COMPLEX(fftkind), DIMENSION(*), INTENT(INOUT)        :: array
!>    INTEGER,          DIMENSION(:), INTENT(IN)           :: shape
!>    INTEGER,          DIMENSION(:), INTENT(IN),  OPTIONAL:: dim
!>    LOGICAL,                        INTENT(IN),  OPTIONAL:: inv
!>    INTEGER,                        INTENT(OUT), OPTIONAL:: stat
!>
!>
!> Formal Parameters:
!>
!>   array    The complex array to be transformed. array can be of arbitrary
!>            rank (i.e. up to seven).
!>
!>   shape    With subroutine fftn, the shape of the array to be transformed
!>            has to be passed separately, since fftradix - the internal trans-
!>            formation routine - will treat array always as one dimensional.
!>            The product of elements in shape must be the number of
!>            elements in array.
!>            Although passing array with assumed shape would have been nicer,
!>            I prefered assumed size in order to prevent the compiler from
!>            using a copy-in-copy-out mechanism. That would generally be
!>            necessary with fftn passing array to fftradix and with fftn
!>            being prepared for accepting non consecutive array sections.
!>            Using assumed size, it's up to the user to pass an array argu-
!>            ment, that can be addressed as continous one dimensional array
!>            without copying. Otherwise, transformation will not really be
!>            performed in place.
!>            On the other hand, since the rank of array and the size of
!>            shape needn't match, fftn is appropriate for handling more than
!>            seven dimensions.
!>            As far as function fft is concerned all this doesn't matter,
!>            because the argument will be copied anyway. Thus no extra
!>            shape argument is needed for fft.
!>
!> Optional Parameters:
!>
!>   dim      One dimensional integer array, containing the dimensions to be
!>            transformed. Default is (/1,...,N/) with N being the rank of
!>            array, i.e. complete transform. dim can restrict transformation
!>            to a subset of available dimensions. Its size must not exceed the
!>            rank of array or the size of shape respectivly.
!>
!>   inv      If .true., inverse transformation will be performed. Default is
!>            .false., i.e. forward transformation.
!>
!>   stat     If present, a system dependent nonzero status value will be
!>            returned in stat, if allocation of temporary storage failed.
!>
!>
!> Scaling:
!>
!>   Transformation results will always be scaled by the square root of the
!>   product of sizes of each dimension in dim. (See examples below)
!>
!>
!> Examples:
!>
!>   Let A be a L*M*N three dimensional complex array. Then
!>
!>     result = fft(A)
!>
!>   will produce a three dimensional transform, scaled by sqrt(L*M*N), while
!>
!>     call fftn(A, SHAPE(A))
!>
!>   will do the same in place.
!>
!>     result = fft(A, dim=(/1,3/))
!>
!>   will transform with respect to the first and the third dimension, scaled
!>   by sqrt(L*N).
!>
!>     result = fft(fft(A), inv=.true.)
!>
!>   should (approximately) reproduce A.
!>   With B having the same shape as A
!>
!>     result = fft(fft(A) * CONJG(fft(B)), inv=.true.)
!>
!>   will correlate A and B.
!>
!>
!> Remarks:
!>
!>   Following changes have been introduced with respect to fftn.c:
!>   - complex arguments and results are of type complex, rather than
!>     real an imaginary part separately.
!>   - increment parameter (magnitude of isign) has been dropped,
!>     inc is always one, direction of transform is given by inv.     
!>   - maxf and maxp have been dropped. The amount of temporary storage
!>     needed is determined by the fftradix routine. Both fftn and fft
!>     can handle any size of array. (Maybe they take a lot of time and
!>     memory, but they will do it)
!>
!>   Redesigning fftradix in a way, that it handles assumed shape arrays
!>   would have been desirable. However, I found it rather hard to do this
!>   in an efficient way. Problems were:
!>   - to prevent stride multiplications when indexing arrays. At least our
!>     compiler was not clever enough to discover that in fact additions
!>     would do the job as well. On the other hand, I haven't been clever
!>     enough to find an implementation using array operations.
!>   - fftradix is rather large and different versions would be necessaray
!>     for each possible rank of array.
!>   Consequently, in place transformation still needs the argument stored
!>   in a consecutive bunch of memory and can't be performed on array
!>   sections like A(100:199:-3, 50:1020). Calling fftn with such sections
!>   will most probably imply copy-in-copy-out. However, the function fft
!>   works with everything it gets and should be convenient to use.
!>
!> Michael Steffens, 09.12.96, <Michael.Steffens@mbox.muk.uni-hannover.de>
!> Restructured fftradix for better optimization. M. Steffens, 4 June 1997
!------------------------------------------------------------------------------!
 MODULE singleton
 

    USE kinds

    IMPLICIT NONE

    PRIVATE
    PUBLIC:: fft, fftn

    REAL(wp), PARAMETER:: sin60 = 0.86602540378443865_wp
    REAL(wp), PARAMETER:: cos72 = 0.30901699437494742_wp
    REAL(wp), PARAMETER:: sin72 = 0.95105651629515357_wp
    REAL(wp), PARAMETER:: pi    = 3.14159265358979323_wp

    INTERFACE fft
       MODULE PROCEDURE fft1d
       MODULE PROCEDURE fft2d
       MODULE PROCEDURE fft3d
       MODULE PROCEDURE fft4d
       MODULE PROCEDURE fft5d
       MODULE PROCEDURE fft6d
       MODULE PROCEDURE fft7d
    END INTERFACE


 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing function description.
!------------------------------------------------------------------------------!
 FUNCTION fft1d(array, dim, inv, stat) RESULT(ft)
!
!-- Formal parameters
    COMPLEX(wp),  DIMENSION(:), INTENT(IN)           :: array
    INTEGER(iwp), DIMENSION(:), INTENT(IN),  OPTIONAL:: dim
    INTEGER(iwp),               INTENT(OUT), OPTIONAL:: stat
    LOGICAL,                    INTENT(IN),  OPTIONAL:: inv
!
!-- Function result
    COMPLEX(wp), DIMENSION(SIZE(array, 1)):: ft

    INTEGER(iwp)::  idum
    INTEGER(iwp)::  ishape(1)

!
!-- Intrinsics used
    INTRINSIC SIZE, SHAPE

    ft = array
    ishape = SHAPE( array )
    CALL fftn(ft, ishape, inv = inv,  stat = stat)
!
!-- Next statement to prevent compiler warning about unused variable
    IF ( PRESENT( dim ) )  idum = 1

 END FUNCTION fft1d


!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing function description.
!------------------------------------------------------------------------------!
 FUNCTION fft2d(array, dim, inv, stat) RESULT(ft)
!
!-- Formal parameters
    COMPLEX(wp),  DIMENSION(:,:), INTENT(IN)           :: array
    INTEGER(iwp), DIMENSION(:),   INTENT(IN),  OPTIONAL:: dim
    INTEGER(iwp),                 INTENT(OUT), OPTIONAL:: stat
    LOGICAL,                      INTENT(IN),  OPTIONAL:: inv
!
!-- Function result
    COMPLEX(wp), DIMENSION(SIZE(array, 1), SIZE(array, 2)):: ft

    INTEGER(iwp) ::  ishape(2)
!
!-- Intrinsics used
    INTRINSIC SIZE, SHAPE

    ft = array
    ishape = SHAPE( array )
    CALL fftn(ft, ishape, dim, inv, stat)

 END FUNCTION fft2d


!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing function description.
!------------------------------------------------------------------------------!
 FUNCTION fft3d(array, dim, inv, stat) RESULT(ft)
!
!-- Formal parameters
    COMPLEX(wp),  DIMENSION(:,:,:), INTENT(IN)           :: array
    INTEGER(iwp), DIMENSION(:),     INTENT(IN),  OPTIONAL:: dim
    INTEGER(iwp),                   INTENT(OUT), OPTIONAL:: stat
    LOGICAL,                        INTENT(IN),  OPTIONAL:: inv
!
!-- Function result
    COMPLEX(wp), &
         DIMENSION(SIZE(array, 1), SIZE(array, 2), SIZE(array, 3)):: ft

    INTEGER(iwp) ::  ishape(3)

!
!-- Intrinsics used
    INTRINSIC SIZE, SHAPE

    ft = array
    ishape = SHAPE( array)
    CALL fftn(ft, ishape, dim, inv, stat)

 END FUNCTION fft3d


!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing function description.
!------------------------------------------------------------------------------!
 FUNCTION fft4d(array, dim, inv, stat) RESULT(ft)
!
!-- Formal parameters
    COMPLEX(wp),  DIMENSION(:,:,:,:), INTENT(IN)           :: array
    INTEGER(iwp), DIMENSION(:),       INTENT(IN),  OPTIONAL:: dim
    INTEGER(iwp),                     INTENT(OUT), OPTIONAL:: stat
    LOGICAL,                          INTENT(IN),  OPTIONAL:: inv
!
!-- Function result
    COMPLEX(wp), DIMENSION( &
         SIZE(array, 1), SIZE(array, 2), SIZE(array, 3), SIZE(array, 4)):: ft

    INTEGER(iwp) ::  ishape(4)
!
!-- Intrinsics used
    INTRINSIC SIZE, SHAPE

    ft = array
    ishape = SHAPE( array )
    CALL fftn(ft, ishape, dim, inv, stat)

 END FUNCTION fft4d


!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing function description.
!------------------------------------------------------------------------------!
 FUNCTION fft5d(array, dim, inv, stat) RESULT(ft)
!
!-- Formal parameters
    COMPLEX(wp),  DIMENSION(:,:,:,:,:), INTENT(IN)           :: array
    INTEGER(iwp), DIMENSION(:),         INTENT(IN),  OPTIONAL:: dim
    INTEGER(iwp),                       INTENT(OUT), OPTIONAL:: stat
    LOGICAL,                            INTENT(IN),  OPTIONAL:: inv
!
!-- Function result
    COMPLEX(wp), DIMENSION( &
         SIZE(array, 1), SIZE(array, 2), SIZE(array, 3), SIZE(array, 4), &
         SIZE(array, 5)):: ft

    INTEGER(iwp) ::  ishape(5)

!
!-- Intrinsics used
    INTRINSIC SIZE, SHAPE

    ft = array
    ishape = SHAPE( array )
    CALL fftn(ft, ishape, dim, inv, stat)

 END FUNCTION fft5d


!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing function description.
!------------------------------------------------------------------------------!
 FUNCTION fft6d(array, dim, inv, stat) RESULT(ft)
!
!-- Formal parameters
    COMPLEX(wp),  DIMENSION(:,:,:,:,:,:), INTENT(IN)           :: array
    INTEGER(iwp), DIMENSION(:),           INTENT(IN),  OPTIONAL:: dim
    INTEGER(iwp),                         INTENT(OUT), OPTIONAL:: stat
    LOGICAL,                              INTENT(IN),  OPTIONAL:: inv
!
!-- Function result
    COMPLEX(wp), DIMENSION( &
         SIZE(array, 1), SIZE(array, 2), SIZE(array, 3), SIZE(array, 4), &
         SIZE(array, 5), SIZE(array, 6)):: ft

    INTEGER(iwp) ::  ishape(6)

!
!-- Intrinsics used
    INTRINSIC SIZE, SHAPE

    ft = array
    ishape = SHAPE( array )
    CALL fftn(ft, ishape, dim, inv, stat)

 END FUNCTION fft6d


!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing function description.
!------------------------------------------------------------------------------!
 FUNCTION fft7d(array, dim, inv, stat) RESULT(ft)
!
!-- Formal parameters
    COMPLEX(wp), DIMENSION(:,:,:,:,:,:,:), INTENT(IN)           :: array
    INTEGER(iwp),          DIMENSION(:),   INTENT(IN),  OPTIONAL:: dim
    INTEGER(iwp),                          INTENT(OUT), OPTIONAL:: stat
    LOGICAL,                               INTENT(IN),  OPTIONAL:: inv
!
!-- Function result
    COMPLEX(wp), DIMENSION( &
         SIZE(array, 1), SIZE(array, 2), SIZE(array, 3), SIZE(array, 4), &
         SIZE(array, 5), SIZE(array, 6), SIZE(array, 7)):: ft

    INTEGER(iwp) ::  ishape(7)

!
!-- Intrinsics used
    INTRINSIC SIZE, SHAPE

    ft = array
    ishape = SHAPE( array )
    CALL fftn(ft, ishape, dim, inv, stat)

 END FUNCTION fft7d


!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
 SUBROUTINE fftn(array, shape, dim, inv, stat)
!
!-- Formal parameters
    COMPLEX(wp),  DIMENSION(*), INTENT(INOUT)        :: array
    INTEGER(iwp), DIMENSION(:), INTENT(IN)           :: shape
    INTEGER(iwp), DIMENSION(:), INTENT(IN),  OPTIONAL:: dim
    INTEGER(iwp),               INTENT(OUT), OPTIONAL:: stat
    LOGICAL,                    INTENT(IN),  OPTIONAL:: inv
!
!-- Local arrays
    INTEGER(iwp), DIMENSION(SIZE(shape)):: d
!
!-- Local scalars
    LOGICAL      :: inverse
    INTEGER(iwp) :: i, ndim, ntotal
    REAL(wp):: scale
!
!-- Intrinsics used
    INTRINSIC PRESENT, MIN, PRODUCT, SIZE, SQRT

!
!-- Optional parameter settings
    IF (PRESENT(inv)) THEN
       inverse = inv
    ELSE
       inverse = .FALSE.
    END IF
    IF (PRESENT(dim)) THEN
       ndim = MIN(SIZE(dim), SIZE(d))
       d(1:ndim) = DIM(1:ndim)
    ELSE
       ndim = SIZE(d)
       d = (/(i, i = 1, SIZE(d))/)
    END IF

    ntotal = PRODUCT(shape)
    scale = SQRT(1.0_wp / PRODUCT(shape(d(1:ndim))))
    DO i = 1, ntotal
       array(i) = CMPLX(REAL(array(i)) * scale, AIMAG(array(i)) * scale, &
            KIND=wp)
    END DO

    DO i = 1, ndim
       CALL fftradix(array, ntotal, shape(d(i)), PRODUCT(shape(1:d(i))), &
            inverse, stat)
       IF (PRESENT(stat)) THEN
          IF (stat /=0) RETURN
       END IF
    END DO

 END SUBROUTINE fftn


!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
 SUBROUTINE fftradix(array, ntotal, npass, nspan, inv, stat)
!
!-- Formal parameters
    COMPLEX(wp),  DIMENSION(*), INTENT(INOUT)        :: array
    INTEGER(iwp),               INTENT(IN)           :: ntotal, npass, nspan
    INTEGER(iwp),               INTENT(OUT), OPTIONAL:: stat
    LOGICAL,                    INTENT(IN)           :: inv
!
!-- Local arrays
    COMPLEX(wp),  DIMENSION(:), ALLOCATABLE  :: ctmp
    INTEGER(iwp), DIMENSION(BIT_SIZE(0))     :: factor
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE  :: perm
    REAL(wp),     DIMENSION(:), ALLOCATABLE  :: sine, cosine
!
!-- Local scalars
    INTEGER(iwp)         :: maxfactor, nfactor, nsquare, nperm
!
!-- Intrinsics used
    INTRINSIC MAXVAL, MOD, PRESENT, ISHFT, BIT_SIZE, SIN, COS, &
         CMPLX, REAL, AIMAG

    IF (npass <= 1) RETURN

    CALL factorize(npass, factor, nfactor, nsquare)

    maxfactor = MAXVAL(factor(:nfactor))
    IF (nfactor - ISHFT(nsquare, 1) > 0) THEN
       nperm = MAX(nfactor + 1, PRODUCT(factor(nsquare+1: nfactor-nsquare)) - 1)
    ELSE
       nperm = nfactor + 1
    END IF

    IF (PRESENT(stat)) THEN
       ALLOCATE(ctmp(maxfactor), sine(maxfactor), cosine(maxfactor), STAT=stat)
       IF (stat /= 0) RETURN
       CALL transform(array, ntotal, npass, nspan, &
            factor, nfactor, ctmp, sine, cosine, inv)
       DEALLOCATE(sine, cosine, STAT=stat)
       IF (stat /= 0) RETURN
       ALLOCATE(perm(nperm), STAT=stat)
       IF (stat /= 0) RETURN
       CALL permute(array, ntotal, npass, nspan, &
            factor, nfactor, nsquare, maxfactor, &
            ctmp, perm)
       DEALLOCATE(perm, ctmp, STAT=stat)
       IF (stat /= 0) RETURN
    ELSE
       ALLOCATE(ctmp(maxfactor), sine(maxfactor), cosine(maxfactor))
       CALL transform(array, ntotal, npass, nspan, &
            factor, nfactor, ctmp, sine, cosine, inv)
       DEALLOCATE(sine, cosine)
       ALLOCATE(perm(nperm))
       CALL permute(array, ntotal, npass, nspan, &
            factor, nfactor, nsquare, maxfactor, &
            ctmp, perm)
       DEALLOCATE(perm, ctmp)
    END IF


  CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE factorize(npass, factor, nfactor, nsquare)
!
!--   Formal parameters
      INTEGER(iwp),               INTENT(IN) :: npass
      INTEGER(iwp), DIMENSION(*), INTENT(OUT):: factor
      INTEGER(iwp),               INTENT(OUT):: nfactor, nsquare
!
!--   Local scalars
      INTEGER(iwp):: j, jj, k

      nfactor = 0
      k = npass
      DO WHILE (MOD(k, 16) == 0) 
         nfactor = nfactor + 1
         factor(nfactor) = 4
         k = k / 16
      END DO
      j = 3
      jj = 9
      DO
         DO WHILE (MOD(k, jj) == 0)
            nfactor = nfactor + 1
            factor(nfactor) = j
            k = k / jj
         END DO
         j = j + 2
         jj = j * j
         IF (jj > k) EXIT
      END DO
      IF (k <= 4) THEN
         nsquare = nfactor
         factor(nfactor + 1) = k
         IF (k /= 1) nfactor = nfactor + 1
      ELSE 
         IF (k - ISHFT(k / 4, 2) == 0) THEN
            nfactor = nfactor + 1
            factor(nfactor) = 2
            k = k / 4
         END IF
         nsquare = nfactor
         j = 2
         DO
            IF (MOD(k, j) == 0) THEN
               nfactor = nfactor + 1
               factor(nfactor) = j
               k = k / j
            END IF
            j = ISHFT((j + 1) / 2, 1) + 1
            IF (j > k) EXIT
         END DO
      END IF
      IF (nsquare > 0) THEN
         j = nsquare
         DO
            nfactor = nfactor + 1
            factor(nfactor) = factor(j)
            j = j - 1
            IF (j==0) EXIT
         END DO
      END IF

    END SUBROUTINE factorize


!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE transform(array, ntotal, npass, nspan, &
         factor, nfactor, ctmp, sine, cosine, inv) !-- compute fourier transform
!
!--   Formal parameters
      COMPLEX(wp),  DIMENSION(*), INTENT(IN OUT):: array
      COMPLEX(wp),  DIMENSION(*), INTENT(OUT)   :: ctmp
      INTEGER(iwp),               INTENT(IN)    :: ntotal, npass, nspan
      INTEGER(iwp), DIMENSION(*), INTENT(IN)    :: factor
      INTEGER(iwp),               INTENT(IN)    :: nfactor
      LOGICAL,                    INTENT(IN)    :: inv
      REAL(wp),     DIMENSION(*), INTENT(OUT)   :: sine, cosine
!
!--   Local scalars
      INTEGER(iwp):: ii, ispan
      INTEGER(iwp):: j, jc, jf, jj
      INTEGER(iwp):: k, kk, kspan, k1, k2, k3, k4
      INTEGER(iwp):: nn, nt
      REAL(wp)    :: s60, c72, s72, pi2, radf
      REAL(wp)    :: c1, s1, c2, s2, c3, s3, cd, sd, ak
      COMPLEX(wp) :: cc, cj, ck, cjp, cjm, ckp, ckm

      c72 = cos72
      IF (inv) THEN
         s72 = sin72
         s60 = sin60
         pi2 = pi
      ELSE
         s72 = -sin72
         s60 = -sin60
         pi2 = -pi
      END IF

      nt = ntotal
      nn = nt - 1
      kspan = nspan
      jc = nspan / npass
      radf = pi2 * jc
      pi2 = pi2 * 2.0_wp !-- use 2 PI from here on

      ii = 0
      jf = 0
      DO
         sd = radf / kspan
         cd = SIN(sd)
         cd = 2.0_wp * cd * cd
         sd = SIN(sd + sd)
         kk = 1
         ii = ii + 1

         SELECT CASE (factor(ii))
         CASE (2)
!
!--         Transform for factor of 2 (including rotation factor)
            kspan = kspan / 2
            k1 = kspan + 2
            DO
               DO
                  k2 = kk + kspan
                  ck = array(k2)
                  array(k2) = array(kk)-ck
                  array(kk) = array(kk) + ck
                  kk = k2 + kspan
                  IF (kk > nn) EXIT
               END DO
               kk = kk - nn
               IF (kk > jc) EXIT
            END DO
            IF (kk > kspan) RETURN
            DO
               c1 = 1.0_wp - cd
               s1 = sd
               DO
                  DO
                     DO
                        k2 = kk + kspan
                        ck = array(kk) - array(k2)
                        array(kk) = array(kk) + array(k2)
                        array(k2) = ck * CMPLX(c1, s1, KIND=wp)
                        kk = k2 + kspan
                        IF (kk >= nt) EXIT
                     END DO
                     k2 = kk - nt
                     c1 = -c1
                     kk = k1 - k2
                     IF (kk <= k2) EXIT
                  END DO
                  ak = c1 - (cd * c1 + sd * s1)
                  s1 = sd * c1 - cd * s1 + s1
                  c1 = 2.0_wp - (ak * ak + s1 * s1)
                  s1 = s1 * c1
                  c1 = c1 * ak
                  kk = kk + jc
                  IF (kk >= k2) EXIT
               END DO
               k1 = k1 + 1 + 1
               kk = (k1 - kspan) / 2 + jc
               IF (kk > jc + jc) EXIT
            END DO

         CASE (4) !-- transform for factor of 4
            ispan = kspan
            kspan = kspan / 4

            DO
               c1 = 1.0_wp
               s1 = 0.0_wp
               DO
                  DO
                     k1 = kk + kspan
                     k2 = k1 + kspan
                     k3 = k2 + kspan
                     ckp = array(kk) + array(k2)
                     ckm = array(kk) - array(k2)
                     cjp = array(k1) + array(k3)
                     cjm = array(k1) - array(k3)
                     array(kk) = ckp + cjp
                     cjp = ckp - cjp
                     IF (inv) THEN
                        ckp = ckm + CMPLX(-AIMAG(cjm), REAL(cjm), KIND=wp)
                        ckm = ckm + CMPLX(AIMAG(cjm), -REAL(cjm), KIND=wp)
                     ELSE
                        ckp = ckm + CMPLX(AIMAG(cjm), -REAL(cjm), KIND=wp)
                        ckm = ckm + CMPLX(-AIMAG(cjm), REAL(cjm), KIND=wp)
                     END IF
!
!--                  Avoid useless multiplies
                     IF (s1 == 0.0_wp) THEN
                        array(k1) = ckp
                        array(k2) = cjp
                        array(k3) = ckm
                     ELSE
                        array(k1) = ckp * CMPLX(c1, s1, KIND=wp)
                        array(k2) = cjp * CMPLX(c2, s2, KIND=wp)
                        array(k3) = ckm * CMPLX(c3, s3, KIND=wp)
                     END IF
                     kk = k3 + kspan
                     IF (kk > nt) EXIT
                  END DO

                  c2 = c1 - (cd * c1 + sd * s1)
                  s1 = sd * c1 - cd * s1 + s1
                  c1 = 2.0_wp - (c2 * c2 + s1 * s1)
                  s1 = s1 * c1
                  c1 = c1 * c2
!
!--               Values of c2, c3, s2, s3 that will get used next time
                  c2 = c1 * c1 - s1 * s1
                  s2 = 2.0_wp * c1 * s1
                  c3 = c2 * c1 - s2 * s1
                  s3 = c2 * s1 + s2 * c1
                  kk = kk - nt + jc
                  IF (kk > kspan) EXIT
               END DO
               kk = kk - kspan + 1
               IF (kk > jc) EXIT
            END DO
            IF (kspan == jc) RETURN

         CASE default
!
!--         Transform for odd factors
            k = factor(ii)
            ispan = kspan
            kspan = kspan / k

            SELECT CASE (k)
            CASE (3) !-- transform for factor of 3 (optional code)
               DO
                  DO
                     k1 = kk + kspan
                     k2 = k1 + kspan
                     ck = array(kk)
                     cj = array(k1) + array(k2)
                     array(kk) = ck + cj
                     ck = ck - CMPLX( &
                          0.5_wp * REAL (cj), &
                          0.5_wp * AIMAG(cj), &
                          KIND=wp)
                     cj = CMPLX( &
                          (REAL (array(k1)) - REAL (array(k2))) * s60, &
                          (AIMAG(array(k1)) - AIMAG(array(k2))) * s60, &
                          KIND=wp)
                     array(k1) = ck + CMPLX(-AIMAG(cj), REAL(cj), KIND=wp)
                     array(k2) = ck + CMPLX(AIMAG(cj), -REAL(cj), KIND=wp)
                     kk = k2 + kspan
                     IF (kk >= nn) EXIT
                  END DO
                  kk = kk - nn
                  IF (kk > kspan) EXIT
               END DO

            CASE (5) !-- transform for factor of 5 (optional code)
               c2 = c72 * c72 - s72 * s72
               s2 = 2.0_wp * c72 * s72
               DO
                  DO
                     k1 = kk + kspan
                     k2 = k1 + kspan
                     k3 = k2 + kspan
                     k4 = k3 + kspan
                     ckp = array(k1) + array(k4)
                     ckm = array(k1) - array(k4)
                     cjp = array(k2) + array(k3)
                     cjm = array(k2) - array(k3)
                     cc = array(kk)
                     array(kk) = cc + ckp + cjp
                     ck = CMPLX(REAL(ckp) * c72, AIMAG(ckp) * c72, &
                          KIND=wp) + &
                          CMPLX(REAL(cjp) * c2,  AIMAG(cjp) * c2,  &
                          KIND=wp) + cc
                     cj = CMPLX(REAL(ckm) * s72, AIMAG(ckm) * s72, &
                          KIND=wp) + &
                          CMPLX(REAL(cjm) * s2,  AIMAG(cjm) * s2,  &
                          KIND=wp)
                     array(k1) = ck + CMPLX(-AIMAG(cj), REAL(cj), KIND=wp)
                     array(k4) = ck + CMPLX(AIMAG(cj), -REAL(cj), KIND=wp)
                     ck = CMPLX(REAL(ckp) * c2,  AIMAG(ckp) * c2,  &
                          KIND=wp) + &
                          CMPLX(REAL(cjp) * c72, AIMAG(cjp) * c72, &
                          KIND=wp) + cc
                     cj = CMPLX(REAL(ckm) * s2,  AIMAG(ckm) * s2,  &
                          KIND=wp) - &
                          CMPLX(REAL(cjm) * s72, AIMAG(cjm) * s72, &
                          KIND=wp)
                     array(k2) = ck + CMPLX(-AIMAG(cj), REAL(cj), KIND=wp)
                     array(k3) = ck + CMPLX(AIMAG(cj), -REAL(cj), KIND=wp)
                     kk = k4 + kspan
                     IF (kk >= nn) EXIT
                  END DO
                  kk = kk - nn
                  IF (kk > kspan) EXIT
               END DO

            CASE default
               IF (k /= jf) THEN
                  jf = k
                  s1 = pi2 / k
                  c1 = COS(s1)
                  s1 = SIN(s1)
                  cosine (jf) = 1.0_wp
                  sine (jf) = 0.0_wp
                  j = 1
                  DO
                     cosine (j) = cosine (k) * c1 + sine (k) * s1
                     sine (j) = cosine (k) * s1 - sine (k) * c1
                     k = k-1
                     cosine (k) = cosine (j)
                     sine (k) = -sine (j)
                     j = j + 1
                     IF (j >= k) EXIT
                  END DO
               END IF
               DO
                  DO
                     k1 = kk
                     k2 = kk + ispan
                     cc = array(kk)
                     ck = cc
                     j = 1
                     k1 = k1 + kspan
                     DO
                        k2 = k2 - kspan
                        j = j + 1
                        ctmp(j) = array(k1) + array(k2)
                        ck = ck + ctmp(j)
                        j = j + 1
                        ctmp(j) = array(k1) - array(k2)
                        k1 = k1 + kspan
                        IF (k1 >= k2) EXIT
                     END DO
                     array(kk) = ck
                     k1 = kk
                     k2 = kk + ispan
                     j = 1
                     DO
                        k1 = k1 + kspan
                        k2 = k2 - kspan
                        jj = j
                        ck = cc
                        cj = (0.0_wp, 0.0_wp)
                        k = 1
                        DO
                           k = k + 1
                           ck = ck + CMPLX( &
                                REAL (ctmp(k)) * cosine(jj), &
                                AIMAG(ctmp(k)) * cosine(jj), KIND=wp)
                           k = k + 1
                           cj = cj + CMPLX( &
                                REAL (ctmp(k)) * sine(jj), &
                                AIMAG(ctmp(k)) * sine(jj), KIND=wp)
                           jj = jj + j
                           IF (jj > jf) jj = jj - jf
                           IF (k >= jf) EXIT
                        END DO
                        k = jf - j
                        array(k1) = ck + CMPLX(-AIMAG(cj), REAL(cj), &
                             KIND=wp)
                        array(k2) = ck + CMPLX(AIMAG(cj), -REAL(cj), &
                             KIND=wp)
                        j = j + 1
                        IF (j >= k) EXIT
                     END DO
                     kk = kk + ispan
                     IF (kk > nn) EXIT
                  END DO
                  kk = kk - nn
                  IF (kk > kspan) EXIT
               END DO

            END SELECT
!
!--         Multiply by rotation factor (except for factors of 2 and 4)
            IF (ii == nfactor) RETURN
            kk = jc + 1
            DO
               c2 = 1.0_wp - cd
               s1 = sd
               DO
                  c1 = c2
                  s2 = s1
                  kk = kk + kspan
                  DO
                     DO
                        array(kk) = CMPLX(c2, s2, KIND=wp) * array(kk)
                        kk = kk + ispan
                        IF (kk > nt) EXIT
                     END DO
                     ak = s1 * s2
                     s2 = s1 * c2 + c1 * s2
                     c2 = c1 * c2 - ak
                     kk = kk - nt + kspan
                     IF (kk > ispan) EXIT
                  END DO
                  c2 = c1 - (cd * c1 + sd * s1)
                  s1 = s1 + sd * c1 - cd * s1
                  c1 = 2.0_wp - (c2 * c2 + s1 * s1)
                  s1 = s1 * c1
                  c2 = c2 * c1
                  kk = kk - ispan + jc
                  IF (kk > kspan) EXIT
               END DO
               kk = kk - kspan + jc + 1
               IF (kk > jc + jc) EXIT
            END DO

         END SELECT
      END DO
    END SUBROUTINE transform


!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE permute(array, ntotal, npass, nspan, &
         factor, nfactor, nsquare, maxfactor, &
         ctmp, perm)
!
!--   Formal parameters
      COMPLEX(wp),  DIMENSION(*), INTENT(IN OUT):: array
      COMPLEX(wp),  DIMENSION(*), INTENT(OUT)   :: ctmp
      INTEGER(iwp),               INTENT(IN)    :: ntotal, npass, nspan
      INTEGER(iwp), DIMENSION(*), INTENT(IN OUT):: factor
      INTEGER(iwp),               INTENT(IN)    :: nfactor, nsquare
      INTEGER(iwp),               INTENT(IN)    :: maxfactor
      INTEGER(iwp), DIMENSION(*), INTENT(OUT)   :: perm
!
!--   Local scalars
      COMPLEX(wp) :: ck
      INTEGER(iwp):: ii, ispan
      INTEGER(iwp):: j, jc, jj
      INTEGER(iwp):: k, kk, kspan, kt, k1, k2, k3
      INTEGER(iwp):: nn, nt
!
!--   Permute the results to normal order---done in two stages
!--   Permutation for square factors of n

      nt = ntotal
      nn = nt - 1
      kt = nsquare
      kspan = nspan
      jc = nspan / npass

      perm (1) = nspan
      IF (kt > 0) THEN
         k = kt + kt + 1
         IF (nfactor < k) k = k - 1
         j = 1
         perm (k + 1) = jc
         DO
            perm (j + 1) = perm (j) / factor(j)
            perm (k) = perm (k + 1) * factor(j)
            j = j + 1
            k = k - 1
            IF (j >= k) EXIT
         END DO
         k3 = perm (k + 1)
         kspan = perm (2)
         kk = jc + 1
         k2 = kspan + 1
         j = 1

         IF (npass /= ntotal) THEN
            permute_multi: DO
               DO
                  DO
                     k = kk + jc
                     DO
!
!--                     Swap array(kk) <> array(k2)
                        ck = array(kk)
                        array(kk) = array(k2)
                        array(k2) = ck
                        kk = kk + 1
                        k2 = k2 + 1
                        IF (kk >= k) EXIT
                     END DO
                     kk = kk + nspan - jc
                     k2 = k2 + nspan - jc
                     IF (kk >= nt) EXIT
                  END DO
                  kk = kk - nt + jc
                  k2 = k2 - nt + kspan
                  IF (k2 >= nspan) EXIT
               END DO
               DO
                  DO
                     k2 = k2 - perm (j)
                     j = j + 1
                     k2 = perm (j + 1) + k2
                     IF (k2 <= perm (j)) EXIT
                  END DO
                  j = 1
                  DO
                     IF (kk < k2) CYCLE permute_multi
                     kk = kk + jc
                     k2 = k2 + kspan
                     IF (k2 >= nspan) EXIT
                  END DO
                  IF (kk >= nspan) EXIT
               END DO
               EXIT
            END DO permute_multi
         ELSE
            permute_single: DO
               DO
!
!--               Swap array(kk) <> array(k2)
                  ck = array(kk)
                  array(kk) = array(k2)
                  array(k2) = ck
                  kk = kk + 1
                  k2 = k2 + kspan
                  IF (k2 >= nspan) EXIT
               END DO
               DO
                  DO
                     k2 = k2 - perm (j)
                     j = j + 1
                     k2 = perm (j + 1) + k2
                     IF (k2 <= perm (j)) EXIT
                  END DO
                  j = 1
                  DO
                     IF (kk < k2) CYCLE permute_single
                     kk = kk + 1
                     k2 = k2 + kspan
                     IF (k2 >= nspan) EXIT
                  END DO
                  IF (kk >= nspan) EXIT
               END DO
               EXIT
            END DO permute_single
         END IF
         jc = k3
      END IF

      IF (ISHFT(kt, 1) + 1 >= nfactor) RETURN

      ispan = perm (kt + 1)
!
!--   Permutation for square-free factors of n
      j = nfactor - kt
      factor(j + 1) = 1
      DO
         factor(j) = factor(j) * factor(j+1)
         j = j - 1
         IF (j == kt) EXIT
      END DO
      kt = kt + 1
      nn = factor(kt) - 1
      j = 0
      jj = 0
      DO
         k = kt + 1
         k2 = factor(kt)
         kk = factor(k)
         j = j + 1
         IF (j > nn) EXIT !-- exit infinite loop
         jj = jj + kk
         DO WHILE (jj >= k2)
            jj = jj - k2
            k2 = kk
            k = k + 1
            kk = factor(k)
            jj = jj + kk
         END DO
         perm (j) = jj
      END DO
!
!--   Determine the permutation cycles of length greater than 1
      j = 0
      DO
         DO
            j = j + 1
            kk = perm(j)
            IF (kk >= 0) EXIT
         END DO
         IF (kk /= j) THEN
            DO
               k = kk
               kk = perm (k)
               perm (k) = -kk
               IF (kk == j) EXIT
            END DO
            k3 = kk
         ELSE
            perm (j) = -j
            IF (j == nn) EXIT !-- exit infinite loop
         END IF
      END DO
!
!--   Reorder a and b, following the permutation cycles
      DO
         j = k3 + 1
         nt = nt - ispan
         ii = nt - 1 + 1
         IF (nt < 0) EXIT !-- exit infinite loop
         DO
            DO
               j = j-1
               IF (perm(j) >= 0) EXIT
            END DO
            jj = jc
            DO
               kspan = jj
               IF (jj > maxfactor) kspan = maxfactor
               jj = jj - kspan
               k = perm(j)
               kk = jc * k + ii + jj
               k1 = kk + kspan
               k2 = 0
               DO
                  k2 = k2 + 1
                  ctmp(k2) = array(k1)
                  k1 = k1 - 1
                  IF (k1 == kk) EXIT
               END DO
               DO
                  k1 = kk + kspan
                  k2 = k1 - jc * (k + perm(k))
                  k = -perm(k)
                  DO
                     array(k1) = array(k2)
                     k1 = k1 - 1
                     k2 = k2 - 1
                     IF (k1 == kk) EXIT
                  END DO
                  kk = k2
                  IF (k == j) EXIT
               END DO
               k1 = kk + kspan
               k2 = 0
               DO
                  k2 = k2 + 1
                  array(k1) = ctmp(k2)
                  k1 = k1 - 1
                  IF (k1 == kk) EXIT
               END DO
               IF (jj == 0) EXIT
            END DO
            IF (j == 1) EXIT
         END DO
      END DO

    END SUBROUTINE permute

 END SUBROUTINE fftradix

 END MODULE singleton
