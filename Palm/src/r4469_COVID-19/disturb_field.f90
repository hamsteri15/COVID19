!> @file disturb_field.f90
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
! $Id: disturb_field.f90 4457 2020-03-11 14:20:43Z raasch $
! use statement for exchange horiz added
! 
! 4360 2020-01-07 11:25:50Z suehring
! Introduction of wall_flags_total_0, which currently sets bits based on static
! topography information used in wall_flags_static_0
! 
! 4329 2019-12-10 15:46:36Z motisi
! Renamed wall_flags_0 to wall_flags_static_0
! 
! 4237 2019-09-25 11:33:42Z knoop
! Added missing OpenMP directives
! 
! 4182 2019-08-22 15:20:23Z scharf
! Corrected "Former revisions" section
! 
! 3849 2019-04-01 16:35:16Z knoop
! Corrected "Former revisions" section
!
! Revision 1.1  1998/02/04 15:40:45  raasch
! Initial revision
!
!
! Description:
! ------------
!> Imposing a random perturbation on a 3D-array.
!> On parallel computers, the random number generator is as well called for all
!> gridpoints of the total domain to ensure, regardless of the number of PEs
!> used, that the elements of the array have the same values in the same
!> order in every case. The perturbation range is steered by dist_range.
!------------------------------------------------------------------------------!
 SUBROUTINE disturb_field( var_char, dist1, field )
 

    USE control_parameters,                                                    &
        ONLY:  dist_nxl, dist_nxr, dist_nyn, dist_nys, dist_range,             &
               disturbance_amplitude, disturbance_created,                     &
               disturbance_level_ind_b, disturbance_level_ind_t, iran,         &
               random_generator, topography
                
    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point
        
    USE exchange_horiz_mod,                                                    &
        ONLY:  exchange_horiz

    USE indices,                                                               &
        ONLY:  nbgp, nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg, nzb, nzb_max, &
               nzt, wall_flags_total_0
        
    USE kinds
    
    USE random_function_mod,                                                   &
        ONLY: random_function
        
    USE random_generator_parallel,                                             &
        ONLY:  random_number_parallel, random_seed_parallel, random_dummy,     &
               seq_random_array

    IMPLICIT NONE

    CHARACTER (LEN = *) ::  var_char !< flag to distinguish betwenn u- and v-component

    INTEGER(iwp) ::  flag_nr !< number of respective flag for u- or v-grid 
    INTEGER(iwp) ::  i       !< index variable
    INTEGER(iwp) ::  j       !< index variable
    INTEGER(iwp) ::  k       !< index variable

    REAL(wp) ::  randomnumber  !<
    
    REAL(wp) ::  dist1(nzb:nzt+1,nysg:nyng,nxlg:nxrg)  !<
    REAL(wp) ::  field(nzb:nzt+1,nysg:nyng,nxlg:nxrg)  !<
    
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  dist2  !<


    CALL cpu_log( log_point(20), 'disturb_field', 'start' )
!
!-- Set flag number, 20 for u-grid, 21 for v-grid, required to mask topography
    flag_nr = MERGE( 20, 21, TRIM(var_char) == 'u' )
!
!-- Create an additional temporary array and initialize the arrays needed
!-- to store the disturbance
    ALLOCATE( dist2(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
!$ACC DATA CREATE(dist2(nzb:nzt+1,nysg:nyng,nxlg:nxrg))

!
!-- dist1 is initialized on the host (see below) and then updated on the device.
    dist1 = 0.0_wp
    !$ACC KERNELS PRESENT(dist2)
    dist2 = 0.0_wp
    !$ACC END KERNELS

!
!-- Create the random perturbation and store it on temporary array
    IF ( random_generator == 'numerical-recipes' )  THEN
       DO  i = dist_nxl(dist_range), dist_nxr(dist_range)
          DO  j = dist_nys(dist_range), dist_nyn(dist_range)
             DO  k = disturbance_level_ind_b, disturbance_level_ind_t
                randomnumber = 3.0_wp * disturbance_amplitude *                &
                               ( random_function( iran ) - 0.5_wp )
                IF ( nxl <= i  .AND.  nxr >= i  .AND.  nys <= j  .AND.         &
                     nyn >= j )                                                &
                THEN
                   dist1(k,j,i) = randomnumber
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ELSEIF ( random_generator == 'random-parallel' )  THEN
       DO  i = dist_nxl(dist_range), dist_nxr(dist_range)
          DO  j = dist_nys(dist_range), dist_nyn(dist_range)
             CALL random_seed_parallel( put=seq_random_array(:, j, i) )
             DO  k = disturbance_level_ind_b, disturbance_level_ind_t
                CALL random_number_parallel( random_dummy )
                randomnumber = 3.0_wp * disturbance_amplitude *                &
                               ( random_dummy - 0.5_wp )
                IF ( nxl <= i  .AND.  nxr >= i  .AND.  nys <= j  .AND.         &
                     nyn >= j )                                                &
                THEN
                   dist1(k,j,i) = randomnumber
                ENDIF
             ENDDO
             CALL random_seed_parallel( get=seq_random_array(:, j, i) )
          ENDDO
       ENDDO
    ELSEIF ( random_generator == 'system-specific' )  THEN
       DO  i = dist_nxl(dist_range), dist_nxr(dist_range)
          DO  j = dist_nys(dist_range), dist_nyn(dist_range)
             DO  k = disturbance_level_ind_b, disturbance_level_ind_t
                CALL RANDOM_NUMBER( randomnumber )
                randomnumber = 3.0_wp * disturbance_amplitude *                &
                                ( randomnumber - 0.5_wp )
                IF ( nxl <= i .AND. nxr >= i .AND. nys <= j .AND. nyn >= j )   &
                THEN
                   dist1(k,j,i) = randomnumber
                ENDIF
             ENDDO
          ENDDO
       ENDDO

    ENDIF

!
!-- Update dist1 on the device, this is expected by exchange_horiz!
    !$ACC UPDATE DEVICE(dist1(nzb:nzt+1,nysg:nyng,nxlg:nxrg))

!
!-- Exchange of ghost points for the random perturbation

    CALL exchange_horiz( dist1, nbgp )
!
!-- Applying the Shuman filter in order to smooth the perturbations.
!-- Neighboured grid points in all three directions are used for the
!-- filter operation.
!-- Loop has been splitted to make runs reproducible on HLRN systems using
!-- compiler option -O3
     !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(i, j, k) PRESENT(dist1, dist2)
     !$OMP PARALLEL DO PRIVATE(i, j, k)
     DO  i = nxl, nxr
        DO  j = nys, nyn
          DO  k = disturbance_level_ind_b-1, disturbance_level_ind_t+1
             dist2(k,j,i) = ( dist1(k,j,i-1) + dist1(k,j,i+1)                  &
                            + dist1(k,j+1,i) + dist1(k+1,j,i)                  &
                            ) / 12.0_wp
          ENDDO
          DO  k = disturbance_level_ind_b-1, disturbance_level_ind_t+1
              dist2(k,j,i) = dist2(k,j,i) + ( dist1(k,j-1,i) + dist1(k-1,j,i)  &
                            + 6.0_wp * dist1(k,j,i)                            &
                            ) / 12.0_wp
          ENDDO
        ENDDO
     ENDDO

!
!-- Exchange of ghost points for the filtered perturbation.
!-- Afterwards, filter operation and exchange of ghost points are repeated.
    CALL exchange_horiz( dist2, nbgp )

    !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(i, j, k) PRESENT(dist1, dist2)
    !$OMP PARALLEL DO PRIVATE(i, j, k)
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = disturbance_level_ind_b-2, disturbance_level_ind_t+2
             dist1(k,j,i) = ( dist2(k,j,i-1) + dist2(k,j,i+1) + dist2(k,j-1,i) &
                            + dist2(k,j+1,i) + dist2(k+1,j,i) + dist2(k-1,j,i) &
                            + 6.0_wp * dist2(k,j,i)                            &
                            ) / 12.0_wp
          ENDDO
       ENDDO
    ENDDO

    CALL exchange_horiz( dist1, nbgp )

!
!-- Remove perturbations below topography (including one gridpoint above it
!-- in order to allow for larger timesteps at the beginning of the simulation
!-- (diffusion criterion))
    IF ( TRIM( topography ) /= 'flat' )  THEN
       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
             DO  k = nzb, nzb_max
                dist1(k,j,i) = MERGE( dist1(k,j,i), 0.0_wp,                    &
                                  BTEST( wall_flags_total_0(k,j,i), flag_nr )  &
                                    )
             ENDDO
          ENDDO
       ENDDO
    ENDIF

!
!-- Random perturbation is added to the array to be disturbed.
    !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i, j, k) PRESENT(field, dist1)
    !$OMP PARALLEL DO PRIVATE(i, j, k)
    DO  i = nxlg, nxrg
       DO  j = nysg, nyng
          DO  k = disturbance_level_ind_b-2, disturbance_level_ind_t+2
             field(k,j,i) = field(k,j,i) + dist1(k,j,i)
          ENDDO
       ENDDO
    ENDDO

!$ACC END DATA

!
!-- Deallocate the temporary array
    DEALLOCATE( dist2 )

!
!-- Set a flag, which indicates that a random perturbation is imposed
    disturbance_created = .TRUE.


    CALL cpu_log( log_point(20), 'disturb_field', 'stop' )


 END SUBROUTINE disturb_field
