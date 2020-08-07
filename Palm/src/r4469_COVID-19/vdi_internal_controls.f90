!> @file vdi_internal_controls.f90
!--------------------------------------------------------------------------------!
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
! Copyright 2019-2019 Leibniz Universitaet Hannover
!--------------------------------------------------------------------------------!
!
! Current revisions:
! -----------------
! 
! 
! Former revisions:
! -----------------
! $Id: vdi_internal_controls.f90 4415 2020-02-20 10:30:33Z raasch $
! missing preprocessor directive added
! 
! 4346 2019-12-18 11:55:56Z motisi
! Introduction of wall_flags_total_0, which currently sets bits based on static
! topography information used in wall_flags_static_0
! 
! 4329 2019-12-10 15:46:36Z motisi
! Renamed wall_flags_0 to wall_flags_static_0
! 
! 4182 2019-08-22 15:20:23Z scharf
! added "Authors" section
! 
! 4175 2019-08-20 13:19:16Z gronemeier
! bugfix: removed unused variables
!
! 4173 2019-08-20 12:04:06Z weniger
! Initial version
!
! Authors:
! --------
! @author Viola Weniger
!
!
! Description:
! ------------
!> According to VDI Guideline 3783 Part 9, internal assessment have to be
!> carried out within the program for the model to be considered as evaluated.
!------------------------------------------------------------------------------!
 MODULE vdi_internal_controls

    USE arrays_3d,                          &
        ONLY:  dzw,                         &
               pt,                          &
               q,                           &
               u,                           &
               u_p,                         &
               v,                           &
               w

    USE control_parameters,                 &
        ONLY:  bc_dirichlet_l,              &
               bc_dirichlet_n,              &
               bc_dirichlet_r,              &
               bc_dirichlet_s,              &
               bc_lr_cyc,                   &
               bc_ns_cyc,                   &
               humidity,                    &
               end_time,                    &
               message_string,              &
               neutral,                     &
               time_since_reference_point

    USE indices,                            &
        ONLY:  nx,                          &
               nxl,                         &
               nxlg,                        &
               nxr,                         &
               nxrg,                        &
               ny,                          &
               nyn,                         &
               nyng,                        &
               nys,                         &
               nysg,                        &
               nzb,                         &
               nzt,                         &
               wall_flags_total_0

    USE kinds

#if defined( __parallel )
    USE pegrid,                             &
        ONLY:  collective_wait,             &
               comm2d,                      &
               ierr,                        &
               MPI_DOUBLE_PRECISION,        &
               MPI_INTEGER,                 &
               MPI_MAX,                     &
               MPI_SUM,                     &
               myid
#else
    USE pegrid,                             &
        ONLY:  myid
#endif


    USE grid_variables,                     &
        ONLY:  dx,                          &
               dy

    USE pmc_interface,                      &
        ONLY: nested_run

    IMPLICIT NONE
    
    INTEGER(iwp) ::  internal_count = 0  !< counts calls to this module

    INTERFACE vdi_2_deltat_wave
       MODULE PROCEDURE vdi_2_deltat_wave
    END INTERFACE vdi_2_deltat_wave

    INTERFACE vdi_standard_differences
       MODULE PROCEDURE vdi_standard_differences
    END INTERFACE vdi_standard_differences

    INTERFACE vdi_domain_averages
       MODULE PROCEDURE vdi_domain_averages
    END INTERFACE vdi_domain_averages

    INTERFACE vdi_plausible_values
       MODULE PROCEDURE vdi_plausible_values
    END INTERFACE vdi_plausible_values

    INTERFACE vdi_actions
       MODULE PROCEDURE vdi_actions
    END INTERFACE vdi_actions

    INTERFACE vdi_conservation_of_mass
       MODULE PROCEDURE vdi_conservation_of_mass
    END INTERFACE vdi_conservation_of_mass

    SAVE

    PRIVATE

!
!-- Public functions
    PUBLIC          &
       vdi_actions


 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!> @todo Add proper description
!------------------------------------------------------------------------------!
 SUBROUTINE vdi_actions( location )

    CHARACTER (LEN=*), INTENT(IN) ::  location  !< call location string


    SELECT CASE ( location )

       CASE ( 'after_integration' )

          internal_count = internal_count + 1
       
          CALL vdi_2_deltat_wave

          CALL vdi_standard_differences

          CALL vdi_domain_averages

          CALL vdi_conservation_of_mass

          CALL vdi_plausible_values

       CASE DEFAULT

          CONTINUE

    END SELECT

 END SUBROUTINE vdi_actions
!------------------------------------------------------------------------------!
! Description:
! ------------
!> At a control grid point in the interior of the model domain,
!> 2 * deltat waves are not to be generated with increasing simulation time.
!------------------------------------------------------------------------------!
 SUBROUTINE vdi_2_deltat_wave

    INTEGER(iwp) ::  count_wave = 0  !< counts the number of consecutive waves
    INTEGER(iwp) ::  count_time = 0  !< counter, so that the waves follow each other without gaps
    INTEGER(iwp) ::  cgp_i = 0       !< x coordinate of the control grid point for testing 2deltat waves
    INTEGER(iwp) ::  cgp_j = 0       !< y coordinate of the control grid point for testing 2deltat waves
    INTEGER(iwp) ::  cgp_k = 0       !< z coordinate of the control grid point for testing 2deltat waves

    INTEGER(iwp), DIMENSION(4) ::  sig_arr = (/ 0, 0, 0, 0/)  !< indicates an increase(1) or a decrease (0) of u in the last four time steps

    REAL(wp) ::  random  !< random number

!
!-- Defining the control grid point
    IF ( internal_count == 1 )  THEN
       cgp_i = INT( nxl + ( nxr - nxl ) / 2 )
       cgp_j = INT( nys + ( nyn - nys ) / 2 )
       cgp_k = INT( nzt / 2 )
!
!--    If the grid point lies in a building, a new point is defined
       DO WHILE ( .NOT. BTEST( wall_flags_total_0(cgp_k,cgp_j,cgp_i), 1 ) )
          CALL RANDOM_NUMBER( random )
          cgp_k = cgp_k + FLOOR( ( nzt - cgp_k ) * random )   !< Random number upon cgp_k
!
!--       If there is topography in the entire grid column, a new x coordinate is chosen
          IF ( cgp_k >= nzt-1 )  THEN
             CALL RANDOM_NUMBER( random )
             cgp_i = nxl + FLOOR( ( nxr + 1 - nxl ) * random )
             cgp_k = INT( nzt / 2 )
          ENDIF
       ENDDO
    ENDIF

    CALL testing_2_deltat_wave( u_p(cgp_k,cgp_j,cgp_i), u(cgp_k,cgp_j,cgp_i), &
                                sig_arr, count_wave, count_time )

 END SUBROUTINE vdi_2_deltat_wave


!------------------------------------------------------------------------------!
! Description:
! ------------
!> In this subroutine a quantity quant is tested for 2 delta t waves.
!> For this, the size must have a wave-shaped course over 4*4 time steps
!> and the amplitude of the wave has to be greater than the change of quant with
!> increasing time.
!------------------------------------------------------------------------------!
SUBROUTINE testing_2_deltat_wave( quant_p_r, quant_r, sig_arr, count_wave, count_time )

    INTEGER(iwp), INTENT(INOUT) ::  count_wave        !< counts the number of consecutive waves
    INTEGER(iwp), INTENT(INOUT) ::  count_time        !< counter, so that the waves follow each other without gaps
    INTEGER(iwp), PARAMETER     ::  number_wave = 10  !< number of consecutive waves that are not allowed

    REAL(wp), INTENT(IN) ::  quant_p_r                !< quantity from the previous time step as a real
    REAL(wp), INTENT(IN) ::  quant_r                  !< quantity as a real
    REAL(wp)             ::  quant_rel = 0.0_wp       !< rel. change of the quantity to the previous time step

    INTEGER(iwp), DIMENSION(4), INTENT(INOUT) ::  sig_arr !< indicates an increase(1) or a decrease (0) of
                                                          !> quantity quant in the last four time steps


    IF ( quant_p_r - quant_r > 0.0 )  THEN
       sig_arr(4) = 0
    ELSE
       sig_arr(4) = 1
    ENDIF

    quant_rel = ABS( ( quant_p_r - quant_r ) / quant_p_r )

!
!-- With this criterion 2 delta t waves are detected if the amplitude of
!-- the wave is greater than the change of quant with increasing time
    IF ( ALL( sig_arr(1:4) == (/ 1, 0, 1, 0 /) )  .AND.  quant_rel > 0.01 )  THEN

       count_wave = count_wave + 1

       IF ( count_wave == number_wave  .AND.  count_time == 4 )  THEN
          message_string = '2 deltat waves are generated '
          CALL message( 'vdi_2_deltat_wave', 'PA0669', 2, 2, myid, 6, 0 )
       ENDIF

       count_time = 0

    ELSE

       IF ( count_time >= 4 )  THEN
          count_wave = 0
       ENDIF

    ENDIF

    sig_arr(1) = sig_arr(2)
    sig_arr(2) = sig_arr(3)
    sig_arr(3) = sig_arr(4)

    count_time = count_time + 1

 END SUBROUTINE testing_2_deltat_wave


!------------------------------------------------------------------------------!
! Description:
! ------------
!> In this internal assessment the maxima of standarddifferences of the
!> meteorological variables, computed layer by layer will be checked.
!> The maxima should not to remain at the open edges of the model or
!> travel from there into the interior of the domain with increasing
!> simulation time.
!> @todo try to reduce repeating code.
!------------------------------------------------------------------------------!
SUBROUTINE vdi_standard_differences

    INTEGER(iwp) ::  position_u_deviation = 0     !< position of the maximum of the standard deviation of u
    INTEGER(iwp) ::  position_u_deviation_p = 0   !< position of the maximum of the standard deviation of u to the previous time step
    INTEGER(iwp) ::  position_u_deviation_pp = 0  !< position of the maximum of the standard deviation of u two time steps ago
    INTEGER(iwp) ::  position_v_deviation = 0     !< position of the maximum of the standard deviation of v
    INTEGER(iwp) ::  position_v_deviation_p = 0   !< position of the maximum of the standard deviation of v to the previous time step
    INTEGER(iwp) ::  position_v_deviation_pp = 0  !< position of the maximum of the standard deviation of v two time steps ago
    INTEGER(iwp) ::  position_w_deviation = 0     !< position of the maximum of the standard deviation of w
    INTEGER(iwp) ::  position_w_deviation_p = 0   !< position of the maximum of the standard deviation of w to the previous time step
    INTEGER(iwp) ::  position_w_deviation_pp = 0  !< position of the maximum of the standard deviation of w two time steps ago
    INTEGER(iwp) ::  position_pt_deviation = 0    !< position of the maximum of the standard deviation of pt
    INTEGER(iwp) ::  position_pt_deviation_p = 0  !< position of the maximum of the standard deviation of pt to the previous time step
    INTEGER(iwp) ::  position_pt_deviation_pp = 0 !< position of the maximum of the standard deviation of pt two time steps ago
    INTEGER(iwp) ::  position_q_deviation = 0     !< position of the maximum of the standard deviation of q
    INTEGER(iwp) ::  position_q_deviation_p = 0   !< position of the maximum of the standard deviation of q to the previous time step
    INTEGER(iwp) ::  position_q_deviation_pp = 0  !< position of the maximum of the standard deviation of q two time steps ago

    REAL(wp), DIMENSION(nzb:nzt+1) ::  u_deviation  !< standard deviation of u depending on k
    REAL(wp), DIMENSION(nzb:nzt+1) ::  v_deviation  !< standard deviation of v depending on k
    REAL(wp), DIMENSION(nzb:nzt+1) ::  w_deviation  !< standard deviation of w depending on k
    REAL(wp), DIMENSION(nzb:nzt+1) ::  pt_deviation !< standard deviation of pt depending on k
    REAL(wp), DIMENSION(nzb:nzt+1) ::  q_deviation  !< standard deviation of q depending on k

!
!-- Calculation of the standard deviation of u
    CALL calc_standard_deviation( u, u_deviation, 1 )

!
!-- Determination of the position of the maximum
    position_u_deviation = MAXLOC( u_deviation, DIM=1 )

!
!-- Check the position of the maximum of the standard deviation of u
    IF ( internal_count > 2 )  THEN
       CALL check_position( position_u_deviation, position_u_deviation_p, position_u_deviation_pp )
    ENDIF

    position_u_deviation_pp = position_u_deviation_p
    position_u_deviation_p = position_u_deviation

!
!-- Calculation of the standard deviation of v
    CALL calc_standard_deviation( v, v_deviation, 2 )

!
!-- Determination of the position of the maximum
    position_v_deviation = MAXLOC( v_deviation, DIM=1 )

!
!-- Check the position of the maximum of the standard deviation of v
    IF ( internal_count > 2 )  THEN
       CALL check_position( position_v_deviation, position_v_deviation_p, position_v_deviation_pp )
    ENDIF

   position_v_deviation_pp = position_v_deviation_p
   position_v_deviation_p = position_v_deviation

!
!-- Calculation of the standard deviation of w
    CALL calc_standard_deviation( w, w_deviation, 3 )

!
!-- Determination of the position of the maximum
    position_w_deviation = MAXLOC( w_deviation, DIM=1 )

!
!-- Check the position of the maximum of the standard deviation of w
    IF ( internal_count > 2 )  THEN
       CALL check_position( position_w_deviation, position_w_deviation_p, position_w_deviation_pp )
    ENDIF

    position_w_deviation_pp = position_w_deviation_p
    position_w_deviation_p = position_w_deviation


!
!-- Calculation of the standard deviation of pt
    IF ( .NOT. neutral )  THEN
       CALL calc_standard_deviation( pt, pt_deviation, 0 )
!
!--    Determination of the position of the maximum
       position_pt_deviation = MAXLOC( pt_deviation, DIM=1 )

!
!--    Check the position of the maximum of the standard deviation of pt
       IF ( internal_count > 2 )  THEN
          CALL check_position( position_pt_deviation,   &
                               position_pt_deviation_p, &
                               position_pt_deviation_pp )
       ENDIF

       position_pt_deviation_pp = position_pt_deviation_p
       position_pt_deviation_p = position_pt_deviation

    ENDIF

!
!-- Calculation of the standard deviation of q
    IF ( humidity )  THEN
       CALL calc_standard_deviation( q, q_deviation, 0 )

!
!--    Determination of the position of the maximum
       position_q_deviation = MAXLOC( q_deviation, DIM=1 )

!
!--    Check the position of the maximum of the standard deviation of q
       IF ( internal_count > 2 )  THEN
          CALL check_position( position_q_deviation,   &
                               position_q_deviation_p, &
                               position_q_deviation_pp )
       ENDIF

       position_q_deviation_pp = position_q_deviation_p
       position_q_deviation_p = position_q_deviation

    ENDIF

END SUBROUTINE vdi_standard_differences


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculation of the standard deviation
!------------------------------------------------------------------------------!
SUBROUTINE calc_standard_deviation( quant, std_deviation, quant_type )

    INTEGER(iwp)             ::  i           !< loop index
    INTEGER(iwp)             ::  j           !< loop index
    INTEGER(iwp)             ::  k           !< loop index
    INTEGER(iwp), INTENT(IN) ::  quant_type  !< bit position (1 for u, 2 for v, 3 for w and 0 for scalar)

    INTEGER(iwp), DIMENSION(nzb:nzt+1) ::  count_2d_l  !< counter for averaging (local)
    INTEGER(iwp), DIMENSION(nzb:nzt+1) ::  count_2d    !< counter for averaging

    REAL(wp) ::  flag  !< flag indicating atmosphere (1) or wall (0) grid point

    REAL(wp), DIMENSION(nzb:nzt+1)              ::  quant_av_k_l     !< Mean of the quantity quant depending on k (local)
    REAL(wp), DIMENSION(nzb:nzt+1)              ::  quant_av_k       !< Mean of the quantity quant depending on k
    REAL(wp), DIMENSION(nzb:nzt+1), INTENT(OUT) ::  std_deviation    !< standard deviation of quant
    REAL(wp), DIMENSION(nzb:nzt+1)              ::  std_deviation_l  !< standard deviation of quant (local)

    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg), INTENT(IN) ::  quant !< Quantity

!
!-- Calculation of the standard deviation
    quant_av_k_l    = 0.0_wp
    quant_av_k      = 0.0_wp
    std_deviation   = 0.0_wp
    std_deviation_l = 0.0_wp
!
!-- Average
    count_2d_l = 0
    count_2d   = 0
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt+1
             flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), quant_type ) )
             quant_av_k_l(k) = quant_av_k_l(k) + quant(k,j,i) * flag
             count_2d_l(k)   = count_2d_l(k) + INT( flag, KIND=iwp )
          ENDDO
       ENDDO
    ENDDO

#if defined( __parallel )
    CALL MPI_ALLREDUCE( quant_av_k_l, quant_av_k, nzt+1-nzb+1, &
                        MPI_REAL, MPI_SUM, comm2d, ierr )

    CALL MPI_ALLREDUCE( count_2d_l, count_2d, nzt+1-nzb+1,     &
                        MPI_INTEGER, MPI_SUM, comm2d, ierr )
#else
    quant_av_k = quant_av_k_l
    count_2d   = count_2d_l
#endif

    DO  k = nzb+1, nzt+1
       quant_av_k(k) = quant_av_k(k) / REAL( count_2d(k), KIND=wp )
    ENDDO

    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt+1
             std_deviation_l(k) = std_deviation_l(k)                  &
                                + ( quant(k,j,i) - quant_av_k(k) )**2 &
                                * MERGE( 1.0_wp, 0.0_wp,              &
                                         BTEST( wall_flags_total_0(k,j,i), quant_type ) )
          ENDDO
       ENDDO
    ENDDO


#if defined( __parallel )
    CALL MPI_ALLREDUCE( std_deviation_l, std_deviation, nzt+1-nzb+1, &
                        MPI_REAL, MPI_SUM, comm2d, ierr )
#else
    std_deviation = std_deviation_l
#endif

    DO  k = nzb+1, nzt+1
       std_deviation(k) = SQRT( std_deviation(k) / REAL( count_2d(k), KIND=wp ) )
    ENDDO

END SUBROUTINE calc_standard_deviation


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Tests for the position of the maxima of the standard deviation.
!> If the maxima remain at the open edges of the model or travel from
!> the open edges into the interior of the domain with increasing
!> simulation time, the simulation should be aborted.
!------------------------------------------------------------------------------!
 SUBROUTINE check_position( position_std_deviation, position_std_deviation_p, &
                            position_std_deviation_pp )

    INTEGER(iwp), INTENT(IN) ::  position_std_deviation     !< position of the maximum of the std
    INTEGER(iwp), INTENT(IN) ::  position_std_deviation_p   !< previous position of std-max
    INTEGER(iwp), INTENT(IN) ::  position_std_deviation_pp  !< prev. prev. position of std-max


    IF ( position_std_deviation == nzt    .AND.  &
         position_std_deviation_p == nzt  .AND.  &
         position_std_deviation_pp == nzt        )  THEN
       message_string = 'The maxima of the standard deviation remain ' // &
                        'at the open edges of the model.'
       CALL message( 'vdi_standard_differences', 'PA0663', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( position_std_deviation == nzt-2    .AND. &
         position_std_deviation_p == nzt-1  .AND. &
         position_std_deviation_pp == nzt         )  THEN
       message_string = 'The maxima of the standard deviation travel ' // &
                        'from the open edges into the interior ' //       &
                        'of the domain with increasing simulation time.'
       CALL message( 'vdi_standard_differences', 'PA0664', 1, 2, 0, 6, 0 )
    ENDIF

END SUBROUTINE check_position


!------------------------------------------------------------------------------!
! Description:
! ------------
!> In this control it will be checked, if the means of the meteorological
!> variables over the model grid are not to exhibit 2 deltat waves or
!> monotonic increase or decrease with increasing simulation time.
!------------------------------------------------------------------------------!
SUBROUTINE vdi_domain_averages

   INTEGER(iwp) ::  mono_count_u = 0   !< counter for monotonic decrease or increase of u
   INTEGER(iwp) ::  mono_count_v = 0   !< counter for monotonic decrease or increase of v
   INTEGER(iwp) ::  mono_count_w = 0   !< counter for monotonic decrease or increase of w
   INTEGER(iwp) ::  mono_count_q = 0   !< counter for monotonic decrease or increase of q
   INTEGER(iwp) ::  mono_count_pt = 0  !< counter for monotonic decrease or increase of pt
   INTEGER(iwp) ::  count_time_u = 0   !< counter, so that the waves of u follow each other without gaps
   INTEGER(iwp) ::  count_time_v = 0   !< counter, so that the waves of v follow each other without gaps
   INTEGER(iwp) ::  count_time_w = 0   !< counter, so that the waves of w follow each other without gaps
   INTEGER(iwp) ::  count_time_q = 0   !< counter, so that the waves of q follow each other without gaps
   INTEGER(iwp) ::  count_time_pt = 0  !< counter, so that the waves of pt follow each other without gaps
   INTEGER(iwp) ::  count_wave_u = 0   !< counts the number of consecutive waves of u
   INTEGER(iwp) ::  count_wave_v = 0   !< counts the number of consecutive waves of v
   INTEGER(iwp) ::  count_wave_w = 0   !< counts the number of consecutive waves of w
   INTEGER(iwp) ::  count_wave_q = 0   !< counts the number of consecutive waves of q
   INTEGER(iwp) ::  count_wave_pt = 0  !< counts the number of consecutive waves of pt

   INTEGER(iwp), DIMENSION(4) ::  sig_u_arr = (/ 0, 0, 0, 0/)   !< indicates an increase(1) or a decrease (0) of u in the last four time steps
   INTEGER(iwp), DIMENSION(4) ::  sig_v_arr = (/ 0, 0, 0, 0/)   !< indicates an increase(1) or a decrease (0) of v in the last four time steps
   INTEGER(iwp), DIMENSION(4) ::  sig_w_arr = (/ 0, 0, 0, 0/)   !< indicates an increase(1) or a decrease (0) of w in the last four time steps
   INTEGER(iwp), DIMENSION(4) ::  sig_q_arr = (/ 0, 0, 0, 0/)   !< indicates an increase(1) or a decrease (0) of q in the last four time steps
   INTEGER(iwp), DIMENSION(4) ::  sig_pt_arr = (/ 0, 0, 0, 0/)  !< indicates an increase(1) or a decrease (0) of pt in the last four time steps

   REAL(wp) ::  u_av = 0.0_wp     !< Mean of u
   REAL(wp) ::  u_av_p = 0.0_wp   !< Mean of u at the previous time step
   REAL(wp) ::  v_av = 0.0_wp     !< Mean of v
   REAL(wp) ::  v_av_p = 0.0_wp   !< Mean of v at the previous time step
   REAL(wp) ::  w_av = 0.0_wp     !< Mean of w
   REAL(wp) ::  w_av_p = 0.0_wp   !< Mean of w at the previous time step
   REAL(wp) ::  q_av = 0.0_wp     !< Mean of q
   REAL(wp) ::  q_av_p = 0.0_wp   !< Mean of q at the previous time step
   REAL(wp) ::  pt_av = 0.0_wp    !< Mean of pt
   REAL(wp) ::  pt_av_p = 0.0_wp  !< Mean of pt at the previous time step

!
!-- Averaging the meteorological variables over the model grid
    CALL calc_average( u, u_av, 1 )
    CALL calc_average( v, v_av, 2 )
    CALL calc_average( w, w_av, 3 )
    IF ( .NOT. neutral )  THEN
       CALL calc_average( pt, pt_av, 0 )
    ENDIF
    IF ( humidity )  THEN
       CALL calc_average( q, q_av, 0 )
    ENDIF

!
!-- Testing the meteorological variables for 2 delta t waves
    IF ( internal_count > 1 )  THEN
       CALL testing_2_deltat_wave( u_av_p, u_av, sig_u_arr, count_wave_u, count_time_u )
       CALL testing_2_deltat_wave( v_av_p, v_av, sig_v_arr, count_wave_v, count_time_v )
       CALL testing_2_deltat_wave( w_av_p, w_av, sig_w_arr, count_wave_w, count_time_w )
       IF ( .NOT. neutral )  THEN
          CALL testing_2_deltat_wave( pt_av_p, pt_av, sig_pt_arr, count_wave_pt, count_time_pt )
       ENDIF
       IF ( humidity )  THEN
          CALL testing_2_deltat_wave( q_av_p, q_av, sig_q_arr, count_wave_q, count_time_q )
       ENDIF
    ENDIF

!
!-- Testing if there is a monotonic increase or decrease with increasing simulation time
    IF ( sig_u_arr(2) /= sig_u_arr(3) )  THEN
       mono_count_u = 0
    ELSE
       mono_count_u = mono_count_u + 1
    ENDIF

    IF ( time_since_reference_point >= end_time  .AND.   &
         mono_count_u > 0.9_wp * internal_count )  THEN

       message_string = 'Monotonic decrease or increase with ' // &
                        'increasing simulation time for u'
       CALL message( 'vdi_domain_averages', 'PA0665', 0, 1, 0, 6, 0 )
    ENDIF

    IF ( sig_v_arr(2) /= sig_v_arr(3) )  THEN
       mono_count_v = 0
    ELSE
       mono_count_v = mono_count_v + 1
    ENDIF

    IF ( time_since_reference_point >= end_time  .AND.   &
         mono_count_v > 0.9_wp * internal_count )  THEN
       message_string = 'Monotonic decrease or increase with ' // &
                        'increasing simulation time for v'
       CALL message( 'vdi_domain_averages', 'PA0665', 0, 1, 0, 6, 0 )
    ENDIF

    IF ( sig_w_arr(2) /= sig_w_arr(3) )  THEN
       mono_count_w = 0
    ELSE
       mono_count_w = mono_count_w + 1
    ENDIF

    IF ( time_since_reference_point >= end_time  .AND.   &
         mono_count_w > 0.9_wp * internal_count )  THEN
       message_string = 'Monotonic decrease or increase with ' // &
                        'increasing simulation time for w'
       CALL message( 'vdi_domain_averages', 'PA0665', 0, 1, 0, 6, 0 )
    ENDIF

    IF ( .NOT. neutral )  THEN
       IF ( sig_pt_arr(2) /= sig_pt_arr(3) )  THEN
          mono_count_pt = 0
       ELSE
          mono_count_pt = mono_count_pt + 1
       ENDIF

       IF ( time_since_reference_point >= end_time  .AND.    &
            mono_count_pt > 0.9_wp * internal_count )  THEN
          message_string = 'Monotonic decrease or increase with ' // &
                           'increasing simulation time for pt'
          CALL message( 'vdi_domain_averages', 'PA0665', 0, 1, 0, 6, 0 )
       ENDIF
    ENDIF

    IF ( humidity )  THEN
       IF ( sig_q_arr(2) /= sig_q_arr(3) )  THEN
          mono_count_q = 0
       ELSE
          mono_count_q = mono_count_q + 1
       ENDIF

       IF ( time_since_reference_point >= end_time  .AND.   &
            mono_count_q > 0.9_wp * internal_count )  THEN
          message_string = 'Monotonic decrease or increase with ' // &
                           'increasing simulation time for q'
          CALL message( 'vdi_domain_averages', 'PA0665', 0, 1, 0, 6, 0 )
       ENDIF
    ENDIF

!
!-- Save the values from the previous time step
    u_av_p = u_av
    v_av_p = v_av
    w_av_p = w_av

    IF ( .NOT. neutral )  THEN
       pt_av_p = pt_av
    ENDIF

    IF ( humidity )  THEN
       q_av_p = q_av
    ENDIF

 END SUBROUTINE vdi_domain_averages


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate the average of a quantity 'quant'.
!------------------------------------------------------------------------------!
 SUBROUTINE calc_average( quant, quant_av, quant_type )

    INTEGER(iwp) ::  average_count = 0    !< counter for averaging
    INTEGER(iwp) ::  average_count_l = 0  !< counter for averaging (local)
    INTEGER      ::  i                    !< loop index
    INTEGER      ::  j                    !< loop index
    INTEGER      ::  k                    !< loop index
    INTEGER(iwp) ::  quant_type           !< bit position (1 for u, 2 for v, 3 for w and 0 for scalar)

    REAL(wp) ::  flag                     !< flag indicating atmosphere (1) or wall (0) grid point
    REAL(wp) ::  quant_av                 !< average of the quantity quant
    REAL(wp) ::  quant_av_l = 0.0_wp      !< average of the quantity quant (local)

    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  quant

!
!-- Averaging the quantity over the model grid
   average_count_l = 0
   quant_av_l = 0.0_wp
   DO  i = nxl, nxr
      DO  j = nys, nyn
         DO  k = nzb, nzt+1
            flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i), quant_type ) )
            quant_av_l = quant_av_l + quant(k,j,i) * flag
            average_count_l = average_count_l + INT( flag, KIND=iwp )
         ENDDO
      ENDDO
   ENDDO

#if defined( __parallel )
    CALL MPI_ALLREDUCE( quant_av_l, quant_av, 1,        &
                        MPI_REAL, MPI_SUM, comm2d, ierr )
    CALL MPI_ALLREDUCE( average_count_l, average_count, 1, &
                        MPI_INTEGER, MPI_SUM, comm2d, ierr )
#else
    quant_av = quant_av_l
    average_count = average_count_l
#endif

    quant_av = quant_av / REAL( average_count, KIND(wp) )

 END SUBROUTINE calc_average


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Testing for conservation of mass.
!------------------------------------------------------------------------------!
 SUBROUTINE vdi_conservation_of_mass

    INTEGER(iwp) ::  i              !< loop index
    INTEGER(iwp) ::  j              !< loop index
    INTEGER(iwp) ::  k              !< loop index

    REAL(wp)     ::  sum_mass_flux  !< sum of the mass flow

    REAL(wp), DIMENSION(1:3) ::  volume_flow_l  !< volume flow (local)
    REAL(wp), DIMENSION(1:3) ::  volume_flow    !< volume flow


    volume_flow   = 0.0_wp
    volume_flow_l = 0.0_wp

!
!-- Left/right:
!-- Sum up the volume flow through the left boundary
    IF ( nxl == 0 )  THEN
       i = 0
       DO  j = nys, nyn
          DO  k = nzb+1, nzt
             volume_flow_l(1) = volume_flow_l(1)                        &
                              + u(k,j,i) * dzw(k) * dy                  &
                              * MERGE( 1.0_wp, 0.0_wp,                  &
                                 BTEST( wall_flags_total_0(k,j,i), 1 )  &
                                     )
          ENDDO
       ENDDO
    ENDIF
! 
!-- Sum up the volume flow through the right boundary
    IF ( nxr == nx )  THEN
       i = nx+1
       DO  j = nys, nyn
          DO  k = nzb+1, nzt
             volume_flow_l(1) = volume_flow_l(1)                        &
                              - u(k,j,i) * dzw(k) * dy                  &
                              * MERGE( 1.0_wp, 0.0_wp,                  &
                                 BTEST( wall_flags_total_0(k,j,i), 1 )  &
                                     )
          ENDDO
       ENDDO
    ENDIF
!
!-- South/north:
!-- Sum up the volume flow through the south boundary
    IF ( nys == 0 )  THEN
       j = 0
       DO  i = nxl, nxr
          DO  k = nzb+1, nzt
             volume_flow_l(2) = volume_flow_l(2)                        &
                              + v(k,j,i) * dzw(k) * dx                  &
                              * MERGE( 1.0_wp, 0.0_wp,                  &
                                 BTEST( wall_flags_total_0(k,j,i), 2 )  &
                                     )
          ENDDO
       ENDDO
    ENDIF
!
!-- Sum up the volume flow through the north boundary
    IF ( nyn == ny )  THEN
       j = ny+1
       DO  i = nxl, nxr
          DO  k = nzb+1, nzt
             volume_flow_l(2) = volume_flow_l(2)                        &
                              - v(k,j,i) * dzw(k) * dx                  &
                              * MERGE( 1.0_wp, 0.0_wp,                  &
                                 BTEST( wall_flags_total_0(k,j,i), 2 )  &
                                     )
          ENDDO
       ENDDO
    ENDIF
!
!-- Top boundary
    k = nzt
    DO  i = nxl, nxr
       DO  j = nys, nyn
          volume_flow_l(3) = volume_flow_l(3) - w(k,j,i) * dx * dy
       ENDDO
    ENDDO

#if defined( __parallel )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( volume_flow_l, volume_flow, 3, MPI_REAL, MPI_SUM, comm2d, ierr )
#else
    volume_flow = volume_flow_l
#endif

    sum_mass_flux = SUM( volume_flow ) / ( ( nx + 1 ) * dx * ( ny + 1 ) * dy )

   IF ( ABS( sum_mass_flux ) > 0.001 )  THEN
      message_string = 'The mass is not conserved. '
      CALL message( 'vdi_conservation_of_mass', 'PA0666', 1, 2, 0, 6, 0 )
   ENDIF

END SUBROUTINE vdi_conservation_of_mass


!------------------------------------------------------------------------------!
! Description:
! ------------
!> The results will be checked for exceedance the specified limits.
!> The controls are performed at every time step and at every grid point.
!> No wind component is allowed to have a magnitude greater than ten times
!> the maximum wind velocity at the approach flow profile (Vdi 3783 part 9).
!> Note, that the supersaturation can not be higher than 10%. Therefore, no
!> test is required.
!------------------------------------------------------------------------------!
SUBROUTINE vdi_plausible_values

    INTEGER(iwp) ::  i          !< loop index
    INTEGER(iwp) ::  j          !< loop index
    INTEGER(iwp) ::  k          !< loop index

    REAL(wp)     :: max_uv_l_l  !< maximum speed at the left edge (local)
    REAL(wp)     :: max_uv_l    !< maximum speed at the left edge
    REAL(wp)     :: max_uv_r_l  !< maximum speed at the right edge (local)
    REAL(wp)     :: max_uv_r    !< maximum speed at the right edge
    REAL(wp)     :: max_uv_s_l  !< maximum speed at the south edge (local)
    REAL(wp)     :: max_uv_s    !< maximum speed at the south edge
    REAL(wp)     :: max_uv_n_l  !< maximum speed at the north edge (local)
    REAL(wp)     :: max_uv_n    !< maximum speed at the north edge
    REAL(wp)     :: max_uv      !< maximum speed of all edges

    REAL(wp), DIMENSION(4)                 ::  max_arr    !< 
    REAL(wp), DIMENSION(:), ALLOCATABLE    ::  uv         !< wind velocity at the approach flow
    REAL(wp), DIMENSION(:), ALLOCATABLE    ::  uv_l       !< wind velocity at the approach flow (local)

    REAL(wp), DIMENSION(nzb:nzt+1,nys:nyn) ::  uv_l_nest  !< wind profile at the left edge (nesting)
    REAL(wp), DIMENSION(nzb:nzt+1,nys:nyn) ::  uv_r_nest  !< wind profile at the right edge (nesting)
    REAL(wp), DIMENSION(nzb:nzt+1,nxl:nxr) ::  uv_s_nest  !< wind profile at the south edge (nesting)
    REAL(wp), DIMENSION(nzb:nzt+1,nxl:nxr) ::  uv_n_nest  !< wind profile at the north edge (nesting)


    IF ( .NOT. ALLOCATED( uv ) )  THEN
      ALLOCATE( uv(nzb:nzt+1)   )
      ALLOCATE( uv_l(nzb:nzt+1) )

      uv   = 0.0_wp
      uv_l = 0.0_wp
    ENDIF

!
!-- Determination of the approach flow profile
    IF ( nested_run )  THEN

       uv_l_nest = 0.0_wp
       uv_r_nest = 0.0_wp
       uv_s_nest = 0.0_wp
       uv_n_nest = 0.0_wp
!
!--    Left boundary
       IF ( nxl == 0 )  THEN
          i = nxl
          DO j = nys, nyn
             DO k = nzb, nzt+1
                uv_l_nest(k,j) = SQRT( ( 0.5_wp * ( u(k,j,i-1) + u(k,j,i) ) )**2  &
                                     + ( 0.5_wp * ( v(k,j-1,i) + v(k,j,i) ) )**2  )
             ENDDO
          ENDDO
          max_uv_l_l = MAXVAL(uv_l_nest)
       ENDIF
!
!--    Right boundary
       IF( nxr == nx )  THEN
          i = nxr
          DO j = nys, nyn
             DO k = nzb, nzt+1
                uv_r_nest(k,j) = SQRT( ( 0.5_wp * ( u(k,j,i-1) + u(k,j,i) ) )**2  &
                                     + ( 0.5_wp * ( v(k,j-1,i) + v(k,j,i) ) )**2  )

             ENDDO
          ENDDO
          max_uv_r_l = MAXVAL(uv_r_nest)
       ENDIF
!
!--    South boundary
       IF ( nys == 0 )  THEN
          j = nys
          DO i = nxl, nxr
             DO k = nzb, nzt+1
                uv_s_nest(k,i) = SQRT( ( 0.5_wp * ( u(k,j,i-1) + u(k,j,i) ) )**2  &
                                     + ( 0.5_wp * ( v(k,j-1,i) + v(k,j,i) ) )**2  )
             ENDDO
          ENDDO
          max_uv_s_l = MAXVAL(uv_s_nest)
       ENDIF
!
!--    North boundary
       IF ( nyn == ny )  THEN
          j = nyn
          DO i = nxl, nxr
             DO k = nzb, nzt+1
                uv_n_nest(k,i) = SQRT( ( 0.5_wp * ( u(k,j,i-1) + u(k,j,i) ) )**2  &
                                     + ( 0.5_wp * ( v(k,j-1,i) + v(k,j,i) ) )**2  )

             ENDDO
          ENDDO
          max_uv_n_l = MAXVAL(uv_n_nest)
       ENDIF

#if defined( __parallel )
       CALL MPI_ALLREDUCE( max_uv_l_l, max_uv_l, 1, MPI_REAL, MPI_MAX, comm2d, ierr )
       CALL MPI_ALLREDUCE( max_uv_r_l, max_uv_r, 1, MPI_REAL, MPI_MAX, comm2d, ierr )
       CALL MPI_ALLREDUCE( max_uv_s_l, max_uv_s, 1, MPI_REAL, MPI_MAX, comm2d, ierr )
       CALL MPI_ALLREDUCE( max_uv_n_l, max_uv_n, 1, MPI_REAL, MPI_MAX, comm2d, ierr )
#else
       max_uv_l = max_uv_l_l
       max_uv_r = max_uv_r_l
       max_uv_s = max_uv_s_l
       max_uv_n = max_uv_n_l
#endif

       max_arr = (/ max_uv_r, max_uv_l, max_uv_s, max_uv_n /)
       max_uv = MAXVAl( max_arr )

    ELSE  ! non-nested run

       IF ( bc_lr_cyc  .AND.  bc_ns_cyc )  THEN
          IF ( nxl == 0  .AND.  nys == 0 )  THEN
             DO  k = nzb, nzt+1
                uv_l(k) = SQRT( ( 0.5_wp * ( u(k,0,-1) + u(k,0,0) ) )**2  &
                              + ( 0.5_wp * ( v(k,-1,0) + v(k,0,0) ) )**2  )
             ENDDO
          ENDIF
       ENDIF


       IF ( bc_dirichlet_l )  THEN
          IF ( nxl == 0 .AND. nys == 0 )  THEN
             DO  k = nzb, nzt+1
                uv_l(k) = SQRT( ( 0.5_wp * ( u(k,0,-1) + u(k,0,0) ) )**2  &
                              + ( 0.5_wp * ( v(k,-1,0) + v(k,0,0) ) )**2  )
             ENDDO
          ENDIF

       ELSEIF (bc_dirichlet_r )  THEN
          IF ( nxr == nx .AND. nys == 0 )  THEN
             DO  k = nzb, nzt+1
                uv_l(k) = SQRT( ( 0.5_wp * ( u(k,0,nxr) + u(k,0,nxr+1) ) )**2  &
                              + ( 0.5_wp * ( v(k,-1,nxr) + v(k,0,nxr) ) )**2   )
             ENDDO
          ENDIF
       ENDIF

       IF ( bc_dirichlet_n )  THEN
          IF ( nxl == 0 .AND. nyn == ny )  THEN
             DO  k = nzb, nzt+1
                uv_l(k) = SQRT( ( 0.5_wp * ( u(k,nyn,-1) + u(k,nyn,0) ) )**2  &
                              + ( 0.5_wp * ( v(k,nyn+1,0) + v(k,nyn,0) ) )**2 )
             ENDDO
          ENDIF

       ELSEIF ( bc_dirichlet_s )  THEN
          IF ( nxl == 0 .AND. nys == 0 )  THEN
             DO  k = nzb, nzt+1
                uv_l(k) = SQRT( ( 0.5_wp * ( u(k,0,-1) + u(k,0,0) ) )**2  &
                              + ( 0.5_wp * ( v(k,-1,0) + v(k,0,0) ) )**2  )
             ENDDO
          ENDIF
       ENDIF

#if defined( __parallel )
       CALL MPI_ALLREDUCE( uv_l, uv, nzt+1-nzb+1, MPI_REAL, MPI_MAX, comm2d, ierr )
#else
       uv = uv_l
#endif

       max_uv = MAXVAL( uv )

   ENDIF

!
!-- Test for exceedance the specified limits
    message_string = 'A wind component have a magnitude greater ' //  &
                     'than ten times the maximum wind velocity ' //   &
                     'at the approach flow profile.'

    IF ( MAXVAL( ABS( u ) ) > 10.0_wp * max_uv )  THEN
       CALL message( 'vdi_plausible_values', 'PA0667', 2, 2, myid, 6, 0 )
    ENDIF

    IF ( MAXVAL( ABS( v ) ) > 10.0_wp * max_uv )  THEN
       CALL message( 'vdi_plausible_values', 'PA0667', 2, 2, myid, 6, 0 )
    ENDIF

    IF ( MAXVAL( ABS( w ) ) > 10.0_wp * max_uv )  THEN
       CALL message( 'vdi_plausible_values', 'PA0667', 2, 2, myid, 6, 0 )
    ENDIF

!
!-- Test if the potential temperature lies between 220 K and 330 K
    IF ( MAXVAL( pt ) > 330.0_wp .OR. MAXVAL( pt ) < 220.0_wp )  THEN
       message_string = 'The potential temperature does not lie ' //  &
                        'between 220 K and 330 K.'
       CALL message( 'vdi_plausible_values', 'PA0668', 2, 2, myid, 6, 0 )
    ENDIF

 END SUBROUTINE vdi_plausible_values

 END MODULE vdi_internal_controls
