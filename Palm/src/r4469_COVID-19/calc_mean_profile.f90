!> @file calc_mean_profile.f90
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
! $Id: calc_mean_profile.f90 4360 2020-01-07 11:25:50Z suehring $
! Introduction of wall_flags_total_0, which currently sets bits based on static
! topography information used in wall_flags_static_0
! 
! 4329 2019-12-10 15:46:36Z motisi
! Renamed wall_flags_0 to wall_flags_static_0
! 
! 4182 2019-08-22 15:20:23Z scharf
! Corrected "Former revisions" section
! 
! 3655 2019-01-07 16:51:22Z knoop
! nopointer option removed
! 
! 1365 2014-04-22 15:03:56Z boeske
! Initial revision
!
! Description:
! ------------
!> Calculate the horizontally averaged vertical temperature profile (pr=4 in case
!> of potential temperature, 44 in case of virtual potential temperature, and 64
!> in case of density (ocean runs)).
!------------------------------------------------------------------------------!
 MODULE calc_mean_profile_mod
 

    PRIVATE
    PUBLIC calc_mean_profile

    INTERFACE calc_mean_profile
       MODULE PROCEDURE calc_mean_profile
    END INTERFACE calc_mean_profile

 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE calc_mean_profile( var, pr )

       USE control_parameters,                                                 &
           ONLY:  intermediate_timestep_count

       USE indices,                                                            &
           ONLY:  ngp_2dh_s_inner, nxl, nxr, nyn, nys, nzb, nzb, nzt,          &
                  wall_flags_total_0

       USE kinds

       USE pegrid

       USE statistics,                                                         &
           ONLY:  flow_statistics_called, hom, sums, sums_l


       IMPLICIT NONE
       
       INTEGER(iwp) ::  i                  !<
       INTEGER(iwp) ::  j                  !<
       INTEGER(iwp) ::  k                  !<
       INTEGER(iwp) ::  pr                 !< 
!$     INTEGER(iwp) ::  omp_get_thread_num !<
       INTEGER(iwp) ::  tn                 !<
       
       REAL(wp), DIMENSION(:,:,:), POINTER ::  var

!
!--    Computation of the horizontally averaged profile of variable var, unless
!--    already done by the relevant call from flow_statistics. The calculation
!--    is done only for the first respective intermediate timestep in order to
!--    spare communication time and to produce identical model results with jobs
!--    which are calling flow_statistics at different time intervals. At 
!--    initialization, intermediate_timestep_count = 0 is considered as well.

       IF ( .NOT. flow_statistics_called  .AND.                                &
            intermediate_timestep_count <= 1 )  THEN

!
!--       Horizontal average of variable var
          tn           =   0  ! Default thread number in case of one thread
          !$OMP PARALLEL PRIVATE( i, j, k, tn )
!$        tn = omp_get_thread_num()
          sums_l(:,pr,tn) = 0.0_wp
          !$OMP DO
          DO  i = nxl, nxr
             DO  j =  nys, nyn
                DO  k = nzb, nzt+1
                   sums_l(k,pr,tn) = sums_l(k,pr,tn) + var(k,j,i)              &
                                     * MERGE( 1.0_wp, 0.0_wp,                  &
                                          BTEST( wall_flags_total_0(k,j,i), 22 ) )
                ENDDO
             ENDDO
          ENDDO
          !$OMP END PARALLEL

          DO  i = 1, threads_per_task-1
             sums_l(:,pr,0) = sums_l(:,pr,0) + sums_l(:,pr,i)
          ENDDO

#if defined( __parallel )

          IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
          CALL MPI_ALLREDUCE( sums_l(nzb,pr,0), sums(nzb,pr), nzt+2-nzb,       &
                              MPI_REAL, MPI_SUM, comm2d, ierr )

#else

          sums(:,pr) = sums_l(:,pr,0)

#endif

          DO  k = nzb, nzt+1
             IF ( ngp_2dh_s_inner(k,0) /= 0 )  THEN
                hom(k,1,pr,0) = sums(k,pr) / ngp_2dh_s_inner(k,0)
             ENDIF
          ENDDO

       ENDIF


    END SUBROUTINE calc_mean_profile

 END MODULE calc_mean_profile_mod
