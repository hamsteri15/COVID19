!> @file buoyancy.f90
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
! $Id: buoyancy.f90 4360 2020-01-07 11:25:50Z suehring $
! Corrected "Former revisions" section
! 
! 3725 2019-02-07 10:11:02Z raasch
! unused variables removed
! 
! 3655 2019-01-07 16:51:22Z knoop
! nopointer option removed
!
! Revision 1.1  1997/08/29 08:56:48  raasch
! Initial revision
!
!
! Description:
! ------------
!> Buoyancy term of the third component of the equation of motion.
!> @attention Humidity is not regarded when using a sloping surface!
!------------------------------------------------------------------------------!
 MODULE buoyancy_mod

    USE basic_constants_and_equations_mod,                                     &
        ONLY:  g

    USE arrays_3d,                                                             &
        ONLY:  pt, pt_slope_ref, ref_state, tend

    USE control_parameters,                                                    &
        ONLY:  atmos_ocean_sign, cos_alpha_surface, message_string, pt_surface,&
               sin_alpha_surface, sloping_surface

    USE kinds

    PRIVATE
    PUBLIC buoyancy

    INTERFACE buoyancy
       MODULE PROCEDURE buoyancy
       MODULE PROCEDURE buoyancy_ij
    END INTERFACE buoyancy

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE buoyancy( var, wind_component )

       USE indices,                                                            &
           ONLY:  nxl, nxlu, nxr, nyn, nys, nzb, nzt


       IMPLICIT NONE

       INTEGER(iwp) ::  i              !<
       INTEGER(iwp) ::  j              !<
       INTEGER(iwp) ::  k              !<
       INTEGER(iwp) ::  wind_component !<
       
       REAL(wp), DIMENSION(:,:,:), POINTER ::  var


       IF ( .NOT. sloping_surface )  THEN
!
!--       Normal case: horizontal surface
          !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(i, j, k) &
          !$ACC PRESENT(var) &
          !$ACC PRESENT(ref_state) &
          !$ACC PRESENT(tend)
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt-1
                   tend(k,j,i) = tend(k,j,i) + atmos_ocean_sign * g * 0.5_wp *  &
                          (                                                     &
                             ( var(k,j,i)   - ref_state(k) )   / ref_state(k) + &
                             ( var(k+1,j,i) - ref_state(k+1) ) / ref_state(k+1) &
                          )
                ENDDO
             ENDDO
          ENDDO

       ELSE
!
!--       Buoyancy term for a surface with a slope in x-direction. The equations
!--       for both the u and w velocity-component contain proportionate terms.
!--       Temperature field at time t=0 serves as environmental temperature.
!--       Reference temperature (pt_surface) is the one at the lower left corner
!--       of the total domain.
          IF ( wind_component == 1 )  THEN

             DO  i = nxlu, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt-1
                      tend(k,j,i) = tend(k,j,i) + g * sin_alpha_surface *         &
                           0.5_wp * ( ( pt(k,j,i-1)         + pt(k,j,i)         ) &
                                    - ( pt_slope_ref(k,i-1) + pt_slope_ref(k,i) ) &
                                    ) / pt_surface
                   ENDDO
                ENDDO
             ENDDO

          ELSEIF ( wind_component == 3 )  THEN

             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt-1
                      tend(k,j,i) = tend(k,j,i) + g * cos_alpha_surface *         &
                           0.5_wp * ( ( pt(k,j,i)         + pt(k+1,j,i)         ) &
                                    - ( pt_slope_ref(k,i) + pt_slope_ref(k+1,i) ) &
                                    ) / pt_surface
                   ENDDO
                ENDDO
            ENDDO

          ELSE
             
             WRITE( message_string, * ) 'no term for component "',             &
                                       wind_component,'"'
             CALL message( 'buoyancy', 'PA0159', 1, 2, 0, 6, 0 )

          ENDIF

       ENDIF

    END SUBROUTINE buoyancy


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!> @attention PGI-compiler creates SIGFPE if opt>1 is used! Therefore, opt=1 is
!>            forced by compiler-directive.
!------------------------------------------------------------------------------!
!pgi$r opt=1
    SUBROUTINE buoyancy_ij( i, j, var, wind_component )


       USE indices,                                                            &
           ONLY:  nzb, nzt

       IMPLICIT NONE

       INTEGER(iwp) ::  i              !<
       INTEGER(iwp) ::  j              !<
       INTEGER(iwp) ::  k              !<
       INTEGER(iwp) ::  wind_component !<
       
       REAL(wp), DIMENSION(:,:,:), POINTER ::  var


       IF ( .NOT. sloping_surface )  THEN
!
!--       Normal case: horizontal surface
          DO  k = nzb+1, nzt-1
              tend(k,j,i) = tend(k,j,i) + atmos_ocean_sign * g * 0.5_wp * (    &
                        ( var(k,j,i)   - ref_state(k)   ) / ref_state(k)   +   &
                        ( var(k+1,j,i) - ref_state(k+1) ) / ref_state(k+1)     &
                                                                          )
          ENDDO

       ELSE
!
!--       Buoyancy term for a surface with a slope in x-direction. The equations
!--       for both the u and w velocity-component contain proportionate terms.
!--       Temperature field at time t=0 serves as environmental temperature.
!--       Reference temperature (pt_surface) is the one at the lower left corner
!--       of the total domain.
          IF ( wind_component == 1 )  THEN

             DO  k = nzb+1, nzt-1
                tend(k,j,i) = tend(k,j,i) + g * sin_alpha_surface *               &
                           0.5_wp * ( ( pt(k,j,i-1)         + pt(k,j,i)         ) &
                                    - ( pt_slope_ref(k,i-1) + pt_slope_ref(k,i) ) &
                                    ) / pt_surface
             ENDDO

          ELSEIF ( wind_component == 3 )  THEN

             DO  k = nzb+1, nzt-1
                tend(k,j,i) = tend(k,j,i) + g * cos_alpha_surface *               &
                           0.5_wp * ( ( pt(k,j,i)         + pt(k+1,j,i)         ) &
                                    - ( pt_slope_ref(k,i) + pt_slope_ref(k+1,i) ) &
                                    ) / pt_surface
             ENDDO

          ELSE

             WRITE( message_string, * ) 'no term for component "',             &
                                       wind_component,'"'
             CALL message( 'buoyancy', 'PA0159', 1, 2, 0, 6, 0 )

          ENDIF

       ENDIF

    END SUBROUTINE buoyancy_ij

 END MODULE buoyancy_mod
