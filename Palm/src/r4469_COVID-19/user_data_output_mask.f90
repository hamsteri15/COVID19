!> @file user_data_output_mask.f90
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
! $Id: user_data_output_mask.f90 4360 2020-01-07 11:25:50Z suehring $
! Corrected "Former revisions" section
! 
! 4168 2019-08-16 13:50:17Z suehring
! Remove dependency on surface_mod + example for terrain-following output 
! adjusted
! 
! 4069 2019-07-01 14:05:51Z Giersch
! Masked output running index mid has been introduced as a local variable to 
! avoid runtime error (Loop variable has been modified) in time_integration
! 
! 3768 2019-02-27 14:35:58Z raasch
! variables commented + statement added to avoid compiler warnings about unused variables
! 
! 3655 2019-01-07 16:51:22Z knoop
! Add terrain-following output
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! Description:
! ------------
!> Resorts the user-defined output quantity with indices (k,j,i) to a
!> temporary array with indices (i,j,k) for masked data output.
!------------------------------------------------------------------------------!
 SUBROUTINE user_data_output_mask( av, variable, found, local_pf, mid )
 

    USE control_parameters
        
    USE indices
    
    USE kinds
    
    USE user

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  variable  !< 

    INTEGER(iwp) ::  av             !< 
    INTEGER(iwp) ::  mid            !< masked output running index
!    INTEGER(iwp) ::  i              !<
!    INTEGER(iwp) ::  j              !<
!    INTEGER(iwp) ::  k              !<
!    INTEGER(iwp) ::  topo_top_index !< k index of highest horizontal surface

    LOGICAL ::  found               !< 

    REAL(wp),                                                                  &
       DIMENSION(mask_size_l(mid,1),mask_size_l(mid,2),mask_size_l(mid,3)) ::  &
          local_pf   !< 

!
!-- Next line is to avoid compiler warning about unused variables. Please remove.
    IF ( av == 0  .OR.                                                                             &
         local_pf(mask_size_l(mid,1),mask_size_l(mid,2),mask_size_l(mid,3)) == 0.0_wp )  CONTINUE


    found = .TRUE.

    SELECT CASE ( TRIM( variable ) )

!--    Uncomment and extend the following lines, if necessary.
!--    The arrays for storing the user defined quantities (here u2 and u2_av)
!--    have to be declared and defined by the user!
!--    Sample for user-defined output:
!       CASE ( 'u2' )
!          IF ( av == 0 )  THEN
!             IF ( .NOT. mask_surface(mid) )  THEN
!!
!!--             Default masked output
!                DO  i = 1, mask_size_l(mid,1)
!                   DO  j = 1, mask_size_l(mid,2)
!                      DO  k = 1, mask_size_l(mid,3)
!                         local_pf(i,j,k) = u2(mask_k(mid,k),                  &
!                                              mask_j(mid,j),                  &
!                                              mask_i(mid,i))
!                      ENDDO
!                   ENDDO
!                ENDDO
!             ELSE
!!
!!--             Terrain-following masked output
!                DO  i = 1, mask_size_l(mid,1)
!                   DO  j = 1, mask_size_l(mid,2)
!!
!!--                   Get k index of highest horizontal surface
!                      topo_top_index = topo_top_ind( &
!                                        mask_j(mid,j), &
!                                        mask_i(mid,i), &
!                                        1          )
!!
!!--                   Save output array
!                      DO  k = 1, mask_size_l(mid,3)
!                         local_pf(i,j,k) = u2(MIN( topo_top_index+mask_k(mid,k),&
!                                                   nzt+1 ),                     &
!                                              mask_j(mid,j),                    &
!                                              mask_i(mid,i)                   )
!                      ENDDO
!                   ENDDO
!                ENDDO
!             ENDIF
!          ELSE
!             IF ( .NOT. mask_surface(mid) )  THEN
!!
!!--             Default masked output
!                DO  i = 1, mask_size_l(mid,1)
!                   DO  j = 1, mask_size_l(mid,2)
!                      DO  k = 1, mask_size_l(mid,3)
!                          local_pf(i,j,k) = u2_av(mask_k(mid,k),              &
!                                                  mask_j(mid,j),              &
!                                                  mask_i(mid,i) )
!                       ENDDO
!                    ENDDO
!                 ENDDO
!             ELSE
!!
!!--             Terrain-following masked output
!                DO  i = 1, mask_size_l(mid,1)
!                   DO  j = 1, mask_size_l(mid,2)
!!
!!--                   Get k index of highest horizontal surface
!                      topo_top_index = topo_top_ind(   &
!                                        mask_j(mid,j), &
!                                        mask_i(mid,i), &
!                                        1 )
!!
!!--                   Save output array
!                      DO  k = 1, mask_size_l(mid,3)
!                         local_pf(i,j,k) = u2_av(                               &
!                                              MIN( topo_top_index+mask_k(mid,k),&
!                                                   nzt+1 ),                     &
!                                              mask_j(mid,j),                    &
!                                              mask_i(mid,i)                   )
!                      ENDDO
!                   ENDDO
!                ENDDO
!             ENDIF
!          ENDIF

       CASE DEFAULT
          found = .FALSE.

    END SELECT


 END SUBROUTINE user_data_output_mask
