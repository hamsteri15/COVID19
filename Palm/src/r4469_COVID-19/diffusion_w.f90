!> @file diffusion_w.f90
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
! $Id: diffusion_w.f90 4360 2020-01-07 11:25:50Z suehring $
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
! OpenACC port for SPEC
!
! Revision 1.1  1997/09/12 06:24:11  raasch
! Initial revision
!
!
! Description:
! ------------
!> Diffusion term of the w-component
!------------------------------------------------------------------------------!
 MODULE diffusion_w_mod
 

    PRIVATE
    PUBLIC diffusion_w

    INTERFACE diffusion_w
       MODULE PROCEDURE diffusion_w
       MODULE PROCEDURE diffusion_w_ij
    END INTERFACE diffusion_w

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE diffusion_w

       USE arrays_3d,                                                          &          
           ONLY :  ddzu, ddzw, km, tend, u, v, w, drho_air_zw, rho_air
           
       USE grid_variables,                                                     &     
           ONLY :  ddx, ddy
           
       USE indices,                                                            &            
           ONLY :  nxl, nxr, nyn, nys, nzb, nzt, wall_flags_total_0
           
       USE kinds

       USE surface_mod,                                                        &
           ONLY :  surf_def_v, surf_lsm_v, surf_usm_v

       IMPLICIT NONE

       INTEGER(iwp) ::  i             !< running index x direction
       INTEGER(iwp) ::  j             !< running index y direction
       INTEGER(iwp) ::  k             !< running index z direction
       INTEGER(iwp) ::  l             !< running index of surface type, south- or north-facing wall
       INTEGER(iwp) ::  m             !< running index surface elements
       INTEGER(iwp) ::  surf_e        !< End index of surface elements at (j,i)-gridpoint
       INTEGER(iwp) ::  surf_s        !< Start index of surface elements at (j,i)-gridpoint
       
       REAL(wp) ::  flag              !< flag to mask topography grid points
       REAL(wp) ::  kmxm              !< diffusion coefficient on leftward side of the w-gridbox - interpolated onto xu-y grid
       REAL(wp) ::  kmxp              !<diffusion coefficient on rightward side of the w-gridbox - interpolated onto xu-y grid
       REAL(wp) ::  kmym              !< diffusion coefficient on southward side of the w-gridbox - interpolated onto x-yv grid
       REAL(wp) ::  kmyp              !< diffusion coefficient on northward side of the w-gridbox - interpolated onto x-yv grid
       REAL(wp) ::  mask_west         !< flag to mask vertical wall west of the grid point
       REAL(wp) ::  mask_east         !< flag to mask vertical wall east of the grid point 
       REAL(wp) ::  mask_south        !< flag to mask vertical wall south of the grid point 
       REAL(wp) ::  mask_north        !< flag to mask vertical wall north of the grid point 



       !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(i, j, k, l, m) &
       !$ACC PRIVATE(surf_e, surf_s, flag, kmxm, kmxp, kmym, kmyp) &
       !$ACC PRIVATE(mask_west, mask_east, mask_south, mask_north) &
       !$ACC PRESENT(wall_flags_total_0, km) &
       !$ACC PRESENT(u, v, w) &
       !$ACC PRESENT(ddzu, ddzw, rho_air, drho_air_zw) &
       !$ACC PRESENT(surf_def_v(0:3)) &
       !$ACC PRESENT(surf_lsm_v(0:3)) &
       !$ACC PRESENT(surf_usm_v(0:3)) &
       !$ACC PRESENT(tend)
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt-1
!
!--             Predetermine flag to mask topography and wall-bounded grid points. 
                flag       = MERGE( 1.0_wp, 0.0_wp,                            &
                                    BTEST( wall_flags_total_0(k,j,i),   3 ) ) 
                mask_east  = MERGE( 1.0_wp, 0.0_wp,                            &
                                    BTEST( wall_flags_total_0(k,j,i+1), 3 ) )
                mask_west  = MERGE( 1.0_wp, 0.0_wp,                            &
                                    BTEST( wall_flags_total_0(k,j,i-1), 3 ) )
                mask_south = MERGE( 1.0_wp, 0.0_wp,                            &
                                    BTEST( wall_flags_total_0(k,j-1,i), 3 ) )
                mask_north = MERGE( 1.0_wp, 0.0_wp,                            &
                                    BTEST( wall_flags_total_0(k,j+1,i), 3 ) )
!
!--             Interpolate eddy diffusivities on staggered gridpoints
                kmxp = 0.25_wp * ( km(k,j,i)   + km(k,j,i+1)   +               &
                                   km(k+1,j,i) + km(k+1,j,i+1) )
                kmxm = 0.25_wp * ( km(k,j,i)   + km(k,j,i-1)   +               &
                                   km(k+1,j,i) + km(k+1,j,i-1) )
                kmyp = 0.25_wp * ( km(k,j,i)   + km(k+1,j,i)   +               &
                                   km(k,j+1,i) + km(k+1,j+1,i) )
                kmym = 0.25_wp * ( km(k,j,i)   + km(k+1,j,i)   +               &
                                   km(k,j-1,i) + km(k+1,j-1,i) )

                tend(k,j,i) = tend(k,j,i)                                      &
                       + ( mask_east *  kmxp * (                               &
                                   ( w(k,j,i+1)   - w(k,j,i)   ) * ddx         &
                                 + ( u(k+1,j,i+1) - u(k,j,i+1) ) * ddzu(k+1)   &
                                               )                               &
                         - mask_west * kmxm *  (                               &
                                   ( w(k,j,i)     - w(k,j,i-1) ) * ddx         &
                                 + ( u(k+1,j,i)   - u(k,j,i)   ) * ddzu(k+1)   &
                                               )                               &
                         ) * ddx                                 * flag        &
                       + ( mask_north * kmyp * (                               &
                                   ( w(k,j+1,i)   - w(k,j,i)   ) * ddy         &
                                 + ( v(k+1,j+1,i) - v(k,j+1,i) ) * ddzu(k+1)   &
                                               )                               &
                         - mask_south * kmym * (                               &
                                   ( w(k,j,i)     - w(k,j-1,i) ) * ddy         &
                                 + ( v(k+1,j,i)   - v(k,j,i)   ) * ddzu(k+1)   &
                                               )                               &
                         ) * ddy                                 * flag        &
                       + 2.0_wp * (                                            &
                         km(k+1,j,i) * ( w(k+1,j,i) - w(k,j,i) )   * ddzw(k+1) &
                                     * rho_air(k+1)                            &
                       - km(k,j,i)   * ( w(k,j,i)   - w(k-1,j,i) ) * ddzw(k)   &
                                     * rho_air(k)                              &
                                  ) * ddzu(k+1) * drho_air_zw(k) * flag
             ENDDO

!
!--          Add horizontal momentum flux v'w' at north- (l=0) and south-facing (l=1)
!--          surfaces. Note, in the the flat case, loops won't be entered as 
!--          start_index > end_index. Furtermore, note, no vertical natural surfaces
!--          so far.           
!--          Default-type surfaces
             DO  l = 0, 1
                surf_s = surf_def_v(l)%start_index(j,i)
                surf_e = surf_def_v(l)%end_index(j,i)
                DO  m = surf_s, surf_e
                   k           = surf_def_v(l)%k(m)
                   tend(k,j,i) = tend(k,j,i) +                                 &
                                     surf_def_v(l)%mom_flux_w(m) * ddy
                ENDDO   
             ENDDO
!
!--          Natural-type surfaces
             DO  l = 0, 1
                surf_s = surf_lsm_v(l)%start_index(j,i)
                surf_e = surf_lsm_v(l)%end_index(j,i)
                DO  m = surf_s, surf_e
                   k           = surf_lsm_v(l)%k(m)
                   tend(k,j,i) = tend(k,j,i) +                                 &
                                     surf_lsm_v(l)%mom_flux_w(m) * ddy
                ENDDO   
             ENDDO
!
!--          Urban-type surfaces
             DO  l = 0, 1
                surf_s = surf_usm_v(l)%start_index(j,i)
                surf_e = surf_usm_v(l)%end_index(j,i)
                DO  m = surf_s, surf_e
                   k           = surf_usm_v(l)%k(m)
                   tend(k,j,i) = tend(k,j,i) +                                 &
                                     surf_usm_v(l)%mom_flux_w(m) * ddy
                ENDDO   
             ENDDO
!
!--          Add horizontal momentum flux u'w' at east- (l=2) and west-facing (l=3)
!--          surface.
!--          Default-type surfaces
             DO  l = 2, 3
                surf_s = surf_def_v(l)%start_index(j,i)
                surf_e = surf_def_v(l)%end_index(j,i)
                DO  m = surf_s, surf_e
                   k           = surf_def_v(l)%k(m)
                   tend(k,j,i) = tend(k,j,i) +                                 &
                                     surf_def_v(l)%mom_flux_w(m) * ddx
                ENDDO   
             ENDDO
!
!--          Natural-type surfaces
             DO  l = 2, 3
                surf_s = surf_lsm_v(l)%start_index(j,i)
                surf_e = surf_lsm_v(l)%end_index(j,i)
                DO  m = surf_s, surf_e
                   k           = surf_lsm_v(l)%k(m)
                   tend(k,j,i) = tend(k,j,i) +                                 &
                                     surf_lsm_v(l)%mom_flux_w(m) * ddx
                ENDDO   
             ENDDO
!
!--          Urban-type surfaces
             DO  l = 2, 3
                surf_s = surf_usm_v(l)%start_index(j,i)
                surf_e = surf_usm_v(l)%end_index(j,i)
                DO  m = surf_s, surf_e
                   k           = surf_usm_v(l)%k(m)
                   tend(k,j,i) = tend(k,j,i) +                                 &
                                     surf_usm_v(l)%mom_flux_w(m) * ddx
                ENDDO   
             ENDDO

          ENDDO
       ENDDO

    END SUBROUTINE diffusion_w


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE diffusion_w_ij( i, j )

       USE arrays_3d,                                                          &          
           ONLY :  ddzu, ddzw, km, tend, u, v, w, drho_air_zw, rho_air
           
       USE grid_variables,                                                     &     
           ONLY :  ddx, ddy
           
       USE indices,                                                            &            
           ONLY :  nzb, nzt, wall_flags_total_0
           
       USE kinds

       USE surface_mod,                                                        &
           ONLY :  surf_def_v, surf_lsm_v, surf_usm_v

       IMPLICIT NONE


       INTEGER(iwp) ::  i             !< running index x direction
       INTEGER(iwp) ::  j             !< running index y direction
       INTEGER(iwp) ::  k             !< running index z direction
       INTEGER(iwp) ::  l             !< running index of surface type, south- or north-facing wall
       INTEGER(iwp) ::  m             !< running index surface elements
       INTEGER(iwp) ::  surf_e        !< End index of surface elements at (j,i)-gridpoint
       INTEGER(iwp) ::  surf_s        !< Start index of surface elements at (j,i)-gridpoint
       
       REAL(wp) ::  flag              !< flag to mask topography grid points
       REAL(wp) ::  kmxm              !< diffusion coefficient on leftward side of the w-gridbox - interpolated onto xu-y grid
       REAL(wp) ::  kmxp              !< diffusion coefficient on rightward side of the w-gridbox - interpolated onto xu-y grid
       REAL(wp) ::  kmym              !< diffusion coefficient on southward side of the w-gridbox - interpolated onto x-yv grid
       REAL(wp) ::  kmyp              !< diffusion coefficient on northward side of the w-gridbox - interpolated onto x-yv grid
       REAL(wp) ::  mask_west         !< flag to mask vertical wall west of the grid point
       REAL(wp) ::  mask_east         !< flag to mask vertical wall east of the grid point 
       REAL(wp) ::  mask_south        !< flag to mask vertical wall south of the grid point 
       REAL(wp) ::  mask_north        !< flag to mask vertical wall north of the grid point 


       DO  k = nzb+1, nzt-1
!
!--       Predetermine flag to mask topography and wall-bounded grid points. 
          flag       = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i),   3 ) ) 
          mask_east  = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i+1), 3 ) )
          mask_west  = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j,i-1), 3 ) )
          mask_south = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j-1,i), 3 ) )
          mask_north = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_total_0(k,j+1,i), 3 ) )
!
!--       Interpolate eddy diffusivities on staggered gridpoints
          kmxp = 0.25_wp * ( km(k,j,i)+km(k,j,i+1)+km(k+1,j,i)+km(k+1,j,i+1) )
          kmxm = 0.25_wp * ( km(k,j,i)+km(k,j,i-1)+km(k+1,j,i)+km(k+1,j,i-1) )
          kmyp = 0.25_wp * ( km(k,j,i)+km(k+1,j,i)+km(k,j+1,i)+km(k+1,j+1,i) )
          kmym = 0.25_wp * ( km(k,j,i)+km(k+1,j,i)+km(k,j-1,i)+km(k+1,j-1,i) )

          tend(k,j,i) = tend(k,j,i)                                            &
                       + ( mask_east *  kmxp * (                               &
                                   ( w(k,j,i+1)   - w(k,j,i)   ) * ddx         &
                                 + ( u(k+1,j,i+1) - u(k,j,i+1) ) * ddzu(k+1)   &
                                               )                               &
                         - mask_west * kmxm *  (                               &
                                   ( w(k,j,i)     - w(k,j,i-1) ) * ddx         &
                                 + ( u(k+1,j,i)   - u(k,j,i)   ) * ddzu(k+1)   &
                                               )                               &
                         ) * ddx                                 * flag        &
                       + ( mask_north * kmyp * (                               &
                                   ( w(k,j+1,i)   - w(k,j,i)   ) * ddy         &
                                 + ( v(k+1,j+1,i) - v(k,j+1,i) ) * ddzu(k+1)   &
                                               )                               &
                         - mask_south * kmym * (                               &
                                   ( w(k,j,i)     - w(k,j-1,i) ) * ddy         &
                                 + ( v(k+1,j,i)   - v(k,j,i)   ) * ddzu(k+1)   &
                                               )                               &
                         ) * ddy                                 * flag        &
                       + 2.0_wp * (                                            &
                         km(k+1,j,i) * ( w(k+1,j,i) - w(k,j,i) ) * ddzw(k+1)   &
                                     * rho_air(k+1)                            &
                       - km(k,j,i)   * ( w(k,j,i)   - w(k-1,j,i) ) * ddzw(k)   &
                                     * rho_air(k)                              &
                                  ) * ddzu(k+1) * drho_air_zw(k) * flag
       ENDDO
!
!--    Add horizontal momentum flux v'w' at north- (l=0) and south-facing (l=1)
!--    surfaces. Note, in the the flat case, loops won't be entered as 
!--    start_index > end_index. Furtermore, note, no vertical natural surfaces
!--    so far.           
!--    Default-type surfaces
       DO  l = 0, 1
          surf_s = surf_def_v(l)%start_index(j,i)
          surf_e = surf_def_v(l)%end_index(j,i)
          DO  m = surf_s, surf_e
             k           = surf_def_v(l)%k(m)
             tend(k,j,i) = tend(k,j,i) +                                       &
                                     surf_def_v(l)%mom_flux_w(m) * ddy
          ENDDO   
       ENDDO
!
!--    Natural-type surfaces
       DO  l = 0, 1
          surf_s = surf_lsm_v(l)%start_index(j,i)
          surf_e = surf_lsm_v(l)%end_index(j,i)
          DO  m = surf_s, surf_e
             k           = surf_lsm_v(l)%k(m)
             tend(k,j,i) = tend(k,j,i) +                                       &
                                     surf_lsm_v(l)%mom_flux_w(m) * ddy
          ENDDO   
       ENDDO
!
!--    Urban-type surfaces
       DO  l = 0, 1
          surf_s = surf_usm_v(l)%start_index(j,i)
          surf_e = surf_usm_v(l)%end_index(j,i)
          DO  m = surf_s, surf_e
             k           = surf_usm_v(l)%k(m)
             tend(k,j,i) = tend(k,j,i) +                                       &
                                     surf_usm_v(l)%mom_flux_w(m) * ddy
          ENDDO   
       ENDDO
!
!--    Add horizontal momentum flux u'w' at east- (l=2) and west-facing (l=3)
!--    surfaces. 
!--    Default-type surfaces
       DO  l = 2, 3
          surf_s = surf_def_v(l)%start_index(j,i)
          surf_e = surf_def_v(l)%end_index(j,i)
          DO  m = surf_s, surf_e
             k           = surf_def_v(l)%k(m)
             tend(k,j,i) = tend(k,j,i) +                                       &
                                     surf_def_v(l)%mom_flux_w(m) * ddx
          ENDDO   
       ENDDO
!
!--    Natural-type surfaces
       DO  l = 2, 3
          surf_s = surf_lsm_v(l)%start_index(j,i)
          surf_e = surf_lsm_v(l)%end_index(j,i)
          DO  m = surf_s, surf_e
             k           = surf_lsm_v(l)%k(m)
             tend(k,j,i) = tend(k,j,i) +                                       &
                                     surf_lsm_v(l)%mom_flux_w(m) * ddx
          ENDDO   
       ENDDO
!
!--    Urban-type surfaces
       DO  l = 2, 3
          surf_s = surf_usm_v(l)%start_index(j,i)
          surf_e = surf_usm_v(l)%end_index(j,i)
          DO  m = surf_s, surf_e
             k           = surf_usm_v(l)%k(m)
             tend(k,j,i) = tend(k,j,i) +                                       &
                                     surf_usm_v(l)%mom_flux_w(m) * ddx
          ENDDO   
       ENDDO


    END SUBROUTINE diffusion_w_ij

 END MODULE diffusion_w_mod
