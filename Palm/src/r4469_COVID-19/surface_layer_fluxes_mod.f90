!> @file surface_layer_fluxes_mod.f90
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
!
!------------------------------------------------------------------------------!
!
! Current revisions:
! ------------------
! 
! 
! Former revisions:
! -----------------
! $Id: surface_layer_fluxes_mod.f90 4370 2020-01-10 14:00:44Z raasch $
! bugfix: openacc porting for vector version of OL calculation added
! 
! 4366 2020-01-09 08:12:43Z raasch
! vector version for calculation of Obukhov length via Newton iteration added
! 
! 4360 2020-01-07 11:25:50Z suehring
! Calculation of diagnostic-only 2-m potential temperature moved to 
! diagnostic_output_quantities.
! 
! 4298 2019-11-21 15:59:16Z suehring
! Calculation of 2-m temperature adjusted to the case the 2-m level is above 
! the first grid point.
! 
! 4258 2019-10-07 13:29:08Z suehring
! Initialization of Obukhov lenght also at vertical surfaces (if allocated).
! 
! 4237 2019-09-25 11:33:42Z knoop
! Added missing OpenMP directives
! 
! 4186 2019-08-23 16:06:14Z suehring
! - To enable limitation of Obukhov length, move it before exit-cycle construct. 
!   Further, change the limit to 10E-5 in order to get rid-off unrealistic 
!   peaks in the heat fluxes during nighttime
! - Unused variable removed
! 
! 4182 2019-08-22 15:20:23Z scharf
! Corrected "Former revisions" section
! 
! 3987 2019-05-22 09:52:13Z kanani
! Introduce alternative switch for debug output during timestepping
! 
! 3885 2019-04-11 11:29:34Z kanani
! Changes related to global restructuring of location messages and introduction 
! of additional debug messages
! 
! 3881 2019-04-10 09:31:22Z suehring
! Assure that Obukhov length does not become zero
! 
! 3834 2019-03-28 15:40:15Z forkel
! added USE chem_gasphase_mod 
! 
! 3787 2019-03-07 08:43:54Z raasch
! unused variables removed
! 
! 3745 2019-02-15 18:57:56Z suehring
! Bugfix, missing calculation of 10cm temperature at vertical building walls, 
! required for indoor model
! 
! 3744 2019-02-15 18:38:58Z suehring
! Some interface calls moved to module_interface + cleanup
! 
! 3668 2019-01-14 12:49:24Z maronga
! Removed methods "circular" and "lookup"
! 
! 3655 2019-01-07 16:51:22Z knoop
! OpenACC port for SPEC
!
! Revision 1.1  1998/01/23 10:06:06  raasch
! Initial revision
!
!
! Description:
! ------------
!> Diagnostic computation of vertical fluxes in the constant flux layer from the
!> values of the variables at grid point k=1 based on Newton iteration
!>
!> @todo (re)move large_scale_forcing actions
!> @todo check/optimize OpenMP directives
!> @todo simplify if conditions (which flux need to be computed in which case)
!------------------------------------------------------------------------------!
 MODULE surface_layer_fluxes_mod

    USE arrays_3d,                                                             &
        ONLY:  e, kh, nc, nr, pt, q, ql, qc, qr, s, u, v, vpt, w, zu, zw,      &
               drho_air_zw, rho_air_zw, d_exner

    USE basic_constants_and_equations_mod,                                     &
        ONLY:  g, kappa, lv_d_cp, pi, rd_d_rv

    USE chem_gasphase_mod,                                                     &
        ONLY:  nvar

    USE chem_modules,                                                          &
        ONLY:  constant_csflux

    USE cpulog

    USE control_parameters,                                                    &
        ONLY:  air_chemistry, cloud_droplets,                                  &
               constant_heatflux, constant_scalarflux,                         &
               constant_waterflux, coupling_mode,                              &
               debug_output_timestep,                                          &
               humidity, loop_optimization,                                    &
               ibc_e_b, ibc_pt_b, indoor_model,                                &
               land_surface, large_scale_forcing, lsf_surf, message_string,    &
               neutral, passive_scalar, pt_surface, q_surface,                 &
               run_coupled, surface_pressure, simulated_time,                  &
               time_since_reference_point, urban_surface,                      &
               use_free_convection_scaling, zeta_max, zeta_min

    USE grid_variables,                                                        &
        ONLY:  dx, dy  

    USE indices,                                                               &
        ONLY:  nzt

    USE kinds

    USE bulk_cloud_model_mod,                                                  &
        ONLY: bulk_cloud_model, microphysics_morrison, microphysics_seifert

    USE pegrid

    USE land_surface_model_mod,                                                &
        ONLY:  aero_resist_kray, skip_time_do_lsm

    USE surface_mod,                                                           &
        ONLY :  surf_def_h, surf_def_v, surf_lsm_h, surf_lsm_v, surf_type,     &
                surf_usm_h, surf_usm_v
        

    IMPLICIT NONE

    INTEGER(iwp) ::  i              !< loop index x direction
    INTEGER(iwp) ::  j              !< loop index y direction
    INTEGER(iwp) ::  k              !< loop index z direction
    INTEGER(iwp) ::  l              !< loop index for surf type

    LOGICAL      ::  coupled_run       !< Flag for coupled atmosphere-ocean runs
    LOGICAL      ::  downward = .FALSE.!< Flag indicating downward-facing horizontal surface
    LOGICAL      ::  mom_uv  = .FALSE. !< Flag indicating calculation of usvs and vsus at vertical surfaces
    LOGICAL      ::  mom_w   = .FALSE. !< Flag indicating calculation of wsus and wsvs at vertical surfaces
    LOGICAL      ::  mom_tke = .FALSE. !< Flag indicating calculation of momentum fluxes at vertical surfaces used for TKE production 
    LOGICAL      ::  surf_vertical     !< Flag indicating vertical surfaces

    REAL(wp)     ::  e_s,               & !< Saturation water vapor pressure
                     ol_max = 1.0E6_wp, & !< Maximum Obukhov length
                     z_mo                 !< Height of the constant flux layer where MOST is assumed

    TYPE(surf_type), POINTER ::  surf     !< surf-type array, used to generalize subroutines


    SAVE

    PRIVATE

    PUBLIC init_surface_layer_fluxes,                                          &
           phi_m,                                                              &
           psi_h,                                                              &
           psi_m,                                                              &
           surface_layer_fluxes

    INTERFACE init_surface_layer_fluxes
       MODULE PROCEDURE init_surface_layer_fluxes
    END INTERFACE init_surface_layer_fluxes

    INTERFACE phi_m
       MODULE PROCEDURE phi_m
    END INTERFACE phi_m

    INTERFACE psi_h
       MODULE PROCEDURE psi_h
    END INTERFACE psi_h

    INTERFACE psi_m
       MODULE PROCEDURE psi_m
    END INTERFACE psi_m

    INTERFACE surface_layer_fluxes
       MODULE PROCEDURE surface_layer_fluxes
    END INTERFACE surface_layer_fluxes


 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Main routine to compute the surface fluxes
!------------------------------------------------------------------------------!
    SUBROUTINE surface_layer_fluxes

       IMPLICIT NONE


       IF ( debug_output_timestep )  CALL debug_message( 'surface_layer_fluxes', 'start' )

       surf_vertical = .FALSE. !< flag indicating vertically orientated surface elements
       downward      = .FALSE. !< flag indicating downward-facing surface elements
!
!--    Derive potential temperature and specific humidity at first grid level 
!--    from the fields pt and q
!
!--    First call for horizontal default-type surfaces (l=0 - upward facing, 
!--    l=1 - downward facing)
       DO  l = 0, 1
          IF ( surf_def_h(l)%ns >= 1  )  THEN
             surf => surf_def_h(l)
             CALL calc_pt_q
             IF ( .NOT. neutral )  THEN
                CALL calc_pt_surface
                IF ( humidity )  THEN
                   CALL calc_q_surface
                   CALL calc_vpt_surface
                ENDIF
             ENDIF
          ENDIF
       ENDDO
!
!--    Call for natural-type horizontal surfaces
       IF ( surf_lsm_h%ns >= 1  )  THEN
          surf => surf_lsm_h
          CALL calc_pt_q
       ENDIF

!
!--    Call for urban-type horizontal surfaces
       IF ( surf_usm_h%ns >= 1  )  THEN
          surf => surf_usm_h
          CALL calc_pt_q
       ENDIF

!
!--    Call for natural-type vertical surfaces
       DO  l = 0, 3
          IF ( surf_lsm_v(l)%ns >= 1  )  THEN
             surf => surf_lsm_v(l)
             CALL calc_pt_q
          ENDIF

!--       Call for urban-type vertical surfaces
          IF ( surf_usm_v(l)%ns >= 1  )  THEN
             surf => surf_usm_v(l)
             CALL calc_pt_q
          ENDIF
       ENDDO

!
!--    First, calculate the new Obukhov length, then new friction velocity,
!--    followed by the new scaling parameters (th*, q*, etc.), and the new
!--    surface fluxes if required. Note, each routine is called for different surface types. 
!--    First call for default-type horizontal surfaces, for natural- and 
!--    urban-type surfaces. Note, at this place only upward-facing horizontal
!--    surfaces are treated. 

!
!--    Default-type upward-facing horizontal surfaces
       IF ( surf_def_h(0)%ns >= 1  )  THEN
          surf => surf_def_h(0)
          CALL calc_uvw_abs
          IF ( .NOT. neutral )  CALL calc_ol
          CALL calc_us
          CALL calc_scaling_parameters
          CALL calc_surface_fluxes
       ENDIF
!
!--    Natural-type horizontal surfaces
       IF ( surf_lsm_h%ns >= 1  )  THEN
          surf => surf_lsm_h
          CALL calc_uvw_abs
          IF ( .NOT. neutral )  CALL calc_ol
          CALL calc_us
          CALL calc_scaling_parameters
          CALL calc_surface_fluxes
       ENDIF
!
!--    Urban-type horizontal surfaces
       IF ( surf_usm_h%ns >= 1  )  THEN
          surf => surf_usm_h
          CALL calc_uvw_abs
          IF ( .NOT. neutral )  CALL calc_ol
          CALL calc_us
          CALL calc_scaling_parameters
          CALL calc_surface_fluxes
!
!--       Calculate 10cm temperature, required in indoor model
          IF ( indoor_model )  CALL calc_pt_near_surface ( '10cm' )
       ENDIF

!
!--    Treat downward-facing horizontal surfaces. Note, so far, these are 
!--    always default type. Stratification is not considered
!--    in this case, hence, no further distinction between different 
!--    most_method is required.  
       IF ( surf_def_h(1)%ns >= 1  )  THEN
          downward = .TRUE.
          surf => surf_def_h(1)
          CALL calc_uvw_abs
          CALL calc_us
          CALL calc_surface_fluxes
          downward = .FALSE.
       ENDIF
!
!--    Calculate surfaces fluxes at vertical surfaces for momentum 
!--    and subgrid-scale TKE.
!--    No stability is considered. Therefore, scaling parameters and Obukhov-
!--    length do not need to be calculated and no distinction in 'circular',
!--    'Newton' or 'lookup' is necessary so far. 
!--    Note, this will change if stability is once considered.
       surf_vertical = .TRUE.
!
!--    Calculate horizontal momentum fluxes at north- and south-facing 
!--    surfaces(usvs).
!--    For default-type surfaces
       mom_uv = .TRUE. 
       DO  l = 0, 1
          IF ( surf_def_v(l)%ns >= 1  )  THEN
             surf => surf_def_v(l)
!
!--          Compute surface-parallel velocity
             CALL calc_uvw_abs_v_ugrid
!
!--          Compute respective friction velocity on staggered grid
             CALL calc_us
!
!--          Compute respective surface fluxes for momentum and TKE
             CALL calc_surface_fluxes
          ENDIF
       ENDDO
!
!--    For natural-type surfaces. Please note, even though stability is not
!--    considered for the calculation of momentum fluxes at vertical surfaces,
!--    scaling parameters and Obukov length are calculated nevertheless in this
!--    case. This is due to the requirement of ts in parameterization of heat
!--    flux in land-surface model in case of aero_resist_kray is not true.
       IF ( .NOT. aero_resist_kray )  THEN
          DO  l = 0, 1
             IF ( surf_lsm_v(l)%ns >= 1  )  THEN
                surf => surf_lsm_v(l)
!
!--             Compute surface-parallel velocity
                CALL calc_uvw_abs_v_ugrid
!
!--             Compute Obukhov length
                IF ( .NOT. neutral )  CALL calc_ol
!
!--             Compute respective friction velocity on staggered grid
                CALL calc_us
!
!--             Compute scaling parameters
                CALL calc_scaling_parameters
!
!--             Compute respective surface fluxes for momentum and TKE
                CALL calc_surface_fluxes
             ENDIF
          ENDDO
!
!--    No ts is required, so scaling parameters and Obukhov length do not need
!--    to be computed.
       ELSE
          DO  l = 0, 1
             IF ( surf_lsm_v(l)%ns >= 1  )  THEN
                surf => surf_lsm_v(l)
!   
!--             Compute surface-parallel velocity
                CALL calc_uvw_abs_v_ugrid
!
!--             Compute respective friction velocity on staggered grid
                CALL calc_us
!
!--             Compute respective surface fluxes for momentum and TKE
                CALL calc_surface_fluxes
             ENDIF
          ENDDO
       ENDIF
!
!--    For urban-type surfaces
       DO  l = 0, 1
          IF ( surf_usm_v(l)%ns >= 1  )  THEN
             surf => surf_usm_v(l)
!
!--          Compute surface-parallel velocity
             CALL calc_uvw_abs_v_ugrid
!
!--          Compute respective friction velocity on staggered grid
             CALL calc_us
!
!--          Compute respective surface fluxes for momentum and TKE
             CALL calc_surface_fluxes
!
!--          Calculate 10cm temperature, required in indoor model
             IF ( indoor_model )  CALL calc_pt_near_surface ( '10cm' )
          ENDIF
       ENDDO
!
!--    Calculate horizontal momentum fluxes at east- and west-facing 
!--    surfaces (vsus).
!--    For default-type surfaces
       DO  l = 2, 3
          IF ( surf_def_v(l)%ns >= 1  )  THEN
             surf => surf_def_v(l)
!
!--          Compute surface-parallel velocity
             CALL calc_uvw_abs_v_vgrid
!
!--          Compute respective friction velocity on staggered grid
             CALL calc_us
!
!--          Compute respective surface fluxes for momentum and TKE
             CALL calc_surface_fluxes
             
          ENDIF
       ENDDO
!
!--    For natural-type surfaces. Please note, even though stability is not
!--    considered for the calculation of momentum fluxes at vertical surfaces,
!--    scaling parameters and Obukov length are calculated nevertheless in this
!--    case. This is due to the requirement of ts in parameterization of heat
!--    flux in land-surface model in case of aero_resist_kray is not true.
       IF ( .NOT. aero_resist_kray )  THEN
          DO  l = 2, 3
             IF ( surf_lsm_v(l)%ns >= 1  )  THEN
                surf => surf_lsm_v(l)
!
!--             Compute surface-parallel velocity
                CALL calc_uvw_abs_v_vgrid
!
!--             Compute Obukhov length
                IF ( .NOT. neutral )  CALL calc_ol
!
!--             Compute respective friction velocity on staggered grid
                CALL calc_us
!
!--             Compute scaling parameters
                CALL calc_scaling_parameters
!
!--             Compute respective surface fluxes for momentum and TKE
                CALL calc_surface_fluxes
             ENDIF
          ENDDO
       ELSE
          DO  l = 2, 3
             IF ( surf_lsm_v(l)%ns >= 1  )  THEN
                surf => surf_lsm_v(l)
!   
!--             Compute surface-parallel velocity
                CALL calc_uvw_abs_v_vgrid
!
!--             Compute respective friction velocity on staggered grid
                CALL calc_us
!
!--             Compute respective surface fluxes for momentum and TKE
                CALL calc_surface_fluxes
             ENDIF
          ENDDO
       ENDIF
!
!--    For urban-type surfaces
       DO  l = 2, 3
          IF ( surf_usm_v(l)%ns >= 1  )  THEN
             surf => surf_usm_v(l)
!
!--          Compute surface-parallel velocity
             CALL calc_uvw_abs_v_vgrid
!
!--          Compute respective friction velocity on staggered grid
             CALL calc_us
!
!--          Compute respective surface fluxes for momentum and TKE
             CALL calc_surface_fluxes
!
!--          Calculate 10cm temperature, required in indoor model             
             IF ( indoor_model )  CALL calc_pt_near_surface ( '10cm' )
          ENDIF
       ENDDO
       mom_uv = .FALSE.
!
!--    Calculate horizontal momentum fluxes of w (wsus and wsvs) at vertial 
!--    surfaces.
       mom_w = .TRUE.
!
!--    Default-type surfaces
       DO  l = 0, 3
          IF ( surf_def_v(l)%ns >= 1  )  THEN
             surf => surf_def_v(l)
             CALL calc_uvw_abs_v_wgrid
             CALL calc_us
             CALL calc_surface_fluxes
          ENDIF
       ENDDO 
!
!--    Natural-type surfaces
       DO  l = 0, 3
          IF ( surf_lsm_v(l)%ns >= 1  )  THEN
             surf => surf_lsm_v(l)
             CALL calc_uvw_abs_v_wgrid
             CALL calc_us
             CALL calc_surface_fluxes
          ENDIF
       ENDDO 
!
!--    Urban-type surfaces
       DO  l = 0, 3
          IF ( surf_usm_v(l)%ns >= 1  )  THEN
             surf => surf_usm_v(l)
             CALL calc_uvw_abs_v_wgrid
             CALL calc_us
             CALL calc_surface_fluxes
          ENDIF
       ENDDO 
       mom_w = .FALSE.   
!
!--    Calculate momentum fluxes usvs, vsus, wsus and wsvs at vertical 
!--    surfaces for TKE production. Note, here, momentum fluxes are defined 
!--    at grid center and are not staggered as before.
       mom_tke = .TRUE.
!
!--    Default-type surfaces
       DO  l = 0, 3
          IF ( surf_def_v(l)%ns >= 1  )  THEN
             surf => surf_def_v(l)
             CALL calc_uvw_abs_v_sgrid
             CALL calc_us
             CALL calc_surface_fluxes
          ENDIF
       ENDDO 
!
!--    Natural-type surfaces
       DO  l = 0, 3
          IF ( surf_lsm_v(l)%ns >= 1  )  THEN
             surf => surf_lsm_v(l)
             CALL calc_uvw_abs_v_sgrid
             CALL calc_us
             CALL calc_surface_fluxes
          ENDIF
       ENDDO 
!
!--    Urban-type surfaces
       DO  l = 0, 3
          IF ( surf_usm_v(l)%ns >= 1  )  THEN
             surf => surf_usm_v(l)
             CALL calc_uvw_abs_v_sgrid
             CALL calc_us
             CALL calc_surface_fluxes
          ENDIF
       ENDDO 
       mom_tke = .FALSE.

       IF ( debug_output_timestep )  CALL debug_message( 'surface_layer_fluxes', 'end' )

    END SUBROUTINE surface_layer_fluxes


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initializing actions for the surface layer routine.
!------------------------------------------------------------------------------!
    SUBROUTINE init_surface_layer_fluxes

       IMPLICIT NONE

       INTEGER(iwp) ::  l  !< running index for vertical surface orientation

       CALL location_message( 'initializing surface layer', 'start' )

!
!--    In case of runs with neutral statification, set Obukhov length to a
!--    large value
       IF ( neutral )  THEN
          IF ( surf_def_h(0)%ns >= 1 )  surf_def_h(0)%ol = 1.0E10_wp
          IF ( surf_lsm_h%ns    >= 1 )  surf_lsm_h%ol    = 1.0E10_wp
          IF ( surf_usm_h%ns    >= 1 )  surf_usm_h%ol    = 1.0E10_wp
          
          DO  l = 0, 3
             IF ( surf_def_v(l)%ns >= 1  .AND.                                 &
                  ALLOCATED( surf_def_v(l)%ol ) )  surf_def_v(l)%ol = 1.0E10_wp
             IF ( surf_lsm_v(l)%ns >= 1  .AND.                                 &
                  ALLOCATED( surf_lsm_v(l)%ol ) )  surf_lsm_v(l)%ol = 1.0E10_wp
             IF ( surf_usm_v(l)%ns >= 1  .AND.                                 &
                  ALLOCATED( surf_usm_v(l)%ol ) )  surf_usm_v(l)%ol = 1.0E10_wp  
          ENDDO
          
       ENDIF

       CALL location_message( 'initializing surface layer', 'finished' )

    END SUBROUTINE init_surface_layer_fluxes


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the absolute value of the horizontal velocity (relative to the
!> surface) for horizontal surface elements. This is required by all methods.
!------------------------------------------------------------------------------!
    SUBROUTINE calc_uvw_abs
    
       IMPLICIT NONE

       INTEGER(iwp) ::  i             !< running index x direction
       INTEGER(iwp) ::  ibit          !< flag to mask computation of relative velocity in case of downward-facing surfaces
       INTEGER(iwp) ::  j             !< running index y direction
       INTEGER(iwp) ::  k             !< running index z direction
       INTEGER(iwp) ::  m             !< running index surface elements

       REAL(wp)     :: w_lfc          !< local free convection velocity scale
!
!--    ibit is 1 for upward-facing surfaces, zero for downward-facing surfaces.
       ibit = MERGE( 1, 0, .NOT. downward )

       !$OMP PARALLEL DO PRIVATE(i, j, k, w_lfc)
       !$ACC PARALLEL LOOP PRIVATE(i, j, k, w_lfc) &
       !$ACC PRESENT(surf, u, v)
       DO  m = 1, surf%ns

          i   = surf%i(m)            
          j   = surf%j(m)
          k   = surf%k(m)
       
!
!--       Calculate free convection velocity scale w_lfc is 
!--       use_free_convection_scaling = .T.. This will maintain a horizontal
!--       velocity even for very weak wind convective conditions. SIGN is used
!--       to set w_lfc to zero under stable conditions.
          IF ( use_free_convection_scaling )  THEN 
             w_lfc = ABS(g / surf%pt1(m) * surf%z_mo(m) * surf%shf(m))
             w_lfc = ( 0.5_wp * ( w_lfc + SIGN(w_lfc,surf%shf(m)) ) )**(0.33333_wp)
          ELSE
             w_lfc = 0.0_wp
          ENDIF

!
!--       Compute the absolute value of the horizontal velocity.
!--       (relative to the surface in case the lower surface is the ocean).
!--       Please note, in new surface modelling concept the index values changed,
!--       i.e. the reference grid point is not the surface-grid point itself but
!--       the first grid point outside of the topography. 
!--       Note, in case of coupled ocean-atmosphere simulations relative velocity
!--       with respect to the ocean surface is used, hence, (k-1,j,i) values
!--       are used to calculate the absolute velocity. 
!--       However, this do not apply for downward-facing walls. To mask this, 
!--       use ibit, which checks for upward/downward-facing surfaces. 
          surf%uvw_abs(m) = SQRT(                                              &
                              ( 0.5_wp * (   u(k,j,i)   + u(k,j,i+1)           &
                                        -  ( u(k-1,j,i) + u(k-1,j,i+1)         &
                                           ) * ibit                            &
                                         )                                     &
                              )**2 +                                           &
                              ( 0.5_wp * (   v(k,j,i)   + v(k,j+1,i)           &
                                        -  ( v(k-1,j,i) + v(k-1,j+1,i)         &
                                           ) * ibit                            &
                                         )                                     &
                              )**2  + w_lfc**2                                 &
                                )

                                

       ENDDO

    END SUBROUTINE calc_uvw_abs


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the absolute value of the horizontal velocity (relative to the
!> surface) for horizontal surface elements. This is required by all methods.
!------------------------------------------------------------------------------!
    SUBROUTINE calc_uvw_abs_v_ugrid

       IMPLICIT NONE

       INTEGER(iwp) ::  i   !< running index x direction
       INTEGER(iwp) ::  j   !< running index y direction
       INTEGER(iwp) ::  k   !< running index z direction
       INTEGER(iwp) ::  m   !< running index surface elements

       REAL(wp)     ::  u_i !< u-component on xu-grid
       REAL(wp)     ::  w_i !< w-component on xu-grid


       DO  m = 1, surf%ns
          i   = surf%i(m)            
          j   = surf%j(m)
          k   = surf%k(m)
!
!--       Compute the absolute value of the surface parallel velocity on u-grid.
          u_i  = u(k,j,i)
          w_i  = 0.25_wp * ( w(k-1,j,i-1) + w(k-1,j,i) +                       &
                             w(k,j,i-1)   + w(k,j,i) ) 

          surf%uvw_abs(m) = SQRT( u_i**2 + w_i**2 ) 

       ENDDO

    END SUBROUTINE calc_uvw_abs_v_ugrid

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the absolute value of the horizontal velocity (relative to the
!> surface) for horizontal surface elements. This is required by all methods.
!------------------------------------------------------------------------------!
    SUBROUTINE calc_uvw_abs_v_vgrid

       IMPLICIT NONE

       INTEGER(iwp) ::  i   !< running index x direction
       INTEGER(iwp) ::  j   !< running index y direction
       INTEGER(iwp) ::  k   !< running index z direction
       INTEGER(iwp) ::  m   !< running index surface elements

       REAL(wp)     ::  v_i !< v-component on yv-grid
       REAL(wp)     ::  w_i !< w-component on yv-grid


       DO  m = 1, surf%ns
          i   = surf%i(m)            
          j   = surf%j(m)
          k   = surf%k(m)

          v_i  = u(k,j,i)
          w_i  = 0.25_wp * ( w(k-1,j-1,i) + w(k-1,j,i) +                       &
                             w(k,j-1,i)   + w(k,j,i) ) 

          surf%uvw_abs(m) = SQRT( v_i**2 + w_i**2 ) 

       ENDDO

    END SUBROUTINE calc_uvw_abs_v_vgrid

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the absolute value of the horizontal velocity (relative to the
!> surface) for horizontal surface elements. This is required by all methods.
!------------------------------------------------------------------------------!
    SUBROUTINE calc_uvw_abs_v_wgrid

       IMPLICIT NONE

       INTEGER(iwp) ::  i   !< running index x direction
       INTEGER(iwp) ::  j   !< running index y direction
       INTEGER(iwp) ::  k   !< running index z direction
       INTEGER(iwp) ::  m   !< running index surface elements

       REAL(wp)     ::  u_i !< u-component on x-zw-grid
       REAL(wp)     ::  v_i !< v-component on y-zw-grid
       REAL(wp)     ::  w_i !< w-component on zw-grid
!
!--    North- (l=0) and south-facing (l=1) surfaces 
       IF ( l == 0  .OR.  l == 1 )  THEN
          DO  m = 1, surf%ns
             i   = surf%i(m)            
             j   = surf%j(m)
             k   = surf%k(m)

             u_i  = 0.25_wp * ( u(k+1,j,i+1) + u(k+1,j,i) +                    &
                                u(k,j,i+1)   + u(k,j,i) )
             v_i  = 0.0_wp
             w_i  = w(k,j,i)

             surf%uvw_abs(m) = SQRT( u_i**2 + v_i**2 + w_i**2 ) 
          ENDDO
!
!--    East- (l=2) and west-facing (l=3) surfaces
       ELSE
          DO  m = 1, surf%ns
             i   = surf%i(m)            
             j   = surf%j(m)
             k   = surf%k(m)

             u_i  = 0.0_wp
             v_i  = 0.25_wp * ( v(k+1,j+1,i) + v(k+1,j,i) +                    &
                                v(k,j+1,i)   + v(k,j,i) )
             w_i  = w(k,j,i)

             surf%uvw_abs(m) = SQRT( u_i**2 + v_i**2 + w_i**2 ) 
          ENDDO
       ENDIF            

    END SUBROUTINE calc_uvw_abs_v_wgrid

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the absolute value of the horizontal velocity (relative to the
!> surface) for horizontal surface elements. This is required by all methods.
!------------------------------------------------------------------------------!
    SUBROUTINE calc_uvw_abs_v_sgrid

       IMPLICIT NONE

       INTEGER(iwp) ::  i   !< running index x direction
       INTEGER(iwp) ::  j   !< running index y direction
       INTEGER(iwp) ::  k   !< running index z direction
       INTEGER(iwp) ::  m   !< running index surface elements

       REAL(wp)     ::  u_i !< u-component on scalar grid
       REAL(wp)     ::  v_i !< v-component on scalar grid
       REAL(wp)     ::  w_i !< w-component on scalar grid

!
!--    North- (l=0) and south-facing (l=1) walls 
       IF ( l == 0  .OR.  l == 1 )  THEN
          DO  m = 1, surf%ns
             i   = surf%i(m)            
             j   = surf%j(m)
             k   = surf%k(m)

             u_i = 0.5_wp * ( u(k,j,i) + u(k,j,i+1) )
             v_i = 0.0_wp
             w_i = 0.5_wp * ( w(k,j,i) + w(k-1,j,i) )

             surf%uvw_abs(m) = SQRT( u_i**2 + v_i**2 + w_i**2 ) 
          ENDDO
!
!--    East- (l=2) and west-facing (l=3) walls 
       ELSE
          DO  m = 1, surf%ns
             i   = surf%i(m)            
             j   = surf%j(m)
             k   = surf%k(m)

             u_i = 0.0_wp
             v_i = 0.5_wp * ( v(k,j,i) + v(k,j+1,i) )
             w_i = 0.5_wp * ( w(k,j,i) + w(k-1,j,i) )

             surf%uvw_abs(m) = SQRT( u_i**2 + v_i**2 + w_i**2 ) 
          ENDDO
       ENDIF  

    END SUBROUTINE calc_uvw_abs_v_sgrid


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate the Obukhov length (L) and Richardson flux number (z/L)
!------------------------------------------------------------------------------!
    SUBROUTINE calc_ol

       IMPLICIT NONE

       INTEGER(iwp) ::  iter    !< Newton iteration step
       INTEGER(iwp) ::  m       !< loop variable over all horizontal wall elements 

       LOGICAL, DIMENSION(surf%ns) ::  convergence_reached  !< convergence switch for vectorization
       !$ACC DECLARE CREATE( convergence_reached )

       REAL(wp)     :: f,      & !< Function for Newton iteration: f = Ri - [...]/[...]^2 = 0
                       f_d_ol, & !< Derivative of f
                       ol_l,   & !< Lower bound of L for Newton iteration
                       ol_m,   & !< Previous value of L for Newton iteration
                       ol_old, & !< Previous time step value of L
                       ol_u      !< Upper bound of L for Newton iteration

       REAL(wp), DIMENSION(surf%ns) ::  ol_old_vec  !< temporary array required for vectorization
       REAL(wp), DIMENSION(surf%ns) ::  z_mo_vec    !< temporary array required for vectorization
       !$ACC DECLARE CREATE( ol_old_vec, z_mo_vec )

!
!--    Evaluate bulk Richardson number (calculation depends on
!--    definition based on setting of boundary conditions
       IF ( ibc_pt_b /= 1 )  THEN
          IF ( humidity )  THEN
             !$OMP PARALLEL DO PRIVATE( z_mo )
             DO  m = 1, surf%ns

                z_mo = surf%z_mo(m)

                surf%rib(m) = g * z_mo                                         &
                              * ( surf%vpt1(m) - surf%vpt_surface(m) )         &
                              / ( surf%uvw_abs(m)**2 * surf%vpt1(m)            &
                              + 1.0E-20_wp )
             ENDDO
          ELSE
             !$OMP PARALLEL DO PRIVATE( z_mo )
             DO  m = 1, surf%ns

                z_mo = surf%z_mo(m)

                surf%rib(m) = g * z_mo                                         &
                              * ( surf%pt1(m) - surf%pt_surface(m) )           &
                              / ( surf%uvw_abs(m)**2 * surf%pt1(m) + 1.0E-20_wp )
             ENDDO
          ENDIF
       ELSE
          IF ( humidity )  THEN
             !$OMP PARALLEL DO PRIVATE( k, z_mo )
             DO  m = 1, surf%ns

                k   = surf%k(m)

                z_mo = surf%z_mo(m)

                surf%rib(m) = - g * z_mo * ( ( 1.0_wp + 0.61_wp                &
                           * surf%qv1(m) ) * surf%shf(m) + 0.61_wp             &
                           * surf%pt1(m) * surf%qsws(m) ) *                    &
                             drho_air_zw(k-1)                /                 &
                         ( surf%uvw_abs(m)**3 * surf%vpt1(m) * kappa**2        &
                           + 1.0E-20_wp )
             ENDDO
          ELSE
             !$OMP PARALLEL DO PRIVATE( k, z_mo )
             !$ACC PARALLEL LOOP PRIVATE(k, z_mo) &
             !$ACC PRESENT(surf, drho_air_zw)
             DO  m = 1, surf%ns

                k   = surf%k(m)

                z_mo = surf%z_mo(m)

                surf%rib(m) = - g * z_mo * surf%shf(m) *                    &
                                     drho_air_zw(k-1)            /          &
                        ( surf%uvw_abs(m)**3 * surf%pt1(m) * kappa**2       &
                        + 1.0E-20_wp )
             ENDDO
          ENDIF
       ENDIF

       IF ( loop_optimization == 'cache' )  THEN
!
!--       Calculate the Obukhov length using Newton iteration
          !$OMP PARALLEL DO PRIVATE(i, j, z_mo) &
          !$OMP PRIVATE(ol_old, ol_m, ol_l, ol_u, f, f_d_ol)
          !$ACC PARALLEL LOOP PRIVATE(i, j, z_mo) &
          !$ACC PRIVATE(ol_old, ol_m, ol_l, ol_u, f, f_d_ol) &
          !$ACC PRESENT(surf)
          DO  m = 1, surf%ns

             i   = surf%i(m)
             j   = surf%j(m)

             z_mo = surf%z_mo(m)

!
!--          Store current value in case the Newton iteration fails
             ol_old = surf%ol(m)

!
!--          Ensure that the bulk Richardson number and the Obukhov
!--          length have the same sign
             IF ( surf%rib(m) * surf%ol(m) < 0.0_wp  .OR.  ABS( surf%ol(m) ) == ol_max )  THEN
                IF ( surf%rib(m) > 1.0_wp ) surf%ol(m) =  0.01_wp
                IF ( surf%rib(m) < 0.0_wp ) surf%ol(m) = -0.01_wp
             ENDIF
!
!--          Iteration to find Obukhov length
             iter = 0
             DO
                iter = iter + 1
!
!--             In case of divergence, use the value of the previous time step
                IF ( iter > 1000 )  THEN
                   surf%ol(m) = ol_old
                   EXIT
                ENDIF

                ol_m = surf%ol(m)
                ol_l = ol_m - 0.001_wp * ol_m
                ol_u = ol_m + 0.001_wp * ol_m


                IF ( ibc_pt_b /= 1 )  THEN
!
!--                Calculate f = Ri - [...]/[...]^2 = 0
                   f = surf%rib(m) - ( z_mo / ol_m ) * ( LOG( z_mo / surf%z0h(m) )                 &
                                                       - psi_h( z_mo / ol_m )                      &
                                                       + psi_h( surf%z0h(m) / ol_m )               &
                                                       ) /                                         &
                                                       ( LOG( z_mo / surf%z0(m) )                  &
                                                      - psi_m( z_mo / ol_m )                       &
                                                      + psi_m( surf%z0(m) /  ol_m )                &
                                                       )**2

!
!--                Calculate df/dL
                   f_d_ol = ( - ( z_mo / ol_u ) * ( LOG( z_mo / surf%z0h(m) )                      &
                                                  - psi_h( z_mo / ol_u )                           &
                                                  + psi_h( surf%z0h(m) / ol_u )                    &
                                                  ) /                                              &
                                                  ( LOG( z_mo / surf%z0(m) )                       &
                                                  - psi_m( z_mo / ol_u )                           &
                                                  + psi_m( surf%z0(m) / ol_u )                     &
                                                  )**2                                             &
                              + ( z_mo / ol_l ) * ( LOG( z_mo / surf%z0h(m) )                      &
                                                  - psi_h( z_mo / ol_l )                           &
                                                  + psi_h( surf%z0h(m) / ol_l )                    &
                                                  ) /                                              &
                                                  ( LOG( z_mo / surf%z0(m) )                       &
                                                  - psi_m( z_mo / ol_l )                           &
                                                  + psi_m( surf%z0(m) / ol_l )                     &
                                                  )**2                                             &
                           ) / ( ol_u - ol_l )
                ELSE
!
!--                Calculate f = Ri - 1 /[...]^3 = 0
                   f = surf%rib(m) - ( z_mo / ol_m ) / ( LOG( z_mo / surf%z0(m) )                  &
                                                       - psi_m( z_mo / ol_m )                      &
                                                       + psi_m( surf%z0(m) / ol_m )                &
                                                       )**3

!
!--                Calculate df/dL
                   f_d_ol = ( - ( z_mo / ol_u ) / ( LOG( z_mo / surf%z0(m) )                       &
                                                  - psi_m( z_mo / ol_u )                           &
                                                  + psi_m( surf%z0(m) / ol_u )                     &
                                                  )**3                                             &
                              + ( z_mo / ol_l ) / ( LOG( z_mo / surf%z0(m) )                       &
                                                  - psi_m( z_mo / ol_l )                           &
                                                  + psi_m( surf%z0(m) / ol_l )                     &
                                                  )**3                                             &
                             ) / ( ol_u - ol_l )
                ENDIF
!
!--             Calculate new L
                surf%ol(m) = ol_m - f / f_d_ol

!
!--             Ensure that the bulk Richardson number and the Obukhov
!--             length have the same sign and ensure convergence.
                IF ( surf%ol(m) * ol_m < 0.0_wp )  surf%ol(m) = ol_m * 0.5_wp
!
!--             If unrealistic value occurs, set L to the maximum
!--             value that is allowed
                IF ( ABS( surf%ol(m) ) > ol_max )  THEN
                   surf%ol(m) = ol_max
                   EXIT
                ENDIF
!
!--             Assure that Obukhov length does not become zero. If the limit is
!--             reached, exit the loop.
                IF ( ABS( surf%ol(m) ) < 1E-5_wp )  THEN
                   surf%ol(m) = SIGN( 1E-5_wp, surf%ol(m) )
                   EXIT
                ENDIF
!
!--             Check for convergence
                IF ( ABS( ( surf%ol(m) - ol_m ) /  surf%ol(m) ) < 1.0E-4_wp )  EXIT

             ENDDO
          ENDDO

!
!--    Vector Version
       ELSE
!
!--       Calculate the Obukhov length using Newton iteration
!--       First set arrays required for vectorization
          !$ACC PARALLEL LOOP &
          !$ACC PRESENT(surf)
          DO  m = 1, surf%ns

             z_mo_vec(m) = surf%z_mo(m)

!
!--          Store current value in case the Newton iteration fails
             ol_old_vec(m) = surf%ol(m)

!
!--          Ensure that the bulk Richardson number and the Obukhov length have the same sign
             IF ( surf%rib(m) * surf%ol(m) < 0.0_wp  .OR.  ABS( surf%ol(m) ) == ol_max )  THEN
                IF ( surf%rib(m) > 1.0_wp ) surf%ol(m) =  0.01_wp
                IF ( surf%rib(m) < 0.0_wp ) surf%ol(m) = -0.01_wp
             ENDIF
!
!--          Initialize convergence flag
             convergence_reached(m) = .FALSE.

          ENDDO

!
!--       Iteration to find Obukhov length
          iter = 0
          DO

             iter = iter + 1
!
!--          In case of divergence, use the value(s) of the previous time step
             IF ( iter > 1000 )  THEN
                !$ACC PARALLEL LOOP &
                !$ACC PRESENT(surf)
                DO  m = 1, surf%ns
                   IF ( .NOT. convergence_reached(m) )  surf%ol(m) = ol_old_vec(m)
                ENDDO
                EXIT
             ENDIF

             !$ACC PARALLEL LOOP PRIVATE(ol_m, ol_l, ol_u, f, f_d_ol) &
             !$ACC PRESENT(surf)
             DO  m = 1, surf%ns

                IF ( convergence_reached(m) )  CYCLE

                ol_m = surf%ol(m)
                ol_l = ol_m - 0.001_wp * ol_m
                ol_u = ol_m + 0.001_wp * ol_m


                IF ( ibc_pt_b /= 1 )  THEN
!
!--                Calculate f = Ri - [...]/[...]^2 = 0
                   f = surf%rib(m) - ( z_mo_vec(m) / ol_m ) * ( LOG( z_mo_vec(m) / surf%z0h(m) )   &
                                                              - psi_h( z_mo_vec(m) / ol_m )        &
                                                              + psi_h( surf%z0h(m) / ol_m )        &
                                                              ) /                                  &
                                                              ( LOG( z_mo_vec(m) / surf%z0(m) )    &
                                                             - psi_m( z_mo_vec(m) / ol_m )         &
                                                             + psi_m( surf%z0(m) /  ol_m )         &
                                                              )**2

!
!--                Calculate df/dL
                   f_d_ol = ( - ( z_mo_vec(m) / ol_u ) * ( LOG( z_mo_vec(m) / surf%z0h(m) )        &
                                                         - psi_h( z_mo_vec(m) / ol_u )             &
                                                         + psi_h( surf%z0h(m) / ol_u )             &
                                                         ) /                                       &
                                                         ( LOG( z_mo_vec(m) / surf%z0(m) )         &
                                                         - psi_m( z_mo_vec(m) / ol_u )             &
                                                         + psi_m( surf%z0(m) / ol_u )              &
                                                         )**2                                      &
                              + ( z_mo_vec(m) / ol_l ) * ( LOG( z_mo_vec(m) / surf%z0h(m) )        &
                                                         - psi_h( z_mo_vec(m) / ol_l )             &
                                                         + psi_h( surf%z0h(m) / ol_l )             &
                                                         ) /                                       &
                                                         ( LOG( z_mo_vec(m) / surf%z0(m) )         &
                                                         - psi_m( z_mo_vec(m) / ol_l )             &
                                                         + psi_m( surf%z0(m) / ol_l )              &
                                                         )**2                                      &
                            ) / ( ol_u - ol_l )
                ELSE
!
!--                Calculate f = Ri - 1 /[...]^3 = 0
                   f = surf%rib(m) - ( z_mo_vec(m) / ol_m ) / ( LOG( z_mo_vec(m) / surf%z0(m) )    &
                                                              - psi_m( z_mo_vec(m) / ol_m )        &
                                                              + psi_m( surf%z0(m) / ol_m )         &
                                                              )**3

!
!--                Calculate df/dL
                   f_d_ol = ( - ( z_mo_vec(m) / ol_u ) / ( LOG( z_mo_vec(m) / surf%z0(m) )         &
                                                         - psi_m( z_mo_vec(m) / ol_u )             &
                                                         + psi_m( surf%z0(m) / ol_u )              &
                                                         )**3                                      &
                              + ( z_mo_vec(m) / ol_l ) / ( LOG( z_mo_vec(m) / surf%z0(m) )         &
                                                         - psi_m( z_mo_vec(m) / ol_l )             &
                                                         + psi_m( surf%z0(m) / ol_l )              &
                                                         )**3                                      &
                            ) / ( ol_u - ol_l )
                ENDIF
!
!--             Calculate new L
                surf%ol(m) = ol_m - f / f_d_ol

!
!--             Ensure that the bulk Richardson number and the Obukhov
!--             length have the same sign and ensure convergence.
                IF ( surf%ol(m) * ol_m < 0.0_wp )  surf%ol(m) = ol_m * 0.5_wp

!
!--             Check for convergence
!--             This check does not modify surf%ol, therefore this is done first
                IF ( ABS( ( surf%ol(m) - ol_m ) /  surf%ol(m) ) < 1.0E-4_wp )  THEN
                   convergence_reached(m) = .TRUE.
                ENDIF
!
!--             If unrealistic value occurs, set L to the maximum allowed value
                IF ( ABS( surf%ol(m) ) > ol_max )  THEN
                   surf%ol(m) = ol_max
                   convergence_reached(m) = .TRUE.
                ENDIF

             ENDDO
!
!--          Assure that Obukhov length does not become zero
             !$ACC PARALLEL LOOP &
             !$ACC PRESENT(surf)
             DO  m = 1, surf%ns
                IF ( convergence_reached(m) )  CYCLE
                IF ( ABS( surf%ol(m) ) < 1E-5_wp )  THEN
                   surf%ol(m) = SIGN( 10E-6_wp, surf%ol(m) )
                   convergence_reached(m) = .TRUE.
                ENDIF
             ENDDO

             IF ( ALL( convergence_reached ) )  EXIT

          ENDDO  ! end of iteration loop

       ENDIF  ! end of vector branch

    END SUBROUTINE calc_ol

!
!-- Calculate friction velocity u*
    SUBROUTINE calc_us

       IMPLICIT NONE

       INTEGER(iwp) ::  m       !< loop variable over all horizontal surf elements 

!
!--    Compute u* at horizontal surfaces at the scalars' grid points
       IF ( .NOT. surf_vertical )  THEN
!
!--       Compute u* at upward-facing surfaces
          IF ( .NOT. downward )  THEN
             !$OMP PARALLEL  DO PRIVATE( z_mo )
             !$ACC PARALLEL LOOP PRIVATE(z_mo) &
             !$ACC PRESENT(surf)
             DO  m = 1, surf%ns

                z_mo = surf%z_mo(m)
!
!--             Compute u* at the scalars' grid points
                surf%us(m) = kappa * surf%uvw_abs(m) /                         &
                            ( LOG( z_mo / surf%z0(m) )                         &
                           - psi_m( z_mo / surf%ol(m) )                        &
                           + psi_m( surf%z0(m) / surf%ol(m) ) )
   
             ENDDO
!
!--       Compute u* at downward-facing surfaces. This case, do not consider
!--       any stability. 
          ELSE
             !$OMP PARALLEL  DO PRIVATE( z_mo )
             !$ACC PARALLEL LOOP PRIVATE(z_mo) &
             !$ACC PRESENT(surf)
             DO  m = 1, surf%ns

                z_mo = surf%z_mo(m)
!
!--             Compute u* at the scalars' grid points
                surf%us(m) = kappa * surf%uvw_abs(m) / LOG( z_mo / surf%z0(m) )
   
             ENDDO
          ENDIF
!
!--    Compute u* at vertical surfaces at the u/v/v grid, respectively. 
!--    No stability is considered in this case.
       ELSE
          !$OMP PARALLEL DO PRIVATE( z_mo )
          !$ACC PARALLEL LOOP PRIVATE(z_mo) &
          !$ACC PRESENT(surf)
          DO  m = 1, surf%ns
             z_mo = surf%z_mo(m)

             surf%us(m) = kappa * surf%uvw_abs(m) / LOG( z_mo / surf%z0(m) )
          ENDDO
       ENDIF

    END SUBROUTINE calc_us

!
!-- Calculate potential temperature, specific humidity, and virtual potential
!-- temperature at first grid level
    SUBROUTINE calc_pt_q

       IMPLICIT NONE

       INTEGER(iwp) ::  m       !< loop variable over all horizontal surf elements 

       !$OMP PARALLEL DO PRIVATE( i, j, k )
       !$ACC PARALLEL LOOP PRIVATE(i, j, k) &
       !$ACC PRESENT(surf, pt)
       DO  m = 1, surf%ns 

          i   = surf%i(m)            
          j   = surf%j(m)
          k   = surf%k(m)

#ifndef _OPENACC
          IF ( bulk_cloud_model ) THEN
             surf%pt1(m) = pt(k,j,i) + lv_d_cp * d_exner(k) * ql(k,j,i)
             surf%qv1(m) = q(k,j,i) - ql(k,j,i)
          ELSEIF( cloud_droplets ) THEN
             surf%pt1(m) = pt(k,j,i) + lv_d_cp * d_exner(k) * ql(k,j,i)
             surf%qv1(m) = q(k,j,i) 
          ELSE
#endif
             surf%pt1(m) = pt(k,j,i)
#ifndef _OPENACC
             IF ( humidity )  THEN
                surf%qv1(m) = q(k,j,i)
             ELSE
#endif
                surf%qv1(m) = 0.0_wp
#ifndef _OPENACC
             ENDIF
          ENDIF 

          IF ( humidity )  THEN
             surf%vpt1(m) = pt(k,j,i) * ( 1.0_wp + 0.61_wp * q(k,j,i) )
          ENDIF
#endif
          
       ENDDO

    END SUBROUTINE calc_pt_q


!
!-- Set potential temperature at surface grid level.
!-- ( only for upward-facing surfs )
    SUBROUTINE calc_pt_surface

       IMPLICIT NONE

       INTEGER(iwp) ::  k_off    !< index offset between surface and atmosphere grid point (-1 for upward-, +1 for downward-facing walls) 
       INTEGER(iwp) ::  m       !< loop variable over all horizontal surf elements 
       
       k_off = surf%koff
       !$OMP PARALLEL DO PRIVATE( i, j, k )
       !$ACC PARALLEL LOOP PRIVATE(i, j, k) &
       !$ACC PRESENT(surf, pt)
       DO  m = 1, surf%ns 

          i   = surf%i(m)            
          j   = surf%j(m)
          k   = surf%k(m)

          surf%pt_surface(m) = pt(k+k_off,j,i)

       ENDDO

    END SUBROUTINE calc_pt_surface

!
!-- Set mixing ratio at surface grid level. ( Only for upward-facing surfs. )
    SUBROUTINE calc_q_surface

       IMPLICIT NONE

       INTEGER(iwp) ::  k_off   !< index offset between surface and atmosphere grid point (-1 for upward-, +1 for downward-facing walls) 
       INTEGER(iwp) ::  m       !< loop variable over all horizontal surf elements 
       
       k_off = surf%koff
       !$OMP PARALLEL DO PRIVATE( i, j, k )
       DO  m = 1, surf%ns 

          i   = surf%i(m)            
          j   = surf%j(m)
          k   = surf%k(m)

          surf%q_surface(m) = q(k+k_off,j,i)

       ENDDO

    END SUBROUTINE calc_q_surface
    
!
!-- Set virtual potential temperature at surface grid level.
!-- ( only for upward-facing surfs )
    SUBROUTINE calc_vpt_surface

       IMPLICIT NONE

       INTEGER(iwp) ::  k_off    !< index offset between surface and atmosphere grid point (-1 for upward-, +1 for downward-facing walls) 
       INTEGER(iwp) ::  m       !< loop variable over all horizontal surf elements 
       
       k_off = surf%koff
       !$OMP PARALLEL DO PRIVATE( i, j, k )
       DO  m = 1, surf%ns 

          i   = surf%i(m)            
          j   = surf%j(m)
          k   = surf%k(m)

          surf%vpt_surface(m) = vpt(k+k_off,j,i)

       ENDDO

    END SUBROUTINE calc_vpt_surface
    
!
!-- Calculate the other MOST scaling parameters theta*, q*, (qc*, qr*, nc*, nr*)
    SUBROUTINE calc_scaling_parameters

       IMPLICIT NONE


       INTEGER(iwp)  ::  m       !< loop variable over all horizontal surf elements 
       INTEGER(iwp)  ::  lsp     !< running index for chemical species
! 
!--    Compute theta* at horizontal surfaces
       IF ( constant_heatflux  .AND.  .NOT. surf_vertical )  THEN
!
!--       For a given heat flux in the surface layer:

          !$OMP PARALLEL DO PRIVATE( i, j, k )
          !$ACC PARALLEL LOOP PRIVATE(i, j, k) &
          !$ACC PRESENT(surf, drho_air_zw)
          DO  m = 1, surf%ns 

             i   = surf%i(m)            
             j   = surf%j(m)
             k   = surf%k(m)

             surf%ts(m) = -surf%shf(m) * drho_air_zw(k-1) /                    &
                                  ( surf%us(m) + 1E-30_wp )

!
!--          ts must be limited, because otherwise overflow may occur in case
!--          of us=0 when computing ol further below
             IF ( surf%ts(m) < -1.05E5_wp )  surf%ts(m) = -1.0E5_wp
             IF ( surf%ts(m) >  1.0E5_wp  )  surf%ts(m) =  1.0E5_wp

          ENDDO

       ELSEIF ( .NOT. surf_vertical ) THEN
!
!--       For a given surface temperature:
          IF ( large_scale_forcing  .AND.  lsf_surf )  THEN

             !$OMP PARALLEL DO PRIVATE( i, j, k )
             DO  m = 1, surf%ns 
                i   = surf%i(m)            
                j   = surf%j(m)
                k   = surf%k(m)

                pt(k-1,j,i) = pt_surface
             ENDDO
          ENDIF

          !$OMP PARALLEL DO PRIVATE( z_mo )
          DO  m = 1, surf%ns   

             z_mo = surf%z_mo(m)

             surf%ts(m) = kappa * ( surf%pt1(m) - surf%pt_surface(m) )      &
                                  / ( LOG( z_mo / surf%z0h(m) )             &
                                      - psi_h( z_mo / surf%ol(m) )          &
                                      + psi_h( surf%z0h(m) / surf%ol(m) ) )

          ENDDO

       ENDIF
! 
!--    Compute theta* at vertical surfaces. This is only required in case of 
!--    land-surface model, in order to compute aerodynamical resistance.
       IF ( surf_vertical )  THEN
          !$OMP PARALLEL DO PRIVATE( i, j )
          DO  m = 1, surf%ns 

             i          =  surf%i(m)            
             j          =  surf%j(m)
             surf%ts(m) = -surf%shf(m) / ( surf%us(m) + 1E-30_wp )
!
!--          ts must be limited, because otherwise overflow may occur in case
!--          of us=0 when computing ol further below
             IF ( surf%ts(m) < -1.05E5_wp )  surf%ts(m) = -1.0E5_wp
             IF ( surf%ts(m) >  1.0E5_wp  )  surf%ts(m) =  1.0E5_wp

          ENDDO
       ENDIF

!
!--    If required compute q* at horizontal surfaces
       IF ( humidity )  THEN
          IF ( constant_waterflux  .AND.  .NOT. surf_vertical )  THEN
!
!--          For a given water flux in the surface layer
             !$OMP PARALLEL DO PRIVATE( i, j, k )
             DO  m = 1, surf%ns  

                i   = surf%i(m)            
                j   = surf%j(m)
                k   = surf%k(m)
                surf%qs(m) = -surf%qsws(m) * drho_air_zw(k-1) /                &
                                               ( surf%us(m) + 1E-30_wp )

             ENDDO

          ELSEIF ( .NOT. surf_vertical )  THEN 
             coupled_run = ( coupling_mode == 'atmosphere_to_ocean'  .AND.     &
                             run_coupled )

             IF ( large_scale_forcing  .AND.  lsf_surf )  THEN
                !$OMP PARALLEL DO PRIVATE( i, j, k )
                DO  m = 1, surf%ns  

                   i   = surf%i(m)           
                   j   = surf%j(m)
                   k   = surf%k(m)
                   q(k-1,j,i) = q_surface
                   
                ENDDO
             ENDIF

!
!--          Assume saturation for atmosphere coupled to ocean (but not
!--          in case of precursor runs)
             IF ( coupled_run )  THEN
                !$OMP PARALLEL DO PRIVATE( i, j, k, e_s )
                DO  m = 1, surf%ns   
                   i   = surf%i(m)            
                   j   = surf%j(m)
                   k   = surf%k(m)
                   e_s = 6.1_wp * &
                              EXP( 0.07_wp * ( MIN(pt(k-1,j,i),pt(k,j,i))      &
                                               - 273.15_wp ) )
                   q(k-1,j,i) = rd_d_rv * e_s / ( surface_pressure - e_s )
                ENDDO
             ENDIF

             IF ( bulk_cloud_model  .OR.  cloud_droplets )  THEN
               !$OMP PARALLEL DO PRIVATE( i, j, k, z_mo )
                DO  m = 1, surf%ns   

                   i   = surf%i(m)            
                   j   = surf%j(m)
                   k   = surf%k(m)
   
                   z_mo = surf%z_mo(m)

                   surf%qs(m) = kappa * ( surf%qv1(m) - surf%q_surface(m) )    &
                                        / ( LOG( z_mo / surf%z0q(m) )          &
                                            - psi_h( z_mo / surf%ol(m) )       &
                                            + psi_h( surf%z0q(m) /             &
                                                     surf%ol(m) ) )
                ENDDO
             ELSE
                !$OMP PARALLEL DO PRIVATE( i, j, k, z_mo )
                DO  m = 1, surf%ns   

                   i   = surf%i(m)            
                   j   = surf%j(m)
                   k   = surf%k(m)
   
                   z_mo = surf%z_mo(m)

                   surf%qs(m) = kappa * ( q(k,j,i) - q(k-1,j,i) )              &
                                        / ( LOG( z_mo / surf%z0q(m) )          &
                                            - psi_h( z_mo / surf%ol(m) )       &
                                            + psi_h( surf%z0q(m) /             &
                                                     surf%ol(m) ) )
                ENDDO
             ENDIF
          ENDIF
! 
!--       Compute q* at vertical surfaces
          IF ( surf_vertical )  THEN
             !$OMP PARALLEL DO PRIVATE( i, j )
             DO  m = 1, surf%ns  

                i          =  surf%i(m)            
                j          =  surf%j(m)
                surf%qs(m) = -surf%qsws(m) / ( surf%us(m) + 1E-30_wp )

             ENDDO
          ENDIF
       ENDIF
       
!
!--    If required compute s*
       IF ( passive_scalar )  THEN
!
!--       At horizontal surfaces
          IF ( constant_scalarflux  .AND.  .NOT. surf_vertical )  THEN
!
!--          For a given scalar flux in the surface layer
             !$OMP PARALLEL DO PRIVATE( i, j )
             DO  m = 1, surf%ns  
                i   = surf%i(m)            
                j   = surf%j(m)
                surf%ss(m) = -surf%ssws(m) / ( surf%us(m) + 1E-30_wp )
             ENDDO
          ENDIF
!
!--       At vertical surfaces
          IF ( surf_vertical )  THEN
             !$OMP PARALLEL DO PRIVATE( i, j )
             DO  m = 1, surf%ns  
                i          =  surf%i(m)            
                j          =  surf%j(m)
                surf%ss(m) = -surf%ssws(m) / ( surf%us(m) + 1E-30_wp )
             ENDDO
          ENDIF
       ENDIF

!
!--    If required compute cs* (chemical species)
       IF ( air_chemistry  )  THEN  
!
!--       At horizontal surfaces                             
          DO  lsp = 1, nvar
             IF ( constant_csflux(lsp)  .AND.  .NOT.  surf_vertical )  THEN
!--             For a given chemical species' flux in the surface layer
                !$OMP PARALLEL DO PRIVATE( i, j )
                DO  m = 1, surf%ns  
                   i   = surf%i(m)            
                   j   = surf%j(m)
                   surf%css(lsp,m) = -surf%cssws(lsp,m) / ( surf%us(m) + 1E-30_wp )
                ENDDO
             ENDIF
          ENDDO
!
!--       At vertical surfaces
          IF ( surf_vertical )  THEN
             DO  lsp = 1, nvar
                !$OMP PARALLEL DO PRIVATE( i, j )
                DO  m = 1, surf%ns  
                   i   = surf%i(m)            
                   j   = surf%j(m)
                   surf%css(lsp,m) = -surf%cssws(lsp,m) / ( surf%us(m) + 1E-30_wp )
                ENDDO
             ENDDO
          ENDIF
       ENDIF

!
!--    If required compute qc* and nc*
       IF ( bulk_cloud_model  .AND.  microphysics_morrison  .AND.              &
            .NOT. surf_vertical )  THEN
          !$OMP PARALLEL DO PRIVATE( i, j, k, z_mo )
          DO  m = 1, surf%ns    

             i    = surf%i(m)            
             j    = surf%j(m)
             k    = surf%k(m)

             z_mo = surf%z_mo(m)

             surf%qcs(m) = kappa * ( qc(k,j,i) - qc(k-1,j,i) )                 &
                                 / ( LOG( z_mo / surf%z0q(m) )                 &
                                 - psi_h( z_mo / surf%ol(m) )                  &
                                 + psi_h( surf%z0q(m) / surf%ol(m) ) )

             surf%ncs(m) = kappa * ( nc(k,j,i) - nc(k-1,j,i) )                 &
                                 / ( LOG( z_mo / surf%z0q(m) )                 &
                                 - psi_h( z_mo / surf%ol(m) )                  &
                                 + psi_h( surf%z0q(m) / surf%ol(m) ) )
          ENDDO

       ENDIF

!
!--    If required compute qr* and nr*
       IF ( bulk_cloud_model  .AND.  microphysics_seifert  .AND.               &
            .NOT. surf_vertical )  THEN
          !$OMP PARALLEL DO PRIVATE( i, j, k, z_mo )
          DO  m = 1, surf%ns    

             i    = surf%i(m)            
             j    = surf%j(m)
             k    = surf%k(m)

             z_mo = surf%z_mo(m)

             surf%qrs(m) = kappa * ( qr(k,j,i) - qr(k-1,j,i) )                 &
                                 / ( LOG( z_mo / surf%z0q(m) )                 &
                                 - psi_h( z_mo / surf%ol(m) )                  &
                                 + psi_h( surf%z0q(m) / surf%ol(m) ) )

             surf%nrs(m) = kappa * ( nr(k,j,i) - nr(k-1,j,i) )                 &
                                 / ( LOG( z_mo / surf%z0q(m) )                 &
                                 - psi_h( z_mo / surf%ol(m) )                  &
                                 + psi_h( surf%z0q(m) / surf%ol(m) ) )
          ENDDO

       ENDIF

    END SUBROUTINE calc_scaling_parameters



!
!-- Calculate surface fluxes usws, vsws, shf, qsws, (qcsws, qrsws, ncsws, nrsws)
    SUBROUTINE calc_surface_fluxes

       IMPLICIT NONE

       INTEGER(iwp)  ::  m       !< loop variable over all horizontal surf elements
       INTEGER(iwp)  ::  lsp     !< running index for chemical species

       REAL(wp)                            ::  dum     !< dummy to precalculate logarithm
       REAL(wp)                            ::  flag_u  !< flag indicating u-grid, used for calculation of horizontal momentum fluxes at vertical surfaces
       REAL(wp)                            ::  flag_v  !< flag indicating v-grid, used for calculation of horizontal momentum fluxes at vertical surfaces
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  u_i     !< u-component interpolated onto scalar grid point, required for momentum fluxes at vertical surfaces 
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  v_i     !< v-component interpolated onto scalar grid point, required for momentum fluxes at vertical surfaces 
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  w_i     !< w-component interpolated onto scalar grid point, required for momentum fluxes at vertical surfaces 

!
!--    Calcuate surface fluxes at horizontal walls
       IF ( .NOT. surf_vertical )  THEN
!
!--       Compute u'w' for the total model domain at upward-facing surfaces.
!--       First compute the corresponding component of u* and square it.
          IF ( .NOT. downward )  THEN
             !$OMP PARALLEL DO PRIVATE( i, j, k, z_mo )
             !$ACC PARALLEL LOOP PRIVATE(i, j, k, z_mo) &
             !$ACC PRESENT(surf, u, rho_air_zw)
             DO  m = 1, surf%ns  
   
                i = surf%i(m)            
                j = surf%j(m)
                k = surf%k(m)

                z_mo = surf%z_mo(m)

                surf%usws(m) = kappa * ( u(k,j,i) - u(k-1,j,i) )               &
                              / ( LOG( z_mo / surf%z0(m) )                     &
                                  - psi_m( z_mo / surf%ol(m) )                 &
                                  + psi_m( surf%z0(m) / surf%ol(m) ) )
!
!--             Please note, the computation of usws is not fully accurate. Actually 
!--             a further interpolation of us onto the u-grid, where usws is defined, 
!--             is required. However, this is not done as this would require several
!--             data transfers between 2D-grid and the surf-type. 
!--             The impact of the missing interpolation is negligible as several 
!--             tests had shown. 
!--             Same also for ol.  
                surf%usws(m) = -surf%usws(m) * surf%us(m) * rho_air_zw(k-1)

             ENDDO
!
!--       At downward-facing surfaces
          ELSE
             !$OMP PARALLEL DO PRIVATE( i, j, k, z_mo )
             DO  m = 1, surf%ns  
   
                i = surf%i(m)            
                j = surf%j(m)
                k = surf%k(m)

                z_mo = surf%z_mo(m)

                surf%usws(m) = kappa * u(k,j,i) / LOG( z_mo / surf%z0(m) )
                surf%usws(m) = surf%usws(m) * surf%us(m) * rho_air_zw(k)

             ENDDO     
          ENDIF

!
!--       Compute v'w' for the total model domain.
!--       First compute the corresponding component of u* and square it.
!--       Upward-facing surfaces
          IF ( .NOT. downward )  THEN
             !$OMP PARALLEL DO PRIVATE( i, j, k, z_mo )
             !$ACC PARALLEL LOOP PRIVATE(i, j, k, z_mo) &
             !$ACC PRESENT(surf, v, rho_air_zw)
             DO  m = 1, surf%ns  
                i = surf%i(m)            
                j = surf%j(m)
                k = surf%k(m)

                z_mo = surf%z_mo(m)

                surf%vsws(m) = kappa * ( v(k,j,i) - v(k-1,j,i) )               &
                           / ( LOG( z_mo / surf%z0(m) )                        &
                               - psi_m( z_mo / surf%ol(m) )                    &
                               + psi_m( surf%z0(m) / surf%ol(m) ) )
!
!--             Please note, the computation of vsws is not fully accurate. Actually 
!--             a further interpolation of us onto the v-grid, where vsws is defined, 
!--             is required. However, this is not done as this would require several
!--             data transfers between 2D-grid and the surf-type. 
!--             The impact of the missing interpolation is negligible as several 
!--             tests had shown. 
!--             Same also for ol.  
                surf%vsws(m) = -surf%vsws(m) * surf%us(m) * rho_air_zw(k-1)
             ENDDO
!
!--       Downward-facing surfaces
          ELSE
             !$OMP PARALLEL DO PRIVATE( i, j, k, z_mo )
             DO  m = 1, surf%ns  
                i = surf%i(m)            
                j = surf%j(m)
                k = surf%k(m)

                z_mo = surf%z_mo(m)

                surf%vsws(m) = kappa * v(k,j,i) / LOG( z_mo / surf%z0(m) )
                surf%vsws(m) = surf%vsws(m) * surf%us(m) * rho_air_zw(k)
             ENDDO
          ENDIF
!
!--       Compute the vertical kinematic heat flux
          IF (  .NOT.  constant_heatflux  .AND.  ( ( time_since_reference_point&
               <=  skip_time_do_lsm  .AND. simulated_time > 0.0_wp ) .OR.      &
               .NOT.  land_surface )  .AND.  .NOT. urban_surface  .AND.        &
               .NOT. downward )  THEN
             !$OMP PARALLEL DO PRIVATE( i, j, k )
             DO  m = 1, surf%ns 
                i    = surf%i(m)            
                j    = surf%j(m)
                k    = surf%k(m)
                surf%shf(m) = -surf%ts(m) * surf%us(m) * rho_air_zw(k-1)
             ENDDO
          ENDIF
!
!--       Compute the vertical water flux
          IF (  .NOT.  constant_waterflux  .AND.  humidity  .AND.              &
               ( ( time_since_reference_point <= skip_time_do_lsm  .AND.       &
               simulated_time > 0.0_wp ) .OR.  .NOT.  land_surface  )  .AND.   &
               .NOT.  urban_surface  .AND.  .NOT. downward )  THEN
             !$OMP PARALLEL DO PRIVATE( i, j, k )
             DO  m = 1, surf%ns 
                i    = surf%i(m)            
                j    = surf%j(m)
                k    = surf%k(m)
                surf%qsws(m) = -surf%qs(m) * surf%us(m) * rho_air_zw(k-1)
             ENDDO
          ENDIF
!
!--       Compute the vertical scalar flux
          IF (  .NOT.  constant_scalarflux  .AND.  passive_scalar  .AND.       &
                .NOT.  downward )  THEN
             !$OMP PARALLEL DO PRIVATE( i, j )
             DO  m = 1, surf%ns   

                i    = surf%i(m)            
                j    = surf%j(m)
                surf%ssws(m) = -surf%ss(m) * surf%us(m)

             ENDDO
          ENDIF   
!
!--       Compute the vertical chemical species' flux
          DO  lsp = 1, nvar
             IF (  .NOT.  constant_csflux(lsp)  .AND.  air_chemistry  .AND.    &
                   .NOT.  downward )  THEN
                !$OMP PARALLEL DO PRIVATE( i, j )
                DO  m = 1, surf%ns   
                   i     = surf%i(m)            
                   j     = surf%j(m)
                   surf%cssws(lsp,m) = -surf%css(lsp,m) * surf%us(m)
                ENDDO
             ENDIF   
          ENDDO

!
!--       Compute (turbulent) fluxes of cloud water content and cloud drop conc.
          IF ( bulk_cloud_model  .AND.  microphysics_morrison  .AND.           &
               .NOT. downward)  THEN
             !$OMP PARALLEL DO PRIVATE( i, j )
             DO  m = 1, surf%ns   

                i    = surf%i(m)            
                j    = surf%j(m)

                surf%qcsws(m) = -surf%qcs(m) * surf%us(m)
                surf%ncsws(m) = -surf%ncs(m) * surf%us(m)
             ENDDO
          ENDIF    
!
!--       Compute (turbulent) fluxes of rain water content and rain drop conc.
          IF ( bulk_cloud_model  .AND.  microphysics_seifert  .AND.            &
               .NOT. downward)  THEN
             !$OMP PARALLEL DO PRIVATE( i, j )
             DO  m = 1, surf%ns   

                i    = surf%i(m)            
                j    = surf%j(m)

                surf%qrsws(m) = -surf%qrs(m) * surf%us(m)
                surf%nrsws(m) = -surf%nrs(m) * surf%us(m)
             ENDDO
          ENDIF

!
!--       Bottom boundary condition for the TKE. 
          IF ( ibc_e_b == 2 )  THEN
             !$OMP PARALLEL DO PRIVATE( i, j, k )
             DO  m = 1, surf%ns   

                i    = surf%i(m)            
                j    = surf%j(m)
                k    = surf%k(m)

                e(k,j,i) = ( surf%us(m) / 0.1_wp )**2
!
!--             As a test: cm = 0.4
!               e(k,j,i) = ( us(j,i) / 0.4_wp )**2
                e(k-1,j,i)   = e(k,j,i)

             ENDDO
          ENDIF
!
!--    Calcuate surface fluxes at vertical surfaces. No stability is considered. 
       ELSE
!
!--       Compute usvs l={0,1} and vsus l={2,3}
          IF ( mom_uv )  THEN
!
!--          Generalize computation by introducing flags. At north- and south-
!--          facing surfaces u-component is used, at east- and west-facing
!--          surfaces v-component is used.
             flag_u = MERGE( 1.0_wp, 0.0_wp, l == 0  .OR.  l == 1 )   
             flag_v = MERGE( 1.0_wp, 0.0_wp, l == 2  .OR.  l == 3 )   
             !$OMP PARALLEL  DO PRIVATE( i, j, k, z_mo )
             DO  m = 1, surf%ns  
                i = surf%i(m)            
                j = surf%j(m)
                k = surf%k(m)

                z_mo = surf%z_mo(m)

                surf%mom_flux_uv(m) = kappa *                                  &
                                ( flag_u * u(k,j,i) + flag_v * v(k,j,i) )  /   &
                                                        LOG( z_mo / surf%z0(m) )

               surf%mom_flux_uv(m) =                                           &
                                    - surf%mom_flux_uv(m) * surf%us(m)
             ENDDO
          ENDIF
!
!--       Compute wsus l={0,1} and wsvs l={2,3}
          IF ( mom_w )  THEN
             !$OMP PARALLEL  DO PRIVATE( i, j, k, z_mo )
             DO  m = 1, surf%ns  
                i = surf%i(m)            
                j = surf%j(m)
                k = surf%k(m)

                z_mo = surf%z_mo(m)

                surf%mom_flux_w(m) = kappa * w(k,j,i) / LOG( z_mo / surf%z0(m) )

                surf%mom_flux_w(m) =                                           &
                                     - surf%mom_flux_w(m) * surf%us(m)
             ENDDO
          ENDIF
!
!--       Compute momentum fluxes used for subgrid-scale TKE production at 
!--       vertical surfaces. In constrast to the calculated momentum fluxes at 
!--       vertical surfaces before, which are defined on the u/v/w-grid, 
!--       respectively), the TKE fluxes are defined at the scalar grid. 
!--       
          IF ( mom_tke )  THEN
!
!--          Precalculate velocity components at scalar grid point. 
             ALLOCATE( u_i(1:surf%ns) )
             ALLOCATE( v_i(1:surf%ns) )
             ALLOCATE( w_i(1:surf%ns) )

             IF ( l == 0  .OR.  l == 1 )  THEN
                !$OMP PARALLEL DO PRIVATE( i, j, k )
                DO  m = 1, surf%ns  
                   i = surf%i(m)            
                   j = surf%j(m)
                   k = surf%k(m)

                   u_i(m) = 0.5_wp * ( u(k,j,i) + u(k,j,i+1) )
                   v_i(m) = 0.0_wp
                   w_i(m) = 0.5_wp * ( w(k,j,i) + w(k-1,j,i) )
                ENDDO
             ELSE
                !$OMP PARALLEL DO PRIVATE( i, j, k )
                DO  m = 1, surf%ns  
                   i = surf%i(m)            
                   j = surf%j(m)
                   k = surf%k(m)

                   u_i(m) = 0.0_wp
                   v_i(m) = 0.5_wp * ( v(k,j,i) + v(k,j+1,i) )
                   w_i(m) = 0.5_wp * ( w(k,j,i) + w(k-1,j,i) )
                ENDDO
             ENDIF

             !$OMP PARALLEL DO PRIVATE( i, j, dum, z_mo )
             DO  m = 1, surf%ns  
                i = surf%i(m)            
                j = surf%j(m)

                z_mo = surf%z_mo(m)

                dum = kappa / LOG( z_mo / surf%z0(m) )
!
!--             usvs (l=0,1) and vsus (l=2,3)
                surf%mom_flux_tke(0,m) = dum * ( u_i(m) + v_i(m) )
!
!--             wsvs (l=0,1) and wsus (l=2,3)
                surf%mom_flux_tke(1,m) = dum * w_i(m)

                surf%mom_flux_tke(0:1,m) =                                     &
                               - surf%mom_flux_tke(0:1,m) * surf%us(m)
             ENDDO
!
!--          Deallocate temporary arrays
             DEALLOCATE( u_i )             
             DEALLOCATE( v_i )  
             DEALLOCATE( w_i )  
          ENDIF
       ENDIF

    END SUBROUTINE calc_surface_fluxes

    
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculates temperature near surface (10 cm) for indoor model or 2 m 
!> temperature for output
!------------------------------------------------------------------------------!
    SUBROUTINE calc_pt_near_surface ( z_char )

       IMPLICIT NONE

       CHARACTER (LEN = *), INTENT(IN) :: z_char !< string identifier to identify z level
       INTEGER(iwp)                    :: i      !< grid index x-dimension
       INTEGER(iwp)                    :: j      !< grid index y-dimension 
       INTEGER(iwp)                    :: k      !< grid index z-dimension
       INTEGER(iwp)                    :: m      !< running index for surface elements

       
       SELECT CASE ( z_char)
            
       
          CASE ( '10cm' )

             DO  m = 1, surf%ns

                i = surf%i(m)
                j = surf%j(m)
                k = surf%k(m)

                surf%pt_10cm(m) = surf%pt_surface(m) + surf%ts(m) / kappa      &
                                   * ( LOG( 0.1_wp /  surf%z0h(m) )            &
                                     - psi_h( 0.1_wp / surf%ol(m) )            &
                                     + psi_h( surf%z0h(m) / surf%ol(m) ) )

             ENDDO

       END SELECT

    END SUBROUTINE calc_pt_near_surface
    

!
!-- Integrated stability function for momentum
    FUNCTION psi_m( zeta ) 
       !$ACC ROUTINE SEQ
       
       USE kinds

       IMPLICIT NONE 

       REAL(wp)            :: psi_m !< Integrated similarity function result
       REAL(wp)            :: zeta  !< Stability parameter z/L
       REAL(wp)            :: x     !< dummy variable

       REAL(wp), PARAMETER :: a = 1.0_wp            !< constant
       REAL(wp), PARAMETER :: b = 0.66666666666_wp  !< constant
       REAL(wp), PARAMETER :: c = 5.0_wp            !< constant
       REAL(wp), PARAMETER :: d = 0.35_wp           !< constant
       REAL(wp), PARAMETER :: c_d_d = c / d         !< constant
       REAL(wp), PARAMETER :: bc_d_d = b * c / d    !< constant


       IF ( zeta < 0.0_wp )  THEN
          x = SQRT( SQRT( 1.0_wp  - 16.0_wp * zeta ) )
          psi_m = pi * 0.5_wp - 2.0_wp * ATAN( x ) + LOG( ( 1.0_wp + x )**2    &
                  * ( 1.0_wp + x**2 ) * 0.125_wp )
       ELSE

          psi_m = - b * ( zeta - c_d_d ) * EXP( -d * zeta ) - a * zeta         &
                   - bc_d_d
!
!--       Old version for stable conditions (only valid for z/L < 0.5)
!--       psi_m = - 5.0_wp * zeta

       ENDIF

    END FUNCTION psi_m


!
!-- Integrated stability function for heat and moisture
    FUNCTION psi_h( zeta ) 
       !$ACC ROUTINE SEQ
       
       USE kinds

       IMPLICIT NONE 

       REAL(wp)            :: psi_h !< Integrated similarity function result
       REAL(wp)            :: zeta  !< Stability parameter z/L
       REAL(wp)            :: x     !< dummy variable

       REAL(wp), PARAMETER :: a = 1.0_wp            !< constant
       REAL(wp), PARAMETER :: b = 0.66666666666_wp  !< constant
       REAL(wp), PARAMETER :: c = 5.0_wp            !< constant
       REAL(wp), PARAMETER :: d = 0.35_wp           !< constant
       REAL(wp), PARAMETER :: c_d_d = c / d         !< constant
       REAL(wp), PARAMETER :: bc_d_d = b * c / d    !< constant


       IF ( zeta < 0.0_wp )  THEN
          x = SQRT( 1.0_wp  - 16.0_wp * zeta )
          psi_h = 2.0_wp * LOG( (1.0_wp + x ) / 2.0_wp )
       ELSE
          psi_h = - b * ( zeta - c_d_d ) * EXP( -d * zeta ) - (1.0_wp          &
                  + 0.66666666666_wp * a * zeta )**1.5_wp - bc_d_d             &
                  + 1.0_wp
!
!--       Old version for stable conditions (only valid for z/L < 0.5)
!--       psi_h = - 5.0_wp * zeta
       ENDIF

    END FUNCTION psi_h


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculates stability function for momentum
!>
!> @author Hauke Wurps
!------------------------------------------------------------------------------!
    FUNCTION phi_m( zeta ) 
       !$ACC ROUTINE SEQ
   
       IMPLICIT NONE 
   
       REAL(wp)            :: phi_m         !< Value of the function
       REAL(wp)            :: zeta          !< Stability parameter z/L
   
       REAL(wp), PARAMETER :: a = 16.0_wp   !< constant
       REAL(wp), PARAMETER :: c = 5.0_wp    !< constant
   
       IF ( zeta < 0.0_wp )  THEN
          phi_m = 1.0_wp / SQRT( SQRT( 1.0_wp - a * zeta ) )
       ELSE
          phi_m = 1.0_wp + c * zeta
       ENDIF
   
    END FUNCTION phi_m

 END MODULE surface_layer_fluxes_mod
