!> @file data_output_mask.f90
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
! $Id: data_output_mask.f90 4457 2020-03-11 14:20:43Z raasch $
! use statement for exchange horiz added
! 
! 4444 2020-03-05 15:59:50Z raasch
! bugfix: cpp-directives for serial mode added
!
! 4377 2020-01-15 11:10:51Z gronemeier
! bugfix: set fill value for output according to wall_flags_total_0 for
!         non-terrain following output
!
! 4360 2020-01-07 11:25:50Z suehring
! Introduction of wall_flags_total_0, which currently sets bits based on static
! topography information used in wall_flags_static_0
!
! 4331 2019-12-10 18:25:02Z suehring
! Formatting adjustment
!
! 4329 2019-12-10 15:46:36Z motisi
! Renamed wall_flags_0 to wall_flags_static_0
!
! 4246 2019-09-30 09:27:52Z pavelkrc
! Corrected "Former revisions" section
!
! 4168 2019-08-16 13:50:17Z suehring
! Remove variable grid
!
! 4167 2019-08-16 11:01:48Z suehring
! Changed behaviour of masked output over surface to follow terrain and ignore
! buildings (J.Resler, T.Gronemeier)
!
! 4069 2019-07-01 14:05:51Z Giersch
! Masked output running index mid has been introduced as a local variable to
! avoid runtime error (Loop variable has been modified) in time_integration
!
! 4039 2019-06-18 10:32:41Z suehring
! Modularize diagnostic output
!
! 3994 2019-05-22 18:08:09Z suehring
! output of turbulence intensity added
!
! 3665 2019-01-10 08:28:24Z raasch
! unused variables removed
!
! 3655 2019-01-07 16:51:22Z knoop
! Fix output time levels (use time_since_reference_point)
!
! 410 2009-12-04 17:05:40Z letzel
! Initial version
!
! Description:
! ------------
!> Masked data output in netCDF format for current mask (current value of mid).
!------------------------------------------------------------------------------!
 SUBROUTINE data_output_mask( av, mid )



#if defined( __netcdf )
    USE arrays_3d,                                                             &
        ONLY:  e, nc, nr, p, pt, q, qc, ql, ql_c, ql_v, qr, rho_ocean, s, sa,  &
               tend, u, v, vpt, w, d_exner

    USE averaging,                                                             &
        ONLY:  e_av, lpt_av, nc_av, nr_av, p_av, pc_av, pr_av, pt_av, q_av,    &
               qc_av, ql_av, ql_c_av, ql_v_av, ql_vp_av, qv_av, qr_av,         &
               rho_ocean_av, s_av, sa_av, u_av, v_av, vpt_av, w_av

    USE basic_constants_and_equations_mod,                                     &
        ONLY:  lv_d_cp

    USE chemistry_model_mod,                                                   &
        ONLY:  chem_data_output_mask

    USE control_parameters,                                                    &
        ONLY:  air_chemistry, domask, domask_no, domask_time_count, mask_i,    &
               mask_j, mask_k, mask_size_l, mask_surface,                                                   &
               max_masks, message_string, nz_do3d, salsa,                      &
               time_since_reference_point

#if defined( __parallel )
    USE control_parameters,                                                    &
        ONLY:  mask_size, mask_start_l
#endif

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point

    USE diagnostic_output_quantities_mod,                                      &
        ONLY:  doq_output_mask

    USE exchange_horiz_mod,                                                    &
        ONLY:  exchange_horiz

    USE indices,                                                               &
        ONLY:  nbgp, nxl, nxr, nyn, nys, nzb, nzt, wall_flags_total_0

    USE kinds

    USE bulk_cloud_model_mod,                                                  &
        ONLY:  bulk_cloud_model

    USE NETCDF

    USE netcdf_interface,                                                      &
        ONLY:  fill_value, id_set_mask, id_var_domask, id_var_time_mask,       &
               nc_stat, netcdf_data_format, netcdf_handle_error

    USE particle_attributes,                                                   &
        ONLY:  grid_particles, number_of_particles, particles,                 &
               particle_advection_start, prt_count

    USE pegrid

    USE radiation_model_mod,                                                   &
        ONLY:  radiation, radiation_data_output_mask

    USE salsa_mod,                                                             &
        ONLY:  salsa_data_output_mask


    IMPLICIT NONE

    INTEGER(iwp) ::  av                      !< flag for (non-)average output
    INTEGER(iwp) ::  flag_nr                 !< number of masking flag
    INTEGER(iwp) ::  i                       !< loop index
    INTEGER(iwp) ::  ivar                    !< variable index
    INTEGER(iwp) ::  j                       !< loop index
    INTEGER(iwp) ::  k                       !< loop index
    INTEGER(iwp) ::  im                      !< loop index for masked variables
    INTEGER(iwp) ::  jm                      !< loop index for masked variables
    INTEGER(iwp) ::  kk                      !< vertical index
    INTEGER(iwp) ::  mid                     !< masked output running index
    INTEGER(iwp) ::  n                       !< loop index
    INTEGER(iwp) ::  netcdf_data_format_save !< value of netcdf_data_format
    INTEGER(iwp) ::  ktt                     !< k index of highest terrain surface
#if defined( __parallel )
    INTEGER(iwp) ::  ngp                     !< number of grid points of an output slice
    INTEGER(iwp) ::  sender                  !< PE id of sending PE
    INTEGER(iwp) ::  ind(6)                  !< index limits (lower/upper bounds) of array 'local_2d'
#endif

    LOGICAL ::  found      !< true if output variable was found
    LOGICAL ::  resorted   !< true if variable is resorted

    REAL(wp) ::  mean_r    !< mean particle radius
    REAL(wp) ::  s_r2      !< sum( particle-radius**2 )
    REAL(wp) ::  s_r3      !< sum( particle-radius**3 )

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  local_pf    !< output array
#if defined( __parallel )
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  total_pf    !< collected output array
#endif
    REAL(wp), DIMENSION(:,:,:), POINTER ::  to_be_resorted  !< points to array which shall be output

!
!-- Return, if nothing to output
    IF ( domask_no(mid,av) == 0 )  RETURN

    CALL cpu_log (log_point(49),'data_output_mask','start')

!
!-- Parallel netcdf output is not tested so far for masked data, hence
!-- netcdf_data_format is switched back to non-paralell output.
    netcdf_data_format_save = netcdf_data_format
    IF ( netcdf_data_format == 5 ) netcdf_data_format = 3
    IF ( netcdf_data_format == 6 ) netcdf_data_format = 4

!
!-- Open output file.
    IF ( myid == 0  .OR.  netcdf_data_format > 4 )  THEN
       CALL check_open( 200+mid+av*max_masks )
    ENDIF

!
!-- Allocate total and local output arrays.
#if defined( __parallel )
    IF ( myid == 0 )  THEN
       ALLOCATE( total_pf(mask_size(mid,1),mask_size(mid,2),mask_size(mid,3)) )
    ENDIF
#endif
    ALLOCATE( local_pf(mask_size_l(mid,1),mask_size_l(mid,2), &
                       mask_size_l(mid,3)) )

!
!-- Update the netCDF time axis.
    domask_time_count(mid,av) = domask_time_count(mid,av) + 1
    IF ( myid == 0  .OR.  netcdf_data_format > 4 )  THEN
       nc_stat = NF90_PUT_VAR( id_set_mask(mid,av), id_var_time_mask(mid,av), &
                               (/ time_since_reference_point /),              &
                               start = (/ domask_time_count(mid,av) /),       &
                               count = (/ 1 /) )
       CALL netcdf_handle_error( 'data_output_mask', 460 )
    ENDIF

!
!-- Loop over all variables to be written.
    ivar = 1

    DO  WHILE ( domask(mid,av,ivar)(1:1) /= ' ' )
!
!--    Reallocate local_pf on PE 0 since its shape changes during MPI exchange
       IF ( netcdf_data_format < 5   .AND.  myid == 0  .AND.  ivar > 1 )  THEN
          DEALLOCATE( local_pf )
          ALLOCATE( local_pf(mask_size_l(mid,1),mask_size_l(mid,2), &
                             mask_size_l(mid,3)) )
       ENDIF
!
!--    Set masking flag for topography for not resorted arrays
       flag_nr = 0
!
!--    Store the variable chosen.
       resorted = .FALSE.
       SELECT CASE ( TRIM( domask(mid,av,ivar) ) )

          CASE ( 'e' )
             IF ( av == 0 )  THEN
                to_be_resorted => e
             ELSE
                to_be_resorted => e_av
             ENDIF

          CASE ( 'thetal' )
             IF ( av == 0 )  THEN
                to_be_resorted => pt
             ELSE
                to_be_resorted => lpt_av
             ENDIF

          CASE ( 'nc' )
             IF ( av == 0 )  THEN
                to_be_resorted => nc
             ELSE
                to_be_resorted => nc_av
             ENDIF

          CASE ( 'nr' )
             IF ( av == 0 )  THEN
                to_be_resorted => nr
             ELSE
                to_be_resorted => nr_av
             ENDIF

          CASE ( 'p' )
             IF ( av == 0 )  THEN
                to_be_resorted => p
             ELSE
                to_be_resorted => p_av
             ENDIF

          CASE ( 'pc' )  ! particle concentration (requires ghostpoint exchange)
             IF ( av == 0 )  THEN
                tend = prt_count
                CALL exchange_horiz( tend, nbgp )
                IF ( .NOT. mask_surface(mid) )  THEN
                   DO  i = 1, mask_size_l(mid,1)
                      DO  j = 1, mask_size_l(mid,2)
                         DO  k = 1, mask_size_l(mid,3)
                            local_pf(i,j,k) =  tend(mask_k(mid,k), &
                                      mask_j(mid,j),mask_i(mid,i))
                         ENDDO
                      ENDDO
                   ENDDO
                ELSE
!
!--                Terrain-following masked output
                   DO  i = 1, mask_size_l(mid,1)
                      DO  j = 1, mask_size_l(mid,2)
!--                      Get k index of the highest terraing surface
                         im = mask_i(mid,i)
                         jm = mask_j(mid,j)
                         ktt = MINLOC( MERGE( 1, 0, BTEST( wall_flags_total_0(:,jm,im), 5 )),&
                                       DIM = 1 ) - 1
                         DO  k = 1, mask_size_l(mid,3)
                            kk = MIN( ktt+mask_k(mid,k), nzt+1 )
!--                         Set value if not in building
                            IF ( BTEST( wall_flags_total_0(kk,jm,im), 6 ) )  THEN
                               local_pf(i,j,k) = fill_value
                            ELSE
                               local_pf(i,j,k) =  tend(kk,jm,im)
                            ENDIF
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF
                resorted = .TRUE.
             ELSE
                CALL exchange_horiz( pc_av, nbgp )
                to_be_resorted => pc_av
             ENDIF

          CASE ( 'pr' )  ! mean particle radius (effective radius)
             IF ( av == 0 )  THEN
                IF ( time_since_reference_point >= particle_advection_start )  THEN
                   DO  i = nxl, nxr
                      DO  j = nys, nyn
                         DO  k = nzb, nz_do3d
                            number_of_particles = prt_count(k,j,i)
                            IF (number_of_particles <= 0)  CYCLE
                            particles => grid_particles(k,j,i)%particles(1:number_of_particles)
                            s_r2 = 0.0_wp
                            s_r3 = 0.0_wp
                            DO  n = 1, number_of_particles
                               IF ( particles(n)%particle_mask )  THEN
                                  s_r2 = s_r2 + grid_particles(k,j,i)%particles(n)%radius**2 * &
                                         grid_particles(k,j,i)%particles(n)%weight_factor
                                  s_r3 = s_r3 + grid_particles(k,j,i)%particles(n)%radius**3 * &
                                         grid_particles(k,j,i)%particles(n)%weight_factor
                               ENDIF
                            ENDDO
                            IF ( s_r2 > 0.0_wp )  THEN
                               mean_r = s_r3 / s_r2
                            ELSE
                               mean_r = 0.0_wp
                            ENDIF
                            tend(k,j,i) = mean_r
                         ENDDO
                      ENDDO
                   ENDDO
                   CALL exchange_horiz( tend, nbgp )
                ELSE
                   tend = 0.0_wp
                ENDIF
                IF ( .NOT. mask_surface(mid) )  THEN
                   DO  i = 1, mask_size_l(mid,1)
                      DO  j = 1, mask_size_l(mid,2)
                         DO  k = 1, mask_size_l(mid,3)
                            local_pf(i,j,k) =  tend(mask_k(mid,k), &
                                      mask_j(mid,j),mask_i(mid,i))
                         ENDDO
                      ENDDO
                   ENDDO
                ELSE
!
!--                Terrain-following masked output
                   DO  i = 1, mask_size_l(mid,1)
                      DO  j = 1, mask_size_l(mid,2)
!--                      Get k index of the highest terraing surface
                         im = mask_i(mid,i)
                         jm = mask_j(mid,j)
                         ktt = MINLOC( MERGE( 1, 0, BTEST( wall_flags_total_0(:,jm,im), 5 )), &
                                         DIM = 1 ) - 1
                         DO  k = 1, mask_size_l(mid,3)
                            kk = MIN( ktt+mask_k(mid,k), nzt+1 )
!--                         Set value if not in building
                            IF ( BTEST( wall_flags_total_0(kk,jm,im), 6 ) )  THEN
                               local_pf(i,j,k) = fill_value
                            ELSE
                               local_pf(i,j,k) =  tend(kk,jm,im)
                            ENDIF
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF
                resorted = .TRUE.
             ELSE
                CALL exchange_horiz( pr_av, nbgp )
                to_be_resorted => pr_av
             ENDIF

          CASE ( 'theta' )
             IF ( av == 0 )  THEN
                IF ( .NOT. bulk_cloud_model ) THEN
                   to_be_resorted => pt
                ELSE
                   IF ( .NOT. mask_surface(mid) )  THEN
                      DO  i = 1, mask_size_l(mid,1)
                         DO  j = 1, mask_size_l(mid,2)
                            DO  k = 1, mask_size_l(mid,3)
                               local_pf(i,j,k) =  &
                                  pt(mask_k(mid,k),mask_j(mid,j),mask_i(mid,i)) &
                                  + lv_d_cp * d_exner(mask_k(mid,k)) *          &
                                    ql(mask_k(mid,k),mask_j(mid,j),mask_i(mid,i))
                            ENDDO
                         ENDDO
                      ENDDO
                   ELSE
!
!--                   Terrain-following masked output
                      DO  i = 1, mask_size_l(mid,1)
                         DO  j = 1, mask_size_l(mid,2)
!--                         Get k index of the highest terraing surface
                            im = mask_i(mid,i)
                            jm = mask_j(mid,j)
                            ktt = MINLOC( MERGE( 1, 0, BTEST( wall_flags_total_0(:,jm,im), 5 )), &
                                             DIM = 1 ) - 1
                            DO  k = 1, mask_size_l(mid,3)
                               kk = MIN( ktt+mask_k(mid,k), nzt+1 )
!--                            Set value if not in building
                               IF ( BTEST( wall_flags_total_0(kk,jm,im), 6 ) )  THEN
                                  local_pf(i,j,k) = fill_value
                               ELSE
                                  local_pf(i,j,k) = pt(kk,jm,im) + lv_d_cp * d_exner(kk) * ql(kk,jm,im)
                               ENDIF
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDIF
                   resorted = .TRUE.
                ENDIF
             ELSE
                to_be_resorted => pt_av
             ENDIF

          CASE ( 'q' )
             IF ( av == 0 )  THEN
                to_be_resorted => q
             ELSE
                to_be_resorted => q_av
             ENDIF

          CASE ( 'qc' )
             IF ( av == 0 )  THEN
                to_be_resorted => qc
             ELSE
                to_be_resorted => qc_av
             ENDIF

          CASE ( 'ql' )
             IF ( av == 0 )  THEN
                to_be_resorted => ql
             ELSE
                to_be_resorted => ql_av
             ENDIF

          CASE ( 'ql_c' )
             IF ( av == 0 )  THEN
                to_be_resorted => ql_c
             ELSE
                to_be_resorted => ql_c_av
             ENDIF

          CASE ( 'ql_v' )
             IF ( av == 0 )  THEN
                to_be_resorted => ql_v
             ELSE
                to_be_resorted => ql_v_av
             ENDIF

          CASE ( 'ql_vp' )
             IF ( av == 0 )  THEN
                IF ( time_since_reference_point >= particle_advection_start )  THEN
                   DO  i = nxl, nxr
                      DO  j = nys, nyn
                         DO  k = nzb, nz_do3d
                            number_of_particles = prt_count(k,j,i)
                            IF (number_of_particles <= 0)  CYCLE
                            particles => grid_particles(k,j,i)%particles(1:number_of_particles)
                            DO  n = 1, number_of_particles
                               IF ( particles(n)%particle_mask )  THEN
                                  tend(k,j,i) = tend(k,j,i) + &
                                          particles(n)%weight_factor / &
                                          prt_count(k,j,i)
                               ENDIF
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDDO
                   CALL exchange_horiz( tend, nbgp )
                ELSE
                   tend = 0.0_wp
                ENDIF
                IF ( .NOT. mask_surface(mid) )  THEN
                   DO  i = 1, mask_size_l(mid,1)
                      DO  j = 1, mask_size_l(mid,2)
                         DO  k = 1, mask_size_l(mid,3)
                            local_pf(i,j,k) =  tend(mask_k(mid,k), &
                                      mask_j(mid,j),mask_i(mid,i))
                         ENDDO
                      ENDDO
                   ENDDO
                ELSE
!
!--                Terrain-following masked output
                   DO  i = 1, mask_size_l(mid,1)
                      DO  j = 1, mask_size_l(mid,2)
!--                      Get k index of the highest terraing surface
                         im = mask_i(mid,i)
                         jm = mask_j(mid,j)
                         ktt = MINLOC( MERGE( 1, 0, BTEST( wall_flags_total_0(:,jm,im), 5 )), &
                                          DIM = 1 ) - 1
                         DO  k = 1, mask_size_l(mid,3)
                            kk = MIN( ktt+mask_k(mid,k), nzt+1 )
!--                         Set value if not in building
                            IF ( BTEST( wall_flags_total_0(kk,jm,im), 6 ) )  THEN
                               local_pf(i,j,k) = fill_value
                            ELSE
                               local_pf(i,j,k) = tend(kk,jm,im)
                            ENDIF
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF
                resorted = .TRUE.
             ELSE
                CALL exchange_horiz( ql_vp_av, nbgp )
                to_be_resorted => ql_vp_av
             ENDIF

          CASE ( 'qv' )
             IF ( av == 0 )  THEN
                IF ( .NOT. mask_surface(mid) )  THEN
                   DO  i = 1, mask_size_l(mid,1)
                      DO  j = 1, mask_size_l(mid,2)
                         DO  k = 1, mask_size_l(mid,3)
                            local_pf(i,j,k) =  &
                                 q(mask_k(mid,k),mask_j(mid,j),mask_i(mid,i)) -  &
                                 ql(mask_k(mid,k),mask_j(mid,j),mask_i(mid,i))
                         ENDDO
                      ENDDO
                   ENDDO
                ELSE
!
!--                Terrain-following masked output
                   DO  i = 1, mask_size_l(mid,1)
                      DO  j = 1, mask_size_l(mid,2)
!--                      Get k index of the highest terraing surface
                         im = mask_i(mid,i)
                         jm = mask_j(mid,j)
                         ktt = MINLOC( MERGE( 1, 0, BTEST( wall_flags_total_0(:,jm,im), 5 )), &
                                          DIM = 1 ) - 1
                         DO  k = 1, mask_size_l(mid,3)
                            kk = MIN( ktt+mask_k(mid,k), nzt+1 )
!--                         Set value if not in building
                            IF ( BTEST( wall_flags_total_0(kk,jm,im), 6 ) )  THEN
                               local_pf(i,j,k) = fill_value
                            ELSE
                               local_pf(i,j,k) = q(kk,jm,im) - ql(kk,jm,im)
                            ENDIF
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF
                resorted = .TRUE.
             ELSE
                to_be_resorted => qv_av
             ENDIF

          CASE ( 'qr' )
             IF ( av == 0 )  THEN
                to_be_resorted => qr
             ELSE
                to_be_resorted => qr_av
             ENDIF

          CASE ( 'rho_sea_water' )
             IF ( av == 0 )  THEN
                to_be_resorted => rho_ocean
             ELSE
                to_be_resorted => rho_ocean_av
             ENDIF

          CASE ( 's' )
             IF ( av == 0 )  THEN
                to_be_resorted => s
             ELSE
                to_be_resorted => s_av
             ENDIF

          CASE ( 'sa' )
             IF ( av == 0 )  THEN
                to_be_resorted => sa
             ELSE
                to_be_resorted => sa_av
             ENDIF

          CASE ( 'u' )
             flag_nr = 1
             IF ( av == 0 )  THEN
                to_be_resorted => u
             ELSE
                to_be_resorted => u_av
             ENDIF

          CASE ( 'v' )
             flag_nr = 2
             IF ( av == 0 )  THEN
                to_be_resorted => v
             ELSE
                to_be_resorted => v_av
             ENDIF

          CASE ( 'thetav' )
             IF ( av == 0 )  THEN
                to_be_resorted => vpt
             ELSE
                to_be_resorted => vpt_av
             ENDIF

          CASE ( 'w' )
             flag_nr = 3
             IF ( av == 0 )  THEN
                to_be_resorted => w
             ELSE
                to_be_resorted => w_av
             ENDIF

          CASE DEFAULT
!
!--          Set flag to steer output of radiation, land-surface, or user-defined
!--          quantities
             found = .FALSE.
!
!--          Radiation quantity
             IF ( .NOT. found  .AND. radiation )  THEN
                CALL radiation_data_output_mask(av, domask(mid,av,ivar), found,&
                                                local_pf, mid )
             ENDIF

             IF ( .NOT. found  .AND. air_chemistry )  THEN
                CALL chem_data_output_mask(av, domask(mid,av,ivar), found,     &
                                           local_pf, mid )
             ENDIF
!
!--          Check for diagnostic quantities
             IF ( .NOT. found )  THEN
                CALL doq_output_mask( av, domask(mid,av,ivar), found, local_pf,   &
                                      mid)
             ENDIF
!
!--          SALSA quantities
             IF ( .NOT. found .AND. salsa )  THEN
                CALL salsa_data_output_mask( av, domask(mid,av,ivar), found,   &
                                             local_pf, mid )
             ENDIF
!
!--          User defined quantity
             IF ( .NOT. found )  THEN
                CALL user_data_output_mask(av, domask(mid,av,ivar), found,     &
                                           local_pf, mid )
             ENDIF

             resorted = .TRUE.

             IF ( .NOT. found )  THEN
                WRITE ( message_string, * ) 'no masked output available for: ',&
                                            TRIM( domask(mid,av,ivar) )
                CALL message( 'data_output_mask', 'PA0327', 0, 0, 0, 6, 0 )
             ENDIF

       END SELECT

!
!--    Resort the array to be output, if not done above
       IF ( .NOT. resorted )  THEN
          IF ( .NOT. mask_surface(mid) )  THEN
!
!--          Default masked output
             DO  i = 1, mask_size_l(mid,1)
                DO  j = 1, mask_size_l(mid,2)
                   DO  k = 1, mask_size_l(mid,3)
                      local_pf(i,j,k) = MERGE( to_be_resorted(mask_k(mid,k),  &
                                                              mask_j(mid,j),  &
                                                              mask_i(mid,i)), &
                                               REAL( fill_value, KIND = wp ), &
                                               BTEST( wall_flags_total_0(     &
                                                              mask_k(mid,k),  &
                                                              mask_j(mid,j),  &
                                                              mask_i(mid,i)), &
                                                      flag_nr ) )
                   ENDDO
                ENDDO
             ENDDO

          ELSE
!
!--          Terrain-following masked output
             DO  i = 1, mask_size_l(mid,1)
                DO  j = 1, mask_size_l(mid,2)
!--                Get k index of the highest terraing surface
                   im = mask_i(mid,i)
                   jm = mask_j(mid,j)
                   ktt = MINLOC( MERGE( 1, 0, BTEST( wall_flags_total_0(:,jm,im), 5 )), &
                                    DIM = 1 ) - 1
                   DO  k = 1, mask_size_l(mid,3)
                      kk = MIN( ktt+mask_k(mid,k), nzt+1 )
!--                   Set value if not in building
                      IF ( BTEST( wall_flags_total_0(kk,jm,im), 6 ) )  THEN
                         local_pf(i,j,k) = fill_value
                      ELSE
                         local_pf(i,j,k) = to_be_resorted(kk,jm,im)
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO

          ENDIF
       ENDIF

!
!--    I/O block. I/O methods are implemented
!--    (1) for parallel execution
!--     a. with netCDF 4 parallel I/O-enabled library
!--     b. with netCDF 3 library
!--    (2) for serial execution.
!--    The choice of method depends on the correct setting of preprocessor
!--    directives __parallel and __netcdf4_parallel as well as on the parameter
!--    netcdf_data_format.
#if defined( __parallel )
#if defined( __netcdf4_parallel )
       IF ( netcdf_data_format > 4 )  THEN
!
!--       (1) a. Parallel I/O using netCDF 4 (not yet tested)
          nc_stat = NF90_PUT_VAR( id_set_mask(mid,av),                         &
               id_var_domask(mid,av,ivar), local_pf,                           &
               start = (/ mask_start_l(mid,1), mask_start_l(mid,2),            &
                          mask_start_l(mid,3), domask_time_count(mid,av) /),   &
               count = (/ mask_size_l(mid,1), mask_size_l(mid,2),              &
                          mask_size_l(mid,3), 1 /) )
          CALL netcdf_handle_error( 'data_output_mask', 461 )
       ELSE
#endif
!
!--       (1) b. Conventional I/O only through PE0
!--       PE0 receives partial arrays from all processors of the respective mask
!--       and outputs them. Here a barrier has to be set, because otherwise
!--       "-MPI- FATAL: Remote protocol queue full" may occur.
          CALL MPI_BARRIER( comm2d, ierr )

          ngp = mask_size_l(mid,1) * mask_size_l(mid,2) * mask_size_l(mid,3)
          IF ( myid == 0 )  THEN
!
!--          Local array can be relocated directly.
             total_pf( &
               mask_start_l(mid,1):mask_start_l(mid,1)+mask_size_l(mid,1)-1, &
               mask_start_l(mid,2):mask_start_l(mid,2)+mask_size_l(mid,2)-1, &
               mask_start_l(mid,3):mask_start_l(mid,3)+mask_size_l(mid,3)-1 ) &
               = local_pf
!
!--          Receive data from all other PEs.
             DO  n = 1, numprocs-1
!
!--             Receive index limits first, then array.
!--             Index limits are received in arbitrary order from the PEs.
                CALL MPI_RECV( ind(1), 6, MPI_INTEGER, MPI_ANY_SOURCE, 0,  &
                     comm2d, status, ierr )
!
!--             Not all PEs have data for the mask
                IF ( ind(1) /= -9999 )  THEN
                   ngp = ( ind(2)-ind(1)+1 ) * (ind(4)-ind(3)+1 ) *  &
                         ( ind(6)-ind(5)+1 )
                   sender = status(MPI_SOURCE)
                   DEALLOCATE( local_pf )
                   ALLOCATE(local_pf(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)))
                   CALL MPI_RECV( local_pf(ind(1),ind(3),ind(5)), ngp,  &
                        MPI_REAL, sender, 1, comm2d, status, ierr )
                   total_pf(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)) &
                        = local_pf
                ENDIF
             ENDDO

             nc_stat = NF90_PUT_VAR( id_set_mask(mid,av),                      &
                  id_var_domask(mid,av,ivar), total_pf,                        &
                  start = (/ 1, 1, 1, domask_time_count(mid,av) /),            &
                  count = (/ mask_size(mid,1), mask_size(mid,2),               &
                             mask_size(mid,3), 1 /) )
             CALL netcdf_handle_error( 'data_output_mask', 462 )

          ELSE
!
!--          If at least part of the mask resides on the PE, send the index
!--          limits for the target array, otherwise send -9999 to PE0.
             IF ( mask_size_l(mid,1) > 0 .AND.  mask_size_l(mid,2) > 0 .AND. &
                  mask_size_l(mid,3) > 0  ) &
                  THEN
                ind(1) = mask_start_l(mid,1)
                ind(2) = mask_start_l(mid,1) + mask_size_l(mid,1) - 1
                ind(3) = mask_start_l(mid,2)
                ind(4) = mask_start_l(mid,2) + mask_size_l(mid,2) - 1
                ind(5) = mask_start_l(mid,3)
                ind(6) = mask_start_l(mid,3) + mask_size_l(mid,3) - 1
             ELSE
                ind(1) = -9999; ind(2) = -9999
                ind(3) = -9999; ind(4) = -9999
                ind(5) = -9999; ind(6) = -9999
             ENDIF
             CALL MPI_SEND( ind(1), 6, MPI_INTEGER, 0, 0, comm2d, ierr )
!
!--          If applicable, send data to PE0.
             IF ( ind(1) /= -9999 )  THEN
                CALL MPI_SEND( local_pf(1,1,1), ngp, MPI_REAL, 0, 1, comm2d, &
                     ierr )
             ENDIF
          ENDIF
!
!--       A barrier has to be set, because otherwise some PEs may proceed too
!--       fast so that PE0 may receive wrong data on tag 0.
          CALL MPI_BARRIER( comm2d, ierr )
#if defined( __netcdf4_parallel )
       ENDIF
#endif
#else
!
!--    (2) For serial execution of PALM, the single processor (PE0) holds all
!--    data and writes them directly to file.
       nc_stat = NF90_PUT_VAR( id_set_mask(mid,av),                            &
                               id_var_domask(mid,av,ivar), local_pf,           &
                             start = (/ 1, 1, 1, domask_time_count(mid,av) /), &
                             count = (/ mask_size_l(mid,1), mask_size_l(mid,2),&
                               mask_size_l(mid,3), 1 /) )
       CALL netcdf_handle_error( 'data_output_mask', 463 )
#endif

       ivar = ivar + 1

    ENDDO

!
!-- Deallocate temporary arrays.
    DEALLOCATE( local_pf )
#if defined( __parallel )
    IF ( myid == 0 )  THEN
       DEALLOCATE( total_pf )
    ENDIF
#endif

!
!-- Switch back to original format given by user (see beginning of this routine)
    netcdf_data_format = netcdf_data_format_save

    CALL cpu_log( log_point(49), 'data_output_mask', 'stop' )
#endif


 END SUBROUTINE data_output_mask
