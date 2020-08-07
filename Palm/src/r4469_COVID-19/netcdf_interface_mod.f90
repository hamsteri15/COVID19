!> @file netcdf_interface_mod.f90
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
! $Id: netcdf_interface_mod.f90 4455 2020-03-11 12:20:29Z Giersch $
! Axis attribute added to netcdf output
! 
! 4400 2020-02-10 20:32:41Z suehring
! Move routine to transform coordinates from netcdf_interface_mod to 
! basic_constants_and_equations_mod
! 
! 4360 2020-01-07 11:25:50Z suehring
! Adjusted output of multi-agent system for biometeorology
! 
! 4227 2019-09-10 18:04:34Z gronemeier
! Replace function date_time_string by call to get_date_time
! 
! 4223 2019-09-10 09:20:47Z gronemeier
! replaced rotation angle from input-netCDF file
! by namelist parameter 'rotation_angle'
! 
! 4182 2019-08-22 15:20:23Z scharf
! Corrected "Former revisions" section
! 
! 4127 2019-07-30 14:47:10Z suehring
! -Introduce new vertical dimension for plant-canopy output. 
! -Temporarlily disable masked output for soil (merge from branch resler)
! 
! 4069 2019-07-01 14:05:51Z Giersch
! Masked output running index mid has been introduced as a local variable to 
! avoid runtime error (Loop variable has been modified) in time_integration
! 
! 4046 2019-06-21 17:32:04Z knoop
! removal of special treatment for usm_define_netcdf_grid call
! 
! 4039 2019-06-18 10:32:41Z suehring
! Rename subroutines in module for diagnostic quantities
!
! 4029 2019-06-14 14:04:35Z raasch
! netcdf variable NF90_NOFILL is used as argument instead of "1" in calls to NF90_DEF_VAR_FILL
!
! 3995 2019-05-22 18:59:54Z suehring
! output of turbulence intensity added
!
! 3994 2019-05-22 18:08:09Z suehring
! remove origin time from time unit, compose origin_time_string within
! subroutine netcdf_create_global_atts
!
! 3954 2019-05-06 12:49:42Z gronemeier
! bugfix: corrected format for date_time_string
!
! 3953 2019-05-06 12:11:55Z gronemeier
! bugfix: set origin_time and starting point of time coordinate according to
!         day_of_year_init and time_utc_init
!
! 3942 2019-04-30 13:08:30Z kanani
! Add specifier to netcdf_handle_error to simplify identification of attribute
! causing the error
!
! 3766 2019-02-26 16:23:41Z raasch
! bugfix in im_define_netcdf_grid argument list
!
! 3745 2019-02-15 18:57:56Z suehring
! Add indoor model
!
! 3744 2019-02-15 18:38:58Z suehring
! Bugfix: - initialize return values to ensure they are set before returning
!           (routine define_geo_coordinates)
!         - change order of dimensions for some variables
!
! 3727 2019-02-08 14:52:10Z gronemeier
! make several routines publicly available
!
! 3701 2019-01-26 18:57:21Z knoop
! Statement added to prevent compiler warning about unused variable
!
! 3655 2019-01-07 16:51:22Z knoop
! Move the control parameter "salsa" from salsa_mod to control_parameters
! (M. Kurppa)
!
! Revision 1.1  2005/05/18 15:37:16  raasch
! Initial revision
!
!
! Description:
! ------------
!> In case of extend = .FALSE.:
!> Define all necessary dimensions, axes and variables for the different
!> netCDF datasets. This subroutine is called from check_open after a new
!> dataset is created. It leaves the open netCDF files ready to write.
!>
!> In case of extend = .TRUE.:
!> Find out if dimensions and variables of an existing file match the values
!> of the actual run. If so, get all necessary information (ids, etc.) from
!> this file.
!>
!> Parameter av can assume values 0 (non-averaged data) and 1 (time averaged
!> data)
!>
!> @todo calculation of output time levels for parallel NetCDF still does not
!>       cover every exception (change of dt_do, end_time in restart)
!> @todo timeseries and profile output still needs to be rewritten to allow
!>       modularization
!> @todo output 2d UTM coordinates without global arrays
!> @todo output longitude/latitude also with non-parallel output (3d and xy)
!------------------------------------------------------------------------------!
 MODULE netcdf_interface

    USE control_parameters,                                                    &
        ONLY:  biometeorology, fl_max,                                         &
               max_masks, multi_agent_system_end,                              &
               multi_agent_system_start,                                       &
               rotation_angle,                                                 &
               var_fl_max, varnamelength
    USE kinds
#if defined( __netcdf )
    USE NETCDF
#endif
    USE mas_global_attributes,                                                 &
        ONLY:  dim_size_agtnum

    USE netcdf_data_input_mod,                                                 &
        ONLY: coord_ref_sys,                                                   &
              crs_list,                                                        &
              init_model

    PRIVATE

    CHARACTER (LEN=16), DIMENSION(13) ::  agt_var_names =                      &
          (/ 'ag_id           ', 'ag_x            ', 'ag_y            ',       &
             'ag_wind         ', 'ag_temp         ', 'ag_group        ',       &
             'ag_iPT          ', 'ag_PM10         ', 'ag_PM25         ',       &
             'not_used        ', 'not_used        ', 'not_used        ',       &
             'not_used        ' /)

    CHARACTER (LEN=16), DIMENSION(13) ::  agt_var_units = &
          (/ 'dim_less        ', 'meters          ', 'meters          ',       &
             'm/s             ', 'K               ', 'dim_less        ',       &
             'C               ', 'tbd             ', 'tbd             ',       &
             'tbd             ', 'not_used        ', 'not_used        ',       &
             'not_used        ' /)

    INTEGER(iwp), PARAMETER ::  dopr_norm_num = 7, dopts_num = 29, dots_max = 100

    CHARACTER (LEN=7), DIMENSION(dopr_norm_num) ::  dopr_norm_names =          &
         (/ 'wtheta0', 'ws2    ', 'tsw2   ', 'ws3    ', 'ws2tsw ', 'wstsw2 ',  &
            'z_i    ' /)

    CHARACTER (LEN=7), DIMENSION(dopr_norm_num) ::  dopr_norm_longnames =      &
         (/ 'wtheta0', 'w*2    ', 't*w2   ', 'w*3    ', 'w*2t*w ', 'w*t*w2 ',  &
            'z_i    ' /)

    CHARACTER (LEN=7), DIMENSION(dopts_num) :: dopts_label =                   &
          (/ 'tnpt   ', 'x_     ', 'y_     ', 'z_     ', 'z_abs  ', 'u      ', &
             'v      ', 'w      ', 'u"     ', 'v"     ', 'w"     ', 'npt_up ', &
             'w_up   ', 'w_down ', 'radius ', 'r_min  ', 'r_max  ', 'npt_max', &
             'npt_min', 'x*2    ', 'y*2    ', 'z*2    ', 'u*2    ', 'v*2    ', &
             'w*2    ', 'u"2    ', 'v"2    ', 'w"2    ', 'npt*2  ' /)

    CHARACTER (LEN=7), DIMENSION(dopts_num) :: dopts_unit =                    &
          (/ 'number ', 'm      ', 'm      ', 'm      ', 'm      ', 'm/s    ', &
             'm/s    ', 'm/s    ', 'm/s    ', 'm/s    ', 'm/s    ', 'number ', &
             'm/s    ', 'm/s    ', 'm      ', 'm      ', 'm      ', 'number ', &
             'number ', 'm2     ', 'm2     ', 'm2     ', 'm2/s2  ', 'm2/s2  ', &
             'm2/s2  ', 'm2/s2  ', 'm2/s2  ', 'm2/s2  ', 'number2' /)

    INTEGER(iwp) ::  dots_num  = 25  !< number of timeseries defined by default
    INTEGER(iwp) ::  dots_soil = 26  !< starting index for soil-timeseries
    INTEGER(iwp) ::  dots_rad  = 32  !< starting index for radiation-timeseries

    CHARACTER (LEN=13), DIMENSION(dots_max) :: dots_label =                    &
          (/ 'E            ', 'E*           ', 'dt           ',                &
             'us*          ', 'th*          ', 'umax         ',                &
             'vmax         ', 'wmax         ', 'div_new      ',                &
             'div_old      ', 'zi_wtheta    ', 'zi_theta     ',                &
             'w*           ', 'w"theta"0    ', 'w"theta"     ',                &
             'wtheta       ', 'theta(0)     ', 'theta(z_mo)  ',                &
             'w"u"0        ', 'w"v"0        ', 'w"q"0        ',                &
             'ol           ', 'q*           ', 'w"s"         ',                &
             's*           ', 'ghf          ', 'qsws_liq     ',                &
             'qsws_soil    ', 'qsws_veg     ', 'r_a          ',                &
             'r_s          ',                                                  &
             'rad_net      ', 'rad_lw_in    ', 'rad_lw_out   ',                &
             'rad_sw_in    ', 'rad_sw_out   ', 'rrtm_aldif   ',                &
             'rrtm_aldir   ', 'rrtm_asdif   ', 'rrtm_asdir   ',                &
             ( 'unknown      ', i9 = 1, dots_max-40 ) /)

    CHARACTER (LEN=13), DIMENSION(dots_max) :: dots_unit =                     &
          (/ 'm2/s2        ', 'm2/s2        ', 's            ',                &
             'm/s          ', 'K            ', 'm/s          ',                &
             'm/s          ', 'm/s          ', 's-1          ',                &
             's-1          ', 'm            ', 'm            ',                &
             'm/s          ', 'K m/s        ', 'K m/s        ',                &
             'K m/s        ', 'K            ', 'K            ',                &
             'm2/s2        ', 'm2/s2        ', 'kg m/s       ',                &
             'm            ', 'kg/kg        ', 'kg m/(kg s)  ',                &
             'kg/kg        ', 'W/m2         ', 'W/m2         ',                &
             'W/m2         ', 'W/m2         ', 's/m          ',                &
             's/m          ',                                                  &
             'W/m2         ', 'W/m2         ', 'W/m2         ',                &
             'W/m2         ', 'W/m2         ', '             ',                &
             '             ', '             ', '             ',                &
             ( 'unknown      ', i9 = 1, dots_max-40 ) /)

    CHARACTER (LEN=16) :: heatflux_output_unit     !< unit for heatflux output
    CHARACTER (LEN=16) :: waterflux_output_unit    !< unit for waterflux output
    CHARACTER (LEN=16) :: momentumflux_output_unit !< unit for momentumflux output

    CHARACTER (LEN=9), DIMENSION(300) ::  dopr_unit = 'unknown'

    CHARACTER (LEN=7), DIMENSION(0:1,500) ::  do2d_unit, do3d_unit

!    CHARACTER (LEN=16), DIMENSION(25) ::  prt_var_names = &
!          (/ 'pt_age          ', 'pt_dvrp_size    ', 'pt_origin_x     ', &
!             'pt_origin_y     ', 'pt_origin_z     ', 'pt_radius       ', &
!             'pt_speed_x      ', 'pt_speed_y      ', 'pt_speed_z      ', &
!             'pt_weight_factor', 'pt_x            ', 'pt_y            ', &
!             'pt_z            ', 'pt_color        ', 'pt_group        ', &
!             'pt_tailpoints   ', 'pt_tail_id      ', 'pt_density_ratio', &
!             'pt_exp_arg      ', 'pt_exp_term     ', 'not_used        ', &
!             'not_used        ', 'not_used        ', 'not_used        ', &
!             'not_used        ' /)

!    CHARACTER (LEN=16), DIMENSION(25) ::  prt_var_units = &
!          (/ 'seconds         ', 'meters          ', 'meters          ', &
!             'meters          ', 'meters          ', 'meters          ', &
!             'm/s             ', 'm/s             ', 'm/s             ', &
!             'factor          ', 'meters          ', 'meters          ', &
!             'meters          ', 'none            ', 'none            ', &
!             'none            ', 'none            ', 'ratio           ', &
!             'none            ', 'none            ', 'not_used        ', &
!             'not_used        ', 'not_used        ', 'not_used        ', &
!             'not_used        ' /)

    CHARACTER(LEN=20), DIMENSION(11) ::  netcdf_precision = ' '
    CHARACTER(LEN=40) ::  netcdf_data_format_string

    INTEGER(iwp) ::  id_dim_agtnum, id_dim_time_agt,                           &
                     id_dim_time_fl, id_dim_time_pr,                           &
                     id_dim_time_pts, id_dim_time_sp, id_dim_time_ts,          &
                     id_dim_x_sp, id_dim_y_sp, id_dim_zu_sp, id_dim_zw_sp,     &
                     id_set_agt, id_set_fl, id_set_pr, id_set_prt, id_set_pts, &
                     id_set_sp, id_set_ts, id_var_agtnum, id_var_time_agt,     &
                     id_var_time_fl, id_var_rnoa_agt, id_var_time_pr,          &
                     id_var_time_pts, id_var_time_sp, id_var_time_ts,          &
                     id_var_x_sp, id_var_y_sp, id_var_zu_sp, id_var_zw_sp,     &
                     nc_stat


    INTEGER(iwp), DIMENSION(0:1) ::  id_dim_time_xy, id_dim_time_xz, &
                    id_dim_time_yz, id_dim_time_3d, id_dim_x_xy, id_dim_xu_xy, &
                    id_dim_x_xz, id_dim_xu_xz, id_dim_x_yz, id_dim_xu_yz, &
                    id_dim_x_3d, id_dim_xu_3d, id_dim_y_xy, id_dim_yv_xy, &
                    id_dim_y_xz, id_dim_yv_xz, id_dim_y_yz, id_dim_yv_yz, &
                    id_dim_y_3d, id_dim_yv_3d, id_dim_zs_xy, id_dim_zs_xz, &
                    id_dim_zs_yz, id_dim_zs_3d, id_dim_zpc_3d, &
                    id_dim_zu_xy, id_dim_zu1_xy, &
                    id_dim_zu_xz, id_dim_zu_yz, id_dim_zu_3d, id_dim_zw_xy, &
                    id_dim_zw_xz, id_dim_zw_yz, id_dim_zw_3d, id_set_xy, &
                    id_set_xz, id_set_yz, id_set_3d, id_var_ind_x_yz, &
                    id_var_ind_y_xz, id_var_ind_z_xy, id_var_time_xy, &
                    id_var_time_xz, id_var_time_yz, id_var_time_3d, id_var_x_xy, &
                    id_var_xu_xy, id_var_x_xz, id_var_xu_xz, id_var_x_yz, &
                    id_var_xu_yz, id_var_x_3d, id_var_xu_3d, id_var_y_xy, &
                    id_var_yv_xy, id_var_y_xz, id_var_yv_xz, id_var_y_yz, &
                    id_var_yv_yz, id_var_y_3d, id_var_yv_3d, id_var_zs_xy, &
                    id_var_zs_xz, id_var_zs_yz, id_var_zs_3d, id_var_zpc_3d, &
                    id_var_zusi_xy, id_var_zusi_3d, id_var_zu_xy, id_var_zu1_xy, id_var_zu_xz, &
                    id_var_zu_yz, id_var_zu_3d, id_var_zwwi_xy, id_var_zwwi_3d, &
                    id_var_zw_xy, id_var_zw_xz, id_var_zw_yz, id_var_zw_3d

    INTEGER(iwp), DIMENSION(0:2,0:1) ::  id_var_eutm_3d, id_var_nutm_3d, &
                                         id_var_eutm_xy, id_var_nutm_xy, &
                                         id_var_eutm_xz, id_var_nutm_xz, &
                                         id_var_eutm_yz, id_var_nutm_yz

    INTEGER(iwp), DIMENSION(0:2,0:1) ::  id_var_lat_3d, id_var_lon_3d, &
                                         id_var_lat_xy, id_var_lon_xy, &
                                         id_var_lat_xz, id_var_lon_xz, &
                                         id_var_lat_yz, id_var_lon_yz

    INTEGER ::  netcdf_data_format = 2  !< NetCDF3 64bit offset format
    INTEGER ::  netcdf_deflate = 0      !< NetCDF compression, default: no
                                        !< compression

    INTEGER(iwp)                 ::  dofl_time_count
    INTEGER(iwp), DIMENSION(10)  ::  id_var_dospx, id_var_dospy
    INTEGER(iwp), DIMENSION(20)  ::  id_var_agt
!    INTEGER(iwp), DIMENSION(20)  ::  id_var_prt
    INTEGER(iwp), DIMENSION(11)  ::  nc_precision
    INTEGER(iwp), DIMENSION(dopr_norm_num) ::  id_var_norm_dopr

    INTEGER(iwp), DIMENSION(fl_max) ::  id_dim_x_fl, id_dim_y_fl, id_dim_z_fl
    INTEGER(iwp), DIMENSION(fl_max) ::  id_var_x_fl, id_var_y_fl, id_var_z_fl

    CHARACTER (LEN=20), DIMENSION(fl_max*var_fl_max) :: dofl_label
    CHARACTER (LEN=20), DIMENSION(fl_max*var_fl_max) :: dofl_unit
    CHARACTER (LEN=20), DIMENSION(fl_max) :: dofl_dim_label_x
    CHARACTER (LEN=20), DIMENSION(fl_max) :: dofl_dim_label_y
    CHARACTER (LEN=20), DIMENSION(fl_max) :: dofl_dim_label_z

    INTEGER(iwp), DIMENSION(fl_max*var_fl_max) :: id_var_dofl

    INTEGER(iwp), DIMENSION(dopts_num,0:10) ::  id_var_dopts
    INTEGER(iwp), DIMENSION(0:1,500)        ::  id_var_do2d, id_var_do3d
    INTEGER(iwp), DIMENSION(100,0:99)       ::  id_dim_z_pr, id_var_dopr, &
                                                id_var_z_pr
    INTEGER(iwp), DIMENSION(dots_max,0:99)  ::  id_var_dots

!
!-- Masked output
    CHARACTER (LEN=7), DIMENSION(max_masks,0:1,100) ::  domask_unit

    LOGICAL ::  output_for_t0 = .FALSE.

    INTEGER(iwp), DIMENSION(1:max_masks,0:1) ::  id_dim_time_mask, id_dim_x_mask, &
                   id_dim_xu_mask, id_dim_y_mask, id_dim_yv_mask, id_dim_zs_mask, &
                   id_dim_zu_mask, id_dim_zw_mask, &
                   id_set_mask, &
                   id_var_time_mask, id_var_x_mask, id_var_xu_mask, &
                   id_var_y_mask, id_var_yv_mask, id_var_zs_mask, &
                   id_var_zu_mask, id_var_zw_mask, &
                   id_var_zusi_mask, id_var_zwwi_mask

    INTEGER(iwp), DIMENSION(0:2,1:max_masks,0:1) ::  id_var_eutm_mask, &
                                                     id_var_nutm_mask

    INTEGER(iwp), DIMENSION(0:2,1:max_masks,0:1) ::  id_var_lat_mask, &
                                                     id_var_lon_mask

    INTEGER(iwp), DIMENSION(1:max_masks,0:1,100) ::  id_var_domask

    REAL(wp) ::  fill_value = -9999.0_wp    !< value for the _FillValue attribute


    PUBLIC  dofl_dim_label_x, dofl_dim_label_y, dofl_dim_label_z, dofl_label,  &
            dofl_time_count, dofl_unit, domask_unit, dopr_unit, dopts_num,     &
            dots_label, dots_max, dots_num, dots_rad, dots_soil, dots_unit,    &
            do2d_unit, do3d_unit, fill_value, id_set_agt, id_set_fl,           &
            id_set_mask, id_set_pr, id_set_prt, id_set_pts, id_set_sp,         &
            id_set_ts, id_set_xy, id_set_xz, id_set_yz, id_set_3d, id_var_agt, &
            id_var_domask, id_var_dofl, id_var_dopr, id_var_dopts,             &
            id_var_dospx, id_var_dospy, id_var_dots, id_var_do2d, id_var_do3d, &
            id_var_norm_dopr, id_var_time_agt, id_var_time_fl,                 &
            id_var_time_mask, id_var_time_pr, id_var_rnoa_agt, id_var_time_pts,&
            id_var_time_sp, id_var_time_ts,                                    &
            id_var_time_xy, id_var_time_xz, id_var_time_yz, id_var_time_3d,    &
            id_var_x_fl, id_var_y_fl, id_var_z_fl,  nc_stat,                   &
            netcdf_data_format, netcdf_data_format_string, netcdf_deflate,     &
            netcdf_precision, output_for_t0, heatflux_output_unit,             &
            waterflux_output_unit, momentumflux_output_unit

    SAVE

    INTERFACE netcdf_create_dim
       MODULE PROCEDURE netcdf_create_dim
    END INTERFACE netcdf_create_dim

    INTERFACE netcdf_create_file
       MODULE PROCEDURE netcdf_create_file
    END INTERFACE netcdf_create_file

    INTERFACE netcdf_create_global_atts
       MODULE PROCEDURE netcdf_create_global_atts
    END INTERFACE netcdf_create_global_atts

    INTERFACE netcdf_create_var
       MODULE PROCEDURE netcdf_create_var
    END INTERFACE netcdf_create_var

    INTERFACE netcdf_create_att
       MODULE PROCEDURE netcdf_create_att_string
    END INTERFACE netcdf_create_att

    INTERFACE netcdf_define_header
       MODULE PROCEDURE netcdf_define_header
    END INTERFACE netcdf_define_header

    INTERFACE netcdf_handle_error
       MODULE PROCEDURE netcdf_handle_error
    END INTERFACE netcdf_handle_error

    INTERFACE netcdf_open_write_file
       MODULE PROCEDURE netcdf_open_write_file
    END INTERFACE netcdf_open_write_file

    PUBLIC netcdf_create_att, netcdf_create_dim, netcdf_create_file,           &
           netcdf_create_global_atts, netcdf_create_var, netcdf_define_header, &
           netcdf_handle_error, netcdf_open_write_file

 CONTAINS

 SUBROUTINE netcdf_define_header( callmode, extend, av )

#if defined( __netcdf )

    USE arrays_3d,                                                             &
        ONLY:  zu, zw

    USE biometeorology_mod,                                                    &
        ONLY:  bio_define_netcdf_grid

    USE chemistry_model_mod,                                                   &
        ONLY:  chem_define_netcdf_grid

    USE basic_constants_and_equations_mod,                                     &
        ONLY:  convert_utm_to_geographic,                                      &
               pi

    USE control_parameters,                                                    &
        ONLY:  agent_time_unlimited, air_chemistry, averaging_interval,        &
               averaging_interval_pr, data_output_pr, domask, dopr_n,          &
               dopr_time_count, dopts_time_count, dots_time_count,             &
               do2d, do2d_at_begin, do2d_xz_time_count, do3d, do3d_at_begin,   &
               do2d_yz_time_count, dt_data_output_av, dt_do2d_xy, dt_do2d_xz,  &
               dt_do2d_yz, dt_do3d, dt_write_agent_data, mask_size,            &
               do2d_xy_time_count, do3d_time_count, domask_time_count,         &
               end_time, indoor_model, land_surface,                           &
               mask_size_l, mask_i, mask_i_global, mask_j, mask_j_global,      &
               mask_k_global, mask_surface,                                    &
               message_string, ntdim_2d_xy, ntdim_2d_xz,                       &
               ntdim_2d_yz, ntdim_3d, nz_do3d, ocean_mode, plant_canopy,       &
               run_description_header, salsa, section, simulated_time,         &
               simulated_time_at_begin, skip_time_data_output_av,              &
               skip_time_do2d_xy, skip_time_do2d_xz, skip_time_do2d_yz,        &
               skip_time_do3d, topography, num_leg, num_var_fl,                &
               urban_surface

    USE diagnostic_output_quantities_mod,                                      &
        ONLY:  doq_define_netcdf_grid

    USE grid_variables,                                                        &
        ONLY:  dx, dy, zu_s_inner, zw_w_inner

    USE gust_mod,                                                              &
        ONLY: gust_define_netcdf_grid, gust_module_enabled

    USE indices,                                                               &
        ONLY:  nx, nxl, nxr, ny, nys, nyn, nz ,nzb, nzt

    USE kinds

    USE indoor_model_mod,                                                      &
        ONLY: im_define_netcdf_grid

    USE land_surface_model_mod,                                                &
        ONLY: lsm_define_netcdf_grid, nzb_soil, nzt_soil, nzs, zs

    USE ocean_mod,                                                             &
        ONLY:  ocean_define_netcdf_grid

    USE pegrid

    USE particle_attributes,                                                   &
        ONLY:  number_of_particle_groups

    USE plant_canopy_model_mod,                                                &
        ONLY:  pch_index, pcm_define_netcdf_grid

    USE profil_parameter,                                                      &
        ONLY:  crmax, cross_profiles, dopr_index, profile_columns, profile_rows

    USE radiation_model_mod,                                                   &
        ONLY: radiation, radiation_define_netcdf_grid

    USE salsa_mod,                                                             &
        ONLY:  salsa_define_netcdf_grid

    USE spectra_mod,                                                           &
        ONLY:  averaging_interval_sp, comp_spectra_level, data_output_sp, dosp_time_count, spectra_direction

    USE statistics,                                                            &
        ONLY:  hom, statistic_regions

    USE turbulence_closure_mod,                                                &
        ONLY:  tcm_define_netcdf_grid

    USE urban_surface_mod,                                                     &
        ONLY:  usm_define_netcdf_grid

    USE user,                                                                  &
        ONLY:  user_module_enabled, user_define_netcdf_grid



    IMPLICIT NONE

    CHARACTER (LEN=3)              ::  suffix                !<
    CHARACTER (LEN=2), INTENT (IN) ::  callmode              !<
    CHARACTER (LEN=4)              ::  grid_x                !<
    CHARACTER (LEN=4)              ::  grid_y                !<
    CHARACTER (LEN=4)              ::  grid_z                !<
    CHARACTER (LEN=6)              ::  mode                  !<
    CHARACTER (LEN=10)             ::  precision             !<
    CHARACTER (LEN=10)             ::  var                   !<
    CHARACTER (LEN=20)             ::  netcdf_var_name       !<
    CHARACTER (LEN=varnamelength)  ::  trimvar               !< TRIM of output-variable string
    CHARACTER (LEN=80)             ::  time_average_text     !<
    CHARACTER (LEN=4000)           ::  char_cross_profiles   !<
    CHARACTER (LEN=4000)           ::  var_list              !<
    CHARACTER (LEN=4000)           ::  var_list_old          !<

    CHARACTER (LEN=100), DIMENSION(1:crmax) ::  cross_profiles_adj   !<
    CHARACTER (LEN=100), DIMENSION(1:crmax) ::  cross_profiles_char  !<

    INTEGER(iwp) ::  av                                      !<
    INTEGER(iwp) ::  cross_profiles_count                    !<
    INTEGER(iwp) ::  cross_profiles_maxi                     !<
    INTEGER(iwp) ::  delim                                   !<
    INTEGER(iwp) ::  delim_old                               !<
    INTEGER(iwp) ::  file_id                                 !<
    INTEGER(iwp) ::  i                                       !<
    INTEGER(iwp) ::  id_last                                 !<
    INTEGER(iwp) ::  id_x                                    !<
    INTEGER(iwp) ::  id_y                                    !<
    INTEGER(iwp) ::  id_z                                    !<
    INTEGER(iwp) ::  j                                       !<
    INTEGER(iwp) ::  k                                       !<
    INTEGER(iwp) ::  kk                                      !<
    INTEGER(iwp) ::  mid                                     !< masked output running index
    INTEGER(iwp) ::  ns                                      !<
    INTEGER(iwp) ::  ns_do                                   !< actual value of ns for soil model data
    INTEGER(iwp) ::  ns_old                                  !<
    INTEGER(iwp) ::  ntime_count                             !< number of time levels found in file
    INTEGER(iwp) ::  nz_old                                  !<
    INTEGER(iwp) ::  l                                       !<

    INTEGER(iwp), SAVE ::  oldmode                           !<

    INTEGER(iwp), DIMENSION(1) ::  id_dim_time_old           !<
    INTEGER(iwp), DIMENSION(1) ::  id_dim_x_yz_old           !<
    INTEGER(iwp), DIMENSION(1) ::  id_dim_y_xz_old           !<
    INTEGER(iwp), DIMENSION(1) ::  id_dim_zu_sp_old          !<
    INTEGER(iwp), DIMENSION(1) ::  id_dim_zu_xy_old          !<
    INTEGER(iwp), DIMENSION(1) ::  id_dim_zu_3d_old          !<
    INTEGER(iwp), DIMENSION(1) ::  id_dim_zu_mask_old        !<


    INTEGER(iwp), DIMENSION(1:crmax) ::  cross_profiles_numb !<

    LOGICAL ::  found                                        !<

    LOGICAL, INTENT (INOUT) ::  extend                       !<

    LOGICAL, SAVE ::  init_netcdf = .FALSE.                  !<

    REAL(wp) ::  cos_rot_angle                               !< cosine of rotation_angle
    REAL(wp) ::  eutm                                        !< easting (UTM)
    REAL(wp) ::  nutm                                        !< northing (UTM)
    REAL(wp) ::  shift_x                                     !< shift of x coordinate
    REAL(wp) ::  shift_y                                     !< shift of y coordinate
    REAL(wp) ::  sin_rot_angle                               !< sine of rotation_angle

    REAL(wp), DIMENSION(1) ::  last_time_coordinate          !< last time value in file

    REAL(wp), DIMENSION(:), ALLOCATABLE   ::  netcdf_data    !<
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  netcdf_data_2d !<
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  lat            !< latitude
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  lon            !< longitude


!
!-- Initializing actions
    IF ( .NOT. init_netcdf )  THEN
!
!--    Check and set accuracy for netCDF output. First set default value
       nc_precision = NF90_REAL4

       i = 1
       DO  WHILE ( netcdf_precision(i) /= ' ' )
          j = INDEX( netcdf_precision(i), '_' )
          IF ( j == 0 )  THEN
             WRITE ( message_string, * ) 'netcdf_precision must contain a ', &
                                         '"_"netcdf_precision(', i, ')="',   &
                                         TRIM( netcdf_precision(i) ),'"'
             CALL message( 'netcdf_define_header', 'PA0241', 2, 2, 0, 6, 0 )
          ENDIF

          var       = netcdf_precision(i)(1:j-1)
          precision = netcdf_precision(i)(j+1:)

          IF ( precision == 'NF90_REAL4' )  THEN
             j = NF90_REAL4
          ELSEIF ( precision == 'NF90_REAL8' )  THEN
             j = NF90_REAL8
          ELSE
             WRITE ( message_string, * ) 'illegal netcdf precision: ',  &
                                         'netcdf_precision(', i, ')="', &
                                         TRIM( netcdf_precision(i) ),'"'
             CALL message( 'netcdf_define_header', 'PA0242', 1, 2, 0, 6, 0 )
          ENDIF

          SELECT CASE ( var )
             CASE ( 'xy' )
                nc_precision(1) = j
             CASE ( 'xz' )
                nc_precision(2) = j
             CASE ( 'yz' )
                nc_precision(3) = j
             CASE ( '2d' )
                nc_precision(1:3) = j
             CASE ( '3d' )
                nc_precision(4) = j
             CASE ( 'pr' )
                nc_precision(5) = j
             CASE ( 'ts' )
                nc_precision(6) = j
             CASE ( 'sp' )
                nc_precision(7) = j
             CASE ( 'prt' )
                nc_precision(8) = j
             CASE ( 'masks' )
                nc_precision(11) = j
             CASE ( 'fl' )
                nc_precision(9) = j
             CASE ( 'all' )
                nc_precision    = j

             CASE DEFAULT
                WRITE ( message_string, * ) 'unknown variable in ' //          &
                                  'initialization_parameters ',                &
                                  'assignment: netcdf_precision(', i, ')="',   &
                                            TRIM( netcdf_precision(i) ),'"'
                CALL message( 'netcdf_define_header', 'PA0243', 1, 2, 0, 6, 0 )

          END SELECT

          i = i + 1
          IF ( i > 50 )  EXIT
       ENDDO

!
!--    Check for allowed parameter range
       IF ( netcdf_deflate < 0  .OR.  netcdf_deflate > 9 )  THEN
          WRITE ( message_string, '(A,I3,A)' ) 'netcdf_deflate out of ' //     &
                                      'range & given value: ', netcdf_deflate, &
                                      ', allowed range: 0-9'
          CALL message( 'netcdf_define_header', 'PA0355', 2, 2, 0, 6, 0 )
       ENDIF
!
!--    Data compression does not work with parallel NetCDF/HDF5
       IF ( netcdf_deflate > 0  .AND.  netcdf_data_format /= 3 )  THEN
          message_string = 'netcdf_deflate reset to 0'
          CALL message( 'netcdf_define_header', 'PA0356', 0, 1, 0, 6, 0 )

          netcdf_deflate = 0
       ENDIF

       init_netcdf = .TRUE.

    ENDIF

!
!-- Determine the mode to be processed
    IF ( extend )  THEN
       mode = callmode // '_ext'
    ELSE
       mode = callmode // '_new'
    ENDIF

!
!-- Select the mode to be processed. Possibilities are 3d, ma (mask), xy, xz,
!-- yz, pr (profiles), ps (particle timeseries), fl (flight data), ts
!-- (timeseries) or sp (spectra)
    SELECT CASE ( mode )

       CASE ( 'ma_new' )

!
!--       decompose actual parameter file_id (=formal parameter av) into
!--       mid and av
          file_id = av
          IF ( file_id <= 200+max_masks )  THEN
             mid = file_id - 200
             av = 0
          ELSE
             mid = file_id - (200+max_masks)
             av = 1
          ENDIF

!
!--       Define some global attributes of the dataset
          IF ( av == 0 )  THEN
             CALL netcdf_create_global_atts( id_set_mask(mid,av), 'podsmasked', TRIM( run_description_header ), 464 )
             time_average_text = ' '
          ELSE
             CALL netcdf_create_global_atts( id_set_mask(mid,av), 'podsmasked', TRIM( run_description_header ), 464 )
             WRITE ( time_average_text,'(F7.1,'' s avg'')' )  averaging_interval
             nc_stat = NF90_PUT_ATT( id_set_mask(mid,av), NF90_GLOBAL, 'time_avg',   &
                                     TRIM( time_average_text ) )
             CALL netcdf_handle_error( 'netcdf_define_header', 466 )
          ENDIF

!
!--       Define time coordinate for volume data (unlimited dimension)
          CALL netcdf_create_dim( id_set_mask(mid,av), 'time', NF90_UNLIMITED, &
                                  id_dim_time_mask(mid,av), 467 )
          CALL netcdf_create_var( id_set_mask(mid,av),                         &
                                  (/ id_dim_time_mask(mid,av) /), 'time',      &
                                  NF90_DOUBLE, id_var_time_mask(mid,av),       &
                                 'seconds', 'time', 468, 469, 000 )
          CALL netcdf_create_att( id_set_mask(mid,av), id_var_time_mask(mid,av), 'standard_name', 'time', 000)
          CALL netcdf_create_att( id_set_mask(mid,av), id_var_time_mask(mid,av), 'axis', 'T', 000)

!
!--       Define spatial dimensions and coordinates:
          IF ( mask_surface(mid) )  THEN
!
!--          In case of terrain-following output, the vertical dimensions are
!--          indices, not meters
             CALL netcdf_create_dim( id_set_mask(mid,av), 'ku_above_surf',     &
                                     mask_size(mid,3), id_dim_zu_mask(mid,av), &
                                     470 )
             CALL netcdf_create_var( id_set_mask(mid,av),                      &
                                     (/ id_dim_zu_mask(mid,av) /),             &
                                     'ku_above_surf',                          &
                                     NF90_DOUBLE, id_var_zu_mask(mid,av),      &
                                     '1', 'grid point above terrain',          &
                                     471, 472, 000 )
             CALL netcdf_create_att( id_set_mask(mid,av),                      &
                                     id_var_zu_mask(mid,av), 'axis', 'Z', 000)
                                     
             CALL netcdf_create_dim( id_set_mask(mid,av), 'kw_above_surf',     &
                                     mask_size(mid,3), id_dim_zw_mask(mid,av), &
                                     473 )
             CALL netcdf_create_var( id_set_mask(mid,av),                      &
                                     (/ id_dim_zw_mask(mid,av) /),             &
                                     'kw_above_surf',                          &
                                     NF90_DOUBLE, id_var_zw_mask(mid,av),      &
                                    '1', 'grid point above terrain',           &
                                    474, 475, 000 )
             CALL netcdf_create_att( id_set_mask(mid,av),                      &
                                     id_var_zw_mask(mid,av), 'axis', 'Z', 000)
          ELSE
!
!--          Define vertical coordinate grid (zu grid)
             CALL netcdf_create_dim( id_set_mask(mid,av), 'zu_3d',             &
                                     mask_size(mid,3), id_dim_zu_mask(mid,av), &
                                     470 )
             CALL netcdf_create_var( id_set_mask(mid,av),                      &
                                     (/ id_dim_zu_mask(mid,av) /), 'zu_3d',    &
                                     NF90_DOUBLE, id_var_zu_mask(mid,av),      &
                                     'meters', '', 471, 472, 000 )
             CALL netcdf_create_att( id_set_mask(mid,av),                      &
                                     id_var_zu_mask(mid,av), 'axis', 'Z', 000)
!
!--          Define vertical coordinate grid (zw grid)
             CALL netcdf_create_dim( id_set_mask(mid,av), 'zw_3d',             &
                                     mask_size(mid,3), id_dim_zw_mask(mid,av), &
                                     473 )
             CALL netcdf_create_var( id_set_mask(mid,av),                      &
                                     (/ id_dim_zw_mask(mid,av) /), 'zw_3d',    &
                                     NF90_DOUBLE, id_var_zw_mask(mid,av),      &
                                    'meters', '', 474, 475, 000 )
             CALL netcdf_create_att( id_set_mask(mid,av),                      &
                                     id_var_zw_mask(mid,av), 'axis', 'Z', 000)
          ENDIF
!
!--       Define x-axis (for scalar position)
          CALL netcdf_create_dim( id_set_mask(mid,av), 'x', mask_size(mid,1),  &
                                  id_dim_x_mask(mid,av), 476 )
          CALL netcdf_create_var( id_set_mask(mid,av),                         &
                                  (/ id_dim_x_mask(mid,av) /), 'x',            &
                                  NF90_DOUBLE, id_var_x_mask(mid,av),          &
                                  'meters', '', 477, 478, 000 )
          CALL netcdf_create_att( id_set_mask(mid,av), id_var_x_mask(mid,av),  &
                                  'axis', 'X', 000)
!
!--       Define x-axis (for u position)
          CALL netcdf_create_dim( id_set_mask(mid,av), 'xu', mask_size(mid,1), &
                                  id_dim_xu_mask(mid,av), 479 )
          CALL netcdf_create_var( id_set_mask(mid,av),                         &
                                  (/ id_dim_xu_mask(mid,av) /), 'xu',          &
                                  NF90_DOUBLE, id_var_xu_mask(mid,av),         &
                                  'meters', '', 480, 481, 000 )
          CALL netcdf_create_att( id_set_mask(mid,av), id_var_xu_mask(mid,av), & 
                                  'axis', 'X', 000)
!
!--       Define y-axis (for scalar position)
          CALL netcdf_create_dim( id_set_mask(mid,av), 'y', mask_size(mid,2),  &
                                  id_dim_y_mask(mid,av), 482 )
          CALL netcdf_create_var( id_set_mask(mid,av),                         &
                                  (/ id_dim_y_mask(mid,av) /), 'y',            &
                                  NF90_DOUBLE, id_var_y_mask(mid,av),          &
                                  'meters', '', 483, 484, 000 )
          CALL netcdf_create_att( id_set_mask(mid,av), id_var_y_mask(mid,av),  &
                                  'axis', 'Y', 000)
!
!--       Define y-axis (for v position)
          CALL netcdf_create_dim( id_set_mask(mid,av), 'yv', mask_size(mid,2), &
                                  id_dim_yv_mask(mid,av), 485 )
          CALL netcdf_create_var( id_set_mask(mid,av),                         &
                                  (/ id_dim_yv_mask(mid,av) /),                &
                                  'yv', NF90_DOUBLE, id_var_yv_mask(mid,av),   &
                                  'meters', '', 486, 487, 000 )
          CALL netcdf_create_att( id_set_mask(mid,av), id_var_yv_mask(mid,av), & 
                                  'axis', 'Y', 000)
!
!--       Define UTM and geographic coordinates
          CALL define_geo_coordinates( id_set_mask(mid,av),               &
                  (/ id_dim_x_mask(mid,av), id_dim_xu_mask(mid,av) /),    &
                  (/ id_dim_y_mask(mid,av), id_dim_yv_mask(mid,av) /),    &
                  id_var_eutm_mask(:,mid,av), id_var_nutm_mask(:,mid,av), &
                  id_var_lat_mask(:,mid,av), id_var_lon_mask(:,mid,av)    )
!
!--       Define coordinate-reference system
          CALL netcdf_create_crs( id_set_mask(mid,av), 000 )
!
!--       In case of non-flat topography define 2d-arrays containing the height
!--       information. Only for parallel netcdf output.
          IF ( TRIM( topography ) /= 'flat'  .AND.                             &
               netcdf_data_format > 4 )  THEN
!
!--          Define zusi = zu(nzb_s_inner)
             CALL netcdf_create_var( id_set_mask(mid,av),                      &
                                     (/ id_dim_x_mask(mid,av),                 &
                                        id_dim_y_mask(mid,av) /), 'zusi',      &
                                     NF90_DOUBLE, id_var_zusi_mask(mid,av),    &
                                     'meters', 'zu(nzb_s_inner)', 488, 489,    &
                                     490 )
!
!--          Define zwwi = zw(nzb_w_inner)
             CALL netcdf_create_var( id_set_mask(mid,av),                      &
                                     (/ id_dim_x_mask(mid,av),                 &
                                        id_dim_y_mask(mid,av) /), 'zwwi',      &
                                     NF90_DOUBLE, id_var_zwwi_mask(mid,av),    &
                                     'meters', 'zw(nzb_w_inner)', 491, 492,    &
                                     493 )
          ENDIF

          IF ( land_surface )  THEN
!
!--          Define vertical coordinate grid (zw grid)
             CALL netcdf_create_dim( id_set_mask(mid,av), 'zs_3d',             &
                                     mask_size(mid,3), id_dim_zs_mask(mid,av), &
                                     536 )
             CALL netcdf_create_var( id_set_mask(mid,av),                      &
                                     (/ id_dim_zs_mask(mid,av) /), 'zs_3d',    &
                                     NF90_DOUBLE, id_var_zs_mask(mid,av),      &
                                     'meters', '', 537, 555, 000 )
             CALL netcdf_create_att( id_set_mask(mid,av),                      &
                                     id_var_zs_mask(mid,av), 'axis', 'Z', 000)
          ENDIF

!
!--       Define the variables
          var_list = ';'
          i = 1

          DO WHILE ( domask(mid,av,i)(1:1) /= ' ' )

             trimvar = TRIM( domask(mid,av,i) )
!
!--          Check for the grid
             found = .FALSE.
             SELECT CASE ( trimvar )
!
!--             Most variables are defined on the scalar grid
                CASE ( 'e', 'nc', 'nr', 'p', 'pc', 'pr', 'prr',                &
                       'q', 'qc', 'ql', 'ql_c', 'ql_v', 'ql_vp', 'qr', 'qv',   &
                       's', 'theta', 'thetal', 'thetav' )

                   grid_x = 'x'
                   grid_y = 'y'
                   grid_z = 'zu'
!
!--             u grid
                CASE ( 'u' )

                   grid_x = 'xu'
                   grid_y = 'y'
                   grid_z = 'zu'
!
!--             v grid
                CASE ( 'v' )

                   grid_x = 'x'
                   grid_y = 'yv'
                   grid_z = 'zu'
!
!--             w grid
                CASE ( 'w' )

                   grid_x = 'x'
                   grid_y = 'y'
                   grid_z = 'zw'


                CASE DEFAULT
!
!--                Check for quantities defined in other modules
                   CALL tcm_define_netcdf_grid( trimvar, found, grid_x, grid_y, grid_z )

                   IF ( .NOT. found  .AND.  air_chemistry )  THEN
                      CALL chem_define_netcdf_grid( trimvar, found, grid_x, grid_y, grid_z )
                   ENDIF

                   IF ( .NOT. found )                                          &
                      CALL doq_define_netcdf_grid( trimvar, found, grid_x, grid_y, grid_z )

                   IF ( .NOT. found  .AND.  gust_module_enabled )  THEN
                      CALL gust_define_netcdf_grid( trimvar, found, grid_x, grid_y, grid_z )
                   ENDIF

                   IF ( .NOT. found  .AND. land_surface )  THEN
                      CALL lsm_define_netcdf_grid( trimvar, found, grid_x, grid_y, grid_z )
                   ENDIF

                   IF ( .NOT. found  .AND.  ocean_mode )  THEN
                      CALL ocean_define_netcdf_grid( trimvar, found, grid_x, grid_y, grid_z )
                   ENDIF

                   IF ( .NOT. found  .AND.  plant_canopy )  THEN
                      CALL pcm_define_netcdf_grid( trimvar, found, grid_x, grid_y, grid_z )
                   ENDIF

                   IF ( .NOT. found  .AND.  radiation )  THEN
                      CALL radiation_define_netcdf_grid( trimvar, found, grid_x, grid_y, grid_z )
                   ENDIF
!
!--                Check for SALSA quantities
                   IF ( .NOT. found  .AND.  salsa )  THEN
                      CALL salsa_define_netcdf_grid( trimvar, found, grid_x, grid_y, grid_z )
                   ENDIF

                   IF ( .NOT. found  .AND.  urban_surface )  THEN
                      CALL usm_define_netcdf_grid( trimvar, found, grid_x, grid_y, grid_z )
                   ENDIF
!
!--                Now check for user-defined quantities
                   IF ( .NOT. found  .AND.  user_module_enabled )  THEN
                      CALL user_define_netcdf_grid( trimvar, found, grid_x, grid_y, grid_z )
                   ENDIF

                   IF ( .NOT. found )  THEN
                      WRITE ( message_string, * ) 'no grid defined for variable ', TRIM( trimvar )
                      CALL message( 'define_netcdf_header', 'PA0244', 0, 1, 0, 6, 0 )
                   ENDIF

             END SELECT

!
!--          Select the respective dimension ids
             IF ( grid_x == 'x' )  THEN
                id_x = id_dim_x_mask(mid,av)
             ELSEIF ( grid_x == 'xu' )  THEN
                id_x = id_dim_xu_mask(mid,av)
             ENDIF

             IF ( grid_y == 'y' )  THEN
                id_y = id_dim_y_mask(mid,av)
             ELSEIF ( grid_y == 'yv' )  THEN
                id_y = id_dim_yv_mask(mid,av)
             ENDIF

             IF ( grid_z == 'zu' )  THEN
                id_z = id_dim_zu_mask(mid,av)
             ELSEIF ( grid_z == 'zw' )  THEN
                id_z = id_dim_zw_mask(mid,av)
             ELSEIF ( grid_z == "zs" )  THEN
                id_z = id_dim_zs_mask(mid,av)
             ENDIF

!
!--          Define the grid
             CALL netcdf_create_var( id_set_mask(mid,av), (/ id_x, id_y, id_z, &
                                     id_dim_time_mask(mid,av) /),              &
                                     domask(mid,av,i), nc_precision(11),       &
                                     id_var_domask(mid,av,i),                  &
                                     TRIM( domask_unit(mid,av,i) ),            &
                                     domask(mid,av,i), 494, 495, 496, .TRUE. )

             var_list = TRIM( var_list ) // TRIM( domask(mid,av,i) ) // ';'

             i = i + 1

          ENDDO

!
!--       No arrays to output
          IF ( i == 1 )  RETURN

!
!--       Write the list of variables as global attribute (this is used by
!--       restart runs and by combine_plot_fields)
          nc_stat = NF90_PUT_ATT( id_set_mask(mid,av), NF90_GLOBAL, &
                                  'VAR_LIST', var_list )
          CALL netcdf_handle_error( 'netcdf_define_header', 497 )

!
!--       Leave netCDF define mode
          nc_stat = NF90_ENDDEF( id_set_mask(mid,av) )
          CALL netcdf_handle_error( 'netcdf_define_header', 498 )

!
!--       Write data for x (shifted by +dx/2) and xu axis
          ALLOCATE( netcdf_data(mask_size(mid,1)) )

          netcdf_data = ( mask_i_global(mid,:mask_size(mid,1)) + 0.5_wp ) * dx

          nc_stat = NF90_PUT_VAR( id_set_mask(mid,av), id_var_x_mask(mid,av), &
                                  netcdf_data, start = (/ 1 /),               &
                                  count = (/ mask_size(mid,1) /) )
          CALL netcdf_handle_error( 'netcdf_define_header', 499 )

          netcdf_data = mask_i_global(mid,:mask_size(mid,1)) * dx

          nc_stat = NF90_PUT_VAR( id_set_mask(mid,av), id_var_xu_mask(mid,av),&
                                  netcdf_data, start = (/ 1 /),               &
                                  count = (/ mask_size(mid,1) /) )
          CALL netcdf_handle_error( 'netcdf_define_header', 500 )

          DEALLOCATE( netcdf_data )

!
!--       Write data for y (shifted by +dy/2) and yv axis
          ALLOCATE( netcdf_data(mask_size(mid,2)) )

          netcdf_data = ( mask_j_global(mid,:mask_size(mid,2)) + 0.5_wp ) * dy

          nc_stat = NF90_PUT_VAR( id_set_mask(mid,av), id_var_y_mask(mid,av), &
                                  netcdf_data, start = (/ 1 /),               &
                                  count = (/ mask_size(mid,2) /))
          CALL netcdf_handle_error( 'netcdf_define_header', 501 )

          netcdf_data = mask_j_global(mid,:mask_size(mid,2)) * dy

          nc_stat = NF90_PUT_VAR( id_set_mask(mid,av), id_var_yv_mask(mid,av), &
                                  netcdf_data, start = (/ 1 /),    &
                                  count = (/ mask_size(mid,2) /))
          CALL netcdf_handle_error( 'netcdf_define_header', 502 )

          DEALLOCATE( netcdf_data )

!
!--       Write UTM coordinates
          IF ( rotation_angle == 0.0_wp )  THEN
!
!--          1D in case of no rotation
             cos_rot_angle = COS( rotation_angle * pi / 180.0_wp )
!
!--          x coordinates
             ALLOCATE( netcdf_data(mask_size(mid,1)) )
             DO  k = 0, 2
!
!--             Scalar grid points
                IF ( k == 0 )  THEN
                   shift_x = 0.5
!
!--             u grid points
                ELSEIF ( k == 1 )  THEN
                   shift_x = 0.0
!
!--             v grid points
                ELSEIF ( k == 2 )  THEN
                   shift_x = 0.5
                ENDIF

                netcdf_data = init_model%origin_x + cos_rot_angle              &
                       * ( mask_i_global(mid,:mask_size(mid,1)) + shift_x ) * dx

                nc_stat = NF90_PUT_VAR( id_set_mask(mid,av), &
                                        id_var_eutm_mask(k,mid,av), &
                                        netcdf_data, start = (/ 1 /), &
                                        count = (/ mask_size(mid,1) /) )
                CALL netcdf_handle_error( 'netcdf_define_header', 555 )

             ENDDO
             DEALLOCATE( netcdf_data )
!
!--          y coordinates
             ALLOCATE( netcdf_data(mask_size(mid,2)) )
             DO  k = 0, 2
!
!--             Scalar grid points
                IF ( k == 0 )  THEN
                   shift_y = 0.5
!
!--             u grid points
                ELSEIF ( k == 1 )  THEN
                   shift_y = 0.5
!
!--             v grid points
                ELSEIF ( k == 2 )  THEN
                   shift_y = 0.0
                ENDIF

                netcdf_data = init_model%origin_y + cos_rot_angle              &
                       * ( mask_j_global(mid,:mask_size(mid,2)) + shift_y ) * dy

                nc_stat = NF90_PUT_VAR( id_set_mask(mid,av), &
                                        id_var_nutm_mask(k,mid,av), &
                                        netcdf_data, start = (/ 1 /), &
                                        count = (/ mask_size(mid,2) /) )
                CALL netcdf_handle_error( 'netcdf_define_header', 556 )

             ENDDO
             DEALLOCATE( netcdf_data )

          ELSE
!
!--          2D in case of rotation
             ALLOCATE( netcdf_data_2d(mask_size(mid,1),mask_size(mid,2)) )
             cos_rot_angle = COS( rotation_angle * pi / 180.0_wp )
             sin_rot_angle = SIN( rotation_angle * pi / 180.0_wp )

             DO  k = 0, 2
!
!--            Scalar grid points
               IF ( k == 0 )  THEN
                  shift_x = 0.5 ; shift_y = 0.5
!
!--            u grid points
               ELSEIF ( k == 1 )  THEN
                  shift_x = 0.0 ; shift_y = 0.5
!
!--            v grid points
               ELSEIF ( k == 2 )  THEN
                  shift_x = 0.5 ; shift_y = 0.0
               ENDIF

               DO  j = 1, mask_size(mid,2)
                  DO  i = 1, mask_size(mid,1)
                     netcdf_data_2d(i,j) = init_model%origin_x                        &
                           + cos_rot_angle * ( mask_i_global(mid,i) + shift_x ) * dx  &
                           + sin_rot_angle * ( mask_j_global(mid,j) + shift_y ) * dy
                  ENDDO
               ENDDO

               nc_stat = NF90_PUT_VAR( id_set_mask(mid,av), &
                                       id_var_eutm_mask(k,mid,av), &
                                       netcdf_data_2d, start = (/ 1, 1 /), &
                                       count = (/ mask_size(mid,1), &
                                                  mask_size(mid,2) /) )
               CALL netcdf_handle_error( 'netcdf_define_header', 555 )

               DO  j = 1, mask_size(mid,2)
                  DO  i = 1, mask_size(mid,1)
                     netcdf_data_2d(i,j) = init_model%origin_y                        &
                           - sin_rot_angle * ( mask_i_global(mid,i) + shift_x ) * dx  &
                           + cos_rot_angle * ( mask_j_global(mid,j) + shift_y ) * dy
                  ENDDO
               ENDDO

               nc_stat = NF90_PUT_VAR( id_set_mask(mid,av), &
                                       id_var_nutm_mask(k,mid,av), &
                                       netcdf_data_2d, start = (/ 1, 1 /), &
                                       count = (/ mask_size(mid,1), &
                                                  mask_size(mid,2) /) )
               CALL netcdf_handle_error( 'netcdf_define_header', 556 )

             ENDDO
             DEALLOCATE( netcdf_data_2d )
          ENDIF
!
!--       Write lon and lat data.
          ALLOCATE( lat(mask_size(mid,1),mask_size(mid,2)) )
          ALLOCATE( lon(mask_size(mid,1),mask_size(mid,2)) )
          cos_rot_angle = COS( rotation_angle * pi / 180.0_wp )
          sin_rot_angle = SIN( rotation_angle * pi / 180.0_wp )

          DO  k = 0, 2
!
!--          Scalar grid points
             IF ( k == 0 )  THEN
                shift_x = 0.5 ; shift_y = 0.5
!
!--          u grid points
             ELSEIF ( k == 1 )  THEN
                shift_x = 0.0 ; shift_y = 0.5
!
!--          v grid points
             ELSEIF ( k == 2 )  THEN
                shift_x = 0.5 ; shift_y = 0.0
             ENDIF

             DO  j = 1, mask_size(mid,2)
                DO  i = 1, mask_size(mid,1)
                   eutm = init_model%origin_x                                      &
                        + cos_rot_angle * ( mask_i_global(mid,i) + shift_x ) * dx  &
                        + sin_rot_angle * ( mask_j_global(mid,j) + shift_y ) * dy
                   nutm = init_model%origin_y                                      &
                        - sin_rot_angle * ( mask_i_global(mid,i) + shift_x ) * dx  &
                        + cos_rot_angle * ( mask_j_global(mid,j) + shift_y ) * dy

                   CALL  convert_utm_to_geographic( crs_list,          &
                                                    eutm, nutm,        &
                                                    lon(i,j), lat(i,j) )
                ENDDO
             ENDDO

             nc_stat = NF90_PUT_VAR( id_set_mask(mid,av),           &
                                     id_var_lon_mask(k,mid,av),     &
                                     lon, start = (/ 1, 1 /),       &
                                     count = (/ mask_size(mid,1),   &
                                                mask_size(mid,2) /) )
             CALL netcdf_handle_error( 'netcdf_define_header', 556 )

             nc_stat = NF90_PUT_VAR( id_set_mask(mid,av),           &
                                     id_var_lat_mask(k,mid,av),     &
                                     lat, start = (/ 1, 1 /),       &
                                     count = (/ mask_size(mid,1),   &
                                                mask_size(mid,2) /) )
             CALL netcdf_handle_error( 'netcdf_define_header', 556 )
          ENDDO

          DEALLOCATE( lat )
          DEALLOCATE( lon )
!
!--       Write zu and zw data (vertical axes)
          ALLOCATE( netcdf_data(mask_size(mid,3)) )

          IF ( mask_surface(mid) )  THEN

             netcdf_data = mask_k_global(mid,:mask_size(mid,3))

             nc_stat = NF90_PUT_VAR( id_set_mask(mid,av), id_var_zu_mask(mid,av), &
                                     netcdf_data, start = (/ 1 /), &
                                     count = (/ mask_size(mid,3) /) )
             CALL netcdf_handle_error( 'netcdf_define_header', 503 )

             netcdf_data = mask_k_global(mid,:mask_size(mid,3))

             nc_stat = NF90_PUT_VAR( id_set_mask(mid,av), id_var_zw_mask(mid,av), &
                                     netcdf_data, start = (/ 1 /), &
                                     count = (/ mask_size(mid,3) /) )
             CALL netcdf_handle_error( 'netcdf_define_header', 504 )

          ELSE

             netcdf_data = zu( mask_k_global(mid,:mask_size(mid,3)) )

             nc_stat = NF90_PUT_VAR( id_set_mask(mid,av), id_var_zu_mask(mid,av), &
                                     netcdf_data, start = (/ 1 /), &
                                     count = (/ mask_size(mid,3) /) )
             CALL netcdf_handle_error( 'netcdf_define_header', 503 )

             netcdf_data = zw( mask_k_global(mid,:mask_size(mid,3)) )

             nc_stat = NF90_PUT_VAR( id_set_mask(mid,av), id_var_zw_mask(mid,av), &
                                     netcdf_data, start = (/ 1 /), &
                                     count = (/ mask_size(mid,3) /) )
             CALL netcdf_handle_error( 'netcdf_define_header', 504 )

          ENDIF

          DEALLOCATE( netcdf_data )

!
!--       In case of non-flat topography write height information
          IF ( TRIM( topography ) /= 'flat'  .AND.                             &
               netcdf_data_format > 4 )  THEN

             ALLOCATE( netcdf_data_2d(mask_size_l(mid,1),mask_size_l(mid,2)) )
             netcdf_data_2d = zu_s_inner( mask_i(mid,:mask_size_l(mid,1)),     &
                                          mask_j(mid,:mask_size_l(mid,2)) )

             nc_stat = NF90_PUT_VAR( id_set_mask(mid,av),                      &
                                     id_var_zusi_mask(mid,av),                 &
                                     netcdf_data_2d,                           &
                                     start = (/ 1, 1 /),                       &
                                     count = (/ mask_size_l(mid,1),            &
                                                mask_size_l(mid,2) /) )
             CALL netcdf_handle_error( 'netcdf_define_header', 505 )

             netcdf_data_2d = zw_w_inner( mask_i(mid,:mask_size_l(mid,1)),     &
                                          mask_j(mid,:mask_size_l(mid,2)) )

             nc_stat = NF90_PUT_VAR( id_set_mask(mid,av),                      &
                                     id_var_zwwi_mask(mid,av),                 &
                                     netcdf_data_2d,                           &
                                     start = (/ 1, 1 /),                       &
                                     count = (/ mask_size_l(mid,1),            &
                                                mask_size_l(mid,2) /) )
             CALL netcdf_handle_error( 'netcdf_define_header', 506 )

             DEALLOCATE( netcdf_data_2d )

          ENDIF
!
!--       soil is not in masked output for now - disable temporary this block
!          IF ( land_surface )  THEN
!
!--          Write zs data (vertical axes for soil model), use negative values
!--          to indicate soil depth
!             ALLOCATE( netcdf_data(mask_size(mid,3)) )
!
!             netcdf_data = zs( mask_k_global(mid,:mask_size(mid,3)) )
!
!             nc_stat = NF90_PUT_VAR( id_set_mask(mid,av), id_var_zs_mask(mid,av), &
!                                     netcdf_data, start = (/ 1 /), &
!                                     count = (/ mask_size(mid,3) /) )
!             CALL netcdf_handle_error( 'netcdf_define_header', 538 )
!
!             DEALLOCATE( netcdf_data )
!
!          ENDIF

!
!--       restore original parameter file_id (=formal parameter av) into av
          av = file_id


       CASE ( 'ma_ext' )

!
!--       decompose actual parameter file_id (=formal parameter av) into
!--       mid and av
          file_id = av
          IF ( file_id <= 200+max_masks )  THEN
             mid = file_id - 200
             av = 0
          ELSE
             mid = file_id - (200+max_masks)
             av = 1
          ENDIF

!
!--       Get the list of variables and compare with the actual run.
!--       First var_list_old has to be reset, since GET_ATT does not assign
!--       trailing blanks.
          var_list_old = ' '
          nc_stat = NF90_GET_ATT( id_set_mask(mid,av), NF90_GLOBAL, 'VAR_LIST',&
                                  var_list_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 507 )

          var_list = ';'
          i = 1
          DO WHILE ( domask(mid,av,i)(1:1) /= ' ' )
             var_list = TRIM(var_list) // TRIM( domask(mid,av,i) ) // ';'
             i = i + 1
          ENDDO

          IF ( av == 0 )  THEN
             var = '(mask)'
          ELSE
             var = '(mask_av)'
          ENDIF

          IF ( TRIM( var_list ) /= TRIM( var_list_old ) )  THEN
             WRITE ( message_string, * ) 'netCDF file for ', TRIM( var ),       &
                  ' data for mask', mid, ' from previous run found,',           &
                  '&but this file cannot be extended due to variable ',         &
                  'mismatch.&New file is created instead.'
             CALL message( 'define_netcdf_header', 'PA0335', 0, 1, 0, 6, 0 )
             extend = .FALSE.
             RETURN
          ENDIF

!
!--       Get and compare the number of vertical gridpoints
          nc_stat = NF90_INQ_VARID( id_set_mask(mid,av), 'zu_3d', &
                                    id_var_zu_mask(mid,av) )
          CALL netcdf_handle_error( 'netcdf_define_header', 508 )

          nc_stat = NF90_INQUIRE_VARIABLE( id_set_mask(mid,av),     &
                                           id_var_zu_mask(mid,av),  &
                                           dimids = id_dim_zu_mask_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 509 )
          id_dim_zu_mask(mid,av) = id_dim_zu_mask_old(1)

          nc_stat = NF90_INQUIRE_DIMENSION( id_set_mask(mid,av),               &
                                            id_dim_zu_mask(mid,av),            &
                                            len = nz_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 510 )

          IF ( mask_size(mid,3) /= nz_old )  THEN
             WRITE ( message_string, * ) 'netCDF file for ', TRIM( var ),      &
                  '&data for mask', mid, ' from previous run found,',          &
                  ' but this file cannot be extended due to mismatch in ',     &
                  ' number of vertical grid points.',                          &
                  '&New file is created instead.'
             CALL message( 'define_netcdf_header', 'PA0336', 0, 1, 0, 6, 0 )
             extend = .FALSE.
             RETURN
          ENDIF

!
!--       Get the id of the time coordinate (unlimited coordinate) and its
!--       last index on the file. The next time level is plmask..count+1.
!--       The current time must be larger than the last output time
!--       on the file.
          nc_stat = NF90_INQ_VARID( id_set_mask(mid,av), 'time',               &
                                    id_var_time_mask(mid,av) )
          CALL netcdf_handle_error( 'netcdf_define_header', 511 )

          nc_stat = NF90_INQUIRE_VARIABLE( id_set_mask(mid,av),                &
                                           id_var_time_mask(mid,av),           &
                                           dimids = id_dim_time_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 512 )
          id_dim_time_mask(mid,av) = id_dim_time_old(1)

          nc_stat = NF90_INQUIRE_DIMENSION( id_set_mask(mid,av),               &
                                            id_dim_time_mask(mid,av),          &
                                            len = domask_time_count(mid,av) )
          CALL netcdf_handle_error( 'netcdf_define_header', 513 )

          nc_stat = NF90_GET_VAR( id_set_mask(mid,av),                         &
                                  id_var_time_mask(mid,av),                    &
                                  last_time_coordinate,                        &
                                  start = (/ domask_time_count(mid,av) /),     &
                                  count = (/ 1 /) )
          CALL netcdf_handle_error( 'netcdf_define_header', 514 )

          IF ( last_time_coordinate(1) >= simulated_time )  THEN
             WRITE ( message_string, * ) 'netCDF file for ', TRIM( var ),      &
                  ' data for mask', mid, ' from previous run found,',          &
                  '&but this file cannot be extended because the current ',    &
                  'output time is less or equal than the last output time ',   &
                  'on this file.&New file is created instead.'
             CALL message( 'define_netcdf_header', 'PA0337', 0, 1, 0, 6, 0 )
             domask_time_count(mid,av) = 0
             extend = .FALSE.
             RETURN
          ENDIF

!
!--       Dataset seems to be extendable.
!--       Now get the variable ids.
          i = 1
          DO WHILE ( domask(mid,av,i)(1:1) /= ' ' )
             nc_stat = NF90_INQ_VARID( id_set_mask(mid,av), &
                                       TRIM( domask(mid,av,i) ), &
                                       id_var_domask(mid,av,i) )
             CALL netcdf_handle_error( 'netcdf_define_header', 515 )
             i = i + 1
          ENDDO

!
!--       Update the title attribute on file
!--       In order to avoid 'data mode' errors if updated attributes are larger
!--       than their original size, NF90_PUT_ATT is called in 'define mode'
!--       enclosed by NF90_REDEF and NF90_ENDDEF calls. This implies a possible
!--       performance loss due to data copying; an alternative strategy would be
!--       to ensure equal attribute size in a job chain. Maybe revise later.
          IF ( av == 0 )  THEN
             time_average_text = ' '
          ELSE
             WRITE (time_average_text, '('', '',F7.1,'' s average'')')         &
                                                            averaging_interval
          ENDIF
          nc_stat = NF90_REDEF( id_set_mask(mid,av) )
          CALL netcdf_handle_error( 'netcdf_define_header', 516 )
          nc_stat = NF90_PUT_ATT( id_set_mask(mid,av), NF90_GLOBAL, 'title',   &
                                  TRIM( run_description_header ) //            &
                                  TRIM( time_average_text ) )
          CALL netcdf_handle_error( 'netcdf_define_header', 517 )
          nc_stat = NF90_ENDDEF( id_set_mask(mid,av) )
          CALL netcdf_handle_error( 'netcdf_define_header', 518 )
          WRITE ( message_string, * ) 'netCDF file for ', TRIM( var ),         &
               ' data for mask', mid, ' from previous run found.',             &
               ' &This file will be extended.'
          CALL message( 'define_netcdf_header', 'PA0338', 0, 0, 0, 6, 0 )
!
!--       restore original parameter file_id (=formal parameter av) into av
          av = file_id


       CASE ( '3d_new' )

!
!--       Define some global attributes of the dataset
          IF ( av == 0 )  THEN
             CALL netcdf_create_global_atts( id_set_3d(av), '3d', TRIM( run_description_header ), 62 )
             time_average_text = ' '
          ELSE
             CALL netcdf_create_global_atts( id_set_3d(av), '3d_av', TRIM( run_description_header ), 62 )
             WRITE ( time_average_text,'(F7.1,'' s avg'')' )  averaging_interval
             nc_stat = NF90_PUT_ATT( id_set_3d(av), NF90_GLOBAL, 'time_avg',   &
                                     TRIM( time_average_text ) )
             CALL netcdf_handle_error( 'netcdf_define_header', 63 )
          ENDIF

!
!--       Define time coordinate for volume data.
!--       For parallel output the time dimensions has to be limited, otherwise
!--       the performance drops significantly.
          IF ( netcdf_data_format < 5 )  THEN
             CALL netcdf_create_dim( id_set_3d(av), 'time', NF90_UNLIMITED,    &
                                     id_dim_time_3d(av), 64 )
          ELSE
             CALL netcdf_create_dim( id_set_3d(av), 'time', ntdim_3d(av),      &
                                     id_dim_time_3d(av), 523 )
          ENDIF

          CALL netcdf_create_var( id_set_3d(av), (/ id_dim_time_3d(av) /),     &
                                  'time', NF90_DOUBLE, id_var_time_3d(av),     &
                                  'seconds', 'time', 65, 66, 00 )
          CALL netcdf_create_att( id_set_3d(av), id_var_time_3d(av), 'standard_name', 'time', 000)
          CALL netcdf_create_att( id_set_3d(av), id_var_time_3d(av), 'axis', 'T', 000)
!
!--       Define spatial dimensions and coordinates:
!--       Define vertical coordinate grid (zu grid)
          CALL netcdf_create_dim( id_set_3d(av), 'zu_3d', nz_do3d-nzb+1,       &
                                  id_dim_zu_3d(av), 67 )
          CALL netcdf_create_var( id_set_3d(av), (/ id_dim_zu_3d(av) /),       &
                                  'zu_3d', NF90_DOUBLE, id_var_zu_3d(av),      &
                                  'meters', '', 68, 69, 00 )
          CALL netcdf_create_att( id_set_3d(av), id_var_zu_3d(av), 'axis',     &
                                  'Z', 000)
!
!--       Define vertical coordinate grid (zw grid)
          CALL netcdf_create_dim( id_set_3d(av), 'zw_3d', nz_do3d-nzb+1,       &
                                  id_dim_zw_3d(av), 70 )
          CALL netcdf_create_var( id_set_3d(av), (/ id_dim_zw_3d(av) /),       &
                                  'zw_3d', NF90_DOUBLE, id_var_zw_3d(av),      &
                                  'meters', '', 71, 72, 00 )
          CALL netcdf_create_att( id_set_3d(av), id_var_zw_3d(av), 'axis',     &
                                  'Z', 000)
!
!--       Define x-axis (for scalar position)
          CALL netcdf_create_dim( id_set_3d(av), 'x', nx+1, id_dim_x_3d(av),   &
                                  73 )
          CALL netcdf_create_var( id_set_3d(av), (/ id_dim_x_3d(av) /), 'x',   &
                                  NF90_DOUBLE, id_var_x_3d(av), 'meters', '',  &
                                  74, 75, 00 )
          CALL netcdf_create_att( id_set_3d(av), id_var_x_3d(av), 'axis',     &
                                  'X', 000)
!
!--       Define x-axis (for u position)
          CALL netcdf_create_dim( id_set_3d(av), 'xu', nx+1, id_dim_xu_3d(av), &
                                  358 )
          CALL netcdf_create_var( id_set_3d(av), (/ id_dim_xu_3d(av) /), 'xu', &
                                  NF90_DOUBLE, id_var_xu_3d(av), 'meters', '', &
                                  359, 360, 000 )
          CALL netcdf_create_att( id_set_3d(av), id_var_xu_3d(av), 'axis',     &
                                  'X', 000)
!
!--       Define y-axis (for scalar position)
          CALL netcdf_create_dim( id_set_3d(av), 'y', ny+1, id_dim_y_3d(av),   &
                                  76 )
          CALL netcdf_create_var( id_set_3d(av), (/ id_dim_y_3d(av) /), 'y',   &
                                  NF90_DOUBLE, id_var_y_3d(av), 'meters', '',  &
                                  77, 78, 00 )
          CALL netcdf_create_att( id_set_3d(av), id_var_y_3d(av), 'axis',     &
                                  'Y', 000)
!
!--       Define y-axis (for v position)
          CALL netcdf_create_dim( id_set_3d(av), 'yv', ny+1, id_dim_yv_3d(av), &
                                  361 )
          CALL netcdf_create_var( id_set_3d(av), (/ id_dim_yv_3d(av) /), 'yv', &
                                  NF90_DOUBLE, id_var_yv_3d(av), 'meters', '', &
                                  362, 363, 000 )
          CALL netcdf_create_att( id_set_3d(av), id_var_yv_3d(av), 'axis',     &
                                  'Y', 000)
!
!--       Define UTM and geographic coordinates
          CALL define_geo_coordinates( id_set_3d(av),         &
                  (/ id_dim_x_3d(av), id_dim_xu_3d(av) /),    &
                  (/ id_dim_y_3d(av), id_dim_yv_3d(av) /),    &
                  id_var_eutm_3d(:,av), id_var_nutm_3d(:,av), &
                  id_var_lat_3d(:,av), id_var_lon_3d(:,av)    )
!
!--       Define coordinate-reference system
          CALL netcdf_create_crs( id_set_3d(av), 000 )
!
!--       In case of non-flat topography define 2d-arrays containing the height
!--       information. Only output 2d topography information in case of parallel
!--       output.
          IF ( TRIM( topography ) /= 'flat'  .AND.                             &
               netcdf_data_format > 4 )  THEN
!
!--          Define zusi = zu(nzb_s_inner)
             CALL netcdf_create_var( id_set_3d(av), (/ id_dim_x_3d(av),        &
                                     id_dim_y_3d(av) /), 'zusi', NF90_DOUBLE,  &
                                     id_var_zusi_3d(av), 'meters',             &
                                     'zu(nzb_s_inner)', 413, 414, 415 )
!
!--          Define zwwi = zw(nzb_w_inner)
             CALL netcdf_create_var( id_set_3d(av), (/ id_dim_x_3d(av),        &
                                     id_dim_y_3d(av) /), 'zwwi', NF90_DOUBLE,  &
                                     id_var_zwwi_3d(av), 'meters',             &
                                     'zw(nzb_w_inner)', 416, 417, 418 )

          ENDIF

          IF ( land_surface )  THEN
!
!--          Define vertical coordinate grid (zs grid)
             CALL netcdf_create_dim( id_set_3d(av), 'zs_3d',                   &
                                     nzt_soil-nzb_soil+1, id_dim_zs_3d(av), 70 )
             CALL netcdf_create_var( id_set_3d(av), (/ id_dim_zs_3d(av) /),    &
                                     'zs_3d', NF90_DOUBLE, id_var_zs_3d(av),   &
                                     'meters', '', 71, 72, 00 )
             CALL netcdf_create_att( id_set_3d(av), id_var_zs_3d(av), 'axis',  &
                                     'Z', 000)

          ENDIF

          IF ( plant_canopy )  THEN
!
!--          Define vertical coordinate grid (zpc grid)
             CALL netcdf_create_dim( id_set_3d(av), 'zpc_3d',                  &
                                     pch_index+1, id_dim_zpc_3d(av), 70 )
             !netcdf_create_dim(ncid, dim_name, ncdim_type, ncdim_id, error_no)
             CALL netcdf_create_var( id_set_3d(av), (/ id_dim_zpc_3d(av) /),   &
                                     'zpc_3d', NF90_DOUBLE, id_var_zpc_3d(av), &
                                     'meters', '', 71, 72, 00 )
             CALL netcdf_create_att( id_set_3d(av), id_var_zpc_3d(av), 'axis', &
                                     'Z', 000)

          ENDIF

!
!--       Define the variables
          var_list = ';'
          i = 1

          DO WHILE ( do3d(av,i)(1:1) /= ' ' )
!
!--          Temporary solution to account for data output within the new urban
!--          surface model (urban_surface_mod.f90), see also SELECT CASE ( trimvar )
             trimvar = TRIM( do3d(av,i) )
             IF ( urban_surface  .AND.  trimvar(1:4) == 'usm_' )  THEN
                trimvar = 'usm_output'
             ENDIF
!
!--          Check for the grid
             found = .FALSE.
             SELECT CASE ( trimvar )
!
!--             Most variables are defined on the scalar grid
                CASE ( 'e', 'nc', 'nr', 'p', 'pc', 'pr', 'prr',   &
                       'q', 'qc', 'ql', 'ql_c', 'ql_v', 'ql_vp', 'qr', 'qv',   &
                       's', 'theta', 'thetal', 'thetav' )

                   grid_x = 'x'
                   grid_y = 'y'
                   grid_z = 'zu'
!
!--             u grid
                CASE ( 'u' )

                   grid_x = 'xu'
                   grid_y = 'y'
                   grid_z = 'zu'
!
!--             v grid
                CASE ( 'v' )

                   grid_x = 'x'
                   grid_y = 'yv'
                   grid_z = 'zu'
!
!--             w grid
                CASE ( 'w' )

                   grid_x = 'x'
                   grid_y = 'y'
                   grid_z = 'zw'

!
!--             Block of urban surface model outputs
                CASE ( 'usm_output' )
                   CALL usm_define_netcdf_grid( do3d(av,i), found, &
                                                   grid_x, grid_y, grid_z )

                CASE DEFAULT

                   CALL tcm_define_netcdf_grid( do3d(av,i), found, &
                                                   grid_x, grid_y, grid_z )

!
!--                Check for land surface quantities
                   IF ( .NOT. found .AND. land_surface )  THEN
                      CALL lsm_define_netcdf_grid( do3d(av,i), found, grid_x,  &
                                                   grid_y, grid_z )
                   ENDIF
!
!--                Check for ocean quantities
                   IF ( .NOT. found  .AND.  ocean_mode )  THEN
                      CALL ocean_define_netcdf_grid( do3d(av,i), found,  &
                                                     grid_x, grid_y, grid_z )
                   ENDIF

!
!--                Check for plant canopy quantities
                   IF ( .NOT. found  .AND.  plant_canopy )  THEN
                      CALL pcm_define_netcdf_grid( do3d(av,i), found, grid_x,  &
                                                   grid_y, grid_z )
                   ENDIF

!
!--                Check for radiation quantities
                   IF ( .NOT. found  .AND.  radiation )  THEN
                      CALL radiation_define_netcdf_grid( do3d(av,i), found,    &
                                                         grid_x, grid_y,       &
                                                         grid_z )
                   ENDIF

!--                Check for gust module quantities
                   IF ( .NOT. found  .AND.  gust_module_enabled )  THEN
                      CALL gust_define_netcdf_grid( do3d(av,i), found, grid_x, &
                                                    grid_y, grid_z )
                   ENDIF
!
!--                Check for indoor model quantities
                   IF ( .NOT. found .AND. indoor_model ) THEN
                      CALL im_define_netcdf_grid( do3d(av,i), found,           &
                                                  grid_x, grid_y, grid_z )
                   ENDIF

!
!--                Check for biometeorology quantities
                   IF ( .NOT. found  .AND.  biometeorology )  THEN
                      CALL bio_define_netcdf_grid( do3d(av,i), found,          &
                                                   grid_x, grid_y, grid_z )
                   ENDIF

!
!--                Check for chemistry quantities
                   IF ( .NOT. found  .AND.  air_chemistry )  THEN
                      CALL chem_define_netcdf_grid( do3d(av,i), found,         &
                                                    grid_x, grid_y, grid_z )
                   ENDIF

!
!--                Check for SALSA quantities
                   IF ( .NOT. found  .AND.  salsa )  THEN
                      CALL salsa_define_netcdf_grid( do3d(av,i), found, grid_x,&
                                                     grid_y, grid_z )
                   ENDIF
!
!--                Check for user-defined quantities
                   IF ( .NOT. found  .AND.  user_module_enabled )  THEN
                      CALL user_define_netcdf_grid( do3d(av,i), found, grid_x, &
                                                    grid_y, grid_z )
                   ENDIF

                   IF ( .NOT. found )                                          &
                      CALL doq_define_netcdf_grid( do3d(av,i), found, grid_x,  &
                                                   grid_y, grid_z        )

                   IF ( .NOT. found )  THEN
                      WRITE ( message_string, * ) 'no grid defined for varia', &
                                                  'ble ', TRIM( do3d(av,i) )
                      CALL message( 'define_netcdf_header', 'PA0244', 0, 1, 0, &
                                    6, 0 )
                   ENDIF

             END SELECT

!
!--          Select the respective dimension ids
             IF ( grid_x == 'x' )  THEN
                id_x = id_dim_x_3d(av)
             ELSEIF ( grid_x == 'xu' )  THEN
                id_x = id_dim_xu_3d(av)
             ENDIF

             IF ( grid_y == 'y' )  THEN
                id_y = id_dim_y_3d(av)
             ELSEIF ( grid_y == 'yv' )  THEN
                id_y = id_dim_yv_3d(av)
             ENDIF

             IF ( grid_z == 'zu' )  THEN
                id_z = id_dim_zu_3d(av)
             ELSEIF ( grid_z == 'zw' )  THEN
                id_z = id_dim_zw_3d(av)
             ELSEIF ( grid_z == 'zs' )  THEN
                id_z = id_dim_zs_3d(av)
             ELSEIF ( grid_z == 'zpc' )  THEN
                id_z = id_dim_zpc_3d(av)
             ENDIF

!
!--          Define the grid
             CALL netcdf_create_var( id_set_3d(av),(/ id_x, id_y, id_z,        &
                                     id_dim_time_3d(av) /), do3d(av,i),        &
                                     nc_precision(4), id_var_do3d(av,i),       &
                                     TRIM( do3d_unit(av,i) ), do3d(av,i), 79,  &
                                     80, 357, .TRUE. )
#if defined( __netcdf4_parallel )
             IF ( netcdf_data_format > 4 )  THEN
!
!--             Set no fill for every variable to increase performance.
                nc_stat = NF90_DEF_VAR_FILL( id_set_3d(av),     &
                                             id_var_do3d(av,i), &
                                             NF90_NOFILL, 0 )
                CALL netcdf_handle_error( 'netcdf_define_header', 532 )
!
!--             Set collective io operations for parallel io
                nc_stat = NF90_VAR_PAR_ACCESS( id_set_3d(av),     &
                                               id_var_do3d(av,i), &
                                               NF90_COLLECTIVE )
                CALL netcdf_handle_error( 'netcdf_define_header', 445 )
             ENDIF
#endif
             var_list = TRIM( var_list ) // TRIM( do3d(av,i) ) // ';'

             i = i + 1

          ENDDO

!
!--       No arrays to output
          IF ( i == 1 )  RETURN

!
!--       Write the list of variables as global attribute (this is used by
!--       restart runs and by combine_plot_fields)
          nc_stat = NF90_PUT_ATT( id_set_3d(av), NF90_GLOBAL, 'VAR_LIST', &
                                  var_list )
          CALL netcdf_handle_error( 'netcdf_define_header', 81 )

!
!--       Set general no fill, otherwise the performance drops significantly for
!--       parallel output.
          nc_stat = NF90_SET_FILL( id_set_3d(av), NF90_NOFILL, oldmode )
          CALL netcdf_handle_error( 'netcdf_define_header', 528 )

!
!--       Leave netCDF define mode
          nc_stat = NF90_ENDDEF( id_set_3d(av) )
          CALL netcdf_handle_error( 'netcdf_define_header', 82 )

!
!--       These data are only written by PE0 for parallel output to increase
!--       the performance.
          IF ( myid == 0  .OR.  netcdf_data_format < 5 )  THEN
!
!--          Write data for x (shifted by +dx/2) and xu axis
             ALLOCATE( netcdf_data(0:nx) )

             DO  i = 0, nx
                netcdf_data(i) = ( i + 0.5 ) * dx
             ENDDO

             nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_x_3d(av),  &
                                     netcdf_data, start = (/ 1 /),    &
                                     count = (/ nx+1 /) )
             CALL netcdf_handle_error( 'netcdf_define_header', 83 )

             DO  i = 0, nx
                netcdf_data(i) = i * dx
             ENDDO

             nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_xu_3d(av), &
                                     netcdf_data, start = (/ 1 /),    &
                                     count = (/ nx+1 /) )
             CALL netcdf_handle_error( 'netcdf_define_header', 385 )

             DEALLOCATE( netcdf_data )

!
!--          Write data for y (shifted by +dy/2) and yv axis
             ALLOCATE( netcdf_data(0:ny) )

             DO  i = 0, ny
                netcdf_data(i) = ( i + 0.5_wp ) * dy
             ENDDO

             nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_y_3d(av),  &
                                     netcdf_data, start = (/ 1 /),    &
                                     count = (/ ny+1 /) )
             CALL netcdf_handle_error( 'netcdf_define_header', 84 )

             DO  i = 0, ny
                netcdf_data(i) = i * dy
             ENDDO

             nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_yv_3d(av), &
                                     netcdf_data, start = (/ 1 /),    &
                                     count = (/ ny+1 /))
             CALL netcdf_handle_error( 'netcdf_define_header', 387 )

             DEALLOCATE( netcdf_data )

!
!--          Write UTM coordinates
             IF ( rotation_angle == 0.0_wp )  THEN
!
!--             1D in case of no rotation
                cos_rot_angle = COS( rotation_angle * pi / 180.0_wp )
!
!--             x coordinates
                ALLOCATE( netcdf_data(0:nx) )
                DO  k = 0, 2
!
!--                Scalar grid points
                   IF ( k == 0 )  THEN
                      shift_x = 0.5
!
!--                u grid points
                   ELSEIF ( k == 1 )  THEN
                      shift_x = 0.0
!
!--                v grid points
                   ELSEIF ( k == 2 )  THEN
                      shift_x = 0.5
                   ENDIF

                   DO  i = 0, nx
                     netcdf_data(i) = init_model%origin_x                      &
                                    + cos_rot_angle * ( i + shift_x ) * dx
                   ENDDO

                   nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_eutm_3d(k,av),&
                                           netcdf_data, start = (/ 1 /),       &
                                           count = (/ nx+1 /) )
                   CALL netcdf_handle_error( 'netcdf_define_header', 555 )

                ENDDO
                DEALLOCATE( netcdf_data )
!
!--             y coordinates
                ALLOCATE( netcdf_data(0:ny) )
                DO  k = 0, 2
!
!--                Scalar grid points
                   IF ( k == 0 )  THEN
                      shift_y = 0.5
!
!--                u grid points
                   ELSEIF ( k == 1 )  THEN
                      shift_y = 0.5
!
!--                v grid points
                   ELSEIF ( k == 2 )  THEN
                      shift_y = 0.0
                   ENDIF

                   DO  j = 0, ny
                      netcdf_data(j) = init_model%origin_y                     &
                                     + cos_rot_angle * ( j + shift_y ) * dy
                   ENDDO

                   nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_nutm_3d(k,av),&
                                           netcdf_data, start = (/ 1 /),       &
                                           count = (/ ny+1 /) )
                   CALL netcdf_handle_error( 'netcdf_define_header', 556 )

                ENDDO
                DEALLOCATE( netcdf_data )

             ELSE
!
!--             2D in case of rotation
                ALLOCATE( netcdf_data_2d(0:nx,0:ny) )
                cos_rot_angle = COS( rotation_angle * pi / 180.0_wp )
                sin_rot_angle = SIN( rotation_angle * pi / 180.0_wp )

                DO  k = 0, 2
!
!--               Scalar grid points
                  IF ( k == 0 )  THEN
                     shift_x = 0.5 ; shift_y = 0.5
!
!--               u grid points
                  ELSEIF ( k == 1 )  THEN
                     shift_x = 0.0 ; shift_y = 0.5
!
!--               v grid points
                  ELSEIF ( k == 2 )  THEN
                     shift_x = 0.5 ; shift_y = 0.0
                  ENDIF

                  DO  j = 0, ny
                     DO  i = 0, nx
                        netcdf_data_2d(i,j) = init_model%origin_x                   &
                                            + cos_rot_angle * ( i + shift_x ) * dx  &
                                            + sin_rot_angle * ( j + shift_y ) * dy
                     ENDDO
                  ENDDO

                  nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_eutm_3d(k,av),  &
                                          netcdf_data_2d, start = (/ 1, 1 /),   &
                                          count = (/ nx+1, ny+1 /) )
                  CALL netcdf_handle_error( 'netcdf_define_header', 555 )

                  DO  j = 0, ny
                     DO  i = 0, nx
                        netcdf_data_2d(i,j) = init_model%origin_y                   &
                                            - sin_rot_angle * ( i + shift_x ) * dx  &
                                            + cos_rot_angle * ( j + shift_y ) * dy
                     ENDDO
                  ENDDO

                  nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_nutm_3d(k,av),  &
                                          netcdf_data_2d, start = (/ 1, 1 /),   &
                                          count = (/ nx+1, ny+1 /) )
                  CALL netcdf_handle_error( 'netcdf_define_header', 556 )

                ENDDO
                DEALLOCATE( netcdf_data_2d )
             ENDIF
!
!--          Write zu and zw data (vertical axes)
             nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_zu_3d(av),  &
                                     zu(nzb:nz_do3d), start = (/ 1 /), &
                                     count = (/ nz_do3d-nzb+1 /) )
             CALL netcdf_handle_error( 'netcdf_define_header', 85 )


             nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_zw_3d(av),  &
                                     zw(nzb:nz_do3d), start = (/ 1 /), &
                                     count = (/ nz_do3d-nzb+1 /) )
             CALL netcdf_handle_error( 'netcdf_define_header', 86 )

             IF ( land_surface )  THEN
!
!--             Write zs grid
                nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_zs_3d(av),  &
                                        - zs(nzb_soil:nzt_soil), start = (/ 1 /), &
                                        count = (/ nzt_soil-nzb_soil+1 /) )
                CALL netcdf_handle_error( 'netcdf_define_header', 86 )
             ENDIF

             IF ( plant_canopy )  THEN
!
!--             Write zpc grid
                nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_zpc_3d(av),  &
                                        zu(nzb:nzb+pch_index), start = (/ 1 /), &
                                        count = (/ pch_index+1 /) )
                CALL netcdf_handle_error( 'netcdf_define_header', 86 )
             ENDIF

          ENDIF
!
!--       Write lon and lat data. Only for parallel output.
          IF ( netcdf_data_format > 4 )  THEN

             ALLOCATE( lat(nxl:nxr,nys:nyn) )
             ALLOCATE( lon(nxl:nxr,nys:nyn) )
             cos_rot_angle = COS( rotation_angle * pi / 180.0_wp )
             sin_rot_angle = SIN( rotation_angle * pi / 180.0_wp )

             DO  k = 0, 2
!
!--             Scalar grid points
                IF ( k == 0 )  THEN
                   shift_x = 0.5 ; shift_y = 0.5
!
!--             u grid points
                ELSEIF ( k == 1 )  THEN
                   shift_x = 0.0 ; shift_y = 0.5
!
!--             v grid points
                ELSEIF ( k == 2 )  THEN
                   shift_x = 0.5 ; shift_y = 0.0
                ENDIF

                DO  j = nys, nyn
                   DO  i = nxl, nxr
                      eutm = init_model%origin_x                   &
                           + cos_rot_angle * ( i + shift_x ) * dx  &
                           + sin_rot_angle * ( j + shift_y ) * dy
                      nutm = init_model%origin_y                   &
                           - sin_rot_angle * ( i + shift_x ) * dx  &
                           + cos_rot_angle * ( j + shift_y ) * dy

                      CALL  convert_utm_to_geographic( crs_list,          &
                                                       eutm, nutm,        &
                                                       lon(i,j), lat(i,j) )
                   ENDDO
                ENDDO

                nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_lon_3d(k,av), &
                                     lon, start = (/ nxl+1, nys+1 /),       &
                                     count = (/ nxr-nxl+1, nyn-nys+1 /) )
                CALL netcdf_handle_error( 'netcdf_define_header', 556 )

                nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_lat_3d(k,av), &
                                     lat, start = (/ nxl+1, nys+1 /),       &
                                     count = (/ nxr-nxl+1, nyn-nys+1 /) )
                CALL netcdf_handle_error( 'netcdf_define_header', 556 )
             ENDDO

             DEALLOCATE( lat )
             DEALLOCATE( lon )

          ENDIF
!
!--       In case of non-flat topography write height information. Only for
!--       parallel netcdf output.
          IF ( TRIM( topography ) /= 'flat'  .AND.                             &
               netcdf_data_format > 4 )  THEN

!             IF ( nxr == nx  .AND.  nyn /= ny )  THEN
!                nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_zusi_3d(av),     &
!                                        zu_s_inner(nxl:nxr+1,nys:nyn),         &
!                                        start = (/ nxl+1, nys+1 /),            &
!                                        count = (/ nxr-nxl+2, nyn-nys+1 /) )
!             ELSEIF ( nxr /= nx  .AND.  nyn == ny )  THEN
!                nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_zusi_3d(av),     &
!                                        zu_s_inner(nxl:nxr,nys:nyn+1),         &
!                                        start = (/ nxl+1, nys+1 /),            &
!                                        count = (/ nxr-nxl+1, nyn-nys+2 /) )
!             ELSEIF ( nxr == nx  .AND.  nyn == ny )  THEN
!                nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_zusi_3d(av),     &
!                                        zu_s_inner(nxl:nxr+1,nys:nyn+1),       &
!                                        start = (/ nxl+1, nys+1 /),            &
!                                        count = (/ nxr-nxl+2, nyn-nys+2 /) )
!             ELSE
                nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_zusi_3d(av),     &
                                        zu_s_inner(nxl:nxr,nys:nyn),           &
                                        start = (/ nxl+1, nys+1 /),            &
                                        count = (/ nxr-nxl+1, nyn-nys+1 /) )
!             ENDIF
             CALL netcdf_handle_error( 'netcdf_define_header', 419 )

!             IF ( nxr == nx  .AND.  nyn /= ny )  THEN
!                nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_zwwi_3d(av),     &
!                                        zw_w_inner(nxl:nxr+1,nys:nyn),         &
!                                        start = (/ nxl+1, nys+1 /),            &
!                                        count = (/ nxr-nxl+2, nyn-nys+1 /) )
!             ELSEIF ( nxr /= nx  .AND.  nyn == ny )  THEN
!                nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_zwwi_3d(av),     &
!                                        zw_w_inner(nxl:nxr,nys:nyn+1),         &
!                                        start = (/ nxl+1, nys+1 /),            &
!                                        count = (/ nxr-nxl+1, nyn-nys+2 /) )
!             ELSEIF ( nxr == nx  .AND.  nyn == ny )  THEN
!                nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_zwwi_3d(av),     &
!                                        zw_w_inner(nxl:nxr+1,nys:nyn+1),       &
!                                        start = (/ nxl+1, nys+1 /),            &
!                                        count = (/ nxr-nxl+2, nyn-nys+2 /) )
!             ELSE
                nc_stat = NF90_PUT_VAR( id_set_3d(av), id_var_zwwi_3d(av),     &
                                        zw_w_inner(nxl:nxr,nys:nyn),           &
                                        start = (/ nxl+1, nys+1 /),            &
                                        count = (/ nxr-nxl+1, nyn-nys+1 /) )
!             ENDIF
             CALL netcdf_handle_error( 'netcdf_define_header', 420 )

          ENDIF

       CASE ( '3d_ext' )

!
!--       Get the list of variables and compare with the actual run.
!--       First var_list_old has to be reset, since GET_ATT does not assign
!--       trailing blanks.
          var_list_old = ' '
          nc_stat = NF90_GET_ATT( id_set_3d(av), NF90_GLOBAL, 'VAR_LIST', &
                                  var_list_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 87 )

          var_list = ';'
          i = 1
          DO WHILE ( do3d(av,i)(1:1) /= ' ' )
             var_list = TRIM(var_list) // TRIM( do3d(av,i) ) // ';'
             i = i + 1
          ENDDO

          IF ( av == 0 )  THEN
             var = '(3d)'
          ELSE
             var = '(3d_av)'
          ENDIF

          IF ( TRIM( var_list ) /= TRIM( var_list_old ) )  THEN
             message_string = 'netCDF file for volume data ' //             &
                              TRIM( var ) // ' from previous run found,' // &
                              '&but this file cannot be extended due to' // &
                              ' variable mismatch.' //                      &
                              '&New file is created instead.'
             CALL message( 'define_netcdf_header', 'PA0245', 0, 1, 0, 6, 0 )
             extend = .FALSE.
             RETURN
          ENDIF

!
!--       Get and compare the number of vertical gridpoints
          nc_stat = NF90_INQ_VARID( id_set_3d(av), 'zu_3d', id_var_zu_3d(av) )
          CALL netcdf_handle_error( 'netcdf_define_header', 88 )

          nc_stat = NF90_INQUIRE_VARIABLE( id_set_3d(av), id_var_zu_3d(av), &
                                           dimids = id_dim_zu_3d_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 89 )
          id_dim_zu_3d(av) = id_dim_zu_3d_old(1)

          nc_stat = NF90_INQUIRE_DIMENSION( id_set_3d(av), id_dim_zu_3d(av), &
                                            len = nz_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 90 )

          IF ( nz_do3d-nzb+1 /= nz_old )  THEN
              message_string = 'netCDF file for volume data ' //             &
                               TRIM( var ) // ' from previous run found,' // &
                               '&but this file cannot be extended due to' // &
                               ' mismatch in number of' //                   &
                               ' vertical grid points (nz_do3d).' //         &
                               '&New file is created instead.'
             CALL message( 'define_netcdf_header', 'PA0246', 0, 1, 0, 6, 0 )
             extend = .FALSE.
             RETURN
          ENDIF

!
!--       Get the id of the time coordinate (unlimited coordinate) and its
!--       last index on the file. The next time level is pl3d..count+1.
!--       The current time must be larger than the last output time
!--       on the file.
          nc_stat = NF90_INQ_VARID( id_set_3d(av), 'time', id_var_time_3d(av) )
          CALL netcdf_handle_error( 'netcdf_define_header', 91 )

          nc_stat = NF90_INQUIRE_VARIABLE( id_set_3d(av), id_var_time_3d(av), &
                                           dimids = id_dim_time_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 92 )

          id_dim_time_3d(av) = id_dim_time_old(1)

          nc_stat = NF90_INQUIRE_DIMENSION( id_set_3d(av), id_dim_time_3d(av), &
                                            len = ntime_count )
          CALL netcdf_handle_error( 'netcdf_define_header', 93 )

!
!--       For non-parallel output use the last output time level of the netcdf
!--       file because the time dimension is unlimited. In case of parallel
!--       output the variable ntime_count could get the value of 9*10E36 because
!--       the time dimension is limited.
          IF ( netcdf_data_format < 5 ) do3d_time_count(av) = ntime_count

          nc_stat = NF90_GET_VAR( id_set_3d(av), id_var_time_3d(av), &
                                  last_time_coordinate,              &
                                  start = (/ do3d_time_count(av) /), &
                                  count = (/ 1 /) )
          CALL netcdf_handle_error( 'netcdf_define_header', 94 )

          IF ( last_time_coordinate(1) >= simulated_time )  THEN
             message_string = 'netCDF file for volume data ' //             &
                              TRIM( var ) // ' from previous run found,' // &
                              '&but this file cannot be extended becaus' // &
                              'e the current output time' //                &
                              '&is less or equal than the last output t' // &
                              'ime on this file.' //                        &
                              '&New file is created instead.'
             CALL message( 'define_netcdf_header', 'PA0247', 0, 1, 0, 6, 0 )
             do3d_time_count(av) = 0
             extend = .FALSE.
             RETURN
          ENDIF

          IF ( netcdf_data_format > 4 )  THEN
!
!--          Check if the needed number of output time levels is increased
!--          compared to the number of time levels in the existing file.
             IF ( ntdim_3d(av) > ntime_count )  THEN
                message_string = 'netCDF file for volume data ' // &
                                 TRIM( var ) // ' from previous run found,' // &
                                 '&but this file cannot be extended becaus' // &
                                 'e the number of output time levels has b' // &
                                 'een increased compared to the previous s' // &
                                 'imulation.' //                               &
                                 '&New file is created instead.'
                CALL message( 'define_netcdf_header', 'PA0388', 0, 1, 0, 6, 0 )
                do3d_time_count(av) = 0
                extend = .FALSE.
!
!--             Recalculate the needed time levels for the new file.
                IF ( av == 0 )  THEN
                   ntdim_3d(0) = CEILING(                               &
                           ( end_time - MAX( skip_time_do3d,            &
                                             simulated_time_at_begin )  &
                           ) / dt_do3d )
                   IF ( do3d_at_begin )  ntdim_3d(0) = ntdim_3d(0) + 1
                ELSE
                   ntdim_3d(1) = CEILING(                               &
                           ( end_time - MAX( skip_time_data_output_av,  &
                                             simulated_time_at_begin )  &
                           ) / dt_data_output_av )
                ENDIF
                RETURN
             ENDIF
          ENDIF

!
!--       Dataset seems to be extendable.
!--       Now get the variable ids.
          i = 1
          DO WHILE ( do3d(av,i)(1:1) /= ' ' )
             nc_stat = NF90_INQ_VARID( id_set_3d(av), TRIM( do3d(av,i) ), &
                                       id_var_do3d(av,i) )
             CALL netcdf_handle_error( 'netcdf_define_header', 95 )
#if defined( __netcdf4_parallel )
!
!--          Set collective io operations for parallel io
             IF ( netcdf_data_format > 4 )  THEN
                nc_stat = NF90_VAR_PAR_ACCESS( id_set_3d(av),     &
                                               id_var_do3d(av,i), &
                                               NF90_COLLECTIVE )
                CALL netcdf_handle_error( 'netcdf_define_header', 453 )
             ENDIF
#endif
             i = i + 1
          ENDDO

!
!--       Update the title attribute on file
!--       In order to avoid 'data mode' errors if updated attributes are larger
!--       than their original size, NF90_PUT_ATT is called in 'define mode'
!--       enclosed by NF90_REDEF and NF90_ENDDEF calls. This implies a possible
!--       performance loss due to data copying; an alternative strategy would be
!--       to ensure equal attribute size. Maybe revise later.
          IF ( av == 0 )  THEN
             time_average_text = ' '
          ELSE
             WRITE (time_average_text, '('', '',F7.1,'' s average'')') &
                                                            averaging_interval
          ENDIF
          nc_stat = NF90_REDEF( id_set_3d(av) )
          CALL netcdf_handle_error( 'netcdf_define_header', 429 )
          nc_stat = NF90_PUT_ATT( id_set_3d(av), NF90_GLOBAL, 'title', &
                                  TRIM( run_description_header ) //    &
                                  TRIM( time_average_text ) )
          CALL netcdf_handle_error( 'netcdf_define_header', 96 )
          nc_stat = NF90_ENDDEF( id_set_3d(av) )
          CALL netcdf_handle_error( 'netcdf_define_header', 430 )
          message_string = 'netCDF file for volume data ' //             &
                           TRIM( var ) // ' from previous run found.' // &
                           '&This file will be extended.'
          CALL message( 'define_netcdf_header', 'PA0248', 0, 0, 0, 6, 0 )


       CASE ( 'ag_new' )

!
!--       Define some global attributes of the dataset
          nc_stat = NF90_PUT_ATT( id_set_agt, NF90_GLOBAL, 'title', &
                                  TRIM( run_description_header ) )
          CALL netcdf_handle_error( 'netcdf_define_header', 330 )
!
!--       Switch for unlimited time dimension
          IF ( agent_time_unlimited ) THEN
             CALL netcdf_create_dim( id_set_agt, 'time', NF90_UNLIMITED,       &
                                     id_dim_time_agt, 331 )
          ELSE
             CALL netcdf_create_dim( id_set_agt, 'time',                       &
                                     INT( ( MIN( multi_agent_system_end,       &
                                                 end_time ) -                  &
                                            multi_agent_system_start ) /       &
                                            dt_write_agent_data * 1.1 ),       &
                                     id_dim_time_agt, 331 )
          ENDIF

          CALL netcdf_create_var( id_set_agt, (/ id_dim_time_agt /), 'time',   &
                                  NF90_REAL4, id_var_time_agt, 'seconds', 'time',  &
                                  332, 333, 000 )
          CALL netcdf_create_att( id_set_agt, id_var_time_agt, 'standard_name', 'time', 000)
          CALL netcdf_create_att( id_set_agt, id_var_time_agt, 'axis', 'T', 000)

          CALL netcdf_create_dim( id_set_agt, 'agent_number',                  &
                                  dim_size_agtnum, id_dim_agtnum, 334 )

          CALL netcdf_create_var( id_set_agt, (/ id_dim_agtnum /),             &
                                  'agent_number', NF90_REAL4,                  &
                                  id_var_agtnum, 'agent number', '', 335,      &
                                  336, 000 )
!
!--       Define variable which contains the real number of agents in use
          CALL netcdf_create_var( id_set_agt, (/ id_dim_time_agt /),           &
                                  'real_num_of_agt', NF90_REAL4,               &
                                  id_var_rnoa_agt, 'agent number', '', 337,    &
                                  338, 000 )
          i = 1
          CALL netcdf_create_var( id_set_agt, (/ id_dim_agtnum,                &
                                  id_dim_time_agt /), agt_var_names(i),        &
                                  NF90_DOUBLE, id_var_agt(i),                  &
                                  TRIM( agt_var_units(i) ),                    &
                                  TRIM( agt_var_names(i) ), 339, 340, 341 )
!
!--       Define the variables
          DO  i = 2, 6
             CALL netcdf_create_var( id_set_agt, (/ id_dim_agtnum,             &
                                     id_dim_time_agt /), agt_var_names(i),     &
                                     NF90_REAL4, id_var_agt(i),                &
                                     TRIM( agt_var_units(i) ),                 &
                                     TRIM( agt_var_names(i) ), 339, 340, 341 )

          ENDDO
!
!--       Define vars for biometeorology
          IF ( biometeorology )  THEN
             CALL netcdf_create_var( id_set_agt, (/ id_dim_agtnum,             &
                                     id_dim_time_agt /), agt_var_names(7),     &
                                     nc_precision(8), id_var_agt(7),           &
                                     TRIM( agt_var_units(7) ),                 &
                                     TRIM( agt_var_names(7) ), 339, 340, 341 )                      
          ENDIF
!
!--       Leave netCDF define mode
          nc_stat = NF90_ENDDEF( id_set_agt )
          CALL netcdf_handle_error( 'netcdf_define_header', 342 )


!        CASE ( 'ag_ext' )
! !+?agent extend output for restart runs has to be adapted
!
! !
! !--       Get the id of the time coordinate (unlimited coordinate) and its
! !--       last index on the file. The next time level is prt..count+1.
! !--       The current time must be larger than the last output time
! !--       on the file.
!           nc_stat = NF90_INQ_VARID( id_set_agt, 'time', id_var_time_agt )
!           CALL netcdf_handle_error( 'netcdf_define_header', 343 )
!
!           nc_stat = NF90_INQUIRE_VARIABLE( id_set_agt, id_var_time_agt, &
!                                            dimids = id_dim_time_old )
!           CALL netcdf_handle_error( 'netcdf_define_header', 344 )
!           id_dim_time_agt = id_dim_time_old(1)
!
!           nc_stat = NF90_INQUIRE_DIMENSION( id_set_agt, id_dim_time_agt, &
!                                             len = agt_time_count )
!           CALL netcdf_handle_error( 'netcdf_define_header', 345 )
!
!           nc_stat = NF90_GET_VAR( id_set_agt, id_var_time_agt,  &
!                                   last_time_coordinate,         &
!                                   start = (/ agt_time_count /), &
!                                   count = (/ 1 /) )
!           CALL netcdf_handle_error( 'netcdf_define_header', 346 )
!
!           IF ( last_time_coordinate(1) >= simulated_time )  THEN
!              message_string = 'netCDF file for agents ' //                  &
!                               'from previous run found,' //                 &
!                               '&but this file cannot be extended becaus' // &
!                               'e the current output time' //                &
!                               '&is less or equal than the last output t' // &
!                               'ime on this file.' //                        &
!                               '&New file is created instead.'
!              CALL message( 'define_netcdf_header', 'PA0265', 0, 1, 0, 6, 0 )
!              agt_time_count = 0
!              extend = .FALSE.
!              RETURN
!           ENDIF
!
! !
! !--       Dataset seems to be extendable.
! !--       Now get the variable ids.
!           nc_stat = NF90_INQ_VARID( id_set_agt, 'real_num_of_agt', &
!                                     id_var_rnoa_agt )
!           CALL netcdf_handle_error( 'netcdf_define_header', 347 )
!
!           DO  i = 1, 17
!
!              nc_stat = NF90_INQ_VARID( id_set_agt, agt_var_names(i), &
!                                        id_var_prt(i) )
!              CALL netcdf_handle_error( 'netcdf_define_header', 348 )
!
!           ENDDO
!
!           message_string = 'netCDF file for particles ' // &
!                            'from previous run found.' //   &
!                            '&This file will be extended.'
!           CALL message( 'define_netcdf_header', 'PA0266', 0, 0, 0, 6, 0 )


       CASE ( 'xy_new' )

!
!--       Define some global attributes of the dataset
          IF ( av == 0 )  THEN
             CALL netcdf_create_global_atts( id_set_xy(av), 'xy', TRIM( run_description_header ), 97 )
             time_average_text = ' '
          ELSE
             CALL netcdf_create_global_atts( id_set_xy(av), 'xy_av', TRIM( run_description_header ), 97 )
             WRITE ( time_average_text,'(F7.1,'' s avg'')' )  averaging_interval
             nc_stat = NF90_PUT_ATT( id_set_xy(av), NF90_GLOBAL, 'time_avg',   &
                                     TRIM( time_average_text ) )
             CALL netcdf_handle_error( 'netcdf_define_header', 98 )
          ENDIF

!
!--       Define time coordinate for xy sections.
!--       For parallel output the time dimensions has to be limited, otherwise
!--       the performance drops significantly.
          IF ( netcdf_data_format < 5 )  THEN
             CALL netcdf_create_dim( id_set_xy(av), 'time', NF90_UNLIMITED,    &
                                     id_dim_time_xy(av), 99 )
          ELSE
             CALL netcdf_create_dim( id_set_xy(av), 'time', ntdim_2d_xy(av),   &
                                     id_dim_time_xy(av), 524 )
          ENDIF

          CALL netcdf_create_var( id_set_xy(av), (/ id_dim_time_xy(av) /),     &
                                  'time', NF90_DOUBLE, id_var_time_xy(av),     &
                                  'seconds', 'time', 100, 101, 000 )
          CALL netcdf_create_att( id_set_xy(av), id_var_time_xy(av), 'standard_name', 'time', 000)
          CALL netcdf_create_att( id_set_xy(av), id_var_time_xy(av), 'axis', 'T', 000)
!
!--       Define the spatial dimensions and coordinates for xy-sections.
!--       First, determine the number of horizontal sections to be written.
          IF ( section(1,1) == -9999 )  THEN
             RETURN
          ELSE
             ns = 1
             DO WHILE ( section(ns,1) /= -9999  .AND.  ns <= 100 )
                ns = ns + 1
             ENDDO
             ns = ns - 1
          ENDIF

!
!--       Define vertical coordinate grid (zu grid)
          CALL netcdf_create_dim( id_set_xy(av), 'zu_xy', ns,                  &
                                  id_dim_zu_xy(av), 102 )
          CALL netcdf_create_var( id_set_xy(av), (/ id_dim_zu_xy(av) /),       &
                                  'zu_xy', NF90_DOUBLE, id_var_zu_xy(av),      &
                                  'meters', '', 103, 104, 000 )
          CALL netcdf_create_att( id_set_xy(av), id_var_zu_xy(av), 'axis',     &
                                  'Z', 000)
!
!--       Define vertical coordinate grid (zw grid)
          CALL netcdf_create_dim( id_set_xy(av), 'zw_xy', ns,                  &
                                  id_dim_zw_xy(av), 105 )
          CALL netcdf_create_var( id_set_xy(av), (/ id_dim_zw_xy(av) /),       &
                                  'zw_xy', NF90_DOUBLE, id_var_zw_xy(av),      &
                                  'meters', '', 106, 107, 000 )
          CALL netcdf_create_att( id_set_xy(av), id_var_zw_xy(av), 'axis',     &
                                  'Z', 000)

          IF ( land_surface )  THEN

             ns_do = 1
             DO WHILE ( section(ns_do,1) /= -9999  .AND.  ns_do < nzs )
                ns_do = ns_do + 1
             ENDDO
!
!--          Define vertical coordinate grid (zs grid)
             CALL netcdf_create_dim( id_set_xy(av), 'zs_xy', ns_do,            &
                                     id_dim_zs_xy(av), 539 )
             CALL netcdf_create_var( id_set_xy(av), (/ id_dim_zs_xy(av) /),    &
                                     'zs_xy', NF90_DOUBLE, id_var_zs_xy(av),   &
                                     'meters', '', 540, 541, 000 )
             CALL netcdf_create_att( id_set_xy(av), id_var_zs_xy(av), 'axis',  &
                                     'Z', 000)

          ENDIF

!
!--       Define a pseudo vertical coordinate grid for the surface variables
!--       u* and t* to store their height level
          CALL netcdf_create_dim( id_set_xy(av), 'zu1_xy', 1,                  &
                                  id_dim_zu1_xy(av), 108 )
          CALL netcdf_create_var( id_set_xy(av), (/ id_dim_zu1_xy(av) /),      &
                                  'zu1_xy', NF90_DOUBLE, id_var_zu1_xy(av),    &
                                  'meters', '', 109, 110, 000 )
          CALL netcdf_create_att( id_set_xy(av), id_var_zu1_xy(av), 'axis',    &
                                  'Z', 000)
!
!--       Define a variable to store the layer indices of the horizontal cross
!--       sections, too
          CALL netcdf_create_var( id_set_xy(av), (/ id_dim_zu_xy(av) /),       &
                                  'ind_z_xy', NF90_DOUBLE,                     &
                                  id_var_ind_z_xy(av), 'gridpoints', '', 111,  &
                                  112, 000 )
!
!--       Define x-axis (for scalar position)
          CALL netcdf_create_dim( id_set_xy(av), 'x', nx+1, id_dim_x_xy(av),   &
                                  113 )
          CALL netcdf_create_var( id_set_xy(av), (/ id_dim_x_xy(av) /), 'x',   &
                                  NF90_DOUBLE, id_var_x_xy(av), 'meters', '',  &
                                  114, 115, 000 )
          CALL netcdf_create_att( id_set_xy(av), id_var_x_xy(av), 'axis',      &
                                  'X', 000)
!
!--       Define x-axis (for u position)
          CALL netcdf_create_dim( id_set_xy(av), 'xu', nx+1,                   &
                                  id_dim_xu_xy(av), 388 )
          CALL netcdf_create_var( id_set_xy(av), (/ id_dim_xu_xy(av) /), 'xu', &
                                  NF90_DOUBLE, id_var_xu_xy(av), 'meters', '', &
                                  389, 390, 000 )
          CALL netcdf_create_att( id_set_xy(av), id_var_xu_xy(av), 'axis',     &
                                  'X', 000)
!
!--       Define y-axis (for scalar position)
          CALL netcdf_create_dim( id_set_xy(av), 'y', ny+1, id_dim_y_xy(av),   &
                                  116 )
          CALL netcdf_create_var( id_set_xy(av), (/ id_dim_y_xy(av) /), 'y',   &
                                  NF90_DOUBLE, id_var_y_xy(av), 'meters', '',  &
                                  117, 118, 000 )
          CALL netcdf_create_att( id_set_xy(av), id_var_y_xy(av), 'axis',      &
                                  'Y', 000)
!
!--       Define y-axis (for scalar position)
          CALL netcdf_create_dim( id_set_xy(av), 'yv', ny+1,                   &
                                  id_dim_yv_xy(av), 364 )
          CALL netcdf_create_var( id_set_xy(av), (/ id_dim_yv_xy(av) /), 'yv', &
                                  NF90_DOUBLE, id_var_yv_xy(av), 'meters', '', &
                                  365, 366, 000 )
          CALL netcdf_create_att( id_set_xy(av), id_var_yv_xy(av), 'axis',     &
                                  'Y', 000)
!
!--       Define UTM and geographic coordinates
          CALL define_geo_coordinates( id_set_xy(av),         &
                  (/ id_dim_x_xy(av), id_dim_xu_xy(av) /),    &
                  (/ id_dim_y_xy(av), id_dim_yv_xy(av) /),    &
                  id_var_eutm_xy(:,av), id_var_nutm_xy(:,av), &
                  id_var_lat_xy(:,av), id_var_lon_xy(:,av)    )
!
!--       Define coordinate-reference system
          CALL netcdf_create_crs( id_set_xy(av), 000 )
!
!--       In case of non-flat topography define 2d-arrays containing the height
!--       information. Only for parallel netcdf output.
          IF ( TRIM( topography ) /= 'flat'  .AND.                             &
               netcdf_data_format > 4  )  THEN
!
!--          Define zusi = zu(nzb_s_inner)
             CALL netcdf_create_var( id_set_xy(av), (/ id_dim_x_xy(av),        &
                                     id_dim_y_xy(av) /), 'zusi', NF90_DOUBLE,  &
                                     id_var_zusi_xy(av), 'meters',             &
                                     'zu(nzb_s_inner)', 421, 422, 423 )
!
!--          Define zwwi = zw(nzb_w_inner)
             CALL netcdf_create_var( id_set_xy(av), (/ id_dim_x_xy(av),        &
                                     id_dim_y_xy(av) /), 'zwwi', NF90_DOUBLE,  &
                                     id_var_zwwi_xy(av), 'meters',             &
                                     'zw(nzb_w_inner)', 424, 425, 426 )

          ENDIF

!
!--       Define the variables
          var_list = ';'
          i = 1

          DO WHILE ( do2d(av,i)(1:1) /= ' ' )

             IF ( INDEX( do2d(av,i), 'xy' ) /= 0 )  THEN
!
!--             If there is a star in the variable name (u* or t*), it is a
!--             surface variable. Define it with id_dim_zu1_xy.
                IF ( INDEX( do2d(av,i), '*' ) /= 0 )  THEN

                   CALL netcdf_create_var( id_set_xy(av), (/ id_dim_x_xy(av),  &
                                           id_dim_y_xy(av), id_dim_zu1_xy(av), &
                                           id_dim_time_xy(av) /), do2d(av,i),  &
                                           nc_precision(1), id_var_do2d(av,i), &
                                           TRIM( do2d_unit(av,i) ),            &
                                           do2d(av,i), 119, 120, 354, .TRUE. )

                ELSE

!
!--                Check for the grid
                   found = .FALSE.
                   SELECT CASE ( do2d(av,i) )
!
!--                   Most variables are defined on the zu grid
                      CASE ( 'e_xy', 'nc_xy', 'nr_xy', 'p_xy',                 &
                             'pc_xy', 'pr_xy', 'prr_xy', 'q_xy',               &
                             'qc_xy', 'ql_xy', 'ql_c_xy', 'ql_v_xy',           &
                             'ql_vp_xy', 'qr_xy', 'qv_xy',                     &
                             's_xy',                                           &
                             'theta_xy', 'thetal_xy', 'thetav_xy' )

                         grid_x = 'x'
                         grid_y = 'y'
                         grid_z = 'zu'
!
!--                   u grid
                      CASE ( 'u_xy' )

                         grid_x = 'xu'
                         grid_y = 'y'
                         grid_z = 'zu'
!
!--                   v grid
                      CASE ( 'v_xy' )

                         grid_x = 'x'
                         grid_y = 'yv'
                         grid_z = 'zu'
!
!--                   w grid
                      CASE ( 'w_xy' )

                         grid_x = 'x'
                         grid_y = 'y'
                         grid_z = 'zw'


                      CASE DEFAULT
!
!--                      Check for land surface quantities
                         IF ( land_surface )  THEN
                            CALL lsm_define_netcdf_grid( do2d(av,i), found,    &
                                                   grid_x, grid_y, grid_z )
                         ENDIF

                         IF ( .NOT. found )  THEN
                            CALL tcm_define_netcdf_grid( do2d(av,i), found,    &
                                                         grid_x, grid_y,       &
                                                         grid_z )
                         ENDIF

!
!--                      Check for ocean quantities
                         IF ( .NOT. found  .AND.  ocean_mode )  THEN
                            CALL ocean_define_netcdf_grid( do2d(av,i), found,  &
                                                           grid_x, grid_y,     &
                                                           grid_z )
                         ENDIF
!
!--                      Check for radiation quantities
                         IF ( .NOT. found  .AND.  radiation )  THEN
                            CALL radiation_define_netcdf_grid( do2d(av,i),     &
                                                         found, grid_x, grid_y,&
                                                         grid_z )
                         ENDIF

!
!--                      Check for SALSA quantities
                         IF ( .NOT. found  .AND.  salsa )  THEN
                            CALL salsa_define_netcdf_grid( do2d(av,i), found,  &
                                                           grid_x, grid_y,     &
                                                           grid_z )
                         ENDIF

!
!--                      Check for gust module quantities
                         IF ( .NOT. found  .AND.  gust_module_enabled )  THEN
                            CALL gust_define_netcdf_grid( do2d(av,i), found,   &
                                                          grid_x, grid_y,      &
                                                          grid_z )
                         ENDIF
!
!--                      Check for biometeorology quantities
                         IF ( .NOT. found  .AND.  biometeorology )  THEN
                            CALL bio_define_netcdf_grid( do2d( av, i), found,  &
                                                         grid_x, grid_y,       &
                                                         grid_z )
                         ENDIF
!
!--                      Check for chemistry quantities
                         IF ( .NOT. found  .AND.  air_chemistry )  THEN
                            CALL chem_define_netcdf_grid( do2d(av,i), found,   &
                                                          grid_x, grid_y,      &
                                                          grid_z )
                         ENDIF

                         IF ( .NOT. found )                                    &
                            CALL doq_define_netcdf_grid(                       &
                                                    do2d(av,i), found, grid_x, &
                                                    grid_y, grid_z              )
!
!--                      Check for user-defined quantities
                         IF ( .NOT. found  .AND.  user_module_enabled )  THEN
                            CALL user_define_netcdf_grid( do2d(av,i), found,   &
                                                          grid_x, grid_y,      &
                                                          grid_z )
                         ENDIF

                         IF ( .NOT. found )  THEN
                            WRITE ( message_string, * ) 'no grid defined for', &
                                                ' variable ', TRIM( do2d(av,i) )
                            CALL message( 'define_netcdf_header', 'PA0244',    &
                                          0, 1, 0, 6, 0 )
                         ENDIF

                   END SELECT

!
!--                Select the respective dimension ids
                   IF ( grid_x == 'x' )  THEN
                      id_x = id_dim_x_xy(av)
                   ELSEIF ( grid_x == 'xu' )  THEN
                      id_x = id_dim_xu_xy(av)
                   ENDIF

                   IF ( grid_y == 'y' )  THEN
                      id_y = id_dim_y_xy(av)
                   ELSEIF ( grid_y == 'yv' )  THEN
                      id_y = id_dim_yv_xy(av)
                   ENDIF

                   IF ( grid_z == 'zu' )  THEN
                      id_z = id_dim_zu_xy(av)
                   ELSEIF ( grid_z == 'zw' )  THEN
                      id_z = id_dim_zw_xy(av)
                   ELSEIF ( grid_z == 'zs' )  THEN
                      id_z = id_dim_zs_xy(av)
                   ELSEIF ( grid_z == 'zu1' )  THEN
                      id_z = id_dim_zu1_xy(av)
                   ENDIF

!
!--                Define the grid
                   CALL netcdf_create_var( id_set_xy(av), (/ id_x, id_y, id_z, &
                                           id_dim_time_xy(av) /), do2d(av,i),  &
                                           nc_precision(1), id_var_do2d(av,i), &
                                           TRIM( do2d_unit(av,i) ),            &
                                           do2d(av,i), 119, 120, 354, .TRUE. )

                ENDIF

#if defined( __netcdf4_parallel )
                IF ( netcdf_data_format > 4 )  THEN
!
!--                Set no fill for every variable to increase performance.
                   nc_stat = NF90_DEF_VAR_FILL( id_set_xy(av),     &
                                                id_var_do2d(av,i), &
                                                NF90_NOFILL, 0 )
                   CALL netcdf_handle_error( 'netcdf_define_header', 533 )
!
!--                Set collective io operations for parallel io
                   nc_stat = NF90_VAR_PAR_ACCESS( id_set_xy(av),     &
                                                  id_var_do2d(av,i), &
                                                  NF90_COLLECTIVE )
                   CALL netcdf_handle_error( 'netcdf_define_header', 448 )
                ENDIF
#endif
                var_list = TRIM( var_list) // TRIM( do2d(av,i) ) // ';'

             ENDIF

             i = i + 1

          ENDDO

!
!--       No arrays to output. Close the netcdf file and return.
          IF ( i == 1 )  RETURN

!
!--       Write the list of variables as global attribute (this is used by
!--       restart runs and by combine_plot_fields)
          nc_stat = NF90_PUT_ATT( id_set_xy(av), NF90_GLOBAL, 'VAR_LIST', &
                                  var_list )
          CALL netcdf_handle_error( 'netcdf_define_header', 121 )

!
!--       Set general no fill, otherwise the performance drops significantly for
!--       parallel output.
          nc_stat = NF90_SET_FILL( id_set_xy(av), NF90_NOFILL, oldmode )
          CALL netcdf_handle_error( 'netcdf_define_header', 529 )

!
!--       Leave netCDF define mode
          nc_stat = NF90_ENDDEF( id_set_xy(av) )
          CALL netcdf_handle_error( 'netcdf_define_header', 122 )

!
!--       These data are only written by PE0 for parallel output to increase
!--       the performance.
          IF ( myid == 0  .OR.  netcdf_data_format < 5 )  THEN

!
!--          Write axis data: z_xy, x, y
             ALLOCATE( netcdf_data(1:ns) )

!
!--          Write zu data
             DO  i = 1, ns
                IF( section(i,1) == -1 )  THEN
                   netcdf_data(i) = -1.0_wp  ! section averaged along z
                ELSE
                   netcdf_data(i) = zu( section(i,1) )
                ENDIF
             ENDDO
             nc_stat = NF90_PUT_VAR( id_set_xy(av), id_var_zu_xy(av), &
                                     netcdf_data, start = (/ 1 /),    &
                                     count = (/ ns /) )
             CALL netcdf_handle_error( 'netcdf_define_header', 123 )

!
!--          Write zw data
             DO  i = 1, ns
                IF( section(i,1) == -1 )  THEN
                   netcdf_data(i) = -1.0_wp  ! section averaged along z
                ELSE
                   netcdf_data(i) = zw( section(i,1) )
                ENDIF
             ENDDO
             nc_stat = NF90_PUT_VAR( id_set_xy(av), id_var_zw_xy(av), &
                                     netcdf_data, start = (/ 1 /),    &
                                     count = (/ ns /) )
             CALL netcdf_handle_error( 'netcdf_define_header', 124 )

!
!--          Write zs data
             IF ( land_surface )  THEN
                ns_do = 0
                DO  i = 1, ns
                   IF( section(i,1) == -1 )  THEN
                      netcdf_data(i) = 1.0_wp  ! section averaged along z
                      ns_do = ns_do + 1
                   ELSEIF ( section(i,1) < nzs )  THEN
                      netcdf_data(i) = - zs( section(i,1) )
                      ns_do = ns_do + 1
                   ENDIF
                ENDDO

                nc_stat = NF90_PUT_VAR( id_set_xy(av), id_var_zs_xy(av), &
                                        netcdf_data(1:ns_do), start = (/ 1 /),    &
                                        count = (/ ns_do /) )
                CALL netcdf_handle_error( 'netcdf_define_header', 124 )

             ENDIF

!
!--          Write gridpoint number data
             netcdf_data(1:ns) = section(1:ns,1)
             nc_stat = NF90_PUT_VAR( id_set_xy(av), id_var_ind_z_xy(av), &
                                     netcdf_data, start = (/ 1 /),       &
                                     count = (/ ns /) )
             CALL netcdf_handle_error( 'netcdf_define_header', 125 )

             DEALLOCATE( netcdf_data )

!
!--          Write the cross section height u*, t*
             nc_stat = NF90_PUT_VAR( id_set_xy(av), id_var_zu1_xy(av), &
                                     (/ zu(nzb+1) /), start = (/ 1 /), &
                                     count = (/ 1 /) )
             CALL netcdf_handle_error( 'netcdf_define_header', 126 )

!
!--          Write data for x (shifted by +dx/2) and xu axis
             ALLOCATE( netcdf_data(0:nx) )

             DO  i = 0, nx
                netcdf_data(i) = ( i + 0.5_wp ) * dx
             ENDDO

             nc_stat = NF90_PUT_VAR( id_set_xy(av), id_var_x_xy(av), &
                                     netcdf_data, start = (/ 1 /),   &
                                     count = (/ nx+1 /) )
             CALL netcdf_handle_error( 'netcdf_define_header', 127 )

             DO  i = 0, nx
                netcdf_data(i) = i * dx
             ENDDO

             nc_stat = NF90_PUT_VAR( id_set_xy(av), id_var_xu_xy(av), &
                                     netcdf_data, start = (/ 1 /),    &
                                     count = (/ nx+1 /) )
             CALL netcdf_handle_error( 'netcdf_define_header', 367 )

             DEALLOCATE( netcdf_data )

!
!--          Write data for y (shifted by +dy/2) and yv axis
             ALLOCATE( netcdf_data(0:ny+1) )

             DO  i = 0, ny
                netcdf_data(i) = ( i + 0.5_wp ) * dy
             ENDDO

             nc_stat = NF90_PUT_VAR( id_set_xy(av), id_var_y_xy(av), &
                                     netcdf_data, start = (/ 1 /),   &
                                     count = (/ ny+1 /))
             CALL netcdf_handle_error( 'netcdf_define_header', 128 )

             DO  i = 0, ny
                netcdf_data(i) = i * dy
             ENDDO

             nc_stat = NF90_PUT_VAR( id_set_xy(av), id_var_yv_xy(av), &
                                     netcdf_data, start = (/ 1 /),    &
                                     count = (/ ny+1 /))
             CALL netcdf_handle_error( 'netcdf_define_header', 368 )

             DEALLOCATE( netcdf_data )
!
!--          Write UTM coordinates
             IF ( rotation_angle == 0.0_wp )  THEN
!
!--             1D in case of no rotation
                cos_rot_angle = COS( rotation_angle * pi / 180.0_wp )
!
!--             x coordinates
                ALLOCATE( netcdf_data(0:nx) )
                DO  k = 0, 2
!
!--                Scalar grid points
                   IF ( k == 0 )  THEN
                      shift_x = 0.5
!
!--                u grid points
                   ELSEIF ( k == 1 )  THEN
                      shift_x = 0.0
!
!--                v grid points
                   ELSEIF ( k == 2 )  THEN
                      shift_x = 0.5
                   ENDIF

                   DO  i = 0, nx
                     netcdf_data(i) = init_model%origin_x                      &
                                    + cos_rot_angle * ( i + shift_x ) * dx
                   ENDDO

                   nc_stat = NF90_PUT_VAR( id_set_xy(av), id_var_eutm_xy(k,av),&
                                           netcdf_data, start = (/ 1 /),       &
                                           count = (/ nx+1 /) )
                   CALL netcdf_handle_error( 'netcdf_define_header', 555 )

                ENDDO
                DEALLOCATE( netcdf_data )
!
!--             y coordinates
                ALLOCATE( netcdf_data(0:ny) )
                DO  k = 0, 2
!
!--                Scalar grid points
                   IF ( k == 0 )  THEN
                      shift_y = 0.5
!
!--                u grid points
                   ELSEIF ( k == 1 )  THEN
                      shift_y = 0.5
!
!--                v grid points
                   ELSEIF ( k == 2 )  THEN
                      shift_y = 0.0
                   ENDIF

                   DO  j = 0, ny
                      netcdf_data(j) = init_model%origin_y                     &
                                     + cos_rot_angle * ( j + shift_y ) * dy
                   ENDDO

                   nc_stat = NF90_PUT_VAR( id_set_xy(av), id_var_nutm_xy(k,av),&
                                           netcdf_data, start = (/ 1 /),       &
                                           count = (/ ny+1 /) )
                   CALL netcdf_handle_error( 'netcdf_define_header', 556 )

                ENDDO
                DEALLOCATE( netcdf_data )

             ELSE
!
!--             2D in case of rotation
                ALLOCATE( netcdf_data_2d(0:nx,0:ny) )
                cos_rot_angle = COS( rotation_angle * pi / 180.0_wp )
                sin_rot_angle = SIN( rotation_angle * pi / 180.0_wp )

                DO  k = 0, 2
!
!--               Scalar grid points
                  IF ( k == 0 )  THEN
                     shift_x = 0.5 ; shift_y = 0.5
!
!--               u grid points
                  ELSEIF ( k == 1 )  THEN
                     shift_x = 0.0 ; shift_y = 0.5
!
!--               v grid points
                  ELSEIF ( k == 2 )  THEN
                     shift_x = 0.5 ; shift_y = 0.0
                  ENDIF

                  DO  j = 0, ny
                     DO  i = 0, nx
                        netcdf_data_2d(i,j) = init_model%origin_x                   &
                                            + cos_rot_angle * ( i + shift_x ) * dx  &
                                            + sin_rot_angle * ( j + shift_y ) * dy
                     ENDDO
                  ENDDO

                  nc_stat = NF90_PUT_VAR( id_set_xy(av), id_var_eutm_xy(k,av),  &
                                          netcdf_data_2d, start = (/ 1, 1 /),   &
                                          count = (/ nx+1, ny+1 /) )
                  CALL netcdf_handle_error( 'netcdf_define_header', 555 )

                  DO  j = 0, ny
                     DO  i = 0, nx
                        netcdf_data_2d(i,j) = init_model%origin_y                   &
                                            - sin_rot_angle * ( i + shift_x ) * dx  &
                                            + cos_rot_angle * ( j + shift_y ) * dy
                     ENDDO
                  ENDDO

                  nc_stat = NF90_PUT_VAR( id_set_xy(av), id_var_nutm_xy(k,av),  &
                                          netcdf_data_2d, start = (/ 1, 1 /),   &
                                          count = (/ nx+1, ny+1 /) )
                  CALL netcdf_handle_error( 'netcdf_define_header', 556 )

                ENDDO
                DEALLOCATE( netcdf_data_2d )
             ENDIF

          ENDIF
!
!--       Write lon and lat data. Only for parallel output.
          IF ( netcdf_data_format > 4 )  THEN

             ALLOCATE( lat(nxl:nxr,nys:nyn) )
             ALLOCATE( lon(nxl:nxr,nys:nyn) )
             cos_rot_angle = COS( rotation_angle * pi / 180.0_wp )
             sin_rot_angle = SIN( rotation_angle * pi / 180.0_wp )

             DO  k = 0, 2
!
!--             Scalar grid points
                IF ( k == 0 )  THEN
                   shift_x = 0.5 ; shift_y = 0.5
!
!--             u grid points
                ELSEIF ( k == 1 )  THEN
                   shift_x = 0.0 ; shift_y = 0.5
!
!--             v grid points
                ELSEIF ( k == 2 )  THEN
                   shift_x = 0.5 ; shift_y = 0.0
                ENDIF

                DO  j = nys, nyn
                   DO  i = nxl, nxr
                      eutm = init_model%origin_x                   &
                           + cos_rot_angle * ( i + shift_x ) * dx  &
                           + sin_rot_angle * ( j + shift_y ) * dy
                      nutm = init_model%origin_y                   &
                           - sin_rot_angle * ( i + shift_x ) * dx  &
                           + cos_rot_angle * ( j + shift_y ) * dy

                      CALL  convert_utm_to_geographic( crs_list,          &
                                                       eutm, nutm,        &
                                                       lon(i,j), lat(i,j) )
                   ENDDO
                ENDDO

                nc_stat = NF90_PUT_VAR( id_set_xy(av), id_var_lon_xy(k,av), &
                                     lon, start = (/ nxl+1, nys+1 /),       &
                                     count = (/ nxr-nxl+1, nyn-nys+1 /) )
                CALL netcdf_handle_error( 'netcdf_define_header', 556 )

                nc_stat = NF90_PUT_VAR( id_set_xy(av), id_var_lat_xy(k,av), &
                                     lat, start = (/ nxl+1, nys+1 /),       &
                                     count = (/ nxr-nxl+1, nyn-nys+1 /) )
                CALL netcdf_handle_error( 'netcdf_define_header', 556 )
             ENDDO

             DEALLOCATE( lat )
             DEALLOCATE( lon )

          ENDIF
!
!--       In case of non-flat topography write height information. Only for
!--       parallel netcdf output.
          IF ( TRIM( topography ) /= 'flat'  .AND.                             &
               netcdf_data_format > 4  )  THEN

!             IF ( nxr == nx  .AND.  nyn /= ny )  THEN
!                nc_stat = NF90_PUT_VAR( id_set_xy(av), id_var_zusi_xy(av),     &
!                                        zu_s_inner(nxl:nxr+1,nys:nyn),         &
!                                        start = (/ nxl+1, nys+1 /),            &
!                                        count = (/ nxr-nxl+2, nyn-nys+1 /) )
!             ELSEIF ( nxr /= nx  .AND.  nyn == ny )  THEN
!                nc_stat = NF90_PUT_VAR( id_set_xy(av), id_var_zusi_xy(av),     &
!                                        zu_s_inner(nxl:nxr,nys:nyn+1),         &
!                                        start = (/ nxl+1, nys+1 /),            &
!                                        count = (/ nxr-nxl+1, nyn-nys+2 /) )
!             ELSEIF ( nxr == nx  .AND.  nyn == ny )  THEN
!                nc_stat = NF90_PUT_VAR( id_set_xy(av), id_var_zusi_xy(av),     &
!                                        zu_s_inner(nxl:nxr+1,nys:nyn+1),       &
!                                        start = (/ nxl+1, nys+1 /),            &
!                                        count = (/ nxr-nxl+2, nyn-nys+2 /) )
!             ELSE
                nc_stat = NF90_PUT_VAR( id_set_xy(av), id_var_zusi_xy(av),     &
                                        zu_s_inner(nxl:nxr,nys:nyn),           &
                                        start = (/ nxl+1, nys+1 /),            &
                                        count = (/ nxr-nxl+1, nyn-nys+1 /) )
!             ENDIF
             CALL netcdf_handle_error( 'netcdf_define_header', 427 )

!             IF ( nxr == nx  .AND.  nyn /= ny )  THEN
!                nc_stat = NF90_PUT_VAR( id_set_xy(av), id_var_zwwi_xy(av),     &
!                                        zw_w_inner(nxl:nxr+1,nys:nyn),         &
!                                        start = (/ nxl+1, nys+1 /),            &
!                                        count = (/ nxr-nxl+2, nyn-nys+1 /) )
!             ELSEIF ( nxr /= nx  .AND.  nyn == ny )  THEN
!                nc_stat = NF90_PUT_VAR( id_set_xy(av), id_var_zwwi_xy(av),     &
!                                        zw_w_inner(nxl:nxr,nys:nyn+1),         &
!                                        start = (/ nxl+1, nys+1 /),            &
!                                        count = (/ nxr-nxl+1, nyn-nys+2 /) )
!             ELSEIF ( nxr == nx  .AND.  nyn == ny )  THEN
!                nc_stat = NF90_PUT_VAR( id_set_xy(av), id_var_zwwi_xy(av),     &
!                                        zw_w_inner(nxl:nxr+1,nys:nyn+1),       &
!                                        start = (/ nxl+1, nys+1 /),            &
!                                        count = (/ nxr-nxl+2, nyn-nys+2 /) )
!             ELSE
                nc_stat = NF90_PUT_VAR( id_set_xy(av), id_var_zwwi_xy(av),     &
                                        zw_w_inner(nxl:nxr,nys:nyn),           &
                                        start = (/ nxl+1, nys+1 /),            &
                                        count = (/ nxr-nxl+1, nyn-nys+1 /) )
!             ENDIF
             CALL netcdf_handle_error( 'netcdf_define_header', 428 )

          ENDIF

       CASE ( 'xy_ext' )

!
!--       Get the list of variables and compare with the actual run.
!--       First var_list_old has to be reset, since GET_ATT does not assign
!--       trailing blanks.
          var_list_old = ' '
          nc_stat = NF90_GET_ATT( id_set_xy(av), NF90_GLOBAL, 'VAR_LIST', &
                                  var_list_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 129 )

          var_list = ';'
          i = 1
          DO WHILE ( do2d(av,i)(1:1) /= ' ' )
             IF ( INDEX( do2d(av,i), 'xy' ) /= 0 )  THEN
                var_list = TRIM( var_list ) // TRIM( do2d(av,i) ) // ';'
             ENDIF
             i = i + 1
          ENDDO

          IF ( av == 0 )  THEN
             var = '(xy)'
          ELSE
             var = '(xy_av)'
          ENDIF

          IF ( TRIM( var_list ) /= TRIM( var_list_old ) )  THEN
             message_string = 'netCDF file for cross-sections ' //           &
                              TRIM( var ) // ' from previous run found,' //  &
                              '&but this file cannot be extended due to' //  &
                              ' variable mismatch.' //                       &
                              '&New file is created instead.'
             CALL message( 'define_netcdf_header', 'PA0249', 0, 1, 0, 6, 0 )
             extend = .FALSE.
             RETURN
          ENDIF

!
!--       Calculate the number of current sections
          ns = 1
          DO WHILE ( section(ns,1) /= -9999  .AND.  ns <= 100 )
             ns = ns + 1
          ENDDO
          ns = ns - 1

!
!--       Get and compare the number of horizontal cross sections
          nc_stat = NF90_INQ_VARID( id_set_xy(av), 'zu_xy', id_var_zu_xy(av) )
          CALL netcdf_handle_error( 'netcdf_define_header', 130 )

          nc_stat = NF90_INQUIRE_VARIABLE( id_set_xy(av), id_var_zu_xy(av), &
                                           dimids = id_dim_zu_xy_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 131 )
          id_dim_zu_xy(av) = id_dim_zu_xy_old(1)

          nc_stat = NF90_INQUIRE_DIMENSION( id_set_xy(av), id_dim_zu_xy(av), &
                                            len = ns_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 132 )

          IF ( ns /= ns_old )  THEN
             message_string = 'netCDF file for cross-sections ' //          &
                              TRIM( var ) // ' from previous run found,' // &
                              '&but this file cannot be extended due to' // &
                              ' mismatch in number of' //                   &
                              ' cross sections.' //                         &
                              '&New file is created instead.'
             CALL message( 'define_netcdf_header', 'PA0250', 0, 1, 0, 6, 0 )
             extend = .FALSE.
             RETURN
          ENDIF

!
!--       Get and compare the heights of the cross sections
          ALLOCATE( netcdf_data(1:ns_old) )

          nc_stat = NF90_GET_VAR( id_set_xy(av), id_var_zu_xy(av), netcdf_data )
          CALL netcdf_handle_error( 'netcdf_define_header', 133 )

          DO  i = 1, ns
             IF ( section(i,1) /= -1 )  THEN
                IF ( zu(section(i,1)) /= netcdf_data(i) )  THEN
                   message_string = 'netCDF file for cross-sections ' //       &
                               TRIM( var ) // ' from previous run found,' //   &
                               ' but this file cannot be extended' //          &
                               ' due to mismatch in cross' //                  &
                               ' section levels.' //                           &
                               ' New file is created instead.'
                   CALL message( 'define_netcdf_header', 'PA0251',             &
                                                                 0, 1, 0, 6, 0 )
                   extend = .FALSE.
                   RETURN
                ENDIF
             ELSE
                IF ( -1.0_wp /= netcdf_data(i) )  THEN
                   message_string = 'netCDF file for cross-sections ' //       &
                               TRIM( var ) // ' from previous run found,' //   &
                               ' but this file cannot be extended' //          &
                               ' due to mismatch in cross' //                  &
                               ' section levels.' //                           &
                               ' New file is created instead.'
                   CALL message( 'define_netcdf_header', 'PA0251',             &
                                                                 0, 1, 0, 6, 0 )
                   extend = .FALSE.
                   RETURN
                ENDIF
             ENDIF
          ENDDO

          DEALLOCATE( netcdf_data )

!
!--       Get the id of the time coordinate (unlimited coordinate) and its
!--       last index on the file. The next time level is do2d..count+1.
!--       The current time must be larger than the last output time
!--       on the file.
          nc_stat = NF90_INQ_VARID( id_set_xy(av), 'time', id_var_time_xy(av) )
          CALL netcdf_handle_error( 'netcdf_define_header', 134 )

          nc_stat = NF90_INQUIRE_VARIABLE( id_set_xy(av), id_var_time_xy(av), &
                                           dimids = id_dim_time_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 135 )
          id_dim_time_xy(av) = id_dim_time_old(1)

          nc_stat = NF90_INQUIRE_DIMENSION( id_set_xy(av), id_dim_time_xy(av), &
                                            len = ntime_count )
          CALL netcdf_handle_error( 'netcdf_define_header', 136 )

!
!--       For non-parallel output use the last output time level of the netcdf
!--       file because the time dimension is unlimited. In case of parallel
!--       output the variable ntime_count could get the value of 9*10E36 because
!--       the time dimension is limited.
          IF ( netcdf_data_format < 5 ) do2d_xy_time_count(av) = ntime_count

          nc_stat = NF90_GET_VAR( id_set_xy(av), id_var_time_xy(av),           &
                                  last_time_coordinate,                        &
                                  start = (/ do2d_xy_time_count(av) /),        &
                                  count = (/ 1 /) )
          CALL netcdf_handle_error( 'netcdf_define_header', 137 )

          IF ( last_time_coordinate(1) >= simulated_time )  THEN
             message_string = 'netCDF file for cross sections ' //             &
                              TRIM( var ) // ' from previous run found,' //    &
                              '&but this file cannot be extended becaus' //    &
                              'e the current output time' //                   &
                              '&is less or equal than the last output t' //    &
                              'ime on this file.' //                           &
                              '&New file is created instead.'
             CALL message( 'define_netcdf_header', 'PA0252', 0, 1, 0, 6, 0 )
             do2d_xy_time_count(av) = 0
             extend = .FALSE.
             RETURN
          ENDIF

          IF ( netcdf_data_format > 4 )  THEN
!
!--          Check if the needed number of output time levels is increased
!--          compared to the number of time levels in the existing file.
             IF ( ntdim_2d_xy(av) > ntime_count )  THEN
                message_string = 'netCDF file for cross sections ' //          &
                                 TRIM( var ) // ' from previous run found,' // &
                                 '&but this file cannot be extended becaus' // &
                                 'e the number of output time levels has b' // &
                                 'een increased compared to the previous s' // &
                                 'imulation.' //                               &
                                 '&New file is created instead.'
                CALL message( 'define_netcdf_header', 'PA0389', 0, 1, 0, 6, 0 )
                do2d_xy_time_count(av) = 0
                extend = .FALSE.
!
!--             Recalculate the needed time levels for the new file.
                IF ( av == 0 )  THEN
                   ntdim_2d_xy(0) = CEILING(                            &
                           ( end_time - MAX( skip_time_do2d_xy,         &
                                             simulated_time_at_begin )  &
                           ) / dt_do2d_xy )
                   IF ( do2d_at_begin )  ntdim_2d_xy(0) = ntdim_2d_xy(0) + 1
                ELSE
                   ntdim_2d_xy(1) = CEILING(                            &
                           ( end_time - MAX( skip_time_data_output_av,  &
                                             simulated_time_at_begin )  &
                           ) / dt_data_output_av )
                ENDIF
                RETURN
             ENDIF
          ENDIF

!
!--       Dataset seems to be extendable.
!--       Now get the variable ids.
          i = 1
          DO WHILE ( do2d(av,i)(1:1) /= ' ' )
             IF ( INDEX( do2d(av,i), 'xy' ) /= 0 )  THEN
                nc_stat = NF90_INQ_VARID( id_set_xy(av), do2d(av,i), &
                                          id_var_do2d(av,i) )
                CALL netcdf_handle_error( 'netcdf_define_header', 138 )
#if defined( __netcdf4_parallel )
!
!--             Set collective io operations for parallel io
                IF ( netcdf_data_format > 4 )  THEN
                   nc_stat = NF90_VAR_PAR_ACCESS( id_set_xy(av),     &
                                                  id_var_do2d(av,i), &
                                                  NF90_COLLECTIVE )
                   CALL netcdf_handle_error( 'netcdf_define_header', 454 )
                ENDIF
#endif
             ENDIF
             i = i + 1
          ENDDO

!
!--       Update the title attribute on file
!--       In order to avoid 'data mode' errors if updated attributes are larger
!--       than their original size, NF90_PUT_ATT is called in 'define mode'
!--       enclosed by NF90_REDEF and NF90_ENDDEF calls. This implies a possible
!--       performance loss due to data copying; an alternative strategy would be
!--       to ensure equal attribute size in a job chain. Maybe revise later.
          IF ( av == 0 )  THEN
             time_average_text = ' '
          ELSE
             WRITE (time_average_text, '('', '',F7.1,'' s average'')') &
                                                            averaging_interval
          ENDIF
          nc_stat = NF90_REDEF( id_set_xy(av) )
          CALL netcdf_handle_error( 'netcdf_define_header', 431 )
          nc_stat = NF90_PUT_ATT( id_set_xy(av), NF90_GLOBAL, 'title',         &
                                  TRIM( run_description_header ) //            &
                                  TRIM( time_average_text ) )
          CALL netcdf_handle_error( 'netcdf_define_header', 139 )
          nc_stat = NF90_ENDDEF( id_set_xy(av) )
          CALL netcdf_handle_error( 'netcdf_define_header', 432 )
          message_string = 'netCDF file for cross-sections ' //                &
                            TRIM( var ) // ' from previous run found.' //      &
                           '&This file will be extended.'
          CALL message( 'define_netcdf_header', 'PA0253', 0, 0, 0, 6, 0 )


       CASE ( 'xz_new' )

!
!--       Define some global attributes of the dataset
          IF ( av == 0 )  THEN
             CALL netcdf_create_global_atts( id_set_xz(av), 'xz', TRIM( run_description_header ), 140 )
             time_average_text = ' '
          ELSE
             CALL netcdf_create_global_atts( id_set_xz(av), 'xz_av', TRIM( run_description_header ), 140 )
             WRITE ( time_average_text,'(F7.1,'' s avg'')' )  averaging_interval
             nc_stat = NF90_PUT_ATT( id_set_xz(av), NF90_GLOBAL, 'time_avg',   &
                                     TRIM( time_average_text ) )
             CALL netcdf_handle_error( 'netcdf_define_header', 141 )
          ENDIF

!
!--       Define time coordinate for xz sections.
!--       For parallel output the time dimensions has to be limited, otherwise
!--       the performance drops significantly.
          IF ( netcdf_data_format < 5 )  THEN
             CALL netcdf_create_dim( id_set_xz(av), 'time', NF90_UNLIMITED,    &
                                     id_dim_time_xz(av), 142 )
          ELSE
             CALL netcdf_create_dim( id_set_xz(av), 'time', ntdim_2d_xz(av),   &
                                     id_dim_time_xz(av), 525 )
          ENDIF

          CALL netcdf_create_var( id_set_xz(av), (/ id_dim_time_xz(av) /),     &
                                  'time', NF90_DOUBLE, id_var_time_xz(av),     &
                                  'seconds', 'time', 143, 144, 000 )
          CALL netcdf_create_att( id_set_xz(av), id_var_time_xz(av), 'standard_name', 'time', 000)
          CALL netcdf_create_att( id_set_xz(av), id_var_time_xz(av), 'axis', 'T', 000)
!
!--       Define the spatial dimensions and coordinates for xz-sections.
!--       First, determine the number of vertical sections to be written.
          IF ( section(1,2) == -9999 )  THEN
             RETURN
          ELSE
             ns = 1
             DO WHILE ( section(ns,2) /= -9999  .AND.  ns <= 100 )
                ns = ns + 1
             ENDDO
             ns = ns - 1
          ENDIF

!
!--       Define y-axis (for scalar position)
          CALL netcdf_create_dim( id_set_xz(av), 'y_xz', ns, id_dim_y_xz(av),  &
                                  145 )
          CALL netcdf_create_var( id_set_xz(av), (/ id_dim_y_xz(av) /),        &
                                  'y_xz', NF90_DOUBLE, id_var_y_xz(av),        &
                                  'meters', '', 146, 147, 000 )
          CALL netcdf_create_att( id_set_xz(av), id_var_y_xz(av), 'axis',      &
                                  'Y', 000)
!
!--       Define y-axis (for v position)
          CALL netcdf_create_dim( id_set_xz(av), 'yv_xz', ns,                  &
                                  id_dim_yv_xz(av), 369 )
          CALL netcdf_create_var( id_set_xz(av), (/ id_dim_yv_xz(av) /),       &
                                  'yv_xz', NF90_DOUBLE, id_var_yv_xz(av),      &
                                  'meters', '', 370, 371, 000 )
          CALL netcdf_create_att( id_set_xz(av), id_var_yv_xz(av), 'axis',     &
                                  'Y', 000)
!
!--       Define a variable to store the layer indices of the vertical cross
!--       sections
          CALL netcdf_create_var( id_set_xz(av), (/ id_dim_y_xz(av) /),        &
                                  'ind_y_xz', NF90_DOUBLE,                     &
                                  id_var_ind_y_xz(av), 'gridpoints', '', 148,  &
                                  149, 000 )
!
!--       Define x-axis (for scalar position)
          CALL netcdf_create_dim( id_set_xz(av), 'x', nx+1, id_dim_x_xz(av),   &
                                  150 )
          CALL netcdf_create_var( id_set_xz(av), (/ id_dim_x_xz(av) /), 'x',   &
                                  NF90_DOUBLE, id_var_x_xz(av), 'meters', '',  &
                                  151, 152, 000 )
          CALL netcdf_create_att( id_set_xz(av), id_var_x_xz(av), 'axis',      &
                                  'X', 000)
!
!--       Define x-axis (for u position)
          CALL netcdf_create_dim( id_set_xz(av), 'xu', nx+1, id_dim_xu_xz(av), &
                                  372 )
          CALL netcdf_create_var( id_set_xz(av), (/ id_dim_xu_xz(av) /), 'xu', &
                                  NF90_DOUBLE, id_var_xu_xz(av), 'meters', '', &
                                  373, 374, 000 )
          CALL netcdf_create_att( id_set_xz(av), id_var_xu_xz(av), 'axis',     &
                                  'X', 000)
                                  
!
!--       Define the three z-axes (zu, zw, and zs)
          CALL netcdf_create_dim( id_set_xz(av), 'zu', nz+2, id_dim_zu_xz(av), &
                                  153 )
          CALL netcdf_create_var( id_set_xz(av), (/ id_dim_zu_xz(av) /), 'zu', &
                                  NF90_DOUBLE, id_var_zu_xz(av), 'meters', '', &
                                  154, 155, 000 )
          CALL netcdf_create_att( id_set_xz(av), id_var_zu_xz(av), 'axis',     &
                                  'Z', 000)
                                  
          CALL netcdf_create_dim( id_set_xz(av), 'zw', nz+2, id_dim_zw_xz(av), &
                                  156 )
          CALL netcdf_create_var( id_set_xz(av), (/ id_dim_zw_xz(av) /), 'zw', &
                                  NF90_DOUBLE, id_var_zw_xz(av), 'meters', '', &
                                  157, 158, 000 )
          CALL netcdf_create_att( id_set_xz(av), id_var_zw_xz(av), 'axis',     &
                                  'Z', 000)
!
!--       Define UTM and geographic coordinates
          CALL define_geo_coordinates( id_set_xz(av),         &
                  (/ id_dim_x_xz(av), id_dim_xu_xz(av) /),    &
                  (/ id_dim_y_xz(av), id_dim_yv_xz(av) /),    &
                  id_var_eutm_xz(:,av), id_var_nutm_xz(:,av), &
                  id_var_lat_xz(:,av), id_var_lon_xz(:,av)    )
!
!--       Define coordinate-reference system
          CALL netcdf_create_crs( id_set_xz(av), 000 )

          IF ( land_surface )  THEN

             CALL netcdf_create_dim( id_set_xz(av), 'zs', nzs,                 &
                                     id_dim_zs_xz(av), 542 )
             CALL netcdf_create_var( id_set_xz(av), (/ id_dim_zs_xz(av) /),    &
                                     'zs', NF90_DOUBLE, id_var_zs_xz(av),      &
                                     'meters', '', 543, 544, 000 )
             CALL netcdf_create_att( id_set_xz(av), id_var_zs_xz(av), 'axis',  &
                                     'Z', 000)

          ENDIF

!
!--       Define the variables
          var_list = ';'
          i = 1

          DO WHILE ( do2d(av,i)(1:1) /= ' ' )

             IF ( INDEX( do2d(av,i), 'xz' ) /= 0 )  THEN

!
!--             Check for the grid
                found = .FALSE.
                SELECT CASE ( do2d(av,i) )
!
!--                Most variables are defined on the zu grid
                   CASE ( 'e_xz', 'nc_xz', 'nr_xz', 'p_xz', 'pc_xz',           &
                          'pr_xz', 'prr_xz', 'q_xz', 'qc_xz',                  &
                          'ql_xz', 'ql_c_xz', 'ql_v_xz', 'ql_vp_xz', 'qr_xz',  &
                          'qv_xz', 's_xz',                                     &
                          'theta_xz', 'thetal_xz', 'thetav_xz' )

                      grid_x = 'x'
                      grid_y = 'y'
                      grid_z = 'zu'
!
!--                u grid
                   CASE ( 'u_xz' )

                      grid_x = 'xu'
                      grid_y = 'y'
                      grid_z = 'zu'
!
!--                v grid
                   CASE ( 'v_xz' )

                      grid_x = 'x'
                      grid_y = 'yv'
                      grid_z = 'zu'
!
!--                w grid
                   CASE ( 'w_xz' )

                      grid_x = 'x'
                      grid_y = 'y'
                      grid_z = 'zw'

                   CASE DEFAULT

!
!--                   Check for land surface quantities
                      IF ( land_surface )  THEN
                         CALL lsm_define_netcdf_grid( do2d(av,i), found,       &
                                                      grid_x, grid_y, grid_z )
                      ENDIF

                      IF ( .NOT. found )  THEN
                         CALL tcm_define_netcdf_grid( do2d(av,i), found,       &
                                                      grid_x, grid_y, grid_z )
                      ENDIF

!
!--                   Check for ocean quantities
                      IF ( .NOT. found  .AND.  ocean_mode )  THEN
                         CALL ocean_define_netcdf_grid( do2d(av,i), found,  &
                                                        grid_x, grid_y, grid_z )
                      ENDIF
!
!--                   Check for radiation quantities
                      IF ( .NOT. found  .AND.  radiation )  THEN
                         CALL radiation_define_netcdf_grid( do2d(av,i), found, &
                                                            grid_x, grid_y,    &
                                                            grid_z )
                      ENDIF
!
!--                   Check for SALSA quantities
                      IF ( .NOT. found  .AND.  salsa )  THEN
                         CALL salsa_define_netcdf_grid( do2d(av,i), found,     &
                                                        grid_x, grid_y, grid_z )
                      ENDIF

!
!--                   Check for gust module quantities
                      IF ( .NOT. found  .AND.  gust_module_enabled )  THEN
                         CALL gust_define_netcdf_grid( do2d(av,i), found,      &
                                                       grid_x, grid_y, grid_z )
                      ENDIF

!
!--                   Check for chemistry quantities
                      IF ( .NOT. found  .AND.  air_chemistry )  THEN
                         CALL chem_define_netcdf_grid( do2d(av,i), found,      &
                                                       grid_x, grid_y,         &
                                                       grid_z )
                      ENDIF

                      IF ( .NOT. found )                                       &
                         CALL doq_define_netcdf_grid( do2d(av,i), found,       &
                                                      grid_x, grid_y, grid_z )

!
!--                   Check for user-defined quantities
                      IF ( .NOT. found  .AND.  user_module_enabled )  THEN
                         CALL user_define_netcdf_grid( do2d(av,i), found,      &
                                                       grid_x, grid_y, grid_z )
                      ENDIF

                      IF ( .NOT. found )  THEN
                         WRITE ( message_string, * ) 'no grid defined for',    &
                                                ' variable ', TRIM( do2d(av,i) )
                         CALL message( 'define_netcdf_header', 'PA0244',       &
                                       0, 1, 0, 6, 0 )
                      ENDIF

                END SELECT

!
!--             Select the respective dimension ids
                IF ( grid_x == 'x' )  THEN
                   id_x = id_dim_x_xz(av)
                ELSEIF ( grid_x == 'xu' )  THEN
                   id_x = id_dim_xu_xz(av)
                ENDIF

                IF ( grid_y == 'y' )  THEN
                   id_y = id_dim_y_xz(av)
                ELSEIF ( grid_y == 'yv' )  THEN
                   id_y = id_dim_yv_xz(av)
                ENDIF

                IF ( grid_z == 'zu' )  THEN
                   id_z = id_dim_zu_xz(av)
                ELSEIF ( grid_z == 'zw' )  THEN
                   id_z = id_dim_zw_xz(av)
                ELSEIF ( grid_z == 'zs' )  THEN
                   id_z = id_dim_zs_xz(av)
                ENDIF

!
!--             Define the grid
                CALL netcdf_create_var( id_set_xz(av), (/ id_x, id_y, id_z,    &
                                        id_dim_time_xz(av) /), do2d(av,i),     &
                                        nc_precision(2), id_var_do2d(av,i),    &
                                        TRIM( do2d_unit(av,i) ), do2d(av,i),   &
                                        159, 160, 355, .TRUE. )

#if defined( __netcdf4_parallel )

                IF ( netcdf_data_format > 4 )  THEN
!
!--                Set no fill for every variable to increase performance.
                   nc_stat = NF90_DEF_VAR_FILL( id_set_xz(av),     &
                                                id_var_do2d(av,i), &
                                                NF90_NOFILL, 0 )
                   CALL netcdf_handle_error( 'netcdf_define_header', 534 )
!
!--                Set independent io operations for parallel io. Collective io
!--                is only allowed in case of a 1d-decomposition along x,
!--                because otherwise, not all PEs have output data.
                   IF ( npey == 1 )  THEN
                      nc_stat = NF90_VAR_PAR_ACCESS( id_set_xz(av),     &
                                                     id_var_do2d(av,i), &
                                                     NF90_COLLECTIVE )
                   ELSE
!
!--                   Test simulations showed that the output of cross sections
!--                   by all PEs in data_output_2d using NF90_COLLECTIVE is
!--                   faster than the output by the first row of PEs in
!--                   x-direction using NF90_INDEPENDENT.
                      nc_stat = NF90_VAR_PAR_ACCESS( id_set_xz(av),    &
                                                    id_var_do2d(av,i), &
                                                    NF90_COLLECTIVE )
!                      nc_stat = NF90_VAR_PAR_ACCESS( id_set_xz(av),     &
!                                                     id_var_do2d(av,i), &
!                                                     NF90_INDEPENDENT )
                   ENDIF
                   CALL netcdf_handle_error( 'netcdf_define_header', 449 )
                ENDIF
#endif
                var_list = TRIM( var_list ) // TRIM( do2d(av,i) ) // ';'

             ENDIF

             i = i + 1

          ENDDO

!
!--       No arrays to output. Close the netcdf file and return.
          IF ( i == 1 )  RETURN

!
!--       Write the list of variables as global attribute (this is used by
!--       restart runs and by combine_plot_fields)
          nc_stat = NF90_PUT_ATT( id_set_xz(av), NF90_GLOBAL, 'VAR_LIST', &
                                  var_list )
          CALL netcdf_handle_error( 'netcdf_define_header', 161 )

!
!--       Set general no fill, otherwise the performance drops significantly for
!--       parallel output.
          nc_stat = NF90_SET_FILL( id_set_xz(av), NF90_NOFILL, oldmode )
          CALL netcdf_handle_error( 'netcdf_define_header', 530 )

!
!--       Leave netCDF define mode
          nc_stat = NF90_ENDDEF( id_set_xz(av) )
          CALL netcdf_handle_error( 'netcdf_define_header', 162 )

!
!--       These data are only written by PE0 for parallel output to increase
!--       the performance.
          IF ( myid == 0  .OR.  netcdf_data_format < 5 )  THEN

!
!--          Write axis data: y_xz, x, zu, zw
             ALLOCATE( netcdf_data(1:ns) )

!
!--          Write y_xz data (shifted by +dy/2)
             DO  i = 1, ns
                IF( section(i,2) == -1 )  THEN
                   netcdf_data(i) = -1.0_wp  ! section averaged along y
                ELSE
                   netcdf_data(i) = ( section(i,2) + 0.5_wp ) * dy
                ENDIF
             ENDDO
             nc_stat = NF90_PUT_VAR( id_set_xz(av), id_var_y_xz(av), &
                                     netcdf_data, start = (/ 1 /),   &
                                     count = (/ ns /) )
             CALL netcdf_handle_error( 'netcdf_define_header', 163 )

!
!--          Write yv_xz data
             DO  i = 1, ns
                IF( section(i,2) == -1 )  THEN
                   netcdf_data(i) = -1.0_wp  ! section averaged along y
                ELSE
                   netcdf_data(i) = section(i,2) * dy
                ENDIF
             ENDDO
             nc_stat = NF90_PUT_VAR( id_set_xz(av), id_var_yv_xz(av), &
                                     netcdf_data, start = (/ 1 /),    &
                                     count = (/ ns /) )
             CALL netcdf_handle_error( 'netcdf_define_header', 375 )

!
!--          Write gridpoint number data
             netcdf_data(1:ns) = section(1:ns,2)
             nc_stat = NF90_PUT_VAR( id_set_xz(av), id_var_ind_y_xz(av), &
                                     netcdf_data, start = (/ 1 /),       &
                                     count = (/ ns /) )
             CALL netcdf_handle_error( 'netcdf_define_header', 164 )


             DEALLOCATE( netcdf_data )

!
!--          Write data for x (shifted by +dx/2) and xu axis
             ALLOCATE( netcdf_data(0:nx) )

             DO  i = 0, nx
                netcdf_data(i) = ( i + 0.5_wp ) * dx
             ENDDO

             nc_stat = NF90_PUT_VAR( id_set_xz(av), id_var_x_xz(av), &
                                     netcdf_data, start = (/ 1 /),   &
                                     count = (/ nx+1 /) )
             CALL netcdf_handle_error( 'netcdf_define_header', 165 )

             DO  i = 0, nx
                netcdf_data(i) = i * dx
             ENDDO

             nc_stat = NF90_PUT_VAR( id_set_xz(av), id_var_xu_xz(av), &
                                     netcdf_data, start = (/ 1 /),    &
                                     count = (/ nx+1 /) )
             CALL netcdf_handle_error( 'netcdf_define_header', 377 )

             DEALLOCATE( netcdf_data )

!
!--          Write zu and zw data (vertical axes)
             ALLOCATE( netcdf_data(0:nz+1) )

             netcdf_data(0:nz+1) = zu(nzb:nzt+1)
             nc_stat = NF90_PUT_VAR( id_set_xz(av), id_var_zu_xz(av), &
                                     netcdf_data, start = (/ 1 /),    &
                                     count = (/ nz+2 /) )
             CALL netcdf_handle_error( 'netcdf_define_header', 166 )

             netcdf_data(0:nz+1) = zw(nzb:nzt+1)
             nc_stat = NF90_PUT_VAR( id_set_xz(av), id_var_zw_xz(av), &
                                     netcdf_data, start = (/ 1 /),    &
                                     count = (/ nz+2 /) )
             CALL netcdf_handle_error( 'netcdf_define_header', 167 )

!
!--          Write zs data
             IF ( land_surface )  THEN
                netcdf_data(0:nzs-1) = - zs(nzb_soil:nzt_soil)
                nc_stat = NF90_PUT_VAR( id_set_xz(av), id_var_zs_xz(av), &
                                        netcdf_data(0:nzs), start = (/ 1 /),    &
                                        count = (/ nzt_soil-nzb_soil+1 /) )
               CALL netcdf_handle_error( 'netcdf_define_header', 548 )
             ENDIF

             DEALLOCATE( netcdf_data )
!
!--          Write UTM coordinates
             IF ( rotation_angle == 0.0_wp )  THEN
!
!--             1D in case of no rotation
                cos_rot_angle = COS( rotation_angle * pi / 180.0_wp )
!
!--             x coordinates
                ALLOCATE( netcdf_data(0:nx) )
                DO  k = 0, 2
!
!--                Scalar grid points
                   IF ( k == 0 )  THEN
                      shift_x = 0.5
!
!--                u grid points
                   ELSEIF ( k == 1 )  THEN
                      shift_x = 0.0
!
!--                v grid points
                   ELSEIF ( k == 2 )  THEN
                      shift_x = 0.5
                   ENDIF

                   DO  i = 0, nx
                     netcdf_data(i) = init_model%origin_x                      &
                                    + cos_rot_angle * ( i + shift_x ) * dx
                   ENDDO

                   nc_stat = NF90_PUT_VAR( id_set_xz(av), id_var_eutm_xz(k,av),&
                                           netcdf_data, start = (/ 1 /),       &
                                           count = (/ nx+1 /) )
                   CALL netcdf_handle_error( 'netcdf_define_header', 555 )

                ENDDO
                DEALLOCATE( netcdf_data )
!
!--             y coordinates
                ALLOCATE( netcdf_data(1:ns) )
                DO  k = 0, 2
!
!--                Scalar grid points
                   IF ( k == 0 )  THEN
                      shift_y = 0.5
!
!--                u grid points
                   ELSEIF ( k == 1 )  THEN
                      shift_y = 0.5
!
!--                v grid points
                   ELSEIF ( k == 2 )  THEN
                      shift_y = 0.0
                   ENDIF

                   DO  i = 1, ns
                      IF( section(i,2) == -1 )  THEN
                         netcdf_data(i) = -1.0_wp  ! section averaged along y
                      ELSE
                         netcdf_data(i) = init_model%origin_y &
                                     + cos_rot_angle * ( section(i,2) + shift_y ) * dy
                      ENDIF
                   ENDDO

                   nc_stat = NF90_PUT_VAR( id_set_xz(av), id_var_nutm_xz(k,av),&
                                           netcdf_data, start = (/ 1 /),   &
                                           count = (/ ns /) )
                   CALL netcdf_handle_error( 'netcdf_define_header', 556 )

                ENDDO
                DEALLOCATE( netcdf_data )

             ELSE
!
!--             2D in case of rotation
                ALLOCATE( netcdf_data_2d(0:nx,1:ns) )
                cos_rot_angle = COS( rotation_angle * pi / 180.0_wp )
                sin_rot_angle = SIN( rotation_angle * pi / 180.0_wp )

                DO  k = 0, 2
!
!--                Scalar grid points
                   IF ( k == 0 )  THEN
                      shift_x = 0.5 ; shift_y = 0.5
!
!--                u grid points
                   ELSEIF ( k == 1 )  THEN
                      shift_x = 0.0 ; shift_y = 0.5
!
!--                v grid points
                   ELSEIF ( k == 2 )  THEN
                      shift_x = 0.5 ; shift_y = 0.0
                   ENDIF

                   DO  j = 1, ns
                      IF( section(j,2) == -1 )  THEN
                         netcdf_data_2d(:,j) = -1.0_wp  ! section averaged along y
                      ELSE
                         DO  i = 0, nx
                            netcdf_data_2d(i,j) = init_model%origin_x                 &
                                    + cos_rot_angle * ( i + shift_x ) * dx            &
                                    + sin_rot_angle * ( section(j,2) + shift_y ) * dy
                         ENDDO
                      ENDIF
                   ENDDO

                   nc_stat = NF90_PUT_VAR( id_set_xz(av), id_var_eutm_xz(k,av),  &
                                           netcdf_data_2d, start = (/ 1, 1 /),   &
                                           count = (/ nx+1, ns /) )
                   CALL netcdf_handle_error( 'netcdf_define_header', 555 )

                   DO  j = 1, ns
                      IF( section(j,2) == -1 )  THEN
                         netcdf_data_2d(:,j) = -1.0_wp  ! section averaged along y
                      ELSE
                         DO  i = 0, nx
                            netcdf_data_2d(i,j) = init_model%origin_y                 &
                                    - sin_rot_angle * ( i + shift_x ) * dx            &
                                    + cos_rot_angle * ( section(j,2) + shift_y ) * dy
                         ENDDO
                      ENDIF
                   ENDDO

                   nc_stat = NF90_PUT_VAR( id_set_xz(av), id_var_nutm_xz(k,av),  &
                                           netcdf_data_2d, start = (/ 1, 1 /),   &
                                           count = (/ nx+1, ns /) )
                   CALL netcdf_handle_error( 'netcdf_define_header', 556 )

                ENDDO
                DEALLOCATE( netcdf_data_2d )
             ENDIF
!
!--          Write lon and lat data
             ALLOCATE( lat(0:nx,1:ns) )
             ALLOCATE( lon(0:nx,1:ns) )
             cos_rot_angle = COS( rotation_angle * pi / 180.0_wp )
             sin_rot_angle = SIN( rotation_angle * pi / 180.0_wp )

             DO  k = 0, 2
!
!--             Scalar grid points
                IF ( k == 0 )  THEN
                   shift_x = 0.5 ; shift_y = 0.5
!
!--             u grid points
                ELSEIF ( k == 1 )  THEN
                   shift_x = 0.0 ; shift_y = 0.5
!
!--             v grid points
                ELSEIF ( k == 2 )  THEN
                   shift_x = 0.5 ; shift_y = 0.0
                ENDIF

                DO  j = 1, ns
                   IF( section(j,2) == -1 )  THEN
                      lat(:,j) = -90.0_wp  ! section averaged along y
                      lon(:,j) = -180.0_wp  ! section averaged along y
                   ELSE
                      DO  i = 0, nx
                         eutm = init_model%origin_x                   &
                              + cos_rot_angle * ( i + shift_x ) * dx  &
                              + sin_rot_angle * ( section(j,2) + shift_y ) * dy
                         nutm = init_model%origin_y                   &
                              - sin_rot_angle * ( i + shift_x ) * dx  &
                              + cos_rot_angle * ( section(j,2) + shift_y ) * dy

                         CALL  convert_utm_to_geographic( crs_list,          &
                                                          eutm, nutm,        &
                                                          lon(i,j), lat(i,j) )
                      ENDDO
                   ENDIF
                ENDDO

                nc_stat = NF90_PUT_VAR( id_set_xz(av), id_var_lon_xz(k,av), &
                                     lon, start = (/ 1, 1 /),       &
                                     count = (/ nx+1, ns /) )
                CALL netcdf_handle_error( 'netcdf_define_header', 556 )

                nc_stat = NF90_PUT_VAR( id_set_xz(av), id_var_lat_xz(k,av), &
                                     lat, start = (/ 1, 1 /),       &
                                     count = (/ nx+1, ns /) )
                CALL netcdf_handle_error( 'netcdf_define_header', 556 )
             ENDDO

             DEALLOCATE( lat )
             DEALLOCATE( lon )

          ENDIF


       CASE ( 'xz_ext' )

!
!--       Get the list of variables and compare with the actual run.
!--       First var_list_old has to be reset, since GET_ATT does not assign
!--       trailing blanks.
          var_list_old = ' '
          nc_stat = NF90_GET_ATT( id_set_xz(av), NF90_GLOBAL, 'VAR_LIST', &
                                  var_list_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 168 )

          var_list = ';'
          i = 1
          DO WHILE ( do2d(av,i)(1:1) /= ' ' )
             IF ( INDEX( do2d(av,i), 'xz' ) /= 0 )  THEN
                var_list = TRIM( var_list ) // TRIM( do2d(av,i) ) // ';'
             ENDIF
             i = i + 1
          ENDDO

          IF ( av == 0 )  THEN
             var = '(xz)'
          ELSE
             var = '(xz_av)'
          ENDIF

          IF ( TRIM( var_list ) /= TRIM( var_list_old ) )  THEN
             message_string = 'netCDF file for cross-sections ' //           &
                              TRIM( var ) // ' from previous run found,' //  &
                              '&but this file cannot be extended due to' //  &
                              ' variable mismatch.' //                       &
                              '&New file is created instead.'
             CALL message( 'define_netcdf_header', 'PA0249', 0, 1, 0, 6, 0 )
             extend = .FALSE.
             RETURN
          ENDIF

!
!--       Calculate the number of current sections
          ns = 1
          DO WHILE ( section(ns,2) /= -9999  .AND.  ns <= 100 )
             ns = ns + 1
          ENDDO
          ns = ns - 1

!
!--       Get and compare the number of vertical cross sections
          nc_stat = NF90_INQ_VARID( id_set_xz(av), 'y_xz', id_var_y_xz(av) )
          CALL netcdf_handle_error( 'netcdf_define_header', 169 )

          nc_stat = NF90_INQUIRE_VARIABLE( id_set_xz(av), id_var_y_xz(av), &
                                           dimids = id_dim_y_xz_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 170 )
          id_dim_y_xz(av) = id_dim_y_xz_old(1)

          nc_stat = NF90_INQUIRE_DIMENSION( id_set_xz(av), id_dim_y_xz(av), &
                                            len = ns_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 171 )

          IF ( ns /= ns_old )  THEN
             message_string = 'netCDF file for cross-sections ' //          &
                              TRIM( var ) // ' from previous run found,' // &
                              '&but this file cannot be extended due to' // &
                              ' mismatch in number of' //                   &
                              ' cross sections.' //                         &
                              '&New file is created instead.'
             CALL message( 'define_netcdf_header', 'PA0250', 0, 1, 0, 6, 0 )
             extend = .FALSE.
             RETURN
          ENDIF

!
!--       Get and compare the heights of the cross sections
          ALLOCATE( netcdf_data(1:ns_old) )

          nc_stat = NF90_GET_VAR( id_set_xz(av), id_var_y_xz(av), netcdf_data )
          CALL netcdf_handle_error( 'netcdf_define_header', 172 )

          DO  i = 1, ns
             IF ( section(i,2) /= -1 )  THEN
                IF ( ( ( section(i,2) + 0.5 ) * dy ) /= netcdf_data(i) )  THEN
                   message_string = 'netCDF file for cross-sections ' //       &
                               TRIM( var ) // ' from previous run found,' //   &
                               ' but this file cannot be extended' //          &
                               ' due to mismatch in cross' //                  &
                               ' section levels.' //                           &
                               ' New file is created instead.'
                   CALL message( 'define_netcdf_header', 'PA0251',             &
                                                                 0, 1, 0, 6, 0 )
                   extend = .FALSE.
                   RETURN
                ENDIF
             ELSE
                IF ( -1.0_wp /= netcdf_data(i) )  THEN
                   message_string = 'netCDF file for cross-sections ' //       &
                               TRIM( var ) // ' from previous run found,' //   &
                               ' but this file cannot be extended' //          &
                               ' due to mismatch in cross' //                  &
                               ' section levels.' //                           &
                               ' New file is created instead.'
                   CALL message( 'define_netcdf_header', 'PA0251',             &
                                                                 0, 1, 0, 6, 0 )
                   extend = .FALSE.
                   RETURN
                ENDIF
             ENDIF
          ENDDO

          DEALLOCATE( netcdf_data )

!
!--       Get the id of the time coordinate (unlimited coordinate) and its
!--       last index on the file. The next time level is do2d..count+1.
!--       The current time must be larger than the last output time
!--       on the file.
          nc_stat = NF90_INQ_VARID( id_set_xz(av), 'time', id_var_time_xz(av) )
          CALL netcdf_handle_error( 'netcdf_define_header', 173 )

          nc_stat = NF90_INQUIRE_VARIABLE( id_set_xz(av), id_var_time_xz(av), &
                                           dimids = id_dim_time_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 174 )
          id_dim_time_xz(av) = id_dim_time_old(1)

          nc_stat = NF90_INQUIRE_DIMENSION( id_set_xz(av), id_dim_time_xz(av), &
                                            len = ntime_count )
          CALL netcdf_handle_error( 'netcdf_define_header', 175 )

!
!--       For non-parallel output use the last output time level of the netcdf
!--       file because the time dimension is unlimited. In case of parallel
!--       output the variable ntime_count could get the value of 9*10E36 because
!--       the time dimension is limited.
          IF ( netcdf_data_format < 5 ) do2d_xz_time_count(av) = ntime_count

          nc_stat = NF90_GET_VAR( id_set_xz(av), id_var_time_xz(av),           &
                                  last_time_coordinate,                        &
                                  start = (/ do2d_xz_time_count(av) /),        &
                                  count = (/ 1 /) )
          CALL netcdf_handle_error( 'netcdf_define_header', 176 )

          IF ( last_time_coordinate(1) >= simulated_time )  THEN
             message_string = 'netCDF file for cross sections ' //             &
                              TRIM( var ) // ' from previous run found,' //    &
                              '&but this file cannot be extended becaus' //    &
                              'e the current output time' //                   &
                              '&is less or equal than the last output t' //    &
                              'ime on this file.' //                           &
                              '&New file is created instead.'
             CALL message( 'define_netcdf_header', 'PA0252', 0, 1, 0, 6, 0 )
             do2d_xz_time_count(av) = 0
             extend = .FALSE.
             RETURN
          ENDIF

          IF ( netcdf_data_format > 4 )  THEN
!
!--          Check if the needed number of output time levels is increased
!--          compared to the number of time levels in the existing file.
             IF ( ntdim_2d_xz(av) > ntime_count )  THEN
                message_string = 'netCDF file for cross sections ' // &
                                 TRIM( var ) // ' from previous run found,' // &
                                 '&but this file cannot be extended becaus' // &
                                 'e the number of output time levels has b' // &
                                 'een increased compared to the previous s' // &
                                 'imulation.' //                               &
                                 '&New file is created instead.'
                CALL message( 'define_netcdf_header', 'PA0390', 0, 1, 0, 6, 0 )
                do2d_xz_time_count(av) = 0
                extend = .FALSE.
!
!--             Recalculate the needed time levels for the new file.
                IF ( av == 0 )  THEN
                   ntdim_2d_xz(0) = CEILING(                            &
                           ( end_time - MAX( skip_time_do2d_xz,         &
                                             simulated_time_at_begin )  &
                           ) / dt_do2d_xz )
                   IF ( do2d_at_begin )  ntdim_2d_xz(0) = ntdim_2d_xz(0) + 1
                ELSE
                   ntdim_2d_xz(1) = CEILING(                            &
                           ( end_time - MAX( skip_time_data_output_av,  &
                                             simulated_time_at_begin )  &
                           ) / dt_data_output_av )
                ENDIF
                RETURN
             ENDIF
          ENDIF

!
!--       Dataset seems to be extendable.
!--       Now get the variable ids.
          i = 1
          DO WHILE ( do2d(av,i)(1:1) /= ' ' )
             IF ( INDEX( do2d(av,i), 'xz' ) /= 0 )  THEN
                nc_stat = NF90_INQ_VARID( id_set_xz(av), do2d(av,i), &
                                          id_var_do2d(av,i) )
                CALL netcdf_handle_error( 'netcdf_define_header', 177 )
#if defined( __netcdf4_parallel )
!
!--             Set independent io operations for parallel io. Collective io
!--             is only allowed in case of a 1d-decomposition along x, because
!--             otherwise, not all PEs have output data.
                IF ( netcdf_data_format > 4 )  THEN
                   IF ( npey == 1 )  THEN
                      nc_stat = NF90_VAR_PAR_ACCESS( id_set_xz(av),     &
                                                     id_var_do2d(av,i), &
                                                     NF90_COLLECTIVE )
                   ELSE
!
!--                   Test simulations showed that the output of cross sections
!--                   by all PEs in data_output_2d using NF90_COLLECTIVE is
!--                   faster than the output by the first row of PEs in
!--                   x-direction using NF90_INDEPENDENT.
                      nc_stat = NF90_VAR_PAR_ACCESS( id_set_xz(av),     &
                                                     id_var_do2d(av,i), &
                                                     NF90_COLLECTIVE )
!                      nc_stat = NF90_VAR_PAR_ACCESS( id_set_xz(av),     &
!                                                     id_var_do2d(av,i), &
!                                                     NF90_INDEPENDENT )
                   ENDIF
                   CALL netcdf_handle_error( 'netcdf_define_header', 455 )
                ENDIF
#endif
             ENDIF
             i = i + 1
          ENDDO

!
!--       Update the title attribute on file
!--       In order to avoid 'data mode' errors if updated attributes are larger
!--       than their original size, NF90_PUT_ATT is called in 'define mode'
!--       enclosed by NF90_REDEF and NF90_ENDDEF calls. This implies a possible
!--       performance loss due to data copying; an alternative strategy would be
!--       to ensure equal attribute size in a job chain. Maybe revise later.
          IF ( av == 0 )  THEN
             time_average_text = ' '
          ELSE
             WRITE (time_average_text, '('', '',F7.1,'' s average'')') &
                                                            averaging_interval
          ENDIF
          nc_stat = NF90_REDEF( id_set_xz(av) )
          CALL netcdf_handle_error( 'netcdf_define_header', 433 )
          nc_stat = NF90_PUT_ATT( id_set_xz(av), NF90_GLOBAL, 'title',         &
                                  TRIM( run_description_header ) //            &
                                  TRIM( time_average_text ) )
          CALL netcdf_handle_error( 'netcdf_define_header', 178 )
          nc_stat = NF90_ENDDEF( id_set_xz(av) )
          CALL netcdf_handle_error( 'netcdf_define_header', 434 )
          message_string = 'netCDF file for cross-sections ' //                &
                            TRIM( var ) // ' from previous run found.' //      &
                           '&This file will be extended.'
          CALL message( 'define_netcdf_header', 'PA0253', 0, 0, 0, 6, 0 )


       CASE ( 'yz_new' )

!
!--       Define some global attributes of the dataset
          IF ( av == 0 )  THEN
             CALL netcdf_create_global_atts( id_set_yz(av), 'yz', TRIM( run_description_header ), 179 )
             time_average_text = ' '
          ELSE
             CALL netcdf_create_global_atts( id_set_yz(av), 'yz_av', TRIM( run_description_header ), 179 )
             WRITE ( time_average_text,'(F7.1,'' s avg'')' )  averaging_interval
             nc_stat = NF90_PUT_ATT( id_set_yz(av), NF90_GLOBAL, 'time_avg',   &
                                     TRIM( time_average_text ) )
             CALL netcdf_handle_error( 'netcdf_define_header', 180 )
          ENDIF

!
!--       Define time coordinate for yz sections.
!--       For parallel output the time dimensions has to be limited, otherwise
!--       the performance drops significantly.
          IF ( netcdf_data_format < 5 )  THEN
             CALL netcdf_create_dim( id_set_yz(av), 'time', NF90_UNLIMITED,    &
                                     id_dim_time_yz(av), 181 )
          ELSE
             CALL netcdf_create_dim( id_set_yz(av), 'time', ntdim_2d_yz(av),   &
                                     id_dim_time_yz(av), 526 )
          ENDIF

          CALL netcdf_create_var( id_set_yz(av), (/ id_dim_time_yz(av) /),     &
                                  'time', NF90_DOUBLE, id_var_time_yz(av),     &
                                  'seconds', 'time', 182, 183, 000 )
          CALL netcdf_create_att( id_set_yz(av), id_var_time_yz(av), 'standard_name', 'time', 000)
          CALL netcdf_create_att( id_set_yz(av), id_var_time_yz(av), 'axis', 'T', 000)
!
!--       Define the spatial dimensions and coordinates for yz-sections.
!--       First, determine the number of vertical sections to be written.
          IF ( section(1,3) == -9999 )  THEN
             RETURN
          ELSE
             ns = 1
             DO WHILE ( section(ns,3) /= -9999  .AND.  ns <= 100 )
                ns = ns + 1
             ENDDO
             ns = ns - 1
          ENDIF

!
!--       Define x axis (for scalar position)
          CALL netcdf_create_dim( id_set_yz(av), 'x_yz', ns, id_dim_x_yz(av),  &
                                  184 )
          CALL netcdf_create_var( id_set_yz(av), (/ id_dim_x_yz(av) /),        &
                                  'x_yz', NF90_DOUBLE, id_var_x_yz(av),        &
                                  'meters', '', 185, 186, 000 )
          CALL netcdf_create_att( id_set_yz(av), id_var_x_yz(av), 'axis',      &
                                  'X', 000)
!
!--       Define x axis (for u position)
          CALL netcdf_create_dim( id_set_yz(av), 'xu_yz', ns,                  &
                                  id_dim_xu_yz(av), 377 )
          CALL netcdf_create_var( id_set_yz(av), (/ id_dim_xu_yz(av) /),       &
                                  'xu_yz', NF90_DOUBLE, id_var_xu_yz(av),      &
                                  'meters', '', 378, 379, 000 )
          CALL netcdf_create_att( id_set_yz(av), id_var_xu_yz(av), 'axis',     &
                                  'X', 000)
!
!--       Define a variable to store the layer indices of the vertical cross
!--       sections
          CALL netcdf_create_var( id_set_yz(av), (/ id_dim_x_yz(av) /),        &
                                  'ind_x_yz', NF90_DOUBLE,                     &
                                  id_var_ind_x_yz(av), 'gridpoints', '', 187,  &
                                  188, 000 )
!
!--       Define y-axis (for scalar position)
          CALL netcdf_create_dim( id_set_yz(av), 'y', ny+1, id_dim_y_yz(av),   &
                                  189 )
          CALL netcdf_create_var( id_set_yz(av), (/ id_dim_y_yz(av) /), 'y',   &
                                  NF90_DOUBLE, id_var_y_yz(av), 'meters', '',  &
                                  190, 191, 000 )
          CALL netcdf_create_att( id_set_yz(av), id_var_y_yz(av), 'axis',      &
                                  'Y', 000)
!
!--       Define y-axis (for v position)
          CALL netcdf_create_dim( id_set_yz(av), 'yv', ny+1, id_dim_yv_yz(av), &
                                  380 )
          CALL netcdf_create_var( id_set_yz(av), (/ id_dim_yv_yz(av) /), 'yv', &
                                  NF90_DOUBLE, id_var_yv_yz(av), 'meters', '', &
                                  381, 382, 000 )
          CALL netcdf_create_att( id_set_yz(av), id_var_yv_yz(av), 'axis',     &
                                  'Y', 000)
!
!--       Define the two z-axes (zu and zw)
          CALL netcdf_create_dim( id_set_yz(av), 'zu', nz+2, id_dim_zu_yz(av), &
                                  192 )
          CALL netcdf_create_var( id_set_yz(av), (/ id_dim_zu_yz(av) /), 'zu', &
                                  NF90_DOUBLE, id_var_zu_yz(av), 'meters', '', &
                                  193, 194, 000 )
          CALL netcdf_create_att( id_set_yz(av), id_var_zu_yz(av), 'axis',     &
                                  'Z', 000)

          CALL netcdf_create_dim( id_set_yz(av), 'zw', nz+2, id_dim_zw_yz(av), &
                                  195 )
          CALL netcdf_create_var( id_set_yz(av), (/ id_dim_zw_yz(av) /), 'zw', &
                                  NF90_DOUBLE, id_var_zw_yz(av), 'meters', '', &
                                  196, 197, 000 )
          CALL netcdf_create_att( id_set_yz(av), id_var_zw_yz(av), 'axis',     &
                                  'Z', 000)
!
!--       Define UTM and geographic coordinates
          CALL define_geo_coordinates( id_set_yz(av),         &
                  (/ id_dim_x_yz(av), id_dim_xu_yz(av) /),    &
                  (/ id_dim_y_yz(av), id_dim_yv_yz(av) /),    &
                  id_var_eutm_yz(:,av), id_var_nutm_yz(:,av), &
                  id_var_lat_yz(:,av), id_var_lon_yz(:,av)    )
!
!--       Define coordinate-reference system
          CALL netcdf_create_crs( id_set_yz(av), 000 )

          IF ( land_surface )  THEN

             CALL netcdf_create_dim( id_set_yz(av), 'zs', nzs,                 &
                                     id_dim_zs_yz(av), 545 )
             CALL netcdf_create_var( id_set_yz(av), (/ id_dim_zs_yz(av) /),    &
                                     'zs', NF90_DOUBLE, id_var_zs_yz(av),      &
                                     'meters', '', 546, 547, 000 )
             CALL netcdf_create_att( id_set_yz(av), id_var_zs_yz(av), 'axis',  &
                                     'Z', 000)

          ENDIF

!
!--       Define the variables
          var_list = ';'
          i = 1

          DO WHILE ( do2d(av,i)(1:1) /= ' ' )

             IF ( INDEX( do2d(av,i), 'yz' ) /= 0 )  THEN

!
!--             Check for the grid
                found = .FALSE.
                SELECT CASE ( do2d(av,i) )
!
!--                Most variables are defined on the zu grid
                   CASE ( 'e_yz', 'nc_yz', 'nr_yz', 'p_yz', 'pc_yz',           &
                          'pr_yz','prr_yz', 'q_yz', 'qc_yz', 'ql_yz',          &
                          'ql_c_yz', 'ql_v_yz', 'ql_vp_yz', 'qr_yz', 'qv_yz',  &
                          's_yz',                                              &
                          'theta_yz', 'thetal_yz', 'thetav_yz', 'ti_yz' )

                      grid_x = 'x'
                      grid_y = 'y'
                      grid_z = 'zu'
!
!--                u grid
                   CASE ( 'u_yz' )

                      grid_x = 'xu'
                      grid_y = 'y'
                      grid_z = 'zu'
!
!--                v grid
                   CASE ( 'v_yz' )

                      grid_x = 'x'
                      grid_y = 'yv'
                      grid_z = 'zu'
!
!--                w grid
                   CASE ( 'w_yz' )

                      grid_x = 'x'
                      grid_y = 'y'
                      grid_z = 'zw'


                   CASE DEFAULT
!
!--                   Check for land surface quantities
                      IF ( land_surface )  THEN
                         CALL lsm_define_netcdf_grid( do2d(av,i), found,       &
                                                      grid_x, grid_y, grid_z )
                      ENDIF

                      IF ( .NOT. found )  THEN
                         CALL tcm_define_netcdf_grid( do2d(av,i), found,       &
                                                      grid_x, grid_y, grid_z )
                      ENDIF

!
!--                   Check for ocean quantities
                      IF ( .NOT. found  .AND.  ocean_mode )  THEN
                         CALL ocean_define_netcdf_grid( do2d(av,i), found,     &
                                                       grid_x, grid_y, grid_z )
                      ENDIF
!
!--                   Check for radiation quantities
                      IF ( .NOT. found  .AND.  radiation )  THEN
                         CALL radiation_define_netcdf_grid( do2d(av,i), found, &
                                                            grid_x, grid_y,    &
                                                            grid_z )
                      ENDIF
!
!--                   Check for SALSA quantities
                      IF ( .NOT. found  .AND.  salsa )  THEN
                         CALL salsa_define_netcdf_grid( do2d(av,i), found,     &
                                                        grid_x, grid_y, grid_z )
                      ENDIF
!
!--                   Check for gust module quantities
                      IF ( .NOT. found  .AND.  gust_module_enabled )  THEN
                         CALL gust_define_netcdf_grid( do2d(av,i), found,      &
                                                       grid_x, grid_y, grid_z )
                      ENDIF

!
!--                   Check for chemistry quantities
                      IF ( .NOT. found  .AND.  air_chemistry )  THEN
                         CALL chem_define_netcdf_grid( do2d(av,i), found,      &
                                                       grid_x, grid_y,         &
                                                       grid_z )
                      ENDIF

                      IF ( .NOT. found )                                       &
                         CALL doq_define_netcdf_grid(                          &
                                                    do2d(av,i), found, grid_x, &
                                                    grid_y, grid_z           )
!
!--                   Check for user-defined quantities
                      IF ( .NOT. found  .AND.  user_module_enabled )  THEN
                         CALL user_define_netcdf_grid( do2d(av,i), found,      &
                                                       grid_x, grid_y, grid_z )
                      ENDIF

                      IF ( .NOT. found )  THEN
                         WRITE ( message_string, * ) 'no grid defined for',    &
                                                ' variable ', TRIM( do2d(av,i) )
                         CALL message( 'define_netcdf_header', 'PA0244',       &
                                       0, 1, 0, 6, 0 )
                      ENDIF

                END SELECT

!
!--             Select the respective dimension ids
                IF ( grid_x == 'x' )  THEN
                   id_x = id_dim_x_yz(av)
                ELSEIF ( grid_x == 'xu' )  THEN
                   id_x = id_dim_xu_yz(av)
                ENDIF

                IF ( grid_y == 'y' )  THEN
                   id_y = id_dim_y_yz(av)
                ELSEIF ( grid_y == 'yv' )  THEN
                   id_y = id_dim_yv_yz(av)
                ENDIF

                IF ( grid_z == 'zu' )  THEN
                   id_z = id_dim_zu_yz(av)
                ELSEIF ( grid_z == 'zw' )  THEN
                   id_z = id_dim_zw_yz(av)
                ELSEIF ( grid_z == 'zs' )  THEN
                   id_z = id_dim_zs_yz(av)
                ENDIF

!
!--             Define the grid
                CALL netcdf_create_var( id_set_yz(av),  (/ id_x, id_y, id_z,   &
                                        id_dim_time_yz(av) /), do2d(av,i),     &
                                        nc_precision(3), id_var_do2d(av,i),    &
                                        TRIM( do2d_unit(av,i) ), do2d(av,i),   &
                                        198, 199, 356, .TRUE. )

#if defined( __netcdf4_parallel )
                IF ( netcdf_data_format > 4 )  THEN
!
!--                Set no fill for every variable to increase performance.
                   nc_stat = NF90_DEF_VAR_FILL( id_set_yz(av),     &
                                                id_var_do2d(av,i), &
                                                NF90_NOFILL, 0 )
                   CALL netcdf_handle_error( 'netcdf_define_header', 535 )
!
!--                Set independent io operations for parallel io. Collective io
!--                is only allowed in case of a 1d-decomposition along y,
!--                because otherwise, not all PEs have output data.
                   IF ( npex == 1 )  THEN
                      nc_stat = NF90_VAR_PAR_ACCESS( id_set_yz(av),     &
                                                     id_var_do2d(av,i), &
                                                     NF90_COLLECTIVE )
                   ELSE
!
!--                   Test simulations showed that the output of cross sections
!--                   by all PEs in data_output_2d using NF90_COLLECTIVE is
!--                   faster than the output by the first row of PEs in
!--                   y-direction using NF90_INDEPENDENT.
                      nc_stat = NF90_VAR_PAR_ACCESS( id_set_yz(av),     &
                                                     id_var_do2d(av,i), &
                                                     NF90_COLLECTIVE )
!                      nc_stat = NF90_VAR_PAR_ACCESS( id_set_yz(av),     &
!                                                     id_var_do2d(av,i), &
!                                                     NF90_INDEPENDENT )
                   ENDIF
                   CALL netcdf_handle_error( 'netcdf_define_header', 450 )
                ENDIF
#endif
                var_list = TRIM( var_list ) // TRIM( do2d(av,i) ) // ';'

             ENDIF

             i = i + 1

          ENDDO

!
!--       No arrays to output. Close the netcdf file and return.
          IF ( i == 1 )  RETURN

!
!--       Write the list of variables as global attribute (this is used by
!--       restart runs and by combine_plot_fields)
          nc_stat = NF90_PUT_ATT( id_set_yz(av), NF90_GLOBAL, 'VAR_LIST', &
                                  var_list )
          CALL netcdf_handle_error( 'netcdf_define_header', 200 )

!
!--       Set general no fill, otherwise the performance drops significantly for
!--       parallel output.
          nc_stat = NF90_SET_FILL( id_set_yz(av), NF90_NOFILL, oldmode )
          CALL netcdf_handle_error( 'netcdf_define_header', 531 )

!
!--       Leave netCDF define mode
          nc_stat = NF90_ENDDEF( id_set_yz(av) )
          CALL netcdf_handle_error( 'netcdf_define_header', 201 )

!
!--       These data are only written by PE0 for parallel output to increase
!--       the performance.
          IF ( myid == 0  .OR.  netcdf_data_format < 5 )  THEN

!
!--          Write axis data: x_yz, y, zu, zw
             ALLOCATE( netcdf_data(1:ns) )

!
!--          Write x_yz data (shifted by +dx/2)
             DO  i = 1, ns
                IF( section(i,3) == -1 )  THEN
                   netcdf_data(i) = -1.0_wp  ! section averaged along x
                ELSE
                   netcdf_data(i) = ( section(i,3) + 0.5_wp ) * dx
                ENDIF
             ENDDO
             nc_stat = NF90_PUT_VAR( id_set_yz(av), id_var_x_yz(av), &
                                     netcdf_data, start = (/ 1 /),   &
                                     count = (/ ns /) )
             CALL netcdf_handle_error( 'netcdf_define_header', 202 )

!
!--          Write x_yz data (xu grid)
             DO  i = 1, ns
                IF( section(i,3) == -1 )  THEN
                   netcdf_data(i) = -1.0_wp  ! section averaged along x
                ELSE
                   netcdf_data(i) = section(i,3) * dx
                ENDIF
             ENDDO
             nc_stat = NF90_PUT_VAR( id_set_yz(av), id_var_xu_yz(av), &
                                     netcdf_data, start = (/ 1 /),    &
                                     count = (/ ns /) )
             CALL netcdf_handle_error( 'netcdf_define_header', 383 )

!
!--          Write gridpoint number data
             netcdf_data(1:ns) = section(1:ns,3)
             nc_stat = NF90_PUT_VAR( id_set_yz(av), id_var_ind_x_yz(av), &
                                     netcdf_data, start = (/ 1 /),       &
                                     count = (/ ns /) )
             CALL netcdf_handle_error( 'netcdf_define_header', 203 )

             DEALLOCATE( netcdf_data )

!
!--          Write data for y (shifted by +dy/2) and yv axis
             ALLOCATE( netcdf_data(0:ny) )

             DO  j = 0, ny
                netcdf_data(j) = ( j + 0.5_wp ) * dy
             ENDDO

             nc_stat = NF90_PUT_VAR( id_set_yz(av), id_var_y_yz(av), &
                                     netcdf_data, start = (/ 1 /),   &
                                     count = (/ ny+1 /) )
             CALL netcdf_handle_error( 'netcdf_define_header', 204 )

             DO  j = 0, ny
                netcdf_data(j) = j * dy
             ENDDO

             nc_stat = NF90_PUT_VAR( id_set_yz(av), id_var_yv_yz(av), &
                                     netcdf_data, start = (/ 1 /),    &
                                     count = (/ ny+1 /) )
             CALL netcdf_handle_error( 'netcdf_define_header', 384 )

             DEALLOCATE( netcdf_data )

!
!--          Write zu and zw data (vertical axes)
             ALLOCATE( netcdf_data(0:nz+1) )

             netcdf_data(0:nz+1) = zu(nzb:nzt+1)
             nc_stat = NF90_PUT_VAR( id_set_yz(av), id_var_zu_yz(av), &
                                     netcdf_data, start = (/ 1 /),    &
                                     count = (/ nz+2 /) )
             CALL netcdf_handle_error( 'netcdf_define_header', 205 )

             netcdf_data(0:nz+1) = zw(nzb:nzt+1)
             nc_stat = NF90_PUT_VAR( id_set_yz(av), id_var_zw_yz(av), &
                                     netcdf_data, start = (/ 1 /),    &
                                     count = (/ nz+2 /) )
             CALL netcdf_handle_error( 'netcdf_define_header', 206 )

             DEALLOCATE( netcdf_data )
!
!--          Write UTM coordinates
             IF ( rotation_angle == 0.0_wp )  THEN
!
!--             1D in case of no rotation
                cos_rot_angle = COS( rotation_angle * pi / 180.0_wp )
!
!--             x coordinates
                ALLOCATE( netcdf_data(1:ns) )
                DO  k = 0, 2
!
!--                Scalar grid points
                   IF ( k == 0 )  THEN
                      shift_x = 0.5
!
!--                u grid points
                   ELSEIF ( k == 1 )  THEN
                      shift_x = 0.0
!
!--                v grid points
                   ELSEIF ( k == 2 )  THEN
                      shift_x = 0.5
                   ENDIF

                   DO  i = 1, ns
                      IF( section(i,3) == -1 )  THEN
                         netcdf_data(i) = -1.0_wp  ! section averaged along x
                      ELSE
                         netcdf_data(i) = init_model%origin_x &
                                     + cos_rot_angle * ( section(i,3) + shift_x ) * dx
                      ENDIF
                   ENDDO

                   nc_stat = NF90_PUT_VAR( id_set_yz(av), id_var_eutm_yz(k,av),&
                                           netcdf_data, start = (/ 1 /),   &
                                           count = (/ ns /) )
                   CALL netcdf_handle_error( 'netcdf_define_header', 555 )

                ENDDO
                DEALLOCATE( netcdf_data )
!
!--             y coordinates
                ALLOCATE( netcdf_data(0:ny) )
                DO  k = 0, 2
!
!--                Scalar grid points
                   IF ( k == 0 )  THEN
                      shift_y = 0.5
!
!--                u grid points
                   ELSEIF ( k == 1 )  THEN
                      shift_y = 0.5
!
!--                v grid points
                   ELSEIF ( k == 2 )  THEN
                      shift_y = 0.0
                   ENDIF

                   DO  i = 0, ny
                     netcdf_data(i) = init_model%origin_y                      &
                                    + cos_rot_angle * ( i + shift_y ) * dy
                   ENDDO

                   nc_stat = NF90_PUT_VAR( id_set_yz(av), id_var_nutm_yz(k,av),&
                                           netcdf_data, start = (/ 1 /),   &
                                           count = (/ ny+1 /) )
                   CALL netcdf_handle_error( 'netcdf_define_header', 556 )

                ENDDO
                DEALLOCATE( netcdf_data )

             ELSE
!
!--             2D in case of rotation
                ALLOCATE( netcdf_data_2d(1:ns,0:ny) )
                cos_rot_angle = COS( rotation_angle * pi / 180.0_wp )
                sin_rot_angle = SIN( rotation_angle * pi / 180.0_wp )

                DO  k = 0, 2
!
!--                Scalar grid points
                   IF ( k == 0 )  THEN
                      shift_x = 0.5 ; shift_y = 0.5
!
!--                u grid points
                   ELSEIF ( k == 1 )  THEN
                      shift_x = 0.0 ; shift_y = 0.5
!
!--                v grid points
                   ELSEIF ( k == 2 )  THEN
                      shift_x = 0.5 ; shift_y = 0.0
                   ENDIF

                   DO  j = 0, ny
                      DO  i = 1, ns
                         IF( section(i,3) == -1 )  THEN
                            netcdf_data_2d(i,:) = -1.0_wp !section averaged along x
                         ELSE
                            netcdf_data_2d(i,j) = init_model%origin_x                 &
                                    + cos_rot_angle * ( section(i,3) + shift_x ) * dx &
                                    + sin_rot_angle * ( j + shift_y ) * dy
                         ENDIF
                      ENDDO
                   ENDDO

                   nc_stat = NF90_PUT_VAR( id_set_yz(av), id_var_eutm_yz(k,av),  &
                                           netcdf_data_2d, start = (/ 1, 1 /),   &
                                           count = (/ ns, ny+1 /) )
                   CALL netcdf_handle_error( 'netcdf_define_header', 555 )

                   DO  j = 0, ny
                      DO  i = 1, ns
                         IF( section(i,3) == -1 )  THEN
                            netcdf_data_2d(i,:) = -1.0_wp !section averaged along x
                         ELSE
                            netcdf_data_2d(i,j) = init_model%origin_y                 &
                                    - sin_rot_angle * ( section(i,3) + shift_x ) * dx &
                                    + cos_rot_angle * ( j + shift_y ) * dy
                         ENDIF
                      ENDDO
                   ENDDO

                   nc_stat = NF90_PUT_VAR( id_set_yz(av), id_var_nutm_yz(k,av),  &
                                           netcdf_data_2d, start = (/ 1, 1 /),   &
                                           count = (/ ns, ny+1 /) )
                   CALL netcdf_handle_error( 'netcdf_define_header', 556 )

                ENDDO
                DEALLOCATE( netcdf_data_2d )
             ENDIF
!
!--          Write lon and lat data
             ALLOCATE( lat(1:ns,0:ny) )
             ALLOCATE( lon(1:ns,0:ny) )
             cos_rot_angle = COS( rotation_angle * pi / 180.0_wp )
             sin_rot_angle = SIN( rotation_angle * pi / 180.0_wp )

             DO  k = 0, 2
!
!--             Scalar grid points
                IF ( k == 0 )  THEN
                   shift_x = 0.5 ; shift_y = 0.5
!
!--             u grid points
                ELSEIF ( k == 1 )  THEN
                   shift_x = 0.0 ; shift_y = 0.5
!
!--             v grid points
                ELSEIF ( k == 2 )  THEN
                   shift_x = 0.5 ; shift_y = 0.0
                ENDIF

                DO  j = 0, ny
                   DO  i = 1, ns
                      IF( section(i,3) == -1 )  THEN
                         lat(i,:) = -90.0_wp   ! section averaged along x
                         lon(i,:) = -180.0_wp  ! section averaged along x
                      ELSE
                         eutm = init_model%origin_x                              &
                              + cos_rot_angle * ( section(i,3) + shift_x ) * dx  &
                              + sin_rot_angle * ( j + shift_y ) * dy
                         nutm = init_model%origin_y                              &
                              - sin_rot_angle * ( section(i,3) + shift_x ) * dx  &
                              + cos_rot_angle * ( j + shift_y ) * dy

                         CALL  convert_utm_to_geographic( crs_list,          &
                                                          eutm, nutm,        &
                                                          lon(i,j), lat(i,j) )
                      ENDIF
                   ENDDO
                ENDDO

                nc_stat = NF90_PUT_VAR( id_set_yz(av), id_var_lon_yz(k,av), &
                                     lon, start = (/ 1, 1 /),       &
                                     count = (/ ns, ny+1 /) )
                CALL netcdf_handle_error( 'netcdf_define_header', 556 )

                nc_stat = NF90_PUT_VAR( id_set_yz(av), id_var_lat_yz(k,av), &
                                     lat, start = (/ 1, 1 /),       &
                                     count = (/ ns, ny+1 /) )
                CALL netcdf_handle_error( 'netcdf_define_header', 556 )
             ENDDO

             DEALLOCATE( lat )
             DEALLOCATE( lon )

          ENDIF


       CASE ( 'yz_ext' )

!
!--       Get the list of variables and compare with the actual run.
!--       First var_list_old has to be reset, since GET_ATT does not assign
!--       trailing blanks.
          var_list_old = ' '
          nc_stat = NF90_GET_ATT( id_set_yz(av), NF90_GLOBAL, 'VAR_LIST', &
                                  var_list_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 207 )

          var_list = ';'
          i = 1
          DO WHILE ( do2d(av,i)(1:1) /= ' ' )
             IF ( INDEX( do2d(av,i), 'yz' ) /= 0 )  THEN
                var_list = TRIM( var_list ) // TRIM( do2d(av,i) ) // ';'
             ENDIF
             i = i + 1
          ENDDO

          IF ( av == 0 )  THEN
             var = '(yz)'
          ELSE
             var = '(yz_av)'
          ENDIF

          IF ( TRIM( var_list ) /= TRIM( var_list_old ) )  THEN
             message_string = 'netCDF file for cross-sections ' //           &
                              TRIM( var ) // ' from previous run found,' //  &
                              '&but this file cannot be extended due to' //  &
                              ' variable mismatch.' //                       &
                              '&New file is created instead.'
             CALL message( 'define_netcdf_header', 'PA0249', 0, 1, 0, 6, 0 )
             extend = .FALSE.
             RETURN
          ENDIF

!
!--       Calculate the number of current sections
          ns = 1
          DO WHILE ( section(ns,3) /= -9999  .AND.  ns <= 100 )
             ns = ns + 1
          ENDDO
          ns = ns - 1

!
!--       Get and compare the number of vertical cross sections
          nc_stat = NF90_INQ_VARID( id_set_yz(av), 'x_yz', id_var_x_yz(av) )
          CALL netcdf_handle_error( 'netcdf_define_header', 208 )

          nc_stat = NF90_INQUIRE_VARIABLE( id_set_yz(av), id_var_x_yz(av), &
                                           dimids = id_dim_x_yz_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 209 )
          id_dim_x_yz(av) = id_dim_x_yz_old(1)

          nc_stat = NF90_INQUIRE_DIMENSION( id_set_yz(av), id_dim_x_yz(av), &
                                            len = ns_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 210 )

          IF ( ns /= ns_old )  THEN
             message_string = 'netCDF file for cross-sections ' //          &
                              TRIM( var ) // ' from previous run found,' // &
                              '&but this file cannot be extended due to' // &
                              ' mismatch in number of' //                   &
                              ' cross sections.' //                         &
                              '&New file is created instead.'
             CALL message( 'define_netcdf_header', 'PA0250', 0, 1, 0, 6, 0 )
             extend = .FALSE.
             RETURN
          ENDIF

!
!--       Get and compare the heights of the cross sections
          ALLOCATE( netcdf_data(1:ns_old) )

          nc_stat = NF90_GET_VAR( id_set_yz(av), id_var_x_yz(av), netcdf_data )
          CALL netcdf_handle_error( 'netcdf_define_header', 211 )

          DO  i = 1, ns
             IF ( section(i,3) /= -1 )  THEN
                IF ( ( ( section(i,3) + 0.5 ) * dx ) /= netcdf_data(i) )  THEN
                   message_string = 'netCDF file for cross-sections ' //       &
                              TRIM( var ) // ' from previous run found,' //    &
                              ' but this file cannot be extended' //           &
                              ' due to mismatch in cross' //                   &
                              ' section levels.' //                            &
                              ' New file is created instead.'
                   CALL message( 'define_netcdf_header', 'PA0251',             &
                                                                 0, 1, 0, 6, 0 )
                   extend = .FALSE.
                   RETURN
                ENDIF
             ELSE
                IF ( -1.0_wp /= netcdf_data(i) )  THEN
                   message_string = 'netCDF file for cross-sections ' //       &
                              TRIM( var ) // ' from previous run found,' //    &
                              ' but this file cannot be extended' //           &
                              ' due to mismatch in cross' //                   &
                              ' section levels.' //                            &
                              ' New file is created instead.'
                   CALL message( 'define_netcdf_header', 'PA0251',             &
                                                                 0, 1, 0, 6, 0 )
                   extend = .FALSE.
                   RETURN
                ENDIF
             ENDIF
          ENDDO

          DEALLOCATE( netcdf_data )

!
!--       Get the id of the time coordinate (unlimited coordinate) and its
!--       last index on the file. The next time level is pl2d..count+1.
!--       The current time must be larger than the last output time
!--       on the file.
          nc_stat = NF90_INQ_VARID( id_set_yz(av), 'time', id_var_time_yz(av) )
          CALL netcdf_handle_error( 'netcdf_define_header', 212 )

          nc_stat = NF90_INQUIRE_VARIABLE( id_set_yz(av), id_var_time_yz(av), &
                                           dimids = id_dim_time_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 213 )
          id_dim_time_yz(av) = id_dim_time_old(1)

          nc_stat = NF90_INQUIRE_DIMENSION( id_set_yz(av), id_dim_time_yz(av), &
                                            len = ntime_count )
          CALL netcdf_handle_error( 'netcdf_define_header', 214 )

!
!--       For non-parallel output use the last output time level of the netcdf
!--       file because the time dimension is unlimited. In case of parallel
!--       output the variable ntime_count could get the value of 9*10E36 because
!--       the time dimension is limited.
          IF ( netcdf_data_format < 5 ) do2d_yz_time_count(av) = ntime_count

          nc_stat = NF90_GET_VAR( id_set_yz(av), id_var_time_yz(av),           &
                                  last_time_coordinate,                        &
                                  start = (/ do2d_yz_time_count(av) /),        &
                                  count = (/ 1 /) )
          CALL netcdf_handle_error( 'netcdf_define_header', 215 )

          IF ( last_time_coordinate(1) >= simulated_time )  THEN
             message_string = 'netCDF file for cross sections ' //             &
                              TRIM( var ) // ' from previous run found,' //    &
                              '&but this file cannot be extended becaus' //    &
                              'e the current output time' //                   &
                              '&is less or equal than the last output t' //    &
                              'ime on this file.' //                           &
                              '&New file is created instead.'
             CALL message( 'define_netcdf_header', 'PA0252', 0, 1, 0, 6, 0 )
             do2d_yz_time_count(av) = 0
             extend = .FALSE.
             RETURN
          ENDIF

          IF ( netcdf_data_format > 4 )  THEN
!
!--          Check if the needed number of output time levels is increased
!--          compared to the number of time levels in the existing file.
             IF ( ntdim_2d_yz(av) > ntime_count )  THEN
                message_string = 'netCDF file for cross sections ' //          &
                                 TRIM( var ) // ' from previous run found,' // &
                                 '&but this file cannot be extended becaus' // &
                                 'e the number of output time levels has b' // &
                                 'een increased compared to the previous s' // &
                                 'imulation.' //                               &
                                 '&New file is created instead.'
                CALL message( 'define_netcdf_header', 'PA0391', 0, 1, 0, 6, 0 )
                do2d_yz_time_count(av) = 0
                extend = .FALSE.
!
!--             Recalculate the needed time levels for the new file.
                IF ( av == 0 )  THEN
                   ntdim_2d_yz(0) = CEILING(                            &
                           ( end_time - MAX( skip_time_do2d_yz,         &
                                             simulated_time_at_begin )  &
                           ) / dt_do2d_yz )
                   IF ( do2d_at_begin )  ntdim_2d_yz(0) = ntdim_2d_yz(0) + 1
                ELSE
                   ntdim_2d_yz(1) = CEILING(                            &
                           ( end_time - MAX( skip_time_data_output_av,  &
                                             simulated_time_at_begin )  &
                           ) / dt_data_output_av )
                ENDIF
                RETURN
             ENDIF
          ENDIF

!
!--       Dataset seems to be extendable.
!--       Now get the variable ids.
          i = 1
          DO WHILE ( do2d(av,i)(1:1) /= ' ' )
             IF ( INDEX( do2d(av,i), 'yz' ) /= 0 )  THEN
                nc_stat = NF90_INQ_VARID( id_set_yz(av), do2d(av,i), &
                                          id_var_do2d(av,i) )
                CALL netcdf_handle_error( 'netcdf_define_header', 216 )
#if defined( __netcdf4_parallel )
!
!--             Set independent io operations for parallel io. Collective io
!--             is only allowed in case of a 1d-decomposition along y, because
!--             otherwise, not all PEs have output data.
                IF ( netcdf_data_format > 4 )  THEN
                   IF ( npex == 1 )  THEN
                      nc_stat = NF90_VAR_PAR_ACCESS( id_set_yz(av),     &
                                                     id_var_do2d(av,i), &
                                                     NF90_COLLECTIVE )
                   ELSE
!
!--                   Test simulations showed that the output of cross sections
!--                   by all PEs in data_output_2d using NF90_COLLECTIVE is
!--                   faster than the output by the first row of PEs in
!--                   y-direction using NF90_INDEPENDENT.
                      nc_stat = NF90_VAR_PAR_ACCESS( id_set_yz(av),     &
                                                     id_var_do2d(av,i), &
                                                     NF90_COLLECTIVE )
!                      nc_stat = NF90_VAR_PAR_ACCESS( id_set_yz(av),     &
!                                                     id_var_do2d(av,i), &
!                                                     NF90_INDEPENDENT )
                   ENDIF
                   CALL netcdf_handle_error( 'netcdf_define_header', 450 )
                ENDIF
#endif
             ENDIF
             i = i + 1
          ENDDO

!
!--       Update the title attribute on file
!--       In order to avoid 'data mode' errors if updated attributes are larger
!--       than their original size, NF90_PUT_ATT is called in 'define mode'
!--       enclosed by NF90_REDEF and NF90_ENDDEF calls. This implies a possible
!--       performance loss due to data copying; an alternative strategy would be
!--       to ensure equal attribute size in a job chain. Maybe revise later.
          IF ( av == 0 )  THEN
             time_average_text = ' '
          ELSE
             WRITE (time_average_text, '('', '',F7.1,'' s average'')') &
                                                            averaging_interval
          ENDIF
          nc_stat = NF90_REDEF( id_set_yz(av) )
          CALL netcdf_handle_error( 'netcdf_define_header', 435 )
          nc_stat = NF90_PUT_ATT( id_set_yz(av), NF90_GLOBAL, 'title',         &
                                  TRIM( run_description_header ) //            &
                                  TRIM( time_average_text ) )
          CALL netcdf_handle_error( 'netcdf_define_header', 217 )
          nc_stat = NF90_ENDDEF( id_set_yz(av) )
          CALL netcdf_handle_error( 'netcdf_define_header', 436 )
          message_string = 'netCDF file for cross-sections ' //                &
                            TRIM( var ) // ' from previous run found.' //      &
                           '&This file will be extended.'
          CALL message( 'define_netcdf_header', 'PA0253', 0, 0, 0, 6, 0 )


       CASE ( 'pr_new' )

!
!--       Define some global attributes of the dataset

          IF ( averaging_interval_pr /= 0.0_wp )  THEN
             CALL netcdf_create_global_atts( id_set_pr, 'podsprav', TRIM( run_description_header ), 451 )
             WRITE ( time_average_text,'(F7.1,'' s avg'')' ) averaging_interval_pr
             nc_stat = NF90_PUT_ATT( id_set_pr, NF90_GLOBAL, 'time_avg',       &
                                     TRIM( time_average_text ) )
          ELSE
             CALL netcdf_create_global_atts( id_set_pr, 'podspr', TRIM( run_description_header ), 451 )
          ENDIF
          CALL netcdf_handle_error( 'netcdf_define_header', 219 )
!
!--       Write number of columns and rows of coordinate systems to be plotted
!--       on one page to the netcdf header.
!--       This information can be used by palmplot.
          nc_stat = NF90_PUT_ATT( id_set_pr, NF90_GLOBAL,                     &
                                  'no_rows',                                  &
                                  profile_rows )
          CALL netcdf_handle_error( 'netcdf_define_header', 519 )

          nc_stat = NF90_PUT_ATT( id_set_pr, NF90_GLOBAL,                     &
                                  'no_columns',                               &
                                  profile_columns )
          CALL netcdf_handle_error( 'netcdf_define_header', 520 )


          cross_profiles_adj  = ADJUSTL( cross_profiles )
          cross_profiles_numb = 999999
          cross_profiles_char = ''

!
!--       Each profile defined in cross_profiles is written to an array
!--       (cross_profiles_char). The number of the respective coordinate
!--       system is assigned in a second array (cross_profiles_numb).
          k = 1

          DO  i = 1, crmax

             IF ( TRIM( cross_profiles_adj(i) ) == ' ' )  EXIT
             delim_old = 0

             DO   j = 1, crmax
                delim = INDEX( cross_profiles_adj(i)(delim_old+1:), ' ' )
                IF ( delim == 1 )  EXIT
                kk = MIN( crmax, k )
                cross_profiles_char(kk) = cross_profiles_adj(i)(delim_old+1: &
                                                              delim_old+delim-1)
                cross_profiles_numb(kk) = i
                k = k + 1
                cross_profiles_maxi  = i
                delim_old = delim_old + delim
             ENDDO

          ENDDO

          cross_profiles_count = MIN( crmax, k-1 )
!
!--       Check if all profiles defined in cross_profiles are defined in
!--       data_output_pr. If not, they will be skipped.
          DO  i = 1, cross_profiles_count
             DO  j = 1, dopr_n

                IF ( TRIM(cross_profiles_char(i)) == TRIM(data_output_pr(j)) ) &
                THEN
                   EXIT
                ENDIF

                IF ( j == dopr_n )  THEN
                   cross_profiles_numb(i) = 999999
                ENDIF

             ENDDO
          ENDDO

          DO i = 1, crmax
             IF ( cross_profiles_numb(i) == 999999 ) THEN
                DO j = i + 1, crmax
                   IF ( cross_profiles_numb(j) /= 999999 ) THEN
                      cross_profiles_char(i) = cross_profiles_char(j)
                      cross_profiles_numb(i) = cross_profiles_numb(j)
                      cross_profiles_numb(j) = 999999
                      EXIT
                   ENDIF
                ENDDO
             ENDIF
          ENDDO

          DO i = 1, crmax-1
             IF ( cross_profiles_numb(i + 1) == 999999 ) THEN
                cross_profiles_count = i
                EXIT
             ENDIF
          ENDDO
!
!--       Check if all profiles defined in data_output_pr are defined in
!--       cross_profiles. If not, they will be added to cross_profiles.
          DO  i = 1, dopr_n
             DO  j = 1, cross_profiles_count

                IF ( TRIM(cross_profiles_char(j)) == TRIM(data_output_pr(i)))  &
                THEN
                   EXIT
                ENDIF

                IF (( j == cross_profiles_count ) .AND.                        &
                    ( cross_profiles_count <= crmax - 1))  THEN
                   cross_profiles_count = cross_profiles_count + 1
                   cross_profiles_maxi  = cross_profiles_maxi  + 1
                   cross_profiles_char(MIN( crmax, cross_profiles_count )) =   &
                                                      TRIM( data_output_pr(i) )
                   cross_profiles_numb(MIN( crmax, cross_profiles_count )) =   &
                                                      cross_profiles_maxi
                ENDIF

             ENDDO
          ENDDO

          IF ( cross_profiles_count >= crmax )  THEN
             message_string = 'It is not allowed to arrange more than '        &
                              // '100 profiles with & cross_profiles. Apart'   &
                              // ' from that, all profiles are saved & to '    &
                              // 'the netCDF file.'
             CALL message( 'define_netcdf_header', 'PA0354', 0, 0, 0, 6, 0 )
          ENDIF

!
!--       Writing cross_profiles to netcdf header. This information can be
!--       used by palmplot. Each profile is separated by ",", each cross is
!--       separated by ";".
          char_cross_profiles = ';'
          id_last = 1
          cross_profiles_count = MIN( cross_profiles_count, crmax )

          DO  i = 1, cross_profiles_count

             IF ( cross_profiles_numb(i) /= 999999 )  THEN
                IF ( TRIM( char_cross_profiles ) == ';' )  THEN
                   char_cross_profiles = TRIM( char_cross_profiles ) // &
                                         TRIM( cross_profiles_char(i) )
                ELSEIF ( id_last == cross_profiles_numb(i) )  THEN
                   char_cross_profiles = TRIM( char_cross_profiles ) // &
                                         ',' // TRIM( cross_profiles_char(i) )
                ELSE
                   char_cross_profiles = TRIM( char_cross_profiles ) // &
                                         ';' // TRIM( cross_profiles_char(i) )
                ENDIF
                id_last = cross_profiles_numb(i)
             ENDIF

          ENDDO

          char_cross_profiles = TRIM( char_cross_profiles ) // ';'

          nc_stat = NF90_PUT_ATT( id_set_pr, NF90_GLOBAL, 'cross_profiles',   &
                                  TRIM( char_cross_profiles ) )
          CALL netcdf_handle_error( 'netcdf_define_header', 521 )

!
!--       Define time coordinate for profiles (unlimited dimension)
          CALL netcdf_create_dim( id_set_pr, 'time', NF90_UNLIMITED,           &
                                  id_dim_time_pr, 220 )
          CALL netcdf_create_var( id_set_pr, (/ id_dim_time_pr /), 'time',     &
                                  NF90_DOUBLE, id_var_time_pr, 'seconds', 'time',  &
                                  221, 222, 000 )
          CALL netcdf_create_att( id_set_pr, id_var_time_pr, 'standard_name', 'time', 000)
          CALL netcdf_create_att( id_set_pr, id_var_time_pr, 'axis', 'T', 000)
!
!--       Define the variables
          var_list = ';'
          DO  i = 1, dopr_n

             IF ( statistic_regions == 0 )  THEN

!
!--             Define the z-axes (each variable gets its own z-axis)
                CALL netcdf_create_dim( id_set_pr,                             &
                                        'z' // TRIM( data_output_pr(i) ),      &
                                        nzt+2-nzb, id_dim_z_pr(i,0), 223 )
                CALL netcdf_create_var( id_set_pr, (/ id_dim_z_pr(i,0) /),     &
                                        'z' // TRIM( data_output_pr(i) ),      &
                                       NF90_DOUBLE, id_var_z_pr(i,0),          &
                                       'meters', '', 224, 225, 000 )
                CALL netcdf_create_att( id_set_pr, id_var_z_pr(i,0), 'axis',   &
                                        'Z', 000)
!
!--             Define the variable
                CALL netcdf_create_var( id_set_pr, (/ id_dim_z_pr(i,0),        &
                                        id_dim_time_pr /), data_output_pr(i),  &
                                        nc_precision(5), id_var_dopr(i,0),     &
                                        TRIM( dopr_unit(i) ),                  &
                                        TRIM( data_output_pr(i) ), 226, 227,   &
                                        228 )

                var_list = TRIM( var_list ) // TRIM( data_output_pr(i) ) //  ';'

             ELSE
!
!--             If statistic regions are defined, add suffix _SR+#SR to the
!--             names
                DO  j = 0, statistic_regions
                   WRITE ( suffix, '(''_'',I2.2)' )  j

!
!--                Define the z-axes (each variable gets it own z-axis)
                   CALL netcdf_create_dim( id_set_pr, 'z' //                   &
                                           TRIM(data_output_pr(i)) // suffix,  &
                                           nzt+2-nzb, id_dim_z_pr(i,j), 229 )
                   CALL netcdf_create_var( id_set_pr, (/ id_dim_z_pr(i,j) /),  &
                                           'z' // TRIM(data_output_pr(i)) //   &
                                           suffix, NF90_DOUBLE,                &
                                           id_var_z_pr(i,j), 'meters', '',     &
                                           230, 231, 000 )
                   CALL netcdf_create_att( id_set_pr, id_var_z_pr(i,j), 'axis',&
                                           'Z', 000)
!
!--                Define the variable
                   CALL netcdf_create_var( id_set_pr, (/ id_dim_z_pr(i,j),     &
                                           id_dim_time_pr /),                  &
                                           TRIM(data_output_pr(i)) // suffix,  &
                                           nc_precision(5), id_var_dopr(i,j),  &
                                           TRIM( dopr_unit(i) ),               &
                                           TRIM( data_output_pr(i) ) //        &
                                           ' SR ', 232, 233, 234 )

                   var_list = TRIM( var_list ) // TRIM( data_output_pr(i) ) // &
                              suffix // ';'

                ENDDO

             ENDIF

          ENDDO

!
!--       Write the list of variables as global attribute (this is used by
!--       restart runs)
          nc_stat = NF90_PUT_ATT( id_set_pr, NF90_GLOBAL, 'VAR_LIST', var_list )
          CALL netcdf_handle_error( 'netcdf_define_header', 235 )

!
!--       Define normalization variables (as time series)
          DO  i = 1, dopr_norm_num

             CALL netcdf_create_var( id_set_pr, (/ id_dim_time_pr /),          &
                                     'NORM_' // TRIM( dopr_norm_names(i) ),    &
                                     nc_precision(5), id_var_norm_dopr(i),     &
                                     '', TRIM( dopr_norm_longnames(i) ), 236,  &
                                     000, 237 )

          ENDDO

!
!--       Leave netCDF define mode
          nc_stat = NF90_ENDDEF( id_set_pr )
          CALL netcdf_handle_error( 'netcdf_define_header', 238 )

!
!--       Write z-axes data
          DO  i = 1, dopr_n
             DO  j = 0, statistic_regions

                nc_stat = NF90_PUT_VAR( id_set_pr, id_var_z_pr(i,j),      &
                                        hom(nzb:nzt+1,2,dopr_index(i),0), &
                                        start = (/ 1 /),                  &
                                        count = (/ nzt-nzb+2 /) )
                CALL netcdf_handle_error( 'netcdf_define_header', 239 )

             ENDDO
          ENDDO


       CASE ( 'pr_ext' )

!
!--       Get the list of variables and compare with the actual run.
!--       First var_list_old has to be reset, since GET_ATT does not assign
!--       trailing blanks.
          var_list_old = ' '
          nc_stat = NF90_GET_ATT( id_set_pr, NF90_GLOBAL, 'VAR_LIST', &
                                  var_list_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 240 )

          var_list = ';'
          DO  i = 1, dopr_n

             IF ( statistic_regions == 0 )  THEN
                var_list = TRIM( var_list ) // TRIM( data_output_pr(i) ) // ';'
             ELSE
                DO  j = 0, statistic_regions
                   WRITE ( suffix, '(''_'',I2.2)' )  j
                   var_list = TRIM( var_list ) // TRIM( data_output_pr(i) ) // &
                              suffix // ';'
                ENDDO
             ENDIF

          ENDDO

          IF ( TRIM( var_list ) /= TRIM( var_list_old ) )  THEN
             message_string = 'netCDF file for vertical profiles ' //          &
                              'from previous run found,' //                    &
                              '&but this file cannot be extended due to' //    &
                              ' variable mismatch.' //                         &
                              '&New file is created instead.'
             CALL message( 'define_netcdf_header', 'PA0254', 0, 1, 0, 6, 0 )
             extend = .FALSE.
             RETURN
          ENDIF

!
!--       Get the id of the time coordinate (unlimited coordinate) and its
!--       last index on the file. The next time level is dopr..count+1.
!--       The current time must be larger than the last output time
!--       on the file.
          nc_stat = NF90_INQ_VARID( id_set_pr, 'time', id_var_time_pr )
          CALL netcdf_handle_error( 'netcdf_define_header', 241 )

          nc_stat = NF90_INQUIRE_VARIABLE( id_set_pr, id_var_time_pr, &
                                           dimids = id_dim_time_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 242 )
          id_dim_time_pr = id_dim_time_old(1)

          nc_stat = NF90_INQUIRE_DIMENSION( id_set_pr, id_dim_time_pr, &
                                            len = dopr_time_count )
          CALL netcdf_handle_error( 'netcdf_define_header', 243 )

          nc_stat = NF90_GET_VAR( id_set_pr, id_var_time_pr,        &
                                  last_time_coordinate,             &
                                  start = (/ dopr_time_count /), &
                                  count = (/ 1 /) )
          CALL netcdf_handle_error( 'netcdf_define_header', 244 )

          IF ( last_time_coordinate(1) >= simulated_time )  THEN
             message_string = 'netCDF file for vertical profiles ' //          &
                              'from previous run found,' //                    &
                              '&but this file cannot be extended becaus' //    &
                              'e the current output time' //                   &
                              '&is less or equal than the last output t' //    &
                              'ime on this file.' //                           &
                              '&New file is created instead.'
             CALL message( 'define_netcdf_header', 'PA0255', 0, 1, 0, 6, 0 )
             dopr_time_count = 0
             extend = .FALSE.
             RETURN
          ENDIF

!
!--       Dataset seems to be extendable.
!--       Now get the variable ids.
          i = 1
          DO  i = 1, dopr_n

             IF ( statistic_regions == 0 )  THEN
                nc_stat = NF90_INQ_VARID( id_set_pr, data_output_pr(i),        &
                                          id_var_dopr(i,0) )
                CALL netcdf_handle_error( 'netcdf_define_header', 245 )
             ELSE
                DO  j = 0, statistic_regions
                   WRITE ( suffix, '(''_'',I2.2)' )  j
                   netcdf_var_name = TRIM( data_output_pr(i) ) // suffix
                   nc_stat = NF90_INQ_VARID( id_set_pr, netcdf_var_name,       &
                                             id_var_dopr(i,j) )
                   CALL netcdf_handle_error( 'netcdf_define_header', 246 )
                ENDDO
             ENDIF

          ENDDO

!
!--       Get ids of the normalization variables
          DO  i = 1, dopr_norm_num
             nc_stat = NF90_INQ_VARID( id_set_pr,                             &
                                       'NORM_' // TRIM( dopr_norm_names(i) ), &
                                       id_var_norm_dopr(i) )
             CALL netcdf_handle_error( 'netcdf_define_header', 247 )
          ENDDO

!
!--       Update the title attribute on file
!--       In order to avoid 'data mode' errors if updated attributes are larger
!--       than their original size, NF90_PUT_ATT is called in 'define mode'
!--       enclosed by NF90_REDEF and NF90_ENDDEF calls. This implies a possible
!--       performance loss due to data copying; an alternative strategy would be
!--       to ensure equal attribute size in a job chain. Maybe revise later.
          IF ( averaging_interval_pr == 0.0_wp )  THEN
             time_average_text = ' '
          ELSE
             WRITE (time_average_text, '('', '',F7.1,'' s average'')') &
                                                            averaging_interval_pr
          ENDIF
          nc_stat = NF90_REDEF( id_set_pr )
          CALL netcdf_handle_error( 'netcdf_define_header', 437 )
          nc_stat = NF90_PUT_ATT( id_set_pr, NF90_GLOBAL, 'title',             &
                                  TRIM( run_description_header ) //            &
                                  TRIM( time_average_text ) )
          CALL netcdf_handle_error( 'netcdf_define_header', 248 )

          nc_stat = NF90_ENDDEF( id_set_pr )
          CALL netcdf_handle_error( 'netcdf_define_header', 438 )
          message_string = 'netCDF file for vertical profiles ' //             &
                           'from previous run found.' //                       &
                           '&This file will be extended.'
          CALL message( 'define_netcdf_header', 'PA0256', 0, 0, 0, 6, 0 )


       CASE ( 'ts_new' )

!
!--       Define some global attributes of the dataset
          CALL netcdf_create_global_atts( id_set_ts, 'podsts', TRIM(run_description_header), 329 )

          ! nc_stat = NF90_PUT_ATT( id_set_ts, NF90_GLOBAL, 'title',             &
          !                         TRIM( run_description_header ) )
          ! CALL netcdf_handle_error( 'netcdf_define_header', 249 )

!
!--       Define time coordinate for time series (unlimited dimension)
          CALL netcdf_create_dim( id_set_ts, 'time', NF90_UNLIMITED,           &
                                  id_dim_time_ts, 250 )
          CALL netcdf_create_var( id_set_ts, (/ id_dim_time_ts /), 'time',     &
                                  NF90_DOUBLE, id_var_time_ts, 'seconds', 'time',  &
                                  251, 252, 000 )
          CALL netcdf_create_att( id_set_ts, id_var_time_ts, 'standard_name', 'time', 000)
          CALL netcdf_create_att( id_set_ts, id_var_time_ts, 'axis', 'T', 000)
!
!--       Define the variables
          var_list = ';'
          DO  i = 1, dots_num

             IF ( statistic_regions == 0 )  THEN

                CALL netcdf_create_var( id_set_ts, (/ id_dim_time_ts /),       &
                                        dots_label(i), nc_precision(6),        &
                                        id_var_dots(i,0),                      &
                                        TRIM( dots_unit(i) ),                  &
                                        TRIM( dots_label(i) ), 253, 254, 255 )

                var_list = TRIM( var_list ) // TRIM( dots_label(i) ) // ';'

             ELSE
!
!--             If statistic regions are defined, add suffix _SR+#SR to the
!--             names
                DO  j = 0, statistic_regions
                   WRITE ( suffix, '(''_'',I2.2)' )  j

                   CALL netcdf_create_var( id_set_ts, (/ id_dim_time_ts /),    &
                                           TRIM( dots_label(i) ) // suffix,    &
                                           nc_precision(6), id_var_dots(i,j),  &
                                           TRIM( dots_unit(i) ),               &
                                           TRIM( dots_label(i) ) // ' SR ' //  &
                                           suffix(2:2), 256, 257, 347)

                   var_list = TRIM( var_list ) // TRIM( dots_label(i) ) // &
                              suffix // ';'

                ENDDO

             ENDIF

          ENDDO

!
!--       Write the list of variables as global attribute (this is used by
!--       restart runs)
          nc_stat = NF90_PUT_ATT( id_set_ts, NF90_GLOBAL, 'VAR_LIST', var_list )
          CALL netcdf_handle_error( 'netcdf_define_header', 258 )

!
!--       Leave netCDF define mode
          nc_stat = NF90_ENDDEF( id_set_ts )
          CALL netcdf_handle_error( 'netcdf_define_header', 259 )


       CASE ( 'ts_ext' )

!
!--       Get the list of variables and compare with the actual run.
!--       First var_list_old has to be reset, since GET_ATT does not assign
!--       trailing blanks.
          var_list_old = ' '
          nc_stat = NF90_GET_ATT( id_set_ts, NF90_GLOBAL, 'VAR_LIST', &
                                  var_list_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 260 )

          var_list = ';'
          i = 1
          DO  i = 1, dots_num

             IF ( statistic_regions == 0 )  THEN
                var_list = TRIM( var_list ) // TRIM( dots_label(i) ) // ';'
             ELSE
                DO  j = 0, statistic_regions
                   WRITE ( suffix, '(''_'',I2.2)' )  j
                   var_list = TRIM( var_list ) // TRIM( dots_label(i) ) //     &
                              suffix // ';'
                ENDDO
             ENDIF

          ENDDO

          IF ( TRIM( var_list ) /= TRIM( var_list_old ) )  THEN
             message_string = 'netCDF file for time series ' //                &
                              'from previous run found,' //                    &
                              '&but this file cannot be extended due to' //    &
                              ' variable mismatch.' //                         &
                              '&New file is created instead.'
             CALL message( 'define_netcdf_header', 'PA0257', 0, 1, 0, 6, 0 )
             extend = .FALSE.
             RETURN
          ENDIF

!
!--       Get the id of the time coordinate (unlimited coordinate) and its
!--       last index on the file. The next time level is dots..count+1.
!--       The current time must be larger than the last output time
!--       on the file.
          nc_stat = NF90_INQ_VARID( id_set_ts, 'time', id_var_time_ts )
          CALL netcdf_handle_error( 'netcdf_define_header', 261 )

          nc_stat = NF90_INQUIRE_VARIABLE( id_set_ts, id_var_time_ts,          &
                                           dimids = id_dim_time_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 262 )
          id_dim_time_ts = id_dim_time_old(1)

          nc_stat = NF90_INQUIRE_DIMENSION( id_set_ts, id_dim_time_ts,         &
                                            len = dots_time_count )
          CALL netcdf_handle_error( 'netcdf_define_header', 263 )

          nc_stat = NF90_GET_VAR( id_set_ts, id_var_time_ts,                   &
                                  last_time_coordinate,                        &
                                  start = (/ dots_time_count /),               &
                                  count = (/ 1 /) )
          CALL netcdf_handle_error( 'netcdf_define_header', 264 )

          IF ( last_time_coordinate(1) >= simulated_time )  THEN
             message_string = 'netCDF file for time series ' //                &
                              'from previous run found,' //                    &
                              '&but this file cannot be extended becaus' //    &
                              'e the current output time' //                   &
                              '&is less or equal than the last output t' //    &
                              'ime on this file.' //                           &
                              '&New file is created instead.'
             CALL message( 'define_netcdf_header', 'PA0258', 0, 1, 0, 6, 0 )
             dots_time_count = 0
             extend = .FALSE.
             RETURN
          ENDIF

!
!--       Dataset seems to be extendable.
!--       Now get the variable ids
          i = 1
          DO  i = 1, dots_num

             IF ( statistic_regions == 0 )  THEN
                nc_stat = NF90_INQ_VARID( id_set_ts, dots_label(i), &
                                          id_var_dots(i,0) )
                CALL netcdf_handle_error( 'netcdf_define_header', 265 )
             ELSE
                DO  j = 0, statistic_regions
                   WRITE ( suffix, '(''_'',I2.2)' )  j
                   netcdf_var_name = TRIM( dots_label(i) ) // suffix
                   nc_stat = NF90_INQ_VARID( id_set_ts, netcdf_var_name, &
                                             id_var_dots(i,j) )
                   CALL netcdf_handle_error( 'netcdf_define_header', 266 )
                ENDDO
             ENDIF

          ENDDO

!
!--       Update the title attribute on file
!--       In order to avoid 'data mode' errors if updated attributes are larger
!--       than their original size, NF90_PUT_ATT is called in 'define mode'
!--       enclosed by NF90_REDEF and NF90_ENDDEF calls. This implies a possible
!--       performance loss due to data copying; an alternative strategy would be
!--       to ensure equal attribute size in a job chain. Maybe revise later.
          nc_stat = NF90_REDEF( id_set_ts )
          CALL netcdf_handle_error( 'netcdf_define_header', 439 )
          nc_stat = NF90_PUT_ATT( id_set_ts, NF90_GLOBAL, 'title',             &
                                  TRIM( run_description_header ) )
          CALL netcdf_handle_error( 'netcdf_define_header', 267 )
          nc_stat = NF90_ENDDEF( id_set_ts )
          CALL netcdf_handle_error( 'netcdf_define_header', 440 )
          message_string = 'netCDF file for time series ' //                   &
                           'from previous run found.' //                       &
                           '&This file will be extended.'
          CALL message( 'define_netcdf_header', 'PA0259', 0, 0, 0, 6, 0 )


       CASE ( 'sp_new' )

!
!--       Define some global attributes of the dataset
          IF ( averaging_interval_sp /= 0.0_wp )  THEN
             WRITE (time_average_text,'('', '',F7.1,'' s average'')')          &
                                                            averaging_interval_sp
             nc_stat = NF90_PUT_ATT( id_set_sp, NF90_GLOBAL, 'title',          &
                                     TRIM( run_description_header ) //         &
                                     TRIM( time_average_text ) )
             CALL netcdf_handle_error( 'netcdf_define_header', 268 )

             WRITE ( time_average_text,'(F7.1,'' s avg'')' )  averaging_interval_sp
             nc_stat = NF90_PUT_ATT( id_set_sp, NF90_GLOBAL, 'time_avg', &
                                     TRIM( time_average_text ) )
          ELSE
             nc_stat = NF90_PUT_ATT( id_set_sp, NF90_GLOBAL, 'title', &
                                     TRIM( run_description_header ) )
          ENDIF
          CALL netcdf_handle_error( 'netcdf_define_header', 269 )

!
!--       Define time coordinate for spectra (unlimited dimension)
          CALL netcdf_create_dim( id_set_sp, 'time', NF90_UNLIMITED,           &
                                  id_dim_time_sp, 270 )
          CALL netcdf_create_var( id_set_sp, (/ id_dim_time_sp /), 'time',     &
                                  NF90_DOUBLE, id_var_time_sp, 'seconds', 'time',  &
                                  271, 272, 000 )
          CALL netcdf_create_att( id_set_sp, id_var_time_sp, 'standard_name', 'time', 000)
          CALL netcdf_create_att( id_set_sp, id_var_time_sp, 'axis', 'T', 000)
!
!--       Define the spatial dimensions and coordinates for spectra.
!--       First, determine the number of vertical levels for which spectra
!--       are to be output.
          ns = 1
          DO WHILE ( comp_spectra_level(ns) /= 999999  .AND.  ns <= 100 )
             ns = ns + 1
          ENDDO
          ns = ns - 1

!
!--       Define vertical coordinate grid (zu grid)
          CALL netcdf_create_dim( id_set_sp, 'zu_sp', ns, id_dim_zu_sp, 273 )
          CALL netcdf_create_var( id_set_sp, (/ id_dim_zu_sp /), 'zu_sp',      &
                                  NF90_DOUBLE, id_var_zu_sp, 'meters', '',     &
                                  274, 275, 000 )
          CALL netcdf_create_att( id_set_sp, id_var_zu_sp, 'axis', 'Z', 000)
!
!--       Define vertical coordinate grid (zw grid)
          CALL netcdf_create_dim( id_set_sp, 'zw_sp', ns, id_dim_zw_sp, 276 )
          CALL netcdf_create_var( id_set_sp, (/ id_dim_zw_sp /), 'zw_sp',      &
                                  NF90_DOUBLE, id_var_zw_sp, 'meters', '',     &
                                  277, 278, 000 )
          CALL netcdf_create_att( id_set_sp, id_var_zw_sp, 'axis', 'Z', 000)
!
!--       Define x-axis
          CALL netcdf_create_dim( id_set_sp, 'k_x', nx/2, id_dim_x_sp, 279 )
          CALL netcdf_create_var( id_set_sp, (/ id_dim_x_sp /), 'k_x',         &
                                  NF90_DOUBLE, id_var_x_sp, 'm-1', '', 280,    &
                                  281, 000 )
          CALL netcdf_create_att( id_set_sp, id_var_x_sp, 'axis', 'X', 000)
!
!--       Define y-axis
          CALL netcdf_create_dim( id_set_sp, 'k_y', ny/2, id_dim_y_sp, 282 )
          CALL netcdf_create_var( id_set_sp, (/ id_dim_y_sp /), 'k_y',         &
                                  NF90_DOUBLE, id_var_y_sp, 'm-1', '', 283,    &
                                  284, 000 )
          CALL netcdf_create_att( id_set_sp, id_var_y_sp, 'axis', 'Y', 000)
!
!--       Define the variables
          var_list = ';'
          i = 1
          DO WHILE ( data_output_sp(i) /= ' '  .AND.  i <= 10 )
!
!--          First check for the vertical grid
             found = .FALSE.
             SELECT CASE ( data_output_sp(i) )
!
!--             Most variables are defined on the zu levels
                CASE ( 'e', 'nc', 'nr', 'p', 'pc', 'pr', 'prr',   &
                       'q', 'qc', 'ql', 'ql_c', 'ql_v', 'ql_vp', 'qr', 'qv',   &
                       'rho_sea_water', 's', 'sa', &
                       'theta', 'thetal', 'thetav', 'u', 'v' )

                   grid_z = 'zu'
!
!--             zw levels
                CASE ( 'w' )

                   grid_z = 'zw'

                CASE DEFAULT
!
!--                Check for user-defined quantities (found, grid_x and grid_y
!--                are dummies)
                   IF ( user_module_enabled )  THEN
                      CALL user_define_netcdf_grid( data_output_sp(i), found,  &
                                                    grid_x, grid_y, grid_z )
                   ENDIF

             END SELECT

             IF ( INDEX( spectra_direction(i), 'x' ) /= 0 )  THEN

!
!--             Define the variable
                netcdf_var_name = TRIM( data_output_sp(i) ) // '_x'
                IF ( TRIM( grid_z ) == 'zw' )  THEN
                   CALL netcdf_create_var( id_set_sp, (/ id_dim_x_sp,          &
                                           id_dim_zw_sp, id_dim_time_sp /),    &
                                           netcdf_var_name, nc_precision(7),   &
                                           id_var_dospx(i), 'unknown',         &
                                           netcdf_var_name, 285, 286, 287 )
                ELSE
                   CALL netcdf_create_var( id_set_sp, (/ id_dim_x_sp,          &
                                           id_dim_zu_sp, id_dim_time_sp /),    &
                                           netcdf_var_name, nc_precision(7),   &
                                           id_var_dospx(i), 'unknown',         &
                                           netcdf_var_name, 285, 286, 287 )
                ENDIF

                var_list = TRIM( var_list ) // TRIM( netcdf_var_name ) // ';'

             ENDIF

             IF ( INDEX( spectra_direction(i), 'y' ) /= 0 )  THEN

!
!--             Define the variable
                netcdf_var_name = TRIM( data_output_sp(i) ) // '_y'
                IF ( TRIM( grid_z ) == 'zw' )  THEN
                   CALL netcdf_create_var( id_set_sp, (/ id_dim_y_sp,          &
                                           id_dim_zw_sp, id_dim_time_sp /),    &
                                           netcdf_var_name, nc_precision(7),   &
                                           id_var_dospy(i), 'unknown',         &
                                           netcdf_var_name, 288, 289, 290 )
                ELSE
                   CALL netcdf_create_var( id_set_sp, (/ id_dim_y_sp,          &
                                           id_dim_zu_sp, id_dim_time_sp /),    &
                                           netcdf_var_name, nc_precision(7),   &
                                           id_var_dospy(i), 'unknown',         &
                                           netcdf_var_name, 288, 289, 290 )
                ENDIF

                var_list = TRIM( var_list ) // TRIM( netcdf_var_name ) // ';'

             ENDIF

             i = i + 1

          ENDDO

!
!--       Write the list of variables as global attribute (this is used by
!--       restart runs)
          nc_stat = NF90_PUT_ATT( id_set_sp, NF90_GLOBAL, 'VAR_LIST', var_list )
          CALL netcdf_handle_error( 'netcdf_define_header', 291 )

!
!--       Leave netCDF define mode
          nc_stat = NF90_ENDDEF( id_set_sp )
          CALL netcdf_handle_error( 'netcdf_define_header', 292 )

!
!--       Write axis data: zu_sp, zw_sp, k_x, k_y
          ALLOCATE( netcdf_data(1:ns) )

!
!--       Write zu data
          netcdf_data(1:ns) = zu( comp_spectra_level(1:ns) )
          nc_stat = NF90_PUT_VAR( id_set_sp, id_var_zu_sp, netcdf_data, &
                                  start = (/ 1 /), count = (/ ns /) )
          CALL netcdf_handle_error( 'netcdf_define_header', 293 )

!
!--       Write zw data
          netcdf_data(1:ns) = zw( comp_spectra_level(1:ns) )
          nc_stat = NF90_PUT_VAR( id_set_sp, id_var_zw_sp, netcdf_data, &
                                  start = (/ 1 /), count = (/ ns /) )
          CALL netcdf_handle_error( 'netcdf_define_header', 294 )

          DEALLOCATE( netcdf_data )

!
!--       Write data for x and y axis (wavenumbers)
          ALLOCATE( netcdf_data(nx/2) )
          DO  i = 1, nx/2
             netcdf_data(i) = 2.0_wp * pi * i / ( dx * ( nx + 1 ) )
          ENDDO

          nc_stat = NF90_PUT_VAR( id_set_sp, id_var_x_sp, netcdf_data, &
                                  start = (/ 1 /), count = (/ nx/2 /) )
          CALL netcdf_handle_error( 'netcdf_define_header', 295 )

          DEALLOCATE( netcdf_data )

          ALLOCATE( netcdf_data(ny/2) )
          DO  i = 1, ny/2
             netcdf_data(i) = 2.0_wp * pi * i / ( dy * ( ny + 1 ) )
          ENDDO

          nc_stat = NF90_PUT_VAR( id_set_sp, id_var_y_sp, netcdf_data, &
                                  start = (/ 1 /), count = (/ ny/2 /) )
          CALL netcdf_handle_error( 'netcdf_define_header', 296 )

          DEALLOCATE( netcdf_data )


       CASE ( 'sp_ext' )

!
!--       Get the list of variables and compare with the actual run.
!--       First var_list_old has to be reset, since GET_ATT does not assign
!--       trailing blanks.
          var_list_old = ' '
          nc_stat = NF90_GET_ATT( id_set_sp, NF90_GLOBAL, 'VAR_LIST', &
                                  var_list_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 297 )
          var_list = ';'
          i = 1
          DO WHILE ( data_output_sp(i) /= ' '  .AND.  i <= 10 )

             IF ( INDEX( spectra_direction(i), 'x' ) /= 0 )  THEN
                netcdf_var_name = TRIM( data_output_sp(i) ) // '_x'
                var_list = TRIM( var_list ) // TRIM( netcdf_var_name ) // ';'
             ENDIF

             IF ( INDEX( spectra_direction(i), 'y' ) /= 0 )  THEN
                netcdf_var_name = TRIM( data_output_sp(i) ) // '_y'
                var_list = TRIM( var_list ) // TRIM( netcdf_var_name ) // ';'
             ENDIF

             i = i + 1

          ENDDO

          IF ( TRIM( var_list ) /= TRIM( var_list_old ) )  THEN
             message_string = 'netCDF file for spectra  ' //                   &
                              'from previous run found,' //                    &
                              '&but this file cannot be extended due to' //    &
                              ' variable mismatch.' //                         &
                              '&New file is created instead.'
             CALL message( 'define_netcdf_header', 'PA0260', 0, 1, 0, 6, 0 )
             extend = .FALSE.
             RETURN
          ENDIF

!
!--       Determine the number of current vertical levels for which spectra
!--       shall be output
          ns = 1
          DO WHILE ( comp_spectra_level(ns) /= 999999  .AND.  ns <= 100 )
             ns = ns + 1
          ENDDO
          ns = ns - 1

!
!--       Get and compare the number of vertical levels
          nc_stat = NF90_INQ_VARID( id_set_sp, 'zu_sp', id_var_zu_sp )
          CALL netcdf_handle_error( 'netcdf_define_header', 298 )

          nc_stat = NF90_INQUIRE_VARIABLE( id_set_sp, id_var_zu_sp, &
                                           dimids = id_dim_zu_sp_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 299 )
          id_dim_zu_sp = id_dim_zu_sp_old(1)

          nc_stat = NF90_INQUIRE_DIMENSION( id_set_sp, id_dim_zu_sp, &
                                            len = ns_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 300 )

          IF ( ns /= ns_old )  THEN
             message_string = 'netCDF file for spectra ' //                    &
                              ' from previous run found,' //                   &
                              '&but this file cannot be extended due to' //    &
                              ' mismatch in number of' //                      &
                              ' vertical levels.' //                           &
                              '&New file is created instead.'
             CALL message( 'define_netcdf_header', 'PA0261', 0, 1, 0, 6, 0 )
             extend = .FALSE.
             RETURN
          ENDIF

!
!--       Get and compare the heights of the cross sections
          ALLOCATE( netcdf_data(1:ns_old) )

          nc_stat = NF90_GET_VAR( id_set_sp, id_var_zu_sp, netcdf_data )
          CALL netcdf_handle_error( 'netcdf_define_header', 301 )

          DO  i = 1, ns
             IF ( zu(comp_spectra_level(i)) /= netcdf_data(i) )  THEN
                message_string = 'netCDF file for spectra ' //                 &
                                 ' from previous run found,' //                &
                                 '&but this file cannot be extended due to' // &
                                 ' mismatch in heights of' //                  &
                                 ' vertical levels.' //                        &
                                 '&New file is created instead.'
                CALL message( 'define_netcdf_header', 'PA0262', 0, 1, 0, 6, 0 )
                extend = .FALSE.
                RETURN
             ENDIF
          ENDDO

          DEALLOCATE( netcdf_data )

!
!--       Get the id of the time coordinate (unlimited coordinate) and its
!--       last index on the file. The next time level is plsp..count+1.
!--       The current time must be larger than the last output time
!--       on the file.
          nc_stat = NF90_INQ_VARID( id_set_sp, 'time', id_var_time_sp )
          CALL netcdf_handle_error( 'netcdf_define_header', 302 )

          nc_stat = NF90_INQUIRE_VARIABLE( id_set_sp, id_var_time_sp, &
                                           dimids = id_dim_time_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 303 )
          id_dim_time_sp = id_dim_time_old(1)

          nc_stat = NF90_INQUIRE_DIMENSION( id_set_sp, id_dim_time_sp, &
                                            len = dosp_time_count )
          CALL netcdf_handle_error( 'netcdf_define_header', 304 )

          nc_stat = NF90_GET_VAR( id_set_sp, id_var_time_sp,        &
                                  last_time_coordinate,             &
                                  start = (/ dosp_time_count /), &
                                  count = (/ 1 /) )
          CALL netcdf_handle_error( 'netcdf_define_header', 305 )

          IF ( last_time_coordinate(1) >= simulated_time )  THEN
             message_string = 'netCDF file for spectra ' //                    &
                              'from previous run found,' //                    &
                              '&but this file cannot be extended becaus' //    &
                              'e the current output time' //                   &
                              '&is less or equal than the last output t' //    &
                              'ime on this file.' //                           &
                              '&New file is created instead.'
             CALL message( 'define_netcdf_header', 'PA0263', 0, 1, 0, 6, 0 )
             dosp_time_count = 0
             extend = .FALSE.
             RETURN
          ENDIF

!
!--       Dataset seems to be extendable.
!--       Now get the variable ids.
          i = 1
          DO WHILE ( data_output_sp(i) /= ' '  .AND.  i <= 10 )

             IF ( INDEX( spectra_direction(i), 'x' ) /= 0 )  THEN
                netcdf_var_name = TRIM( data_output_sp(i) ) // '_x'
                nc_stat = NF90_INQ_VARID( id_set_sp, netcdf_var_name, &
                                          id_var_dospx(i) )
                CALL netcdf_handle_error( 'netcdf_define_header', 306 )
             ENDIF

             IF ( INDEX( spectra_direction(i), 'y' ) /= 0 )  THEN
                netcdf_var_name = TRIM( data_output_sp(i) ) // '_y'
                nc_stat = NF90_INQ_VARID( id_set_sp, netcdf_var_name, &
                                          id_var_dospy(i) )
                CALL netcdf_handle_error( 'netcdf_define_header', 307 )
             ENDIF

             i = i + 1

          ENDDO

!
!--       Update the title attribute on file
!--       In order to avoid 'data mode' errors if updated attributes are larger
!--       than their original size, NF90_PUT_ATT is called in 'define mode'
!--       enclosed by NF90_REDEF and NF90_ENDDEF calls. This implies a possible
!--       performance loss due to data copying; an alternative strategy would be
!--       to ensure equal attribute size in a job chain. Maybe revise later.
          nc_stat = NF90_REDEF( id_set_sp )
          CALL netcdf_handle_error( 'netcdf_define_header', 441 )
          IF ( averaging_interval_sp /= 0.0_wp )  THEN
             WRITE (time_average_text,'('', '',F7.1,'' s average'')') &
                                                           averaging_interval_sp
             nc_stat = NF90_PUT_ATT( id_set_sp, NF90_GLOBAL, 'title',  &
                                     TRIM( run_description_header ) // &
                                     TRIM( time_average_text ) )
             CALL netcdf_handle_error( 'netcdf_define_header', 308 )

             WRITE ( time_average_text,'(F7.1,'' s avg'')' )  averaging_interval_sp
             nc_stat = NF90_PUT_ATT( id_set_sp, NF90_GLOBAL, 'time_avg',       &
                                     TRIM( time_average_text ) )
          ELSE
             nc_stat = NF90_PUT_ATT( id_set_sp, NF90_GLOBAL, 'title',          &
                                     TRIM( run_description_header ) )
          ENDIF
          CALL netcdf_handle_error( 'netcdf_define_header', 309 )
          nc_stat = NF90_ENDDEF( id_set_sp )
          CALL netcdf_handle_error( 'netcdf_define_header', 442 )
          message_string = 'netCDF file for spectra ' //                       &
                           'from previous run found.' //                       &
                           '&This file will be extended.'
          CALL message( 'define_netcdf_header', 'PA0264', 0, 0, 0, 6, 0 )

!
!--     Currently disabled (DATA_PRT_NETCDF)
!       CASE ( 'pt_new' )

!
!--       Define some global attributes of the dataset
!          nc_stat = NF90_PUT_ATT( id_set_prt, NF90_GLOBAL, 'title',            &
!                                  TRIM( run_description_header ) )
!          CALL netcdf_handle_error( 'netcdf_define_header', 310 )

!
!--       Define time coordinate for particles (unlimited dimension)
!          CALL netcdf_create_dim( id_set_prt, 'time', NF90_UNLIMITED,          &
!                                  id_dim_time_prt, 311 )
!          CALL netcdf_create_var( id_set_prt, (/ id_dim_time_prt /), 'time',   &
!                                  NF90_DOUBLE, id_var_time_prt, 'seconds', '', &
!                                  312, 313, 000 )
!
!--       netCDF4 allows more than one unlimited dimension
!          CALL netcdf_create_dim( id_set_prt, 'particle_number',            &
!                                  NF90_UNLIMITED, id_dim_prtnum, 314 )

!          CALL netcdf_create_var( id_set_prt, (/ id_dim_prtnum /),             &
!                                  'particle_number', NF90_DOUBLE,              &
!                                  id_var_prtnum, 'particle number', '', 315,   &
!                                  316, 000 )
!
!--       Define variable which contains the real number of particles in use
!          CALL netcdf_create_var( id_set_prt, (/ id_dim_time_prt /),           &
!                                  'real_num_of_prt', NF90_DOUBLE,              &
!                                  id_var_rnop_prt, 'particle number', '', 317, &
!                                  318, 000 )
!
!--       Define the variables
!          DO  i = 1, 17

!             CALL netcdf_create_var( id_set_prt, (/ id_dim_prtnum,             &
!                                     id_dim_time_prt /), prt_var_names(i),     &
!                                     nc_precision(8), id_var_prt(i),           &
!                                     TRIM( prt_var_units(i) ),                 &
!                                     TRIM( prt_var_names(i) ), 319, 320, 321 )

!          ENDDO

!
!--       Leave netCDF define mode
!          nc_stat = NF90_ENDDEF( id_set_prt )
!          CALL netcdf_handle_error( 'netcdf_define_header', 322 )

!
!--     Currently disabled (DATA_PRT_NETCDF)
!       CASE ( 'pt_ext' )

!
!--       Get the id of the time coordinate (unlimited coordinate) and its
!--       last index on the file. The next time level is prt..count+1.
!--       The current time must be larger than the last output time
!--       on the file.
!          nc_stat = NF90_INQ_VARID( id_set_prt, 'time', id_var_time_prt )
!          CALL netcdf_handle_error( 'netcdf_define_header', 323 )

!          nc_stat = NF90_INQUIRE_VARIABLE( id_set_prt, id_var_time_prt, &
!                                           dimids = id_dim_time_old )
!          CALL netcdf_handle_error( 'netcdf_define_header', 324 )
!          id_dim_time_prt = id_dim_time_old(1)

!          nc_stat = NF90_INQUIRE_DIMENSION( id_set_prt, id_dim_time_prt, &
!                                            len = prt_time_count )
!          CALL netcdf_handle_error( 'netcdf_define_header', 325 )

!          nc_stat = NF90_GET_VAR( id_set_prt, id_var_time_prt,  &
!                                  last_time_coordinate,         &
!                                  start = (/ prt_time_count /), &
!                                  count = (/ 1 /) )
!          CALL netcdf_handle_error( 'netcdf_define_header', 326 )

!          IF ( last_time_coordinate(1) >= simulated_time )  THEN
!             message_string = 'netCDF file for particles ' //                  &
!                              'from previous run found,' //                    &
!                              '&but this file cannot be extended becaus' //    &
!                              'e the current output time' //                   &
!                              '&is less or equal than the last output t' //    &
!                              'ime on this file.' //                           &
!                              '&New file is created instead.'
!             CALL message( 'define_netcdf_header', 'PA0265', 0, 1, 0, 6, 0 )
!             prt_time_count = 0
!             extend = .FALSE.
!             RETURN
!          ENDIF

!
!--       Dataset seems to be extendable.
!--       Now get the variable ids.
!         nc_stat = NF90_INQ_VARID( id_set_prt, 'real_num_of_prt',             &
!                                   id_var_rnop_prt )
!         CALL netcdf_handle_error( 'netcdf_define_header', 327 )

!          DO  i = 1, 17

!             nc_stat = NF90_INQ_VARID( id_set_prt, prt_var_names(i),           &
!                                       id_var_prt(i) )
!             CALL netcdf_handle_error( 'netcdf_define_header', 328 )

!          ENDDO

!          message_string = 'netCDF file for particles ' //                     &
!                           'from previous run found.' //                       &
!                           '&This file will be extended.'
!          CALL message( 'define_netcdf_header', 'PA0266', 0, 0, 0, 6, 0 )



       CASE ( 'ps_new' )

!
!--       Define some global attributes of the dataset
          nc_stat = NF90_PUT_ATT( id_set_pts, NF90_GLOBAL, 'title',            &
                                  TRIM( run_description_header ) )
          CALL netcdf_handle_error( 'netcdf_define_header', 396 )

!
!--       Define time coordinate for particle time series (unlimited dimension)
          CALL netcdf_create_dim( id_set_pts, 'time', NF90_UNLIMITED,          &
                                  id_dim_time_pts, 397 )
          CALL netcdf_create_var( id_set_pts, (/ id_dim_time_pts /), 'time',   &
                                  NF90_DOUBLE, id_var_time_pts, 'seconds', 'time', &
                                  398, 399, 000 )
          CALL netcdf_create_att( id_set_pts, id_var_time_pts, 'standard_name', 'time', 000)
          CALL netcdf_create_att( id_set_pts, id_var_time_pts, 'axis', 'T', 000)
!
!--       Define the variables. If more than one particle group is defined,
!--       define seperate variables for each group
          var_list = ';'
          DO  i = 1, dopts_num

             DO  j = 0, number_of_particle_groups

                IF ( j == 0 )  THEN
                   suffix = ''
                ELSE
                   WRITE ( suffix, '(''_'',I2.2)' )  j
                ENDIF

                IF ( j == 0 )  THEN
                   CALL netcdf_create_var( id_set_pts, (/ id_dim_time_pts /),  &
                                           TRIM( dopts_label(i) ) // suffix,  &
                                           nc_precision(6), id_var_dopts(i,j), &
                                           TRIM( dopts_unit(i) ),              &
                                           TRIM( dopts_label(i) ), 400, 401,   &
                                           402 )
                ELSE
                   CALL netcdf_create_var( id_set_pts, (/ id_dim_time_pts /),  &
                                           TRIM( dopts_label(i) ) // suffix,  &
                                           nc_precision(6), id_var_dopts(i,j), &
                                           TRIM( dopts_unit(i) ),              &
                                           TRIM( dopts_label(i) ) // ' PG ' // &
                                           suffix(2:3), 400, 401, 402 )
                ENDIF

                var_list = TRIM( var_list ) // TRIM( dopts_label(i) ) // &
                           suffix // ';'

                IF ( number_of_particle_groups == 1 )  EXIT

             ENDDO

          ENDDO

!
!--       Write the list of variables as global attribute (this is used by
!--       restart runs)
          nc_stat = NF90_PUT_ATT( id_set_pts, NF90_GLOBAL, 'VAR_LIST', &
                                  var_list )
          CALL netcdf_handle_error( 'netcdf_define_header', 403 )


!
!--       Leave netCDF define mode
          nc_stat = NF90_ENDDEF( id_set_pts )
          CALL netcdf_handle_error( 'netcdf_define_header', 404 )


       CASE ( 'ps_ext' )

!
!--       Get the list of variables and compare with the actual run.
!--       First var_list_old has to be reset, since GET_ATT does not assign
!--       trailing blanks.
          var_list_old = ' '
          nc_stat = NF90_GET_ATT( id_set_pts, NF90_GLOBAL, 'VAR_LIST', &
                                  var_list_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 405 )

          var_list = ';'
          i = 1
          DO  i = 1, dopts_num

             DO  j = 0, number_of_particle_groups

                IF ( j == 0 )  THEN
                   suffix = ''
                ELSE
                   WRITE ( suffix, '(''_'',I2.2)' )  j
                ENDIF

                var_list = TRIM( var_list ) // TRIM( dopts_label(i) ) // &
                           suffix // ';'

                IF ( number_of_particle_groups == 1 )  EXIT

             ENDDO

          ENDDO

          IF ( TRIM( var_list ) /= TRIM( var_list_old ) )  THEN
             message_string = 'netCDF file for particle time series ' //       &
                              'from previous run found,' //                    &
                              '&but this file cannot be extended due to' //    &
                              ' variable mismatch.' //                         &
                              '&New file is created instead.'
             CALL message( 'define_netcdf_header', 'PA0267', 0, 1, 0, 6, 0 )
             extend = .FALSE.
             RETURN
          ENDIF

!
!--       Get the id of the time coordinate (unlimited coordinate) and its
!--       last index on the file. The next time level is dots..count+1.
!--       The current time must be larger than the last output time
!--       on the file.
          nc_stat = NF90_INQ_VARID( id_set_pts, 'time', id_var_time_pts )
          CALL netcdf_handle_error( 'netcdf_define_header', 406 )

          nc_stat = NF90_INQUIRE_VARIABLE( id_set_pts, id_var_time_pts, &
                                           dimids = id_dim_time_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 407 )
          id_dim_time_pts = id_dim_time_old(1)

          nc_stat = NF90_INQUIRE_DIMENSION( id_set_pts, id_dim_time_pts, &
                                            len = dopts_time_count )
          CALL netcdf_handle_error( 'netcdf_define_header', 408 )

          nc_stat = NF90_GET_VAR( id_set_pts, id_var_time_pts,    &
                                  last_time_coordinate,           &
                                  start = (/ dopts_time_count /), &
                                  count = (/ 1 /) )
          CALL netcdf_handle_error( 'netcdf_define_header', 409 )

          IF ( last_time_coordinate(1) >= simulated_time )  THEN
             message_string = 'netCDF file for particle time series ' //       &
                              'from previous run found,' //                    &
                              '&but this file cannot be extended becaus' //    &
                              'e the current output time' //                   &
                              '&is less or equal than the last output t' //    &
                              'ime on this file.' //                           &
                              '&New file is created instead.'
             CALL message( 'define_netcdf_header', 'PA0268', 0, 1, 0, 6, 0 )
             dopts_time_count = 0
             extend = .FALSE.
             RETURN
          ENDIF

!
!--       Dataset seems to be extendable.
!--       Now get the variable ids
          i = 1
          DO  i = 1, dopts_num

             DO  j = 0, number_of_particle_groups

                IF ( j == 0 )  THEN
                   suffix = ''
                ELSE
                   WRITE ( suffix, '(''_'',I2.2)' )  j
                ENDIF

                netcdf_var_name = TRIM( dopts_label(i) ) // suffix

                nc_stat = NF90_INQ_VARID( id_set_pts, netcdf_var_name, &
                                          id_var_dopts(i,j) )
                CALL netcdf_handle_error( 'netcdf_define_header', 410 )

                IF ( number_of_particle_groups == 1 )  EXIT

             ENDDO

          ENDDO

!
!--       Update the title attribute on file
!--       In order to avoid 'data mode' errors if updated attributes are larger
!--       than their original size, NF90_PUT_ATT is called in 'define mode'
!--       enclosed by NF90_REDEF and NF90_ENDDEF calls. This implies a possible
!--       performance loss due to data copying; an alternative strategy would be
!--       to ensure equal attribute size in a job chain. Maybe revise later.
          nc_stat = NF90_REDEF( id_set_pts )
          CALL netcdf_handle_error( 'netcdf_define_header', 443 )
          nc_stat = NF90_PUT_ATT( id_set_pts, NF90_GLOBAL, 'title',            &
                                  TRIM( run_description_header ) )
          CALL netcdf_handle_error( 'netcdf_define_header', 411 )
          nc_stat = NF90_ENDDEF( id_set_pts )
          CALL netcdf_handle_error( 'netcdf_define_header', 444 )
          message_string = 'netCDF file for particle time series ' //          &
                           'from previous run found.' //                       &
                           '&This file will be extended.'
          CALL message( 'netcdf_define_header', 'PA0269', 0, 0, 0, 6, 0 )

!
!--    Flight data
       CASE ( 'fl_new' )
!
!--       Define some global attributes of the dataset
          nc_stat = NF90_PUT_ATT( id_set_fl, NF90_GLOBAL, 'title',             &
                                  TRIM( run_description_header ) )
          CALL netcdf_handle_error( 'netcdf_define_header', 249 )

!
!--       Define time and location coordinates for flight space-time series
!--       (unlimited dimension)
!--       Error number must still be set appropriately.
          CALL netcdf_create_dim( id_set_fl, 'time', NF90_UNLIMITED,           &
                                  id_dim_time_fl, 250 )
          CALL netcdf_create_var( id_set_fl, (/ id_dim_time_fl /), 'time',     &
                                  NF90_DOUBLE, id_var_time_fl, 'seconds', 'time',  &
                                  251, 252, 000 )
          CALL netcdf_create_att( id_set_fl, id_var_time_fl, 'standard_name', 'time', 000)
          CALL netcdf_create_att( id_set_fl, id_var_time_fl, 'axis', 'T', 000)

          DO l = 1, num_leg
             CALL netcdf_create_dim( id_set_fl, dofl_dim_label_x(l),           &
                                     NF90_UNLIMITED, id_dim_x_fl(l), 250 )
             CALL netcdf_create_dim( id_set_fl, dofl_dim_label_y(l),           &
                                     NF90_UNLIMITED, id_dim_y_fl(l), 250 )
             CALL netcdf_create_dim( id_set_fl, dofl_dim_label_z(l),           &
                                     NF90_UNLIMITED, id_dim_z_fl(l), 250 )

             CALL netcdf_create_var( id_set_fl, (/ id_dim_x_fl(l) /),          &
                                     dofl_dim_label_x(l), NF90_DOUBLE,         &
                                     id_var_x_fl(l), 'm', '', 251, 252, 000 )
             CALL netcdf_create_var( id_set_fl, (/ id_dim_y_fl(l) /),          &
                                     dofl_dim_label_y(l), NF90_DOUBLE,         &
                                     id_var_y_fl(l), 'm', '', 251, 252, 000 )
             CALL netcdf_create_var( id_set_fl, (/ id_dim_z_fl(l) /),          &
                                     dofl_dim_label_z(l), NF90_DOUBLE,         &
                                     id_var_z_fl(l), 'm', '', 251, 252, 000 )
          ENDDO
!
!--       Define the variables
          var_list = ';'
          k = 1
          DO  l = 1, num_leg
             DO i = 1, num_var_fl

                CALL netcdf_create_var( id_set_fl, (/ id_dim_time_fl /),       &
                                        dofl_label(k), nc_precision(9),        &
                                        id_var_dofl(k),                        &
                                        TRIM( dofl_unit(k) ),                  &
                                        TRIM( dofl_label(k) ), 253, 254, 255 )

                k = k + 1

             ENDDO

          ENDDO

!
!--       Write the list of variables as global attribute (this is used by
!--       restart runs)
          nc_stat = NF90_PUT_ATT( id_set_fl, NF90_GLOBAL, 'VAR_LIST', var_list )
          CALL netcdf_handle_error( 'netcdf_define_header', 258 )

!
!--       Leave netCDF define mode
          nc_stat = NF90_ENDDEF( id_set_fl )
          CALL netcdf_handle_error( 'netcdf_define_header', 259 )


       CASE ( 'fl_ext' )

!
!--       Get the list of variables and compare with the actual run.
!--       First var_list_old has to be reset, since GET_ATT does not assign
!--       trailing blanks.
          var_list_old = ' '
          nc_stat = NF90_GET_ATT( id_set_fl, NF90_GLOBAL, 'VAR_LIST',          &
                                  var_list_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 260 )

          var_list = ';'
          i = 1
          DO  i = 1, num_leg * num_var_fl

             var_list = TRIM( var_list ) // TRIM( dofl_label(i) ) // ';'

          ENDDO

          IF ( TRIM( var_list ) /= TRIM( var_list_old ) )  THEN
             message_string = 'netCDF file for flight time series ' //         &
                              'from previous run found,' //                    &
                              '&but this file cannot be extended due to' //    &
                              ' variable mismatch.' //                         &
                              '&New file is created instead.'
             CALL message( 'define_netcdf_header', 'PA0257', 0, 1, 0, 6, 0 )
             extend = .FALSE.
             RETURN
          ENDIF

!
!--       Get the id of the time coordinate (unlimited coordinate) and its
!--       last index on the file. The next time level is dofl_time_count+1.
!--       The current time must be larger than the last output time
!--       on the file.
          nc_stat = NF90_INQ_VARID( id_set_fl, 'time', id_var_time_fl )
          CALL netcdf_handle_error( 'netcdf_define_header', 261 )

          nc_stat = NF90_INQUIRE_VARIABLE( id_set_fl, id_var_time_fl,          &
                                           dimids = id_dim_time_old )
          CALL netcdf_handle_error( 'netcdf_define_header', 262 )
          id_dim_time_fl = id_dim_time_old(1)

          nc_stat = NF90_INQUIRE_DIMENSION( id_set_fl, id_dim_time_fl,         &
                                            len = dofl_time_count )
          CALL netcdf_handle_error( 'netcdf_define_header', 263 )

          nc_stat = NF90_GET_VAR( id_set_fl, id_var_time_fl,        &
                                  last_time_coordinate,             &
                                  start = (/ dofl_time_count /), &
                                  count = (/ 1 /) )
          CALL netcdf_handle_error( 'netcdf_define_header', 264 )

          IF ( last_time_coordinate(1) >= simulated_time )  THEN
             message_string = 'netCDF file for flight-time series ' //         &
                              'from previous run found,' //                    &
                              '&but this file cannot be extended becaus' //    &
                              'e the current output time' //                   &
                              '&is less or equal than the last output t' //    &
                              'ime on this file.' //                           &
                              '&New file is created instead.'
             CALL message( 'define_netcdf_header', 'PA0258', 0, 1, 0, 6, 0 )
             dofl_time_count = 0
             extend = .FALSE.
             RETURN
          ENDIF

!
!--       Dataset seems to be extendable.
!--       Now get the remaining dimension and variable ids
          DO l = 1, num_leg
             nc_stat = NF90_INQ_VARID( id_set_fl, dofl_dim_label_x(l),         &
                                       id_var_x_fl(l) )
             CALL netcdf_handle_error( 'netcdf_define_header', 265 )
             nc_stat = NF90_INQ_VARID( id_set_fl, dofl_dim_label_y(l),         &
                                       id_var_y_fl(l) )
             CALL netcdf_handle_error( 'netcdf_define_header', 265 )
             nc_stat = NF90_INQ_VARID( id_set_fl, dofl_dim_label_z(l),         &
                                       id_var_z_fl(l) )
             CALL netcdf_handle_error( 'netcdf_define_header', 265 )

          ENDDO


          DO  i = 1, num_leg * num_var_fl

            nc_stat = NF90_INQ_VARID( id_set_fl, dofl_label(i), &
                                       id_var_dofl(i) )
            CALL netcdf_handle_error( 'netcdf_define_header', 265 )

          ENDDO

!
!--       Update the title attribute on file
!--       In order to avoid 'data mode' errors if updated attributes are larger
!--       than their original size, NF90_PUT_ATT is called in 'define mode'
!--       enclosed by NF90_REDEF and NF90_ENDDEF calls. This implies a possible
!--       performance loss due to data copying; an alternative strategy would be
!--       to ensure equal attribute size in a job chain. Maybe revise later.
          nc_stat = NF90_REDEF( id_set_fl )
          CALL netcdf_handle_error( 'netcdf_define_header', 439 )
          nc_stat = NF90_PUT_ATT( id_set_fl, NF90_GLOBAL, 'title',             &
                                  TRIM( run_description_header ) )
          CALL netcdf_handle_error( 'netcdf_define_header', 267 )
          nc_stat = NF90_ENDDEF( id_set_fl )
          CALL netcdf_handle_error( 'netcdf_define_header', 440 )
          message_string = 'netCDF file for flight-time series ' //            &
                           'from previous run found.' //                       &
                           '&This file will be extended.'
          CALL message( 'define_netcdf_header', 'PA0259', 0, 0, 0, 6, 0 )


       CASE DEFAULT

          message_string = 'mode "' // TRIM( mode) // '" not supported'
          CALL message( 'netcdf_define_header', 'PA0270', 0, 0, 0, 6, 0 )

    END SELECT

#endif
 END SUBROUTINE netcdf_define_header


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Creates a netCDF file and give back the id. The parallel flag has to be TRUE
!> for parallel netCDF output support.
!------------------------------------------------------------------------------!

 SUBROUTINE netcdf_create_file( filename , id, parallel, errno )
#if defined( __netcdf )

    USE pegrid

    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT(IN) :: filename
    INTEGER, INTENT(IN)           :: errno
    INTEGER, INTENT(OUT)          :: id
    INTEGER                       :: idum  !< dummy variable used to avoid compiler warnings about unused variables
    LOGICAL, INTENT(IN)           :: parallel

!
!-- Next line is just to avoid compiler warning about unused variable
    IF ( parallel )  idum = 0

!
!-- Create a new netCDF output file with requested netCDF format
    IF ( netcdf_data_format == 1 )  THEN
!
!--    Classic netCDF format
       nc_stat = NF90_CREATE( filename, NF90_NOCLOBBER, id )

    ELSEIF ( netcdf_data_format == 2 )  THEN
!
!--    64bit-offset format
       nc_stat = NF90_CREATE( filename,                                        &
                              IOR( NF90_NOCLOBBER, NF90_64BIT_OFFSET ), id )

#if defined( __netcdf4 )
    ELSEIF ( netcdf_data_format == 3  .OR.                                     &
             ( .NOT. parallel  .AND.  netcdf_data_format == 5 ) )  THEN
!
!--    netCDF4/HDF5 format
       nc_stat = NF90_CREATE( filename, IOR( NF90_NOCLOBBER, NF90_NETCDF4 ), id )

    ELSEIF ( netcdf_data_format == 4  .OR.                                     &
             ( .NOT. parallel  .AND.  netcdf_data_format == 6 ) )  THEN
!
!--    netCDF4/HDF5 format with classic model flag
       nc_stat = NF90_CREATE( filename,                                        &
                              IOR( NF90_NOCLOBBER,                             &
                              IOR( NF90_CLASSIC_MODEL, NF90_HDF5 ) ), id )

#if defined( __netcdf4_parallel )
    ELSEIF ( netcdf_data_format == 5  .AND.  parallel )  THEN
!
!--    netCDF4/HDF5 format, parallel
       nc_stat = NF90_CREATE( filename,                                        &
                              IOR( NF90_NOCLOBBER,                             &
                              IOR( NF90_NETCDF4, NF90_MPIIO ) ),               &
                              id, COMM = comm2d, INFO = MPI_INFO_NULL )

    ELSEIF ( netcdf_data_format == 6  .AND.  parallel )  THEN
!
!--    netCDF4/HDF5 format with classic model flag, parallel
       nc_stat = NF90_CREATE( filename,                                        &
                              IOR( NF90_NOCLOBBER,                             &
                              IOR( NF90_MPIIO,                                 &
                              IOR( NF90_CLASSIC_MODEL, NF90_HDF5 ) ) ),        &
                              id, COMM = comm2d, INFO = MPI_INFO_NULL )

#endif
#endif
    ENDIF

    CALL netcdf_handle_error( 'netcdf_create_file', errno )
#endif
 END SUBROUTINE netcdf_create_file

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Opens an existing netCDF file for writing and gives back the id.
!> The parallel flag has to be TRUE for parallel netCDF output support.
!------------------------------------------------------------------------------!
 SUBROUTINE netcdf_open_write_file( filename, id, parallel, errno )
#if defined( __netcdf )

    USE pegrid

    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT(IN) :: filename
    INTEGER, INTENT(IN)           :: errno
    INTEGER, INTENT(OUT)          :: id
    LOGICAL, INTENT(IN)           :: parallel


    IF ( netcdf_data_format < 5  .OR.  .NOT. parallel )  THEN
       nc_stat = NF90_OPEN( filename, NF90_WRITE, id )
#if defined( __netcdf4 )
#if defined( __netcdf4_parallel )
    ELSEIF ( netcdf_data_format > 4  .AND.  parallel )  THEN
       nc_stat = NF90_OPEN( filename, IOR( NF90_WRITE, NF90_MPIIO ), id,  &
                            COMM = comm2d, INFO = MPI_INFO_NULL )
#endif
#endif
    ENDIF

    CALL netcdf_handle_error( 'netcdf_open_write_file', errno )
#endif
 END SUBROUTINE netcdf_open_write_file


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Prints out a text message corresponding to the current status.
!------------------------------------------------------------------------------!

 SUBROUTINE netcdf_handle_error( routine_name, errno )
#if defined( __netcdf )


    USE control_parameters,                                                    &
        ONLY:  message_string

    IMPLICIT NONE

    CHARACTER(LEN=6) ::  message_identifier
    CHARACTER(LEN=*) ::  routine_name

    INTEGER(iwp) ::  errno

    IF ( nc_stat /= NF90_NOERR )  THEN

       WRITE( message_identifier, '(''NC'',I4.4)' )  errno

       message_string = TRIM( NF90_STRERROR( nc_stat ) )

       CALL message( routine_name, message_identifier, 2, 2, 0, 6, 1 )

    ENDIF

#endif
 END SUBROUTINE netcdf_handle_error


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Create a dimension in NetCDF file
!------------------------------------------------------------------------------!

 SUBROUTINE netcdf_create_dim(ncid, dim_name, ncdim_type, ncdim_id, error_no)

#if defined( __netcdf )

    USE kinds

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) ::  dim_name

    INTEGER, INTENT(IN)  ::  error_no
    INTEGER, INTENT(IN)  ::  ncid
    INTEGER, INTENT(OUT) ::  ncdim_id
    INTEGER, INTENT(IN)  ::  ncdim_type

!
!-- Define time coordinate for volume data (unlimited dimension)
    nc_stat = NF90_DEF_DIM( ncid, dim_name, ncdim_type, ncdim_id )
    CALL netcdf_handle_error( 'netcdf_create_dim', error_no )

#endif

 END SUBROUTINE netcdf_create_dim


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Create a one dimensional variable in specific units in NetCDF file
!------------------------------------------------------------------------------!

 SUBROUTINE netcdf_create_var( ncid, dim_id, var_name, var_type, var_id,       &
                               unit_name, long_name, error_no1, error_no2,     &
                               error_no3, fill )

#if defined( __netcdf )
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) ::  long_name
    CHARACTER(LEN=*), INTENT(IN) ::  unit_name
    CHARACTER(LEN=*), INTENT(IN) ::  var_name

    LOGICAL, OPTIONAL ::  fill  !< indicates setting of _FillValue attribute

    INTEGER, INTENT(IN)  ::  error_no1
    INTEGER, INTENT(IN)  ::  error_no2
    INTEGER, INTENT(IN)  ::  error_no3
    INTEGER, INTENT(IN)  ::  ncid
    INTEGER, INTENT(OUT) ::  var_id
    INTEGER, INTENT(IN)  ::  var_type

    INTEGER, DIMENSION(:), INTENT(IN) ::  dim_id

!
!-- Define variable
    nc_stat = NF90_DEF_VAR( ncid, var_name, var_type, dim_id, var_id )
    CALL netcdf_handle_error( 'netcdf_create_var', error_no1 )

#if defined( __netcdf4 )
!
!-- Check if variable should be deflate (including shuffling)
!-- and if it is possible (only NetCDF4 with HDF5 supports compression)
    IF ( netcdf_data_format > 2  .AND.  netcdf_deflate > 0 )  THEN
       nc_stat = NF90_DEF_VAR_DEFLATE( ncid, var_id, 1, 1, netcdf_deflate )
       CALL netcdf_handle_error( 'netcdf_create_var_deflate', error_no1 )
    ENDIF
#endif
!
!-- Set unit name if set
    IF ( unit_name /= '' )  THEN
       nc_stat = NF90_PUT_ATT( ncid, var_id, 'units', unit_name )
       CALL netcdf_handle_error( 'netcdf_create_var', error_no2 )
    ENDIF

!
!-- Set long name if set
    IF ( long_name /= '' )  THEN
       nc_stat = NF90_PUT_ATT( ncid, var_id, 'long_name', long_name )
       CALL netcdf_handle_error( 'netcdf_create_var', error_no3 )
    ENDIF

!
!-- Set _FillValue for all variables, except for dimension variables.
!-- Set the fill values accordingly to the corresponding output precision.
    IF ( PRESENT( fill ) )  THEN
       IF ( var_type == NF90_REAL4 )  THEN
          nc_stat = NF90_PUT_ATT( ncid, var_id, '_FillValue',                  &
                                  REAL( fill_value, KIND = 4 ) )
          CALL netcdf_handle_error( 'netcdf_create_var', 0 )
       ELSE
          nc_stat = NF90_PUT_ATT( ncid, var_id, '_FillValue',                  &
                                  REAL( fill_value, KIND = 8 ) )
          CALL netcdf_handle_error( 'netcdf_create_var', 0 )
       ENDIF
    ENDIF

#endif
 END SUBROUTINE netcdf_create_var


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Write attributes to file.
!------------------------------------------------------------------------------!
 SUBROUTINE netcdf_create_att_string( ncid, varid, name, value, err )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) ::  name    !< attribute name
    CHARACTER(LEN=*), INTENT(IN) ::  value   !< attribute value

    INTEGER, INTENT(IN) ::  err              !< error id
    INTEGER, INTENT(IN) ::  ncid             !< file id

    INTEGER, INTENT(IN), OPTIONAL ::  varid  !< variable id

#if defined( __netcdf )
    IF ( PRESENT( varid ) )  THEN
       nc_stat = NF90_PUT_ATT( ncid, varid, TRIM( name ), TRIM( value ) )
    ELSE
       nc_stat = NF90_PUT_ATT( ncid, NF90_GLOBAL, TRIM( name ), TRIM( value ) )
    ENDIF
    CALL netcdf_handle_error( 'netcdf_create_att_string', err )
#endif

 END SUBROUTINE netcdf_create_att_string


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Write a set of global attributes to file.
!------------------------------------------------------------------------------!
 SUBROUTINE netcdf_create_global_atts( ncid, data_content, title, error_no )

    USE control_parameters,                                                    &
        ONLY:  revision, run_date, run_time, run_zone, runnr, version

    USE netcdf_data_input_mod,                                                 &
        ONLY:  input_file_atts

    USE palm_date_time_mod,                                                    &
        ONLY:  date_time_str_len, get_date_time

    IMPLICIT NONE

    CHARACTER(LEN=date_time_str_len) ::  origin_time_string  !< string containing date-time of origin

    CHARACTER(LEN=*), INTENT(IN)  ::  data_content  !< describes the type of data in file
    CHARACTER(LEN=*), INTENT(IN)  ::  title         !< file title

    INTEGER, INTENT(IN)  ::  error_no  !< error number
    INTEGER, INTENT(IN)  ::  ncid      !< file id
!
!-- Get date-time string for origin_time
    CALL get_date_time( 0.0_wp, date_time_str=origin_time_string )

#if defined( __netcdf )
    nc_stat = NF90_PUT_ATT( ncid, NF90_GLOBAL, 'title', TRIM( title ) )
    CALL netcdf_handle_error( 'netcdf_create_global_atts 1', error_no )
    nc_stat = NF90_PUT_ATT( ncid, NF90_GLOBAL, 'Conventions', 'CF-1.7' )
    CALL netcdf_handle_error( 'netcdf_create_global_atts 2', error_no )
    nc_stat = NF90_PUT_ATT( ncid, NF90_GLOBAL, 'creation_time', TRIM( run_date )//' '//TRIM( run_time )//' '//run_zone(1:3) )
    CALL netcdf_handle_error( 'netcdf_create_global_atts 3', error_no )
    nc_stat = NF90_PUT_ATT( ncid, NF90_GLOBAL, 'data_content', TRIM(data_content) )
    CALL netcdf_handle_error( 'netcdf_create_global_atts 4', error_no )
    nc_stat = NF90_PUT_ATT( ncid, NF90_GLOBAL, 'version', runnr+1 )
    CALL netcdf_handle_error( 'netcdf_create_global_atts 5', error_no )

    nc_stat = NF90_PUT_ATT( ncid, NF90_GLOBAL, 'origin_time', origin_time_string )
    CALL netcdf_handle_error( 'netcdf_create_global_atts 6', error_no )
    nc_stat = NF90_PUT_ATT( ncid, NF90_GLOBAL, 'origin_lat', init_model%latitude )
    CALL netcdf_handle_error( 'netcdf_create_global_atts 7', error_no )
    nc_stat = NF90_PUT_ATT( ncid, NF90_GLOBAL, 'origin_lon', init_model%longitude )
    CALL netcdf_handle_error( 'netcdf_create_global_atts 8', error_no )
    nc_stat = NF90_PUT_ATT( ncid, NF90_GLOBAL, 'origin_x', init_model%origin_x )
    CALL netcdf_handle_error( 'netcdf_create_global_atts 9', error_no )
    nc_stat = NF90_PUT_ATT( ncid, NF90_GLOBAL, 'origin_y', init_model%origin_y )
    CALL netcdf_handle_error( 'netcdf_create_global_atts 10', error_no )
    nc_stat = NF90_PUT_ATT( ncid, NF90_GLOBAL, 'origin_z', init_model%origin_z )
    CALL netcdf_handle_error( 'netcdf_create_global_atts 11', error_no )
    nc_stat = NF90_PUT_ATT( ncid, NF90_GLOBAL, 'rotation_angle', rotation_angle )
    CALL netcdf_handle_error( 'netcdf_create_global_atts 12', error_no )

    nc_stat = NF90_PUT_ATT( ncid, NF90_GLOBAL, 'dependencies', '' )
    CALL netcdf_handle_error( 'netcdf_create_global_atts 13', error_no )
    nc_stat = NF90_PUT_ATT( ncid, NF90_GLOBAL, 'history', '' )
    CALL netcdf_handle_error( 'netcdf_create_global_atts 14', error_no )

    nc_stat = NF90_PUT_ATT( ncid, NF90_GLOBAL, TRIM( input_file_atts%author_char ), TRIM( input_file_atts%author ) )
    CALL netcdf_handle_error( 'netcdf_create_global_atts 15', error_no )
    nc_stat = NF90_PUT_ATT( ncid, NF90_GLOBAL, TRIM( input_file_atts%contact_person_char ), TRIM( input_file_atts%contact_person ) )
    CALL netcdf_handle_error( 'netcdf_create_global_atts 16', error_no )
    nc_stat = NF90_PUT_ATT( ncid, NF90_GLOBAL, TRIM( input_file_atts%institution_char ), TRIM( input_file_atts%institution ) )
    CALL netcdf_handle_error( 'netcdf_create_global_atts 17', error_no )
    nc_stat = NF90_PUT_ATT( ncid, NF90_GLOBAL, TRIM( input_file_atts%acronym_char ), TRIM( input_file_atts%acronym ) )
    CALL netcdf_handle_error( 'netcdf_create_global_atts 18', error_no )

    nc_stat = NF90_PUT_ATT( ncid, NF90_GLOBAL, TRIM( input_file_atts%campaign_char ), TRIM( input_file_atts%campaign ) )
    CALL netcdf_handle_error( 'netcdf_create_global_atts 19', error_no )
    nc_stat = NF90_PUT_ATT( ncid, NF90_GLOBAL, TRIM( input_file_atts%location_char ), TRIM( input_file_atts%location ) )
    CALL netcdf_handle_error( 'netcdf_create_global_atts 20', error_no )
    nc_stat = NF90_PUT_ATT( ncid, NF90_GLOBAL, TRIM( input_file_atts%site_char ), TRIM( input_file_atts%site ) )
    CALL netcdf_handle_error( 'netcdf_create_global_atts 21', error_no )

    nc_stat = NF90_PUT_ATT( ncid, NF90_GLOBAL, 'source', TRIM( version )//' '//TRIM( revision ) )
    CALL netcdf_handle_error( 'netcdf_create_global_atts 22', error_no )
    nc_stat = NF90_PUT_ATT( ncid, NF90_GLOBAL, TRIM( input_file_atts%references_char ), TRIM( input_file_atts%references ) )
    CALL netcdf_handle_error( 'netcdf_create_global_atts 23', error_no )
    nc_stat = NF90_PUT_ATT( ncid, NF90_GLOBAL, TRIM( input_file_atts%keywords_char ), TRIM( input_file_atts%keywords ) )
    CALL netcdf_handle_error( 'netcdf_create_global_atts 24', error_no )
    nc_stat = NF90_PUT_ATT( ncid, NF90_GLOBAL, TRIM( input_file_atts%licence_char ), TRIM( input_file_atts%licence ) )
    CALL netcdf_handle_error( 'netcdf_create_global_atts 25', error_no )
    nc_stat = NF90_PUT_ATT( ncid, NF90_GLOBAL, TRIM( input_file_atts%comment_char ), TRIM( input_file_atts%comment ) )
    CALL netcdf_handle_error( 'netcdf_create_global_atts 26', error_no )

#endif

 END SUBROUTINE netcdf_create_global_atts

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Create a variable holding the coordinate-reference-system information.
!------------------------------------------------------------------------------!
 SUBROUTINE netcdf_create_crs( ncid, error_no )

    IMPLICIT NONE

    INTEGER, INTENT(IN)  ::  error_no  !< error number
    INTEGER, INTENT(IN)  ::  ncid      !< file id
    INTEGER              ::  var_id    !< variable id

#if defined( __netcdf )
!
!-- Define variable
    nc_stat = NF90_DEF_VAR( ncid, 'crs', NF90_INT, VARID = var_id )
    CALL netcdf_handle_error( 'netcdf_create_crs', error_no )
!
!-- Set attributes
    nc_stat = NF90_PUT_ATT( ncid, var_id, 'epsg_code', &
                            coord_ref_sys%epsg_code )
    CALL netcdf_handle_error( 'netcdf_create_crs', error_no )

    nc_stat = NF90_PUT_ATT( ncid, var_id, 'false_easting', &
                            coord_ref_sys%false_easting )
    CALL netcdf_handle_error( 'netcdf_create_crs', error_no )

    nc_stat = NF90_PUT_ATT( ncid, var_id, 'false_northing', &
                            coord_ref_sys%false_northing )
    CALL netcdf_handle_error( 'netcdf_create_crs', error_no )

    nc_stat = NF90_PUT_ATT( ncid, var_id, 'grid_mapping_name', &
                            coord_ref_sys%grid_mapping_name )
    CALL netcdf_handle_error( 'netcdf_create_crs', error_no )

    nc_stat = NF90_PUT_ATT( ncid, var_id, 'inverse_flattening', &
                            coord_ref_sys%inverse_flattening )
    CALL netcdf_handle_error( 'netcdf_create_crs', error_no )

    nc_stat = NF90_PUT_ATT( ncid, var_id, 'latitude_of_projection_origin', &
                            coord_ref_sys%latitude_of_projection_origin )
    CALL netcdf_handle_error( 'netcdf_create_crs', error_no )

    nc_stat = NF90_PUT_ATT( ncid, var_id, 'long_name', &
                            coord_ref_sys%long_name )
    CALL netcdf_handle_error( 'netcdf_create_crs', error_no )

    nc_stat = NF90_PUT_ATT( ncid, var_id, 'longitude_of_central_meridian', &
                            coord_ref_sys%longitude_of_central_meridian )
    CALL netcdf_handle_error( 'netcdf_create_crs', error_no )

    nc_stat = NF90_PUT_ATT( ncid, var_id, 'longitude_of_prime_meridian', &
                            coord_ref_sys%longitude_of_prime_meridian )
    CALL netcdf_handle_error( 'netcdf_create_crs', error_no )

    nc_stat = NF90_PUT_ATT( ncid, var_id, 'scale_factor_at_central_meridian', &
                            coord_ref_sys%scale_factor_at_central_meridian )
    CALL netcdf_handle_error( 'netcdf_create_crs', error_no )

    nc_stat = NF90_PUT_ATT( ncid, var_id, 'semi_major_axis', &
                            coord_ref_sys%semi_major_axis )
    CALL netcdf_handle_error( 'netcdf_create_crs', error_no )

    nc_stat = NF90_PUT_ATT( ncid, var_id, 'units', &
                            coord_ref_sys%units )
    CALL netcdf_handle_error( 'netcdf_create_crs', error_no )

#endif
 END SUBROUTINE netcdf_create_crs


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Define UTM coordinates and longitude and latitude in file.
!------------------------------------------------------------------------------!
 SUBROUTINE define_geo_coordinates( id_set, id_dim_x, id_dim_y, id_var_eutm, id_var_nutm, id_var_lat, id_var_lon )

    IMPLICIT NONE

    INTEGER ::  i                    !< loop index
    INTEGER, INTENT(IN)  ::  id_set  !< file id

    INTEGER(iwp), DIMENSION(0:1), INTENT(IN) ::  id_dim_x  !< dimension id of x and xu
    INTEGER(iwp), DIMENSION(0:1), INTENT(IN) ::  id_dim_y  !< dimension id of y and yv

    INTEGER(iwp), DIMENSION(0:2), INTENT(OUT) ::  id_var_eutm  !< variable id for E_UTM coordinates
    INTEGER(iwp), DIMENSION(0:2), INTENT(OUT) ::  id_var_lat   !< variable id for latitude coordinates
    INTEGER(iwp), DIMENSION(0:2), INTENT(OUT) ::  id_var_lon   !< variable id for longitude coordinates
    INTEGER(iwp), DIMENSION(0:2), INTENT(OUT) ::  id_var_nutm  !< variable id for N_UTM coordinates

!
!-- Initialize return values
    id_var_lat  = -1
    id_var_lon  = -1
    id_var_eutm = -1
    id_var_nutm = -1

#if defined( __netcdf )
!
!-- Define UTM coordinates
    IF ( rotation_angle == 0.0_wp )  THEN
       CALL netcdf_create_var( id_set, (/ id_dim_x(0) /), 'E_UTM', NF90_DOUBLE, id_var_eutm(0), 'm', 'easting', 000, 000, 000 )
       CALL netcdf_create_var( id_set, (/ id_dim_y(0) /), 'N_UTM', NF90_DOUBLE, id_var_nutm(0), 'm', 'northing', 000, 000, 000 )
       CALL netcdf_create_var( id_set, (/ id_dim_x(1) /), 'Eu_UTM', NF90_DOUBLE, id_var_eutm(1), 'm', 'easting', 000, 000, 000 )
       CALL netcdf_create_var( id_set, (/ id_dim_y(0) /), 'Nu_UTM', NF90_DOUBLE, id_var_nutm(1), 'm', 'northing', 000, 000, 000 )
       CALL netcdf_create_var( id_set, (/ id_dim_x(0) /), 'Ev_UTM', NF90_DOUBLE, id_var_eutm(2), 'm', 'easting', 000, 000, 000 )
       CALL netcdf_create_var( id_set, (/ id_dim_y(1) /), 'Nv_UTM', NF90_DOUBLE, id_var_nutm(2), 'm', 'northing', 000, 000, 000 )
    ELSE
       CALL netcdf_create_var( id_set, (/ id_dim_x(0), id_dim_y(0) /), &
                               'E_UTM', NF90_DOUBLE, id_var_eutm(0), 'm', 'easting', 000, 000, 000 )
       CALL netcdf_create_var( id_set, (/ id_dim_x(0), id_dim_y(0) /), &
                               'N_UTM', NF90_DOUBLE, id_var_nutm(0), 'm', 'northing', 000, 000, 000 )
       CALL netcdf_create_var( id_set, (/ id_dim_x(1), id_dim_y(0) /), &
                               'Eu_UTM', NF90_DOUBLE, id_var_eutm(1), 'm', 'easting', 000, 000, 000 )
       CALL netcdf_create_var( id_set, (/ id_dim_x(1), id_dim_y(0) /), &
                               'Nu_UTM', NF90_DOUBLE, id_var_nutm(1), 'm', 'northing', 000, 000, 000 )
       CALL netcdf_create_var( id_set, (/ id_dim_x(0), id_dim_y(1) /), &
                               'Ev_UTM', NF90_DOUBLE, id_var_eutm(2), 'm', 'easting', 000, 000, 000 )
       CALL netcdf_create_var( id_set, (/ id_dim_x(0), id_dim_y(1) /), &
                               'Nv_UTM', NF90_DOUBLE, id_var_nutm(2), 'm', 'northing', 000, 000, 000 )
    ENDIF
!
!-- Define geographic coordinates
    CALL netcdf_create_var( id_set, (/ id_dim_x(0), id_dim_y(0) /), 'lon', NF90_DOUBLE, id_var_lon(0), &
                            'degrees_east', 'longitude', 000, 000, 000 )
    CALL netcdf_create_var( id_set, (/ id_dim_x(0), id_dim_y(0) /), 'lat', NF90_DOUBLE, id_var_lat(0), &
                            'degrees_north', 'latitude', 000, 000, 000 )
    CALL netcdf_create_var( id_set, (/ id_dim_x(1), id_dim_y(0) /), 'lonu', NF90_DOUBLE, id_var_lon(1), &
                            'degrees_east', 'longitude', 000, 000, 000 )
    CALL netcdf_create_var( id_set, (/ id_dim_x(1), id_dim_y(0) /), 'latu', NF90_DOUBLE, id_var_lat(1), &
                            'degrees_north', 'latitude', 000, 000, 000 )
    CALL netcdf_create_var( id_set, (/ id_dim_x(0), id_dim_y(1) /), 'lonv', NF90_DOUBLE, id_var_lon(2), &
                            'degrees_east', 'longitude', 000, 000, 000 )
    CALL netcdf_create_var( id_set, (/ id_dim_x(0), id_dim_y(1) /), 'latv', NF90_DOUBLE, id_var_lat(2), &
                            'degrees_north', 'latitude', 000, 000, 000 )

    DO  i = 0, 2
       CALL netcdf_create_att( id_set, id_var_eutm(i), 'standard_name', 'projection_x_coordinate', 000)
       CALL netcdf_create_att( id_set, id_var_nutm(i), 'standard_name', 'projection_y_coordinate', 000)

       CALL netcdf_create_att( id_set, id_var_lat(i), 'standard_name', 'latitude', 000)
       CALL netcdf_create_att( id_set, id_var_lon(i), 'standard_name', 'longitude', 000)
    ENDDO

#endif
 END SUBROUTINE define_geo_coordinates


 END MODULE netcdf_interface
