!> @file basic_constants_and_equations_mod.f90
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
! $Id: basic_constants_and_equations_mod.f90 4400 2020-02-10 20:32:41Z suehring $
! Move routine to transform coordinates from netcdf_interface_mod to 
! basic_constants_and_equations_mod
! 
! 4360 2020-01-07 11:25:50Z suehring
! Corrected "Former revisions" section
! 
! 4088 2019-07-11 13:57:56Z Giersch
! Comment of barometric formula improved, function for ideal gas law revised 
! 
! 4084 2019-07-10 17:09:11Z knoop
! Changed precomputed fractions to be variable based
! 
! 4055 2019-06-27 09:47:29Z suehring
! Added rgas_univ (universal gas constant) (E.C. Chan)
! 
! 
! 3655 2019-01-07 16:51:22Z knoop
! OpenACC port for SPEC
! 3361 2018-10-16 20:39:37Z knoop
! New module (introduced with modularization of bulk cloud physics model)
!
! 
!
!
! Description:
! ------------
!> This module contains all basic (physical) constants
!> and
!> functions for the calculation of diagnostic quantities.
!------------------------------------------------------------------------------!
 MODULE basic_constants_and_equations_mod


    USE kinds

    IMPLICIT NONE

    REAL(wp), PARAMETER ::  c_p = 1005.0_wp                           !< heat capacity of dry air (J kg-1 K-1)
    REAL(wp), PARAMETER ::  degc_to_k = 273.15_wp                     !< temperature (in K) of 0 deg C (K)
    REAL(wp), PARAMETER ::  g = 9.81_wp                               !< gravitational acceleration (m s-2)
    REAL(wp), PARAMETER ::  kappa = 0.4_wp                            !< von Karman constant
    REAL(wp), PARAMETER ::  l_m = 0.33E+06_wp                         !< latent heat of water melting (J kg-1)
    REAL(wp), PARAMETER ::  l_v = 2.5E+06_wp                          !< latent heat of water vaporization (J kg-1)
    REAL(wp), PARAMETER ::  l_s = l_m + l_v                           !< latent heat of water sublimation (J kg-1)
    REAL(wp), PARAMETER ::  molecular_weight_of_nacl = 0.05844_wp     !< mol. m. NaCl (kg mol-1)
    REAL(wp), PARAMETER ::  molecular_weight_of_c3h4o4 = 0.10406_wp   !< mol. m. malonic acid (kg mol-1)
    REAL(wp), PARAMETER ::  molecular_weight_of_nh4no3 = 0.08004_wp   !< mol. m. ammonium sulfate (kg mol-1)
    REAL(wp), PARAMETER ::  molecular_weight_of_water = 0.01801528_wp !< mol. m. H2O (kg mol-1)
    REAL(wp), PARAMETER ::  pi = 3.141592654_wp                       !< PI
    !$ACC DECLARE COPYIN(pi)
    REAL(wp), PARAMETER ::  rgas_univ = 8.31446261815324_wp           !< universal gas constant (J K-1 mol-1)
    REAL(wp), PARAMETER ::  rho_l = 1.0E3_wp                          !< density of water (kg m-3)
    REAL(wp), PARAMETER ::  rho_nacl = 2165.0_wp                      !< density of NaCl (kg m-3)
    REAL(wp), PARAMETER ::  rho_c3h4o4 = 1600.0_wp                    !< density of malonic acid (kg m-3)
    REAL(wp), PARAMETER ::  rho_nh4no3 = 1720.0_wp                    !< density of ammonium sulfate (kg m-3)
    REAL(wp), PARAMETER ::  r_d = 287.0_wp                            !< sp. gas const. dry air (J kg-1 K-1)
    REAL(wp), PARAMETER ::  r_v = 461.51_wp                           !< sp. gas const. water vapor (J kg-1 K-1)
    REAL(wp), PARAMETER ::  sigma_sb = 5.67037E-08_wp                 !< Stefan-Boltzmann constant
    REAL(wp), PARAMETER ::  solar_constant = 1368.0_wp                !< solar constant at top of atmosphere
    REAL(wp), PARAMETER ::  vanthoff_nacl = 2.0_wp                    !< van't Hoff factor for NaCl
    REAL(wp), PARAMETER ::  vanthoff_c3h4o4 = 1.37_wp                 !< van't Hoff factor for malonic acid
    REAL(wp), PARAMETER ::  vanthoff_nh4no3 = 2.31_wp                 !< van't Hoff factor for ammonium sulfate

    REAL(wp), PARAMETER ::  p_0 = 100000.0_wp                         !< standard pressure reference state

    REAL(wp), PARAMETER ::  g_d_cp  = g   / c_p   !< precomputed g / c_p
    REAL(wp), PARAMETER ::  lv_d_cp = l_v / c_p   !< precomputed l_v / c_p
    REAL(wp), PARAMETER ::  lv_d_rd = l_v / r_d   !< precomputed l_v / r_d
    REAL(wp), PARAMETER ::  rd_d_rv = r_d / r_v   !< precomputed r_d / r_v
    REAL(wp), PARAMETER ::  rd_d_cp = r_d / c_p   !< precomputed r_d / c_p
    REAL(wp), PARAMETER ::  cp_d_rd = c_p / r_d   !< precomputed c_p / r_d

    REAL(wp) ::  molecular_weight_of_solute = molecular_weight_of_nacl  !< mol. m. NaCl (kg mol-1)
    REAL(wp) ::  rho_s = rho_nacl                                       !< density of NaCl (kg m-3)
    REAL(wp) ::  vanthoff = vanthoff_nacl                               !< van't Hoff factor for NaCl


    SAVE

    PRIVATE magnus_0d, &
            magnus_1d, &
            ideal_gas_law_rho_0d, &
            ideal_gas_law_rho_1d, &
            ideal_gas_law_rho_pt_0d, &
            ideal_gas_law_rho_pt_1d, &
            exner_function_0d, &
            exner_function_1d, &
            exner_function_invers_0d, &
            exner_function_invers_1d, &
            barometric_formula_0d, &
            barometric_formula_1d


    INTERFACE convert_utm_to_geographic
       MODULE PROCEDURE convert_utm_to_geographic
    END INTERFACE convert_utm_to_geographic

    INTERFACE magnus
       MODULE PROCEDURE magnus_0d
       MODULE PROCEDURE magnus_1d
    END INTERFACE magnus

    INTERFACE ideal_gas_law_rho
       MODULE PROCEDURE ideal_gas_law_rho_0d
       MODULE PROCEDURE ideal_gas_law_rho_1d
    END INTERFACE ideal_gas_law_rho

    INTERFACE ideal_gas_law_rho_pt
       MODULE PROCEDURE ideal_gas_law_rho_pt_0d
       MODULE PROCEDURE ideal_gas_law_rho_pt_1d
    END INTERFACE ideal_gas_law_rho_pt

    INTERFACE exner_function
       MODULE PROCEDURE exner_function_0d
       MODULE PROCEDURE exner_function_1d
    END INTERFACE exner_function

    INTERFACE exner_function_invers
       MODULE PROCEDURE exner_function_invers_0d
       MODULE PROCEDURE exner_function_invers_1d
    END INTERFACE exner_function_invers

    INTERFACE barometric_formula
       MODULE PROCEDURE barometric_formula_0d
       MODULE PROCEDURE barometric_formula_1d
    END INTERFACE barometric_formula
!
!-- Public routines
    PUBLIC convert_utm_to_geographic

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Convert UTM coordinates into geographic latitude and longitude. Conversion
!> is based on the work of Krüger (1912) DOI: 10.2312/GFZ.b103-krueger28
!> and Karney (2013) DOI: 10.1007/s00190-012-0578-z
!> Based on a JavaScript of the geodesy function library written by chrisveness
!> https://github.com/chrisveness/geodesy
!------------------------------------------------------------------------------!
    SUBROUTINE convert_utm_to_geographic( crs, eutm, nutm, lon, lat )

       INTEGER(iwp) ::  j   !< loop index

       REAL(wp), INTENT(in)  ::  eutm !< easting (UTM)
       REAL(wp), INTENT(out) ::  lat  !< geographic latitude in degree
       REAL(wp), INTENT(out) ::  lon  !< geographic longitude in degree
       REAL(wp), INTENT(in)  ::  nutm !< northing (UTM)

       REAL(wp) ::  a           !< 2*pi*a is the circumference of a meridian
       REAL(wp) ::  cos_eta_s   !< cos(eta_s)
       REAL(wp) ::  delta_i     !<
       REAL(wp) ::  delta_tau_i !<
       REAL(wp) ::  e           !< eccentricity
       REAL(wp) ::  eta         !<
       REAL(wp) ::  eta_s       !<
       REAL(wp) ::  n           !< 3rd flattening
       REAL(wp) ::  n2          !< n^2
       REAL(wp) ::  n3          !< n^3
       REAL(wp) ::  n4          !< n^4
       REAL(wp) ::  n5          !< n^5
       REAL(wp) ::  n6          !< n^6
       REAL(wp) ::  nu          !<
       REAL(wp) ::  nu_s        !<
       REAL(wp) ::  sin_eta_s   !< sin(eta_s)
       REAL(wp) ::  sinh_nu_s   !< sinush(nu_s)
       REAL(wp) ::  tau_i       !<
       REAL(wp) ::  tau_i_s     !<
       REAL(wp) ::  tau_s       !<
       REAL(wp) ::  x           !< adjusted easting
       REAL(wp) ::  y           !< adjusted northing

       REAL(wp), DIMENSION(6) ::  beta !< 6th order Krüger expressions

       REAL(wp), DIMENSION(8), INTENT(in) ::  crs !< coordinate reference system, consists of
                                                  !< (/semi_major_axis,
                                                  !<   inverse_flattening,
                                                  !<   longitude_of_prime_meridian,
                                                  !<   longitude_of_central_meridian,
                                                  !<   scale_factor_at_central_meridian,
                                                  !<   latitude_of_projection_origin,
                                                  !<   false_easting,
                                                  !<   false_northing /)

       x = eutm - crs(7)  ! remove false easting
       y = nutm - crs(8)  ! remove false northing
!
!--    from Karney 2011 Eq 15-22, 36:
       e = SQRT( 1.0_wp / crs(2) * ( 2.0_wp - 1.0_wp / crs(2) ) )
       n = 1.0_wp / crs(2) / ( 2.0_wp - 1.0_wp / crs(2) )
       n2 = n * n
       n3 = n * n2
       n4 = n * n3
       n5 = n * n4
       n6 = n * n5

       a = crs(1) / ( 1.0_wp + n ) * ( 1.0_wp + 0.25_wp * n2       &
                                              + 0.015625_wp * n4   &
                                              + 3.90625E-3_wp * n6 )

       nu  = x / ( crs(5) * a )
       eta = y / ( crs(5) * a )

!--    According to Krüger (1912), eq. 26*
       beta(1) =        0.5_wp                  * n  &
               -        2.0_wp /         3.0_wp * n2 &
               +       37.0_wp /        96.0_wp * n3 &
               -        1.0_wp /       360.0_wp * n4 &
               -       81.0_wp /       512.0_wp * n5 &
               +    96199.0_wp /    604800.0_wp * n6

       beta(2) =        1.0_wp /        48.0_wp * n2 &
               +        1.0_wp /        15.0_wp * n3 &
               -      437.0_wp /      1440.0_wp * n4 &
               +       46.0_wp /       105.0_wp * n5 &
               -  1118711.0_wp /   3870720.0_wp * n6

       beta(3) =       17.0_wp /       480.0_wp * n3 &
               -       37.0_wp /       840.0_wp * n4 &
               -      209.0_wp /      4480.0_wp * n5 &
               +     5569.0_wp /     90720.0_wp * n6

       beta(4) =     4397.0_wp /    161280.0_wp * n4 &
               -       11.0_wp /       504.0_wp * n5 &
               -   830251.0_wp /   7257600.0_wp * n6

       beta(5) =     4583.0_wp /    161280.0_wp * n5 &
               -   108847.0_wp /   3991680.0_wp * n6

       beta(6) = 20648693.0_wp / 638668800.0_wp * n6

       eta_s = eta
       nu_s  = nu
       DO  j = 1, 6
         eta_s = eta_s - beta(j) * SIN(2.0_wp * j * eta) * COSH(2.0_wp * j * nu)
         nu_s  = nu_s  - beta(j) * COS(2.0_wp * j * eta) * SINH(2.0_wp * j * nu)
       ENDDO

       sinh_nu_s = SINH( nu_s )
       sin_eta_s = SIN( eta_s )
       cos_eta_s = COS( eta_s )

       tau_s = sin_eta_s / SQRT( sinh_nu_s**2 + cos_eta_s**2 )

       tau_i = tau_s
       delta_tau_i = 1.0_wp

       DO WHILE ( ABS( delta_tau_i ) > 1.0E-12_wp )

         delta_i = SINH( e * ATANH( e * tau_i / SQRT( 1.0_wp + tau_i**2 ) ) )

         tau_i_s = tau_i   * SQRT( 1.0_wp + delta_i**2 )  &
                  - delta_i * SQRT( 1.0_wp + tau_i**2 )

         delta_tau_i = ( tau_s - tau_i_s ) / SQRT( 1.0_wp + tau_i_s**2 )  &
                      * ( 1.0_wp + ( 1.0_wp - e**2 ) * tau_i**2 )          &
                      / ( ( 1.0_wp - e**2 ) * SQRT( 1.0_wp + tau_i**2 ) )

         tau_i = tau_i + delta_tau_i

       ENDDO

       lat = ATAN( tau_i ) / pi * 180.0_wp
       lon = ATAN2( sinh_nu_s, cos_eta_s ) / pi * 180.0_wp + crs(4)

    END SUBROUTINE convert_utm_to_geographic

!------------------------------------------------------------------------------!
! Description:
! ------------
!> This function computes the magnus formula (Press et al., 1992).
!> The magnus formula is needed to calculate the saturation vapor pressure
!------------------------------------------------------------------------------!
    FUNCTION magnus_0d( t )

       IMPLICIT NONE

       REAL(wp), INTENT(IN) ::  t  !< temperature (K)

       REAL(wp) ::  magnus_0d
!
!--    Saturation vapor pressure for a specific temperature:
       magnus_0d =  611.2_wp * EXP( 17.62_wp * ( t - degc_to_k ) /             &
                                               ( t - 29.65_wp  ) )

    END FUNCTION magnus_0d

!------------------------------------------------------------------------------!
! Description:
! ------------
!> This function computes the magnus formula (Press et al., 1992).
!> The magnus formula is needed to calculate the saturation vapor pressure
!------------------------------------------------------------------------------!
    FUNCTION magnus_1d( t )

       IMPLICIT NONE

       REAL(wp), INTENT(IN), DIMENSION(:) ::  t  !< temperature (K)

       REAL(wp), DIMENSION(size(t)) ::  magnus_1d
!
!--    Saturation vapor pressure for a specific temperature:
       magnus_1d =  611.2_wp * EXP( 17.62_wp * ( t - degc_to_k ) /             &
                                               ( t - 29.65_wp  ) )

    END FUNCTION magnus_1d

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the ideal gas law for scalar arguments.
!------------------------------------------------------------------------------!
    FUNCTION ideal_gas_law_rho_0d( p, t )

       IMPLICIT NONE

       REAL(wp), INTENT(IN) ::  p  !< pressure (Pa)
       REAL(wp), INTENT(IN) ::  t  !< temperature (K)

       REAL(wp) ::  ideal_gas_law_rho_0d
!
!--    compute density according to ideal gas law:
       ideal_gas_law_rho_0d = p / (r_d * t)

    END FUNCTION ideal_gas_law_rho_0d

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the ideal gas law for 1-D array arguments.
!------------------------------------------------------------------------------!
    FUNCTION ideal_gas_law_rho_1d( p, t )

       IMPLICIT NONE

       REAL(wp), INTENT(IN), DIMENSION(:) ::  p  !< pressure (Pa)
       REAL(wp), INTENT(IN), DIMENSION(:) ::  t  !< temperature (K)

       REAL(wp), DIMENSION(size(p)) ::  ideal_gas_law_rho_1d
!
!--    compute density according to ideal gas law:
       ideal_gas_law_rho_1d = p / (r_d * t)

    END FUNCTION ideal_gas_law_rho_1d

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the ideal gas law for scalar arguments.
!------------------------------------------------------------------------------!
    FUNCTION ideal_gas_law_rho_pt_0d( p, t )

       IMPLICIT NONE

       REAL(wp), INTENT(IN) ::  p  !< pressure (Pa)
       REAL(wp), INTENT(IN) ::  t  !< temperature (K)

       REAL(wp) ::  ideal_gas_law_rho_pt_0d
!
!--    compute density according to ideal gas law:
       ideal_gas_law_rho_pt_0d = p / (r_d * exner_function(p) * t)

    END FUNCTION ideal_gas_law_rho_pt_0d

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the ideal gas law for 1-D array arguments.
!------------------------------------------------------------------------------!
    FUNCTION ideal_gas_law_rho_pt_1d( p, t )

       IMPLICIT NONE

       REAL(wp), INTENT(IN), DIMENSION(:) ::  p  !< pressure (Pa)
       REAL(wp), INTENT(IN), DIMENSION(:) ::  t  !< temperature (K)

       REAL(wp), DIMENSION(size(p)) ::  ideal_gas_law_rho_pt_1d
!
!--    compute density according to ideal gas law:
       ideal_gas_law_rho_pt_1d = p / (r_d * exner_function(p) * t)

    END FUNCTION ideal_gas_law_rho_pt_1d

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the exner function for scalar arguments.
!------------------------------------------------------------------------------!
    FUNCTION exner_function_0d( p )

       IMPLICIT NONE

       REAL(wp), INTENT(IN) ::  p    !< pressure (Pa)

       REAL(wp) ::  exner_function_0d
!
!--    compute exner function:
       exner_function_0d = ( p / p_0 )**( rd_d_cp )

    END FUNCTION exner_function_0d

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the exner function for 1-D array arguments.
!------------------------------------------------------------------------------!
    FUNCTION exner_function_1d( p )

       IMPLICIT NONE

       REAL(wp), INTENT(IN), DIMENSION(:) ::  p  !< pressure (Pa)

       REAL(wp), DIMENSION(size(p)) ::  exner_function_1d
!
!--    compute exner function:
       exner_function_1d = ( p / p_0 )**( rd_d_cp )

    END FUNCTION exner_function_1d

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the exner function for scalar arguments.
!------------------------------------------------------------------------------!
    FUNCTION exner_function_invers_0d( p )

       IMPLICIT NONE

       REAL(wp), INTENT(IN) ::  p    !< pressure (Pa)

       REAL(wp) ::  exner_function_invers_0d
!
!--    compute exner function:
       exner_function_invers_0d = ( p_0 / p )**( rd_d_cp )

    END FUNCTION exner_function_invers_0d

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the exner function for 1-D array arguments.
!------------------------------------------------------------------------------!
    FUNCTION exner_function_invers_1d( p )

       IMPLICIT NONE

       REAL(wp), INTENT(IN), DIMENSION(:) ::  p  !< pressure (Pa)

       REAL(wp), DIMENSION(size(p)) ::  exner_function_invers_1d
!
!--    compute exner function:
       exner_function_invers_1d = ( p_0 / p )**( rd_d_cp )

    END FUNCTION exner_function_invers_1d

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the barometric formula for scalar arguments. The calculation is 
!> based on the assumption of a polytropic atmosphere and neutral 
!> stratification, where the temperature lapse rate is g/cp. 
!------------------------------------------------------------------------------!
    FUNCTION barometric_formula_0d( z, t_0, p_0)

       IMPLICIT NONE

       REAL(wp), INTENT(IN) ::  z    !< height (m)
       REAL(wp), INTENT(IN) ::  t_0  !< temperature reference state (K)
       REAL(wp), INTENT(IN) ::  p_0  !< surface pressure (Pa)

       REAL(wp) ::  barometric_formula_0d
!
!--    compute barometric formula:
       barometric_formula_0d =  p_0 * ( (t_0 - g_d_cp * z) / t_0 )**( cp_d_rd )

    END FUNCTION barometric_formula_0d

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the barometric formula for 1-D array arguments. The calculation is 
!> based on the assumption of a polytropic atmosphere and neutral 
!> stratification, where the temperature lapse rate is g/cp.
!------------------------------------------------------------------------------!
    FUNCTION barometric_formula_1d( z, t_0, p_0)

       IMPLICIT NONE

       REAL(wp), INTENT(IN), DIMENSION(:) ::  z  !< height (m)
       REAL(wp), INTENT(IN) ::  t_0              !< temperature reference state (K)
       REAL(wp), INTENT(IN) ::  p_0              !< surface pressure (Pa)

       REAL(wp), DIMENSION(size(z)) ::  barometric_formula_1d
!
!--    compute barometric formula:
       barometric_formula_1d =  p_0 * ( (t_0 - g_d_cp * z) / t_0 )**( cp_d_rd )

    END FUNCTION barometric_formula_1d

 END MODULE basic_constants_and_equations_mod
