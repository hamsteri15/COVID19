!> @file mod_kinds.f90
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
! $Id: mod_kinds.f90 4360 2020-01-07 11:25:50Z suehring $
! Corrected "Former revisions" section
! 
! 4000 2019-05-24 07:20:44Z raasch
! preprocessor switch added for choosing the real precision
! 
! 3655 2019-01-07 16:51:22Z knoop
! Corrected "Former revisions" section
! 
! 1306 2014-03-13 14:30:59Z raasch
! Initial revision
!
! Description:
! ------------
!> Standard kind definitions
!> wp (working precision) and iwp (integer working precision) are the kinds
!> used by default in all variable declarations.
!> By default, PALM is using wp = dp (64bit), and iwp = isp (32bit).
!> If you like to switch to other precision, then please set wp/iwp
!> appropriately by assigning other kinds below.
!------------------------------------------------------------------------------!
 MODULE kinds
 

    IMPLICIT NONE

!
!-- Floating point kinds
    INTEGER, PARAMETER ::  sp = 4           !< single precision (32 bit)
    INTEGER, PARAMETER ::  dp = 8           !< double precision (64 bit)

!
!-- Integer kinds
    INTEGER, PARAMETER ::  isp = SELECTED_INT_KIND(  9 )   !< single precision (32 bit)
    INTEGER, PARAMETER ::  idp = SELECTED_INT_KIND( 14 )   !< double precision (64 bit)

!
!-- Set kinds to be used as defaults
#if defined( __single_precision )
    INTEGER, PARAMETER ::   wp =  sp          !< default real kind
#else
    INTEGER, PARAMETER ::   wp =  dp          !< default real kind
#endif
    INTEGER, PARAMETER ::  iwp = isp          !< default integer kind

    SAVE

 END MODULE kinds
