!> @file palm_date_time_mod.f90
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
! $Id: palm_date_time_mod.f90 4360 2020-01-07 11:25:50Z suehring $
! Add days of northward- and southward equinox
! 
! 4227 2019-09-10 18:04:34Z gronemeier
! Complete rework of module date_and_time_mod:
!  - renamed module to prevent confusion with
!    FORTRAN Standard routine date_and_time
!  - introduce date_time_type
!  - add set_reference_date_time
!  - add get_date_time
!  - capsule whole calculation of date and time variables within this routine
!  - removed all variables/routines not belonging to this module
!
!
! Authors:
! --------
!> @author Tobias Gronemeier (LUH)
!
! Description:
! ------------
!> This routine calculates all needed information on date and time used by
!> other modules
!>
!> @todo Consider leap seconds
!> @note Time_zone only supports full-hour time zones, i.e., time zones like
!>       Australian Central Standard Time (UTC+9.5) are not possible
!------------------------------------------------------------------------------!
 MODULE palm_date_time_mod

    USE control_parameters,                                                    &
         ONLY:  message_string

    USE kinds

    IMPLICIT NONE

!
!-- Parameter Definition
    INTEGER(iwp), PARAMETER ::  date_time_str_len  = 23_iwp                                 !< length of date_time strings
    INTEGER(iwp), PARAMETER ::  days_per_week      = 7_iwp                                  !< days in a week
    INTEGER(iwp), PARAMETER ::  hours_per_day      = 24_iwp                                 !< hours in a day
    INTEGER(iwp), PARAMETER ::  minutes_per_hour   = 60_iwp                                 !< minutes in an hour
    INTEGER(iwp), PARAMETER ::  months_per_year    = 12_iwp                                 !< months in a year
!
!-- Definition of mean northward and southward equinox (summer and winter half year) 
!-- in days of year. For simplicity, March 21 and September 21 is assumed.
    INTEGER(iwp), PARAMETER ::  northward_equinox  = 80_iwp
    INTEGER(iwp), PARAMETER ::  southward_equinox  = 264_iwp

    REAL(wp),     PARAMETER ::  seconds_per_minute = 60.0_wp                                !< seconds in a minute
    REAL(wp),     PARAMETER ::  seconds_per_hour   = seconds_per_minute * minutes_per_hour  !< seconds in an hour
    REAL(wp),     PARAMETER ::  seconds_per_day    = seconds_per_hour * hours_per_day       !< seconds in a day

    CHARACTER(LEN=3), DIMENSION(days_per_week), PARAMETER ::  &
       weekdays = (/"Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"/)                       !< names of weekdays

    INTEGER, DIMENSION(months_per_year), PARAMETER ::  &
       days_per_month_noleapyear = (/31,28,31,30,31,30,31,31,30,31,30,31/)                  !< days for each month (no leap year)

    INTEGER, DIMENSION(months_per_year), PARAMETER ::  &
       days_per_month_leapyear = (/31,29,31,30,31,30,31,31,30,31,30,31/)                    !< days for each month (leap year)

    INTEGER, DIMENSION(121), PARAMETER ::  leap_years = &                                   !< list of leap years
        (/1804_iwp, 1808_iwp, 1812_iwp, 1816_iwp, 1820_iwp, 1824_iwp, 1828_iwp, 1832_iwp, &
          1836_iwp, 1840_iwp, 1844_iwp, 1848_iwp, 1852_iwp, 1856_iwp, 1860_iwp, 1864_iwp, &
          1868_iwp, 1872_iwp, 1876_iwp, 1880_iwp, 1884_iwp, 1888_iwp, 1892_iwp, 1896_iwp, &
          1904_iwp, 1908_iwp, 1912_iwp, 1916_iwp, 1920_iwp, 1924_iwp, 1928_iwp, 1932_iwp, &
          1936_iwp, 1940_iwp, 1944_iwp, 1948_iwp, 1952_iwp, 1956_iwp, 1960_iwp, 1964_iwp, &
          1968_iwp, 1972_iwp, 1976_iwp, 1980_iwp, 1984_iwp, 1988_iwp, 1992_iwp, 1996_iwp, &
          2000_iwp, 2004_iwp, 2008_iwp, 2012_iwp, 2016_iwp, 2020_iwp, 2024_iwp, 2028_iwp, &
          2032_iwp, 2036_iwp, 2040_iwp, 2044_iwp, 2048_iwp, 2052_iwp, 2056_iwp, 2060_iwp, &
          2064_iwp, 2068_iwp, 2072_iwp, 2076_iwp, 2080_iwp, 2084_iwp, 2088_iwp, 2092_iwp, &
          2096_iwp, 2104_iwp, 2108_iwp, 2112_iwp, 2116_iwp, 2120_iwp, 2124_iwp, 2128_iwp, &
          2132_iwp, 2136_iwp, 2140_iwp, 2144_iwp, 2148_iwp, 2152_iwp, 2156_iwp, 2160_iwp, &
          2164_iwp, 2168_iwp, 2172_iwp, 2176_iwp, 2180_iwp, 2184_iwp, 2188_iwp, 2192_iwp, &
          2196_iwp, 2204_iwp, 2208_iwp, 2212_iwp, 2216_iwp, 2220_iwp, 2224_iwp, 2228_iwp, &
          2232_iwp, 2236_iwp, 2240_iwp, 2244_iwp, 2248_iwp, 2252_iwp, 2256_iwp, 2260_iwp, &
          2264_iwp, 2268_iwp, 2272_iwp, 2276_iwp, 2280_iwp, 2284_iwp, 2288_iwp, 2292_iwp, &
          2296_iwp /)

!
!-- Type Definition
    TYPE date_time_type
       INTEGER(iwp)                        ::  year           = -HUGE(0_iwp)               !< year
       INTEGER(iwp)                        ::  month          = -HUGE(0_iwp)               !< month of year
       INTEGER(iwp)                        ::  day            = -HUGE(0_iwp)               !< day of month
       INTEGER(iwp)                        ::  hour           = -HUGE(0_iwp)               !< hour of day
       INTEGER(iwp)                        ::  minute         = -HUGE(0_iwp)               !< minute of hour
       INTEGER(iwp)                        ::  zone           = -HUGE(0_iwp)               !< time zone

       REAL(wp)                            ::  second         = -HUGE(0.0_wp)              !< second of minute
       REAL(wp)                            ::  second_of_year = -HUGE(0.0_wp)              !< second of year

       INTEGER(iwp)                        ::  days_per_year  = -HUGE(0_iwp)               !< days within a year

       INTEGER, DIMENSION(months_per_year) ::  days_per_month = days_per_month_noleapyear  !< list of total days per month
    END TYPE date_time_type

!
!-- Variable Declaration
    LOGICAL              ::  reference_date_time_is_set = .FALSE.  !< true if reference_date_time is set

    TYPE(date_time_type) ::  reference_date_time                   !< reference date-time

    SAVE

    PRIVATE
!
!-- Set reference date and time
    INTERFACE set_reference_date_time
        MODULE PROCEDURE set_reference_date_time
    END INTERFACE set_reference_date_time
!
!-- Return date and time information
    INTERFACE get_date_time
       MODULE PROCEDURE get_date_time
    END INTERFACE get_date_time
!
!-- Public Interfaces
    PUBLIC &
       get_date_time, &
       set_reference_date_time
!
!-- Public variables
    PUBLIC                 &
       date_time_str_len,  &
       days_per_week,      &
       hours_per_day,      &
       minutes_per_hour,   &
       months_per_year,    &
       northward_equinox,  &
       seconds_per_minute, &
       seconds_per_hour,   &
       seconds_per_day,    &
       southward_equinox,  &
       weekdays

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set reference date-time.
!> Only a single call is allowed to this routine during execution.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE set_reference_date_time( date_time_str )

    CHARACTER(LEN=date_time_str_len), INTENT(IN) ::  date_time_str  !< string containing date-time information

!
!-- Check if date and time are already set
    IF ( reference_date_time_is_set )  THEN
       !> @note This error should never be observed by a user.
       !>       It can only appear if the code was modified.
       WRITE( message_string, * ) 'Multiple calls to set_reference_date_time detected.&' //  &
                                  'The reference date-time must be set only once.'
       CALL message( 'set_reference_date_time', 'PA0680', 2, 2, 0, 6, 0 )
       RETURN

    ELSE

       reference_date_time = convert_string_to_date_time( date_time_str )

       reference_date_time_is_set = .TRUE.

    ENDIF

 END SUBROUTINE set_reference_date_time


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Return requested date-time information of the reference time + time_since_reference.
!> An alternative reference date-time string can be specified via 'reference_date_time_str'.
!> Call to this routine is only possible if a reference time is either specified in the call itself
!> via 'reference_date_time_str' or previously set by calling routine 'set_reference_date_time'.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE get_date_time( time_since_reference, reference_date_time_str,    &
                           year, month, day, hour, minute, second, zone,     &
                           second_of_day, second_of_year,                    &
                           day_of_year, day_of_week, weekday, date_time_str, &
                           days_per_month, days_per_year                     )

    CHARACTER(LEN=date_time_str_len), INTENT(OUT), OPTIONAL ::  date_time_str            !< date-time as string
    CHARACTER(LEN=1)                                        ::  plus_minus               !< either '+' or '-'
    CHARACTER(LEN=date_time_str_len), INTENT(IN),  OPTIONAL ::  reference_date_time_str  !< alternative reference date-time
    CHARACTER(LEN=3),                 INTENT(OUT), OPTIONAL ::  weekday                  !< weekday

    INTEGER(iwp),                             INTENT(OUT), OPTIONAL ::  day              !< day of month
    INTEGER(iwp),                             INTENT(OUT), OPTIONAL ::  day_of_week      !< number of weekday
    INTEGER(iwp),                             INTENT(OUT), OPTIONAL ::  day_of_year      !< day of the year
    INTEGER(iwp),                             INTENT(OUT), OPTIONAL ::  hour             !< hour of day
    INTEGER(iwp),                             INTENT(OUT), OPTIONAL ::  minute           !< minute of hour
    INTEGER(iwp),                             INTENT(OUT), OPTIONAL ::  month            !< month of year
    INTEGER(iwp),                             INTENT(OUT), OPTIONAL ::  year             !< year
    INTEGER(iwp),                             INTENT(OUT), OPTIONAL ::  zone             !< time zone
    INTEGER(iwp),                             INTENT(OUT), OPTIONAL ::  days_per_year    !< days per year
    INTEGER(iwp), DIMENSION(months_per_year), INTENT(OUT), OPTIONAL ::  days_per_month   !< days per year

    REAL(wp), INTENT(OUT), OPTIONAL ::  second                        !< second of minute
    REAL(wp), INTENT(OUT), OPTIONAL ::  second_of_day                 !< second of day
    REAL(wp), INTENT(OUT), OPTIONAL ::  second_of_year                !< second of year
    REAL(wp), INTENT(IN)            ::  time_since_reference          !< seconds between reference time and current time

    TYPE(date_time_type)            ::  date_time                     !< date-time which to return
    TYPE(date_time_type)            ::  internal_reference_date_time  !< internal reference date-time

!
!-- Check if a reference date-time is given
    IF ( .NOT. reference_date_time_is_set  .AND.  .NOT. PRESENT( reference_date_time_str ) )  THEN
       !> @note This error should never be observed by a user.
       !>       It can only appear if the code was modified.
       WRITE( message_string, * ) 'No reference date-time defined. '//                   &
                                  'Returning date-time information is not possible. ' // &
                                  'Either specify reference_date_time_str ' //           &
                                  'or set a reference via set_reference_date_time.'
       CALL message( 'get_date_time', 'PA0677', 2, 2, 0, 6, 0 )
       RETURN
    ENDIF
!
!-- Set internal reference date-time
    IF ( PRESENT( reference_date_time_str ) )  THEN
       internal_reference_date_time = convert_string_to_date_time( reference_date_time_str )
    ELSE
       internal_reference_date_time = reference_date_time
    ENDIF
!
!-- Add time to reference time
    date_time = add_date_time( time_since_reference, internal_reference_date_time )
!
!-- Set requested return values
    IF ( PRESENT( year           ) )  year           = date_time%year
    IF ( PRESENT( month          ) )  month          = date_time%month
    IF ( PRESENT( day            ) )  day            = date_time%day
    IF ( PRESENT( hour           ) )  hour           = date_time%hour
    IF ( PRESENT( minute         ) )  minute         = date_time%minute
    IF ( PRESENT( second         ) )  second         = date_time%second
    IF ( PRESENT( zone           ) )  zone           = date_time%zone
    IF ( PRESENT( second_of_year ) )  second_of_year = date_time%second_of_year
    IF ( PRESENT( second_of_day  ) )  second_of_day  = get_second_of_day( date_time )
    IF ( PRESENT( day_of_year    ) )  day_of_year    = get_day_of_year( date_time )
    IF ( PRESENT( day_of_week    ) )  day_of_week    = get_day_of_week( date_time )
    IF ( PRESENT( weekday        ) )  weekday        = weekdays( get_day_of_week( date_time ) )
    IF ( PRESENT( days_per_month ) )  days_per_month = date_time%days_per_month
    IF ( PRESENT( days_per_year  ) )  days_per_year  = date_time%days_per_year

    IF ( PRESENT( date_time_str ) )  THEN
       IF ( date_time%zone < 0_iwp )  THEN
          plus_minus = '-'
       ELSE
          plus_minus = '+'
       ENDIF
       WRITE( UNIT = date_time_str,                                                 &
              FMT = '(I4,"-",I2.2,"-",I2.2,1X,I2.2,":",I2.2,":",I2.2,1X,A1,I2.2)' ) &
          date_time%year, date_time%month, date_time%day,                           &
          date_time%hour, date_time%minute, INT( date_time%second ),                &
          plus_minus, ABS( date_time%zone )
    ENDIF

 END SUBROUTINE get_date_time


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Convert a date-time string into a date_time object.
!--------------------------------------------------------------------------------------------------!
 FUNCTION convert_string_to_date_time( date_time_str ) RESULT( date_time )

    CHARACTER(LEN=date_time_str_len), INTENT(IN) ::  date_time_str  !< date-time as string

    INTEGER(iwp)                                 ::  read_status    !< returned status of read command

    TYPE(date_time_type)                         ::  date_time      !< requested date-time object

!
!-- Decompose string to date-time information
    READ( UNIT = date_time_str( 1: 4), IOSTAT = read_status, FMT = '(I4)'   )  date_time%year
    IF ( read_status == 0 )  &
       READ( UNIT = date_time_str( 6: 7), IOSTAT = read_status, FMT = '(I2)'   )  date_time%month
    IF ( read_status == 0 )  &
       READ( UNIT = date_time_str( 9:10), IOSTAT = read_status, FMT = '(I2)'   )  date_time%day
    IF ( read_status == 0 )  &
       READ( UNIT = date_time_str(12:13), IOSTAT = read_status, FMT = '(I2)'   )  date_time%hour
    IF ( read_status == 0 )  &
       READ( UNIT = date_time_str(15:16), IOSTAT = read_status, FMT = '(I2)'   )  date_time%minute
    IF ( read_status == 0 )  &
       READ( UNIT = date_time_str(18:19), IOSTAT = read_status, FMT = '(F2.0)' )  date_time%second
    IF ( read_status == 0 )  &
       READ( UNIT = date_time_str(21:23), IOSTAT = read_status, FMT = '(I3)'   )  date_time%zone

    IF ( read_status /= 0 )  THEN
       WRITE( message_string, * ) 'Error while converting date-time string. ' //  &
                                  'Please check format of date-time string: "' // &
                                  TRIM( date_time_str ) // '". ' //               &
                                  'Format must be "YYYY-MM-DD hh:mm:ss ZZZ".'
       CALL message( 'convert_string_to_date_time', 'PA0678', 2, 2, 0, 6, 0 )
       RETURN
    ENDIF

    date_time = update_leapyear_setting( date_time )

    IF ( check_date( date_time, date_time_str ) == 0 )  THEN

       date_time%second_of_year = get_second_of_year( date_time )

!
!--    Shift time to UTC and set zone to UTC
       date_time = add_date_time( REAL( -1 * date_time%zone, KIND = wp ) &
                                  * seconds_per_hour,                    &
                                  date_time )
       date_time%zone = 0_iwp
    ENDIF

 END FUNCTION convert_string_to_date_time


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Add time increment (in seconds) to given date-time and return shifted date-time
!--------------------------------------------------------------------------------------------------!
 FUNCTION add_date_time( inc_seconds, date_time_base ) RESULT( date_time )

    INTEGER(iwp)                      ::  i                 !< loop index

    REAL(wp)                          ::  seconds_left      !< seconds which must still be added to new date-time
    REAL(wp)                          ::  seconds_per_year  !< number of seconds in a year

    REAL(wp),             INTENT(IN)  ::  inc_seconds       !< seconds to be added to date-time

    TYPE(date_time_type)              ::  date_time         !< shifted date-time
    TYPE(date_time_type), INTENT(IN)  ::  date_time_base    !< date-time to be shifted

!
!-- Define some parameters
    date_time = date_time_base
    seconds_per_year = REAL( date_time%days_per_year,  KIND = wp ) * seconds_per_day
!
!-- Shift time
    date_time%second_of_year = date_time%second_of_year + inc_seconds
!
!-- Check if year changes
!-- First, if year is reduced
    DO WHILE ( date_time%second_of_year < 0.0_wp )
       date_time%year = date_time%year - 1_iwp
       date_time = update_leapyear_setting( date_time )
       seconds_per_year = REAL( date_time%days_per_year * seconds_per_day, KIND = wp )
       date_time%second_of_year = date_time%second_of_year + seconds_per_year
    ENDDO
!
!-- Now, if year is increased
    DO WHILE ( date_time%second_of_year > seconds_per_year )
       date_time%year = date_time%year + 1_iwp
       date_time = update_leapyear_setting( date_time )
       date_time%second_of_year = date_time%second_of_year - seconds_per_year
       seconds_per_year = REAL( date_time%days_per_year * seconds_per_day, KIND = wp )
    ENDDO
!
!-- Based on updated year and second_of_year, update month, day, hour, minute, second
    DO  i = 1, months_per_year
       IF ( date_time%second_of_year < SUM( date_time%days_per_month(1:i) ) * seconds_per_day ) &
       THEN
          date_time%month  = i
          seconds_left     = date_time%second_of_year                                &
                           - REAL( SUM( date_time%days_per_month(1:i-1) ), KIND=wp ) &
                           * seconds_per_day
          date_time%day    = INT( seconds_left / seconds_per_day, KIND = iwp ) + 1_iwp
          seconds_left     = seconds_left &
                           - REAL( date_time%day - 1_iwp, KIND = wp ) * seconds_per_day
          date_time%hour   = INT( seconds_left / seconds_per_hour, KIND = iwp )
          seconds_left     = seconds_left - REAL( date_time%hour, KIND = wp ) * seconds_per_hour
          date_time%minute = INT( seconds_left / seconds_per_minute, KIND = iwp )
          date_time%second = seconds_left - REAL( date_time%minute, KIND = wp ) * seconds_per_minute
          EXIT
       ENDIF
    ENDDO

 END FUNCTION add_date_time


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Return day of year of given date.
!--------------------------------------------------------------------------------------------------!
 FUNCTION get_day_of_year( date_time ) RESULT( day_of_year )

    INTEGER(iwp)                     ::  day_of_year         !< day of the year

    TYPE(date_time_type), INTENT(IN) ::  date_time           !< date of which to calculate day of year
    TYPE(date_time_type)             ::  date_time_internal  !< internal copy of input date-time


    date_time_internal = update_leapyear_setting( date_time )

    day_of_year = date_time_internal%day &
                + SUM( date_time_internal%days_per_month(:date_time_internal%month-1) )

 END FUNCTION get_day_of_year


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Return second of day of given time.
!--------------------------------------------------------------------------------------------------!
 FUNCTION get_second_of_day( date_time ) RESULT( second_of_day )

    REAL(wp)                         ::  second_of_day  !< second of the day

    TYPE(date_time_type), INTENT(IN) ::  date_time      !< date of which to calculate second of the day


    second_of_day = date_time%second                                                            &
                  + REAL( ( date_time%hour * minutes_per_hour ) + date_time%minute, KIND = wp ) &
                  * seconds_per_minute

 END FUNCTION get_second_of_day


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Return second of year of given date-time.
!--------------------------------------------------------------------------------------------------!
 FUNCTION get_second_of_year( date_time ) RESULT( second_of_year )

    REAL(wp)                         ::  second_of_year  !< second of the year

    TYPE(date_time_type), INTENT(IN) ::  date_time       !< date of which to calculate second of the year


    second_of_year = get_second_of_day( date_time ) &
                   + REAL( get_day_of_year( date_time ) - 1_iwp, KIND = wp ) * seconds_per_day

 END FUNCTION get_second_of_year


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Return index of the day of the week of the given date-time.
!--------------------------------------------------------------------------------------------------!
 FUNCTION get_day_of_week( date_time_in ) RESULT( day_of_week )

    INTEGER(iwp)                     ::  date_time_internal_reference_day_of_week  !< day of week of reference date
    INTEGER(iwp)                     ::  day_difference                            !< day between given date and reference
    INTEGER(iwp)                     ::  day_of_week                               !< day of the week

    TYPE(date_time_type), INTENT(IN) ::  date_time_in                              !< date of which to get the weekday
    TYPE(date_time_type)             ::  date_time_internal                        !< internal date-time

!
!-- Define reference date from which on the current day of week can be determined
    date_time_internal%year                  = 2000_iwp
    date_time_internal%month                 = 1_iwp
    date_time_internal%day                   = 1_iwp
    date_time_internal_reference_day_of_week = 6_iwp

!
!-- First, get the difference if both dates would be in the same year
    day_difference = get_day_of_year( date_time_in ) - get_day_of_year( date_time_internal )
!
!-- Now, shift the year and add the corresponding number of days to the difference
    IF ( date_time_internal%year < date_time_in%year )  THEN

       DO WHILE ( date_time_internal%year /= date_time_in%year )

          date_time_internal = update_leapyear_setting( date_time_internal )
          day_difference = day_difference + date_time_internal%days_per_year

          date_time_internal%year = date_time_internal%year + 1_iwp

       ENDDO

    ELSEIF ( date_time_internal%year > date_time_in%year )  THEN

       DO WHILE ( date_time_internal%year /= date_time_in%year )

          date_time_internal%year = date_time_internal%year - 1_iwp

          date_time_internal = update_leapyear_setting( date_time_internal )
          day_difference = day_difference - date_time_internal%days_per_year

       ENDDO

    ENDIF
!
!-- Remove full weeks of day_difference and shift day_of_week of reference by
!-- remaining days
    day_of_week = date_time_internal_reference_day_of_week + MOD( day_difference, days_per_week )

    IF ( day_of_week > days_per_week )  THEN
!
!--    Shift index again if it is next week (above days_per_week)...
       day_of_week = day_of_week - days_per_week
    ELSEIF ( day_of_week <= 0_iwp )  THEN
!
!--    ...or if it is last week (below 1)
       day_of_week = day_of_week + days_per_week
    ENDIF

 END FUNCTION get_day_of_week


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check if given year is a leap year and update days per month accordingly.
!--------------------------------------------------------------------------------------------------!
 FUNCTION update_leapyear_setting( date_time_in ) RESULT( date_time_out )

    TYPE(date_time_type), INTENT(IN) ::  date_time_in   !< input date-time
    TYPE(date_time_type)             ::  date_time_out  !< return date-time


    date_time_out = date_time_in

    IF ( ANY( date_time_in%year == leap_years ) )  THEN
       date_time_out%days_per_month = days_per_month_leapyear
    ELSE
      date_time_out%days_per_month = days_per_month_noleapyear
    ENDIF

    date_time_out%days_per_year = SUM( date_time_out%days_per_month )

 END FUNCTION update_leapyear_setting


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check if given date and time are valid. Returns 0 if all checks are passed.
!> @todo Revise error message. ATM, gives only last errorneous value even if
!>       multiple values violate the bounds.
!--------------------------------------------------------------------------------------------------!
 FUNCTION check_date( date_time, date_time_str ) RESULT( error_code )


    CHARACTER(LEN=6), DIMENSION(7), PARAMETER ::  error_str_list =  &  !< string used for error message
       (/'year  ', 'month ', 'day   ', 'hour  ', 'minute', 'second', 'zone  '/)

    CHARACTER(LEN=date_time_str_len), INTENT(IN) ::  date_time_str     !< date-time as string

    INTEGER(iwp)                                 ::  error_code        !< error code

    TYPE(date_time_type),             INTENT(IN) ::  date_time         !< date-time to be checked


    error_code = 0
!
!-- Check date
    IF ( date_time%month < 1_iwp  .OR.     &
         date_time%month > months_per_year )  THEN
       error_code = 2
    ELSE
       IF ( date_time%day < 1_iwp  .OR.                               &
            date_time%day > date_time%days_per_month(date_time%month) )  THEN
          error_code = 3
       ENDIF
    ENDIF
!
!-- Check time
    IF ( date_time%hour < 0_iwp  .OR.   &
         date_time%hour > hours_per_day )  THEN
       error_code = 4
    ELSE
        IF ( date_time%minute < 0_iwp  .OR.      &
             date_time%minute > minutes_per_hour )  THEN
           error_code = 5
        ELSE
           IF ( date_time%second < 0.0_wp  .OR.        &
                date_time%second >= seconds_per_minute )  THEN
              error_code = 6
           ENDIF
        ENDIF
    ENDIF
!
!-- Check time zone
!-- Bounds defined by maximum and minimum time zone present on earth
    IF ( date_time%zone < -12_iwp  .OR.  &
         date_time%zone > 14_iwp )  THEN
       error_code = 7
    ENDIF
!
!-- Raise error if any check is marked invalid
    IF ( error_code /= 0 )  THEN
       WRITE( message_string, * ) 'Date-time values out of bounds: "' //                    &
                                  TRIM( error_str_list(error_code) ) //                     &
                                  '" is out of bounds. Please check date-time string: "' // &
                                  TRIM( date_time_str ) // '". ' //                         &
                                  'Format must be "YYYY-MM-DD hh:mm:ss ZZZ".'
       CALL message( 'check_date', 'PA0679', 2, 2, 0, 6, 0 )
       RETURN
    ENDIF

 END FUNCTION check_date

 END MODULE palm_date_time_mod
