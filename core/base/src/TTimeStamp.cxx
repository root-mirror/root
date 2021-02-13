// @(#)root/base:$Id$
// Author: R. Hatcher   30/9/2001

/*************************************************************************
 * Copyright (C) 1995-2002, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TTimeStamp
\ingroup Base

The TTimeStamp encapsulates seconds and ns since EPOCH

This extends (and isolates) struct timespec
~~~ {.cpp}
   struct timespec
      {
         time_t   tv_sec;   // seconds
         long     tv_nsec;  // nanoseconds
      }
   time_t seconds is relative to Jan 1, 1970 00:00:00 UTC
~~~
No accounting of leap seconds is made.

Due to ROOT/CINT limitations TTimeStamp does not explicitly
hold a timespec struct; attempting to do so means the Streamer
must be hand written.  Instead we have chosen to simply contain
similar fields within the private area of this class.

NOTE: the use of time_t (and its default implementation as a 32 int)
      implies overflow conditions occurs somewhere around
      `Jan 18, 19:14:07, 2038`.
      If this experiment is still going when it becomes significant
      someone will have to deal with it.
*/

#include "TTimeStamp.h"
#include "TString.h"
#include "TError.h"
#include <iostream>
#ifdef R__WIN32
#include "Windows4Root.h"
#else
#include <unistd.h>
#include <sys/time.h>
#endif
#include "TVirtualMutex.h"

ClassImp(TTimeStamp);


TVirtualMutex *gTimeMutex = 0; // local mutex

////////////////////////////////////////////////////////////////////////////////
/// Write time stamp to std::ostream.

std::ostream& operator<<(std::ostream& os, const TTimeStamp& ts)
{
   if (os.good()) {
      if (os.tie()) os.tie()->flush(); // instead of opfx
      os << ts.AsString("c");
   }
   // instead of os.osfx()
   if (os.flags() & std::ios::unitbuf) os.flush();
   return os;
}

////////////////////////////////////////////////////////////////////////////////
/// Read time stamp from TBuffer.

TBuffer &operator>>(TBuffer &buf, TTimeStamp &ts)
{
   ts.Streamer(buf);
   return buf;
}

////////////////////////////////////////////////////////////////////////////////
/// Write time stamp to TBuffer.

TBuffer &operator<<(TBuffer &buf, const TTimeStamp &ts)
{
   ((TTimeStamp&)ts).Streamer(buf);
   return buf;
}

////////////////////////////////////////////////////////////////////////////////
/// Default ctor. Create a TTimeStamp and set it to the current time
/// (as best possible). The nanosecond part is faked so that subsequent
/// calls simply add 1 to ensure that sequential calls are distinct
/// (and sortable).

TTimeStamp::TTimeStamp()
{
   Set();
}

////////////////////////////////////////////////////////////////////////////////
/// Create a TTimeStamp and set it to the specified year, month,
/// day, time, hour, minute, second and nanosec.
/// If !isUTC then it is assumed to be the standard local time zone.
///
/// If local time is PST then one can use
/// ~~~ {.cpp}
///    TTimeStamp(year,month,day,hour,min,sec,nsec,kFALSE,0);
/// ~~~
/// or
/// ~~~ {.cpp}
///    Int_t secOffset = 8*60*60;
///    TTimeStamp timeStamp(year,month,day,hour,min,sec,nsec,kTRUE,8*60*60);
/// ~~~

TTimeStamp::TTimeStamp(UInt_t year, UInt_t month,
                       UInt_t day,  UInt_t hour,
                       UInt_t min,  UInt_t sec,
                       UInt_t nsec,
                       Bool_t isUTC, Int_t secOffset)
{
   Set(year, month, day, hour, min, sec, nsec, isUTC, secOffset);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a TTimeStamp and set it to the specified date, time, nanosec.
/// If !isUTC then it is assumed to be the standard local time zone.
///
/// \warning Watch out! C++ overload resolution often chooses the constructor
/// `TTimeStamp(UInt_t tloc, Bool_t isUTC, Int_t secOffset, Bool_t dosDate)`
/// instead of this one. Your best bet is to explicitly pass UInt_t values instead
/// of Int_t values. When calling with integer literals, pass for instance
/// ~~~ {.cpp}
/// TTimeStamp timeStamp(20150610u,80448u,0u)
/// ~~~
/// to disambiguate.

TTimeStamp::TTimeStamp(UInt_t date, UInt_t time,
                       UInt_t nsec,
                       Bool_t isUTC, Int_t secOffset)
{
   Set(date, time, nsec, isUTC, secOffset);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a TTimeStamp and set it to tloc which must be a time_t value
/// returned by time(). This value is the number of seconds since the EPOCH
/// (i.e. 00:00:00 on Jan 1m 1970). If dosDate is true then the input
/// is a dosDate value.

TTimeStamp::TTimeStamp(UInt_t tloc, Bool_t isUTC, Int_t secOffset, Bool_t dosDate)
{
   Set(tloc, isUTC, secOffset, dosDate);
}

////////////////////////////////////////////////////////////////////////////////
/// Return Greenwich mean sidereal time (GMST) in hour-angle. Return value
/// will always be between 0 and 24 (hours). Sidereal time is most accurately
/// calculated from UT1. If fSec and fNanoSec are in UTC (which they are by
/// default), the optional argument UT1Offset can be supplied (in
/// milliseconds). If UT1Offset is not supplied, conversion has maximum error
/// of 1s. If offset is supplied error can be reduced to us level. Values for
/// UT1Offset can be found in IERS Bulletin B:
/// ftp://ftp.iers.org/products/eop/bulletinb/format_2009/
/// The conversion to sidereal time used here is given by
/// Aoki et. al. Astron. Astrophys. 105, 359-362 (1982)
/// http://adsabs.harvard.edu/abs/1982A%26A...105..359A

Double_t TTimeStamp::AsGMST(Double_t UT1Offset) const
{
   Double_t D = (AsJulianDate() + UT1Offset/86400000.0) - 2451545.0;
   Double_t D_r = D - fmod(2.0*D+1.0, 2.0)/2.0;
   Double_t T = D_r/36525.0;
   Double_t sidereal = (24110.54841 + 8640184.812866*T + 0.093142*T*T
      - 0.0000062*T*T*T + (D - D_r)*86400.0*1.002737909350795)/3600.0;
   Double_t rval = fmod(sidereal, 24.0);
   return rval < 0 ? rval + 24.0 : rval;
}

////////////////////////////////////////////////////////////////////////////////
/// Return Greenwich apparent sidereal time (GAST) in hour-angle. Return
/// value will always be between 0 and 24 (hours). Sidereal time is most
/// accurately calculated from UT1. If fSec and fNanoSec are in UTC (which
/// they are by default), the optional argument UT1Offset can be supplied (in
/// milliseconds). If UT1Offset is not supplied, conversion has maximum error
/// of 1s. If offset is supplied error can be reduced to us level. Values for
/// UT1Offset can be found in IERS Bulletin B:
/// ftp://ftp.iers.org/products/eop/bulletinb/format_2009/
/// Equation of the equinoxes is given by USNO:
/// http://aa.usno.navy.mil/faq/docs/GAST.php

Double_t TTimeStamp::AsGAST(Double_t UT1Offset) const
{
   Double_t Pi = 3.14159265358979323846;
   Double_t D = (AsJulianDate() + UT1Offset/86400000.0) - 2451545.0;
   Double_t epsilon = (23.4393 - 0.0000004 * D) * Pi / 180.0;
   Double_t L = (280.47 + 0.98565 * D) * Pi / 180.0;
   Double_t Omega = (125.04 - 0.052954 * D) * Pi / 180.0;
   Double_t Deltapsi = -0.000319 * std::sin(Omega) - 0.000024 * std::sin(2.0 * L);
   Double_t eqeq = Deltapsi * std::cos(epsilon);
   Double_t rval = fmod(AsGMST(UT1Offset) + eqeq, 24.0);
   return rval < 0 ? rval + 24.0 : rval;
}

////////////////////////////////////////////////////////////////////////////////
/// Return local mean sidereal time (LMST) in hour-angle, given a longitude
/// in degrees. Return value will always be between 0 and 24 (hours).
/// Sidereal time is most accurately calculated from UT1. If fSec and
/// fNanoSec are in UTC (which they are by default), the optional argument
/// UT1Offset can be supplied (in milliseconds). If UT1Offset is not
/// supplied, conversion has maximum error of 1s. If offset is supplied error
/// can be reduced to us level. Values for UT1Offset can be found in IERS
/// Bulletin B: ftp://ftp.iers.org/products/eop/bulletinb/format_2009/

Double_t TTimeStamp::AsLMST(Double_t Longitude, Double_t UT1Offset) const
{
   Double_t rval = fmod(AsGMST(UT1Offset) + Longitude/15.0, 24.0);
   return rval < 0 ? rval + 24.0 : rval;
}

////////////////////////////////////////////////////////////////////////////////
/// Return local apparent sidereal time (LAST) in hour-angle, given a
/// longitude in degrees. Return value will always be between 0 and 24
/// (hours). Sidereal time is most accurately calculated from UT1. If fSec
/// and fNanoSec are in UTC (which they are by default), the optional
/// argument UT1Offset can be supplied (in milliseconds). If UT1Offset is not
/// supplied, conversion has maximum error of 1s. If offset is supplied error
/// can be reduced to us level. Values for UT1Offset can be found in IERS
/// Bulletin B: ftp://ftp.iers.org/products/eop/bulletinb/format_2009/

Double_t TTimeStamp::AsLAST(Double_t Longitude, Double_t UT1Offset) const
{
   Double_t rval = fmod(AsGAST(UT1Offset) + Longitude/15.0, 24.0);
   return rval < 0 ? rval + 24.0 : rval;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the date & time as a string.
///
/// Result is pointer to a statically allocated string.
/// User should copy this into their own buffer before calling
/// this method again.
///
/// Option "l" returns it in local zone format
/// (can be applied to default or compact format).
///
/// Default format is RFC822 compliant:
/// ~~~ {.cpp}
///   "Mon, 02 Jan 2001 18:11:12 +0000 (GMT) +999999999 nsec"
///   "Mon, 02 Jan 2001 10:11:12 -0800 (PST) +999999999 nsec"
/// ~~~
///
/// Option "c" compact is (almost) ISO 8601 compliant:
///  - "2001-01-02 18:11:12.9999999999Z"
///  - "2001-01-02 10:11:12.9999999999-0800"  if PST
///      * uses "-" as date separator as specified in ISO 8601
///      * uses "." rather than preferred "," for decimal separator
///      * -HHMM is the difference between local and UTC (if behind, + if ahead).
///
///   The "-HHMM" is replaced with "Z" if given as UTC.
///   To be strictly conforming it should use "T" instead of the
///   blank separating the date and time.
///
/// Option "2" returns as {sec,nsec} integers.
///
/// Option "s" returns "2001-01-02 18:11:12" with an implied UTC,
/// overrides "l" option.
///
/// Internally uses a circular list of buffers to avoid problems
/// using AsString multiple times in a single statement.

const Char_t *TTimeStamp::AsString(Option_t *option) const
{
   const Int_t nbuffers = 8;     // # of buffers

   static Char_t formatted[nbuffers][64];  // strftime fields substituted
   static Char_t formatted2[nbuffers][64]; // nanosec field substituted

   static Int_t ibuffer = nbuffers;

   R__LOCKGUARD2(gTimeMutex);

   ibuffer = (ibuffer+1)%nbuffers; // each call moves to next buffer

   TString opt = option;
   opt.ToLower();

   if (opt.Contains("2")) {
      // return string formatted as integer {sec,nsec}
      sprintf(formatted[ibuffer], "{%d,%d}", fSec, fNanoSec);
      return formatted[ibuffer];
   }

#ifdef R__LINUX
   // under linux %z is the hour offset and %Z is the timezone name
   const Char_t *kRFC822   = "%a, %d %b %Y %H:%M:%S %z (%Z) +#9ld nsec";
   const Char_t *kISO8601  = "%Y-%m-%d %H:%M:%S.#9.9ld%z";
   const Char_t *kISO8601Z = "%Y-%m-%d %H:%M:%S.#9.9ldZ";
#else
   // otherwise only %Z is guaranteed to be defined
   const Char_t *kRFC822   = "%a, %d %b %Y %H:%M:%S %Z +#9ld nsec";
   const Char_t *kISO8601  = "%Y-%m-%d %H:%M:%S.#9.9ld%Z";
   const Char_t *kISO8601Z = "%Y-%m-%d %H:%M:%S.#9.9ldZ";
#endif
   const Char_t *kSQL = "%Y-%m-%d %H:%M:%S";

   Bool_t asLocal = opt.Contains("l");
   Bool_t asSQL   = opt.Contains("s");
   if (asSQL) asLocal = kFALSE;

   const Char_t *format = kRFC822;
   if (opt.Contains("c")) {
      format = kISO8601;
      if (!asLocal) format = kISO8601Z;
   }
   if (asSQL) format = kSQL;

   // get the components into a tm struct
   time_t seconds = (time_t) fSec;
#ifdef _REENTRANT
   struct tm buf;
   struct tm *ptm = (asLocal) ? localtime_r(&seconds, &buf) : gmtime_r(&seconds, &buf);
#else
   struct tm *ptm = (asLocal) ? localtime(&seconds) : gmtime(&seconds);
#endif

   // format all but the nsec field
   // size_t length =
   strftime(formatted[ibuffer], sizeof(formatted[ibuffer]), format, ptm);

   if (asSQL) return formatted[ibuffer];

   // hack in the nsec part
   Char_t *ptr = strrchr(formatted[ibuffer], '#');
   if (ptr) *ptr = '%';    // substitute % for #
   sprintf(formatted2[ibuffer], formatted[ibuffer], fNanoSec);

   return formatted2[ibuffer];
}

////////////////////////////////////////////////////////////////////////////////
/// Copy this to ts.

void TTimeStamp::Copy(TTimeStamp &ts) const
{
   ts.fSec     = fSec;
   ts.fNanoSec = fNanoSec;
}

////////////////////////////////////////////////////////////////////////////////
/// Return date in form of 19971224 (i.e. 24/12/1997),
/// if non-zero pointers supplied for year, month, day fill those as well.

UInt_t TTimeStamp::GetDate(Bool_t inUTC, Int_t secOffset,
                           UInt_t *year, UInt_t *month, UInt_t *day) const
{
   time_t atime = fSec + secOffset;
#ifdef _REENTRANT
   struct tm buf;
   struct tm *ptm = (inUTC) ? gmtime_r(&atime, &buf) : localtime_r(&atime, &buf);
#else
   struct tm *ptm = (inUTC) ? gmtime(&atime) : localtime(&atime);
#endif

   if (day)   *day   = ptm->tm_mday;
   if (month) *month = ptm->tm_mon + 1;
   if (year)  *year  = ptm->tm_year + 1900;

   return (1900+ptm->tm_year)*10000 + (1+ptm->tm_mon)*100 + ptm->tm_mday;
}

////////////////////////////////////////////////////////////////////////////////
/// Return time in form of 123623 (i.e. 12:36:23),
/// if non-zero pointers supplied for hour, min, sec fill those as well.

UInt_t TTimeStamp::GetTime(Bool_t inUTC, Int_t secOffset,
                           UInt_t *hour, UInt_t *min, UInt_t *sec) const
{
   time_t atime = fSec + secOffset;
#ifdef _REENTRANT
   struct tm buf;
   struct tm *ptm = (inUTC) ? gmtime_r(&atime, &buf) : localtime_r(&atime, &buf);
#else
   struct tm *ptm = (inUTC) ? gmtime(&atime) : localtime(&atime);
#endif

   if (hour) *hour = ptm->tm_hour;
   if (min)  *min  = ptm->tm_min;
   if (sec)  *sec  = ptm->tm_sec;

   return ptm->tm_hour*10000 + ptm->tm_min*100 + ptm->tm_sec;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the day of the year represented by this time stamp value.
/// Valid return values range between 1 and 366, where January 1 = 1.

Int_t TTimeStamp::GetDayOfYear(Bool_t inUTC, Int_t secOffset) const
{
   time_t atime = fSec + secOffset;
#ifdef _REENTRANT
   struct tm buf;
   struct tm *ptm = (inUTC) ? gmtime_r(&atime, &buf) : localtime_r(&atime, &buf);
#else
   struct tm *ptm = (inUTC) ? gmtime(&atime) : localtime(&atime);
#endif

   Int_t day   = ptm->tm_mday;
   Int_t month = ptm->tm_mon + 1;
   Int_t year  = ptm->tm_year + 1900;

   return GetDayOfYear(day, month, year);
}

////////////////////////////////////////////////////////////////////////////////
/// Method is using Zeller's formula for calculating the day number.
/// Valid return values range between 1 and 7, where Monday = 1.

Int_t TTimeStamp::GetDayOfWeek(Bool_t inUTC, Int_t secOffset) const
{
   time_t atime = fSec + secOffset;
#ifdef _REENTRANT
   struct tm buf;
   struct tm *ptm = (inUTC) ? gmtime_r(&atime, &buf) : localtime_r(&atime, &buf);
#else
   struct tm *ptm = (inUTC) ? gmtime(&atime) : localtime(&atime);
#endif

   Int_t day   = ptm->tm_mday;
   Int_t month = ptm->tm_mon + 1;
   Int_t year  = ptm->tm_year + 1900;

   return GetDayOfWeek(day, month, year);
}

////////////////////////////////////////////////////////////////////////////////
/// Get the month of the year. Valid return values are between 1 and 12.

Int_t TTimeStamp::GetMonth(Bool_t inUTC, Int_t secOffset) const
{
   time_t atime = fSec + secOffset;
#ifdef _REENTRANT
   struct tm buf;
   struct tm *ptm = (inUTC) ? gmtime_r(&atime, &buf) : localtime_r(&atime, &buf);
#else
   struct tm *ptm = (inUTC) ? gmtime(&atime) : localtime(&atime);
#endif

   return ptm->tm_mon + 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the week of the year. Valid week values are between 1 and 53.
/// The return value is the year*100+week (1 Jan may be in the last
/// week of the previous year so the year must be returned too).

Int_t TTimeStamp::GetWeek(Bool_t inUTC, Int_t secOffset) const
{
   time_t atime = fSec + secOffset;
#ifdef _REENTRANT
   struct tm buf;
   struct tm *ptm = (inUTC) ? gmtime_r(&atime, &buf) : localtime_r(&atime, &buf);
#else
   struct tm *ptm = (inUTC) ? gmtime(&atime) : localtime(&atime);
#endif

   Int_t day   = ptm->tm_mday;
   Int_t month = ptm->tm_mon + 1;
   Int_t year  = ptm->tm_year + 1900;

   return GetWeek(day, month, year);
}

////////////////////////////////////////////////////////////////////////////////
/// Is the year a leap year.
/// The calendar year is 365 days long, unless the year is exactly divisible
/// by 4, in which case an extra day is added to February to make the year
/// 366 days long. If the year is the last year of a century, eg. 1700, 1800,
/// 1900, 2000, then it is only a leap year if it is exactly divisible by
/// 400. Therefore, 1900 wasn't a leap year but 2000 was. The reason for
/// these rules is to bring the average length of the calendar year into
/// line with the length of the Earth's orbit around the Sun, so that the
/// seasons always occur during the same months each year.

Bool_t TTimeStamp::IsLeapYear(Bool_t inUTC, Int_t secOffset) const
{
   time_t atime = fSec + secOffset;
#ifdef _REENTRANT
   struct tm buf;
   struct tm *ptm = (inUTC) ? gmtime_r(&atime, &buf) : localtime_r(&atime, &buf);
#else
   struct tm *ptm = (inUTC) ? gmtime(&atime) : localtime(&atime);
#endif

   Int_t year = ptm->tm_year + 1900;

   return IsLeapYear(year);
}

////////////////////////////////////////////////////////////////////////////////
/// Static method returning local (current) time zone offset from UTC.
/// This is the value in seconds one must add to the local time to arrive at
/// Coordinated Universal Time, so it is negative east of the Prime Meridian.

Int_t TTimeStamp::GetZoneOffset()
{
   // ?? should tzset (_tzset) be called?
#ifndef R__WIN32
   tzset();
#if defined(R__WINGCC)
   return _timezone;
#else
#if !defined(R__FBSD) && !defined(R__OBSD)
   return  timezone;   // unix has extern long int
#else
   time_t tp = 0;
   time(&tp);
#ifdef _REENTRANT
   struct tm buf;
   return -localtime_r(&tp, &buf)->tm_gmtoff;
#else
   return -localtime(&tp)->tm_gmtoff;
#endif
#endif
#endif
#else
   _tzset();
   return _timezone;   // Win32 prepends "_"
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Add "offset" as a delta time.

void TTimeStamp::Add(const TTimeStamp &offset)
{
   fSec     += offset.fSec;
   fNanoSec += offset.fNanoSec;
   NormalizeNanoSec();
}

////////////////////////////////////////////////////////////////////////////////
/// Print date and time.

void TTimeStamp::Print(Option_t *option) const
{
   printf("Date/Time = %s\n", AsString(option));
}

////////////////////////////////////////////////////////////////////////////////
/// Set Date/Time to current time as reported by the system.
/// No accounting for nanoseconds with std ANSI functions,
/// ns part faked so that subsequent calls simply add 1 to it
/// this ensures that calls within the same second come back
/// distinct (and sortable). Time is since Jan 1, 1970.

void TTimeStamp::Set()
{
#ifdef R__WIN32
   ULARGE_INTEGER time;
   GetSystemTimeAsFileTime((FILETIME *)&time);
   // NT keeps time in FILETIME format which is 100ns ticks since
   // Jan 1, 1601. TTimeStamp uses time in 100ns ticks since Jan 1, 1970,
   // the difference is 134774 days.
   fNanoSec = Int_t((time.QuadPart * (unsigned __int64) 100) %
                    (unsigned __int64) 1000000000);
   time.QuadPart -=
            (unsigned __int64) (1000*1000*10)       // seconds
          * (unsigned __int64) (60 * 60 * 24)       // days
          * (unsigned __int64) (134774);            // # of days

   fSec     = Int_t(time.QuadPart/(unsigned __int64) (1000*1000*10));
#else
   struct timeval tp;
   gettimeofday(&tp, 0);
   fSec     = tp.tv_sec;
   fNanoSec = tp.tv_usec*1000;
#endif

   static Int_t sec = 0, nsec = 0, fake_ns = 0;

   R__LOCKGUARD2(gTimeMutex);

   if (fSec == sec && fNanoSec == nsec)
      fNanoSec += ++fake_ns;
   else {
      fake_ns = 0;
      sec     = fSec;
      nsec    = fNanoSec;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set Date/Time from components.
///
/// Month & day both use normal 1..12 and 1..31 counting,
/// hours, min, sec run from 0 to 23, 59, 59 respectively,
/// secOffset provides method for adjusting for alternative timezones
///
/// ~~~ {.cpp}
/// "year"  |    0    1 ... 37 | 38...69   |   70 .. 100  101 ..  137
/// true    | 2000 2001   2037 | undefined | 1970   2000 2001 .. 2037
///
/// "year"  | 138...1969 | 1970 .. 2037 | ...
/// true    | undefined  | 1970 .. 2037 | undefined
/// ~~~

void TTimeStamp::Set(Int_t year, Int_t month, Int_t day,
                     Int_t hour, Int_t min, Int_t sec,
                     Int_t nsec, Bool_t isUTC, Int_t secOffset)
{
   // deal with special formats of year
   if (year <= 37)                year += 2000;
   if (year >= 70 && year <= 137) year += 1900;
   // tm.tm_year is years since 1900
   if (year >= 1900)              year -= 1900;

   struct tm tmstruct;
   tmstruct.tm_year  = year;    // years since 1900
   tmstruct.tm_mon   = month-1; // months since Jan [0,11]
   tmstruct.tm_mday  = day;     // day of the month [1,31]
   tmstruct.tm_hour  = hour;    // hours since midnight [0,23]
   tmstruct.tm_min   = min;     // minutes after the hour [0,59]
   tmstruct.tm_sec   = sec + secOffset;  // seconds after the minute [0,59]
   tmstruct.tm_isdst = -1;     // let "mktime" determine DST setting

   const time_t bad_time_t = (time_t) -1;
   // convert tm struct to time_t, if values are given in UTC then
   // no standard routine exists and we'll have to use our homegrown routine,
   // if values are given in local time then use "mktime"
   // which also normalizes the tm struct as a byproduct
   time_t utc_sec = (isUTC) ? MktimeFromUTC(&tmstruct) : mktime(&tmstruct);

   if (utc_sec == bad_time_t)
      Error("TTimeStamp::Set","mktime returned -1");

   fSec     = utc_sec;
   fNanoSec = nsec;

   NormalizeNanoSec();
}

////////////////////////////////////////////////////////////////////////////////
/// Set date/time from integers of the form [yy]YYMMDD and HHMMSS,
/// assume UTC (UTC) components:
///
/// ~~~ {.cpp}
///  MM: 01=January .. 12=December
///  DD: 01 .. 31
///
///  HH: 00=midnight .. 23
///  MM: 00 .. 59
///  SS: 00 .. 69
/// ~~~
///
///  - Date must be in format 980418 or 19980418
///                       1001127 or 20001127  (i.e. year 100 = 2000),
///  - time must be in format 224512 (second precision),
///  - date must be >= 700101.

void TTimeStamp::Set(Int_t date, Int_t time, Int_t nsec,
                     Bool_t isUTC, Int_t secOffset)
{
   Int_t year  = date/10000;
   Int_t month = (date-year*10000)/100;
   Int_t day   = date%100;

   // protect against odd attempts at time offsets
   const Int_t oneday = 240000;
   while (time < 0) {
      time += oneday;
      day  -= 1;
   }
   while (time > oneday) {
      time -= oneday;
      day  += 1;
   }
   Int_t hour  = time/10000;
   Int_t min   = (time-hour*10000)/100;
   Int_t sec   = time%100;

   Set(year, month, day, hour, min, sec, nsec, isUTC, secOffset);
}

////////////////////////////////////////////////////////////////////////////////
/// The input arg is a time_t value returned by time() or a value
/// returned by Convert(). This value is the number of seconds since
/// the EPOCH (i.e. 00:00:00 on Jan 1m 1970). If dosDate is true then
/// the input is a dosDate value.

void TTimeStamp::Set(UInt_t tloc, Bool_t isUTC, Int_t secOffset, Bool_t dosDate)
{
   struct tm localtm;
   memset (&localtm, 0, sizeof (localtm));

   if (dosDate) {
      localtm.tm_year  = ((tloc >> 25) & 0x7f) + 80;
      localtm.tm_mon   = ((tloc >> 21) & 0xf);
      localtm.tm_mday  = (tloc >> 16) & 0x1f;
      localtm.tm_hour  = (tloc >> 11) & 0x1f;
      localtm.tm_min   = (tloc >> 5) & 0x3f;
      localtm.tm_sec   = (tloc & 0x1f) * 2 + secOffset;
      localtm.tm_isdst = -1;
   } else {
      time_t t = (time_t) tloc;
#ifdef _REENTRANT
      struct tm tpa;
      struct tm *tp = localtime_r(&t, &tpa);
#else
      struct tm *tp = localtime(&t);
#endif
      localtm.tm_year  = tp->tm_year;
      localtm.tm_mon   = tp->tm_mon;
      localtm.tm_mday  = tp->tm_mday;
      localtm.tm_hour  = tp->tm_hour;
      localtm.tm_min   = tp->tm_min;
      localtm.tm_sec   = tp->tm_sec + secOffset;
      localtm.tm_isdst = -1;
   }

   const time_t bad_time_t = (time_t) -1;
   // convert tm struct to time_t, if values are given in UTC then
   // no standard routine exists and we'll have to use our homegrown routine,
   // if values are given in local time then use "mktime"
   // which also normalizes the tm struct as a byproduct
   time_t utc_sec = (isUTC && dosDate) ? MktimeFromUTC(&localtm) : mktime(&localtm);

   if (utc_sec == bad_time_t)
      Error("TTimeStamp::Set","mktime returned -1");

   fSec     = utc_sec;
   fNanoSec = 0;  //nsec;

   NormalizeNanoSec();
}

////////////////////////////////////////////////////////////////////////////////
/// Ensure that the fNanoSec field is in range [0,999999999].

void TTimeStamp::NormalizeNanoSec()
{
   const Int_t kNsPerSec = 1000000000;
   // deal with negative values
   while (fNanoSec < 0) {
      fNanoSec += kNsPerSec;
      fSec -= 1;
   }

   // deal with values inf fNanoSec greater than one sec
   while (fNanoSec >= kNsPerSec) {
      fNanoSec -= kNsPerSec;
      fSec += 1;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Equivalent of standard routine "mktime" but
/// using the assumption that tm struct is filled with UTC, not local, time.
///
/// This version *ISN'T* configured to handle every possible
/// weirdness of out-of-range values in the case of normalizing
/// the tm struct.
///
/// This version *DOESN'T* correctly handle values that can't be
/// fit into a time_t (i.e. beyond year 2038-01-18 19:14:07, or
/// before the start of Epoch).

time_t TTimeStamp::MktimeFromUTC(tm_t *tmstruct)
{
   Int_t daysInMonth[] = { 31, 0, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
   Int_t year = tmstruct->tm_year + 1900;
   daysInMonth[1] = IsLeapYear(year) ? 29 : 28;

   // fill in tmstruct->tm_yday

   Int_t &ref_tm_mon = tmstruct->tm_mon;
   Int_t &ref_tm_mday = tmstruct->tm_mday;
   // count days in months past
   tmstruct->tm_yday = 0;
   for (Int_t imonth = 0; imonth < ref_tm_mon; imonth++) {
      tmstruct->tm_yday += daysInMonth[imonth];
   }
   tmstruct->tm_yday += ref_tm_mday - 1;  // day [1-31] but yday [0-365]

   // adjust if day in this month is more than the month has
   while (ref_tm_mday > daysInMonth[ref_tm_mon]) {
      ref_tm_mday -= daysInMonth[ref_tm_mon];
      ref_tm_mon++;
   }

   // *should* calculate tm_wday (0-6) here ...

   // UTC is never DST
   tmstruct->tm_isdst = 0;

   // Calculate seconds since the Epoch based on formula in
   // POSIX  IEEEE Std 1003.1b-1993 pg 22

   Int_t utc_sec = tmstruct->tm_sec +
                   tmstruct->tm_min*60 +
                   tmstruct->tm_hour*3600 +
                   tmstruct->tm_yday*86400 +
                   (tmstruct->tm_year-70)*31536000 +
                   ((tmstruct->tm_year-69)/4)*86400;

   return utc_sec;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the day of the year represented by day, month and year.
/// Valid return values range between 1 and 366, where January 1 = 1.

Int_t TTimeStamp::GetDayOfYear(Int_t day, Int_t month, Int_t year)
{
   Int_t dayOfYear = 0;
   Int_t daysInMonth[] = { 31, 0, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
   daysInMonth[1] = IsLeapYear(year) ? 29 : 28;

   for (Int_t i = 0; i < (month - 1); i++)
      dayOfYear += daysInMonth[i];
   dayOfYear += day;

   return dayOfYear;
}

////////////////////////////////////////////////////////////////////////////////
/// Method is using Zeller's formula for calculating the day number.
/// Valid return values range between 1 and 7, where Monday = 1.

Int_t TTimeStamp::GetDayOfWeek(Int_t day, Int_t month, Int_t year)
{
   Int_t dayno;

   if (month < 3) {
      year--;
      month += 12;
   }

   dayno = 1 + day + 2*month + 3*(month + 1)/5 + year + year/4 - year/100 + year/400;
   dayno %= 7;

   // make monday first day of week
   return ((dayno == 0) ? 7 : dayno);
}

////////////////////////////////////////////////////////////////////////////////
/// Get the week of the year. Valid week values are between 1 and 53.
/// The return value is the year*100+week (1 Jan may be in the last
/// week of the previous year so the year must be returned too).

Int_t TTimeStamp::GetWeek(Int_t day, Int_t month, Int_t year)
{
   Int_t dayOfYear = GetDayOfYear(day, month, year);
   Int_t dayJan1st = GetDayOfWeek(1, 1, year);
   Int_t week = (dayOfYear + dayJan1st - 2) / 7 + 1;

   if (dayJan1st > 4)
      week--;

   if (week == 53) {
      Int_t dayNextJan1st = GetDayOfWeek(1, 1, year + 1);
      if (dayNextJan1st > 1 && dayNextJan1st < 5) {
         year++;
         week = 1;
      }
   } else if (week == 0) {
      Int_t dayPrevJan1st = GetDayOfWeek(1, 1, year - 1);
      week = (dayPrevJan1st < 5 && dayJan1st > 4) ? 53 : 52;
      year--;
   }
   return year * 100 + week;
}

////////////////////////////////////////////////////////////////////////////////
/// Is the given year a leap year.
/// The calendar year is 365 days long, unless the year is exactly divisible
/// by 4, in which case an extra day is added to February to make the year
/// 366 days long. If the year is the last year of a century, eg. 1700, 1800,
/// 1900, 2000, then it is only a leap year if it is exactly divisible by
/// 400. Therefore, 1900 wasn't a leap year but 2000 was. The reason for
/// these rules is to bring the average length of the calendar year into
/// line with the length of the Earth's orbit around the Sun, so that the
/// seasons always occur during the same months each year.

Bool_t TTimeStamp::IsLeapYear(Int_t year)
{
   return (year % 4 == 0) && !((year % 100 == 0) && (year % 400 > 0));
}

////////////////////////////////////////////////////////////////////////////////
/// Print out the "tm" structure:
/// ~~~ {.cpp}
/// tmstruct.tm_year = year;    // years since 1900
/// tmstruct.tm_mon  = month-1; // months since Jan [0,11]
/// tmstruct.tm_mday = day;     // day of the month [1,31]
/// tmstruct.tm_hour = hour;    // hours since midnight [0,23]
/// tmstruct.tm_min  = min;     // minutes after the hour [0,59]
/// tmstruct.tm_sec  = sec;     // seconds after the minute [0,59]
/// tmstruct.tm_wday            // day of week [0,6]
/// tmstruct.tm_yday            // days in year [0,365]
/// tmstruct.tm_isdst           // DST [-1/0/1]  (unknown,false,true)
/// ~~~

void TTimeStamp::DumpTMStruct(const tm_t &tmstruct)
{
   printf(" tm { year %4d, mon   %2d, day   %2d,\n",
          tmstruct.tm_year,
          tmstruct.tm_mon,
          tmstruct.tm_mday);
   printf("      hour   %2d, min   %2d, sec   %2d,\n",
          tmstruct.tm_hour,
          tmstruct.tm_min,
          tmstruct.tm_sec);
   printf("      wday   %2d, yday %3d, isdst %2d",
          tmstruct.tm_wday,
          tmstruct.tm_yday,
          tmstruct.tm_isdst);
#if (defined(linux) && !defined(R__WINGCC)) || defined(R__MACOSX)
   printf(",\n      tm_gmtoff %6ld, tm_zone \"%s\"",
#if defined(__USE_BSD) || defined(R__MACOSX) || defined(__USE_MISC)
   tmstruct.tm_gmtoff, tmstruct.tm_zone);
#else
   tmstruct.__tm_gmtoff, tmstruct.__tm_zone);
#endif
#endif
   printf(" }\n");
}
