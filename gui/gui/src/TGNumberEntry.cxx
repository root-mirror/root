// @(#)root/gui:$Id$
// Author: Daniel Sigg   03/09/2001

/*************************************************************************
 * Copyright (C) 1995-2001, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


/** \class TGNumberEntry
    \ingroup guiwidgets

TGNumberEntry is a number entry input widget with up/down buttons.
TGNumberEntryField is a number entry input widget.
TGNumberFormat contains enum types to specify the numeric format.

The number entry widget is based on TGTextEntry but allows only
numerical input. The widget support numerous formats including
integers, hex numbers, real numbers, fixed fraction reals and
time/date formats. The widget also allows to restrict input values
to non-negative or positive numbers and to specify explicit limits.

The following styles are supported:
  - kNESInteger:        integer number
  - kNESRealOne:        real number with one digit (no exponent)
  - kNESRealTwo:        real number with two digits (no exponent)
  - kNESRealThree:      real number with three digits (no exponent)
  - kNESRealFour:       real number with four digits (no exponent)
  - kNESReal:           arbitrary real number
  - kNESDegree:         angle in degree:minutes:seconds format
  - kNESMinSec:         time in minutes:seconds format
  - kNESHourMin:        time in hour:minutes format
  - kNESHourMinSec:     time in hour:minutes:seconds format
  - kNESDayMYear:       date in day/month/year format
  - kNESMDayYear:       date in month/day/year format
  - kNESHex:            hex number

The following attributes can be specified:
  - kNEAAnyNumber:      any number is allowed
  - kNEANonNegative:    only non-negative numbers are allowed
  - kNEAPositive:       only positive numbers are allowed

Explicit limits can be specified individually:
  - kNELNoLimits:       no limits
  - kNELLimitMin:       lower limit only
  - kNELLimitMax        upper limit only
  - kNELLimitMinMax     both lower and upper limits

TGNumberEntryField is a plain vanilla entry field, whereas
TGNumberEntry adds two small buttons to increase and decrease the
numerical value in the field. The number entry widgets also support
using the up and down cursor keys to change the numerical values.
The step size can be selected with control and shift keys:
  - --                  small step (1 unit/factor of 3)
  - shift               medium step (10 units/factor of 10)
  - control             large step (100 units/factor of 30)
  - shift-control       huge step (1000 units/factor of 100)

The steps are either linear or logarithmic. The default behaviour
is set when the entry field is created, but it can be changed by
pressing the alt key at the same time.

Changing the number in the widget will generate the event:
  - kC_TEXTENTRY, kTE_TEXTCHANGED, widget id, 0.
Hitting the enter key will generate:
  - kC_TEXTENTRY, kTE_ENTER, widget id, 0.
Hitting the tab key will generate:
  - kC_TEXTENTRY, kTE_TAB, widget id, 0.

*/


#include "TGNumberEntry.h"
#include "KeySymbols.h"
#include "TTimer.h"
#include "TSystem.h"
#include "TGToolTip.h"
#include "TMath.h"
#include "TVirtualX.h"
#include "strlcpy.h"
#include "snprintf.h"

#include <cctype>
#include <iostream>


ClassImp(TGNumberFormat);
ClassImp(TGNumberEntryField);
ClassImp(TGNumberEntryLayout);
ClassImp(TGNumberEntry);



//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Miscellaneous routines for handling numeric values <-> strings       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
enum ERealStyle {         // Style of real
   kRSInt = 0,            // Integer
   kRSFrac = 1,           // Fraction only
   kRSExpo = 2,           // Exponent only
   kRSFracExpo = 3        // Fraction and Exponent
};

////////////////////////////////////////////////////////////////////////////////

struct RealInfo_t {
   ERealStyle fStyle{kRSInt};  // Style of real
   Int_t fFracDigits{0};       // Number of fractional digits
   Int_t fFracBase{0};         // Base of fractional digits
   Int_t fIntNum{0};           // Integer number
   Int_t fFracNum{0};          // Fraction
   Int_t fExpoNum{0};          // Exponent
   Int_t fSign{0};             // Sign
};

////////////////////////////////////////////////////////////////////////////////

const Double_t kEpsilon = 1E-12;

////////////////////////////////////////////////////////////////////////////////

const Int_t kDays[13] =
    { 0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

////////////////////////////////////////////////////////////////////////////////

static Long_t Round(Double_t x)
{
   if (x > 0) {
      return (Long_t) (x + 0.5);
   } else if (x < 0) {
      return (Long_t) (x - 0.5);
   } else {
      return 0;
   }
}

////////////////////////////////////////////////////////////////////////////////

static Long_t Truncate(Double_t x)
{
   if (x > 0) {
      return (Long_t) (x + kEpsilon);
   } else if (x < 0) {
      return (Long_t) (x - kEpsilon);
   } else {
      return 0;
   }
}

////////////////////////////////////////////////////////////////////////////////

static Bool_t IsLeapYear(Int_t year)
{
   return ((year % 4 == 0) && ((year % 100 != 0) || (year % 400 == 0)));
}

////////////////////////////////////////////////////////////////////////////////

static Bool_t IsGoodChar(char c, TGNumberFormat::EStyle style,
                         TGNumberFormat::EAttribute attr)
{
   if (isdigit(c)) {
      return kTRUE;
   }
   if (isxdigit(c) && (style == TGNumberFormat::kNESHex)) {
      return kTRUE;
   }
   if ((c == '-') && (style == TGNumberFormat::kNESInteger) &&
       (attr == TGNumberFormat::kNEAAnyNumber)) {
      return kTRUE;
   }
   if ((c == '-') &&
       ((style == TGNumberFormat::kNESRealOne) ||
        (style == TGNumberFormat::kNESRealTwo) ||
        (style == TGNumberFormat::kNESRealThree) ||
        (style == TGNumberFormat::kNESRealFour) ||
        (style == TGNumberFormat::kNESReal) ||
        (style == TGNumberFormat::kNESDegree) ||
        (style == TGNumberFormat::kNESMinSec)) &&
       (attr == TGNumberFormat::kNEAAnyNumber)) {
      return kTRUE;
   }
   if ((c == '-') && (style == TGNumberFormat::kNESReal)) {
      return kTRUE;
   }
   if (((c == '.') || (c == ',')) &&
       ((style == TGNumberFormat::kNESRealOne) ||
        (style == TGNumberFormat::kNESRealTwo) ||
        (style == TGNumberFormat::kNESRealThree) ||
        (style == TGNumberFormat::kNESRealFour) ||
        (style == TGNumberFormat::kNESReal) ||
        (style == TGNumberFormat::kNESDegree) ||
        (style == TGNumberFormat::kNESMinSec) ||
        (style == TGNumberFormat::kNESHourMin) ||
        (style == TGNumberFormat::kNESHourMinSec) ||
        (style == TGNumberFormat::kNESDayMYear) ||
        (style == TGNumberFormat::kNESMDayYear))) {
      return kTRUE;
   }
   if ((c == ':') &&
       ((style == TGNumberFormat::kNESDegree) ||
        (style == TGNumberFormat::kNESMinSec) ||
        (style == TGNumberFormat::kNESHourMin) ||
        (style == TGNumberFormat::kNESHourMinSec) ||
        (style == TGNumberFormat::kNESDayMYear) ||
        (style == TGNumberFormat::kNESMDayYear))) {
      return kTRUE;
   }
   if ((c == '/') &&
       ((style == TGNumberFormat::kNESDayMYear) ||
        (style == TGNumberFormat::kNESMDayYear))) {
      return kTRUE;
   }
   if (((c == 'e') || (c == 'E')) && (style == TGNumberFormat::kNESReal)) {
      return kTRUE;
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////

static char *EliminateGarbage(char *text,
                              TGNumberFormat::EStyle style,
                              TGNumberFormat::EAttribute attr)
{
   if (text == 0) {
      return 0;
   }
   for (Int_t i = strlen(text) - 1; i >= 0; i--) {
      if (!IsGoodChar(text[i], style, attr)) {
         memmove(text + i, text + i + 1, strlen(text) - i);
      }
   }
   return text;
}

////////////////////////////////////////////////////////////////////////////////

static Long_t IntStr(const char *text)
{
   Long_t l = 0;
   Int_t sign = 1;
   for (UInt_t i = 0; i < strlen(text); i++) {
      if (text[i] == '-') {
         sign = -1;
      } else if ((isdigit(text[i])) && (l < kMaxLong)) {
         l = 10 * l + (text[i] - '0');
      }
   }
   return sign * l;
}

////////////////////////////////////////////////////////////////////////////////

static char *StrInt(char *text, Long_t i, Int_t digits)
{
   snprintf(text, 250, "%li", TMath::Abs(i));
   TString s = text;
   while (digits > s.Length()) {
      s = "0" + s;
   }
   if (i < 0) {
      s = "-" + s;
   }
   strlcpy(text, (const char *) s, 250);
   return text;
}

////////////////////////////////////////////////////////////////////////////////

static TString StringInt(Long_t i, Int_t digits)
{
   char text[256];
   StrInt(text, i, digits);
   return TString(text);
}

////////////////////////////////////////////////////////////////////////////////

static char *RealToStr(char *text, const RealInfo_t & ri)
{
   char *p = text;
   if (text == 0) {
      return 0;
   }
   strlcpy(p, "", 256);
   if (ri.fSign < 0) {
      strlcpy(p, "-", 256);
      p++;
   }
   StrInt(p, TMath::Abs(ri.fIntNum), 0);
   p += strlen(p);
   if ((ri.fStyle == kRSFrac) || (ri.fStyle == kRSFracExpo)) {
      strlcpy(p, ".", 256-strlen(p));
      p++;
      StrInt(p, TMath::Abs(ri.fFracNum), ri.fFracDigits);
      p += strlen(p);
   }
   if ((ri.fStyle == kRSExpo) || (ri.fStyle == kRSFracExpo)) {
      strlcpy(p, "e", 256-strlen(p));
      p++;
      StrInt(p, ri.fExpoNum, 0);
      p += strlen(p);
   }
   return text;
}

////////////////////////////////////////////////////////////////////////////////

static Double_t StrToReal(const char *text, RealInfo_t & ri)
{
   char *s;
   char *frac;
   char *expo;
   char *minus;
   char buf[256];

   if ((text == 0) || (!text[0])) {
      ri.fStyle = kRSInt;
      ri.fIntNum = 0;
      ri.fSign = 1;
      return 0.0;
   }
   strlcpy(buf, text, sizeof(buf));
   s = buf;
   frac = strchr(s, '.');
   if (frac == 0) {
      frac = strchr(s, ',');
   }
   expo = strchr(s, 'e');
   minus = strchr(s, '-');
   if (expo == 0) {
      expo = strchr(s, 'E');
   }
   if ((frac != 0) && (expo != 0) && (frac > expo)) {
      frac = 0;
   }
   if ((minus != 0) && ((expo == 0) || (minus < expo))) {
      ri.fSign = -1;
   } else {
      ri.fSign = 1;
   }
   if ((frac == 0) && (expo == 0)) {
      ri.fStyle = kRSInt;
   } else if (frac == 0) {
      ri.fStyle = kRSExpo;
   } else if (expo == 0) {
      ri.fStyle = kRSFrac;
   } else {
      ri.fStyle = kRSFracExpo;
   }
   if (frac != 0) {
      *frac = 0;
      frac++;
   }
   if (expo != 0) {
      *expo = 0;
      expo++;
   }
   ri.fIntNum = TMath::Abs(IntStr(s));
   if (expo != 0) {
      ri.fExpoNum = IntStr(expo);
   } else {
      ri.fExpoNum = 0;
   }
   if (ri.fExpoNum > 999) {
      ri.fExpoNum = 999;
   }
   if (ri.fExpoNum < -999) {
      ri.fExpoNum = -999;
   }
   ri.fFracDigits = 0;
   ri.fFracBase = 1;
   ri.fFracNum = 0;
   if (frac != 0) {
      for (UInt_t i = 0; i < strlen(frac); i++) {
         if (isdigit(frac[i])) {
            if (ri.fFracNum + 9 < kMaxInt / 10) {
               ri.fFracNum = 10 * ri.fFracNum + (frac[i] - '0');
               ri.fFracDigits++;
               ri.fFracBase *= 10;
            }
         }
      }
   }
   if ((ri.fFracDigits == 0) && (ri.fStyle == kRSFrac)) {
      ri.fStyle = kRSInt;
   }
   if ((ri.fFracDigits == 0) && (ri.fStyle == kRSFracExpo)) {
      ri.fStyle = kRSExpo;
   }
   switch (ri.fStyle) {
   case kRSInt:
      return ri.fSign * ri.fIntNum;
   case kRSFrac:
      return ri.fSign *
          (ri.fIntNum + (Double_t) ri.fFracNum / ri.fFracBase);
   case kRSExpo:
      return ri.fSign * (ri.fIntNum * TMath::Power(10, ri.fExpoNum));
   case kRSFracExpo:
      return ri.fSign * (ri.fIntNum +
                         (Double_t) ri.fFracNum / ri.fFracBase) *
          TMath::Power(10, ri.fExpoNum);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////

static ULong_t HexStrToInt(const char *s)
{
   ULong_t w = 0;
   for (UInt_t i = 0; i < strlen(s); i++) {
      if ((s[i] >= '0') && (s[i] <= '9')) {
         w = 16 * w + (s[i] - '0');
      } else if ((toupper(s[i]) >= 'A') && (toupper(s[i]) <= 'F')) {
         w = 16 * w + (toupper(s[i]) - 'A' + 10);
      }
   }
   return w;
}

////////////////////////////////////////////////////////////////////////////////

static char *IntToHexStr(char *text, ULong_t l)
{
   const char *const digits = "0123456789ABCDEF";
   char buf[64];
   char *p = buf + 62;
   // coverity[secure_coding]
   strcpy(p, "");
   while (l > 0) {
      *(--p) = digits[l % 16];
      l /= 16;
   }
   if (!p[0]) {
      // coverity[secure_coding]
      strcpy(text, "0");
   } else {
      // coverity[secure_coding]
      strcpy(text, p);
   }
   return text;
}

////////////////////////////////////////////////////////////////////////////////

static char *MIntToStr(char *text, Long_t l, Int_t digits)
{
   TString s;
   Int_t base;
   switch (digits) {
   case 0:
      base = 1;
      break;
   case 1:
      base = 10;
      break;
   case 2:
      base = 100;
      break;
   case 3:
      base = 1000;
      break;
   default:
   case 4:
      base = 10000;
      break;
   }
   s = StringInt(TMath::Abs(l) / base, 0) + "." +
       StringInt(TMath::Abs(l) % base, digits);
   if (l < 0) {
      s = "-" + s;
   }
   strlcpy(text, (const char *) s, 256);
   return text;
}

////////////////////////////////////////////////////////////////////////////////

static char *DIntToStr(char *text, Long_t l, Bool_t Sec, char Del)
{
   TString s;
   if (Sec) {
      s = StringInt(TMath::Abs(l) / 3600, 0) + Del +
          StringInt((TMath::Abs(l) % 3600) / 60, 2) + Del +
          StringInt(TMath::Abs(l) % 60, 2);
   } else {
      s = StringInt(TMath::Abs(l) / 60, 0) + Del +
          StringInt(TMath::Abs(l) % 60, 2);
   }
   if (l < 0) {
      s = "-" + s;
   }
   strlcpy(text, (const char *) s, 256);
   return text;
}

////////////////////////////////////////////////////////////////////////////////

static void GetNumbers(const char *s, Int_t & Sign,
                       Long_t & n1, Int_t maxd1,
                       Long_t & n2, Int_t maxd2,
                       Long_t & n3, Int_t maxd3, const char *Delimiters)
{
   Long_t n;
   Long_t d = 0;
   Sign = +1;
   n1 = 0;
   n2 = 0;
   n3 = 0;
   if (*s == '-') {
      Sign = -1;
      s++;
   }
   if (!isdigit(*s) && !strchr(Delimiters, *s)) {
      return;
   }
   while ((*s != 0) && ((strchr(Delimiters, *s) == 0) || (maxd2 == 0))) {
      if (isdigit(*s) && (d < maxd1)) {
         if (n1 < kMaxLong) {
            n1 = 10 * n1 + (*s - '0');
         }
         d++;
      }
      s++;
   }
   if (strcspn(s, Delimiters) == strlen(s)) {
      return;
   }
   Int_t dummy = 0;
   GetNumbers(s + 1, dummy, n2, maxd2, n3, maxd3, n, d, Delimiters);
}

////////////////////////////////////////////////////////////////////////////////

static Long_t GetSignificant(Long_t l, Int_t Max)
{
   while (TMath::Abs(l) >= Max) {
      l /= 10;
   }
   return l;
}

////////////////////////////////////////////////////////////////////////////////

static void AppendFracZero(char *text, Int_t digits)
{
   char *p;
   Int_t found = 0;
   p = strchr(text, '.');
   if (p == 0) {
      p = strchr(text, ',');
   }
   if (p == 0) {
      return;
   }
   p++;
   for (UInt_t i = 0; i < strlen(p); i++) {
      if (isdigit(*p)) {
         found++;
      }
   }
   while (found < digits) {
      // coverity[secure_coding]
      strcpy(p + strlen(p), "0");
      found++;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create a number entry with year/month/day information.

static Long_t MakeDateNumber(const char * /*text*/, Long_t Day,
                             Long_t Month, Long_t Year)
{
   Day = TMath::Abs(Day);
   Month = TMath::Abs(Month);
   Year = TMath::Abs(Year);
   if (Year < 100) {
      Year += 2000;
   }
   Month = GetSignificant(Month, 100);
   if (Month > 12)
      Month = 12;
   if (Month == 0)
      Month = 1;
   Day = GetSignificant(Day, 100);
   if (Day == 0)
      Day = 1;
   if (Day > kDays[Month])
      Day = kDays[Month];
   if ((Month == 2) && (Day > 28) && !IsLeapYear(Year))
      Day = 28;
   return 10000 * Year + 100 * Month + Day;
}

////////////////////////////////////////////////////////////////////////////////
/// Translate a string to a number value.

static Long_t TranslateToNum(const char *text,
                             TGNumberFormat::EStyle style, RealInfo_t & ri)
{
   Long_t n1;
   Long_t n2;
   Long_t n3;
   Int_t sign;
   switch (style) {
   case TGNumberFormat::kNESInteger:
      GetNumbers(text, sign, n1, 12, n2, 0, n3, 0, "");
      return sign * n1;
   case TGNumberFormat::kNESRealOne:
      GetNumbers(text, sign, n1, 12, n2, 1, n3, 0, ".,");
      return sign * (10 * n1 + GetSignificant(n2, 10));
   case TGNumberFormat::kNESRealTwo:
      {
         char buf[256];
         strlcpy(buf, text, sizeof(buf));
         AppendFracZero(buf, 2);
         GetNumbers(buf, sign, n1, 12, n2, 2, n3, 0, ".,");
         return sign * (100 * n1 + GetSignificant(n2, 100));
      }
   case TGNumberFormat::kNESRealThree:
      {
         char buf[256];
         strlcpy(buf, text, sizeof(buf));
         AppendFracZero(buf, 3);
         GetNumbers(buf, sign, n1, 12, n2, 3, n3, 0, ".,");
         return sign * (1000 * n1 + GetSignificant(n2, 1000));
      }
   case TGNumberFormat::kNESRealFour:
      {
         char buf[256];
         strlcpy(buf, text, sizeof(buf));
         AppendFracZero(buf, 4);
         GetNumbers(buf, sign, n1, 12, n2, 4, n3, 0, ".,");
         return sign * (10000 * n1 + GetSignificant(n2, 10000));
      }
   case TGNumberFormat::kNESReal:
      return (Long_t) StrToReal(text, ri);
   case TGNumberFormat::kNESDegree:
      GetNumbers(text, sign, n1, 12, n2, 2, n3, 2, ".,:");
      return sign * (3600 * n1 + 60 * GetSignificant(n2, 60) +
                     GetSignificant(n3, 60));
   case TGNumberFormat::kNESHourMinSec:
      GetNumbers(text, sign, n1, 12, n2, 2, n3, 2, ".,:");
      return 3600 * n1 + 60 * GetSignificant(n2, 60) +
          GetSignificant(n3, 60);
   case TGNumberFormat::kNESMinSec:
      GetNumbers(text, sign, n1, 12, n2, 2, n3, 0, ".,:");
      return sign * (60 * n1 + GetSignificant(n2, 60));
   case TGNumberFormat::kNESHourMin:
      GetNumbers(text, sign, n1, 12, n2, 2, n3, 0, ".,:");
      return 60 * n1 + GetSignificant(n2, 60);
   case TGNumberFormat::kNESDayMYear:
      GetNumbers(text, sign, n1, 2, n2, 2, n3, 4, ".,/");
      return MakeDateNumber(text, n1, n2, n3);
   case TGNumberFormat::kNESMDayYear:
      GetNumbers(text, sign, n2, 2, n1, 2, n3, 4, ".,/");
      return MakeDateNumber(text, n1, n2, n3);
   case TGNumberFormat::kNESHex:
      return HexStrToInt(text);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Translate a number value to a string.

static char *TranslateToStr(char *text, Long_t l,
                            TGNumberFormat::EStyle style, const RealInfo_t & ri)
{
   switch (style) {
   case TGNumberFormat::kNESInteger:
      return StrInt(text, l, 0);
   case TGNumberFormat::kNESRealOne:
      return MIntToStr(text, l, 1);
   case TGNumberFormat::kNESRealTwo:
      return MIntToStr(text, l, 2);
   case TGNumberFormat::kNESRealThree:
      return MIntToStr(text, l, 3);
   case TGNumberFormat::kNESRealFour:
      return MIntToStr(text, l, 4);
   case TGNumberFormat::kNESReal:
      return RealToStr(text, ri);
   case TGNumberFormat::kNESDegree:
      return DIntToStr(text, l, kTRUE, '.');
   case TGNumberFormat::kNESHourMinSec:
      return DIntToStr(text, l % (24 * 3600), kTRUE, ':');
   case TGNumberFormat::kNESMinSec:
      return DIntToStr(text, l, kFALSE, ':');
   case TGNumberFormat::kNESHourMin:
      return DIntToStr(text, l % (24 * 60), kFALSE, ':');
   case TGNumberFormat::kNESDayMYear:
      {
         TString date =
             StringInt(TMath::Abs(l) % 100, 0) + "/" +
             StringInt((TMath::Abs(l) / 100) % 100, 0) + "/" +
             StringInt(TMath::Abs(l) / 10000, 0);
         strlcpy(text, (const char *) date, 256);
         return text;
      }
   case TGNumberFormat::kNESMDayYear:
      {
         TString date =
             StringInt((TMath::Abs(l) / 100) % 100, 0) + "/" +
             StringInt(TMath::Abs(l) % 100, 0) + "/" +
             StringInt(TMath::Abs(l) / 10000, 0);
         strlcpy(text, (const char *) date, 256);
         return text;
      }
   case TGNumberFormat::kNESHex:
      return IntToHexStr(text, (ULong_t) l);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Convert to double format.

static Double_t RealToDouble(const RealInfo_t ri)
{
   switch (ri.fStyle) {
      // Integer type real
   case kRSInt:
      return (Double_t) ri.fSign * ri.fIntNum;
      // Fraction type real
   case kRSFrac:
      return (Double_t) ri.fSign * ((Double_t) TMath::Abs(ri.fIntNum) +
                                    (Double_t) ri.fFracNum / ri.fFracBase);
      // Exponent only
   case kRSExpo:
      return (Double_t) ri.fSign * ri.fIntNum *
          TMath::Power(10, ri.fExpoNum);
      // Fraction and exponent
   case kRSFracExpo:
      return (Double_t) ri.fSign * ((Double_t) TMath::Abs(ri.fIntNum) +
                                    (Double_t) ri.fFracNum /
                                    ri.fFracBase) * TMath::Power(10,
                                                                 ri.fExpoNum);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Check min/max limits for the set value.

static void CheckMinMax(Long_t & l, TGNumberFormat::EStyle style,
                        TGNumberFormat::ELimit limits,
                        Double_t min, Double_t max)
{
   if ((limits == TGNumberFormat::kNELNoLimits) ||
       (style == TGNumberFormat::kNESReal)) {
      return;
   }
   // check min
   if ((limits == TGNumberFormat::kNELLimitMin) ||
       (limits == TGNumberFormat::kNELLimitMinMax)) {
      Long_t lower;
      switch (style) {
      case TGNumberFormat::kNESRealOne:
         lower = Round(10.0 * min);
         break;
      case TGNumberFormat::kNESRealTwo:
         lower = Round(100.0 * min);
         break;
      case TGNumberFormat::kNESRealThree:
         lower = Round(1000.0 * min);
         break;
      case TGNumberFormat::kNESRealFour:
         lower = Round(10000.0 * min);
         break;
      case TGNumberFormat::kNESHex:
         lower = (ULong_t) Round(min);
         break;
      default:
         lower = Round(min);
         break;
      }
      if (style != TGNumberFormat::kNESHex) {
         if (l < lower)
            l = lower;
      } else {
         if (lower < 0)
            lower = 0;
         if ((ULong_t) l < (ULong_t) lower)
            l = lower;
      }
   }
   // check max
   if ((limits == TGNumberFormat::kNELLimitMax) ||
       (limits == TGNumberFormat::kNELLimitMinMax)) {
      Long_t upper;
      switch (style) {
      case TGNumberFormat::kNESRealOne:
         upper = Round(10.0 * max);
         break;
      case TGNumberFormat::kNESRealTwo:
         upper = Round(100.0 * max);
         break;
      case TGNumberFormat::kNESRealThree:
         upper = Round(1000.0 * max);
         break;
      case TGNumberFormat::kNESRealFour:
         upper = Round(10000.0 * max);
         break;
      case TGNumberFormat::kNESHex:
         upper = (ULong_t) Round(max);
         break;
      default:
         upper = Round(max);
         break;
      }
      if (style != TGNumberFormat::kNESHex) {
         if (l > upper)
            l = upper;
      } else {
         if (upper < 0)
            upper = 0;
         if ((ULong_t) l > (ULong_t) upper)
            l = upper;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Convert to double format.

static void IncreaseReal(RealInfo_t & ri, Double_t mag, Bool_t logstep,
                         TGNumberFormat::ELimit limits =
                         TGNumberFormat::kNELNoLimits, Double_t min = 0,
                         Double_t max = 1)
{
   Double_t x = RealToDouble(ri);

   // apply step
   if (logstep) {
      x *= mag;
   } else {
      switch (ri.fStyle) {
      case kRSInt:
         x = x + mag;
         break;
      case kRSFrac:
         x = x + mag / ri.fFracBase;
         break;
      case kRSExpo:
         x = x + mag * TMath::Power(10, ri.fExpoNum);
         break;
      case kRSFracExpo:
         x = x + (mag / ri.fFracBase) * TMath::Power(10, ri.fExpoNum);
         break;
      }
   }
   // check min
   if ((limits == TGNumberFormat::kNELLimitMin) ||
       (limits == TGNumberFormat::kNELLimitMinMax)) {
      if (x < min)
         x = min;
   }
   // check max
   if ((limits == TGNumberFormat::kNELLimitMax) ||
       (limits == TGNumberFormat::kNELLimitMinMax)) {
      if (x > max)
         x = max;
   }
   // check format after log step
   if ((x != 0) && logstep && (TMath::Abs(mag) > kEpsilon)) {
      for (int j = 0; j < 10; j++) {
         // Integer: special case
         if ((ri.fStyle == kRSInt) && (TMath::Abs(x) < 1) &&
             (TMath::Abs(x) > kEpsilon)) {
            ri.fStyle = kRSFrac;
            ri.fFracDigits = 1;
            ri.fFracBase = 10;
            continue;
         }
         if ((ri.fStyle == kRSInt) && (TMath::Abs(x) > 10000)) {
            ri.fStyle = kRSFracExpo;
            ri.fExpoNum = 4;
            ri.fFracDigits = 4;
            ri.fFracBase = 10000;
            Long_t rest = Round(TMath::Abs(x)) % 10000;
            for (int k = 0; k < 4; k++) {
               if (rest % 10 != 0) {
                  break;
               }
               ri.fFracDigits--;
               ri.fFracBase /= 10;
               rest /= 10;
            }
            if (ri.fFracDigits == 0) {
               ri.fStyle = kRSExpo;
            }
            continue;
         }
         if (ri.fStyle == kRSInt)
            break;

         // calculate first digit
         Double_t y;
         if ((ri.fStyle == kRSExpo) || (ri.fStyle == kRSFracExpo)) {
            y = TMath::Abs(x) * TMath::Power(10, -ri.fExpoNum);
         } else {
            y = TMath::Abs(x);
         }
         // adjust exponent if num < 1
         if ((Truncate(y) == 0) && (y > 0.001)) {
            if ((ri.fStyle == kRSExpo) || (ri.fStyle == kRSFracExpo)) {
               ri.fExpoNum--;
            } else {
               ri.fStyle = kRSFracExpo;
               ri.fExpoNum = -1;
            }
            continue;
         }
         // adjust exponent if num > 10
         if (Truncate(y) >= 10) {
            if ((ri.fStyle == kRSExpo) || (ri.fStyle == kRSFracExpo)) {
               ri.fExpoNum++;
            } else {
               ri.fStyle = kRSFracExpo;
               ri.fExpoNum = 1;
            }
            continue;
         }
         break;
      }
   }
   // convert back to RealInfo_t
   switch (ri.fStyle) {
      // Integer type real
   case kRSInt:
      {
         ri.fSign = (x < 0) ? -1 : 1;
         ri.fIntNum = Round(TMath::Abs(x));
         break;
      }
      // Fraction type real
   case kRSFrac:
      {
         ri.fSign = (x < 0) ? -1 : 1;
         ri.fIntNum = Truncate(TMath::Abs(x));
         ri.fFracNum = Round((TMath::Abs(x) - TMath::Abs(ri.fIntNum)) * ri.fFracBase);
         break;
      }
      // Exponent only
   case kRSExpo:
      {
         ri.fSign = (x < 0) ? -1 : 1;
         ri.fIntNum = Round(TMath::Abs(x) * TMath::Power(10, -ri.fExpoNum));
         if (ri.fIntNum == 0) {
            ri.fStyle = kRSInt;
         }
         break;
      }
      // Fraction and exponent
   case kRSFracExpo:
      {
         ri.fSign = (x < 0) ? -1 : 1;
         Double_t y = TMath::Abs(x) * TMath::Power(10, -ri.fExpoNum);
         ri.fIntNum = Truncate(y);
         ri.fFracNum = Round((y - TMath::Abs(ri.fIntNum)) * ri.fFracBase);
         if ((ri.fIntNum == 0) && (ri.fFracNum == 0)) {
            ri.fStyle = kRSFrac;
         }
         break;
      }
   }

   // check if the back conversion violated limits
   if (limits != TGNumberFormat::kNELNoLimits) {
      x = RealToDouble(ri);
      // check min
      if ((limits == TGNumberFormat::kNELLimitMin) ||
          (limits == TGNumberFormat::kNELLimitMinMax)) {
         if (x < min) {
            char text[256];
            snprintf(text, 255, "%g", min);
            StrToReal(text, ri);
         }
      }
      // check max
      if ((limits == TGNumberFormat::kNELLimitMax) ||
          (limits == TGNumberFormat::kNELLimitMinMax)) {
         if (x > max) {
            char text[256];
            snprintf(text, 255, "%g", max);
            StrToReal(text, ri);
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Change year/month/day format.

static void IncreaseDate(Long_t & l, TGNumberFormat::EStepSize step, Int_t sign)
{
   Long_t year;
   Long_t month;
   Long_t day;

   // get year/month/day format
   year = l / 10000;
   month = (TMath::Abs(l) / 100) % 100;
   if (month > 12)
      month = 12;
   if (month == 0)
      month = 1;
   day = TMath::Abs(l) % 100;
   if (day > kDays[month])
      day = kDays[month];
   if ((month == 2) && (day > 28) && !IsLeapYear(year)) {
      day = 28;
   }
   if (day == 0)
      day = 0;

   // apply step
   if (step == TGNumberFormat::kNSSHuge) {
      year += sign * 10;
   } else if (step == TGNumberFormat::kNSSLarge) {
      year += sign;
   } else if (step == TGNumberFormat::kNSSMedium) {
      month += sign;
      if (month > 12) {
         month = 1;
         year++;
      }
      if (month < 1) {
         month = 12;
         year--;
      }
   } else if (step == TGNumberFormat::kNSSSmall) {
      day += sign;
      if ((sign > 0) &&
          ((day > kDays[month]) ||
           ((month == 2) && (day > 28) && !IsLeapYear(year)))) {
         day = 1;
         month++;
         if (month > 12) {
            month = 1;
            year++;
         }
      }
      if ((sign < 0) && (day == 0)) {
         month--;
         if (month < 1) {
            month = 12;
            year--;
         }
         day = kDays[month];
      }
   }
   // check again for valid date
   if (year < 0)
      year = 0;
   if (day > kDays[month])
      day = kDays[month];
   if ((month == 2) && (day > 28) && !IsLeapYear(year)) {
      day = 28;
   }
   l = 10000 * year + 100 * month + day;
}




////////////////////////////////////////////////////////////////////////////////
/// Constructs a number entry field.

TGNumberEntryField::TGNumberEntryField(const TGWindow * p, Int_t id,
                                       Double_t val, GContext_t norm,
                                       FontStruct_t font, UInt_t option,
                                       ULong_t back)
   : TGTextEntry(p, new TGTextBuffer(), id, norm, font, option, back),
     fNeedsVerification(kFALSE), fNumStyle(kNESReal), fNumAttr(kNEAAnyNumber),
     fNumLimits(kNELNoLimits), fNumMin(0.0), fNumMax(1.0)
{
   fStepLog = kFALSE;
   SetAlignment(kTextRight);
   SetNumber(val);
   fEditDisabled = kEditDisable | kEditDisableGrab;
}

////////////////////////////////////////////////////////////////////////////////
/// Constructs a number entry field.

TGNumberEntryField::TGNumberEntryField(const TGWindow * parent,
                                       Int_t id, Double_t val,
                                       EStyle style, EAttribute attr,
                                       ELimit limits, Double_t min,
                                       Double_t max)
   : TGTextEntry(parent, "", id), fNeedsVerification(kFALSE), fNumStyle(style),
     fNumAttr(attr), fNumLimits(limits), fNumMin(min), fNumMax(max)
{
   fStepLog = kFALSE;
   SetAlignment(kTextRight);
   SetNumber(val);
   fEditDisabled = kEditDisable | kEditDisableGrab;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the numeric value (floating point representation).

void TGNumberEntryField::SetNumber(Double_t val, Bool_t emit)
{
   switch (fNumStyle) {
   case kNESInteger:
      SetIntNumber(Round(val), emit);
      break;
   case kNESRealOne:
      SetIntNumber(Round(10.0 * val), emit);
      break;
   case kNESRealTwo:
      SetIntNumber(Round(100.0 * val), emit);
      break;
   case kNESRealThree:
      SetIntNumber(Round(1000.0 * val), emit);
      break;
   case kNESRealFour:
      SetIntNumber(Round(10000.0 * val), emit);
      break;
   case kNESReal:
      {
         char text[256];
         snprintf(text, 255, "%g", val);
         SetText(text, emit);
         break;
      }
   case kNESDegree:
      SetIntNumber(Round(val), emit);
      break;
   case kNESHourMinSec:
      SetIntNumber(Round(val), emit);
      break;
   case kNESMinSec:
      SetIntNumber(Round(val), emit);
      break;
   case kNESHourMin:
      SetIntNumber(Round(val), emit);
      break;
   case kNESDayMYear:
      SetIntNumber(Round(val), emit);
      break;
   case kNESMDayYear:
      SetIntNumber(Round(val), emit);
      break;
   case kNESHex:
      SetIntNumber((UInt_t) (TMath::Abs(val) + 0.5), emit);
      break;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the numeric value (integer representation).

void TGNumberEntryField::SetIntNumber(Long_t val, Bool_t emit)
{
   char text[256];
   RealInfo_t ri;
   if (fNumStyle == kNESReal) {
      TranslateToStr(text, val, kNESInteger, ri);
   } else {
      TranslateToStr(text, val, fNumStyle, ri);
   }
   SetText(text, emit);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the numeric value (time format).

void TGNumberEntryField::SetTime(Int_t hour, Int_t min, Int_t sec, Bool_t emit)
{
   switch (fNumStyle) {
   case kNESHourMinSec:
      SetIntNumber(3600 * TMath::Abs(hour) + 60 * TMath::Abs(min) +
                   TMath::Abs(sec), emit);
      break;
   case kNESMinSec:
      {
         SetIntNumber(60 * min + sec, emit);
         break;
      }
   case kNESHourMin:
      SetIntNumber(60 * TMath::Abs(hour) + TMath::Abs(min), emit);
      break;
   default:
      break;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the numeric value (date format).

void TGNumberEntryField::SetDate(Int_t year, Int_t month, Int_t day, Bool_t emit)
{
   switch (fNumStyle) {
   case kNESDayMYear:
   case kNESMDayYear:
      {
         SetIntNumber(10000 * TMath::Abs(year) + 100 * TMath::Abs(month) +
                      TMath::Abs(day), emit);
      }
   default:
      {
         break;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the numeric value (hex format).

void TGNumberEntryField::SetHexNumber(ULong_t val, Bool_t emit)
{
   SetIntNumber((Long_t) val, emit);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the value (text format).

void TGNumberEntryField::SetText(const char *text, Bool_t emit)
{
   char buf[256];
   strlcpy(buf, text, sizeof(buf));
   EliminateGarbage(buf, fNumStyle, fNumAttr);
   TGTextEntry::SetText(buf, emit);
   fNeedsVerification = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the numeric value (floating point representation).

Double_t TGNumberEntryField::GetNumber() const
{
   switch (fNumStyle) {
   case kNESInteger:
      return (Double_t) GetIntNumber();
   case kNESRealOne:
      return (Double_t) GetIntNumber() / 10.0;
   case kNESRealTwo:
      return (Double_t) GetIntNumber() / 100.0;
   case kNESRealThree:
      return (Double_t) GetIntNumber() / 1000.0;
   case kNESRealFour:
      return (Double_t) GetIntNumber() / 10000.0;
   case kNESReal:
      {
         char text[256];
         RealInfo_t ri;
         strlcpy(text, GetText(), sizeof(text));
         return StrToReal(text, ri);
      }
   case kNESDegree:
      return (Double_t) GetIntNumber();
   case kNESHourMinSec:
      return (Double_t) GetIntNumber();
   case kNESMinSec:
      return (Double_t) GetIntNumber();
   case kNESHourMin:
      return (Double_t) GetIntNumber();
   case kNESDayMYear:
      return (Double_t) GetIntNumber();
   case kNESMDayYear:
      return (Double_t) GetIntNumber();
   case kNESHex:
      return (Double_t) (ULong_t) GetIntNumber();
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the numeric value (integer representation).

Long_t TGNumberEntryField::GetIntNumber() const
{
   RealInfo_t ri;
   return TranslateToNum(GetText(), fNumStyle, ri);
}

////////////////////////////////////////////////////////////////////////////////
/// Get the numeric value (time format).

void TGNumberEntryField::GetTime(Int_t & hour, Int_t & min, Int_t & sec) const
{
   switch (fNumStyle) {
   case kNESHourMinSec:
      {
         Long_t l = GetIntNumber();
         hour = TMath::Abs(l) / 3600;
         min = (TMath::Abs(l) % 3600) / 60;
         sec = TMath::Abs(l) % 60;
         break;
      }
   case kNESMinSec:
      {
         Long_t l = GetIntNumber();
         hour = 0;
         min = TMath::Abs(l) / 60;
         sec = TMath::Abs(l) % 60;
         if (l < 0) {
            min *= -1;
            sec *= -1;
         }
         break;
      }
   case kNESHourMin:
      {
         Long_t l = GetIntNumber();
         hour = TMath::Abs(l) / 60;
         min = TMath::Abs(l) % 60;
         sec = 0;
         break;
      }
   default:
      {
         hour = 0;
         min = 0;
         sec = 0;
         break;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get the numeric value (date format).

void TGNumberEntryField::GetDate(Int_t & year, Int_t & month, Int_t & day) const
{
   switch (fNumStyle) {
   case kNESDayMYear:
   case kNESMDayYear:
      {
         Long_t l = GetIntNumber();
         year = l / 10000;
         month = (l % 10000) / 100;
         day = l % 100;
         break;
      }
   default:
      {
         year = 0;
         month = 0;
         day = 0;
         break;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get the numeric value (hex format).

ULong_t TGNumberEntryField::GetHexNumber() const
{
   return (ULong_t) GetIntNumber();
}

////////////////////////////////////////////////////////////////////////////////
/// Get the text width in pixels.

Int_t TGNumberEntryField::GetCharWidth(const char *text) const
{
   return gVirtualX->TextWidth(fFontStruct, text, strlen(text));
}

////////////////////////////////////////////////////////////////////////////////
/// Increase the number value.

void TGNumberEntryField::IncreaseNumber(EStepSize step,
                                        Int_t stepsign, Bool_t logstep)
{
   Long_t l = 0;
   RealInfo_t ri;
   Long_t mag = 0;
   Double_t rmag = 0.0;
   Int_t sign = stepsign;

   // save old text field
   TString oldtext = GetText();
   // Get number
   if (fNumStyle != kNESReal) {
      l = GetIntNumber();
   } else {
      StrToReal(oldtext, ri);
   }

   // magnitude of step
   if ((fNumStyle == kNESDegree) || (fNumStyle == kNESHourMinSec) ||
       (fNumStyle == kNESMinSec) || (fNumStyle == kNESHourMin) ||
       (fNumStyle == kNESDayMYear) || (fNumStyle == kNESMDayYear) ||
       (fNumStyle == kNESHex)) {
      logstep = kFALSE;
      switch (step) {
      case kNSSSmall:
         mag = 1;
         break;
      case kNSSMedium:
         mag = 10;
         break;
      case kNSSLarge:
         mag = 100;
         break;
      case kNSSHuge:
         mag = 1000;
         break;
      }
   } else {
      Int_t msd = TMath::Abs((fNumStyle == kNESReal) ? ri.fIntNum : l);
      while (msd >= 10)
         msd /= 10;
      Bool_t odd = (msd < 3);
      if (sign < 0)
         odd = !odd;
      switch (step) {
      case kNSSSmall:
         rmag = (!logstep) ? 1. : (odd ? 3. : 10. / 3.);
         break;
      case kNSSMedium:
         rmag = 10.;
         break;
      case kNSSLarge:
         rmag = (!logstep) ? 100. : (odd ? 30. : 100. / 3.);
         break;
      case kNSSHuge:
         rmag = (!logstep) ? 1000. : 100.;
         break;
      }
      if (sign < 0)
         rmag = logstep ? 1. / rmag : -rmag;
   }

   // sign of step
   if (sign == 0) {
      logstep = kFALSE;
      rmag = 0;
      mag = 0;
   } else {
      sign = (sign > 0) ? 1 : -1;
   }
   // add/multiply step
   switch (fNumStyle) {
   case kNESInteger:
   case kNESRealOne:
   case kNESRealTwo:
   case kNESRealThree:
   case kNESRealFour:
      {
         l = logstep ? Round(l * rmag) : Round(l + rmag);
         CheckMinMax(l, fNumStyle, fNumLimits, fNumMin, fNumMax);
         if ((l < 0) && (fNumAttr == kNEANonNegative))
            l = 0;
         if ((l <= 0) && (fNumAttr == kNEAPositive))
            l = 1;
         break;
      }
   case kNESReal:
      {
         IncreaseReal(ri, rmag, logstep, fNumLimits, fNumMin, fNumMax);
         if (((fNumAttr == kNEANonNegative) ||
              (fNumAttr == kNEAPositive)) && (ri.fSign < 0)) {
            ri.fIntNum = 0;
            ri.fFracNum = 0;
            ri.fExpoNum = 0;
            ri.fSign = 1;
         }
         break;
      }
   case kNESDegree:
      {
         if (mag > 60)
            l += sign * 36 * mag;
         else if (mag > 6)
            l += sign * 6 * mag;
         else
            l += sign * mag;
         CheckMinMax(l, fNumStyle, fNumLimits, fNumMin, fNumMax);
         if ((l < 0) && (fNumAttr == kNEANonNegative))
            l = 0;
         if ((l <= 0) && (fNumAttr == kNEAPositive))
            l = 1;
         break;
      }
   case kNESHourMinSec:
      {
         if (mag > 60)
            l += sign * 36 * mag;
         else if (mag > 6)
            l += sign * 6 * mag;
         else
            l += sign * mag;
         CheckMinMax(l, fNumStyle, fNumLimits, fNumMin, fNumMax);
         if (l < 0)
            l = (24 * 3600) - ((-l) % (24 * 3600));
         if (l > 0)
            l = l % (24 * 3600);
         break;
      }
   case kNESMinSec:
      {
         if (mag > 6)
            l += sign * 6 * mag;
         else
            l += sign * mag;
         CheckMinMax(l, fNumStyle, fNumLimits, fNumMin, fNumMax);
         if ((l < 0) && (fNumAttr == kNEANonNegative))
            l = 0;
         if ((l <= 0) && (fNumAttr == kNEAPositive))
            l = 1;
         break;
      }
   case kNESHourMin:
      {
         if (mag > 6)
            l += sign * 6 * mag;
         else
            l += sign * mag;
         CheckMinMax(l, fNumStyle, fNumLimits, fNumMin, fNumMax);
         if (l < 0)
            l = (24 * 60) - ((-l) % (24 * 60));
         if (l > 0)
            l = l % (24 * 60);
         break;
      }
   case kNESDayMYear:
   case kNESMDayYear:
      {
         IncreaseDate(l, step, sign);
         CheckMinMax(l, fNumStyle, fNumLimits, fNumMin, fNumMax);
         break;
      }
   case kNESHex:
      {
         ULong_t ll = (ULong_t) l;
         if (mag > 500)
            ll += sign * 4096 * mag / 1000;
         else if (mag > 50)
            ll += sign * 256 * mag / 100;
         else if (mag > 5)
            ll += sign * 16 * mag / 10;
         else
            ll += sign * mag;
         l = (Long_t) ll;
         CheckMinMax(l, fNumStyle, fNumLimits, fNumMin, fNumMax);
         break;
      }
   }
   if (fNumStyle != kNESReal) {
      SetIntNumber(l);
   } else {
      char buf[256];
      RealToStr(buf, ri);
      SetText(buf);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the numerical format.

void TGNumberEntryField::SetFormat(EStyle style, EAttribute attr)
{
   Double_t val = GetNumber();
   fNumStyle = style;
   fNumAttr = attr;
   if ((fNumAttr != kNEAAnyNumber) && (val < 0))
      val = 0;
   SetNumber(val);
   // make sure we have a valid number by increasing it by 0
   IncreaseNumber(kNSSSmall, 0, kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the numerical limits.

void TGNumberEntryField::SetLimits(ELimit limits,
                                   Double_t min, Double_t max)
{
   Double_t val = GetNumber();
   fNumLimits = limits;
   fNumMin = min;
   fNumMax = max;
   SetNumber(val);
   // make sure we have a valid number by increasing it by 0
   IncreaseNumber(kNSSSmall, 0, kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the active state.

void TGNumberEntryField::SetState(Bool_t state)
{
   if (!state && fNeedsVerification) {
      // make sure we have a valid number by increasing it by 0
      IncreaseNumber(kNSSSmall, 0, kFALSE);
   }
   TGTextEntry::SetState(state);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle keys.

Bool_t TGNumberEntryField::HandleKey(Event_t * event)
{
   if (!IsEnabled()) {
      return TGTextEntry::HandleKey(event);
   }

   Int_t n;
   char tmp[10];
   UInt_t keysym;
   gVirtualX->LookupString(event, tmp, sizeof(tmp), keysym);
   n = strlen(tmp);

   // intercept up key
   if ((EKeySym) keysym == kKey_Up) {
      // Get log step / alt key
      Bool_t logstep = fStepLog;
      if (event->fState & kKeyMod1Mask)
         logstep = !logstep;
      // shift-cntrl-up
      if ((event->fState & kKeyShiftMask) &&
          (event->fState & kKeyControlMask)) {
         IncreaseNumber(kNSSHuge, 1, logstep);
      }
      // cntrl-up
      else if (event->fState & kKeyControlMask) {
         IncreaseNumber(kNSSLarge, 1, logstep);

      }
      // shift-up
      else if (event->fState & kKeyShiftMask) {
         IncreaseNumber(kNSSMedium, 1, logstep);
      }

      // up
      else {
         IncreaseNumber(kNSSSmall, 1, logstep);
      }
      return kTRUE;
   }
   // intercept down key
   else if ((EKeySym) keysym == kKey_Down) {
      // Get log step / alt key
      Bool_t logstep = fStepLog;
      if (event->fState & kKeyMod1Mask)
         logstep = !logstep;
      // shift-cntrl-down
      if ((event->fState & kKeyShiftMask) &&
          (event->fState & kKeyControlMask)) {
         IncreaseNumber(kNSSHuge, -1, logstep);
      }
      // cntrl-down
      else if (event->fState & kKeyControlMask) {
         IncreaseNumber(kNSSLarge, -1, logstep);
      }
      // shift-down
      else if (event->fState & kKeyShiftMask) {
         IncreaseNumber(kNSSMedium, -1, logstep);
      }
      // down
      else {
         IncreaseNumber(kNSSSmall, -1, logstep);
      }
      return kTRUE;
   }
   // intercept printable characters
   else if (n && (keysym < 127) && (keysym >= 32) &&
            ((EKeySym) keysym != kKey_Delete) &&
            ((EKeySym) keysym != kKey_Backspace) &&
            ((event->fState & kKeyControlMask) == 0)) {
      if (IsGoodChar(tmp[0], fNumStyle, fNumAttr)) {
         return TGTextEntry::HandleKey(event);
      } else {
         return kTRUE;
      }
   }
   // otherwise use default behaviour
   else {
      return TGTextEntry::HandleKey(event);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handle focus change.

Bool_t TGNumberEntryField::HandleFocusChange(Event_t * event)
{
   if (IsEnabled() && fNeedsVerification &&
       (event->fCode == kNotifyNormal) &&
       (event->fState != kNotifyPointer) && (event->fType == kFocusOut)) {
      // make sure we have a valid number by increasing it by 0
      IncreaseNumber(kNSSSmall, 0, kFALSE);
   }

   return TGTextEntry::HandleFocusChange(event);
}

////////////////////////////////////////////////////////////////////////////////
/// Text has changed message.

void TGNumberEntryField::TextChanged(const char *text)
{
   TGTextEntry::TextChanged(text);
   fNeedsVerification = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return was pressed.

void TGNumberEntryField::ReturnPressed()
{
   TString instr, outstr;
   instr = TGTextEntry::GetBuffer()->GetString();

   if (fNeedsVerification) {
      // make sure we have a valid number by increasing it by 0
      IncreaseNumber(kNSSSmall, 0, kFALSE);
   }
   outstr = TGTextEntry::GetBuffer()->GetString();
   if (instr != outstr) {
      InvalidInput(instr);
      gVirtualX->Bell(0);
   }
   TGTextEntry::ReturnPressed();
}

////////////////////////////////////////////////////////////////////////////////
/// Layout.

void TGNumberEntryField::Layout()
{
   if (GetAlignment() == kTextRight) {
      End(kFALSE);
   } else {
      Home(kFALSE);
   }
}

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGNumberEntryLayout                                                  //
//                                                                      //
// Layout manager for number entry widget                               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
/// Layout the internal GUI elements in use.

void TGNumberEntryLayout::Layout()
{
   if (fBox == 0) {
      return;
   }
   UInt_t w = fBox->GetWidth();
   UInt_t h = fBox->GetHeight();
   UInt_t upw = 2 * h / 3;
   UInt_t uph = h / 2;
   Int_t upx = (w > h) ? (Int_t) w - (Int_t) upw : -1000;
   Int_t upy = 0;
   Int_t downx = (w > h) ? (Int_t) w - (Int_t) upw : -1000;
   Int_t downy = h / 2;
   UInt_t downw = upw;
   UInt_t downh = h - downy;
   UInt_t numw = (w > h) ? w - upw : w;
   UInt_t numh = h;
   if (fBox->GetNumberEntry())
      fBox->GetNumberEntry()->MoveResize(0, 0, numw, numh);
   if (fBox->GetButtonUp())
      fBox->GetButtonUp()->MoveResize(upx, upy, upw, uph);
   if (fBox->GetButtonDown())
      fBox->GetButtonDown()->MoveResize(downx, downy, downw, downh);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the default size of the numeric control box.

TGDimension TGNumberEntryLayout::GetDefaultSize() const
{
   return fBox->GetSize();
}



//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TRepeatTimer                                                         //
//                                                                      //
// Timer for numeric control box buttons.                               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

class TGRepeatFireButton;

////////////////////////////////////////////////////////////////////////////////

class TRepeatTimer : public TTimer {
private:
   TGRepeatFireButton *fButton;  // Fire button

public:
   TRepeatTimer(TGRepeatFireButton * button, Long_t ms)
    : TTimer(ms, kTRUE), fButton(button) { }
   virtual Bool_t Notify();
};



//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TRepeatFireButton                                                    //
//                                                                      //
// Picture button which fires repeatedly as long as the button is pressed //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

class TGRepeatFireButton : public TGPictureButton {
protected:
   TRepeatTimer             *fTimer;           // the timer
   Int_t                     fIgnoreNextFire;  // flag for skipping next
   TGNumberFormat::EStepSize fStep;            // increment/decrement step
   Bool_t                    fStepLog;         // logarithmic step flag
   Bool_t                    fDoLogStep;       // flag for using logarithmic step

   Bool_t IsEditableParent();

public:
   TGRepeatFireButton(const TGWindow *p, const TGPicture *pic,
                      Int_t id, Bool_t logstep)
    : TGPictureButton(p, pic, id), fTimer(0), fIgnoreNextFire(0),
       fStep(TGNumberFormat::kNSSSmall), fStepLog(logstep), fDoLogStep(logstep)
       { fEditDisabled = kEditDisable | kEditDisableGrab; }
   virtual ~TGRepeatFireButton() { delete fTimer; }

   virtual  Bool_t HandleButton(Event_t *event);
            void   FireButton();
   virtual  void   SetLogStep(Bool_t on = kTRUE) { fStepLog = on; }
};

////////////////////////////////////////////////////////////////////////////////
/// Return kTRUE if one of the parents is in edit mode.

Bool_t TGRepeatFireButton::IsEditableParent()
{
   TGWindow *parent = (TGWindow*)GetParent();

   while (parent && (parent != fClient->GetDefaultRoot())) {
      if (parent->IsEditable()) {
         return kTRUE;
      }
      parent = (TGWindow*)parent->GetParent();
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle messages for number entry widget according to the user input.

Bool_t TGRepeatFireButton::HandleButton(Event_t * event)
{
   const Int_t t0 = 200;
   if (fTip)
      fTip->Hide();

   // disable button handling while gui building
   if (IsEditableParent()) {
      return kTRUE;
   }

   if (fState == kButtonDisabled)
      return kTRUE;

   if (event->fType == kButtonPress) {
      // Get log step / alt key
      fDoLogStep = fStepLog;
      if (event->fState & kKeyMod1Mask)
         fDoLogStep = !fDoLogStep;
      if ((event->fState & kKeyShiftMask) &&
          (event->fState & kKeyControlMask)) {
         fStep = TGNumberFormat::kNSSHuge;
      } else if (event->fState & kKeyControlMask) {
         fStep = TGNumberFormat::kNSSLarge;
      } else if (event->fState & kKeyShiftMask) {
         fStep = TGNumberFormat::kNSSMedium;
      } else {
         fStep = TGNumberFormat::kNSSSmall;
      }
      SetState(kButtonDown);
      fIgnoreNextFire = 0;
      FireButton();
      fIgnoreNextFire = 2;

      if (fTimer == 0) {
         fTimer = new TRepeatTimer(this, t0);
      }
      fTimer->Reset();
      gSystem->AddTimer(fTimer);
   } else {
      SetState(kButtonUp);
      if (fTimer != 0) {
         fTimer->Remove();
         fTimer->SetTime(t0);
      }
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Process messages for fire button.

void TGRepeatFireButton::FireButton()
{
   if (fIgnoreNextFire <= 0) {
      SendMessage(fMsgWindow, MK_MSG(kC_COMMAND, kCM_BUTTON),
                  fWidgetId, (Long_t) fStep + (fDoLogStep ? 100 : 0));
   } else {
      fIgnoreNextFire--;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Notify when timer times out and reset the timer.

Bool_t TRepeatTimer::Notify()
{
   fButton->FireButton();
   Reset();
   if ((Long64_t)fTime > 20) fTime -= 10;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Constructs a numeric entry widget.

TGNumberEntry::TGNumberEntry(const TGWindow *parent,
                             Double_t val, Int_t wdigits, Int_t id,
                             EStyle style,
                             EAttribute attr,
                             ELimit limits, Double_t min, Double_t max)
 : TGCompositeFrame(parent, 10 * wdigits, 25), fButtonToNum(kTRUE)
{
   fWidgetId  = id;
   fMsgWindow = parent;
   fPicUp = fClient->GetPicture("arrow_up.xpm");
   if (!fPicUp)
      Error("TGNumberEntry", "arrow_up.xpm not found");
   fPicDown = fClient->GetPicture("arrow_down.xpm");
   if (!fPicDown)
      Error("TGNumberEntry", "arrow_down.xpm not found");

   // create gui elements
   fNumericEntry = new TGNumberEntryField(this, id, val, style, attr,
                                          limits, min, max);
   fNumericEntry->Connect("ReturnPressed()", "TGNumberEntry", this,
                          "ValueSet(Long_t=0)");
   fNumericEntry->Connect("ReturnPressed()", "TGNumberEntry", this,
                          "Modified()");
   fNumericEntry->Associate(fMsgWindow);
   AddFrame(fNumericEntry, 0);
   fButtonUp = new TGRepeatFireButton(this, fPicUp, 1,
                                      fNumericEntry->IsLogStep());
   fButtonUp->Associate(this);
   AddFrame(fButtonUp, 0);
   fButtonDown = new TGRepeatFireButton(this, fPicDown, 2,
                                        fNumericEntry->IsLogStep());
   fButtonDown->Associate(this);
   AddFrame(fButtonDown, 0);

   // resize
   Int_t h = fNumericEntry->GetDefaultHeight();
   Int_t charw = fNumericEntry->GetCharWidth("0123456789");
   Int_t w = charw * TMath::Abs(wdigits) / 10 + 8 + 2 * h / 3;
   SetLayoutManager(new TGNumberEntryLayout(this));
   MapSubwindows();
   Resize(w, h);
   fEditDisabled = kEditDisableLayout | kEditDisableHeight;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructs a numeric entry widget.

TGNumberEntry::~TGNumberEntry()
{
   gClient->FreePicture(fPicUp);
   gClient->FreePicture(fPicDown);

   Cleanup();
}

////////////////////////////////////////////////////////////////////////////////
/// Make w the window that will receive the generated messages.

void TGNumberEntry::Associate(const TGWindow *w)
{
   TGWidget::Associate(w);
   fNumericEntry->Associate(w);
}

////////////////////////////////////////////////////////////////////////////////
/// Set log steps.

void TGNumberEntry::SetLogStep(Bool_t on)
{
   fNumericEntry->SetLogStep(on);
   ((TGRepeatFireButton *) fButtonUp)->SetLogStep(fNumericEntry->IsLogStep());
   ((TGRepeatFireButton *) fButtonDown)->SetLogStep(fNumericEntry->IsLogStep());
}

////////////////////////////////////////////////////////////////////////////////
/// Set the active state.

void TGNumberEntry::SetState(Bool_t enable)
{
   if (enable) {
      fButtonUp->SetState(kButtonUp);
      fButtonDown->SetState(kButtonUp);
      fNumericEntry->SetState(kTRUE);
   } else {
      fButtonUp->SetState(kButtonDisabled);
      fButtonDown->SetState(kButtonDisabled);
      fNumericEntry->SetState(kFALSE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Send button messages to the number field (true) or parent widget (false).
/// When the message is sent to the parent widget, it is responsible to change
/// the numerical value accordingly. This can be useful to implement cursors
/// which move from data point to data point. For the message being sent
/// see ProcessMessage().

void TGNumberEntry::SetButtonToNum(Bool_t state)
{
   fButtonToNum = state;
}

////////////////////////////////////////////////////////////////////////////////
/// Process the up/down button messages. If fButtonToNum is false the
/// following message is sent: kC_COMMAND, kCM_BUTTON, widget id, param
/// param % 100 is the step size
/// param % 10000 / 100 != 0 indicates log step
/// param / 10000 != 0 indicates button down

Bool_t TGNumberEntry::ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2)
{
   switch (GET_MSG(msg)) {
   case kC_COMMAND:
      {
         if ((GET_SUBMSG(msg) == kCM_BUTTON) &&
             (parm1 >= 1) && (parm1 <= 2)) {
            if (fButtonToNum) {
               Int_t sign = (parm1 == 1) ? 1 : -1;
               EStepSize step = (EStepSize) (parm2 % 100);
               Bool_t logstep = (parm2 >= 100);
               fNumericEntry->IncreaseNumber(step, sign, logstep);
            } else {
               SendMessage(fMsgWindow, msg, fWidgetId,
                           10000 * (parm1 - 1) + parm2);
               ValueChanged(10000 * (parm1 - 1) + parm2);
            }
         // Emit a signal needed by pad editor
         ValueSet(10000 * (parm1 - 1) + parm2);
         Modified();
         }
         break;
      }
   }
   return kTRUE;
}


////////////////////////////////////////////////////////////////////////////////
/// Return layout manager.

TGLayoutManager *TGNumberEntry::GetLayoutManager() const
{
   TGNumberEntry *entry = (TGNumberEntry*)this;

   if (entry->fLayoutManager->IsA() != TGNumberEntryLayout::Class()) {
      entry->SetLayoutManager(new TGNumberEntryLayout(entry));
   }

   return entry->fLayoutManager;
}

////////////////////////////////////////////////////////////////////////////////
/// Emit ValueChanged(Long_t) signal. This signal is emitted when
/// fButtonToNum is false. The val has the following meaning:
/// val % 100 is the step size
/// val % 10000 / 100 != 0 indicates log step
/// val / 10000 != 0 indicates button down

void TGNumberEntry::ValueChanged(Long_t val)
{
   Emit("ValueChanged(Long_t)", val);
}

////////////////////////////////////////////////////////////////////////////////
/// Emit ValueSet(Long_t) signal. This signal is emitted when the
/// number entry value is changed. The val has the following meaning:
/// val % 100 is the step size
/// val % 10000 / 100 != 0 indicates log step
/// val / 10000 != 0 indicates button down

void TGNumberEntry::ValueSet(Long_t val)
{
   Emit("ValueSet(Long_t)", val);
}

////////////////////////////////////////////////////////////////////////////////
/// Emit Modified() signal. This signal is emitted when the
/// number entry value is changed.

void TGNumberEntry::Modified()
{
   Emit("Modified()");
}

////////////////////////////////////////////////////////////////////////////////
/// Save a number entry widget as a C++ statement(s) on output stream out.

void TGNumberEntry::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   char quote = '"';

   // to calculate the digits parameter
   Int_t w = fNumericEntry->GetWidth();
   Int_t h = fNumericEntry->GetHeight();
   Int_t charw  = fNumericEntry->GetCharWidth("0123456789");
   Int_t digits = (30*w - 240 -20*h)/(3*charw) + 3;

   // for time format
   Int_t hour, min, sec;
   GetTime(hour, min, sec);

   // for date format
   Int_t yy, mm, dd;
   GetDate(yy, mm, dd);

   out << "   TGNumberEntry *";
   out << GetName() << " = new TGNumberEntry(" << fParent->GetName() << ", (Double_t) ";
   switch (GetNumStyle()){
      case kNESInteger:
         out << GetIntNumber() << "," << digits << "," << WidgetId()
             << ",(TGNumberFormat::EStyle) " << GetNumStyle();
         break;
      case kNESRealOne:
         out << GetNumber() << "," << digits << "," << WidgetId()
             << ",(TGNumberFormat::EStyle) " << GetNumStyle();
         break;
      case kNESRealTwo:
         out << GetNumber() << "," << digits << "," << WidgetId()
             << ",(TGNumberFormat::EStyle) " << GetNumStyle();
         break;
      case kNESRealThree:
         out << GetNumber() << "," << digits << "," << WidgetId()
             << ",(TGNumberFormat::EStyle) " << GetNumStyle();
         break;
      case kNESRealFour:
         out << GetNumber() << "," << digits << "," << WidgetId()
             << ",(TGNumberFormat::EStyle) " << GetNumStyle();
         break;
      case kNESReal:
         out << GetNumber() << "," << digits << "," << WidgetId()
             << ",(TGNumberFormat::EStyle) " << GetNumStyle();
         break;
      case kNESDegree:
         out << GetIntNumber() << "," << digits << "," << WidgetId()
             << ",(TGNumberFormat::EStyle) " << GetNumStyle();
         break;
      case kNESMinSec:
         out << min*60 + sec << "," << digits << "," << WidgetId()
             << ",(TGNumberFormat::EStyle) " << GetNumStyle();
         break;
      case kNESHourMin:
         out << hour*60 + min << "," << digits << "," << WidgetId()
             << ",(TGNumberFormat::EStyle) " << GetNumStyle();
         break;
      case kNESHourMinSec:
         out << hour*3600 + min*60 + sec << "," << digits << "," << WidgetId()
             << ",(TGNumberFormat::EStyle) " << GetNumStyle();
         break;
      case kNESDayMYear:
         out << yy << mm << dd << "," << digits << "," << WidgetId()
             << ",(TGNumberFormat::EStyle) " << GetNumStyle();
         break;
      case kNESMDayYear:
         out << yy << mm << dd << "," << digits << "," << WidgetId()
             << ",(TGNumberFormat::EStyle) " << GetNumStyle();
         break;
      case kNESHex:
      {  char hex[256];
         ULong_t l = GetHexNumber();
         IntToHexStr(hex, l);
         std::ios::fmtflags f = out.flags(); // store flags
         out << "0x" << std::hex << "U," << digits << "," << WidgetId()
             << ",(TGNumberFormat::EStyle) " << GetNumStyle();
         out.flags( f ); // restore flags (reset std::hex)
         break;
      }
   }
   if (GetNumMax() ==1) {
      if (GetNumMin() == 0) {
         if (GetNumLimits() == kNELNoLimits) {
            if (GetNumAttr() == kNEAAnyNumber) {
               out << ");" << std::endl;
            } else {
               out << ",(TGNumberFormat::EAttribute) " << GetNumAttr() << ");" << std::endl;
            }
         } else {
            out << ",(TGNumberFormat::EAttribute) " << GetNumAttr()
                << ",(TGNumberFormat::ELimit) " << GetNumLimits() << ");" << std::endl;
         }
      } else {
         out << ",(TGNumberFormat::EAttribute) " << GetNumAttr()
             << ",(TGNumberFormat::ELimit) " << GetNumLimits()
             << "," << GetNumMin() << ");" << std::endl;
      }
   } else {
         out << ",(TGNumberFormat::EAttribute) " << GetNumAttr()
             << ",(TGNumberFormat::ELimit) " << GetNumLimits()
             << "," << GetNumMin() << "," << GetNumMax() << ");" << std::endl;
   }
   if (option && strstr(option, "keep_names"))
      out << "   " << GetName() << "->SetName(\"" << GetName() << "\");" << std::endl;
   if (fButtonDown->GetState() == kButtonDisabled)
      out << "   " << GetName() << "->SetState(kFALSE);" << std::endl;

   TGToolTip *tip = GetNumberEntry()->GetToolTip();
   if (tip) {
      TString tiptext = tip->GetText()->GetString();
      tiptext.ReplaceAll("\n", "\\n");
      out << "   ";
      out << GetName() << "->GetNumberEntry()->SetToolTipText(" << quote
          << tiptext << quote << ");"  << std::endl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Save a number entry widget as a C++ statement(s) on output stream out.

void TGNumberEntryField::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   char quote = '"';

   // for time format
   Int_t hour, min, sec;
   GetTime(hour, min, sec);

   // for date format
   Int_t yy, mm, dd;
   GetDate(yy, mm, dd);

   out << "   TGNumberEntryField *";
   out << GetName() << " = new TGNumberEntryField(" << fParent->GetName()
       << ", " << WidgetId() << ", (Double_t) ";
   switch (GetNumStyle()){
      case kNESInteger:
         out << GetIntNumber()
             << ",(TGNumberFormat::EStyle) " << GetNumStyle();
         break;
      case kNESRealOne:
         out << GetNumber()
             << ",(TGNumberFormat::EStyle) " << GetNumStyle();
         break;
      case kNESRealTwo:
         out << GetNumber()
             << ",(TGNumberFormat::EStyle) " << GetNumStyle();
         break;
      case kNESRealThree:
         out << GetNumber()
             << ",(TGNumberFormat::EStyle) " << GetNumStyle();
         break;
      case kNESRealFour:
         out << GetNumber()
             << ",(TGNumberFormat::EStyle) " << GetNumStyle();
         break;
      case kNESReal:
         out << GetNumber()
             << ",(TGNumberFormat::EStyle) " << GetNumStyle();
         break;
      case kNESDegree:
         out << GetIntNumber()
             << ",(TGNumberFormat::EStyle) " << GetNumStyle();
         break;
      case kNESMinSec:
         out << min*60 + sec
             << ",(TGNumberFormat::EStyle) " << GetNumStyle();
         break;
      case kNESHourMin:
         out << hour*60 + min
             << ",(TGNumberFormat::EStyle) " << GetNumStyle();
         break;
      case kNESHourMinSec:
         out << hour*3600 + min*60 + sec
             << ",(TGNumberFormat::EStyle) " << GetNumStyle();
         break;
      case kNESDayMYear:
         out << yy << mm << dd
             << ",(TGNumberFormat::EStyle) " << GetNumStyle();
         break;
      case kNESMDayYear:
         out << yy << mm << dd
             << ",(TGNumberFormat::EStyle) " << GetNumStyle();
         break;
      case kNESHex:
      {  char hex[256];
         ULong_t l = GetHexNumber();
         IntToHexStr(hex, l);
         std::ios::fmtflags f = out.flags(); // store flags
         out << "0x" << std::hex << "U"
             << ",(TGNumberFormat::EStyle) " << GetNumStyle();
         out.flags( f ); // restore flags (reset std::hex)
         break;
      }
   }
   if (GetNumMax() ==1) {
      if (GetNumMin() == 0) {
         if (GetNumLimits() == kNELNoLimits) {
            if (GetNumAttr() == kNEAAnyNumber) {
               out << ");" << std::endl;
            } else {
               out << ",(TGNumberFormat::EAttribute) " << GetNumAttr() << ");" << std::endl;
            }
         } else {
            out << ",(TGNumberFormat::EAttribute) " << GetNumAttr()
                << ",(TGNumberFormat::ELimit) " << GetNumLimits() << ");" << std::endl;
         }
      } else {
         out << ",(TGNumberFormat::EAttribute) " << GetNumAttr()
             << ",(TGNumberFormat::ELimit) " << GetNumLimits()
             << "," << GetNumMin() << ");" << std::endl;
      }
   } else {
         out << ",(TGNumberFormat::EAttribute) " << GetNumAttr()
             << ",(TGNumberFormat::ELimit) " << GetNumLimits()
             << "," << GetNumMin() << "," << GetNumMax() << ");" << std::endl;
   }
   if (option && strstr(option, "keep_names"))
      out << "   " << GetName() << "->SetName(\"" << GetName() << "\");" << std::endl;
   if (!IsEnabled())
      out << "   " << GetName() << "->SetState(kFALSE);" << std::endl;

   out << "   " << GetName() << "->Resize("<< GetWidth() << "," << GetName()
       << "->GetDefaultHeight());" << std::endl;

   TGToolTip *tip = GetToolTip();
   if (tip) {
      TString tiptext = tip->GetText()->GetString();
      tiptext.ReplaceAll("\n", "\\n");
      out << "   ";
      out << GetName() << "->SetToolTipText(" << quote
          << tiptext << quote << ");"  << std::endl;
   }
}
