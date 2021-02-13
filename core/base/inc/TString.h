// @(#)root/base:$Id$
// Author: Fons Rademakers   04/08/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TString
#define ROOT_TString


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TString                                                              //
//                                                                      //
// Basic string class.                                                  //
//                                                                      //
// Cannot be stored in a TCollection... use TObjString instead.         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "Rtypes.h"

#include "TMathBase.h"
#include "ROOT/RStringView.hxx"
#include "ROOT/TypeTraits.hxx"
#include "snprintf.h"

#include <iosfwd>
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <string>

class TRegexp;
class TPRegexp;
class TString;
class TSubString;
class TObjArray;
class TVirtualMutex;
class TBuffer;
class TClass;
class TBufferFile;

TString operator+(const TString &s1, const TString &s2);
TString operator+(const TString &s,  const char *cs);
TString operator+(const char *cs, const TString &s);
TString operator+(const TString &s, char c);
TString operator+(char c, const TString &s);
Bool_t  operator==(const TString &s1, const TString &s2);
Bool_t  operator==(const TString &s1, const char *s2);
Bool_t  operator==(const TSubString &s1, const TSubString &s2);
Bool_t  operator==(const TSubString &s1, const TString &s2);
Bool_t  operator==(const TSubString &s1, const char *s2);
/*
template<class T>
struct is_signed_numeral : std::integral_constant<bool,
   std::is_integral<T>::value && std::is_signed<T>::value
> {};

template<class T>
struct is_unsigned_numeral : std::integral_constant<bool,
   std::is_integral<T>::value && !std::is_signed<T>::value
> {};

template<class T>
using is_float_numeral = std::is_floating_point<T>;
*/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  TSubString                                                          //
//                                                                      //
//  The TSubString class allows selected elements to be addressed.      //
//  There are no public constructors.                                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class TSubString {

friend class TStringLong;
friend class TString;

friend Bool_t operator==(const TSubString &s1, const TSubString &s2);
friend Bool_t operator==(const TSubString &s1, const TString &s2);
friend Bool_t operator==(const TSubString &s1, const char *s2);

private:
   TString      &fStr;           // Referenced string
   Ssiz_t        fBegin;         // Index of starting character
   Ssiz_t        fExtent;        // Length of TSubString

   // NB: the only constructor is private
   TSubString(const TString &s, Ssiz_t start, Ssiz_t len);

protected:
   void          SubStringError(Ssiz_t, Ssiz_t, Ssiz_t) const;
   void          AssertElement(Ssiz_t i) const;  // Verifies i is valid index

public:
   TSubString(const TSubString &s)
     : fStr(s.fStr), fBegin(s.fBegin), fExtent(s.fExtent) { }

   TSubString   &operator=(const char *s);       // Assignment from a char*
   TSubString   &operator=(const TString &s);    // Assignment from a TString
   TSubString   &operator=(const TSubString &s); // Assignment from a TSubString
   char         &operator()(Ssiz_t i);           // Index with optional bounds checking
   char         &operator[](Ssiz_t i);           // Index with bounds checking
   char          operator()(Ssiz_t i) const;     // Index with optional bounds checking
   char          operator[](Ssiz_t i) const;     // Index with bounds checking

   operator std::string_view() const { return std::string_view(Data(),fExtent); }
   operator std::string() const { return std::string(Data(),fExtent); }

   const char   *Data() const;
   Ssiz_t        Length() const          { return fExtent; }
   Ssiz_t        Start() const           { return fBegin; }
   TString&      String()                { return fStr; }
   void          ToLower();              // Convert self to lower-case
   void          ToUpper();              // Convert self to upper-case

   // For detecting null substrings
   Bool_t        IsNull() const          { return fBegin == kNPOS; }
   int           operator!() const       { return fBegin == kNPOS; }
};


//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  TString                                                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class TString {

friend class TStringLong;
friend class TSubString;
friend class TBufferFile;

friend TString operator+(const TString &s1, const TString &s2);
friend TString operator+(const TString &s,  const char *cs);
friend TString operator+(const char *cs, const TString &s);
friend TString operator+(const TString &s, char c);
friend TString operator+(char c, const TString &s);

template<class T>
friend typename std::enable_if<ROOT::TypeTraits::IsSignedNumeral<T>::value,TString>::type
operator+(TString s, T i);
template<class T>
friend typename std::enable_if<ROOT::TypeTraits::IsUnsignedNumeral<T>::value,TString>::type
operator+(TString s, T u);
template<class T>
friend typename std::enable_if<ROOT::TypeTraits::IsFloatNumeral<T>::value,TString>::type
operator+(TString s, T f);
template<class T>
friend typename std::enable_if<ROOT::TypeTraits::IsSignedNumeral<T>::value,TString>::type
operator+(T i, const TString &s);
template<class T>
friend typename std::enable_if<ROOT::TypeTraits::IsUnsignedNumeral<T>::value,TString>::type
operator+(T u, const TString &s);
template<class T>
friend typename std::enable_if<ROOT::TypeTraits::IsFloatNumeral<T>::value,TString>::type
operator+(T f, const TString &s);

friend Bool_t  operator==(const TString &s1, const TString &s2);
friend Bool_t  operator==(const TString &s1, const char *s2);

private:
#ifdef R__BYTESWAP
   enum { kShortMask = 0x01, kLongMask  = 0x1 };
#else
   enum { kShortMask = 0x80, kLongMask  = 0x80000000 };
#endif

   struct LongStr_t
   {
      Ssiz_t    fCap;    // Max string length (including null)
      Ssiz_t    fSize;   // String length (excluding null)
      char     *fData;   // Long string data
   };

   enum { kMinCap = (sizeof(LongStr_t) - 1)/sizeof(char) > 2 ?
                    (sizeof(LongStr_t) - 1)/sizeof(char) : 2 };

   struct ShortStr_t
   {
      unsigned char fSize;           // String length (excluding null)
      char          fData[kMinCap];  // Short string data
   };

   union UStr_t { LongStr_t fL; ShortStr_t fS; };

   enum { kNwords = sizeof(UStr_t) / sizeof(Ssiz_t)};

   struct RawStr_t
   {
      Ssiz_t fWords[kNwords];
   };

   struct Rep_t
   {
      union
      {
         LongStr_t  fLong;
         ShortStr_t fShort;
         RawStr_t   fRaw;
      };
   };

protected:
#ifndef __CINT__
   Rep_t          fRep;           //! String data
#endif

   // Special concatenation constructor
   TString(const char *a1, Ssiz_t n1, const char *a2, Ssiz_t n2);
   void           AssertElement(Ssiz_t nc) const; // Index in range
   void           Clobber(Ssiz_t nc);             // Remove old contents
   void           InitChar(char c);               // Initialize from char

   enum { kAlignment = 16 };
   static Ssiz_t  Align(Ssiz_t s) { return (s + (kAlignment-1)) & ~(kAlignment-1); }
   static Ssiz_t  Recommend(Ssiz_t s) { return (s < kMinCap ? kMinCap : Align(s+1)) - 1; }
   static Ssiz_t  AdjustCapacity(Ssiz_t oldCap, Ssiz_t newCap);

private:
   Bool_t         IsLong() const { return Bool_t(fRep.fShort.fSize & kShortMask); }
#ifdef R__BYTESWAP
   void           SetShortSize(Ssiz_t s) { fRep.fShort.fSize = (unsigned char)(s << 1); }
   Ssiz_t         GetShortSize() const { return fRep.fShort.fSize >> 1; }
#else
   void           SetShortSize(Ssiz_t s) { fRep.fShort.fSize = (unsigned char)s; }
   Ssiz_t         GetShortSize() const { return fRep.fShort.fSize; }
#endif
   void           SetLongSize(Ssiz_t s) { fRep.fLong.fSize = s; }
   Ssiz_t         GetLongSize() const { return fRep.fLong.fSize; }
   void           SetSize(Ssiz_t s) { IsLong() ? SetLongSize(s) : SetShortSize(s); }
   void           SetLongCap(Ssiz_t s) { fRep.fLong.fCap = kLongMask | s; }
   Ssiz_t         GetLongCap() const { return fRep.fLong.fCap & ~kLongMask; }
   void           SetLongPointer(char *p) { fRep.fLong.fData = p; }
   char          *GetLongPointer() { return fRep.fLong.fData; }
   const char    *GetLongPointer() const { return fRep.fLong.fData; }
   char          *GetShortPointer() { return fRep.fShort.fData; }
   const char    *GetShortPointer() const { return fRep.fShort.fData; }
   char          *GetPointer() { return IsLong() ? GetLongPointer() : GetShortPointer(); }
   const char    *GetPointer() const { return IsLong() ? GetLongPointer() : GetShortPointer(); }
#ifdef R__BYTESWAP
   static Ssiz_t  MaxSize() { return kMaxInt - 1; }
#else
   static Ssiz_t  MaxSize() { return (kMaxInt >> 1) - 1; }
#endif
   void           UnLink() const { if (IsLong()) delete [] fRep.fLong.fData; }
   void           Zero() {
      Ssiz_t (&a)[kNwords] = fRep.fRaw.fWords;
      for (UInt_t i = 0; i < kNwords; ++i)
         a[i] = 0;
   }
   char          *Init(Ssiz_t capacity, Ssiz_t nchar);
   void           Clone(Ssiz_t nc); // Make self a distinct copy w. capacity nc
   void           FormImp(const char *fmt, va_list ap);
   UInt_t         HashCase() const;
   UInt_t         HashFoldCase() const;

public:
   enum EStripType   { kLeading = 0x1, kTrailing = 0x2, kBoth = 0x3 };
   enum ECaseCompare { kExact, kIgnoreCase };
   static const Ssiz_t kNPOS = ::kNPOS;

   TString();                           // Null string
   explicit TString(Ssiz_t ic);         // Suggested capacity
   TString(const TString &s);           // Copy constructor
   TString(TString &&s) noexcept;       // Move constructor
   TString(const char *s);              // Copy to embedded null
   TString(const char *s, Ssiz_t n);    // Copy past any embedded nulls
   TString(const std::string &s);
   TString(char c);
   TString(char c, Ssiz_t s);
   explicit TString(const std::string_view &sub);
   TString(const TSubString &sub);

   virtual ~TString();

   // ROOT I/O interface
   virtual void     FillBuffer(char *&buffer) const;
   virtual void     ReadBuffer(char *&buffer);
   virtual Int_t    Sizeof() const;

   static TString  *ReadString(TBuffer &b, const TClass *clReq);
   static void      WriteString(TBuffer &b, const TString *a);

   friend TBuffer &operator<<(TBuffer &b, const TString *obj);

   // C I/O interface
   Bool_t   Gets(FILE *fp, Bool_t chop=kTRUE);
   void     Puts(FILE *fp);

   // Type conversion
   operator const char*() const { return GetPointer(); }
#if (__cplusplus >= 201700L) && (!defined(__clang_major__) || __clang_major__ > 5)
   // Clang 5.0 support for explicit conversion is still inadequate even in c++17 mode.
   // (It leads to extraneous ambiguous overload errors)
   explicit operator std::string() const { return std::string(GetPointer(),Length()); }
   explicit operator ROOT::Internal::TStringView() const { return ROOT::Internal::TStringView(GetPointer(),Length()); }
   operator std::string_view() const { return std::string_view(GetPointer(),Length()); }
#else
   operator ROOT::Internal::TStringView() const { return ROOT::Internal::TStringView(GetPointer(),Length()); }
#endif

   // Assignment
   TString    &operator=(char s);                // Replace string
   TString    &operator=(const char *s);
   TString    &operator=(const TString &s);
   TString    &operator=(TString &&s) noexcept;
   TString    &operator=(const std::string &s);
   TString    &operator=(const std::string_view &s);
   TString    &operator=(const TSubString &s);
   TString    &operator+=(const char *s);        // Append string
   TString    &operator+=(const TString &s);
   TString    &operator+=(char c);

   template<class T>
   typename std::enable_if<ROOT::TypeTraits::IsSignedNumeral<T>::value,TString>::type
              &operator+=(T i);
   template<class T>
   typename std::enable_if<ROOT::TypeTraits::IsUnsignedNumeral<T>::value,TString>::type
              &operator+=(T u);
   template<class T>
   typename std::enable_if<ROOT::TypeTraits::IsFloatNumeral<T>::value,TString>::type
              &operator+=(T f);

   // Indexing operators
   char         &operator[](Ssiz_t i);         // Indexing with bounds checking
   char         &operator()(Ssiz_t i);         // Indexing with optional bounds checking
   char          operator[](Ssiz_t i) const;
   char          operator()(Ssiz_t i) const;
   TSubString    operator()(Ssiz_t start, Ssiz_t len) const;   // Sub-string operator
   TSubString    operator()(const TRegexp &re) const;          // Match the RE
   TSubString    operator()(const TRegexp &re, Ssiz_t start) const;
   TSubString    operator()(TPRegexp &re) const;               // Match the Perl compatible Regular Expression
   TSubString    operator()(TPRegexp &re, Ssiz_t start) const;
   TSubString    SubString(const char *pat, Ssiz_t start = 0,
                           ECaseCompare cmp = kExact) const;

   // Non-static member functions
   TString     &Append(const char *cs);
   TString     &Append(const char *cs, Ssiz_t n);
   TString     &Append(const TString &s);
   TString     &Append(const TString &s, Ssiz_t n);
   TString     &Append(char c, Ssiz_t rep = 1);   // Append c rep times
   Int_t        Atoi() const;
   Long64_t     Atoll() const;
   Double_t     Atof() const;
   Bool_t       BeginsWith(const char *s,      ECaseCompare cmp = kExact) const;
   Bool_t       BeginsWith(const TString &pat, ECaseCompare cmp = kExact) const;
   Ssiz_t       Capacity() const { return (IsLong() ? GetLongCap() : kMinCap) - 1; }
   Ssiz_t       Capacity(Ssiz_t n);
   TString     &Chop();
   void         Clear();
   int          CompareTo(const char *cs,    ECaseCompare cmp = kExact) const;
   int          CompareTo(const TString &st, ECaseCompare cmp = kExact) const;
   Bool_t       Contains(const char *pat,    ECaseCompare cmp = kExact) const;
   Bool_t       Contains(const TString &pat, ECaseCompare cmp = kExact) const;
   Bool_t       Contains(const TRegexp &pat) const;
   Bool_t       Contains(TPRegexp &pat) const;
   Int_t        CountChar(Int_t c) const;
   TString      Copy() const;
   const char  *Data() const { return GetPointer(); }
   Bool_t       EndsWith(const char *pat, ECaseCompare cmp = kExact) const;
   Bool_t       EqualTo(const char *cs,    ECaseCompare cmp = kExact) const;
   Bool_t       EqualTo(const TString &st, ECaseCompare cmp = kExact) const;
   Ssiz_t       First(char c) const;
   Ssiz_t       First(const char *cs) const;
   void         Form(const char *fmt, ...)
#if defined(__GNUC__) && !defined(__CINT__)
   __attribute__((format(printf, 2, 3)))   /* 1 is the this pointer */
#endif
   ;
   UInt_t       Hash(ECaseCompare cmp = kExact) const;
   Ssiz_t       Index(const char *pat, Ssiz_t i = 0,
                      ECaseCompare cmp = kExact) const;
   Ssiz_t       Index(const TString &s, Ssiz_t i = 0,
                      ECaseCompare cmp = kExact) const;
   Ssiz_t       Index(const char *pat, Ssiz_t patlen, Ssiz_t i,
                      ECaseCompare cmp) const;
   Ssiz_t       Index(const TString &s, Ssiz_t patlen, Ssiz_t i,
                      ECaseCompare cmp) const;
   Ssiz_t       Index(const TRegexp &pat, Ssiz_t i = 0) const;
   Ssiz_t       Index(const TRegexp &pat, Ssiz_t *ext, Ssiz_t i = 0) const;
   Ssiz_t       Index(TPRegexp &pat, Ssiz_t i = 0) const;
   Ssiz_t       Index(TPRegexp &pat, Ssiz_t *ext, Ssiz_t i = 0) const;
   TString     &Insert(Ssiz_t pos, const char *s);
   TString     &Insert(Ssiz_t pos, const char *s, Ssiz_t extent);
   TString     &Insert(Ssiz_t pos, const TString &s);
   TString     &Insert(Ssiz_t pos, const TString &s, Ssiz_t extent);
   Bool_t       IsAscii() const;
   Bool_t       IsAlpha() const;
   Bool_t       IsAlnum() const;
   Bool_t       IsDigit() const;
   Bool_t       IsFloat() const;
   Bool_t       IsHex() const;
   Bool_t       IsBin() const;
   Bool_t       IsOct() const;
   Bool_t       IsDec() const;
   Bool_t       IsInBaseN(Int_t base) const;
   Bool_t       IsNull() const         { return Length() == 0; }
   Bool_t       IsWhitespace() const   { return (Length() == CountChar(' ')); }
   Ssiz_t       Last(char c) const;
   Ssiz_t       Length() const         { return IsLong() ? GetLongSize() : GetShortSize(); }
   Bool_t       MaybeRegexp() const;
   Bool_t       MaybeWildcard() const;
   TString      MD5() const;
   TString     &Prepend(const char *cs);     // Prepend a character string
   TString     &Prepend(const char *cs, Ssiz_t n);
   TString     &Prepend(const TString &s);
   TString     &Prepend(const TString &s, Ssiz_t n);
   TString     &Prepend(char c, Ssiz_t rep = 1);  // Prepend c rep times
   std::istream     &ReadFile(std::istream &str);      // Read to EOF or null character
   std::istream     &ReadLine(std::istream &str,
                         Bool_t skipWhite = kTRUE);   // Read to EOF or newline
   std::istream     &ReadString(std::istream &str);             // Read to EOF or null character
   std::istream     &ReadToDelim(std::istream &str, char delim = '\n'); // Read to EOF or delimitor
   std::istream     &ReadToken(std::istream &str);                // Read separated by white space
   TString     &Remove(Ssiz_t pos);                     // Remove pos to end of string
   TString     &Remove(Ssiz_t pos, Ssiz_t n);           // Remove n chars starting at pos
   TString     &Remove(EStripType s, char c);           // Like Strip() but changing string directly
   TString     &Replace(Ssiz_t pos, Ssiz_t n, const char *s);
   TString     &Replace(Ssiz_t pos, Ssiz_t n, const char *s, Ssiz_t ns);
   TString     &Replace(Ssiz_t pos, Ssiz_t n, const TString &s);
   TString     &Replace(Ssiz_t pos, Ssiz_t n1, const TString &s, Ssiz_t n2);
   TString     &ReplaceAll(const TString &s1, const TString &s2); // Find&Replace all s1 with s2 if any
   TString     &ReplaceAll(const TString &s1, const char *s2);    // Find&Replace all s1 with s2 if any
   TString     &ReplaceAll(const    char *s1, const TString &s2); // Find&Replace all s1 with s2 if any
   TString     &ReplaceAll(const char *s1, const char *s2);       // Find&Replace all s1 with s2 if any
   TString     &ReplaceAll(const char *s1, Ssiz_t ls1, const char *s2, Ssiz_t ls2);  // Find&Replace all s1 with s2 if any
   void         Resize(Ssiz_t n);                       // Truncate or add blanks as necessary
   TSubString   Strip(EStripType s = kTrailing, char c = ' ') const;
   TString     &Swap(TString &other); // Swap the contents of this and other without reallocation
   void         ToLower();                              // Change self to lower-case
   void         ToUpper();                              // Change self to upper-case
   TObjArray   *Tokenize(const TString &delim) const;
   Bool_t       Tokenize(TString &tok, Ssiz_t &from, const char *delim = " ") const;
   std::string_view View() const { return std::string_view(GetPointer(),Length()); }

   // Static member functions
   static UInt_t  Hash(const void *txt, Int_t ntxt);    // Calculates hash index from any char string.
   static Ssiz_t  InitialCapacity(Ssiz_t ic = 15);      // Initial allocation capacity
   static Ssiz_t  MaxWaste(Ssiz_t mw = 15);             // Max empty space before reclaim
   static Ssiz_t  ResizeIncrement(Ssiz_t ri = 16);      // Resizing increment
   static Ssiz_t  GetInitialCapacity();
   static Ssiz_t  GetResizeIncrement();
   static Ssiz_t  GetMaxWaste();
   static TString Itoa   (    Int_t value, Int_t base);  // Converts int to string with respect to the base specified (2-36)
   static TString UItoa  (   UInt_t value, Int_t base);
   static TString LLtoa  ( Long64_t value, Int_t base);
   static TString ULLtoa (ULong64_t value, Int_t base);
   static TString BaseConvert(const TString& s_in, Int_t base_in, Int_t base_out);  // Converts string from base base_in to base base_out (supported bases 2-36)
   static TString Format(const char *fmt, ...)
#if defined(__GNUC__) && !defined(__CINT__)
   __attribute__((format(printf, 1, 2)))
#endif
   ;

   ClassDef(TString,2)  //Basic string class
};

// Related global functions
std::istream  &operator>>(std::istream &str,       TString &s);
std::ostream  &operator<<(std::ostream &str, const TString &s);
#if defined(R__TEMPLATE_OVERLOAD_BUG)
template <>
#endif
TBuffer  &operator>>(TBuffer &buf,       TString *&sp);
TBuffer  &operator<<(TBuffer &buf, const TString * sp);

// Conversion operator (per se).
inline std::string& operator+=(std::string &left, const TString &right)
{
   return left.append(right.Data());
}

TString ToLower(const TString &s);    // Return lower-case version of argument
TString ToUpper(const TString &s);    // Return upper-case version of argument

inline UInt_t Hash(const TString &s) { return s.Hash(); }
inline UInt_t Hash(const TString *s) { return s->Hash(); }
       UInt_t Hash(const char *s);

extern char *Form(const char *fmt, ...)      // format in circular buffer
#if defined(__GNUC__) && !defined(__CINT__)
__attribute__((format(printf, 1, 2)))
#endif
;
extern void  Printf(const char *fmt, ...)    // format and print
#if defined(__GNUC__) && !defined(__CINT__)
__attribute__((format(printf, 1, 2)))
#endif
;
extern char *Strip(const char *str, char c = ' '); // strip c off str, free with delete []
extern char *StrDup(const char *str);        // duplicate str, free with delete []
extern char *Compress(const char *str);      // remove blanks from string, free with delele []
extern int   EscChar(const char *src, char *dst, int dstlen, char *specchars,
                     char escchar);          // copy from src to dst escaping specchars by escchar
extern int   UnEscChar(const char *src, char *dst, int dstlen, char *specchars,
                       char escchar);        // copy from src to dst removing escchar from specchars

#ifdef NEED_STRCASECMP
extern int strcasecmp(const char *str1, const char *str2);
extern int strncasecmp(const char *str1, const char *str2, Ssiz_t n);
#endif

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  Inlines                                                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

template<class T>
inline typename std::enable_if<ROOT::TypeTraits::IsSignedNumeral<T>::value,TString>::type
operator+(TString s, T i)
{ return s += i; }

template<class T>
inline typename std::enable_if<ROOT::TypeTraits::IsUnsignedNumeral<T>::value,TString>::type
operator+(TString s, T u)
{ return s += u; }

template<class T>
inline typename std::enable_if<ROOT::TypeTraits::IsFloatNumeral<T>::value,TString>::type
operator+(TString s, T f)
{ return s += f; }

template<class T>
inline typename std::enable_if<ROOT::TypeTraits::IsSignedNumeral<T>::value,TString>::type
operator+(T i, const TString &s)
{
    char buffer[32];
    // coverity[secure_coding] Buffer is large enough (2^64 = 20 digits).
    snprintf(buffer, sizeof(buffer), "%lld", static_cast<Long64_t>(i));
    return TString(buffer, strlen(buffer), s.Data(), s.Length());
}

template<class T>
inline typename std::enable_if<ROOT::TypeTraits::IsUnsignedNumeral<T>::value,TString>::type
operator+(T u, const TString &s)
{
    char buffer[32];
    // coverity[secure_coding] Buffer is large enough (2^64 = 20 digits).
    snprintf(buffer, sizeof(buffer), "%llu", static_cast<ULong64_t>(u));
    return TString(buffer, strlen(buffer), s.Data(), s.Length());
}

template<class T>
inline typename std::enable_if<ROOT::TypeTraits::IsFloatNumeral<T>::value,TString>::type
operator+(T f, const TString &s)
{
    char buffer[32];
    // coverity[secure_coding] Buffer is large enough: width specified in format
    snprintf(buffer, sizeof(buffer), "%.17Lg", static_cast<LongDouble_t>(f));
    return TString(buffer, strlen(buffer), s.Data(), s.Length());
}

inline TString &TString::Append(const char *cs)
{ return Replace(Length(), 0, cs, cs ? strlen(cs) : 0); }

inline TString &TString::Append(const char *cs, Ssiz_t n)
{ return Replace(Length(), 0, cs, n); }

inline TString &TString::Append(const TString &s)
{ return Replace(Length(), 0, s.Data(), s.Length()); }

inline TString &TString::Append(const TString &s, Ssiz_t n)
{ return Replace(Length(), 0, s.Data(), TMath::Min(n, s.Length())); }

inline TString &TString::operator+=(const char *cs)
{ return Append(cs, cs ? strlen(cs) : 0); }

inline TString &TString::operator+=(const TString &s)
{ return Append(s.Data(), s.Length()); }

inline TString &TString::operator+=(char c)
{ return Append(c); }

template<class T>
inline typename std::enable_if<ROOT::TypeTraits::IsSignedNumeral<T>::value,TString>::type
&TString::operator+=(T i)
{
   char buffer[32];
   // coverity[secure_coding] Buffer is large enough (2^64 = 20 digits).
   snprintf(buffer, sizeof(buffer), "%lld", static_cast<Long64_t>(i));
   return operator+=(buffer);
}

template<class T>
inline typename std::enable_if<ROOT::TypeTraits::IsUnsignedNumeral<T>::value,TString>::type
&TString::operator+=(T u)
{
   char buffer[32];
   // coverity[secure_coding] Buffer is large enough (2^64 = 20 digits).
   snprintf(buffer, sizeof(buffer), "%llu", static_cast<ULong64_t>(u));
   return operator+=(buffer);
}

template<class T>
inline typename std::enable_if<ROOT::TypeTraits::IsFloatNumeral<T>::value,TString>::type
&TString::operator+=(T f)
{
   char buffer[32];
   // coverity[secure_coding] Buffer is large enough: width specified in format
   snprintf(buffer, sizeof(buffer), "%.17Lg", static_cast<LongDouble_t>(f));
   return operator+=(buffer);
}

inline Bool_t TString::BeginsWith(const char *s, ECaseCompare cmp) const
{ return Index(s, s ? strlen(s) : (Ssiz_t)0, (Ssiz_t)0, cmp) == 0; }

inline Bool_t TString::BeginsWith(const TString &pat, ECaseCompare cmp) const
{ return Index(pat.Data(), pat.Length(), (Ssiz_t)0, cmp) == 0; }

inline Bool_t TString::Contains(const TString &pat, ECaseCompare cmp) const
{ return Index(pat.Data(), pat.Length(), (Ssiz_t)0, cmp) != kNPOS; }

inline Bool_t TString::Contains(const char *s, ECaseCompare cmp) const
{ return Index(s, s ? strlen(s) : 0, (Ssiz_t)0, cmp) != kNPOS; }

inline Bool_t TString::Contains(const TRegexp &pat) const
{ return Index(pat, (Ssiz_t)0) != kNPOS; }

inline Bool_t TString::Contains(TPRegexp &pat) const
{ return Index(pat, (Ssiz_t)0) != kNPOS; }

inline Bool_t TString::EqualTo(const char *cs, ECaseCompare cmp) const
{ return (CompareTo(cs, cmp) == 0) ? kTRUE : kFALSE; }

inline Bool_t TString::EqualTo(const TString &st, ECaseCompare cmp) const
{ return (CompareTo(st, cmp) == 0) ? kTRUE : kFALSE; }

inline Ssiz_t TString::Index(const char *s, Ssiz_t i, ECaseCompare cmp) const
{ return Index(s, s ? strlen(s) : 0, i, cmp); }

inline Ssiz_t TString::Index(const TString &s, Ssiz_t i, ECaseCompare cmp) const
{ return Index(s.Data(), s.Length(), i, cmp); }

inline Ssiz_t TString::Index(const TString &pat, Ssiz_t patlen, Ssiz_t i,
                             ECaseCompare cmp) const
{ return Index(pat.Data(), patlen, i, cmp); }

inline TString &TString::Insert(Ssiz_t pos, const char *cs)
{ return Replace(pos, 0, cs, cs ? strlen(cs) : 0); }

inline TString &TString::Insert(Ssiz_t pos, const char *cs, Ssiz_t n)
{ return Replace(pos, 0, cs, n); }

inline TString &TString::Insert(Ssiz_t pos, const TString &s)
{ return Replace(pos, 0, s.Data(), s.Length()); }

inline TString &TString::Insert(Ssiz_t pos, const TString &s, Ssiz_t n)
{ return Replace(pos, 0, s.Data(), TMath::Min(n, s.Length())); }

inline TString &TString::Prepend(const char *cs)
{ return Replace(0, 0, cs, cs ? strlen(cs) : 0); }

inline TString &TString::Prepend(const char *cs, Ssiz_t n)
{ return Replace(0, 0, cs, n); }

inline TString &TString::Prepend(const TString &s)
{ return Replace(0, 0, s.Data(), s.Length()); }

inline TString &TString::Prepend(const TString &s, Ssiz_t n)
{ return Replace(0, 0, s.Data(), TMath::Min(n, s.Length())); }

inline TString &TString::Remove(Ssiz_t pos)
{ return Replace(pos, TMath::Max(0, Length()-pos), 0, 0); }

inline TString &TString::Remove(Ssiz_t pos, Ssiz_t n)
{ return Replace(pos, n, 0, 0); }

inline TString &TString::Chop()
{ return Remove(TMath::Max(0, Length()-1)); }

inline TString &TString::Replace(Ssiz_t pos, Ssiz_t n, const char *cs)
{ return Replace(pos, n, cs, cs ? strlen(cs) : 0); }

inline TString &TString::Replace(Ssiz_t pos, Ssiz_t n, const TString& s)
{ return Replace(pos, n, s.Data(), s.Length()); }

inline TString &TString::Replace(Ssiz_t pos, Ssiz_t n1, const TString &s,
                                 Ssiz_t n2)
{ return Replace(pos, n1, s.Data(), TMath::Min(s.Length(), n2)); }

inline TString &TString::ReplaceAll(const TString &s1, const TString &s2)
{ return ReplaceAll(s1.Data(), s1.Length(), s2.Data(), s2.Length()) ; }

inline TString &TString::ReplaceAll(const TString &s1, const char *s2)
{ return ReplaceAll(s1.Data(), s1.Length(), s2, s2 ? strlen(s2) : 0); }

inline TString &TString::ReplaceAll(const char *s1, const TString &s2)
{ return ReplaceAll(s1, s1 ? strlen(s1) : 0, s2.Data(), s2.Length()); }

inline TString &TString::ReplaceAll(const char *s1,const char *s2)
{ return ReplaceAll(s1, s1 ? strlen(s1) : 0, s2, s2 ? strlen(s2) : 0); }

inline TString &TString::Swap(TString &other) {
   // Swap the contents of other and this without reallocation.
#ifndef __CINT__
   Rep_t tmp = other.fRep;
   other.fRep = fRep;
   fRep = tmp;
#endif
   return *this;
}

inline char &TString::operator()(Ssiz_t i)
{ return GetPointer()[i]; }

inline char TString::operator()(Ssiz_t i) const
{ return GetPointer()[i]; }

inline char &TString::operator[](Ssiz_t i)
{ AssertElement(i); return GetPointer()[i]; }

inline char TString::operator[](Ssiz_t i) const
{ AssertElement(i); return GetPointer()[i]; }

inline const char *TSubString::Data() const
{
   // Return a pointer to the beginning of the substring. Note that the
   // terminating null is in the same place as for the original
   // TString, so this method is not appropriate for converting the
   // TSubString to a string. To do that, construct a TString from the
   // TSubString. For example:
   //
   //   root [0] TString s("hello world")
   //   root [1] TSubString sub=s(0, 5)
   //   root [2] sub.Data()
   //   (const char* 0x857c8b8)"hello world"
   //   root [3] TString substr(sub)
   //   root [4] substr
   //   (class TString)"hello"

   return fStr.Data() + fBegin;
}

// Access to elements of sub-string with bounds checking
inline char TSubString::operator[](Ssiz_t i) const
{ AssertElement(i); return fStr.GetPointer()[fBegin+i]; }

inline char TSubString::operator()(Ssiz_t i) const
{ return fStr.GetPointer()[fBegin+i]; }

inline TSubString &TSubString::operator=(const TSubString &s)
{ fStr = s.fStr; fBegin = s.fBegin; fExtent = s.fExtent; return *this; }


// String Logical operators
inline Bool_t operator==(const TString &s1, const TString &s2)
{
   return ((s1.Length() == s2.Length()) &&
            !memcmp(s1.Data(), s2.Data(), s1.Length()));
}

inline Bool_t operator!=(const TString &s1, const TString &s2)
{ return !(s1 == s2); }

inline Bool_t operator<(const TString &s1, const TString &s2)
{ return s1.CompareTo(s2) < 0; }

inline Bool_t operator>(const TString &s1, const TString &s2)
{ return s1.CompareTo(s2) > 0; }

inline Bool_t operator<=(const TString &s1, const TString &s2)
{ return s1.CompareTo(s2) <= 0; }

inline Bool_t operator>=(const TString &s1, const TString &s2)
{ return s1.CompareTo(s2) >= 0; }

//     Bool_t operator==(const TString &s1, const char *s2);
inline Bool_t operator!=(const TString &s1, const char *s2)
{ return !(s1 == s2); }

inline Bool_t operator<(const TString &s1, const char *s2)
{ return s1.CompareTo(s2) < 0; }

inline Bool_t operator>(const TString &s1, const char *s2)
{ return s1.CompareTo(s2) > 0; }

inline Bool_t operator<=(const TString &s1, const char *s2)
{ return s1.CompareTo(s2) <= 0; }

inline Bool_t operator>=(const TString &s1, const char *s2)
{ return s1.CompareTo(s2) >= 0; }

inline Bool_t operator==(const char *s1, const TString &s2)
{ return (s2 == s1); }

inline Bool_t operator!=(const char *s1, const TString &s2)
{ return !(s2 == s1); }

inline Bool_t operator<(const char *s1, const TString &s2)
{ return s2.CompareTo(s1) > 0; }

inline Bool_t operator>(const char *s1, const TString &s2)
{ return s2.CompareTo(s1) < 0; }

inline Bool_t operator<=(const char *s1, const TString &s2)
{ return s2.CompareTo(s1) >= 0; }

inline Bool_t operator>=(const char *s1, const TString &s2)
{ return s2.CompareTo(s1) <= 0; }

// SubString Logical operators
//     Bool_t operator==(const TSubString &s1, const TSubString &s2);
//     Bool_t operator==(const TSubString &s1, const char *s2);
//     Bool_t operator==(const TSubString &s1, const TString &s2);
inline Bool_t operator==(const TString &s1, const TSubString &s2)
{ return (s2 == s1); }

inline Bool_t operator==(const char *s1, const TSubString &s2)
{ return (s2 == s1); }

inline Bool_t operator!=(const TSubString &s1, const char *s2)
{ return !(s1 == s2); }

inline Bool_t operator!=(const TSubString &s1, const TString &s2)
{ return !(s1 == s2); }

inline Bool_t operator!=(const TSubString &s1, const TSubString &s2)
{ return !(s1 == s2); }

inline Bool_t operator!=(const TString &s1, const TSubString &s2)
{ return !(s2 == s1); }

inline Bool_t operator!=(const char *s1, const TSubString &s2)
{ return !(s2 == s1); }

#ifndef WIN32
// To avoid ambiguities.
inline Bool_t operator==(const char *s1, const std::string_view &s2)
{
  return std::string_view(s1) == s2;
}

inline Bool_t operator==(const std::string_view &s1, const char *s2)
{
  return s1 == std::string_view(s2);
}
#endif

namespace llvm {
   class raw_ostream;
}

namespace cling {
  std::string printValue(const TString* val);
  std::string printValue(const TSubString* val);
  std::string printValue(const std::string_view* val);
}

#endif
