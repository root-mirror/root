// @(#)root/base:$Id$
// Author: Eddy Offermann   24/06/05

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/* \class TPRegexp
\ingroup Base

C++ Wrapper for the "Perl Compatible Regular Expressions" library
 The PCRE lib can be found at: http://www.pcre.org/

Extensive documentation about Regular expressions in Perl can be
found at : http://perldoc.perl.org/perlre.html
*/

#include <iostream>
#include "TPRegexp.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TError.h"

#ifdef R__WIN32
#define PCRE_STATIC
#endif
#include <pcre.h>

#include <vector>
#include <stdexcept>

struct PCREPriv_t {
   pcre       *fPCRE;
   pcre_extra *fPCREExtra;

   PCREPriv_t() { fPCRE = 0; fPCREExtra = 0; }
};


ClassImp(TPRegexp);

Bool_t TPRegexp::fgThrowAtCompileError = kFALSE;

////////////////////////////////////////////////////////////////////////////////
/// Default ctor.

TPRegexp::TPRegexp()
{
   fPriv     = new PCREPriv_t;
   fPCREOpts = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Create and initialize with pat.

TPRegexp::TPRegexp(const TString &pat)
{
   fPattern  = pat;
   fPriv     = new PCREPriv_t;
   fPCREOpts = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy ctor.

TPRegexp::TPRegexp(const TPRegexp &p)
{
   fPattern  = p.fPattern;
   fPriv     = new PCREPriv_t;
   fPCREOpts = p.fPCREOpts;
}

////////////////////////////////////////////////////////////////////////////////
/// Cleanup.

TPRegexp::~TPRegexp()
{
   if (fPriv->fPCRE)
      pcre_free(fPriv->fPCRE);
   if (fPriv->fPCREExtra)
      pcre_free(fPriv->fPCREExtra);
   delete fPriv;
}

////////////////////////////////////////////////////////////////////////////////
/// Assignment operator.

TPRegexp &TPRegexp::operator=(const TPRegexp &p)
{
   if (this != &p) {
      fPattern = p.fPattern;
      if (fPriv->fPCRE)
         pcre_free(fPriv->fPCRE);
      fPriv->fPCRE = 0;
      if (fPriv->fPCREExtra)
         pcre_free(fPriv->fPCREExtra);
      fPriv->fPCREExtra = 0;
      fPCREOpts  = p.fPCREOpts;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Translate Perl modifier flags into pcre flags.
/// The supported modStr characters are: g, i, m, o, s, x, and the
/// special d for debug. The meaning of the letters is:
/// - m
///   Treat string as multiple lines. That is, change "^" and "$" from
///   matching the start or end of the string to matching the start or
///   end of any line anywhere within the string.
/// - s
///   Treat string as single line. That is, change "." to match any
///   character whatsoever, even a newline, which normally it would not match.
///   Used together, as /ms, they let the "." match any character whatsoever,
///   while still allowing "^" and "$" to match, respectively, just after and
///   just before newlines within the string.
/// - i
///   Do case-insensitive pattern matching.
/// - x
///   Extend your pattern's legibility by permitting whitespace and comments.
/// - p
///   Preserve the string matched such that ${^PREMATCH}, ${^MATCH},
///   and ${^POSTMATCH} are available for use after matching.
/// - g and c
///   Global matching, and keep the Current position after failed matching.
///   Unlike i, m, s and x, these two flags affect the way the regex is used
///   rather than the regex itself. See Using regular expressions in Perl in
///   perlretut for further explanation of the g and c modifiers.
/// For more detail see: http://perldoc.perl.org/perlre.html#Modifiers.

UInt_t TPRegexp::ParseMods(const TString &modStr) const
{
   UInt_t opts = 0;

   if (modStr.Length() <= 0)
      return fPCREOpts;

   //translate perl flags into pcre flags
   const char *m = modStr;
   while (*m) {
      switch (*m) {
         case 'g':
            opts |= kPCRE_GLOBAL;
            break;
         case 'i':
            opts |= PCRE_CASELESS;
            break;
         case 'm':
            opts |= PCRE_MULTILINE;
            break;
         case 'o':
            opts |= kPCRE_OPTIMIZE;
            break;
         case 's':
            opts |= PCRE_DOTALL;
            break;
         case 'x':
            opts |= PCRE_EXTENDED;
            break;
         case 'd': // special flag to enable debug printing (not Perl compat.)
            opts |= kPCRE_DEBUG_MSGS;
            break;
         default:
            Error("ParseMods", "illegal pattern modifier: %c", *m);
            opts = 0;
      }
      ++m;
   }
   return opts;
}

////////////////////////////////////////////////////////////////////////////////
/// Return PCRE modifier options as string.
/// For meaning of mods see ParseMods().

TString TPRegexp::GetModifiers() const
{
   TString ret;

   if (fPCREOpts & kPCRE_GLOBAL)     ret += 'g';
   if (fPCREOpts & PCRE_CASELESS)    ret += 'i';
   if (fPCREOpts & PCRE_MULTILINE)   ret += 'm';
   if (fPCREOpts & PCRE_DOTALL)      ret += 's';
   if (fPCREOpts & PCRE_EXTENDED)    ret += 'x';
   if (fPCREOpts & kPCRE_OPTIMIZE)   ret += 'o';
   if (fPCREOpts & kPCRE_DEBUG_MSGS) ret += 'd';

   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Compile the fPattern.

void TPRegexp::Compile()
{
   if (fPriv->fPCRE)
      pcre_free(fPriv->fPCRE);

   if (fPCREOpts & kPCRE_DEBUG_MSGS)
      Info("Compile", "PREGEX compiling %s", fPattern.Data());

   const char *errstr;
   Int_t patIndex;
   fPriv->fPCRE = pcre_compile(fPattern.Data(), fPCREOpts & kPCRE_INTMASK,
                               &errstr, &patIndex, 0);

   if (!fPriv->fPCRE) {
      if (fgThrowAtCompileError) {
         throw std::runtime_error
            (TString::Format("TPRegexp::Compile() compilation of TPRegexp(%s) failed at: %d because %s",
                             fPattern.Data(), patIndex, errstr).Data());
      } else {
         Error("Compile", "compilation of TPRegexp(%s) failed at: %d because %s",
               fPattern.Data(), patIndex, errstr);
         return;
      }
   }

   if (fPriv->fPCREExtra || (fPCREOpts & kPCRE_OPTIMIZE))
      Optimize();
}

////////////////////////////////////////////////////////////////////////////////
/// Send the pattern through the optimizer.

void TPRegexp::Optimize()
{
   if (fPriv->fPCREExtra)
      pcre_free(fPriv->fPCREExtra);

   if (fPCREOpts & kPCRE_DEBUG_MSGS)
      Info("Optimize", "PREGEX studying %s", fPattern.Data());

   const char *errstr;
   // pcre_study allows less options - see pcre_internal.h PUBLIC_STUDY_OPTIONS.
   fPriv->fPCREExtra = pcre_study(fPriv->fPCRE, 0, &errstr);

   if (!fPriv->fPCREExtra && errstr) {
      Error("Optimize", "Optimization of TPRegexp(%s) failed: %s",
            fPattern.Data(), errstr);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the number of expanded '$' constructs.

Int_t TPRegexp::ReplaceSubs(const TString &s, TString &final,
                            const TString &replacePattern,
                            Int_t *offVec, Int_t nrMatch) const
{
   Int_t nrSubs = 0;
   const char *p = replacePattern;

   Int_t state = 0;
   Int_t subnum = 0;
   while (state != -1) {
      switch (state) {
         case 0:
            if (!*p) {
               state = -1;
               break;
            }
            if (*p == '$') {
               state = 1;
               subnum = 0;
               if (p[1] == '&') {
                  p++;
                  if (isdigit(p[1]))
                     p++;
               } else if (!isdigit(p[1])) {
                  Error("ReplaceSubs", "badly formed replacement pattern: %s",
                        replacePattern.Data());
               }
            } else
               final += *p;
            break;
         case 1:
            if (isdigit(*p)) {
               subnum *= 10;
               subnum += (*p)-'0';
            } else {
               if (fPCREOpts & kPCRE_DEBUG_MSGS)
                  Info("ReplaceSubs", "PREGEX appending substr #%d", subnum);
               if (subnum < 0 || subnum > nrMatch-1) {
                  Error("ReplaceSubs","bad string number: %d",subnum);
               } else {
                  const TString subStr = s(offVec[2*subnum],offVec[2*subnum+1]-offVec[2*subnum]);
                  final += subStr;
                  nrSubs++;
               }
               state = 0;
               continue;  // send char to start state
            }
      }
      p++;
   }
   return nrSubs;
}

////////////////////////////////////////////////////////////////////////////////
/// Perform the actual matching - protected method.

Int_t TPRegexp::MatchInternal(const TString &s, Int_t start,
                              Int_t nMaxMatch, TArrayI *pos) const
{
   Int_t *offVec = new Int_t[3*nMaxMatch];
   // pcre_exec allows less options - see pcre_internal.h PUBLIC_EXEC_OPTIONS.
   Int_t nrMatch = pcre_exec(fPriv->fPCRE, fPriv->fPCREExtra, s.Data(),
                             s.Length(), start, 0,
                             offVec, 3*nMaxMatch);

   if (nrMatch == PCRE_ERROR_NOMATCH)
      nrMatch = 0;
   else if (nrMatch <= 0) {
      Error("Match","pcre_exec error = %d", nrMatch);
      delete [] offVec;
      return 0;
   }

   if (pos)
      pos->Set(2*nrMatch, offVec);
   delete [] offVec;

   return nrMatch;
}

////////////////////////////////////////////////////////////////////////////////
/// The number of matches is returned, this equals the full match +
/// sub-pattern matches.
/// nMaxMatch is the maximum allowed number of matches.
/// pos contains the string indices of the matches. Its usage is
/// shown in the routine MatchS.
/// For meaning of mods see ParseMods().

Int_t TPRegexp::Match(const TString &s, const TString &mods, Int_t start,
                      Int_t nMaxMatch, TArrayI *pos)
{
   UInt_t opts = ParseMods(mods);

   if (!fPriv->fPCRE || opts != fPCREOpts) {
      fPCREOpts = opts;
      Compile();
   }

   return MatchInternal(s, start, nMaxMatch, pos);
}


////////////////////////////////////////////////////////////////////////////////
/// Returns a TObjArray of matched substrings as TObjString's.
/// The TObjArray is owner of the objects and must be deleted by the user.
/// The first entry is the full matched pattern, followed by the sub-patterns.
/// If a pattern was not matched, it will return an empty substring:
/// ~~~ {.cpp}
/// TObjArray *subStrL = TPRegexp("(a|(z))(bc)").MatchS("abc");
/// for (Int_t i = 0; i < subStrL->GetLast()+1; i++) {
///    const TString subStr = ((TObjString *)subStrL->At(i))->GetString();
///    std::cout << "\"" << subStr << "\" ";
/// }
/// std::cout << subStr << std::endl;
/// ~~~
/// produces:  "abc" "a" "" "bc"
///
/// For meaning of mods see ParseMods().

TObjArray *TPRegexp::MatchS(const TString &s, const TString &mods,
                            Int_t start, Int_t nMaxMatch)
{
   TArrayI pos;
   Int_t nrMatch = Match(s, mods, start, nMaxMatch, &pos);

   TObjArray *subStrL = new TObjArray();
   subStrL->SetOwner();

   for (Int_t i = 0; i < nrMatch; i++) {
      Int_t startp = pos[2*i];
      Int_t stopp  = pos[2*i+1];
      if (startp >= 0 && stopp >= 0) {
         const TString subStr = s(pos[2*i], pos[2*i+1]-pos[2*i]);
         subStrL->Add(new TObjString(subStr));
      } else
         subStrL->Add(new TObjString());
   }

   return subStrL;
}

////////////////////////////////////////////////////////////////////////////////
/// Perform pattern substitution with optional back-ref replacement
/// - protected method.

Int_t TPRegexp::SubstituteInternal(TString &s, const TString &replacePattern,
                                   Int_t start, Int_t nMaxMatch,
                                   Bool_t doDollarSubst) const
{
   Int_t *offVec = new Int_t[3*nMaxMatch];

   TString final;
   Int_t nrSubs = 0;
   Int_t offset = start;
   Int_t last = 0;

   while (kTRUE) {

      // find next matching subs
      // pcre_exec allows less options - see pcre_internal.h PUBLIC_EXEC_OPTIONS.
      Int_t nrMatch = pcre_exec(fPriv->fPCRE, fPriv->fPCREExtra, s.Data(),
                                s.Length(), offset, 0,
                                offVec, 3*nMaxMatch);

      if (nrMatch == PCRE_ERROR_NOMATCH) {
         nrMatch = 0;
         break;
      } else if (nrMatch <= 0) {
         Error("Substitute", "pcre_exec error = %d", nrMatch);
         break;
      }

      // append anything previously unmatched, but not substituted
      if (last <= offVec[0]) {
         final += s(last,offVec[0]-last);
         last = offVec[1];
      }

      // replace stuff in s
      if (doDollarSubst) {
         ReplaceSubs(s, final, replacePattern, offVec, nrMatch);
      } else {
         final += replacePattern;
      }
      ++nrSubs;

      // if global gotta check match at every pos
      if (!(fPCREOpts & kPCRE_GLOBAL))
         break;

      if (offVec[0] != offVec[1])
         offset = offVec[1];
      else {
         // matched empty string
         if (offVec[1] == s.Length())
         break;
         offset = offVec[1]+1;
      }
   }

   delete [] offVec;

   final += s(last,s.Length()-last);
   s = final;

   return nrSubs;
}

////////////////////////////////////////////////////////////////////////////////
/// Substitute replaces the string s by a new string in which matching
/// patterns are replaced by the replacePattern string. The number of
/// substitutions are returned.
/// ~~~ {.cpp}
/// TString s("aap noot mies");
/// const Int_t nrSub = TPRegexp("(\\w*) noot (\\w*)").Substitute(s,"$2 noot $1");
/// std::cout << nrSub << " \"" << s << "\"" <<std::endl;
/// ~~~
/// produces: 2 "mies noot aap"
///
/// For meaning of mods see ParseMods().

Int_t TPRegexp::Substitute(TString &s, const TString &replacePattern,
                           const TString &mods, Int_t start, Int_t nMaxMatch)
{
   UInt_t opts = ParseMods(mods);

   if (!fPriv->fPCRE || opts != fPCREOpts) {
      fPCREOpts = opts;
      Compile();
   }

   return SubstituteInternal(s, replacePattern, start, nMaxMatch, kTRUE);
}


////////////////////////////////////////////////////////////////////////////////
/// Returns true if underlying PCRE structure has been successfully
/// generated via regexp compilation.

Bool_t TPRegexp::IsValid() const
{
   return fPriv->fPCRE != 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get value of static flag controlling whether exception should be thrown upon an
/// error during regular expression compilation by the PCRE engine.

Bool_t TPRegexp::GetThrowAtCompileError()
{
   return fgThrowAtCompileError;
}

////////////////////////////////////////////////////////////////////////////////
/// Set static flag controlling whether exception should be thrown upon an
/// error during regular expression compilation by the PCRE engine.

void TPRegexp::SetThrowAtCompileError(Bool_t throwp)
{
   fgThrowAtCompileError = throwp;
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// TString member functions, put here so the linker will include              //
// them only if regular expressions are used.                                 //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Find the first occurrence of the regexp in string and return the position.
/// Start is the offset at which the search should start.

Ssiz_t TString::Index(TPRegexp& r, Ssiz_t start) const
{
   TArrayI pos;
   Int_t nrMatch = r.Match(*this,"",start,10,&pos);
   if (nrMatch > 0)
      return pos[0];
   else
      return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Find the first occurrence of the regexp in string and return the position.
/// Extent is length of the matched string and start is the offset at which
/// the matching should start.

Ssiz_t TString::Index(TPRegexp& r, Ssiz_t* extent, Ssiz_t start) const
{
   TArrayI pos;
   const Int_t nrMatch = r.Match(*this,"",start,10,&pos);
   if (nrMatch > 0) {
      *extent = pos[1]-pos[0];
      return pos[0];
   } else {
      *extent = 0;
      return -1;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return the substring found by applying the regexp starting at start.

TSubString TString::operator()(TPRegexp& r, Ssiz_t start) const
{
   Ssiz_t len;
   Ssiz_t begin = Index(r, &len, start);
   return TSubString(*this, begin, len);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the substring found by applying the regexp.

TSubString TString::operator()(TPRegexp& r) const
{
   return (*this)(r, 0);
}


/** \class TPMERegexp

Wrapper for PCRE library (Perl Compatible Regular Expressions).
Based on PME - PCRE Made Easy by Zachary Hansen.

Supports main Perl operations using regular expressions (Match,
Substitute and Split). To retrieve the results one can simply use
operator[] returning a TString.

See $ROOTSYS/tutorials/regexp_pme.C for examples.
*/

ClassImp(TPMERegexp);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor. This regexp will match an empty string.

TPMERegexp::TPMERegexp() :
   TPRegexp(),
   fNMaxMatches(10),
   fNMatches(0),
   fAddressOfLastString(0),
   fLastGlobalPosition(0)
{
   Compile();
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.
///
/// \param[in] s      string to compile into regular expression
/// \param[in] opts   perl-style character flags to be set on TPME object
/// \param[in] nMatchMax  maximum number of matches

TPMERegexp::TPMERegexp(const TString& s, const TString& opts, Int_t nMatchMax) :
   TPRegexp(s),
   fNMaxMatches(nMatchMax),
   fNMatches(0),
   fAddressOfLastString(0),
   fLastGlobalPosition(0)
{
   fPCREOpts = ParseMods(opts);
   Compile();
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.
///
/// \param[in] s          string to compile into regular expression
/// \param[in] opts       PCRE-style option flags to be set on TPME object
/// \param[in] nMatchMax  maximum number of matches

TPMERegexp::TPMERegexp(const TString& s, UInt_t opts, Int_t nMatchMax) :
   TPRegexp(s),
   fNMaxMatches(nMatchMax),
   fNMatches(0),
   fAddressOfLastString(0),
   fLastGlobalPosition(0)
{
   fPCREOpts = opts;
   Compile();
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.
/// Only PCRE specifics are copied, not last-match or global-match
/// information.

TPMERegexp::TPMERegexp(const TPMERegexp& r) :
   TPRegexp(r),
   fNMaxMatches(r.fNMaxMatches),
   fNMatches(0),
   fAddressOfLastString(0),
   fLastGlobalPosition(0)
{
   Compile();
}

////////////////////////////////////////////////////////////////////////////////
/// Reset the pattern and options.
/// If 'nMatchMax' other than -1 (the default) is passed, it is also set.

void TPMERegexp::Reset(const TString& s, const TString& opts, Int_t nMatchMax)
{
   Reset(s, ParseMods(opts), nMatchMax);
}

////////////////////////////////////////////////////////////////////////////////
/// Reset the pattern and options.
/// If 'nMatchMax' other than -1 (the default) is passed, it is also set.

void TPMERegexp::Reset(const TString& s, UInt_t opts, Int_t nMatchMax)
{
   fPattern = s;
   fPCREOpts = opts;
   Compile();

   if (nMatchMax != -1)
      fNMatches = nMatchMax;
   fNMatches = 0;
   fLastGlobalPosition = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy global-match state from 're; so that this regexp can continue
/// parsing the string from where 're' left off.
///
/// Alternatively, GetGlobalPosition() get be used to retrieve the
/// last match position so that it can passed to Match().
///
/// Ideally, as it is done in PERL, the last match position would be
/// stored in the TString itself.

void TPMERegexp::AssignGlobalState(const TPMERegexp& re)
{
   fLastStringMatched  = re.fLastStringMatched;
   fLastGlobalPosition = re.fLastGlobalPosition;
}

////////////////////////////////////////////////////////////////////////////////
/// Reset state of global match.
/// This happens automatically when a new string is passed for matching.
/// But be carefull, as the address of last TString object is used
/// to make this decision.

void TPMERegexp::ResetGlobalState()
{
   fLastGlobalPosition = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Runs a match on s against the regex 'this' was created with.
///
/// \param[in] s         string to match against
/// \param[in] start     offset at which to start matching
/// \return              number of matches found

Int_t TPMERegexp::Match(const TString& s, UInt_t start)
{
   // If we got a new string, reset the global position counter.
   if (fAddressOfLastString != (void*) &s) {
      fLastGlobalPosition = 0;
   }

   if (fPCREOpts & kPCRE_GLOBAL) {
      start += fLastGlobalPosition;
   }

   //fprintf(stderr, "string: '%s' length: %d offset: %d\n", s.Data(), s.length(), offset);
   fNMatches = MatchInternal(s, start, fNMaxMatches, &fMarkers);

   //fprintf(stderr, "MatchInternal_exec result = %d\n", fNMatches);

   fLastStringMatched   = s;
   fAddressOfLastString = (void*) &s;

   if (fPCREOpts & kPCRE_GLOBAL) {
      if (fNMatches == PCRE_ERROR_NOMATCH) {
         // fprintf(stderr, "TPME RESETTING: reset for no match\n");
         fLastGlobalPosition = 0; // reset the position for next match (perl does this)
      } else if (fNMatches > 0) {
         // fprintf(stderr, "TPME RESETTING: setting to %d\n", marks[0].second);
         fLastGlobalPosition = fMarkers[1]; // set to the end of the match
      } else {
         // fprintf(stderr, "TPME RESETTING: reset for no unknown\n");
         fLastGlobalPosition = 0;
      }
   }

   return fNMatches;
}

////////////////////////////////////////////////////////////////////////////////
/// Splits into at most maxfields. If maxfields is unspecified or
/// 0, trailing empty matches are discarded. If maxfields is
/// positive, no more than maxfields fields will be returned and
/// trailing empty matches are preserved. If maxfields is empty,
/// all fields (including trailing empty ones) are returned. This
/// *should* be the same as the perl behaviour.
///
/// If pattern produces sub-matches, these are also stored in
/// the result.
///
/// A pattern matching the null string will split the value of EXPR
/// into separate characters at each point it matches that way.
///
/// \param[in] s           string to split
/// \param[in] maxfields   maximum number of fields to be split out.  0 means
///                        split all fields, but discard any trailing empty bits.
///                        Negative means split all fields and keep trailing empty bits.
///                        Positive means keep up to N fields including any empty fields
///                        less than N. Anything remaining is in the last field.
/// \return                number of fields found

Int_t TPMERegexp::Split(const TString& s, Int_t maxfields)
{
   typedef std::pair<int, int>   MarkerLoc_t;
   typedef std::vector<MarkerLoc_t> MarkerLocVec_t;

   // stores the marks for the split
   MarkerLocVec_t oMarks;

   // this is a list of current trailing empty matches if maxfields is
   //   unspecified or 0.  If there is stuff in it and a non-empty match
   //   is found, then everything in here is pushed into oMarks and then
   //   the new match is pushed on.  If the end of the string is reached
   //   and there are empty matches in here, they are discarded.
   MarkerLocVec_t oCurrentTrailingEmpties;

   Int_t nOffset = 0;
   Int_t nMatchesFound = 0;

   // while we are still finding matches and maxfields is 0 or negative
   //   (meaning we get all matches), or we haven't gotten to the number
   //   of specified matches
   Int_t matchRes;
   while ((matchRes = Match(s, nOffset)) &&
          ((maxfields < 1) || nMatchesFound < maxfields)) {
      ++nMatchesFound;

      if (fMarkers[1] - fMarkers[0] == 0) {
         oMarks.push_back(MarkerLoc_t(nOffset, nOffset + 1));
         ++nOffset;
         if (nOffset >= s.Length())
            break;
         else
            continue;
      }

      // match can be empty
      if (nOffset != fMarkers[0]) {
         if (!oCurrentTrailingEmpties.empty()) {
            oMarks.insert(oMarks.end(),
                          oCurrentTrailingEmpties.begin(),
                          oCurrentTrailingEmpties.end());
            oCurrentTrailingEmpties.clear();
         }
         oMarks.push_back(MarkerLoc_t(nOffset, fMarkers[0]));
      } else {
         // empty match
         if (maxfields == 0) {
            // store for possible later inclusion
            oCurrentTrailingEmpties.push_back(MarkerLoc_t(nOffset, nOffset));
         } else {
            oMarks.push_back(MarkerLoc_t(nOffset, nOffset));
         }
      }

      nOffset = fMarkers[1];

      if (matchRes > 1) {
         for (Int_t i = 1; i < matchRes; ++i)
            oMarks.push_back(MarkerLoc_t(fMarkers[2*i], fMarkers[2*i + 1]));
      }
   }


   // if there were no matches found, push the whole thing on
   if (nMatchesFound == 0) {
      oMarks.push_back(MarkerLoc_t(0, s.Length()));
   }
   // if we ran out of matches, then append the rest of the string
   //   onto the end of the last split field
   else if (maxfields > 0 && nMatchesFound >= maxfields) {
      oMarks[oMarks.size() - 1].second = s.Length();
   }
   // else we have to add another entry for the end of the string
   else {
      Bool_t last_empty = (nOffset == s.Length());
      if (!last_empty || maxfields < 0) {
         if (!oCurrentTrailingEmpties.empty()) {
            oMarks.insert(oMarks.end(),
                          oCurrentTrailingEmpties.begin(),
                          oCurrentTrailingEmpties.end());
         }
         oMarks.push_back(MarkerLoc_t(nOffset, s.Length()));
      }
   }

   fNMatches = oMarks.size();
   fMarkers.Set(2*fNMatches);
   for (Int_t i = 0; i < fNMatches; ++i) {
      fMarkers[2*i]     = oMarks[i].first;
      fMarkers[2*i + 1] = oMarks[i].second;
   }

   // fprintf(stderr, "match returning %d\n", fNMatches);
   return fNMatches;
}

////////////////////////////////////////////////////////////////////////////////
/// Substitute matching part of s with r, dollar back-ref
/// substitution is performed if doDollarSubst is true (default).
/// Returns the number of substitutions made.
///
/// After the substitution, another pass is made over the resulting
/// string and the following special tokens are interpreted:
///  - `\l`   lowercase next char,
///  - `\u`   uppercase next char,
///  - `\L`   lowercase till `\E`,
///  - `\U`   uppercase till `\E`, and
///  - `\E`   end case modification.

Int_t TPMERegexp::Substitute(TString& s, const TString& r, Bool_t doDollarSubst)
{
   Int_t cnt = SubstituteInternal(s, r, 0, fNMaxMatches, doDollarSubst);

   TString ret;
   Int_t   state = 0;
   Ssiz_t  pos = 0, len = s.Length();
   const Char_t *data = s.Data();
   while (pos < len) {
      Char_t c = data[pos];
      if (c == '\\') {
         c = data[pos+1]; // Rely on string-data being null-terminated.
         switch (c) {
            case  0 : ret += '\\'; break;
            case 'l': state = 1;   break;
            case 'u': state = 2;   break;
            case 'L': state = 3;   break;
            case 'U': state = 4;   break;
            case 'E': state = 0;   break;
            default : ret += '\\'; ret += c; break;
         }
         pos += 2;
      } else {
         switch (state) {
            case 0:  ret += c; break;
            case 1:  ret += (Char_t) tolower(c); state = 0; break;
            case 2:  ret += (Char_t) toupper(c); state = 0; break;
            case 3:  ret += (Char_t) tolower(c); break;
            case 4:  ret += (Char_t) toupper(c); break;
            default: Error("TPMERegexp::Substitute", "invalid state.");
         }
         ++pos;
      }
   }

   s = ret;

   return cnt;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the sub-string from the internal fMarkers vector.
/// Requires having run match or split first.

TString TPMERegexp::operator[](int index)
{
   if (index >= fNMatches)
      return "";

   Int_t begin = fMarkers[2*index];
   Int_t end   = fMarkers[2*index + 1];
   return fLastStringMatched(begin, end-begin);
}

////////////////////////////////////////////////////////////////////////////////
/// Print the regular expression and modifier options.
/// If 'option' contains "all", prints also last string match and
/// match results.

void TPMERegexp::Print(Option_t* option)
{
   TString opt = option;
   opt.ToLower();

   Printf("Regexp='%s', Opts='%s'", fPattern.Data(), GetModifiers().Data());
   if (opt.Contains("all")) {
      Printf("  last string='%s'", fLastStringMatched.Data());
      Printf("  number of matches = %d", fNMatches);
      for (Int_t i=0; i<fNMatches; ++i)
         Printf("  %d - %s", i, operator[](i).Data());
   }
}


/** \class TStringToken
Provides iteration through tokens of a given string.

  - fFullStr     stores the string to be split. It is never modified.
  - fSplitRe     is the perl-re that is used to separate the tokens.
  - fReturnVoid  if true, empty strings will be returned.

Current token is stored in the TString base-class.
During construction no match is done, use NextToken() to get the first
and all subsequent tokens.
*/

ClassImp(TStringToken);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TStringToken::TStringToken(const TString& fullStr, const TString& splitRe, Bool_t retVoid) :
   fFullStr    (fullStr),
   fSplitRe    (splitRe),
   fReturnVoid (retVoid),
   fPos        (0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Get the next token, it is stored in this TString.
/// Returns true if new token is available, false otherwise.

Bool_t TStringToken::NextToken()
{
   TArrayI x;
   while (fPos < fFullStr.Length()) {
      if (fSplitRe.Match(fFullStr, "", fPos, 2, &x)) {
         TString::operator=(fFullStr(fPos, x[0] - fPos));
         fPos = x[1];
      } else {
         TString::operator=(fFullStr(fPos, fFullStr.Length() - fPos));
         fPos = fFullStr.Length() + 1;
      }
      if (Length() || fReturnVoid)
         return kTRUE;
   }

   // Special case: void-strings are requested and the full-string
   // ends with the separator. Thus we return another empty string.
   if (fPos == fFullStr.Length() && fReturnVoid) {
      TString::operator=("");
      fPos = fFullStr.Length() + 1;
      return kTRUE;
   }

   return kFALSE;
}
