// @(#)root/base:$Id$
// Author: Eddy Offermann   24/06/05

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TPRegexp
#define ROOT_TPRegexp

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TPRegexp                                                             //
//                                                                      //
// C++ Wrapper for the "Perl Compatible Regular Expressions" library    //
//  The PCRE lib can be found at:                                       //
//              http://www.pcre.org/                                    //
//                                                                      //
// Extensive documentation about Regular expressions in Perl can be     //
// found at :                                                           //
//              http://perldoc.perl.org/perlre.html                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "Rtypes.h"
#include "TString.h"
#include "TArrayI.h"

struct PCREPriv_t;


class TPRegexp {

protected:
   enum {
      kPCRE_GLOBAL     = 0x80000000,
      kPCRE_OPTIMIZE   = 0x40000000,
      kPCRE_DEBUG_MSGS = 0x20000000,
      kPCRE_INTMASK    = 0x0FFF
   };

   TString     fPattern;
   PCREPriv_t *fPriv;
   UInt_t      fPCREOpts;

   static Bool_t fgThrowAtCompileError;

   void     Compile();
   void     Optimize();
   UInt_t   ParseMods(const TString &mods) const;
   Int_t    ReplaceSubs(const TString &s, TString &final,
                        const TString &replacePattern,
                        Int_t *ovec, Int_t nmatch) const;

   Int_t    MatchInternal(const TString& s, Int_t start,
                          Int_t nMaxMatch, TArrayI *pos=0) const;

   Int_t    SubstituteInternal(TString &s, const TString &replace,
                               Int_t start, Int_t nMaxMatch0,
                               Bool_t doDollarSubst) const;

public:
   TPRegexp();
   TPRegexp(const TString &pat);
   TPRegexp(const TPRegexp &p);
   virtual ~TPRegexp();

   Bool_t     IsValid() const;

   Int_t      Match(const TString &s, const TString &mods="",
                    Int_t start=0, Int_t nMaxMatch=10, TArrayI *pos=0);
   TObjArray *MatchS(const TString &s, const TString &mods="",
                     Int_t start=0, Int_t nMaxMatch=10);
   Bool_t     MatchB(const TString &s, const TString &mods="",
                     Int_t start=0, Int_t nMaxMatch=10) {
                           return (Match(s,mods,start,nMaxMatch) > 0); }
   Int_t      Substitute(TString &s, const TString &replace,
                         const TString &mods="", Int_t start=0,
                         Int_t nMatchMax=10);

   TString GetPattern()   const { return fPattern; }
   TString GetModifiers() const;

   TPRegexp &operator=(const TPRegexp &p);

   static Bool_t GetThrowAtCompileError();
   static void   SetThrowAtCompileError(Bool_t throwp);

   ClassDef(TPRegexp,0)  // Perl Compatible Regular Expression Class
};


class TPMERegexp : protected TPRegexp {

private:
   TPMERegexp& operator=(const TPMERegexp&) = delete;

protected:
   Int_t    fNMaxMatches;         // maximum number of matches
   Int_t    fNMatches;            // number of matches returned from last pcre_exec call
   TArrayI  fMarkers;             // last set of indexes of matches

   TString  fLastStringMatched;   // copy of the last TString matched
   void    *fAddressOfLastString; // used for checking for change of TString in global match

   Int_t    fLastGlobalPosition;  // end of last match when kPCRE_GLOBAL is set

public:
   TPMERegexp();
   TPMERegexp(const TString& s, const TString& opts = "", Int_t nMatchMax = 10);
   TPMERegexp(const TString& s, UInt_t opts, Int_t nMatchMax = 10);
   TPMERegexp(const TPMERegexp& r);

   virtual ~TPMERegexp() {}

   void    Reset(const TString& s, const TString& opts = "", Int_t nMatchMax = -1);
   void    Reset(const TString& s, UInt_t opts, Int_t nMatchMax = -1);

   Int_t   GetNMaxMatches()   const { return fNMaxMatches; }
   void    SetNMaxMatches(Int_t nm) { fNMaxMatches = nm; }

   Int_t   GetGlobalPosition() const { return fLastGlobalPosition; }
   void    AssignGlobalState(const TPMERegexp& re);
   void    ResetGlobalState();

   Int_t   Match(const TString& s, UInt_t start = 0);
   Int_t   Split(const TString& s, Int_t maxfields = 0);
   Int_t   Substitute(TString& s, const TString& r, Bool_t doDollarSubst=kTRUE);

   Int_t   NMatches() const { return fNMatches; }
   TString operator[](Int_t);

   virtual void Print(Option_t* option="");

   ClassDef(TPMERegexp, 0); // Wrapper for Perl-like regular expression matching.
};


class TStringToken : public TString {

protected:
   const TString fFullStr;
   TPRegexp      fSplitRe;
   Bool_t        fReturnVoid;
   Int_t         fPos;

public:
   TStringToken(const TString& fullStr, const TString& splitRe, Bool_t retVoid=kFALSE);
   virtual ~TStringToken() {}

   Bool_t NextToken();
   Bool_t AtEnd() const { return fPos >= fFullStr.Length(); }

   ClassDef(TStringToken,0) // String tokenizer using PCRE for finding next tokens.
};

#endif
