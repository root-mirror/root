// @(#)root/eve:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TEveUtil
#define ROOT_TEveUtil

#include "TObject.h"
#include "TString.h"
#include "TError.h"

#include "GuiTypes.h"

#include <map>
#include <string>
#include <exception>

class TVirtualPad;
class TGeoManager;

class TEveElement;

/******************************************************************************/
// TEveUtil
/******************************************************************************/

class TEveUtil
{
private:
   static TObjArray* fgDefaultColors;

public:
   virtual ~TEveUtil() {}

   // Environment, Macro functions

   static void   SetupEnvironment();
   static void   SetupGUI();

   static Bool_t CheckMacro(const char* mac);
   static void   AssertMacro(const char* mac);
   static void   Macro(const char* mac);
   static void   LoadMacro(const char* mac);

   // Color management

   static void     ColorFromIdx(Color_t ci, UChar_t col[4], Bool_t alpha=kTRUE);
   static void     ColorFromIdx(Color_t ci, UChar_t col[4], Char_t transparency);
   static void     ColorFromIdx(Float_t f1, Color_t c1, Float_t f2, Color_t c2,
                                UChar_t col[4], Bool_t alpha=kTRUE);
   static Color_t* FindColorVar(TObject* obj, const char* varname);

   static void     SetColorBrightness(Float_t value, Bool_t full_redraw=kFALSE);


   // Math utilities

   static Bool_t IsU1IntervalContainedByMinMax  (Float_t minM, Float_t maxM,
                                                 Float_t minQ, Float_t maxQ);
   static Bool_t IsU1IntervalOverlappingByMinMax(Float_t minM, Float_t maxM,
                                                 Float_t minQ, Float_t maxQ);

   static Bool_t IsU1IntervalContainedByMeanDelta  (Float_t meanM, Float_t deltaM,
                                                    Float_t meanQ, Float_t deltaQ);
   static Bool_t IsU1IntervalOverlappingByMeanDelta(Float_t meanM, Float_t deltaM,
                                                    Float_t meanQ, Float_t deltaQ);

   static Float_t GetFraction(Float_t minM, Float_t maxM, Float_t minQ, Float_t maxQ);


   ClassDef(TEveUtil, 0); // Standard utility functions for Reve.
};

inline Bool_t TEveUtil::IsU1IntervalContainedByMeanDelta(Float_t meanM, Float_t deltaM,
                                                         Float_t meanQ, Float_t deltaQ)
{
   return IsU1IntervalContainedByMinMax(meanM - deltaM, meanM + deltaM, meanQ - deltaQ, meanQ + deltaQ);
}

inline Bool_t TEveUtil::IsU1IntervalOverlappingByMeanDelta(Float_t meanM, Float_t deltaM,
                                                             Float_t meanQ, Float_t deltaQ)
{
   return IsU1IntervalContainedByMinMax(meanM - deltaM, meanM + deltaM, meanQ - deltaQ, meanQ + deltaQ);
}


/******************************************************************************/
// Exceptions, string functions
/******************************************************************************/

bool operator==(const TString& t, const std::string& s);
bool operator==(const std::string& s, const TString& t);

class TEveException : public std::exception, public TString
{
public:
   TEveException() {}
   TEveException(const TString& s) : TString(s) {}
   TEveException(const char* s)    : TString(s) {}
   TEveException(const std::string& s);

   virtual ~TEveException() noexcept {}

   virtual const char* what() const noexcept { return Data(); }

   ClassDef(TEveException, 1); // Exception-type thrown by Eve classes.
};

TEveException operator+(const TEveException &s1, const std::string  &s2);
TEveException operator+(const TEveException &s1, const TString &s2);
TEveException operator+(const TEveException &s1, const char    *s2);


/******************************************************************************/
// Exception-safe global variable holders
/******************************************************************************/

class TEvePadHolder
{
private:
   TVirtualPad *fOldPad;
   Bool_t       fModifyUpdateP;

   TEvePadHolder(const TEvePadHolder&);            // Not implemented
   TEvePadHolder& operator=(const TEvePadHolder&); // Not implemented

public:
   TEvePadHolder(Bool_t modify_update_p, TVirtualPad* new_pad=0, Int_t subpad=0);
   virtual ~TEvePadHolder();

   ClassDef(TEvePadHolder, 0); // Exception-safe wrapper for temporary setting of gPad variable.
};

class TEveGeoManagerHolder
{
private:
   TGeoManager *fManager;
   Int_t        fNSegments;

   TEveGeoManagerHolder(const TEveGeoManagerHolder&);            // Not implemented
   TEveGeoManagerHolder& operator=(const TEveGeoManagerHolder&); // Not implemented

public:
   TEveGeoManagerHolder(TGeoManager* new_gmgr=0, Int_t n_seg=0);
   virtual ~TEveGeoManagerHolder();

   ClassDef(TEveGeoManagerHolder, 0); // Exception-safe wrapper for temporary setting of gGeoManager variable.
};


/******************************************************************************/
// TEveRefCnt base-class (interface)
/******************************************************************************/

class TEveRefCnt
{
protected:
   Int_t fRefCount;

public:
   TEveRefCnt() : fRefCount(0) {}
   virtual ~TEveRefCnt() {}

   TEveRefCnt(const TEveRefCnt&) : fRefCount(0) {}
   TEveRefCnt& operator=(const TEveRefCnt&) { return *this; }

   void IncRefCount() { ++fRefCount; }
   void DecRefCount() { if(--fRefCount <= 0) OnZeroRefCount(); }

   virtual void OnZeroRefCount() { delete this; }

   ClassDef(TEveRefCnt, 0); // Base-class for reference-counted objects.
};

/******************************************************************************/
// TEveRefBackPtr reference-count with back pointers
/******************************************************************************/

class TEveRefBackPtr : public TEveRefCnt
{
protected:
   typedef std::map<TEveElement*, Int_t> RefMap_t;
   typedef RefMap_t::iterator            RefMap_i;

   RefMap_t fBackRefs;

public:
   TEveRefBackPtr();
   virtual ~TEveRefBackPtr();

   TEveRefBackPtr(const TEveRefBackPtr&);
   TEveRefBackPtr& operator=(const TEveRefBackPtr&);

   using TEveRefCnt::IncRefCount;
   using TEveRefCnt::DecRefCount;
   virtual void IncRefCount(TEveElement* re);
   virtual void DecRefCount(TEveElement* re);

   virtual void StampBackPtrElements(UChar_t stamps);

   ClassDef(TEveRefBackPtr, 0); // Base-class for reference-counted objects with reverse references to TEveElement objects.
};

#endif
