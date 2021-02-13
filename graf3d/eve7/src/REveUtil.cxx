// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/REveUtil.hxx>
#include <ROOT/REveElement.hxx>
#include <ROOT/REveManager.hxx>

#include <ROOT/RLogger.hxx>

#include "TError.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TClass.h"
#include "TMath.h"

#include "TColor.h"

#include "TROOT.h"
#include "TInterpreter.h"

#include <list>
#include <algorithm>
#include <string>

using namespace ROOT::Experimental;
namespace REX = ROOT::Experimental;

REX::RLogChannel &REX::EveLog() {
   static RLogChannel sLog("ROOT.Eve");
   return sLog;
}

/** \class REveUtil
\ingroup REve
Standard utility functions for Eve.
*/

TObjArray *REveUtil::fgDefaultColors = nullptr;

namespace
{
////////////////////////////////////////////////////////////////////////////////
/// Remove last part of string 's', starting from the last
/// occurrence of character 'c'.
/// Remove directory part -- everything until the last '/'.

void ChompTailAndDir(TString& s, char c='.')
{
   Ssiz_t p = s.Last(c);
   if (p != kNPOS)
      s.Remove(p);

   Ssiz_t ls = s.Last('/');
   if (ls != kNPOS)
      s.Remove(0, ls + 1);
}
}

////////////////////////////////////////////////////////////////////////////////
/// Checks if macro 'mac' is loaded.

Bool_t REveUtil::CheckMacro(const char* mac)
{
   // Axel's advice; now sth seems slow, using old method below for test.
   // return gROOT->GetInterpreter()->IsLoaded(mac);

   // Previous version expected function with same name and used ROOT's
   // list of global functions.

   TString foo(mac); ChompTailAndDir(foo);
   if (gROOT->GetGlobalFunction(foo.Data(), 0, kFALSE) != 0)
      return kTRUE;
   else
      return (gROOT->GetGlobalFunction(foo.Data(), 0, kTRUE) != 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Load and execute macro 'mac' if it has not been loaded yet.

void REveUtil::AssertMacro(const char* mac)
{
   if( CheckMacro(mac) == kFALSE) {
      gROOT->Macro(mac);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Execute macro 'mac'. Do not reload the macro.

void REveUtil::Macro(const char* mac)
{
   if (CheckMacro(mac) == kFALSE) {
      gROOT->LoadMacro(mac);
   }
   TString foo(mac); ChompTailAndDir(foo); foo += "()";
   gROOT->ProcessLine(foo.Data());
}

////////////////////////////////////////////////////////////////////////////////
/// Makes sure that macro 'mac' is loaded, but do not reload it.

void REveUtil::LoadMacro(const char* mac)
{
   if (CheckMacro(mac) == kFALSE) {
      gROOT->LoadMacro(mac);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fill col with RGBA values corresponding to index ci. If alpha
/// is true, set alpha component of col to 255.
/// ROOT's indexed color palette does not support transparency.

void REveUtil::ColorFromIdx(Color_t ci, UChar_t col[4], Bool_t alpha)
{
   TColor* c = gROOT->GetColor(ci);
   if (c)
   {
      col[0] = (UChar_t)(255*c->GetRed());
      col[1] = (UChar_t)(255*c->GetGreen());
      col[2] = (UChar_t)(255*c->GetBlue());
      if (alpha) col[3] = 255;
   }
   else
   {
      // Set to magenta.
      col[0] = 255; col[1] = 0; col[2] = 255;
      if (alpha) col[3] = 255;
      return;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fill col with RGBA values corresponding to index ci and transparency.
/// ROOT's indexed color palette does not support transparency.

void REveUtil::ColorFromIdx(Color_t ci, UChar_t col[4], Char_t transparency)
{
   UChar_t alpha = (255*(100 - transparency))/100;

   TColor* c = gROOT->GetColor(ci);
   if (c)
   {
      col[0] = (UChar_t)(255*c->GetRed());
      col[1] = (UChar_t)(255*c->GetGreen());
      col[2] = (UChar_t)(255*c->GetBlue());
      col[3] = alpha;
   }
   else
   {
      // Set to magenta.
      col[0] = 255; col[1] = 0; col[2] = 255; col[3] = alpha;
      return;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fill col with weighted RGBA values corresponding to
/// color-indices c1 and c2. If alpha is true, set alpha component
/// of col to 255.

void REveUtil::ColorFromIdx(Float_t f1, Color_t c1, Float_t f2, Color_t c2,
                            UChar_t col[4], Bool_t alpha)
{
   TColor* t1 = gROOT->GetColor(c1);
   TColor* t2 = gROOT->GetColor(c2);
   if(t1 && t2) {
      col[0] = (UChar_t)(255*(f1*t1->GetRed()   + f2*t2->GetRed()));
      col[1] = (UChar_t)(255*(f1*t1->GetGreen() + f2*t2->GetGreen()));
      col[2] = (UChar_t)(255*(f1*t1->GetBlue()  + f2*t2->GetBlue()));
      if (alpha) col[3] = 255;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Find address of Color_t data-member with name varname in object
/// obj.
///
/// This is used to access color information for TGListTreeItem
/// coloration from visualization macros that wrap TObjects into
/// REveElementObjectPtr instances.

Color_t* REveUtil::FindColorVar(TObject* obj, const char* varname)
{
   static const REveException eh("REveUtil::FindColorVar");

   Int_t off = obj->IsA()->GetDataMemberOffset(varname);
   if(off == 0)
      throw(eh + "could not find member '" + varname + "' in class " + obj->IsA()->GetName() + ".");
   return (Color_t*) (((char*)obj) + off);
}

////////////////////////////////////////////////////////////////////////////////
/// Tweak all ROOT colors to become brighter (if value > 0) or
/// darker (value < 0). Reasonable values for the value argument are
/// from -2.5 to 2.5 (error will be printed otherwise).
/// If value is zero, the original colors are restored.
///
/// You should call REveManager::FullRedraw3D() afterwards or set
/// the argument full_redraw to true (default is false).

void REveUtil::SetColorBrightness(Float_t value, Bool_t full_redraw)
{
   if (value < -2.5 || value > 2.5)
   {
      Error("REveUtil::SetColorBrightness", "value '%f' out of range [-0.5, 0.5].", value);
      return;
   }

   TObjArray *colors = (TObjArray*) gROOT->GetListOfColors();

   if (fgDefaultColors == 0)
   {
      const Int_t n_col = colors->GetEntriesFast();
      fgDefaultColors = new TObjArray(n_col);
      for (Int_t i = 0; i < n_col; ++i)
      {
         TColor* c = (TColor*) colors->At(i);
         if (c)
            fgDefaultColors->AddAt(new TColor(*c), i);
      }
   }

   const Int_t n_col = fgDefaultColors->GetEntriesFast();
   for (Int_t i = 0; i < n_col; ++i)
   {
      TColor* cdef = (TColor*) fgDefaultColors->At(i);
      if (cdef)
      {
         TColor* croot = (TColor*) colors->At(i);
         if (!croot)
         {
            croot = new TColor(*cdef);
            colors->AddAt(croot, i);
         }
         else
         {
            cdef->Copy(*croot);
         }

         Float_t r, g, b;
         croot->GetRGB(r, g, b);
         r = TMath::Power( r, (2.5 - value)/2.5);
         g = TMath::Power(g, (2.5 - value)/2.5);
         b = TMath::Power(b, (2.5 - value)/2.5);

         r = TMath::Min(r, 1.0f);
         g = TMath::Min(g, 1.0f);
         b = TMath::Min(b, 1.0f);

         croot->SetRGB(r, g, b);
      }
      else
      {
         delete colors->RemoveAt(i);
      }
   }

   if (full_redraw && REX::gEve)
      REX::gEve->FullRedraw3D();
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if interval Q is contained within interval M for U1 variables.
/// It is assumed that all values are within the [-2pi, 2pi] interval and
/// minM <= maxM & minQ <= maxQ.

Bool_t REveUtil::IsU1IntervalContainedByMinMax(Float_t minM, Float_t maxM,
                                               Float_t minQ, Float_t maxQ)
{
   using namespace TMath;

   if (maxQ < minM)
   {
      minQ += TwoPi(); maxQ += TwoPi();
   }
   else if (minQ > maxM)
   {
      minQ -= TwoPi(); maxQ -= TwoPi();
   }
   return minQ >= minM && maxQ <= maxM;
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if interval Q is overlapping within interval M for U1 variables.
/// It is assumed that all values are within the [-2pi, 2pi] interval and
/// minM <= maxM & minQ <= maxQ.

Bool_t REveUtil::IsU1IntervalOverlappingByMinMax(Float_t minM, Float_t maxM,
                                                 Float_t minQ, Float_t maxQ)
{
   using namespace TMath;

   if (maxQ < minM)
   {
      minQ += TwoPi(); maxQ += TwoPi();
   }
   else if (minQ > maxM)
   {
      minQ -= TwoPi(); maxQ -= TwoPi();
   }
   return maxQ >= minM && minQ <= maxM;
}

////////////////////////////////////////////////////////////////////////////////
/// Get fraction of interval [minQ, maxQ] in [minM, maxM]

Float_t REveUtil::GetFraction(Float_t minM, Float_t maxM, Float_t minQ, Float_t maxQ)
{
   if (minQ>=minM && maxQ<=maxM)
      return 1;

   else if (minQ<minM && maxQ>maxM)
      return (maxM-minM)/(maxQ-minQ);

   else if (minQ>=minM && maxQ>maxM)
      return (maxM-minQ)/(maxQ-minQ);

   else if (minQ<minM && maxQ<=maxM)
      return (maxQ-minM)/(maxQ-minQ);

   return 0;
}


/** \class REveGeoManagerHolder
\ingroup REve
Exception safe wrapper for setting gGeoManager.
Functionality to lock-unlock via setting of a static lock in
TGeoManager should be added (new feature of TGeoManager).
*/

////////////////////////////////////////////////////////////////////////////////
/// Constructor.
/// If n_seg is specified and larger than 2, the new geo-manager's
/// NSegments is set to this value.

REveGeoManagerHolder::REveGeoManagerHolder(TGeoManager* new_gmgr, Int_t n_seg) :
   fManager   (gGeoManager),
   fNSegments (0)
{
   gGeoManager = new_gmgr;
   if (gGeoManager) {
      gGeoIdentity = (TGeoIdentity *)gGeoManager->GetListOfMatrices()->At(0);
      if (n_seg > 2) {
         fNSegments = gGeoManager->GetNsegments();
         gGeoManager->SetNsegments(n_seg);
      }
   } else {
      gGeoIdentity = nullptr;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

REveGeoManagerHolder::~REveGeoManagerHolder()
{
   if (gGeoManager && fNSegments > 2) {
      gGeoManager->SetNsegments(fNSegments);
   }
   gGeoManager = fManager;
   if (gGeoManager) {
      gGeoIdentity = (TGeoIdentity *)gGeoManager->GetListOfMatrices()->At(0);
   } else {
      gGeoIdentity = nullptr;
   }
}

/** \class REveRefCnt
\ingroup REve
Base-class for reference-counted objects.
By default the object is destroyed when zero reference-count is reached.
*/

/** \class REveRefBackPtr
\ingroup REve
Base-class for reference-counted objects with reverse references to
REveElement objects.
*/

////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

REveRefBackPtr::REveRefBackPtr() :
   REveRefCnt(),
   fBackRefs()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor. Noop, should complain if back-ref list is not empty.

REveRefBackPtr::~REveRefBackPtr()
{
   // !!! Complain if list not empty.
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor. New copy starts with zero reference count and
/// empty back-reference list.

REveRefBackPtr::REveRefBackPtr(const REveRefBackPtr&) :
   REveRefCnt(),
   fBackRefs()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Assignment operator. Reference count and back-reference
/// information is not assigned as these object hold pointers to a
/// specific object.

REveRefBackPtr& REveRefBackPtr::operator=(const REveRefBackPtr&)
{
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Increase reference count and add re to the list of back-references.

void REveRefBackPtr::IncRefCount(REveElement* re)
{
   REveRefCnt::IncRefCount();
   ++fBackRefs[re];
}

////////////////////////////////////////////////////////////////////////////////
/// Decrease reference count and remove re from the list of back-references.

void REveRefBackPtr::DecRefCount(REveElement *re)
{
   auto i = fBackRefs.find(re);
   if (i != fBackRefs.end()) {
      if (--(i->second) <= 0)
         fBackRefs.erase(i);
      REveRefCnt::DecRefCount();
   } else {
      Warning("REveRefBackPtr::DecRefCount", "element '%s' not found in back-refs.", re->GetCName());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add given stamps to elements in the list of reverse references.

void REveRefBackPtr::StampBackPtrElements(UChar_t stamps)
{
   for (auto &i: fBackRefs)
      i.first->AddStamp(stamps);
}
