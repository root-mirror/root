// @(#)root/eve:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TEveUtil.h"
#include "TEveElement.h"
#include "TEveManager.h"

#include "TError.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TClass.h"
#include "TMath.h"
#include "TColor.h"

#include "TROOT.h"
#include "TInterpreter.h"
#include "TSystem.h"

#include "TGClient.h"
#include "TGMimeTypes.h"

#include <list>
#include <algorithm>
#include <string>

/** \class TEveUtil
\ingroup TEve
Standard utility functions for Eve.
*/

ClassImp(TEveUtil);

TObjArray* TEveUtil::fgDefaultColors = 0;

////////////////////////////////////////////////////////////////////////////////
/// Setup Include and Macro paths.
/// Since inclusion into ROOT this does nothing but could
/// potentially be reactivated if some common macros are established
/// and shipped with binary ROOT (in macros/eve). For example, these
/// might be used to spawn specific GUI / GL configurations.

void TEveUtil::SetupEnvironment()
{
   static const TEveException eh("TEveUtil::SetupEnvironment");
   static Bool_t setupDone = kFALSE;

   if (setupDone) {
      Info(eh.Data(), "has already been run.");
      return;
   }

   // Old initialization for ALICE.
   // Left as an example.
   /*
   // Check if REVESYS exists, try fallback to $ALICE_ROOT/EVE.
   if(gSystem->Getenv("REVESYS") == 0) {
      if(gSystem->Getenv("ALICE_ROOT") != 0) {
         Info(eh.Data(), "setting REVESYS from ALICE_ROOT.");
         gSystem->Setenv("REVESYS", Form("%s/EVE", gSystem->Getenv("ALICE_ROOT")));
      } else {
         Error(eh.Data(), "REVESYS not defined, neither is ALICE_ROOT.");
         gSystem->Exit(1);
      }
   }
   if(gSystem->AccessPathName(gSystem->Getenv("REVESYS")) == kTRUE) {
      Error(eh.Data(), "REVESYS '%s' does not exist.", gSystem->Getenv("REVESYS"));
      gSystem->Exit(1);
   }

   TString macPath(gROOT->GetMacroPath());
   macPath += Form(":%s/macros", gSystem->Getenv("REVESYS"));
   gInterpreter->AddIncludePath(gSystem->Getenv("REVESYS"));
   if(gSystem->Getenv("ALICE_ROOT") != 0) {
      macPath += Form(":%s/alice-macros", gSystem->Getenv("REVESYS"));
      gInterpreter->AddIncludePath(Form("%s/include", gSystem->Getenv("ALICE_ROOT")));
      gInterpreter->AddIncludePath(gSystem->Getenv("ALICE_ROOT"));
   }
   gROOT->SetMacroPath(macPath);
   */
}

////////////////////////////////////////////////////////////////////////////////
/// Setup icon pictures and mime-types.

void TEveUtil::SetupGUI()
{
   TEveElement::fgRnrIcons[0] = gClient->GetPicture("eve_rnr00_t.xpm");
   TEveElement::fgRnrIcons[1] = gClient->GetPicture("eve_rnr01_t.xpm");
   TEveElement::fgRnrIcons[2] = gClient->GetPicture("eve_rnr10_t.xpm");
   TEveElement::fgRnrIcons[3] = gClient->GetPicture("eve_rnr11_t.xpm");

   TEveElement::fgListTreeIcons[0] = gClient->GetPicture("folder_t.xpm");
   TEveElement::fgListTreeIcons[1] = gClient->GetPicture("eve_viewer.xpm");
   TEveElement::fgListTreeIcons[2] = gClient->GetPicture("eve_scene.xpm");
   TEveElement::fgListTreeIcons[3] = gClient->GetPicture("eve_pointset.xpm");
   TEveElement::fgListTreeIcons[4] = gClient->GetPicture("eve_track.xpm");
   TEveElement::fgListTreeIcons[5] = gClient->GetPicture("eve_text.gif");
   TEveElement::fgListTreeIcons[6] = gClient->GetPicture("eve_axes.xpm");
   TEveElement::fgListTreeIcons[7] = gClient->GetPicture("ofolder_t.xpm");
   TEveElement::fgListTreeIcons[8] = gClient->GetPicture("eve_line.xpm");

   gClient->GetMimeTypeList()->AddType("root/tmacro", "TEveMacro",
                                       "tmacro_s.xpm", "tmacro_t.xpm", "");
}

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

Bool_t TEveUtil::CheckMacro(const char* mac)
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

void TEveUtil::AssertMacro(const char* mac)
{
   if( CheckMacro(mac) == kFALSE) {
      gROOT->Macro(mac);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Execute macro 'mac'. Do not reload the macro.

void TEveUtil::Macro(const char* mac)
{
   if (CheckMacro(mac) == kFALSE) {
      gROOT->LoadMacro(mac);
   }
   TString foo(mac); ChompTailAndDir(foo); foo += "()";
   gROOT->ProcessLine(foo.Data());
}

////////////////////////////////////////////////////////////////////////////////
/// Makes sure that macro 'mac' is loaded, but do not reload it.

void TEveUtil::LoadMacro(const char* mac)
{
   if (CheckMacro(mac) == kFALSE) {
      gROOT->LoadMacro(mac);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fill col with RGBA values corresponding to index ci. If alpha
/// is true, set alpha component of col to 255.
/// ROOT's indexed color palette does not support transparency.

void TEveUtil::ColorFromIdx(Color_t ci, UChar_t col[4], Bool_t alpha)
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

void TEveUtil::ColorFromIdx(Color_t ci, UChar_t col[4], Char_t transparency)
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

void TEveUtil::ColorFromIdx(Float_t f1, Color_t c1, Float_t f2, Color_t c2,
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
/// TEveElementObjectPtr instances.

Color_t* TEveUtil::FindColorVar(TObject* obj, const char* varname)
{
   static const TEveException eh("TEveUtil::FindColorVar");

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
/// You should call TEveManager::FullRedraw3D() afterwards or set
/// the argument full_redraw to true (default is false).

void TEveUtil::SetColorBrightness(Float_t value, Bool_t full_redraw)
{
   if (value < -2.5 || value > 2.5)
   {
      Error("TEveUtil::SetColorBrightness", "value '%f' out of range [-0.5, 0.5].", value);
      return;
   }

   TObjArray   *colors = (TObjArray*) gROOT->GetListOfColors();

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
         TColor* croot = (TColor*)  colors->At(i);
         if (croot == 0)
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

   if (full_redraw && gEve != 0)
      gEve->FullRedraw3D();
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if interval Q is contained within interval M for U1 variables.
/// It is assumed that all values are within the [-2pi, 2pi] interval and
/// minM <= maxM & minQ <= maxQ.

Bool_t TEveUtil::IsU1IntervalContainedByMinMax(Float_t minM, Float_t maxM,
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

Bool_t TEveUtil::IsU1IntervalOverlappingByMinMax(Float_t minM, Float_t maxM,
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

Float_t TEveUtil::GetFraction(Float_t minM, Float_t maxM, Float_t minQ, Float_t maxQ)
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


/** \class TEveException
\ingroup TEve
Exception class thrown by TEve classes and macros.
*/

ClassImp(TEveException);

////////////////////////////////////////////////////////////////////////////////

bool operator==(const TString& t, const std::string& s)
{ return (s == t.Data()); }

bool operator==(const std::string&  s, const TString& t)
{ return (s == t.Data()); }

// Exc

TEveException::TEveException(const std::string& s) : TString(s.c_str())
{
   // Constructor.
}

// Exc + ops

TEveException operator+(const TEveException &s1, const std::string &s2)
{ TEveException r(s1); r += s2; return r; }

TEveException operator+(const TEveException &s1, const TString &s2)
{ TEveException r(s1); r += s2; return r; }

TEveException operator+(const TEveException &s1,  const char *s2)
{ TEveException r(s1); r += s2; return r; }


/** \class TEvePadHolder
\ingroup TEve
Exception safe wrapper for setting gPad.
Optionally calls gPad->Modified()/Update() in destructor.
*/

ClassImp(TEvePadHolder);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEvePadHolder::TEvePadHolder(Bool_t modify_update_p, TVirtualPad* new_pad, Int_t subpad) :
   fOldPad        (gPad),
   fModifyUpdateP (modify_update_p)
{
   if (new_pad != 0)
      new_pad->cd(subpad);
   else
      gPad = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TEvePadHolder::~TEvePadHolder()
{
   if (fModifyUpdateP && gPad != 0) {
      gPad->Modified();
      gPad->Update();
   }
   gPad = fOldPad;
}

/** \class TEveGeoManagerHolder
\ingroup TEve
Exception safe wrapper for setting gGeoManager.
Functionality to lock-unlock via setting of a static lock in
TGeoManager should be added (new feature of TGeoManager).
*/

ClassImp(TEveGeoManagerHolder);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.
/// If n_seg is specified and larger than 2, the new geo-manager's
/// NSegments is set to this value.

TEveGeoManagerHolder::TEveGeoManagerHolder(TGeoManager* new_gmgr, Int_t n_seg) :
   fManager   (gGeoManager),
   fNSegments (0)
{
   gGeoManager = new_gmgr;
   if (gGeoManager)
   {
      gGeoIdentity = (TGeoIdentity*) gGeoManager->GetListOfMatrices()->At(0);
      if (n_seg > 2)
      {
         fNSegments = gGeoManager->GetNsegments();
         gGeoManager->SetNsegments(n_seg);
      }
   }
   else
   {
      gGeoIdentity = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TEveGeoManagerHolder::~TEveGeoManagerHolder()
{
   if (gGeoManager && fNSegments > 2)
   {
      gGeoManager->SetNsegments(fNSegments);
   }
   gGeoManager = fManager;
   if (gGeoManager)
   {
      gGeoIdentity = (TGeoIdentity*) gGeoManager->GetListOfMatrices()->At(0);
   }
   else
   {
      gGeoIdentity = 0;
   }
}

/** \class TEveRefCnt
\ingroup TEve
Base-class for reference-counted objects.
By default the object is destroyed when zero reference-count is reached.
*/

ClassImp(TEveRefCnt);

/** \class TEveRefBackPtr
\ingroup TEve
Base-class for reference-counted objects with reverse references to
TEveElement objects.
*/

ClassImp(TEveRefBackPtr);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

TEveRefBackPtr::TEveRefBackPtr() :
   TEveRefCnt(),
   fBackRefs()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor. Noop, should complain if back-ref list is not empty.

TEveRefBackPtr::~TEveRefBackPtr()
{
   // !!! Complain if list not empty.
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor. New copy starts with zero reference count and
/// empty back-reference list.

TEveRefBackPtr::TEveRefBackPtr(const TEveRefBackPtr&) :
   TEveRefCnt(),
   fBackRefs()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Assignment operator. Reference count and back-reference
/// information is not assigned as these object hold pointers to a
/// specific object.

TEveRefBackPtr& TEveRefBackPtr::operator=(const TEveRefBackPtr&)
{
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Increase reference count and add re to the list of back-references.

void TEveRefBackPtr::IncRefCount(TEveElement* re)
{
   TEveRefCnt::IncRefCount();
   ++fBackRefs[re];
}

////////////////////////////////////////////////////////////////////////////////
/// Decrease reference count and remove re from the list of back-references.

void TEveRefBackPtr::DecRefCount(TEveElement* re)
{
   static const TEveException eh("TEveRefBackPtr::DecRefCount ");

   RefMap_i i = fBackRefs.find(re);
   if (i != fBackRefs.end()) {
      if (--(i->second) <= 0)
         fBackRefs.erase(i);
      TEveRefCnt::DecRefCount();
   } else {
      Warning(eh, "render element '%s' not found in back-refs.",
                  re->GetObject(eh)->GetName());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add given stamps to elements in the list of reverse references.

void TEveRefBackPtr::StampBackPtrElements(UChar_t stamps)
{
   RefMap_i i = fBackRefs.begin();
   while (i != fBackRefs.end())
   {
      i->first->AddStamp(stamps);
      ++i;
   }
}
