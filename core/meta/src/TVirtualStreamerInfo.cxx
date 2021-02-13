// @(#)root/base:$Id$
// Author: Rene Brun   05/02/2007
/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TVirtualStreamerInfo
Abstract Interface class describing Streamer information for one class.
*/

#include "TROOT.h"
#include "TSystem.h"
#include "TClass.h"
#include "TVirtualMutex.h"
#include "TInterpreter.h"
#include "TVirtualStreamerInfo.h"
#include "TPluginManager.h"
#include "TStreamerElement.h"
#include "TError.h"
#include "TObjArray.h"

TVirtualStreamerInfo *TVirtualStreamerInfo::fgInfoFactory    = 0;

Bool_t  TVirtualStreamerInfo::fgCanDelete        = kTRUE;
Bool_t  TVirtualStreamerInfo::fgOptimize         = kTRUE;
Bool_t  TVirtualStreamerInfo::fgStreamMemberWise = kTRUE;

ClassImp(TVirtualStreamerInfo);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

TVirtualStreamerInfo::TVirtualStreamerInfo() : fOptimized(kFALSE), fIsBuilt(kFALSE), fIsCompiled(kFALSE)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

TVirtualStreamerInfo::TVirtualStreamerInfo(TClass *cl)
   : TNamed(cl->GetName(),""), fOptimized(kFALSE), fIsBuilt(kFALSE), fIsCompiled(kFALSE)
{
}

////////////////////////////////////////////////////////////////////////////////
///copy constructor

TVirtualStreamerInfo::TVirtualStreamerInfo(const TVirtualStreamerInfo& info)
  : TNamed(info), fOptimized(kFALSE), fIsBuilt(kFALSE), fIsCompiled(kFALSE)
{
}

////////////////////////////////////////////////////////////////////////////////
///assignment operator

TVirtualStreamerInfo& TVirtualStreamerInfo::operator=(const TVirtualStreamerInfo& info)
{
   if(this!=&info) {
      TNamed::operator=(info);
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TVirtualStreamerInfo::~TVirtualStreamerInfo()
{
}

////////////////////////////////////////////////////////////////////////////////
/// static function returning true if ReadBuffer can delete object

Bool_t TVirtualStreamerInfo::CanDelete()
{
   return fgCanDelete;
}

////////////////////////////////////////////////////////////////////////////////
/// static function returning true if optimization can be on

Bool_t TVirtualStreamerInfo::CanOptimize()
{
   return fgOptimize;
}

////////////////////////////////////////////////////////////////////////////////
/// Given a comment/title declaring an array counter, for example:
/// ~~~ {.cpp}
///    //[fArraySize] array of size fArraySize
/// ~~~
/// return the start of the array dimension declaration start in the string
/// (so the location of the 'f'.

const char *TVirtualStreamerInfo::GetElementCounterStart(const char *dmTitle)
{
   for (const char *lbracket = dmTitle; *lbracket; ++lbracket) {
      // = ::strchr(dmTitle, '[');
      if ( (*lbracket) == '[' ) return lbracket;
      if ( (*lbracket) != '/' && !isspace(*lbracket) ) {
         // Allow only comment delimiters and white spaces
         // before the array information.
         return 0;
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get pointer to a TStreamerBasicType in TClass *cl
///static function

TStreamerBasicType *TVirtualStreamerInfo::GetElementCounter(const char *countName, TClass *cl)
{
   TVirtualStreamerInfo *info;
   {
      R__LOCKGUARD(gInterpreterMutex);
      const TObjArray *sinfos = cl->GetStreamerInfos();
      info = (TVirtualStreamerInfo *)sinfos->At(cl->GetClassVersion());
   }

   if (!info || !info->IsCompiled()) {
      // Even if the streamerInfo exist, it could still need to be 'build'
      // It is important to figure this out, because
      //   a) if it is not build, we need to build
      //   b) if is build, we should not build it (or we could end up in an
      //      infinite loop, if the element and its counter are in the same
      //      class!
      // Checking IsCompiled is sufficient here even-though it is set only at
      // the end of the call to Build as this function has an
      // internal recursion prevention (setting and testing kBuildRunning).
      info = cl->GetStreamerInfo();
   }
   if (!info) return 0;
   TStreamerElement *element = (TStreamerElement *)info->GetElements()->FindObject(countName);
   if (!element) return 0;
   if (element->IsA() == TStreamerBasicType::Class()) return (TStreamerBasicType*)element;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return whether the TStreamerInfos will save the collections in
/// "member-wise" order whenever possible.    The default is to store member-wise.
///  - kTRUE indicates member-wise storing
///  - kFALSE inddicates object-wise storing
///
/// A collection can be saved member wise when it contain is guaranteed to be
/// homogeneous.  For example std::vector<THit> can be stored member wise,
/// while std::vector<THit*> can not (possible use of polymorphism).

Bool_t TVirtualStreamerInfo::GetStreamMemberWise()
{
   return fgStreamMemberWise;
}

////////////////////////////////////////////////////////////////////////////////
///  This is a static function.
///  Set optimization option.
///  When this option is activated (default), consecutive data members
///  of the same type are merged into an array (faster).
///  Optimization must be off in TTree split mode.

void TVirtualStreamerInfo::Optimize(Bool_t opt)
{
   fgOptimize = opt;
}

////////////////////////////////////////////////////////////////////////////////
/// Static function returning a pointer to a new TVirtualStreamerInfo object.
/// If the Info factory does not exist, it is created via the plugin manager.
/// In reality the factory is an empty TStreamerInfo object.

TVirtualStreamerInfo *TVirtualStreamerInfo::Factory()
{
   if (!fgInfoFactory) {
      R__LOCKGUARD(gInterpreterMutex);
      TPluginHandler *h;
      if ((h = gROOT->GetPluginManager()->FindHandler("TVirtualStreamerInfo","TStreamerInfo"))) {
         if (h->LoadPlugin() == -1) {
            ::Fatal("TVirtualStreamerInfo::Factory",
                    "The plugin handler for TVirtualStreamerInfo was found but failed to load!");
         }
         fgInfoFactory = (TVirtualStreamerInfo*) h->ExecPlugin(0);
         if (fgInfoFactory == 0) {
            ::Fatal("TVirtualStreamerInfo::Factory",
                    "The plugin handler for TVirtualStreamerInfo was found but failed to create the factory object!");
         }
      } else {
         TString filename("$ROOTSYS/etc/plugins/TVirtualStreamerInfo");
         gSystem->ExpandPathName(filename);
         if (gSystem->AccessPathName(filename)) {
            ::Fatal("TVirtualStreamerInfo::Factory",
                    "Cannot find the plugin handler for TVirtualStreamerInfo! "
                    "$ROOTSYS/etc/plugins/TVirtualStreamerInfo does not exist "
                    "or is inaccessible.");
         } else {
            ::Fatal("TVirtualStreamerInfo::Factory",
                    "Cannot find the plugin handler for TVirtualStreamerInfo! "
                    "However $ROOTSYS/etc/plugins/TVirtualStreamerInfo is accessible, "
                    "Check the content of this directory!");
         }
      }
   }

   return fgInfoFactory;
}

////////////////////////////////////////////////////////////////////////////////
///  This is a static function.
///  Set object delete option.
///  When this option is activated (default), ReadBuffer automatically
///  delete objects when a data member is a pointer to an object.
///  If your constructor is not presetting pointers to 0, you must
///  call this static function TStreamerInfo::SetCanDelete(kFALSE);

void TVirtualStreamerInfo::SetCanDelete(Bool_t opt)
{
   fgCanDelete = opt;
}

////////////////////////////////////////////////////////////////////////////////
///static function: Set the StreamerInfo factory

void TVirtualStreamerInfo::SetFactory(TVirtualStreamerInfo *factory)
{
   R__LOCKGUARD(gInterpreterMutex);
   auto old = fgInfoFactory;
   fgInfoFactory = factory;
   if (old) delete old;
}

////////////////////////////////////////////////////////////////////////////////
/// Set whether the TStreamerInfos will save the collections in
/// "member-wise" order whenever possible.  The default is to store member-wise.
///  - kTRUE indicates member-wise storing
///  - kFALSE indicates object-wise storing
/// This function returns the previous value of fgStreamMemberWise.

Bool_t TVirtualStreamerInfo::SetStreamMemberWise(Bool_t enable)
{
   // A collection can be saved member wise when it contain is guaranteed to be
   // homogeneous.  For example std::vector<THit> can be stored member wise,
   // while std::vector<THit*> can not (possible use of polymorphism).

   Bool_t prev = fgStreamMemberWise;
   fgStreamMemberWise = enable;
   return prev;
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TVirtualStreamerInfo.

void TVirtualStreamerInfo::Streamer(TBuffer &R__b)
{
   TNamed::Streamer(R__b);
}
