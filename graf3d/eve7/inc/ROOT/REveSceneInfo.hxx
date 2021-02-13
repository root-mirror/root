// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_REveSceneInfo
#define ROOT7_REveSceneInfo

#include <ROOT/REveElement.hxx>

namespace ROOT {
namespace Experimental {

class REveViewer;
class REveScene;

/////////////////////////////////////////////////////////////////////////////
/// REveSceneInfo
/// Scene in a viewer.
/////////////////////////////////////////////////////////////////////////////

class REveSceneInfo : public REveElement
{
private:
   REveSceneInfo(const REveSceneInfo &) = delete;
   REveSceneInfo &operator=(const REveSceneInfo &) = delete;

protected:
   REveViewer *fViewer{nullptr};   ///<!
   REveScene *fScene{nullptr};     ///<!

public:
   REveSceneInfo(REveViewer *viewer, REveScene *scene);
   virtual ~REveSceneInfo() {}

   Int_t WriteCoreJson(nlohmann::json &j, Int_t rnr_offset) override;

   REveViewer *GetViewer() const { return fViewer; }
   REveScene *GetScene() const { return fScene; }

   Bool_t SingleRnrState() const override { return kTRUE; }

   Bool_t AcceptElement(REveElement *el) override;
};

} // namespace Experimental
} // namespace ROOT

#endif
