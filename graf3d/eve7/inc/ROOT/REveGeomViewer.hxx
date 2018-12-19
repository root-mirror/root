// @(#)root/eve:$Id$
// Author: Sergey Linev, 13.12.2018

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_REveGeomViewer
#define ROOT7_REveGeomViewer

#include <ROOT/RWebDisplayArgs.hxx>
#include <ROOT/REveGeomData.hxx>

#include <memory>

class TGeoManager;

namespace ROOT {
namespace Experimental {

class RWebWindow;
class REveManager;

class REveGeomViewer {

   friend class REveManager;

protected:

   TGeoManager *fGeoManager{nullptr};        ///<! geometry to show
   REveGeomDescription fDesc;                ///<! description

   std::shared_ptr<RWebWindow> fWebWindow;   ///<! web window to show geometry

   void WebWindowCallback(unsigned connid, const std::string &arg);

public:

   REveGeomViewer(TGeoManager *mgr);
   virtual ~REveGeomViewer();

   void SelectVolume(const std::string &volname);

   void SetDrawOptions(const std::string &opt) { fDesc.SetDrawOptions(opt); }

   void Show(const RWebDisplayArgs &args = "");

};

}}

#endif
