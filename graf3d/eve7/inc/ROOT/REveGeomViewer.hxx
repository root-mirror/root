// @(#)root/eve7:$Id$
// Author: Sergey Linev, 13.12.2018

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
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
   std::string fSelectedVolume;              ///<! name of selected volume
   REveGeomDescription fDesc;                ///<! geometry description, send to the client as first message
   bool fShowHierarchy{true};                ///<! if hierarchy visible by default
   std::string fTitle;                       ///<! title of geometry viewer

   std::shared_ptr<RWebWindow> fWebWindow;   ///<! web window to show geometry

   void WebWindowCallback(unsigned connid, const std::string &arg);

   std::vector<int> GetStackFromJson(const std::string &json, bool node_ids = false);

   void SendGeometry(unsigned connid);

public:

   REveGeomViewer(TGeoManager *mgr = nullptr, const std::string &volname = "");
   virtual ~REveGeomViewer();

   void SetTitle(const std::string &title) { fTitle = title; }
   const std::string &GetTitle() const { return fTitle; }

   std::string GetWindowAddr() const;

   void SetGeometry(TGeoManager *mgr, const std::string &volname = "");

   void SelectVolume(const std::string &volname);

   /** Configures maximal number of visible nodes and faces */
   void SetLimits(int nnodes = 5000, int nfaces = 100000)
   {
      fDesc.SetMaxVisNodes(nnodes);
      fDesc.SetMaxVisFaces(nfaces);
   }

   /** Configures maximal visible level */
   void SetVisLevel(int lvl = 3)
   {
      fDesc.SetVisLevel(lvl);
   }

   /** Configures default hierarchy browser visibility, only has effect before showing web window */
   void SetShowHierarchy(bool on = true) { fShowHierarchy = on; }

   /** Returns default hierarchy browser visibility */
   bool GetShowHierarchy() const { return fShowHierarchy; }

   void SetDrawOptions(const std::string &opt);

   void Show(const RWebDisplayArgs &args = "", bool always_start_new_browser = false);

   void Update();

   void SaveImage(const std::string &fname = "geometry.png");

};

} // namespace Experimental
} // namespace ROOT

#endif
