// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007, 2018

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_REveViewer
#define ROOT7_REveViewer

#include <ROOT/REveElement.hxx>

namespace ROOT {
namespace Experimental {

class REveScene;

////////////////////////////////////////////////////////////////////////////////
/// REveViewer
/// Reve representation of TGLViewer.
////////////////////////////////////////////////////////////////////////////////

class REveViewer : public REveElement
{
private:
   REveViewer(const REveViewer&) = delete;
   REveViewer& operator=(const REveViewer&) = delete;

public:
   REveViewer(const std::string &n="REveViewer", const std::string &t="");
   virtual ~REveViewer();

   void Redraw(Bool_t resetCameras=kFALSE);

   virtual void AddScene(REveScene* scene);
   // XXX Missing RemoveScene() ????

   void RemoveElementLocal(REveElement *el) override;
   void RemoveElementsLocal() override;
};


////////////////////////////////////////////////////////////////////////////////
/// REveViewerList
/// List of Viewers providing common operations on REveViewer collections.
////////////////////////////////////////////////////////////////////////////////

class REveViewerList : public REveElement
{
private:
   REveViewerList(const REveViewerList&) = delete;
   REveViewerList& operator=(const REveViewerList&) = delete;

protected:
   Bool_t        fShowTooltip;

   Float_t       fBrightness;
   Bool_t        fUseLightColorSet;

   void HandleTooltip();

public:
   REveViewerList(const std::string &n="REveViewerList", const std::string &t="");
   virtual ~REveViewerList();

   void AddElement(REveElement* el) override;
   void RemoveElementLocal(REveElement* el) override;
   void RemoveElementsLocal() override;

   // --------------------------------

   virtual void Connect();
   virtual void Disconnect();

   void RepaintChangedViewers(Bool_t resetCameras, Bool_t dropLogicals);
   void RepaintAllViewers(Bool_t resetCameras, Bool_t dropLogicals);
   void DeleteAnnotations();

   void SceneDestructing(REveScene* scene);

   // --------------------------------

   void OnMouseOver(TObject* obj, UInt_t state);
   void OnReMouseOver(TObject* obj, UInt_t state);
   void OnUnMouseOver(TObject* obj, UInt_t state);
   void OnClicked(TObject *obj, UInt_t button, UInt_t state);
   void OnReClicked(TObject *obj, UInt_t button, UInt_t state);
   void OnUnClicked(TObject *obj, UInt_t button, UInt_t state);

   // --------------------------------

   Bool_t  GetShowTooltip()     const { return fShowTooltip; }
   void    SetShowTooltip(Bool_t x)   { fShowTooltip = x; }

   Float_t GetColorBrightness() const { return fBrightness; }
   void    SetColorBrightness(Float_t b);

   Bool_t  UseLightColorSet()   const { return fUseLightColorSet; }
   void    SwitchColorSet();
};

} // namespace Experimental
} // namespace ROOT

#endif
