// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007, 2018

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_REveManager
#define ROOT7_REveManager

#include <ROOT/REveElement.hxx>

#include <ROOT/RWebDisplayArgs.hxx>

#include "TSysEvtHandler.h"
#include "TTimer.h"

#include <memory>
#include <unordered_map>

class TMap;
class TExMap;
class TMacro;
class TFolder;
class TGeoManager;

namespace ROOT {
namespace Experimental {

class REveSelection;
class REveViewer;
class REveViewerList;
class REveScene;
class REveSceneList;

class RWebWindow;
class REveGeomViewer;

class REveManager
{
   REveManager(const REveManager&) = delete;
   REveManager& operator=(const REveManager&) = delete;

public:
   class RRedrawDisabler {
   private:
      RRedrawDisabler(const RRedrawDisabler &) = delete;
      RRedrawDisabler &operator=(const RRedrawDisabler &) = delete;

      REveManager *fMgr{nullptr};

   public:
      RRedrawDisabler(REveManager *m) : fMgr(m)
      {
         if (fMgr)
            fMgr->DisableRedraw();
      }
      virtual ~RRedrawDisabler()
      {
         if (fMgr)
            fMgr->EnableRedraw();
      }
   };

   class RExceptionHandler : public TStdExceptionHandler {
   public:
      RExceptionHandler() : TStdExceptionHandler() { Add(); }
      virtual ~RExceptionHandler()                 { Remove(); }

      virtual EStatus  Handle(std::exception& exc);

      ClassDef(RExceptionHandler, 0);
   };

   struct Conn
   {
      unsigned fId{0};

      Conn() = default;
      Conn(unsigned int cId) : fId(cId) {}
   };

protected:
   RExceptionHandler        *fExcHandler{nullptr};   //!< exception handler

   TMap                     *fVizDB{nullptr};
   Bool_t                    fVizDBReplace{kFALSE};
   Bool_t                    fVizDBUpdate{kFALSE};

   TMap                     *fGeometries{nullptr};      //  TODO: use std::map<std::string, std::unique_ptr<TGeoManager>>
   TMap                     *fGeometryAliases{nullptr}; //  TODO: use std::map<std::string, std::string>

   TFolder                  *fMacroFolder{nullptr};

   REveScene                *fWorld{nullptr};

   REveViewerList           *fViewers{nullptr};
   REveSceneList            *fScenes{nullptr};

   REveScene                *fGlobalScene{nullptr};
   REveScene                *fEventScene{nullptr};

   Int_t                     fRedrawDisabled{0};
   Bool_t                    fFullRedraw{kFALSE};
   Bool_t                    fResetCameras{kFALSE};
   Bool_t                    fDropLogicals{kFALSE};
   Bool_t                    fKeepEmptyCont{kFALSE};
   Bool_t                    fTimerActive{kFALSE};
   TTimer                    fRedrawTimer;

   // ElementId management
   std::unordered_map<ElementId_t, REveElement*> fElementIdMap;
   ElementId_t                                   fLastElementId{0};
   ElementId_t                                   fNumElementIds{0};
   ElementId_t                                   fMaxElementIds{std::numeric_limits<ElementId_t>::max()};

   // Selection / highlight elements
   REveElement              *fSelectionList{nullptr};
   REveSelection            *fSelection{nullptr};
   REveSelection            *fHighlight{nullptr};

   std::shared_ptr<ROOT::Experimental::RWebWindow>  fWebWindow;
   std::vector<Conn>                                fConnList;

   void WindowConnect(unsigned connid);
   void WindowData(unsigned connid, const std::string &arg);
   void WindowDisconnect(unsigned connid);

public:
   REveManager(); // (Bool_t map_window=kTRUE, Option_t* opt="FI");
   virtual ~REveManager();

   RExceptionHandler *GetExcHandler() const { return fExcHandler; }

   REveSelection *GetSelection() const { return fSelection; }
   REveSelection *GetHighlight() const { return fHighlight; }

   REveSceneList *GetScenes() const { return fScenes; }
   REveViewerList *GetViewers() const { return fViewers; }

   REveScene *GetGlobalScene() const { return fGlobalScene; }
   REveScene *GetEventScene() const { return fEventScene; }

   REveScene *GetWorld() const { return fWorld; }

   REveViewer *SpawnNewViewer(const char *name, const char *title = "");
   REveScene *SpawnNewScene(const char *name, const char *title = "");

   TFolder *GetMacroFolder() const { return fMacroFolder; }
   TMacro *GetMacro(const char *name) const;

   void DisableRedraw() { ++fRedrawDisabled; }
   void EnableRedraw()  { --fRedrawDisabled; if (fRedrawDisabled <= 0) Redraw3D(); }

   void Redraw3D(Bool_t resetCameras = kFALSE, Bool_t dropLogicals = kFALSE)
   {
      if (fRedrawDisabled <= 0 && !fTimerActive) RegisterRedraw3D();
      if (resetCameras) fResetCameras = kTRUE;
      if (dropLogicals) fDropLogicals = kTRUE;
   }
   void RegisterRedraw3D();
   void DoRedraw3D();
   void FullRedraw3D(Bool_t resetCameras = kFALSE, Bool_t dropLogicals = kFALSE);

   void ClearAllSelections();

   Bool_t GetKeepEmptyCont() const   { return fKeepEmptyCont; }
   void   SetKeepEmptyCont(Bool_t k) { fKeepEmptyCont = k; }

   void AddElement(REveElement *element, REveElement *parent = nullptr);
   void AddGlobalElement(REveElement *element, REveElement *parent = nullptr);

   void RemoveElement(REveElement* element, REveElement *parent);

   REveElement *FindElementById (ElementId_t id) const;
   void         AssignElementId (REveElement* element);
   void         PreDeleteElement(REveElement* element);
   void         BrowseElement(ElementId_t id);

   // VizDB - Visualization-parameter data-base.
   Bool_t       InsertVizDBEntry(const TString& tag, REveElement* model,
                                 Bool_t replace, Bool_t update);
   Bool_t       InsertVizDBEntry(const TString& tag, REveElement* model);
   REveElement *FindVizDBEntry  (const TString& tag);

   void         LoadVizDB(const TString& filename, Bool_t replace, Bool_t update);
   void         LoadVizDB(const TString& filename);
   void         SaveVizDB(const TString& filename);

   Bool_t GetVizDBReplace()   const { return fVizDBReplace; }
   Bool_t GetVizDBUpdate ()   const { return fVizDBUpdate;  }
   void   SetVizDBReplace(Bool_t r) { fVizDBReplace = r; }
   void   SetVizDBUpdate (Bool_t u) { fVizDBUpdate  = u; }


   // Geometry management.
   TGeoManager *GetGeometry(const TString& filename);
   TGeoManager *GetGeometryByAlias(const TString& alias);
   TGeoManager *GetDefaultGeometry();
   void         RegisterGeometryAlias(const TString& alias, const TString& filename);

   void ClearROOTClassSaved();

   void AddLocation(const std::string& name, const std::string& path);
   void SetDefaultHtmlPage(const std::string& path);
   void SetClientVersion(const std::string& version);

   static REveManager* Create();
   static void         Terminate();

   // Access to internals, needed for low-level control in advanced
   // applications.

   void EnforceTimerActive (Bool_t ta) { fTimerActive = ta; }

   std::shared_ptr<RWebWindow> GetWebWindow() const { return fWebWindow; }

   // void Send(void* buff, unsigned connid);
   void Send(unsigned connid, const std::string &data);
   void SendBinary(unsigned connid, const void *data, std::size_t len);

   void DestroyElementsOf(REveElement::List_t &els);

   void BroadcastElementsOf(REveElement::List_t &els);

   void Show(const RWebDisplayArgs &args = "");

   std::shared_ptr<REveGeomViewer> ShowGeometry(const RWebDisplayArgs &args = "");
};

R__EXTERN REveManager* gEve;

}}

#endif
