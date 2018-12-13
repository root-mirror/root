// @(#)root/eve:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/REveManager.hxx>

#include <ROOT/REveSelection.hxx>
#include <ROOT/REveViewer.hxx>
#include <ROOT/REveScene.hxx>
#include <ROOT/REveClient.hxx>
#include <ROOT/REveGeomViewer.hxx>
#include <ROOT/RWebWindowsManager.hxx>

#include "TGeoManager.h"
#include "TObjString.h"
#include "TROOT.h"
#include "TFile.h"
#include "TMap.h"
#include "TExMap.h"
#include "TMacro.h"
#include "TFolder.h"
#include "TSystem.h"
#include "TRint.h"
#include "TEnv.h"
#include "TColor.h"
#include "TPluginManager.h"
#include "TPRegexp.h"
#include "TClass.h"
#include "THttpServer.h"

#include "Riostream.h"

#include "json.hpp"


using namespace ROOT::Experimental;
namespace REX = ROOT::Experimental;

REveManager* REX::gEve = 0;




/** \class REveManager
\ingroup REve
Central application manager for Eve.
Manages elements, GUI, GL scenes and GL viewers.
*/

////////////////////////////////////////////////////////////////////////////////

REveManager::REveManager() : // (Bool_t map_window, Option_t* opt) :
   fExcHandler  (nullptr),
   fVizDB       (nullptr),
   fVizDBReplace(kTRUE),
   fVizDBUpdate(kTRUE),
   fGeometries  (nullptr),
   fGeometryAliases (nullptr),

   fMacroFolder (nullptr),

   fViewers        (nullptr),
   fScenes         (nullptr),
   fGlobalScene    (nullptr),
   fEventScene     (nullptr),

   fRedrawDisabled (0),
   fResetCameras   (kFALSE),
   fDropLogicals   (kFALSE),
   fKeepEmptyCont  (kFALSE),
   fTimerActive    (kFALSE),
   fRedrawTimer    (),

   fStampedElements(0),
   fSelection      (0),
   fHighlight      (0),

   fOrphanage      (0),
   fUseOrphanage   (kFALSE)
{
   // Constructor.

   static const REveException eh("REveManager::REveManager ");

   if (REX::gEve != 0)
      throw eh + "There can be only one!";

   REX::gEve = this;

   fExcHandler = new RExceptionHandler;

   fGeometries      = new TMap; fGeometries->SetOwnerKeyValue();
   fGeometryAliases = new TMap; fGeometryAliases->SetOwnerKeyValue();
   fVizDB           = new TMap; fVizDB->SetOwnerKeyValue();

   fElementIdMap[0] = nullptr; // do not increase count for null element.

   fStampedElements = new TExMap;

   fSelection = new REveSelection("Global Selection");
   fSelection->IncDenyDestroy();
   AssignElementId(fSelection);
   fHighlight = new REveSelection("Global Highlight");
   fHighlight->SetHighlightMode();
   fHighlight->IncDenyDestroy();
   AssignElementId(fHighlight);

   fOrphanage = new REveElementList("Global Orphanage");
   fOrphanage->IncDenyDestroy();
   AssignElementId(fOrphanage);

   fRedrawTimer.Connect("Timeout()", "ROOT::Experimental::REveManager", this, "DoRedraw3D()");
   fMacroFolder = new TFolder("EVE", "Visualization macros");
   gROOT->GetListOfBrowsables()->Add(fMacroFolder);

   fWorld = new REveScene("EveWorld", "Top-level Eve Scene");
   fWorld->IncDenyDestroy();
   AssignElementId(fWorld);

   fViewers = new REveViewerList("Viewers");
   fViewers->IncDenyDestroy();
   fWorld->AddElement(fViewers);

   fScenes  = new REveSceneList ("Scenes");
   fScenes->IncDenyDestroy();
   fWorld->AddElement(fScenes);

   fGlobalScene = new REveScene("Geometry scene");
   fGlobalScene->IncDenyDestroy();
   fScenes->AddElement(fGlobalScene);

   fEventScene = new REveScene("Event scene");
   fEventScene->IncDenyDestroy();
   fScenes->AddElement(fEventScene);

   {
      REveViewer *v = SpawnNewViewer("Default Viewer");
      v->AddScene(fGlobalScene);
      v->AddScene(fEventScene);
   }

   // !!! AMT increase threshold to enable color pick on client
   TColor::SetColorThreshold(0.1);

   fWebWindow =  ROOT::Experimental::RWebWindowsManager::Instance()->CreateWindow();

   TString evedir = TString::Format("%s/eve7", TROOT::GetEtcDir().Data());
   if (gSystem->ExpandPathName(evedir)) {
      Warning("REveManager", "problems resolving %s for HTML sources", evedir.Data());
      evedir = ".";
   }

   fWebWindow->GetServer()->AddLocation("/evedir/",  evedir.Data());
   fWebWindow->SetDefaultPage(Form("file:%s/index.html", evedir.Data()));

   // this is call-back, invoked when message received via websocket
   fWebWindow->SetDataCallBack([this](unsigned connid, const std::string &arg) { this->HttpServerCallback(connid, arg); });
   fWebWindow->SetGeometry(900, 700); // configure predefined window geometry
   fWebWindow->SetConnLimit(100); // maximal number of connections
   fWebWindow->SetMaxQueueLength(30); // number of allowed entries in the window queue
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

REveManager::~REveManager()
{
   // Stop timer and deny further redraw requests.
   fRedrawTimer.Stop();
   fTimerActive = kTRUE;

   fGlobalScene->DecDenyDestroy();
   fEventScene->DecDenyDestroy();
   fScenes->DestroyScenes();
   fScenes->DecDenyDestroy();
   // Not needed - no more top-items: fScenes->Destroy();
   fScenes = nullptr;

   fViewers->DestroyElements();
   fViewers->DecDenyDestroy();
   // Not needed - no more top-items: fViewers->Destroy();
   fViewers = nullptr;

   // fWindowManager->DestroyWindows();
   // fWindowManager->DecDenyDestroy();
   // fWindowManager->Destroy();
   // fWindowManager = 0;

   fOrphanage->DecDenyDestroy();
   fHighlight->DecDenyDestroy();
   fSelection->DecDenyDestroy();

   gROOT->GetListOfBrowsables()->Remove(fMacroFolder);
   delete fMacroFolder;

   delete fGeometryAliases;
   delete fGeometries;
   delete fVizDB;
   delete fExcHandler;
   delete fStampedElements;
}

////////////////////////////////////////////////////////////////////////////////
/// Clear the orphanage.

void REveManager::ClearOrphanage()
{
   Bool_t old_state = fUseOrphanage;
   fUseOrphanage = kFALSE;
   fOrphanage->DestroyElements();
   fUseOrphanage = old_state;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a new GL viewer.

REveViewer* REveManager::SpawnNewViewer(const char* name, const char* title)
{
   REveViewer* v = new REveViewer(name, title);
   fViewers->AddElement(v);
   return v;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a new scene.

REveScene* REveManager::SpawnNewScene(const char* name, const char* title)
{
   REveScene* s = new REveScene(name, title);
   AddElement(s, fScenes);
   return s;
}

////////////////////////////////////////////////////////////////////////////////
/// Find macro in fMacroFolder by name.

TMacro* REveManager::GetMacro(const char* name) const
{
   return dynamic_cast<TMacro*>(fMacroFolder->FindObject(name));
}

////////////////////////////////////////////////////////////////////////////////
/// Show element in default editor.

void REveManager::EditElement(REveElement* /*element*/)
{
   static const REveException eh("REveManager::EditElement ");

   // GetEditor()->DisplayElement(element);
}

////////////////////////////////////////////////////////////////////////////////
/// Register a request for 3D redraw.

void REveManager::RegisterRedraw3D()
{
   fRedrawTimer.Start(0, kTRUE);
   fTimerActive = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Perform 3D redraw of scenes and viewers whose contents has
/// changed.

void REveManager::DoRedraw3D()
{
   static const REveException eh("REveManager::DoRedraw3D ");

   // Process changes in scenes.
   fScenes ->ProcessSceneChanges();


   fResetCameras = kFALSE;
   fDropLogicals = kFALSE;

   fTimerActive = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Perform 3D redraw of all scenes and viewers.

void REveManager::FullRedraw3D(Bool_t /*resetCameras*/, Bool_t /*dropLogicals*/)
{
   // XXXX fScenes ->RepaintAllScenes (dropLogicals);
   // XXXX fViewers->RepaintAllViewers(resetCameras, dropLogicals);
}

////////////////////////////////////////////////////////////////////////////////
/// Element was changed, perform framework side action.
/// Called from REveElement::ElementChanged().

void REveManager::ElementChanged(REveElement* element, Bool_t update_scenes, Bool_t redraw)
{
   static const REveException eh("REveElement::ElementChanged ");

   if (update_scenes) {
      REveElement::List_t scenes;
      element->CollectSceneParents(scenes);
      ScenesChanged(scenes);
   }

   if (redraw)
      Redraw3D();
}

////////////////////////////////////////////////////////////////////////////////
/// Mark all scenes from the given list as changed.

void REveManager::ScenesChanged(REveElement::List_t& scenes)
{
   for (REveElement::List_i s=scenes.begin(); s!=scenes.end(); ++s)
      ((REveScene*)*s)->Changed();
}

////////////////////////////////////////////////////////////////////////////////
/// Mark element as changed -- it will be processed on next redraw.

void REveManager::ElementStamped(REveElement* element)
{
   UInt_t slot;
   if (fStampedElements->GetValue((ULong64_t) element, (Long64_t) element, slot) == 0)
   {
      fStampedElements->AddAt(slot, (ULong64_t) element, (Long64_t) element, 1);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add an element. If parent is not specified it is added into
/// current event (which is created if does not exist).

void REveManager::AddElement(REveElement* element, REveElement* parent)
{
   if (parent == 0) {
      // XXXX
   }

   parent->AddElement(element);
}

////////////////////////////////////////////////////////////////////////////////
/// Add a global element, i.e. one that does not change on each
/// event, like geometry or projection manager.
/// If parent is not specified it is added to a global scene.

void REveManager::AddGlobalElement(REveElement* element, REveElement* parent)
{
   if (parent == 0)
      parent = fGlobalScene;

   parent->AddElement(element);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove element from parent.

void REveManager::RemoveElement(REveElement* element,
                                REveElement* parent)
{
   parent->RemoveElement(element);
}

////////////////////////////////////////////////////////////////////////////////
/// Lookup ElementId in element map and return corresponding REveElement*.
/// Returns nullptr if the id is not found

REveElement* REveManager::FindElementById(ElementId_t id) const
{
   static const REveException eh("REveManager::FindElementById ");

   auto it = fElementIdMap.find(id);
   return (it != fElementIdMap.end()) ? it->second : nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Assign a unique ElementId to given element.

void REveManager::AssignElementId(REveElement* element)
{
   static const REveException eh("REveManager::AssignElementId ");

   if (fNumElementIds == fMaxElementIds)
      throw eh + "ElementId map is full.";

   while (fElementIdMap.find(++fLastElementId) != fElementIdMap.end());

   element->fElementId = fLastElementId;
   fElementIdMap.insert(std::make_pair(fLastElementId, element));
   ++fNumElementIds;
}

////////////////////////////////////////////////////////////////////////////////
/// Called from REveElement prior to its destruction so the
/// framework components (like object editor) can unreference it.

void REveManager::PreDeleteElement(REveElement* element)
{
   // XXXX if (fScenes) fScenes->DestroyElementRenderers(element);

   if (fStampedElements->GetValue((ULong64_t) element, (Long64_t) element) != 0)
      fStampedElements->Remove((ULong64_t) element, (Long64_t) element);

   if (element->fImpliedSelected > 0)
      fSelection->RemoveImpliedSelected(element);
   if (element->fImpliedHighlighted > 0)
      fHighlight->RemoveImpliedSelected(element);

   if (element->fElementId != 0)
   {
      auto it = fElementIdMap.find(element->fElementId);
      if (it != fElementIdMap.end())
      {
         if (it->second == element)
         {
            fElementIdMap.erase(it);
            --fNumElementIds;
         }
         else Error("PreDeleteElement", "element ptr in ElementIdMap does not match the argument element.");
      }
      else Error("PreDeleteElement", "element id %u was not registered in ElementIdMap.", element->fElementId);
   }
   else Error("PreDeleteElement", "element with 0 ElementId passed in.");
}

////////////////////////////////////////////////////////////////////////////////
/// Select an element.
/// Now it only calls EditElement() - should also update selection state.

void REveManager::ElementSelect(REveElement* element)
{
   if (element != 0)
      EditElement(element);
}

////////////////////////////////////////////////////////////////////////////////
/// Paste has been called.

Bool_t REveManager::ElementPaste(REveElement* element)
{
   // The object to paste is taken from the editor (this is not
   // exactly right) and handed to 'element' for pasting.

   REveElement* src = 0; // GetEditor()->GetEveElement();
   if (src)
      return element->HandleElementPaste(src);
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Insert a new visualization-parameter database entry. Returns
/// true if the element is inserted successfully.
/// If entry with the same key already exists the behaviour depends on the
/// 'replace' flag:
///  - true  - The old model is deleted and new one is inserted (default).
///            Clients of the old model are transferred to the new one and
///            if 'update' flag is true (default), the new model's parameters
///            are assigned to all clients.
///  - false - The old model is kept, false is returned.
///
/// If insert is successful, the ownership of the model-element is
/// transferred to the manager.

Bool_t REveManager::InsertVizDBEntry(const TString& tag, REveElement* model,
                                     Bool_t replace, Bool_t update)
{
   TPair* pair = (TPair*) fVizDB->FindObject(tag);
   if (pair)
   {
      if (replace)
      {
         model->IncDenyDestroy();
         model->SetRnrChildren(kFALSE);

         REveElement* old_model = dynamic_cast<REveElement*>(pair->Value());
         if (old_model)
         {
            while (old_model->HasChildren())
            {
               REveElement *el = old_model->FirstChild();
               el->SetVizModel(model);
               if (update)
               {
                  el->CopyVizParams(model);
                  el->PropagateVizParamsToProjecteds();
               }
            }
            old_model->DecDenyDestroy();
         }
         pair->SetValue(dynamic_cast<TObject*>(model));
         return kTRUE;
      }
      else
      {
         return kFALSE;
      }
   }
   else
   {
      model->IncDenyDestroy();
      model->SetRnrChildren(kFALSE);
      fVizDB->Add(new TObjString(tag), dynamic_cast<TObject*>(model));
      return kTRUE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Insert a new visualization-parameter database entry with the default
/// parameters for replace and update, as specified by members
/// fVizDBReplace(default=kTRUE) and fVizDBUpdate(default=kTRUE).
/// See docs of the above function.

Bool_t REveManager::InsertVizDBEntry(const TString& tag, REveElement* model)
{
   return InsertVizDBEntry(tag, model, fVizDBReplace, fVizDBUpdate);
}

////////////////////////////////////////////////////////////////////////////////
/// Find a visualization-parameter database entry corresponding to tag.
/// If the entry is not found 0 is returned.

REveElement* REveManager::FindVizDBEntry(const TString& tag)
{
   return dynamic_cast<REveElement*>(fVizDB->GetValue(tag));
}

////////////////////////////////////////////////////////////////////////////////
/// Load visualization-parameter database from file filename. The
/// replace, update arguments replace the values of fVizDBReplace
/// and fVizDBUpdate members for the duration of the macro
/// execution.

void REveManager::LoadVizDB(const TString& filename, Bool_t replace, Bool_t update)
{
   Bool_t ex_replace = fVizDBReplace;
   Bool_t ex_update  = fVizDBUpdate;
   fVizDBReplace = replace;
   fVizDBUpdate  = update;

   LoadVizDB(filename);

   fVizDBReplace = ex_replace;
   fVizDBUpdate  = ex_update;
}

////////////////////////////////////////////////////////////////////////////////
/// Load visualization-parameter database from file filename.
/// State of data-members fVizDBReplace and fVizDBUpdate determine
/// how the registered entries are handled.

void REveManager::LoadVizDB(const TString& filename)
{
   REveUtil::Macro(filename);
   Redraw3D();
}

////////////////////////////////////////////////////////////////////////////////
/// Save visualization-parameter database to file filename.

void REveManager::SaveVizDB(const TString& filename)
{
   TPMERegexp re("(.+)\\.\\w+");
   if (re.Match(filename) != 2) {
      Error("SaveVizDB", "filename does not match required format '(.+)\\.\\w+'.");
      return;
   }

   TString exp_filename(filename);
   gSystem->ExpandPathName(exp_filename);

   std::ofstream out(exp_filename, std::ios::out | std::ios::trunc);
   out << "void " << re[1] << "()\n";
   out << "{\n";
   out << "   REveManager::Create();\n";

   ClearROOTClassSaved();

   Int_t       var_id = 0;
   TString     var_name;
   TIter       next(fVizDB);
   TObjString *key;
   while ((key = (TObjString*)next()))
   {
      REveElement* mdl = dynamic_cast<REveElement*>(fVizDB->GetValue(key));
      if (mdl)
      {
         var_name.Form("x%03d", var_id++);
         mdl->SaveVizParams(out, key->String(), var_name);
      }
      else
      {
         Warning("SaveVizDB", "Saving failed for key '%s'.", key->String().Data());
      }
   }

   out << "}\n";
   out.close();
}

////////////////////////////////////////////////////////////////////////////////
/// Get geometry with given filename.
/// This is cached internally so the second time this function is
/// called with the same argument the same geo-manager is returned.
/// gGeoManager is set to the return value.

TGeoManager* REveManager::GetGeometry(const TString& filename)
{
   static const REveException eh("REveManager::GetGeometry ");

   TString exp_filename = filename;
   gSystem->ExpandPathName(exp_filename);
   printf("%s loading: '%s' -> '%s'.\n", eh.Data(),
          filename.Data(), exp_filename.Data());

   gGeoManager = (TGeoManager*) fGeometries->GetValue(filename);
   if (gGeoManager)
   {
      gGeoIdentity = (TGeoIdentity*) gGeoManager->GetListOfMatrices()->At(0);
   }
   else
   {
      Bool_t locked = TGeoManager::IsLocked();
      if (locked) {
         Warning(eh, "TGeoManager is locked ... unlocking it.");
         TGeoManager::UnlockGeometry();
      }
      if (TGeoManager::Import(filename) == 0) {
         throw eh + "TGeoManager::Import() failed for '" + exp_filename + "'.";
      }
      if (locked) {
         TGeoManager::LockGeometry();
      }

      gGeoManager->GetTopVolume()->VisibleDaughters(1);

      // Import colors exported by Gled, if they exist.
      {
         TFile f(exp_filename, "READ");
         TObjArray* collist = (TObjArray*) f.Get("ColorList");
         f.Close();
         if (collist != 0) {
            TIter next(gGeoManager->GetListOfVolumes());
            TGeoVolume* vol;
            while ((vol = (TGeoVolume*) next()) != 0)
            {
               Int_t oldID = vol->GetLineColor();
               TColor* col = (TColor*)collist->At(oldID);
               Float_t r, g, b;
               col->GetRGB(r, g, b);
               Int_t  newID = TColor::GetColor(r,g,b);
               vol->SetLineColor(newID);
            }
         }
      }

      fGeometries->Add(new TObjString(filename), gGeoManager);
   }
   return gGeoManager;
}

////////////////////////////////////////////////////////////////////////////////
/// Get geometry with given alias.
/// The alias must be registered via RegisterGeometryAlias().

TGeoManager* REveManager::GetGeometryByAlias(const TString& alias)
{
   static const REveException eh("REveManager::GetGeometry ");

   TObjString* full_name = (TObjString*) fGeometryAliases->GetValue(alias);
   if (!full_name)
      throw eh + "geometry alias '" + alias + "' not registered.";
   return GetGeometry(full_name->String());
}

////////////////////////////////////////////////////////////////////////////////
/// Get the default geometry.
/// It should be registered via RegisterGeometryName("Default", <URL>).

TGeoManager* REveManager::GetDefaultGeometry()
{
   return GetGeometryByAlias("Default");
}

////////////////////////////////////////////////////////////////////////////////
/// Register 'name' as an alias for geometry file 'filename'.
/// The old aliases are silently overwritten.
/// After that the geometry can be retrieved also by calling:
///   REX::gEve->GetGeometryByName(name);

void REveManager::RegisterGeometryAlias(const TString& alias, const TString& filename)
{
   fGeometryAliases->Add(new TObjString(alias), new TObjString(filename));
}

////////////////////////////////////////////////////////////////////////////////
/// Work-around uber ugly hack used in SavePrimitive and co.

void REveManager::ClearROOTClassSaved()
{
   TIter   nextcl(gROOT->GetListOfClasses());
   TClass *cls;
   while((cls = (TClass *)nextcl()))
   {
      cls->ResetBit(TClass::kClassSaved);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// If global REveManager* REX::gEve is not set initialize it.
/// Returns REX::gEve.

REveManager* REveManager::Create()
{
   static const REveException eh("REveManager::Create ");

   if (!REX::gEve)
   {
      // XXXX Initialize some server stuff ???

      REX::gEve = new REveManager();
   }
   return REX::gEve;
}

////////////////////////////////////////////////////////////////////////////////
/// Properly terminate global REveManager.

void REveManager::Terminate()
{
   if (!REX::gEve) return;

   delete REX::gEve;
   REX::gEve = nullptr;
}

/** \class REveManager::RExceptionHandler
\ingroup REve
Exception handler for Eve exceptions.
*/


////////////////////////////////////////////////////////////////////////////////
/// Handle exceptions deriving from REveException.

TStdExceptionHandler::EStatus
REveManager::RExceptionHandler::Handle(std::exception& exc)
{
   REveException* ex = dynamic_cast<REveException*>(&exc);
   if (ex) {
      Info("Handle", "%s", ex->Data());
      // REX::gEve->SetStatusLine(ex->Data());
      gSystem->Beep();
      return kSEHandled;
   } else {
      return kSEProceed;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Callback from THttpServer

void REveManager::HttpServerCallback(unsigned connid, const std::string &arg)
{
   if (arg == "CONN_READY")
   {
      fConnList.emplace_back(connid);
      printf("connection established %u\n", connid);

      printf("\nEVEMNG ............. streaming the world scene.\n");

      // This prepares core and render data buffers.
      fWorld->StreamElements();

      printf("   sending json, len = %d\n", (int) fWorld->fOutputJson.size());
      Send(connid, fWorld->fOutputJson);
      printf("   for now assume world-scene has no render data, binary-size=%d\n", fWorld->fTotalBinarySize);
      assert(fWorld->fTotalBinarySize == 0);

      for (REveElement::List_i it = fScenes->BeginChildren(); it != fScenes->EndChildren(); ++it)
      {
         REveScene* scene = dynamic_cast<REveScene*>(*it);

         scene->AddSubscriber(std::make_unique<REveClient>(connid, fWebWindow));
         printf("\nEVEMNG ............. streaming scene %s [%s]\n",
                scene->GetElementTitle() ,scene->GetElementName());

         // This prepares core and render data buffers.
         scene->StreamElements();

         printf("   sending json, len = %d\n", (int) scene->fOutputJson.size());
         Send(connid, scene->fOutputJson);
         printf("   sending binary, len = %d\n", scene->fTotalBinarySize);
         SendBinary(connid, &scene->fOutputBinary[0], scene->fTotalBinarySize);
      }
      return;
   }

   // find connection object
   std::vector<Conn>::iterator conn = fConnList.end();
   for (auto i = fConnList.begin(); i != fConnList.end(); ++i)
   {
      if (i->fId == connid)
      {
         conn = i;
         break;
      }
   }
   // this should not happen, just check
   if (conn == fConnList.end()) {
      printf("error, connection not found!");
      return;
   }

   if (arg == "CONN_CLOSED") {
      printf("connection closed\n");
      fConnList.erase(conn);
      for (auto i = fScenes->BeginChildren(); i != fScenes->EndChildren(); ++i)
      {
         REveScene* scene = dynamic_cast<REveScene*>(*i);
         scene->RemoveSubscriber(connid);
      }

      return;
   }
   else
   {
      fScenes->AcceptChanges(true);

      // MIR
      nlohmann::json cj =  nlohmann::json::parse(arg.c_str());
      printf("MIR test %s \n", cj.dump().c_str());
      std::string mir =  cj["mir"];
      std::string ctype =  cj["class"];
      int id = cj["fElementId"];

      auto el =  FindElementById(id);
      char cmd[128];
      sprintf(cmd, "((%s*)%p)->%s;", ctype.c_str(), el, mir.c_str());
      printf("MIR cmd %s\n", cmd);
      gROOT->ProcessLine(cmd);

      fScenes->AcceptChanges(false);
      Redraw3D();


      /*
      nlohmann::json resp;
      resp["function"] = "replaceElement";
      //el->SetCoreJson(resp);
      for (auto i = fConnList.begin(); i != fConnList.end(); ++i)
      {
         fWebWindow->Send(i->fId, resp.dump());
      }
      */
   }
}

void REveManager::Send(unsigned connid, const std::string &data)
{
   fWebWindow->Send(connid, data);
}


void REveManager::SendBinary(unsigned connid, const void *data, std::size_t len)
{
   fWebWindow->SendBinary(connid, data, len);
}

//------------------------------------------------------------------------------

void REveManager::DestroyElementsOf(REveElement::List_t& els)
{
   fScenes->AcceptChanges(false);

   nlohmann::json jarr = nlohmann::json::array();

   nlohmann::json jhdr = {};
   jhdr["content"]  = "REveManager::DestroyElementsOf";

   nlohmann::json jels = nlohmann::json::array();

   for (auto & ep : els)
   {
      jels.push_back(ep->GetElementId());

      ep->DestroyElements();
   }

   jhdr["element_ids"] = jels;

   jarr.push_back(jhdr);

   std::string msg = jarr.dump();

   // XXXX Do we have broadcast?

   for (auto && conn: fConnList)
   {
      fWebWindow->Send(conn.fId, msg);
   }
}

void REveManager::BroadcastElementsOf(REveElement::List_t& els)
{
   for (auto & ep : els)
   {
      REveScene* scene = dynamic_cast<REveScene*>(ep);
      assert (scene != 0);

      printf("\nEVEMNG ............. streaming scene %s [%s]\n",
             scene->GetElementTitle(), scene->GetElementName());

      // This prepares core and render data buffers.
      scene->StreamElements();

      for (auto i = fConnList.begin(); i != fConnList.end(); ++i)
      {
         printf("   sending json, len = %d --> to conn_id = %d\n", (int) scene->fOutputJson.size(), i->fId);
         fWebWindow->Send(i->fId, scene->fOutputJson);
         printf("   sending binary, len = %d --> to conn_id = %d\n", scene->fTotalBinarySize, i->fId);
         fWebWindow->SendBinary(i->fId, &scene->fOutputBinary[0], scene->fTotalBinarySize);
      }
   }

   // AMT: This call may not be necessary
   fScenes->AcceptChanges(true);
}

//////////////////////////////////////////////////////////////////
/// Show eve manager in specified browser
/// If rootrc variable WebEve.DisableShow configured, just HTTP server will be started and URL printout

void REveManager::Show(const RWebDisplayArgs &args)
{
   if (gEnv->GetValue("WebEve.DisableShow", 0) != 0) {
      std::string url = fWebWindow->GetUrl(true);
      printf("EVE URL %s\n", url.c_str());
   } else {
      fWebWindow->Show(args);
   }
}

//////////////////////////////////////////////////////////////////
/// Show current geometry in web browser

std::shared_ptr<REveGeomViewer> REveManager::ShowGeometry(const RWebDisplayArgs &args)
{
   if (!gGeoManager) {
      Error("ShowGeometry", "No geometry is loaded");
      return nullptr;
   }

   auto viewer = std::make_shared<REveGeomViewer>(gGeoManager);

   viewer->Show(args);

   return viewer;
}
