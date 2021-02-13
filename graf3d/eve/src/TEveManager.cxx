// @(#)root/eve:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TEveManager.h"

#include "TEveSelection.h"
#include "TEveViewer.h"
#include "TEveScene.h"
#include "TEveEventManager.h"
#include "TEveWindowManager.h"

#include "TEveBrowser.h"
#include "TEveGedEditor.h"

#include "TGStatusBar.h"

#include "TGLSAViewer.h"

#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TObjString.h"
#include "TROOT.h"
#include "TFile.h"
#include "TMap.h"
#include "TExMap.h"
#include "TMacro.h"
#include "TFolder.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TColor.h"
#include "TPluginManager.h"
#include "TPRegexp.h"
#include "TClass.h"

#include <iostream>
#include <fstream>

TEveManager *gEve = nullptr;

/** \class TEveManager
\ingroup TEve
Central application manager for Eve.
Manages elements, GUI, GL scenes and GL viewers.
*/

ClassImp(TEveManager);

////////////////////////////////////////////////////////////////////////////////

TEveManager::TEveManager(UInt_t w, UInt_t h, Bool_t map_window, Option_t* opt) :
   fExcHandler  (0),
   fVizDB       (0), fVizDBReplace(kTRUE), fVizDBUpdate(kTRUE),
   fGeometries  (0),
   fGeometryAliases (0),
   fBrowser     (0),

   fMacroFolder (0),

   fWindowManager  (0),
   fViewers        (0),
   fScenes         (0),
   fGlobalScene    (0),
   fEventScene     (0),
   fCurrentEvent   (0),

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
   // If map_window is true, the TEveBrowser window is mapped.
   //
   // Option string is first parsed for the following characters:
   //   V - spawn a default GL viewer.
   //
   // The consumed characters are removed from the options and they
   // are passed to TEveBrowser for creation of additional plugins.
   //
   // Default options: "FIV" - file-browser, command-line, GL-viewer.


   static const TEveException eh("TEveManager::TEveManager ");

   if (gEve != 0)
      throw(eh + "There can be only one!");

   gEve = this;

   fExcHandler = new TExceptionHandler;

   fGeometries      = new TMap; fGeometries->SetOwnerKeyValue();
   fGeometryAliases = new TMap; fGeometryAliases->SetOwnerKeyValue();
   fVizDB           = new TMap; fVizDB->SetOwnerKeyValue();

   fStampedElements = new TExMap;

   fSelection = new TEveSelection("Global Selection");
   fSelection->IncDenyDestroy();
   fHighlight = new TEveSelection("Global Highlight");
   fHighlight->SetHighlightMode();
   fHighlight->IncDenyDestroy();

   fOrphanage = new TEveElementList("Global Orphanage");
   fOrphanage->IncDenyDestroy();

   fRedrawTimer.Connect("Timeout()", "TEveManager", this, "DoRedraw3D()");
   fMacroFolder = new TFolder("EVE", "Visualization macros");
   gROOT->GetListOfBrowsables()->Add(fMacroFolder);


   fWindowManager = new TEveWindowManager("WindowManager", "Manager of EVE windows");

   // Build GUI
   fBrowser   = new TEveBrowser(w, h);

   // ListTreeEditor
   fBrowser->StartEmbedding(0);
   fLTEFrame = new TEveGListTreeEditorFrame;
   fBrowser->StopEmbedding("Eve");
   fLTEFrame->ConnectSignals();

   // See how many GL viewers are requested, remove from options.
   TString str_opt(opt);
   TPMERegexp viewer_re("V", "g");
   Int_t viewer_count = viewer_re.Substitute(str_opt, "", kFALSE);

   // Create the main window / browse.
   fBrowser->InitPlugins(str_opt);
   if (map_window)
      fBrowser->MapWindow();

   // --------------------------------

   fWindowManager->IncDenyDestroy();
   AddToListTree(fWindowManager, kFALSE);

   fViewers = new TEveViewerList("Viewers");
   fViewers->IncDenyDestroy();
   AddToListTree(fViewers, kFALSE);

   fScenes  = new TEveSceneList ("Scenes");
   fScenes->IncDenyDestroy();
   AddToListTree(fScenes, kFALSE);

   fGlobalScene = new TEveScene("Geometry scene");
   fGlobalScene->IncDenyDestroy();
   fScenes->AddElement(fGlobalScene);

   fEventScene = new TEveScene("Event scene");
   fEventScene->IncDenyDestroy();
   fScenes->AddElement(fEventScene);

   for (Int_t vc = 0; vc < viewer_count; ++vc)
   {
      TEveViewer* v = SpawnNewViewer(Form("Viewer %d", vc+1));
      v->AddScene(fGlobalScene);
      v->AddScene(fEventScene);
   }

   if (GetDefaultViewer())
   {
      EditElement(GetDefaultViewer());
   }

   gSystem->ProcessEvents();
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TEveManager::~TEveManager()
{
   // Stop timer and deny further redraw requests.
   fRedrawTimer.Stop();
   fTimerActive = kTRUE;

   delete fCurrentEvent;
   fCurrentEvent = 0;

   fGlobalScene->DecDenyDestroy();
   fEventScene->DecDenyDestroy();
   fScenes->DestroyScenes();
   fScenes->DecDenyDestroy();
   fScenes->Destroy();
   fScenes = 0;

   fViewers->DestroyElements();
   fViewers->DecDenyDestroy();
   fViewers->Destroy();
   fViewers = 0;

   fWindowManager->DestroyWindows();
   fWindowManager->DecDenyDestroy();
   fWindowManager->Destroy();
   fWindowManager = 0;

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

   fLTEFrame->DeleteWindow();

   fBrowser->DontCallClose();
   fBrowser->TRootBrowser::CloseWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Clear the orphanage.

void TEveManager::ClearOrphanage()
{
   Bool_t old_state = fUseOrphanage;
   fUseOrphanage = kFALSE;
   fOrphanage->DestroyElements();
   fUseOrphanage = old_state;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the main window, i.e. EVE-browser.

TGWindow* TEveManager::GetMainWindow() const
{
   return fBrowser;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the default viewer - the first one in the fViewers list.

TEveViewer* TEveManager::GetDefaultViewer() const
{
   return dynamic_cast<TEveViewer*>(fViewers->FirstChild());
}

////////////////////////////////////////////////////////////////////////////////
/// Get TGLViewer of the default TEveViewer.

TGLViewer* TEveManager::GetDefaultGLViewer() const
{
   TEveViewer *ev = GetDefaultViewer();
   return ev ? ev->GetGLViewer() : 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns main object editor.

TEveGedEditor* TEveManager::GetEditor() const
{
   return fLTEFrame->GetEditor();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns main window status bar.

TGStatusBar* TEveManager::GetStatusBar() const
{
   return fBrowser->GetStatusBar();
}

////////////////////////////////////////////////////////////////////////////////
/// Add a new canvas tab.

TCanvas* TEveManager::AddCanvasTab(const char* name)
{
   fBrowser->StartEmbedding(1, -1);
   TCanvas* c = new TCanvas;
   fBrowser->StopEmbedding(name);

   return c;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a new GL viewer.

TEveViewer* TEveManager::SpawnNewViewer(const char* name, const char* title,
                                        Bool_t embed)
{
   TEveWindowSlot* slot = 0;
   if (embed)
   {
      slot = fWindowManager->GetCurrentWindowAsSlot();
      if (slot == 0)
      {
         // In principle should have some default/current container
         // in TEveWindowManager.
         // Also to store closed windows.
         slot = TEveWindow::CreateWindowInTab(fBrowser->GetTabRight());
         fBrowser->SanitizeTabCounts();
      }
   }
   else
   {
      slot = TEveWindow::CreateWindowMainFrame();
   }

   TEveViewer* v = new TEveViewer(name, title);
   v->SpawnGLViewer(embed ? GetEditor() : 0);

   slot->ReplaceWindow(v);

   fViewers->AddElement(v);

   return v;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a new scene.

TEveScene* TEveManager::SpawnNewScene(const char* name, const char* title)
{
   TEveScene* s = new TEveScene(name, title);
   AddElement(s, fScenes);
   return s;
}

////////////////////////////////////////////////////////////////////////////////
/// Find macro in fMacroFolder by name.

TMacro* TEveManager::GetMacro(const char* name) const
{
   return dynamic_cast<TMacro*>(fMacroFolder->FindObject(name));
}

////////////////////////////////////////////////////////////////////////////////
/// Show element in default editor.

void TEveManager::EditElement(TEveElement* element)
{
   static const TEveException eh("TEveManager::EditElement ");

   GetEditor()->DisplayElement(element);
}

////////////////////////////////////////////////////////////////////////////////
/// Register a request for 3D redraw.

void TEveManager::RegisterRedraw3D()
{
   fRedrawTimer.Start(0, kTRUE);
   fTimerActive = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Perform 3D redraw of scenes and viewers whose contents has
/// changed.

void TEveManager::DoRedraw3D()
{
   static const TEveException eh("TEveManager::DoRedraw3D ");

   // printf("TEveManager::DoRedraw3D redraw triggered\n");

   // Process element visibility changes, mark relevant scenes as changed.
   {
      TEveElement::List_t scenes;
      Long64_t   key, value;
      TExMapIter stamped_elements(fStampedElements);
      while (stamped_elements.Next(key, value))
      {
         TEveElement *el = reinterpret_cast<TEveElement*>(key);
         if (el->GetChangeBits() & TEveElement::kCBVisibility)
         {
            el->CollectSceneParents(scenes);
         }
      }
      ScenesChanged(scenes);
   }

   // Process changes in scenes.
   fScenes ->ProcessSceneChanges(fDropLogicals, fStampedElements);
   fViewers->RepaintChangedViewers(fResetCameras, fDropLogicals);

   // Process changed elements again, update GUI (just editor so far,
   // but more can come).
   {
      Long64_t   key, value;
      TExMapIter stamped_elements(fStampedElements);
      while (stamped_elements.Next(key, value))
      {
         TEveElement *el = reinterpret_cast<TEveElement*>(key);
         if (GetEditor()->GetModel() == el->GetEditorObject(eh))
            EditElement(el);
         TEveGedEditor::ElementChanged(el);

         el->ClearStamps();
      }
   }
   fStampedElements->Delete();
   GetListTree()->ClearViewPort(); // Fix this when several list-trees can be added.

   fResetCameras = kFALSE;
   fDropLogicals = kFALSE;

   fTimerActive = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Perform 3D redraw of all scenes and viewers.

void TEveManager::FullRedraw3D(Bool_t resetCameras, Bool_t dropLogicals)
{
   fScenes ->RepaintAllScenes (dropLogicals);
   fViewers->RepaintAllViewers(resetCameras, dropLogicals);
}

////////////////////////////////////////////////////////////////////////////////
/// Element was changed, perform framework side action.
/// Called from TEveElement::ElementChanged().

void TEveManager::ElementChanged(TEveElement* element, Bool_t update_scenes, Bool_t redraw)
{
   static const TEveException eh("TEveElement::ElementChanged ");

   if (GetEditor()->GetModel() == element->GetEditorObject(eh))
      EditElement(element);
   TEveGedEditor::ElementChanged(element);

   if (update_scenes) {
      TEveElement::List_t scenes;
      element->CollectSceneParents(scenes);
      ScenesChanged(scenes);
   }

   if (redraw)
      Redraw3D();
}

////////////////////////////////////////////////////////////////////////////////
/// Mark all scenes from the given list as changed.

void TEveManager::ScenesChanged(TEveElement::List_t& scenes)
{
   for (TEveElement::List_i s=scenes.begin(); s!=scenes.end(); ++s)
      ((TEveScene*)*s)->Changed();
}

////////////////////////////////////////////////////////////////////////////////
/// Mark element as changed -- it will be processed on next redraw.

void TEveManager::ElementStamped(TEveElement* element)
{
   UInt_t slot;
   if (fStampedElements->GetValue((ULong64_t) element, (Long64_t) element, slot) == 0)
   {
      fStampedElements->AddAt(slot, (ULong64_t) element, (Long64_t) element, 1);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get default list-tree widget.

TGListTree* TEveManager::GetListTree() const
{
   return fLTEFrame->fListTree;
}

TGListTreeItem*
TEveManager::AddToListTree(TEveElement* re, Bool_t open, TGListTree* lt)
{
   // Add element as a top-level to a list-tree.
   // Only add a single copy of a render-element as a top level.

   if (lt == 0) lt = GetListTree();
   TGListTreeItem* lti = re->AddIntoListTree(lt, (TGListTreeItem*)0);
   if (open) lt->OpenItem(lti);
   return lti;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove top-level element from list-tree with specified tree-item.

void TEveManager::RemoveFromListTree(TEveElement* element,
                                     TGListTree* lt, TGListTreeItem* lti)
{
   static const TEveException eh("TEveManager::RemoveFromListTree ");

   if (lti->GetParent())
      throw(eh + "not a top-level item.");

   element->RemoveFromListTree(lt, 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Add a new event and make it the current event.
/// It is added into the event-scene and as a top-level list-tree
/// item.

TGListTreeItem* TEveManager::AddEvent(TEveEventManager* event)
{
   fCurrentEvent = event;
   fCurrentEvent->IncDenyDestroy();
   AddElement(fCurrentEvent, fEventScene);
   return AddToListTree(event, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Add an element. If parent is not specified it is added into
/// current event (which is created if does not exist).

void TEveManager::AddElement(TEveElement* element, TEveElement* parent)
{
   if (parent == 0) {
      if (fCurrentEvent == 0)
         AddEvent(new TEveEventManager("Event", "Auto-created event directory"));
      parent = fCurrentEvent;
   }

   parent->AddElement(element);
}

////////////////////////////////////////////////////////////////////////////////
/// Add a global element, i.e. one that does not change on each
/// event, like geometry or projection manager.
/// If parent is not specified it is added to a global scene.

void TEveManager::AddGlobalElement(TEveElement* element, TEveElement* parent)
{
   if (parent == 0)
      parent = fGlobalScene;

   parent->AddElement(element);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove element from parent.

void TEveManager::RemoveElement(TEveElement* element,
                                TEveElement* parent)
{
   parent->RemoveElement(element);
}

////////////////////////////////////////////////////////////////////////////////
/// Called from TEveElement prior to its destruction so the
/// framework components (like object editor) can unreference it.

void TEveManager::PreDeleteElement(TEveElement* element)
{
   if (GetEditor()->GetEveElement() == element)
      EditElement(0);
   TEveGedEditor::ElementDeleted(element);

   if (fScenes)
      fScenes->DestroyElementRenderers(element);

   if (fStampedElements->GetValue((ULong64_t) element, (Long64_t) element) != 0)
      fStampedElements->Remove((ULong64_t) element, (Long64_t) element);

   if (element->fImpliedSelected > 0)
      fSelection->RemoveImpliedSelected(element);
   if (element->fImpliedHighlighted > 0)
      fHighlight->RemoveImpliedSelected(element);
}

////////////////////////////////////////////////////////////////////////////////
/// Select an element.
/// Now it only calls EditElement() - should also update selection state.

void TEveManager::ElementSelect(TEveElement* element)
{
   if (element != 0)
      EditElement(element);
}

////////////////////////////////////////////////////////////////////////////////
/// Paste has been called.

Bool_t TEveManager::ElementPaste(TEveElement* element)
{
   // The object to paste is taken from the editor (this is not
   // exactly right) and handed to 'element' for pasting.

   TEveElement* src = GetEditor()->GetEveElement();
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

Bool_t TEveManager::InsertVizDBEntry(const TString& tag, TEveElement* model,
                                     Bool_t replace, Bool_t update)
{
   TPair* pair = (TPair*) fVizDB->FindObject(tag);
   if (pair)
   {
      if (replace)
      {
         model->IncDenyDestroy();
         model->SetRnrChildren(kFALSE);

         TEveElement* old_model = dynamic_cast<TEveElement*>(pair->Value());
         if (old_model)
         {
            while (old_model->HasChildren())
            {
               TEveElement *el = old_model->FirstChild();
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

Bool_t TEveManager::InsertVizDBEntry(const TString& tag, TEveElement* model)
{
   return InsertVizDBEntry(tag, model, fVizDBReplace, fVizDBUpdate);
}

////////////////////////////////////////////////////////////////////////////////
/// Find a visualization-parameter database entry corresponding to tag.
/// If the entry is not found 0 is returned.

TEveElement* TEveManager::FindVizDBEntry(const TString& tag)
{
   return dynamic_cast<TEveElement*>(fVizDB->GetValue(tag));
}

////////////////////////////////////////////////////////////////////////////////
/// Load visualization-parameter database from file filename. The
/// replace, update arguments replace the values of fVizDBReplace
/// and fVizDBUpdate members for the duration of the macro
/// execution.

void TEveManager::LoadVizDB(const TString& filename, Bool_t replace, Bool_t update)
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

void TEveManager::LoadVizDB(const TString& filename)
{
   TEveUtil::Macro(filename);
   Redraw3D();
}

////////////////////////////////////////////////////////////////////////////////
/// Save visualization-parameter database to file filename.

void TEveManager::SaveVizDB(const TString& filename)
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
   out << "   TEveManager::Create();\n";

   ClearROOTClassSaved();

   Int_t       var_id = 0;
   TString     var_name;
   TIter       next(fVizDB);
   TObjString *key;
   while ((key = (TObjString*)next()))
   {
      TEveElement* mdl = dynamic_cast<TEveElement*>(fVizDB->GetValue(key));
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

TGeoManager* TEveManager::GetGeometry(const TString& filename)
{
   static const TEveException eh("TEveManager::GetGeometry ");

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
         throw(eh + "TGeoManager::Import() failed for '" + exp_filename + "'.");
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

TGeoManager* TEveManager::GetGeometryByAlias(const TString& alias)
{
   static const TEveException eh("TEveManager::GetGeometry ");

   TObjString* full_name = (TObjString*) fGeometryAliases->GetValue(alias);
   if (!full_name)
      throw(eh + "geometry alias '" + alias + "' not registered.");
   return GetGeometry(full_name->String());
}

////////////////////////////////////////////////////////////////////////////////
/// Get the default geometry.
/// It should be registered via RegisterGeometryName("Default", <URL>).

TGeoManager* TEveManager::GetDefaultGeometry()
{
   return GetGeometryByAlias("Default");
}

////////////////////////////////////////////////////////////////////////////////
/// Register 'name' as an alias for geometry file 'filename'.
/// The old aliases are silently overwritten.
/// After that the geometry can be retrieved also by calling:
///   gEve->GetGeometryByName(name);

void TEveManager::RegisterGeometryAlias(const TString& alias, const TString& filename)
{
   fGeometryAliases->Add(new TObjString(alias), new TObjString(filename));
}

////////////////////////////////////////////////////////////////////////////////
/// Set the text in the right side of browser's status bar.

void TEveManager::SetStatusLine(const char* text)
{
   fBrowser->SetStatusText(text, 1);
}

////////////////////////////////////////////////////////////////////////////////
/// Work-around uber ugly hack used in SavePrimitive and co.

void TEveManager::ClearROOTClassSaved()
{
   TIter   nextcl(gROOT->GetListOfClasses());
   TClass *cls;
   while((cls = (TClass *)nextcl()))
   {
      cls->ResetBit(TClass::kClassSaved);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Close button has been clicked on EVE main window (browser).
/// Cleanup and terminate application.

void TEveManager::CloseEveWindow()
{
   TGMainFrame *mf = (TGMainFrame*) gTQSender;
   TEveBrowser *eb = dynamic_cast<TEveBrowser*>(mf);
   if (eb == fBrowser)
   {
      mf->DontCallClose();
      Terminate();
      gApplication->Terminate();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// If global TEveManager* gEve is not set initialize it.
/// Returns gEve.

TEveManager* TEveManager::Create(Bool_t map_window, Option_t* opt)
{
   static const TEveException eh("TEveManager::Create ");

   if (gEve == 0)
   {
      // Make sure that the GUI system is initialized.
      if (gROOT->IsBatch())
      {
         throw eh + "ROOT is running in batch mode.";
      }
      TApplication::NeedGraphicsLibs();
      gApplication->InitializeGraphics();
      if (gROOT->IsBatch() || gClient == 0 || gClient->IsZombie())
      {
         throw eh + "window system not initialized.";
      }

      Int_t w = 1024;
      Int_t h =  768;

      TEveUtil::SetupEnvironment();
      TEveUtil::SetupGUI();
      gEve = new TEveManager(w, h, map_window, opt);
   }
   return gEve;
}

////////////////////////////////////////////////////////////////////////////////
/// Properly terminate global TEveManager.

void TEveManager::Terminate()
{
   if (!gEve) return;

   TEveGedEditor::DestroyEditors();

   delete gEve;
   gEve = 0;
}

/** \class TEveManager::TExceptionHandler
\ingroup TEve
Exception handler for Eve exceptions.
*/

ClassImp(TEveManager::TExceptionHandler);

////////////////////////////////////////////////////////////////////////////////
/// Handle exceptions deriving from TEveException.

TStdExceptionHandler::EStatus
TEveManager::TExceptionHandler::Handle(std::exception& exc)
{
   TEveException* ex = dynamic_cast<TEveException*>(&exc);
   if (ex) {
      Info("Handle", "%s", ex->Data());
      gEve->SetStatusLine(ex->Data());
      gSystem->Beep();
      return kSEHandled;
   } else {
      return kSEProceed;
   }
}
