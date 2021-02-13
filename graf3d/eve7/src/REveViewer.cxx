// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/REveViewer.hxx>

#include <ROOT/REveUtil.hxx>
#include <ROOT/REveScene.hxx>
#include <ROOT/REveSceneInfo.hxx>
#include <ROOT/REveManager.hxx>
#include <ROOT/REveSelection.hxx>

using namespace ROOT::Experimental;
namespace REX = ROOT::Experimental;

/** \class REveViewer
\ingroup REve
Eve representation of a GL view. In a gist, it's a camera + a list of scenes.

*/

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

REveViewer::REveViewer(const std::string& n, const std::string& t) :
   REveElement(n, t)
{
   // SetChildClass(TClass::GetClass<REveSceneInfo>());
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

REveViewer::~REveViewer()
{}

////////////////////////////////////////////////////////////////////////////////
/// Redraw viewer immediately.

void REveViewer::Redraw(Bool_t /*resetCameras*/)
{
   // if (resetCameras) fGLViewer->PostSceneBuildSetup(kTRUE);
   // fGLViewer->RequestDraw(TGLRnrCtx::kLODHigh);
}

////////////////////////////////////////////////////////////////////////////////
/// Add 'scene' to the list of scenes.

void REveViewer::AddScene(REveScene *scene)
{
   static const REveException eh("REveViewer::AddScene ");

   for (auto &c: RefChildren()) {
      auto sinfo = dynamic_cast<REveSceneInfo*>(c);

      if (sinfo && sinfo->GetScene() == scene)
      {
         throw eh + "scene already in the viewer.";
      }
   }

   auto si = new REveSceneInfo(this, scene);
   AddElement(si);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove element 'el' from the list of children and also remove
/// appropriate GLScene from GLViewer's list of scenes.
/// Virtual from REveElement.

void REveViewer::RemoveElementLocal(REveElement* /*el*/)
{
   // fGLViewer->RemoveScene(((REveSceneInfo*)el)->GetGLScene());

   // XXXXX Notify clients !!! Or will this be automatic?
}

////////////////////////////////////////////////////////////////////////////////
/// Remove all children, forwarded to GLViewer.
/// Virtual from REveElement.

void REveViewer::RemoveElementsLocal()
{
   // fGLViewer->RemoveAllScenes();

   // XXXXX Notify clients !!! Or will this be automatic?
}


/** \class REveViewerList
\ingroup REve
List of Viewers providing common operations on REveViewer collections.
*/

////////////////////////////////////////////////////////////////////////////////

REveViewerList::REveViewerList(const std::string &n, const std::string &t) :
   REveElement  (n, t),
   fShowTooltip (kTRUE),

   fBrightness(0),
   fUseLightColorSet(kFALSE)
{
   // Constructor.

   SetChildClass(TClass::GetClass<REveViewer>());
   Connect();
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

REveViewerList::~REveViewerList()
{
   Disconnect();
}

////////////////////////////////////////////////////////////////////////////////
/// Call base-class implementation.
/// If compound is open and compound of the new element is not set,
/// the el's compound is set to this.

void REveViewerList::AddElement(REveElement* el)
{
   REveElement::AddElement(el);
}

////////////////////////////////////////////////////////////////////////////////
/// Decompoundofy el, call base-class version.

void REveViewerList::RemoveElementLocal(REveElement* el)
{
   // This was needed as viewer was in EveWindowManager hierarchy, too.
   // el->DecParentIgnoreCnt();

   REveElement::RemoveElementLocal(el);
}

////////////////////////////////////////////////////////////////////////////////
/// Decompoundofy children, call base-class version.

void REveViewerList::RemoveElementsLocal()
{
   // This was needed as viewer was in EveWindowManager hierarchy, too.
   // el->DecParentIgnoreCnt();
   // for (auto &c: fChildren)
   // {
   //    c->DecParentIgnoreCnt();
   // }

   REveElement::RemoveElementsLocal();
}

////////////////////////////////////////////////////////////////////////////////
/// Connect to TGLViewer class-signals.

void REveViewerList::Connect()
{
   // TQObject::Connect("TGLViewer", "MouseOver(TObject*,UInt_t)",
   //                   "REveViewerList", this, "OnMouseOver(TObject*,UInt_t)");

   // TQObject::Connect("TGLViewer", "ReMouseOver(TObject*,UInt_t)",
   //                   "REveViewerList", this, "OnReMouseOver(TObject*,UInt_t)");

   // TQObject::Connect("TGLViewer", "UnMouseOver(TObject*,UInt_t)",
   //                   "REveViewerList", this, "OnUnMouseOver(TObject*,UInt_t)");

   // TQObject::Connect("TGLViewer", "Clicked(TObject*,UInt_t,UInt_t)",
   //                   "REveViewerList", this, "OnClicked(TObject*,UInt_t,UInt_t)");

   // TQObject::Connect("TGLViewer", "ReClicked(TObject*,UInt_t,UInt_t)",
   //                   "REveViewerList", this, "OnReClicked(TObject*,UInt_t,UInt_t)");

   // TQObject::Connect("TGLViewer", "UnClicked(TObject*,UInt_t,UInt_t)",
   //                   "REveViewerList", this, "OnUnClicked(TObject*,UInt_t,UInt_t)");
}

////////////////////////////////////////////////////////////////////////////////
/// Disconnect from TGLViewer class-signals.

void REveViewerList::Disconnect()
{
   // TQObject::Disconnect("TGLViewer", "MouseOver(TObject*,UInt_t)",
   //                      this, "OnMouseOver(TObject*,UInt_t)");

   // TQObject::Disconnect("TGLViewer", "ReMouseOver(TObject*,UInt_t)",
   //                      this, "OnReMouseOver(TObject*,UInt_t)");

   // TQObject::Disconnect("TGLViewer", "UnMouseOver(TObject*,UInt_t)",
   //                      this, "OnUnMouseOver(TObject*,UInt_t)");

   // TQObject::Disconnect("TGLViewer", "Clicked(TObject*,UInt_t,UInt_t)",
   //                      this, "OnClicked(TObject*,UInt_t,UInt_t)");

   // TQObject::Disconnect("TGLViewer", "ReClicked(TObject*,UInt_t,UInt_t)",
   //                      this, "OnReClicked(TObject*,UInt_t,UInt_t)");

   // TQObject::Disconnect("TGLViewer", "UnClicked(TObject*,UInt_t,UInt_t)",
   //                      this, "OnUnClicked(TObject*,UInt_t,UInt_t)");
}

////////////////////////////////////////////////////////////////////////////////
/// Repaint viewers that are tagged as changed.

void REveViewerList::RepaintChangedViewers(Bool_t /*resetCameras*/, Bool_t /*dropLogicals*/)
{
   //for (auto &c: fChildren)  {
      // TGLViewer* glv = ((REveViewer*)c)->GetGLViewer();
      // if (glv->IsChanged())
      // {
      //    if (resetCameras) glv->PostSceneBuildSetup(kTRUE);
      //    if (dropLogicals) glv->SetSmartRefresh(kFALSE);

      //    glv->RequestDraw(TGLRnrCtx::kLODHigh);

      //    if (dropLogicals) glv->SetSmartRefresh(kTRUE);
      // }
   //}
}

////////////////////////////////////////////////////////////////////////////////
/// Repaint all viewers.

void REveViewerList::RepaintAllViewers(Bool_t /*resetCameras*/, Bool_t /*dropLogicals*/)
{
   // for (auto &c: fChildren) {
      // TGLViewer* glv = ((REveViewer *)c)->GetGLViewer();

      // if (resetCameras) glv->PostSceneBuildSetup(kTRUE);
      // if (dropLogicals) glv->SetSmartRefresh(kFALSE);

      // glv->RequestDraw(TGLRnrCtx::kLODHigh);

      // if (dropLogicals) glv->SetSmartRefresh(kTRUE);
   // }
}

////////////////////////////////////////////////////////////////////////////////
/// Delete annotations from all viewers.

void REveViewerList::DeleteAnnotations()
{
   // for (auto &c: fChildren) {
      // TGLViewer* glv = ((REveViewer *)c)->GetGLViewer();
      // glv->DeleteOverlayAnnotations();
  // }
}

////////////////////////////////////////////////////////////////////////////////
/// Callback done from a REveScene destructor allowing proper
/// removal of the scene from affected viewers.

void REveViewerList::SceneDestructing(REveScene* scene)
{
   for (auto &viewer: fChildren) {
      for (auto &j: viewer->RefChildren()) {
         REveSceneInfo* sinfo = (REveSceneInfo *) j;
         if (sinfo->GetScene() == scene)
            viewer->RemoveElement(sinfo);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Show / hide tooltip for various MouseOver events.
/// Must be called from slots where sender is TGLEventHandler.

void REveViewerList::HandleTooltip()
{
   if (fShowTooltip)
   {
      // TGLViewer       *glw = dynamic_cast<TGLViewer*>((TQObject*) gTQSender);
      // TGLEventHandler *glh = (TGLEventHandler*) glw->GetEventHandler();
      // if (REX::gEve->GetHighlight()->NumChildren() == 1)
      // {
      //    TString title(REX::gEve->GetHighlight()->FirstChild()->GetHighlightTooltip());
      //    if ( ! title.IsNull())
      //       glh->TriggerTooltip(title);
      // }
      // else
      // {
      //    glh->RemoveTooltip();
      // }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for global TGLViewer::MouseOver() signal.
///
/// The attempt is made to determine the REveElement being
/// represented by the physical shape and global highlight is updated
/// accordingly.
///
/// If REveElement::IsPickable() returns false, the element is not
/// highlighted.
///
/// Highlight is always in single-selection mode.

void REveViewerList::OnMouseOver(TObject *, UInt_t /*state*/)
{
   // REveElement *el = dynamic_cast<REveElement*>(obj);
   // if (el && !el->IsPickable())
   //   el = nullptr;

   // void *qsender = gTQSender;
   // REX::gEve->GetHighlight()->UserPickedElement(el, kFALSE);
   // gTQSender = qsender;

   HandleTooltip();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for global TGLViewer::ReMouseOver().
///
/// The obj is dyn-casted to the REveElement and global selection is
/// updated accordingly.
///
/// If REveElement::IsPickable() returns false, the element is not
/// selected.

void REveViewerList::OnReMouseOver(TObject *obj, UInt_t /*state*/)
{
   REveElement* el = dynamic_cast<REveElement*>(obj);
   if (el && !el->IsPickable())
      el = 0;

   // void *qsender = gTQSender;
   // REX::gEve->GetHighlight()->UserRePickedElement(el);
   // gTQSender = qsender;

   HandleTooltip();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for global TGLViewer::UnMouseOver().
///
/// The obj is dyn-casted to the REveElement and global selection is
/// updated accordingly.
///
/// If REveElement::IsPickable() returns false, the element is not
/// selected.

void REveViewerList::OnUnMouseOver(TObject *obj, UInt_t /*state*/)
{
   REveElement* el = dynamic_cast<REveElement*>(obj);
   if (el && !el->IsPickable())
      el = 0;

   // void *qsender = gTQSender;
   // REX::gEve->GetHighlight()->UserUnPickedElement(el);
   // gTQSender = qsender;

   HandleTooltip();
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for global TGLViewer::Clicked().
///
/// The obj is dyn-casted to the REveElement and global selection is
/// updated accordingly.
///
/// If REveElement::IsPickable() returns false, the element is not
/// selected.

void REveViewerList::OnClicked(TObject *obj, UInt_t /*button*/, UInt_t state)
{
   REveElement* el = dynamic_cast<REveElement*>(obj);
   if (el && !el->IsPickable())
      el = 0;
   REX::gEve->GetSelection()->UserPickedElement(el, state & kKeyControlMask);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for global TGLViewer::ReClicked().
///
/// The obj is dyn-casted to the REveElement and global selection is
/// updated accordingly.
///
/// If REveElement::IsPickable() returns false, the element is not
/// selected.

void REveViewerList::OnReClicked(TObject *obj, UInt_t /*button*/, UInt_t /*state*/)
{
   REveElement* el = dynamic_cast<REveElement*>(obj);
   if (el && !el->IsPickable())
      el = 0;
   REX::gEve->GetSelection()->UserRePickedElement(el);
}

////////////////////////////////////////////////////////////////////////////////
/// Slot for global TGLViewer::UnClicked().
///
/// The obj is dyn-casted to the REveElement and global selection is
/// updated accordingly.
///
/// If REveElement::IsPickable() returns false, the element is not
/// selected.

void REveViewerList::OnUnClicked(TObject *obj, UInt_t /*button*/, UInt_t /*state*/)
{
   REveElement* el = dynamic_cast<REveElement*>(obj);
   if (el && !el->IsPickable())
      el = 0;
   REX::gEve->GetSelection()->UserUnPickedElement(el);
}

////////////////////////////////////////////////////////////////////////////////
/// Set color brightness.

void REveViewerList::SetColorBrightness(Float_t b)
{
   REveUtil::SetColorBrightness(b, 1);
}

////////////////////////////////////////////////////////////////////////////////
/// Switch background color.

void REveViewerList::SwitchColorSet()
{
   fUseLightColorSet = ! fUseLightColorSet;
   // for (auto &c: fChildren) {
      // TGLViewer* glv = ((REveViewer *)c)->GetGLViewer();
      // if ( fUseLightColorSet)
      //    glv->UseLightColorSet();
      // else
      //    glv->UseDarkColorSet();

      // glv->RequestDraw(TGLRnrCtx::kLODHigh);
   // }
}
